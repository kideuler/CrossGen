// E3 -- Benchmark sweep over the single-material corpus (Sec. 6, 0.85 page).
//
// All three methods on every model in data/meshes/singlemat, with the aggregate
// table the paper prints and the per-model table the supplement carries.
//
// The column the outline calls the paper's key number is the **conversion cost**
// of B1, and it is two things measured separately:
//
//   holonomy_lost_edges   interior edges where the resampled face field's jump
//                         is not the P1 field's own transport across that edge.
//                         The P1 field is continuous, so the transport is
//                         determined; the face field's jump is that transport
//                         reduced to (-pi, pi]. Where they differ by a full
//                         spin turn, the matching the quantizer reads on that
//                         edge is not the field's holonomy.
//   topology_changed      singularities created or annihilated by the
//                         resampling, matched greedily at equal index within
//                         three local edge lengths.
//
// If either is materially non-zero across the corpus, the motivation of Sec. 1
// is not an argument but a measurement. If it is zero, the honest thing is to
// say so and lean on E4 instead -- so this experiment prints the number either
// way and asserts nothing about which way it comes out.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#include "common/Corpus.hxx"
#include "common/Methods.hxx"
#include "common/Metrics.hxx"
#include "common/Report.hxx"

using namespace paper;

namespace {

constexpr double kGammaEval = 10.0;
constexpr double kTol = 1e-9;   // see E1(d); the field-quality experiments run tight

struct Row {
    std::string model, method;
    bool ok = false;
    int vertices = 0, triangles = 0, chi = 0;
    int singularities = 0, index4Sum = 0, phResidual = 0, phGeom = 0;
    int boundaryIndex4Sum = 0, boundaryAnomalies = 0, ambiguousCorners = 0;
    int boundaryAdjacent = 0, closePairs = 0, closeOpposite = 0;
    double minBoundaryDistanceH = 0.0;
    double maxSpinJump = 0.0;
    int starsPastCondition = 0;
    double energy = 0.0;
    double boundaryAlignMaxDeg = 0.0, boundaryAlignP95Deg = 0.0;
    int iterations = 0;
    double assembleSeconds = 0.0, solveSeconds = 0.0, totalSeconds = 0.0;
    bool converged = false;
    metrics::ConversionReport conv;
    bool hasConv = false;
};

double mean(const std::vector<double> &v) {
    return v.empty() ? 0.0 : std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}
double median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

} // namespace

int main(int argc, char **argv) {
    std::string outDir = "results";
    std::string sub = "singlemat";
    std::string vtkDir;
    std::vector<std::string> only;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) outDir = argv[++i];
        else if (a == "--set" && i + 1 < argc) sub = argv[++i];
        else if (a == "--vtk" && i + 1 < argc) vtkDir = argv[++i];
        else if (a == "--help") {
            std::cout << "Usage: " << argv[0]
                      << " [--out DIR] [--set singlemat|multimat] [--vtk DIR] [model ...]\n";
            return 0;
        } else only.push_back(a);
    }
    ensureDir(outDir);
    if (!vtkDir.empty()) ensureDir(vtkDir);
    Verdicts v;

    banner("E3  Benchmark sweep: data/meshes/" + sub);

    std::vector<Model> models = select(corpus(sub), only);
    if (models.empty()) {
        std::cout << kFail << " no models found under " << PAPER_MESH_DIR << "/" << sub << "\n";
        return 2;
    }
    std::cout << models.size() << " model(s); MBO convergence tolerance " << kTol
              << ", energies at gamma = " << kGammaEval << ".\n";

    Csv csv(outDir + "/E3_" + sub + "_permodel.csv",
            {"model", "method", "vertices", "triangles", "chi",
             "singularities", "index4_sum", "ph_residual", "ph_residual_geometric",
             "boundary_index4_sum", "boundary_anomalies", "ambiguous_corners",
             "boundary_adjacent",
             "min_boundary_distance_h", "close_pairs", "close_opposite_pairs",
             "max_spin_jump_rad", "stars_past_condition",
             "energy", "boundary_align_max_deg", "boundary_align_p95_deg",
             "iterations", "converged", "assemble_s", "solve_s", "total_s",
             "conv_degenerate_faces", "conv_holonomy_lost_edges", "conv_untrackable_edges",
             "conv_singularities_before", "conv_singularities_after",
             "conv_created", "conv_annihilated", "conv_topology_changed"});

    std::vector<Row> rows;
    int loadFailures = 0;

    MethodOptions base;
    base.convergenceTol = kTol;

    for (const Model &mm : models) {
        std::string why;
        std::shared_ptr<Mesh> mesh = load(mm, why);
        if (!mesh) { v.warn(mm.name + ": " + why); ++loadFailures; continue; }

        const metrics::EdgeWeights w = metrics::edgeWeights(*mesh, kGammaEval);
        const std::vector<char> skip = interfaceEdgeFlags(*mesh);
        const int chi = metrics::eulerCharacteristic(*mesh);

        std::cout << "\n" << mm.name << "  (" << mesh->vertices.size() << " v, "
                  << mesh->triangles.size() << " t, chi = " << chi << ")\n";
        Table t({"method", "#sing", "indices", "sum 4I", "PH", "PH geom", "dS anom",
                 "near dS", "min d/h", "pairs<2h", "energy", "dS align max", "iters", "s"});

        for (const std::string &name : methodNames()) {
            const FieldRun run = runMethod(name, mesh, base);
            Row r;
            r.model = mm.name;
            r.method = name;
            r.vertices = static_cast<int>(mesh->vertices.size());
            r.triangles = static_cast<int>(mesh->triangles.size());
            r.chi = chi;
            if (!run.ok) {
                v.warn(mm.name + "/" + name + ": " + run.error);
                rows.push_back(r);
                continue;
            }
            const metrics::SingularityReport s = metrics::singularities(*mesh, run.u);
            const metrics::AlignmentError ba = metrics::boundaryAlignmentError(*mesh, run.u);

            r.ok = true;
            r.singularities = static_cast<int>(s.interior.size());
            r.index4Sum = s.interiorIndex4Sum;
            r.phResidual = s.poincareHopfResidual;
            r.phGeom = s.geometricPoincareHopfResidual;
            r.boundaryIndex4Sum = s.boundaryIndex4Sum;
            r.boundaryAnomalies = s.boundaryAnomalies;
            r.ambiguousCorners = s.ambiguousCorners;
            r.boundaryAdjacent = s.boundaryAdjacent;
            r.closePairs = s.closePairs;
            r.closeOpposite = s.closeOppositePairs;
            r.minBoundaryDistanceH = s.minBoundaryDistanceH;
            r.maxSpinJump = s.maxSpinJump;
            r.starsPastCondition = s.starsPastCondition;
            r.energy = metrics::commonEnergy(*mesh, w, run.u, skip);
            r.boundaryAlignMaxDeg = ba.maxAngle * 180.0 / M_PI;
            r.boundaryAlignP95Deg = ba.p95Angle * 180.0 / M_PI;
            r.iterations = run.iterations;
            r.converged = run.converged;
            r.assembleSeconds = run.assembleSeconds;
            r.solveSeconds = run.solveSeconds;
            r.totalSeconds = run.totalSeconds;
            r.conv = run.conversion;
            r.hasConv = run.hasConversion;
            rows.push_back(r);

            t.row({name, num(r.singularities), s.indexMultiset(), num(r.index4Sum),
                   num(r.phResidual), num(r.phGeom), num(r.boundaryAnomalies),
                   num(r.boundaryAdjacent),
                   num(r.minBoundaryDistanceH, 3), num(r.closePairs),
                   num(r.energy, 6), num(r.boundaryAlignMaxDeg, 3),
                   num(r.iterations), num(r.totalSeconds, 3)});

            csv.row({{"model", r.model}, {"method", r.method},
                     {"vertices", num(r.vertices)}, {"triangles", num(r.triangles)},
                     {"chi", num(r.chi)},
                     {"singularities", num(r.singularities)},
                     {"index4_sum", num(r.index4Sum)},
                     {"ph_residual", num(r.phResidual)},
                     {"ph_residual_geometric", num(r.phGeom)},
                     {"boundary_index4_sum", num(r.boundaryIndex4Sum)},
                     {"boundary_anomalies", num(r.boundaryAnomalies)},
                     {"ambiguous_corners", num(r.ambiguousCorners)},
                     {"boundary_adjacent", num(r.boundaryAdjacent)},
                     {"min_boundary_distance_h", num(r.minBoundaryDistanceH, 6)},
                     {"close_pairs", num(r.closePairs)},
                     {"close_opposite_pairs", num(r.closeOpposite)},
                     {"max_spin_jump_rad", num(r.maxSpinJump, 6)},
                     {"stars_past_condition", num(r.starsPastCondition)},
                     {"energy", num(r.energy, 10)},
                     {"boundary_align_max_deg", num(r.boundaryAlignMaxDeg, 6)},
                     {"boundary_align_p95_deg", num(r.boundaryAlignP95Deg, 6)},
                     {"iterations", num(r.iterations)},
                     {"converged", num(r.converged)},
                     {"assemble_s", num(r.assembleSeconds, 5)},
                     {"solve_s", num(r.solveSeconds, 5)},
                     {"total_s", num(r.totalSeconds, 5)},
                     {"conv_degenerate_faces", r.hasConv ? num(r.conv.degenerateFaces) : ""},
                     {"conv_holonomy_lost_edges", r.hasConv ? num(r.conv.holonomyLostEdges) : ""},
                     {"conv_untrackable_edges", r.hasConv ? num(r.conv.untrackableEdges) : ""},
                     {"conv_singularities_before", r.hasConv ? num(r.conv.singularitiesBefore) : ""},
                     {"conv_singularities_after", r.hasConv ? num(r.conv.singularitiesAfter) : ""},
                     {"conv_created", r.hasConv ? num(r.conv.created) : ""},
                     {"conv_annihilated", r.hasConv ? num(r.conv.annihilated) : ""},
                     {"conv_topology_changed", r.hasConv ? num(r.conv.topologyChanged) : ""}});

            if (!vtkDir.empty())
                writeFieldVTK(vtkDir + "/E3_" + mm.name + "_" + name + ".vtk", *mesh, run.u);
        }
        t.print();
    }

    // -----------------------------------------------------------------------
    // Aggregates
    // -----------------------------------------------------------------------
    banner("E3  Aggregate over " + num((int)models.size() - loadFailures) + " model(s)");

    Csv agg(outDir + "/E3_" + sub + "_aggregate.csv",
            {"method", "models", "solved", "ph_ok", "singularities_mean",
             "ph_geom_ok", "boundary_adjacent_total", "close_pairs_total",
             "energy_ratio_vs_dualmbo_mean", "energy_ratio_vs_dualmbo_median",
             "boundary_align_max_deg_worst", "iterations_mean", "total_s_mean"});

    // The energy ratio is per model, so a model nobody solved drops out of all
    // three rather than out of one.
    std::map<std::string, double> dualMBOEnergy;
    for (const Row &r : rows) if (r.ok && r.method == "DualMBO") dualMBOEnergy[r.model] = r.energy;

    Table at({"method", "solved", "PH ok", "PH geom ok", "#sing mean", "near-dS total", "pairs<2h total",
              "E/E_DualMBO mean", "median", "worst dS align", "iters mean", "s mean"});

    for (const std::string &name : methodNames()) {
        int solved = 0, phOk = 0, phGeomOk = 0, nearTotal = 0, pairTotal = 0;
        std::vector<double> sings, ratios, iters, secs;
        double worstAlign = 0.0;
        for (const Row &r : rows) {
            if (r.method != name || !r.ok) continue;
            ++solved;
            if (r.phResidual == 0) ++phOk;
            if (r.phGeom == 0) ++phGeomOk;
            nearTotal += r.boundaryAdjacent;
            pairTotal += r.closePairs;
            sings.push_back(r.singularities);
            iters.push_back(r.iterations);
            secs.push_back(r.totalSeconds);
            worstAlign = std::max(worstAlign, r.boundaryAlignMaxDeg);
            auto it = dualMBOEnergy.find(r.model);
            if (it != dualMBOEnergy.end() && it->second > 1e-12)
                ratios.push_back(r.energy / it->second);
        }
        at.row({name, num(solved), num(phOk), num(phGeomOk), num(mean(sings), 4), num(nearTotal),
                num(pairTotal), num(mean(ratios), 4), num(median(ratios), 4),
                num(worstAlign, 3), num(mean(iters), 4), num(mean(secs), 4)});
        agg.row({{"method", name}, {"models", num((int)models.size() - loadFailures)},
                 {"solved", num(solved)}, {"ph_ok", num(phOk)}, {"ph_geom_ok", num(phGeomOk)},
                 {"singularities_mean", num(mean(sings), 6)},
                 {"boundary_adjacent_total", num(nearTotal)},
                 {"close_pairs_total", num(pairTotal)},
                 {"energy_ratio_vs_dualmbo_mean", num(mean(ratios), 6)},
                 {"energy_ratio_vs_dualmbo_median", num(median(ratios), 6)},
                 {"boundary_align_max_deg_worst", num(worstAlign, 6)},
                 {"iterations_mean", num(mean(iters), 6)},
                 {"total_s_mean", num(mean(secs), 6)}});
    }
    at.print();

    // -----------------------------------------------------------------------
    // The conversion cost
    // -----------------------------------------------------------------------
    heading("The cost of the P1 -> face conversion (B1)");
    {
        int modelsWithLoss = 0, modelsWithTopoChange = 0;
        int lostTotal = 0, createdTotal = 0, annihilatedTotal = 0, degenerateTotal = 0;
        int interiorEdgeTotal = 0;
        Table ct({"model", "interior edges", "holonomy lost", "%", "sing before", "after",
                  "created", "annihilated", "degenerate faces"});
        for (const Row &r : rows) {
            if (r.method != "B1" || !r.ok || !r.hasConv) continue;
            // Interior edges: 3T minus the boundary edges, halved.
            const int approxInterior = (3 * r.triangles - (2 * r.triangles + 2 - 2 * r.vertices)) / 2;
            interiorEdgeTotal += approxInterior;
            lostTotal += r.conv.holonomyLostEdges;
            createdTotal += r.conv.created;
            annihilatedTotal += r.conv.annihilated;
            degenerateTotal += r.conv.degenerateFaces;
            if (r.conv.holonomyLostEdges > 0) ++modelsWithLoss;
            if (r.conv.topologyChanged > 0) ++modelsWithTopoChange;
            if (r.conv.holonomyLostEdges > 0 || r.conv.topologyChanged > 0) {
                ct.row({r.model, num(approxInterior), num(r.conv.holonomyLostEdges),
                        num(100.0 * r.conv.holonomyLostEdges / std::max(1, approxInterior), 3),
                        num(r.conv.singularitiesBefore), num(r.conv.singularitiesAfter),
                        num(r.conv.created), num(r.conv.annihilated),
                        num(r.conv.degenerateFaces)});
            }
        }
        if (modelsWithLoss || modelsWithTopoChange) ct.print();
        else std::cout << "  (no model lost a matching or changed its singularity set)\n";

        std::cout << "\n  " << modelsWithLoss << " model(s) lost at least one edge matching; "
                  << lostTotal << " edge(s) in total.\n"
                  << "  " << modelsWithTopoChange
                  << " model(s) changed singularity set; " << createdTotal << " created, "
                  << annihilatedTotal << " annihilated.\n"
                  << "  " << degenerateTotal
                  << " face(s) had no direction to average at all.\n";

        Csv cc(outDir + "/E3_" + sub + "_conversion.csv",
               {"models_with_holonomy_loss", "holonomy_lost_edges_total",
                "models_with_topology_change", "created_total", "annihilated_total",
                "degenerate_faces_total", "interior_edges_total"});
        cc.row({num(modelsWithLoss), num(lostTotal), num(modelsWithTopoChange),
                num(createdTotal), num(annihilatedTotal), num(degenerateTotal),
                num(interiorEdgeTotal)});
    }

    // -----------------------------------------------------------------------
    // Verdicts. Only DualMBO is asserted about: it is our method and Prop. 3 is
    // ours to hold. The baselines are reported.
    // -----------------------------------------------------------------------
    heading("Verdicts");
    for (const Row &r : rows) {
        if (r.method != "DualMBO") continue;
        v.check(r.ok, r.model + ": DualMBO produced a field");
        if (r.ok) v.check(r.phResidual == 0, r.model + ": Poincare-Hopf holds for DualMBO");
    }
    {
        int b1Bad = 0, b2Bad = 0;
        for (const Row &r : rows) {
            if (!r.ok) continue;
            if (r.method == "B1" && r.phResidual != 0) ++b1Bad;
            if (r.method == "B2" && r.phResidual != 0) ++b2Bad;
        }
        std::cout << "  (reported, not asserted) Poincare-Hopf fails on " << b1Bad
                  << " B1 run(s) and " << b2Bad << " B2 run(s).\n";
    }

    banner("E3 summary");
    std::cout << (v.failures ? kFail : kPass) << " " << v.failures << " failure(s), "
              << v.warnings << " warning(s). CSVs in " << outDir << "/\n";
    return v.failures ? 1 : 0;
}
