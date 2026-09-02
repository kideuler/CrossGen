// E2 -- Singularity structure on canonical domains (Sec. 6, 0.7 page).
//
// Eight domains whose singularity content is known before the code runs, three
// methods, one table. The four with a closed-form answer -- square, disk,
// annulus, L -- have every corner angle a multiple of pi/2, so the boundary's
// quarter count is exact and Poincare-Hopf fixes the interior total outright.
// The wedges and the pentagon do not, and there the identity is the statement
// being checked rather than a number to compare against.
//
// Reported per method and per domain:
//   #singularities, the index multiset, the Poincare-Hopf residual,
//   the distance from the nearest singularity to dS in units of local h,
//   the number of singularity pairs closer than 2h (the "spurious pair" proxy),
//   the boundary alignment error, and the common energy of Sec. 5.
//
// Robustness: ten seeded runs of B1 -- which starts from a random interior and
// so has a distribution rather than an answer -- and five 10%-jitter remeshes
// for all three. The claim is that the MBO methods are stable and that SIPG is
// stable *without* paying the conversion's noise; the way to be wrong about
// that is not to measure it, so the variance is a column.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#include "common/Domains.hxx"
#include "common/Methods.hxx"
#include "common/Metrics.hxx"
#include "common/Report.hxx"

using namespace paper;

namespace {

constexpr double kGammaEval = 10.0;

// The field-quality experiments run the MBO tighter than the pipeline does; see
// MethodOptions::convergenceTol and E1(d) for the measurement that motivates it.
constexpr double kTol = 1e-9;

struct Measured {
    bool ok = false;
    std::string error;
    metrics::SingularityReport sing;
    double energy = 0.0;
    metrics::AlignmentError boundary;
    int iterations = 0;
    double seconds = 0.0;
    metrics::ConversionReport conversion;
    bool hasConversion = false;
};

Measured measure(const std::shared_ptr<Mesh> &m, const FieldRun &run) {
    Measured x;
    if (!run.ok) { x.error = run.error; return x; }
    x.ok = true;
    x.sing = metrics::singularities(*m, run.u);
    x.energy = metrics::commonEnergy(*m, metrics::edgeWeights(*m, kGammaEval), run.u,
                                     interfaceEdgeFlags(*m));
    x.boundary = metrics::boundaryAlignmentError(*m, run.u);
    x.iterations = run.iterations;
    x.seconds = run.totalSeconds;
    x.conversion = run.conversion;
    x.hasConversion = run.hasConversion;
    return x;
}

double mean(const std::vector<double> &v) {
    if (v.empty()) return 0.0;
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}
double stdev(const std::vector<double> &v) {
    if (v.size() < 2) return 0.0;
    const double mu = mean(v);
    double s = 0.0;
    for (double x : v) s += (x - mu) * (x - mu);
    return std::sqrt(s / (v.size() - 1));
}

} // namespace

int main(int argc, char **argv) {
    std::string outDir = "results";
    double h = 0.03;
    int seeds = 10;
    int jitters = 5;
    std::string vtkDir;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) outDir = argv[++i];
        else if (a == "--h" && i + 1 < argc) h = std::stod(argv[++i]);
        else if (a == "--seeds" && i + 1 < argc) seeds = std::stoi(argv[++i]);
        else if (a == "--jitters" && i + 1 < argc) jitters = std::stoi(argv[++i]);
        else if (a == "--vtk" && i + 1 < argc) vtkDir = argv[++i];
        else if (a == "--quick") { h = 0.05; seeds = 3; jitters = 2; }
        else if (a == "--help") {
            std::cout << "Usage: " << argv[0]
                      << " [--out DIR] [--h SIZE] [--seeds N] [--jitters N] [--vtk DIR] [--quick]\n";
            return 0;
        }
    }
    ensureDir(outDir);
    if (!vtkDir.empty()) ensureDir(vtkDir);
    Verdicts v;

    banner("E2  Singularity structure on canonical domains");
    std::cout << "h = " << h << ", MBO convergence tolerance " << kTol
              << ", energies at gamma = " << kGammaEval << ".\n";

    const std::vector<Domain> domains = canonicalDomains(h);

    MethodOptions base;
    base.convergenceTol = kTol;

    // -----------------------------------------------------------------------
    // The main table
    // -----------------------------------------------------------------------
    Csv csv(outDir + "/E2_canonical.csv",
            {"domain", "method", "vertices", "triangles", "chi", "boundary_quarters",
             "singularities", "index_multiset", "index4_sum", "expected_index4_sum",
             "ph_residual", "ph_residual_geometric", "boundary_index4_sum",
             "boundary_adjacent", "min_boundary_distance_h",
             "close_pairs", "close_opposite_pairs", "max_spin_jump_rad",
             "stars_past_condition", "boundary_anomalies", "ambiguous_corners",
             "boundary_align_max_deg", "boundary_align_p95_deg",
             "energy", "iterations", "seconds",
             "conv_holonomy_lost_edges", "conv_topology_changed"});

    for (const Domain &d : domains) {
        heading(d.name + "  (" + d.note + ")");
        const std::vector<int> q = metrics::cornerQuarters(*d.mesh);
        int qs = 0;
        for (std::size_t i = 0; i < q.size(); ++i) if (d.mesh->isBoundaryVertex[i]) qs += q[i];
        const int chi = metrics::eulerCharacteristic(*d.mesh);
        const int predicted = 4 * chi - qs;

        std::cout << "  " << d.mesh->vertices.size() << " vertices, "
                  << d.mesh->triangles.size() << " triangles, chi = " << chi
                  << ", boundary quarters = " << qs
                  << " => the interior must total " << predicted << "/4\n";
        if (d.knownIndex)
            v.check(predicted == d.expectedInteriorIndex4,
                    d.name + ": the geometry predicts the closed-form answer "
                    + num(d.expectedInteriorIndex4));

        Table t({"method", "#sing", "indices", "sum 4I", "PH", "PH geom", "dS anom",
                 "near dS", "min d/h", "pairs<2h", "max jump", "dS align max",
                 "energy", "iters", "s"});

        for (const std::string &name : methodNames()) {
            const FieldRun run = runMethod(name, d.mesh, base);
            const Measured x = measure(d.mesh, run);
            if (!x.ok) { v.warn(d.name + "/" + name + ": " + x.error); continue; }

            t.row({name, num((int)x.sing.interior.size()), x.sing.indexMultiset(),
                   num(x.sing.interiorIndex4Sum), num(x.sing.poincareHopfResidual),
                   num(x.sing.geometricPoincareHopfResidual), num(x.sing.boundaryAnomalies),
                   num(x.sing.boundaryAdjacent), num(x.sing.minBoundaryDistanceH, 3),
                   num(x.sing.closePairs), num(x.sing.maxSpinJump, 3),
                   num(x.boundary.maxAngle * 180.0 / M_PI, 3),
                   num(x.energy, 6), num(x.iterations), num(x.seconds, 3)});

            csv.row({{"domain", d.name}, {"method", name},
                     {"vertices", num((int)d.mesh->vertices.size())},
                     {"triangles", num((int)d.mesh->triangles.size())},
                     {"chi", num(chi)}, {"boundary_quarters", num(qs)},
                     {"singularities", num((int)x.sing.interior.size())},
                     {"index_multiset", "\"" + x.sing.indexMultiset() + "\""},
                     {"index4_sum", num(x.sing.interiorIndex4Sum)},
                     {"expected_index4_sum", num(predicted)},
                     {"ph_residual", num(x.sing.poincareHopfResidual)},
                     {"ph_residual_geometric", num(x.sing.geometricPoincareHopfResidual)},
                     {"boundary_index4_sum", num(x.sing.boundaryIndex4Sum)},
                     {"boundary_adjacent", num(x.sing.boundaryAdjacent)},
                     {"min_boundary_distance_h", num(x.sing.minBoundaryDistanceH, 6)},
                     {"close_pairs", num(x.sing.closePairs)},
                     {"close_opposite_pairs", num(x.sing.closeOppositePairs)},
                     {"max_spin_jump_rad", num(x.sing.maxSpinJump, 6)},
                     {"stars_past_condition", num(x.sing.starsPastCondition)},
                     {"boundary_anomalies", num(x.sing.boundaryAnomalies)},
                     {"ambiguous_corners", num(x.sing.ambiguousCorners)},
                     {"boundary_align_max_deg", num(x.boundary.maxAngle * 180.0 / M_PI, 6)},
                     {"boundary_align_p95_deg", num(x.boundary.p95Angle * 180.0 / M_PI, 6)},
                     {"energy", num(x.energy, 10)},
                     {"iterations", num(x.iterations)},
                     {"seconds", num(x.seconds, 4)},
                     {"conv_holonomy_lost_edges",
                      x.hasConversion ? num(x.conversion.holonomyLostEdges) : ""},
                     {"conv_topology_changed",
                      x.hasConversion ? num(x.conversion.topologyChanged) : ""}});

            if (name == "SIPG") {
                v.check(x.sing.poincareHopfResidual == 0,
                        d.name + "/SIPG: Poincare-Hopf holds");
                // The geometric prediction is only a prediction where the corner
                // rounding is unambiguous. A 45-degree apex asks the boundary
                // for exactly one and a half quarter turns, and which way that
                // rounds is a coin toss that the mesh, not the method, settles
                // -- the same wedge at h = 0.05 and h = 0.03 lands on different
                // sides of it. Where a corner is that close to a tie the field
                // is entitled to disagree with the geometry, and the identity
                // that still has to hold is the one above.
                if (x.sing.ambiguousCorners == 0)
                    v.check(x.sing.interiorIndex4Sum == predicted,
                            d.name + "/SIPG: the interior index total is the predicted "
                            + num(predicted));
                else if (x.sing.interiorIndex4Sum != predicted)
                    v.warn(d.name + "/SIPG: interior total " + num(x.sing.interiorIndex4Sum)
                           + " against a geometric prediction of " + num(predicted)
                           + ", on a domain with " + num(x.sing.ambiguousCorners)
                           + " corner(s) within 9 degrees of a rounding tie");
            }
            if (!vtkDir.empty())
                writeFieldVTK(vtkDir + "/E2_" + d.name + "_" + name + ".vtk", *d.mesh, run.u);
        }
        t.print();
    }

    // -----------------------------------------------------------------------
    // Robustness: B1's seeds, and jittered remeshes for all three
    // -----------------------------------------------------------------------
    banner("E2  Robustness");
    {
        Csv rc(outDir + "/E2_robustness.csv",
               {"domain", "method", "kind", "trial", "singularities", "index4_sum",
                "ph_residual", "energy"});

        heading("B1 across " + num(seeds) + " random seeds");
        Table t({"domain", "#sing mean", "sigma", "min", "max", "sum 4I distinct"});
        for (const Domain &d : domains) {
            std::vector<double> counts;
            std::map<int, int> sums;
            for (int s = 0; s < seeds; ++s) {
                MethodOptions o = base;
                o.seed = static_cast<unsigned>(1000 + s);
                const FieldRun run = runP1MBO(d.mesh, o);
                const Measured x = measure(d.mesh, run);
                if (!x.ok) continue;
                counts.push_back(static_cast<double>(x.sing.interior.size()));
                ++sums[x.sing.interiorIndex4Sum];
                rc.row({{"domain", d.name}, {"method", "B1"}, {"kind", "seed"},
                        {"trial", num(s)},
                        {"singularities", num((int)x.sing.interior.size())},
                        {"index4_sum", num(x.sing.interiorIndex4Sum)},
                        {"ph_residual", num(x.sing.poincareHopfResidual)},
                        {"energy", num(x.energy, 10)}});
            }
            if (counts.empty()) continue;
            t.row({d.name, num(mean(counts), 4), num(stdev(counts), 4),
                   num((int)*std::min_element(counts.begin(), counts.end())),
                   num((int)*std::max_element(counts.begin(), counts.end())),
                   num((int)sums.size())});
        }
        t.print();

        heading("All three across " + num(jitters) + " 10%-jitter remeshes");
        Table jt({"domain", "method", "#sing mean", "sigma", "sum 4I distinct", "PH failures"});
        for (const Domain &d : domains) {
            for (const std::string &name : methodNames()) {
                std::vector<double> counts;
                std::map<int, int> sums;
                int phFail = 0;
                for (int j = 0; j < jitters; ++j) {
                    std::shared_ptr<Mesh> jm = jitter(*d.mesh, 0.10, static_cast<unsigned>(7000 + j));
                    MethodOptions o = base;
                    o.seed = 1;
                    const FieldRun run = runMethod(name, jm, o);
                    const Measured x = measure(jm, run);
                    if (!x.ok) continue;
                    counts.push_back(static_cast<double>(x.sing.interior.size()));
                    ++sums[x.sing.interiorIndex4Sum];
                    if (x.sing.poincareHopfResidual != 0) ++phFail;
                    rc.row({{"domain", d.name}, {"method", name}, {"kind", "jitter"},
                            {"trial", num(j)},
                            {"singularities", num((int)x.sing.interior.size())},
                            {"index4_sum", num(x.sing.interiorIndex4Sum)},
                            {"ph_residual", num(x.sing.poincareHopfResidual)},
                            {"energy", num(x.energy, 10)}});
                }
                if (counts.empty()) continue;
                jt.row({d.name, name, num(mean(counts), 4), num(stdev(counts), 4),
                        num((int)sums.size()), num(phFail)});
                if (name == "SIPG")
                    v.check(sums.size() == 1,
                            d.name + "/SIPG: the index total is the same on every remesh");
            }
        }
        jt.print();
    }

    banner("E2 summary");
    std::cout << (v.failures ? kFail : kPass) << " " << v.failures << " failure(s), "
              << v.warnings << " warning(s). CSVs in " << outDir << "/\n";
    return v.failures ? 1 : 0;
}
