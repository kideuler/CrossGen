// E4 -- End-to-end quad layout (Sec. 6, 1.15 pages). The money experiment.
//
// The whole pipeline, twice, with one thing changed: the cross field. Stages 0b
// to 11 of TORSION -- interface network, cone set, cut, frames, integration,
// layout energies, separatrices, arrangement, spline fit, meshing -- are the
// same code, on the same mesh, at the same target edge length, with the same
// seeding tolerance. What differs is the vector of per-face crosses Stage 0
// would otherwise have solved for, and what is measured is the layout that came
// out the far end.
//
// Reported per model and method:
//   whether the pipeline reached a mesh, and which stage it stopped at,
//   #patches, #irregular vertices (interior and boundary),
//   min/mean scaled Jacobian and the angle-deviation summary,
//   the minimum zone dimension as a CFL proxy (R3/R4), and the wall time.
//
// Two things this experiment is careful about.
//
//   * It is *conservative by construction*. Stage 1 does not take the field's
//     singularities as given: it prescribes cones on the interface network,
//     cancels dipoles inside a region and rebalances to satisfy Gauss-Bonnet.
//     A field whose indices are wrong is therefore partly repaired before the
//     layout sees them, so a difference that survives to the mesh understates
//     what the field cost rather than overstating it.
//   * Not every model reaches a mesh on either field. That is stated in the
//     outline and it is a column here, not an omission: "reached" says where a
//     run stopped, and the aggregate counts successes rather than assuming them.
//
// The last section lists **figure candidates** -- the models where the two
// fields disagree most, split by kind -- because choosing the figures out of a
// 35-row table by eye is the slowest part of writing Sec. 6.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "common/Corpus.hxx"
#include "common/Layout.hxx"
#include "common/Methods.hxx"
#include "common/Metrics.hxx"
#include "common/Report.hxx"

using namespace paper;

namespace {

struct Outcome {
    std::string model, method;
    bool fieldOk = false;
    LayoutResult layout;
};

double mean(const std::vector<double> &v) {
    return v.empty() ? 0.0 : std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

std::vector<std::string> split(const std::string &s, char sep) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, sep)) if (!tok.empty()) out.push_back(tok);
    return out;
}

} // namespace

int main(int argc, char **argv) {
    std::string outDir = "results";
    std::string sub = "singlemat";
    std::string objDir;
    std::vector<std::string> methods{"SIPG", "B1"};
    std::vector<std::string> only;
    double target = 0.05;
    bool diskTemplates = false;
    bool b1Continuation = false;
    // The p=0 penalty weight, which at this order *is* the discrete Laplacian.
    // "orth" is the two-point/finite-volume weight; see SIPG::PenaltyWeight.
    std::string weight = "min";
    int limit = 0;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) outDir = argv[++i];
        else if (a == "--set" && i + 1 < argc) sub = argv[++i];
        else if (a == "--methods" && i + 1 < argc) methods = split(argv[++i], ',');
        else if (a == "--target" && i + 1 < argc) target = std::stod(argv[++i]);
        else if (a == "--obj" && i + 1 < argc) objDir = argv[++i];
        else if (a == "--limit" && i + 1 < argc) limit = std::stoi(argv[++i]);
        else if (a == "--disk-templates") diskTemplates = true;
        else if (a == "--b1-continuation") b1Continuation = true;
        else if (a == "--weight" && i + 1 < argc) weight = argv[++i];
        else if (a == "--help") {
            std::cout << "Usage: " << argv[0]
                      << " [--out DIR] [--set singlemat|multimat] [--methods SIPG,B1,B2]\n"
                         "       [--target EDGE] [--obj DIR] [--limit N] [--disk-templates]\n"
                         "       [--b1-continuation] [--weight min|harm|orth] [model ...]\n";
            return 0;
        } else only.push_back(a);
    }
    ensureDir(outDir);
    if (!objDir.empty()) ensureDir(objDir);
    Verdicts v;

    banner("E4  End-to-end quad layout through the pipeline: data/meshes/" + sub);

    std::vector<Model> models = select(corpus(sub), only);
    if (limit > 0 && static_cast<int>(models.size()) > limit) models.resize(limit);
    if (models.empty()) {
        std::cout << kFail << " no models found under " << PAPER_MESH_DIR << "/" << sub << "\n";
        return 2;
    }
    std::cout << models.size() << " model(s), methods:";
    for (const std::string &m : methods) std::cout << " " << m;
    std::cout << "; target edge " << target << "\n"
              << "The MBO runs at the tolerance the pipeline ships with (1e-5), because what\n"
              << "is being measured here is the pipeline.\n";

    Csv csv(outDir + "/E4_" + sub + "_permodel.csv",
            {"model", "method", "vertices", "triangles",
             "field_ok", "pipeline_ran", "valid", "reached",
             "frames_valid", "immersion_valid", "layout_valid", "arrangement_valid",
             "splines_valid", "mesh_valid", "mesh_conforming",
             "interior_cones", "boundary_cones", "integration_flipped_faces",
             "separatrices", "separatrices_unresolved",
             "patches", "coverage", "unmeshed_patches",
             "mesh_vertices", "mesh_quads",
             "irregular_interior", "irregular_boundary", "inverted_quads",
             "min_scaled_jacobian", "mean_scaled_jacobian",
             "pipeline_min_scaled_jacobian",
             "max_angle_deviation_deg", "mean_angle_deviation_deg",
             "min_zone_dimension", "p01_zone_dimension", "mean_zone_dimension",
             "mixed_quads", "pipeline_mixed_quads", "seconds", "error"});

    std::vector<Outcome> outcomes;

    for (const Model &mm : models) {
        std::string why;
        std::shared_ptr<Mesh> mesh = load(mm, why);
        if (!mesh) { v.warn(mm.name + ": " + why); continue; }

        std::cout << "\n" << mm.name << "  (" << mesh->vertices.size() << " v, "
                  << mesh->triangles.size() << " t)" << std::flush;
        Table t({"method", "reached", "valid", "patches", "quads", "irr int", "irr bnd",
                 "min SJ", "mean SJ", "inv", "min zone", "1% zone", "s"});

        for (const std::string &name : methods) {
            MethodOptions fo;   // the shipped tolerance: this measures the pipeline
            // The tau-continuation is part of our method; --b1-continuation
            // gives the baseline the same ladder, which is the fair-comparison
            // variant the protocol paragraph has to report beside the headline.
            fo.b1TauContinuation = b1Continuation;
            fo.penaltyWeight = weight == "orth" ? SIPG::PenaltyWeight::Orthogonal
                             : weight == "harm" ? SIPG::PenaltyWeight::HarmonicHeight
                                                : SIPG::PenaltyWeight::MinHeight;
            const FieldRun run = runMethod(name, mesh, fo);

            Outcome o;
            o.model = mm.name;
            o.method = name;
            o.fieldOk = run.ok;
            if (!run.ok) {
                v.warn(mm.name + "/" + name + ": field failed: " + run.error);
                outcomes.push_back(o);
                continue;
            }

            LayoutOptions lo;
            lo.quadTargetEdge = target;
            lo.diskTemplates = diskTemplates;
            o.layout = runLayout(mesh, run.u, lo, !objDir.empty());
            outcomes.push_back(o);

            const LayoutResult &L = o.layout;
            const metrics::QuadMetrics &Q = L.quality;
            t.row({name, L.reachedStage, L.valid ? "yes" : "no", num(L.patches),
                   num(Q.quads), num(Q.irregularInterior), num(Q.irregularBoundary),
                   num(Q.minScaledJacobian, 4), num(Q.meanScaledJacobian, 4),
                   num(Q.invertedQuads), num(Q.minZoneDimension, 4),
                   num(Q.p01ZoneDimension, 4), num(L.seconds, 2)});

            csv.row({{"model", mm.name}, {"method", name},
                     {"vertices", num((int)mesh->vertices.size())},
                     {"triangles", num((int)mesh->triangles.size())},
                     {"field_ok", num(true)}, {"pipeline_ran", num(L.ran)},
                     {"valid", num(L.valid)}, {"reached", L.reachedStage},
                     {"frames_valid", num(L.framesValid)},
                     {"immersion_valid", num(L.immersionValid)},
                     {"layout_valid", num(L.layoutValid)},
                     {"arrangement_valid", num(L.arrangementValid)},
                     {"splines_valid", num(L.splinesValid)},
                     {"mesh_valid", num(L.meshValid)},
                     {"mesh_conforming", num(L.meshConforming)},
                     {"interior_cones", num(L.interiorCones)},
                     {"boundary_cones", num(L.boundaryCones)},
                     {"integration_flipped_faces", num(L.integrationFlippedFaces)},
                     {"separatrices", num(L.separatrices)},
                     {"separatrices_unresolved", num(L.separatricesUnresolved)},
                     {"patches", num(L.patches)}, {"coverage", num(L.coverage, 6)},
                     {"unmeshed_patches", num(L.unmeshedPatches)},
                     {"mesh_vertices", num(L.meshVertices)}, {"mesh_quads", num(L.meshQuads)},
                     {"irregular_interior", num(Q.irregularInterior)},
                     {"irregular_boundary", num(Q.irregularBoundary)},
                     {"inverted_quads", num(Q.invertedQuads)},
                     {"min_scaled_jacobian", num(Q.minScaledJacobian, 6)},
                     {"mean_scaled_jacobian", num(Q.meanScaledJacobian, 6)},
                     {"pipeline_min_scaled_jacobian", num(L.pipelineMinScaledJacobian, 6)},
                     {"max_angle_deviation_deg", num(Q.maxAngleDeviation, 6)},
                     {"mean_angle_deviation_deg", num(Q.meanAngleDeviation, 6)},
                     {"min_zone_dimension", num(Q.minZoneDimension, 6)},
                     {"p01_zone_dimension", num(Q.p01ZoneDimension, 6)},
                     {"mean_zone_dimension", num(Q.meanZoneDimension, 6)},
                     {"mixed_quads", num(Q.mixedQuads)},
                     {"pipeline_mixed_quads", num(L.pipelineMixedQuads)},
                     {"seconds", num(L.seconds, 3)},
                     {"error", L.error}});

            if (!objDir.empty() && !L.quadCells.empty())
                writeQuadOBJ(objDir + "/E4_" + mm.name + "_" + name + ".obj", L);
        }
        std::cout << "\n";
        t.print();
    }

    // -----------------------------------------------------------------------
    // Aggregate: where runs stop, and the paired comparison where both finished
    // -----------------------------------------------------------------------
    banner("E4  Aggregate");

    heading("Where the pipeline stopped");
    {
        std::map<std::string, std::map<std::string, int>> byStage;
        for (const Outcome &o : outcomes) byStage[o.method][o.layout.reachedStage]++;
        // The stages in pipeline order, so the table reads as a funnel.
        const std::vector<std::string> order{"none", "field", "cones", "cut", "frames",
                                             "psi_0", "layout(partial)", "layout",
                                             "separatrices", "arrangement", "splines",
                                             "mesh(partial)", "mesh"};
        std::vector<std::string> header{"method"};
        for (const std::string &s : order) header.push_back(s);
        Table t(header);
        for (const std::string &name : methods) {
            std::vector<std::string> row{name};
            for (const std::string &s : order) row.push_back(num(byStage[name][s]));
            t.row(row);
        }
        t.print();
    }

    heading("Success rate");
    Csv agg(outDir + "/E4_" + sub + "_aggregate.csv",
            {"method", "models", "reached_mesh", "valid_layout",
             "patches_mean", "irregular_interior_mean", "inverted_total",
             "min_scaled_jacobian_worst", "min_scaled_jacobian_mean",
             "min_zone_dimension_worst", "p01_zone_dimension_mean", "seconds_mean"});
    {
        Table t({"method", "reached mesh", "valid layout", "patches mean",
                 "irr int mean", "inverted total", "worst min SJ", "mean min SJ",
                 "worst min zone", "mean 1% zone", "s mean"});
        for (const std::string &name : methods) {
            int reached = 0, valid = 0, inverted = 0, n = 0;
            std::vector<double> patches, irr, sj, zone, p01, secs;
            double worstSJ = 1.0, worstZone = 1e30;
            for (const Outcome &o : outcomes) {
                if (o.method != name) continue;
                ++n;
                secs.push_back(o.layout.seconds);
                if (!o.layout.meshValid && o.layout.quality.quads == 0) continue;
                ++reached;
                if (o.layout.valid) ++valid;
                inverted += o.layout.quality.invertedQuads;
                patches.push_back(o.layout.patches);
                irr.push_back(o.layout.quality.irregularInterior);
                sj.push_back(o.layout.quality.minScaledJacobian);
                zone.push_back(o.layout.quality.minZoneDimension);
                p01.push_back(o.layout.quality.p01ZoneDimension);
                worstSJ = std::min(worstSJ, o.layout.quality.minScaledJacobian);
                worstZone = std::min(worstZone, o.layout.quality.minZoneDimension);
            }
            if (worstZone > 1e29) worstZone = 0.0;
            t.row({name, num(reached) + "/" + num(n), num(valid) + "/" + num(n),
                   num(mean(patches), 4), num(mean(irr), 4), num(inverted),
                   num(worstSJ, 4), num(mean(sj), 4), num(worstZone, 4),
                   num(mean(p01), 4), num(mean(secs), 3)});
            agg.row({{"method", name}, {"models", num(n)},
                     {"reached_mesh", num(reached)}, {"valid_layout", num(valid)},
                     {"patches_mean", num(mean(patches), 6)},
                     {"irregular_interior_mean", num(mean(irr), 6)},
                     {"inverted_total", num(inverted)},
                     {"min_scaled_jacobian_worst", num(worstSJ, 6)},
                     {"min_scaled_jacobian_mean", num(mean(sj), 6)},
                     {"min_zone_dimension_worst", num(worstZone, 6)},
                     {"p01_zone_dimension_mean", num(mean(p01), 6)},
                     {"seconds_mean", num(mean(secs), 6)}});
        }
        t.print();
    }

    // -----------------------------------------------------------------------
    // Figure candidates
    // -----------------------------------------------------------------------
    if (methods.size() >= 2) {
        const std::string A = methods[0], B = methods[1];
        heading("Figure candidates: " + A + " against " + B);

        std::map<std::string, const Outcome *> a, b;
        for (const Outcome &o : outcomes) {
            if (o.method == A) a[o.model] = &o;
            if (o.method == B) b[o.model] = &o;
        }

        auto meshed = [](const Outcome *o) { return o && o->layout.quality.quads > 0; };

        Table only({"model", A + " reached", B + " reached", "note"});
        Table better({"model", A + " min SJ", B + " min SJ", "gain",
                      A + " irr", B + " irr", A + " min zone", B + " min zone"});
        int bothMeshed = 0, aOnly = 0, bOnly = 0, neither = 0, same = 0;

        for (const auto &[model, oa] : a) {
            auto it = b.find(model);
            const Outcome *ob = it == b.end() ? nullptr : it->second;
            const bool ma = meshed(oa), mb = meshed(ob);
            if (ma && !mb) { ++aOnly; only.row({model, oa->layout.reachedStage,
                                                ob ? ob->layout.reachedStage : "-",
                                                A + " only"}); }
            else if (!ma && mb) { ++bOnly; only.row({model, oa->layout.reachedStage,
                                                     ob->layout.reachedStage, B + " only"}); }
            else if (!ma && !mb) { ++neither; only.row({model, oa->layout.reachedStage,
                                                        ob ? ob->layout.reachedStage : "-",
                                                        "neither"}); }
            else {
                ++bothMeshed;
                const double sa = oa->layout.quality.minScaledJacobian;
                const double sb = ob->layout.quality.minScaledJacobian;
                const double za = oa->layout.quality.minZoneDimension;
                const double zb = ob->layout.quality.minZoneDimension;
                const int ia = oa->layout.quality.irregularInterior;
                const int ib = ob->layout.quality.irregularInterior;
                const bool interesting =
                    std::fabs(sa - sb) > 0.05 || ia != ib ||
                    std::fabs(za - zb) > 0.1 * std::max(za, zb);
                if (!interesting) { ++same; continue; }
                better.row({model, num(sa, 4), num(sb, 4), num(sa - sb, 4),
                            num(ia), num(ib), num(za, 4), num(zb, 4)});
            }
        }

        std::cout << "  " << bothMeshed << " model(s) meshed on both fields, " << aOnly
                  << " on " << A << " only, " << bOnly << " on " << B << " only, "
                  << neither << " on neither.\n";
        if (aOnly || bOnly || neither) { std::cout << "\n"; only.print(); }
        std::cout << "\n  Of the " << bothMeshed << " that meshed on both, " << same
                  << " are indistinguishable and the rest are:\n";
        better.print();
        std::cout << "\n  Honest reading: where the two columns agree, the conversion is\n"
                  << "  harmless on that model, and the paper should say so.\n";

        // -------------------------------------------------------------------
        // The paired scoreboard.
        //
        // The aggregate above is a mean over a corpus whose models differ by
        // orders of magnitude, and a mean can be carried by one bad row. This
        // is the paired statement instead: on the models both fields meshed,
        // how often is each metric better, equal, or worse. It is the honest
        // form of "equal-or-better layouts", and if the answer is that they are
        // not, the paper needs to know that before a reviewer finds it.
        // -------------------------------------------------------------------
        heading("Paired scoreboard on the models both fields meshed");
        {
            struct Metric {
                std::string name;
                double (*get)(const metrics::QuadMetrics &);
                bool higherIsBetter;
                double tol;
            };
            const std::vector<Metric> mets{
                {"min scaled Jacobian",  [](const metrics::QuadMetrics &q) { return q.minScaledJacobian; },  true,  0.02},
                {"mean scaled Jacobian", [](const metrics::QuadMetrics &q) { return q.meanScaledJacobian; }, true,  0.002},
                {"inverted quads",       [](const metrics::QuadMetrics &q) { return (double)q.invertedQuads; }, false, 0.5},
                {"irregular interior",   [](const metrics::QuadMetrics &q) { return (double)q.irregularInterior; }, false, 0.5},
                {"min zone dimension",   [](const metrics::QuadMetrics &q) { return q.minZoneDimension; },   true,  0.01},
                {"1% zone dimension",    [](const metrics::QuadMetrics &q) { return q.p01ZoneDimension; },   true,  0.01},
                {"max angle deviation",  [](const metrics::QuadMetrics &q) { return q.maxAngleDeviation; },  false, 1.0},
            };

            Table t({"metric", A + " better", "tie", B + " better"});
            for (const Metric &me : mets) {
                int win = 0, tie = 0, loss = 0;
                for (const auto &[model, oa] : a) {
                    auto it = b.find(model);
                    if (it == b.end()) continue;
                    const Outcome *ob = it->second;
                    if (!meshed(oa) || !meshed(ob)) continue;
                    const double xa = me.get(oa->layout.quality);
                    const double xb = me.get(ob->layout.quality);
                    if (std::fabs(xa - xb) <= me.tol) { ++tie; continue; }
                    const bool aBetter = me.higherIsBetter ? (xa > xb) : (xa < xb);
                    if (aBetter) ++win; else ++loss;
                }
                t.row({me.name, num(win), num(tie), num(loss)});
            }
            t.print();

            // Patch count is not a quality metric in either direction -- fewer
            // patches is a simpler layout and more patches can be the layout the
            // geometry actually needs -- so it is reported and not scored.
            int pa = 0, pb = 0, npair = 0;
            for (const auto &[model, oa] : a) {
                auto it = b.find(model);
                if (it == b.end() || !meshed(oa) || !meshed(it->second)) continue;
                pa += oa->layout.patches;
                pb += it->second->layout.patches;
                ++npair;
            }
            if (npair)
                std::cout << "\n  Patches over the same " << npair << " models: " << A << " "
                          << pa << ", " << B << " " << pb
                          << " (reported, not scored -- a layout with more patches is a more"
                          << " detailed\n  decomposition, not a worse one).\n";

            std::cout <<
                "\n  How to read this table before writing Sec. 6. It scores the *final\n"
                "  element quality*, and that is not where a face-native field is claimed to\n"
                "  help. What E3 and the funnel above measure -- whether the pipeline reaches\n"
                "  a mesh at all, whether the field's topological content survives to the\n"
                "  quantizer, whether Poincare-Hopf holds -- is the claim. If this scoreboard\n"
                "  does not come out in our favour, the sentence in Sec. 1 has to be about\n"
                "  robustness and the absence of a conversion step, not about better\n"
                "  elements: a claim of 'equal-or-better layouts' that this table contradicts\n"
                "  is the kind a reviewer checks first.\n";
        }
    }

    banner("E4 summary");
    std::cout << (v.failures ? kFail : kPass) << " " << v.failures << " failure(s), "
              << v.warnings << " warning(s). CSVs in " << outDir << "/\n";
    return v.failures ? 1 : 0;
}
