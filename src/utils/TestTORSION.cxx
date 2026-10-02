// Utility to run TORSION -- Pipeline B of docs/cf_flow_pipeline.md -- on a
// mesh: the same Stages 1 to 10 as MERIDIAN, with the initial map psi_0 built
// by integrating the DualMBO cross field instead of by unfolding a Ricci metric.
//
//   TestTORSION <mesh.obj> [options]
//   TestTORSION --selftest
//   TestTORSION --ref-test <mesh.obj>
//
// Reports each stage and exits non-zero if the pipeline did not reach a valid
// layout. The two other modes are the plan's own validation steps, in the
// order it puts them:
//
//   --selftest   the four new pieces on inputs whose answers are known in
//                closed form -- a square whose cross field is constant, where
//                the integration has an exact affine solution, and where
//                Reference::Field must report det J = 1 on every triangle.
//                Sec. 9 step 1 and Sec. 13's first six checks.
//
//   --ref-test   Sec. 4's standalone test, and Sec. 10's step 1: run *Pipeline
//                A* unchanged and swap only E1's reference, so that Stage 6
//                starts from psi_R and measures against the field. C4 is
//                answered if maxConeAngleResidual does not regress. This
//                isolates the reference metric from every other source of
//                error before the integration exists to be blamed for it.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "MERIDIAN/MERIDIAN.hxx"
#include "MERIDIAN/QuadMesh.hxx"
#include "MERIDIAN/SplineFit.hxx"
#include "dualmbo/DualMBO.hxx"
#include "TORSION/ConeMetric.hxx"
#include "TORSION/MaterialLayout.hxx"
#include "TORSION/TORSION.hxx"
#include "TestHelper.hxx"
#include "mesh/QuadMesh.hxx"
#include "mesh/TMOP.hxx"
#include "mesh/Pillow.hxx"

namespace {

const char *kPass = "\033[32m[PASS]\033[0m";
const char *kFail = "\033[31m[FAIL]\033[0m";
const char *kWarn = "\033[33m[WARN]\033[0m";

int failures = 0;

void verdict(bool ok, const std::string &what) {
    std::cout << "  " << (ok ? kPass : kFail) << " " << what << "\n";
    if (!ok) ++failures;
}

void check(bool ok, const std::string &what) { verdict(ok, what); }

void heading(const std::string &title) {
    std::cout << "\n" << title << "\n";
    std::cout << std::string(title.size(), '-') << "\n";
}

void warn(const std::string &what) {
    std::cout << "  " << kWarn << " " << what << "\n";
}

// ---------------------------------------------------------------------------
// selfTest()
//
// The four pieces this pipeline adds, on inputs whose answers are known before
// the code runs.
//
// The square is the whole point of the first case. Its boundary is axis
// aligned, so the DualMBO field's Dirichlet data is the same cross on every
// boundary triangle and the smoothest field satisfying it is *constant*; a
// constant cross field is integrable exactly, its four corners take +1 each so
// sum I = 4 = 4 chi with no interior cone at all, and G is therefore empty and
// Omega is S. Every quantity downstream then has a closed form:
//
//     every matching is 0 and the comb is the identity
//     the integration's answer is the affine map x -> x/h, exactly
//     its residual against the field is 0 and nothing is inverted
//     psi_0 is an isometry of the field metric, so Immersion's metric residual
//         -- the one the field route stops treating as a verdict -- is 0 here,
//         because on this model the field really is integrable
//     Reference::Field, handed that map, reports det J = 1 on every triangle
//
// A failure in any of them is a failure in the new code and not in the model,
// which is what makes this worth having in front of the corpus.
// ---------------------------------------------------------------------------
int selfTest() {
    std::cout << "TORSION self-test (Pipeline B, docs/cf_flow_pipeline.md)\n";

    // -----------------------------------------------------------------
    heading("Case 1  Square, constant cross field: the integration is exact");
    // -----------------------------------------------------------------
    {
        std::shared_ptr<Mesh> mesh = TestHelper::createBox(0.0, 0.0, 1.0, 1.0, 0.05);
        const double h = 0.25;

        TORSION::Options opts;
        opts.targetEdge = h;
        opts.runLayout = false;
        TORSION pipe(mesh, opts);
        const bool ran = pipe.run();
        const TORSION::Status &st = pipe.getStatus();
        (void)ran;

        std::cout << "  " << mesh->vertices.size() << " vertices, "
                  << mesh->triangles.size() << " triangles; h = " << h << "\n";

        check(pipe.hasFrames(), "Stage 3F ran");
        if (!pipe.hasFrames()) return failures;

        const FieldFrames &ff = pipe.getFrames();
        const FieldFrames::Report &fr = ff.getReport();
        std::cout << "  Comb: " << fr.combedFaces << " of " << fr.faces
                  << " faces from seed " << fr.seedFace
                  << "; worst frame jump " << std::scientific << std::setprecision(3)
                  << fr.maxFrameJump << std::defaultfloat << " rad\n";

        int nonZeroMatchings = 0;
        for (int p : ff.matchings()) if (p != 0) ++nonZeroMatchings;
        std::cout << "  Matchings: " << nonZeroMatchings << " of "
                  << ff.matchings().size() << " non-zero\n";

        check(fr.unreachedFaces == 0, "The comb reached every face of Omega");
        check(fr.combingDefects == 0, "Every loop of Omega closes: a_g - a_f = -p_fg");
        check(st.secondSeedAgrees, "Combing from a second seed gives the same branch");
        check(nonZeroMatchings == 0, "A constant field has every matching at zero");
        check(fr.maxFrameJump < 1e-6, "A constant field turns by nothing between faces");
        check(fr.indexMismatches == 0, "Sec. 5.1: the index from the matchings is the field's");
        check(fr.admissible, "Sec. 5.1: sum I(v) = 4 chi(S)");
        check(fr.leftHandedFrames == 0, "det J*_t > 0 on every triangle");
        check(st.interiorCones == 0, "A square takes no interior cone");

        check(pipe.hasIntegration(), "Stage 4F ran");
        if (!pipe.hasIntegration()) return failures;
        const FieldIntegration::Report &ir = pipe.getIntegration().getReport();
        std::cout << "  Integration: " << ir.constraintRows << " constraint row(s) over "
                  << ir.seamPairs << " seam pair(s); fit residual " << std::scientific
                  << std::setprecision(3) << ir.maxFitResidual << " worst, "
                  << ir.meanFitResidual << " mean" << std::defaultfloat << "\n";
        check(ir.solved, "The constrained least-squares system solved");
        check(ir.solvedWithLDLT, "The saddle system was quasi-definite (no LU fallback)");
        check(ir.maxFitResidual < 1e-9,
              "grad u = X and grad v = Y exactly: a constant field is integrable");
        check(ir.flippedFaces == 0, "Sec. 7.1: nothing inverted");
        check(std::fabs(ir.minAreaRatio - 1.0) < 1e-9,
              "Every image triangle has exactly the area det J* asks for");

        // The map is an isometry of the model scaled by 1/h. Read straight off
        // the edges rather than off any report, so that the report and the
        // measurement are independent of each other.
        const Mesh &omega = pipe.getCut().getCutMesh();
        const std::vector<Point> &psi = pipe.getIntegratedMap();
        double worstEdge = 0.0;
        for (const auto &e : omega.edges) {
            const double want = normP(omega.vertices[e[1]] - omega.vertices[e[0]]) / h;
            const double got = normP(psi[e[1]] - psi[e[0]]);
            if (want > 0.0) worstEdge = std::max(worstEdge, std::fabs(got - want) / want);
        }
        std::cout << "  Image edge lengths off |e|/h by at most " << std::scientific
                  << std::setprecision(3) << worstEdge << std::defaultfloat << " relative\n";
        check(worstEdge < 1e-9, "psi_0 is an exact isometry of the field metric here");

        check(pipe.hasImmersion(), "Stage 4 accepted psi_0");
        if (pipe.hasImmersion()) {
            const Immersion::Report &imr = pipe.getImmersion().getReport();
            check(imr.flippedFaces == 0, "Q1: no face of psi_0 is inverted");
            check(imr.frameKConflicts == 0, "Every arc's quarter turn is read consistently");
            check(imr.maxMetricResidual < 1e-9,
                  "Immersion agrees: the field metric is realised on this model");
        }

        // ---- Tutte, on the same Omega --------------------------------
        TutteEmbedding tutte(omega, TutteEmbedding::Options());
        const TutteEmbedding::Report &tr = tutte.getReport();
        std::cout << "  Tutte: " << tr.boundaryVertices << " boundary + "
                  << tr.interiorVertices << " interior vertices onto a circle of radius "
                  << std::fixed << std::setprecision(4) << tr.radius << std::defaultfloat
                  << "; smallest signed area " << std::scientific << std::setprecision(3)
                  << tr.minSignedArea << std::defaultfloat << "\n";
        check(tr.solved, "The Tutte system solved");
        check(tr.flippedFaces == 0, "Tutte's theorem: the embedding is bijective");
        check(tr.minSignedArea > 0.0, "Every Tutte triangle is positively oriented");

        // ---- Reference::Field, on the map it was designed for ---------
        //
        // refM = (J* E)^-1, so a map whose Jacobian is J* has J = I and det J = 1
        // on every triangle. That identity is the whole of what Sec. 4 asks for,
        // and getting it wrong -- composing the two the other way round, or
        // building E from the canonical side-length triangle instead of the real
        // one -- leaves a J that is a rotation of the identity on every face and
        // an E1 that is minimised somewhere else entirely.
        SubdomainLabels::Options lopts;
        lopts.seedTopoConstraints = false;
        lopts.interfaceCorners = false;
        SubdomainLabels labels(pipe.getImmersion(), lopts);

        LayoutEnergy::Options eopts;
        eopts.reference = LayoutEnergy::Reference::Field;
        eopts.fieldFrames = ff.frames();
        LayoutEnergy le(pipe.getImmersion(), labels, eopts);
        const std::vector<double> dets = le.determinants();
        double worstDet = 0.0;
        for (double d : dets) worstDet = std::max(worstDet, std::fabs(d - 1.0));
        std::cout << "  Reference::Field on psi_0: det J is 1 to " << std::scientific
                  << std::setprecision(3) << worstDet << std::defaultfloat
                  << " over " << dets.size() << " triangles\n";
        check(le.getReport().leftHandedFrames == 0, "No left-handed frame reached the reference");
        check(worstDet < 1e-9, "Sec. 4: J = I when the map's Jacobian is J*");

        const double gradErr = le.checkGradient(16);
        std::cout << "  Assembled gradient vs. central difference: " << std::scientific
                  << std::setprecision(3) << gradErr << std::defaultfloat << " relative\n";
        check(gradErr < 1e-4, "The energy's gradient is right against the field reference");
    }

    // -----------------------------------------------------------------
    heading("Case 2  Disk: cones, a cutting graph, and a non-trivial branch");
    // -----------------------------------------------------------------
    {
        // A circle has no axis-aligned boundary anywhere, so the field is not
        // constant, the cut is not empty and the branch has something to
        // select. What is still known in advance is the arithmetic: chi = 1, so
        // sum I = 4, every arc carries one whole number of quarter turns, and
        // the branch closes around every loop of Omega.
        std::shared_ptr<Mesh> mesh = TestHelper::createCircle(0.0, 0.0, 1.0, 0.05);

        TORSION::Options opts;
        opts.runLayout = false;
        TORSION pipe(mesh, opts);
        pipe.run();
        const TORSION::Status &st = pipe.getStatus();

        std::cout << "  " << mesh->vertices.size() << " vertices, "
                  << mesh->triangles.size() << " triangles; "
                  << st.interiorCones << " interior + " << st.boundaryCones
                  << " boundary cone(s)\n";

        check(pipe.hasFrames(), "Stage 3F ran");
        if (!pipe.hasFrames()) return failures;
        const FieldFrames::Report &fr = pipe.getFrames().getReport();
        std::cout << "  Comb: " << fr.combedFaces << " of " << fr.faces
                  << " faces; worst frame jump " << std::scientific << std::setprecision(3)
                  << fr.maxFrameJump << std::defaultfloat << " rad; "
                  << fr.indexMismatches << " index mismatch(es)\n";
        check(fr.unreachedFaces == 0, "The comb reached every face of Omega");
        check(fr.combingDefects == 0, "Every loop of Omega closes");
        check(st.secondSeedAgrees, "Combing from a second seed gives the same branch");
        check(fr.admissible, "sum I(v) = 4 chi(S)");
        check(fr.indexMismatches == 0, "The index from the matchings is the field's own");

        if (pipe.hasIntegration()) {
            const FieldIntegration::Report &ir = pipe.getIntegration().getReport();
            std::cout << "  Integration: seam residual " << std::scientific
                      << std::setprecision(3) << ir.maxSeamResidual
                      << " of the extent; fit residual " << ir.maxFitResidual << " worst, "
                      << ir.meanFitResidual << " mean" << std::defaultfloat << "; "
                      << ir.flippedFaces << " inverted face(s)\n";
            check(ir.solved, "The constrained least-squares system solved");
            check(ir.maxSeamResidual < 1e-9,
                  "C2: the seam constraints hold to rounding, not to a tolerance");
        }
        if (pipe.hasImmersion()) {
            const Immersion::Report &imr = pipe.getImmersion().getReport();
            std::cout << "  psi_0: " << imr.arcs << " arc(s), " << imr.seamEdgePairs
                      << " seam pair(s); Gamma_Hol_0..3 = " << imr.holonomyCount[0] << ", "
                      << imr.holonomyCount[1] << ", " << imr.holonomyCount[2] << ", "
                      << imr.holonomyCount[3] << "\n";
            check(imr.frameKConflicts == 0,
                  "Every seam edge of an arc agrees with it about the quarter turn");
            check(imr.flippedFaces == 0, "Q1: psi_0 is locally injective");
        } else {
            warn("Stage 4 did not accept psi_0 on the disk; see the messages above.");
        }
        for (const std::string &m : st.messages) {
            if (m.find("Stage 4") != std::string::npos) std::cout << "     " << m << "\n";
        }
    }

    // -----------------------------------------------------------------
    heading("Case 3  The cone metric, where its answer is known before it runs");
    // -----------------------------------------------------------------
    //
    // Sec. 4's Newton solve, on two models whose conformal factor is known in
    // closed form.
    //
    // The square is the sharper of the two: its four corners are the four
    // cones, each carries I = +1, and pi - (pi/2)(+1) = pi/2 is *already* the
    // angle the model has there. Every other vertex wants what it already has
    // too, so K = Kbar everywhere before a step is taken, the metric is the
    // model's own, and u is constant. A solve that moves anything here is
    // moving it for a reason that does not exist.
    //
    // The disk has no such luck -- its boundary carries 2 pi of turning that
    // has to go into the cones -- so what is known there is not u but what u
    // must achieve: K = Kbar at every vertex, to rounding, and a metric that is
    // still a metric when it gets there.
    {
        struct Case {
            const char *what;
            std::shared_ptr<Mesh> mesh;
            bool flat;      // is u constant?
        };
        std::vector<Case> cases = {
            {"square", TestHelper::createBox(0.0, 0.0, 1.0, 1.0, 0.05), true},
            {"disk", TestHelper::createCircle(0.0, 0.0, 1.0, 0.05), false},
        };
        for (Case &c : cases) {
            DualMBO field(c.mesh, 500, 10.0);
            field.initialize();
            const double nTris = static_cast<double>(c.mesh->triangles.size());
            for (int i = 0; i < 500; ++i) {
                field.step();
                if (field.error < 2.0 * nTris * 1e-5) break;
            }
            field.computeSingularities();
            ConeSingularities cones(field);
            ConeMetric metric(*c.mesh, cones);
            const ConeMetric::Report &r = metric.getReport();
            std::cout << "  " << c.what << ": ||K - Kbar||_inf " << std::scientific
                      << std::setprecision(3) << r.initialError << " -> " << r.finalError
                      << " rad in " << r.newtonIterations << " Newton step(s); exp(u) in ["
                      << std::fixed << std::setprecision(6) << r.minScale << ", " << r.maxScale
                      << "]" << std::defaultfloat << "\n";
            check(r.solved, std::string(c.what) + ": every face of the metric is realisable");
            check(r.finalError < 1e-6,
                  std::string(c.what) + ": K = Kbar at every vertex");
            check(r.coneResidual < 1e-6,
                  std::string(c.what) + ": C4, the angle sum at every cone is 2 pi - (pi/2) I");
            if (c.flat) {
                check(std::fabs(r.maxScale / r.minScale - 1.0) < 1e-9,
                      "square: u is constant, because the model already is the flat cone metric");
            } else {
                check(r.maxScale / r.minScale > 1.001,
                      "disk: u is not constant, because 2 pi of boundary turning had to move");
            }
        }
    }

    heading("Self-test result");
    if (failures == 0) {
        std::cout << "  " << kPass << " Every closed-form check held.\n";
    } else {
        std::cout << "  " << kFail << " " << failures << " check(s) failed.\n";
    }
    return failures == 0 ? 0 : 1;
}

// ---------------------------------------------------------------------------
// referenceTest()  --  Sec. 4's standalone test, Sec. 10's step 1
//
// Pipeline A, unchanged, four times: same psi_R, same labelling, same schedule,
// one difference, which is what E1 measures J against. Sec. 4 is explicit about
// what this is for -- "that isolates C4 from every other source of error before
// the integration exists to be blamed" -- and about the number to read:
// maxConeAngleResidual out of Stage 6 must not regress against Reference::Ricci.
//
// The four:
//
//   Ricci      Pipeline A's flat cone metric. The yardstick.
//   Euclidean  the model's own geometry. LayoutEnergy.hxx measures what this
//              does on geom003 and it is the failure C4 is a warning about.
//   Field      Sec. 4's construction: the Euclidean reference composed with the
//              field's per-triangle frame.
//   Induced    the metric psi_R itself induces, which is what Pipeline B uses
//              (there, from psi_0) and is here to show the mechanism works from
//              either map.
//
// **What it finds.** Field and Euclidean are the same reference. The frame is a
// scaled rotation, so composing with it multiplies J on the right by a
// rotation, and neither ||J||_F nor ||J^-1||_F can see one; at h = 1 the two
// energies are identical term for term. The first verdict below asserts that
// identity directly on det J, because it is the kind of fact that is easy to
// state and easy to stop being true, and because Sec. 4's whole argument turns
// on the frame carrying the cones into E1 -- which it does not.
// ---------------------------------------------------------------------------
int referenceTest(const std::string &path) {
    std::shared_ptr<Mesh> mesh;
    try {
        mesh = std::make_shared<Mesh>(path);
    } catch (const std::exception &e) {
        std::cout << kFail << " Failed to load mesh: " << e.what() << "\n";
        return 2;
    }

    std::cout << "Reference metrics for E1, standalone on psi_R (Sec. 4, Sec. 10 step 1)\n";
    std::cout << "Mesh: " << path << "\n  " << mesh->vertices.size() << " vertices, "
              << mesh->triangles.size() << " triangles\n";

    // Stages 0 to 4 once. Every run below starts from this psi_R, so the only
    // thing that differs between them is the reference.
    MERIDIAN::Options opts;
    MERIDIAN::Options stopAt4 = opts;
    stopAt4.runLayout = false;
    MERIDIAN front(mesh, stopAt4);
    if (!front.run() || !front.hasImmersion()) {
        std::cout << "  " << kFail << " Pipeline A did not reach psi_R on this model.\n";
        return 4;
    }
    FieldFrames frames(front.getField(), front.getCut(), front.getCones(),
                       FieldFrames::Options());

    // psi_R's own induced metric, per edge of S. The same construction TORSION
    // uses on psi_0; see TORSION::run().
    std::vector<double> induced(mesh->edges.size(), 0.0);
    {
        const Mesh &om = front.getCut().getCutMesh();
        const auto &c2o = front.getCut().getCutVertexToOriginal();
        const std::vector<Point> &psi = front.getImmersion().getUV();
        std::unordered_map<MeshEdgeKey, int, MeshEdgeKeyHash> idx;
        for (int e = 0; e < static_cast<int>(mesh->edges.size()); ++e) {
            idx.emplace(MeshEdgeKey(mesh->edges[e][0], mesh->edges[e][1]), e);
        }
        for (int e = 0; e < static_cast<int>(om.edges.size()); ++e) {
            auto it = idx.find(MeshEdgeKey(c2o[om.edges[e][0]], c2o[om.edges[e][1]]));
            if (it == idx.end() || induced[it->second] > 0.0) continue;
            induced[it->second] = normP(psi[om.edges[e][1]] - psi[om.edges[e][0]]);
        }
        for (size_t e = 0; e < induced.size(); ++e) {
            if (!(induced[e] > 0.0)) {
                induced[e] = normP(mesh->vertices[mesh->edges[e][1]] -
                                   mesh->vertices[mesh->edges[e][0]]);
            }
        }
    }

    // The angle sum the reference metric itself puts at each cone, against the
    // 2 pi - (pi/2) I the layout needs. This is C4 stated as a number, on the
    // reference alone, before any continuation is allowed to use it.
    auto coneAnglesOf = [&](const std::vector<double> &len) {
        std::vector<double> sum(mesh->vertices.size(), 0.0);
        for (int t = 0; t < static_cast<int>(mesh->triangles.size()); ++t) {
            const double l[3] = {len[mesh->triangleEdges[t][0]],
                                 len[mesh->triangleEdges[t][1]],
                                 len[mesh->triangleEdges[t][2]]};
            for (int i = 0; i < 3; ++i) {
                const double p = l[i], q = l[(i + 2) % 3], r = l[(i + 1) % 3];
                if (!(p > 0.0) || !(q > 0.0)) continue;
                double c = (p * p + q * q - r * r) / (2.0 * p * q);
                c = std::max(-1.0, std::min(1.0, c));
                sum[mesh->triangles[t][i]] += std::acos(c);
            }
        }
        const std::vector<int> &I = front.getCones().getIndices();
        double worst = 0.0;
        for (const auto &c : front.getCones().getCones()) {
            const double full = mesh->isBoundaryVertex[c.vertex] ? M_PI : 2.0 * M_PI;
            worst = std::max(worst, std::fabs(sum[c.vertex] - (full - M_PI_2 * I[c.vertex])));
        }
        return worst;
    };
    std::vector<double> euclidLen(mesh->edges.size(), 0.0);
    for (size_t e = 0; e < euclidLen.size(); ++e) {
        euclidLen[e] = normP(mesh->vertices[mesh->edges[e][1]] -
                             mesh->vertices[mesh->edges[e][0]]);
    }

    struct Run {
        const char *name;
        LayoutEnergy::Reference ref;
        double refConeAngle = 0.0;
        double coneResidual = 0.0, regularResidual = 0.0;
        double minDetJ = 0.0, boundary = 0.0, seam = 0.0, topo = 0.0;
        int coneValenceChanges = 0, inverted = 0, outer = 0;
        double minScaledJacobian = 0.0;
        int invertedQuads = 0, meshQuads = 0, patches = 0, quadPatches = 0;
        std::vector<double> dets;
    };
    Run runs[4] = {{"Ricci", LayoutEnergy::Reference::Ricci},
                   {"Euclidean", LayoutEnergy::Reference::Euclidean},
                   {"Field", LayoutEnergy::Reference::Field},
                   {"Induced", LayoutEnergy::Reference::Induced}};

    for (Run &r : runs) {
        SubdomainLabels::Options lopts;
        lopts.nearMissTolerance = opts.topoNearMiss;
        lopts.maxTraceSteps = opts.separatrixMaxSteps;
        SubdomainLabels labels(front.getImmersion(), lopts,
                               front.hasInterfaces() ? &front.getInterfaces() : nullptr);

        LayoutEnergy::Options eopts;
        eopts.lambdaInit = opts.lambdaInit;
        eopts.lambdaGrowth = opts.lambdaGrowth;
        eopts.outerSteps = opts.outerSteps;
        eopts.innerIterations = opts.innerIterations;
        eopts.reference = r.ref;
        if (r.ref == LayoutEnergy::Reference::Field) eopts.fieldFrames = frames.frames();
        if (r.ref == LayoutEnergy::Reference::Induced) eopts.referenceLengths = induced;

        LayoutEnergy le(front.getImmersion(), labels, eopts);
        // Before any minimisation: what the reference is, per triangle.
        r.dets = le.determinants();
        switch (r.ref) {
            case LayoutEnergy::Reference::Ricci:
                r.refConeAngle = coneAnglesOf(front.getImmersion().getFlatEdgeLengths()); break;
            case LayoutEnergy::Reference::Induced:
                r.refConeAngle = coneAnglesOf(induced); break;
            default:
                r.refConeAngle = coneAnglesOf(euclidLen); break;
        }

        le.run();
        Separatrices::Options sopts;
        sopts.coneSnapTolerance = opts.separatrixSnap;
        sopts.maxSteps = opts.separatrixMaxSteps;
        sopts.coneSnapRings = opts.separatrixSnapRings;
        sopts.nearMissWindow = opts.repairGapLimit;
        MERIDIAN::RepairOptions ropts;
        ropts.passes = opts.repairPasses;
        ropts.maxPerPass = opts.repairMaxPerPass;
        ropts.gapLimit = opts.repairGapLimit;
        ropts.lambdaBoost = opts.repairLambdaBoost;
        ropts.outerSteps = opts.repairOuterSteps;
        ropts.patience = opts.repairPatience;
        MERIDIAN::RepairResult rep =
            MERIDIAN::traceAndRepair(labels, le, sopts, ropts, [](const std::string &) {});

        const LayoutEnergy::Report &er = le.getReport();
        r.coneResidual = er.maxConeAngleResidual;
        r.regularResidual = er.maxRegularAngleResidual;
        r.minDetJ = er.minDetJ;
        r.boundary = er.maxBoundaryResidual;
        r.seam = er.maxSeamResidual;
        r.topo = er.maxTopoResidual;
        r.coneValenceChanges = er.coneValenceChanges;
        r.inverted = er.invertedTriangles;
        r.outer = er.outerSteps;

        if (rep.separatrices) {
            try {
                Arrangement arr(*rep.separatrices, labels, ropts.arrangement);
                r.patches = arr.getReport().patches;
                r.quadPatches = arr.getReport().simpleQuads;
                SplineFit::Options sfo;
                sfo.segments = opts.splineSegments;
                sfo.samples = opts.splineSamples;
                SplineFit fit(arr, sfo);
                QuadMesh::Options qo;
                qo.targetEdgeLength = opts.quadTargetEdge;
                qo.collapseSpan = opts.quadCollapseSpan;
                qo.smoothingPasses = opts.quadSmoothingPasses;
                QuadMesh qm(fit, qo);
                r.minScaledJacobian = qm.getReport().minScaledJacobian;
                r.invertedQuads = qm.getReport().invertedQuads;
                r.meshQuads = qm.getReport().quads;
            } catch (const std::exception &) {
                r.meshQuads = -1;
            }
        }
    }

    heading("The same psi_R, four references for E1");
    std::cout << "  " << std::left << std::setw(20) << "";
    for (const Run &r : runs) std::cout << std::right << std::setw(13) << r.name;
    std::cout << "\n";
    auto row = [&](const char *what, double Run::*m) {
        std::cout << "  " << std::left << std::setw(20) << what << std::right
                  << std::scientific << std::setprecision(3);
        for (const Run &r : runs) std::cout << std::setw(13) << r.*m;
        std::cout << std::defaultfloat << "\n";
    };
    auto irow = [&](const char *what, int Run::*m) {
        std::cout << "  " << std::left << std::setw(20) << what << std::right;
        for (const Run &r : runs) std::cout << std::setw(13) << r.*m;
        std::cout << "\n";
    };
    row("C4 ref cone angle", &Run::refConeAngle);
    row("Q2 at the cones", &Run::coneResidual);
    row("Q2 elsewhere", &Run::regularResidual);
    row("Q3 boundary", &Run::boundary);
    row("Q4 seam", &Run::seam);
    row("Q5 topo", &Run::topo);
    row("min det J", &Run::minDetJ);
    irow("cone valence moves", &Run::coneValenceChanges);
    irow("inverted triangles", &Run::inverted);
    irow("outer steps", &Run::outer);
    irow("patches", &Run::patches);
    irow("four-sided patches", &Run::quadPatches);
    row("Stage 10 min sJ", &Run::minScaledJacobian);
    irow("inverted quads", &Run::invertedQuads);
    irow("quads", &Run::meshQuads);

    // Field against Euclidean, on the reference itself rather than on anything
    // the continuation did with it.
    double fieldVsEuclid = 0.0;
    for (size_t i = 0; i < runs[1].dets.size() && i < runs[2].dets.size(); ++i) {
        const double d = std::fabs(runs[1].dets[i]);
        fieldVsEuclid = std::max(fieldVsEuclid,
                                 std::fabs(runs[1].dets[i] - runs[2].dets[i]) /
                                     std::max(d, 1e-30));
    }
    double ricciVsEuclid = 0.0;
    for (size_t i = 0; i < runs[0].dets.size() && i < runs[1].dets.size(); ++i) {
        const double d = std::fabs(runs[1].dets[i]);
        ricciVsEuclid = std::max(ricciVsEuclid,
                                 std::fabs(runs[0].dets[i] - runs[1].dets[i]) /
                                     std::max(d, 1e-30));
    }

    heading("Verdict");
    std::cout << "  det J at psi_R: Field vs Euclidean differ by " << std::scientific
              << std::setprecision(3) << fieldVsEuclid << ", Ricci vs Euclidean by "
              << ricciVsEuclid << std::defaultfloat << " relative\n";
    // The finding, asserted rather than described: Sec. 4's reference is the
    // Euclidean one. J_field = h J_true R(theta), a rotation on the right, and
    // the symmetric Dirichlet energy cannot see one.
    verdict(fieldVsEuclid < 1e-12,
            "Reference::Field is Reference::Euclidean: the frame is a scaled rotation and "
            "E1 is blind to a rotation on the right of J");
    verdict(ricciVsEuclid > 1e-3,
            "Reference::Ricci is genuinely a different reference (the control for the above)");
    verdict(runs[1].refConeAngle > 1.0,
            "C4 restated: the Euclidean reference has no cones at all");
    verdict(runs[3].refConeAngle < 1e-3,
            "C4 answered: the induced reference has the cones Q2 asks for");
    verdict(runs[3].coneResidual <= std::max(10.0 * runs[0].coneResidual, 1e-9),
            "Q2 out of Stage 6 does not regress against Reference::Ricci");
    verdict(runs[3].inverted == 0, "Q1 survives the reference switch");
    // Reported, not asserted. Sec. 10 names this number on geom003 and it is
    // the right one to watch, but what moves it is the *boundary*, not the
    // cones: a reference built from a map whose boundary is not geodesic keeps
    // that boundary's corners, and Stage 10 meets them as reflex patch corners.
    // That is a property of the map the reference was taken from, and saying so
    // is more use than failing a run over it.
    if (runs[3].invertedQuads > runs[0].invertedQuads) {
        warn("Stage 10 inverts " + std::to_string(runs[3].invertedQuads) +
             " element(s) under the induced reference against " +
             std::to_string(runs[0].invertedQuads) +
             " under Ricci. The induced metric keeps whatever boundary corners the map it "
             "was taken from had, and only a stage that makes dS geodesic removes them.");
    }
    return failures == 0 ? 0 : 1;
}

void usage(const char *prog) {
    std::cout << "Usage: " << prog << " <mesh.obj> [options]\n"
              << "       " << prog << " --selftest\n"
              << "       " << prog << " --ref-test <mesh.obj>\n\n"
              << "Stage 3F / 4F / 4R (the field route)\n"
              << "  --h <len>          target edge length h of the frame     (default 1)\n"
              << "  --no-untangle      leave psi_0 as the integration left it, flips and all\n"
              << "  --untangle-outer <n>  outer steps of the target fit      (default 8)\n"
              << "  --untangle-mu <m>  starting weight on E4 in that fit     (default 1e-2)\n"
              << "  --reg <e>          saddle-block regularisation           (default 1e-10)\n"
              << "  --no-seed-check    skip the second-seed spot check\n"
              << "  --no-align         leave Q3 and E3 entirely to Stage 6 (Sec. 6.4 off)\n"
              << "  --no-align-fallback  keep the aligned solve even when it inverts more\n"
              << "  --no-seam-align    vote an interface branch the cut crosses as one chart\n"
              << "  --choose-by-flips  pick among the alignment attempts by raw flip count,\n"
              << "                     not by what Sec. 7.2a's ladder leaves of each\n"
              << "  --no-regularised-untangle  Sec. 7.2a's kernel pass alone, no Sec. 7.2c\n"
              << "  --no-pull          hand Stage 5 psi_0 unaligned when Stage 4F dropped the\n"
              << "                     alignment, instead of pulling it back on (Sec. 7.2b)\n"
              << "  --release-rounds <n>  rounds of releasing only the chains a tangle touches\n"
              << "                     before dropping the alignment              (default 4)\n"
              << "  --no-sector-reconcile  move the field's indices vertex by vertex, dS\n"
              << "                     absorbing anything, and re-smooth holding one face\n"
              << "  --no-conformal-sizing  leave the frame unscaled: J* = (1/h) R(-theta)\n"
              << "  --ref <r>          E1's reference metric, C4:                (default cone)\n"
              << "                       cone      the flat cone metric of the cone set, from\n"
              << "                                 ConeMetric's Newton solve on the model's\n"
              << "                                 own triangles. The answer to C4\n"
              << "                       induced   the metric psi_0 induces on Omega, which\n"
              << "                                 has the cone angles the seam constraints\n"
              << "                                 already forced into it\n"
              << "                       field     Sec. 4's frame construction -- measured,\n"
              << "                                 and identical to euclidean, see --ref-test\n"
              << "                       euclidean the model's own geometry, which on a planar\n"
              << "                                 input has no cones at all\n"
              << "                       ricci     Pipeline A's flat cone metric, computed for\n"
              << "                                 the reference alone; the known-good yardstick\n\n"
              << "Per-material mode (a multi-material model; docs/cf_flow_pipeline.md Sec. 15)\n"
              << "  --whole-model      lay the whole model out at once rather than one material\n"
              << "                     region at a time (Options::perMaterial off)\n"
              << "  --pm-rounds <n>    rounds of matching at most                  (default 10)\n"
              << "  --pm-threads <n>   regions laid out at once, 0 = one per hardware thread\n"
              << "  --pm-tolerance <e> how many edges apart two layout vertices across an\n"
              << "                     interface may be and still be matched       (default 3)\n"
              << "  --pm-no-fallback   keep the glued layout even when it is not valid, rather\n"
              << "                     than running the whole model as well\n"
              << "  --pm-cancel-always keep a region's +1/-1 cancellation even when it leaves\n"
              << "                     faces that are not simple quadrilaterals\n"
              << "  --pm-most-pairs    match by most pairs alone, wherever the glued nodes land\n"
              << "  --pm-move-ends     move a matched curve's last point onto its glued node\n"
              << "                     instead of bending the curve into it\n\n"
              << "Stages 5 to 10 (MERIDIAN's, unchanged)\n"
              << "  --gamma <g>        edge penalty                          (default 10)\n"
              << "  --steps <n>        dual-mesh MBO steps                        (default 500)\n"
              << "  --weight <w>       p=0 penalty weight: min | harm | orth   (default orth)\n"
              << "                       At p=0 the weight is the discrete Laplacian itself,\n"
              << "                       not a stabilisation parameter; orth is the two-point\n"
              << "                       (finite-volume) weight and the consistent one\n"
              << "  --no-continuation  one tau = D^2/10 instead of the annealed ladder\n"
              << "  --tau-ratio <r>    ratio between continuation levels      (default 0.25)\n"
              << "  --tau-floor <c>    stop once the diffusion length reaches c mean edges\n"
              << "                                                            (default 20)\n"
              << "  --cut-to-graph     let a cone arc stop on an earlier arc\n"
              << "  --no-interfaces    ignore the material tags\n"
              << "  --no-field-interfaces  do not align the field to the interfaces\n"
              << "  --no-cancel-dipoles    keep the +1/-1 pairs on curved interfaces\n"
              << "  --no-dipole-retry  do not rerun Stages 1-4R keeping the pairs when\n"
              << "                     cancelling them cost Stage 4F the alignment\n"
              << "  --keep-flat-cones  leave a +1 cone where Stage 1 put it on straight dS\n"
              << "                     (a patch corner of pi; Stage 1 before 2026-10-02)\n"
              << "  --no-e6            leave the interfaces to E3 alone\n"
              << "  --no-seam-turns    give every chain of an interface branch the branch's\n"
              << "                     label, ignoring the quarter turns of the cuts it crosses\n"
              << "  --outer <n>        penalty continuation steps            (default 16)\n"
              << "  --inner <n>        inner iterations per step             (default 60)\n"
              << "  --lambda <l>       initial lambda_2..lambda_5            (default 1e-2)\n"
              << "  --growth <g>       lambda growth per outer step          (default 10)\n"
              << "  --align-factor <f> factor on lambda_2, lambda_3          (default 0.1)\n"
              << "  --seam-factor <f>  factor on lambda_4                    (default 10)\n"
              << "  --no-topo          skip the Gamma_topo seeding (E5 off)\n"
              << "  --near-miss <f>    Gamma_topo seeding tolerance          (default 0.15)\n"
              << "  --retry <f>        tolerance to re-seed from when the first layout\n"
              << "                     leaves a piece of S no grid covers      (default 0.04)\n"
              << "  --last-retry <f>   a third, tighter rung of that retry           (default 0.01)\n"
              << "  --no-unseeded-retry  stop the retry before its last rung: no Gamma_topo\n"
              << "                     seeding at all, the repair loop alone\n"
              << "  --no-retry         one attempt only, at --near-miss\n"
              << "  --no-layout        stop after Stage 4\n"
              << "  --no-trace         skip Stage 7\n"
              << "  --no-arrange       skip Stage 8\n"
              << "  --no-splines       skip Stage 9\n"
              << "  --no-mesh          skip Stage 10\n"
              << "  --target <h>       Stage 10 target edge, model units     (default 0.05)\n"
              << "  --collapse-span <f>  contract a chord whose every patch is thinner\n"
              << "                     than this times the target                 (default 0.5)\n"
              << "  --no-collapse      keep every chord, however thin its patches\n"
              << "  --sample-materials give each element the material its samples find, not\n"
              << "                     its layout face's (Stage 10 before 2026-10-01)\n"
              << "  --contract-at-root place a contracted run of nodes where the union-find\n"
              << "                     left its root, on dS or an interface or not\n"
              << "  --disk-templates   excise every circular inclusion (Stage 0c) and fill it\n"
              << "                     back in with an O-grid template (Stage 11)\n"
              << "  --disk-squareness <w>  how square the template's core is  (default 0.55)\n"
              << "  --disk-ring <n>    rows of elements in the ring, 0 = auto (default 0)\n"
              << "  --disk-smooth <n>  smoothing sweeps over the template     (default 300)\n"
              << "  --repair <n>       rounds of Sec. 3.3's repair           (default 6)\n"
              << "  --cones <n>        list at most n cones                  (default 20)\n\n"
              << "Output\n"
              << "  --psi <file.obj>   write psi_0, the map Stage 4 accepted\n"
              << "  --raw <file.obj>   write the integration's map before any untangling\n"
              << "  --tutte <file.obj> write the Tutte embedding\n"
              << "  --layout <file.obj> write Psi, the Stage 6 result\n"
              << "  --cut <file.obj>   write Omega\n"
              << "  --arr <file.obj>   write the Stage 8 arrangement, on the model\n"
              << "  --quads <file.obj> write the Stage 10 quadrilateral mesh\n"
              << "  --mfem <f.mesh>    write the final mesh for MFEM, material id per element\n"
              << "  --tmop <n>         smooth the final mesh with n TMOP sweeps (metric 007,\n"
              << "                     shape and size) before writing it, 0 = off (default 0)\n"
              << "  --tmop-gauss       sample TMOP at the 2x2 Gauss points, not the element\n"
              << "                     corners (Stage 12 before 2026-10-02)\n"
              << "  --no-pillow        TMOP straight away, without first pillowing the feature\n"
              << "                     corners an element spans flat (Stage 12 before 2026-10-02)\n";
}

} // namespace

// ---------------------------------------------------------------------------
// Stages 1 to 6 as the whole-model route reports them. Returns -1 to go on to
// Stage 7, or main()'s exit code where the pipeline stopped short.
// ---------------------------------------------------------------------------
struct WholeModelOutputs {
    std::string cutOut, rawOut, tutteOut, psiOut, layoutOut;
    int coneListLimit = 20;
};

static int reportWholeModel(const TORSION &pipeline, const std::shared_ptr<Mesh> &mesh,
                            const TORSION::Options &opts, const WholeModelOutputs &outputs) {
    const TORSION::Status &st = pipeline.getStatus();
    auto stageMessages = [&](const std::string &prefix) {
        for (const std::string &m : st.messages) {
            if (m.rfind(prefix, 0) == 0) warn(m.substr(prefix.size()));
        }
    };
    const int coneListLimit = outputs.coneListLimit;
    const std::string &cutOut = outputs.cutOut;
    const std::string &rawOut = outputs.rawOut;
    const std::string &tutteOut = outputs.tutteOut;
    const std::string &psiOut = outputs.psiOut;
    const std::string &layoutOut = outputs.layoutOut;

    // ---------------------------------------------------------------------
    heading("Stage 1  Cone singularities (Sec. 3.1)");
    std::cout << "  " << st.interiorCones << " interior + " << st.boundaryCones
              << " boundary cone(s)";
    if (st.coneDipoleUnits > 0) std::cout << "; " << st.coneDipoleUnits << " dipole unit(s) cancelled";
    if (st.dipolesKept) std::cout << "; the same-region pairs kept (see below)";
    std::cout << "\n";
    verdict(st.conesAdmissible, "sum I(v) = 4 chi(S), Eq. (4)");
    if (st.dipoleRetryRan) {
        for (const std::string &m : st.messages) {
            if (m.rfind("Stage 1: Stage 1 cancelled", 0) == 0) warn(m.substr(9));
        }
    }
    if (!st.conesAdmissible) {
        stageMessages("Stage 1: ");
        heading("Result");
        std::cout << "  " << kFail << " No seamless map with these holonomies exists.\n";
        return 4;
    }
    if (coneListLimit > 0 && pipeline.hasCones()) {
        int shown = 0;
        for (const auto &c : pipeline.getCones().getCones()) {
            if (shown++ >= coneListLimit) break;
            const Point &p = mesh->vertices[c.vertex];
            std::cout << "     v" << c.vertex << " (" << std::fixed << std::setprecision(4)
                      << p[0] << ", " << p[1] << ")" << std::defaultfloat
                      << " I = " << std::showpos << c.index << std::noshowpos
                      << (c.onBoundary ? "  on dS" : "")
                      << (c.prescribed ? "  prescribed" : "")
                      << (c.fromRebalance ? "  rebalanced" : "") << "\n";
        }
    }

    // ---------------------------------------------------------------------
    heading("Stage 2  Cutting graph and cut disk (Sec. 3.2.2)");
    const ConeCut::Report &cr = pipeline.getCut().getReport();
    std::cout << "  Omega: " << pipeline.getCut().getCutMesh().vertices.size()
              << " vertices, chi = " << cr.eulerCharacteristic << ", "
              << cr.boundaryComponents << " boundary component(s); G has "
              << pipeline.getCut().getCutEdges().size() << " edge(s)\n";
    verdict(cr.isDisk, "Omega = S - G is a topological disk");
    verdict(cr.allConesOnBoundary, "P is contained in G union dS");
    // Sec. 6.1: with every cone arc running to dS the seam graph is a forest,
    // which is what lets the elimination be one Cholesky rather than a saddle
    // system. TORSION keeps the saddle system either way, so a junction costs
    // conditioning and not correctness -- but it is still the case the option
    // was meant to avoid.
    verdict(cr.interiorJunctions == 0,
            "Sec. 6.1: the seam graph is a forest (no interior junction of G)");
    if (!cutOut.empty()) {
        if (pipeline.getCut().writeOBJ(cutOut)) std::cout << "  Wrote Omega to " << cutOut << "\n";
        else warn("Failed to write " + cutOut);
    }

    if (!pipeline.hasFrames()) {
        stageMessages("Stopping");
        heading("Result");
        std::cout << "  " << kFail << " The pipeline stopped before the combing.\n";
        return 5;
    }

    // ---------------------------------------------------------------------
    if (st.coneMetricRan) {
        heading("Stage 4C  The flat cone metric of the cone set (Sec. 4, C4)");
        std::cout << "  ||K - Kbar||_inf " << std::scientific << std::setprecision(3)
                  << st.coneMetricInitialError << " rad on the model, " << st.coneMetricLinearError
                  << " after the linear solve, " << st.coneMetricFinalError << " after "
                  << st.coneMetricNewtonSteps << " Newton step(s)" << std::defaultfloat << "\n";
        std::cout << "  exp(u) in [" << std::fixed << std::setprecision(4) << st.coneMetricMinScale
                  << ", " << st.coneMetricMaxScale << "]" << std::defaultfloat << "; "
                  << st.coneMetricNonRealisable << " face(s) fail the triangle inequality\n";
        verdict(st.coneMetricSolved, "The cone metric is a metric: every face realisable");
        verdict(st.coneMetricConeResidual < 1e-6,
                "Its angle sum at every cone is the 2 pi - (pi/2) I the layout is held to");
        if (st.conformalSizingUsed) {
            std::cout << "  The frame is scaled by it: J*_t = (exp(u_t)/h) R(-theta_t)\n";
        }
        stageMessages("Stage 4C: ");
    }

    // ---------------------------------------------------------------------
    heading("Stage 3F  Combing, matchings and the index audit (Sec. 5)");
    const FieldFrames &ff = pipeline.getFrames();
    const FieldFrames::Report &fr = ff.getReport();
    int nonZero = 0;
    for (int p : ff.matchings()) if (p != 0) ++nonZero;
    std::cout << "  Combed " << fr.combedFaces << " of " << fr.faces
              << " face(s) from seed " << fr.seedFace << "; " << nonZero
              << " non-zero matching(s) over " << ff.matchings().size() << " edge(s)\n";
    std::cout << "  Worst turn between neighbouring faces: " << std::scientific
              << std::setprecision(3) << fr.maxFrameJump << std::defaultfloat
              << " rad (this is where Sec. 7's flips will be)\n";
    std::cout << "  Field metric: h = " << std::fixed << std::setprecision(4)
              << opts.targetEdge << ", per-edge disagreement at most " << std::scientific
              << std::setprecision(3) << fr.maxMetricDisagreement << std::defaultfloat
              << " relative\n";
    verdict(fr.unreachedFaces == 0, "The comb reached every face of Omega");
    verdict(fr.combingDefects == 0, "Every loop of Omega closes: a_g - a_f = -p_fg");
    verdict(st.secondSeedAgrees, "Combing from a second seed selects the same branch");
    verdict(fr.indexMismatches == 0,
            "Sec. 5.1: I from the matchings agrees with the field at every interior vertex");
    verdict(fr.leftHandedFrames == 0, "det J*_t > 0 on every triangle");
    std::cout << "  dS: the frame turns a quarter at " << fr.boundaryTurns
              << " vertex/vertices, " << fr.boundaryTurnMismatches
              << " of them where the cone set says it should not (worst rounding "
              << std::scientific << std::setprecision(3) << fr.maxBoundaryTurnResidual
              << " rad)" << std::defaultfloat << "\n";
    // Not a verdict on the new code, and stated as a prediction rather than as
    // a finding, because that is what it is: each of these is a corner the
    // finished layout will have at a vertex Stage 1 called regular, and the
    // count here is the count the arrangement will report as extra nodes and Q2
    // as non-cone vertices off pi by a quarter turn. Reading it now, before the
    // integration, is the whole point of Sec. 5.1.
    if (fr.boundaryTurnMismatches == 0) {
        std::cout << "  " << kPass
                  << " The frame's boundary staircase has its steps at cones and nowhere else\n";
    } else {
        warn("The layout's boundary will turn at " + std::to_string(fr.boundaryTurnMismatches) +
             " vertex/vertices Stage 1 called regular. See the message below: this is the "
             "structural gap between the two routes, not a defect of the combing.");
    }
    if (fr.clusteredConePairs > 0) {
        warn(std::to_string(fr.clusteredConePairs) +
             " triangle(s) carry two or more cones (Sec. 3.1's clustering)");
    }
    if (fr.highIndexCones > 0) {
        warn(std::to_string(fr.highIndexCones) + " cone(s) of index |I| > 2");
    }
    // Sec. 6.4. The three numbers that matter are the last three: how many
    // edges Stage 4F will hold outright, how many steps of the field's own
    // staircase were taken out to do it, and how far off an axis the field was
    // on the worst of them.
    std::cout << "  Alignment: " << fr.alignmentChains << " chain(s) of dS - G and of the "
              << "interface network -> " << fr.alignedBoundaryEdges << " boundary + "
              << fr.alignedInterfaceEdges << " interface edge(s) held; "
              << fr.alignmentOverrides << " staircase step(s) removed over "
              << fr.alignmentSteppedChains << " chain(s); " << fr.alignmentSeamCrossings
              << " crossing(s) of G read across in the branch's own chart\n";
    std::cout << "  Left free: " << fr.alignmentClosedChains << " closed, "
              << fr.alignmentReversedChains << " that turn back on themselves, "
              << fr.alignmentAbstained << " the field is not on an axis of\n";
    if (fr.alignedBoundaryEdges + fr.alignedInterfaceEdges > 0) {
        std::cout << "  Worst edge held is " << std::scientific << std::setprecision(3)
                  << fr.maxAlignmentResidual << std::defaultfloat
                  << " rad off the axis the frame reads for it\n";
    }
    stageMessages("Stage 3F: ");

    if (!pipeline.hasIntegration()) {
        heading("Result");
        std::cout << "  " << kFail << " The pipeline stopped before the integration.\n";
        return 5;
    }

    // ---------------------------------------------------------------------
    heading("Stage 4F  Seamless integration on Omega (Sec. 6)");
    const FieldIntegration::Report &ir = pipeline.getIntegration().getReport();
    std::cout << "  " << ir.vertices << " vertices, " << ir.constraintRows
              << " constraint row(s) over " << ir.seamPairs << " seam pair(s); solved with "
              << (ir.solvedWithLDLT ? "LDL^T" : "the LU fallback") << "\n";
    std::cout << "  Seam residual " << std::scientific << std::setprecision(3)
              << ir.maxSeamResidual << " of the extent" << std::defaultfloat << "\n";
    if (ir.alignedEdges > 0) {
        std::cout << "  Sec. 6.4: " << ir.alignedEdges << " edge(s) held to an axis; residual "
                  << std::scientific << std::setprecision(3) << ir.maxAlignResidual
                  << " of the extent, strain " << ir.alignmentStrain << " relative"
                  << std::defaultfloat << "\n";
    } else if (st.alignmentWasDropped) {
        std::cout << "  Sec. 6.4: the alignment was dropped -- see the message below\n";
    }
    std::cout << "  Non-integrability, per triangle: " << std::scientific << std::setprecision(3)
              << ir.maxFitResidual << " worst, " << ir.meanFitResidual << " mean"
              << std::defaultfloat;
    if (ir.worstFitFace >= 0) std::cout << " (worst at face " << ir.worstFitFace << ")";
    std::cout << "\n";
    std::cout << "  Sec. 7.1 census: " << ir.flippedFaces << " inverted face(s), "
              << std::fixed << std::setprecision(3)
              << (100.0 * st.integrationFlippedAreaFraction) << "% of the model by area; "
              << ir.flipsAdjacentToCone << " in a cone's one ring"
              << std::defaultfloat << "\n";
    if (ir.flippedFaces > 0) {
        std::cout << "     nearest to a cone " << std::scientific << std::setprecision(3)
                  << ir.nearestFlipToCone << " of the diagonal, farthest "
                  << ir.farthestFlipFromCone << std::defaultfloat << "\n";
    }
    verdict(ir.solved, "The constrained least-squares system solved");
    verdict(ir.maxSeamResidual < 1e-9, "C2: the seam holds to rounding, not to a tolerance");
    // Not a verdict. A discrete cross field is generically non-integrable, so
    // flips here are the expected cost of the substitution -- Sec. 7 is written
    // to pay it, and Sec. 7.1 says in as many words not to assume zero.
    if (ir.flippedFaces == 0) {
        std::cout << "  " << kPass << " Nothing inverted: no untangling needed on this model\n";
    } else {
        warn("Sec. 7: the field is not integrable here, so the least-squares map inverted "
             + std::to_string(ir.flippedFaces) + " face(s). Stage 4R is what pays for that.");
    }
    stageMessages("Stage 4F: ");
    if (!rawOut.empty()) {
        std::ofstream out(rawOut);
        if (out) {
            for (const Point &p : pipeline.getIntegratedMap()) out << "v " << p[0] << " " << p[1] << " 0\n";
            for (const Triangle &t : pipeline.getCut().getCutMesh().triangles) {
                out << "f " << (t[0] + 1) << " " << (t[1] + 1) << " " << (t[2] + 1) << "\n";
            }
            std::cout << "  Wrote the raw integrated map to " << rawOut << "\n";
        } else warn("Failed to write " + rawOut);
    }

    // ---------------------------------------------------------------------
    if (st.untangleRan || pipeline.hasTutte() || st.localUntangleRan || st.pullRan) {
        heading("Stage 4R  Untangling: the local kernel pass, then Tutte (Secs. 7.2a, 7.2)");
        if (st.localUntangleRan) {
            std::cout << "  Sec. 7.2a: " << st.localUntangleFlippedFaces
                      << " inverted face(s) left after moving the interior vertices around "
                      << "the tangle to the centres of their one-ring kernels\n";
            verdict(st.localUntangleFlippedFaces == 0 || pipeline.hasTutte(),
                    "The tangle was cleared locally, or handed on to Sec. 7.2");
        }
        if (pipeline.hasTutte()) {
            const TutteEmbedding::Report &tr = pipeline.getTutte().getReport();
            std::cout << "  Tutte: " << tr.boundaryVertices << " boundary + "
                      << tr.interiorVertices << " interior vertices onto a circle of radius "
                      << std::fixed << std::setprecision(4) << tr.radius << std::defaultfloat
                      << "; " << tr.flippedFaces << " inverted face(s)\n";
            verdict(tr.valid, "Tutte's theorem: the starting map is bijective");
            if (!tutteOut.empty()) {
                std::ofstream out(tutteOut);
                if (out) {
                    for (const Point &p : pipeline.getTutte().getUV()) out << "v " << p[0] << " " << p[1] << " 0\n";
                    for (const Triangle &t : pipeline.getCut().getCutMesh().triangles) {
                        out << "f " << (t[0] + 1) << " " << (t[1] + 1) << " " << (t[2] + 1) << "\n";
                    }
                    std::cout << "  Wrote the Tutte embedding to " << tutteOut << "\n";
                } else warn("Failed to write " + tutteOut);
            }
        }
        if (st.untangleRan) {
            std::cout << "  Target fit: " << st.untangleOuterSteps
                      << " outer step(s), " << st.untangleFlippedFaces
                      << " inverted face(s) left, seam " << std::scientific
                      << std::setprecision(3) << st.untangleSeamResidual
                      << " of the extent" << std::defaultfloat << "\n";
            verdict(st.untangleFlippedFaces == 0,
                    "The barrier's flip cap kept every triangle positively oriented");
        }
        if (st.pullRan) {
            std::cout << "  Sec. 7.2b: pulled back onto the alignment to Q3 " << std::scientific
                      << std::setprecision(3) << st.pullBoundaryResidual << ", features "
                      << st.pullFeatureResidual << std::defaultfloat
                      << (st.pullProjected ? ", then projected onto it exactly\n" : "\n");
            verdict(st.pullKept, "Sec. 7.2b: the pull kept every triangle positively oriented");
        }
        stageMessages("Stage 4R: ");
    }

    // ---------------------------------------------------------------------
    if (!pipeline.hasImmersion()) {
        heading("Result");
        std::cout << "  " << kFail << " Stage 4 could not accept psi_0.\n";
        return 6;
    }
    heading("Stage 4  psi_0 handed to the layout energies");
    const Immersion::Report &imr = pipeline.getImmersion().getReport();
    std::cout << "  Image " << std::fixed << std::setprecision(3)
              << (imr.uvMax[0] - imr.uvMin[0]) << " x " << (imr.uvMax[1] - imr.uvMin[1])
              << std::defaultfloat << "; " << imr.arcs << " arc(s), " << imr.seamEdgePairs
              << " seam pair(s); Gamma_Hol_0..3 = " << imr.holonomyCount[0] << ", "
              << imr.holonomyCount[1] << ", " << imr.holonomyCount[2] << ", "
              << imr.holonomyCount[3] << "\n";
    std::cout << "  Field metric realised to " << std::scientific << std::setprecision(3)
              << imr.maxMetricResidual << " relative (this is the non-integrability, per edge)"
              << std::defaultfloat << "\n";
    std::cout << "  Q4 as a check rather than a rounding: the fit disagrees with the "
              << "matchings' quarter turn by " << std::scientific << std::setprecision(3)
              << imr.maxSnapError << " rad at worst; frame residual " << imr.maxFrameKResidual
              << std::defaultfloat << "\n";
    std::cout << "  Angle sums off target by " << std::scientific << std::setprecision(3)
              << imr.maxConeAngleResidual << " rad at the cones, "
              << imr.maxRegularAngleResidual << " elsewhere" << std::defaultfloat << "\n";
    verdict(imr.flippedFaces == 0, "Q1: no face of psi_0 is inverted");
    verdict(imr.frameKConflicts == 0,
            "C2: every seam edge of an arc agrees with it about the quarter turn");
    verdict(imr.degenerateArcs == 0, "Every arc of G paired into two sides");
    if (!psiOut.empty()) {
        if (pipeline.getImmersion().writeOBJ(psiOut)) std::cout << "  Wrote psi_0 to " << psiOut << "\n";
        else warn("Failed to write " + psiOut);
    }
    stageMessages("Stage 4: ");

    if (!pipeline.hasLayout()) {
        heading("Result");
        const bool immOk = imr.flippedFaces == 0;
        std::cout << "  " << (immOk ? kPass : kFail)
                  << " psi_0 " << (immOk ? "satisfies Q1 and Q4" : "does not satisfy Q1")
                  << "; the layout energies were not run.\n";
        return immOk ? 0 : 6;
    }

    // ---------------------------------------------------------------------
    heading("Stage 5  Subdomain labelling (Sec. 3.3)");
    std::cout << "  Gamma_u " << st.boundaryEdgesU << " edge(s), Gamma_v "
              << st.boundaryEdgesV << "; " << st.featureChains << " feature chain(s), "
              << st.topoPaths << " Gamma_topo path(s)\n";
    if (st.materials > 1) {
        std::cout << "  " << st.interfaceCorners << " interface sector(s) for E6; "
                  << st.interfaceLabelsCorrected << " chain label(s) where the flux said otherwise, "
                  << st.interfaceChainsSeamFlipped << " turned by a crossing of G\n";
    }
    std::cout << "  Near-miss tolerance " << st.topoNearMissUsed
              << (st.topoNearMissRetried ? "  (after a retry)" : "") << "\n";
    stageMessages("Stage 5: ");

    // ---------------------------------------------------------------------
    heading("Stage 6  Layout-inducing energies against the field reference");
    const LayoutEnergy::Report &er = pipeline.getLayout().getReport();
    const char *refName = "induced (psi_0's own metric)";
    switch (opts.reference) {
        case TORSION::Options::Reference::Cone:      refName = "the flat cone metric (Sec. 4)"; break;
        case TORSION::Options::Reference::Field:     refName = "the field frame (Sec. 4)"; break;
        case TORSION::Options::Reference::Euclidean: refName = "the model's Euclidean geometry"; break;
        case TORSION::Options::Reference::Ricci:     refName = "a Ricci flat cone metric"; break;
        default: break;
    }
    std::cout << "  Reference: " << refName << "; " << er.outerSteps
              << " outer step(s), " << er.innerIterations << " inner iteration(s)\n";
    // C4 in one number, taken on the reference metric itself rather than
    // inferred from what the continuation produced: how far the reference's own
    // angle sum at each cone is from the 2 pi - (pi/2) I the layout needs. A
    // reference with no cones reports (pi/2)|I| here.
    std::cout << "  C4: the reference's own cone angles are off by "
              << std::scientific << std::setprecision(3) << st.referenceConeResidual
              << " rad" << std::defaultfloat << "\n";
    verdict(st.referenceConeResidual < 1e-3,
            "C4: E1's reference metric has the cones Q2 asks for");
    std::cout << "  Residuals: Q3 " << std::scientific << std::setprecision(3)
              << er.maxBoundaryResidual << ", features " << er.maxFeatureResidual
              << ", Q4 " << er.maxSeamResidual << ", Q5 " << er.maxTopoResidual
              << std::defaultfloat << "\n";
    std::cout << "  Q2: " << std::scientific << std::setprecision(3)
              << er.maxConeAngleResidual << " rad at the cones, " << er.maxRegularAngleResidual
              << " elsewhere; " << er.coneValenceChanges << " cone(s) changed valence"
              << std::defaultfloat << "\n";
    std::cout << "  min det J " << std::scientific << std::setprecision(3) << er.minDetJ
              << std::defaultfloat << "; " << er.invertedTriangles << " inverted triangle(s)\n";
    verdict(er.injective, "Q1: det J > 0 on every triangle");
    verdict(er.constraintsMet, "Q3, Q4, Q5 and the features are under tolerance");
    verdict(er.anglesHeld, "Q2: cone angle sums are the prescribed multiples of pi/2");
    stageMessages("Stage 6: ");
    if (er.leftHandedFrames > 0) {
        warn(std::to_string(er.leftHandedFrames) +
             " triangle(s) fell back to the Euclidean reference for want of a valid frame");
    }
    if (!layoutOut.empty()) {
        if (pipeline.getLayout().writeOBJ(layoutOut)) std::cout << "  Wrote Psi to " << layoutOut << "\n";
        else warn("Failed to write " + layoutOut);
    }

    return -1;
}

// ---------------------------------------------------------------------------
// Stages 1 to 8 as the per-material mode ran them: one line per material
// region, then what the matching across the interfaces came to.
// ---------------------------------------------------------------------------
static void reportPerMaterial(const TORSION &pipeline) {
    const TORSION::Status &st = pipeline.getStatus();
    heading("Stages 1 to 8, one material region at a time (Sec. 15)");
    if (!pipeline.hasMaterialLayout()) {
        std::cout << "  " << kFail << " The regions were not laid out.\n";
        return;
    }
    const MaterialLayout &ml = pipeline.getMaterialLayout();
    const MaterialLayout::Report &mr = ml.getReport();
    std::cout << "  " << mr.regions << " region(s); " << mr.regionRuns << " layout(s) from Stage 1 and "
              << mr.regionRelayouts << " from Stage 5 again over " << mr.rounds
              << " round(s) of matching, " << std::fixed << std::setprecision(1) << mr.seconds
              << " s" << std::defaultfloat << "\n";
    if (!mr.askedPerRound.empty()) {
        std::cout << "  Layout edges asked per round:";
        for (int a : mr.askedPerRound) std::cout << " " << a;
        std::cout << "\n";
    }
    std::cout << "  " << mr.emitters << " layout edge(s) asked of a region by its neighbours; "
              << mr.matched << " layout vertex/vertices matched across an interface, "
              << mr.unmatched << " left single; worst match "
              << std::fixed << std::setprecision(2) << mr.maxMatchOffset << " edge(s)"
              << std::defaultfloat << "\n";
    const std::vector<MaterialLayout::Region> &regions = ml.getRegions();
    for (size_t r = 0; r < regions.size(); ++r) {
        const MaterialLayout::Region &R = regions[r];
        std::cout << "    r" << r << "  material " << R.material << ", " << R.triangles.size()
                  << " triangle(s)";
        if (!R.layout) {
            std::cout << ": not laid out\n";
            continue;
        }
        const TORSION::Status &rs = R.layout->getStatus();
        std::cout << "; " << rs.interiorCones << " + " << rs.boundaryCones << " cone(s)";
        if (R.laidOut()) {
            const Arrangement::Report &ar = R.layout->getArrangement().getReport();
            std::cout << "; " << ar.patches << " patch(es), " << ar.simpleQuads << " simple";
        }
        std::cout << "; " << R.emitters.size() << " emitter(s)"
                  << (R.dipolesKept ? "; +1/-1 pairs kept" : "")
                  << (rs.layoutValid ? "" : "; Stage 6 short of Definition 2.1") << "; "
                  << std::fixed << std::setprecision(2) << R.seconds << " s" << std::defaultfloat
                  << "\n";
    }
    verdict(mr.regionsValid == mr.regions,
            "Every region's own layout reaches Definition 2.1 with every face a simple quadrilateral");
    verdict(mr.converged, "The matching settled: no region is asked for another layout edge");
    verdict(mr.unmatched == 0, "Every layout vertex on an interface is matched on the other side");
    for (const std::string &m : st.messages) {
        if (m.rfind("Per material: ", 0) == 0) warn(m.substr(14));
    }
    if (st.perMaterialFallbackRan) {
        std::cout << "  The whole model was laid out as well; the "
                  << (st.perMaterialKept ? "glued" : "whole-model") << " layout is the one kept.\n";
    }
}

int main(int argc, char **argv) {
    if (argc < 2) { usage(argv[0]); return 1; }
    const std::string first = argv[1];
    if (first == "--selftest") return selfTest();
    if (first == "--ref-test") {
        if (argc < 3) { usage(argv[0]); return 1; }
        return referenceTest(argv[2]);
    }

    const std::string path = first;
    TORSION::Options opts;
    std::string psiOut, rawOut, tutteOut, layoutOut, cutOut, quadOut, arrOut, mfemOut;
    int coneListLimit = 20;
    int tmopSweeps = 0;
    bool tmopGauss = false;
    bool tmopPillow = true;

    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--h" && i + 1 < argc)                 opts.targetEdge = std::stod(argv[++i]);
        else if (a == "--no-untangle")                  opts.untangle = false;
        else if (a == "--untangle-outer" && i + 1 < argc) opts.untangleOuterSteps = std::stoi(argv[++i]);
        else if (a == "--untangle-mu" && i + 1 < argc)  opts.untangleSeamWeight = std::stod(argv[++i]);
        else if (a == "--reg" && i + 1 < argc)          opts.integrationRegularisation = std::stod(argv[++i]);
        else if (a == "--no-seed-check")                opts.checkSecondSeed = false;
        else if (a == "--no-align")                     opts.alignInIntegration = false;
        else if (a == "--no-conformal-sizing")          opts.conformalSizing = false;
        else if (a == "--no-align-fallback")            opts.alignmentFallback = false;
        else if (a == "--no-seam-align")                opts.alignAcrossSeams = false;
        else if (a == "--choose-by-flips")              opts.alignmentChooseByLadder = false;
        else if (a == "--no-regularised-untangle")      opts.regularisedUntangleRings = 0;
        else if (a == "--no-sector-reconcile")          opts.reconcileSectors = false;
        else if (a == "--release-rounds" && i + 1 < argc) opts.alignmentReleaseRounds = std::stoi(argv[++i]);
        else if (a == "--no-pull")                      opts.pullOntoAlignment = false;
        else if (a == "--whole-model")                  opts.perMaterial = false;
        else if (a == "--pm-rounds" && i + 1 < argc)    opts.perMaterialRounds = std::stoi(argv[++i]);
        else if (a == "--pm-threads" && i + 1 < argc)   opts.perMaterialThreads = std::stoi(argv[++i]);
        else if (a == "--pm-tolerance" && i + 1 < argc) opts.perMaterialMatchTolerance = std::stod(argv[++i]);
        else if (a == "--pm-no-fallback")               opts.perMaterialFallback = false;
        else if (a == "--pm-cancel-always")             opts.perMaterialKeepDipoles = false;
        else if (a == "--pm-most-pairs")                opts.perMaterialSpacedMatching = false;
        else if (a == "--pm-move-ends")                 opts.perMaterialBendEnds = false;
        else if (a == "--ref" && i + 1 < argc) {
            const std::string r = argv[++i];
            if (r == "cone")           opts.reference = TORSION::Options::Reference::Cone;
            else if (r == "induced")   opts.reference = TORSION::Options::Reference::Induced;
            else if (r == "field")     opts.reference = TORSION::Options::Reference::Field;
            else if (r == "euclidean") opts.reference = TORSION::Options::Reference::Euclidean;
            else if (r == "ricci")     opts.reference = TORSION::Options::Reference::Ricci;
            else { std::cerr << "Unknown reference: " << r << "\n"; return 1; }
        }
        else if (a == "--gamma" && i + 1 < argc)        opts.dualMBOGamma = std::stod(argv[++i]);
        else if (a == "--steps" && i + 1 < argc)        opts.dualMBOMaxSteps = std::stoi(argv[++i]);
        else if (a == "--weight" && i + 1 < argc) {
            const std::string w = argv[++i];
            if      (w == "min")  opts.dualMBOWeight = DualMBO::PenaltyWeight::MinHeight;
            else if (w == "harm") opts.dualMBOWeight = DualMBO::PenaltyWeight::HarmonicHeight;
            else if (w == "orth") opts.dualMBOWeight = DualMBO::PenaltyWeight::Orthogonal;
            else { std::cerr << "Unknown --weight '" << w << "' (expected min|harm|orth)\n"; return 1; }
        }
        else if (a == "--no-continuation")              opts.dualMBOTauContinuation = false;
        else if (a == "--tau-ratio" && i + 1 < argc)     opts.dualMBOTauRatio = std::stod(argv[++i]);
        else if (a == "--tau-floor" && i + 1 < argc)     opts.dualMBOTauFloorEdges = std::stod(argv[++i]);
        else if (a == "--cut-to-graph")                 opts.coneCutsToBoundary = false;
        else if (a == "--no-interfaces")                opts.materialInterfaces = false;
        else if (a == "--no-field-interfaces")          opts.alignFieldToInterfaces = false;
        else if (a == "--no-cancel-dipoles")            opts.cancelInterfaceDipoles = false;
        else if (a == "--no-dipole-retry")              opts.retryKeepingDipoles = false;
        else if (a == "--keep-flat-cones")              opts.relocateFlatCones = false;
        else if (a == "--no-e6")                        opts.interfaceCorners = false;
        else if (a == "--no-seam-turns")                opts.seamTurnInterfaceLabels = false;
        else if (a == "--outer" && i + 1 < argc)        opts.outerSteps = std::stoi(argv[++i]);
        else if (a == "--inner" && i + 1 < argc)        opts.innerIterations = std::stoi(argv[++i]);
        else if (a == "--lambda" && i + 1 < argc)       opts.lambdaInit = std::stod(argv[++i]);
        else if (a == "--growth" && i + 1 < argc)       opts.lambdaGrowth = std::stod(argv[++i]);
        else if (a == "--align-factor" && i + 1 < argc) opts.lambdaAlignmentFactor = std::stod(argv[++i]);
        else if (a == "--seam-factor" && i + 1 < argc)  opts.lambdaSeamFactor = std::stod(argv[++i]);
        else if (a == "--no-topo")                      opts.seedTopoConstraints = false;
        else if (a == "--no-retry") {
            opts.topoNearMissRetry = 0.0;
            opts.topoNearMissLastRetry = 0.0;
            opts.topoRetryUnseeded = false;
        }
        else if (a == "--retry" && i + 1 < argc)       opts.topoNearMissRetry = std::stod(argv[++i]);
        else if (a == "--last-retry" && i + 1 < argc)  opts.topoNearMissLastRetry = std::stod(argv[++i]);
        else if (a == "--no-unseeded-retry")           opts.topoRetryUnseeded = false;
        else if (a == "--near-miss" && i + 1 < argc)    opts.topoNearMiss = std::stod(argv[++i]);
        else if (a == "--no-layout")                    opts.runLayout = false;
        else if (a == "--no-trace")                     opts.runSeparatrices = false;
        else if (a == "--no-arrange")                   opts.runArrangement = false;
        else if (a == "--no-splines")                   opts.runSplines = false;
        else if (a == "--no-mesh")                      opts.runQuadMesh = false;
        else if (a == "--collapse-span" && i + 1 < argc) opts.quadCollapseSpan = std::stod(argv[++i]);
        else if (a == "--no-collapse")                 opts.quadCollapseSpan = 0.0;
        else if (a == "--sample-materials")            opts.quadMaterialsFromFaces = false;
        else if (a == "--contract-at-root")            opts.quadContractOntoFeatures = false;
        else if (a == "--target" && i + 1 < argc)       opts.quadTargetEdge = std::stod(argv[++i]);
        else if (a == "--disk-templates")               opts.diskTemplates = true;
        else if (a == "--disk-squareness" && i + 1 < argc)
            opts.diskCoreSquareness = std::stod(argv[++i]);
        else if (a == "--disk-ring" && i + 1 < argc)    opts.diskRingDepth = std::stoi(argv[++i]);
        else if (a == "--disk-smooth" && i + 1 < argc)
            opts.diskSmoothingPasses = std::stoi(argv[++i]);
        else if (a == "--repair" && i + 1 < argc)       opts.repairPasses = std::stoi(argv[++i]);
        else if (a == "--cones" && i + 1 < argc)        coneListLimit = std::stoi(argv[++i]);
        else if (a == "--psi" && i + 1 < argc)          psiOut = argv[++i];
        else if (a == "--raw" && i + 1 < argc)          rawOut = argv[++i];
        else if (a == "--tutte" && i + 1 < argc)        tutteOut = argv[++i];
        else if (a == "--layout" && i + 1 < argc)       layoutOut = argv[++i];
        else if (a == "--cut" && i + 1 < argc)          cutOut = argv[++i];
        else if (a == "--quads" && i + 1 < argc)        quadOut = argv[++i];
        else if (a == "--mfem" && i + 1 < argc)         mfemOut = argv[++i];
        else if (a == "--arr" && i + 1 < argc)          arrOut = argv[++i];
        else if (a == "--tmop" && i + 1 < argc)         tmopSweeps = std::stoi(argv[++i]);
        else if (a == "--tmop-gauss")                    tmopGauss = true;
        else if (a == "--no-pillow")                     tmopPillow = false;
        else { std::cerr << "Unknown option: " << a << "\n"; usage(argv[0]); return 1; }
    }

    std::shared_ptr<Mesh> mesh;
    try {
        mesh = std::make_shared<Mesh>(path);
    } catch (const std::exception &e) {
        std::cout << kFail << " Failed to load mesh: " << e.what() << "\n";
        return 2;
    }

    std::cout << "TORSION -- Pipeline B: the layout of Shepherd, Gu and Hughes (2022) with\n"
              << "psi_0 integrated from the DualMBO cross field (docs/cf_flow_pipeline.md)\n";
    std::cout << "Mesh: " << path << "\n";
    std::cout << "  " << mesh->vertices.size() << " vertices, "
              << mesh->edges.size() << " edges, "
              << mesh->triangles.size() << " triangles, "
              << mesh->boundaryEdges.size() << " boundary edges\n";

    TORSION pipeline(mesh, opts);
    bool ok = false;
    try {
        ok = pipeline.run();
    } catch (const std::exception &e) {
        std::cout << kFail << " Pipeline threw: " << e.what() << "\n";
        return 3;
    }
    const TORSION::Status &st = pipeline.getStatus();

    auto stageMessages = [&](const std::string &prefix) {
        for (const std::string &m : st.messages) {
            if (m.rfind(prefix, 0) == 0) warn(m.substr(prefix.size()));
        }
    };

    // ---------------------------------------------------------------------
    heading("Stage 0  Cross field (p=0 dual-mesh MBO)");
    std::cout << "  MBO steps: " << st.mboSteps << ", residual " << std::scientific
              << std::setprecision(3) << pipeline.getField().error << std::defaultfloat << "\n";
    verdict(st.fieldConverged, "Field converged");
    if (st.materials > 1) {
        std::cout << "  " << st.materials << " material(s) in " << st.regions
                  << " region(s); " << st.interfaceEdges << " interface edge(s) in "
                  << st.interfaceBranches << " branch(es), " << st.interfaceNodes << " node(s)\n";
        verdict(st.fieldAlignedToInterfaces, "The interfaces are Dirichlet data for the field");
        verdict(st.regionsBalanced == st.regions,
                "Every material region satisfies its own Eq. (4)");
    }

    // ---------------------------------------------------------------------
    // Stages 1 to 6 (or 1 to 8), by the route that produced the layout kept.
    if (st.perMaterialRan) reportPerMaterial(pipeline);
    if (!st.perMaterialRan || !st.perMaterialKept) {
        WholeModelOutputs outputs;
        outputs.cutOut = cutOut;
        outputs.rawOut = rawOut;
        outputs.tutteOut = tutteOut;
        outputs.psiOut = psiOut;
        outputs.layoutOut = layoutOut;
        outputs.coneListLimit = coneListLimit;
        const int rc = reportWholeModel(pipeline, mesh, opts, outputs);
        if (rc >= 0) return rc;
    }

    // ---------------------------------------------------------------------
    if (st.separatricesRan) {
        heading("Stage 7  Separatrices (Sec. 4)");
        std::cout << "  " << st.separatrices << " curve(s): " << st.separatricesToCone
                  << " to a cone, " << st.separatricesToBoundary << " out through dS, "
                  << st.separatricesUnresolved << " terminating nowhere, "
                  << st.separatricesNearMisses << " grazing a cone\n";
        std::cout << "  Repair: " << st.repairPasses << " pass(es), "
                  << st.repairConstraintsAdded << " constraint(s) added\n";
        verdict(st.q5Verified, "Q5: every curve ends at a cone or leaves through dS");
        stageMessages("Stage 7: ");
    }
    if (st.arrangementRan) {
        heading("Stage 8  Arrangement (Sec. 4)");
        std::cout << "  " << st.layoutNodes << " node(s), " << st.layoutArcs << " arc(s), "
                  << st.layoutPatches << " patch(es) of which " << st.layoutQuads
                  << " are simple quadrilaterals; coverage " << std::fixed
                  << std::setprecision(4) << st.layoutCoverage << std::defaultfloat << "\n";
        verdict(st.arrangementValid, "The curves partition S into four-sided faces");
        if (!arrOut.empty()) {
            if (pipeline.getArrangement().writeOBJ(arrOut)) std::cout << "  Wrote the arrangement to " << arrOut << "\n";
            else warn("Failed to write " + arrOut);
        }
        stageMessages("Stage 8: ");
    }
    if (st.splinesRan) {
        heading("Stage 9  Spline reconstruction (Sec. 5)");
        const SplineFit::Report &sr = pipeline.getSplines().getReport();
        std::cout << "  " << st.splinePatches << " bicubic patch(es), "
                  << st.splineControlPoints << " control point(s) per arc; worst deviation "
                  << std::scientific << std::setprecision(3) << st.splineMaxDeviation
                  << std::defaultfloat << "\n";
        std::cout << "  " << sr.fittedArcs << " arc(s) fitted, " << sr.exactArcs
                  << " carried exactly; worst sampled cell " << std::fixed
                  << std::setprecision(4) << sr.minCellRatio << " of the mean"
                  << std::defaultfloat << "\n";
        std::cout << "  B-rep: " << sr.brepFaces << " face(s), " << sr.brepEdges << " edge(s) -- "
                  << sr.brepSharedEdges << " shared by two faces, " << sr.brepFreeEdges
                  << " bounding one; the kernel calls it " << (sr.brepValid ? "valid" : "INVALID")
                  << "\n";
        verdict(st.splinesWatertight, "Watertight: both patches on an arc share its control points");
        verdict(sr.foldedPatches == 0, "No patch folds under its Coons blend");
        if (sr.foldedPatches > 0) {
            warn(std::to_string(sr.foldedPatches) +
                 " patch(es) fold: the four arcs they are built on do not bound a "
                 "quadrilateral the blend can fill, so Stage 10 meshes a fold.");
        }
        stageMessages("Stage 9: ");
    }
    if (st.quadMeshRan) {
        heading("Stage 10  Quadrilateral mesh");
        std::cout << "  " << st.meshVertices << " vertices, " << st.meshQuads
                  << " quadrilateral(s) over " << st.meshChords << " chord(s); "
                  << st.meshUnmeshedPatches << " patch(es) unmeshed\n";
        const QuadMesh::Report &qr2 = pipeline.getQuadMesh().getReport();
        std::cout << "  Scaled Jacobian " << std::fixed << std::setprecision(4)
                  << st.meshMinScaledJacobian << " worst, " << qr2.meanScaledJacobian
                  << " mean" << std::defaultfloat << "\n";
        std::cout << "  " << qr2.invertedQuads << " element(s) of non-positive area, "
                  << qr2.reflexCorners << " reflex corner(s) of the layout\n";
        std::cout << "  Edge length in [" << std::fixed << std::setprecision(4)
                  << qr2.minEdge << ", " << qr2.maxEdge << "], rms log ratio "
                  << qr2.edgeRatioRms << std::defaultfloat << "\n";
        if (qr2.collapsedChords > 0) {
            std::cout << "  Contracted " << qr2.collapsedChords << " chord(s), which merged "
                      << qr2.collapsedPatches << " face(s) -- " << std::fixed
                      << std::setprecision(2) << (100.0 * qr2.collapsedArea)
                      << "% of S -- into their neighbours on " << qr2.weldedVertices
                      << " welded vertex/vertices" << std::defaultfloat << "\n";
        }
        verdict(st.meshConforming, "Conforming: every edge is shared by two elements or bounds the mesh");
        verdict(st.meshValid, "No element has non-positive area");
        // An element with a positive area can still have a reversed corner, and
        // a solver will not take it. The area test is what `valid` is written
        // on; this is the test the mesh is actually used under.
        verdict(st.meshMinScaledJacobian > 0.0,
                "Every corner of every element turns the right way (scaled Jacobian > 0)");
        stageMessages("Stage 10: ");
    }
    if (pipeline.hasDiskTemplate()) {
        heading("Stage 11  O-grid templates on the excised inclusions");
        const DiskTemplate::Report &dr = pipeline.getDiskTemplate().getReport();
        std::cout << "  " << dr.filled << " of " << dr.inclusions
                  << " inclusion(s) templated: " << dr.blocks << " block(s), "
                  << dr.quads << " element(s), " << dr.vertices << " new vertex/vertices\n";
        std::cout << "  Merged mesh: " << dr.mergedVertices << " vertices, "
                  << dr.mergedQuads << " quadrilateral(s); scaled Jacobian "
                  << std::fixed << std::setprecision(4) << dr.minScaledJacobian
                  << " worst, " << dr.meanScaledJacobian << " mean ("
                  << dr.templateMinScaledJacobian << " worst on the templates)"
                  << std::defaultfloat << "\n";
        verdict(dr.filled == dr.inclusions, "Every circular inclusion was templated");
        verdict(dr.nonManifoldEdges == 0 && dr.cracks == 0,
                "Conforming and watertight across the rims");
        verdict(dr.invertedQuads == 0, "No element of the merged mesh is inverted");
        stageMessages("Stage 11: ");
    }
    // Whatever came last -- the Stage 10 mesh, or the matrix plus its filled
    // inclusions -- adopted into the standalone quad mesh class of
    // src/mesh/QuadMesh.hxx, the same conversion TestMERIDIAN does. Its
    // topology is rebuilt from the elements alone, so this is a second opinion
    // on the mesh rather than a restatement of the pipeline's own report.
    if (st.quadMeshRan && pipeline.hasQuadMesh()) {
        heading("Final mesh (mesh::QuadMesh)");
        mesh::QuadMesh finalMesh =
            pipeline.hasDiskTemplate() ? mesh::QuadMesh::from(pipeline.getDiskTemplate())
                                       : mesh::QuadMesh::from(pipeline.getQuadMesh());
        const mesh::QuadMesh::Quality &fq = finalMesh.quality;
        std::cout << "  " << fq.vertices << " vertices, " << fq.quads
                  << " quadrilateral(s), " << finalMesh.edges.size() << " edge(s)\n";
        std::cout << "  Nodes: " << fq.freeNodes << " free, " << fq.slidingNodes
                  << " sliding on a feature, " << fq.fixedNodes << " fixed\n";
        std::cout << "  Scaled Jacobian: worst " << std::fixed << std::setprecision(3)
                  << fq.minScaledJacobian << ", mean " << fq.meanScaledJacobian
                  << std::defaultfloat << "\n";
        std::cout << "  Boundary: " << fq.boundaryEdges << " edge(s) in " << fq.boundaryLoops
                  << " loop(s)\n";
        std::cout << "  Materials:";
        for (const std::pair<int, int> &kv : finalMesh.materialCounts())
            std::cout << " " << kv.first << " (" << kv.second << " element(s))";
        std::cout << "\n";
        verdict(fq.allCounterClockwise && fq.invertedQuads == 0,
                "Every element is counter-clockwise with a positive Jacobian at every corner");
        verdict(fq.nonManifoldEdges == 0, "Every edge is shared by at most two elements");

        // Stage 12 starts with a pillow, as TestMERIDIAN's does: one layer of
        // quads along each stretch of feature on which a block corner spans
        // 180 degrees, so that its node has two elements and the corner is
        // half an element inside, where TMOP can move it. src/mesh/Pillow.hxx.
        if (tmopSweeps > 0 && tmopPillow) {
            heading("Pillowing flat feature corners (mesh::Pillow)");
            mesh::Pillow pillow(finalMesh);
            pillow.run();
            const mesh::Pillow::Report &pr = pillow.getReport();
            std::cout << "  " << pr.defectsBefore << " element corner(s) with both sides on a "
                      << "feature, opening past " << pillow.getOptions().flatAngle << " degrees";
            if (pr.defectsBefore > 0)
                std::cout << " (the flattest " << std::fixed << std::setprecision(1)
                          << pr.flattestBefore << std::defaultfloat << ")";
            std::cout << "\n";
            if (pr.layers > 0) {
                std::cout << "  " << pr.layers << " layer(s)";
                if (pr.closedLayers > 0)
                    std::cout << ", " << pr.closedLayers << " round a closed loop";
                std::cout << ": " << pr.layerQuads << " quad(s) along the features, "
                          << pr.splitQuads << " split carrying the layers' chords on to dS; "
                          << pr.quadsBefore << " -> " << pr.quadsAfter << " quadrilateral(s)\n";
                std::cout << "  Scaled Jacobian before smoothing: worst " << std::fixed
                          << std::setprecision(4) << pr.minScaledJacobianBefore << " -> "
                          << pr.minScaledJacobianAfter << std::defaultfloat << "\n";
            }
            for (const std::string &m : pr.messages) warn(m);
            verdict(pr.defectsAfter == 0,
                    "No element spans a flat corner with both of its sides on a feature");
        }

        if (tmopSweeps > 0) {
            heading("TMOP smoothing (mesh::TMOP)");
            mesh::TMOP::Options topt;
            topt.metric = mesh::TMOP::ShapeSize007;
            topt.maxSweeps = tmopSweeps;
            // Sampled at the corners, as the viewer's Stage 12: at the 2x2
            // Gauss points a corner can fold without the barrier seeing it.
            topt.quadrature = tmopGauss ? mesh::TMOP::Gauss2x2 : mesh::TMOP::Corners;
            mesh::TMOP smoother(finalMesh, topt);
            smoother.run();
            const mesh::TMOP::Report &tr = smoother.getReport();
            std::cout << "  " << tr.sweeps << " sweep(s)";
            if (tr.untangleSweeps > 0)
                std::cout << " after " << tr.untangleSweeps << " untangling sweep(s)";
            std::cout << ", " << tr.colors << " colour(s), " << tr.threads << " thread(s)"
                      << (tr.openMP ? "" : " (no OpenMP)") << ", " << std::fixed
                      << std::setprecision(3) << tr.seconds << " s\n";
            std::cout << "  Scaled Jacobian: worst " << std::fixed << std::setprecision(4)
                      << tr.minScaledJacobianBefore << " -> " << tr.minScaledJacobianAfter
                      << ", mean " << tr.meanScaledJacobianBefore << " -> "
                      << tr.meanScaledJacobianAfter << std::defaultfloat << "\n";
            for (const std::string &m : tr.messages) warn(m);
            verdict(tr.invertedAfter == 0, "No element of the smoothed mesh is inverted");
        }

        if (!mfemOut.empty()) {
            mesh::QuadMesh::MFEMOptions mo;
            mesh::QuadMesh::MFEMReport mr;
            if (finalMesh.writeMFEM(mfemOut, mo, &mr)) {
                std::cout << "  Wrote the MFEM mesh to " << mfemOut << ": " << mr.elements
                          << " element(s), " << mr.boundaryElements << " boundary segment(s), "
                          << mr.vertices << " vertices\n";
                for (const std::pair<int, int> &kv : mr.attributeCounts)
                    std::cout << "    attribute " << kv.first << ": " << kv.second
                              << " element(s)\n";
                if (mr.attributesRemapped)
                    warn("material ids were shifted so the MFEM attribute is positive");
                if (mr.unusedVertices > 0)
                    warn(std::to_string(mr.unusedVertices) +
                         " vertex/vertices no element uses were dropped");
            } else {
                warn("Failed to write " + mfemOut);
            }
        }
    }

    if (!quadOut.empty() && st.quadMeshRan) {
        const bool wrote = pipeline.hasDiskTemplate()
                               ? pipeline.getDiskTemplate().writeOBJ(quadOut)
                               : pipeline.getQuadMesh().writeOBJ(quadOut);
        if (wrote) std::cout << "  Wrote the quad mesh to " << quadOut << "\n";
        else warn("Failed to write " + quadOut);
    }

    heading("Result");
    for (const std::string &m : st.messages) {
        if (m.rfind("Stopping", 0) == 0) warn(m);
    }
    if (ok && st.perMaterialRan && st.perMaterialKept) {
        std::cout << "  " << kPass << " A layout satisfying Q1-Q5 on every material region, "
                  << "each integrated out of the cross field, glued across the interfaces.\n";
    } else if (ok) {
        std::cout << "  " << kPass << " A layout satisfying Q1-Q5, from a map integrated out of "
                  << "the cross field rather than unfolded from a metric.\n";
    } else {
        std::cout << "  " << kFail << " The field-integrated map did not reach Definition 2.1.\n";
    }
    return ok ? 0 : 7;
}
