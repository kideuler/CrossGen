// Utility to run Stages 1-6 of Shepherd, Gu and Hughes (2022) on a mesh: cone
// singularities from the SIPG cross field, the cutting graph and the cut disk,
// the discrete Gauss-Bonnet check, discrete surface Ricci flow, the metric
// immersion psi_R, the subdomain labelling, and the penalty continuation
// against the layout-inducing energies.
//
//   TestMERIDIAN <mesh.obj> [options]
//   TestMERIDIAN --selftest
//
// Reports each stage and exits non-zero if the pipeline did not reach a valid
// layout. --selftest drives Stage 3 alone on prescribed cone sets that the
// models in data/meshes never produce; see selfTest() for why that is worth
// having.

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "MERIDIAN/MERIDIAN.hxx"
#include "SIPG/SIPG.hxx"
#include "TestHelper.hxx"

namespace {

const char *kPass = "\033[32m[PASS]\033[0m";
const char *kFail = "\033[31m[FAIL]\033[0m";
const char *kWarn = "\033[33m[WARN]\033[0m";

void verdict(bool ok, const std::string &what) {
    std::cout << "  " << (ok ? kPass : kFail) << " " << what << "\n";
}

void heading(const std::string &title) {
    std::cout << "\n" << title << "\n";
    std::cout << std::string(title.size(), '-') << "\n";
}

// ---------------------------------------------------------------------------
// selfTest()
//
// Stage 3 on cone sets chosen by hand rather than read off a field. The point
// is coverage the corpus cannot give: every model in data/meshes lands on a
// cone set mild enough that the metric barely moves, so the weighted-Delaunay
// flipping of Sec. 3.2.1 -- the robustness fix the paper takes from [87], and
// the one thing standing between Newton and a triangle that stops closing --
// fires once or not at all across all thirty-seven of them. A prescribed cone
// of index -4 is a cone angle of 4 pi at a single vertex, and that does
// exercise it.
//
// Each case also re-reads the answer independently of the solver's own stopping
// test: the flat cone metric is *defined* by its angle sums, so those are
// measured directly at the cones and away from them.
int selfTest() {
    struct Case {
        const char *name;
        std::shared_ptr<Mesh> mesh;
        std::vector<std::pair<int, int>> cones; // (vertex, index)
    };

    auto interiorNear = [](const Mesh &m, double x, double y) {
        int best = -1;
        double bd = std::numeric_limits<double>::infinity();
        for (int v = 0; v < static_cast<int>(m.vertices.size()); ++v) {
            if (m.isBoundaryVertex[v]) continue;
            const double d = std::hypot(m.vertices[v][0] - x, m.vertices[v][1] - y);
            if (d < bd) { bd = d; best = v; }
        }
        return best;
    };
    auto boundarySpread = [](const Mesh &m, int count, int index,
                             std::vector<std::pair<int, int>> &out) {
        const int n = static_cast<int>(m.boundaryVertices.size());
        int placed = 0;
        for (int i = 0; i < n && placed < count; ++i) {
            const int v = m.boundaryVertices[(i * 7919) % n]; // spread, not clustered
            bool clash = false;
            for (const auto &c : out) if (c.first == v) clash = true;
            if (clash) continue;
            out.push_back({v, index});
            ++placed;
        }
    };

    std::vector<Case> cases;
    {
        // Already flat: four +1 corners on a square. u must not move at all.
        Case c{"square, four +1 corners (already flat)", TestHelper::createBox(0, 0, 1, 1, 0.05), {}};
        double xs[4] = {0, 1, 1, 0}, ys[4] = {0, 0, 1, 1};
        for (int k = 0; k < 4; ++k) {
            int best = -1;
            double bd = std::numeric_limits<double>::infinity();
            for (int v : c.mesh->boundaryVertices) {
                const double d = std::hypot(c.mesh->vertices[v][0] - xs[k],
                                            c.mesh->vertices[v][1] - ys[k]);
                if (d < bd) { bd = d; best = v; }
            }
            c.cones.push_back({best, 1});
        }
        cases.push_back(c);
    }
    {
        // One interior cone of index -4: a 4 pi cone angle, eight boundary +1
        // cones to balance it. Sum: 8 - 4 = 4 = 4 chi.
        Case c{"square + interior I=-4 (valence 8)", TestHelper::createBox(0, 0, 1, 1, 0.05), {}};
        c.cones.push_back({interiorNear(*c.mesh, 0.5, 0.5), -4});
        boundarySpread(*c.mesh, 8, 1, c.cones);
        cases.push_back(c);
    }
    {
        // The clustering failure mode of Sec. 3.1: two index -2 cones a few
        // edges apart. Sum: 8 - 4 = 4.
        Case c{"disk + two adjacent I=-2 cones", TestHelper::createCircle(0, 0, 1, 0.04), {}};
        c.cones.push_back({interiorNear(*c.mesh, -0.03, 0.0), -2});
        c.cones.push_back({interiorNear(*c.mesh, 0.09, 0.0), -2});
        boundarySpread(*c.mesh, 8, 1, c.cones);
        cases.push_back(c);
    }
    {
        // Many cones at once. Sum: 8 - 4 = 4.
        Case c{"disk + eight interior I=+1 cones", TestHelper::createCircle(0, 0, 1, 0.03), {}};
        for (int k = 0; k < 8; ++k) {
            const double a = 2.0 * M_PI * k / 8.0;
            c.cones.push_back({interiorNear(*c.mesh, 0.55 * std::cos(a), 0.55 * std::sin(a)), 1});
        }
        boundarySpread(*c.mesh, 4, -1, c.cones);
        cases.push_back(c);
    }

    std::cout << "MERIDIAN Stage 3 self-test (prescribed cone sets)\n";
    int failures = 0;

    for (const Case &c : cases) {
        const int nV = static_cast<int>(c.mesh->vertices.size());
        std::vector<double> Kbar(nV, 0.0);
        int sumI = 0;
        for (const auto &p : c.cones) { Kbar[p.first] = M_PI_2 * p.second; sumI += p.second; }
        const int chi = nV - static_cast<int>(c.mesh->edges.size()) +
                        static_cast<int>(c.mesh->triangles.size());

        heading(c.name);
        std::cout << "  " << nV << " vertices, chi = " << chi
                  << ", sum I = " << sumI << " (need " << 4 * chi << ")\n";
        if (sumI != 4 * chi) {
            std::cout << "  " << kFail << " Test case is inadmissible; the case is wrong, "
                      << "not the solver.\n";
            ++failures;
            continue;
        }

        RicciFlow flow(c.mesh, Kbar);
        flow.setTolerance(1e-8);
        flow.setMaxIterations(200);
        const bool ok = flow.solve();
        const RicciFlow::Report &r = flow.getReport();

        // Independent read of the answer: a flat cone metric is one whose angle
        // sums are 2pi away from the cones and 2pi - Kbar at them.
        const std::vector<double> &K = flow.getCurvature();
        double worstCone = 0.0, worstFlat = 0.0;
        for (const auto &p : c.cones) worstCone = std::max(worstCone, std::fabs(K[p.first] - Kbar[p.first]));
        for (int v = 0; v < nV; ++v) {
            if (Kbar[v] == 0.0) worstFlat = std::max(worstFlat, std::fabs(K[v]));
        }

        std::cout << "  Newton steps: " << r.newtonIterations << ", edge flips: " << r.flips << "\n";
        std::cout << "  ||K - Kbar||_inf: " << std::scientific << std::setprecision(2)
                  << r.initialError << "  ->  " << r.finalError
                  << "    cone angles off by " << worstCone
                  << ", flat vertices off by " << worstFlat << std::defaultfloat << "\n";

        const double hess = flow.checkHessian(40);
        verdict(ok, "Converged to the flat cone metric");
        verdict(worstCone < 1e-8, "Cone angle sums are 2pi - (pi/2) I");
        verdict(worstFlat < 1e-8, "Every non-cone vertex is flat");
        verdict(r.nonRealisableFaces == 0, "Every face closes");
        verdict(r.minInteriorEdgeWeight >= -1e-9, "Weighted-Delaunay (interior w_ij >= 0)");
        verdict(hess < 1e-4, "Hessian matches a finite difference of K");

        if (!(ok && worstCone < 1e-8 && worstFlat < 1e-8 && r.nonRealisableFaces == 0 &&
              r.minInteriorEdgeWeight >= -1e-9 && hess < 1e-4)) {
            ++failures;
        }
        for (const std::string &m : r.messages) std::cout << "  " << kWarn << " " << m << "\n";
    }

    heading("Result");
    std::cout << "  " << (failures == 0 ? kPass : kFail) << " " << (cases.size() - failures)
              << " of " << cases.size() << " self-test case(s) passed\n";
    return failures == 0 ? 0 : 6;
}

void usage(const char *argv0) {
    std::cerr << "Usage: " << argv0 << " <mesh.obj> [options]\n"
              << "       " << argv0 << " --selftest\n"
              << "Options:\n"
              << "  --gamma <g>        SIPG penalty parameter        (default 10)\n"
              << "  --steps <n>        max MBO iterations            (default 500)\n"
              << "  --tol <t>          Ricci ||K - Kbar||_inf target (default 1e-8)\n"
              << "  --newton <n>       max Newton iterations         (default 100)\n"
              << "  --no-flips         disable weighted-Delaunay flipping\n"
              << "  --no-rebalance     do not repair Eq. (4) automatically\n"
              << "  --cut <file.obj>   write the cut disk Omega\n"
              << "  --cones <n>        list at most n cones          (default 20)\n"
              << "  --outer <n>        penalty continuation steps    (default 10)\n"
              << "  --inner <n>        inner iterations per step     (default 40)\n"
              << "  --lambda <l>       initial lambda_2..lambda_5    (default 1e-2)\n"
              << "  --growth <g>       lambda growth per outer step  (default 10)\n"
              << "  --near-miss <f>    Gamma_topo seeding tolerance      (default 0.15)\n"
              << "  --no-topo          skip the Gamma_topo seeding (E5 off)\n"
              << "  --no-layout        stop after Stage 4\n"
              << "  --psi <file.obj>   write psi_R, the Stage 4 immersion\n"
              << "  --layout <file.obj> write Psi, the Stage 6 result\n";
}

} // namespace

int main(int argc, char **argv) {
    if (argc < 2) { usage(argv[0]); return 1; }
    if (std::string(argv[1]) == "--selftest") return selfTest();

    const std::string path = argv[1];
    MERIDIAN::Options opts;
    std::string cutOut, psiOut, layoutOut;
    int coneListLimit = 20;

    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--gamma" && i + 1 < argc)        opts.sipgGamma = std::stod(argv[++i]);
        else if (a == "--steps" && i + 1 < argc)   opts.sipgMaxSteps = std::stoi(argv[++i]);
        else if (a == "--tol" && i + 1 < argc)     opts.ricciTolerance = std::stod(argv[++i]);
        else if (a == "--newton" && i + 1 < argc)  opts.ricciMaxIterations = std::stoi(argv[++i]);
        else if (a == "--no-flips")                opts.delaunayFlips = false;
        else if (a == "--no-rebalance")            opts.autoRebalance = false;
        else if (a == "--cut" && i + 1 < argc)     cutOut = argv[++i];
        else if (a == "--cones" && i + 1 < argc)   coneListLimit = std::stoi(argv[++i]);
        else if (a == "--outer" && i + 1 < argc)   opts.outerSteps = std::stoi(argv[++i]);
        else if (a == "--inner" && i + 1 < argc)   opts.innerIterations = std::stoi(argv[++i]);
        else if (a == "--lambda" && i + 1 < argc)  opts.lambdaInit = std::stod(argv[++i]);
        else if (a == "--growth" && i + 1 < argc)  opts.lambdaGrowth = std::stod(argv[++i]);
        else if (a == "--no-topo")                 opts.seedTopoConstraints = false;
        else if (a == "--near-miss" && i + 1 < argc) opts.topoNearMiss = std::stod(argv[++i]);
        else if (a == "--no-layout")               opts.runLayout = false;
        else if (a == "--psi" && i + 1 < argc)     psiOut = argv[++i];
        else if (a == "--layout" && i + 1 < argc)  layoutOut = argv[++i];
        else { std::cerr << "Unknown option: " << a << "\n"; usage(argv[0]); return 1; }
    }

    std::shared_ptr<Mesh> mesh;
    try {
        mesh = std::make_shared<Mesh>(path);
    } catch (const std::exception &e) {
        std::cout << kFail << " Failed to load mesh: " << e.what() << "\n";
        return 2;
    }

    std::cout << "MERIDIAN -- Shepherd, Gu and Hughes (2022), Stages 1-6\n";
    std::cout << "Mesh: " << path << "\n";
    std::cout << "  " << mesh->vertices.size() << " vertices, "
              << mesh->edges.size() << " edges, "
              << mesh->triangles.size() << " triangles, "
              << mesh->boundaryEdges.size() << " boundary edges\n";

    MERIDIAN pipeline(mesh, opts);
    bool ok = false;
    try {
        ok = pipeline.run();
    } catch (const std::exception &e) {
        std::cout << kFail << " Pipeline threw: " << e.what() << "\n";
        return 3;
    }
    const MERIDIAN::Status &st = pipeline.getStatus();

    // ---------------------------------------------------------------------
    // Stage 0 -- cross field
    // ---------------------------------------------------------------------
    heading("Stage 0  Cross field (p=0 DG/SIPG MBO)");
    std::cout << "  MBO steps: " << st.mboSteps
              << ", residual " << std::scientific << std::setprecision(3)
              << pipeline.getField().error << std::defaultfloat << "\n";
    verdict(st.fieldConverged, "Field converged");

    // ---------------------------------------------------------------------
    // Stage 1 -- cone singularities and Gauss-Bonnet
    // ---------------------------------------------------------------------
    heading("Stage 1  Cone singularities (Sec. 3.1)");
    const ConeSingularities &cones = pipeline.getCones();
    const auto gb = cones.gaussBonnet();

    std::cout << "  V - E + F = " << gb.V << " - " << gb.E << " + " << gb.F
              << " = " << gb.eulerCharacteristic
              << "   (" << gb.boundaryLoops << " boundary loop"
              << (gb.boundaryLoops == 1 ? "" : "s");
    if (gb.isolatedVertices > 0) {
        std::cout << ", " << gb.isolatedVertices << " isolated vertices excluded";
    }
    std::cout << ")\n";
    std::cout << "  Cones: " << cones.interiorCones().size() << " interior, "
              << cones.boundaryCones().size() << " boundary\n";
    if (gb.rebalanceUnits != 0) {
        std::cout << "  Rebalanced " << gb.rebalanceUnits << " index unit(s) at a cost of "
                  << std::fixed << std::setprecision(3) << gb.rebalanceCost
                  << std::defaultfloat << " quarter turns\n";
    }

    if (coneListLimit > 0) {
        const auto &list = cones.getCones();
        std::cout << "     vertex     I   valence  where      raw\n";
        int shown = 0;
        for (const auto &c : list) {
            if (shown++ >= coneListLimit) {
                std::cout << "     ... " << (list.size() - coneListLimit) << " more\n";
                break;
            }
            std::cout << "  " << std::setw(9) << c.vertex
                      << std::setw(6) << std::showpos << c.index << std::noshowpos
                      << std::setw(8) << c.valence
                      << "  " << (c.onBoundary ? "boundary" : "interior")
                      << std::setw(9) << std::fixed << std::setprecision(3) << c.raw
                      << std::defaultfloat
                      << (c.fromRebalance ? "  (rebalanced)" : "") << "\n";
        }
    }

    heading("         Discrete Gauss-Bonnet");
    std::cout << "  Eq. (4)  sum I(v)  = " << std::showpos << gb.indexSum << std::noshowpos
              << "   4 chi = " << std::showpos << gb.indexTarget << std::noshowpos
              << "   (unrounded " << std::fixed << std::setprecision(6) << gb.rawIndexSum
              << std::defaultfloat << ")\n";
    std::cout << "  Eq. (7)  sum K_v   = " << std::fixed << std::setprecision(9)
              << gb.curvatureSum << "   2 pi chi = " << gb.curvatureTarget
              << "   residual " << std::scientific << std::setprecision(3)
              << gb.curvatureResidual << std::defaultfloat << "\n";
    std::cout << "           sum Kbar  = " << std::fixed << std::setprecision(9)
              << gb.targetCurvatureSum << std::defaultfloat << "\n";
    verdict(gb.metricConsistent, "Input metric satisfies Eq. (7)");
    verdict(gb.admissible, "Cone set is admissible, Eq. (4): sum I(v) = 4 chi(S)");

    // The pipeline collects every stage's notes into one list, tagged by stage;
    // the later stages print their own below, so only the untagged ones (the
    // field, the rebalance) and Stage 1's belong here.
    for (const std::string &m : st.messages) {
        if (m.rfind("Stage 2: ", 0) == 0 || m.rfind("Stage 3: ", 0) == 0) continue;
        std::cout << "  " << kWarn << " " << m << "\n";
    }

    if (!st.conesAdmissible) {
        std::cout << "\n" << kFail << " Stopped before Ricci flow: the flat cone metric "
                  << "asked for does not exist.\n";
        return 4;
    }

    // ---------------------------------------------------------------------
    // Stage 2 -- cutting graph
    // ---------------------------------------------------------------------
    heading("Stage 2  Cutting graph and cut disk (Sec. 3.2.2)");
    const ConeCut &cut = pipeline.getCut();
    const ConeCut::Report &cr = cut.getReport();
    std::cout << "  Voids: " << cr.voids << ", HarmonicCut arcs: " << cr.harmonicCuts
              << ", cone arcs: " << cr.conesRouted << " of " << cr.interiorCones << "\n";
    std::cout << "  Omega: " << cut.getCutMesh().vertices.size() << " vertices, "
              << "chi = " << cr.eulerCharacteristic << ", "
              << cr.boundaryComponents << " boundary component"
              << (cr.boundaryComponents == 1 ? "" : "s") << ", "
              << cr.triangleComponents << " triangle component"
              << (cr.triangleComponents == 1 ? "" : "s") << "\n";

    double totalCutLength = 0.0;
    for (const auto &c : cut.getVoidCuts()) totalCutLength += c.length;
    for (const auto &c : cut.getConePaths()) totalCutLength += c.length;
    std::cout << "  Cutting graph: " << cut.getCutEdges().size() << " edges, total length "
              << std::fixed << std::setprecision(4) << totalCutLength
              << std::defaultfloat << "\n";

    verdict(cr.isDisk, "Omega = S - G is a topological disk");
    verdict(cr.allConesOnBoundary, "P is contained in G union dS (all cones on the boundary)");
    // Not a verdict: Sec. 3.2.2 prefers cuts that stop at a cone but does not
    // require it, and Fig. 7 is a cone split into three children.
    if (cr.conesSplitByCut > 0) {
        std::cout << "  " << kWarn << " " << cr.conesSplitByCut
                  << " cone(s) split by the cutting graph rather than left as leaves "
                  << "(legal; Sec. 3.2.2 prefers otherwise)\n";
    } else {
        std::cout << "  " << kPass << " Every cone is a leaf of G (Sec. 3.2.2's preferred form)\n";
    }

    for (const std::string &m : cr.messages) std::cout << "  " << kWarn << " " << m << "\n";

    if (!cutOut.empty()) {
        if (cut.writeOBJ(cutOut)) std::cout << "  Wrote cut mesh to " << cutOut << "\n";
        else std::cout << "  " << kWarn << " Failed to write " << cutOut << "\n";
    }

    // ---------------------------------------------------------------------
    // Stage 3 -- Ricci flow
    // ---------------------------------------------------------------------
    heading("Stage 3  Discrete surface Ricci flow (Sec. 3.2.1)");
    const RicciFlow &ricci = pipeline.getRicci();
    const RicciFlow::Report &rr = ricci.getReport();

    std::cout << "  Newton steps: " << rr.newtonIterations
              << ", edge flips: " << rr.flips
              << " over " << rr.flipPasses << " pass"
              << (rr.flipPasses == 1 ? "" : "es") << "\n";
    std::cout << "  ||K - Kbar||_inf: " << std::scientific << std::setprecision(3)
              << rr.initialError << "  ->  " << rr.finalError
              << "   (tolerance " << opts.ricciTolerance << ")" << std::defaultfloat << "\n";
    std::cout << "  Inversive distance cos(phi) in ["
              << std::fixed << std::setprecision(4) << rr.minInversiveDistance << ", "
              << rr.maxInversiveDistance << "]"
              << (rr.minInversiveDistance >= 1.0
                      ? "   (all >= 1: the generalised flow of [87] is doing the work)"
                      : "")
              << std::defaultfloat << "\n";
    std::cout << "  Smallest edge weight w_ij: " << std::scientific << std::setprecision(3)
              << rr.minInteriorEdgeWeight << " interior, "
              << rr.minBoundaryEdgeWeight << " boundary" << std::defaultfloat << "\n";

    const auto &uu = ricci.getU();
    if (!uu.empty()) {
        const auto mm = std::minmax_element(uu.begin(), uu.end());
        std::cout << "  Conformal factor u in [" << std::fixed << std::setprecision(4)
                  << *mm.first << ", " << *mm.second << "]"
                  << "  (metric scaled by up to " << std::setprecision(2)
                  << std::exp(*mm.second - *mm.first) << "x)" << std::defaultfloat << "\n";
    }
    if (rr.flippedAwayOriginalEdges > 0) {
        std::cout << "  " << rr.flippedAwayOriginalEdges
                  << " input edge(s) no longer exist in the flow triangulation\n";
    }

    const double hessErr = ricci.checkHessian(24);
    std::cout << "  Hessian vs. central difference: " << std::scientific
              << std::setprecision(3) << hessErr << std::defaultfloat << " relative\n";

    verdict(rr.gaussBonnetResidual < 1e-6, "Newton system is consistent (sum Kbar = 2 pi chi)");
    verdict(rr.converged, "Flat cone metric reached: ||K - Kbar||_inf < tolerance");
    verdict(rr.nonRealisableFaces == 0, "Every face of the flow metric closes");
    // Only interior edges: a boundary edge has one power height rather than two
    // and goes negative on an obtuse triangle, with no quad to flip.
    verdict(rr.minInteriorEdgeWeight >= -1e-9,
            "Flow triangulation is weighted-Delaunay (interior w_ij >= 0)");
    verdict(hessErr < 1e-4, "Hessian matches a finite difference of K");

    for (const std::string &m : rr.messages) std::cout << "  " << kWarn << " " << m << "\n";

    if (!st.readyForImmersion) {
        heading("Result");
        std::cout << "  " << kFail
                  << " Pipeline did not reach a usable flat cone metric; see above.\n";
        return 5;
    }

    // ---------------------------------------------------------------------
    // Stage 4 -- the metric immersion psi_R
    // ---------------------------------------------------------------------
    heading("Stage 4  Metric immersion psi_R (Sec. 3.2.2)");
    const Immersion &psi = pipeline.getImmersion();
    const Immersion::Report &pr = psi.getReport();

    std::cout << "  Omega laid out: " << pr.placedVertices << " of " << pr.cutVertices
              << " vertices, image " << std::fixed << std::setprecision(3)
              << (pr.uvMax[0] - pr.uvMin[0]) << " x " << (pr.uvMax[1] - pr.uvMin[1])
              << std::defaultfloat << (pr.mirrored ? "  (reflected to positive orientation)" : "")
              << "\n";
    std::cout << "  Metric realised to " << std::scientific << std::setprecision(3)
              << pr.maxMetricResidual << " relative; path-independence gap "
              << pr.maxClosureGap << std::defaultfloat << "\n";
    std::cout << "  Angle sums off target by " << std::scientific << std::setprecision(3)
              << pr.maxConeAngleResidual << " rad at the cones, "
              << pr.maxRegularAngleResidual << " elsewhere" << std::defaultfloat << "\n";
    std::cout << "  Cutting graph: " << pr.arcs << " arc(s), " << pr.seamEdgePairs
              << " seam edge pair(s); Gamma_Hol_0..3 = " << pr.holonomyCount[0] << ", "
              << pr.holonomyCount[1] << ", " << pr.holonomyCount[2] << ", "
              << pr.holonomyCount[3] << "\n";
    std::cout << "  Transition fit: rotation snapped by at most " << std::scientific
              << std::setprecision(3) << pr.maxSnapError << " rad, residual "
              << pr.maxFitResidual << std::defaultfloat << "\n";
    if (pr.alignedBoundaryEdge >= 0) {
        std::cout << "  Rotated by " << std::fixed << std::setprecision(4) << pr.globalRotation
                  << " rad to put boundary edge " << pr.alignedBoundaryEdge
                  << " on an axis" << std::defaultfloat << "\n";
    }

    verdict(pr.unplacedVertices == 0, "Every vertex of Omega was placed");
    verdict(pr.maxMetricResidual < 1e-6, "The unfolding realises the flat metric");
    verdict(pr.maxClosureGap < 1e-6,
            "The unfolding is path-independent (no cone enclosed by a loop)");
    verdict(pr.flippedFaces == 0, "Q1: no face of psi_R is inverted");
    verdict(pr.maxConeAngleResidual < 1e-6, "Q2: cone angle sums are multiples of pi/2");
    verdict(pr.degenerateArcs == 0 && pr.maxSnapError < 1e-6,
            "Q4: every arc's transition is exactly a quarter-turn rotation");
    if (pr.clusteredConePairs > 0) {
        std::cout << "  " << kWarn << " " << pr.clusteredConePairs
                  << " clustered cone pair(s); closest " << std::scientific
                  << std::setprecision(3) << pr.minConeSeparation
                  << " of the model apart" << std::defaultfloat << "\n";
    }
    for (const std::string &m : pr.messages) std::cout << "  " << kWarn << " " << m << "\n";

    if (!psiOut.empty()) {
        if (psi.writeOBJ(psiOut)) std::cout << "  Wrote psi_R to " << psiOut << "\n";
        else std::cout << "  " << kWarn << " Failed to write " << psiOut << "\n";
    }

    if (!pipeline.hasLabels()) {
        heading("Result");
        const bool immOk = pr.valid;
        std::cout << "  " << (immOk ? kPass : kFail) << " "
                  << (immOk ? "psi_R is a valid metric immersion (Stages 5 and 6 not run)."
                            : "psi_R is not usable; the later stages were not run.")
                  << "\n";
        return immOk ? 0 : 5;
    }

    // ---------------------------------------------------------------------
    // Stage 5 -- subdomain labelling
    // ---------------------------------------------------------------------
    heading("Stage 5  Subdomain labelling (Sec. 3.3)");
    const SubdomainLabels &lab = pipeline.getLabels();
    const SubdomainLabels::Report &sr = lab.getReport();

    std::cout << "  dS: " << sr.boundaryEdges << " edge(s) -> Gamma_u " << sr.boundaryEdgesU
              << ", Gamma_v " << sr.boundaryEdgesV << ", in " << sr.boundaryChains
              << " chain(s)";
    if (sr.ambiguousBoundaryEdges > 0) {
        std::cout << "  (" << sr.ambiguousBoundaryEdges << " near-tied)";
    }
    std::cout << "\n";
    std::cout << "  Features: " << sr.featureEdges << " edge(s) in " << sr.featureChains
              << " chain(s) -> Gamma_u^feat " << sr.featureChainsU
              << ", Gamma_v^feat " << sr.featureChainsV << "\n";
    std::cout << "  Separatrices traced: " << sr.separatrices << " -> "
              << sr.separatricesToCone << " to a cone, " << sr.separatricesToBoundary
              << " out through dS, " << sr.separatricesCapped << " unresolved\n";
    std::cout << "  Gamma_topo: " << sr.topoPaths << " path(s)";
    if (sr.topoPaths > 0) {
        std::cout << ", worst advance " << std::scientific << std::setprecision(3)
                  << sr.maxTopoResidual << std::defaultfloat;
    }
    std::cout << "   (mean cone spacing " << std::fixed << std::setprecision(4)
              << sr.meanConeSpacing << std::defaultfloat << ")\n";

    verdict(sr.boundaryEdges > 0, "Every curve of dS carries a label");
    for (const std::string &m : sr.messages) std::cout << "  " << kWarn << " " << m << "\n";

    if (!pipeline.hasLayout()) {
        heading("Result");
        std::cout << "  " << kWarn << " Stage 6 was not run.\n";
        return 0;
    }

    // ---------------------------------------------------------------------
    // Stage 6 -- the layout-inducing energies
    // ---------------------------------------------------------------------
    heading("Stage 6  Layout-inducing energies (Sec. 3.3)");
    const LayoutEnergy &lay = pipeline.getLayout();
    const LayoutEnergy::Report &er = lay.getReport();

    std::cout << "  Outer steps: " << er.outerSteps << " of " << opts.outerSteps
              << ", inner iterations " << er.innerIterations
              << ", line-search stops " << er.lineSearchFailures << "\n";
    std::cout << "  lambda_2..5 reached " << std::scientific << std::setprecision(2)
              << er.lambdaFinal[0] << ", " << er.lambdaFinal[1] << ", "
              << er.lambdaFinal[2] << ", " << er.lambdaFinal[3] << std::defaultfloat;
    if (er.referenceSwitches || er.relabels) {
        std::cout << "   (" << er.referenceSwitches << " reference switch(es), "
                  << er.relabels << " relabel(s))";
    }
    std::cout << "\n";
    // The two totals are each measured at their own lambda, so the pair is not
    // a decrease and is not meant to read as one: a rise means the penalties
    // outgrew the distortion, which is what a working continuation looks like.
    std::cout << "  Energy at lambda_init " << std::scientific << std::setprecision(4)
              << er.energyStart << ", at the end " << er.energyEnd
              << "   E1 " << er.e1 << ", E2 " << er.e2
              << ", E3 " << er.e3 << ", E4 " << er.e4 << ", E5 " << er.e5
              << std::defaultfloat << "\n";
    std::cout << "  Residuals (relative to the image extent)\n";
    std::cout << "     Q3 boundary  " << std::scientific << std::setprecision(3)
              << er.initialBoundaryResidual << "  ->  " << er.maxBoundaryResidual << "\n";
    std::cout << "     features     " << "        " << "  ->  " << er.maxFeatureResidual << "\n";
    std::cout << "     Q4 seam      " << er.initialSeamResidual << "  ->  "
              << er.maxSeamResidual << "\n";
    std::cout << "     Q5 topo      " << er.initialTopoResidual << "  ->  "
              << er.maxTopoResidual << std::defaultfloat << "\n";
    std::cout << "  det J in [" << std::fixed << std::setprecision(4) << er.minDetJ
              << ", ...], " << er.invertedTriangles << " inverted triangle(s)"
              << std::defaultfloat << "\n";
    std::cout << "  Q2 angle sums off the prescribed value by " << std::scientific
              << std::setprecision(3) << er.maxConeAngleResidual << " rad at the cones, "
              << er.maxRegularAngleResidual << " elsewhere" << std::defaultfloat;
    if (er.coneValenceChanges > 0) {
        std::cout << "; " << er.coneValenceChanges << " cone(s) changed valence";
    }
    std::cout << "\n";

    const double gradErr = lay.checkGradient(32);
    std::cout << "  Gradient vs. central difference: " << std::scientific
              << std::setprecision(3) << gradErr << std::defaultfloat << " relative\n";

    verdict(gradErr < 1e-4, "The assembled gradient matches a finite difference");
    verdict(er.injective, "Q1: det J > 0 on every triangle");
    verdict(er.maxSeamResidual < 1e-6, "Q4: seam transitions are exactly R_k");
    verdict(er.maxBoundaryResidual < 1e-6, "Q3: every curve of dS is on a coordinate line");
    verdict(er.maxFeatureResidual < 1e-6, "Feature chains are layout edges");
    verdict(sr.topoPaths == 0 || er.maxTopoResidual < 1e-6,
            "Q5: the connectivity constraints are met");
    verdict(er.anglesHeld, "Q2: every cone kept the angle Stage 1 prescribed for it");
    for (const std::string &m : er.messages) std::cout << "  " << kWarn << " " << m << "\n";

    if (!layoutOut.empty()) {
        if (lay.writeOBJ(layoutOut)) std::cout << "  Wrote Psi to " << layoutOut << "\n";
        else std::cout << "  " << kWarn << " Failed to write " << layoutOut << "\n";
    }

    // ---------------------------------------------------------------------
    heading("Result");
    if (ok) {
        std::cout << "  " << kPass
                  << " Psi satisfies Q1-Q5: a quadrilateral layout in the sense of "
                  << "Definition 2.1. Ready for Stage 7 (separatrix tracing).\n";
    } else {
        std::cout << "  " << kFail
                  << " The continuation did not reach a valid layout; see the residuals "
                  << "above.\n";
    }
    return ok ? 0 : 5;
}
