// Utility to run Stages 1-3 of Shepherd, Gu and Hughes (2022) on a mesh:
// cone singularities from the SIPG cross field, the cutting graph and the cut
// disk, the discrete Gauss-Bonnet check, and discrete surface Ricci flow.
//
//   TestMERIDIAN <mesh.obj> [options]
//   TestMERIDIAN --selftest
//
// Reports each stage and exits non-zero if the pipeline did not reach a flat
// cone metric on a disk. --selftest drives Stage 3 alone on prescribed cone
// sets that the models in data/meshes never produce; see selfTest() for why
// that is worth having.

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
              << "  --cones <n>        list at most n cones          (default 20)\n";
}

} // namespace

int main(int argc, char **argv) {
    if (argc < 2) { usage(argv[0]); return 1; }
    if (std::string(argv[1]) == "--selftest") return selfTest();

    const std::string path = argv[1];
    MERIDIAN::Options opts;
    std::string cutOut;
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
        else { std::cerr << "Unknown option: " << a << "\n"; usage(argv[0]); return 1; }
    }

    std::shared_ptr<Mesh> mesh;
    try {
        mesh = std::make_shared<Mesh>(path);
    } catch (const std::exception &e) {
        std::cout << kFail << " Failed to load mesh: " << e.what() << "\n";
        return 2;
    }

    std::cout << "MERIDIAN -- Shepherd, Gu and Hughes (2022), Stages 1-3\n";
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

    // ---------------------------------------------------------------------
    heading("Result");
    if (ok) {
        std::cout << "  " << kPass
                  << " Flat cone metric on a topological disk. Ready for Stage 4 "
                  << "(metric immersion).\n";
    } else {
        std::cout << "  " << kFail
                  << " Pipeline did not reach a usable flat cone metric; see above.\n";
    }
    return ok ? 0 : 5;
}
