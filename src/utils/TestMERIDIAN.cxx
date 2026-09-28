// Utility to run Stages 1-9 of Shepherd, Gu and Hughes (2022) on a mesh: cone
// singularities from the DualMBO cross field, the cutting graph and the cut disk,
// the discrete Gauss-Bonnet check, discrete surface Ricci flow, the metric
// immersion psi_R, the subdomain labelling, the penalty continuation against
// the layout-inducing energies, the separatrices of the layout that comes out,
// the arrangement they cut S into, and the bicubic patches fitted to it.
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
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "MERIDIAN/MERIDIAN.hxx"
#include "MERIDIAN/ParallelLDLT.hxx"
#include "dualmbo/DualMBO.hxx"
#include "mesh/QuadMesh.hxx"
#include "mesh/TMOP.hxx"
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

    // ParallelLDLT against Eigen's SimplicialLDLT, which it must reproduce to
    // the last bit on any number of threads (see its class comment). The
    // matrix has Stage 6's shape: two unknowns per vertex of a disk, a 6x6
    // block per triangle coupling them, and a handful of long-range terms
    // like E5's; the shift after the first round is innerSolve()'s Tikhonov
    // fallback, new values on the same pattern.
    {
        heading("ParallelLDLT == SimplicialLDLT, bit for bit");
        const std::shared_ptr<Mesh> disk = TestHelper::createCircle(0, 0, 1, 0.03);
        const int nV = static_cast<int>(disk->vertices.size());
        std::vector<Eigen::Triplet<double>> trips;
        for (int t = 0; t < static_cast<int>(disk->triangles.size()); ++t) {
            const Triangle &tri = disk->triangles[t];
            const Point &a = disk->vertices[tri[0]], &b = disk->vertices[tri[1]], &c = disk->vertices[tri[2]];
            const Point e[3] = {c - b, a - c, b - a};
            const double area = 0.5 * std::fabs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
            const double w[3] = {2.0 + std::sin(0.7 * t), 1.0 + 0.5 * std::cos(1.3 * t), 0.3 * std::sin(0.1 * t)};
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    const double g = (e[i][0] * e[j][0] + e[i][1] * e[j][1]) / (4.0 * area);
                    const double W[2][2] = {{w[0], w[2]}, {w[2], w[1]}};
                    for (int p = 0; p < 2; ++p)
                        for (int q = 0; q < 2; ++q)
                            trips.emplace_back(2 * tri[i] + p, 2 * tri[j] + q, g * W[p][q]);
                }
            }
        }
        for (int d = 0; d < 2 * nV; ++d) trips.emplace_back(d, d, 1e-3);
        for (int k = 0; k + 1 < static_cast<int>(disk->boundaryVertices.size()); k += 37) {
            const int u = 2 * disk->boundaryVertices[k], v = 2 * disk->boundaryVertices[(k * 7 + 11) % disk->boundaryVertices.size()];
            if (u == v) continue;
            trips.emplace_back(u, u, 5.0); trips.emplace_back(v, v, 5.0);
            trips.emplace_back(u, v, -5.0); trips.emplace_back(v, u, -5.0);
        }
        Eigen::SparseMatrix<double> A(2 * nV, 2 * nV);
        A.setFromTriplets(trips.begin(), trips.end());
        A.makeCompressed();
        Eigen::VectorXd rhs(2 * nV);
        for (int d = 0; d < 2 * nV; ++d) rhs[d] = std::sin(0.37 * d) + 0.1;

        bool allSame = true;
        for (int round = 0; round < 2; ++round) {
            if (round == 1) {
                for (int d = 0; d < 2 * nV; ++d) A.coeffRef(d, d) += 1e-6;
            }
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ref;
            ref.analyzePattern(A);
            ref.factorize(A);
            const Eigen::VectorXd xRef = ref.solve(rhs);
            const auto &LRef = ref.matrixL().nestedExpression();
            for (int threads : {1, 2, 4}) {
                ParallelLDLT par;
                par.setThreads(threads);
                par.analyzePattern(A);
                par.factorize(A);
                par.factorize(A);   // twice: the second runs on the first's state
                const Eigen::VectorXd x = par.solve(rhs);
                const auto &L = par.matrixL().nestedExpression();
                const bool same =
                    par.info() == Eigen::Success && L.nonZeros() == LRef.nonZeros() &&
                    std::equal(L.innerIndexPtr(), L.innerIndexPtr() + L.nonZeros(), LRef.innerIndexPtr()) &&
                    std::memcmp(L.valuePtr(), LRef.valuePtr(), sizeof(double) * L.nonZeros()) == 0 &&
                    std::memcmp(par.vectorD().data(), ref.vectorD().data(), sizeof(double) * ref.vectorD().size()) == 0 &&
                    std::memcmp(x.data(), xRef.data(), sizeof(double) * x.size()) == 0;
                verdict(same, std::string(round ? "shifted, " : "") + std::to_string(threads) +
                                  " thread(s) asked, " + std::to_string(par.lastThreads()) +
                                  " used: L, D and a solve identical to Eigen's (" +
                                  std::to_string(2 * nV) + " unknowns)");
                allSame = allSame && same;
            }
        }
        if (!allSame) ++failures;
    }

    heading("Result");
    const size_t total = cases.size() + 1;
    std::cout << "  " << (failures == 0 ? kPass : kFail) << " " << (total - failures)
              << " of " << total << " self-test case(s) passed\n";
    return failures == 0 ? 0 : 6;
}

void usage(const char *argv0) {
    std::cerr << "Usage: " << argv0 << " <mesh.obj> [options]\n"
              << "       " << argv0 << " --selftest\n"
              << "Options:\n"
              << "  --gamma <g>        edge penalty parameter        (default 10)\n"
              << "  --steps <n>        max MBO iterations            (default 500)\n"
              << "  --tol <t>          Ricci ||K - Kbar||_inf target (default 1e-8)\n"
              << "  --newton <n>       max Newton iterations         (default 100)\n"
              << "  --no-flips         disable weighted-Delaunay flipping\n"
              << "  --no-rebalance     do not repair Eq. (4) automatically\n"
              << "  --cut-to-graph     cone cuts may stop on an earlier cut, not only dS\n"
              << "  --cut-through-interfaces  route the cutting graph as if the material\n"
              << "                     interfaces were not there (Stage 2 before it knew)\n"
              << "  --cut <file.obj>   write the cut disk Omega\n"
              << "  --no-interfaces    ignore the material tags: no Stage 0b, no E6\n"
              << "  --no-field-interfaces  leave the Stage 0 field aligned to dS alone\n"
              << "  --no-cancel-dipoles    keep the +1/-1 cone pairs on curved interfaces\n"
              << "  --no-e6            keep the interface network but drop E6, so the\n"
              << "                     interfaces are held by E3 alone\n"
              << "  --no-prescribe     leave the cone index at an interface node to the\n"
              << "                     cross field rather than to the geometry\n"
              << "  --no-propagate     label each interface chain by its own flux rather\n"
              << "                     than by the node quantisation\n"
              << "  --no-seam-turns    give every chain of a branch the branch's label,\n"
              << "                     ignoring the quarter turns of the cuts it crosses\n"
              << "  --lag-e6           re-read E6's tangent lengths from the map between\n"
              << "                     outer steps (off; it costs more than it buys)\n"
              << "  --fit-features     mesh dS and the interfaces on their spline fits like\n"
              << "                     every other arc, rather than on the traced curves\n"
              << "  --kink <deg>       interface corner threshold        (default 45)\n"
              << "  --loop-splits <n>  arcs a closed interface is cut into (default 4)\n"
              << "  --split-disks      cut a disk rim into arcs, as any other closed loop\n"
              << "  --disk-templates   excise every circular inclusion (Stage 0c) and fill it\n"
              << "                     back in with an O-grid template (Stage 11)\n"
              << "  --disk-squareness <w>  how square the template's core is  (default 0.55)\n"
              << "  --disk-ring <n>    rows of elements in the ring, 0 = auto  (default 0)\n"
              << "  --disk-smooth <n>  smoothing sweeps over the template     (default 300)\n"
              << "  --nodes <n>        list at most n interface nodes    (default 12)\n"
              << "  --interfaces <f.obj> write the interface network as polylines\n"
              << "  --cones <n>        list at most n cones          (default 20)\n"
              << "  --outer <n>        penalty continuation steps    (default 10)\n"
              << "  --inner <n>        inner iterations per step     (default 40)\n"
              << "  --lambda <l>       initial lambda_2..lambda_5    (default 1e-2)\n"
              << "  --growth <g>       lambda growth per outer step  (default 10)\n"
              << "  --near-miss <f>    Gamma_topo seeding tolerance      (default 0.15)\n"
              << "  --retry <f>        tolerance to re-seed from when the first layout\n"
              << "                     leaves a piece of S no grid covers  (default 0.04)\n"
              << "  --no-retry         one attempt only, at --near-miss\n"
              << "  --alt-ref          let a stalled Stage 6 swap E1's reference to\n"
              << "                     the model's Euclidean metric (off; it has no\n"
              << "                     cones, so it shreds the cone one-rings)\n"
              << "  --no-topo          skip the Gamma_topo seeding (E5 off)\n"
              << "  --no-layout        stop after Stage 4\n"
              << "  --no-trace         skip Stage 7 (separatrix tracing)\n"
              << "  --snap <t>         cone snap tolerance, of the extent (default 1e-6)\n"
              << "  --snap-rings <n>   cone snap search radius, in faces (default 2)\n"
              << "  --max-steps <n>    separatrix step cap             (default 50000)\n"
              << "  --no-cycles        do not stop a curve once it is provably winding\n"
              << "  --no-self-returns  do not seed Gamma_topo for a curve back to its own cone\n"
              << "  --one-per-pair     keep only one Gamma_topo path per pair of cones\n"
              << "  --repair <n>       rounds of Sec. 3.3's repair       (default 6)\n"
              << "  --repair-max <n>   constraints per round, 0 = all    (default 4)\n"
              << "  --repair-gap <f>   near-miss window, of the extent   (default 1e-3)\n"
              << "  --curves <n>       list at most n separatrices     (default 10)\n"
              << "  --psi <file.obj>   write psi_R, the Stage 4 immersion\n"
              << "  --layout <file.obj> write Psi, the Stage 6 result\n"
              << "  --sep <file.obj>   write the separatrices on the model\n"
              << "  --sep-uv <file.obj> write the separatrices in the layout image\n"
              << "  --no-arrange       skip Stage 8 (arrangement)\n"
              << "  --merge <t>        node merge tolerance, of the model  (default 1e-7)\n"
              << "  --corner <r>       corner tolerance, radians           (default 0.35)\n"
              << "  --collapse <t>     sliver-arc tolerance, of the model  (default 1e-4)\n"
              << "  --no-trim          leave a winding separatrix's loose end where it stopped\n"
              << "  --no-arr-score     judge a repair round on the curves alone, not the layout\n"
              << "  --patience <n>     repair rounds that may fail to improve  (default 2)\n"
              << "  --arcs <file.obj>  write the layout arcs on the model\n"
              << "  --faces <file.obj> write the layout faces as polylines\n"
              << "  --no-splines       skip Stage 9 (spline reconstruction)\n"
              << "  --segments <n>     cubic segments per arc              (default 3)\n"
              << "  --spline-boundary  fit the arcs on dS instead of carrying them as the\n"
              << "                     polylines of mesh edges they are\n"
              << "  --spline-interfaces  the same for the material interface network\n"
              << "  --samples <n>      patch sampling grid                 (default 8)\n"
              << "  --fit <file.obj>   write the fitted arc curves\n"
              << "  --net <file.obj>   write the patch control nets\n"
              << "  --surf <file.obj>  write the reconstructed patches\n"
              << "  --step <file.step> write the reconstructed model as a STEP B-rep\n"
              << "  --brep <file.brep> ... and in OpenCASCADE's native BREP format\n"
              << "  --no-mesh          skip Stage 10 (quadrilateral meshing)\n"
              << "  --target <h>       target edge length, model units      (default 0.05)\n"
              << "  --collapse-span <f>  contract a chord whose every patch is thinner\n"
              << "                     than this times the target             (default 0.5)\n"
              << "  --no-collapse      keep every chord, however thin its patches\n"
              << "  --min-edges <n>    fewest edges per chord               (default 1)\n"
              << "  --max-edges <n>    most edges per chord, 0 = no cap     (default 0)\n"
              << "  --polyline-mesh    mesh the traced arcs, not the spline fits\n"
              << "  --smooth <n>       Winslow sweeps per block, 0 = off      (default 500)\n"
              << "  --smooth-below <j> smooth blocks worse than this Jacobian (default 0)\n"
              << "  --chords <n>       list at most n chords                (default 0)\n"
              << "  --tmop <n>         smooth the final mesh with n TMOP sweeps, 0 = off\n"
              << "  --tmop-power <p>   TMOP exponent: 1 optimises the average, 2 (default)\n"
              << "                     and up chase the worst element\n"
              << "  --pin-features     TMOP pins every feature node instead of sliding\n"
              << "  --mesh <file.obj>  write the quadrilateral mesh\n"
              << "  --mesh-vtu <f.vtu> write it as a VTK unstructured grid\n"
              << "  --mfem <f.mesh>   write it as an MFEM mesh, material id per element\n";
}

} // namespace

int main(int argc, char **argv) {
    if (argc < 2) { usage(argv[0]); return 1; }
    if (std::string(argv[1]) == "--selftest") return selfTest();

    const std::string path = argv[1];
    MERIDIAN::Options opts;
    std::string cutOut, psiOut, layoutOut, sepOut, sepUVOut;
    std::string arcsOut, facesOut, fitOut, netOut, surfOut, meshOut, meshVTUOut, mfemOut;
    std::string stepOut, brepOut;
    int tmopSweeps = 0;
    double tmopPower = 2.0;
    bool tmopPinFeatures = false;
    int coneListLimit = 20;
    int curveListLimit = 10;
    int chordListLimit = 0;
    int interfaceListLimit = 12;
    std::string interfacesOut;

    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--gamma" && i + 1 < argc)        opts.dualMBOGamma = std::stod(argv[++i]);
        else if (a == "--steps" && i + 1 < argc)   opts.dualMBOMaxSteps = std::stoi(argv[++i]);
        else if (a == "--tol" && i + 1 < argc)     opts.ricciTolerance = std::stod(argv[++i]);
        else if (a == "--newton" && i + 1 < argc)  opts.ricciMaxIterations = std::stoi(argv[++i]);
        else if (a == "--no-flips")                opts.delaunayFlips = false;
        else if (a == "--no-rebalance")            opts.autoRebalance = false;
        else if (a == "--cut-to-graph")            opts.coneCutsToBoundary = false;
        else if (a == "--cut-through-interfaces")  opts.coneCutInterfaceAvoidance = 0.0;
        else if (a == "--cut" && i + 1 < argc)     cutOut = argv[++i];
        else if (a == "--no-interfaces")           opts.materialInterfaces = false;
        else if (a == "--no-field-interfaces")     opts.alignFieldToInterfaces = false;
        else if (a == "--no-cancel-dipoles")       opts.cancelInterfaceDipoles = false;
        else if (a == "--no-e6")                   opts.interfaceCorners = false;
        else if (a == "--no-prescribe")            opts.prescribeInterfaceCones = false;
        else if (a == "--no-propagate")            opts.propagateInterfaceLabels = false;
        else if (a == "--no-seam-turns")           opts.seamTurnInterfaceLabels = false;
        else if (a == "--lag-e6")                  opts.lagInterfaceScales = true;
        else if (a == "--fit-features")            opts.quadFeaturesOnTracedArcs = false;
        else if (a == "--kink" && i + 1 < argc)
            opts.interfaceKinkAngle = std::stod(argv[++i]) * M_PI / 180.0;
        else if (a == "--split-disks")             opts.splitCircleLoops = true;
        else if (a == "--disk-templates")          opts.diskTemplates = true;
        else if (a == "--disk-squareness" && i + 1 < argc)
            opts.diskCoreSquareness = std::stod(argv[++i]);
        else if (a == "--disk-ring" && i + 1 < argc)
            opts.diskRingDepth = std::stoi(argv[++i]);
        else if (a == "--disk-smooth" && i + 1 < argc)
            opts.diskSmoothingPasses = std::stoi(argv[++i]);
        else if (a == "--loop-splits" && i + 1 < argc)
            opts.interfaceLoopSplits = std::stoi(argv[++i]);
        else if (a == "--nodes" && i + 1 < argc)   interfaceListLimit = std::stoi(argv[++i]);
        else if (a == "--interfaces" && i + 1 < argc) interfacesOut = argv[++i];
        else if (a == "--cones" && i + 1 < argc)   coneListLimit = std::stoi(argv[++i]);
        else if (a == "--outer" && i + 1 < argc)   opts.outerSteps = std::stoi(argv[++i]);
        else if (a == "--inner" && i + 1 < argc)   opts.innerIterations = std::stoi(argv[++i]);
        else if (a == "--lambda" && i + 1 < argc)  opts.lambdaInit = std::stod(argv[++i]);
        else if (a == "--growth" && i + 1 < argc)  opts.lambdaGrowth = std::stod(argv[++i]);
        else if (a == "--alt-ref")                 opts.alternateReference = true;
        else if (a == "--no-topo")                 opts.seedTopoConstraints = false;
        else if (a == "--no-retry")                opts.topoNearMissRetry = 0.0;
        else if (a == "--retry" && i + 1 < argc)   opts.topoNearMissRetry = std::stod(argv[++i]);
        else if (a == "--near-miss" && i + 1 < argc) opts.topoNearMiss = std::stod(argv[++i]);
        else if (a == "--no-layout")               opts.runLayout = false;
        else if (a == "--no-trace")                opts.runSeparatrices = false;
        else if (a == "--snap" && i + 1 < argc)    opts.separatrixSnap = std::stod(argv[++i]);
        else if (a == "--snap-rings" && i + 1 < argc) opts.separatrixSnapRings = std::stoi(argv[++i]);
        else if (a == "--max-steps" && i + 1 < argc) opts.separatrixMaxSteps = std::stoi(argv[++i]);
        else if (a == "--no-cycles")               opts.separatrixDetectCycles = false;
        else if (a == "--no-self-returns")         opts.seedSelfReturns = false;
        else if (a == "--one-per-pair")            opts.seedAllConnections = false;
        else if (a == "--repair" && i + 1 < argc)  opts.repairPasses = std::stoi(argv[++i]);
        else if (a == "--repair-max" && i + 1 < argc) opts.repairMaxPerPass = std::stoi(argv[++i]);
        else if (a == "--repair-gap" && i + 1 < argc) opts.repairGapLimit = std::stod(argv[++i]);
        else if (a == "--curves" && i + 1 < argc)  curveListLimit = std::stoi(argv[++i]);
        else if (a == "--psi" && i + 1 < argc)     psiOut = argv[++i];
        else if (a == "--layout" && i + 1 < argc)  layoutOut = argv[++i];
        else if (a == "--sep" && i + 1 < argc)     sepOut = argv[++i];
        else if (a == "--sep-uv" && i + 1 < argc)  sepUVOut = argv[++i];
        else if (a == "--no-arrange")              opts.runArrangement = false;
        else if (a == "--merge" && i + 1 < argc)   opts.arrangementMerge = std::stod(argv[++i]);
        else if (a == "--corner" && i + 1 < argc)  opts.arrangementCorner = std::stod(argv[++i]);
        else if (a == "--collapse" && i + 1 < argc) opts.arrangementCollapse = std::stod(argv[++i]);
        else if (a == "--no-trim")                 opts.arrangementTrim = false;
        else if (a == "--no-arr-score")            opts.repairScoresArrangement = false;
        else if (a == "--patience" && i + 1 < argc) opts.repairPatience = std::stoi(argv[++i]);
        else if (a == "--arcs" && i + 1 < argc)    arcsOut = argv[++i];
        else if (a == "--faces" && i + 1 < argc)   facesOut = argv[++i];
        else if (a == "--no-splines")              opts.runSplines = false;
        else if (a == "--segments" && i + 1 < argc) opts.splineSegments = std::stoi(argv[++i]);
        else if (a == "--spline-boundary")         opts.fitBoundaryArcs = true;
        else if (a == "--spline-interfaces")       opts.fitInterfaceArcs = true;
        else if (a == "--samples" && i + 1 < argc) opts.splineSamples = std::stoi(argv[++i]);
        else if (a == "--fit" && i + 1 < argc)     fitOut = argv[++i];
        else if (a == "--net" && i + 1 < argc)     netOut = argv[++i];
        else if (a == "--surf" && i + 1 < argc)    surfOut = argv[++i];
        else if (a == "--step" && i + 1 < argc)    stepOut = argv[++i];
        else if (a == "--brep" && i + 1 < argc)    brepOut = argv[++i];
        else if (a == "--no-mesh")                 opts.runQuadMesh = false;
        else if (a == "--target" && i + 1 < argc)  opts.quadTargetEdge = std::stod(argv[++i]);
        else if (a == "--collapse-span" && i+1 < argc) opts.quadCollapseSpan = std::stod(argv[++i]);
        else if (a == "--no-collapse")             opts.quadCollapseSpan = 0.0;
        else if (a == "--min-edges" && i + 1 < argc) opts.quadMinIntervals = std::stoi(argv[++i]);
        else if (a == "--max-edges" && i + 1 < argc) opts.quadMaxIntervals = std::stoi(argv[++i]);
        else if (a == "--polyline-mesh")           opts.quadUseSplines = false;
        else if (a == "--smooth" && i + 1 < argc)  opts.quadSmoothingPasses = std::stoi(argv[++i]);
        else if (a == "--smooth-below" && i + 1 < argc) opts.quadSmoothingThreshold = std::stod(argv[++i]);
        else if (a == "--chords" && i + 1 < argc)  chordListLimit = std::stoi(argv[++i]);
        else if (a == "--tmop" && i + 1 < argc)    tmopSweeps = std::stoi(argv[++i]);
        else if (a == "--tmop-power" && i + 1 < argc) tmopPower = std::stod(argv[++i]);
        else if (a == "--pin-features")            tmopPinFeatures = true;
        else if (a == "--mesh" && i + 1 < argc)    meshOut = argv[++i];
        else if (a == "--mesh-vtu" && i + 1 < argc) meshVTUOut = argv[++i];
        else if (a == "--mfem" && i + 1 < argc)    mfemOut = argv[++i];
        else { std::cerr << "Unknown option: " << a << "\n"; usage(argv[0]); return 1; }
    }

    std::shared_ptr<Mesh> mesh;
    try {
        mesh = std::make_shared<Mesh>(path);
    } catch (const std::exception &e) {
        std::cout << kFail << " Failed to load mesh: " << e.what() << "\n";
        return 2;
    }

    std::cout << "MERIDIAN -- Shepherd, Gu and Hughes (2022), Stages 1-9, then meshing\n";
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
    heading("Stage 0  Cross field (p=0 dual-mesh MBO)");
    std::cout << "  MBO steps: " << st.mboSteps
              << ", residual " << std::scientific << std::setprecision(3)
              << pipeline.getField().error << std::defaultfloat << "\n";
    verdict(st.fieldConverged, "Field converged");

    // ---------------------------------------------------------------------
    // Stage 0b -- the material interface network
    // ---------------------------------------------------------------------
    if (pipeline.hasInterfaces() && pipeline.getInterfaces().multiMaterial()) {
        heading("Stage 0b  Material interfaces (feature lines from the material tags)");
        const Interfaces &itf = pipeline.getInterfaces();
        const Interfaces::Report &fr = itf.getReport();

        std::cout << "  " << fr.materials << " material(s) in " << fr.regions
                  << " connected region(s); " << fr.interfaceEdges
                  << " interface edge(s) in " << fr.branches << " branch(es)\n";
        std::cout << "  Nodes: " << fr.nodes << " -- " << fr.junctions << " junction, "
                  << fr.landings << " on dS, " << fr.kinks << " kink, "
                  << fr.loopSplits << " loop split";
        if (fr.dangling) std::cout << ", " << fr.dangling << " dangling";
        std::cout << "\n";
        std::cout << "  Sectors: " << fr.illPosedNodes << " of " << fr.nodes
                  << " node(s) are more than "
                  << std::fixed << std::setprecision(0) << (10.0)
                  << " degrees from whole quarter turns; the worst is "
                  << std::setprecision(1) << (fr.worstSectorResidual * 180.0 / M_PI)
                  << " degrees" << std::defaultfloat;
        if (fr.worstSectorNode >= 0 && fr.worstSectorNode < static_cast<int>(itf.nodes().size())) {
            std::cout << " (at vertex " << itf.nodes()[fr.worstSectorNode].vertex << ")";
        }
        std::cout << "\n";
        std::cout << "  Cones prescribed from the geometry: " << fr.prescribedCones
                  << ", summing to " << std::showpos << fr.prescribedIndexSum << std::noshowpos
                  << "; the field was " << st.interfacePrescriptionShift
                  << " index unit(s) away from them\n";
        std::cout << "  Region balance: " << fr.regionsBalanced << " of " << fr.regions
                  << " satisfy sum I + sum (2 - q) = 4 chi";
        if (fr.quartersMoved || fr.cornersInserted) {
            std::cout << ", after moving " << fr.quartersMoved << " quarter turn(s) and "
                      << "inserting " << fr.cornersInserted << " corner(s)";
        }
        std::cout << "\n";

        if (interfaceListLimit > 0) {
            std::cout << "     vertex  kind       I   sectors (deg -> quarters)\n";
            int shown = 0;
            for (const Interfaces::Node &n : itf.nodes()) {
                if (shown++ >= interfaceListLimit) {
                    std::cout << "     ... " << (itf.nodes().size() - shown + 1) << " more\n";
                    break;
                }
                const char *kind = n.kind == Interfaces::NodeKind::Junction  ? "junction"
                                 : n.kind == Interfaces::NodeKind::Landing   ? "landing "
                                 : n.kind == Interfaces::NodeKind::Kink      ? "kink    "
                                 : n.kind == Interfaces::NodeKind::LoopSplit ? "loopsplit"
                                 : n.kind == Interfaces::NodeKind::Balance   ? "balance "
                                                                             : "dangling";
                std::cout << "     " << std::setw(6) << n.vertex << "  " << kind << " "
                          << std::showpos << std::setw(3) << n.index << std::noshowpos << "   ";
                for (size_t k = 0; k < n.sector.size(); ++k) {
                    std::cout << std::fixed << std::setprecision(0)
                              << (n.sector[k] * 180.0 / M_PI) << "->" << n.quarters[k]
                              << (k + 1 < n.sector.size() ? ", " : "") << std::defaultfloat;
                }
                if (n.rebalanced) std::cout << "  (rebalanced)";
                if (!n.wellPosed) std::cout << "  [ill-posed]";
                std::cout << "\n";
            }
        }

        verdict(fr.dangling == 0, "Every interface branch runs between two nodes");
        verdict(fr.regionsBalanced == fr.regions,
                "Every material region can be a whole number of quadrilaterals");
        if (fr.illPosedNodes > 0) {
            std::cout << "  " << kWarn << " " << fr.illPosedNodes
                      << " node(s) are ill-posed for a vertex-based cross field: no single "
                      << "cross is tangent to all the interfaces meeting there. The layout "
                      << "absorbs the turn, the field cannot represent it.\n";
        }
        if (!interfacesOut.empty()) {
            const bool w = itf.writeOBJ(interfacesOut);
            std::cout << "  " << (w ? "Wrote " : "Could not write ") << interfacesOut << "\n";
        }
    }

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
    // every stage from the second on prints its own below, so only the untagged
    // ones (the field, the rebalance) and Stage 1's belong here.
    for (const std::string &m : st.messages) {
        bool later = false;
        for (int s = 2; s <= 9; ++s) {
            if (m.rfind("Stage " + std::to_string(s) + ": ", 0) == 0) later = true;
        }
        if (later) continue;
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
              << ", cone arcs: " << cr.conesRouted << " of " << cr.interiorCones
              << " (" << cr.conesToBoundary << " reached dS";
    if (cr.conesFellBack > 0) std::cout << ", " << cr.conesFellBack << " fell back";
    std::cout << ")\n";
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

    std::cout << "  Junctions of G away from dS: " << cr.interiorJunctions
              << (opts.coneCutsToBoundary ? "  (expected 0: every arc runs to dS)" : "")
              << "\n";

    // What the cut has in common with the material interface network. Nothing
    // is what E3 and E6 need of it; see ConeCut's class comment.
    if (st.materials > 1) {
        std::cout << "  Interface network: " << cr.interfaceVertsOnCut
                  << " vertex/vertices of G on it (" << cr.interfaceNodesOnCut
                  << " node(s)), " << cr.interfaceEdgesOnCut << " arc(s) along it";
        if (cr.conesOnInterfaceNodes > 0) {
            std::cout << "; " << cr.conesOnInterfaceNodes
                      << " cone(s) sit on a node and must leave from one";
        }
        std::cout << "\n";
        const Mesh &om = cut.getOriginalMesh();
        for (const auto &cp : cut.getConePaths()) {
            if (cp.interfaceVerts == 0) continue;
            const Point &a = om.vertices[cp.cone];
            std::cout << "     cone " << cp.cone << " (" << std::fixed << std::setprecision(4)
                      << a[0] << ", " << a[1] << ")" << std::defaultfloat << " I = " << cp.index
                      << ": " << cp.interfaceVerts << " interface vertex/vertices ("
                      << cp.interfaceNodes << " node(s))"
                      << (cp.startsOnNode ? ", the cone itself among them" : "")
                      << " on a " << cp.path.size() << "-vertex arc\n";
        }
    }

    verdict(cr.isDisk, "Omega = S - G is a topological disk");
    if (st.materials > 1) {
        // A verdict on the half of it that is always reachable. Running an arc
        // *along* an interface is never forced -- a single transversal crossing
        // is cheaper than two vertices of a branch under any triangulation --
        // and it is the contact that costs E3 a whole chain, so it is a
        // failure. Crossing is the other half, and a cone inside an inclusion
        // has no interface-free route to dS at all, so that one is counted on
        // the line above and warned about below rather than failed.
        verdict(cr.interfaceEdgesOnCut == 0, "No arc of G runs along a material interface");
    }
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
    if (sr.interfaceCorners > 0 || sr.featureChainsPropagated > 0) {
        std::cout << "  Interfaces: " << sr.featureChainsPropagated
                  << " chain(s) labelled from the node quantisation ("
                  << sr.featureLabelsCorrected << " where the flux said otherwise), "
                  << sr.interfaceCorners << " sector(s) for E6";
        if (sr.featureChainsPastSeam > 0 || sr.branchesUnwalked > 0) {
            std::cout << "; " << sr.featureChainsPastSeam << " chain(s) past a crossing of G, "
                      << sr.featureChainsSeamFlipped << " turned an odd number of quarters by it";
            if (sr.branchesUnwalked > 0) {
                std::cout << ", " << sr.branchesUnwalked << " branch(es) not walked";
            }
        }
        if (sr.interfaceCornersSpanningCut > 0) {
            std::cout << "; " << sr.interfaceCornersSpanningCut
                      << " dropped where the cutting graph runs through the sector";
        }
        std::cout << "\n";
    }
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
    std::cout << "  Near-miss tolerance " << st.topoNearMissUsed
              << (st.topoNearMissRetried ? "  (after a retry)" : "") << "\n";

    verdict(sr.boundaryEdges > 0, "Every curve of dS carries a label");
    for (const std::string &m : sr.messages) std::cout << "  " << kWarn << " " << m << "\n";
    // Stage 5 and 6 notes the *repair loop* logged through MERIDIAN::Status.
    // Every other stage prints its own Report, but traceAndRepair has none of
    // its own -- it writes through a callback -- so without this its rollbacks
    // and its pruning happen silently.
    for (const std::string &m : st.messages) {
        if (m.rfind("Stage 5: ", 0) == 0 || m.rfind("Stage 6: ", 0) == 0) {
            std::cout << "  " << kWarn << " " << m << "\n";
        }
    }

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
    std::cout << "  lambda_2..6 reached " << std::scientific << std::setprecision(2)
              << er.lambdaFinal[0] << ", " << er.lambdaFinal[1] << ", "
              << er.lambdaFinal[2] << ", " << er.lambdaFinal[3] << ", "
              << er.lambdaFinal[4] << std::defaultfloat;
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
              << ", E6 " << er.e6 << std::defaultfloat << "\n";
    std::cout << "  Residuals (relative to the image extent)\n";
    std::cout << "     Q3 boundary  " << std::scientific << std::setprecision(3)
              << er.initialBoundaryResidual << "  ->  " << er.maxBoundaryResidual << "\n";
    std::cout << "     features     " << "        " << "  ->  " << er.maxFeatureResidual << "\n";
    std::cout << "     Q4 seam      " << er.initialSeamResidual << "  ->  "
              << er.maxSeamResidual << "\n";
    std::cout << "     Q5 topo      " << er.initialTopoResidual << "  ->  "
              << er.maxTopoResidual << std::defaultfloat << "\n";
    if (er.interfaceCorners > 0) {
        std::cout << "  E6: the worst interface sector is " << std::scientific
                  << std::setprecision(3) << er.maxInterfaceResidual
                  << " rad from the quarter turns the model has there"
                  << std::defaultfloat;
        if (er.interfaceCornerChanges > 0) {
            std::cout << "; " << er.interfaceCornerChanges
                      << " turned a different whole number of quarters";
        }
        std::cout << "\n";
    }
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
    if (er.interfaceCorners > 0) {
        verdict(er.interfaceCornerChanges == 0 && er.maxInterfaceResidual < 1e-3,
                "E6: every interface sector turns the quarters the model has");
    }
    verdict(sr.topoPaths == 0 || er.maxTopoResidual < 1e-6,
            "Q5: the connectivity constraints are met");
    verdict(er.anglesHeld, "Q2: every cone kept the angle Stage 1 prescribed for it");
    for (const std::string &m : er.messages) std::cout << "  " << kWarn << " " << m << "\n";

    if (!layoutOut.empty()) {
        if (lay.writeOBJ(layoutOut)) std::cout << "  Wrote Psi to " << layoutOut << "\n";
        else std::cout << "  " << kWarn << " Failed to write " << layoutOut << "\n";
    }

    if (!pipeline.hasSeparatrices()) {
        heading("Result");
        std::cout << "  " << (ok ? kPass : kFail) << " "
                  << (ok ? "Psi satisfies Q1-Q5 (Stage 7 not run)."
                         : "The continuation did not reach a valid layout; see above.")
                  << "\n";
        return ok ? 0 : 5;
    }

    // ---------------------------------------------------------------------
    // Stage 7 -- the separatrices of Psi
    // ---------------------------------------------------------------------
    heading("Stage 7  Separatrix tracing (Sec. 4, Q5 of Definition 2.1)");
    const Separatrices &sep = pipeline.getSeparatrices();
    const Separatrices::Report &tr = sep.getReport();

    std::cout << "  Emitted " << tr.emitted << " separatrix/ces from " << tr.cones
              << " cone(s); the indices prescribe " << tr.prescribed
              << " (4 - I interior, 1 - I on dS)\n";
    if (tr.designatedCone >= 0) {
        std::cout << "  Footnote 3: no cones, so vertex " << tr.designatedCone
                  << " was called the surface's only singularity\n";
    }
    std::cout << "  Ends: " << tr.endedAtCone << " at a cone, " << tr.endedAtBoundary
              << " out through dS, " << tr.cycled << " on a closed orbit, " << tr.capped
              << " at the step cap";
    if (tr.stuck > 0 || tr.degenerate > 0) {
        std::cout << ", " << tr.stuck << " stuck, " << tr.degenerate << " degenerate";
    }
    std::cout << "\n";
    if (tr.nearMisses > 0 || tr.grazes > 0) {
        std::cout << "  Near misses within " << std::scientific << std::setprecision(2)
                  << tr.nearMissWindow << " of the image extent: " << tr.nearMisses
                  << " unterminated, " << tr.grazes
                  << " that grazed a cone and left through dS anyway" << std::defaultfloat
                  << "\n";
    }
    std::cout << "  " << tr.triangleSteps << " triangle crossing(s), longest curve "
              << tr.maxTriangleSteps << "; " << tr.seamCrossings
              << " seam crossing(s), most on one curve " << tr.maxSeamCrossings << "\n";
    // Both gaps are lengths in the image, so they are reported against its
    // extent -- the same scale Stage 6 measured its constraint residuals in.
    std::cout << "  Snap tolerance " << std::scientific << std::setprecision(2)
              << opts.separatrixSnap << " of the image; worst snap actually taken "
              << (tr.extent > 0.0 ? tr.maxSnapGap / tr.extent : 0.0) << std::defaultfloat << "\n";
    if (tr.worstMissCurve >= 0) {
        std::cout << "  Worst unresolved curve passed " << std::scientific
                  << std::setprecision(2)
                  << (tr.extent > 0.0 ? tr.worstMissGap / tr.extent : 0.0)
                  << " of the image from the nearest cone it crossed"
                  << std::defaultfloat << "\n";
    }
    std::cout << "  Pullback continuity across the cuts: " << std::scientific
              << std::setprecision(3) << tr.maxPullbackGap << " of the model"
              << std::defaultfloat << "\n";
    std::cout << "  Fan sweep: cone angles off the prescribed value by " << std::scientific
              << std::setprecision(3) << tr.maxConeAngleResidual << " rad"
              << std::defaultfloat;
    if (tr.worstConeAngleSlot >= 0 &&
        tr.worstConeAngleSlot < static_cast<int>(sep.coneVertices().size())) {
        std::cout << " (worst at vertex " << sep.coneVertices()[tr.worstConeAngleSlot] << ")";
    }
    std::cout << "\n";

    if (curveListLimit > 0 && !sep.curves().empty()) {
        std::cout << "      from     I  dir   tris  seams  ends            gap\n";
        int shown = 0;
        for (const Separatrices::Curve &c : sep.curves()) {
            if (shown++ >= curveListLimit) {
                std::cout << "     ... " << (sep.curves().size() - curveListLimit) << " more\n";
                break;
            }
            static const char *kDir[4] = {"+u", "+v", "-u", "-v"};
            std::cout << "  " << std::setw(9) << sep.coneVertices()[c.cone]
                      << std::setw(6) << std::showpos << c.index << std::noshowpos
                      << std::setw(5) << kDir[c.dir]
                      << std::setw(7) << c.steps.size()
                      << std::setw(7) << c.seamCrossings << "  ";
            switch (c.end) {
                case Separatrices::End::Cone:
                    std::cout << "cone " << std::setw(9) << std::left
                              << sep.coneVertices()[c.toCone] << std::right;
                    break;
                case Separatrices::End::Boundary:
                    std::cout << std::setw(14) << std::left << "dS" << std::right;
                    break;
                case Separatrices::End::Capped:
                    std::cout << std::setw(14) << std::left << "step cap" << std::right;
                    break;
                case Separatrices::End::Cycle:
                    std::cout << std::setw(14) << std::left << "closed orbit" << std::right;
                    break;
                case Separatrices::End::Stuck:
                    std::cout << std::setw(14) << std::left << "stuck" << std::right;
                    break;
                default:
                    std::cout << std::setw(14) << std::left << "degenerate" << std::right;
                    break;
            }
            const double g = (c.end == Separatrices::End::Cone) ? c.gap : c.nearestConeGap;
            if (std::isfinite(g) && tr.extent > 0.0) {
                std::cout << std::scientific << std::setprecision(1) << g / tr.extent
                          << std::defaultfloat;
            } else {
                std::cout << "      -";
            }
            std::cout << "\n";
        }
    }

    verdict(tr.fanFailures == 0, "Every cone's one-ring fan in Omega was swept");
    verdict(tr.emitted == tr.prescribed,
            "Every cone emitted the number of separatrices its index prescribes");
    verdict(tr.stuck == 0 && tr.degenerate == 0, "Every ray was traceable");
    verdict(tr.maxPullbackGap < 1e-9,
            "Every curve is continuous on S across the cutting graph");
    // Not a verdict. A separatrix that runs on past the cap is Q5 failing to
    // hold on a curve E5 was never given a constraint for, and the remedy is
    // upstream of this stage: away from the cones the flat metric makes the
    // integral curves geodesics of a translation surface, and a geodesic in a
    // direction nothing has quantised does not close -- it fills the model and
    // leaves through dS eventually or not at all. Remark 3.1 is this same fact
    // stated as patch counts.
    if (tr.capped + tr.cycled == 0) {
        std::cout << "  " << kPass
                  << " Q5: every separatrix terminates at a cone or leaves through dS\n";
    } else {
        std::cout << "  " << kWarn << " " << (tr.capped + tr.cycled) << " of " << tr.emitted
                  << " separatrices terminated at neither: Q5 holds on the Gamma_topo paths "
                  << "E5 was given, not on every integral curve\n";
    }
    if (tr.grazes > 0) {
        std::cout << "  " << kWarn << " " << tr.grazes
                  << " separatrix/ces passed within " << std::scientific << std::setprecision(2)
                  << tr.nearMissWindow << " of the image extent of a cone and left through dS "
                  << "instead of stopping at it" << std::defaultfloat
                  << ": Q5 is satisfied but each is a quadrilateral of poor aspect ratio "
                  << "(Remark 3.1), and a connectivity constraint is what closes it\n";
    }
    for (const std::string &m : tr.messages) std::cout << "  " << kWarn << " " << m << "\n";

    if (!sepOut.empty()) {
        if (sep.writeOBJ(sepOut, Separatrices::Space::Model)) {
            std::cout << "  Wrote the separatrices on the model to " << sepOut << "\n";
        } else {
            std::cout << "  " << kWarn << " Failed to write " << sepOut << "\n";
        }
    }
    if (!sepUVOut.empty()) {
        if (sep.writeOBJ(sepUVOut, Separatrices::Space::Image)) {
            std::cout << "  Wrote the separatrices in the image to " << sepUVOut << "\n";
        } else {
            std::cout << "  " << kWarn << " Failed to write " << sepUVOut << "\n";
        }
    }

    if (!pipeline.hasArrangement()) {
        heading("Result");
        std::cout << "  " << (ok ? kPass : kFail) << " "
                  << (ok ? "Psi satisfies Q1-Q5 (Stage 8 not run)."
                         : "The continuation did not reach a valid layout; see above.")
                  << "\n";
        return ok ? 0 : 5;
    }

    // ---------------------------------------------------------------------
    // Stage 8 -- the arrangement, and the layout read off it
    // ---------------------------------------------------------------------
    heading("Stage 8  Arrangement and layout extraction (Sec. 4)");
    const Arrangement &arrangement = pipeline.getArrangement();
    const Arrangement::Report &arep = arrangement.getReport();

    std::cout << "  Nodes: " << arep.nodes << " -- " << arep.coneNodes << " cone, "
              << arep.crossingNodes << " separatrix crossing, " << arep.boundaryHitNodes
              << " on dS, " << arep.boundaryCornerNodes << " boundary corner";
    if (arep.interfaceNodes > 0) {
        std::cout << ", " << arep.interfaceNodes << " interface node";
    }
    if (arep.interfaceHitNodes > 0) {
        std::cout << ", " << arep.interfaceHitNodes << " separatrix on an interface";
    }
    if (arep.danglingNodes > 0) std::cout << ", " << arep.danglingNodes << " dangling";
    if (arep.mergedNodes > 0) std::cout << "  (" << arep.mergedNodes << " coincided)";
    std::cout << "\n";
    // Not a defect, and worth saying plainly: a cone-to-cone arc is a
    // separatrix of both its cones, so Stage 7 traced every one of them twice.
    std::cout << "  Curves traced from both ends: " << arep.duplicateCurves
              << " second copy/ies dropped\n";
    std::cout << "  Arcs: " << arep.arcs << " -- " << arep.separatrixArcs << " separatrix, "
              << arep.boundaryArcs << " on dS";
    if (arep.interfaceArcs > 0) std::cout << ", " << arep.interfaceArcs << " on an interface";
    if (arep.degenerateArcs > 0) std::cout << ", " << arep.degenerateArcs << " degenerate";
    if (arep.danglingArcs > 0) std::cout << ", " << arep.danglingArcs << " dangling";
    std::cout << "\n";
    std::cout << "  Faces: " << arep.faces << " -- " << arep.patches << " inside S, of which "
              << arep.quads << " have four corners and " << arep.simpleQuads
              << " one arc on each side\n";
    std::cout << "  Patch area from " << std::scientific << std::setprecision(3)
              << arep.minPatchArea << " to " << arep.maxPatchArea << "; they cover "
              << std::fixed << std::setprecision(6) << arep.areaCoverage << " of S"
              << std::defaultfloat << "\n";
    std::cout << "  Sectors off a whole quarter turn by at most " << std::scientific
              << std::setprecision(3) << arep.maxSectorResidual << " rad";
    if (arep.worstSectorNode >= 0) std::cout << " (at node " << arep.worstSectorNode << ")";
    if (arep.ambiguousSectors > 0) {
        std::cout << "; " << arep.ambiguousSectors << " farther out than "
                  << opts.arrangementCorner << " rad";
    }
    std::cout << std::defaultfloat << "\n";
    if (arep.selfCrossings > 0 || arep.parallelOverlaps > 0 || arep.collapsedArcs > 0 ||
        arep.sliversKept > 0 || arep.clusteredCones > 0) {
        std::cout << "  " << arep.collapsedArcs << " sliver arc(s) collapsed, "
                  << arep.sliversKept << " left as they were, "
                  << arep.selfCrossings << " self-crossing(s), " << arep.parallelOverlaps
                  << " segment pair(s) too nearly collinear to cut";
        if (arep.clusteredCones > 0) {
            std::cout << ", " << arep.clusteredCones << " clustered cone pair(s) left alone";
        }
        std::cout << "\n";
    }
    if (arep.trimmedCurves > 0) {
        std::cout << "  " << arep.trimmedCurves << " separatrix/ces that Q5 did not close were "
                  << "ended at the first layout edge they met, as T-junctions rather than "
                  << "loose ends\n";
    }
    if (arep.featureChains > 0) {
        std::cout << "  Feature chains: " << arep.featureChainsCovered << " of "
                  << arep.featureChains << " are a union of arcs; worst gap "
                  << std::scientific << std::setprecision(3) << arep.maxFeatureGap
                  << " of the model" << std::defaultfloat << "\n";
    }
    if (arep.maxConeSnapGap > 0.0) {
        std::cout << "  Curve ends moved onto their cone by at most " << std::scientific
                  << std::setprecision(3) << arep.maxConeSnapGap << " of the model"
                  << std::defaultfloat << "\n";
    }
    if (arep.interpolatedAngles > 0) {
        std::cout << "  " << kWarn << " " << arep.interpolatedAngles
                  << " arc end(s) reach their cone from outside its one-ring, so the fan "
                  << "placed them by interpolation rather than exactly\n";
    }

    verdict(arep.danglingNodes == 0, "Every curve of the layout ends at a node");
    verdict(arep.wrongCornerFaces == 0, "Every patch has exactly four corners");
    verdict(arep.tJunctions == 0, "No T-junction: every node a patch passes turns a quarter");
    verdict(arep.simpleQuads == arep.patches, "Every patch has one arc on each side");
    verdict(arep.unsharedArcs == 0, "Every arc separates two patches, or a patch from dS");
    if (st.materials > 1) {
        verdict(arep.mixedPatches == 0, "Every patch lies inside one material");
    }
    verdict(arep.coneValenceErrors == 0 && arep.isolatedCones == 0,
            "Every cone has the arcs its index prescribes");
    verdict(arep.areaCoverage > 0.999 && arep.areaCoverage < 1.001, "The patches tile S");
    if (arep.featureChains > 0) {
        verdict(arep.featureChainsCovered == arep.featureChains,
                "Every feature chain is a union of arcs");
    }
    for (const std::string &m : arep.messages) std::cout << "  " << kWarn << " " << m << "\n";

    // The same layout, as the class ATLAS's Stage 5 hands back too
    // (mesh/BlockDecomposition.hxx) -- checked against Arrangement's own
    // counts here, so a driver run that already exercises Stage 8 on real
    // models is also what exercises the shared representation, rather than
    // leaving that to the viewer alone.
    {
        const BlockDecomposition D = arrangement.blockDecomposition(
            pipeline.hasSplines() ? &pipeline.getSplines() : nullptr);
        bool sidesOk = true;
        for (size_t b = 0; b < D.blocks.size() && sidesOk; ++b) {
            for (int s = 0; s < 4; ++s) {
                if (D.blocks[b].edges[s] < 0 || D.sidePolyline(static_cast<int>(b), s).size() < 2) {
                    sidesOk = false;
                    break;
                }
            }
        }
        verdict(static_cast<int>(D.blocks.size()) == arep.simpleQuads,
                "blockDecomposition() carries every simple quad patch as a block");
        verdict(sidesOk, "... with all four sides valid macro edges");
    }

    if (!arcsOut.empty()) {
        if (arrangement.writeOBJ(arcsOut)) {
            std::cout << "  Wrote the layout arcs to " << arcsOut << "\n";
        } else {
            std::cout << "  " << kWarn << " Failed to write " << arcsOut << "\n";
        }
    }
    if (!facesOut.empty()) {
        if (arrangement.writePatchOBJ(facesOut)) {
            std::cout << "  Wrote the patch outlines to " << facesOut << "\n";
        } else {
            std::cout << "  " << kWarn << " Failed to write " << facesOut << "\n";
        }
    }

    if (!pipeline.hasSplines()) {
        heading("Result");
        std::cout << "  " << (ok ? kPass : kFail) << " "
                  << (ok ? "Psi satisfies Q1-Q5 (Stage 9 not run)."
                         : "The continuation did not reach a valid layout; see above.")
                  << "\n";
        return ok ? 0 : 5;
    }

    // ---------------------------------------------------------------------
    // Stage 9 -- the spline reconstruction
    // ---------------------------------------------------------------------
    heading("Stage 9  Spline reconstruction (Sec. 5)");
    const SplineFit &fit = pipeline.getSplines();
    const SplineFit::Report &sf = fit.getReport();

    std::cout << "  Fitted " << sf.fittedArcs << " of " << sf.arcs << " arc(s) to cubic "
              << "B-splines: " << sf.controlPointsPerArc << " control points each, "
              << fit.getOptions().segments << " Bezier segment(s)\n";
    std::cout << "  Deviation from the traced curves: " << std::scientific
              << std::setprecision(3) << sf.maxDeviation << " worst, " << sf.rmsDeviation
              << " rms, of the model" << std::defaultfloat;
    if (sf.worstArc >= 0) std::cout << " (arc " << sf.worstArc << ")";
    std::cout << "\n";
    if (sf.exactArcs > 0) {
        std::cout << "  Carried exactly, as the polylines the mesh has them as: "
                  << sf.exactArcs << " arc(s) on dS or an interface; "
                  << sf.correctedPatches << " patch(es) blend one\n";
        std::cout << "  Their control nets alone would have been out by " << std::scientific
                  << std::setprecision(3) << sf.maxNetDeviation << " of the model"
                  << std::defaultfloat;
        if (sf.worstNetArc >= 0) std::cout << " (arc " << sf.worstNetArc << ")";
        std::cout << ", which is what building those patches on the polylines takes out\n";
    }
    std::cout << "  Patches: " << sf.patches << " of " << sf.faces << " face(s), "
              << sf.controlPointsPerArc << " x " << sf.controlPointsPerArc
              << " control points each";
    if (sf.skipped > 0) std::cout << "; " << sf.skipped << " skipped";
    std::cout << "\n";
    std::cout << "  Watertightness: " << sf.sharedArcs
              << " arc(s) shared by two patches, control points apart by "
              << std::scientific << std::setprecision(3) << sf.maxSeamGap
              << ", corners off their node by " << sf.maxCornerGap
              << ", patch boundaries off their arc by " << sf.maxBoundaryGap
              << std::defaultfloat << "\n";
    std::cout << "  B-rep: " << sf.brepFaces << " face(s), " << sf.brepEdges << " edge(s) -- "
              << sf.brepSharedEdges << " shared by two faces, " << sf.brepFreeEdges
              << " bounding one; the kernel calls it " << (sf.brepValid ? "valid" : "INVALID")
              << "\n";
    std::cout << "  Worst sampled cell is " << std::fixed << std::setprecision(4)
              << sf.minCellRatio << " of the mean; " << sf.foldedPatches
              << " folded patch(es)" << std::defaultfloat << "\n";
    std::cout << "  Patch area " << std::scientific << std::setprecision(6) << sf.patchArea
              << " against " << sf.faceArea << " for the same faces of the arrangement"
              << std::defaultfloat << "\n";

    verdict(sf.curves == sf.arcs, "Every arc came out with a curve");
    verdict(sf.underdetermined == 0, "Every arc had more sample points than control points");
    verdict(sf.skipped == 0, "Every patch of the layout became a Coons patch");
    verdict(sf.maxBoundaryGap <= SplineFit::Report::boundaryTolerance,
            "Every patch lies on the four arcs it was built from");
    verdict(sf.watertight, "Watertight: both patches on an arc carry the same control points");
    verdict(sf.brepFaces == sf.patches && sf.brepSharedEdges == sf.sharedArcs && sf.brepValid,
            "The B-rep is closed: every patch is a face, every shared arc one edge of two");
    verdict(sf.foldedPatches == 0, "No patch folds");
    for (const std::string &m : sf.messages) std::cout << "  " << kWarn << " " << m << "\n";

    if (!fitOut.empty()) {
        if (fit.writeCurvesOBJ(fitOut)) {
            std::cout << "  Wrote the fitted curves to " << fitOut << "\n";
        } else {
            std::cout << "  " << kWarn << " Failed to write " << fitOut << "\n";
        }
    }
    if (!netOut.empty()) {
        if (fit.writeNetOBJ(netOut)) {
            std::cout << "  Wrote the control nets to " << netOut << "\n";
        } else {
            std::cout << "  " << kWarn << " Failed to write " << netOut << "\n";
        }
    }
    if (!surfOut.empty()) {
        if (fit.writeSurfaceOBJ(surfOut)) {
            std::cout << "  Wrote the patches to " << surfOut << "\n";
        } else {
            std::cout << "  " << kWarn << " Failed to write " << surfOut << "\n";
        }
    }
    if (!stepOut.empty()) {
        if (fit.writeSTEP(stepOut)) {
            std::cout << "  Wrote the B-rep to " << stepOut << "\n";
        } else {
            std::cout << "  " << kWarn << " Failed to write " << stepOut << "\n";
        }
    }
    if (!brepOut.empty()) {
        if (fit.writeBREP(brepOut)) {
            std::cout << "  Wrote the B-rep to " << brepOut << "\n";
        } else {
            std::cout << "  " << kWarn << " Failed to write " << brepOut << "\n";
        }
    }

    if (!pipeline.hasQuadMesh()) {
        heading("Result");
        std::cout << "  " << (ok ? kPass : kFail) << " "
                  << (ok ? "Psi satisfies Q1-Q5 (Stage 10 not run)."
                         : "The continuation did not reach a valid layout; see above.")
                  << "\n";
        return ok ? 0 : 5;
    }

    // ---------------------------------------------------------------------
    // Stage 10 -- the quadrilateral mesh
    // ---------------------------------------------------------------------
    heading("Stage 10  Quadrilateral meshing (chords and interval assignment)");
    const QuadMesh &qm = pipeline.getQuadMesh();
    const QuadMesh::Report &qr = qm.getReport();

    std::cout << "  Target edge length " << std::fixed << std::setprecision(4) << qr.target
              << " on a model " << qr.modelExtent << " across" << std::defaultfloat << "\n";
    std::cout << "  Chords: " << qr.chords << " over " << qr.arcsAssigned
              << " arc(s); " << qr.minIntervals << " to " << qr.maxIntervals
              << " edge(s) each, mean " << std::fixed << std::setprecision(2)
              << qr.meanIntervals << std::defaultfloat;
    if (qr.clampedChords > 0) {
        std::cout << "  (" << qr.clampedChords << " held at the floor or the cap)";
    }
    std::cout << "\n";
    std::cout << "  Mesh: " << qr.vertices << " vertices, " << qr.quads
              << " quadrilateral(s) in " << qr.blocks << " structured block(s)";
    if (qr.unmeshedPatches > 0) {
        std::cout << "; " << qr.unmeshedPatches << " patch(es) left unmeshed, "
                  << std::fixed << std::setprecision(2) << (100.0 * qr.unmeshedArea)
                  << "% of S" << std::defaultfloat;
    }
    std::cout << "\n";
    if (qr.collapsedChords > 0) {
        std::cout << "  Contracted " << qr.collapsedChords << " chord(s), which merged "
                  << qr.collapsedPatches << " face(s) -- " << std::fixed
                  << std::setprecision(2) << (100.0 * qr.collapsedArea)
                  << "% of S -- into their neighbours on " << qr.weldedVertices
                  << " welded vertex/vertices" << std::defaultfloat << "\n";
    }
    // The spread is the price of the integer constraint: a chord that runs
    // through patches of different sizes has one count for all of them, so an
    // arc at either end of its length range is cut into edges away from the
    // target. The rms log ratio is what the assignment minimised.
    std::cout << "  Edge length in [" << std::fixed << std::setprecision(4) << qr.minEdge
              << ", " << qr.maxEdge << "], mean " << qr.meanEdge
              << "; worst is " << std::setprecision(2) << qr.worstEdgeRatio
              << "x the target, rms log ratio " << std::setprecision(3) << qr.edgeRatioRms
              << std::defaultfloat << "\n";
    std::cout << "  Scaled Jacobian " << std::fixed << std::setprecision(4)
              << qr.minScaledJacobian << " worst, " << qr.meanScaledJacobian << " mean; "
              << qr.invertedQuads << " inverted element(s)" << std::defaultfloat << "\n";
    if (opts.quadSmoothingPasses > 0 && qr.smoothedBlocks > 0) {
        std::cout << "  Winslow smoothing: " << qr.smoothedBlocks << " block(s), at most "
                  << qr.smoothingSweeps << " sweep(s) on one; before it the worst scaled "
                  << "Jacobian was " << std::fixed << std::setprecision(4)
                  << qr.minScaledJacobianBefore << " with " << qr.invertedBefore
                  << " inverted element(s)" << std::defaultfloat << "\n";
    }
    if (qr.reflexCorners > 0) {
        // Not a defect of this stage and worth saying so where it is reported:
        // a layout corner the model turns more than a half turn at reverses the
        // element there whatever the grid does, and the remedy is upstream.
        std::cout << "  " << kWarn << " " << qr.reflexCorners
                  << " corner(s) of the layout turn through more than a half turn on the "
                  << "model; a structured grid reverses its corner element at each of them "
                  << "and no smoothing can reach it\n";
    }
    std::cout << "  Element area from " << std::scientific << std::setprecision(3)
              << qr.minQuadArea << " to " << qr.maxQuadArea << "; the mesh covers "
              << std::fixed << std::setprecision(6)
              << (qr.patchArea > 0.0 ? qr.meshArea / qr.patchArea : 0.0)
              << " of the patches it was built on" << std::defaultfloat << "\n";
    std::cout << "  Edges: " << qr.interiorEdges << " shared by two elements, "
              << qr.boundaryEdges << " on the boundary of the mesh";
    if (qr.nonManifoldEdges > 0) std::cout << ", " << qr.nonManifoldEdges << " used by more";
    if (qr.cracks > 0) std::cout << "; " << qr.cracks << " coincident vertex pair(s)";
    std::cout << "\n";
    if (st.materials > 1) {
        std::cout << "  Materials: " << qr.materials << " on the elements, "
                  << qr.interfaceEdges << " element edge(s) shared by two of them, "
                  << qr.mixedQuads << " element(s) straddling an interface";
        if (qr.unlocatedQuads > 0) {
            std::cout << " (" << qr.unlocatedQuads
                      << " took theirs from the nearest triangle)";
        }
        std::cout << "\n";
    }

    if (chordListLimit > 0 && !qm.chords().empty()) {
        std::cout << "     chord   arcs  edges   ideal    shortest arc   longest arc\n";
        int shown = 0;
        for (size_t c = 0; c < qm.chords().size(); ++c) {
            if (shown++ >= chordListLimit) {
                std::cout << "     ... " << (qm.chords().size() - chordListLimit) << " more\n";
                break;
            }
            const QuadMesh::Chord &ch = qm.chords()[c];
            std::cout << "  " << std::setw(8) << c
                      << std::setw(7) << ch.arcs.size()
                      << std::setw(7) << ch.intervals
                      << std::setw(9) << std::fixed << std::setprecision(2) << ch.idealIntervals
                      << std::setw(12) << std::setprecision(4) << ch.minLength
                      << std::setw(14) << ch.maxLength << std::defaultfloat
                      << (ch.clamped ? "  (clamped)" : "") << "\n";
        }
    }

    verdict(qr.quads > 0, "The layout produced a mesh");
    verdict(qr.nonManifoldEdges == 0,
            "Conforming: every edge is shared by two elements or bounds the mesh");
    verdict(qr.cracks == 0, "Watertight: no two vertices sit at the same point");
    verdict(qr.invertedQuads == 0, "No element is inverted");
    verdict(qr.unmeshedPatches == 0, "Every patch of the layout was meshed");
    if (st.materials > 1) {
        verdict(qr.mixedQuads == 0, "No element straddles a material interface");
    }
    for (const std::string &m : qr.messages) std::cout << "  " << kWarn << " " << m << "\n";

    // ---------------------------------------------------------------------
    // Stage 11 -- the O-grid templates
    // ---------------------------------------------------------------------
    bool templatesOK = true;
    if (pipeline.hasDiskTemplate()) {
        heading("Stage 11  O-grid templates on the excised inclusions");
        const DiskTemplate &dt = pipeline.getDiskTemplate();
        const DiskTemplate::Report &dr = dt.getReport();
        templatesOK = dr.valid;

        std::cout << "  " << dr.filled << " of " << dr.inclusions
                  << " inclusion(s) templated: " << dr.blocks << " block(s), "
                  << dr.quads << " element(s), " << dr.vertices << " new vertex/vertices";
        if (dr.smoothingSweeps > 0) {
            std::cout << "; at most " << dr.smoothingSweeps << " smoothing sweep(s) on one";
        }
        std::cout << "\n";
        if (dr.refusedOdd || dr.refusedShort || dr.refusedOpen) {
            std::cout << "  Refused: " << dr.refusedOdd << " rim(s) with an odd edge count, "
                      << dr.refusedShort << " too short, " << dr.refusedOpen
                      << " that Stage 10 did not close\n";
        }
        std::cout << "  Merged mesh: " << dr.mergedVertices << " vertices, "
                  << dr.mergedQuads << " quadrilateral(s); scaled Jacobian "
                  << std::fixed << std::setprecision(4) << dr.minScaledJacobian
                  << " worst, " << dr.meanScaledJacobian << " mean ("
                  << dr.templateMinScaledJacobian << " worst on the templates)"
                  << std::defaultfloat << "\n";
        std::cout << "  Edge length in [" << std::fixed << std::setprecision(4) << dr.minEdge
                  << ", " << dr.maxEdge << "]" << std::defaultfloat << "\n";
        std::cout << "  Edges: " << dr.interiorEdges << " shared by two elements, "
                  << dr.boundaryEdges << " on the boundary of the mesh\n";

        verdict(dr.filled == dr.inclusions, "Every circular inclusion was templated");
        verdict(dr.nonManifoldEdges == 0,
                "Conforming: every edge is shared by two elements or bounds the mesh");
        verdict(dr.cracks == 0, "Watertight: the templates share the rim rather than repeat it");
        verdict(dr.invertedQuads == 0, "No element of the merged mesh is inverted");
        for (const std::string &m : dr.messages) std::cout << "  " << kWarn << " " << m << "\n";
    }

    // ---------------------------------------------------------------------
    // Whatever came last -- the Stage 10 mesh, or the matrix plus its filled
    // inclusions -- is adopted into the standalone quad mesh class of
    // src/mesh/QuadMesh.hxx. That is the form the .mesh writer and any
    // smoother work on, and its topology is rebuilt from the elements alone,
    // so what it reports here is a second opinion on the same mesh rather
    // than a restatement of the pipeline's own bookkeeping.
    heading("Final mesh (mesh::QuadMesh)");
    mesh::QuadMesh::Options finalOpts;
    finalOpts.fixAllFeatureNodes = tmopPinFeatures;
    mesh::QuadMesh finalMesh = pipeline.hasDiskTemplate()
                                   ? mesh::QuadMesh::from(pipeline.getDiskTemplate(), finalOpts)
                                   : mesh::QuadMesh::from(qm, finalOpts);
    {
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
    }

    // ---------------------------------------------------------------------
    // Stage 12 -- TMOP. Node-local Newton on the TMOP energy over the mesh
    // above, with the boundary and the material interfaces free to slide along
    // themselves. Everything written below this point is the smoothed mesh, so
    // --mesh, --mesh-vtu and --mfem all agree with what is reported here.
    bool tmopOK = true;
    if (tmopSweeps > 0) {
        heading("Stage 12  TMOP smoothing (mesh::TMOP)");
        mesh::TMOP::Options topt;
        topt.metric = mesh::TMOP::ShapeSize007;
        topt.maxSweeps = tmopSweeps;
        topt.exponent = tmopPower;
        mesh::TMOP smoother(finalMesh, topt);
        tmopOK = smoother.run();
        const mesh::TMOP::Report &tr2 = smoother.getReport();
        std::cout << "  " << tr2.sweeps << " sweep(s)";
        if (tr2.untangleSweeps > 0)
            std::cout << " after " << tr2.untangleSweeps << " untangling sweep(s)";
        std::cout << ", " << tr2.colors << " colour(s), " << tr2.threads << " thread(s)"
                  << (tr2.openMP ? "" : " (no OpenMP)") << ", " << std::fixed
                  << std::setprecision(3) << tr2.seconds << " s\n";
        std::cout << "  Nodes: " << tr2.freeNodes << " free, " << tr2.slidingNodes
                  << " sliding, " << tr2.fixedNodes << " fixed\n";
        std::cout << "  Energy per unit target area: " << std::scientific
                  << std::setprecision(3) << tr2.energyBefore << " -> " << tr2.energyAfter
                  << std::defaultfloat << "\n";
        std::cout << "  Scaled Jacobian: worst " << std::fixed << std::setprecision(4)
                  << tr2.minScaledJacobianBefore << " -> " << tr2.minScaledJacobianAfter
                  << ", mean " << tr2.meanScaledJacobianBefore << " -> "
                  << tr2.meanScaledJacobianAfter << "\n";
        std::cout << "  Worst aspect ratio: " << tr2.worstAspectBefore << " -> "
                  << tr2.worstAspectAfter << std::defaultfloat << "\n";
        for (const std::string &m : tr2.messages) std::cout << "  " << kWarn << " " << m << "\n";
        verdict(tr2.invertedAfter == 0, "No element of the smoothed mesh is inverted");
        verdict(tr2.minScaledJacobianAfter >= tr2.minScaledJacobianBefore - 1e-12,
                "The worst element is no worse than it was");
        // Sliding moves nodes along the chord through their two feature
        // neighbours, which leaves the triangle those neighbours span -- and so
        // the polygon's area -- unchanged. A difference here means the domain
        // itself moved, which no amount of smoothing is worth.
        verdict(std::fabs(tr2.areaAfter - tr2.areaBefore) <= 1e-9 * std::fabs(tr2.areaBefore),
                "The domain has exactly the area it started with");
        if (!tr2.converged)
            std::cout << "  " << kWarn << " stopped on the sweep cap, still moving "
                      << std::scientific << std::setprecision(2) << tr2.lastSweepMove
                      << " mean edge(s) per sweep" << std::defaultfloat << "\n";
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
                std::cout << "  " << kWarn
                          << " material ids were shifted so the MFEM attribute is positive\n";
            if (mr.unusedVertices > 0)
                std::cout << "  " << kWarn << " " << mr.unusedVertices
                          << " vertex/vertices no element uses were dropped\n";
        } else {
            std::cout << "  " << kWarn << " Failed to write " << mfemOut << "\n";
        }
    }

    // Without --tmop these go out through the pipeline's own writers, which
    // carry a little more than the adopted mesh does; with it, the smoothed
    // positions are only in finalMesh, so that is what has to be written.
    if (!meshOut.empty()) {
        const bool wrote = tmopSweeps > 0
                               ? finalMesh.writeOBJ(meshOut)
                               : (pipeline.hasDiskTemplate()
                                      ? pipeline.getDiskTemplate().writeOBJ(meshOut)
                                      : qm.writeOBJ(meshOut));
        if (wrote) std::cout << "  Wrote the quad mesh to " << meshOut << "\n";
        else std::cout << "  " << kWarn << " Failed to write " << meshOut << "\n";
    }
    if (!meshVTUOut.empty()) {
        const bool wrote = tmopSweeps > 0
                               ? finalMesh.writeVTU(meshVTUOut)
                               : (pipeline.hasDiskTemplate()
                                      ? pipeline.getDiskTemplate().writeVTU(meshVTUOut)
                                      : qm.writeVTU(meshVTUOut));
        if (wrote) std::cout << "  Wrote the quad mesh to " << meshVTUOut << "\n";
        else std::cout << "  " << kWarn << " Failed to write " << meshVTUOut << "\n";
    }

    // ---------------------------------------------------------------------
    heading("Result");
    if (ok && tr.valid && arep.valid && sf.valid && qr.valid && templatesOK && tmopOK) {
        std::cout << "  " << kPass << " " << sf.patches
                  << " watertight bicubic patch(es), C0 across their shared curves and C2 "
                  << "inside, from a layout satisfying Q1-Q5, meshed into " << qr.quads
                  << " conforming quadrilateral(s) at a target edge length of " << qr.target;
        if (pipeline.hasDiskTemplate()) {
            const DiskTemplate::Report &dr = pipeline.getDiskTemplate().getReport();
            std::cout << ", plus " << dr.filled << " templated inclusion(s) for a merged "
                      << "mesh of " << dr.mergedQuads << " element(s)";
        }
        if (tmopSweeps > 0)
            std::cout << ", smoothed to a worst scaled Jacobian of " << std::fixed
                      << std::setprecision(4) << finalMesh.quality.minScaledJacobian
                      << std::defaultfloat;
        std::cout << ".\n";
        return 0;
    }
    if (ok && tr.valid && arep.valid && sf.valid && qr.valid && templatesOK && !tmopOK) {
        std::cout << "  " << kWarn << " The layout and the mesh are clean, but TMOP left "
                  << "the mesh no better than it found it; see Stage 12 above.\n";
        return 5;
    }
    if (ok && tr.valid && arep.valid && sf.valid && qr.valid && !templatesOK) {
        const DiskTemplate::Report &dr = pipeline.getDiskTemplate().getReport();
        std::cout << "  " << kWarn << " The layout of the excised model is clean -- "
                  << sf.patches << " watertight patch(es), " << qr.quads
                  << " conforming element(s) -- but " << (dr.inclusions - dr.filled)
                  << " of " << dr.inclusions << " inclusion(s) could not be templated.\n";
        return 4;
    }
    if (ok && tr.valid && arep.valid && sf.valid) {
        std::cout << "  " << kWarn << " " << sf.patches
                  << " watertight bicubic patch(es) from a layout satisfying Q1-Q5, but the "
                  << "mesh on them is not clean: " << qr.invertedQuads << " inverted element(s), "
                  << qr.nonManifoldEdges << " over-used edge(s), " << qr.cracks
                  << " coincident vertex pair(s).\n";
        return 0;
    }
    if (ok && tr.valid) {
        std::cout << "  " << kWarn
                  << " Psi is a quadrilateral layout in the sense of Definition 2.1, but the "
                  << "arrangement is not yet a set of usable patches: "
                  << (arep.patches - arep.simpleQuads) << " of " << arep.patches
                  << " face(s) are not quadrilaterals with one arc a side"
                  << (sf.foldedPatches > 0 ? ", and some of the rest fold when blended" : "")
                  << ". They were left out of the mesh, which covers the other "
                  << std::fixed << std::setprecision(2) << (100.0 * (1.0 - qr.unmeshedArea))
                  << "% of S in " << std::defaultfloat << qr.quads
                  << " quadrilateral(s). Sec. 4's remedy is Sec. 3.3's -- raise lambda_5, add "
                  << "the Gamma_topo constraint for the pair that nearly met, re-run Stage 6 "
                  << "from the current phi.\n";
        return 0;
    }
    if (ok && tr.valid) {
        std::cout << "  " << kPass
                  << " Psi satisfies Q1-Q5: a quadrilateral layout in the sense of "
                  << "Definition 2.1, and its " << tr.emitted
                  << " separatrices are all finite. Ready for Stage 8 (arrangement).\n";
    } else if (ok) {
        std::cout << "  " << kWarn
                  << " Psi satisfies Q1-Q5 as the energies measure them, but "
                  << (tr.capped + tr.stuck + tr.degenerate)
                  << " separatrix/ces did not terminate: Q5 holds on the Gamma_topo paths "
                  << "E5 was given and not on every integral curve.\n";
    } else {
        std::cout << "  " << kFail
                  << " The continuation did not reach a valid layout; see the residuals "
                  << "above.\n";
    }
    return ok ? 0 : 5;
}
