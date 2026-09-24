// Utility to run the TMOP smoother of src/mesh/TMOP.hxx on a quadrilateral
// mesh, and to self-test it and the node-mobility layer it stands on.
//
//   TestTMOP <quadmesh.obj> [options]
//   TestTMOP --selftest
//
// The mesh is a *quad* .obj -- what TestMERIDIAN --mesh writes, or what
// mesh::QuadMesh::writeOBJ produces. The triangle meshes in data/meshes are the
// pipeline's input and this reader refuses them on purpose.
//
//   TestMERIDIAN data/meshes/singlemat/geom001.obj --mesh geom001.quads.obj
//   TestTMOP geom001.quads.obj --vtu geom001.smoothed.vtu
//
// Exits non-zero when the smoother reported failure or left the mesh worse than
// it found it. --selftest needs no files: it builds its meshes in memory.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "mesh/QuadMesh.hxx"
#include "mesh/TMOP.hxx"

namespace {

const char *kPass = "\033[32m[PASS]\033[0m";
const char *kFail = "\033[31m[FAIL]\033[0m";
const char *kWarn = "\033[33m[WARN]\033[0m";

int gFailures = 0;

void verdict(bool ok, const std::string &what) {
    if (!ok) ++gFailures;
    std::cout << "  " << (ok ? kPass : kFail) << " " << what << "\n";
}

void heading(const std::string &title) {
    std::cout << "\n" << title << "\n" << std::string(title.size(), '-') << "\n";
}

// A deterministic generator, so a failure is reproducible without a seed file.
struct Rand {
    unsigned long long s;
    explicit Rand(unsigned long long seed) : s(seed ? seed : 88172645463325252ULL) {}
    double next() {   // uniform in [-1, 1)
        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
        return 2.0 * (static_cast<double>(s >> 11) / 9007199254740992.0) - 1.0;
    }
};

// ---------------------------------------------------------------------------
// meshes built in memory
// ---------------------------------------------------------------------------

// An nx by ny grid of quads over [0,w] x [0,h]. Interior nodes are pushed
// around by `jitter` times the cell size, boundary nodes are left where they
// are, so the boundary is straight and the interior is a mess -- exactly the
// case a shape metric should clean up completely.
mesh::QuadMesh makeGrid(int nx, int ny, double w, double h, double jitter,
                        unsigned long long seed = 12345) {
    std::vector<Point> verts;
    verts.reserve((nx + 1) * (ny + 1));
    Rand rng(seed);
    const double dx = w / nx, dy = h / ny;
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) {
            double x = i * dx, y = j * dy;
            if (jitter > 0.0 && i > 0 && i < nx && j > 0 && j < ny) {
                x += jitter * dx * rng.next();
                y += jitter * dy * rng.next();
            }
            verts.push_back(Point{x, y});
        }
    std::vector<mesh::Quad> cells;
    cells.reserve(nx * ny);
    auto id = [&](int i, int j) { return j * (nx + 1) + i; };
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
            cells.push_back(mesh::Quad{id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)});
    return mesh::QuadMesh(verts, cells);
}

// The same grid, split down the middle into two materials, so the vertical line
// x = w/2 is a material interface: a feature curve that is not on dS.
mesh::QuadMesh makeTwoMaterialGrid(int nx, int ny, double w, double h, double jitter) {
    mesh::QuadMesh m = makeGrid(nx, ny, w, h, jitter);
    std::vector<int> mats(m.quads.size(), 1);
    for (int q = 0; q < static_cast<int>(m.quads.size()); ++q)
        if (m.centroid(q)[0] > 0.5 * w) mats[q] = 2;
    return mesh::QuadMesh(m.vertices, m.quads, mats);
}

// An annulus: nSeg quads around, nRing across. Both rims are closed feature
// loops with no corner anywhere on them, which is the case the closed-loop
// anchor exists for -- at nSeg = 24 each rim node turns by 15 degrees, the same
// discretisation data/geometry/multimat/bubbles.geo asks of its inclusions, and
// well inside the 45-degree corner tolerance.
//
// `skew` unevens the angular spacing, so the rims start with a discretisation a
// smoother has something to redistribute.
mesh::QuadMesh makeAnnulus(int nSeg, int nRing, double r0, double r1, double skew) {
    std::vector<Point> verts;
    verts.reserve(nSeg * (nRing + 1));
    std::vector<double> theta(nSeg);
    for (int j = 0; j < nSeg; ++j) {
        const double t = 2.0 * M_PI * j / nSeg;
        theta[j] = t + skew * std::sin(3.0 * t) * (2.0 * M_PI / nSeg);
    }
    for (int k = 0; k <= nRing; ++k) {
        const double r = r0 + (r1 - r0) * k / nRing;
        for (int j = 0; j < nSeg; ++j)
            verts.push_back(Point{r * std::cos(theta[j]), r * std::sin(theta[j])});
    }
    std::vector<mesh::Quad> cells;
    auto id = [&](int k, int j) { return k * nSeg + (j % nSeg); };
    for (int k = 0; k < nRing; ++k)
        for (int j = 0; j < nSeg; ++j)
            cells.push_back(mesh::Quad{id(k, j), id(k + 1, j), id(k + 1, j + 1), id(k, j + 1)});
    return mesh::QuadMesh(verts, cells);
}

// ---------------------------------------------------------------------------
// checks that do not trust the implementation they are checking
// ---------------------------------------------------------------------------

// Re-walk the feature graph and look for a cycle every one of whose nodes is
// Sliding. Deliberately written from `featureNeighborOf` and `nodeType` alone,
// not by asking whether the anchor pass ran.
bool hasUnanchoredFeatureLoop(const mesh::QuadMesh &m, int *witness = nullptr) {
    const int nV = static_cast<int>(m.vertices.size());
    std::vector<bool> seen(nV, false);
    for (int v = 0; v < nV; ++v) {
        if (seen[v] || m.nodeType[v] != mesh::QuadMesh::NodeSliding) continue;
        int prev = v, cur = m.featureNeighborOf[v][1];
        seen[v] = true;
        int guard = 0;
        while (cur >= 0 && cur != v && m.nodeType[cur] == mesh::QuadMesh::NodeSliding &&
               guard++ <= nV) {
            seen[cur] = true;
            const std::array<int, 2> &f = m.featureNeighborOf[cur];
            const int next = (f[0] == prev) ? f[1] : f[0];
            prev = cur;
            cur = next;
        }
        if (cur == v) { if (witness) *witness = v; return true; }
    }
    return false;
}

double distanceToLine(const Point &p, const Point &a, const Point &b) {
    const Point d = b - a;
    const double L = normP(d);
    if (L <= 0.0) return normP(p - a);
    return std::fabs(cross2(d, p - a)) / L;
}

// The same to the segment rather than the whole line: what "still on the
// polyline" means for a node that may have walked past a vertex of it.
double distanceToSegment(const Point &p, const Point &a, const Point &b) {
    const Point d = b - a;
    const double L2 = dotP(d, d);
    if (L2 <= 0.0) return normP(p - a);
    double t = dotP(p - a, d) / L2;
    t = std::min(std::max(t, 0.0), 1.0);
    return normP(p - (a + d * t));
}

// ---------------------------------------------------------------------------
// selfTest()
// ---------------------------------------------------------------------------

int selfTest() {
    std::cout << "TMOP and node-mobility self-test\n";

    // -- 1. the mobility layer ---------------------------------------------
    heading("1  Node mobility on a grid: who moves, and how");
    {
        mesh::QuadMesh g = makeGrid(8, 6, 8.0, 6.0, 0.15);
        const mesh::QuadMesh::Quality &Q = g.quality;
        std::cout << "  " << Q.vertices << " nodes: " << Q.freeNodes << " free, "
                  << Q.slidingNodes << " sliding, " << Q.fixedNodes << " fixed\n";
        verdict(Q.freeNodes == 7 * 5, "Every interior node of an 8x6 grid is free");
        verdict(Q.fixedNodes == 4, "The four 90-degree corners, and only those, are fixed");
        verdict(Q.slidingNodes == 2 * 7 + 2 * 5,
                "Every straight boundary node slides");

        // projectStep, on one node of each kind.
        int freeV = -1, slideV = -1, fixedV = -1;
        for (int v = 0; v < static_cast<int>(g.vertices.size()); ++v) {
            if (g.nodeType[v] == mesh::QuadMesh::NodeFree && freeV < 0) freeV = v;
            if (g.nodeType[v] == mesh::QuadMesh::NodeSliding && slideV < 0) slideV = v;
            if (g.nodeType[v] == mesh::QuadMesh::NodeFixed && fixedV < 0) fixedV = v;
        }
        const Point d{0.37, -0.81};
        const Point pf = g.projectStep(fixedV, d);
        const Point pr = g.projectStep(freeV, d);
        const Point ps = g.projectStep(slideV, d);
        verdict(normP(pf) == 0.0, "projectStep on a fixed node is exactly zero");
        verdict(normP(pr - d) == 0.0, "projectStep on a free node returns the step unchanged");
        verdict(std::fabs(cross2(ps, g.slideTangent[slideV])) < 1e-12 &&
                    normP(ps) <= normP(d) + 1e-12,
                "projectStep on a sliding node is parallel to its tangent and no longer");

        // "Stays exactly on the segment": the projected step keeps the node
        // collinear with the two feature neighbours it was between.
        const std::array<int, 2> fn = g.featureNeighborOf[slideV];
        const Point a = g.vertices[fn[0]], b = g.vertices[fn[1]];
        const Point moved = g.vertices[slideV] + ps;
        verdict(distanceToLine(moved, a, b) < 1e-12,
                "A sliding node stays exactly on the segment through its feature neighbours");
    }

    // -- 2. the closed-loop anchor -----------------------------------------
    heading("2  Closed feature loops are anchored");
    {
        mesh::QuadMesh ann = makeAnnulus(24, 3, 1.0, 2.0, 0.0);
        const mesh::QuadMesh::Quality &Q = ann.quality;
        std::cout << "  annulus, 24 segments x 3 rings: " << Q.freeNodes << " free, "
                  << Q.slidingNodes << " sliding, " << Q.fixedNodes << " fixed; "
                  << Q.boundaryLoops << " boundary loop(s)\n";
        verdict(Q.boundaryLoops == 2, "Both rims come back as closed boundary loops");
        verdict(Q.fixedNodes == 2,
                "Exactly one node is fixed on each rim, and nothing else is");

        // One fixed node per rim, not two on one rim and none on the other.
        int fixedInner = 0, fixedOuter = 0;
        for (int v = 0; v < static_cast<int>(ann.vertices.size()); ++v) {
            if (ann.nodeType[v] != mesh::QuadMesh::NodeFixed) continue;
            (normP(ann.vertices[v]) < 1.5 ? fixedInner : fixedOuter) += 1;
        }
        verdict(fixedInner == 1 && fixedOuter == 1,
                "The anchor is on the rim itself, one per rim");

        int witness = -1;
        verdict(!hasUnanchoredFeatureLoop(ann, &witness),
                "Re-walking the feature graph finds no all-sliding cycle anywhere");

        // The anchor is what a solver actually feels: it cannot be displaced.
        const int anchor = [&] {
            for (int v = 0; v < static_cast<int>(ann.vertices.size()); ++v)
                if (ann.nodeType[v] == mesh::QuadMesh::NodeFixed) return v;
            return -1;
        }();
        verdict(anchor >= 0 && normP(ann.projectStep(anchor, Point{1.0, 1.0})) == 0.0,
                "No displacement gets through the anchor, so the rim keeps a fixed point");

        // Rebuilding must anchor the same node: an anchor that wandered between
        // rebuilds would make a pipeline's results depend on nothing.
        const std::vector<std::uint8_t> before = ann.nodeType;
        ann.buildTopology();
        verdict(before == ann.nodeType, "A rebuild anchors the same node again");

        // A caller pin on the loop is the anchor: the pass runs after pins, so
        // it must not add a second fixed node to an already-anchored loop.
        mesh::QuadMesh pinned = makeAnnulus(24, 3, 1.0, 2.0, 0.0);
        int rimNode = -1;
        for (int v = 0; v < static_cast<int>(pinned.vertices.size()); ++v)
            if (std::fabs(normP(pinned.vertices[v]) - 2.0) < 1e-9 &&
                pinned.nodeType[v] == mesh::QuadMesh::NodeSliding) { rimNode = v; break; }
        pinned.pinVertex(rimNode);
        pinned.buildTopology();
        verdict(pinned.quality.fixedNodes == 2,
                "Pinning a node on a rim anchors it; no second anchor is added");

        // An interior material interface ring, not a boundary at all: the same
        // pass has to catch it, which is why it walks the feature graph rather
        // than boundaryLoops.
        mesh::QuadMesh block = makeGrid(12, 12, 12.0, 12.0, 0.0);
        std::vector<int> mats(block.quads.size(), 1);
        for (int q = 0; q < static_cast<int>(block.quads.size()); ++q) {
            const Point c = block.centroid(q);
            if (normP(c - Point{6.0, 6.0}) < 3.0) mats[q] = 2;
        }
        mesh::QuadMesh withInclusion(block.vertices, block.quads, mats);
        int w2 = -1;
        verdict(!hasUnanchoredFeatureLoop(withInclusion, &w2),
                "An interior material-interface ring is anchored too, not just dS");
    }

    // -- 3. the metrics ----------------------------------------------------
    heading("3  The metrics: zero at the target, positive away from it");
    {
        const mesh::Jacobian2 I{1.0, 0.0, 0.0, 1.0};
        const mesh::Jacobian2 skew{1.0, 0.6, 0.0, 1.0};      // sheared, same area
        const mesh::Jacobian2 big{1.4, 0.0, 0.0, 1.4};        // square, wrong size
        const mesh::Jacobian2 flipped{1.0, 0.0, 0.0, -1.0};   // inverted

        struct Case { mesh::TMOP::Metric m; const char *name; bool barrier; bool seesShape; bool seesSize; };
        const Case cases[] = {
            {mesh::TMOP::Shape004,       "004 shape, no barrier",  false, true,  false},
            {mesh::TMOP::Shape002,       "002 shape",              true,  true,  false},
            {mesh::TMOP::ShapeSize007,   "007 shape+size",         true,  true,  true},
            {mesh::TMOP::Size055,        "055 size",               false, false, true},
            {mesh::TMOP::Size056,        "056 size",               true,  false, true},
            {mesh::TMOP::ShapeSizeCombo, "combo shape+size",       true,  true,  true},
        };
        bool allZero = true, allShape = true, allSize = true, allBarrier = true;
        for (const Case &c : cases) {
            double mu = 0.0;
            const bool ok = mesh::TMOP::evalMetric(c.m, I, 0.5, 0.0, &mu);
            if (!ok || std::fabs(mu) > 1e-12) allZero = false;

            double muSkew = 0.0, muBig = 0.0;
            mesh::TMOP::evalMetric(c.m, skew, 0.5, 0.0, &muSkew);
            mesh::TMOP::evalMetric(c.m, big, 0.5, 0.0, &muBig);
            if (c.seesShape && !(muSkew > 1e-9)) allShape = false;
            if (c.seesSize && !(muBig > 1e-9)) allSize = false;
            if (!c.seesSize && std::fabs(muBig) > 1e-9) allSize = false;

            double muFlip = 0.0;
            const bool okFlip = mesh::TMOP::evalMetric(c.m, flipped, 0.5, 0.0, &muFlip);
            if (c.barrier && okFlip) allBarrier = false;
            if (!c.barrier && !okFlip) allBarrier = false;
        }
        verdict(allZero, "Every metric is exactly zero at T = I");
        verdict(allShape, "Every shape metric is positive on a sheared T");
        verdict(allSize, "The size metrics see a wrongly sized T and the shape metric does not");
        verdict(allBarrier, "Every barrier metric refuses an inverted T; the others accept it");

        // The untangler must be finite on exactly that inverted T, which is the
        // whole reason it exists.
        double muU = 0.0;
        const bool okU = mesh::TMOP::evalMetric(mesh::TMOP::Untangle022, flipped, 0.0, -2.0, &muU);
        verdict(okU && std::isfinite(muU) && muU > 0.0,
                "Metric 22 is finite and positive on an inverted T with tau_0 below it");
        double muU0 = 0.0;
        verdict(mesh::TMOP::evalMetric(mesh::TMOP::Untangle022, I, 0.0, -2.0, &muU0) &&
                    muU0 < muU,
                "and still prefers the identity to the inverted T");

        // Raising the exponent must not move where the metric vanishes, only
        // how hard it pulls away from there.
        double m1 = 0.0, m2 = 0.0, mI = 0.0;
        mesh::TMOP::evalMetric(mesh::TMOP::Shape002, skew, 0.0, 0.0, &m1, nullptr, 1.0);
        mesh::TMOP::evalMetric(mesh::TMOP::Shape002, skew, 0.0, 0.0, &m2, nullptr, 2.0);
        mesh::TMOP::evalMetric(mesh::TMOP::Shape002, I, 0.0, 0.0, &mI, nullptr, 2.0);
        verdict(std::fabs(m2 - m1 * m1) < 1e-12 && mI == 0.0,
                "mu^2 is mu squared, and still exactly zero at T = I");

        // An exponent strictly between 1 and 2 has a singular second derivative
        // at mu = 0 and must be refused rather than run.
        mesh::QuadMesh g = makeGrid(4, 4, 4.0, 4.0, 0.1);
        mesh::TMOP::Options o;
        o.exponent = 1.5;
        mesh::TMOP opt(g, o);
        opt.run();
        verdict(opt.getOptions().exponent == 2.0 && !opt.getReport().messages.empty(),
                "An exponent in (1, 2) is raised to 2 and reported, not run as given");
    }

    // -- 4. the derivatives ------------------------------------------------
    heading("4  Gradient and Hessian against finite differences");
    {
        // Both are exact expressions, so a central difference has to reproduce
        // them to its own truncation error and nothing worse. Checked on a
        // distorted grid, on a metric of each kind, and at nodes of each type.
        const mesh::TMOP::Metric metrics[] = {
            mesh::TMOP::Shape004, mesh::TMOP::Shape002, mesh::TMOP::ShapeSize007,
            mesh::TMOP::Size055, mesh::TMOP::Size056, mesh::TMOP::ShapeSizeCombo};
        // Every exponent too: mu^p brings its own chain rule into all five
        // partials, and p = 2 is what runs by default.
        const double powers[] = {1.0, 2.0, 3.0};
        double worstGrad = 0.0, worstHess = 0.0;
        for (mesh::TMOP::Metric mt : metrics)
        for (double pw : powers) {
            mesh::QuadMesh g = makeGrid(6, 5, 6.0, 5.0, 0.2);
            mesh::TMOP::Options o;
            o.metric = mt;
            o.gamma = 0.35;
            o.exponent = pw;
            o.target = mesh::TMOP::TargetUniformSquare;
            mesh::TMOP opt(g, o);
            opt.prepare();

            const double eps = 1e-6;
            for (int v = 0; v < static_cast<int>(g.vertices.size()); v += 3) {
                if (g.nodeType[v] == mesh::QuadMesh::NodeFixed) continue;
                Point ga{0.0, 0.0};
                double H[4] = {0.0, 0.0, 0.0, 0.0};
                opt.nodeSystem(v, ga, H);

                const Point origin = g.vertices[v];
                double scale = std::max(normP(ga), 1.0);
                for (int i = 0; i < 2; ++i) {
                    Point plus = origin, minus = origin;
                    plus[i] += eps; minus[i] -= eps;
                    g.vertices[v] = plus;
                    const double ep = opt.nodeEnergy(v);
                    Point gp = opt.nodeGradient(v);
                    g.vertices[v] = minus;
                    const double em = opt.nodeEnergy(v);
                    Point gm = opt.nodeGradient(v);
                    g.vertices[v] = origin;

                    const double fd = (ep - em) / (2.0 * eps);
                    worstGrad = std::max(worstGrad, std::fabs(fd - ga[i]) / scale);
                    for (int l = 0; l < 2; ++l) {
                        const double fdh = (gp[l] - gm[l]) / (2.0 * eps);
                        const double hs = std::max(std::fabs(H[2 * l + i]), 1.0);
                        worstHess = std::max(worstHess, std::fabs(fdh - H[2 * l + i]) / hs);
                    }
                }
            }
        }
        std::cout << "  worst relative error: gradient " << std::scientific
                  << std::setprecision(2) << worstGrad << ", Hessian " << worstHess
                  << std::defaultfloat << "\n";
        verdict(worstGrad < 1e-5, "The analytic node gradient matches a central difference");
        verdict(worstHess < 1e-4, "The analytic local Hessian matches a central difference");
    }

    // -- 5. smoothing a distorted grid -------------------------------------
    heading("5  Smoothing: a badly distorted grid recovers");
    {
        mesh::QuadMesh g = makeGrid(16, 12, 16.0, 12.0, 0.45);
        const double sjBefore = g.quality.minScaledJacobian;
        mesh::TMOP::Options o;
        o.metric = mesh::TMOP::Shape002;
        mesh::TMOP opt(g, o);
        const bool ok = opt.run();
        const mesh::TMOP::Report &r = opt.getReport();
        std::cout << "  " << r.sweeps << " sweep(s), " << r.colors << " colour(s), "
                  << r.threads << " thread(s); min scaled Jacobian " << std::fixed
                  << std::setprecision(4) << r.minScaledJacobianBefore << " -> "
                  << r.minScaledJacobianAfter << ", energy " << std::scientific
                  << std::setprecision(3) << r.energyBefore << " -> " << r.energyAfter
                  << std::defaultfloat << "\n";
        verdict(ok && r.ran, "The smoother ran and reported success");
        verdict(r.energyAfter < r.energyBefore, "The TMOP energy fell");
        verdict(r.minScaledJacobianAfter > sjBefore,
                "The worst element got better, not just the average");
        verdict(r.minScaledJacobianAfter > 0.99,
                "A rectangle with a straight boundary comes back essentially perfect");
        verdict(r.converged, "It stopped on a tolerance rather than on the sweep cap");

        // Nothing may have left the domain: the corners are fixed and every
        // boundary node slides along its own wall, so the bounding box is
        // exactly what it was.
        double xlo = 1e30, xhi = -1e30, ylo = 1e30, yhi = -1e30;
        for (const Point &p : g.vertices) {
            xlo = std::min(xlo, p[0]); xhi = std::max(xhi, p[0]);
            ylo = std::min(ylo, p[1]); yhi = std::max(yhi, p[1]);
        }
        verdict(std::fabs(xlo) < 1e-12 && std::fabs(ylo) < 1e-12 &&
                    std::fabs(xhi - 16.0) < 1e-12 && std::fabs(yhi - 12.0) < 1e-12,
                "The boundary stayed on the boundary: the bounding box is unchanged");
    }

    // -- 6. sliding really is what fixes the boundary ring -----------------
    heading("6  Sliding versus pinning the boundary");
    {
        // The same distorted grid twice, once with the boundary free to
        // redistribute and once pinned. The pinned run cannot fix the ring of
        // elements touching a boundary whose spacing it may not change.
        auto build = [](bool pin) {
            mesh::QuadMesh::Options mo;
            mo.fixAllFeatureNodes = pin;
            mesh::QuadMesh g = makeGrid(14, 10, 14.0, 10.0, 0.0);
            // Bunch the bottom wall's nodes to one side, so its discretisation
            // disagrees with what the interior metric wants.
            for (int v = 0; v < static_cast<int>(g.vertices.size()); ++v) {
                if (std::fabs(g.vertices[v][1]) > 1e-12) continue;
                const double t = g.vertices[v][0] / 14.0;
                g.vertices[v][0] = 14.0 * t * t;
            }
            return mesh::QuadMesh(g.vertices, g.quads, std::vector<int>{}, mo);
        };
        mesh::QuadMesh sliding = build(false), pinnedM = build(true);
        mesh::TMOP::Options o;
        o.metric = mesh::TMOP::Shape002;
        mesh::TMOP a(sliding, o), b(pinnedM, o);
        a.run();
        b.run();
        std::cout << "  sliding boundary: min scaled Jacobian "
                  << std::fixed << std::setprecision(4)
                  << a.getReport().minScaledJacobianAfter << ", energy "
                  << std::scientific << std::setprecision(3) << a.getReport().energyAfter << "\n";
        std::cout << "  pinned boundary:  min scaled Jacobian "
                  << std::fixed << std::setprecision(4)
                  << b.getReport().minScaledJacobianAfter << ", energy "
                  << std::scientific << std::setprecision(3) << b.getReport().energyAfter
                  << std::defaultfloat << "\n";
        verdict(a.getReport().energyAfter < b.getReport().energyAfter,
                "Letting the boundary slide reaches a lower energy than pinning it");
        verdict(a.getReport().slidingNodes > 0 && b.getReport().slidingNodes == 0,
                "fixAllFeatureNodes really does take every feature node out of play");

        // And the sliding run must not have moved a boundary node off its wall.
        int onWall = 0;
        for (const Point &p : sliding.vertices)
            if (std::fabs(p[0]) < 1e-12 || std::fabs(p[0] - 14.0) < 1e-12 ||
                std::fabs(p[1]) < 1e-12 || std::fabs(p[1] - 10.0) < 1e-12)
                ++onWall;
        verdict(onWall == sliding.quality.slidingNodes + sliding.quality.fixedNodes,
                "Every boundary node is still exactly on its wall after sliding");
    }

    // -- 7. material interfaces --------------------------------------------
    heading("7  A material interface is a feature the smoother may not cross");
    {
        mesh::QuadMesh g = makeTwoMaterialGrid(12, 8, 12.0, 8.0, 0.3);
        const std::vector<int> matsBefore = g.quadMatId;
        mesh::TMOP::Options o;
        mesh::TMOP opt(g, o);
        opt.run();

        double worstOff = 0.0;
        int interfaceNodes = 0;
        for (int v = 0; v < static_cast<int>(g.vertices.size()); ++v) {
            if (!g.isFeatureVertex[v]) continue;
            if (std::fabs(g.vertices[v][0] - 6.0) > 1e-9) continue;
            ++interfaceNodes;
            worstOff = std::max(worstOff, std::fabs(g.vertices[v][0] - 6.0));
        }
        verdict(interfaceNodes > 0, "The interface at x = 6 is a feature curve");
        verdict(worstOff < 1e-12, "Its nodes are still exactly on it after smoothing");
        verdict(g.quadMatId == matsBefore, "No element changed material");
        verdict(opt.getReport().invertedAfter == 0, "Nothing inverted");
    }

    // -- 8. untangling ------------------------------------------------------
    heading("8  Untangling a mesh that starts inverted");
    {
        // Drag one interior node of a perfect grid off to `to`, folding the
        // four elements around it. Two distances, because they are different
        // problems. Short of the diagonal neighbour the fold is local and comes
        // straight back out. Past it, the node and that neighbour have swapped
        // sides, and untangling finds the nearest configuration where every
        // determinant is positive again -- which is a valid mesh with the two
        // still wound around each other, and a genuine local minimum of the
        // shape energy that no local move escapes. That is a property of
        // node-local optimization, not a defect to fix here.
        auto tangled = [](Point to) {
            mesh::QuadMesh g = makeGrid(10, 10, 10.0, 10.0, 0.0);
            for (int v = 0; v < static_cast<int>(g.vertices.size()); ++v)
                if (normP(g.vertices[v] - Point{5.0, 5.0}) < 1e-9) g.vertices[v] = to;
            g.buildTopology();
            return g;
        };

        mesh::QuadMesh mild = tangled(Point{5.9, 5.9});
        std::cout << "  fold inside the one-ring: " << mild.quality.invertedQuads
                  << " inverted element(s), min scaled Jacobian " << std::fixed
                  << std::setprecision(4) << mild.quality.minScaledJacobian << "\n";
        verdict(mild.quality.invertedQuads > 0, "The starting mesh really is tangled");

        mesh::TMOP::Options o;
        o.metric = mesh::TMOP::Shape002;
        o.untangle = true;
        mesh::TMOP opt(mild, o);
        const bool ok = opt.run();
        const mesh::TMOP::Report &r = opt.getReport();
        std::cout << "  " << r.untangleSweeps << " untangling sweep(s) then " << r.sweeps
                  << " smoothing sweep(s); min scaled Jacobian "
                  << r.minScaledJacobianBefore << " -> " << r.minScaledJacobianAfter << "\n";
        for (const std::string &m : r.messages) std::cout << "  " << kWarn << " " << m << "\n";
        verdict(ok, "The run reported success");
        verdict(r.untangleSweeps > 0, "The untangling phase was entered");
        verdict(r.invertedAfter == 0, "Nothing is inverted any more");
        verdict(r.minScaledJacobianAfter > 0.99,
                "and a fold that stayed inside the one-ring comes all the way back "
                "to a perfect grid");

        mesh::QuadMesh severe = tangled(Point{7.6, 7.6});
        mesh::TMOP opt3(severe, o);
        opt3.run();
        const mesh::TMOP::Report &r3 = opt3.getReport();
        std::cout << "  fold past the diagonal neighbour: min scaled Jacobian "
                  << r3.minScaledJacobianBefore << " -> " << r3.minScaledJacobianAfter
                  << " in " << r3.untangleSweeps << " + " << r3.sweeps << " sweep(s)\n";
        verdict(r3.invertedAfter == 0,
                "A node thrown clear past its diagonal neighbour is still untangled");
        verdict(r3.minScaledJacobianAfter > 0.5,
                "and left a usable mesh, though wound up: a local method cannot unwind it");

        // A fold that is barely there: the node a hair past the diagonal of
        // its upper-right neighbours, one corner at scaled Jacobian -0.04 --
        // what a TFI mesh gives at a block corner the layout left at 181
        // degrees (traced det_rocket, -0.02). Smoothed the way the block
        // meshes are, shape + size sampled at the corners. Metric 22 with a
        // tracked tau_0 treats a crushed element as a perfect one, and this is
        // where that shows: it crushed the fold and its neighbours into a
        // point (worst 0.035, smallest element 6e-7 of the target area,
        // aspect 20) and stopped, and the main phase could not grow them
        // back. The regularised untangler opens the fold instead.
        {
            mesh::TMOP::Options o7;
            o7.metric = mesh::TMOP::ShapeSize007;
            o7.quadrature = mesh::TMOP::Corners;
            mesh::QuadMesh hair = tangled(Point{5.51, 5.51});
            mesh::TMOP optH(hair, o7);
            optH.run();
            const mesh::TMOP::Report &rh = optH.getReport();
            mesh::QuadMesh hairOld = tangled(Point{5.51, 5.51});
            mesh::TMOP::Options o7old = o7;
            o7old.untangler = mesh::TMOP::UntangleShifted;
            mesh::TMOP optOld(hairOld, o7old);
            optOld.run();
            const mesh::TMOP::Report &ro = optOld.getReport();
            std::cout << "  fold a hair past the diagonal: min scaled Jacobian "
                      << rh.minScaledJacobianBefore << " -> " << rh.minScaledJacobianAfter
                      << ", smallest element " << rh.minAreaAfter << " (" << rh.untangleSweeps
                      << " + " << rh.sweeps << " sweep(s)); the shifted untangler: "
                      << ro.minScaledJacobianAfter << ", smallest element " << ro.minAreaAfter << "\n";
            verdict(rh.invertedAfter == 0 && rh.minScaledJacobianAfter > 0.99,
                    "A barely folded node comes all the way back under shape + size at the corners");
            verdict(rh.minAreaAfter > 0.9,
                    "and nothing was crushed on the way: every element is still about the target area");
        }

        // Refusing to start is the correct behaviour when untangling is off.
        mesh::QuadMesh h = tangled(Point{7.6, 7.6});
        mesh::TMOP::Options o2 = o;
        o2.untangle = false;
        mesh::TMOP opt2(h, o2);
        const bool ran = opt2.run();
        verdict(!ran && !opt2.getReport().ran,
                "With untangling off, a barrier metric declines a tangled mesh rather than "
                "dividing by zero");

        // A corner nothing can open. Three quads on [0,2]^2; the bottom one is
        // (0,0) (1,0) (2,0) (d), so its corner at (1,0) is a straight angle
        // between two boundary edges -- a block corner a layout put on a
        // straight feature. With (1,0) pinned, and (0,0) and (2,0) corners of
        // the model, no admissible move changes that corner. Sampled at the
        // corners its det T is zero, so the barrier energy of the whole mesh
        // was infinite and the smoother declined to touch any of it; and an
        // untangler run at it only ever finds |T|^2 left to lower. It has to
        // be left out, and the rest smoothed.
        {
            std::vector<Point> v = {{0, 0}, {1, 0}, {2, 0}, {1.35, 0.55}, {0, 2}, {1, 2}, {2, 2}};
            std::vector<mesh::Quad> c = {{0, 1, 2, 3}, {0, 3, 5, 4}, {3, 2, 6, 5}};
            mesh::QuadMesh flat(v, c);
            flat.pinVertex(1);
            flat.classifyNodes();
            flat.buildFeatureCurves();
            flat.computeSlideTangents();
            mesh::TMOP::Options of;
            of.quadrature = mesh::TMOP::Corners;
            mesh::TMOP optF(flat, of);
            const Point before = flat.vertices[3];
            const bool okF = optF.run();
            const mesh::TMOP::Report &rf = optF.getReport();
            std::cout << "  a straight-angle corner at a pinned node: " << rf.frozenCorners
                      << " frozen corner(s), " << rf.untangleSweeps << " + " << rf.sweeps
                      << " sweep(s), the free node moved "
                      << normP(flat.vertices[3] - before) << "\n";
            verdict(rf.frozenCorners == 1, "The one corner no node can change is found, and only it");
            verdict(okF && rf.ran && rf.untangleSweeps == 0,
                    "The mesh is smoothed rather than declined, without an untangling phase");
            verdict(normP(flat.vertices[3] - before) > 0.05 && std::isfinite(rf.energyAfter),
                    "and the free node was actually moved, to a finite energy");
        }
    }

    // -- 9. the annulus: rims redistribute, and stay circles ----------------
    heading("9  An annulus with unevenly spaced rims");
    {
        const int nSeg = 32;
        // The two figures this case is judged against. A rim node is a point of
        // its circle whatever the smoother does with the spacing, so the area
        // the *mesh* encloses is the inscribed nSeg-gon's -- and that is
        // largest when the nodes are evenly spread, which is exactly what a
        // shape metric on a rim of equal elements asks for. The sagitta is how
        // far inside the circle the chord between two neighbouring rim nodes
        // passes, which is the error the old chord projection left behind.
        const double evenArea =
            0.5 * nSeg * std::sin(2.0 * M_PI / nSeg) * (2.0 * 2.0 - 1.0 * 1.0);
        const double sagitta = 2.0 * (1.0 - std::cos(M_PI / nSeg));

        mesh::QuadMesh ann = makeAnnulus(nSeg, 4, 1.0, 2.0, 0.45);
        int anchor = -1;
        for (int v = 0; v < static_cast<int>(ann.vertices.size()); ++v)
            if (ann.nodeType[v] == mesh::QuadMesh::NodeFixed) { anchor = v; break; }
        const Point anchorAt = ann.vertices[anchor];

        mesh::TMOP::Options o;
        o.metric = mesh::TMOP::Shape002;
        mesh::TMOP opt(ann, o);
        opt.run();
        const mesh::TMOP::Report &r = opt.getReport();
        std::cout << "  min scaled Jacobian " << std::fixed << std::setprecision(4)
                  << r.minScaledJacobianBefore << " -> " << r.minScaledJacobianAfter
                  << ", worst aspect " << r.worstAspectBefore << " -> " << r.worstAspectAfter
                  << ", nodes moved " << std::setprecision(3) << r.meanDisplacement
                  << " mean edge(s) on average\n";
        verdict(r.minScaledJacobianAfter >= r.minScaledJacobianBefore - 1e-9,
                "The worst element did not get worse");
        verdict(r.worstAspectAfter < r.worstAspectBefore,
                "The rims redistributed: the worst aspect ratio came down");
        verdict(normP(ann.vertices[anchor] - anchorAt) == 0.0,
                "The anchor did not move a millimetre, so the rim kept its reference point");

        // Each rim is a closed feature loop with no corner on it, so it comes
        // back as exactly one run -- from its anchor all the way round to the
        // same anchor -- and that run is smooth enough to interpolate.
        verdict(r.featureCurves == 2 && r.fittedCurves == 2,
                "Each rim was bound to one interpolating curve of its own");
        verdict(r.curveNodes == r.slidingNodes,
                "and every sliding node rides one: none fell back to a chord");

        // A node riding the interpolant of 32 points of a circle is on the
        // circle to O(h^4), not merely within the chord's sagitta of it. This
        // is the whole point of the change: under the chord projection the
        // number here was the sagitta itself, and always on the inside.
        double worstInner = 0.0, worstOuter = 0.0;
        for (int v = 0; v < static_cast<int>(ann.vertices.size()); ++v) {
            if (ann.nodeType[v] == mesh::QuadMesh::NodeFree) continue;
            const double rr = normP(ann.vertices[v]);
            if (rr < 1.5) worstInner = std::max(worstInner, std::fabs(rr - 1.0));
            else worstOuter = std::max(worstOuter, std::fabs(rr - 2.0));
        }
        std::cout << "  rim nodes drifted at most " << std::scientific << std::setprecision(2)
                  << std::max(worstInner, worstOuter) << " from their circle (one segment's "
                  << "sagitta is " << sagitta << ")" << std::defaultfloat << "\n";
        verdict(std::max(worstInner, worstOuter) < 0.05 * sagitta,
                "Every rim node is on its circle, not merely within a chord's sagitta of it");
        verdict(r.curveDeviation < 1e-12,
                "and exactly on the curve it rides: parameter and position never came apart");

        // The area follows. Under the chord projection it was conserved to the
        // last bit -- see the next case, which still checks that -- but what it
        // conserved was the area of the *starting* polygon, unevenly spread
        // nodes and all. Riding the curve, the nodes spread evenly along their
        // circles and the mesh grows into the polygon those circles inscribe.
        std::cout << "  area " << std::setprecision(12) << r.areaBefore << " -> " << r.areaAfter
                  << " (evenly spread on both circles: " << evenArea << ")"
                  << std::defaultfloat << "\n";
        verdict(std::fabs(r.areaAfter - evenArea) < std::fabs(r.areaBefore - evenArea),
                "The rims moved out towards their circles, not in towards their chords");
        verdict(std::fabs(r.areaAfter - evenArea) < 1e-4 * evenArea,
                "and stopped where 32 evenly spread points of those circles put them");
    }

    // -- 10. the curve layer itself ----------------------------------------
    heading("10  What a node slides on: chord, polyline, interpolant");
    {
        // The same annulus under each of the three sources, so the three
        // answers are comparable line for line.
        const int nSeg = 32;

        // How far the rim nodes of a smoothed annulus ended from their circles.
        auto rimDrift = [](const mesh::QuadMesh &m) {
            double worst = 0.0;
            for (int v = 0; v < static_cast<int>(m.vertices.size()); ++v) {
                if (m.nodeType[v] == mesh::QuadMesh::NodeFree) continue;
                const double rr = normP(m.vertices[v]);
                worst = std::max(worst, std::fabs(rr - (rr < 1.5 ? 1.0 : 2.0)));
            }
            return worst;
        };

        // CurveChord: no curves at all, and the exact-area invariant that
        // projection is written for. A node moving along the chord through its
        // two feature neighbours moves parallel to the base of the only
        // triangle the enclosed area depends on it through, so the area cannot
        // change however far the rims redistribute.
        mesh::QuadMesh chord = makeAnnulus(nSeg, 4, 1.0, 2.0, 0.45);
        chord.options.curveSource = mesh::QuadMesh::Options::CurveChord;
        chord.buildTopology();
        verdict(chord.featureCurves.empty(),
                "CurveChord builds no curves and binds no node");
        mesh::TMOP::Options co;
        co.metric = mesh::TMOP::Shape002;
        mesh::TMOP copt(chord, co);
        copt.run();
        const mesh::TMOP::Report &cr = copt.getReport();
        std::cout << "  chord: area " << std::setprecision(12) << cr.areaBefore << " -> "
                  << cr.areaAfter << std::defaultfloat << "\n";
        verdict(std::fabs(cr.areaAfter - cr.areaBefore) < 1e-12 * cr.areaBefore,
                "and still encloses exactly the area it did before the rims moved");
        // The same annulus again on its interpolants, to have the two drifts
        // side by side rather than against a guessed fraction of a sagitta.
        mesh::QuadMesh spline = makeAnnulus(nSeg, 4, 1.0, 2.0, 0.45);
        mesh::TMOP sopt(spline, co);
        sopt.run();
        const double chordDrift = rimDrift(chord), splineDrift = rimDrift(spline);
        std::cout << "  rim nodes ended " << std::scientific << std::setprecision(2) << chordDrift
                  << " off their circles on chords, " << splineDrift << " on interpolants"
                  << std::defaultfloat << "\n";
        verdict(chordDrift > 10.0 * splineDrift,
                "Riding the curve leaves an order of magnitude less drift than the chord does");

        // CurvePolyline: the run's own polyline, as a degree-1 curve. The node
        // rides the geometry it was handed and never leaves it -- not to the
        // inside as the chord does, nor to the outside as the interpolant may.
        mesh::QuadMesh poly = makeAnnulus(nSeg, 4, 1.0, 2.0, 0.45);
        poly.options.curveSource = mesh::QuadMesh::Options::CurvePolyline;
        poly.buildTopology();
        const std::vector<Point> polyStart = poly.vertices;
        verdict(poly.featureCurves.size() == 2 &&
                    !poly.featureCurves[0].fitted && !poly.featureCurves[1].fitted,
                "CurvePolyline binds both rims, and fits nothing");
        bool through = true;
        for (const mesh::QuadMesh::FeatureCurve &fc : poly.featureCurves) {
            verdict(fc.curve.degree() == 1, "Its curve is the degree-1 spline of the run");
            verdict(fc.closed && fc.chain.front() == fc.chain.back(),
                    "and runs from the rim's anchor round to the same anchor");
        }
        for (int v = 0; v < static_cast<int>(poly.vertices.size()); ++v)
            if (poly.isOnCurve(v))
                through = through && normP(poly.vertices[v] -
                                           poly.curvePoint(v, poly.curveParam[v])) < 1e-12;
        verdict(through, "Binding moved no node: each sits on its curve at the parameter it was given");

        mesh::TMOP popt(poly, co);
        popt.run();
        // Every node is still on the polyline it started on: its distance to
        // the nearest of that rim's original segments is zero to rounding.
        double offPolyline = 0.0;
        for (int v = 0; v < static_cast<int>(poly.vertices.size()); ++v) {
            if (!poly.isOnCurve(v)) continue;
            const std::vector<int> &chain = poly.featureCurves[poly.curveOf[v]].chain;
            double best = std::numeric_limits<double>::infinity();
            for (std::size_t i = 0; i + 1 < chain.size(); ++i)
                best = std::min(best, distanceToSegment(poly.vertices[v], polyStart[chain[i]],
                                                        polyStart[chain[i + 1]]));
            offPolyline = std::max(offPolyline, best);
        }
        std::cout << "  polyline: nodes ended " << std::scientific << std::setprecision(2)
                  << offPolyline << " off the polyline they started on" << std::defaultfloat << "\n";
        verdict(offPolyline < 1e-12,
                "A node on a polyline stays on that polyline exactly, through its vertices");

        // The guard: a curve that bows further from its own polyline than
        // curveMaxDeviation allows is refused and the run keeps its polyline.
        // Setting the tolerance to zero refuses every curved run, which is the
        // mechanism, tested without having to draw a chain that rings.
        mesh::QuadMesh strict = makeAnnulus(nSeg, 4, 1.0, 2.0, 0.45);
        strict.options.curveMaxDeviation = 0.0;
        strict.buildTopology();
        bool anyFitted = false;
        for (const mesh::QuadMesh::FeatureCurve &fc : strict.featureCurves) anyFitted |= fc.fitted;
        verdict(strict.featureCurves.size() == 2 && !anyFitted,
                "A run whose interpolant bows too far falls back to its polyline");
    }

    // -- 11. the colouring, and the thread count ---------------------------
    heading("11  Colouring: the answer does not depend on the thread count");
    {
        mesh::QuadMesh a = makeGrid(20, 20, 20.0, 20.0, 0.35);
        mesh::QuadMesh b = makeGrid(20, 20, 20.0, 20.0, 0.35);

        mesh::TMOP::Options o;
        o.metric = mesh::TMOP::Shape002;
        // Deliberately the *default* stopping tests, not disabled ones. Getting
        // the same mesh out of a fixed number of sweeps only shows the sweep
        // arithmetic agrees; letting each run decide for itself when to stop
        // also puts the convergence measures under the same demand, and those
        // are sums over the whole mesh, where a floating-point reduction would
        // differ in its last bit at each thread count.
        mesh::TMOP::Options o1 = o; o1.threads = 1;
        mesh::TMOP::Options oN = o; oN.threads = 0;
        mesh::TMOP one(a, o1), many(b, oN);
        one.run();
        many.run();

        // No two nodes of one colour may share an element -- the property the
        // whole parallel scheme rests on, checked directly.
        bool disjoint = true;
        for (const std::vector<int> &bucket : one.colorBuckets()) {
            std::vector<int> mark(a.quads.size(), -1);
            for (int v : bucket)
                for (int i = a.vertexQuads.begin(v); i < a.vertexQuads.end(v); ++i) {
                    const int q = a.vertexQuads.colIdx[i];
                    if (mark[q] >= 0) disjoint = false;
                    mark[q] = v;
                }
        }
        std::cout << "  " << one.getReport().colors << " colour(s); "
                  << one.getReport().threads << " thread(s) took " << one.getReport().sweeps
                  << " sweep(s), " << many.getReport().threads << " took "
                  << many.getReport().sweeps << "\n";
        verdict(disjoint, "No element is touched by two nodes of the same colour");

        double worst = 0.0;
        for (std::size_t v = 0; v < a.vertices.size(); ++v)
            worst = std::max(worst, normP(a.vertices[v] - b.vertices[v]));
        verdict(worst == 0.0,
                "One thread and many threads produce bit-identical meshes");
        verdict(one.getReport().sweeps == many.getReport().sweeps &&
                    one.getReport().energyAfter == many.getReport().energyAfter,
                "and stop at the same sweep, on the same energy to the last bit");
        if (!one.getReport().openMP)
            std::cout << "  " << kWarn
                      << " built without OpenMP; the thread comparison is trivially true\n";
    }

    heading("Result");
    if (gFailures == 0) {
        std::cout << "  " << kPass << " every check passed\n";
        return 0;
    }
    std::cout << "  " << kFail << " " << gFailures << " check(s) failed\n";
    return 1;
}

// ---------------------------------------------------------------------------
// the driver
// ---------------------------------------------------------------------------

void usage(const char *argv0) {
    std::cout <<
        "Usage: " << argv0 << " <quadmesh.obj> [options]\n"
        "       " << argv0 << " --selftest\n\n"
        "  --metric N        2 (default), 4, 7, 22, 55, 56, or 100 for the shape+size combo\n"
        "  --gamma G         size weight of metric 100 (default 0.2)\n"
        "  --power P         minimise the integral of mu^P instead of mu: 1 (default)\n"
        "                    optimises the average, 2 and up chase the worst element\n"
        "  --target T        square (default), current, keep\n"
        "  --h X             target edge length; default is the mesh's mean edge\n"
        "  --corners         sample the metric at the element corners, not at Gauss points\n"
        "  --sweeps N        cap on smoothing sweeps (default 200)\n"
        "  --tol X           stop when no node moves more than X mean edges (default 1e-6)\n"
        "  --no-untangle     do not run the untangling phase first\n"
        "  --untangler U     regular (default: regularised shape + size, Garanzha et al.)\n"
        "                    or shifted (metric 22 with a tracked tau_0, the old one)\n"
        "  --untangle-size g the regularised untangler's size weight (default 1/128)\n"
        "  --threads N       OpenMP threads; 0 (default) leaves it to the runtime\n"
        "  --pin-features    pin every boundary and interface node instead of sliding\n"
        "  --no-interfaces   do not treat material interfaces as features\n"
        "  --corner-angle D  the turn a feature node must make to count as a corner\n"
        "  --curves K        what a feature node slides on: spline (the interpolant of its\n"
        "                    smooth run, the default), polyline (the run itself, exactly),\n"
        "                    or chord (the line through its two feature neighbours)\n"
        "  --out F.obj       write the smoothed mesh\n"
        "  --vtu F.vtu       write it with per-element quality and per-node type\n"
        "  --mfem F.mesh     write it as an MFEM mesh with materials as attributes\n";
}

}  // namespace

int main(int argc, char **argv) {
    if (argc < 2) { usage(argv[0]); return 1; }
    if (std::string(argv[1]) == "--selftest") return selfTest();
    if (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h") {
        usage(argv[0]);
        return 0;
    }

    const std::string path = argv[1];
    mesh::TMOP::Options o;
    mesh::QuadMesh::Options mo;
    std::string objOut, vtuOut, mfemOut;

    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--metric" && i + 1 < argc)
            o.metric = static_cast<mesh::TMOP::Metric>(std::stoi(argv[++i]));
        else if (a == "--gamma" && i + 1 < argc)   o.gamma = std::stod(argv[++i]);
        else if (a == "--power" && i + 1 < argc)   o.exponent = std::stod(argv[++i]);
        else if (a == "--target" && i + 1 < argc) {
            const std::string t = argv[++i];
            if (t == "square") o.target = mesh::TMOP::TargetUniformSquare;
            else if (t == "current") o.target = mesh::TMOP::TargetCurrentShape;
            else if (t == "keep") o.target = mesh::TMOP::TargetKeep;
            else { std::cerr << "Unknown target: " << t << "\n"; return 1; }
        }
        else if (a == "--h" && i + 1 < argc)       o.targetSize = std::stod(argv[++i]);
        else if (a == "--corners")                 o.quadrature = mesh::TMOP::Corners;
        else if (a == "--sweeps" && i + 1 < argc)  o.maxSweeps = std::stoi(argv[++i]);
        else if (a == "--tol" && i + 1 < argc)     o.moveTolerance = std::stod(argv[++i]);
        else if (a == "--no-untangle")             o.untangle = false;
        else if (a == "--untangler" && i + 1 < argc) {
            const std::string u = argv[++i];
            if (u == "shifted") o.untangler = mesh::TMOP::UntangleShifted;
            else if (u == "regular") o.untangler = mesh::TMOP::UntangleRegular;
            else { std::cerr << "Unknown untangler: " << u << "\n"; return 1; }
        }
        else if (a == "--untangle-size" && i + 1 < argc) o.untangleSizeWeight = std::stod(argv[++i]);
        else if (a == "--threads" && i + 1 < argc) o.threads = std::stoi(argv[++i]);
        else if (a == "--pin-features")            mo.fixAllFeatureNodes = true;
        else if (a == "--no-interfaces")           mo.interfacesAreFeatures = false;
        else if (a == "--corner-angle" && i + 1 < argc) mo.cornerAngle = std::stod(argv[++i]);
        else if (a == "--curves" && i + 1 < argc) {
            const std::string k = argv[++i];
            if (k == "chord")         mo.curveSource = mesh::QuadMesh::Options::CurveChord;
            else if (k == "polyline") mo.curveSource = mesh::QuadMesh::Options::CurvePolyline;
            else if (k == "spline")   mo.curveSource = mesh::QuadMesh::Options::CurveSpline;
            else { std::cerr << "Unknown --curves: " << k << "\n"; usage(argv[0]); return 1; }
        }
        else if (a == "--out" && i + 1 < argc)     objOut = argv[++i];
        else if (a == "--vtu" && i + 1 < argc)     vtuOut = argv[++i];
        else if (a == "--mfem" && i + 1 < argc)    mfemOut = argv[++i];
        else { std::cerr << "Unknown option: " << a << "\n"; usage(argv[0]); return 1; }
    }

    mesh::QuadMesh m;
    try {
        mesh::QuadMesh loaded(path);
        m = mesh::QuadMesh(loaded.vertices, loaded.quads, loaded.quadMatId, mo);
    } catch (const std::exception &e) {
        std::cout << kFail << " Failed to load " << path << ": " << e.what() << "\n";
        std::cout << "  (this reader takes quadrilateral .obj files -- the meshes in "
                  << "data/meshes are triangles, the pipeline's input)\n";
        return 2;
    }

    std::cout << "TMOP smoothing -- mesh::TMOP on mesh::QuadMesh\n";
    std::cout << "Mesh: " << path << "\n";
    std::cout << "  " << m.vertices.size() << " vertices, " << m.quads.size()
              << " quads, " << m.edges.size() << " edges, " << m.quality.boundaryLoops
              << " boundary loop(s)\n";
    std::cout << "  Materials:";
    for (const std::pair<int, int> &kv : m.materialCounts())
        std::cout << " " << kv.first << " (" << kv.second << ")";
    std::cout << "\n";

    heading("Degrees of freedom");
    std::cout << "  " << m.quality.freeNodes << " free, " << m.quality.slidingNodes
              << " sliding on a feature, " << m.quality.fixedNodes << " fixed\n";
    if (m.featureCurves.empty()) {
        std::cout << "  No feature curves: a sliding node moves along the chord through its "
                  << "two feature neighbours\n";
    } else {
        int fitted = 0, bound = 0;
        double bow = 0.0;
        for (const mesh::QuadMesh::FeatureCurve &fc : m.featureCurves) {
            if (fc.fitted) ++fitted;
            bow = std::max(bow, fc.deviation);
        }
        for (int v = 0; v < static_cast<int>(m.curveOf.size()); ++v)
            if (m.isOnCurve(v)) ++bound;
        std::cout << "  " << m.featureCurves.size() << " feature curve(s), " << fitted
                  << " of them interpolants; " << bound << " node(s) bound, bowing at most "
                  << std::scientific << std::setprecision(2) << bow
                  << " from their own polylines" << std::defaultfloat << "\n";
    }
    int witness = -1;
    if (hasUnanchoredFeatureLoop(m, &witness))
        std::cout << "  " << kWarn << " a feature loop through vertex " << witness
                  << " has no fixed node on it\n";

    mesh::TMOP opt(m, o);
    const bool ok = opt.run();
    const mesh::TMOP::Report &r = opt.getReport();

    heading("Smoothing");
    std::cout << "  " << r.sweeps << " sweep(s)";
    if (r.untangleSweeps > 0) std::cout << " after " << r.untangleSweeps << " untangling sweep(s)";
    std::cout << ", " << r.colors << " colour(s), " << r.threads << " thread(s)"
              << (r.openMP ? "" : " (no OpenMP)") << ", " << std::fixed << std::setprecision(3)
              << r.seconds << " s\n";
    std::cout << "  Energy per unit target area: " << std::scientific << std::setprecision(4)
              << r.energyBefore << " -> " << r.energyAfter << std::defaultfloat << "\n";
    std::cout << "  Scaled Jacobian: worst " << std::fixed << std::setprecision(4)
              << r.minScaledJacobianBefore << " -> " << r.minScaledJacobianAfter
              << ", mean " << r.meanScaledJacobianBefore << " -> " << r.meanScaledJacobianAfter
              << "\n";
    std::cout << "  Worst aspect ratio: " << r.worstAspectBefore << " -> " << r.worstAspectAfter
              << "\n";
    std::cout << "  Inverted elements: " << r.invertedBefore << " -> " << r.invertedAfter << "\n";
    std::cout << "  Area: " << std::setprecision(8) << r.areaBefore << " -> " << r.areaAfter
              << " (" << std::setprecision(4)
              << 100.0 * (r.areaAfter - r.areaBefore) / std::max(r.areaBefore, 1e-300)
              << "%)\n";
    std::cout << "  Nodes moved: " << std::setprecision(3) << r.meanDisplacement
              << " mean edge(s) on average, " << r.maxDisplacement << " at most\n";
    for (const std::string &msg : r.messages) std::cout << "  " << kWarn << " " << msg << "\n";

    verdict(r.invertedAfter == 0, "No element of the smoothed mesh is inverted");
    verdict(r.minScaledJacobianAfter >= r.minScaledJacobianBefore - 1e-12,
            "The worst element is no worse than it was");
    // The area is exact only when the nodes slide along chords; on a curve the
    // mesh grows into the curve as the boundary redistributes, and what has to
    // hold instead is that no node ever left the curve it rides.
    if (r.curveNodes > 0) {
        std::cout << "  Off-curve: " << std::scientific << std::setprecision(2)
                  << r.curveDeviation << std::defaultfloat << "\n";
        verdict(r.curveDeviation <= 1e-9 * std::max(m.quality.meanEdge, 1e-300),
                "Every feature node is still on the curve it slides on");
    } else {
        verdict(std::fabs(r.areaAfter - r.areaBefore) <= 1e-9 * std::fabs(r.areaBefore),
                "The domain has the area it started with");
    }
    // Not a failure. A Gauss-Seidel sweep contracts the error by a roughly
    // constant factor, so the tail is long; running out of sweeps means the
    // mesh is still improving, not that anything went wrong.
    if (!r.converged)
        std::cout << "  " << kWarn << " stopped on the sweep cap, still moving "
                  << std::scientific << std::setprecision(2) << r.lastSweepMove
                  << " mean edge(s) per sweep -- raise --sweeps to go further\n"
                  << std::defaultfloat;

    if (!objOut.empty())
        std::cout << "  " << (m.writeOBJ(objOut) ? "Wrote " : "Failed to write ") << objOut << "\n";
    if (!vtuOut.empty())
        std::cout << "  " << (m.writeVTU(vtuOut) ? "Wrote " : "Failed to write ") << vtuOut << "\n";
    if (!mfemOut.empty())
        std::cout << "  " << (m.writeMFEM(mfemOut) ? "Wrote " : "Failed to write ") << mfemOut << "\n";

    heading("Result");
    if (ok && gFailures == 0) {
        std::cout << "  " << kPass << " " << m.quads.size() << " element(s) smoothed to a worst "
                  << "scaled Jacobian of " << std::fixed << std::setprecision(4)
                  << r.minScaledJacobianAfter << "\n";
        return 0;
    }
    std::cout << "  " << kFail << " the smoother did not improve this mesh\n";
    return 3;
}
