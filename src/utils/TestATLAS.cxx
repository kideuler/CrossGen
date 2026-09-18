// Utility to run ATLAS -- square-transport blocking, Stages 1 to 6 of
// docs/square_transport_2d_theory_and_implementation.md -- on a mesh.
//
//   TestATLAS <mesh.obj> [options]
//   TestATLAS --selftest
//
// Reports each stage and exits non-zero unless a valid conforming blocking came
// out. A valid blocking may still be the fine fallback of Sec. 14 -- every
// carrier cell its own block -- and the report says so rather than calling it
// a failure: the spec's guarantee is validity, and coarseness is what the
// search is for.
//
//   --selftest   the regression geometries of Sec. 13.3 whose answer is known
//                before the code runs: the integer transport algebra; one
//                triangle, where Sec. 3.1's corner determinants are exactly
//                2A (1/4, 1/6, 1/12, 1/6) and the blocking is three quads; a
//                rectangle (one block); a disk (a five-block O-grid); an
//                annulus (four sectors, no centre); an L (three blocks, the
//                reentrant corner a macrovertex); two materials side by side
//                (two blocks sharing the interface as a macro edge); the
//                certificate refusing a periodic band and an L-tromino; and
//                Stage 6 alone lowering the valence defect of a raw carrier
//                without ever leaving it invalid.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "ATLAS/ATLAS.hxx"
#include "TestHelper.hxx"
#include "mesh/QuadMesh.hxx"

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
// Meshes for the self-test
// ---------------------------------------------------------------------------

// Triangulate polygon loops with Triangle. Each loop is a list of corners;
// every side is subdivided at spacing about h so the boundary is resolved
// before Triangle sees it. types: 0 outer, 1 hole.
std::shared_ptr<Mesh> meshLoops(const std::vector<std::vector<Point>> &loops,
                                const std::vector<int> &types, double h) {
    using namespace triangle_wrapper;
    TriangleMesher2D::MeshInput input;
    for (size_t l = 0; l < loops.size(); ++l) {
        const std::vector<Point> &L = loops[l];
        const int first = static_cast<int>(input.vertlist.size());
        for (size_t k = 0; k < L.size(); ++k) {
            const Point a = L[k], b = L[(k + 1) % L.size()];
            const int n = std::max(1, static_cast<int>(std::lround(normP(b - a) / h)));
            for (int i = 0; i < n; ++i) {
                const Point p = a + (b - a) * (static_cast<double>(i) / n);
                input.vertlist.push_back({p[0], p[1]});
            }
        }
        const int last = static_cast<int>(input.vertlist.size());
        std::vector<std::array<int, 2>> segs;
        for (int i = first; i < last; ++i) segs.push_back({i, i + 1 < last ? i + 1 : first});
        input.segment_loops.push_back(segs);
        input.type.push_back(types[l]);
    }
    input.h = h;
    TriangleMesher2D::Options o;
    o.min_angle_degrees = 20.0;
    TriangleMesher2D mesher(o);
    auto result = mesher.triangulate(input);
    std::vector<Point> verts(result.verts.size());
    for (size_t i = 0; i < result.verts.size(); ++i) verts[i] = {result.verts[i][0], result.verts[i][1]};
    std::vector<Triangle> tris(result.triangles.size());
    for (size_t i = 0; i < result.triangles.size(); ++i) {
        tris[i] = {result.triangles[i][0], result.triangles[i][1], result.triangles[i][2]};
    }
    return std::make_shared<Mesh>(verts, tris);
}

std::vector<Point> circlePoints(double cx, double cy, double r, int n) {
    std::vector<Point> p;
    for (int i = 0; i < n; ++i) {
        const double t = 2.0 * M_PI * i / n;
        p.push_back({cx + r * std::cos(t), cy + r * std::sin(t)});
    }
    return p;
}

// A structured triangulation of [0,W] x [0,H], nx by ny squares, each cut on
// alternating diagonals; material 1 left of x = split, 2 right of it.
std::shared_ptr<Mesh> gridMesh(int nx, int ny, double W, double H, double split) {
    std::vector<Point> v;
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) v.push_back({W * i / nx, H * j / ny});
    }
    auto id = [&](int i, int j) { return i + (nx + 1) * j; };
    std::vector<Triangle> t;
    std::vector<int> mat;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const int a = id(i, j), b = id(i + 1, j), c = id(i + 1, j + 1), d = id(i, j + 1);
            const int m = (W * (i + 0.5) / nx) < split ? 1 : 2;
            if ((i + j) % 2 == 0) {
                t.push_back({a, b, c});
                t.push_back({a, c, d});
            } else {
                t.push_back({a, b, d});
                t.push_back({b, c, d});
            }
            mat.push_back(m);
            mat.push_back(m);
        }
    }
    return std::make_shared<Mesh>(v, t, mat);
}

ATLAS::Options quietOptions() {
    ATLAS::Options o;
    o.rewrite.timeBudget = 20.0;
    o.timeBudget = 60.0;
    return o;
}

void summary(const ATLAS &pipe) {
    const ATLAS::Status &st = pipe.getStatus();
    std::cout << "  " << st.triangles << " triangle(s) -> " << st.initialCells << " carrier cell(s)";
    if (pipe.hasTemplates()) {
        const ExplicitTemplates::Report &tr = pipe.getTemplateReport();
        std::cout << "; Stage 3 templated " << tr.templated << " of " << tr.regions << " region(s)";
    }
    if (pipe.hasCover()) {
        const BlockCover::Report &cr = pipe.getCover().getReport();
        std::cout << "; " << cr.blocks << " block(s) over " << pipe.getCarrier().numCells() << " cell(s)";
    }
    std::cout << "\n";
    if (pipe.hasTemplates()) {
        for (const ExplicitTemplates::Attempt &a : pipe.getTemplateReport().attempts) {
            std::cout << "     region of " << a.cells << " cell(s), " << a.corners << " corner(s), "
                      << a.holes << " hole(s): " << (a.accepted ? a.family : std::string("rejected"));
            if (a.accepted) std::cout << ", " << a.blocks << " block(s), " << a.maps;
            else std::cout << " -- " << a.reason;
            std::cout << "\n";
        }
    }
}

// ---------------------------------------------------------------------------
// selfTest()
// ---------------------------------------------------------------------------
int selfTest() {
    std::cout << "ATLAS self-test (docs/square_transport_2d_theory_and_implementation.md)\n";

    // -----------------------------------------------------------------
    heading("Case 1  The integer transport algebra (Sec. 2.1)");
    // -----------------------------------------------------------------
    {
        // Every pair of glued sides: g_qr o g_rq = I, and g_qr carries r's
        // copy of the shared edge onto q's, endpoints swapped.
        int bad = 0;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                const SquareTransport g = transportAcross(i, j), h = transportAcross(j, i);
                if (!g.compose(h).isIdentity() || !h.compose(g).isIdentity()) ++bad;
                if (g.apply(squareCorner(j)) != squareCorner(i + 1)) ++bad;
                if (g.apply(squareCorner(j + 1)) != squareCorner(i)) ++bad;
                if (g != h.inverse()) ++bad;
                // r's square lands outside q's: its centre (doubled) is not
                // inside (0,2)^2.
                const IPoint c = g.apply({0, 0}), d = g.apply({1, 1});
                const long long cx = c[0] + d[0], cy = c[1] + d[1];
                if (cx > 0 && cx < 2 && cy > 0 && cy < 2) ++bad;
            }
        }
        check(bad == 0, "g_rq = g_qr^-1 and the shared edge maps onto itself, for all 16 side pairs");
        // Two squares side by side in one frame: q's side 1 meets r's side 3,
        // and the transport is the translation by (1, 0).
        const SquareTransport g = transportAcross(1, 3);
        check(g.rot == 0 && g.tx == 1 && g.ty == 0, "q's side 1 against r's side 3 is the translation (1, 0)");
        // Associativity and the inverse on a few compositions.
        const SquareTransport a(1, 2, -3), b(3, -1, 5), c(2, 4, 0);
        check(a.compose(b).compose(c) == a.compose(b.compose(c)), "composition is associative");
        check(a.compose(a.inverse()).isIdentity(), "a o a^-1 = I");
    }

    // -----------------------------------------------------------------
    heading("Case 2  One triangle: Sec. 3.1's determinants, three blocks");
    // -----------------------------------------------------------------
    {
        auto mesh = std::make_shared<Mesh>(std::vector<Point>{{0.0, 0.0}, {2.0, 0.0}, {0.3, 1.5}},
                                           std::vector<Triangle>{{0, 1, 2}});
        ATLAS pipe(mesh, quietOptions());
        const bool ok = pipe.run();
        summary(pipe);
        const SquareCarrier::Report &cr = pipe.getInitialCarrierReport();
        std::cout << "  corner det / 2A in [" << std::setprecision(10) << cr.minCornerDetRatio << ", "
                  << cr.maxCornerDetRatio << "]" << std::setprecision(6) << "\n";
        check(cr.cells == 3, "N_Q = 3 N_T");
        check(cr.jacobianBoundViolations == 0,
              "every corner determinant is exactly 2A (1/4, 1/6, 1/12, 1/6)");
        check(std::fabs(cr.minCornerDetRatio - 1.0 / 12.0) < 1e-12 &&
                  std::fabs(cr.maxCornerDetRatio - 0.25) < 1e-12,
              "det DX spans exactly [1/12, 1/4] of 2A");
        check(cr.eulerHolds && cr.eulerLHS == 4, "sum (4 - q) + sum (2 - q) = 4 on a disk");
        check(cr.boundaryParityEven, "the split doubles every boundary edge: even total");
        check(ok, "the pipeline returns a valid blocking");
        if (pipe.hasCover()) check(pipe.getCover().getReport().blocks == 3, "three blocks, round one centre");
    }

    // -----------------------------------------------------------------
    heading("Case 3  Rectangle: one block (Sec. 7.1)");
    // -----------------------------------------------------------------
    std::unique_ptr<ATLAS> rect;
    {
        auto mesh = TestHelper::createBox(0.0, 0.0, 2.0, 1.0, 0.1);
        rect = std::make_unique<ATLAS>(mesh, quietOptions());
        const bool ok = rect->run();
        summary(*rect);
        check(rect->getInitialCarrierReport().valid, "Stage 2: the carrier validates");
        check(rect->getInitialCarrierReport().jacobianBoundViolations == 0, "Sec. 3.1's bound on every cell");
        check(rect->getTemplateReport().templated == 1, "Stage 3 replaced the region");
        check(ok, "the pipeline returns a valid blocking");
        if (rect->hasCover()) {
            const BlockCover::Report &br = rect->getCover().getReport();
            check(br.blocks == 1, "one block");
            check(br.macroVertices == 4, "four macrovertices, the four corners");
            const RectangleCertifier::Certificate &c = rect->getCover().getBlocks()[0].cert;
            check(c.nu * c.nv == rect->getCarrier().numCells(), "the block is the whole carrier, n_u x n_v");
        }
    }

    // -----------------------------------------------------------------
    heading("Case 4  Disk: a five-block O-grid (Sec. 7.3)");
    // -----------------------------------------------------------------
    {
        auto mesh = TestHelper::createCircle(0.0, 0.0, 1.0, 0.08);
        ATLAS pipe(mesh, quietOptions());
        const bool ok = pipe.run();
        summary(pipe);
        check(ok, "the pipeline returns a valid blocking");
        if (pipe.hasCover()) {
            const BlockCover::Report &br = pipe.getCover().getReport();
            check(br.blocks == 5, "five blocks: one core, four shells");
            check(br.irregularMacroVertices == 4, "four three-valent macrovertices, the core's corners");
        }
    }

    // -----------------------------------------------------------------
    heading("Case 5  Annulus: four sectors, the hole kept (Sec. 7.4)");
    // -----------------------------------------------------------------
    std::unique_ptr<ATLAS> ring;
    {
        auto mesh = meshLoops({circlePoints(0.0, 0.0, 1.0, 72), circlePoints(0.05, 0.0, 0.45, 36)},
                              {0, 1}, 0.09);
        ring = std::make_unique<ATLAS>(mesh, quietOptions());
        const bool ok = ring->run();
        summary(*ring);
        check(ring->getInitialCarrierReport().holes == 1, "the carrier has one hole");
        check(ring->getInitialCarrierReport().eulerRHS == 0, "4 (1 - h) = 0");
        check(ok, "the pipeline returns a valid blocking");
        if (ring->hasTemplates() && !ring->getTemplateReport().attempts.empty()) {
            const ExplicitTemplates::Attempt &a = ring->getTemplateReport().attempts.back();
            check(a.accepted && a.family == "annulus" && a.blocks == 4, "Stage 3 lays out four sectors");
        }
        if (ring->hasCover()) {
            // Stage 5 then merges two neighbouring sectors: the half-annulus
            // certifies and the objective prefers it. Three is the minimum --
            // two blocks would share both radial sides, which Sec. 6 forbids.
            const BlockCover::Report &br = ring->getCover().getReport();
            check(br.blocks == 3, "Stage 5 reaches the minimum, three blocks");
            check(br.multipleContacts == 0, "no two blocks share two sides");
            check(br.irregularMacroVertices == 0, "no singular macrovertex anywhere");
        }
    }

    // -----------------------------------------------------------------
    heading("Case 6  L: the reentrant corner is a macrovertex (Sec. 13.3)");
    // -----------------------------------------------------------------
    {
        auto mesh = meshLoops({{{0.0, 0.0}, {2.0, 0.0}, {2.0, 1.0}, {1.0, 1.0}, {1.0, 2.0}, {0.0, 2.0}}},
                              {0}, 0.1);
        ATLAS pipe(mesh, quietOptions());
        const bool ok = pipe.run();
        summary(pipe);
        check(ok, "the pipeline returns a valid blocking");
        if (pipe.hasCover()) {
            const BlockCover::Report &br = pipe.getCover().getReport();
            check(br.blocks == 3, "three blocks");
            check(br.protectedNotMacro == 0, "every protected corner is a macrovertex");
        }
    }

    // -----------------------------------------------------------------
    heading("Case 7  Two materials: the interface is a shared macro edge");
    // -----------------------------------------------------------------
    {
        auto mesh = gridMesh(20, 10, 2.0, 1.0, 1.0);
        ATLAS pipe(mesh, quietOptions());
        const bool ok = pipe.run();
        summary(pipe);
        check(pipe.getDomain().getReport().interfaceLandings == 2, "Stage 1 finds the interface's two landings");
        check(ok, "the pipeline returns a valid blocking");
        if (pipe.hasCover()) {
            const BlockCover::Report &br = pipe.getCover().getReport();
            check(br.blocks == 2, "two blocks, one per material");
            check(br.interfaceInside == 0, "no block straddles the interface");
            int interfaceEdges = 0;
            for (const BlockCover::MacroEdge &me : pipe.getCover().getMacroEdges()) interfaceEdges += me.interface;
            check(interfaceEdges == 1, "the interface is exactly one shared macro edge");
        }
    }

    // -----------------------------------------------------------------
    heading("Case 8  The certificate refuses what weaker tests accept (Sec. 5.2)");
    // -----------------------------------------------------------------
    {
        // The annulus after Stage 3 is a perfectly regular grid everywhere --
        // every interior vertex four-valent, zero holonomy on every cell --
        // and the band is still not a block.
        if (ring && ring->hasTemplates()) {
            const SquareCarrier &C = ring->getTemplatedCarrier();
            RectangleCertifier R(C, RectangleCertifier::Options());
            std::vector<int> all(C.numCells());
            for (int q = 0; q < C.numCells(); ++q) all[q] = q;
            RectangleCertifier::Certificate cert;
            RectangleCertifier::Witness w;
            const bool accepted = R.certify(all, cert, w);
            std::cout << "  the whole band: " << (accepted ? "accepted" : RectangleCertifier::failureName(w.kind)) << "\n";
            check(!accepted, "a periodic band is refused");
            check(C.getReport().irregularInterior == 0, "... although no interior vertex of it is irregular");
        }
        // An L-tromino out of the rectangle's single block: three cells whose
        // bounding box has four squares.
        if (rect && rect->hasCover() && rect->getCover().getBlocks().size() == 1) {
            const SquareCarrier &C = rect->getCarrier();
            const RectangleCertifier::Certificate &b = rect->getCover().getBlocks()[0].cert;
            RectangleCertifier R(C, RectangleCertifier::Options());
            std::vector<int> tromino{b.cells[0], b.cells[1], b.cells[b.nu]};
            RectangleCertifier::Certificate cert;
            RectangleCertifier::Witness w;
            const bool accepted = R.certify(tromino, cert, w);
            std::cout << "  an L-tromino: " << (accepted ? "accepted" : RectangleCertifier::failureName(w.kind)) << "\n";
            check(!accepted && w.kind == RectangleCertifier::Failure::MissingSquare, "an L-tromino is refused: missing square");
            std::vector<int> two{b.cells[0], b.cells[1], b.cells[b.nu], b.cells[b.nu + 1]};
            check(R.certify(two, cert, w) && cert.nu == 2 && cert.nv == 2, "the 2 x 2 box it came from is accepted");
        }
    }

    // -----------------------------------------------------------------
    heading("Case 9  Stage 6 alone, on a raw carrier");
    // -----------------------------------------------------------------
    {
        // No templates: every irregular vertex the triangulation left is still
        // there, and Stage 6 is what removes them. What is known in advance is
        // not the answer but its invariants: the carrier validates after every
        // round, and neither the defect nor the best objective goes up.
        auto mesh = TestHelper::createCircle(0.0, 0.0, 1.0, 0.15);
        ATLAS::Options o = quietOptions();
        o.runTemplates = false;
        o.rewriteRounds = 6;
        ATLAS pipe(mesh, o);
        const bool ok = pipe.run();
        const ATLAS::Status &st = pipe.getStatus();
        for (const ATLAS::Round &r : st.history) {
            std::cout << "  round " << r.round << ": " << r.cells << " cells, defect " << r.defect
                      << ", " << r.blocks << " block(s), " << r.committed << " rewrite(s)\n";
        }
        bool monotone = true;
        for (size_t k = 1; k < st.history.size(); ++k) {
            if (st.history[k].defect > st.history[k - 1].defect) monotone = false;
        }
        check(ok, "the pipeline returns a valid blocking");
        check(st.rewritesValid, "the carrier validated after every Stage 6 round");
        check(monotone, "the valence defect never rises");
        check(st.history.size() > 1 && st.history.back().defect < st.history.front().defect,
              "Stage 6 lowered the defect");
        check(st.history.size() > 1 && st.bestBlocks < st.history.front().blocks,
              "... and the best blocking is coarser than the first");
    }

    heading("Self-test result");
    if (failures == 0) {
        std::cout << "  " << kPass << " Every check held.\n";
    } else {
        std::cout << "  " << kFail << " " << failures << " check(s) failed.\n";
    }
    return failures == 0 ? 0 : 1;
}

// ---------------------------------------------------------------------------
// Reporting for a single model
// ---------------------------------------------------------------------------
void printCarrier(const SquareCarrier::Report &r) {
    std::cout << "  " << r.cells << " cells, " << r.vertices << " vertices, " << r.edges << " edges ("
              << r.boundaryEdges << " on dS, " << r.interfaceEdges << " on interfaces)\n";
    std::cout << "  " << r.components << " component(s), " << r.holes << " hole(s); "
              << r.protectedVertices << " protected, " << r.designatedVertices << " designated vertex/vertices\n";
    std::cout << "  Scaled Jacobian " << std::fixed << std::setprecision(4) << r.minScaledJacobian
              << " worst, " << r.meanScaledJacobian << " mean" << std::defaultfloat << "\n";
    std::cout << "  Irregular: " << r.irregularInterior << " interior, " << r.irregularBoundary
              << " on dS; total valence defect " << r.totalDefect << "\n";
    std::cout << "  Sec. 4: sum_int (4 - q) + sum_bdy (2 - q) = " << r.eulerLHS << ", 4 sum (1 - h) = "
              << r.eulerRHS << "\n";
}

void usage(const char *prog) {
    std::cout << "Usage: " << prog << " <mesh.obj> [options]\n"
              << "       " << prog << " --selftest\n\n"
              << "Stage 1  the domain\n"
              << "  --corner-angle <deg>  a boundary turn past this is a protected corner (default 30)\n"
              << "  --kink-angle <deg>    the same along an interface                   (default 30)\n"
              << "  --no-interfaces       ignore the material tags\n\n"
              << "Stage 3  explicit replacements\n"
              << "  --no-templates        skip Stage 3 entirely\n"
              << "  --no-sections         no sectioned quadrilaterals (Secs. 7.1, 7.2)\n"
              << "  --no-stars            no three- or five-block stars\n"
              << "  --no-ogrid            no O-grids (Sec. 7.3)\n"
              << "  --no-annulus          no annuli (Sec. 7.4)\n"
              << "  --reflex-angle <deg>  a corner past pi by this emits sections      (default 20)\n"
              << "  --core <f>            O-grid core at this fraction of the radius   (default 0.5)\n"
              << "  --no-split            never subdivide a boundary segment (Sec. 8.4)\n\n"
              << "Stage 4  rectangles\n"
              << "  --no-merges           no unions of adjacent base patches\n"
              << "  --no-growth           no boxes grown row by row\n"
              << "  --max-candidates <n>  candidate budget                             (default 200000)\n\n"
              << "Stage 5  cover\n"
              << "  --alpha <a>           weight on block distortion D(P)               (default 0.05)\n"
              << "  --beta <b>            weight on side complexity C(P)                (default 0.05)\n"
              << "  --multi-contact       allow two blocks to share two sides\n"
              << "  --sweeps <n>          improvement sweeps                            (default 20)\n\n"
              << "Stage 6  cavity rewrites\n"
              << "  --no-rewrite          skip Stage 6\n"
              << "  --rounds <n>          rounds of Stages 4-6                          (default 12)\n"
              << "  --rings <n>           starting cavity radius, in rings              (default 3)\n"
              << "  --max-rings <n>       widest radius a stalled round may grow to     (default 8)\n"
              << "  --rewrite-floor <s>   scaled-Jacobian floor for a rewrite           (default 0.2)\n"
              << "  --no-grid-rewrites    no single-grid refills\n"
              << "  --no-star-rewrites    no star refills\n"
              << "  --time <s>            budget for the whole Stage 4-6 loop           (default 300)\n\n"
              << "Output\n"
              << "  --initial <f.obj>     the Stage 2 carrier\n"
              << "  --templated <f.obj>   the carrier after Stage 3\n"
              << "  --carrier <f.obj>     the carrier the best blocking lives on\n"
              << "  --blocks <f.obj>      the macro edges of the best blocking, as polylines\n"
              << "  --mfem <f.mesh>       the carrier of the best blocking for MFEM, material per element\n";
}

} // namespace

int main(int argc, char **argv) {
    if (argc < 2) { usage(argv[0]); return 1; }
    const std::string first = argv[1];
    if (first == "--selftest") return selfTest();
    if (first == "-h" || first == "--help") { usage(argv[0]); return 0; }

    const std::string path = first;
    ATLAS::Options opts;
    std::string initialOut, templatedOut, carrierOut, blocksOut, mfemOut;

    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--corner-angle" && i + 1 < argc)        opts.domain.cornerAngle = std::stod(argv[++i]);
        else if (a == "--kink-angle" && i + 1 < argc)     opts.domain.kinkAngle = std::stod(argv[++i]);
        else if (a == "--no-interfaces")                  opts.domain.materialInterfaces = false;
        else if (a == "--no-templates")                   opts.runTemplates = false;
        else if (a == "--no-sections")                    opts.templates.sections = false;
        else if (a == "--no-stars")                       opts.templates.stars = false;
        else if (a == "--no-ogrid")                       opts.templates.ogrids = false;
        else if (a == "--no-annulus")                     opts.templates.annuli = false;
        else if (a == "--reflex-angle" && i + 1 < argc)   opts.templates.reflexAngle = std::stod(argv[++i]);
        else if (a == "--core" && i + 1 < argc)           opts.templates.coreFraction = std::stod(argv[++i]);
        else if (a == "--no-split")                       opts.templates.splitBoundary = false;
        else if (a == "--no-merges")                      opts.rectangles.merges = false;
        else if (a == "--no-growth")                      opts.rectangles.growth = false;
        else if (a == "--max-candidates" && i + 1 < argc) opts.rectangles.maxCandidates = std::stoi(argv[++i]);
        else if (a == "--alpha" && i + 1 < argc)          opts.cover.alpha = std::stod(argv[++i]);
        else if (a == "--beta" && i + 1 < argc)           opts.cover.beta = std::stod(argv[++i]);
        else if (a == "--multi-contact")                  opts.cover.singleContact = false;
        else if (a == "--sweeps" && i + 1 < argc)         opts.cover.maxSweeps = std::stoi(argv[++i]);
        else if (a == "--no-rewrite")                     opts.runRewrite = false;
        else if (a == "--rounds" && i + 1 < argc)         opts.rewriteRounds = std::stoi(argv[++i]);
        else if (a == "--rings" && i + 1 < argc)          opts.rewrite.maxRings = std::stoi(argv[++i]);
        else if (a == "--max-rings" && i + 1 < argc)      opts.maxRings = std::stoi(argv[++i]);
        else if (a == "--rewrite-floor" && i + 1 < argc)  opts.rewrite.minScaledJacobian = std::stod(argv[++i]);
        else if (a == "--no-grid-rewrites")               opts.rewrite.grids = false;
        else if (a == "--no-star-rewrites")               opts.rewrite.stars = false;
        else if (a == "--time" && i + 1 < argc)           opts.timeBudget = std::stod(argv[++i]);
        else if (a == "--initial" && i + 1 < argc)        initialOut = argv[++i];
        else if (a == "--templated" && i + 1 < argc)      templatedOut = argv[++i];
        else if (a == "--carrier" && i + 1 < argc)        carrierOut = argv[++i];
        else if (a == "--blocks" && i + 1 < argc)         blocksOut = argv[++i];
        else if (a == "--mfem" && i + 1 < argc)           mfemOut = argv[++i];
        else { std::cerr << "Unknown option: " << a << "\n"; usage(argv[0]); return 1; }
    }

    std::shared_ptr<Mesh> mesh;
    try {
        mesh = std::make_shared<Mesh>(path);
    } catch (const std::exception &e) {
        std::cout << kFail << " Failed to load mesh: " << e.what() << "\n";
        return 2;
    }

    std::cout << "ATLAS -- square-transport blocking (docs/square_transport_2d_theory_and_implementation.md)\n";
    std::cout << "Mesh: " << path << "\n";
    std::cout << "  " << mesh->vertices.size() << " vertices, " << mesh->edges.size() << " edges, "
              << mesh->triangles.size() << " triangles, " << mesh->boundaryEdges.size()
              << " boundary edges\n";

    ATLAS pipe(mesh, opts);
    bool ok = false;
    try {
        ok = pipe.run();
    } catch (const std::exception &e) {
        std::cout << kFail << " Pipeline threw: " << e.what() << "\n";
        return 3;
    }
    const ATLAS::Status &st = pipe.getStatus();

    // ---------------------------------------------------------------------
    heading("Stage 1  Validate and tag the planar domain (Sec. 1.1)");
    {
        const PlanarDomain::Report &dr = pipe.getDomain().getReport();
        std::cout << "  " << dr.components << " component(s), " << dr.loops << " boundary loop(s), "
                  << dr.holes << " hole(s); area " << dr.area << ", boundary length "
                  << dr.boundaryLength << "\n";
        std::cout << "  Smallest triangle angle " << std::fixed << std::setprecision(2)
                  << dr.minAngleDegrees << " deg" << std::defaultfloat << "\n";
        std::cout << "  " << dr.materials << " material(s), " << dr.interfaceEdges << " interface edge(s)\n";
        std::cout << "  Protected: " << dr.boundaryCorners << " boundary corner(s), " << dr.interfaceJunctions
                  << " junction(s), " << dr.interfaceLandings << " landing(s), " << dr.interfaceKinks
                  << " kink(s), " << dr.interfaceDangling << " dangling end(s)\n";
        verdict(dr.degenerateTriangles == 0 && dr.invertedTriangles == 0,
                "Every triangle is nondegenerate and counter-clockwise");
        verdict(dr.nonManifoldEdges == 0 && dr.nonManifoldVertices == 0,
                "Manifold: no edge in three triangles, no pinched vertex");
        verdict(dr.componentsWithoutOuterLoop == 0 && dr.eulerMismatches == 0,
                "Each component has one outer loop and V - E + F = 1 - h");
        for (const std::string &m : dr.messages) warn(m);
        if (!st.domainValid) {
            heading("Result");
            std::cout << "  " << kFail << " The domain does not satisfy Sec. 1.1's hypotheses.\n";
            return 4;
        }
    }

    // ---------------------------------------------------------------------
    heading("Stage 2  The guaranteed carrier (Sec. 3)");
    {
        const SquareCarrier::Report &r = pipe.getInitialCarrierReport();
        printCarrier(r);
        std::cout << "  Corner determinants / 2A_T in [" << std::setprecision(8) << r.minCornerDetRatio
                  << ", " << r.maxCornerDetRatio << "]" << std::setprecision(6) << "\n";
        verdict(r.cells == 3 * r.sourceTriangles, "N_Q = 3 N_T");
        verdict(r.jacobianBoundViolations == 0 && r.nonConvexCells == 0,
                "Sec. 3.1: every corner determinant is 2A (1/4, 1/6, 1/12, 1/6)");
        verdict(r.transportMismatches == 0, "Sec. 2.1: g_rq = g_qr^-1 on every shared edge");
        verdict(r.angleSumViolations == 0 && r.areaError < 1e-9,
                "Sec. 9.2: every one-ring winds once and the cells cover the domain");
        verdict(r.eulerHolds, "Sec. 4: the Euler identity holds");
        verdict(r.boundaryParityEven, "Sec. 4: an even number of boundary edges");
        verdict(r.boundaryPreservationErrors == 0 && r.missingBoundaryVertices == 0,
                "Sec. 1.1: every boundary segment is preserved");
        std::cout << "  The singleton blocking -- " << r.cells << " blocks -- is the incumbent.\n";
        if (!initialOut.empty()) {
            if (pipe.getInitialCarrier().writeOBJ(initialOut)) std::cout << "  Wrote the carrier to " << initialOut << "\n";
            else warn("Failed to write " + initialOut);
        }
        if (!st.carrierValid) {
            heading("Result");
            std::cout << "  " << kFail << " The initial carrier does not validate.\n";
            return 5;
        }
    }

    // ---------------------------------------------------------------------
    if (pipe.hasTemplates()) {
        heading("Stage 3  Explicit replacements (Sec. 7)");
        const ExplicitTemplates::Report &tr = pipe.getTemplateReport();
        std::cout << "  " << tr.templated << " of " << tr.regions << " region(s) replaced, "
                  << std::fixed << std::setprecision(1) << 100.0 * tr.templatedArea << "% of the area"
                  << std::defaultfloat << "; " << tr.blocks << " block(s) planned; "
                  << tr.boundarySplits << " boundary point(s) inserted\n";
        std::cout << "  Cells " << tr.cellsBefore << " -> " << tr.cellsAfter << "\n";
        int shown = 0;
        for (const ExplicitTemplates::Attempt &a : tr.attempts) {
            if (shown++ >= 40) { std::cout << "     ...\n"; break; }
            std::cout << "     mat " << a.material << ", " << a.cells << " cell(s), " << a.corners
                      << " corner(s), " << a.holes << " hole(s): ";
            if (a.accepted) {
                std::cout << a.family << ", " << a.blocks << " block(s), " << a.newCells << " cell(s), "
                          << a.maps << "; worst scaled Jacobian " << std::fixed << std::setprecision(3)
                          << a.minScaledJacobian << std::defaultfloat;
                if (a.splits > 0) std::cout << "; " << a.splits << " split(s)";
            } else {
                std::cout << "kept (" << a.reason << ")";
            }
            std::cout << "\n";
        }
        verdict(st.templatesValid, "The carrier still validates after every replacement");
        if (!templatedOut.empty()) {
            if (pipe.getTemplatedCarrier().writeOBJ(templatedOut)) std::cout << "  Wrote the templated carrier to " << templatedOut << "\n";
            else warn("Failed to write " + templatedOut);
        }
        if (!st.templatesValid) {
            for (const std::string &m : st.messages) warn(m);
            heading("Result");
            std::cout << "  " << kFail << " A replacement broke the carrier.\n";
            return 6;
        }
    }

    if (!pipe.hasCover()) {
        heading("Result");
        std::cout << "  " << kFail << " No cover was selected.\n";
        return 7;
    }

    // ---------------------------------------------------------------------
    heading("Stages 4-5  First round: certified rectangles and a conforming cover (Secs. 5, 6)");
    {
        const RectangleCertifier::Report &rr = pipe.getFirstRectangles().getReport();
        const BlockCover::Report &br = pipe.getFirstCover().getReport();
        std::cout << "  Forced macrovertices " << rr.forcedVertices << ", designated " << rr.softVertices
                  << ", cuts added " << rr.cutVertices << "\n";
        std::cout << "  Base complex: " << rr.basePatchesAll << " patch(es) with the designated vertices, "
                  << rr.basePatchesHard << " without\n";
        std::cout << "  Candidates: " << rr.candidates << " (" << rr.fromBaseAll << " + " << rr.fromBaseHard
                  << " base, " << rr.fromMerges << " merged, " << rr.fromGrowth << " grown); "
                  << rr.certifyCalls << " certification(s), largest " << rr.largestCandidate << " cell(s)\n";
        if (rr.failures > 0) {
            std::cout << "  Failure witnesses:";
            for (int k = 0; k < static_cast<int>(rr.failuresByKind.size()); ++k) {
                if (rr.failuresByKind[k] > 0) {
                    std::cout << " " << rr.failuresByKind[k] << " "
                              << RectangleCertifier::failureName(static_cast<RectangleCertifier::Failure>(k)) << ";";
                }
            }
            std::cout << "\n";
        }
        verdict(rr.cutsConverged, "Every base patch certified after cutting");
        std::cout << "  Cover: " << br.blocks << " block(s) from '" << br.incumbent << "', " << br.swaps
                  << " swap(s) in " << br.sweeps << " sweep(s)\n";
        verdict(br.valid, "Sec. 13.2: the first cover is conforming");
    }

    // ---------------------------------------------------------------------
    if (pipe.hasRewrite()) {
        heading("Stage 6  Cavity rewrites, Stages 4-6 repeated (Sec. 8)");
        const CavityRewrite::Report &wr = pipe.getRewriteReport();
        std::cout << "  " << std::setw(6) << "round" << std::setw(9) << "cells" << std::setw(9) << "defect"
                  << std::setw(11) << "irregular" << std::setw(12) << "candidates" << std::setw(9)
                  << "blocks" << std::setw(10) << "rewrites" << std::setw(7) << "rings" << std::setw(9)
                  << "sec" << "\n";
        for (const ATLAS::Round &r : st.history) {
            std::cout << "  " << std::setw(6) << r.round << std::setw(9) << r.cells << std::setw(9) << r.defect
                      << std::setw(11) << r.irregular << std::setw(12) << r.candidates << std::setw(9)
                      << r.blocks << std::setw(10) << r.committed << std::setw(7) << r.rings
                      << std::setw(9) << std::fixed
                      << std::setprecision(2) << r.seconds << std::defaultfloat
                      << (r.round == st.bestRound ? "  <- best" : "") << "\n";
        }
        std::cout << "  " << wr.witnesses << " witness(es), " << wr.cavities << " cavity/cavities, "
                  << wr.proposals << " proposal(s), " << wr.realised << " realised, " << wr.certified
                  << " certified, " << wr.committed << " committed (" << wr.gridFills << " grid, "
                  << wr.starFills << " star)\n";
        std::cout << "  Valence defect " << wr.defectBefore << " -> " << wr.defectAfter << ", irregular vertices "
                  << wr.irregularBefore << " -> " << wr.irregularAfter << ", cells " << wr.cellsBefore
                  << " -> " << wr.cellsAfter << "\n";
        if (!wr.rejections.empty()) {
            std::cout << "  Rejected:";
            for (const auto &kv : wr.rejections) std::cout << " " << kv.second << " " << kv.first << ";";
            std::cout << "\n";
        }
        if (wr.timedOut) warn("a round ran out of its time budget");
        verdict(st.rewritesValid, "The carrier validated after every round of rewrites");
    }

    // ---------------------------------------------------------------------
    heading("Best blocking (round " + std::to_string(st.bestRound) + ")");
    const SquareCarrier &C = pipe.getCarrier();
    const BlockCover::Report &br = pipe.getCover().getReport();
    {
        std::cout << "  " << br.blocks << " block(s) over " << C.numCells() << " carrier cell(s) ("
                  << br.singletonBlocks << " single-cell, largest " << br.largestBlock << ")\n";
        std::cout << "  " << br.macroVertices << " macrovertices (" << br.irregularMacroVertices
                  << " irregular), " << br.macroEdges << " macro edges\n";
        std::cout << "  Objective " << std::fixed << std::setprecision(3) << br.objective << "; distortion "
                  << br.meanDistortion << " mean, " << br.maxDistortion << " worst; side turning "
                  << br.meanComplexity << " quarter turns per block" << std::defaultfloat << "\n";
        std::cout << "  Sec. 4 on the macro complex: " << br.eulerLHS << " = " << br.eulerRHS << "\n";
        verdict(br.uncovered == 0 && br.overcovered == 0, "Every carrier cell is in exactly one block");
        verdict(br.sideMismatches == 0, "Every block side is dS or one full side of one neighbour");
        verdict(br.hangingVertices == 0, "No macrovertex inside a block side (no T-junction)");
        verdict(br.multipleContacts == 0 || !opts.cover.singleContact, "No two blocks share two sides");
        verdict(br.protectedNotMacro == 0, "Every protected vertex is a macrovertex");
        verdict(br.interfaceInside == 0, "Every interface lies on block sides");
        verdict(br.eulerHolds, "The macro complex satisfies Sec. 4's identity");
        for (const std::string &m : br.messages) warn(m);

        if (!carrierOut.empty()) {
            if (C.writeOBJ(carrierOut)) std::cout << "  Wrote the carrier to " << carrierOut << "\n";
            else warn("Failed to write " + carrierOut);
        }
        if (!blocksOut.empty()) {
            if (pipe.getCover().writeOBJ(blocksOut)) std::cout << "  Wrote the macro edges to " << blocksOut << "\n";
            else warn("Failed to write " + blocksOut);
        }
        if (!mfemOut.empty()) {
            mesh::QuadMesh qm(C.vertices, C.cells, C.cellMaterial);
            mesh::QuadMesh::MFEMReport mr;
            if (qm.writeMFEM(mfemOut, mesh::QuadMesh::MFEMOptions(), &mr)) {
                std::cout << "  Wrote the MFEM mesh to " << mfemOut << ": " << mr.elements << " element(s), "
                          << mr.boundaryElements << " boundary segment(s)\n";
            } else {
                warn("Failed to write " + mfemOut);
            }
        }
    }

    heading("Result");
    for (const std::string &m : st.messages) {
        if (m.rfind("Stopping", 0) == 0 || m.rfind("Stages 4-6", 0) == 0) warn(m);
    }
    if (ok) {
        if (br.singletonBlocks == br.blocks) {
            std::cout << "  " << kPass << " A valid conforming blocking -- the fine fallback: every carrier "
                      << "cell is its own block.\n";
        } else {
            std::cout << "  " << kPass << " A valid conforming blocking of " << br.blocks << " block(s) ("
                      << st.initialCells << " in the singleton fallback).\n";
        }
    } else {
        std::cout << "  " << kFail << " No valid conforming blocking came out.\n";
    }
    return ok ? 0 : 7;
}
