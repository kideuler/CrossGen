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
//                certificate refusing a periodic band and an L-tromino;
//                Stage 6 alone lowering the valence defect of a raw carrier
//                without ever leaving it invalid; a plate with two holes,
//                which no template covers, blocked by a coarse search and
//                realised on the input; Sec. 13.3's last row, two meshes of
//                one boundary giving blockings of about the same size; and the
//                TFI mesh on the blocks (BlockMesh), whose counts and
//                boundary are known on a rectangle and two materials.
//
// --mesh <h> then meshes the chosen blocking at target edge length h, as the
// viewer's ATLAS mode does (Sec. 11.2's counts, transfinite interpolation per
// block, Winslow where a block folds), and --mesh-tmop runs mesh::TMOP on it.
//
// A model is searched on its own carrier (Sec. 3's split of it) and on coarse
// re-triangulations of its domain (CoarseDomain), whose layouts are carried
// back onto it (Realisation) and certified there; the report shows every
// search and which one won. Exit codes: 0 valid, 2 the mesh did not load, 3
// the pipeline threw, 4 the domain failed Stage 1, 5 the carrier failed Stage
// 2, 7 no valid cover, 8 valid but coarser than --max-blocks allows, 9 the
// --mesh asked for folds (after TMOP, when --mesh-tmop asked for that too).

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "ATLAS/ATLAS.hxx"
#include "ATLAS/BlockMesh.hxx"
#include "TestHelper.hxx"
#include "mesh/QuadMesh.hxx"
#include "mesh/TMOP.hxx"

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

// The fine search: Stage 2's own carrier, the incumbent every run has.
const ATLAS::Search &fineSearch(const ATLAS &pipe) { return pipe.getSearch(0); }

void summary(const ATLAS &pipe) {
    const ATLAS::Status &st = pipe.getStatus();
    std::cout << "  " << st.triangles << " triangle(s) -> " << st.initialCells << " carrier cell(s)";
    if (pipe.hasCover()) {
        const BlockCover::Report &cr = pipe.getCover().getReport();
        std::cout << "; " << cr.blocks << " block(s) over " << pipe.getCarrier().numCells() << " cell(s), from the "
                  << pipe.getChosen().name << " search";
    }
    std::cout << "\n";
    for (int i = 0; i < pipe.numSearches(); ++i) {
        const ATLAS::Search &s = pipe.getSearch(i);
        std::cout << "     " << s.name << ": ";
        if (s.succeeded()) std::cout << s.finalCover()->getReport().blocks << " block(s)";
        else std::cout << "no valid cover";
        if (s.templates) {
            for (const ExplicitTemplates::Attempt &a : s.templates->getReport().attempts) {
                std::cout << "; region of " << a.cells << " cell(s), " << a.corners << " corner(s), " << a.holes
                          << " hole(s): " << (a.accepted ? a.family : std::string("kept"));
                if (a.accepted) std::cout << " (" << a.blocks << " block(s))";
            }
        }
        for (const std::string &m : s.messages) std::cout << "; " << m;
        std::cout << "\n";
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
        check(fineSearch(*rect).templates && fineSearch(*rect).templates->getReport().templated == 1,
              "Stage 3 replaced the region");
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
        if (fineSearch(*ring).templates && !fineSearch(*ring).templates->getReport().attempts.empty()) {
            const ExplicitTemplates::Attempt &a = fineSearch(*ring).templates->getReport().attempts.back();
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
    std::unique_ptr<ATLAS> twoMat;
    {
        auto mesh = gridMesh(20, 10, 2.0, 1.0, 1.0);
        twoMat = std::make_unique<ATLAS>(mesh, quietOptions());
        ATLAS &pipe = *twoMat;
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
        if (ring && fineSearch(*ring).templated) {
            const SquareCarrier &C = *fineSearch(*ring).templated;
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
        o.coarse = false;
        o.fineRewrite = true;
        o.rewriteRounds = 6;
        ATLAS pipe(mesh, o);
        const bool ok = pipe.run();
        const ATLAS::Search &s = fineSearch(pipe);
        for (const ATLAS::Round &r : s.history) {
            std::cout << "  round " << r.round << ": " << r.cells << " cells, defect " << r.defect
                      << ", " << r.blocks << " block(s), " << r.committed << " rewrite(s)\n";
        }
        bool monotone = true;
        for (size_t k = 1; k < s.history.size(); ++k) {
            if (s.history[k].defect > s.history[k - 1].defect) monotone = false;
        }
        check(ok, "the pipeline returns a valid blocking");
        check(s.rewritesValid, "the carrier validated after every Stage 6 round");
        check(monotone, "the valence defect never rises");
        check(s.history.size() > 1 && s.history.back().defect < s.history.front().defect,
              "Stage 6 lowered the defect");
        check(s.history.size() > 1 && pipe.getStatus().bestBlocks < s.history.front().blocks,
              "... and the best blocking is coarser than the first");
    }

    // -----------------------------------------------------------------
    heading("Case 10  A coarse search realised on the input (CoarseDomain, Realisation)");
    // -----------------------------------------------------------------
    {
        // A plate with two holes: no single template covers it (Sec. 7.4),
        // so the fine search can only aggregate the triangulation's carrier.
        // The coarse search must find a layout, carry it back onto the input's
        // own boundary, and have that pass the same validate() as Stage 2.
        auto mesh = meshLoops({{{0.0, 0.0}, {3.0, 0.0}, {3.0, 1.5}, {0.0, 1.5}},
                               circlePoints(0.8, 0.75, 0.35, 48), circlePoints(2.1, 0.7, 0.4, 48)},
                              {0, 1, 1}, 0.05);
        ATLAS pipe(mesh, quietOptions());
        const bool ok = pipe.run();
        summary(pipe);
        int realised = 0, validRealised = 0;
        for (int i = 1; i < pipe.numSearches(); ++i) {
            const ATLAS::Search &s = pipe.getSearch(i);
            if (!s.realised) continue;
            ++realised;
            if (s.realisedReport.valid) ++validRealised;
            const Realisation::Report &rr = s.realisation->getReport();
            std::cout << "  " << s.name << ": coarse " << s.coarseDomain->getReport().triangles << " triangle(s), "
                      << s.best->numCells() << " cell(s) -> realised " << rr.cells << " cell(s), "
                      << rr.inverted << " inverted, worst scaled Jacobian " << std::fixed << std::setprecision(3)
                      << rr.minScaledJacobian << std::defaultfloat << "\n";
        }
        // A coarse layout may be valid on the coarse proxy and still fold when
        // it is carried onto the true curves; that search then offers no
        // answer, and the others compete without it. What must hold is that
        // one of them does realise, and that what is chosen is certified on
        // the input by the same validate() as Stage 2's carrier.
        check(ok, "the pipeline returns a valid blocking");
        check(realised > 0 && validRealised > 0, "at least one coarse layout realises on the input");
        check(pipe.getChosen().coarse && pipe.getCarrier().getReport().valid,
              "a coarse search wins, and its carrier passes validate() on the input domain");
        const int fineFirst = fineSearch(pipe).firstCover ? fineSearch(pipe).firstCover->getReport().blocks : 0;
        if (pipe.hasCover()) {
            const BlockCover::Report &br = pipe.getCover().getReport();
            std::cout << "  " << br.blocks << " block(s); the fine carrier's first cover had " << fineFirst << "\n";
            check(br.blocks <= 60, "a coarse blocking: at most 60 blocks");
            check(br.protectedNotMacro == 0 && br.hangingVertices == 0, "every corner a macrovertex, no T-junction");
        }
    }

    // -----------------------------------------------------------------
    heading("Case 11  Same boundary, different triangulations (Sec. 13.3)");
    // -----------------------------------------------------------------
    {
        // The last regression row of Sec. 13.3: whether the blocking depends
        // on the interior triangle density. The coarse domain is sampled from
        // the boundary alone, so two meshes of one boundary differing only
        // inside should give blockings of about the same size.
        std::vector<int> blocks;
        for (double h : {0.06, 0.03}) {
            auto mesh = meshLoops({{{0.0, 0.0}, {2.0, 0.0}, {2.0, 1.0}, {0.0, 1.0}},
                                   circlePoints(0.7, 0.5, 0.25, 40), circlePoints(1.4, 0.45, 0.2, 40)},
                                  {0, 1, 1}, h);
            ATLAS pipe(mesh, quietOptions());
            pipe.run();
            const int b = pipe.hasCover() ? pipe.getCover().getReport().blocks : -1;
            std::cout << "  interior spacing " << h << ": " << mesh->triangles.size() << " triangle(s) -> " << b
                      << " block(s) from the " << (pipe.hasCover() ? pipe.getChosen().name : std::string("-"))
                      << " search\n";
            blocks.push_back(b);
        }
        check(blocks[0] > 0 && blocks[1] > 0, "both meshes give a valid blocking");
        check(blocks[0] > 0 && blocks[1] > 0 && std::abs(blocks[0] - blocks[1]) <= std::max(2, blocks[0] / 4),
              "their block counts agree to within a quarter");
    }

    // -----------------------------------------------------------------
    heading("Case 12  A TFI mesh on the blocks (Secs. 11.2, 11.3)");
    // -----------------------------------------------------------------
    {
        // Known answers first. The rectangle is one 2 x 1 block, so a target
        // of 0.1 is 20 x 10 squares whether the interior comes from the chart
        // or from the Coons blend, and a straight boundary is reproduced
        // exactly. Then the two materials: two blocks, the interface one
        // shared macro edge, so its nodes are shared and no element has two
        // materials. Then the annulus, whose curved sides the elements chord.
        if (rect && rect->hasCover()) {
            for (bool chart : {true, false}) {
                BlockMesh::Options mo;
                mo.targetEdgeLength = 0.1;
                mo.useChart = chart;
                BlockMesh bm(rect->getCover(), mo);
                const BlockMesh::Report &r = bm.getReport();
                std::cout << "  rectangle, " << (chart ? "chart" : "Coons") << ": " << r.quads << " quads over "
                          << r.blocks << " block(s), worst scaled Jacobian " << std::fixed << std::setprecision(6)
                          << r.minScaledJacobian << ", boundary deviation " << r.boundaryDeviation
                          << std::defaultfloat << "\n";
                check(r.valid && r.quads == 200 && bm.blocks().size() == 1 && bm.blocks()[0].ns * bm.blocks()[0].nt == 200,
                      std::string("rectangle: 20 x 10 elements at h = 0.1 (") + (chart ? "chart)" : "Coons)"));
                check(r.minScaledJacobian > 0.999 && r.boundaryDeviation < 1e-12,
                      "... square, and the boundary reproduced exactly");
            }
        }
        if (twoMat && twoMat->hasCover()) {
            BlockMesh::Options mo;
            mo.targetEdgeLength = 0.1;
            BlockMesh bm(twoMat->getCover(), mo);
            const BlockMesh::Report &r = bm.getReport();
            std::cout << "  two materials: " << r.quads << " quads, " << r.materials << " material(s), "
                      << r.interfaceEdges << " element edge(s) on the interface\n";
            check(r.valid && r.conforming, "two materials: a valid, conforming mesh");
            check(r.materials == 2 && r.interfaceEdges == 10 && r.interfaceDeviation < 1e-12,
                  "... two materials meeting along the interface's ten shared edges");
        }
        if (ring && ring->hasCover()) {
            const double h = 0.1;
            BlockMesh::Options mo;
            mo.targetEdgeLength = h;
            BlockMesh bm(ring->getCover(), mo);
            const BlockMesh::Report &r = bm.getReport();
            std::cout << "  annulus: " << r.quads << " quads over " << r.blocks << " block(s), " << r.chords
                      << " chord(s), worst scaled Jacobian " << std::fixed << std::setprecision(3)
                      << r.minScaledJacobian << ", boundary deviation " << std::setprecision(4)
                      << r.boundaryDeviation << ", edge ratio rms " << std::setprecision(3) << r.edgeRatioRms
                      << std::defaultfloat << "\n";
            check(r.valid && r.invertedQuads == 0, "annulus: a valid mesh, nothing folded");
            check(r.boundaryDeviation < 0.1 * h, "... chording the circles by under a tenth of the target");
            check(r.edgeRatioRms < 0.35, "... at edges near the target (rms log ratio under 0.35)");
            // The adaptor TMOP is handed the mesh through.
            mesh::QuadMesh qm = mesh::QuadMesh::from(bm);
            check(static_cast<int>(qm.quads.size()) == r.quads && qm.nonManifoldEdges.empty() &&
                      qm.boundaryLoops.size() == 2,
                  "mesh::QuadMesh::from takes it whole: every element, two boundary loops, manifold");
        }
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

// Stages 3-6 of one search, and for a coarse one its domain and realisation.
void printSearch(const ATLAS::Search &s, bool chosen) {
    heading("Search: " + s.name + (chosen ? "  (chosen)" : ""));
    if (s.coarse) {
        if (!s.coarseDomain) return;
        const CoarseDomain::Report &cr = s.coarseDomain->getReport();
        std::cout << "  CoarseDomain: " << cr.samples << " boundary sample(s) on " << cr.chains << " chain(s), "
                  << cr.triangles << " triangle(s), spacing " << std::setprecision(3) << cr.spacingMin << " to "
                  << cr.spacingMax << std::setprecision(6) << ", " << std::fixed << std::setprecision(3)
                  << cr.seconds << " s" << std::defaultfloat << "\n";
        if (!cr.valid) {
            warn("no coarse domain: " + cr.reason);
            return;
        }
        std::cout << "  Coarse carrier: " << s.initialReport.cells << " cells, " << s.initialReport.irregularInterior
                  << " + " << s.initialReport.irregularBoundary << " irregular, defect " << s.initialReport.totalDefect
                  << "\n";
        verdict(s.carrierValid, "The coarse carrier validates (Sec. 3 on the coarse triangulation)");
    }
    if (s.templates) {
        const ExplicitTemplates::Report &tr = s.templates->getReport();
        std::cout << "  Stage 3: " << tr.templated << " of " << tr.regions << " region(s) replaced, " << std::fixed
                  << std::setprecision(1) << 100.0 * tr.templatedArea << "% of the area" << std::defaultfloat
                  << ", " << tr.blocks << " block(s) planned; cells " << tr.cellsBefore << " -> " << tr.cellsAfter
                  << "\n";
        int shown = 0;
        for (const ExplicitTemplates::Attempt &a : tr.attempts) {
            if (shown++ >= 12) { std::cout << "     ...\n"; break; }
            std::cout << "     mat " << a.material << ", " << a.cells << " cell(s), " << a.corners << " corner(s), "
                      << a.holes << " hole(s): ";
            if (a.accepted) {
                std::cout << a.family << ", " << a.blocks << " block(s), " << a.maps << "; worst scaled Jacobian "
                          << std::fixed << std::setprecision(3) << a.minScaledJacobian << std::defaultfloat;
            } else {
                std::cout << "kept (" << a.reason << ")";
            }
            std::cout << "\n";
        }
        verdict(s.templatesValid, "The carrier still validates after every replacement");
    }
    if (!s.history.empty()) {
        std::cout << "  " << std::setw(6) << "round" << std::setw(9) << "cells" << std::setw(9) << "defect"
                  << std::setw(11) << "irregular" << std::setw(12) << "candidates" << std::setw(9) << "blocks"
                  << std::setw(10) << "rewrites" << std::setw(7) << "rings" << std::setw(9) << "sec" << "\n";
        const size_t n = s.history.size();
        for (size_t k = 0; k < n; ++k) {
            // Long searches: the first rounds, the best one and the last.
            const ATLAS::Round &r = s.history[k];
            if (n > 16 && k >= 6 && k + 4 < n && r.round != s.bestRound) {
                if (k == 6) std::cout << "  " << std::setw(6) << "..." << "\n";
                continue;
            }
            std::cout << "  " << std::setw(6) << r.round << std::setw(9) << r.cells << std::setw(9) << r.defect
                      << std::setw(11) << r.irregular << std::setw(12) << r.candidates << std::setw(9) << r.blocks
                      << std::setw(10) << r.committed << std::setw(7) << r.rings << std::setw(9) << std::fixed
                      << std::setprecision(2) << r.seconds << std::defaultfloat
                      << (r.round == s.bestRound ? "  <- best" : "") << "\n";
        }
    }
    if (s.rewrite) {
        const CavityRewrite::Report &wr = s.rewrite->getReport();
        std::cout << "  Stage 6: " << wr.committed << " rewrite(s) committed of " << wr.proposals << " proposal(s) ("
                  << wr.gridFills << " grid, " << wr.starFills << " star); defect " << wr.defectBefore << " -> "
                  << wr.defectAfter << ", irregular " << wr.irregularBefore << " -> " << wr.irregularAfter
                  << ", cells " << wr.cellsBefore << " -> " << wr.cellsAfter << "\n";
        if (!wr.rejections.empty()) {
            std::cout << "  Rejected:";
            for (const auto &kv : wr.rejections) std::cout << " " << kv.second << " " << kv.first << ";";
            std::cout << "\n";
        }
        if (s.annealed) {
            const CavityRewrite::AnnealReport &ar = s.rewrite->getAnnealReport();
            std::cout << "  Annealing (Sec. 8.3's energy, N_B included): " << ar.moves << " move(s), " << ar.certified
                      << " certified, " << ar.accepted << " accepted (" << ar.uphill << " uphill), "
                      << ar.improvements << " new best; energy " << std::fixed << std::setprecision(2)
                      << ar.energyBefore << " -> " << ar.energyAfter << ", base patches " << ar.blocksBefore
                      << " -> " << ar.blocksAfter << "; " << ar.seconds << " s" << std::defaultfloat
                      << (ar.timedOut ? " (hit the time cap)" : "") << "\n";
        }
        verdict(s.rewritesValid, "The carrier validated after every round of rewrites");
    }
    if (s.firstCover && s.firstRects) {
        const RectangleCertifier::Report &rr = s.firstRects->getReport();
        std::cout << "  First cover: " << s.firstCover->getReport().blocks << " block(s) (" << rr.forcedVertices
                  << " forced macrovertices, " << rr.candidates << " candidates)";
        if (s.bestCover) {
            std::cout << "; best: " << s.bestCover->getReport().blocks << " block(s) over " << s.best->numCells()
                      << " cells, round " << s.bestRound;
        }
        std::cout << "\n";
    }
    if (s.coarse && s.realisation) {
        const Realisation::Report &rr = s.realisation->getReport();
        std::cout << "  Realisation on the input: " << rr.classes << " count class(es), " << rr.cells << " cell(s), "
                  << rr.snapped << " snapped, " << rr.boundarySplits << " boundary point(s) inserted\n";
        std::cout << "  Interpolated: " << rr.invertedInterpolated << " inverted, worst scaled Jacobian " << std::fixed
                  << std::setprecision(3) << rr.minScaledJacobianInterpolated << "; "
                  << (rr.invertedHarmonic >= 0 ? "harmonic map " + std::to_string(rr.invertedHarmonic) + " inverted" +
                                                     (rr.harmonic ? " (kept); " : "; ")
                                               : std::string())
                  << (rr.smoothed ? "after TMOP" : "TMOP not kept") << ": " << rr.inverted << " inverted, worst "
                  << rr.minScaledJacobian << ", mean " << rr.meanScaledJacobian << "; " << rr.seconds << " s"
                  << std::defaultfloat << "\n";
        if (s.realised) {
            verdict(s.realisedReport.valid, "The realised carrier passes Stage 2's validate() on the input domain");
        }
        if (s.realisedCover) {
            const BlockCover::Report &cr = s.realisedCover->getReport();
            std::cout << "  On the input: " << cr.blocks << " block(s) over " << s.realised->numCells() << " cells, "
                      << cr.irregularMacroVertices << " irregular macrovertices\n";
            verdict(cr.valid, "Sec. 13.2: the realised cover is conforming");
        }
    }
    for (const std::string &m : s.messages) warn(m);
    std::cout << "  " << std::fixed << std::setprecision(2) << s.seconds << " s" << std::defaultfloat << "\n";
}

// The TFI mesh on the chosen blocks, in the words the viewer's Mesh phase uses.
void printMesh(const BlockMesh &bm) {
    const BlockMesh::Report &r = bm.getReport();
    std::cout << "  " << r.quads << " quad(s) on " << r.vertices << " vertices over " << r.blocks << " block(s); "
              << r.chords << " chord(s) over " << r.arcsAssigned << " macro edge(s), " << r.minIntervals << " to "
              << r.maxIntervals << " edges each (mean " << std::fixed << std::setprecision(2) << r.meanIntervals
              << ")" << std::defaultfloat;
    if (r.clampedChords > 0) std::cout << ", " << r.clampedChords << " clamped by a bound";
    std::cout << "\n";
    std::cout << "  Edges " << std::fixed << std::setprecision(4) << r.minEdge << " to " << r.maxEdge
              << " against a target of " << r.target << " (worst " << std::setprecision(2) << r.worstEdgeRatio
              << "x, rms log ratio " << std::setprecision(3) << r.edgeRatioRms << ")\n";
    std::cout << "  Scaled Jacobian " << std::setprecision(4) << r.minScaledJacobian << " worst, "
              << r.meanScaledJacobian << " mean; before the Winslow pass " << r.minScaledJacobianBefore << " worst, "
              << r.invertedBefore << " folded (" << r.smoothedBlocks << " block(s) smoothed)\n";
    std::cout << "  Boundary: chords within " << std::setprecision(5) << r.boundaryDeviation << " of dS ("
              << std::setprecision(3) << (r.target > 0.0 ? r.boundaryDeviation / r.target : 0.0)
              << " of the target)";
    if (r.materials > 1) {
        std::cout << "; interfaces within " << std::setprecision(5) << r.interfaceDeviation << ", "
                  << r.interfaceEdges << " element edge(s) on them";
    }
    std::cout << "; area " << std::setprecision(6) << r.meshArea << " of the domain's " << r.domainArea
              << std::defaultfloat << "\n";
    verdict(r.conforming, "Conforming: no edge used three times, no crack (Sec. 11.3's incidence)");
    verdict(r.invertedQuads == 0, "No element folds (Sec. 11.3's four corner Jacobians)");
    verdict(r.unmeshedBlocks == 0, "Every block meshed");
    for (const std::string &m : r.messages) warn(m);
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
              << "  --time <s>            budget for one search's Stage 4-6 loop        (default 300)\n"
              << "  --fine-rewrite        run Stage 6 on the fine carrier even when a coarse search succeeds\n\n"
              << "Coarse searches (CoarseDomain, Realisation)\n"
              << "  --no-coarse           search on the fine carrier only\n"
              << "  --coarse-spacing <h>  a coarse search at this fraction of the diagonal; repeatable\n"
              << "                        (default: 0.125 and 0.07)\n"
              << "  --coarse-rounds <n>   rounds of Stages 4-6 on a coarse carrier       (default 30)\n"
              << "  --anneal <n>          moves of annealed search per coarse carrier, 0 = off (default 15000)\n"
              << "  --seed <n>            the annealer's random seed\n"
              << "  --realise-size <h>    fine edge length of the realisation (default: the input's mean)\n"
              << "  --realise-sweeps <n>  TMOP sweeps on the realisation, 0 = none       (default 100)\n\n"
              << "Output\n"
              << "  --initial <f.obj>     the Stage 2 carrier\n"
              << "  --templated <f.obj>   the chosen search's carrier after Stage 3\n"
              << "  --coarse <f.obj>      the chosen coarse search's best coarse carrier\n"
              << "  --coarse-blocks <f.obj> its macro edges\n"
              << "  --dump-dir <dir>      every search's coarse carrier, coarse blocks and realisation\n"
              << "  --max-blocks <n>      fail (exit 8) unless the blocking has at most n blocks\n"
              << "  --cover-model <f>     Sec. 6's exact-cover program for the layout the answer came from\n"
              << "                        (the coarse carrier's for a coarse search), for an external solver\n"
              << "  --carrier <f.obj>     the carrier the best blocking lives on\n"
              << "  --blocks <f.obj>      the macro edges of the best blocking, as polylines\n"
              << "  --mfem <f.mesh>       the carrier of the best blocking for MFEM, material per element\n\n"
              << "A quadrilateral mesh on the blocks (BlockMesh; Secs. 11.2, 11.3)\n"
              << "  --mesh <h>            mesh the best blocking at target edge length h, transfinite\n"
              << "                        interpolation per block; exit 9 unless it validates\n"
              << "  --mesh-coons          interiors by the Coons blend of the sides, not through the chart\n"
              << "  --mesh-min <n>        fewest edges per chord                         (default 1)\n"
              << "  --mesh-max <n>        most edges per chord, 0 = none                 (default 0)\n"
              << "  --mesh-no-smooth      no Winslow pass on the folded blocks\n"
              << "  --mesh-tmop <n>       then TMOP (mesh::TMOP's defaults, sampled at the corners) for at\n"
              << "                        most n sweeps\n"
              << "  --mesh-obj <f.obj>    write the mesh (after TMOP when that ran)\n";
}

} // namespace

int main(int argc, char **argv) {
    if (argc < 2) { usage(argv[0]); return 1; }
    const std::string first = argv[1];
    if (first == "--selftest") return selfTest();
    if (first == "-h" || first == "--help") { usage(argv[0]); return 0; }

    const std::string path = first;
    ATLAS::Options opts;
    std::string initialOut, templatedOut, carrierOut, blocksOut, mfemOut, coarseOut, coarseBlocksOut, dumpDir;
    std::string modelOut;
    int maxBlocks = -1;
    std::vector<double> spacings;
    BlockMesh::Options meshOpts;
    meshOpts.targetEdgeLength = 0.0;   // no mesh unless --mesh
    int meshTMOP = 0;
    std::string meshOut;

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
        else if (a == "--fine-rewrite")                   opts.fineRewrite = true;
        else if (a == "--no-coarse")                      opts.coarse = false;
        else if (a == "--coarse-spacing" && i + 1 < argc) spacings.push_back(std::stod(argv[++i]));
        else if (a == "--coarse-rounds" && i + 1 < argc)  opts.coarseRewriteRounds = std::stoi(argv[++i]);
        else if (a == "--anneal" && i + 1 < argc)         opts.anneal.maxMoves = std::stoi(argv[++i]);
        else if (a == "--seed" && i + 1 < argc)           opts.anneal.seed = static_cast<unsigned>(std::stoul(argv[++i]));
        else if (a == "--realise-size" && i + 1 < argc)   opts.realisation.size = std::stod(argv[++i]);
        else if (a == "--realise-sweeps" && i + 1 < argc) opts.realisation.smoothingSweeps = std::stoi(argv[++i]);
        else if (a == "--coarse" && i + 1 < argc)         coarseOut = argv[++i];
        else if (a == "--coarse-blocks" && i + 1 < argc)  coarseBlocksOut = argv[++i];
        else if (a == "--dump-dir" && i + 1 < argc)       dumpDir = argv[++i];
        else if (a == "--cover-model" && i + 1 < argc)    modelOut = argv[++i];
        else if (a == "--max-blocks" && i + 1 < argc)     maxBlocks = std::stoi(argv[++i]);
        else if (a == "--initial" && i + 1 < argc)        initialOut = argv[++i];
        else if (a == "--templated" && i + 1 < argc)      templatedOut = argv[++i];
        else if (a == "--carrier" && i + 1 < argc)        carrierOut = argv[++i];
        else if (a == "--blocks" && i + 1 < argc)         blocksOut = argv[++i];
        else if (a == "--mfem" && i + 1 < argc)           mfemOut = argv[++i];
        else if (a == "--mesh" && i + 1 < argc)           meshOpts.targetEdgeLength = std::stod(argv[++i]);
        else if (a == "--mesh-coons")                     meshOpts.useChart = false;
        else if (a == "--mesh-min" && i + 1 < argc)       meshOpts.minIntervals = std::stoi(argv[++i]);
        else if (a == "--mesh-max" && i + 1 < argc)       meshOpts.maxIntervals = std::stoi(argv[++i]);
        else if (a == "--mesh-no-smooth")                 meshOpts.smoothingPasses = 0;
        else if (a == "--mesh-tmop" && i + 1 < argc)      meshTMOP = std::stoi(argv[++i]);
        else if (a == "--mesh-obj" && i + 1 < argc)       meshOut = argv[++i];
        else { std::cerr << "Unknown option: " << a << "\n"; usage(argv[0]); return 1; }
    }

    if (!spacings.empty()) opts.coarseSpacings = spacings;

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
    for (int i = 0; i < pipe.numSearches(); ++i) printSearch(pipe.getSearch(i), i == st.chosen);
    if (!dumpDir.empty()) {
        for (int i = 0; i < pipe.numSearches(); ++i) {
            const ATLAS::Search &s = pipe.getSearch(i);
            const std::string base = dumpDir + "/search" + std::to_string(i);
            if (s.initial) s.initial->writeOBJ(base + "_initial.obj");
            if (s.best) s.best->writeOBJ(base + "_best.obj");
            if (s.bestCover) s.bestCover->writeOBJ(base + "_bestblocks.obj");
            if (s.realised) s.realised->writeOBJ(base + "_realised.obj");
            if (s.realisedCover) s.realisedCover->writeOBJ(base + "_realisedblocks.obj");
        }
        std::cout << "\n  Dumped every search to " << dumpDir << "\n";
    }

    if (!pipe.hasCover()) {
        heading("Result");
        std::cout << "  " << kFail << " No cover was selected.\n";
        return 7;
    }

    // ---------------------------------------------------------------------
    const ATLAS::Search &chosen = pipe.getChosen();
    heading("Best blocking: the " + chosen.name + " search" +
            (chosen.coarse ? std::string(", realised on the input") : ", round " + std::to_string(chosen.bestRound)));
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
        std::cout << "  Carrier scaled Jacobian " << std::fixed << std::setprecision(4)
                  << C.getReport().minScaledJacobian << " worst, " << C.getReport().meanScaledJacobian << " mean"
                  << std::defaultfloat << "\n";
        std::cout << "  Sec. 4 on the macro complex: " << br.eulerLHS << " = " << br.eulerRHS << "\n";
        verdict(C.getReport().valid, "The carrier the blocks live on passes Stage 2's validate() on the input");
        verdict(br.uncovered == 0 && br.overcovered == 0, "Every carrier cell is in exactly one block");
        verdict(br.sideMismatches == 0, "Every block side is dS or one full side of one neighbour");
        verdict(br.hangingVertices == 0, "No macrovertex inside a block side (no T-junction)");
        verdict(br.multipleContacts == 0 || !opts.cover.singleContact, "No two blocks share two sides");
        verdict(br.protectedNotMacro == 0, "Every protected vertex is a macrovertex");
        verdict(br.interfaceInside == 0, "Every interface lies on block sides");
        verdict(br.eulerHolds, "The macro complex satisfies Sec. 4's identity");
        for (const std::string &m : br.messages) warn(m);

        if (!templatedOut.empty() && chosen.templated) {
            if (chosen.templated->writeOBJ(templatedOut)) std::cout << "  Wrote the templated carrier to " << templatedOut << "\n";
            else warn("Failed to write " + templatedOut);
        }
        if (!coarseOut.empty()) {
            if (chosen.coarse && chosen.best && chosen.best->writeOBJ(coarseOut)) {
                std::cout << "  Wrote the coarse carrier to " << coarseOut << "\n";
            } else {
                warn("No coarse carrier to write to " + coarseOut);
            }
        }
        if (!coarseBlocksOut.empty()) {
            if (chosen.coarse && chosen.bestCover && chosen.bestCover->writeOBJ(coarseBlocksOut)) {
                std::cout << "  Wrote the coarse macro edges to " << coarseBlocksOut << "\n";
            } else {
                warn("No coarse cover to write to " + coarseBlocksOut);
            }
        }
        if (!modelOut.empty()) {
            const BlockCover *mc = chosen.bestCover.get();
            if (mc && mc->writeModel(modelOut)) std::cout << "  Wrote Sec. 6's exact-cover program to " << modelOut << "\n";
            else warn("Failed to write " + modelOut);
        }
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

    // ---------------------------------------------------------------------
    bool meshValid = true;
    if (meshOpts.targetEdgeLength > 0.0) {
        heading("A quadrilateral mesh on the blocks (Secs. 11.2, 11.3)");
        BlockMesh bm(pipe.getCover(), meshOpts);
        printMesh(bm);
        meshValid = bm.getReport().valid;
        mesh::QuadMesh qm = mesh::QuadMesh::from(bm);
        if (meshTMOP > 0) {
            // Sampled at the corners, as the viewer's ATLAS mode does: the
            // barrier then guards the corner Jacobians Sec. 9.1 judges an
            // element by. At the 2x2 Gauss points a corner can turn over
            // unseen, and on six corpus models it did.
            mesh::TMOP::Options to;
            to.maxSweeps = meshTMOP;
            to.quadrature = mesh::TMOP::Corners;
            mesh::TMOP smoother(qm, to);
            const bool improved = smoother.run();
            const mesh::TMOP::Report &tr = smoother.getReport();
            std::cout << "  TMOP: " << tr.sweeps << " sweep(s)";
            if (tr.untangleSweeps > 0) std::cout << " after " << tr.untangleSweeps << " untangling";
            std::cout << "; scaled Jacobian " << std::fixed << std::setprecision(4) << tr.minScaledJacobianBefore
                      << " -> " << tr.minScaledJacobianAfter << " worst, " << tr.meanScaledJacobianBefore << " -> "
                      << tr.meanScaledJacobianAfter << " mean, folds " << tr.invertedBefore << " -> "
                      << tr.invertedAfter << std::defaultfloat << "\n";
            verdict(improved && tr.invertedAfter == 0, "TMOP left the mesh better than it found it, nothing folded");
            // The mesh handed on is the smoothed one: it is what must not fold.
            meshValid = tr.invertedAfter == 0;
        }
        if (!meshOut.empty()) {
            if (qm.writeOBJ(meshOut)) std::cout << "  Wrote the mesh to " << meshOut << "\n";
            else warn("Failed to write " + meshOut);
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
            std::cout << "  " << kPass << " A valid conforming blocking of " << br.blocks << " block(s) from the "
                      << chosen.name << " search (" << st.initialCells << " in the singleton fallback).\n";
        }
    } else {
        std::cout << "  " << kFail << " No valid conforming blocking came out.\n";
    }
    if (ok && maxBlocks >= 0) {
        const bool coarse = br.blocks <= maxBlocks;
        verdict(coarse, "At most " + std::to_string(maxBlocks) + " blocks (--max-blocks)");
        if (!coarse) return 8;
    }
    if (ok && !meshValid) return 9;
    return ok ? 0 : 7;
}
