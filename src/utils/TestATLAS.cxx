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
// --mesh asked for folds (after TMOP, when --mesh-tmop asked for that too), 10
// that mesh's worst element is below --min-sj.
//
// The searches are scored against a DualMBO reference cross field by default
// (docs/atlas_crossfield_guidance.md); --field harmonic --no-signed-defects
// --no-patch-topology is ATLAS as it was before it.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
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

// Print the reference field's diagnostics (--field-report, or a DualMBO run).
bool showField = false;

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

    // -----------------------------------------------------------------
    heading("Case 13  Half disk: a four-block half O-grid");
    // -----------------------------------------------------------------
    {
        // Every body that crosses the axis of an (r, z) model is one of these.
        // Its two corners are real quarter turns, and Sec. 7.3's O-grid would
        // make them split points and cut each into two cells; half of an
        // O-grid gives each one cell and needs two three-valent vertices
        // inside instead of four.
        std::vector<Point> loop;
        const int n = 36;
        for (int i = 0; i <= n; ++i) loop.push_back({std::cos(M_PI * i / n), std::sin(M_PI * i / n)});
        auto mesh = meshLoops({loop}, {0}, 0.08);
        ATLAS pipe(mesh, quietOptions());
        const bool ok = pipe.run();
        summary(pipe);
        check(ok && pipe.hasCover(), "the pipeline returns a valid blocking");
        if (fineSearch(pipe).templates && !fineSearch(pipe).templates->getReport().attempts.empty()) {
            const ExplicitTemplates::Attempt &a = fineSearch(pipe).templates->getReport().attempts.back();
            check(a.accepted && a.family == "half O-grid" && a.blocks == 4,
                  "Stage 3 lays a core on the diameter and three shells round it");
        }
        if (ok && pipe.hasCover()) {
            const BlockCover &cover = pipe.getCover();
            const SquareCarrier &C = pipe.getCarrier();
            check(cover.getReport().blocks == 4, "four blocks");
            // Block corners at each carrier vertex: one at each end of the
            // diameter, three at the core's two top corners, none irregular
            // anywhere else.
            std::vector<int> corners(C.vertices.size(), 0);
            for (const BlockCover::Block &b : cover.getBlocks()) for (int v : b.cert.corners) ++corners[v];
            int endsOneBlock = 0, threeValent = 0, otherIrregular = 0;
            for (size_t v = 0; v < C.vertices.size(); ++v) {
                if (corners[v] == 0) continue;
                const bool end = std::fabs(std::fabs(C.vertices[v][0]) - 1.0) < 1e-9 && std::fabs(C.vertices[v][1]) < 1e-9;
                if (end && corners[v] == 1) ++endsOneBlock;
                else if (!C.boundaryVertex[v] && corners[v] == 3) ++threeValent;
                else if (C.boundaryVertex[v] ? corners[v] != 2 : corners[v] != 4) ++otherIrregular;
            }
            check(endsOneBlock == 2, "each end of the diameter is the corner of one block");
            check(threeValent == 2 && otherIrregular == 0,
                  "two three-valent macrovertices, the core's top corners, and nothing else irregular");
            BlockMesh::Options mo;
            mo.targetEdgeLength = 0.1;
            BlockMesh bm(cover, mo);
            std::cout << "  TFI mesh at h = 0.1: worst scaled Jacobian " << std::fixed << std::setprecision(3)
                      << bm.getReport().minScaledJacobian << std::defaultfloat << "\n";
            // 1/sqrt(2) is the O-grid's own figure: at a core corner the core
            // takes a quarter turn and the two shells 135 degrees each.
            check(bm.getReport().valid && bm.getReport().minScaledJacobian > 0.65,
                  "the TFI mesh on it is valid, no element below 0.65");
        }
    }

    // -----------------------------------------------------------------
    heading("Case 14  The reference cross field (docs/atlas_crossfield_guidance.md)");
    // -----------------------------------------------------------------
    {
        ATLAS::Options fo = quietOptions();
        fo.field.reference = ATLAS::Options::Field::Reference::DualMBO;
        const ReferenceField::SingularityWeights &w = fo.field.singularity;

        // A rectangle turned 30 degrees: the field is its axes everywhere, so
        // the angle conventions -- the field's, the cell's own cross, the
        // average -- agree end to end or this fails.
        const double a = M_PI / 6.0, ca = std::cos(a), sa = std::sin(a);
        std::vector<Point> box;
        for (const Point &p : std::vector<Point>{{0.0, 0.0}, {2.0, 0.0}, {2.0, 1.0}, {0.0, 1.0}}) {
            box.push_back({ca * p[0] - sa * p[1], sa * p[0] + ca * p[1]});
        }
        auto rmesh = meshLoops({box}, {0}, 0.08);
        ReferenceField F(rmesh, {}, ReferenceField::Options());
        double rho = 0.0;
        const double th = F.crossAt({ca * 1.0 - sa * 0.5, sa * 1.0 + ca * 0.5}, &rho);
        std::cout << "  rotated rectangle: " << F.cones().size() << " cone(s), cross at the centre "
                  << std::setprecision(4) << th * 180.0 / M_PI << " deg, coherence " << rho << std::setprecision(6)
                  << ", " << F.getReport().levels << " tau level(s), " << F.getReport().seconds << " s\n";
        check(F.built() && F.cones().empty(), "the rectangle's field has no cones");
        check(std::fabs(th - a) < M_PI / 180.0 && rho > 0.99, "its cross at the centre is the rectangle's, 30 degrees");
        {
            ATLAS pipe(rmesh, fo);
            const bool ok = pipe.run();
            summary(pipe);
            check(ok && pipe.hasCover() && pipe.getCover().getReport().blocks == 1, "one block with the field on");
            if (pipe.hasCover()) {
                const double e = F.directionEnergy(pipe.getCarrier());
                std::cout << "  E_dir of the block " << e << "\n";
                check(e < 1e-3, "E_dir of the block is below 1e-3: its cells follow the field");
            }
        }

        // E_sing's matching, on hand-placed units (r0 = 0.1 of the diagonal).
        {
            const double r0 = w.r0 * F.diagonal();
            auto unit = [](double x, double y, int sign, bool boundary) {
                ReferenceField::Singularity s;
                s.x = {x, y};
                s.sign = sign;
                s.boundary = boundary;
                return s;
            };
            const std::vector<ReferenceField::Singularity> none, conePlus{unit(0.0, 0.0, 1, false)};
            const double e1 = F.singularityEnergy({unit(0.2 * r0, 0.0, 1, false)}, conePlus, w);
            const double e2 = F.singularityEnergy({unit(0.0, 0.0, 1, true)}, conePlus, w);
            const double e3 = F.singularityEnergy(none, conePlus, w);
            const double e4 = F.singularityEnergy({unit(0.0, 0.0, 1, false), unit(0.3 * r0, 0.0, -1, false)}, none, w);
            const double e5 = F.singularityEnergy({unit(0.0, 0.0, -1, false)}, conePlus, w);
            const double e6 = F.singularityEnergy({unit(0.0, 0.0, -1, true)}, none, w);
            std::cout << "  E_sing: near " << e1 << ", pushed to dS " << e2 << ", missing " << e3 << ", dipole " << e4
                      << ", wrong sign " << e5 << ", -1/4 on dS " << e6 << "\n";
            check(std::fabs(e1 - 0.2 * w.wPos) < 1e-12, "a singularity 0.2 r0 from its cone costs 0.2 wPos");
            check(std::fabs(e2 - w.wEdgePlus) < 1e-12, "a cone pushed onto dS costs wEdge+ and is not also missing");
            check(std::fabs(e3 - w.wMissing) < 1e-12, "a cone with nothing at it costs wMissing");
            check(std::fabs(e4 - 0.3 * w.wPos) < 1e-12, "a carrier dipole 0.3 r0 wide cancels at 0.3 wPos");
            check(std::fabs(e5 - w.wExtra - w.wMissing) < 1e-12, "opposite signs never match");
            check(std::fabs(e6 - w.wEdgeMinus) < 1e-12, "a -1/4 on dS costs wEdge-");
        }

        // A disk: the four cones on the diagonals (DualMBO pins the centre),
        // the O-grid laid on them, and Sec. 4.3's invariance.
        {
            auto dmesh = TestHelper::createCircle(0.0, 0.0, 1.0, 0.08);
            ATLAS pipe(dmesh, fo);
            const bool ok = pipe.run();
            summary(pipe);
            const ReferenceField *D = pipe.getField();
            check(D != nullptr, "ATLAS built the field after Stage 1");
            if (D) {
                int plus = 0, minus = 0;
                double worst = 0.0, radius = 0.0;
                for (const ReferenceField::Singularity &s : D->cones()) {
                    (s.sign > 0 ? plus : minus)++;
                    const double ang = std::atan2(s.x[1], s.x[0]);
                    const double k = std::round((ang - M_PI_4) / M_PI_2);
                    worst = std::max(worst, std::fabs(ang - M_PI_4 - k * M_PI_2));
                    radius += normP(s.x) / 4.0;
                }
                std::cout << "  disk: cones +" << plus << "/-" << minus << ", worst " << worst * 180.0 / M_PI
                          << " deg off a diagonal, mean radius " << radius << "\n";
                check(plus == 4 && minus == 0, "four +1/4 cones");
                check(worst < 5.0 * M_PI / 180.0, "each within 5 degrees of a diagonal");
            }
            check(ok && pipe.hasCover() && pipe.getCover().getReport().blocks == 5, "a five-block O-grid");
            if (fineSearch(pipe).templates && !fineSearch(pipe).templates->getReport().attempts.empty()) {
                const ExplicitTemplates::Attempt &at = fineSearch(pipe).templates->getReport().attempts.back();
                check(at.accepted && at.family == "O-grid" && at.conePlaced, "Stage 3 placed the O-grid on the cones");
                check(at.singPlus == 4 && at.edgePlus == 0 && at.edgeMinus == 0,
                      "its four three-valent vertices are the singular ones, nothing on dS");
            }
            if (D && ok && pipe.hasCover()) {
                const SquareCarrier &C = pipe.getCarrier();
                ReferenceField::SingularityReport sr;
                const double es = D->singularityEnergy(C, fo.field.singularity, &sr);
                std::cout << "  E_sing " << es << ": " << sr.matched << " matched, mean " << sr.meanDistance << " r0\n";
                check(sr.matched == 4 && sr.missing == 0 && sr.extra == 0 && es < 2.0 * w.wPos * 0.25,
                      "every cone matched, each within a quarter of r0");
                // The split rays on the diagonals: the four block corners on the circle.
                double off = 0.0;
                int onCircle = 0;
                for (int v : pipe.getCover().getMacroVertices()) {
                    if (!C.boundaryVertex[v]) continue;
                    ++onCircle;
                    const double ang = std::atan2(C.vertices[v][1], C.vertices[v][0]);
                    const double k = std::round((ang - M_PI_4) / M_PI_2);
                    off = std::max(off, std::fabs(ang - M_PI_4 - k * M_PI_2));
                }
                std::cout << "  split rays: " << onCircle << " on the circle, worst " << off * 180.0 / M_PI
                          << " deg off a diagonal\n";
                check(onCircle == 4 && off < 8.0 * M_PI / 180.0, "the split rays run on the diagonals");
                // Sec. 4.3: sum over blocks of their cells' misalignment is the
                // carrier's, whatever the cover -- which is why it is not in cost(P).
                double byBlock = 0.0;
                for (const BlockCover::Block &b : pipe.getCover().getBlocks()) {
                    for (int q : b.cert.cells) byBlock += D->misalignment(C, q);
                }
                const double whole = D->directionEnergy(C) * D->area();
                std::cout << "  misalignment by block " << byBlock << ", of the carrier " << whole << "\n";
                check(std::fabs(byBlock - whole) <= 1e-9 * std::max(1.0, whole),
                      "cover invariance: the blocks' sum is the carrier's (Sec. 4.3)");
            }
        }

        // The half disk again: with the field, Stage 3 builds both the half
        // O-grid and the O-grid and keeps the half O-grid on score -- the
        // dispatch rule of Case 13, reached without being written down.
        {
            std::vector<Point> loop;
            const int n = 36;
            for (int i = 0; i <= n; ++i) loop.push_back({std::cos(M_PI * i / n), std::sin(M_PI * i / n)});
            // At h = 0.1 Sec. 7.3's O-grid certifies on it too (at 0.08 one
            // of its rays misses the core), so there is a choice to make.
            auto hmesh = meshLoops({loop}, {0}, 0.1);
            ATLAS pipe(hmesh, fo);
            const bool ok = pipe.run();
            summary(pipe);
            check(ok && pipe.hasCover(), "the pipeline returns a valid blocking");
            if (fineSearch(pipe).templates && !fineSearch(pipe).templates->getReport().attempts.empty()) {
                const ExplicitTemplates::Attempt &at = fineSearch(pipe).templates->getReport().attempts.back();
                std::cout << "  " << at.family << ": score " << at.score << " (E_sing " << at.eSing << ", cones +"
                          << at.conesPlus << "/-" << at.conesMinus << "); also built: " << at.alternatives << "\n";
                check(at.accepted && at.family == "half O-grid" && at.alternatives.find("O-grid") != std::string::npos,
                      "both families built, the half O-grid kept on score");
            }
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
            if (a.scored) {
                std::cout << "       field: cones +" << a.conesPlus << "/-" << a.conesMinus << "; template +"
                          << a.singPlus << "/-" << a.singMinus << " inside, +" << a.edgePlus << "/-" << a.edgeMinus
                          << " on dS; E_dir " << std::fixed << std::setprecision(4) << a.eDir << ", E_sing "
                          << std::setprecision(2) << a.eSing << ", score " << a.score << std::defaultfloat
                          << (a.conePlaced ? "; placed on the cones" : "");
                if (!a.alternatives.empty()) std::cout << "; also built: " << a.alternatives;
                std::cout << "\n";
            }
        }
        verdict(s.templatesValid, "The carrier still validates after every replacement");
    }
    if (!s.history.empty()) {
        std::cout << "  " << std::setw(6) << "round" << std::setw(9) << "cells" << std::setw(9) << "defect"
                  << std::setw(11) << "irregular" << std::setw(12) << "candidates" << std::setw(9) << "blocks"
                  << std::setw(10) << "rewrites" << std::setw(7) << "rings" << std::setw(9) << "sec";
        if (showField) std::cout << std::setw(9) << "score" << std::setw(8) << "E_dir" << std::setw(8) << "E_sing";
        std::cout << "\n";
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
                      << std::setprecision(2) << r.seconds << std::defaultfloat;
            if (showField) {
                std::cout << std::fixed << std::setprecision(2) << std::setw(9) << r.score << std::setprecision(3)
                          << std::setw(8) << r.eDir << std::setprecision(2) << std::setw(8) << r.eSing
                          << std::defaultfloat;
            }
            std::cout << (r.round == s.bestRound ? "  <- best" : "") << "\n";
        }
    }
    if (s.rewrite) {
        const CavityRewrite::Report &wr = s.rewrite->getReport();
        std::cout << "  Stage 6: " << wr.committed << " rewrite(s) committed of " << wr.proposals << " proposal(s) ("
                  << wr.gridFills << " grid, " << wr.starFills << " star); defect " << wr.defectBefore << " -> "
                  << wr.defectAfter << ", irregular " << wr.irregularBefore << " -> " << wr.irregularAfter
                  << ", cells " << wr.cellsBefore << " -> " << wr.cellsAfter << "\n";
        if (wr.fieldChoices > 0 || wr.conePlacedStars > 0) {
            std::cout << "  Field in the greedy rounds: " << wr.fieldChoices << " fill(s) chosen over the first "
                      << "certified, " << wr.conePlacedStars << " star(s) centred on a cone\n";
        }
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
            if (ar.directed > 0 || ar.conePlacedStars > 0) {
                std::cout << "  Field-directed: " << ar.directed << " witness(es) where the field disagrees, "
                          << ar.conePlacedStars << " star fill(s) centred on a cone\n";
            }
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
              << "  --no-half-ogrid       no half O-grids (two corners on a straight side, e.g. an axis)\n"
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
              << "  --realise-sweeps <n>  TMOP sweeps on the realisation, 0 = none       (default 100)\n"
              << "  --realisations <n>    a coarse search's layouts tried, best first, until one\n"
              << "                        realises                                        (default 4)\n"
              << "  --last-resort <n>     ... and up to n when no coarse search realised in those (default 16)\n\n"
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
              << "  --mesh-curves <k>     what a feature node slides on: spline (the interpolant of its\n"
              << "                        run, the default), polyline (the run itself), or chord (the\n"
              << "                        line through its two feature neighbours)\n"
              << "  --mesh-obj <f.obj>    write the mesh (after TMOP when that ran)\n"
              << "  --mesh-frac <f>       as --mesh, at f times the bounding-box diagonal\n\n"
              << "The reference cross field (docs/atlas_crossfield_guidance.md)\n"
              << "  --field <k>           dualmbo (the default): a DualMBO field on the input, scored by\n"
              << "                        Stages 3, 5 and 6; harmonic: the annealer's old term only; none\n"
              << "  --field-report        E_dir, E_sing and the cone match of the answer, per-search scores,\n"
              << "                        and a SUMMARY line (a field is solved for it if the run had none)\n"
              << "  --w-dir <w>           blocks per unit of E_dir                        (default 10)\n"
              << "  --w-edge <w>          blocks per +1/4 defect on dS                    (default 1.5)\n"
              << "  --w-edge-minus <w>    blocks per -1/4 defect on dS                    (default 0.25)\n"
              << "  --w-pos <w>           an interior singularity r0 or more from its cone (default 0.5)\n"
              << "  --w-extra <w>         an interior singularity with no cone            (default 1)\n"
              << "  --w-missing <w>       a cone with no singularity                      (default 1)\n"
              << "  --r0 <f>              E_sing's length scale, fraction of the diagonal (default 0.1)\n"
              << "  --no-field-templates  Stage 3 without the field\n"
              << "  --no-cone-placement   Stage 3: no template parameters from the cones\n"
              << "  --no-field-choice     Stage 3: the first family that validates, not the best score\n"
              << "  --no-field-anneal     Stage 6's annealer without the field\n"
              << "  --no-field-greedy     Stage 6's greedy rounds without the field\n"
              << "  --field-greedy-fine   ... and with it on the fine carrier too (default: coarse only)\n"
              << "  --no-directed         no field-directed witnesses or cone-centred stars\n"
              << "  --directed <f>        fraction of the annealer's witnesses the field picks (default 0.3)\n"
              << "  --no-field-arbiter    rounds and searches compared on Stage 5's objective alone\n"
              << "  --signed-defects <p> <m>  the annealer's weight of a +1/4 and a -1/4 boundary defect\n"
              << "                        (default 4.5 1.5)\n"
              << "  --greedy-signed <p> <m>  Stage 6's greedy rounds on coarse carriers: a +1/4 and a -1/4\n"
              << "                        boundary defect against an interior one's 1 (default 1.5 0.5)\n"
              << "  --no-signed-defects   both of the above back to the symmetric weights\n"
              << "  --arbiter-shape <w>   the arbiter's weight on the worst carrier cell below SJ 0.35\n"
              << "  --no-patch-topology   the annealer counts a base patch round a hole as one block, not\n"
              << "                        the three or more it needs (as before 2026-09-19)\n"
              << "  --min-sj <s>          exit 10 unless the --mesh (after TMOP, when that ran) has worst\n"
              << "                        scaled Jacobian at least s\n";
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
    mesh::QuadMesh::Options meshNodeOpts;
    std::string meshOut;
    bool fieldReport = false;
    double meshFrac = 0.0;
    double minSJ = -2.0;

    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--corner-angle" && i + 1 < argc)        opts.domain.cornerAngle = std::stod(argv[++i]);
        else if (a == "--kink-angle" && i + 1 < argc)     opts.domain.kinkAngle = std::stod(argv[++i]);
        else if (a == "--no-interfaces")                  opts.domain.materialInterfaces = false;
        else if (a == "--no-templates")                   opts.runTemplates = false;
        else if (a == "--no-sections")                    opts.templates.sections = false;
        else if (a == "--no-stars")                       opts.templates.stars = false;
        else if (a == "--no-ogrid")                       opts.templates.ogrids = false;
        else if (a == "--no-half-ogrid")                  opts.templates.halfOGrids = false;
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
        else if (a == "--realisations" && i + 1 < argc)   opts.maxRealisations = std::stoi(argv[++i]);
        else if (a == "--last-resort" && i + 1 < argc)    opts.maxRealisationsLastResort = std::stoi(argv[++i]);
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
        else if (a == "--mesh-curves" && i + 1 < argc) {
            const std::string k = argv[++i];
            if (k == "chord")         meshNodeOpts.curveSource = mesh::QuadMesh::Options::CurveChord;
            else if (k == "polyline") meshNodeOpts.curveSource = mesh::QuadMesh::Options::CurvePolyline;
            else if (k == "spline")   meshNodeOpts.curveSource = mesh::QuadMesh::Options::CurveSpline;
            else { std::cerr << "Unknown --mesh-curves: " << k << "\n"; usage(argv[0]); return 1; }
        }
        else if (a == "--mesh-obj" && i + 1 < argc)       meshOut = argv[++i];
        else if (a == "--mesh-frac" && i + 1 < argc)      meshFrac = std::stod(argv[++i]);
        else if (a == "--field" && i + 1 < argc) {
            const std::string k = argv[++i];
            if (k == "harmonic")     opts.field.reference = ATLAS::Options::Field::Reference::Harmonic;
            else if (k == "dualmbo") opts.field.reference = ATLAS::Options::Field::Reference::DualMBO;
            else if (k == "none")    opts.field.reference = ATLAS::Options::Field::Reference::None;
            else { std::cerr << "Unknown --field: " << k << "\n"; usage(argv[0]); return 1; }
        }
        else if (a == "--field-report")                   fieldReport = true;
        else if (a == "--w-dir" && i + 1 < argc)          opts.field.wDir = std::stod(argv[++i]);
        else if (a == "--w-edge" && i + 1 < argc)         opts.field.singularity.wEdgePlus = std::stod(argv[++i]);
        else if (a == "--w-edge-minus" && i + 1 < argc)   opts.field.singularity.wEdgeMinus = std::stod(argv[++i]);
        else if (a == "--w-pos" && i + 1 < argc)          opts.field.singularity.wPos = std::stod(argv[++i]);
        else if (a == "--w-extra" && i + 1 < argc)        opts.field.singularity.wExtra = std::stod(argv[++i]);
        else if (a == "--w-missing" && i + 1 < argc)      opts.field.singularity.wMissing = std::stod(argv[++i]);
        else if (a == "--r0" && i + 1 < argc)             opts.field.singularity.r0 = std::stod(argv[++i]);
        else if (a == "--no-field-templates")             opts.field.templates = false;
        else if (a == "--no-cone-placement")              opts.templates.conePlacement = false;
        else if (a == "--no-field-choice")                opts.templates.chooseByField = false;
        else if (a == "--no-field-anneal")                opts.field.anneal = false;
        else if (a == "--no-field-greedy")                opts.field.greedy = false;
        else if (a == "--field-greedy-fine")              opts.field.greedyOnFine = true;
        else if (a == "--no-field-arbiter")               opts.field.arbiter = false;
        else if (a == "--no-directed")                    opts.field.directedMoves = false;
        else if (a == "--directed" && i + 1 < argc)       opts.field.directedFraction = std::stod(argv[++i]);
        else if (a == "--signed-defects" && i + 2 < argc) {
            opts.anneal.wBoundaryDefectPlus = std::stod(argv[++i]);
            opts.anneal.wBoundaryDefectMinus = std::stod(argv[++i]);
        }
        else if (a == "--arbiter-shape" && i + 1 < argc)  opts.arbiterShape = std::stod(argv[++i]);
        else if (a == "--patch-topology")                 opts.anneal.patchTopology = true;
        else if (a == "--no-patch-topology")              opts.anneal.patchTopology = false;
        else if (a == "--no-signed-defects") {
            opts.anneal.wBoundaryDefectPlus = opts.anneal.wBoundaryDefectMinus = -1.0;
            opts.rewrite.boundaryPlusWeight = opts.rewrite.boundaryMinusWeight = 1.0;
        }
        else if (a == "--min-sj" && i + 1 < argc)         minSJ = std::stod(argv[++i]);
        else if (a == "--greedy-signed" && i + 2 < argc) {
            opts.rewrite.boundaryPlusWeight = std::stod(argv[++i]);
            opts.rewrite.boundaryMinusWeight = std::stod(argv[++i]);
        }
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

    showField = fieldReport || opts.field.reference == ATLAS::Options::Field::Reference::DualMBO ||
                opts.arbiterShape > 0.0;
    if (meshFrac > 0.0 && !mesh->vertices.empty()) {
        Point lo = mesh->vertices[0], hi = lo;
        for (const Point &p : mesh->vertices) {
            lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
            hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
        }
        meshOpts.targetEdgeLength = meshFrac * std::hypot(hi[0] - lo[0], hi[1] - lo[1]);
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
    // The reference field against the chosen answer (the guidance note's
    // Sec. 2 probe). With --field harmonic|none a field is solved here for
    // the report only; the searches never saw it.
    std::ostringstream summary;
    if (showField) {
        heading("The reference cross field (docs/atlas_crossfield_guidance.md)");
        std::unique_ptr<ReferenceField> own;
        const ReferenceField *F = pipe.getField();
        if (!F) {
            std::vector<int> interfaces;
            const PlanarDomain &D = pipe.getDomain();
            for (int e = 0; e < static_cast<int>(D.interfaceEdge.size()); ++e) if (D.interfaceEdge[e]) interfaces.push_back(e);
            own = std::make_unique<ReferenceField>(mesh, interfaces, opts.field.solve);
            F = own.get();
            std::cout << "  Solved for this report only: the searches ran without it\n";
        }
        const ReferenceField::Report &fr = F->getReport();
        std::cout << "  DualMBO: " << fr.levels << " tau level(s), " << fr.steps << " MBO step(s), "
                  << (fr.converged ? "converged" : "not converged") << ", " << std::fixed << std::setprecision(2)
                  << fr.seconds << " s" << std::defaultfloat << "; lookup grid " << fr.gridNx << " x " << fr.gridNy << "\n";
        std::cout << "  Cones: +" << fr.rawPlus << "/-" << fr.rawMinus << ", " << fr.dipoleUnits
                  << " dipole(s) cancelled -> +" << fr.conesPlus << "/-" << fr.conesMinus << "\n";
        ReferenceField::SingularityReport sr;
        const double eDir = F->directionEnergy(C);
        const double eSing = F->singularityEnergy(C, opts.field.singularity, &sr);
        std::cout << "  The answer: E_dir " << std::fixed << std::setprecision(4) << eDir << " (x w_dir "
                  << std::setprecision(1) << opts.field.wDir << " = " << std::setprecision(2) << opts.field.wDir * eDir
                  << " block(s)), E_sing " << eSing << std::defaultfloat << "\n";
        std::cout << "  Singular vertices: +" << sr.interiorPlus << "/-" << sr.interiorMinus << " inside, +"
                  << sr.boundaryPlus << "/-" << sr.boundaryMinus << " on dS; " << sr.matched
                  << " matched to a cone (mean " << std::fixed << std::setprecision(2) << sr.meanDistance
                  << " r0), " << sr.absorbed << " cone(s) pushed to dS, " << sr.missing << " missing, "
                  << sr.extra << " extra, " << sr.dipoles << " carrier dipole(s)" << std::defaultfloat << "\n";
        for (int i = 0; i < pipe.numSearches(); ++i) {
            const ATLAS::Search &sx = pipe.getSearch(i);
            if (!sx.succeeded()) continue;
            std::cout << "  " << sx.name << ": " << sx.finalCover()->getReport().blocks << " block(s), objective "
                      << std::fixed << std::setprecision(2) << sx.finalCover()->getReport().objective << ", score "
                      << sx.finalScore << " (E_dir " << std::setprecision(3) << sx.finalDir << ", E_sing "
                      << std::setprecision(2) << sx.finalSing << ", shape " << sx.finalShape << ")"
                      << std::defaultfloat << (i == st.chosen ? "  <- chosen" : "") << "\n";
        }
        summary << std::fixed << std::setprecision(4) << " eDir=" << eDir << " eSing=" << eSing
                << " iPlus=" << sr.interiorPlus << " iMinus=" << sr.interiorMinus << " bPlus=" << sr.boundaryPlus
                << " bMinus=" << sr.boundaryMinus << " matched=" << sr.matched << " missing=" << sr.missing
                << " extra=" << sr.extra << " conesPlus=" << fr.conesPlus << " conesMinus=" << fr.conesMinus
                << " fieldSec=" << fr.seconds;
    }

    // ---------------------------------------------------------------------
    bool meshValid = true;
    double tfiWorst = 0.0, tmopWorst = 0.0;
    int tfiFolded = -1, tmopFolded = -1;
    if (meshOpts.targetEdgeLength > 0.0) {
        heading("A quadrilateral mesh on the blocks (Secs. 11.2, 11.3)");
        BlockMesh bm(pipe.getCover(), meshOpts);
        printMesh(bm);
        meshValid = bm.getReport().valid;
        tfiWorst = bm.getReport().minScaledJacobian;
        tfiFolded = bm.getReport().invertedQuads;
        mesh::QuadMesh qm = mesh::QuadMesh::from(bm, meshNodeOpts);
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
            std::cout << "  Feature nodes: " << tr.curveNodes << " of " << tr.slidingNodes
                      << " sliding on " << tr.featureCurves << " curve(s), " << tr.fittedCurves
                      << " of them interpolants (bow " << std::scientific << std::setprecision(2)
                      << tr.curveBow << ", off-curve " << tr.curveDeviation << ")"
                      << std::defaultfloat << "\n";
            verdict(tr.curveDeviation < 1e-9 * std::max(1.0, bm.getReport().modelExtent),
                    "Every feature node ended on the curve it was bound to");
            verdict(improved && tr.invertedAfter == 0, "TMOP left the mesh better than it found it, nothing folded");
            // The mesh handed on is the smoothed one: it is what must not fold.
            meshValid = tr.invertedAfter == 0;
            tmopWorst = tr.minScaledJacobianAfter;
            tmopFolded = tr.invertedAfter;
        }
        if (!meshOut.empty()) {
            if (qm.writeOBJ(meshOut)) std::cout << "  Wrote the mesh to " << meshOut << "\n";
            else warn("Failed to write " + meshOut);
        }
    }

    if (showField) {
        // One line for corpus sweeps.
        std::cout << "\nSUMMARY model=" << path << " blocks=" << br.blocks << std::fixed << std::setprecision(4)
                  << " objective=" << br.objective << " score=" << st.bestScore << " search=\"" << chosen.name
                  << "\" carrierSJ=" << C.getReport().minScaledJacobian << summary.str();
        if (tfiFolded >= 0) std::cout << " tfiSJ=" << tfiWorst << " tfiFolds=" << tfiFolded;
        if (tmopFolded >= 0) std::cout << " tmopSJ=" << tmopWorst << " tmopFolds=" << tmopFolded;
        std::cout << " searchSec=" << st.secondsSearch << " fieldStageSec=" << st.secondsField
                  << std::defaultfloat << "\n";
    }

    heading("Result");
    for (const std::string &m : st.messages) {
        if (m.rfind("Stopping", 0) == 0 || m.rfind("Stages 4-6", 0) == 0 || m.rfind("Stage 1b", 0) == 0) warn(m);
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
    if (ok && minSJ > -2.0) {
        const double worst = tmopFolded >= 0 ? tmopWorst : tfiWorst;
        const bool good = tfiFolded >= 0 && worst >= minSJ;
        std::ostringstream os;
        os << "The mesh's worst scaled Jacobian " << std::fixed << std::setprecision(3) << worst << " is at least "
           << minSJ << " (--min-sj)";
        verdict(good, os.str());
        if (!good) return 10;
    }
    return ok ? 0 : 7;
}
