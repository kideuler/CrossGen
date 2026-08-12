// Checks on the QGP quantizer (Campen et al. 2015, Sec. 6, interval-
// assignment variant): a hand-built T-mesh with a T-junction, the boundary
// (source-to-sink) generating vectors, Stage II shrinking toward xIdeal = 1,
// and the geometric block welder that feeds MedialAxisTMesh blocks in.
// Usage: TestQuantization

#include <iostream>
#include <string>

#include "medialaxis/MedialAxisTMesh.hxx"
#include "quantization/QuantTMesh.hxx"
#include "quantization/QuantTMeshConvert.hxx"
#include "quantization/TMeshQuantizer.hxx"

#define GREEN "\033[32m"
#define RED   "\033[31m"
#define RESET "\033[0m"

static int failures = 0;

static void check(bool ok, const std::string &what) {
    std::cout << (ok ? GREEN "  PASS " : RED "  FAIL ") << RESET << what
              << "\n";
    if (!ok) ++failures;
}

// The smallest genuine T-mesh: a unit-square domain of three faces, the top
// face C spanning the two bottom faces A and B, so the node between them is
// a T-junction on C's lower side.
//
//   +-----t-----+
//   |     C     |
//   m1----+----m2
//   | A   |  B  |
//   a1----+----a2
//
static QuantTMesh buildTJunction(double idealTop) {
    QuantTMesh tm;
    const int a1 = tm.addEdge(), a2 = tm.addEdge();   // bottom boundary
    const int m1 = tm.addEdge(), m2 = tm.addEdge();   // interior mid line
    const int t = tm.addEdge(idealTop);               // top boundary
    const int L0 = tm.addEdge(), M0 = tm.addEdge(), R0 = tm.addEdge();
    const int L1 = tm.addEdge(), R1 = tm.addEdge();

    tm.addFace({a1}, {M0}, {m1}, {L0});       // A
    tm.addFace({a2}, {R0}, {m2}, {M0});       // B
    tm.addFace({m1, m2}, {R1}, {t}, {L1});    // C
    return tm;
}

static void testTJunction() {
    std::cout << "T-junction block decomposition (xIdeal = 1):\n";
    QuantTMesh tm = buildTJunction(1.0);
    std::string error;
    check(tm.finalize(&error), "finalize: " + error);
    check(tm.edges[0].onBoundary() && !tm.edges[2].onBoundary(),
          "boundary/interior classification");

    TMeshQuantizer::Report report = TMeshQuantizer(tm).run();
    check(report.consistent, "quantization consistent (Ax = 0)");
    bool positive = true;
    for (const auto &e : tm.edges) positive = positive && e.x >= 1;
    check(positive, "every edge >= 1");
    check(tm.sideSum(2, 0) == tm.sideSum(2, 2),
          "T-junction side balances the top");
    // Minimal quantization: the split sides force the top edge to 2 and
    // everything else to 1, objective 1.
    check(tm.edges[4].x == 2, "top edge quantized to 2, got " +
                                  std::to_string(tm.edges[4].x));
    check(report.objective <= 1.0 + 1e-12,
          "Stage II reached the minimal objective");
}

static void testStageIIGrows() {
    std::cout << "Stage II grows toward a sizing field (h = 1/2):\n";
    // Same topology with every ideal doubled, as a target element size half
    // the block size would set them: Stage II must reach the exact sizing,
    // objective 0.
    QuantTMesh tm;
    for (int e = 0; e < 10; ++e) tm.addEdge(e == 4 ? 4.0 : 2.0);
    tm.addFace({0}, {6}, {2}, {5});       // A
    tm.addFace({1}, {7}, {3}, {6});       // B
    tm.addFace({2, 3}, {9}, {4}, {8});    // C
    std::string error;
    check(tm.finalize(&error), "finalize: " + error);
    TMeshQuantizer::Report report = TMeshQuantizer(tm).run();
    check(report.consistent, "quantization consistent (Ax = 0)");
    check(tm.edges[4].x == 4,
          "top edge reaches its ideal 4, got " + std::to_string(tm.edges[4].x));
    check(report.objective == 0.0, "sizing met exactly, objective 0");
}

// The same three-face configuration as geometry only, plus corner counting:
// the welder has to find the T-junction on C's lower side by itself.
static void testBlockWelder() {
    std::cout << "Geometric welder on three blocks:\n";
    auto quad = [](Point p0, Point p1, Point p2, Point p3) {
        TMeshBlock b;
        b.outline = {p0, p1, p2, p3};
        b.corners = b.outline;
        return b;
    };
    std::vector<TMeshBlock> blocks;
    blocks.push_back(quad({0, 0}, {1, 0}, {1, 1}, {0, 1}));   // A
    blocks.push_back(quad({1, 0}, {2, 0}, {2, 1}, {1, 1}));   // B
    blocks.push_back(quad({0, 1}, {2, 1}, {2, 2}, {0, 2}));   // C

    BlockQuant welded = makeQuantTMesh(blocks, 1e-9);
    check(welded.ok, "welded T-mesh finalizes: " + welded.error);
    check(welded.skippedBlocks == 0, "no block skipped");
    check(welded.nodes.size() == 8, "8 welded nodes, got " +
                                        std::to_string(welded.nodes.size()));
    // 10 edges: C's lower side must have been split at (1,1).
    check(welded.tmesh.edges.size() == 10,
          "10 edges after splitting at the T-junction, got " +
              std::to_string(welded.tmesh.edges.size()));

    TMeshQuantizer::Report report = TMeshQuantizer(welded.tmesh).run();
    check(report.consistent, "welded quantization consistent");
    long long top = 0;
    for (size_t e = 0; e < welded.tmesh.edges.size(); ++e) {
        // The top side of C is the single edge from (0,2) to (2,2).
        const auto &g = welded.edgeGeometry[e];
        if (g.size() >= 2 && g.front()[1] == 2.0 && g.back()[1] == 2.0) {
            top = welded.tmesh.edges[e].x;
        }
    }
    check(top == 2, "top edge of C quantized to 2, got " + std::to_string(top));
}

// A triangular block is not a T-mesh cell, but dropping it would leave a
// hole in the decomposition. The converter halves its longest side to make
// a fourth, the way the paper's templates carve a cap.
static void testTrianglePromotion() {
    std::cout << "Triangular block promoted to a four-sided cell:\n";
    TMeshBlock tri;
    tri.outline = {{0, 0}, {4, 0}, {0, 3}};
    tri.corners = tri.outline;
    TMeshBlock quad;
    quad.outline = {{4, 0}, {4, 3}, {0, 3}};   // shares the hypotenuse
    quad.corners = {{4, 0}, {4, 3}, {0, 3}};

    BlockQuant welded = makeQuantTMesh({tri, quad}, 1e-9);
    check(welded.ok, "welded T-mesh finalizes: " + welded.error);
    check(welded.skippedBlocks == 0,
          "no block dropped, got " + std::to_string(welded.skippedBlocks));
    check(welded.tmesh.faces.size() == 2,
          "both blocks became cells, got " +
              std::to_string(welded.tmesh.faces.size()));
    // The split lands on the longest side, the 5-long hypotenuse, so its
    // midpoint joins the 4 corners already present.
    check(welded.nodes.size() == 5,
          "5 nodes after the split, got " + std::to_string(welded.nodes.size()));

    TMeshQuantizer::Report report = TMeshQuantizer(welded.tmesh).run();
    check(report.consistent, "quantization consistent (Ax = 0)");
    bool positive = true;
    for (const auto &e : welded.tmesh.edges) positive = positive && e.x >= 1;
    check(positive, "every edge >= 1");
    for (size_t f = 0; f < welded.tmesh.faces.size(); ++f) {
        check(welded.tmesh.sideSum(f, 0) == welded.tmesh.sideSum(f, 2) &&
                  welded.tmesh.sideSum(f, 1) == welded.tmesh.sideSum(f, 3),
              "face " + std::to_string(f) + " has matching opposite sides");
        // What the renderer needs: opposite sides carry the same tick count.
        check(sideCurve(welded, f, 0).cells() == sideCurve(welded, f, 2).cells() &&
                  sideCurve(welded, f, 1).cells() == sideCurve(welded, f, 3).cells(),
              "face " + std::to_string(f) + " grids weld");
    }
}

int main() {
    testTJunction();
    testStageIIGrows();
    testBlockWelder();
    testTrianglePromotion();
    std::cout << (failures == 0 ? GREEN "\nAll checks passed\n" RESET
                                : RED "\nFAILURES\n" RESET);
    return failures == 0 ? 0 : 1;
}
