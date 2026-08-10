// Structural checks on the discrete medial axis built from an .obj mesh:
// dual-complex integrity, boundary orientation, the inside filter, cell
// merging, and branch extraction.
// Usage: TestMedialAxis <mesh.obj>

#include <iostream>
#include <iomanip>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <cmath>

#include "mesh/Mesh.hxx"
#include "triangle/TriangleMesher.hpp"
#include "medialaxis/MedialAxis.hxx"

// ANSI colour helpers
#define GREEN "\033[32m"
#define RED   "\033[31m"
#define YELLOW "\033[33m"
#define RESET "\033[0m"

// ─── Delaunay re-triangulation (same logic as ViewerMain.cxx) ───────────────
static std::shared_ptr<Mesh> buildDelaunayMesh(const std::shared_ptr<Mesh>& mesh) {
    // 1) Collect ordered boundary loops
    std::unordered_map<int, std::vector<int>> bAdj;
    for (int beIdx : mesh->boundaryEdges) {
        int a = mesh->edges[beIdx][0];
        int b = mesh->edges[beIdx][1];
        bAdj[a].push_back(b);
        bAdj[b].push_back(a);
    }

    std::unordered_set<int> visited;
    std::vector<std::vector<int>> loops;
    for (int bv : mesh->boundaryVertices) {
        if (visited.count(bv)) continue;
        std::vector<int> loop;
        int prev = -1, curr = bv;
        while (true) {
            visited.insert(curr);
            loop.push_back(curr);
            int next = -1;
            for (int nb : bAdj[curr]) {
                if (nb != prev && !visited.count(nb)) { next = nb; break; }
            }
            if (next == -1) break;
            prev = curr;
            curr = next;
        }
        if (loop.size() >= 3) loops.push_back(std::move(loop));
    }

    if (loops.empty()) {
        std::cerr << "No boundary loops found!\n";
        return nullptr;
    }

    // 2) Build TriangleMesher input
    std::unordered_map<int, int> vertRemap;
    std::vector<std::array<double, 2>> vertlist;
    for (const auto& loop : loops) {
        for (int vi : loop) {
            if (!vertRemap.count(vi)) {
                int newIdx = static_cast<int>(vertlist.size());
                vertRemap[vi] = newIdx;
                vertlist.push_back(mesh->vertices[vi]);
            }
        }
    }

    std::vector<std::vector<std::array<int, 2>>> segment_loops;
    std::vector<int> loopTypes;
    for (size_t li = 0; li < loops.size(); ++li) {
        const auto& loop = loops[li];
        std::vector<std::array<int, 2>> segments;
        for (size_t i = 0; i < loop.size(); ++i) {
            int a = vertRemap[loop[i]];
            int b = vertRemap[loop[(i + 1) % loop.size()]];
            segments.push_back({a, b});
        }
        segment_loops.push_back(std::move(segments));
        loopTypes.push_back(li == 0 ? 0 : 1);
    }

    // 3) Triangulate
    triangle_wrapper::TriangleMesher2D::Options opts;
    opts.just_delaunay = true;
    triangle_wrapper::TriangleMesher2D mesher(opts);

    triangle_wrapper::TriangleMesher2D::MeshInput input;
    input.vertlist = vertlist;
    input.segment_loops = segment_loops;
    input.type = loopTypes;
    input.h = 0.0;

    auto output = mesher.triangulate(input);
    return std::make_shared<Mesh>(output.verts, output.triangles);
}

// ─── Check harness ──────────────────────────────────────────────────────────
namespace {

int g_fails = 0;

void check(bool ok, const std::string& label, const std::string& detail = "") {
    std::cout << (ok ? GREEN "[PASS]" RESET : RED "[FAIL]" RESET) << " " << label;
    if (!detail.empty()) std::cout << " (" << detail << ")";
    std::cout << "\n";
    if (!ok) ++g_fails;
}

// Every structural invariant the dual complex is supposed to maintain. Run
// before and after simplification, since a merge is the operation most likely
// to break one.
void checkComplexIntegrity(const MedialAxis& ma, const std::string& when) {
    int badRing = 0, badKey = 0, asymmetric = 0, badCount = 0, repeatedVert = 0;
    std::unordered_map<int, int> observedVertexCells;

    for (int ci = 0; ci < static_cast<int>(ma.cells.size()); ++ci) {
        const auto& cell = ma.cells[ci];
        if (!cell.active) continue;

        const int n = static_cast<int>(cell.verts.size());
        if (n < 3) { ++badRing; continue; }

        std::unordered_set<int> seen;
        for (int v : cell.verts) {
            if (!seen.insert(v).second) ++repeatedVert;
            ++observedVertexCells[v];
        }

        for (int i = 0; i < n; ++i) {
            const auto e = ma.ringEdge(ci, i);
            if (!ma.edgeCells.count(MeshEdgeKey(e[0], e[1]))) { ++badKey; continue; }

            const int nb = ma.neighborAcross(ci, i);
            if (nb < 0) continue;
            // The neighbour must name this cell back across the same edge.
            bool mutual = false;
            const int nn = static_cast<int>(ma.cells[nb].verts.size());
            for (int j = 0; j < nn; ++j) {
                if (ma.neighborAcross(nb, j) == ci) { mutual = true; break; }
            }
            if (!mutual) ++asymmetric;
        }
    }

    for (const auto& kv : observedVertexCells) {
        if (ma.vertexCellCount(kv.first) != kv.second) ++badCount;
    }

    check(badRing == 0, when + ": every active cell has >= 3 ring vertices",
          std::to_string(badRing) + " bad");
    check(repeatedVert == 0, when + ": cell rings are simple (no repeated vertex)",
          std::to_string(repeatedVert) + " repeats");
    check(badKey == 0, when + ": every ring edge is present in edgeCells",
          std::to_string(badKey) + " missing");
    check(asymmetric == 0, when + ": cell adjacency is symmetric",
          std::to_string(asymmetric) + " one-sided");
    check(badCount == 0, when + ": vertexCellCount matches the rings",
          std::to_string(badCount) + " stale");
}

void checkAxisConsistency(const MedialAxis& ma, const std::string& when) {
    const int nv = static_cast<int>(ma.medialVertices.size());
    int oob = 0, nonFinite = 0, badBackref = 0, degreeMismatch = 0, selfLoop = 0;

    for (const auto& e : ma.medialEdges) {
        if (e[0] < 0 || e[0] >= nv || e[1] < 0 || e[1] >= nv) ++oob;
        else if (e[0] == e[1]) ++selfLoop;
    }

    for (int i = 0; i < nv; ++i) {
        const auto& mv = ma.medialVertices[i];
        if (!std::isfinite(mv.coord[0]) || !std::isfinite(mv.coord[1])) ++nonFinite;
        if (mv.degree != static_cast<int>(mv.incidentEdges.size())) ++degreeMismatch;

        // Each incident edge must actually name this vertex, and the dual cell
        // must point back at it.
        for (int e : mv.incidentEdges) {
            if (e < 0 || e >= static_cast<int>(ma.medialEdges.size())) { ++badBackref; continue; }
            if (ma.medialEdges[e][0] != i && ma.medialEdges[e][1] != i) ++badBackref;
        }
        if (mv.cell < 0 || mv.cell >= static_cast<int>(ma.cells.size()) ||
            ma.cells[mv.cell].medialVertex != i) {
            ++badBackref;
        }
    }

    check(oob == 0, when + ": medial edge indices in range", std::to_string(oob) + " bad");
    check(selfLoop == 0, when + ": no self-loop medial edges", std::to_string(selfLoop) + " found");
    check(nonFinite == 0, when + ": medial vertex coordinates finite",
          std::to_string(nonFinite) + " bad");
    check(degreeMismatch == 0, when + ": degree == incidentEdges.size()",
          std::to_string(degreeMismatch) + " bad");
    check(badBackref == 0, when + ": incidentEdges / cell back-references agree",
          std::to_string(badBackref) + " bad");
}

} // namespace

// ─── Main ───────────────────────────────────────────────────────────────────
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <mesh.obj>\n";
        return 1;
    }

    const std::string path = argv[1];
    std::cout << "Loading mesh: " << path << "\n";

    auto originalMesh = std::make_shared<Mesh>(path);
    std::cout << "  Original mesh: " << originalMesh->triangles.size() << " triangles, "
              << originalMesh->vertices.size() << " vertices, "
              << originalMesh->boundaryEdges.size() << " boundary edges\n";

    auto delaunayMesh = buildDelaunayMesh(originalMesh);
    if (!delaunayMesh) {
        std::cerr << RED "[FAIL]" RESET " Could not build Delaunay mesh.\n";
        return 1;
    }
    std::cout << "  Delaunay mesh:  " << delaunayMesh->triangles.size() << " triangles, "
              << delaunayMesh->vertices.size() << " vertices, "
              << delaunayMesh->boundaryEdges.size() << " boundary edges\n";

    MedialAxis ma(delaunayMesh);
    std::cout << "  Dual cells:      " << ma.cells.size() << "\n";
    std::cout << "  Medial vertices: " << ma.medialVertices.size() << "\n";
    std::cout << "  Medial edges:    " << ma.medialEdges.size() << "\n";
    std::cout << "  Stats: " << ma.stats.degenerateCells << " degenerate, "
              << ma.stats.centersOutsideDomain << " centers outside, "
              << ma.stats.edgesCrossingBoundary << " edges crossing boundary, "
              << ma.stats.nonManifoldBoundaryVertices << " pinched boundary vertices\n";

    std::cout << "\n── Raw complex ──\n";

    check(ma.cells.size() == delaunayMesh->triangles.size(),
          "one dual cell per triangle",
          std::to_string(ma.cells.size()) + " vs " + std::to_string(delaunayMesh->triangles.size()));

    check(ma.medialVertices.size() == ma.cells.size(),
          "one medial vertex per active cell",
          std::to_string(ma.medialVertices.size()) + " vs " + std::to_string(ma.cells.size()));

    // The inside filter can only remove candidates, never invent them.
    {
        int numInternal = 0;
        for (const auto& et : delaunayMesh->edgeTriangles) {
            if (et[0] != -1 && et[1] != -1) ++numInternal;
        }
        const int kept = static_cast<int>(ma.medialEdges.size());
        check(kept + ma.stats.edgesCrossingBoundary == numInternal,
              "kept edges + filtered edges == internal Delaunay edges",
              std::to_string(kept) + " + " + std::to_string(ma.stats.edgesCrossingBoundary) +
              " vs " + std::to_string(numInternal));
    }

    // Cells are normalised to CCW, which is what makes ring order, incident
    // edge order and boundary orientation agree.
    {
        int cw = 0;
        for (const auto& cell : ma.cells) {
            if (!cell.active) continue;
            double area2 = 0.0;
            const int n = static_cast<int>(cell.verts.size());
            for (int i = 0; i < n; ++i) {
                const Point& p = delaunayMesh->vertices[cell.verts[i]];
                const Point& q = delaunayMesh->vertices[cell.verts[(i + 1) % n]];
                area2 += cross2(p, q);
            }
            if (area2 <= 0.0) ++cw;
        }
        check(cw == 0, "every cell ring is CCW", std::to_string(cw) + " reversed");
    }

    checkComplexIntegrity(ma, "raw");
    checkAxisConsistency(ma, "raw");

    // ── Boundary orientation ──
    // next/prev must keep the interior on the left. That makes the outer loop
    // CCW and hole loops CW, and it makes every interior angle land in
    // (0, 2*pi) with sum (k-2)*pi over an outer loop.
    {
        std::vector<bool> seen(delaunayMesh->vertices.size(), false);
        int loops = 0, ccwLoops = 0, badAngleSum = 0, unlinked = 0;
        double totalArea = 0.0;

        for (int bv : delaunayMesh->boundaryVertices) {
            if (ma.boundaryLinks[bv].nextBoundaryVertex == -1 ||
                ma.boundaryLinks[bv].prevBoundaryVertex == -1) {
                ++unlinked;
            }
        }

        for (int bv : delaunayMesh->boundaryVertices) {
            if (seen[bv] || ma.boundaryLinks[bv].nextBoundaryVertex == -1) continue;
            ++loops;

            double area2 = 0.0, angleSum = 0.0;
            int count = 0, cur = bv;
            do {
                seen[cur] = true;
                const int nxt = ma.boundaryLinks[cur].nextBoundaryVertex;
                if (nxt == -1) break;
                area2 += cross2(delaunayMesh->vertices[cur], delaunayMesh->vertices[nxt]);
                angleSum += ma.interiorAngle[cur];
                ++count;
                cur = nxt;
            } while (cur != bv && count < static_cast<int>(delaunayMesh->vertices.size()) + 1);

            totalArea += 0.5 * area2;
            if (area2 > 0.0) ++ccwLoops;
            // The domain-side angles of a CCW outer loop sum to (k-2)*pi. Around
            // a hole the domain is on the outside, so its k reflex angles sum
            // to (k+2)*pi instead.
            const double expected = (area2 > 0.0 ? (count - 2.0) : (count + 2.0)) * M_PI;
            if (std::abs(angleSum - expected) > 1e-6 * std::max(1.0, std::abs(expected))) {
                ++badAngleSum;
            }
        }

        check(unlinked == 0, "every boundary vertex has next and prev",
              std::to_string(unlinked) + " unlinked");
        check(loops > 0, "boundary loops found", std::to_string(loops) + " loops");
        check(ccwLoops == 1, "exactly one CCW loop (the outer one; holes run CW)",
              std::to_string(ccwLoops) + " CCW of " + std::to_string(loops));
        check(totalArea > 0.0, "signed area of all loops is positive (interior on the left)",
              std::to_string(totalArea));
        check(badAngleSum == 0, "interior angles sum to (k-2)*pi per loop",
              std::to_string(badAngleSum) + " loops off");
    }

    // Medial vertices must lie inside the domain to be medial at all.
    {
        int outside = 0;
        for (const auto& cell : ma.cells) {
            if (cell.active && !cell.centerInside) ++outside;
        }
        check(outside == ma.stats.centersOutsideDomain,
              "centersOutsideDomain matches the cell flags",
              std::to_string(outside) + " vs " + std::to_string(ma.stats.centersOutsideDomain));
        std::cout << (outside == 0 ? GREEN "[INFO]" RESET : YELLOW "[WARN]" RESET)
                  << " " << outside << " circumcenters lie outside the domain"
                  << (outside ? "; the boundary sampling is too coarse for a faithful axis" : "")
                  << "\n";
    }

    // ── Simplification ──
    std::cout << "\n── After deduplication ──\n";
    const size_t rawVertices = ma.medialVertices.size();
    ma.deduplicateMedialVertices();

    int activeCells = 0;
    int polygonalCells = 0;
    double worstDeviation = 0.0;
    for (const auto& cell : ma.cells) {
        if (!cell.active) continue;
        ++activeCells;
        if (!cell.dualIsTriangle()) ++polygonalCells;
        worstDeviation = std::max(worstDeviation, cell.maxRadiusDeviation);
    }
    std::cout << "  Active cells:    " << activeCells << " (" << polygonalCells
              << " non-triangular)\n";
    std::cout << "  Medial vertices: " << rawVertices << " -> " << ma.medialVertices.size() << "\n";
    std::cout << "  Medial edges:    " << ma.medialEdges.size() << "\n";
    std::cout << "  Merges:          " << ma.stats.mergedCells << " applied, "
              << ma.stats.mergesRejectedMultiEdge << " rejected\n";
    std::cout << "  Worst cocircularity drift: " << worstDeviation << " of the medial radius\n";

    check(static_cast<int>(ma.medialVertices.size()) == activeCells,
          "one medial vertex per surviving cell");
    check(static_cast<int>(ma.cells.size()) - activeCells == ma.stats.mergedCells,
          "deactivated cells == merges applied",
          std::to_string(static_cast<int>(ma.cells.size()) - activeCells) + " vs " +
          std::to_string(ma.stats.mergedCells));

    checkComplexIntegrity(ma, "post-dedup");
    checkAxisConsistency(ma, "post-dedup");

    // Merging must not fuse two circumcenters that were meant to stay apart,
    // nor leave two survivors sitting on top of each other.
    {
        const double minSep = 1e-9 * ma.boundingBoxDiagonal();
        int coincident = 0;
        for (const auto& e : ma.medialEdges) {
            if (normP(ma.medialVertices[e[0]].coord - ma.medialVertices[e[1]].coord) < minSep) {
                ++coincident;
            }
        }
        check(coincident == 0, "no medial edge has coincident endpoints",
              std::to_string(coincident) + " degenerate");
    }

    // ── Branches ──
    std::cout << "\n── Branches ──\n";
    ma.createPolylines();

    int cyclic = 0, pureCycles = 0;
    for (size_t i = 0; i < ma.polyLines.size(); ++i) {
        if (!ma.polyLineIsCycle[i]) continue;
        ++cyclic;
        // A cycle made only of regular vertices has no seed of degree != 2, so
        // it is reachable only by the dedicated cycle pass.
        if (ma.medialVertices[ma.polyLines[i].front()].degree == 2) ++pureCycles;
    }
    std::cout << "  Polylines: " << ma.polyLines.size() << " (" << cyclic << " cyclic, "
              << pureCycles << " of them pure degree-2 cycles)\n";

    check(ma.polyLines.size() == ma.polyLineIsCycle.size(),
          "polyLineIsCycle is parallel to polyLines");

    // Every medial edge must appear in exactly one branch. This is the check
    // that catches both a branch traced twice and a cycle never traced at all.
    {
        std::map<std::pair<int, int>, int> edgeUse;
        for (const auto& e : ma.medialEdges) {
            edgeUse[{std::min(e[0], e[1]), std::max(e[0], e[1])}] = 0;
        }
        int unknown = 0;
        for (const auto& pl : ma.polyLines) {
            for (size_t j = 1; j < pl.size(); ++j) {
                const auto key = std::make_pair(std::min(pl[j - 1], pl[j]),
                                                std::max(pl[j - 1], pl[j]));
                auto it = edgeUse.find(key);
                if (it == edgeUse.end()) ++unknown;
                else ++it->second;
            }
        }
        int never = 0, twice = 0;
        for (const auto& kv : edgeUse) {
            if (kv.second == 0) ++never;
            else if (kv.second > 1) ++twice;
        }
        check(unknown == 0, "branches only use real medial edges",
              std::to_string(unknown) + " phantom");
        check(never == 0, "every medial edge is covered by a branch",
              std::to_string(never) + " uncovered");
        check(twice == 0, "no medial edge is covered twice",
              std::to_string(twice) + " duplicated");
    }

    // Interior vertices of an open branch must be regular; a cycle is closed.
    {
        int badInterior = 0, badCycle = 0;
        for (size_t i = 0; i < ma.polyLines.size(); ++i) {
            const auto& pl = ma.polyLines[i];
            for (size_t j = 1; j + 1 < pl.size(); ++j) {
                if (ma.medialVertices[pl[j]].degree != 2) ++badInterior;
            }
            if (ma.polyLineIsCycle[i]) {
                if (pl.size() < 3 || pl.front() != pl.back()) ++badCycle;
            } else if (pl.size() > 1 && pl.front() == pl.back()) {
                ++badCycle;
            }
        }
        check(badInterior == 0, "branch interiors are all degree-2 vertices",
              std::to_string(badInterior) + " bad");
        check(badCycle == 0, "cycle flag matches closure", std::to_string(badCycle) + " bad");
    }

    // ── Summary ──
    std::cout << "\n";
    if (g_fails == 0) std::cout << GREEN "[ALL PASSED]" RESET "\n";
    else              std::cout << RED "[" << g_fails << " FAILURE(S)]" RESET "\n";
    return g_fails > 0 ? 1 : 0;
}
