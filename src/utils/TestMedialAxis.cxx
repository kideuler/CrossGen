// Utility to construct the medial axis from an .obj mesh and debug constructPhiInverseMapping.
// Usage: TestMedialAxis <mesh.obj>

#include <iostream>
#include <iomanip>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
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

    // ── Delaunay re-triangulation ──
    auto delaunayMesh = buildDelaunayMesh(originalMesh);
    if (!delaunayMesh) {
        std::cerr << RED "[FAIL]" RESET " Could not build Delaunay mesh.\n";
        return 1;
    }
    std::cout << "  Delaunay mesh:  " << delaunayMesh->triangles.size() << " triangles, "
              << delaunayMesh->vertices.size() << " vertices, "
              << delaunayMesh->boundaryEdges.size() << " boundary edges\n";

    // ── Construct medial axis ──
    MedialAxis ma(delaunayMesh);
    std::cout << "  Medial vertices: " << ma.medialVertices.size() << "\n";
    std::cout << "  Medial edges:    " << ma.medialEdges.size() << "\n";

    // ── Sanity checks ──
    int totalFails = 0;

    // Check 1: There should be one medial vertex per triangle
    {
        bool ok = (ma.medialVertices.size() == delaunayMesh->triangles.size());
        std::cout << (ok ? GREEN "[PASS]" RESET : RED "[FAIL]" RESET)
                  << " medialVertices.size() == triangles.size() ("
                  << ma.medialVertices.size() << " vs " << delaunayMesh->triangles.size() << ")\n";
        if (!ok) ++totalFails;
    }

    // Check 2: All medial edges should reference valid vertex indices
    {
        int pass = 0, fail = 0;
        int nv = static_cast<int>(ma.medialVertices.size());
        for (const auto& edge : ma.medialEdges) {
            if (edge[0] < 0 || edge[0] >= nv || edge[1] < 0 || edge[1] >= nv) {
                if (fail < 5) {
                    std::cout << YELLOW "  medial edge (" << edge[0] << ", " << edge[1]
                              << ") out of range [0, " << nv << ")" RESET "\n";
                }
                ++fail;
            } else {
                ++pass;
            }
        }
        std::cout << (fail == 0 ? GREEN "[PASS]" RESET : RED "[FAIL]" RESET)
                  << " Medial edge index validity: " << pass << " pass, " << fail << " fail\n";
        totalFails += fail;
    }

    // Check 3: Medial edges should only come from internal mesh edges (no boundary duals)
    {
        int numInternal = 0;
        for (const auto& et : delaunayMesh->edgeTriangles) {
            if (et[0] != -1 && et[1] != -1) ++numInternal;
        }
        bool ok = (static_cast<int>(ma.medialEdges.size()) == numInternal);
        std::cout << (ok ? GREEN "[PASS]" RESET : RED "[FAIL]" RESET)
                  << " medialEdges.size() == internal edges ("
                  << ma.medialEdges.size() << " vs " << numInternal << ")\n";
        if (!ok) ++totalFails;
    }

    // Check 4: Medial vertices (circumcenters) should be finite
    {
        int pass = 0, fail = 0;
        for (size_t i = 0; i < ma.medialVertices.size(); ++i) {
            const auto& v = ma.medialVertices[i];
            if (!std::isfinite(v.coord[0]) || !std::isfinite(v.coord[1])) {
                if (fail < 5) {
                    std::cout << YELLOW "  medial vertex " << i
                              << " is not finite: (" << v.coord[0] << ", " << v.coord[1] << ")" RESET "\n";
                }
                ++fail;
            } else {
                ++pass;
            }
        }
        std::cout << (fail == 0 ? GREEN "[PASS]" RESET : RED "[FAIL]" RESET)
                  << " Medial vertex finiteness: " << pass << " pass, " << fail << " fail\n";
        totalFails += fail;
    }

    // ── Run deduplication ──
    ma.deduplicateMedialVertices();
    {
        int activeCount = 0;
        for (const auto& mv : ma.medialVertices) {
            if (mv.active) ++activeCount;
        }
        std::cout << "  After deduplication: " << activeCount << " active medial vertices\n";
    }

    // Check 5: After deduplication, no active degree-1 vertex should be within
    // 0.01 of any other active vertex.
    {
        const double minSepDist = 1e-8;
        int fail = 0;
        for (size_t i = 0; i < ma.medialVertices.size(); ++i) {
            const auto& mv = ma.medialVertices[i];
            if (!mv.active) continue;

            for (size_t j = 0; j < ma.medialVertices.size(); ++j) {
                if (j == i) continue;
                const auto& other = ma.medialVertices[j];
                if (!other.active) continue;

                double dist = normP(mv.coord - other.coord);
                if (dist < minSepDist) {
                    if (fail < 5) {
                        std::cout << YELLOW "  vertex " << i
                                  << " (deg=" << mv.degree << ") is " << dist
                                  << " from vertex " << j
                                  << " (deg=" << other.degree << ")" RESET "\n";
                        // Print neighbors of both
                        std::cout << "    neighbors of " << i << ":";
                        for (int nb : mv.neighbors) std::cout << " " << nb;
                        std::cout << "\n";
                        std::cout << "    neighbors of " << j << ":";
                        for (int nb : other.neighbors) std::cout << " " << nb;
                        std::cout << "\n";
                        // Check if they share a neighbor
                        bool isNeighbor = mv.neighbors.count(static_cast<int>(j)) > 0;
                        std::cout << "    directly connected: " << (isNeighbor ? "YES" : "NO") << "\n";
                        std::cout << "    radius_i=" << mv.radius << " radius_j=" << other.radius
                                  << " tol_i=" << MEDIAL_VERTEX_MERGE_TOLERANCE * mv.radius
                                  << " tol_j=" << MEDIAL_VERTEX_MERGE_TOLERANCE * other.radius << "\n";
                        // Walk path from i to j via neighbors (BFS, max 10 hops)
                        {
                            std::unordered_map<int, int> parent;
                            std::vector<int> queue = {static_cast<int>(i)};
                            parent[static_cast<int>(i)] = -1;
                            bool found = false;
                            for (size_t qi = 0; qi < queue.size() && qi < 1000; ++qi) {
                                int cur = queue[qi];
                                if (cur == static_cast<int>(j)) { found = true; break; }
                                for (int nb : ma.medialVertices[cur].neighbors) {
                                    if (!ma.medialVertices[nb].active) continue;
                                    if (parent.count(nb)) continue;
                                    parent[nb] = cur;
                                    queue.push_back(nb);
                                }
                            }
                            if (found) {
                                std::vector<int> path;
                                for (int c = static_cast<int>(j); c != -1; c = parent[c])
                                    path.push_back(c);
                                std::reverse(path.begin(), path.end());
                                std::cout << "    path (" << path.size() - 1 << " hops):";
                                for (int p : path) std::cout << " " << p;
                                std::cout << "\n";
                            } else {
                                std::cout << "    no path found (disconnected components)\n";
                            }
                        }
                    }
                    ++fail;
                    break;
                }
            }
        }
        std::cout << (fail == 0 ? GREEN "[PASS]" RESET : RED "[FAIL]" RESET)
                  << " Post-dedup separation (>= " << minSepDist << "): "
                  << fail << " violations\n";
        totalFails += fail;
    }

    // ── Create polylines ──
    ma.createPolylines();
    std::cout << "  Polylines:      " << ma.polyLines.size() << "\n";

    // Check 6: Every polyline has at least 2 vertices, all indices in bounds, and no duplicates (forwards/backwards)
    {
        int fail = 0;
        std::set<std::vector<int>> uniquePolys;
        for (size_t i = 0; i < ma.polyLines.size(); ++i) {
            const auto& pl = ma.polyLines[i];
            if (pl.size() < 2) {
                std::cout << YELLOW "  polyline " << i << " has <2 vertices\n" RESET;
                ++fail;
            }
            for (int vi : pl) {
                if (vi < 0 || vi >= static_cast<int>(ma.medialVertices.size())) {
                    std::cout << YELLOW "  polyline " << i << " has out-of-bounds vertex " << vi << "\n" RESET;
                    ++fail;
                }
            }
            // Check for duplicate (forwards/backwards)
            std::vector<int> rev = pl;
            std::reverse(rev.begin(), rev.end());
            if (uniquePolys.count(pl) || uniquePolys.count(rev)) {
                std::cout << YELLOW "  polyline " << i << " is a duplicate (forwards/backwards)\n" RESET;
                ++fail;
            }
            uniquePolys.insert(pl);
        }
        std::cout << (fail == 0 ? GREEN "[PASS]" RESET : RED "[FAIL]" RESET)
                  << " Polyline structure: " << ma.polyLines.size() << " polylines, " << fail << " issues\n";
        totalFails += fail;
    }

    // Check 7: All polyline edges exist in the medial axis
    {
        int fail = 0;
        std::unordered_set<uint64_t> edgeSet;
        for (const auto& e : ma.medialEdges) {
            int a = std::min(e[0], e[1]);
            int b = std::max(e[0], e[1]);
            edgeSet.insert((static_cast<uint64_t>(a) << 32) | static_cast<uint32_t>(b));
        }
        for (size_t i = 0; i < ma.polyLines.size(); ++i) {
            const auto& pl = ma.polyLines[i];
            for (size_t j = 1; j < pl.size(); ++j) {
                int a = std::min(pl[j-1], pl[j]);
                int b = std::max(pl[j-1], pl[j]);
                if (a == b) continue;
                uint64_t key = (static_cast<uint64_t>(a) << 32) | static_cast<uint32_t>(b);
                if (!edgeSet.count(key)) {
                    std::cout << YELLOW "  polyline " << i << " has non-medial edge (" << pl[j-1] << ", " << pl[j] << ")\n" RESET;
                    ++fail;
                }
            }
        }
        std::cout << (fail == 0 ? GREEN "[PASS]" RESET : RED "[FAIL]" RESET)
                  << " Polyline edges exist in medial axis: " << fail << " issues\n";
        totalFails += fail;
    }

    // Check 8: All non-degree-2 vertices are endpoints of at least one polyline
    {
        int fail = 0;
        for (size_t i = 0; i < ma.medialVertices.size(); ++i) {
            const auto& mv = ma.medialVertices[i];
            if (!mv.active) continue;
            if (mv.degree == 2) continue;

            // Count how many polylines have this vertex as an endpoint
            int endpointCount = 0;
            for (const auto& pl : ma.polyLines) {
                if (pl.empty()) continue;
                if (pl.front() == static_cast<int>(i) || pl.back() == static_cast<int>(i)) {
                    ++endpointCount;
                }
            }

            if (endpointCount == 0) {
                std::cout << YELLOW "  non-degree-2 vertex " << i << " is not a polyline endpoint\n" RESET;
                ++fail;
                continue;
            }

            // Warn if the number of polylines attached doesn't match the degree
            if (endpointCount != mv.degree) {
                std::cout << YELLOW "  vertex " << i << " has degree " << mv.degree
                          << " but is an endpoint of " << endpointCount << " polylines\n" RESET;
            }
        }
        std::cout << (fail == 0 ? GREEN "[PASS]" RESET : RED "[FAIL]" RESET)
                  << " All non-degree-2 vertices are polyline endpoints: " << fail << " issues\n";
        totalFails += fail;
    }

    // Check 9: All degree-2 vertices are interior to exactly one polyline (except for pure cycles)
    {
        int fail = 0;
        for (size_t i = 0; i < ma.medialVertices.size(); ++i) {
            const auto& mv = ma.medialVertices[i];
            if (!mv.active) continue;
            if (mv.degree != 2) continue;
            int endpointCount = 0;
            for (const auto& pl : ma.polyLines) {
                if (pl.front() == static_cast<int>(i) || pl.back() == static_cast<int>(i)) {
                    ++endpointCount;
                }
            }
            if (endpointCount > 0) continue; // degree-2 vertex is an endpoint (cycle)
            int count = 0;
            for (const auto& pl : ma.polyLines) {
                for (size_t j = 1; j + 1 < pl.size(); ++j) {
                    if (pl[j] == static_cast<int>(i)) ++count;
                }
            }
            if (count != 1) {
                std::cout << YELLOW "  degree-2 vertex " << i << " is interior to " << count << " polylines\n" RESET;
                ++fail;
            }
        }
        std::cout << (fail == 0 ? GREEN "[PASS]" RESET : RED "[FAIL]" RESET)
                  << " All degree-2 vertices are interior to one polyline: " << fail << " issues\n";
        totalFails += fail;
    }

    // ── Summary ──
    std::cout << "\n";
    if (totalFails == 0) {
        std::cout << GREEN "[ALL PASSED]" RESET "\n";
    } else {
        std::cout << RED "[" << totalFails << " FAILURE(S)]" RESET "\n";
    }

    return totalFails > 0 ? 1 : 0;
}
