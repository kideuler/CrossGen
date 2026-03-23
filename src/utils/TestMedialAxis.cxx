// Utility to construct the medial axis from an .obj mesh and debug constructPhiInverseMapping.
// Usage: TestMedialAxis <mesh.obj>

#include <iostream>
#include <iomanip>
#include <memory>
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
    std::cout << "  Medial nodes:   " << ma.medialNodes.size() << "\n";
    std::cout << "  Medial edges:   " << ma.medialEdges.size() << "\n";

    // ── Count degrees ──
    int nDeg1 = 0, nDeg2 = 0, nDeg3 = 0;
    for (const auto& node : ma.medialNodes) {
        if (node.degree == 1) ++nDeg1;
        else if (node.degree == 2) ++nDeg2;
        else if (node.degree == 3) ++nDeg3;
    }
    std::cout << "  Degree counts:  deg1=" << nDeg1 << "  deg2=" << nDeg2 << "  deg3=" << nDeg3 << "\n";

    // ── Run constructPhiInverseMapping ──
    ma.constructMappingPhase1();
    ma.constructMappingPhase2();

    // ── Sanity checks ──
    int totalFails = 0;

    // Check 1: Every degree-1 node should have exactly 2 preImage entry (from the initial pass)
    {
        int pass = 0, fail = 0;
        for (const auto& node : ma.medialNodes) {
            if (node.degree != 1) continue;
            int expected = 2; // corner triangle: 2 boundary edges
            if ((int)node.preImage.size() != expected) {
                if (fail < 10) {
                    std::cout << YELLOW "  [deg1] node " << node.id
                              << ": expected " << expected << " preImage entries, got "
                              << node.preImage.size() << RESET "\n";
                }
                ++fail;
            } else {
                ++pass;
            }
        }
        std::cout << (fail == 0 ? GREEN "[PASS]" RESET : RED "[FAIL]" RESET)
                  << " Degree-1 preImage count: " << pass << " pass, " << fail << " fail\n";
        totalFails += fail;
    }

    // Check 2: Every degree-2 node should have exactly 2 preImage entries
    //   (1 from the boundary edge of the triangle + 1 from the adjacent boundary edges through intersection)
    {
        int pass = 0, fail = 0;
        for (const auto& node : ma.medialNodes) {
            if (node.degree != 2) continue;
            int expected = 2;
            if ((int)node.preImage.size() != expected) {
                if (fail < 10) {
                    std::cout << YELLOW "  [deg2] node " << node.id
                              << ": expected " << expected << " preImage entries, got "
                              << node.preImage.size() << "\n";

                    // Detailed debug info for this node
                    int triIndex = node.id;
                    std::cout << "    triangle vertices:";
                    for (int k = 0; k < 3; ++k) {
                        int vi = delaunayMesh->triangles[triIndex][k];
                        std::cout << " v" << vi << "(" << delaunayMesh->vertices[vi][0]
                                  << ", " << delaunayMesh->vertices[vi][1] << ")";
                    }
                    std::cout << "\n";
                    std::cout << "    node coord: (" << node.coord[0] << ", " << node.coord[1] << ")\n";

                    // Print which boundary edges this triangle has
                    std::cout << "    triangle edges:";
                    for (int e = 0; e < 3; ++e) {
                        int edgeIdx = delaunayMesh->triangleEdges[triIndex][e];
                        bool isBdy = delaunayMesh->isBoundaryEdge[edgeIdx];
                        std::cout << " e" << edgeIdx << (isBdy ? "(bdy)" : "(int)");
                    }
                    std::cout << "\n";

                    // Print the preImage entries we do have
                    for (int pi = 0; pi < (int)node.preImage.size(); ++pi) {
                        const auto& bmp = node.preImage[pi];
                        Point pt = bmp.evaluate(delaunayMesh);
                        std::cout << "    preImage[" << pi << "]: edge=" << bmp.edgeIndex
                                  << " lid=" << (int)bmp.lid << " t=" << bmp.t
                                  << " -> (" << pt[0] << ", " << pt[1] << ")\n";
                    }

                    // Print the vertex opposite the boundary edge and its adjacent bdy edges
                    if (!node.preImage.empty()) {
                        int vlid = edgeToOppositeVertex[node.preImage[0].lid];
                        int vIdx = delaunayMesh->triangles[triIndex][vlid];
                        std::cout << "    opposite vertex: v" << vIdx
                                  << " (" << delaunayMesh->vertices[vIdx][0]
                                  << ", " << delaunayMesh->vertices[vIdx][1] << ")\n";

                        if (ma.bdyNodeToEdges.count(vIdx)) {
                            int be0 = ma.bdyNodeToEdges[vIdx][0];
                            int be1 = ma.bdyNodeToEdges[vIdx][1];
                            int v0 = delaunayMesh->edges[be0][0];
                            int v1 = delaunayMesh->edges[be1][1];
                            std::cout << "    be0=" << be0 << " (v" << delaunayMesh->edges[be0][0]
                                      << "->v" << delaunayMesh->edges[be0][1] << ")\n";
                            std::cout << "    be1=" << be1 << " (v" << delaunayMesh->edges[be1][0]
                                      << "->v" << delaunayMesh->edges[be1][1] << ")\n";

                            // Compute and print the intersection point c
                            int mn0 = ma.EdgeToMedialNode[be0][0];
                            int mn1 = ma.EdgeToMedialNode[be1][0];
                            if (mn0 >= 0 && mn1 >= 0) {
                                Point p0 = ma.medialNodes[mn0].coord;
                                Point p1 = (delaunayMesh->vertices[v0] + delaunayMesh->vertices[vIdx]) * 0.5;
                                Point p3 = ma.medialNodes[mn1].coord;
                                Point p4 = (delaunayMesh->vertices[v1] + delaunayMesh->vertices[vIdx]) * 0.5;

                                std::cout << "    line0: (" << p0[0] << "," << p0[1] << ") -> ("
                                          << p1[0] << "," << p1[1] << ")\n";
                                std::cout << "    line1: (" << p3[0] << "," << p3[1] << ") -> ("
                                          << p4[0] << "," << p4[1] << ")\n";

                                // Ray from c through node.coord
                                // Compute c as the intersection of the two lines
                                // (same logic as in constructPhiInverseMapping)
                                // We need to replicate it here since it's private
                                // Instead, print what the ray direction would be
                                std::cout << "    (use debug output above to check ray/edge intersection)\n";
                            } else {
                                std::cout << "    EdgeToMedialNode not set for be0 or be1!\n";
                                if (mn0 < 0) std::cout << "      be0 (edge " << be0 << ") has no medial node mapping\n";
                                if (mn1 < 0) std::cout << "      be1 (edge " << be1 << ") has no medial node mapping\n";
                            }
                        } else {
                            std::cout << "    vertex v" << vIdx << " not in bdyNodeToEdges!\n";
                        }
                    }
                }
                std::cout << RESET;
                ++fail;
            } else {
                ++pass;
            }
        }
        std::cout << (fail == 0 ? GREEN "[PASS]" RESET : RED "[FAIL]" RESET)
                  << " Degree-2 preImage count: " << pass << " pass, " << fail << " fail\n";
        totalFails += fail;
    }

    // Check 3: Every degree-3 node should have exactly 0 preImage entries (no boundary edge)
    {
        int pass = 0, fail = 0;
        for (const auto& node : ma.medialNodes) {
            if (node.degree != 3) continue;
            if (node.preImage.size() != 3) {
                if (fail < 5) {
                    std::cout << YELLOW "  [deg3] node " << node.id
                              << ": expected 3 preImage entries, got "
                              << node.preImage.size() << RESET "\n";
                }
                ++fail;
            } else {
                ++pass;
            }
        }
        std::cout << (fail == 0 ? GREEN "[PASS]" RESET : RED "[FAIL]" RESET)
                  << " Degree-3 preImage count: " << pass << " pass, " << fail << " fail\n";
        totalFails += fail;
    }

    // Check 4: All preImage t-values should be in [0, 1]
    {
        int pass = 0, fail = 0;
        for (const auto& node : ma.medialNodes) {
            for (const auto& bmp : node.preImage) {
                if (bmp.t < 0.0 || bmp.t > 1.0) {
                    if (fail < 5) {
                        std::cout << YELLOW "  node " << node.id << ": preImage t=" << bmp.t
                                  << " out of [0,1] on edge " << bmp.edgeIndex << RESET "\n";
                    }
                    ++fail;
                } else {
                    ++pass;
                }
            }
        }
        std::cout << (fail == 0 ? GREEN "[PASS]" RESET : RED "[FAIL]" RESET)
                  << " PreImage t-values in [0,1]: " << pass << " pass, " << fail << " fail\n";
        totalFails += fail;
    }

    // Check 5: bdyNodeToEdges should be populated for all boundary vertices
    {
        int pass = 0, fail = 0;
        for (int bv : delaunayMesh->boundaryVertices) {
            if (!ma.bdyNodeToEdges.count(bv)) {
                if (fail < 5) {
                    std::cout << YELLOW "  boundary vertex v" << bv << " missing from bdyNodeToEdges" RESET "\n";
                }
                ++fail;
            } else {
                int be0 = ma.bdyNodeToEdges[bv][0];
                int be1 = ma.bdyNodeToEdges[bv][1];
                if (be0 < 0 || be1 < 0) {
                    if (fail < 5) {
                        std::cout << YELLOW "  boundary vertex v" << bv
                                  << ": bdyNodeToEdges has invalid entry ("
                                  << be0 << ", " << be1 << ")" RESET "\n";
                    }
                    ++fail;
                } else {
                    ++pass;
                }
            }
        }
        std::cout << (fail == 0 ? GREEN "[PASS]" RESET : RED "[FAIL]" RESET)
                  << " bdyNodeToEdges completeness: " << pass << " pass, " << fail << " fail\n";
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
