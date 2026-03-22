#include "medialaxis/MedialAxis.hxx"
#include <cmath>
#include <stdexcept>



MedialAxis::MedialAxis(std::shared_ptr<Mesh> mesh) : mesh(mesh) {
    if (!mesh) {
        throw std::runtime_error("Mesh pointer cannot be null.");
    }

    int numTriangles = static_cast<int>(mesh->triangles.size());
    medialNodes.resize(numTriangles);

    // Step 1: Create medial nodes from Voronoi vertices
    // Each triangle's circumcenter is a discrete medial axis vertex.
    // The index of the node in medialNodes directly matches the triangle index.
    for (int t = 0; t < numTriangles; ++t) {
        medialNodes[t].id = t;
        medialNodes[t].coord = computeCircumcenter(t);
    }

    // Step 2: Loop through all triangles and extract neigborship information to build the medial edges and compute node degrees.
    for (int t = 0; t < numTriangles; ++t) {
        const auto& adj = mesh->triangleAdjacency[t];
        for (int e = 0; e < 3; ++e) {
            int neighborTri = adj[e];
            if (neighborTri != -1) {
                // This is an internal edge, so we add a medial edge between the current triangle's medial node and the neighbor triangle's medial node.
                medialNodes[t].neighbors.push_back(neighborTri);
            }
        }
     }

    // store the degree of each medial node for convenience (can be derived from neighbors.size())
    for (auto& node : medialNodes) {
        int deg = static_cast<int>(node.neighbors.size());
        if (deg < 1 || deg > 3) {
            throw std::runtime_error("Invalid medial node degree: " + std::to_string(deg));
        }
        node.degree = deg;
     }

    // create the list of medial edges from the neighbor information. We only add an edge (t, neighbor) if t < neighbor to avoid duplicates since the adjacency is symmetric.
    for (int t = 0; t < numTriangles; ++t) {
        const auto& neighbors = medialNodes[t].neighbors;
        for (int neighbor : neighbors) {
            if (t < neighbor) {
                medialEdges.push_back({t, neighbor});
            }
        }
    }

    // Step 3: Build bdyNodeToEdges map for boundary vertices.
    // For each boundary vertex v, find the two boundary edges incident on it.
    // bdyNodeToEdges[v] = {edge0, edge1} where both are boundary edge indices
    // incident on v. We fill slot [0] first, then [1].
    // Note: edges are stored in canonical (min,max) order, so we cannot rely
    // on vertex position within the edge to distinguish the two edges.
    // Pre-initialize entries for all boundary vertices to {-1, -1}.
    for (int bv : mesh->boundaryVertices) {
        bdyNodeToEdges[bv] = {-1, -1};
    }
    for (int beIdx : mesh->boundaryEdges) {
        int v0 = mesh->edges[beIdx][0];
        int v1 = mesh->edges[beIdx][1];

        // For each endpoint, fill the first available slot
        auto& entry0 = bdyNodeToEdges[v0];
        if (entry0[0] == -1) {
            entry0[0] = beIdx;
        } else {
            entry0[1] = beIdx;
        }

        auto& entry1 = bdyNodeToEdges[v1];
        if (entry1[0] == -1) {
            entry1[0] = beIdx;
        } else {
            entry1[1] = beIdx;
        }
    }
}

Point MedialAxis::computeCircumcenter(int triIndex) {
    const Triangle& tri = mesh->triangles[triIndex];
    const Point& a = mesh->vertices[tri[0]];
    const Point& b = mesh->vertices[tri[1]];
    const Point& c = mesh->vertices[tri[2]];

    double ax = a[0], ay = a[1];
    double bx = b[0], by = b[1];
    double cx = c[0], cy = c[1];

    // Denominator for the circumcenter Cartesian coordinate formula
    double d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));

    // Handle degenerate/collinear triangles to prevent division by zero
    if (std::abs(d) < 1e-9) {
        // Fallback to the triangle centroid
        return {(ax + bx + cx) / 3.0, (ay + by + cy) / 3.0};
    }

    double aSq = ax * ax + ay * ay;
    double bSq = bx * bx + by * by;
    double cSq = cx * cx + cy * cy;

    double ux = (aSq * (by - cy) + bSq * (cy - ay) + cSq * (ay - by)) / d;
    double uy = (aSq * (cx - bx) + bSq * (ax - cx) + cSq * (bx - ax)) / d;

    return {ux, uy};
}

void MedialAxis::constructPhiInverseMapping() {
    // for each medial node, we want to find the corresponding point(s) on the boundary that map to it under the phi mapping.
    // first we loop through all medial nodes and identify which ones correspond to boundary triangles (these will have degree 2) and which correspond to corner triangles (degree 1).
    int n = 0;
    EdgeToMedialNode.resize(mesh->edges.size(), {-1, -1});
    for (auto& node : medialNodes) {
        int triIndex = node.id;

        // for each edge of this triangle, check if it's a boundary edge. If so, we add the corresponding BoundaryMappedPoint to the node's preImage.
        int ld = 0;
        for (int e = 0; e < 3; ++e) {
            int edgeIdx = mesh->triangleEdges[triIndex][e];
            if (mesh->isBoundaryEdge[edgeIdx]) {
                // This edge is a boundary edge, so we create a BoundaryMappedPoint for it.
                BoundaryMappedPoint bmp;
                bmp.edgeIndex = edgeIdx;
                bmp.t = 0.5; // Map to the midpoint for now
                bmp.lid = e; // local edge index in the triangle
                node.preImage.push_back(bmp);
                EdgeToMedialNode[edgeIdx][0] = n;
                EdgeToMedialNode[edgeIdx][1] = ld;
                ld++;
            }
        }
        n++;
    }

    n = 0;
    for (auto& node : medialNodes) {
        int triIndex = node.id;
        int deg = node.degree;
        if (deg == 2) {
            int vlid = edgeToOppositeVertex[node.preImage[0].lid]; // local vertex index opposite the first boundary edge
            int vIdx = mesh->triangles[triIndex][vlid]; // global vertex
            // get the two boundary edges connected to this vertex from bdyNodeToEdges
            int be0 = bdyNodeToEdges[vIdx][0];
            int be1 = bdyNodeToEdges[vIdx][1];
            // find the "other" vertex of each boundary edge (the one that isn't vIdx)
            int v0 = (mesh->edges[be0][0] == vIdx) ? mesh->edges[be0][1] : mesh->edges[be0][0];
            int v1 = (mesh->edges[be1][0] == vIdx) ? mesh->edges[be1][1] : mesh->edges[be1][0];

            // check if 3 points v0, vIdx, v1 are collinear
            if (areCollinear(mesh->vertices[v0], mesh->vertices[vIdx], mesh->vertices[v1])) {
                // If they are colinear we simply project the medial node onto the line defined by v0 and v1. This is a simple fallback to handle degenerate cases where the medial node lies on a straight boundary segment.
                auto [t0, valid0] = computePointEdgeProjection(node.coord, mesh->vertices[mesh->edges[be0][0]], mesh->vertices[mesh->edges[be0][1]]);
                auto [t1, valid1] = computePointEdgeProjection(node.coord, mesh->vertices[mesh->edges[be1][0]], mesh->vertices[mesh->edges[be1][1]]);
                if (valid0) {
                    BoundaryMappedPoint bmp;
                    bmp.edgeIndex = be0;
                    bmp.t = t0;
                    node.preImage.push_back(bmp);
                }
                if (valid1) {
                    BoundaryMappedPoint bmp;
                    bmp.edgeIndex = be1;
                    bmp.t = t1;
                    node.preImage.push_back(bmp);
                }

            } else {
                // if they are not colinear we compute the intersection between the lines coming from edges be0 and be1.
                Point p0 = medialNodes[EdgeToMedialNode[be0][0]].coord;
                Point p1 = (mesh->vertices[v0] + mesh->vertices[vIdx]) * 0.5;
                Point p3 = medialNodes[EdgeToMedialNode[be1][0]].coord;
                Point p4 = (mesh->vertices[v1] + mesh->vertices[vIdx]) * 0.5;

                Point c = computeLineIntersection(p0, p1, p3, p4);

                // Ray direction from c through node.coord
                double rdx = node.coord[0] - c[0];
                double rdy = node.coord[1] - c[1];
                double rayLenSq = rdx * rdx + rdy * rdy;

                if (rayLenSq < 1e-20) {
                    // Degenerate case: c ≈ node.coord (e.g. near-circular geometry
                    // where all circumcenters collapse to the same point).
                    // Fall back to orthogonal projection.
                    auto [t0, valid0] = computePointEdgeProjection(node.coord, mesh->vertices[mesh->edges[be0][0]], mesh->vertices[mesh->edges[be0][1]]);
                    auto [t1, valid1] = computePointEdgeProjection(node.coord, mesh->vertices[mesh->edges[be1][0]], mesh->vertices[mesh->edges[be1][1]]);
                    if (valid0) {
                        BoundaryMappedPoint bmp;
                        bmp.edgeIndex = be0;
                        bmp.t = t0;
                        node.preImage.push_back(bmp);
                    }
                    if (valid1) {
                        BoundaryMappedPoint bmp;
                        bmp.edgeIndex = be1;
                        bmp.t = t1;
                        node.preImage.push_back(bmp);
                    }
                } else {
                    // try intersection for both edges be0 and be1
                    auto [t0, valid0] = computeLineEdgeIntersection(c, node.coord, mesh->vertices[mesh->edges[be0][0]], mesh->vertices[mesh->edges[be0][1]]);
                    auto [t1, valid1] = computeLineEdgeIntersection(c, node.coord, mesh->vertices[mesh->edges[be1][0]], mesh->vertices[mesh->edges[be1][1]]);

                    if (valid0) {
                        BoundaryMappedPoint bmp;
                        bmp.edgeIndex = be0;
                        bmp.t = t0;
                        node.preImage.push_back(bmp);
                    }
                    if (valid1) {
                        BoundaryMappedPoint bmp;
                        bmp.edgeIndex = be1;
                        bmp.t = t1;
                        node.preImage.push_back(bmp);
                    }
                }
            }
        }
    }  
}

bool MedialAxis::areCollinear(const Point& a, const Point& b, const Point& c) {
    // Compute the area of the triangle formed by points a, b, c using the determinant method
    double area = 0.5 * std::abs(a[0] * (b[1] - c[1]) + 
                                  b[0] * (c[1] - a[1]) + 
                                  c[0] * (a[1] - b[1]));
    return area < 1e-9; // Consider collinear if area is very small
}

Point MedialAxis::computeLineIntersection(const Point& p1, const Point& p2, const Point& p3, const Point& p4) {
    // Line 1: p1 -> p2,  direction d1 = p2 - p1
    // Line 2: p3 -> p4,  direction d2 = p4 - p3
    //
    // Parametric form:  L1(s) = p1 + s * d1
    //                   L2(t) = p3 + t * d2
    //
    // Solve:  p1 + s * d1 = p3 + t * d2
    //   =>    s * d1 - t * d2 = p3 - p1
    //
    // In 2D this is the 2×2 system:
    //   | d1x  -d2x | | s |   | p3x - p1x |
    //   | d1y  -d2y | | t | = | p3y - p1y |
    //
    // Cramer's rule:  det = d1x * (-d2y) - (-d2x) * d1y
    //                     = d2x * d1y - d1x * d2y
    //                     = cross2(d2, d1)

    double d1x = p2[0] - p1[0];
    double d1y = p2[1] - p1[1];
    double d2x = p4[0] - p3[0];
    double d2y = p4[1] - p3[1];

    double denom = d2x * d1y - d1x * d2y; // cross2(d2, d1)

    // Near-parallel guard: fall back to midpoint of the two reference points
    if (std::abs(denom) < 1e-14) {
        return {0.5 * (p1[0] + p3[0]), 0.5 * (p1[1] + p3[1])};
    }

    double rx = p3[0] - p1[0];
    double ry = p3[1] - p1[1];

    // s = (rx * (-d2y) - (-d2x) * ry) / denom = (d2x * ry - rx * d2y) / denom
    double s = (d2x * ry - rx * d2y) / denom;

    return {p1[0] + s * d1x, p1[1] + s * d1y};
}

std::pair<double, bool> MedialAxis::computeLineEdgeIntersection(const Point& p1, const Point& p2, const Point& e0, const Point& e1) {
    // Line:  L(s) = p1 + s * d1,        d1 = p2 - p1   (unbounded)
    // Edge:  E(t) = e0 + t * d2,         d2 = e1 - e0   (t ∈ [0, 1])
    //
    // Solve the same 2×2 system as computeLineIntersection but also
    // recover t to check whether the hit lies on the segment.

    double d1x = p2[0] - p1[0];
    double d1y = p2[1] - p1[1];
    double d2x = e1[0] - e0[0];
    double d2y = e1[1] - e0[1];

    double denom = d2x * d1y - d1x * d2y; // cross2(d2, d1)

    // Near-parallel: no valid intersection
    if (std::abs(denom) < 1e-14) {
        return {-1.0, false};
    }

    double rx = e0[0] - p1[0];
    double ry = e0[1] - p1[1];

    // s for the line, t for the edge
    double s = (d2x * ry - rx * d2y) / denom;
    double t = (d1x * ry - rx * d1y) / denom;

    // Clamp t to [0, 1] — intersection must lie on the edge segment
    if (t < 0.0 || t > 1.0) {
        return {-1.0, false};
    }

    return {t, true};
}

std::pair<double, bool> MedialAxis::computePointEdgeProjection(const Point& p, const Point& e0, const Point& e1) {
    // Project point p onto the line through e0 → e1.
    // E(t) = e0 + t * d,  d = e1 - e0
    // The closest point on the line satisfies:  dot(p - E(t), d) = 0
    //   dot(p - e0 - t*d, d) = 0
    //   dot(p - e0, d) = t * dot(d, d)
    //   t = dot(p - e0, d) / dot(d, d)

    double dx = e1[0] - e0[0];
    double dy = e1[1] - e0[1];
    double lenSq = dx * dx + dy * dy;

    if (lenSq < 1e-28) {
        // Degenerate edge (zero length)
        return {-1.0, false};
    }

    double t = ((p[0] - e0[0]) * dx + (p[1] - e0[1]) * dy) / lenSq;

    if (t < 0.0 || t > 1.0) {
        return {t, false};
    }

    return {t, true};
}