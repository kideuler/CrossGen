#include "medialaxis/MedialAxis.hxx"
#include <cmath>
#include <limits>
#include <stdexcept>

MedialAxis::MedialAxis(std::shared_ptr<Mesh> mesh) : mesh(mesh) {
    if (!mesh) {
        throw std::runtime_error("Mesh pointer cannot be null.");
    }

    int numTriangles = static_cast<int>(mesh->triangles.size());
    medialVertices.reserve(numTriangles);

    // Step 1: Calculate Voronoi vertices
    // Each triangle's circumcenter is a discrete medial axis vertex.
    // The index of the circumcenter in medialVertices directly matches the triangle index.
    for (int t = 0; t < numTriangles; ++t) {
        Point cc = computeCircumcenter(t);
        const Triangle& tri = mesh->triangles[t];

        // Compute radius as distance from circumcenter to any triangle vertex
        const Point& v0 = mesh->vertices[tri[0]];
        double radius = normP(cc - v0);

        // Collect the touch points (triangle vertices equidistant to the circumcenter)
        std::unordered_set<int> touchPoints;
        touchPoints.insert(tri[0]);
        touchPoints.insert(tri[1]);
        touchPoints.insert(tri[2]);

        MedialVertex mv;
        mv.coord = cc;
        mv.triangleIndex = t;
        mv.touchPoints = std::move(touchPoints);
        mv.degree = 0;       // will be set after edges are built
        mv.radius = radius;
        mv.active = true;
        mv.nodeType = TopMakerNodeType::Normal;

        medialVertices.push_back(std::move(mv));
    }

    // Step 2: Extract internal Voronoi edges
    // The dual of any internal Delaunay edge is a medial axis edge 
    // connecting the circumcenters of the two adjacent triangles.
    for (const auto& edgeTri : mesh->edgeTriangles) {
        int t0 = edgeTri[0];
        int t1 = edgeTri[1];

        // An edge is strictly internal if it is shared by two valid triangles 
        // (boundary edges have -1 for one of the adjacent triangles).
        if (t0 != -1 && t1 != -1) {
            medialEdges.push_back({t0, t1});
            // Update neighbor sets and degree for both endpoints
            medialVertices[t0].neighbors.insert(t1);
            medialVertices[t1].neighbors.insert(t0);
        }
    }

    // Step 3: Set the degree of each medial vertex from its neighbor set
    for (auto& mv : medialVertices) {
        mv.degree = static_cast<int>(mv.neighbors.size());
    }

    // Step 4: Build boundary vertex connectivity (ordered boundary loops)
    // For each boundary edge, determine the oriented direction using the
    // adjacent triangle's winding, then build next/prev maps directly.
    int numVertices = static_cast<int>(mesh->vertices.size());

    // Initialize boundaryVertices array indexed by vertex index
    boundaryVertices.resize(numVertices);
    for (int i = 0; i < numVertices; ++i) {
        boundaryVertices[i].vertexIndex = i;
        boundaryVertices[i].nextBoundaryVertex = -1;
        boundaryVertices[i].prevBoundaryVertex = -1;
    }

    // For each boundary edge, use the adjacent triangle's winding to determine
    // the oriented boundary direction: orientedFrom -> orientedTo.
    // This directly gives us next[orientedFrom] = orientedTo for that edge.
    // For each boundary edge, use the adjacent triangle's winding to determine
    // the oriented boundary direction: orientedFrom -> orientedTo.
    // This directly gives us next[orientedFrom] = orientedTo for that edge.
    for (int eIdx : mesh->boundaryEdges) {
        // Find the (single) adjacent triangle
        int triIdx = (mesh->edgeTriangles[eIdx][0] != -1)
                         ? mesh->edgeTriangles[eIdx][0]
                         : mesh->edgeTriangles[eIdx][1];
        const Triangle& tri = mesh->triangles[triIdx];

        // Find which local edge of the triangle this is
        int localEdge = -1;
        for (int le = 0; le < 3; ++le) {
            if (mesh->triangleEdges[triIdx][le] == eIdx) {
                localEdge = le;
                break;
            }
        }

        // Local edge le in CCW winding: tri[le] -> tri[(le+1)%3].
        // The boundary runs in the opposite direction so that the mesh
        // interior stays on the left: orientedFrom -> orientedTo.
        int orientedFrom = tri[(localEdge + 1) % 3];
        int orientedTo   = tri[localEdge];

        boundaryVertices[orientedFrom].nextBoundaryVertex = orientedTo;
        boundaryVertices[orientedTo].prevBoundaryVertex = orientedFrom;
    }

    // Step 5: Initialize sharpVertices based on the angle at each boundary vertex.
    // A boundary vertex is "sharp" if the angle between its two adjacent boundary
    // edges is below a threshold (i.e. the boundary makes a sharp turn).
    const double sharpAngleThreshold = M_PI * 0.75; // 135 degrees
    sharpVertices.resize(numVertices, false);
    for (int bv : mesh->boundaryVertices) {
        int prev = boundaryVertices[bv].prevBoundaryVertex;
        int next = boundaryVertices[bv].nextBoundaryVertex;
        if (prev == -1 || next == -1) continue;

        Point eToPrev = mesh->vertices[prev] - mesh->vertices[bv];
        Point eToNext = mesh->vertices[next] - mesh->vertices[bv];

        double lenPrev = normP(eToPrev);
        double lenNext = normP(eToNext);
        if (lenPrev < 1e-12 || lenNext < 1e-12) continue;

        double cosAngle = dotP(eToPrev, eToNext) / (lenPrev * lenNext);
        cosAngle = std::max(-1.0, std::min(1.0, cosAngle)); // clamp for numerical safety
        double angle = std::acos(cosAngle);

        if (angle < sharpAngleThreshold) {
            sharpVertices[bv] = true;
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

void MedialAxis::deduplicateMedialVertices() {
    
    // Pass 1: Deduplicate degree-1 medial vertices.
    // Loop through degree 1 medial vertices and check if the connected vertex
    // is close within tolerance with respect to radius. If so, absorb the
    // degree-1 vertex into its neighbor: transfer touch points, remove the
    // leaf, and update the neighbor's degree.

    bool done = false;
    while (!done) {
        done = true;

        for (size_t i = 0; i < medialVertices.size(); ++i) {
            MedialVertex& mv = medialVertices[i];
            if (!mv.active || mv.degree != 1) continue;

            // Get the single neighbor
            int neighborIdx = *mv.neighbors.begin();
            MedialVertex& neighbor = medialVertices[neighborIdx];

            // Check distance between mv and neighbor relative to radius
            double dist = normP(mv.coord - neighbor.coord);
            double localRadius = std::max(mv.radius, neighbor.radius);
            if (dist < MEDIAL_VERTEX_MERGE_TOLERANCE * localRadius) {
                // Transfer touch points to the neighbor
                neighbor.touchPoints.insert(mv.touchPoints.begin(), mv.touchPoints.end());

                // Remove mv from neighbor's neighbor set
                neighbor.neighbors.erase(static_cast<int>(i));

                // Mark mv as inactive and clear its connections
                mv.active = false;
                mv.neighbors.clear();
                mv.degree = 0;

                // Update neighbor's degree
                neighbor.degree = static_cast<int>(neighbor.neighbors.size());

                done = false;
                break;
            }
        }
    }

    // Pass 2: Deduplicate degree-2 medial vertices.
    // A degree-2 vertex sits in the middle of a chain. If it is close to one
    // of its two neighbors, absorb it into that neighbor and reconnect the
    // other neighbor to the merge target.

    done = false;
    while (!done) {
        done = true;

        for (size_t i = 0; i < medialVertices.size(); ++i) {
            MedialVertex& mv = medialVertices[i];
            if (!mv.active || mv.degree != 2) continue;

            // Find the closest neighbor within tolerance
            int mergeTarget = -1;
            double bestDist = std::numeric_limits<double>::max();
            for (int nb : mv.neighbors) {
                double dist = normP(mv.coord - medialVertices[nb].coord);
                double localRadius = std::max(mv.radius, medialVertices[nb].radius);
                if (dist < MEDIAL_VERTEX_MERGE_TOLERANCE * localRadius && dist < bestDist) {
                    bestDist = dist;
                    mergeTarget = nb;
                }
            }
            if (mergeTarget == -1) continue;

            MedialVertex& target = medialVertices[mergeTarget];

            // Transfer touch points
            target.touchPoints.insert(mv.touchPoints.begin(), mv.touchPoints.end());

            // Reconnect all other neighbors of mv to the merge target
            for (int nb : mv.neighbors) {
                if (nb == mergeTarget) continue;
                MedialVertex& other = medialVertices[nb];
                other.neighbors.erase(static_cast<int>(i));
                other.neighbors.insert(mergeTarget);
                target.neighbors.insert(nb);
            }

            // Remove mv from the merge target's neighbor set
            target.neighbors.erase(static_cast<int>(i));

            // Deactivate mv
            mv.active = false;
            mv.neighbors.clear();
            mv.degree = 0;

            // Update target's degree
            target.degree = static_cast<int>(target.neighbors.size());

            // Update degrees of reconnected neighbors
            for (int nb : target.neighbors) {
                medialVertices[nb].degree = static_cast<int>(medialVertices[nb].neighbors.size());
            }

            done = false;
            break;
        }
    }

    // Pass 3: Deduplicate degree-3+ medial vertices.
    // A high-degree vertex is a junction. If it is close to one of its
    // neighbors, absorb it into that neighbor and rewire all other neighbors.

    done = false;
    while (!done) {
        done = true;

        for (size_t i = 0; i < medialVertices.size(); ++i) {
            MedialVertex& mv = medialVertices[i];
            if (!mv.active || mv.degree < 3) continue;

            // Find the closest neighbor within tolerance
            int mergeTarget = -1;
            double bestDist = std::numeric_limits<double>::max();
            for (int nb : mv.neighbors) {
                double dist = normP(mv.coord - medialVertices[nb].coord);
                double localRadius = std::max(mv.radius, medialVertices[nb].radius);
                if (dist < MEDIAL_VERTEX_MERGE_TOLERANCE * localRadius && dist < bestDist) {
                    bestDist = dist;
                    mergeTarget = nb;
                }
            }
            if (mergeTarget == -1) continue;

            MedialVertex& target = medialVertices[mergeTarget];

            // Transfer touch points
            target.touchPoints.insert(mv.touchPoints.begin(), mv.touchPoints.end());

            // Reconnect all other neighbors of mv to the merge target
            for (int nb : mv.neighbors) {
                if (nb == mergeTarget) continue;
                MedialVertex& other = medialVertices[nb];
                other.neighbors.erase(static_cast<int>(i));
                other.neighbors.insert(mergeTarget);
                target.neighbors.insert(nb);
            }

            // Remove mv from the merge target's neighbor set
            target.neighbors.erase(static_cast<int>(i));

            // Deactivate mv
            mv.active = false;
            mv.neighbors.clear();
            mv.degree = 0;

            // Update target's degree
            target.degree = static_cast<int>(target.neighbors.size());

            // Update degrees of reconnected neighbors
            for (int nb : target.neighbors) {
                medialVertices[nb].degree = static_cast<int>(medialVertices[nb].neighbors.size());
            }

            done = false;
            break;
        }
    }

    // ── Compaction: remove inactive vertices and rebuild edges with new indices ──

    // Build old-to-new index mapping for active vertices only
    std::vector<int> oldToNew(medialVertices.size(), -1);
    int newIdx = 0;
    for (size_t i = 0; i < medialVertices.size(); ++i) {
        if (medialVertices[i].active) {
            oldToNew[i] = newIdx++;
        }
    }

    // Build compacted vertex list with remapped neighbor sets
    std::vector<MedialVertex> compacted;
    compacted.reserve(newIdx);
    for (size_t i = 0; i < medialVertices.size(); ++i) {
        if (!medialVertices[i].active) continue;
        MedialVertex mv = std::move(medialVertices[i]);

        // Remap neighbor indices
        std::unordered_set<int> newNeighbors;
        for (int nb : mv.neighbors) {
            int mapped = oldToNew[nb];
            if (mapped != -1) {
                newNeighbors.insert(mapped);
            }
        }
        mv.neighbors = std::move(newNeighbors);
        mv.degree = static_cast<int>(mv.neighbors.size());

        compacted.push_back(std::move(mv));
    }

    // Rebuild medial edges from neighbor sets (each edge stored once, with i < j)
    std::vector<Edge> newEdges;
    for (int i = 0; i < static_cast<int>(compacted.size()); ++i) {
        for (int nb : compacted[i].neighbors) {
            if (i < nb) {
                newEdges.push_back({i, nb});
            }
        }
    }

    medialVertices = std::move(compacted);
    medialEdges = std::move(newEdges);
}