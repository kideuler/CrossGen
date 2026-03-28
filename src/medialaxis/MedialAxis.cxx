#include "medialaxis/MedialAxis.hxx"
#include <cmath>
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