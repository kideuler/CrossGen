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
        medialVertices.push_back(computeCircumcenter(t));
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