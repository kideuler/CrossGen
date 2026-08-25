#ifndef __MESH_HXX__
#define __MESH_HXX__

#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cctype>
#include <stdexcept>
#include "VertexTriangleCSR.hxx"


typedef std::array<double, 2> Point;
typedef std::array<int, 3> Triangle;
typedef std::array<int, 2> Edge;

// operators and functions for Point
inline Point operator+(const Point &a, const Point &b) {
    return {a[0] + b[0], a[1] + b[1]};
}

inline Point operator-(const Point &a, const Point &b) {
    return {a[0] - b[0], a[1] - b[1]};
}

inline Point operator*(const Point &a, double s) {
    return {a[0] * s, a[1] * s};
}

inline Point operator/(const Point &a, double s) {
    return {a[0] / s, a[1] / s};
}

inline double dotP(const Point &a, const Point &b) {
    return a[0] * b[0] + a[1] * b[1];
}

inline double cross2(const Point &a, const Point &b) {
    return a[0] * b[1] - a[1] * b[0];
}

inline double normP(const Point &a) {
    return std::sqrt(dotP(a, a));
}

inline Point normalizeP(const Point &a) {
    double n = normP(a);
    if (n <= 0.0) return {0.0, 0.0};
    return a / n;
}

// compute angle from vector (Point)
inline double computeAngle(const Point &p) {
    return std::atan2(p[1], p[0]);
}

// simple wrap to (-pi, pi]
inline double wrap_pi(double a) {
  a = std::fmod(a + M_PI, 2.0*M_PI);
  if (a < 0) a += 2.0*M_PI;
  return a - M_PI; // now in (-pi, pi]
}

// find rotation matrix index k in {0,1,2,3} that minimizes |(theta_new + k*pi/2) - theta_old|
inline int find_rotation_matrix(double theta_new, double theta_old) {
    int k = -1;
    double max = 1e34;
    for (int i = 0; i < 4; ++i) {
        double angle = i * M_PI_2;
        double diff = std::fabs(wrap_pi((theta_new + angle) - theta_old));
        if (diff < max) {
            max = diff;
            k = i;
        }
    }
    return k;
}

inline Point rotateVector(const Point &u, int k) {
    switch (k) {
        case 0: return u;
        case 1: return Point{ -u[1], u[0] };
        case 2: return Point{ -u[0], -u[1] };
        case 3: return Point{ u[1], -u[0] };
        default: return u;
    }
}


// Undirected edge key, normalized to (min, max) so that the two orientations
// of an edge hash and compare equal. Shared by everything that keys a set or
// map on mesh edges -- see CutMesh::EdgeKey and HarmonicCut::EdgeKey, which
// are aliases of it, so cut edge sets from either are interchangeable.
struct MeshEdgeKey {
    int a = -1;
    int b = -1;
    MeshEdgeKey() = default;
    MeshEdgeKey(int u, int v) {
        if (u < v) { a = u; b = v; }
        else       { a = v; b = u; }
    }
    bool operator==(const MeshEdgeKey &o) const { return a == o.a && b == o.b; }
};

struct MeshEdgeKeyHash {
    std::size_t operator()(const MeshEdgeKey &k) const {
        return static_cast<std::size_t>(k.a) * 73856093u ^ static_cast<std::size_t>(k.b) * 19349663u;
    }
};

class Mesh {
 public:
    std::vector<Point> vertices; // List of 2D points
    std::vector<Triangle> triangles; // List of triangles defined by vertex indices
    std::vector<std::array<int, 3>> triangleAdjacency; // Adjacency info for triangles (0: left, 1: right, 2: below, -1 indicates boundary)
    std::vector<std::array<int, 2>> boundaryTriangles; // List of boundary triangle indices and their corresponding edge (0,1,2)
    std::vector<std::array<int,3>> cornerTriangles; // List of corner triangle indices and their corresponding boundary edges
    std::vector<int> boundaryVertices; // List of vertex indices that lie on the boundary
    std::vector<bool> isBoundaryVertex; // Boolean flag per vertex indicating if it's a boundary vertex
    // Material ID per triangle, one entry per triangle. Single-material meshes
    // carry all 1s. An .obj sets these from its `usemtl mat<id>` lines, which
    // Mesh2Dgmsh writes from the Physical Surface tags of the .geo.
    std::vector<int> triangleMatId;

    // Edge data structures
    std::vector<std::array<int, 2>> edges; // Unique edges: edge index -> [v0, v1] vertex ids
    std::vector<std::array<int, 3>> triangleEdges; // Triangle -> edge map: tri index -> [e0, e1, e2] edge indices
    std::vector<std::array<int, 2>> edgeTriangles; // Edge -> triangle map: edge index -> [t0, t1] triangle indices (-1 for boundary)
    std::vector<int> boundaryEdges; // List of boundary edge indices
    std::vector<bool> isBoundaryEdge; // Boolean flag per edge indicating if it's a boundary edge

    // CSR mapping: vertex -> incident triangles in CCW order
    VertexTriangleCSR vertexTriangles;

    Mesh() = default;

    Mesh(const std::string &filename); // Load mesh from an .obj file

    Mesh(const std::vector<Point> &vertices, const std::vector<Triangle> &triangles); // Construct mesh from given vertices and triangles

    // As above, with an explicit material id per triangle. An empty matIds is
    // taken as a single-material mesh and fills triangleMatId with 1s.
    Mesh(const std::vector<Point> &vertices, const std::vector<Triangle> &triangles,
         const std::vector<int> &matIds);

    int findTriangleContainingPoint(const Point &p) const; // Find the triangle index that contains point p, or -1 if not found
};

#endif // __MESH_HXX__
