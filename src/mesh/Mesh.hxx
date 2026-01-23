#ifndef __MESH_HXX__
#define __MESH_HXX__

#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cctype>
#include <stdexcept>
#include "VertexTriangleCSR.hxx"


typedef std::array<double, 2> Point;
typedef std::array<int, 3> Triangle;

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


class Mesh {
 public:
    std::vector<Point> vertices; // List of 2D points
    std::vector<Triangle> triangles; // List of triangles defined by vertex indices
    std::vector<std::array<int, 3>> triangleAdjacency; // Adjacency info for triangles (0: left, 1: right, 2: below, -1 indicates boundary)
    std::vector<std::array<int, 2>> boundaryTriangles; // List of boundary triangle indices and their corresponding edge (0,1,2)
    std::vector<std::array<int,3>> cornerTriangles; // List of corner triangle indices and their corresponding boundary edges

    // CSR mapping: vertex -> incident triangles in CCW order
    VertexTriangleCSR vertexTriangles;

    Mesh() = default;

    Mesh(const std::string &filename); // Load mesh from an .obj file
};

#endif // __MESH_HXX__
