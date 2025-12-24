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