#ifndef __VERTEX_TRIANGLE_CSR_HXX__
#define __VERTEX_TRIANGLE_CSR_HXX__

#include <vector>
#include <utility>
#include <cmath>

// Forward declaration to avoid circular include with Mesh.hxx
class Mesh;

// Minimal CSR-like structure mapping vertex -> incident triangles (CCW order)
// rowPtr has size |V|+1; entries for vertex i are in colIdx[rowPtr[i]..rowPtr[i+1]-1]
struct VertexTriangleCSR {
    std::vector<int> rowPtr;   // offsets per vertex (size: nVertices+1)
    std::vector<int> colIdx;   // triangle indices

    // Build mapping from a Mesh; triangles are ordered CCW around each vertex
    static VertexTriangleCSR buildFromMesh(const Mesh &mesh);

    // Access helpers
    inline int vertexDegree(int v) const { return rowPtr[v+1] - rowPtr[v]; }
    inline std::pair<const int*, const int*> trianglesForVertex(int v) const {
        return { &colIdx[rowPtr[v]], &colIdx[rowPtr[v+1]] };
    }
};

#endif // __VERTEX_TRIANGLE_CSR_HXX__
