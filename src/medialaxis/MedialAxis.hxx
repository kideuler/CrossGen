#ifndef _MEDIAL_AXIS_HXX_
#define _MEDIAL_AXIS_HXX_

#include "mesh/Mesh.hxx"
#include <unordered_set>

enum class TopMakerNodeType {
    Normal,
    Corner,
    Dangle,
};

struct BoundaryVertex {
    int vertexIndex; // Index of the vertex in the original mesh
    int nextBoundaryVertex; // Index of the next boundary vertex in the same boundary loop
    int prevBoundaryVertex; // Index of the previous boundary vertex in the same boundary loop
};

struct MedialVertex {
    Point coord; // 2D position of the medial vertex
    int triangleIndex; // Index of the corresponding triangle in the mesh
    std::unordered_set<int> neighbors; // Set of neighboring medial vertex indices (connected by medial edges)
    std::unordered_set<int> touchPoints; // Set of triangle vertex indices that are equidistant to this medial vertex (the "touch points" on the original mesh)
    int degree; // Number of connected medial edges (degree of the vertex in the medial graph)
    double radius; // Distance from the medial vertex to any of the triangle's vertices.
    bool active; // Whether this medial vertex is active (not pruned) in the current iteration of TopMaker
    TopMakerNodeType nodeType; // Type of the node for TopMaker (Normal, Corner, Dangle)
};

class MedialAxis {
    public:
        std::shared_ptr<Mesh> mesh;
        std::vector<MedialVertex> medialVertices; // List of medial axis vertices
        std::vector<Edge> medialEdges; // List of medial axis edges
        std::vector<BoundaryVertex> boundaryVertices; // List of boundary vertices and their connectivity
        std::vector<bool> sharpVertices; // Whether each boundary vertex is a "sharp" vertex geometrically.

        MedialAxis(std::shared_ptr<Mesh> mesh);
    private:
        Point computeCircumcenter(int triIndex);
        
};

#endif // _MEDIAL_AXIS_HXX_