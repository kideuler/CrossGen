#ifndef _MEDIAL_AXIS_HXX_
#define _MEDIAL_AXIS_HXX_

#include "mesh/Mesh.hxx"
#include <unordered_set>

const double MEDIAL_VERTEX_MERGE_TOLERANCE = 0.01; // Tolerance for merging close medial vertices in deduplication

enum class TopMakerNodeType {
    Normal,
    Corner,
    Dangle,
    NONE = -1
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
    std::unordered_set<int> polyLines; // Set of polyline indices that this medial vertex belongs to (for later use in TopMaker)
    int degree; // Number of connected medial edges (degree of the vertex in the medial graph)
    double radius; // Distance from the medial vertex to any of the triangle's vertices.
    bool active; // Whether this medial vertex is active (not pruned) in the current iteration of TopMaker
    TopMakerNodeType nodeType = TopMakerNodeType::NONE; // Type of the node for TopMaker (Normal, Corner, Dangle)
    int cornerIndex = -1; // If this is a corner node, the index of the corresponding sharp corner vertex in the original mesh; otherwise -1
    bool partOfPolyline = false; // Whether this medial vertex is part of any polyline
    int timesFused = 0; // Number of times this vertex has been fused into another during deduplication (if more than 1 it is a sign this vertex is a dangle)
};

class MedialAxis {
    public:
        std::shared_ptr<Mesh> mesh;
        std::vector<MedialVertex> medialVertices; // List of medial axis vertices
        std::vector<Edge> medialEdges; // List of medial axis edges
        std::vector<BoundaryVertex> boundaryVertices; // List of boundary vertices and their connectivity
        std::vector<bool> sharpVertices; // Whether each boundary vertex is a "sharp" vertex geometrically.
        std::vector<std::vector<int>> polyLines; // List of polylines (each represented as a vector of vertex indices)

        MedialAxis(std::shared_ptr<Mesh> mesh);

        void deduplicateMedialVertices();
        void createPolylines();
        void classifyMedialVertices();
    private:
        Point computeCircumcenter(int triIndex);
        
};

#endif // _MEDIAL_AXIS_HXX_