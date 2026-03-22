#ifndef _MEDIAL_AXIS_HXX_
#define _MEDIAL_AXIS_HXX_

#include "mesh/Mesh.hxx"

constexpr int edgeToOppositeVertex[3] = {2, 0, 1}; // For triangle edge i, the opposite vertex is at index edgeToOppositeVertex[i]

// Represents a precise location on the domain boundary
struct BoundaryMappedPoint {
    int edgeIndex; // Index into mesh->edges (must be a boundary edge)
    int8_t lid; // local edge index to triangle (0, 1, or 2)
    double t;      // Parametric value [0.0, 1.0]. exactly 0.5 is the midpoint.

    // Evaluates the physical 2D coordinate on demand
    Point evaluate(const std::shared_ptr<Mesh>& mesh) const {
        const auto& edge = mesh->edges[edgeIndex];
        const Point& v0 = mesh->vertices[edge[0]];
        const Point& v1 = mesh->vertices[edge[1]];
        return v0 + (v1 - v0) * t;
    }
};

struct MedialNode {
    int id; // Matches the original mesh->triangles index initially
    Point coord; // The circumcenter (or resampled location)
    
    // The phi^-1 mapping. 
    // - Size 2 for corner triangles (degree 1)
    // - Size 2 for boundary triangles (degree 2)
    // - Size 3 for internal triangles (degree 3)
    std::vector<BoundaryMappedPoint> preImage; 
    
    // Adjacent MedialNode IDs. Replaces std::vector<Edge> medialEdges.
    std::vector<int> neighbors; 

    int degree; // 1 for end, 2 for regular, 3 for branch. Can be derived from neighbors.size() but stored for convenience.
};

class MedialAxis {
    public:
        std::shared_ptr<Mesh> mesh;
        std::vector<MedialNode> medialNodes; // List of medial axis vertices
        std::vector<Edge> medialEdges; // List of medial axis edges
        std::unordered_map<int, std::array<int, 2>> bdyNodeToEdges; // For boundary nodes: the ordered boundary edges connected to this node
        std::vector<std::array<int,2>> EdgeToMedialNode; // Maps each internal edge index to the corresponding medial node index (size = number of edges, -1 for boundary edges)

        MedialAxis(std::shared_ptr<Mesh> mesh);

        void constructPhiInverseMapping();
    private:
        Point computeCircumcenter(int triIndex);
        bool areCollinear(const Point& a, const Point& b, const Point& c);
        // compute the intersection of two lines defined by (p1, p2) and (p3, p4). Assumes lines are not parallel.
        Point computeLineIntersection(const Point& p1, const Point& p2, const Point& p3, const Point& p4);
        std::pair<double, bool> computeLineEdgeIntersection(const Point& p1, const Point& p2, const Point& e0, const Point& e1);
        // compute the projection of point p onto the line defined by edge (e0, e1). Returns the parametric value t along the edge and a boolean indicating if the projection is valid (i.e., if it falls within the edge segment).
        std::pair<double, bool> computePointEdgeProjection(const Point& p, const Point& e0, const Point& e1);
};

#endif // _MEDIAL_AXIS_HXX_