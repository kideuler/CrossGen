#ifndef _MEDIAL_AXIS_HXX_
#define _MEDIAL_AXIS_HXX_

#include "mesh/Mesh.hxx"
#include <deque>
#include <stdexcept>

constexpr int edgeToOppositeVertex[3] = {2, 0, 1}; // For triangle edge i, the opposite vertex is at index edgeToOppositeVertex[i]

// An ordered chain of boundary edges around a boundary vertex.
// Edges are maintained in counterclockwise order so that consecutive
// edges share a vertex: the "back vertex" of edge i equals the
// "front vertex" of edge i+1.  Querying front() and back() always
// returns the first and last edges of the current chain bounds.
//
// Usage:  chain.insert(edgeIdx, v0, v1, vertices)
//   where v0, v1 are the two endpoint vertex-indices of the edge and
//   vertices is the coordinate array (needed for the CCW orientation test).
struct BoundaryEdgeChain {
    std::deque<int> edges;    // edge indices in chain order
    std::deque<int> verts;    // directed vertex sequence (size = edges.size()+1)
                              // verts[i] is the "start" vertex of edges[i]
                              // verts.back() is the free end after the last edge

    // Insert an edge into the chain, automatically placing it at the
    // correct end so the chain stays connected.  For the very first
    // two edges the CCW test (cross product at the shared vertex) is
    // used to decide the ordering.
    void insert(int edgeIdx, int v0, int v1,
                const std::vector<Point>& vertices) {
        if (edges.empty()) {
            // First edge — just store it with an arbitrary direction;
            // orientation will be fixed when the second edge arrives.
            edges.push_back(edgeIdx);
            verts.push_back(v0);
            verts.push_back(v1);
            return;
        }

        int frontV = verts.front();   // free vertex at the start of the chain
        int backV  = verts.back();    // free vertex at the end of the chain

        // Determine which endpoint of the new edge is the shared vertex
        // with an existing chain end.
        if (v0 == backV || v1 == backV) {
            // Connects to the back of the chain
            int shared = backV;
            int other  = (v0 == shared) ? v1 : v0;
            edges.push_back(edgeIdx);
            verts.push_back(other);
        } else if (v0 == frontV || v1 == frontV) {
            // Connects to the front of the chain
            int shared = frontV;
            int other  = (v0 == shared) ? v1 : v0;
            edges.push_front(edgeIdx);
            verts.push_front(other);
        } else {
            // The edge doesn't connect to either end — this means the
            // chain currently has only one edge and its direction may
            // need flipping, OR something is wrong.
            // Try reversing the chain direction and retry.
            std::reverse(verts.begin(), verts.end());
            frontV = verts.front();
            backV  = verts.back();
            if (v0 == backV || v1 == backV) {
                int shared = backV;
                int other  = (v0 == shared) ? v1 : v0;
                edges.push_back(edgeIdx);
                verts.push_back(other);
            } else if (v0 == frontV || v1 == frontV) {
                int shared = frontV;
                int other  = (v0 == shared) ? v1 : v0;
                edges.push_front(edgeIdx);
                verts.push_front(other);
            } else {
                throw std::runtime_error(
                    "BoundaryEdgeChain::insert: edge does not connect to chain");
            }
        }

        // After inserting the second edge, enforce CCW ordering.
        // The "owner" vertex (the boundary vertex this chain belongs to)
        // sits in the middle of verts when size==3: verts = {A, owner, B}.
        // CCW means cross(owner→A, owner→B) > 0  ⟹  A is before owner,
        // B is after owner in the CCW sense.  If cross < 0 we reverse.
        if (edges.size() == 2 && verts.size() == 3) {
            int A     = verts[0];
            int owner = verts[1]; // shared vertex = the boundary node
            int B     = verts[2];
            const Point& pA = vertices[A];
            const Point& pO = vertices[owner];
            const Point& pB = vertices[B];
            double cross = (pA[0] - pO[0]) * (pB[1] - pO[1])
                         - (pA[1] - pO[1]) * (pB[0] - pO[0]);
            if (cross < 0.0) {
                // Reverse the whole chain to restore CCW order
                std::reverse(edges.begin(), edges.end());
                std::reverse(verts.begin(), verts.end());
            }
        }
    }

    // Number of edges in the chain
    std::size_t size() const { return edges.size(); }

    // Indexed access (0-based)
    int  operator[](std::size_t i) const { return edges[i]; }
    int& operator[](std::size_t i)       { return edges[i]; }

    // First / last edge in the CCW chain
    int front() const { return edges.front(); }
    int back()  const { return edges.back();  }
};

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
        std::unordered_map<int, BoundaryEdgeChain> bdyNodeToEdges; // For boundary nodes: CCW-ordered chain of incident boundary edges
        std::vector<std::array<int,2>> EdgeToMedialNode; // Maps each internal edge index to the corresponding medial node index (size = number of edges, -1 for boundary edges)

        MedialAxis(std::shared_ptr<Mesh> mesh);

        void constructMappingPhase1();
        void constructMappingPhase2();
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