#ifndef __CUTMESH_HXX__
#define __CUTMESH_HXX__

#include <array>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "Mesh.hxx"

// Forward declaration (definition in PolyVectors.hxx)
class PolyField;

// Cut a triangle mesh (represented by PolyField::getMesh()) into a topological disk.
//
// Implements the MIQ-style cutting strategy:
//  1) Build a spanning tree of the dual graph and cut the complement primal edges.
//  2) Ensure all singularities lie on the boundary by cutting shortest paths from
//     each singular vertex to the current cut graph / boundary.
//  3) Materialize an explicit cut mesh by duplicating vertices along cut edges.
//
// The resulting mesh is stored in `getCutMesh()` as a regular Mesh with updated
// adjacency and boundary lists.
class CutMesh {
public:
    // Undirected edge key (min,max) suitable for hash containers.
    struct EdgeKey {
        int a = -1;
        int b = -1;
        EdgeKey() = default;
        EdgeKey(int u, int v) {
            if (u < v) { a = u; b = v; }
            else       { a = v; b = u; }
        }
        bool operator==(const EdgeKey &o) const { return a == o.a && b == o.b; }
    };

    struct EdgeKeyHash {
        std::size_t operator()(const EdgeKey &k) const {
            return static_cast<std::size_t>(k.a) * 73856093u ^ static_cast<std::size_t>(k.b) * 19349663u;
        }
    };

    struct SanityReport {
        bool trianglesConnected = false;
        int triangleComponents = 0;
        int boundaryComponents = 0;
        int eulerCharacteristic = 0;
        bool looksLikeDisk = false;
        bool allSingularitiesOnBoundary = false;
        std::vector<std::string> messages;
    };

    explicit CutMesh(const PolyField &field);

    const Mesh& getOriginalMesh() const { return orig; }
    const Mesh& getCutMesh() const { return cut; }

    const std::unordered_set<EdgeKey, EdgeKeyHash>& getCutEdges() const { return cutEdges; }

    // Subset of cut edges that were added specifically to connect singularities
    // to the boundary / existing cut graph (step 2 in the MIQ-style strategy).
    // These are expressed in original-mesh vertex indices, same as getCutEdges().
    const std::unordered_set<EdgeKey, EdgeKeyHash>& getSingularityPathCutEdges() const { return singularityPathCutEdges; }

    // Mapping: cut-mesh vertex -> original-mesh vertex.
    const std::vector<int>& getCutVertexToOriginal() const { return cutVertToOrig; }
    // Mapping: original vertex -> list of cut-mesh vertices.
    const std::vector<std::vector<int>>& getOriginalToCutVertices() const { return origToCutVerts; }

    // Singularities (original vertex id, index4)
    const std::vector<std::pair<int,int>>& getSingularities() const { return singularities; }

    // Write the cut mesh as an OBJ file (z=0).
    bool writeOBJ(const std::string &filename) const;

    // Basic topological sanity checks for the cut mesh.
    SanityReport sanityCheck() const;

private:
    Mesh orig;
    Mesh cut;

    std::vector<std::pair<int,int>> singularities;
    std::unordered_set<EdgeKey, EdgeKeyHash> cutEdges;
    std::unordered_set<EdgeKey, EdgeKeyHash> singularityPathCutEdges;

    std::vector<int> cutVertToOrig;
    std::vector<std::vector<int>> origToCutVerts;

    void buildDualSpanningTreeCuts();
    void connectSingularitiesWithShortestPaths();
    void buildExplicitCutMesh();
};

#endif // __CUTMESH_HXX__
