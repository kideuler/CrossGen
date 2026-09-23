#ifndef __HARMONICCUT_HXX__
#define __HARMONICCUT_HXX__

#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "mesh/Mesh.hxx"

// Cut generation of
//   Wang, Ren, Fang, Lin, Xu, Bao and Huang, "IGA-suitable planar
//   parameterization with patch structure simplification of closed-form
//   polysquare", CMAME 392 (2022) 114678, Section 4.1.
//
// A planar mesh M_beta with beta voids carries beta non-trivial harmonic
// 1-forms, and those are the extra degrees of freedom that promote a common
// (exact-form) polysquare into a closed-form one -- the difference between
// putting four forced corners on an annulus and putting none. Section 4.1
// exposes them as transitions on beta cuts: one cut per void, each connecting
// an inner boundary loop to the outer one, breaking the void open. Applied
// recursively, M_beta -> M_{beta-1} -> ... -> M_0, this ends at a disk.
//
// The cuts must be disjoint -- from each other and from the boundary, except
// at their two endpoints -- because that is what makes the beta harmonic forms
// independently representable, one per cut. That requirement is the reason the
// paper rules out the usual spanning-tree cut graph (Bommes et al. 2009), and
// it is the difference between this class and CutMesh: CutMesh cuts with a
// dual spanning tree plus shortest paths dragging each singularity to the
// boundary, which yields a cut *tree* that touches the boundary all over. Both
// produce a disk; only this one produces the beta independent transitions
// Sec. 4.1 asks for, and it never depends on where the singularities of some
// input field happen to sit.
//
// Sec. 4.1 also prefers short cuts, to keep down the number of transition
// constraints (Eq. 8) later, so the endpoints are the nearest vertex pair
// between the two loops and the path between them is a shortest path through
// interior edges. The result is independent of the specific choice of cut
// anyway -- a standard property of seamless parameterization (see the
// discussion around Fig. 4 of Aigerman et al. 2015) -- so "short" is about
// cost, not correctness.
class HarmonicCut {
public:
    using EdgeKey = MeshEdgeKey;
    using EdgeKeyHash = MeshEdgeKeyHash;

    // One cut C_gamma: the vertex path from p on the already-merged boundary
    // to q on the void being opened, and its edges. Downstream stages need the
    // cuts *individually*, not as one edge set -- Eq. (7) extracts one
    // rotation Pi_gamma per cut -- which is what this grouping is for.
    struct Cut {
        std::vector<int> path;      // original vertex ids, p first, q last
        std::vector<EdgeKey> edges; // path edges, all interior
        int fromLoop = -1;          // boundary loop of p (already merged)
        int toLoop = -1;            // boundary loop of q (the void being opened)
        double length = 0.0;        // geometric length
    };

    struct Report {
        int boundaryLoops = 0;          // 1 + beta on a well-formed planar mesh
        int voids = 0;                  // beta
        int cutsMade = 0;               // beta when every void was opened
        int eulerCharacteristic = 0;    // of the cut mesh; 1 for a disk
        int boundaryComponents = 0;     // of the cut mesh; 1 for a disk
        bool trianglesConnected = false;
        bool isDisk = false;
        std::vector<std::string> messages;
    };

    // `openVoids` false makes every void stay shut: no cuts, and a "cut mesh"
    // that is the input mesh vertex for vertex. Everything downstream still
    // works -- there are simply no seams, no transitions and no harmonic
    // degrees of freedom -- so it is the switch between a *closed-form*
    // polysquare and a *common* one, and it is a real choice rather than a
    // degradation.
    //
    // Sec. 4.1 wants the cuts because its goal is the closed form: an annulus
    // then maps to a ring with no corners on it at all, which is one spline
    // patch instead of four and is what an IGA solver wants to be handed.
    // A *block decomposition* wants the opposite. Its blocks are
    // four-cornered, so a cornerless ring is not a block it can represent, and
    // the ring's image is refused whole (`BlockLayout::DecompositionReport`
    // counts it under notFourSided); the four reflex corners a common
    // polysquare is forced to put on the hole are exactly the four the
    // motorcycle graph needs to send rays from. Measured over
    // data/meshes/singlemat, every model that lost coverage but one was a
    // holed one, and the worst of them -- geom026, geom032 -- were the ones
    // whose cuts carried a 180-degree transition, which folds the image over
    // itself. See [[umber-no-cut-parameterization]].
    //
    // The default stays true because MERIDIAN::ConeCut uses this class for
    // Stage 2, where a disk is exactly what is wanted.
    explicit HarmonicCut(std::shared_ptr<Mesh> mesh, bool openVoids = true);

    // Whether the voids were opened at all -- see the constructor.
    bool getOpenVoids() const { return openVoids; }

    const Mesh& getOriginalMesh() const { return *orig; }
    std::shared_ptr<Mesh> getOriginalMeshPtr() const { return orig; }

    // The disk M_0: the same triangles, with vertices duplicated along the
    // cuts.
    const Mesh& getCutMesh() const { return cut; }

    // All cut edges at once, in original-mesh vertex indices. This is C of
    // Eq. (2) for the frame field optimization.
    const std::unordered_set<EdgeKey, EdgeKeyHash>& getCutEdges() const { return cutEdges; }

    // The cuts one by one, in the order they were made.
    const std::vector<Cut>& getCuts() const { return cuts; }

    // Boundary loops of the *input* mesh, each an ordered vertex ring, and the
    // index of the outer one (the loop of largest enclosed area).
    const std::vector<std::vector<int>>& getBoundaryLoops() const { return boundaryLoops; }
    int getOuterLoop() const { return outerLoop; }

    // Mapping: cut-mesh vertex -> original vertex, and its inverse.
    const std::vector<int>& getCutVertexToOriginal() const { return cutVertToOrig; }
    const std::vector<std::vector<int>>& getOriginalToCutVertices() const { return origToCutVerts; }

    // Whether the cut mesh came out a disk, plus what went wrong if not.
    const Report& getReport() const { return report; }

    bool writeOBJ(const std::string &filename) const;

private:
    void findBoundaryLoops();
    void buildVertexAdjacency();

    // The Sec. 4.1 recursion: open one void per iteration, nearest loop first.
    void generateCuts();

    // Shortest path from any vertex in `sources` to any vertex in `targets`,
    // through interior edges only, visiting no vertex marked in `blocked` --
    // which holds every boundary vertex and every vertex already touched by a
    // cut, so the result satisfies the disjointness requirement by
    // construction. Endpoints are exempt from the block. Empty if no such path
    // exists.
    std::vector<int> shortestInteriorPath(const std::vector<int> &sources,
                                          const std::unordered_set<int> &targets,
                                          const std::vector<char> &blocked) const;

    void applyCut(const std::vector<int> &path, int fromLoop, int toLoop);

    void buildExplicitCutMesh();
    Report checkCutMesh() const;

    std::shared_ptr<Mesh> orig;
    bool openVoids = true;
    Mesh cut;
    Report report;

    std::vector<std::vector<int>> boundaryLoops; // ordered rings, original ids
    std::vector<int> vertexLoop;                 // vertex -> loop id, -1 interior
    int outerLoop = -1;

    std::vector<std::vector<std::pair<int, double>>> adjacency; // primal graph, weighted
    std::vector<char> blockedVertex;  // boundary vertices and vertices on cuts

    std::vector<Cut> cuts;
    std::unordered_set<EdgeKey, EdgeKeyHash> cutEdges;

    std::vector<int> cutVertToOrig;
    std::vector<std::vector<int>> origToCutVerts;
};

#endif // __HARMONICCUT_HXX__
