#ifndef __CONECUT_HXX__
#define __CONECUT_HXX__

#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include "MERIDIAN/ConeSingularities.hxx"
#include "Parameterization/HarmonicCut.hxx"
#include "mesh/Mesh.hxx"

// Stage 2 of Shepherd, Gu and Hughes (2022), Section 3.2.2: the cutting graph
// G, and the cut surface Omega = S - G it defines.
//
// Definition 2.1 asks for a G along edges of S such that S - G is a topological
// disk and P is contained in G union dS -- every cone either already sits on
// the boundary or is dragged onto it by a cut. Those are two separate jobs on a
// planar model and they are done by two separate mechanisms here:
//
//   * The voids. A planar mesh with beta holes is not a disk, and each hole
//     needs one arc joining it to the outer boundary. That is exactly the
//     recursion of Wang et al. Sec. 4.1, so HarmonicCut does it: one short cut
//     per void, disjoint from the others and meeting the boundary only at its
//     endpoints. Using it rather than a fresh spanning-tree cut also keeps this
//     stage's arcs interchangeable with the ones UMBER's frame field is built
//     on, since both are keyed on original-mesh vertex indices.
//
//   * The interior cones. HarmonicCut has no notion of a singularity and will
//     not route anything to one, so the cones are dragged out afterwards by the
//     Sec. 3.2.2 procedure -- shortest path in the edge graph from each cone to
//     whatever is already boundary or cut, nearest cone first. Options::
//     conesToBoundary narrows "whatever" to dS alone; see it for why.
//
// Sec. 3.2.2 adds that "it is often preferable that these cuts go to, but not
// through singular points (meaning that a small neighborhood of every singular
// point is a single connected component under the cutting operation)". Note
// what that is and is not: a preference, for computational simplicity, and not
// a condition of Definition 2.1. A cone the graph runs across is split into
// several children -- Fig. 7 is exactly that picture, a valence-five cone in
// three pieces -- and Q2's angle sum in Eq. (1) is written to be read over all
// of them, so such a cut is legal, just more bookkeeping for the seam
// transitions of Stage 4.
//
// The preference is honoured where this class does the routing: a cone arc
// neither passes through nor terminates on another cone, which leaves every
// cone a leaf of G and singly represented in Omega. It cannot be honoured in
// HarmonicCut's void arcs, which pick their endpoints by proximity between
// boundary loops and know nothing about P, so an arc landing on a boundary cone
// is possible. That is counted in Report::conesSplitByCut rather than assumed
// away.
//
// Note that Stage 3 does not consume any of this. Ricci flow (Sec. 3.2.1) runs
// on the uncut triangulation, because a conformal factor per vertex has nothing
// to say about seams; Omega is what Sec. 3.2.2's immersion is laid out on, one
// stage later. The cut is built here, and validated here, so that the metric
// coming out of Stage 3 has somewhere to go.
class ConeCut {
public:
    using EdgeKey = MeshEdgeKey;
    using EdgeKeyHash = MeshEdgeKeyHash;

    struct Options {
        // Where a cone arc is allowed to stop.
        //
        // Sec. 3.2.2 only asks that P end up inside G union dS, so the cheapest
        // arc -- the one to whatever is nearest, which is usually an arc some
        // earlier cone already laid -- satisfies it. That is what this class did
        // unconditionally, and it has a cost the definition does not see. An arc
        // that stops on an earlier arc puts a degree-three junction in the
        // middle of G, and a junction is a vertex of Omega whose one-ring is
        // split into three sectors by two different seams. The immersion of
        // Stage 4 then has to carry two transitions past a single point, and the
        // cone at the far end of the tree is separated from dS by every arc
        // between it and the boundary rather than by its own: the seam it lives
        // on is a chain of arcs, each with its own quarter-turn, and the layout
        // of Stage 6 sees the whole chain as one rigid stack of constraints.
        //
        // With this on, every cone arc runs to dS itself and is vertex-disjoint
        // from every other arc, so G becomes a set of independent slits, each
        // one cone deep. The cut is longer -- that is the trade -- but each cone
        // reaches the boundary through a seam of its own.
        bool conesToBoundary = true;
    };

    // One arc of G routed from an interior cone to the rest of the graph.
    struct ConePath {
        int cone = -1;              // the singular vertex, path.front()
        int index = 0;              // I(cone), carried through for convenience
        std::vector<int> path;      // original vertex ids, cone first
        std::vector<EdgeKey> edges;
        double length = 0.0;
        bool reachedBoundary = false; // ended on dS rather than on an earlier arc
    };

    struct Report {
        int voids = 0;                 // beta, holes in the input
        int harmonicCuts = 0;          // arcs HarmonicCut made, beta when it worked
        int interiorCones = 0;
        int conesRouted = 0;           // interior cones that reached G union dS
        // With Options::conesToBoundary: arcs that made it to dS, and arcs that
        // could not and had to stop on an earlier arc after all. A fallback is
        // not a failure -- the cut is still valid -- but it is the case the
        // option was meant to avoid, so it is counted rather than hidden.
        int conesToBoundary = 0;
        int conesFellBack = 0;
        // Degree-three (or higher) junctions of G away from dS: one per place
        // where an arc stopped on another. Zero is what conesToBoundary buys.
        int interiorJunctions = 0;
        // Cones the graph runs across rather than stopping at, so that they
        // have more than one child in Omega. Legal but not preferred; see the
        // class comment.
        int conesSplitByCut = 0;

        int eulerCharacteristic = 0;   // of Omega; 1 for a disk
        int boundaryComponents = 0;    // of Omega; 1 for a disk
        int triangleComponents = 0;
        bool trianglesConnected = false;
        bool isDisk = false;
        bool allConesOnBoundary = false; // P subset of G union dS, as Def. 2.1 asks

        std::vector<std::string> messages;
    };

    ConeCut(std::shared_ptr<Mesh> mesh, const ConeSingularities &cones,
            const Options &options);
    // Options{} cannot be spelled as a default argument here: the member
    // initialiser above is not yet usable inside the class body.
    ConeCut(std::shared_ptr<Mesh> mesh, const ConeSingularities &cones)
        : ConeCut(std::move(mesh), cones, Options()) {}

    const Mesh& getOriginalMesh() const { return *orig; }
    std::shared_ptr<Mesh> getOriginalMeshPtr() const { return orig; }

    // Omega: the same triangles, with vertices duplicated along G.
    const Mesh& getCutMesh() const { return cut; }

    // G as one edge set, in original-mesh vertex indices: the HarmonicCut arcs
    // and the cone paths together.
    const std::unordered_set<EdgeKey, EdgeKeyHash>& getCutEdges() const { return cutEdges; }

    // The two halves of G separately. Stage 4 fits one transition per arc
    // (Sec. 3.2.2), and the two kinds of arc behave differently -- a void arc
    // has both endpoints on dS, a cone arc has one endpoint at a cone -- so
    // they are kept apart rather than merged into one edge soup.
    const std::vector<HarmonicCut::Cut>& getVoidCuts() const { return harmonic->getCuts(); }
    const std::vector<ConePath>& getConePaths() const { return conePaths; }

    const HarmonicCut& getHarmonicCut() const { return *harmonic; }

    // Omega vertex -> S vertex, and its inverse. A cone at a leaf of G has one
    // child; a vertex the graph passes through has one per angular sector.
    const std::vector<int>& getCutVertexToOriginal() const { return cutVertToOrig; }
    const std::vector<std::vector<int>>& getOriginalToCutVertices() const { return origToCutVerts; }

    const Report& getReport() const { return report; }

    bool writeOBJ(const std::string &filename) const;

private:
    // Sec. 3.2.2, cones nearest the existing graph first. "Nearest" is measured
    // in the same weighted edge graph the paths are routed in, so the ordering
    // and the routing agree with each other.
    void routeConePaths(const ConeSingularities &cones);

    // Dijkstra from `source` to any vertex of `targets`, weight l_ij, never
    // passing through a vertex marked in `forbidden` (other cones) and never
    // through a target (which would run the cut along the boundary).
    std::vector<int> shortestPathToSet(int source,
                                       const std::vector<char> &targets,
                                       const std::vector<char> &forbidden) const;

    void buildExplicitCutMesh();
    void check(const ConeSingularities &cones);

    std::shared_ptr<Mesh> orig;
    Options opts;
    std::unique_ptr<HarmonicCut> harmonic;

    Mesh cut;
    Report report;

    std::unordered_set<EdgeKey, EdgeKeyHash> cutEdges;
    std::vector<ConePath> conePaths;

    std::vector<int> cutVertToOrig;
    std::vector<std::vector<int>> origToCutVerts;
};

#endif // __CONECUT_HXX__
