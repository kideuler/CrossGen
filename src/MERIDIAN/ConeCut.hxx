#ifndef __CONECUT_HXX__
#define __CONECUT_HXX__

#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include "MERIDIAN/ConeSingularities.hxx"
#include "MERIDIAN/Interfaces.hxx"
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
// ### The material interface network
//
// On a multi-material domain there is a second set of curves G must keep away
// from, and Definition 2.1 says nothing about it because the paper's features
// are an input to Stage 0 and its cuts are chosen long before. An interface is
// a curve the layout has to keep: E3 of Sec. 3.3 holds each branch of it on a
// coordinate line, and E6 holds the sectors at each of its nodes to whole
// quarter turns. Both of those are statements about Psi, and Psi lives on
// Omega, so both of them are read through the cut -- and a cut that touches the
// network breaks them in three different ways:
//
//   * An arc running *along* a branch makes that branch a seam. Its two sides
//     are then two chains of Omega joined by a quarter-turn transition, and
//     Stage 5 labels both of them from the one direction the node quantisation
//     gave the branch. E3 is then asked to hold u constant on a curve and on
//     its own image under a quarter turn, which is a map that does not exist.
//     On geom002 this is visible before Stage 5 even runs: a 56-edge interface
//     arrives at the labelling as 66 feature edges, the ten the cut ran along
//     counted twice.
//
//   * An arc *crossing* a branch splits it into two chains the same way, with
//     the same contradiction between their shared label and the transition
//     between them.
//
//   * An arc through, or ending at, a *node* splits the node's fan. Two rays
//     either side of the arc are then tangents at two different children of the
//     node, which are not in the same frame, and the sector between them is
//     dropped from E6 (SubdomainLabels::InterfaceCorner::spansCut). What the
//     sector was asserting -- that the layout turns q_s quarters there -- is
//     then asserted by nothing.
//
// Over the fourteen models of data/meshes/multimat that is the whole of the
// difference between a layout and a near-miss: the five whose cut touched the
// network are exactly the five that failed Stage 6, and the nine whose cut
// happened to miss it are exactly the nine that passed.
//
// So the network is routed around, by Options::interfaceAvoidance: every vertex
// of it is charged more than the longest path in the mesh could ever cost, so
// Dijkstra takes any interface-free detour, however long, over any route that
// touches the network at all. It is a penalty and not a wall because the detour
// need not exist -- a cone inside an inclusion is fenced in by a closed
// interface loop and *has* to cross it to reach dS -- and a cut that has to
// cross should cross once, transversally, at the cheapest vertex available,
// rather than fail. What contact is left is measured in Report::
// interfaceEdgesOnCut and friends rather than assumed away.
//
// The obvious next move on a forced crossing was tried and is not an
// improvement. All four of geom009's interior cones sit inside one ellipse, and
// under conesToBoundary each pays for a crossing of its own; letting the last
// three stop on the first one's arc instead -- still inside the ellipse, so
// they cross nothing -- takes the crossings from four to one and the cut from
// 4.09 to 3.09. It also puts three interior junctions in G, and Stage 6 comes
// out worse for it: 17 failed checks against 14, having lost E6 and Q2 to the
// stacked quarter-turns at the junctions and gained nothing on the feature
// chains, since one split chain is as unlabellable as four. The junction is the
// more expensive of the two after all, which is what Options::conesToBoundary
// says in the first place. A forced crossing is a Stage 5 problem and not a
// Stage 2 one: what a chain split by a seam needs is the label on its far side
// turned by the seam's own transition, and until Stage 5 does that no routing
// choice here rescues it.
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

        // Keep G off the material interface network. See the class comment for
        // what touching it costs; this is the multiplier on the penalty, in
        // units of "every edge of the mesh, end to end", which is an upper
        // bound on any path a cone arc could take. At 1 a route that avoids the
        // network always beats one that does not, no matter how much longer it
        // is, and among routes with the same contact the shortest still wins --
        // the penalty is a lexicographic order, not a weighting to tune. A node
        // of the network costs twice a plain vertex of a branch, so a crossing
        // forced by an inclusion prefers the middle of a branch to a junction.
        //
        // Zero restores the interface-blind routing, which is what the
        // single-material pipeline does in any case: with no network there are
        // no penalised vertices and this option has no effect at all. Below 1
        // it stops being an order and becomes a weighting, where a long enough
        // detour is worth a crossing after all; there is no model in the corpus
        // that wants one, so 1 is both the default and the only value measured.
        double interfaceAvoidance = 1.0;
    };

    // One arc of G routed from an interior cone to the rest of the graph.
    struct ConePath {
        int cone = -1;              // the singular vertex, path.front()
        int index = 0;              // I(cone), carried through for convenience
        std::vector<int> path;      // original vertex ids, cone first
        std::vector<EdgeKey> edges;
        double length = 0.0;
        bool reachedBoundary = false; // ended on dS rather than on an earlier arc
        // Vertices of the arc that lie on the material interface network, and
        // how many of those are nodes of it. `startsOnNode` separates the one
        // kind of contact routing cannot remove -- a cone that *is* a node of
        // the network has to leave from one -- from the crossings, which are
        // forced only by a cone with no interface-free route to dS.
        int interfaceVerts = 0;
        int interfaceNodes = 0;
        bool startsOnNode = false;
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

        // What G still has in common with the material interface network, all
        // zero on a single-material mesh and all zero on every model of
        // data/meshes/multimat. Non-zero is legal and is not a failure -- an
        // inclusion leaves a cone no interface-free route to dS -- but each one
        // costs E3 a feature chain or E6 a sector, so they are counted.
        int interfaceEdgesOnCut = 0;    // arcs of G running *along* an interface
        int interfaceVertsOnCut = 0;    // vertices of G on the network at all
        int interfaceNodesOnCut = 0;    // ... of those, nodes of it
        int conesTouchingInterface = 0; // cone arcs with any of the above on them
        // Cone arcs whose *source* is itself a node of the network. Unavoidable
        // -- a junction carrying an index is a cone and Def. 2.1 drags it to dS
        // like any other -- and counted separately because it is the one kind of
        // contact no routing can remove.
        int conesOnInterfaceNodes = 0;

        int eulerCharacteristic = 0;   // of Omega; 1 for a disk
        int boundaryComponents = 0;    // of Omega; 1 for a disk
        int triangleComponents = 0;
        bool trianglesConnected = false;
        bool isDisk = false;
        bool allConesOnBoundary = false; // P subset of G union dS, as Def. 2.1 asks

        std::vector<std::string> messages;
    };

    // `interfaces` is the Stage 0b network, or null on a single-material mesh
    // (or to route as if there were none). It is read and not kept.
    ConeCut(std::shared_ptr<Mesh> mesh, const ConeSingularities &cones,
            const Options &options, const Interfaces *interfaces = nullptr);
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
    // in the same weighted edge graph the paths are routed in -- penalties and
    // all -- so the ordering and the routing agree with each other.
    void routeConePaths(const ConeSingularities &cones);

    // The cost of arriving at each vertex, over and above the edge length that
    // got there: zero everywhere except on the interface network. Empty when
    // there is no network or Options::interfaceAvoidance is zero, in which case
    // every lookup below is a no-op.
    void buildVertexPenalty(const Interfaces *interfaces);
    double penaltyAt(int v) const {
        return vertexPenalty.empty() ? 0.0 : vertexPenalty[v];
    }

    // Dijkstra from `source` to any vertex of `targets`, weight l_ij plus the
    // arrival penalty of the vertex reached, never passing through a vertex
    // marked in `forbidden` (other cones) and never through a target (which
    // would run the cut along the boundary).
    std::vector<int> shortestPathToSet(int source,
                                       const std::vector<char> &targets,
                                       const std::vector<char> &forbidden) const;

    void buildExplicitCutMesh();
    void check(const ConeSingularities &cones, const Interfaces *interfaces);

    std::shared_ptr<Mesh> orig;
    Options opts;
    std::unique_ptr<HarmonicCut> harmonic;

    Mesh cut;
    Report report;

    std::unordered_set<EdgeKey, EdgeKeyHash> cutEdges;
    std::vector<ConePath> conePaths;

    std::vector<int> cutVertToOrig;
    std::vector<std::vector<int>> origToCutVerts;

    std::vector<double> vertexPenalty;
};

#endif // __CONECUT_HXX__
