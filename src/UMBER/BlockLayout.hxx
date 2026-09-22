#ifndef __BLOCK_LAYOUT_HXX__
#define __BLOCK_LAYOUT_HXX__

#include <string>
#include <vector>

#include "UMBER/MotorcycleGraph.hxx"
#include "mesh/BlockDecomposition.hxx"
#include "tracing/QuadLayout.hxx"

// ---------------------------------------------------------------------------
// The meta-blocks of MotorcycleGraph as a graph rather than as a colouring.
//
// MotorcycleGraph says which block each triangle belongs to, which is all that
// is needed to look at the partition and nothing like enough to operate on it.
// Everything after this point asks questions of the block structure itself --
// which block is on the other side of this side, which four nodes are this
// block's corners, which arc joins those two nodes, whether this side is an
// interface -- and none of them can be answered from a per-triangle label.
// This assembles the structure the questions are about: the nodes, the arcs
// between them, and the blocks as cycles of arcs.
//
// The pieces are already there. The rays were traced as polylines, the corners
// and crossings and exits were found while tracing, and each node knows where
// it sits along the lines through it (MotorcycleGraph::Node::param). So the
// arcs are had by cutting each ray at the nodes on it, in order along it, and
// each boundary loop at the nodes on it -- and the blocks by walking the plane
// graph that leaves. That last part is QuadLayout's, which does exactly this
// for the separatrix partition of Viertel, Osting and Staten and takes an
// arrangement handed to it directly, so it is used here rather than repeated.
//
// What comes out should be a quad layout in the strict sense -- every block
// four-sided, every side a single arc, no T-junctions -- because that is what
// the motorcycle graph is built to give: Sec. 5 of Wang et al. opens by
// claiming the meta-block structure "has neither T-junction nor internal
// singularity", and MotorcycleGraph carries every iso-line through to the
// boundary rather than stopping it on another one, which is what makes the
// claim true. Every interior node is therefore a crossing of two lines, with
// four arcs and four right angles, and every side of every block runs from one
// corner to the next without interruption.
//
// Where that fails it is worth knowing about, so it is measured rather than
// assumed: QuadLayout::Report::badFaces counts the blocks that did not come
// out four-sided, and DecompositionReport turns that into the number that
// matters -- how much of the model is left with no block on it at all.
// ---------------------------------------------------------------------------
class BlockLayout {
public:
    struct Report {
        int boundaryLoops = 0;
        // Rays that stopped in the interior -- MotorcycleGraph::Report::ranOut
        // seen from here. Each leaves a free end in the arrangement, and the
        // blocks around it are not four-sided, so this being nonzero is why a
        // layout can come out with bad faces in it.
        int danglingEnds = 0;
        // Nodes that turned out to be the same place: a crossing found twice
        // because it sits on the edge between two triangles, or one that landed
        // on a corner.
        int mergedNodes = 0;
        // Pieces of a line shorter than the node tolerance, which cannot carry
        // a direction and so are not arcs. Dropping one leaves a hole in the
        // arrangement, so this should be zero.
        int degenerateArcs = 0;
        // Nodes of the motorcycle graph that could not be placed on a boundary
        // loop. Also a hole in the arrangement.
        int unplacedNodes = 0;
    };

    // What became of the layout's faces when it was handed over as the shared
    // block decomposition, which is a different question from whether the
    // layout is sound: a face has to be four-sided *and* have one arc per side
    // to be a block there, and a T-junction costs the second without costing
    // the first.
    struct DecompositionReport {
        int faces = 0;            // faces of the layout offered
        int blocks = 0;           // ... that became blocks
        int notFourSided = 0;     // corners != 4
        int multiArcSides = 0;    // four corners, a side made of several arcs
        int unmatchedSides = 0;   // macro edges with a block on one side only
        int materials = 0;        // distinct materials the blocks landed in
        int straddlingBlocks = 0; // blocks whose samples disagreed on it

        // How much of the model the blocks actually cover, by area of the
        // layout's faces: `covered` over `total`.
        //
        // This is the number to watch, and the counts above are only ever the
        // explanation for it. A face that does not qualify is not a face drawn
        // in a different colour, it is a piece of the model with no block on it
        // -- so it gets no elements, and the mesh has a hole in it exactly
        // there. One refused face can be half the model: on
        // data/meshes/multimat/rocket a single face carrying one T-junction was
        // 51% of it.
        double coveredArea = 0.0;
        double totalArea = 0.0;
    };

    // The fraction of the model the blocks cover, 0 to 1, or 1 when there is
    // nothing to cover.
    static double coverageOf(const DecompositionReport &r) {
        return (r.totalArea > 0.0) ? r.coveredArea / r.totalArea : 1.0;
    }

    // The graph must have been built.
    explicit BlockLayout(const MotorcycleGraph &graph);

    void build();

    // The layout as the one representation ATLAS, MERIDIAN and TORSION also
    // hand their blocks back as (mesh/BlockDecomposition.hxx), so that drawing
    // the blocks, meshing them or asking anything of them is the same code
    // whichever method produced them.
    //
    // Static, and taking the layout rather than reading this object's own, so
    // that any layout can be read as a decomposition and not only the one this
    // object holds. `blockDecomposition()` below is the convenience for the
    // ordinary case of wanting this object's.
    //
    // Only a face with four corners and one arc per side becomes a block. The
    // rest are counted in `out` and simply left out -- the same discipline
    // Arrangement::blockDecomposition follows, and for the same reason: this
    // class holds a decomposition once it already is one, and a face the
    // tracing left open is not one.
    //
    // `mesh` is what the material of each block is read from, by locating a
    // point inside the block's outline. On a single-material model every block
    // comes out material 1 and nothing is searched for twice.
    static BlockDecomposition blockDecompositionOf(const QuadLayout &layout, const Mesh &mesh,
                                                   DecompositionReport *out = nullptr,
                                                   const std::string &source = "UMBER");

    BlockDecomposition blockDecomposition(DecompositionReport *out = nullptr,
                                          const std::string &source = "UMBER") const {
        return blockDecompositionOf(layout_, *mesh, out, source);
    }

    const QuadLayout &getLayout() const { return layout_; }
    QuadLayout &getLayout() { return layout_; }
    const Report &getReport() const { return report_; }

    // The boundary loops of the model, as vertex indices, each with the
    // material on its left -- the orientation QuadLayout's face walk reads the
    // unbounded side off.
    const std::vector<std::vector<int>> &getBoundaryLoops() const { return loops; }

    // The material interface chains, as vertex indices, cut at every vertex
    // the network does not run straight through. Empty on a single-material
    // model. These are the curves MotorcycleGraph::launchFeatureSectors has to
    // agree with about where the network has a node, so they are exposed for
    // checking that it does.
    const std::vector<std::vector<int>> &getInterfaceChains() const { return ichains; }

    // One flag per arc of getLayout(): whether that arc is a piece of a
    // material interface.
    //
    // The decomposition answers the same question with MacroEdge::interface,
    // and better, since it derives it from the materials on the two sides
    // rather than carrying a flag. This is the answer for the *layout*, which
    // is a wider question: it covers the arcs of faces the decomposition
    // refused, and those are exactly the places worth looking at when the
    // coverage is not all of the model.
    const std::vector<char> &getInterfaceArcs() const { return arcIsInterface; }

private:
    // A place along a polyline: p = (index of the point before it) + (fraction
    // of the way along that step), so everything on one line sorts on one
    // number. The same convention MotorcycleGraph::Node::param uses.
    struct Split {
        double p = 0.0;
        int node = -1;
    };

    static double toleranceOf(const MotorcycleGraph &graph);

    void buildBoundaryLoops();
    // The interface network cut into chains at every vertex the network does
    // not simply carry on through -- a junction, an end, a landing on dS.
    // These are sides of the block structure in exactly the sense the boundary
    // loops are: a ray ends on one, a block never crosses one, and a block
    // beside one has it for a side. Without them a multi-material model has
    // rays running to nodes the arrangement cannot place and blocks whose
    // fourth side is not there.
    void buildInterfaceChains();
    void collectNodes();
    void buildRayArcs();
    void buildBoundaryArcs();
    void buildInterfaceArcs();

    // With merging of nodes that are the same place, as QuadLayout does when it
    // builds an arrangement of its own.
    int addNode(const Point &pos, QuadLayout::NodeKind kind);
    int findNode(const Point &pos) const;
    int addArc(std::vector<Point> pts, int a, int b, int ray, bool onBoundary);

    const MotorcycleGraph *graph = nullptr;
    const Mesh *mesh = nullptr;
    double tol = 1e-9;

    QuadLayout layout_;
    Report report_;

    std::vector<QuadLayout::Node> nodes;
    std::vector<QuadLayout::Arc> arcs;

    std::vector<std::vector<int>> loops;
    // vertex -> (loop, index in it), for the boundary vertices.
    std::vector<std::pair<int, int>> loopOfVertex;

    // The interface chains, and mesh edge -> (chain, index of its first vertex
    // in that chain). Keyed by edge rather than by vertex because a vertex can
    // be on two chains -- that is what a junction is -- while an edge is on
    // exactly one.
    std::vector<std::vector<int>> ichains;
    std::vector<char> ichainClosed;
    std::vector<std::pair<int, int>> ichainOfEdge;
    // vertex -> every (chain, index in it) it appears at. A junction is on
    // several, which is why this is a list and the boundary's is a pair.
    std::vector<std::vector<std::pair<int, int>>> ichainPosOfVertex;

    std::vector<std::vector<Split>> splitsOfRay;
    std::vector<std::vector<Split>> splitsOfLoop;
    std::vector<std::vector<Split>> splitsOfChain;
    std::vector<char> arcIsInterface;   // per arc of layout_
};

#endif // __BLOCK_LAYOUT_HXX__
