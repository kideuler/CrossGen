#ifndef __BLOCK_LAYOUT_HXX__
#define __BLOCK_LAYOUT_HXX__

#include <vector>

#include "UMBER/MotorcycleGraph.hxx"
#include "tracing/QuadLayout.hxx"

// ---------------------------------------------------------------------------
// The meta-blocks of MotorcycleGraph as a graph rather than as a colouring.
//
// MotorcycleGraph says which block each triangle belongs to, which is all that
// is needed to look at the partition and nothing like enough to operate on it.
// A chord collapse asks questions of the block structure itself -- which block
// is on the other side of this side, which four nodes are this block's corners,
// which arc joins those two nodes -- and none of them can be answered from a
// per-triangle label. This assembles the structure the questions are about: the
// nodes, the arcs between them, and the blocks as cycles of arcs.
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
// assumed: QuadLayout::Report::badFaces counts the blocks that did not come out
// four-sided and ChordCollapse refuses to run a chord through one.
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

    // The graph must have been built.
    explicit BlockLayout(const MotorcycleGraph &graph);

    void build();

    const QuadLayout &getLayout() const { return layout_; }
    QuadLayout &getLayout() { return layout_; }
    const Report &getReport() const { return report_; }

    // The boundary loops of the model, as vertex indices, each with the
    // material on its left -- the orientation QuadLayout's face walk reads the
    // unbounded side off.
    const std::vector<std::vector<int>> &getBoundaryLoops() const { return loops; }

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
    void collectNodes();
    void buildRayArcs();
    void buildBoundaryArcs();

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

    std::vector<std::vector<Split>> splitsOfRay;
    std::vector<std::vector<Split>> splitsOfLoop;
};

#endif // __BLOCK_LAYOUT_HXX__
