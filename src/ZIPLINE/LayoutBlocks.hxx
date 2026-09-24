#ifndef __LAYOUT_BLOCKS_HXX__
#define __LAYOUT_BLOCKS_HXX__

#include <array>
#include <string>
#include <vector>

#include "geom/BSpline.hxx"
#include "geom/Polyline.hxx"
#include "geom/Topology.hxx"
#include "mesh/BlockDecomposition.hxx"
#include "mesh/Mesh.hxx"
#include "ZIPLINE/QuadLayout.hxx"

// ---------------------------------------------------------------------------
// The last step of Viertel, Osting and Staten (IMR 2019) as this codebase
// finishes every method: the simplified partition of Sec. 4 read as the shared
// block decomposition (mesh/BlockDecomposition.hxx), on spline geometry built
// in src/geom -- the construction MERIDIAN's Stage 9 (MERIDIAN/SplineFit.hxx)
// makes of an arrangement, made here of a QuadLayout.
//
// ### Which faces are blocks
//
// The paper ends on a T-layout: a partition into four-sided components in
// which a side may be several arcs, because a T-junction may sit on it. This
// codebase does not mesh T-junctions yet, so a component becomes a block only
// when it is a *conforming* quadrilateral:
//
//   * its boundary is one simple cycle (a disk with no arc visited twice and no
//     node visited twice);
//   * it turns exactly four corners (QuadLayout's corner rule, docs/
//     viertel_2019.md Sec. 5);
//   * no node along any of its four sides is anything but a pass-through --
//     a node where two arcs meet and nothing else does, which is what a
//     heteroclinic join and the leftovers of a chord collapse are. A node of
//     valence three or more that the component runs straight past is a
//     T-junction, and the component carrying it is refused.
//
// A refused component is not drawn in another colour and meshed anyway: it
// has no block, so it gets no elements, and the fraction of the model the
// blocks cover (Report::coveredArea over Report::totalArea) is the number that
// says how close the method came. The T-junctions themselves are listed
// (tJunctionNodes()) so the viewer can mark each one.
//
// Sides through pass-through nodes are one side. A chain of arcs joined at
// valence-two nodes is merged into one *macro arc* before anything else is
// done with it, so that a heteroclinic joint in the middle of a side does not
// cost the block it is on.
//
// ### The geometry
//
// One rule, Stage 9's: **fit each macro arc exactly once, and give the curve to
// both blocks that share it.** Two blocks sharing a side then hold the same
// geom::Edge, the B-rep counts it once, and the mesh built on the
// decomposition is conforming because both blocks read the one polyline.
//
//   * A macro arc on the boundary of the model, or on a material interface,
//     is carried exactly, as the polyline of mesh edges it is
//     (geom::Polyline). It is geometry the input gave, and a least-squares
//     cubic through it would round off the model's own features (SplineFit's
//     header has the numbers for why).
//   * Every other macro arc is a traced separatrix, or a blend of two made by a
//     zip collapse, and is fitted by least squares as a cubic B-spline with its
//     two ends pinned to the nodes (geom::fitCurve). Unlike Stage 9 the number
//     of segments is not global: it starts at Options::segments and doubles
//     until the fit is within Options::fitTolerance mean mesh edges of the
//     traced polyline. Stage 9 needs one knot vector everywhere because it
//     blends control nets; here the patch is geom::coonsSurface() of the four
//     side curves, which makes opposite sides compatible itself
//     (geom::makeCompatible), so a long side may carry more segments than a
//     short one without the patch caring.
//   * Each block is the Coons surface of its four side curves, bounded by
//     their four edges: a geom::Face, and together a geom::Shape that can be
//     written as STEP or BREP.
//
// The decomposition handed on (decomposition()) carries each fitted macro arc
// sampled densely enough that its polyline is the spline to plotting
// accuracy, and each exact one as its own polyline, so anything reading the
// decomposition -- the viewer, mesh/BlockQuadMesh -- reads the spline geometry.
//
// ### Materials
//
// With the interfaces honoured (SeparatrixTrace given the network, the field
// aligned to it) every interface is a union of layout arcs, so no face
// straddles one. Each block's material is still read off the model -- a vote
// over sample points located in the triangulation -- rather than trusted, and
// a block whose samples disagree is counted in Report::straddlingBlocks: zero
// is the check that the interfaces really were followed.
// ---------------------------------------------------------------------------
class LayoutBlocks {
public:
    struct Options {
        // Cubic segments a fitted macro arc starts with, and the most it may
        // be given; the count doubles from the first until the fit is within
        // fitTolerance (in mean mesh edges) of the traced polyline.
        int segments = 3;
        int maxSegments = 48;
        double fitTolerance = 0.1;
        // Tikhonov pull towards the chord, as Stage 9 uses it
        // (geom::FitOptions::regularisation).
        double regularisation = 1e-6;
        // Fit boundary macro arcs too, instead of carrying them exactly.
        bool fitBoundaryArcs = false;
        // And interface ones: off for the same reason (SplineFit's header has
        // the numbers -- a cubic through an interface rounds off what the
        // input drew and puts elements astride it).
        bool fitInterfaceArcs = false;
        // Points per cubic segment when a fitted macro arc is written into the
        // decomposition as a polyline.
        int samplesPerSegment = 24;
        // Grid for the fold check of each Coons patch.
        int patchSamples = 8;
    };

    // A maximal chain of layout arcs joined at pass-through nodes.
    struct MacroArc {
        int from = -1, to = -1;        // layout nodes at its two ends
        std::vector<int> darts;        // layout darts, from -> to
        std::vector<Point> points;     // the traced polyline, from -> to
        bool onBoundary = false;
        bool onInterface = false;      // a piece of a material interface
        double length = 0.0;

        // Geometry, for the macro arcs some block uses.
        bool used = false;
        bool exact = false;            // carried as `poly`, not fitted
        geom::Polyline<2> poly;        // when exact
        geom::BSplineCurve<2> spline;  // when fitted
        int segments = 0;
        double maxDeviation = 0.0;     // of the curve from `points`
        geom::Edge edge;               // null if the kernel refused it
        int decompositionEdge = -1;    // index into decomposition().edges
    };

    struct Block {
        int face = -1;                          // QuadLayout face
        std::array<int, 4> corners{{-1, -1, -1, -1}};   // layout nodes, cyclic
        std::array<int, 4> sides{{-1, -1, -1, -1}};     // macro arcs
        // Side k runs corners[k] -> corners[k+1]; forward when that is the
        // macro arc's own from -> to.
        std::array<bool, 4> forward{{true, true, true, true}};
        geom::BSplineSurface<2> surface;        // the Coons patch
        geom::Face brepFace;                    // null if the kernel refused it
        double area = 0.0;                      // of the layout face
        double minCellRatio = 0.0;              // of the sampled patch; <= 0 folded
        int material = 0;
    };

    struct Report {
        int faces = 0;              // components of the layout
        int blocks = 0;             // ... that became blocks
        int notFourSided = 0;       // corners != 4
        int notDisks = 0;           // an arc or a node visited twice
        int tJunctionFaces = 0;     // four corners, a T-junction along a side
        int tJunctions = 0;         // distinct T-junction nodes
        double coveredArea = 0.0;   // area of the blocks' faces
        double totalArea = 0.0;     // area of every face
        // The uncovered area by why: which of the three refusals above cost
        // how much of the model.
        double tJunctionArea = 0.0, notFourSidedArea = 0.0, notDiskArea = 0.0;

        int macroArcs = 0;          // used by some block
        int fittedArcs = 0;
        int exactArcs = 0;
        int maxSegmentsUsed = 0;
        int unconverged = 0;        // fits that hit maxSegments above tolerance
        double maxDeviation = 0.0;  // worst fitted arc, in mean mesh edges
        double meanDeviation = 0.0; // over the fitted arcs, in mean mesh edges

        int foldedPatches = 0;
        // Corner interpolation, read back off the patches: how far a Coons
        // patch's corner is from the node it was built on. Zero, since every
        // fitted curve is pinned at its ends; a check, not a tolerance.
        double maxCornerGap = 0.0;
        // How far a patch's edge is from the macro arc's curve, sampled: the
        // statement that two blocks sharing a side meet along it.
        double maxSideGap = 0.0;

        int materials = 0;
        int straddlingBlocks = 0;

        int brepFaces = 0, brepEdges = 0, brepSharedEdges = 0, brepFreeEdges = 0;
        int brepFailures = 0;
        bool brepValid = false;

        std::vector<std::string> messages;
    };

    // `mesh` is the model the layout was traced on, for the materials and the
    // mean edge the tolerances are in; the layout is read and not kept.
    LayoutBlocks(const QuadLayout &layout, const Mesh &mesh, const Options &opts);
    LayoutBlocks(const QuadLayout &layout, const Mesh &mesh)
        : LayoutBlocks(layout, mesh, Options()) {}

    const Report &getReport() const { return report_; }
    const Options &getOptions() const { return opts_; }

    // The fraction of the model the blocks cover, 0 to 1.
    double coverage() const {
        return report_.totalArea > 0.0 ? report_.coveredArea / report_.totalArea : 1.0;
    }

    const std::vector<MacroArc> &macroArcs() const { return arcs_; }
    const std::vector<Block> &blocks() const { return blocks_; }

    // The shared representation, on the spline geometry. Block i of it is
    // blocks()[i]; its edges are the used macro arcs.
    const BlockDecomposition &decomposition() const { return decomp_; }

    // Nodes of the layout where a side runs straight past a third arc, and
    // where they are.
    const std::vector<int> &tJunctionNodes() const { return tNodes_; }
    const std::vector<Point> &tJunctionPoints() const { return tPoints_; }
    // Per face of the layout: whether it became a block.
    const std::vector<char> &faceIsBlock() const { return faceIsBlock_; }

    // The B-rep: a face per block on its Coons patch, sharing the edges of
    // the macro arcs, plus the macro arcs no block uses as loose edges.
    const geom::Shape &shape() const { return shape_; }
    bool writeSTEP(const std::string &filename) const { return shape_.writeSTEP(filename); }
    bool writeBREP(const std::string &filename) const { return shape_.writeBREP(filename); }

private:
    void classifyFaces(const QuadLayout &layout);
    void buildMacroArcs(const QuadLayout &layout);
    void buildBlocks(const QuadLayout &layout);
    void fitArcs();
    void buildPatches();
    void buildShape();
    void buildDecomposition(const Mesh &mesh);

    Options opts_;
    Report report_;
    double meanEdge_ = 1.0;

    std::vector<Point> nodePos_;          // per layout node
    std::vector<char> nodeOnBoundary_;
    std::vector<int> macroOfArc_;         // layout arc -> macro arc
    std::vector<MacroArc> arcs_;
    std::vector<Block> blocks_;
    std::vector<char> faceIsBlock_;
    std::vector<int> tNodes_;
    std::vector<Point> tPoints_;
    BlockDecomposition decomp_;
    geom::Shape shape_;
};

#endif // __LAYOUT_BLOCKS_HXX__
