#ifndef __MESH_PILLOW_HXX__
#define __MESH_PILLOW_HXX__

#include <string>
#include <vector>
#include "QuadMesh.hxx"

// Pillowing for mesh::QuadMesh: one layer of quads inserted along a stretch of
// boundary or material interface, so that no element is left with two sides
// on a feature meeting at a flat node.
//
// ## The defect
//
// A quad whose two sides at a corner both lie on features, with the corner
// opening to about 180 degrees: a block corner the layout put on a straight
// stretch of dS or of an interface. Stage 1's flat +1 cones are where these
// come from (docs/cf_flow_pipeline.md Sec. 15.9); a +1 sitting on a straight
// feature is a patch corner of pi, and Stage 10 meshes that patch with a single
// element there. The element is a triangle with a fourth node on one side.
//
// TMOP cannot repair it, whatever metric or quadrature it runs. The corner node
// is a feature node -- Fixed, or Sliding along the very line both sides lie on --
// so every admissible move keeps the two sides collinear and the corner at 180
// degrees. TMOP's frozen-corner rule (TMOP::findFrozenCorners) exists so that
// such a corner neither blocks the smoothing of everything else nor drags its
// neighbours onto itself. It leaves the corner as it found it, with a scaled
// Jacobian of zero.
//
// ## The repair, and why it is a whole layer
//
// The node has one element where its angle wants two. The only operation on a
// conforming quad mesh that changes how many elements meet at an existing node
// is pillowing (Mitchell and Tautges, "Pillowing doublets", 1995): duplicate
// the nodes of a chain of edges, give the elements on one side the copies, and
// fill the gap with a new quad per chain edge. Splitting a chord, or any other
// refinement, adds nodes but never changes the valence of an old one. After a
// pillow along the feature through the flat node, the node has two new quads
// at about 90 degrees each. Its copy is an interior node with three -- the
// +1 the layout put on the feature, moved half an element inside, where TMOP
// is free to move it in both directions.
//
// A layer of quads is a dual chord, and a chord cannot stop in the middle of
// the mesh: it closes up or both ends leave through dS. That fixes how far the
// layer has to run. From the flat node it follows the feature both ways, on the
// side the defect is on, through every node where the feature is roughly
// straight (Options::passAngle). Each such node gets two new quads and its copy
// gets the old ones. That is a straight split of the boundary row where the
// node was regular, and it moves a defect into the interior where it was not.
// At a corner the layer stops, and its copy goes onto the first edge of the
// corner's fan: the next feature side when one element fills the corner, or
// the edge between its first two elements when several do. That edge is split
// there. If dS lies beyond the edge, the chord ends. Otherwise -- an interface,
// or an interior edge at a reflex corner -- the chord carries on as a split of
// the quads beyond, entering each through one side and leaving through the
// opposite one, until it reaches dS or meets itself. Every node keeps its
// valence through all of this except the flat ones the layer passed.
//
// The alternative is to close the layer through the interior around a few
// elements. That keeps it local, but every place it turns off the feature
// becomes a feature node with three elements at 60 degrees and an interior node
// with three. Ending at corners leaves no such pair anywhere, at the price of a
// row of elements split in two along the chord.
//
// A closed feature loop with no corner on it -- an inclusion's rim -- is
// pillowed all the way round and needs no ending.
//
// A layer has to start valid. No element may be turned inside out -- the
// rebuild would re-wind it against its neighbours -- and no corner folded that
// was not folded before; elements Stage 10 already left folded, and the pieces
// a chord cuts them into, are not held against it. Half and a quarter of the
// depth are tried before a corner is given up on, and the report says where
// and why. A chord carried past an end that runs round a closed loop of
// elements, never reaching dS, is the other way a corner is given up on.
//
// ## What it does not touch
//
// A corner quad at a convex feature corner of about 90 degrees also has two
// sides on features. That is the ideal there. A layer around the corner would
// give the corner two elements at 45 degrees, so the walk ends at such corners
// and never pillows across one. Acute corners are left alone for the same
// reason; a second element would halve an angle that is already too small.
//
// The mesh is changed in place and its topology rebuilt (QuadMesh::buildTopology),
// so the feature curves, node types and targets of a smoother run afterwards
// are those of the pillowed mesh. Existing vertices keep their indices and
// their positions; new ones are appended, so anything that indexes the old
// vertex array -- the viewer's block walls -- still finds them.
namespace mesh {

class Pillow {
public:
    struct Options {
        // A corner whose two sides both lie on features is a defect when it
        // opens wider than this, in degrees. 150 is a scaled Jacobian of 0.5
        // at that corner; anything flatter is what TMOP cannot improve.
        double flatAngle = 150.0;

        // The walk along the feature passes through a node whose fan, on the
        // pillowed side, spans 180 degrees give or take this much -- the range
        // in which two elements is the right number there -- and ends at the
        // first node outside it. A node with one element and an angle past
        // flatAngle is passed through as well, since that is a defect too.
        double passAngle = 45.0;

        // Where a copy is put: this fraction of the way from the feature node
        // to the mean centroid of the elements it takes over. At 1 a regular
        // node's copy lands half an element in, so the boundary row is split
        // down the middle. At a flat corner the element is a triangle and its
        // centroid sits a quarter of the way in, so the copy starts below its
        // neighbours' copies. That is the convex side, and the start is valid.
        double depth = 1.0;

        // Where the copy of an end node goes on its exit edge, and where a
        // continuing chord crosses each edge after that, as a fraction of the
        // edge from the side the layer is on.
        double exitFraction = 0.5;

        // At most this many layers on one mesh. Each one fixes every defect it
        // passes, so this is far more than the corpus needs.
        int maxLayers = 32;

        // Accept a layer that folds an element corner Stage 10 had not folded,
        // at the first depth that turns no element inside out, and leave the
        // fold to TMOP's untangler. Off, such a layer is tried thinner and then
        // refused, and the corner stays flat.
        bool allowFolds = false;
    };

    struct Report {
        bool ran = false;

        int defectsBefore = 0, defectsAfter = 0;
        // The flattest defect corner, in degrees, before and after (0 when
        // there are none).
        double flattestBefore = 0.0, flattestAfter = 0.0;

        int layers = 0;            // layers inserted
        int closedLayers = 0;      // of which went round a closed feature loop
        int skipped = 0;           // flat corners left that no layer could be built for
        int layerQuads = 0;        // quads added along the features
        int splitQuads = 0;        // quads split by chords carried past an end
        int featureSplits = 0;     // feature edges a chord ended or crossed on

        int quadsBefore = 0, quadsAfter = 0;
        int verticesBefore = 0, verticesAfter = 0;

        // Over every corner of every element, before and after.
        double minScaledJacobianBefore = 0.0, minScaledJacobianAfter = 0.0;
        int invertedBefore = 0, invertedAfter = 0;

        std::vector<std::string> messages;
    };

    // One defect corner: corner `corner` of quad `quad`, at `vertex`, whose two
    // sides lie on features and open to `angle` degrees.
    struct Defect {
        int quad = -1;
        int corner = -1;
        int vertex = -1;
        double angle = 0.0;
    };

    explicit Pillow(QuadMesh &mesh);
    Pillow(QuadMesh &mesh, const Options &opts);

    // Pillow every defect corner on the mesh, one layer at a time, rebuilding
    // the topology after each. Returns true when at least one layer went in.
    bool run();

    const Report &getReport() const { return report; }
    const Options &getOptions() const { return options; }

    // Every defect corner on the mesh as it is now, by quad and corner. What
    // run() works through, and what a caller can count without changing
    // anything.
    static std::vector<Defect> findDefects(const QuadMesh &mesh, double flatAngle);

private:
    // The elements around a feature node on one side of the feature: from the
    // feature edge the walk arrived on, turning away from it until the next
    // feature edge.
    struct Fan {
        std::vector<int> quads;    // in turning order
        std::vector<int> spokes;   // far end of the edge after each quad
        int next = -1;             // far end of the feature edge that closes it
        double angle = 0.0;        // radians, measured through the fan
    };

    // A node of the layer's path. `fan` holds the elements the copy takes over:
    // all of them where the layer passes through, the first only at an end.
    struct PathNode {
        int v = -1;
        std::vector<int> fan;
        bool end = false;
        int exitTo = -1;           // at an end: far node of the exit edge
    };

    bool turn(int u, int from, int start, Fan &fan) const;
    // Walk from v along the feature edge (v, first), whose element on the
    // pillowed side is `startQuad`. `nodes` gets the nodes after v, `sideQuads`
    // the pillowed-side element of each edge walked. False when the walk had to
    // give up; `closed` when it came back round to v.
    bool walk(int v, int first, int startQuad, std::vector<PathNode> &nodes,
              std::vector<int> &sideQuads, bool &closed) const;
    // A point on the edge (a, b) at fraction s from a: on the feature curve
    // the edge belongs to when it is a feature edge and has one, on the chord
    // otherwise.
    Point onEdge(int a, int b, double s, bool feature) const;
    // Build and insert one layer through the defect. False when no layer
    // could be built or none of the depths tried gave a valid start.
    bool insertLayer(const Defect &d);

    QuadMesh &mesh;
    Options options;
    Report report;
    // Why the last layer that could not be built was refused, for the report.
    mutable std::string refusal;
};

}  // namespace mesh

#endif  // __MESH_PILLOW_HXX__
