#ifndef __QUAD_LAYOUT_HXX__
#define __QUAD_LAYOUT_HXX__

#include <string>
#include <vector>

#include "tracing/SeparatrixTrace.hxx"

// ---------------------------------------------------------------------------
// The quad layout with T-junctions of Viertel, Osting and Staten, IMR 2019,
// Sec. 3: the partition of the model cut out by the traced separatrices, which
// is step 2 of their Algorithm 1. The simplification of Sec. 4 is deliberately
// not here.
//
// The separatrices and the boundary together form a plane graph and the
// components of the partition are its faces. Building it is three steps:
//
//   1. Nodes. Singularities, boundary corners, the points where a separatrix
//      runs into the boundary, the points where two separatrices cross, and the
//      ends of separatrices that were cut off in the interior. Crossings and
//      cuts were found while tracing, so nothing is searched for twice.
//
//   2. Arcs. Each separatrix is cut at the nodes lying on it, in order along
//      it, and each boundary loop at the nodes lying on it. An arc is a
//      polyline between two nodes with no node in between.
//
//   3. Faces. Sort the arcs around each node by the direction they leave in and
//      walk: having arrived at a node, leave by the arc one step clockwise from
//      the one just arrived along. Every such walk closes, and the walks that
//      close the wrong way round are the unbounded face and are dropped.
//
// A component of a proper quad layout has four sides, and a side is one arc or
// several when T-junctions land on it. So the count that decides whether a face
// is a quad is not how many arcs bound it but at how many of its nodes it
// actually turns a corner -- about ninety degrees rather than about straight
// on. That is what `corners` reports, and the layout is valid when every face
// has four.
// ---------------------------------------------------------------------------

class QuadLayout {
public:
    enum class NodeKind {
        Singularity,     // an irregular node of the layout, valence 3 or 5
        BoundaryCorner,  // a corner of the model
        BoundaryHit,     // where a separatrix runs into the boundary
        Crossing,        // where two separatrices cross
        Heteroclinic,    // where two separatrices were joined head-on into one curve
        TJunction,       // where a separatrix stopped on another one
        Dangling,        // a free end: a separatrix that stopped on nothing at all
        InterfaceNode,   // a node of the material interface network: a corner of its regions
        InterfaceHit     // where a separatrix crosses a material interface
    };

    struct Node {
        Point pos{0.0, 0.0};
        NodeKind kind = NodeKind::Crossing;
        int source = -1;             // singularity / corner / separatrix index, by kind
        std::vector<int> darts;      // outgoing darts, sorted counter-clockwise
        std::vector<double> angles;  // the direction each of them leaves in
    };

    // A dart is an arc with a direction: 2 * arc + 0 runs a -> b, 2 * arc + 1
    // runs b -> a. Faces are cycles of darts.
    struct Arc {
        std::vector<Point> pts;  // pts.front() sits on node a, pts.back() on node b
        int a = -1, b = -1;
        int separatrix = -1;     // -1 for a piece of the boundary or of an interface
        bool onBoundary = false;
        // A piece of a material interface: a curve of the model, like dS, and
        // so never moved or deleted by partition simplification -- but with a
        // component of the layout on both sides, which dS does not have.
        bool onInterface = false;
        double length = 0.0;
    };

    struct Face {
        std::vector<int> darts;
        std::vector<int> nodes;      // nodes[i] is where darts[i] starts
        std::vector<char> isCorner;  // whether the face turns a corner at nodes[i]
        std::vector<double> turn;    // how far it turns there, signed, in radians
        int corners = 0;
        double area = 0.0;
        // The component's sides: the runs of darts between one corner and the
        // next, in the same order the darts are in. A four-sided component has
        // four of them however many arcs each is made of -- a side carrying
        // T-junctions is several arcs and still one side, which is the whole
        // point of allowing T-junctions.
        std::vector<std::vector<int>> sides;
    };

    struct Report {
        int nodes = 0;
        int arcs = 0;
        int faces = 0;
        int singularities = 0;
        int tJunctions = 0;
        int danglingEnds = 0;   // separatrices that stopped on nothing: always a defect
        int quadFaces = 0;      // faces with exactly four corners
        int badFaces = 0;
        int unboundedCycles = 0;  // should be one per connected piece of the complement
        // Places two arcs cross without a node between them. The whole face
        // walk assumes the graph is embedded, so this has to be zero: it is the
        // one check that tells a wrong layout from an ugly one.
        int arcCrossings = 0;
        // Nodes where the arcs do not meet anywhere near squarely: tangential
        // crossings, and separatrices running into the boundary well off
        // normal. Not a validity failure -- the layout is still four-sided
        // there -- but the count of places its geometry is degenerate.
        int skewNodes = 0;
        double smallestArea = 0.0;
        double largestArea = 0.0;
        double totalArea = 0.0;
        double meshArea = 0.0;  // for checking the faces actually cover the model
    };

    explicit QuadLayout(const SeparatrixTrace &trace);

    // For an arrangement handed over directly rather than read off a trace.
    // `tol` is the distance below which two nodes are the same place.
    explicit QuadLayout(double tol);

    void build();

    // Take an arrangement as given and work out its faces. Partition
    // simplification rewrites nodes and arcs and needs the faces of what it
    // has rewritten; everything from the sorting of darts onwards is the same
    // work as build() does.
    void rebuild(std::vector<Node> newNodes, std::vector<Arc> newArcs);

    // The mesh the layout was traced on; null for one handed over directly to
    // rebuild(), which keeps no reference to a model.
    const Mesh *getMesh() const { return mesh; }
    // The trace it was built from, likewise; null for one handed over directly.
    const SeparatrixTrace *getTrace() const { return trace; }

    static int arcOfDart(int dart) { return dart >> 1; }

    const std::vector<Node> &getNodes() const { return nodes; }
    const std::vector<Arc> &getArcs() const { return arcs; }
    const std::vector<Face> &getFaces() const { return faces; }
    const Report &getReport() const { return report_; }

    // The layout as one poly-line cell per arc.
    bool writeArcsVTU(const std::string &filename) const;
    // The faces as polygons, coloured by how many corners each turned out to
    // have, which is what makes a bad one easy to find in a picture.
    bool writeFacesVTU(const std::string &filename) const;
    // The nodes as vertex cells, coloured by kind.
    bool writeNodesVTU(const std::string &filename) const;

private:
    // A place along a polyline: p = (segment - 1) + t, so the path's own point i
    // sits at p = i and everything sorts on one number.
    struct Split { double p; int node; };

    int addNode(const Point &pos, NodeKind kind, int source);
    void collectNodes();
    void buildSeparatrixArcs();
    void buildBoundaryArcs();
    void buildInterfaceArcs();
    void sortAroundNodes();
    void traceFaces();
    void checkEmbedding();
    void finish();   // the tail shared by build() and rebuild()
    int addArc(std::vector<Point> pts, int a, int b, int separatrix, bool onBoundary);

    const SeparatrixTrace *trace = nullptr;
    const Mesh *mesh = nullptr;

    std::vector<Node> nodes;
    std::vector<Arc> arcs;
    std::vector<Face> faces;
    Report report_;

    std::vector<std::vector<Split>> splitsOfSeparatrix;
    // Node ids on each boundary loop, as (edge index in the loop + parameter
    // along it, node).
    std::vector<std::vector<Split>> splitsOfLoop;
    // And on each interface branch (Interfaces::Branch::verts), the same way.
    std::vector<std::vector<Split>> splitsOfBranch;

    // Position lookup for merging coincident nodes.
    std::vector<std::pair<long long, int>> nodeHash;
    double tol = 1e-9;
};

#endif // __QUAD_LAYOUT_HXX__
