#ifndef __INTERFACES_HXX__
#define __INTERFACES_HXX__

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "mesh/Mesh.hxx"

// Stage 0b: the material interface network of a multi-material domain, read as
// the feature graph the rest of the pipeline is written against.
//
// Shepherd, Gu and Hughes (2022) take features as an input to Stage 0: trim
// curves, dihedral creases and whatever the user designated, marked on the
// triangulation before anything else runs, and thereafter honoured by E3 of
// Sec. 3.3 and checked by Sec. 4's "every feature chain is a union of arcs".
// The paper never says where they come from because for its models they come
// from the CAD.
//
// Here they come from the material tags. These are planar midsurfaces meshed
// conformally across the physical surfaces of the .geo, so the interface
// between two materials is a chain of mesh edges whose two triangles disagree
// on Mesh::triangleMatId, and it is a feature in exactly the paper's sense: a
// curve of the input that the output has to keep, because an element straddling
// it would carry two materials and no analysis code can integrate that.
//
// ### Why this is a stage and not a predicate
//
// Extracting the edges is three lines and SubdomainLabels used to do it inline.
// What is not three lines, and what that inline version got wrong on every
// model in the corpus, is everything the *graph* of those edges says:
//
//   * Where it branches. Three materials meeting along a line meet at a point,
//     and that point is a singular point of any layout containing the
//     interfaces -- the layout has to turn a whole number of quarters in each
//     of the sectors between them, and the number is not free.
//
//   * Where it lands on dS. An interface reaching the boundary splits the
//     boundary vertex's fan into sectors with the same condition on each.
//
//   * Where it kinks. An interface that is continuous but not smooth -- the
//     inner face of the buttering layer in geom013 is the corpus's honest
//     example -- turns through an angle the layout must also turn, and turning
//     it anywhere other than at the kink means an element straddling it.
//
//   * Where it closes. A closed interface loop (an inclusion) has no node at
//     all, and a layout edge that is a closed curve through no singular point
//     is not a layout edge: the loop has to be cut into at least four arcs, one
//     per quarter turn its image makes, or it bounds a face with no corners.
//
// Each of those is a *node*, and each node carries a prescribed cone index that
// follows from its sectors and from nothing else -- see nodeIndex() below. That
// is the whole content of this stage: cone indices for Stage 1 that come from
// the geometry rather than from the cross field, at the places where the field
// cannot be trusted to supply them because the field does not know the
// interfaces are there.
//
// ### The index at a node, and what "ill-posed" means
//
// Around an interior node the incident branches cut the disk into sectors
// summing to 2 pi; around a boundary node the branches and the two boundary
// edges cut the fan into sectors summing to the interior angle. In a
// quadrilateral layout every sector is a whole number of quarter turns, at
// least one, and the total is the cone angle:
//
//     sum_s q_s = 4 - I(v)    interior            (cone angle 2pi - (pi/2) I)
//     sum_s q_s = 2 - I(v)    on dS               (angle    pi - (pi/2) I)
//
// so taking q_s as the nearest positive integer to sector_s / (pi/2) fixes I(v)
// with no freedom left. A junction of three interfaces at 120 degrees gives
// q = (1,1,1) and I = +1; an orthogonal T gives q = (1,1,2) and I = 0, which is
// to say a regular vertex, which is right -- a layout can run straight through
// an orthogonal T and nothing singular happens there.
//
// The residual max_s |sector_s - q_s pi/2| is the measure of how much the
// layout has to distort to make the sectors integral, and it is what the
// well-posed/ill-posed distinction in data/geometry/multimat is about. It is
// reported per node and never acted on: a large residual is a statement about
// the domain, not an error, and the domains that have one were built to have
// it.
class Interfaces {
public:
    enum class NodeKind {
        Junction,   // three or more interface branches meet inside S
        Landing,    // an interface reaches dS
        Kink,       // two branches meet at an angle the layout has to turn
        LoopSplit,  // a closed interface loop, cut so that it has corners at all
        Balance,    // a corner balance() had to put on a smooth interface
        Dangling    // one branch and nothing else: a tag error, reported
    };

    // One curve leaving a node: an interface branch, or -- at a node on dS --
    // one of the two boundary edges, which bound the fan and are sector
    // boundaries in exactly the same way.
    //
    // The two faces are the triangles immediately either side of the ray, and
    // they are what makes the sector condition writable as an energy. A tangent
    // in the image is a difference of two values of Psi, and Psi lives on Omega,
    // where the node may have been split into several children by the cutting
    // graph; naming the face says which child, and taking the face on the side
    // of the ray that the sector is on means the two tangents of one sector are
    // read at the same child whenever the cut does not run through that sector.
    struct Ray {
        int branch = -1;      // index into branches(); -1 for a boundary ray
        int neighbour = -1;   // the vertex it points at
        int edge = -1;        // the mesh edge it runs along
        int faceCCW = -1;     // the triangle counter-clockwise of it at the node
        int faceCW = -1;      // and the one clockwise of it
        double dir = 0.0;     // its direction on the model, radians
    };

    // One node of the interface graph, with the cone Stage 1 is asked to put
    // there.
    //
    // `rays` runs counter-clockwise. At an interior node it is cyclic and
    // `sector` has one entry per ray, the last wrapping round; at a node on dS
    // it starts and ends with the two boundary rays and `sector` has one fewer,
    // spanning the interior angle and no more.
    struct Node {
        int vertex = -1;
        NodeKind kind = NodeKind::Junction;
        bool onBoundary = false;

        std::vector<Ray> rays;
        std::vector<double> sector;      // between consecutive rays, CCW
        std::vector<int> quarters;       // q_s, the nearest positive integer
        std::vector<int> sectorRegion;   // which material region each sector is in
        int index = 0;                   // I(v) = (4 or 2) - sum q_s
        double residual = 0.0;           // max_s |sector_s - q_s pi/2|, radians
        bool wellPosed = true;           // residual under Options::wellPosedAngle
        // Quarters balance() moved between this node's sectors, which changes
        // where the layout turns without changing I(v).
        int rebalanced = 0;

        // Interface branches only, in the order they appear in `rays`.
        std::vector<int> branches;
    };

    // One branch: a maximal run of interface edges between two nodes, or a
    // closed loop when Options::splitLoops is off and there is no node on it.
    struct Branch {
        std::vector<int> verts;          // original mesh vertices, node to node
        std::vector<int> edges;          // original mesh edges, in the same order
        int node0 = -1, node1 = -1;      // indices into nodes(), -1 if closed
        int matLeft = 0, matRight = 0;   // the two materials it separates
        double length = 0.0;             // on the model
        double turning = 0.0;            // total signed turn along it, radians
        bool closed = false;
    };

    // One connected piece of one material: a component of S cut along the
    // interfaces. This is what has to come out of the pipeline as a whole
    // number of quadrilateral patches, and it is what balance() is written
    // against.
    struct Region {
        std::vector<int> triangles;
        int material = 0;
        int chi = 1;                     // V - E + F of its own triangles
        double area = 0.0;
        // Quarter turns its boundary makes, plus the indices of the cones
        // strictly inside it. A quadrilateral layout of a surface with boundary
        // has this equal to 4 chi; see balance().
        int turning = 0;
        int deficit = 0;                 // 4 chi - turning
        std::vector<int> branches;       // interface branches on its boundary
    };

    struct Options {
        // Two branches meeting at a degree-2 vertex are one branch through it
        // unless the direction turns by more than this, in which case the
        // vertex is a Kink node. The layout has to turn a quarter somewhere
        // along a corner and the only place it can turn one without an element
        // straddling the interface is at the corner itself.
        //
        // A quarter of a turn is pi/2; half of that is the natural split, and
        // it is loose enough not to fire on the discretisation of a smooth arc
        // (a 70-segment quarter circle turns 1.3 degrees per vertex).
        double kinkAngle = M_PI / 4.0;

        // Cut a closed interface loop that carries no node into arcs. A closed
        // curve through no singular point cannot be a union of layout edges:
        // the face inside it would have no corners. Four is the fewest that
        // works and is what an inclusion wants.
        bool splitLoops = true;
        int loopSplits = 4;

        // A node's sectors are called well-posed when every one of them is
        // within this of a whole number of quarter turns. Purely diagnostic.
        // 10 degrees: comfortably above the discretisation of a smooth
        // interface and far below the 28 and 30 degree residuals the corpus's
        // deliberately oblique junctions carry.
        double wellPosedAngle = 10.0 * M_PI / 180.0;

        // Prescribe an index at Landing and Kink nodes as well as at
        // Junctions. Off leaves the field's own reading at those vertices,
        // which is what the single-material pipeline does everywhere.
        bool prescribeLandings = true;
        bool prescribeKinks = true;

        // Balance each material region's own Gauss-Bonnet count, moving a
        // quarter turn across an existing node where that will do and putting a
        // new corner on a smooth interface where it will not. Off is the way to
        // see what the sector reading alone comes to; it is not a usable
        // setting, because E3 and E6 are then asking for a map that does not
        // exist. See balance().
        bool balanceRegions = true;
        // A cap, so that a badly tagged model reports rather than grinds. Two
        // per region is more than any model in the corpus needs.
        int maxBalanceMoves = 64;
    };

    struct Report {
        int materials = 0;
        int interfaceEdges = 0;
        int interfaceVertices = 0;
        int branches = 0;
        int closedLoops = 0;

        int nodes = 0;
        int junctions = 0, landings = 0, kinks = 0, loopSplits = 0, dangling = 0;
        int illPosedNodes = 0;
        double worstSectorResidual = 0.0;   // radians
        int worstSectorNode = -1;

        int prescribedCones = 0;            // nodes with I != 0
        int prescribedIndexSum = 0;

        // balance()
        int regions = 0;
        int regionsBalanced = 0;            // of those, deficit 0 when it finished
        int quartersMoved = 0;              // across an existing node
        int cornersInserted = 0;            // new nodes on a smooth interface
        int worstRegionDeficit = 0;         // what was left, largest in magnitude

        // Interface edges whose two triangles are the same material after all,
        // and vertices where the graph is neither a path nor a recognised node.
        int nonManifoldVertices = 0;

        std::vector<std::string> messages;
    };

    // Two forms rather than a default argument, matching the rest of MERIDIAN:
    // Options carries default member initialisers.
    explicit Interfaces(std::shared_ptr<Mesh> mesh);
    Interfaces(std::shared_ptr<Mesh> mesh, const Options &opts);

    // True when the mesh carries more than one material at all. Everything
    // below is empty when it does not, and the single-material pipeline is then
    // bit-for-bit what it was.
    bool multiMaterial() const { return report.materials > 1; }

    const std::vector<Node>& nodes() const { return nodeList; }
    const std::vector<Branch>& branches() const { return branchList; }
    const std::vector<Region>& regions() const { return regionList; }
    // Which region each triangle of the input mesh belongs to, or -1 before
    // balance() has run.
    const std::vector<int>& triangleRegion() const { return triRegion; }

    // Make each material region's own quadrilateral layout possible.
    //
    // The sector reading of quantiseNode() is local: it says what the layout
    // does where two interfaces meet, and it is right there. It says nothing
    // about a region as a whole, and a region as a whole is subject to a
    // condition of its own -- the discrete Gauss-Bonnet count of Eq. (4)
    // restricted to it,
    //
    //     sum_{v inside}  I(v)  +  sum_{v on its boundary} (2 - q_v)  =  4 chi
    //
    // -- which the sectors satisfy only by accident. Where it fails the failure
    // is not a slack constraint, it is a contradiction: geom001 is a quarter
    // disk inside a square, its interface is one smooth arc meeting the bottom
    // edge and the left edge at right angles, and E3 asks that arc to be one
    // straight line in the image while E6 asks it to leave one end vertically
    // and the other end horizontally. There is no such line. The continuation
    // does not fail to converge, it converges with det J at 1e-4 and both
    // residuals stuck at 1e-5, which is a layout-shaped object that is not a
    // layout. What is missing is one corner in the middle of the arc, and the
    // count above is what says so -- the quarter disk's boundary makes three
    // quarter turns and needs four.
    //
    // Two operations, cheaper first:
    //
    //   * move a quarter across an existing node, from the sector on the side
    //     with too many to the sector on the side with too few. I(v) is
    //     unchanged, since it is the total; what changes is where the layout
    //     turns. This is what an inclusion needs -- the four points a closed
    //     loop was cut at read as two straight-through sectors apiece, and the
    //     layout wants one quarter inside and three outside.
    //
    //   * put a new node on a branch between the two regions, with one quarter
    //     on the side that is short and three on the side that is over. It goes
    //     where the branch has turned half of whatever it turns, so that a
    //     smooth arc is cut where a coordinate line can follow both halves.
    //
    // Takes Stage 1's index per vertex, because the cones on dS are half the
    // count and they are not known until Stage 1 has run. Regions whose deficit
    // it cannot place are left alone and reported: the remedy for those is
    // upstream in Stage 1, which chose where the boundary cones went.
    void balance(const std::vector<int> &coneIndex);

    // Interface edges of the *original* mesh, as a set and as a flag per edge.
    const std::vector<int>& interfaceEdges() const { return edgeList; }
    bool isInterfaceEdge(int e) const {
        return e >= 0 && e < static_cast<int>(edgeIsInterface.size()) && edgeIsInterface[e];
    }
    // Node number at an original vertex, or -1.
    int nodeAt(int v) const {
        return (v >= 0 && v < static_cast<int>(vertNode.size())) ? vertNode[v] : -1;
    }
    // Branch number an original *vertex* lies on, or -1. A node lies on
    // several, so this is -1 there.
    int branchAt(int v) const {
        return (v >= 0 && v < static_cast<int>(vertBranch.size())) ? vertBranch[v] : -1;
    }
    // Branch number an original *edge* belongs to, or -1 if it is not an
    // interface edge. Every interface edge belongs to exactly one branch.
    int branchOfEdge(int e) const {
        return (e >= 0 && e < static_cast<int>(edgeBranch.size())) ? edgeBranch[e] : -1;
    }
    // How many interface branches leave an original vertex.
    int degreeAt(int v) const {
        return (v >= 0 && v < static_cast<int>(vertDegree.size())) ? vertDegree[v] : 0;
    }

    // (vertex, I(v)) for every node this stage asks Stage 1 to prescribe --
    // that is, every node whose kind Options allows and whose sectors are
    // determined. Nodes with I = 0 are included: a node the field reads as a
    // cone and the geometry says is regular has to be told so, or the layout
    // acquires a singularity in the middle of a straight interface.
    std::vector<std::pair<int, int>> prescription() const;

    // The nodes that are singular points of the *layout*, as the vertices they
    // sit on: the ones with a sector wider than a single quadrilateral, plus
    // any whose index makes them a cone outright.
    //
    // This is the emitter set Stages 5 and 7 trace from, and the reason it is a
    // subset rather than every node is what an emitter is also allowed to be: a
    // point a separatrix may *end* at. The layout edges at a node with sectors
    // q_1..q_k are the k interface branches and boundary edges already there
    // plus sum_s (q_s - 1) more, so a node whose every sector is exactly one
    // quarter has none to spare -- nothing to emit, and no room to receive.
    // Ending a separatrix there splits a sector that is already one
    // quadrilateral wide into two, which is not a layout; on geom001 it is what
    // pinned both of the -1 cone's flanking rays to the two points where the
    // arc meets dS and turned two of the eight faces into triangles.
    //
    // The same rule already holds for cones on dS -- one of index +1 emits
    // 1 - I = 0 rays and Separatrices will not snap a curve to it -- so this is
    // that rule, restated for a node whose ray count comes from its sectors
    // rather than from its index.
    std::vector<int> emitterNodes() const;

    // Which region a vertex lies inside, or -1 on an interface (where it is in
    // two of them) or before the regions exist.
    int regionAt(int v) const;

    // The material on each side of an interface edge, in no particular order.
    std::pair<int, int> materialsAcross(int edge) const;

    const Mesh& getMesh() const { return *mesh; }
    std::shared_ptr<Mesh> getMeshPtr() const { return mesh; }
    const Report& getReport() const { return report; }

    // The interface network as a polyline .obj, one line per branch. What the
    // stage found, in a form that can be laid over the model.
    bool writeOBJ(const std::string &filename) const;

private:
    void collectEdges();
    void findNodes();
    void buildBranches();
    void splitClosedLoops();
    void measureNodes();
    void buildRegions();
    void measureRegions(const std::vector<int> &coneIndex);
    // Put a corner on `branch`, at the point it has turned half its total, with
    // one quarter on the side of `region` and three on the other. Returns false
    // when the branch has no interior vertex left to put one on.
    bool insertCorner(int branch, int region, const std::vector<int> &coneIndex);

    // Sectors, quarters and I(v) at one node, from the directions already in
    // it. Interior nodes close on 2 pi; boundary nodes run from one boundary
    // edge to the other across the interior angle.
    void quantiseNode(Node &n) const;

    std::shared_ptr<Mesh> mesh;
    Options options;

    std::vector<int> edgeList;            // interface edges of the input mesh
    std::vector<char> edgeIsInterface;    // per input edge
    std::vector<int> vertDegree;          // interface edges at each vertex
    std::vector<int> vertNode;            // node number at a vertex, or -1
    std::vector<int> vertBranch;          // branch number at a vertex, or -1
    std::vector<int> edgeBranch;          // branch number at an interface edge, or -1
    std::vector<std::vector<int>> vertEdges;  // interface edges at each vertex

    std::vector<Node> nodeList;
    std::vector<Branch> branchList;
    std::vector<Region> regionList;
    std::vector<int> triRegion;           // per triangle of the input mesh, or -1

    // Per-sector adjustments balance() made, keyed by the node's vertex so that
    // they survive the rebuild of the node list that inserting a corner forces.
    // quantiseNode() applies them, which is what keeps a moved quarter moved.
    std::unordered_map<int, std::vector<int>> sectorAdjust;

    Report report;
};

#endif // __INTERFACES_HXX__
