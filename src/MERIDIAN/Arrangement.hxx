#ifndef __ARRANGEMENT_HXX__
#define __ARRANGEMENT_HXX__

#include <array>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "MERIDIAN/Immersion.hxx"
#include "MERIDIAN/Separatrices.hxx"
#include "MERIDIAN/SubdomainLabels.hxx"
#include "mesh/BlockDecomposition.hxx"
#include "mesh/Mesh.hxx"

// Stage 9's fit is who supplies a smooth curve for a separatrix arc; only a
// pointer to it crosses into this header, so it stays a forward declaration
// and Arrangement.cxx is the one file that needs SplineFit.hxx -- the same
// way blockDecomposition() below takes both without either header including
// the other.
class SplineFit;

// Stage 8 of Shepherd, Gu and Hughes (2022): the arrangement of the traced
// curves, and the quadrilateral layout read off it. (docs/shepherd2022.pdf,
// Sec. 4; docs/ricci_flow_pipeline.md Sec. 10.)
//
// Stage 7 handed over a bundle of curves. Nothing about a bundle of curves is a
// layout: the layout is the *partition of S* they cut it into, and that is a
// planar subdivision -- nodes, arcs, faces -- which has to be built. The recipe
// is short and the whole of Sec. 4 fits in it:
//
//     NODES = cones u boundary corners
//           u separatrix x separatrix crossings
//           u separatrix x boundary hits
//     split every curve at NODES                          -> ARCS
//     at each node sort the outgoing arcs by direction
//     twin(a) = the reversed arc
//     next(a) = the cyclic predecessor of twin(a) at its origin
//     follow next() to enumerate the faces
//
// ### Where the arrangement lives
//
// On S, not in the image. Psi is locally injective and globally not -- that is
// Q1, and the self-overlap is the magenta hatching of Figs. 10b and 10c -- so
// two sheets of the image lie on top of each other and a 2-D arrangement
// computed there would invent intersections between curves that are nowhere
// near each other on the model. Sec. 4 says as much and prescribes the remedy:
// compute the crossings **triangle-locally**. A separatrix is stored as a list
// of (face, barycentric entry, barycentric exit), so two curves cross where two
// of those segments cross inside one shared face, and that question is asked in
// barycentric coordinates, which are the same numbers in both spaces. No global
// arrangement is ever formed.
//
// The angles, on the other hand, are read in the image and only there. "The
// incident arc directions turn by pi/2" is a statement about Psi: on S the same
// corner is whatever angle the triangulation happens to make. So every node
// carries a *layout angle* per incident arc, accumulated around the node's fan
// of faces using Psi, and the sector between two consecutive arcs is what gets
// compared against a quarter turn.
//
// Accumulating rather than taking atan2 of the image direction matters at a
// cone and nowhere else. A cone of index I has cone angle 2pi - (pi/2) I, which
// for I < 0 is more than a full turn, so two of its arcs can leave along the
// same image direction and only an accumulated angle tells them apart. It is
// the same argument, for the same reason, as the fan sweep of Stage 7.
//
// ### What a face is, and which faces are patches
//
// Following next() enumerates every face of the subdivision of the *plane*,
// which is more than the patches: the unbounded face is one, and so is each
// hole of the model. Both are excluded by the same test, and it is exact rather
// than a point-in-polygon guess -- a boundary arc of the arrangement lies on dS,
// which has the model on exactly one side, so the half-edge running the other
// way belongs to a face that is not part of S. A face holding such a half-edge
// is not a patch.
//
// ### The validation, and what each check is really asking
//
// Sec. 4's list, with the reason each one is worth running:
//
//   * Every face has exactly **four** corners, a corner being a node where the
//     boundary of the face turns by a quarter. A turn of 0 is a T-junction and
//     says the quantisation has not converged; a turn of a half or three
//     quarters says a separatrix is missing or went to the wrong cone.
//   * Every arc bounds exactly two faces, or one face and dS.
//   * Every cone has as many incident arcs as its index prescribes: 4 - I in
//     the interior, 3 - I on dS. Fewer means a separatrix that Stage 7 could
//     not trace, or two that were snapped to the same place.
//   * Every feature chain is a union of arcs -- that is what E3 was for.
//
// Failure is not repaired here. Sec. 4's own remedy is upstream and is already
// implemented: raise lambda_5, add the Gamma_topo constraint for the pair that
// nearly met, and re-run Stage 6 *from the current phi*. That is
// MERIDIAN::traceAndRepair, and by the time this stage runs it has already been
// applied until it stopped helping. What is left here is to say precisely which
// face failed and how, which is the thing a caller can act on.
class Arrangement {
public:
    using Align = SubdomainLabels::Align;

    enum class NodeKind {
        Cone,            // a singularity of P
        Crossing,        // two separatrices met inside a face
        BoundaryHit,     // a separatrix left transversely through dS
        BoundaryCorner,  // dS changes its label here without a cone to explain it
        InterfaceNode,   // a node of the material interface network, Stage 0b
        InterfaceHit,    // a separatrix crossed an interface
        Dangling         // a separatrix that Q5 did not close simply stopped
    };

    // Boundary and Interface arcs are both runs of mesh edges and are built the
    // same way; the difference is which faces they bound. A boundary arc has
    // the model on one side and nothing on the other, and that is how a patch
    // is told from a hole. An interface arc has a patch on both sides -- which
    // is the whole point of it, since the two are different materials and no
    // element may straddle them.
    enum class ArcKind { Separatrix, Boundary, Interface };

    // One node of the arrangement. `angle` is the layout angle of each outgoing
    // half-edge, accumulated around the node's fan in the image and therefore
    // running to `totalAngle` -- 2pi in the interior, pi on dS, and the cone
    // angle 2pi - (pi/2) I or pi - (pi/2) I at a singularity.
    struct Node {
        Point p{0.0, 0.0};       // where it is on S
        NodeKind kind = NodeKind::Crossing;

        int cone = -1;           // slot into Separatrices::coneVertices()
        int index = 0;           // I(v), at a cone
        int vertex = -1;         // the mesh vertex it sits on, or -1
        int face = -1;           // a face of S containing it
        int edge = -1;           // the mesh edge it lies on, or -1
        double t = 0.0;          // its parameter along that edge
        bool onBoundary = false;

        std::vector<int> out;    // outgoing half-edges, by ascending layout angle
        std::vector<double> angle;
        double totalAngle = 0.0;
        int valence = 0;         // what the index prescribes; 0 where nothing does
    };

    // One arc: a piece of a separatrix between two nodes, or a piece of dS
    // between two nodes. Both carry their polyline on S, which is what Stage 9
    // fits; a separatrix arc also keeps the (face, barycentric) steps it was cut
    // from, so the fit can be re-taken on a refined triangulation without
    // re-tracing.
    struct Arc {
        ArcKind kind = ArcKind::Separatrix;
        int from = -1, to = -1;

        int curve = -1;                        // separatrix id, for a traced arc
        std::vector<Separatrices::Step> steps; // clipped to this arc
        std::vector<int> verts;                // mesh vertices, for a dS arc

        // Where each end leaves from: the face of S it starts and finishes in,
        // and -- for a boundary arc -- which directed edge of dS. Both ends
        // need it, because the layout angle of an arc at a node is read in the
        // image of the face the arc is actually in and in no other.
        int faceFrom = -1, faceTo = -1;
        int bdFrom = -1, bdTo = -1;
        // For an interface arc, the directed interface edge each end leaves
        // along, and which branch of the network it came from.
        int ifFrom = -1, ifTo = -1;
        int branch = -1;
        int matLeft = 0, matRight = 0;

        std::vector<Point> points;             // the polyline on S, from -> to
        double modelLength = 0.0;
        double imageLength = 0.0;

        // Which coordinate the arc holds constant. Read off the direction of
        // travel at the `from` end; a separatrix that crosses an arc of
        // Gamma_Hol_k with k odd swaps the two, and `familyTurns` says so.
        Align family = Align::None;
        int dirFrom = 0, dirTo = 0;
        bool familyTurns = false;

        bool dangling = false;   // its far end is a node Q5 does not allow
        bool degenerate = false; // shorter than the merge tolerance
    };

    struct HalfEdge {
        int arc = -1;
        bool forward = true;     // runs arc.from -> arc.to
        int origin = -1;
        int twin = -1;
        int next = -1;
        int face = -1;
        double angle = 0.0;      // the layout angle at `origin`
        int slot = 0;            // its position in origin's `out`
    };

    // One face of the subdivision. `sides` groups the arcs between consecutive
    // corners: a patch in the paper's sense has four corners and one arc on
    // each side, and it is the second half of that which a T-junction breaks.
    struct Face {
        std::vector<int> half;                 // in next() order
        std::vector<int> corners;              // node ids, in the same order
        std::vector<std::vector<int>> sides;   // half-edges between corners
        std::vector<int> turns;                // quarter turns at each node visited

        double area = 0.0;       // signed, on S; negative for the unbounded face
        bool patch = false;      // lies inside S
        int material = 0;        // the material of the region it lies in, or 0
        bool mixed = false;      // sampled points found more than one material
        bool quad = false;       // exactly four corners
        bool simple = false;     // and exactly one arc per side
        int tJunctions = 0;      // nodes it passes through without turning
    };

    struct Options {
        // How near two nodes have to be, as a fraction of the diagonal of S,
        // before they are the same node. Two curves that cross exactly on a
        // mesh edge are found twice, once in each incident face, and this is
        // what makes them one crossing rather than two joined by an arc of
        // length 1e-16.
        double mergeTolerance = 1e-7;

        // How far a sector may be from a whole number of quarter turns before
        // the reading is called into question, in radians, and how far the
        // image angle at a point of dS may be from pi before that point is a
        // corner of the layout. It is not a tight tolerance and should not be:
        // what it settles is whether a turn is a quarter or a half, and those
        // are 1.57 rad apart. A sector is still *classified* by the nearest
        // multiple; what exceeding this does is get it counted, because a
        // layout whose corners are 20 degrees out is one whose quantisation has
        // not finished even though every face still has four of them.
        double cornerTolerance = 0.35;

        // Keep the arcs of separatrices that Q5 did not close -- the ones that
        // ran out at the step cap or onto a closed orbit. They leave a dangling
        // node in the middle of a face, which is a defect and is reported as
        // one, but dropping them silently would hide it and would also merge
        // two faces that the finished layout keeps apart. Off is the way to see
        // what the layout looks like without them.
        bool keepUnresolved = true;

        // End such a curve at the first layout edge it meets rather than
        // wherever the step cap or the cycle test happened to stop it.
        //
        // Away from the cones Psi is flat and a separatrix is one of its
        // geodesics, so a curve Q5 does not close does not run off anywhere: it
        // winds, and where it stops is an artefact of the cap. The part of it
        // that means something is the part before it first crossed a curve that
        // *is* closed, because up to there it was cutting a face in two and
        // after there it is only cutting pieces of a face it has already left.
        // Ending it at that crossing turns a dangling node -- an edge with a
        // loose end, which is not a subdivision of S at all and which nothing
        // downstream can mesh -- into a T-junction, which is a defect the paper
        // names, which is reported as one, and which a quad mesher can at least
        // resolve with a transition.
        //
        // It does not make Q5 hold and is not meant to; Sec. 3.3's repair,
        // which MERIDIAN::traceAndRepair runs to exhaustion first, is what does
        // that. This is what to do with what is left.
        bool trimUnresolvedAtCrossings = true;

        // Add a node at any point of dS where the boundary turns a quarter in
        // the image without a cone there to account for it. Q3 puts each
        // component of dS - G on a coordinate line and Stage 3 drove the
        // boundary's curvature to zero away from the cones, so a regular
        // boundary vertex has image angle pi exactly; anything else is a corner
        // of the layout that Stage 1's cone set does not know about. Better to
        // carry it as a node and report it than to let a patch come out with
        // three corners and no explanation.
        bool boundaryCornerNodes = true;

        // An arc shorter than this fraction of the diagonal of S is not an arc
        // of the layout, and its two nodes are one node. Two things make them:
        // Stage 7 snapped a curve's end onto its cone, so the curve's own
        // polyline still passes a hair to one side of it and anything crossing
        // there crosses beside the cone rather than at it; and two crossings
        // that fall either side of a triangulation edge can survive the merge
        // tolerance. Both leave a sliver arc, and a sliver arc puts a node in
        // the middle of a patch side, which costs the patch its Coons fit for
        // no reason the layout knows about. The tolerance sits well above the
        // snap gap Stage 7 reports and well below any arc that means anything.
        //
        // Only slivers whose collapse cannot change what the model is are
        // taken: two nodes of dS joined by a piece of dS, and two interior
        // nodes joined by a piece of a separatrix. See collapseShortArcs() for
        // what the other two cases do to a face. Set to 0 to keep every arc.
        double collapseTolerance = 1e-4;

        // How near two curves have to be, as a fraction of their own length,
        // before they are the same curve traced from both ends. Every arc of a
        // finished layout that joins two cones is traced twice -- once from
        // each -- because Q5 makes it a separatrix of both, and the second copy
        // is not a second layout edge. It has to be recognised as the same one
        // or every cone comes out with twice the arcs its index prescribes and
        // the pair bounds a face of zero area between them.
        double duplicateTolerance = 1e-2;

        // Check that every feature chain is a union of arcs, and how near an
        // arc has to pass a feature edge to count as covering it, as a fraction
        // of the diagonal of S.
        bool checkFeatures = true;
        double featureTolerance = 1e-3;
    };

    struct Report {
        int nodes = 0;
        int coneNodes = 0, crossingNodes = 0, boundaryHitNodes = 0;
        int boundaryCornerNodes = 0, danglingNodes = 0;
        int mergedNodes = 0;         // crossings that turned out to coincide
        int isolatedCones = 0;       // cones with no arc at all

        int arcs = 0, separatrixArcs = 0, boundaryArcs = 0, interfaceArcs = 0;
        int interfaceNodes = 0, interfaceHitNodes = 0;
        // Patches whose interior is not all one material. Zero is what an
        // interface-respecting layout looks like and is the property the whole
        // multi-material extension exists to produce.
        int mixedPatches = 0;
        int degenerateArcs = 0;
        int duplicateCurves = 0;     // traced from both ends; the second dropped
        int collapsedArcs = 0;       // slivers whose two ends were one node
        int sliversKept = 0;         // slivers a collapse would have pinched
        int clusteredCones = 0;      // two cones a sliver apart, left alone
        int selfCrossings = 0;       // a separatrix that crossed itself
        int parallelOverlaps = 0;    // two segments too nearly collinear to cut
        // Curves Q5 did not close that were ended at a layout edge instead of
        // being left with a loose end. See Options::trimUnresolvedAtCrossings.
        int trimmedCurves = 0;

        int faces = 0;
        int patches = 0;             // faces inside S
        int quads = 0;               // of those, with exactly four corners
        int simpleQuads = 0;         // ... and one arc per side
        int annularFaces = 0;        // a face whose boundary is not one cycle
        // A face with dS on both sides of it: it is partly inside the model and
        // partly outside, which a planar subdivision cannot be. It means two
        // arcs were sorted the wrong way round at some node -- in practice two
        // that leave it along very nearly the same direction, which is a layout
        // whose quantisation has not converged rather than an arrangement that
        // could have been built differently.
        int mixedBoundaryFaces = 0;

        // Sec. 4's validation, as counts of what failed.
        int tJunctions = 0;          // face corners where the turn is zero
        int wrongCornerFaces = 0;    // patches without exactly four corners
        int coneValenceErrors = 0;   // cones whose degree is not 4 - I / 3 - I
        int danglingArcs = 0;
        int unsharedArcs = 0;        // a separatrix arc not between two patches

        // How far the sectors are from whole quarter turns -- the same question
        // Stage 6 answered as a residual, asked again of the finished layout.
        double maxSectorResidual = 0.0;
        int worstSectorNode = -1;
        int ambiguousSectors = 0;    // farther from a quarter turn than the tolerance

        // How far a curve's snapped end had to be moved to sit on its cone,
        // measured on S rather than in the image, relative to the diagonal.
        double maxConeSnapGap = 0.0;
        // Nodes whose fan could not place an arc exactly, so its layout angle
        // was interpolated across the fan instead. Non-zero only when Stage 7
        // snapped a curve to a cone from outside that cone's one-ring.
        int interpolatedAngles = 0;

        double domainArea = 0.0;
        double patchArea = 0.0;      // the patches together
        double areaCoverage = 0.0;   // patchArea / domainArea; 1 when they tile S
        double minPatchArea = 0.0, maxPatchArea = 0.0;

        int featureChains = 0, featureChainsCovered = 0;
        double maxFeatureGap = 0.0;

        // Every patch is a quadrilateral with one arc per side, every cone has
        // the arcs its index asks for, every arc is shared, and the patches
        // tile S. This is Sec. 4's validation passing.
        bool valid = false;

        std::vector<std::string> messages;
    };

    Arrangement(const Separatrices &separatrices, const SubdomainLabels &labels);
    Arrangement(const Separatrices &separatrices, const SubdomainLabels &labels,
                const Options &opts);

    // A layout put together rather than traced: TORSION's per-material mode
    // (MaterialLayout) lays each material region of S out on its own, and what
    // comes back from each is an arrangement of that region alone. They are
    // glued here -- the nodes two regions share on an interface made one node,
    // the two copies of an interface arc between them made one arc with a
    // patch on either side -- by the caller, which is the only party that
    // knows which node of one region is which node of the next. What arrives
    // is therefore already a subdivision: every face with its half-edge cycle,
    // its corners, its sides and its turns, every half-edge with its twin and
    // its successor, the two half-edges of arc a at 2a (running from -> to)
    // and 2a + 1. This constructor measures S, collects the patches and runs
    // check() on them, the same validation a traced arrangement gets, and
    // nothing else: there is no map behind the result, so hasMap() is false
    // and getSeparatrices(), getLabels() and getImmersion() throw. Stages 9 to
    // 11 read none of the three.
    struct Assembly {
        std::vector<Node> nodes;
        std::vector<Arc> arcs;
        std::vector<HalfEdge> halves;
        std::vector<Face> faces;
    };
    Arrangement(const Mesh &mesh, Assembly parts, const Options &opts);

    const std::vector<Node>& getNodes() const { return nodes; }
    const std::vector<Arc>& getArcs() const { return arcs; }
    const std::vector<HalfEdge>& getHalfEdges() const { return halves; }
    const std::vector<Face>& getFaces() const { return faces; }
    const Report& getReport() const { return report; }
    const Options& getOptions() const { return options; }

    // Whether the arrangement was traced on a map of S (the two-argument
    // constructors) or assembled from pieces laid out on their own.
    bool hasMap() const { return sep != nullptr; }
    const Separatrices& getSeparatrices() const {
        if (!sep) throw std::logic_error("Arrangement: an assembled layout has no separatrices");
        return *sep;
    }
    const SubdomainLabels& getLabels() const {
        if (!lab) throw std::logic_error("Arrangement: an assembled layout has no labels");
        return *lab;
    }
    const Immersion& getImmersion() const {
        if (!imm) throw std::logic_error("Arrangement: an assembled layout has no immersion");
        return *imm;
    }
    const Mesh& getMesh() const { return *orig; }

    // The patches, in the order Stage 9 will fit them.
    const std::vector<int>& patchFaces() const { return patches; }

    // The four sides of a patch as arcs, in cyclic order starting at its first
    // corner, each flagged with whether the arc runs forwards along the side.
    // Empty unless the face is a simple quadrilateral.
    struct Side {
        int arc = -1;
        bool forward = true;
    };
    std::vector<Side> patchSides(int face) const;

    // The layout's quadrilaterals as the shared block-decomposition
    // representation (mesh/BlockDecomposition.hxx): only the faces Sec. 4's
    // validation already accepts as blocks -- flagged `patch` and `simple` --
    // become one, exactly the faces `patchFaces()`/`patchSides()` describe. A
    // face that did not come out four-sided or one-arc-per-side is still
    // counted in Report::wrongCornerFaces; it is simply not a block here.
    //
    // `fit` supplies the fitted spline for a side where Stage 9 fitted one,
    // sampled at `samples` steps; a side with no fit behind it -- carried
    // exactly, or fit skipped -- is the arc's own polyline. `fit` may be
    // null, for a decomposition read straight off Stage 8 before Stage 9 has
    // run. `source` is carried through unchanged, for TORSION to pass
    // "TORSION" where it shares this stage with MERIDIAN.
    BlockDecomposition blockDecomposition(const SplineFit *fit = nullptr, int samples = 16,
                                          const std::string &source = "MERIDIAN") const;

    // The arcs as polylines on S: the layout drawn on the model, which is the
    // left half of the paper's Fig. 12.
    bool writeOBJ(const std::string &filename) const;
    // The patches as closed polylines, one per face.
    bool writePatchOBJ(const std::string &filename) const;

private:
    struct Seg {
        int curve = -1;
        int step = 0;
        Point a{0.0, 0.0}, b{0.0, 0.0};
    };
    struct Event {
        double param = 0.0;   // step index plus the fraction across it
        int node = -1;
    };
    // A boundary edge of S, directed so that the model lies on its left. This
    // is the orientation everything downstream is written against: the face to
    // the left of a boundary arc is a patch, the face to its right is the
    // unbounded face or a hole.
    struct DirEdge {
        int a = -1, b = -1;   // mesh vertices, in the direction of travel
        int edge = -1;        // index into Mesh::edges
        int face = -1;        // the one triangle it belongs to
    };

    void buildFan(int v);
    // How big S is and its area: the part of buildDomain() an assembled
    // arrangement needs too.
    void measureDomain();
    void buildDomain();
    void dedupeCurves();
    void collapseShortArcs();
    void collectSegments();
    void findCrossings();
    void trimUnresolved();
    void buildEndNodes();
    void buildBoundaryNodes();
    void buildInterfaceNodes();
    void findInterfaceCrossings();
    void splitCurves();
    void buildBoundaryArcs();
    void buildInterfaceArcs();
    void classifyMaterials();
    void buildAngles();
    void linkHalfEdges();
    void extractFaces();
    void classifyFaces();
    // The patches among the classified faces, and the report's area and
    // corner counts over them: the end of classifyFaces(), which an assembled
    // arrangement arrives already past.
    void collectPatches();
    void checkFeatures();
    void check();

    // The faces of S incident to vertex v, in counter-clockwise order; for a
    // boundary vertex the walk starts at the boundary edge that leaves v with
    // the model on its left, so index 0 is always the start of the fan.
    bool vertexFan(int v, std::vector<int> &fan) const;
    // The corners of face f in counter-clockwise order, as local indices.
    std::array<int, 3> ccw(int f) const;
    // The interior angle of face f at vertex v, measured in the image.
    double imageAngle(int f, int v) const;
    // The image direction of a directed edge of dS, read in its own face.
    Point dirEdgeImage(int be) const;
    // The same for a directed interface edge, read in the face named on it.
    Point dirInterfaceImage(int ie) const;
    // The layout angle at node n of a direction leaving it into face f.
    double layoutAngle(int n, int f, const Point &imageDir, const Point &modelDir);

    int addNode(const Point &p, NodeKind kind);
    int findOrAddCrossing(const Point &p, int face);
    long long cellOf(const Point &p) const;

    const Separatrices *sep = nullptr;
    const SubdomainLabels *lab = nullptr;
    const Immersion *imm = nullptr;
    const Mesh *orig = nullptr;
    const Mesh *cut = nullptr;
    const std::vector<Point> *uv = nullptr;
    Options options;

    std::vector<Node> nodes;
    std::vector<Arc> arcs;
    std::vector<HalfEdge> halves;
    std::vector<Face> faces;
    std::vector<int> patches;

    // Scratch shared between the passes.
    std::vector<Seg> segs;
    std::vector<std::vector<int>> faceSegs;      // face of S -> indices into segs
    std::vector<std::vector<Event>> curveEvents; // per traced curve
    std::vector<char> curveDropped;              // the second copy of a curve
    std::vector<int> coneNode;                   // cone slot -> node id
    std::vector<int> vertexNode;                 // mesh vertex -> node id, or -1
    std::vector<DirEdge> dirBoundary;            // ordered dS edges, model on left
    std::vector<int> boundaryOut;                // vertex -> its outgoing dS edge
    std::vector<std::vector<Event>> boundaryEvents; // per entry of dirBoundary
    // The interface network, walked branch by branch: dirInterface holds every
    // interface edge directed the way its branch runs, branchStart says where
    // each branch's run begins, and interfaceEvents carries the separatrix
    // crossings on each edge in the same way boundaryEvents does for dS.
    const Interfaces *itf = nullptr;
    std::vector<DirEdge> dirInterface;
    std::vector<std::vector<Event>> interfaceEvents;
    std::vector<int> branchStart;                // branch -> first dirInterface
    std::vector<int> branchCount;                // ... and how many
    // Gamma_u / Gamma_v per edge of dS, carried over from Stage 5.
    std::unordered_map<MeshEdgeKey, int, MeshEdgeKeyHash> edgeLabel;
    std::vector<char> fanBuilt;                  // per vertex
    std::vector<std::vector<int>> fanCache;      // its faces, counter-clockwise
    std::vector<std::vector<double>> fanOffset;  // accumulated image angle
    std::vector<double> fanTotal;

    // The reference direction a boundary node measures its layout angles from:
    // forward along dS, so that the sector inside the model is [0, totalAngle].
    std::vector<Point> nodeRef;
    // Every node, bucketed by position, so that a crossing found twice -- once
    // in each face sharing the edge it landed on -- is found to be one node.
    std::unordered_map<long long, std::vector<int>> nodeGrid;
    double cellSize = 1.0;

    double modelExtent = 1.0;
    double mergeTol = 0.0;

    Report report;
};

#endif // __ARRANGEMENT_HXX__
