#ifndef __MESH_QUADMESH_HXX__
#define __MESH_QUADMESH_HXX__

#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <cstdint>
#include "Mesh.hxx"
#include "geom/BSpline.hxx"

// A standalone quadrilateral mesh, the quad counterpart of `Mesh`: positions,
// cells, the topology derived from them, and the per-node bookkeeping a
// variational smoother needs. It is deliberately *not* the Stage 10 class of
// the same name in src/MERIDIAN/QuadMesh.hxx -- that one is a stage of the
// layout pipeline and owns a SplineFit, an arrangement, chords and intervals;
// this one owns nothing but a mesh and can be built from an .obj, from arrays,
// or from the elements the pipeline produced. It lives in `namespace mesh` so
// both can be included in the same translation unit; refer to it as
// `mesh::QuadMesh`.
//
// It holds no smoothing of any kind. Everything here is data: what a TMOP
// (target-matrix optimization paradigm) smoother has to look up per node and
// per element, precomputed once so the optimizer's inner loop is index
// arithmetic rather than searching.
//
// ## Conventions
//
// A quad is four vertex indices in **counter-clockwise** order, so its signed
// area is positive and every corner Jacobian determinant is positive on a
// convex element. Side `i` of a quad runs from corner `i` to corner
// `(i+1) % 4`, and `quadAdjacency[q][i]` is the quad across that side, so
// sides, edges and neighbours all share one index.
//
// The reference element is the **unit square [0,1]^2**, with corner `c` at
//
//     (0,0), (1,0), (1,1), (0,1)   for c = 0, 1, 2, 3
//
// and the bilinear map x(xi, eta) = sum_c N_c(xi, eta) x_c. This matches
// MFEM's quad reference element, so target matrices and quadrature rules
// carry over without a change of variables.
//
// ## What a TMOP smoother needs, and where it is here
//
// | need | member |
// |---|---|
// | element Jacobian A at a quadrature point | `jacobianAt`, `cornerJacobian` |
// | target matrix W per element | `targetJacobian`, `targetAt` |
// | which nodes an element touches | `quads[q]` |
// | which elements a node touches, and *at which corner* | `vertexQuads` |
// | which nodes may move, and how | `nodeType`, `slideTangent` |
// | the curve a feature node is confined to | `featureCurves`, `curveOf`, `curveParam` |
// | a step legal for that node | `projectStep` |
// | one-ring stencil for a local/Gauss-Seidel sweep | `vertexNeighbors` |
// | current quality, to accept or reject a step | `scaledJacobian`, `quality` |
//
// The corner index in `vertexQuads` is the load-bearing one. A TMOP gradient
// contribution needs dA/dx for the *specific* node, which depends on which
// corner of the element that node occupies; storing it alongside the element
// index removes a four-way search from the innermost loop.
class FeatureFrame;

namespace mesh {

typedef std::array<int, 4> Quad;

// A 2x2 Jacobian, row-major: {dx/dxi, dx/deta, dy/dxi, dy/deta}. Column 0 is
// the image of the reference xi direction, column 1 of eta -- the same layout
// TMOP's A and W matrices use, so a target can be written straight in.
typedef std::array<double, 4> Jacobian2;

inline double det2(const Jacobian2 &J) { return J[0] * J[3] - J[1] * J[2]; }

// Frobenius norm, |J|_F -- the quantity most TMOP shape metrics are built from.
inline double frob2(const Jacobian2 &J) {
    return std::sqrt(J[0] * J[0] + J[1] * J[1] + J[2] * J[2] + J[3] * J[3]);
}

inline Jacobian2 mul2(const Jacobian2 &A, const Jacobian2 &B) {
    return { A[0] * B[0] + A[1] * B[2], A[0] * B[1] + A[1] * B[3],
             A[2] * B[0] + A[3] * B[2], A[2] * B[1] + A[3] * B[3] };
}

// Inverse; returns the identity and sets `ok` false when singular, so a caller
// forming T = A W^{-1} does not have to guard every call site itself.
inline Jacobian2 inv2(const Jacobian2 &J, bool &ok) {
    const double d = det2(J);
    if (std::fabs(d) < 1e-300) { ok = false; return {1.0, 0.0, 0.0, 1.0}; }
    ok = true;
    return { J[3] / d, -J[1] / d, -J[2] / d, J[0] / d };
}

// vertex -> incident quads, with the corner of the quad the vertex occupies.
// The counterpart of VertexTriangleCSR, plus that corner. Entries for vertex
// `v` are colIdx[rowPtr[v] .. rowPtr[v+1]-1] and corner[] runs in lockstep.
//
// The quads around an interior vertex are ordered counter-clockwise, walked
// side to side through `quadAdjacency`; around a boundary vertex the walk is
// rewound to the clockwise-most quad first, so the list starts on one boundary
// side and ends on the other. A vertex whose one-ring is non-manifold gets its
// incident quads in index order and `manifoldRing[v]` false -- the ordering is
// what breaks there, not the membership, so a smoother that only needs the set
// can still use it.
struct VertexQuadCSR {
    std::vector<int> rowPtr;      // size nVertices + 1
    std::vector<int> colIdx;      // quad indices
    std::vector<int> corner;      // corner 0..3 of that quad, same length as colIdx

    inline int vertexDegree(int v) const { return rowPtr[v + 1] - rowPtr[v]; }
    inline int begin(int v) const { return rowPtr[v]; }
    inline int end(int v) const { return rowPtr[v + 1]; }
};

class QuadMesh {
public:
    // How a node may move. A TMOP solve assembles only over Free and Sliding
    // nodes; Sliding ones project their step onto `slideTangent` so the node
    // stays on the curve it was placed on.
    enum NodeType : std::uint8_t {
        NodeFree    = 0,  // interior, moves in both directions
        NodeSliding = 1,  // on a smooth stretch of boundary or interface
        NodeFixed   = 2   // corner, feature junction, or pinned by the caller
    };

    struct Options {
        // Boundary and interface nodes are classified by the turn the feature
        // makes at them: a node whose two incident feature edges meet at less
        // than (180 - cornerAngle) degrees is a corner and is fixed, anything
        // straighter slides. 45 degrees is deliberately generous -- a node on a
        // discretised circle turns by 360/nSeg, which is 15 degrees on the
        // interface discretisation data/geometry/multimat/bubbles.geo asks for,
        // and none of those may be mistaken for corners.
        double cornerAngle = 45.0;

        // Pin every feature node instead of letting the smooth ones slide.
        //
        // A debug and isolation knob, not the setting to start from. Pinning
        // the whole boundary freezes its discretisation density at whatever an
        // earlier stage left it while the interior metric pulls towards its own
        // ideal spacing, and the mismatch is absorbed by the one ring of
        // elements touching the boundary -- the "bubbly" boundary ring. Sliding
        // (the default) lets that ring redistribute instead, and it is exact
        // rather than approximate here: dS and the material interfaces are
        // never splined, so the chord through a node's two feature neighbours
        // *is* the polyline the node was placed on. Turn this on to compare a
        // suspect result against a fully pinned one and tell whether a
        // regression came from the sliding or from the interior metric.
        bool fixAllFeatureNodes = false;

        // Treat material interfaces as features. Off, only dS is a feature and
        // a smoother is free to drag an element across an interface, which
        // destroys the per-element material assignment the mesh carries.
        bool interfacesAreFeatures = true;

        // What a sliding node slides *on*.
        //
        // The three answers differ in where the node is after a step, and the
        // difference is the whole of this setting:
        //
        //   * Chord takes the direction of the chord through the node's two
        //     feature neighbours and moves along it. The node is then off the
        //     feature by that chord's sagitta -- on a straight run, nothing; on
        //     a discretised circle, r(1 - cos(theta/2)) *inwards*, every step,
        //     never recovered. A rim of such nodes erodes towards its own
        //     chords, which is the polygonalisation seen on ATLAS's meshes of
        //     the disk models.
        //   * Polyline binds each smooth run of the feature to its own
        //     polyline, carried as the degree-1 B-spline geom::Polyline builds,
        //     and moves the node in that curve's parameter. The node lands on
        //     the input polyline exactly, however many steps it takes, and
        //     passes through the input's own vertices rather than cutting them.
        //     The geometry is unchanged from the mesh it was handed.
        //   * Spline binds the run to a C2 cubic interpolating its nodes
        //     (geom::interpolateCurve, so the curve passes through every one of
        //     them and binding moves nothing). Between them the curve is the
        //     smooth curve the discretisation was sampled from, so a node
        //     sliding on it rounds a discretised circle out towards the circle
        //     instead of in towards the polygon. **The default.**
        //
        // Feature corners, junctions, ends and pins are NodeFixed and are the
        // ends of these runs, so no curve is ever fitted across a corner and
        // none of them moves. A run whose interpolant bows further from its own
        // polyline than curveMaxDeviation falls back to Polyline: that is the
        // guard against a cubic ringing through a turn that cornerAngle let
        // past.
        enum CurveSource : std::uint8_t {
            CurveChord    = 0,
            CurvePolyline = 1,
            CurveSpline   = 2
        };
        CurveSource curveSource = CurveSpline;

        // How far a fitted curve may bow from the polyline of the same nodes,
        // as a fraction of that run's mean segment length, before the run falls
        // back to its polyline. A circle cut into n segments bows by about
        // pi^2 / (2 n^2) of a segment -- 0.03 at the 15 degrees per vertex
        // bubbles.geo asks for, and less on anything finer -- so 0.25 passes
        // every genuinely curved run in the corpus and still refuses the
        // overshoot a cubic makes at a 45 degree turn.
        double curveMaxDeviation = 0.25;
    };

    // Element quality, computed on demand by `computeQuality()`. Nothing here
    // is used by the class itself; it is what a smoother reports before and
    // after, and what an acceptance test reads.
    struct Quality {
        int vertices = 0;
        int quads = 0;
        int freeNodes = 0, slidingNodes = 0, fixedNodes = 0;

        double minScaledJacobian = 0.0;   // over every corner of every quad
        double meanScaledJacobian = 0.0;
        int invertedQuads = 0;            // a corner with non-positive determinant
        int nonConvexQuads = 0;           // reflex corner, but still positive area

        double minArea = 0.0, maxArea = 0.0, totalArea = 0.0;
        double minEdge = 0.0, maxEdge = 0.0, meanEdge = 0.0;

        // Aspect ratio per quad, the longer of the two mid-side spans over the
        // shorter, worst over the mesh.
        double worstAspect = 1.0;

        int nonManifoldEdges = 0;         // used by three or more quads
        int boundaryEdges = 0;
        int boundaryLoops = 0;
        bool allCounterClockwise = false;
    };

    // One smooth run of a feature, as a curve: the geometry the sliding nodes
    // on it are confined to. `chain` is the run's vertices in parameter order,
    // its first and last being the NodeFixed ends that bound it -- the same
    // vertex twice when the run is a closed loop, which is exactly the case
    // anchorFreeFeatureLoops() guarantees by fixing one node on any feature
    // loop that had none. So no curve here is periodic and no parameter wraps.
    struct FeatureCurve {
        geom::BSplineCurve<2> curve;
        std::vector<int> chain;
        bool closed = false;      // chain.front() == chain.back()
        bool fitted = false;      // an interpolating cubic, not the polyline
        double deviation = 0.0;   // how far it bows from the chain's polyline
        double length = 0.0;      // of that polyline
    };

    QuadMesh() = default;

    // Load from an .obj. Faces with four vertices become quads; a face with
    // three is refused rather than silently split, and one with more is
    // refused, because either would change the element count the caller thinks
    // it has. `usemtl mat<id>` sets the material id, as in Mesh.
    explicit QuadMesh(const std::string &filename);

    QuadMesh(const std::vector<Point> &verts, const std::vector<Quad> &cells);

    // As above with an explicit material id per quad. An empty matIds is taken
    // as a single-material mesh and fills quadMatId with 1s.
    QuadMesh(const std::vector<Point> &verts, const std::vector<Quad> &cells,
             const std::vector<int> &matIds);

    QuadMesh(const std::vector<Point> &verts, const std::vector<Quad> &cells,
             const std::vector<int> &matIds, const Options &opts);

    // Adopt the elements of anything that exposes `vertices()`, `quads()` and
    // `quadMaterials()`. Both of the pipeline's final meshes do: the Stage 10
    // `::QuadMesh` of src/MERIDIAN/QuadMesh.hxx, and the merged matrix +
    // O-grid mesh of `DiskTemplate` -- so
    //
    //     mesh::QuadMesh out = pipeline.hasDiskTemplate()
    //         ? mesh::QuadMesh::from(pipeline.getDiskTemplate())
    //         : mesh::QuadMesh::from(pipeline.getQuadMesh());
    //
    // is the whole of the conversion, for MERIDIAN and for TORSION alike.
    //
    // It is a template rather than two overloads deliberately: this header
    // must not include the pipeline's, because the pipeline's Stage 10 class
    // is also called QuadMesh and including both here would put the two names
    // in front of every user of this file. Duck typing keeps the dependency
    // pointing one way, from the pipeline to the mesh library and never back.
    //
    // The elements are taken as they are; only orientation is normalised (see
    // orientQuads) and the topology rebuilt. Material ids come across
    // unchanged, which is what the .mesh writer below turns into element
    // attributes.
    template <class Source>
    static QuadMesh from(const Source &src, const Options &opts) {
        return QuadMesh(src.vertices(), src.quads(), src.quadMaterials(), opts);
    }
    // Two overloads rather than a defaulted argument: a default of `Options{}`
    // would need Options' member initializers inside the class body, which is
    // not allowed before the class is complete.
    template <class Source>
    static QuadMesh from(const Source &src) {
        return QuadMesh(src.vertices(), src.quads(), src.quadMaterials(), Options());
    }

    // ---- the mesh itself -------------------------------------------------

    std::vector<Point> vertices;
    std::vector<Quad> quads;            // CCW, four vertex indices
    std::vector<int> quadMatId;         // one per quad; all 1 when unspecified

    // ---- topology, filled by buildTopology() -----------------------------

    // Unique undirected edges, stored (min, max) as in Mesh.
    std::vector<Edge> edges;
    // quad -> its four edges, side i joining corners i and (i+1)%4.
    std::vector<std::array<int, 4>> quadEdges;
    // edge -> the quads on it, -1 for the missing second one on a boundary
    // edge. A third user makes the edge non-manifold: it is recorded in
    // nonManifoldEdges and dropped here rather than overwriting a neighbour.
    std::vector<std::array<int, 2>> edgeQuads;
    std::vector<int> boundaryEdges;
    std::vector<bool> isBoundaryEdge;
    std::vector<int> nonManifoldEdges;

    // quad -> the quad across each of its four sides, -1 on the boundary.
    std::vector<std::array<int, 4>> quadAdjacency;

    std::vector<int> boundaryVertices;
    std::vector<bool> isBoundaryVertex;

    // Boundary loops as ordered vertex rings, each closed implicitly (the
    // first vertex is not repeated at the end). The interior is on the left of
    // the direction of travel, so the outer loop is counter-clockwise and a
    // hole is clockwise -- which also means `boundaryLoops[i]` can be handed
    // straight to a polygon area test to tell the two apart.
    std::vector<std::vector<int>> boundaryLoops;
    // vertex -> (loop, position in it), or (-1,-1) for an interior vertex.
    std::vector<std::array<int, 2>> boundaryLoopOf;

    VertexQuadCSR vertexQuads;
    // vertex -> its one-ring neighbour vertices, counter-clockwise where the
    // ring is manifold. Every vertex joined to it by a quad side appears once.
    std::vector<std::vector<int>> vertexNeighbors;
    // false where the incident quads could not be walked as a single fan.
    std::vector<bool> manifoldRing;

    // ---- features and degrees of freedom ---------------------------------

    // An edge with different material ids on its two sides. Empty in a
    // single-material mesh.
    std::vector<bool> isInterfaceEdge;
    // Boundary, interface (when Options::interfacesAreFeatures), or marked by
    // the caller through `markFeatureEdge`. A feature edge is one a smoother
    // may not move a node off.
    std::vector<bool> isFeatureEdge;
    std::vector<bool> isFeatureVertex;

    std::vector<std::uint8_t> nodeType;   // NodeFree / NodeSliding / NodeFixed
    // Unit tangent of the feature curve at a sliding node, {0,0} otherwise.
    std::vector<Point> slideTangent;

    // The two feature-curve neighbours of a vertex that has exactly two, and
    // {-1,-1} for every other vertex -- an interior node, a curve end, a
    // junction where three feature curves meet. This is what makes the tangent
    // refresh below O(1); it is also the feature graph itself, walkable as an
    // adjacency list.
    //
    // It is filled for *any* vertex with two feature neighbours, not only the
    // sliding ones, so a walk along a feature polyline passes through the
    // corners on it rather than stopping at them. The node type is what says
    // whether a vertex may move; this says only what it is joined to.
    std::vector<std::array<int, 2>> featureNeighborOf;

    // The feature curves, one per smooth run, and the binding of the sliding
    // nodes to them: curveOf[v] indexes featureCurves and is -1 for every node
    // that is not bound (interior, fixed, or Options::curveSource == CurveChord),
    // curveParam[v] is that node's parameter on its curve.
    //
    // The parameter is the node's degree of freedom, not a cache of its
    // position: a bound node is moved by changing it and evaluating the curve
    // there (setVertexParameter), which is what keeps the node on the curve to
    // the last bit however far it travels. `vertices[v]` and
    // `curve.evaluate(curveParam[v])` therefore agree by construction, and a
    // caller that writes `vertices[v]` directly is expected to rebind or to
    // stop using the curve.
    std::vector<FeatureCurve> featureCurves;
    std::vector<int> curveOf;
    std::vector<double> curveParam;

    // ---- TMOP targets ----------------------------------------------------

    // Target (ideal) Jacobian W per quad, row-major as Jacobian2. The identity
    // scaled by the target edge length is "an axis-aligned square of side h";
    // `setUniformTargets` and `setTargetsFromCurrentShape` build the two
    // ordinary cases. Empty means "no targets set", and `targetAt` then returns
    // the identity so a caller can run shape-only metrics without filling it.
    //
    // One target per element, not per quadrature point. Every target this
    // codebase has wanted so far -- an ideal square, a per-element size from
    // the layout, a shape read off the initial mesh -- is constant within an
    // element, and `targetAt(q, c)` is the only place that assumption lives:
    // making targets vary within an element means changing that accessor and
    // its storage, and nothing else.
    std::vector<Jacobian2> targetJacobian;

    Options options;
    Quality quality;

    // ---- construction ----------------------------------------------------

    // Rebuild everything derived from vertices/quads/quadMatId: edges,
    // adjacency, boundary, loops, the CSR, features, node types and the feature
    // curves. Safe to call again after the connectivity changes. A smoother
    // that only moves vertices does not need it -- no member here depends on
    // position except `slideTangent`, `quality` and the feature curves, each of
    // which has its own refresh call.
    //
    // Note that the curves are built from the positions as they are *now*: they
    // are the geometry the mesh was handed, and rebuilding them from a mesh a
    // smoother has already moved would let the boundary walk away one rebuild
    // at a time. Smooth, then rebuild only if the connectivity changed.
    void buildTopology();

    // Reorder any quad whose signed area is negative so every element is CCW.
    // Called by buildTopology(); exposed because a caller assembling quads by
    // hand may want to fix orientation before anything else reads them.
    int orientQuads();

    // Recompute slideTangent from the current positions. Cheap, and the one
    // topology-independent thing a smoother invalidates when it moves a node.
    // A thin loop over updateSlideTangent().
    void computeSlideTangents();

    // Recompute slideTangent[v] alone, in O(1): the unit tangent of v's feature
    // curve at its current parameter when it is bound to one, and otherwise the
    // chord through the current positions of v's two feature neighbours.
    //
    // A Gauss-Seidel sweep has to call this immediately before each node's own
    // local solve: the neighbours may already have moved earlier in the same
    // sweep, and a stale tangent no longer lies along the current polyline, so
    // the "the node stays exactly on the segment it was placed on" guarantee
    // that projectStep() gives would quietly stop holding.
    void updateSlideTangent(int v);

    // The displacement node `v` is actually allowed to take, given a
    // displacement a solver would like it to take: zero for NodeFixed,
    // unchanged for NodeFree, and for NodeSliding the part of it the node's
    // constraint permits.
    //
    // For a node bound to a feature curve that means travelling along the
    // curve: the tangential part of the displacement is read as an arc length,
    // turned into a parameter step by dividing by |C'(u)|, and the node goes to
    // C(u + du) -- so the returned displacement is a chord of the curve, and
    // the node lands **on** the curve rather than on its tangent. `newParam`,
    // when given, comes back with that u + du, which is what a caller applying
    // the step must write into curveParam[v]; setVertexParameter() does both.
    //
    // For an unbound sliding node (Options::CurveChord, or a run no curve could
    // be built for) it is the component along slideTangent[v], the chord through
    // the node's two feature neighbours. The node is then off the feature by
    // that chord's sagitta, but the enclosed **area** is exact: a node moving
    // parallel to the line through its neighbours leaves the triangle they span
    // unchanged, and that triangle is the polygon's only dependence on the
    // node. A node on a curve trades that exactness for staying on the
    // geometry -- the area it encloses then follows the curve's, not the
    // polyline's, which on a discretised circle is the point.
    //
    // This is the one place the mobility constraint lives, so an optimizer's
    // inner loop never re-derives it per call site. It is a post-hoc
    // projection of a full 2-D step rather than a reduced local system solved
    // in the tangent direction -- the same accept-a-displacement shape the
    // hand-rolled smoothers in DiskTemplate::smooth() and
    // TORSION::relaxToKernel() already use.
    Point projectStep(int v, const Point &displacement, double *newParam = nullptr) const;

    // ---- feature curves --------------------------------------------------

    // Cut the feature graph into its smooth runs, build a curve for each per
    // Options::curveSource, and bind the sliding nodes on it. Returns the
    // number of curves built. Called by buildTopology() after classifyNodes(),
    // and safe to call again -- it rebuilds from scratch. Does nothing but
    // clear the binding when curveSource is CurveChord.
    int buildFeatureCurves();

    // Whether v's steps are taken on a curve rather than along a chord.
    bool isOnCurve(int v) const {
        return v >= 0 && v < static_cast<int>(curveOf.size()) && curveOf[v] >= 0;
    }
    // The point of v's curve at parameter u, and the unit tangent there with
    // (optionally) |C'(u)|, the speed the parameter runs at. The node's own
    // position and {0,0} for a vertex that is not bound.
    Point curvePoint(int v, double u) const;
    Point curveTangent(int v, double u, double *speed = nullptr) const;
    // u brought into the domain of v's curve.
    double clampParameter(int v, double u) const;
    // Put a bound node at parameter u: sets curveParam[v] and vertices[v] to
    // the curve's own evaluation there, so the node sits on the curve exactly
    // rather than at the far end of an accumulated displacement.
    void setVertexParameter(int v, double u);
    // How far each bound node is from the curve it is bound to, at most. Zero
    // to rounding on a mesh nothing has moved off its curve; this is the
    // invariant a smoother is checked against.
    double maxCurveDeviation() const;

    // Classify nodes into Free / Sliding / Fixed from the feature flags and
    // Options::cornerAngle, fill featureNeighborOf, and anchor any feature loop
    // that would otherwise be free to drift (see the .cxx). Called by
    // buildTopology().
    void classifyNodes();

    // Mark an edge (by its index in `edges`) as a feature, and both its
    // vertices as feature vertices. Re-run classifyNodes() afterwards to have
    // the node types follow. Returns false for an out-of-range index.
    bool markFeatureEdge(int edge);
    // The same by vertex pair, for a caller that has the polyline and not the
    // edge indices. Returns false when the two vertices share no edge.
    bool markFeatureEdge(int va, int vb);

    // Pin a node regardless of where it sits. The escape hatch for a caller
    // with a constraint this class knows nothing about -- a symmetry plane, a
    // node the physics fixes, the centre of a templated disk.
    void pinVertex(int v);
    // Release every pin and reclassify. Does not forget marked feature edges.
    void unpinAll();

    // ---- targets ---------------------------------------------------------

    // W = h * I for every element: the ideal square of side h. Passing h <= 0
    // uses the mean edge length of the current mesh, which is the usual
    // "keep the size you have, fix the shape" target.
    void setUniformTargets(double h = 0.0);
    // W per element from its own current shape, so the metric measures
    // departure from where the mesh started rather than from a square. The
    // Jacobian at the element centre is used, which is the mean of the four
    // corner Jacobians for a bilinear map.
    void setTargetsFromCurrentShape();
    // W = h_q * I with a per-element size, for a target coming from the layout
    // (a chord's interval length, say) rather than from one global number.
    void setSizeTargets(const std::vector<double> &hPerQuad);

    inline Jacobian2 targetAt(int q, int /*corner*/) const {
        if (q < 0 || q >= static_cast<int>(targetJacobian.size()))
            return {1.0, 0.0, 0.0, 1.0};
        return targetJacobian[q];
    }

    // ---- geometry --------------------------------------------------------

    // Jacobian of the bilinear map of quad `q` at reference (xi, eta) in
    // [0,1]^2. Column 0 is dx/dxi, column 1 is dx/deta.
    Jacobian2 jacobianAt(int q, double xi, double eta) const;
    // The same at corner `c`, where it reduces to the two element sides
    // leaving that corner -- so its determinant is the corner's cross product
    // and its sign is the usual convexity test.
    Jacobian2 cornerJacobian(int q, int c) const;

    // Derivative of the bilinear shape function of corner `c` at (xi, eta),
    // returned as {dN/dxi, dN/deta}. This is the whole of what a TMOP gradient
    // needs beyond the Jacobians: dA/dx_c is the outer product of that pair
    // with the coordinate direction.
    static Point shapeGrad(int c, double xi, double eta);
    static double shapeFn(int c, double xi, double eta);

    double signedArea(int q) const;
    Point centroid(int q) const;
    // Scaled Jacobian at corner `c`: det / (|column 0| |column 1|), i.e. the
    // sine of the corner angle. Negative means the corner has turned over.
    double scaledJacobian(int q, int c) const;
    // The worst of the four. This is the number the corpus reports as
    // "min scaled Jacobian" and the one an untangler drives above zero.
    double minScaledJacobian(int q) const;
    // Mean of the two mid-side spans in each direction, as a pair -- the
    // element's size in xi and in eta. What a size target is compared against.
    void elementSpans(int q, double &sXi, double &sEta) const;

    // Which quad contains `p`, or -1. Brute force over the elements, splitting
    // each into two triangles; the counterpart of Mesh::findTriangleContainingPoint
    // and no faster.
    int findQuadContainingPoint(const Point &p) const;

    void computeQuality();

    // BlockDecomposition::alignmentQuality() read on the elements themselves:
    // each element's two mid-side directions held to `frame`'s cross at its
    // centre, by area and the frame's |u|, graded by FeatureFrame::Alignment
    // (docs/block_decomposition_metrics.md Sec. 7). 1 for a mesh whose
    // element sides run along or across the walls' cross wherever it has a
    // direction, 1/2 for one at a random angle to it. Not in `quality`, which
    // needs nothing but the mesh: this needs the model's frame.
    double alignmentQuality(const FeatureFrame &frame) const;

    // ---- conversion and I/O ----------------------------------------------

    // Split each quad into two triangles on its shorter diagonal and return
    // the result as a Mesh, so anything already written against Mesh -- point
    // location, the CSR, circle fitting -- can be run on a quad mesh. The
    // vertex indices are preserved, so a per-vertex field maps across
    // unchanged.
    Mesh toTriangleMesh() const;

    bool writeOBJ(const std::string &filename) const;
    // Legacy-VTK unstructured grid with VTK_QUAD cells, carrying the material
    // id per cell and the node type per point so a smoother's degrees of
    // freedom can be seen in ParaView.
    bool writeVTU(const std::string &filename) const;

    // ---- MFEM ------------------------------------------------------------

    // How the .mesh writer turns this mesh's materials and boundary into the
    // attributes MFEM reads.
    struct MFEMOptions {
        // Boundary attribute for every boundary segment when `boundaryByAxis`
        // is off. MFEM requires a positive attribute; 1 is the usual "one wall,
        // one tag" answer.
        int boundaryAttribute = 1;

        // Tag boundary segments by which coordinate the wall holds fixed:
        // **1** for a segment lying in a vertical wall (x fixed), **2** for one
        // in a horizontal wall (y fixed), `boundaryAttribute` for anything
        // else. This is MFEM's own 1/2/3 = fixed-x/y/z convention and what a
        // solver imposing v.n = 0 per wall expects. Off by default because it
        // is meaningful only on an axis-aligned domain -- a curved boundary
        // gets one tag and the convention says nothing useful about it.
        bool boundaryByAxis = false;
        // How far off axis-aligned a segment may be, as |cos| or |sin| of its
        // direction, and still count as a wall for the rule above.
        double axisTolerance = 1e-9;

        bool writeBoundary = true;

        // Digits written per coordinate. 17 round-trips a double exactly, and
        // that matters here: a mesh that is watertight in memory and not on
        // disk is a mesh a solver will find cracks in.
        int precision = 17;

        // Emit the `# material <id>: <n> elements` comment block. Cheap, and
        // it is the thing to read first when a hydro run comes back
        // single-material.
        bool annotate = true;
    };

    struct MFEMReport {
        int elements = 0;
        int boundaryElements = 0;
        int vertices = 0;          // written, after unused ones are dropped
        int unusedVertices = 0;    // present in this mesh, not in the file
        int attributes = 0;        // distinct element attributes written
        // (attribute, element count), ascending by attribute. The per-material
        // zone census a solver prints at startup, available before the file is
        // even handed over.
        std::vector<std::pair<int, int>> attributeCounts;
        // True when an attribute had to be shifted to make it positive; see
        // the note in writeMFEM.
        bool attributesRemapped = false;
    };

    // Write MFEM's `MFEM mesh v1.0` format: a 2-D unstructured grid of SQUARE
    // elements, each carrying **its material id as its element attribute**.
    //
    // That is the whole point of this writer, and the reason it is worth
    // having one rather than exporting through VTK. MFEM addresses material
    // regions by element attribute, not by position -- Laghos' problem 8, for
    // one, reads `T.Attribute` inside a Coefficient to give the inclusions a
    // different density from the matrix -- so a mesh whose attributes are all
    // 1 silently runs single-material no matter how carefully the interface
    // was meshed. Elements generated after the layout, like the O-grid inside
    // a templated inclusion, are exactly the ones at risk of losing the tag,
    // which is why `quadMatId` is carried through the conversion above rather
    // than re-derived from geometry here.
    //
    // MFEM requires strictly positive attributes. Ids are written unchanged
    // whenever every one of them is already >= 1 -- the common case, since the
    // pipeline's ids come from the .geo's Physical Surface tags and start at 1
    // -- and otherwise every id is shifted by one constant so the smallest
    // becomes 1, which preserves the partition and the gaps. `report`
    // (optional) says whether that happened, along with the zone census.
    //
    // Only the domain boundary becomes boundary elements. Material interfaces
    // deliberately do not: in MFEM an interior interface is the boundary
    // between two attributes and is found from them, and emitting interior
    // segments as boundary elements would make the mesh non-conforming to
    // MFEM's own reading of it.
    //
    // Vertices no element references are dropped and the rest renumbered, so
    // the file is always compact even when this mesh is not.
    bool writeMFEM(const std::string &filename) const;
    bool writeMFEM(const std::string &filename, const MFEMOptions &mopts,
                   MFEMReport *report = nullptr) const;

    // (material id, number of quads) ascending by id. What the writer reports,
    // available without writing anything.
    std::vector<std::pair<int, int>> materialCounts() const;

private:
    void buildEdges();
    void buildBoundary();
    void buildVertexQuads();
    void buildBoundaryLoops();
    void markInterfaceEdges();
    // Fix one node on every closed feature loop that has no fixed node on it
    // already. Run at the end of classifyNodes(); see the .cxx for why it is
    // unconditional rather than an Options flag.
    void anchorFreeFeatureLoops();
    // The feature graph cut at every node that cannot slide: each run is
    // [fixed, sliding..., fixed], and carries at least one sliding node.
    std::vector<std::vector<int>> featureChains() const;

    // Caller-pinned vertices, kept across a rebuild so classifyNodes() can
    // re-apply them.
    std::vector<bool> pinned;
    // Caller-marked feature edges, by vertex pair so they survive a rebuild
    // that renumbers edges.
    std::vector<Edge> markedFeatures;
};

}  // namespace mesh

#endif  // __MESH_QUADMESH_HXX__
