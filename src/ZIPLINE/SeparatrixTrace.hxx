#ifndef __SEPARATRIX_TRACE_HXX__
#define __SEPARATRIX_TRACE_HXX__

#include <memory>
#include <unordered_map>
#include <vector>

#include "ZIPLINE/FieldTracer.hxx"

class Interfaces;
class ConformalField;

// ---------------------------------------------------------------------------
// The separatrices of a cross field, traced until the stopping conditions of
// Viertel, Osting and Staten, IMR 2019, Sec. 3.3 apply. What comes out is the
// input to QuadLayout, which turns the arrangement of these curves into a quad
// layout with T-junctions.
//
// Streamlines are launched from two kinds of place:
//
//   * every interior singularity, one along each of its 4 - d ports;
//
//   * every corner of the boundary, one along each axis direction of the field
//     that points into the domain. A corner whose interior angle is m pi/2 has
//     m - 1 of them: none at a convex 90-degree corner, whose two boundary
//     edges are already sides of the layout, two at a 270-degree reflex corner,
//     three at a slit. That count comes from the corner index of Table 1 and
//     not from testing directions against the geometry, for the reason set out
//     in MotorcycleGraph::launch(): a corner that measures 275 degrees rather
//     than 270 admits a spurious third ray which then runs alongside the
//     boundary for the length of the model.
//
// By the Poincare-Bendixson theorem a streamline on a bounded surface either
// leaves through the boundary, joins singularities, or winds onto a limit
// cycle; on a discrete field the middle case never happens exactly, so the
// tracing has to be cut off by hand. Three conditions do it (Sec. 3.3):
//
//   1. it leaves through the boundary;
//   2. it crosses the same separatrix a second time -- which is what a limit
//      cycle does immediately, and what an honest streamline does not;
//   3. it crosses a separatrix of a singularity orthogonally inside that
//      singularity's own triangle, which is the near-miss that would otherwise
//      leave two streamlines running side by side down the whole model.
//
// The last two end a separatrix in the interior, i.e. at a T-junction. Removing
// those is what Sec. 4 is for and is deliberately not done here.
//
// Two meetings the conditions do not cover are settled the way
// docs/viertel_2019.md Sec. 4.4 settles them: two separatrices meeting head-on
// at a shallow angle are one curve and are joined (HETEROCLINIC), and two
// meeting at a shallow angle running the same way are one curve too, so the
// arriving one stops on the other (TANGENTIAL). A join cuts the other curve
// back, and whatever had stopped on the part that went is resumed.
// ---------------------------------------------------------------------------

enum class TerminationReason {
    RUNNING,
    EXIT_BOUNDARY,      // condition 1
    CROSSED_TWICE,      // condition 2
    CUT_AT_SINGULARITY, // condition 3
    HETEROCLINIC,       // met another separatrix head-on and was joined to it
    LIMIT_CYCLE,        // ran past the step budget without ever meeting anything
    STUCK,              // the walk could not be continued; a bug or a broken mesh
    TANGENTIAL,         // met another separatrix running the same way, at a shallow angle
    // Settings::captureAtSingularities (docs/quad_layout_pipeline.md Secs. 7-9):
    AT_SINGULARITY,     // entered a singular point's critical patch and was joined to it
    REDUNDANT,          // the other copy of a separatrix traced from both ends: removed
    CUT_LIMIT_CYCLE     // a limit cycle, cut at its first orthogonal crossing: a T-junction
};

enum class SeparatrixOrigin { Singularity, BoundaryCorner, InterfaceNode };

struct Separatrix {
    std::vector<TracePoint> path;   // path[0] is the origin; segment i runs path[i-1] -> path[i]
    int id = -1;
    SeparatrixOrigin originKind = SeparatrixOrigin::Singularity;
    int origin_singularity_id = -1;  // index into singularities, or into boundaryCorners
    int origin_singularity_port = -1;
    bool active = true;
    TerminationReason termination_reason = TerminationReason::RUNNING;

    // Where it stopped, when it stopped on another separatrix (conditions 2
    // and 3). `endOnSegment` indexes into that separatrix's path the same way.
    int endOnSeparatrix = -1;
    int endOnSegment = -1;

    // Where it stopped, when it left through the boundary.
    int endBoundaryEdge = -1;    // global mesh edge id, -1 if it left at a vertex
    int endBoundaryVertex = -1;

    // AT_SINGULARITY: the singular point it was joined to -- an interior
    // singularity, a corner of dS or an interface node, by `endKind`, indexed
    // as origin_singularity_id is -- the port of that point it came in along
    // (Sec. 7.6's arrival slot), and how far off that port's direction it
    // arrived, in radians.
    SeparatrixOrigin endKind = SeparatrixOrigin::Singularity;
    int endSingularity = -1;
    int endPort = -1;
    double arrivalMisfit = 0.0;
};

// Two separatrices meeting transversally, recorded as they are traced so the
// layout does not have to rediscover them.
struct SeparatrixCrossing {
    int sepA = -1, segA = -1;
    int sepB = -1, segB = -1;
    Point pos{0.0, 0.0};
};

// A separatrix crossing a material interface. The field is aligned to the
// interface from both sides, so a streamline meets it square and carries on
// into the other material; the crossing is a node of the layout, on the
// interface arc as well as on the separatrix. It always sits at a hand-off
// between two triangles, so it is a point of the path, `pathIndex`, and lies
// on mesh edge `edge` or -- where the walk went through a vertex -- on
// `vertex`.
struct InterfaceCrossing {
    int sep = -1;
    int pathIndex = -1;
    Point pos{0.0, 0.0};
    int edge = -1;
    int vertex = -1;
    double angle = 0.0;   // between the separatrix and the interface, in [0, pi/2]
};

// A node of the material interface network (MERIDIAN's Stage 0b Interfaces):
// a junction, a landing on dS, a kink. A node of the layout whether or not
// anything is launched from it.
struct InterfaceNodeInfo {
    int vertex = -1;
    int node = -1;           // index into Interfaces::nodes()
    bool onBoundary = false;
    int index = 0;           // I(v) from the sectors, in quarter turns
    std::vector<int> quarters;   // per sector of Interfaces::Node, as launched from
    std::vector<int> separatrixIds;
};

// A corner of the boundary: a node of the layout whether or not anything is
// launched from it.
struct BoundaryCorner {
    int vertex = -1;
    int quarters = 2;          // interior angle in units of pi/2, rounded (Table 1)
    double interiorAngle = 0.0;
    std::vector<int> separatrixIds;
};

class SeparatrixTrace {
public:
    // The knobs the stopping conditions leave open. They belong to the
    // constructor because the launching and the per-singularity radii are
    // settled there; setting them on a built object would be too late to have
    // any effect.
    struct Settings {
        // Budget per separatrix, in triangles crossed. Only a streamline
        // winding onto a limit cycle in a region no other separatrix reaches
        // gets near it.
        int maxStepsPerSeparatrix = 20000;

        // Third stopping condition. Turning it off traces streamlines straight
        // past singularities, which is worth seeing when judging whether a
        // T-junction was earned.
        bool cutAtSingularities = true;

        // How far from a singularity that condition still applies, in units of
        // the singular triangle's own mean edge.
        //
        // The paper states it inside the singular triangle, which is where the
        // sweep of Sec. 3.2.2 can watch the crossing happen. But whether a near
        // miss lands inside that one triangle or just outside it is an accident
        // of the mesh, and outside it the streamline runs on and ends up
        // alongside the next separatrix of the same singularity -- a sliver of
        // a component with three corners and its neighbour with five. What the
        // condition is really about is a separatrix coming within
        // discretisation error of a singularity, since in the continuum the two
        // curves would have met, and discretisation error is measured in edge
        // lengths.
        //
        // One edge is the smallest radius that leaves every model in
        // data/meshes sound; at zero -- the condition exactly as the paper
        // states it -- two of them keep a separatrix that winds onto a limit
        // cycle no other separatrix ever reaches, and it ends nowhere. Above
        // one the partition gets steadily coarser without the number of
        // components that are not four-sided changing, so there is nothing to
        // be had by going further.
        double singularityCutRadius = 1.0;

        // In the continuum two streamlines of a cross field can only cross each
        // other orthogonally. A crossing at a shallow angle is therefore not a
        // crossing at all but two separatrices following what is really one
        // curve, pushed apart by the discretisation. Anything within this many
        // radians of anti-parallel counts as one of those and is joined up.
        double tangentialAngle = M_PI_4;

        // A separatrix that reaches the boundary this close to running along it
        // has not really arrived there.
        //
        // The field is boundary aligned, so in the continuum a streamline meets
        // the boundary square or runs parallel to it; there is no such thing as
        // a five degree arrival. One happens where the mesh is too coarse to
        // resolve a curved boundary -- Sec. 3.3's "few crosses that are
        // actually aligned with the discrete boundary of the triangle mesh" --
        // and the paper assumes it away ("assuming a sufficiently fine triangle
        // mesh along the boundary such that no separatrices exit tangentially").
        //
        // It cannot be assumed away here, and it is not harmless: the landing
        // splits the boundary at a place that means nothing, giving the
        // component on one side a corner where the boundary in fact runs
        // straight past, and the one on the other side a splinter of no area.
        // That is one pentagon and one sliver triangle from every occurrence,
        // and it is the whole of what stops these layouts being four-sided.
        //
        // Over the sixteen models, 420 of 424 landings arrive between 77 and 90
        // degrees and four arrive between 5 and 20, so there is nothing in
        // between for the threshold to cut through.
        //
        // Off by default, because it does not work. At pi/6 it moves all four
        // and leaves no landing under 30 degrees, and the layout is no better
        // for it: 1116 of 1124 components are four-sided either way, and after
        // Sec. 4 it is 818 of 824 against 820 of 825. The degeneracy does not
        // go, it changes shape. The region between where the streamline landed
        // and the corner it was heading for has no area whichever end the node
        // sits at: planting it at the landing makes that region a splinter
        // triangle, planting it on the corner makes it a two-sided lens, and
        // meanwhile the corner collects every streamline that was running
        // alongside the boundary -- three of them on one 124-degree corner of
        // geom013, splitting its wedge into four 31-degree sectors.
        //
        // What the region needs is contracting, not relabelling, and that is
        // PartitionSimplify's sliver removal -- which skips it only because its
        // short side lies on the boundary and shortening the model is not
        // allowed. Absorbing that side into the neighbouring boundary arcs, the
        // way chord collapse already handles a rung on the boundary, is the fix
        // this was standing in for.
        double boundaryTangentialAngle = 0.0;

        // Where such a landing belongs instead. Each of those four is within an
        // element of a corner of the model, which is what a streamline running
        // alongside the boundary is heading for: at a convex corner the field
        // turns through ninety degrees, so the streamline that came in parallel
        // to one edge leaves along the other, and the corner is where it goes.
        // In mean mesh edges.
        double boundaryCornerSnap = 1.5;

        // A separatrix landing squarely on the boundary within this many mean
        // edges of a corner that emits a separatrix straight back along it
        // ends on that corner, and the corner's own separatrix -- the same
        // streamline traced from the other end -- is cut back to nothing
        // (joinCornerLanding()). Zero turns it off.
        double cornerJoinRadius = 1.0;

        // Divide a corner's wedge into `quarters` equal parts rather than into
        // exact right angles measured off the boundary. On a corner that really
        // is 270 degrees the two agree; on one that measures 226 they do not,
        // and right angles leave a 46-degree scrap of wedge whose component is
        // a sliver running along the boundary.
        bool evenCornerRays = true;

        // Grow every separatrix at the same rate in arc length rather than in
        // triangles (docs/viertel_2019.md Sec. 4.4, "scheduling"). Condition
        // 2 depends on which of two separatrices reaches their second
        // crossing first, and the paper does not say how they are ordered;
        // one triangle each per round orders them by how finely the mesh
        // happens to be cut where each one is, which on a graded mesh is not
        // a property of the field at all. On: the front advances by one mean
        // edge per stepAndCheck(), and within a round the separatrix that is
        // furthest behind always moves next -- the spec's priority queue,
        // cut into rounds so that the viewer can still animate it.
        //
        // Off by default, on measurement. The spec's ordering is its own
        // implementation decision ([I]), not the paper's, and on data/meshes
        // the two orders give the same layout on 45 of the 46 models; the
        // 46th, geom031, is covered 95.6% by blocks growing a triangle at a
        // time and 86.0% growing by arc length (a ring of separatrices round
        // the gear stops on different crossings). The corpus meshes are close
        // to uniform, which is where the two orders agree; on a strongly
        // graded mesh this is the setting to try.
        bool growByArcLength = false;

        // A crossing within `tangentialAngle` of *parallel* is two
        // separatrices running the same way along what in the continuum is
        // one streamline. Counted as an ordinary crossing it leaves a strip
        // between them whose two long sides meet at the crossing at a few
        // degrees -- a component with a cusp for a corner. The spec (Sec.
        // 4.4) stops the arriving one there as a T-junction on the other and
        // flags it, which is what this does; off, such a crossing counts
        // towards condition 2 like any other.
        bool stopParallelTangential = true;

        // A heteroclinic join truncates the separatrix it joins. Anything
        // that had stopped on the part that goes -- a second crossing of it,
        // or a parallel meeting -- is then stopped on nothing, and becomes a
        // loose end of the layout. Such separatrices are resumed from where
        // they stopped (Sec. 4.4, "retroactive truncation"), and the
        // crossings recorded on the removed part are taken back.
        bool resumeAfterTruncation = true;

        // Look for singularities in triangles that touch the boundary too
        // (FieldTracer's `boundaryTriangles`). Off reproduces the older
        // behaviour of reading CrossField::singularTriangles, which skips them.
        bool singularitiesAtBoundary = true;

        // Honour the material interfaces handed to the constructor: their
        // network nodes emit separatrices as boundary corners do, and every
        // separatrix crossing an interface is a node on it. Off traces as if
        // the model had one material, which is the ablation.
        bool respectInterfaces = true;

        // ---- the tracing of docs/quad_layout_pipeline.md Secs. 7-9 ----------
        //
        // Jezdimirovic et al.'s rules instead of Viertel's (P1 Sec. 5, P2
        // Sec. 6.1). Viertel stops a separatrix where it meets what is
        // already there -- a second crossing, a near miss of a singularity
        // -- and leaves a T-junction, because on a field that misses its
        // targets by a fraction of an edge that is what keeps the partition
        // finite. On the corrected field (ConformalField) the targets are
        // where the field says they are, and the near miss can be taken for
        // what it is in the continuum, a connection:
        //
        //   * Sec. 7.5 rule 2: a separatrix entering the *critical patch* of a
        //     singular point -- a disk of captureRadius local mean edges round
        //     an interior singularity, a corner of dS that emits, or an
        //     interface node -- is joined to that point by a straight last
        //     segment (P1 Fig. 13), arriving along the port of it nearest its
        //     arrival (the slot of Sec. 7.6).
        //   * Sec. 7.1: separatrices do *not* stop on one another
        //     (stopAtCrossings off). Every crossing is recorded and nothing
        //     else is done with it.
        //   * Sec. 8 (P1 Algorithm 1 with the slot test of Sec. 8.7c): a
        //     connection between two singular points is traced from both
        //     ends; one copy is kept -- the smaller arrival misfit, Sec. 8.7d
        //     -- and the other removed. A capture whose slot is already
        //     taken, or whose own copy is a connection to somewhere else, is
        //     undone and the separatrix traced on past the point: P1 Sec. 10's
        //     re-propagation (Sec. 12.4 here).
        //   * Sec. 9: what is left crossing one separatrix twice, or reaching
        //     dS far from square, or running out of length, is a limit cycle
        //     and is cut at its first square crossing -- a T-junction, which
        //     Sec. 4 of Viertel and StemExtension then deal with as before.
        //
        // The critical patch trades two failures against each other (Sec.
        // 7.2): too small and a near miss runs on, too large and a separatrix
        // that should pass a point is captured by it.
        //
        // Off here, so SeparatrixTrace on its own traces as Viertel does; on in
        // ZIPLINE::Options, which applies it only when the trace follows the
        // corrected field. Measured over data/meshes/singlemat (h = 0.05,
        // TMOP 200 sweeps), with Viertel's stops kept (stopAtCrossings) and
        // captureRadius at each value, against Viertel alone on the same
        // corrected field:
        //
        //   radius   cover%  full  T-jct  raw T  blocks  folded  TMOP worst  < 0.3
        //   off      99.50    30    14     438    1300     0       0.817       0
        //   0.5      99.67    31    10     144    1280     0       0.817       0
        //   0.75     99.67    31    10      52    1288     1       0.818       0
        //   1.0      99.60    30    12      26    1266     1       0.814       1
        //   1.5      99.72    31     8      24    1171     1       0.820       1
        //
        // ("folded": models with an inverted element before TMOP.) Half an
        // edge is the largest patch that folds nothing: at a full edge the
        // last segment of a capture bends a connection into its point over a
        // whole element, and that hook, a side of two blocks, is where the
        // transfinite grid turns over (geom015). Without Viertel's stops --
        // P2 as written -- the same radius covers 98.8% and leaves one model
        // unsound: what the patch misses winds on and is cut as a limit cycle
        // where nothing can remove the T-junction (geom026, a hole ringed by
        // closed trajectories).
        bool captureAtSingularities = false;
        double captureRadius = 0.5;
        // A capture stands only when the separatrix was running back along a
        // port of the point -- within this fraction of the half-gap between
        // that port and the direction beside it (buildTargets() has the
        // numbers). Further off, it is passing the point and is traced on.
        double captureMisfit = 0.5;
        // Viertel's stopping rules 2 and 3, the joins and the tangential stop.
        // Off is Sec. 7.1: every separatrix runs to dS, to a critical patch,
        // or out of length.
        bool stopAtCrossings = true;
        // Sec. 8's reduction and Sec. 9's limit-cycle cut, with capture on.
        bool reduceSeparatrices = true;
        bool cutLimitCycles = true;
        // Sec. 9.3: within this of square is a crossing a limit cycle may be
        // cut at, and an arrival on dS that is not a limit cycle (the doc's
        // 15-20 degrees).
        double orthogonalTolerance = M_PI / 9.0;
        // With stopAtCrossings off nothing else stops a separatrix winding
        // onto a cycle, so it gets a length as well as a step budget: this
        // many bounding-box diagonals (Sec. 14: "a few times").
        double maxLengthDiagonals = 4.0;

        // A hole no separatrix touches once tracing stops -- an annulus with
        // no corner above all, which has no singularity to trace from (docs/
        // quad_layout_pipeline.md Sec. 6.9) -- is cut by four streamlines
        // square to it (seedUntouchedHoles()). Off leaves the face round it a
        // face that is not a disc, and no blocks.
        bool seedUntouchedHoles = true;
    };

    // What the trace did that is not visible in the separatrices themselves.
    struct Report {
        // Poincare-Hopf, docs/viertel_2019.md Sec. 3.2: the interior indices
        // plus the boundary ones of Table 1 add up to the Euler characteristic
        // of the model, four times over so that all three are integers. It
        // failing means a singularity was missed or a corner misclassified,
        // and the layout cannot close up however the tracing goes.
        int eulerCharacteristic = 0;      // 2 - boundary loops: a planar domain
        int interiorIndexSum4 = 0;        // sum d over singular triangles
        int boundaryIndexSum4 = 0;        // sum (2 - quarters) over corners
        bool poincareHopf = false;
        int multipleSingularities = 0;    // |d| >= 2
        int droppedSingularities = 0;     // no ports: d >= 3 or d <= -5

        int heteroclinicJoins = 0;
        int cornerJoins = 0;              // Settings::cornerJoinRadius
        int resumed = 0;                  // separatrices restarted after a truncation
        int crossingsRetracted = 0;       // crossings that were on a removed tail
        int tangentialParallel = 0;       // TerminationReason::TANGENTIAL
        // Boundary landings arriving within tangentialAngle of the boundary,
        // the spec's `tangential_boundary_exit`: counted, not acted on (see
        // Settings::boundaryTangentialAngle for what acting on them does).
        int tangentialBoundaryExits = 0;

        // Material interfaces (all zero on a single-material model).
        int interfaceNodes = 0;
        int interfaceIndexSum4 = 0;       // sum I(v) over them, part of the check above
        int interfaceEmitted = 0;         // separatrices launched from them
        int absorbedSingularities = 0;    // singular triangles at a node, not traced from
        int interfaceCrossings = 0;
        // Crossings shallower than tangentialAngle: a separatrix running
        // alongside an interface rather than across it.
        int interfaceGrazes = 0;

        // Settings::captureAtSingularities. Captures made (counting each
        // re-made one again), kept as connections, and undone (re-propagated);
        // copies removed, of which `reciprocal` were the true other copy
        // (Sec. 8.7c: start slot of one = arrival slot of the other) and the
        // rest a capture whose other copy had run on elsewhere; limit cycles
        // cut, and ones Sec. 9 found but could not cut (no square crossing).
        int captures = 0;
        int connections = 0;
        int released = 0;
        int redundant = 0;
        int reciprocal = 0;
        int limitCycles = 0;
        int limitCyclesCut = 0;
        // Captures undone because, with a connection already kept between
        // the same two points, they made a loop round nothing (Sec. 8.2).
        int lenses = 0;
        // Sec. 8.6 step 4: singular points with a slot used by no separatrix
        // or by two, after the reduction. Zero is the check.
        int slotAuditFailures = 0;
        // Settings::seedUntouchedHoles: streamlines launched from holes.
        int holeSeeds = 0;
    };

    // `useActualSingularityCoordinates` is passed through to FieldTracer: solve
    // for the point where the representation vector vanishes rather than taking
    // the barycentre of the singular triangle.
    // `interfaces` is the material interface network of the model, or null.
    // It is read, not kept beyond the object's life, and must outlive it: the
    // layout built on this trace reads the branches off it too.
    // `conformal`, when given, is the corrected field of docs/
    // quad_layout_pipeline.md Sec. 6 (ConformalField) and is what is traced:
    // its singularities, its closed-form ports, its exact singular part. `cf`
    // is then that field as a CrossField. With interfaces it is traced too, as
    // long as it was solved for the same network (ConformalField's multi-
    // material constructor: region by region, tangent to the interfaces, the
    // nodes' sector quarters its own, and launched from by those).
    SeparatrixTrace(std::shared_ptr<CrossField> cf, bool useActualSingularityCoordinates,
                    const Settings &settings, const Interfaces *interfaces = nullptr,
                    std::shared_ptr<const ConformalField> conformal = nullptr);
    explicit SeparatrixTrace(std::shared_ptr<CrossField> cf,
                             bool useActualSingularityCoordinates = true)
        : SeparatrixTrace(std::move(cf), useActualSingularityCoordinates, Settings()) {}

    // One round: every live separatrix advanced to the next front -- one mean
    // edge further in arc length (Settings::growByArcLength), or one triangle
    // further -- testing each step against what is already there. Tracing
    // them together rather than one after another is what makes the "crossed
    // the same separatrix twice" test independent of the order the
    // separatrices happen to be numbered in.
    void stepAndCheck();

    // stepAndCheck() until nothing is live.
    void run();

    // How many landings the tangential rule moved onto a corner.
    int tangentialLandings = 0;

    const Report &getReport() const { return report_; }

    std::vector<Separatrix> separatrices;
    std::vector<Singularity> singularities;    // mirror of the tracer's, with ports filled in
    std::vector<BoundaryCorner> boundaryCorners;
    std::vector<SeparatrixCrossing> crossings;
    std::vector<InterfaceNodeInfo> interfaceNodes;
    std::vector<InterfaceCrossing> interfaceCrossings;

    // The interface network being honoured, or null on a single-material
    // model or with Settings::respectInterfaces off.
    const Interfaces *getInterfaces() const { return interfaces; }

    // Boundary loops as vertex rings, oriented with the interior on the left.
    std::vector<std::vector<int>> boundaryLoops;

    std::shared_ptr<CrossField> crossField;
    bool finishedTracing = false;
    int steps = 0;

    const FieldTracer &getTracer() const { return *tracer; }
    const Settings &getSettings() const { return settings; }

    int countActive() const;

private:
    struct SegRef { int sep; int seg; };

    void buildBoundary();
    void launchFromSingularities();
    void launchFromCorners();
    void launchFromInterfaceNodes();
    // Whether the step separatrix k just took handed it across an interface,
    // and if so the crossing.
    void recordInterfaceCrossing(int k);
    // Launch one separatrix from mesh vertex v in direction dir.
    int launchFromVertex(int v, double dir, SeparatrixOrigin kind, int origin, int port);

    // Register the segments a step just added and act on what they crossed.
    // Returns false when the separatrix was terminated part way through them.
    // `stopRules` false records the crossings and applies none of Viertel's
    // stopping rules to them (a capture's last segment).
    bool registerNewSegments(Separatrix &sep, int firstNewIndex, bool stopRules = true);

    // A separatrix that has just reached the boundary: if it arrived tangentially
    // and a corner of the model is within reach, move its last point onto that
    // corner, so no node is planted where the arrival happened to land.
    bool snapTangentialLanding(Separatrix &sep);

    // A separatrix that has just reached the boundary squarely, next to a
    // corner of the model with a separatrix of its own leaving straight back
    // along it: the two are one corner-to-corner (or singularity-to-corner)
    // streamline, traced from both ends. This one is ended on the corner and
    // the corner's is cut back to it, with whatever had stopped on it resumed.
    bool joinCornerLanding(Separatrix &sep);

    void terminateAt(Separatrix &sep, int segIndex, const Point &at, TerminationReason why,
                     int onSep, int onSeg);

    // ---- Settings::captureAtSingularities --------------------------------
    // A singular point a separatrix can be joined to: its position, the
    // radius of its critical patch, what it is, and its ports -- the
    // directions its own separatrices leave in, and which they are.
    struct Target {
        Point pos{0.0, 0.0};
        double radius = 0.0;
        SeparatrixOrigin kind = SeparatrixOrigin::Singularity;
        int index = -1;              // into singularities / boundaryCorners / interfaceNodes
        int triangle = -1;           // the singular triangle; -1 for a mesh vertex
        int vertex = -1;             // the mesh vertex of a corner or a node
        std::vector<double> portDirs;
        std::vector<int> portSeps;
        double maxMisfit = M_PI;     // Settings::captureMisfit, in radians, for this point
        std::vector<double> sideDirs; // the ports and, at a vertex, the walls of its wedge
    };
    void buildTargets();
    // The first place along the segments from `firstNew` on where separatrix
    // k enters a critical patch it may be captured by, and the patch: false
    // when there is none. On success `seg` is the segment and `t` the
    // parameter along it.
    bool findCapture(int k, int firstNew, int &seg, double &t, int &target) const;
    // The straight approach from `from` (in triangle `tri`) to a target, as
    // trace points one per triangle crossed; false if it leaves the mesh.
    bool approach(int tri, const Point &from, const Target &g, std::vector<TracePoint> &out,
                  const Point *via = nullptr) const;
    // Join separatrix k to target g at segment `seg`, parameter t.
    bool capture(int k, int seg, double t, int target, int firstNew);
    // Once nothing is running: keep one copy of each connection, undo the
    // captures that cannot stand, and say whether anything was restarted.
    bool resolveCaptures();
    // Sec. 9, after the reduction.
    void cutLimitCycles();
    // Remove the segments of separatrix k from `fromSeg` on: from the
    // triangle buckets and from the crossings, then cut the path back.
    void dropSegments(int k, int fromSeg);
    // The port of target g nearest direction `dir` (leaving the target), and
    // how far off it is.
    int nearestPort(const Target &g, double dir, double &misfit) const;
    int targetOf(SeparatrixOrigin kind, int index) const;

    // One triangle of separatrix k, and what it ran into. Returns whether it
    // is still live.
    bool advanceOne(int k);

    // Separatrix `host` has just been cut back to end at segment `seg`, point
    // `at`: take back the crossings on the part that went and resume what had
    // stopped against it.
    void retractTail(int host, int seg, const Point &at);

    // Whether point q, known to lie on segment `seg` of separatrix `host` as
    // it was before a truncation, is still on it after it was cut at (seg, at).
    bool survivesTruncation(int host, int seg, const Point &q, int cutSeg, const Point &cut) const;

    Settings settings;
    Report report_;
    const Interfaces *interfaces = nullptr;
    std::vector<char> interfaceNodeVertex;   // per mesh vertex
    std::unique_ptr<FieldTracer> tracer;
    const Mesh *mesh = nullptr;

    std::vector<Walker> walkers;                     // one per separatrix
    std::vector<int> stepsTaken;
    std::vector<double> arcLength;                   // per separatrix, as traced so far
    double front = 0.0;                              // Settings::growByArcLength
    std::vector<int> resumeQueue;                    // restarted by retractTail()
    std::vector<double> cutRadius;                   // per singularity, in model units
    std::vector<std::vector<SegRef>> segmentsOfTriangle;
    std::vector<std::unordered_map<int, int>> crossCount;  // per separatrix: other id -> times met

    // Settings::captureAtSingularities
    std::vector<Target> targets;
    std::vector<std::vector<int>> targetsOfTriangle;   // triangle -> targets whose patch reaches it
    std::vector<int> originTarget;                     // per separatrix
    std::vector<char> leftOrigin;                      // ... has been outside its origin's patch
    std::vector<std::vector<int>> excludedTargets;     // ... released from these
    struct Resume { int tri = -1; Point pos{0.0, 0.0}; double dir = 0.0; int seg = -1; };
    std::vector<Resume> captureResume;                 // ... where a capture can be undone to
    double maxLength = 0.0;
    bool capturesResolved = false;
    int resolvePasses = 0;
    std::vector<char> connectionKept_;                 // ... kept as a connection by a resolution
    // Restart whatever had stopped on separatrix `host` -- a T-junction on it,
    // a join with it -- now that it is gone; true if anything was.
    bool resumeStoppedOn(int host);
    // Settings::seedUntouchedHoles; how many it launched.
    int seedUntouchedHoles();
    bool holesSeeded = false;
    // Whether a closed curve winds round a hole of the model, or round a
    // singularity other than targets skipA and skipB.
    bool enclosesSomething(const std::vector<Point> &loop, int skipA, int skipB) const;
};

#endif // __SEPARATRIX_TRACE_HXX__
