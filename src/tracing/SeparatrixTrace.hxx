#ifndef __SEPARATRIX_TRACE_HXX__
#define __SEPARATRIX_TRACE_HXX__

#include <memory>
#include <unordered_map>
#include <vector>

#include "tracing/FieldTracer.hxx"

class Interfaces;

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
    TANGENTIAL          // met another separatrix running the same way, at a shallow angle
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
    };

    // `useActualSingularityCoordinates` is passed through to FieldTracer: solve
    // for the point where the representation vector vanishes rather than taking
    // the barycentre of the singular triangle.
    // `interfaces` is the material interface network of the model, or null.
    // It is read, not kept beyond the object's life, and must outlive it: the
    // layout built on this trace reads the branches off it too.
    SeparatrixTrace(std::shared_ptr<CrossField> cf, bool useActualSingularityCoordinates,
                    const Settings &settings, const Interfaces *interfaces = nullptr);
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
    bool registerNewSegments(Separatrix &sep, int firstNewIndex);

    // A separatrix that has just reached the boundary: if it arrived tangentially
    // and a corner of the model is within reach, move its last point onto that
    // corner, so no node is planted where the arrival happened to land.
    bool snapTangentialLanding(Separatrix &sep);

    void terminateAt(Separatrix &sep, int segIndex, const Point &at, TerminationReason why,
                     int onSep, int onSeg);

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
};

#endif // __SEPARATRIX_TRACE_HXX__
