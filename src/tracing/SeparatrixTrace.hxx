#ifndef __SEPARATRIX_TRACE_HXX__
#define __SEPARATRIX_TRACE_HXX__

#include <memory>
#include <unordered_map>
#include <vector>

#include "tracing/FieldTracer.hxx"

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
// ---------------------------------------------------------------------------

enum class TerminationReason {
    RUNNING,
    EXIT_BOUNDARY,      // condition 1
    CROSSED_TWICE,      // condition 2
    CUT_AT_SINGULARITY, // condition 3
    HETEROCLINIC,       // met another separatrix head-on and was joined to it
    LIMIT_CYCLE,        // ran past the step budget without ever meeting anything
    STUCK               // the walk could not be continued; a bug or a broken mesh
};

enum class SeparatrixOrigin { Singularity, BoundaryCorner };

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

        // Divide a corner's wedge into `quarters` equal parts rather than into
        // exact right angles measured off the boundary. On a corner that really
        // is 270 degrees the two agree; on one that measures 226 they do not,
        // and right angles leave a 46-degree scrap of wedge whose component is
        // a sliver running along the boundary.
        bool evenCornerRays = true;
    };

    // `useActualSingularityCoordinates` is passed through to FieldTracer: solve
    // for the point where the representation vector vanishes rather than taking
    // the barycentre of the singular triangle.
    SeparatrixTrace(std::shared_ptr<CrossField> cf, bool useActualSingularityCoordinates,
                    const Settings &settings);
    explicit SeparatrixTrace(std::shared_ptr<CrossField> cf,
                             bool useActualSingularityCoordinates = true)
        : SeparatrixTrace(std::move(cf), useActualSingularityCoordinates, Settings()) {}

    // Advance every live separatrix by one triangle and test what it met.
    // Tracing them together rather than one after another is what makes the
    // "crossed the same separatrix twice" test independent of the order the
    // separatrices happen to be numbered in.
    void stepAndCheck();

    // stepAndCheck() until nothing is live.
    void run();

    std::vector<Separatrix> separatrices;
    std::vector<Singularity> singularities;    // mirror of the tracer's, with ports filled in
    std::vector<BoundaryCorner> boundaryCorners;
    std::vector<SeparatrixCrossing> crossings;

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

    // Register the segments a step just added and act on what they crossed.
    // Returns false when the separatrix was terminated part way through them.
    bool registerNewSegments(Separatrix &sep, int firstNewIndex);

    void terminateAt(Separatrix &sep, int segIndex, const Point &at, TerminationReason why,
                     int onSep, int onSeg);

    Settings settings;
    std::unique_ptr<FieldTracer> tracer;
    const Mesh *mesh = nullptr;

    std::vector<Walker> walkers;                     // one per separatrix
    std::vector<int> stepsTaken;
    std::vector<double> cutRadius;                   // per singularity, in model units
    std::vector<std::vector<SegRef>> segmentsOfTriangle;
    std::vector<std::unordered_map<int, int>> crossCount;  // per separatrix: other id -> times met
};

#endif // __SEPARATRIX_TRACE_HXX__
