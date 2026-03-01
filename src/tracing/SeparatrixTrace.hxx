#ifndef __SEPARATRIX_TRACE_HXX__
#define __SEPARATRIX_TRACE_HXX__

#include "crossfield/CrossField.hxx"
#include <deque>
#include <queue>
#include <set>
#include <unordered_map>
#include <utility>
#include <iostream>

// Robustness tolerances.
static const double EPS_INTERSECT_T = 1e-10;      // param epsilon for segment intersections
static const double EPS_BARY = 1e-9;              // barycentric inside tolerance
static const double MIN_ADVANCE_REL = 1e-7;       // minimum accepted segment length as fraction of avg edge
static const int MAX_STEPS = 50000;


// TracePoint represents a point along a separatrix trace, storing the triangle it is in, its barycentric coordinates, and its global position.
struct TracePoint {
    int face_id;
    std::array<double, 3> barycentric;
    Point global_pos;
    double field_angle; // interpolated cross field angle at this point from the current triangle
    double trace_direction; // actual traced direction (angle from entry to exit), used for phase continuity
    int edge_id = -1; // ID of the edge crossed to be used in the next iteration
    int local_edge_index = -1; // local edge index (0,1,2) of the triangle for the edge crossed to be used in the next iteration
    double edge_crossing_t = 0.0; // parametric t along the edge where the crossing occurred, used for interpolation
};

// TerminationReason enumerates the possible reasons for a separatrix trace to terminate, such as hitting a boundary, self-intersection, or reaching a maximum step count.
enum class TerminationReason {
    RUNNING, // still running
    EXIT_BOUNDARY, // Hit the mesh boundary
    CONNECT_TANGENTIAL_PRIMARY, // Connected tangentially to another separatrix and is being kept
    CONNECT_TANGENTIAL_SECONDARY, // Connected tangentially to another separatrix and is being removed
    LIMIT_CYCLE,
    ORTHOGONAL_TO_SINGULARITY_SEPARATRIX,
    MAX_STEPS_REACHED
};

struct Separatrix {
    std::deque<TracePoint> path; // The sequence of points along the separatrix trace.
    int id; // Unique identifier for the separatrix.
    int origin_singularity_id; // ID of the singularity from which this separatrix originates.
    int origin_singularity_port; // Port index at the singularity (0-4) corresponding to the initial direction.
    bool active = true; // Indicates whether the trace is still active (not terminated).
    bool in_singularity_zone = false; // Whether the trace is currently within the "singularity zone" of any singularity.
    TerminationReason termination_reason = TerminationReason::RUNNING; // Reason for termination if not active.

    std::set<int> visited_edges; // Set of edge IDs visited by this separatrix, used for self-intersection detection.
};

class SeparatrixTrace {
public:
    std::vector<Separatrix> separatrices; // List of all separatrices being traced
    std::shared_ptr<CrossField> crossField; // Shared pointer to the cross field 
    
    SeparatrixTrace(std::shared_ptr<CrossField> cf);

    // Convert global coordinates to barycentric coordinates for a given triangle
    std::array<double, 3> globalToBarycentric(int triangleIndex, const Point& p);

    // Compute the intersection of a ray from a point inside a triangle with one of its edges
    // Returns: (intersection point, local edge index, edge parameter t, neighbor triangle index)
    // excludeEdge: optional edge index to exclude (e.g., the entry edge when tracing)
    std::tuple<Point, int, double, int> rayEdgeIntersection(int triangleIndex, const Point& origin, double direction, int excludeEdge = -1);

    // sister method which takes a direction vector instead of an angle
    std::tuple<Point, int, double, int> rayEdgeIntersection(int triangleIndex, const Point& origin, const Point& direction, int excludeEdge = -1);

    // For a given reference angle and another angle, compute the angle that is closest to the reference angle but still matches the cross field direction at that point
    double makeAngleSamePhase(double referenceAngle, double angleCandidate);

    // Find phase difference between two angles and return k such that angleCandidate + k*(pi/2) is closest to referenceAngle
    int findPhaseDifference(double referenceAngle, double angleCandidate);

    // Perform one tracing step using Heun's method, returning the next TracePoint
    TracePoint stepHeuns(const TracePoint& current);

private:
    // Private members for internal use during tracing
    std::vector<bool> isSingularTriangle; // Precomputed lookup for whether a triangle is singular
    std::unordered_map<int, std::pair<std::vector<int>, bool>> triangleSeparatrixMap; // triangle index -> (list of separatrix IDs passing through, is an intersection present)
    std::queue<int> Intersections; // Queue of triangle indices where intersections have been detected, to be processed.
};


#endif // __SEPARATRIX_TRACE_HXX__