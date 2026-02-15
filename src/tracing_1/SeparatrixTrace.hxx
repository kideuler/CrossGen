#ifndef __SEPARATRIX_TRACE_HXX__
#define __SEPARATRIX_TRACE_HXX__

#include "crossfield/CrossField.hxx"
#include <deque>
#include <set>
#include <unordered_map>

// Data structures assumed from previous context
struct TracePoint {
    int face_id;
    std::array<double, 3> barycentric;
    Point global_pos;
};

enum class TerminationReason {
    RUNNING,
    HIT_BOUNDARY,
    SELF_INTERSECTION,
    SINGULARITY_CUTOFF,
    MAX_STEPS_REACHED,
    MERGED_TANGENTIAL,
    LIMIT_CYCLE,
    REPEATED_CROSSING  // Hit the same other separatrix twice (phase 2 retracing)
};

struct Separatrix {
    std::deque<TracePoint> path;
    int id;
    int origin_singularity_id;
    bool active = true;
    TerminationReason termination_reason = TerminationReason::RUNNING;
    int merged_with_id = -1; // ID of separatrix this was merged with (-1 if not merged)

    std::set<int> visited_edges;
};

class SeparatrixTrace {
public:
    std::vector<Separatrix> separatrices;
    std::shared_ptr<CrossField> crossField;

    SeparatrixTrace(std::shared_ptr<CrossField> cf);

    // Initialize separatrices from all singularities in the cross field
    void initializeSeparatrices();

    void trace();

private:
    // Helper: Interpolate the field angle at a barycentric coordinate.
    // Matches the result to 'prevAngle' to ensure continuity.
    double interpolateFieldDirection(int triangleIndex, const std::array<double, 3>& bary, double prevAngle);

    // Compute the starting point and directions for a single singularity
    std::vector<std::pair<TracePoint, double>> computeSingularityPorts(int triangleIndex, double crossFieldIndex);
    
    // Get the cross field direction at a vertex, rotated into the triangle's plane
    double getFieldAlignmentAngle(int triangleIndex, int vertexLocalIndex);
    
    // Compute the intersection of a ray from a point inside a triangle with one of its edges
    std::tuple<Point, int, double, int> rayEdgeIntersection(
        int triangleIndex, const Point& origin, double direction);
    
    // Convert global coordinates to barycentric coordinates for a given triangle
    std::array<double, 3> globalToBarycentric(int triangleIndex, const Point& p);

    // Cache for quick singularity lookups during tracing
    std::unordered_map<int, double> singularityMap;
};

#endif // __SEPARATRIX_TRACE_HXX__
