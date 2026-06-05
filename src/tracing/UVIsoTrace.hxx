#ifndef __UV_ISO_TRACE_HXX__
#define __UV_ISO_TRACE_HXX__

#include "Parameterization/UVGParam.hxx"
#include "IntervalTree.hxx"
#include <iostream>
#include <deque>

const int MAX_TRACE_STEPS = 5000; // Maximum number of steps to trace an isoline before termination

enum class IsoTraceType {
    U_ISOLINE,
    V_ISOLINE
};

enum class TraceTerminationReason {
    RUNNING, // still running
    EXIT_BOUNDARY, // Hit the mesh boundary
    HIT_INVALID_TRIANGLE, // Entered an invalid triangle (e.g., flipped or near singularity)
    MAX_STEPS_REACHED, // Reached maximum step count without termination
};

struct IntegralCurve {
    std::deque<Point> points;          // 2D spatial coordinates of the integral curve
    std::deque<Point> uv_points;       // Corresponding (u,v) coordinates along the curve
    std::deque<int> traversedFaces;   // Sequential indices of triangles crossed

    IsoTraceType traceType;            // Type of isoline (u or v)
    double coordinateValue;            // The constant parameter value tracked (e.g., u = C or v = C)
    TraceTerminationReason terminationReason = TraceTerminationReason::RUNNING; // Reason for termination
};


class UVIsoTrace {
public:
    UVIsoTrace(std::shared_ptr<UVGParam> uvParam, int nU, int nV);

    void printQueries() const; // For debugging: print the interval trees

private:
    std::shared_ptr<UVGParam> uvParam_;
    int nU_;
    int nV_;
    double deltaU_;
    double deltaV_;
    Point uvMin_;
    Point uvMax_;
    std::vector<bool> isValidTriangle_; // Flags for triangles that are valid for tracing (not near singularities or flipped)
    IntervalTree uIntervalTree_; // Interval tree for fast triangle lookup along u-isolines
    IntervalTree vIntervalTree_; // Interval tree for fast triangle lookup along v-isolines

    std::unordered_set<int> TrianglesInterstedByUIsolines; // Set of triangle indices intersected by any u-isoline (for quick validity checks to make sure we got all branches)
};

#endif // __UV_ISO_TRACE_HXX__