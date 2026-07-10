#ifndef __UV_ISO_TRACE_HXX__
#define __UV_ISO_TRACE_HXX__

#include "Parameterization/UVGParam.hxx"
#include "IntervalTree.hxx"
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_set>

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

    void traceIsolines(); // Main method to trace all u and v isolines and populate integralCurves_

    const std::vector<IntegralCurve>& getIntegralCurves() const { return integralCurves_; }

    // Write the cut mesh + all traced isocurves to a VTK legacy unstructured grid
    // file for viewing in ParaView.  Returns false on I/O error.
    //
    // The file contains:
    //   - Mesh triangles as VTK_TRIANGLE (type 5) cells.
    //   - Each isocurve as a VTK_POLY_LINE (type 4) cell.
    //
    // Cell data arrays:
    //   cell_type  (int)    : 0 = mesh triangle, 1 = u-isoline, 2 = v-isoline
    //   coord_value (double): constant u or v value of the isoline (0 for triangles)
    //
    // Point data arrays (on mesh vertices and curve sample points):
    //   u_coord, v_coord (double): UV parameterization values at each point
    bool writeVTK(const std::string& filename) const;

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


    std::vector<IntegralCurve> integralCurves_; // List of all traced integral curves (u and v isolines)
};

#endif // __UV_ISO_TRACE_HXX__