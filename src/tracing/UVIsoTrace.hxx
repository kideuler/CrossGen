#ifndef __UV_ISO_TRACE_HXX__
#define __UV_ISO_TRACE_HXX__

#include "Parameterization/UVGParam.hxx"

enum class IsoTraceType {
    U_ISOLINE,
    V_ISOLINE
};

struct IntegralCurve {
    std::vector<Point> points;          // 2D spatial coordinates of the integral curve
    std::vector<Point> uv_points;       // Corresponding (u,v) coordinates along the curve
    std::vector<int> traversedFaces;   // Sequential indices of triangles crossed

    IsoTraceType traceType;            // Type of isoline (u or v)
    double coordinateValue;            // The constant parameter value tracked (e.g., u = C or v = C)
};


class UVIsoTrace {
public:
    UVIsoTrace(std::shared_ptr<UVGParam> uvParam, int nU, int nV);

private:
    std::shared_ptr<UVGParam> uvParam_;
    int nU_;
    int nV_;
    double deltaU_;
    double deltaV_;
    Point uvMin_;
    Point uvMax_;
    std::vector<bool> isValidTriangle_; // Flags for triangles that are valid for tracing (not near singularities or flipped)
};

#endif // __UV_ISO_TRACE_HXX__