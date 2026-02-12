#ifndef __SEPARATRIX_TRACE_HXX__
#define __SEPARATRIX_TRACE_HXX__

#include "crossfield/CrossField.hxx"
#include <deque>

// Data structures assumed from previous context
struct TracePoint {
    int face_id;
    std::array<double, 3> barycentric;
    Point global_pos;
};

struct Separatrix {
    std::deque<TracePoint> path;
    int id;
    int origin_singularity_id;
    bool active = true;
};

class SeparatrixTrace {
public:
    std::vector<Separatrix> separatrices;
    std::shared_ptr<CrossField> crossField;

    SeparatrixTrace(std::shared_ptr<CrossField> cf);

    // Initialize separatrices from all singularities in the cross field
    void initializeSeparatrices();

private:
    // Compute the starting point and directions for a single singularity
    // Returns pairs of (TracePoint, direction_angle) for each port
    std::vector<std::pair<TracePoint, double>> computeSingularityPorts(int triangleIndex, double crossFieldIndex);
    
    // Get the cross field direction at a vertex, rotated into the triangle's plane
    // Returns the angle of the nearest cross component relative to the reference axis
    double getFieldAlignmentAngle(int triangleIndex, int vertexLocalIndex);
    
    // Compute the intersection of a ray from a point inside a triangle with one of its edges
    // Returns: (intersection point, edge index 0/1/2, parametric t along edge, neighbor triangle index)
    // Edge i connects vertex i to vertex (i+1)%3
    std::tuple<Point, int, double, int> rayEdgeIntersection(
        int triangleIndex, const Point& origin, double direction);
    
    // Convert global coordinates to barycentric coordinates for a given triangle
    std::array<double, 3> globalToBarycentric(int triangleIndex, const Point& p);
};

#endif // __SEPARATRIX_TRACE_HXX__
