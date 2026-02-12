#include "tracing/SeparatrixTrace.hxx"
#include <cmath>
#include <tuple>

SeparatrixTrace::SeparatrixTrace(std::shared_ptr<CrossField> cf) : crossField(cf) {
    separatrices.clear();
}

std::array<double, 3> SeparatrixTrace::globalToBarycentric(int triangleIndex, const Point& p) {
    const Mesh& mesh = crossField->getMesh();
    const Triangle& tri = mesh.triangles[triangleIndex];
    
    const Point& v0 = mesh.vertices[tri[0]];
    const Point& v1 = mesh.vertices[tri[1]];
    const Point& v2 = mesh.vertices[tri[2]];
    
    // Compute barycentric coordinates using the standard formula
    double denom = (v1[1] - v2[1]) * (v0[0] - v2[0]) + (v2[0] - v1[0]) * (v0[1] - v2[1]);
    
    double lambda0 = ((v1[1] - v2[1]) * (p[0] - v2[0]) + (v2[0] - v1[0]) * (p[1] - v2[1])) / denom;
    double lambda1 = ((v2[1] - v0[1]) * (p[0] - v2[0]) + (v0[0] - v2[0]) * (p[1] - v2[1])) / denom;
    double lambda2 = 1.0 - lambda0 - lambda1;
    
    return {lambda0, lambda1, lambda2};
}

std::tuple<Point, int, double, int> SeparatrixTrace::rayEdgeIntersection(
    int triangleIndex, const Point& origin, double direction) {
    
    const Mesh& mesh = crossField->getMesh();
    const Triangle& tri = mesh.triangles[triangleIndex];
    
    // Direction vector from angle
    Point dir = {std::cos(direction), std::sin(direction)};
    
    // Get triangle vertices
    std::array<Point, 3> verts = {
        mesh.vertices[tri[0]],
        mesh.vertices[tri[1]],
        mesh.vertices[tri[2]]
    };
    
    // Check intersection with each edge
    // Edge i connects vertex i to vertex (i+1)%3
    double bestT = std::numeric_limits<double>::max();
    int bestEdge = -1;
    double bestEdgeParam = 0.0;
    Point bestIntersection = origin;
    
    for (int i = 0; i < 3; ++i) {
        int j = (i + 1) % 3;
        
        const Point& e0 = verts[i];
        const Point& e1 = verts[j];
        
        // Edge direction
        Point edgeDir = {e1[0] - e0[0], e1[1] - e0[1]};
        
        // Solve: origin + t * dir = e0 + s * edgeDir
        // This gives us:
        // origin.x + t * dir.x = e0.x + s * edgeDir.x
        // origin.y + t * dir.y = e0.y + s * edgeDir.y
        
        double denom = dir[0] * edgeDir[1] - dir[1] * edgeDir[0];
        
        if (std::abs(denom) < 1e-12) {
            // Ray is parallel to edge
            continue;
        }
        
        // Vector from edge start to ray origin
        double dx = origin[0] - e0[0];
        double dy = origin[1] - e0[1];
        
        // Parametric value along the ray
        double t = (edgeDir[0] * dy - edgeDir[1] * dx) / denom;
        
        // Parametric value along the edge
        double s = (dir[0] * dy - dir[1] * dx) / denom;
        
        // Check if intersection is valid:
        // - t > epsilon (ray goes forward, not backward)
        // - s in [0, 1] (intersection is on the edge)
        const double eps = 1e-10;
        if (t > eps && s >= -eps && s <= 1.0 + eps) {
            // Clamp s to [0, 1]
            s = std::max(0.0, std::min(1.0, s));
            
            if (t < bestT) {
                bestT = t;
                bestEdge = i;
                bestEdgeParam = s;
                bestIntersection = {
                    e0[0] + s * edgeDir[0],
                    e0[1] + s * edgeDir[1]
                };
            }
        }
    }
    
    // Get the neighbor triangle across this edge
    int neighborTriangle = -1;
    if (bestEdge >= 0) {
        neighborTriangle = mesh.triangleAdjacency[triangleIndex][bestEdge];
    }
    
    return {bestIntersection, bestEdge, bestEdgeParam, neighborTriangle};
}

double SeparatrixTrace::getFieldAlignmentAngle(int triangleIndex, int vertexLocalIndex) {
    const Mesh& mesh = crossField->getMesh();
    const Triangle& tri = mesh.triangles[triangleIndex];
    
    // Get global vertex index for the chosen node
    int vertexIndex = tri[vertexLocalIndex];
    
    // Get the complex field value at this vertex (u^4 representation)
    std::complex<double> u4 = crossField->u_k[vertexIndex];
    
    // Extract the base angle: u = e^(i*theta) means u^4 = e^(i*4*theta)
    // So theta = arg(u^4) / 4
    double fieldAngle = std::arg(u4) / 4.0;
    
    // Compute the reference axis: ray from barycenter to vertex n
    Point v0 = mesh.vertices[tri[0]];
    Point v1 = mesh.vertices[tri[1]];
    Point v2 = mesh.vertices[tri[2]];
    
    // Barycenter
    Point barycenter = {
        (v0[0] + v1[0] + v2[0]) / 3.0,
        (v0[1] + v1[1] + v2[1]) / 3.0
    };
    
    // Vertex position
    Point vertexPos = mesh.vertices[vertexIndex];
    
    // Reference axis direction (from barycenter to vertex)
    Point refAxis = {
        vertexPos[0] - barycenter[0],
        vertexPos[1] - barycenter[1]
    };
    
    // Angle of reference axis
    double refAngle = std::atan2(refAxis[1], refAxis[0]);
    
    // alpha is the angle the field makes relative to the reference axis
    // We want the smallest such angle (modulo pi/2 symmetry of cross field)
    double alpha = fieldAngle - refAngle;
    
    // Normalize to [-pi/4, pi/4) to get the nearest cross component
    while (alpha > M_PI_4) alpha -= M_PI_2;
    while (alpha <= -M_PI_4) alpha += M_PI_2;
    
    return alpha;
}

std::vector<std::pair<TracePoint, double>> SeparatrixTrace::computeSingularityPorts(int triangleIndex, double crossFieldIndex) {
    std::vector<std::pair<TracePoint, double>> ports;
    
    const Mesh& mesh = crossField->getMesh();
    const Triangle& tri = mesh.triangles[triangleIndex];
    
    // Step 1: Compute the starting point (barycenter of the triangle)
    Point v0 = mesh.vertices[tri[0]];
    Point v1 = mesh.vertices[tri[1]];
    Point v2 = mesh.vertices[tri[2]];
    
    Point barycenter = {
        (v0[0] + v1[0] + v2[0]) / 3.0,
        (v0[1] + v1[1] + v2[1]) / 3.0
    };
    
    // Barycentric coordinates of barycenter are (1/3, 1/3, 1/3)
    std::array<double, 3> baryCoords = {1.0/3.0, 1.0/3.0, 1.0/3.0};
    
    // Step 2: Compute alpha (field alignment angle)
    // Use vertex 0 as the reference node
    double alpha = getFieldAlignmentAngle(triangleIndex, 0);
    
    // Step 3: Compute d from cross field index
    // Cross field index = d/4, so d = 4 * index
    // For index +1/4 (3-valent): d = 1
    // For index -1/4 (5-valent): d = -1
    int d = static_cast<int>(std::round(4.0 * crossFieldIndex));
    
    // Number of separatrices = |4 - d|
    // For d=1: 3 separatrices
    // For d=-1: 5 separatrices
    int numPorts = std::abs(4 - d);
    
    // Step 4: Compute port angles using the formula:
    // theta_port = alpha + (2*pi*k) / (4 - d)
    // But we need to express these as absolute angles (not relative to reference)
    
    // Get reference axis angle (from barycenter to vertex 0)
    Point refAxis = {
        v0[0] - barycenter[0],
        v0[1] - barycenter[1]
    };
    double refAngle = std::atan2(refAxis[1], refAxis[0]);
    
    for (int k = 0; k < numPorts; ++k) {
        double theta_port = alpha + (2.0 * M_PI * k) / static_cast<double>(4 - d);
        
        // Convert to absolute angle
        double absoluteAngle = refAngle + theta_port;
        
        // Normalize to [-pi, pi)
        while (absoluteAngle > M_PI) absoluteAngle -= 2.0 * M_PI;
        while (absoluteAngle <= -M_PI) absoluteAngle += 2.0 * M_PI;
        
        // Create the trace point
        TracePoint tp;
        tp.face_id = triangleIndex;
        tp.barycentric = baryCoords;
        tp.global_pos = barycenter;
        
        ports.emplace_back(tp, absoluteAngle);
    }
    
    return ports;
}

void SeparatrixTrace::initializeSeparatrices() {
    separatrices.clear();
    
    const auto& singularities = crossField->singularTriangles;
    
    int separatrixId = 0;
    
    for (size_t singIdx = 0; singIdx < singularities.size(); ++singIdx) {
        int triangleIndex = singularities[singIdx].first;
        double crossFieldIndex = singularities[singIdx].second;
        
        // Compute the starting ports for this singularity
        std::vector<std::pair<TracePoint, double>> ports = computeSingularityPorts(triangleIndex, crossFieldIndex);
        
        // Create a separatrix for each port
        for (const auto& [startPoint, direction] : ports) {
            Separatrix sep;
            sep.id = separatrixId++;
            sep.origin_singularity_id = static_cast<int>(singIdx);
            sep.active = true;
            
            // First point: barycenter of singular triangle
            sep.path.push_back(startPoint);
            
            // Second point: intersection with triangle edge
            auto [intersection, edgeIdx, edgeParam, neighborTri] = 
                rayEdgeIntersection(triangleIndex, startPoint.global_pos, direction);
            
            if (edgeIdx >= 0) {
                // Valid intersection found
                TracePoint edgePoint;
                edgePoint.global_pos = intersection;
                
                // The edge point lies on the boundary between two triangles
                // We store it as belonging to the neighbor triangle (where tracing continues)
                // If neighborTri == -1, we've hit the domain boundary
                if (neighborTri >= 0) {
                    edgePoint.face_id = neighborTri;
                    edgePoint.barycentric = globalToBarycentric(neighborTri, intersection);
                } else {
                    // Hit boundary - store in current triangle, mark separatrix as inactive
                    edgePoint.face_id = triangleIndex;
                    edgePoint.barycentric = globalToBarycentric(triangleIndex, intersection);
                    sep.active = false;
                }
                
                sep.path.push_back(edgePoint);
            } else {
                // No valid intersection (shouldn't happen for a point inside a triangle)
                sep.active = false;
            }
            
            separatrices.push_back(sep);
        }
    }
}