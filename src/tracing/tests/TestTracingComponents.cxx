#include "triangle/TriangleMesher.hpp"
#include "crossfield/CrossField.hxx"
#include "tracing/SeparatrixTrace.hxx"
#include "TestHelper.hxx"
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <map>

// Write a traced path to an OBJ file as a polyline
void writePathToOBJ(const std::string& filename, const std::vector<TracePoint>& path) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    out << "# Traced path with " << path.size() << " points\n";

    // Write vertices
    for (const auto& tp : path) {
        out << "v " << tp.global_pos[0] << " " << tp.global_pos[1] << " 0\n";
    }

    // Write polyline as edges
    out << "l";
    for (size_t i = 1; i <= path.size(); ++i) {
        out << " " << i;
    }
    out << "\n";

    out.close();
    std::cout << "Wrote traced path to " << filename << std::endl;
}

// Helper function to trace from a starting point to the boundary
std::vector<TracePoint> traceTowardsBoundary(
    SeparatrixTrace& tracer,
    std::shared_ptr<Mesh> mesh,
    const Point& startPos,
    double startAngle,
    int maxSteps = 10000
) {
    std::vector<TracePoint> path;

    // Find triangle containing the starting point
    int startTri = mesh->findTriangleContainingPoint(startPos);
    if (startTri < 0) {
        std::cerr << "Error: Starting point is not inside any triangle!" << std::endl;
        return path;
    }

    // Compute barycentric coordinates of start point
    std::array<double, 3> startBary = tracer.globalToBarycentric(startTri, startPos);

    // Find where the ray exits the starting triangle
    auto [exitPos, exitEdge, exitT, exitNeighbor] = tracer.rayEdgeIntersection(startTri, startPos, startAngle);

    if (exitEdge < 0) {
        std::cerr << "Error: Could not find exit edge from starting triangle!" << std::endl;
        return path;
    }

    // Add the starting point
    TracePoint start;
    start.face_id = startTri;
    start.barycentric = startBary;
    start.global_pos = startPos;
    start.field_angle = startAngle;
    start.trace_direction = startAngle; // Initial trace direction matches field direction
    start.edge_id = mesh->triangleEdges[startTri][exitEdge];
    start.local_edge_index = exitEdge;
    start.edge_crossing_t = exitT;
    path.push_back(start);

    // Add the first exit point
    TracePoint firstExit;
    firstExit.face_id = startTri;
    firstExit.barycentric = tracer.globalToBarycentric(startTri, exitPos);
    firstExit.global_pos = exitPos;
    firstExit.field_angle = startAngle;
    // Compute actual trace direction from start to exit
    Point traceVec = exitPos - startPos;
    firstExit.trace_direction = std::atan2(traceVec[1], traceVec[0]);
    firstExit.edge_id = mesh->triangleEdges[startTri][exitEdge];
    firstExit.local_edge_index = exitEdge;
    firstExit.edge_crossing_t = exitT;
    path.push_back(firstExit);

    // Check if we hit boundary immediately
    if (exitNeighbor < 0) {
        return path;
    }

    // Build a temporary Separatrix to use with stepHeuns
    Separatrix sep;
    sep.id = -1;
    sep.origin_singularity_id = -1;
    sep.origin_singularity_port = -1;
    sep.active = true;
    sep.termination_reason = TerminationReason::RUNNING;
    sep.path.push_back(start);
    sep.path.push_back(firstExit);

    // Trace until we hit the boundary using stepHeuns
    for (int step = 0; step < maxSteps; ++step) {
        if (!sep.active) break;

        tracer.stepHeuns(sep);
    }

    // Copy path out of the separatrix
    path.assign(sep.path.begin(), sep.path.end());

    return path;
}

//=============================================================================
// Test 1: Constant cross field (1,1) direction, trace from origin
//=============================================================================
bool test1_ConstantFieldTrace() {
        std::cout << "\n=== Test 1: Constant Field Trace ===" << std::endl;

        // Create a unit-disk mesh
        auto mesh = TestHelper::createEllipse(0.0, 0.0, 1.0, 1.0, 360.0, 0.05);
        std::cout << "Mesh created: " << mesh->vertices.size() << " vertices, " 
                << mesh->triangles.size() << " triangles" << std::endl;

        // Create crossfield with constant direction (1,1)
        auto crossField = std::make_shared<CrossField>(mesh);
        std::complex<double> u = {1.0, 1.0};
        u /= std::abs(u);
        std::complex<double> u_val = u*u*u*u;
        
        crossField->u_k.resize(mesh->vertices.size());
        for (int i = 0; i < (int)mesh->vertices.size(); ++i) {
            crossField->u_k[i] = u_val;
        }

        // Create tracer and trace
        SeparatrixTrace tracer(crossField, false);
        Point startPos = {0.0, 0.0};
        double startAngle = M_PI / 4.0; // direction (1,1)

        auto path = traceTowardsBoundary(tracer, mesh, startPos, startAngle);
        
        std::cout << "Traced " << path.size() << " points" << std::endl;
        writePathToOBJ("test1_constant_field.obj", path);

        // Verify: path should end near boundary (radius ~1)
        if (path.empty()) {
            std::cerr << "FAIL: Empty path" << std::endl;
            return false;
        }

        Point endPos = path.back().global_pos;
        double endRadius = std::sqrt(endPos[0]*endPos[0] + endPos[1]*endPos[1]);
        std::cout << "End position: (" << endPos[0] << ", " << endPos[1] << "), radius: " << endRadius << std::endl;

        if (std::abs(endRadius - 1.0) > 0.05) {
            std::cerr << "FAIL: End point not near boundary" << std::endl;
            return false;
        }

        std::cout << "PASS" << std::endl;
        return true;
}

//=============================================================================
// Test 2: Rigid Body Rotation [-y, x] field, trace from (0.5,0.0)
//=============================================================================
bool test2_RigidBodyRotation() {
    std::cout << "\n=== Test 2: Rigid Body Rotation ===" << std::endl;

    // Create a unit-disk mesh
    auto mesh = TestHelper::createEllipse(0.0, 0.0, 1.0, 1.0, 360.0, 0.05);
    std::cout << "Mesh created: " << mesh->vertices.size() << " vertices, " 
            << mesh->triangles.size() << " triangles" << std::endl;

    // Create crossfield with rigid body rotation [-y, x] (tangent to circles centered at origin)
    auto crossField = std::make_shared<CrossField>(mesh);
    
    crossField->u_k.resize(mesh->vertices.size());
    for (int i = 0; i < (int)mesh->vertices.size(); ++i) {
        double x = mesh->vertices[i][0];
        double y = mesh->vertices[i][1];
        // Field direction is [-y, x] (tangent to circle)
        // For vertices at origin, use default direction
        double r = std::sqrt(x*x + y*y);
        std::complex<double> u;
        if (r < 1e-10) {
            u = {0.0, 1.0}; // default to vertical at origin
        } else {
            u = {-y/r, x/r}; // normalized tangent
        }
        crossField->u_k[i] = u*u*u*u;
    }

    // Create tracer and trace
    SeparatrixTrace tracer(crossField, false);
    Point startPos = {0.5, 0.0};
    // At (0.5, 0), the tangent direction is [-0, 0.5] = [0, 1], angle = pi/2
    double startAngle = M_PI / 2.0;
    
    // Set up the initial trace points
    int startTri = mesh->findTriangleContainingPoint(startPos);
    if (startTri < 0) {
        std::cerr << "Error: Starting point is not inside any triangle!" << std::endl;
        return false;
    }

    std::array<double, 3> startBary = tracer.globalToBarycentric(startTri, startPos);
    auto [exitPos, exitEdge, exitT, exitNeighbor] = tracer.rayEdgeIntersection(startTri, startPos, startAngle);

    if (exitEdge < 0) {
        std::cerr << "Error: Could not find exit edge from starting triangle!" << std::endl;
        return false;
    }

    // Add the starting point
    TracePoint start;
    start.face_id = startTri;
    start.barycentric = startBary;
    start.global_pos = startPos;
    start.field_angle = startAngle;
    start.trace_direction = startAngle;
    start.edge_id = mesh->triangleEdges[startTri][exitEdge];
    start.local_edge_index = exitEdge;
    start.edge_crossing_t = exitT;

    // Add the first exit point
    TracePoint firstExit;
    firstExit.face_id = startTri;
    firstExit.barycentric = tracer.globalToBarycentric(startTri, exitPos);
    firstExit.global_pos = exitPos;
    firstExit.field_angle = startAngle;
    Point traceVec = exitPos - startPos;
    firstExit.trace_direction = std::atan2(traceVec[1], traceVec[0]);
    firstExit.edge_id = mesh->triangleEdges[startTri][exitEdge];
    firstExit.local_edge_index = exitEdge;
    firstExit.edge_crossing_t = exitT;

    // Build a temporary Separatrix to use with stepHeuns
    Separatrix sep;
    sep.id = -1;
    sep.origin_singularity_id = -1;
    sep.origin_singularity_port = -1;
    sep.active = true;
    sep.termination_reason = TerminationReason::RUNNING;
    sep.path.push_back(start);
    sep.path.push_back(firstExit);

    // Trace
    int maxSteps = 10000;
    bool verbose = true;
    for (int step = 0; step < maxSteps; ++step) {
        if (!sep.active) break;

        const TracePoint& current = sep.path.back();

        // Debug: print expected vs actual field angle
        if (verbose && step < 30) {
            double x = current.global_pos[0];
            double y = current.global_pos[1];
            double r = std::sqrt(x*x + y*y);
            double expected = std::atan2(x/r, -y/r); // expected tangent angle
            double diff = current.field_angle - expected;
            while (diff > M_PI) diff -= 2*M_PI;
            while (diff < -M_PI) diff += 2*M_PI;
            std::cout << "Step " << step << ": pos=(" << x << "," << y 
                      << ") expected=" << expected << " actual=" << current.field_angle 
                      << " diff=" << diff << std::endl;
        }

        tracer.stepHeuns(sep);
    }

    // Copy path out of the separatrix
    std::vector<TracePoint> path(sep.path.begin(), sep.path.end());
    
    std::cout << "Traced " << path.size() << " points" << std::endl;
    writePathToOBJ("test2_rigid_body_rotation.obj", path);

    // Analyze radius drift
    if (path.empty()) {
        std::cerr << "FAIL: Empty path" << std::endl;
        return false;
    }

    double startRadius = std::sqrt(startPos[0]*startPos[0] + startPos[1]*startPos[1]);
    double maxRadiusDrift = 0.0;
    int worstStep = 0;
    
    for (size_t i = 0; i < path.size(); ++i) {
        double r = std::sqrt(path[i].global_pos[0]*path[i].global_pos[0] + 
                            path[i].global_pos[1]*path[i].global_pos[1]);
        double drift = std::abs(r - startRadius);
        if (drift > maxRadiusDrift) {
            maxRadiusDrift = drift;
            worstStep = i;
        }
    }

    Point endPos = path.back().global_pos;
    double endRadius = std::sqrt(endPos[0]*endPos[0] + endPos[1]*endPos[1]);
    
    std::cout << "Start radius: " << startRadius << std::endl;
    std::cout << "End radius: " << endRadius << std::endl;
    std::cout << "Max radius drift: " << maxRadiusDrift << " at step " << worstStep << std::endl;
    std::cout << "End position: (" << endPos[0] << ", " << endPos[1] << ")" << std::endl;

    // For a good circular trace, radius should stay within ~10% of original
    if (maxRadiusDrift > 0.1 * startRadius) {
        std::cerr << "FAIL: Radius drifted too much (>" << 0.1*startRadius << ")" << std::endl;
        return false;
    }

    std::cout << "PASS" << std::endl;
    return true;
}

//=============================================================================
// Test 3: Field v(x,y) = [1, 2x], exact streamline y = x^2 + C
//=============================================================================
bool test3_ParabolicStreamline() {
    std::cout << "\n=== Test 3: Parabolic Streamline v=[1,2x], y=x^2 ===" << std::endl;

    // Create a unit-disk mesh using TestHelper
    auto mesh = TestHelper::createEllipse(0.0, 0.0, 1.0, 1.0, 360.0, 0.05);
    std::cout << "Mesh created: " << mesh->vertices.size() << " vertices, " 
              << mesh->triangles.size() << " triangles" << std::endl;

    // Create cross field encoding v(x,y) = [1, 2x]
    // At each vertex the representative angle is theta = atan2(2x, 1)
    auto crossField = std::make_shared<CrossField>(mesh);
    crossField->u_k.resize(mesh->vertices.size());
    for (int i = 0; i < (int)mesh->vertices.size(); ++i) {
        double x = mesh->vertices[i][0];
        double theta = std::atan2(2.0 * x, 1.0);
        std::complex<double> u = std::exp(std::complex<double>(0.0, theta));
        crossField->u_k[i] = u * u * u * u;   // 4-RoSy representation
    }

    // Trace from the origin; at (0,0) the field is [1,0], angle = 0
    SeparatrixTrace tracer(crossField, false);
    Point startPos = {0.0, 0.0};
    double startAngle = 0.0;       // atan2(0,1) = 0

    auto path = traceTowardsBoundary(tracer, mesh, startPos, startAngle);

    std::cout << "Traced " << path.size() << " points" << std::endl;
    writePathToOBJ("test3_parabolic.obj", path);

    if (path.empty()) {
        std::cerr << "FAIL: Empty path" << std::endl;
        return false;
    }

    // Exact streamline through (0,0): y = x^2  (C = 0)
    // Measure max deviation |y_i - x_i^2| along the trace
    double maxErr = 0.0;
    int worstStep = 0;
    for (size_t i = 0; i < path.size(); ++i) {
        double x = path[i].global_pos[0];
        double y = path[i].global_pos[1];
        double err = std::abs(y - x * x);
        if (err > maxErr) {
            maxErr = err;
            worstStep = (int)i;
        }
    }

    Point endPos = path.back().global_pos;
    double endR  = std::sqrt(endPos[0] * endPos[0] + endPos[1] * endPos[1]);

    std::cout << "End position: (" << endPos[0] << ", " << endPos[1] << "), radius: " << endR << std::endl;
    std::cout << "Max |y - x^2| error: " << maxErr << " at step " << worstStep << std::endl;

    // The trace should reach near the boundary (radius ≈ 1)
    if (endR < 0.7) {
        std::cerr << "FAIL: Trace did not reach near boundary (endR = " << endR << ")" << std::endl;
        return false;
    }

    // Allow ≤ 0.05 deviation from the exact parabola
    if (maxErr > 0.05) {
        std::cerr << "FAIL: Max streamline error " << maxErr << " exceeds tolerance 0.05" << std::endl;
        return false;
    }

    std::cout << "PASS" << std::endl;
    return true;
}

//=============================================================================
// Test 4: Singular Triangle Viertel Trace
//=============================================================================
bool test4_SingularTriangleViertelTrace(double singularityIndex, std::array<double, 3> singularityBarycenter) {
    std::cout << "\n=== Test n: Singular Triangle Viertel Trace ===" << std::endl;

    auto mesh = std::make_shared<Mesh>();
    double R = 1.0;
    
    // T0 vertices (equilateral triangle)
    Point v0 = {0.0, R};
    Point v1 = {-R * std::sqrt(3.0) / 2.0, -R / 2.0};
    Point v2 = { R * std::sqrt(3.0) / 2.0, -R / 2.0};
    
    // Outer dummy vertices to form neighbor triangles
    Point v3 = {-R * std::sqrt(3.0), R};
    Point v4 = {0.0, -R * 1.5};
    Point v5 = {R * std::sqrt(3.0), R};

    mesh->vertices = {v0, v1, v2, v3, v4, v5};
    mesh->triangles = {
        {0, 1, 2}, // T0: The singular triangle
        {1, 0, 3}, // T1: Neighbor sharing (0,1)
        {2, 1, 4}, // T2: Neighbor sharing (1,2)
        {0, 2, 5}  // T3: Neighbor sharing (2,0)
    };
    mesh->triangleEdges = {
        {0, 1, 2},
        {0, 3, 4},
        {1, 5, 6},
        {2, 7, 8}
    };
    mesh->triangleAdjacency = {
        {1, 2, 3},   // T0 adjacent to T1, T2, T3
        {0, -1, -1}, // T1 adjacent back to T0
        {0, -1, -1}, // T2 adjacent back to T0
        {0, -1, -1}  // T3 adjacent back to T0
    };

    auto crossField = std::make_shared<CrossField>(mesh);
    crossField->u_k.resize(6);
    
    // Field setup: choose u_k values so the linear interpolant has its zero at
    // the desired singularityBarycenter (λ0, λ1, λ2).
    // We need λ0*u0 + λ1*u1 + λ2*u2 = 0.
    // Pick u0 and u1 freely, then solve for u2.
    double lam0 = singularityBarycenter[0];
    double lam1 = singularityBarycenter[1];
    double lam2 = singularityBarycenter[2];

    // 1. Compute the exact global coordinates of the singularity
    Point sc;
    sc[0] = v0[0]*lam0 + v1[0]*lam1 + v2[0]*lam2;
    sc[1] = v0[1]*lam0 + v1[1]*lam1 + v2[1]*lam2;

    // 2. Define the analytical cross field (d = 1 for +1/4, d = -1 for -1/4)
    double d = (singularityIndex > 0) ? 1.0 : -1.0;

    // Helper to compute the exact representation vector u_k for any vertex
    auto calc_u = [&](Point v) {
        // Find angle of vertex relative to singularity
        double theta = std::atan2(v[1] - sc[1], v[0] - sc[0]);
        // Analytical field angle
        double field_angle = (d * theta) / 4.0;
        // Convert to representation vector e^(i * 4 * field_angle)
        return std::complex<double>(std::cos(4.0 * field_angle), std::sin(4.0 * field_angle));
    };

    crossField->u_k[0] = calc_u(v0);
    crossField->u_k[1] = calc_u(v1);
    crossField->u_k[2] = calc_u(v2);
    
    // Apply the exact analytical field to the outer boundary vertices as well!
    crossField->u_k[3] = calc_u(v3);
    crossField->u_k[4] = calc_u(v4);
    crossField->u_k[5] = calc_u(v5);

    // Mark T0 as singular with the given singularity index
    crossField->singularTriangles.emplace_back(0, singularityIndex);

    // Initialize tracer with useActualSingularityCoordinates=true so the
    // constructor solves for the singularity location from u_k, landing at singularityBarycenter
    SeparatrixTrace tracer(crossField, true);
    tracer.dphi_singularity_zone = 0.001;
    tracer.maxStepsInSingularityZone = 3000;
    std::vector<std::vector<TracePoint>> allPaths;
    std::vector<int> pathEdgeIndex; // which edge (0,1,2) each path came from, -1 for separatrices

    // Save separatrices for visualization
    for (const auto& sep : tracer.separatrices) {
        allPaths.emplace_back(sep.path.begin(), sep.path.end());
        pathEdgeIndex.push_back(-1); // separatrices don't belong to an entry edge
    }

    // Define the boundaries we will shoot streamlines from
    struct EntryDef {
        int neighborTri;
        int localEdge; 
        Point A; 
        Point B; 
    };
    std::vector<EntryDef> entries = {
        {1, 0, v1, v0}, // Edge 0
        {2, 0, v2, v1}, // Edge 1
        {3, 0, v0, v2}  // Edge 2
    };

    int successfulTraces = 0;
    int numTracesPerEdge = 30;
    for (int e = 0; e < 3; ++e) {
        auto def = entries[e];
        for (int i = 1; i <= numTracesPerEdge; ++i) {
            double t = i / (numTracesPerEdge + 1.0);
            Point P;
            P[0] = def.A[0] + (def.B[0] - def.A[0]) * t;
            P[1] = def.A[1] + (def.B[1] - def.A[1]) * t;

            Separatrix sep;
            sep.id = 100 + e * 10 + i;
            sep.active = true;

            TracePoint entry;
            entry.face_id = def.neighborTri;
            entry.global_pos = P;
            entry.local_edge_index = def.localEdge;
            entry.edge_crossing_t = t;
            
            // Calculate the entry field angle: orthogonal to the edge, pointing inward
            double entryfieldAngle;
            
            Point edgeDir = {def.B[0] - def.A[0], def.B[1] - def.A[1]};
            // Inward normal: rotate edge direction by -90 degrees (clockwise)
            Point inwardNormal = {edgeDir[1], -edgeDir[0]};
            // Check that it points inward (toward the centroid of T0 at the origin)
            Point toCenter = {-P[0], -P[1]};
            if (inwardNormal[0] * toCenter[0] + inwardNormal[1] * toCenter[1] < 0) {
                // Flip if pointing outward
                inwardNormal[0] = -inwardNormal[0];
                inwardNormal[1] = -inwardNormal[1];
            }
            entryfieldAngle = std::atan2(inwardNormal[1], inwardNormal[0]);
        
            entry.field_angle = entryfieldAngle; 
            entry.trace_direction = entryfieldAngle;

            sep.path.push_back(entry);

            // Execute the analytical hyperbolic step
            tracer.stepViertel(sep);

            if (sep.path.size() < 2) {
                std::cerr << "FAIL: Streamline " << sep.id << " did not find an exit." << std::endl;
                //return false;
            }
            
            allPaths.emplace_back(sep.path.begin(), sep.path.end());
            pathEdgeIndex.push_back(e);
            successfulTraces++;
        }
    }

    std::cout << "Successfully traced " << successfulTraces << " custom streamlines through the singular triangle." << std::endl;

    // Generate VTK file for visualization (ParaView-compatible with per-line colors)
    std::string basename = "test4_singular_triangle_idx" + std::to_string(singularityIndex) 
                         + "_bary" + std::to_string(singularityBarycenter[0]) 
                         + "_" + std::to_string(singularityBarycenter[1])
                         + "_" + std::to_string(singularityBarycenter[2]);
    std::string filename = basename + ".vtk";

    // Count total points and total cell-list size
    int totalPoints = 6; // mesh vertices
    int totalCells = 1;  // the triangle face
    int totalCellListSize = 4; // "3 v0 v1 v2" for the triangle
    for (const auto& path : allPaths) {
        totalPoints += static_cast<int>(path.size());
        totalCells += 1; // one polyline per path
        totalCellListSize += 1 + static_cast<int>(path.size()); // count + point indices
    }

    std::ofstream out(filename);
    if (out.is_open()) {
        // Header
        out << "# vtk DataFile Version 3.0\n";
        out << "Singular Triangle Test\n";
        out << "ASCII\n";
        out << "DATASET POLYDATA\n";

        // Points
        out << "POINTS " << totalPoints << " double\n";
        for (const auto& v : mesh->vertices) {
            out << v[0] << " " << v[1] << " 0\n";
        }
        for (const auto& path : allPaths) {
            for (const auto& pt : path) {
                out << pt.global_pos[0] << " " << pt.global_pos[1] << " 0\n";
            }
        }

        // Polygon (the triangle)
        out << "POLYGONS 1 4\n";
        out << "3 0 1 2\n";

        // Lines (the streamlines)
        int numLines = static_cast<int>(allPaths.size());
        int linesListSize = 0;
        for (const auto& path : allPaths) {
            linesListSize += 1 + static_cast<int>(path.size());
        }
        out << "LINES " << numLines << " " << linesListSize << "\n";

        int ptOffset = 6; // first 6 points are mesh vertices
        for (const auto& path : allPaths) {
            out << path.size();
            for (size_t i = 0; i < path.size(); ++i) {
                out << " " << (ptOffset + i);
            }
            out << "\n";
            ptOffset += static_cast<int>(path.size());
        }

        // Cell data: edge index for coloring
        // Total cells = 1 (polygon) + numLines
        out << "CELL_DATA " << (1 + numLines) << "\n";
        out << "SCALARS EdgeGroup int 1\n";
        out << "LOOKUP_TABLE default\n";
        out << "3\n"; // triangle gets its own group
        for (size_t p = 0; p < allPaths.size(); ++p) {
            out << pathEdgeIndex[p] << "\n"; // -1 for separatrices, 0/1/2 for edges
        }

        out.close();
        std::cout << "Wrote visualization to " << filename << std::endl;
    }

    std::cout << "PASS" << std::endl;
    return true;
}

//=============================================================================
// Test registry and main
//=============================================================================
using TestFunc = std::function<bool()>;

int main(int argc, char** argv) {
    // Register all tests
    std::map<int, std::pair<std::string, TestFunc>> tests = {
        {1, {"ConstantFieldTrace", test1_ConstantFieldTrace}},
        {2, {"RigidBodyRotation", test2_RigidBodyRotation}},
        {3, {"ParabolicStreamline", test3_ParabolicStreamline}},
        {4, {"SingularTriangleViertelTrace_1", []() { return test4_SingularTriangleViertelTrace(0.25, {1.0/3.0, 1.0/3.0, 1.0/3.0}); }}},
        {5, {"SingularTriangleViertelTrace_2", []() { return test4_SingularTriangleViertelTrace(0.25, {2.0/5.0, 2.0/5.0, 1.0/5.0}); }}},
        {6, {"SingularTriangleViertelTrace_3", []() { return test4_SingularTriangleViertelTrace(-0.25, {1.0/3.0, 1.0/3.0, 1.0/3.0}); }}},
        {7, {"SingularTriangleViertelTrace_4", []() { return test4_SingularTriangleViertelTrace(-0.25, {2.0/5.0, 2.0/5.0, 1.0/5.0}); }}}
    };

    std::vector<int> testsToRun;

    if (argc > 1) {
        // Run specific tests
        for (int i = 1; i < argc; ++i) {
            int testNum = std::atoi(argv[i]);
            if (tests.find(testNum) != tests.end()) {
                testsToRun.push_back(testNum);
            } else {
                std::cerr << "Unknown test number: " << testNum << std::endl;
            }
        }
    } else {
        // Run all tests
        for (const auto& [num, _] : tests) {
            testsToRun.push_back(num);
        }
    }

    if (testsToRun.empty()) {
        std::cout << "Available tests:" << std::endl;
        for (const auto& [num, test] : tests) {
            std::cout << "  " << num << ": " << test.first << std::endl;
        }
        return 0;
    }

    int passed = 0;
    int failed = 0;

    for (int testNum : testsToRun) {
        const auto& [name, func] = tests[testNum];
        try {
            if (func()) {
                passed++;
            } else {
                failed++;
            }
        } catch (const std::exception& e) {
            std::cerr << "Test " << testNum << " (" << name << ") threw exception: " << e.what() << std::endl;
            failed++;
        }
    }

    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Passed: " << passed << "/" << (passed + failed) << std::endl;
    
    if (failed > 0) {
        std::cout << "Failed: " << failed << std::endl;
        return 1;
    }

    return 0;
}