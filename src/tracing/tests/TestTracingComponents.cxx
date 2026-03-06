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
// Test registry and main
//=============================================================================
using TestFunc = std::function<bool()>;

int main(int argc, char** argv) {
    // Register all tests
    std::map<int, std::pair<std::string, TestFunc>> tests = {
        {1, {"ConstantFieldTrace", test1_ConstantFieldTrace}},
        {2, {"RigidBodyRotation", test2_RigidBodyRotation}},
        {3, {"ParabolicStreamline", test3_ParabolicStreamline}},
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