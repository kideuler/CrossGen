// Viewer executable entry point (split out from the old monolithic Viewer.cxx).

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <thread>

#include "viewer/Geometry.hxx"
#include "viewer/Interaction.hxx"
#include "viewer/Render.hxx"

#include "IGM/CutMesh.hxx"
#include "IGM/MIQ.hxx"
#include "polyvector/PolyVectors.hxx"
#include "crossfield/CrossField.hxx"
#include "tracing/SeparatrixTrace.hxx"
#include "triangle/TriangleMesher.hpp"
#include "medialaxis/MedialAxis.hxx"

namespace {

enum class Mode {
    Unselected = 0,
    PolyVector = 1,
    MBO = 2,
    MedialAxis = 3,
};

enum class Phase {
    MeshOnly = 1,
    CrossField = 2,
    Singularities = 3,
    CutSeams = 4,
    UVMesh = 5,
};

enum class MBOPhase {
    MeshOnly = 1,
    CrossField = 2,
    Stepping = 3,
    Separatrices = 4,
    Trace = 5,
};

enum class MedialAxisPhase {
    MeshOnly = 0,
    DelaunayMesh = 1,
    MedialAxis = 2,
};

Phase nextPhase(Phase p) {
    switch (p) {
        case Phase::MeshOnly: return Phase::CrossField;
        case Phase::CrossField: return Phase::Singularities;
        case Phase::Singularities: return Phase::CutSeams;
        case Phase::CutSeams: return Phase::UVMesh;
        case Phase::UVMesh: return Phase::UVMesh;
    }
    return Phase::UVMesh;
}

MBOPhase nextMBOPhase(MBOPhase p) {
    switch (p) {
        case MBOPhase::MeshOnly: return MBOPhase::CrossField;
        case MBOPhase::CrossField: return MBOPhase::Stepping;
        case MBOPhase::Stepping: return MBOPhase::Separatrices;
        case MBOPhase::Separatrices: return MBOPhase::Trace;
        case MBOPhase::Trace: return MBOPhase::Trace;
    }
    return MBOPhase::Trace;
}

MedialAxisPhase nextMedialAxisPhase(MedialAxisPhase p) {
    switch (p) {
        case MedialAxisPhase::MeshOnly: return MedialAxisPhase::DelaunayMesh;
        case MedialAxisPhase::DelaunayMesh: return MedialAxisPhase::MedialAxis;
        case MedialAxisPhase::MedialAxis: return MedialAxisPhase::MedialAxis;
    }
    return MedialAxisPhase::MedialAxis;
}

const char *phaseName(Phase p) {
    switch (p) {
        case Phase::MeshOnly: return "1) mesh";
        case Phase::CrossField: return "2) crossfield";
        case Phase::Singularities: return "3) singularities";
        case Phase::CutSeams: return "4) cut seams";
        case Phase::UVMesh: return "5) UV mesh (MIQ)";
    }
    return "?";
}

const char *mboPhaseName(MBOPhase p) {
    switch (p) {
        case MBOPhase::MeshOnly: return "1) mesh";
        case MBOPhase::CrossField: return "2) MBO crossfield";
        case MBOPhase::Stepping: return "3) MBO stepping";
        case MBOPhase::Separatrices: return "4) separatrices";
        case MBOPhase::Trace: return "5) trace";
    }
    return "?";
}

const char *medialAxisPhaseName(MedialAxisPhase p) {
    switch (p) {
        case MedialAxisPhase::MeshOnly: return "1) mesh";
        case MedialAxisPhase::DelaunayMesh: return "2) Delaunay re-triangulation";
        case MedialAxisPhase::MedialAxis: return "3) Medial axis";
    }
    return "?";
}

const char *modeName(Mode m) {
    switch (m) {
        case Mode::Unselected: return "unselected";
        case Mode::PolyVector: return "PolyVector";
        case Mode::MBO: return "MBO";
        case Mode::MedialAxis: return "Medial Axis";
    }
    return "?";
}

// Format duration in milliseconds with 2 decimal places
std::string formatMs(double ms) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << ms << " ms";
    return oss.str();
}

// Global MBO stepping parameters
const int MBO_MAX_STEPS = 500;

} // namespace

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: Viewer <mesh->obj>\n";
        return 1;
    }

    std::string path = argv[1];
    std::shared_ptr<Mesh> mesh;
    try {
        mesh = std::make_shared<Mesh>(path);
    } catch (const std::exception &e) {
        std::cerr << "Failed to load mesh: " << e.what() << "\n";
        return 2;
    }

    std::optional<PolyField> field;
    std::optional<CutMesh> cutMesh;
    std::optional<MIQSolver> miqSolver;
    std::optional<CrossField> crossField;
    std::shared_ptr<SeparatrixTrace> separatrixTrace;
    std::shared_ptr<Mesh> delaunayMesh;
    std::shared_ptr<MedialAxis> medialAxis;
    Mode mode = Mode::Unselected;
    Phase phase = Phase::MeshOnly;
    MBOPhase mboPhase = MBOPhase::MeshOnly;
    MedialAxisPhase maPhase = MedialAxisPhase::MeshOnly;
    bool cWasDown = false;
    bool rWasDown = false;
    bool oneWasDown = false;
    bool twoWasDown = false;
    bool threeWasDown = false;
    bool singularitiesLogged = false;
    bool mboSteppingStarted = false;
    bool mboConverged = false;
    bool mboTracingStarted = false;
    bool mboTracingFinished = false;
    int mboStepCount = 0;
    
    // Console for timing output
    viewer::Console console;
    console.setMaxLines(8);

    // Log initial mesh info
    {
        std::ostringstream oss;
        oss << "Loaded mesh: " << mesh->triangles.size() << " triangles, " 
            << mesh->vertices.size() << " vertices";
        console.log(oss.str());
    }

    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        return 3;
    }

    glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);
    glfwWindowHint(GLFW_SAMPLES, 8);
#ifdef __APPLE__
    glfwWindowHint(GLFW_COCOA_RETINA_FRAMEBUFFER, GLFW_TRUE);
#endif

    GLFWwindow *window = glfwCreateWindow(1200, 900, "CrossGen Viewer", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        std::cerr << "Failed to create window\n";
        return 4;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    viewer::Bounds B = viewer::computeBounds(*mesh);
    double dx = B.maxx - B.minx;
    double dy = B.maxy - B.miny;
    double ext = std::max(dx, dy);
    if (ext <= 0) ext = 1.0;
    double pad = 0.1 * ext;

    viewer::ViewState view;
    view.cx = 0.5 * (B.minx + B.maxx);
    view.cy = 0.5 * (B.miny + B.maxy);
    view.baseW = (B.maxx - B.minx) + 2.0 * pad;
    view.baseH = (B.maxy - B.miny) + 2.0 * pad;
    if (view.baseW <= 0.0) view.baseW = 1.0;
    if (view.baseH <= 0.0) view.baseH = 1.0;
    view.zoom = 1.0;

    glfwGetFramebufferSize(window, &view.fbw, &view.fbh);
    glfwSetWindowUserPointer(window, &view);
    viewer::applyOrtho(view);
    viewer::installInteractionCallbacks(window);

    double avgEdge = viewer::averageTriangleEdgeLength(*mesh);
    double scale = 0.7 * avgEdge;

    std::cerr << "[Viewer] Phase " << phaseName(phase) << " (press '1' for PolyVector mode, '2' for MBO mode, '3' for Medial Axis mode)\n";

    glEnable(GL_MULTISAMPLE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glShadeModel(GL_SMOOTH);

    using Clock = std::chrono::high_resolution_clock;

    while (!glfwWindowShouldClose(window)) {
        // Reset (edge-triggered) - press 'r' to restart as if freshly loaded
        {
            bool rDown = (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS);
            if (rDown && !rWasDown) {
                // Reset all computed data
                field.reset();
                cutMesh.reset();
                miqSolver.reset();
                crossField.reset();
                separatrixTrace.reset();
                delaunayMesh.reset();
                medialAxis.reset();

                // Reset mode and phase state
                mode = Mode::Unselected;
                phase = Phase::MeshOnly;
                mboPhase = MBOPhase::MeshOnly;
                maPhase = MedialAxisPhase::MeshOnly;

                // Reset flags
                singularitiesLogged = false;
                mboSteppingStarted = false;
                mboConverged = false;
                mboTracingStarted = false;
                mboTracingFinished = false;
                mboStepCount = 0;

                // Reset view to original mesh bounds
                view.cx = 0.5 * (B.minx + B.maxx);
                view.cy = 0.5 * (B.miny + B.maxy);
                view.baseW = (B.maxx - B.minx) + 2.0 * pad;
                view.baseH = (B.maxy - B.miny) + 2.0 * pad;
                if (view.baseW <= 0.0) view.baseW = 1.0;
                if (view.baseH <= 0.0) view.baseH = 1.0;
                view.zoom = 1.0;
                viewer::applyOrtho(view);

                // Reset console
                console.clear();
                console.setMaxLines(8);
                {
                    std::ostringstream oss;
                    oss << "Loaded mesh: " << mesh->triangles.size() << " triangles, "
                        << mesh->vertices.size() << " vertices";
                    console.log(oss.str());
                }
                console.log("[Reset] Restarted viewer.");
                std::cerr << "[Viewer] Reset. Phase " << phaseName(phase)
                          << " (press '1' for PolyVector mode, '2' for MBO mode, '3' for Medial Axis mode)\n";
            }
            rWasDown = rDown;
        }

        // Mode selection (edge-triggered) - only when mode is unselected and in MeshOnly phase
        if (mode == Mode::Unselected && phase == Phase::MeshOnly) {
            bool oneDown = (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS);
            bool twoDown = (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS);
            bool threeDown = (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS);
            
            if (oneDown && !oneWasDown) {
                mode = Mode::PolyVector;
                std::cerr << "[Viewer] Selected mode: " << modeName(mode) << " (press 'c' to advance)\n";
                console.log("Selected mode: PolyVector");
            }
            if (twoDown && !twoWasDown) {
                mode = Mode::MBO;
                std::cerr << "[Viewer] Selected mode: " << modeName(mode) << " (press 'c' to advance)\n";
                console.log("Selected mode: MBO");
            }
            if (threeDown && !threeWasDown) {
                mode = Mode::MedialAxis;
                std::cerr << "[Viewer] Selected mode: " << modeName(mode) << " (press 'c' to advance)\n";
                console.log("Selected mode: Medial Axis");
            }
            
            oneWasDown = oneDown;
            twoWasDown = twoDown;
            threeWasDown = threeDown;
        }

        // Phase progression (edge-triggered) - only when mode is selected
        bool cDown = (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS);
        if (cDown && !cWasDown && mode != Mode::Unselected) {
            if (mode == Mode::PolyVector) {
                Phase old = phase;
                phase = nextPhase(phase);
                if (phase != old) {
                    std::cerr << "[Viewer] Phase " << phaseName(phase) << "\n";
                }
            } else if (mode == Mode::MBO) {
                MBOPhase old = mboPhase;
                mboPhase = nextMBOPhase(mboPhase);
                if (mboPhase != old) {
                    std::cerr << "[Viewer] MBO Phase " << mboPhaseName(mboPhase) << "\n";
                }
            } else if (mode == Mode::MedialAxis) {
                MedialAxisPhase old = maPhase;
                maPhase = nextMedialAxisPhase(maPhase);
                if (maPhase != old) {
                    std::cerr << "[Viewer] Medial Axis Phase " << medialAxisPhaseName(maPhase) << "\n";
                }
            }
        }
        cWasDown = cDown;

        // MBO mode: Initialize CrossField when entering CrossField phase
        if (mode == Mode::MBO && mboPhase >= MBOPhase::CrossField && !crossField.has_value()) {
            auto t0 = Clock::now();
            crossField.emplace(mesh);
            crossField->initialize(1); // method 1: random initialization for non-boundary vertices
            auto t1 = Clock::now();
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            console.log("[MBO] Initialized CrossField: " + formatMs(ms));
        }

        // MBO mode: Run stepping iterations when in Stepping phase
        if (mode == Mode::MBO && mboPhase == MBOPhase::Stepping && crossField.has_value() && !mboSteppingStarted) {
            mboSteppingStarted = true;
            mboStepCount = 0;
            console.log("[MBO] Starting " + std::to_string(MBO_MAX_STEPS) + " iterations...");
        }
        
        // Run 10 steps at a time and redraw
        if (mode == Mode::MBO && mboPhase == MBOPhase::Stepping && mboSteppingStarted && mboStepCount < MBO_MAX_STEPS && !mboConverged) {
            // Run 10 steps
            double nv = static_cast<double>(crossField->u_k.size());
            // stop if error < 2*n * 1e-4 (this is a very loose threshold just to prevent unnecessary stepping after convergence, since the viewer is not meant for precise timing/benchmarking)
            for (int i = 0; i < 2 && mboStepCount < MBO_MAX_STEPS; ++i) {
                crossField->step();
                mboStepCount++;
                if (crossField->error < 2.0 * nv * 1e-7) {
                    console.log("[MBO] Convergence reached at step " + std::to_string(mboStepCount) + " with error " + std::to_string(crossField->error));
                    mboConverged = true;
                    break;
                }
            }

            crossField->computeSingularities(); // update singularities for visualization during stepping (not just at the end)
            
            // Log progress
            std::ostringstream oss;
            oss << "[MBO] Step " << mboStepCount << "/" << MBO_MAX_STEPS;
            console.log(oss.str());
            
            // Render and sleep to show progress
            glClearColor(0.1f, 0.1f, 0.12f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glDisable(GL_DEPTH_TEST);
            
            viewer::drawMesh(*mesh);
            viewer::drawVertexCrossFieldUK(*mesh, *crossField, scale);
            
            // Draw singularities at triangle centroids
            double ballRadius = 0.5 * avgEdge;
            for (const auto &sig : crossField->singularTriangles) {
                int triIdx = sig.first;
                double crossIndex = sig.second;
                if (triIdx < 0 || triIdx >= static_cast<int>(mesh->triangles.size())) continue;
                
                // Compute triangle centroid
                const Triangle &tri = mesh->triangles[triIdx];
                const Point &p0 = mesh->vertices[tri[0]];
                const Point &p1 = mesh->vertices[tri[1]];
                const Point &p2 = mesh->vertices[tri[2]];
                Point centroid = {(p0[0] + p1[0] + p2[0]) / 3.0,
                                  (p0[1] + p1[1] + p2[1]) / 3.0};
                
                // Color based on index: blue for +1/4, red for -1/4
                if (crossIndex > 0) {
                    viewer::drawDisk3D(centroid, ballRadius, 0.2f, 0.2f, 0.95f);
                } else {
                    viewer::drawDisk3D(centroid, ballRadius, 0.95f, 0.2f, 0.2f);
                }
            }
            
            console.draw(window, 55.0f);
            viewer::drawTextOverlay(window, "MBO stepping in progress...\npress 'q' to quit", 10.0f, 20.0f, 0.8f, 0.8f, 0.8f);
            
            glfwSwapBuffers(window);
            glfwPollEvents();
            
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
            
            // Skip normal rendering this frame
            continue;
        }

        // MBO mode: Lazily construct SeparatrixTrace when entering Separatrices or Trace phase
        if (mode == Mode::MBO && mboPhase >= MBOPhase::Separatrices && crossField.has_value() && !separatrixTrace) {
            auto t0 = Clock::now();
            // Wrap the existing CrossField in a shared_ptr with a no-op deleter (ownership stays with the optional)
            auto cfPtr = std::shared_ptr<CrossField>(&*crossField, [](CrossField*){});
            separatrixTrace = std::make_shared<SeparatrixTrace>(cfPtr, false);
            auto t1 = Clock::now();
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            std::ostringstream oss;
            oss << "[Separatrices] Initialized " << separatrixTrace->separatrices.size()
                << " separatrices from " << separatrixTrace->singularities.size()
                << " singularities: " << formatMs(ms);
            console.log(oss.str());
        }

        // MBO mode: Trace phase - step separatrices one iteration at a time with animation
        if (mode == Mode::MBO && mboPhase == MBOPhase::Trace && separatrixTrace && !mboTracingFinished) {
            if (!mboTracingStarted) {
                mboTracingStarted = true;
                console.log("[Trace] Starting separatrix tracing...");
            }

            // Run one stepAndCheck iteration
            separatrixTrace->stepAndCheck();

            if (separatrixTrace->finishedTracing) {
                mboTracingFinished = true;
                console.log("[Trace] Tracing complete.");
            }

            // Render: mesh + separatrices
            glClearColor(0.1f, 0.1f, 0.12f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glDisable(GL_DEPTH_TEST);

            viewer::drawMesh(*mesh);

            // Draw separatrices: red if active, green if inactive
            glLineWidth(3.0f);
            for (const auto &sep : separatrixTrace->separatrices) {
                if (sep.path.size() < 2) continue;
                if (sep.active) {
                    glColor3f(0.95f, 0.1f, 0.1f);   // red for active
                } else {
                    glColor3f(0.1f, 0.9f, 0.2f);     // green for finished
                }
                glBegin(GL_LINE_STRIP);
                for (const auto &tp : sep.path) {
                    glVertex2d(tp.global_pos[0], tp.global_pos[1]);
                }
                glEnd();
            }
            glLineWidth(1.0f);

            console.draw(window, 55.0f);
            viewer::drawTextOverlay(window, "Tracing separatrices...\npress 'q' to quit", 10.0f, 20.0f, 0.8f, 0.8f, 0.8f);

            glfwSwapBuffers(window);
            glfwPollEvents();

            std::this_thread::sleep_for(std::chrono::milliseconds(10));

            // Skip normal rendering this frame
            continue;
        }

        // Medial Axis mode: Lazily re-triangulate boundary with Delaunay
        if (mode == Mode::MedialAxis && maPhase >= MedialAxisPhase::DelaunayMesh && !delaunayMesh) {
            auto t0 = Clock::now();

            // 1) Collect ordered boundary loops from the original mesh boundary edges
            //    Build adjacency for boundary vertices along boundary edges
            std::unordered_map<int, std::vector<int>> bAdj;
            for (int beIdx : mesh->boundaryEdges) {
                int a = mesh->edges[beIdx][0];
                int b = mesh->edges[beIdx][1];
                bAdj[a].push_back(b);
                bAdj[b].push_back(a);
            }

            // Walk boundary loops
            std::unordered_set<int> visited;
            std::vector<std::vector<int>> loops; // vertex index loops
            for (int bv : mesh->boundaryVertices) {
                if (visited.count(bv)) continue;
                std::vector<int> loop;
                int prev = -1;
                int curr = bv;
                while (true) {
                    visited.insert(curr);
                    loop.push_back(curr);
                    int next = -1;
                    for (int nb : bAdj[curr]) {
                        if (nb != prev && !visited.count(nb)) {
                            next = nb;
                            break;
                        }
                    }
                    if (next == -1) break; // loop closed or dead end
                    prev = curr;
                    curr = next;
                }
                if (loop.size() >= 3) {
                    loops.push_back(std::move(loop));
                }
            }

            if (loops.empty()) {
                console.log("[MedialAxis] No boundary loops found!");
            } else {
                // 2) Build TriangleMesher input from boundary loops
                //    Collect unique boundary vertices and remap indices
                std::unordered_map<int, int> vertRemap; // old index -> new index
                std::vector<std::array<double, 2>> vertlist;

                for (const auto &loop : loops) {
                    for (int vi : loop) {
                        if (!vertRemap.count(vi)) {
                            int newIdx = static_cast<int>(vertlist.size());
                            vertRemap[vi] = newIdx;
                            vertlist.push_back(mesh->vertices[vi]);
                        }
                    }
                }

                std::vector<std::vector<std::array<int, 2>>> segment_loops;
                std::vector<int> loopTypes;

                for (size_t li = 0; li < loops.size(); ++li) {
                    const auto &loop = loops[li];
                    std::vector<std::array<int, 2>> segments;
                    for (size_t i = 0; i < loop.size(); ++i) {
                        int a = vertRemap[loop[i]];
                        int b = vertRemap[loop[(i + 1) % loop.size()]];
                        segments.push_back({a, b});
                    }
                    segment_loops.push_back(std::move(segments));
                    // First loop is exterior, subsequent loops are holes
                    loopTypes.push_back(li == 0 ? 0 : 1);
                }

                // 3) Triangulate with just_delaunay = true
                triangle_wrapper::TriangleMesher2D::Options opts;
                opts.just_delaunay = true;
                triangle_wrapper::TriangleMesher2D mesher(opts);

                triangle_wrapper::TriangleMesher2D::MeshInput input;
                input.vertlist = vertlist;
                input.segment_loops = segment_loops;
                input.type = loopTypes;
                input.h = 0.0; // no refinement

                auto output = mesher.triangulate(input);

                // 4) Build new Mesh from the Delaunay output
                delaunayMesh = std::make_shared<Mesh>(output.verts, output.triangles);

                auto t1 = Clock::now();
                double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
                std::ostringstream oss;
                oss << "[MedialAxis] Delaunay re-triangulation: "
                    << delaunayMesh->triangles.size() << " triangles, "
                    << delaunayMesh->vertices.size() << " vertices: " << formatMs(ms);
                console.log(oss.str());
            }
        }

        // Medial Axis mode: Lazily compute medial axis
        if (mode == Mode::MedialAxis && maPhase >= MedialAxisPhase::MedialAxis && delaunayMesh && !medialAxis) {
            auto t0 = Clock::now();
            medialAxis = std::make_shared<MedialAxis>(delaunayMesh);
            int rawVerts = static_cast<int>(medialAxis->medialVertices.size());
            int rawEdges = static_cast<int>(medialAxis->medialEdges.size());
            medialAxis->deduplicateMedialVertices();
            auto t1 = Clock::now();
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            std::ostringstream oss;
            oss << "[MedialAxis] Computed medial axis: "
                << rawVerts << " -> " << medialAxis->medialVertices.size() << " vertices, "
                << rawEdges << " -> " << medialAxis->medialEdges.size() << " edges: " << formatMs(ms);
            console.log(oss.str());
        }

        // PolyVector mode: Lazily compute data when entering phases (with timing).
        if (mode == Mode::PolyVector && phase >= Phase::CrossField && !field.has_value()) {
            auto t0 = Clock::now();
            field.emplace(mesh);
            field->solveForPolyCoeffs();
            auto t1 = Clock::now();
            // convertToFieldVectors() also computes singularities internally
            field->convertToFieldVectors();
            auto t2 = Clock::now();
            double msCoeffs = std::chrono::duration<double, std::milli>(t1 - t0).count();
            double msField = std::chrono::duration<double, std::milli>(t2 - t1).count();
            console.log("[CrossField] Solved poly-coeffs: " + formatMs(msCoeffs));
            console.log("[CrossField] Converted to field vectors: " + formatMs(msField));
        }
        if (mode == Mode::PolyVector && phase >= Phase::Singularities && field.has_value() && !singularitiesLogged) {
            // Re-compute singularities to get accurate timing (they were computed in convertToFieldVectors)
            auto t0 = Clock::now();
            field->computeUSingularities();
            auto t1 = Clock::now();
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            std::ostringstream oss;
            oss << "[Singularities] Found " << field->uSingularities.size() 
                << " singularities: " << formatMs(ms);
            console.log(oss.str());
            singularitiesLogged = true;
        }
        if (mode == Mode::PolyVector && phase >= Phase::CutSeams && field.has_value() && !cutMesh.has_value()) {
            auto t0 = Clock::now();
            cutMesh.emplace(*field);
            auto t1 = Clock::now();
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

            std::ostringstream oss;
            oss << "[CutSeams] Generated " << cutMesh->getCutEdges().size() 
                << " cut edges: " << formatMs(ms);
            console.log(oss.str());

            std::cerr << "[Viewer] #tri=" << mesh->triangles.size() << " #vtx=" << mesh->vertices.size()
                      << " | uSingularities=" << field->uSingularities.size() << " | cutEdges="
                      << cutMesh->getCutEdges().size() << " | singularityPathCutEdges="
                      << cutMesh->getSingularityPathCutEdges().size() << "\n";

            if (!field->uSingularities.empty() && cutMesh->getSingularityPathCutEdges().empty()) {
                std::cerr
                    << "[Viewer] Note: no singularity->boundary path cuts were added (they may already lie on the boundary/cut graph)."
                    << " Falling back to showing all cut edges.\n";
            }
        }
        if (mode == Mode::PolyVector && phase >= Phase::UVMesh && cutMesh.has_value() && !miqSolver.has_value()) {
            auto t0 = Clock::now();
            miqSolver.emplace(*cutMesh);
            // Parameters: gradientSize, stiffness, directRound, iter, localIter, doRound, singularityRound, boundaryFeatures
            miqSolver->solve(100.0, 5.0, false, 10, 5000, true, true, true);
            auto t1 = Clock::now();
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

            int flips = miqSolver->numFlips();
            const auto &UV = miqSolver->getUV();

            std::ostringstream oss;
            oss << "[MIQ] Computed UV mesh: " << UV.rows() << " vertices, "
                << flips << " flips: " << formatMs(ms);
            console.log(oss.str());

            std::cerr << "[Viewer] MIQ parametrization: " << UV.rows() << " UV vertices, "
                      << flips << " flipped triangles\n";

            // Initialize view for UV space (reset pan/zoom for new coordinate system)
            viewer::computeUVMeshBounds(*miqSolver, view.cx, view.cy, view.baseW, view.baseH);
            view.zoom = 1.0;
            viewer::applyOrtho(view);
        }

        glClearColor(0.1f, 0.1f, 0.12f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glDisable(GL_DEPTH_TEST);

        // Phase 5: UV mesh (MIQ) - special case: clear and draw only UV mesh
        if (mode == Mode::PolyVector && phase == Phase::UVMesh && miqSolver.has_value()) {
            viewer::drawUVMesh(*miqSolver);
            // Draw singularities on UV mesh with same coloring
            if (field.has_value() && cutMesh.has_value()) {
                // Use a radius proportional to the UV mesh extent
                double uvRadius = 0.8; // fixed size that looks good on integer grid
                viewer::drawSingularitiesOnUV(*miqSolver, *cutMesh, *field, uvRadius);
            }
        } else if (mode == Mode::MBO) {
            // MBO mode rendering
            viewer::drawMesh(*mesh);
            
            // Draw crossfield on vertices if initialized
            if (mboPhase >= MBOPhase::CrossField && crossField.has_value()) {
                if (mboPhase >= MBOPhase::Stepping && mboStepCount > 0) {
                    // Use u_k after stepping has started
                    viewer::drawVertexCrossFieldUK(*mesh, *crossField, scale);
                } else {
                    // Use u_k_prev for initial display
                    viewer::drawVertexCrossField(*mesh, *crossField, scale);
                }
                
                // Draw singularities at triangle centroids (not in Separatrices phase)
                if (mboPhase < MBOPhase::Separatrices) {
                double ballRadius = 0.5 * avgEdge;
                for (const auto &sig : crossField->singularTriangles) {
                    int triIdx = sig.first;
                    double crossIndex = sig.second;
                    if (triIdx < 0 || triIdx >= static_cast<int>(mesh->triangles.size())) continue;
                    
                    // Compute triangle centroid
                    const Triangle &tri = mesh->triangles[triIdx];
                    const Point &p0 = mesh->vertices[tri[0]];
                    const Point &p1 = mesh->vertices[tri[1]];
                    const Point &p2 = mesh->vertices[tri[2]];
                    Point centroid = {(p0[0] + p1[0] + p2[0]) / 3.0,
                                      (p0[1] + p1[1] + p2[1]) / 3.0};
                    
                    // Color based on index: blue for +1/4, red for -1/4
                    if (crossIndex > 0) {
                        viewer::drawDisk3D(centroid, ballRadius, 0.2f, 0.2f, 0.95f);
                    } else {
                        viewer::drawDisk3D(centroid, ballRadius, 0.95f, 0.2f, 0.2f);
                    }
                }
                }
            }

            // Draw separatrices: red if active, green if inactive
            if (mboPhase >= MBOPhase::Separatrices && separatrixTrace) {
                glLineWidth(3.0f);
                for (const auto &sep : separatrixTrace->separatrices) {
                    if (sep.path.size() < 2) continue;
                    if (sep.active) {
                        glColor3f(0.95f, 0.1f, 0.1f);   // red for active
                    } else {
                        glColor3f(0.1f, 0.9f, 0.2f);     // green for finished
                    }
                    glBegin(GL_LINE_STRIP);
                    for (const auto &tp : sep.path) {
                        glVertex2d(tp.global_pos[0], tp.global_pos[1]);
                    }
                    glEnd();
                }
                glLineWidth(1.0f);
            }
        } else if (mode == Mode::MedialAxis) {
            // Medial Axis mode rendering
            // Draw full mesh (Delaunay if available, otherwise original)
            if (delaunayMesh) {
                viewer::drawMesh(*delaunayMesh);
            } else {
                viewer::drawMesh(*mesh);
            }

            // Draw medial axis overlay
            if (maPhase >= MedialAxisPhase::MedialAxis && medialAxis) {
                double ballRadius_ma = avgEdge / 5.0;
                viewer::drawMedialAxis(*medialAxis, ballRadius_ma);

                // Draw red circles at sharp corner vertices
                for (int i = 0; i < static_cast<int>(medialAxis->sharpVertices.size()); ++i) {
                    if (medialAxis->sharpVertices[i]) {
                        const Point &p = medialAxis->mesh->vertices[i];
                        viewer::drawDisk3D(p, ballRadius_ma, 0.95f, 0.2f, 0.2f);
                    }
                }
            }
        } else {
            // PolyVector mode rendering (phases 1-4)
            // Phase 1: mesh
            if (phase >= Phase::MeshOnly) {
                viewer::drawMesh(*mesh);
            }

            // Phase 2: crossfield (only in phases 2 and 3, not in phase 4)
            if (phase >= Phase::CrossField && phase < Phase::CutSeams && field.has_value()) {
                viewer::drawField(*mesh, *field, scale);
            }

            // Phase 3: singularities
            if (phase >= Phase::Singularities && field.has_value()) {
                double ballRadius = 0.5 * avgEdge; // large, relative to mesh scale
                for (const auto &sig : field->uSingularities) {
                    int vid = sig.first;
                    int index4 = sig.second;
                    if (vid < 0 || vid >= static_cast<int>(mesh->vertices.size())) continue;
                    const Point &c = mesh->vertices[vid];

                    if (index4 == 1) {
                        viewer::drawDisk3D(c, ballRadius, 0.2f, 0.2f, 0.95f);
                    } else if (index4 == -1) {
                        viewer::drawDisk3D(c, ballRadius, 0.95f, 0.2f, 0.2f);
                    }
                }
            }

            // Phase 4: seam cuts
            if (phase >= Phase::CutSeams && cutMesh.has_value()) {
                // Draw U field (green) and V field (red) - single direction per triangle
                viewer::drawUField(*mesh, cutMesh->getUField(), scale);
                viewer::drawVField(*mesh, cutMesh->getVField(), scale);

                if (!cutMesh->getSingularityPathCutEdges().empty()) {
                    viewer::drawEdgeSetOnMesh(*mesh, cutMesh->getCutEdges(), 1.0f, 0.75f, 0.1f, 4.0f);
                } else {
                    viewer::drawEdgeSetOnMesh(*mesh, cutMesh->getCutEdges(), 1.0f, 0.2f, 0.9f, 3.5f);
                }
            }
        }

        // Draw console at top (below help text)
        console.draw(window, 55.0f);

        // Draw help text overlay based on mode
        if (mode == Mode::Unselected) {
            viewer::drawTextOverlay(window, "press '1' for PolyVector mode\npress '2' for MBO mode\npress '3' for Medial Axis mode\npress 'r' to restart\npress 'q' to quit", 10.0f, 20.0f, 0.8f, 0.8f, 0.8f);
        } else {
            viewer::drawTextOverlay(window, "press 'c' to continue\npress 'r' to restart\npress 'q' to quit", 10.0f, 20.0f, 0.8f, 0.8f, 0.8f);
        }

        glfwSwapBuffers(window);
        glfwPollEvents();

        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS ||
            glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) {
            glfwSetWindowShouldClose(window, 1);
        }
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
