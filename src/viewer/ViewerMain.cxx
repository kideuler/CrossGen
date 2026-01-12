// Viewer executable entry point (split out from the old monolithic Viewer.cxx).

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>

#include "viewer/Geometry.hxx"
#include "viewer/Interaction.hxx"
#include "viewer/Render.hxx"

#include "polyvec/CutMesh.hxx"
#include "polyvec/IntegerGridMap.hxx"
#include "polyvec/PolyVectors.hxx"

namespace {

enum class Phase {
    MeshOnly = 1,
    CrossField = 2,
    Singularities = 3,
    CutSeams = 4,
    UVMesh = 5,
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

const char *phaseName(Phase p) {
    switch (p) {
        case Phase::MeshOnly: return "1) mesh";
        case Phase::CrossField: return "2) crossfield";
        case Phase::Singularities: return "3) singularities";
        case Phase::CutSeams: return "4) cut seams";
        case Phase::UVMesh: return "5) uv mesh";
    }
    return "?";
}

// Format duration in milliseconds with 2 decimal places
std::string formatMs(double ms) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << ms << " ms";
    return oss.str();
}

} // namespace

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: Viewer <mesh.obj>\n";
        return 1;
    }

    std::string path = argv[1];
    Mesh mesh;
    try {
        mesh = Mesh(path);
    } catch (const std::exception &e) {
        std::cerr << "Failed to load mesh: " << e.what() << "\n";
        return 2;
    }

    std::optional<PolyField> field;
    std::optional<CutMesh> cutMesh;
    std::optional<IntegerGridMap> igm;
    Phase phase = Phase::MeshOnly;
    bool cWasDown = false;
    bool singularitiesLogged = false;
    
    // Console for timing output
    viewer::Console console;
    console.setMaxLines(8);

    // Log initial mesh info
    {
        std::ostringstream oss;
        oss << "Loaded mesh: " << mesh.triangles.size() << " triangles, " 
            << mesh.vertices.size() << " vertices";
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

    viewer::Bounds B = viewer::computeBounds(mesh);
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

    double avgEdge = viewer::averageTriangleEdgeLength(mesh);
    double scale = 0.7 * avgEdge;

    std::cerr << "[Viewer] Phase " << phaseName(phase) << " (press 'c' to advance)\n";

    glEnable(GL_MULTISAMPLE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glShadeModel(GL_SMOOTH);

    using Clock = std::chrono::high_resolution_clock;

    while (!glfwWindowShouldClose(window)) {
        // Phase progression (edge-triggered)
        bool cDown = (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS);
        if (cDown && !cWasDown) {
            Phase old = phase;
            phase = nextPhase(phase);
            if (phase != old) {
                std::cerr << "[Viewer] Phase " << phaseName(phase) << "\n";
            }
        }
        cWasDown = cDown;

        // Lazily compute data when entering phases (with timing).
        if (phase >= Phase::CrossField && !field.has_value()) {
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
        if (phase >= Phase::Singularities && field.has_value() && !singularitiesLogged) {
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
        if (phase >= Phase::CutSeams && field.has_value() && !cutMesh.has_value()) {
            auto t0 = Clock::now();
            cutMesh.emplace(*field);
            auto t1 = Clock::now();
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

            std::ostringstream oss;
            oss << "[CutSeams] Generated " << cutMesh->getCutEdges().size() 
                << " cut edges: " << formatMs(ms);
            console.log(oss.str());

            std::cerr << "[Viewer] #tri=" << mesh.triangles.size() << " #vtx=" << mesh.vertices.size()
                      << " | uSingularities=" << field->uSingularities.size() << " | cutEdges="
                      << cutMesh->getCutEdges().size() << " | singularityPathCutEdges="
                      << cutMesh->getSingularityPathCutEdges().size() << "\n";

            if (!field->uSingularities.empty() && cutMesh->getSingularityPathCutEdges().empty()) {
                std::cerr
                    << "[Viewer] Note: no singularity->boundary path cuts were added (they may already lie on the boundary/cut graph)."
                    << " Falling back to showing all cut edges.\n";
            }
        }
        if (phase >= Phase::UVMesh && cutMesh.has_value() && !igm.has_value()) {
            auto t0 = Clock::now();
            igm.emplace(*cutMesh);
            IntegerGridMap::Options opt;
            // Use auto-computed h: set h=0 to let it be computed based on mesh size.
            // target_grid_lines=100 means UV coordinates will span roughly [-50, 50],
            // giving many integer positions for singularity snapping.
            opt.h = 0.0;  // auto-compute based on target_grid_lines
            opt.target_grid_lines = 100;
            opt.stiffening_iterations = 10;
            opt.do_rounding = true;
            bool success = igm->solve(opt);
            auto t1 = Clock::now();
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

            if (success) {
                auto st = igm->computeStats();  // uses effective h from solve
                std::ostringstream oss;
                oss << "[UVMesh] Solved IGM: " << formatMs(ms) 
                    << " (flipped=" << st.flipped_triangles << ")";
                console.log(oss.str());

                // Print UV coordinates of singularities to terminal
                const auto &sings = cutMesh->getSingularities();
                const auto &uvCoords = igm->uv();
                const auto &origToCut = cutMesh->getOriginalToCutVertices();
                
                if (!sings.empty()) {
                    std::cerr << "\n[UVMesh] Singularity UV coordinates:\n";
                    std::cerr << std::fixed << std::setprecision(4);
                    for (size_t i = 0; i < sings.size(); ++i) {
                        int origVid = sings[i].first;
                        int index4 = sings[i].second;
                        const char* typeStr = (index4 == 1) ? "+1" : (index4 == -1) ? "-1" : "??";
                        
                        // Get UV from first cut vertex mapping to this original vertex
                        if (origVid >= 0 && origVid < static_cast<int>(origToCut.size()) && !origToCut[origVid].empty()) {
                            int cutVid = origToCut[origVid][0];
                            if (cutVid >= 0 && cutVid < static_cast<int>(uvCoords.size())) {
                                const Point &uv = uvCoords[cutVid];
                                std::cerr << "  Singularity " << i << " (index=" << typeStr << "): "
                                          << "UV = (" << uv[0] << ", " << uv[1] << ")\n";
                            }
                        }
                    }
                    std::cerr << std::endl;
                }
            } else {
                console.log("[UVMesh] IGM solve failed: " + formatMs(ms));
            }
        }

        glClearColor(0.1f, 0.1f, 0.12f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glDisable(GL_DEPTH_TEST);

        // Phase 5: UV mesh (clears everything and shows only UV coordinates)
        if (phase == Phase::UVMesh && igm.has_value()) {
            // Compute bounds of UV coordinates for proper view
            const auto &uvCoords = igm->uv();
            if (!uvCoords.empty()) {
                double uvMinX = uvCoords[0][0], uvMaxX = uvCoords[0][0];
                double uvMinY = uvCoords[0][1], uvMaxY = uvCoords[0][1];
                for (const auto &uv : uvCoords) {
                    uvMinX = std::min(uvMinX, uv[0]);
                    uvMaxX = std::max(uvMaxX, uv[0]);
                    uvMinY = std::min(uvMinY, uv[1]);
                    uvMaxY = std::max(uvMaxY, uv[1]);
                }
                double uvDx = uvMaxX - uvMinX;
                double uvDy = uvMaxY - uvMinY;
                double uvExt = std::max(uvDx, uvDy);
                if (uvExt <= 0) uvExt = 1.0;
                double uvPad = 0.1 * uvExt;

                // Set up orthographic projection for UV space
                glMatrixMode(GL_PROJECTION);
                glLoadIdentity();
                int fbw, fbh;
                glfwGetFramebufferSize(window, &fbw, &fbh);
                double aspect = static_cast<double>(fbw) / static_cast<double>(fbh);
                double uvCx = 0.5 * (uvMinX + uvMaxX);
                double uvCy = 0.5 * (uvMinY + uvMaxY);
                double uvW = (uvDx + 2.0 * uvPad);
                double uvH = (uvDy + 2.0 * uvPad);
                if (uvW <= 0.0) uvW = 1.0;
                if (uvH <= 0.0) uvH = 1.0;
                double halfW, halfH;
                if (aspect > uvW / uvH) {
                    halfH = uvH / 2.0;
                    halfW = halfH * aspect;
                } else {
                    halfW = uvW / 2.0;
                    halfH = halfW / aspect;
                }
                glOrtho(uvCx - halfW, uvCx + halfW, uvCy - halfH, uvCy + halfH, -1.0, 1.0);
                glMatrixMode(GL_MODELVIEW);
                glLoadIdentity();

                viewer::drawUVMesh(*igm, 0.3f, 0.85f, 0.95f, 1.25f);
            }
        } else {
            // Phases 1-4: draw on original mesh
            // Phase 1: mesh
            if (phase >= Phase::MeshOnly) {
                viewer::drawMesh(mesh);
            }

            // Phase 2: crossfield
            if (phase >= Phase::CrossField && field.has_value()) {
                viewer::drawField(mesh, *field, scale);
            }

            // Phase 3: singularities
            if (phase >= Phase::Singularities && field.has_value()) {
                double ballRadius = 0.5 * avgEdge; // large, relative to mesh scale
                for (const auto &sig : field->uSingularities) {
                    int vid = sig.first;
                    int index4 = sig.second;
                    if (vid < 0 || vid >= static_cast<int>(mesh.vertices.size())) continue;
                    const Point &c = mesh.vertices[vid];

                    if (index4 == 1) {
                        viewer::drawDisk3D(c, ballRadius, 0.2f, 0.2f, 0.95f);
                    } else if (index4 == -1) {
                        viewer::drawDisk3D(c, ballRadius, 0.95f, 0.2f, 0.2f);
                    }
                }
            }

            // Phase 4: seam cuts
            if (phase >= Phase::CutSeams && cutMesh.has_value()) {
                if (!cutMesh->getSingularityPathCutEdges().empty()) {
                    viewer::drawEdgeSetOnMesh(mesh, cutMesh->getCutEdges(), 1.0f, 0.75f, 0.1f, 4.0f);
                } else {
                    viewer::drawEdgeSetOnMesh(mesh, cutMesh->getCutEdges(), 1.0f, 0.2f, 0.9f, 3.5f);
                }
            }
        }

        // Draw console at top (below help text)
        console.draw(window, 55.0f);

        // Draw help text overlay
        viewer::drawTextOverlay(window, "press 'c' to continue\npress 'q' to quit", 10.0f, 20.0f, 0.8f, 0.8f, 0.8f);

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
