// Viewer executable entry point (split out from the old monolithic Viewer.cxx).

#include <algorithm>
#include <iostream>
#include <optional>
#include <string>

#include "viewer/Geometry.hxx"
#include "viewer/Interaction.hxx"
#include "viewer/Render.hxx"

#include "polyvec/CutMesh.hxx"
#include "polyvec/PolyVectors.hxx"

namespace {

enum class Phase {
    MeshOnly = 1,
    CrossField = 2,
    Singularities = 3,
    CutSeams = 4,
};

Phase nextPhase(Phase p) {
    switch (p) {
        case Phase::MeshOnly: return Phase::CrossField;
        case Phase::CrossField: return Phase::Singularities;
        case Phase::Singularities: return Phase::CutSeams;
        case Phase::CutSeams: return Phase::CutSeams;
    }
    return Phase::CutSeams;
}

const char *phaseName(Phase p) {
    switch (p) {
        case Phase::MeshOnly: return "1) mesh";
        case Phase::CrossField: return "2) crossfield";
        case Phase::Singularities: return "3) singularities";
        case Phase::CutSeams: return "4) cut seams";
    }
    return "?";
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
    Phase phase = Phase::MeshOnly;
    bool cWasDown = false;

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

        // Lazily compute data when entering phases.
        if (phase >= Phase::CrossField && !field.has_value()) {
            field.emplace(mesh);
            field->solveForPolyCoeffs();
            field->convertToFieldVectors();
        }
        if (phase >= Phase::Singularities && field.has_value() && field->uSingularities.empty()) {
            field->computeUSingularities();
        }
        if (phase >= Phase::CutSeams && field.has_value() && !cutMesh.has_value()) {
            cutMesh.emplace(*field);

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

        glClearColor(0.1f, 0.1f, 0.12f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glDisable(GL_DEPTH_TEST);

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
                viewer::drawEdgeSetOnMesh(mesh, cutMesh->getSingularityPathCutEdges(), 1.0f, 0.75f, 0.1f, 4.0f);
            } else {
                viewer::drawEdgeSetOnMesh(mesh, cutMesh->getCutEdges(), 1.0f, 0.2f, 0.9f, 3.5f);
            }
        }

        glfwSwapBuffers(window);
        glfwPollEvents();

        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
            glfwSetWindowShouldClose(window, 1);
        }
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
