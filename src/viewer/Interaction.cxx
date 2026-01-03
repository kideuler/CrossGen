#include "viewer/Interaction.hxx"

#include <algorithm>
#include <cmath>

namespace viewer {

static void computeWorldBox(const ViewState &vs, double &worldW, double &worldH) {
    double viewAspect = static_cast<double>(vs.fbw) / static_cast<double>(std::max(1, vs.fbh));
    worldW = vs.baseW * vs.zoom;
    worldH = vs.baseH * vs.zoom;
    double worldAspect = worldW / worldH;

    // Preserve aspect by expanding the shorter dimension
    if (viewAspect > worldAspect) {
        worldW = worldH * viewAspect;
    } else {
        worldH = worldW / viewAspect;
    }
}

void applyOrtho(const ViewState &view) {
    ViewState tmp = view;
    tmp.fbw = std::max(1, tmp.fbw);
    tmp.fbh = std::max(1, tmp.fbh);

    double worldW = 1.0, worldH = 1.0;
    computeWorldBox(tmp, worldW, worldH);

    double left = tmp.cx - 0.5 * worldW;
    double right = tmp.cx + 0.5 * worldW;
    double bottom = tmp.cy - 0.5 * worldH;
    double top = tmp.cy + 0.5 * worldH;

    glViewport(0, 0, tmp.fbw, tmp.fbh);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(left, right, bottom, top, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void resizeAndApplyOrtho(ViewState &view, int fbw, int fbh) {
    view.fbw = std::max(1, fbw);
    view.fbh = std::max(1, fbh);
    applyOrtho(view);
}

void installInteractionCallbacks(GLFWwindow *window) {
    glfwSetFramebufferSizeCallback(window, [](GLFWwindow *win, int w, int h) {
        auto *vs = reinterpret_cast<ViewState *>(glfwGetWindowUserPointer(win));
        if (!vs) return;
        resizeAndApplyOrtho(*vs, w, h);
    });

    glfwSetMouseButtonCallback(window, [](GLFWwindow *win, int button, int action, int /*mods*/) {
        auto *vs = reinterpret_cast<ViewState *>(glfwGetWindowUserPointer(win));
        if (!vs) return;
        if (button == GLFW_MOUSE_BUTTON_RIGHT) {
            if (action == GLFW_PRESS) {
                vs->drag.active = true;
                double x, y;
                glfwGetCursorPos(win, &x, &y);
                vs->drag.lastX = x;
                vs->drag.lastY = y;
            } else if (action == GLFW_RELEASE) {
                vs->drag.active = false;
            }
        }
    });

    glfwSetCursorPosCallback(window, [](GLFWwindow *win, double x, double y) {
        auto *vs = reinterpret_cast<ViewState *>(glfwGetWindowUserPointer(win));
        if (!vs || !vs->drag.active) return;

        double dxPx = x - vs->drag.lastX;
        double dyPx = y - vs->drag.lastY;
        vs->drag.lastX = x;
        vs->drag.lastY = y;

        // Convert to world delta using current world-per-pixel scale
        int fbw = std::max(1, vs->fbw);
        int fbh = std::max(1, vs->fbh);

        double worldW = 1.0, worldH = 1.0;
        ViewState tmp = *vs;
        tmp.fbw = fbw;
        tmp.fbh = fbh;
        computeWorldBox(tmp, worldW, worldH);

        double sx = worldW / static_cast<double>(fbw);
        double sy = worldH / static_cast<double>(fbh);

        // Move center opposite to mouse movement; flip Y (window coords grow down)
        vs->cx -= dxPx * sx;
        vs->cy += dyPx * sy;

        applyOrtho(*vs);
    });

    glfwSetScrollCallback(window, [](GLFWwindow *win, double /*xoff*/, double yoff) {
        auto *vs = reinterpret_cast<ViewState *>(glfwGetWindowUserPointer(win));
        if (!vs) return;

        // Zoom around cursor position
        double mx, my;
        glfwGetCursorPos(win, &mx, &my);
        int fbw, fbh;
        glfwGetFramebufferSize(win, &fbw, &fbh);
        fbw = std::max(1, fbw);
        fbh = std::max(1, fbh);

        ViewState tmp = *vs;
        tmp.fbw = fbw;
        tmp.fbh = fbh;

        double worldW = 1.0, worldH = 1.0;
        computeWorldBox(tmp, worldW, worldH);

        double left = vs->cx - 0.5 * worldW;
        double bottom = vs->cy - 0.5 * worldH;
        double wx = left + (mx / fbw) * worldW;
        double wy = bottom + ((1.0 - (my / fbh)) * worldH);

        double zoomFactor = std::pow(0.9, yoff); // scroll up => zoom in
        double newZoom = std::clamp(vs->zoom * zoomFactor, 0.02, 50.0);

        // Keep cursor position stable in world coords: recompute center
        tmp.zoom = newZoom;
        computeWorldBox(tmp, worldW, worldH);

        double nx = mx / fbw;
        double ny = 1.0 - (my / fbh);
        double newLeft = wx - nx * worldW;
        double newBottom = wy - ny * worldH;
        vs->cx = newLeft + 0.5 * worldW;
        vs->cy = newBottom + 0.5 * worldH;
        vs->zoom = newZoom;

        resizeAndApplyOrtho(*vs, fbw, fbh);
    });
}

} // namespace viewer
