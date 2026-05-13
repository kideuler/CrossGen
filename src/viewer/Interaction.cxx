#include "viewer/Interaction.hxx"
#include "viewer/GL.hxx"

#include <algorithm>
#include <cmath>

namespace viewer {

void computeWorldBox(const ViewState &vs, double &worldW, double &worldH) {
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

    double left   = tmp.cx - 0.5 * worldW;
    double right  = tmp.cx + 0.5 * worldW;
    double bottom = tmp.cy - 0.5 * worldH;
    double top    = tmp.cy + 0.5 * worldH;

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

void panView(ViewState &view, double dxPx, double dyPx) {
    int fbw = std::max(1, view.fbw);
    int fbh = std::max(1, view.fbh);

    double worldW = 1.0, worldH = 1.0;
    computeWorldBox(view, worldW, worldH);

    double sx = worldW / static_cast<double>(fbw);
    double sy = worldH / static_cast<double>(fbh);

    // Move center opposite to mouse movement; flip Y (window coords grow down)
    view.cx -= dxPx * sx;
    view.cy += dyPx * sy;
    // Caller must invoke applyOrtho() when the GL context is current (i.e. in paintGL).
}

void zoomView(ViewState &view, double scrollDelta, double cursorX, double cursorY) {
    int fbw = std::max(1, view.fbw);
    int fbh = std::max(1, view.fbh);

    double worldW = 1.0, worldH = 1.0;
    computeWorldBox(view, worldW, worldH);

    double left   = view.cx - 0.5 * worldW;
    double bottom = view.cy - 0.5 * worldH;
    // World position under cursor
    double wx = left   + (cursorX / fbw) * worldW;
    double wy = bottom + ((1.0 - cursorY / fbh)) * worldH;

    double zoomFactor = std::pow(0.9, scrollDelta); // positive delta => zoom in
    double newZoom = std::clamp(view.zoom * zoomFactor, 0.02, 50.0);

    ViewState tmp = view;
    tmp.zoom = newZoom;
    computeWorldBox(tmp, worldW, worldH);

    double nx = cursorX / fbw;
    double ny = 1.0 - (cursorY / fbh);
    double newLeft   = wx - nx * worldW;
    double newBottom = wy - ny * worldH;
    view.cx   = newLeft   + 0.5 * worldW;
    view.cy   = newBottom + 0.5 * worldH;
    view.zoom = newZoom;
    // Caller must invoke applyOrtho() when the GL context is current (i.e. in paintGL).
}

} // namespace viewer
