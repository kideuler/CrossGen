#pragma once

// Viewer-specific simple types (no OpenGL/GLFW headers here).

namespace viewer {

struct Bounds {
    double minx{0}, miny{0}, maxx{0}, maxy{0};
};

struct ViewState {
    // World-space center and zoom
    double cx{0.0}, cy{0.0};
    double baseW{1.0}, baseH{1.0}; // world box derived from mesh + padding
    double zoom{1.0};              // 1.0 => fit to baseW/baseH; smaller => zoom in
    int fbw{1}, fbh{1};

    struct Drag {
        bool active{false};
        double lastX{0.0}, lastY{0.0}; // last cursor position in screen pixels
    } drag;
};

} // namespace viewer
