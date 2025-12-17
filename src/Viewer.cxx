// Minimal OpenGL viewer for Mesh + PolyField
// - Loads an OBJ via Mesh
// - Builds PolyField, converts to field vectors
// - Renders triangles and per-triangle arrows using GLFW + OpenGL

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

#include <Eigen/Dense>

#include "polyvec/PolyVectors.hxx"

// GLFW header: do not include GL headers before it when using GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

struct Bounds {
    double minx{0}, miny{0}, maxx{0}, maxy{0};
};

struct ViewState {
    // world-space center and zoom
    double cx{0.0}, cy{0.0};
    double baseW{1.0}, baseH{1.0}; // world box derived from mesh + padding
    double zoom{1.0};              // 1.0 => fit to baseW/baseH; smaller => zoom in
    int fbw{1}, fbh{1};
    struct Drag {
        bool active{false};
        double lastWX{0.0}, lastWY{0.0}; // legacy fields (unused after change)
        double lastX{0.0}, lastY{0.0};   // last cursor position in screen pixels
    } drag;
};

static Bounds computeBounds(const Mesh &m) {
    Bounds b; bool first = true;
    for (const auto &p : m.vertices) {
        if (first) { b.minx = b.maxx = p[0]; b.miny = b.maxy = p[1]; first = false; continue; }
        b.minx = std::min(b.minx, p[0]); b.maxx = std::max(b.maxx, p[0]);
        b.miny = std::min(b.miny, p[1]); b.maxy = std::max(b.maxy, p[1]);
    }
    if (first) { b.minx = b.miny = 0; b.maxx = b.maxy = 1; }
    return b;
}

static void drawMesh(const Mesh &m) {
    glColor3f(0.85f, 0.85f, 0.85f);
    glLineWidth(1.25f);
    glBegin(GL_LINES);
    for (const auto &tri : m.triangles) {
        const Point &a = m.vertices[tri[0]];
        const Point &b = m.vertices[tri[1]];
        const Point &c = m.vertices[tri[2]];
        glVertex2d(a[0], a[1]); glVertex2d(b[0], b[1]);
        glVertex2d(b[0], b[1]); glVertex2d(c[0], c[1]);
        glVertex2d(c[0], c[1]); glVertex2d(a[0], a[1]);
    }
    glEnd();
}

static Point triangleCentroid(const Mesh &m, const Triangle &t) {
    const Point &a = m.vertices[t[0]];
    const Point &b = m.vertices[t[1]];
    const Point &c = m.vertices[t[2]];
    return Point{ (a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0 };
}

static void drawArrow(const Point &p, const Point &dir, double scale, float r, float g, float b) {
    Point d{ dir[0] * scale, dir[1] * scale };
    Point q{ p[0] + d[0], p[1] + d[1] };
    glColor3f(r, g, b);
    glBegin(GL_LINES);
        glVertex2d(p[0], p[1]);
        glVertex2d(q[0], q[1]);
    glEnd();
    // tiny head
    Point ortho{ -d[1], d[0] };
    double head = 0.2 * scale;
    double olen = std::sqrt(ortho[0]*ortho[0] + ortho[1]*ortho[1]);
    if (olen > 1e-12) { ortho[0] /= olen; ortho[1] /= olen; }
    Point h1{ q[0] - 0.3 * d[0] + head * ortho[0], q[1] - 0.3 * d[1] + head * ortho[1] };
    Point h2{ q[0] - 0.3 * d[0] - head * ortho[0], q[1] - 0.3 * d[1] - head * ortho[1] };
    glBegin(GL_LINES);
        glVertex2d(q[0], q[1]); glVertex2d(h1[0], h1[1]);
        glVertex2d(q[0], q[1]); glVertex2d(h2[0], h2[1]);
    glEnd();
}

static void drawField(const Mesh &m, const PolyField &field, double scale) {
    glLineWidth(1.5f);
    const float br = 0.2f, bg = 0.2f, bb = 0.95f; // unified blue color
    for (size_t i = 0; i < m.triangles.size(); ++i) {
        Point c = triangleCentroid(m, m.triangles[i]);
        if (i < field.field.size()) {
            const Point &u = field.field[i].u;
            const Point &v = field.field[i].v;
            // draw u and its opposite
            drawArrow(c, u, scale, br, bg, bb);
            Point uopp{ -u[0], -u[1] };
            drawArrow(c, uopp, scale, br, bg, bb);
            // draw v and its opposite
            drawArrow(c, v, scale, br, bg, bb);
            Point vopp{ -v[0], -v[1] };
            drawArrow(c, vopp, scale, br, bg, bb);
        }
    }
}

static double averageTriangleEdgeLength(const Mesh &m) {
    if (m.triangles.empty()) return 1.0;
    auto elen = [&](const Point &p, const Point &q){
        double dx = q[0] - p[0];
        double dy = q[1] - p[1];
        return std::sqrt(dx*dx + dy*dy);
    };
    double sum = 0.0;
    for (const auto &t : m.triangles) {
        const Point &a = m.vertices[t[0]];
        const Point &b = m.vertices[t[1]];
        const Point &c = m.vertices[t[2]];
        double l0 = elen(a,b);
        double l1 = elen(b,c);
        double l2 = elen(c,a);
        sum += (l0 + l1 + l2) / 3.0;
    }
    return sum / static_cast<double>(m.triangles.size());
}

int main(int argc, char** argv) {
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

    PolyField field(mesh);
    field.solveForPolyCoeffs();
    field.convertToFieldVectors();

    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        return 3;
    }
    glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);
    // Anti-aliasing via MSAA for smoother lines
    glfwWindowHint(GLFW_SAMPLES, 8);
    // Prefer Retina framebuffer on macOS for higher resolution
    #ifdef __APPLE__
    glfwWindowHint(GLFW_COCOA_RETINA_FRAMEBUFFER, GLFW_TRUE);
    #endif
    GLFWwindow* window = glfwCreateWindow(1200, 900, "CrossGen Viewer", nullptr, nullptr);
    if (!window) { glfwTerminate(); std::cerr << "Failed to create window\n"; return 4; }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    Bounds B = computeBounds(mesh);
    double dx = B.maxx - B.minx;
    double dy = B.maxy - B.miny;
    double ext = std::max(dx, dy);
    if (ext <= 0) ext = 1.0;
    double pad = 0.1 * ext;

    ViewState view;
    view.cx = 0.5 * (B.minx + B.maxx);
    view.cy = 0.5 * (B.miny + B.maxy);
    view.baseW = (B.maxx - B.minx) + 2.0 * pad;
    view.baseH = (B.maxy - B.miny) + 2.0 * pad;
    if (view.baseW <= 0.0) view.baseW = 1.0;
    if (view.baseH <= 0.0) view.baseH = 1.0;
    view.zoom = 1.0; // start fitted

    auto updateOrtho = [&](int w, int h){
        view.fbw = (w > 0) ? w : 1;
        view.fbh = (h > 0) ? h : 1;
        glViewport(0, 0, view.fbw, view.fbh);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        double viewAspect = static_cast<double>(view.fbw) / static_cast<double>(view.fbh);
        double worldW = view.baseW * view.zoom;
        double worldH = view.baseH * view.zoom;
        double worldAspect = worldW / worldH;
        // Preserve aspect by expanding the shorter dimension
        if (viewAspect > worldAspect) {
            worldW = worldH * viewAspect;
        } else {
            worldH = worldW / viewAspect;
        }
        double left = view.cx - 0.5 * worldW;
        double right = view.cx + 0.5 * worldW;
        double bottom = view.cy - 0.5 * worldH;
        double top = view.cy + 0.5 * worldH;
        glOrtho(left, right, bottom, top, -1, 1);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    };

    glfwGetFramebufferSize(window, &view.fbw, &view.fbh);
    updateOrtho(view.fbw, view.fbh);

    glfwSetWindowUserPointer(window, &view);
    glfwSetFramebufferSizeCallback(window, [](GLFWwindow* win, int w, int h){
        auto *vs = reinterpret_cast<ViewState*>(glfwGetWindowUserPointer(win));
        if (!vs) return;
        // Reconstruct a compatible update function
        glViewport(0, 0, std::max(1, w), std::max(1, h));
        vs->fbw = std::max(1, w);
        vs->fbh = std::max(1, h);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        double viewAspect = static_cast<double>(vs->fbw) / static_cast<double>(vs->fbh);
        double worldW = vs->baseW * vs->zoom;
        double worldH = vs->baseH * vs->zoom;
        double worldAspect = worldW / worldH;
        if (viewAspect > worldAspect) worldW = worldH * viewAspect; else worldH = worldW / viewAspect;
        double left = vs->cx - 0.5 * worldW;
        double right = vs->cx + 0.5 * worldW;
        double bottom = vs->cy - 0.5 * worldH;
        double top = vs->cy + 0.5 * worldH;
        glOrtho(left, right, bottom, top, -1, 1);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    });

    // Helper to get world coords under screen point for current view
    auto screenToWorld = [&](double sx, double sy){
        double viewAspect = static_cast<double>(view.fbw) / std::max(1, view.fbh);
        double worldW = view.baseW * view.zoom;
        double worldH = view.baseH * view.zoom;
        double worldAspect = worldW / worldH;
        if (viewAspect > worldAspect) worldW = worldH * viewAspect; else worldH = worldW / viewAspect;
        double left = view.cx - 0.5 * worldW;
        double bottom = view.cy - 0.5 * worldH;
        double nx = sx / static_cast<double>(std::max(1, view.fbw));
        double ny = 1.0 - (sy / static_cast<double>(std::max(1, view.fbh)));
        double wx = left + nx * worldW;
        double wy = bottom + ny * worldH;
        return std::pair<double,double>(wx, wy);
    };
    glfwSetMouseButtonCallback(window, [](GLFWwindow* win, int button, int action, int /*mods*/){
        auto *vs = reinterpret_cast<ViewState*>(glfwGetWindowUserPointer(win));
        if (!vs) return;
        if (button == GLFW_MOUSE_BUTTON_RIGHT) {
            if (action == GLFW_PRESS) {
                vs->drag.active = true;
                double x, y; glfwGetCursorPos(win, &x, &y);
                vs->drag.lastX = x; vs->drag.lastY = y;
            } else if (action == GLFW_RELEASE) {
                vs->drag.active = false;
            }
        }
    });
    glfwSetCursorPosCallback(window, [](GLFWwindow* win, double x, double y){
        auto *vs = reinterpret_cast<ViewState*>(glfwGetWindowUserPointer(win));
        if (!vs || !vs->drag.active) return;
        // Compute pixel delta
        double dxPx = x - vs->drag.lastX;
        double dyPx = y - vs->drag.lastY;
        vs->drag.lastX = x; vs->drag.lastY = y;

        // Convert to world delta using current world-per-pixel scale
        int fbw = std::max(1, vs->fbw);
        int fbh = std::max(1, vs->fbh);
        double viewAspect = static_cast<double>(fbw) / static_cast<double>(fbh);
        double worldW = vs->baseW * vs->zoom;
        double worldH = vs->baseH * vs->zoom;
        double worldAspect = worldW / worldH;
        if (viewAspect > worldAspect) worldW = worldH * viewAspect; else worldH = worldW / viewAspect;
        double sx = worldW / static_cast<double>(fbw);
        double sy = worldH / static_cast<double>(fbh);
        // Move center opposite to mouse movement; flip Y (window coords grow down)
        vs->cx -= dxPx * sx;
        vs->cy += dyPx * sy;
        // update projection
        glViewport(0, 0, fbw, fbh);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        double leftN = vs->cx - 0.5 * worldW;
        double rightN = vs->cx + 0.5 * worldW;
        double bottomN = vs->cy - 0.5 * worldH;
        double topN = vs->cy + 0.5 * worldH;
        glOrtho(leftN, rightN, bottomN, topN, -1, 1);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    });

    glfwSetScrollCallback(window, [](GLFWwindow* win, double /*xoff*/, double yoff){
        auto *vs = reinterpret_cast<ViewState*>(glfwGetWindowUserPointer(win));
        if (!vs) return;
        // zoom around cursor position
        double mx, my; glfwGetCursorPos(win, &mx, &my);
        int fbw, fbh; glfwGetFramebufferSize(win, &fbw, &fbh);
        double viewAspect = static_cast<double>(fbw) / std::max(1, fbh);
        double worldW = vs->baseW * vs->zoom;
        double worldH = vs->baseH * vs->zoom;
        double worldAspect = worldW / worldH;
        if (viewAspect > worldAspect) worldW = worldH * viewAspect; else worldH = worldW / viewAspect;
        double left = vs->cx - 0.5 * worldW;
        double bottom = vs->cy - 0.5 * worldH;
        double wx = left + (mx / std::max(1, fbw)) * worldW;
        double wy = bottom + ((1.0 - (my / std::max(1, fbh))) * worldH);

        // adjust zoom factor
        double zoomFactor = std::pow(0.9, yoff); // scroll up (positive yoff) => zoom in
        double newZoom = std::clamp(vs->zoom * zoomFactor, 0.02, 50.0);
        // keep cursor position stable in world coords: recompute center
        double newWorldW = vs->baseW * newZoom;
        double newWorldH = vs->baseH * newZoom;
        if (viewAspect > (newWorldW / newWorldH)) newWorldW = newWorldH * viewAspect; else newWorldH = newWorldW / viewAspect;
        double nx = mx / std::max(1, fbw);
        double ny = 1.0 - (my / std::max(1, fbh));
        double newLeft = wx - nx * newWorldW;
        double newBottom = wy - ny * newWorldH;
        vs->cx = newLeft + 0.5 * newWorldW;
        vs->cy = newBottom + 0.5 * newWorldH;
        vs->zoom = newZoom;

        // update projection
        glViewport(0, 0, std::max(1, fbw), std::max(1, fbh));
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        double leftN = vs->cx - 0.5 * newWorldW;
        double rightN = vs->cx + 0.5 * newWorldW;
        double bottomN = vs->cy - 0.5 * newWorldH;
        double topN = vs->cy + 0.5 * newWorldH;
        glOrtho(leftN, rightN, bottomN, topN, -1, 1);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    });

    // Estimate vector scale from average triangle size (edge length)
    double avgEdge = averageTriangleEdgeLength(mesh);
    // Use a fraction so arrows fit within typical triangles
    double scale = 0.7*avgEdge;

    // Smoothing and anti-aliasing for better visuals
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

    while (!glfwWindowShouldClose(window)) {
        glClearColor(0.1f, 0.1f, 0.12f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glDisable(GL_DEPTH_TEST);

    drawMesh(mesh);
    drawField(mesh, field, scale);

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
