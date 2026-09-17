#include "viewer/Render.hxx"

#include "viewer/Export.hxx"
#include "viewer/GL.hxx"
#include "viewer/Geometry.hxx"
#include "viewer/Interaction.hxx"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#include <Eigen/Dense>

namespace viewer {

// ============================================================================
// Console implementation
// ============================================================================

void Console::log(const std::string &msg) {
    lines_.push_back(msg);
    // Keep only the last maxLines_ entries
    if (static_cast<int>(lines_.size()) > maxLines_) {
        lines_.erase(lines_.begin());
    }
}

void Console::clear() {
    lines_.clear();
}

void Console::draw(int fbw, int fbh, float startY) const {
    if (lines_.empty()) return;

    // Save current projection/modelview and set up screen-space orthographic
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, fbw, fbh, 0, -1, 1); // top-left origin

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    float scale = 2.0f * textScale();
    float charH = 8.0f * scale;
    float padding = 8.0f;
    float lineSpacing = charH + 2.0f;
    
    // Calculate console dimensions
    float consoleHeight = lines_.size() * lineSpacing + padding * 2;
    float consoleWidth = fbw - 20.0f;  // nearly full width with margins
    
    // Draw semi-transparent background
    viewer::color4f(0.0f, 0.0f, 0.0f, 0.6f);
    glBegin(GL_QUADS);
    glVertex2f(10.0f, startY - padding);
    glVertex2f(10.0f + consoleWidth, startY - padding);
    glVertex2f(10.0f + consoleWidth, startY + consoleHeight - padding);
    glVertex2f(10.0f, startY + consoleHeight - padding);
    glEnd();

    // Draw border
    viewer::color3f(0.3f, 0.6f, 0.3f);
    viewer::lineWidth(1.0f);
    glBegin(GL_LINE_LOOP);
    glVertex2f(10.0f, startY - padding);
    glVertex2f(10.0f + consoleWidth, startY - padding);
    glVertex2f(10.0f + consoleWidth, startY + consoleHeight - padding);
    glVertex2f(10.0f, startY + consoleHeight - padding);
    glEnd();

    // Restore matrices for text drawing
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();

    // Draw each line using drawTextOverlay (it sets up its own projection)
    float y = startY;
    for (const auto &line : lines_) {
        drawTextOverlay(fbw, fbh, line.c_str(), 18.0f, y, 0.4f, 0.9f, 0.4f);
        y += lineSpacing;
    }
}

// ============================================================================
// Mesh drawing
// ============================================================================

// Color for material id 0..9 from a fixed qualitative palette; anything
// beyond that gets a color hashed from the id, so it is stable across
// frames (no flicker) without needing a palette entry for every material.
static void materialColor(int matId, float &r, float &g, float &b) {
    static const float palette[10][3] = {
        {0.10f, 0.80f, 0.80f}, // 0 cyan
        {0.85f, 0.15f, 0.15f}, // 1 red
        {0.15f, 0.75f, 0.15f}, // 2 green
        {0.55f, 0.20f, 0.80f}, // 3 purple
        {0.95f, 0.55f, 0.10f}, // 4 orange
        {0.90f, 0.85f, 0.10f}, // 5 yellow
        {0.10f, 0.15f, 0.60f}, // 6 dark blue
        {0.55f, 0.05f, 0.05f}, // 7 dark red
        {0.05f, 0.40f, 0.05f}, // 8 dark green
        {0.30f, 0.05f, 0.45f}, // 9 dark purple
    };
    if (matId >= 0 && matId < 10) {
        r = palette[matId][0];
        g = palette[matId][1];
        b = palette[matId][2];
        return;
    }
    std::size_t h = std::hash<int>{}(matId);
    float hue = static_cast<float>(h % 360u) / 360.0f;
    float s = 0.65f, v = 0.9f;
    float i = std::floor(hue * 6.0f);
    float f = hue * 6.0f - i;
    float p = v * (1.0f - s);
    float q = v * (1.0f - f * s);
    float t = v * (1.0f - (1.0f - f) * s);
    switch (static_cast<int>(i) % 6) {
        case 0: r = v; g = t; b = p; break;
        case 1: r = q; g = v; b = p; break;
        case 2: r = p; g = v; b = t; break;
        case 3: r = p; g = q; b = v; break;
        case 4: r = t; g = p; b = v; break;
        default: r = v; g = p; b = q; break;
    }
}

void drawMesh(const Mesh &m) {
    viewer::lineWidth(1.25f);
    glBegin(GL_LINES);
    for (std::size_t e = 0; e < m.edges.size(); ++e) {
        int t0 = m.edgeTriangles[e][0];
        int t1 = m.edgeTriangles[e][1];
        float r, g, b;
        if (t1 < 0) {
            // Boundary edge: colored by its one incident triangle's material.
            materialColor(m.triangleMatId[t0], r, g, b);
        } else if (m.triangleMatId[t0] != m.triangleMatId[t1]) {
            // Interface edge, touching more than one material: highlight white.
            r = g = b = 1.0f;
        } else {
            materialColor(m.triangleMatId[t0], r, g, b);
        }
        viewer::color3f(r, g, b);
        const Point &a = m.vertices[m.edges[e][0]];
        const Point &b_ = m.vertices[m.edges[e][1]];
        glVertex2d(a[0], a[1]);
        glVertex2d(b_[0], b_[1]);
    }
    glEnd();
}

void drawMeshOverlay(const Mesh &m, float r, float g, float b, float a, float lineWidth) {
    viewer::color4f(r, g, b, a);
    viewer::lineWidth(lineWidth);
    glBegin(GL_LINES);
    for (const auto &tri : m.triangles) {
        const Point &p0 = m.vertices[tri[0]];
        const Point &p1 = m.vertices[tri[1]];
        const Point &p2 = m.vertices[tri[2]];
        glVertex2d(p0[0], p0[1]);
        glVertex2d(p1[0], p1[1]);
        glVertex2d(p1[0], p1[1]);
        glVertex2d(p2[0], p2[1]);
        glVertex2d(p2[0], p2[1]);
        glVertex2d(p0[0], p0[1]);
    }
    glEnd();
    viewer::lineWidth(1.0f);
}

void drawEdgeSetOnMesh(const Mesh &m,
                       const std::unordered_set<CutMesh::EdgeKey, CutMesh::EdgeKeyHash> &edges,
                       float r, float g, float b,
                       float lineWidth) {
    viewer::color3f(r, g, b);
    viewer::lineWidth(lineWidth);
    glBegin(GL_LINES);
    for (const auto &e : edges) {
        if (e.a < 0 || e.b < 0 || e.a >= static_cast<int>(m.vertices.size()) ||
            e.b >= static_cast<int>(m.vertices.size())) {
            continue;
        }
        const Point &pa = m.vertices[e.a];
        const Point &pb = m.vertices[e.b];
        glVertex2d(pa[0], pa[1]);
        glVertex2d(pb[0], pb[1]);
    }
    glEnd();
}

void drawArrow(const Point &p, const Point &dir, double scale, float r, float g, float b) {
    Point d{dir[0] * scale, dir[1] * scale};
    Point q{p[0] + d[0], p[1] + d[1]};
    viewer::color3f(r, g, b);
    glBegin(GL_LINES);
    glVertex2d(p[0], p[1]);
    glVertex2d(q[0], q[1]);
    glEnd();

    // Tiny head
    Point ortho{-d[1], d[0]};
    double head = 0.2 * scale;
    double olen = std::sqrt(ortho[0] * ortho[0] + ortho[1] * ortho[1]);
    if (olen > 1e-12) {
        ortho[0] /= olen;
        ortho[1] /= olen;
    }
    Point h1{q[0] - 0.3 * d[0] + head * ortho[0], q[1] - 0.3 * d[1] + head * ortho[1]};
    Point h2{q[0] - 0.3 * d[0] - head * ortho[0], q[1] - 0.3 * d[1] - head * ortho[1]};
    glBegin(GL_LINES);
    glVertex2d(q[0], q[1]);
    glVertex2d(h1[0], h1[1]);
    glVertex2d(q[0], q[1]);
    glVertex2d(h2[0], h2[1]);
    glEnd();
}

void drawField(const Mesh &m, const PolyField &field, double scale) {
    viewer::lineWidth(2.5f);
    const float br = 0.2f, bg = 0.2f, bb = 0.95f; // unified blue color

    for (size_t i = 0; i < m.triangles.size(); ++i) {
        Point c = triangleCentroid(m, m.triangles[i]);
        if (i < field.field.size()) {
            const Point &u = field.field[i].u;
            const Point &v = field.field[i].v;

            // draw u and its opposite
            drawArrow(c, u, scale, br, bg, bb);
            Point uopp{-u[0], -u[1]};
            drawArrow(c, uopp, scale, br, bg, bb);

            // draw v and its opposite
            drawArrow(c, v, scale, br, bg, bb);
            Point vopp{-v[0], -v[1]};
            drawArrow(c, vopp, scale, br, bg, bb);
        }
    }
}

void drawUField(const Mesh &m, const std::vector<Point> &uField, double scale) {
    viewer::lineWidth(2.5f);
    const float br = 0.2f, bg = 0.7f, bb = 0.2f; // green color for U field

    for (size_t i = 0; i < m.triangles.size(); ++i) {
        Point c = triangleCentroid(m, m.triangles[i]);
        if (i < uField.size()) {
            const Point &u = uField[i];
            drawArrow(c, u, scale, br, bg, bb);
        }
    }
}

void drawVField(const Mesh &m, const std::vector<Point> &vField, double scale) {
    viewer::lineWidth(2.5f);
    const float vr = 0.9f, vg = 0.2f, vb = 0.2f; // red color for V field

    for (size_t i = 0; i < m.triangles.size(); ++i) {
        Point c = triangleCentroid(m, m.triangles[i]);
        if (i < vField.size()) {
            const Point &v = vField[i];
            drawArrow(c, v, scale, vr, vg, vb);
        }
    }
}

void drawVertexCrossField(const Mesh &m, const CrossField &cf, double scale) {
    viewer::lineWidth(2.5f);
    const float br = 0.2f, bg = 0.2f, bb = 0.95f; // unified blue color

    const Eigen::VectorXcd &u_k_prev = cf.u_k_prev;
    
    for (size_t i = 0; i < m.vertices.size(); ++i) {
        if (static_cast<int>(i) >= u_k_prev.size()) continue;
        
        const Point &c = m.vertices[i];
        std::complex<double> u4 = u_k_prev[i];
        
        // Skip if magnitude is too small
        if (std::abs(u4) < 1e-14) continue;
        
        // Compute fourth root: u = u4^(1/4)
        // u4 = r * e^(i*theta), so u4^(1/4) = r^(1/4) * e^(i*theta/4)
        double r = std::abs(u4);
        double theta = std::arg(u4);
        double r4 = std::pow(r, 0.25);
        double theta4 = theta / 4.0;
        
        // The four directions are at angles: theta4, theta4 + pi/2, theta4 + pi, theta4 + 3*pi/2
        for (int k = 0; k < 4; ++k) {
            double angle = theta4 + k * M_PI_2;
            Point dir{r4 * std::cos(angle), r4 * std::sin(angle)};
            drawArrow(c, dir, scale, br, bg, bb);
        }
    }
}

void drawVertexCrossFieldUK(const Mesh &m, const CrossField &cf, double scale) {
    viewer::lineWidth(2.5f);
    const float br = 0.2f, bg = 0.2f, bb = 0.95f; // unified blue color

    const Eigen::VectorXcd &u_k = cf.u_k;
    
    for (size_t i = 0; i < m.vertices.size(); ++i) {
        if (static_cast<int>(i) >= u_k.size()) continue;
        
        const Point &c = m.vertices[i];
        std::complex<double> u4 = u_k[i];
        
        // Skip if magnitude is too small
        if (std::abs(u4) < 1e-14) continue;
        
        // Compute fourth root: u = u4^(1/4)
        // u4 = r * e^(i*theta), so u4^(1/4) = r^(1/4) * e^(i*theta/4)
        double r = std::abs(u4);
        double theta = std::arg(u4);
        double r4 = std::pow(r, 0.25);
        double theta4 = theta / 4.0;
        
        // The four directions are at angles: theta4, theta4 + pi/2, theta4 + pi, theta4 + 3*pi/2
        for (int k = 0; k < 4; ++k) {
            double angle = theta4 + k * M_PI_2;
            Point dir{r4 * std::cos(angle), r4 * std::sin(angle)};
            drawArrow(c, dir, scale, br, bg, bb);
        }
    }
}

void drawTriangleCrossField(const Mesh &m, const DualMBO &dualMBO, double scale) {
    viewer::lineWidth(2.5f);
    const float br = 0.45f, bg = 0.05f, bb = 0.55f; // dark purple

    const Eigen::VectorXcd &u_k = dualMBO.u_k;

    for (int t = 0; t < static_cast<int>(m.triangles.size()); ++t) {
        if (t >= u_k.size()) continue;

        std::complex<double> u4 = u_k[t];
        if (std::abs(u4) < 1e-14) continue;

        // Compute centroid of triangle
        const Triangle &tri = m.triangles[t];
        const Point &p0 = m.vertices[tri[0]];
        const Point &p1 = m.vertices[tri[1]];
        const Point &p2 = m.vertices[tri[2]];
        Point c = {(p0[0] + p1[0] + p2[0]) / 3.0,
                   (p0[1] + p1[1] + p2[1]) / 3.0};

        // Decode spin-4 angle: u4 = exp(4i*theta)
        double theta4 = std::arg(u4);
        double theta  = theta4 / 4.0;

        // Draw four arms of the cross
        for (int k = 0; k < 4; ++k) {
            double angle = theta + k * M_PI_2;
            Point dir{std::cos(angle), std::sin(angle)};
            drawArrow(c, dir, scale, br, bg, bb);
        }
    }
}

void drawDisk3D(const Point &center, double radius, float baseR, float baseG, float baseB, int segments) {
    // Fake light direction in view space (towards viewer, slightly to top-right)
    Eigen::Vector3d L = Eigen::Vector3d(0.4, 0.4, 0.8).normalized();

    glBegin(GL_TRIANGLE_FAN);
    // Center normal pointing out of screen (0,0,1)
    double ndotl_center = std::max(0.0, L[2]);
    double shade_center = 0.3 + 0.7 * ndotl_center; // ambient + diffuse
    viewer::color3f(baseR * shade_center, baseG * shade_center, baseB * shade_center);
    glVertex2d(center[0], center[1]);

    // Rim vertices: compute per-vertex shading by mapping disk to sphere cap
    for (int i = 0; i <= segments; ++i) {
        double ang = (static_cast<double>(i) / segments) * 2.0 * M_PI;
        double dx = std::cos(ang);
        double dy = std::sin(ang);
        double x = center[0] + radius * dx;
        double y = center[1] + radius * dy;

        // Map (dx,dy) on unit disk to hemisphere normal: z = sqrt(max(0, 1 - r^2))
        double r2 = dx * dx + dy * dy; // equals 1 on rim
        double z = std::sqrt(std::max(0.0, 1.0 - r2));
        Eigen::Vector3d N(dx, dy, z);
        N.normalize();
        double ndotl = std::max(0.0, N.dot(L));
        double shade = 0.25 + 0.75 * ndotl; // ambient + diffuse
        viewer::color3f(baseR * shade, baseG * shade, baseB * shade);
        glVertex2d(x, y);
    }
    glEnd();

    // Subtle outline to enhance 3D look
    viewer::lineWidth(1.0f);
    glBegin(GL_LINE_LOOP);
    for (int i = 0; i < segments; ++i) {
        double ang = (static_cast<double>(i) / segments) * 2.0 * M_PI;
        double x = center[0] + radius * std::cos(ang);
        double y = center[1] + radius * std::sin(ang);
        viewer::color3f(baseR * 0.5f, baseG * 0.5f, baseB * 0.5f);
        glVertex2d(x, y);
    }
    glEnd();
}

// Simple 5x7 bitmap font for basic ASCII characters (space through ~)
// Each character is stored as 7 bytes, one per row, with 5 bits per row.
static const unsigned char kFont5x7[95][7] = {
    {0x00,0x00,0x00,0x00,0x00,0x00,0x00}, // ' '
    {0x04,0x04,0x04,0x04,0x04,0x00,0x04}, // '!'
    {0x0A,0x0A,0x00,0x00,0x00,0x00,0x00}, // '"'
    {0x0A,0x0A,0x1F,0x0A,0x1F,0x0A,0x0A}, // '#'
    {0x04,0x0F,0x14,0x0E,0x05,0x1E,0x04}, // '$'
    {0x18,0x19,0x02,0x04,0x08,0x13,0x03}, // '%'
    {0x08,0x14,0x14,0x08,0x15,0x12,0x0D}, // '&'
    {0x04,0x04,0x00,0x00,0x00,0x00,0x00}, // '\''
    {0x02,0x04,0x08,0x08,0x08,0x04,0x02}, // '('
    {0x08,0x04,0x02,0x02,0x02,0x04,0x08}, // ')'
    {0x00,0x04,0x15,0x0E,0x15,0x04,0x00}, // '*'
    {0x00,0x04,0x04,0x1F,0x04,0x04,0x00}, // '+'
    {0x00,0x00,0x00,0x00,0x00,0x04,0x08}, // ','
    {0x00,0x00,0x00,0x1F,0x00,0x00,0x00}, // '-'
    {0x00,0x00,0x00,0x00,0x00,0x00,0x04}, // '.'
    {0x00,0x01,0x02,0x04,0x08,0x10,0x00}, // '/'
    {0x0E,0x11,0x13,0x15,0x19,0x11,0x0E}, // '0'
    {0x04,0x0C,0x04,0x04,0x04,0x04,0x0E}, // '1'
    {0x0E,0x11,0x01,0x06,0x08,0x10,0x1F}, // '2'
    {0x0E,0x11,0x01,0x06,0x01,0x11,0x0E}, // '3'
    {0x02,0x06,0x0A,0x12,0x1F,0x02,0x02}, // '4'
    {0x1F,0x10,0x1E,0x01,0x01,0x11,0x0E}, // '5'
    {0x06,0x08,0x10,0x1E,0x11,0x11,0x0E}, // '6'
    {0x1F,0x01,0x02,0x04,0x08,0x08,0x08}, // '7'
    {0x0E,0x11,0x11,0x0E,0x11,0x11,0x0E}, // '8'
    {0x0E,0x11,0x11,0x0F,0x01,0x02,0x0C}, // '9'
    {0x00,0x00,0x04,0x00,0x00,0x04,0x00}, // ':'
    {0x00,0x00,0x04,0x00,0x00,0x04,0x08}, // ';'
    {0x02,0x04,0x08,0x10,0x08,0x04,0x02}, // '<'
    {0x00,0x00,0x1F,0x00,0x1F,0x00,0x00}, // '='
    {0x08,0x04,0x02,0x01,0x02,0x04,0x08}, // '>'
    {0x0E,0x11,0x01,0x06,0x04,0x00,0x04}, // '?'
    {0x0E,0x11,0x17,0x15,0x17,0x10,0x0E}, // '@'
    {0x0E,0x11,0x11,0x1F,0x11,0x11,0x11}, // 'A'
    {0x1E,0x11,0x11,0x1E,0x11,0x11,0x1E}, // 'B'
    {0x0E,0x11,0x10,0x10,0x10,0x11,0x0E}, // 'C'
    {0x1E,0x11,0x11,0x11,0x11,0x11,0x1E}, // 'D'
    {0x1F,0x10,0x10,0x1E,0x10,0x10,0x1F}, // 'E'
    {0x1F,0x10,0x10,0x1E,0x10,0x10,0x10}, // 'F'
    {0x0E,0x11,0x10,0x17,0x11,0x11,0x0E}, // 'G'
    {0x11,0x11,0x11,0x1F,0x11,0x11,0x11}, // 'H'
    {0x0E,0x04,0x04,0x04,0x04,0x04,0x0E}, // 'I'
    {0x07,0x02,0x02,0x02,0x02,0x12,0x0C}, // 'J'
    {0x11,0x12,0x14,0x18,0x14,0x12,0x11}, // 'K'
    {0x10,0x10,0x10,0x10,0x10,0x10,0x1F}, // 'L'
    {0x11,0x1B,0x15,0x15,0x11,0x11,0x11}, // 'M'
    {0x11,0x19,0x15,0x13,0x11,0x11,0x11}, // 'N'
    {0x0E,0x11,0x11,0x11,0x11,0x11,0x0E}, // 'O'
    {0x1E,0x11,0x11,0x1E,0x10,0x10,0x10}, // 'P'
    {0x0E,0x11,0x11,0x11,0x15,0x12,0x0D}, // 'Q'
    {0x1E,0x11,0x11,0x1E,0x14,0x12,0x11}, // 'R'
    {0x0E,0x11,0x10,0x0E,0x01,0x11,0x0E}, // 'S'
    {0x1F,0x04,0x04,0x04,0x04,0x04,0x04}, // 'T'
    {0x11,0x11,0x11,0x11,0x11,0x11,0x0E}, // 'U'
    {0x11,0x11,0x11,0x11,0x11,0x0A,0x04}, // 'V'
    {0x11,0x11,0x11,0x15,0x15,0x1B,0x11}, // 'W'
    {0x11,0x11,0x0A,0x04,0x0A,0x11,0x11}, // 'X'
    {0x11,0x11,0x0A,0x04,0x04,0x04,0x04}, // 'Y'
    {0x1F,0x01,0x02,0x04,0x08,0x10,0x1F}, // 'Z'
    {0x0E,0x08,0x08,0x08,0x08,0x08,0x0E}, // '['
    {0x00,0x10,0x08,0x04,0x02,0x01,0x00}, // '\\'
    {0x0E,0x02,0x02,0x02,0x02,0x02,0x0E}, // ']'
    {0x04,0x0A,0x11,0x00,0x00,0x00,0x00}, // '^'
    {0x00,0x00,0x00,0x00,0x00,0x00,0x1F}, // '_'
    {0x08,0x04,0x00,0x00,0x00,0x00,0x00}, // '`'
    {0x00,0x00,0x0E,0x01,0x0F,0x11,0x0F}, // 'a'
    {0x10,0x10,0x1E,0x11,0x11,0x11,0x1E}, // 'b'
    {0x00,0x00,0x0E,0x11,0x10,0x11,0x0E}, // 'c'
    {0x01,0x01,0x0F,0x11,0x11,0x11,0x0F}, // 'd'
    {0x00,0x00,0x0E,0x11,0x1F,0x10,0x0E}, // 'e'
    {0x06,0x08,0x1E,0x08,0x08,0x08,0x08}, // 'f'
    {0x00,0x00,0x0F,0x11,0x0F,0x01,0x0E}, // 'g'
    {0x10,0x10,0x1E,0x11,0x11,0x11,0x11}, // 'h'
    {0x04,0x00,0x0C,0x04,0x04,0x04,0x0E}, // 'i'
    {0x02,0x00,0x06,0x02,0x02,0x12,0x0C}, // 'j'
    {0x10,0x10,0x12,0x14,0x18,0x14,0x12}, // 'k'
    {0x0C,0x04,0x04,0x04,0x04,0x04,0x0E}, // 'l'
    {0x00,0x00,0x1A,0x15,0x15,0x11,0x11}, // 'm'
    {0x00,0x00,0x1E,0x11,0x11,0x11,0x11}, // 'n'
    {0x00,0x00,0x0E,0x11,0x11,0x11,0x0E}, // 'o'
    {0x00,0x00,0x1E,0x11,0x1E,0x10,0x10}, // 'p'
    {0x00,0x00,0x0F,0x11,0x0F,0x01,0x01}, // 'q'
    {0x00,0x00,0x16,0x19,0x10,0x10,0x10}, // 'r'
    {0x00,0x00,0x0F,0x10,0x0E,0x01,0x1E}, // 's'
    {0x08,0x08,0x1E,0x08,0x08,0x09,0x06}, // 't'
    {0x00,0x00,0x11,0x11,0x11,0x11,0x0F}, // 'u'
    {0x00,0x00,0x11,0x11,0x11,0x0A,0x04}, // 'v'
    {0x00,0x00,0x11,0x11,0x15,0x15,0x0A}, // 'w'
    {0x00,0x00,0x11,0x0A,0x04,0x0A,0x11}, // 'x'
    {0x00,0x00,0x11,0x11,0x0F,0x01,0x0E}, // 'y'
    {0x00,0x00,0x1F,0x02,0x04,0x08,0x1F}, // 'z'
    {0x02,0x04,0x04,0x08,0x04,0x04,0x02}, // '{'
    {0x04,0x04,0x04,0x04,0x04,0x04,0x04}, // '|'
    {0x08,0x04,0x04,0x02,0x04,0x04,0x08}, // '}'
    {0x00,0x00,0x08,0x15,0x02,0x00,0x00}, // '~'
};

void drawTextOverlay(int fbw, int fbh, const char *text, float x, float y, float r, float g, float b) {

    // Save current projection/modelview and set up screen-space orthographic
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, fbw, fbh, 0, -1, 1); // top-left origin

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    viewer::color3f(r, g, b);
    
    float scale = 2.0f * textScale(); // readability, times the export scale
    float charW = 6.0f * scale;  // 5 pixels + 1 spacing
    float charH = 8.0f * scale;  // 7 pixels + 1 spacing
    
    float curX = x;
    float curY = y;

    glBegin(GL_QUADS);
    for (const char *p = text; *p; ++p) {
        if (*p == '\n') {
            curX = x;
            curY += charH;
            continue;
        }

        int ch = static_cast<unsigned char>(*p);
        if (ch < 32 || ch > 126) ch = '?';
        int idx = ch - 32;

        // Draw character as filled quads for each pixel
        for (int row = 0; row < 7; ++row) {
            unsigned char rowBits = kFont5x7[idx][row];
            for (int col = 0; col < 5; ++col) {
                if (rowBits & (0x10 >> col)) {
                    float px = curX + col * scale;
                    float py = curY + row * scale;
                    glVertex2f(px, py);
                    glVertex2f(px + scale, py);
                    glVertex2f(px + scale, py + scale);
                    glVertex2f(px, py + scale);
                }
            }
        }
        curX += charW;
    }
    glEnd();

    // Restore matrices
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
}

#ifdef CROSSGEN_WITH_COMISO
void computeUVMeshBounds(const MIQSolver &miq, double &cx, double &cy, double &baseW, double &baseH) {
    const Eigen::MatrixXd &UV = miq.getUV();

    if (UV.rows() == 0) {
        cx = 0.0;
        cy = 0.0;
        baseW = 1.0;
        baseH = 1.0;
        return;
    }

    // Compute UV bounds
    double minU = UV.col(0).minCoeff();
    double maxU = UV.col(0).maxCoeff();
    double minV = UV.col(1).minCoeff();
    double maxV = UV.col(1).maxCoeff();

    double du = maxU - minU;
    double dv = maxV - minV;
    double ext = std::max(du, dv);
    if (ext <= 0) ext = 1.0;
    double pad = 0.1 * ext;

    cx = 0.5 * (minU + maxU);
    cy = 0.5 * (minV + maxV);
    baseW = du + 2.0 * pad;
    baseH = dv + 2.0 * pad;
    if (baseW <= 0.0) baseW = 1.0;
    if (baseH <= 0.0) baseH = 1.0;
}

void drawUVMesh(const MIQSolver &miq) {
    const Eigen::MatrixXd &UV = miq.getUV();
    const Eigen::MatrixXi &FUV = miq.getFUV();

    if (UV.rows() == 0 || FUV.rows() == 0) return;

    // Compute UV bounds for grid drawing
    double minU = UV.col(0).minCoeff();
    double maxU = UV.col(0).maxCoeff();
    double minV = UV.col(1).minCoeff();
    double maxV = UV.col(1).maxCoeff();

    double du = maxU - minU;
    double dv = maxV - minV;
    double ext = std::max(du, dv);
    if (ext <= 0) ext = 1.0;
    double pad = 0.1 * ext;

    // First pass: draw flipped triangles as filled red polygons
    viewer::color4f(0.9f, 0.2f, 0.2f, 0.7f); // red with some transparency
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < FUV.rows(); ++i) {
        if (!miq.isFlipped(i)) continue;

        int v0 = FUV(i, 0);
        int v1 = FUV(i, 1);
        int v2 = FUV(i, 2);

        if (v0 < 0 || v0 >= UV.rows() ||
            v1 < 0 || v1 >= UV.rows() ||
            v2 < 0 || v2 >= UV.rows()) {
            continue;
        }

        double u0 = UV(v0, 0), w0 = UV(v0, 1);
        double u1 = UV(v1, 0), w1 = UV(v1, 1);
        double u2 = UV(v2, 0), w2 = UV(v2, 1);

        glVertex2d(u0, w0);
        glVertex2d(u1, w1);
        glVertex2d(u2, w2);
    }
    glEnd();

    // Draw UV mesh edges
    viewer::color3f(0.3f, 0.8f, 0.9f); // cyan color for UV mesh
    viewer::lineWidth(1.5f);
    glBegin(GL_LINES);
    for (int i = 0; i < FUV.rows(); ++i) {
        int v0 = FUV(i, 0);
        int v1 = FUV(i, 1);
        int v2 = FUV(i, 2);

        if (v0 < 0 || v0 >= UV.rows() ||
            v1 < 0 || v1 >= UV.rows() ||
            v2 < 0 || v2 >= UV.rows()) {
            continue;
        }

        double u0 = UV(v0, 0), w0 = UV(v0, 1);
        double u1 = UV(v1, 0), w1 = UV(v1, 1);
        double u2 = UV(v2, 0), w2 = UV(v2, 1);

        // Edge 0-1
        glVertex2d(u0, w0);
        glVertex2d(u1, w1);
        // Edge 1-2
        glVertex2d(u1, w1);
        glVertex2d(u2, w2);
        // Edge 2-0
        glVertex2d(u2, w2);
        glVertex2d(u0, w0);
    }
    glEnd();

    // Draw integer grid lines for reference
    int gridMinU = static_cast<int>(std::floor(minU));
    int gridMaxU = static_cast<int>(std::ceil(maxU));
    int gridMinV = static_cast<int>(std::floor(minV));
    int gridMaxV = static_cast<int>(std::ceil(maxV));

    viewer::color4f(0.4f, 0.4f, 0.4f, 0.5f);
    viewer::lineWidth(1.0f);
    glBegin(GL_LINES);
    // Vertical lines (constant U)
    for (int u = gridMinU; u <= gridMaxU; ++u) {
        glVertex2d(static_cast<double>(u), minV - pad);
        glVertex2d(static_cast<double>(u), maxV + pad);
    }
    // Horizontal lines (constant V)
    for (int v = gridMinV; v <= gridMaxV; ++v) {
        glVertex2d(minU - pad, static_cast<double>(v));
        glVertex2d(maxU + pad, static_cast<double>(v));
    }
    glEnd();
}

void drawSingularitiesOnUV(const MIQSolver &miq, const CutMesh &cutMesh,
                           const PolyField &field, double radius) {
    const Eigen::MatrixXd &UV = miq.getUV();
    const auto &origToCut = cutMesh.getOriginalToCutVertices();
    const auto &singularities = field.uSingularities;

    if (UV.rows() == 0) return;

    for (const auto &sig : singularities) {
        int origVid = sig.first;
        int index4 = sig.second;

        // Skip if original vertex index is out of range
        if (origVid < 0 || origVid >= static_cast<int>(origToCut.size())) continue;

        // Get all cut mesh vertices corresponding to this original vertex
        const auto &cutVerts = origToCut[origVid];
        if (cutVerts.empty()) continue;

        // Use the first cut vertex to get UV coordinates
        // (all copies of a singularity should map to the same UV location in a valid parametrization)
        int cutVid = cutVerts[0];
        if (cutVid < 0 || cutVid >= UV.rows()) continue;

        double u = UV(cutVid, 0);
        double v = UV(cutVid, 1);
        Point center{u, v};

        // Draw with same coloring as 3D view: blue for +1, red for -1
        if (index4 == 1) {
            drawDisk3D(center, radius, 0.2f, 0.2f, 0.95f);
        } else if (index4 == -1) {
            drawDisk3D(center, radius, 0.95f, 0.2f, 0.2f);
        }
    }
}
#endif // CROSSGEN_WITH_COMISO

void drawMedialAxis(const MedialAxis &ma, double vertexRadius) {
    // Draw medial axis edges (Voronoi dual edges) in orange
    viewer::color3f(1.0f, 0.6f, 0.1f);
    viewer::lineWidth(2.5f);
    glBegin(GL_LINES);
    for (const auto &edge : ma.medialEdges) {
        const Point &a = ma.medialVertices[edge[0]].coord;
        const Point &b = ma.medialVertices[edge[1]].coord;
        glVertex2d(a[0], a[1]);
        glVertex2d(b[0], b[1]);
    }
    glEnd();
    viewer::lineWidth(1.0f);

    // Draw medial axis vertices (circumcenters) as small cyan disks. Vertices
    // whose circumcenter fell outside the domain are drawn red instead: they
    // are not on the axis, and seeing them is the point.
    for (const auto &v : ma.medialVertices) {
        if (v.insideDomain) drawDisk3D(v.coord, vertexRadius, 0.1f, 0.85f, 0.85f);
        else                drawDisk3D(v.coord, vertexRadius, 0.95f, 0.1f, 0.1f);
    }
}

namespace {

void medialClassColor(MedialColor c, float alpha) {
    switch (c) {
        case MedialColor::Green:  viewer::color4f(0.15f, 0.85f, 0.25f, alpha); return;
        case MedialColor::Red:    viewer::color4f(0.95f, 0.25f, 0.20f, alpha); return;
        case MedialColor::Blue:   viewer::color4f(0.30f, 0.50f, 0.95f, alpha); return;
        case MedialColor::Purple: viewer::color4f(0.80f, 0.30f, 0.95f, alpha); return;
    }
}

} // namespace

void drawMedialTMesh(const MedialAxisTMesh &tm, double cornerRadius) {
    const MedialAxis &ma = *tm.axis;

    // ── Block fills ──
    //
    // One fan per block, from its centroid. A template block is star-shaped
    // about its centroid -- its four sides are two straight template edges
    // and at most two boundary runs, none of which turns back on itself --
    // so the fan is a valid triangulation and needs no ear clipping.
    for (const TMeshBlock &block : tm.blocks) {
        if (block.outline.size() < 3) continue;
        Point c{0.0, 0.0};
        for (const Point &p : block.outline) c = c + p;
        c = c / static_cast<double>(block.outline.size());

        medialClassColor(block.color, 0.28f);
        glBegin(GL_TRIANGLE_FAN);
        glVertex2d(c[0], c[1]);
        for (const Point &p : block.outline) glVertex2d(p[0], p[1]);
        glVertex2d(block.outline.front()[0], block.outline.front()[1]);
        glEnd();
    }

    // ── The downsampled axis in class colour, under the block walls ──
    viewer::lineWidth(3.0f);
    for (const MedialZone &zone : tm.zones) {
        if (zone.chain.size() < 2) continue;
        medialClassColor(zone.color, 1.0f);
        glBegin(GL_LINE_STRIP);
        for (int m : zone.chain) {
            const Point &p = ma.medialVertices[m].coord;
            glVertex2d(p[0], p[1]);
        }
        glEnd();
    }

    // ── The T-mesh blocking: every template block's closed outline ──
    viewer::lineWidth(2.2f);
    viewer::color4f(0.94f, 0.94f, 0.94f, 0.95f);
    for (const TMeshBlock &block : tm.blocks) {
        glBegin(GL_LINE_LOOP);
        for (const Point &p : block.outline) glVertex2d(p[0], p[1]);
        glEnd();
    }
    viewer::lineWidth(1.0f);

    // ── Block corners ──
    for (const TMeshBlock &block : tm.blocks) {
        for (const Point &p : block.corners) {
            drawDisk3D(p, cornerRadius, 0.95f, 0.9f, 0.85f);
        }
    }
}

namespace {

// Line segments per cell along a grid curve, as a floor rather than the
// measure: see gridSamples() below for why the cell count alone is not
// enough to decide how finely a grid line has to be walked.
constexpr int QUANT_CURVE_SUBDIV = 8;

// Ceiling on the samples one grid line may take, so a block sitting on a
// pathologically dense stretch of boundary cannot cost unbounded time.
constexpr int QUANT_CURVE_MAX_SAMPLES = 2000;

// Width of every line of the quantized picture -- grid lines, the outlines
// of blocks with no grid, and the outlines of components that never made it
// into the T-mesh. Kept in one place so those three cannot drift apart:
// they are all meant to read as the same decomposition.
constexpr float QUANT_LINE_WIDTH = 3.0f;

// Transfinite (Coons) interpolation between the four sides of a block.
// Along each side the formula collapses to that side's own geometry, so a
// block wall is drawn exactly as it curves and two blocks sharing it agree
// on it to the last point.
Point coons(const SideCurve &bottom, const SideCurve &top,
            const SideCurve &left, const SideCurve &right,
            double u, double v) {
    const Point p = bottom.at(u) * (1.0 - v) + top.at(u) * v +
                    left.at(v) * (1.0 - u) + right.at(v) * u;
    return p - (bottom.at(0.0) * ((1.0 - u) * (1.0 - v)) +
                bottom.at(1.0) * (u * (1.0 - v)) +
                top.at(0.0)    * ((1.0 - u) * v) +
                top.at(1.0)    * (u * v));
}

// Total turning of a side's polyline, unsigned: how far its direction
// rotates from one end to the other. Near zero for a straight or gently
// curved side; approaching 2*pi for a side that wraps a hole.
double sideWinding(const std::vector<Point> &pts) {
    double total = 0.0;
    Point prev{0.0, 0.0};
    bool havePrev = false;
    for (size_t i = 1; i < pts.size(); ++i) {
        const Point d = pts[i] - pts[i - 1];
        if (normP(d) <= 0.0) continue;
        if (havePrev) {
            total += std::atan2(prev[0] * d[1] - prev[1] * d[0],
                                prev[0] * d[0] + prev[1] * d[1]);
        }
        prev = d;
        havePrev = true;
    }
    return std::fabs(total);
}

// Above this winding, a pair of sides no longer looks anything like the
// straight rails of a quad and the Coons shape terms taken from the other
// pair stop being a correction and start being a lie (see below).
constexpr double QUANT_WINDING_LIMIT = 2.0 * M_PI / 3.0;

// The arc length behind SideCurve::at(t): the tick-interpolated position,
// as a fraction of the side's full length.
double arcFractionAt(const SideCurve &sc, double t) {
    const int n = sc.cells();
    if (n <= 0 || sc.length <= 0.0) return 0.0;
    const double s = std::min(std::max(t, 0.0), 1.0) * n;
    int i = static_cast<int>(s);
    if (i >= n) i = n - 1;
    const double f = s - i;
    return (sc.tickAt[i] + (sc.tickAt[i + 1] - sc.tickAt[i]) * f) / sc.length;
}

// One point of a ruled interior line between sides A and B of a winding
// pair: at parameter w in [0, 1] the line has slid from A's endpoint tick
// (arc fraction fa of A) to B's (fb of B), and sits at the blend of the two
// sides sampled at the *interpolated* fraction. A plain chord from tick to
// tick is not enough: on a block that wraps a hole the two ticks can sit at
// different angles around it, and the chord then cuts straight across.
// Sliding the fraction keeps the line inside the band the two sides bound,
// and it still lands exactly on both ticks.
Point ruledPoint(const SideCurve &a, const SideCurve &b, double fa, double fb,
                 double w) {
    const double f = fa * (1.0 - w) + fb * w;
    return a.atArc(f * a.length) * (1.0 - w) + b.atArc(f * b.length) * w;
}

// How finely to walk a grid line that spans `cells` cells and runs
// alongside the two sides `a` and `b`.
//
// The cell count on its own is the wrong measure. With xIdeal = 1 Stage II
// drives every edge toward its minimum, so `cells` is routinely 1, while a
// single block side can be a long run of the model boundary carrying
// hundreds of polyline points. Scaling by cells alone would then draw that
// side as QUANT_CURVE_SUBDIV straight segments -- and since the Coons
// formula collapses to the side's own geometry at u = 0 and u = 1, the
// outermost grid line *is* that side, so the block visibly stops following
// the boundary it was cut from and reads as a coarse polygon.
//
// The sides know their own true geometry, so ask them: enough samples to
// resolve every point they carry, and never fewer than the per-cell floor.
int gridSamples(int cells, const SideCurve &a, const SideCurve &b) {
    const int byCells = std::max(1, cells) * QUANT_CURVE_SUBDIV;
    const int byGeometry = static_cast<int>(a.pts.size() + b.pts.size());
    return std::clamp(std::max(byCells, byGeometry), 1, QUANT_CURVE_MAX_SAMPLES);
}

}  // namespace

namespace {

// Shared by drawQuantizedBlocks() and drawQuantizedLayout(): both quant
// structs (BlockQuant, QuadLayoutQuant) carry a QuantTMesh plus per-edge
// geometry and per-side reversal flags, and sideCurve() is overloaded on
// both, so the grid-drawing walk below only needs the type as a parameter.
//
// `vertexRadius` > 0 additionally marks every grid vertex -- every point
// where two grid lines cross, block corners and interior cell corners
// alike -- with a disk of that radius. They are all vertices of the final
// quad decomposition: a block-structure T-junction included, since with
// every edge quantized >= 1 the grids of the faces meeting there weld
// flush and the point is a regular vertex of the result.
template <typename Quant>
void drawQuantizedTMesh(const Quant &bq, double vertexRadius = 0.0) {
    const QuantTMesh &q = bq.tmesh;

    viewer::lineWidth(QUANT_LINE_WIDTH);
    viewer::color4f(0.35f, 0.62f, 0.98f, 0.95f);

    // Every edge of a face at its own full geometry: what to fall back on
    // whenever the face has no grid to draw. Drawing it beats leaving a
    // hole, and since it is the edge polylines themselves it follows the
    // boundary exactly.
    const auto outlineFace = [&](size_t f) {
        for (int s = 0; s < 4; ++s) {
            for (int e : q.faces[f].sides[s]) {
                glBegin(GL_LINE_STRIP);
                for (const Point &p : bq.edgeGeometry[e]) glVertex2d(p[0], p[1]);
                glEnd();
            }
        }
    };

    for (size_t f = 0; f < q.faces.size(); ++f) {
        const int nu = static_cast<int>(q.sideSum(static_cast<int>(f), 0));
        const int nv = static_cast<int>(q.sideSum(static_cast<int>(f), 1));

        // A block the quantization collapsed has no cells to draw, but it
        // is still part of the decomposition: outline it, so the region
        // reads as a block with no subdivision rather than as a hole.
        if (nu <= 0 || nv <= 0) {
            outlineFace(f);
            continue;
        }

        // The four sides on a common (u, v) frame: bottom and left run away
        // from corner 0, top and right toward corner 2, so opposite sides
        // are parametrized alike.
        const SideCurve bottom = sideCurve(bq, static_cast<int>(f), 0);
        const SideCurve right  = sideCurve(bq, static_cast<int>(f), 1);
        const SideCurve top    = sideCurve(bq, static_cast<int>(f), 2).reversed();
        const SideCurve left   = sideCurve(bq, static_cast<int>(f), 3).reversed();

        // An inconsistent quantization would leave opposite sides with
        // different cell counts, and then there is no grid to lay down --
        // but the block is still there, so outline it rather than dropping
        // it out of the picture without a word.
        if (bottom.cells() != nu || top.cells()  != nu ||
            right.cells()  != nv || left.cells() != nv) {
            outlineFace(f);
            continue;
        }

        const int sampleU = gridSamples(nu, bottom, top);
        const int sampleV = gridSamples(nv, left, right);

        // Which interpolant each family of interior lines gets.
        //
        // An iso-u line hangs its endpoints on bottom and top and takes its
        // shape through the interior from left and right (and vice versa
        // for iso-v). On a block that wraps a hole -- left and right two
        // near-full circles, bottom and top the short cut between them --
        // that borrowing inverts: an iso-v line's endpoints travel around
        // the ring with left and right, while its shape terms stay pinned
        // to the cut's position, and the line is dragged clear across the
        // hole and out of the block. So when a family's *endpoint* pair
        // winds far more than its *shape* pair, drop the shape terms and
        // rule straight between the endpoints: for the ring that is a
        // radial chord, which is exactly right. The other family keeps
        // Coons -- its lines run along the winding and the winding sides
        // are then its shape terms, which is what carries them around.
        const double windB = sideWinding(bottom.pts), windT = sideWinding(top.pts);
        const double windL = sideWinding(left.pts), windR = sideWinding(right.pts);
        const bool ruledU = std::max(windB, windT) > QUANT_WINDING_LIMIT &&
                            std::max(windB, windT) > std::max(windL, windR);
        const bool ruledV = std::max(windL, windR) > QUANT_WINDING_LIMIT &&
                            std::max(windL, windR) > std::max(windB, windT);

        // The outermost grid lines are the block's own sides: the Coons
        // formula collapses to left/right at u = 0/1 and to bottom/top at
        // v = 0/1. So draw those four from their polylines verbatim rather
        // than resampling them.
        //
        // This is not just cheaper, it is the only way to get them right.
        // Sampling walks a side at even steps of arc length, and even
        // steps need not land on the side's own vertices -- so wherever a
        // side turns a genuine corner between two of them, the drawn line
        // cuts it off, and no amount of extra sampling fixes it because
        // the corner is a real kink and not a resolution artifact. Traced
        // layouts have such corners: a block side that runs along the
        // model boundary inherits every corner of it. Drawn from the
        // polyline the side is exact, kinks included, so a block that sits
        // on the boundary conforms to it exactly.
        const auto strip = [](const std::vector<Point> &pts) {
            glBegin(GL_LINE_STRIP);
            for (const Point &p : pts) glVertex2d(p[0], p[1]);
            glEnd();
        };

        for (int i = 0; i <= nu; ++i) {
            if (i == 0)  { strip(left.pts);  continue; }
            if (i == nu) { strip(right.pts); continue; }
            const double u = static_cast<double>(i) / nu;
            const double fb = arcFractionAt(bottom, u), ft = arcFractionAt(top, u);
            glBegin(GL_LINE_STRIP);
            for (int k = 0; k <= sampleV; ++k) {
                const double v = static_cast<double>(k) / sampleV;
                const Point p = ruledU ? ruledPoint(bottom, top, fb, ft, v)
                                       : coons(bottom, top, left, right, u, v);
                glVertex2d(p[0], p[1]);
            }
            glEnd();
        }
        for (int j = 0; j <= nv; ++j) {
            if (j == 0)  { strip(bottom.pts); continue; }
            if (j == nv) { strip(top.pts);    continue; }
            const double v = static_cast<double>(j) / nv;
            const double fl = arcFractionAt(left, v), fr = arcFractionAt(right, v);
            glBegin(GL_LINE_STRIP);
            for (int k = 0; k <= sampleU; ++k) {
                const double u = static_cast<double>(k) / sampleU;
                const Point p = ruledV ? ruledPoint(left, right, fl, fr, u)
                                       : coons(bottom, top, left, right, u, v);
                glVertex2d(p[0], p[1]);
            }
            glEnd();
        }

        // The grid vertices, on the same interpolants as the lines so the
        // disks land exactly on the crossings. Where a family is ruled its
        // line is the trustworthy one, so its point wins; ticks on the
        // sides come from the sides themselves either way. Shared ticks are
        // drawn once per incident face, which is harmless overdraw.
        if (vertexRadius > 0.0) {
            for (int i = 0; i <= nu; ++i) {
                const double u = static_cast<double>(i) / nu;
                const double fb = arcFractionAt(bottom, u);
                const double ft = arcFractionAt(top, u);
                for (int j = 0; j <= nv; ++j) {
                    const double v = static_cast<double>(j) / nv;
                    Point p{0.0, 0.0};
                    if (i == 0)       p = left.at(v);
                    else if (i == nu) p = right.at(v);
                    else if (j == 0)  p = bottom.at(u);
                    else if (j == nv) p = top.at(u);
                    else if (ruledU)  p = ruledPoint(bottom, top, fb, ft, v);
                    else if (ruledV)  p = ruledPoint(left, right,
                                                     arcFractionAt(left, v),
                                                     arcFractionAt(right, v), u);
                    else              p = coons(bottom, top, left, right, u, v);
                    drawDisk3D(p, vertexRadius, 0.98f, 0.85f, 0.1f);
                }
            }
            // The disks painted over the line color; restore it for the
            // next face's grid.
            viewer::color4f(0.35f, 0.62f, 0.98f, 0.95f);
        }
    }
    viewer::lineWidth(1.0f);
}

}  // namespace

void drawQuantizedBlocks(const BlockQuant &bq) { drawQuantizedTMesh(bq); }

void drawQuantizedLayout(const QuadLayoutQuant &lq, double vertexRadius) {
    drawQuantizedTMesh(lq, vertexRadius);

    // Components that were not four-sided took no part in the quantization,
    // but they are still ground the partition covers -- outline them in the
    // same blue, unsubdivided, so the picture has no hole where one sits and
    // the domain boundary stays unbroken across its share of it.
    viewer::lineWidth(QUANT_LINE_WIDTH);
    viewer::color4f(0.35f, 0.62f, 0.98f, 0.95f);
    for (const auto &outline : lq.skippedOutlines) {
        glBegin(GL_LINE_STRIP);
        for (const Point &p : outline) glVertex2d(p[0], p[1]);
        glEnd();
    }
    viewer::lineWidth(1.0f);
}

namespace {

// Fewer decimals for a coarse step (ticks land on whole numbers) and more
// for a fine one, so a label never reads as more precise than the spacing
// it names, e.g. step=0.05 -> "0.05" not "0.050000" or "0".
std::string formatAxisTick(double v, double step) {
    int decimals = 0;
    if (step > 0.0 && step < 1.0) {
        decimals = static_cast<int>(std::ceil(-std::log10(step)));
        decimals = std::clamp(decimals, 0, 6);
    }
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(decimals) << v;
    return oss.str();
}

} // namespace

void drawAxis(const ViewState &vs) {
    // The switch is here and not at the fifteen call sites, for the same reason
    // the theme's is inside color3f: one place to keep in step.
    if (!axisVisible()) return;

    double worldW = 1.0, worldH = 1.0;
    computeWorldBox(vs, worldW, worldH);

    const double left   = vs.cx - 0.5 * worldW;
    const double right  = vs.cx + 0.5 * worldW;
    const double bottom = vs.cy - 0.5 * worldH;
    const double top    = vs.cy + 0.5 * worldH;

    // Round the tick spacing to 1/2/5 * 10^n so it reads as a scale rather
    // than an arbitrary fraction of the view, with roughly 8 ticks across
    // the wider extent.
    const double ext = std::max(worldW, worldH);
    const double rawStep = ext / 8.0;
    const double mag = std::pow(10.0, std::floor(std::log10(std::max(rawStep, 1e-12))));
    const double norm = rawStep / mag;
    const double niceNorm = (norm < 1.5) ? 1.0 : (norm < 3.0) ? 2.0 : (norm < 7.0) ? 5.0 : 10.0;
    const double step = niceNorm * mag;
    const double tick = 0.01 * ext;

    viewer::lineWidth(1.5f);

    viewer::color3f(0.85f, 0.3f, 0.3f);
    glBegin(GL_LINES);
    glVertex2d(left, 0.0);
    glVertex2d(right, 0.0);
    for (double x = std::ceil(left / step) * step; x <= right; x += step) {
        if (std::fabs(x) < 1e-9) continue; // the axes already cross at the origin
        glVertex2d(x, -tick);
        glVertex2d(x, tick);
    }
    glEnd();

    viewer::color3f(0.3f, 0.8f, 0.35f);
    glBegin(GL_LINES);
    glVertex2d(0.0, bottom);
    glVertex2d(0.0, top);
    for (double y = std::ceil(bottom / step) * step; y <= top; y += step) {
        if (std::fabs(y) < 1e-9) continue;
        glVertex2d(-tick, y);
        glVertex2d(tick, y);
    }
    glEnd();

    viewer::lineWidth(1.0f);

    // Numeric labels in screen space. drawTextOverlay sets up its own
    // screen-space projection (0..fbw, 0..fbh, top-left origin) rather than
    // reusing the world ortho above, so each tick's world position is mapped
    // through the same left/right/bottom/top box by hand; vs.fbw/fbh is the
    // physical size of whatever viewport is currently bound (the full window,
    // or a split-screen half via applyHalfOrtho), so this lines up in both.
    const int fbw = std::max(1, vs.fbw);
    const int fbh = std::max(1, vs.fbh);
    auto worldToScreen = [&](double wx, double wy, float &sx, float &sy) {
        sx = static_cast<float>((wx - left) / (right - left) * fbw);
        sy = static_cast<float>((1.0 - (wy - bottom) / (top - bottom)) * fbh);
    };

    for (double x = std::ceil(left / step) * step; x <= right; x += step) {
        if (std::fabs(x) < 1e-9) continue;
        float sx, sy;
        worldToScreen(x, 0.0, sx, sy);
        drawTextOverlay(fbw, fbh, formatAxisTick(x, step).c_str(), sx + 4.0f, sy + 4.0f,
                        0.95f, 0.55f, 0.55f);
    }
    for (double y = std::ceil(bottom / step) * step; y <= top; y += step) {
        if (std::fabs(y) < 1e-9) continue;
        float sx, sy;
        worldToScreen(0.0, y, sx, sy);
        drawTextOverlay(fbw, fbh, formatAxisTick(y, step).c_str(), sx + 4.0f, sy - 14.0f,
                        0.55f, 0.9f, 0.6f);
    }
    {
        float sx, sy;
        worldToScreen(0.0, 0.0, sx, sy);
        drawTextOverlay(fbw, fbh, "0", sx + 4.0f, sy + 4.0f, 0.75f, 0.75f, 0.75f);
    }
}

namespace {

// Diverging ramp: two hues with a neutral gray midpoint, stepped for a dark
// surface. A quasi-eigenfunction is signed and oscillates about zero, so the
// job is polarity, not magnitude — zero must land exactly on the neutral, and
// the two arms must be opposite hues (never a rainbow, never a hue at the
// midpoint). Blue and red separate under all three CVD types.
void divergingColor(double t, float &r, float &g, float &b) {
    static const float kNeg[3] = {0.224f, 0.529f, 0.898f}; // #3987e5 blue
    static const float kMid[3] = {0.220f, 0.220f, 0.208f}; // #383835 neutral
    static const float kPos[3] = {0.902f, 0.404f, 0.404f}; // #e66767 red

    if (t < -1.0) t = -1.0;
    if (t > 1.0) t = 1.0;

    const float *end = (t < 0.0) ? kNeg : kPos;
    const float a = static_cast<float>(std::fabs(t));
    r = kMid[0] + a * (end[0] - kMid[0]);
    g = kMid[1] + a * (end[1] - kMid[1]);
    b = kMid[2] + a * (end[2] - kMid[2]);
}

} // namespace

void drawScalarField(const Mesh &m, const Eigen::VectorXd &f, double vmax) {
    if (!(vmax > 0.0)) vmax = 1.0;
    const Eigen::Index n = f.size();

    glBegin(GL_TRIANGLES);
    for (const auto &tri : m.triangles) {
        for (int k = 0; k < 3; ++k) {
            const int v = tri[k];
            if (v < 0 || v >= n) continue;
            float r, g, b;
            divergingColor(f[v] / vmax, r, g, b);
            viewer::color3f(r, g, b);
            glVertex2d(m.vertices[v][0], m.vertices[v][1]);
        }
    }
    glEnd();
}

void drawScalarFieldLegend(int fbw, int fbh, double vmin, double vmax, const char *title) {
    const float barW = 220.0f;
    const float barH = 18.0f;
    const float x0 = 20.0f;
    const float y0 = static_cast<float>(fbh) - 70.0f;

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, fbw, fbh, 0, -1, 1); // top-left origin

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    // The ramp itself, drawn as a strip so the neutral midpoint is visible.
    const int steps = 64;
    glBegin(GL_QUAD_STRIP);
    for (int i = 0; i <= steps; ++i) {
        const double t = -1.0 + 2.0 * static_cast<double>(i) / steps;
        float r, g, b;
        divergingColor(t, r, g, b);
        viewer::color3f(r, g, b);
        const float x = x0 + barW * static_cast<float>(i) / steps;
        glVertex2f(x, y0);
        glVertex2f(x, y0 + barH);
    }
    glEnd();

    // Border
    viewer::color3f(0.55f, 0.55f, 0.55f);
    viewer::lineWidth(1.0f);
    glBegin(GL_LINE_LOOP);
    glVertex2f(x0, y0);
    glVertex2f(x0 + barW, y0);
    glVertex2f(x0 + barW, y0 + barH);
    glVertex2f(x0, y0 + barH);
    glEnd();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    // Labels: the extremes and the zero level, in ink rather than ramp colors.
    std::ostringstream lo, hi;
    lo << std::scientific << std::setprecision(2) << vmin;
    hi << std::scientific << std::setprecision(2) << vmax;

    drawTextOverlay(fbw, fbh, title, x0, y0 - 22.0f, 0.8f, 0.8f, 0.8f);
    drawTextOverlay(fbw, fbh, lo.str().c_str(), x0, y0 + barH + 6.0f, 0.8f, 0.8f, 0.8f);
    drawTextOverlay(fbw, fbh, "0", x0 + 0.5f * barW - 6.0f, y0 + barH + 6.0f, 0.8f, 0.8f, 0.8f);
    drawTextOverlay(fbw, fbh, hi.str().c_str(), x0 + barW - 60.0f, y0 + barH + 6.0f, 0.8f, 0.8f, 0.8f);
}

void drawBoundaryEdges(const Mesh &m) {
    viewer::color3f(0.7f, 0.7f, 0.7f);
    viewer::lineWidth(2.0f);
    glBegin(GL_LINES);
    for (int beIdx : m.boundaryEdges) {
        const Point &a = m.vertices[m.edges[beIdx][0]];
        const Point &b = m.vertices[m.edges[beIdx][1]];
        glVertex2d(a[0], a[1]);
        glVertex2d(b[0], b[1]);
    }
    glEnd();
    viewer::lineWidth(1.0f);
}

void computeUVGParamBounds(const UVGParam &uvp, double &cx, double &cy, double &baseW, double &baseH) {
    const Eigen::VectorXd &u = uvp.getU();
    const Eigen::VectorXd &v = uvp.getV();

    if (u.size() == 0) {
        cx = 0.0; cy = 0.0; baseW = 1.0; baseH = 1.0;
        return;
    }

    double minU = u.minCoeff(), maxU = u.maxCoeff();
    double minV = v.minCoeff(), maxV = v.maxCoeff();

    double du = maxU - minU;
    double dv = maxV - minV;
    double ext = std::max(du, dv);
    if (ext <= 0) ext = 1.0;
    double pad = 0.1 * ext;

    cx = 0.5 * (minU + maxU);
    cy = 0.5 * (minV + maxV);
    baseW = du + 2.0 * pad;
    baseH = dv + 2.0 * pad;
    if (baseW <= 0.0) baseW = 1.0;
    if (baseH <= 0.0) baseH = 1.0;
}

void drawUVGParam(const UVGParam &uvp) {
    const Eigen::VectorXd &u = uvp.getU();
    const Eigen::VectorXd &v = uvp.getV();
    const Mesh &cutMesh = uvp.getCutMesh().getCutMesh();

    if (u.size() == 0) return;

    int nV = static_cast<int>(u.size());

    // Draw triangle edges in UV space
    viewer::color3f(0.3f, 0.8f, 0.9f);
    viewer::lineWidth(1.5f);
    glBegin(GL_LINES);
    for (const auto &tri : cutMesh.triangles) {
        for (int e = 0; e < 3; ++e) {
            int a = tri[e];
            int b = tri[(e + 1) % 3];
            if (a < 0 || a >= nV || b < 0 || b >= nV) continue;
            glVertex2d(u(a), v(a));
            glVertex2d(u(b), v(b));
        }
    }
    glEnd();
    viewer::lineWidth(1.0f);
}

void drawFlippedUVTriangles(const UVGParam &uvp) {
    const Eigen::VectorXd &u = uvp.getU();
    const Eigen::VectorXd &v = uvp.getV();
    const Mesh &cutMesh = uvp.getCutMesh().getCutMesh();

    if (u.size() == 0) return;
    int nV = static_cast<int>(u.size());

    viewer::color4f(0.9f, 0.1f, 0.1f, 0.45f);
    glBegin(GL_TRIANGLES);
    for (const auto &tri : cutMesh.triangles) {
        int i = tri[0], j = tri[1], k = tri[2];
        if (i < 0 || i >= nV || j < 0 || j >= nV || k < 0 || k >= nV) continue;
        double signedArea2 = (u(j) - u(i)) * (v(k) - v(i)) - (u(k) - u(i)) * (v(j) - v(i));
        if (signedArea2 < 0.0) {
            glVertex2d(u(i), v(i));
            glVertex2d(u(j), v(j));
            glVertex2d(u(k), v(k));
        }
    }
    glEnd();
}

void drawSingularitiesOnUVG(const UVGParam &uvp,
                             const std::vector<std::pair<int, double>> &singularVertices,
                             double radius) {
    const Eigen::VectorXd &u = uvp.getU();
    const Eigen::VectorXd &v = uvp.getV();
    const auto &origToCut = uvp.getCutMesh().getOriginalToCutVertices();
    int nV = static_cast<int>(u.size());

    for (const auto &[origVtx, crossIndex] : singularVertices) {
        if (origVtx < 0 || origVtx >= static_cast<int>(origToCut.size())) continue;
        const auto &cutVerts = origToCut[origVtx];
        if (cutVerts.empty()) continue;

        int cutVid = cutVerts[0];
        if (cutVid < 0 || cutVid >= nV) continue;

        Point center{u(cutVid), v(cutVid)};
        if (crossIndex > 0)
            drawDisk3D(center, radius, 0.2f, 0.2f, 0.95f);
        else
            drawDisk3D(center, radius, 0.95f, 0.2f, 0.2f);
    }
}

// ── UMBER polysquare ─────────────────────────────────────────────────────────

void computePolysquareBounds(const Polysquare &ps, double &cx, double &cy,
                             double &baseW, double &baseH) {
    const auto &uv = ps.getUV();
    if (uv.empty()) {
        cx = 0.0; cy = 0.0; baseW = 1.0; baseH = 1.0;
        return;
    }

    double minU = uv[0][0], maxU = uv[0][0];
    double minV = uv[0][1], maxV = uv[0][1];
    for (const auto &p : uv) {
        minU = std::min(minU, p[0]); maxU = std::max(maxU, p[0]);
        minV = std::min(minV, p[1]); maxV = std::max(maxV, p[1]);
    }

    const double du = maxU - minU;
    const double dv = maxV - minV;
    double ext = std::max(du, dv);
    if (ext <= 0) ext = 1.0;
    const double pad = 0.1 * ext;

    cx = 0.5 * (minU + maxU);
    cy = 0.5 * (minV + maxV);
    baseW = du + 2.0 * pad;
    baseH = dv + 2.0 * pad;
    if (baseW <= 0.0) baseW = 1.0;
    if (baseH <= 0.0) baseH = 1.0;
}

void drawPolysquare(const Polysquare &ps) {
    const auto &uv = ps.getUV();
    const Mesh &cm = ps.getCutMesh();
    if (uv.empty()) return;
    const int nV = static_cast<int>(uv.size());

    viewer::color3f(0.3f, 0.8f, 0.9f);
    viewer::lineWidth(1.0f);
    glBegin(GL_LINES);
    for (const auto &tri : cm.triangles) {
        for (int e = 0; e < 3; ++e) {
            const int a = tri[e], b = tri[(e + 1) % 3];
            if (a < 0 || a >= nV || b < 0 || b >= nV) continue;
            glVertex2d(uv[a][0], uv[a][1]);
            glVertex2d(uv[b][0], uv[b][1]);
        }
    }
    glEnd();
}

void drawFlippedPolysquareTriangles(const Polysquare &ps) {
    const auto &uv = ps.getUV();
    const Mesh &cm = ps.getCutMesh();
    if (uv.empty()) return;
    const int nV = static_cast<int>(uv.size());

    viewer::color4f(0.9f, 0.1f, 0.1f, 0.45f);
    glBegin(GL_TRIANGLES);
    for (const auto &tri : cm.triangles) {
        const int i = tri[0], j = tri[1], k = tri[2];
        if (i < 0 || i >= nV || j < 0 || j >= nV || k < 0 || k >= nV) continue;
        const double area2 = (uv[j][0] - uv[i][0]) * (uv[k][1] - uv[i][1]) -
                             (uv[k][0] - uv[i][0]) * (uv[j][1] - uv[i][1]);
        if (area2 < 0.0) {
            glVertex2d(uv[i][0], uv[i][1]);
            glVertex2d(uv[j][0], uv[j][1]);
            glVertex2d(uv[k][0], uv[k][1]);
        }
    }
    glEnd();
}

void drawPolysquareStructure(const Polysquare &ps, const HarmonicCut &hc) {
    const auto &uv = ps.getUV();
    const Mesh &cm = ps.getCutMesh();
    if (uv.empty()) return;
    const int nV = static_cast<int>(uv.size());
    const auto &toOrig = hc.getCutVertexToOriginal();
    const auto &cuts = hc.getCutEdges();

    // The banks first, so the boundary of the model draws over them.
    for (int pass = 0; pass < 2; ++pass) {
        if (pass == 0) { viewer::color3f(1.0f, 0.2f, 0.9f); viewer::lineWidth(2.5f); }
        else           { viewer::color3f(1.0f, 0.85f, 0.3f); viewer::lineWidth(3.5f); }

        glBegin(GL_LINES);
        for (int e : cm.boundaryEdges) {
            const int a = cm.edges[e][0], b = cm.edges[e][1];
            if (a < 0 || a >= nV || b < 0 || b >= nV) continue;
            if (a >= static_cast<int>(toOrig.size()) || b >= static_cast<int>(toOrig.size())) continue;

            const bool onCut = cuts.count(MeshEdgeKey(toOrig[a], toOrig[b])) > 0;
            if (onCut != (pass == 0)) continue;

            glVertex2d(uv[a][0], uv[a][1]);
            glVertex2d(uv[b][0], uv[b][1]);
        }
        glEnd();
    }
    viewer::lineWidth(1.0f);
}

void drawPolysquareCorners(const Polysquare &ps, const HarmonicCut &hc,
                           const std::vector<std::pair<int, int>> &corners, double radius) {
    const auto &uv = ps.getUV();
    if (uv.empty()) return;
    const int nV = static_cast<int>(uv.size());
    const auto &origToCut = hc.getOriginalToCutVertices();

    for (const auto &[origVtx, k] : corners) {
        if (origVtx < 0 || origVtx >= static_cast<int>(origToCut.size())) continue;
        for (int cv : origToCut[origVtx]) {
            if (cv < 0 || cv >= nV) continue;
            if (k == 1)       drawDisk3D(uv[cv], radius, 0.2f, 0.2f, 0.95f);
            else if (k == -1) drawDisk3D(uv[cv], radius, 0.95f, 0.2f, 0.2f);
            else              drawDisk3D(uv[cv], radius, 0.95f, 0.85f, 0.1f);
        }
    }
}

// ── UMBER block structure ────────────────────────────────────────────────────

namespace {
// One colour for every edge of the block decomposition, traced or boundary.
inline void blockEdgeColour() { viewer::color3f(1.0f, 0.75f, 0.15f); }
} // namespace

void drawBlockEdges(const MotorcycleGraph &mg, bool parameterDomain, float lineWidth) {
    const auto &segs = mg.getSegments();
    if (segs.empty()) return;

    blockEdgeColour();
    viewer::lineWidth(lineWidth);
    glBegin(GL_LINES);
    for (const auto &s : segs) {
        const Point &a = parameterDomain ? s.ua : s.a;
        const Point &b = parameterDomain ? s.ub : s.b;
        glVertex2d(a[0], a[1]);
        glVertex2d(b[0], b[1]);
    }
    glEnd();
    viewer::lineWidth(1.0f);
}

void drawBlockBoundary(const Polysquare &ps, const HarmonicCut &hc,
                       bool parameterDomain, float lineWidth) {
    const auto &uv = ps.getUV();
    const Mesh &cm = ps.getCutMesh();
    if (uv.empty()) return;
    const int nV = static_cast<int>(uv.size());
    const auto &toOrig = hc.getCutVertexToOriginal();
    const auto &cuts = hc.getCutEdges();

    blockEdgeColour();
    viewer::lineWidth(lineWidth);
    glBegin(GL_LINES);
    for (int e : cm.boundaryEdges) {
        const int a = cm.edges[e][0], b = cm.edges[e][1];
        if (a < 0 || a >= nV || b < 0 || b >= nV) continue;
        if (a >= static_cast<int>(toOrig.size()) || b >= static_cast<int>(toOrig.size())) continue;
        if (cuts.count(MeshEdgeKey(toOrig[a], toOrig[b]))) continue; // a seam, not an edge

        const Point &pa = parameterDomain ? uv[a] : cm.vertices[a];
        const Point &pb = parameterDomain ? uv[b] : cm.vertices[b];
        glVertex2d(pa[0], pa[1]);
        glVertex2d(pb[0], pb[1]);
    }
    glEnd();
    viewer::lineWidth(1.0f);
}

void drawBlockNodes(const MotorcycleGraph &mg, bool parameterDomain, double radius) {
    for (const auto &n : mg.getNodes()) {
        const Point &p = parameterDomain ? n.uv : n.xy;
        switch (n.kind) {
            case MotorcycleGraph::Node::Corner:
                drawDisk3D(p, radius, 0.95f, 0.85f, 0.1f);
                break;
            case MotorcycleGraph::Node::BoundaryEnd:
                drawDisk3D(p, radius, 0.1f, 0.9f, 0.3f);
                break;
            case MotorcycleGraph::Node::Crossing:
                drawDisk3D(p, radius, 0.1f, 0.85f, 0.95f);
                break;
        }
    }
}

void drawQuadLayoutArcs(const QuadLayout &layout, float lineWidth, float r, float g, float b) {
    viewer::color3f(r, g, b);
    viewer::lineWidth(lineWidth);
    for (const auto &arc : layout.getArcs()) {
        if (arc.pts.size() < 2) continue;
        glBegin(GL_LINE_STRIP);
        for (const Point &p : arc.pts) glVertex2d(p[0], p[1]);
        glEnd();
    }
    viewer::lineWidth(1.0f);
}

void drawQuadLayoutNodes(const QuadLayout &layout, double radius) {
    for (const auto &n : layout.getNodes()) {
        switch (n.kind) {
            case QuadLayout::NodeKind::Singularity:
                drawDisk3D(n.pos, radius, 0.95f, 0.85f, 0.1f);
                break;
            case QuadLayout::NodeKind::BoundaryCorner:
                drawDisk3D(n.pos, radius, 0.95f, 0.55f, 0.1f);
                break;
            case QuadLayout::NodeKind::BoundaryHit:
                drawDisk3D(n.pos, radius, 0.1f, 0.9f, 0.3f);
                break;
            case QuadLayout::NodeKind::Crossing:
                drawDisk3D(n.pos, radius, 0.1f, 0.85f, 0.95f);
                break;
            case QuadLayout::NodeKind::Heteroclinic:
                drawDisk3D(n.pos, radius, 0.85f, 0.1f, 0.85f);
                break;
            case QuadLayout::NodeKind::TJunction:
                drawDisk3D(n.pos, radius, 0.2f, 0.4f, 0.95f);
                break;
            case QuadLayout::NodeKind::Dangling:
                drawDisk3D(n.pos, radius, 0.9f, 0.1f, 0.1f);
                break;
        }
    }
}

// ============================================================================
// MERIDIAN -- cones, the cutting graph, and the flat cone metric
// ============================================================================

namespace {

// Blue at +1, red at -1, cyan at -2, yellow beyond. The first two follow the
// convention the rest of the viewer already uses for a cross field's index; the
// cyan is Fig. 11's colour for the valence-six cones it merges pairs of
// valence-five ones into.
void coneColor(int index, float &r, float &g, float &b) {
    if (index == 1)       { r = 0.25f; g = 0.45f; b = 0.95f; }
    else if (index == -1) { r = 0.95f; g = 0.30f; b = 0.25f; }
    else if (index == -2) { r = 0.20f; g = 0.85f; b = 0.90f; }
    else                  { r = 0.95f; g = 0.85f; b = 0.15f; }
}

void drawPolylineOnMesh(const Mesh &m, const std::vector<int> &path) {
    if (path.size() < 2) return;
    glBegin(GL_LINE_STRIP);
    for (const int v : path) {
        if (v < 0 || v >= static_cast<int>(m.vertices.size())) continue;
        glVertex2d(m.vertices[v][0], m.vertices[v][1]);
    }
    glEnd();
}

} // namespace

// ── Stage 0b: the material interfaces ────────────────────────────────────────

namespace {

// One colour per kind of node, in the order Interfaces::NodeKind declares them.
// The kinds are not degrees of the same thing -- a junction is a property of
// three materials meeting, a kink of one interface bending, a balance corner of
// a region's Gauss-Bonnet count -- so they are given colours that do not read
// as a scale.
void interfaceNodeColor(Interfaces::NodeKind kind, float &r, float &g, float &b) {
    switch (kind) {
        case Interfaces::NodeKind::Junction:  r = 0.97f; g = 0.97f; b = 0.97f; break;
        case Interfaces::NodeKind::Landing:   r = 0.25f; g = 0.85f; b = 0.35f; break;
        case Interfaces::NodeKind::Kink:      r = 0.98f; g = 0.60f; b = 0.10f; break;
        case Interfaces::NodeKind::LoopSplit: r = 0.72f; g = 0.35f; b = 0.95f; break;
        case Interfaces::NodeKind::Balance:   r = 0.20f; g = 0.85f; b = 0.90f; break;
        default:                              r = 0.95f; g = 0.15f; b = 0.15f; break;
    }
}

} // namespace

void drawMaterialFill(const Mesh &m, float alpha) {
    if (m.triangleMatId.size() != m.triangles.size()) return;

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBegin(GL_TRIANGLES);
    for (std::size_t t = 0; t < m.triangles.size(); ++t) {
        float r, g, b;
        materialColor(m.triangleMatId[t], r, g, b);
        viewer::color4f(r, g, b, alpha);
        const auto &tri = m.triangles[t];
        for (int k = 0; k < 3; ++k) {
            const Point &p = m.vertices[tri[k]];
            glVertex2d(p[0], p[1]);
        }
    }
    glEnd();
    glDisable(GL_BLEND);
}

void drawInterfaceNetwork(const Interfaces &itf, double nodeRadius, float lineWidth) {
    if (!itf.multiMaterial()) return;
    const Mesh &m = itf.getMesh();

    // The branches, in the white the wireframe already gives an interface edge,
    // so the network and the edges it is made of read as one thing. A closed
    // branch -- an inclusion the loop splitting left whole -- is drawn in the
    // same colour: what makes it different is that it carries no node, which
    // the picture shows by there being none on it.
    viewer::color3f(0.97f, 0.97f, 0.99f);
    viewer::lineWidth(lineWidth);
    for (const Interfaces::Branch &br : itf.branches()) drawPolylineOnMesh(m, br.verts);
    viewer::lineWidth(1.0f);

    if (!(nodeRadius > 0.0)) return;
    for (const Interfaces::Node &n : itf.nodes()) {
        if (n.vertex < 0 || n.vertex >= static_cast<int>(m.vertices.size())) continue;
        const Point &p = m.vertices[n.vertex];
        // The halo is drawn under the disk, so an ill-posed node is a ring
        // rather than a different colour: the kind still has to be readable.
        if (!n.wellPosed) drawDisk3D(p, 1.6 * nodeRadius, 0.35f, 0.10f, 0.10f);
        float r, g, b;
        interfaceNodeColor(n.kind, r, g, b);
        drawDisk3D(p, nodeRadius, r, g, b);
    }
}

void drawInterfaceLegend(int fbw, int fbh, bool separatrixLegendShown) {
    struct Row { Interfaces::NodeKind kind; const char *label; };
    static const Row kRows[] = {
        { Interfaces::NodeKind::Junction,  "junction (3+ branches)" },
        { Interfaces::NodeKind::Landing,   "landing on dS"          },
        { Interfaces::NodeKind::Kink,      "kink"                   },
        { Interfaces::NodeKind::LoopSplit, "loop split"             },
        { Interfaces::NodeKind::Balance,   "balance corner"         },
        { Interfaces::NodeKind::Dangling,  "dangling (tag error)"   },
    };
    const int kRowCount = static_cast<int>(sizeof(kRows) / sizeof(kRows[0]));

    const float x0 = 20.0f;
    const float sw = 16.0f;
    const float lineH = 22.0f;
    // Above the cone legend, and above the separatrix legend when that one is
    // between them.
    float y0 = static_cast<float>(fbh) - (separatrixLegendShown ? 466.0f : 334.0f);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, fbw, fbh, 0, -1, 1); // top-left origin
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glBegin(GL_QUADS);
    for (int i = 0; i < kRowCount; ++i) {
        float r, g, b;
        interfaceNodeColor(kRows[i].kind, r, g, b);
        viewer::color3f(r, g, b);
        const float y = y0 + i * lineH;
        glVertex2f(x0, y);
        glVertex2f(x0 + sw, y);
        glVertex2f(x0 + sw, y + sw);
        glVertex2f(x0, y + sw);
    }
    glEnd();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    for (int i = 0; i < kRowCount; ++i) {
        drawTextOverlay(fbw, fbh, kRows[i].label, x0 + sw + 8.0f, y0 + i * lineH + 2.0f,
                        0.8f, 0.8f, 0.8f);
    }
}

void drawCones(const Mesh &m, const ConeSingularities &cones, double radius) {
    if (!(radius > 0.0)) return;
    for (const auto &c : cones.getCones()) {
        if (c.vertex < 0 || c.vertex >= static_cast<int>(m.vertices.size())) continue;
        const Point &p = m.vertices[c.vertex];
        // An interior cone is the kind that still needs an arc of the cutting
        // graph run out to the boundary, so it is haloed and a boundary cone,
        // which is already there, is not.
        if (!c.onBoundary) drawDisk3D(p, 1.45 * radius, 0.85f, 0.85f, 0.88f);
        float r, g, b;
        coneColor(c.index, r, g, b);
        drawDisk3D(p, radius, r, g, b);
    }
}

void drawConeLegend(int fbw, int fbh) {
    struct Row { int index; const char *label; };
    static const Row kRows[] = {
        {  1, "I=+1  valence 3" },
        { -1, "I=-1  valence 5" },
        { -2, "I=-2  valence 6" },
        { -3, "|I|>2 higher"    },
    };

    const float x0 = 20.0f;
    const float sw = 16.0f;                       // swatch side
    const float lineH = 22.0f;
    float y0 = static_cast<float>(fbh) - 190.0f;  // above the scalar-field legend

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, fbw, fbh, 0, -1, 1); // top-left origin
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glBegin(GL_QUADS);
    for (int i = 0; i < 4; ++i) {
        float r, g, b;
        coneColor(kRows[i].index, r, g, b);
        viewer::color3f(r, g, b);
        const float y = y0 + i * lineH;
        glVertex2f(x0, y);
        glVertex2f(x0 + sw, y);
        glVertex2f(x0 + sw, y + sw);
        glVertex2f(x0, y + sw);
    }
    glEnd();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    for (int i = 0; i < 4; ++i) {
        drawTextOverlay(fbw, fbh, kRows[i].label, x0 + sw + 8.0f, y0 + i * lineH + 2.0f,
                        0.8f, 0.8f, 0.8f);
    }
}

// ── TORSION Stage 3F: the combed frames ──────────────────────────────────────

namespace {

// The palette a_t is cycled through. Eight entries, chosen to be
// distinguishable next to each other rather than to form a ramp: a_t is an
// integer label and not a magnitude, and a ramp would invite reading a
// difference of three as bigger than a difference of one when what matters is
// only whether two neighbouring faces share it.
void branchColor(int a, float &r, float &g, float &b) {
    static const float kPal[8][3] = {
        {0.20f, 0.65f, 0.95f},  // blue
        {0.95f, 0.60f, 0.15f},  // orange
        {0.35f, 0.80f, 0.40f},  // green
        {0.90f, 0.35f, 0.55f},  // pink
        {0.75f, 0.70f, 0.20f},  // olive
        {0.55f, 0.45f, 0.90f},  // violet
        {0.25f, 0.80f, 0.78f},  // teal
        {0.85f, 0.30f, 0.25f},  // red
    };
    const int k = ((a % 8) + 8) % 8;
    r = kPal[k][0]; g = kPal[k][1]; b = kPal[k][2];
}

} // namespace

void drawCombedFrames(const Mesh &m, const FieldFrames &ff, double scale) {
    const std::vector<double> &th = ff.combedAngle();
    const std::vector<int> &a = ff.branch();
    if (th.size() != m.triangles.size() || a.size() != th.size()) return;

    viewer::lineWidth(2.5f);
    for (int t = 0; t < static_cast<int>(m.triangles.size()); ++t) {
        const Triangle &tri = m.triangles[t];
        const Point &p0 = m.vertices[tri[0]];
        const Point &p1 = m.vertices[tri[1]];
        const Point &p2 = m.vertices[tri[2]];
        const Point c = {(p0[0] + p1[0] + p2[0]) / 3.0,
                         (p0[1] + p1[1] + p2[1]) / 3.0};

        float r, g, b;
        branchColor(a[t], r, g, b);

        // X_t full length, Y_t at two thirds and dimmed. The two rows of J*_t
        // are not interchangeable -- the integration fits grad u to X and
        // grad v to Y -- so the picture has to say which is which, and a cross
        // of four equal arms cannot.
        const Point x{std::cos(th[t]), std::sin(th[t])};
        const Point y{-std::sin(th[t]), std::cos(th[t])};
        drawArrow(c, x, scale, r, g, b);
        drawArrow(c, y, scale * 0.66, 0.45f * r + 0.15f, 0.45f * g + 0.15f,
                  0.45f * b + 0.15f);
    }
}

void drawCombedFrameLegend(int fbw, int fbh, const FieldFrames &ff) {
    const std::vector<int> &a = ff.branch();
    int lo = 0, hi = 0;
    for (int v : a) { lo = std::min(lo, v); hi = std::max(hi, v); }
    // More than eight distinct branches and the palette repeats, so the swatch
    // strip would be claiming a distinction it cannot draw; the range is still
    // worth saying, and the picture is still worth reading for where the steps
    // are rather than for which value each patch has.
    const int rows = std::min(hi - lo + 1, 8);

    const float x0 = 20.0f;
    const float sw = 16.0f;
    const float lineH = 22.0f;
    float y0 = static_cast<float>(fbh) - 190.0f - rows * lineH;

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, fbw, fbh, 0, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glBegin(GL_QUADS);
    for (int i = 0; i < rows; ++i) {
        float r, g, b;
        branchColor(lo + i, r, g, b);
        viewer::color3f(r, g, b);
        const float y = y0 + i * lineH;
        glVertex2f(x0, y);
        glVertex2f(x0 + sw, y);
        glVertex2f(x0 + sw, y + sw);
        glVertex2f(x0, y + sw);
    }
    glEnd();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    for (int i = 0; i < rows; ++i) {
        std::string label = "a = " + std::to_string(lo + i);
        if (i == 0) label += "   (branch of the comb)";
        drawTextOverlay(fbw, fbh, label.c_str(), x0 + sw + 8.0f,
                        y0 + i * lineH + 2.0f, 0.8f, 0.8f, 0.8f);
    }
}

void drawCuttingGraph(const ConeCut &cut, float lineWidth) {
    const Mesh &m = cut.getOriginalMesh();

    viewer::lineWidth(lineWidth);

    // Void arcs: one per hole, both ends on the boundary (Wang et al. Sec. 4.1).
    viewer::color3f(1.0f, 0.2f, 0.9f);
    for (const auto &c : cut.getVoidCuts()) drawPolylineOnMesh(m, c.path);

    // Cone arcs: one per interior cone, one end *at* the cone (Sec. 3.2.2).
    viewer::color3f(1.0f, 0.65f, 0.1f);
    for (const auto &c : cut.getConePaths()) drawPolylineOnMesh(m, c.path);

    viewer::lineWidth(1.0f);
}

// ── the flat cone metric ─────────────────────────────────────────────────────

FlatMetric buildFlatMetric(const RicciFlow &flow) {
    FlatMetric fm;
    const Mesh &m = flow.getMesh();
    const auto flat = flow.flowEdgeLengths();
    if (flat.empty()) return fm;

    auto inputLength = [&](int a, int b) {
        const Point &pa = m.vertices[a];
        const Point &pb = m.vertices[b];
        return std::hypot(pa[0] - pb[0], pa[1] - pb[1]);
    };

    std::unordered_set<MeshEdgeKey, MeshEdgeKeyHash> inputEdges;
    inputEdges.reserve(m.edges.size() * 2);
    for (const auto &e : m.edges) inputEdges.emplace(e[0], e[1]);

    fm.edges.reserve(flat.size());
    double sum = 0.0;
    for (const auto &[key, len] : flat) {
        const double l0 = inputLength(key.a, key.b);
        if (!(len > 0.0) || !(l0 > 0.0)) continue;
        FlatMetric::Edge e;
        e.a = key.a;
        e.b = key.b;
        e.t = std::log(len / l0);
        e.newDiagonal = (inputEdges.find(key) == inputEdges.end());
        sum += e.t;
        fm.edges.push_back(e);
    }
    if (fm.edges.empty()) return fm;

    // Remove the mean: u -> u + c rescales the whole metric and is not part of
    // the answer, so what is drawn is the stretch relative to the average one.
    const double mean = sum / static_cast<double>(fm.edges.size());
    double lo = 0.0, hi = 0.0;
    for (auto &e : fm.edges) {
        e.t -= mean;
        lo = std::min(lo, e.t);
        hi = std::max(hi, e.t);
        fm.absMax = std::max(fm.absMax, std::fabs(e.t));
    }
    fm.absMax = std::max(fm.absMax, 1e-12);
    fm.minRatio = std::exp(lo);
    fm.maxRatio = std::exp(hi);

    // The input edges the flipping took out. Their flat length is undefined
    // without unfolding, but where they *were* is worth seeing next to the
    // diagonals that replaced them.
    for (const auto &e : m.edges) {
        if (flat.find(MeshEdgeKey(e[0], e[1])) == flat.end()) {
            fm.replaced.push_back({e[0], e[1]});
        }
    }
    return fm;
}

void drawFlatMetric(const Mesh &m, const FlatMetric &fm, float lineWidth) {
    if (fm.edges.empty()) return;
    const int nV = static_cast<int>(m.vertices.size());
    auto valid = [&](int a, int b) { return a >= 0 && b >= 0 && a < nV && b < nV; };

    // The removed input edges go down first, dim, so the diagonals that took
    // their place read as an overlay rather than as more of the metric.
    if (!fm.replaced.empty()) {
        viewer::lineWidth(std::max(1.0f, lineWidth * 0.6f));
        viewer::color4f(0.40f, 0.45f, 0.55f, 0.75f);
        glBegin(GL_LINES);
        for (const auto &e : fm.replaced) {
            if (!valid(e[0], e[1])) continue;
            glVertex2d(m.vertices[e[0]][0], m.vertices[e[0]][1]);
            glVertex2d(m.vertices[e[1]][0], m.vertices[e[1]][1]);
        }
        glEnd();
    }

    viewer::lineWidth(lineWidth);
    glBegin(GL_LINES);
    for (const auto &e : fm.edges) {
        if (e.newDiagonal || !valid(e.a, e.b)) continue;
        float r, g, b;
        divergingColor(e.t / fm.absMax, r, g, b);
        viewer::color3f(r, g, b);
        glVertex2d(m.vertices[e.a][0], m.vertices[e.a][1]);
        glVertex2d(m.vertices[e.b][0], m.vertices[e.b][1]);
    }
    glEnd();

    // The diagonals the weighted-Delaunay flipping introduced, on top and in a
    // colour the ramp never produces.
    viewer::lineWidth(lineWidth + 1.5f);
    viewer::color3f(0.15f, 0.95f, 0.35f);
    glBegin(GL_LINES);
    for (const auto &e : fm.edges) {
        if (!e.newDiagonal || !valid(e.a, e.b)) continue;
        glVertex2d(m.vertices[e.a][0], m.vertices[e.a][1]);
        glVertex2d(m.vertices[e.b][0], m.vertices[e.b][1]);
    }
    glEnd();
    viewer::lineWidth(1.0f);
}

// ── the cone angles, unfolded ────────────────────────────────────────────────

std::vector<ConeFan> buildConeFans(const RicciFlow &flow, const ConeSingularities &cones) {
    std::vector<ConeFan> fans;
    const Mesh &m = flow.getMesh();
    const auto &faces = flow.getFlowFaces();
    const auto flat = flow.flowEdgeLengths();
    if (faces.empty() || flat.empty() || cones.getCones().empty()) return fans;

    auto lengthOf = [&](int a, int b) -> double {
        auto it = flat.find(MeshEdgeKey(a, b));
        return (it == flat.end()) ? -1.0 : it->second;
    };

    std::vector<std::vector<int>> incident(m.vertices.size());
    for (int f = 0; f < static_cast<int>(faces.size()); ++f) {
        for (int k = 0; k < 3; ++k) {
            const int v = faces[f][k];
            if (v >= 0 && v < static_cast<int>(incident.size())) incident[v].push_back(f);
        }
    }

    for (const auto &cone : cones.getCones()) {
        const int v = cone.vertex;
        if (v < 0 || v >= static_cast<int>(incident.size())) continue;
        if (incident[v].empty()) continue;

        // In a face (a, b, c) written counter-clockwise, the counter-clockwise
        // sweep at a runs from the edge a->b towards a->c. Chaining those pairs
        // walks the star in order without needing any adjacency structure.
        std::unordered_map<int, std::pair<int, double>> nextOf; // nbr -> (next nbr, angle)
        std::unordered_set<int> isSuccessor;
        bool ok = true;

        for (const int f : incident[v]) {
            int i0 = -1;
            for (int k = 0; k < 3; ++k) if (faces[f][k] == v) i0 = k;
            if (i0 < 0) { ok = false; break; }
            const int a = faces[f][(i0 + 1) % 3];
            const int b = faces[f][(i0 + 2) % 3];

            const double la = lengthOf(v, a), lb = lengthOf(v, b), lab = lengthOf(a, b);
            if (!(la > 0.0) || !(lb > 0.0) || !(lab > 0.0)) { ok = false; break; }
            double cosA = (la * la + lb * lb - lab * lab) / (2.0 * la * lb);
            cosA = std::max(-1.0, std::min(1.0, cosA));

            if (!nextOf.emplace(a, std::make_pair(b, std::acos(cosA))).second) { ok = false; break; }
            isSuccessor.insert(b);
        }
        if (!ok || nextOf.empty()) continue;

        // A boundary fan starts at the neighbour nothing leads to; an interior
        // one is a cycle, so any start will do.
        int start = -1;
        for (const auto &[from, to] : nextOf) {
            if (!isSuccessor.count(from)) { start = from; break; }
        }
        const bool closed = (start < 0);
        if (closed) start = nextOf.begin()->first;

        ConeFan fan;
        fan.vertex = v;
        fan.index = cone.index;
        fan.onBoundary = cone.onBoundary;
        fan.closed = closed;

        double theta = 0.0;
        double maxR = 0.0;
        std::vector<std::array<double, 2>> polar; // (radius, angle)

        int cur = start;
        for (size_t step = 0; step <= nextOf.size(); ++step) {
            const double r = lengthOf(v, cur);
            if (!(r > 0.0)) { ok = false; break; }
            polar.push_back({{r, theta}});
            maxR = std::max(maxR, r);

            auto it = nextOf.find(cur);
            if (it == nextOf.end()) break;          // open fan, walked to the far edge
            theta += it->second.second;
            cur = it->second.first;
            if (cur == start) {
                // Closed fan: record the first neighbour a second time, at the
                // angle the walk arrived back at it. The gap between the two is
                // the cone angle's departure from 2pi, which is the picture.
                const double r2 = lengthOf(v, cur);
                if (r2 > 0.0) polar.push_back({{r2, theta}});
                break;
            }
        }
        // Two points is a legitimate fan, not a degenerate one: a convex
        // boundary corner sharp enough to be filled by a single triangle has
        // exactly two neighbours, and its one wedge is the whole cone angle.
        if (!ok || polar.size() < 2 || !(maxR > 0.0)) continue;

        fan.angleSum = theta;
        fan.ring.reserve(polar.size());
        for (const auto &p : polar) {
            const double r = p[0] / maxR;
            fan.ring.push_back({{r * std::cos(p[1]), r * std::sin(p[1])}});
        }
        fans.push_back(std::move(fan));
    }

    // Lay the gallery out in a grid, interior cones first (they are the ones a
    // cone angle other than 2pi was asked for in the first place).
    std::stable_sort(fans.begin(), fans.end(), [](const ConeFan &a, const ConeFan &b) {
        if (a.onBoundary != b.onBoundary) return !a.onBoundary;
        return a.vertex < b.vertex;
    });

    const int n = static_cast<int>(fans.size());
    const int cols = std::max(1, static_cast<int>(std::ceil(std::sqrt(static_cast<double>(n)))));
    const double pitch = 2.8;
    for (int i = 0; i < n; ++i) {
        fans[i].center = {{(i % cols) * pitch, -(i / cols) * pitch}};
    }
    return fans;
}

void computeConeFanBounds(const std::vector<ConeFan> &fans,
                          double &cx, double &cy, double &baseW, double &baseH) {
    cx = cy = 0.0;
    baseW = baseH = 1.0;
    if (fans.empty()) return;

    double minx = 1e300, maxx = -1e300, miny = 1e300, maxy = -1e300;
    for (const auto &f : fans) {
        minx = std::min(minx, f.center[0] - 1.4);
        maxx = std::max(maxx, f.center[0] + 1.4);
        miny = std::min(miny, f.center[1] - 1.4);
        maxy = std::max(maxy, f.center[1] + 1.4);
    }
    cx = 0.5 * (minx + maxx);
    cy = 0.5 * (miny + maxy);
    baseW = std::max(1e-6, (maxx - minx) * 1.05);
    baseH = std::max(1e-6, (maxy - miny) * 1.05);
}

void drawConeFans(const std::vector<ConeFan> &fans) {
    for (const auto &f : fans) {
        if (f.ring.size() < 2) continue;
        const double cx = f.center[0], cy = f.center[1];
        float r, g, b;
        coneColor(f.index, r, g, b);

        // The unit circle the fan is measured against: a full turn of the
        // one-ring, so the gap or the overlap against it is the cone angle.
        viewer::lineWidth(1.0f);
        viewer::color4f(0.45f, 0.45f, 0.50f, 0.8f);
        glBegin(GL_LINE_LOOP);
        for (int i = 0; i < 64; ++i) {
            const double a = 2.0 * M_PI * i / 64.0;
            glVertex2d(cx + std::cos(a), cy + std::sin(a));
        }
        glEnd();

        // The triangles, translucent, so that where a fan of more than 2pi laps
        // itself the overlap shows up as a brighter wedge.
        viewer::color4f(r, g, b, 0.28f);
        glBegin(GL_TRIANGLES);
        for (size_t k = 0; k + 1 < f.ring.size(); ++k) {
            glVertex2d(cx, cy);
            glVertex2d(cx + f.ring[k][0], cy + f.ring[k][1]);
            glVertex2d(cx + f.ring[k + 1][0], cy + f.ring[k + 1][1]);
        }
        glEnd();

        // Their edges: the spokes to each neighbour, and the chords between.
        viewer::lineWidth(1.0f);
        viewer::color4f(r, g, b, 0.75f);
        glBegin(GL_LINES);
        for (size_t k = 0; k < f.ring.size(); ++k) {
            glVertex2d(cx, cy);
            glVertex2d(cx + f.ring[k][0], cy + f.ring[k][1]);
        }
        glEnd();
        glBegin(GL_LINE_STRIP);
        for (const auto &p : f.ring) glVertex2d(cx + p[0], cy + p[1]);
        glEnd();

        // The two spokes that carry the answer. On an interior fan they are the
        // same mesh edge reached from either side of the star, so the angle
        // between them is exactly the cone angle less 2pi; on a boundary fan
        // they are the two boundary edges, and the angle between them is the
        // cone angle itself.
        viewer::lineWidth(3.0f);
        glBegin(GL_LINES);
        viewer::color3f(0.95f, 0.95f, 0.95f);
        glVertex2d(cx, cy);
        glVertex2d(cx + f.ring.front()[0], cy + f.ring.front()[1]);
        viewer::color3f(0.98f, 0.85f, 0.15f);
        glVertex2d(cx, cy);
        glVertex2d(cx + f.ring.back()[0], cy + f.ring.back()[1]);
        glEnd();
        viewer::lineWidth(1.0f);

        drawDisk3D(Point{cx, cy}, 0.11, r, g, b);
    }
}

// ============================================================================
// MERIDIAN -- the layout Psi on Omega
// ============================================================================

void computeLayoutBounds(const std::vector<Point> &uv,
                         double &cx, double &cy, double &baseW, double &baseH) {
    cx = cy = 0.0;
    baseW = baseH = 1.0;
    if (uv.empty()) return;

    double minx = 1e300, maxx = -1e300, miny = 1e300, maxy = -1e300;
    for (const Point &p : uv) {
        minx = std::min(minx, p[0]);  maxx = std::max(maxx, p[0]);
        miny = std::min(miny, p[1]);  maxy = std::max(maxy, p[1]);
    }
    const double du = maxx - minx, dv = maxy - miny;
    double ext = std::max(du, dv);
    if (!(ext > 0.0)) ext = 1.0;
    const double pad = 0.08 * ext;

    cx = 0.5 * (minx + maxx);
    cy = 0.5 * (miny + maxy);
    baseW = std::max(1e-9, du + 2.0 * pad);
    baseH = std::max(1e-9, dv + 2.0 * pad);
}

void drawLayoutUV(const Immersion &imm, const SubdomainLabels *labels,
                  const std::vector<Point> &uv, double coneRadius) {
    const Mesh &cm = imm.getCutMesh();
    const int nV = static_cast<int>(uv.size());
    if (nV == 0) return;

    auto ok = [&](int v) { return v >= 0 && v < nV; };

    // Q1 first, underneath everything: a triangle the map turned over. It is
    // filled rather than outlined because a fold is usually a handful of
    // triangles in a crease and an outline of one is lost among the wireframe.
    viewer::color4f(0.9f, 0.1f, 0.1f, 0.45f);
    glBegin(GL_TRIANGLES);
    for (const auto &tri : cm.triangles) {
        const int i = tri[0], j = tri[1], k = tri[2];
        if (!ok(i) || !ok(j) || !ok(k)) continue;
        const double area2 = (uv[j][0] - uv[i][0]) * (uv[k][1] - uv[i][1]) -
                             (uv[k][0] - uv[i][0]) * (uv[j][1] - uv[i][1]);
        if (area2 < 0.0) {
            glVertex2d(uv[i][0], uv[i][1]);
            glVertex2d(uv[j][0], uv[j][1]);
            glVertex2d(uv[k][0], uv[k][1]);
        }
    }
    glEnd();

    // The triangulation, in the same cyan the other parameter domains use.
    viewer::color3f(0.3f, 0.8f, 0.9f);
    viewer::lineWidth(1.0f);
    glBegin(GL_LINES);
    for (const auto &tri : cm.triangles) {
        for (int e = 0; e < 3; ++e) {
            const int a = tri[e], b = tri[(e + 1) % 3];
            if (!ok(a) || !ok(b)) continue;
            glVertex2d(uv[a][0], uv[a][1]);
            glVertex2d(uv[b][0], uv[b][1]);
        }
    }
    glEnd();

    // The two banks of each arc of G. Drawn apart because Q4 is a statement
    // about the pair: they are the same curve up to R_k, and an arc whose banks
    // are not congruent is where E4 still has work.
    viewer::lineWidth(2.5f);
    for (const Immersion::Arc &arc : imm.getArcs()) {
        for (int side = 0; side < 2; ++side) {
            const std::vector<int> &chain = side == 0 ? arc.plusChain : arc.minusChain;
            if (chain.size() < 2) continue;
            if (side == 0) viewer::color4f(0.98f, 0.70f, 0.15f, 0.95f);
            else           viewer::color4f(0.90f, 0.25f, 0.90f, 0.95f);
            glBegin(GL_LINE_STRIP);
            for (const int v : chain) {
                if (!ok(v)) continue;
                glVertex2d(uv[v][0], uv[v][1]);
            }
            glEnd();
        }
    }

    // dS, coloured by the label Stage 5 gave it. Q3 is read straight off this:
    // every blue run should be vertical and every green one horizontal.
    viewer::lineWidth(3.0f);
    glBegin(GL_LINES);
    if (labels) {
        for (const auto &be : labels->boundaryEdges()) {
            if (!ok(be.a) || !ok(be.b)) continue;
            if (be.label == SubdomainLabels::Align::U)      viewer::color3f(0.35f, 0.55f, 0.95f);
            else if (be.label == SubdomainLabels::Align::V) viewer::color3f(0.20f, 0.85f, 0.40f);
            else                                            viewer::color3f(0.70f, 0.70f, 0.72f);
            glVertex2d(uv[be.a][0], uv[be.a][1]);
            glVertex2d(uv[be.b][0], uv[be.b][1]);
        }
    } else {
        viewer::color3f(0.70f, 0.70f, 0.72f);
        for (const int be : cm.boundaryEdges) {
            const int a = cm.edges[be][0], b = cm.edges[be][1];
            if (!ok(a) || !ok(b)) continue;
            glVertex2d(uv[a][0], uv[a][1]);
            glVertex2d(uv[b][0], uv[b][1]);
        }
    }
    glEnd();
    viewer::lineWidth(1.0f);

    // The feature chains -- on a multi-material model, the interfaces -- in the
    // image, which is where E3 and E6 are read off. Q3's rule applies to them
    // in the same form as to dS: a chain labelled u should come out vertical
    // and one labelled v horizontal, and a chain that is neither is where the
    // E3 residual is. At a node of the network the picture says something E3
    // cannot: the angle between two chains meeting there should be a whole
    // number of right angles, and E6 is the term that puts it there.
    //
    // Drawn in the label colours with a white core, because that is what the
    // model panel draws an interface in: the two halves have to be the same
    // curve to the eye.
    if (labels) {
        for (int pass = 0; pass < 2; ++pass) {
            viewer::lineWidth(pass == 0 ? 4.5f : 1.5f);
            for (const auto &fc : labels->featureChains()) {
                if (fc.verts.size() < 2) continue;
                if (pass == 0) {
                    if (fc.label == SubdomainLabels::Align::U)      viewer::color3f(0.35f, 0.55f, 0.95f);
                    else if (fc.label == SubdomainLabels::Align::V) viewer::color3f(0.20f, 0.85f, 0.40f);
                    else                                            viewer::color3f(0.70f, 0.70f, 0.72f);
                } else {
                    viewer::color3f(0.97f, 0.97f, 0.99f);
                }
                glBegin(GL_LINE_STRIP);
                for (const int v : fc.verts) {
                    if (!ok(v)) continue;
                    glVertex2d(uv[v][0], uv[v][1]);
                }
                glEnd();
            }
        }
        viewer::lineWidth(1.0f);
    }

    // The cones, at every child the cut left them with, in the index colours
    // the model panel uses -- so a cone can be found in both halves at once.
    if (coneRadius > 0.0) {
        const std::vector<std::vector<int>> &children = imm.getConeChildren();
        const std::vector<int> &idx = imm.getConeIndices();
        for (size_t c = 0; c < children.size(); ++c) {
            float r, g, b;
            coneColor(c < idx.size() ? idx[c] : 0, r, g, b);
            for (const int v : children[c]) {
                if (!ok(v)) continue;
                drawDisk3D(Point{uv[v][0], uv[v][1]}, coneRadius, r, g, b);
            }
        }
    }
}

// ============================================================================
// MERIDIAN -- the separatrices of Psi
// ============================================================================

namespace {

// How a curve ended, as the picture wants to say it. Q5's two allowed endings
// are the first two; the third is not an ending at all but the thing the eye
// most needs pointed out -- a curve that passed a cone close enough to have
// been meant to stop at it and went on out through the boundary instead. It
// looks exactly like a healthy boundary-bound separatrix until it is coloured
// differently, and it is the commonest symptom of a missing connectivity
// constraint (Remark 3.1: the patch it leaves behind is a sliver).
enum class SepClass { Cone = 0, Boundary = 1, Graze = 2, Unterminated = 3, Untraceable = 4 };

SepClass classifySeparatrix(const Separatrices::Curve &c, double window) {
    switch (c.end) {
        case Separatrices::End::Cone: return SepClass::Cone;
        case Separatrices::End::Boundary:
            return (c.nearestCone >= 0 && std::isfinite(c.nearestConeGap) &&
                    c.nearestConeGap <= window) ? SepClass::Graze : SepClass::Boundary;
        case Separatrices::End::Capped:
        case Separatrices::End::Cycle: return SepClass::Unterminated;
        default: return SepClass::Untraceable;
    }
}

void separatrixColor(SepClass k, float &r, float &g, float &b) {
    switch (k) {
        case SepClass::Cone:         r = 0.25f; g = 0.90f; b = 0.45f; break;
        case SepClass::Boundary:     r = 0.35f; g = 0.60f; b = 0.98f; break;
        case SepClass::Graze:        r = 0.98f; g = 0.75f; b = 0.15f; break;
        case SepClass::Unterminated: r = 0.95f; g = 0.20f; b = 0.20f; break;
        default:                     r = 0.90f; g = 0.35f; b = 0.90f; break;
    }
}

} // namespace

void drawSeparatrices(const Separatrices &sep, Separatrices::Space space,
                      double endRadius, float lineWidth) {
    const std::vector<Separatrices::Curve> &curves = sep.curves();
    if (curves.empty()) return;

    std::vector<Point> pts;
    std::vector<int> breaks;

    // The near-miss window in absolute image units, so a boundary-bound curve
    // that grazed a cone can be told from one that did not.
    const Separatrices::Report &rep = sep.getReport();
    const double window = sep.getOptions().nearMissWindow * rep.extent;

    // The curves first, then the termini, so a disk is never buried under the
    // line of a curve that happens to pass over it.
    viewer::lineWidth(lineWidth);
    for (const Separatrices::Curve &c : curves) {
        float r, g, b;
        separatrixColor(classifySeparatrix(c, window), r, g, b);
        viewer::color3f(r, g, b);

        pts = sep.polyline(c, space, &breaks);
        if (pts.size() < 2) continue;

        // breaks[k] is the index the run after the k-th break starts at, so the
        // runs are [0, b0), [b0, b1), ..., [blast, end). In the image a break is
        // a seam crossing and the two runs sit on opposite banks of the cut; on
        // the model there are none and the whole curve is one strip.
        size_t start = 0;
        for (size_t k = 0; k <= breaks.size(); ++k) {
            const size_t stop = (k < breaks.size()) ? static_cast<size_t>(breaks[k]) : pts.size();
            if (stop > start + 1) {
                glBegin(GL_LINE_STRIP);
                for (size_t i = start; i < stop; ++i) glVertex2d(pts[i][0], pts[i][1]);
                glEnd();
            }
            start = stop;
        }
    }
    viewer::lineWidth(1.0f);

    if (endRadius <= 0.0) return;
    for (const Separatrices::Curve &c : curves) {
        if (c.end == Separatrices::End::Cone || c.steps.empty()) continue;
        const SepClass k = classifySeparatrix(c, window);
        float r, g, b;
        separatrixColor(k, r, g, b);
        drawDisk3D(sep.point(c.steps.back(), true, space), endRadius, r, g, b);

        // A grazing curve is also marked at the cone it passed, not only where
        // it stopped. That cone is the far end of the connectivity constraint
        // Sec. 3.3 would add, and it is generally nowhere near the end of the
        // curve -- which is the whole reason the miss is hard to see.
        if (k == SepClass::Graze && c.nearestChild >= 0) {
            const Point at = (space == Separatrices::Space::Image)
                                 ? sep.getUV()[c.nearestChild]
                                 : sep.getCutMesh().vertices[c.nearestChild];
            drawDisk3D(at, endRadius * 1.6, r, g, b);
        }
    }
}

void drawLayoutPatches(const Arrangement &arr, const SplineFit *fit,
                       double nodeRadius, float lineWidth, int samples) {
    const std::vector<Arrangement::Node> &nodes = arr.getNodes();
    const std::vector<Arrangement::Arc> &arcs = arr.getArcs();
    const std::vector<Arrangement::Face> &faces = arr.getFaces();
    const std::vector<Arrangement::HalfEdge> &halves = arr.getHalfEdges();
    if (arcs.empty() || halves.empty()) return;

    // Every arc of a finished layout bounds two patches, so walking the faces
    // would draw each side twice -- harmless on screen but not on the nodes,
    // which are drawn as disks and would then be laid one over another. One
    // pass, marking as it goes.
    std::vector<char> drawn(arcs.size(), 0);
    std::vector<char> seen(nodes.size(), 0);

    const int steps = (samples < 2) ? 2 : samples;
    std::vector<Point> poly;

    viewer::color3f(0.42f, 0.74f, 1.0f);
    viewer::lineWidth(lineWidth);
    for (int f : arr.patchFaces()) {
        if (f < 0 || f >= static_cast<int>(faces.size())) continue;
        for (int h : faces[f].half) {
            if (h < 0 || h >= static_cast<int>(halves.size())) continue;
            const int a = halves[h].arc;
            if (a < 0 || a >= static_cast<int>(arcs.size()) || drawn[a]) continue;
            drawn[a] = 1;

            const Arrangement::Arc &arc = arcs[a];
            if (arc.from >= 0 && arc.from < static_cast<int>(nodes.size())) seen[arc.from] = 1;
            if (arc.to   >= 0 && arc.to   < static_cast<int>(nodes.size())) seen[arc.to] = 1;

            // The fitted spline where Stage 9 produced one, and the polyline it
            // would have been fitted to where it did not.
            poly.clear();
            const SplineFit::Curve *c =
                (fit && a < static_cast<int>(fit->curves().size()) &&
                 fit->curves()[a].spline.size() >= 2)
                    ? &fit->curves()[a]
                    : nullptr;
            if (c) {
                poly.reserve(steps + 1);
                for (int i = 0; i <= steps; ++i)
                    poly.push_back(fit->evaluate(*c, static_cast<double>(i) / steps));
            } else {
                poly = arc.points;
            }
            if (poly.size() < 2) continue;

            glBegin(GL_LINE_STRIP);
            for (const Point &p : poly) glVertex2d(p[0], p[1]);
            glEnd();
        }
    }
    viewer::lineWidth(1.0f);

    // The nodes on top, so a disk is never buried under the side of the
    // neighbouring patch that runs into it.
    if (nodeRadius <= 0.0) return;
    for (size_t n = 0; n < nodes.size(); ++n) {
        if (!seen[n]) continue;
        drawDisk3D(nodes[n].p, nodeRadius, 0.15f, 0.88f, 0.30f);
    }
}

namespace {

// One structured block, as the two counts and the vertex array both QuadMesh
// and DiskTemplate store it in. The two types are the same grid and differ only
// in which stage built it, so the drawing takes this rather than either.
struct GridRef {
    int ns = 0, nt = 0;
    const std::vector<int> *vert = nullptr;
};

// The body of drawQuadMesh, written against the arrays rather than against the
// stage that produced them, so that Stage 10's mesh and the merged mesh of
// Stage 11 are drawn by one routine and not by two that could drift apart.
void drawQuadMeshArrays(const std::vector<Point> &V,
                        const std::vector<std::array<int, 4>> &Q,
                        const std::vector<int> &mat,
                        const std::vector<GridRef> &grids,
                        float lineWidth, float blockLineWidth, bool materialFill) {
    if (V.empty() || Q.empty()) return;

    auto ok = [&](int v) { return v >= 0 && v < static_cast<int>(V.size()); };
    auto signedArea = [&](const std::array<int, 4> &q) {
        double a = 0.0;
        for (int k = 0; k < 4; ++k) {
            const Point &p = V[q[k]];
            const Point &n = V[q[(k + 1) & 3]];
            a += p[0] * n[1] - n[0] * p[1];
        }
        return 0.5 * a;
    };

    // The materials, filled, under everything -- including under the folds, so
    // a folded element in the middle of a region still reads as red.
    if (materialFill && mat.size() == Q.size()) {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glBegin(GL_QUADS);
        for (std::size_t i = 0; i < Q.size(); ++i) {
            const auto &q = Q[i];
            if (!ok(q[0]) || !ok(q[1]) || !ok(q[2]) || !ok(q[3])) continue;
            float r, g, b;
            materialColor(mat[i], r, g, b);
            viewer::color4f(r, g, b, 0.30f);
            for (int k = 0; k < 4; ++k) glVertex2d(V[q[k]][0], V[q[k]][1]);
        }
        glEnd();
        glDisable(GL_BLEND);
    }

    // The folds, filled, over the material tint and under everything else.
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    viewer::color4f(0.92f, 0.22f, 0.22f, 0.55f);
    glBegin(GL_QUADS);
    for (const auto &q : Q) {
        if (!ok(q[0]) || !ok(q[1]) || !ok(q[2]) || !ok(q[3])) continue;
        if (signedArea(q) > 0.0) continue;
        for (int k = 0; k < 4; ++k) glVertex2d(V[q[k]][0], V[q[k]][1]);
    }
    glEnd();
    glDisable(GL_BLEND);

    // Every edge once. An interior edge is shared by two quads and would
    // otherwise be laid down twice, which on a translucent line reads as a
    // darker line and makes the grid look like it has a pattern in it.
    //
    // On a multi-material mesh the edge carries the material with it, exactly
    // as the triangulation's wireframe does: an edge inside one material takes
    // that material's colour, an edge between two takes white. That makes the
    // interface network readable in the quads without needing the fill on.
    const bool haveMat = mat.size() == Q.size();
    bool multiMat = false;
    if (haveMat)
        for (std::size_t i = 1; i < mat.size() && !multiMat; ++i)
            multiMat = mat[i] != mat[0];

    // -1 marks "not yet seen"; -2 marks "seen with two different materials".
    std::vector<std::pair<std::pair<int, int>, int>> edges;
    edges.reserve(Q.size() * 4);
    for (std::size_t i = 0; i < Q.size(); ++i) {
        const auto &q = Q[i];
        for (int k = 0; k < 4; ++k) {
            int a = q[k], b = q[(k + 1) & 3];
            if (!ok(a) || !ok(b) || a == b) continue;
            edges.emplace_back(std::make_pair(std::min(a, b), std::max(a, b)),
                               haveMat ? mat[i] : 0);
        }
    }
    std::sort(edges.begin(), edges.end());
    // Collapse the duplicates, folding the two incident materials together.
    std::size_t out = 0;
    for (std::size_t i = 0; i < edges.size();) {
        std::size_t j = i;
        int m = edges[i].second;
        while (j < edges.size() && edges[j].first == edges[i].first) {
            if (edges[j].second != m) m = -2;
            ++j;
        }
        edges[out++] = { edges[i].first, m };
        i = j;
    }
    edges.resize(out);

    viewer::lineWidth(lineWidth);
    glBegin(GL_LINES);
    if (multiMat) {
        for (const auto &e : edges) {
            float r, g, b;
            if (e.second == -2) {
                r = g = b = 1.0f;  // interface edge
            } else {
                materialColor(e.second, r, g, b);
            }
            viewer::color3f(r, g, b);
            glVertex2d(V[e.first.first][0],   V[e.first.first][1]);
            glVertex2d(V[e.first.second][0],  V[e.first.second][1]);
        }
    } else {
        viewer::color3f(0.78f, 0.80f, 0.84f);
        for (const auto &e : edges) {
            glVertex2d(V[e.first.first][0],   V[e.first.first][1]);
            glVertex2d(V[e.first.second][0],  V[e.first.second][1]);
        }
    }
    glEnd();
    viewer::lineWidth(1.0f);

    if (blockLineWidth <= 0.0f) return;

    // The block walls: the four sides of each structured grid, taken off the
    // block's own vertex array rather than off the arrangement, so a block that
    // was skipped leaves a gap here exactly as it does in the mesh.
    viewer::color3f(0.42f, 0.74f, 1.0f);
    viewer::lineWidth(blockLineWidth);
    for (const GridRef &b : grids) {
        const int ns = b.ns, nt = b.nt;
        if (!b.vert || ns < 1 || nt < 1) continue;
        if (static_cast<int>(b.vert->size()) != (ns + 1) * (nt + 1)) continue;
        auto at = [&](int i, int j) { return (*b.vert)[j * (ns + 1) + i]; };
        auto strip = [&](auto next) {
            glBegin(GL_LINE_STRIP);
            next();
            glEnd();
        };
        strip([&] { for (int i = 0; i <= ns; ++i) if (ok(at(i, 0)))  glVertex2d(V[at(i, 0)][0],  V[at(i, 0)][1]); });
        strip([&] { for (int i = 0; i <= ns; ++i) if (ok(at(i, nt))) glVertex2d(V[at(i, nt)][0], V[at(i, nt)][1]); });
        strip([&] { for (int j = 0; j <= nt; ++j) if (ok(at(0, j)))  glVertex2d(V[at(0, j)][0],  V[at(0, j)][1]); });
        strip([&] { for (int j = 0; j <= nt; ++j) if (ok(at(ns, j))) glVertex2d(V[at(ns, j)][0], V[at(ns, j)][1]); });
    }
    viewer::lineWidth(1.0f);
}

std::vector<GridRef> gridsOf(const QuadMesh &qm) {
    std::vector<GridRef> g;
    g.reserve(qm.blocks().size());
    for (const QuadMesh::Block &b : qm.blocks()) g.push_back({ b.ns, b.nt, &b.vert });
    return g;
}

} // namespace

void drawQuadMesh(const QuadMesh &qm, float lineWidth, float blockLineWidth,
                  bool materialFill) {
    drawQuadMeshArrays(qm.vertices(), qm.quads(), qm.quadMaterials(), gridsOf(qm),
                       lineWidth, blockLineWidth, materialFill);
}

void drawQuadMesh(const DiskTemplate &dt, const QuadMesh *qm, float lineWidth,
                  float blockLineWidth, bool materialFill) {
    std::vector<GridRef> grids = qm ? gridsOf(*qm) : std::vector<GridRef>();
    for (const DiskTemplate::Block &b : dt.blocks()) grids.push_back({ b.ns, b.nt, &b.vert });
    drawQuadMeshArrays(dt.vertices(), dt.quads(), dt.quadMaterials(), grids,
                       lineWidth, blockLineWidth, materialFill);
}

void drawQuadMesh(const mesh::QuadMesh &sm, const QuadMesh *qm, const DiskTemplate *dt,
                  float lineWidth, float blockLineWidth, bool materialFill) {
    std::vector<GridRef> grids = qm ? gridsOf(*qm) : std::vector<GridRef>();
    if (dt)
        for (const DiskTemplate::Block &b : dt->blocks()) grids.push_back({ b.ns, b.nt, &b.vert });
    drawQuadMeshArrays(sm.vertices, sm.quads, sm.quadMatId, grids,
                       lineWidth, blockLineWidth, materialFill);
}

// The fitted circle of every inclusion, sampled. Dashed would say "removed"
// more plainly, but a dash pattern at this radius reads as a coarse polygon, so
// it is drawn whole and dim instead and the hole underneath does the saying.
void drawInclusionCircles(const std::vector<DiskTemplate::Inclusion> &inclusions,
                          float lineWidth, int samples) {
    if (inclusions.empty() || samples < 3) return;
    viewer::lineWidth(lineWidth);
    viewer::color3f(0.55f, 0.62f, 0.72f);
    for (const DiskTemplate::Inclusion &inc : inclusions) {
        if (inc.circle.radius <= 0.0) continue;
        glBegin(GL_LINE_LOOP);
        for (int k = 0; k < samples; ++k) {
            const double t = 2.0 * M_PI * k / samples;
            glVertex2d(inc.circle.center[0] + inc.circle.radius * std::cos(t),
                       inc.circle.center[1] + inc.circle.radius * std::sin(t));
        }
        glEnd();
    }
    viewer::lineWidth(1.0f);
}

void drawSeparatrixLegend(int fbw, int fbh) {
    struct Row { SepClass kind; const char *label; };
    static const Row kRows[] = {
        { SepClass::Cone,         "ends at a cone"          },
        { SepClass::Boundary,     "ends through dS"         },
        { SepClass::Graze,        "grazed a cone, went on"  },
        { SepClass::Unterminated, "ended nowhere"           },
        { SepClass::Untraceable,  "not traceable"           },
    };
    const int kRowCount = static_cast<int>(sizeof(kRows) / sizeof(kRows[0]));

    const float x0 = 20.0f;
    const float sw = 16.0f;
    const float lineH = 22.0f;
    float y0 = static_cast<float>(fbh) - 312.0f;  // above the cone legend

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, fbw, fbh, 0, -1, 1); // top-left origin
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glBegin(GL_QUADS);
    for (int i = 0; i < kRowCount; ++i) {
        float r, g, b;
        separatrixColor(kRows[i].kind, r, g, b);
        viewer::color3f(r, g, b);
        const float y = y0 + i * lineH;
        glVertex2f(x0, y);
        glVertex2f(x0 + sw, y);
        glVertex2f(x0 + sw, y + sw);
        glVertex2f(x0, y + sw);
    }
    glEnd();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    for (int i = 0; i < kRowCount; ++i) {
        drawTextOverlay(fbw, fbh, kRows[i].label, x0 + sw + 8.0f, y0 + i * lineH + 2.0f,
                        0.8f, 0.8f, 0.8f);
    }
}

} // namespace viewer
