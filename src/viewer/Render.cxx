#include "viewer/Render.hxx"

#include "viewer/GL.hxx"
#include "viewer/Geometry.hxx"
#include "viewer/Interaction.hxx"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

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

    float scale = 2.0f;
    float charH = 8.0f * scale;
    float padding = 8.0f;
    float lineSpacing = charH + 2.0f;
    
    // Calculate console dimensions
    float consoleHeight = lines_.size() * lineSpacing + padding * 2;
    float consoleWidth = fbw - 20.0f;  // nearly full width with margins
    
    // Draw semi-transparent background
    glColor4f(0.0f, 0.0f, 0.0f, 0.6f);
    glBegin(GL_QUADS);
    glVertex2f(10.0f, startY - padding);
    glVertex2f(10.0f + consoleWidth, startY - padding);
    glVertex2f(10.0f + consoleWidth, startY + consoleHeight - padding);
    glVertex2f(10.0f, startY + consoleHeight - padding);
    glEnd();

    // Draw border
    glColor3f(0.3f, 0.6f, 0.3f);
    glLineWidth(1.0f);
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

void drawMesh(const Mesh &m) {
    glColor3f(0.85f, 0.85f, 0.85f);
    glLineWidth(1.25f);
    glBegin(GL_LINES);
    for (const auto &tri : m.triangles) {
        const Point &a = m.vertices[tri[0]];
        const Point &b = m.vertices[tri[1]];
        const Point &c = m.vertices[tri[2]];
        glVertex2d(a[0], a[1]);
        glVertex2d(b[0], b[1]);
        glVertex2d(b[0], b[1]);
        glVertex2d(c[0], c[1]);
        glVertex2d(c[0], c[1]);
        glVertex2d(a[0], a[1]);
    }
    glEnd();
}

void drawMeshOverlay(const Mesh &m, float r, float g, float b, float a, float lineWidth) {
    glColor4f(r, g, b, a);
    glLineWidth(lineWidth);
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
    glLineWidth(1.0f);
}

void drawEdgeSetOnMesh(const Mesh &m,
                       const std::unordered_set<CutMesh::EdgeKey, CutMesh::EdgeKeyHash> &edges,
                       float r, float g, float b,
                       float lineWidth) {
    glColor3f(r, g, b);
    glLineWidth(lineWidth);
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
    glColor3f(r, g, b);
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
    glLineWidth(2.5f);
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
    glLineWidth(2.5f);
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
    glLineWidth(2.5f);
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
    glLineWidth(2.5f);
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
    glLineWidth(2.5f);
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

void drawTriangleCrossField(const Mesh &m, const SIPG &sipg, double scale) {
    glLineWidth(2.5f);
    const float br = 0.45f, bg = 0.05f, bb = 0.55f; // dark purple

    const Eigen::VectorXcd &u_k = sipg.u_k;

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
    glColor3f(baseR * shade_center, baseG * shade_center, baseB * shade_center);
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
        glColor3f(baseR * shade, baseG * shade, baseB * shade);
        glVertex2d(x, y);
    }
    glEnd();

    // Subtle outline to enhance 3D look
    glLineWidth(1.0f);
    glBegin(GL_LINE_LOOP);
    for (int i = 0; i < segments; ++i) {
        double ang = (static_cast<double>(i) / segments) * 2.0 * M_PI;
        double x = center[0] + radius * std::cos(ang);
        double y = center[1] + radius * std::sin(ang);
        glColor3f(baseR * 0.5f, baseG * 0.5f, baseB * 0.5f);
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

    glColor3f(r, g, b);
    
    float scale = 2.0f; // scale factor for readability
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
    glColor4f(0.9f, 0.2f, 0.2f, 0.7f); // red with some transparency
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
    glColor3f(0.3f, 0.8f, 0.9f); // cyan color for UV mesh
    glLineWidth(1.5f);
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

    glColor4f(0.4f, 0.4f, 0.4f, 0.5f);
    glLineWidth(1.0f);
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

void drawMedialAxis(const MedialAxis &ma, double vertexRadius) {
    // Draw medial axis edges (Voronoi dual edges) in orange
    glColor3f(1.0f, 0.6f, 0.1f);
    glLineWidth(2.5f);
    glBegin(GL_LINES);
    for (const auto &edge : ma.medialEdges) {
        const Point &a = ma.medialVertices[edge[0]].coord;
        const Point &b = ma.medialVertices[edge[1]].coord;
        glVertex2d(a[0], a[1]);
        glVertex2d(b[0], b[1]);
    }
    glEnd();
    glLineWidth(1.0f);

    // Draw medial axis vertices (circumcenters) as small cyan disks
    for (const auto &v : ma.medialVertices) {
        drawDisk3D(v.coord, vertexRadius, 0.1f, 0.85f, 0.85f);
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
            glColor3f(r, g, b);
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
        glColor3f(r, g, b);
        const float x = x0 + barW * static_cast<float>(i) / steps;
        glVertex2f(x, y0);
        glVertex2f(x, y0 + barH);
    }
    glEnd();

    // Border
    glColor3f(0.55f, 0.55f, 0.55f);
    glLineWidth(1.0f);
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
    glColor3f(0.7f, 0.7f, 0.7f);
    glLineWidth(2.0f);
    glBegin(GL_LINES);
    for (int beIdx : m.boundaryEdges) {
        const Point &a = m.vertices[m.edges[beIdx][0]];
        const Point &b = m.vertices[m.edges[beIdx][1]];
        glVertex2d(a[0], a[1]);
        glVertex2d(b[0], b[1]);
    }
    glEnd();
    glLineWidth(1.0f);
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
    glColor3f(0.3f, 0.8f, 0.9f);
    glLineWidth(1.5f);
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
    glLineWidth(1.0f);
}

void drawFlippedUVTriangles(const UVGParam &uvp) {
    const Eigen::VectorXd &u = uvp.getU();
    const Eigen::VectorXd &v = uvp.getV();
    const Mesh &cutMesh = uvp.getCutMesh().getCutMesh();

    if (u.size() == 0) return;
    int nV = static_cast<int>(u.size());

    glColor4f(0.9f, 0.1f, 0.1f, 0.45f);
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

    glColor3f(0.3f, 0.8f, 0.9f);
    glLineWidth(1.0f);
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

    glColor4f(0.9f, 0.1f, 0.1f, 0.45f);
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
        if (pass == 0) { glColor3f(1.0f, 0.2f, 0.9f); glLineWidth(2.5f); }
        else           { glColor3f(1.0f, 0.85f, 0.3f); glLineWidth(3.5f); }

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
    glLineWidth(1.0f);
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
inline void blockEdgeColour() { glColor3f(1.0f, 0.75f, 0.15f); }
} // namespace

void drawBlockEdges(const MotorcycleGraph &mg, bool parameterDomain, float lineWidth) {
    const auto &segs = mg.getSegments();
    if (segs.empty()) return;

    blockEdgeColour();
    glLineWidth(lineWidth);
    glBegin(GL_LINES);
    for (const auto &s : segs) {
        const Point &a = parameterDomain ? s.ua : s.a;
        const Point &b = parameterDomain ? s.ub : s.b;
        glVertex2d(a[0], a[1]);
        glVertex2d(b[0], b[1]);
    }
    glEnd();
    glLineWidth(1.0f);
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
    glLineWidth(lineWidth);
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
    glLineWidth(1.0f);
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

void drawQuadLayoutArcs(const QuadLayout &layout, float lineWidth) {
    glColor3f(0.25f, 0.55f, 1.0f);
    glLineWidth(lineWidth);
    for (const auto &arc : layout.getArcs()) {
        if (arc.pts.size() < 2) continue;
        glBegin(GL_LINE_STRIP);
        for (const Point &p : arc.pts) glVertex2d(p[0], p[1]);
        glEnd();
    }
    glLineWidth(1.0f);
}

} // namespace viewer

