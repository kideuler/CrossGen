#include "viewer/Render.hxx"

#include "viewer/GL.hxx"
#include "viewer/Geometry.hxx"

#include <algorithm>
#include <cmath>

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

void Console::draw(GLFWwindow *window, float startY) const {
    if (lines_.empty()) return;

    int fbw, fbh;
    glfwGetFramebufferSize(window, &fbw, &fbh);

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
        drawTextOverlay(window, line.c_str(), 18.0f, y, 0.4f, 0.9f, 0.4f);
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
    glLineWidth(1.5f);
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

void drawTextOverlay(GLFWwindow *window, const char *text, float x, float y, float r, float g, float b) {
    int fbw, fbh;
    glfwGetFramebufferSize(window, &fbw, &fbh);

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

void drawUVMesh(const IntegerGridMap &igm, float r, float g, float b, float lineWidth) {
    const Mesh &cut = igm.getMesh();
    const std::vector<Point> &uv = igm.uv();

    if (uv.size() != cut.vertices.size()) return;

    glColor3f(r, g, b);
    glLineWidth(lineWidth);
    glBegin(GL_LINES);
    for (const auto &tri : cut.triangles) {
        const Point &a = uv[tri[0]];
        const Point &b = uv[tri[1]];
        const Point &c = uv[tri[2]];
        glVertex2d(a[0], a[1]);
        glVertex2d(b[0], b[1]);
        glVertex2d(b[0], b[1]);
        glVertex2d(c[0], c[1]);
        glVertex2d(c[0], c[1]);
        glVertex2d(a[0], a[1]);
    }
    glEnd();
}

} // namespace viewer