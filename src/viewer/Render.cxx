#include "viewer/Render.hxx"

#include "viewer/GL.hxx"
#include "viewer/Geometry.hxx"

#include <algorithm>
#include <cmath>

#include <Eigen/Dense>

namespace viewer {

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

} // namespace viewer
