#include "viewer/Geometry.hxx"

#include <algorithm>
#include <cmath>

namespace viewer {

Bounds computeBounds(const Mesh &m) {
    Bounds b;
    bool first = true;
    for (const auto &p : m.vertices) {
        if (first) {
            b.minx = b.maxx = p[0];
            b.miny = b.maxy = p[1];
            first = false;
            continue;
        }
        b.minx = std::min(b.minx, p[0]);
        b.maxx = std::max(b.maxx, p[0]);
        b.miny = std::min(b.miny, p[1]);
        b.maxy = std::max(b.maxy, p[1]);
    }
    if (first) {
        b.minx = b.miny = 0;
        b.maxx = b.maxy = 1;
    }
    return b;
}

Point triangleCentroid(const Mesh &m, const Triangle &t) {
    const Point &a = m.vertices[t[0]];
    const Point &b = m.vertices[t[1]];
    const Point &c = m.vertices[t[2]];
    return Point{(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0};
}

double averageTriangleEdgeLength(const Mesh &m) {
    if (m.triangles.empty()) return 1.0;

    auto elen = [&](const Point &p, const Point &q) {
        double dx = q[0] - p[0];
        double dy = q[1] - p[1];
        return std::sqrt(dx * dx + dy * dy);
    };

    double sum = 0.0;
    for (const auto &t : m.triangles) {
        const Point &a = m.vertices[t[0]];
        const Point &b = m.vertices[t[1]];
        const Point &c = m.vertices[t[2]];
        double l0 = elen(a, b);
        double l1 = elen(b, c);
        double l2 = elen(c, a);
        sum += (l0 + l1 + l2) / 3.0;
    }

    return sum / static_cast<double>(m.triangles.size());
}

} // namespace viewer
