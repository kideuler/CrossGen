#include "mesh/BlockDecomposition.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>

namespace {

// Cumulative arc length of a polyline, and the point at a fraction of it.
// Equal arc length rather than equal parameter is the whole reason this is
// here: the sides of a block are traced or fitted curves whose points are as
// close together as the tracing happened to put them, and stepping over them
// by index gives a row of mesh nodes that bunches wherever the tracing did.
std::vector<double> arcLengths(const std::vector<Point> &p) {
    std::vector<double> s(p.size(), 0.0);
    for (size_t i = 1; i < p.size(); ++i) s[i] = s[i - 1] + normP(p[i] - p[i - 1]);
    return s;
}

Point atFraction(const std::vector<Point> &p, const std::vector<double> &s, double t) {
    if (p.empty()) return Point{0.0, 0.0};
    if (p.size() == 1) return p[0];
    const double total = s.back();
    if (!(total > 0.0)) return p[0];
    const double target = std::min(std::max(t, 0.0), 1.0) * total;
    // The first point at or past the target; the segment before it is the one
    // the point sits on.
    const size_t hi = static_cast<size_t>(
        std::lower_bound(s.begin(), s.end(), target) - s.begin());
    if (hi == 0) return p.front();
    if (hi >= p.size()) return p.back();
    const double seg = s[hi] - s[hi - 1];
    const double a = (seg > 0.0) ? (target - s[hi - 1]) / seg : 0.0;
    return p[hi - 1] * (1.0 - a) + p[hi] * a;
}

// The winding-number test rather than the crossing-number one: a block outline
// is a closed polyline that may well have a corner exactly on the ray a
// crossing count would cast, and winding gets that case right.
bool insidePolygon(const std::vector<Point> &poly, const Point &q) {
    if (poly.size() < 3) return false;
    int wind = 0;
    for (size_t i = 0; i < poly.size(); ++i) {
        const Point &a = poly[i];
        const Point &b = poly[(i + 1) % poly.size()];
        if (a[1] <= q[1]) {
            if (b[1] > q[1] && cross2(b - a, q - a) > 0.0) ++wind;
        } else if (b[1] <= q[1] && cross2(b - a, q - a) < 0.0) {
            --wind;
        }
    }
    return wind != 0;
}

} // namespace

std::vector<Point> BlockDecomposition::sidePolyline(int block, int side) const {
    if (block < 0 || block >= static_cast<int>(blocks.size())) return {};
    if (side < 0 || side > 3) return {};
    const Block &b = blocks[block];
    const int e = b.edges[side];
    if (e < 0 || e >= static_cast<int>(edges.size())) return {};
    const MacroEdge &me = edges[e];
    if (!b.flip[side]) return me.points;
    return std::vector<Point>(me.points.rbegin(), me.points.rend());
}

std::vector<Point> BlockDecomposition::sampleSide(int block, int side, int n) const {
    const std::vector<Point> p = sidePolyline(block, side);
    if (p.size() < 2) return {};
    if (n < 1) n = 1;
    const std::vector<double> s = arcLengths(p);
    std::vector<Point> out;
    out.reserve(n + 1);
    for (int i = 0; i <= n; ++i) out.push_back(atFraction(p, s, static_cast<double>(i) / n));
    // The ends are the corners themselves and not an interpolation of them,
    // so that a node placed on a corner by one side and by another is one
    // node to the bit rather than to the tolerance.
    out.front() = p.front();
    out.back() = p.back();
    return out;
}

Point BlockDecomposition::pointOnSide(int block, int side, double t) const {
    const std::vector<Point> p = sidePolyline(block, side);
    if (p.empty()) return Point{0.0, 0.0};
    return atFraction(p, arcLengths(p), t);
}

Point BlockDecomposition::coonsPoint(int block, double u, double v) const {
    u = std::min(std::max(u, 0.0), 1.0);
    v = std::min(std::max(v, 0.0), 1.0);
    // Side 0 runs corner 0 -> 1 (the u direction at v = 0), side 1 corner
    // 1 -> 2 (v at u = 1), side 2 corner 2 -> 3 (u backwards at v = 1), side 3
    // corner 3 -> 0 (v backwards at u = 0). So the top and left curves are
    // their sides read the other way.
    const Point b = pointOnSide(block, 0, u);
    const Point r = pointOnSide(block, 1, v);
    const Point t = pointOnSide(block, 2, 1.0 - u);
    const Point l = pointOnSide(block, 3, 1.0 - v);

    if (block < 0 || block >= static_cast<int>(blocks.size())) return b;
    const Point c00 = pointOnSide(block, 0, 0.0);
    const Point c10 = pointOnSide(block, 1, 0.0);
    const Point c11 = pointOnSide(block, 2, 0.0);
    const Point c01 = pointOnSide(block, 3, 0.0);

    const Point ruled = b * (1.0 - v) + t * v + l * (1.0 - u) + r * u;
    const Point bilinear = c00 * ((1.0 - u) * (1.0 - v)) + c10 * (u * (1.0 - v)) +
                           c11 * (u * v) + c01 * ((1.0 - u) * v);
    return ruled - bilinear;
}

bool BlockDecomposition::interiorPoint(int block, Point &out) const {
    if (block < 0 || block >= static_cast<int>(blocks.size())) return false;

    // The outline, at whatever resolution its sides carry, so that "inside"
    // means inside the block as drawn and not inside the quadrilateral of its
    // four corners.
    std::vector<Point> poly;
    for (int s = 0; s < 4; ++s) {
        const std::vector<Point> side = sidePolyline(block, s);
        for (size_t k = 0; k + 1 < side.size(); ++k) poly.push_back(side[k]);
    }
    if (poly.size() < 3) return false;

    // The Coons centre, which is the answer for anything remotely convex, then
    // a coarse sweep of the parameter square for a block that curls.
    const Point mid = coonsPoint(block, 0.5, 0.5);
    if (insidePolygon(poly, mid)) { out = mid; return true; }
    for (int j = 1; j <= 7; ++j) {
        for (int i = 1; i <= 7; ++i) {
            const Point q = coonsPoint(block, i / 8.0, j / 8.0);
            if (insidePolygon(poly, q)) { out = q; return true; }
        }
    }

    // Last resort, and the one that cannot fail on a simple polygon: at the
    // lowest-leftmost vertex -- necessarily convex -- the chord between its
    // two neighbours is interior unless some other vertex lies in the ear, and
    // then the deepest such vertex gives an interior chord instead.
    size_t v = 0;
    for (size_t i = 1; i < poly.size(); ++i) {
        if (poly[i][1] < poly[v][1] ||
            (poly[i][1] == poly[v][1] && poly[i][0] < poly[v][0])) v = i;
    }
    const Point &a = poly[(v + poly.size() - 1) % poly.size()];
    const Point &c = poly[(v + 1) % poly.size()];
    const Point &p = poly[v];
    size_t deepest = poly.size();
    double best = 0.0;
    for (size_t i = 0; i < poly.size(); ++i) {
        if (i == v || i == (v + 1) % poly.size() ||
            i == (v + poly.size() - 1) % poly.size()) continue;
        // Inside the ear (a, p, c), which is a positively oriented triangle
        // because p is the lowest-leftmost vertex of a counter-clockwise loop.
        const double s0 = cross2(p - a, poly[i] - a);
        const double s1 = cross2(c - p, poly[i] - p);
        const double s2 = cross2(a - c, poly[i] - c);
        if (!((s0 >= 0.0 && s1 >= 0.0 && s2 >= 0.0) ||
              (s0 <= 0.0 && s1 <= 0.0 && s2 <= 0.0))) continue;
        const Point d = c - a;
        const double dd = normP(d);
        const double depth = (dd > 0.0) ? std::fabs(cross2(d, poly[i] - a)) / dd : 0.0;
        if (deepest == poly.size() || depth > best) { deepest = i; best = depth; }
    }
    out = (deepest == poly.size()) ? (a + c) * 0.5 : (p + poly[deepest]) * 0.5;
    return true;
}

bool BlockDecomposition::writeEdgesOBJ(const std::string &path) const {
    std::ofstream out(path);
    if (!out) return false;
    out.precision(17);
    out << "# " << source << " block decomposition: " << blocks.size() << " blocks, "
        << edges.size() << " macro edges, " << vertices.size() << " macrovertices\n";
    int base = 1;
    for (const MacroEdge &me : edges) {
        if (me.points.size() < 2) continue;
        for (const Point &p : me.points) out << "v " << p[0] << " " << p[1] << " 0\n";
        out << "l";
        for (size_t i = 0; i < me.points.size(); ++i) out << " " << (base + static_cast<int>(i));
        out << "\n";
        base += static_cast<int>(me.points.size());
    }
    return static_cast<bool>(out);
}

bool BlockDecomposition::writeBlocksOBJ(const std::string &path) const {
    std::ofstream out(path);
    if (!out) return false;
    out.precision(17);
    out << "# " << source << " blocks: " << blocks.size() << "\n";
    int base = 1;
    for (size_t b = 0; b < blocks.size(); ++b) {
        std::vector<Point> loop;
        for (int s = 0; s < 4; ++s) {
            const std::vector<Point> side = sidePolyline(static_cast<int>(b), s);
            for (size_t k = 0; k + 1 < side.size(); ++k) loop.push_back(side[k]);
        }
        if (loop.size() < 3) continue;
        for (const Point &p : loop) out << "v " << p[0] << " " << p[1] << " 0\n";
        out << "l";
        for (size_t i = 0; i < loop.size(); ++i) out << " " << (base + static_cast<int>(i));
        out << " " << base << "\n";
        base += static_cast<int>(loop.size());
    }
    return static_cast<bool>(out);
}
