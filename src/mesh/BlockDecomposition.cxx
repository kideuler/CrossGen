#include "mesh/BlockDecomposition.hxx"

#include "mesh/BoundaryFeatures.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <unordered_map>

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

// ── For the comparison metrics ──────────────────────────────────────────────

// The same relation BlockQuadMesh::assignIntervals() builds its chords from,
// kept private to each file rather than exported for two users.
struct UnionFind {
    std::vector<int> parent;
    explicit UnionFind(size_t n) : parent(n) { std::iota(parent.begin(), parent.end(), 0); }
    int find(int x) {
        while (parent[x] != x) x = parent[x] = parent[parent[x]];
        return x;
    }
    void unite(int a, int b) { parent[find(b)] = find(a); }
};

// Every side has a polyline to measure. BlockQuadMesh meshes nothing else, so
// nothing else is scored.
bool completeBlock(const BlockDecomposition &D, int b) {
    for (int e : D.blocks[b].edges) {
        if (e < 0 || e >= static_cast<int>(D.edges.size())) return false;
        if (D.edges[e].points.size() < 2) return false;
    }
    return true;
}

// Positive for a counter-clockwise outline. The same shoelace the coverage
// uses: each side repeats the corner the last one ended on, which adds nothing.
double signedBlockArea(const BlockDecomposition &D, int b) {
    double twice = 0.0;
    for (int s = 0; s < 4; ++s) {
        const std::vector<Point> side = D.sidePolyline(b, s);
        for (size_t i = 0; i + 1 < side.size(); ++i) twice += cross2(side[i], side[i + 1]);
    }
    return 0.5 * twice;
}

// The unit tangent into a polyline at one of its ends, read from the first
// tenth of it (B1: not from the first segment, which on a traced curve
// measures the tracing step). The chord to the tenth would lean inward on a
// bending side by half the angle the side turns in that tenth -- 2.25 degrees
// on a 45-degree arc, enough to leave the quarter disk's three right-angled
// corners owing 1.05 quarter turns in regularity()'s count where it is exactly 1.
// So the tangent is the linear term of the least-squares quadratic
// p0 + t s + c s^2 through the points of that tenth -- the polyline's own, and
// two resampled at 1/20 and 1/10 so that a side of few points still has two --
// which is exact on a circle to second order and averages tracing noise
// instead of taking it from one point. Measured over the corpus, the quarter
// disk then scores 1.001 and geom006 4.01 (against 1.05 and 4.24 by chord). A
// side that really hooks within its first tenth is reported at its angle at
// the vertex, not at the tenth: UMBER's geom006 has two sides leaving a point
// of the hole 7 degrees apart and 20 degrees apart a tenth later, and scores
// the corner at the first. A fit that points backwards, which only a polyline
// doubling back on itself could give, falls back to the chord.
Point endTangent(const std::vector<Point> &p, bool atBack) {
    const std::vector<double> s = arcLengths(p);
    const double L = s.back();
    if (!(L > 0.0)) return Point{0.0, 0.0};
    const Point p0 = atBack ? p.back() : p.front();
    const Point tenth = atFraction(p, s, atBack ? 0.90 : 0.10);

    // Normal equations of min sum |q - p0 - t u - c u^2|^2 in the fraction
    // u = s / L, which leaves t's direction alone and keeps the sums of order 1.
    double u2 = 0.0, u3 = 0.0, u4 = 0.0;
    Point r1{0.0, 0.0}, r2{0.0, 0.0};
    auto add = [&](double u, const Point &q) {
        const Point d = q - p0;
        u2 += u * u;
        u3 += u * u * u;
        u4 += u * u * u * u;
        r1 = r1 + d * u;
        r2 = r2 + d * (u * u);
    };
    for (size_t i = 0; i < p.size(); ++i) {
        const double u = (atBack ? L - s[i] : s[i]) / L;
        if (u > 0.0 && u < 0.1) add(u, p[i]);
    }
    add(0.05, atFraction(p, s, atBack ? 0.95 : 0.05));
    add(0.10, tenth);

    const Point chord = tenth - p0;
    const double det = u2 * u4 - u3 * u3;
    const Point t = (det > 0.0) ? (r1 * u4 - r2 * u3) / det : chord;
    return normalizeP(dotP(t, chord) > 0.0 ? t : chord);
}

// The sectors of the header comment. Block corner c = 4 b + k is corner k of
// block b; `of[c]` is its sector, or -1 for a corner of a block not scored.
struct Sectors {
    std::vector<double> alpha;   // per corner, in [0, 2 pi)
    std::vector<int> of;         // per corner
    std::vector<double> area;    // per block, unsigned
    std::vector<double> theta;   // per sector: 2 pi for a closed ring, else the sum of its corners
    std::vector<int> count;      // per sector: its block corners
    std::vector<bool> ring;      // per sector: joined all the way round, so an interior vertex
};

Sectors findSectors(const BlockDecomposition &D) {
    const int nB = static_cast<int>(D.blocks.size());
    const int nE = static_cast<int>(D.edges.size());

    // One tangent per end of each macro edge, shared by the two corners that
    // edge separates there, so that the corners of a sector add up to the
    // angle between its two bounding tangents exactly, whatever the estimate.
    std::vector<Point> tFrom(nE, Point{0.0, 0.0}), tTo(nE, Point{0.0, 0.0});
    for (int e = 0; e < nE; ++e) {
        if (D.edges[e].points.size() < 2) continue;
        tFrom[e] = endTangent(D.edges[e].points, false);
        tTo[e] = endTangent(D.edges[e].points, true);
    }

    Sectors S;
    S.alpha.assign(4 * static_cast<size_t>(nB), 0.0);
    S.of.assign(4 * static_cast<size_t>(nB), -1);
    S.area.assign(nB, 0.0);

    // Each corner meets two edge ends: the start of its own side k and the end
    // of side k - 1. An end is keyed 2 e + (0 at `from`, 1 at `to`), and the
    // two corners holding the same key are neighbours round the vertex.
    std::unordered_map<int, std::vector<int>> holders;
    std::vector<bool> scored(4 * static_cast<size_t>(nB), false);
    for (int b = 0; b < nB; ++b) {
        if (!completeBlock(D, b)) continue;
        const BlockDecomposition::Block &B = D.blocks[b];
        const double a = signedBlockArea(D, b);
        S.area[b] = std::fabs(a);
        const double turn = (a < 0.0) ? -1.0 : 1.0;   // the side the block lies on
        for (int k = 0; k < 4; ++k) {
            const int km = (k + 3) % 4;
            // Side k leaves the corner: read forwards it starts at its edge's
            // `from`. Side k - 1 arrives at it: read forwards it ends at `to`.
            const int eo = B.edges[k], endO = B.flip[k] ? 1 : 0;
            const int ei = B.edges[km], endI = B.flip[km] ? 0 : 1;
            const Point to = endO ? tTo[eo] : tFrom[eo];
            const Point ti = endI ? tTo[ei] : tFrom[ei];
            // Swept from the leaving side to the arriving one through the
            // block: counter-clockwise for a counter-clockwise outline.
            double alpha = std::atan2(turn * cross2(to, ti), dotP(to, ti));
            if (alpha < 0.0) alpha += 2.0 * M_PI;
            const int c = 4 * b + k;
            S.alpha[c] = alpha;
            scored[c] = true;
            holders[2 * eo + endO].push_back(c);
            holders[2 * ei + endI].push_back(c);
        }
    }

    // Joined across an edge only where the edge is neither dS nor an interface
    // and both its blocks are of one material -- the last because not every
    // adapter flags an interface edge that borders a gap.
    UnionFind uf(4 * static_cast<size_t>(nB));
    std::vector<int> joins(4 * static_cast<size_t>(nB), 0);
    for (const auto &kv : holders) {
        const std::vector<int> &cs = kv.second;
        if (cs.size() != 2) continue;
        const BlockDecomposition::MacroEdge &me = D.edges[kv.first / 2];
        if (me.boundary || me.interface) continue;
        if (D.blocks[cs[0] / 4].material != D.blocks[cs[1] / 4].material) continue;
        uf.unite(cs[0], cs[1]);
        ++joins[cs[0]];
        ++joins[cs[1]];
    }

    std::unordered_map<int, int> sectorOf;
    for (int c = 0; c < 4 * nB; ++c) {
        if (!scored[c]) continue;
        const int root = uf.find(c);
        auto it = sectorOf.find(root);
        if (it == sectorOf.end()) {
            it = sectorOf.emplace(root, static_cast<int>(S.theta.size())).first;
            S.theta.push_back(0.0);
            S.count.push_back(0);
            S.ring.push_back(true);
        }
        const int s = it->second;
        S.of[c] = s;
        S.theta[s] += S.alpha[c];
        ++S.count[s];
        if (joins[c] < 2) S.ring[s] = false;
    }
    // A closed ring is a full turn by topology. Its corners sum to 2 pi anyway
    // unless a block is folded, and a fold is not a reason to expect a sixth
    // block at the vertex.
    for (size_t s = 0; s < S.theta.size(); ++s)
        if (S.ring[s]) S.theta[s] = 2.0 * M_PI;
    return S;
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

bool BlockDecomposition::covers() const {
    if (blocks.empty()) return false;
    for (int b = 0; b < static_cast<int>(blocks.size()); ++b) {
        if (!completeBlock(*this, b)) return false;
        for (int e : blocks[b].edges)
            if (!edges[e].boundary && (edges[e].blockA < 0 || edges[e].blockB < 0)) return false;
    }
    return true;
}

double BlockDecomposition::regularity(const BoundaryFeatures &model) const {
    if (!covers()) return std::numeric_limits<double>::quiet_NaN();
    // The smooth band is the model's corner rule, so that a boundary point is
    // a corner to the layout exactly when it is one to minimumDefect().
    const double smooth = model.options().cornerAngle * M_PI / 180.0;
    const Sectors S = findSectors(*this);
    double sum = 0.0;
    for (size_t s = 0; s < S.theta.size(); ++s) {
        double theta = S.theta[s];
        if (!S.ring[s] && std::fabs(theta - M_PI) < smooth) theta = M_PI;
        sum += std::fabs(2.0 * theta / M_PI - S.count[s]);
    }
    // A model corner with no macrovertex on it lies inside a side. Half its
    // shorter boundary edge is the tolerance: a macrovertex a whole mesh edge
    // away has put the corner inside a side just the same.
    for (const BoundaryFeatures::Corner &c : model.corners()) {
        bool onVertex = false;
        for (const MacroVertex &v : vertices) {
            if (normP(v.p - c.p) <= 0.5 * c.step) {
                onVertex = true;
                break;
            }
        }
        if (!onVertex) sum += std::fabs(2.0 * c.angle / M_PI - 2.0);
    }
    // W_min is a lower bound on W for any conforming layout; the tangents
    // this class measures sector angles with can leave W a hair under it, and
    // no layout is more regular than the model allows.
    return std::min(1.0, (1.0 + model.minimumDefect()) / (1.0 + sum));
}

double BlockDecomposition::angleQuality() const {
    if (!covers()) return std::numeric_limits<double>::quiet_NaN();
    const Sectors S = findSectors(*this);
    double num = 0.0, den = 0.0;
    for (size_t c = 0; c < S.of.size(); ++c) {
        const int s = S.of[c];
        if (s < 0) continue;
        const double dev = S.alpha[c] - S.theta[s] / S.count[s];
        const double w = S.area[c / 4];
        num += w * dev * dev;
        den += w;
    }
    if (!(den > 0.0)) return 0.0;
    // Against a right angle, the deviation at which a corner has folded flat.
    return std::max(1.0 - std::sqrt(num / den) / M_PI_2, 0.0);
}

double BlockDecomposition::chordQuality() const {
    if (!covers()) return std::numeric_limits<double>::quiet_NaN();
    const int nB = static_cast<int>(blocks.size());
    const int nE = static_cast<int>(edges.size());
    UnionFind uf(static_cast<size_t>(nE));
    std::vector<double> blockArea(nB, 0.0);
    for (int b = 0; b < nB; ++b) {
        if (!completeBlock(*this, b)) continue;
        blockArea[b] = std::fabs(signedBlockArea(*this, b));
        uf.unite(blocks[b].edges[0], blocks[b].edges[2]);
        uf.unite(blocks[b].edges[1], blocks[b].edges[3]);
    }

    // A chord's weight is the area it runs through, one block-direction at a
    // time, so its pair of sides 0/2 and its pair 1/3 each give the block's
    // area to the chord that crosses them.
    std::vector<double> area(nE, 0.0);
    for (int b = 0; b < nB; ++b) {
        if (!(blockArea[b] > 0.0)) continue;
        area[uf.find(blocks[b].edges[0])] += blockArea[b];
        area[uf.find(blocks[b].edges[1])] += blockArea[b];
    }
    std::vector<double> shortest(nE, std::numeric_limits<double>::infinity()), longest(nE, 0.0);
    for (int e = 0; e < nE; ++e) {
        const int r = uf.find(e);
        if (!(area[r] > 0.0)) continue;
        const std::vector<double> s = arcLengths(edges[e].points);
        shortest[r] = std::min(shortest[r], s.back());
        longest[r] = std::max(longest[r], s.back());
    }

    double num = 0.0, den = 0.0;
    for (int r = 0; r < nE; ++r) {
        if (!(area[r] > 0.0)) continue;
        num += area[r] * std::log(longest[r] / shortest[r]);
        den += area[r];
    }
    if (!(den > 0.0)) return 0.0;
    // 1/S of the weighted geometric mean S: no free scale, and 1/2 is a chord
    // whose elements differ in size by a factor of two end to end.
    return std::exp(-num / den);
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
