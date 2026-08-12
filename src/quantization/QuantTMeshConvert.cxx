#include "quantization/QuantTMeshConvert.hxx"

#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace {

// ─── QuadLayout ──────────────────────────────────────────────────────────────

// Lazily give a layout arc a tmesh edge.
int edgeFor(QuadLayoutQuant &out, const std::vector<QuadLayout::Arc> &arcs,
            int arc, double h) {
    int &e = out.edgeOfArc[arc];
    if (e < 0) {
        e = out.tmesh.addEdge(h > 0.0 ? arcs[arc].length / h : 1.0);
    }
    return e;
}

// ─── Geometric welding of blocks ─────────────────────────────────────────────

double sq(double v) { return v * v; }

// Distance from p to segment ab, and the arc-length offset of the closest
// point from a.
double pointSegment(const Point &p, const Point &a, const Point &b,
                    double &along) {
    const Point d = b - a;
    const double len2 = sq(d[0]) + sq(d[1]);
    double t = 0.0;
    if (len2 > 0.0) {
        t = ((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1]) / len2;
        t = std::max(0.0, std::min(1.0, t));
    }
    const Point c = a + d * t;
    along = t * std::sqrt(len2);
    return normP(p - c);
}

double polylineLength(const std::vector<Point> &pts) {
    double len = 0.0;
    for (size_t i = 1; i < pts.size(); ++i) len += normP(pts[i] - pts[i - 1]);
    return len;
}

// The point at arc length `target` along the polyline.
Point pointAt(const std::vector<Point> &pts, double target) {
    double acc = 0.0;
    for (size_t i = 1; i < pts.size(); ++i) {
        const double seg = normP(pts[i] - pts[i - 1]);
        if (acc + seg >= target && seg > 0.0) {
            return pts[i - 1] + (pts[i] - pts[i - 1]) * ((target - acc) / seg);
        }
        acc += seg;
    }
    return pts.back();
}

// Welds coincident points into node ids with a hash grid.
class NodeSet {
public:
    explicit NodeSet(double tol) : tol_(tol) {}

    int findOrAdd(const Point &p, std::vector<Point> &nodes) {
        const int found = find(p, nodes);
        if (found >= 0) return found;
        nodes.push_back(p);
        const int id = static_cast<int>(nodes.size()) - 1;
        grid_.insert({cell(p), id});
        return id;
    }

    int find(const Point &p, const std::vector<Point> &nodes) const {
        const long long cx = coord(p[0]), cy = coord(p[1]);
        for (long long dx = -1; dx <= 1; ++dx) {
            for (long long dy = -1; dy <= 1; ++dy) {
                auto range = grid_.equal_range(key(cx + dx, cy + dy));
                for (auto it = range.first; it != range.second; ++it) {
                    if (normP(p - nodes[it->second]) < tol_) return it->second;
                }
            }
        }
        return -1;
    }

private:
    long long coord(double v) const {
        return static_cast<long long>(std::floor(v / (2.0 * tol_)));
    }
    long long key(long long cx, long long cy) const {
        return cx * 2000003LL + cy;
    }
    long long cell(const Point &p) const { return key(coord(p[0]), coord(p[1])); }

    double tol_;
    std::unordered_multimap<long long, int> grid_;
};

// The side polylines of a block, cyclically between its corners, one per
// corner. Returns false when a corner cannot be located on the outline.
bool cornerSides(const TMeshBlock &b, double tol,
                 std::vector<std::vector<Point>> &sides) {
    const int nc = static_cast<int>(b.corners.size());
    std::vector<int> idx(nc, -1);
    for (int c = 0; c < nc; ++c) {
        double best = tol;
        for (size_t i = 0; i < b.outline.size(); ++i) {
            const double d = normP(b.corners[c] - b.outline[i]);
            if (d < best) {
                best = d;
                idx[c] = static_cast<int>(i);
            }
        }
        if (idx[c] < 0) return false;
    }
    const int n = static_cast<int>(b.outline.size());
    sides.assign(nc, {});
    for (int c = 0; c < nc; ++c) {
        for (int i = idx[c];; i = (i + 1) % n) {
            sides[c].push_back(b.outline[i]);
            if (i == idx[(c + 1) % nc] && sides[c].size() > 1) break;
        }
    }
    return true;
}

bool blockSides(const TMeshBlock &b, double tol,
                std::array<std::vector<Point>, 4> &sides) {
    std::vector<std::vector<Point>> s;
    if (b.corners.size() != 4 || !cornerSides(b, tol, s)) return false;
    for (int c = 0; c < 4; ++c) sides[c] = std::move(s[c]);
    return true;
}

// Turn a three-sided block into a four-sided cell by cutting its longest
// side in half.
//
// The medial decomposition leaves triangles behind: every cap zone is one,
// and so are the odd pieces of the red and purple templates. A T-mesh cell
// must have four sides, and the alternative to splitting -- treating the
// triangle as a quad with one side of length zero -- would force that side
// to quantize to zero and collapse the cell. Halving the longest side is
// what the paper's own templates do to the caps they consume (Fig. 17,
// where the polar cap is split at its midpoint W and joined to the
// flanks), so the two agree on how a cap is carved.
//
// The new corner is a plain node, and the neighbour across that side picks
// it up as a T-junction through the ordinary splitting pass.
TMeshBlock splitTriangle(const TMeshBlock &b, double tol) {
    std::vector<std::vector<Point>> sides;
    if (b.corners.size() != 3 || !cornerSides(b, tol, sides)) return b;

    int longest = 0;
    double best = -1.0;
    for (int c = 0; c < 3; ++c) {
        const double len = polylineLength(sides[c]);
        if (len > best) {
            best = len;
            longest = c;
        }
    }
    if (!(best > 0.0)) return b;
    const Point mid = pointAt(sides[longest], 0.5 * best);

    TMeshBlock out = b;
    // Place the cut in the outline, unless a vertex already sits there.
    int at = -1;
    for (size_t i = 0; i < out.outline.size(); ++i) {
        if (normP(out.outline[i] - mid) < tol) { at = static_cast<int>(i); break; }
    }
    if (at < 0) {
        // The outline runs corner[longest] -> corner[longest + 1]; the cut
        // goes on the segment of that run which straddles it.
        double acc = 0.0;
        const double half = 0.5 * best;
        size_t k = 1;
        for (; k < sides[longest].size(); ++k) {
            const double step = normP(sides[longest][k] - sides[longest][k - 1]);
            if (acc + step >= half) break;
            acc += step;
        }
        const Point &after = sides[longest][std::min(k, sides[longest].size() - 1)];
        for (size_t i = 0; i < out.outline.size(); ++i) {
            if (normP(out.outline[i] - after) < tol) { at = static_cast<int>(i); break; }
        }
        if (at < 0) return b;
        out.outline.insert(out.outline.begin() + at, mid);
    }
    out.corners.insert(out.corners.begin() + longest + 1, mid);
    return out;
}

// Deduplicates edges shared by two blocks: same endpoint nodes and the same
// curve between them.
//
// Two distinct edges can join the same pair of nodes -- the two long sides
// of a thin sliver do exactly that -- so the endpoints alone cannot decide
// it. The curve is compared at three interior stations rather than at the
// midpoint alone, because a single sample is easily fooled: the two sides
// of a sliver have near-coincident midpoints, and welding them collapses a
// healthy block into a non-cell. Two blocks that genuinely share a side
// trace it from the same source geometry, so their samples agree to
// rounding, far inside the tolerance.
class EdgeSet {
public:
    EdgeSet(double tol, double h) : tol_(tol), h_(h) {}

    // `reversed` reports whether the stored geometry runs against the given
    // one, which happens for the second of the two faces sharing the edge.
    int findOrAdd(int a, int b, std::vector<Point> geometry, BlockQuant &out,
                  bool &reversed) {
        const Signature sig = sample(geometry);
        const long long k = key(a, b);
        auto range = made_.equal_range(k);
        for (auto it = range.first; it != range.second; ++it) {
            if (matches(sig, it->second.first)) {
                const int e = it->second.second;
                reversed = normP(out.edgeGeometry[e].front() - geometry.front()) >
                           normP(out.edgeGeometry[e].back() - geometry.front());
                return e;
            }
        }
        reversed = false;
        const double len = polylineLength(geometry);
        const int e = out.tmesh.addEdge(h_ > 0.0 ? len / h_ : 1.0);
        out.edgeNodes.push_back({std::min(a, b), std::max(a, b)});
        out.edgeGeometry.push_back(std::move(geometry));
        made_.insert({k, {sig, e}});
        return e;
    }

private:
    // The curve at a quarter, a half and three quarters of its length.
    struct Signature { std::array<Point, 3> p; };

    static Signature sample(const std::vector<Point> &geometry) {
        const double len = polylineLength(geometry);
        Signature s;
        for (int i = 0; i < 3; ++i) s.p[i] = pointAt(geometry, len * (i + 1) / 4.0);
        return s;
    }

    // Either traversal direction counts as the same curve: the two blocks
    // sharing a side walk it opposite ways.
    bool matches(const Signature &a, const Signature &b) const {
        bool fwd = true, rev = true;
        for (int i = 0; i < 3; ++i) {
            if (normP(a.p[i] - b.p[i]) >= tol_) fwd = false;
            if (normP(a.p[i] - b.p[2 - i]) >= tol_) rev = false;
        }
        return fwd || rev;
    }

    long long key(int a, int b) const {
        return (static_cast<long long>(std::min(a, b)) << 32) |
               static_cast<long long>(std::max(a, b));
    }

    double tol_;
    double h_;
    std::unordered_multimap<long long, std::pair<Signature, int>> made_;
};

}  // namespace

QuadLayoutQuant makeQuantTMesh(const QuadLayout &layout, double h) {
    QuadLayoutQuant out;
    const auto &arcs = layout.getArcs();
    const auto &faces = layout.getFaces();
    out.edgeOfArc.assign(arcs.size(), -1);
    out.faceOfLayout.assign(faces.size(), -1);

    for (size_t fi = 0; fi < faces.size(); ++fi) {
        const QuadLayout::Face &f = faces[fi];
        if (f.sides.size() != 4) {
            ++out.skippedFaces;
            continue;
        }
        std::array<std::vector<int>, 4> sides;
        for (int s = 0; s < 4; ++s) {
            for (int dart : f.sides[s]) {
                sides[s].push_back(
                    edgeFor(out, arcs, QuadLayout::arcOfDart(dart), h));
            }
        }
        out.faceOfLayout[fi] = out.tmesh.addFace(sides[0], sides[1],
                                                 sides[2], sides[3]);
    }

    out.ok = out.tmesh.finalize(&out.error);
    return out;
}

BlockQuant makeQuantTMesh(const std::vector<TMeshBlock> &rawBlocks,
                          double weldTol, double h) {
    BlockQuant out;
    out.faceOfBlock.assign(rawBlocks.size(), -1);
    NodeSet nodeSet(weldTol);

    // Pass 0: give every triangle a fourth side, so it can be a cell like
    // any other rather than a hole in the decomposition. One in, one out,
    // so block indices still line up with the caller's.
    std::vector<TMeshBlock> blocks;
    blocks.reserve(rawBlocks.size());
    for (const TMeshBlock &b : rawBlocks) {
        blocks.push_back(b.corners.size() == 3 ? splitTriangle(b, weldTol) : b);
    }

    // Pass 1: every corner of every four-sided block becomes a node. All
    // nodes must exist before any side is split, since the T-junctions on a
    // side are other blocks' corners.
    for (const TMeshBlock &b : blocks) {
        if (b.corners.size() != 4) continue;
        for (const Point &c : b.corners) nodeSet.findOrAdd(c, out.nodes);
    }

    // Pass 2: split each side at the nodes on it and weld the pieces.
    //
    // The edge tolerance is far tighter than the node one on purpose. It
    // does not express geometric slack: two blocks sharing a side trace it
    // from the same points, so their samples agree to rounding. It only has
    // to separate distinct curves that run between the same two nodes, and
    // those can be very close indeed -- the spokes of neighbouring cap
    // slivers converge on a common centre -- so anything as loose as the
    // node tolerance welds them together and lands the edge in three faces.
    EdgeSet edgeSet(1e-3 * weldTol, h);
    std::vector<int> edgeUses;  // incident faces committed so far, per edge
    for (size_t bi = 0; bi < blocks.size(); ++bi) {
        const TMeshBlock &b = blocks[bi];
        std::array<std::vector<Point>, 4> sidePts;
        if (b.corners.size() != 4) {
            ++out.skippedBlocks;
            ++out.skippedNonQuad;
            continue;
        }
        if (!blockSides(b, weldTol, sidePts)) {
            ++out.skippedBlocks;
            ++out.skippedUnlocatable;
            continue;
        }

        std::array<std::vector<int>, 4> sides;
        std::array<std::vector<char>, 4> reversed;
        for (int s = 0; s < 4; ++s) {
            const std::vector<Point> &pts = sidePts[s];
            const double total = polylineLength(pts);

            // Every node on the side, by its arc-length position. The
            // endpoints land at 0 and total and are kept; duplicates from
            // adjacent cells of the hash are welded by the sort-unique.
            std::vector<std::pair<double, int>> cuts;
            for (int ni = 0; ni < static_cast<int>(out.nodes.size()); ++ni) {
                double acc = 0.0, bestD = weldTol, bestAt = -1.0;
                for (size_t i = 1; i < pts.size(); ++i) {
                    double along;
                    const double d =
                        pointSegment(out.nodes[ni], pts[i - 1], pts[i], along);
                    if (d < bestD) {
                        bestD = d;
                        bestAt = acc + along;
                    }
                    acc += normP(pts[i] - pts[i - 1]);
                }
                if (bestAt >= 0.0) cuts.push_back({bestAt, ni});
            }
            std::sort(cuts.begin(), cuts.end());

            // Walk the cuts and emit one edge per piece.
            std::vector<Point> piece;
            double prevAt = 0.0;
            int prevNode = cuts.empty() ? -1 : cuts.front().second;
            size_t pi = 0;
            double acc = 0.0;
            for (size_t ci = 1; ci < cuts.size(); ++ci) {
                const auto [at, node] = cuts[ci];
                if (at - prevAt < weldTol) continue;  // welded duplicate
                piece.clear();
                piece.push_back(pointAt(pts, prevAt));
                while (pi + 1 < pts.size() && acc + normP(pts[pi + 1] - pts[pi]) <
                                                  at - weldTol) {
                    acc += normP(pts[pi + 1] - pts[pi]);
                    ++pi;
                    piece.push_back(pts[pi]);
                }
                piece.push_back(pointAt(pts, std::min(at, total)));
                bool rev = false;
                sides[s].push_back(
                    edgeSet.findOrAdd(prevNode, node, piece, out, rev));
                reversed[s].push_back(rev ? 1 : 0);
                prevAt = at;
                prevNode = node;
            }
            if (sides[s].empty()) {
                // No interior cut and endpoints only: the side is one edge.
                const int a = nodeSet.find(pts.front(), out.nodes);
                const int c = nodeSet.find(pts.back(), out.nodes);
                bool rev = false;
                sides[s].push_back(edgeSet.findOrAdd(a, c, pts, out, rev));
                reversed[s].push_back(rev ? 1 : 0);
            }
        }
        // A block whose corners very nearly coincide can have two of its
        // sides collapse onto the same welded edge. That is not a T-mesh
        // cell, and admitting it would make the whole system unusable, so
        // the block is dropped the way a non-quad one is.
        bool degenerate = false;
        std::vector<int> used;
        for (int s = 0; s < 4 && !degenerate; ++s) {
            for (int e : sides[s]) {
                if (std::find(used.begin(), used.end(), e) != used.end()) {
                    degenerate = true;
                    break;
                }
                used.push_back(e);
            }
        }
        if (degenerate) {
            ++out.skippedBlocks;
            ++out.skippedDegenerate;
            continue;
        }

        // An edge borders at most two cells. A block that would claim a
        // third share of one overlaps a pair already in place, so it is
        // dropped rather than allowed to make the whole system unusable.
        bool overlapping = false;
        for (int e : used) {
            if (e < static_cast<int>(edgeUses.size()) && edgeUses[e] >= 2) {
                overlapping = true;
                break;
            }
        }
        if (overlapping) {
            ++out.skippedBlocks;
            ++out.skippedOverlapping;
            continue;
        }
        edgeUses.resize(out.tmesh.edges.size(), 0);
        for (int e : used) ++edgeUses[e];

        const int face = out.tmesh.addFace(sides[0], sides[1],
                                           sides[2], sides[3]);
        out.faceOfBlock[bi] = face;
        out.blockOfFace.push_back(static_cast<int>(bi));
        out.sideReversed.push_back(std::move(reversed));
    }

    out.ok = out.tmesh.finalize(&out.error);
    return out;
}

Point SideCurve::atArc(double s) const {
    if (pts.empty()) return Point{0.0, 0.0};
    if (s <= 0.0) return pts.front();
    if (s >= length) return pts.back();
    const size_t i = static_cast<size_t>(
        std::upper_bound(arc.begin(), arc.end(), s) - arc.begin());
    if (i == 0 || i >= pts.size()) return pts.back();
    const double seg = arc[i] - arc[i - 1];
    const double t = seg > 0.0 ? (s - arc[i - 1]) / seg : 0.0;
    return pts[i - 1] + (pts[i] - pts[i - 1]) * t;
}

Point SideCurve::at(double u) const {
    const int n = cells();
    if (n <= 0) return pts.empty() ? Point{0.0, 0.0} : pts.front();
    // Walk the ticks, not the raw arc length: cell i must span exactly the
    // stretch between tick i and tick i + 1, whatever its length.
    const double s = std::min(std::max(u, 0.0), 1.0) * n;
    int i = static_cast<int>(s);
    if (i >= n) i = n - 1;
    const double f = s - i;
    return atArc(tickAt[i] + (tickAt[i + 1] - tickAt[i]) * f);
}

SideCurve SideCurve::reversed() const {
    SideCurve r;
    r.pts.assign(pts.rbegin(), pts.rend());
    r.length = length;
    r.arc.resize(arc.size());
    for (size_t i = 0; i < arc.size(); ++i) r.arc[i] = length - arc[arc.size() - 1 - i];
    r.tickAt.resize(tickAt.size());
    for (size_t i = 0; i < tickAt.size(); ++i)
        r.tickAt[i] = length - tickAt[tickAt.size() - 1 - i];
    return r;
}

SideCurve sideCurve(const BlockQuant &bq, int face, int side) {
    SideCurve out;
    const std::vector<int> &edges = bq.tmesh.faces[face].sides[side];
    const std::vector<char> &rev = bq.sideReversed[face][side];
    for (size_t i = 0; i < edges.size(); ++i) {
        const int e = edges[i];
        std::vector<Point> geometry = bq.edgeGeometry[e];
        if (i < rev.size() && rev[i]) {
            std::reverse(geometry.begin(), geometry.end());
        }
        const double base = out.length;
        if (out.pts.empty()) {
            out.pts.push_back(geometry.front());
            out.arc.push_back(0.0);
            out.tickAt.push_back(0.0);
        }
        for (size_t k = 1; k < geometry.size(); ++k) {
            const double step = normP(geometry[k] - geometry[k - 1]);
            if (step <= 0.0) continue;  // the joint the previous edge left
            out.length += step;
            out.pts.push_back(geometry[k]);
            out.arc.push_back(out.length);
        }
        // A zero edge contributes no tick: the grid has no cell across it.
        const double len = out.length - base;
        const int n = std::max(0, bq.tmesh.edges[e].x);
        for (int k = 1; k <= n; ++k) out.tickAt.push_back(base + len * k / n);
    }
    return out;
}
