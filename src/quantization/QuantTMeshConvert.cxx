#include "quantization/QuantTMeshConvert.hxx"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

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
            // Not a T-mesh cell, but still real ground: walk every dart of
            // it in face order, true geometry, so the outline is exactly
            // the boundary this component actually has -- including its
            // share of the domain boundary, if any.
            std::vector<Point> outline;
            for (int dart : f.darts) {
                const int arc = QuadLayout::arcOfDart(dart);
                std::vector<Point> pts = arcs[arc].pts;
                if (dart & 1) std::reverse(pts.begin(), pts.end());
                const size_t start = outline.empty() ? 0 : 1;
                outline.insert(outline.end(), pts.begin() + start, pts.end());
            }
            out.skippedOutlines.push_back(std::move(outline));
            continue;
        }
        std::array<std::vector<int>, 4> sides;
        std::array<std::vector<char>, 4> reversed;
        for (int s = 0; s < 4; ++s) {
            for (int dart : f.sides[s]) {
                const int arc = QuadLayout::arcOfDart(dart);
                const int e = edgeFor(out, arcs, arc, h);
                if (e >= static_cast<int>(out.edgeGeometry.size())) {
                    out.edgeGeometry.resize(e + 1);
                    out.edgeGeometry[e] = arcs[arc].pts;
                }
                sides[s].push_back(e);
                // Dart parity: 2*arc+0 runs a -> b, the direction pts is
                // stored in; 2*arc+1 is the reverse.
                reversed[s].push_back((dart & 1) ? 1 : 0);
            }
        }
        out.faceOfLayout[fi] = out.tmesh.addFace(sides[0], sides[1],
                                                 sides[2], sides[3]);
        out.sideReversed.push_back(std::move(reversed));
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

namespace {

// Shared by both sideCurve() overloads below: the two quant structs carry
// their per-edge geometry and per-side reversal flags in the same shape,
// just under different names, so only this needs to know the layout.
SideCurve sideCurveImpl(const QuantTMesh &tmesh,
                        const std::vector<std::vector<Point>> &edgeGeometry,
                        const std::vector<std::array<std::vector<char>, 4>> &sideReversed,
                        int face, int side) {
    SideCurve out;
    const std::vector<int> &edges = tmesh.faces[face].sides[side];
    const std::vector<char> &rev = sideReversed[face][side];
    for (size_t i = 0; i < edges.size(); ++i) {
        const int e = edges[i];
        std::vector<Point> geometry = edgeGeometry[e];
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
        const int n = std::max(0, tmesh.edges[e].x);
        for (int k = 1; k <= n; ++k) out.tickAt.push_back(base + len * k / n);
    }
    // A side that *ends* in zero edges has geometry past its last tick, and
    // at(1.0) would stop there rather than at the corner -- the Coons
    // corner terms then subtract the wrong point and every interior grid
    // line of the face inherits the offset. Snap the last tick to the full
    // length: the zero edge is parametrically zero wide, so its geometry
    // belongs to the final cell -- the mirror image of what already happens
    // at the start of a side, where tick 0 sits at the true corner and a
    // leading zero edge falls into the first cell.
    if (out.tickAt.size() > 1) out.tickAt.back() = out.length;
    return out;
}

}  // namespace

SideCurve sideCurve(const BlockQuant &bq, int face, int side) {
    return sideCurveImpl(bq.tmesh, bq.edgeGeometry, bq.sideReversed, face, side);
}

SideCurve sideCurve(const QuadLayoutQuant &lq, int face, int side) {
    return sideCurveImpl(lq.tmesh, lq.edgeGeometry, lq.sideReversed, face, side);
}

// ─── The medial-axis T-mesh and its Prop. 1 requirements ────────────────────

namespace {

double distanceToSegment(const Point &p, const Point &a, const Point &b) {
    double along = 0.0;
    return pointSegment(p, a, b, along);
}

// Corner and side incidences of the T-mesh nodes: how many faces have a node
// as a corner, and how many have it inside one of their sides (a T-junction).
void nodeIncidences(const BlockQuant &bq, std::vector<int> &asCorner,
                    std::vector<int> &insideSide) {
    asCorner.assign(bq.nodes.size(), 0);
    insideSide.assign(bq.nodes.size(), 0);
    for (const QuantTMesh::Face &face : bq.tmesh.faces) {
        std::unordered_map<int, int> corner;
        for (int side = 0; side < 4; ++side) {
            std::unordered_map<int, int> seen;
            for (int e : face.sides[side]) {
                ++seen[bq.edgeNodes[e][0]];
                ++seen[bq.edgeNodes[e][1]];
            }
            for (const auto &kv : seen) {
                if (kv.second == 1) corner[kv.first] = 1;
                else ++insideSide[kv.first];
            }
        }
        for (const auto &kv : corner) ++asCorner[kv.first];
    }
}

}  // namespace

BlockQuant makeQuantTMesh(const MedialAxisTMesh &tm, double h) {
    if (!(h > 0.0)) h = tm.targetSize();
    const double tol = 1e-4 * tm.targetSize();
    BlockQuant bq = makeQuantTMesh(tm.blocks, tol, h);
    if (!bq.ok) return bq;

    std::set<std::vector<int>> seen;
    auto add = [&](std::vector<int> req, int atLeast, int &counter) {
        if (req.empty()) return;
        std::sort(req.begin(), req.end());
        if (!seen.insert(req).second) return;
        bq.requirements.push_back({std::move(req), atLeast});
        ++counter;
    };

    // (1) The integer distance between two boundary components stays
    // positive: every radius chain -- all T-mesh edges along one radius
    // [m, phi^-1(m)] -- keeps a length of at least 1, so no half-region
    // loses its thickness. A radius a template did not keep as an edge (the
    // two ends of a purple cap) has no chain and needs none.
    for (const std::array<Point, 2> &r : tm.radii()) {
        if (normP(r[1] - r[0]) <= tol) continue;
        std::vector<int> chain;
        for (size_t e = 0; e < bq.edgeGeometry.size(); ++e) {
            bool along = true;
            for (const Point &p : bq.edgeGeometry[e]) {
                if (distanceToSegment(p, r[0], r[1]) > tol) {
                    along = false;
                    break;
                }
            }
            if (along) chain.push_back(static_cast<int>(e));
        }
        add(std::move(chain), 1, bq.radiusRequirements);
    }

    // (2) Every boundary feature -- the arc between two boundary nodes of
    // valence 2, a node only one cell has a corner at -- keeps a length of at
    // least 1, and (3) every boundary loop one of at least 3. [R] Prop. 1
    // asks a loop for 1 only, but a loop of one quad edge is a vertex joined
    // to itself and one of two a digon, neither a boundary a quad mesh can
    // have: a hole narrower than an element would close up in Step 10. Three
    // is the least a polygon has.
    std::vector<int> asCorner, insideSide;
    nodeIncidences(bq, asCorner, insideSide);
    std::unordered_map<int, std::vector<int>> boundaryAt;
    for (size_t e = 0; e < bq.tmesh.edges.size(); ++e) {
        if (!bq.tmesh.edges[e].onBoundary()) continue;
        boundaryAt[bq.edgeNodes[e][0]].push_back(static_cast<int>(e));
        boundaryAt[bq.edgeNodes[e][1]].push_back(static_cast<int>(e));
    }
    auto isFeatureNode = [&](int n) { return asCorner[n] == 1 && insideSide[n] == 0; };
    std::vector<char> used(bq.tmesh.edges.size(), 0);
    for (size_t e0 = 0; e0 < bq.tmesh.edges.size(); ++e0) {
        if (used[e0] || !bq.tmesh.edges[e0].onBoundary()) continue;
        // Walk the loop through e0 into its sequence of (edge, node ahead).
        std::vector<std::pair<int, int>> loop;
        int e = static_cast<int>(e0);
        int node = bq.edgeNodes[e][1];
        while (e >= 0 && !used[e]) {
            used[e] = 1;
            loop.push_back({e, node});
            int next = -1;
            for (int f : boundaryAt[node]) {
                if (!used[f]) { next = f; break; }
            }
            if (next < 0) break;
            node = bq.edgeNodes[next][0] == node ? bq.edgeNodes[next][1] : bq.edgeNodes[next][0];
            e = next;
        }
        // Start the walk just after a feature node, so arcs do not wrap.
        int start = -1;
        for (size_t k = 0; k < loop.size(); ++k) {
            if (isFeatureNode(loop[k].second)) { start = static_cast<int>(k); break; }
        }
        std::vector<int> all;
        for (const auto &step : loop) all.push_back(step.first);
        add(std::move(all), 3, bq.loopRequirements);
        if (start < 0) continue;
        std::vector<int> arc;
        for (size_t k = 1; k <= loop.size(); ++k) {
            const auto &step = loop[(start + k) % loop.size()];
            arc.push_back(step.first);
            if (isFeatureNode(step.second)) {
                add(std::move(arc), 1, bq.featureRequirements);
                arc.clear();
            }
        }
    }
    return bq;
}

// ─── Step 10 ────────────────────────────────────────────────────────────────

namespace {

struct UnionFind {
    std::vector<int> parent;
    int add() {
        parent.push_back(static_cast<int>(parent.size()));
        return parent.back();
    }
    int find(int a) {
        while (parent[a] != a) a = parent[a] = parent[parent[a]];
        return a;
    }
    void unite(int a, int b) {
        a = find(a);
        b = find(b);
        if (a != b) parent[std::max(a, b)] = std::min(a, b);
    }
};

Point coonsPoint(const SideCurve &bottom, const SideCurve &top, const SideCurve &left,
                 const SideCurve &right, double u, double v) {
    const Point p = bottom.at(u) * (1.0 - v) + top.at(u) * v +
                    left.at(v) * (1.0 - u) + right.at(v) * u;
    return p - (bottom.at(0.0) * ((1.0 - u) * (1.0 - v)) + bottom.at(1.0) * (u * (1.0 - v)) +
                top.at(0.0) * ((1.0 - u) * v) + top.at(1.0) * (u * v));
}

QuadExtraction extract(const BlockQuant &bq) {
    QuadExtraction out;
    const QuantTMesh &tm = bq.tmesh;
    if (!tm.consistent()) {
        out.error = "the quantization is not consistent";
        return out;
    }

    UnionFind uf;
    std::vector<Point> pos;
    std::vector<char> boundary;
    auto vertex = [&](const Point &p, bool onBoundary) {
        pos.push_back(p);
        boundary.push_back(onBoundary ? 1 : 0);
        return uf.add();
    };
    for (const Point &p : bq.nodes) vertex(p, false);

    // Each edge's vertices in the direction of its geometry, x + 1 of them;
    // a zero edge keeps its two ends, now one vertex.
    std::vector<std::vector<int>> edgeVerts(tm.edges.size());
    for (size_t e = 0; e < tm.edges.size(); ++e) {
        const std::vector<Point> &g = bq.edgeGeometry[e];
        int a = bq.edgeNodes[e][0];
        int b = bq.edgeNodes[e][1];
        if (normP(bq.nodes[a] - g.front()) > normP(bq.nodes[b] - g.front())) std::swap(a, b);
        const bool onBoundary = tm.edges[e].onBoundary();
        if (onBoundary) boundary[a] = boundary[b] = 1;
        const int x = tm.edges[e].x;
        if (x <= 0) {
            uf.unite(a, b);
            edgeVerts[e] = {a, b};
            continue;
        }
        double length = 0.0;
        for (size_t i = 1; i < g.size(); ++i) length += normP(g[i] - g[i - 1]);
        edgeVerts[e].push_back(a);
        for (int k = 1; k < x; ++k) {
            edgeVerts[e].push_back(vertex(pointAt(g, length * k / x), onBoundary));
        }
        edgeVerts[e].push_back(b);
    }

    // The vertices along one side of a face, corner to corner: sideSum + 1.
    auto sideVerts = [&](int f, int side) {
        std::vector<int> verts;
        const std::vector<int> &edges = tm.faces[f].sides[side];
        for (size_t i = 0; i < edges.size(); ++i) {
            std::vector<int> seq = edgeVerts[edges[i]];
            if (i < bq.sideReversed[f][side].size() && bq.sideReversed[f][side][i]) {
                std::reverse(seq.begin(), seq.end());
            }
            if (verts.empty()) verts.push_back(seq.front());
            if (tm.edges[edges[i]].x <= 0) continue;
            verts.insert(verts.end(), seq.begin() + 1, seq.end());
        }
        return verts;
    };

    std::vector<std::array<int, 4>> quads;
    for (size_t fi = 0; fi < tm.faces.size(); ++fi) {
        const int f = static_cast<int>(fi);
        const int a = static_cast<int>(tm.sideSum(f, 0));
        const int b = static_cast<int>(tm.sideSum(f, 1));
        const std::vector<int> s0 = sideVerts(f, 0), s1 = sideVerts(f, 1);
        const std::vector<int> s2 = sideVerts(f, 2), s3 = sideVerts(f, 3);

        if (a <= 0 || b <= 0) {
            // The face collapses: S1 onto reversed S3 and S2 onto reversed
            // S4, vertex by vertex -- everything onto one point if both are
            // zero.
            ++out.collapsedFaces;
            if (a > 0) {
                for (int i = 0; i <= a; ++i) uf.unite(s0[i], s2[a - i]);
            } else if (b > 0) {
                for (int j = 0; j <= b; ++j) uf.unite(s1[j], s3[b - j]);
            } else {
                for (const auto *s : {&s0, &s1, &s2, &s3}) {
                    for (int v : *s) uf.unite(v, s0.front());
                }
            }
            continue;
        }

        const SideCurve bottom = sideCurve(bq, f, 0);
        const SideCurve right = sideCurve(bq, f, 1);
        const SideCurve top = sideCurve(bq, f, 2).reversed();
        const SideCurve left = sideCurve(bq, f, 3).reversed();

        // grid[j][i], i along S1 (corner 0 -> 1), j along S4 reversed
        // (corner 0 -> 3).
        std::vector<std::vector<int>> grid(b + 1, std::vector<int>(a + 1, -1));
        for (int i = 0; i <= a; ++i) {
            grid[0][i] = s0[i];
            grid[b][i] = s2[a - i];
        }
        for (int j = 0; j <= b; ++j) {
            grid[j][a] = s1[j];
            grid[j][0] = s3[b - j];
        }
        for (int j = 1; j < b; ++j) {
            for (int i = 1; i < a; ++i) {
                grid[j][i] = vertex(coonsPoint(bottom, top, left, right,
                                               static_cast<double>(i) / a,
                                               static_cast<double>(j) / b),
                                    false);
            }
        }
        // The corners as placed decide which way the face turns; every quad
        // of it is emitted CCW.
        const Point c0 = bottom.at(0.0), c1 = bottom.at(1.0), c2 = top.at(1.0), c3 = top.at(0.0);
        const double area2 = cross2(c1 - c0, c2 - c0) + cross2(c2 - c0, c3 - c0);
        for (int j = 0; j < b; ++j) {
            for (int i = 0; i < a; ++i) {
                std::array<int, 4> q{grid[j][i], grid[j][i + 1], grid[j + 1][i + 1], grid[j + 1][i]};
                if (area2 < 0.0) std::swap(q[1], q[3]);
                quads.push_back(q);
            }
        }
    }

    // Resolve the welds: one vertex per class, at the mean of its boundary
    // members if it has any.
    const int nv = static_cast<int>(pos.size());
    std::vector<Point> sum(nv, Point{0.0, 0.0}), sumB(nv, Point{0.0, 0.0});
    std::vector<int> count(nv, 0), countB(nv, 0);
    for (int v = 0; v < nv; ++v) {
        const int r = uf.find(v);
        sum[r] = sum[r] + pos[v];
        ++count[r];
        if (boundary[v]) {
            sumB[r] = sumB[r] + pos[v];
            ++countB[r];
        }
    }
    std::vector<int> index(nv, -1);
    for (const auto &q : quads) {
        std::array<int, 4> r;
        for (int k = 0; k < 4; ++k) r[k] = uf.find(q[k]);
        bool degenerate = false;
        for (int k = 0; k < 4 && !degenerate; ++k) {
            for (int l = k + 1; l < 4; ++l) degenerate = degenerate || r[k] == r[l];
        }
        if (degenerate) {
            ++out.droppedQuads;
            continue;
        }
        for (int k = 0; k < 4; ++k) {
            if (index[r[k]] < 0) {
                index[r[k]] = static_cast<int>(out.vertices.size());
                out.vertices.push_back(countB[r[k]] > 0 ? sumB[r[k]] / countB[r[k]]
                                                         : sum[r[k]] / count[r[k]]);
                out.onBoundary.push_back(countB[r[k]] > 0 ? 1 : 0);
            }
            r[k] = index[r[k]];
        }
        out.quads.push_back(r);
    }

    // Doublets: an interior vertex in exactly two quads, which share both of
    // its edges. Collapsed faces leave one now and then (geom003, multimat/
    // geom004), and no placement of its vertices makes both quads convex --
    // their two angles at it add to 360 degrees -- so one of them is always
    // folded at a corner, however long an untangler works. The two quads
    // [v, b, p, d] and [v, d, q, b] are one quad [b, p, d, q] without it,
    // and Sigma (4 - deg) is unchanged: v's +2 goes, b and d gain +1 each.
    for (bool again = true; again;) {
        again = false;
        std::vector<std::vector<int>> at(out.vertices.size());
        for (size_t qi = 0; qi < out.quads.size(); ++qi) {
            for (const int v : out.quads[qi]) at[v].push_back(static_cast<int>(qi));
        }
        std::vector<char> gone(out.quads.size(), 0);
        for (size_t v = 0; v < out.vertices.size(); ++v) {
            if (out.onBoundary[v] || at[v].size() != 2) continue;
            const int q1 = at[v][0], q2 = at[v][1];
            if (gone[q1] || gone[q2]) continue;
            auto from = [&](int qi) {   // the quad rotated to start at v
                std::array<int, 4> r = out.quads[qi];
                while (r[0] != static_cast<int>(v)) std::rotate(r.begin(), r.begin() + 1, r.end());
                return r;
            };
            const std::array<int, 4> a = from(q1), b = from(q2);
            const bool sameRing = (a[1] == b[1] && a[3] == b[3]) || (a[1] == b[3] && a[3] == b[1]);
            if (!sameRing || a[2] == b[2]) continue;
            out.quads[q1] = {a[1], a[2], a[3], b[2]};
            gone[q2] = 1;
            ++out.doublets;
            again = true;
        }
        if (!again) break;
        std::vector<std::array<int, 4>> kept;
        for (size_t qi = 0; qi < out.quads.size(); ++qi) {
            if (!gone[qi]) kept.push_back(out.quads[qi]);
        }
        out.quads = std::move(kept);
    }
    if (out.doublets > 0) {
        // The vertices no quad uses any more.
        std::vector<int> index(out.vertices.size(), -1);
        std::vector<Point> vertices;
        std::vector<char> onBoundary;
        for (auto &q : out.quads) {
            for (int &v : q) {
                if (index[v] < 0) {
                    index[v] = static_cast<int>(vertices.size());
                    vertices.push_back(out.vertices[v]);
                    onBoundary.push_back(out.onBoundary[v]);
                }
                v = index[v];
            }
        }
        out.vertices = std::move(vertices);
        out.onBoundary = std::move(onBoundary);
    }

    // One orientation for the whole mesh, read off its connectivity. Each
    // face's quads were turned by the sign of the face's own area, and a face
    // Step 10 placed folded -- an exact disc's caps, a sliver the welds turned
    // over -- came out wound against its neighbours: a mesh no placement of
    // its vertices can untangle, since its quads disagree about which way is
    // up. Across every edge two quads share, they must run it opposite ways;
    // each connected piece is walked from one quad and turned to agree, and
    // then turned over as a whole if its signed area says it was wound
    // clockwise.
    {
        std::map<std::pair<int, int>, std::vector<int>> byEdge;
        for (size_t qi = 0; qi < out.quads.size(); ++qi) {
            const auto &q = out.quads[qi];
            for (int k = 0; k < 4; ++k)
                byEdge[std::minmax(q[k], q[(k + 1) % 4])].push_back(static_cast<int>(qi));
        }
        auto runs = [&](int qi, int a, int b) {   // true when quad qi has the edge as a -> b
            const auto &q = out.quads[qi];
            for (int k = 0; k < 4; ++k)
                if (q[k] == a && q[(k + 1) % 4] == b) return true;
            return false;
        };
        auto flip = [&](int qi) { std::swap(out.quads[qi][1], out.quads[qi][3]); };
        std::vector<char> seen(out.quads.size(), 0);
        for (size_t start = 0; start < out.quads.size(); ++start) {
            if (seen[start]) continue;
            std::vector<int> piece{static_cast<int>(start)}, stack{static_cast<int>(start)};
            seen[start] = 1;
            int turned = 0;
            while (!stack.empty()) {
                const int qi = stack.back();
                stack.pop_back();
                for (int k = 0; k < 4; ++k) {
                    const int a = out.quads[qi][k], b = out.quads[qi][(k + 1) % 4];
                    const auto &sharing = byEdge[std::minmax(a, b)];
                    if (sharing.size() != 2) continue;   // dS, or non-manifold
                    const int other = sharing[0] == qi ? sharing[1] : sharing[0];
                    if (seen[other]) continue;
                    seen[other] = 1;
                    if (runs(other, a, b)) {
                        flip(other);
                        ++turned;
                    }
                    piece.push_back(other);
                    stack.push_back(other);
                }
            }
            double area = 0.0;
            for (const int qi : piece) {
                const auto &q = out.quads[qi];
                area += cross2(out.vertices[q[2]] - out.vertices[q[0]],
                               out.vertices[q[3]] - out.vertices[q[1]]);
            }
            if (area < 0.0) {
                for (const int qi : piece) flip(qi);
                turned = static_cast<int>(piece.size()) - turned;
            }
            out.reorientedQuads += turned;
        }
    }

    // Topology: every mesh edge in one or two quads, the boundary edges
    // closing into loops, and V - E + F.
    std::map<std::pair<int, int>, int> edgeUse;
    for (const auto &q : out.quads) {
        for (int k = 0; k < 4; ++k) {
            const int u = q[k], w = q[(k + 1) % 4];
            ++edgeUse[{std::min(u, w), std::max(u, w)}];
        }
    }
    std::unordered_map<int, std::vector<int>> boundaryNext;
    for (const auto &kv : edgeUse) {
        if (kv.second > 2) ++out.nonManifoldEdges;
        if (kv.second == 1) {
            boundaryNext[kv.first.first].push_back(kv.first.second);
            boundaryNext[kv.first.second].push_back(kv.first.first);
        }
    }
    for (const auto &kv : boundaryNext) {
        if (kv.second.size() != 2) ++out.nonManifoldVertices;
    }
    std::unordered_set<int> visited;
    for (const auto &kv : boundaryNext) {
        if (visited.count(kv.first)) continue;
        ++out.boundaryLoops;
        std::vector<int> stack{kv.first};
        while (!stack.empty()) {
            const int v = stack.back();
            stack.pop_back();
            if (!visited.insert(v).second) continue;
            for (int w : boundaryNext[v]) stack.push_back(w);
        }
    }
    out.eulerCharacteristic = static_cast<int>(out.vertices.size()) -
                              static_cast<int>(edgeUse.size()) +
                              static_cast<int>(out.quads.size());
    out.ok = out.nonManifoldEdges == 0 && out.nonManifoldVertices == 0;
    if (!out.ok) {
        out.error = std::to_string(out.nonManifoldEdges) + " non-manifold edges, " +
                    std::to_string(out.nonManifoldVertices) + " non-manifold vertices";
    }
    return out;
}

}  // namespace

QuadExtraction extractQuadMesh(const BlockQuant &bq, int subdivisions) {
    if (subdivisions <= 1) return extract(bq);
    // Every length times k is as consistent as the lengths were (Ax = 0 is
    // homogeneous), and it is the same mesh with each quad cut k x k: the
    // ticks of an edge at multiples of k are where its ticks were.
    BlockQuant fine = bq;
    for (auto &e : fine.tmesh.edges) e.x *= subdivisions;
    return extract(fine);
}
