#include "medialaxis/MedialAxisMap.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

// origin + t*dir  =  a + s*(b - a). False on (near) parallel lines. The
// tolerance is relative so the test does not depend on the model's scale.
bool lineSegmentIntersect(const Point &origin, const Point &dir,
                          const Point &a, const Point &b,
                          Point &hit, double &s) {
    const Point e = b - a;
    const double denom = cross2(dir, e);
    const double scale = normP(dir) * normP(e);
    if (scale <= 0.0 || std::abs(denom) <= 1e-12 * scale) return false;
    const Point w = a - origin;
    s = cross2(w, dir) / denom;
    hit = a + e * s;
    return true;
}

Point closestOnSegment(const Point &p, const Point &a, const Point &b,
                       double &dist2) {
    const Point e = b - a;
    const double len2 = dotP(e, e);
    const double t = len2 > 0.0 ? std::clamp(dotP(p - a, e) / len2, 0.0, 1.0) : 0.0;
    const Point q = a + e * t;
    dist2 = dotP(p - q, p - q);
    return q;
}

// Slack on the segment parameter: a polar ray through a chain vertex passes
// through a segment endpoint exactly, and must not be lost to roundoff.
const double SEGMENT_PARAM_EPS = 1e-6;

} // namespace

// ─── Construction ───────────────────────────────────────────────────────────

MedialAxisMap::MedialAxisMap(std::shared_ptr<MedialAxis> axisIn)
    : axis(std::move(axisIn)) {
    if (!axis) {
        throw std::runtime_error("MedialAxis pointer cannot be null.");
    }

    const int numVertices = static_cast<int>(axis->mesh->vertices.size());
    fanOfVertex_.assign(numVertices, -1);
    for (int v = 0; v < numVertices; ++v) {
        buildFan(v);
    }
    stats_.fans = static_cast<int>(fans_.size());
    buildSpokes();
}

const MedialFan *MedialAxisMap::fanOf(int boundaryVertex) const {
    if (boundaryVertex < 0 ||
        boundaryVertex >= static_cast<int>(fanOfVertex_.size()))
        return nullptr;
    const int f = fanOfVertex_[boundaryVertex];
    return f < 0 ? nullptr : &fans_[f];
}

// ─── Fans ───────────────────────────────────────────────────────────────────

void MedialAxisMap::buildFan(int v) {
    const MedialAxis &ma = *axis;
    const int prevB = ma.boundaryLinks[v].prevBoundaryVertex;
    const int nextB = ma.boundaryLinks[v].nextBoundaryVertex;
    if (prevB < 0 && nextB < 0) return;   // interior vertex
    if (prevB < 0 || nextB < 0) {
        // Half-linked: a pinched boundary vertex rebuildBoundaryLinks gave up
        // on. No single fan can describe it.
        ++stats_.skippedVertices;
        return;
    }

    // The fan starts at the cell carrying the constrained edge (prev, v),
    // whose dual Voronoi edge crosses the boundary there.
    const auto it = ma.edgeCells.find(MeshEdgeKey(prevB, v));
    if (it == ma.edgeCells.end()) {
        ++stats_.skippedVertices;
        return;
    }
    int cur = -1;
    for (int c : it->second) {
        if (c >= 0 && ma.cells[c].active) cur = c;
    }
    if (cur < 0) {
        ++stats_.skippedVertices;
        return;
    }

    MedialFan fan;
    fan.boundaryVertex = v;
    const Point &pv = ma.mesh->vertices[v];
    const Point &pp = ma.mesh->vertices[prevB];
    const Point &pn = ma.mesh->vertices[nextB];
    fan.midPrev = (pp + pv) * 0.5;
    fan.midNext = (pv + pn) * 0.5;

    // Walk the cells around v. Rings are CCW with the interior on the left of
    // constrained edges, so crossing the ring edge that *starts* at v rotates
    // the fan from the (prev, v) side towards (v, next); the walk ends when
    // that edge is the constrained edge (v, next).
    int guard = static_cast<int>(ma.cells.size()) + 1;
    while (cur >= 0 && guard-- > 0) {
        fan.cells.push_back(cur);
        fan.chain.push_back(ma.cells[cur].medialVertex);
        const std::vector<int> &ring = ma.cells[cur].verts;
        int idx = -1;
        for (int i = 0; i < static_cast<int>(ring.size()); ++i) {
            if (ring[i] == v) { idx = i; break; }
        }
        if (idx < 0) break;   // inconsistent complex; keep what we have
        cur = ma.neighborAcross(fan.cells.back(), idx);
    }
    if (fan.chain.empty() || fan.chain.front() < 0) {
        ++stats_.skippedVertices;
        return;
    }

    // Polar center: intersection of the two boundary-crossing Voronoi edges.
    // Each runs along the perpendicular bisector of its boundary edge, so the
    // midpoint is a numerically safe origin; the outward normal (the interior
    // lies left of prev -> v -> next) picks the branch that leaves the domain.
    const Point ePrev = pv - pp;
    const Point eNext = pn - pv;
    const Point u1 = normalizeP(Point{ePrev[1], -ePrev[0]});
    const Point u2 = normalizeP(Point{eNext[1], -eNext[0]});
    const double denom = cross2(u1, u2);
    if (std::abs(denom) > 1e-9) {
        const double t = cross2(fan.midNext - fan.midPrev, u2) / denom;
        fan.center = fan.midPrev + u1 * t;
        fan.centerFinite = true;
    } else {
        fan.centerFinite = false;
        fan.parallelDir = u1;
    }

    const int k = static_cast<int>(fan.chain.size());
    fan.collapsed = (k == 1);
    fan.footpoints.resize(k);
    if (fan.collapsed) {
        // The whole sub-polyline maps to the one medial vertex: a polar
        // section of the map, the discrete medial axis endpoint.
        fan.footpoints[0] = pv;
        fan.image = ma.medialVertices[fan.chain[0]].coord;
        ++stats_.collapsedFans;
    } else {
        fan.footpoints[0] = fan.midPrev;
        fan.footpoints[k - 1] = fan.midNext;
        for (int i = 1; i + 1 < k; ++i) {
            bool missed = false;
            fan.footpoints[i] = footpointOnBoundary(
                fan, ma.medialVertices[fan.chain[i]].coord, pv, missed);
            if (missed) ++stats_.fallbackProjections;
        }
        bool missed = false;
        fan.image = mapOnFan(fan, pv, missed);
        if (missed) ++stats_.fallbackProjections;
    }

    fanOfVertex_[v] = static_cast<int>(fans_.size());
    fans_.push_back(std::move(fan));
}

// ─── Evaluation ─────────────────────────────────────────────────────────────

Point MedialAxisMap::mapToAxis(int boundaryVertex, const Point &p) const {
    const MedialFan *fan = fanOf(boundaryVertex);
    if (!fan) return p;
    bool missed = false;
    return mapOnFan(*fan, p, missed);
}

Point MedialAxisMap::mapOnFan(const MedialFan &fan, const Point &p,
                              bool &missed) const {
    missed = false;
    const std::vector<MedialVertex> &mv = axis->medialVertices;
    if (fan.collapsed || fan.chain.size() < 2) {
        return mv[fan.chain.front()].coord;
    }

    // The polar line through p: through the center when it is finite, in the
    // common ray direction when the pencil is parallel.
    Point origin, dir;
    if (fan.centerFinite) {
        origin = fan.center;
        dir = p - origin;
    } else {
        origin = p;
        dir = fan.parallelDir;
    }

    // Monotonicity puts a single chain crossing on the line; if roundoff
    // yields more than one candidate, the closest to p is the right one.
    Point best{0.0, 0.0};
    double bestDist = std::numeric_limits<double>::max();
    bool found = false;
    for (size_t i = 0; i + 1 < fan.chain.size(); ++i) {
        const Point &a = mv[fan.chain[i]].coord;
        const Point &b = mv[fan.chain[i + 1]].coord;
        Point hit;
        double s = 0.0;
        if (!lineSegmentIntersect(origin, dir, a, b, hit, s)) continue;
        if (s < -SEGMENT_PARAM_EPS || s > 1.0 + SEGMENT_PARAM_EPS) continue;
        const double d = dotP(hit - p, hit - p);
        if (d < bestDist) {
            bestDist = d;
            best = hit;
            found = true;
        }
    }
    if (found) return best;

    // The ray missed the chain — a degenerate center or a chain vertex whose
    // circumcenter fell outside the domain. Fall back to the closest point.
    missed = true;
    for (size_t i = 0; i + 1 < fan.chain.size(); ++i) {
        double d = 0.0;
        const Point q = closestOnSegment(p, mv[fan.chain[i]].coord,
                                         mv[fan.chain[i + 1]].coord, d);
        if (d < bestDist) {
            bestDist = d;
            best = q;
        }
    }
    return best;
}

Point MedialAxisMap::footpointOnBoundary(const MedialFan &fan, const Point &m,
                                         const Point &vertexPos,
                                         bool &missed) const {
    missed = false;

    Point origin, dir;
    if (fan.centerFinite) {
        origin = fan.center;
        dir = m - origin;
    } else {
        origin = m;
        dir = fan.parallelDir;
    }

    const std::array<std::array<Point, 2>, 2> segs{{
        {fan.midPrev, vertexPos},
        {vertexPos, fan.midNext},
    }};

    Point best{0.0, 0.0};
    double bestDist = std::numeric_limits<double>::max();
    bool found = false;
    for (const auto &seg : segs) {
        Point hit;
        double s = 0.0;
        if (!lineSegmentIntersect(origin, dir, seg[0], seg[1], hit, s)) continue;
        if (s < -SEGMENT_PARAM_EPS || s > 1.0 + SEGMENT_PARAM_EPS) continue;
        const double d = dotP(hit - m, hit - m);
        if (d < bestDist) {
            bestDist = d;
            best = hit;
            found = true;
        }
    }
    if (found) return best;

    missed = true;
    for (const auto &seg : segs) {
        double d = 0.0;
        const Point q = closestOnSegment(m, seg[0], seg[1], d);
        if (d < bestDist) {
            bestDist = d;
            best = q;
        }
    }
    return best;
}

// ─── Spokes ─────────────────────────────────────────────────────────────────

void MedialAxisMap::buildSpokes() {
    spokes_.clear();
    const std::vector<MedialVertex> &mv = axis->medialVertices;

    for (const MedialFan &fan : fans_) {
        const Point &pv = axis->mesh->vertices[fan.boundaryVertex];
        if (fan.collapsed) {
            const Point &m0 = mv[fan.chain[0]].coord;
            spokes_.push_back({fan.midPrev, m0, true});
            spokes_.push_back({pv, m0, true});
            continue;
        }
        // Everything but the midNext radius, which the next fan emits as its
        // midPrev one (the two fans share that chain vertex).
        const int k = static_cast<int>(fan.chain.size());
        for (int i = 0; i + 1 < k; ++i) {
            spokes_.push_back({fan.footpoints[i], mv[fan.chain[i]].coord, false});
        }
        spokes_.push_back({pv, fan.image, false});
    }
}
