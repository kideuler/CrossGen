#include "medialaxis/MedialAxis.hxx"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <stdexcept>

namespace {

// ─── Circumcircle ───────────────────────────────────────────────────────────

// Circumcircle of a triangle. Returns false for a degenerate (collinear)
// triangle, in which case the diametral circle of the longest edge is written
// out instead: the true circumcircle has unbounded radius there, and a bounded
// placeholder keeps the complex well-formed while the inside filter discards
// whatever Voronoi edges it would have produced.
bool circumcircle(const Point &a, const Point &b, const Point &c,
                  Point &center, double &radius) {
    const Point ab = b - a;
    const Point ac = c - a;
    const double d = 2.0 * cross2(ab, ac);

    // Relative test: |d| is 2 * |ab| * |ac| * sin(angle at a), so comparing it
    // against the edge lengths makes the threshold scale-invariant, unlike a
    // bare |d| < eps.
    const double scale = normP(ab) * normP(ac);
    if (std::abs(d) > 1e-12 * scale && scale > 0.0) {
        const double ab2 = dotP(ab, ab);
        const double ac2 = dotP(ac, ac);
        const Point off{(ac[1] * ab2 - ab[1] * ac2) / d,
                        (ab[0] * ac2 - ac[0] * ab2) / d};
        center = a + off;
        radius = normP(off);
        return true;
    }

    // Degenerate: fall back to the diametral circle of the longest edge.
    const Point *p0 = &a;
    const Point *p1 = &b;
    double best = dotP(ab, ab);
    if (dotP(ac, ac) > best) { best = dotP(ac, ac); p0 = &a; p1 = &c; }
    const Point bc = c - b;
    if (dotP(bc, bc) > best) { p0 = &b; p1 = &c; }
    center = (*p0 + *p1) * 0.5;
    radius = normP(*p1 - *p0) * 0.5;
    return false;
}

// ─── Segment intersection ───────────────────────────────────────────────────

inline double orient(const Point &a, const Point &b, const Point &c) {
    return cross2(b - a, c - a);
}

// Proper (transversal) intersection only. Touching and collinear overlaps read
// as "no crossing", which is what the callers want: a Voronoi vertex sitting
// exactly on the boundary, or a query segment ending on a boundary vertex,
// must not be reported as leaving the domain.
bool properlyIntersect(const Point &p1, const Point &p2,
                       const Point &q1, const Point &q2, double eps) {
    const double d1 = orient(q1, q2, p1);
    const double d2 = orient(q1, q2, p2);
    const double d3 = orient(p1, p2, q1);
    const double d4 = orient(p1, p2, q2);
    const bool straddleP = (d1 > eps && d2 < -eps) || (d1 < -eps && d2 > eps);
    const bool straddleQ = (d3 > eps && d4 < -eps) || (d3 < -eps && d4 > eps);
    return straddleP && straddleQ;
}

// ─── Uniform grid over the boundary polyline ────────────────────────────────
//
// Answers "does this segment leave the domain?" without an O(n) sweep per
// query. Segments are bucketed by their bounding box, which is fine here
// because both the boundary edges and the Voronoi edges being tested are
// short relative to the domain.
class BoundarySegmentGrid {
public:
    void build(std::vector<std::array<Point, 2>> segments, double diagonal) {
        segments_ = std::move(segments);
        eps_ = 1e-12 * diagonal * diagonal;
        if (segments_.empty()) return;

        lo_ = hi_ = segments_[0][0];
        for (const auto &s : segments_) {
            for (const Point &p : s) {
                lo_[0] = std::min(lo_[0], p[0]);
                lo_[1] = std::min(lo_[1], p[1]);
                hi_[0] = std::max(hi_[0], p[0]);
                hi_[1] = std::max(hi_[1], p[1]);
            }
        }
        // Roughly one segment per grid cell.
        const int target = std::max(1, static_cast<int>(
            std::sqrt(static_cast<double>(segments_.size()))));
        nx_ = std::clamp(target, 1, 512);
        ny_ = std::clamp(target, 1, 512);
        cw_ = std::max((hi_[0] - lo_[0]) / nx_, 1e-300);
        ch_ = std::max((hi_[1] - lo_[1]) / ny_, 1e-300);

        buckets_.assign(static_cast<size_t>(nx_) * ny_, {});
        for (int i = 0; i < static_cast<int>(segments_.size()); ++i) {
            int x0, y0, x1, y1;
            cellRange(segments_[i][0], segments_[i][1], x0, y0, x1, y1);
            for (int y = y0; y <= y1; ++y)
                for (int x = x0; x <= x1; ++x)
                    buckets_[static_cast<size_t>(y) * nx_ + x].push_back(i);
        }
        stamp_.assign(segments_.size(), 0);
    }

    bool crossesBoundary(const Point &a, const Point &b) const {
        if (segments_.empty()) return false;

        int x0, y0, x1, y1;
        cellRange(a, b, x0, y0, x1, y1);
        ++epoch_;
        for (int y = y0; y <= y1; ++y) {
            for (int x = x0; x <= x1; ++x) {
                for (int id : buckets_[static_cast<size_t>(y) * nx_ + x]) {
                    if (stamp_[id] == epoch_) continue;
                    stamp_[id] = epoch_;
                    if (properlyIntersect(a, b, segments_[id][0], segments_[id][1], eps_))
                        return true;
                }
            }
        }
        return false;
    }

private:
    void cellRange(const Point &a, const Point &b,
                   int &x0, int &y0, int &x1, int &y1) const {
        const double minx = std::min(a[0], b[0]);
        const double maxx = std::max(a[0], b[0]);
        const double miny = std::min(a[1], b[1]);
        const double maxy = std::max(a[1], b[1]);
        x0 = std::clamp(static_cast<int>((minx - lo_[0]) / cw_), 0, nx_ - 1);
        x1 = std::clamp(static_cast<int>((maxx - lo_[0]) / cw_), 0, nx_ - 1);
        y0 = std::clamp(static_cast<int>((miny - lo_[1]) / ch_), 0, ny_ - 1);
        y1 = std::clamp(static_cast<int>((maxy - lo_[1]) / ch_), 0, ny_ - 1);
    }

    std::vector<std::array<Point, 2>> segments_;
    std::vector<std::vector<int>> buckets_;
    Point lo_{0.0, 0.0}, hi_{0.0, 0.0};
    int nx_ = 1, ny_ = 1;
    double cw_ = 1.0, ch_ = 1.0;
    double eps_ = 0.0;
    mutable std::vector<int> stamp_;
    mutable int epoch_ = 0;
};

} // namespace

// ─── Construction ───────────────────────────────────────────────────────────

MedialAxis::MedialAxis(std::shared_ptr<Mesh> meshIn) : mesh(std::move(meshIn)) {
    if (!mesh) {
        throw std::runtime_error("Mesh pointer cannot be null.");
    }

    Point lo = mesh->vertices.empty() ? Point{0.0, 0.0} : mesh->vertices[0];
    Point hi = lo;
    for (const Point &p : mesh->vertices) {
        lo[0] = std::min(lo[0], p[0]);
        lo[1] = std::min(lo[1], p[1]);
        hi[0] = std::max(hi[0], p[0]);
        hi[1] = std::max(hi[1], p[1]);
    }
    bboxDiagonal_ = std::max(normP(hi - lo), 1e-300);

    buildComplex();
    rebuildBoundaryLinks();
    computeInteriorAngles();

    // Whether each circumcenter is interior is decided once, here, while every
    // cell is still a triangle: a triangle's centroid is guaranteed to lie
    // inside the domain, which gives the ray a known-good starting point. A
    // merge never moves the survivor's center, so the answer stays valid.
    BoundarySegmentGrid grid;
    grid.build(boundarySegments(), bboxDiagonal_);
    for (DualCell &c : cells) {
        const Point centroid = (mesh->vertices[c.verts[0]] +
                                mesh->vertices[c.verts[1]] +
                                mesh->vertices[c.verts[2]]) / 3.0;
        c.centerInside = !grid.crossesBoundary(centroid, c.center);
    }

    buildAxis();
}

void MedialAxis::buildComplex() {
    cells.clear();
    edgeCells.clear();
    stats = MedialAxisStats{};

    const int numTriangles = static_cast<int>(mesh->triangles.size());
    cells.reserve(numTriangles);

    for (int t = 0; t < numTriangles; ++t) {
        const Triangle &tri = mesh->triangles[t];
        DualCell cell;
        cell.verts = {tri[0], tri[1], tri[2]};

        // Normalise the ring to CCW so that ring order, incident-edge order
        // and boundary orientation all agree downstream.
        const Point &a = mesh->vertices[cell.verts[0]];
        const Point &b = mesh->vertices[cell.verts[1]];
        const Point &c = mesh->vertices[cell.verts[2]];
        if (cross2(b - a, c - a) < 0.0) {
            std::swap(cell.verts[1], cell.verts[2]);
        }

        if (!circumcircle(mesh->vertices[cell.verts[0]],
                          mesh->vertices[cell.verts[1]],
                          mesh->vertices[cell.verts[2]],
                          cell.center, cell.radius)) {
            ++stats.degenerateCells;
        }
        cells.push_back(std::move(cell));
    }

    edgeCells.reserve(cells.size() * 2);
    vertexCellCount_.assign(mesh->vertices.size(), 0);
    for (int ci = 0; ci < static_cast<int>(cells.size()); ++ci) {
        const DualCell &cell = cells[ci];
        const int n = static_cast<int>(cell.verts.size());
        for (int i = 0; i < n; ++i) {
            ++vertexCellCount_[cell.verts[i]];
            const MeshEdgeKey key(cell.verts[i], cell.verts[(i + 1) % n]);
            auto &slots = edgeCells.try_emplace(key, std::array<int, 2>{-1, -1}).first->second;
            if (slots[0] == -1)      slots[0] = ci;
            else if (slots[1] == -1) slots[1] = ci;
            // A third incident cell means the input is not a manifold surface;
            // ignore it rather than corrupt the adjacency.
        }
    }
}

std::vector<std::array<Point, 2>> MedialAxis::boundarySegments() const {
    std::vector<std::array<Point, 2>> segments;
    for (const auto &kv : edgeCells) {
        const std::array<int, 2> &slots = kv.second;
        const bool aActive = slots[0] >= 0 && cells[slots[0]].active;
        const bool bActive = slots[1] >= 0 && cells[slots[1]].active;
        if (aActive != bActive) {
            segments.push_back({mesh->vertices[kv.first.a], mesh->vertices[kv.first.b]});
        }
    }
    return segments;
}

void MedialAxis::rebuildBoundaryLinks() {
    const int numVertices = static_cast<int>(mesh->vertices.size());
    boundaryLinks.assign(numVertices, BoundaryVertex{});
    for (int i = 0; i < numVertices; ++i) {
        boundaryLinks[i].vertexIndex = i;
    }
    stats.nonManifoldBoundaryVertices = 0;

    // A cell's ring is CCW, so the domain lies to the LEFT of every ring edge
    // v[i] -> v[i+1]. Orienting constrained ring edges that way makes the
    // outer loop CCW and hole loops CW, with the interior on the left of
    // `next` in both cases.
    for (int ci = 0; ci < static_cast<int>(cells.size()); ++ci) {
        const DualCell &cell = cells[ci];
        if (!cell.active) continue;
        const int n = static_cast<int>(cell.verts.size());
        for (int i = 0; i < n; ++i) {
            if (!isConstrained(ci, i)) continue;
            const int u = cell.verts[i];
            const int v = cell.verts[(i + 1) % n];
            if (boundaryLinks[u].nextBoundaryVertex != -1 ||
                boundaryLinks[v].prevBoundaryVertex != -1) {
                // The vertex already carries a link, so it is pinched between
                // two boundary loops and a single next/prev pair cannot
                // describe it. Report rather than silently overwrite.
                ++stats.nonManifoldBoundaryVertices;
                continue;
            }
            boundaryLinks[u].nextBoundaryVertex = v;
            boundaryLinks[v].prevBoundaryVertex = u;
        }
    }
}

void MedialAxis::computeInteriorAngles() {
    const int numVertices = static_cast<int>(mesh->vertices.size());
    interiorAngle.assign(numVertices, std::numeric_limits<double>::quiet_NaN());

    for (int v = 0; v < numVertices; ++v) {
        const int prev = boundaryLinks[v].prevBoundaryVertex;
        const int next = boundaryLinks[v].nextBoundaryVertex;
        if (prev == -1 || next == -1) continue;

        const Point in  = mesh->vertices[v] - mesh->vertices[prev];
        const Point out = mesh->vertices[next] - mesh->vertices[v];
        if (normP(in) < 1e-300 || normP(out) < 1e-300) continue;

        // `next` keeps the interior on the left, so the turn is positive at a
        // convex corner, and the interior angle is pi minus that turn. The
        // result lands in (0, 2*pi): below pi convex, above pi reflex.
        const double turn = std::atan2(cross2(in, out), dotP(in, out));
        interiorAngle[v] = M_PI - turn;
    }
}

// ─── Axis extraction ────────────────────────────────────────────────────────

void MedialAxis::buildAxis() {
    medialVertices.clear();
    medialEdges.clear();
    polyLines.clear();
    polyLineIsCycle.clear();
    stats.centersOutsideDomain = 0;
    stats.edgesCrossingBoundary = 0;

    for (DualCell &c : cells) c.medialVertex = -1;

    for (int ci = 0; ci < static_cast<int>(cells.size()); ++ci) {
        DualCell &cell = cells[ci];
        if (!cell.active) continue;
        if (!cell.centerInside) ++stats.centersOutsideDomain;

        cell.medialVertex = static_cast<int>(medialVertices.size());
        MedialVertex mv;
        mv.coord = cell.center;
        mv.radius = cell.radius;
        mv.cell = ci;
        mv.dualIsTriangle = cell.dualIsTriangle();
        mv.insideDomain = cell.centerInside;
        medialVertices.push_back(std::move(mv));
    }

    BoundarySegmentGrid grid;
    grid.build(boundarySegments(), bboxDiagonal_);

    // The axis is the part of the Voronoi diagram that lies inside the domain,
    // so a candidate edge survives only when both of its endpoints are
    // interior and the segment between them does not cross the boundary.
    // Dropping an edge here means the boundary sampling was too coarse for the
    // dual to be a faithful medial axis; the count is reported in `stats`.
    std::unordered_map<MeshEdgeKey, int, MeshEdgeKeyHash> medialEdgeOf;
    medialEdgeOf.reserve(edgeCells.size());

    // Walking each cell's own ring in order gives incidentEdges in CCW order.
    for (int ci = 0; ci < static_cast<int>(cells.size()); ++ci) {
        const DualCell &cell = cells[ci];
        if (!cell.active) continue;
        const int m = cell.medialVertex;
        const int n = static_cast<int>(cell.verts.size());

        for (int i = 0; i < n; ++i) {
            const int nb = neighborAcross(ci, i);
            if (nb < 0) continue;

            const MeshEdgeKey key(cell.verts[i], cell.verts[(i + 1) % n]);
            auto it = medialEdgeOf.find(key);
            if (it == medialEdgeOf.end()) {
                const bool inside = cell.centerInside && cells[nb].centerInside &&
                                    !grid.crossesBoundary(cell.center, cells[nb].center);
                if (!inside) {
                    medialEdgeOf.emplace(key, -1);
                    ++stats.edgesCrossingBoundary;
                    continue;
                }
                const int id = static_cast<int>(medialEdges.size());
                medialEdges.push_back({m, cells[nb].medialVertex});
                medialEdgeOf.emplace(key, id);
                medialVertices[m].incidentEdges.push_back(id);
            } else if (it->second >= 0) {
                medialVertices[m].incidentEdges.push_back(it->second);
            }
        }
    }

    for (MedialVertex &mv : medialVertices) {
        mv.degree = static_cast<int>(mv.incidentEdges.size());
    }
    // Neighbour sets index into medialVertices, so they can only be filled
    // once every medial vertex exists.
    for (int e = 0; e < static_cast<int>(medialEdges.size()); ++e) {
        medialVertices[medialEdges[e][0]].neighbors.insert(medialEdges[e][1]);
        medialVertices[medialEdges[e][1]].neighbors.insert(medialEdges[e][0]);
    }
}

// ─── Complex queries ────────────────────────────────────────────────────────

std::array<int, 2> MedialAxis::ringEdge(int cell, int edgeIdx) const {
    const std::vector<int> &verts = cells[cell].verts;
    const int n = static_cast<int>(verts.size());
    return {verts[edgeIdx], verts[(edgeIdx + 1) % n]};
}

int MedialAxis::neighborAcross(int cell, int edgeIdx) const {
    const std::array<int, 2> e = ringEdge(cell, edgeIdx);
    const auto it = edgeCells.find(MeshEdgeKey(e[0], e[1]));
    if (it == edgeCells.end()) return -1;
    for (int c : it->second) {
        if (c >= 0 && c != cell && cells[c].active) return c;
    }
    return -1;
}

bool MedialAxis::isConstrained(int cell, int edgeIdx) const {
    return neighborAcross(cell, edgeIdx) < 0;
}

int MedialAxis::cellDegree(int cell) const {
    const int n = static_cast<int>(cells[cell].verts.size());
    int degree = 0;
    for (int i = 0; i < n; ++i) {
        if (neighborAcross(cell, i) >= 0) ++degree;
    }
    return degree;
}

int MedialAxis::vertexCellCount(int v) const {
    if (v < 0 || v >= static_cast<int>(vertexCellCount_.size())) return 0;
    return vertexCellCount_[v];
}

bool MedialAxis::isSharpCorner(int v, double threshold) const {
    const double a = interiorAngle[v];
    return std::isfinite(a) && a < threshold;
}

bool MedialAxis::isConcaveCorner(int v) const {
    const double a = interiorAngle[v];
    return std::isfinite(a) && a > M_PI;
}

// ─── Cell merging ───────────────────────────────────────────────────────────

bool MedialAxis::mergeCells(int survivor, int absorbed) {
    DualCell &S = cells[survivor];
    DualCell &A = cells[absorbed];
    if (!S.active || !A.active || survivor == absorbed) return false;

    const int nA = static_cast<int>(A.verts.size());
    const int nS = static_cast<int>(S.verts.size());

    // Exactly one shared ring edge, otherwise the union of the two rings is
    // not a simple polygon.
    int p = -1;
    int sharedCount = 0;
    for (int i = 0; i < nA; ++i) {
        if (neighborAcross(absorbed, i) == survivor) { p = i; ++sharedCount; }
    }
    if (sharedCount != 1) {
        ++stats.mergesRejectedMultiEdge;
        return false;
    }

    const int a = A.verts[p];
    const int b = A.verts[(p + 1) % nA];

    // The shared edge runs a -> b in A, hence b -> a in the CCW ring of the
    // cell on the other side.
    int q = -1;
    for (int i = 0; i < nS; ++i) {
        if (S.verts[i] == b && S.verts[(i + 1) % nS] == a) { q = i; break; }
    }
    if (q == -1) {
        ++stats.mergesRejectedMultiEdge;
        return false;
    }

    // Any vertex shared beyond the two edge endpoints pinches the union, so
    // the merged ring would visit it twice.
    {
        std::unordered_set<int> sVerts(S.verts.begin(), S.verts.end());
        for (int v : A.verts) {
            if (v != a && v != b && sVerts.count(v)) {
                ++stats.mergesRejectedMultiEdge;
                return false;
            }
        }
    }

    // Splice A's chain from b round to a into S's ring in place of edge b -> a.
    std::vector<int> merged;
    merged.reserve(nS + nA - 2);
    for (int j = 0; j < nS; ++j) {
        const int idx = (q + j) % nS;
        merged.push_back(S.verts[idx]);
        if (j == 0) {
            for (int t = 2; t < nA; ++t) {
                merged.push_back(A.verts[(p + t) % nA]);
            }
        }
    }

    // Re-point A's remaining ring edges at S, and drop the shared one.
    for (int i = 0; i < nA; ++i) {
        const MeshEdgeKey key(A.verts[i], A.verts[(i + 1) % nA]);
        if (i == p) {
            edgeCells.erase(key);
            continue;
        }
        auto it = edgeCells.find(key);
        if (it == edgeCells.end()) continue;
        for (int &slot : it->second) {
            if (slot == absorbed) slot = survivor;
        }
    }

    // a and b lose one incident cell; every other vertex of A simply changes
    // owner, so its count is unaffected.
    --vertexCellCount_[a];
    --vertexCellCount_[b];

    S.verts = std::move(merged);
    S.timesFused += A.timesFused + 1;

    // The survivor keeps its circle, so the merged ring is only approximately
    // cocircular. Record the drift: a real simplification pass would instead
    // project the absorbed vertices onto this circle and keep it exact.
    double deviation = S.maxRadiusDeviation;
    if (S.radius > 0.0) {
        for (int v : S.verts) {
            const double d = std::abs(normP(mesh->vertices[v] - S.center) - S.radius) / S.radius;
            deviation = std::max(deviation, d);
        }
    }
    S.maxRadiusDeviation = deviation;

    A.active = false;
    A.medialVertex = -1;
    A.verts.clear();
    return true;
}

void MedialAxis::deduplicateMedialVertices(double tolerance) {
    stats.mergedCells = 0;
    stats.mergesRejectedMultiEdge = 0;

    // Leaves first, then chains, then junctions: pruning noise leaves before
    // touching junctions keeps a spurious leaf from dragging a real branch
    // point with it.
    const int phases[3][2] = {{1, 1}, {2, 2}, {3, std::numeric_limits<int>::max()}};

    std::vector<char> queued(cells.size(), 0);
    std::deque<int> queue;

    auto enqueue = [&](int ci) {
        if (ci < 0 || !cells[ci].active || queued[ci]) return;
        queued[ci] = 1;
        queue.push_back(ci);
    };

    for (const auto &phase : phases) {
        const int lo = phase[0];
        const int hi = phase[1];

        std::fill(queued.begin(), queued.end(), 0);
        queue.clear();
        for (int ci = 0; ci < static_cast<int>(cells.size()); ++ci) {
            if (cells[ci].active) enqueue(ci);
        }

        // A worklist rather than a rescan-from-zero: each merge only disturbs
        // the survivor and its neighbours, so only those need revisiting.
        while (!queue.empty()) {
            const int ci = queue.front();
            queue.pop_front();
            queued[ci] = 0;
            if (!cells[ci].active) continue;

            const int degree = cellDegree(ci);
            if (degree < lo || degree > hi) continue;

            // Closest neighbour within tolerance of the local medial radius.
            int target = -1;
            double bestDist = std::numeric_limits<double>::max();
            const int n = static_cast<int>(cells[ci].verts.size());
            for (int i = 0; i < n; ++i) {
                const int nb = neighborAcross(ci, i);
                if (nb < 0) continue;
                const double dist = normP(cells[ci].center - cells[nb].center);
                const double localRadius = std::max(cells[ci].radius, cells[nb].radius);
                if (dist < tolerance * localRadius && dist < bestDist) {
                    bestDist = dist;
                    target = nb;
                }
            }
            if (target < 0) continue;

            if (!mergeCells(target, ci)) continue;
            ++stats.mergedCells;

            enqueue(target);
            const int tn = static_cast<int>(cells[target].verts.size());
            for (int i = 0; i < tn; ++i) enqueue(neighborAcross(target, i));
        }
    }

    rebuildBoundaryLinks();
    computeInteriorAngles();
    buildAxis();
}

// ─── Branches ───────────────────────────────────────────────────────────────

void MedialAxis::createPolylines() {
    polyLines.clear();
    polyLineIsCycle.clear();

    const int numEdges = static_cast<int>(medialEdges.size());
    std::vector<char> used(numEdges, 0);

    auto otherEnd = [&](int edge, int from) {
        return medialEdges[edge][0] == from ? medialEdges[edge][1] : medialEdges[edge][0];
    };

    // Walk from `start` along `edge` until reaching a vertex whose degree is
    // not 2. Marking edges rather than vertices is what keeps a branch from
    // being traced once from each end.
    auto walk = [&](int start, int edge) {
        std::vector<int> path{start};
        int cur = start;
        int ce = edge;
        while (true) {
            used[ce] = 1;
            const int next = otherEnd(ce, cur);
            path.push_back(next);
            if (medialVertices[next].degree != 2) break;

            int following = -1;
            for (int e : medialVertices[next].incidentEdges) {
                if (e != ce) { following = e; break; }
            }
            if (following == -1 || used[following]) break;
            cur = next;
            ce = following;
        }
        return path;
    };

    for (int i = 0; i < static_cast<int>(medialVertices.size()); ++i) {
        const MedialVertex &mv = medialVertices[i];
        if (mv.degree == 2) continue;
        if (mv.degree == 0) {
            // An isolated medial vertex still belongs in the output, otherwise
            // it silently disappears from anything iterating branches.
            polyLines.push_back({i});
            polyLineIsCycle.push_back(false);
            continue;
        }
        for (int e : mv.incidentEdges) {
            if (used[e]) continue;
            polyLines.push_back(walk(i, e));
            // A branch can leave a junction and come straight back to it, so
            // closure has to be tested here too, not only in the cycle pass.
            polyLineIsCycle.push_back(polyLines.back().size() > 1 &&
                                      polyLines.back().front() == polyLines.back().back());
        }
    }

    // Whatever is left is a closed cycle of degree-2 vertices. These are real
    // medial structure -- every hole in the domain produces one -- and seeding
    // only from vertices of degree != 2 would never reach them.
    for (int e = 0; e < numEdges; ++e) {
        if (used[e]) continue;
        std::vector<int> path = walk(medialEdges[e][0], e);
        polyLines.push_back(std::move(path));
        polyLineIsCycle.push_back(polyLines.back().size() > 1 &&
                                  polyLines.back().front() == polyLines.back().back());
    }
}
