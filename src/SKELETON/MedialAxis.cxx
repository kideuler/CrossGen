#include "SKELETON/MedialAxis.hxx"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <queue>
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

double distToSegment(const Point &p, const Point &a, const Point &b) {
    const Point e = b - a;
    const double len2 = dotP(e, e);
    const double t = len2 > 0.0 ? std::clamp(dotP(p - a, e) / len2, 0.0, 1.0) : 0.0;
    return normP(p - (a + e * t));
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

    // The input boundary, frozen before simplify() moves anything: x0 for the
    // displacement budget of Step 4, and the loop order for walking the
    // original polyline through vertices 4.2 later removes.
    originalPosition = mesh->vertices;
    originalNext.assign(mesh->vertices.size(), -1);
    for (int v = 0; v < static_cast<int>(mesh->vertices.size()); ++v) {
        originalNext[v] = boundaryLinks[v].nextBoundaryVertex;
    }
    locked.assign(mesh->vertices.size(), 0);

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
    medialEdgeDual.clear();
    polyLines.clear();
    polyLineEdges.clear();
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
                medialEdgeDual.push_back({key.a, key.b});
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
        return false;
    }

    // Any vertex shared beyond the two edge endpoints pinches the union, so
    // the merged ring would visit it twice.
    {
        std::unordered_set<int> sVerts(S.verts.begin(), S.verts.end());
        for (int v : A.verts) {
            if (v != a && v != b && sVerts.count(v)) {
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

    // The survivor keeps its circle. Every caller has put the absorbed
    // cell's private vertices on it first (Step 4.1 projects them, the
    // cocircular merge finds them there), so this only records rounding.
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

// ─── Step 4: simplification ─────────────────────────────────────────────────
//
// Topological noise -- endpoints and intersection points crowded together by
// a dense or slightly irregular sampling (Fig. 10, 15) -- each puts a polar
// region into phi, so it has to go before the map is built. The axis is
// never edited directly: boundary vertices are moved or removed so that the
// complex, and hence the axis and phi, stay consistent. Every move is
// charged against x0, so no boundary vertex ever ends up farther than
// delta_max from where the input put it.

void MedialAxis::simplify(const SimplifyOptions &opts) {
    simplifyStats = SimplifyStats{};
    const double h = opts.targetSize > 0.0 ? opts.targetSize : defaultTargetSize();
    const double budget = opts.maxDisplacement >= 0.0 ? opts.maxDisplacement : 0.05 * h;
    simplifyStats.budget = budget;

    // Locked features (Step 1.4), from the input geometry: a turning angle
    // above alpha_feat either way is a CAD corner, convex or concave, and
    // Step 4 may neither move it nor flatten it away. A pinched vertex has no
    // single loop to slide along, so it is locked as well. The boundary is
    // only boundary loops here -- no interfaces, hence no junctions, tips or
    // contacts to add.
    const int nv = static_cast<int>(mesh->vertices.size());
    std::vector<int> originalPrev(nv, -1);
    for (int v = 0; v < nv; ++v) {
        if (originalNext[v] >= 0) originalPrev[originalNext[v]] = v;
    }
    locked.assign(nv, 0);
    for (int v = 0; v < nv; ++v) {
        const int p = originalPrev[v];
        const int n = originalNext[v];
        if (p < 0 && n < 0) continue;
        if (p < 0 || n < 0 || boundaryLinks[v].nextBoundaryVertex < 0 ||
            boundaryLinks[v].prevBoundaryVertex < 0) {
            locked[v] = 1;
            continue;
        }
        const Point in  = originalPosition[v] - originalPosition[p];
        const Point out = originalPosition[n] - originalPosition[v];
        const double turn = std::atan2(cross2(in, out), dotP(in, out));
        if (std::abs(turn) > opts.featureAngle) locked[v] = 1;
    }

    // Regions are the connected components of the complex. "Merging would
    // leave the region with a single face" is checked against these counts:
    // keeping two faces keeps at least one medial edge, without which there
    // is no axis to lay templates on.
    regionOf_.assign(cells.size(), -1);
    std::vector<int> regionCells;
    for (int c0 = 0; c0 < static_cast<int>(cells.size()); ++c0) {
        if (!cells[c0].active || regionOf_[c0] >= 0) continue;
        const int region = static_cast<int>(regionCells.size());
        regionCells.push_back(0);
        std::deque<int> queue{c0};
        regionOf_[c0] = region;
        while (!queue.empty()) {
            const int c = queue.front();
            queue.pop_front();
            ++regionCells[region];
            for (int i = 0; i < static_cast<int>(cells[c].verts.size()); ++i) {
                const int nb = neighborAcross(c, i);
                if (nb >= 0 && regionOf_[nb] < 0) {
                    regionOf_[nb] = region;
                    queue.push_back(nb);
                }
            }
        }
    }

    // 4.1 to a fixpoint, then 4.2 (the paper's order), and round again while
    // either changes anything: deleting a 2-gon turns its neighbour into a
    // new endpoint that 4.1 has not seen. The cocircular merge leads each
    // round, since a cell that 4.1 grew can come to share its circle with a
    // neighbour.
    for (int round = 0; round < opts.maxRounds; ++round) {
        ++simplifyStats.rounds;
        const int changed = mergeCocircular(opts.cocircularTolerance * budget, regionCells) +
                            removeEndpoints(budget, regionCells) +
                            flattenArcs(budget, regionCells);
        if (changed == 0) break;
    }

    for (int c = 0; c < static_cast<int>(cells.size()); ++c) {
        if (cells[c].active && cellDegree(c) == 1) ++simplifyStats.endpointsKept;
    }

    rebuildBoundaryLinks();
    computeInteriorAngles();
    for (int v = 0; v < nv; ++v) {
        if (boundaryLinks[v].nextBoundaryVertex < 0) continue;
        simplifyStats.maxDisplacement = std::max(
            simplifyStats.maxDisplacement, normP(mesh->vertices[v] - originalPosition[v]));
    }
    buildAxis();
}

bool MedialAxis::crossesBoundary(const Point &a, const Point &b,
                                 const std::unordered_set<int> &ignore) const {
    const double eps = 1e-12 * bboxDiagonal_ * bboxDiagonal_;
    for (int u = 0; u < static_cast<int>(boundaryLinks.size()); ++u) {
        const int w = boundaryLinks[u].nextBoundaryVertex;
        if (w < 0 || ignore.count(u) || ignore.count(w)) continue;
        if (properlyIntersect(a, b, mesh->vertices[u], mesh->vertices[w], eps)) return true;
    }
    return false;
}

// Step 2.3: a cocircular set of samples has no unique Delaunay
// triangulation, and whichever one Triangle picked splits a single medial
// vertex into several sitting on top of each other, joined by medial edges
// of no length. Merging the cells back costs nothing -- every vertex is
// already on the survivor's circle to within `tolerance`, so nothing moves.
int MedialAxis::mergeCocircular(double tolerance, std::vector<int> &regionCells) {
    std::vector<char> queued(cells.size(), 0);
    std::deque<int> queue;
    auto enqueue = [&](int c) {
        if (c < 0 || !cells[c].active || queued[c]) return;
        queued[c] = 1;
        queue.push_back(c);
    };
    for (int c = 0; c < static_cast<int>(cells.size()); ++c) enqueue(c);

    int merged = 0;
    while (!queue.empty()) {
        const int f = queue.front();
        queue.pop_front();
        queued[f] = 0;
        if (!cells[f].active) continue;

        for (int i = 0; i < static_cast<int>(cells[f].verts.size()); ++i) {
            const int g = neighborAcross(f, i);
            if (g < 0) continue;
            const DualCell &F = cells[f];
            const DualCell &G = cells[g];
            // Rounding is all the slack an exactly cocircular input needs.
            const double tol = std::max(tolerance, 1e-12 * std::max(F.radius, G.radius));
            if (std::abs(F.radius - G.radius) > tol) continue;
            if (regionCells[regionOf_[f]] <= 2) continue;
            // The samples themselves decide it, not the two centres: a
            // degenerate cell's placeholder circle agrees with nothing, and
            // two nearly equal circles through the same samples can have
            // centres much farther apart than the samples are off either.
            bool onCircle = true;
            for (int v : F.verts) {
                if (std::abs(normP(mesh->vertices[v] - G.center) - G.radius) > tol) {
                    onCircle = false;
                    break;
                }
            }
            if (!onCircle || !mergeCells(g, f)) continue;
            --regionCells[regionOf_[f]];
            ++merged;
            enqueue(g);
            break;
        }
    }
    simplifyStats.cocircularMerges += merged;
    return merged;
}

// Step 4.1 (Fig. 11, 12): an endpoint F is absorbed into its one interior
// neighbour G by projecting F's private vertices radially onto G's circle.
// Those vertices belong to F alone -- both their boundary edges are F's -- so
// the move changes no other part of the axis, and the merged ring is exactly
// cocircular on G's circle, which G keeps. Absorbing an endpoint can make G
// one, so merges cascade along a branch until the budget stops them: noisy
// branches are pruned and the polar sections along a round boundary merge
// into one (Fig. 15).
int MedialAxis::removeEndpoints(double budget, std::vector<int> &regionCells) {
    // The unique interior neighbour of an endpoint, and the ring index of
    // the edge shared with it.
    auto soleNeighbour = [&](int f, int &edge) {
        const int n = static_cast<int>(cells[f].verts.size());
        int g = -1;
        for (int i = 0; i < n; ++i) {
            const int nb = neighborAcross(f, i);
            if (nb < 0) continue;
            if (g >= 0) return -1;
            g = nb;
            edge = i;
        }
        return g;
    };

    // Where F's private vertices land on G's circle, and the largest
    // cumulative displacement that costs. Negative when F cannot be absorbed
    // at all: a locked feature among them, or a vertex some other cell also
    // uses (a pinch).
    auto plan = [&](int f, int g, int edge, std::vector<int> &Q, std::vector<Point> &xNew) {
        Q.clear();
        xNew.clear();
        const DualCell &F = cells[f];
        const DualCell &G = cells[g];
        const int n = static_cast<int>(F.verts.size());
        double cost = 0.0;
        for (int t = 2; t < n; ++t) {
            const int q = F.verts[(edge + t) % n];
            if (locked[q] || vertexCellCount(q) != 1) return -1.0;
            const Point d = mesh->vertices[q] - G.center;
            const double len = normP(d);
            if (!(len > 0.0)) return -1.0;
            const Point x = G.center + d * (G.radius / len);
            cost = std::max(cost, normP(x - originalPosition[q]));
            Q.push_back(q);
            xNew.push_back(x);
        }
        return cost;
    };

    using Item = std::pair<double, int>;
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> heap;
    std::vector<int> Q;
    std::vector<Point> xNew;
    auto push = [&](int f) {
        int edge = -1;
        const int g = soleNeighbour(f, edge);
        if (g < 0) return;
        const double cost = plan(f, g, edge, Q, xNew);
        if (cost >= 0.0 && cost <= budget) heap.push({cost, f});
    };
    for (int f = 0; f < static_cast<int>(cells.size()); ++f) {
        if (cells[f].active) push(f);
    }

    const double eps = 1e-12 * bboxDiagonal_ * bboxDiagonal_;
    int merges = 0;
    while (!heap.empty()) {
        const int f = heap.top().second;
        heap.pop();
        if (!cells[f].active) continue;
        int edge = -1;
        const int g = soleNeighbour(f, edge);
        if (g < 0) continue;
        if (regionCells[regionOf_[f]] <= 2) continue;
        // Re-planned on pop: F may have grown since it was queued, and the
        // queue order is only a heuristic ([R]: least displacement first).
        const double cost = plan(f, g, edge, Q, xNew);
        if (cost < 0.0 || cost > budget) continue;

        std::vector<Point> xOld;
        for (size_t k = 0; k < Q.size(); ++k) {
            xOld.push_back(mesh->vertices[Q[k]]);
            mesh->vertices[Q[k]] = xNew[k];
        }
        const std::unordered_set<int> moved(Q.begin(), Q.end());

        // The merged ring F u G, spliced the way mergeCells will splice it,
        // must stay a simple polygon. Only an edge touching a moved vertex
        // can have become a crossing.
        const DualCell &F = cells[f];
        const DualCell &G = cells[g];
        const int nF = static_cast<int>(F.verts.size());
        const int a = F.verts[edge];
        const int b = F.verts[(edge + 1) % nF];
        std::vector<int> ring;
        for (size_t j = 0; j < G.verts.size(); ++j) {
            if (G.verts[j] != b) continue;
            for (size_t k = 0; k < G.verts.size(); ++k) {
                ring.push_back(G.verts[(j + k) % G.verts.size()]);
                if (k == 0) {
                    for (int t = 2; t < nF; ++t) ring.push_back(F.verts[(edge + t) % nF]);
                }
            }
            break;
        }
        bool ok = !ring.empty();
        const int nr = static_cast<int>(ring.size());
        double area2 = 0.0;
        for (int i = 0; i < nr && ok; ++i) {
            const Point &p0 = mesh->vertices[ring[i]];
            const Point &p1 = mesh->vertices[ring[(i + 1) % nr]];
            area2 += cross2(p0, p1);
            if (!moved.count(ring[i]) && !moved.count(ring[(i + 1) % nr])) continue;
            for (int j = 0; j < nr; ++j) {
                if (j == i || j == (i + 1) % nr || (j + 1) % nr == i) continue;
                if (properlyIntersect(p0, p1, mesh->vertices[ring[j]],
                                      mesh->vertices[ring[(j + 1) % nr]], eps)) {
                    ok = false;
                    break;
                }
            }
        }
        ok = ok && area2 > 0.0;

        // Nor may a moved boundary edge cut through the rest of the boundary.
        for (size_t k = 0; k < Q.size() && ok; ++k) {
            const int q = Q[k];
            for (int nbr : {boundaryLinks[q].prevBoundaryVertex,
                            boundaryLinks[q].nextBoundaryVertex}) {
                if (nbr < 0) continue;
                const std::unordered_set<int> ignore{q, nbr};
                if (crossesBoundary(mesh->vertices[q], mesh->vertices[nbr], ignore)) {
                    ok = false;
                    break;
                }
            }
        }

        // The optional check of Step 4.1: a moved vertex must not enter the
        // circle of another cell next to G, or that cell's centre would stop
        // being a medial vertex of the new boundary.
        for (int i = 0; i < static_cast<int>(G.verts.size()) && ok; ++i) {
            const int h = neighborAcross(g, i);
            if (h < 0 || h == f) continue;
            for (int q : Q) {
                if (normP(mesh->vertices[q] - cells[h].center) <
                    cells[h].radius * (1.0 - 1e-9)) {
                    ok = false;
                    break;
                }
            }
        }

        if (!ok || !mergeCells(g, f)) {
            for (size_t k = 0; k < Q.size(); ++k) mesh->vertices[Q[k]] = xOld[k];
            continue;
        }
        --regionCells[regionOf_[f]];
        ++merges;
        push(g);
    }
    simplifyStats.endpointMerges += merges;
    return merges;
}

// Step 4.2 (Fig. 13): a chain of two or more boundary edges in a merged cell
// is a polar section -- phi collapses the whole chain onto the cell's
// medial vertex. Where the input boundary under the chain is straight to
// within delta_max, the chain is a circular arc only because the sampling
// put it on one: its interior vertices are deleted and the chord takes its
// place, so the contact becomes a single boundary edge and the polar region
// disappears.
int MedialAxis::flattenArcs(double budget, std::vector<int> &regionCells) {
    int changes = 0;
    for (int f = 0; f < static_cast<int>(cells.size()); ++f) {
        bool again = true;
        while (again && cells[f].active && !cells[f].dualIsTriangle()) {
            again = false;
            const std::vector<int> ring = cells[f].verts;
            const int n = static_cast<int>(ring.size());
            std::vector<char> constrained(n, 0);
            int start = -1;
            for (int i = 0; i < n; ++i) {
                constrained[i] = isConstrained(f, i) ? 1 : 0;
                if (!constrained[i]) start = i;
            }
            if (start < 0) break;   // a cell with no interior edge: nothing to keep

            // Maximal runs of constrained ring edges, scanned from just after
            // an interior one so that no run straddles the start.
            for (int k = 1; k <= n && !again; ++k) {
                const int s = (start + k) % n;
                if (!constrained[s] || constrained[(s + n - 1) % n]) continue;
                int len = 0;
                while (len < n && constrained[(s + len) % n]) ++len;
                if (len < 2) continue;

                const int first = ring[s];
                const int last = ring[(s + len) % n];
                std::vector<int> interior;
                bool admissible = true;
                for (int t = 1; t < len; ++t) {
                    const int v = ring[(s + t) % n];
                    if (locked[v] || vertexCellCount(v) != 1) admissible = false;
                    interior.push_back(v);
                }
                if (!admissible) continue;

                // The deviation the chord would cost, against the input
                // boundary it replaces ([R]: original positions and the
                // original polyline between the chain ends) and against the
                // chain as it now stands.
                const Point &pa = mesh->vertices[first];
                const Point &pb = mesh->vertices[last];
                double d = 0.0;
                for (const Point &p : originalPolyline(first, last)) {
                    d = std::max(d, distToSegment(p, pa, pb));
                }
                for (int v : interior) {
                    d = std::max(d, distToSegment(mesh->vertices[v], pa, pb));
                }
                if (d > budget) continue;

                // The whole boundary of an endpoint: the cell degenerates to
                // the 2-gon (first, last) and is deleted, its interior edge
                // becoming a boundary edge of the neighbour ([R]).
                const bool twoGon = (len == n - 1);
                if (twoGon && regionCells[regionOf_[f]] <= 2) continue;
                const MeshEdgeKey chord(first, last);
                if (!twoGon && edgeCells.count(chord)) continue;

                std::unordered_set<int> ignore(interior.begin(), interior.end());
                ignore.insert(first);
                ignore.insert(last);
                if (crossesBoundary(pa, pb, ignore)) continue;

                // Commit.
                for (int t = 0; t < len; ++t) {
                    edgeCells.erase(MeshEdgeKey(ring[(s + t) % n], ring[(s + t + 1) % n]));
                }
                for (int v : interior) {
                    vertexCellCount_[v] = 0;
                    boundaryLinks[v].nextBoundaryVertex = -1;
                    boundaryLinks[v].prevBoundaryVertex = -1;
                }
                boundaryLinks[first].nextBoundaryVertex = last;
                boundaryLinks[last].prevBoundaryVertex = first;

                if (twoGon) {
                    for (int &slot : edgeCells[chord]) {
                        if (slot == f) slot = -1;
                    }
                    --vertexCellCount_[first];
                    --vertexCellCount_[last];
                    cells[f].active = false;
                    cells[f].medialVertex = -1;
                    cells[f].verts.clear();
                    --regionCells[regionOf_[f]];
                    ++simplifyStats.deletedTwoGons;
                } else {
                    edgeCells[chord] = {f, -1};
                    std::vector<int> kept;
                    for (int v : ring) {
                        if (!ignore.count(v) || v == first || v == last) kept.push_back(v);
                    }
                    cells[f].verts = std::move(kept);
                }
                ++simplifyStats.flattenedChains;
                simplifyStats.removedVertices += static_cast<int>(interior.size());
                ++changes;
                again = true;
            }
        }
    }
    return changes;
}

// ─── Branches ───────────────────────────────────────────────────────────────

void MedialAxis::createPolylines() {
    polyLines.clear();
    polyLineEdges.clear();
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
        std::vector<int> pathEdges;
        int cur = start;
        int ce = edge;
        while (true) {
            used[ce] = 1;
            pathEdges.push_back(ce);
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
        polyLines.push_back(std::move(path));
        polyLineEdges.push_back(std::move(pathEdges));
        // A branch can leave a junction and come straight back to it, so
        // closure has to be tested for every branch, not only in the cycle
        // pass.
        polyLineIsCycle.push_back(polyLines.back().size() > 1 &&
                                  polyLines.back().front() == polyLines.back().back());
    };

    for (int i = 0; i < static_cast<int>(medialVertices.size()); ++i) {
        const MedialVertex &mv = medialVertices[i];
        if (mv.degree == 2) continue;
        if (mv.degree == 0) {
            // An isolated medial vertex still belongs in the output, otherwise
            // it silently disappears from anything iterating branches.
            polyLines.push_back({i});
            polyLineEdges.push_back({});
            polyLineIsCycle.push_back(false);
            continue;
        }
        for (int e : mv.incidentEdges) {
            if (!used[e]) walk(i, e);
        }
    }

    // Whatever is left is a closed cycle of degree-2 vertices. These are real
    // medial structure -- every hole in the domain produces one -- and seeding
    // only from vertices of degree != 2 would never reach them.
    for (int e = 0; e < numEdges; ++e) {
        if (!used[e]) walk(medialEdges[e][0], e);
    }
}

std::vector<Point> MedialAxis::originalPolyline(int a, int b) const {
    std::vector<Point> pts{originalPosition[a]};
    int v = a;
    int guard = static_cast<int>(originalNext.size()) + 1;
    while (guard-- > 0) {
        v = originalNext[v];
        if (v < 0) break;
        pts.push_back(originalPosition[v]);
        if (v == b) return pts;
    }
    // b is not downstream of a on any original loop: the two are not the
    // ends of one current boundary edge. The straight segment is all there
    // is to say about them.
    return {originalPosition[a], originalPosition[b]};
}
