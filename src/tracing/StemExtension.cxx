#include "tracing/StemExtension.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>

#include "tracing/TriangleLocator.hxx"

namespace {

// Segment p0->p1 against q0->q1, half-open in the first (a meeting at p1
// belongs to the next segment along the continuation) and closed in the
// second. Returns the parameters along each.
bool segmentHit(const Point &p0, const Point &p1, const Point &q0, const Point &q1, double &t,
                double &u) {
    const Point r = p1 - p0, s = q1 - q0;
    const double den = cross2(r, s);
    if (std::fabs(den) <= 1e-12 * normP(r) * normP(s)) return false;   // parallel, to rounding
    const Point w = q0 - p0;
    t = cross2(w, s) / den;
    u = cross2(w, r) / den;
    const double eps = 1e-12;
    return t >= -eps && t < 1.0 - eps && u >= -eps && u <= 1.0 + eps;
}

// The arcs' segments bucketed on a uniform grid, for the one question the
// continuation asks after every step: which arcs does this segment cross.
class SegmentGrid {
public:
    struct Ref { int arc, seg; };

    SegmentGrid(const std::vector<QuadLayout::Arc> &arcs, double cellSize) : cell(cellSize) {
        lo = Point{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
        Point hi{-lo[0], -lo[1]};
        for (const auto &a : arcs)
            for (const Point &p : a.pts) {
                lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
                hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
            }
        if (!std::isfinite(lo[0]) || !(cell > 0.0)) return;
        nx = std::max(1, static_cast<int>((hi[0] - lo[0]) / cell) + 1);
        ny = std::max(1, static_cast<int>((hi[1] - lo[1]) / cell) + 1);
        // A grid far finer than the layout is only slower; cap the cell count.
        while (static_cast<long long>(nx) * ny > 4000000LL) {
            cell *= 2.0;
            nx = std::max(1, static_cast<int>((hi[0] - lo[0]) / cell) + 1);
            ny = std::max(1, static_cast<int>((hi[1] - lo[1]) / cell) + 1);
        }
        cells.assign(static_cast<size_t>(nx) * ny, {});
        for (size_t a = 0; a < arcs.size(); ++a)
            for (size_t k = 0; k + 1 < arcs[a].pts.size(); ++k)
                forCells(arcs[a].pts[k], arcs[a].pts[k + 1], [&](std::vector<Ref> &c) {
                    c.push_back(Ref{static_cast<int>(a), static_cast<int>(k)});
                });
    }

    std::vector<Ref> query(const Point &p0, const Point &p1) {
        std::vector<Ref> out;
        forCells(p0, p1, [&](std::vector<Ref> &c) { out.insert(out.end(), c.begin(), c.end()); });
        std::sort(out.begin(), out.end(), [](const Ref &x, const Ref &y) {
            return x.arc != y.arc ? x.arc < y.arc : x.seg < y.seg;
        });
        out.erase(std::unique(out.begin(), out.end(),
                              [](const Ref &x, const Ref &y) { return x.arc == y.arc && x.seg == y.seg; }),
                  out.end());
        return out;
    }

private:
    template <class F>
    void forCells(const Point &a, const Point &b, F &&f) {
        if (cells.empty()) return;
        const int i0 = cx(std::min(a[0], b[0])), i1 = cx(std::max(a[0], b[0]));
        const int j0 = cy(std::min(a[1], b[1])), j1 = cy(std::max(a[1], b[1]));
        for (int j = j0; j <= j1; ++j)
            for (int i = i0; i <= i1; ++i) f(cells[static_cast<size_t>(j) * nx + i]);
    }
    int cx(double x) const { return std::min(nx - 1, std::max(0, static_cast<int>((x - lo[0]) / cell))); }
    int cy(double y) const { return std::min(ny - 1, std::max(0, static_cast<int>((y - lo[1]) / cell))); }

    double cell = 1.0;
    Point lo{0.0, 0.0};
    int nx = 0, ny = 0;
    std::vector<std::vector<Ref>> cells;
};

double polylineLength(const std::vector<Point> &p) {
    double L = 0.0;
    for (size_t i = 1; i < p.size(); ++i) L += normP(p[i] - p[i - 1]);
    return L;
}

// Faces that are not disks: an arc with the same face on both sides, or a
// face visiting a node twice. Neither may be made by a continuation.
int nonDisks(const QuadLayout &L) {
    int n = 0;
    std::vector<int> faceOf(2 * L.getArcs().size(), -1);
    for (size_t f = 0; f < L.getFaces().size(); ++f)
        for (const int d : L.getFaces()[f].darts) faceOf[d] = static_cast<int>(f);
    std::vector<char> bad(L.getFaces().size(), 0);
    for (size_t a = 0; a < L.getArcs().size(); ++a)
        if (faceOf[2 * a] >= 0 && faceOf[2 * a] == faceOf[2 * a + 1]) bad[faceOf[2 * a]] = 1;
    for (size_t f = 0; f < L.getFaces().size(); ++f) {
        std::unordered_set<int> seen;
        for (const int nd : L.getFaces()[f].nodes)
            if (!seen.insert(nd).second) bad[f] = 1;
        n += bad[f];
    }
    return n;
}

// Move a polyline's last point to `to`, and the points before it by the same
// displacement falling smoothly to nothing `over` back along the line.
void easeOnto(std::vector<Point> &pts, const Point &to, double over) {
    if (pts.size() < 2) { if (!pts.empty()) pts.back() = to; return; }
    const Point d = to - pts.back();
    double s = 0.0;
    pts.back() = to;
    for (size_t i = pts.size() - 1; i-- > 1;) {
        s += normP(pts[i + 1] - pts[i]);
        if (s >= over) break;
        const double t = s / over;
        pts[i] = pts[i] + d * (1.0 - t * t * (3.0 - 2.0 * t));
    }
}

} // namespace

StemExtension::StemExtension(const QuadLayout &layout, const FieldTracer &tracer,
                             const Settings &settings)
    : layout_(layout), tracer_(&tracer), settings_(settings) {
    report_.tJunctionsBefore = layout_.getReport().tJunctions;
    report_.tJunctionsAfter = report_.tJunctionsBefore;
}

// ---------------------------------------------------------------------------
// findJunctions()
//
// Structurally, as LayoutBlocks finds them: a face that does not turn a corner
// at a node where three arcs meet, away from the boundary. The two darts that
// face runs along through the node are the bar; the third is the stem.
// ---------------------------------------------------------------------------
std::vector<StemExtension::Junction> StemExtension::findJunctions() const {
    const auto &nodes = layout_.getNodes();
    const auto &arcs = layout_.getArcs();
    std::vector<Junction> out;
    std::unordered_set<int> seen;
    for (const QuadLayout::Face &f : layout_.getFaces()) {
        const size_t k = f.darts.size();
        for (size_t i = 0; i < k; ++i) {
            if (f.isCorner[i]) continue;
            const int n = f.nodes[i];
            const auto &nd = nodes[n];
            if (nd.darts.size() != 3 || nd.kind == QuadLayout::NodeKind::Singularity ||
                nd.kind == QuadLayout::NodeKind::InterfaceNode)
                continue;
            bool boundary = false;
            for (const int d : nd.darts) if (arcs[QuadLayout::arcOfDart(d)].onBoundary) boundary = true;
            if (boundary || !seen.insert(n).second) continue;
            const int out1 = f.darts[i];
            const int out2 = f.darts[(i + k - 1) % k] ^ 1;
            Junction j;
            j.node = n;
            for (const int d : nd.darts) if (d != out1 && d != out2) j.stemDart = d;
            if (j.stemDart >= 0) out.push_back(j);
        }
    }
    return out;
}

void StemExtension::run() {
    // Each success changes the layout the next one is traced against, so the
    // junctions are found again after every one; a junction that failed once
    // is not retried, since nothing it would be traced through has changed in
    // a way that makes it likelier to succeed.
    std::unordered_set<long long> tried;
    auto key = [&](const Point &p) {
        const double h = std::max(1e-12, 1e-7 * tracer_->averageEdgeLength());
        return static_cast<long long>(std::llround(p[0] / h)) * 1000003LL +
               static_cast<long long>(std::llround(p[1] / h));
    };
    for (int guard = 0; guard < 10000; ++guard) {
        bool progressed = false;
        for (const Junction &j : findJunctions()) {
            const Point at = layout_.getNodes()[j.node].pos;
            if (!tried.insert(key(at)).second) continue;
            ++report_.attempted;
            if (extend(j)) {
                ++report_.extended;
                progressed = true;
                break;
            }
        }
        if (!progressed) break;
    }
    report_.tJunctionsAfter = layout_.getReport().tJunctions;
}

// ---------------------------------------------------------------------------
// extend()  --  one continuation
//
// Traced with the same FieldTracer the separatrices were, from the junction
// along the direction the stem arrived in, which lifts the field onto the
// branch the stem was following. After every triangle the new segments are
// tested against every arc of the layout; the first boundary arc met ends it,
// any other arc met squarely is recorded as a crossing, anything else refuses
// it. The layout is then rebuilt with the crossed arcs split and the
// continuation added, and checked.
// ---------------------------------------------------------------------------
bool StemExtension::extend(const Junction &j) {
    const auto &nodes = layout_.getNodes();
    const auto &arcs = layout_.getArcs();
    const Mesh &mesh = tracer_->getMesh();
    const double h = tracer_->averageEdgeLength();
    const double nodeTol = 1e-4 * h;

    const Point t = nodes[j.node].pos;
    const QuadLayout::Arc &stem = arcs[QuadLayout::arcOfDart(j.stemDart)];
    std::vector<Point> sp = stem.pts;
    if (j.stemDart & 1) std::reverse(sp.begin(), sp.end());   // now from the junction outwards
    Point dir{0.0, 0.0};
    for (size_t i = 1; i < sp.size() && normP(dir) < 1e-9 * h; ++i) dir = t - sp[i];
    if (normP(dir) < 1e-12 * h) return false;
    dir = dir / normP(dir);

    if (!locator_) locator_ = std::make_unique<TriangleLocator>(mesh);
    // Just ahead of the junction, so that a junction on a mesh edge is placed
    // in the triangle the continuation leaves into and not the one behind it.
    const int tri = locator_->locate(t + dir * (1e-6 * h));
    if (tri < 0) { ++report_.noBoundary; return false; }

    SegmentGrid grid(arcs, 2.0 * h);

    // The other T-junctions, by node, with the direction each one's stem
    // leaves it in: the ones this continuation may end on (see below).
    std::vector<Point> stemOut(nodes.size(), Point{0.0, 0.0});
    std::vector<char> joinable(nodes.size(), 0);
    if (settings_.joinRadius > 0.0) {
        for (const Junction &o : findJunctions()) {
            if (o.node == j.node) continue;
            std::vector<Point> op = arcs[QuadLayout::arcOfDart(o.stemDart)].pts;
            if (o.stemDart & 1) std::reverse(op.begin(), op.end());
            const Point on = nodes[o.node].pos;
            Point od{0.0, 0.0};
            for (size_t i = 1; i < op.size() && normP(od) < 1e-9 * h; ++i) od = op[i] - on;
            if (normP(od) < 1e-12 * h) continue;
            stemOut[o.node] = od / normP(od);
            joinable[o.node] = 1;
        }
    }
    int joinedAt = -1;

    struct Crossing {
        int arc = -1;
        double param = 0.0;      // seg + u along the arc
        int pathSeg = 0;         // the continuation's segment it is on
        double along = 0.0;      // t along that segment
        Point at{0.0, 0.0};
        bool boundary = false;
    };
    std::vector<Crossing> xs;
    std::vector<Point> path{t};

    Walker w = tracer_->startAt(tri, t, std::atan2(dir[1], dir[0]));
    bool done = false;
    for (int step = 0; step < settings_.maxSteps && !done; ++step) {
        std::vector<TracePoint> seg;
        const FieldTracer::Status st = tracer_->advance(w, seg);
        for (const TracePoint &tp : seg) {
            const Point p0 = path.back(), p1 = tp.global_pos;
            if (normP(p1 - p0) < 1e-14 * h) continue;
            const int pathSeg = static_cast<int>(path.size()) - 1;
            std::vector<Crossing> hits;
            for (const SegmentGrid::Ref &r : grid.query(p0, p1)) {
                const auto &a = arcs[r.arc];
                double s = 0.0, u = 0.0;
                if (!segmentHit(p0, p1, a.pts[r.seg], a.pts[r.seg + 1], s, u)) continue;
                const Point at = p0 + (p1 - p0) * s;
                // Where it starts: on the bar and the stem, at the junction.
                if (normP(at - t) < nodeTol && (a.a == j.node || a.b == j.node)) continue;
                Crossing c;
                c.arc = r.arc;
                c.param = r.seg + std::min(1.0, std::max(0.0, u));
                c.pathSeg = pathSeg;
                c.along = s;
                c.at = at;
                c.boundary = a.onBoundary;
                hits.push_back(c);
            }
            std::sort(hits.begin(), hits.end(),
                      [](const Crossing &x, const Crossing &y) { return x.along < y.along; });
            // A crossing through a vertex of an arc is found on both segments
            // meeting there; it is one crossing.
            hits.erase(std::unique(hits.begin(), hits.end(),
                                   [&](const Crossing &x, const Crossing &y) {
                                       return x.arc == y.arc && std::fabs(x.param - y.param) < 1e-9;
                                   }),
                       hits.end());
            for (const Crossing &c : hits) {
                const auto &a = arcs[c.arc];
                // Nose to nose with another T-junction: the continuation
                // crosses the side that one stopped on, within joinRadius of
                // it, heading the way its stem leaves it. The two stems are
                // one streamline the discretisation split -- one stopped on a
                // curve the other then stopped on a fraction of an element
                // further along -- so the continuation ends there, on that
                // junction, and both become crossings. Carrying on instead
                // would run alongside the other stem (refused as tangential)
                // or through its node (refused as through a node), and the
                // pair would stay.
                if (!c.boundary && settings_.joinRadius > 0.0) {
                    int J = -1;
                    double best = settings_.joinRadius * h;
                    for (const int e : {a.a, a.b}) {
                        if (e < 0 || !joinable[e]) continue;
                        const double dd = normP(c.at - nodes[e].pos);
                        const Point d = p1 - p0;
                        const double cosd = dotP(stemOut[e], d) / std::max(1e-300, normP(d));
                        if (dd < best && cosd > std::cos(settings_.minCrossingAngle)) { best = dd; J = e; }
                    }
                    if (J >= 0) {
                        joinedAt = J;
                        xs.push_back(c);   // where it crossed; eased onto J below
                        done = true;
                        break;
                    }
                }
                if (normP(c.at - nodes[a.a].pos) < nodeTol || normP(c.at - nodes[a.b].pos) < nodeTol) {
                    ++report_.throughNode;
                    return false;
                }
                const int k = std::min(static_cast<int>(a.pts.size()) - 2,
                                       static_cast<int>(std::floor(c.param)));
                const Point e = a.pts[k + 1] - a.pts[k];
                const Point d = p1 - p0;
                const double cosw = std::fabs(dotP(e, d)) / std::max(1e-300, normP(e) * normP(d));
                if (!c.boundary && std::acos(std::min(1.0, cosw)) < settings_.minCrossingAngle) {
                    ++report_.tangential;
                    return false;
                }
                xs.push_back(c);
                if (c.boundary) { done = true; break; }
                if (static_cast<int>(xs.size()) > settings_.maxCrossings) {
                    ++report_.tooManyCrossings;
                    return false;
                }
            }
            if (done) break;
            path.push_back(p1);
        }
        if (done) break;
        if (st == FieldTracer::Status::Boundary) {
            // Left the mesh through a boundary edge, at the last point traced:
            // it lies on a boundary arc, which is split there.
            const Point end = path.back();
            int best = -1, bestSeg = 0;
            double bestD = std::numeric_limits<double>::max(), bestU = 0.0;
            for (size_t a = 0; a < arcs.size(); ++a) {
                if (!arcs[a].onBoundary) continue;
                for (size_t k = 0; k + 1 < arcs[a].pts.size(); ++k) {
                    const Point q0 = arcs[a].pts[k], e = arcs[a].pts[k + 1] - q0;
                    const double ee = dotP(e, e);
                    const double u = ee > 0.0 ? std::min(1.0, std::max(0.0, dotP(end - q0, e) / ee)) : 0.0;
                    const double dd = normP(end - (q0 + e * u));
                    if (dd < bestD) { bestD = dd; best = static_cast<int>(a); bestSeg = static_cast<int>(k); bestU = u; }
                }
            }
            if (best < 0 || bestD > 1e-6 * h) { ++report_.noBoundary; return false; }
            const auto &a = arcs[best];
            if (normP(end - nodes[a.a].pos) < nodeTol || normP(end - nodes[a.b].pos) < nodeTol) {
                ++report_.throughNode;
                return false;
            }
            Crossing c;
            c.arc = best;
            c.param = bestSeg + bestU;
            c.pathSeg = static_cast<int>(path.size()) - 1;
            c.along = 0.0;
            c.at = end;
            c.boundary = true;
            path.pop_back();   // the end point is the node itself
            c.pathSeg = static_cast<int>(path.size()) - 1;
            c.along = 1.0;
            xs.push_back(c);
            done = true;
            break;
        }
        if (st != FieldTracer::Status::Ok) { ++report_.noBoundary; return false; }
    }
    if (!done || xs.empty() || (!xs.back().boundary && joinedAt < 0)) { ++report_.noBoundary; return false; }

    // --- the new layout ---------------------------------------------------
    std::vector<QuadLayout::Node> newNodes = nodes;
    for (auto &n : newNodes) { n.darts.clear(); n.angles.clear(); }
    newNodes[j.node].kind = QuadLayout::NodeKind::Crossing;
    std::vector<int> nodeOf(xs.size());
    for (size_t i = 0; i < xs.size(); ++i) {
        if (joinedAt >= 0 && i + 1 == xs.size()) {
            nodeOf[i] = joinedAt;
            newNodes[joinedAt].kind = QuadLayout::NodeKind::Crossing;
            continue;
        }
        QuadLayout::Node n;
        n.pos = xs[i].at;
        n.kind = xs[i].boundary ? QuadLayout::NodeKind::BoundaryHit : QuadLayout::NodeKind::Crossing;
        nodeOf[i] = static_cast<int>(newNodes.size());
        newNodes.push_back(n);
    }

    auto makeArc = [&](std::vector<Point> pts, int a, int b, int separatrix, bool boundary,
                       bool interface) {
        pts.front() = newNodes[a].pos;
        pts.back() = newNodes[b].pos;
        std::vector<Point> clean;
        for (const Point &p : pts)
            if (clean.empty() || normP(p - clean.back()) > 1e-12 * h) clean.push_back(p);
        if (clean.size() < 2) clean.push_back(newNodes[b].pos);
        clean.back() = newNodes[b].pos;
        QuadLayout::Arc arc;
        arc.pts = std::move(clean);
        arc.a = a;
        arc.b = b;
        arc.separatrix = separatrix;
        arc.onBoundary = boundary;
        arc.onInterface = interface;
        arc.length = polylineLength(arc.pts);
        return arc;
    };

    // Every crossed arc cut at its crossings, in order along it.
    std::vector<std::vector<std::pair<double, int>>> cuts(arcs.size());
    for (size_t i = 0; i < xs.size(); ++i)
        if (!(joinedAt >= 0 && i + 1 == xs.size())) cuts[xs[i].arc].push_back({xs[i].param, nodeOf[i]});
    std::vector<QuadLayout::Arc> newArcs;
    for (size_t a = 0; a < arcs.size(); ++a) {
        if (cuts[a].empty()) { newArcs.push_back(arcs[a]); continue; }
        auto c = cuts[a];
        std::sort(c.begin(), c.end());
        const auto &pts = arcs[a].pts;
        double prevP = 0.0;
        int prevNode = arcs[a].a;
        auto pointAt = [&](double p) {
            const int k = std::min(static_cast<int>(pts.size()) - 2, static_cast<int>(std::floor(p)));
            const double u = p - k;
            return pts[k] * (1.0 - u) + pts[k + 1] * u;
        };
        c.push_back({static_cast<double>(pts.size() - 1), arcs[a].b});
        for (const auto &[p, nd] : c) {
            std::vector<Point> piece{pointAt(prevP)};
            for (int k = static_cast<int>(std::floor(prevP)) + 1; k <= static_cast<int>(std::ceil(p)); ++k)
                if (k > prevP + 1e-12 && k < p - 1e-12) piece.push_back(pts[k]);
            piece.push_back(pointAt(p));
            newArcs.push_back(makeArc(piece, prevNode, nd, arcs[a].separatrix, arcs[a].onBoundary,
                                      arcs[a].onInterface));
            prevP = p;
            prevNode = nd;
        }
    }

    // The continuation, cut at the same crossings.
    {
        int prevNode = j.node;
        std::vector<Point> piece{t};
        size_t xi = 0;
        for (size_t s = 0; s < path.size(); ++s) {
            if (s > 0) piece.push_back(path[s]);
            while (xi < xs.size() && xs[xi].pathSeg == static_cast<int>(s)) {
                piece.push_back(xs[xi].at);
                if (joinedAt >= 0 && xi + 1 == xs.size()) {
                    // Bent onto the junction over the last stretch rather than
                    // in the last segment: the continuation passed it at up to
                    // joinRadius to one side.
                    easeOnto(piece, nodes[joinedAt].pos, 4.0 * settings_.joinRadius * h);
                }
                newArcs.push_back(makeArc(piece, prevNode, nodeOf[xi], -1, false, false));
                prevNode = nodeOf[xi];
                piece.assign(1, xs[xi].at);
                ++xi;
            }
        }
        if (xi != xs.size()) { ++report_.invalid; return false; }
    }

    QuadLayout next = layout_;
    next.rebuild(std::move(newNodes), std::move(newArcs));
    const auto &rb = layout_.getReport();
    const auto &ra = next.getReport();
    double lenBefore = 0.0, lenAfter = 0.0;
    for (const auto &a : layout_.getArcs()) if (a.onInterface) lenBefore += a.length;
    for (const auto &a : next.getArcs()) if (a.onInterface) lenAfter += a.length;
    const bool ok = ra.arcCrossings == 0 && ra.danglingEnds == 0 &&
                    std::fabs(lenAfter - lenBefore) <= 1e-9 * std::max(1.0, lenBefore) &&
                    std::fabs(ra.totalArea - rb.totalArea) <= 1e-9 * std::max(1.0, rb.totalArea) &&
                    ra.tJunctions < rb.tJunctions && ra.badFaces <= rb.badFaces &&
                    nonDisks(next) <= nonDisks(layout_);
    if (!ok) { ++report_.invalid; return false; }
    layout_ = std::move(next);
    if (joinedAt >= 0) ++report_.joined;
    return true;
}
