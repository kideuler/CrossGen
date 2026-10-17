#include "SKELETON/MedialAxisTMesh.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <unordered_set>

namespace {

// Polyline helpers for the walk and the template construction.

// Append, dropping a duplicated joint point.
void appendPts(std::vector<Point> &dst, const std::vector<Point> &src) {
    for (const Point &p : src) {
        if (!dst.empty() && normP(p - dst.back()) < 1e-14) continue;
        dst.push_back(p);
    }
}

std::vector<double> cumulativeArc(const std::vector<Point> &pts) {
    std::vector<double> arc(pts.size(), 0.0);
    for (size_t i = 1; i < pts.size(); ++i) {
        arc[i] = arc[i - 1] + normP(pts[i] - pts[i - 1]);
    }
    return arc;
}

// The point at arc-length fraction t of a polyline whose cumulative arc
// length is `arc`.
Point pointAtFraction(const std::vector<Point> &pts, const std::vector<double> &arc,
                      double t) {
    if (pts.size() == 1 || !(arc.back() > 0.0)) return pts.front();
    const double target = std::clamp(t, 0.0, 1.0) * arc.back();
    for (size_t i = 1; i < pts.size(); ++i) {
        if (arc[i] >= target) {
            const double seg = arc[i] - arc[i - 1];
            const double f = seg > 0.0 ? (target - arc[i - 1]) / seg : 0.0;
            return pts[i - 1] + (pts[i] - pts[i - 1]) * f;
        }
    }
    return pts.back();
}

// Split a polyline in two at arc-length fraction t; both pieces share the
// cut point.
void splitAtFraction(const std::vector<Point> &pts, double t,
                     std::vector<Point> &head, std::vector<Point> &tail) {
    const std::vector<double> arc = cumulativeArc(pts);
    const double target = std::clamp(t, 0.0, 1.0) * arc.back();
    head = {pts.front()};
    size_t i = 1;
    while (i < pts.size() && arc[i] < target) head.push_back(pts[i++]);
    const Point cut = pointAtFraction(pts, arc, t);
    appendPts(head, {cut});
    tail = {cut};
    for (; i < pts.size(); ++i) appendPts(tail, {pts[i]});
}

double distToSegment(const Point &p, const Point &a, const Point &b) {
    const Point e = b - a;
    const double len2 = dotP(e, e);
    const double t = len2 > 0.0 ? std::clamp(dotP(p - a, e) / len2, 0.0, 1.0) : 0.0;
    const Point q = a + e * t;
    return normP(p - q);
}

// Arc-length parameter of p along the sub-polyline midPrev -> v -> midNext,
// deciding which leg p sits on by proximity.
double subPolylineParam(const Point &p, const Point &midPrev, const Point &v,
                        const Point &midNext) {
    const double lenPrev = normP(v - midPrev);
    if (distToSegment(p, midPrev, v) <= distToSegment(p, v, midNext)) {
        return normP(p - midPrev);
    }
    return lenPrev + normP(p - v);
}

double signedArea(const std::vector<Point> &pts) {
    double a = 0.0;
    for (size_t i = 0; i < pts.size(); ++i) {
        a += cross2(pts[i], pts[(i + 1) % pts.size()]);
    }
    return 0.5 * a;
}

// Unsigned angle between two vectors, in [0, pi].
double angleBetween(const Point &a, const Point &b) {
    return std::atan2(std::abs(cross2(a, b)), dotP(a, b));
}

long long pairKey(int a, int b) {
    return (static_cast<long long>(std::min(a, b)) << 32) |
           static_cast<long long>(std::max(a, b));
}

MedialAxisTMesh::Options withSize(double targetSize) {
    MedialAxisTMesh::Options o;
    o.targetSize = targetSize;
    return o;
}

} // namespace

// ─── Construction ───────────────────────────────────────────────────────────

MedialAxisTMesh::MedialAxisTMesh(const MedialAxisMap &map, double targetSize)
    : MedialAxisTMesh(map, withSize(targetSize)) {}

MedialAxisTMesh::MedialAxisTMesh(const MedialAxisMap &map, const Options &opts)
    : axis(map.axis), opts_(opts) {
    if (!axis) {
        throw std::runtime_error("MedialAxisMap carries no axis.");
    }
    targetSize_ = opts_.targetSize > 0.0 ? opts_.targetSize : axis->defaultTargetSize();

    // Branches are the unit of downsampling, so make sure they exist.
    if (axis->polyLines.empty()) axis->createPolylines();
    branchScale_.assign(axis->polyLines.size(), 1.0);

    const int nv = static_cast<int>(axis->mesh->vertices.size());
    originalPrev_.assign(nv, -1);
    for (int v = 0; v < nv; ++v) {
        const int n = axis->originalNext[v];
        if (n >= 0) originalPrev_[n] = v;
    }

    // Step 8.6: a face that is not a usable cell means h was too coarse
    // there -- Sec. 5.3 finds invalid embeddings only when it is -- so h is
    // halved on the branches that produced one and everything from the
    // downsampling on is rebuilt. A round is kept only if it leaves fewer bad
    // faces than the best so far: not every failure is one of h -- a purple
    // cap on a lobe whose neck narrows away from it folds at any size, as
    // does the collocated fan at a concave corner -- and refining around one
    // of those only multiplies the faces that fail.
    build(map);
    std::vector<int> bad = checkGeometry();
    std::vector<double> bestScale = branchScale_;
    int best = stats_.badFaces;
    int refinements = 0;
    while (!bad.empty() && refinements < opts_.maxRefinements) {
        for (int b : bad) branchScale_[b] *= 0.5;
        build(map);
        bad = checkGeometry();
        if (stats_.badFaces >= best) break;
        best = stats_.badFaces;
        bestScale = branchScale_;
        ++refinements;
    }
    if (branchScale_ != bestScale) {
        branchScale_ = bestScale;
        build(map);
        checkGeometry();
    }
    stats_.refinements = refinements;
}

void MedialAxisTMesh::build(const MedialAxisMap &map) {
    zones.clear();
    edges.clear();
    blocks.clear();
    loops_.clear();
    stats_ = TMeshStats{};

    selectNodes();
    gatherLoops(map);
    emitZones();
    buildEdges();
    classify();
    buildBlocks();

    stats_.zones = static_cast<int>(zones.size());
    stats_.coarseEdges = static_cast<int>(edges.size());
    stats_.blocks = static_cast<int>(blocks.size());
}

int MedialAxisTMesh::nodeCount() const {
    return static_cast<int>(axis->medialVertices.size() + inserted_.size());
}

Point MedialAxisTMesh::nodeCoord(int node) const {
    const int nAxis = static_cast<int>(axis->medialVertices.size());
    return node < nAxis ? axis->medialVertices[node].coord : inserted_[node - nAxis].coord;
}

double MedialAxisTMesh::nodeRadius(int node) const {
    const int nAxis = static_cast<int>(axis->medialVertices.size());
    return node < nAxis ? axis->medialVertices[node].radius : inserted_[node - nAxis].radius;
}

int MedialAxisTMesh::nodeDegree(int node) const {
    const int nAxis = static_cast<int>(axis->medialVertices.size());
    return node < nAxis ? axis->medialVertices[node].degree : 2;
}

bool MedialAxisTMesh::isKept(int node) const {
    const int nAxis = static_cast<int>(axis->medialVertices.size());
    if (node < 0) return false;
    if (node >= nAxis) return node < nodeCount();
    return node < static_cast<int>(kept_.size()) && kept_[node] != 0;
}

int MedialAxisTMesh::contactCount(int node) const {
    return node >= 0 && node < static_cast<int>(runCount_.size()) ? runCount_[node] : 0;
}

std::vector<std::array<Point, 2>> MedialAxisTMesh::radii() const {
    std::vector<std::array<Point, 2>> out;
    out.reserve(spokes_.size());
    for (const Spoke &sp : spokes_) out.push_back({nodeCoord(sp.node), sp.foot});
    return out;
}

double MedialAxisTMesh::originalTurn(int v) const {
    const int p = originalPrev_[v];
    const int n = axis->originalNext[v];
    if (p < 0 || n < 0) return 0.0;
    const Point in  = axis->originalPosition[v] - axis->originalPosition[p];
    const Point out = axis->originalPosition[n] - axis->originalPosition[v];
    return std::atan2(cross2(in, out), dotP(in, out));
}

// ─── Step 6.1: downsampling ─────────────────────────────────────────────────

void MedialAxisTMesh::selectNodes() {
    const int nAxis = static_cast<int>(axis->medialVertices.size());
    kept_.assign(nAxis, 0);
    inserted_.clear();
    insertedOn_.clear();
    coarseBranch_.clear();

    const double tiny = 1e-12 * axis->boundingBoxDiagonal();
    for (size_t b = 0; b < axis->polyLines.size(); ++b) {
        const std::vector<int> &pl = axis->polyLines[b];
        const std::vector<int> &pe = axis->polyLineEdges[b];
        // A branch of one vertex is a medial vertex the inside filter left
        // isolated: its circle is not inside the domain, so it is not on the
        // axis and must not cut the boundary. Its passages still lie in the
        // walk, inside whichever zone they fall in.
        if (pl.size() < 2) continue;
        // Branch ends are always kept: junctions and endpoints are where the
        // coarse axis has to have its vertices.
        kept_[pl.front()] = 1;
        kept_[pl.back()] = 1;

        std::vector<double> arc(pl.size(), 0.0);
        for (size_t i = 1; i < pl.size(); ++i) {
            arc[i] = arc[i - 1] + normP(axis->medialVertices[pl[i]].coord -
                                        axis->medialVertices[pl[i - 1]].coord);
        }
        const double L = arc.back();

        // n = max(n_min, round(L / h)). n_min [R] keeps the endpoint
        // templates of a branch from overlapping: a branch ending in an
        // endpoint gets two coarse edges, one ending in two endpoints or
        // closing on itself three.
        const bool cyclic = axis->polyLineIsCycle[b];
        const bool endA = !cyclic && axis->medialVertices[pl.front()].degree == 1;
        const bool endB = !cyclic && axis->medialVertices[pl.back()].degree == 1;
        const int nMin = (cyclic || (endA && endB)) ? 3 : ((endA || endB) ? 2 : 1);
        const double hb = targetSize_ * branchScale_[b];
        int n = std::max(nMin, static_cast<int>(std::lround(L / hb)));
        // A branch of no length has nowhere to put samples apart.
        if (!(L > tiny)) n = 1;
        const double spacing = L / n;

        std::vector<int> nodes{pl.front()};
        size_t prev = 0;
        for (int k = 1; k < n; ++k) {
            const double target = spacing * k;
            // The regular vertex closest to k L / n, beyond the last pick and
            // within half a spacing of the target: anything farther is
            // closer to a neighbouring target, and the sampling would bunch.
            size_t best = 0;
            double bestD = 0.5 * spacing;
            for (size_t i = prev + 1; i + 1 < pl.size(); ++i) {
                const double d = std::abs(arc[i] - target);
                if (d <= bestD) {
                    bestD = d;
                    best = i;
                }
                if (arc[i] > target + 0.5 * spacing) break;
            }
            if (best > 0) {
                kept_[pl[best]] = 1;
                nodes.push_back(pl[best]);
                prev = best;
                continue;
            }

            // Too few fine vertices: insert a sample on the fine edge that
            // holds the target. It sits on the Voronoi edge between the two
            // samples that edge is dual to, which gives its radius and, in
            // the walk, its preimage through each of their cells.
            size_t i = prev;
            while (i + 2 < pl.size() && arc[i + 1] < target) ++i;
            const double seg = arc[i + 1] - arc[i];
            const double f = seg > 0.0 ? std::clamp((target - arc[i]) / seg, 0.0, 1.0) : 0.0;
            const Point &a = axis->medialVertices[pl[i]].coord;
            const Point &c = axis->medialVertices[pl[i + 1]].coord;
            const Point q = a + (c - a) * f;
            const Edge &dual = axis->medialEdgeDual[pe[i]];
            InsertedNode node;
            node.coord = q;
            node.radius = normP(q - axis->mesh->vertices[dual[0]]);
            const int id = nAxis + static_cast<int>(inserted_.size());
            inserted_.push_back(node);
            insertedOn_[MeshEdgeKey(dual[0], dual[1])].push_back(id);
            nodes.push_back(id);
            prev = i;
        }
        nodes.push_back(pl.back());
        for (size_t k = 0; k + 1 < nodes.size(); ++k) {
            coarseBranch_[pairKey(nodes[k], nodes[k + 1])] = static_cast<int>(b);
        }
    }

    stats_.keptVertices = 0;
    for (char k : kept_) stats_.keptVertices += k;
    stats_.insertedVertices = static_cast<int>(inserted_.size());
}

// ─── The walk ───────────────────────────────────────────────────────────────

// The passages a fan contributes, in fan order: its chain vertices 0..k-2
// (the last belongs to the next fan, which starts with it) with the samples
// inserted on the medial edge between each consecutive pair. `feet` gets the
// preimage of each on the fan's sub-polyline. A collapsed fan contributes
// its one vertex at the interval entry, so consecutive collapsed fans build
// up the run that is a polar section.
std::vector<int> MedialAxisTMesh::fanSequenceNodes(const MedialAxisMap &map, int v,
                                                   std::vector<Point> &feet) const {
    feet.clear();
    const MedialFan *fan = map.fanOf(v);
    if (!fan || fan->chain.empty()) return {};
    if (fan->collapsed) {
        feet.push_back(fan->midPrev);
        return {fan->chain[0]};
    }

    std::vector<int> nodes;
    const int k = static_cast<int>(fan->chain.size());
    for (int i = 0; i + 1 < k; ++i) {
        nodes.push_back(fan->chain[i]);
        feet.push_back(fan->footpoints[i]);
        if (i >= static_cast<int>(fan->crossed.size())) continue;
        const auto it = insertedOn_.find(MeshEdgeKey(v, fan->crossed[i]));
        if (it == insertedOn_.end()) continue;
        std::vector<int> ids = it->second;
        const Point from = nodeCoord(fan->chain[i]);
        std::sort(ids.begin(), ids.end(), [&](int x, int y) {
            return normP(nodeCoord(x) - from) < normP(nodeCoord(y) - from);
        });
        for (int id : ids) {
            bool missed = false;
            nodes.push_back(id);
            feet.push_back(map.preimage(v, nodeCoord(id), missed));
        }
    }
    return nodes;
}

void MedialAxisTMesh::gatherLoops(const MedialAxisMap &map) {
    const std::vector<Point> &x = axis->mesh->vertices;
    const std::vector<Point> &x0 = axis->originalPosition;
    const int nv = static_cast<int>(x.size());
    std::vector<char> visited(nv, 0);
    runCount_.assign(nodeCount(), 0);

    // A passage as a parameter on a current boundary edge: t in [0, 1/2] on
    // the half beside the edge's start, [1/2, 1] beside its end.
    struct Passage {
        int node;
        double t;
    };

    for (int start = 0; start < nv; ++start) {
        if (visited[start] || axis->boundaryLinks[start].nextBoundaryVertex < 0) continue;
        ++stats_.loops;

        std::vector<int> loopVerts;
        int v = start;
        int guard = nv + 1;
        bool closed = false;
        while (guard-- > 0) {
            visited[v] = 1;
            loopVerts.push_back(v);
            v = axis->boundaryLinks[v].nextBoundaryVertex;
            if (v == start) { closed = true; break; }
            // Broken or pinched links: the loop never closes, and zones cut
            // from it would not either.
            if (v < 0 || visited[v]) break;
        }
        if (!closed) { ++stats_.skippedLoops; continue; }
        const int n = static_cast<int>(loopVerts.size());

        // Each fan's passages split at its vertex: those before it lie on
        // the edge coming in, those after it on the edge going out.
        // Monotonicity orders them; the clamp only absorbs roundoff.
        std::vector<std::vector<Passage>> before(n), after(n);
        for (int i = 0; i < n; ++i) {
            const int u = loopVerts[i];
            const MedialFan *fan = map.fanOf(u);
            if (!fan) continue;
            std::vector<Point> feet;
            const std::vector<int> nodes = fanSequenceNodes(map, u, feet);
            const double vparam = normP(x[u] - fan->midPrev);
            const double nextLen = normP(fan->midNext - x[u]);
            double last = 0.0;
            for (size_t j = 0; j < nodes.size(); ++j) {
                double p = fan->collapsed
                               ? 0.0
                               : subPolylineParam(feet[j], fan->midPrev, x[u], fan->midNext);
                p = std::max(p, last);
                last = p;
                if (fan->collapsed || p <= vparam) {
                    before[i].push_back({nodes[j], 0.5 + 0.5 * (vparam > 0.0 ? p / vparam : 1.0)});
                } else {
                    after[i].push_back(
                        {nodes[j], 0.5 * (nextLen > 0.0 ? (p - vparam) / nextLen : 0.0)});
                }
            }
        }

        // The events of each current edge (u, w): u itself, then the
        // passages after u and before w, interleaved with whatever input
        // vertices Step 4.2 removed from between them -- every position
        // read off the original polyline the edge stands for, at the same
        // arc-length fraction.
        LoopWalk loop;
        for (int i = 0; i < n; ++i) {
            const int u = loopVerts[i];
            const int w = loopVerts[(i + 1) % n];
            const std::vector<Point> orig = axis->originalPolyline(u, w);
            const std::vector<double> arc = cumulativeArc(orig);
            loop.events.push_back({x0[u], -1, u});

            std::vector<Passage> on = after[i];
            const std::vector<Passage> &tail = before[(i + 1) % n];
            on.insert(on.end(), tail.begin(), tail.end());
            for (size_t j = 1; j < on.size(); ++j) on[j].t = std::max(on[j].t, on[j - 1].t);

            size_t a = 0;
            size_t r = 1;   // removed vertices are orig[1 .. size - 2]
            while (a < on.size() || r + 1 < orig.size()) {
                const bool takeRemoved =
                    r + 1 < orig.size() &&
                    (a >= on.size() || arc[r] < on[a].t * arc.back());
                if (takeRemoved) {
                    loop.events.push_back({orig[r], -1, -1});
                    ++r;
                } else {
                    loop.events.push_back({pointAtFraction(orig, arc, on[a].t), on[a].node, -1});
                    ++a;
                }
            }
        }

        // Rotate the stream so it starts at a run boundary -- the first medial
        // event whose node differs from the previous medial event's -- so
        // runs never straddle the wrap-around.
        const int ne = static_cast<int>(loop.events.size());
        int lastNode = -1;
        for (int i = ne - 1; i >= 0; --i) {
            if (loop.events[i].node >= 0) { lastNode = loop.events[i].node; break; }
        }
        if (lastNode < 0) {
            // A loop no fan reached at all.
            ++stats_.skippedLoops;
            continue;
        }
        int rot = -1;
        int prevNode = lastNode;
        for (int i = 0; i < ne; ++i) {
            if (loop.events[i].node < 0) continue;
            if (loop.events[i].node != prevNode) { rot = i; break; }
            prevNode = loop.events[i].node;
        }

        if (rot < 0) {
            // Every passage is of the same node: the whole loop is one
            // contact (a domain reduced to a single inscribed circle).
            int first = -1, lastEv = -1;
            for (int i = 0; i < ne; ++i) {
                if (loop.events[i].node < 0) continue;
                if (first < 0) first = i;
                lastEv = i;
            }
            loop.runs.push_back({lastNode, first, lastEv});
        } else {
            std::rotate(loop.events.begin(), loop.events.begin() + rot, loop.events.end());
            for (int i = 0; i < ne; ++i) {
                const int m = loop.events[i].node;
                if (m < 0) continue;
                if (!loop.runs.empty() && loop.runs.back().m == m) {
                    loop.runs.back().last = i;
                } else {
                    loop.runs.push_back({m, i, i});
                }
            }
        }
        for (const MedialRun &run : loop.runs) ++runCount_[run.m];
        loops_.push_back(std::move(loop));
    }
}

// ─── Step 6.2/6.3: zones ────────────────────────────────────────────────────

void MedialAxisTMesh::emitZones() {
    // A place on a loop's event stream: event i itself, or a point strictly
    // between events i and i + 1.
    struct Loc {
        int i = -1;
        Point p{0.0, 0.0};
        bool exact = true;
    };
    // Where one zone hands over to the next. A single preimage hands over at
    // a point (enter == exit). An endpoint's polar section hands over across
    // its span, and the piece between enter and exit is its cap zone.
    struct Cut {
        int run = -1;
        Loc enter, exit;
        bool cap = false;
        int enterSpoke = -1, exitSpoke = -1;
    };
    spokes_.clear();

    for (const LoopWalk &loop : loops_) {
        const std::vector<WalkEvent> &events = loop.events;
        const std::vector<MedialRun> &runs = loop.runs;
        const int n = static_cast<int>(events.size());
        const int nr = static_cast<int>(runs.size());
        if (nr == 0) continue;

        std::vector<double> s(n, 0.0);
        for (int i = 1; i < n; ++i) s[i] = s[i - 1] + normP(events[i].pos - events[i - 1].pos);

        auto at = [&](int i) { return Loc{i, events[i].pos, true}; };

        std::vector<Cut> cuts;
        for (int r = 0; r < nr; ++r) {
            const MedialRun &run = runs[r];
            if (!isKept(run.m)) continue;
            if (run.first == run.last) {
                cuts.push_back({r, at(run.first), at(run.first), false});
            } else if (nodeDegree(run.m) == 1) {
                cuts.push_back({r, at(run.first), at(run.last), true});
            } else {
                // A polar section away from an endpoint, left by a pruned
                // branch or an arc 4.2 could not flatten. Its arc-length
                // midpoint is the single preimage of both subdomains beside
                // it (Step 5.3, [R]). Runs do not straddle the wrap, so the
                // span is first .. last in index order.
                const double mid = 0.5 * (s[run.first] + s[run.last]);
                int j = run.first;
                while (j + 1 < run.last && s[j + 1] <= mid) ++j;
                Loc loc;
                if (s[j] == mid) {
                    loc = at(j);
                } else {
                    const double seg = s[j + 1] - s[j];
                    const double f = seg > 0.0 ? (mid - s[j]) / seg : 0.0;
                    loc = Loc{j, events[j].pos + (events[j + 1].pos - events[j].pos) * f, false};
                }
                cuts.push_back({r, loc, loc, false});
                ++stats_.midpointCuts;
            }
        }
        if (cuts.empty()) {
            // No kept vertex on this loop -- degenerate sampling. Nothing
            // sound to cut, so report rather than fabricate.
            ++stats_.skippedLoops;
            continue;
        }
        // One radius per preimage: a cut at a point is one spoke, shared by
        // the zones on either side of it; a cap's two ends are two.
        for (Cut &c : cuts) {
            const int m = runs[c.run].m;
            c.enterSpoke = static_cast<int>(spokes_.size());
            spokes_.push_back({m, c.enter.p});
            c.exitSpoke = c.enterSpoke;
            if (c.cap) {
                c.exitSpoke = static_cast<int>(spokes_.size());
                spokes_.push_back({m, c.exit.p});
            }
        }

        // The stretch of the stream from a forward to b, both ends included:
        // a's point, every event after it up to b's index, and b's own point
        // when it falls between events. Two cuts never share a segment --
        // the next cut belongs to a later run -- so the count below is the
        // whole story, and a == b (a loop with one cut) walks it all.
        auto slice = [&](const Loc &a, const Loc &b, MedialZone &zone) {
            zone.boundarySide = {a.p};
            zone.boundaryVertexAt = {a.exact ? events[a.i].boundaryVertex : -1};
            const int first = (a.i + 1) % n;
            const int count = (b.i - first + n) % n + 1;
            for (int k = 0, j = first; k < count; ++k, j = (j + 1) % n) {
                zone.boundarySide.push_back(events[j].pos);
                zone.boundaryVertexAt.push_back(events[j].boundaryVertex);
            }
            if (!b.exact) {
                zone.boundarySide.push_back(b.p);
                zone.boundaryVertexAt.push_back(-1);
            }
        };

        for (const Cut &c : cuts) {
            if (!c.cap) continue;
            MedialZone zone;
            slice(c.enter, c.exit, zone);
            zone.spokes = {c.enterSpoke, c.exitSpoke};
            zone.chain = {runs[c.run].m};
            zone.cap = true;
            zone.capVertex = runs[c.run].m;
            zones.push_back(std::move(zone));
            ++stats_.capZones;
        }
        // A loop that is a single contact has nothing to cut between: its
        // cap, if its node is an endpoint, is all of it.
        if (nr == 1) continue;

        // The interval between each pair of consecutive cuts, with the nodes
        // passed in between as its chain.
        const int nc = static_cast<int>(cuts.size());
        for (int j = 0; j < nc; ++j) {
            const Cut &a = cuts[j];
            const Cut &b = cuts[(j + 1) % nc];
            MedialZone zone;
            slice(a.exit, b.enter, zone);
            zone.spokes = {a.exitSpoke, b.enterSpoke};
            zone.chain.push_back(runs[a.run].m);
            for (int ri = (a.run + 1) % nr; ri != b.run; ri = (ri + 1) % nr) {
                zone.chain.push_back(runs[ri].m);
            }
            zone.chain.push_back(runs[b.run].m);
            zones.push_back(std::move(zone));
        }
    }
}

// ─── Edges and Step 7 coloring ──────────────────────────────────────────────

void MedialAxisTMesh::buildEdges() {
    // The two zones flanking a chain see it in opposite directions, so the
    // canonical signature is the lexicographically smaller of the chain and
    // its reverse. Two parallel branches with identical signatures (a double
    // edge between two junctions) are paired in encounter order.
    std::map<std::vector<int>, std::vector<int>> bySignature;
    for (int z = 0; z < static_cast<int>(zones.size()); ++z) {
        if (zones[z].cap) continue;
        std::vector<int> sig = zones[z].chain;
        std::vector<int> rev(sig.rbegin(), sig.rend());
        if (rev < sig) sig = std::move(rev);
        bySignature[sig].push_back(z);
    }

    for (auto &kv : bySignature) {
        const std::vector<int> &members = kv.second;
        for (size_t i = 0; i < members.size(); i += 2) {
            CoarseEdge edge;
            edge.chain = kv.first;
            const auto it = coarseBranch_.find(pairKey(edge.chain.front(), edge.chain.back()));
            if (it != coarseBranch_.end()) edge.branch = it->second;
            edge.zones[0] = members[i];
            zones[members[i]].coarseEdge = static_cast<int>(edges.size());
            if (i + 1 < members.size()) {
                edge.zones[1] = members[i + 1];
                zones[members[i + 1]].coarseEdge = static_cast<int>(edges.size());
            } else {
                ++stats_.unpairedEdges;
            }
            edges.push_back(std::move(edge));
        }
    }
}

void MedialAxisTMesh::classify() {
    // Caps first: a blue cap needs its corner p, which is also what makes a
    // merged endpoint blue rather than purple. A triangle's section is two
    // boundary edges with the corner between them; a merged cell's section
    // is blue only if it holds a locked convex corner ([R]), the sharpest of
    // them serving as p.
    std::unordered_map<int, int> capOf;
    for (int z = 0; z < static_cast<int>(zones.size()); ++z) {
        MedialZone &zone = zones[z];
        if (!zone.cap) continue;
        capOf[zone.capVertex] = z;
        const int cell = axis->medialVertices[zone.capVertex].cell;
        const bool triangle = cell >= 0 && axis->cells[cell].dualIsTriangle();
        int corner = -1;
        double sharpest = 0.0;
        const int np = static_cast<int>(zone.boundarySide.size());
        for (int i = 1; i + 1 < np; ++i) {
            const int v = zone.boundaryVertexAt[i];
            if (v < 0) continue;
            if (triangle) { corner = i; break; }
            const double turn = originalTurn(v);
            if (axis->locked[v] && turn > sharpest) {
                sharpest = turn;
                corner = i;
            }
        }
        zone.capCorner = corner;
        zone.color = corner >= 0 ? MedialColor::Blue : MedialColor::Purple;
        if (corner >= 0 && !triangle) ++stats_.blueByCorner;
    }

    for (CoarseEdge &edge : edges) {
        const int m1 = edge.chain.front();
        const int m2 = edge.chain.back();
        int endpoint = -1;
        if (nodeDegree(m1) == 1) endpoint = m1;
        else if (nodeDegree(m2) == 1) endpoint = m2;

        if (endpoint >= 0) {
            // An edge incident to an endpoint takes the endpoint's class.
            const auto it = capOf.find(endpoint);
            edge.color = (it != capOf.end()) ? zones[it->second].color : MedialColor::Purple;
        } else if (edge.zones[1] >= 0) {
            // theta_{m,E} at each end: the angle between the radii of the
            // two contacts beside E, which are the spokes of E's two zones
            // there. A zone's boundary side runs with its chain, so a zone
            // walking m1 -> m2 has its m1 spoke at the front.
            for (int end = 0; end < 2; ++end) {
                const int m = end == 0 ? m1 : m2;
                const Point c = nodeCoord(m);
                Point feet[2];
                for (int k = 0; k < 2; ++k) {
                    const MedialZone &z = zones[edge.zones[k]];
                    feet[k] = (z.chain.front() == m) ? z.boundarySide.front()
                                                     : z.boundarySide.back();
                }
                edge.endAngle[end] = angleBetween(feet[0] - c, feet[1] - c);
            }
            edge.meanAngle = 0.5 * (edge.endAngle[0] + edge.endAngle[1]);
            edge.color = edge.meanAngle >= MEDIAL_COLOR_ANGLE_THRESHOLD ? MedialColor::Green
                                                                        : MedialColor::Red;
        }
        ++stats_.colorCounts[static_cast<int>(edge.color)];

        for (int z : edge.zones) {
            if (z >= 0) zones[z].color = edge.color;
        }
    }
}

// ─── Step 8: templates (Fig. 17) ────────────────────────────────────────────
//
// Each class of coarse edge remeshes its subdomain with one quadrilateral
// template, written exactly as docs/medial_axis.md Step 8 lists them: E =
// (m1, m2) with m1 the end the template treats specially, L and R the sides
// on the left and right walking m1 -> m2, b1/b2 the preimages at m1/m2, and
// s_sigma the midpoint of the radius [m2, b2_sigma]. Every quad comes out
// CCW, so a cell with non-positive area is a folded one (Step 8.6).

namespace {

// Assembles one block from its sides. `corner` adds a single vertex; `run`
// appends a polyline and marks both its ends as corners. Points repeated
// across consecutive calls are welded, so sharing an endpoint between two
// sides costs nothing at the call site.
class BlockBuilder {
public:
    BlockBuilder(std::vector<TMeshBlock> &out, MedialColor color, int edge, double weld)
        : out_(out), color_(color), edge_(edge), weld_(weld) {}

    BlockBuilder &corner(const Point &p) {
        push(p, true);
        return *this;
    }

    BlockBuilder &run(const std::vector<Point> &pts) {
        for (size_t i = 0; i < pts.size(); ++i) {
            push(pts[i], i == 0 || i + 1 == pts.size());
        }
        return *this;
    }

    // Index of the emitted block, or -1 if it had too few points to be one.
    int emit() {
        // The outline is a closed loop, so a last point on top of the first
        // is a duplicate rather than a side.
        if (block_.outline.size() > 1 &&
            normP(block_.outline.front() - block_.outline.back()) < weld_) {
            block_.outline.pop_back();
        }
        if (block_.corners.size() > 1 &&
            normP(block_.corners.front() - block_.corners.back()) < weld_) {
            block_.corners.pop_back();
        }
        int index = -1;
        if (block_.outline.size() >= 3) {
            block_.color = color_;
            block_.edge = edge_;
            index = static_cast<int>(out_.size());
            out_.push_back(std::move(block_));
        }
        block_ = TMeshBlock{};
        return index;
    }

private:
    void push(const Point &p, bool isCorner) {
        if (block_.outline.empty() || normP(p - block_.outline.back()) >= weld_) {
            block_.outline.push_back(p);
        }
        if (isCorner && (block_.corners.empty() ||
                         normP(p - block_.corners.back()) >= weld_)) {
            block_.corners.push_back(p);
        }
    }

    std::vector<TMeshBlock> &out_;
    MedialColor color_;
    int edge_;
    double weld_;
    TMeshBlock block_;
};

// A zone kept whole, for what no template claimed, CCW like the rest.
TMeshBlock zoneAsBlock(const MedialAxisTMesh &tm, const MedialZone &zone, double weld) {
    TMeshBlock b;
    b.outline = zone.boundarySide;
    std::vector<Point> chainPts;
    for (int m : zone.chain) chainPts.push_back(tm.nodeCoord(m));
    // A zone that leaves a vertex and comes back to it -- around a tip whose
    // branch the inside filter dropped, so nothing pairs with it, or beside
    // the axis of a disc, whose two ends are one point -- is shaped like a
    // cap, and is closed like one.
    if (zone.cap || normP(chainPts.front() - chainPts.back()) <= weld) {
        // Open caps close through their vertex; a cap spanning a whole closed
        // loop is already a cycle and needs no spoke at all.
        if (normP(zone.boundarySide.front() - zone.boundarySide.back()) > weld) {
            b.outline.push_back(chainPts.front());
        }
        b.corners = {zone.boundarySide.front(), zone.boundarySide.back(), chainPts.front()};
    } else {
        appendPts(b.outline, {chainPts.back(), chainPts.front()});
        b.corners = {zone.boundarySide.front(), zone.boundarySide.back(),
                     chainPts.back(), chainPts.front()};
    }
    if (signedArea(b.outline) < 0.0) {
        std::reverse(b.outline.begin(), b.outline.end());
        std::reverse(b.corners.begin(), b.corners.end());
    }
    b.color = zone.color;
    return b;
}

} // namespace

void MedialAxisTMesh::buildBlocks() {
    blocks.clear();

    std::unordered_map<int, int> capOf;
    for (int z = 0; z < static_cast<int>(zones.size()); ++z) {
        if (zones[z].cap) capOf[zones[z].capVertex] = z;
    }
    std::vector<char> consumed(zones.size(), 0);
    const double weld = 1e-12 * axis->boundingBoxDiagonal();

    // Orientation of every template: m1 is the end it treats specially --
    // the smaller medial radius for red ([R]: along red edges r is strictly
    // monotone, App. A, so the inner end is well defined), the endpoint for
    // blue and purple, either end for green. R walks m1 -> m2 (a zone lies
    // right of its own chain) and L walks m2 -> m1; anything else means the
    // pairing went wrong and the zones are kept as they are.
    struct Plan {
        int m1 = -1, m2 = -1;
        const MedialZone *R = nullptr;
        const MedialZone *L = nullptr;
    };
    std::vector<Plan> plans(edges.size());
    for (int e = 0; e < static_cast<int>(edges.size()); ++e) {
        const CoarseEdge &edge = edges[e];
        if (edge.zones[1] < 0) continue;
        int m1 = edge.chain.front();
        int m2 = edge.chain.back();
        if (edge.color == MedialColor::Red) {
            if (nodeRadius(m2) < nodeRadius(m1)) std::swap(m1, m2);
        } else if (edge.color == MedialColor::Blue || edge.color == MedialColor::Purple) {
            const bool endFront = nodeDegree(m1) == 1 && capOf.count(m1);
            const bool endBack = nodeDegree(m2) == 1 && capOf.count(m2);
            if (!endFront && endBack) std::swap(m1, m2);
        }
        const MedialZone &A = zones[edge.zones[0]];
        const MedialZone &B = zones[edge.zones[1]];
        Plan &plan = plans[e];
        plan.m1 = m1;
        plan.m2 = m2;
        if (A.chain.front() == m1 && A.chain.back() == m2 &&
            B.chain.front() == m2 && B.chain.back() == m1) {
            plan.R = &A;
            plan.L = &B;
        } else if (B.chain.front() == m1 && B.chain.back() == m2 &&
                   A.chain.front() == m2 && A.chain.back() == m1) {
            plan.R = &B;
            plan.L = &A;
        }
    }

    // The points s_sigma on the radii at m2. Assembly item 1 creates one per
    // radius, shared by every template that asks for it; where two ask --
    // the two sides of a radius at a junction both outer, or the valence-2
    // case below -- the nearer to m2 serves both.
    //
    // [R] The guide puts s_sigma at the radius midpoint, read off its
    // canonical wedge where r(m1) = r(m2) / 2. There the midpoint is also the
    // foot of the perpendicular from m1 onto the radius, and for other ratios
    // the midpoint can fold the diamond D = (m1, sR, m2, sL) at m1: a point
    // at distance t from m2 along a radius that leans back towards m1 by
    // theta/2 keeps m1 outside the triangle (m2, sR, sL) only while
    // t cos(theta/2) < |m1 - m2|, i.e. t * foot < |m1 - m2|^2 with foot the
    // projection of m1 - m2 on the radius -- which every red edge more than
    // three samples from its apex breaks. So the midpoint stays while it is
    // within half that bound, and the foot takes its place beyond: with the
    // feet D has right angles at sR and sL, is inscribed in the circle on
    // m1-m2 and so convex, and a wedge comes out as a grid of squares. A
    // radius at right angles to the edge -- a purple cap on a strip of
    // constant width -- has no foot to speak of and puts no bound on t, so
    // its midpoint always stays.
    std::vector<double> spokeT(spokes_.size(), std::numeric_limits<double>::infinity());
    auto request = [&](int spoke, int m1) {
        if (spoke < 0) return;
        const Spoke &sp = spokes_[spoke];
        const Point c = nodeCoord(sp.node);
        const double len = normP(sp.foot - c);
        if (!(len > 0.0)) return;
        const Point u = (sp.foot - c) / len;
        const Point w = nodeCoord(m1) - c;
        const double foot = dotP(w, u);
        double t = 0.5 * len;
        if (foot > 0.0 && t * foot > 0.5 * dotP(w, w)) t = std::min(t, foot);
        spokeT[spoke] = std::min(spokeT[spoke], t);
    };
    auto spokePoint = [&](int spoke) {
        const Spoke &sp = spokes_[spoke];
        const Point c = nodeCoord(sp.node);
        const double len = normP(sp.foot - c);
        const double t = std::isfinite(spokeT[spoke]) ? spokeT[spoke] : 0.5 * len;
        return len > 0.0 ? c + (sp.foot - c) * (t / len) : c;
    };
    for (int e = 0; e < static_cast<int>(edges.size()); ++e) {
        const Plan &plan = plans[e];
        if (!plan.R || edges[e].color == MedialColor::Green) continue;
        request(plan.R->spokes[1], plan.m1);
        request(plan.L->spokes[0], plan.m1);
    }

    // The diamond D of every non-green template, for the valence-2 pass.
    std::vector<int> diamondOf(edges.size(), -1);
    for (int e = 0; e < static_cast<int>(edges.size()); ++e) {
        const CoarseEdge &edge = edges[e];
        const Plan &plan = plans[e];
        if (!plan.R) continue;

        const std::vector<Point> &arcR = plan.R->boundarySide;   // b1R -> b2R
        const std::vector<Point> &arcL = plan.L->boundarySide;   // b2L -> b1L
        const Point P1 = nodeCoord(plan.m1);
        const Point P2 = nodeCoord(plan.m2);
        // An edge of no length has no subdomain to remesh: the axis of a
        // disc, both its endpoints on the centre. Its caps stay as they are.
        if (normP(P2 - P1) <= weld) {
            ++stats_.pointEdges;
            continue;
        }
        consumed[edge.zones[0]] = 1;
        consumed[edge.zones[1]] = 1;

        BlockBuilder b(blocks, edge.color, e, weld);
        if (edge.color == MedialColor::Green) {
            // Keep the initial decomposition: QL, QR, the medial edge a shared
            // mesh edge, straight between its coarse ends (Step 6.4, [R]).
            b.corner(P1).corner(P2).run(arcL).emit();                // (m1, m2, b2L, b1L)
            b.corner(P1).run(arcR).corner(P2).emit();                // (m1, b1R, b2R, m2)
            continue;
        }
        const Point sR = spokePoint(plan.R->spokes[1]);
        const Point sL = spokePoint(plan.L->spokes[0]);

        // The polar region at an endpoint m1, running aL -> aR with the walk:
        // it starts where L's boundary ends and ends where R's begins.
        std::vector<Point> cap;
        int capCorner = -1;
        if (edge.color != MedialColor::Red) {
            const auto it = capOf.find(plan.m1);
            if (it != capOf.end()) {
                const MedialZone &C = zones[it->second];
                if (normP(C.boundarySide.front() - arcL.back()) <= weld &&
                    normP(C.boundarySide.back() - arcR.front()) <= weld) {
                    cap = C.boundarySide;
                    capCorner = C.capCorner;
                    consumed[it->second] = 1;
                }
            }
        }

        if (edge.color == MedialColor::Blue && !cap.empty() && capCorner > 0) {
            // C = (m_e, aL, p, aR): the endpoint stays regular (valence 4)
            // and the corner p is a boundary node of valence 2.
            const std::vector<Point> toCorner(cap.begin(), cap.begin() + capCorner + 1);
            const std::vector<Point> fromCorner(cap.begin() + capCorner, cap.end());
            b.corner(P1).run(toCorner).run(fromCorner).emit();
            b.corner(P1).corner(sL).run(arcL).emit();               // (m_e, sL, b2L, aL)
            b.corner(P1).run(arcR).corner(sR).emit();               // (m_e, aR, b2R, sR)
        } else if (edge.color == MedialColor::Purple && !cap.empty()) {
            // Blue with aL and aR merged into one apex q, the arc-length
            // midpoint of the section: C disappears, m_e has valence 3 and
            // q boundary valence 3, so the cap adds no boundary singularity.
            std::vector<Point> toQ, fromQ;
            splitAtFraction(cap, 0.5, toQ, fromQ);
            std::vector<Point> left = arcL;
            appendPts(left, toQ);
            std::vector<Point> right = fromQ;
            appendPts(right, arcR);
            b.corner(P1).corner(sL).run(left).emit();               // (m_e, sL, b2L, q)
            b.corner(P1).run(right).corner(sR).emit();              // (m_e, q, b2R, sR)
        } else {
            // Red, m1 inner and m2 outer -- or an endpoint whose cap could
            // not be matched, meshed the same way around whatever contacts
            // it has.
            b.corner(P1).corner(sL).run(arcL).emit();               // (m1, sL, b2L, b1L)
            b.corner(P1).run(arcR).corner(sR).emit();               // (m1, b1R, b2R, sR)
        }
        // The diamond: the medial edge m1-m2 is its diagonal, not a mesh edge.
        diamondOf[e] = b.corner(P1).corner(sR).corner(P2).corner(sL).emit();
    }

    // Whatever the templates did not absorb -- zones of unpaired edges, caps
    // no edge matched -- still tiles part of the domain, so it stays a block,
    // unless it encloses nothing (the two sides of an edge of no length).
    const double minArea = 1e-12 * targetSize_ * targetSize_;
    for (int z = 0; z < static_cast<int>(zones.size()); ++z) {
        if (consumed[z]) continue;
        TMeshBlock b = zoneAsBlock(*this, zones[z], weld);
        if (signedArea(b.outline) > minArea) blocks.push_back(std::move(b));
    }

    dissolveValenceTwo(diamondOf);
}

// Step 8.5: an interior vertex where two templates both put their outer end
// -- two non-green edges meeting at a regular vertex, each with it as m2 --
// is touched only by the two diamonds, along the two radii they share: a
// vertex of valence 2. Merging the two diamonds removes it ([R]); the
// merged cell is the quad on the two diamonds' far corners.
void MedialAxisTMesh::dissolveValenceTwo(const std::vector<int> &diamondOf) {
    std::unordered_map<int, std::vector<int>> outerAt;
    for (int e = 0; e < static_cast<int>(edges.size()); ++e) {
        const int d = diamondOf[e];
        if (d < 0 || blocks[d].corners.size() != 4) continue;
        // D = (m1, sR, m2, sL): the outer end is corner 2.
        const Point outer = blocks[d].corners[2];
        const int m = normP(nodeCoord(edges[e].chain.front()) - outer) <
                              normP(nodeCoord(edges[e].chain.back()) - outer)
                          ? edges[e].chain.front()
                          : edges[e].chain.back();
        outerAt[m].push_back(d);
    }

    const double weld = 1e-12 * axis->boundingBoxDiagonal();
    std::vector<char> drop(blocks.size(), 0);
    for (const auto &kv : outerAt) {
        if (kv.second.size() != 2 || nodeDegree(kv.first) != 2) continue;
        // Rotated so the shared vertex comes last: D1 = (a1, b1, c1, m) and
        // D2 = (a2, b2, c2, m) share m-a1 = m-c2 and m-c1 = m-a2.
        std::array<std::array<Point, 4>, 2> q;
        for (int k = 0; k < 2; ++k) {
            const std::vector<Point> &c = blocks[kv.second[k]].corners;
            for (int i = 0; i < 4; ++i) q[k][i] = c[(i + 3) % 4];
        }
        if (normP(q[0][0] - q[1][2]) > weld || normP(q[0][2] - q[1][0]) > weld) continue;
        TMeshBlock &merged = blocks[kv.second[0]];
        merged.corners = {q[0][0], q[0][1], q[0][2], q[1][1]};
        merged.outline = merged.corners;
        drop[kv.second[1]] = 1;
        ++stats_.dissolvedVertices;
    }
    std::vector<TMeshBlock> kept;
    kept.reserve(blocks.size());
    for (size_t i = 0; i < blocks.size(); ++i) {
        if (!drop[i]) kept.push_back(std::move(blocks[i]));
    }
    blocks = std::move(kept);
}

// Step 8.6: every face has positive area and is not severely non-convex. Two
// corners closer than the T-mesh welder's tolerance (1e-4 h, see
// QuantTMeshConvert.hxx) count as neither: they become one node, that side
// has no length, and the cell is the collocated fan Step 11 exists to repair
// at a concave corner.
std::vector<int> MedialAxisTMesh::checkGeometry() {
    std::vector<int> branches;
    stats_.badFaces = 0;
    const double minArea = 1e-10 * targetSize_ * targetSize_;
    const double minSide = 1e-4 * targetSize_;
    for (const TMeshBlock &b : blocks) {
        if (b.edge < 0 || b.corners.size() != 4) continue;
        bool bad = signedArea(b.outline) <= minArea;
        for (int i = 0; i < 4 && !bad; ++i) {
            const Point toNext = b.corners[(i + 1) % 4] - b.corners[i];
            const Point toPrev = b.corners[(i + 3) % 4] - b.corners[i];
            if (normP(toNext) <= minSide) {
                bad = true;
                break;
            }
            double interior = std::atan2(cross2(toNext, toPrev), dotP(toNext, toPrev));
            if (interior < 0.0) interior += 2.0 * M_PI;
            if (interior > opts_.maxCornerAngle) bad = true;
        }
        if (!bad) continue;
        ++stats_.badFaces;
        const int branch = edges[b.edge].branch;
        if (branch >= 0 && std::find(branches.begin(), branches.end(), branch) == branches.end()) {
            branches.push_back(branch);
        }
    }
    return branches;
}
