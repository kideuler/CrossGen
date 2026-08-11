#include "medialaxis/MedialAxisTMesh.hxx"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <unordered_map>

namespace {

// Polyline helpers for the template construction.

std::vector<Point> reversed(std::vector<Point> pts) {
    std::reverse(pts.begin(), pts.end());
    return pts;
}

// Append, dropping a duplicated joint point.
void appendPts(std::vector<Point> &dst, const std::vector<Point> &src) {
    for (const Point &p : src) {
        if (!dst.empty() && normP(p - dst.back()) < 1e-14) continue;
        dst.push_back(p);
    }
}

// Split a polyline at the given increasing arc-length fractions. Fills
// `pieces` (fractions.size() + 1 of them, each sharing its cut points with
// its neighbours) and `cuts` (the split points themselves).
void splitPolyline(const std::vector<Point> &pts,
                   const std::vector<double> &fractions,
                   std::vector<std::vector<Point>> &pieces,
                   std::vector<Point> &cuts) {
    pieces.clear();
    cuts.clear();
    std::vector<double> arc(pts.size(), 0.0);
    for (size_t i = 1; i < pts.size(); ++i) {
        arc[i] = arc[i - 1] + normP(pts[i] - pts[i - 1]);
    }
    const double total = arc.back();

    std::vector<Point> cur{pts.front()};
    size_t i = 1;
    for (double f : fractions) {
        const double target = total * f;
        while (i < pts.size() && arc[i] < target) {
            cur.push_back(pts[i]);
            ++i;
        }
        Point cut = pts.back();
        if (i < pts.size()) {
            const double seg = arc[i] - arc[i - 1];
            const double t = seg > 0.0 ? (target - arc[i - 1]) / seg : 0.0;
            cut = pts[i - 1] + (pts[i] - pts[i - 1]) * t;
        }
        cur.push_back(cut);
        cuts.push_back(cut);
        pieces.push_back(std::move(cur));
        cur = {cut};
    }
    for (; i < pts.size(); ++i) cur.push_back(pts[i]);
    pieces.push_back(std::move(cur));
}

// One step of the boundary walk: either a plain boundary vertex (medial < 0)
// or the footpoint of a medial vertex passage.
struct WalkEvent {
    Point pos{0.0, 0.0};
    int medial = -1;
};

// A maximal group of consecutive medial events sharing one medial vertex,
// plain boundary events in between included. A run longer than one event is
// a stretch of boundary collapsing onto that vertex -- a polar section.
struct MedialRun {
    int m = -1;
    int first = -1;   // event index of the first passage footpoint
    int last = -1;    // event index of the last one (== first for a point)
};

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

// The fan's contribution to the walk, in boundary order. Each fan emits its
// chain vertices 0..k-2 (the last belongs to the next fan, which starts with
// it), plus the boundary vertex itself slotted in by arc length. A collapsed
// fan emits its single vertex at the interval entry, so consecutive collapsed
// fans build up the run that becomes a cap zone.
void appendFanEvents(const MedialFan &fan, const Point &vpos,
                     std::vector<WalkEvent> &events) {
    const int k = static_cast<int>(fan.chain.size());
    if (k == 1) {
        events.push_back({fan.midPrev, fan.chain[0]});
        events.push_back({vpos, -1});
        return;
    }

    const double vparam = normP(vpos - fan.midPrev);
    bool vEmitted = false;
    for (int i = 0; i + 1 < k; ++i) {
        const double p = subPolylineParam(fan.footpoints[i], fan.midPrev, vpos,
                                          fan.midNext);
        if (!vEmitted && p > vparam) {
            events.push_back({vpos, -1});
            vEmitted = true;
        }
        events.push_back({fan.footpoints[i], fan.chain[i]});
    }
    if (!vEmitted) events.push_back({vpos, -1});
}

} // namespace

// ─── Construction ───────────────────────────────────────────────────────────

MedialAxisTMesh::MedialAxisTMesh(const MedialAxisMap &map, double targetSize)
    : axis(map.axis) {
    if (!axis) {
        throw std::runtime_error("MedialAxisMap carries no axis.");
    }

    targetSize_ = targetSize > 0.0 ? targetSize : 0.1 * axis->boundingBoxDiagonal();
    passages_.assign(axis->medialVertices.size(), {});

    selectKeptVertices();

    std::vector<char> visited(axis->mesh->vertices.size(), 0);
    for (int v = 0; v < static_cast<int>(axis->mesh->vertices.size()); ++v) {
        if (visited[v] || axis->boundaryLinks[v].nextBoundaryVertex < 0) continue;
        walkLoop(map, v, visited);
    }

    buildEdges();
    classify();
    buildBlocks();

    stats_.zones = static_cast<int>(zones.size());
    stats_.coarseEdges = static_cast<int>(edges.size());
    stats_.blocks = static_cast<int>(blocks.size());
}

// ─── Downsampling ───────────────────────────────────────────────────────────

void MedialAxisTMesh::selectKeptVertices() {
    kept_.assign(axis->medialVertices.size(), 0);

    // Branches are the unit of downsampling, so make sure they exist.
    if (axis->polyLines.empty()) axis->createPolylines();

    for (size_t b = 0; b < axis->polyLines.size(); ++b) {
        const std::vector<int> &pl = axis->polyLines[b];
        kept_[pl.front()] = 1;
        kept_[pl.back()] = 1;
        if (pl.size() < 2) continue;

        std::vector<double> arc(pl.size(), 0.0);
        for (size_t i = 1; i < pl.size(); ++i) {
            arc[i] = arc[i - 1] + normP(axis->medialVertices[pl[i]].coord -
                                        axis->medialVertices[pl[i - 1]].coord);
        }
        const double length = arc.back();

        int segments = std::max(1, static_cast<int>(std::lround(length / targetSize_)));
        // A cycle's ends coincide, so fewer than three samples would leave a
        // hole bounded by fewer than two zones per side.
        if (axis->polyLineIsCycle[b]) segments = std::max(3, segments);

        size_t cursor = 1;
        for (int s = 1; s < segments; ++s) {
            const double target = length * s / segments;
            while (cursor + 1 < pl.size() - 1 && arc[cursor] < target) ++cursor;
            // Nearest of the two candidates around the target, clamped to the
            // branch interior (the ends are already kept).
            size_t pick = cursor;
            if (cursor > 1 &&
                target - arc[cursor - 1] < arc[cursor] - target) pick = cursor - 1;
            if (pick >= 1 && pick + 1 <= pl.size() - 1) kept_[pl[pick]] = 1;
        }
    }

    for (char k : kept_) stats_.keptVertices += k;
}

// ─── The walk ───────────────────────────────────────────────────────────────

void MedialAxisTMesh::walkLoop(const MedialAxisMap &map, int startVertex,
                               std::vector<char> &visited) {
    ++stats_.loops;

    // Gather the loop's event stream.
    std::vector<WalkEvent> events;
    int v = startVertex;
    int guard = static_cast<int>(axis->mesh->vertices.size()) + 1;
    while (guard-- > 0) {
        visited[v] = 1;
        const MedialFan *fan = map.fanOf(v);
        if (fan) {
            appendFanEvents(*fan, axis->mesh->vertices[v], events);
        } else {
            events.push_back({axis->mesh->vertices[v], -1});
        }
        v = axis->boundaryLinks[v].nextBoundaryVertex;
        if (v == startVertex) break;
        if (v < 0 || visited[v]) {
            // Broken or pinched links: the loop never closes, and zones cut
            // from it would not either.
            ++stats_.skippedLoops;
            return;
        }
    }

    // Rotate the stream so it starts at a run boundary -- the first medial
    // event whose vertex differs from the previous medial event's -- so runs
    // never straddle the wrap-around.
    const int n = static_cast<int>(events.size());
    int lastMedial = -1;
    for (int i = n - 1; i >= 0; --i) {
        if (events[i].medial >= 0) { lastMedial = events[i].medial; break; }
    }
    if (lastMedial < 0) {
        // A loop no fan reached at all.
        ++stats_.skippedLoops;
        return;
    }
    int start = -1;
    int prevMedial = lastMedial;
    for (int i = 0; i < n; ++i) {
        if (events[i].medial < 0) continue;
        if (events[i].medial != prevMedial) { start = i; break; }
        prevMedial = events[i].medial;
    }

    std::vector<MedialRun> runs;
    if (start < 0) {
        // Every passage hits the same medial vertex: the whole loop is one
        // polar section (a domain reduced to a single inscribed circle).
        runs.push_back({lastMedial, 0, n - 1});
    } else {
        std::rotate(events.begin(), events.begin() + start, events.end());
        for (int i = 0; i < n; ++i) {
            if (events[i].medial < 0) continue;
            if (!runs.empty() && runs.back().m == events[i].medial) {
                runs.back().last = i;
            } else {
                runs.push_back({events[i].medial, i, i});
            }
        }
    }

    // Record the passages: one footpoint per run, at the middle of the run's
    // interval. These are the discrete medial radii the classification needs.
    for (const MedialRun &r : runs) {
        passages_[r.m].push_back((events[r.first].pos + events[r.last].pos) * 0.5);
    }

    std::vector<int> cuts;   // indices into runs
    for (int i = 0; i < static_cast<int>(runs.size()); ++i) {
        if (kept_[runs[i].m]) cuts.push_back(i);
    }
    if (cuts.empty()) {
        // No kept vertex on this loop -- degenerate sampling. Nothing sound
        // to cut, so report rather than fabricate a zone.
        ++stats_.skippedLoops;
        return;
    }

    // Inclusive cyclic slice of event positions.
    auto slice = [&](int from, int to) {
        std::vector<Point> pts;
        for (int i = from;; i = (i + 1) % n) {
            pts.push_back(events[i].pos);
            if (i == to) break;
        }
        return pts;
    };

    // Cap zones: every kept run that covers more than a point.
    for (int c : cuts) {
        const MedialRun &r = runs[c];
        if (r.first == r.last) continue;
        MedialZone zone;
        zone.boundarySide = slice(r.first, r.last);
        // A run spanning the whole loop is a domain reduced to one inscribed
        // circle; close the polyline so the cap covers the full perimeter.
        if (runs.size() == 1) zone.boundarySide.push_back(events[r.first].pos);
        zone.chain = {r.m};
        zone.cap = true;
        zone.capVertex = r.m;
        zones.push_back(std::move(zone));
        ++stats_.capZones;
    }

    // Quad zones: the interval between each pair of consecutive cuts, with
    // the medial vertices passed in between as the chain. A single cut on the
    // loop (a domain whose axis is one kept vertex) leaves no interval.
    if (cuts.size() < 2 && runs.size() == 1) return;
    const int nc = static_cast<int>(cuts.size());
    for (int j = 0; j < nc; ++j) {
        const MedialRun &a = runs[cuts[j]];
        const MedialRun &b = runs[cuts[(j + 1) % nc]];
        MedialZone zone;
        zone.boundarySide = slice(a.last, b.first);
        zone.chain.push_back(a.m);
        for (int r = cuts[j] + 1;; ++r) {
            const int ri = r % static_cast<int>(runs.size());
            if (ri == cuts[(j + 1) % nc]) break;
            zone.chain.push_back(runs[ri].m);
        }
        zone.chain.push_back(b.m);
        zones.push_back(std::move(zone));
    }
}

// ─── Edges and classification ───────────────────────────────────────────────

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
    // Blue/purple caps: the polar region around an endpoint, split by what
    // its dual cell is -- a triangle marks a tip that stayed sharp, anything
    // larger a round cap the simplification merged together.
    for (MedialZone &zone : zones) {
        if (!zone.cap) continue;
        zone.color = axis->medialVertices[zone.capVertex].dualIsTriangle
                         ? MedialColor::Blue
                         : MedialColor::Purple;
        ++stats_.colorCounts[static_cast<int>(zone.color)];
    }

    for (CoarseEdge &edge : edges) {
        // Sec. 4.2: an edge incident to an endpoint takes the endpoint's
        // class; every other edge splits green/red on the medial angle.
        int endpoint = -1;
        if (axis->medialVertices[edge.chain.front()].degree <= 1)
            endpoint = edge.chain.front();
        else if (axis->medialVertices[edge.chain.back()].degree <= 1)
            endpoint = edge.chain.back();

        if (endpoint >= 0) {
            edge.color = axis->medialVertices[endpoint].dualIsTriangle
                             ? MedialColor::Blue
                             : MedialColor::Purple;
        } else {
            // The medial angle at each chain vertex passed exactly twice: the
            // angle its two radii (per-side passage footpoints) span.
            edge.minAngle = M_PI;
            double sum = 0.0;
            int measured = 0;
            for (int m : edge.chain) {
                if (passages_[m].size() != 2) continue;
                const Point &c = axis->medialVertices[m].coord;
                const Point a = passages_[m][0] - c;
                const Point b = passages_[m][1] - c;
                const double theta = std::atan2(std::abs(cross2(a, b)), dotP(a, b));
                edge.minAngle = std::min(edge.minAngle, theta);
                sum += theta;
                ++measured;
            }
            edge.meanAngle = measured > 0 ? sum / measured : M_PI;
            edge.color = edge.meanAngle >= MEDIAL_COLOR_ANGLE_THRESHOLD
                             ? MedialColor::Green
                             : MedialColor::Red;
        }
        ++stats_.colorCounts[static_cast<int>(edge.color)];

        for (int z : edge.zones) {
            if (z >= 0) zones[z].color = edge.color;
        }
    }
}

// ─── Templates (Fig. 17) ────────────────────────────────────────────────────
//
// Each class of coarse edge remeshes its subdomain with one quadrilateral
// template. The machinery below keeps that table-driven: a template is a
// function over a TemplateSubdomain -- the subdomain's geometry, oriented so
// that the end the template treats specially comes first -- which emits its
// blocks through a BlockBuilder.
//
// To add a template: write one such function and add a row to kTemplates
// naming the class it serves, which end it wants first, and whether it
// swallows the polar cap at that end. Nothing else here needs to change.

namespace {

// Which end of the coarse edge a template wants as `inner`.
enum class Orient {
    Any,       // symmetric template; take the edge as it comes
    Narrow,    // the smaller medial radius, i.e. the subdomain's concave corner
    Endpoint,  // the end that terminates the axis
};

// The subdomain of one coarse edge, oriented for its template.
struct TemplateSubdomain {
    Point inner{0.0, 0.0};   // the end the template treats specially
    Point outer{0.0, 0.0};   // the other one

    std::vector<Point> chain;         // every medial vertex, inner -> outer
    std::vector<Point> sideA, sideB;  // the two boundary runs, inner -> outer

    // The polar cap at `inner`, split at its midpoint W and joined to the
    // flanks: arcUp runs outer's contact on side A to W, and arcLo runs W to
    // outer's contact on side B, so together they are the entire boundary of
    // a capped subdomain. Only filled for templates that asked for the cap.
    std::vector<Point> arcUp, arcLo;

    // The cap zone folded in above, for the caller to mark consumed.
    int capZone = -1;

    // Midpoints of the two radii at `outer` -- where a template that leaves a
    // T-junction on those radii puts it.
    Point midRadiusA() const { return (sideA.back() + outer) * 0.5; }
    Point midRadiusB() const { return (sideB.back() + outer) * 0.5; }
};

// Assembles one block from its sides. `corner` adds a single vertex; `run`
// appends a polyline and marks both its ends as corners. Points repeated
// across consecutive calls are welded, so sharing an endpoint between two
// sides costs nothing at the call site.
class BlockBuilder {
public:
    BlockBuilder(std::vector<TMeshBlock> &out, MedialColor color, double weld)
        : out_(out), color_(color), weld_(weld) {}

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

    void emit() {
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
        if (block_.outline.size() >= 3) {
            block_.color = color_;
            out_.push_back(std::move(block_));
        }
        block_ = TMeshBlock{};
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
    double weld_;
    TMeshBlock block_;
};

// ── The four templates ──

// Green: the axis runs along quad edges, so the two zones already are the
// template -- one quad per side, sharing the medial chain.
void templateGreen(const TemplateSubdomain &s, BlockBuilder &b) {
    b.run(s.sideA).run(reversed(s.chain)).emit();
    b.run(s.sideB).run(reversed(s.chain)).emit();
}

// Red: three quads around the concave corner.
//
// At the narrow medial vertex both radii point back toward the apex, so the
// subdomain is reflex there -- that is its one concave corner. Splitting it
// means running an edge from that corner to each of the two opposite sides,
// which here are the radii at `outer`. Each lands mid-radius and so leaves a
// T-junction, the subdomain across that radius having no vertex there. The
// medial edge ends up as the diagonal of the middle quad, which is what red
// asks for: the axis follows quad diagonals.
void templateRed(const TemplateSubdomain &s, BlockBuilder &b) {
    const Point tA = s.midRadiusA();
    const Point tB = s.midRadiusB();
    b.corner(s.inner).run(s.sideA).corner(tA).emit();
    b.corner(s.inner).corner(tA).corner(s.outer).corner(tB).emit();
    b.corner(s.inner).corner(tB).run(reversed(s.sideB)).emit();
}

// Blue: one quad spanning the wedge at a sharp tip, the axis on its diagonal.
void templateBlue(const TemplateSubdomain &s, BlockBuilder &b) {
    b.corner(s.outer).run(s.arcUp).run(s.arcLo).emit();
}

// Purple: red's construction at a round cap. The T-junctions and the
// half-radius running from each one out to the boundary are the same; what
// differs is the far end, where an endpoint has no pair of inner radii to
// serve as a concave corner, so its polar cap stands in for them -- split at
// W, the single boundary point the endpoint maps to.
void templatePurple(const TemplateSubdomain &s, BlockBuilder &b) {
    const Point tA = s.midRadiusA();
    const Point tB = s.midRadiusB();
    b.run(s.arcUp).corner(s.inner).corner(tA).emit();
    b.run(s.arcLo).corner(tB).corner(s.inner).emit();
    b.corner(s.inner).corner(tA).corner(s.outer).corner(tB).emit();
}

struct TemplateEntry {
    MedialColor color;
    Orient orient;
    bool wantsCap;
    void (*apply)(const TemplateSubdomain &, BlockBuilder &);
};

const TemplateEntry kTemplates[] = {
    {MedialColor::Green,  Orient::Any,      false, templateGreen},
    {MedialColor::Red,    Orient::Narrow,   false, templateRed},
    {MedialColor::Blue,   Orient::Endpoint, true,  templateBlue},
    {MedialColor::Purple, Orient::Endpoint, true,  templatePurple},
};

const TemplateEntry &templateFor(MedialColor c) {
    for (const TemplateEntry &t : kTemplates) {
        if (t.color == c) return t;
    }
    return kTemplates[0];   // unreachable while every class has a row
}

// A zone kept whole, for what no template claimed.
TMeshBlock zoneAsBlock(const MedialAxis &axis, const MedialZone &zone) {
    TMeshBlock b;
    b.outline = zone.boundarySide;
    std::vector<Point> chainPts;
    for (int m : zone.chain) chainPts.push_back(axis.medialVertices[m].coord);
    if (zone.cap) {
        // Open caps close through their vertex; a cap spanning a whole closed
        // loop is already a cycle and needs no spoke at all.
        if (normP(zone.boundarySide.front() - zone.boundarySide.back()) >
            1e-9 * axis.boundingBoxDiagonal()) {
            b.outline.push_back(chainPts.front());
        }
        b.corners = {zone.boundarySide.front(), zone.boundarySide.back(),
                     chainPts.front()};
    } else {
        appendPts(b.outline, reversed(chainPts));
        b.corners = {zone.boundarySide.front(), zone.boundarySide.back(),
                     chainPts.back(), chainPts.front()};
    }
    b.color = zone.color;
    return b;
}

// Orient the two flanking zones of `edge` into the frame `t` expects.
TemplateSubdomain makeSubdomain(const MedialAxis &axis,
                                const std::vector<MedialZone> &zones,
                                const std::unordered_map<int, int> &capOf,
                                const CoarseEdge &edge,
                                const TemplateEntry &t) {
    const MedialZone &A = zones[edge.zones[0]];
    const MedialZone &B = zones[edge.zones[1]];

    // The two zones walk the shared chain in opposite directions, so side B
    // comes onto side A's orientation before anything else.
    std::vector<int> chainIdx = A.chain;
    std::vector<Point> sideA = A.boundarySide;
    std::vector<Point> sideB =
        (B.chain == A.chain) ? B.boundarySide : reversed(B.boundarySide);

    bool innerAtFront = true;
    if (t.orient == Orient::Narrow) {
        innerAtFront = axis.medialVertices[chainIdx.front()].radius <=
                       axis.medialVertices[chainIdx.back()].radius;
    } else if (t.orient == Orient::Endpoint) {
        const bool frontIsEnd = axis.medialVertices[chainIdx.front()].degree <= 1;
        const bool backIsEnd  = axis.medialVertices[chainIdx.back()].degree <= 1;
        // Both ends can terminate the axis (a lone branch), so prefer the one
        // that actually carries a cap for the template to swallow.
        if (frontIsEnd && backIsEnd) {
            innerAtFront = capOf.count(chainIdx.front()) ||
                           !capOf.count(chainIdx.back());
        } else {
            innerAtFront = frontIsEnd;
        }
    }
    if (!innerAtFront) {
        std::reverse(chainIdx.begin(), chainIdx.end());
        std::reverse(sideA.begin(), sideA.end());
        std::reverse(sideB.begin(), sideB.end());
    }

    TemplateSubdomain s;
    s.chain.reserve(chainIdx.size());
    for (int m : chainIdx) s.chain.push_back(axis.medialVertices[m].coord);
    s.inner = s.chain.front();
    s.outer = s.chain.back();
    s.sideA = std::move(sideA);
    s.sideB = std::move(sideB);
    if (!t.wantsCap) return s;

    // The polar cap at `inner`, oriented to leave side A and arrive at side B.
    std::vector<Point> cap;
    const auto capIt = capOf.find(chainIdx.front());
    if (capIt != capOf.end()) {
        s.capZone = capIt->second;
        const std::vector<Point> &c = zones[s.capZone].boundarySide;
        const bool forward = normP(c.front() - s.sideA.front()) <=
                             normP(c.back() - s.sideA.front());
        cap = forward ? c : reversed(c);
    }

    s.arcUp = reversed(s.sideA);                 // outer's contact -> inner's
    if (cap.size() >= 2) {
        std::vector<std::vector<Point>> pieces;
        std::vector<Point> cuts;
        splitPolyline(cap, {0.5}, pieces, cuts);
        appendPts(s.arcUp, pieces[0]);           // ... -> W
        s.arcLo = pieces[1];                     // W -> ...
    } else {
        // No cap: the two inner contacts coincide, and that point is W.
        s.arcLo = {s.sideA.front()};
    }
    appendPts(s.arcLo, s.sideB);
    return s;
}

} // namespace

void MedialAxisTMesh::buildBlocks() {
    blocks.clear();

    // Caps by endpoint, so a template that wants one can find it.
    std::unordered_map<int, int> capOf;
    for (int z = 0; z < static_cast<int>(zones.size()); ++z) {
        if (zones[z].cap) capOf[zones[z].capVertex] = z;
    }
    std::vector<char> consumed(zones.size(), 0);
    const double weld = 1e-12 * axis->boundingBoxDiagonal();

    for (const CoarseEdge &edge : edges) {
        // Without both flanks the subdomain is incomplete; keep what exists.
        if (edge.zones[1] < 0) {
            if (edge.zones[0] >= 0) {
                blocks.push_back(zoneAsBlock(*axis, zones[edge.zones[0]]));
                consumed[edge.zones[0]] = 1;
            }
            continue;
        }
        consumed[edge.zones[0]] = 1;
        consumed[edge.zones[1]] = 1;

        const TemplateEntry &t = templateFor(edge.color);
        const TemplateSubdomain s = makeSubdomain(*axis, zones, capOf, edge, t);
        if (s.capZone >= 0) consumed[s.capZone] = 1;

        BlockBuilder builder(blocks, edge.color, weld);
        t.apply(s, builder);
    }

    // Whatever the templates did not absorb -- standalone caps, zones of
    // unpaired edges -- still tiles part of the domain, so it stays a block.
    for (int z = 0; z < static_cast<int>(zones.size()); ++z) {
        if (!consumed[z]) blocks.push_back(zoneAsBlock(*axis, zones[z]));
    }
}
