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

void MedialAxisTMesh::buildBlocks() {
    blocks.clear();

    // Caps by endpoint, so the blue and purple templates can absorb the polar
    // region at their endpoint into the subdomain they remesh.
    std::unordered_map<int, int> capOf;
    for (int z = 0; z < static_cast<int>(zones.size()); ++z) {
        if (zones[z].cap) capOf[zones[z].capVertex] = z;
    }
    std::vector<char> consumed(zones.size(), 0);

    auto zoneAsBlock = [&](const MedialZone &zone) {
        TMeshBlock b;
        b.outline = zone.boundarySide;
        std::vector<Point> chainPts;
        for (int m : zone.chain) chainPts.push_back(axis->medialVertices[m].coord);
        if (zone.cap) {
            // Open caps close through their vertex; a cap spanning a whole
            // closed loop is already a cycle and needs no spoke at all.
            if (normP(zone.boundarySide.front() - zone.boundarySide.back()) >
                1e-9 * axis->boundingBoxDiagonal()) {
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
    };

    for (const CoarseEdge &edge : edges) {
        // Without both flanks the subdomain is incomplete; keep what exists.
        if (edge.zones[1] < 0) {
            if (edge.zones[0] >= 0) {
                blocks.push_back(zoneAsBlock(zones[edge.zones[0]]));
                consumed[edge.zones[0]] = 1;
            }
            continue;
        }
        const MedialZone &A = zones[edge.zones[0]];
        const MedialZone &B = zones[edge.zones[1]];
        consumed[edge.zones[0]] = 1;
        consumed[edge.zones[1]] = 1;

        if (edge.color == MedialColor::Green) {
            // Two quads, one per side: the zones already are the template.
            blocks.push_back(zoneAsBlock(A));
            blocks.push_back(zoneAsBlock(B));
            continue;
        }

        // ── A common frame for the directional templates ──
        //
        // chain / sideA / sideB oriented so that the template's special end
        // (the narrow end for red, the endpoint for blue and purple) is at
        // the front. Side B is stored against its own walk direction, i.e.
        // reversed relative to A, unless the pairing matched equal chains.
        std::vector<int> chainIdx = A.chain;
        std::vector<Point> sideA = A.boundarySide;
        std::vector<Point> sideB =
            (B.chain == A.chain) ? B.boundarySide : reversed(B.boundarySide);

        bool specialAtFront;
        if (edge.color == MedialColor::Red) {
            specialAtFront = axis->medialVertices[chainIdx.front()].radius <=
                             axis->medialVertices[chainIdx.back()].radius;
        } else {
            const bool frontIsEnd =
                axis->medialVertices[chainIdx.front()].degree <= 1;
            const bool backIsEnd =
                axis->medialVertices[chainIdx.back()].degree <= 1;
            // Both ends can be endpoints (a single-branch domain); prefer the
            // one whose cap exists so the template absorbs it.
            if (frontIsEnd && backIsEnd)
                specialAtFront = capOf.count(chainIdx.front()) || !capOf.count(chainIdx.back());
            else
                specialAtFront = frontIsEnd;
        }
        if (!specialAtFront) {
            std::reverse(chainIdx.begin(), chainIdx.end());
            std::reverse(sideA.begin(), sideA.end());
            std::reverse(sideB.begin(), sideB.end());
        }
        std::vector<Point> chainPts;
        for (int m : chainIdx) chainPts.push_back(axis->medialVertices[m].coord);

        const Point L = chainPts.front();
        const Point R = chainPts.back();

        if (edge.color == MedialColor::Red) {
            // Three quads around the concave corner (Fig. 17, column 2).
            //
            // At the narrow medial vertex L both radii point back toward the
            // apex, so the subdomain hexagon is reflex there -- L is its one
            // concave corner. It is split by running an edge from L to each
            // of the two opposite sides, which here are the radii at the wide
            // vertex R. Those edges land mid-side, so each leaves a
            // T-junction: the neighbouring subdomain across that radius has
            // no vertex there. The medial edge L-R is left as the diagonal of
            // the middle quad rather than an edge of it, which is exactly
            // what red asks for -- the axis follows quad diagonals.
            const Point A = (sideA.back() + R) * 0.5;   // T-junction on R's radius
            const Point B = (sideB.back() + R) * 0.5;

            TMeshBlock up;
            up.outline = {L};
            appendPts(up.outline, sideA);
            appendPts(up.outline, {A});
            up.corners = {L, sideA.front(), sideA.back(), A};
            up.color = edge.color;
            blocks.push_back(std::move(up));

            TMeshBlock mid;
            mid.outline = {L, A, R, B};
            mid.corners = mid.outline;
            mid.color = edge.color;
            blocks.push_back(std::move(mid));

            TMeshBlock lo;
            lo.outline = {L, B};
            appendPts(lo.outline, reversed(sideB));
            lo.corners = {L, B, sideB.back(), sideB.front()};
            lo.color = edge.color;
            blocks.push_back(std::move(lo));
            continue;
        }

        // Blue and purple own the polar cap at their endpoint, oriented to
        // run from side A's inner corner around the tip to side B's.
        std::vector<Point> cap;
        const auto capIt = capOf.find(chainIdx.front());
        if (capIt != capOf.end()) {
            consumed[capIt->second] = 1;
            const std::vector<Point> &c = zones[capIt->second].boundarySide;
            const bool forward = normP(c.front() - sideA.front()) <=
                                 normP(c.back() - sideA.front());
            cap = forward ? c : reversed(c);
        }

        if (edge.color == MedialColor::Blue) {
            // One quad spanning the whole wedge, the axis on its diagonal:
            // the far vertex, its two contacts, and the tip between them.
            std::vector<Point> arc = reversed(sideA);
            appendPts(arc, cap);
            appendPts(arc, sideB);

            TMeshBlock quad;
            quad.outline = {R};
            appendPts(quad.outline, arc);
            quad.corners = {R, sideA.back(),
                            cap.empty() ? sideA.front() : cap[cap.size() / 2],
                            sideB.back()};
            quad.color = edge.color;
            blocks.push_back(std::move(quad));
            continue;
        }

        // Purple: three quads, the same construction as red (Fig. 17, col. 4).
        //
        // The T-junctions sit at the midpoints of the two radii at R, and the
        // half-radius from each one out to the boundary is a template edge,
        // exactly as in the red case. What differs is the far end: an
        // endpoint has no pair of inner radii to act as the corners of a
        // concave vertex, so its polar cap stands in for them, split at its
        // midpoint W -- the single boundary point the endpoint maps to. The
        // medial edge L-R is again the diagonal of the middle quad.
        const Point Tu = (sideA.back() + R) * 0.5;   // T-junction on R's radius
        const Point Tl = (sideB.back() + R) * 0.5;

        // The subdomain's boundary, from R's contact on side A around the
        // tip to its contact on side B, cut at W.
        std::vector<Point> arcUp = reversed(sideA);   // uR -> gA
        std::vector<Point> arcLo;                     // W  -> lR
        if (cap.size() >= 2) {
            std::vector<std::vector<Point>> capPieces;
            std::vector<Point> capCuts;
            splitPolyline(cap, {0.5}, capPieces, capCuts);
            appendPts(arcUp, capPieces[0]);           // gA -> W
            arcLo = capPieces[1];                     // W  -> gB
        } else {
            // No cap: the two inner contacts coincide and W is that point.
            arcLo = {sideA.front()};
        }
        appendPts(arcLo, sideB);                      // gB -> lR

        TMeshBlock up;
        up.outline = arcUp;
        appendPts(up.outline, {L, Tu});
        up.corners = {sideA.back(), arcUp.back(), L, Tu};
        up.color = edge.color;
        blocks.push_back(std::move(up));

        TMeshBlock lo;
        lo.outline = arcLo;
        appendPts(lo.outline, {Tl, L});
        lo.corners = {arcLo.front(), sideB.back(), Tl, L};
        lo.color = edge.color;
        blocks.push_back(std::move(lo));

        TMeshBlock diamond;
        diamond.outline = {L, Tu, R, Tl};
        diamond.corners = diamond.outline;
        diamond.color = edge.color;
        blocks.push_back(std::move(diamond));
    }

    // Whatever the templates did not absorb -- standalone caps, zones of
    // unpaired edges -- still tiles part of the domain, so it stays a block.
    for (int z = 0; z < static_cast<int>(zones.size()); ++z) {
        if (!consumed[z]) blocks.push_back(zoneAsBlock(zones[z]));
    }
}
