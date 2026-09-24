#include "tracing/SeparatrixTrace.hxx"

#include "MERIDIAN/Interfaces.hxx"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <unordered_map>

namespace {

// Intersection of segments p0->p1 and q0->q1, counted on the half-open
// convention: a meeting at the far end of either segment belongs to the next
// segment along, not to this one.
//
// Half-open rather than strictly interior because a streamline that passes
// through a mesh vertex is snapped onto it, so two of them can share a polyline
// vertex exactly. Excluding both ends there loses the crossing twice over --
// once at t = 1 on the way in and once at u = 1 -- and it reappears as a place
// two arcs of the finished layout cross with no node between them. The one
// meeting this convention does count that should not is two separatrices
// leaving the same singularity, which share their very first point; the caller
// excludes that pair by origin.
bool properIntersection(const Point &p0, const Point &p1, const Point &q0, const Point &q1,
                        double &t, double &u, Point &at) {
    const Point r = p1 - p0;
    const Point s = q1 - q0;
    const double den = cross2(r, s);
    // Parallel to within rounding, relative to the lengths: see
    // QuadLayout::checkEmbedding for what an absolute threshold lets through.
    if (std::fabs(den) <= 1e-12 * normP(r) * normP(s)) return false;

    const Point w = q0 - p0;
    t = cross2(w, s) / den;
    u = cross2(w, r) / den;
    const double eps = 1e-9;
    if (t < -eps || t >= 1.0 - eps) return false;
    if (u < -eps || u >= 1.0 - eps) return false;
    at = p0 + r * t;
    return true;
}

} // namespace

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
SeparatrixTrace::SeparatrixTrace(std::shared_ptr<CrossField> cf,
                                 bool useActualSingularityCoordinates, const Settings &s,
                                 const Interfaces *itf)
    : crossField(cf), settings(s) {
    // The interface network's nodes are the layout's singular points wherever
    // interfaces meet, land or kink -- they have their own sectors, and the
    // field's singular triangles there are theirs, not singularities to trace
    // from (FieldTracer's `absorbAt`).
    if (itf && settings.respectInterfaces && itf->multiMaterial()) interfaces = itf;
    std::vector<char> absorbAt;
    if (interfaces) {
        interfaceNodeVertex.assign(cf->mesh->vertices.size(), 0);
        for (const Interfaces::Node &n : interfaces->nodes()) {
            if (n.kind == Interfaces::NodeKind::Dangling || n.vertex < 0) continue;
            interfaceNodeVertex[n.vertex] = 1;
        }
        absorbAt = interfaceNodeVertex;
    }
    tracer = std::make_unique<FieldTracer>(cf, useActualSingularityCoordinates,
                                           settings.singularitiesAtBoundary, absorbAt);
    mesh = &tracer->getMesh();

    segmentsOfTriangle.assign(mesh->triangles.size(), {});

    buildBoundary();
    launchFromSingularities();
    launchFromCorners();
    launchFromInterfaceNodes();

    singularities = tracer->getSingularities();

    cutRadius.assign(singularities.size(), 0.0);
    for (size_t i = 0; i < singularities.size(); ++i) {
        const Triangle &tri = mesh->triangles[singularities[i].triangleIndex];
        double h = 0.0;
        for (int k = 0; k < 3; ++k)
            h += normP(mesh->vertices[tri[(k + 1) % 3]] - mesh->vertices[tri[k]]);
        cutRadius[i] = settings.singularityCutRadius * h / 3.0;
    }

    arcLength.assign(separatrices.size(), 0.0);

    // Poincare-Hopf, before anything is traced: a failure here is a missed
    // singularity or a misread corner, and it is worth knowing about before
    // blaming the tracing for the components it cannot close.
    //
    // With interfaces, a node's index is the one its sectors prescribe
    // (Interfaces::Node::index) and stands in both for the field's singular
    // triangles absorbed there and, at a landing, for the corner reading of
    // the boundary vertex it sits on.
    report_.eulerCharacteristic = 2 - static_cast<int>(boundaryLoops.size());
    report_.interiorIndexSum4 = tracer->interiorIndexSum4();
    for (const BoundaryCorner &c : boundaryCorners) {
        if (interfaces && interfaceNodeVertex[c.vertex]) continue;
        report_.boundaryIndexSum4 += 2 - c.quarters;
    }
    for (const InterfaceNodeInfo &n : interfaceNodes) report_.interfaceIndexSum4 += n.index;
    report_.interfaceNodes = static_cast<int>(interfaceNodes.size());
    report_.absorbedSingularities = tracer->absorbedSingularityCount();
    report_.poincareHopf = report_.interiorIndexSum4 + report_.boundaryIndexSum4 +
                               report_.interfaceIndexSum4 ==
                           4 * report_.eulerCharacteristic;
    report_.multipleSingularities = tracer->multipleSingularityCount();
    report_.droppedSingularities = tracer->droppedSingularityCount();

    finishedTracing = separatrices.empty();
}

// ---------------------------------------------------------------------------
// buildBoundary()  --  loops, interior angles, corners
//
// A boundary edge is walked in the orientation its own triangle gives it, so
// each loop comes out with the interior on its left without anything having to
// be tested for being the outer one.
// ---------------------------------------------------------------------------
void SeparatrixTrace::buildBoundary() {
    const int nV = static_cast<int>(mesh->vertices.size());

    // Directed boundary edges, interior on the left.
    std::unordered_map<int, int> nextOf;   // vertex -> vertex
    nextOf.reserve(mesh->boundaryEdges.size() * 2);
    for (const int e : mesh->boundaryEdges) {
        const int f = (mesh->edgeTriangles[e][0] >= 0) ? mesh->edgeTriangles[e][0]
                                                       : mesh->edgeTriangles[e][1];
        if (f < 0) continue;
        for (int i = 0; i < 3; ++i) {
            if (mesh->triangleEdges[f][i] != e) continue;
            nextOf[mesh->triangles[f][i]] = mesh->triangles[f][(i + 1) % 3];
            break;
        }
    }

    std::vector<char> seen(nV, 0);
    for (const auto &kv : nextOf) {
        if (seen[kv.first]) continue;
        std::vector<int> loop;
        int v = kv.first;
        while (!seen[v]) {
            seen[v] = 1;
            loop.push_back(v);
            auto it = nextOf.find(v);
            if (it == nextOf.end()) { loop.clear(); break; }
            v = it->second;
        }
        if (loop.size() >= 3) boundaryLoops.push_back(std::move(loop));
    }

    // Interior angle at each boundary vertex, as the sum of the tip angles of
    // the triangles around it.
    std::vector<double> interior(nV, 0.0);
    for (int f = 0; f < static_cast<int>(mesh->triangles.size()); ++f) {
        const Triangle &tri = mesh->triangles[f];
        for (int i = 0; i < 3; ++i) {
            const int v = tri[i];
            if (!mesh->isBoundaryVertex[v]) continue;
            const Point a = mesh->vertices[tri[(i + 1) % 3]] - mesh->vertices[v];
            const Point b = mesh->vertices[tri[(i + 2) % 3]] - mesh->vertices[v];
            const double na = normP(a), nb = normP(b);
            if (na < 1e-18 || nb < 1e-18) continue;
            interior[v] += std::acos(std::max(-1.0, std::min(1.0, dotP(a, b) / (na * nb))));
        }
    }

    // Table 1: the corner index is the interior angle rounded to a multiple of
    // pi/2. A vertex that rounds to pi is not a corner at all and carries no
    // node; every other one does, whether or not it launches anything.
    //
    // Read with Table 1's own intervals -- (0, 3pi/4), [3pi/4, 5pi/4],
    // (5pi/4, 7pi/4], (7pi/4, 2pi) -- and not by rounding to the nearest
    // quarter, because CrossField::initialize sets the boundary data with
    // exactly these, and the corners the layout turns at have to be the ones
    // the field was told to turn at. The two readings differ at 225 and 315
    // degrees, where rounding goes up and Table 1 goes down: multimat/rocket
    // has five such corners, and read the other way its Poincare-Hopf count
    // did not close.
    auto quartersOf = [](double a) {
        if (a < 0.75 * M_PI) return 1;
        if (a <= 1.25 * M_PI) return 2;
        if (a <= 1.75 * M_PI) return 3;
        return std::max(4, static_cast<int>(std::lround(a / M_PI_2)));
    };
    for (const auto &loop : boundaryLoops) {
        for (const int v : loop) {
            const int q = quartersOf(interior[v]);
            if (q == 2) continue;
            BoundaryCorner c;
            c.vertex = v;
            c.quarters = q;
            c.interiorAngle = interior[v];
            boundaryCorners.push_back(std::move(c));
        }
    }
}

// ---------------------------------------------------------------------------
// Launching
// ---------------------------------------------------------------------------
void SeparatrixTrace::launchFromSingularities() {
    auto &sings = tracer->getSingularities();
    for (int si = 0; si < static_cast<int>(sings.size()); ++si) {
        Singularity &s = sings[si];
        for (int port = 0; port < s.numPorts(); ++port) {
            Separatrix sep;
            sep.id = static_cast<int>(separatrices.size());
            sep.originKind = SeparatrixOrigin::Singularity;
            sep.origin_singularity_id = si;
            sep.origin_singularity_port = port;

            TracePoint tp;
            tp.global_pos = s.coordinates;
            tp.face_id = s.triangleIndex;
            tp.theta = s.portAngles[port];
            sep.path.push_back(tp);

            s.portSeparatrixIds[port] = sep.id;
            walkers.push_back(tracer->startAt(s.triangleIndex, s.coordinates, s.portAngles[port]));
            stepsTaken.push_back(0);
            crossCount.emplace_back();
            separatrices.push_back(std::move(sep));
        }
    }
}

int SeparatrixTrace::launchFromVertex(int v, double dir, SeparatrixOrigin kind, int origin,
                                      int port) {
    const auto range = mesh->vertexTriangles.trianglesForVertex(v);
    if (range.first == range.second) return -1;
    dir = wrap_pi(dir);

    Separatrix sep;
    sep.id = static_cast<int>(separatrices.size());
    sep.originKind = kind;
    sep.origin_singularity_id = origin;
    sep.origin_singularity_port = port;

    Walker w;
    w.tri = *range.first;
    w.pos = mesh->vertices[v];
    w.entryEdge = -1;
    w.atVertex = v;
    w.dir = dir;
    w.crossDir = dir;

    TracePoint tp;
    tp.global_pos = w.pos;
    tp.face_id = w.tri;
    tp.theta = dir;
    sep.path.push_back(tp);

    walkers.push_back(w);
    stepsTaken.push_back(0);
    crossCount.emplace_back();
    separatrices.push_back(std::move(sep));
    return static_cast<int>(separatrices.size()) - 1;
}

void SeparatrixTrace::launchFromCorners() {
    // Direction from a boundary vertex to the next one round the loop; the
    // interior lies to its left, so the axis directions that point inwards are
    // at +pi/2, +pi, ... from it.
    std::unordered_map<int, int> nextOnLoop;
    for (const auto &loop : boundaryLoops)
        for (size_t i = 0; i < loop.size(); ++i)
            nextOnLoop[loop[i]] = loop[(i + 1) % loop.size()];

    for (int ci = 0; ci < static_cast<int>(boundaryCorners.size()); ++ci) {
        BoundaryCorner &c = boundaryCorners[ci];
        const int rays = c.quarters - 1;
        if (rays <= 0) continue;
        // An interface landing here cuts the corner's wedge into sectors of
        // its own, and launchFromInterfaceNodes() launches into those.
        if (interfaces && interfaceNodeVertex[c.vertex]) continue;

        auto it = nextOnLoop.find(c.vertex);
        if (it == nextOnLoop.end()) continue;
        const Point along = mesh->vertices[it->second] - mesh->vertices[c.vertex];
        if (normP(along) < 1e-18) continue;
        const double base = computeAngle(along);

        const double step = settings.evenCornerRays ? (c.interiorAngle / c.quarters) : M_PI_2;
        for (int j = 1; j <= rays; ++j) {
            const int id = launchFromVertex(c.vertex, base + j * step,
                                            SeparatrixOrigin::BoundaryCorner, ci, j - 1);
            if (id >= 0) c.separatrixIds.push_back(id);
        }
    }
}

// ---------------------------------------------------------------------------
// launchFromInterfaceNodes()  --  the corners the materials make
//
// A node of the interface network is a boundary corner of every region it
// touches: its rays (the interface branches, and at a landing the two edges of
// dS) cut its fan into sectors, and a sector of q quarter turns (Interfaces::
// quantiseNode) wants q - 1 separatrices in it, exactly as a corner of dS does.
// MERIDIAN's Stage 7 emits the same count from the same nodes
// (Separatrices::Options::extraEmitters); the rays running along the
// interfaces themselves are never emitted, since the interface is already
// that edge of the layout. Each sector is divided evenly, which is what
// Settings::evenCornerRays does at a corner of dS, and for the same reason:
// the ray only chooses the branch, and the field is followed from the first
// triangle on.
// ---------------------------------------------------------------------------
void SeparatrixTrace::launchFromInterfaceNodes() {
    if (!interfaces) return;
    std::unordered_map<int, int> cornerQuarters;   // boundary vertex -> q, where q != 2
    for (const BoundaryCorner &c : boundaryCorners) cornerQuarters[c.vertex] = c.quarters;

    const auto &nodes = interfaces->nodes();
    for (int ni = 0; ni < static_cast<int>(nodes.size()); ++ni) {
        const Interfaces::Node &n = nodes[ni];
        if (n.kind == Interfaces::NodeKind::Dangling || n.vertex < 0) continue;
        InterfaceNodeInfo info;
        info.vertex = n.vertex;
        info.node = ni;
        info.onBoundary = n.onBoundary;
        info.quarters = n.quarters;

        // How many quarter turns the node's fan makes in the layout is the
        // field's to say, not the geometry's: it is the field the separatrices
        // follow, and a node that turns a quarter more or less than the field
        // round it leaves a component with a corner too many or too few. The
        // field's count is the fan's own -- 4 inside, the corner reading of
        // Table 1 on dS (buildBoundary), the same boundary data the field was
        // solved with -- less the index of the singular triangles the node
        // absorbed, which is where the field turns by more or less than that.
        //
        // The geometry says only how the count is shared out: by largest
        // remainder over the sectors' own quarter turns, each keeping at least
        // one. Interfaces::quantiseNode rounds each sector on its own instead,
        // and the two disagree exactly where it matters -- an interface
        // bisecting a 270-degree corner leaves two 135-degree sectors that
        // round to two apiece, four in a corner of three. MERIDIAN settles
        // that region by region (Interfaces::balance); here the field does.
        if (!n.sector.empty()) {
            int fan = 4;
            if (n.onBoundary) {
                auto it = cornerQuarters.find(n.vertex);
                fan = (it != cornerQuarters.end()) ? it->second : 2;
            }
            const int total = std::max(static_cast<int>(n.sector.size()),
                                       fan - tracer->absorbedIndexAt(n.vertex));
            std::vector<int> q(n.sector.size(), 1);
            int used = static_cast<int>(q.size());
            std::vector<std::pair<double, int>> want;
            for (size_t k = 0; k < n.sector.size(); ++k) want.push_back({n.sector[k] / M_PI_2, static_cast<int>(k)});
            while (used < total) {
                int best = -1;
                double bestGap = -1e300;
                for (const auto &[w, k] : want) {
                    const double gap = w - q[k];
                    if (gap > bestGap) { bestGap = gap; best = k; }
                }
                ++q[best];
                ++used;
            }
            info.quarters = q;
        }
        int total = 0;
        for (const int q : info.quarters) total += q;
        info.index = (n.onBoundary ? 2 : 4) - total;

        const int slot = static_cast<int>(interfaceNodes.size());
        int port = 0;
        for (size_t k = 0; k < n.sector.size() && k < info.quarters.size(); ++k) {
            const int q = info.quarters[k];
            if (q < 2 || k >= n.rays.size()) continue;
            for (int j = 1; j < q; ++j) {
                const int id = launchFromVertex(n.vertex, n.rays[k].dir + j * n.sector[k] / q,
                                                SeparatrixOrigin::InterfaceNode, slot, port++);
                if (id >= 0) {
                    info.separatrixIds.push_back(id);
                    ++report_.interfaceEmitted;
                }
            }
        }
        interfaceNodes.push_back(std::move(info));
    }
}

// ---------------------------------------------------------------------------
// recordInterfaceCrossing()
//
// A step ends on the edge (or vertex) it hands the walk across, so a
// separatrix crosses an interface exactly when the triangle it just crossed
// and the one it is about to cross are different materials, and the point
// where it does is the last point of its path.
// ---------------------------------------------------------------------------
void SeparatrixTrace::recordInterfaceCrossing(int k) {
    const Separatrix &sep = separatrices[k];
    const Walker &w = walkers[k];
    if (sep.path.size() < 2 || w.tri < 0) return;
    const int f0 = sep.path.back().face_id;
    if (f0 < 0 || f0 == w.tri) return;
    if (mesh->triangleMatId[f0] == mesh->triangleMatId[w.tri]) return;

    InterfaceCrossing x;
    x.sep = k;
    x.pathIndex = static_cast<int>(sep.path.size()) - 1;
    x.pos = sep.path.back().global_pos;
    Point along{0.0, 0.0};
    if (w.entryEdge >= 0) {
        x.edge = mesh->triangleEdges[w.tri][w.entryEdge];
        along = mesh->vertices[mesh->edges[x.edge][1]] - mesh->vertices[mesh->edges[x.edge][0]];
    } else {
        double best = std::numeric_limits<double>::max();
        for (int i = 0; i < 3; ++i) {
            const int v = mesh->triangles[w.tri][i];
            const double d = normP(mesh->vertices[v] - x.pos);
            if (d < best) { best = d; x.vertex = v; }
        }
        if (best > 1e-9 * tracer->averageEdgeLength()) return;
    }
    const Point in = x.pos - sep.path[sep.path.size() - 2].global_pos;
    if (normP(along) > 0.0 && normP(in) > 0.0) {
        const double c = std::fabs(dotP(along, in)) / (normP(along) * normP(in));
        x.angle = std::acos(std::min(1.0, c));
        if (x.angle < settings.tangentialAngle) ++report_.interfaceGrazes;
    } else {
        x.angle = M_PI_2;
    }
    interfaceCrossings.push_back(x);
    ++report_.interfaceCrossings;
}

// ---------------------------------------------------------------------------
// Stepping
// ---------------------------------------------------------------------------
int SeparatrixTrace::countActive() const {
    int n = 0;
    for (const auto &s : separatrices) if (s.active) ++n;
    return n;
}

void SeparatrixTrace::terminateAt(Separatrix &sep, int segIndex, const Point &at,
                                  TerminationReason why, int onSep, int onSeg) {
    sep.path.resize(segIndex + 1);
    sep.path.back().global_pos = at;
    sep.active = false;
    sep.termination_reason = why;
    sep.endOnSeparatrix = onSep;
    sep.endOnSegment = onSeg;
    double L = 0.0;
    for (size_t i = 1; i < sep.path.size(); ++i)
        L += normP(sep.path[i].global_pos - sep.path[i - 1].global_pos);
    if (sep.id >= 0 && sep.id < static_cast<int>(arcLength.size())) arcLength[sep.id] = L;
}

// ---------------------------------------------------------------------------
// retractTail()  --  what a heteroclinic join takes away from the other curve
//
// The join cuts `host` back to the meeting point, so every event recorded on
// the part beyond it is now an event on nothing. A crossing there is simply
// withdrawn, with its count. A separatrix that had *stopped* there -- on its
// second crossing of the host, or on meeting it tangentially -- has lost the
// thing it stopped against and would be a loose end of the layout, so it is
// resumed from where it stopped (docs/viertel_2019.md Sec. 4.4). The spec's
// remark that the removed tail is short because both fronts arrive at about
// the same time holds because every separatrix is grown at once, and is why
// this is rare -- on data/meshes it never happens -- not why it is unnecessary.
//
// One case needs care: a separatrix stopped by condition 2 whose *first*
// crossing of the host was on the removed part. Its stop is still on the host
// but is now its only crossing of it, so condition 2 no longer applies. It is
// resumed too, with the count for the stopping crossing withdrawn; the first
// segment it traces starts on the host, and the half-open intersection test
// counts that as the crossing it now is.
// ---------------------------------------------------------------------------
bool SeparatrixTrace::survivesTruncation(int host, int seg, const Point &q, int cutSeg,
                                         const Point &cut) const {
    if (seg < cutSeg) return true;
    if (seg > cutSeg) return false;
    const Point &a = separatrices[host].path[cutSeg - 1].global_pos;
    return normP(q - a) <= normP(cut - a) + 1e-12 * tracer->averageEdgeLength();
}

void SeparatrixTrace::retractTail(int host, int seg, const Point &at) {
    std::vector<SeparatrixCrossing> kept;
    kept.reserve(crossings.size());
    for (const SeparatrixCrossing &x : crossings) {
        const bool gone = (x.sepA == host && !survivesTruncation(host, x.segA, x.pos, seg, at)) ||
                          (x.sepB == host && !survivesTruncation(host, x.segB, x.pos, seg, at));
        if (!gone) { kept.push_back(x); continue; }
        --crossCount[x.sepA][x.sepB];
        --crossCount[x.sepB][x.sepA];
        ++report_.crossingsRetracted;
    }
    crossings.swap(kept);

    for (Separatrix &z : separatrices) {
        if (z.active || z.id == host || z.endOnSeparatrix != host || z.path.size() < 2) continue;
        const bool twice = z.termination_reason == TerminationReason::CROSSED_TWICE;
        if (!twice && z.termination_reason != TerminationReason::TANGENTIAL &&
            z.termination_reason != TerminationReason::CUT_AT_SINGULARITY)
            continue;
        const bool stopSurvives =
            survivesTruncation(host, z.endOnSegment, z.path.back().global_pos, seg, at);
        // The stopping crossing of condition 2 was counted when it happened.
        const bool stillTwice = twice && crossCount[z.id][host] >= 2;
        if (stopSurvives && (!twice || stillTwice)) continue;
        if (twice) {
            --crossCount[z.id][host];
            --crossCount[host][z.id];
        }
        z.active = true;
        z.termination_reason = TerminationReason::RUNNING;
        z.endOnSeparatrix = -1;
        z.endOnSegment = -1;
        const TracePoint &tp = z.path.back();
        walkers[z.id] = tracer->startAt(tp.face_id, tp.global_pos, tp.theta);
        resumeQueue.push_back(z.id);
        ++report_.resumed;
    }
}

// ---------------------------------------------------------------------------
// registerNewSegments()  --  what the last step ran into
//
// Every segment lies inside one triangle, so two of them can only meet if they
// share a triangle, and the whole crossing test is against the handful of
// segments already filed under that triangle.
// ---------------------------------------------------------------------------
bool SeparatrixTrace::registerNewSegments(Separatrix &sep, int firstNewIndex) {
    for (int i = firstNewIndex; i < static_cast<int>(sep.path.size()); ++i) {
        const int f = sep.path[i].face_id;
        if (f < 0 || f >= static_cast<int>(segmentsOfTriangle.size())) continue;

        const Point &a0 = sep.path[i - 1].global_pos;
        const Point &a1 = sep.path[i].global_pos;

        // Everything this segment meets, in the order it meets it. Acting on
        // each as it is found instead would let the order the triangle's
        // segments happen to be filed in decide the outcome: a stopping
        // condition found early returns, and a crossing further along the same
        // segment -- which is still there, on the part of it that survives --
        // is never looked at and ends up in the layout with no node on it.
        struct Hit { double t; SegRef other; Point at; };
        std::vector<Hit> hits;
        for (const SegRef &other : segmentsOfTriangle[f]) {
            // Consecutive segments of one separatrix share an endpoint by
            // construction; a separatrix crossing itself further along is a
            // genuine event and is left in.
            if (other.sep == sep.id && std::abs(other.seg - i) <= 1) continue;

            const Separatrix &osep = separatrices[other.sep];
            // Two separatrices out of one singularity, or one boundary corner,
            // start at the same point: their first segments meet there and it
            // is not a crossing.
            if (i == 1 && other.seg == 1 && osep.originKind == sep.originKind &&
                osep.origin_singularity_id == sep.origin_singularity_id)
                continue;
            if (other.seg <= 0 || other.seg >= static_cast<int>(osep.path.size())) continue;

            double t = 0.0, u = 0.0;
            Point at{0.0, 0.0};
            if (properIntersection(a0, a1, osep.path[other.seg - 1].global_pos,
                                   osep.path[other.seg].global_pos, t, u, at))
                hits.push_back(Hit{t, other, at});
        }
        std::sort(hits.begin(), hits.end(),
                  [](const Hit &x, const Hit &y) { return x.t < y.t; });

        for (const Hit &h : hits) {
            const Separatrix &osep = separatrices[h.other.sep];
            const Point &b0 = osep.path[h.other.seg - 1].global_pos;
            const Point &b1 = osep.path[h.other.seg].global_pos;

            // Sec. 3.3, first case: two separatrices crossing at a shallow
            // angle head-on are one curve the discretisation has split in two,
            // and cutting both at the meeting joins them into the heteroclinic
            // orbit their continuous counterparts form. The layout needs
            // nothing more than that -- both ending at the same point makes one
            // node of valence two, which a face walks straight through.
            //
            // Measured on the chords rather than on the field angles: what
            // decides whether the join reads as one curve is the direction the
            // two polylines actually run in at the meeting.
            const double phi = std::fabs(wrap_pi(computeAngle(a1 - a0) - computeAngle(b1 - b0)));
            if (phi > M_PI - settings.tangentialAngle && osep.active && h.other.sep != sep.id) {
                const int host = h.other.sep, hostSeg = h.other.seg;
                terminateAt(sep, i, h.at, TerminationReason::HETEROCLINIC, host, hostSeg);
                terminateAt(separatrices[host], hostSeg, h.at, TerminationReason::HETEROCLINIC,
                            sep.id, i);
                segmentsOfTriangle[f].push_back(SegRef{sep.id, i});
                ++report_.heteroclinicJoins;
                if (settings.resumeAfterTruncation) retractTail(host, hostSeg, h.at);
                return false;
            }

            // Sec. 3.3, third stopping condition. A streamline passing close to
            // a singularity crosses one of its separatrices orthogonally at its
            // closest approach -- that is what the hyperbolas of Prop. 1 do --
            // and this is that crossing whenever it happens near enough to the
            // singularity for the two curves to be the same curve as far as the
            // discretisation can tell.
            if (settings.cutAtSingularities && osep.originKind == SeparatrixOrigin::Singularity &&
                osep.origin_singularity_id >= 0 &&
                normP(h.at - singularities[osep.origin_singularity_id].coordinates) <=
                    cutRadius[osep.origin_singularity_id]) {
                terminateAt(sep, i, h.at, TerminationReason::CUT_AT_SINGULARITY, h.other.sep,
                            h.other.seg);
                segmentsOfTriangle[f].push_back(SegRef{sep.id, i});
                return false;
            }

            // The same argument as the heteroclinic join, for two separatrices
            // running the same way: in the continuum they are one streamline,
            // so this is not a crossing of two families but the discretisation
            // letting one drift across the other. The arriving one stops on
            // the other as a T-junction (docs/viertel_2019.md Sec. 4.4, the
            // row the paper leaves open). A separatrix drifting across its own
            // earlier path this way is winding onto a limit cycle, and stops
            // on the first such meeting rather than the second.
            if (settings.stopParallelTangential && phi < settings.tangentialAngle) {
                terminateAt(sep, i, h.at, TerminationReason::TANGENTIAL, h.other.sep,
                            h.other.seg);
                segmentsOfTriangle[f].push_back(SegRef{sep.id, i});
                ++report_.tangentialParallel;
                return false;
            }

            crossings.push_back(SeparatrixCrossing{sep.id, i, h.other.sep, h.other.seg, h.at});

            const int before = crossCount[sep.id][h.other.sep]++;
            crossCount[h.other.sep][sep.id]++;
            if (before >= 1) {
                // Sec. 3.3, second stopping condition. Cutting only the
                // separatrix whose step produced the crossing keeps the rule
                // local: the other one keeps its own count and stops on its own
                // terms.
                crossings.pop_back();
                terminateAt(sep, i, h.at, TerminationReason::CROSSED_TWICE, h.other.sep,
                            h.other.seg);
                segmentsOfTriangle[f].push_back(SegRef{sep.id, i});
                return false;
            }
        }
        segmentsOfTriangle[f].push_back(SegRef{sep.id, i});
    }
    return true;
}

// ---------------------------------------------------------------------------
// snapTangentialLanding()  --  Sec. 3.3's tangential rule, at the boundary
//
// The paper applies it to two separatrices meeting head on: in the continuum
// they could only cross squarely, so a shallow crossing is one curve the
// discretisation split in two, and it joins them. A separatrix meeting the
// *boundary* is the same statement -- the field is boundary aligned, so a
// streamline either arrives square or runs alongside -- and the paper does not
// apply it there, it assumes a mesh fine enough that the case cannot arise.
//
// Where it does arise the streamline is running alongside the boundary towards
// a corner. At a convex corner the field turns through ninety degrees, so a
// streamline that came in parallel to one edge leaves along the other: the
// corner is where it is going, and where it lands otherwise is an accident of
// which triangle the polyline happened to cross out of.
// ---------------------------------------------------------------------------
bool SeparatrixTrace::snapTangentialLanding(Separatrix &sep) {
    if (sep.path.size() < 2) return false;
    const Point &end = sep.path.back().global_pos;
    const Point in = end - sep.path[sep.path.size() - 2].global_pos;
    const double nin = normP(in);
    if (nin < 1e-18) return false;

    // The boundary it landed on: the edge it left through, or the two edges
    // meeting at the vertex it left at.
    std::vector<int> edges;
    if (sep.endBoundaryEdge >= 0) {
        edges.push_back(sep.endBoundaryEdge);
    } else if (sep.endBoundaryVertex >= 0) {
        for (const int e : mesh->boundaryEdges)
            if (mesh->edges[e][0] == sep.endBoundaryVertex || mesh->edges[e][1] == sep.endBoundaryVertex)
                edges.push_back(e);
    }
    if (edges.empty()) return false;

    // How square the arrival is, taken against whichever boundary edge it is
    // most nearly parallel to: running along one of the two edges meeting at a
    // vertex is running along the boundary.
    double bestSin = 1.0;
    for (const int e : edges) {
        const Point d = mesh->vertices[mesh->edges[e][1]] - mesh->vertices[mesh->edges[e][0]];
        const double nd = normP(d);
        if (nd < 1e-18) continue;
        bestSin = std::min(bestSin, std::fabs(cross2(in, d)) / (nin * nd));
    }
    if (bestSin < std::sin(settings.tangentialAngle)) ++report_.tangentialBoundaryExits;
    if (bestSin > std::sin(settings.boundaryTangentialAngle)) return false;   // arrived squarely

    // The corner it was heading for.
    const double reach = settings.boundaryCornerSnap * tracer->averageEdgeLength();
    int bestVertex = -1;
    double bestDist = reach;
    for (const BoundaryCorner &c : boundaryCorners) {
        const double d = normP(mesh->vertices[c.vertex] - end);
        if (d < bestDist) { bestDist = d; bestVertex = c.vertex; }
    }
    if (bestVertex < 0) return false;

    sep.path.back().global_pos = mesh->vertices[bestVertex];
    sep.endBoundaryVertex = bestVertex;
    sep.endBoundaryEdge = -1;
    return true;
}

// ---------------------------------------------------------------------------
// joinCornerLanding()  --  the corner-to-corner heteroclinic connection
//
// Where a neck meets the body of a model, the reflex corners at its two sides
// each send a separatrix straight across it, and in the continuum the two are
// one streamline running corner to corner. Traced, each lands a fraction of an
// element beside the other corner, and the layout gets two curves across the
// neck a ten-thousandth of an edge apart: a strip between two boundaries that
// Sec. 4 may not collapse (a rung from a corner to the boundary), which cuts
// every chord running along the neck in two, and leaves a node where the two
// crossings of each longitudinal separatrix merge. On geom022 that is what
// keeps the neck, and with it a fifth of the model, from ever becoming blocks.
//
// It is the heteroclinic join of Sec. 4.4 with the corner standing in for the
// other separatrix's origin: square arrival (a tangential one is
// snapTangentialLanding()'s), within Settings::cornerJoinRadius of the corner,
// heading within tangentialAngle of straight back along one of the corner's
// own separatrices. The arrival is moved onto the corner and that separatrix
// is cut back to its origin through the same retractTail() a join uses, so
// crossings on it are withdrawn and anything that had stopped on it resumes.
// ---------------------------------------------------------------------------
bool SeparatrixTrace::joinCornerLanding(Separatrix &sep) {
    if (!(settings.cornerJoinRadius > 0.0) || sep.path.size() < 2) return false;
    const double h = tracer->averageEdgeLength();
    const Point end = sep.path.back().global_pos;

    // Directions measured over about an edge, not over the last or first
    // segment alone: a segment can be a sliver of a triangle.
    auto pointBack = [&](const std::vector<TracePoint> &p, bool fromEnd) {
        const Point o = fromEnd ? p.back().global_pos : p.front().global_pos;
        const int n = static_cast<int>(p.size());
        for (int i = 1; i < n; ++i) {
            const Point &q = p[fromEnd ? n - 1 - i : i].global_pos;
            if (normP(q - o) >= h) return q;
        }
        return fromEnd ? p.front().global_pos : p.back().global_pos;
    };
    Point in = end - pointBack(sep.path, true);
    if (normP(in) < 1e-12 * h) return false;
    in = in / normP(in);

    int bestSep = -1, bestCorner = -1;
    double bestDist = settings.cornerJoinRadius * h;
    for (size_t ci = 0; ci < boundaryCorners.size(); ++ci) {
        const Point c = mesh->vertices[boundaryCorners[ci].vertex];
        const double d = normP(c - end);
        if (d >= bestDist) continue;
        if (sep.originKind == SeparatrixOrigin::BoundaryCorner &&
            sep.origin_singularity_id == static_cast<int>(ci))
            continue;   // back to where it started: not a connection
        for (const Separatrix &o : separatrices) {
            if (o.id == sep.id || o.originKind != SeparatrixOrigin::BoundaryCorner ||
                o.origin_singularity_id != static_cast<int>(ci) || o.path.size() < 2)
                continue;
            Point out = pointBack(o.path, false) - o.path.front().global_pos;
            if (normP(out) < 1e-12 * h) continue;
            out = out / normP(out);
            if (dotP(in, out) > -std::cos(settings.tangentialAngle)) continue;
            bestSep = o.id;
            bestCorner = static_cast<int>(ci);
            bestDist = d;
        }
    }
    if (bestSep < 0) return false;
    const Point c = mesh->vertices[boundaryCorners[bestCorner].vertex];

    // One that crossed another of the corner's separatrices on its way in, near
    // the corner, did not arrive end on: it came in across the corner's own
    // wedge (on geom035 it grazes the corner a hundred-thousandth of an edge
    // out), and ending it on the corner would leave that crossing sitting on
    // top of the corner with the two curves still crossing beside it. Those
    // are left as they landed.
    for (const SeparatrixCrossing &x : crossings) {
        const int other = (x.sepA == sep.id) ? x.sepB : (x.sepB == sep.id) ? x.sepA : -1;
        if (other < 0 || separatrices[other].originKind != SeparatrixOrigin::BoundaryCorner ||
            separatrices[other].origin_singularity_id != bestCorner)
            continue;
        if (normP(x.pos - c) <= settings.cornerJoinRadius * h) return false;
    }

    sep.path.back().global_pos = c;
    sep.endBoundaryVertex = boundaryCorners[bestCorner].vertex;
    sep.endBoundaryEdge = -1;

    Separatrix &o = separatrices[bestSep];
    retractTail(o.id, 1, c);
    interfaceCrossings.erase(std::remove_if(interfaceCrossings.begin(), interfaceCrossings.end(),
                                            [&](const InterfaceCrossing &x) { return x.sep == o.id; }),
                             interfaceCrossings.end());
    terminateAt(o, 0, c, TerminationReason::HETEROCLINIC, sep.id,
                static_cast<int>(sep.path.size()) - 1);
    return true;
}

// ---------------------------------------------------------------------------
// advanceOne()  --  one triangle of one separatrix, and what it ran into
// ---------------------------------------------------------------------------
bool SeparatrixTrace::advanceOne(int k) {
    Separatrix &sep = separatrices[k];
    if (!sep.active) return false;

    const int firstNew = static_cast<int>(sep.path.size());
    FieldTracer::CutInfo cut;
    const FieldTracer::Status st = tracer->advance(walkers[k], sep.path, &cut);
    for (size_t i = static_cast<size_t>(std::max(firstNew, 1)); i < sep.path.size(); ++i)
        arcLength[k] += normP(sep.path[i].global_pos - sep.path[i - 1].global_pos);

    if (static_cast<int>(sep.path.size()) > firstNew) {
        if (!registerNewSegments(sep, firstNew)) return false;
    }
    if (interfaces && st == FieldTracer::Status::Ok) recordInterfaceCrossing(k);

    switch (st) {
        case FieldTracer::Status::Boundary:
            sep.active = false;
            sep.termination_reason = TerminationReason::EXIT_BOUNDARY;
            sep.endBoundaryVertex = walkers[k].atVertex;
            sep.endBoundaryEdge =
                (walkers[k].atVertex < 0 && walkers[k].entryEdge >= 0)
                    ? mesh->triangleEdges[walkers[k].tri][walkers[k].entryEdge]
                    : -1;
            if (snapTangentialLanding(sep)) ++tangentialLandings;
            else if (joinCornerLanding(sep)) ++report_.cornerJoins;
            break;

        case FieldTracer::Status::Cut: {
            if (!settings.cutAtSingularities) {
                // Asked not to apply condition 3: nothing sensible is left
                // to do inside the singular triangle, so stop without
                // claiming a T-junction.
                sep.active = false;
                sep.termination_reason = TerminationReason::STUCK;
                break;
            }
            const int target = (cut.singularity >= 0)
                                   ? singularities[cut.singularity].portSeparatrixIds[cut.port]
                                   : -1;
            sep.active = false;
            sep.termination_reason = TerminationReason::CUT_AT_SINGULARITY;
            sep.endOnSeparatrix = target;
            sep.endOnSegment = 1;  // the port's first segment, out of the singular triangle
            break;
        }

        case FieldTracer::Status::Stuck:
            sep.active = false;
            sep.termination_reason = TerminationReason::STUCK;
            break;

        case FieldTracer::Status::Ok:
            if (++stepsTaken[k] >= settings.maxStepsPerSeparatrix) {
                sep.active = false;
                sep.termination_reason = TerminationReason::LIMIT_CYCLE;
            }
            break;
    }
    return sep.active;
}

// ---------------------------------------------------------------------------
// stepAndCheck()  --  one round
//
// With Settings::growByArcLength the round is "every live separatrix up to the
// next front": the front moves on by one mean edge, and until every live one
// has reached it the one furthest behind is advanced by a triangle. That is the
// spec's priority queue on arc length, in slices a frame of the viewer can
// show. Without it a round is the older rule of one triangle each, in id order.
// ---------------------------------------------------------------------------
void SeparatrixTrace::stepAndCheck() {
    ++steps;

    if (!settings.growByArcLength) {
        for (size_t k = 0; k < separatrices.size(); ++k) advanceOne(static_cast<int>(k));
        resumeQueue.clear();   // resumed ones are live, and next round takes them
    } else {
        front += tracer->averageEdgeLength();
        using Item = std::pair<double, int>;
        std::priority_queue<Item, std::vector<Item>, std::greater<Item>> queue;
        for (size_t k = 0; k < separatrices.size(); ++k)
            if (separatrices[k].active && arcLength[k] < front)
                queue.push({arcLength[k], static_cast<int>(k)});
        while (!queue.empty()) {
            const int k = queue.top().second;
            queue.pop();
            if (advanceOne(k) && arcLength[k] < front) queue.push({arcLength[k], k});
            for (const int r : resumeQueue)
                if (separatrices[r].active && arcLength[r] < front) queue.push({arcLength[r], r});
            resumeQueue.clear();
        }
    }

    finishedTracing = countActive() == 0;

    // Ports keep their identity through the trace, so the public copy of the
    // singularity table is refreshed rather than rebuilt.
    singularities = tracer->getSingularities();
}

void SeparatrixTrace::run() {
    while (!finishedTracing) stepAndCheck();
}
