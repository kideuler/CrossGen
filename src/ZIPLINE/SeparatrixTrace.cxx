#include "ZIPLINE/SeparatrixTrace.hxx"

#include "MERIDIAN/Interfaces.hxx"
#include "ZIPLINE/ConformalField.hxx"

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
                                 const Interfaces *itf,
                                 std::shared_ptr<const ConformalField> conformal)
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
    if (conformal && !interfaces)
        tracer = std::make_unique<FieldTracer>(cf, std::move(conformal));
    else
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

    if (settings.captureAtSingularities) {
        buildTargets();
        double lo[2] = {1e300, 1e300}, hi[2] = {-1e300, -1e300};
        for (const Point &p : mesh->vertices)
            for (int k = 0; k < 2; ++k) { lo[k] = std::min(lo[k], p[k]); hi[k] = std::max(hi[k], p[k]); }
        maxLength = settings.maxLengthDiagonals * std::hypot(hi[0] - lo[0], hi[1] - lo[1]);
    }

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

    // Nothing to trace is finished -- unless there is a hole nothing will
    // touch, which the first round seeds (seedUntouchedHoles()).
    finishedTracing = separatrices.empty() && !(settings.seedUntouchedHoles && boundaryLoops.size() > 1);
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
    // On the corrected field, the corners it was solved with, which are
    // Table 1's unless the first field turned one the other way
    // (ConformalField's constructor has the case).
    const ConformalField *cf = tracer->getConformal();
    for (const auto &loop : boundaryLoops) {
        for (const int v : loop) {
            const int q = (cf && v < static_cast<int>(cf->cornerQuarters().size()))
                              ? cf->cornerQuarters()[v]
                              : quartersOf(interior[v]);
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
bool SeparatrixTrace::registerNewSegments(Separatrix &sep, int firstNewIndex, bool stopRules) {
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
        // The triangles a meeting with this segment can be filed under. Its
        // own, and -- when it starts on the boundary of its triangle -- every
        // other triangle holding its start point: the one it was handed over
        // from across an edge, or the whole fan of the vertex it left. Two
        // separatrices meeting exactly at a hand-off point count the meeting
        // on the two segments that *start* there (the half-open convention
        // above), and those can lie in different triangles: one crossing an
        // edge from T to T', the other from T' to T at the same point, or two
        // leaving one vertex into different wedges. A field the mesh resolves
        // exactly makes that common -- the corrected field of a staircase
        // with no singularity is the constant cross, its separatrices are the
        // mesh's own grid lines, and they cross at its vertices (geom027) --
        // and a crossing missed there is two arcs touching with no node.
        int extra[64];
        int nExtra = 0;
        {
            const int pf = sep.path[i - 1].face_id;
            if (pf >= 0 && pf != f && pf < static_cast<int>(segmentsOfTriangle.size())) extra[nExtra++] = pf;
            for (int k = 0; k < 3; ++k) {
                const int v = mesh->triangles[f][k];
                if (a0[0] != mesh->vertices[v][0] || a0[1] != mesh->vertices[v][1]) continue;
                const auto range = mesh->vertexTriangles.trianglesForVertex(v);
                for (const int *it = range.first; it != range.second && nExtra < 64; ++it)
                    if (*it != f && *it != pf) extra[nExtra++] = *it;
            }
        }
        std::vector<SegRef> candidates = segmentsOfTriangle[f];
        for (int x = 0; x < nExtra; ++x)
            candidates.insert(candidates.end(), segmentsOfTriangle[extra[x]].begin(),
                              segmentsOfTriangle[extra[x]].end());
        for (const SegRef &other : candidates) {
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

        // docs/quad_layout_pipeline.md Sec. 7.1: separatrices cross one another
        // freely, and a crossing is a node of the layout and nothing more.
        // Always so for the last segment of a capture, inside the critical
        // patch, where Viertel's rules would stop a separatrix half way
        // into the point it is being joined to.
        if (!settings.stopAtCrossings || !stopRules) {
            for (const Hit &h : hits) {
                crossings.push_back(SeparatrixCrossing{sep.id, i, h.other.sep, h.other.seg, h.at});
                crossCount[sep.id][h.other.sep]++;
                if (h.other.sep != sep.id) crossCount[h.other.sep][sep.id]++;
            }
            segmentsOfTriangle[f].push_back(SegRef{sep.id, i});
            continue;
        }

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
    FieldTracer::Status st = tracer->advance(walkers[k], sep.path, &cut);
    for (size_t i = static_cast<size_t>(std::max(firstNew, 1)); i < sep.path.size(); ++i)
        arcLength[k] += normP(sep.path[i].global_pos - sep.path[i - 1].global_pos);

    // docs/quad_layout_pipeline.md Sec. 7.5 rule 2: into a critical patch, and
    // joined to its point. Looked for before the crossings of this step are
    // registered, since none past the patch's edge is on the separatrix any
    // more. Viertel's condition 3 inside a singular triangle (Status::Cut) is
    // the same event seen by the sweep, and becomes the same capture.
    if (settings.captureAtSingularities && static_cast<int>(sep.path.size()) > firstNew) {
        int seg = -1, g = -1;
        double t = 0.0;
        if (findCapture(k, firstNew, seg, t, g) && capture(k, seg, t, g, firstNew)) return false;
        if (st == FieldTracer::Status::Cut) {
            const int gc = targetOf(SeparatrixOrigin::Singularity, cut.singularity);
            const bool excluded = gc < 0 ||
                std::find(excludedTargets[k].begin(), excludedTargets[k].end(), gc) != excludedTargets[k].end() ||
                (gc == originTarget[k] && !leftOrigin[k]);
            if (!excluded &&
                capture(k, static_cast<int>(sep.path.size()) - 1, 1.0, gc, firstNew))
                return false;
            // Not taken: with Viertel's rules on as well it is his condition 3
            // and the T-junction stands; without them it carries on past.
            if (!settings.stopAtCrossings) st = FieldTracer::Status::Ok;
        }
        if (!leftOrigin[k] && originTarget[k] >= 0) {
            const Target &o = targets[originTarget[k]];
            for (size_t i = static_cast<size_t>(firstNew); i < sep.path.size(); ++i)
                if (normP(sep.path[i].global_pos - o.pos) > o.radius) leftOrigin[k] = 1;
        }
    }

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
            // The corner join first. A separatrix that runs exactly along the
            // continuation of an edge of a reflex corner -- which on the
            // corrected field of a model with only right-angled corners and
            // no singularity is every corner separatrix, since that field is
            // the constant cross -- reaches the opposite corner's vertex
            // parallel to that corner's own edge, and so reads as a tangential
            // landing. It is the corner-to-corner separatrix traced from its
            // other end (docs/quad_layout_pipeline.md Sec. 8.1, the duplicate
            // in its purest form), and snapping it instead of joining it left
            // the two copies lying exactly on top of each other: an arc
            // crossing no segment test can see, and no faces at all
            // (geom018). On the first field the order made no difference: an
            // arrival that exactly parallel never happened there.
            // With capture on, a landing beside a corner that emits has
            // already been taken by that corner's critical patch, and a join
            // here would decide the copy question before Sec. 8 can.
            if (!settings.captureAtSingularities && joinCornerLanding(sep)) ++report_.cornerJoins;
            else if (snapTangentialLanding(sep)) ++tangentialLandings;
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
            if (++stepsTaken[k] >= settings.maxStepsPerSeparatrix ||
                (settings.captureAtSingularities && maxLength > 0.0 && arcLength[k] > maxLength)) {
                sep.active = false;
                sep.termination_reason = TerminationReason::LIMIT_CYCLE;
            }
            break;
    }
    return sep.active;
}

// ---------------------------------------------------------------------------
// seedUntouchedHoles()  --  docs/quad_layout_pipeline.md Sec. 6.9's annulus
//
// A hole that no separatrix leaves from or lands on, once tracing has stopped,
// sits inside one face of the layout, and that face is not a disc: it has the
// hole for a second boundary, and no block can be made of it. The commonest
// case is the whole model -- an annulus with no corner, where Poincare-Hopf
// asks for no singularity and none is traced (Sec. 6.9: "with no singularities
// there are no separatrices, so no layout"; Viertel's spec lists it as a known
// limitation, Sec. 11 test 10). The doc's remedy is to impose a 3-5 pair. But
// an annulus needs no singularity to be quad-meshed: the polar grid has none,
// and its radial lines are streamlines of the field's other family. So the
// hole is cut along four of them instead, launched from it square to it at
// four vertices a quarter of its length apart. Each launch vertex becomes a
// corner of its loop with quarters = 2 -- flat, adding nothing to
// Poincare-Hopf -- so the layout puts a node there, and each face meeting it
// gets one corner on either side of the streamline. They are traced like any
// separatrix, under the same rules, and land on whatever bounds the face. Four
// rather than two: two would leave two faces each a half annulus, sound but
// bent through a half turn.
// ---------------------------------------------------------------------------
int SeparatrixTrace::seedUntouchedHoles() {
    if (boundaryLoops.empty()) return 0;
    // Which loop each boundary vertex is on.
    std::vector<int> loopOf(mesh->vertices.size(), -1);
    for (size_t l = 0; l < boundaryLoops.size(); ++l)
        for (const int v : boundaryLoops[l]) loopOf[v] = static_cast<int>(l);
    std::vector<char> touched(boundaryLoops.size(), 0);
    for (const BoundaryCorner &c : boundaryCorners)
        if (c.vertex >= 0 && loopOf[c.vertex] >= 0) touched[loopOf[c.vertex]] = 1;
    for (const InterfaceNodeInfo &n : interfaceNodes)
        if (n.vertex >= 0 && loopOf[n.vertex] >= 0) touched[loopOf[n.vertex]] = 1;
    for (const Separatrix &sep : separatrices) {
        if (sep.termination_reason != TerminationReason::EXIT_BOUNDARY) continue;
        int v = sep.endBoundaryVertex;
        if (v < 0 && sep.endBoundaryEdge >= 0) v = mesh->edges[sep.endBoundaryEdge][0];
        if (v >= 0 && loopOf[v] >= 0) touched[loopOf[v]] = 1;
    }
    int launched = 0;
    for (size_t l = 0; l < boundaryLoops.size(); ++l) {
        const auto &loop = boundaryLoops[l];
        const int n = static_cast<int>(loop.size());
        if (touched[l] || n < 8) continue;
        // Holes only: the loop turns by -2 pi with the domain on its left.
        double turning = 0.0, length = 0.0;
        for (int i = 0; i < n; ++i) {
            const Point &a = mesh->vertices[loop[(i + n - 1) % n]], &b = mesh->vertices[loop[i]];
            const Point &c = mesh->vertices[loop[(i + 1) % n]];
            turning += std::atan2(cross2(b - a, c - b), dotP(b - a, c - b));
            length += normP(c - b);
        }
        if (turning > 0.0) continue;
        // Start at the vertex furthest along +x, so the seeds do not depend on
        // where the loop happens to begin.
        int start = 0;
        for (int i = 1; i < n; ++i)
            if (mesh->vertices[loop[i]][0] > mesh->vertices[loop[start]][0]) start = i;
        double along = 0.0;
        int next = 0;
        for (int k = 0; k < n && next < 4; ++k) {
            const int i = (start + k) % n;
            if (along >= next * 0.25 * length - 1e-12 * length) {
                const int v = loop[i];
                const Point &a = mesh->vertices[loop[(i + n - 1) % n]], &c = mesh->vertices[loop[(i + 1) % n]];
                // Square to the loop, into the domain: the left normal of its
                // tangent, the domain being on its left.
                const Point t = normalizeP(c - a);
                const double dir = computeAngle(Point{-t[1], t[0]});
                BoundaryCorner bc;
                bc.vertex = v;
                bc.quarters = 2;
                const int ci = static_cast<int>(boundaryCorners.size());
                boundaryCorners.push_back(bc);
                const int id = launchFromVertex(v, dir, SeparatrixOrigin::BoundaryCorner, ci, 0);
                if (id >= 0) {
                    boundaryCorners.back().separatrixIds.push_back(id);
                    arcLength.push_back(0.0);
                    if (settings.captureAtSingularities) {
                        originTarget.push_back(-1);
                        leftOrigin.push_back(1);
                        excludedTargets.emplace_back();
                        captureResume.push_back(Resume{});
                    }
                    ++launched;
                }
                ++next;
            }
            along += normP(mesh->vertices[loop[(i + 1) % n]] - mesh->vertices[loop[i]]);
        }
    }
    report_.holeSeeds += launched;
    return launched;
}

// ===========================================================================
// Settings::captureAtSingularities: docs/quad_layout_pipeline.md Secs. 7-9
// ===========================================================================

// ---------------------------------------------------------------------------
// buildTargets()  --  the singular points, their critical patches, their ports
//
// Every interior singularity, and every corner of dS and node of the interface
// network that emits: the places a separatrix can end other than dS, which is
// where Sec. 3.3 of the doc puts the singularities of the pattern (corners
// included, Sec. 2.4). A corner that emits nothing -- a convex one -- is not a
// place a streamline of a boundary-aligned field can arrive anywhere but along
// dS, so it is no target.
//
// The patch radius is captureRadius in mean edges *of the point's own
// neighbourhood* -- the singular triangle's edges, or the boundary edges at a
// corner -- since a near miss is a matter of the local discretisation.
// ---------------------------------------------------------------------------
void SeparatrixTrace::buildTargets() {
    targets.clear();
    const double h = tracer->averageEdgeLength();
    auto addTarget = [&](Target g) {
        if (g.portSeps.empty()) return;
        targets.push_back(std::move(g));
    };
    for (size_t i = 0; i < singularities.size(); ++i) {
        const Singularity &s = singularities[i];
        Target g;
        g.pos = s.coordinates;
        g.kind = SeparatrixOrigin::Singularity;
        g.index = static_cast<int>(i);
        g.triangle = s.triangleIndex;
        const Triangle &tri = mesh->triangles[s.triangleIndex];
        double e = 0.0;
        for (int k = 0; k < 3; ++k) e += normP(mesh->vertices[tri[(k + 1) % 3]] - mesh->vertices[tri[k]]);
        g.radius = settings.captureRadius * e / 3.0;
        for (int p = 0; p < s.numPorts(); ++p) {
            if (p >= static_cast<int>(s.portSeparatrixIds.size()) || s.portSeparatrixIds[p] < 0) continue;
            g.portDirs.push_back(s.portAngles[p]);
            g.portSeps.push_back(s.portSeparatrixIds[p]);
        }
        addTarget(std::move(g));
    }
    auto vertexTarget = [&](int v, SeparatrixOrigin kind, int index, const std::vector<int> &seps) {
        Target g;
        g.pos = mesh->vertices[v];
        g.kind = kind;
        g.index = index;
        g.vertex = v;
        double e = 0.0;
        int n = 0;
        for (const int be : mesh->boundaryEdges)
            if (mesh->edges[be][0] == v || mesh->edges[be][1] == v) {
                e += normP(mesh->vertices[mesh->edges[be][1]] - mesh->vertices[mesh->edges[be][0]]);
                ++n;
            }
        g.radius = settings.captureRadius * (n > 0 ? e / n : h);
        for (const int id : seps) {
            if (id < 0 || separatrices[id].path.empty()) continue;
            g.portDirs.push_back(separatrices[id].path[0].theta);
            g.portSeps.push_back(id);
        }
        addTarget(std::move(g));
    };
    for (size_t ci = 0; ci < boundaryCorners.size(); ++ci)
        vertexTarget(boundaryCorners[ci].vertex, SeparatrixOrigin::BoundaryCorner, static_cast<int>(ci),
                     boundaryCorners[ci].separatrixIds);
    for (size_t ni = 0; ni < interfaceNodes.size(); ++ni)
        vertexTarget(interfaceNodes[ni].vertex, SeparatrixOrigin::InterfaceNode, static_cast<int>(ni),
                     interfaceNodes[ni].separatrixIds);

    // How far off a port an arrival may be and still be running back along
    // it: captureMisfit of the half-gap between that port and the nearest
    // direction beside it -- the next port, or at a corner or node a wall of
    // its wedge. A quarter of the sector at a singularity (30 degrees at a
    // valence 3, 18 at a valence 5), 22.5 degrees at a right-angled reflex
    // corner.
    for (Target &g : targets) {
        std::vector<double> sides = g.portDirs;
        if (g.vertex >= 0)
            for (const int be : mesh->boundaryEdges) {
                const int a = mesh->edges[be][0], b = mesh->edges[be][1];
                if (a != g.vertex && b != g.vertex) continue;
                sides.push_back(computeAngle(mesh->vertices[a == g.vertex ? b : a] - g.pos));
            }
        double gap = 2.0 * M_PI;
        for (size_t p = 0; p < g.portDirs.size(); ++p)
            for (size_t q = 0; q < sides.size(); ++q) {
                if (q == p) continue;
                const double d = std::fabs(wrap_pi(sides[q] - g.portDirs[p]));
                if (d > 1e-9) gap = std::min(gap, d);
            }
        g.maxMisfit = settings.captureMisfit * 0.5 * gap;
        g.sideDirs = sides;
    }

    // Which triangles each patch reaches: a walk out from the point's own
    // triangle(s) over every triangle within the radius of it.
    targetsOfTriangle.assign(mesh->triangles.size(), {});
    auto distToTriangle = [&](int f, const Point &p) {
        const std::array<double, 3> b = tracer->barycentric(f, p);
        if (b[0] >= 0.0 && b[1] >= 0.0 && b[2] >= 0.0) return 0.0;
        double best = std::numeric_limits<double>::max();
        const Triangle &tri = mesh->triangles[f];
        for (int k = 0; k < 3; ++k) {
            const Point &a = mesh->vertices[tri[k]], &c = mesh->vertices[tri[(k + 1) % 3]];
            const Point d = c - a;
            const double dd = dotP(d, d);
            const double t = dd > 0.0 ? std::min(1.0, std::max(0.0, dotP(p - a, d) / dd)) : 0.0;
            best = std::min(best, normP(p - (a + d * t)));
        }
        return best;
    };
    for (size_t gi = 0; gi < targets.size(); ++gi) {
        const Target &g = targets[gi];
        std::vector<int> seeds;
        if (g.triangle >= 0) seeds.push_back(g.triangle);
        if (g.vertex >= 0) {
            const auto range = mesh->vertexTriangles.trianglesForVertex(g.vertex);
            for (const int *it = range.first; it != range.second; ++it) seeds.push_back(*it);
        }
        std::vector<int> stack = seeds;
        std::unordered_map<int, char> seen;
        for (const int f : seeds) seen[f] = 1;
        while (!stack.empty()) {
            const int f = stack.back();
            stack.pop_back();
            if (distToTriangle(f, g.pos) > g.radius) continue;
            targetsOfTriangle[f].push_back(static_cast<int>(gi));
            for (int k = 0; k < 3; ++k) {
                const int n = mesh->triangleAdjacency[f][k];
                if (n >= 0 && !seen[n]) { seen[n] = 1; stack.push_back(n); }
            }
        }
    }

    const size_t n = separatrices.size();
    originTarget.assign(n, -1);
    leftOrigin.assign(n, 0);
    excludedTargets.assign(n, {});
    captureResume.assign(n, Resume{});
    for (size_t gi = 0; gi < targets.size(); ++gi)
        for (const int id : targets[gi].portSeps)
            if (id >= 0 && id < static_cast<int>(n)) originTarget[id] = static_cast<int>(gi);
}

int SeparatrixTrace::targetOf(SeparatrixOrigin kind, int index) const {
    for (size_t gi = 0; gi < targets.size(); ++gi)
        if (targets[gi].kind == kind && targets[gi].index == index) return static_cast<int>(gi);
    return -1;
}

int SeparatrixTrace::nearestPort(const Target &g, double dir, double &misfit) const {
    int best = -1;
    misfit = std::numeric_limits<double>::max();
    for (size_t p = 0; p < g.portDirs.size(); ++p) {
        const double d = std::fabs(wrap_pi(dir - g.portDirs[p]));
        if (d < misfit) { misfit = d; best = static_cast<int>(p); }
    }
    return best;
}

// ---------------------------------------------------------------------------
// findCapture()
//
// "An arriving separatrix stops on entering it" (Sec. 7.2): the first point of
// the new segments inside a patch the separatrix may end in -- not its own
// origin's until it has been out of it once, not one a resolution already
// released it from.
// ---------------------------------------------------------------------------
bool SeparatrixTrace::findCapture(int k, int firstNew, int &seg, double &t, int &target) const {
    const Separatrix &sep = separatrices[k];
    for (int i = std::max(firstNew, 1); i < static_cast<int>(sep.path.size()); ++i) {
        const int f = sep.path[i].face_id;
        if (f < 0 || f >= static_cast<int>(targetsOfTriangle.size())) continue;
        const Point &a = sep.path[i - 1].global_pos;
        const Point &b = sep.path[i].global_pos;
        double bestT = 2.0;
        int bestG = -1;
        for (const int g : targetsOfTriangle[f]) {
            if (g == originTarget[k] && !leftOrigin[k]) continue;
            if (std::find(excludedTargets[k].begin(), excludedTargets[k].end(), g) != excludedTargets[k].end())
                continue;
            const Target &tg = targets[g];
            // Where |a + s (b - a) - X| first reaches the radius.
            const Point d = b - a, w = a - tg.pos;
            const double A = dotP(d, d), B = 2.0 * dotP(w, d), C = dotP(w, w) - tg.radius * tg.radius;
            double s = -1.0;
            if (C <= 0.0) s = 0.0;
            else if (A > 0.0) {
                const double disc = B * B - 4.0 * A * C;
                if (disc >= 0.0) {
                    const double r = (-B - std::sqrt(disc)) / (2.0 * A);
                    if (r >= 0.0 && r <= 1.0) s = r;
                }
            }
            if (s >= 0.0 && s < bestT) { bestT = s; bestG = g; }
        }
        if (bestG >= 0) {
            seg = i;
            t = bestT;
            target = bestG;
            return true;
        }
    }
    return false;
}

// ---------------------------------------------------------------------------
// approach()  --  P1 Fig. 13's last segment, straight to the point
//
// Walked through the triangles it crosses so the path keeps one triangle per
// segment, which everything that intersects paths relies on.
// ---------------------------------------------------------------------------
bool SeparatrixTrace::approach(int tri, const Point &from, const Target &g,
                               std::vector<TracePoint> &out, const Point *via) const {
    out.clear();
    int f = tri;
    Point p = from;
    int exclude = -1;
    if (via) {
        // First to the waypoint, as a target of its own that lies in
        // whichever triangle the walk finds it in.
        Target w = g;
        w.pos = *via;
        w.triangle = -2;
        w.vertex = -1;
        const Point d = *via - from;
        const double dir0 = computeAngle(d);
        for (int guard = 0; guard < 256; ++guard) {
            const double L = normP(*via - p);
            const std::array<double, 3> b = tracer->barycentric(f, *via);
            if (L <= 1e-12 * tracer->averageEdgeLength() ||
                (b[0] >= -1e-12 && b[1] >= -1e-12 && b[2] >= -1e-12)) {
                out.push_back(TracePoint{*via, f, dir0});
                p = *via;
                break;
            }
            Point hit{0.0, 0.0};
            int edge = -1;
            double along = 0.0;
            if (!tracer->exitOfRay(f, p, (*via - p) / L, exclude, hit, edge, along)) return false;
            out.push_back(TracePoint{hit, f, dir0});
            const int nb = mesh->triangleAdjacency[f][edge];
            if (nb < 0) return false;
            const int va = mesh->triangles[f][edge], vb = mesh->triangles[f][(edge + 1) % 3];
            exclude = -1;
            for (int i = 0; i < 3; ++i) {
                const int a = mesh->triangles[nb][i], b2 = mesh->triangles[nb][(i + 1) % 3];
                if ((a == va && b2 == vb) || (a == vb && b2 == va)) exclude = i;
            }
            f = nb;
            p = hit;
            if (guard == 255) return false;
        }
        exclude = -1;
    }
    const double dir = computeAngle(g.pos - p);
    auto holds = [&](int t) {
        if (g.triangle >= 0) return t == g.triangle;
        const Triangle &tr = mesh->triangles[t];
        return tr[0] == g.vertex || tr[1] == g.vertex || tr[2] == g.vertex;
    };
    for (int guard = 0; guard < 256 && f >= 0; ++guard) {
        if (holds(f)) {
            out.push_back(TracePoint{g.pos, f, dir});
            return true;
        }
        const Point d = g.pos - p;
        const double L = normP(d);
        if (L <= 1e-12 * tracer->averageEdgeLength()) {
            out.push_back(TracePoint{g.pos, f, dir});
            return true;
        }
        Point hit{0.0, 0.0};
        int edge = -1;
        double along = 0.0;
        if (!tracer->exitOfRay(f, p, d / L, exclude, hit, edge, along)) return false;
        if (normP(hit - p) >= L) {
            out.push_back(TracePoint{g.pos, f, dir});
            return true;
        }
        out.push_back(TracePoint{hit, f, dir});
        const int nb = mesh->triangleAdjacency[f][edge];
        if (nb < 0) return false;
        const int va = mesh->triangles[f][edge], vb = mesh->triangles[f][(edge + 1) % 3];
        exclude = -1;
        for (int i = 0; i < 3; ++i) {
            const int a = mesh->triangles[nb][i], b = mesh->triangles[nb][(i + 1) % 3];
            if ((a == va && b == vb) || (a == vb && b == va)) exclude = i;
        }
        f = nb;
        p = hit;
    }
    return false;
}

// ---------------------------------------------------------------------------
// capture()  --  end separatrix k on target g
//
// The separatrix is cut where it entered the patch, the crossings up to there
// registered as any step's are, and the straight approach appended and
// registered too. Where it was cut, and which way the field ran there, is kept,
// so that resolveCaptures() can undo the capture and trace on from that point.
// ---------------------------------------------------------------------------
bool SeparatrixTrace::capture(int k, int seg, double t, int target, int firstNew) {
    Separatrix &sep = separatrices[k];
    const Target &g = targets[target];
    const Point q = sep.path[seg - 1].global_pos * (1.0 - t) + sep.path[seg].global_pos * t;
    const int tri = sep.path[seg].face_id;
    if (tri < 0) return false;

    // The arrival slot (Sec. 7.6): the port whose separatrix this one is the
    // other copy of, which is the port it was travelling back along -- read
    // off the last stretch of the path, about a patch radius of it, before it
    // entered. Not the direction of the entry point from the point, which is
    // what Sec. 7.4's pseudocode stores: for a separatrix running past a
    // point a fraction of the radius off -- what the concentric arcs round a
    // hole do, geom009 -- the entry point lies abeam of it, the radial
    // direction is a port or more away from the one it was running along,
    // and the copy removed is a separatrix that was no copy at all while the
    // real one runs on alongside.
    const double arrival = [&] {
        for (int i = seg - 1; i >= 0; --i)
            if (normP(sep.path[i].global_pos - q) >= g.radius) return computeAngle(sep.path[i].global_pos - q);
        return computeAngle(sep.path.front().global_pos - q);
    }();
    double misfit = 0.0;
    const int port = nearestPort(g, arrival, misfit);
    // Running back along no port at all is passing the point, not reaching
    // it (Sec. 7.2's false capture): a corner's separatrix sliding into the
    // patch of the next corner along dS, a streamline swinging round a
    // singularity on its hyperbola.
    if (port < 0 || misfit > g.maxMisfit) return false;
    // And it must enter between the rays of the directions beside that port
    // -- the neighbouring ports, or a wall of a corner's wedge -- and not
    // beyond one, for its last segment to reach the point without crossing
    // the point's own separatrices. The last segment is then not P1 Fig.
    // 13's straight one but two: from the entry to half way down the port's
    // ray, then along the ray. Both lie in the one convex sector between
    // the port and the neighbour on the entry's side, so they meet nothing,
    // and the separatrix comes in along the slot it claims. Straight from an
    // entry abeam of the next port, it ran down that port's own separatrix
    // instead (geom033).
    const double radialDir = computeAngle(q - g.pos);
    {
        const double off = wrap_pi(radialDir - g.portDirs[port]);
        // The nearest side direction on the entry's side of the port.
        double side = 2.0 * M_PI;
        for (const double a : g.sideDirs) {
            const double d = wrap_pi(a - g.portDirs[port]);
            if (std::fabs(d) < 1e-9) continue;
            if ((d > 0.0) == (off > 0.0)) side = std::min(side, std::fabs(d));
        }
        if (std::fabs(off) >= side) return false;
    }
    const Point mid = g.pos + Point{std::cos(g.portDirs[port]), std::sin(g.portDirs[port])} *
                                  (0.5 * normP(q - g.pos));

    std::vector<TracePoint> tail;
    if (!approach(tri, q, g, tail, &mid)) return false;
    // Nor may the last segment cross a separatrix already there to reach the
    // point: arriving *behind* another curve that passes the point closer
    // than this one entered, it is that curve's capture to make, if
    // anyone's, and joining this one through it makes a crossing at
    // whatever angle the two happen to meet -- ten degrees on geom032, and a
    // two-cornered sliver either side. Its own point's separatrices are no
    // obstacle: the legs lie in one sector of them and meet them only at it.
    {
        Point prev = q;
        for (const TracePoint &tp : tail) {
            const Point &a = prev, &b = tp.global_pos;
            prev = b;
            const int f = tp.face_id;
            if (f < 0 || f >= static_cast<int>(segmentsOfTriangle.size())) continue;
            for (const SegRef &other : segmentsOfTriangle[f]) {
                if (other.sep == k) continue;
                const Separatrix &o = separatrices[other.sep];
                if (other.seg <= 0 || other.seg >= static_cast<int>(o.path.size())) continue;
                if (o.originKind == g.kind && o.origin_singularity_id == g.index) continue;
                double tt = 0.0, uu = 0.0;
                Point at{0.0, 0.0};
                if (properIntersection(a, b, o.path[other.seg - 1].global_pos, o.path[other.seg].global_pos,
                                       tt, uu, at) &&
                    normP(at - g.pos) > 1e-9 * tracer->averageEdgeLength())
                    return false;
            }
        }
    }

    const double along = computeAngle(sep.path[seg].global_pos - sep.path[seg - 1].global_pos);
    sep.path.resize(seg + 1);
    sep.path.back().global_pos = q;
    captureResume[k] = Resume{tri, q, tracer->angleAt(tri, q, along), seg};
    if (!registerNewSegments(sep, std::max(firstNew, 1))) return true;   // a Viertel rule stopped it first

    const int before = static_cast<int>(sep.path.size());
    sep.path.insert(sep.path.end(), tail.begin(), tail.end());
    registerNewSegments(sep, before, false);

    sep.active = false;
    sep.termination_reason = TerminationReason::AT_SINGULARITY;
    sep.endKind = g.kind;
    sep.endSingularity = g.index;
    sep.endPort = port;
    sep.arrivalMisfit = misfit;
    double L = 0.0;
    for (size_t i = 1; i < sep.path.size(); ++i) L += normP(sep.path[i].global_pos - sep.path[i - 1].global_pos);
    arcLength[k] = L;
    ++report_.captures;
    return true;
}

// ---------------------------------------------------------------------------
// dropSegments()  --  separatrix k from segment `fromSeg` on, taken back
// ---------------------------------------------------------------------------
void SeparatrixTrace::dropSegments(int k, int fromSeg) {
    Separatrix &sep = separatrices[k];
    fromSeg = std::max(fromSeg, 1);
    for (int i = fromSeg; i < static_cast<int>(sep.path.size()); ++i) {
        const int f = sep.path[i].face_id;
        if (f < 0 || f >= static_cast<int>(segmentsOfTriangle.size())) continue;
        auto &bucket = segmentsOfTriangle[f];
        bucket.erase(std::remove_if(bucket.begin(), bucket.end(),
                                    [&](const SegRef &r) { return r.sep == k && r.seg == i; }),
                     bucket.end());
    }
    std::vector<SeparatrixCrossing> kept;
    kept.reserve(crossings.size());
    for (const SeparatrixCrossing &x : crossings) {
        const bool gone = (x.sepA == k && x.segA >= fromSeg) || (x.sepB == k && x.segB >= fromSeg);
        if (!gone) { kept.push_back(x); continue; }
        --crossCount[x.sepA][x.sepB];
        if (x.sepA != x.sepB) --crossCount[x.sepB][x.sepA];
    }
    crossings.swap(kept);
    interfaceCrossings.erase(std::remove_if(interfaceCrossings.begin(), interfaceCrossings.end(),
                                            [&](const InterfaceCrossing &x) {
                                                return x.sep == k && x.pathIndex >= fromSeg;
                                            }),
                             interfaceCrossings.end());
    if (static_cast<int>(sep.path.size()) > fromSeg) sep.path.resize(fromSeg);
    double L = 0.0;
    for (size_t i = 1; i < sep.path.size(); ++i) L += normP(sep.path[i].global_pos - sep.path[i - 1].global_pos);
    arcLength[k] = L;
}

// ---------------------------------------------------------------------------
// resolveCaptures()  --  Sec. 8, P1's Algorithm 1 with Sec. 8.7's slot test
//
// A singular point's port is a *slot* (Sec. 7.6): one separatrix leaves by it,
// and in the finished layout exactly one curve may use it. A capture claims the
// slot it arrived in. When separatrix s, leaving slot i of A, is captured in
// slot j of B, then s *is* B's separatrix j traced from the other end, and
// that separatrix -- call it u -- is s's copy:
//
//   * if u was captured in slot i of A, the two are the true pair of Sec. 8.1
//     (start slot of one = arrival slot of the other, Sec. 8.7c), and one is
//     removed -- the one that arrived further off its slot (Sec. 8.7d),
//     where P1 takes one at random;
//   * if u ran on somewhere else, it missed A by more than A's patch while s
//     came within B's -- the same near miss, measured at the two ends in two
//     different local scales (exp H differs at A and B; Sec. 6.2) -- and u is
//     removed all the same: the slot belongs to the connection, and what u
//     did after passing A is the discretisation's.
//
// Done greedily, the smallest misfit first, as a matching on the slots: a
// capture whose slot is already taken, or whose copy u is itself a connection
// kept to somewhere else, cannot stand, and is undone -- the separatrix traced
// on past the point, which may not catch it again. That is P1 Sec. 10's
// re-propagation (Sec. 12.4 of the doc). Connections kept in an earlier pass go
// first in the next, so a later capture cannot unseat one.
//
// Returns whether anything was restarted; tracing then carries on, and this
// runs again when it stops.
// ---------------------------------------------------------------------------
bool SeparatrixTrace::resolveCaptures() {
    const int n = static_cast<int>(separatrices.size());
    std::vector<int> caps;
    for (int s = 0; s < n; ++s)
        if (separatrices[s].termination_reason == TerminationReason::AT_SINGULARITY) caps.push_back(s);
    std::vector<char> &kept = connectionKept_;
    kept.resize(n, 0);
    std::sort(caps.begin(), caps.end(), [&](int a, int b) {
        if (kept[a] != kept[b]) return kept[a] > kept[b];
        return separatrices[a].arrivalMisfit < separatrices[b].arrivalMisfit;
    });
    auto slotOf = [&](int s) -> std::pair<int, int> {   // origin slot (target, port)
        const int g = originTarget[s];
        if (g < 0) return {-1, -1};
        const auto &ps = targets[g].portSeps;
        for (size_t p = 0; p < ps.size(); ++p) if (ps[p] == s) return {g, static_cast<int>(p)};
        return {g, -1};
    };
    std::vector<char> dead(n, 0);
    std::vector<char> keep(n, 0);
    std::unordered_map<long long, int> claimed;
    auto key = [](int g, int p) { return static_cast<long long>(g) * 1024 + p; };
    std::vector<int> release;
    // A connection kept by an earlier pass already removed its copy, which
    // cannot come back: it holds its two slots outright, and only what has
    // been captured since is weighed against it.
    for (const int s : caps) {
        if (!kept[s]) continue;
        const Separatrix &S = separatrices[s];
        const int g = targetOf(S.endKind, S.endSingularity);
        keep[s] = 1;
        claimed[key(g, S.endPort)] = s;
        const auto origin = slotOf(s);
        if (origin.first >= 0) claimed[key(origin.first, origin.second)] = s;
    }
    for (const int s : caps) {
        if (dead[s] || kept[s]) continue;
        const Separatrix &S = separatrices[s];
        const int g = targetOf(S.endKind, S.endSingularity);
        const int p = S.endPort;
        if (g < 0 || p < 0 || p >= static_cast<int>(targets[g].portSeps.size())) { release.push_back(s); continue; }
        const int u = targets[g].portSeps[p];
        const auto origin = slotOf(s);
        if (claimed.count(key(g, p)) || u == s || (u >= 0 && keep[u]) ||
            (origin.first >= 0 && claimed.count(key(origin.first, origin.second)))) {
            release.push_back(s);
            continue;
        }
        // P1's own redundancy test (Sec. 8.2: same ends, same separatrices
        // crossed), in its geometric form. Two connections between the same
        // pair of points that leave and arrive by *different* slots pass the
        // slot test, and are still one connection traced twice when the loop
        // they make together encloses nothing -- no hole, no other singular
        // point: the region between them is then P1 Fig. 16's strip of
        // slivers ending in two triangles, a two-cornered face in the layout.
        // Enclosing a hole, they are the torus case of Sec. 8.3 and both are
        // real. One connection that leaves a point and comes back to it is
        // tested against itself the same way. The worse of a lens is undone.
        bool lens = false;
        for (int o = 0; o < n && !lens; ++o) {
            if (!keep[o] && o != s) continue;
            const Separatrix &O = separatrices[o];
            const int oa = originTarget[o], ob = targetOf(O.endKind, O.endSingularity);
            const int sa = originTarget[s];
            if (o == s ? (sa != g) : !((oa == sa && ob == g) || (oa == g && ob == sa))) continue;
            std::vector<Point> loop;
            for (const TracePoint &tp : S.path) loop.push_back(tp.global_pos);
            if (o != s) {
                if (oa == sa)
                    for (auto it = O.path.rbegin(); it != O.path.rend(); ++it) loop.push_back(it->global_pos);
                else
                    for (const TracePoint &tp : O.path) loop.push_back(tp.global_pos);
            }
            lens = !enclosesSomething(loop, sa, g);
        }
        if (lens) {
            release.push_back(s);
            ++report_.lenses;
            continue;
        }
        keep[s] = 1;
        claimed[key(g, p)] = s;
        if (origin.first >= 0) claimed[key(origin.first, origin.second)] = s;
        if (u >= 0 && !dead[u]) {
            dead[u] = 1;
            const Separatrix &U = separatrices[u];
            if (U.termination_reason == TerminationReason::AT_SINGULARITY && origin.first >= 0 &&
                targetOf(U.endKind, U.endSingularity) == origin.first && U.endPort == origin.second)
                ++report_.reciprocal;
        }
    }

    bool restarted = false;
    for (int u = 0; u < n; ++u) {
        if (!dead[u] || separatrices[u].termination_reason == TerminationReason::REDUNDANT) continue;
        Separatrix &U = separatrices[u];
        dropSegments(u, 1);
        U.path.resize(1);
        U.active = false;
        U.termination_reason = TerminationReason::REDUNDANT;
        U.endOnSeparatrix = U.endOnSegment = -1;
        ++report_.redundant;
        if (resumeStoppedOn(u)) restarted = true;
    }
    for (int s = 0; s < n; ++s)
        if (keep[s] && !kept[s]) { kept[s] = 1; ++report_.connections; }

    for (const int s : release) {
        // Released earlier in this pass and then taken as the copy of a
        // connection kept after it: it is gone, not restarted.
        if (dead[s]) continue;
        Separatrix &S = separatrices[s];
        const Resume &r = captureResume[s];
        const int g = targetOf(S.endKind, S.endSingularity);
        if (g >= 0) excludedTargets[s].push_back(g);
        if (r.tri < 0 || r.seg < 1) continue;
        dropSegments(s, r.seg + 1);
        S.path.back().global_pos = r.pos;
        S.active = true;
        S.termination_reason = TerminationReason::RUNNING;
        S.endSingularity = S.endPort = -1;
        S.arrivalMisfit = 0.0;
        walkers[s] = tracer->startAt(r.tri, r.pos, r.dir);
        kept[s] = 0;
        restarted = true;
        ++report_.released;
    }
    return restarted;
}

// ---------------------------------------------------------------------------
// enclosesSomething()  --  whether a closed curve has topology inside it
//
// A hole is inside when one of its boundary vertices is (no separatrix crosses
// dS, so a hole is wholly in or wholly out); a singular point when its
// position is, the two the curve joins excepted. By winding number, so a curve
// that crosses itself is read the way a face walk would read it.
// ---------------------------------------------------------------------------
bool SeparatrixTrace::enclosesSomething(const std::vector<Point> &loop, int skipA, int skipB) const {
    if (loop.size() < 3) return false;
    auto winding = [&](const Point &q) {
        double w = 0.0;
        for (size_t i = 0; i < loop.size(); ++i) {
            const Point a = loop[i] - q, b = loop[(i + 1) % loop.size()] - q;
            w += std::atan2(cross2(a, b), dotP(a, b));
        }
        return std::lround(w / (2.0 * M_PI));
    };
    for (const auto &bl : boundaryLoops) {
        if (bl.empty()) continue;
        // Holes only: the outer loop cannot be inside anything.
        double turning = 0.0;
        for (size_t i = 0; i < bl.size(); ++i) {
            const Point &a = mesh->vertices[bl[(i + bl.size() - 1) % bl.size()]];
            const Point &b = mesh->vertices[bl[i]];
            const Point &c = mesh->vertices[bl[(i + 1) % bl.size()]];
            turning += std::atan2(cross2(b - a, c - b), dotP(b - a, c - b));
        }
        if (turning > 0.0) continue;
        if (winding(mesh->vertices[bl[0]]) != 0) return true;
    }
    for (size_t g = 0; g < targets.size(); ++g) {
        if (static_cast<int>(g) == skipA || static_cast<int>(g) == skipB) continue;
        if (targets[g].kind != SeparatrixOrigin::Singularity) continue;
        if (winding(targets[g].pos) != 0) return true;
    }
    return false;
}

// ---------------------------------------------------------------------------
// resumeStoppedOn()  --  what stopped against a removed copy, restarted
//
// Only with Viertel's stops on as well (Settings::stopAtCrossings): a
// separatrix can then have ended on the copy just removed -- a second crossing
// of it, a near miss of a singularity on it, a tangential meeting, a
// head-on join -- and would be a loose end. It is restarted from where it
// stopped, as retractTail() restarts the ones a join cuts the ground from
// under, with the stopping crossing's count taken back.
// ---------------------------------------------------------------------------
bool SeparatrixTrace::resumeStoppedOn(int host) {
    bool any = false;
    for (Separatrix &z : separatrices) {
        if (z.active || z.id == host || z.endOnSeparatrix != host || z.path.size() < 2) continue;
        const TerminationReason r = z.termination_reason;
        if (r != TerminationReason::CROSSED_TWICE && r != TerminationReason::TANGENTIAL &&
            r != TerminationReason::CUT_AT_SINGULARITY && r != TerminationReason::HETEROCLINIC &&
            r != TerminationReason::CUT_LIMIT_CYCLE)
            continue;
        if (r == TerminationReason::CROSSED_TWICE && crossCount[z.id][host] > 0) {
            --crossCount[z.id][host];
            --crossCount[host][z.id];
        }
        z.active = true;
        z.termination_reason = TerminationReason::RUNNING;
        z.endOnSeparatrix = z.endOnSegment = -1;
        const TracePoint &tp = z.path.back();
        walkers[z.id] = tracer->startAt(tp.face_id, tp.global_pos, tp.theta);
        ++report_.resumed;
        any = true;
    }
    return any;
}

// ---------------------------------------------------------------------------
// cutLimitCycles()  --  Sec. 9, after the reduction
//
// Sec. 9.2's reading: a pair crossing each other at least twice is a possible
// limit cycle, and of the two the authentic one is the one whose crossings with
// the other are spread furthest along it; so, outright, is a separatrix that
// ran out of length or reached dS far from square. Each authentic one is cut at
// its first crossing, counted from its origin and outside its origin's patch,
// that is within orthogonalTolerance of square (Sec. 9.3's "closest orthogonal
// intersection", read as the doc reads it): it ends there, a T-junction on the
// separatrix it was crossing.
// ---------------------------------------------------------------------------
void SeparatrixTrace::cutLimitCycles() {
    const int n = static_cast<int>(separatrices.size());
    // Arc length along each path, per point.
    std::vector<std::vector<double>> s(n);
    for (int k = 0; k < n; ++k) {
        const auto &p = separatrices[k].path;
        s[k].assign(p.size(), 0.0);
        for (size_t i = 1; i < p.size(); ++i) s[k][i] = s[k][i - 1] + normP(p[i].global_pos - p[i - 1].global_pos);
    }
    auto param = [&](int k, int seg, const Point &at) {
        const auto &p = separatrices[k].path;
        if (seg < 1 || seg >= static_cast<int>(p.size())) return 0.0;
        return s[k][seg - 1] + normP(at - p[seg - 1].global_pos);
    };
    struct Span { double lo = 1e300, hi = -1e300; int count = 0; };
    std::unordered_map<long long, Span> span;   // (a, b): a's parameters at crossings with b
    for (const SeparatrixCrossing &x : crossings) {
        if (x.sepA == x.sepB) continue;
        for (int side = 0; side < 2; ++side) {
            const int a = side ? x.sepB : x.sepA, b = side ? x.sepA : x.sepB;
            const double t = side ? param(x.sepB, x.segB, x.pos) : param(x.sepA, x.segA, x.pos);
            Span &sp = span[static_cast<long long>(a) * n + b];
            sp.lo = std::min(sp.lo, t);
            sp.hi = std::max(sp.hi, t);
            ++sp.count;
        }
    }
    std::vector<char> authentic(n, 0);
    for (const auto &kv : span) {
        const int a = static_cast<int>(kv.first / n), b = static_cast<int>(kv.first % n);
        if (a > b || kv.second.count < 2) continue;
        const auto it = span.find(static_cast<long long>(b) * n + a);
        if (it == span.end()) continue;
        // A connection, joined to a singular point at both ends, is no
        // limit cycle whatever it crosses: its pair is the one that winds.
        // Cutting it would also free its arrival slot while the copy it
        // replaced stays removed (Sec. 8), and empty that slot.
        const bool ca = separatrices[a].termination_reason == TerminationReason::AT_SINGULARITY;
        const bool cb = separatrices[b].termination_reason == TerminationReason::AT_SINGULARITY;
        if (ca && cb) continue;
        const double sa = kv.second.hi - kv.second.lo, sb = it->second.hi - it->second.lo;
        authentic[ca ? b : cb ? a : (sa >= sb ? a : b)] = 1;
    }
    const double tol = settings.orthogonalTolerance;
    for (int k = 0; k < n; ++k) {
        const Separatrix &sep = separatrices[k];
        if (sep.termination_reason == TerminationReason::LIMIT_CYCLE ||
            sep.termination_reason == TerminationReason::STUCK)
            authentic[k] = 1;
        if (sep.termination_reason == TerminationReason::EXIT_BOUNDARY && sep.path.size() >= 2) {
            // How square it arrived, against the boundary edge(s) it left by.
            const Point in = sep.path.back().global_pos - sep.path[sep.path.size() - 2].global_pos;
            double worst = 0.0;
            bool any = false;
            auto off = [&](int e) {
                const Point d = mesh->vertices[mesh->edges[e][1]] - mesh->vertices[mesh->edges[e][0]];
                const double c = std::fabs(dotP(in, d)) / std::max(1e-300, normP(in) * normP(d));
                return M_PI_2 - std::acos(std::min(1.0, c));   // 0 when square
            };
            if (sep.endBoundaryEdge >= 0) { worst = off(sep.endBoundaryEdge); any = true; }
            else if (sep.endBoundaryVertex >= 0) {
                double best = 1e300;
                for (const int e : mesh->boundaryEdges)
                    if (mesh->edges[e][0] == sep.endBoundaryVertex || mesh->edges[e][1] == sep.endBoundaryVertex)
                        best = std::min(best, off(e));
                if (best < 1e300) { worst = best; any = true; }
            }
            if (any && worst > tol) authentic[k] = 1;
        }
    }

    for (int k = 0; k < n; ++k) {
        if (!authentic[k]) continue;
        ++report_.limitCycles;
        Separatrix &sep = separatrices[k];
        const Target *o = (originTarget[k] >= 0) ? &targets[originTarget[k]] : nullptr;
        // The first square crossing outside the origin's patch.
        int best = -1;
        double bestT = 1e300;
        for (size_t xi = 0; xi < crossings.size(); ++xi) {
            const SeparatrixCrossing &x = crossings[xi];
            if (x.sepA != k && x.sepB != k) continue;
            if (x.sepA == x.sepB) continue;
            const bool mine = (x.sepA == k);
            const int seg = mine ? x.segA : x.segB;
            const int other = mine ? x.sepB : x.sepA, oseg = mine ? x.segB : x.segA;
            // On a separatrix that stays whole: one that is cut itself may
            // lose the stretch this one would end on.
            if (authentic[other]) continue;
            if (o && normP(x.pos - o->pos) <= o->radius) continue;
            const auto &op = separatrices[other].path;
            if (seg < 1 || oseg < 1 || oseg >= static_cast<int>(op.size())) continue;
            const Point a = sep.path[seg].global_pos - sep.path[seg - 1].global_pos;
            const Point b = op[oseg].global_pos - op[oseg - 1].global_pos;
            const double c = std::fabs(dotP(a, b)) / std::max(1e-300, normP(a) * normP(b));
            if (std::acos(std::min(1.0, c)) < M_PI_2 - tol) continue;   // not square
            const double t = param(k, seg, x.pos);
            if (t < bestT) { bestT = t; best = static_cast<int>(xi); }
        }
        if (best < 0) continue;
        const SeparatrixCrossing x = crossings[best];
        const bool mine = (x.sepA == k);
        const int seg = mine ? x.segA : x.segB;
        const int other = mine ? x.sepB : x.sepA, oseg = mine ? x.segB : x.segA;
        // Everything after the cut goes, including the crossing itself, which
        // is now where the separatrix ends on the other.
        dropSegments(k, seg + 1);
        std::vector<SeparatrixCrossing> keepX;
        for (const SeparatrixCrossing &y : crossings) {
            const bool onCut = (y.sepA == k && y.segA == seg) || (y.sepB == k && y.segB == seg);
            if (onCut) {
                const double ty = (y.sepA == k) ? param(k, y.segA, y.pos) : param(k, y.segB, y.pos);
                if (ty >= bestT - 1e-12 * tracer->averageEdgeLength()) {
                    --crossCount[y.sepA][y.sepB];
                    --crossCount[y.sepB][y.sepA];
                    continue;
                }
            }
            keepX.push_back(y);
        }
        crossings.swap(keepX);
        sep.path.back().global_pos = x.pos;
        {
            double L = 0.0;
            for (size_t i = 1; i < sep.path.size(); ++i)
                L += normP(sep.path[i].global_pos - sep.path[i - 1].global_pos);
            arcLength[k] = L;
        }
        sep.active = false;
        sep.termination_reason = TerminationReason::CUT_LIMIT_CYCLE;
        sep.endOnSeparatrix = other;
        sep.endOnSegment = oseg;
        ++report_.limitCyclesCut;
    }
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

    // A hole nothing has touched gets streamlines of its own, once, before
    // anything else is settled (seedUntouchedHoles()).
    if (finishedTracing && settings.seedUntouchedHoles && !holesSeeded) {
        holesSeeded = true;
        if (seedUntouchedHoles() > 0) finishedTracing = false;
    }

    // With capture on, stopping is not the end: Secs. 8 and 9 run once
    // nothing is moving, and the reduction may restart some of what stopped.
    // Bounded, should captures keep being undone and remade: after the last
    // pass a restarted separatrix may be caught by nothing but dS.
    if (finishedTracing && settings.captureAtSingularities && !capturesResolved) {
        bool restarted = false;
        if (settings.reduceSeparatrices) {
            ++resolvePasses;
            restarted = resolveCaptures();
            if (restarted && resolvePasses >= 16)
                for (size_t k = 0; k < separatrices.size(); ++k)
                    if (separatrices[k].active)
                        for (size_t g = 0; g < targets.size(); ++g) excludedTargets[k].push_back(static_cast<int>(g));
        }
        if (restarted) {
            finishedTracing = false;
        } else {
            // With Viertel's stops on, a separatrix that winds has already
            // stopped at its second crossing (his condition 2, which is P2's
            // "possible limit cycle" acted on as it happens): there is nothing
            // left for Sec. 9 to cut.
            if (settings.cutLimitCycles && !settings.stopAtCrossings) cutLimitCycles();
            // Sec. 8.6 step 4: every slot used exactly once.
            for (size_t g = 0; g < targets.size(); ++g) {
                bool bad = false;
                for (size_t p = 0; p < targets[g].portSeps.size(); ++p) {
                    int uses = 0;
                    const int own = targets[g].portSeps[p];
                    if (own >= 0 && separatrices[own].termination_reason != TerminationReason::REDUNDANT) ++uses;
                    for (const Separatrix &z : separatrices)
                        if (z.termination_reason == TerminationReason::AT_SINGULARITY &&
                            targetOf(z.endKind, z.endSingularity) == static_cast<int>(g) &&
                            z.endPort == static_cast<int>(p))
                            ++uses;
                    if (uses != 1) bad = true;
                }
                if (bad) ++report_.slotAuditFailures;
            }
            capturesResolved = true;
        }
    }

    // Ports keep their identity through the trace, so the public copy of the
    // singularity table is refreshed rather than rebuilt.
    singularities = tracer->getSingularities();
}

void SeparatrixTrace::run() {
    while (!finishedTracing) stepAndCheck();
}
