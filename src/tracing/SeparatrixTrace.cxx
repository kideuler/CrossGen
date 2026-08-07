#include "tracing/SeparatrixTrace.hxx"

#include <algorithm>
#include <cmath>
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
    if (std::fabs(den) < 1e-30) return false;

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
                                 bool useActualSingularityCoordinates, const Settings &s)
    : crossField(cf), settings(s) {
    tracer = std::make_unique<FieldTracer>(cf, useActualSingularityCoordinates);
    mesh = &tracer->getMesh();

    segmentsOfTriangle.assign(mesh->triangles.size(), {});

    buildBoundary();
    launchFromSingularities();
    launchFromCorners();

    singularities = tracer->getSingularities();

    cutRadius.assign(singularities.size(), 0.0);
    for (size_t i = 0; i < singularities.size(); ++i) {
        const Triangle &tri = mesh->triangles[singularities[i].triangleIndex];
        double h = 0.0;
        for (int k = 0; k < 3; ++k)
            h += normP(mesh->vertices[tri[(k + 1) % 3]] - mesh->vertices[tri[k]]);
        cutRadius[i] = settings.singularityCutRadius * h / 3.0;
    }

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
    for (const auto &loop : boundaryLoops) {
        for (const int v : loop) {
            const int q = std::max(1, static_cast<int>(std::lround(interior[v] / M_PI_2)));
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

        auto it = nextOnLoop.find(c.vertex);
        if (it == nextOnLoop.end()) continue;
        const Point along = mesh->vertices[it->second] - mesh->vertices[c.vertex];
        if (normP(along) < 1e-18) continue;
        const double base = computeAngle(along);

        const auto range = mesh->vertexTriangles.trianglesForVertex(c.vertex);
        if (range.first == range.second) continue;

        const double step = settings.evenCornerRays ? (c.interiorAngle / c.quarters) : M_PI_2;
        for (int j = 1; j <= rays; ++j) {
            const double dir = wrap_pi(base + j * step);

            Separatrix sep;
            sep.id = static_cast<int>(separatrices.size());
            sep.originKind = SeparatrixOrigin::BoundaryCorner;
            sep.origin_singularity_id = ci;
            sep.origin_singularity_port = j - 1;

            Walker w;
            w.tri = *range.first;
            w.pos = mesh->vertices[c.vertex];
            w.entryEdge = -1;
            w.atVertex = c.vertex;
            w.dir = dir;
            w.crossDir = dir;

            TracePoint tp;
            tp.global_pos = w.pos;
            tp.face_id = w.tri;
            tp.theta = dir;
            sep.path.push_back(tp);

            c.separatrixIds.push_back(sep.id);
            walkers.push_back(w);
            stepsTaken.push_back(0);
            crossCount.emplace_back();
            separatrices.push_back(std::move(sep));
        }
    }
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
                terminateAt(sep, i, h.at, TerminationReason::HETEROCLINIC, h.other.sep,
                            h.other.seg);
                terminateAt(separatrices[h.other.sep], h.other.seg, h.at,
                            TerminationReason::HETEROCLINIC, sep.id, i);
                segmentsOfTriangle[f].push_back(SegRef{sep.id, i});
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

void SeparatrixTrace::stepAndCheck() {
    finishedTracing = true;
    ++steps;

    for (size_t k = 0; k < separatrices.size(); ++k) {
        Separatrix &sep = separatrices[k];
        if (!sep.active) continue;

        const int firstNew = static_cast<int>(sep.path.size());
        FieldTracer::CutInfo cut;
        const FieldTracer::Status st = tracer->advance(walkers[k], sep.path, &cut);

        if (static_cast<int>(sep.path.size()) > firstNew) {
            if (!registerNewSegments(sep, firstNew)) continue;
        }

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

        if (sep.active) finishedTracing = false;
    }

    // Ports keep their identity through the trace, so the public copy of the
    // singularity table is refreshed rather than rebuilt.
    singularities = tracer->getSingularities();
}

void SeparatrixTrace::run() {
    while (!finishedTracing) stepAndCheck();
}
