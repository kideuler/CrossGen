#include "MotorcycleGraph.hxx"

#include <algorithm>
#include <cmath>
#include <deque>
#include <fstream>
#include <iostream>
#include <limits>
#include <array>
#include <queue>
#include <stdexcept>

namespace {

// Solve J d = e for d, i.e. the model-space direction whose image is e.
inline Point solve2(const double J[2][2], const Point &e) {
    const double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    if (std::fabs(det) < 1e-18) return Point{0.0, 0.0};
    return Point{( J[1][1] * e[0] - J[0][1] * e[1]) / det,
                 (-J[1][0] * e[0] + J[0][0] * e[1]) / det};
}

} // namespace

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
MotorcycleGraph::MotorcycleGraph(const Polysquare &ps) {
    poly = &ps;
    mesh = &ps.getMesh();
    cut = &ps.getCutMesh();
    uv = &ps.getUV();
    if (uv->empty()) throw std::runtime_error("MotorcycleGraph: the polysquare is not solved");
}

// ---------------------------------------------------------------------------
// triangleJacobian()  --  grad phi on one triangle
// ---------------------------------------------------------------------------
void MotorcycleGraph::triangleJacobian(int f, double J[2][2]) const {
    const Triangle &t = mesh->triangles[f];
    const Point &p0 = mesh->vertices[t[0]];
    const Point &p1 = mesh->vertices[t[1]];
    const Point &p2 = mesh->vertices[t[2]];
    const double area2 = cross2(p1 - p0, p2 - p0);

    J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;
    if (std::fabs(area2) < 1e-18) return;

    const Point p[3] = {p0, p1, p2};
    for (int i = 0; i < 3; ++i) {
        const Point e = p[(i + 2) % 3] - p[(i + 1) % 3];
        const Point g{-e[1] / area2, e[0] / area2}; // grad lambda_i
        const Point &q = (*uv)[cut->triangles[f][i]];
        for (int r = 0; r < 2; ++r)
            for (int c = 0; c < 2; ++c) J[r][c] += q[r] * (c == 0 ? g[0] : g[1]);
    }
}

// ---------------------------------------------------------------------------
// exitPoint()  --  where a straight step leaves a triangle
// ---------------------------------------------------------------------------
bool MotorcycleGraph::exitPoint(int f, const Point &from, const Point &meshDir,
                                Point &hit, int &edge, double &along) const {
    const Triangle &t = mesh->triangles[f];
    double best = std::numeric_limits<double>::max();
    edge = -1;

    for (int i = 0; i < 3; ++i) {
        const int va = t[i], vb = t[(i + 1) % 3];
        const Point &A = mesh->vertices[va];
        const Point &B = mesh->vertices[vb];
        const Point e = B - A;

        // from + s*dir = A + u*e
        const double den = cross2(meshDir, e);
        if (std::fabs(den) < 1e-18) continue; // parallel to this edge
        const Point w = A - from;
        const double s = cross2(w, e) / den;
        const double u = cross2(w, meshDir) / den;
        if (s <= 1e-12) continue;                 // behind, or where we started
        if (u < -1e-9 || u > 1.0 + 1e-9) continue; // misses the segment

        if (s < best) {
            best = s;
            edge = i;
            along = std::min(1.0, std::max(0.0, u));
            hit = Point{from[0] + s * meshDir[0], from[1] + s * meshDir[1]};
        }
    }
    return edge >= 0;
}

// ---------------------------------------------------------------------------
// launch()  --  one ray per axis direction that points into the domain
//
// How many, and which way, is decided by the corner index rather than by
// testing directions for pointing into the domain. A corner of index k has an
// image wedge of 180 - k*90 degrees, so it holds 1 - k axis directions between
// its two boundary rays -- two at a reflex corner, three at a 360-degree one,
// none at a convex one -- and they sit at exact multiples of 90 degrees from
// the boundary.
//
// Reading them off the geometry instead does not survive contact with a real
// mesh. A corner whose image wedge comes out at 275 degrees rather than 270
// admits a third axis direction, and that third ray leaves within a few
// degrees of the boundary and runs alongside it for the length of the model,
// shaving the boundary layer off into slivers. On geom014 that turned 24 rays
// into 26 and 27 blocks into 41.
// ---------------------------------------------------------------------------
void MotorcycleGraph::launch() {
    bikes.clear();
    const auto &corners = poly->getBoundaryCorners();

    const int nV = static_cast<int>(mesh->vertices.size());
    std::vector<std::array<int, 2>> incidentBoundary(nV, std::array<int, 2>{-1, -1});
    std::vector<int> boundaryDegree(nV, 0);
    for (int e : mesh->boundaryEdges) {
        for (int k = 0; k < 2; ++k) {
            const int v = mesh->edges[e][k];
            if (boundaryDegree[v] < 2) incidentBoundary[v][boundaryDegree[v]] = e;
            ++boundaryDegree[v];
        }
    }

    auto imageAt = [&](int f, int origVertex) -> Point {
        const Triangle &t = mesh->triangles[f];
        for (int i = 0; i < 3; ++i) if (t[i] == origVertex) return (*uv)[cut->triangles[f][i]];
        return Point{0.0, 0.0};
    };

    for (int v : mesh->boundaryVertices) {
        if (v >= static_cast<int>(corners.size())) continue;
        const int k = corners[v];
        if (k >= 0 || boundaryDegree[v] != 2) continue;

        // The star, walked from one boundary edge to the other.
        int entry = incidentBoundary[v][0];
        int cur = (mesh->edgeTriangles[entry][0] >= 0) ? mesh->edgeTriangles[entry][0]
                                                       : mesh->edgeTriangles[entry][1];
        std::vector<int> fan;
        const int guard = mesh->vertexTriangles.vertexDegree(v) + 1;
        int firstEdge = entry;
        while (cur >= 0) {
            fan.push_back(cur);
            int next = -1;
            for (int i = 0; i < 3; ++i) {
                const int e = mesh->triangleEdges[cur][i];
                if (e < 0 || e == entry) continue;
                if (mesh->edges[e][0] == v || mesh->edges[e][1] == v) { next = e; break; }
            }
            if (next < 0) { fan.clear(); break; }
            if (mesh->isBoundaryEdge[next]) break;
            const int nb = (mesh->edgeTriangles[next][0] == cur) ? mesh->edgeTriangles[next][1]
                                                                 : mesh->edgeTriangles[next][0];
            if (nb < 0 || static_cast<int>(fan.size()) > guard) { fan.clear(); break; }
            entry = next;
            cur = nb;
        }
        if (fan.empty()) continue;

        // Walk it counter-clockwise, the sense the corner index is measured in.
        {
            const Triangle &tri = mesh->triangles[fan.front()];
            int i0 = -1;
            for (int i = 0; i < 3; ++i) if (tri[i] == v) i0 = i;
            if (i0 < 0) continue;
            const int va = tri[(i0 + 1) % 3], vb = tri[(i0 + 2) % 3];
            const double sgn = cross2(mesh->vertices[va] - mesh->vertices[v],
                                      mesh->vertices[vb] - mesh->vertices[v]);
            const int w = (mesh->edges[firstEdge][0] == v) ? mesh->edges[firstEdge][1]
                                                           : mesh->edges[firstEdge][0];
            const bool ccw = (w == va) ? (sgn > 0.0) : (sgn < 0.0);
            if (!ccw) { std::reverse(fan.begin(), fan.end()); firstEdge = incidentBoundary[v][1]; }
        }

        // The boundary ray the wedge starts from, snapped to an axis: after
        // snapBoundary() it is on one, to within a rounding error.
        const int w0 = (mesh->edges[firstEdge][0] == v) ? mesh->edges[firstEdge][1]
                                                        : mesh->edges[firstEdge][0];
        const Point startImage = imageAt(fan.front(), w0) - imageAt(fan.front(), v);
        if (normP(startImage) < 1e-15) continue;
        const int startAxis = static_cast<int>(
            std::lround(computeAngle(startImage) / M_PI_2) + 8) % 4;

        // How many rays and which way is settled; all that is left is which
        // triangle of the star each one starts in. That is a containment
        // question, so it is answered by testing -- the triangle whose wedge
        // the direction sits furthest inside -- rather than by adding up image
        // angles, which puts a ray in the neighbouring triangle whenever the
        // wedges do not sum quite as expected and leaves it pointing out of
        // the triangle it was given, with nowhere to go on its first step.
        const Point axes[4] = {{1.0, 0.0}, {0.0, 1.0}, {-1.0, 0.0}, {0.0, -1.0}};
        for (int j = 1; j <= 1 - k; ++j) {
            const Point e = axes[(startAxis + j) % 4];

            int host = -1;
            double bestMargin = -std::numeric_limits<double>::max();
            Point bestDir{0.0, 0.0};

            for (int f : fan) {
                const Triangle &tri = mesh->triangles[f];
                int i0 = -1;
                for (int i = 0; i < 3; ++i) if (tri[i] == v) i0 = i;
                if (i0 < 0) continue;

                double J[2][2];
                triangleJacobian(f, J);
                const Point d = solve2(J, e);
                const double dn = normP(d);
                if (dn < 1e-18) continue;

                const Point a = mesh->vertices[tri[(i0 + 1) % 3]] - mesh->vertices[v];
                const Point b = mesh->vertices[tri[(i0 + 2) % 3]] - mesh->vertices[v];
                const double na = normP(a), nb = normP(b);
                if (na < 1e-18 || nb < 1e-18) continue;

                // Sines of the angles from each wedge edge to the direction,
                // signed so that both are positive exactly when it is inside.
                const double s = (cross2(a, b) > 0.0) ? 1.0 : -1.0;
                const double margin = std::min(s * cross2(a, d) / (na * dn),
                                               s * cross2(d, b) / (nb * dn));
                if (margin > bestMargin) { bestMargin = margin; host = f; bestDir = d; }
            }
            // A direction that points into no triangle of the star is not a
            // direction into the domain: the corner index says the wedge holds
            // it and the parameterization there says otherwise. That happens
            // where Sec. 4.3 left the image inconsistent -- a folded triangle,
            // or a corner whose wedge did not come out at the angle its index
            // claims -- and launching the ray anyway only produces one that
            // cannot take a step. Skipped and counted instead.
            if (host < 0 || bestMargin < -1e-9) { ++report_.skippedCorners; continue; }

            Motorcycle m;
            m.tri = host;
            m.originVertex = v;
            m.pos = mesh->vertices[v];
            m.dir = e;
            m.id = static_cast<int>(bikes.size());
            bikes.push_back(m);

            // Seal the corner the ray leaves from. A ray's wall is the edges
            // it crosses, and the first of those is the far side of the
            // triangle it starts in -- so without this the flood walks from
            // one side of the ray to the other straight through that triangle,
            // around the vertex, and the wall separates nothing. Both edges of
            // the launch triangle at the corner are walled instead, which
            // closes the gap on either side; the triangle itself is then
            // isolated and folded back in as an unresolved block.
            for (int i = 0; i < 3; ++i) {
                const int ed = mesh->triangleEdges[host][i];
                if (ed < 0) continue;
                if (mesh->edges[ed][0] == v || mesh->edges[ed][1] == v) {
                    tracedEdges.insert(EdgeKey(mesh->edges[ed][0], mesh->edges[ed][1]));
                }
            }
        }
    }
    report_.motorcycles = static_cast<int>(bikes.size());
}

// ---------------------------------------------------------------------------
// run()  --  advance whichever ray has travelled least
// ---------------------------------------------------------------------------
void MotorcycleGraph::run() {
    traces.assign(bikes.size(), {});
    segments.clear();
    exitEdge.assign(bikes.size(), -1);
    exitAlong.assign(bikes.size(), 0.0);
    onWall.assign(mesh->triangles.size(), 0);
    for (size_t i = 0; i < bikes.size(); ++i) {
        traces[i].push_back(bikes[i].pos);
        onWall[bikes[i].tri] = 1;
    }

    // A ray cannot need more triangles than the mesh has, a few times over, so
    // anything past that is going in circles rather than getting anywhere.
    const int perRayCap = 4 * static_cast<int>(mesh->triangles.size()) + 64;

    // The scale the tolerances in here are quoted in.
    double meanEdge = 0.0;
    for (const auto &e : mesh->edges) meanEdge += normP(mesh->vertices[e[1]] - mesh->vertices[e[0]]);
    if (!mesh->edges.empty()) meanEdge /= static_cast<double>(mesh->edges.size());
    const double arrivalDist = arrivalTol * meanEdge;

    // How far a ray has to have gone before any of that applies.
    //
    // In a polysquare that closes up, a ray is a straight line in the parameter
    // domain from a corner to a wall, so it cannot be longer than the domain is
    // wide; the map is near isometric, so nor can its image in the model. A ray
    // past twice the bounding box diagonal is therefore not a long ray, it is a
    // ray in a parameterization that is not a polysquare -- and that, not the
    // distance to a corner, is what separates the ray that should be stopped
    // from the ones that are merely passing close to one.
    Point lo{1e300, 1e300}, hi{-1e300, -1e300};
    for (const auto &p : mesh->vertices) {
        lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
        hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
    }
    const double lostAfter = 2.0 * normP(hi - lo);

    // The boundary vertex a ray is standing on, if it is standing on one. The
    // tolerance is against the triangle it is in rather than the model, since
    // that is the scale the arithmetic that got it here worked at; a ray that
    // lands this near a vertex arrived at it exactly and the distance is
    // rounding.
    auto boundaryVertexAt = [&](int f, const Point &p) -> int {
        const Triangle &t = mesh->triangles[f];
        double longest = 0.0;
        for (int i = 0; i < 3; ++i)
            longest = std::max(longest,
                               normP(mesh->vertices[t[(i + 1) % 3]] - mesh->vertices[t[i]]));
        for (int i = 0; i < 3; ++i) {
            if (!mesh->isBoundaryVertex[t[i]]) continue;
            if (normP(mesh->vertices[t[i]] - p) <= 1e-6 * longest) return t[i];
        }
        return -1;
    };

    // A boundary edge at that vertex, and which end of it the vertex is in the
    // edge's own orientation -- the same pair a crossing records, so that the
    // block structure can place the node on its loop either way.
    auto boundaryEdgeAt = [&](int v, double &alongOut) -> int {
        const auto &vt = mesh->vertexTriangles;
        for (int k = vt.rowPtr[v]; k < vt.rowPtr[v + 1]; ++k) {
            const int f = vt.colIdx[k];
            for (int i = 0; i < 3; ++i) {
                const int e = mesh->triangleEdges[f][i];
                if (e < 0 || !mesh->isBoundaryEdge[e]) continue;
                if (mesh->edges[e][0] == v) { alongOut = 0.0; return e; }
                if (mesh->edges[e][1] == v) { alongOut = 1.0; return e; }
            }
        }
        return -1;
    };

    for (size_t id = 0; id < bikes.size(); ++id) {
        Motorcycle &m = bikes[id];
        while (m.alive) {
        if (++m.steps > perRayCap) { m.alive = false; ++report_.ranOut; break; }

        double J[2][2];
        triangleJacobian(m.tri, J);
        const Point meshDir = normalizeP(solve2(J, m.dir));
        if (normP(meshDir) < 1e-12) { m.alive = false; ++report_.ranOut; break; }

        Point hit;
        int localEdge = -1;
        double along = 0.0;
        if (!exitPoint(m.tri, m.pos, meshDir, hit, localEdge, along)) {
            // A ray standing exactly on a boundary vertex has left the model,
            // through a vertex rather than through an edge. exitPoint() cannot
            // say so: the ray is on a corner of the triangle it is in, so the
            // two edges there hold the point it would be starting from and the
            // third is behind it, and there is no crossing to return.
            //
            // That is what an iso-line landing on the corner it was aimed at
            // looks like from in here, and it is the *intended* end of a trace:
            // in a polysquare the line out of a reflex corner runs into another
            // corner, and Polysquare::snapCorners() is what makes the two share
            // a coordinate so that it lands rather than slips past. Counting it
            // as a ray that got nowhere would report the good case as the bad
            // one -- and worse, leave the trace with a free end, which costs
            // the block structure a four-sided block right where it was most
            // nearly right.
            const int landed = boundaryVertexAt(m.tri, m.pos);
            if (landed >= 0) {
                double landedAlong = 0.0;
                const int be = boundaryEdgeAt(landed, landedAlong);
                if (be >= 0) {
                    exitEdge[id] = be;
                    exitAlong[id] = landedAlong;
                }
                // Exactly on the vertex, so the node the block structure makes
                // here coincides with the corner instead of sitting a rounding
                // error off it.
                m.pos = mesh->vertices[landed];
                if (!traces[id].empty()) traces[id].back() = m.pos;

                // Seal the corner it stopped at, for the same reason launch()
                // seals the one it started from: the ray's wall is the edges it
                // crossed, and it crossed none at this end, so without this the
                // flood walks from one side of the ray to the other through the
                // triangle it stopped in.
                for (int i = 0; i < 3; ++i) {
                    const int ed = mesh->triangleEdges[m.tri][i];
                    if (ed < 0) continue;
                    if (mesh->edges[ed][0] == landed || mesh->edges[ed][1] == landed)
                        tracedEdges.insert(EdgeKey(mesh->edges[ed][0], mesh->edges[ed][1]));
                }

                m.alive = false;
                ++report_.reachedBoundary;
                break;
            }
            m.alive = false;
            ++report_.ranOut;
            break;
        }

        // Arriving at the corner another ray was launched from.
        //
        // In an exact polysquare the iso-line out of one reflex corner runs
        // into another one, and once Polysquare::snapCorners() has put the two
        // on a common coordinate that is what happens: the ray lands on the
        // corner and the vertex test above ends it. Where the parameterization
        // is not exact it lands a little to one side instead, slips past, and
        // keeps going -- on geom035 for 7428 steps and 126 units of a model two
        // units across, crossing the other lines 216 times on the way.
        //
        // So a ray that passes close enough to one of those corners is taken to
        // have arrived at it. The catch is that "close enough" does not
        // separate the two cases on its own, and it is worth recording why,
        // because the test looks obvious until it is measured: rays that are
        // doing exactly what they should also pass launch corners in mid
        // flight, and closer than the lost ones do. geom008's ray 5 passes one
        // at 3.2% of a mesh edge with a fifth of its length still to run,
        // geom009's ray 9 at 40%, geom011's ray 0 at 35% -- while the first
        // useful pass on geom035 is at 12.8%. The populations overlap the wrong
        // way round. Direction does not help either: a corner launches two rays
        // ninety degrees apart, so one of them is aligned with whatever is
        // arriving, whatever that is.
        //
        // What does separate them is how far the ray has already come, which is
        // what `lostAfter` is. A ray in a polysquare that closes up cannot be
        // longer than the domain is wide; one that is has established that it
        // is not in a polysquare, and nothing it does afterwards is worth
        // preserving. So the snap is only offered to rays past that point, and
        // the ones that are merely passing a corner never reach it.
        if (arrivalTol > 0.0 && m.distance > lostAfter) {
            int landedOn = -1;
            Point landedAt{0.0, 0.0};
            double bestT = 2.0;
            for (const auto &other : bikes) {
                if (other.originVertex < 0 || other.originVertex == m.originVertex) continue;
                const Point &p = mesh->vertices[other.originVertex];
                const Point d = hit - m.pos;
                const double dd = dotP(d, d);
                if (dd < 1e-30) continue;
                double t = dotP(p - m.pos, d) / dd;
                t = std::min(1.0, std::max(0.0, t));
                if (normP(p - (m.pos + d * t)) > arrivalDist) continue;
                if (t < bestT) { bestT = t; landedOn = other.originVertex; landedAt = p; }
            }
            if (landedOn >= 0) {
                Segment last;
                last.a = m.pos;
                last.b = landedAt;
                last.ua = imageOfPoint(m.tri, m.pos);
                last.ub = imageOfPoint(m.tri, landedAt);
                last.tri = m.tri;
                last.ray = static_cast<int>(id);
                last.step = m.steps - 1;
                if (normP(last.b - last.a) > 1e-15) segments.push_back(last);

                m.distance += normP(landedAt - m.pos);
                m.pos = landedAt;
                traces[id].push_back(landedAt);

                double landedAlong = 0.0;
                const int be = boundaryEdgeAt(landedOn, landedAlong);
                if (be >= 0) {
                    exitEdge[id] = be;
                    exitAlong[id] = landedAlong;
                }
                // Seal the corner it stopped at, where the triangle it stopped
                // in has it -- the same closing launch() does at the other end.
                for (int i = 0; i < 3; ++i) {
                    const int ed = mesh->triangleEdges[m.tri][i];
                    if (ed < 0) continue;
                    if (mesh->edges[ed][0] == landedOn || mesh->edges[ed][1] == landedOn)
                        tracedEdges.insert(EdgeKey(mesh->edges[ed][0], mesh->edges[ed][1]));
                }

                m.alive = false;
                ++report_.reachedBoundary;
                ++report_.arrivals;
                break;
            }
        }

        const int e = mesh->triangleEdges[m.tri][localEdge];
        const EdgeKey key(mesh->edges[e][0], mesh->edges[e][1]);

        Segment seg;
        seg.a = m.pos;
        seg.b = hit;
        seg.ua = imageOfPoint(m.tri, m.pos);
        seg.ub = imageOfPoint(m.tri, hit);
        seg.tri = m.tri;
        seg.ray = static_cast<int>(id);
        seg.step = m.steps - 1;
        segments.push_back(seg);

        m.distance += normP(hit - m.pos);
        m.pos = hit;
        traces[id].push_back(hit);

        // An iso-line does not stop where it meets another one. It is drawn
        // from the corner it belongs to all the way out of the model, and the
        // two simply cross. That is what leaves the block structure free of
        // T-junctions, which is the property Sec. 5 opens by claiming.
        tracedEdges.insert(key);

        const int next = mesh->triangleAdjacency[m.tri][localEdge];
        if (next < 0) {                 // out at the boundary of the model
            // Where on the boundary, in the edge's own orientation: exitPoint()
            // measured `along` from the triangle's corner localEdge, which need
            // not be the edge's first vertex.
            exitEdge[id] = e;
            exitAlong[id] = (mesh->edges[e][0] == mesh->triangles[m.tri][localEdge])
                                ? along : 1.0 - along;
            m.alive = false;
            ++report_.reachedBoundary;
            break;
        }

        // Carry the direction into the next triangle. The shared edge has an
        // image in each of them; across a cut those differ by the transition,
        // which is the rotation taking one to the other.
        const int va = mesh->edges[e][0], vb = mesh->edges[e][1];
        auto imageOf = [&](int f, int origVertex) -> Point {
            const Triangle &t = mesh->triangles[f];
            for (int i = 0; i < 3; ++i) if (t[i] == origVertex) return (*uv)[cut->triangles[f][i]];
            return Point{0.0, 0.0};
        };
        const Point eHere = imageOf(m.tri, vb) - imageOf(m.tri, va);
        const Point eThere = imageOf(next, vb) - imageOf(next, va);
        if (normP(eHere) > 1e-15 && normP(eThere) > 1e-15) {
            const double turn = computeAngle(eThere) - computeAngle(eHere);
            const double c = std::cos(turn), s = std::sin(turn);
            m.dir = Point{c * m.dir[0] - s * m.dir[1], s * m.dir[0] + c * m.dir[1]};
        }

        onWall[next] = 1;
        m.tri = next;
        }
    }

    // An iso-line drawn from both of its ends is one iso-line.
    //
    // Two reflex corners that face each other across the domain each send a ray
    // at the other, and they sit on a common coordinate, so the two rays run
    // down the same line in opposite directions and each finishes on the corner
    // the other started from. Before Polysquare::snapCorners() that could only
    // happen by luck -- the two corners were a rounding error apart and the
    // rays slid past each other -- so both were kept and neither was right.
    // Now it is how the pair is supposed to end, and keeping both would put two
    // coincident walls in the block structure: a pair of sides that lie on top
    // of each other with a block of no width in between. So the second one goes.
    //
    // The test is that the two run between the same two points the opposite way
    // round, are the same length, and pass through the same midpoint. Two
    // different curves can share their endpoints -- they would bound a lens
    // between them -- and the third condition is what tells that apart from one
    // curve seen twice.
    const double dupTol = 1e-6 * std::max(meanEdge, 1e-300);

    auto traceLength = [](const std::vector<Point> &t) {
        double l = 0.0;
        for (size_t i = 1; i < t.size(); ++i) l += normP(t[i] - t[i - 1]);
        return l;
    };
    auto traceMidpoint = [](const std::vector<Point> &t, double half) {
        double l = 0.0;
        for (size_t i = 1; i < t.size(); ++i) {
            const double s = normP(t[i] - t[i - 1]);
            if (l + s >= half) {
                const double f = (s > 0.0) ? (half - l) / s : 0.0;
                return t[i - 1] * (1.0 - f) + t[i] * f;
            }
            l += s;
        }
        return t.back();
    };

    std::vector<char> duplicate(traces.size(), 0);
    for (size_t i = 0; i < traces.size(); ++i) {
        if (duplicate[i] || traces[i].size() < 2) continue;
        const double li = traceLength(traces[i]);
        for (size_t j = i + 1; j < traces.size(); ++j) {
            if (duplicate[j] || traces[j].size() < 2) continue;
            if (normP(traces[i].front() - traces[j].back()) > dupTol) continue;
            if (normP(traces[i].back() - traces[j].front()) > dupTol) continue;
            const double lj = traceLength(traces[j]);
            if (std::fabs(li - lj) > dupTol) continue;
            if (normP(traceMidpoint(traces[i], 0.5 * li) - traceMidpoint(traces[j], 0.5 * lj)) >
                dupTol) continue;

            duplicate[j] = 1;
            ++report_.duplicates;
            if (exitEdge[j] >= 0) --report_.reachedBoundary; else --report_.ranOut;
        }
    }

    if (report_.duplicates > 0) {
        for (size_t j = 0; j < traces.size(); ++j) {
            if (!duplicate[j]) continue;
            traces[j].clear();
            exitEdge[j] = -1;
        }
        // The walls the dropped ray put down stay: its twin crossed the same
        // edges, so tracedEdges is already what it should be.
        std::vector<Segment> keep;
        keep.reserve(segments.size());
        for (const auto &s : segments)
            if (s.ray < 0 || s.ray >= static_cast<int>(duplicate.size()) || !duplicate[s.ray])
                keep.push_back(s);
        segments.swap(keep);
    }

    report_.tracedEdges = static_cast<int>(tracedEdges.size());
}

// ---------------------------------------------------------------------------
// floodBlocks()  --  the blocks are what is left between the trails
// ---------------------------------------------------------------------------
void MotorcycleGraph::floodBlocks() {
    const int nT = static_cast<int>(mesh->triangles.size());
    blockOfTriangle.assign(nT, -1);

    int block = 0;
    std::vector<int> sizes;
    for (int seed = 0; seed < nT; ++seed) {
        if (blockOfTriangle[seed] >= 0) continue;

        int count = 0;
        std::deque<int> queue{seed};
        blockOfTriangle[seed] = block;
        while (!queue.empty()) {
            const int f = queue.front();
            queue.pop_front();
            ++count;

            for (int i = 0; i < 3; ++i) {
                const int nb = mesh->triangleAdjacency[f][i];
                if (nb < 0 || blockOfTriangle[nb] >= 0) continue;
                const int e = mesh->triangleEdges[f][i];
                if (e < 0) continue;
                // A cut is not a wall: it is interior to the model, and the
                // parameterization is seamless across it.
                if (tracedEdges.count(EdgeKey(mesh->edges[e][0], mesh->edges[e][1]))) continue;
                blockOfTriangle[nb] = block;
                queue.push_back(nb);
            }
        }
        (void)count;
        ++block;
    }

    // The walls are mesh edges but the rays are not, so the flood pinches off
    // regions that are not blocks of the structure. Two kinds: where two rays
    // cross, the triangle holding the crossing has every edge walled off and
    // comes out as a block of one; and where a ray runs along a boundary
    // segment it happens to be collinear with, the strip between them is a
    // block of zero width that the mesh renders as one or two triangles.
    //
    // Both are the same statement -- a block the mesh cannot resolve -- and
    // both are caught by measuring the block in the parameter domain rather
    // than counting its triangles. On data/meshes the two populations are far
    // apart: the unresolved ones are 1.0 to 1.9 image edges across and the
    // real ones start at 3.4, so the cut at two edges is not delicate.
    double avgImageEdge = 0.0;
    int edgeCount = 0;
    for (int f = 0; f < nT; ++f) {
        for (int i = 0; i < 3; ++i) {
            const Point &a = (*uv)[cut->triangles[f][i]];
            const Point &b = (*uv)[cut->triangles[f][(i + 1) % 3]];
            avgImageEdge += normP(b - a);
            ++edgeCount;
        }
    }
    if (edgeCount > 0) avgImageEdge /= edgeCount;
    const double thinLimit = 2.0 * avgImageEdge;

    std::vector<int> remap(block);
    for (int b = 0; b < block; ++b) remap[b] = b;
    auto root = [&](int b) { while (remap[b] != b) b = remap[b] = remap[remap[b]]; return b; };

    for (int pass = 0; pass < 8; ++pass) {
        // The extent of each block in the parameter domain.
        std::vector<double> minU(block, std::numeric_limits<double>::max());
        std::vector<double> maxU(block, -std::numeric_limits<double>::max());
        std::vector<double> minV(block, std::numeric_limits<double>::max());
        std::vector<double> maxV(block, -std::numeric_limits<double>::max());
        for (int f = 0; f < nT; ++f) {
            const int b = root(blockOfTriangle[f]);
            for (int i = 0; i < 3; ++i) {
                const Point &q = (*uv)[cut->triangles[f][i]];
                minU[b] = std::min(minU[b], q[0]); maxU[b] = std::max(maxU[b], q[0]);
                minV[b] = std::min(minV[b], q[1]); maxV[b] = std::max(maxV[b], q[1]);
            }
        }

        bool changed = false;
        for (int f = 0; f < nT; ++f) {
            const int b = root(blockOfTriangle[f]);
            if (std::min(maxU[b] - minU[b], maxV[b] - minV[b]) >= thinLimit) continue;

            int best = -1;
            double bestSize = -1.0;
            for (int i = 0; i < 3; ++i) {
                const int nb = mesh->triangleAdjacency[f][i];
                if (nb < 0) continue;
                const int ob = root(blockOfTriangle[nb]);
                if (ob == b) continue;
                const double size = std::min(maxU[ob] - minU[ob], maxV[ob] - minV[ob]);
                if (size > bestSize) { bestSize = size; best = ob; }
            }
            if (best < 0) continue;
            remap[b] = best;
            changed = true;
            ++report_.mergedSlivers;
        }
        if (!changed) break;
    }

    // Re-index so the blocks are 0..n-1 again.
    std::vector<int> compact(block, -1);
    int nBlocks = 0;
    for (int f = 0; f < nT; ++f) {
        const int b = root(blockOfTriangle[f]);
        if (compact[b] < 0) compact[b] = nBlocks++;
        blockOfTriangle[f] = compact[b];
    }

    std::vector<int> finalSizes(nBlocks, 0);
    for (int f = 0; f < nT; ++f) ++finalSizes[blockOfTriangle[f]];

    report_.blocks = nBlocks;
    report_.smallestBlock =
        finalSizes.empty() ? 0 : *std::min_element(finalSizes.begin(), finalSizes.end());
}

// ---------------------------------------------------------------------------
// imageOfPoint()
// ---------------------------------------------------------------------------
Point MotorcycleGraph::imageOfPoint(int f, const Point &p) const {
    const Triangle &t = mesh->triangles[f];
    const Point &p0 = mesh->vertices[t[0]];
    const Point &p1 = mesh->vertices[t[1]];
    const Point &p2 = mesh->vertices[t[2]];
    const double area2 = cross2(p1 - p0, p2 - p0);
    if (std::fabs(area2) < 1e-18) return (*uv)[cut->triangles[f][0]];

    const double l0 = cross2(p1 - p, p2 - p) / area2;
    const double l1 = cross2(p2 - p, p0 - p) / area2;
    const double l2 = 1.0 - l0 - l1;

    const Point &q0 = (*uv)[cut->triangles[f][0]];
    const Point &q1 = (*uv)[cut->triangles[f][1]];
    const Point &q2 = (*uv)[cut->triangles[f][2]];
    return Point{l0 * q0[0] + l1 * q1[0] + l2 * q2[0],
                 l0 * q0[1] + l1 * q1[1] + l2 * q2[1]};
}

// ---------------------------------------------------------------------------
// findNodes()  --  where the edges of the block structure meet
//
// Two iso-lines that cross do so inside one triangle, so the crossings are
// found by intersecting the steps that share a triangle rather than by any
// search over the plane. That also keeps the answer right where the image
// overlaps itself: two lines that pass over each other on different sheets
// share no triangle and are not counted.
// ---------------------------------------------------------------------------
void MotorcycleGraph::findNodes() {
    nodes.clear();

    // Corners of the polysquare, including the convex ones, which launch no
    // ray but bound a block just the same.
    const auto &corners = poly->getBoundaryCorners();
    for (int v : mesh->boundaryVertices) {
        if (v >= static_cast<int>(corners.size()) || corners[v] == 0) continue;
        const auto &vt = mesh->vertexTriangles;
        if (vt.rowPtr[v] >= vt.rowPtr[v + 1]) continue;
        Node n;
        n.xy = mesh->vertices[v];
        n.uv = imageOfPoint(vt.colIdx[vt.rowPtr[v]], mesh->vertices[v]);
        n.kind = Node::Corner;
        n.vertex = v;
        nodes.push_back(n);
    }

    // Where each line leaves the model. A ray that never got out has no node
    // here: its trace stops in the middle of the model, and pretending its last
    // point is a boundary end would hide that.
    for (size_t r = 0; r < traces.size(); ++r) {
        if (traces[r].size() < 2 || r >= exitEdge.size() || exitEdge[r] < 0) continue;
        Node n;
        n.xy = traces[r].back();
        n.kind = Node::BoundaryEnd;
        n.ray[0] = static_cast<int>(r);
        n.param[0] = static_cast<double>(traces[r].size() - 1);
        n.boundaryEdge = exitEdge[r];
        n.alongEdge = exitAlong[r];
        for (const auto &seg : segments)
            if (seg.ray == static_cast<int>(r) && seg.step == n.param[0] - 1) { n.uv = seg.ub; break; }
        nodes.push_back(n);
    }

    // Crossings, triangle by triangle.
    std::unordered_map<int, std::vector<int>> perTriangle;
    for (size_t i = 0; i < segments.size(); ++i) perTriangle[segments[i].tri].push_back(static_cast<int>(i));

    for (const auto &[f, list] : perTriangle) {
        for (size_t i = 0; i < list.size(); ++i) {
            for (size_t j = i + 1; j < list.size(); ++j) {
                const Segment &s1 = segments[list[i]];
                const Segment &s2 = segments[list[j]];
                if (s1.ray == s2.ray) continue;
                // Two rays out of one corner meet at that corner, which is a
                // node already and not a crossing of anything.
                if (normP(s1.a - s2.a) < 1e-12 || normP(s1.a - s2.b) < 1e-12 ||
                    normP(s1.b - s2.a) < 1e-12 || normP(s1.b - s2.b) < 1e-12) continue;

                const Point d1 = s1.b - s1.a, d2 = s2.b - s2.a;
                const double den = cross2(d1, d2);
                if (std::fabs(den) < 1e-18) continue;
                const Point w = s2.a - s1.a;
                const double t1 = cross2(w, d2) / den;
                const double t2 = cross2(w, d1) / den;
                if (t1 < -1e-9 || t1 > 1.0 + 1e-9 || t2 < -1e-9 || t2 > 1.0 + 1e-9) continue;

                Node n;
                n.xy = Point{s1.a[0] + t1 * d1[0], s1.a[1] + t1 * d1[1]};
                n.uv = Point{s1.ua[0] + t1 * (s1.ub[0] - s1.ua[0]),
                             s1.ua[1] + t1 * (s1.ub[1] - s1.ua[1])};
                n.kind = Node::Crossing;
                n.ray[0] = s1.ray;
                n.param[0] = s1.step + t1;
                n.ray[1] = s2.ray;
                n.param[1] = s2.step + t2;
                nodes.push_back(n);
            }
        }
    }

    report_.crossings = 0;
    for (const auto &n : nodes) if (n.kind == Node::Crossing) ++report_.crossings;
    report_.nodes = static_cast<int>(nodes.size());
}

void MotorcycleGraph::build() {
    launch();
    run();
    findNodes();
    floodBlocks();
}

// ---------------------------------------------------------------------------
// VTK output
// ---------------------------------------------------------------------------
bool MotorcycleGraph::writeBlocksVTU(const std::string &filename) const {
    if (blockOfTriangle.empty()) return false;
    std::ofstream out(filename);
    if (!out) return false;

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << mesh->vertices.size()
        << "\" NumberOfCells=\"" << mesh->triangles.size() << "\">\n";

    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto &v : mesh->vertices) out << "          " << v[0] << " " << v[1] << " 0.0\n";
    out << "        </DataArray>\n";
    out << "      </Points>\n";

    out << "      <CellData Scalars=\"block\">\n";
    out << "        <DataArray type=\"Int32\" Name=\"block\" format=\"ascii\">\n";
    for (int b : blockOfTriangle) out << "          " << b << "\n";
    out << "        </DataArray>\n";
    out << "      </CellData>\n";

    out << "      <Cells>\n";
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto &t : mesh->triangles) out << "          " << t[0] << " " << t[1] << " " << t[2] << "\n";
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (size_t i = 1; i <= mesh->triangles.size(); ++i) out << "          " << i * 3 << "\n";
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < mesh->triangles.size(); ++i) out << "          5\n";
    out << "        </DataArray>\n";
    out << "      </Cells>\n";

    out << "    </Piece>\n";
    out << "  </UnstructuredGrid>\n";
    out << "</VTKFile>\n";
    return true;
}

bool MotorcycleGraph::writeTracesVTU(const std::string &filename) const {
    std::ofstream out(filename);
    if (!out) return false;

    size_t nPoints = 0, nCells = 0;
    for (const auto &tr : traces) {
        if (tr.size() < 2) continue;
        nPoints += tr.size();
        ++nCells;
    }

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << nCells << "\">\n";

    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto &tr : traces) {
        if (tr.size() < 2) continue;
        for (const auto &p : tr) out << "          " << p[0] << " " << p[1] << " 0.0\n";
    }
    out << "        </DataArray>\n";
    out << "      </Points>\n";

    out << "      <CellData Scalars=\"motorcycle\">\n";
    out << "        <DataArray type=\"Int32\" Name=\"motorcycle\" format=\"ascii\">\n";
    for (size_t i = 0; i < traces.size(); ++i) {
        if (traces[i].size() < 2) continue;
        out << "          " << i << "\n";
    }
    out << "        </DataArray>\n";
    out << "      </CellData>\n";

    out << "      <Cells>\n";
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    size_t base = 0;
    for (const auto &tr : traces) {
        if (tr.size() < 2) continue;
        out << "         ";
        for (size_t i = 0; i < tr.size(); ++i) out << " " << base + i;
        out << "\n";
        base += tr.size();
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    size_t off = 0;
    for (const auto &tr : traces) {
        if (tr.size() < 2) continue;
        off += tr.size();
        out << "          " << off << "\n";
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < nCells; ++i) out << "          4\n"; // VTK_POLY_LINE
    out << "        </DataArray>\n";
    out << "      </Cells>\n";

    out << "    </Piece>\n";
    out << "  </UnstructuredGrid>\n";
    out << "</VTKFile>\n";
    return true;
}
