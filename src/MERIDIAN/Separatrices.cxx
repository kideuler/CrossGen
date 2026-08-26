#include "Separatrices.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#include "MERIDIAN/LayoutEnergy.hxx"

namespace {

inline int localOf(const Triangle &t, int v) {
    if (t[0] == v) return 0;
    if (t[1] == v) return 1;
    if (t[2] == v) return 2;
    return -1;
}

// Distance from q to the segment [a, b], and where along it that was.
double segmentClosest(const Point &q, const Point &a, const Point &b, double &s) {
    const Point ab = b - a;
    const double denom = dotP(ab, ab);
    s = 0.0;
    if (denom > 0.0) s = std::max(0.0, std::min(1.0, dotP(q - a, ab) / denom));
    return normP(q - (a + ab * s));
}

} // namespace

// ---------------------------------------------------------------------------
Separatrices::Separatrices(const Immersion &immersion, const std::vector<Point> &map)
    : Separatrices(immersion, map, Options()) {}

Separatrices::Separatrices(const Immersion &immersion, const std::vector<Point> &map,
                           const Options &opts)
    : imm(&immersion), options(opts), uv(map) {
    if (uv.size() != immersion.getCutMesh().vertices.size()) {
        throw std::runtime_error("Separatrices: the map has a different number of points "
                                 "than Omega has vertices");
    }
    buildTables();
    buildEmitters();
    measureExtent();
    traceAll();
    check();
}

Separatrices::Separatrices(const LayoutEnergy &layout)
    : Separatrices(layout.getImmersion(), layout.getUV(), Options()) {}

Separatrices::Separatrices(const LayoutEnergy &layout, const Options &opts)
    : Separatrices(layout.getImmersion(), layout.getUV(), opts) {}

// ---------------------------------------------------------------------------
// buildTables()
//
// The same three lookups Stage 5 builds -- Omega's edges by their endpoints,
// which of its boundary edges came from dS rather than from a cut, and the seam
// pairing of the ones that came from a cut -- plus one this stage needs and
// Stage 5 did without: the seam partner of a *vertex*.
//
// A curve that meets an arc of G part way along an edge crosses it by the edge
// pairing. One that runs into the arc head-on, at a vertex, has no edge to
// cross and has to be carried over by the pairing of the vertex itself. Stage 5
// gave up there, which is fair for a seeder looking for near-misses; a
// separatrix that stopped at a seam vertex would be reported as violating Q5
// when it does no such thing.
// ---------------------------------------------------------------------------
void Separatrices::buildTables() {
    const Mesh &cm = imm->getCutMesh();
    const Mesh &om = imm->getOriginalMesh();
    const std::vector<int> &c2o = imm->getCut().getCutVertexToOriginal();

    cutEdgeIndex.reserve(cm.edges.size() * 2);
    for (int e = 0; e < static_cast<int>(cm.edges.size()); ++e) {
        cutEdgeIndex.emplace(EdgeKey(cm.edges[e][0], cm.edges[e][1]), e);
    }

    std::unordered_map<EdgeKey, int, EdgeKeyHash> origEdgeIndex;
    origEdgeIndex.reserve(om.edges.size() * 2);
    for (int e = 0; e < static_cast<int>(om.edges.size()); ++e) {
        origEdgeIndex.emplace(EdgeKey(om.edges[e][0], om.edges[e][1]), e);
    }

    parentOnBoundary.assign(cm.edges.size(), 0);
    vertexOnRealBoundary.assign(cm.vertices.size(), 0);
    for (int e = 0; e < static_cast<int>(cm.edges.size()); ++e) {
        auto it = origEdgeIndex.find(EdgeKey(c2o[cm.edges[e][0]], c2o[cm.edges[e][1]]));
        if (it != origEdgeIndex.end() && om.isBoundaryEdge[it->second]) {
            parentOnBoundary[e] = 1;
            vertexOnRealBoundary[cm.edges[e][0]] = 1;
            vertexOnRealBoundary[cm.edges[e][1]] = 1;
        }
    }

    const auto &pairs = imm->getSeamPairs();
    seamSide.reserve(pairs.size() * 4);
    for (int i = 0; i < static_cast<int>(pairs.size()); ++i) {
        seamSide[EdgeKey(pairs[i][0], pairs[i][1])] = 2 * i + 0;   // the plus child
        seamSide[EdgeKey(pairs[i][2], pairs[i][3])] = 2 * i + 1;   // the minus child
    }

    // The vertex pairing, read off the arcs rather than the edge pairs, so that
    // the k that goes with it is the arc's. A vertex carried by two arcs at
    // once is a junction of G: there is no single rotation to cross it by, and
    // it is marked rather than guessed at.
    for (const Immersion::Arc &a : imm->getArcs()) {
        if (a.degenerate) continue;
        const int k = ((a.k % 4) + 4) % 4;
        for (size_t i = 0; i < a.plusChain.size() && i < a.minusChain.size(); ++i) {
            const int pv = a.plusChain[i], mv = a.minusChain[i];
            if (pv == mv) continue;
            for (int side = 0; side < 2; ++side) {
                const int from = side == 0 ? pv : mv;
                const int to = side == 0 ? mv : pv;
                auto it = vertexSeam.find(from);
                if (it == vertexSeam.end()) {
                    VertexSeam vs;
                    vs.partner = to;
                    vs.k = k;
                    vs.sign = side == 0 ? +1 : -1;
                    vertexSeam.emplace(from, vs);
                } else if (it->second.partner != to || it->second.k != k) {
                    it->second.ambiguous = true;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// buildEmitters()
//
// Stage 4's cones, or -- when there are none, which is footnote 3's annulus --
// one regular vertex pressed into service as the surface's only singularity.
// ---------------------------------------------------------------------------
void Separatrices::buildEmitters() {
    const Mesh &cm = imm->getCutMesh();
    const Mesh &om = imm->getOriginalMesh();

    coneVertex = imm->getConeVertices();
    coneIndex = imm->getConeIndices();
    coneChildren = imm->getConeChildren();
    vertCone = imm->getCutVertexCone();
    if (vertCone.size() != cm.vertices.size()) vertCone.assign(cm.vertices.size(), -1);

    coneOnBoundary.assign(coneVertex.size(), 0);
    for (size_t i = 0; i < coneVertex.size(); ++i) {
        const int v = coneVertex[i];
        coneOnBoundary[i] = (v >= 0 && v < static_cast<int>(om.isBoundaryVertex.size()) &&
                             om.isBoundaryVertex[v]) ? 1 : 0;
    }
    if (!coneVertex.empty()) return;

    // Footnote 3. The point has to be one whose fan is a full 2pi in Omega --
    // interior to S and untouched by the cutting graph -- or the four rays it
    // is supposed to emit are not all there.
    Point centre{0.0, 0.0};
    for (const Point &p : uv) centre = centre + p;
    if (!uv.empty()) centre = centre / static_cast<double>(uv.size());

    int best = -1;
    double bestD = std::numeric_limits<double>::infinity();
    for (int v = 0; v < static_cast<int>(cm.vertices.size()); ++v) {
        if (vertexOnRealBoundary[v] || vertexSeam.count(v)) continue;
        const int ov = imm->getCut().getCutVertexToOriginal()[v];
        if (ov >= 0 && ov < static_cast<int>(om.isBoundaryVertex.size()) &&
            om.isBoundaryVertex[ov]) continue;
        const double d = normP(uv[v] - centre);
        if (d < bestD) { bestD = d; best = v; }
    }
    if (best < 0) {
        report.messages.push_back(
            "There are no cones and no regular interior point to stand in for one, so no "
            "separatrix could be emitted (footnote 3 does not apply to this mesh).");
        return;
    }

    coneVertex.push_back(imm->getCut().getCutVertexToOriginal()[best]);
    coneIndex.push_back(0);
    coneChildren.push_back(std::vector<int>{best});
    coneOnBoundary.push_back(0);
    vertCone[best] = 0;
    report.designatedCone = coneVertex.back();
    report.messages.push_back(
        "The cone set is empty, so footnote 3 was applied: regular vertex " +
        std::to_string(coneVertex.back()) +
        " was called the surface's only singularity and emitted four separatrices.");
}

// ---------------------------------------------------------------------------
void Separatrices::measureExtent() {
    Point lo{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    Point hi{-lo[0], -lo[1]};
    for (const Point &p : uv) {
        lo[0] = std::min(lo[0], p[0]); hi[0] = std::max(hi[0], p[0]);
        lo[1] = std::min(lo[1], p[1]); hi[1] = std::max(hi[1], p[1]);
    }
    extent = uv.empty() ? 1.0 : std::hypot(hi[0] - lo[0], hi[1] - lo[1]);
    if (!(extent > 0.0)) extent = 1.0;
    snapTol = options.coneSnapTolerance * extent;
    report.extent = extent;
    report.snapTolerance = snapTol;

    // The same for S itself. The Ricci metric is only defined up to a scale, so
    // the image and the model can be any factor apart, and a tolerance measured
    // on one of them means nothing on the other.
    Point mlo{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    Point mhi{-mlo[0], -mlo[1]};
    for (const Point &p : imm->getCutMesh().vertices) {
        mlo[0] = std::min(mlo[0], p[0]); mhi[0] = std::max(mhi[0], p[0]);
        mlo[1] = std::min(mlo[1], p[1]); mhi[1] = std::max(mhi[1], p[1]);
    }
    modelExtent = std::hypot(mhi[0] - mlo[0], mhi[1] - mlo[1]);
    if (!(modelExtent > 0.0)) modelExtent = 1.0;
}

// ---------------------------------------------------------------------------
// sweepCone()
//
// The one-ring fan of a cone in Omega, in order, as a list of angular sectors.
//
// The walk is counter-clockwise, which for a positively oriented image triangle
// means: at the corner where the cone sits, the sector runs from the next
// vertex round to the one after, and the sector after this one is across the
// second of those two edges. Stage 4 normalises the immersion so that every
// image triangle is positively oriented and Stage 6's barrier keeps it that
// way, so the two are the same walk.
//
// Three things end it, and only one of them is an end:
//
//   * the next edge is a boundary edge of Omega whose parent is in dS. The cone
//     is on the boundary of the surface and its fan is the interior angle
//     between two boundary curves; the sweep stops.
//   * the next edge is a seam. The fan continues on the far side of the cut, at
//     the partner child of this vertex, and the accumulated angle runs on
//     through the crossing. This is the case Fig. 7 draws: one cone, three
//     children, one fan.
//   * the walk comes back to the sector it started in. The cone is interior to
//     S and its fan has closed.
// ---------------------------------------------------------------------------
bool Separatrices::sweepCone(int slot, std::vector<Sector> &out, bool &closed) const {
    const Mesh &cm = imm->getCutMesh();
    out.clear();
    closed = false;

    // Where the sweep starts: the sector with no counter-clockwise predecessor.
    // At a cone on dS there are two such sectors when the cutting graph also
    // ends there, and the one wanted is the one whose *first* edge is the
    // boundary of the surface.
    int startFace = -1, startChild = -1;
    for (int child : coneChildren[slot]) {
        if (child < 0 || child >= static_cast<int>(cm.vertices.size())) continue;
        const auto &vt = cm.vertexTriangles;
        for (int i = vt.rowPtr[child]; i < vt.rowPtr[child + 1]; ++i) {
            const int f = vt.colIdx[i];
            const int lc = localOf(cm.triangles[f], child);
            if (lc < 0) continue;
            if (cm.triangleAdjacency[f][lc] >= 0) continue;   // has a predecessor
            const int e = cm.triangleEdges[f][lc];
            if (coneOnBoundary[slot] && !parentOnBoundary[e]) continue;
            startFace = f;
            startChild = child;
            break;
        }
        if (startFace >= 0) break;
    }
    if (startFace < 0) {
        // No boundary edge at all in the fan: an interior vertex whose fan is a
        // closed 2pi ring, which is footnote 3's designated point.
        for (int child : coneChildren[slot]) {
            const auto &vt = cm.vertexTriangles;
            if (child < 0 || child >= static_cast<int>(cm.vertices.size())) continue;
            if (vt.rowPtr[child + 1] > vt.rowPtr[child]) {
                startFace = vt.colIdx[vt.rowPtr[child]];
                startChild = child;
                break;
            }
        }
    }
    if (startFace < 0) return false;

    int f = startFace, child = startChild;
    double acc = 0.0;
    const int cap = 4096;

    for (int guard = 0; guard < cap; ++guard) {
        const Triangle &t = cm.triangles[f];
        const int lc = localOf(t, child);
        if (lc < 0) return false;

        const int fromV = t[(lc + 1) % 3];
        const int toV = t[(lc + 2) % 3];
        const Point A = uv[fromV] - uv[child];
        const Point B = uv[toV] - uv[child];
        const double cr = cross2(A, B), dt = dotP(A, B);
        if (normP(A) <= 0.0 || normP(B) <= 0.0) return false;

        Sector s;
        s.face = f;
        s.child = child;
        s.lc = lc;
        s.startAngle = std::atan2(A[1], A[0]);
        // A positively oriented image triangle has cr > 0 and an angle in
        // (0, pi). A non-positive cr means Stage 6 came back with an inverted
        // triangle, which is Q1 failing; the magnitude is still the best
        // available reading of the angle, and the caller counts the failure.
        s.angle = std::fabs(std::atan2(cr, dt));
        s.acc = acc;
        acc += s.angle;
        out.push_back(s);

        const int endLocal = (lc + 2) % 3;
        const int nb = cm.triangleAdjacency[f][endLocal];
        if (nb >= 0) {
            f = nb;
            if (f == startFace && child == startChild) { closed = true; return true; }
            continue;
        }

        const int e = cm.triangleEdges[f][endLocal];
        if (parentOnBoundary[e]) return true;   // a cone on dS: the sweep is done

        // A seam. Cross to the partner child and carry on; the accumulated
        // angle does not care that the sheet changed, and the direction the
        // rays come out at is read from the sector they land in, so the
        // rotation by k quarter turns is picked up for free.
        auto it = seamSide.find(EdgeKey(child, toV));
        if (it == seamSide.end()) return false;
        const auto &pr = imm->getSeamPairs()[it->second / 2];
        const bool onPlus = (it->second % 2) == 0;

        int partner = -1, partnerOther = -1;
        if (onPlus) {
            if (child == pr[0])      { partner = pr[2]; partnerOther = pr[3]; }
            else if (child == pr[1]) { partner = pr[3]; partnerOther = pr[2]; }
        } else {
            if (child == pr[2])      { partner = pr[0]; partnerOther = pr[1]; }
            else if (child == pr[3]) { partner = pr[1]; partnerOther = pr[0]; }
        }
        if (partner < 0) return false;

        auto ie = cutEdgeIndex.find(EdgeKey(partner, partnerOther));
        if (ie == cutEdgeIndex.end()) return false;
        const int te = ie->second;
        const int tf = cm.edgeTriangles[te][0] >= 0 ? cm.edgeTriangles[te][0]
                                                    : cm.edgeTriangles[te][1];
        if (tf < 0) return false;
        const int lp = localOf(cm.triangles[tf], partner);
        // The partner edge has to be the *first* edge of its sector, or the two
        // sides of the cut are glued with opposite orientations and the sweep
        // is not a sweep.
        if (lp < 0 || cm.triangles[tf][(lp + 1) % 3] != partnerOther) return false;

        f = tf;
        child = partner;
        if (f == startFace && child == startChild) { closed = true; return true; }
    }
    return false;
}

// ---------------------------------------------------------------------------
// traceAll()
//
// Sec. 3.4's ray count and ray placement, then the march.
//
//   interior cone   valence 4 - I, and all 4 - I rays are separatrices
//   boundary cone   valence 3 - I, of which two run along dS, leaving 1 - I
//
// The rays sit at the accumulated fan angles that put the direction of travel
// on an axis. Only the first one has to be solved for: the fan angle a and the
// image angle differ by a constant within a sheet and by a multiple of pi/2
// across a seam, so once the first axis direction is found the rest are a
// quarter turn apart all the way round.
// ---------------------------------------------------------------------------
void Separatrices::traceAll() {
    traced.clear();

    for (int slot = 0; slot < static_cast<int>(coneVertex.size()); ++slot) {
        std::vector<Sector> fan;
        bool closed = false;
        if (!sweepCone(slot, fan, closed) || fan.empty()) {
            ++report.fanFailures;
            std::ostringstream oss;
            oss << "The one-ring fan of the cone at vertex " << coneVertex[slot]
                << " could not be swept, so its separatrices were not emitted.";
            report.messages.push_back(oss.str());
            continue;
        }

        const int I = coneIndex[slot];
        const bool onBoundary = coneOnBoundary[slot] != 0;

        // A cone interior to S has a fan that closes; one on dS has a fan that
        // runs from one boundary curve to the other and does not. Anything else
        // means the sweep followed the cutting graph somewhere it should not
        // have, and the ray directions that come out of it are not to be
        // trusted.
        if (closed == onBoundary) {
            ++report.fanFailures;
            std::ostringstream oss;
            oss << "The fan of the cone at vertex " << coneVertex[slot] << " "
                << (closed ? "closed although the cone lies on dS"
                           : "did not close although the cone is interior to S")
                << "; no separatrices were emitted from it.";
            report.messages.push_back(oss.str());
            continue;
        }

        const double total = fan.back().acc + fan.back().angle;
        const double target = onBoundary ? (M_PI - M_PI_2 * I) : (2.0 * M_PI - M_PI_2 * I);
        const double residual = std::fabs(total - target);
        if (residual > report.maxConeAngleResidual) {
            report.maxConeAngleResidual = residual;
            report.worstConeAngleSlot = slot;
        }

        int rays = onBoundary ? (1 - I) : (4 - I);
        report.prescribed += std::max(0, rays);

        // Q2 is what makes the prescribed count the right count. Where the map
        // did not realise the prescribed angle, the fan is the only thing that
        // can be believed, so the count is taken from it instead and the
        // disagreement is named.
        if (residual > 0.05) {
            const int measured = static_cast<int>(std::lround(total / M_PI_2)) -
                                 (onBoundary ? 1 : 0);
            std::ostringstream oss;
            oss << "The cone at vertex " << coneVertex[slot] << " has index " << I
                << ", which asks for a cone angle of " << target << " rad, but its fan in "
                << "the layout measures " << total << " rad; ";
            if (measured == rays) {
                oss << "the ray count is unaffected at " << std::max(0, rays) << ".";
            } else {
                oss << std::max(0, measured) << " separatrices were emitted rather than the "
                    << std::max(0, rays) << " the index asks for.";
            }
            report.messages.push_back(oss.str());
            rays = measured;
        }
        if (rays <= 0) continue;   // a convex boundary corner emits nothing
        ++report.cones;

        // The first ray: the smallest fan angle at which the direction of
        // travel lands on an axis. At a boundary cone the two ends of the fan
        // are the two curves of dS that meet at the corner, and Q3 has put both
        // on an axis, so the ray at a = 0 has to be stepped over.
        double a0 = std::fmod(-fan.front().startAngle, M_PI_2);
        if (a0 < 0.0) a0 += M_PI_2;
        if (onBoundary) {
            while (a0 <= options.fanEndTolerance) a0 += M_PI_2;
        }

        for (int j = 0; j < rays; ++j) {
            const double a = a0 + j * M_PI_2;

            size_t si = 0;
            while (si + 1 < fan.size() && fan[si + 1].acc <= a) ++si;

            const Sector &s = fan[si];
            const double abs = s.startAngle + (a - s.acc);
            const int dir = ((static_cast<int>(std::lround(abs / M_PI_2)) % 4) + 4) % 4;

            Curve c = trace(slot, s.child, s.face, dir, j, a);
            ++report.emitted;
            switch (c.end) {
                case End::Cone:     ++report.endedAtCone; break;
                case End::Boundary: ++report.endedAtBoundary; break;
                case End::Capped:   ++report.capped; break;
                case End::Stuck:    ++report.stuck; break;
                case End::Degenerate: ++report.degenerate; break;
            }
            report.seamCrossings += c.seamCrossings;
            report.maxSeamCrossings = std::max(report.maxSeamCrossings, c.seamCrossings);
            report.triangleSteps += static_cast<long long>(c.steps.size());
            report.maxTriangleSteps = std::max(report.maxTriangleSteps,
                                               static_cast<int>(c.steps.size()));
            if (c.end == End::Cone) report.maxSnapGap = std::max(report.maxSnapGap, c.gap);
            traced.push_back(std::move(c));
        }
    }
}

// ---------------------------------------------------------------------------
// trace()  --  the marching of Sec. 3.4
//
// phi is affine on a triangle, so the isoline of the held coordinate is a
// straight segment and the exit point is a linear interpolation along one edge.
// The walk keeps the coordinate it entered the triangle with, which is what
// makes the curve exact rather than merely close: the exit point of one
// triangle is the entry point of the next and they share the constant
// coordinate, so nothing accumulates.
//
// Four things can happen at the far side of a triangle. Three of them are
// Q5's -- a cone, the boundary, or an arc of G to be continued across -- and
// the fourth is a vertex, which is not in the paper's list because in the
// smooth setting it does not arise. It arises here constantly: Q3 puts the
// boundary of the image on the axes, so an axis-parallel isoline runs into
// boundary vertices head-on. The walk turns through the fan at such a vertex,
// which is the discrete reading of "continue in the same direction".
// ---------------------------------------------------------------------------
Separatrices::Curve Separatrices::trace(int slot, int child, int face, int dirT,
                                        int ray, double fanAngle) const {
    const Mesh &cm = imm->getCutMesh();
    const auto &pairs = imm->getSeamPairs();
    const auto &pairArc = imm->getSeamPairArc();
    const auto &arcList = imm->getArcs();

    Curve cv;
    cv.cone = slot;
    cv.child = child;
    cv.index = coneIndex[slot];
    cv.boundaryCone = coneOnBoundary[slot] != 0;
    cv.ray = ray;
    cv.dir = ((dirT % 4) + 4) % 4;
    cv.fanAngle = fanAngle;
    cv.end = End::Capped;
    cv.nearestConeGap = std::numeric_limits<double>::infinity();

    const double tiny = 1e-12 * extent;
    const double vertexEps = 1e-9 * extent;

    int f = face;
    int entryEdge = -1;
    int dir = cv.dir;
    Point p = uv[child];
    std::array<double, 3> entryB = baryUnit(f, child);
    int fanSpins = 0;

    auto pushStep = [&](const std::array<double, 3> &a, const std::array<double, 3> &b) {
        Step st;
        st.face = f;
        st.entry = a;
        st.exit = b;
        cv.imageLength += normP(point(st, true, Space::Image) - point(st, false, Space::Image));
        cv.modelLength += normP(point(st, true, Space::Model) - point(st, false, Space::Model));
        cv.steps.push_back(st);
    };

    for (int step = 0; step < options.maxSteps; ++step) {
        const int tc = dir % 2;        // the coordinate the curve advances in
        const int cc = 1 - tc;         // the one it holds constant
        const double sgn = (dir < 2) ? 1.0 : -1.0;
        const double c0 = p[cc];

        const Triangle &t = cm.triangles[f];
        int best = -1;
        double bestT = 0.0, bestAdv = tiny;
        Point bestP{0.0, 0.0};

        for (int q = 0; q < 3; ++q) {
            if (q == entryEdge) continue;
            const int m = q, n = (q + 1) % 3;
            const double A = uv[t[m]][cc];
            const double B = uv[t[n]][cc];
            const double den = B - A;
            if (std::fabs(den) < 1e-300) continue;   // the edge lies on the isoline
            const double tt = (c0 - A) / den;
            if (tt < -1e-9 || tt > 1.0 + 1e-9) continue;
            const Point P = uv[t[m]] * (1.0 - tt) + uv[t[n]] * tt;
            const double adv = (P[tc] - p[tc]) * sgn;
            if (adv > bestAdv) {
                bestAdv = adv;
                best = q;
                bestT = std::min(1.0, std::max(0.0, tt));
                bestP = P;
            }
        }

        if (best >= 0) {
            // A crossing that lands on an endpoint is put exactly on it. The
            // curve then leaves through that vertex's fan rather than through
            // the edge, and the two have to agree about where the vertex is:
            // the recorded exit of this triangle and the recorded entry of the
            // next are the same point of S only if both are the vertex itself.
            const int m = best, n = (best + 1) % 3;
            if (normP(bestP - uv[t[m]]) <= vertexEps)      bestT = 0.0;
            else if (normP(bestP - uv[t[n]]) <= vertexEps) bestT = 1.0;
            bestP = uv[t[m]] * (1.0 - bestT) + uv[t[n]] * bestT;
        }

        if (best < 0) {
            // The curve is at a vertex of this triangle and leaves it through a
            // different triangle of that vertex's fan.
            int w = -1;
            for (int q = 0; q < 3; ++q) {
                if (normP(p - uv[t[q]]) <= vertexEps) { w = t[q]; break; }
            }
            if (w >= 0 && vertCone[w] >= 0 && !cv.steps.empty()) {
                cv.end = End::Cone;
                cv.toCone = vertCone[w];
                cv.toChild = w;
                cv.gap = 0.0;
                cv.endDir = dir;
                cv.nearestConeGap = 0.0;
                cv.nearestCone = vertCone[w];
                return cv;
            }

            const int nf = (w >= 0 && fanSpins < 64) ? faceContaining(w, dir, f) : -1;
            if (nf >= 0) {
                ++fanSpins;
                f = nf;
                entryEdge = -1;
                p = uv[w];
                entryB = baryUnit(f, w);
                continue;
            }
            if (w >= 0 && vertexOnRealBoundary[w]) {
                cv.end = End::Boundary;   // Q5 case 2, met at a vertex
                cv.endDir = dir;
                return cv;
            }
            // The fan ends at a seam: the curve met the cutting graph head-on
            // rather than crossing one of its edges. It carries over to the
            // partner child, turned by that arc's k, exactly as an edge
            // crossing would.
            if (w >= 0 && fanSpins < 64) {
                auto it = vertexSeam.find(w);
                if (it != vertexSeam.end() && !it->second.ambiguous) {
                    const int k = it->second.k;
                    const int nd = it->second.sign > 0 ? (((dir - k) % 4 + 4) % 4)
                                                       : (((dir + k) % 4 + 4) % 4);
                    const int g = faceContaining(it->second.partner, nd, -1);
                    if (g >= 0) {
                        ++fanSpins;
                        ++cv.seamCrossings;
                        f = g;
                        dir = nd;
                        entryEdge = -1;
                        p = uv[it->second.partner];
                        entryB = baryUnit(f, it->second.partner);
                        continue;
                    }
                }
            }
            cv.end = cv.steps.empty() ? End::Degenerate : End::Stuck;
            cv.endDir = dir;
            return cv;
        }
        fanSpins = 0;

        // Q5 case 1. Only the cones at the corners of the triangle being
        // crossed are candidates: Psi overlaps itself, so image distance alone
        // would match cones the curve is nowhere near on S, and a cone within
        // the snap tolerance on S is one whose one-ring the curve is inside.
        const Point seg = bestP - p;
        const double segLen = normP(seg);
        int snapSlot = -1, snapChild = -1;
        double snapGap = snapTol, snapS = 0.0;
        for (int q = 0; q < 3; ++q) {
            const int w = t[q];
            const int other = vertCone[w];
            if (other < 0) continue;
            double ss = 0.0;
            const double d = segmentClosest(uv[w], p, bestP, ss);
            // The cone it left from, before it has gone anywhere, is not a
            // termination. Once the curve is away it is a candidate like any
            // other -- a separatrix returning to its own singularity is Q5's
            // "possibly identical" case and the curve of the paper's Fig. 9.
            if (w == child && cv.imageLength + ss * segLen <= snapTol) continue;
            if (d < cv.nearestConeGap) { cv.nearestConeGap = d; cv.nearestCone = other; }
            if (d < snapGap) { snapGap = d; snapSlot = other; snapChild = w; snapS = ss; }
        }
        if (snapSlot >= 0) {
            const Point hit = p + seg * snapS;
            pushStep(entryB, baryOfPoint(f, hit));
            cv.end = End::Cone;
            cv.toCone = snapSlot;
            cv.toChild = snapChild;
            cv.gap = snapGap;
            cv.endDir = dir;
            return cv;
        }

        const int e = cm.triangleEdges[f][best];
        const int nb = cm.triangleAdjacency[f][best];
        const int x = t[best], y = t[(best + 1) % 3];

        pushStep(entryB, baryFromEdge(f, x, y, bestT));

        if (nb >= 0) {
            const Triangle &tn = cm.triangles[nb];
            const int lx = localOf(tn, x), ly = localOf(tn, y);
            entryEdge = (lx < 0 || ly < 0) ? -1 : localEdgeBetween(lx, ly);
            entryB = baryFromEdge(nb, x, y, bestT);
            f = nb;
            p = bestP;
            continue;
        }

        // Q5 case 2: out through dS, transversely.
        if (parentOnBoundary[e]) {
            cv.end = End::Boundary;
            cv.endDir = dir;
            return cv;
        }

        // Q5 case 3: across an arc of G. The point does not have to be carried
        // by the fitted transition -- the two child edges are the same edge of
        // S, so the same parameter along it is the same point, and using the
        // parameter rather than R_k x + t keeps the crossing exact however far
        // Stage 6 moved the map. Only the direction is turned, by the arc's k.
        auto it = seamSide.find(EdgeKey(x, y));
        if (it == seamSide.end()) { cv.end = End::Stuck; cv.endDir = dir; return cv; }
        const int pairIdx = it->second / 2;
        const bool onPlus = (it->second % 2) == 0;
        const auto &pr = pairs[pairIdx];
        const int k = arcList[pairArc[pairIdx]].k;

        int ta = -1, tb = -1, newDir = dir;
        double tPar = bestT;
        if (onPlus) {
            tPar = (x == pr[0]) ? bestT : 1.0 - bestT;
            ta = pr[2]; tb = pr[3];
            newDir = ((dir - k) % 4 + 4) % 4;
        } else {
            tPar = (x == pr[2]) ? bestT : 1.0 - bestT;
            ta = pr[0]; tb = pr[1];
            newDir = ((dir + k) % 4 + 4) % 4;
        }

        auto ie = cutEdgeIndex.find(EdgeKey(ta, tb));
        if (ie == cutEdgeIndex.end()) { cv.end = End::Stuck; cv.endDir = dir; return cv; }
        const int te = ie->second;
        const int tf = cm.edgeTriangles[te][0] >= 0 ? cm.edgeTriangles[te][0]
                                                    : cm.edgeTriangles[te][1];
        if (tf < 0) { cv.end = End::Stuck; cv.endDir = dir; return cv; }

        ++cv.seamCrossings;
        const Triangle &tt = cm.triangles[tf];
        const int la = localOf(tt, ta), lb = localOf(tt, tb);
        entryEdge = (la < 0 || lb < 0) ? -1 : localEdgeBetween(la, lb);
        entryB = baryFromEdge(tf, ta, tb, tPar);
        f = tf;
        dir = newDir;
        p = uv[ta] * (1.0 - tPar) + uv[tb] * tPar;
    }

    cv.end = End::Capped;
    cv.endDir = dir;
    return cv;
}

// ---------------------------------------------------------------------------
int Separatrices::faceContaining(int w, int dirT, int exclude) const {
    const Mesh &cm = imm->getCutMesh();
    if (w < 0 || w >= static_cast<int>(cm.vertices.size())) return -1;
    const Point e = axis(dirT);
    const auto &vt = cm.vertexTriangles;
    for (int i = vt.rowPtr[w]; i < vt.rowPtr[w + 1]; ++i) {
        const int g = vt.colIdx[i];
        if (g == exclude) continue;
        const int lw = localOf(cm.triangles[g], w);
        if (lw < 0) continue;
        const Point A = normalizeP(uv[cm.triangles[g][(lw + 1) % 3]] - uv[w]);
        const Point B = normalizeP(uv[cm.triangles[g][(lw + 2) % 3]] - uv[w]);
        // Half-open, closed at A and open at B, so a direction lying exactly
        // along a shared edge belongs to one of the two faces and not to both.
        if (cross2(A, e) > -1e-9 && cross2(e, B) > 1e-9) return g;
    }
    return -1;
}

// ---------------------------------------------------------------------------
std::array<double, 3> Separatrices::baryUnit(int face, int v) const {
    std::array<double, 3> b{{0.0, 0.0, 0.0}};
    const int l = localOf(imm->getCutMesh().triangles[face], v);
    if (l >= 0) b[l] = 1.0;
    return b;
}

std::array<double, 3> Separatrices::baryFromEdge(int face, int x, int y, double t) const {
    std::array<double, 3> b{{0.0, 0.0, 0.0}};
    const Triangle &tri = imm->getCutMesh().triangles[face];
    const int lx = localOf(tri, x), ly = localOf(tri, y);
    if (lx < 0 || ly < 0) return b;
    b[lx] = 1.0 - t;
    b[ly] = t;
    return b;
}

std::array<double, 3> Separatrices::baryOfPoint(int face, const Point &p) const {
    const Triangle &tri = imm->getCutMesh().triangles[face];
    const Point &a = uv[tri[0]], &b = uv[tri[1]], &c = uv[tri[2]];
    const Point v0 = b - a, v1 = c - a, v2 = p - a;
    const double den = cross2(v0, v1);
    std::array<double, 3> out{{1.0, 0.0, 0.0}};
    if (std::fabs(den) < 1e-300) return out;
    const double l1 = cross2(v2, v1) / den;
    const double l2 = cross2(v0, v2) / den;
    out[1] = std::max(0.0, std::min(1.0, l1));
    out[2] = std::max(0.0, std::min(1.0, l2));
    out[0] = std::max(0.0, 1.0 - out[1] - out[2]);
    const double s = out[0] + out[1] + out[2];
    if (s > 0.0) { out[0] /= s; out[1] /= s; out[2] /= s; }
    return out;
}

// ---------------------------------------------------------------------------
Point Separatrices::point(const Step &s, bool exit, Space space) const {
    const Mesh &cm = imm->getCutMesh();
    const Triangle &t = cm.triangles[s.face];
    const std::array<double, 3> &b = exit ? s.exit : s.entry;
    if (space == Space::Image) {
        return uv[t[0]] * b[0] + uv[t[1]] * b[1] + uv[t[2]] * b[2];
    }
    // Omega carries the positions of S, so the pullback is the same
    // barycentric combination of the same corners.
    return cm.vertices[t[0]] * b[0] + cm.vertices[t[1]] * b[1] + cm.vertices[t[2]] * b[2];
}

std::vector<Point> Separatrices::polyline(const Curve &c, Space space,
                                          std::vector<int> *breaks) const {
    std::vector<Point> out;
    if (breaks) breaks->clear();
    if (c.steps.empty()) return out;

    const double eps = 1e-9 * (space == Space::Image ? extent : modelExtent);
    out.push_back(point(c.steps.front(), false, space));
    for (size_t i = 0; i < c.steps.size(); ++i) {
        out.push_back(point(c.steps[i], true, space));
        if (i + 1 >= c.steps.size()) break;
        const Point next = point(c.steps[i + 1], false, space);
        // In the image the curve jumps to the far side of the cut at every seam
        // crossing; on the model the two sides are the same point of S and the
        // polyline runs on unbroken.
        if (normP(next - out.back()) > eps) {
            if (breaks) breaks->push_back(static_cast<int>(out.size()));
            out.push_back(next);
        }
    }
    return out;
}

bool Separatrices::writeOBJ(const std::string &filename, Space space) const {
    std::ofstream out(filename);
    if (!out) return false;

    int base = 1;
    std::vector<int> breaks;
    for (const Curve &c : traced) {
        const std::vector<Point> pts = polyline(c, space, &breaks);
        if (pts.size() < 2) continue;
        for (const Point &p : pts) out << "v " << p[0] << " " << p[1] << " 0\n";

        size_t bi = 0;
        out << "l";
        for (size_t i = 0; i < pts.size(); ++i) {
            if (bi < breaks.size() && static_cast<int>(i) == breaks[bi]) {
                out << "\nl";
                ++bi;
            }
            out << " " << (base + static_cast<int>(i));
        }
        out << "\n";
        base += static_cast<int>(pts.size());
    }
    return true;
}

// ---------------------------------------------------------------------------
// check()
// ---------------------------------------------------------------------------
void Separatrices::check() {
    report.worstMissGap = 0.0;
    report.worstMissCurve = -1;

    // Continuity of the pullback. Consecutive triangle crossings meet at a
    // point of S: within a sheet trivially, and across an arc of G because the
    // two child edges are the same edge of S at the same parameter. Measuring
    // it is how the seam handling is checked.
    double worstJoin = 0.0;
    for (const Curve &c : traced) {
        for (size_t i = 1; i < c.steps.size(); ++i) {
            const Point a = point(c.steps[i - 1], true, Space::Model);
            const Point b = point(c.steps[i], false, Space::Model);
            worstJoin = std::max(worstJoin, normP(b - a));
        }
    }
    report.maxPullbackGap = worstJoin / modelExtent;

    for (size_t i = 0; i < traced.size(); ++i) {
        const Curve &c = traced[i];
        if (c.end == End::Cone || c.end == End::Boundary) continue;
        const double g = std::isfinite(c.nearestConeGap) ? c.nearestConeGap : extent;
        if (report.worstMissCurve < 0 || g > report.worstMissGap) {
            report.worstMissGap = g;
            report.worstMissCurve = static_cast<int>(i);
        }
    }

    report.valid = report.fanFailures == 0 && report.capped == 0 && report.stuck == 0 &&
                   report.degenerate == 0 && report.emitted == report.prescribed &&
                   report.maxPullbackGap < 1e-9;

    if (report.maxPullbackGap >= 1e-9) {
        std::ostringstream oss;
        oss << "A traced curve is discontinuous on S by " << report.maxPullbackGap
            << " of the model between one triangle and the next, which means a seam was "
            << "crossed with the wrong pairing or the wrong parameter.";
        report.messages.push_back(oss.str());
    }

    if (report.capped > 0) {
        std::ostringstream oss;
        oss << report.capped << " separatrix/ces reached the step cap of " << options.maxSteps
            << " without terminating at a cone or leaving through dS: Q5 does not hold on "
            << "this map. Sec. 3.3's remedy is upstream -- raise lambda_5, add a Gamma_topo "
            << "constraint for the pair that nearly met, and re-run Stage 6 from the current "
            << "phi rather than from psi_R.";
        report.messages.push_back(oss.str());
    }
    if (report.stuck > 0) {
        std::ostringstream oss;
        oss << report.stuck << " separatrix/ces could not be continued: the walk met a "
            << "junction of the cutting graph head-on, where no single quarter turn carries "
            << "it across.";
        report.messages.push_back(oss.str());
    }
    if (report.degenerate > 0) {
        std::ostringstream oss;
        oss << report.degenerate << " ray/s never entered their starting triangle; the fan "
            << "sector they were placed in does not contain them.";
        report.messages.push_back(oss.str());
    }
    if (report.emitted != report.prescribed) {
        std::ostringstream oss;
        oss << "The cones emitted " << report.emitted << " separatrices where their indices "
            << "prescribe " << report.prescribed << " (4 - I in the interior, 1 - I on the "
            << "boundary).";
        report.messages.push_back(oss.str());
    }
}
