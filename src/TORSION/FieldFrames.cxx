#include "FieldFrames.hxx"

#include <algorithm>
#include <cmath>
#include <complex>
#include <deque>
#include <sstream>
#include <stdexcept>

#include "SIPG/SIPG.hxx"

namespace {

// The representative of d modulo pi/2, in [-pi/4, pi/4). Identical in
// intent to ConeSingularities' crossTurn, and written on the angles rather than
// on the representations because the angles are what the comb carries.
inline double wrapQuarter(double d) {
    const double quarter = M_PI_2;      // the period: the cross has four arms
    const double half = M_PI_4;         // half of it
    double r = std::fmod(d + half, quarter);
    if (r < 0.0) r += quarter;
    return r - half;
}

} // namespace

FieldFrames::FieldFrames(const SIPG &field, const ConeCut &cut,
                         const ConeSingularities &cones, const Options &options)
    : mesh(&cut.getOriginalMesh()), opts(options) {
    if (&cones.getMesh() != mesh) {
        throw std::runtime_error("FieldFrames: the cones were measured on a different mesh "
                                 "from the one the cut was built on");
    }
    if (&field.getMesh() != mesh) {
        throw std::runtime_error("FieldFrames: the field was solved on a different mesh from "
                                 "the one the cut was built on");
    }
    if (field.u_k.size() != static_cast<Eigen::Index>(mesh->triangles.size())) {
        throw std::runtime_error("FieldFrames: the field has one value per face and this one "
                                 "does not");
    }
    if (!opts.sizing.empty() && opts.sizing.size() != mesh->triangles.size()) {
        throw std::runtime_error("FieldFrames: the sizing field has one h per face and this "
                                 "one does not");
    }

    readField(field);
    comb(cut);
    buildFrames();
    audit(cut, cones);

    report.valid = report.unreachedFaces == 0 && report.combingDefects == 0 &&
                   report.leftHandedFrames == 0 && report.indexMismatches == 0;
}

// ---------------------------------------------------------------------------
// readField()
//
// theta_f off u_k, and the per-edge decomposition of the difference into a
// quarter turn and a remainder. Both are properties of the field alone -- no
// cut, no branch -- which is why they are taken first and why the comb below
// has nothing left to decide except which faces to visit.
// ---------------------------------------------------------------------------
void FieldFrames::readField(const SIPG &field) {
    const int nT = static_cast<int>(mesh->triangles.size());
    const int nE = static_cast<int>(mesh->edges.size());
    report.faces = nT;

    rawTheta.assign(nT, 0.0);
    for (int t = 0; t < nT; ++t) rawTheta[t] = std::arg(field.u_k[t]) * 0.25;

    matching.assign(nE, 0);
    delta.assign(nE, 0.0);
    for (int e = 0; e < nE; ++e) {
        const int f = mesh->edgeTriangles[e][0];
        const int g = mesh->edgeTriangles[e][1];
        if (f < 0 || g < 0) continue;
        const double d = rawTheta[g] - rawTheta[f];
        delta[e] = wrapQuarter(d);
        matching[e] = static_cast<int>(std::lround((d - delta[e]) / M_PI_2));
    }
}

// ---------------------------------------------------------------------------
// comb()
//
// BFS over the faces of Omega. Two faces of S are adjacent in Omega exactly
// when the edge between them is interior and is not an edge of G, and that test
// is the whole of what makes this traversal path-independent -- so it is made
// against ConeCut's own edge set rather than against anything rederived here.
//
// The loop check is the second half of the same statement and costs one more
// pass: every dual edge of Omega that the tree did not use closes a cycle, and
// on a disk with no cone inside it the branch has to come back to itself.
// ---------------------------------------------------------------------------
void FieldFrames::comb(const ConeCut &cut) {
    const int nT = static_cast<int>(mesh->triangles.size());
    const auto &cutEdges = cut.getCutEdges();

    auto crossesG = [&](int e) {
        return cutEdges.count(MeshEdgeKey(mesh->edges[e][0], mesh->edges[e][1])) > 0;
    };

    period.assign(nT, 0);
    std::vector<char> seen(nT, 0);

    int seed = opts.seedFace;
    if (seed < 0 || seed >= nT) seed = 0;
    report.seedFace = seed;

    std::deque<int> q;
    seen[seed] = 1;
    period[seed] = 0;
    q.push_back(seed);
    int reached = 1;

    while (!q.empty()) {
        const int f = q.front();
        q.pop_front();
        for (int le = 0; le < 3; ++le) {
            const int e = mesh->triangleEdges[f][le];
            const int a = mesh->edgeTriangles[e][0];
            const int b = mesh->edgeTriangles[e][1];
            if (a < 0 || b < 0) continue;
            const int g = (a == f) ? b : a;
            if (crossesG(e)) continue;
            if (seen[g]) continue;
            // matching[] is oriented from edgeTriangles[e][0] to [1]; the walk
            // may be going the other way, and p_gf = -p_fg.
            const int p = (a == f) ? matching[e] : -matching[e];
            period[g] = period[f] - p;
            seen[g] = 1;
            ++reached;
            q.push_back(g);
        }
    }

    report.combedFaces = reached;
    report.unreachedFaces = nT - reached;

    theta.assign(nT, 0.0);
    for (int t = 0; t < nT; ++t) theta[t] = rawTheta[t] + M_PI_2 * period[t];

    // The loops. Every interior, non-G edge is a dual edge of Omega; the tree
    // used |Omega faces| - 1 of them and the rest close cycles.
    for (int e = 0; e < static_cast<int>(mesh->edges.size()); ++e) {
        const int f = mesh->edgeTriangles[e][0];
        const int g = mesh->edgeTriangles[e][1];
        if (f < 0 || g < 0) continue;
        if (crossesG(e)) continue;
        if (!seen[f] || !seen[g]) continue;
        if (period[g] - period[f] != -matching[e]) ++report.combingDefects;
        report.maxFrameJump = std::max(report.maxFrameJump, std::fabs(delta[e]));
    }

    if (report.unreachedFaces > 0) {
        std::ostringstream oss;
        oss << report.unreachedFaces << " face(s) of Omega were never reached from the seed: "
            << "the dual graph of the cut mesh is not connected, so there is no single branch "
            << "of the field over it and no map to integrate.";
        report.messages.push_back(oss.str());
    }
    if (report.combingDefects > 0) {
        std::ostringstream oss;
        oss << report.combingDefects << " dual edge(s) of Omega close a loop the branch does "
            << "not: a_g - a_f disagrees with -p_fg. Omega is a disk with every cone on its "
            << "boundary, so the traversal leaked across G somewhere -- the arc membership "
            << "test is what to look at, not the field.";
        report.messages.push_back(oss.str());
    }
}

// ---------------------------------------------------------------------------
// buildFrames()
//
// J*_t and the field metric's edge lengths. The metric is (1/h_t^2) I because
// the frame is a scaled rotation, so an edge's length under it is its Euclidean
// length over h -- one number from each of its two triangles, averaged, with
// the disagreement between them reported. Sec. 6.3 asks for that number; what
// it measures here is the *sizing* field's variation across the edge, since a
// uniform h makes the two answers identical. The non-integrability that section
// is really after is Report::maxFrameJump above, and the per-triangle fit
// residual FieldIntegration reports after the solve.
// ---------------------------------------------------------------------------
void FieldFrames::buildFrames() {
    const int nT = static_cast<int>(mesh->triangles.size());
    h.assign(nT, opts.targetEdge);
    if (!opts.sizing.empty()) {
        for (int t = 0; t < nT; ++t) {
            if (opts.sizing[t] > 0.0 && std::isfinite(opts.sizing[t])) h[t] = opts.sizing[t];
            else ++report.sizingFallbacks;
        }
    }
    for (int t = 0; t < nT; ++t) {
        if (!(h[t] > 0.0) || !std::isfinite(h[t])) { h[t] = 1.0; ++report.sizingFallbacks; }
    }

    frame.assign(nT, {0.0, 0.0, 0.0, 0.0});
    for (int t = 0; t < nT; ++t) {
        const double c = std::cos(theta[t]) / h[t];
        const double s = std::sin(theta[t]) / h[t];
        // Rows are X_t^T and Y_t^T: the gradients of u and v the map is to have.
        frame[t] = {c, s, -s, c};
        const double det = frame[t][0] * frame[t][3] - frame[t][1] * frame[t][2];
        if (!(det > 0.0) || !std::isfinite(det)) ++report.leftHandedFrames;
    }
    if (report.leftHandedFrames > 0) {
        std::ostringstream oss;
        oss << report.leftHandedFrames << " triangle(s) carry a frame with det J* <= 0. "
            << "det J* is 1/h^2 for any real combed angle, so this is a NaN in the field or a "
            << "zero in the sizing field; until it is fixed the flip cap's identity "
            << "\"det J > 0 iff the image triangle is not inverted\" is void on them.";
        report.messages.push_back(oss.str());
    }

    const int nE = static_cast<int>(mesh->edges.size());
    metricLen.assign(nE, 0.0);
    for (int e = 0; e < nE; ++e) {
        const double len = normP(mesh->vertices[mesh->edges[e][1]] -
                                 mesh->vertices[mesh->edges[e][0]]);
        const int f = mesh->edgeTriangles[e][0];
        const int g = mesh->edgeTriangles[e][1];
        if (f >= 0 && g >= 0) {
            const double lf = len / h[f], lg = len / h[g];
            metricLen[e] = 0.5 * (lf + lg);
            if (metricLen[e] > 0.0) {
                report.maxMetricDisagreement =
                    std::max(report.maxMetricDisagreement,
                             std::fabs(lf - lg) / metricLen[e]);
            }
        } else {
            const int only = (f >= 0) ? f : g;
            metricLen[e] = (only >= 0) ? len / h[only] : len;
        }
    }
}

// ---------------------------------------------------------------------------
// audit()  --  Sec. 5.1
// ---------------------------------------------------------------------------
void FieldFrames::audit(const ConeCut &cut, const ConeSingularities &cones) {
    const int nV = static_cast<int>(mesh->vertices.size());
    const auto &vt = mesh->vertexTriangles;
    const std::vector<int> &current = cones.getIndices();
    const std::vector<int> &want =
        opts.referenceIndex.size() == static_cast<size_t>(nV) ? opts.referenceIndex : current;

    // The shared edge of two faces of one ring, for reading p between them.
    auto matchingBetween = [&](int f, int g) -> int {
        for (int le = 0; le < 3; ++le) {
            const int e = mesh->triangleEdges[f][le];
            const int a = mesh->edgeTriangles[e][0];
            const int b = mesh->edgeTriangles[e][1];
            if ((a == f && b == g)) return matching[e];
            if ((b == f && a == g)) return -matching[e];
        }
        return 0;   // not adjacent; the ring is not a fan of shared edges
    };

    indexFromMatchings.assign(nV, 0);
    double worst = -1.0;
    for (int v = 0; v < nV; ++v) {
        if (mesh->isBoundaryVertex[v]) {
            // Not a closed loop, so the ring sum is not the index. The value
            // ConeSingularities read from the boundary fan
            // (I ~ (Theta + pi - Omega)/(pi/2), rounded per boundary loop) is
            // the one that applies and is carried through unaudited.
            indexFromMatchings[v] = current[v];
            continue;
        }
        const int begin = vt.rowPtr[v], end = vt.rowPtr[v + 1];
        const int star = end - begin;
        if (star < 3) { indexFromMatchings[v] = current[v]; continue; }

        int sum = 0;
        for (int k = begin; k < end; ++k) {
            const int tCur = vt.colIdx[k];
            const int tNext = vt.colIdx[(k - begin + 1) % star + begin];
            sum += matchingBetween(tCur, tNext);
        }
        // sum(delta) + (pi/2) sum(p) = 0 around a closed ring, and the field's
        // index is sum(delta) / (pi/2).
        indexFromMatchings[v] = -sum;

        if (indexFromMatchings[v] != want[v]) {
            ++report.indexMismatches;
            const double d = std::fabs(static_cast<double>(indexFromMatchings[v] - want[v]));
            if (d > worst) { worst = d; report.worstMismatchVertex = v; }
        }
    }

    // --- the boundary half of the audit -------------------------------
    //
    // At a vertex of dS with boundary edges (a, v) and (v, b) walked with the
    // interior on the left, the image tangents are J*_{f1}(v - a) and
    // J*_{f2}(b - v) -- both axis aligned, because the field is tangent to dS
    // and the frame is that tangent rotated onto an axis. The turn between them
    // is the quarter turn the layout has to make at v, and pi - (pi/2) I is the
    // angle Q2 then asks of it, so the implied index is that number of quarters.
    frameBoundaryIndex.assign(nV, 0);
    {
        // The boundary edges at each vertex, and the face each belongs to.
        std::vector<std::array<int, 2>> incident(nV, {-1, -1});
        std::vector<int> degree(nV, 0);
        for (int e : mesh->boundaryEdges) {
            for (int i = 0; i < 2; ++i) {
                const int v = mesh->edges[e][i];
                ++degree[v];
                if (incident[v][0] < 0) incident[v][0] = e;
                else if (incident[v][1] < 0) incident[v][1] = e;
            }
        }
        auto imageDir = [&](int e, int from, int to) -> double {
            const int f = (mesh->edgeTriangles[e][0] >= 0) ? mesh->edgeTriangles[e][0]
                                                           : mesh->edgeTriangles[e][1];
            if (f < 0) return 0.0;
            const Point d = mesh->vertices[to] - mesh->vertices[from];
            const std::array<double, 4> &J = frame[f];
            return std::atan2(J[2] * d[0] + J[3] * d[1], J[0] * d[0] + J[1] * d[1]);
        };
        // Which way round the two boundary edges at v run is settled by the
        // face they sit in: the one that traverses (a, v) in its own winding is
        // the incoming edge, since the mesh is positively oriented and dS is
        // therefore walked with the interior on the left.
        auto traversedAs = [&](int e, int from, int to) -> bool {
            const int f = (mesh->edgeTriangles[e][0] >= 0) ? mesh->edgeTriangles[e][0]
                                                           : mesh->edgeTriangles[e][1];
            if (f < 0) return false;
            const Triangle &t = mesh->triangles[f];
            for (int i = 0; i < 3; ++i) {
                if (t[i] == from && t[(i + 1) % 3] == to) return true;
            }
            return false;
        };

        for (int v = 0; v < nV; ++v) {
            if (!mesh->isBoundaryVertex[v]) continue;
            // A pinch -- more than two boundary edges at one vertex -- has no
            // single corner to read, and reading two of the several would be
            // reporting a number about a vertex that does not have one.
            if (degree[v] != 2) continue;
            const int e1 = incident[v][0], e2 = incident[v][1];
            if (e1 < 0 || e2 < 0) continue;
            const int a = (mesh->edges[e1][0] == v) ? mesh->edges[e1][1] : mesh->edges[e1][0];
            const int b = (mesh->edges[e2][0] == v) ? mesh->edges[e2][1] : mesh->edges[e2][0];

            int inE = e1, outE = e2, inFrom = a, outTo = b;
            if (!traversedAs(e1, a, v)) {
                if (traversedAs(e2, b, v)) { inE = e2; outE = e1; inFrom = b; outTo = a; }
                else continue;   // neither runs into v; the winding is not readable here
            }
            const double before = imageDir(inE, inFrom, v);
            const double after = imageDir(outE, v, outTo);
            const double turn = wrap_pi(after - before) / M_PI_2;
            const long q = std::lround(turn);
            const double residual = std::fabs(turn - static_cast<double>(q)) * M_PI_2;
            if (residual > report.maxBoundaryTurnResidual) {
                report.maxBoundaryTurnResidual = residual;
            }
            frameBoundaryIndex[v] = static_cast<int>(q);
            if (q != 0) ++report.boundaryTurns;
            if (static_cast<int>(q) != current[v]) {
                ++report.boundaryTurnMismatches;
                if (report.worstBoundaryTurnVertex < 0) report.worstBoundaryTurnVertex = v;
            }
        }
    }

    // Eq. (4) on the set as it stands, which is the set Stage 2 cut to.
    const ConeSingularities::GaussBonnetReport gb = cones.gaussBonnet();
    report.indexSum = gb.indexSum;
    report.indexTarget = gb.indexTarget;
    report.admissible = gb.admissible;

    // Two cones in one triangle. Immersion measures the same thing a stage
    // later as a distance in the image; here it is exact and combinatorial, and
    // finding it now is what Sec. 5.1 asks for.
    for (const Triangle &t : mesh->triangles) {
        int carried = 0;
        for (int i = 0; i < 3; ++i) if (current[t[i]] != 0) ++carried;
        if (carried >= 2) ++report.clusteredConePairs;
    }
    for (const auto &c : cones.getCones()) {
        if (std::abs(c.index) > 2) ++report.highIndexCones;
    }

    if (report.boundaryTurnMismatches > 0) {
        std::ostringstream oss;
        oss << report.boundaryTurnMismatches << " vertex/vertices of dS where the frame turns "
            << "by a different number of quarters from the index the cone set carries (first at "
            << "vertex " << report.worstBoundaryTurnVertex << "; " << report.boundaryTurns
            << " vertices turn at all, worst rounding " << report.maxBoundaryTurnResidual
            << " rad from a whole quarter). Each of these is a corner the layout's boundary "
            << "will have and Q2 will report as a non-cone vertex off pi by a quarter turn, and "
            << "each becomes a node of the arrangement that Pipeline A does not have. It is not "
            << "a defect of the combing: a field tangent to a curving boundary quantises it into "
            << "a staircase, and only a stage that makes the boundary geodesic -- which is what "
            << "Ricci flow does and what this route does not have -- removes the steps.";
        report.messages.push_back(oss.str());
    }
    if (report.indexMismatches > 0) {
        std::ostringstream oss;
        oss << report.indexMismatches << " interior vertex/vertices where the index recomputed "
            << "from the matchings disagrees with the one the field was read for (worst at "
            << "vertex " << report.worstMismatchVertex << ": "
            << indexFromMatchings[report.worstMismatchVertex] << " against "
            << want[report.worstMismatchVertex] << "). Either the ring sum and the field's own "
            << "measurement have parted company, or the cone set was moved after it was read "
            << "-- prescribe(), rebalance() and cancelDipoles() all do that on purpose, and "
            << "Options::referenceIndex is how to tell the audit which set it is auditing.";
        report.messages.push_back(oss.str());
    }
    if (!report.admissible) {
        std::ostringstream oss;
        oss << "sum I(v) = " << report.indexSum << " against the 4 chi(S) = "
            << report.indexTarget << " of Eq. (4): the flat cone metric this set asks for does "
            << "not exist, and neither does a seamless map with these holonomies.";
        report.messages.push_back(oss.str());
    }
    if (report.clusteredConePairs > 0) {
        std::ostringstream oss;
        oss << report.clusteredConePairs << " triangle(s) carry two or more cones. Sec. 3.1 "
            << "calls this clustering and merges such a cluster into one cone of the summed "
            << "index; left alone it comes back as a Q2 residual out of Stage 6.";
        report.messages.push_back(oss.str());
    }
    if (report.highIndexCones > 0) {
        std::ostringstream oss;
        oss << report.highIndexCones << " cone(s) of index |I| > 2. Legitimate -- the paper "
            << "wants high-valence cones -- and rare by accident, so it is reported and not "
            << "rejected.";
        report.messages.push_back(oss.str());
    }
    if (cut.getReport().interiorJunctions > 0) {
        std::ostringstream oss;
        oss << cut.getReport().interiorJunctions << " interior junction(s) in G. The seam "
            << "elimination of Sec. 6.1 assumes a forest; the constrained solve falls back to "
            << "the saddle system, which does not, so this costs conditioning and not "
            << "correctness.";
        report.messages.push_back(oss.str());
    }
}
