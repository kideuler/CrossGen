#include "FieldFrames.hxx"

#include <algorithm>
#include <cmath>
#include <complex>
#include <deque>
#include <unordered_map>
#include <sstream>
#include <stdexcept>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

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
                         const ConeSingularities &cones, const Options &options,
                         const Interfaces *interfaces)
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
    if (opts.reconcileIndices) reconcile(cut, cones);
    comb(cut);
    if (report.reconciledEdges > 0) smoothBranch(cut);
    buildFrames();
    audit(cut, cones);
    if (opts.buildAlignment) buildAlignment(cut, cones, interfaces);

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
// ringSign()  --  which way an edge enters a vertex's ring sum
//
// The index at an interior vertex is minus the sum of the matchings taken
// between consecutive faces of its one ring, and each edge at the vertex
// contributes exactly one of those terms -- with a sign, because matching[] is
// oriented from edgeTriangles[e][0] to [1] and the ring may cross it either
// way. The two ends of an edge get opposite signs, which is what makes adding a
// constant to one matching a *transport* of index from one end to the other
// rather than a change to both.
//
// Read off Mesh::vertexTriangles rather than derived from the winding, so that
// it is the same traversal audit() sums over and cannot disagree with it.
// ---------------------------------------------------------------------------
int FieldFrames::ringSign(int v, int e) const {
    const auto &vt = mesh->vertexTriangles;
    const int begin = vt.rowPtr[v], end = vt.rowPtr[v + 1];
    const int star = end - begin;
    if (star < 3) return 0;
    const int a = mesh->edgeTriangles[e][0], b = mesh->edgeTriangles[e][1];
    if (a < 0 || b < 0) return 0;
    for (int k = begin; k < end; ++k) {
        const int cur = vt.colIdx[k];
        const int next = vt.colIdx[(k - begin + 1) % star + begin];
        if (cur == a && next == b) return +1;
        if (cur == b && next == a) return -1;
    }
    return 0;
}

// ---------------------------------------------------------------------------
// indexRing()
//
// I(v) = - sum over the one ring of p, at every interior vertex. Boundary
// vertices are left at zero: their fan is not a closed loop and their index is
// ConeSingularities::computeBoundaryIndices' business.
// ---------------------------------------------------------------------------
void FieldFrames::indexRing(std::vector<int> &out) const {
    const int nV = static_cast<int>(mesh->vertices.size());
    const auto &vt = mesh->vertexTriangles;
    out.assign(nV, 0);
    for (int v = 0; v < nV; ++v) {
        if (mesh->isBoundaryVertex[v]) continue;
        const int begin = vt.rowPtr[v], end = vt.rowPtr[v + 1];
        const int star = end - begin;
        if (star < 3) continue;
        int sum = 0;
        for (int k = begin; k < end; ++k) {
            const int cur = vt.colIdx[k];
            const int next = vt.colIdx[(k - begin + 1) % star + begin];
            for (int le = 0; le < 3; ++le) {
                const int e = mesh->triangleEdges[cur][le];
                const int a = mesh->edgeTriangles[e][0], b = mesh->edgeTriangles[e][1];
                if (a == cur && b == next) { sum += matching[e]; break; }
                if (b == cur && a == next) { sum -= matching[e]; break; }
            }
        }
        out[v] = -sum;
    }
}

// ---------------------------------------------------------------------------
// reconcile()  --  the field's singularities against Stage 1's cone set
//
// Stage 1 does not hand the cut the indices the field read. prescribe() writes
// the interface network's own indices over them, cancelDipoles() removes the
// +1/-1 pairs an aligned field carries on a curved interface, and rebalance()
// moves whatever is left to restore Eq. (4). All three are right, and all three
// are stated about the *layout*.
//
// Pipeline A can take them at face value because the flow is driven by the cone
// set and the field is never looked at again. This pipeline integrates the
// field, and there the difference is not a matter of taste:
//
//     Stage 2 cuts to the cones Stage 1 ended with. A vertex the field is
//     singular at and Stage 1 called regular is therefore *not* cut to, a loop
//     of Omega encircles it, and the branch the comb selects is not single
//     valued -- which is not a degradation but the end of the route. On
//     multimat/geom011, where cancelDipoles() removes three pairs, 208 dual
//     edges of Omega close a loop the branch does not and the pipeline stops
//     before the integration.
//
// So the field is moved to the cone set rather than the cone set to the field.
// The move is exact and combinatorial. I(v) is minus the ring sum of the
// matchings, and each edge enters the ring sums of its two endpoints with
// opposite signs, so adding one to a matching moves one unit of index from one
// end of that edge to the other. A shortest path in Omega from a vertex with
// too much to one with too little -- or to dS, whose index is not a ring sum
// and which absorbs anything -- transports it the whole way, changing nothing
// in between.
//
// The path is not allowed to cross G, because the quarter turn of each arc is
// read off the frames either side of it and moving a matching there would move
// Q4 with it.
//
// What the transport costs is smoothness: across each edge it touched, the
// combed angle now jumps by a quarter that the field does not. smoothBranch()
// is what pays that back.
// ---------------------------------------------------------------------------
void FieldFrames::reconcile(const ConeCut &cut, const ConeSingularities &cones) {
    const int nV = static_cast<int>(mesh->vertices.size());
    const std::vector<int> &want = cones.getIndices();
    if (static_cast<int>(want.size()) != nV) return;

    std::vector<int> have;
    indexRing(have);

    std::vector<int> excess(nV, 0);
    int outstanding = 0;
    for (int v = 0; v < nV; ++v) {
        if (mesh->isBoundaryVertex[v]) continue;
        excess[v] = have[v] - want[v];
        if (excess[v] != 0) { ++report.reconciledVertices; outstanding += std::abs(excess[v]); }
    }
    if (outstanding == 0) return;

    const auto &cutEdges = cut.getCutEdges();
    auto crossesG = [&](int e) {
        return cutEdges.count(MeshEdgeKey(mesh->edges[e][0], mesh->edges[e][1])) > 0;
    };

    // The vertices each vertex can reach without crossing G. Built once as an
    // edge list per vertex, since the search below runs a handful of times.
    std::vector<std::vector<int>> incident(nV);
    for (int e = 0; e < static_cast<int>(mesh->edges.size()); ++e) {
        if (crossesG(e)) continue;
        incident[mesh->edges[e][0]].push_back(e);
        incident[mesh->edges[e][1]].push_back(e);
    }

    std::vector<int> cameBy(nV, -1), cameFrom(nV, -1);
    std::deque<int> q;
    int guard = 0;
    const int guardLimit = 4 * outstanding + 16;

    while (outstanding > 0 && guard++ < guardLimit) {
        int source = -1;
        for (int v = 0; v < nV && source < 0; ++v) if (excess[v] != 0) source = v;
        if (source < 0) break;
        const int sign = (excess[source] > 0) ? +1 : -1;

        // Nearest sink: a vertex wanting the opposite, or dS, whose index is
        // read from its boundary fan and not from a ring sum.
        std::fill(cameBy.begin(), cameBy.end(), -1);
        std::fill(cameFrom.begin(), cameFrom.end(), -1);
        q.clear();
        q.push_back(source);
        cameFrom[source] = source;
        int sink = -1;
        while (!q.empty() && sink < 0) {
            const int v = q.front();
            q.pop_front();
            for (int e : incident[v]) {
                const int w = (mesh->edges[e][0] == v) ? mesh->edges[e][1] : mesh->edges[e][0];
                if (cameFrom[w] >= 0) continue;
                cameFrom[w] = v;
                cameBy[w] = e;
                if (mesh->isBoundaryVertex[w] ||
                    (excess[w] != 0 && (excess[w] > 0) != (sign > 0))) { sink = w; break; }
                q.push_back(w);
            }
        }
        if (sink < 0) { ++report.reconcileFailures; excess[source] = 0; --outstanding; continue; }

        for (int w = sink; w != source; w = cameFrom[w]) {
            const int e = cameBy[w];
            const int from = cameFrom[w];
            const int s = ringSign(from, e);
            if (s == 0) continue;   // `from` is on dS: nothing to transport out of
            // index(from) changes by -s * delta, so delta = s * sign removes
            // one unit of `sign` from it and puts it at the other end.
            matching[e] += s * sign;
            ++report.reconciledEdges;
        }
        excess[source] -= sign;
        if (!mesh->isBoundaryVertex[sink]) excess[sink] += sign;
        --outstanding;
        ++report.reconciledUnits;
    }

    // delta is the field's own variation less the quarter turn the matching now
    // claims, and the matchings just moved, so it is re-read rather than left
    // saying something about the set they used to be.
    for (int e = 0; e < static_cast<int>(mesh->edges.size()); ++e) {
        const int f = mesh->edgeTriangles[e][0], g = mesh->edgeTriangles[e][1];
        if (f < 0 || g < 0) continue;
        delta[e] = (rawTheta[g] - rawTheta[f]) - M_PI_2 * matching[e];
    }

    std::ostringstream oss;
    oss << "The field is singular at " << report.reconciledVertices
        << " interior vertex/vertices where Stage 1's cone set is not, or carries a different "
        << "index: " << report.reconciledUnits << " index unit(s) were transported along "
        << report.reconciledEdges << " edge(s) of Omega so that the two agree. Without it the "
        << "cut does not enclose the field's singularities and the branch over Omega is not "
        << "single valued -- prescribe(), cancelDipoles() and rebalance() all move indices on "
        << "purpose, and Pipeline A can take that at face value because it never looks at the "
        << "field again.";
    report.messages.push_back(oss.str());
    if (report.reconcileFailures > 0) {
        std::ostringstream f;
        f << report.reconcileFailures << " index unit(s) had nowhere to go: no vertex wanting "
          << "the opposite and no boundary reachable without crossing G. The comb will report "
          << "the loops they leave open.";
        report.messages.push_back(f.str());
    }
}

// ---------------------------------------------------------------------------
// smoothBranch()
//
// reconcile() moved singularities by rewriting matchings, and a rewritten
// matching is a quarter turn the combed angle now takes across an edge that the
// field does not take. Left alone that is a crease in J* along the transport
// path, and the integration answers a crease with inverted triangles.
//
// The repair is the same statement the comb already is, solved rather than
// walked. theta_hat is the function on the faces of Omega whose difference
// across every dual edge is delta_fg, and the comb finds it by following a
// spanning tree; the edges reconcile() touched are the ones where that
// difference is now a quarter out. So ask for the function whose differences
// are the field's *original* variation everywhere and zero across the touched
// edges,
//
//     min_theta  sum_e w_e ( theta_g - theta_f - t_e )^2 ,
//     t_e = delta_e originally, 0 where a matching was moved
//
// which is one Poisson solve on the dual graph of Omega with theta pinned at
// the seed face. Where nothing was moved t_e is exactly the comb's own
// difference and the answer is the comb, to rounding; where something was, the
// quarter is spread over the neighbourhood instead of standing on one edge.
//
// The target is not integrable -- that is the whole point, the index has to go
// somewhere -- so this is a least-squares fit and not an interpolation, and
// Report::maxBranchCorrection is how far it had to move the branch to make it.
// ---------------------------------------------------------------------------
void FieldFrames::smoothBranch(const ConeCut &cut) {
    const int nT = static_cast<int>(mesh->triangles.size());
    if (nT == 0 || theta.size() != static_cast<size_t>(nT)) return;

    const auto &cutEdges = cut.getCutEdges();
    auto crossesG = [&](int e) {
        return cutEdges.count(MeshEdgeKey(mesh->edges[e][0], mesh->edges[e][1])) > 0;
    };

    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(static_cast<size_t>(mesh->edges.size()) * 4 + 1);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(nT);

    for (int e = 0; e < static_cast<int>(mesh->edges.size()); ++e) {
        const int f = mesh->edgeTriangles[e][0], g = mesh->edgeTriangles[e][1];
        if (f < 0 || g < 0) continue;
        if (crossesG(e)) continue;
        // The comb's own difference across this dual edge, and the one it
        // should have had: they differ only where reconcile() moved a matching,
        // and there by exactly the quarter it moved.
        const double now = theta[g] - theta[f];
        const double t = wrapQuarter(now);
        const double w = 1.0;
        trip.emplace_back(f, f, w);
        trip.emplace_back(g, g, w);
        trip.emplace_back(f, g, -w);
        trip.emplace_back(g, f, -w);
        rhs[f] -= w * t;
        rhs[g] += w * t;
    }

    int seed = report.seedFace;
    if (seed < 0 || seed >= nT) seed = 0;
    // The pin, as a large diagonal rather than as an eliminated row, so that the
    // pattern stays symmetric positive definite and the seed keeps the value the
    // comb gave it.
    double scale = 0.0;
    for (const auto &t : trip) if (t.row() == t.col()) scale = std::max(scale, t.value());
    const double pin = std::max(1.0, scale) * 1e6;
    trip.emplace_back(seed, seed, pin);
    rhs[seed] += pin * theta[seed];

    Eigen::SparseMatrix<double> L(nT, nT);
    L.setFromTriplets(trip.begin(), trip.end());
    L.makeCompressed();
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
    ldlt.compute(L);
    if (ldlt.info() != Eigen::Success) {
        report.messages.push_back(
            "The branch could not be re-smoothed after the index transport; J* keeps its "
            "quarter-turn crease along the transported path and the integration will invert "
            "there.");
        return;
    }
    const Eigen::VectorXd sol = ldlt.solve(rhs);
    if (ldlt.info() != Eigen::Success || !sol.allFinite()) {
        report.messages.push_back("The branch re-smoothing did not solve.");
        return;
    }
    for (int t = 0; t < nT; ++t) {
        report.maxBranchCorrection =
            std::max(report.maxBranchCorrection, std::fabs(sol[t] - theta[t]));
        theta[t] = sol[t];
    }
    report.branchResmoothed = true;

    std::ostringstream oss;
    oss << "The branch was re-smoothed after the transport: the quarter turn each moved "
        << "matching put on one dual edge is spread over its neighbourhood instead, and the "
        << "branch moved by at most " << report.maxBranchCorrection << " rad doing it.";
    report.messages.push_back(oss.str());
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
    if (static_cast<int>(opts.metricLengths.size()) == nE) {
        // Supplied. maxMetricDisagreement then reads how far it is from the
        // |e|/h the frame itself would have given, face by face -- which is the
        // difference between a vertex scaling and a per-face one, and is the
        // number that says whether the frame and the metric are telling the
        // integration the same thing.
        metricLen = opts.metricLengths;
        for (int e = 0; e < nE; ++e) {
            if (!(metricLen[e] > 0.0)) continue;
            const double len = normP(mesh->vertices[mesh->edges[e][1]] -
                                     mesh->vertices[mesh->edges[e][0]]);
            for (int i = 0; i < 2; ++i) {
                const int f = mesh->edgeTriangles[e][i];
                if (f < 0) continue;
                report.maxMetricDisagreement =
                    std::max(report.maxMetricDisagreement,
                             std::fabs(len / h[f] - metricLen[e]) / metricLen[e]);
            }
        }
        return;
    }
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


// ---------------------------------------------------------------------------
// buildAlignment()  --  Sec. 6.4
//
// Which coordinate each edge of Omega is to hold constant, for the integration
// to impose as an equation rather than for Stage 6 to ask for with a penalty.
//
// The statement is short. The field is pinned tangent to dS and to every
// interface, so J*_t applied to an edge of either is axis aligned already, and
// the axis is the whole number of quarter turns
//
//     q(e) = round( (dir(e) - theta_hat_f) / (pi/2) )
//
// with dir(e) taken in the direction the curve is walked. An even q puts the
// image of e on a horizontal, whose v is constant; an odd one on a vertical.
// One equation per edge, u_i = u_j or v_i = v_j, and Q3 holds the moment the
// solve returns instead of after a continuation has driven it there.
//
// ### Why the axis is decided per chain and not per edge
//
// A cross field following a curving boundary quantises it into a staircase:
// theta_hat tracks the tangent until the tracking slips a quarter, and it slips
// where the field's own smoothness puts the slip rather than where Stage 1 put
// a cone. Every step inside a chain is a corner in the image that Stage 8 finds
// as a node, that Q2 reports as a non-cone vertex off pi by a quarter turn, and
// that no amount of lambda_2 removes, because E2 is written per edge and a
// staircase satisfies it exactly.
//
// So the chains are cut at the points the layout is *allowed* to turn -- the
// cones, the nodes of the interface network, and wherever the cutting graph
// interrupts dS -- and each chain gets one axis, by the vote of its edges
// weighted by their length. That is the same thing Pipeline A gets from Stage
// 3: prescribing zero curvature at every non-cone boundary vertex is exactly
// "no step between two cones", and the flow straightens the boundary to match.
//
// ### Why the vote is on q and not on its parity
//
// The constraint can only express the parity -- u_i = u_j says the image is
// vertical and says nothing about which way up -- so the parity is all that is
// written down. It is not all that is *read*. Two edges of one chain with
// q = 1 and q = 3 have the same parity and opposite directions: holding both to
// u = const is satisfiable and folds the chain back along itself, and the
// interior folds with it. Voting on q mod 4 is what sees that, and a chain
// where it happens is left free rather than folded -- the field is saying the
// boundary turns half way round inside a run Stage 1 called straight, and the
// place to fix that is the cone set.
//
// Reading q mod 4 at all needs the edges *oriented*, since reversing an edge
// changes q by two; that is why the walk carries the direction it is going
// rather than reading Mesh::edges as stored.
//
// ### The two arrays
//
// `alignAxis` holds every chain that was not folded. `strictAxis` holds only
// the chains whose edges all agreed -- no staircase step to override -- and is
// what Stage 4F falls back to when holding the full set inverted something: a
// step is a right angle the map is being asked to unbend, and unbending it is
// where the least-squares fit runs out of room. Between them and the free solve
// the integration has three attempts, and takes the first that inverts nothing.
// ---------------------------------------------------------------------------
void FieldFrames::buildAlignment(const ConeCut &cut, const ConeSingularities &cones,
                                 const Interfaces *interfaces) {
    const Mesh &om = cut.getCutMesh();
    const auto &c2o = cut.getCutVertexToOriginal();
    const std::vector<int> &index = cones.getIndices();
    const Interfaces *itf = (interfaces && interfaces->multiMaterial()) ? interfaces : nullptr;

    const int nCutE = static_cast<int>(om.edges.size());
    alignAxis.assign(nCutE, -1);
    strictAxis.assign(nCutE, -1);
    if (nCutE == 0 || theta.empty()) return;

    // --- the parent of every edge of Omega --------------------------------
    std::unordered_map<MeshEdgeKey, int, MeshEdgeKeyHash> origEdge;
    origEdge.reserve(mesh->edges.size() * 2);
    for (int e = 0; e < static_cast<int>(mesh->edges.size()); ++e) {
        origEdge.emplace(MeshEdgeKey(mesh->edges[e][0], mesh->edges[e][1]), e);
    }
    std::vector<int> parent(nCutE, -1);
    for (int e = 0; e < nCutE; ++e) {
        const int a = om.edges[e][0], b = om.edges[e][1];
        if (a < 0 || b < 0) continue;
        auto it = origEdge.find(MeshEdgeKey(c2o[a], c2o[b]));
        if (it != origEdge.end()) parent[e] = it->second;
    }

    // One step of a chain: the edge of S it runs along, walked from `from` to
    // `to`, and the edges of Omega that are copies of it. An interface edge the
    // cutting graph ran along has two, and both sides want the same axis.
    struct Step {
        int orig = -1;
        int from = -1, to = -1;      // original vertex ids, in walk order
        std::vector<int> cutEdges;
    };

    std::unordered_map<int, std::vector<int>> children;
    children.reserve(nCutE * 2);
    for (int e = 0; e < nCutE; ++e) {
        if (parent[e] >= 0) children[parent[e]].push_back(e);
    }

    // q, its rounding residual and the edge's length, read in the direction the
    // walk is going. The face is the one the frame is taken from: a boundary
    // edge has one, an interface edge has two that agree to the field's own
    // variation across it.
    auto read = [&](const Step &st, long &q, double &residual, double &len) {
        q = 0; residual = M_PI; len = 0.0;
        if (st.orig < 0 || st.from < 0 || st.to < 0) return;
        const int f = (mesh->edgeTriangles[st.orig][0] >= 0) ? mesh->edgeTriangles[st.orig][0]
                                                             : mesh->edgeTriangles[st.orig][1];
        if (f < 0) return;
        const Point d = mesh->vertices[st.to] - mesh->vertices[st.from];
        len = normP(d);
        if (!(len > 0.0)) { len = 0.0; return; }
        // The image direction under J*_f = R(-theta_hat_f)/h, which is a
        // rotation: the edge's own angle less the combed angle.
        const double img = std::atan2(d[1], d[0]) - theta[f];
        q = std::lround(img / M_PI_2);
        residual = std::fabs(img - static_cast<double>(q) * M_PI_2);
    };

    // --- what one chain does once its steps are known ---------------------
    auto commit = [&](const std::vector<Step> &chain, bool onInterface) {
        if (chain.empty()) return;
        double vote[4] = {0.0, 0.0, 0.0, 0.0};
        double total = 0.0;
        for (const Step &st : chain) {
            long q; double res, len;
            read(st, q, res, len);
            if (len <= 0.0) continue;
            total += len;
            if (res > opts.alignmentVoteTolerance) continue;   // too near a coin toss to vote
            vote[((q % 4) + 4) % 4] += len;
        }
        if (!(total > 0.0)) return;

        int winner = -1;
        double best = 0.0;
        for (int k = 0; k < 4; ++k) if (vote[k] > best) { best = vote[k]; winner = k; }
        if (winner < 0) {
            // Every edge abstained: the field is nowhere near an axis on this
            // curve, and there is no chain-wide statement to make about it.
            ++report.alignmentAbstained;
            return;
        }

        // A chain the field turns half way round inside. See the header: the
        // parity is satisfiable and the geometry is a fold.
        const double reversed = vote[(winner + 2) % 4];
        if (reversed > 0.02 * total) {
            ++report.alignmentReversedChains;
            return;
        }
        const double stepped = vote[(winner + 1) % 4] + vote[(winner + 3) % 4];

        const int hold = (winner % 2 == 0) ? 1 : 0;   // even q -> horizontal -> hold v
        ++report.alignmentChains;
        for (const Step &st : chain) {
            long q; double res, len;
            read(st, q, res, len);
            if (len > 0.0 && res <= opts.alignmentVoteTolerance &&
                ((q % 4) + 4) % 4 != winner) {
                ++report.alignmentOverrides;
            }
            if (len > 0.0 && res < M_PI) {
                report.maxAlignmentResidual = std::max(report.maxAlignmentResidual, res);
            }
            for (int ce : st.cutEdges) {
                if (alignAxis[ce] < 0) {
                    alignAxis[ce] = hold;
                    if (onInterface) ++report.alignedInterfaceEdges;
                    else ++report.alignedBoundaryEdges;
                }
                if (stepped <= 0.0 && strictAxis[ce] < 0) strictAxis[ce] = hold;
            }
        }
        if (stepped > 0.0) ++report.alignmentSteppedChains;
    };

    auto stepFor = [&](int orig, int from, int to) {
        Step st;
        st.orig = orig; st.from = from; st.to = to;
        auto it = children.find(orig);
        if (it != children.end()) st.cutEdges = it->second;
        return st;
    };

    // --- the chains of dS - G, walked on Omega ----------------------------
    {
        std::vector<char> onDS(nCutE, 0);
        std::vector<std::vector<int>> at(om.vertices.size());
        for (int e : om.boundaryEdges) {
            const int oe = parent[e];
            if (oe < 0 || !mesh->isBoundaryEdge[oe]) continue;
            onDS[e] = 1;
            at[om.edges[e][0]].push_back(e);
            at[om.edges[e][1]].push_back(e);
        }

        // A vertex the chain has to stop at: where the layout may turn (a cone,
        // a node of the interface network) and where the walk has no single way
        // on (a cut landing, which leaves one dS edge at that child of it, or a
        // pinch, which leaves three).
        auto isBreak = [&](int v) {
            if (at[v].size() != 2) return true;
            const int ov = c2o[v];
            if (ov < 0 || ov >= static_cast<int>(index.size())) return true;
            if (index[ov] != 0) return true;
            if (itf && itf->nodeAt(ov) >= 0) return true;
            return false;
        };

        std::vector<char> used(nCutE, 0);
        auto walkFrom = [&](int startV, int startE) {
            std::vector<Step> chain;
            int v = startV, e = startE;
            while (true) {
                used[e] = 1;
                const int w = (om.edges[e][0] == v) ? om.edges[e][1] : om.edges[e][0];
                if (parent[e] >= 0) chain.push_back(stepFor(parent[e], c2o[v], c2o[w]));
                if (isBreak(w)) break;
                int next = -1;
                for (int f : at[w]) if (f != e && !used[f]) next = f;
                if (next < 0) break;
                v = w; e = next;
            }
            commit(chain, false);
        };

        for (int v = 0; v < static_cast<int>(om.vertices.size()); ++v) {
            if (at[v].empty() || !isBreak(v)) continue;
            for (int e : at[v]) if (!used[e]) walkFrom(v, e);
        }
        // Whatever is left runs round a boundary component carrying no cone and
        // no node. Holding a closed curve to one coordinate line collapses it,
        // so it is left free; a boundary component with no cone on it is Stage
        // 1's business and E2 still has it.
        for (int e = 0; e < nCutE; ++e) {
            if (onDS[e] && !used[e]) { ++report.alignmentClosedChains; used[e] = 1; }
        }
    }

    // --- and on the interface network -------------------------------------
    //
    // A branch already runs node to node, which is the decomposition the walk
    // above had to construct: Interfaces ends a branch at every junction,
    // landing, kink, loop split and balance corner, and those are exactly the
    // points a layout may turn on an interface. So there is nothing to walk --
    // the branches are the chains, and Branch::verts is already in order.
    if (itf) {
        for (const Interfaces::Branch &br : itf->branches()) {
            if (br.closed) { ++report.alignmentClosedChains; continue; }
            if (br.edges.size() + 1 != br.verts.size()) continue;
            std::vector<Step> chain;
            chain.reserve(br.edges.size());
            for (size_t i = 0; i < br.edges.size(); ++i) {
                chain.push_back(stepFor(br.edges[i], br.verts[i], br.verts[i + 1]));
            }
            commit(chain, true);
        }
    }

    if (report.alignmentOverrides > 0) {
        std::ostringstream oss;
        oss << "The frame reads a different axis from its own chain on "
            << report.alignmentOverrides << " edge(s) of dS or of the interface network, over "
            << report.alignmentSteppedChains << " chain(s). That is the staircase of the "
            << "boundary turn audit above, seen per edge: the field slips a quarter where its "
            << "own smoothness puts the slip and not where Stage 1 put a cone. Each is held "
            << "to its chain's axis instead, which is what stops it becoming a node of the "
            << "arrangement Pipeline A does not have -- and it is a right angle the map is "
            << "being asked to unbend, which is what Stage 4F's fallback is for.";
        report.messages.push_back(oss.str());
    }
    if (report.alignmentReversedChains > 0) {
        std::ostringstream oss;
        oss << report.alignmentReversedChains << " chain(s) of dS or of the interface network "
            << "turn half way round inside themselves and were left free. Two edges of one "
            << "chain pointing opposite ways have the same parity, so holding both to one "
            << "coordinate line is satisfiable and folds the chain along itself; the field is "
            << "saying there is a corner of two quarters inside a run Stage 1 called straight, "
            << "and the cone set is where that is fixed.";
        report.messages.push_back(oss.str());
    }
    if (report.alignmentAbstained > 0) {
        std::ostringstream oss;
        oss << report.alignmentAbstained << " chain(s) had no edge within "
            << opts.alignmentVoteTolerance << " rad of an axis and were left free. The field "
            << "is pinned tangent to dS and to the interfaces, so this is a pin that did not "
            << "take -- the usual cause is a triangle carrying two boundary edges that meet at "
            << "an angle that is not a multiple of pi/2, where exp(4 i theta) of the two "
            << "tangents pull against each other and the average is not either of them.";
        report.messages.push_back(oss.str());
    }
}
