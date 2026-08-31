#include "TutteEmbedding.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <vector>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

namespace {

inline double signedArea(const Point &a, const Point &b, const Point &c) {
    return 0.5 * cross2(b - a, c - a);
}

} // namespace

// ---------------------------------------------------------------------------
// boundaryLoop()
//
// Omega's boundary edges walked into one cycle. A disk has exactly one; more
// than one means ConeCut did not close the holes and the theorem does not
// apply, which is reported rather than worked around.
// ---------------------------------------------------------------------------
std::vector<int> TutteEmbedding::boundaryLoop(const Mesh &m) const {
    std::unordered_map<int, std::vector<int>> nbr;
    for (int e : m.boundaryEdges) {
        nbr[m.edges[e][0]].push_back(m.edges[e][1]);
        nbr[m.edges[e][1]].push_back(m.edges[e][0]);
    }
    std::vector<int> loop;
    if (nbr.empty()) return loop;
    for (const auto &kv : nbr) {
        if (kv.second.size() != 2) return loop;   // a pinch; not a simple cycle
    }

    const int start = nbr.begin()->first;
    int prev = -1, cur = start;
    do {
        loop.push_back(cur);
        const std::vector<int> &n = nbr[cur];
        const int next = (n[0] == prev) ? n[1] : n[0];
        prev = cur;
        cur = next;
        if (loop.size() > nbr.size()) return {};   // did not close
    } while (cur != start);

    if (loop.size() != nbr.size()) return {};

    // Orient it the way the model is oriented. The walk above starts at
    // whichever vertex the hash table happened to hand back first and leaves in
    // whichever of the two directions the adjacency listed first, so its sense
    // is arbitrary -- and getting it wrong is not a cosmetic error. A boundary
    // laid onto the circle backwards is a globally reflected embedding, which
    // is still injective and still satisfies the theorem, and in which *every*
    // triangle is negatively oriented. The untangling pass would then start
    // from a map the barrier reads as entirely inverted. Omega's faces are
    // positively oriented, so the outer boundary walked the same way has
    // positive signed area, and that is the test.
    double area2 = 0.0;
    for (size_t i = 0; i < loop.size(); ++i) {
        const Point &a = m.vertices[loop[i]];
        const Point &b = m.vertices[loop[(i + 1) % loop.size()]];
        area2 += cross2(a, b);
    }
    if (area2 < 0.0) std::reverse(loop.begin(), loop.end());
    return loop;
}

// ---------------------------------------------------------------------------
TutteEmbedding::TutteEmbedding(const Mesh &omega, const Options &opts) {
    const int nV = static_cast<int>(omega.vertices.size());
    uv.assign(nV, Point{0.0, 0.0});
    if (nV == 0) return;

    const std::vector<int> loop = boundaryLoop(omega);
    if (loop.size() < 3) {
        report.messages.push_back(
            "The boundary of Omega is not one simple cycle, so Tutte's theorem does not "
            "apply and no embedding was built.");
        return;
    }
    report.boundaryVertices = static_cast<int>(loop.size());
    report.interiorVertices = nV - report.boundaryVertices;
    report.boundaryLoops = 1;

    // --- the convex target ------------------------------------------------
    double modelArea = 0.0;
    for (const Triangle &t : omega.triangles) {
        modelArea += signedArea(omega.vertices[t[0]], omega.vertices[t[1]], omega.vertices[t[2]]);
    }
    modelArea = std::fabs(modelArea);
    double radius = opts.radius;
    if (!(radius > 0.0)) {
        const double h = (opts.targetEdge > 0.0) ? opts.targetEdge : 1.0;
        radius = std::sqrt(std::max(modelArea, 1e-30) / M_PI) / h;
    }
    report.radius = radius;

    // Arc length around dOmega, so that the boundary is stretched as little as
    // the circle allows. Any convex target satisfies the theorem; this one
    // costs the barrier least to unpick afterwards.
    const int nB = static_cast<int>(loop.size());
    std::vector<double> cum(nB + 1, 0.0);
    for (int i = 0; i < nB; ++i) {
        cum[i + 1] = cum[i] + normP(omega.vertices[loop[(i + 1) % nB]] - omega.vertices[loop[i]]);
    }
    const double total = cum[nB] > 0.0 ? cum[nB] : 1.0;

    std::vector<char> onBoundary(nV, 0);
    for (int i = 0; i < nB; ++i) {
        const double a = 2.0 * M_PI * cum[i] / total;
        uv[loop[i]] = Point{radius * std::cos(a), radius * std::sin(a)};
        onBoundary[loop[i]] = 1;
    }

    // --- mean-value weights ----------------------------------------------
    // w_ij = (tan(alpha/2) + tan(beta/2)) / |p_i - p_j|, with alpha and beta the
    // two angles at i between the edge (i, j) and its neighbours in the two
    // incident triangles. Strictly positive on any triangulation, which is what
    // the convex-combination half of the theorem needs and what the cotangent
    // weights do not give on an obtuse triangle.
    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(static_cast<size_t>(nV) * 8);
    Eigen::VectorXd bu = Eigen::VectorXd::Zero(nV), bv = Eigen::VectorXd::Zero(nV);
    std::vector<std::unordered_map<int, double>> w(nV);

    for (const Triangle &t : omega.triangles) {
        for (int i = 0; i < 3; ++i) {
            const int vi = t[i], vj = t[(i + 1) % 3], vk = t[(i + 2) % 3];
            if (onBoundary[vi]) continue;
            const Point a = omega.vertices[vj] - omega.vertices[vi];
            const Point b = omega.vertices[vk] - omega.vertices[vi];
            const double la = normP(a), lb = normP(b);
            if (!(la > 0.0) || !(lb > 0.0)) continue;
            // The angle at vi between the two rays, halved.
            const double ang = std::atan2(std::fabs(cross2(a, b)), dotP(a, b));
            const double th = std::tan(0.5 * std::max(0.0, std::min(M_PI - 1e-12, ang)));
            // Each of the triangle's two rays at vi picks up tan(angle/2) once
            // from this triangle; summed over both incident triangles that is
            // Floater's (tan(alpha/2) + tan(beta/2)).
            w[vi][vj] += th / la;
            w[vi][vk] += th / lb;
        }
    }

    for (int v = 0; v < nV; ++v) {
        if (onBoundary[v]) {
            trip.emplace_back(v, v, 1.0);
            bu[v] = uv[v][0];
            bv[v] = uv[v][1];
            continue;
        }
        double s = 0.0;
        for (const auto &kv : w[v]) s += kv.second;
        if (!(s > 0.0)) {
            // An isolated interior vertex; hold it where it is rather than
            // leaving a zero row.
            trip.emplace_back(v, v, 1.0);
            bu[v] = omega.vertices[v][0];
            bv[v] = omega.vertices[v][1];
            continue;
        }
        trip.emplace_back(v, v, s);
        for (const auto &kv : w[v]) trip.emplace_back(v, kv.first, -kv.second);
    }

    Eigen::SparseMatrix<double> A(nV, nV);
    A.setFromTriplets(trip.begin(), trip.end());
    A.makeCompressed();

    // Not symmetric -- the boundary rows are identities and the mean-value
    // weights are not symmetric anyway -- so an LU, not a Cholesky.
    Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
    lu.compute(A);
    if (lu.info() != Eigen::Success) {
        report.messages.push_back("The Tutte system would not factorise.");
        return;
    }
    const Eigen::VectorXd su = lu.solve(bu);
    const Eigen::VectorXd sv = lu.solve(bv);
    if (lu.info() != Eigen::Success) {
        report.messages.push_back("The Tutte system would not solve.");
        return;
    }
    for (int v = 0; v < nV; ++v) uv[v] = Point{su[v], sv[v]};
    report.solved = true;

    // --- what the theorem promised, measured ------------------------------
    report.minSignedArea = std::numeric_limits<double>::infinity();
    for (const Triangle &t : omega.triangles) {
        const double a = signedArea(uv[t[0]], uv[t[1]], uv[t[2]]);
        report.totalArea += a;
        report.minSignedArea = std::min(report.minSignedArea, a);
        if (!(a > 0.0)) ++report.flippedFaces;
    }
    if (!std::isfinite(report.minSignedArea)) report.minSignedArea = 0.0;

    if (report.flippedFaces > 0) {
        std::ostringstream oss;
        oss << report.flippedFaces << " face(s) came out inverted, which Tutte's theorem says "
            << "cannot happen on a disk with positive convex weights. Either Omega is not the "
            << "disk it was taken for or the solve did not converge; in either case the "
            << "untangling pass has nothing injective to start from.";
        report.messages.push_back(oss.str());
    }
    report.valid = report.solved && report.flippedFaces == 0;
}
