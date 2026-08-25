#include "RicciFlow.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "MERIDIAN/ConeSingularities.hxx"

namespace {

inline double dist(const Mesh &m, int a, int b) {
    const Point &pa = m.vertices[a];
    const Point &pb = m.vertices[b];
    const double dx = pa[0] - pb[0];
    const double dy = pa[1] - pb[1];
    return std::sqrt(dx * dx + dy * dy);
}

inline double safeAcos(double c) {
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    return std::acos(c);
}

// Four-point Gauss-Legendre on [0, 1]: nodes and weights.
const double kGLNode[4] = {0.0694318442029737, 0.3300094782075719,
                           0.6699905217924281, 0.9305681557970263};
const double kGLWeight[4] = {0.1739274225687269, 0.3260725774312731,
                             0.3260725774312731, 0.1739274225687269};

} // namespace

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
RicciFlow::RicciFlow(std::shared_ptr<Mesh> m, const ConeSingularities &cones)
    : RicciFlow(std::move(m), cones.getTargetCurvature()) {}

RicciFlow::RicciFlow(std::shared_ptr<Mesh> m, const std::vector<double> &targetCurvature)
    : mesh(std::move(m)) {
    if (!mesh) throw std::runtime_error("RicciFlow: null mesh");
    if (mesh->triangles.empty()) throw std::runtime_error("RicciFlow: empty mesh");
    if (targetCurvature.size() != mesh->vertices.size()) {
        throw std::runtime_error("RicciFlow: target curvature has " +
                                 std::to_string(targetCurvature.size()) + " entries for " +
                                 std::to_string(mesh->vertices.size()) + " vertices");
    }
    target = targetCurvature;
    faces = mesh->triangles;
    initialiseMetric();
}

// ---------------------------------------------------------------------------
// initialiseMetric()
//
// gamma_i^{ijk} = ( l_ij + l_ik - l_jk ) / 2 on each face, gamma_i the smallest
// over the faces around v_i; then cos(phi_ij) from Eq. (11) inverted. Taking
// the minimum is what guarantees gamma_i + gamma_j <= l_ij on every face, which
// in turn puts cos(phi_ij) at or above 1 everywhere: the circles overlap in the
// inversive sense from the very first iterate. That is exactly the regime
// Thurston's packing excludes and Yang et al. [87] admit, so nothing is clamped
// here -- clamping would silently replace the input metric with a different one
// before the flow ever starts.
// ---------------------------------------------------------------------------
void RicciFlow::initialiseMetric() {
    const int nV = static_cast<int>(mesh->vertices.size());

    active.assign(nV, 0);
    for (const Triangle &t : faces) {
        active[t[0]] = 1;
        active[t[1]] = 1;
        active[t[2]] = 1;
    }
    activeCount = 0;
    for (int v = 0; v < nV; ++v) if (active[v]) ++activeCount;
    if (activeCount != nV) {
        std::ostringstream oss;
        oss << (nV - activeCount) << " vertex/vertices carry no triangle and are held "
            << "out of the flow.";
        report.messages.push_back(oss.str());
    }

    std::vector<double> gamma(nV, std::numeric_limits<double>::infinity());
    for (const Triangle &t : faces) {
        const double lab = dist(*mesh, t[0], t[1]);
        const double lbc = dist(*mesh, t[1], t[2]);
        const double lca = dist(*mesh, t[2], t[0]);
        gamma[t[0]] = std::min(gamma[t[0]], 0.5 * (lab + lca - lbc));
        gamma[t[1]] = std::min(gamma[t[1]], 0.5 * (lab + lbc - lca));
        gamma[t[2]] = std::min(gamma[t[2]], 0.5 * (lbc + lca - lab));
    }

    u.assign(nV, 0.0);
    for (int v = 0; v < nV; ++v) {
        if (!active[v]) { gamma[v] = 1.0; u[v] = 0.0; continue; }
        if (!(gamma[v] > 0.0) || !std::isfinite(gamma[v])) {
            // A degenerate face: the triangle inequality is an equality there,
            // so the inscribed radius collapses. Fall back to something
            // positive so the solve can start and the failure is reported by
            // the realisability count rather than by a NaN.
            gamma[v] = 1e-12;
            std::ostringstream oss;
            oss << "Vertex " << v << " has a non-positive circle radius; the input "
                << "triangulation has a degenerate face there.";
            report.messages.push_back(oss.str());
        }
        u[v] = std::log(gamma[v]);
    }

    cosPhi.clear();
    for (const Triangle &t : faces) {
        for (int p = 0; p < 3; ++p) {
            const int a = t[p];
            const int b = t[(p + 1) % 3];
            const EdgeKey key(a, b);
            if (cosPhi.count(key)) continue;
            const double l = dist(*mesh, a, b);
            cosPhi[key] = (l * l - gamma[a] * gamma[a] - gamma[b] * gamma[b]) /
                          (2.0 * gamma[a] * gamma[b]);
        }
    }

    initialCosPhi = cosPhi;

    curvature.assign(nV, 0.0);
    angleSum.assign(nV, 0.0);

    rebuildEdges();
    rebuildLengths();
    rebuildAnglesAndCurvature();
    rebuildWeights();
}

// ---------------------------------------------------------------------------
// rebuildEdges()
//
// The edge table is derived from `faces` rather than maintained through the
// flips, so a flip only has to fix up the two face triples and drop one
// cos(phi) in for the new diagonal. (i, j) is stored in the orientation face f0
// sees, which is what lets tryFlip() name the quad's corners without redoing
// the orientation test.
// ---------------------------------------------------------------------------
void RicciFlow::rebuildEdges() {
    const int nF = static_cast<int>(faces.size());
    edges.clear();
    edgeIndex.clear();
    faceEdges.assign(nF, std::array<int, 3>{-1, -1, -1});
    edges.reserve(static_cast<size_t>(nF) * 2);

    for (int f = 0; f < nF; ++f) {
        const Triangle &t = faces[f];
        for (int p = 0; p < 3; ++p) {
            const int a = t[p];
            const int b = t[(p + 1) % 3];
            const int c = t[(p + 2) % 3];
            const EdgeKey key(a, b);

            auto it = edgeIndex.find(key);
            if (it == edgeIndex.end()) {
                FlowEdge fe;
                fe.i = a;
                fe.j = b;
                fe.f0 = f;
                fe.k = c;
                auto cp = cosPhi.find(key);
                fe.cosPhi = (cp == cosPhi.end()) ? 1.0 : cp->second;
                const int idx = static_cast<int>(edges.size());
                edges.push_back(fe);
                edgeIndex.emplace(key, idx);
                faceEdges[f][p] = idx;
            } else {
                FlowEdge &fe = edges[it->second];
                fe.f1 = f;
                fe.l = c;
                faceEdges[f][p] = it->second;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// lengthsAt() / rebuildLengths()  --  Eq. (11)
// ---------------------------------------------------------------------------
std::vector<double> RicciFlow::lengthsAt(const std::vector<double> &uTrial) const {
    std::vector<double> len(edges.size(), 0.0);
    for (size_t e = 0; e < edges.size(); ++e) {
        const FlowEdge &fe = edges[e];
        const double gi = std::exp(uTrial[fe.i]);
        const double gj = std::exp(uTrial[fe.j]);
        const double sq = gi * gi + gj * gj + 2.0 * gi * gj * fe.cosPhi;
        len[e] = (sq > 0.0) ? std::sqrt(sq) : 0.0;
    }
    return len;
}

void RicciFlow::rebuildLengths() {
    const std::vector<double> len = lengthsAt(u);
    for (size_t e = 0; e < edges.size(); ++e) edges[e].length = len[e];
}

// ---------------------------------------------------------------------------
// curvatureFrom()  --  Eqs. (5) and (6)
// ---------------------------------------------------------------------------
void RicciFlow::curvatureFrom(const std::vector<double> &len,
                              std::vector<double> &K,
                              std::vector<double> &angles) const {
    const int nV = static_cast<int>(mesh->vertices.size());
    angles.assign(nV, 0.0);

    for (int f = 0; f < static_cast<int>(faces.size()); ++f) {
        const Triangle &t = faces[f];
        for (int m = 0; m < 3; ++m) {
            // At local vertex m the two incident sides are local edges m and
            // (m+2)%3; the side opposite it is local edge (m+1)%3.
            const double a = len[faceEdges[f][m]];
            const double b = len[faceEdges[f][(m + 2) % 3]];
            const double c = len[faceEdges[f][(m + 1) % 3]];
            if (a <= 0.0 || b <= 0.0) continue;
            angles[t[m]] += safeAcos((a * a + b * b - c * c) / (2.0 * a * b));
        }
    }

    K.assign(nV, 0.0);
    for (int v = 0; v < nV; ++v) {
        if (!active[v]) continue; // not on the surface, so no curvature
        const double full = mesh->isBoundaryVertex[v] ? M_PI : 2.0 * M_PI;
        K[v] = full - angles[v];
    }
}

void RicciFlow::rebuildAnglesAndCurvature() {
    std::vector<double> len(edges.size());
    for (size_t e = 0; e < edges.size(); ++e) len[e] = edges[e].length;
    curvatureFrom(len, curvature, angleSum);
}

// ---------------------------------------------------------------------------
// powerHeight()
//
// Put a at the origin and b on the positive x-axis. The power (radical) centre
// O of the three circles satisfies |O - v|^2 - gamma_v^2 equal for all three,
// which pins it to
//
//     x = ( l_ab^2 + ga^2 - gb^2 ) / ( 2 l_ab ),
//     x cos(alpha_a) + y sin(alpha_a) = ( l_ca^2 + ga^2 - gc^2 ) / ( 2 l_ca ),
//
// and h_ab^c is that y -- positive when the power centre lies on the same side
// of ab as c. With all three radii equal the power centre is the circumcentre
// and h_ab^c = (l_ab / 2) cot(alpha_c), so w_ab collapses to the cotangent
// weight (cot alpha_c + cot alpha_d) / 2 and the weighted-Delaunay test
// h_ab^c + h_ab^d >= 0 collapses to alpha_c + alpha_d <= pi. Those are the two
// checks to reach for if this formula is ever in doubt.
// ---------------------------------------------------------------------------
double RicciFlow::powerHeight(double lab, double lbc, double lca,
                              double ga, double gb, double gc) {
    if (lab <= 0.0 || lca <= 0.0) return 0.0;
    const double cosA = (lab * lab + lca * lca - lbc * lbc) / (2.0 * lab * lca);
    const double alpha = safeAcos(cosA);
    const double sinA = std::sin(alpha);
    if (sinA < 1e-14) return 0.0; // degenerate corner: no usable power centre

    const double dab = (lab * lab + ga * ga - gb * gb) / (2.0 * lab);
    const double dac = (lca * lca + ga * ga - gc * gc) / (2.0 * lca);
    return (dac - dab * std::cos(alpha)) / sinA;
}

// ---------------------------------------------------------------------------
// rebuildWeights()
//
//     w_ij = ( h_ij^k + h_ij^l ) / l_ij
//
// with only the one term on a boundary edge, exactly as the cotangent Laplacian
// keeps only one cotangent there.
// ---------------------------------------------------------------------------
void RicciFlow::rebuildWeights() {
    double minInterior = std::numeric_limits<double>::infinity();
    double minBoundary = std::numeric_limits<double>::infinity();

    auto lengthOf = [&](int a, int b) -> double {
        auto it = edgeIndex.find(EdgeKey(a, b));
        return (it == edgeIndex.end()) ? 0.0 : edges[it->second].length;
    };

    for (FlowEdge &fe : edges) {
        const double gi = std::exp(u[fe.i]);
        const double gj = std::exp(u[fe.j]);
        double h = 0.0;

        if (fe.k >= 0) {
            const double gk = std::exp(u[fe.k]);
            h += powerHeight(fe.length, lengthOf(fe.j, fe.k), lengthOf(fe.k, fe.i), gi, gj, gk);
        }
        if (fe.f1 >= 0 && fe.l >= 0) {
            const double gl = std::exp(u[fe.l]);
            h += powerHeight(fe.length, lengthOf(fe.j, fe.l), lengthOf(fe.l, fe.i), gi, gj, gl);
        }
        fe.weight = (fe.length > 0.0) ? h / fe.length : 0.0;
        double &slot = (fe.f1 >= 0) ? minInterior : minBoundary;
        slot = std::min(slot, fe.weight);
    }

    report.minInteriorEdgeWeight = std::isfinite(minInterior) ? minInterior : 0.0;
    report.minBoundaryEdgeWeight = std::isfinite(minBoundary) ? minBoundary : 0.0;
}

// ---------------------------------------------------------------------------
// realisable()
// ---------------------------------------------------------------------------
bool RicciFlow::realisable(const std::vector<double> &len, int *badFace) const {
    for (int f = 0; f < static_cast<int>(faces.size()); ++f) {
        const double a = len[faceEdges[f][0]];
        const double b = len[faceEdges[f][1]];
        const double c = len[faceEdges[f][2]];
        if (!(a > 0.0 && b > 0.0 && c > 0.0) ||
            a + b <= c || b + c <= a || c + a <= b) {
            if (badFace) *badFace = f;
            return false;
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// tryFlip()  --  Sec. 3.2.1 / [87]
//
// The two triangles either side of a non-Delaunay edge are unfolded into the
// plane, the other diagonal of the resulting quad is measured, and that length
// becomes the new edge -- from which a new cos(phi) follows, since the radii do
// not change. The flip is refused if the quad is not convex (the new diagonal
// would fall outside it) or if either new triangle would not close, both of
// which mean the metric is too distorted for this particular edge to be the
// culprit.
// ---------------------------------------------------------------------------
bool RicciFlow::tryFlip(int edgeIdx) {
    FlowEdge fe = edges[edgeIdx];
    if (fe.f1 < 0 || fe.k < 0 || fe.l < 0) return false;  // boundary edge
    if (fe.k == fe.l) return false;                       // would fold the mesh
    // cosPhi is keyed on exactly the current edge set (erased and inserted as
    // flips happen), so it is the live test for "that diagonal already exists"
    // even part-way through a pass, when the edge table is deliberately stale.
    if (cosPhi.count(EdgeKey(fe.k, fe.l))) return false;

    auto lengthOf = [&](int a, int b) -> double {
        auto it = edgeIndex.find(EdgeKey(a, b));
        return (it == edgeIndex.end()) ? -1.0 : edges[it->second].length;
    };

    const double lij = fe.length;
    const double lik = lengthOf(fe.i, fe.k), ljk = lengthOf(fe.j, fe.k);
    const double lil = lengthOf(fe.i, fe.l), ljl = lengthOf(fe.j, fe.l);
    if (lij <= 0.0 || lik <= 0.0 || ljk <= 0.0 || lil <= 0.0 || ljl <= 0.0) return false;

    // Unfold: i at the origin, j on the x-axis, k above, l below.
    const double xk = (lik * lik + lij * lij - ljk * ljk) / (2.0 * lij);
    const double yk2 = lik * lik - xk * xk;
    const double xl = (lil * lil + lij * lij - ljl * ljl) / (2.0 * lij);
    const double yl2 = lil * lil - xl * xl;
    if (yk2 <= 0.0 || yl2 <= 0.0) return false;
    const double yk = std::sqrt(yk2);
    const double yl = -std::sqrt(yl2);

    // Convexity: the segment kl has to cross the segment ij strictly inside it.
    const double cross = xk + (xl - xk) * (yk / (yk - yl));
    if (!(cross > 1e-12 && cross < lij - 1e-12)) return false;

    const double lkl = std::hypot(xk - xl, yk - yl);
    if (!(lkl > 0.0)) return false;

    // The two triangles the flip produces must themselves close.
    auto closes = [](double a, double b, double c) {
        return a + b > c && b + c > a && c + a > b;
    };
    if (!closes(lik, lil, lkl) || !closes(ljk, ljl, lkl)) return false;

    const double gk = std::exp(u[fe.k]);
    const double gl = std::exp(u[fe.l]);
    const double newCos = (lkl * lkl - gk * gk - gl * gl) / (2.0 * gk * gl);
    // Eq. (11) has to stay solvable for a positive length at any u, which needs
    // cos(phi) > -1. Anything at or below that is a metric this flip cannot
    // represent, so the edge is left alone and the line search deals with it.
    if (!(newCos > -1.0 + 1e-9) || !std::isfinite(newCos)) return false;

    // (i, j, k) and (j, i, l) become (k, i, l) and (k, l, j), which is the
    // orientation-preserving relabelling of the quad i -> l -> j -> k.
    faces[fe.f0] = Triangle{fe.k, fe.i, fe.l};
    faces[fe.f1] = Triangle{fe.k, fe.l, fe.j};

    cosPhi.erase(EdgeKey(fe.i, fe.j));
    cosPhi[EdgeKey(fe.k, fe.l)] = newCos;
    return true;
}

// ---------------------------------------------------------------------------
// flipPass()
//
// Flip every non-Delaunay edge whose two triangles this pass has not already
// disturbed, then rebuild once. Rebuilding after each individual flip would be
// correct too, and quadratic; the restriction to untouched faces is what makes
// a whole sweep safe on one snapshot, since every length tryFlip() reads
// belongs to the two faces it is about to replace. Passes repeat until nothing
// moves, which is what the flip algorithm terminating means in practice.
// ---------------------------------------------------------------------------
int RicciFlow::flipPass() {
    int total = 0;
    const int cap = 8 * static_cast<int>(edges.size()) + 64;

    for (int sweep = 0; sweep < 64; ++sweep) {
        std::vector<char> touched(faces.size(), 0);
        int flips = 0;

        for (int e = 0; e < static_cast<int>(edges.size()); ++e) {
            const FlowEdge &fe = edges[e];
            if (fe.weight >= 0.0) continue;
            if (fe.f0 < 0 || fe.f1 < 0) continue;
            if (touched[fe.f0] || touched[fe.f1]) continue;

            const int f0 = fe.f0, f1 = fe.f1;
            if (!tryFlip(e)) continue;
            touched[f0] = 1;
            touched[f1] = 1;
            ++flips;
            if (total + flips >= cap) break;
        }

        if (flips == 0) break;
        total += flips;

        rebuildEdges();
        rebuildLengths();
        rebuildWeights();

        if (total >= cap) {
            report.messages.push_back(
                "Weighted-Delaunay flipping hit its cap of " + std::to_string(cap) +
                " flips; the metric may still have negative edge weights.");
            break;
        }
    }
    return total;
}

// ---------------------------------------------------------------------------
// energyDelta()
// ---------------------------------------------------------------------------
double RicciFlow::energyDelta(const std::vector<double> &du, double lambda) const {
    const int nV = static_cast<int>(u.size());
    std::vector<double> uT(nV), K, ang;
    double acc = 0.0;

    for (int q = 0; q < 4; ++q) {
        const double t = lambda * kGLNode[q];
        for (int v = 0; v < nV; ++v) uT[v] = u[v] + t * du[v];
        curvatureFrom(lengthsAt(uT), K, ang);

        double dot = 0.0;
        for (int v = 0; v < nV; ++v) dot += (K[v] - target[v]) * du[v];
        acc += kGLWeight[q] * dot;
    }
    return lambda * acc;
}

// ---------------------------------------------------------------------------
// newtonDirection()
//
// Assemble Delta, drop the pinned row and column, solve Delta du = Kbar - K.
//
// The sign is worth being explicit about, because getting it backwards still
// looks plausible. Raising u_i lengthens every edge at v_i without touching the
// opposite sides, which narrows the angle at v_i and so *raises* K_i: dK/du is
// the positive-diagonal Laplacian Delta, not its negative. The Newton step for
// driving K to Kbar is therefore Delta du = Kbar - K, and du carries the same
// sign as the curvature deficit. checkHessian() is what pins this down
// numerically.
//
// Delta annihilates the constants -- u -> u + c is a global rescaling of the
// metric and changes no angle -- so one vertex is pinned. That also means the
// residual has to be orthogonal to the constants for the system to be
// consistent, i.e. sum(Kbar) = sum(K) = 2 pi chi. Gauss-Bonnet is not a
// formality here; it is the solvability condition, which is why solve() checks
// it before starting.
// ---------------------------------------------------------------------------
bool RicciFlow::newtonDirection(std::vector<double> &du) {
    const int nV = static_cast<int>(u.size());
    du.assign(nV, 0.0);

    std::vector<int> reduced(nV, -1);
    int n = 0;
    for (int v = 0; v < nV; ++v) if (active[v] && v != pinnedVertex) reduced[v] = n++;
    if (n == 0) return false;

    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(edges.size() * 4 + nV);
    std::vector<double> diag(nV, 0.0);

    for (const FlowEdge &fe : edges) {
        const double w = fe.weight;
        diag[fe.i] += w;
        diag[fe.j] += w;
        if (reduced[fe.i] >= 0 && reduced[fe.j] >= 0) {
            trips.emplace_back(reduced[fe.i], reduced[fe.j], -w);
            trips.emplace_back(reduced[fe.j], reduced[fe.i], -w);
        }
    }
    for (int v = 0; v < nV; ++v) {
        if (reduced[v] >= 0) trips.emplace_back(reduced[v], reduced[v], diag[v]);
    }

    Eigen::SparseMatrix<double> L(n, n);
    L.setFromTriplets(trips.begin(), trips.end());

    Eigen::VectorXd rhs(n);
    for (int v = 0; v < nV; ++v) {
        if (reduced[v] >= 0) rhs[reduced[v]] = target[v] - curvature[v];
    }

    Eigen::VectorXd sol;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
    ldlt.compute(L);
    if (ldlt.info() == Eigen::Success) {
        sol = ldlt.solve(rhs);
    }
    if (ldlt.info() != Eigen::Success || !sol.allFinite()) {
        // A negative weight survived the flip pass, so Delta is not positive
        // definite and the Cholesky gives up. LU still solves it, and the line
        // search below is what catches a direction that is not a descent one.
        Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
        lu.compute(L);
        if (lu.info() != Eigen::Success) return false;
        sol = lu.solve(rhs);
        if (!sol.allFinite()) return false;
    }

    for (int v = 0; v < nV; ++v) if (reduced[v] >= 0) du[v] = sol[reduced[v]];
    return true;
}

// ---------------------------------------------------------------------------
// solve()
//
//   fix u_0 = 0 on one vertex
//   repeat:
//       rebuild l from u; rebuild angles; rebuild K
//       if ||K - Kbar||_inf < tol: stop
//       restore weighted-Delaunay by edge flips
//       assemble Delta; solve Delta du = Kbar - K
//       backtrack lambda until u + lambda du is realisable and the energy drops
//       u <- u + lambda du
// ---------------------------------------------------------------------------
bool RicciFlow::solve() {
    const int nV = static_cast<int>(u.size());
    if (pinnedVertex < 0 || pinnedVertex >= nV || !active[pinnedVertex]) {
        pinnedVertex = 0;
        while (pinnedVertex < nV && !active[pinnedVertex]) ++pinnedVertex;
        if (pinnedVertex >= nV) {
            report.messages.push_back("No vertex carries a triangle.");
            return false;
        }
    }

    // Solvability: sum(Kbar) must equal 2 pi chi, Eq. (7). chi counts only the
    // vertices that are actually on the surface.
    const int chi = activeCount -
                    static_cast<int>(mesh->edges.size()) +
                    static_cast<int>(mesh->triangles.size());
    double targetSum = 0.0;
    for (int v = 0; v < nV; ++v) if (active[v]) targetSum += target[v];
    report.gaussBonnetResidual = std::fabs(targetSum - 2.0 * M_PI * chi);
    if (report.gaussBonnetResidual > 1e-6) {
        std::ostringstream oss;
        oss << "sum(Kbar) = " << targetSum << " but 2 pi chi = " << (2.0 * M_PI * chi)
            << " (residual " << report.gaussBonnetResidual
            << "); the cone set is not admissible and the Newton system is inconsistent.";
        report.messages.push_back(oss.str());
        report.converged = false;
        report.finalError = maxCurvatureError();
        report.initialError = report.finalError;
        return false;
    }

    rebuildLengths();
    rebuildAnglesAndCurvature();
    rebuildWeights();
    report.initialError = maxCurvatureError();

    std::vector<double> du;
    bool converged = false;

    for (int iter = 0; iter < maxIterations; ++iter) {
        rebuildLengths();
        rebuildAnglesAndCurvature();

        if (maxCurvatureError() < tolerance) { converged = true; break; }

        if (delaunayFlips) {
            const int flips = flipPass();
            if (flips > 0) {
                ++report.flipPasses;
                report.flips += flips;
                // A flip leaves the metric and every K alone, but the Hessian
                // is a different operator on a different graph, so the angles
                // are refreshed before it is assembled.
                rebuildAnglesAndCurvature();
            }
        }
        rebuildWeights();

        if (!newtonDirection(du)) {
            report.messages.push_back("Newton system could not be factorised.");
            break;
        }

        // The gradient of Eq. (10) is K - Kbar, so a descent direction has
        // (K - Kbar) . du < 0. Delta is positive semi-definite whenever the
        // triangulation is weighted-Delaunay, which makes the Newton direction
        // a descent one; if a negative weight survived, fall back to the flow
        // of Eq. (9) itself, du = Kbar - K, which always is.
        double slope = 0.0;
        for (int v = 0; v < nV; ++v) slope += (curvature[v] - target[v]) * du[v];
        if (!(slope < 0.0)) {
            for (int v = 0; v < nV; ++v) du[v] = active[v] ? target[v] - curvature[v] : 0.0;
            du[pinnedVertex] = 0.0;
            slope = 0.0;
            for (int v = 0; v < nV; ++v) slope += (curvature[v] - target[v]) * du[v];
            if (!(slope < 0.0)) {
                report.messages.push_back("No descent direction at iteration " +
                                          std::to_string(iter) + ".");
                break;
            }
        }

        double lambda = 1.0;
        std::vector<double> uTrial(nV);
        bool stepped = false;
        for (int back = 0; back < 60; ++back) {
            for (int v = 0; v < nV; ++v) uTrial[v] = u[v] + lambda * du[v];
            if (realisable(lengthsAt(uTrial)) &&
                energyDelta(du, lambda) <= 1e-4 * lambda * slope) {
                u.swap(uTrial);
                stepped = true;
                break;
            }
            lambda *= 0.5;
        }

        if (!stepped) {
            report.messages.push_back(
                "Line search failed at iteration " + std::to_string(iter) +
                ": no realisable step with a decrease in the Ricci energy.");
            break;
        }
        ++report.newtonIterations;
    }

    // One last restoration pass. The loop breaks on convergence *before* it
    // flips, so the final u is otherwise left on whatever triangulation the
    // previous step happened to use, and an edge that only went negative on the
    // last step would never be fixed. A flip changes neither the metric nor any
    // K -- it re-cuts one planar quad along its other diagonal -- so this
    // cannot undo the convergence it comes after; it only leaves Stage 4 a
    // weighted-Delaunay triangulation to lay out.
    rebuildLengths();
    rebuildAnglesAndCurvature();
    rebuildWeights();
    if (delaunayFlips) {
        const int flips = flipPass();
        if (flips > 0) {
            ++report.flipPasses;
            report.flips += flips;
            rebuildAnglesAndCurvature();
            rebuildWeights();
        }
    }
    if (!converged) converged = maxCurvatureError() < tolerance;

    // u is only determined up to an additive constant; pinning it makes the
    // reported factors reproducible.
    const double shift = u[pinnedVertex];
    for (int v = 0; v < nV; ++v) u[v] = active[v] ? u[v] - shift : 0.0;
    rebuildLengths();
    rebuildWeights();

    report.converged = converged;
    report.finalError = maxCurvatureError();

    int bad = -1;
    report.nonRealisableFaces = 0;
    {
        std::vector<double> len(edges.size());
        for (size_t e = 0; e < edges.size(); ++e) len[e] = edges[e].length;
        for (int f = 0; f < static_cast<int>(faces.size()); ++f) {
            const double a = len[faceEdges[f][0]];
            const double b = len[faceEdges[f][1]];
            const double c = len[faceEdges[f][2]];
            if (!(a > 0.0 && b > 0.0 && c > 0.0) || a + b <= c || b + c <= a || c + a <= b) {
                ++report.nonRealisableFaces;
                bad = f;
            }
        }
    }
    if (report.nonRealisableFaces > 0) {
        std::ostringstream oss;
        oss << report.nonRealisableFaces << " face(s) violate the triangle inequality "
            << "in the final metric (e.g. face " << bad << ").";
        report.messages.push_back(oss.str());
    }

    report.minInversiveDistance = std::numeric_limits<double>::infinity();
    report.maxInversiveDistance = -std::numeric_limits<double>::infinity();
    for (const FlowEdge &fe : edges) {
        report.minInversiveDistance = std::min(report.minInversiveDistance, fe.cosPhi);
        report.maxInversiveDistance = std::max(report.maxInversiveDistance, fe.cosPhi);
    }
    if (edges.empty()) { report.minInversiveDistance = 0.0; report.maxInversiveDistance = 0.0; }

    report.flippedAwayOriginalEdges = 0;
    for (const auto &e : mesh->edges) {
        if (!edgeIndex.count(EdgeKey(e[0], e[1]))) ++report.flippedAwayOriginalEdges;
    }

    if (!report.converged) {
        std::ostringstream oss;
        oss << "Ricci flow stopped at ||K - Kbar||_inf = " << report.finalError
            << " after " << report.newtonIterations << " Newton step(s).";
        report.messages.push_back(oss.str());
    }
    return report.converged;
}

// ---------------------------------------------------------------------------
// Accessors
// ---------------------------------------------------------------------------
std::vector<double> RicciFlow::getRadii() const {
    std::vector<double> g(u.size());
    for (size_t v = 0; v < u.size(); ++v) g[v] = std::exp(u[v]);
    return g;
}

double RicciFlow::maxCurvatureError() const {
    double e = 0.0;
    for (size_t v = 0; v < curvature.size(); ++v) {
        if (!active[v]) continue;
        e = std::max(e, std::fabs(curvature[v] - target[v]));
    }
    return e;
}

std::unordered_map<RicciFlow::EdgeKey, double, RicciFlow::EdgeKeyHash>
RicciFlow::flowEdgeLengths() const {
    std::unordered_map<EdgeKey, double, EdgeKeyHash> out;
    out.reserve(edges.size());
    for (const FlowEdge &fe : edges) out.emplace(EdgeKey(fe.i, fe.j), fe.length);
    return out;
}

std::vector<double> RicciFlow::originalEdgeLengths() const {
    std::vector<double> out(mesh->edges.size(), -1.0);
    for (size_t e = 0; e < mesh->edges.size(); ++e) {
        auto it = edgeIndex.find(EdgeKey(mesh->edges[e][0], mesh->edges[e][1]));
        if (it != edgeIndex.end()) out[e] = edges[it->second].length;
    }
    return out;
}

std::vector<double> RicciFlow::originalEdgeLengthsCompleted(int *recovered) const {
    std::vector<double> out = originalEdgeLengths();
    int filled = 0;
    for (size_t e = 0; e < out.size(); ++e) {
        if (out[e] >= 0.0) continue;
        const EdgeKey key(mesh->edges[e][0], mesh->edges[e][1]);
        auto it = initialCosPhi.find(key);
        if (it == initialCosPhi.end()) continue;   // not an edge of this mesh at all
        // Eq. (11) on the edge's own inversive distance and the converged u.
        const double gi = std::exp(u[key.a]);
        const double gj = std::exp(u[key.b]);
        const double sq = gi * gi + gj * gj + 2.0 * gi * gj * it->second;
        if (!(sq > 0.0)) continue;
        out[e] = std::sqrt(sq);
        ++filled;
    }
    if (recovered) *recovered = filled;
    return out;
}

// ---------------------------------------------------------------------------
// checkHessian()
//
// dK_i/du_j against a central difference. Sampling a spread of vertices rather
// than all of them keeps this usable inside a test on a 20k-vertex mesh; the
// error reported is the largest relative disagreement seen, normalised by the
// size of the row so that a vertex whose whole row is tiny does not dominate.
// ---------------------------------------------------------------------------
double RicciFlow::checkHessian(int samples, double h) const {
    const int nV = static_cast<int>(u.size());
    if (nV == 0 || samples <= 0) return 0.0;

    // The analytic row for v: Delta_vv = sum_j w_vj, Delta_vj = -w_vj.
    std::vector<std::vector<std::pair<int, double>>> incident(nV);
    std::vector<double> diag(nV, 0.0);
    for (const FlowEdge &fe : edges) {
        incident[fe.i].push_back({fe.j, -fe.weight});
        incident[fe.j].push_back({fe.i, -fe.weight});
        diag[fe.i] += fe.weight;
        diag[fe.j] += fe.weight;
    }

    const int stride = std::max(1, nV / samples);
    double worst = 0.0;
    std::vector<double> uP(nV), uM(nV), Kp, Km, ang;

    for (int v = 0; v < nV; v += stride) {
        if (!active[v]) continue;
        uP = u; uM = u;
        uP[v] += h;
        uM[v] -= h;
        curvatureFrom(lengthsAt(uP), Kp, ang);
        curvatureFrom(lengthsAt(uM), Km, ang);

        // Column v of dK/du: (K(u + h e_v) - K(u - h e_v)) / 2h, compared with
        // column v of Delta, which by symmetry is row v.
        double scale = std::fabs(diag[v]);
        for (const auto &p : incident[v]) scale = std::max(scale, std::fabs(p.second));
        if (scale < 1e-12) continue;

        auto compare = [&](int i, double analytic) {
            const double numeric = (Kp[i] - Km[i]) / (2.0 * h);
            worst = std::max(worst, std::fabs(numeric - analytic) / scale);
        };
        compare(v, diag[v]);
        for (const auto &p : incident[v]) compare(p.first, p.second);
    }
    return worst;
}
