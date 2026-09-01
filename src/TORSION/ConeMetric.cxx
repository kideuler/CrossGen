#include "ConeMetric.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include "MERIDIAN/ConeSingularities.hxx"

namespace {

inline double safeAcos(double c) {
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    return std::acos(c);
}

// The three angles of a triangle from its three sides, with side[i] joining
// corner i and corner i+1 -- Mesh::triangleEdges' own convention.
inline void anglesOf(const double l[3], double a[3]) {
    for (int i = 0; i < 3; ++i) {
        const double p = l[i], q = l[(i + 2) % 3], r = l[(i + 1) % 3];
        if (!(p > 0.0) || !(q > 0.0)) { a[i] = 0.0; continue; }
        a[i] = safeAcos((p * p + q * q - r * r) / (2.0 * p * q));
    }
}

// Strictly realisable: every side shorter than the sum of the other two, by a
// margin taken relative to the longest side so the test is scale free.
inline bool realisable(const double l[3], double margin) {
    const double m = std::max(l[0], std::max(l[1], l[2]));
    if (!(m > 0.0)) return false;
    for (int i = 0; i < 3; ++i) {
        if (l[i] >= l[(i + 1) % 3] + l[(i + 2) % 3] - margin * m) return false;
    }
    return true;
}

} // namespace

ConeMetric::ConeMetric(const Mesh &mesh, const ConeSingularities &cones, const Options &options)
    : opts(options) {
    const int nV = static_cast<int>(mesh.vertices.size());
    const int nE = static_cast<int>(mesh.edges.size());
    u.assign(nV, 0.0);
    baseLen.assign(nE, 0.0);
    for (int e = 0; e < nE; ++e) {
        baseLen[e] = normP(mesh.vertices[mesh.edges[e][1]] - mesh.vertices[mesh.edges[e][0]]);
    }
    rebuild(mesh);
    solve(mesh, cones);
}

void ConeMetric::rebuild(const Mesh &mesh) {
    const int nE = static_cast<int>(mesh.edges.size());
    len.assign(nE, 0.0);
    for (int e = 0; e < nE; ++e) {
        const int a = mesh.edges[e][0], b = mesh.edges[e][1];
        len[e] = std::exp(0.5 * (u[a] + u[b])) * baseLen[e];
    }
    const int nT = static_cast<int>(mesh.triangles.size());
    fscale.assign(nT, 1.0);
    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = mesh.triangles[t];
        fscale[t] = std::exp((u[tri[0]] + u[tri[1]] + u[tri[2]]) / 3.0);
    }
    report.minScale = std::numeric_limits<double>::infinity();
    report.maxScale = 0.0;
    for (int v = 0; v < static_cast<int>(u.size()); ++v) {
        const double s = std::exp(u[v]);
        report.minScale = std::min(report.minScale, s);
        report.maxScale = std::max(report.maxScale, s);
    }
    if (!std::isfinite(report.minScale)) report.minScale = 1.0;
}

std::vector<double> ConeMetric::sizingField(double targetEdge) const {
    std::vector<double> h(fscale.size(), targetEdge);
    for (size_t t = 0; t < fscale.size(); ++t) {
        if (fscale[t] > 0.0 && std::isfinite(fscale[t])) h[t] = targetEdge / fscale[t];
    }
    return h;
}

// ---------------------------------------------------------------------------
// solve()
//
// Newton on the vertex scaling, from u = 0. Each step assembles the cotangent
// Laplacian of the *current* metric -- which is the exact Jacobian dK/du, see
// the header -- and solves L du = Kbar - K with one vertex of each connected
// component pinned to remove the constants. The step is
// halved until every face is still realisable, and the iteration stops when
// ||K - Kbar||_inf is under tolerance or a step can no longer be taken.
// ---------------------------------------------------------------------------
void ConeMetric::solve(const Mesh &mesh, const ConeSingularities &cones) {
    const int nV = static_cast<int>(mesh.vertices.size());
    const int nT = static_cast<int>(mesh.triangles.size());
    const std::vector<double> &target = cones.getTargetCurvature();
    const std::vector<char> &active = cones.getActiveVertices();
    if (static_cast<int>(target.size()) != nV || static_cast<int>(active.size()) != nV) {
        report.messages.push_back(
            "the cone set and the mesh disagree about how many vertices there are, so no cone "
            "metric was built and E1 falls back to the map's own induced lengths");
        return;
    }

    // The curvature of a metric given by edge lengths, and the angles it needs.
    std::vector<double> K(nV, 0.0);
    std::vector<double> angle(3 * nT, 0.0);
    auto measure = [&]() {
        for (int v = 0; v < nV; ++v) {
            K[v] = active[v] ? (mesh.isBoundaryVertex[v] ? M_PI : 2.0 * M_PI) : 0.0;
        }
        report.nonRealisable = 0;
        for (int t = 0; t < nT; ++t) {
            const std::array<int, 3> &te = mesh.triangleEdges[t];
            const double l[3] = {len[te[0]], len[te[1]], len[te[2]]};
            if (!realisable(l, opts.realisableMargin)) ++report.nonRealisable;
            double a[3];
            anglesOf(l, a);
            for (int i = 0; i < 3; ++i) {
                angle[3 * t + i] = a[i];
                K[mesh.triangles[t][i]] -= a[i];
            }
        }
    };
    auto error = [&]() {
        double worst = 0.0;
        for (int v = 0; v < nV; ++v) {
            if (!active[v]) continue;
            worst = std::max(worst, std::fabs(K[v] - target[v]));
        }
        return worst;
    };

    measure();
    report.initialError = error();

    // Eq. (4) as a solvability statement rather than as an admissibility one:
    // L du = Kbar - K has a solution only if the right-hand side is orthogonal
    // to the constants, which is exactly sum(Kbar) = sum(K) = 2 pi chi.
    {
        double s = 0.0;
        for (int v = 0; v < nV; ++v) if (active[v]) s += target[v] - K[v];
        report.gaussBonnetResidual = std::fabs(s);
        if (report.gaussBonnetResidual > 1e-6) {
            std::ostringstream oss;
            oss << "sum(Kbar) - sum(K) = " << report.gaussBonnetResidual
                << " rather than zero, so no metric conformal to the model carries this cone "
                << "set. That is Eq. (4) failing, and Stage 1 is where it is answered.";
            report.messages.push_back(oss.str());
            return;
        }
    }

    // The pins. Constants are a global scaling of the metric and change no
    // angle, so L has one null vector per connected component and each of them
    // has to be removed -- pinning a single vertex leaves the other components
    // singular, and a singular block does not announce itself: LDL^T returns a
    // number. So the components are found and one vertex of each is pinned.
    std::vector<char> pinned(nV, 0);
    {
        std::vector<int> comp(nV, -1), stack;
        int nComp = 0;
        for (int seed = 0; seed < nV; ++seed) {
            if (!active[seed] || comp[seed] >= 0) continue;
            comp[seed] = nComp;
            pinned[seed] = 1;
            stack.assign(1, seed);
            while (!stack.empty()) {
                const int v = stack.back();
                stack.pop_back();
                const auto &vt = mesh.vertexTriangles;
                for (int k = vt.rowPtr[v]; k < vt.rowPtr[v + 1]; ++k) {
                    const Triangle &tri = mesh.triangles[vt.colIdx[k]];
                    for (int i = 0; i < 3; ++i) {
                        const int w = tri[i];
                        if (!active[w] || comp[w] >= 0) continue;
                        comp[w] = nComp;
                        stack.push_back(w);
                    }
                }
            }
            ++nComp;
        }
        if (nComp == 0) return;
        if (nComp > 1) {
            std::ostringstream oss;
            oss << "the model has " << nComp << " connected components, so the cone metric is "
                << "pinned once in each. Eq. (4) is checked over the whole model and each "
                << "component needs it separately, which Stage 1 does not test.";
            report.messages.push_back(oss.str());
        }
    }

    std::vector<double> du(nV, 0.0), prev(nV, 0.0);
    std::vector<Eigen::Triplet<double>> trip;
    Eigen::VectorXd rhs(nV), sol(nV);

    for (int it = 0; it < opts.newtonSteps; ++it) {
        const double err = error();
        if (it == 1) report.linearError = err;
        if (err < opts.tolerance) { report.converged = true; break; }

        // dK = L du with w_ij = (1/2) sum_t cot(alpha_ij^t): the cotangent
        // Laplacian of the metric as it now stands.
        trip.clear();
        trip.reserve(nT * 12 + nV);
        std::vector<double> diag(nV, 0.0);
        for (int t = 0; t < nT; ++t) {
            const Triangle &tri = mesh.triangles[t];
            for (int i = 0; i < 3; ++i) {
                // The edge joining tri[i] and tri[i+1] is opposite corner i+2.
                const int a = tri[i], b = tri[(i + 1) % 3];
                if (!active[a] || !active[b]) continue;
                const double op = angle[3 * t + (i + 2) % 3];
                const double s = std::sin(op);
                if (!(std::fabs(s) > 1e-12)) continue;
                const double w = 0.5 * std::cos(op) / s;
                if (!std::isfinite(w)) continue;
                trip.emplace_back(a, b, -w);
                trip.emplace_back(b, a, -w);
                diag[a] += w;
                diag[b] += w;
            }
        }
        for (int v = 0; v < nV; ++v) {
            if (!active[v]) { trip.emplace_back(v, v, 1.0); continue; }
            trip.emplace_back(v, v, diag[v]);
        }
        // The pins, as rows and columns, so the matrix stays symmetric.
        std::vector<Eigen::Triplet<double>> kept;
        kept.reserve(trip.size() + nV);
        for (const auto &t : trip) {
            if (pinned[t.row()] || pinned[t.col()]) continue;
            kept.push_back(t);
        }
        for (int v = 0; v < nV; ++v) if (pinned[v]) kept.emplace_back(v, v, 1.0);

        Eigen::SparseMatrix<double> L(nV, nV);
        L.setFromTriplets(kept.begin(), kept.end());
        L.makeCompressed();

        for (int v = 0; v < nV; ++v) {
            rhs[v] = (active[v] && !pinned[v]) ? (target[v] - K[v]) : 0.0;
        }

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
        ldlt.compute(L);
        if (ldlt.info() != Eigen::Success) {
            report.messages.push_back(
                "the Laplacian of the current metric would not factor, so the cone metric "
                "stopped at the iterate it had");
            break;
        }
        sol = ldlt.solve(rhs);
        if (ldlt.info() != Eigen::Success) break;

        bool finite = true;
        for (int v = 0; v < nV; ++v) {
            du[v] = active[v] ? sol[v] : 0.0;
            if (!std::isfinite(du[v])) finite = false;
        }
        if (!finite) {
            report.messages.push_back("the Newton step came back non-finite; the cone metric "
                                      "stopped at the iterate it had");
            break;
        }

        // Halve until the whole metric is realisable *and* the error has not
        // grown. A full step is what is wanted and what is taken on almost
        // every model; the halving is for the cone strong enough to pull a
        // triangle through its own opposite side.
        prev = u;
        double frac = 1.0;
        bool took = false;
        while (frac >= opts.minStepFraction) {
            for (int v = 0; v < nV; ++v) u[v] = prev[v] + frac * du[v];
            rebuild(mesh);
            measure();
            if (report.nonRealisable == 0 && error() < err) { took = true; break; }
            frac *= 0.5;
            ++report.halvings;
        }
        if (!took) {
            u = prev;
            rebuild(mesh);
            measure();
            std::ostringstream oss;
            oss << "the Newton step could not be shortened into the realisable set at "
                << "||K - Kbar||_inf = " << err << " rad; the metric stands at that iterate";
            report.messages.push_back(oss.str());
            break;
        }
        ++report.newtonIterations;
    }

    measure();
    report.finalError = error();
    if (report.newtonIterations <= 1) report.linearError = report.finalError;
    if (report.finalError < opts.tolerance) report.converged = true;
    report.solved = report.nonRealisable == 0 && std::isfinite(report.finalError);

    // C4 in the number Stage 6 reports: the worst cone whose angle sum here is
    // not the 2 pi - (pi/2) I the layout is going to be held to.
    const std::vector<int> &idx = cones.getIndices();
    for (const auto &c : cones.getCones()) {
        const int v = c.vertex;
        if (v < 0 || v >= nV || !active[v]) continue;
        const double full = mesh.isBoundaryVertex[v] ? M_PI : 2.0 * M_PI;
        const double want = full - M_PI_2 * idx[v];
        report.coneResidual = std::max(report.coneResidual, std::fabs((full - K[v]) - want));
    }

    if (!report.converged && report.solved) {
        std::ostringstream oss;
        oss << "the cone metric stopped at ||K - Kbar||_inf = " << report.finalError
            << " rad after " << report.newtonIterations << " Newton step(s), from "
            << report.initialError << ". E1 measures against it as it stands, and C4 is the "
            << "number that says whether that is good enough.";
        report.messages.push_back(oss.str());
    }
    if (report.nonRealisable > 0) {
        std::ostringstream oss;
        oss << report.nonRealisable << " face(s) of the cone metric fail the triangle "
            << "inequality, so it is not a metric at all and E1 falls back to the map's own "
            << "induced lengths. A cone strong enough to do that wants the flow's edge flips, "
            << "which is what --ref ricci is.";
        report.messages.push_back(oss.str());
    }
}
