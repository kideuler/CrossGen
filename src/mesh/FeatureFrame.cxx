#include "mesh/FeatureFrame.hxx"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>

#include <algorithm>
#include <cmath>
#include <limits>

std::complex<double> FeatureFrame::crossOf(const Point &d) {
    const double n2 = dotP(d, d);
    if (!(n2 > 0.0)) return {0.0, 0.0};
    const std::complex<double> z(d[0], d[1]);
    const std::complex<double> z2 = z * z;
    return z2 * z2 / (n2 * n2);
}

FeatureFrame::FeatureFrame(const Mesh &model) : FeatureFrame(model, Options()) {}

FeatureFrame::FeatureFrame(const Mesh &model, const Options &opts)
    : opts_(opts), vertices_(model.vertices), triangles_(model.triangles) {
    const int nV = static_cast<int>(vertices_.size());
    const int nT = static_cast<int>(triangles_.size());
    report_.vertices = nV;
    report_.triangles = nT;
    u_.assign(nV, {0.0, 0.0});
    if (nT == 0) {
        report_.reason = "the model has no triangles";
        return;
    }

    // A vertex no triangle uses has no equation; mesh::Mesh keeps every vertex
    // of its file, and the selection dataset's meshes have had such strays.
    std::vector<char> used(nV, 0);
    for (const Triangle &t : triangles_)
        for (int v : t) used[v] = 1;

    // Each feature edge's cross, summed onto both its ends.
    std::vector<std::complex<double>> fixedValue(nV, {0.0, 0.0});
    std::vector<int> featureEdges(nV, 0);
    for (size_t e = 0; e < model.edges.size(); ++e) {
        const bool boundary = e < model.isBoundaryEdge.size() && model.isBoundaryEdge[e];
        bool interface = false;
        if (!boundary && opts_.interfaces && e < model.edgeTriangles.size()) {
            const std::array<int, 2> &et = model.edgeTriangles[e];
            interface = et[0] >= 0 && et[1] >= 0 &&
                        model.triangleMatId[et[0]] != model.triangleMatId[et[1]];
        }
        if (!boundary && !interface) continue;
        const int a = model.edges[e][0], b = model.edges[e][1];
        const std::complex<double> c = crossOf(vertices_[b] - vertices_[a]);
        fixedValue[a] += c;
        fixedValue[b] += c;
        ++featureEdges[a];
        ++featureEdges[b];
        ++(boundary ? report_.boundaryEdges : report_.interfaceEdges);
    }

    std::vector<int> unknown(nV, -1);
    int n = 0;
    for (int v = 0; v < nV; ++v) {
        if (featureEdges[v] > 0) {
            u_[v] = fixedValue[v] / static_cast<double>(featureEdges[v]);
            ++report_.fixedVertices;
        } else if (used[v]) {
            unknown[v] = n++;
        }
    }

    // The cotangent Laplacian, assembled on the unknowns, with the fixed
    // vertices' values moved to the right-hand side.
    std::vector<Eigen::Triplet<double>> K;
    K.reserve(9 * static_cast<size_t>(nT));
    Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(n, 2);
    for (const Triangle &t : triangles_) {
        for (int k = 0; k < 3; ++k) {
            const int a = t[k], b = t[(k + 1) % 3], c = t[(k + 2) % 3];
            const Point ea = vertices_[a] - vertices_[c];
            const Point eb = vertices_[b] - vertices_[c];
            const double twiceArea = std::fabs(cross2(ea, eb));
            if (!(twiceArea > 0.0)) continue;
            // Half the cotangent of the angle at c, the weight of edge (a, b).
            const double w = 0.5 * dotP(ea, eb) / twiceArea;
            const int ia = unknown[a], ib = unknown[b];
            if (ia >= 0) K.emplace_back(ia, ia, w);
            if (ib >= 0) K.emplace_back(ib, ib, w);
            if (ia >= 0 && ib >= 0) {
                K.emplace_back(ia, ib, -w);
                K.emplace_back(ib, ia, -w);
            } else if (ia >= 0) {
                rhs(ia, 0) += w * u_[b].real();
                rhs(ia, 1) += w * u_[b].imag();
            } else if (ib >= 0) {
                rhs(ib, 0) += w * u_[a].real();
                rhs(ib, 1) += w * u_[a].imag();
            }
        }
    }

    if (n > 0) {
        Eigen::SparseMatrix<double> A(n, n);
        A.setFromTriplets(K.begin(), K.end());
        Eigen::MatrixXd x;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(A);
        if (ldlt.info() == Eigen::Success) x = ldlt.solve(rhs);
        if (ldlt.info() != Eigen::Success || !x.allFinite()) {
            // A stiffness matrix with Dirichlet rows is positive definite, so
            // this is for a mesh degenerate enough to defeat the factorisation.
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg(A);
            cg.setTolerance(1e-10);
            x = cg.solve(rhs);
            if (cg.info() != Eigen::Success || !x.allFinite()) {
                report_.reason = "the Laplace solve did not converge";
                x = Eigen::MatrixXd::Zero(n, 2);
            }
        }
        for (int v = 0; v < nV; ++v)
            if (unknown[v] >= 0) u_[v] = {x(unknown[v], 0), x(unknown[v], 1)};
    }
    report_.solved = report_.reason.empty();

    double area = 0.0, weighted = 0.0;
    for (const Triangle &t : triangles_) {
        const double a = 0.5 * std::fabs(cross2(vertices_[t[1]] - vertices_[t[0]],
                                                vertices_[t[2]] - vertices_[t[0]]));
        area += a;
        weighted += a * std::abs((u_[t[0]] + u_[t[1]] + u_[t[2]]) / 3.0);
    }
    report_.meanCoherence = area > 0.0 ? weighted / area : 0.0;

    buildBuckets();
}

void FeatureFrame::buildBuckets() {
    const int nT = static_cast<int>(triangles_.size());
    Point hi{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    lo_ = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    for (const Triangle &t : triangles_) {
        for (int v : t) {
            lo_ = {std::min(lo_[0], vertices_[v][0]), std::min(lo_[1], vertices_[v][1])};
            hi = {std::max(hi[0], vertices_[v][0]), std::max(hi[1], vertices_[v][1])};
        }
    }
    const double w = std::max(hi[0] - lo_[0], 1e-300);
    const double h = std::max(hi[1] - lo_[1], 1e-300);
    const double cell = std::sqrt(2.0 * w * h / std::max(nT, 1));
    nx_ = std::min(std::max(static_cast<int>(std::ceil(w / cell)), 1), 1024);
    ny_ = std::min(std::max(static_cast<int>(std::ceil(h / cell)), 1), 1024);
    cellW_ = w / nx_;
    cellH_ = h / ny_;

    auto range = [&](const Triangle &t, int &x0, int &x1, int &y0, int &y1) {
        double a0 = vertices_[t[0]][0], a1 = a0, b0 = vertices_[t[0]][1], b1 = b0;
        for (int k = 1; k < 3; ++k) {
            a0 = std::min(a0, vertices_[t[k]][0]);
            a1 = std::max(a1, vertices_[t[k]][0]);
            b0 = std::min(b0, vertices_[t[k]][1]);
            b1 = std::max(b1, vertices_[t[k]][1]);
        }
        x0 = std::min(std::max(static_cast<int>((a0 - lo_[0]) / cellW_), 0), nx_ - 1);
        x1 = std::min(std::max(static_cast<int>((a1 - lo_[0]) / cellW_), 0), nx_ - 1);
        y0 = std::min(std::max(static_cast<int>((b0 - lo_[1]) / cellH_), 0), ny_ - 1);
        y1 = std::min(std::max(static_cast<int>((b1 - lo_[1]) / cellH_), 0), ny_ - 1);
    };
    bucketStart_.assign(static_cast<size_t>(nx_) * ny_ + 1, 0);
    for (const Triangle &t : triangles_) {
        int x0, x1, y0, y1;
        range(t, x0, x1, y0, y1);
        for (int y = y0; y <= y1; ++y)
            for (int x = x0; x <= x1; ++x) ++bucketStart_[static_cast<size_t>(y) * nx_ + x + 1];
    }
    for (size_t i = 1; i < bucketStart_.size(); ++i) bucketStart_[i] += bucketStart_[i - 1];
    bucketTriangles_.assign(bucketStart_.back(), -1);
    std::vector<int> fill(bucketStart_.begin(), bucketStart_.end() - 1);
    for (int i = 0; i < nT; ++i) {
        int x0, x1, y0, y1;
        range(triangles_[i], x0, x1, y0, y1);
        for (int y = y0; y <= y1; ++y)
            for (int x = x0; x <= x1; ++x) bucketTriangles_[fill[static_cast<size_t>(y) * nx_ + x]++] = i;
    }
}

std::complex<double> FeatureFrame::at(const Point &p) const {
    if (triangles_.empty() || nx_ == 0) return {0.0, 0.0};
    const int cx = std::min(std::max(static_cast<int>(std::floor((p[0] - lo_[0]) / cellW_)), 0), nx_ - 1);
    const int cy = std::min(std::max(static_cast<int>(std::floor((p[1] - lo_[1]) / cellH_)), 0), ny_ - 1);

    // The triangle whose smallest barycentric coordinate at p is largest:
    // the one containing p when one does, and otherwise one p is just outside
    // of. A triangle containing p is in p's own bucket, so ring 0 settles
    // every point of the model; ring 1 and on serve the points outside it.
    int best = -1;
    double bestMin = -std::numeric_limits<double>::infinity();
    std::array<double, 3> bestL{0.0, 0.0, 0.0};
    const int rings = std::max(nx_, ny_);
    for (int r = 0; r <= rings; ++r) {
        for (int y = cy - r; y <= cy + r; ++y) {
            if (y < 0 || y >= ny_) continue;
            for (int x = cx - r; x <= cx + r; ++x) {
                if (x < 0 || x >= nx_) continue;
                if (std::max(std::abs(x - cx), std::abs(y - cy)) != r) continue;   // the ring only
                const size_t bucket = static_cast<size_t>(y) * nx_ + x;
                for (int k = bucketStart_[bucket]; k < bucketStart_[bucket + 1]; ++k) {
                    const Triangle &t = triangles_[bucketTriangles_[k]];
                    const Point &a = vertices_[t[0]];
                    const double det = cross2(vertices_[t[1]] - a, vertices_[t[2]] - a);
                    if (std::fabs(det) < 1e-300) continue;
                    const double l1 = cross2(p - a, vertices_[t[2]] - a) / det;
                    const double l2 = cross2(vertices_[t[1]] - a, p - a) / det;
                    const double l0 = 1.0 - l1 - l2;
                    const double m = std::min(l0, std::min(l1, l2));
                    if (m > bestMin) {
                        bestMin = m;
                        best = bucketTriangles_[k];
                        bestL = {l0, l1, l2};
                    }
                }
            }
        }
        if (best >= 0 && (bestMin >= -1e-12 || r >= 1)) break;
    }
    if (best < 0) return {0.0, 0.0};

    // Outside: clamp to the triangle, which leaves a point on its boundary.
    double s = 0.0;
    for (double &l : bestL) s += (l = std::max(l, 0.0));
    if (!(s > 0.0)) return {0.0, 0.0};
    const Triangle &t = triangles_[best];
    return (bestL[0] * u_[t[0]] + bestL[1] * u_[t[1]] + bestL[2] * u_[t[2]]) / s;
}

void FeatureFrame::Alignment::add(const Point &centre, const Point &du, const Point &dv, double area) {
    if (!(area > 0.0)) return;
    const std::complex<double> u = frame_.at(centre);
    const double rho = std::abs(u);
    if (!(rho > 0.0)) return;
    // sin^2 2 Delta = (1 - cos 4 Delta) / 2, cos 4 Delta = Re(c conj(u)) / |u|.
    const double off = 0.5 * (1.0 - std::real(crossOf(du) * std::conj(u)) / rho) +
                       0.5 * (1.0 - std::real(crossOf(dv) * std::conj(u)) / rho);
    num_ += area * rho * 0.5 * off;
    den_ += area * rho;
}

double FeatureFrame::Alignment::quality() const {
    if (!(den_ > 0.0)) return 1.0;
    const double E = std::min(std::max(num_ / den_, 0.0), 1.0);
    return 1.0 - (2.0 / M_PI) * std::asin(std::sqrt(E));
}
