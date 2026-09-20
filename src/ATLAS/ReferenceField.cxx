#include "ATLAS/ReferenceField.hxx"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <sstream>

namespace {

typedef std::complex<double> cplx;

// Minimum-cost assignment of n rows to distinct columns among m >= n
// (Kuhn-Munkres with potentials, O(n^2 m)). Returns the column of each row.
std::vector<int> assignRows(const std::vector<std::vector<double>> &a, int n, int m) {
    const double INF = std::numeric_limits<double>::infinity();
    std::vector<double> u(n + 1, 0.0), v(m + 1, 0.0);
    std::vector<int> p(m + 1, 0), way(m + 1, 0);
    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;
        std::vector<double> minv(m + 1, INF);
        std::vector<char> used(m + 1, 0);
        do {
            used[j0] = 1;
            const int i0 = p[j0];
            int j1 = 0;
            double delta = INF;
            for (int j = 1; j <= m; ++j) {
                if (used[j]) continue;
                const double cur = a[i0 - 1][j - 1] - u[i0] - v[j];
                if (cur < minv[j]) { minv[j] = cur; way[j] = j0; }
                if (minv[j] < delta) { delta = minv[j]; j1 = j; }
            }
            for (int j = 0; j <= m; ++j) {
                if (used[j]) { u[p[j]] += delta; v[j] -= delta; }
                else minv[j] -= delta;
            }
            j0 = j1;
        } while (p[j0] != 0);
        do {
            const int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0);
    }
    std::vector<int> col(n, -1);
    for (int j = 1; j <= m; ++j) if (p[j] > 0) col[p[j] - 1] = j - 1;
    return col;
}

bool insidePolygon(const std::vector<Point> &poly, const Point &x) {
    bool in = false;
    const size_t n = poly.size();
    for (size_t i = 0, j = n - 1; i < n; j = i++) {
        const Point &a = poly[i], &b = poly[j];
        if ((a[1] > x[1]) != (b[1] > x[1])) {
            const double xc = a[0] + (x[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1]);
            if (x[0] < xc) in = !in;
        }
    }
    return in;
}

} // namespace

ReferenceField::ReferenceField(std::shared_ptr<Mesh> mesh, const std::vector<int> &interfaceEdges,
                               const Options &opts)
    : mesh_(std::move(mesh)), opts_(opts) {
    const auto t0 = std::chrono::steady_clock::now();
    if (!mesh_ || mesh_->triangles.empty() || mesh_->vertices.empty()) {
        report_.reason = "no mesh";
        return;
    }
    report_.triangles = static_cast<int>(mesh_->triangles.size());
    Point lo = mesh_->vertices[0], hi = lo;
    for (const Point &p : mesh_->vertices) {
        lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
        hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
    }
    diag_ = std::hypot(hi[0] - lo[0], hi[1] - lo[1]);
    area_ = 0.0;
    for (const Triangle &t : mesh_->triangles) {
        area_ += 0.5 * std::fabs(cross2(mesh_->vertices[t[1]] - mesh_->vertices[t[0]],
                                        mesh_->vertices[t[2]] - mesh_->vertices[t[0]]));
    }
    if (!(diag_ > 0.0) || !(area_ > 0.0)) {
        report_.reason = "the mesh has no extent";
        return;
    }
    solve(interfaceEdges);
    buildGrid();
    report_.built = true;
    report_.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// ---------------------------------------------------------------------------
// The solve: TORSION::runField's tau-continuation, unchanged.
// ---------------------------------------------------------------------------
void ReferenceField::solve(const std::vector<int> &interfaceEdges) {
    // One level: a fresh operator at this tau with the same Dirichlet data,
    // started from the free triangles of the level before.
    auto buildLevel = [&](double tauScale, const Eigen::VectorXcd &carried) {
        auto f = std::make_unique<DualMBO>(mesh_, opts_.maxSteps, opts_.gamma);
        f->setPenaltyWeight(opts_.weight);
        f->setTauScale(tauScale);
        f->setPinDiskCenters(opts_.pinDiskCenters);
        if (opts_.alignToInterfaces && !interfaceEdges.empty()) f->setAlignedInteriorEdges(interfaceEdges);
        f->initialize();
        if (carried.size() == f->u_k_prev.size()) {
            Eigen::VectorXcd start = carried;
            for (const auto &[ti, bc] : f->getBoundaryData()) start[ti] = bc;
            f->u_k_prev = start;
            f->u_k = start;
        }
        return f;
    };
    std::unique_ptr<DualMBO> field = buildLevel(1.0, Eigen::VectorXcd());

    // Largest tau first, down to the floor read off the assembled operator.
    std::vector<double> ladder{1.0};
    if (opts_.tauContinuation && opts_.tauRatio > 0.0 && opts_.tauRatio < 1.0) {
        const double rate = field->medianDiffusionRate();
        if (rate > 0.0) {
            const double tau0 = diag_ * diag_ / 10.0;
            const double c = opts_.tauFloorEdges;
            const double tauMin = c * c / rate;
            double tau = tau0;
            while (tau * opts_.tauRatio > tauMin && ladder.size() < 40) {
                tau *= opts_.tauRatio;
                ladder.push_back(ladder.back() * opts_.tauRatio);
            }
        }
    }
    const double nTris = static_cast<double>(mesh_->triangles.size());
    Eigen::VectorXcd carried;
    for (std::size_t level = 0; level < ladder.size(); ++level) {
        if (level > 0) field = buildLevel(ladder[level], carried);
        const int cap = ladder.size() == 1 ? opts_.maxSteps : std::min(opts_.maxSteps, opts_.tauLevelSteps);
        report_.converged = false;
        for (int i = 0; i < cap; ++i) {
            field->step();
            ++report_.steps;
            if (field->error < 2.0 * nTris * 1e-5) { report_.converged = true; break; }
        }
        carried = field->u_k_prev;
    }
    report_.levels = static_cast<int>(ladder.size());
    if (!report_.converged) {
        std::ostringstream os;
        os << "the reference field did not converge in " << opts_.maxSteps << " MBO steps at its last level";
        report_.messages.push_back(os.str());
    }
    field->computeSingularities();
    u_ = field->u_k_prev;
    for (Eigen::Index t = 0; t < u_.size(); ++t) {
        const double m = std::abs(u_[t]);
        u_[t] = m > 0.0 ? u_[t] / m : cplx(0.0, 0.0);
    }
    findCones(*field);
}

void ReferenceField::findCones(const DualMBO &f) {
    const Mesh &M = *mesh_;
    cones_.clear();
    for (const auto &[v, idx] : f.singularVertices) {
        if (M.isBoundaryVertex[v]) continue;
        const int units = static_cast<int>(std::lround(std::fabs(idx) * 4.0));
        if (units == 0) continue;
        // The material component the cone sits in; -1 on an interface, where
        // it is no one region's to cancel.
        int region = -2;
        for (int k = M.vertexTriangles.rowPtr[v]; k < M.vertexTriangles.rowPtr[v + 1]; ++k) {
            const int t = M.vertexTriangles.colIdx[k];
            const int c = t < static_cast<int>(M.triangleComponent.size()) ? M.triangleComponent[t] : 0;
            if (region == -2) region = c;
            else if (region != c) region = -1;
        }
        Singularity s;
        s.x = M.vertices[v];
        s.sign = idx > 0.0 ? 1 : -1;
        s.vertex = v;
        s.region = region < 0 ? -1 : region;
        for (int k = 0; k < units; ++k) cones_.push_back(s);
        (s.sign > 0 ? report_.rawPlus : report_.rawMinus) += units;
    }
    if (opts_.cancelDipoles) {
        const double radius = opts_.dipoleRadius > 0.0 ? opts_.dipoleRadius * diag_
                                                       : std::numeric_limits<double>::infinity();
        std::vector<char> alive(cones_.size(), 1);
        while (true) {
            int bi = -1, bj = -1;
            double bd = radius;
            for (size_t i = 0; i < cones_.size(); ++i) {
                if (!alive[i] || cones_[i].region < 0) continue;
                for (size_t j = i + 1; j < cones_.size(); ++j) {
                    if (!alive[j] || cones_[j].region != cones_[i].region) continue;
                    if (cones_[j].sign == cones_[i].sign) continue;
                    const double d = normP(cones_[i].x - cones_[j].x);
                    if (d < bd) { bd = d; bi = static_cast<int>(i); bj = static_cast<int>(j); }
                }
            }
            if (bi < 0) break;
            alive[bi] = alive[bj] = 0;
            report_.dipoleUnits += 1;
        }
        std::vector<Singularity> kept;
        for (size_t i = 0; i < cones_.size(); ++i) if (alive[i]) kept.push_back(cones_[i]);
        cones_.swap(kept);
    }
    for (const Singularity &s : cones_) (s.sign > 0 ? report_.conesPlus : report_.conesMinus)++;
}

// ---------------------------------------------------------------------------
// The lookup grid: u averaged to the mesh vertices (so it falls to zero at a
// cone, whose star it winds round), linear on each triangle, sampled at the
// grid nodes; nodes outside the domain take their known neighbours' mean.
// ---------------------------------------------------------------------------
void ReferenceField::buildGrid() {
    const Mesh &M = *mesh_;
    const int NV = static_cast<int>(M.vertices.size());
    std::vector<cplx> uv(NV, cplx(0.0, 0.0));
    std::vector<double> wv(NV, 0.0);
    double edgeSum = 0.0;
    int edgeCount = 0;
    for (size_t t = 0; t < M.triangles.size(); ++t) {
        const Triangle &T = M.triangles[t];
        const double A = 0.5 * std::fabs(cross2(M.vertices[T[1]] - M.vertices[T[0]], M.vertices[T[2]] - M.vertices[T[0]]));
        for (int k = 0; k < 3; ++k) {
            uv[T[k]] += A * u_[static_cast<Eigen::Index>(t)];
            wv[T[k]] += A;
            edgeSum += normP(M.vertices[T[(k + 1) % 3]] - M.vertices[T[k]]);
            ++edgeCount;
        }
    }
    for (int v = 0; v < NV; ++v) if (wv[v] > 0.0) uv[v] /= wv[v];

    Point lo = M.vertices[0], hi = lo;
    for (const Point &p : M.vertices) {
        lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
        hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
    }
    const double meanEdge = edgeCount > 0 ? edgeSum / edgeCount : diag_ / 100.0;
    step_ = std::max(opts_.gridStep * meanEdge, 1e-12 * diag_);
    const double W = hi[0] - lo[0], H = hi[1] - lo[1];
    while ((W / step_ + 5.0) * (H / step_ + 5.0) > static_cast<double>(opts_.maxGridNodes)) step_ *= 1.25;
    const int pad = 2;
    lo_ = {lo[0] - pad * step_, lo[1] - pad * step_};
    nx_ = static_cast<int>(std::ceil(W / step_)) + 2 * pad + 1;
    ny_ = static_cast<int>(std::ceil(H / step_)) + 2 * pad + 1;
    report_.gridNx = nx_;
    report_.gridNy = ny_;
    grid_.assign(static_cast<size_t>(nx_) * ny_, cplx(0.0, 0.0));
    std::vector<char> known(grid_.size(), 0);
    for (const Triangle &t : M.triangles) {
        const Point &a = M.vertices[t[0]], &b = M.vertices[t[1]], &c = M.vertices[t[2]];
        const double det = cross2(b - a, c - a);
        if (!(std::fabs(det) > 0.0)) continue;
        const int i0 = std::max(0, static_cast<int>(std::floor((std::min({a[0], b[0], c[0]}) - lo_[0]) / step_)));
        const int i1 = std::min(nx_ - 1, static_cast<int>(std::ceil((std::max({a[0], b[0], c[0]}) - lo_[0]) / step_)));
        const int j0 = std::max(0, static_cast<int>(std::floor((std::min({a[1], b[1], c[1]}) - lo_[1]) / step_)));
        const int j1 = std::min(ny_ - 1, static_cast<int>(std::ceil((std::max({a[1], b[1], c[1]}) - lo_[1]) / step_)));
        for (int j = j0; j <= j1; ++j) {
            for (int i = i0; i <= i1; ++i) {
                const Point p{lo_[0] + i * step_, lo_[1] + j * step_};
                const double l1 = cross2(c - b, p - b) / det, l2 = cross2(a - c, p - c) / det;
                const double l3 = 1.0 - l1 - l2;
                if (l1 < -1e-12 || l2 < -1e-12 || l3 < -1e-12) continue;
                const size_t k = static_cast<size_t>(i) + static_cast<size_t>(nx_) * j;
                grid_[k] = l1 * uv[t[0]] + l2 * uv[t[1]] + l3 * uv[t[2]];
                known[k] = 1;
            }
        }
    }
    // A coarse chord across a concave arc can stray a few percent of the
    // diagonal outside; fill that far and a little more.
    const int sweeps = std::max(8, static_cast<int>(std::ceil(0.06 * diag_ / step_)));
    for (int sweep = 0; sweep < sweeps; ++sweep) {
        std::vector<char> next = known;
        bool changed = false;
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                const size_t k = static_cast<size_t>(i) + static_cast<size_t>(nx_) * j;
                if (known[k]) continue;
                cplx s(0.0, 0.0);
                int n = 0;
                for (int dj = -1; dj <= 1; ++dj) {
                    for (int di = -1; di <= 1; ++di) {
                        const int ii = i + di, jj = j + dj;
                        if (ii < 0 || jj < 0 || ii >= nx_ || jj >= ny_) continue;
                        const size_t kk = static_cast<size_t>(ii) + static_cast<size_t>(nx_) * jj;
                        if (!known[kk]) continue;
                        s += grid_[kk];
                        ++n;
                    }
                }
                if (n > 0) {
                    grid_[k] = s / static_cast<double>(n);
                    next[k] = 1;
                    changed = true;
                }
            }
        }
        known.swap(next);
        if (!changed) break;
    }
}

std::complex<double> ReferenceField::valueAt(const Point &p) const {
    if (grid_.empty()) return cplx(0.0, 0.0);
    const double x = (p[0] - lo_[0]) / step_, y = (p[1] - lo_[1]) / step_;
    const int i = std::max(0, std::min(nx_ - 2, static_cast<int>(std::floor(x))));
    const int j = std::max(0, std::min(ny_ - 2, static_cast<int>(std::floor(y))));
    const double fx = std::max(0.0, std::min(1.0, x - i)), fy = std::max(0.0, std::min(1.0, y - j));
    auto at = [&](int a, int b) { return grid_[static_cast<size_t>(a) + static_cast<size_t>(nx_) * b]; };
    return (1 - fx) * (1 - fy) * at(i, j) + fx * (1 - fy) * at(i + 1, j) + (1 - fx) * fy * at(i, j + 1) +
           fx * fy * at(i + 1, j + 1);
}

double ReferenceField::crossAt(const Point &p, double *coherence) const {
    const cplx u = valueAt(p);
    if (coherence) *coherence = std::min(1.0, std::abs(u));
    double th = std::arg(u) / 4.0;
    th = std::fmod(th, M_PI_2);
    if (th < 0.0) th += M_PI_2;
    return th;
}

std::complex<double> ReferenceField::average(const std::array<Point, 4> &q) const {
    static const double s3[3] = {1.0 / 6.0, 0.5, 5.0 / 6.0};
    cplx sum(0.0, 0.0), plain(0.0, 0.0);
    double W = 0.0;
    for (double s : s3) {
        for (double t : s3) {
            const Point x = q[0] * ((1 - s) * (1 - t)) + q[1] * (s * (1 - t)) + q[2] * (s * t) + q[3] * ((1 - s) * t);
            const Point xs = (q[1] - q[0]) * (1 - t) + (q[2] - q[3]) * t;
            const Point xt = (q[3] - q[0]) * (1 - s) + (q[2] - q[1]) * s;
            const double w = std::max(0.0, cross2(xs, xt));
            const cplx u = valueAt(x);
            sum += w * u;
            plain += u;
            W += w;
        }
    }
    return W > 0.0 ? sum / W : plain / 9.0;
}

std::complex<double> ReferenceField::cellCross(const std::array<Point, 4> &q) {
    const Point du = (q[1] - q[0]) + (q[2] - q[3]), dv = (q[3] - q[0]) + (q[2] - q[1]);
    const double a = 4.0 * std::atan2(du[1], du[0]);
    const double b = 4.0 * (std::atan2(dv[1], dv[0]) - M_PI_2);
    const cplx c = std::polar(1.0, a) + std::polar(1.0, b);
    const double m = std::abs(c);
    return m > 1e-12 ? c / m : cplx(0.0, 0.0);
}

double ReferenceField::quadArea(const std::array<Point, 4> &q) {
    return 0.5 * cross2(q[2] - q[0], q[3] - q[1]);
}

double ReferenceField::misalignment(const std::array<Point, 4> &q) const {
    const cplx ub = average(q);
    const double rho = std::min(1.0, std::abs(ub));
    const cplx c = cellCross(q);
    const double m = 0.5 * (rho - (std::conj(c) * ub).real());
    return std::fabs(quadArea(q)) * std::max(0.0, m);
}

double ReferenceField::misalignment(const SquareCarrier &C, int q) const {
    const auto &c = C.cells[q];
    return misalignment({C.vertices[c[0]], C.vertices[c[1]], C.vertices[c[2]], C.vertices[c[3]]});
}

double ReferenceField::directionEnergy(const SquareCarrier &C) const {
    double s = 0.0;
    for (int q = 0; q < C.numCells(); ++q) s += misalignment(C, q);
    return s / area_;
}

std::vector<ReferenceField::Singularity> ReferenceField::conesIn(const std::vector<Point> &outer,
                                                                 const std::vector<std::vector<Point>> &holes) const {
    std::vector<Singularity> out;
    for (const Singularity &s : cones_) {
        if (!insidePolygon(outer, s.x)) continue;
        bool inHole = false;
        for (const auto &h : holes) if (insidePolygon(h, s.x)) { inHole = true; break; }
        if (!inHole) out.push_back(s);
    }
    return out;
}

std::vector<ReferenceField::Singularity> ReferenceField::carrierSingularities(const SquareCarrier &C,
                                                                              const std::vector<char> *only) {
    std::vector<Singularity> S;
    for (int v = 0; v < C.numVertices(); ++v) {
        if (C.valence[v] == 0) continue;
        if (only && !(*only)[v]) continue;
        const bool bd = C.boundaryVertex[v] != 0;
        // Index (target - valence) / 4: too few cells is +, too many -.
        const int d = (bd ? C.targetValence(v) : 4) - C.valence[v];
        if (d == 0) continue;
        Singularity s;
        s.x = C.vertices[v];
        s.sign = d > 0 ? 1 : -1;
        s.boundary = bd;
        s.vertex = v;
        for (int k = 0; k < std::abs(d); ++k) S.push_back(s);
    }
    return S;
}

// ---------------------------------------------------------------------------
// E_sing: every boundary defect pays its flat cost and may take a same-signed
// cone off the table (it is where the cone went); every interior singularity
// is matched to a same-signed cone at wPos min(1, d / r0), or pays wExtra;
// every cone left over pays wMissing. Opposite-signed interior pairs of the
// carrier within r0 -- a dislocation the field has no counterpart for -- are
// cancelled first, closest first, at wPos d / r0.
// ---------------------------------------------------------------------------
double ReferenceField::singularityEnergy(const std::vector<Singularity> &S, const std::vector<Singularity> &F,
                                         const SingularityWeights &w, SingularityReport *rep) const {
    SingularityReport R;
    R.cones = static_cast<int>(F.size());
    const double r0 = std::max(1e-300, w.r0 * diag_);
    std::vector<int> I, B;
    for (int i = 0; i < static_cast<int>(S.size()); ++i) {
        if (S[i].boundary) {
            B.push_back(i);
            (S[i].sign > 0 ? R.boundaryPlus : R.boundaryMinus)++;
        } else {
            I.push_back(i);
            (S[i].sign > 0 ? R.interiorPlus : R.interiorMinus)++;
        }
    }
    double E = 0.0;
    for (int b : B) E += S[b].sign > 0 ? w.wEdgePlus : w.wEdgeMinus;

    // The carrier's own dipoles.
    {
        struct Pair { double d; int a, b; };
        std::vector<Pair> pairs;
        for (size_t x = 0; x < I.size(); ++x) {
            for (size_t y = x + 1; y < I.size(); ++y) {
                const Singularity &a = S[I[x]], &b = S[I[y]];
                if (a.sign == b.sign) continue;
                const double d = normP(a.x - b.x);
                if (d < r0) pairs.push_back({d, static_cast<int>(x), static_cast<int>(y)});
            }
        }
        std::sort(pairs.begin(), pairs.end(), [](const Pair &p, const Pair &q) { return p.d < q.d; });
        std::vector<char> gone(I.size(), 0);
        for (const Pair &p : pairs) {
            if (gone[p.a] || gone[p.b]) continue;
            gone[p.a] = gone[p.b] = 1;
            E += w.wPos * p.d / r0;
            ++R.dipoles;
        }
        std::vector<int> keep;
        for (size_t x = 0; x < I.size(); ++x) if (!gone[x]) keep.push_back(I[x]);
        I.swap(keep);
    }

    // Rows: interior singularities, then boundary defects. Columns: the
    // cones, then one "unmatched" column per row. Costs are relative to the
    // cone being missed, so wMissing |F| is added back after.
    const int nI = static_cast<int>(I.size()), nB = static_cast<int>(B.size());
    const int n = nI + nB, nF = static_cast<int>(F.size());
    E += w.wMissing * nF;
    const double BIG = 1e6;
    auto rowCost = [&](int r, int j) {
        const Singularity &s = S[r < nI ? I[r] : B[r - nI]];
        if (j >= nF) return r < nI ? w.wExtra : 0.0;
        if (F[j].sign != s.sign) return BIG;
        if (r >= nI) return -w.wMissing;
        return w.wPos * std::min(1.0, normP(s.x - F[j].x) / r0) - w.wMissing;
    };
    std::vector<int> col(n, -1);
    if (n > 0 && nF > 0) {
        if (n <= 200 && nF <= 400) {
            const int m = nF + n;
            std::vector<std::vector<double>> a(n, std::vector<double>(m));
            for (int r = 0; r < n; ++r) for (int j = 0; j < m; ++j) a[r][j] = rowCost(r, j);
            col = assignRows(a, n, m);
        } else {
            // Too many for the exact assignment (a fine carrier): greedy by
            // saving, which is exact whenever no two rows compete.
            struct Opt { double saving; int r, j; };
            std::vector<Opt> opts;
            for (int r = 0; r < n; ++r) {
                const double idle = rowCost(r, nF);
                for (int j = 0; j < nF; ++j) {
                    const double c = rowCost(r, j);
                    if (c < idle) opts.push_back({idle - c, r, j});
                }
            }
            std::sort(opts.begin(), opts.end(), [](const Opt &x, const Opt &y) { return x.saving > y.saving; });
            std::vector<char> takenF(nF, 0);
            for (const Opt &o : opts) {
                if (col[o.r] >= 0 || takenF[o.j]) continue;
                col[o.r] = o.j;
                takenF[o.j] = 1;
            }
        }
    }
    double dist = 0.0;
    std::vector<char> close(nF, 0);
    const bool loose = rep != nullptr;
    for (int r = 0; r < n; ++r) {
        const int j = col[r];
        const Singularity &s = S[r < nI ? I[r] : B[r - nI]];
        if (j >= 0 && j < nF && F[j].sign == s.sign) {
            E += rowCost(r, j);
            if (r < nI) {
                ++R.matched;
                const double dd = std::min(1.0, normP(s.x - F[j].x) / r0);
                dist += dd;
                if (dd <= 0.5) close[j] = 1;
                else if (loose) R.looseVertices.push_back(s.vertex);
            } else {
                ++R.absorbed;
                if (loose) R.looseVertices.push_back(s.vertex);
            }
        } else {
            if (r < nI) {
                E += w.wExtra;
                ++R.extra;
            }
            if (loose) R.looseVertices.push_back(s.vertex);
        }
    }
    if (loose) for (int j = 0; j < nF; ++j) if (!close[j]) R.looseCones.push_back(F[j].x);
    R.missing = nF - R.matched - R.absorbed;
    R.meanDistance = R.matched > 0 ? dist / R.matched : 0.0;
    R.energy = E;
    if (rep) *rep = R;
    return E;
}

double ReferenceField::singularityEnergy(const SquareCarrier &C, const SingularityWeights &w,
                                         SingularityReport *rep) const {
    return singularityEnergy(carrierSingularities(C), cones_, w, rep);
}
