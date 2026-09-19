#include "ATLAS/CoarseDomain.hxx"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <unordered_map>

#include <Eigen/Sparse>

#include "triangle/TriangleMesher.hpp"

namespace {

double wrap(double s, double L) {
    if (L <= 0.0) return 0.0;
    s = std::fmod(s, L);
    return s < 0.0 ? s + L : s;
}

// Integral of 1/h over a segment of length l whose spacing runs linearly from
// h0 to h1: the number of coarse segments it asks for.
double density(double l, double h0, double h1) {
    if (std::fabs(h1 - h0) <= 1e-12 * std::max(h0, h1)) return l / h0;
    return l * std::log(h1 / h0) / (h1 - h0);
}

double minAngleOf(const Point &a, const Point &b, const Point &c) {
    auto ang = [](const Point &p, const Point &q, const Point &r) {
        const Point u = q - p, w = r - p;
        return std::fabs(std::atan2(cross2(u, w), dotP(u, w)));
    };
    return std::min({ang(a, b, c), ang(b, c, a), ang(c, a, b)});
}

// In the three-quad split a boundary vertex's valence is its number of
// triangles, so a quality triangulation that cuts a right-angled corner into
// two triangles hands the layout a corner with two cells -- whose separatrix
// then crosses the domain on the diagonal -- and Stage 6 has to undo it
// later, which it does badly (a corner's fan is only in a cavity whole). An
// edge flip at the corner does it for nothing. So: while some boundary vertex
// has more triangles than its layout angle asks for, flip the interior edge at
// it that lowers the total boundary excess and leaves the best smallest
// angle, never below minAngle.
int flipBoundaryValences(const std::vector<Point> &V, std::vector<Triangle> &T, const std::vector<int> &target,
                         double minAngle) {
    int flips = 0;
    auto key = [](int a, int b) { return (static_cast<long long>(std::min(a, b)) << 32) | std::max(a, b); };
    for (int pass = 0; pass < 64; ++pass) {
        std::unordered_map<long long, std::vector<int>> edgeTris;
        std::vector<int> deg(V.size(), 0);
        for (int t = 0; t < static_cast<int>(T.size()); ++t) {
            for (int k = 0; k < 3; ++k) {
                edgeTris[key(T[t][k], T[t][(k + 1) % 3])].push_back(t);
                ++deg[T[t][k]];
            }
        }
        auto excess = [&](int v, int d) { return target[v] >= 1 ? std::abs(d - target[v]) : 0; };
        bool changed = false;
        for (int v = 0; v < static_cast<int>(V.size()) && !changed; ++v) {
            if (target[v] < 1 || deg[v] <= target[v]) continue;
            int bestT1 = -1, bestT2 = -1;
            Triangle n1{}, n2{};
            double bestQ = minAngle;
            for (int t1 = 0; t1 < static_cast<int>(T.size()); ++t1) {
                int k = -1;
                for (int j = 0; j < 3; ++j) if (T[t1][j] == v) k = j;
                if (k < 0) continue;
                const int w = T[t1][(k + 1) % 3], a = T[t1][(k + 2) % 3];   // T1 = (v, w, a)
                const auto &sh = edgeTris[key(v, w)];
                if (sh.size() != 2) continue;                                // a boundary edge
                const int t2 = sh[0] == t1 ? sh[1] : sh[0];
                int b = -1;
                for (int j = 0; j < 3; ++j) if (T[t2][j] != v && T[t2][j] != w) b = T[t2][j];
                if (b < 0 || edgeTris.count(key(a, b))) continue;
                // Quad v, b, w, a; new triangles (v, b, a) and (b, w, a).
                if (!(cross2(V[b] - V[v], V[a] - V[v]) > 0.0) || !(cross2(V[w] - V[b], V[a] - V[b]) > 0.0)) continue;
                const int gain = excess(v, deg[v] - 1) + excess(w, deg[w] - 1) + excess(a, deg[a] + 1) +
                                 excess(b, deg[b] + 1) - excess(v, deg[v]) - excess(w, deg[w]) - excess(a, deg[a]) -
                                 excess(b, deg[b]);
                if (gain >= 0) continue;
                const double q = std::min(minAngleOf(V[v], V[b], V[a]), minAngleOf(V[b], V[w], V[a]));
                if (q > bestQ) {
                    bestQ = q;
                    bestT1 = t1;
                    bestT2 = t2;
                    n1 = {v, b, a};
                    n2 = {b, w, a};
                }
            }
            if (bestT1 >= 0) {
                T[bestT1] = n1;
                T[bestT2] = n2;
                ++flips;
                changed = true;
            }
        }
        if (!changed) break;
    }
    return flips;
}

} // namespace

CoarseDomain::CoarseDomain(const PlanarDomain &fine, const Options &opts) : fine_(fine), opts_(opts) {
    const auto t0 = std::chrono::steady_clock::now();
    if (!fine_.getReport().valid) {
        report_.reason = "the fine domain did not pass Stage 1";
        return;
    }
    buildArcs();
    std::vector<std::vector<int>> samples;
    if (sample(samples) && triangulate(samples)) {
        buildField();
        report_.valid = true;
    }
    report_.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

void CoarseDomain::buildArcs() {
    const Mesh &M = fine_.getMesh();
    arcIndex_.assign(M.vertices.size(), -1);
    for (const PlanarDomain::Loop &L : fine_.loops) {
        Arc a;
        a.vertices = L.vertices;
        a.edges = L.edges;
        const int n = static_cast<int>(L.vertices.size());
        a.s.resize(n);
        double acc = 0.0;
        for (int i = 0; i < n; ++i) {
            a.s[i] = acc;
            acc += normP(M.vertices[L.vertices[(i + 1) % n]] - M.vertices[L.vertices[i]]);
            arcIndex_[L.vertices[i]] = i;
        }
        a.length = acc;
        arcs_.push_back(std::move(a));
    }
}

// ---------------------------------------------------------------------------
// The spacing along one loop: min of the global bound, the gap and the
// curvature, graded.
// ---------------------------------------------------------------------------
std::vector<double> CoarseDomain::spacing(int l) const {
    const Mesh &M = fine_.getMesh();
    const Arc &A = arcs_[l];
    const int n = static_cast<int>(A.vertices.size());
    const double scale = fine_.getReport().scale;
    const double hmax = opts_.maxSpacing * scale;
    // Below this the fine boundary itself is the limit, and a sharp corner --
    // where both bounds below go to zero at the apex -- would otherwise ask
    // for an unbounded number of samples.
    const double hmin = std::max(0.01 * scale, 1e-12);
    auto P = [&](int i) { return M.vertices[A.vertices[((i % n) + n) % n]]; };
    auto seg = [&](int i) { return normP(P(i + 1) - P(i)); };

    std::vector<double> h(n, hmax);

    // Curvature: the turning of the loop per unit length over a window, the
    // protected corners excluded (they are chain ends, not curvature).
    std::vector<double> turn(n, 0.0), share(n, 0.0);
    for (int i = 0; i < n; ++i) {
        const int v = A.vertices[i];
        turn[i] = fine_.protectedVertex[v] ? 0.0 : std::fabs(M_PI - fine_.targetAngle[v]);
        share[i] = 0.5 * (seg(i - 1) + seg(i));
    }
    const double W = std::min(0.5 * hmax, 0.25 * A.length);
    for (int i = 0; i < n; ++i) {
        double T = turn[i], S = share[i];
        double d = 0.0;
        for (int k = 1; k < n; ++k) {           // forward
            d += seg(i + k - 1);
            if (d > W) break;
            const int j = (i + k) % n;
            T += turn[j];
            S += share[j];
        }
        d = 0.0;
        for (int k = 1; k < n; ++k) {           // backward
            d += seg(i - k);
            if (d > W) break;
            const int j = ((i - k) % n + n) % n;
            T += turn[j];
            S += share[j];
        }
        if (T > 1e-9 && S > 0.0) h[i] = std::min(h[i], opts_.curvatureFraction * S / T);
    }

    // Gap: the radius of the largest disk inside the domain that touches dS
    // at this vertex, |w|^2 / (2 w.n) minimised over the other boundary
    // vertices in front of it -- across a narrow passage, half its width.
    // Only vertices far from this one *along* the boundary count (another
    // loop, or an arc more than three times the chord): near a convex corner
    // the disk shrinks to nothing, but a corner is not a passage -- a block
    // fits a right angle as it is -- and counting it would grade the samples
    // down to the input's own boundary spacing there, making the coarse
    // domain depend on the input triangulation (the dependence Sec. 13.3's
    // last row asks to be rid of).
    for (int i = 0; i < n; ++i) {
        const Point x = P(i);
        Point t = P(i + 1) - P(i - 1);
        const double tl = normP(t);
        if (tl <= 0.0) continue;
        t = t / tl;
        const Point nrm{-t[1], t[0]};
        double r = std::numeric_limits<double>::infinity();
        for (int m = 0; m < static_cast<int>(arcs_.size()); ++m) {
            for (int f : arcs_[m].vertices) {
                if (f == A.vertices[i]) continue;
                const Point w = M.vertices[f] - x;
                if (m == l) {
                    const double ds = std::fabs(arcs_[m].s[arcIndex_[f]] - A.s[i]);
                    const double arc = std::min(ds, A.length - ds);
                    if (arc <= 3.0 * normP(w)) continue;
                }
                const double dn = dotP(w, nrm);
                if (dn <= 1e-12 * scale) continue;
                r = std::min(r, dotP(w, w) / (2.0 * dn));
            }
        }
        if (std::isfinite(r)) h[i] = std::min(h[i], opts_.gapFraction * 2.0 * r);
    }

    for (double &x : h) x = std::max(x, hmin);
    // Grading, round the loop twice each way.
    for (int pass = 0; pass < 2; ++pass) {
        for (int k = 0; k < 2 * n; ++k) {
            const int i = k % n, j = (k + 1) % n;
            h[j] = std::min(h[j], h[i] + opts_.grading * seg(i));
        }
        for (int k = 2 * n; k > 0; --k) {
            const int i = k % n, j = (k - 1) % n;
            h[j] = std::min(h[j], h[i] + opts_.grading * seg(j));
        }
    }
    return h;
}

// ---------------------------------------------------------------------------
// Samples: every protected vertex, then each chain between two of them cut at
// equal steps of the integrated density, snapped to input vertices.
// ---------------------------------------------------------------------------
bool CoarseDomain::sample(std::vector<std::vector<int>> &samples) {
    const Mesh &M = fine_.getMesh();
    samples.assign(arcs_.size(), {});
    report_.spacingMin = std::numeric_limits<double>::infinity();
    report_.spacingMax = 0.0;
    for (int l = 0; l < static_cast<int>(arcs_.size()); ++l) {
        const Arc &A = arcs_[l];
        const int n = static_cast<int>(A.vertices.size());
        if (n < 3) {
            report_.reason = "a boundary loop has fewer than three vertices";
            return false;
        }
        const std::vector<double> h = spacing(l);
        for (double x : h) {
            report_.spacingMin = std::min(report_.spacingMin, x);
            report_.spacingMax = std::max(report_.spacingMax, x);
        }
        auto seg = [&](int i) {
            return normP(M.vertices[A.vertices[(i + 1) % n]] - M.vertices[A.vertices[i % n]]);
        };

        std::vector<int> corners;
        for (int i = 0; i < n; ++i) if (fine_.protectedVertex[A.vertices[i]]) corners.push_back(i);
        struct Chain { int a, len; bool closed; std::vector<double> dens; int segs; };
        std::vector<Chain> chains;
        if (corners.empty()) {
            chains.push_back({0, n, true, {}, 0});
        } else {
            for (size_t k = 0; k < corners.size(); ++k) {
                const int a = corners[k];
                const int b = corners[(k + 1) % corners.size()];
                const int len = corners.size() == 1 ? n : ((b - a) % n + n) % n;
                chains.push_back({a, len, corners.size() == 1, {}, 0});
            }
        }
        int total = 0;
        for (Chain &c : chains) {
            c.dens.assign(c.len + 1, 0.0);
            for (int k = 0; k < c.len; ++k) {
                const int i = (c.a + k) % n;
                c.dens[k + 1] = c.dens[k] + density(seg(i), h[i], h[(i + 1) % n]);
            }
            const int want = static_cast<int>(std::lround(c.dens[c.len]));
            const int least = c.closed ? opts_.minClosedSamples : 1;
            c.segs = std::min(c.len, std::max(least, want));
            total += c.segs;
        }
        // A loop needs three samples to be a polygon at all.
        while (total < 3) {
            Chain *best = nullptr;
            for (Chain &c : chains) {
                if (c.segs >= c.len) continue;
                if (!best || c.dens[c.len] / c.segs > best->dens[best->len] / best->segs) best = &c;
            }
            if (!best) {
                report_.reason = "a boundary loop is too short to sample";
                return false;
            }
            ++best->segs;
            ++total;
        }
        for (const Chain &c : chains) {
            samples[l].push_back(A.vertices[c.a]);
            int prev = 0;
            const double D = c.dens[c.len];
            for (int j = 1; j < c.segs; ++j) {
                const double target = D * j / c.segs;
                // The vertex whose density is nearest the target, leaving room
                // for the samples still to come.
                const int lo = prev + 1, hi = c.len - (c.segs - j);
                int k = static_cast<int>(std::lower_bound(c.dens.begin() + lo, c.dens.begin() + hi + 1, target) -
                                         c.dens.begin());
                if (k > hi) k = hi;
                if (k > lo && std::fabs(c.dens[k - 1] - target) < std::fabs(c.dens[k] - target)) --k;
                samples[l].push_back(A.vertices[(c.a + k) % n]);
                prev = k;
            }
        }
        report_.chains += static_cast<int>(chains.size());
        report_.samples += static_cast<int>(samples[l].size());
    }
    return true;
}

// ---------------------------------------------------------------------------
// Triangle on the sampled loops, then Stage 1 on the result with the layout
// angles of the input.
// ---------------------------------------------------------------------------
bool CoarseDomain::triangulate(const std::vector<std::vector<int>> &samples) {
    using namespace triangle_wrapper;
    const Mesh &M = fine_.getMesh();
    TriangleMesher2D::MeshInput in;
    std::vector<int> fineOf;
    for (size_t l = 0; l < samples.size(); ++l) {
        const int first = static_cast<int>(in.vertlist.size());
        const int m = static_cast<int>(samples[l].size());
        std::vector<std::array<int, 2>> segs;
        for (int k = 0; k < m; ++k) {
            const Point &p = M.vertices[samples[l][k]];
            in.vertlist.push_back({p[0], p[1]});
            fineOf.push_back(samples[l][k]);
            segs.push_back({first + k, first + (k + 1) % m});
        }
        in.segment_loops.push_back(std::move(segs));
        in.type.push_back(fine_.loops[l].outer ? 0 : 1);
    }
    in.h = opts_.maxSpacing * fine_.getReport().scale;
    TriangleMesher2D::Options o;
    o.min_angle_degrees = opts_.minAngle;
    o.suppress_boundary_splitting = true;
    TriangleMesher2D mesher(o);
    TriangleMesher2D::MeshOutput out;
    try {
        out = mesher.triangulate(in);
    } catch (const std::exception &e) {
        report_.reason = std::string("Triangle failed: ") + e.what();
        return false;
    }

    std::vector<Point> verts(out.verts.size());
    for (size_t i = 0; i < out.verts.size(); ++i) verts[i] = {out.verts[i][0], out.verts[i][1]};
    std::vector<Triangle> tris(out.triangles.size());
    for (size_t t = 0; t < out.triangles.size(); ++t) {
        tris[t] = {out.triangles[t][0], out.triangles[t][1], out.triangles[t][2]};
        // Triangle writes counter-clockwise triangles; make sure of it.
        const double A = cross2(verts[tris[t][1]] - verts[tris[t][0]], verts[tris[t][2]] - verts[tris[t][0]]);
        if (A < 0.0) std::swap(tris[t][1], tris[t][2]);
    }
    {
        std::vector<int> target(verts.size(), -1);
        for (size_t k = 0; k < fineOf.size() && k < verts.size(); ++k) {
            const long t = std::lround(fine_.targetAngle[fineOf[k]] / M_PI_2);
            target[k] = static_cast<int>(std::max(1L, std::min(4L, t)));
        }
        report_.flips = flipBoundaryValences(verts, tris, target, opts_.minFlipAngle * M_PI / 180.0);
    }
    const int material = M.triangleMatId.empty() ? 1 : M.triangleMatId[0];
    mesh_ = std::make_shared<Mesh>(verts, tris, std::vector<int>(tris.size(), material));

    const int NV = static_cast<int>(verts.size());
    fineVertex_.assign(NV, -1);
    // Triangle keeps the input points first and in order (switch z).
    for (size_t k = 0; k < fineOf.size() && static_cast<int>(k) < NV; ++k) fineVertex_[k] = fineOf[k];

    PlanarDomain::Options po = fine_.getOptions();
    po.angleOverride.assign(NV, std::numeric_limits<double>::quiet_NaN());
    po.extraCorners.clear();
    for (int k = 0; k < NV; ++k) {
        const int f = fineVertex_[k];
        if (f < 0) continue;
        po.angleOverride[k] = fine_.targetAngle[f];
        if (fine_.protectedVertex[f]) po.extraCorners.push_back(k);
    }
    domain_ = std::make_unique<PlanarDomain>(*mesh_, po);
    report_.vertices = NV;
    report_.triangles = static_cast<int>(tris.size());
    if (!domain_->getReport().valid) {
        report_.reason = "the coarse triangulation does not pass Stage 1";
        return false;
    }
    if (domain_->loops.size() != fine_.loops.size()) {
        report_.reason = "the coarse triangulation has a different number of boundary loops";
        return false;
    }
    // Every coarse boundary vertex must be a sample: a point Triangle added on
    // a boundary segment would have no place on the fine boundary.
    for (int v : mesh_->boundaryVertices) {
        if (fineVertex_[v] < 0) {
            report_.reason = "Triangle put a point on a boundary segment";
            return false;
        }
    }
    edgeFrom_.assign(mesh_->edges.size(), -1);
    for (const PlanarDomain::Loop &L : domain_->loops) {
        for (size_t i = 0; i < L.edges.size(); ++i) edgeFrom_[L.edges[i]] = L.vertices[i];
        // The coarse loop must run the way its fine loop does.
        const int f0 = fineVertex_[L.vertices[0]];
        if (f0 < 0 || arcIndex_[f0] < 0) {
            report_.reason = "a coarse boundary loop has no fine counterpart";
            return false;
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
double CoarseDomain::forward(int loop, double a, double b) const {
    const double L = arcs_[loop].length;
    double d = wrap(b - a, L);
    if (d <= 1e-14 * L) d = L;
    return d;
}

Point CoarseDomain::pointAt(int loop, double s, int *edge, int *index) const {
    const Arc &A = arcs_[loop];
    const Mesh &M = fine_.getMesh();
    const int n = static_cast<int>(A.vertices.size());
    s = wrap(s, A.length);
    int i = static_cast<int>(std::upper_bound(A.s.begin(), A.s.end(), s) - A.s.begin()) - 1;
    i = std::max(0, std::min(n - 1, i));
    const double s0 = A.s[i];
    const double s1 = (i + 1 < n) ? A.s[i + 1] : A.length;
    const double t = s1 > s0 ? (s - s0) / (s1 - s0) : 0.0;
    if (edge) *edge = A.edges[i];
    if (index) *index = i;
    const Point &p = M.vertices[A.vertices[i]], &q = M.vertices[A.vertices[(i + 1) % n]];
    return p + (q - p) * t;
}

CoarseDomain::Location CoarseDomain::locate(const SquareCarrier &C, int v) const {
    Location out;
    const Mesh &M = fine_.getMesh();
    if (C.sourceVertex[v] >= 0) {
        const int f = fineVertex_[C.sourceVertex[v]];
        if (f < 0) return out;
        out.loop = fine_.loopOf[f];
        if (out.loop < 0) return out;
        out.s = arcs_[out.loop].s[arcIndex_[f]];
        return out;
    }
    const int e = C.sourceEdge[v];
    if (e < 0 || e >= static_cast<int>(edgeFrom_.size()) || edgeFrom_[e] < 0) return out;
    const int a = edgeFrom_[e];
    const int b = mesh_->edges[e][0] == a ? mesh_->edges[e][1] : mesh_->edges[e][0];
    const int fa = fineVertex_[a], fb = fineVertex_[b];
    if (fa < 0 || fb < 0) return out;
    const int loop = fine_.loopOf[fa];
    if (loop < 0 || fine_.loopOf[fb] != loop) return out;
    const Point pa = M.vertices[fa], pb = M.vertices[fb];
    const Point d = pb - pa;
    const double L2 = dotP(d, d);
    if (!(L2 > 0.0)) return out;
    const double t = std::max(0.0, std::min(1.0, dotP(C.vertices[v] - pa, d) / L2));
    const double sa = arcs_[loop].s[arcIndex_[fa]];
    out.loop = loop;
    out.s = wrap(sa + t * forward(loop, sa, arcs_[loop].s[arcIndex_[fb]]), arcs_[loop].length);
    return out;
}

// ---------------------------------------------------------------------------
// The reference cross field: harmonic u = exp(4 i theta), then a lookup grid.
// ---------------------------------------------------------------------------
void CoarseDomain::buildField() {
    const Mesh &M = fine_.getMesh();
    const Mesh &Cm = *mesh_;
    const int NV = static_cast<int>(Cm.vertices.size());
    std::vector<std::array<double, 2>> u(NV, {0.0, 0.0});
    std::vector<int> unknown(NV, -1);
    int nu = 0;
    for (int v = 0; v < NV; ++v) {
        const int f = fineVertex_[v];
        if (f >= 0 && !fine_.protectedVertex[f] && arcIndex_[f] >= 0) {
            const Arc &A = arcs_[fine_.loopOf[f]];
            const int n = static_cast<int>(A.vertices.size()), i = arcIndex_[f];
            const Point t = M.vertices[A.vertices[(i + 1) % n]] - M.vertices[A.vertices[(i + n - 1) % n]];
            const double th = 4.0 * std::atan2(t[1], t[0]);
            u[v] = {std::cos(th), std::sin(th)};
        } else {
            unknown[v] = nu++;
        }
    }
    if (nu > 0) {
        std::vector<Eigen::Triplet<double>> T;
        Eigen::VectorXd bx = Eigen::VectorXd::Zero(nu), by = Eigen::VectorXd::Zero(nu);
        std::vector<double> diag(nu, 1e-9);
        for (const auto &e : Cm.edges) {
            for (int k = 0; k < 2; ++k) {
                const int a = e[k], b = e[1 - k];
                const int i = unknown[a];
                if (i < 0) continue;
                diag[i] += 1.0;
                if (unknown[b] >= 0) {
                    T.emplace_back(i, unknown[b], -1.0);
                } else {
                    bx[i] += u[b][0];
                    by[i] += u[b][1];
                }
            }
        }
        for (int i = 0; i < nu; ++i) T.emplace_back(i, i, diag[i]);
        Eigen::SparseMatrix<double> L(nu, nu);
        L.setFromTriplets(T.begin(), T.end());
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(L);
        if (solver.info() == Eigen::Success) {
            const Eigen::VectorXd x = solver.solve(bx), y = solver.solve(by);
            for (int v = 0; v < NV; ++v) if (unknown[v] >= 0) u[v] = {x[unknown[v]], y[unknown[v]]};
        }
    }

    // Lookup grid: barycentric in the coarse triangle over each node, and
    // nodes outside the domain filled from their neighbours so that cells
    // along dS read the boundary's own direction.
    Point lo = Cm.vertices[0], hi = Cm.vertices[0];
    for (const Point &p : Cm.vertices) {
        lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
        hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
    }
    fieldStep_ = std::max(hi[0] - lo[0], hi[1] - lo[1]) / 128.0;
    if (!(fieldStep_ > 0.0)) fieldStep_ = 1.0;
    fieldLo_ = {lo[0] - fieldStep_, lo[1] - fieldStep_};
    fieldNx_ = static_cast<int>(std::ceil((hi[0] - lo[0]) / fieldStep_)) + 3;
    fieldNy_ = static_cast<int>(std::ceil((hi[1] - lo[1]) / fieldStep_)) + 3;
    fieldGrid_.assign(static_cast<size_t>(fieldNx_) * fieldNy_, {0.0, 0.0});
    std::vector<char> known(fieldGrid_.size(), 0);
    for (const Triangle &t : Cm.triangles) {
        const Point &a = Cm.vertices[t[0]], &b = Cm.vertices[t[1]], &c = Cm.vertices[t[2]];
        const double det = cross2(b - a, c - a);
        if (!(std::fabs(det) > 0.0)) continue;
        const int i0 = std::max(0, static_cast<int>(std::floor((std::min({a[0], b[0], c[0]}) - fieldLo_[0]) / fieldStep_)));
        const int i1 = std::min(fieldNx_ - 1, static_cast<int>(std::ceil((std::max({a[0], b[0], c[0]}) - fieldLo_[0]) / fieldStep_)));
        const int j0 = std::max(0, static_cast<int>(std::floor((std::min({a[1], b[1], c[1]}) - fieldLo_[1]) / fieldStep_)));
        const int j1 = std::min(fieldNy_ - 1, static_cast<int>(std::ceil((std::max({a[1], b[1], c[1]}) - fieldLo_[1]) / fieldStep_)));
        for (int j = j0; j <= j1; ++j) {
            for (int i = i0; i <= i1; ++i) {
                const Point p{fieldLo_[0] + i * fieldStep_, fieldLo_[1] + j * fieldStep_};
                const double l1 = cross2(c - b, p - b) / det, l2 = cross2(a - c, p - c) / det;
                const double l3 = 1.0 - l1 - l2;
                if (l1 < -1e-12 || l2 < -1e-12 || l3 < -1e-12) continue;
                const size_t k = static_cast<size_t>(i) + static_cast<size_t>(fieldNx_) * j;
                for (int d = 0; d < 2; ++d) fieldGrid_[k][d] = l1 * u[t[0]][d] + l2 * u[t[1]][d] + l3 * u[t[2]][d];
                known[k] = 1;
            }
        }
    }
    for (int sweep = 0; sweep < 4; ++sweep) {
        std::vector<char> next = known;
        for (int j = 0; j < fieldNy_; ++j) {
            for (int i = 0; i < fieldNx_; ++i) {
                const size_t k = static_cast<size_t>(i) + static_cast<size_t>(fieldNx_) * j;
                if (known[k]) continue;
                std::array<double, 2> s{0.0, 0.0};
                int n = 0;
                for (int dj = -1; dj <= 1; ++dj) {
                    for (int di = -1; di <= 1; ++di) {
                        const int ii = i + di, jj = j + dj;
                        if (ii < 0 || jj < 0 || ii >= fieldNx_ || jj >= fieldNy_) continue;
                        const size_t kk = static_cast<size_t>(ii) + static_cast<size_t>(fieldNx_) * jj;
                        if (!known[kk]) continue;
                        s[0] += fieldGrid_[kk][0];
                        s[1] += fieldGrid_[kk][1];
                        ++n;
                    }
                }
                if (n > 0) {
                    fieldGrid_[k] = {s[0] / n, s[1] / n};
                    next[k] = 1;
                }
            }
        }
        known.swap(next);
    }
    fieldMax_ = 1.0;
}

double CoarseDomain::crossAngle(const Point &p, double *weight) const {
    if (fieldGrid_.empty()) {
        if (weight) *weight = 0.0;
        return 0.0;
    }
    const double x = (p[0] - fieldLo_[0]) / fieldStep_, y = (p[1] - fieldLo_[1]) / fieldStep_;
    const int i = std::max(0, std::min(fieldNx_ - 2, static_cast<int>(std::floor(x))));
    const int j = std::max(0, std::min(fieldNy_ - 2, static_cast<int>(std::floor(y))));
    const double fx = std::max(0.0, std::min(1.0, x - i)), fy = std::max(0.0, std::min(1.0, y - j));
    auto at = [&](int a, int b) { return fieldGrid_[static_cast<size_t>(a) + static_cast<size_t>(fieldNx_) * b]; };
    std::array<double, 2> u{0.0, 0.0};
    for (int d = 0; d < 2; ++d) {
        u[d] = (1 - fx) * (1 - fy) * at(i, j)[d] + fx * (1 - fy) * at(i + 1, j)[d] + (1 - fx) * fy * at(i, j + 1)[d] +
               fx * fy * at(i + 1, j + 1)[d];
    }
    if (weight) *weight = std::min(1.0, std::hypot(u[0], u[1]) / fieldMax_);
    double th = std::atan2(u[1], u[0]) / 4.0;
    th = std::fmod(th, M_PI_2);
    if (th < 0.0) th += M_PI_2;
    return th;
}
