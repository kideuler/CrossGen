#include "ZIPLINE/ConformalField.hxx"

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <queue>
#include <sstream>

namespace {

// The representative of a mod pi/2 nearest ref.
inline double liftNear(double ref, double a) { return a + M_PI_2 * std::round((ref - a) / M_PI_2); }
// Into (-pi/4, pi/4]: how far a cross angle is from zero.
inline double wrapQuarter(double a) { return a - M_PI_2 * std::round(a / M_PI_2); }

// Table 1 of Viertel et al. with its own intervals, exactly as CrossField::
// initialize sets the boundary data and SeparatrixTrace::buildBoundary reads the
// corners back: the corners the second field turns at have to be the ones the
// first was told to turn at and the tracing will launch from.
int quartersOf(double interiorAngle) {
    if (interiorAngle < 0.75 * M_PI) return 1;
    if (interiorAngle <= 1.25 * M_PI) return 2;
    if (interiorAngle <= 1.75 * M_PI) return 3;
    return std::max(4, static_cast<int>(std::lround(interiorAngle / M_PI_2)));
}

inline Point midpoint(const Mesh &m, int t, int k) {
    const Triangle &tri = m.triangles[t];
    return (m.vertices[tri[k]] + m.vertices[tri[(k + 1) % 3]]) * 0.5;
}

} // namespace

ConformalField::ConformalField(std::shared_ptr<Mesh> mesh, std::vector<Singular> pattern,
                               const Options &opts, const CrossField *first)
    : mesh_(std::move(mesh)), opts_(opts), first_(first), sing_(std::move(pattern)) {
    build();
    first_ = nullptr;   // read during construction only
}

double ConformalField::singularPart(const Point &x) const {
    double s = 0.0;
    for (const Singular &q : sing_) s += 0.25 * q.d * std::atan2(x[1] - q.pos[1], x[0] - q.pos[0]);
    return s;
}

double ConformalField::angleAt(int t, const Point &x, double ref) const {
    const Mesh &m = *mesh_;
    const Triangle &tri = m.triangles[t];
    const Point &a = m.vertices[tri[0]], &b = m.vertices[tri[1]], &c = m.vertices[tri[2]];
    const double den = cross2(b - a, c - a);
    double l1 = 1.0 / 3.0, l2 = 1.0 / 3.0;
    if (std::fabs(den) > 1e-300) {
        l1 = cross2(x - a, c - a) / den;
        l2 = cross2(b - a, x - a) / den;
    }
    const double r0 = regular_[tri[0]];
    const double r1 = liftNear(r0, regular_[tri[1]]);
    const double r2 = liftNear(r0, regular_[tri[2]]);
    return liftNear(ref, (1.0 - l1 - l2) * r0 + l1 * r1 + l2 * r2 + singularPart(x));
}

// ---------------------------------------------------------------------------
// build()
// ---------------------------------------------------------------------------
void ConformalField::build() {
    const Mesh &M = *mesh_;
    const int nV = static_cast<int>(M.vertices.size());
    const int nT = static_cast<int>(M.triangles.size());
    report_ = Report{};
    auto fail = [&](const std::string &why) { report_.built = false; report_.failure = why; };
    if (nT == 0) return fail("no triangles");

    double meanEdge = 0.0;
    for (const auto &e : M.edges) meanEdge += normP(M.vertices[e[1]] - M.vertices[e[0]]);
    meanEdge = M.edges.empty() ? 1.0 : meanEdge / static_cast<double>(M.edges.size());

    std::vector<char> used(nV, 0);
    for (const Triangle &t : M.triangles) used[t[0]] = used[t[1]] = used[t[2]] = 1;

    // ---- the boundary, as loops with the domain on their left --------------
    //
    // Read off the geometry rather than the vertex order of the owning
    // triangle, so nothing below depends on the mesh being counter-clockwise.
    std::vector<int> nextOf(nV, -1), prevOf(nV, -1), edgeAfter(nV, -1);
    for (const int e : M.boundaryEdges) {
        const int t = (M.edgeTriangles[e][0] >= 0) ? M.edgeTriangles[e][0] : M.edgeTriangles[e][1];
        if (t < 0) continue;
        int a = M.edges[e][0], b = M.edges[e][1];
        int c = -1;
        for (int k = 0; k < 3; ++k)
            if (M.triangles[t][k] != a && M.triangles[t][k] != b) c = M.triangles[t][k];
        if (c < 0) continue;
        if (cross2(M.vertices[b] - M.vertices[a], M.vertices[c] - M.vertices[a]) < 0.0) std::swap(a, b);
        if (nextOf[a] >= 0 || prevOf[b] >= 0)
            return fail("the boundary is not a disjoint union of simple loops (a pinched vertex)");
        nextOf[a] = b;
        prevOf[b] = a;
        edgeAfter[a] = e;
    }
    std::vector<std::vector<int>> loops;
    std::vector<int> loopOf(nV, -1);
    for (int v = 0; v < nV; ++v) {
        if (nextOf[v] < 0 || loopOf[v] >= 0) continue;
        std::vector<int> loop;
        int u = v;
        while (u >= 0 && loopOf[u] < 0) {
            loopOf[u] = static_cast<int>(loops.size());
            loop.push_back(u);
            u = nextOf[u];
        }
        if (u != v) return fail("a boundary loop does not close");
        loops.push_back(std::move(loop));
    }
    report_.loops = static_cast<int>(loops.size());

    // Interior angle at each boundary vertex, and from it the turning angle
    // tau_b = pi - alpha_b (Sec. 2.3) and the Table 1 quarters q_b.
    std::vector<double> interior(nV, 0.0);
    for (const Triangle &tri : M.triangles) {
        for (int k = 0; k < 3; ++k) {
            const int v = tri[k];
            if (loopOf[v] < 0) continue;
            const Point a = M.vertices[tri[(k + 1) % 3]] - M.vertices[v];
            const Point b = M.vertices[tri[(k + 2) % 3]] - M.vertices[v];
            const double na = normP(a), nb = normP(b);
            if (na < 1e-300 || nb < 1e-300) continue;
            interior[v] += std::acos(std::max(-1.0, std::min(1.0, dotP(a, b) / (na * nb))));
        }
    }
    std::vector<int> quarters(nV, 2);
    std::vector<double> loopTurning(loops.size(), 0.0);
    for (size_t l = 0; l < loops.size(); ++l)
        for (const int b : loops[l]) {
            quarters[b] = quartersOf(interior[b]);
            if (quarters[b] != 2) ++report_.corners;
            loopTurning[l] += M_PI - interior[b];
        }

    // ---- components, and the outer loop of each -----------------------------
    std::vector<int> compOf(nT, -1);
    int nComp = 0;
    for (int t0 = 0; t0 < nT; ++t0) {
        if (compOf[t0] >= 0) continue;
        std::queue<int> q;
        q.push(t0);
        compOf[t0] = nComp;
        while (!q.empty()) {
            const int t = q.front();
            q.pop();
            for (int k = 0; k < 3; ++k) {
                const int n = M.triangleAdjacency[t][k];
                if (n >= 0 && compOf[n] < 0) { compOf[n] = nComp; q.push(n); }
            }
        }
        ++nComp;
    }
    report_.components = nComp;
    std::vector<int> compOfVertex(nV, -1);
    for (int t = 0; t < nT; ++t)
        for (int k = 0; k < 3; ++k) compOfVertex[M.triangles[t][k]] = compOf[t];
    std::vector<std::vector<int>> loopsOfComp(nComp);
    for (size_t l = 0; l < loops.size(); ++l) loopsOfComp[compOfVertex[loops[l][0]]].push_back(static_cast<int>(l));
    std::vector<int> outerOf(nComp, -1);
    for (int c = 0; c < nComp; ++c) {
        for (const int l : loopsOfComp[c]) {
            // +2 pi round the outer loop, -2 pi round a hole (Sec. 2.3).
            if (loopTurning[l] > 0.0) {
                if (outerOf[c] >= 0) return fail("a component with two outer loops");
                outerOf[c] = l;
            } else {
                ++report_.holes;
            }
        }
        if (outerOf[c] < 0) return fail("a component with no outer loop");
    }

    // ---- Sec. 5.4's automatic edits -----------------------------------------
    //
    // Locate each singularity's triangle afresh from its position if the one
    // given does not hold it, merge singularities sharing a triangle (their
    // local models cannot both be the one the sweep uses), then cancel dipoles.
    auto bary = [&](int t, const Point &x) {
        const Triangle &tri = M.triangles[t];
        const Point &a = M.vertices[tri[0]], &b = M.vertices[tri[1]], &c = M.vertices[tri[2]];
        const double den = cross2(b - a, c - a);
        const double l1 = cross2(x - a, c - a) / den, l2 = cross2(b - a, x - a) / den;
        return std::array<double, 3>{{1.0 - l1 - l2, l1, l2}};
    };
    for (Singular &s : sing_) {
        if (s.triangle < 0 || s.triangle >= nT) return fail("a singularity with no triangle");
        const auto l = bary(s.triangle, s.pos);
        if (std::min({l[0], l[1], l[2]}) < -1e-9) return fail("a singularity outside its triangle");
    }
    {
        std::vector<Singular> kept;
        std::vector<int> slot(nT, -1);
        for (const Singular &s : sing_) {
            if (s.d == 0) continue;
            int &k = slot[s.triangle];
            if (k < 0) { k = static_cast<int>(kept.size()); kept.push_back(s); continue; }
            // Two in one triangle: one singularity of the summed index, at
            // the index-weighted mean (or simply the mean when they cancel).
            Singular &o = kept[k];
            const int wsum = std::abs(o.d) + std::abs(s.d);
            o.pos = (o.pos * static_cast<double>(std::abs(o.d)) + s.pos * static_cast<double>(std::abs(s.d))) /
                    static_cast<double>(wsum);
            o.d += s.d;
            ++report_.merged;
        }
        sing_.clear();
        for (const Singular &s : kept)
            if (s.d != 0) sing_.push_back(s);
    }
    if (opts_.cancelDipoles && opts_.dipoleRadius > 0.0) {
        const double R = opts_.dipoleRadius * meanEdge;
        struct Pair { double dist; int a, b; };
        std::vector<Pair> pairs;
        for (size_t i = 0; i < sing_.size(); ++i)
            for (size_t j = i + 1; j < sing_.size(); ++j) {
                if (sing_[i].d + sing_[j].d != 0) continue;
                const double d = normP(sing_[i].pos - sing_[j].pos);
                if (d < R) pairs.push_back({d, static_cast<int>(i), static_cast<int>(j)});
            }
        std::sort(pairs.begin(), pairs.end(), [](const Pair &x, const Pair &y) { return x.dist < y.dist; });
        std::vector<char> gone(sing_.size(), 0);
        for (const Pair &p : pairs) {
            if (gone[p.a] || gone[p.b]) continue;
            gone[p.a] = gone[p.b] = 1;
            ++report_.dipolesCancelled;
        }
        std::vector<Singular> kept;
        for (size_t i = 0; i < sing_.size(); ++i)
            if (!gone[i]) kept.push_back(sing_[i]);
        sing_.swap(kept);
    }
    // Sec. 12.3: a +1/4 singularity hugging a convex corner, moved into the
    // patch (Options::relocateCornerSingularities).
    if (opts_.relocateCornerSingularities && opts_.cornerRelocation > 0.0) {
        // Boundary length from b to the next corner along the loop, either way.
        auto runTo = [&](int b, bool forward) {
            double L = 0.0;
            int v = b;
            for (size_t guard = 0; guard < M.vertices.size(); ++guard) {
                const int n = forward ? nextOf[v] : prevOf[v];
                if (n < 0) break;
                L += normP(M.vertices[n] - M.vertices[v]);
                v = n;
                if (v == b || quarters[v] != 2) break;
            }
            return L;
        };
        auto locate = [&](const Point &x, int hint) {
            // A walk from the hint towards x, then a scan if it fails.
            int t = hint;
            for (int guard = 0; guard < 4096 && t >= 0; ++guard) {
                const auto l = bary(t, x);
                int worst = -1;
                double w = -1e-12;
                for (int k = 0; k < 3; ++k) if (l[k] < w) { w = l[k]; worst = k; }
                if (worst < 0) return t;
                // The edge opposite the most negative coordinate is local edge
                // (worst + 1) % 3.
                t = M.triangleAdjacency[t][(worst + 1) % 3];
            }
            for (int u = 0; u < nT; ++u) {
                const auto l = bary(u, x);
                if (l[0] >= -1e-12 && l[1] >= -1e-12 && l[2] >= -1e-12) return u;
            }
            return -1;
        };
        std::vector<char> holds(nT, 0);
        for (const Singular &q : sing_) holds[q.triangle] = 1;
        for (Singular &q : sing_) {
            if (q.d <= 0) continue;
            int corner = -1;
            double best = 1e300;
            for (const auto &loop : loops)
                for (const int b : loop) {
                    if (quarters[b] != 1) continue;
                    const double local = 0.5 * (normP(M.vertices[nextOf[b]] - M.vertices[b]) +
                                                normP(M.vertices[prevOf[b]] - M.vertices[b]));
                    const double d = normP(q.pos - M.vertices[b]);
                    const Triangle &tr = M.triangles[q.triangle];
                    const bool inFan = tr[0] == b || tr[1] == b || tr[2] == b;
                    if ((inFan || d < opts_.relocationRadius * local) && d < best) { best = d; corner = b; }
                }
            if (corner < 0) continue;
            const Point &c = M.vertices[corner];
            const Point e1 = normalizeP(M.vertices[nextOf[corner]] - c);
            const Point e2 = normalizeP(M.vertices[prevOf[corner]] - c);
            const Point bis = normalizeP(e1 + e2);
            if (normP(bis) < 0.5) continue;
            const double L = std::min(runTo(corner, true), runTo(corner, false));
            const double halfAngle = 0.5 * interior[corner];
            double D = opts_.cornerRelocation * (2.0 / 3.0) * L * std::cos(halfAngle);
            holds[q.triangle] = 0;
            int t = -1;
            Point x{0.0, 0.0};
            for (int tries = 0; tries < 8 && D > 0.25 * best; ++tries, D *= 0.5) {
                x = c + bis * D;
                t = locate(x, q.triangle);
                if (t >= 0 && !holds[t]) {
                    const auto l = bary(t, x);
                    if (std::min({l[0], l[1], l[2]}) > 0.05) break;
                    // Too near an edge of its triangle: its centroid instead.
                    const Triangle &tr = M.triangles[t];
                    x = (M.vertices[tr[0]] + M.vertices[tr[1]] + M.vertices[tr[2]]) / 3.0;
                    break;
                }
                t = -1;
            }
            if (t < 0) { holds[q.triangle] = 1; continue; }
            q.pos = x;
            q.triangle = t;
            holds[t] = 1;
            ++report_.relocated;
        }
    }
    report_.singularities = static_cast<int>(sing_.size());

    // ---- Poincare-Hopf, per component (Sec. 3.3) -----------------------------
    auto poincareHopf = [&](std::string *why) {
        std::vector<int> sum4(nComp, 0);
        for (const Singular &s : sing_) sum4[compOf[s.triangle]] += s.d;
        for (size_t l = 0; l < loops.size(); ++l)
            for (const int b : loops[l]) sum4[compOfVertex[b]] += 2 - quarters[b];
        for (int c = 0; c < nComp; ++c) {
            const int chi = 2 - static_cast<int>(loopsOfComp[c].size());
            if (sum4[c] != 4 * chi) {
                if (why) {
                    std::ostringstream os;
                    os << "Poincare-Hopf fails: indices add to " << sum4[c] << "/4 against chi = " << chi;
                    *why = os.str();
                }
                return false;
            }
        }
        return true;
    };
    {
        std::string why;
        if (!poincareHopf(&why) && first_ && first_->u_k.size() == nV) {
            // Re-read the corners off the first field (see the constructor):
            // the cross turns by the principal matchings along the boundary
            // edges, and across a stretch of boundary the turn Delta and the
            // tangent's turn tau fix the stretch's index, 2 pi k = tau - Delta
            // (Sec. 2.4). Exact over any stretch that begins and ends at
            // vertices carrying the tangent, which every vertex of a straight
            // or smooth run does -- but *not* per corner when two corners are
            // neighbours on the loop (a chamfer: two 225-degree corners a
            // short edge apart, geom033), since each would claim the edge
            // between them. So the stretch is a maximal run of consecutive
            // corners, its index is read whole, and Table 1's split of it is
            // kept and adjusted at the corner nearest its Table 1 threshold
            // until it adds up.
            std::vector<int> reread = quarters;
            int changed = 0;
            auto th = [&](int v) { return 0.25 * std::arg(first_->u_k[v]); };
            auto isCorner = [&](int v) {
                if (quarters[v] != 2) return true;
                const double delta = wrapQuarter(th(v) - th(prevOf[v])) + wrapQuarter(th(nextOf[v]) - th(v));
                return std::lround((M_PI - interior[v] - delta) / M_PI_2) != 0;
            };
            // The nearest Table 1 threshold to an interior angle, in the
            // direction that raises (up) or lowers q.
            auto toThreshold = [&](double a, int q, bool up) {
                const double t = up ? (q == 1 ? 0.75 * M_PI : q == 2 ? 1.25 * M_PI : q == 3 ? 1.75 * M_PI : 1e9)
                                    : (q == 2 ? 0.75 * M_PI : q == 3 ? 1.25 * M_PI : q >= 4 ? 1.75 * M_PI : -1e9);
                return std::fabs(a - t);
            };
            for (const auto &loop : loops) {
                const int m = static_cast<int>(loop.size());
                // Start the walk at a vertex that is no corner, so no run
                // straddles the seam; a loop of nothing but corners is read
                // vertex by vertex.
                int start = -1;
                for (int i = 0; i < m && start < 0; ++i) if (!isCorner(loop[i])) start = i;
                if (start < 0) continue;
                for (int i = 0; i < m;) {
                    const int v0 = loop[(start + i) % m];
                    if (!isCorner(v0)) { ++i; continue; }
                    std::vector<int> run;
                    while (i < m && isCorner(loop[(start + i) % m])) { run.push_back(loop[(start + i) % m]); ++i; }
                    const int before = prevOf[run.front()], after = nextOf[run.back()];
                    double delta = wrapQuarter(th(run.front()) - th(before));
                    double tau = 0.0;
                    for (size_t r = 0; r < run.size(); ++r) {
                        tau += M_PI - interior[run[r]];
                        const int nx = (r + 1 < run.size()) ? run[r + 1] : after;
                        delta += wrapQuarter(th(nx) - th(run[r]));
                    }
                    const int want = static_cast<int>(std::lround((tau - delta) / M_PI_2));   // sum (2 - q)
                    int have = 0;
                    for (const int v : run) have += 2 - reread[v];
                    while (have != want) {
                        const bool raise = have > want;   // more q, less (2 - q)
                        int best = -1;
                        double bestD = 1e300;
                        for (const int v : run) {
                            if (!raise && reread[v] <= 1) continue;
                            const double d = toThreshold(interior[v], reread[v], raise);
                            if (d < bestD) { bestD = d; best = v; }
                        }
                        if (best < 0) break;
                        reread[best] += raise ? 1 : -1;
                        have += raise ? -1 : 1;
                        ++changed;
                    }
                }
            }
            if (changed > 0) {
                std::swap(quarters, reread);
                if (poincareHopf(nullptr)) {
                    report_.cornersReread = changed;
                    report_.corners = 0;
                    for (const auto &loop : loops)
                        for (const int b : loop) if (quarters[b] != 2) ++report_.corners;
                } else {
                    std::swap(quarters, reread);
                }
            }
        }
        if (!poincareHopf(&why)) return fail(why);
        report_.poincareHopf = true;
    }
    quarters_ = quarters;

    // ---- the right-hand side of h (Sec. 6.3, with the singular part out) -----
    //
    // The change of the singular part along a straight segment that misses
    // every X_j, exactly: each term turns by the angle the segment subtends.
    auto singularChange = [&](const Point &p, const Point &q) {
        double s = 0.0;
        for (const Singular &j : sing_) {
            const double a0 = std::atan2(p[1] - j.pos[1], p[0] - j.pos[0]);
            const double a1 = std::atan2(q[1] - j.pos[1], q[0] - j.pos[0]);
            s += 0.25 * j.d * wrap_pi(a1 - a0);
        }
        return s;
    };
    std::vector<double> load(nV, 0.0);
    double loadScale = 0.0;
    for (const auto &loop : loops) {
        for (const int b : loop) {
            const int a = prevOf[b], c = nextOf[b];
            const Point &xb = M.vertices[b];
            const Point mIn = (M.vertices[a] + xb) * 0.5;
            const Point mOut = (xb + M.vertices[c]) * 0.5;
            const double tau = M_PI - interior[b];
            const double cornerTurn = M_PI_2 * (2 - quarters[b]);   // 2 pi k_b
            load[b] = tau - cornerTurn - (singularChange(mIn, xb) + singularChange(xb, mOut));
            loadScale = std::max(loadScale, std::fabs(load[b]));
        }
    }
    {
        std::vector<double> sum(nComp, 0.0);
        for (int v = 0; v < nV; ++v)
            if (loopOf[v] >= 0) sum[compOfVertex[v]] += load[v];
        for (int c = 0; c < nComp; ++c)
            report_.loadSum = std::max(report_.loadSum, std::fabs(sum[c]) / std::max(1.0, loadScale));
    }

    // ---- the stiffness, and h ------------------------------------------------
    std::vector<Point> gradH(nT, Point{0.0, 0.0});
    std::vector<std::array<Point, 3>> grads(nT);
    std::vector<double> area(nT, 0.0);
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(static_cast<size_t>(nT) * 9);
    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = M.triangles[t];
        const Point &p0 = M.vertices[tri[0]], &p1 = M.vertices[tri[1]], &p2 = M.vertices[tri[2]];
        const double a2 = cross2(p1 - p0, p2 - p0);
        area[t] = 0.5 * std::fabs(a2);
        if (std::fabs(a2) < 1e-300) continue;
        grads[t][0] = Point{(p1[1] - p2[1]) / a2, (p2[0] - p1[0]) / a2};
        grads[t][1] = Point{(p2[1] - p0[1]) / a2, (p0[0] - p2[0]) / a2};
        grads[t][2] = Point{(p0[1] - p1[1]) / a2, (p1[0] - p0[0]) / a2};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                trips.emplace_back(tri[i], tri[j], area[t] * dotP(grads[t][i], grads[t][j]));
    }
    // One vertex per component pinned to 0 (h is defined up to a constant),
    // and every vertex no triangle uses.
    std::vector<char> pinned(nV, 0);
    for (int v = 0; v < nV; ++v) if (!used[v]) pinned[v] = 1;
    for (int c = 0; c < nComp; ++c) pinned[loops[outerOf[c]][0]] = 1;
    std::vector<int> freeIdx(nV, -1);
    int nFree = 0;
    for (int v = 0; v < nV; ++v) if (!pinned[v]) freeIdx[v] = nFree++;
    {
        std::vector<Eigen::Triplet<double>> red;
        red.reserve(trips.size());
        for (const auto &tr : trips)
            if (freeIdx[tr.row()] >= 0 && freeIdx[tr.col()] >= 0)
                red.emplace_back(freeIdx[tr.row()], freeIdx[tr.col()], tr.value());
        Eigen::SparseMatrix<double> K(nFree, nFree);
        K.setFromTriplets(red.begin(), red.end());
        Eigen::VectorXd rhs(nFree);
        for (int v = 0; v < nV; ++v) if (freeIdx[v] >= 0) rhs[freeIdx[v]] = load[v];
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(K);
        if (ldlt.info() != Eigen::Success) return fail("the Neumann problem for h would not factorise");
        const Eigen::VectorXd sol = ldlt.solve(rhs);
        if (ldlt.info() != Eigen::Success) return fail("the Neumann problem for h would not solve");
        h_.assign(nV, 0.0);
        for (int v = 0; v < nV; ++v) if (freeIdx[v] >= 0) h_[v] = sol[freeIdx[v]];
    }
    std::vector<Point> w(nT, Point{0.0, 0.0});   // J grad h: the gradient of theta_h
    for (int t = 0; t < nT; ++t) {
        Point g{0.0, 0.0};
        for (int i = 0; i < 3; ++i) g = g + grads[t][i] * h_[M.triangles[t][i]];
        gradH[t] = g;
        w[t] = Point{-g[1], g[0]};
        report_.maxRegularTurnPerEdge = std::max(report_.maxRegularTurnPerEdge, normP(w[t]) * meanEdge);
    }

    // ---- theta_h, by Crouzeix-Raviart integration (Sec. 6.5) ----------------
    std::vector<std::array<double, 3>> cr(nT, std::array<double, 3>{{0.0, 0.0, 0.0}});
    std::vector<char> seen(nT, 0);
    auto fill = [&](int t, int k, double value) {
        const Point m = midpoint(M, t, k);
        for (int kk = 0; kk < 3; ++kk) cr[t][kk] = value + dotP(w[t], midpoint(M, t, kk) - m);
    };
    auto localEdgeOf = [&](int t, int e) {
        for (int k = 0; k < 3; ++k) if (M.triangleEdges[t][k] == e) return k;
        return -1;
    };
    for (int c = 0; c < nComp; ++c) {
        const auto &outer = loops[outerOf[c]];
        const int v0 = outer[0], v1 = nextOf[v0];
        const int e0 = edgeAfter[v0];
        const int t0 = (M.edgeTriangles[e0][0] >= 0) ? M.edgeTriangles[e0][0] : M.edgeTriangles[e0][1];
        const int k0 = localEdgeOf(t0, e0);
        if (k0 < 0) return fail("mesh edge tables disagree");
        const Point m0 = (M.vertices[v0] + M.vertices[v1]) * 0.5;
        // The anchor: the cross at the midpoint of the first edge of the outer
        // loop has a branch along it (P2 eq. 15's theta(P) = theta_P).
        const double phi0 = computeAngle(M.vertices[v1] - M.vertices[v0]);
        fill(t0, k0, phi0 - singularPart(m0));
        seen[t0] = 1;
        std::queue<int> q;
        q.push(t0);
        while (!q.empty()) {
            const int t = q.front();
            q.pop();
            for (int k = 0; k < 3; ++k) {
                const int n = M.triangleAdjacency[t][k];
                if (n < 0 || seen[n]) continue;
                const int kn = localEdgeOf(n, M.triangleEdges[t][k]);
                if (kn < 0) return fail("mesh edge tables disagree");
                fill(n, kn, cr[t][k]);
                seen[n] = 1;
                q.push(n);
            }
        }
    }
    for (int i = 0; i < static_cast<int>(M.edges.size()); ++i) {
        const int t = M.edgeTriangles[i][0], u = M.edgeTriangles[i][1];
        if (t < 0 || u < 0) continue;
        const int kt = localEdgeOf(t, i), ku = localEdgeOf(u, i);
        if (kt < 0 || ku < 0) continue;
        report_.maxMidpointMismatch =
            std::max(report_.maxMidpointMismatch, std::fabs(wrapQuarter(cr[t][kt] - cr[u][ku])));
    }

    // The tangency of the CR field at the midpoints of each loop's edges: the
    // offset delta of Sec. 6.8, which on the anchor's loop is zero by
    // construction.
    auto loopOffsets = [&](const std::vector<int> &loop, double &mean, double &spread) {
        std::complex<double> acc(0.0, 0.0);
        std::vector<double> d;
        for (const int a : loop) {
            const int e = edgeAfter[a], b = nextOf[a];
            const int t = (M.edgeTriangles[e][0] >= 0) ? M.edgeTriangles[e][0] : M.edgeTriangles[e][1];
            const int k = localEdgeOf(t, e);
            const Point m = (M.vertices[a] + M.vertices[b]) * 0.5;
            const double phi = computeAngle(M.vertices[b] - M.vertices[a]);
            const double off = wrapQuarter(cr[t][k] + singularPart(m) - phi);
            d.push_back(off);
            acc += std::polar(1.0, 4.0 * off);
        }
        mean = 0.25 * std::arg(acc);
        spread = 0.0;
        for (const double x : d) spread = std::max(spread, std::fabs(wrapQuarter(x - mean)));
    };

    // ---- the regular part per vertex -----------------------------------------
    //
    // Each triangle's affine theta_h read at the vertex, lifted round the fan
    // one triangle at a time (a sharp corner's fan turns by up to pi/4 in all,
    // more than a single lift against the first value could take), averaged.
    regular_.assign(nV, 0.0);
    {
        std::vector<std::vector<int>> trisOf(nV);
        for (int t = 0; t < nT; ++t)
            for (int k = 0; k < 3; ++k) trisOf[M.triangles[t][k]].push_back(t);
        auto valueAt = [&](int t, int v) {
            for (int k = 0; k < 3; ++k)
                if (M.triangles[t][k] == v) return cr[t][k] + dotP(w[t], M.vertices[v] - midpoint(M, t, k));
            return 0.0;
        };
        for (int v = 0; v < nV; ++v) {
            const auto &ts = trisOf[v];
            if (ts.empty()) continue;
            // Order the fan by walking across the edges at v: start at a
            // triangle with a boundary edge at v when there is one.
            auto edgesAtV = [&](int t, int out[2]) {
                int n = 0;
                for (int k = 0; k < 3; ++k) {
                    const int a = M.triangles[t][k], b = M.triangles[t][(k + 1) % 3];
                    if (a == v || b == v) out[n++] = k;
                }
            };
            int start = ts[0], startEdge = -1;
            for (const int t : ts) {
                int ek[2];
                edgesAtV(t, ek);
                for (const int k : ek)
                    if (M.triangleAdjacency[t][k] < 0) { start = t; startEdge = k; }
                if (startEdge >= 0) break;
            }
            std::vector<int> fan;
            {
                int t = start, came = startEdge;
                if (came < 0) { int ek[2]; edgesAtV(t, ek); came = ek[0]; }
                for (size_t guard = 0; guard <= ts.size() && t >= 0; ++guard) {
                    fan.push_back(t);
                    int ek[2];
                    edgesAtV(t, ek);
                    const int leave = (ek[0] == came) ? ek[1] : ek[0];
                    const int n = M.triangleAdjacency[t][leave];
                    if (n < 0 || n == start) break;
                    came = localEdgeOf(n, M.triangleEdges[t][leave]);
                    t = n;
                }
                if (fan.size() != ts.size()) fan = ts;   // non-manifold: take them as they come
            }
            double prev = valueAt(fan[0], v), sum = prev, lo = prev, hi = prev;
            for (size_t i = 1; i < fan.size(); ++i) {
                const double x = liftNear(prev, valueAt(fan[i], v));
                sum += x;
                lo = std::min(lo, x);
                hi = std::max(hi, x);
                prev = x;
            }
            regular_[v] = sum / static_cast<double>(fan.size());
            if (loopOf[v] < 0) report_.vertexSpread = std::max(report_.vertexSpread, 0.5 * (hi - lo));
        }
    }

    // ---- holes (Sec. 6.8) ----------------------------------------------------
    std::vector<double> offsetOfLoop(loops.size(), 0.0);
    for (size_t l = 0; l < loops.size(); ++l) {
        double mean = 0.0, spread = 0.0;
        loopOffsets(loops[l], mean, spread);
        offsetOfLoop[l] = mean;
        if (loopTurning[l] > 0.0) {
            report_.anchorTangency = std::max(report_.anchorTangency, std::fabs(mean) + spread);
        } else {
            report_.holeOffsets.push_back(mean);
            report_.holeSpread = std::max(report_.holeSpread, spread);
        }
    }
    if (opts_.alignHoles && report_.holes > 0) {
        // Harmonic measure of each hole: Dirichlet 1 on it and 0 on every
        // other loop of its component, harmonic in between. One factorisation
        // serves them all.
        std::vector<int> idx(nV, -1);
        int n = 0;
        for (int v = 0; v < nV; ++v) if (used[v] && loopOf[v] < 0) idx[v] = n++;
        std::vector<Eigen::Triplet<double>> red;
        red.reserve(trips.size());
        for (const auto &tr : trips)
            if (idx[tr.row()] >= 0 && idx[tr.col()] >= 0) red.emplace_back(idx[tr.row()], idx[tr.col()], tr.value());
        Eigen::SparseMatrix<double> Kii(n, n);
        Kii.setFromTriplets(red.begin(), red.end());
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
        if (n > 0) {
            ldlt.compute(Kii);
            if (ldlt.info() != Eigen::Success) return fail("the harmonic measures would not factorise");
        }
        for (size_t l = 0; l < loops.size(); ++l) {
            if (loopTurning[l] > 0.0 || offsetOfLoop[l] == 0.0) continue;
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);
            for (const auto &tr : trips)
                if (idx[tr.row()] >= 0 && idx[tr.col()] < 0 && loopOf[tr.col()] == static_cast<int>(l))
                    rhs[idx[tr.row()]] -= tr.value();
            Eigen::VectorXd omega = (n > 0) ? Eigen::VectorXd(ldlt.solve(rhs)) : Eigen::VectorXd();
            const int comp = compOfVertex[loops[l][0]];
            for (int v = 0; v < nV; ++v) {
                if (compOfVertex[v] != comp) continue;
                const double om = (idx[v] >= 0) ? omega[idx[v]] : (loopOf[v] == static_cast<int>(l) ? 1.0 : 0.0);
                regular_[v] -= offsetOfLoop[l] * om;
            }
        }
    }
    // Tangency of the per-vertex field the tracer reads, at every boundary
    // edge midpoint, on the holes.
    for (size_t l = 0; l < loops.size(); ++l) {
        if (loopTurning[l] > 0.0) continue;
        for (const int a : loops[l]) {
            const int b = nextOf[a];
            const Point m = (M.vertices[a] + M.vertices[b]) * 0.5;
            const double ra = regular_[a], rb = liftNear(ra, regular_[b]);
            const double th = 0.5 * (ra + rb) + singularPart(m);
            report_.holeTangency = std::max(
                report_.holeTangency, std::fabs(wrapQuarter(th - computeAngle(M.vertices[b] - M.vertices[a]))));
        }
    }

    // ---- the field as a CrossField -------------------------------------------
    cross_ = std::make_shared<CrossField>(mesh_);
    cross_->u_k.resize(nV);
    for (int v = 0; v < nV; ++v)
        cross_->u_k[v] = std::polar(1.0, 4.0 * (regular_[v] + singularPart(M.vertices[v])));
    cross_->u_k_prev = cross_->u_k;
    cross_->error = 0.0;

    report_.built = true;
}
