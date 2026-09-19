#include "ATLAS/CavityFill.hxx"

#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace {

inline long long pairKey(long long a, long long b) { return (a << 32) | b; }

double cornerAngleOf(const Point &p, const Point &next, const Point &prev) {
    const Point u = next - p, w = prev - p;
    double a = std::atan2(cross2(u, w), dotP(u, w));
    if (a < 0.0) a += 2.0 * M_PI;
    return a;
}

std::vector<double> arcFractions(const std::vector<Point> &pts) {
    std::vector<double> f(pts.size(), 0.0);
    for (size_t k = 1; k < pts.size(); ++k) f[k] = f[k - 1] + normP(pts[k] - pts[k - 1]);
    const double L = f.empty() ? 0.0 : f.back();
    for (size_t k = 0; k < f.size(); ++k) {
        f[k] = L > 0.0 ? f[k] / L : static_cast<double>(k) / std::max<size_t>(1, f.size() - 1);
    }
    return f;
}

} // namespace

// ---------------------------------------------------------------------------
CavityFill::Boundary CavityFill::boundaryOf(const SquareCarrier &C, const std::vector<int> &cells,
                                            const std::unordered_set<int> &inCavity) {
    Boundary B;
    std::unordered_map<int, std::pair<int, int>> next;
    next.reserve(cells.size() * 2);
    std::vector<int> starts;
    for (int q : cells) {
        for (int i = 0; i < 4; ++i) {
            const int r = C.neighbor[q][i];
            if (r >= 0 && inCavity.count(r)) continue;
            const int a = C.cells[q][i], b = C.cells[q][(i + 1) & 3];
            if (next.count(a)) {
                B.manifold = false;
                B.reason = "the cavity is pinched at a vertex";
                return B;
            }
            next.emplace(a, std::make_pair(b, C.cellEdges[q][i]));
            starts.push_back(a);
        }
    }
    std::unordered_set<int> onLoop;
    for (int s : starts) {
        if (onLoop.count(s)) continue;
        std::vector<int> loop, loopE;
        int v = s;
        size_t guard = 0;
        while (true) {
            auto it = next.find(v);
            if (it == next.end() || onLoop.count(v) || ++guard > next.size() + 1) {
                B.manifold = false;
                B.reason = "the cavity boundary does not close";
                return B;
            }
            onLoop.insert(v);
            loop.push_back(v);
            loopE.push_back(it->second.second);
            v = it->second.first;
            if (v == s) break;
        }
        B.loops.push_back(std::move(loop));
        B.loopEdges.push_back(std::move(loopE));
    }
    for (const auto &loop : B.loops) {
        std::vector<double> ang(loop.size(), 0.0);
        std::vector<int> ins(loop.size(), 0), outs(loop.size(), 0);
        double area = 0.0;
        for (size_t k = 0; k < loop.size(); ++k) {
            const int v = loop[k];
            for (int r = C.ringPtr[v]; r < C.ringPtr[v + 1]; ++r) {
                const int q = C.ringCell[r];
                if (inCavity.count(q)) {
                    ang[k] += C.cornerAngle(q, C.ringCorner[r]);
                    ++ins[k];
                } else {
                    ++outs[k];
                }
            }
            area += cross2(C.vertices[v], C.vertices[loop[(k + 1) % loop.size()]]);
        }
        B.angle.push_back(std::move(ang));
        B.inside.push_back(std::move(ins));
        B.outside.push_back(std::move(outs));
        B.signedArea.push_back(0.5 * area);
    }
    std::unordered_set<int> seen;
    for (int q : cells) {
        for (int v : C.cells[q]) {
            if (onLoop.count(v) || seen.count(v)) continue;
            seen.insert(v);
            B.interior.push_back(v);
        }
    }
    return B;
}

// ---------------------------------------------------------------------------
bool CavityFill::droppable(const SquareCarrier &C, int v, const std::unordered_set<int> &inCavity) {
    if (!C.boundaryVertex[v] || C.protectedVertex[v] || C.designatedVertex[v]) return false;
    if (C.vertexOrigin[v] != Origin::EdgeMidpoint && C.vertexOrigin[v] != Origin::BoundarySplit) return false;
    const int e = C.sourceEdge[v];
    if (e < 0 || !C.getDomain().boundaryEdge[e]) return false;
    for (int r = C.ringPtr[v]; r < C.ringPtr[v + 1]; ++r) {
        if (!inCavity.count(C.ringCell[r])) return false;
    }
    return true;
}

int CavityFill::domainEdge(const SquareCarrier &C, int a, int b) {
    const PlanarDomain &D = C.getDomain();
    for (int v : {a, b}) {
        if (v < 0) continue;
        const int e = C.sourceEdge[v];
        if (e >= 0 && D.boundaryEdge[e] && C.vertexOrigin[v] != Origin::MeshVertex) return e;
    }
    if (a < 0 || b < 0) return -1;
    const int u = C.sourceVertex[a], w = C.sourceVertex[b];
    if (u < 0 || w < 0) return -1;
    const Mesh &M = D.getMesh();
    const auto range = M.vertexTriangles.trianglesForVertex(u);
    for (const int *t = range.first; t != range.second; ++t) {
        for (int k = 0; k < 3; ++k) {
            const int e = M.triangleEdges[*t][k];
            const auto &ev = M.edges[e];
            if ((ev[0] == u && ev[1] == w) || (ev[0] == w && ev[1] == u)) return D.boundaryEdge[e] ? e : -1;
        }
    }
    return -1;
}

std::vector<int> CavityFill::insertOnBoundary(const SquareCarrier &C, Patch &P, const std::vector<int> &ids,
                                              int extra) {
    if (extra <= 0) return ids;
    struct Seg { int a, b, meshEdge; };
    std::vector<Seg> segs;
    for (size_t k = 0; k + 1 < ids.size(); ++k) {
        const int a = ids[k], b = ids[k + 1];
        int e = -1;
        if (a >= 0 && b >= 0) {
            // Only a segment that is on dS: either a carrier boundary edge,
            // or two boundary points a drop has made neighbours.
            const int ce = C.edgeBetween(a, b);
            if (ce < 0 || C.boundaryEdge[ce]) e = domainEdge(C, a, b);
        } else if (a < 0 && P.newOrigin[-1 - a] == Origin::BoundarySplit) {
            e = P.newSourceEdge[-1 - a];
        } else if (b < 0 && P.newOrigin[-1 - b] == Origin::BoundarySplit) {
            e = P.newSourceEdge[-1 - b];
        }
        segs.push_back({a, b, e});
    }
    for (int x = 0; x < extra; ++x) {
        int best = -1;
        double bestL = -1.0;
        for (size_t k = 0; k < segs.size(); ++k) {
            if (segs[k].meshEdge < 0) continue;
            const double L = normP(P.at(C, segs[k].b) - P.at(C, segs[k].a));
            if (L > bestL) { bestL = L; best = static_cast<int>(k); }
        }
        if (best < 0) return {};
        const Seg s = segs[best];
        const int m = P.addVertex((P.at(C, s.a) + P.at(C, s.b)) * 0.5, Origin::BoundarySplit, s.meshEdge);
        segs[best] = Seg{s.a, m, s.meshEdge};
        segs.insert(segs.begin() + best + 1, Seg{m, s.b, s.meshEdge});
    }
    std::vector<int> out{segs.front().a};
    for (const Seg &s : segs) out.push_back(s.b);
    return out;
}

std::vector<int> CavityFill::dropFromBoundary(const SquareCarrier &C, const std::vector<int> &ids,
                                              const std::vector<char> &mayDrop, int count) {
    std::vector<int> out = ids;
    std::vector<char> may = mayDrop;
    for (int x = 0; x < count; ++x) {
        int best = -1;
        double bestL = 1e300;
        for (size_t k = 1; k + 1 < out.size(); ++k) {
            if (!may[k]) continue;
            const double L = normP(C.vertices[out[k]] - C.vertices[out[k - 1]]) +
                             normP(C.vertices[out[k + 1]] - C.vertices[out[k]]);
            if (L < bestL) { bestL = L; best = static_cast<int>(k); }
        }
        if (best < 0) return {};
        out.erase(out.begin() + best);
        may.erase(may.begin() + best);
    }
    return out;
}

// ---------------------------------------------------------------------------
bool CavityFill::grid(const SquareCarrier &C, Patch &P, const std::array<std::vector<int>, 4> &sides,
                      Origin interiorOrigin, std::vector<int> *nodesOut) {
    const int n = static_cast<int>(sides[0].size()) - 1;
    const int m = static_cast<int>(sides[1].size()) - 1;
    if (n < 1 || m < 1) return false;
    if (static_cast<int>(sides[2].size()) != n + 1 || static_cast<int>(sides[3].size()) != m + 1) return false;
    if (sides[0].back() != sides[1].front() || sides[1].back() != sides[2].front() ||
        sides[2].back() != sides[3].front() || sides[3].back() != sides[0].front()) {
        return false;
    }
    const int W = n + 1;
    std::vector<int> g(static_cast<size_t>(W) * (m + 1), 0);
    auto at = [&](int i, int j) -> int & { return g[i + W * j]; };
    for (int i = 0; i <= n; ++i) { at(i, 0) = sides[0][i]; at(i, m) = sides[2][n - i]; }
    for (int j = 0; j <= m; ++j) { at(n, j) = sides[1][j]; at(0, j) = sides[3][m - j]; }

    std::vector<Point> Bp(n + 1), Tp(n + 1), Lp(m + 1), Rp(m + 1);
    for (int i = 0; i <= n; ++i) { Bp[i] = P.at(C, at(i, 0)); Tp[i] = P.at(C, at(i, m)); }
    for (int j = 0; j <= m; ++j) { Lp[j] = P.at(C, at(0, j)); Rp[j] = P.at(C, at(n, j)); }
    const std::vector<double> fb = arcFractions(Bp), ft = arcFractions(Tp);
    const std::vector<double> fl = arcFractions(Lp), fr = arcFractions(Rp);
    const Point P00 = Bp[0], P10 = Bp[n], P11 = Tp[n], P01 = Tp[0];

    for (int j = 1; j < m; ++j) {
        for (int i = 1; i < n; ++i) {
            const double ub = fb[i], ut = ft[i], vl = fl[j], vr = fr[j];
            const double den = 1.0 - (ut - ub) * (vr - vl);
            const double u = (ub + vl * (ut - ub)) / den;
            const double v = vl + u * (vr - vl);
            const Point X = Bp[i] * (1.0 - v) + Tp[i] * v + Lp[j] * (1.0 - u) + Rp[j] * u -
                            (P00 * ((1.0 - u) * (1.0 - v)) + P10 * (u * (1.0 - v)) +
                             P11 * (u * v) + P01 * ((1.0 - u) * v));
            at(i, j) = P.addVertex(X, interiorOrigin);
        }
    }
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < n; ++i) {
            P.cells.push_back({at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)});
        }
    }
    if (nodesOut) *nodesOut = std::move(g);
    return true;
}

std::vector<int> CavityFill::segment(const SquareCarrier &C, Patch &P, int a, int b, int n, Origin o) {
    std::vector<int> ids{a};
    const Point pa = P.at(C, a), pb = P.at(C, b);
    for (int k = 1; k < n; ++k) ids.push_back(P.addVertex(pa + (pb - pa) * (static_cast<double>(k) / n), o));
    ids.push_back(b);
    return ids;
}

bool CavityFill::star(const SquareCarrier &C, Patch &P, const std::vector<std::vector<int>> &arcs,
                      const std::vector<int> &sigma, const Point &centre, Origin interiorOrigin,
                      std::vector<int> *splitPoints, int *centreId) {
    const int K = static_cast<int>(arcs.size());
    if (K < 3 || static_cast<int>(sigma.size()) != K) return false;
    for (int k = 0; k < K; ++k) {
        if (static_cast<int>(arcs[k].size()) - 1 != sigma[(k + K - 1) % K] + sigma[(k + 1) % K]) return false;
    }
    const int c = P.addVertex(centre, interiorOrigin);
    // m_k sits x_k = s_{k-1} edges along arc k.
    std::vector<int> mid(K);
    for (int k = 0; k < K; ++k) mid[k] = arcs[k][sigma[(k + K - 1) % K]];
    std::vector<std::vector<int>> spoke(K);
    for (int k = 0; k < K; ++k) spoke[k] = segment(C, P, mid[k], c, sigma[k], interiorOrigin);
    for (int k = 0; k < K; ++k) {
        const int km = (k + K - 1) % K;
        const int xPrev = sigma[(km + K - 1) % K];
        std::array<std::vector<int>, 4> sides;
        sides[0].assign(arcs[km].begin() + xPrev, arcs[km].end());           // m_{k-1} -> corner k
        sides[1].assign(arcs[k].begin(), arcs[k].begin() + sigma[km] + 1);   // corner k -> m_k
        sides[2] = spoke[k];                                                 // m_k -> centre
        sides[3] = spoke[km];
        std::reverse(sides[3].begin(), sides[3].end());                      // centre -> m_{k-1}
        if (!grid(C, P, sides, interiorOrigin)) return false;
    }
    if (splitPoints) *splitPoints = mid;
    if (centreId) *centreId = c;
    return true;
}

bool CavityFill::solveStar(const std::vector<int> &n, std::vector<int> &s) {
    const int K = static_cast<int>(n.size());
    if (K % 2 == 0 || K < 3) return false;
    long long S = 0;
    for (int x : n) S += x;
    if (S % 2 != 0) return false;
    s.assign(K, 0);
    if (K == 3) {
        for (int i = 0; i < 3; ++i) {
            s[i] = static_cast<int>(S / 2 - n[i]);
            if (s[i] < 1) return false;
        }
        return true;
    }
    // Gaussian elimination on the K x K circulant, then an exact integer check.
    std::vector<std::vector<double>> M(K, std::vector<double>(K + 1, 0.0));
    for (int i = 0; i < K; ++i) {
        M[i][(i + K - 1) % K] += 1.0;
        M[i][(i + 1) % K] += 1.0;
        M[i][K] = n[i];
    }
    for (int c = 0; c < K; ++c) {
        int piv = c;
        for (int r = c + 1; r < K; ++r) if (std::fabs(M[r][c]) > std::fabs(M[piv][c])) piv = r;
        if (std::fabs(M[piv][c]) < 1e-12) return false;
        std::swap(M[c], M[piv]);
        for (int r = 0; r < K; ++r) {
            if (r == c) continue;
            const double f = M[r][c] / M[c][c];
            for (int k = c; k <= K; ++k) M[r][k] -= f * M[c][k];
        }
    }
    for (int i = 0; i < K; ++i) {
        s[i] = static_cast<int>(std::lround(M[i][K] / M[i][i]));
        if (s[i] < 1) return false;
    }
    for (int i = 0; i < K; ++i) {
        if (s[(i + K - 1) % K] + s[(i + 1) % K] != n[i]) return false;
    }
    return true;
}

bool CavityFill::repairStar(const std::vector<int> &n, const std::vector<char> &splittable,
                            std::vector<int> &s) {
    const int K = static_cast<int>(n.size());
    if (K != 3 && K != 5) return false;
    // The real solution, as the tie-break target.
    std::vector<double> ideal(K, 1.0);
    {
        std::vector<std::vector<double>> M(K, std::vector<double>(K + 1, 0.0));
        for (int i = 0; i < K; ++i) {
            M[i][(i + K - 1) % K] += 1.0;
            M[i][(i + 1) % K] += 1.0;
            M[i][K] = n[i];
        }
        for (int c = 0; c < K; ++c) {
            int piv = c;
            for (int r = c + 1; r < K; ++r) if (std::fabs(M[r][c]) > std::fabs(M[piv][c])) piv = r;
            std::swap(M[c], M[piv]);
            for (int r = 0; r < K; ++r) {
                if (r == c) continue;
                const double f = M[r][c] / M[c][c];
                for (int k = c; k <= K; ++k) M[r][k] -= f * M[c][k];
            }
        }
        for (int i = 0; i < K; ++i) ideal[i] = std::max(1.0, M[i][K] / M[i][i]);
    }
    int maxN = 1;
    for (int x : n) maxN = std::max(maxN, x);

    // A side's constraint on the two spokes either side of it.
    auto lower = [&](int side, int other) { return n[side] - other; };
    long long bestD = -1;
    double bestDev = 0.0;
    std::vector<int> cur(K);
    auto consider = [&]() {
        long long D = 0;
        double dev = 0.0;
        for (int i = 0; i < K; ++i) {
            if (cur[i] < 1) return;
            const int have = cur[(i + K - 1) % K] + cur[(i + 1) % K];
            if (have < n[i] || (!splittable[i] && have != n[i])) return;
            D += have - n[i];
            dev += (cur[i] - ideal[i]) * (cur[i] - ideal[i]);
        }
        if (bestD < 0 || D < bestD || (D == bestD && dev < bestDev)) {
            bestD = D;
            bestDev = dev;
            s = cur;
        }
    };
    // Given the enumerated spokes, the remaining one or two are each bounded
    // below by the sides they touch; take the least value that satisfies
    // them, exactly where a side is fixed.
    auto settleSpoke = [&](int j) {
        // Sides touching spoke j: j - 1 (with spoke j - 2) and j + 1 (with spoke j + 2).
        const int sa = (j + K - 1) % K, oa = (j + K - 2) % K;
        const int sb = (j + 1) % K, ob = (j + 2) % K;
        int v = std::max({1, lower(sa, cur[oa]), lower(sb, cur[ob])});
        if (!splittable[sa]) v = lower(sa, cur[oa]);
        if (!splittable[sb]) {
            const int w = lower(sb, cur[ob]);
            if (!splittable[sa] && w != v) return false;
            v = w;
        }
        cur[j] = v;
        return v >= 1;
    };
    if (K == 3) {
        for (cur[0] = 1; cur[0] <= maxN; ++cur[0]) {
            for (cur[1] = 1; cur[1] <= maxN; ++cur[1]) {
                // Spoke 2 touches sides 0 (with spoke 1) and 1 (with spoke 0).
                int v = std::max({1, lower(0, cur[1]), lower(1, cur[0])});
                if (!splittable[0]) v = lower(0, cur[1]);
                if (!splittable[1]) {
                    if (!splittable[0] && lower(1, cur[0]) != v) continue;
                    v = lower(1, cur[0]);
                }
                cur[2] = v;
                consider();
            }
        }
    } else {
        // Spokes 3 and 4 only meet sides that also hold one of spokes 0..2,
        // so they settle independently once those three are fixed.
        const int W = 60;
        int lo[3], hi[3];
        for (int k = 0; k < 3; ++k) {
            const int c = static_cast<int>(std::lround(ideal[k]));
            lo[k] = std::max(1, c - W);
            hi[k] = std::min(maxN, c + W);
        }
        for (cur[0] = lo[0]; cur[0] <= hi[0]; ++cur[0]) {
            for (cur[1] = lo[1]; cur[1] <= hi[1]; ++cur[1]) {
                for (cur[2] = lo[2]; cur[2] <= hi[2]; ++cur[2]) {
                    if (!settleSpoke(3) || !settleSpoke(4)) continue;
                    consider();
                }
            }
        }
    }
    return bestD >= 0;
}

namespace {

// The integer x in [L, H], x >= 1, minimising |x - a| + |x - b|, nearest c
// among the minimisers. False when [L, H] has no integer >= 1.
bool bestSpoke(long long L, long long H, long long a, long long b, double c, long long &x) {
    L = std::max(1LL, L);
    if (L > H) return false;
    const long long p = std::min(a, b), q = std::max(a, b);
    long long lo = std::max(L, p), hi = std::min(H, q);
    if (lo > hi) {
        x = (H < p) ? H : L;
        return true;
    }
    x = std::max(lo, std::min(hi, static_cast<long long>(std::llround(c))));
    return true;
}

} // namespace

bool CavityFill::rangeStar(const std::vector<int> &n, const std::vector<int> &lo, const std::vector<int> &hi,
                           std::vector<int> &s, std::vector<int> &have, int window) {
    const int K = static_cast<int>(n.size());
    if (K != 3 && K != 5) return false;
    std::vector<double> ideal(K, 1.0);
    {
        std::vector<std::vector<double>> M(K, std::vector<double>(K + 1, 0.0));
        for (int i = 0; i < K; ++i) {
            M[i][(i + K - 1) % K] += 1.0;
            M[i][(i + 1) % K] += 1.0;
            M[i][K] = n[i];
        }
        for (int c = 0; c < K; ++c) {
            int piv = c;
            for (int r = c + 1; r < K; ++r) if (std::fabs(M[r][c]) > std::fabs(M[piv][c])) piv = r;
            std::swap(M[c], M[piv]);
            for (int r = 0; r < K; ++r) {
                if (r == c) continue;
                const double f = M[r][c] / M[c][c];
                for (int k = c; k <= K; ++k) M[r][k] -= f * M[c][k];
            }
        }
        for (int i = 0; i < K; ++i) ideal[i] = std::max(1.0, M[i][K] / M[i][i]);
    }
    int maxH = 1;
    for (int x : hi) maxH = std::max(maxH, x);
    long long bestD = -1;
    double bestDev = 0.0;
    std::vector<long long> cur(K, 1);
    auto consider = [&]() {
        long long D = 0;
        double dev = 0.0;
        for (int i = 0; i < K; ++i) {
            if (cur[i] < 1) return;
            const long long h = cur[(i + K - 1) % K] + cur[(i + 1) % K];
            if (h < lo[i] || h > hi[i]) return;
            D += std::llabs(h - n[i]);
            dev += (cur[i] - ideal[i]) * (cur[i] - ideal[i]);
        }
        if (bestD < 0 || D < bestD || (D == bestD && dev < bestDev)) {
            bestD = D;
            bestDev = dev;
            s.assign(cur.begin(), cur.end());
        }
    };
    // Spoke j, given its two neighbours' partners: side j-1 = s_{j-2} + s_j and
    // side j+1 = s_j + s_{j+2}.
    auto settle = [&](int j) {
        const int sa = (j + K - 1) % K, oa = (j + K - 2) % K;
        const int sb = (j + 1) % K, ob = (j + 2) % K;
        const long long L = std::max(static_cast<long long>(lo[sa]) - cur[oa], static_cast<long long>(lo[sb]) - cur[ob]);
        const long long H = std::min(static_cast<long long>(hi[sa]) - cur[oa], static_cast<long long>(hi[sb]) - cur[ob]);
        long long x = 0;
        if (!bestSpoke(L, H, n[sa] - cur[oa], n[sb] - cur[ob], ideal[j], x)) return false;
        cur[j] = x;
        return true;
    };
    // The enumerated spokes range over a window round the real solution;
    // the rest settle exactly.
    long long a[3], b[3];
    for (int k = 0; k < 3; ++k) {
        const long long c = std::llround(ideal[k]);
        a[k] = std::max(1LL, c - window);
        b[k] = std::min(static_cast<long long>(maxH), c + window);
    }
    if (K == 3) {
        for (cur[0] = a[0]; cur[0] <= b[0]; ++cur[0]) {
            for (cur[1] = a[1]; cur[1] <= b[1]; ++cur[1]) {
                if (settle(2)) consider();
            }
        }
    } else {
        for (cur[0] = a[0]; cur[0] <= b[0]; ++cur[0]) {
            for (cur[1] = a[1]; cur[1] <= b[1]; ++cur[1]) {
                for (cur[2] = a[2]; cur[2] <= b[2]; ++cur[2]) {
                    if (settle(3) && settle(4)) consider();
                }
            }
        }
    }
    if (bestD < 0) return false;
    have.assign(K, 0);
    for (int i = 0; i < K; ++i) have[i] = s[(i + K - 1) % K] + s[(i + 1) % K];
    return true;
}

// ---------------------------------------------------------------------------
void CavityFill::smooth(const SquareCarrier &C, Patch &P, int iterations) {
    if (iterations <= 0 || P.newVertices.empty()) return;
    const int K = static_cast<int>(P.newVertices.size());
    // Neighbour lists of the new vertices only; carrier vertices never move.
    std::vector<std::vector<int>> nbr(K);
    for (const auto &q : P.cells) {
        for (int c = 0; c < 4; ++c) {
            const int a = q[c], b = q[(c + 1) & 3];
            if (a < 0) nbr[-1 - a].push_back(b);
            if (b < 0) nbr[-1 - b].push_back(a);
        }
    }
    for (auto &l : nbr) { std::sort(l.begin(), l.end()); l.erase(std::unique(l.begin(), l.end()), l.end()); }
    std::vector<Point> next(K);
    for (int it = 0; it < iterations; ++it) {
        for (int k = 0; k < K; ++k) {
            next[k] = P.newVertices[k];
            if (P.newOrigin[k] == Origin::BoundarySplit || nbr[k].empty()) continue;
            Point s{0.0, 0.0};
            for (int id : nbr[k]) s = s + P.at(C, id);
            next[k] = s / static_cast<double>(nbr[k].size());
        }
        P.newVertices = next;
    }
}

CavityFill::Verdict CavityFill::settle(const SquareCarrier &C, const std::vector<int> &cavity, Patch &P,
                                       double minScaledJacobian, int smoothingIterations) {
    Verdict best = validate(C, cavity, P, minScaledJacobian);
    // Good enough already: an explicit map that certifies does not get
    // smoothed away from itself.
    if (smoothingIterations <= 0 || (best.valid && best.minScaledJacobian >= 0.5)) return best;
    Patch work = P;
    Patch bestPatch = P;
    const int chunk = 5;
    for (int done = 0; done < smoothingIterations; done += chunk) {
        smooth(C, work, chunk);
        const Verdict v = validate(C, cavity, work, minScaledJacobian);
        if (v.valid && (!best.valid || v.minScaledJacobian > best.minScaledJacobian + 1e-12)) {
            best = v;
            bestPatch = work;
        }
    }
    if (best.valid) P = bestPatch;
    return best;
}

// ---------------------------------------------------------------------------
// Sec. 9.2, locally. See the class comment for why these tests suffice.
// ---------------------------------------------------------------------------
CavityFill::Verdict CavityFill::validate(const SquareCarrier &C, const std::vector<int> &cavity,
                                         const Patch &P, double minScaledJacobian) {
    Verdict V;
    auto fail = [&](const std::string &why) { V.valid = false; V.reason = why; return V; };
    if (P.cells.empty()) return fail("the patch has no cells");

    const int NV = C.numVertices();
    auto key = [&](int id) -> long long { return id >= 0 ? id : NV + (-1 - id); };

    std::unordered_set<int> inCav(cavity.begin(), cavity.end());

    // The cavity's boundary half-edges (cavity on the left) and the interior
    // angle it has at each boundary vertex.
    struct Half { int a, b; bool domain; };
    std::vector<Half> oldHalf;
    std::unordered_map<int, double> oldAngle;
    double oldArea = 0.0;
    for (int q : cavity) {
        oldArea += C.cellArea(q);
        for (int i = 0; i < 4; ++i) {
            const int r = C.neighbor[q][i];
            if (r >= 0 && inCav.count(r)) continue;
            oldHalf.push_back({C.cells[q][i], C.cells[q][(i + 1) & 3], r < 0});
        }
    }
    for (const Half &h : oldHalf) oldAngle.emplace(h.a, 0.0);
    for (auto &kv : oldAngle) {
        const int v = kv.first;
        for (int r = C.ringPtr[v]; r < C.ringPtr[v + 1]; ++r) {
            if (inCav.count(C.ringCell[r])) kv.second += C.cornerAngle(C.ringCell[r], C.ringCorner[r]);
        }
    }

    // Edge census of the patch: every undirected edge in at most two cells, and
    // in opposite directions when in two.
    std::unordered_map<long long, int> uses;   // undirected -> count
    std::unordered_map<long long, int> dirUse; // directed -> count
    uses.reserve(P.cells.size() * 4);
    dirUse.reserve(P.cells.size() * 4);
    for (const auto &q : P.cells) {
        for (int c = 0; c < 4; ++c) {
            const long long a = key(q[c]), b = key(q[(c + 1) & 3]);
            if (a == b) return fail("a cell repeats a vertex");
            ++uses[pairKey(std::min(a, b), std::max(a, b))];
            if (++dirUse[pairKey(a, b)] > 1) return fail("two cells run an edge the same way");
        }
    }
    std::unordered_map<long long, int> nextOf;   // patch boundary: start -> end
    int patchBoundary = 0;
    for (const auto &q : P.cells) {
        for (int c = 0; c < 4; ++c) {
            const long long a = key(q[c]), b = key(q[(c + 1) & 3]);
            const int u = uses[pairKey(std::min(a, b), std::max(a, b))];
            if (u > 2) return fail("an edge of the patch is in three cells");
            if (u == 1) {
                if (nextOf.count(a)) return fail("the patch boundary is pinched");
                nextOf[a] = q[(c + 1) & 3];
                ++patchBoundary;
            }
        }
    }

    // Match the cavity's boundary. Between two consecutive vertices the
    // patch must keep (anchors), the old boundary is either one edge -- kept
    // as it is, or, on dS, subdivided by inserted points -- or a run of
    // droppable points along one domain segment, which the patch may walk
    // with any in-order subset of them plus inserted points, all collinear.
    std::unordered_map<int, std::pair<int, bool>> oldNext;
    for (const Half &h : oldHalf) oldNext.emplace(h.a, std::make_pair(h.b, h.domain));
    std::unordered_map<int, char> drop;
    for (const Half &h : oldHalf) drop.emplace(h.a, droppable(C, h.a, inCav) ? 1 : 0);
    int consumed = 0;
    for (const Half &h0 : oldHalf) {
        const int a = h0.a;
        if (drop[a]) continue;
        // The old run from a to the next anchor b.
        std::unordered_set<int> run;
        bool allDomain = true;
        int b = a;
        int guard = 0;
        do {
            const auto &nx = oldNext[b];
            allDomain = allDomain && nx.second;
            b = nx.first;
            if (drop[b]) run.insert(b);
            if (++guard > static_cast<int>(oldHalf.size()) + 1) return fail("the cavity boundary does not close");
        } while (drop[b]);
        auto it = nextOf.find(a);
        if (it == nextOf.end()) return fail("a cavity boundary vertex is not on the patch boundary");
        int w = it->second;
        ++consumed;
        const Point pa = C.vertices[a], pb = C.vertices[b];
        const Point d = pb - pa;
        const double L2 = dotP(d, d);
        double lastT = 0.0;
        guard = 0;
        while (w != b) {
            if (!allDomain) return fail("a point was inserted on an edge shared with a neighbour (hanging node)");
            if (w >= 0) {
                if (!run.count(w)) return fail("the patch boundary leaves the cavity boundary");
            } else if (P.newOrigin[-1 - w] != Origin::BoundarySplit) {
                return fail("an inserted boundary point is not marked as one");
            }
            const Point pw = P.at(C, w);
            const double t = dotP(pw - pa, d) / L2;
            const double off = std::fabs(cross2(pw - pa, d)) / std::sqrt(L2);
            if (!(t > lastT) || !(t < 1.0) || off > 1e-9 * std::sqrt(L2)) {
                return fail("an inserted boundary point is off its segment or out of order");
            }
            lastT = t;
            auto nt = nextOf.find(key(w));
            if (nt == nextOf.end() || ++guard > 100000) return fail("the patch boundary breaks off");
            w = nt->second;
            ++consumed;
        }
    }
    if (consumed != patchBoundary) return fail("the patch has boundary edges the cavity does not");

    // Only boundary vertices of the cavity may be reused; its interior
    // vertices are deleted with it.
    for (const auto &q : P.cells) {
        for (int v : q) if (v >= 0 && !oldAngle.count(v)) return fail("the patch uses a vertex outside the cavity boundary");
    }

    // Cells: Sec. 9.1's four corners, plus the requested quality floor.
    double sjMin = 1.0, sjSum = 0.0, newArea = 0.0;
    std::unordered_map<long long, double> angleSum;
    for (const auto &q : P.cells) {
        double cellMin = 1.0;
        for (int c = 0; c < 4; ++c) {
            const Point p = P.at(C, q[c]);
            const Point nx = P.at(C, q[(c + 1) & 3]), pv = P.at(C, q[(c + 3) & 3]);
            const Point u = nx - p, w = pv - p;
            const double det = cross2(u, w);
            const double nn = normP(u) * normP(w);
            if (!(det > 0.0) || !(nn > 0.0)) return fail("a new cell has a non-positive corner Jacobian");
            cellMin = std::min(cellMin, det / nn);
            angleSum[key(q[c])] += cornerAngleOf(p, nx, pv);
            newArea += 0.5 * cross2(p, nx);
        }
        sjMin = std::min(sjMin, cellMin);
        sjSum += cellMin;
    }
    V.minScaledJacobian = sjMin;
    V.meanScaledJacobian = sjSum / P.cells.size();
    if (sjMin < minScaledJacobian) {
        V.reason = "a new cell is below the scaled-Jacobian floor";
        return V;
    }

    // One-rings: 2 pi inside, the cavity's own angle on its boundary, pi at an
    // inserted boundary point.
    const double tol = 1e-7;
    for (const auto &kv : angleSum) {
        const long long k = kv.first;
        double want;
        if (k < NV) {
            want = oldAngle[static_cast<int>(k)];
        } else {
            const int idx = static_cast<int>(k - NV);
            want = P.newOrigin[idx] == Origin::BoundarySplit ? M_PI : 2.0 * M_PI;
        }
        if (std::fabs(kv.second - want) > tol) return fail("a one-ring does not wind exactly once");
    }
    for (int k = 0; k < static_cast<int>(P.newVertices.size()); ++k) {
        if (!angleSum.count(NV + k)) return fail("a new vertex is used by no cell");
    }
    if (std::fabs(newArea - oldArea) > 1e-9 * std::max(1e-300, std::fabs(oldArea))) {
        return fail("the patch does not have the cavity's area");
    }
    V.valid = true;
    V.reason.clear();
    return V;
}

// ---------------------------------------------------------------------------
double CavityFill::signedArea(const std::vector<Point> &loop) {
    double A = 0.0;
    for (size_t k = 0; k < loop.size(); ++k) A += cross2(loop[k], loop[(k + 1) % loop.size()]);
    return 0.5 * A;
}

Point CavityFill::areaCentroid(const std::vector<Point> &loop) {
    double A = 0.0;
    Point c{0.0, 0.0};
    for (size_t k = 0; k < loop.size(); ++k) {
        const Point &p = loop[k], &q = loop[(k + 1) % loop.size()];
        const double w = cross2(p, q);
        A += w;
        c = c + (p + q) * w;
    }
    if (std::fabs(A) < 1e-300) {
        Point s{0.0, 0.0};
        for (const Point &p : loop) s = s + p;
        return loop.empty() ? s : s / static_cast<double>(loop.size());
    }
    return c / (3.0 * A);
}

double CavityFill::kernelDepth(const std::vector<Point> &loop, const Point &c) {
    double depth = 1e300;
    const size_t n = loop.size();
    for (size_t k = 0; k < n; ++k) {
        const Point &a = loop[k], &b = loop[(k + 1) % n];
        const Point d = b - a;
        const double L = normP(d);
        if (L <= 0.0) continue;
        depth = std::min(depth, cross2(d, c - a) / L);
    }
    return depth;
}

// Ascent on the concave function c -> min_e dist(c, line_e): step along the
// inward normal of the edge that currently attains the minimum, growing the
// step on success and shrinking it on failure. Not the exact Chebyshev centre
// of the kernel, but it only has to find *a* point well inside it.
Point CavityFill::kernelCenter(const std::vector<Point> &loop, const Point &start) {
    if (loop.size() < 3) return start;
    Point lo = loop[0], hi = loop[0];
    for (const Point &p : loop) {
        lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
        hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
    }
    const double scale = normP(hi - lo);
    double step = 0.05 * scale;
    Point c = start;
    double depth = kernelDepth(loop, c);
    const size_t n = loop.size();
    for (int it = 0; it < 400 && step > 1e-10 * scale; ++it) {
        // Inward normal of the minimising edge.
        size_t best = 0;
        double bestD = 1e300;
        for (size_t k = 0; k < n; ++k) {
            const Point &a = loop[k], &b = loop[(k + 1) % n];
            const Point d = b - a;
            const double L = normP(d);
            if (L <= 0.0) continue;
            const double dist = cross2(d, c - a) / L;
            if (dist < bestD) { bestD = dist; best = k; }
        }
        const Point d = loop[(best + 1) % n] - loop[best];
        const Point nrm = normalizeP(Point{-d[1], d[0]});
        const Point trial = c + nrm * step;
        const double td = kernelDepth(loop, trial);
        if (td > depth) { c = trial; depth = td; step *= 1.3; }
        else step *= 0.5;
    }
    return c;
}
