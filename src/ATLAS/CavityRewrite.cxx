#include "ATLAS/CavityRewrite.hxx"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <unordered_set>

namespace {

typedef CavityFill::Patch Patch;
typedef SquareCarrier::Origin Origin;

} // namespace

CavityRewrite::CavityRewrite(SquareCarrier &carrier, const Options &opts) : C_(carrier), opts_(opts) {
    report_.defectBefore = report_.defectAfter = totalDefect();
    report_.irregularBefore = report_.irregularAfter = irregularCount();
    report_.cellsBefore = report_.cellsAfter = C_.numCells();
}

int CavityRewrite::totalDefect() const {
    int d = 0;
    for (int v = 0; v < C_.numVertices(); ++v) d += C_.defect(v);
    return d;
}

int CavityRewrite::irregularCount() const {
    int n = 0;
    for (int v = 0; v < C_.numVertices(); ++v) if (C_.defect(v) != 0) ++n;
    return n;
}

// ---------------------------------------------------------------------------
// Corner sets for one cavity loop, scored by the defect they leave on it.
// ---------------------------------------------------------------------------
void CavityRewrite::propose(const std::vector<int> &loop, const std::vector<double> &angle,
                            const std::vector<int> &outside, std::vector<Proposal> &out) {
    out.clear();
    const int N = static_cast<int>(loop.size());
    const double maxCorner = opts_.maxCornerAngle * M_PI / 180.0;
    std::vector<int> delta(N);
    std::vector<double> cornerGeo(N), sideGeo(N);
    int base = 0;
    double geoBase = 0.0;
    std::vector<char> okCorner(N);
    for (int u = 0; u < N; ++u) {
        const int t = C_.targetValence(loop[u]);
        const int f1 = std::abs(outside[u] + 1 - t);
        const int f2 = std::abs(outside[u] + 2 - t);
        base += f2;
        delta[u] = f1 - f2;
        cornerGeo[u] = (angle[u] - M_PI_2) * (angle[u] - M_PI_2);
        sideGeo[u] = (angle[u] - M_PI) * (angle[u] - M_PI);
        geoBase += sideGeo[u];
        okCorner[u] = angle[u] < maxCorner;
    }

    if (opts_.grids && N >= 4 && N % 2 == 0) {
        const int half = N / 2;
        for (int i0 = 0; i0 < half; ++i0) {
            for (int a = 1; a < half; ++a) {
                const int c[4] = {i0, (i0 + a) % N, (i0 + half) % N, (i0 + half + a) % N};
                if (!okCorner[c[0]] || !okCorner[c[1]] || !okCorner[c[2]] || !okCorner[c[3]]) continue;
                Proposal p;
                p.kind = 0;
                p.corners.assign(c, c + 4);
                p.defect = base;
                p.geo = geoBase;
                for (int k = 0; k < 4; ++k) {
                    p.defect += delta[c[k]];
                    p.geo += cornerGeo[c[k]] - sideGeo[c[k]];
                }
                out.push_back(std::move(p));
            }
        }
    }

    if (opts_.stars && N >= 6) {
        // A pool of the likeliest corners, in loop order.
        std::vector<int> pool;
        for (int u = 0; u < N; ++u) if (okCorner[u]) pool.push_back(u);
        std::stable_sort(pool.begin(), pool.end(), [&](int x, int y) {
            if (delta[x] != delta[y]) return delta[x] < delta[y];
            return cornerGeo[x] < cornerGeo[y];
        });
        if (pool.size() > 9) pool.resize(9);
        std::sort(pool.begin(), pool.end());
        const int P = static_cast<int>(pool.size());
        for (int K : {3, 5}) {
            if (P < K) continue;
            std::vector<int> pick(K);
            std::function<void(int, int)> rec = [&](int start, int depth) {
                if (depth == K) {
                    std::vector<int> n(K), s;
                    for (int k = 0; k < K; ++k) {
                        n[k] = ((pick[(k + 1) % K] - pick[k]) % N + N) % N;
                        if (n[k] == 0) return;
                    }
                    if (!CavityFill::solveStar(n, s)) return;
                    Proposal p;
                    p.kind = 1;
                    p.corners = pick;
                    p.sigma = s;
                    p.defect = base + 1;   // the centre has valence K, one off
                    p.geo = geoBase;
                    for (int c : pick) {
                        p.defect += delta[c];
                        p.geo += cornerGeo[c] - sideGeo[c];
                    }
                    out.push_back(std::move(p));
                    return;
                }
                for (int i = start; i < P; ++i) {
                    pick[depth] = pool[i];
                    rec(i + 1, depth + 1);
                }
            };
            rec(0, 0);
        }
    }
    report_.proposals += static_cast<int>(out.size());
}

bool CavityRewrite::realise(const std::vector<int> &loop, const Proposal &p, Patch &P) const {
    const int N = static_cast<int>(loop.size());
    const int K = static_cast<int>(p.corners.size());
    std::vector<std::vector<int>> arcs(K);
    for (int k = 0; k < K; ++k) {
        int i = p.corners[k];
        const int e = p.corners[(k + 1) % K];
        arcs[k].push_back(loop[i]);
        do { i = (i + 1) % N; arcs[k].push_back(loop[i]); } while (i != e);
    }
    if (p.kind == 0) {
        std::array<std::vector<int>, 4> sides{arcs[0], arcs[1], arcs[2], arcs[3]};
        P.kind = "grid";
        P.blocks = 1;
        return CavityFill::grid(C_, P, sides, Origin::Rewrite);
    }
    std::vector<Point> poly;
    for (int v : loop) poly.push_back(C_.vertices[v]);
    Point c = CavityFill::kernelCenter(poly, CavityFill::areaCentroid(poly));
    if (!(CavityFill::kernelDepth(poly, c) > 0.0)) return false;
    P.kind = "star";
    P.blocks = K;
    return CavityFill::star(C_, P, arcs, p.sigma, c, Origin::Rewrite);
}

// ---------------------------------------------------------------------------
// One round of Stage 6.
// ---------------------------------------------------------------------------
int CavityRewrite::round(const std::vector<int> &priority) {
    const auto t0 = std::chrono::steady_clock::now();
    ++report_.rounds;
    const int NV = C_.numVertices(), NQ = C_.numCells();

    // Witnesses: Stage 4's first, then every vertex of wrong valence, worst
    // first (Sec. 8.3, "prioritize obstruction witnesses").
    std::vector<int> order;
    std::vector<char> queued(NV, 0);
    for (int v : priority) {
        if (v >= 0 && v < NV && !queued[v] && C_.valence[v] > 0) { queued[v] = 1; order.push_back(v); }
    }
    std::vector<int> rest;
    for (int v = 0; v < NV; ++v) if (!queued[v] && C_.defect(v) > 0) rest.push_back(v);
    std::stable_sort(rest.begin(), rest.end(), [&](int a, int b) { return C_.defect(a) > C_.defect(b); });
    order.insert(order.end(), rest.begin(), rest.end());

    std::vector<char> lockedV(NV, 0), lockedQ(NQ, 0);
    std::vector<int> stampQ(NQ, 0), stampV(NV, 0);
    int stamp = 0;
    std::vector<SquareCarrier::Edit> edits;

    struct Option {
        double dE;
        double geo;
        int k;
        Proposal p;
    };

    for (int v : order) {
        const double elapsed =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (elapsed > opts_.timeBudget) { report_.timedOut = true; break; }
        if (lockedV[v] || C_.valence[v] == 0) continue;
        ++report_.witnesses;

        // The material of the cavity: the first non-template cell at v.
        int material = -1;
        for (int r = C_.ringPtr[v]; r < C_.ringPtr[v + 1]; ++r) {
            const int q = C_.ringCell[r];
            if (C_.cellOrigin[q] != SquareCarrier::CellOrigin::Template) { material = C_.cellMaterial[q]; break; }
        }
        if (material < 0 && C_.valence[v] > 0) { reject("witness inside a template"); continue; }

        // Grow the cavity ring by ring, keeping each ring's loop and proposals.
        ++stamp;
        std::vector<int> cav;
        std::vector<int> frontier{v};
        stampV[v] = stamp;
        std::vector<std::vector<int>> cavities, loops;
        std::vector<Option> options;
        bool blocked = false;
        for (int k = 1; k <= opts_.maxRings && !blocked; ++k) {
            std::vector<int> next;
            for (int u : frontier) {
                for (int r = C_.ringPtr[u]; r < C_.ringPtr[u + 1]; ++r) {
                    const int q = C_.ringCell[r];
                    if (stampQ[q] == stamp) continue;
                    if (C_.cellMaterial[q] != material) continue;
                    if (C_.cellOrigin[q] == SquareCarrier::CellOrigin::Template) continue;
                    if (lockedQ[q]) { blocked = true; continue; }
                    stampQ[q] = stamp;
                    cav.push_back(q);
                    for (int w : C_.cells[q]) {
                        if (stampV[w] != stamp) { stampV[w] = stamp; next.push_back(w); }
                    }
                }
            }
            frontier.swap(next);
            if (static_cast<int>(cav.size()) > opts_.maxCavityCells) break;
            ++report_.cavities;

            const std::unordered_set<int> inCav(cav.begin(), cav.end());
            const CavityFill::Boundary B = CavityFill::boundaryOf(C_, cav, inCav);
            if (!B.manifold || B.loops.size() != 1) { reject("cavity is not a disk"); continue; }
            bool ok = true;
            int oldDefect = 0;
            for (int w : B.interior) {
                if (C_.protectedVertex[w] || C_.designatedVertex[w]) { reject("protected vertex inside"); ok = false; break; }
                if (C_.boundaryVertex[w]) { reject("domain-boundary vertex inside"); ok = false; break; }
                if (lockedV[w]) { reject("touches a committed cavity"); ok = false; break; }
                oldDefect += C_.defect(w);
            }
            if (!ok) continue;
            const std::vector<int> &loop = B.loops[0];
            for (int w : loop) {
                if (lockedV[w]) { ok = false; break; }
                oldDefect += C_.defect(w);
            }
            if (!ok) { reject("touches a committed cavity"); continue; }
            if (loop.size() < 4 || loop.size() % 2 != 0) { reject("odd or tiny boundary"); continue; }

            std::vector<Proposal> props;
            propose(loop, B.angle[0], B.outside[0], props);
            for (Proposal &p : props) {
                int cells = 0;
                const int K = static_cast<int>(p.corners.size());
                const int N = static_cast<int>(loop.size());
                if (p.kind == 0) {
                    const int a = ((p.corners[1] - p.corners[0]) % N + N) % N;
                    const int b = ((p.corners[2] - p.corners[1]) % N + N) % N;
                    cells = a * b;
                } else {
                    for (int i = 0; i < K; ++i) cells += p.sigma[i] * p.sigma[(i + K - 1) % K];
                }
                const double dE = opts_.lambdaS * (p.defect - oldDefect) +
                                  opts_.lambdaQ * (cells - static_cast<int>(cav.size()));
                if (dE < -1e-9) {
                    options.push_back({dE, p.geo, static_cast<int>(cavities.size()), std::move(p)});
                }
            }
            cavities.push_back(cav);
            loops.push_back(loop);
        }
        if (options.empty()) continue;
        std::stable_sort(options.begin(), options.end(), [](const Option &a, const Option &b) {
            if (std::fabs(a.dE - b.dE) > 1e-9) return a.dE < b.dE;
            return a.geo < b.geo;
        });

        int tries = 0;
        for (const Option &o : options) {
            if (tries++ >= opts_.realisationsPerWitness) break;
            Patch P;
            if (!realise(loops[o.k], o.p, P)) { reject("no kernel point for a star centre"); continue; }
            ++report_.realised;
            double oldMin = 1.0;
            for (int q : cavities[o.k]) oldMin = std::min(oldMin, C_.minScaledJacobian(q));
            const double floor = std::min(opts_.minScaledJacobian, oldMin);
            const CavityFill::Verdict V =
                CavityFill::settle(C_, cavities[o.k], P, floor, opts_.smoothingIterations);
            if (!V.valid) { reject(V.reason); continue; }
            ++report_.certified;

            SquareCarrier::Edit ed;
            ed.removeCells = cavities[o.k];
            ed.newVertices = P.newVertices;
            ed.newOrigin = P.newOrigin;
            ed.newSourceEdge = P.newSourceEdge;
            ed.cells = P.cells;
            ed.material = material;
            ed.cellOrigin = SquareCarrier::CellOrigin::Rewrite;
            edits.push_back(std::move(ed));
            for (int q : cavities[o.k]) {
                lockedQ[q] = 1;
                for (int w : C_.cells[q]) lockedV[w] = 1;
            }
            ++report_.committed;
            if (o.p.kind == 0) ++report_.gridFills; else ++report_.starFills;
            break;
        }
    }

    if (!edits.empty()) C_.apply(edits);
    report_.defectAfter = totalDefect();
    report_.irregularAfter = irregularCount();
    report_.cellsAfter = C_.numCells();
    return static_cast<int>(edits.size());
}
