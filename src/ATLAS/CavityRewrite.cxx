#include "ATLAS/CavityRewrite.hxx"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <random>
#include <set>
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
void CavityRewrite::sideRange(const Loop &L, int from, int to, int &n, int &lo, int &hi) const {
    const int N = static_cast<int>(L.ids.size());
    n = ((to - from) % N + N) % N;
    if (n == 0) n = N;
    int drop = 0, flex = 0;
    for (int k = 0; k < n; ++k) {
        const int u = (from + k) % N;
        if (k > 0 && L.droppable[u]) ++drop;
        if (L.onBoundary[u]) ++flex;
    }
    lo = n - drop;
    hi = flex > 0 ? n + std::max(2, static_cast<int>(std::ceil(opts_.maxInsertFactor * n))) : n;
}

void CavityRewrite::propose(const Loop &L, std::vector<Proposal> &out) {
    out.clear();
    const int N = static_cast<int>(L.ids.size());
    const double maxCorner = opts_.maxCornerAngle * M_PI / 180.0;
    std::vector<int> delta(N);
    std::vector<double> wdelta(N);
    std::vector<double> cornerGeo(N), sideGeo(N);
    int base = 0;
    double wbase = 0.0;
    double geoBase = 0.0;
    std::vector<char> okCorner(N);
    for (int u = 0; u < N; ++u) {
        const int t = C_.targetValence(L.ids[u]);
        const int f1 = std::abs(L.outside[u] + 1 - t);
        const int f2 = std::abs(L.outside[u] + 2 - t);
        const double w = L.weight.empty() ? 1.0 : L.weight[u];
        base += f2;
        wbase += w * f2;
        delta[u] = f1 - f2;
        wdelta[u] = w * (f1 - f2);
        cornerGeo[u] = (L.angle[u] - M_PI_2) * (L.angle[u] - M_PI_2);
        sideGeo[u] = (L.angle[u] - M_PI) * (L.angle[u] - M_PI);
        geoBase += sideGeo[u];
        okCorner[u] = L.angle[u] < maxCorner && !L.droppable[u];
    }
    // The candidate corners the flexible enumerations draw from: the ones
    // that lower the defect most, then the squarest.
    auto pool = [&](int size) {
        std::vector<int> P;
        for (int u = 0; u < N; ++u) if (okCorner[u]) P.push_back(u);
        std::stable_sort(P.begin(), P.end(), [&](int x, int y) {
            if (delta[x] != delta[y]) return delta[x] < delta[y];
            return cornerGeo[x] < cornerGeo[y];
        });
        if (static_cast<int>(P.size()) > size) P.resize(size);
        std::sort(P.begin(), P.end());
        return P;
    };

    if (opts_.grids && N >= 3) {
        std::set<std::array<int, 4>> seen;
        auto addGrid = [&](std::array<int, 4> c) {
            if (!seen.insert(c).second) return;
            int n[4], lo[4], hi[4];
            for (int k = 0; k < 4; ++k) sideRange(L, c[k], c[(k + 1) % 4], n[k], lo[k], hi[k]);
            int m[2];
            for (int k = 0; k < 2; ++k) {
                const int a = std::max(lo[k], lo[k + 2]), b = std::min(hi[k], hi[k + 2]);
                if (a > b) return;
                // Least total change, then the mean of the two, then fewer.
                const int p = std::min(n[k], n[k + 2]), q = std::max(n[k], n[k + 2]);
                const int x = std::max(a, p), y = std::min(b, q);
                if (x <= y) {
                    m[k] = std::max(x, std::min(y, (n[k] + n[k + 2]) / 2));
                } else {
                    m[k] = (b < p) ? b : a;
                }
            }
            Proposal pr;
            pr.kind = 0;
            pr.corners.assign(c.begin(), c.end());
            pr.counts = {m[0], m[1], m[0], m[1]};
            pr.cells = m[0] * m[1];
            pr.defect = base;
            pr.weighted = wbase;
            pr.geo = geoBase;
            for (int k = 0; k < 4; ++k) {
                pr.defect += delta[c[k]];
                pr.weighted += wdelta[c[k]];
                pr.geo += cornerGeo[c[k]] - sideGeo[c[k]];
            }
            out.push_back(std::move(pr));
        };
        // Every corner set with equal opposite sides as they stand...
        if (N >= 4 && N % 2 == 0) {
            const int half = N / 2;
            for (int i0 = 0; i0 < half; ++i0) {
                for (int a = 1; a < half; ++a) {
                    std::array<int, 4> c{i0, (i0 + a) % N, (i0 + half) % N, (i0 + half + a) % N};
                    if (!okCorner[c[0]] || !okCorner[c[1]] || !okCorner[c[2]] || !okCorner[c[3]]) continue;
                    std::sort(c.begin(), c.end());
                    addGrid(c);
                }
            }
        }
        // ... and every set from the pool, its sides matched on dS.
        const std::vector<int> P = pool(opts_.cornerPool);
        const int np = static_cast<int>(P.size());
        for (int a = 0; a < np; ++a)
            for (int b = a + 1; b < np; ++b)
                for (int c = b + 1; c < np; ++c)
                    for (int d = c + 1; d < np; ++d) addGrid({P[a], P[b], P[c], P[d]});
    }

    if (opts_.stars && N >= 3) {
        const std::vector<int> P = pool(9);
        const int np = static_cast<int>(P.size());
        for (int K : {3, 5}) {
            if (np < K) continue;
            std::vector<int> pick(K);
            std::function<void(int, int)> rec = [&](int start, int depth) {
                if (depth == K) {
                    std::vector<int> n(K), lo(K), hi(K), s, have;
                    for (int k = 0; k < K; ++k) sideRange(L, pick[k], pick[(k + 1) % K], n[k], lo[k], hi[k]);
                    if (!CavityFill::rangeStar(n, lo, hi, s, have, 3)) return;
                    Proposal pr;
                    pr.kind = 1;
                    pr.corners = pick;
                    pr.sigma = s;
                    pr.counts = have;
                    for (int k = 0; k < K; ++k) pr.cells += s[k] * s[(k + K - 1) % K];
                    pr.defect = base + 1;   // the centre has valence K, one off
                    pr.weighted = wbase + 1.0;
                    pr.geo = geoBase;
                    for (int c : pick) {
                        pr.defect += delta[c];
                        pr.weighted += wdelta[c];
                        pr.geo += cornerGeo[c] - sideGeo[c];
                    }
                    out.push_back(std::move(pr));
                    return;
                }
                for (int i = start; i < np; ++i) {
                    pick[depth] = P[i];
                    rec(i + 1, depth + 1);
                }
            };
            rec(0, 0);
        }
    }
    report_.proposals += static_cast<int>(out.size());
}

bool CavityRewrite::realise(const Loop &L, const Proposal &p, Patch &P) const {
    const int N = static_cast<int>(L.ids.size());
    const int K = static_cast<int>(p.corners.size());
    std::vector<std::vector<int>> arcs(K);
    for (int k = 0; k < K; ++k) {
        int i = p.corners[k];
        const int e = p.corners[(k + 1) % K];
        std::vector<int> ids{L.ids[i]};
        std::vector<char> may{0};
        do {
            i = (i + 1) % N;
            ids.push_back(L.ids[i]);
            may.push_back(L.droppable[i]);
        } while (i != e);
        may.back() = 0;
        const int n = static_cast<int>(ids.size()) - 1;
        if (p.counts[k] < n) ids = CavityFill::dropFromBoundary(C_, ids, may, n - p.counts[k]);
        else if (p.counts[k] > n) ids = CavityFill::insertOnBoundary(C_, P, ids, p.counts[k] - n);
        if (ids.empty()) return false;
        arcs[k] = std::move(ids);
    }
    if (p.kind == 0) {
        std::array<std::vector<int>, 4> sides{arcs[0], arcs[1], arcs[2], arcs[3]};
        P.kind = "grid";
        P.blocks = 1;
        return CavityFill::grid(C_, P, sides, Origin::Rewrite);
    }
    std::vector<Point> poly;
    for (const auto &a : arcs) for (size_t j = 0; j + 1 < a.size(); ++j) poly.push_back(P.at(C_, a[j]));
    const Point g = CavityFill::areaCentroid(poly);
    Point c = CavityFill::kernelCenter(poly, g);
    // No kernel point: the centroid, and smoothing and the certificate decide.
    if (!(CavityFill::kernelDepth(poly, c) > 0.0)) c = g;
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
        std::vector<std::vector<int>> cavities;
        std::vector<Loop> loops;
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
            Loop L;
            L.ids = B.loops[0];
            L.angle = B.angle[0];
            L.outside = B.outside[0];
            const int N = static_cast<int>(L.ids.size());
            L.droppable.assign(N, 0);
            L.onBoundary.assign(N, 0);
            for (int u = 0; u < N; ++u) {
                const int w = L.ids[u];
                if (lockedV[w]) { ok = false; break; }
                oldDefect += C_.defect(w);
                L.droppable[u] = CavityFill::droppable(C_, w, inCav) ? 1 : 0;
                L.onBoundary[u] = C_.boundaryEdge[B.loopEdges[0][u]] ? 1 : 0;
            }
            if (!ok) { reject("touches a committed cavity"); continue; }
            if (N < 3) { reject("tiny boundary"); continue; }

            std::vector<Proposal> props;
            propose(L, props);
            for (Proposal &p : props) {
                const double dE = opts_.lambdaS * (p.defect - oldDefect) +
                                  opts_.lambdaQ * (p.cells - static_cast<int>(cav.size()));
                if (dE < -1e-9) {
                    options.push_back({dE, p.geo, static_cast<int>(cavities.size()), std::move(p)});
                }
            }
            cavities.push_back(cav);
            loops.push_back(std::move(L));
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
            if (!realise(loops[o.k], o.p, P)) { reject("the boundary could not be resubdivided"); continue; }
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

// ---------------------------------------------------------------------------
// The annealed search (Sec. 8.3).
// ---------------------------------------------------------------------------
bool CavityRewrite::growCavity(const std::vector<int> &seeds, int rings, const std::unordered_set<int> &exclude,
                               std::vector<int> &cav, Loop &L, int &material, std::unordered_set<int> &inCav) {
    cav.clear();
    inCav.clear();
    material = -1;
    for (int v : seeds) {
        for (int r = C_.ringPtr[v]; r < C_.ringPtr[v + 1] && material < 0; ++r) {
            const int q = C_.ringCell[r];
            if (exclude.count(q)) continue;
            if (C_.cellOrigin[q] != SquareCarrier::CellOrigin::Template) material = C_.cellMaterial[q];
        }
    }
    if (material < 0) return false;
    std::vector<int> frontier = seeds;
    std::unordered_set<int> seenV(seeds.begin(), seeds.end());
    for (int k = 1; k <= rings; ++k) {
        std::vector<int> next;
        for (int u : frontier) {
            for (int r = C_.ringPtr[u]; r < C_.ringPtr[u + 1]; ++r) {
                const int q = C_.ringCell[r];
                if (inCav.count(q) || exclude.count(q) || C_.cellMaterial[q] != material) continue;
                if (C_.cellOrigin[q] == SquareCarrier::CellOrigin::Template) continue;
                inCav.insert(q);
                cav.push_back(q);
                for (int w : C_.cells[q]) if (seenV.insert(w).second) next.push_back(w);
            }
        }
        frontier.swap(next);
    }
    if (cav.empty() || static_cast<int>(cav.size()) > opts_.maxCavityCells) return false;
    const CavityFill::Boundary B = CavityFill::boundaryOf(C_, cav, inCav);
    if (!B.manifold || B.loops.size() != 1) return false;
    for (int w : B.interior) {
        if (C_.protectedVertex[w] || C_.designatedVertex[w]) return false;
    }
    L.ids = B.loops[0];
    L.angle = B.angle[0];
    L.outside = B.outside[0];
    const int N = static_cast<int>(L.ids.size());
    if (N < 3) return false;
    L.droppable.assign(N, 0);
    L.onBoundary.assign(N, 0);
    for (int u = 0; u < N; ++u) {
        L.droppable[u] = CavityFill::droppable(C_, L.ids[u], inCav) ? 1 : 0;
        L.onBoundary[u] = C_.boundaryEdge[B.loopEdges[0][u]] ? 1 : 0;
    }
    return true;
}

double CavityRewrite::energy(const SquareCarrier &C, const AnnealOptions &ao, int *blocks) {
    const int nb = RectangleCertifier::basePatchCount(C);
    double d = 0.0;
    for (int v = 0; v < C.numVertices(); ++v) {
        if (C.valence[v] == 0) continue;
        d += (C.boundaryVertex[v] ? ao.wBoundaryDefect : ao.wDefect) * C.defect(v);
    }
    double shape = 0.0;
    if (ao.wShape > 0.0) {
        for (int q = 0; q < C.numCells(); ++q) {
            const double sj = C.minScaledJacobian(q);
            if (sj < ao.shapeFloor) shape += (ao.shapeFloor - sj) / ao.shapeFloor;
        }
    }
    double align = 0.0;
    if (ao.field && ao.wAlign > 0.0) {
        for (int q = 0; q < C.numCells(); ++q) {
            const auto &c = C.cells[q];
            const Point &p0 = C.vertices[c[0]], &p1 = C.vertices[c[1]], &p2 = C.vertices[c[2]], &p3 = C.vertices[c[3]];
            // The cell's own cross: its two mean edge directions, the second
            // turned back a quarter, averaged as 4-fold representation vectors.
            const Point du = (p1 - p0) + (p2 - p3), dv = (p3 - p0) + (p2 - p1);
            const double a = 4.0 * std::atan2(du[1], du[0]);
            const double b = 4.0 * (std::atan2(dv[1], dv[0]) - M_PI_2);
            const double cellTh = std::atan2(std::sin(a) + std::sin(b), std::cos(a) + std::cos(b)) / 4.0;
            double w = 0.0;
            const double fieldTh = ao.field(C.cellCentroid(q), &w);
            align += w * 0.5 * (1.0 - std::cos(4.0 * (cellTh - fieldTh)));
        }
    }
    if (blocks) *blocks = nb;
    return ao.wBlocks * nb + d + ao.wCells * C.numCells() + ao.wShape * shape + ao.wAlign * align;
}

const CavityRewrite::AnnealReport &CavityRewrite::anneal(const AnnealOptions &ao) {
    const auto t0 = std::chrono::steady_clock::now();
    AnnealReport &R = annealReport_;
    R = AnnealReport();
    std::mt19937 rng(ao.seed);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    int blocks = 0;
    double E = energy(C_, ao, &blocks);
    R.energyBefore = E;
    R.blocksBefore = blocks;
    SquareCarrier best = C_;
    double Ebest = E;
    int blocksBest = blocks;

    std::vector<int> cav;
    std::unordered_set<int> inCav;
    std::vector<Proposal> props;
    checkpoints_.clear();
    double Echeck = 1e300;
    for (int move = 0; move < ao.maxMoves; ++move) {
        // The best state at each quarter of the schedule, kept for a caller
        // that finds the final one cannot be realised (see ATLAS).
        if (move > 0 && move % std::max(1, ao.maxMoves / 4) == 0 && Ebest < Echeck - 1e-9) {
            checkpoints_.push_back(best);
            Echeck = Ebest;
        }
        // The schedule runs on the move count, so a run is reproducible from
        // its seed; the clock is only a safety cap.
        const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (elapsed > ao.seconds) { R.timedOut = true; break; }
        const double progress = static_cast<double>(move) / std::max(1, ao.maxMoves);
        const double T = ao.T0 * std::pow(ao.T1 / ao.T0, progress);
        ++R.moves;

        // A witness: a vertex of wrong valence mostly, any vertex sometimes --
        // aligning two singularities can need a move next to neither.
        int v = -1;
        {
            std::vector<int> W;
            for (int x = 0; x < C_.numVertices(); ++x) if (C_.valence[x] > 0 && C_.defect(x) > 0) W.push_back(x);
            if (!W.empty() && U(rng) < 0.85) {
                v = W[std::uniform_int_distribution<int>(0, static_cast<int>(W.size()) - 1)(rng)];
            } else {
                v = std::uniform_int_distribution<int>(0, C_.numVertices() - 1)(rng);
                if (C_.valence[v] == 0) continue;
            }
        }
        const double pr = U(rng);
        int rings = pr < 0.35 ? 1 : (pr < 0.75 ? 2 : std::min(3, ao.maxRings));
        // The cavity's shape: rings round the witness, round one of its edges,
        // or -- on dS or at a protected vertex -- round part of its fan only.
        // A ring always takes every cell at the witness, so a fill can leave
        // it at most two (the fill's corner or side), which is never enough
        // for a reflex corner; keeping some of its cells outside is what
        // lets a fill set its valence to anything.
        std::vector<int> seeds{v};
        std::unordered_set<int> exclude;
        const double shapePick = U(rng);
        const int q = C_.valence[v];
        if ((C_.boundaryVertex[v] || C_.protectedVertex[v]) && q >= 2 && shapePick < 0.5) {
            const int keep = std::uniform_int_distribution<int>(1, q - 1)(rng);
            const int off = std::uniform_int_distribution<int>(0, q - 1)(rng);
            std::vector<int> kept;
            for (int k = 0; k < q; ++k) {
                const int cell = C_.ringCell[C_.ringPtr[v] + (off + k) % q];
                if (k < keep) kept.push_back(cell);
                else exclude.insert(cell);
            }
            seeds.clear();
            for (int cell : kept) for (int w : C_.cells[cell]) if (w != v) seeds.push_back(w);
            std::sort(seeds.begin(), seeds.end());
            seeds.erase(std::unique(seeds.begin(), seeds.end()), seeds.end());
            rings = std::max(1, rings - 1);
        } else if (shapePick < 0.75) {
            std::vector<int> es;
            C_.vertexEdges(v, es);
            if (!es.empty()) {
                const int e = es[std::uniform_int_distribution<int>(0, static_cast<int>(es.size()) - 1)(rng)];
                seeds.push_back(C_.otherEnd(e, v));
            }
        }
        Loop L;
        int material = -1;
        if (!growCavity(seeds, rings, exclude, cav, L, material, inCav)) continue;
        ++R.cavities;
        const int N = static_cast<int>(L.ids.size());
        L.weight.assign(N, 1.0);
        double oldWeighted = 0.0;
        for (int u = 0; u < N; ++u) {
            L.weight[u] = C_.boundaryVertex[L.ids[u]] ? ao.wBoundaryDefect / ao.wDefect : 1.0;
            oldWeighted += L.weight[u] * C_.defect(L.ids[u]);
        }
        {
            std::unordered_set<int> onLoop(L.ids.begin(), L.ids.end());
            for (int q : cav) {
                for (int w : C_.cells[q]) {
                    if (onLoop.insert(w).second) oldWeighted += C_.defect(w);
                }
            }
        }
        propose(L, props);
        if (props.empty()) continue;
        R.proposals += static_cast<int>(props.size());
        // Local screen in the energy's own units, then a softmax pick among
        // the best few so that equal-defect moves -- the ones that shift a
        // singularity -- get tried as well as the improving ones.
        struct Cand { double d; int k; };
        std::vector<Cand> cand;
        for (int k = 0; k < static_cast<int>(props.size()); ++k) {
            const double d = ao.wDefect * (props[k].weighted - oldWeighted) +
                             ao.wCells * (props[k].cells - static_cast<int>(cav.size()));
            cand.push_back({d, k});
        }
        std::sort(cand.begin(), cand.end(), [](const Cand &a, const Cand &b) { return a.d < b.d; });
        if (cand.size() > 12) cand.resize(12);
        const double Tloc = std::max(0.5, T);
        double Z = 0.0;
        for (const Cand &c : cand) Z += std::exp(-(c.d - cand[0].d) / Tloc);

        const double floor = ao.minScaledJacobian;
        for (int attempt = 0; attempt < ao.realisations; ++attempt) {
            // Draw a candidate.
            double x = U(rng) * Z;
            int pick = 0;
            for (int j = 0; j < static_cast<int>(cand.size()); ++j) {
                x -= std::exp(-(cand[j].d - cand[0].d) / Tloc);
                if (x <= 0.0) { pick = j; break; }
                pick = j;
            }
            const Proposal &p = props[cand[pick].k];
            Patch P;
            if (!realise(L, p, P)) continue;
            const CavityFill::Verdict V = CavityFill::settle(C_, cav, P, floor, opts_.smoothingIterations);
            if (!V.valid) continue;
            ++R.certified;
            SquareCarrier trial = C_;
            SquareCarrier::Edit ed;
            ed.removeCells = cav;
            ed.newVertices = P.newVertices;
            ed.newOrigin = P.newOrigin;
            ed.newSourceEdge = P.newSourceEdge;
            ed.cells = P.cells;
            ed.material = material;
            ed.cellOrigin = SquareCarrier::CellOrigin::Rewrite;
            trial.apply({ed});
            int nb = 0;
            const double E2 = energy(trial, ao, &nb);
            const double dE = E2 - E;
            if (dE <= 0.0 || U(rng) < std::exp(-dE / T)) {
                if (dE > 0.0) ++R.uphill;
                ++R.accepted;
                C_ = std::move(trial);
                E = E2;
                if (E < Ebest - 1e-9) {
                    best = C_;
                    Ebest = E;
                    blocksBest = nb;
                    ++R.improvements;
                }
            }
            break;
        }
    }
    C_ = std::move(best);
    R.energyAfter = Ebest;
    R.blocksAfter = blocksBest;
    R.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    report_.defectAfter = totalDefect();
    report_.irregularAfter = irregularCount();
    report_.cellsAfter = C_.numCells();
    return R;
}
