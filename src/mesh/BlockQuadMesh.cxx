#include "mesh/BlockQuadMesh.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

namespace {

// ---------------------------------------------------------------------------
// One block's four sides, read out once.
//
// Every question asked of a block's geometry -- how long its isoparametric
// lines are, where the Coons blend puts a point -- wants the four side
// polylines and their arc lengths, and asking BlockDecomposition for them each
// time rebuilds and reverses four vectors per call. On a layout of a thousand
// blocks sampled at a few hundred points each that is the whole cost of the
// stage, so the sides are read once per block and the sampling written against
// them.
// ---------------------------------------------------------------------------
struct BlockGeom {
    std::vector<Point> side[4];
    std::vector<double> len[4];   // cumulative arc length, len[s].back() the total

    BlockGeom(const BlockDecomposition &D, int b) {
        for (int s = 0; s < 4; ++s) {
            side[s] = D.sidePolyline(b, s);
            len[s].assign(side[s].size(), 0.0);
            for (size_t i = 1; i < side[s].size(); ++i)
                len[s][i] = len[s][i - 1] + normP(side[s][i] - side[s][i - 1]);
        }
    }

    bool complete() const {
        for (int s = 0; s < 4; ++s) if (side[s].size() < 2) return false;
        return true;
    }

    // The point at arc-length fraction `t` along side `s`, corner to corner.
    Point at(int s, double t) const {
        const std::vector<Point> &p = side[s];
        if (p.empty()) return Point{0.0, 0.0};
        if (p.size() == 1) return p[0];
        const double total = len[s].back();
        if (!(total > 0.0)) return p[0];
        const double target = std::min(std::max(t, 0.0), 1.0) * total;
        const size_t hi = static_cast<size_t>(
            std::lower_bound(len[s].begin(), len[s].end(), target) - len[s].begin());
        if (hi == 0) return p.front();
        if (hi >= p.size()) return p.back();
        const double seg = len[s][hi] - len[s][hi - 1];
        const double a = (seg > 0.0) ? (target - len[s][hi - 1]) / seg : 0.0;
        return p[hi - 1] * (1.0 - a) + p[hi] * a;
    }

    // BlockDecomposition::coonsPoint, on the sides already in hand.
    Point coons(double u, double v) const {
        const Point b = at(0, u), r = at(1, v), t = at(2, 1.0 - u), l = at(3, 1.0 - v);
        const Point c00 = side[0].front(), c10 = side[1].front();
        const Point c11 = side[2].front(), c01 = side[3].front();
        return b * (1.0 - v) + t * v + l * (1.0 - u) + r * u -
               (c00 * ((1.0 - u) * (1.0 - v)) + c10 * (u * (1.0 - v)) +
                c11 * (u * v) + c01 * ((1.0 - u) * v));
    }
};

// The mean length of the block's isoparametric lines in the u direction (the
// one sides 0 and 2 are cut into) when `uDirection`, in the v direction
// otherwise.
//
// Why this rather than the side length: the count on a side does not only cut
// the side, it cuts every row of elements across the block all the way to the
// side facing it. On a block that is nearly a rectangle those are the same
// number. On a lens-shaped one whose two ends taper and whose middle is wide
// they are not close, and judging by the ends puts one element across a middle
// many times the target -- the case Stage 10's header records on geom024.
double isoSpan(const BlockGeom &g, bool uDirection) {
    constexpr int kAcross = 5;    // isolines sampled
    constexpr int kAlong = 16;    // steps along each
    double sum = 0.0;
    for (int k = 0; k < kAcross; ++k) {
        const double w = static_cast<double>(k) / (kAcross - 1);
        double L = 0.0;
        Point prev = uDirection ? g.coons(0.0, w) : g.coons(w, 0.0);
        for (int i = 1; i <= kAlong; ++i) {
            const double a = static_cast<double>(i) / kAlong;
            const Point cur = uDirection ? g.coons(a, w) : g.coons(w, a);
            L += normP(cur - prev);
            prev = cur;
        }
        sum += L;
    }
    return sum / kAcross;
}

double quadArea(const Point c[4]) {
    double a = 0.0;
    for (int k = 0; k < 4; ++k) a += cross2(c[k], c[(k + 1) % 4]);
    return 0.5 * a;
}

// The worst of the four corner Jacobians, normalized by the two incident edge
// lengths -- the standard scaled Jacobian, and the quantity that tells a
// folded element from a merely stretched one.
double scaledJacobian(const Point c[4]) {
    double worst = std::numeric_limits<double>::infinity();
    for (int k = 0; k < 4; ++k) {
        const Point u = c[(k + 1) % 4] - c[k];
        const Point v = c[(k + 3) % 4] - c[k];
        const double nu = normP(u), nv = normP(v);
        if (!(nu > 0.0) || !(nv > 0.0)) return 0.0;
        worst = std::min(worst, cross2(u, v) / (nu * nv));
    }
    return std::isfinite(worst) ? worst : 0.0;
}

double segmentDistance(const Point &p, const Point &a, const Point &b) {
    const Point d = b - a;
    const double dd = dotP(d, d);
    if (!(dd > 0.0)) return normP(p - a);
    double s = dotP(p - a, d) / dd;
    s = std::max(0.0, std::min(1.0, s));
    return normP(p - (a + d * s));
}

// Union-find over the macro edges, for the chords.
struct DisjointSet {
    std::vector<int> parent;
    explicit DisjointSet(size_t n) : parent(n) { std::iota(parent.begin(), parent.end(), 0); }
    int find(int x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    }
    void unite(int a, int b) {
        a = find(a); b = find(b);
        if (a != b) parent[b] = a;
    }
};

} // namespace

// ---------------------------------------------------------------------------
BlockQuadMesh::BlockQuadMesh(const BlockDecomposition &D, const Options &opts) : opts_(opts) {
    report_.target = opts_.targetEdgeLength;

    Point lo{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    Point hi{-lo[0], -lo[1]};
    for (const auto &e : D.edges) {
        for (const Point &p : e.points) {
            lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
            hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
        }
    }
    report_.modelExtent = (std::isfinite(lo[0]) && hi[0] > lo[0]) ? normP(hi - lo) : 1.0;
    if (!(report_.modelExtent > 0.0)) report_.modelExtent = 1.0;

    assignIntervals(D);
    meshEdges(D);
    for (size_t e = 0; e < D.edges.size(); ++e) {
        const BlockDecomposition::MacroEdge &me = D.edges[e];
        if (me.boundary || me.blockB >= 0 || me.blockA < 0) continue;
        for (const int v : edgeNodes_[e]) openVerts_.push_back(v);
    }
    std::sort(openVerts_.begin(), openVerts_.end());
    openVerts_.erase(std::unique(openVerts_.begin(), openVerts_.end()), openVerts_.end());
    meshBlocks(D);
    smooth();
    check(D);
}

// ---------------------------------------------------------------------------
// assignIntervals()
//
// "Opposite sides of a block carry the same count" is an equivalence relation
// on the macro edges, and a class of it is a chord. The classes come from a
// union-find over the edges -- side 0 with side 2 and side 1 with side 3 of
// every block -- rather than from a walk, because that way a chord closing on
// itself and a chord entering one block twice in perpendicular directions,
// both of which real layouts have, come out as the ordinary case.
//
// The count itself minimises F(N) = sum_k (log(S_k / (N h)))^2 over the blocks
// the chord runs through, which is the same objective and the same rounding
// Stage 10 and BlockMesh use: a ratio, so a row twice the target is as wrong
// as one half of it, and whichever of floor(N*) and ceil(N*) has the lower F
// rather than round(N*), which minimises a different objective and is a
// different integer exactly where N* is small.
// ---------------------------------------------------------------------------
void BlockQuadMesh::assignIntervals(const BlockDecomposition &D) {
    const int nE = static_cast<int>(D.edges.size());
    intervals_.assign(nE, 0);
    chordOf_.assign(nE, -1);
    if (nE == 0) return;

    DisjointSet ds(static_cast<size_t>(nE));
    for (const auto &b : D.blocks) {
        for (int pair = 0; pair < 2; ++pair) {
            const int a = b.edges[pair], c = b.edges[pair + 2];
            if (a >= 0 && a < nE && c >= 0 && c < nE) ds.unite(a, c);
        }
    }

    std::unordered_map<int, int> chordIndex;
    for (int e = 0; e < nE; ++e) {
        const int root = ds.find(e);
        auto it = chordIndex.find(root);
        if (it == chordIndex.end()) {
            it = chordIndex.emplace(root, static_cast<int>(chords_.size())).first;
            chords_.emplace_back();
        }
        chordOf_[e] = it->second;
        chords_[it->second].edges.push_back(e);
    }

    // The spans: one per (block, direction) the chord passes through. A block
    // contributes its u span to the chord of sides 0 and 2 and its v span to
    // the chord of sides 1 and 3, so a chord that runs through one block in
    // both directions is counted for both, which is what it has to span.
    std::vector<std::vector<double>> spans(chords_.size());
    for (size_t b = 0; b < D.blocks.size(); ++b) {
        const BlockGeom g(D, static_cast<int>(b));
        if (!g.complete()) continue;
        for (int pair = 0; pair < 2; ++pair) {
            const int e = D.blocks[b].edges[pair];
            if (e < 0 || e >= nE) continue;
            spans[chordOf_[e]].push_back(isoSpan(g, pair == 0));
        }
    }

    const double h = (opts_.targetEdgeLength > 0.0) ? opts_.targetEdgeLength : 1.0;
    long long intervalSum = 0;
    report_.minIntervals = std::numeric_limits<int>::max();
    for (size_t c = 0; c < chords_.size(); ++c) {
        Chord &ch = chords_[c];
        ch.minLength = std::numeric_limits<double>::infinity();
        for (int e : ch.edges) {
            double L = 0.0;
            for (size_t i = 1; i < D.edges[e].points.size(); ++i)
                L += normP(D.edges[e].points[i] - D.edges[e].points[i - 1]);
            ch.minLength = std::min(ch.minLength, L);
            ch.maxLength = std::max(ch.maxLength, L);
        }
        if (!std::isfinite(ch.minLength)) ch.minLength = 0.0;

        // A chord no block claimed -- every one of its edges unmatched -- has
        // no isoline to be judged by, so it falls back to its own length,
        // which is what a side has to span when nothing spans across it.
        std::vector<double> &S = spans[c];
        if (S.empty()) {
            for (int e : ch.edges) {
                double L = 0.0;
                for (size_t i = 1; i < D.edges[e].points.size(); ++i)
                    L += normP(D.edges[e].points[i] - D.edges[e].points[i - 1]);
                if (L > 0.0) S.push_back(L);
            }
        }
        ch.minSpan = S.empty() ? 0.0 : *std::min_element(S.begin(), S.end());
        ch.maxSpan = S.empty() ? 0.0 : *std::max_element(S.begin(), S.end());

        double logSum = 0.0;
        int n = 0;
        for (double s : S) if (s > 0.0) { logSum += std::log(s); ++n; }
        const double ideal = (n > 0) ? std::exp(logSum / n) / h : 1.0;
        ch.idealIntervals = ideal;

        auto F = [&](int N) {
            if (N < 1) return std::numeric_limits<double>::infinity();
            double f = 0.0;
            for (double s : S) {
                if (!(s > 0.0)) continue;
                const double r = std::log(s / (N * h));
                f += r * r;
            }
            return f;
        };
        const int lo = std::max(1, static_cast<int>(std::floor(ideal)));
        const int hiN = std::max(1, static_cast<int>(std::ceil(ideal)));
        int N = (F(lo) <= F(hiN)) ? lo : hiN;

        const int floorN = std::max(1, opts_.minIntervals);
        if (N < floorN) { N = floorN; ch.clamped = true; }
        if (opts_.maxIntervals > 0 && N > opts_.maxIntervals) {
            N = opts_.maxIntervals;
            ch.clamped = true;
        }
        ch.intervals = N;
        if (ch.clamped) ++report_.clampedChords;
        for (int e : ch.edges) intervals_[e] = N;

        report_.minIntervals = std::min(report_.minIntervals, N);
        report_.maxIntervals = std::max(report_.maxIntervals, N);
        intervalSum += static_cast<long long>(N) * static_cast<long long>(ch.edges.size());
    }
    if (report_.minIntervals == std::numeric_limits<int>::max()) report_.minIntervals = 0;
    report_.chords = static_cast<int>(chords_.size());
    report_.edgesAssigned = nE;
    report_.meanIntervals = nE > 0 ? static_cast<double>(intervalSum) / nE : 0.0;
}

// ---------------------------------------------------------------------------
// meshEdges()
//
// Every macro edge meshed once, and the two blocks that share it handed the
// same vertex indices. Conformity is then a property of the construction
// rather than a tolerance to be met: there is no second copy of a shared node
// anywhere for a later check to have to find.
//
// The macrovertices go down first and are shared by every edge that ends
// there, so a corner where four blocks meet is one vertex and not four at the
// same point.
// ---------------------------------------------------------------------------
void BlockQuadMesh::meshEdges(const BlockDecomposition &D) {
    std::vector<int> nodeOfVertex(D.vertices.size(), -1);
    for (size_t v = 0; v < D.vertices.size(); ++v) {
        nodeOfVertex[v] = static_cast<int>(verts_.size());
        verts_.push_back(D.vertices[v].p);
    }

    edgeNodes_.assign(D.edges.size(), {});
    for (size_t e = 0; e < D.edges.size(); ++e) {
        const BlockDecomposition::MacroEdge &me = D.edges[e];
        const int N = std::max(1, intervals_[e]);
        if (me.points.size() < 2) continue;

        std::vector<double> s(me.points.size(), 0.0);
        for (size_t i = 1; i < me.points.size(); ++i)
            s[i] = s[i - 1] + normP(me.points[i] - me.points[i - 1]);
        const double total = s.back();

        std::vector<int> &nodes = edgeNodes_[e];
        nodes.reserve(N + 1);
        // The ends are the macrovertices themselves; only the interior is new.
        nodes.push_back(me.from >= 0 ? nodeOfVertex[me.from]
                                     : static_cast<int>(verts_.size()));
        if (me.from < 0) verts_.push_back(me.points.front());
        for (int k = 1; k < N; ++k) {
            const double target = total * k / N;
            const size_t hi = static_cast<size_t>(
                std::lower_bound(s.begin(), s.end(), target) - s.begin());
            Point p = me.points.back();
            if (hi == 0) {
                p = me.points.front();
            } else if (hi < me.points.size()) {
                const double seg = s[hi] - s[hi - 1];
                const double a = (seg > 0.0) ? (target - s[hi - 1]) / seg : 0.0;
                p = me.points[hi - 1] * (1.0 - a) + me.points[hi] * a;
            }
            nodes.push_back(static_cast<int>(verts_.size()));
            verts_.push_back(p);
        }
        nodes.push_back(me.to >= 0 ? nodeOfVertex[me.to] : static_cast<int>(verts_.size()));
        if (me.to < 0) verts_.push_back(me.points.back());
    }
}

// ---------------------------------------------------------------------------
// meshBlocks()
//
// Each block an ns x nt grid: the four sides' nodes as they were placed, and
// the interior by the discrete Coons blend of them. The blend reproduces the
// four sides exactly on the four edges of the parameter square, so a block's
// grid meets its neighbours' whatever it does inside.
// ---------------------------------------------------------------------------
void BlockQuadMesh::meshBlocks(const BlockDecomposition &D) {
    grids_.reserve(D.blocks.size());
    for (size_t b = 0; b < D.blocks.size(); ++b) {
        const BlockDecomposition::Block &blk = D.blocks[b];

        // The four sides' node lists, each read in the direction this block
        // walks that side.
        std::vector<int> sideNodes[4];
        bool ok = true;
        for (int s = 0; s < 4; ++s) {
            const int e = blk.edges[s];
            if (e < 0 || e >= static_cast<int>(edgeNodes_.size()) || edgeNodes_[e].size() < 2) {
                ok = false;
                break;
            }
            sideNodes[s] = edgeNodes_[e];
            if (blk.flip[s]) std::reverse(sideNodes[s].begin(), sideNodes[s].end());
        }
        if (!ok) { ++report_.unmeshedBlocks; continue; }

        const int ns = static_cast<int>(sideNodes[0].size()) - 1;
        const int nt = static_cast<int>(sideNodes[1].size()) - 1;
        // The chording guarantees this, so a mismatch is a decomposition whose
        // opposite sides were never identified -- a block left out rather than
        // a grid built on a count that does not fit.
        if (ns < 1 || nt < 1 ||
            static_cast<int>(sideNodes[2].size()) - 1 != ns ||
            static_cast<int>(sideNodes[3].size()) - 1 != nt) {
            ++report_.unmeshedBlocks;
            continue;
        }

        Block g;
        g.block = static_cast<int>(b);
        g.ns = ns;
        g.nt = nt;
        g.material = blk.material;
        g.vert.assign(static_cast<size_t>(ns + 1) * (nt + 1), -1);
        const int w = ns + 1;
        auto id = [&](int i, int j) -> int & { return g.vert[static_cast<size_t>(j) * w + i]; };

        // Side 0 runs corner 0 -> 1 along the bottom; side 1 corner 1 -> 2 up
        // the right; side 2 corner 2 -> 3, which is the top backwards; side 3
        // corner 3 -> 0, the left backwards.
        for (int i = 0; i <= ns; ++i) {
            id(i, 0) = sideNodes[0][i];
            id(i, nt) = sideNodes[2][ns - i];
        }
        for (int j = 0; j <= nt; ++j) {
            id(ns, j) = sideNodes[1][j];
            id(0, j) = sideNodes[3][nt - j];
        }

        const Point c00 = verts_[id(0, 0)], c10 = verts_[id(ns, 0)];
        const Point c11 = verts_[id(ns, nt)], c01 = verts_[id(0, nt)];
        for (int j = 1; j < nt; ++j) {
            const double v = static_cast<double>(j) / nt;
            for (int i = 1; i < ns; ++i) {
                const double u = static_cast<double>(i) / ns;
                const Point p = verts_[id(i, 0)] * (1.0 - v) + verts_[id(i, nt)] * v +
                                verts_[id(0, j)] * (1.0 - u) + verts_[id(ns, j)] * u -
                                (c00 * ((1.0 - u) * (1.0 - v)) + c10 * (u * (1.0 - v)) +
                                 c11 * (u * v) + c01 * ((1.0 - u) * v));
                id(i, j) = static_cast<int>(verts_.size());
                verts_.push_back(p);
            }
        }

        for (int j = 0; j < nt; ++j) {
            for (int i = 0; i < ns; ++i) {
                cells_.push_back({id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)});
                cellMaterial_.push_back(blk.material);
            }
        }
        grids_.push_back(std::move(g));
    }
    report_.blocks = static_cast<int>(grids_.size());
}

// ---------------------------------------------------------------------------
// smooth()
//
// Winslow's elliptic system on the interior nodes of the blocks that need it,
// the boundary held -- Stage 10's pass, for Stage 10's reason. A transfinite
// grid inherits whatever its four sides do, and where one side is far more
// curved than the side facing it the blend crowds its rows together and can
// push them through each other. Solving
//
//     alpha x_ss - 2 beta x_st + gamma x_tt = 0,
//     alpha = |x_t|^2,  beta = x_s . x_t,  gamma = |x_s|^2,
//
// by Gauss-Seidel gives the inverse of a harmonic map onto the square, which
// has no interior extremum to fold at. Holding the boundary is what keeps the
// mesh conforming while it does so.
// ---------------------------------------------------------------------------
void BlockQuadMesh::smooth() {
    double before = std::numeric_limits<double>::infinity();
    for (const std::array<int, 4> &q : cells_) {
        const Point c[4] = {verts_[q[0]], verts_[q[1]], verts_[q[2]], verts_[q[3]]};
        const double sj = scaledJacobian(c);
        before = std::min(before, sj);
        if (!(sj > 0.0) || !(quadArea(c) > 0.0)) ++report_.invertedBefore;
    }
    report_.minScaledJacobianBefore = std::isfinite(before) ? before : 0.0;
    if (opts_.smoothingPasses <= 0 || cells_.empty()) return;

    const double tol = opts_.smoothingTolerance * opts_.targetEdgeLength;
    for (const Block &b : grids_) {
        if (b.ns < 2 || b.nt < 2) continue;
        const int w = b.ns + 1;
        auto id = [&](int i, int j) { return b.vert[static_cast<size_t>(j) * w + i]; };
        double worst = std::numeric_limits<double>::infinity();
        for (int j = 0; j < b.nt; ++j) {
            for (int i = 0; i < b.ns; ++i) {
                const Point c[4] = {verts_[id(i, j)], verts_[id(i + 1, j)],
                                    verts_[id(i + 1, j + 1)], verts_[id(i, j + 1)]};
                worst = std::min(worst, scaledJacobian(c));
            }
        }
        if (!(worst < opts_.smoothingThreshold)) continue;
        ++report_.smoothedBlocks;

        int sweep = 0;
        for (; sweep < opts_.smoothingPasses; ++sweep) {
            double moved = 0.0;
            for (int j = 1; j < b.nt; ++j) {
                for (int i = 1; i < b.ns; ++i) {
                    const Point &xe = verts_[id(i + 1, j)];
                    const Point &xw = verts_[id(i - 1, j)];
                    const Point &xn = verts_[id(i, j + 1)];
                    const Point &xs = verts_[id(i, j - 1)];
                    const Point ds = (xe - xw) * 0.5;
                    const Point dt = (xn - xs) * 0.5;
                    const double alpha = dotP(dt, dt);
                    const double beta = dotP(ds, dt);
                    const double gamma = dotP(ds, ds);
                    const double denom = 2.0 * (alpha + gamma);
                    if (!(denom > 0.0)) continue;
                    const Point cross = (verts_[id(i + 1, j + 1)] - verts_[id(i - 1, j + 1)] -
                                         verts_[id(i + 1, j - 1)] + verts_[id(i - 1, j - 1)]) *
                                        0.25;
                    const Point next =
                        ((xe + xw) * alpha + (xn + xs) * gamma - cross * (2.0 * beta)) / denom;
                    Point &here = verts_[id(i, j)];
                    moved = std::max(moved, normP(next - here));
                    here = next;
                }
            }
            if (moved <= tol) { ++sweep; break; }
        }
        report_.smoothingSweeps = std::max(report_.smoothingSweeps, sweep);
    }
}

// ---------------------------------------------------------------------------
// check()
//
// The corner Jacobians, the incidence and the boundary error, read back off
// the result rather than assumed. See the class comment for why each is a
// thing sampling can break.
// ---------------------------------------------------------------------------
void BlockQuadMesh::check(const BlockDecomposition &D) {
    report_.vertices = static_cast<int>(verts_.size());
    report_.quads = static_cast<int>(cells_.size());
    if (cells_.empty()) {
        report_.messages.push_back("no block of the decomposition could be meshed");
        return;
    }

    const long long n = static_cast<long long>(verts_.size());
    std::unordered_map<long long, std::array<int, 3>> edgeUse;  // count, first material, mixed
    edgeUse.reserve(cells_.size() * 4);
    double lenSum = 0.0, logSum = 0.0;
    int lenCount = 0;
    report_.minEdge = std::numeric_limits<double>::infinity();
    const double h = (opts_.targetEdgeLength > 0.0) ? opts_.targetEdgeLength : 1.0;
    for (size_t qi = 0; qi < cells_.size(); ++qi) {
        const std::array<int, 4> &q = cells_[qi];
        for (int k = 0; k < 4; ++k) {
            const int a = q[k], b = q[(k + 1) % 4];
            const long long key = static_cast<long long>(std::min(a, b)) * n + std::max(a, b);
            auto it = edgeUse.find(key);
            if (it == edgeUse.end()) {
                edgeUse.emplace(key, std::array<int, 3>{{1, cellMaterial_[qi], 0}});
                const double L = normP(verts_[b] - verts_[a]);
                report_.minEdge = std::min(report_.minEdge, L);
                report_.maxEdge = std::max(report_.maxEdge, L);
                lenSum += L;
                ++lenCount;
                if (L > 0.0) {
                    const double r = std::log(L / h);
                    logSum += r * r;
                    if (std::fabs(r) > std::fabs(std::log(report_.worstEdgeRatio)))
                        report_.worstEdgeRatio = L / h;
                }
            } else {
                ++it->second[0];
                if (it->second[1] != cellMaterial_[qi]) it->second[2] = 1;
            }
        }
    }
    if (!std::isfinite(report_.minEdge)) report_.minEdge = 0.0;
    report_.meanEdge = lenCount > 0 ? lenSum / lenCount : 0.0;
    report_.edgeRatioRms = lenCount > 0 ? std::sqrt(logSum / lenCount) : 0.0;
    for (const auto &e : edgeUse) {
        if (e.second[0] == 1) ++report_.boundaryEdges;
        else if (e.second[0] == 2) ++report_.interiorEdges;
        else ++report_.nonManifoldEdges;
        if (e.second[2]) ++report_.interfaceEdges;
    }

    std::unordered_set<int> mats;
    for (int m : cellMaterial_) mats.insert(m);
    report_.materials = static_cast<int>(mats.size());

    double sjSum = 0.0;
    report_.minScaledJacobian = std::numeric_limits<double>::infinity();
    report_.minQuadArea = std::numeric_limits<double>::infinity();
    for (const std::array<int, 4> &q : cells_) {
        const Point c[4] = {verts_[q[0]], verts_[q[1]], verts_[q[2]], verts_[q[3]]};
        const double area = quadArea(c);
        const double sj = scaledJacobian(c);
        report_.meshArea += area;
        report_.minQuadArea = std::min(report_.minQuadArea, area);
        report_.maxQuadArea = std::max(report_.maxQuadArea, area);
        report_.minScaledJacobian = std::min(report_.minScaledJacobian, sj);
        sjSum += sj;
        if (!(sj > 0.0) || !(area > 0.0)) ++report_.invertedQuads;
    }
    if (!std::isfinite(report_.minScaledJacobian)) report_.minScaledJacobian = 0.0;
    if (!std::isfinite(report_.minQuadArea)) report_.minQuadArea = 0.0;
    report_.meanScaledJacobian = sjSum / cells_.size();

    // Cracks: two distinct vertices at one point, on a hash grid.
    const double tol = opts_.crackTolerance * report_.modelExtent;
    const double cell = std::max(tol * 4.0, report_.modelExtent * 1e-12);
    std::unordered_map<long long, std::vector<int>> grid;
    auto key = [&](double x, double y) {
        return static_cast<long long>(std::floor(x / cell)) * 73856093LL ^
               static_cast<long long>(std::floor(y / cell)) * 19349663LL;
    };
    for (size_t i = 0; i < verts_.size(); ++i)
        grid[key(verts_[i][0], verts_[i][1])].push_back(static_cast<int>(i));
    for (size_t i = 0; i < verts_.size(); ++i) {
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                auto it = grid.find(key(verts_[i][0] + dx * cell, verts_[i][1] + dy * cell));
                if (it == grid.end()) continue;
                for (int j : it->second) {
                    if (j <= static_cast<int>(i)) continue;
                    if (normP(verts_[j] - verts_[i]) <= tol) ++report_.cracks;
                }
            }
        }
    }

    // Block corners the mesh turns through more than a half turn at: the
    // element there is reversed whatever the interior does, so these are the
    // folds no smoother can take out.
    for (const Block &b : grids_) {
        if (b.ns < 1 || b.nt < 1) continue;
        const int w = b.ns + 1;
        auto id = [&](int i, int j) { return b.vert[static_cast<size_t>(j) * w + i]; };
        const int ci[4] = {0, b.ns, b.ns, 0}, cj[4] = {0, 0, b.nt, b.nt};
        const int ai[4] = {1, b.ns, b.ns - 1, 0}, aj[4] = {0, 1, b.nt, b.nt - 1};
        const int bi[4] = {0, b.ns - 1, b.ns, 1}, bj[4] = {1, 0, b.nt - 1, b.nt};
        for (int k = 0; k < 4; ++k) {
            const Point o = verts_[id(ci[k], cj[k])];
            const Point u = verts_[id(ai[k], aj[k])] - o;
            const Point v = verts_[id(bi[k], bj[k])] - o;
            double a = std::atan2(cross2(u, v), dotP(u, v));
            if (a < 0.0) a += 2.0 * M_PI;
            if (a > M_PI) ++report_.reflexCorners;
        }
    }

    // The boundary error: how far each vertex of a boundary (or interface)
    // macro edge's polyline lies from the chord between the two mesh nodes
    // either side of it. The polyline is the model's own curve between two
    // macrovertices, so this is the whole of what the meshing gave up there.
    for (size_t e = 0; e < D.edges.size(); ++e) {
        const BlockDecomposition::MacroEdge &me = D.edges[e];
        if (!me.boundary && !me.interface) continue;
        const std::vector<int> &nodes = edgeNodes_[e];
        const int N = static_cast<int>(nodes.size()) - 1;
        if (N < 1 || me.points.size() < 3) continue;
        if (me.boundary) report_.boundaryNodes += N;

        std::vector<double> d(me.points.size(), 0.0);
        for (size_t i = 1; i < me.points.size(); ++i)
            d[i] = d[i - 1] + normP(me.points[i] - me.points[i - 1]);
        const double L = d.back();
        double worst = 0.0;
        for (size_t m = 1; m + 1 < me.points.size(); ++m) {
            const int k = (L > 0.0) ? std::min(N - 1, static_cast<int>(std::floor(d[m] / L * N)))
                                    : 0;
            worst = std::max(worst,
                             segmentDistance(me.points[m], verts_[nodes[k]], verts_[nodes[k + 1]]));
        }
        if (me.boundary) report_.boundaryDeviation = std::max(report_.boundaryDeviation, worst);
        else report_.interfaceDeviation = std::max(report_.interfaceDeviation, worst);
    }

    report_.conforming = report_.nonManifoldEdges == 0 && report_.cracks == 0;
    report_.valid = report_.conforming && report_.invertedQuads == 0 &&
                    report_.unmeshedBlocks == 0;

    if (report_.unmeshedBlocks > 0) {
        report_.messages.push_back(
            std::to_string(report_.unmeshedBlocks) +
            " block(s) have a side with no macro edge behind it, or opposite sides that were "
            "given different counts, and were left unmeshed");
    }
    if (report_.invertedQuads > 0) {
        report_.messages.push_back(
            std::to_string(report_.invertedQuads) +
            " element(s) came out with a corner Jacobian that is not positive" +
            (report_.reflexCorners > 0
                 ? ", against " + std::to_string(report_.reflexCorners) +
                       " block corner(s) the mesh turns through more than a half turn at"
                 : std::string("; the block outlines are sound, so these are straight elements "
                               "spanning a curved piece of one")));
    }
    if (report_.nonManifoldEdges > 0) {
        report_.messages.push_back(std::to_string(report_.nonManifoldEdges) +
                                   " edge(s) are used by more than two elements");
    }
    if (report_.cracks > 0) {
        report_.messages.push_back(std::to_string(report_.cracks) +
                                   " pair(s) of distinct vertices are at the same point");
    }
}

// ---------------------------------------------------------------------------
bool BlockQuadMesh::writeOBJ(const std::string &path) const {
    std::ofstream out(path);
    if (!out) return false;
    out.precision(17);
    out << "# block-decomposition quad mesh: " << cells_.size() << " quads on " << verts_.size()
        << " vertices\n";
    for (const Point &p : verts_) out << "v " << p[0] << " " << p[1] << " 0\n";

    // Grouped by material so that `usemtl` appears once per run rather than
    // once per face, which is what mesh::QuadMesh's reader expects.
    std::vector<int> order(cells_.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                     [&](int a, int b) { return cellMaterial_[a] < cellMaterial_[b]; });
    int current = std::numeric_limits<int>::min();
    for (int q : order) {
        if (cellMaterial_[q] != current) {
            current = cellMaterial_[q];
            out << "usemtl mat" << current << "\n";
        }
        out << "f " << cells_[q][0] + 1 << " " << cells_[q][1] + 1 << " " << cells_[q][2] + 1
            << " " << cells_[q][3] + 1 << "\n";
    }
    return static_cast<bool>(out);
}
