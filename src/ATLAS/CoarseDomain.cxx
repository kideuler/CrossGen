#include "ATLAS/CoarseDomain.hxx"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <utility>

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
                         double minAngle, const std::vector<std::array<int, 2>> &frozen) {
    int flips = 0;
    auto key = [](int a, int b) { return (static_cast<long long>(std::min(a, b)) << 32) | std::max(a, b); };
    std::unordered_set<long long> locked;
    for (const auto &s : frozen) locked.insert(key(s[0], s[1]));
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
                if (locked.count(key(v, w))) continue;                       // an interface segment
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
    buildInterfaceArcs();
    std::vector<std::vector<int>> samples;
    if (sample(samples) && triangulate(samples)) {
        buildField();
        report_.valid = true;
    }
    report_.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// Index every vertex of the arc just pushed, and measure it.
void CoarseDomain::indexArc(int id) {
    const Mesh &M = fine_.getMesh();
    Arc &a = arcs_[id];
    const int n = static_cast<int>(a.vertices.size());
    const long long NV = static_cast<long long>(M.vertices.size());
    a.s.assign(n, 0.0);
    double acc = 0.0;
    for (int i = 0; i < n; ++i) {
        a.s[i] = acc;
        if (a.closed || i + 1 < n) {
            const int w = a.vertices[(i + 1) % n];
            const double len = normP(M.vertices[w] - M.vertices[a.vertices[i]]);
            featAdj_[a.vertices[i]].push_back({w, len});
            featAdj_[w].push_back({a.vertices[i], len});
            acc += len;
        }
        arcIndex_.emplace(static_cast<long long>(id) * NV + a.vertices[i], i);
        arcsAt_[a.vertices[i]].push_back(id);
    }
    a.length = acc;
}

void CoarseDomain::buildArcs() {
    const Mesh &M = fine_.getMesh();
    arcIndex_.clear();
    arcsAt_.assign(M.vertices.size(), {});
    featAdj_.assign(M.vertices.size(), {});
    for (const PlanarDomain::Loop &L : fine_.loops) {
        Arc a;
        a.vertices = L.vertices;
        a.edges = L.edges;
        a.closed = true;
        arcs_.push_back(std::move(a));
        indexArc(static_cast<int>(arcs_.size()) - 1);
    }
    numBoundaryArcs_ = static_cast<int>(arcs_.size());
}

// ---------------------------------------------------------------------------
// The interface network as arcs: one chain between each pair of its nodes.
//
// Stage 1 protects every interface vertex that is not a plain two-edge point
// of one chain -- junctions, landings on dS, kinks, dangling ends -- so the
// chains are exactly the runs of unprotected two-edge vertices between them,
// and a component with no protected vertex at all (a circular inclusion that
// neither meets dS nor turns sharply) is one closed chain.
//
// Each chain is walked so that the fine material on its left and on its right
// are well defined; floodMaterials() seeds the coarse triangles from them.
// ---------------------------------------------------------------------------
void CoarseDomain::buildInterfaceArcs() {
    const Mesh &M = fine_.getMesh();
    if (fine_.getReport().interfaceEdges == 0) return;
    const int NV = static_cast<int>(M.vertices.size());

    // Interface edges at each vertex.
    std::vector<std::vector<int>> at(NV);
    std::vector<int> edgeList;
    for (int e = 0; e < static_cast<int>(fine_.interfaceEdge.size()); ++e) {
        if (!fine_.interfaceEdge[e]) continue;
        at[M.edges[e][0]].push_back(e);
        at[M.edges[e][1]].push_back(e);
        edgeList.push_back(e);
    }
    auto isNode = [&](int v) { return at[v].size() != 2 || fine_.protectedVertex[v]; };

    std::vector<char> used(M.edges.size(), 0);
    // Walk from `v` along `e`, stopping at the next node (or back at the
    // start, for a chain with none).
    auto walk = [&](int v0, int e0) {
        Arc a;
        a.closed = false;
        a.onInterface = true;
        a.vertices.push_back(v0);
        int v = v0, e = e0;
        while (true) {
            used[e] = 1;
            a.edges.push_back(e);
            const int w = M.edges[e][0] == v ? M.edges[e][1] : M.edges[e][0];
            if (w == v0) {                       // came back round: a closed chain
                a.closed = true;
                break;
            }
            a.vertices.push_back(w);
            if (isNode(w)) break;
            const int next = at[w][0] == e ? at[w][1] : at[w][0];
            v = w;
            e = next;
        }
        arcs_.push_back(std::move(a));
        indexArc(static_cast<int>(arcs_.size()) - 1);
    };

    for (int v = 0; v < NV; ++v) {
        if (!isNode(v)) continue;
        for (int e : at[v]) if (!used[e]) walk(v, e);
    }
    // Whatever is left is a chain with no node on it: start anywhere.
    for (int e : edgeList) {
        if (!used[e]) walk(M.edges[e][0], e);
    }

    // The materials either side, from the first edge of each chain: the
    // triangle on the left of vertices[0] -> vertices[1] gives leftMaterial.
    for (int a = numBoundaryArcs_; a < static_cast<int>(arcs_.size()); ++a) {
        Arc &A = arcs_[a];
        ++report_.interfaceArcs;
        if (A.edges.empty() || A.vertices.size() < 2) continue;
        const int e = A.edges[0];
        const int t0 = M.edgeTriangles[e][0], t1 = M.edgeTriangles[e][1];
        if (t0 < 0 || t1 < 0) continue;
        // t0 is on the left of vertices[0] -> vertices[1] when the edge runs
        // that way inside t0's counter-clockwise corner order.
        const int p = A.vertices[0], q = A.vertices[1];
        bool t0Left = false;
        for (int k = 0; k < 3; ++k) {
            if (M.triangles[t0][k] == p && M.triangles[t0][(k + 1) % 3] == q) t0Left = true;
        }
        A.leftMaterial = M.triangleMatId[t0Left ? t0 : t1];
        A.rightMaterial = M.triangleMatId[t0Left ? t1 : t0];
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
    // An open arc has no vertex before its first or after its last: the
    // windows below stop there instead of wrapping, and the missing segments
    // have length zero.
    auto wrapi = [&](int i) { return ((i % n) + n) % n; };
    auto inside = [&](int i) { return A.closed || (i >= 0 && i < n); };
    auto P = [&](int i) { return M.vertices[A.vertices[wrapi(i)]]; };
    auto seg = [&](int i) {
        if (!A.closed && (i < 0 || i >= n - 1)) return 0.0;
        return normP(P(i + 1) - P(i));
    };

    std::vector<double> h(n, hmax);

    // Curvature: the turning of the curve per unit length over a window, the
    // protected corners excluded (they are chain ends, not curvature). On dS
    // the turn is the input's own layout angle; along an interface chain it is
    // the polyline's, since an interior interface vertex has no such angle.
    std::vector<double> turn(n, 0.0), share(n, 0.0);
    for (int i = 0; i < n; ++i) {
        const int v = A.vertices[i];
        if (fine_.protectedVertex[v] || !inside(i - 1) || !inside(i + 1)) {
            turn[i] = 0.0;
        } else if (A.onInterface) {
            const Point a = P(i - 1) - P(i), b = P(i + 1) - P(i);
            turn[i] = std::fabs(M_PI - std::fabs(std::atan2(cross2(a, b), dotP(a, b))));
        } else {
            turn[i] = std::fabs(M_PI - fine_.targetAngle[v]);
        }
        share[i] = 0.5 * (seg(i - 1) + seg(i));
    }
    const double W = std::min(0.5 * hmax, 0.25 * A.length);
    for (int i = 0; i < n; ++i) {
        double T = turn[i], S = share[i];
        double d = 0.0;
        for (int k = 1; k < n; ++k) {           // forward
            if (!inside(i + k)) break;
            d += seg(i + k - 1);
            if (d > W) break;
            const int j = wrapi(i + k);
            T += turn[j];
            S += share[j];
        }
        d = 0.0;
        for (int k = 1; k < n; ++k) {           // backward
            if (!inside(i - k)) break;
            d += seg(i - k);
            if (d > W) break;
            const int j = wrapi(i - k);
            T += turn[j];
            S += share[j];
        }
        if (T > 1e-9 && S > 0.0) h[i] = std::min(h[i], opts_.curvatureFraction * S / T);
    }

    // Gap: the radius of the largest disk inside the domain that touches the
    // curve at this vertex, |w|^2 / (2 w.n) minimised over the other feature
    // vertices in front of it -- across a narrow passage, half its width.
    // Only vertices far from this one *along the feature network* count (more
    // than three times the chord away): near a convex corner the disk shrinks
    // to nothing, but a corner is not a passage -- a block fits a right angle
    // as it is -- and counting it would grade the samples down to the input's
    // own spacing there, making the coarse domain depend on the input
    // triangulation (the dependence Sec. 13.3's last row asks to be rid of).
    //
    // The distance is the geodesic through dS *and* the interfaces, not the
    // arc length of one loop, because the network's curves meet: an interface
    // landing is a corner of the material region either side of it, and
    // measured across the gap rather than along the network it would shrink
    // the samples to nothing on both curves at every landing. Within one
    // closed loop the geodesic is the arc length the other way round, so on a
    // single-material domain this is the rule it replaces, unchanged.
    //
    // The walk stops at 3 hmax / gapFraction: past that a candidate cannot
    // bind (it would ask for a spacing above hmax), so being unreachable and
    // being far are the same answer there.
    std::vector<double> dist;
    std::vector<int> touched;
    const double cap = 3.0 * hmax / std::max(opts_.gapFraction, 1e-6);
    for (int i = 0; i < n; ++i) {
        const Point x = P(i);
        // One-sided at the ends of an open chain, where there is no vertex on
        // the other side to difference against.
        Point t = P(inside(i + 1) ? i + 1 : i) - P(inside(i - 1) ? i - 1 : i);
        const double tl = normP(t);
        if (tl <= 0.0) continue;
        t = t / tl;
        const Point nrm{-t[1], t[0]};
        geodesic(A.vertices[i], cap, dist, touched);
        double r = std::numeric_limits<double>::infinity();
        for (const Arc &B : arcs_) {
            for (int f : B.vertices) {
                if (f == A.vertices[i]) continue;
                const Point w = M.vertices[f] - x;
                if (dist[f] <= 3.0 * normP(w)) continue;
                // An interface has the domain on both sides, so the passage
                // may be either way; dS has it on one.
                const double dn = A.onInterface ? std::fabs(dotP(w, nrm)) : dotP(w, nrm);
                if (dn <= 1e-12 * scale) continue;
                r = std::min(r, dotP(w, w) / (2.0 * dn));
            }
        }
        if (std::isfinite(r)) h[i] = std::min(h[i], opts_.gapFraction * 2.0 * r);
    }

    for (double &x : h) x = std::max(x, hmin);
    // Grading, round the loop twice each way -- and on an open chain not past
    // its ends, which are not neighbours.
    for (int pass = 0; pass < 2; ++pass) {
        for (int k = 0; k < 2 * n; ++k) {
            const int i = k % n, j = (k + 1) % n;
            if (!A.closed && i == n - 1) continue;
            h[j] = std::min(h[j], h[i] + opts_.grading * seg(i));
        }
        for (int k = 2 * n; k > 0; --k) {
            const int i = k % n, j = (k - 1) % n;
            if (!A.closed && j == n - 1) continue;
            h[j] = std::min(h[j], h[i] + opts_.grading * seg(j));
        }
    }
    return h;
}

// ---------------------------------------------------------------------------
// Shortest path from one fine feature vertex through the feature network
// (dS and the interfaces), stopping at `cap`. `dist` is infinity everywhere
// else, and `touched` lists what to reset before the next call.
// ---------------------------------------------------------------------------
void CoarseDomain::geodesic(int source, double cap, std::vector<double> &dist,
                            std::vector<int> &touched) const {
    const double inf = std::numeric_limits<double>::infinity();
    if (dist.size() != fine_.getMesh().vertices.size()) {
        dist.assign(fine_.getMesh().vertices.size(), inf);
        touched.clear();
    }
    for (int v : touched) dist[v] = inf;
    touched.clear();

    typedef std::pair<double, int> Item;   // (distance, vertex)
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> q;
    dist[source] = 0.0;
    touched.push_back(source);
    q.push({0.0, source});
    while (!q.empty()) {
        const Item top = q.top();
        q.pop();
        if (top.first > dist[top.second] + 1e-15) continue;
        if (top.first > cap) break;
        for (const auto &e : featAdj_[top.second]) {
            const double d = top.first + e.second;
            if (d > cap || d >= dist[e.first]) continue;
            if (dist[e.first] == inf) touched.push_back(e.first);
            dist[e.first] = d;
            q.push({d, e.first});
        }
    }
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
    // The input's own section arrangements first: where their sections end,
    // the chains are cut, so that each face side is a whole number of chains
    // whose counts harmonise() can set.
    breakpoint_.assign(M.vertices.size(), 0);
    const bool reconcile = opts_.harmoniseCounts && fine_.getReport().interfaceEdges > 0;
    if (reconcile) planSections();
    // The chains of every arc and their counts, so that the counts can be
    // reconciled across the interfaces (harmonise()) before any sample is
    // placed.
    std::vector<std::vector<SampleChain>> arcChains(arcs_.size());
    for (int l = 0; l < static_cast<int>(arcs_.size()); ++l) {
        const Arc &A = arcs_[l];
        const int n = static_cast<int>(A.vertices.size());
        // A closed curve needs three vertices to be a polygon; an open chain
        // can legitimately be one straight segment between two of its nodes.
        if (A.closed ? n < 3 : n < 2) {
            report_.reason = A.closed ? "a closed feature curve has fewer than three vertices"
                                      : "an interface chain has fewer than two vertices";
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
        for (int i = 0; i < n; ++i) {
            const int v = A.vertices[i];
            if (fine_.protectedVertex[v] || breakpoint_[v]) corners.push_back(i);
        }
        std::vector<SampleChain> &chains = arcChains[l];
        auto add = [&](int a, int len, bool closed) {
            SampleChain c;
            c.a = a;
            c.len = len;
            c.closed = closed;
            chains.push_back(std::move(c));
        };
        if (!A.closed) {
            // An open chain runs between two nodes of the interface network,
            // both protected, and any protected vertex between them (a kink)
            // cuts it further. Its ends are not neighbours, so nothing wraps.
            if (corners.size() < 2 || corners.front() != 0 || corners.back() != n - 1) {
                corners.clear();
                corners.push_back(0);
                corners.push_back(n - 1);
            }
            for (size_t k = 0; k + 1 < corners.size(); ++k) add(corners[k], corners[k + 1] - corners[k], false);
        } else if (corners.empty()) {
            add(0, n, true);
        } else {
            for (size_t k = 0; k < corners.size(); ++k) {
                const int a = corners[k];
                const int b = corners[(k + 1) % corners.size()];
                const int len = corners.size() == 1 ? n : ((b - a) % n + n) % n;
                add(a, len, corners.size() == 1);
            }
        }
        int total = 0;
        for (SampleChain &c : chains) {
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
        // A closed curve needs three samples to be a polygon at all; an open
        // chain is already one with its two ends.
        while (A.closed && total < 3) {
            SampleChain *best = nullptr;
            for (SampleChain &c : chains) {
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
    }

    if (reconcile) harmonise(arcChains);

    for (int l = 0; l < static_cast<int>(arcs_.size()); ++l) {
        const Arc &A = arcs_[l];
        const int n = static_cast<int>(A.vertices.size());
        const std::vector<SampleChain> &chains = arcChains[l];
        for (const SampleChain &c : chains) {
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
        // On a closed curve the next chain's start closes the last one; an
        // open chain has to carry its own far end.
        if (!A.closed) samples[l].push_back(A.vertices[n - 1]);
        report_.chains += static_cast<int>(chains.size());
        report_.samples += static_cast<int>(samples[l].size());
    }
    return true;
}

// ---------------------------------------------------------------------------
// Counts across an interface (see the class comment), part one: the input's
// material regions, each one's loops, and the arrangement of the sections its
// corners and flat nodes cast (SectionArrangement), built on the input's own
// curves. Where a section ends, the sampling cuts its chain (breakpoint_), so
// that each side of each face is a whole number of chains.
// ---------------------------------------------------------------------------
void CoarseDomain::planSections() {
    const Mesh &M = fine_.getMesh();
    const int NT = static_cast<int>(M.triangles.size());
    const int NV = static_cast<int>(M.vertices.size());
    region_.assign(NT, -1);
    regions_ = 0;
    for (int s = 0; s < NT; ++s) {
        if (region_[s] >= 0) continue;
        std::vector<int> stack{s};
        region_[s] = regions_;
        while (!stack.empty()) {
            const int t = stack.back();
            stack.pop_back();
            for (int k = 0; k < 3; ++k) {
                const int e = M.triangleEdges[t][k];
                if (fine_.interfaceEdge[e] || fine_.boundaryEdge[e]) continue;
                const int u = M.edgeTriangles[e][0] == t ? M.edgeTriangles[e][1] : M.edgeTriangles[e][0];
                if (u < 0 || region_[u] >= 0) continue;
                region_[u] = regions_;
                stack.push_back(u);
            }
        }
        ++regions_;
    }
    std::vector<std::vector<int>> tris(regions_);
    for (int t = 0; t < NT; ++t) tris[region_[t]].push_back(t);
    plans_.assign(regions_, RegionPlan());
    const double flatTol = 25.0 * M_PI / 180.0;
    std::vector<int> next(NV, -1);
    std::vector<double> angle(NV, 0.0);
    for (int r = 0; r < regions_; ++r) {
        // The boundary half-edges with the region on their left: a triangle's
        // own edge (triangles are counter-clockwise) whose other side is not
        // the region.
        std::vector<int> starts;
        bool pinched = false;
        for (int t : tris[r]) {
            for (int k = 0; k < 3; ++k) {
                const int e = M.triangleEdges[t][k];
                const int u = M.edgeTriangles[e][0] == t ? M.edgeTriangles[e][1] : M.edgeTriangles[e][0];
                if (u >= 0 && region_[u] == r) continue;
                const int a = M.triangles[t][k], b = M.triangles[t][(k + 1) % 3];
                if (next[a] >= 0) pinched = true;
                next[a] = b;
                starts.push_back(a);
            }
        }
        RegionPlan &plan = plans_[r];
        if (!pinched) {
            std::vector<char> seen;
            for (int a : starts) {
                if (next[a] < 0) continue;
                std::vector<int> loop;
                int v = a;
                while (next[v] >= 0) {
                    loop.push_back(v);
                    const int w = next[v];
                    next[v] = -1;
                    v = w;
                }
                if (v != a) { pinched = true; break; }
                plan.loops.push_back(std::move(loop));
            }
        }
        for (int a : starts) next[a] = -1;
        if (pinched || plan.loops.empty()) {
            plan.loops.clear();
            continue;
        }
        // The region's own interior angle at each loop vertex.
        for (int t : tris[r]) {
            for (int k = 0; k < 3; ++k) {
                const int v = M.triangles[t][k];
                const Point p = M.vertices[v];
                const Point u = M.vertices[M.triangles[t][(k + 1) % 3]] - p;
                const Point w = M.vertices[M.triangles[t][(k + 2) % 3]] - p;
                angle[v] += std::fabs(std::atan2(cross2(u, w), dotP(u, w)));
            }
        }
        // Outer loop first: the one of positive area.
        std::vector<SectionArrangement::Loop> L;
        std::vector<std::vector<int>> ordered;
        int outer = -1;
        for (size_t k = 0; k < plan.loops.size(); ++k) {
            double a2 = 0.0;
            const auto &lp = plan.loops[k];
            for (size_t i = 0; i < lp.size(); ++i) a2 += cross2(M.vertices[lp[i]], M.vertices[lp[(i + 1) % lp.size()]]);
            if (a2 > 0.0) {
                if (outer >= 0) { outer = -2; break; }
                outer = static_cast<int>(k);
            }
        }
        if (outer >= 0) {
            ordered.push_back(plan.loops[outer]);
            for (size_t k = 0; k < plan.loops.size(); ++k) if (static_cast<int>(k) != outer) ordered.push_back(plan.loops[k]);
            for (const auto &lp : ordered) {
                SectionArrangement::Loop D;
                for (int v : lp) {
                    D.X.push_back(M.vertices[v]);
                    D.angle.push_back(angle[v]);
                    const bool node = fine_.protectedVertex[v];
                    D.node.push_back(node ? 1 : 0);
                    D.corner.push_back(node && std::fabs(angle[v] - M_PI) >= flatTol ? 1 : 0);
                }
                L.push_back(std::move(D));
            }
            plan.loops = ordered;
            SectionArrangement::Options so;
            so.cross = opts_.cross;
            std::vector<SectionArrangement::Section> S;
            std::string why;
            bool ok = SectionArrangement::cast(L, so, S, why) && SectionArrangement::arrange(L, S, so, plan.arrangement);
            // A plate round one inclusion, the plate four-cornered and the
            // inclusion smooth, casts nothing that reaches the hole; Stage 3
            // gives it an annulus whose four sector rays run through the
            // plate's corners, so those are the sections to plan: from each
            // corner towards the hole's centre, to where they meet it.
            if (!ok && L.size() == 2) {
                std::vector<int> outerCorners;
                bool holeCorner = false;
                for (size_t i = 0; i < L[0].X.size(); ++i) if (L[0].corner[i]) outerCorners.push_back(static_cast<int>(i));
                for (size_t i = 0; i < L[1].X.size(); ++i) holeCorner = holeCorner || L[1].corner[i];
                if (outerCorners.size() == 4 && !holeCorner) {
                    Point c{0.0, 0.0};
                    for (const Point &p : L[1].X) c = c + p;
                    c = c / static_cast<double>(L[1].X.size());
                    S.clear();
                    bool all = true;
                    for (int i : outerCorners) {
                        const Point O = L[0].X[i];
                        const Point d = c - O;
                        // The first hole edge the ray from the corner to c
                        // crosses, and the nearer end of it.
                        int best = -1;
                        double bestT = 1e300;
                        const int NH = static_cast<int>(L[1].X.size());
                        for (int j = 0; j < NH; ++j) {
                            const Point a = L[1].X[j], e = L[1].X[(j + 1) % NH] - a;
                            const double den = cross2(d, e);
                            if (std::fabs(den) < 1e-300) continue;
                            const Point w = a - O;
                            const double t = cross2(w, e) / den, u = cross2(w, d) / den;
                            if (t <= 0.0 || u < 0.0 || u > 1.0 || t >= bestT) continue;
                            bestT = t;
                            best = u < 0.5 ? j : (j + 1) % NH;
                        }
                        SectionArrangement::End a{0, i}, b{1, best};
                        if (best < 0 || !SectionArrangement::clear(L, a, b)) { all = false; break; }
                        S.push_back({a, b, false});
                    }
                    ok = all && SectionArrangement::arrange(L, S, so, plan.arrangement);
                }
            }
            if (ok) {
                for (const SectionArrangement::Section &sc : S) {
                    breakpoint_[plan.loops[sc.b.l][sc.b.i]] = 1;
                    breakpoint_[plan.loops[sc.a.l][sc.a.i]] = 1;
                }
            }
        }
        for (int t : tris[r]) for (int k = 0; k < 3; ++k) angle[M.triangles[t][k]] = 0.0;
    }
    // A breakpoint only matters on a feature curve; the protected ones cut
    // chains anyway.
    for (int v = 0; v < NV; ++v) {
        if (breakpoint_[v] && (arcsAt_[v].empty() || fine_.protectedVertex[v])) breakpoint_[v] = 0;
    }
}

// ---------------------------------------------------------------------------
// Counts across an interface, part two: every face of every region's
// arrangement asks for its opposite sides equal (a grid) or its sides to
// satisfy a star's triangle inequalities, and through the union-find of the
// grids those become classes of sides that must share one count. Only
// interface chains are raised -- a boundary chain is a floor Stage 3 can lift
// -- to a fixpoint, and a constraint the caps leave unmet is dropped and the
// whole thing redone without it, so that no partial raise is left behind.
// ---------------------------------------------------------------------------
void CoarseDomain::harmonise(std::vector<std::vector<SampleChain>> &chains) {
    const Mesh &M = fine_.getMesh();
    struct Ref {
        int arc = -1, k = -1;
        bool iface = false;
        int cap = 0;
        double length = 0.0;
    };
    std::vector<Ref> G;
    std::vector<int> edgeChain(M.edges.size(), -1);
    for (int l = 0; l < static_cast<int>(arcs_.size()); ++l) {
        const Arc &A = arcs_[l];
        const int n = static_cast<int>(A.vertices.size());
        for (int k = 0; k < static_cast<int>(chains[l].size()); ++k) {
            const SampleChain &c = chains[l][k];
            Ref r;
            r.arc = l;
            r.k = k;
            r.iface = A.onInterface;
            const int gid = static_cast<int>(G.size());
            for (int j = 0; j < c.len; ++j) {
                const int e = A.edges[(c.a + j) % n];
                edgeChain[e] = gid;
                r.length += normP(M.vertices[M.edges[e][1]] - M.vertices[M.edges[e][0]]);
            }
            r.cap = std::min(c.len, std::max(c.segs, static_cast<int>(std::floor(opts_.maxHarmoniseFactor * c.segs))));
            G.push_back(r);
        }
    }
    auto segs = [&](int g) -> int & { return chains[G[g].arc][G[g].k].segs; };

    // A side: the chains along one arc of an arrangement, in order.
    struct Side {
        std::vector<int> chains;
        bool loop = false;       // a loop arc (else a section: free)
        double length = 0.0;
    };
    auto edgeBetween = [&](int a, int b) {
        const auto range = M.vertexTriangles.trianglesForVertex(a);
        for (const int *t = range.first; t != range.second; ++t) {
            for (int k = 0; k < 3; ++k) {
                const int e = M.triangleEdges[*t][k];
                if ((M.edges[e][0] == a && M.edges[e][1] == b) || (M.edges[e][0] == b && M.edges[e][1] == a)) return e;
            }
        }
        return -1;
    };
    // Constraints: classes of sides that must share a count, and stars.
    std::vector<std::vector<Side>> classes;
    std::vector<std::array<Side, 3>> stars;
    int constrainedRegions = 0;
    for (const RegionPlan &plan : plans_) {
        const SectionArrangement::Result &ar = plan.arrangement;
        if (!ar.ok) continue;
        const int NA = static_cast<int>(ar.arcs.size());
        std::vector<Side> sides(NA);
        bool broken = false;
        for (int k = 0; k < NA && !broken; ++k) {
            const SectionArrangement::Arc &arc = ar.arcs[k];
            sides[k].length = arc.length;
            if (arc.loop < 0) continue;
            sides[k].loop = true;
            const std::vector<int> &lp = plan.loops[arc.loop];
            for (size_t j = 0; j + 1 < arc.pos.size(); ++j) {
                const int e = edgeBetween(lp[arc.pos[j]], lp[arc.pos[j + 1]]);
                const int g = e >= 0 ? edgeChain[e] : -1;
                if (g < 0) { broken = true; break; }
                if (sides[k].chains.empty() || sides[k].chains.back() != g) sides[k].chains.push_back(g);
            }
        }
        if (broken) continue;
        struct UFs {
            std::vector<int> p;
            explicit UFs(int n) : p(n) { for (int i = 0; i < n; ++i) p[i] = i; }
            int find(int x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }
            void unite(int a, int b) { a = find(a); b = find(b); if (a != b) p[a] = b; }
        } uf(NA);
        bool any = false;
        for (const SectionArrangement::Face &f : ar.faces) {
            if (f.corners.size() == 4) {
                uf.unite(ar.half[f.halves[f.corners[0]]].arc, ar.half[f.halves[f.corners[2]]].arc);
                uf.unite(ar.half[f.halves[f.corners[1]]].arc, ar.half[f.halves[f.corners[3]]].arc);
            } else if (f.corners.size() == 3) {
                std::array<Side, 3> st;
                for (int c = 0; c < 3; ++c) st[c] = sides[ar.half[f.halves[f.corners[c]]].arc];
                stars.push_back(st);
                any = true;
            }
        }
        std::map<int, std::vector<int>> members;
        for (int k = 0; k < NA; ++k) members[uf.find(k)].push_back(k);
        for (const auto &kv : members) {
            std::vector<Side> cls;
            double lo = 1e300, hi = 0.0;
            bool fixed = false;
            for (int k : kv.second) {
                if (!sides[k].loop) continue;
                cls.push_back(sides[k]);
                lo = std::min(lo, sides[k].length);
                hi = std::max(hi, sides[k].length);
                for (int g : sides[k].chains) fixed = fixed || G[g].iface;
            }
            // Only a class of two or more loop sides, one an interface, and
            // none much longer than another: one grid between a 0.65 arc and a
            // 2.9 side (multimat/icf's big region) is not a layout anyone
            // wants, whatever the counts, and equating them raised the thin
            // shells behind it out of step with each other.
            if (cls.size() < 2 || !fixed || hi > opts_.maxHarmoniseAspect * lo) continue;
            classes.push_back(cls);
            any = true;
        }
        if (any) ++constrainedRegions;
    }
    report_.harmonisedRegions = constrainedRegions;
    if (classes.empty() && stars.empty()) return;

    std::vector<int> before(G.size());
    for (int g = 0; g < static_cast<int>(G.size()); ++g) before[g] = segs(g);
    auto count = [&](const Side &sd) {
        int n = 0;
        for (int g : sd.chains) n += segs(g);
        return n;
    };
    auto fixedSide = [&](const Side &sd) {
        for (int g : sd.chains) if (G[g].iface) return true;
        return false;
    };
    // Raise a side's count by `d`, shared among its interface chains by
    // length so that they keep their spacing, each within its cap.
    auto raiseSide = [&](const Side &sd, int d) {
        if (d <= 0) return false;
        double L = 0.0;
        for (int g : sd.chains) if (G[g].iface) L += G[g].length;
        if (!(L > 0.0)) return false;
        bool changed = false;
        int left = d;
        for (size_t k = 0; k < sd.chains.size() && left > 0; ++k) {
            const int g = sd.chains[k];
            if (!G[g].iface) continue;
            bool last = true;
            for (size_t j = k + 1; j < sd.chains.size(); ++j) last = last && !G[sd.chains[j]].iface;
            const int want = last ? left : std::min(left, static_cast<int>(std::lround(d * G[g].length / L)));
            const int add = std::min(want, std::max(0, G[g].cap - segs(g)));
            if (add > 0) {
                segs(g) += add;
                left -= add;
                changed = true;
            }
        }
        return changed;
    };
    // What constraint q still asks: raises of the sides that fall short.
    // Classes are q < classes.size(), stars after.
    const int NQ = static_cast<int>(classes.size() + stars.size());
    auto apply = [&](int q, bool dryRun) {
        bool unmet = false, changed = false;
        if (q < static_cast<int>(classes.size())) {
            const std::vector<Side> &cls = classes[q];
            // Every interface side up to the largest count in the class; a
            // boundary side only sets the floor.
            int target = 0;
            for (const Side &sd : cls) target = std::max(target, count(sd));
            for (const Side &sd : cls) {
                if (!fixedSide(sd)) continue;
                const int d = target - count(sd);
                if (d <= 0) continue;
                unmet = true;
                if (!dryRun) changed = raiseSide(sd, d) || changed;
            }
        } else {
            // Carrier counts are 2 segs; s_i = (n_j + n_k - n_i) / 2 >= 1
            // binds when j and k are both fixed.
            const std::array<Side, 3> &st = stars[q - classes.size()];
            for (int i = 0; i < 3; ++i) {
                const Side &sj = st[(i + 1) % 3], &sk = st[(i + 2) % 3];
                if (!sj.loop || !sk.loop || !st[i].loop || !fixedSide(sj) || !fixedSide(sk)) continue;
                const int deficit = count(st[i]) + 1 - count(sj) - count(sk);
                if (deficit <= 0) continue;
                unmet = true;
                if (dryRun) continue;
                double lj = 0.0, lk = 0.0;
                for (int g : sj.chains) if (G[g].iface) lj += G[g].length;
                for (int g : sk.chains) if (G[g].iface) lk += G[g].length;
                const int dj = static_cast<int>(std::lround(deficit * lj / std::max(1e-300, lj + lk)));
                changed = raiseSide(sj, dj) || changed;
                changed = raiseSide(sk, deficit - dj) || changed;
            }
        }
        return dryRun ? unmet : changed;
    };
    std::vector<char> active(NQ, 1);
    for (int attempt = 0; attempt <= NQ; ++attempt) {
        for (int g = 0; g < static_cast<int>(G.size()); ++g) segs(g) = before[g];
        for (int pass = 0; pass < 100; ++pass) {
            bool changed = false;
            for (int q = 0; q < NQ; ++q) if (active[q]) changed = apply(q, false) || changed;
            if (!changed) break;
        }
        int unmet = -1;
        for (int q = 0; q < NQ && unmet < 0; ++q) if (active[q] && apply(q, true)) unmet = q;
        if (unmet < 0) break;
        active[unmet] = 0;
    }
    for (int g = 0; g < static_cast<int>(G.size()); ++g) {
        if (segs(g) == before[g]) continue;
        ++report_.harmonisedChains;
        report_.harmonisedSegments += segs(g) - before[g];
    }
}

// ---------------------------------------------------------------------------
// The materials of the coarse triangles.
//
// Not by locating each coarse triangle in the input: in the sliver between a
// sampled chord and the arc it cuts that answer is the wrong material, and the
// coarse interface would then not be the polyline that was constrained. The
// chains carry the answer themselves -- the input material on each side of
// each one is known from the input triangles across its first edge -- so the
// triangle on the left of a sampled segment takes the left material, the one
// on its right the right, and the rest of each region follows across the
// unconstrained edges. Every coarse triangle on one side is then of one
// material by construction, whatever the sampling did to the geometry.
//
// A region the chains do not reach is a material component bounded by dS
// alone; it keeps the material of the input triangle at one of its vertices.
// ---------------------------------------------------------------------------
bool CoarseDomain::floodMaterials(const std::vector<Point> &V, const std::vector<Triangle> &T,
                                  const std::vector<std::array<int, 2>> &segments,
                                  const std::vector<int> &segArc, std::vector<int> &matOut) {
    const Mesh &M = fine_.getMesh();
    const int NT = static_cast<int>(T.size());
    const long long numVertices = static_cast<long long>(V.size());
    const int fallback = M.triangleMatId.empty() ? 1 : M.triangleMatId[0];
    matOut.assign(NT, fallback);
    if (segments.empty()) return true;

    auto key = [&](int a, int b) {
        return static_cast<long long>(std::min(a, b)) * numVertices + std::max(a, b);
    };
    // The directed corner each triangle presents to each edge, so "the
    // triangle on the left of a -> b" is a lookup.
    std::unordered_map<long long, std::array<int, 2>> across;   // edge -> its triangles
    std::unordered_map<long long, int> leftOf;                  // directed a -> b: the triangle
    for (int t = 0; t < NT; ++t) {
        for (int k = 0; k < 3; ++k) {
            const int a = T[t][k], b = T[t][(k + 1) % 3];
            leftOf.emplace(static_cast<long long>(a) * numVertices + b, t);
            auto it = across.find(key(a, b));
            if (it == across.end()) across.emplace(key(a, b), std::array<int, 2>{t, -1});
            else it->second[1] = t;
        }
    }
    std::unordered_set<long long> blocked;
    for (const auto &s : segments) blocked.insert(key(s[0], s[1]));

    std::vector<int> seed(NT, 0);
    for (size_t i = 0; i < segments.size(); ++i) {
        const Arc &A = arcs_[segArc[i]];
        const int a = segments[i][0], b = segments[i][1];
        auto L = leftOf.find(static_cast<long long>(a) * numVertices + b);
        auto R = leftOf.find(static_cast<long long>(b) * numVertices + a);
        if (L == leftOf.end() || R == leftOf.end()) {
            report_.reason = "an interface segment is not an edge of the coarse triangulation";
            return false;
        }
        seed[L->second] = A.leftMaterial;
        seed[R->second] = A.rightMaterial;
    }

    // Flood each region from whatever seeds it contains, stopping at the
    // constrained segments.
    std::vector<int> region(NT, -1);
    std::vector<int> stack;
    int regions = 0;
    for (int t0 = 0; t0 < NT; ++t0) {
        if (region[t0] >= 0) continue;
        const int r = regions++;
        int mat = 0;
        std::vector<int> members;
        region[t0] = r;
        stack.assign(1, t0);
        while (!stack.empty()) {
            const int t = stack.back();
            stack.pop_back();
            members.push_back(t);
            if (seed[t] != 0) {
                if (mat != 0 && mat != seed[t]) {
                    report_.reason = "the sampled interface does not separate the materials";
                    return false;
                }
                mat = seed[t];
            }
            for (int k = 0; k < 3; ++k) {
                const int a = T[t][k], b = T[t][(k + 1) % 3];
                if (blocked.count(key(a, b))) continue;
                auto it = across.find(key(a, b));
                if (it == across.end()) continue;
                for (int u : it->second) {
                    if (u < 0 || region[u] >= 0) continue;
                    region[u] = r;
                    stack.push_back(u);
                }
            }
        }
        if (mat == 0) {
            // No chain touches this region: it is a whole material component
            // of the input, bounded by dS alone, so the input's material at
            // any point of it will do.
            mat = fallback;
            const Point c = (V[T[t0][0]] + V[T[t0][1]] + V[T[t0][2]]) / 3.0;
            for (int f = 0; f < static_cast<int>(M.triangles.size()); ++f) {
                const Point &a = M.vertices[M.triangles[f][0]], &b = M.vertices[M.triangles[f][1]],
                            &d = M.vertices[M.triangles[f][2]];
                if (cross2(b - a, c - a) < 0.0 || cross2(d - b, c - b) < 0.0 || cross2(a - d, c - d) < 0.0) continue;
                mat = M.triangleMatId[f];
                break;
            }
        }
        for (int t : members) matOut[t] = mat;
    }
    report_.materialRegions = regions;
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
    // A sample may be on two curves at once -- a landing is on dS and on its
    // interface chain, a junction on each of its branches -- and Triangle
    // needs one point for it, not one per curve.
    std::unordered_map<int, int> pointOf;
    auto point = [&](int f) {
        auto it = pointOf.find(f);
        if (it != pointOf.end()) return it->second;
        const int id = static_cast<int>(in.vertlist.size());
        in.vertlist.push_back({M.vertices[f][0], M.vertices[f][1]});
        fineOf.push_back(f);
        pointOf.emplace(f, id);
        return id;
    };
    // The interface segments, in their arc's own direction: no edge flip may
    // touch one, and floodMaterials() seeds the triangles either side of them.
    std::vector<std::array<int, 2>> ifaceSeg;
    std::vector<int> ifaceArc;
    for (size_t l = 0; l < samples.size(); ++l) {
        const Arc &A = arcs_[l];
        const int m = static_cast<int>(samples[l].size());
        std::vector<std::array<int, 2>> segs;
        // A closed curve's last sample joins back to its first; an open
        // chain's does not, and Triangle takes it as a bare PSLG segment
        // list, which needs no interior seed point (type 0).
        const int last = A.closed ? m : m - 1;
        for (int k = 0; k < last; ++k) {
            const std::array<int, 2> s{point(samples[l][k]), point(samples[l][(k + 1) % m])};
            segs.push_back(s);
            if (A.onInterface) {
                ifaceSeg.push_back(s);
                ifaceArc.push_back(static_cast<int>(l));
            }
        }
        if (!A.closed) point(samples[l][m - 1]);
        in.segment_loops.push_back(std::move(segs));
        if (!A.closed) {
            in.type.push_back(static_cast<int>(TriangleMesher2D::LoopType::Open));
        } else if (!A.onInterface && !fine_.loops[l].outer) {
            in.type.push_back(static_cast<int>(TriangleMesher2D::LoopType::Hole));
        } else {
            in.type.push_back(static_cast<int>(TriangleMesher2D::LoopType::Exterior));
        }
    }
    report_.interfaceSegments = static_cast<int>(ifaceSeg.size());
    in.h = opts_.maxSpacing * fine_.getReport().scale;
    TriangleMesher2D::Options o;
    o.min_angle_degrees = opts_.minAngle;
    // -Y keeps Triangle off dS; with interfaces the same has to hold of them,
    // which is -YY. Every feature vertex of the coarse mesh must be an input
    // vertex, or Realisation has nowhere to put it back.
    o.suppress_boundary_splitting = true;
    o.suppress_all_splitting = report_.interfaceArcs > 0;
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
        // Only vertices on dS have a boundary valence to chase; an interface
        // sample is an interior vertex, whose carrier valence Stage 6 is the
        // one to rewrite.
        std::vector<int> target(verts.size(), -1);
        for (size_t k = 0; k < fineOf.size() && k < verts.size(); ++k) {
            if (fine_.loopOf[fineOf[k]] < 0) continue;
            const long t = std::lround(fine_.targetAngle[fineOf[k]] / M_PI_2);
            target[k] = static_cast<int>(std::max(1L, std::min(4L, t)));
        }
        // An interface segment is constrained: flipping one would take the
        // interface off the curve it is a sampling of.
        std::vector<std::array<int, 2>> frozen = ifaceSeg;
        report_.flips = flipBoundaryValences(verts, tris, target, opts_.minFlipAngle * M_PI / 180.0, frozen);
    }
    std::vector<int> material;
    if (!floodMaterials(verts, tris, ifaceSeg, ifaceArc, material)) return false;
    mesh_ = std::make_shared<Mesh>(verts, tris, material);

    const int NV = static_cast<int>(verts.size());
    fineVertex_.assign(NV, -1);
    // Triangle keeps the input points first and in order (switch z).
    for (size_t k = 0; k < fineOf.size() && static_cast<int>(k) < NV; ++k) fineVertex_[k] = fineOf[k];
    coarseBreakpoint_.clear();
    bool anyBreakpoint = false;
    for (int k = 0; k < NV; ++k) {
        const int f = fineVertex_[k];
        if (f >= 0 && f < static_cast<int>(breakpoint_.size()) && breakpoint_[f]) anyBreakpoint = true;
    }
    if (anyBreakpoint) {
        coarseBreakpoint_.assign(NV, 0);
        for (int k = 0; k < NV; ++k) {
            const int f = fineVertex_[k];
            if (f >= 0 && breakpoint_[f]) coarseBreakpoint_[k] = 1;
        }
    }

    PlanarDomain::Options po = fine_.getOptions();
    po.angleOverride.assign(NV, std::numeric_limits<double>::quiet_NaN());
    po.interfaceAngleOverride.assign(NV, std::numeric_limits<double>::quiet_NaN());
    po.extraCorners.clear();
    for (int k = 0; k < NV; ++k) {
        const int f = fineVertex_[k];
        if (f < 0) continue;
        po.angleOverride[k] = fine_.targetAngle[f];
        // A sample interior to an interface chain is by construction not one
        // of the input's kinks -- those are protected, and so are chain ends
        // -- so the layout sees the chain running straight through it,
        // whatever the chord the sampling drew turns by.
        if (fine_.loopOf[f] < 0 && fine_.interfaceDegree[f] == 2 && !fine_.protectedVertex[f]) {
            po.interfaceAngleOverride[k] = M_PI;
        }
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
    // Every coarse feature vertex must be a sample: a point Triangle added on
    // a boundary or interface segment would have no place on the input's
    // curve, and Realisation could not put it back.
    for (int v : mesh_->boundaryVertices) {
        if (fineVertex_[v] < 0) {
            report_.reason = "Triangle put a point on a boundary segment";
            return false;
        }
    }
    for (int v = 0; v < NV; ++v) {
        if (domain_->interfaceDegree[v] > 0 && fineVertex_[v] < 0) {
            report_.reason = "Triangle put a point on an interface segment";
            return false;
        }
    }
    // The interface the flood fill produced must be the one that was sampled,
    // segment for segment: any other is a material region the sampling moved.
    if (domain_->getReport().interfaceEdges != report_.interfaceSegments) {
        report_.reason = "the coarse interface has " + std::to_string(domain_->getReport().interfaceEdges) +
                         " edge(s), not the " + std::to_string(report_.interfaceSegments) + " that were sampled";
        return false;
    }

    // Which fine curve each coarse feature edge samples, and which way round.
    edgeFrom_.assign(mesh_->edges.size(), -1);
    edgeArc_.assign(mesh_->edges.size(), -1);
    for (size_t l = 0; l < domain_->loops.size(); ++l) {
        const PlanarDomain::Loop &L = domain_->loops[l];
        for (size_t i = 0; i < L.edges.size(); ++i) {
            edgeFrom_[L.edges[i]] = L.vertices[i];
            edgeArc_[L.edges[i]] = static_cast<int>(l);
        }
        // The coarse loop must run the way its fine loop does.
        const int f0 = fineVertex_[L.vertices[0]];
        if (f0 < 0 || arcIndex(static_cast<int>(l), f0) < 0) {
            report_.reason = "a coarse boundary loop has no fine counterpart";
            return false;
        }
    }
    std::unordered_map<long long, int> edgeOf;
    if (!ifaceSeg.empty()) {
        for (size_t e = 0; e < mesh_->edges.size(); ++e) {
            const int a = mesh_->edges[e][0], b = mesh_->edges[e][1];
            edgeOf.emplace(static_cast<long long>(std::min(a, b)) * NV + std::max(a, b), static_cast<int>(e));
        }
    }
    for (size_t i = 0; i < ifaceSeg.size(); ++i) {
        const int a = ifaceSeg[i][0], b = ifaceSeg[i][1];
        auto it = edgeOf.find(static_cast<long long>(std::min(a, b)) * NV + std::max(a, b));
        if (it == edgeOf.end() || !domain_->interfaceEdge[it->second]) {
            report_.reason = "a sampled interface segment is not an interface edge of the coarse mesh";
            return false;
        }
        edgeFrom_[it->second] = a;
        edgeArc_[it->second] = ifaceArc[i];
    }
    return true;
}

// ---------------------------------------------------------------------------
int CoarseDomain::arcIndex(int arc, int fineVertex) const {
    if (arc < 0 || fineVertex < 0) return -1;
    const long long NV = static_cast<long long>(fine_.getMesh().vertices.size());
    auto it = arcIndex_.find(static_cast<long long>(arc) * NV + fineVertex);
    return it == arcIndex_.end() ? -1 : it->second;
}

double CoarseDomain::positionOf(int arc, int fineVertex) const {
    const int i = arcIndex(arc, fineVertex);
    return i < 0 ? -1.0 : arcs_[arc].s[i];
}

double CoarseDomain::forward(int loop, double a, double b) const {
    const Arc &A = arcs_[loop];
    if (!A.closed) return std::max(0.0, std::min(A.length, b) - std::max(0.0, a));
    double d = wrap(b - a, A.length);
    if (d <= 1e-14 * A.length) d = A.length;
    return d;
}

Point CoarseDomain::pointAt(int loop, double s, int *edge, int *index) const {
    const Arc &A = arcs_[loop];
    const Mesh &M = fine_.getMesh();
    const int n = static_cast<int>(A.vertices.size());
    // An open chain does not wrap: past either end is that end.
    s = A.closed ? wrap(s, A.length) : std::max(0.0, std::min(A.length, s));
    int i = static_cast<int>(std::upper_bound(A.s.begin(), A.s.end(), s) - A.s.begin()) - 1;
    i = std::max(0, std::min(A.closed ? n - 1 : n - 2, i));
    const double s0 = A.s[i];
    const double s1 = (i + 1 < n) ? A.s[i + 1] : A.length;
    const double t = s1 > s0 ? (s - s0) / (s1 - s0) : 0.0;
    if (edge) *edge = A.edges[i];
    if (index) *index = i;
    const Point &p = M.vertices[A.vertices[i]], &q = M.vertices[A.vertices[(i + 1) % n]];
    return p + (q - p) * t;
}

void CoarseDomain::locateAll(const SquareCarrier &C, int v, std::vector<Location> &out) const {
    out.clear();
    const Mesh &M = fine_.getMesh();
    // A sample is its own input vertex, and lies on every curve through it.
    if (C.sourceVertex[v] >= 0) {
        const int f = fineVertex_[C.sourceVertex[v]];
        if (f < 0) return;
        for (int a : arcsAt_[f]) {
            Location L;
            L.loop = a;
            L.s = positionOf(a, f);
            out.push_back(L);
        }
        return;
    }
    // Anything else must sit inside one coarse feature edge, which came from
    // exactly one curve: it goes at the same fraction of the input arc
    // between that edge's two samples.
    const int e = C.sourceEdge[v];
    if (e < 0 || e >= static_cast<int>(edgeFrom_.size()) || edgeFrom_[e] < 0) return;
    const int loop = edgeArc_[e];
    if (loop < 0) return;
    const int a = edgeFrom_[e];
    const int b = mesh_->edges[e][0] == a ? mesh_->edges[e][1] : mesh_->edges[e][0];
    const int fa = fineVertex_[a], fb = fineVertex_[b];
    if (fa < 0 || fb < 0) return;
    const double sa = positionOf(loop, fa), sb = positionOf(loop, fb);
    if (sa < 0.0 || sb < 0.0) return;
    const Point pa = M.vertices[fa], pb = M.vertices[fb];
    const Point d = pb - pa;
    const double L2 = dotP(d, d);
    if (!(L2 > 0.0)) return;
    const double t = std::max(0.0, std::min(1.0, dotP(C.vertices[v] - pa, d) / L2));
    Location L;
    L.loop = loop;
    L.s = sa + t * forward(loop, sa, sb);
    if (arcs_[loop].closed) L.s = wrap(L.s, arcs_[loop].length);
    out.push_back(L);
}

CoarseDomain::Location CoarseDomain::locate(const SquareCarrier &C, int v) const {
    std::vector<Location> all;
    locateAll(C, v, all);
    return all.empty() ? Location() : all.front();
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
        // Every feature curve is aligned, dS and the interfaces alike: Sec.
        // 1.1 asks the layout to run along both. A protected vertex is where
        // the tangent jumps, so it is left free.
        const int arc = (f >= 0 && !fine_.protectedVertex[f] && !arcsAt_[f].empty()) ? arcsAt_[f].front() : -1;
        if (arc >= 0) {
            const Arc &A = arcs_[arc];
            const int n = static_cast<int>(A.vertices.size()), i = arcIndex(arc, f);
            const int next = (i + 1 < n || A.closed) ? (i + 1) % n : i;
            const int prev = (i > 0 || A.closed) ? (i + n - 1) % n : i;
            const Point t = M.vertices[A.vertices[next]] - M.vertices[A.vertices[prev]];
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
