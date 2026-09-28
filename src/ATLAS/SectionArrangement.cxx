#include "ATLAS/SectionArrangement.hxx"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace {

inline double wrap2pi(double a) {
    a = std::fmod(a, 2.0 * M_PI);
    if (a < 0.0) a += 2.0 * M_PI;
    return a;
}

// Counter-clockwise angle from direction u to direction w, in [0, 2 pi).
inline double ccwAngle(const Point &u, const Point &w) {
    return wrap2pi(std::atan2(cross2(u, w), dotP(u, w)));
}

inline Point rotate(const Point &d, double phi) {
    const double c = std::cos(phi), s = std::sin(phi);
    return Point{c * d[0] - s * d[1], s * d[0] + c * d[1]};
}

inline double orient(const Point &a, const Point &b, const Point &c) { return cross2(b - a, c - a); }

// Closed-segment intersection with a relative tolerance: touching counts.
bool segmentsTouch(const Point &a, const Point &b, const Point &c, const Point &d, double eps) {
    const double o1 = orient(a, b, c), o2 = orient(a, b, d);
    const double o3 = orient(c, d, a), o4 = orient(c, d, b);
    auto sgn = [&](double x) { return x > eps ? 1 : (x < -eps ? -1 : 0); };
    const int s1 = sgn(o1), s2 = sgn(o2), s3 = sgn(o3), s4 = sgn(o4);
    if (s1 * s2 < 0 && s3 * s4 < 0) return true;
    auto within = [](const Point &p, const Point &q, const Point &r) {
        return std::min(p[0], q[0]) - 1e-15 <= r[0] && r[0] <= std::max(p[0], q[0]) + 1e-15 &&
               std::min(p[1], q[1]) - 1e-15 <= r[1] && r[1] <= std::max(p[1], q[1]) + 1e-15;
    };
    if (s1 == 0 && within(a, b, c)) return true;
    if (s2 == 0 && within(a, b, d)) return true;
    if (s3 == 0 && within(c, d, a)) return true;
    if (s4 == 0 && within(c, d, b)) return true;
    return false;
}

struct UF {
    std::vector<int> p;
    explicit UF(int n) : p(n) { std::iota(p.begin(), p.end(), 0); }
    int find(int x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }
    void unite(int a, int b) { a = find(a); b = find(b); if (a != b) p[a] = b; }
};

double scaleOf(const std::vector<SectionArrangement::Loop> &L) {
    Point lo = L[0].X[0], hi = L[0].X[0];
    for (const auto &D : L) {
        for (const Point &p : D.X) {
            lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
            hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
        }
    }
    return normP(hi - lo);
}

inline double dirAngle(const Point &d) { return std::atan2(d[1], d[0]); }

} // namespace

bool SectionArrangement::adjacent(const std::vector<Loop> &L, int l, int i, int j) {
    const int N = static_cast<int>(L[l].X.size());
    const int d = ((i - j) % N + N) % N;
    return d == 0 || d == 1 || d == N - 1;
}

bool SectionArrangement::clear(const std::vector<Loop> &L, End a, End b) {
    const double scale = scaleOf(L);
    const double eps = 1e-12 * scale * scale;
    const Point Pa = L[a.l].X[a.i], Pb = L[b.l].X[b.i];
    const Point d = Pb - Pa;
    if (!(normP(d) > 0.0)) return false;
    auto inSector = [&](End e, const Point &dir) {
        const Loop &D = L[e.l];
        const int N = static_cast<int>(D.X.size());
        const double x = ccwAngle(D.X[(e.i + 1) % N] - D.X[e.i], dir);
        return x > 1e-3 && x < D.angle[e.i] - 1e-3;
    };
    if (!inSector(a, d) || !inSector(b, d * -1.0)) return false;
    for (int l2 = 0; l2 < static_cast<int>(L.size()); ++l2) {
        const Loop &D = L[l2];
        const int N = static_cast<int>(D.X.size());
        for (int j = 0; j < N; ++j) {
            const int j1 = (j + 1) % N;
            if ((l2 == a.l && (j == a.i || j1 == a.i)) || (l2 == b.l && (j == b.i || j1 == b.i))) continue;
            if (segmentsTouch(Pa, Pb, D.X[j], D.X[j1], eps)) return false;
        }
    }
    return true;
}

bool SectionArrangement::cast(const std::vector<Loop> &L, const Options &o, std::vector<Section> &out,
                              std::string &why) {
    out.clear();
    const int NL = static_cast<int>(L.size());
    if (NL == 0) { why = "no loop"; return false; }
    const double scale = scaleOf(L);
    // The first loop edge a ray from (l, i) meets, and where it ends.
    auto ray = [&](int l, int i, const Point &dir, End &hit, bool &snapped) {
        const Point O = L[l].X[i];
        double bestT = 1e300, bs = 0.0;
        int bl = -1, bj = -1;
        for (int l2 = 0; l2 < NL; ++l2) {
            const Loop &D = L[l2];
            const int N = static_cast<int>(D.X.size());
            for (int j = 0; j < N; ++j) {
                const int j1 = (j + 1) % N;
                if (l2 == l && (j == i || j1 == i)) continue;
                const Point e = D.X[j1] - D.X[j];
                const double den = cross2(dir, e);
                if (std::fabs(den) < 1e-300) continue;
                const Point w = D.X[j] - O;
                const double t = cross2(w, e) / den;
                const double s = cross2(w, dir) / den;
                if (t > 1e-12 * scale && s >= -1e-12 && s <= 1.0 + 1e-12 && t < bestT) {
                    bestT = t; bs = s; bl = l2; bj = j;
                }
            }
        }
        if (bl < 0) return false;
        const Point H = O + dir * bestT;
        // A node, or a preferred vertex, near the hit takes the section.
        int tl = -1, ti = -1;
        double bestD = o.nodeSnap * bestT;
        for (int l2 = 0; l2 < NL; ++l2) {
            const Loop &D = L[l2];
            for (int j = 0; j < static_cast<int>(D.X.size()); ++j) {
                const bool target = D.node[j] || (!D.preferred.empty() && D.preferred[j]);
                if (!target || (l2 == l && j == i)) continue;
                const double dj = normP(D.X[j] - H);
                if (dj < bestD) { bestD = dj; tl = l2; ti = j; }
            }
        }
        snapped = tl >= 0 && L[tl].node[ti];
        if (tl < 0) {
            tl = bl;
            ti = bs < 0.5 ? bj : (bj + 1) % static_cast<int>(L[bl].X.size());
        }
        hit = {tl, ti};
        return true;
    };
    auto same = [](const Section &x, End a, End b) {
        return (x.a.l == a.l && x.a.i == a.i && x.b.l == b.l && x.b.i == b.i) ||
               (x.a.l == b.l && x.a.i == b.i && x.b.l == a.l && x.b.i == a.i);
    };
    for (int l = 0; l < NL; ++l) {
        const Loop &D = L[l];
        const int N = static_cast<int>(D.X.size());
        for (int i = 0; i < N; ++i) {
            if (!D.node[i]) continue;
            // A corner wants round(angle / quarter turn) cells, and sends a
            // section for each past the first.
            const int q = std::max(1, static_cast<int>(std::lround(D.angle[i] / M_PI_2)));
            const bool corner = D.corner[i] && q >= 2;
            const bool flat = !D.corner[i] && o.flatSections;
            if (!corner && !flat) continue;
            std::vector<double> dirs;
            if (corner) {
                for (int k = 1; k < q; ++k) dirs.push_back(k * D.angle[i] / q);
            } else {
                dirs.push_back(0.5 * D.angle[i]);
            }
            const Point d0 = normalizeP(D.X[(i + 1) % N] - D.X[i]);
            for (double phi : dirs) {
                Point dir = rotate(d0, phi);
                if (o.cross) {
                    // The field's direction nearest the rule's, sampled a
                    // little way in along it.
                    double w = 0.0;
                    const double th = o.cross(D.X[i] + dir * (0.01 * scale), &w);
                    if (w >= o.minCoherence) {
                        Point bestDir = dir;
                        double bestTurn = o.maxTurn;
                        for (int k = 0; k < 4; ++k) {
                            const Point c{std::cos(th + k * M_PI_2), std::sin(th + k * M_PI_2)};
                            const double turn = std::acos(std::max(-1.0, std::min(1.0, dotP(c, dir))));
                            const double inside = ccwAngle(d0, c);
                            if (turn < bestTurn && inside > 0.15 && inside < D.angle[i] - 0.15) {
                                bestTurn = turn;
                                bestDir = c;
                            }
                        }
                        dir = bestDir;
                    }
                }
                End hit;
                bool snapped = false;
                const bool ok = ray(l, i, dir, hit, snapped) &&
                                !(hit.l == l && adjacent(L, l, i, hit.i)) && clear(L, {l, i}, hit);
                // A corner's section has to be there; a flat node that cannot
                // send one stays a plain boundary vertex of its face.
                if (!ok) {
                    if (corner) {
                        why = "a section from a corner leaves the region";
                        return false;
                    }
                    continue;
                }
                bool dup = false;
                for (const Section &x : out) dup = dup || same(x, {l, i}, hit);
                if (!dup) out.push_back({{l, i}, hit, !snapped});
            }
        }
    }
    return true;
}

bool SectionArrangement::arrange(const std::vector<Loop> &L, const std::vector<Section> &S, const Options &o,
                                 Result &r) {
    r = Result();
    auto no = [&](const std::string &why) {
        r.why = why;
        return false;
    };
    const int NL = static_cast<int>(L.size());
    const double scale = scaleOf(L);
    std::vector<std::vector<int>> gOf(NL);
    std::vector<std::vector<char>> isEnd(NL);
    for (int l = 0; l < NL; ++l) {
        gOf[l].assign(L[l].X.size(), -1);
        isEnd[l].assign(L[l].X.size(), 0);
    }
    UF joined(NL);
    for (const Section &s : S) {
        isEnd[s.a.l][s.a.i] = isEnd[s.b.l][s.b.i] = 1;
        joined.unite(s.a.l, s.b.l);
    }
    for (int l = 1; l < NL; ++l) {
        if (joined.find(l) != joined.find(0)) return no("a hole is not joined to the outer loop by a section");
    }
    for (int l = 0; l < NL; ++l) {
        for (int i = 0; i < static_cast<int>(L[l].X.size()); ++i) {
            if (!L[l].corner[i] && !isEnd[l][i]) continue;
            gOf[l][i] = static_cast<int>(r.nodes.size());
            r.nodes.push_back({l, i, L[l].X[i]});
        }
    }
    // Crossings.
    struct Hit { double t; int g; };
    std::vector<std::vector<Hit>> on(S.size());
    for (size_t a = 0; a < S.size(); ++a) {
        on[a].push_back({0.0, gOf[S[a].a.l][S[a].a.i]});
        on[a].push_back({1.0, gOf[S[a].b.l][S[a].b.i]});
    }
    for (size_t a = 0; a < S.size(); ++a) {
        for (size_t b = a + 1; b < S.size(); ++b) {
            const Point p = L[S[a].a.l].X[S[a].a.i], rr = L[S[a].b.l].X[S[a].b.i] - p;
            const Point q = L[S[b].a.l].X[S[b].a.i], ss = L[S[b].b.l].X[S[b].b.i] - q;
            const int ga0 = gOf[S[a].a.l][S[a].a.i], ga1 = gOf[S[a].b.l][S[a].b.i];
            const int gb0 = gOf[S[b].a.l][S[b].a.i], gb1 = gOf[S[b].b.l][S[b].b.i];
            if (ga0 == gb0 || ga0 == gb1 || ga1 == gb0 || ga1 == gb1) {
                const int shared = (ga0 == gb0 || ga0 == gb1) ? ga0 : ga1;
                const Point da = (shared == ga0) ? rr : rr * -1.0;
                const Point db = (shared == gb0) ? ss : ss * -1.0;
                const double x = ccwAngle(da, db);
                if (x < 0.2 || 2.0 * M_PI - x < 0.2) return no("two sections leave one node almost together");
                continue;
            }
            const double den = cross2(rr, ss);
            if (std::fabs(den) < 1e-14 * normP(rr) * normP(ss)) {
                // Parallel: on one line they overlap only if their spans do
                // (two sections from opposite corners of a plate, one each
                // side of its hole, are collinear and disjoint).
                if (std::fabs(cross2(q - p, rr)) < 1e-12 * normP(rr) * scale) {
                    const double L2 = dotP(rr, rr);
                    const double u0 = dotP(q - p, rr) / L2, u1 = dotP(q + ss - p, rr) / L2;
                    if (std::max(u0, u1) > 1e-9 && std::min(u0, u1) < 1.0 - 1e-9) return no("two sections overlap");
                }
                continue;
            }
            const double t = cross2(q - p, ss) / den;
            const double u = cross2(q - p, rr) / den;
            if (t <= -1e-9 || t >= 1.0 + 1e-9 || u <= -1e-9 || u >= 1.0 + 1e-9) continue;
            if (t < 0.02 || t > 0.98 || u < 0.02 || u > 0.98) return no("two sections cross next to an end");
            if (std::fabs(den) / (normP(rr) * normP(ss)) < std::sin(15.0 * M_PI / 180.0)) {
                return no("two sections cross at a grazing angle");
            }
            const int g = static_cast<int>(r.nodes.size());
            r.nodes.push_back({-1, -1, p + rr * t});
            on[a].push_back({t, g});
            on[b].push_back({u, g});
        }
    }
    // Arcs.
    for (int l = 0; l < NL; ++l) {
        const Loop &D = L[l];
        const int N = static_cast<int>(D.X.size());
        std::vector<int> gp;
        for (int i = 0; i < N; ++i) if (gOf[l][i] >= 0) gp.push_back(i);
        if (gp.empty()) return no("a loop has no node");
        for (size_t k = 0; k < gp.size(); ++k) {
            Arc arc;
            arc.loop = l;
            const int s0 = gp[k], e0 = gp[(k + 1) % gp.size()];
            arc.a = gOf[l][s0];
            arc.b = gOf[l][e0];
            int i = s0;
            arc.pos.push_back(i);
            do {
                arc.length += normP(D.X[(i + 1) % N] - D.X[i]);
                i = (i + 1) % N;
                arc.pos.push_back(i);
            } while (i != e0);
            r.arcs.push_back(std::move(arc));
        }
    }
    for (size_t a = 0; a < S.size(); ++a) {
        auto &hits = on[a];
        std::sort(hits.begin(), hits.end(), [](const Hit &x, const Hit &y) { return x.t < y.t; });
        for (size_t k = 0; k + 1 < hits.size(); ++k) {
            Arc arc;
            arc.a = hits[k].g;
            arc.b = hits[k + 1].g;
            arc.section = static_cast<int>(a);
            arc.length = normP(r.nodes[arc.b].x - r.nodes[arc.a].x);
            r.arcs.push_back(arc);
        }
    }
    // Half-edges: loop arcs one way (the region on their left), sections both.
    for (int k = 0; k < static_cast<int>(r.arcs.size()); ++k) {
        const Arc &arc = r.arcs[k];
        if (arc.loop >= 0) {
            const Loop &D = L[arc.loop];
            const int s0 = arc.pos.front(), e0 = arc.pos.back();
            const int s1 = arc.pos[1], e1 = arc.pos[arc.pos.size() - 2];
            r.half.push_back({k, false, arc.a, arc.b, dirAngle(D.X[s1] - D.X[s0]), dirAngle(D.X[e1] - D.X[e0])});
        } else {
            const Point d = r.nodes[arc.b].x - r.nodes[arc.a].x;
            r.half.push_back({k, false, arc.a, arc.b, dirAngle(d), dirAngle(d * -1.0)});
            r.half.push_back({k, true, arc.b, arc.a, dirAngle(d * -1.0), dirAngle(d)});
        }
    }
    std::vector<std::vector<int>> outgoing(r.nodes.size());
    for (int h = 0; h < static_cast<int>(r.half.size()); ++h) outgoing[r.half[h].from].push_back(h);
    auto nextHalf = [&](int h) {
        const double ref = r.half[h].inAng;
        int best = -1;
        double bestD = 1e300;
        for (int x : outgoing[r.half[h].to]) {
            double d = wrap2pi(ref - r.half[x].outAng);
            if (d < 1e-12) d = 2.0 * M_PI;
            if (d < bestD) { bestD = d; best = x; }
        }
        return best;
    };
    std::vector<char> used(r.half.size(), 0);
    for (int h0 = 0; h0 < static_cast<int>(r.half.size()); ++h0) {
        if (used[h0]) continue;
        Face f;
        int h = h0, guard = 0;
        do {
            if (used[h] || ++guard > static_cast<int>(r.half.size()) + 1) return no("the arrangement does not close into faces");
            used[h] = 1;
            f.halves.push_back(h);
            h = nextHalf(h);
            if (h < 0) return no("the arrangement has a dead end");
        } while (h != h0);
        const int m = static_cast<int>(f.halves.size());
        for (int k = 0; k < m; ++k) {
            const Half &prev = r.half[f.halves[(k + m - 1) % m]];
            const Half &cur = r.half[f.halves[k]];
            double a = wrap2pi(prev.inAng - cur.outAng);
            if (a < 1e-12) a = 2.0 * M_PI;
            const Node &g = r.nodes[cur.from];
            const bool required = g.l >= 0 && L[g.l].corner[g.i];
            if (required || std::fabs(a - M_PI) > o.cornerTolerance) f.corners.push_back(k);
        }
        const int K = static_cast<int>(f.corners.size());
        if (K < 3 || K > 5) return no("a face has " + std::to_string(K) + " corners");
        for (int c = 0; c < K; ++c) {
            if ((f.corners[(c + 1) % K] - f.corners[c] + m) % m != 1) {
                return no("a node lies inside a face side (it would be a T-junction)");
            }
        }
        r.faces.push_back(std::move(f));
    }
    r.ok = true;
    return true;
}
