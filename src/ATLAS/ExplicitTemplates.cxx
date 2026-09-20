#include "ATLAS/ExplicitTemplates.hxx"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <map>
#include <numeric>
#include <set>
#include <sstream>

namespace {

typedef CavityFill::Patch Patch;
typedef SquareCarrier::Origin Origin;

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

// Union-find with path halving.
struct UF {
    std::vector<int> p;
    explicit UF(int n) : p(n) { std::iota(p.begin(), p.end(), 0); }
    int find(int x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }
    void unite(int a, int b) { a = find(a); b = find(b); if (a != b) p[a] = b; }
};

// Principal axis of a polygon's boundary, weighted by edge length.
double principalAngle(const std::vector<Point> &pts) {
    const size_t n = pts.size();
    double W = 0.0;
    Point m{0.0, 0.0};
    for (size_t k = 0; k < n; ++k) {
        const Point &a = pts[k], &b = pts[(k + 1) % n];
        const double w = normP(b - a);
        m = m + (a + b) * (0.5 * w);
        W += w;
    }
    if (W <= 0.0) return 0.0;
    m = m / W;
    double sxx = 0.0, syy = 0.0, sxy = 0.0;
    for (size_t k = 0; k < n; ++k) {
        const Point &a = pts[k], &b = pts[(k + 1) % n];
        const double w = normP(b - a);
        const Point d = (a + b) * 0.5 - m;
        sxx += w * d[0] * d[0];
        syy += w * d[1] * d[1];
        sxy += w * d[0] * d[1];
    }
    return 0.5 * std::atan2(2.0 * sxy, sxx - syy);
}

bool collinear(const std::vector<Point> &pts, double tol) {
    if (pts.size() < 3) return true;
    const Point a = pts.front(), b = pts.back();
    const Point d = b - a;
    const double L = normP(d);
    if (L <= 0.0) return false;
    for (const Point &p : pts) if (std::fabs(cross2(d, p - a)) / L > tol * L) return false;
    return true;
}

} // namespace

// ---------------------------------------------------------------------------
ExplicitTemplates::ExplicitTemplates(SquareCarrier &carrier, const Options &opts)
    : C_(carrier), opts_(opts) {
    report_.cellsBefore = C_.numCells();
    std::set<int> tried;
    double domainArea = std::fabs(C_.getDomain().getReport().area);
    int guard = 0;
    while (++guard < 100000) {
        std::vector<Region> regions = extractRegions();
        if (report_.regions == 0) report_.regions = static_cast<int>(regions.size());
        // Regions bounded by dS alone first -- nothing couples them to a
        // neighbour -- then the rest, largest first, so that a template's
        // designated corners on an interface are in place before the region
        // across it is tried.
        std::stable_sort(regions.begin(), regions.end(), [](const Region &a, const Region &b) {
            if (a.touchesInterface != b.touchesInterface) return !a.touchesInterface;
            return a.area > b.area;
        });
        bool committed = false;
        for (const Region &R : regions) {
            if (tried.count(R.key)) continue;
            tried.insert(R.key);
            Attempt A;
            A.region = static_cast<int>(report_.attempts.size());
            A.cells = static_cast<int>(R.cells.size());
            A.material = R.material;
            A.holes = static_cast<int>(R.holes.size());
            if (attempt(R, A)) {
                report_.templated++;
                report_.blocks += A.blocks;
                report_.boundarySplits += A.splits;
                report_.templatedArea += domainArea > 0.0 ? R.area / domainArea : 0.0;
                report_.attempts.push_back(A);
                committed = true;
                break;   // cell ids changed; recompute the regions
            }
            report_.attempts.push_back(A);
        }
        if (!committed) break;
    }
    report_.untouched = report_.regions - report_.templated;
    report_.cellsAfter = C_.numCells();
}

bool ExplicitTemplates::isNode(int v) const {
    return C_.protectedVertex[v] || C_.designatedVertex[v];
}

std::vector<ExplicitTemplates::Region> ExplicitTemplates::extractRegions() const {
    std::vector<Region> out;
    const int NQ = C_.numCells();
    std::vector<char> seen(NQ, 0);
    for (int s = 0; s < NQ; ++s) {
        if (seen[s] || C_.cellOrigin[s] == SquareCarrier::CellOrigin::Template) continue;
        Region R;
        R.material = C_.cellMaterial[s];
        std::vector<int> stack{s};
        seen[s] = 1;
        while (!stack.empty()) {
            const int q = stack.back();
            stack.pop_back();
            R.cells.push_back(q);
            for (int i = 0; i < 4; ++i) {
                const int r = C_.neighbor[q][i];
                if (r < 0 || seen[r]) continue;
                if (C_.cellOrigin[r] == SquareCarrier::CellOrigin::Template) continue;
                if (C_.cellMaterial[r] != R.material) continue;
                seen[r] = 1;
                stack.push_back(r);
            }
        }
        R.set.insert(R.cells.begin(), R.cells.end());
        R.key = 1 << 30;
        for (int q : R.cells) {
            R.area += C_.cellArea(q);
            if (C_.cellTriangle[q] >= 0) R.key = std::min(R.key, C_.cellTriangle[q]);
        }
        R.B = CavityFill::boundaryOf(C_, R.cells, R.set);
        double len = 0.0;
        int cnt = 0;
        for (size_t k = 0; k < R.B.loops.size(); ++k) {
            if (R.B.signedArea[k] > 0.0) {
                if (R.outer < 0) R.outer = static_cast<int>(k);
                else R.outer = -2;   // two outer loops: not a region we template
            } else {
                R.holes.push_back(static_cast<int>(k));
            }
            for (size_t i = 0; i < R.B.loops[k].size(); ++i) {
                const int e = R.B.loopEdges[k][i];
                len += normP(C_.vertices[C_.edges[e][1]] - C_.vertices[C_.edges[e][0]]);
                ++cnt;
                if (!C_.boundaryEdge[e]) R.touchesInterface = true;
            }
        }
        R.h = cnt > 0 ? len / cnt : 0.0;
        out.push_back(std::move(R));
    }
    return out;
}

bool ExplicitTemplates::attempt(const Region &R, Attempt &A) {
    if (!R.B.manifold || R.outer < 0) {
        A.family = "none";
        A.reason = R.B.manifold ? "the region has no single outer loop" : R.B.reason;
        return false;
    }
    const std::vector<int> &L = R.B.loops[R.outer];
    int K = 0, reflex = 0;
    const double reflexTol = opts_.reflexAngle * M_PI / 180.0;
    for (size_t i = 0; i < L.size(); ++i) {
        if (!isNode(L[i])) continue;
        ++K;
        if (R.B.angle[R.outer][i] > M_PI + reflexTol) ++reflex;
    }
    for (int h : R.holes) for (int v : R.B.loops[h]) if (isNode(v)) ++K;
    A.corners = K;

    // The region's cones, for the families' parameters and their scores.
    regionCones_.clear();
    if (opts_.field) {
        std::vector<Point> outer;
        for (int v : L) outer.push_back(C_.vertices[v]);
        std::vector<std::vector<Point>> holes;
        for (int h : R.holes) {
            holes.emplace_back();
            for (int v : R.B.loops[h]) holes.back().push_back(C_.vertices[v]);
        }
        regionCones_ = opts_.field->conesIn(outer, holes);
        for (const ReferenceField::Singularity &s : regionCones_) (s.sign > 0 ? A.conesPlus : A.conesMinus)++;
    }
    const bool choose = opts_.field && opts_.chooseByField;

    std::vector<std::string> reasons;
    auto commit = [&](const Patch &P, Attempt &trial, double minSJ) {
        SquareCarrier::Edit ed;
        ed.removeCells = R.cells;
        ed.newVertices = P.newVertices;
        ed.newOrigin = P.newOrigin;
        ed.newSourceEdge = P.newSourceEdge;
        ed.cells = P.cells;
        ed.material = R.material;
        ed.designate = P.designated;
        ed.cellOrigin = SquareCarrier::CellOrigin::Template;
        ed.group = group_++;
        C_.apply({ed});
        trial.accepted = true;
        trial.newCells = static_cast<int>(P.cells.size());
        trial.minScaledJacobian = minSJ;
        trial.blocks = P.blocks;
        A = trial;
    };
    // Without chooseByField a family that certifies is committed at once, in
    // the fixed order; with it, it is kept and the order only breaks ties.
    struct Built {
        Patch P;
        Attempt trial;
        double minSJ = 0.0;
    };
    std::vector<Built> built;
    auto run = [&](const char *family, bool (ExplicitTemplates::*fn)(const Region &, Patch &, Attempt &)) {
        Patch P;
        Attempt trial = A;
        trial.family = family;
        if (!(this->*fn)(R, P, trial)) {
            reasons.push_back(std::string(family) + ": " + trial.reason);
            return false;
        }
        const CavityFill::Verdict V =
            CavityFill::settle(C_, R.cells, P, opts_.minScaledJacobian,
                               (std::string(family) == "O-grid" || std::string(family) == "half O-grid" ||
                                std::string(family) == "annulus")
                                   ? 0 : opts_.smoothingIterations);
        if (!V.valid) {
            reasons.push_back(std::string(family) + ": " + V.reason);
            return false;
        }
        if (opts_.field) score(R, P, trial);
        if (!choose) {
            commit(P, trial, V.minScaledJacobian);
            return true;
        }
        built.push_back({std::move(P), trial, V.minScaledJacobian});
        return false;
    };
    // The least score among the families built; the earliest on a tie.
    auto chooseBest = [&]() {
        if (built.empty()) return false;
        size_t best = 0;
        for (size_t k = 1; k < built.size(); ++k) {
            if (built[k].trial.score < built[best].trial.score - 1e-9) best = k;
        }
        std::ostringstream alt;
        alt << std::setprecision(3);
        for (size_t k = 0; k < built.size(); ++k) {
            if (k == best) continue;
            if (alt.tellp() > 0) alt << "; ";
            alt << built[k].trial.family << " " << built[k].trial.score;
        }
        built[best].trial.alternatives = alt.str();
        commit(built[best].P, built[best].trial, built[best].minSJ);
        return true;
    };

    if (R.holes.empty()) {
        chords_.clear();
        if (opts_.sections && (reflex > 0 || K == 4) && run("sections", &ExplicitTemplates::trySections)) return true;
        if (opts_.sections && reflex == 0 && K >= 6 && K % 2 == 0) {
            // A fan of chords from each corner in turn, until one certifies.
            std::vector<int> npos;
            for (size_t i = 0; i < L.size(); ++i) if (isNode(L[i])) npos.push_back(static_cast<int>(i));
            std::set<std::vector<std::pair<int, int>>> tried;
            const size_t before = reasons.size();
            for (int s0 = 0; s0 < K; ++s0) {
                std::vector<std::pair<int, int>> fan;
                for (int j = 3; j <= K - 3; j += 2) {
                    const int a = npos[s0], b = npos[(s0 + j) % K];
                    fan.push_back({std::min(a, b), std::max(a, b)});
                }
                std::sort(fan.begin(), fan.end());
                if (!tried.insert(fan).second) continue;
                chords_ = fan;
                const size_t before = built.size();
                if (run("sections", &ExplicitTemplates::trySections)) { chords_.clear(); return true; }
                if (built.size() > before) break;   // the first fan that certifies, as without a field
            }
            chords_.clear();
            // One reason for the whole family is enough.
            if (reasons.size() > before + 1) reasons.erase(reasons.begin() + before, reasons.end() - 1);
        }
        if (opts_.stars && reflex == 0 && (K == 3 || K == 5) && run("star", &ExplicitTemplates::tryStar)) return true;
        // Two real corners want a cell each, which the half O-grid gives them
        // and the O-grid cannot; it declines anything else, and the O-grid is
        // still there behind it.
        if (opts_.halfOGrids && K == 2 && run("half O-grid", &ExplicitTemplates::tryHalfOGrid)) return true;
        if (opts_.ogrids && K <= 2 && run("O-grid", &ExplicitTemplates::tryOGrid)) return true;
        if (chooseBest()) return true;
        if (reasons.empty()) {
            std::ostringstream os;
            os << K << " corner(s), " << reflex << " reflex: no family applies";
            reasons.push_back(os.str());
        }
    } else if (R.holes.size() == 1) {
        if (opts_.annuli && run("annulus", &ExplicitTemplates::tryAnnulus)) return true;
        if (chooseBest()) return true;
        if (reasons.empty()) reasons.push_back("annulus: disabled");
    } else {
        std::ostringstream os;
        os << R.holes.size() << " holes: no single template covers a region of this "
           << "connectivity (Sec. 7.4), the carrier is retained";
        reasons.push_back(os.str());
    }
    A.family = "none";
    std::ostringstream os;
    for (size_t k = 0; k < reasons.size(); ++k) os << (k ? "; " : "") << reasons[k];
    A.reason = os.str();
    return false;
}

// ---------------------------------------------------------------------------
// A certified template against the field, on a copy of the carrier with it
// applied. Only the vertices whose every cell is the template's own count
// towards E_sing: a vertex on an interface also has cells of the region across
// it, which may still be the fine split's, and its valence says nothing about
// this template.
void ExplicitTemplates::score(const Region &R, const Patch &P, Attempt &A) const {
    const ReferenceField &F = *opts_.field;
    SquareCarrier T = C_;
    SquareCarrier::Edit ed;
    ed.removeCells = R.cells;
    ed.newVertices = P.newVertices;
    ed.newOrigin = P.newOrigin;
    ed.newSourceEdge = P.newSourceEdge;
    ed.cells = P.cells;
    ed.material = R.material;
    ed.designate = P.designated;
    ed.cellOrigin = SquareCarrier::CellOrigin::Template;
    const int trialGroup = 1 << 29;
    ed.group = trialGroup;
    T.apply({ed});
    double m = 0.0;
    for (int q = 0; q < T.numCells(); ++q) if (T.cellGroup[q] == trialGroup) m += F.misalignment(T, q);
    std::vector<char> own(T.numVertices(), 0);
    for (int v = 0; v < T.numVertices(); ++v) {
        if (T.valence[v] == 0) continue;
        bool all = true;
        for (int r = T.ringPtr[v]; r < T.ringPtr[v + 1] && all; ++r) all = T.cellGroup[T.ringCell[r]] == trialGroup;
        own[v] = all ? 1 : 0;
    }
    ReferenceField::SingularityReport sr;
    A.eSing = F.singularityEnergy(ReferenceField::carrierSingularities(T, &own), regionCones_, opts_.singularity, &sr);
    A.eDir = m / F.area();
    A.singPlus = sr.interiorPlus;
    A.singMinus = sr.interiorMinus;
    A.edgePlus = sr.boundaryPlus;
    A.edgeMinus = sr.boundaryMinus;
    A.score = P.blocks + opts_.wDir * A.eDir + A.eSing;
    A.scored = true;
}

bool ExplicitTemplates::arcSplittable(const std::vector<int> &ids) const {
    if (!opts_.splitBoundary) return false;
    for (size_t k = 0; k + 1 < ids.size(); ++k) {
        if (ids[k] < 0 || ids[k + 1] < 0) continue;
        const int e = C_.edgeBetween(ids[k], ids[k + 1]);
        if (e >= 0 && C_.boundaryEdge[e]) return true;
    }
    return false;
}

std::vector<int> ExplicitTemplates::splitArc(Patch &P, const std::vector<int> &ids, int extra) const {
    if (extra <= 0) return ids;
    if (!opts_.splitBoundary) return {};
    struct Seg { int a, b, meshEdge; bool splittable; };
    std::vector<Seg> segs;
    for (size_t k = 0; k + 1 < ids.size(); ++k) {
        const int a = ids[k], b = ids[k + 1];
        Seg s{a, b, -1, false};
        if (a >= 0 && b >= 0) {
            const int e = C_.edgeBetween(a, b);
            s.splittable = e >= 0 && C_.boundaryEdge[e];
            s.meshEdge = C_.sourceEdge[a] >= 0 ? C_.sourceEdge[a] : C_.sourceEdge[b];
        }
        segs.push_back(s);
    }
    for (int x = 0; x < extra; ++x) {
        int best = -1;
        double bestL = -1.0;
        for (size_t k = 0; k < segs.size(); ++k) {
            if (!segs[k].splittable) continue;
            const double L = normP(P.at(C_, segs[k].b) - P.at(C_, segs[k].a));
            if (L > bestL) { bestL = L; best = static_cast<int>(k); }
        }
        if (best < 0) return {};
        const Seg s = segs[best];
        const int m = P.addVertex((P.at(C_, s.a) + P.at(C_, s.b)) * 0.5, Origin::BoundarySplit, s.meshEdge);
        segs[best] = Seg{s.a, m, s.meshEdge, true};
        segs.insert(segs.begin() + best + 1, Seg{m, s.b, s.meshEdge, true});
    }
    std::vector<int> out{segs.front().a};
    for (const Seg &s : segs) out.push_back(s.b);
    return out;
}

std::string ExplicitTemplates::classifyGrid(const Patch &P, const std::array<std::vector<int>, 4> &sides) const {
    std::array<bool, 4> straight;
    for (int k = 0; k < 4; ++k) {
        std::vector<Point> pts;
        for (int id : sides[k]) pts.push_back(P.at(C_, id));
        straight[k] = collinear(pts, 1e-9);
    }
    if (straight[0] && straight[1] && straight[2] && straight[3]) return "bilinear";
    if ((straight[1] && straight[3]) || (straight[0] && straight[2])) return "ruled";
    return "Coons";
}

// ---------------------------------------------------------------------------
// Sectioned quadrilaterals: Secs. 7.1 and 7.2.
// ---------------------------------------------------------------------------
bool ExplicitTemplates::trySections(const Region &R, Patch &P, Attempt &A) {
    const std::vector<int> &L = R.B.loops[R.outer];
    const std::vector<double> &ang = R.B.angle[R.outer];
    const std::vector<int> &LE = R.B.loopEdges[R.outer];
    const int N = static_cast<int>(L.size());
    auto mod = [N](int i) { return ((i % N) + N) % N; };
    std::vector<Point> X(N);
    for (int i = 0; i < N; ++i) X[i] = C_.vertices[L[i]];
    Point lo = X[0], hi = X[0];
    for (const Point &p : X) {
        lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
        hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
    }
    const double scale = normP(hi - lo);
    const double eps = 1e-12 * scale * scale;
    std::vector<char> node(N, 0);
    for (int i = 0; i < N; ++i) node[i] = isNode(L[i]);
    const double reflexTol = opts_.reflexAngle * M_PI / 180.0;
    const double cornerTol = opts_.cornerTolerance * M_PI / 180.0;

    auto fail = [&](const std::string &why) { A.reason = why; return false; };

    // Is the chord X[i] -> X[k] inside the region? It must leave both ends
    // strictly inside their interior sectors and touch no loop edge except
    // the ones at its ends.
    auto clear = [&](int i, int k) {
        const Point d = X[k] - X[i];
        const double ai = ccwAngle(X[mod(i + 1)] - X[i], d);
        const double ak = ccwAngle(X[mod(k + 1)] - X[k], X[i] - X[k]);
        const double margin = 1e-3;
        if (!(ai > margin && ai < ang[i] - margin)) return false;
        if (!(ak > margin && ak < ang[k] - margin)) return false;
        for (int j = 0; j < N; ++j) {
            const int j1 = mod(j + 1);
            if (j == i || j1 == i || j == k || j1 == k) continue;
            if (segmentsTouch(X[i], X[k], X[j], X[j1], eps)) return false;
        }
        return true;
    };

    // ---- sections from every reflex corner --------------------------------
    std::set<std::pair<int, int>> sectionSet;
    std::vector<std::pair<int, int>> sections;
    for (int i = 0; i < N; ++i) {
        if (!node[i] || ang[i] <= M_PI + reflexTol) continue;
        const int q = std::max(2, static_cast<int>(std::lround(ang[i] / M_PI_2)));
        const Point d0 = normalizeP(X[mod(i + 1)] - X[i]);
        for (int k = 1; k < q; ++k) {
            const Point dir = rotate(d0, k * ang[i] / q);
            int hitEdge = -1;
            double hitT = 1e300, hitS = 0.0;
            for (int j = 0; j < N; ++j) {
                const int j1 = mod(j + 1);
                if (j == i || j1 == i) continue;
                const Point e = X[j1] - X[j];
                const double den = cross2(dir, e);
                if (std::fabs(den) < 1e-300) continue;
                const Point w = X[j] - X[i];
                const double t = cross2(w, e) / den;
                const double s = cross2(w, dir) / den;
                if (t > 1e-12 * scale && s >= -1e-12 && s <= 1.0 + 1e-12 && t < hitT) {
                    hitT = t; hitS = s; hitEdge = j;
                }
            }
            if (hitEdge < 0) return fail("a section from a reflex corner finds no far boundary");
            const Point H = X[i] + dir * hitT;
            int target = -1;
            double bestD = opts_.nodeSnap * hitT;
            for (int j = 0; j < N; ++j) {
                if (!node[j] || j == i) continue;
                const double dj = normP(X[j] - H);
                if (dj < bestD) { bestD = dj; target = j; }
            }
            if (target < 0) target = hitS < 0.5 ? hitEdge : mod(hitEdge + 1);
            if (target == i || mod(target - i) == 1 || mod(i - target) == 1) {
                return fail("a section would end next to its own corner");
            }
            if (!clear(i, target)) return fail("a section leaves the region");
            const std::pair<int, int> key(std::min(i, target), std::max(i, target));
            if (sectionSet.insert(key).second) sections.push_back({i, target});
        }
    }
    for (const auto &ch : chords_) {
        if (!clear(ch.first, ch.second)) return fail("a chord between corners leaves the region");
        const std::pair<int, int> key(std::min(ch.first, ch.second), std::max(ch.first, ch.second));
        if (sectionSet.insert(key).second) sections.push_back(ch);
    }
    int nodeCount = 0;
    for (int i = 0; i < N; ++i) nodeCount += node[i];
    if (sections.empty() && nodeCount != 4) return fail("not four corners and nothing to cut at");

    // ---- the arrangement: graph nodes, arcs, half-edges -------------------
    struct GNode { int pos = -1; int id = 0; Point x; };
    std::vector<GNode> G;
    std::vector<int> gOfPos(N, -1);
    std::vector<char> isEnd(N, 0);
    for (const auto &s : sections) { isEnd[s.first] = 1; isEnd[s.second] = 1; }
    for (int i = 0; i < N; ++i) {
        if (!node[i] && !isEnd[i]) continue;
        gOfPos[i] = static_cast<int>(G.size());
        G.push_back({i, L[i], X[i]});
    }
    if (G.size() < 2) return fail("fewer than two nodes on the boundary");

    // Crossings between sections.
    struct Hit { double t; int g; };
    std::vector<std::vector<Hit>> onSection(sections.size());
    for (size_t a = 0; a < sections.size(); ++a) {
        onSection[a].push_back({0.0, gOfPos[sections[a].first]});
        onSection[a].push_back({1.0, gOfPos[sections[a].second]});
    }
    for (size_t a = 0; a < sections.size(); ++a) {
        for (size_t b = a + 1; b < sections.size(); ++b) {
            const int a0 = sections[a].first, a1 = sections[a].second;
            const int b0 = sections[b].first, b1 = sections[b].second;
            const Point p = X[a0], r = X[a1] - X[a0];
            const Point q = X[b0], s = X[b1] - X[b0];
            if (a0 == b0 || a0 == b1 || a1 == b0 || a1 == b1) {
                // Sharing an end is a node of the loop; they only must not overlap.
                const int shared = (a0 == b0 || a0 == b1) ? a0 : a1;
                const Point da = (shared == a0) ? r : r * -1.0;
                const Point db = (shared == b0) ? s : s * -1.0;
                if (std::fabs(ccwAngle(da, db)) < 0.2 || std::fabs(ccwAngle(da, db) - 2.0 * M_PI) < 0.2) {
                    return fail("two sections leave one corner almost together");
                }
                continue;
            }
            const double den = cross2(r, s);
            if (std::fabs(den) < 1e-14 * normP(r) * normP(s)) {
                if (std::fabs(cross2(q - p, r)) < 1e-12 * normP(r) * scale) return fail("two sections overlap");
                continue;
            }
            const double t = cross2(q - p, s) / den;
            const double u = cross2(q - p, r) / den;
            if (t <= -1e-9 || t >= 1.0 + 1e-9 || u <= -1e-9 || u >= 1.0 + 1e-9) continue;
            if (t < 0.02 || t > 0.98 || u < 0.02 || u > 0.98) return fail("two sections cross next to an end");
            if (std::fabs(den) / (normP(r) * normP(s)) < std::sin(15.0 * M_PI / 180.0)) {
                return fail("two sections cross at a grazing angle");
            }
            const int g = static_cast<int>(G.size());
            G.push_back({-1, 0, p + r * t});
            onSection[a].push_back({t, g});
            onSection[b].push_back({u, g});
        }
    }

    struct Arc {
        int a = -1, b = -1;       // graph nodes
        bool loop = false;
        std::vector<int> pos;     // loop arcs: loop positions a..b
        int count = 0;
        bool splittable = false;
        double length = 0.0;
    };
    std::vector<Arc> arcs;
    // Loop arcs between consecutive graph nodes.
    std::vector<int> gpos;
    for (int i = 0; i < N; ++i) if (gOfPos[i] >= 0) gpos.push_back(i);
    for (size_t k = 0; k < gpos.size(); ++k) {
        Arc arc;
        arc.loop = true;
        const int s = gpos[k], e = gpos[(k + 1) % gpos.size()];
        arc.a = gOfPos[s];
        arc.b = gOfPos[e];
        int i = s;
        arc.pos.push_back(i);
        bool allBoundary = true;
        do {
            if (!C_.boundaryEdge[LE[i]]) allBoundary = false;
            arc.length += normP(X[mod(i + 1)] - X[i]);
            i = mod(i + 1);
            arc.pos.push_back(i);
        } while (i != e);
        arc.count = static_cast<int>(arc.pos.size()) - 1;
        arc.splittable = allBoundary && opts_.splitBoundary;
        arcs.push_back(arc);
    }
    const int loopArcs = static_cast<int>(arcs.size());
    for (size_t a = 0; a < sections.size(); ++a) {
        auto &hits = onSection[a];
        std::sort(hits.begin(), hits.end(), [](const Hit &x, const Hit &y) { return x.t < y.t; });
        for (size_t k = 0; k + 1 < hits.size(); ++k) {
            Arc arc;
            arc.a = hits[k].g;
            arc.b = hits[k + 1].g;
            arc.length = normP(G[arc.b].x - G[arc.a].x);
            arcs.push_back(arc);
        }
    }

    struct Half { int arc; bool rev; int from, to; double outAng, inAng; };
    std::vector<Half> half;
    auto dirAngle = [](const Point &d) { return std::atan2(d[1], d[0]); };
    for (int k = 0; k < static_cast<int>(arcs.size()); ++k) {
        const Arc &arc = arcs[k];
        if (arc.loop) {
            const int s = arc.pos.front(), e = arc.pos.back();
            const int s1 = arc.pos[1], e1 = arc.pos[arc.pos.size() - 2];
            half.push_back({k, false, arc.a, arc.b, dirAngle(X[s1] - X[s]), dirAngle(X[e1] - X[e])});
        } else {
            const Point d = G[arc.b].x - G[arc.a].x;
            half.push_back({k, false, arc.a, arc.b, dirAngle(d), dirAngle(d * -1.0)});
            half.push_back({k, true, arc.b, arc.a, dirAngle(d * -1.0), dirAngle(d)});
        }
    }
    std::vector<std::vector<int>> outgoing(G.size());
    for (int h = 0; h < static_cast<int>(half.size()); ++h) outgoing[half[h].from].push_back(h);

    auto nextHalf = [&](int h) {
        const int w = half[h].to;
        const double ref = half[h].inAng;
        int best = -1;
        double bestD = 1e300;
        for (int o : outgoing[w]) {
            double d = wrap2pi(ref - half[o].outAng);
            if (d < 1e-12) d = 2.0 * M_PI;
            if (d < bestD) { bestD = d; best = o; }
        }
        return best;
    };

    // ---- faces --------------------------------------------------------------
    std::vector<char> used(half.size(), 0);
    std::vector<std::vector<int>> faces;
    for (int h0 = 0; h0 < static_cast<int>(half.size()); ++h0) {
        if (used[h0]) continue;
        std::vector<int> face;
        int h = h0;
        int guard = 0;
        do {
            if (used[h] || ++guard > static_cast<int>(half.size()) + 1) return fail("the arrangement does not close into faces");
            used[h] = 1;
            face.push_back(h);
            h = nextHalf(h);
            if (h < 0) return fail("the arrangement has a dead end");
        } while (h != h0);
        faces.push_back(face);
    }

    // Each face: four corners, one arc per side.
    struct FaceSides { std::array<int, 4> side; };
    std::vector<FaceSides> fsides;
    for (const auto &face : faces) {
        const int m = static_cast<int>(face.size());
        std::vector<int> corners;   // index into face: the half-edge that starts at a corner
        for (int k = 0; k < m; ++k) {
            const Half &prev = half[face[(k + m - 1) % m]];
            const Half &cur = half[face[k]];
            double a = wrap2pi(prev.inAng - cur.outAng);
            if (a < 1e-12) a = 2.0 * M_PI;
            const int g = cur.from;
            const bool required = G[g].pos >= 0 && node[G[g].pos];
            if (required || std::fabs(a - M_PI) > cornerTol) corners.push_back(k);
        }
        if (corners.size() != 4) {
            std::ostringstream os;
            os << "a face of the sectioned region has " << corners.size() << " corners";
            return fail(os.str());
        }
        FaceSides fs;
        for (int c = 0; c < 4; ++c) {
            const int k0 = corners[c], k1 = corners[(c + 1) % 4];
            const int len = (k1 - k0 + m) % m;
            if (len != 1) return fail("a node lies inside a face side (it would be a T-junction)");
            fs.side[c] = face[k0];
        }
        fsides.push_back(fs);
    }

    // ---- counts: Sec. 11.2's union-find, opposite sides equal ---------------
    UF uf(static_cast<int>(arcs.size()));
    for (const FaceSides &fs : fsides) {
        uf.unite(half[fs.side[0]].arc, half[fs.side[2]].arc);
        uf.unite(half[fs.side[1]].arc, half[fs.side[3]].arc);
    }
    std::map<int, std::vector<int>> classes;
    for (int k = 0; k < static_cast<int>(arcs.size()); ++k) classes[uf.find(k)].push_back(k);
    std::vector<int> target(arcs.size(), 1);
    for (const auto &kv : classes) {
        int fixed = -1, maxCur = 0;
        double idealSum = 0.0;
        int nIdeal = 0;
        bool anyLoop = false;
        for (int k : kv.second) {
            const Arc &arc = arcs[k];
            if (arc.loop) {
                anyLoop = true;
                if (!arc.splittable) {
                    if (fixed >= 0 && fixed != arc.count) return fail("two interface arcs that must match have different counts");
                    fixed = arc.count;
                }
                maxCur = std::max(maxCur, arc.count);
            } else {
                idealSum += std::max(1.0, arc.length / std::max(1e-300, R.h));
                ++nIdeal;
            }
        }
        int n;
        if (fixed >= 0) {
            if (maxCur > fixed) return fail("an interface arc is coarser than the boundary arc it must match");
            n = fixed;
        } else if (anyLoop) {
            n = maxCur;
        } else {
            n = std::max(1, static_cast<int>(std::lround(idealSum / std::max(1, nIdeal))));
        }
        for (int k : kv.second) target[k] = n;
    }

    // ---- the patch ----------------------------------------------------------
    std::vector<int> gid(G.size());
    for (size_t g = 0; g < G.size(); ++g) {
        gid[g] = G[g].pos >= 0 ? G[g].id : P.addVertex(G[g].x, Origin::Template);
    }
    std::vector<std::vector<int>> arcIds(arcs.size());
    for (int k = 0; k < static_cast<int>(arcs.size()); ++k) {
        const Arc &arc = arcs[k];
        if (arc.loop) {
            std::vector<int> ids;
            for (int p : arc.pos) ids.push_back(L[p]);
            if (target[k] > arc.count) {
                const std::vector<int> split = splitArc(P, ids, target[k] - arc.count);
                if (split.empty()) return fail("an arc would need subdividing and cannot be");
                A.splits += target[k] - arc.count;
                ids = split;
            }
            arcIds[k] = ids;
        } else {
            arcIds[k] = CavityFill::segment(C_, P, gid[arc.a], gid[arc.b], target[k], Origin::Template);
        }
    }
    (void)loopArcs;
    std::map<std::string, int> mapCount;
    for (const FaceSides &fs : fsides) {
        std::array<std::vector<int>, 4> sides;
        for (int c = 0; c < 4; ++c) {
            const Half &h = half[fs.side[c]];
            sides[c] = arcIds[h.arc];
            if (h.rev) std::reverse(sides[c].begin(), sides[c].end());
        }
        ++mapCount[classifyGrid(P, sides)];
        if (!CavityFill::grid(C_, P, sides, Origin::Template)) return fail("a face's sides do not close");
    }
    for (size_t g = 0; g < G.size(); ++g) P.designated.push_back(gid[g]);
    P.blocks = static_cast<int>(fsides.size());
    P.kind = "sections";
    std::ostringstream os;
    bool first = true;
    for (const auto &kv : mapCount) { os << (first ? "" : ", ") << kv.second << " " << kv.first; first = false; }
    A.maps = os.str();
    return true;
}

// ---------------------------------------------------------------------------
// Stars: K = 3 or 5 corners, one centre of valence K.
// ---------------------------------------------------------------------------
bool ExplicitTemplates::tryStar(const Region &R, Patch &P, Attempt &A) {
    const std::vector<int> &L = R.B.loops[R.outer];
    const int N = static_cast<int>(L.size());
    auto fail = [&](const std::string &why) { A.reason = why; return false; };
    std::vector<int> nodes;
    for (int i = 0; i < N; ++i) if (isNode(L[i])) nodes.push_back(i);
    const int K = static_cast<int>(nodes.size());
    if (K != 3 && K != 5) return fail("not three or five corners");

    std::vector<std::vector<int>> arc(K);
    std::vector<int> n(K);
    std::vector<char> splittable(K);
    for (int k = 0; k < K; ++k) {
        int i = nodes[k];
        const int e = nodes[(k + 1) % K];
        arc[k].push_back(L[i]);
        do { i = (i + 1) % N; arc[k].push_back(L[i]); } while (i != e);
        n[k] = static_cast<int>(arc[k].size()) - 1;
        splittable[k] = arcSplittable(arc[k]);
        // An arc is subdividable only if all of it is domain boundary.
        for (size_t j = 0; j + 1 < arc[k].size() && splittable[k]; ++j) {
            const int ed = C_.edgeBetween(arc[k][j], arc[k][j + 1]);
            if (ed < 0 || !C_.boundaryEdge[ed]) splittable[k] = 0;
        }
    }

    // The least subdivision that makes n_i = s_{i-1} + s_{i+1} solvable.
    std::vector<int> sigma;
    if (!CavityFill::repairStar(n, splittable, sigma)) {
        return fail("no positive integer spoke counts, even subdividing the domain-boundary sides");
    }
    std::vector<int> bestD(K);
    for (int k = 0; k < K; ++k) bestD[k] = sigma[(k + K - 1) % K] + sigma[(k + 1) % K] - n[k];
    for (int k = 0; k < K; ++k) {
        if (bestD[k] > 0) {
            arc[k] = splitArc(P, arc[k], bestD[k]);
            if (arc[k].empty()) return fail("a side could not be subdivided");
            A.splits += bestD[k];
        }
    }

    std::vector<Point> poly;
    for (int v : L) poly.push_back(C_.vertices[v]);
    Point c = CavityFill::kernelCenter(poly, CavityFill::areaCentroid(poly));
    const double depth = CavityFill::kernelDepth(poly, c);
    if (!(depth > 1e-6 * R.h)) return fail("the region is not star-shaped about any point found");
    // The field's one cone of the centre's sign -- +1/4 for three blocks,
    // -1/4 for five -- is where the centre belongs, if it is well inside the
    // kernel (a quarter of the deepest point's depth).
    if (opts_.field && opts_.conePlacement) {
        const int sign = K == 3 ? 1 : -1;
        const ReferenceField::Singularity *only = nullptr;
        int count = 0;
        for (const ReferenceField::Singularity &s : regionCones_) {
            if (s.sign != sign) continue;
            ++count;
            only = &s;
        }
        if (count == 1 && CavityFill::kernelDepth(poly, only->x) > 0.25 * depth) {
            c = only->x;
            A.conePlaced = true;
        }
    }
    std::vector<int> mid;
    int center = 0;
    if (!CavityFill::star(C_, P, arc, sigma, c, Origin::Template, &mid, &center)) {
        return fail("a star block's sides do not close");
    }
    for (int k = 0; k < K; ++k) { P.designated.push_back(L[nodes[k]]); P.designated.push_back(mid[k]); }
    P.designated.push_back(center);
    A.maps = "star of " + std::to_string(K) + " Coons blocks";
    P.blocks = K;
    P.kind = "star";
    return true;
}

// ---------------------------------------------------------------------------
// O-grid: Sec. 7.3.
// ---------------------------------------------------------------------------
bool ExplicitTemplates::tryOGrid(const Region &R, Patch &P, Attempt &A) {
    const std::vector<int> &L = R.B.loops[R.outer];
    const int N = static_cast<int>(L.size());
    auto fail = [&](const std::string &why) { A.reason = why; return false; };
    if (N < 8 || N % 2 != 0) return fail("too few boundary edges for four shells");
    std::vector<Point> X(N);
    for (int i = 0; i < N; ++i) X[i] = C_.vertices[L[i]];
    const Point c = CavityFill::kernelCenter(X, CavityFill::areaCentroid(X));
    const double depth = CavityFill::kernelDepth(X, c);
    if (!(depth > 1e-3 * R.h)) return fail("not star-shaped about any point found (Sec. 7.3 needs a kernel point)");

    std::vector<int> required;
    for (int i = 0; i < N; ++i) if (isNode(L[i])) required.push_back(i);

    // Ideal split directions: the corners of a rectangle on the principal
    // axes with half-widths in proportion to the region's extents.
    const double phi = principalAngle(X);
    const Point e1{std::cos(phi), std::sin(phi)}, e2{-std::sin(phi), std::cos(phi)};
    double w1 = 0.0, w2 = 0.0;
    for (const Point &p : X) {
        w1 = std::max(w1, std::fabs(dotP(p - c, e1)));
        w2 = std::max(w2, std::fabs(dotP(p - c, e2)));
    }
    std::array<double, 4> ideal;
    const int sx[4] = {1, -1, -1, 1}, sy[4] = {1, 1, -1, -1};
    for (int k = 0; k < 4; ++k) {
        const Point d = e1 * (sx[k] * w1) + e2 * (sy[k] * w2);
        ideal[k] = std::atan2(d[1], d[0]);
    }
    // A region holding exactly four +1/4 cones of the field and no -1/4 is the
    // disk-like case, and there the principal axes are degenerate: on a disk
    // they pick the orientation by rounding. The field's cones sit on the
    // diagonals of the cross at the centre (DualMBO::setPinDiskCenters), so
    // they give the split rays, and the core corners go where they are.
    std::array<Point, 4> cone;
    bool fromCones = false;
    if (opts_.field && opts_.conePlacement) {
        std::vector<Point> plus;
        int minus = 0;
        for (const ReferenceField::Singularity &s : regionCones_) {
            if (s.sign > 0) plus.push_back(s.x);
            else ++minus;
        }
        if (plus.size() == 4 && minus == 0) {
            std::sort(plus.begin(), plus.end(), [&](const Point &u, const Point &w) {
                return std::atan2(u[1] - c[1], u[0] - c[0]) < std::atan2(w[1] - c[1], w[0] - c[0]);
            });
            for (int k = 0; k < 4; ++k) {
                cone[k] = plus[k];
                ideal[k] = std::atan2(plus[k][1] - c[1], plus[k][0] - c[0]);
            }
            fromCones = true;
        }
    }
    std::vector<double> theta(N);
    for (int i = 0; i < N; ++i) theta[i] = std::atan2(X[i][1] - c[1], X[i][0] - c[0]);
    auto angDist = [](double a, double b) {
        double d = std::fabs(wrap2pi(a - b));
        return std::min(d, 2.0 * M_PI - d);
    };

    int bestI0 = -1, bestA = -1, bestRot = 0;
    double bestScore = 1e300;
    const int half = N / 2;
    for (int i0 = 0; i0 < half; ++i0) {
        for (int a = 1; a < half; ++a) {
            const int p[4] = {i0, (i0 + a) % N, (i0 + half) % N, (i0 + half + a) % N};
            bool ok = true;
            for (int r : required) {
                if (r != p[0] && r != p[1] && r != p[2] && r != p[3]) { ok = false; break; }
            }
            if (!ok) continue;
            for (int rot = 0; rot < 4; ++rot) {
                double s = 0.0;
                for (int k = 0; k < 4; ++k) {
                    const double dd = angDist(theta[p[k]], ideal[(k + rot) % 4]);
                    s += dd * dd;
                }
                if (s < bestScore) { bestScore = s; bestI0 = i0; bestA = a; bestRot = rot; }
            }
        }
    }
    if (bestI0 < 0) return fail("the protected corners cannot be four O-grid split points");
    const int a = bestA, b = half - bestA;
    const int p[4] = {bestI0, (bestI0 + a) % N, (bestI0 + half) % N, (bestI0 + half + a) % N};
    const int cnt[4] = {a, b, a, b};

    // Core corners and their convexity.
    // The core's corners go on the split rays -- the shells' radial
    // certificate needs the ray through each -- at Sec. 7.3's fraction, or at
    // the cones' mean radius. One fraction for all four: taking each corner
    // from its own cone skews the core, and on bubbles' disks that cost 0.707
    // -> 0.60 in the template's worst cell for nothing E_sing can see.
    double alpha = opts_.coreFraction;
    if (fromCones) {
        double f = 0.0;
        for (int k = 0; k < 4; ++k) {
            const double r = normP(X[p[k]] - c);
            f += r > 0.0 ? normP(cone[(k + bestRot) % 4] - c) / r / 4.0 : alpha / 4.0;
        }
        alpha = std::clamp(f, 0.3, 0.75);
    }
    std::array<Point, 4> q;
    for (int k = 0; k < 4; ++k) q[k] = c + (X[p[k]] - c) * alpha;
    A.conePlaced = fromCones;
    for (int k = 0; k < 4; ++k) {
        if (!(cross2(q[(k + 1) % 4] - q[k], q[(k + 3) % 4] - q[k]) > 0.0)) return fail("the core quadrilateral is not convex");
    }
    std::array<int, 4> qid;
    for (int k = 0; k < 4; ++k) qid[k] = P.addVertex(q[k], Origin::Template);

    // For each shell: the arc's vertices, and their rays' hits on the core side.
    std::array<std::vector<int>, 4> outer, core;
    for (int k = 0; k < 4; ++k) {
        const Point qa = q[k], qb = q[(k + 1) % 4];
        double lastS = -1.0;
        for (int j = 0; j <= cnt[k]; ++j) {
            const int li = (p[k] + j) % N;
            outer[k].push_back(L[li]);
            if (j == 0) { core[k].push_back(qid[k]); lastS = 0.0; continue; }
            if (j == cnt[k]) { core[k].push_back(qid[(k + 1) % 4]); continue; }
            const Point d = X[li] - c, e = qb - qa;
            const double den = cross2(d, e);
            if (std::fabs(den) < 1e-300) return fail("a ray is parallel to the core side");
            const Point w = qa - c;
            const double t = cross2(w, e) / den;
            const double s = cross2(w, d) / den;
            if (!(t > 0.0 && t < 1.0 && s > lastS && s < 1.0)) return fail("a ray misses its core side");
            lastS = s;
            core[k].push_back(P.addVertex(c + d * t, Origin::Template));
        }
    }

    // Radial layers: cells about as deep as the boundary edges are long.
    double radial = 0.0;
    int rc = 0;
    for (int k = 0; k < 4; ++k) {
        for (size_t j = 0; j < outer[k].size(); ++j) {
            radial += normP(C_.vertices[outer[k][j]] - P.at(C_, core[k][j]));
            ++rc;
        }
    }
    const int m = std::max(1, static_cast<int>(std::lround((radial / std::max(1, rc)) / std::max(1e-300, R.h))));

    // Radial lines, shared by neighbouring shells at the four split rays.
    std::array<std::vector<int>, 4> splitRay;
    for (int k = 0; k < 4; ++k) {
        splitRay[k].push_back(qid[k]);
        const Point a0 = q[k], a1 = X[p[k]];
        for (int j = 1; j < m; ++j) splitRay[k].push_back(P.addVertex(a0 + (a1 - a0) * (static_cast<double>(j) / m), Origin::Template));
        splitRay[k].push_back(L[p[k]]);
    }
    for (int k = 0; k < 4; ++k) {
        const int n = cnt[k];
        std::vector<int> nodes(static_cast<size_t>(n + 1) * (m + 1));
        auto at = [&](int i, int j) -> int & { return nodes[i + (n + 1) * j]; };
        for (int i = 0; i <= n; ++i) {
            at(i, 0) = core[k][i];
            at(i, m) = outer[k][i];
        }
        for (int j = 0; j <= m; ++j) {
            at(0, j) = splitRay[k][j];
            at(n, j) = splitRay[(k + 1) % 4][j];
        }
        for (int i = 1; i < n; ++i) {
            const Point a0 = P.at(C_, core[k][i]), a1 = C_.vertices[outer[k][i]];
            for (int j = 1; j < m; ++j) at(i, j) = P.addVertex(a0 + (a1 - a0) * (static_cast<double>(j) / m), Origin::Template);
        }
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < n; ++i) P.cells.push_back({at(i, j), at(i, j + 1), at(i + 1, j + 1), at(i + 1, j)});
        }
    }

    // The core: a straight-line grid, node (i, j) the intersection of the
    // chord bottom_i -> top_i with the chord left_j -> right_j.
    {
        const int n = a, mm = b;
        std::vector<int> nodes(static_cast<size_t>(n + 1) * (mm + 1));
        auto at = [&](int i, int j) -> int & { return nodes[i + (n + 1) * j]; };
        for (int i = 0; i <= n; ++i) { at(i, 0) = core[0][i]; at(i, mm) = core[2][n - i]; }
        for (int j = 0; j <= mm; ++j) { at(n, j) = core[1][j]; at(0, j) = core[3][mm - j]; }
        for (int j = 1; j < mm; ++j) {
            for (int i = 1; i < n; ++i) {
                const Point b0 = P.at(C_, at(i, 0)), t0 = P.at(C_, at(i, mm));
                const Point l0 = P.at(C_, at(0, j)), r0 = P.at(C_, at(n, j));
                const Point r = t0 - b0, s = r0 - l0;
                const double den = cross2(r, s);
                if (std::fabs(den) < 1e-300) return fail("two core chords are parallel");
                const double t = cross2(l0 - b0, s) / den;
                at(i, j) = P.addVertex(b0 + r * t, Origin::Template);
            }
        }
        for (int j = 0; j < mm; ++j) {
            for (int i = 0; i < n; ++i) P.cells.push_back({at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)});
        }
    }

    for (int k = 0; k < 4; ++k) { P.designated.push_back(L[p[k]]); P.designated.push_back(qid[k]); }
    P.blocks = 5;
    P.kind = "O-grid";
    A.maps = "1 straight-line core, 4 radial";
    return true;
}

// ---------------------------------------------------------------------------
// Half O-grid: Sec. 7.3 mirrored across a straight side and cut back.
// ---------------------------------------------------------------------------
//
// Names, walking the loop with the region on the left: the straight side runs
// from corner A to corner B through s1 and s2, the core's feet; the arc runs
// back from B to A through t2 and t1, where the shells split; q2 and q1 are the
// core's top corners, on the rays from c through t2 and t1. The four blocks
// are the core s1 s2 q2 q1, the shell at B (s2 B t2 q2), the top shell
// (q2 t2 t1 q1) and the shell at A (q1 t1 A s1). Opposite sides must carry
// equal counts, so with m the depth of the shells, S the height of the core
// and W its width,
//
//     |A s1| = |s2 B| = m,   |B t2| = |t1 A| = S,   |s1 s2| = |t2 t1| = W,
//
// and the straight side has 2m + W edges against the arc's 2S + W. Their
// totals have one parity -- the loop is a quad region's boundary, so it is
// even -- but on a disk the arc is about pi/2 times longer than the side, so
// with the core at Sec. 7.3's proportions |s1 s2| has fewer edges than
// |t2 t1|. Points are inserted where the difference falls (Sec. 8.4), which on
// the axis of an (r, z) model is dS. Where it cannot be, the only split left
// is the one whose core is as wide as the arc above it; that candidate is
// kept, ranked last.
bool ExplicitTemplates::tryHalfOGrid(const Region &R, Patch &P, Attempt &A) {
    const std::vector<int> &L = R.B.loops[R.outer];
    const int N = static_cast<int>(L.size());
    auto fail = [&](const std::string &why) { A.reason = why; return false; };

    std::vector<int> npos;
    for (int i = 0; i < N; ++i) if (isNode(L[i])) npos.push_back(i);
    if (npos.size() != 2) return fail("not two protected corners");
    const double maxCorner = opts_.halfOGridCorner * M_PI / 180.0;
    for (int i : npos) {
        if (R.B.angle[R.outer][i] > maxCorner) return fail("a corner is flat enough to want two cells");
    }

    // The two sides between the corners, in loop order. One must be straight.
    auto sideIds = [&](int from, int to) {
        std::vector<int> ids;
        for (int i = from;; i = (i + 1) % N) {
            ids.push_back(L[i]);
            if (i == to) break;
        }
        return ids;
    };
    auto positions = [&](const std::vector<int> &ids) {
        std::vector<Point> X;
        X.reserve(ids.size());
        for (int v : ids) X.push_back(C_.vertices[v]);
        return X;
    };
    const std::vector<int> side0 = sideIds(npos[0], npos[1]), side1 = sideIds(npos[1], npos[0]);
    const bool straight0 = collinear(positions(side0), 1e-6), straight1 = collinear(positions(side1), 1e-6);
    if (straight0 == straight1) {
        return fail(straight0 ? "both sides between the corners are straight"
                              : "neither side between the corners is straight");
    }
    const std::vector<int> &base = straight0 ? side0 : side1;   // A -> B
    const std::vector<int> &arc = straight0 ? side1 : side0;    // B -> A
    const int D = static_cast<int>(base.size()) - 1, Ar = static_cast<int>(arc.size()) - 1;
    if (D < 3 || Ar < 3) return fail("too few boundary edges for a core and three shells");
    const Point Apt = C_.vertices[base.front()], Bpt = C_.vertices[base.back()];
    const double len = normP(Bpt - Apt);
    if (!(len > 0.0)) return fail("the straight side has no length");
    const Point e = (Bpt - Apt) / len;
    const Point nrm{-e[1], e[0]};   // into the region
    const std::vector<Point> Xa = positions(arc);

    // c, on the straight side and in the kernel of the arc. The side is a ray
    // from any point of it, so only the arc's edges bound the depth, and their
    // minimum is concave along the side. Sec. 7.3's O-grid starts its kernel
    // point from the area centroid; the mirror image's centroid is this
    // region's dropped onto the side, and c stays there unless that is less
    // than half as deep as the deepest point, when it moves towards it until
    // it is. Not the deepest point itself: on a long flat body (rocket's
    // capsule) the depth is the height all along the middle, the maximum is a
    // plateau, and wherever a search stops on it the core is off-centre.
    auto depthAt = [&](double u) {
        const Point p = Apt + (Bpt - Apt) * u;
        double d = 1e300;
        for (int k = 0; k < Ar; ++k) {
            const Point a = Xa[k], b = Xa[k + 1];
            const double l = normP(b - a);
            if (l > 0.0) d = std::min(d, cross2(b - a, p - a) / l);
        }
        return d;
    };
    double lo = 0.0, hi = 1.0;
    for (int it = 0; it < 80; ++it) {
        const double u1 = lo + (hi - lo) / 3.0, u2 = hi - (hi - lo) / 3.0;
        if (depthAt(u1) < depthAt(u2)) lo = u1;
        else hi = u2;
    }
    const double uDeep = 0.5 * (lo + hi), dMax = depthAt(uDeep);
    if (!(dMax > 1e-3 * R.h)) return fail("not star-shaped about any point of the straight side");
    double uc = std::clamp(dotP(CavityFill::areaCentroid(positions(L)) - Apt, e) / len, 0.0, 1.0);
    if (depthAt(uc) < 0.5 * dMax) {
        // Concave, so the depth rises monotonically from uc to uDeep.
        double a = uc, b = uDeep;
        for (int it = 0; it < 60; ++it) {
            const double m = 0.5 * (a + b);
            if (depthAt(m) < 0.5 * dMax) a = m;
            else b = m;
        }
        uc = b;
    }
    const Point c = Apt + (Bpt - Apt) * uc;

    // Sec. 7.3's ideal split directions for the mirror image: the corners of
    // the rectangle on its axes -- the side and its normal -- with half-widths
    // the extents. On a half disk, 45 and 135 degrees.
    std::vector<double> th(Ar + 1);
    double w1p = 0.0, w1m = 0.0, w2 = 0.0;
    for (int k = 0; k <= Ar; ++k) {
        const Point d = Xa[k] - c;
        th[k] = std::atan2(dotP(d, nrm), dotP(d, e));
        w1p = std::max(w1p, dotP(d, e));
        w1m = std::max(w1m, -dotP(d, e));
        w2 = std::max(w2, dotP(d, nrm));
    }
    double ideal2 = std::atan2(w2, w1p), ideal1 = M_PI - std::atan2(w2, w1m);
    // A half disk's field has two +1/4 cones, the mirror image's four cut
    // back: when the region holds exactly those, the shells split on the rays
    // through them and the core's top corners go where they are.
    Point coneB{0.0, 0.0}, coneA{0.0, 0.0};
    bool fromCones = false;
    if (opts_.field && opts_.conePlacement) {
        std::vector<Point> plus;
        int minus = 0;
        for (const ReferenceField::Singularity &s : regionCones_) {
            if (s.sign > 0) plus.push_back(s.x);
            else ++minus;
        }
        if (plus.size() == 2 && minus == 0) {
            if (dotP(plus[0] - c, e) < dotP(plus[1] - c, e)) std::swap(plus[0], plus[1]);
            const Point dB = plus[0] - c, dA = plus[1] - c;
            if (dotP(dB, e) > 0.0 && dotP(dA, e) < 0.0 && dotP(dB, nrm) > 0.0 && dotP(dA, nrm) > 0.0) {
                coneB = plus[0];
                coneA = plus[1];
                ideal2 = std::atan2(dotP(dB, nrm), dotP(dB, e));
                ideal1 = std::atan2(dotP(dA, nrm), dotP(dA, e));
                fromCones = true;
            }
        }
    }
    A.conePlaced = fromCones;
    std::vector<double> ub(D + 1);
    for (int k = 0; k <= D; ++k) ub[k] = dotP(C_.vertices[base[k]] - c, e);

    // Which stretches may take inserted points: prefix counts of dS segments.
    auto splittablePrefix = [&](const std::vector<int> &ids) {
        std::vector<int> pre(ids.size(), 0);
        for (size_t k = 0; k + 1 < ids.size(); ++k) {
            const int ed = C_.edgeBetween(ids[k], ids[k + 1]);
            pre[k + 1] = pre[k] + ((opts_.splitBoundary && ed >= 0 && C_.boundaryEdge[ed]) ? 1 : 0);
        }
        return pre;
    };
    const std::vector<int> preB = splittablePrefix(base), preA = splittablePrefix(arc);
    auto baseSplits = [&](int i, int j) { return preB[j] - preB[i] > 0; };
    auto arcSplits = [&](int i, int j) { return preA[j] - preA[i] > 0; };

    struct Split { int j2, j1, i1, i2; double cost; };
    std::vector<Split> splits;
    auto nearestArc = [&](double target) {
        int best = 1;
        for (int k = 1; k < Ar; ++k) if (std::fabs(th[k] - target) < std::fabs(th[best] - target)) best = k;
        return best;
    };
    // How far along its ray each top corner sits: Sec. 7.3's fraction, or the
    // two cones' mean radius over the arc's there -- one fraction for both,
    // as the O-grid's, so the core stays the mirror image's.
    double alpha = opts_.coreFraction;
    if (fromCones) {
        const double RB = normP(Xa[nearestArc(ideal2)] - c), RA = normP(Xa[nearestArc(ideal1)] - c);
        if (RB > 0.0 && RA > 0.0) alpha = std::clamp(0.5 * (normP(coneB - c) / RB + normP(coneA - c) / RA), 0.3, 0.75);
    }
    auto alphaAt = [&](const Point &, bool) { return alpha; };
    const int j2c = nearestArc(ideal2), j1c = nearestArc(ideal1);
    const int window = 2;
    for (int j2 = std::max(1, j2c - window); j2 <= std::min(Ar - 2, j2c + window); ++j2) {
        std::set<int> j1s;
        for (int j1 = j1c - window; j1 <= j1c + window; ++j1) j1s.insert(j1);
        j1s.insert(Ar - j2);   // |B t2| = |t1 A| with nothing inserted on the arc
        for (int j1 : j1s) {
            if (j1 <= j2 || j1 > Ar - 1) continue;
            const Point q2 = c + (Xa[j2] - c) * alphaAt(Xa[j2], true), q1 = c + (Xa[j1] - c) * alphaAt(Xa[j1], false);
            const double x2 = dotP(q2 - c, e), x1 = dotP(q1 - c, e);
            if (!(x1 < 0.0 && x2 > 0.0)) continue;
            // The core's feet go under its top corners, as in the mirrored O-grid.
            int i1c = -1, i2c = -1;
            for (int k = 1; k < D; ++k) {
                if (ub[k] < 0.0 && (i1c < 0 || std::fabs(ub[k] - x1) < std::fabs(ub[i1c] - x1))) i1c = k;
                if (ub[k] > 0.0 && (i2c < 0 || std::fabs(ub[k] - x2) < std::fabs(ub[i2c] - x2))) i2c = k;
            }
            auto consider = [&](int i1, int i2) {
                if (i1 < 1 || i2 > D - 1 || i1 >= i2) return;
                if (!(ub[i1] < 0.0 && ub[i2] > 0.0)) return;
                const int nR = j2, nT = j1 - j2, nL = Ar - j1, mA = i1, w0 = i2 - i1, mB = D - i2;
                const int S = std::max(nR, nL), M = std::max(mA, mB), W = std::max(w0, nT);
                if ((S > nR && !arcSplits(0, j2)) || (S > nL && !arcSplits(j1, Ar))) return;
                if ((M > mA && !baseSplits(0, i1)) || (M > mB && !baseSplits(i2, D))) return;
                if ((W > w0 && !baseSplits(i1, i2)) || (W > nT && !arcSplits(j2, j1))) return;
                const double da = th[j2] - ideal2, db = th[j1] - ideal1;
                const double pa = (ub[i1] - x1) / len, pb = (ub[i2] - x2) / len;
                const int extra = (S - nR) + (S - nL) + (M - mA) + (M - mB) + (W - w0) + (W - nT);
                splits.push_back({j2, j1, i1, i2, da * da + db * db + pa * pa + pb * pb + 1e-9 * extra});
            };
            if (i1c >= 0 && i2c >= 0) {
                for (int d1 = -window; d1 <= window; ++d1)
                    for (int d2 = -window; d2 <= window; ++d2) consider(i1c + d1, i2c + d2);
            }
            // Nothing insertable anywhere: the core as wide as the arc above it.
            const int nT = j1 - j2;
            if (D - nT >= 2 && (D - nT) % 2 == 0) consider((D - nT) / 2, D - (D - nT) / 2);
        }
    }
    if (splits.empty()) return fail("no split of the two sides closes the counts");
    std::stable_sort(splits.begin(), splits.end(), [](const Split &a, const Split &b) { return a.cost < b.cost; });

    // Build one split; false with a reason when its geometry does not close.
    auto build = [&](const Split &sp, Patch &Q, int &inserted, std::string &why) {
        auto sub = [](const std::vector<int> &v, int a, int b) {
            return std::vector<int>(v.begin() + a, v.begin() + b + 1);
        };
        std::vector<int> bA = sub(base, 0, sp.i1), bW = sub(base, sp.i1, sp.i2), bB = sub(base, sp.i2, D);
        std::vector<int> aR = sub(arc, 0, sp.j2), aT = sub(arc, sp.j2, sp.j1), aL = sub(arc, sp.j1, Ar);
        const int S = std::max(sp.j2, Ar - sp.j1), M = std::max(sp.i1, D - sp.i2);
        const int W = std::max(sp.i2 - sp.i1, sp.j1 - sp.j2);
        inserted = 0;
        auto grow = [&](std::vector<int> &ids, int to) {
            const int extra = to - (static_cast<int>(ids.size()) - 1);
            if (extra <= 0) return true;
            ids = splitArc(Q, ids, extra);
            inserted += extra;
            return !ids.empty();
        };
        if (!grow(aR, S) || !grow(aL, S) || !grow(bA, M) || !grow(bB, M) || !grow(bW, W) || !grow(aT, W)) {
            why = "a stretch that needs points is not dS";
            return false;
        }
        auto X = [&](int id) { return Q.at(C_, id); };
        const int s1 = bA.back(), s2 = bB.front(), t2 = aR.back(), t1 = aL.front();
        const Point q2p = c + (X(t2) - c) * alphaAt(X(t2), true), q1p = c + (X(t1) - c) * alphaAt(X(t1), false);
        const std::array<Point, 4> quad = {X(s1), X(s2), q2p, q1p};
        for (int k = 0; k < 4; ++k) {
            if (!(cross2(quad[(k + 1) % 4] - quad[k], quad[(k + 3) % 4] - quad[k]) > 0.0)) {
                why = "the core quadrilateral is not convex";
                return false;
            }
        }
        const int q2 = Q.addVertex(q2p, Origin::Template), q1 = Q.addVertex(q1p, Origin::Template);

        // The core side a -> b's nodes: where the rays from c through the
        // shell's arc cross it, a and b themselves at the ends.
        auto raysOnto = [&](const std::vector<int> &outer, int a, int b, std::vector<int> &out) {
            out.assign(1, a);
            const Point pa = X(a), ed = X(b) - pa, w = pa - c;
            double lastS = 0.0;
            const int n = static_cast<int>(outer.size()) - 1;
            for (int k = 1; k < n; ++k) {
                const Point d = X(outer[k]) - c;
                const double den = cross2(d, ed);
                if (std::fabs(den) < 1e-300) return false;
                const double t = cross2(w, ed) / den, s = cross2(w, d) / den;
                if (!(t > 0.0 && t < 1.0 && s > lastS && s < 1.0)) return false;
                lastS = s;
                out.push_back(Q.addVertex(c + d * t, Origin::Template));
            }
            out.push_back(b);
            return true;
        };
        std::vector<int> coreR, coreT, coreL;
        if (!raysOnto(aR, s2, q2, coreR) || !raysOnto(aT, q2, q1, coreT) || !raysOnto(aL, q1, s1, coreL)) {
            why = "a ray misses its core side";
            return false;
        }

        // Radial lines, core end first: uniform on the two split rays, the
        // side's own points on the two stretches of the straight side.
        auto uniformRay = [&](int from, int to) {
            std::vector<int> r{from};
            for (int j = 1; j < M; ++j) {
                r.push_back(Q.addVertex(X(from) + (X(to) - X(from)) * (static_cast<double>(j) / M), Origin::Template));
            }
            r.push_back(to);
            return r;
        };
        const std::vector<int> rayT2 = uniformRay(q2, t2), rayT1 = uniformRay(q1, t1);
        const std::vector<int> rayB = bB, rayA(bA.rbegin(), bA.rend());

        // One shell: core[i] -> outer[i], i = 0..n, with `first` and `last` its
        // radial lines at i = 0 and n. A radial line in between takes the two
        // end lines' fractions blended by i/n, so it meets both exactly and its
        // fractions stay increasing -- all Sec. 7.3's det DX > 0 asks of lambda.
        auto shell = [&](const std::vector<int> &core, const std::vector<int> &outer, const std::vector<int> &first,
                         const std::vector<int> &last) {
            const int n = static_cast<int>(outer.size()) - 1;
            auto fractions = [&](const std::vector<int> &ray) {
                std::vector<double> f(M + 1);
                const Point a = X(ray.front());
                const double l = normP(X(ray.back()) - a);
                for (int j = 0; j <= M; ++j) f[j] = normP(X(ray[j]) - a) / l;
                return f;
            };
            const std::vector<double> f0 = fractions(first), f1 = fractions(last);
            std::vector<int> nodes(static_cast<size_t>(n + 1) * (M + 1));
            auto at = [&](int i, int j) -> int & { return nodes[i + (n + 1) * j]; };
            for (int i = 0; i <= n; ++i) { at(i, 0) = core[i]; at(i, M) = outer[i]; }
            for (int j = 0; j <= M; ++j) { at(0, j) = first[j]; at(n, j) = last[j]; }
            for (int i = 1; i < n; ++i) {
                const double lam = static_cast<double>(i) / n;
                const Point a0 = X(core[i]), a1 = X(outer[i]);
                for (int j = 1; j < M; ++j) {
                    at(i, j) = Q.addVertex(a0 + (a1 - a0) * ((1.0 - lam) * f0[j] + lam * f1[j]), Origin::Template);
                }
            }
            for (int j = 0; j < M; ++j) {
                for (int i = 0; i < n; ++i) Q.cells.push_back({at(i, j), at(i, j + 1), at(i + 1, j + 1), at(i + 1, j)});
            }
        };
        shell(coreR, aR, rayB, rayT2);
        shell(coreT, aT, rayT2, rayT1);
        shell(coreL, aL, rayT1, rayA);

        // The core: a straight-line grid, as the O-grid's, its bottom the
        // straight side between the feet.
        {
            const int n = W, mm = S;
            std::vector<int> nodes(static_cast<size_t>(n + 1) * (mm + 1));
            auto at = [&](int i, int j) -> int & { return nodes[i + (n + 1) * j]; };
            for (int i = 0; i <= n; ++i) { at(i, 0) = bW[i]; at(i, mm) = coreT[n - i]; }
            for (int j = 0; j <= mm; ++j) { at(n, j) = coreR[j]; at(0, j) = coreL[mm - j]; }
            for (int j = 1; j < mm; ++j) {
                for (int i = 1; i < n; ++i) {
                    const Point b0 = X(at(i, 0)), t0 = X(at(i, mm)), l0 = X(at(0, j)), r0 = X(at(n, j));
                    const Point r = t0 - b0, s = r0 - l0;
                    const double den = cross2(r, s);
                    if (std::fabs(den) < 1e-300) {
                        why = "two core chords are parallel";
                        return false;
                    }
                    at(i, j) = Q.addVertex(b0 + r * (cross2(l0 - b0, s) / den), Origin::Template);
                }
            }
            for (int j = 0; j < mm; ++j) {
                for (int i = 0; i < n; ++i) Q.cells.push_back({at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)});
            }
        }
        for (int v : {s1, s2, t1, t2, q1, q2}) Q.designated.push_back(v);
        Q.blocks = 4;
        Q.kind = "half O-grid";
        return true;
    };

    // The best few splits by Sec. 7.3's proportions; the first that certifies.
    std::string why = "no split built";
    const size_t tries = std::min<size_t>(splits.size(), 6);
    for (size_t k = 0; k < tries; ++k) {
        Patch Q;
        int inserted = 0;
        std::string reason;
        if (!build(splits[k], Q, inserted, reason)) {
            why = reason;
            continue;
        }
        const CavityFill::Verdict V = CavityFill::validate(C_, R.cells, Q, opts_.minScaledJacobian);
        if (!V.valid) {
            why = V.reason;
            continue;
        }
        P = std::move(Q);
        A.splits += inserted;
        A.maps = "1 straight-line core, 3 radial";
        return true;
    }
    return fail(why);
}

// ---------------------------------------------------------------------------
// Annulus: Sec. 7.4.
// ---------------------------------------------------------------------------
bool ExplicitTemplates::tryAnnulus(const Region &R, Patch &P, Attempt &A) {
    auto fail = [&](const std::string &why) { A.reason = why; return false; };
    if (R.holes.size() != 1) return fail("not exactly one hole");
    const std::vector<int> &Lo = R.B.loops[R.outer];
    // The hole loop runs clockwise round the hole (region on its left);
    // reversed it runs counter-clockwise round c, like the outer loop.
    std::vector<int> Li = R.B.loops[R.holes[0]];
    std::reverse(Li.begin(), Li.end());
    std::vector<int> nodesO, nodesI;
    for (size_t i = 0; i < Lo.size(); ++i) if (isNode(Lo[i])) nodesO.push_back(static_cast<int>(i));
    for (size_t i = 0; i < Li.size(); ++i) if (isNode(Li[i])) nodesI.push_back(static_cast<int>(i));
    if (!(nodesO.empty() || nodesO.size() == 4) || !(nodesI.empty() || nodesI.size() == 4)) {
        return fail("protected corners on a loop, and not four of them");
    }

    const int No = static_cast<int>(Lo.size()), Ni = static_cast<int>(Li.size());
    std::vector<Point> Xo(No), Xi(Ni);
    for (int i = 0; i < No; ++i) Xo[i] = C_.vertices[Lo[i]];
    for (int i = 0; i < Ni; ++i) Xi[i] = C_.vertices[Li[i]];
    const Point c = CavityFill::kernelCenter(Xi, CavityFill::areaCentroid(Xi));
    if (!(CavityFill::kernelDepth(Xi, c) > 1e-3 * R.h)) return fail("the hole is not a radial graph about any point found");
    if (!(CavityFill::kernelDepth(Xo, c) > 1e-3 * R.h)) return fail("the outer loop is not a radial graph about the hole's centre");

    const double phi = principalAngle(Xo) + M_PI_4;
    auto nearest = [&](const std::vector<Point> &X, double ang) {
        int best = 0;
        double bd = 1e300;
        for (size_t i = 0; i < X.size(); ++i) {
            double d = std::fabs(wrap2pi(std::atan2(X[i][1] - c[1], X[i][0] - c[0]) - ang));
            d = std::min(d, 2.0 * M_PI - d);
            if (d < bd) { bd = d; best = static_cast<int>(i); }
        }
        return best;
    };
    // The four rays: through the protected corners when a loop has four,
    // else on the diagonals of the principal frame.
    std::array<double, 4> rayAngle{};
    for (int k = 0; k < 4; ++k) rayAngle[k] = phi + k * M_PI_2;
    if (nodesO.size() == 4) {
        for (int k = 0; k < 4; ++k) rayAngle[k] = std::atan2(Xo[nodesO[k]][1] - c[1], Xo[nodesO[k]][0] - c[0]);
    } else if (nodesI.size() == 4) {
        for (int k = 0; k < 4; ++k) rayAngle[k] = std::atan2(Xi[nodesI[k]][1] - c[1], Xi[nodesI[k]][0] - c[0]);
    }
    std::array<int, 4> si, so;
    for (int k = 0; k < 4; ++k) {
        si[k] = nodesI.size() == 4 ? nodesI[k] : nearest(Xi, rayAngle[k]);
        so[k] = nodesO.size() == 4 ? nodesO[k] : nearest(Xo, rayAngle[k]);
    }
    // Pair the two loops' split points by angle: rotate the inner list so
    // its first point is the one nearest the outer's first ray.
    {
        int best = 0;
        double bd = 1e300;
        for (int k = 0; k < 4; ++k) {
            double d = std::fabs(wrap2pi(std::atan2(Xi[si[k]][1] - c[1], Xi[si[k]][0] - c[0]) - rayAngle[0]));
            d = std::min(d, 2.0 * M_PI - d);
            if (d < bd) { bd = d; best = k; }
        }
        std::rotate(si.begin(), si.begin() + best, si.end());
    }
    auto ordered = [](const std::array<int, 4> &s, int n) {
        int total = 0;
        for (int k = 0; k < 4; ++k) {
            const int d = ((s[(k + 1) % 4] - s[k]) % n + n) % n;
            if (d == 0) return false;
            total += d;
        }
        return total == n;
    };
    if (!ordered(si, Ni) || !ordered(so, No)) return fail("the split rays do not meet the loops in order");

    auto arcOf = [](const std::vector<int> &L, int a, int b) {
        const int n = static_cast<int>(L.size());
        std::vector<int> ids{L[a]};
        int i = a;
        while (i != b) { i = (i + 1) % n; ids.push_back(L[i]); }
        return ids;
    };
    auto allBoundary = [&](const std::vector<int> &ids) {
        if (!opts_.splitBoundary) return false;
        for (size_t j = 0; j + 1 < ids.size(); ++j) {
            const int e = C_.edgeBetween(ids[j], ids[j + 1]);
            if (e < 0 || !C_.boundaryEdge[e]) return false;
        }
        return true;
    };
    std::array<std::vector<int>, 4> ai, ao;
    for (int k = 0; k < 4; ++k) {
        ai[k] = arcOf(Li, si[k], si[(k + 1) % 4]);
        ao[k] = arcOf(Lo, so[k], so[(k + 1) % 4]);
    }
    // Equalise each sector's two arcs by subdividing the coarser, if it is
    // domain boundary. Failing that, and only when the loops have equal
    // counts, re-place the outer split points by count instead.
    bool needFixed = false;
    for (int k = 0; k < 4; ++k) {
        const int ni = static_cast<int>(ai[k].size()) - 1, no = static_cast<int>(ao[k].size()) - 1;
        if (ni < no && !allBoundary(ai[k])) needFixed = true;
        if (no < ni && !allBoundary(ao[k])) needFixed = true;
    }
    if (needFixed) {
        if (Ni != No) return fail("sector counts differ and the coarser arc is an interface");
        if (!nodesO.empty()) return fail("sector counts differ, the coarser arc is an interface, and the corners pin the outer split points");
        int pos = so[0];
        for (int k = 0; k < 4; ++k) {
            const int ni = static_cast<int>(ai[k].size()) - 1;
            ao[k] = arcOf(Lo, pos, (pos + ni) % No);
            pos = (pos + ni) % No;
        }
    } else {
        for (int k = 0; k < 4; ++k) {
            const int ni = static_cast<int>(ai[k].size()) - 1, no = static_cast<int>(ao[k].size()) - 1;
            if (ni < no) { ai[k] = splitArc(P, ai[k], no - ni); A.splits += no - ni; }
            if (no < ni) { ao[k] = splitArc(P, ao[k], ni - no); A.splits += ni - no; }
            if (ai[k].empty() || ao[k].empty()) return fail("a sector arc could not be subdivided");
        }
    }

    double radial = 0.0;
    int rc = 0;
    for (int k = 0; k < 4; ++k) {
        for (size_t j = 0; j < ai[k].size(); ++j) {
            radial += normP(P.at(C_, ao[k][j]) - P.at(C_, ai[k][j]));
            ++rc;
        }
    }
    const int m = std::max(1, static_cast<int>(std::lround((radial / std::max(1, rc)) / std::max(1e-300, R.h))));

    // Polar interpolation about c between matching inner and outer nodes.
    auto polar = [&](const Point &u, const Point &w, double t) {
        const double tu = std::atan2(u[1] - c[1], u[0] - c[0]);
        double tw = std::atan2(w[1] - c[1], w[0] - c[0]);
        while (tw - tu > M_PI) tw -= 2.0 * M_PI;
        while (tw - tu < -M_PI) tw += 2.0 * M_PI;
        const double th = (1.0 - t) * tu + t * tw;
        const double r = (1.0 - t) * normP(u - c) + t * normP(w - c);
        return c + Point{std::cos(th), std::sin(th)} * r;
    };
    std::array<std::vector<int>, 4> ray;
    for (int k = 0; k < 4; ++k) {
        const int u = ai[k].front(), w = ao[k].front();
        ray[k].push_back(u);
        for (int j = 1; j < m; ++j) ray[k].push_back(P.addVertex(polar(P.at(C_, u), P.at(C_, w), static_cast<double>(j) / m), Origin::Template));
        ray[k].push_back(w);
    }
    for (int k = 0; k < 4; ++k) {
        const int n = static_cast<int>(ai[k].size()) - 1;
        std::vector<int> nodes(static_cast<size_t>(n + 1) * (m + 1));
        auto at = [&](int i, int j) -> int & { return nodes[i + (n + 1) * j]; };
        for (int i = 0; i <= n; ++i) { at(i, 0) = ai[k][i]; at(i, m) = ao[k][i]; }
        for (int j = 0; j <= m; ++j) { at(0, j) = ray[k][j]; at(n, j) = ray[(k + 1) % 4][j]; }
        for (int i = 1; i < n; ++i) {
            const Point u = P.at(C_, ai[k][i]), w = P.at(C_, ao[k][i]);
            for (int j = 1; j < m; ++j) at(i, j) = P.addVertex(polar(u, w, static_cast<double>(j) / m), Origin::Template);
        }
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < n; ++i) P.cells.push_back({at(i, j), at(i, j + 1), at(i + 1, j + 1), at(i + 1, j)});
        }
    }
    for (int k = 0; k < 4; ++k) { P.designated.push_back(ai[k].front()); P.designated.push_back(ao[k].front()); }
    P.blocks = 4;
    P.kind = "annulus";
    A.maps = "4 polar";
    return true;
}
