#include "QUARRY/QUARRY.hxx"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <tuple>

namespace {

// Distance from a point to dS, for points that need it only up to a bound: the boundary edges in a uniform grid of
// cells about four edges wide, searched ring by ring out to the bound.
class BoundaryGrid {
public:
    explicit BoundaryGrid(const Mesh &M) : M_(M) {
        double total = 0.0;
        for (const int e : M.boundaryEdges) total += normP(M.vertices[M.edges[e][1]] - M.vertices[M.edges[e][0]]);
        if (M.boundaryEdges.empty()) return;
        cell_ = std::max(4.0 * total / static_cast<double>(M.boundaryEdges.size()), 1e-12);
        for (const int e : M.boundaryEdges) {
            const Point &a = M.vertices[M.edges[e][0]], &b = M.vertices[M.edges[e][1]];
            const long x0 = key(std::min(a[0], b[0])), x1 = key(std::max(a[0], b[0]));
            const long y0 = key(std::min(a[1], b[1])), y1 = key(std::max(a[1], b[1]));
            for (long x = x0; x <= x1; ++x)
                for (long y = y0; y <= y1; ++y) cells_[{x, y}].push_back(e);
        }
    }
    // The distance to dS, or `bound` if dS is further than that.
    double distance(const Point &p, double bound) const {
        if (cells_.empty()) return bound;
        const long cx = key(p[0]), cy = key(p[1]);
        const long rings = static_cast<long>(std::ceil(bound / cell_)) + 1;
        double best = bound;
        for (long r = 0; r <= rings && (r - 1) * cell_ < best; ++r)
            for (long x = cx - r; x <= cx + r; ++x)
                for (long y = cy - r; y <= cy + r; ++y) {
                    if (std::max(std::labs(x - cx), std::labs(y - cy)) != r) continue;
                    const auto it = cells_.find({x, y});
                    if (it == cells_.end()) continue;
                    for (const int e : it->second) {
                        const Point &a = M_.vertices[M_.edges[e][0]], &b = M_.vertices[M_.edges[e][1]];
                        const Point ab = b - a;
                        const double t = std::clamp(dotP(p - a, ab) / std::max(dotP(ab, ab), 1e-300), 0.0, 1.0);
                        best = std::min(best, normP(p - (a + ab * t)));
                    }
                }
        return best;
    }

private:
    long key(double v) const { return static_cast<long>(std::floor(v / cell_)); }
    const Mesh &M_;
    double cell_ = 1.0;
    std::map<std::pair<long, long>, std::vector<int>> cells_;
};

}  // namespace

QUARRY::QUARRY(std::shared_ptr<const Mesh> mesh, Fields fields, const Options &opts)
    : mesh_(std::move(mesh)), fields_(std::move(fields)), opts_(opts) {
    const size_t nt = mesh_->triangles.size();
    if (fields_.cdf.size() != nt || fields_.dcdf.size() != nt || fields_.offsetCentre.size() != nt ||
        fields_.offsetCorner.size() != nt)
        throw std::invalid_argument("QUARRY: cdf, dcdf, offsetc and offset need one value per triangle");
}

// Stages 1-2, then Stages 3-4 and 6 as often as Stage 5 changes the charts. A chart Stage 4 cannot give four
// corners -- fewer than four junctions and nothing on dS to take the rest from -- is no chart of a quad layout
// however it came about: a sliver the decoder left in a wall, a region a wall gap shrank to a corner, an island.
// It goes to the neighbour it shares the longest boundary with, and the outlines are read again, since its
// neighbours' junctions with it have gone with it. A chart with corners to spare is never dissolved: giving it
// away would only hand its excess to a neighbour.
void QUARRY::run() {
    const auto t0 = std::chrono::steady_clock::now();
    const Mesh &M = *mesh_;
    const int nt = static_cast<int>(M.triangles.size());
    centreTargets_.resize(nt);
    cornerTargets_.resize(nt);
    for (int t = 0; t < nt; ++t) {
        const Point c = (M.vertices[M.triangles[t][0]] + M.vertices[M.triangles[t][1]] +
                         M.vertices[M.triangles[t][2]]) / 3.0;
        centreTargets_[t] = c + fields_.offsetCentre[t];
        cornerTargets_[t] = c + fields_.offsetCorner[t];
    }
    charts_ = std::make_unique<quarry::Basins>(M, fields_.cdf, centreTargets_, opts_.charts);
    duals_ = std::make_unique<quarry::Basins>(M, fields_.dcdf, cornerTargets_, opts_.duals);
    report_.charts = charts_->getReport();
    report_.duals = duals_->getReport();
    features_ = std::make_unique<BoundaryFeatures>(M);
    labels_ = charts_->labels();
    status_.assign(charts_->regions().size(), ChartStatus());

    for (int round = 0;; ++round) {
        std::vector<int> drop = assemble();
        if (drop.empty() || round >= opts_.repairRounds) break;
        ++report_.repairRounds;
        std::vector<double> area(status_.size(), 0.0);
        for (int t = 0; t < nt; ++t) {
            if (labels_[t] < 0) continue;
            const Point &a = M.vertices[M.triangles[t][0]];
            area[labels_[t]] += 0.5 * std::fabs(cross2(M.vertices[M.triangles[t][1]] - a,
                                                      M.vertices[M.triangles[t][2]] - a));
        }
        // The smallest first, and one per round: giving one chart away changes its neighbours' outlines, and a
        // neighbour that was short of a corner may not be any more. A chart with no neighbour to go to (the whole
        // domain one chart short of corners) is passed over.
        std::sort(drop.begin(), drop.end(), [&](int a, int b) { return area[a] < area[b]; });
        bool dissolved = false;
        for (const int k : drop) {
            std::map<int, double> shared;
            for (int t = 0; t < nt; ++t) {
                if (labels_[t] != k) continue;
                for (int e = 0; e < 3; ++e) {
                    const int u = M.triangleAdjacency[t][e];
                    if (u < 0 || labels_[u] < 0 || labels_[u] == k) continue;
                    shared[labels_[u]] +=
                        normP(M.vertices[M.triangles[t][(e + 1) % 3]] - M.vertices[M.triangles[t][e]]);
                }
            }
            int best = -1;
            for (const auto &[nb, len] : shared)
                if (best < 0 || len > shared[best]) best = nb;
            if (best < 0) continue;
            for (int t = 0; t < nt; ++t)
                if (labels_[t] == k) labels_[t] = best;
            status_[k].dissolved = true;
            ++report_.dissolvedCharts;
            dissolved = true;
            break;
        }
        if (!dissolved) break;
    }
    report_.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

int QUARRY::group(int v) {
    while (groupParent_[v] != v) {
        groupParent_[v] = groupParent_[groupParent_[v]];
        v = groupParent_[v];
    }
    return v;
}

// Every chart's boundary as closed curves of half-edges of its own triangles, the chart on their left. From a
// half-edge a -> b the next one is found by turning clockwise round b through the chart's own triangles until an
// edge with something else across it: that keeps the chart on the left through a vertex where it touches itself,
// which a walk by vertex alone would get wrong. The turn always ends, since the triangle across a -> b is not the
// chart's.
std::vector<std::vector<QUARRY::Cycle>> QUARRY::boundaryCycles() const {
    const Mesh &M = *mesh_;
    const std::vector<int> &L = labels_;
    const int nt = static_cast<int>(M.triangles.size());
    std::vector<std::vector<Cycle>> out(status_.size());
    std::vector<char> seen(3 * static_cast<size_t>(nt), 0);
    auto outside = [&](int t, int e) {
        const int u = M.triangleAdjacency[t][e];
        return u < 0 || L[u] != L[t];
    };
    for (int t0 = 0; t0 < nt; ++t0) {
        if (L[t0] < 0) continue;
        for (int e0 = 0; e0 < 3; ++e0) {
            if (seen[3 * t0 + e0] || !outside(t0, e0)) continue;
            Cycle c;
            int t = t0, e = e0;
            bool closed = true;
            do {
                seen[3 * t + e] = 1;
                c.vertex.push_back(M.triangles[t][e]);
                const int u = M.triangleAdjacency[t][e];
                c.across.push_back(u < 0 ? -1 : (L[u] < 0 ? -2 : L[u]));
                c.triangle.push_back(t);
                const int b = M.triangles[t][(e + 1) % 3];
                int s = t, f = (e + 1) % 3;
                int guard = 0;
                while (!outside(s, f) && ++guard < 4096) {
                    const int v = M.triangleAdjacency[s][f];
                    int i = 0;
                    while (i < 3 && M.triangles[v][i] != b) ++i;
                    s = v;
                    f = i;
                }
                if (guard >= 4096) {
                    closed = false;
                    break;
                }
                t = s;
                e = f;
            } while (!(t == t0 && e == e0));
            if (!closed) continue;
            const size_t n = c.vertex.size();
            for (size_t i = 0; i < n; ++i)
                c.signedArea += 0.5 * cross2(M.vertices[c.vertex[i]], M.vertices[c.vertex[(i + 1) % n]]);
            out[L[t0]].push_back(std::move(c));
        }
    }
    return out;
}

// Stages 3, 4 and 6.
//
// Stage 3. A junction is a mesh vertex touched by three or more of {the charts, dS}; everything else on a chart's
// outline lies between exactly two of them. So a chart's outline -- one closed curve if the chart is a disk -- falls
// into arcs at its junctions, each arc against one neighbour or against dS, and an arc between two charts is met
// from both sides: once in each outline, the other way round (its twin).
//
// Stage 4. A conforming quad layout has no T-junction, so every junction on a chart's outline is one of its corners,
// and a chart of the true layout has four of them unless one of two things happened. (1) A vertex where four or
// more charts meet is rarely a single mesh vertex of the extraction mesh: the regions meet at two junctions, three
// charts each, joined by a short arc between the two charts that each see both -- and those two charts then count
// five corners. The arc is collapsed, its two junctions made one corner, when it is short and every chart that
// would lose a corner by it has one to spare; the second condition is what tells a blurred vertex from a short side,
// which no length test can. (2) A corner on dS that no other chart meets -- a corner of the model in one chart, or
// a chart's corner on a smooth stretch of dS -- is no junction at all. A chart short of corners takes them on its
// arcs along dS: first where a dual chart that covers part of it predicts its vertex, since the DCDF and its offset
// describe exactly the layout's vertices, then at the model's corners, sharpest first.
//
// Stage 6. A side runs between two consecutive corners. A side on dS is the boundary loop between them, walked
// with the domain on the left; a side between two charts is their common arc, its ends put at the corners' points
// (a dS corner at a mesh vertex of dS, mesh vertices of dS next to it dropped from the arc, so that no side runs
// along dS for a stretch), then smoothed. Both charts read the one polyline.
std::vector<int> QUARRY::assemble() {
    const Mesh &M = *mesh_;
    const int nt = static_cast<int>(M.triangles.size());
    const int nv = static_cast<int>(M.vertices.size());
    const int K = static_cast<int>(status_.size());
    for (ChartStatus &st : status_) {
        const bool d = st.dissolved;
        st = ChartStatus();
        st.dissolved = d;
    }
    blocks_ = BlockDecomposition();
    blocks_.source = "QUARRY";
    report_.junctions = report_.arcs = report_.collapsedArcs = report_.addedCorners = report_.snappedCorners = 0;
    report_.blocks = report_.failedCharts = report_.unmatchedSides = 0;
    report_.macroVertices = report_.macroEdges = report_.boundaryVertices = 0;
    report_.messages.clear();

    std::vector<double> triArea(nt), area(K, 0.0);
    for (int t = 0; t < nt; ++t) {
        const Point &a = M.vertices[M.triangles[t][0]];
        triArea[t] = 0.5 * std::fabs(cross2(M.vertices[M.triangles[t][1]] - a, M.vertices[M.triangles[t][2]] - a));
        if (labels_[t] >= 0) area[labels_[t]] += triArea[t];
    }

    // ---- Stage 3 ----
    std::vector<std::array<int, 3>> seen(nv, {{INT_MIN, INT_MIN, INT_MIN}});
    std::vector<unsigned char> distinct(nv, 0);
    auto note = [&](int v, int lab) {
        unsigned char &n = distinct[v];
        for (int i = 0; i < n; ++i)
            if (seen[v][i] == lab) return;
        if (n < 3) seen[v][n++] = lab;
    };
    for (int t = 0; t < nt; ++t)
        for (int i = 0; i < 3; ++i) note(M.triangles[t][i], labels_[t] >= 0 ? labels_[t] : -2);
    for (const int v : M.boundaryVertices) note(v, -1);
    std::vector<char> junction(nv, 0);
    for (int v = 0; v < nv; ++v) {
        junction[v] = distinct[v] >= 3;
        report_.junctions += junction[v];
    }

    const std::vector<std::vector<Cycle>> cycles = boundaryCycles();
    std::vector<Arc> arcs;
    std::vector<std::vector<int>> chartArcs(K);
    std::vector<char> outlineOk(K, 0);
    auto edgeLength = [&](const Cycle &c, int i) {
        const int n = static_cast<int>(c.vertex.size());
        return normP(M.vertices[c.vertex[(i + 1) % n]] - M.vertices[c.vertex[i]]);
    };
    for (int k = 0; k < K; ++k) {
        ChartStatus &st = status_[k];
        st.cycles = static_cast<int>(cycles[k].size());
        if (st.dissolved) {
            st.problem = "dissolved";
            continue;
        }
        if (cycles[k].size() != 1) {
            st.problem = cycles[k].empty() ? "no triangles" : std::to_string(cycles[k].size()) + " boundary curves";
            continue;
        }
        const Cycle &c = cycles[k][0];
        if (c.signedArea <= 0.0) {
            st.problem = "outline turns clockwise";
            continue;
        }
        const int n = static_cast<int>(c.vertex.size());
        std::vector<int> J;
        for (int i = 0; i < n; ++i)
            if (junction[c.vertex[i]]) J.push_back(i);
        st.junctions = static_cast<int>(J.size());
        if (J.empty()) {
            Arc a;
            a.chart = k;
            a.across = c.across[0];
            for (int i = 0; i < n; ++i) a.length += edgeLength(c, i);
            a.edges = n;
            chartArcs[k].push_back(static_cast<int>(arcs.size()));
            arcs.push_back(a);
        } else {
            for (size_t p = 0; p < J.size(); ++p) {
                Arc a;
                a.chart = k;
                a.begin = J[p];
                a.end = J[(p + 1) % J.size()];
                a.u = c.vertex[a.begin];
                a.v = c.vertex[a.end];
                a.across = c.across[a.begin];
                int i = a.begin;
                do {
                    a.length += edgeLength(c, i);
                    ++a.edges;
                    i = (i + 1) % n;
                } while (i != a.end);
                chartArcs[k].push_back(static_cast<int>(arcs.size()));
                arcs.push_back(a);
            }
        }
        outlineOk[k] = 1;
    }
    report_.arcs = static_cast<int>(arcs.size());
    {
        std::map<std::tuple<int, int, int, int>, std::vector<int>> byKey;
        for (int i = 0; i < static_cast<int>(arcs.size()); ++i)
            if (arcs[i].across >= 0) byKey[{arcs[i].u, arcs[i].v, arcs[i].chart, arcs[i].across}].push_back(i);
        for (int i = 0; i < static_cast<int>(arcs.size()); ++i) {
            Arc &a = arcs[i];
            if (a.across < 0 || a.twin >= 0) continue;
            const auto it = byKey.find({a.v, a.u, a.across, a.chart});
            if (it == byKey.end()) continue;
            for (const int j : it->second)
                if (arcs[j].twin < 0 && j != i) {
                    a.twin = j;
                    arcs[j].twin = i;
                    break;
                }
        }
    }

    // ---- Stage 4 ----
    groupParent_.resize(nv);
    std::iota(groupParent_.begin(), groupParent_.end(), 0);
    std::vector<std::vector<int>> chartsAt(nv);
    for (const Arc &a : arcs)
        if (a.u >= 0) chartsAt[a.u].push_back(a.chart);
    std::map<int, std::vector<int>> members;   // group -> its junctions
    for (int v = 0; v < nv; ++v)
        if (!chartsAt[v].empty()) members[v] = {v};
    // The corner groups along a chart's outline, a run of one group counted once; `from` read as `to`.
    auto runsOf = [&](int k, int from = -1, int to = -1) {
        std::vector<int> g;
        for (const int a : chartArcs[k]) {
            if (arcs[a].u < 0) continue;
            int x = group(arcs[a].u);
            if (x == from) x = to;
            if (g.empty() || g.back() != x) g.push_back(x);
        }
        while (g.size() > 1 && g.front() == g.back()) g.pop_back();
        return g;
    };
    auto pinched = [](const std::vector<int> &g) { return std::set<int>(g.begin(), g.end()).size() != g.size(); };
    std::vector<int> count(K, 0);
    for (int k = 0; k < K; ++k)
        if (outlineOk[k]) count[k] = static_cast<int>(runsOf(k).size());

    // The corners a chart could still take on dS, best first: where a dual chart that covers part of it predicts
    // a vertex near its stretch of dS (scored 1 + the dual's share of the chart), then the model's corners (scored
    // by their turning over pi). Neither moves when arcs are collapsed, so they are found once.
    std::vector<std::map<int, double>> overlap(K);
    const std::vector<int> &D = duals_->labels();
    for (int t = 0; t < nt; ++t)
        if (labels_[t] >= 0 && D[t] >= 0) overlap[labels_[t]][D[t]] += triArea[t];
    std::vector<double> turnAt(nv, 0.0);
    for (const BoundaryFeatures::Corner &c : features_->corners())
        turnAt[c.vertex] = std::fabs(M_PI - c.angle) / M_PI;
    std::vector<std::vector<std::pair<double, int>>> picks(K);
    for (int k = 0; k < K; ++k) {
        if (!outlineOk[k]) continue;
        const Cycle &c = cycles[k][0];
        const int n = static_cast<int>(c.vertex.size());
        const double scale = std::sqrt(area[k]);
        std::vector<int> onDS;   // cycle indices strictly inside k's arcs along dS
        for (const int a : chartArcs[k]) {
            const Arc &arc = arcs[a];
            if (arc.across != -1) continue;
            if (arc.u < 0) {
                for (int i = 0; i < n; ++i) onDS.push_back(i);
                continue;
            }
            for (int i = (arc.begin + 1) % n; i != arc.end; i = (i + 1) % n) onDS.push_back(i);
        }
        if (onDS.empty()) continue;
        // A vertex a dual chart predicts is one the chart has not got yet only if the nearest point of its dS is
        // nearer to it than every junction of the chart is, by more than cornerGap mean outline edges: otherwise a
        // junction already explains it, as one does wherever the decoder put a vertex a little off dS or two
        // vertices' dual charts merged.
        double perimeter = 0.0;
        for (int i = 0; i < n; ++i) perimeter += edgeLength(c, i);
        const double slack = opts_.cornerGap * perimeter / n;
        for (const auto &[j, share] : overlap[k]) {
            if (share < opts_.cornerShare * area[k]) continue;
            const Point &target = duals_->regions()[j].centre;
            double best = std::numeric_limits<double>::max();
            int at = -1;
            for (const int i : onDS) {
                const double d = normP(M.vertices[c.vertex[i]] - target);
                if (d < best) best = d, at = i;
            }
            bool explained = false;
            for (const int a : chartArcs[k])
                explained |= arcs[a].u >= 0 && normP(M.vertices[arcs[a].u] - target) <= best + slack;
            if (!explained && best <= opts_.dualReach * scale) picks[k].emplace_back(1.0 + share / area[k], at);
        }
        for (const int i : onDS)
            if (turnAt[c.vertex[i]] > 0.0) picks[k].emplace_back(turnAt[c.vertex[i]], i);
        // Last, the DCDF's own maxima along the stretch, scored below every model corner (turning 20 degrees or
        // more scores over 0.11): along dS the DCDF climbs to 1 at each of the chart's corners and falls to 0
        // between them, so a peak there is a corner even where its dual chart merged into a neighbour's and no
        // corner of the model marks the spot -- the case that otherwise costs a real chart in Stage 5.
        std::vector<char> inside(n, 0);
        for (const int i : onDS) inside[i] = 1;
        auto dcdfAt = [&](int i) { return fields_.dcdf[c.triangle[i]]; };
        for (const int i : onDS) {
            bool peak = dcdfAt(i) >= opts_.dcdfCorner;
            for (int d = 1; d <= opts_.cornerGap && peak; ++d)
                for (const int j : {(i + d) % n, (i + n - d) % n})
                    peak &= !inside[j] || dcdfAt(j) <= dcdfAt(i);
            if (peak) picks[k].emplace_back(0.1 * dcdfAt(i), i);
        }
        std::sort(picks[k].begin(), picks[k].end(), std::greater<>());
    }
    // Up to `need` of them, each clear of the chart's junctions and of each other: cornerGap mesh edges away along
    // the outline, or cornerSpacing x sqrt(area) away in the plane. Either says the candidate is a vertex of its
    // own -- a short side is a side, but a corner on top of another is none -- and neither alone is enough: a long
    // thin chart has its short side far below sqrt(area), and a chart a few triangles across has every vertex of
    // its outline within a few edges of a corner. `dualOnly` keeps to the corners a dual chart predicts (score
    // above 1): the model's own corners are everywhere on a polygonal dS -- a small hole drawn with a few dozen
    // chords turns by more than 20 degrees at every one -- so they settle what a chart takes at the end but say
    // nothing about whether a chart could do without a junction.
    auto choose = [&](int k, int need, bool dualOnly) {
        std::vector<int> out;
        if (need <= 0) return out;
        const Cycle &c = cycles[k][0];
        const int n = static_cast<int>(c.vertex.size());
        const double spacing = opts_.cornerSpacing * std::sqrt(area[k]);
        std::vector<int> taken;
        for (const int a : chartArcs[k])
            if (arcs[a].u >= 0) taken.push_back(arcs[a].begin);
        for (const auto &[score, i] : picks[k]) {
            if (static_cast<int>(out.size()) >= need || (dualOnly && score < 1.0)) break;
            bool clear = true;
            for (const int j : taken) {
                const int d = std::min(std::abs(i - j), n - std::abs(i - j));
                clear &= d > 0 && (d >= opts_.cornerGap ||
                                   normP(M.vertices[c.vertex[i]] - M.vertices[c.vertex[j]]) >= spacing);
            }
            if (!clear) continue;
            taken.push_back(i);
            out.push_back(i);
        }
        return out;
    };
    // How far a chart with `corners` junction corners is from being a quad, once it has taken what it can on dS.
    auto badness = [&](int k, int corners) {
        if (corners >= 4) return corners - 4;
        return 4 - corners - static_cast<int>(choose(k, 4 - corners, true).size());
    };

    // Collapses, shortest first among the arcs of at most unresolvedEdges mesh edges or shorter than collapseLength:
    // each taken only if it leaves the charts it touches closer to quads than they were, and never if it would make
    // a chart touch one corner twice. A vertex blurred into two junctions passes (the two charts along the arc each
    // lose the corner too many, or one loses it and takes its own on dS instead); a short true side does not.
    std::vector<std::pair<double, int>> candidates;
    for (int i = 0; i < static_cast<int>(arcs.size()); ++i) {
        const Arc &a = arcs[i];
        if (a.u < 0 || a.u == a.v || (a.across >= 0 && a.across < a.chart)) continue;
        double smaller = area[a.chart];
        if (a.across >= 0) smaller = std::min(smaller, area[a.across]);
        const double rel = a.length / std::sqrt(std::max(smaller, 1e-300));
        if (a.edges <= opts_.unresolvedEdges || rel < opts_.collapseLength) candidates.emplace_back(rel, i);
    }
    std::sort(candidates.begin(), candidates.end());
    // Two passes over them. The first takes a collapse only where it leaves the charts closer to quads (an
    // unresolved arc on a tie too). The second takes, for a chart still with a corner too many, any collapse that
    // leaves things no worse: a chart with an excess corner cannot be repaired, and a neighbour left a corner short
    // can (dissolved, Stage 5), so a tie there is a step towards a layout that covers.
    for (int pass = 0; pass < 2; ++pass)
        for (const auto &[rel, i] : candidates) {
            const int gu = group(arcs[i].u), gv = group(arcs[i].v);
            if (gu == gv) continue;
            std::set<int> atU, affected;
            for (const int m : members[gu]) atU.insert(chartsAt[m].begin(), chartsAt[m].end());
            for (const int m : members[gv])
                for (const int k : chartsAt[m])
                    if (atU.count(k) && outlineOk[k]) affected.insert(k);
            bool ok = true, excess = false;
            int gain = 0;
            std::map<int, int> after;
            for (const int k : affected) {
                const std::vector<int> before = runsOf(k), merged = runsOf(k, gv, gu);
                if (pinched(merged) && !pinched(before)) ok = false;
                after[k] = static_cast<int>(merged.size());
                gain += badness(k, static_cast<int>(before.size())) - badness(k, after[k]);
                excess |= before.size() > 4;
            }
            const int need = pass == 1 ? (excess ? 0 : 1) : (arcs[i].edges <= opts_.unresolvedEdges ? 0 : 1);
            if (!ok || gain < need) continue;
            groupParent_[gv] = gu;
            std::vector<int> &mu = members[gu];
            mu.insert(mu.end(), members[gv].begin(), members[gv].end());
            members.erase(gv);
            for (const auto &[k, c] : after) count[k] = c;
            ++report_.collapsedArcs;
        }

    // Corners on dS.
    std::vector<std::vector<int>> added(K);
    for (int k = 0; k < K; ++k) {
        if (!outlineOk[k] || count[k] >= 4) continue;
        added[k] = choose(k, 4 - count[k], false);
        count[k] += static_cast<int>(added[k].size());
        status_[k].added = static_cast<int>(added[k].size());
        report_.addedCorners += status_[k].added;
    }

    std::vector<int> drop;
    std::vector<char> isBlock(K, 0);
    for (int k = 0; k < K; ++k) {
        if (!outlineOk[k]) continue;
        ChartStatus &st = status_[k];
        const std::vector<int> runs = runsOf(k);
        st.corners = static_cast<int>(runs.size() + added[k].size());
        if (pinched(runs)) {
            st.problem = "a corner twice on its outline";
            continue;
        }
        if (st.corners != 4) {
            st.problem = std::to_string(st.corners) + " corners";
            if (st.corners < 4) drop.push_back(k);
            continue;
        }
        isBlock[k] = 1;
    }

    // ---- Stage 6 ----
    // Macro vertices: one per corner group, at its one junction, or the one on dS nearest the group's mean, or the
    // mean of a group inside the domain; and one per corner added on dS.
    std::map<int, int> macroOfGroup;
    auto macroOf = [&](int g) {
        const auto it = macroOfGroup.find(g);
        if (it != macroOfGroup.end()) return it->second;
        Point mean{0.0, 0.0};
        for (const int m : members[g]) mean = mean + M.vertices[m];
        mean = mean / static_cast<double>(members[g].size());
        int onDS = -1;
        for (const int m : members[g])
            if (M.isBoundaryVertex[m] &&
                (onDS < 0 || normP(M.vertices[m] - mean) < normP(M.vertices[onDS] - mean)))
                onDS = m;
        BlockDecomposition::MacroVertex v;
        v.p = onDS >= 0 ? M.vertices[onDS] : mean;
        v.onBoundary = onDS >= 0;
        const int id = static_cast<int>(blocks_.vertices.size());
        blocks_.vertices.push_back(v);
        macroOfGroup[g] = id;
        return id;
    };
    std::map<int, int> loopVertexOfMacro;   // a dS macro vertex's mesh vertex, for walking dS
    const std::vector<BoundaryFeatures::Loop> &loops = features_->loops();
    std::vector<std::pair<int, int>> loopOf(nv, {-1, -1});
    for (int l = 0; l < static_cast<int>(loops.size()); ++l)
        for (int i = 0; i < static_cast<int>(loops[l].vertices.size()); ++i)
            loopOf[loops[l].vertices[i]] = {l, i};

    std::vector<int> edgeOfArc(arcs.size(), -1);
    std::vector<int> blockOf(K, -1);
    for (int k = 0; k < K; ++k) {
        if (!isBlock[k]) continue;
        blockOf[k] = static_cast<int>(blocks_.blocks.size());
        blocks_.blocks.emplace_back();
    }
    // Corner events along each block's outline: a run of arcs starting in one group is entered at its first arc's
    // start and left at its last arc's start; an added corner is entered and left at its own vertex.
    struct Event {
        int entry, exit, macro;
    };
    std::vector<std::vector<Event>> eventsOf(K);
    for (int k = 0; k < K; ++k) {
        if (!isBlock[k]) continue;
        const Cycle &c = cycles[k][0];
        const int n = static_cast<int>(c.vertex.size());
        std::vector<Event> &events = eventsOf[k];
        const std::vector<int> &ka = chartArcs[k];
        const int m = static_cast<int>(ka.size());
        int first = 0;
        if (arcs[ka[0]].u >= 0) {
            while (first < m && group(arcs[ka[first]].u) == group(arcs[ka[(first + m - 1) % m]].u)) ++first;
            if (first == m) first = 0;
        }
        for (int s = 0; s < m && arcs[ka[0]].u >= 0; ++s) {
            const Arc &a = arcs[ka[(first + s) % m]];
            const int g = group(a.u);
            if (!events.empty() && group(c.vertex[events.back().entry]) == g) {
                events.back().exit = a.begin;
                continue;
            }
            events.push_back({a.begin, a.begin, macroOf(g)});
            const int mv = events.back().macro;
            if (blocks_.vertices[mv].onBoundary) {
                for (const int mm : members[g])
                    if (M.isBoundaryVertex[mm] && normP(M.vertices[mm] - blocks_.vertices[mv].p) == 0.0)
                        loopVertexOfMacro[mv] = mm;
            }
        }
        const int origin = events.empty() ? 0 : events.front().entry;
        for (const int i : added[k]) {
            BlockDecomposition::MacroVertex v;
            v.p = M.vertices[c.vertex[i]];
            v.onBoundary = true;
            const int id = static_cast<int>(blocks_.vertices.size());
            blocks_.vertices.push_back(v);
            loopVertexOfMacro[id] = c.vertex[i];
            events.push_back({i, i, id});
        }
        std::sort(events.begin(), events.end(), [&](const Event &a, const Event &b) {
            return (a.entry - origin + n) % n < (b.entry - origin + n) % n;
        });
    }

    // A layout vertex on dS within cornerSnapEdges mesh edges of a corner of the model that turns by cornerSnapAngle
    // or more, along the boundary loop and with no other layout vertex in between, is moved onto it: the sharpest
    // such corner, a reflex one before any convex one. (The angle keeps the vertices of a polygonal hole out of it:
    // a small hole drawn with a dozen chords turns 30 degrees at every one, and every one is a corner by the
    // 20-degree rule.) A wall the decoder put two edges from a reflex corner otherwise leaves a block side bending a
    // quarter turn the wrong way through it, which no transfinite grid survives, and a convex corner a step off is
    // a block corner a step off.
    if (opts_.cornerSnapEdges > 0) {
        std::vector<int> macroAtLoopVertex(nv, -1);
        for (const auto &[mv, lv] : loopVertexOfMacro) macroAtLoopVertex[lv] = mv;
        std::vector<double> turn(nv, 0.0);
        for (const BoundaryFeatures::Corner &c : features_->corners())
            if (std::fabs(M_PI - c.angle) >= opts_.cornerSnapAngle * M_PI / 180.0)
                turn[c.vertex] = std::fabs(M_PI - c.angle) + (c.angle > M_PI ? M_PI : 0.0);
        for (auto &[mv, lv] : loopVertexOfMacro) {
            const auto [l, i0] = loopOf[lv];
            if (l < 0 || turn[lv] > 0.0) continue;
            const std::vector<int> &lp = loops[l].vertices;
            const int ml = static_cast<int>(lp.size());
            // A reflex corner is looked for twice as far: no quad layout can do without a vertex there, so a wall
            // anywhere near one belongs on it.
            int best = -1;
            for (const int dir : {1, -1})
                for (int step = 1; step <= 2 * opts_.cornerSnapEdges; ++step) {
                    const int w = lp[((i0 + dir * step) % ml + ml) % ml];
                    if (macroAtLoopVertex[w] >= 0) break;
                    const bool reflex = turn[w] > M_PI;
                    if (turn[w] > 0.0 && (reflex || step <= opts_.cornerSnapEdges) && (best < 0 || turn[w] > turn[best]))
                        best = w;
                }
            if (best < 0) continue;
            macroAtLoopVertex[lv] = -1;
            macroAtLoopVertex[best] = mv;
            lv = best;
            blocks_.vertices[mv].p = M.vertices[best];
            ++report_.snappedCorners;
        }
    }

    // Boundary distances, for keeping smoothed sides off dS.
    BoundaryGrid boundary(M);

    for (int k = 0; k < K; ++k) {
        if (!isBlock[k]) continue;
        const Cycle &c = cycles[k][0];
        const int n = static_cast<int>(c.vertex.size());
        const std::vector<int> &ka = chartArcs[k];
        const std::vector<Event> &events = eventsOf[k];
        double outlineEdge = 0.0;
        for (int i = 0; i < n; ++i) outlineEdge += edgeLength(c, i);
        outlineEdge /= n;
        BlockDecomposition::Block &blk = blocks_.blocks[blockOf[k]];
        blk.material = M.triangleMatId.empty() ? 0 : M.triangleMatId[charts_->regions()[k].peak];
        for (int s = 0; s < 4; ++s) {
            const Event &a = events[s], &b = events[(s + 1) % 4];
            blk.corners[s] = a.macro;
            const int across = c.across[a.exit];
            // The arc this side is, for an interior side; sides on dS are never shared.
            int arc = -1;
            for (const int x : ka)
                if (arcs[x].begin == a.exit && arcs[x].u >= 0) arc = x;
            if (across >= 0 && arc >= 0 && arcs[arc].twin >= 0 && edgeOfArc[arcs[arc].twin] >= 0) {
                BlockDecomposition::MacroEdge &e = blocks_.edges[edgeOfArc[arcs[arc].twin]];
                e.blockB = blockOf[k];
                e.sideB = s;
                blk.edges[s] = edgeOfArc[arcs[arc].twin];
                blk.flip[s] = true;
                continue;
            }
            BlockDecomposition::MacroEdge e;
            e.from = a.macro;
            e.to = b.macro;
            e.blockA = blockOf[k];
            e.sideA = s;
            e.boundary = across == -1;
            if (e.boundary) {
                const auto ia = loopVertexOfMacro.find(a.macro), ib = loopVertexOfMacro.find(b.macro);
                bool walked = false;
                if (ia != loopVertexOfMacro.end() && ib != loopVertexOfMacro.end()) {
                    const auto [la, pa] = loopOf[ia->second];
                    const auto [lb, pb] = loopOf[ib->second];
                    if (la >= 0 && la == lb && pa != pb) {
                        const std::vector<int> &lv = loops[la].vertices;
                        const int ml = static_cast<int>(lv.size());
                        for (int i = pa;; i = (i + 1) % ml) {
                            e.points.push_back(M.vertices[lv[i]]);
                            if (i == pb) break;
                        }
                        walked = true;
                    }
                }
                if (!walked) {
                    for (int i = a.exit;; i = (i + 1) % n) {
                        e.points.push_back(M.vertices[c.vertex[i]]);
                        if (i == b.entry) break;
                    }
                    e.points.front() = blocks_.vertices[a.macro].p;
                    e.points.back() = blocks_.vertices[b.macro].p;
                    report_.messages.push_back("chart " + std::to_string(k) + " side " + std::to_string(s) +
                                               ": dS side not on one boundary loop; its outline used instead");
                }
            } else {
                std::vector<int> run;
                for (int i = (a.exit + 1) % n; i != b.entry; i = (i + 1) % n) run.push_back(c.vertex[i]);
                // Next to a corner on dS, the run is cut back to where it leaves dS by endClearance mean outline
                // edges: a decoded wall often meets dS at a grazing angle and runs along it for a stretch first,
                // and the sliver of block between that stretch and dS is folded by any grid laid on it.
                size_t lo = 0, hi = run.size();
                const double clearance = opts_.endClearance * outlineEdge;
                if (blocks_.vertices[a.macro].onBoundary)
                    while (lo < hi && (M.isBoundaryVertex[run[lo]] ||
                                       boundary.distance(M.vertices[run[lo]], clearance) < clearance))
                        ++lo;
                if (blocks_.vertices[b.macro].onBoundary)
                    while (hi > lo && (M.isBoundaryVertex[run[hi - 1]] ||
                                       boundary.distance(M.vertices[run[hi - 1]], clearance) < clearance))
                        --hi;
                e.points.push_back(blocks_.vertices[a.macro].p);
                for (size_t i = lo; i < hi; ++i) e.points.push_back(M.vertices[run[i]]);
                e.points.push_back(blocks_.vertices[b.macro].p);
                // Smoothed, end points held, and no point moved by more than half its distance to dS: a side
                // running a few triangles off a small hole would otherwise be pulled across it, every pass taking
                // the polyline a little further towards the chord.
                const std::vector<Point> raw = e.points;
                for (int pass = 0; pass < opts_.sidePasses; ++pass) {
                    std::vector<Point> q = e.points;
                    for (size_t i = 1; i + 1 < q.size(); ++i)
                        q[i] = e.points[i] * 0.5 + (e.points[i - 1] + e.points[i + 1]) * 0.25;
                    e.points.swap(q);
                }
                for (size_t i = 1; i + 1 < e.points.size(); ++i) {
                    const Point d = e.points[i] - raw[i];
                    const double len = normP(d), limit = 0.5 * boundary.distance(raw[i], len);
                    if (len > limit) e.points[i] = raw[i] + d * (limit / len);
                }
            }
            const int id = static_cast<int>(blocks_.edges.size());
            if (arc >= 0) edgeOfArc[arc] = id;
            blk.edges[s] = id;
            blk.flip[s] = false;
            blocks_.edges.push_back(std::move(e));
        }
        status_[k].block = true;
    }
    for (BlockDecomposition::MacroEdge &e : blocks_.edges) {
        ++blocks_.vertices[e.from].valence;
        ++blocks_.vertices[e.to].valence;
        if (!e.boundary && e.blockB < 0) ++report_.unmatchedSides;
    }
    for (const BlockDecomposition::MacroVertex &v : blocks_.vertices) report_.boundaryVertices += v.onBoundary;
    for (int k = 0; k < K; ++k) report_.failedCharts += !status_[k].block && !status_[k].dissolved;
    report_.blocks = static_cast<int>(blocks_.blocks.size());
    report_.macroVertices = static_cast<int>(blocks_.vertices.size());
    report_.macroEdges = static_cast<int>(blocks_.edges.size());
    report_.covers = blocks_.covers();
    return drop;
}
