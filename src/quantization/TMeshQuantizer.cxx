#include "quantization/TMeshQuantizer.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <utility>

namespace {

const double INF = std::numeric_limits<double>::infinity();

}  // namespace

TMeshQuantizer::TMeshQuantizer(QuantTMesh &tmesh, Options opts)
    : tm_(tmesh), opts_(opts) {
    if (!tm_.finalized()) {
        throw std::runtime_error("TMeshQuantizer: T-mesh not finalized");
    }
    eta_ = static_cast<double>(std::max<size_t>(tm_.edges.size(), 1));
    for (size_t e = 0; e < tm_.edges.size(); ++e) {
        for (int k = 0; k < 2; ++k) {
            if (tm_.edges[e].row[k] == QuantTMesh::PHANTOM) {
                sinks_.push_back(2 * static_cast<int>(e) + k);
            }
        }
    }
}

int TMeshQuantizer::rowOfNode(int n) const {
    return tm_.edges[n >> 1].row[n & 1];
}

// The three-tier weights of Sec. 6.2 with the Sec. 10.1 override: an edge
// already at the lower bound must never be decremented. Strict preference
// goes to edges that stay shorter than ideal after the move, then to edges
// that become longer than ideal, then to edges already too long.
double TMeshQuantizer::weight(int e, int s) const {
    const QuantTMesh::Edge &E = tm_.edges[e];
    if (s < 0 && E.x <= opts_.minLength) return INF;
    const double d = s > 0 ? E.xIdeal - E.x : E.x - E.xIdeal;
    if (d >= 1.0) return 1.0 / (d + 1.0);
    if (d >= 0.0) return eta_ / (d + 1.0);
    return eta_ * eta_ * (1.0 - d);
}

// Arcs out of node (e, k): the row r = edges[e].row[k] is unbalanced by e's
// coefficient there; every edge g of opposite sign in r can rebalance it,
// which in turn unbalances g's other row -- node (g, other slot of r).
template <typename F>
void TMeshQuantizer::forEachOut(int n, F &&f) const {
    const int e = n >> 1, k = n & 1;
    const int r = tm_.edges[e].row[k];
    if (r == QuantTMesh::PHANTOM) return;  // sink
    const QuantTMesh::Row &row = tm_.rows[r];
    for (int g : tm_.edges[e].sign[k] > 0 ? row.neg : row.pos) {
        const int kg = tm_.edges[g].row[0] == r ? 0 : 1;
        f(2 * g + (1 - kg));
    }
}

// Arcs into node (e, k) arrive via e's other row: they come from every edge
// of opposite sign there. A boundary edge's real-row node has no in-arcs --
// those nodes are the sources of D, mirroring the sinks.
template <typename F>
void TMeshQuantizer::forEachIn(int n, F &&f) const {
    const int e = n >> 1, k = n & 1;
    const int r = tm_.edges[e].row[1 - k];
    if (r == QuantTMesh::PHANTOM) return;  // source
    const QuantTMesh::Row &row = tm_.rows[r];
    for (int g : tm_.edges[e].sign[1 - k] > 0 ? row.neg : row.pos) {
        const int kg = tm_.edges[g].row[0] == r ? 0 : 1;
        f(2 * g + kg);
    }
}

void TMeshQuantizer::dijkstra(int n0, int s, int dir, std::vector<double> &dist,
                              std::vector<int> &parent) const {
    const size_t nn = 2 * tm_.edges.size();
    dist.assign(nn, INF);
    parent.assign(nn, -1);

    using Item = std::pair<double, int>;
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> queue;
    dist[n0] = weight(n0 >> 1, s);
    if (dist[n0] == INF) return;
    queue.push({dist[n0], n0});

    while (!queue.empty()) {
        const double d = queue.top().first;
        const int n = queue.top().second;
        queue.pop();
        if (d > dist[n]) continue;
        auto relax = [&](int m) {
            const double w = weight(m >> 1, s);
            if (w == INF) return;
            if (d + w < dist[m]) {
                dist[m] = d + w;
                parent[m] = n;
                queue.push({dist[m], m});
            }
        };
        if (dir > 0) {
            forEachOut(n, relax);
        } else {
            forEachIn(n, relax);
        }
    }
}

bool TMeshQuantizer::findGeneratingVector(int e0, int s,
                                          std::vector<int> &v) const {
    double bestCost = INF;
    std::vector<int> bestNodes;
    std::vector<double> df, db;
    std::vector<int> pf, pb;

    auto trace = [](const std::vector<int> &parent, int n) {
        std::vector<int> walk;
        for (; n >= 0; n = parent[n]) walk.push_back(n);
        return walk;  // from n back to the search root
    };

    for (int k0 = 0; k0 < 2; ++k0) {
        const int n0 = 2 * e0 + k0;
        dijkstra(n0, s, +1, df, pf);
        if (df[n0] == INF) return false;  // e0 itself is not incrementable

        // (a) Elementary circuit through n0: cheapest settled in-neighbour
        // plus the closing arc. Dijkstra paths are node-simple, so the
        // circuit is elementary; an edge may still contribute both its
        // nodes, which is the legitimate v[e] = 2 case.
        forEachIn(n0, [&](int p) {
            if (df[p] < bestCost) {
                bestCost = df[p];
                bestNodes = trace(pf, p);
            }
        });

        // (b) Source-to-sink path through n0, the boundary case of
        // Sec. 9.2: forward to the cheapest sink, backward to the cheapest
        // source, n0 counted once. n0 may itself be the source or sink.
        dijkstra(n0, s, -1, db, pb);
        double toSink = INF;
        int sink = -1;
        for (int t : sinks_) {
            if (df[t] < toSink) {
                toSink = df[t];
                sink = t;
            }
        }
        double fromSource = INF;
        int source = -1;
        for (size_t e = 0; e < tm_.edges.size(); ++e) {
            for (int k = 0; k < 2; ++k) {
                // A source's in-arcs would come via its other row.
                if (tm_.edges[e].row[1 - k] != QuantTMesh::PHANTOM) continue;
                const int u = 2 * static_cast<int>(e) + k;
                if (db[u] < fromSource) {
                    fromSource = db[u];
                    source = u;
                }
            }
        }
        if (sink >= 0 && source >= 0 &&
            toSink + fromSource - df[n0] < bestCost) {
            bestCost = toSink + fromSource - df[n0];
            bestNodes = trace(pf, sink);  // sink .. n0
            std::vector<int> back = trace(pb, source);  // source .. n0
            // n0 closes both traces; keep it once.
            bestNodes.insert(bestNodes.end(), back.begin(), back.end() - 1);
        }
    }

    if (bestCost == INF) return false;
    v.assign(tm_.edges.size(), 0);
    for (int n : bestNodes) ++v[n >> 1];
    return true;
}

// Stage I (Sec. 6.3): grow the trivial solution strip by strip until every
// edge reaches the lower bound. Adding never decreases any x, so at most
// minLength * #edges strips are needed.
void TMeshQuantizer::stageI(Report &report) {
    std::vector<int> v;
    // Edges D has no circuit or source-to-sink path through. Appendix A.2
    // rules these out for a T-mesh traced on the separatrices of a seamless
    // parametrization, but a T-mesh assembled some other way -- a template
    // block decomposition, say -- can contain them. They are recorded once
    // and then left alone: no consistent assignment can lift them off zero.
    std::vector<char> forced(tm_.edges.size(), 0);
    while (true) {
        int e0 = -1;
        double best = INF;
        for (size_t e = 0; e < tm_.edges.size(); ++e) {
            if (tm_.edges[e].x >= opts_.minLength || forced[e]) continue;
            const double w = weight(static_cast<int>(e), +1);
            if (w < best) {
                best = w;
                e0 = static_cast<int>(e);
            }
        }
        if (e0 < 0) return;
        if (!findGeneratingVector(e0, +1, v)) {
            forced[e0] = 1;
            ++report.forcedZeroEdges;
            continue;
        }
        for (size_t e = 0; e < v.size(); ++e) tm_.edges[e].x += v[e];
        ++report.stage1Vectors;
    }
}

// Stage II (Sec. 6.3): repeatedly pick the edge most eager to change, build
// the cheapest strip through it, and keep the move iff no edge drops below
// the bound and the objective improves. The objective strictly decreases
// over integer states, so this terminates.
void TMeshQuantizer::stageII(Report &report) {
    std::vector<std::pair<double, int>> candidates;
    std::vector<int> v;
    double objective = tm_.objective();

    while (opts_.maxStage2Moves < 0 || report.stage2Moves < opts_.maxStage2Moves) {
        candidates.clear();
        for (size_t e = 0; e < tm_.edges.size(); ++e) {
            const int ei = static_cast<int>(e);
            const double w = std::min(weight(ei, +1), weight(ei, -1));
            if (w < INF) candidates.push_back({w, ei});
        }
        std::sort(candidates.begin(), candidates.end());

        bool improved = false;
        for (auto [w, e0] : candidates) {
            const int s = weight(e0, +1) < weight(e0, -1) ? +1 : -1;
            if (!findGeneratingVector(e0, s, v)) continue;
            ++report.stage2Tried;

            for (size_t e = 0; e < v.size(); ++e) tm_.edges[e].x += s * v[e];
            // The INF weight already steers subtractions away from edges at
            // the bound, but a strip crossing an edge twice can still
            // overshoot; this scan is the entire validity check.
            bool valid = true;
            if (s < 0) {
                for (size_t e = 0; e < v.size(); ++e) {
                    if (v[e] > 0 && tm_.edges[e].x < opts_.minLength) {
                        valid = false;
                        break;
                    }
                }
            }
            const double after = tm_.objective();
            if (valid && after < objective) {
                objective = after;
                improved = true;
                ++report.stage2Moves;
                break;
            }
            for (size_t e = 0; e < v.size(); ++e) tm_.edges[e].x -= s * v[e];
        }
        if (!improved) return;
    }
}

TMeshQuantizer::Report TMeshQuantizer::run() {
    for (QuantTMesh::Edge &e : tm_.edges) e.x = 0;
    Report report;
    stageI(report);
    stageII(report);
    report.objective = tm_.objective();
    report.consistent = tm_.consistent();
    return report;
}
