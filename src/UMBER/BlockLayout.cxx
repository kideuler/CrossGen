#include "UMBER/BlockLayout.hxx"

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <unordered_set>

namespace {

// The point of a polyline at p = (index of the point before it) + (fraction of
// the way along that step).
Point pointAtParam(const std::vector<Point> &path, double p) {
    const int n = static_cast<int>(path.size());
    if (n == 0) return Point{0.0, 0.0};
    p = std::min(std::max(p, 0.0), static_cast<double>(n - 1));
    const int i = std::min(n - 2, static_cast<int>(std::floor(p)));
    if (i < 0) return path[0];
    const double t = p - i;
    return path[i] * (1.0 - t) + path[i + 1] * t;
}

// The piece of a polyline between two parameters, ends included.
std::vector<Point> sliceOfPolyline(const std::vector<Point> &path, double p0, double p1) {
    std::vector<Point> pts;
    pts.push_back(pointAtParam(path, p0));
    for (int k = static_cast<int>(std::floor(p0)) + 1; k <= static_cast<int>(std::ceil(p1)); ++k) {
        if (k <= p0 + 1e-12 || k >= p1 - 1e-12) continue;
        if (k >= 0 && k < static_cast<int>(path.size())) pts.push_back(path[k]);
    }
    pts.push_back(pointAtParam(path, p1));
    return pts;
}

} // namespace

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
double BlockLayout::toleranceOf(const MotorcycleGraph &graph) {
    const Mesh &m = graph.getMesh();
    double total = 0.0;
    for (const auto &e : m.edges) total += normP(m.vertices[e[1]] - m.vertices[e[0]]);
    const double mean = m.edges.empty() ? 1.0 : total / static_cast<double>(m.edges.size());
    return 1e-6 * mean;
}

BlockLayout::BlockLayout(const MotorcycleGraph &g)
    : graph(&g), mesh(&g.getMesh()), tol(toleranceOf(g)), layout_(toleranceOf(g)) {}

// ---------------------------------------------------------------------------
// buildBoundaryLoops()
//
// Oriented with the material on the left, which is the one thing the face walk
// needs of them: a cycle that runs along a boundary arc backwards is on the
// outside of the model, and that is how the unbounded face -- and the void
// inside each hole -- is told from a block.
// ---------------------------------------------------------------------------
void BlockLayout::buildBoundaryLoops() {
    loops.clear();
    loopOfVertex.assign(mesh->vertices.size(), {-1, -1});

    std::unordered_map<int, int> nextOf;
    for (int e : mesh->boundaryEdges) {
        const int t = (mesh->edgeTriangles[e][0] >= 0) ? mesh->edgeTriangles[e][0]
                                                       : mesh->edgeTriangles[e][1];
        if (t < 0) continue;
        int a = mesh->edges[e][0], b = mesh->edges[e][1];
        int c = -1;
        for (int i = 0; i < 3; ++i) {
            const int v = mesh->triangles[t][i];
            if (v != a && v != b) c = v;
        }
        if (c < 0) continue;
        if (cross2(mesh->vertices[b] - mesh->vertices[a],
                   mesh->vertices[c] - mesh->vertices[a]) < 0.0) std::swap(a, b);
        nextOf[a] = b;
    }

    std::unordered_set<int> seen;
    for (int e : mesh->boundaryEdges) {
        for (int k = 0; k < 2; ++k) {
            const int start = mesh->edges[e][k];
            if (seen.count(start) || !nextOf.count(start)) continue;
            std::vector<int> loop;
            int v = start;
            while (!seen.count(v)) {
                seen.insert(v);
                loop.push_back(v);
                auto it = nextOf.find(v);
                if (it == nextOf.end()) { loop.clear(); break; }
                v = it->second;
            }
            if (loop.size() < 3) continue;
            const int l = static_cast<int>(loops.size());
            for (size_t i = 0; i < loop.size(); ++i)
                loopOfVertex[loop[i]] = {l, static_cast<int>(i)};
            loops.push_back(std::move(loop));
        }
    }
    report_.boundaryLoops = static_cast<int>(loops.size());
}

// ---------------------------------------------------------------------------
// addNode() / findNode()
// ---------------------------------------------------------------------------
int BlockLayout::findNode(const Point &pos) const {
    for (size_t i = 0; i < nodes.size(); ++i)
        if (normP(nodes[i].pos - pos) <= tol) return static_cast<int>(i);
    return -1;
}

int BlockLayout::addNode(const Point &pos, QuadLayout::NodeKind kind) {
    const int found = findNode(pos);
    if (found >= 0) {
        // The more structural kind wins, the way QuadLayout ranks them: a
        // crossing that lands on a corner of the polysquare is the corner.
        if (static_cast<int>(kind) < static_cast<int>(nodes[found].kind))
            nodes[found].kind = kind;
        ++report_.mergedNodes;
        return found;
    }
    QuadLayout::Node n;
    n.pos = pos;
    n.kind = kind;
    n.source = -1;
    nodes.push_back(std::move(n));
    return static_cast<int>(nodes.size()) - 1;
}

int BlockLayout::addArc(std::vector<Point> pts, int a, int b, int ray, bool onBoundary) {
    // Points closer together than the node tolerance carry no direction, and a
    // direction is exactly what the sort around a node needs.
    std::vector<Point> clean;
    clean.reserve(pts.size());
    for (const Point &p : pts)
        if (clean.empty() || normP(p - clean.back()) > tol * 0.5) clean.push_back(p);
    if (a < 0 || b < 0 || clean.size() < 2) { ++report_.degenerateArcs; return -1; }
    clean.front() = nodes[a].pos;
    clean.back() = nodes[b].pos;
    if (normP(clean.front() - clean.back()) < tol && clean.size() < 3) {
        ++report_.degenerateArcs;
        return -1;
    }

    QuadLayout::Arc arc;
    arc.length = 0.0;
    for (size_t i = 1; i < clean.size(); ++i) arc.length += normP(clean[i] - clean[i - 1]);
    arc.pts = std::move(clean);
    arc.a = a;
    arc.b = b;
    arc.separatrix = ray;
    arc.onBoundary = onBoundary;
    arcs.push_back(std::move(arc));
    return static_cast<int>(arcs.size()) - 1;
}

// ---------------------------------------------------------------------------
// collectNodes()  --  the corners, the exits and the crossings
//
// Corners first, so that anything landing on one merges into it rather than the
// other way round, and each node is filed under the lines it sits on as it is
// made: a node is only of use to the arc builder once it knows where along a
// ray, or around a boundary loop, the node cuts it.
// ---------------------------------------------------------------------------
void BlockLayout::collectNodes() {
    const auto &traces = graph->getTraces();
    const auto &mgNodes = graph->getNodes();

    splitsOfRay.assign(traces.size(), {});
    splitsOfLoop.assign(loops.size(), {});

    auto placeOnLoop = [&](int node, int boundaryEdge, double along) {
        if (boundaryEdge < 0) { ++report_.unplacedNodes; return; }
        const int va = mesh->edges[boundaryEdge][0], vb = mesh->edges[boundaryEdge][1];
        if (va >= static_cast<int>(loopOfVertex.size()) ||
            vb >= static_cast<int>(loopOfVertex.size())) { ++report_.unplacedNodes; return; }
        const auto &ia = loopOfVertex[va];
        const auto &ib = loopOfVertex[vb];
        if (ia.first < 0 || ia.first != ib.first) { ++report_.unplacedNodes; return; }
        const int l = ia.first;
        const int n = static_cast<int>(loops[l].size());
        // Which way round the loop this edge runs.
        if ((ia.second + 1) % n == ib.second)
            splitsOfLoop[l].push_back(Split{ia.second + along, node});
        else if ((ib.second + 1) % n == ia.second)
            splitsOfLoop[l].push_back(Split{ib.second + (1.0 - along), node});
        else
            ++report_.unplacedNodes;
    };

    for (int pass = 0; pass < 3; ++pass) {
        for (const auto &mn : mgNodes) {
            if (pass == 0 && mn.kind != MotorcycleGraph::Node::Corner) continue;
            if (pass == 1 && mn.kind != MotorcycleGraph::Node::BoundaryEnd) continue;
            if (pass == 2 && mn.kind != MotorcycleGraph::Node::Crossing) continue;

            switch (mn.kind) {
                case MotorcycleGraph::Node::Corner: {
                    const int id = addNode(mn.xy, QuadLayout::NodeKind::BoundaryCorner);
                    if (mn.vertex >= 0 && mn.vertex < static_cast<int>(loopOfVertex.size()) &&
                        loopOfVertex[mn.vertex].first >= 0) {
                        const auto &lv = loopOfVertex[mn.vertex];
                        splitsOfLoop[lv.first].push_back(
                            Split{static_cast<double>(lv.second), id});
                    } else {
                        ++report_.unplacedNodes;
                    }
                    break;
                }
                case MotorcycleGraph::Node::BoundaryEnd: {
                    const int id = addNode(mn.xy, QuadLayout::NodeKind::BoundaryHit);
                    placeOnLoop(id, mn.boundaryEdge, mn.alongEdge);
                    if (mn.ray[0] >= 0 && mn.ray[0] < static_cast<int>(splitsOfRay.size()))
                        splitsOfRay[mn.ray[0]].push_back(Split{mn.param[0], id});
                    break;
                }
                case MotorcycleGraph::Node::Crossing: {
                    const int id = addNode(mn.xy, QuadLayout::NodeKind::Crossing);
                    for (int k = 0; k < 2; ++k)
                        if (mn.ray[k] >= 0 && mn.ray[k] < static_cast<int>(splitsOfRay.size()))
                            splitsOfRay[mn.ray[k]].push_back(Split{mn.param[k], id});
                    break;
                }
            }
        }
    }

    // Each ray starts at the corner it was launched from and, if it never got
    // out of the model, ends nowhere. The first is already a node -- a reflex
    // corner is a corner of the polysquare -- and only has to be found; the
    // second is a free end, which the arrangement still has to carry so that
    // what it spoils is visible instead of silently missing.
    for (size_t r = 0; r < traces.size(); ++r) {
        if (traces[r].size() < 2) continue;
        const int start = findNode(traces[r].front());
        if (start >= 0) splitsOfRay[r].push_back(Split{0.0, start});
        else ++report_.unplacedNodes;

        const auto &exit = graph->getRayExitEdge();
        if (r >= exit.size() || exit[r] < 0) {
            const int end = addNode(traces[r].back(), QuadLayout::NodeKind::Dangling);
            splitsOfRay[r].push_back(Split{static_cast<double>(traces[r].size() - 1), end});
            ++report_.danglingEnds;
        }
    }
}

// ---------------------------------------------------------------------------
// buildRayArcs()
// ---------------------------------------------------------------------------
void BlockLayout::buildRayArcs() {
    const auto &traces = graph->getTraces();
    for (size_t r = 0; r < traces.size(); ++r) {
        if (traces[r].size() < 2) continue;
        auto splits = splitsOfRay[r];
        std::sort(splits.begin(), splits.end(),
                  [](const Split &x, const Split &y) { return x.p < y.p; });
        splits.erase(std::unique(splits.begin(), splits.end(),
                                 [](const Split &x, const Split &y) {
                                     return x.node == y.node || std::fabs(x.p - y.p) < 1e-12;
                                 }),
                     splits.end());

        for (size_t i = 1; i < splits.size(); ++i)
            addArc(sliceOfPolyline(traces[r], splits[i - 1].p, splits[i].p),
                   splits[i - 1].node, splits[i].node, static_cast<int>(r), false);
    }
}

// ---------------------------------------------------------------------------
// buildBoundaryArcs()
// ---------------------------------------------------------------------------
void BlockLayout::buildBoundaryArcs() {
    for (size_t l = 0; l < loops.size(); ++l) {
        const auto &loop = loops[l];
        const int n = static_cast<int>(loop.size());

        auto splits = splitsOfLoop[l];
        if (splits.empty()) {
            // A loop nothing lands on -- a hole no iso-line reached and that
            // has no corner of its own, which the polysquare does not produce
            // but which costs nothing to carry. It still has to be an arc or
            // the material beside it is not enclosed.
            splits.push_back(Split{0.0, addNode(mesh->vertices[loop[0]],
                                                QuadLayout::NodeKind::BoundaryHit)});
        }
        std::sort(splits.begin(), splits.end(),
                  [](const Split &x, const Split &y) { return x.p < y.p; });
        splits.erase(std::unique(splits.begin(), splits.end(),
                                 [](const Split &x, const Split &y) {
                                     return x.node == y.node || std::fabs(x.p - y.p) < 1e-12;
                                 }),
                     splits.end());

        auto pointAt = [&](double p) {
            const int i = static_cast<int>(std::floor(p)) % n;
            const double t = p - std::floor(p);
            return mesh->vertices[loop[i]] * (1.0 - t) + mesh->vertices[loop[(i + 1) % n]] * t;
        };

        // Round the loop, wrapping the last piece back to the first split.
        for (size_t i = 0; i < splits.size(); ++i) {
            const double p0 = splits[i].p;
            const double p1 = (i + 1 < splits.size()) ? splits[i + 1].p : splits[0].p + n;
            std::vector<Point> pts;
            pts.push_back(pointAt(p0));
            for (int k = static_cast<int>(std::floor(p0)) + 1;
                 k <= static_cast<int>(std::ceil(p1)); ++k) {
                if (k <= p0 + 1e-12 || k >= p1 - 1e-12) continue;
                pts.push_back(mesh->vertices[loop[k % n]]);
            }
            pts.push_back(pointAt(p1 >= n ? p1 - n : p1));
            addArc(std::move(pts), splits[i].node,
                   splits[(i + 1) % splits.size()].node, -1, true);
        }
    }
}

// ---------------------------------------------------------------------------
// build()
// ---------------------------------------------------------------------------
void BlockLayout::build() {
    nodes.clear();
    arcs.clear();
    report_ = Report{};

    buildBoundaryLoops();
    collectNodes();
    buildRayArcs();
    buildBoundaryArcs();

    layout_.rebuild(nodes, arcs);
}
