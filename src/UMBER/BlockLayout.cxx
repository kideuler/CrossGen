#include "UMBER/BlockLayout.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <unordered_set>

namespace {

// ---------------------------------------------------------------------------
// Which triangle of the model a point of the model is in.
//
// Only one question is ever asked of this -- what material does this block sit
// in -- and it is asked once per block, so a uniform grid over the triangle
// bounding boxes is the right structure: built in one pass, and on a mesh of
// any size the cell a query lands in holds a handful of triangles.
//
// It answers -1 for a point outside every triangle rather than snapping to the
// nearest, because the two cases mean different things to the caller: a sample
// that missed is a sample to ignore, and a block whose every sample missed is
// a block whose outline does not lie on the model, which is worth reporting
// rather than papering over.
// ---------------------------------------------------------------------------
class TriangleGrid {
public:
    explicit TriangleGrid(const Mesh &m) : mesh(&m) {
        lo = Point{std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity()};
        hi = Point{-lo[0], -lo[1]};
        for (const Point &p : m.vertices) {
            lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
            hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
        }
        const int nT = static_cast<int>(m.triangles.size());
        if (nT == 0 || !(hi[0] > lo[0]) || !(hi[1] > lo[1])) return;

        // About one triangle per cell, which is where a uniform grid's build
        // cost and its query cost cross over.
        const double area = (hi[0] - lo[0]) * (hi[1] - lo[1]);
        cell = std::sqrt(area / static_cast<double>(nT));
        if (!(cell > 0.0)) cell = std::max(hi[0] - lo[0], hi[1] - lo[1]);
        nx = std::max(1, static_cast<int>((hi[0] - lo[0]) / cell) + 1);
        ny = std::max(1, static_cast<int>((hi[1] - lo[1]) / cell) + 1);
        cells.assign(static_cast<size_t>(nx) * ny, {});

        for (int t = 0; t < nT; ++t) {
            const Triangle &tri = m.triangles[t];
            Point tlo{std::numeric_limits<double>::infinity(),
                      std::numeric_limits<double>::infinity()};
            Point thi{-tlo[0], -tlo[1]};
            for (int k = 0; k < 3; ++k) {
                const Point &p = m.vertices[tri[k]];
                tlo[0] = std::min(tlo[0], p[0]); tlo[1] = std::min(tlo[1], p[1]);
                thi[0] = std::max(thi[0], p[0]); thi[1] = std::max(thi[1], p[1]);
            }
            const int i0 = clampi(static_cast<int>((tlo[0] - lo[0]) / cell), nx);
            const int i1 = clampi(static_cast<int>((thi[0] - lo[0]) / cell), nx);
            const int j0 = clampi(static_cast<int>((tlo[1] - lo[1]) / cell), ny);
            const int j1 = clampi(static_cast<int>((thi[1] - lo[1]) / cell), ny);
            for (int j = j0; j <= j1; ++j)
                for (int i = i0; i <= i1; ++i)
                    cells[static_cast<size_t>(j) * nx + i].push_back(t);
        }
    }

    int locate(const Point &q) const {
        if (cells.empty()) return -1;
        const int i = clampi(static_cast<int>((q[0] - lo[0]) / cell), nx);
        const int j = clampi(static_cast<int>((q[1] - lo[1]) / cell), ny);
        for (const int t : cells[static_cast<size_t>(j) * nx + i])
            if (inside(t, q)) return t;
        return -1;
    }

private:
    static int clampi(int v, int n) { return std::max(0, std::min(n - 1, v)); }

    bool inside(int t, const Point &q) const {
        const Triangle &tri = mesh->triangles[t];
        const Point &a = mesh->vertices[tri[0]];
        const Point &b = mesh->vertices[tri[1]];
        const Point &c = mesh->vertices[tri[2]];
        // Signs against the triangle's own orientation, so a mesh stored
        // clockwise is not reported empty everywhere.
        const double s0 = cross2(b - a, q - a);
        const double s1 = cross2(c - b, q - b);
        const double s2 = cross2(a - c, q - c);
        return (s0 >= 0.0 && s1 >= 0.0 && s2 >= 0.0) ||
               (s0 <= 0.0 && s1 <= 0.0 && s2 <= 0.0);
    }

    const Mesh *mesh = nullptr;
    Point lo{0.0, 0.0}, hi{0.0, 0.0};
    double cell = 1.0;
    int nx = 0, ny = 0;
    std::vector<std::vector<int>> cells;
};

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
// buildInterfaceChains()
//
// The interface network, cut at every vertex the network does not simply carry
// on through: a junction of three or more branches, a free end, or a landing
// on dS. Each chain is a run of vertices, and each mesh edge belongs to
// exactly one of them -- which is why the lookup below is keyed by edge and
// the boundary's is keyed by vertex. A vertex can be on two chains; that is
// what a junction is.
// ---------------------------------------------------------------------------
void BlockLayout::buildInterfaceChains() {
    ichains.clear();
    ichainClosed.clear();
    ichainOfEdge.assign(mesh->edges.size(), {-1, -1});
    ichainPosOfVertex.assign(mesh->vertices.size(), {});

    const std::vector<int> &featureEdges = graph->getFeatureEdges();
    if (featureEdges.empty()) return;

    const int nV = static_cast<int>(mesh->vertices.size());
    std::vector<std::vector<int>> atVertex(nV);
    for (const int e : featureEdges) {
        if (e < 0 || e >= static_cast<int>(mesh->edges.size())) continue;
        if (mesh->isBoundaryEdge[e]) continue;
        atVertex[mesh->edges[e][0]].push_back(e);
        atVertex[mesh->edges[e][1]].push_back(e);
    }
    auto otherEnd = [&](int e, int v) {
        return (mesh->edges[e][0] == v) ? mesh->edges[e][1] : mesh->edges[e][0];
    };
    auto isChainNode = [&](int v) {
        return atVertex[v].size() != 2 || mesh->isBoundaryVertex[v];
    };

    std::vector<char> used(mesh->edges.size(), 0);
    auto walkFrom = [&](int v, int e) {
        std::vector<int> verts{v};
        std::vector<int> chainEdges;
        int cur = v, ce = e;
        while (ce >= 0 && !used[ce]) {
            used[ce] = 1;
            chainEdges.push_back(ce);
            const int nxt = otherEnd(ce, cur);
            verts.push_back(nxt);
            if (isChainNode(nxt)) break;
            const int e2 = (atVertex[nxt][0] != ce) ? atVertex[nxt][0] : atVertex[nxt][1];
            cur = nxt;
            ce = e2;
        }
        if (chainEdges.empty()) return;
        const int idx = static_cast<int>(ichains.size());
        for (size_t i = 0; i < chainEdges.size(); ++i)
            ichainOfEdge[chainEdges[i]] = {idx, static_cast<int>(i)};
        for (size_t i = 0; i < verts.size(); ++i)
            ichainPosOfVertex[verts[i]].emplace_back(idx, static_cast<int>(i));
        ichainClosed.push_back(verts.front() == verts.back() ? 1 : 0);
        ichains.push_back(std::move(verts));
    };
    for (int v = 0; v < nV; ++v) {
        if (!isChainNode(v)) continue;
        for (const int e : atVertex[v]) if (!used[e]) walkFrom(v, e);
    }
    for (const int e : featureEdges) {
        if (e < 0 || e >= static_cast<int>(mesh->edges.size()) || used[e]) continue;
        if (mesh->isBoundaryEdge[e]) continue;
        walkFrom(mesh->edges[e][0], e);
    }
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
    splitsOfChain.assign(ichains.size(), {});

    // A node on an interface chain: the same job placeOnLoop does for dS, on a
    // run of vertices instead of a ring. The edge says which chain and which
    // step, and whether the chain runs the edge's way.
    auto placeOnChain = [&](int node, int e, double along) {
        if (e < 0 || e >= static_cast<int>(ichainOfEdge.size())) return false;
        const auto [c, i] = ichainOfEdge[e];
        if (c < 0) return false;
        const std::vector<int> &chain = ichains[c];
        if (i + 1 >= static_cast<int>(chain.size())) return false;
        const int va = mesh->edges[e][0];
        const double p = (chain[i] == va) ? i + along : i + (1.0 - along);
        splitsOfChain[c].push_back(Split{p, node});
        return true;
    };

    auto placeOnLoop = [&](int node, int boundaryEdge, double along) {
        if (boundaryEdge < 0) { ++report_.unplacedNodes; return; }
        // A ray that stopped on an interface exits through an interface edge,
        // not a boundary one, and its node belongs on that chain.
        if (!mesh->isBoundaryEdge[boundaryEdge]) {
            if (!placeOnChain(node, boundaryEdge, along)) ++report_.unplacedNodes;
            return;
        }
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
                    bool placed = false;
                    if (mn.vertex >= 0 && mn.vertex < static_cast<int>(loopOfVertex.size()) &&
                        loopOfVertex[mn.vertex].first >= 0) {
                        const auto &lv = loopOfVertex[mn.vertex];
                        splitsOfLoop[lv.first].push_back(
                            Split{static_cast<double>(lv.second), id});
                        placed = true;
                    }
                    // A ray launched from an interface sector starts on the
                    // interface, not on dS, and its corner splits every chain
                    // it sits on -- both of them at a junction.
                    if (mn.vertex >= 0 && mn.vertex < static_cast<int>(ichainPosOfVertex.size())) {
                        for (const auto &[c, i] : ichainPosOfVertex[mn.vertex]) {
                            splitsOfChain[c].push_back(Split{static_cast<double>(i), id});
                            placed = true;
                        }
                    }
                    if (!placed) ++report_.unplacedNodes;
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

    // The ends of every interface chain are nodes of the structure: a landing
    // on dS, or a junction where three or more materials meet. Each goes in as
    // a corner -- and a landing is filed on its boundary loop too, or the two
    // sides of dS there are one arc running straight past the interface and
    // the blocks either side of it have no side between them.
    //
    // "Each is a corner of every region around it" is a claim about the image
    // and it holds only because MotorcycleGraph::launchFeatureSectors() makes
    // it hold: a sector of q quarter turns at one of these gets q - 1 rays, so
    // after the launch every sector here spans exactly one quarter and every
    // one of them is a corner. Where that launch fails the claim fails with it
    // and the face comes out with a corner too many, which is a face the
    // decomposition refuses and the report counts -- the right way round, since
    // the fault is upstream.
    for (size_t c = 0; c < ichains.size(); ++c) {
        if (ichains[c].size() < 2) continue;
        const int ends[2] = {ichains[c].front(), ichains[c].back()};
        const double params[2] = {0.0, static_cast<double>(ichains[c].size() - 1)};
        for (int k = 0; k < 2; ++k) {
            const int v = ends[k];
            const int id = addNode(mesh->vertices[v], QuadLayout::NodeKind::BoundaryCorner);
            splitsOfChain[c].push_back(Split{params[k], id});
            // A closed chain has one vertex at both ends, so the split is
            // filed at both parameters and the ring is cut there once.
            if (v < static_cast<int>(loopOfVertex.size()) && loopOfVertex[v].first >= 0) {
                const auto &lv = loopOfVertex[v];
                splitsOfLoop[lv.first].push_back(Split{static_cast<double>(lv.second), id});
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
// buildInterfaceArcs()
//
// Each chain cut at the nodes on it, exactly as a boundary loop is. The arcs
// are not flagged onBoundary: an interface is interior to the model, and the
// face walk uses that flag to tell the unbounded side from a block. What makes
// an interface a side of the structure is that it is an arc at all.
// ---------------------------------------------------------------------------
void BlockLayout::buildInterfaceArcs() {
    for (size_t c = 0; c < ichains.size(); ++c) {
        const std::vector<int> &chain = ichains[c];
        const int n = static_cast<int>(chain.size());
        if (n < 2) continue;

        auto splits = splitsOfChain[c];
        std::sort(splits.begin(), splits.end(),
                  [](const Split &x, const Split &y) { return x.p < y.p; });
        splits.erase(std::unique(splits.begin(), splits.end(),
                                 [](const Split &x, const Split &y) {
                                     return x.node == y.node || std::fabs(x.p - y.p) < 1e-12;
                                 }),
                     splits.end());
        if (splits.size() < 2) continue;

        auto pointAt = [&](double p) {
            const int i = std::min(n - 2, static_cast<int>(std::floor(p)));
            if (i < 0) return mesh->vertices[chain[0]];
            const double t = p - i;
            return mesh->vertices[chain[i]] * (1.0 - t) + mesh->vertices[chain[i + 1]] * t;
        };

        for (size_t i = 0; i + 1 < splits.size(); ++i) {
            const double p0 = splits[i].p, p1 = splits[i + 1].p;
            std::vector<Point> pts;
            pts.push_back(pointAt(p0));
            for (int k = static_cast<int>(std::floor(p0)) + 1;
                 k <= static_cast<int>(std::ceil(p1)); ++k) {
                if (k <= p0 + 1e-12 || k >= p1 - 1e-12) continue;
                if (k >= 0 && k < n) pts.push_back(mesh->vertices[chain[k]]);
            }
            pts.push_back(pointAt(p1));
            const int id = addArc(std::move(pts), splits[i].node, splits[i + 1].node, -1, false);
            if (id >= 0) {
                if (static_cast<int>(arcIsInterface.size()) <= id)
                    arcIsInterface.resize(id + 1, 0);
                arcIsInterface[id] = 1;
            }
        }
    }
    arcIsInterface.resize(arcs.size(), 0);
}

// ---------------------------------------------------------------------------
// build()
// ---------------------------------------------------------------------------
void BlockLayout::build() {
    nodes.clear();
    arcs.clear();
    arcIsInterface.clear();
    report_ = Report{};

    buildBoundaryLoops();
    buildInterfaceChains();
    collectNodes();
    buildRayArcs();
    buildBoundaryArcs();
    buildInterfaceArcs();

    layout_.rebuild(nodes, arcs);
}

// ---------------------------------------------------------------------------
// blockDecompositionOf()
//
// The layout re-read as the shared representation. Nothing is recomputed: the
// arcs are already the sides, the nodes already the macrovertices, and the
// face walk already decided which sectors are corners. What this does is
// choose which faces qualify, put each arc down once in the direction its
// first claimant walks it, and hand the second claimant the same arc with a
// flip -- which is what makes a side shared full length rather than shared to
// a tolerance.
//
// The material is the one thing that is not already in the layout. It is read
// off the model by locating a point inside each block, and several points
// rather than one, because a block that straddles an interface is exactly the
// defect worth reporting on a multi-material model and a single sample cannot
// see it.
// ---------------------------------------------------------------------------
BlockDecomposition BlockLayout::blockDecompositionOf(const QuadLayout &layout, const Mesh &mesh,
                                                     DecompositionReport *out,
                                                     const std::string &source) {
    BlockDecomposition D;
    D.source = source;
    DecompositionReport rep;

    const std::vector<QuadLayout::Node> &lnodes = layout.getNodes();
    const std::vector<QuadLayout::Arc> &larcs = layout.getArcs();
    const std::vector<QuadLayout::Face> &lfaces = layout.getFaces();
    rep.faces = static_cast<int>(lfaces.size());

    std::vector<int> qualifying;
    qualifying.reserve(lfaces.size());
    for (size_t f = 0; f < lfaces.size(); ++f) {
        const QuadLayout::Face &fc = lfaces[f];
        rep.totalArea += fc.area;
        if (fc.corners != 4 || fc.sides.size() != 4) { ++rep.notFourSided; continue; }
        bool oneArcPerSide = true;
        for (const std::vector<int> &sd : fc.sides) if (sd.size() != 1) oneArcPerSide = false;
        if (!oneArcPerSide) { ++rep.multiArcSides; continue; }
        rep.coveredArea += fc.area;
        qualifying.push_back(static_cast<int>(f));
    }

    std::vector<int> vidx(lnodes.size(), -1);
    auto macroVertex = [&](int n) {
        if (n < 0 || n >= static_cast<int>(lnodes.size())) return -1;
        if (vidx[n] >= 0) return vidx[n];
        BlockDecomposition::MacroVertex mv;
        mv.p = lnodes[n].pos;
        mv.onBoundary = lnodes[n].kind == QuadLayout::NodeKind::BoundaryCorner ||
                        lnodes[n].kind == QuadLayout::NodeKind::BoundaryHit;
        vidx[n] = static_cast<int>(D.vertices.size());
        D.vertices.push_back(mv);
        return vidx[n];
    };

    std::vector<int> edgeOfArc(larcs.size(), -1);
    D.blocks.resize(qualifying.size());
    for (size_t bi = 0; bi < qualifying.size(); ++bi) {
        const QuadLayout::Face &fc = lfaces[qualifying[bi]];
        BlockDecomposition::Block &ob = D.blocks[bi];
        for (int s = 0; s < 4; ++s) {
            const int dart = fc.sides[s][0];
            const int a = QuadLayout::arcOfDart(dart);
            const bool forward = (dart & 1) == 0;
            const QuadLayout::Arc &arc = larcs[a];
            ob.corners[s] = macroVertex(forward ? arc.a : arc.b);

            int eid = edgeOfArc[a];
            if (eid < 0) {
                BlockDecomposition::MacroEdge oe;
                oe.boundary = arc.onBoundary;
                oe.points = arc.pts;
                if (!forward) std::reverse(oe.points.begin(), oe.points.end());
                oe.from = macroVertex(forward ? arc.a : arc.b);
                oe.to = macroVertex(forward ? arc.b : arc.a);
                oe.blockA = static_cast<int>(bi);
                oe.sideA = s;
                eid = static_cast<int>(D.edges.size());
                edgeOfArc[a] = eid;
                D.edges.push_back(std::move(oe));
            } else {
                // The twin dart runs the other way, so the polyline already
                // stored is this side read backwards; only the block needs
                // recording.
                D.edges[eid].blockB = static_cast<int>(bi);
                D.edges[eid].sideB = s;
            }
            ob.edges[s] = eid;
            ob.flip[s] = D.edges[eid].blockB == static_cast<int>(bi);
        }
    }

    // --- materials, read off the model ------------------------------------
    //
    // Nine samples of the parameter square plus a guaranteed interior point.
    // The majority is the block's material and a disagreement is the report's:
    // on a single-material model every sample agrees trivially and nothing
    // below costs anything.
    bool multiMaterial = false;
    for (size_t t = 1; t < mesh.triangleMatId.size(); ++t)
        if (mesh.triangleMatId[t] != mesh.triangleMatId[0]) { multiMaterial = true; break; }

    if (!mesh.triangles.empty()) {
        const TriangleGrid grid(mesh);
        std::unordered_set<int> seenMaterials;
        for (size_t b = 0; b < D.blocks.size(); ++b) {
            std::unordered_map<int, int> votes;
            Point q;
            if (D.interiorPoint(static_cast<int>(b), q)) {
                const int t = grid.locate(q);
                if (t >= 0) votes[mesh.triangleMatId[t]] += 2;   // the sure one counts double
            }
            if (multiMaterial) {
                for (int j = 1; j <= 3; ++j) {
                    for (int i = 1; i <= 3; ++i) {
                        const int t = grid.locate(
                            D.coonsPoint(static_cast<int>(b), i / 4.0, j / 4.0));
                        if (t >= 0) votes[mesh.triangleMatId[t]] += 1;
                    }
                }
            }
            int best = 0, bestVotes = -1;
            for (const auto &[m, v] : votes) if (v > bestVotes) { best = m; bestVotes = v; }
            D.blocks[b].material = best;
            if (votes.size() > 1) ++rep.straddlingBlocks;
            if (bestVotes >= 0) seenMaterials.insert(best);
        }
        rep.materials = static_cast<int>(seenMaterials.size());
    }

    for (BlockDecomposition::MacroEdge &oe : D.edges) {
        oe.matLeft = oe.blockA >= 0 ? D.blocks[oe.blockA].material : 0;
        oe.matRight = oe.blockB >= 0 ? D.blocks[oe.blockB].material : 0;
        // An interface is not a flag the layout carries; it is what an edge
        // with a different material on each side *is*. Reading it back this
        // way rather than plumbing a flag through keeps the two from ever
        // disagreeing.
        oe.interface = !oe.boundary && oe.blockA >= 0 && oe.blockB >= 0 &&
                       oe.matLeft != oe.matRight;
        if (oe.blockB < 0) ++rep.unmatchedSides;
        if (oe.from >= 0) {
            ++D.vertices[oe.from].valence;
            if (oe.boundary) D.vertices[oe.from].onBoundary = true;
            if (oe.interface) D.vertices[oe.from].onInterface = true;
        }
        if (oe.to >= 0) {
            ++D.vertices[oe.to].valence;
            if (oe.boundary) D.vertices[oe.to].onBoundary = true;
            if (oe.interface) D.vertices[oe.to].onInterface = true;
        }
    }

    rep.blocks = static_cast<int>(D.blocks.size());
    if (out) *out = rep;
    return D;
}
