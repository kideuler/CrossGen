#include "tracing/QuadLayout.hxx"

#include "MERIDIAN/Interfaces.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <unordered_map>

namespace {

// The point of a polyline at p = (segment - 1) + t.
Point pointAtParam(const std::vector<TracePoint> &path, double p) {
    const int n = static_cast<int>(path.size());
    if (n == 0) return Point{0.0, 0.0};
    p = std::min(std::max(p, 0.0), static_cast<double>(n - 1));
    const int i = std::min(n - 2, static_cast<int>(std::floor(p)));
    if (i < 0) return path[0].global_pos;
    const double t = p - i;
    return path[i].global_pos * (1.0 - t) + path[i + 1].global_pos * t;
}

// Where on a polyline a point lies: the closest place on it, as a parameter.
// Used for the two events that are known to sit on another separatrix -- a cut
// at a singularity and a second crossing -- where the segment is already known
// but the parameter along it is not.
double paramOfPoint(const std::vector<TracePoint> &path, const Point &q, int hintSeg) {
    const int n = static_cast<int>(path.size());
    double bestD = std::numeric_limits<double>::max();
    double bestP = 0.0;

    auto consider = [&](int i) {
        if (i < 1 || i >= n) return;
        const Point &a = path[i - 1].global_pos;
        const Point &b = path[i].global_pos;
        const Point d = b - a;
        const double dd = dotP(d, d);
        double t = (dd > 1e-30) ? dotP(q - a, d) / dd : 0.0;
        t = std::min(1.0, std::max(0.0, t));
        const Point proj = a + d * t;
        const double dist = normP(q - proj);
        if (dist < bestD) { bestD = dist; bestP = (i - 1) + t; }
    };

    // The hint is right except when a path was truncated under it, so it is
    // checked first and the whole path only if it turns out to be far off.
    consider(hintSeg);
    if (bestD > 1e-7) for (int i = 1; i < n; ++i) consider(i);
    return bestP;
}

inline double wrapTwoPi(double a) {
    a = std::fmod(a, 2.0 * M_PI);
    if (a < 0.0) a += 2.0 * M_PI;
    return a;
}

double polylineArea(const std::vector<Point> &poly) {
    double a = 0.0;
    for (size_t i = 0; i < poly.size(); ++i) {
        const Point &p = poly[i];
        const Point &q = poly[(i + 1) % poly.size()];
        a += cross2(p, q);
    }
    return 0.5 * a;
}

} // namespace

QuadLayout::QuadLayout(const SeparatrixTrace &t) : trace(&t) {
    mesh = &t.getTracer().getMesh();
    tol = 1e-6 * t.getTracer().averageEdgeLength();
}

QuadLayout::QuadLayout(double tolerance) : tol(tolerance) {}

void QuadLayout::rebuild(std::vector<Node> newNodes, std::vector<Arc> newArcs) {
    nodes = std::move(newNodes);
    arcs = std::move(newArcs);
    faces.clear();
    nodeHash.clear();
    report_ = Report{};
    finish();
}

// ---------------------------------------------------------------------------
// addNode()  --  with merging of coincident ones
//
// Three separatrices meeting at what is geometrically one point arrive as three
// separate events, and a separatrix cut on a port ray can land exactly where
// another one crosses it. Merging by position is what keeps those from becoming
// several nodes a hair apart with a degenerate arc between them.
// ---------------------------------------------------------------------------
int QuadLayout::addNode(const Point &pos, NodeKind kind, int source) {
    const double cell = std::max(tol, 1e-300);
    const long long gx = static_cast<long long>(std::floor(pos[0] / cell));
    const long long gy = static_cast<long long>(std::floor(pos[1] / cell));

    for (long long dx = -1; dx <= 1; ++dx) {
        for (long long dy = -1; dy <= 1; ++dy) {
            const long long key = (gx + dx) * 73856093LL ^ (gy + dy) * 19349663LL;
            for (const auto &kv : nodeHash) {
                if (kv.first != key) continue;
                if (normP(nodes[kv.second].pos - pos) <= tol) {
                    // The more structural kind wins: a crossing that lands on a
                    // singularity is the singularity.
                    if (static_cast<int>(kind) < static_cast<int>(nodes[kv.second].kind)) {
                        nodes[kv.second].kind = kind;
                        nodes[kv.second].source = source;
                    }
                    return kv.second;
                }
            }
        }
    }

    Node n;
    n.pos = pos;
    n.kind = kind;
    n.source = source;
    nodes.push_back(std::move(n));
    const long long key = gx * 73856093LL ^ gy * 19349663LL;
    nodeHash.emplace_back(key, static_cast<int>(nodes.size()) - 1);
    return static_cast<int>(nodes.size()) - 1;
}

// ---------------------------------------------------------------------------
// collectNodes()
// ---------------------------------------------------------------------------
void QuadLayout::collectNodes() {
    const auto &seps = trace->separatrices;
    splitsOfSeparatrix.assign(seps.size(), {});
    splitsOfLoop.assign(trace->boundaryLoops.size(), {});

    // Where each boundary vertex sits on which loop, for placing the nodes that
    // land on the boundary.
    std::unordered_map<int, std::pair<int, int>> loopPos;   // vertex -> (loop, index)
    for (size_t l = 0; l < trace->boundaryLoops.size(); ++l)
        for (size_t i = 0; i < trace->boundaryLoops[l].size(); ++i)
            loopPos[trace->boundaryLoops[l][i]] = {static_cast<int>(l), static_cast<int>(i)};

    // Singularities and boundary corners first, so that anything landing on one
    // of them merges into it rather than the other way round.
    std::vector<int> nodeOfSingularity(trace->singularities.size(), -1);
    for (size_t i = 0; i < trace->singularities.size(); ++i)
        nodeOfSingularity[i] =
            addNode(trace->singularities[i].coordinates, NodeKind::Singularity, static_cast<int>(i));

    std::vector<int> nodeOfCorner(trace->boundaryCorners.size(), -1);
    for (size_t i = 0; i < trace->boundaryCorners.size(); ++i) {
        const BoundaryCorner &c = trace->boundaryCorners[i];
        nodeOfCorner[i] = addNode(mesh->vertices[c.vertex], NodeKind::BoundaryCorner,
                                  static_cast<int>(i));
        auto it = loopPos.find(c.vertex);
        if (it != loopPos.end())
            splitsOfLoop[it->second.first].push_back(Split{static_cast<double>(it->second.second),
                                                           nodeOfCorner[i]});
    }

    // The nodes of the material interface network, each where it sits on its
    // branches and, for a landing, on its boundary loop -- a landing on a
    // straight stretch of dS is no corner of the boundary, but it is where
    // the boundary arcs of two different materials meet.
    const Interfaces *itf = trace->getInterfaces();
    std::vector<int> nodeOfInterface(trace->interfaceNodes.size(), -1);
    if (itf) {
        splitsOfBranch.assign(itf->branches().size(), {});
        std::vector<int> nodeOfNetworkNode(itf->nodes().size(), -1);
        for (size_t i = 0; i < trace->interfaceNodes.size(); ++i) {
            const InterfaceNodeInfo &in = trace->interfaceNodes[i];
            nodeOfInterface[i] = addNode(mesh->vertices[in.vertex], NodeKind::InterfaceNode,
                                         static_cast<int>(i));
            nodeOfNetworkNode[in.node] = nodeOfInterface[i];
            auto it = loopPos.find(in.vertex);
            if (in.onBoundary && it != loopPos.end())
                splitsOfLoop[it->second.first].push_back(
                    Split{static_cast<double>(it->second.second), nodeOfInterface[i]});
        }
        for (size_t b = 0; b < itf->branches().size(); ++b) {
            const Interfaces::Branch &br = itf->branches()[b];
            if (br.closed || br.verts.size() < 2) continue;
            for (const int end : {0, 1}) {
                const int nn = end ? br.node1 : br.node0;
                const int v = end ? br.verts.back() : br.verts.front();
                int nd = (nn >= 0) ? nodeOfNetworkNode[nn] : -1;
                // A branch ending on nothing is a tagging error the network
                // already reports; it still needs a node to end on.
                if (nd < 0) nd = addNode(mesh->vertices[v], NodeKind::Dangling, -1);
                splitsOfBranch[b].push_back(
                    Split{end ? static_cast<double>(br.verts.size() - 1) : 0.0, nd});
            }
        }
    }

    // The two ends of every separatrix.
    for (const auto &sep : seps) {
        if (sep.path.size() < 2) continue;

        const int start = (sep.originKind == SeparatrixOrigin::Singularity)
                              ? nodeOfSingularity[sep.origin_singularity_id]
                          : (sep.originKind == SeparatrixOrigin::InterfaceNode)
                              ? nodeOfInterface[sep.origin_singularity_id]
                              : nodeOfCorner[sep.origin_singularity_id];
        splitsOfSeparatrix[sep.id].push_back(Split{0.0, start});

        const double endP = static_cast<double>(sep.path.size() - 1);
        const Point &endPos = sep.path.back().global_pos;

        switch (sep.termination_reason) {
            case TerminationReason::EXIT_BOUNDARY: {
                const int nd = addNode(endPos, NodeKind::BoundaryHit, sep.id);
                splitsOfSeparatrix[sep.id].push_back(Split{endP, nd});

                // Place it on its loop: at a vertex, or along the edge it left
                // through.
                if (sep.endBoundaryVertex >= 0) {
                    auto it = loopPos.find(sep.endBoundaryVertex);
                    if (it != loopPos.end())
                        splitsOfLoop[it->second.first].push_back(
                            Split{static_cast<double>(it->second.second), nd});
                } else if (sep.endBoundaryEdge >= 0) {
                    const int va = mesh->edges[sep.endBoundaryEdge][0];
                    const int vb = mesh->edges[sep.endBoundaryEdge][1];
                    auto ia = loopPos.find(va);
                    auto ib = loopPos.find(vb);
                    if (ia != loopPos.end() && ib != loopPos.end() &&
                        ia->second.first == ib->second.first) {
                        const int l = ia->second.first;
                        const int n = static_cast<int>(trace->boundaryLoops[l].size());
                        // Which of the two orderings this edge has round the loop.
                        const bool aFirst = ((ia->second.second + 1) % n) == ib->second.second;
                        const int i0 = aFirst ? ia->second.second : ib->second.second;
                        const Point &p0 = mesh->vertices[trace->boundaryLoops[l][i0]];
                        const Point &p1 = mesh->vertices[trace->boundaryLoops[l][(i0 + 1) % n]];
                        const Point d = p1 - p0;
                        const double dd = dotP(d, d);
                        double t = (dd > 1e-30) ? dotP(endPos - p0, d) / dd : 0.0;
                        t = std::min(1.0, std::max(0.0, t));
                        splitsOfLoop[l].push_back(Split{i0 + t, nd});
                        // Put the node exactly where its own loop parameter
                        // says it is. The exit point and the parameter agree to
                        // rounding, and where they do not -- a separatrix
                        // leaving within a rounding error of a boundary vertex
                        // -- the two boundary arcs meeting at the node are
                        // ordered one way round the loop and drawn the other,
                        // and cross each other just beside it.
                        if (nodes[nd].kind == NodeKind::BoundaryHit)
                            nodes[nd].pos = p0 * (1.0 - t) + p1 * t;
                    }
                }
                break;
            }

            case TerminationReason::HETEROCLINIC: {
                // Both halves end at the same point, so this merges with its
                // partner into one node of valence two and the two arcs read as
                // the single curve they stand for.
                splitsOfSeparatrix[sep.id].push_back(
                    Split{endP, addNode(endPos, NodeKind::Heteroclinic, sep.id)});
                break;
            }

            case TerminationReason::CUT_AT_SINGULARITY:
            case TerminationReason::CROSSED_TWICE:
            case TerminationReason::TANGENTIAL: {
                const int nd = addNode(endPos, NodeKind::TJunction, sep.id);
                splitsOfSeparatrix[sep.id].push_back(Split{endP, nd});
                const int host = sep.endOnSeparatrix;
                if (host >= 0 && host < static_cast<int>(seps.size()) && seps[host].path.size() > 1)
                    splitsOfSeparatrix[host].push_back(
                        Split{paramOfPoint(seps[host].path, endPos, sep.endOnSegment), nd});
                break;
            }

            default: {
                const int nd = addNode(endPos, NodeKind::Dangling, sep.id);
                splitsOfSeparatrix[sep.id].push_back(Split{endP, nd});
                break;
            }
        }
    }

    // Every crossing of an interface, which splits the separatrix and the
    // interface branch both. A crossing recorded on a stretch of separatrix a
    // later truncation removed is no longer on the path, and is dropped.
    if (itf) {
        for (const InterfaceCrossing &x : trace->interfaceCrossings) {
            if (x.sep < 0 || x.sep >= static_cast<int>(seps.size())) continue;
            const auto &path = seps[x.sep].path;
            if (x.pathIndex <= 0 || x.pathIndex >= static_cast<int>(path.size())) continue;
            if (normP(path[x.pathIndex].global_pos - x.pos) > tol) continue;

            int b = -1;
            double p = 0.0;
            if (x.edge >= 0) {
                b = itf->branchOfEdge(x.edge);
                if (b < 0) continue;
                const auto &br = itf->branches()[b];
                for (size_t i = 0; i < br.edges.size(); ++i) {
                    if (br.edges[i] != x.edge) continue;
                    const Point &A = mesh->vertices[br.verts[i]];
                    const Point d = mesh->vertices[br.verts[i + 1]] - A;
                    const double dd = dotP(d, d);
                    const double t = dd > 0.0 ? std::min(1.0, std::max(0.0, dotP(x.pos - A, d) / dd)) : 0.0;
                    p = i + t;
                    break;
                }
            } else if (x.vertex >= 0) {
                b = itf->branchAt(x.vertex);   // -1 at a node, which is its own split
                if (b >= 0) {
                    const auto &br = itf->branches()[b];
                    for (size_t i = 0; i < br.verts.size(); ++i)
                        if (br.verts[i] == x.vertex) { p = static_cast<double>(i); break; }
                }
            }
            const int nd = addNode(x.pos, NodeKind::InterfaceHit, x.sep);
            splitsOfSeparatrix[x.sep].push_back(Split{static_cast<double>(x.pathIndex), nd});
            if (b >= 0) splitsOfBranch[b].push_back(Split{p, nd});
        }
    }

    // And every crossing, which splits both separatrices through it.
    for (const auto &x : trace->crossings) {
        if (x.sepA < 0 || x.sepB < 0) continue;
        const auto &A = seps[x.sepA];
        const auto &B = seps[x.sepB];
        if (A.path.size() < 2 || B.path.size() < 2) continue;
        const double pa = paramOfPoint(A.path, x.pos, x.segA);
        const double pb = paramOfPoint(B.path, x.pos, x.segB);
        // A crossing recorded before one of the two was truncated can end up
        // past the end of the shortened path; it is no longer on the layout.
        if (pa > static_cast<double>(A.path.size() - 1) - 1e-12) continue;
        if (pb > static_cast<double>(B.path.size() - 1) - 1e-12) continue;

        const int nd = addNode(x.pos, NodeKind::Crossing, -1);
        splitsOfSeparatrix[x.sepA].push_back(Split{pa, nd});
        splitsOfSeparatrix[x.sepB].push_back(Split{pb, nd});
    }
}

// ---------------------------------------------------------------------------
// Arcs
// ---------------------------------------------------------------------------
int QuadLayout::addArc(std::vector<Point> pts, int a, int b, int separatrix, bool onBoundary) {
    // Points closer together than the node tolerance carry no direction, and a
    // direction is exactly what the sort around a node needs.
    std::vector<Point> clean;
    clean.reserve(pts.size());
    for (const Point &p : pts)
        if (clean.empty() || normP(p - clean.back()) > tol * 0.5) clean.push_back(p);
    if (clean.size() < 2) return -1;
    clean.front() = nodes[a].pos;
    clean.back() = nodes[b].pos;
    if (normP(clean.front() - clean.back()) < tol && clean.size() < 3) return -1;

    Arc arc;
    arc.length = 0.0;
    for (size_t i = 1; i < clean.size(); ++i) arc.length += normP(clean[i] - clean[i - 1]);
    arc.pts = std::move(clean);
    arc.a = a;
    arc.b = b;
    arc.separatrix = separatrix;
    arc.onBoundary = onBoundary;
    arcs.push_back(std::move(arc));
    return static_cast<int>(arcs.size()) - 1;
}

void QuadLayout::buildSeparatrixArcs() {
    const auto &seps = trace->separatrices;
    for (const auto &sep : seps) {
        if (sep.path.size() < 2) continue;
        auto splits = splitsOfSeparatrix[sep.id];
        std::sort(splits.begin(), splits.end(),
                  [](const Split &x, const Split &y) { return x.p < y.p; });

        // Two splits at the same place are one node twice.
        splits.erase(std::unique(splits.begin(), splits.end(),
                                 [](const Split &x, const Split &y) {
                                     return x.node == y.node || std::fabs(x.p - y.p) < 1e-12;
                                 }),
                     splits.end());

        for (size_t i = 1; i < splits.size(); ++i) {
            const double p0 = splits[i - 1].p;
            const double p1 = splits[i].p;
            std::vector<Point> pts;
            pts.push_back(pointAtParam(sep.path, p0));
            for (int k = static_cast<int>(std::floor(p0)) + 1; k <= static_cast<int>(std::ceil(p1));
                 ++k) {
                if (k <= p0 + 1e-12 || k >= p1 - 1e-12) continue;
                pts.push_back(sep.path[k].global_pos);
            }
            pts.push_back(pointAtParam(sep.path, p1));
            addArc(std::move(pts), splits[i - 1].node, splits[i].node, sep.id, false);
        }
    }
}

void QuadLayout::buildBoundaryArcs() {
    for (size_t l = 0; l < trace->boundaryLoops.size(); ++l) {
        const auto &loop = trace->boundaryLoops[l];
        const int n = static_cast<int>(loop.size());
        if (n < 3) continue;

        auto splits = splitsOfLoop[l];
        if (splits.empty()) {
            // A loop nothing at all lands on -- an untouched hole. It still has
            // to be an arc or the material next to it is not enclosed, so it
            // gets one node and becomes a single closed arc.
            splits.push_back(Split{0.0, addNode(mesh->vertices[loop[0]], NodeKind::BoundaryHit, -1)});
        }
        std::sort(splits.begin(), splits.end(),
                  [](const Split &x, const Split &y) { return x.p < y.p; });
        splits.erase(std::unique(splits.begin(), splits.end(),
                                 [](const Split &x, const Split &y) {
                                     return x.node == y.node || std::fabs(x.p - y.p) < 1e-12;
                                 }),
                     splits.end());

        const int m = static_cast<int>(splits.size());
        for (int i = 0; i < m; ++i) {
            const Split &s0 = splits[i];
            const Split &s1 = splits[(i + 1) % m];
            const double p0 = s0.p;
            const double p1 = (i + 1 < m) ? s1.p : s1.p + n;  // the last arc wraps

            std::vector<Point> pts;
            pts.push_back(nodes[s0.node].pos);
            for (int k = static_cast<int>(std::floor(p0)) + 1; k <= static_cast<int>(std::floor(p1));
                 ++k) {
                if (k <= p0 + 1e-12 || k >= p1 - 1e-12) continue;
                pts.push_back(mesh->vertices[loop[k % n]]);
            }
            pts.push_back(nodes[s1.node].pos);
            addArc(std::move(pts), s0.node, s1.node, -1, true);
        }
    }
}

// ---------------------------------------------------------------------------
// buildInterfaceArcs()  --  the interface branches, cut at their nodes
//
// Exactly as a boundary loop is cut, with two differences: an interface arc
// has a component on both sides, so it is never the side a face walk leaves
// the model by (onBoundary stays false), and a closed branch -- an inclusion
// the network put no node on -- is cut only where separatrices cross it.
// ---------------------------------------------------------------------------
void QuadLayout::buildInterfaceArcs() {
    const Interfaces *itf = trace ? trace->getInterfaces() : nullptr;
    if (!itf) return;
    for (size_t b = 0; b < itf->branches().size(); ++b) {
        const Interfaces::Branch &br = itf->branches()[b];
        const int n = static_cast<int>(br.verts.size());
        if (n < 2) continue;
        auto splits = splitsOfBranch[b];
        if (br.closed && splits.empty())
            splits.push_back(Split{0.0, addNode(mesh->vertices[br.verts[0]], NodeKind::InterfaceHit, -1)});
        std::sort(splits.begin(), splits.end(),
                  [](const Split &x, const Split &y) { return x.p < y.p; });
        splits.erase(std::unique(splits.begin(), splits.end(),
                                 [](const Split &x, const Split &y) {
                                     return x.node == y.node && std::fabs(x.p - y.p) < 1e-12;
                                 }),
                     splits.end());
        const int m = static_cast<int>(splits.size());
        // A closed branch wraps: its last split joins its first one round the
        // loop, whose length in parameter is n - 1 (verts repeats its start).
        const int arcsHere = br.closed ? m : m - 1;
        for (int i = 0; i < arcsHere; ++i) {
            const Split &s0 = splits[i];
            const Split &s1 = splits[(i + 1) % m];
            const double p0 = s0.p;
            const double p1 = (i + 1 < m) ? s1.p : s1.p + (n - 1);
            std::vector<Point> pts;
            pts.push_back(nodes[s0.node].pos);
            for (int k = static_cast<int>(std::floor(p0)) + 1; k <= static_cast<int>(std::ceil(p1)); ++k) {
                if (k <= p0 + 1e-12 || k >= p1 - 1e-12) continue;
                pts.push_back(mesh->vertices[br.verts[k % (n - 1 > 0 ? n - 1 : 1)]]);
            }
            pts.push_back(nodes[s1.node].pos);
            const int a = addArc(std::move(pts), s0.node, s1.node, -1, false);
            if (a >= 0) arcs[a].onInterface = true;
        }
    }
}

// ---------------------------------------------------------------------------
// sortAroundNodes()
// ---------------------------------------------------------------------------
void QuadLayout::sortAroundNodes() {
    for (auto &n : nodes) { n.darts.clear(); n.angles.clear(); }

    for (int i = 0; i < static_cast<int>(arcs.size()); ++i) {
        const Arc &arc = arcs[i];
        if (arc.a < 0 || arc.b < 0) continue;
        nodes[arc.a].darts.push_back(2 * i);
        nodes[arc.a].angles.push_back(computeAngle(arc.pts[1] - arc.pts[0]));
        nodes[arc.b].darts.push_back(2 * i + 1);
        nodes[arc.b].angles.push_back(
            computeAngle(arc.pts[arc.pts.size() - 2] - arc.pts.back()));
    }

    for (auto &n : nodes) {
        std::vector<int> order(n.darts.size());
        for (size_t i = 0; i < order.size(); ++i) order[i] = static_cast<int>(i);
        std::sort(order.begin(), order.end(),
                  [&](int x, int y) { return n.angles[x] < n.angles[y]; });
        std::vector<int> d(order.size());
        std::vector<double> a(order.size());
        for (size_t i = 0; i < order.size(); ++i) {
            d[i] = n.darts[order[i]];
            a[i] = n.angles[order[i]];
        }
        n.darts = std::move(d);
        n.angles = std::move(a);
    }
}

// ---------------------------------------------------------------------------
// traceFaces()
//
// Having arrived at a node along a dart, the face continues along the dart one
// step clockwise from the reverse of the one just travelled. That is the
// standard walk and it produces every face counter-clockwise and the unbounded
// one clockwise, which is how the unbounded one is told apart and dropped.
// ---------------------------------------------------------------------------
void QuadLayout::traceFaces() {
    const int nDarts = 2 * static_cast<int>(arcs.size());
    if (nDarts == 0) return;

    std::unordered_map<int, int> slotOfDart;
    slotOfDart.reserve(nDarts * 2);
    for (const auto &n : nodes)
        for (size_t i = 0; i < n.darts.size(); ++i) slotOfDart[n.darts[i]] = static_cast<int>(i);

    auto originOf = [&](int dart) {
        const Arc &arc = arcs[dart >> 1];
        return (dart & 1) ? arc.b : arc.a;
    };
    auto nextDart = [&](int dart) {
        const int twin = dart ^ 1;
        const int at = originOf(twin);
        const Node &n = nodes[at];
        auto it = slotOfDart.find(twin);
        if (it == slotOfDart.end() || n.darts.empty()) return -1;
        const int k = static_cast<int>(n.darts.size());
        return n.darts[(it->second - 1 + k) % k];
    };
    auto pointsOf = [&](int dart) {
        const Arc &arc = arcs[dart >> 1];
        if (!(dart & 1)) return arc.pts;
        std::vector<Point> r(arc.pts.rbegin(), arc.pts.rend());
        return r;
    };

    std::vector<char> used(nDarts, 0);
    for (int d0 = 0; d0 < nDarts; ++d0) {
        if (used[d0]) continue;

        std::vector<int> cycle;
        int d = d0;
        bool ok = true;
        for (int guard = 0; guard <= nDarts; ++guard) {
            if (used[d]) { ok = (d == d0 && !cycle.empty()); break; }
            used[d] = 1;
            cycle.push_back(d);
            const int nd = nextDart(d);
            if (nd < 0) { ok = false; break; }
            d = nd;
            if (d == d0) break;
        }
        if (!ok || cycle.empty()) continue;

        // A boundary loop is stored with the material on its left, so a face
        // made of the model runs along it and never against it. That is what
        // separates the unbounded face -- and the void inside each hole, which
        // encloses positive area and would otherwise pass for a component --
        // from the components themselves.
        bool material = true;
        for (const int dd : cycle)
            if ((dd & 1) && arcs[dd >> 1].onBoundary) { material = false; break; }
        if (!material) { ++report_.unboundedCycles; continue; }

        std::vector<Point> poly;
        for (const int dd : cycle) {
            const auto pts = pointsOf(dd);
            for (size_t i = 0; i + 1 < pts.size(); ++i) poly.push_back(pts[i]);
        }
        const double area = polylineArea(poly);
        if (area <= 0.0) { ++report_.unboundedCycles; continue; }

        Face face;
        face.darts = cycle;
        face.area = area;
        face.nodes.resize(cycle.size());
        face.isCorner.assign(cycle.size(), 0);
        face.turn.assign(cycle.size(), 0.0);

        for (size_t i = 0; i < cycle.size(); ++i) {
            const int at = originOf(cycle[i]);
            face.nodes[i] = at;
            const int prev = cycle[(i + cycle.size() - 1) % cycle.size()];
            const auto pin = pointsOf(prev);
            const auto pout = pointsOf(cycle[i]);
            const double in = computeAngle(pin.back() - pin[pin.size() - 2]);
            const double out = computeAngle(pout[1] - pout[0]);
            face.turn[i] = wrap_pi(out - in);

            // Whether the face turns a corner here is a question about the
            // layout's combinatorics, not about the angle: a sector at a node
            // is either one quadrilateral's corner or a side running past a
            // T-junction, and which it is follows from what the node is. Judging
            // by angle instead miscounts every place the model itself is not
            // square -- a 67-degree corner of the geometry is still one corner
            // of one component, and so is a 17-degree spike.
            // What the node is, structurally, rather than what it was labelled
            // when it was created: partition simplification merges nodes, and a
            // label that was right before the merge need not survive it.
            const Node &nd = nodes[at];
            const bool boundary = [&] {
                for (const int d : nd.darts) if (arcs[d >> 1].onBoundary) return true;
                return false;
            }();
            bool corner = true;
            if (nd.kind == NodeKind::BoundaryCorner || nd.kind == NodeKind::InterfaceNode) {
                // A corner of the model is a corner of whatever component sits
                // in it, at any angle: a convex one has only its two boundary
                // arcs and no more, and is still where a component turns.
                corner = true;
            } else if (nd.darts.size() <= 2) {
                // Two arcs meeting anywhere else is one curve carrying on,
                // which is what a heteroclinic joint is; one arc is a loose end
                // the face runs out to and back from.
                corner = false;
            } else if (nd.darts.size() == 3 && !boundary && nd.kind != NodeKind::Singularity) {
                // Three arcs meeting away from the boundary, at a node the
                // field has no irregularity at, is a separatrix stopped against
                // a side. The side it stopped against runs past unbroken, so of
                // the three sectors the widest is not a corner and the other two
                // are.
                int widest = 0;
                double best = -1.0;
                for (size_t k = 0; k < nd.darts.size(); ++k) {
                    const double w = wrapTwoPi(nd.angles[(k + 1) % nd.darts.size()] - nd.angles[k]);
                    if (w > best) { best = w; widest = static_cast<int>(k); }
                }
                auto it = slotOfDart.find(cycle[i]);
                corner = (it == slotOfDart.end()) || (it->second != widest);
            }
            if (corner) {
                face.isCorner[i] = 1;
                ++face.corners;
            }
        }

        // Cut the cycle into sides at the corners, starting at one of them so
        // that a side is never split across the seam.
        int first = -1;
        for (size_t i = 0; i < cycle.size(); ++i) if (face.isCorner[i]) { first = static_cast<int>(i); break; }
        if (first < 0) {
            face.sides.push_back(face.darts);   // a component with no corner at all
        } else {
            std::vector<int> side;
            for (size_t k = 0; k < cycle.size(); ++k) {
                const size_t i = (first + k) % cycle.size();
                if (face.isCorner[i] && !side.empty()) { face.sides.push_back(side); side.clear(); }
                side.push_back(cycle[i]);
            }
            if (!side.empty()) face.sides.push_back(side);
        }
        faces.push_back(std::move(face));
    }
}

// ---------------------------------------------------------------------------
// checkEmbedding()  --  are the arcs actually disjoint?
//
// Everything downstream assumes they are: the face walk reads the plane graph
// off the cyclic order at each node, and if two arcs cross where no node is,
// the walk happily produces two faces that overlap. Their areas then add up to
// more than the model, which is the cheap symptom; this is the direct test.
// ---------------------------------------------------------------------------
void QuadLayout::checkEmbedding() {
    struct Box { Point lo, hi; };
    std::vector<Box> boxes(arcs.size());
    for (size_t i = 0; i < arcs.size(); ++i) {
        Box b{arcs[i].pts[0], arcs[i].pts[0]};
        for (const Point &p : arcs[i].pts) {
            b.lo[0] = std::min(b.lo[0], p[0]); b.lo[1] = std::min(b.lo[1], p[1]);
            b.hi[0] = std::max(b.hi[0], p[0]); b.hi[1] = std::max(b.hi[1], p[1]);
        }
        boxes[i] = b;
    }

    for (size_t i = 0; i < arcs.size(); ++i) {
        for (size_t j = i + 1; j < arcs.size(); ++j) {
            if (boxes[i].hi[0] < boxes[j].lo[0] - tol || boxes[j].hi[0] < boxes[i].lo[0] - tol ||
                boxes[i].hi[1] < boxes[j].lo[1] - tol || boxes[j].hi[1] < boxes[i].lo[1] - tol)
                continue;
            // Two arcs that share a node meet there, and the segments on either
            // side of that meeting are collinear to within the rounding of the
            // node's own position, so a segment-segment test on them reports a
            // crossing a hair from the node. Those two segments are skipped
            // against each other; every other pair still counts, so an arc
            // genuinely doubling back over its neighbour is not missed.
            const bool shareIaJa = arcs[i].a == arcs[j].a, shareIaJb = arcs[i].a == arcs[j].b;
            const bool shareIbJa = arcs[i].b == arcs[j].a, shareIbJb = arcs[i].b == arcs[j].b;
            const size_t ni = arcs[i].pts.size(), nj = arcs[j].pts.size();
            for (size_t x = 1; x < ni; ++x) {
                const bool xAtA = (x == 1), xAtB = (x == ni - 1);
                for (size_t y = 1; y < nj; ++y) {
                    const bool yAtA = (y == 1), yAtB = (y == nj - 1);
                    if ((shareIaJa && xAtA && yAtA) || (shareIaJb && xAtA && yAtB) ||
                        (shareIbJa && xAtB && yAtA) || (shareIbJb && xAtB && yAtB))
                        continue;

                    const Point &p0 = arcs[i].pts[x - 1], &p1 = arcs[i].pts[x];
                    const Point &q0 = arcs[j].pts[y - 1], &q1 = arcs[j].pts[y];
                    const Point r = p1 - p0, s = q1 - q0;
                    const double den = cross2(r, s);
                    // Parallel to within rounding: relative to the two
                    // lengths, not absolute. Two collinear pieces of one
                    // straight edge either side of a node have a cross
                    // product of ~1e-19, which an absolute 1e-30 lets through,
                    // and the parameters that come out of it are noise that
                    // can land inside both -- a crossing that is not there.
                    if (std::fabs(den) <= 1e-12 * normP(r) * normP(s)) continue;
                    const Point w = q0 - p0;
                    const double t = cross2(w, s) / den;
                    const double u = cross2(w, r) / den;
                    if (t <= 0.0 || t >= 1.0 || u <= 0.0 || u >= 1.0) continue;
                    ++report_.arcCrossings;
                    x = ni;
                    break;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// build()
// ---------------------------------------------------------------------------
void QuadLayout::build() {
    nodes.clear();
    arcs.clear();
    faces.clear();
    nodeHash.clear();
    report_ = Report{};

    collectNodes();
    buildSeparatrixArcs();
    buildBoundaryArcs();
    buildInterfaceArcs();
    finish();
}

void QuadLayout::finish() {
    sortAroundNodes();
    traceFaces();
    checkEmbedding();

    report_.nodes = static_cast<int>(nodes.size());
    report_.arcs = static_cast<int>(arcs.size());
    report_.faces = static_cast<int>(faces.size());
    for (const auto &n : nodes) {
        if (n.kind == NodeKind::Singularity) ++report_.singularities;
        if (n.kind == NodeKind::Dangling) ++report_.danglingEnds;
        // Counted by what the node is: three arcs meeting in the interior at
        // something that is not a singularity is a T-junction whatever it was
        // called when it was made.
        if (n.darts.size() != 3 || n.kind == NodeKind::Singularity ||
            n.kind == NodeKind::InterfaceNode)
            continue;
        bool boundary = false;
        for (const int d : n.darts) if (arcs[d >> 1].onBoundary) boundary = true;
        if (!boundary) ++report_.tJunctions;
    }
    // Whether the layout is a valid T-layout is combinatorial, and reported by
    // the corner counts. Whether it is a *good* one is separate: a crossing
    // that is not square, or a separatrix meeting the boundary well off normal,
    // is a component the mesher will not like even though it has four sides.
    // Sec. 4 is what deals with those; this only counts them.
    for (const auto &n : nodes) {
        if (n.darts.size() < 3) continue;
        const double want = 2.0 * M_PI / static_cast<double>(n.darts.size());
        double worst = 0.0;
        for (size_t k = 0; k < n.darts.size(); ++k) {
            const double w = wrapTwoPi(n.angles[(k + 1) % n.darts.size()] - n.angles[k]);
            worst = std::max(worst, std::fabs(w - want));
        }
        bool boundaryHere = false;
        for (const int d : n.darts) if (arcs[d >> 1].onBoundary) boundaryHere = true;
        if (n.darts.size() == 3 && !boundaryHere && n.kind != NodeKind::Singularity &&
            n.kind != NodeKind::InterfaceNode) {
            // Two right angles and a straight side, not three equal sectors.
            worst = 0.0;
            std::vector<double> w(n.darts.size());
            for (size_t k = 0; k < n.darts.size(); ++k)
                w[k] = wrapTwoPi(n.angles[(k + 1) % n.darts.size()] - n.angles[k]);
            std::sort(w.begin(), w.end());
            worst = std::max(std::fabs(w[0] - M_PI_2), std::fabs(w[1] - M_PI_2));
        }
        if (worst > M_PI_4) ++report_.skewNodes;
    }

    report_.smallestArea = std::numeric_limits<double>::max();
    for (const auto &f : faces) {
        if (f.corners == 4) ++report_.quadFaces; else ++report_.badFaces;
        report_.smallestArea = std::min(report_.smallestArea, f.area);
        report_.largestArea = std::max(report_.largestArea, f.area);
        report_.totalArea += f.area;
    }
    if (faces.empty()) report_.smallestArea = 0.0;

    if (mesh)
        for (const auto &t : mesh->triangles)
            report_.meshArea += 0.5 * std::fabs(cross2(mesh->vertices[t[1]] - mesh->vertices[t[0]],
                                                       mesh->vertices[t[2]] - mesh->vertices[t[0]]));
}

// ---------------------------------------------------------------------------
// VTK output
// ---------------------------------------------------------------------------
bool QuadLayout::writeArcsVTU(const std::string &filename) const {
    std::ofstream out(filename);
    if (!out) return false;

    size_t nPoints = 0;
    for (const auto &a : arcs) nPoints += a.pts.size();

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << arcs.size()
        << "\">\n";

    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto &a : arcs)
        for (const auto &p : a.pts) out << "          " << p[0] << " " << p[1] << " 0.0\n";
    out << "        </DataArray>\n";
    out << "      </Points>\n";

    out << "      <CellData Scalars=\"separatrix\">\n";
    out << "        <DataArray type=\"Int32\" Name=\"separatrix\" format=\"ascii\">\n";
    for (const auto &a : arcs) out << "          " << a.separatrix << "\n";
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"onBoundary\" format=\"ascii\">\n";
    for (const auto &a : arcs) out << "          " << (a.onBoundary ? 1 : 0) << "\n";
    out << "        </DataArray>\n";
    out << "      </CellData>\n";

    out << "      <Cells>\n";
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    size_t base = 0;
    for (const auto &a : arcs) {
        out << "         ";
        for (size_t i = 0; i < a.pts.size(); ++i) out << " " << base + i;
        out << "\n";
        base += a.pts.size();
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    size_t off = 0;
    for (const auto &a : arcs) { off += a.pts.size(); out << "          " << off << "\n"; }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < arcs.size(); ++i) out << "          4\n";  // VTK_POLY_LINE
    out << "        </DataArray>\n";
    out << "      </Cells>\n";

    out << "    </Piece>\n";
    out << "  </UnstructuredGrid>\n";
    out << "</VTKFile>\n";
    return true;
}

bool QuadLayout::writeFacesVTU(const std::string &filename) const {
    std::ofstream out(filename);
    if (!out) return false;

    // Each face as its own ring of points, so the colouring is per face.
    std::vector<std::vector<Point>> rings;
    rings.reserve(faces.size());
    for (const auto &f : faces) {
        std::vector<Point> ring;
        for (const int d : f.darts) {
            const Arc &arc = arcs[d >> 1];
            if (d & 1)
                for (size_t i = arc.pts.size() - 1; i > 0; --i) ring.push_back(arc.pts[i]);
            else
                for (size_t i = 0; i + 1 < arc.pts.size(); ++i) ring.push_back(arc.pts[i]);
        }
        rings.push_back(std::move(ring));
    }

    size_t nPoints = 0;
    for (const auto &r : rings) nPoints += r.size();

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << rings.size()
        << "\">\n";

    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto &r : rings)
        for (const auto &p : r) out << "          " << p[0] << " " << p[1] << " 0.0\n";
    out << "        </DataArray>\n";
    out << "      </Points>\n";

    out << "      <CellData Scalars=\"corners\">\n";
    out << "        <DataArray type=\"Int32\" Name=\"corners\" format=\"ascii\">\n";
    for (const auto &f : faces) out << "          " << f.corners << "\n";
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"sides\" format=\"ascii\">\n";
    for (const auto &f : faces) out << "          " << f.darts.size() << "\n";
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Float64\" Name=\"area\" format=\"ascii\">\n";
    for (const auto &f : faces) out << "          " << f.area << "\n";
    out << "        </DataArray>\n";
    out << "      </CellData>\n";

    out << "      <Cells>\n";
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    size_t base = 0;
    for (const auto &r : rings) {
        out << "         ";
        for (size_t i = 0; i < r.size(); ++i) out << " " << base + i;
        out << "\n";
        base += r.size();
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    size_t off = 0;
    for (const auto &r : rings) { off += r.size(); out << "          " << off << "\n"; }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < rings.size(); ++i) out << "          7\n";  // VTK_POLYGON
    out << "        </DataArray>\n";
    out << "      </Cells>\n";

    out << "    </Piece>\n";
    out << "  </UnstructuredGrid>\n";
    out << "</VTKFile>\n";
    return true;
}

bool QuadLayout::writeNodesVTU(const std::string &filename) const {
    std::ofstream out(filename);
    if (!out) return false;

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << nodes.size() << "\" NumberOfCells=\"" << nodes.size()
        << "\">\n";

    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto &n : nodes) out << "          " << n.pos[0] << " " << n.pos[1] << " 0.0\n";
    out << "        </DataArray>\n";
    out << "      </Points>\n";

    out << "      <PointData Scalars=\"kind\">\n";
    out << "        <DataArray type=\"Int32\" Name=\"kind\" format=\"ascii\">\n";
    for (const auto &n : nodes) out << "          " << static_cast<int>(n.kind) << "\n";
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"valence\" format=\"ascii\">\n";
    for (const auto &n : nodes) out << "          " << n.darts.size() << "\n";
    out << "        </DataArray>\n";
    out << "      </PointData>\n";

    out << "      <Cells>\n";
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (size_t i = 0; i < nodes.size(); ++i) out << "          " << i << "\n";
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (size_t i = 1; i <= nodes.size(); ++i) out << "          " << i << "\n";
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < nodes.size(); ++i) out << "          1\n";  // VTK_VERTEX
    out << "        </DataArray>\n";
    out << "      </Cells>\n";

    out << "    </Piece>\n";
    out << "  </UnstructuredGrid>\n";
    out << "</VTKFile>\n";
    return true;
}
