// Pillow.cxx -- see Pillow.hxx.
#include "Pillow.hxx"

#include <algorithm>
#include <cmath>
#include <deque>
#include <map>
#include <iomanip>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace mesh {

namespace {

constexpr double kDeg = M_PI / 180.0;

// Counter-clockwise angle from a to b, in [0, 2 pi).
double ccwAngle(const Point &a, const Point &b) {
    double t = std::atan2(cross2(a, b), dotP(a, b));
    if (t < 0.0) t += 2.0 * M_PI;
    return t;
}

int cornerOf(const Quad &q, int v) {
    for (int k = 0; k < 4; ++k)
        if (q[k] == v) return k;
    return -1;
}

double signedAreaOf(const std::vector<Point> &V, const Quad &q) {
    double a = 0.0;
    for (int k = 0; k < 4; ++k) a += cross2(V[q[k]], V[q[(k + 1) % 4]]);
    return 0.5 * a;
}

// Corners of q whose two sides turn the wrong way or not at all: the corner
// determinants a smoother's untangler would have to fix.
int foldedCorners(const std::vector<Point> &V, const Quad &q) {
    int n = 0;
    for (int k = 0; k < 4; ++k) {
        const Point &p = V[q[k]];
        if (cross2(V[q[(k + 1) % 4]] - p, V[q[(k + 3) % 4]] - p) <= 0.0) ++n;
    }
    return n;
}

}  // namespace

Pillow::Pillow(QuadMesh &m) : Pillow(m, Options()) {}

Pillow::Pillow(QuadMesh &m, const Options &opts) : mesh(m), options(opts) {}

std::vector<Pillow::Defect> Pillow::findDefects(const QuadMesh &m, double flatAngle) {
    std::vector<Defect> out;
    const int nQ = static_cast<int>(m.quads.size());
    if (static_cast<int>(m.quadEdges.size()) != nQ) return out;
    const double limit = flatAngle * kDeg;
    for (int q = 0; q < nQ; ++q) {
        const Quad &Q = m.quads[q];
        for (int c = 0; c < 4; ++c) {
            const int eIn = m.quadEdges[q][(c + 3) % 4];
            const int eOut = m.quadEdges[q][c];
            if (eIn < 0 || eOut < 0) continue;
            if (!m.isFeatureEdge[eIn] || !m.isFeatureEdge[eOut]) continue;
            // The quad is counter-clockwise, so its interior at corner c is
            // swept counter-clockwise from the side to the next corner round to
            // the side to the previous one.
            const Point &p = m.vertices[Q[c]];
            const double a = ccwAngle(m.vertices[Q[(c + 1) % 4]] - p,
                                      m.vertices[Q[(c + 3) % 4]] - p);
            if (a > limit) out.push_back(Defect{q, c, Q[c], a / kDeg});
        }
    }
    return out;
}

// Turning round u from the feature edge (from, u), through the elements on the
// side `start` is on, to the next feature edge. Each element is left through
// its other side at u, and the turn stops on the first of those that is a
// feature. Which way round that is depends only on which corner of `start`
// `from` sits at, so the same code turns either way.
bool Pillow::turn(int u, int from, int start, Fan &fan) const {
    fan = Fan();
    int R = start, w = from;
    bool ccw = false;
    for (int guard = 0; guard < 64; ++guard) {
        const Quad &Q = mesh.quads[R];
        const int i = cornerOf(Q, u);
        if (i < 0) return false;
        const int prev = Q[(i + 3) % 4], next = Q[(i + 1) % 4];
        int x = -1, e = -1;
        if (next == w)      { x = prev; e = mesh.quadEdges[R][(i + 3) % 4]; }
        else if (prev == w) { x = next; e = mesh.quadEdges[R][i]; }
        else return false;
        if (fan.quads.empty()) ccw = (next == w);
        fan.quads.push_back(R);
        fan.spokes.push_back(x);
        if (e < 0) return false;
        if (mesh.isFeatureEdge[e]) { fan.next = x; break; }
        const std::array<int, 2> &eq = mesh.edgeQuads[e];
        const int R2 = (eq[0] == R) ? eq[1] : eq[0];
        if (R2 < 0 || R2 == start) return false;
        R = R2;
        w = x;
    }
    if (fan.next < 0) return false;
    const Point &p = mesh.vertices[u];
    const Point din = mesh.vertices[from] - p, dout = mesh.vertices[fan.next] - p;
    fan.angle = ccw ? ccwAngle(din, dout) : ccwAngle(dout, din);
    return true;
}

bool Pillow::walk(int v, int first, int startQuad, std::vector<PathNode> &nodes,
                  std::vector<int> &sideQuads, bool &closed) const {
    closed = false;
    int prev = v, u = first, inQuad = startQuad;
    sideQuads.push_back(startQuad);
    const int limit = static_cast<int>(mesh.vertices.size());
    for (int guard = 0; guard <= limit; ++guard) {
        if (u == v) { closed = true; return true; }
        Fan fan;
        if (!turn(u, prev, inQuad, fan) || fan.next == prev) {
            refusal = "the elements round node " + std::to_string(u) + " could not be walked";
            return false;
        }
        const double a = fan.angle / kDeg;
        const int k = static_cast<int>(fan.quads.size());
        // Two elements are right for a fan of 135 to 225 degrees, so a layer
        // through such a node leaves the node regular and gives its copy the
        // node's old elements. One element past flatAngle is the defect
        // itself, found again further along.
        const bool pass = std::fabs(a - 180.0) < options.passAngle ||
                          (k == 1 && a > options.flatAngle);
        PathNode n;
        n.v = u;
        if (pass) {
            n.fan = fan.quads;
            nodes.push_back(n);
            prev = u;
            u = fan.next;
            inQuad = fan.quads.back();
            sideQuads.push_back(inQuad);
            continue;
        }
        // A corner. The layer's copy takes over the first element only and
        // goes onto the edge that element leaves u by.
        n.end = true;
        n.fan = {fan.quads.front()};
        n.exitTo = fan.spokes.front();
        nodes.push_back(n);
        return true;
    }
    refusal = "the walk along the feature did not end";
    return false;
}

Point Pillow::onEdge(int a, int b, double s, bool feature) const {
    const Point &pa = mesh.vertices[a], &pb = mesh.vertices[b];
    const Point chord = pa + (pb - pa) * s;
    if (!feature) return chord;

    // A new node on a feature goes onto the curve the feature's sliding nodes
    // are bound to, between the parameters of the edge's two ends, so the
    // curve QuadMesh::buildTopology fits through it afterwards is the one it
    // was cut from. A run's fixed ends are not bound but are the first and last
    // vertices of its chain, at the two ends of the curve's domain.
    const int ci = mesh.isOnCurve(a) ? mesh.curveOf[a] : (mesh.isOnCurve(b) ? mesh.curveOf[b] : -1);
    if (ci < 0) return chord;
    const QuadMesh::FeatureCurve &fc = mesh.featureCurves[ci];
    const double lo = fc.curve.domainBegin(), hi = fc.curve.domainEnd();
    auto bound = [&](int x, double &u) {
        if (mesh.curveOf[x] == ci) { u = mesh.curveParam[x]; return true; }
        return false;
    };
    double ua = 0.0, ub = 0.0;
    const bool ha = bound(a, ua), hb = bound(b, ub);
    auto endOf = [&](int x, double other, double &u) {
        const bool front = !fc.chain.empty() && fc.chain.front() == x;
        const bool back = !fc.chain.empty() && fc.chain.back() == x;
        if (front && back) u = (other - lo < hi - other) ? lo : hi;   // a closed run's anchor
        else if (front) u = lo;
        else if (back) u = hi;
        else return false;
        return true;
    };
    if (!ha && !endOf(a, ub, ua)) return chord;
    if (!hb && !endOf(b, ua, ub)) return chord;
    const Point c = fc.curve.evaluate(ua + s * (ub - ua));
    // A curve point far off the chord means the parameters were not the
    // edge's after all. The chord is then the safer answer.
    if (!(normP(c - chord) <= 0.25 * normP(pb - pa))) return chord;
    return c;
}

bool Pillow::insertLayer(const Defect &d) {
    const int nV0 = static_cast<int>(mesh.vertices.size());
    const int nQ0 = static_cast<int>(mesh.quads.size());
    const Quad &D = mesh.quads[d.quad];
    const int v = d.vertex;
    const int fwd = D[(d.corner + 1) % 4], bwd = D[(d.corner + 3) % 4];

    // ---- the path: from the defect both ways along the feature ------------
    std::vector<PathNode> fnodes, bnodes;
    std::vector<int> fside, bside;
    bool fclosed = false, bclosed = false;
    if (!walk(v, fwd, d.quad, fnodes, fside, fclosed)) return false;

    PathNode centre;
    centre.v = v;
    centre.fan = {d.quad};
    std::vector<PathNode> path;
    std::vector<int> side;        // side[i]: the element of edge (path[i], path[i+1]) being pillowed
    bool closed = false;
    if (fclosed) {
        // Round a closed loop and back: the last edge walked has to be the
        // defect's other side, or the walk went somewhere it should not have.
        if (fnodes.empty() || fnodes.back().v != bwd || fside.back() != d.quad) {
            refusal = "the walk came back round without closing on the corner's other side";
            return false;
        }
        closed = true;
        path.push_back(centre);
        path.insert(path.end(), fnodes.begin(), fnodes.end());
        side = fside;
    } else {
        if (!walk(v, bwd, d.quad, bnodes, bside, bclosed)) return false;
        if (bclosed) {
            refusal = "one walk closed the loop and the other did not";
            return false;
        }
        for (int i = static_cast<int>(bnodes.size()) - 1; i >= 0; --i) path.push_back(bnodes[i]);
        path.push_back(centre);
        path.insert(path.end(), fnodes.begin(), fnodes.end());
        for (int i = static_cast<int>(bside.size()) - 1; i >= 0; --i) side.push_back(bside[i]);
        side.insert(side.end(), fside.begin(), fside.end());
    }
    const int n = static_cast<int>(path.size());
    const int nE = closed ? n : n - 1;
    if (static_cast<int>(side.size()) != nE || n < 2) {
        refusal = "the path and its elements do not match";
        return false;
    }

    // A path that touches itself, or an end whose exit edge leads back onto
    // the path, would have one node copied twice or an exit split crossing
    // the layer. Both happen only on a region whose boundary pinches, and the
    // defect is left alone there.
    std::unordered_set<int> onPath;
    for (const PathNode &p : path)
        if (!onPath.insert(p.v).second) {
            refusal = "the path touches itself at node " + std::to_string(p.v);
            return false;
        }
    for (const PathNode &p : path)
        if (p.end && (p.exitTo < 0 || onPath.count(p.exitTo))) {
            refusal = "the exit at node " + std::to_string(p.v) + " leads back onto the path";
            return false;
        }

    // Each path edge's direction in its pillowed element, read before any
    // corner is re-pointed: the new quad on it is wound the same way.
    std::vector<char> along(nE, 0);
    for (int e = 0; e < nE; ++e) {
        const int a = path[e].v, b = path[(e + 1) % n].v;
        const Quad &S = mesh.quads[side[e]];
        const int ka = cornerOf(S, a);
        if (ka >= 0 && S[(ka + 1) % 4] == b) along[e] = 1;
        else if (ka < 0 || S[(ka + 3) % 4] != b) {
            refusal = "a path edge is not a side of its element";
            return false;
        }
    }

    // Feature edges by vertex pair, for the new nodes that land on one.
    std::unordered_map<MeshEdgeKey, int, MeshEdgeKeyHash> edgeOf;
    edgeOf.reserve(mesh.edges.size());
    for (int e = 0; e < static_cast<int>(mesh.edges.size()); ++e)
        edgeOf.emplace(MeshEdgeKey(mesh.edges[e][0], mesh.edges[e][1]), e);
    auto isFeature = [&](int a, int b) {
        if (a >= nV0 || b >= nV0) return false;
        const auto it = edgeOf.find(MeshEdgeKey(a, b));
        return it != edgeOf.end() && mesh.isFeatureEdge[it->second];
    };

    // Where each copy goes when the layer is at full depth: the end of the
    // vector from the node to the mean centroid of the elements it takes over.
    std::vector<Point> toward(n);
    for (int i = 0; i < n; ++i) {
        if (path[i].end) continue;
        Point c{0.0, 0.0};
        for (int q : path[i].fan) c = c + mesh.centroid(q);
        toward[i] = c / static_cast<double>(path[i].fan.size()) - mesh.vertices[path[i].v];
    }

    // Elements Stage 10 already left folded, and below, the pieces a chord
    // cuts them into. The layer is not what folded them, so they are not held
    // against it -- geom028's chord crosses a block corner Stage 10 turned
    // over, and halving that element doubles its folded corners.
    std::vector<char> foldedAlready(nQ0, 0);
    for (int q = 0; q < nQ0; ++q)
        foldedAlready[q] = foldedCorners(mesh.vertices, mesh.quads[q]) > 0;

    // ---- build it, at the full depth or a thinner one ---------------------
    // A thinner layer starts nearer the mesh it was cut from, so it is tried
    // when the full one would leave an element turned over that was not.
    for (double f : {1.0, 0.5, 0.25}) {
        std::vector<Point> V = mesh.vertices;
        std::vector<Quad> Q = mesh.quads;
        std::vector<int> M = mesh.quadMatId;
        std::vector<char> touched(nQ0, 0);
        std::vector<char> wasFolded = foldedAlready;

        std::vector<int> copy(n);
        for (int i = 0; i < n; ++i) {
            const int u = path[i].v;
            copy[i] = static_cast<int>(V.size());
            if (path[i].end)
                V.push_back(onEdge(u, path[i].exitTo, f * options.exitFraction,
                                   isFeature(u, path[i].exitTo)));
            else
                V.push_back(mesh.vertices[u] + toward[i] * (f * options.depth));
        }

        bool ok = true;
        for (int i = 0; i < n && ok; ++i) {
            for (int R : path[i].fan) {
                const int k = cornerOf(Q[R], path[i].v);
                if (k < 0) { ok = false; break; }
                Q[R][k] = copy[i];
                touched[R] = 1;
            }
        }
        if (!ok) {
            refusal = "an element in a fan does not have the node it was found at";
            return false;
        }

        for (int e = 0; e < nE; ++e) {
            const int i = e, j = (e + 1) % n;
            const int a = path[i].v, b = path[j].v, ca = copy[i], cb = copy[j];
            Q.push_back(along[e] ? Quad{a, b, cb, ca} : Quad{b, a, ca, cb});
            M.push_back(M[side[e]]);
            wasFolded.push_back(0);
        }

        // ---- carry the chord on past the two ends -------------------------
        // Every edge still used whole by an element on one side while the
        // other side already runs through a node in its middle. Splitting the
        // element on the whole side, from that node across to its opposite
        // side, moves the mismatch to the opposite side -- one element along
        // the chord -- until it runs out on dS or meets the other end.
        std::unordered_map<MeshEdgeKey, std::vector<int>, MeshEdgeKeyHash> users;
        users.reserve(Q.size() * 4);
        auto addUses = [&](int q) {
            for (int s = 0; s < 4; ++s)
                users[MeshEdgeKey(Q[q][s], Q[q][(s + 1) % 4])].push_back(q);
        };
        auto dropUses = [&](int q) {
            for (int s = 0; s < 4; ++s) {
                std::vector<int> &L = users[MeshEdgeKey(Q[q][s], Q[q][(s + 1) % 4])];
                L.erase(std::remove(L.begin(), L.end(), q), L.end());
            }
        };
        for (int q = 0; q < static_cast<int>(Q.size()); ++q) addUses(q);

        std::unordered_map<MeshEdgeKey, int, MeshEdgeKeyHash> hanging;
        std::deque<MeshEdgeKey> queue;
        int featureSplits = 0;
        if (!closed) {
            for (int i : {0, n - 1}) {
                const MeshEdgeKey key(path[i].v, path[i].exitTo);
                const auto it = users.find(key);
                if (it != users.end() && !it->second.empty()) {
                    hanging[key] = copy[i];
                    queue.push_back(key);
                }
                if (isFeature(path[i].v, path[i].exitTo)) ++featureSplits;
            }
        }

        // A chord crosses each element at most twice, once each way, and there
        // are two ends. More splits than that is a chord running round a closed
        // loop of elements -- one that never reaches dS -- laying a new copy of
        // itself beside the last on every turn.
        int splits = 0;
        const int guardMax = 4 * nQ0 + 64;
        for (int guard = 0; !queue.empty() && ok; ++guard) {
            if (guard > guardMax) {
                refusal = "the chord carried past an end runs round a closed loop of elements "
                          "and never reaches dS";
                ok = false;
                break;
            }
            const MeshEdgeKey key = queue.front();
            queue.pop_front();
            const auto hit = hanging.find(key);
            if (hit == hanging.end()) continue;
            const int h = hit->second;
            hanging.erase(hit);
            const auto uit = users.find(key);
            if (uit == users.end() || uit->second.empty()) continue;
            if (uit->second.size() != 1) {
                refusal = "the chord reached an edge with more than two elements on it";
                ok = false;
                break;
            }
            const int X = uit->second.front();
            const Quad XQ = Q[X];
            int j = -1;
            for (int s = 0; s < 4; ++s)
                if (MeshEdgeKey(XQ[s], XQ[(s + 1) % 4]) == key) { j = s; break; }
            if (j < 0) {
                refusal = "the chord lost the edge it was crossing";
                ok = false;
                break;
            }
            const int a = XQ[j], b = XQ[(j + 1) % 4], c = XQ[(j + 2) % 4], dd = XQ[(j + 3) % 4];

            // The chord crosses the opposite side at the fraction it crossed
            // this one, measured from the corresponding end.
            const Point ab = V[b] - V[a];
            const double L2 = dotP(ab, ab);
            double s = L2 > 0.0 ? dotP(V[h] - V[a], ab) / L2 : 0.5;
            s = std::min(std::max(s, 0.0), 1.0);

            const MeshEdgeKey okey(dd, c);
            int h2 = -1;
            const auto oh = hanging.find(okey);
            if (oh != hanging.end()) {
                // The other side of the opposite edge is already split: the
                // two ends of the chord have met.
                h2 = oh->second;
                hanging.erase(oh);
            } else {
                const bool feature = isFeature(dd, c);
                h2 = static_cast<int>(V.size());
                V.push_back((dd < nV0 && c < nV0) ? onEdge(dd, c, s, feature)
                                                  : V[dd] + (V[c] - V[dd]) * s);
                if (feature) ++featureSplits;
                if (users[okey].size() >= 2) {
                    hanging[okey] = h2;
                    queue.push_back(okey);
                }
            }
            dropUses(X);
            Q[X] = Quad{a, h, h2, dd};
            Q.push_back(Quad{h, b, c, h2});
            M.push_back(M[X]);
            wasFolded.push_back(wasFolded[X]);
            if (X < nQ0) touched[X] = 1;
            addUses(X);
            addUses(static_cast<int>(Q.size()) - 1);
            ++splits;
        }
        if (!ok || !hanging.empty()) {
            if (ok) refusal = "the chord carried past the ends left an edge split on one side only";
            return false;
        }

        // ---- a valid start, or a thinner layer ----------------------------
        int foldedBefore = 0, foldedAfter = 0;
        bool flipped = false, found = false;
        Point where{0.0, 0.0};
        for (int q = 0; q < static_cast<int>(Q.size()); ++q) {
            const bool isNew = q >= nQ0;
            if (!isNew && !touched[q]) continue;
            // A turned-over element is refused whatever it was before: the
            // rebuild would re-wind it (QuadMesh::orientQuads) against its
            // neighbours. A folded corner only counts where nothing was folded.
            const bool flip = signedAreaOf(V, Q[q]) <= 0.0;
            const int folds = wasFolded[q] ? 0 : foldedCorners(V, Q[q]);
            const int before = (isNew || wasFolded[q]) ? 0 : foldedCorners(mesh.vertices, mesh.quads[q]);
            if ((flip || folds > before) && !found) {
                for (int k = 0; k < 4; ++k) where = where + V[Q[q][k]] * 0.25;
                found = true;
            }
            flipped = flipped || flip;
            foldedAfter += folds;
            foldedBefore += before;
        }
        if (flipped || (foldedAfter > foldedBefore && !options.allowFolds)) {
            std::ostringstream oss;
            oss << "every depth tried turned an element over (the first at (" << std::fixed
                << std::setprecision(4) << where[0] << ", " << where[1] << "))";
            refusal = oss.str();
            continue;
        }

        mesh.vertices = std::move(V);
        mesh.quads = std::move(Q);
        mesh.quadMatId = std::move(M);
        mesh.targetJacobian.clear();
        mesh.buildTopology();

        ++report.layers;
        if (closed) ++report.closedLayers;
        report.layerQuads += nE;
        report.splitQuads += splits;
        report.featureSplits += featureSplits;
        return true;
    }
    return false;
}

bool Pillow::run() {
    report = Report();
    report.ran = true;
    report.quadsBefore = static_cast<int>(mesh.quads.size());
    report.verticesBefore = static_cast<int>(mesh.vertices.size());
    report.minScaledJacobianBefore = mesh.quality.minScaledJacobian;
    report.invertedBefore = mesh.quality.invertedQuads;

    std::vector<Defect> defects = findDefects(mesh, options.flatAngle);
    report.defectsBefore = static_cast<int>(defects.size());
    for (const Defect &d : defects) report.flattestBefore = std::max(report.flattestBefore, d.angle);

    // Flattest first, so that the worst corner has its layer even if the cap
    // is reached. A defect no layer could be built for is remembered by its
    // node and element and not tried again; why it was refused is kept, and
    // reported if the corner is still flat at the end -- a later layer may
    // have passed through it.
    std::set<std::pair<int, int>> refused;
    std::map<int, std::string> whyNot;
    while (report.layers < options.maxLayers) {
        defects = findDefects(mesh, options.flatAngle);
        std::stable_sort(defects.begin(), defects.end(),
                         [](const Defect &a, const Defect &b) { return a.angle > b.angle; });
        const Defect *pick = nullptr;
        for (const Defect &d : defects)
            if (!refused.count({d.vertex, d.quad})) { pick = &d; break; }
        if (!pick) break;
        refusal.clear();
        if (!insertLayer(*pick)) {
            refused.insert({pick->vertex, pick->quad});
            whyNot[pick->vertex] = refusal.empty() ? "no reason recorded" : refusal;
        }
    }

    defects = findDefects(mesh, options.flatAngle);
    report.defectsAfter = static_cast<int>(defects.size());
    std::set<int> reported;
    for (const Defect &d : defects) {
        report.flattestAfter = std::max(report.flattestAfter, d.angle);
        const auto it = whyNot.find(d.vertex);
        if (it == whyNot.end() || !reported.insert(d.vertex).second) continue;
        ++report.skipped;
        const Point &p = mesh.vertices[d.vertex];
        std::ostringstream oss;
        oss << "no layer for the " << std::fixed << std::setprecision(0) << d.angle
            << " degree corner at (" << std::setprecision(4) << p[0] << ", " << p[1]
            << "): " << it->second;
        report.messages.push_back(oss.str());
    }
    report.quadsAfter = static_cast<int>(mesh.quads.size());
    report.verticesAfter = static_cast<int>(mesh.vertices.size());
    report.minScaledJacobianAfter = mesh.quality.minScaledJacobian;
    report.invertedAfter = mesh.quality.invertedQuads;

    if (report.layers >= options.maxLayers && report.defectsAfter > 0) {
        std::ostringstream oss;
        oss << "stopped at the cap of " << options.maxLayers << " layers with "
            << report.defectsAfter << " flat corner(s) left";
        report.messages.push_back(oss.str());
    }
    return report.layers > 0;
}

}  // namespace mesh
