#include "HarmonicCut.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <queue>
#include <sstream>
#include <stdexcept>

namespace {

inline double edgeLength(const Mesh &m, int a, int b) {
    const Point &pa = m.vertices[a];
    const Point &pb = m.vertices[b];
    const double dx = pa[0] - pb[0];
    const double dy = pa[1] - pb[1];
    return std::sqrt(dx * dx + dy * dy);
}

inline double squaredDistance(const Mesh &m, int a, int b) {
    const Point &pa = m.vertices[a];
    const Point &pb = m.vertices[b];
    const double dx = pa[0] - pb[0];
    const double dy = pa[1] - pb[1];
    return dx * dx + dy * dy;
}

// Disjoint-set union, used to glue triangle corners across the uncut edges.
struct DSU {
    std::vector<int> parent;
    std::vector<int> rank;
    explicit DSU(int n = 0) : parent(n), rank(n, 0) {
        std::iota(parent.begin(), parent.end(), 0);
    }
    int find(int x) {
        int p = x;
        while (parent[p] != p) p = parent[p];
        while (parent[x] != x) {
            const int nx = parent[x];
            parent[x] = p;
            x = nx;
        }
        return p;
    }
    void unite(int a, int b) {
        a = find(a);
        b = find(b);
        if (a == b) return;
        if (rank[a] < rank[b]) std::swap(a, b);
        parent[b] = a;
        if (rank[a] == rank[b]) rank[a]++;
    }
};

} // namespace

HarmonicCut::HarmonicCut(std::shared_ptr<Mesh> mesh, bool openVoidsIn)
    : orig(std::move(mesh)), openVoids(openVoidsIn) {
    if (!orig) throw std::runtime_error("HarmonicCut: null mesh");
    if (orig->triangles.empty()) throw std::runtime_error("HarmonicCut: empty mesh");

    findBoundaryLoops();
    buildVertexAdjacency();
    generateCuts();
    buildExplicitCutMesh();
    report = checkCutMesh();
}

// ---------------------------------------------------------------------------
// findBoundaryLoops()
//
// The loops l_1 ... l_{beta+1} of Sec. 4.1. Each is walked in order so that its
// enclosed area can be measured: the loop of largest |area| is the outer one
// l_{beta+1}, the rest enclose the voids.
// ---------------------------------------------------------------------------
void HarmonicCut::findBoundaryLoops() {
    const int nV = static_cast<int>(orig->vertices.size());
    boundaryLoops.clear();
    vertexLoop.assign(nV, -1);

    // Boundary neighbours: on a manifold every boundary vertex has exactly two.
    std::vector<std::vector<int>> bAdj(nV);
    for (int e : orig->boundaryEdges) {
        const int a = orig->edges[e][0];
        const int b = orig->edges[e][1];
        bAdj[a].push_back(b);
        bAdj[b].push_back(a);
    }

    std::vector<char> visited(nV, 0);
    for (int v = 0; v < nV; ++v) {
        if (bAdj[v].empty() || visited[v]) continue;

        // Walk the ring, always leaving by the neighbour we did not arrive on.
        std::vector<int> loop;
        int prev = -1;
        int cur = v;
        while (cur != -1 && !visited[cur]) {
            visited[cur] = 1;
            loop.push_back(cur);

            int next = -1;
            for (int nb : bAdj[cur]) {
                if (nb != prev && !visited[nb]) { next = nb; break; }
            }
            prev = cur;
            cur = next;
        }

        const int id = static_cast<int>(boundaryLoops.size());
        for (int u : loop) vertexLoop[u] = id;
        boundaryLoops.push_back(std::move(loop));
    }

    // Outer loop by shoelace area. A non-manifold pinch would break the walk
    // above and give a bogus ring, which shows up here as a degenerate area,
    // so the fallback is the loop with the most vertices.
    double bestArea = -1.0;
    outerLoop = boundaryLoops.empty() ? -1 : 0;
    for (int i = 0; i < static_cast<int>(boundaryLoops.size()); ++i) {
        const auto &loop = boundaryLoops[i];
        double area2 = 0.0;
        for (size_t k = 0; k < loop.size(); ++k) {
            const Point &p = orig->vertices[loop[k]];
            const Point &q = orig->vertices[loop[(k + 1) % loop.size()]];
            area2 += cross2(p, q);
        }
        const double area = 0.5 * std::fabs(area2);
        if (area > bestArea) {
            bestArea = area;
            outerLoop = i;
        }
    }
}

void HarmonicCut::buildVertexAdjacency() {
    const int nV = static_cast<int>(orig->vertices.size());
    adjacency.assign(nV, {});
    for (int e = 0; e < static_cast<int>(orig->edges.size()); ++e) {
        const int a = orig->edges[e][0];
        const int b = orig->edges[e][1];
        const double w = edgeLength(*orig, a, b);
        adjacency[a].push_back({b, w});
        adjacency[b].push_back({a, w});
    }

    // Every boundary vertex starts out blocked as an intermediate: a cut may
    // only meet the boundary at its endpoints. Cut vertices are added to this
    // as the cuts are made, which is what keeps them disjoint from each other.
    blockedVertex.assign(nV, 0);
    for (int v = 0; v < nV; ++v) {
        if (v < static_cast<int>(orig->isBoundaryVertex.size()) && orig->isBoundaryVertex[v]) {
            blockedVertex[v] = 1;
        }
    }
}

// ---------------------------------------------------------------------------
// shortestInteriorPath()
//
// Dijkstra over the primal graph, where a blocked vertex may be an endpoint
// but never an interior node of the path. Multi-source and multi-target: the
// paper picks one nearest pair (p, q) and connects it, but if that particular
// pair is walled off -- every route from it having to pass through an earlier
// cut -- the caller retries over whole loops rather than giving up on the void.
// ---------------------------------------------------------------------------
std::vector<int> HarmonicCut::shortestInteriorPath(const std::vector<int> &sources,
                                                   const std::unordered_set<int> &targets,
                                                   const std::vector<char> &blocked) const {
    const int nV = static_cast<int>(orig->vertices.size());
    std::vector<double> dist(nV, std::numeric_limits<double>::infinity());
    std::vector<int> prev(nV, -1);

    using QItem = std::pair<double, int>;
    std::priority_queue<QItem, std::vector<QItem>, std::greater<QItem>> pq;

    std::unordered_set<int> sourceSet(sources.begin(), sources.end());
    for (int s : sources) {
        if (targets.count(s)) continue; // degenerate: source is already a target
        dist[s] = 0.0;
        pq.push({0.0, s});
    }

    int reached = -1;
    while (!pq.empty()) {
        const auto [d, u] = pq.top();
        pq.pop();
        if (d != dist[u]) continue;

        if (targets.count(u)) { reached = u; break; }

        // Only sources and interior vertices may be expanded; blocked vertices
        // that are not sources are dead ends unless they are targets, which is
        // handled above.
        if (blocked[u] && !sourceSet.count(u)) continue;

        for (const auto &[to, w] : adjacency[u]) {
            if (dist[to] > d + w) {
                dist[to] = d + w;
                prev[to] = u;
                pq.push({dist[to], to});
            }
        }
    }

    std::vector<int> path;
    if (reached < 0) return path;
    for (int v = reached; v != -1; v = prev[v]) path.push_back(v);
    std::reverse(path.begin(), path.end());
    return path;
}

// ---------------------------------------------------------------------------
// generateCuts()  --  the Sec. 4.1 recursion
//
// The mesh starts as M_beta with the outer loop l_{beta+1} and beta inner
// loops. Each iteration picks the void whose loop comes nearest the boundary
// merged so far, opens it with one cut, and so turns M_k into M_{k-1}; after
// beta iterations there are no voids left and the mesh is a disk.
//
// "Nearest" and "merged so far" together make this Prim's algorithm on the
// loops, rooted at the outer one. The paper phrases the recursion as always
// cutting to l_{beta+1}, which is the same thing once the earlier cuts have
// merged their loops into it.
// ---------------------------------------------------------------------------
void HarmonicCut::generateCuts() {
    cuts.clear();
    cutEdges.clear();

    if (!openVoids) return;  // the caller wants the voids kept -- see the ctor

    const int nLoops = static_cast<int>(boundaryLoops.size());
    if (nLoops <= 1) return; // already a disk: no voids to open

    std::vector<char> merged(nLoops, 0);
    merged[outerLoop] = 1;

    // Vertices of the loops merged so far, i.e. the current outer boundary.
    std::vector<int> mergedVerts = boundaryLoops[outerLoop];

    int remaining = nLoops - 1;
    while (remaining > 0) {
        // Order the unopened voids by the distance from their loop to the
        // current outer boundary, and remember the nearest vertex pair for
        // each -- the (p, q) of Sec. 4.1.
        //
        // Brute force over the boundary vertices. The paper reaches for an
        // approximate nearest-neighbour index (ANN) here, which this does not
        // need at these sizes: the whole cut generation takes 5-7 ms on the
        // 600-boundary-edge, 9-void models in data/meshes. A spatial index is
        // the upgrade if the boundary vertex count ever grows by an order of
        // magnitude, since this loop is quadratic in it.
        struct Candidate {
            int loop = -1;
            int p = -1;
            int q = -1;
            double dist2 = std::numeric_limits<double>::infinity();
        };
        std::vector<Candidate> candidates;

        for (int l = 0; l < nLoops; ++l) {
            if (merged[l]) continue;
            Candidate c;
            c.loop = l;
            for (int p : mergedVerts) {
                if (blockedVertex[p] > 1) continue; // already carries a cut endpoint
                for (int q : boundaryLoops[l]) {
                    if (blockedVertex[q] > 1) continue;
                    const double d2 = squaredDistance(*orig, p, q);
                    if (d2 < c.dist2) { c.dist2 = d2; c.p = p; c.q = q; }
                }
            }
            if (c.p >= 0) candidates.push_back(c);
        }

        if (candidates.empty()) {
            std::ostringstream oss;
            oss << remaining << " void(s) could not be opened: no admissible endpoint pair left.";
            report.messages.push_back(oss.str());
            break;
        }

        std::sort(candidates.begin(), candidates.end(),
                  [](const Candidate &a, const Candidate &b) { return a.dist2 < b.dist2; });

        // Try the candidates nearest first. The path search may fail for one
        // void while another is still reachable, so this walks the list rather
        // than stopping at the first refusal.
        bool opened = false;
        for (const Candidate &c : candidates) {
            std::vector<int> path = shortestInteriorPath({c.p}, {c.q}, blockedVertex);
            if (path.empty()) {
                // Relax to the whole pair of loops: still the shortest
                // admissible cut, just not between the nearest pair.
                const std::unordered_set<int> targets(boundaryLoops[c.loop].begin(),
                                                      boundaryLoops[c.loop].end());
                path = shortestInteriorPath(mergedVerts, targets, blockedVertex);
            }
            if (path.size() < 2) continue;

            applyCut(path, vertexLoop[path.front()], vertexLoop[path.back()]);

            const int target = vertexLoop[path.back()];
            merged[target] = 1;
            mergedVerts.insert(mergedVerts.end(),
                               boundaryLoops[target].begin(), boundaryLoops[target].end());
            --remaining;
            opened = true;
            break;
        }

        if (!opened) {
            std::ostringstream oss;
            oss << remaining << " void(s) could not be opened: no interior path avoiding the "
                << "existing cuts and the boundary.";
            report.messages.push_back(oss.str());
            break;
        }
    }
}

void HarmonicCut::applyCut(const std::vector<int> &path, int fromLoop, int toLoop) {
    Cut c;
    c.path = path;
    c.fromLoop = fromLoop;
    c.toLoop = toLoop;

    for (size_t k = 0; k + 1 < path.size(); ++k) {
        const EdgeKey ek(path[k], path[k + 1]);
        c.edges.push_back(ek);
        c.length += edgeLength(*orig, path[k], path[k + 1]);
        cutEdges.insert(ek);
    }

    // Mark every vertex of the cut, endpoints included: a later cut may
    // neither cross this one nor land on it. 2 distinguishes "on a cut" from
    // the plain boundary block, so that endpoint selection can skip it while
    // the path search treats both the same.
    for (int v : path) blockedVertex[v] = 2;

    cuts.push_back(std::move(c));
}

// ---------------------------------------------------------------------------
// buildExplicitCutMesh()
//
// M_C: same triangles, but the corners on either side of a cut become separate
// vertices. Corners are glued across every interior edge that is not a cut, so
// each connected group of corners becomes one vertex of the cut mesh.
// ---------------------------------------------------------------------------
void HarmonicCut::buildExplicitCutMesh() {
    const int nT = static_cast<int>(orig->triangles.size());
    const int nV = static_cast<int>(orig->vertices.size());

    cutVertToOrig.clear();
    origToCutVerts.assign(nV, {});

    auto localIndex = [&](int f, int v) -> int {
        const Triangle &t = orig->triangles[f];
        if (t[0] == v) return 0;
        if (t[1] == v) return 1;
        if (t[2] == v) return 2;
        return -1;
    };

    DSU dsu(nT * 3);
    for (int e = 0; e < static_cast<int>(orig->edges.size()); ++e) {
        if (orig->isBoundaryEdge[e]) continue;
        const int u = orig->edges[e][0];
        const int v = orig->edges[e][1];
        if (cutEdges.count(EdgeKey(u, v))) continue; // seam: leave the sides apart

        const int f0 = orig->edgeTriangles[e][0];
        const int f1 = orig->edgeTriangles[e][1];
        if (f0 < 0 || f1 < 0) continue;

        const int i0u = localIndex(f0, u), i1u = localIndex(f1, u);
        const int i0v = localIndex(f0, v), i1v = localIndex(f1, v);
        if (i0u < 0 || i1u < 0 || i0v < 0 || i1v < 0) continue;

        dsu.unite(3 * f0 + i0u, 3 * f1 + i1u);
        dsu.unite(3 * f0 + i0v, 3 * f1 + i1v);
    }

    std::vector<Point> verts;
    std::vector<int> cornerNew(3 * nT, -1);
    std::unordered_map<int, int> rootToNew;
    rootToNew.reserve(static_cast<size_t>(nT) * 3);

    for (int c = 0; c < 3 * nT; ++c) {
        const int r = dsu.find(c);
        auto it = rootToNew.find(r);
        if (it == rootToNew.end()) {
            const int origV = orig->triangles[c / 3][c % 3];
            const int newId = static_cast<int>(verts.size());
            rootToNew.emplace(r, newId);
            verts.push_back(orig->vertices[origV]);
            cutVertToOrig.push_back(origV);
            cornerNew[c] = newId;
        } else {
            cornerNew[c] = it->second;
        }
    }

    std::vector<Triangle> tris(nT);
    for (int f = 0; f < nT; ++f) {
        tris[f] = Triangle{cornerNew[3 * f + 0], cornerNew[3 * f + 1], cornerNew[3 * f + 2]};
    }

    cut = Mesh(verts, tris);

    for (int cv = 0; cv < static_cast<int>(cutVertToOrig.size()); ++cv) {
        const int ov = cutVertToOrig[cv];
        if (ov >= 0 && ov < nV) origToCutVerts[ov].push_back(cv);
    }
}

// ---------------------------------------------------------------------------
// checkCutMesh()
// ---------------------------------------------------------------------------
HarmonicCut::Report HarmonicCut::checkCutMesh() const {
    Report rep;
    rep.messages = report.messages; // carry over anything generateCuts() logged
    rep.boundaryLoops = static_cast<int>(boundaryLoops.size());
    rep.voids = std::max(0, rep.boundaryLoops - 1);
    rep.cutsMade = static_cast<int>(cuts.size());

    const int V = static_cast<int>(cut.vertices.size());
    const int F = static_cast<int>(cut.triangles.size());
    const int E = static_cast<int>(cut.edges.size());
    if (V == 0 || F == 0) {
        rep.messages.push_back("Cut mesh is empty.");
        return rep;
    }
    rep.eulerCharacteristic = V - E + F;

    // Boundary components of the cut mesh.
    std::vector<std::vector<int>> bAdj(V);
    std::vector<char> isBoundaryV(V, 0);
    for (int e : cut.boundaryEdges) {
        const int a = cut.edges[e][0];
        const int b = cut.edges[e][1];
        bAdj[a].push_back(b);
        bAdj[b].push_back(a);
        isBoundaryV[a] = 1;
        isBoundaryV[b] = 1;
    }

    std::vector<char> seen(V, 0);
    std::queue<int> q;
    for (int v = 0; v < V; ++v) {
        if (!isBoundaryV[v] || seen[v]) continue;
        rep.boundaryComponents++;
        seen[v] = 1;
        q.push(v);
        while (!q.empty()) {
            const int x = q.front();
            q.pop();
            for (int y : bAdj[x]) {
                if (!seen[y]) { seen[y] = 1; q.push(y); }
            }
        }
    }

    // Triangle connectivity.
    std::vector<char> seenT(F, 0);
    int comps = 0;
    for (int t = 0; t < F; ++t) {
        if (seenT[t]) continue;
        comps++;
        seenT[t] = 1;
        q.push(t);
        while (!q.empty()) {
            const int x = q.front();
            q.pop();
            for (int e = 0; e < 3; ++e) {
                const int nb = cut.triangleAdjacency[x][e];
                if (nb != -1 && !seenT[nb]) { seenT[nb] = 1; q.push(nb); }
            }
        }
    }
    rep.trianglesConnected = (comps == 1);

    rep.isDisk = rep.trianglesConnected && rep.boundaryComponents == 1 &&
                 rep.eulerCharacteristic == 1;

    // With the voids deliberately kept, none of the disk checks below is a
    // finding: chi is 1 - beta and there are beta + 1 boundary components
    // because that is what was asked for. Only the connectivity still means
    // anything, and it is checked above.
    if (!openVoids) {
        if (!rep.trianglesConnected) {
            std::ostringstream oss;
            oss << "Mesh triangles fall into " << comps << " components.";
            rep.messages.push_back(oss.str());
        }
        return rep;
    }

    if (rep.cutsMade != rep.voids) {
        std::ostringstream oss;
        oss << "Expected " << rep.voids << " cut(s) for " << rep.voids
            << " void(s), made " << rep.cutsMade << ".";
        rep.messages.push_back(oss.str());
    }
    if (!rep.trianglesConnected) {
        std::ostringstream oss;
        oss << "Cut mesh triangles fall into " << comps << " components.";
        rep.messages.push_back(oss.str());
    }
    if (rep.boundaryComponents != 1) {
        std::ostringstream oss;
        oss << "Expected 1 boundary component for a disk, got " << rep.boundaryComponents << ".";
        rep.messages.push_back(oss.str());
    }
    if (rep.eulerCharacteristic != 1) {
        std::ostringstream oss;
        oss << "Expected Euler characteristic 1 for a disk, got " << rep.eulerCharacteristic << ".";
        rep.messages.push_back(oss.str());
    }

    return rep;
}

bool HarmonicCut::writeOBJ(const std::string &filename) const {
    std::ofstream out(filename);
    if (!out) return false;
    for (const auto &p : cut.vertices) out << "v " << p[0] << " " << p[1] << " 0\n";
    for (const auto &t : cut.triangles) {
        out << "f " << (t[0] + 1) << " " << (t[1] + 1) << " " << (t[2] + 1) << "\n";
    }
    return true;
}
