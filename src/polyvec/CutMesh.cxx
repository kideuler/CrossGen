#include "CutMesh.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <queue>
#include <sstream>

#include "PolyVectors.hxx" // PolyField

namespace {

static inline double edgeLength(const Mesh &m, int a, int b) {
    const Point &pa = m.vertices[a];
    const Point &pb = m.vertices[b];
    const double dx = pa[0] - pb[0];
    const double dy = pa[1] - pb[1];
    return std::sqrt(dx * dx + dy * dy);
}

// Disjoint-set union for corners
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
            int nx = parent[x];
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

// Fill Mesh::triangleAdjacency, boundaryTriangles, cornerTriangles, vertexTriangles from vertices+triangles.
static void computeTopology(Mesh &m) {
    const int nT = static_cast<int>(m.triangles.size());
    m.triangleAdjacency.assign(nT, std::array<int, 3>{-1, -1, -1});
    m.boundaryTriangles.clear();
    m.cornerTriangles.clear();

    struct EdgeInfo { int tri; int edgeId; };
    std::unordered_map<CutMesh::EdgeKey, EdgeInfo, CutMesh::EdgeKeyHash> edgeMap;
    edgeMap.reserve(static_cast<size_t>(nT) * 3);

    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = m.triangles[t];
        for (int e = 0; e < 3; ++e) {
            const int a = tri[e];
            const int b = tri[(e + 1) % 3];
            const CutMesh::EdgeKey key(a, b);
            auto it = edgeMap.find(key);
            if (it == edgeMap.end()) {
                edgeMap.emplace(key, EdgeInfo{t, e});
            } else {
                // We assume manifold edges (max 2 incident triangles). If we see more,
                // just keep the first pair.
                const int ot = it->second.tri;
                const int oe = it->second.edgeId;
                if (m.triangleAdjacency[t][e] == -1 && m.triangleAdjacency[ot][oe] == -1) {
                    m.triangleAdjacency[t][e] = ot;
                    m.triangleAdjacency[ot][oe] = t;
                }
            }
        }
    }

    // Boundary/corner triangle bookkeeping
    for (int t = 0; t < nT; ++t) {
        int bEdges[3];
        int nb = 0;
        for (int e = 0; e < 3; ++e) {
            if (m.triangleAdjacency[t][e] == -1) bEdges[nb++] = e;
        }
        if (nb == 1) {
            m.boundaryTriangles.push_back(std::array<int, 2>{t, bEdges[0]});
        } else if (nb >= 2) {
            m.cornerTriangles.push_back(std::array<int, 3>{t, bEdges[0], bEdges[1]});
        }
    }

    m.vertexTriangles = VertexTriangleCSR::buildFromMesh(m);
}

} // namespace

CutMesh::CutMesh(const PolyField &field) {
    orig = field.getMesh();
    singularities = field.uSingularities;

    buildDualSpanningTreeCuts();
    connectSingularitiesWithShortestPaths();
    buildExplicitCutMesh();
}

void CutMesh::buildDualSpanningTreeCuts() {
    cutEdges.clear();
    singularityPathCutEdges.clear();
    const int nT = static_cast<int>(orig.triangles.size());
    if (nT == 0) return;

    std::unordered_set<EdgeKey, EdgeKeyHash> interiorEdges;
    interiorEdges.reserve(static_cast<size_t>(nT) * 3);
    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = orig.triangles[t];
        for (int e = 0; e < 3; ++e) {
            if (orig.triangleAdjacency[t][e] == -1) continue;
            interiorEdges.insert(EdgeKey(tri[e], tri[(e + 1) % 3]));
        }
    }

    // Dual spanning forest: triangles are nodes, adjacency across interior edges are dual edges.
    std::unordered_set<EdgeKey, EdgeKeyHash> treeEdges;
    treeEdges.reserve(static_cast<size_t>(nT) * 2);
    std::vector<char> visited(nT, false);
    std::queue<int> q;
    for (int start = 0; start < nT; ++start) {
        if (visited[start]) continue;
        visited[start] = true;
        q.push(start);
        while (!q.empty()) {
            const int t = q.front();
            q.pop();
            const Triangle &tri = orig.triangles[t];
            for (int e = 0; e < 3; ++e) {
                const int nb = orig.triangleAdjacency[t][e];
                if (nb == -1) continue;
                if (!visited[nb]) {
                    visited[nb] = true;
                    q.push(nb);
                    treeEdges.insert(EdgeKey(tri[e], tri[(e + 1) % 3]));
                }
            }
        }
    }

    // Cut all interior edges not in the dual spanning forest.
    cutEdges = std::move(interiorEdges);
    for (const auto &e : treeEdges) {
        cutEdges.erase(e);
    }
}

void CutMesh::connectSingularitiesWithShortestPaths() {
    const int nV = static_cast<int>(orig.vertices.size());
    const int nT = static_cast<int>(orig.triangles.size());
    if (nV == 0 || nT == 0) return;

    // Record only the additional cuts introduced by the singularity-connection paths.
    singularityPathCutEdges.clear();

    // Build unique primal edges with weights.
    std::unordered_map<EdgeKey, double, EdgeKeyHash> wmap;
    wmap.reserve(static_cast<size_t>(nT) * 3);
    auto addEdge = [&](int a, int b) {
        const EdgeKey key(a, b);
        if (wmap.find(key) != wmap.end()) return;
        wmap.emplace(key, edgeLength(orig, a, b));
    };
    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = orig.triangles[t];
        addEdge(tri[0], tri[1]);
        addEdge(tri[1], tri[2]);
        addEdge(tri[2], tri[0]);
    }

    std::vector<std::vector<std::pair<int, double>>> adj(nV);
    adj.reserve(nV);
    for (const auto &kv : wmap) {
        const int a = kv.first.a;
        const int b = kv.first.b;
        const double w = kv.second;
        if (a >= 0 && a < nV && b >= 0 && b < nV) {
            adj[a].push_back({b, w});
            adj[b].push_back({a, w});
        }
    }

    // Boundary markers: ONLY the original mesh boundary.
    //
    // Important: We want paths that truly connect each singularity to the external
    // boundary of the input mesh. If we also treat the current cut graph as
    // boundary, Dijkstra may terminate immediately (or very early), producing no
    // singularity-to-boundary path cuts.
    std::vector<char> isBoundary(nV, false);
    auto markBoundaryEdge = [&](int a, int b) {
        if (a >= 0 && a < nV) isBoundary[a] = true;
        if (b >= 0 && b < nV) isBoundary[b] = true;
    };
    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = orig.triangles[t];
        for (int e = 0; e < 3; ++e) {
            if (orig.triangleAdjacency[t][e] != -1) continue;
            markBoundaryEdge(tri[e], tri[(e + 1) % 3]);
        }
    }

    // Dijkstra from each singularity to the original boundary.
    using QItem = std::pair<double, int>;
    auto runDijkstra = [&](int source, const std::vector<char> &isB, std::vector<int> &prev) -> int {
        prev.assign(nV, -1);
        std::vector<double> dist(nV, std::numeric_limits<double>::infinity());
        std::priority_queue<QItem, std::vector<QItem>, std::greater<QItem>> pq;
        dist[source] = 0.0;
        pq.push({0.0, source});
        while (!pq.empty()) {
            const auto [d, v] = pq.top();
            pq.pop();
            if (d != dist[v]) continue;
            if (isB[v]) return v;
            for (const auto &nw : adj[v]) {
                const int to = nw.first;
                const double w = nw.second;
                if (dist[to] > d + w) {
                    dist[to] = d + w;
                    prev[to] = v;
                    pq.push({dist[to], to});
                }
            }
        }
        return -1;
    };

    std::vector<int> prev;
    for (const auto &s : singularities) {
        const int vSing = s.first;
        if (vSing < 0 || vSing >= nV) continue;
        if (isBoundary[vSing]) continue;

        const int target = runDijkstra(vSing, isBoundary, prev);
        if (target < 0) {
            // Disconnected graph; can't connect this singularity.
            continue;
        }

        // Backtrack target -> source and mark all edges on the path as cuts.
        int cur = target;
        while (cur != vSing) {
            const int p = prev[cur];
            if (p < 0) break;
            const EdgeKey ek(cur, p);
            cutEdges.insert(ek);
            singularityPathCutEdges.insert(ek);
            cur = p;
        }
        // Note: we intentionally don't augment the boundary target set with the
        // newly cut path. Each singularity should connect to the true boundary.
    }
}

void CutMesh::buildExplicitCutMesh() {
    const int nT = static_cast<int>(orig.triangles.size());
    const int nV = static_cast<int>(orig.vertices.size());
    cut = Mesh();
    cut.vertices.clear();
    cut.triangles.clear();
    cutVertToOrig.clear();
    origToCutVerts.clear();

    if (nT == 0 || nV == 0) return;

    DSU dsu(nT * 3);

    // Glue triangles across non-cut interior edges by unifying the corresponding corners.
    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = orig.triangles[t];
        for (int e = 0; e < 3; ++e) {
            const int nb = orig.triangleAdjacency[t][e];
            if (nb == -1) continue; // original boundary

            const int a = tri[e];
            const int b = tri[(e + 1) % 3];
            if (singularityPathCutEdges.find(EdgeKey(a, b)) != singularityPathCutEdges.end()) continue;

            // Find the local indices of a and b in the neighbor triangle.
            const Triangle &triN = orig.triangles[nb];
            int ia = -1, ib = -1;
            for (int lv = 0; lv < 3; ++lv) {
                if (triN[lv] == a) ia = lv;
                if (triN[lv] == b) ib = lv;
            }
            if (ia == -1 || ib == -1) continue;

            dsu.unite(3 * t + e, 3 * nb + ia);
            dsu.unite(3 * t + ((e + 1) % 3), 3 * nb + ib);
        }
    }

    // Assign new vertex ids to DSU roots.
    std::unordered_map<int, int> rootToNew;
    rootToNew.reserve(static_cast<size_t>(nT) * 3);

    cut.triangles.resize(nT);
    for (int t = 0; t < nT; ++t) {
        for (int lv = 0; lv < 3; ++lv) {
            const int corner = 3 * t + lv;
            const int root = dsu.find(corner);
            auto it = rootToNew.find(root);
            int newV = -1;
            if (it == rootToNew.end()) {
                newV = static_cast<int>(cut.vertices.size());
                rootToNew.emplace(root, newV);
                const int origV = orig.triangles[t][lv];
                cut.vertices.push_back(orig.vertices[origV]);
                cutVertToOrig.push_back(origV);
            } else {
                newV = it->second;
            }
            cut.triangles[t][lv] = newV;
        }
    }

    // Build inverse mapping.
    origToCutVerts.assign(static_cast<size_t>(nV), {});
    for (int cv = 0; cv < static_cast<int>(cutVertToOrig.size()); ++cv) {
        const int ov = cutVertToOrig[cv];
        if (ov >= 0 && ov < nV) {
            origToCutVerts[ov].push_back(cv);
        }
    }

    // Recompute topology on the explicit cut mesh.
    computeTopology(cut);
}

bool CutMesh::writeOBJ(const std::string &filename) const {
    std::ofstream out(filename);
    if (!out) return false;

    for (const auto &p : cut.vertices) {
        out << "v " << p[0] << " " << p[1] << " 0\n";
    }
    for (const auto &tri : cut.triangles) {
        out << "f " << (tri[0] + 1) << " " << (tri[1] + 1) << " " << (tri[2] + 1) << "\n";
    }
    return true;
}

CutMesh::SanityReport CutMesh::sanityCheck() const {
    SanityReport rep;
    const int V = static_cast<int>(cut.vertices.size());
    const int F = static_cast<int>(cut.triangles.size());
    if (V == 0 || F == 0) {
        rep.messages.push_back("Cut mesh is empty.");
        return rep;
    }

    // Count unique edges and boundary edges (edges appearing only once).
    std::unordered_map<EdgeKey, int, EdgeKeyHash> edgeCount;
    edgeCount.reserve(static_cast<size_t>(F) * 3);
    auto addEdge = [&](int a, int b) {
        const EdgeKey k(a, b);
        edgeCount[k] += 1;
    };
    for (const auto &tri : cut.triangles) {
        addEdge(tri[0], tri[1]);
        addEdge(tri[1], tri[2]);
        addEdge(tri[2], tri[0]);
    }
    const int E = static_cast<int>(edgeCount.size());
    rep.eulerCharacteristic = V - E + F;

    // Boundary graph.
    std::vector<std::vector<int>> bAdj(V);
    std::vector<char> isBoundaryV(V, false);
    for (const auto &kv : edgeCount) {
        if (kv.second == 1) {
            const int a = kv.first.a;
            const int b = kv.first.b;
            if (a >= 0 && a < V && b >= 0 && b < V) {
                bAdj[a].push_back(b);
                bAdj[b].push_back(a);
                isBoundaryV[a] = true;
                isBoundaryV[b] = true;
            }
        }
    }

    // Count boundary components.
    std::vector<char> visB(V, false);
    int bComp = 0;
    std::queue<int> q;
    for (int v = 0; v < V; ++v) {
        if (!isBoundaryV[v] || visB[v]) continue;
        bComp++;
        visB[v] = true;
        q.push(v);
        while (!q.empty()) {
            int x = q.front();
            q.pop();
            for (int y : bAdj[x]) {
                if (!visB[y]) {
                    visB[y] = true;
                    q.push(y);
                }
            }
        }
    }
    rep.boundaryComponents = bComp;

    // Triangle connectivity components.
    const int nT = static_cast<int>(cut.triangles.size());
    std::vector<char> visT(nT, false);
    int tComp = 0;
    for (int t = 0; t < nT; ++t) {
        if (visT[t]) continue;
        tComp++;
        visT[t] = true;
        q.push(t);
        while (!q.empty()) {
            int x = q.front();
            q.pop();
            const auto &adj = cut.triangleAdjacency[x];
            for (int e = 0; e < 3; ++e) {
                int nb = adj[e];
                if (nb != -1 && !visT[nb]) {
                    visT[nb] = true;
                    q.push(nb);
                }
            }
        }
    }
    rep.triangleComponents = tComp;
    rep.trianglesConnected = (tComp == 1);

    // Singularities on boundary.
    bool singOk = true;
    for (const auto &s : singularities) {
        const int ov = s.first;
        if (ov < 0 || ov >= static_cast<int>(origToCutVerts.size())) continue;
        bool any = false;
        for (int cv : origToCutVerts[ov]) {
            if (cv >= 0 && cv < V && isBoundaryV[cv]) {
                any = true;
                break;
            }
        }
        if (!any) {
            singOk = false;
            std::ostringstream oss;
            oss << "Singularity (orig vertex " << ov << ") is not on the cut mesh boundary.";
            rep.messages.push_back(oss.str());
        }
    }
    rep.allSingularitiesOnBoundary = singOk;

    // Disk heuristic.
    rep.looksLikeDisk = rep.trianglesConnected && rep.boundaryComponents == 1 && rep.eulerCharacteristic == 1;
    if (!rep.trianglesConnected) rep.messages.push_back("Cut mesh triangles are not in a single connected component.");
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
