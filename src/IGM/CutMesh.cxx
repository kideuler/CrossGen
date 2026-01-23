#include "CutMesh.hxx"

#include <algorithm>
#include <cmath>
#include <deque>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>

#include "polyvector/PolyVectors.hxx" // PolyField

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

    // Copy the per-triangle u field directions
    const auto& fieldVecs = field.field;
    uField.resize(fieldVecs.size());
    for (size_t i = 0; i < fieldVecs.size(); ++i) {
        uField[i] = fieldVecs[i].u;
    }

    buildEdgeCuts();
    connectSingularitiesWithShortestPaths();
    buildExplicitCutMesh();
    sanityInfo = sanityCheck();

    if (sanityInfo.looksLikeDisk) {
        combFieldDirections();
    }
    makeVFieldFromUField();
}

void CutMesh::buildEdgeCuts() {
    cutEdges.clear();
    singularityPathCutEdges.clear();

    const int nV = static_cast<int>(orig.vertices.size());
    const int nT = static_cast<int>(orig.triangles.size());
    if (nV == 0 || nT == 0) return;

    // Internal edge structure
    struct InternalEdge {
        int u = -1, v = -1;   // endpoints (u < v)
        int f0 = -1, f1 = -1; // incident faces (f1 == -1 => boundary)
        int bc = -1;          // boundary component id for boundary edges
    };
    struct DualArc { int to; int eid; };

    // -----------------------
    // 1) Build unique undirected edge table (+ incident faces)
    // -----------------------
    std::unordered_map<EdgeKey, int, EdgeKeyHash> edgeMap;
    edgeMap.reserve(static_cast<size_t>(nT) * 3);

    std::vector<InternalEdge> edges;
    edges.reserve(static_cast<size_t>(nT * 3 / 2));

    auto addEdge = [&](int a, int b, int fidx) {
        EdgeKey k(a, b);
        auto it = edgeMap.find(k);
        if (it == edgeMap.end()) {
            int eid = static_cast<int>(edges.size());
            edgeMap.emplace(k, eid);
            InternalEdge e;
            e.u = k.a; e.v = k.b;
            e.f0 = fidx; e.f1 = -1;
            edges.push_back(e);
        } else {
            InternalEdge& e = edges[it->second];
            if (e.f1 == -1) e.f1 = fidx;
        }
    };

    for (int f = 0; f < nT; ++f) {
        const Triangle& t = orig.triangles[f];
        addEdge(t[0], t[1], f);
        addEdge(t[1], t[2], f);
        addEdge(t[2], t[0], f);
    }

    const int nE = static_cast<int>(edges.size());

    // Vertex -> incident edges
    std::vector<std::vector<int>> v2e(nV);
    for (int eid = 0; eid < nE; ++eid) {
        v2e[edges[eid].u].push_back(eid);
        v2e[edges[eid].v].push_back(eid);
    }

    // Track used vertices
    std::vector<uint8_t> usedV(nV, 0);
    for (const auto& t : orig.triangles) {
        usedV[t[0]] = 1; usedV[t[1]] = 1; usedV[t[2]] = 1;
    }
    int Vused = 0;
    for (auto b : usedV) if (b) Vused++;

    // -----------------------
    // 2) Identify boundary components
    // -----------------------
    DSU dsuB(nV);
    for (int eid = 0; eid < nE; ++eid) {
        if (edges[eid].f1 == -1) dsuB.unite(edges[eid].u, edges[eid].v);
    }

    std::vector<int> rootToBC(nV, -1);
    int bCount = 0;
    for (int eid = 0; eid < nE; ++eid) {
        if (edges[eid].f1 != -1) continue;
        int r = dsuB.find(edges[eid].u);
        if (rootToBC[r] == -1) rootToBC[r] = bCount++;
        edges[eid].bc = rootToBC[r];
    }

    std::vector<int> bcRep(bCount, -1);
    for (int eid = 0; eid < nE; ++eid) {
        if (edges[eid].f1 != -1) continue;
        int bc = edges[eid].bc;
        if (bc >= 0 && bc < bCount && bcRep[bc] == -1) bcRep[bc] = edges[eid].u;
    }
    
    // Terminal vertices: one representative per boundary component (for pruning protection)
    std::vector<uint8_t> terminal(nV, 0);
    for (int bc = 0; bc < bCount; ++bc) {
        if (bcRep[bc] >= 0 && bcRep[bc] < nV) terminal[bcRep[bc]] = 1;
    }

    // Compute Euler characteristic
    const int chi = Vused - nE + nT;
    const double g_est = (2.0 - static_cast<double>(bCount) - static_cast<double>(chi)) / 2.0;

    // If already a disk, no cuts needed
    if (bCount == 1 && std::abs(g_est) < 1e-6) {
        return;
    }

    // -----------------------
    // 3) Build augmented dual graph (faces + boundary component nodes)
    // -----------------------
    const int nD = nT + bCount;
    std::vector<std::vector<DualArc>> dualAdj(nD);

    for (int eid = 0; eid < nE; ++eid) {
        const InternalEdge& e = edges[eid];
        if (e.f1 != -1) {
            // Interior edge: connects two faces
            dualAdj[e.f0].push_back({e.f1, eid});
            dualAdj[e.f1].push_back({e.f0, eid});
        } else {
            // Boundary edge: connects face to its boundary component node
            int bd = nT + e.bc;
            dualAdj[e.f0].push_back({bd, eid});
            dualAdj[bd].push_back({e.f0, eid});
        }
    }

    // -----------------------
    // 4) Primal spanning tree T 
    //    Strategy: MAXIMIZE boundary edge usage to leave interior edges for C*
    //    Phase 1: Use boundary edges to create spanning forest on boundary vertices
    //    Phase 2: Connect interior vertices with interior edges
    //    Phase 3: Connect boundary forest to interior tree
    // -----------------------
    std::vector<uint8_t> inT(nE, 0);
    std::vector<uint8_t> visV(nV, 0);
    
    // Track which boundary component each vertex belongs to (-1 for interior)
    std::vector<int> vBC(nV, -1);
    for (int eid = 0; eid < nE; ++eid) {
        if (edges[eid].f1 != -1) continue;
        int bc = edges[eid].bc;
        vBC[edges[eid].u] = bc;
        vBC[edges[eid].v] = bc;
    }
    
    // Phase 1: Span boundary vertices using boundary edges (one tree per boundary component)
    // For each boundary component, do BFS using only boundary edges
    for (int bc = 0; bc < bCount; ++bc) {
        // Find a starting vertex on this boundary
        int startV = -1;
        for (int v = 0; v < nV; ++v) {
            if (vBC[v] == bc && !visV[v]) { startV = v; break; }
        }
        if (startV == -1) continue;
        
        visV[startV] = 1;
        std::deque<int> qb;
        qb.push_back(startV);
        
        while (!qb.empty()) {
            int v = qb.front(); qb.pop_front();
            for (int eid : v2e[v]) {
                if (edges[eid].f1 != -1) continue; // Only boundary edges
                int a = edges[eid].u, b = edges[eid].v;
                int other = (v == a) ? b : a;
                if (!visV[other] && vBC[other] == bc) { // Same boundary component
                    visV[other] = 1;
                    inT[eid] = 1;
                    qb.push_back(other);
                }
            }
        }
    }
    
    // Phase 2 & 3: Connect remaining vertices (interior) using BFS from visited vertices
    std::deque<int> qv;
    for (int v = 0; v < nV; ++v) {
        if (visV[v] && usedV[v]) qv.push_back(v);
    }

    while (!qv.empty()) {
        int v = qv.front(); qv.pop_front();
        for (int eid : v2e[v]) {
            int a = edges[eid].u, b = edges[eid].v;
            int other = (v == a) ? b : a;
            if (!visV[other]) {
                visV[other] = 1;
                inT[eid] = 1;
                qv.push_back(other);
            }
        }
    }
    
    // Check how many vertices were visited
    int visitedCount = 0;
    for (int v = 0; v < nV; ++v) if (visV[v]) visitedCount++;
    
    if (visitedCount != Vused) {
        std::cerr << "WARNING: BFS only reached " << visitedCount << "/" << Vused << " vertices!\n";
    }

    // -----------------------
    // 5) Dual spanning tree C* avoiding T edges (BFS on augmented dual)
    //    IMPORTANT: We must be careful with boundary edges. Using too many
    //    boundary edges in C* will disconnect faces from each other.
    //    Strategy: First BFS through face-to-face edges, then add minimal
    //    boundary edges to reach boundary nodes.
    // -----------------------
    std::vector<uint8_t> inC(nE, 0);
    std::vector<uint8_t> visD(nD, 0);
    std::deque<int> qd;

    // Phase 1: BFS through interior (face-to-face) edges only
    int rootD = 0;
    visD[rootD] = 1;
    qd.push_back(rootD);

    while (!qd.empty()) {
        int d = qd.front(); qd.pop_front();
        if (d >= nT) continue; // Skip if we somehow queued a boundary node
        
        for (const auto& arc : dualAdj[d]) {
            if (inT[arc.eid]) continue; // Edge-disjoint from primal tree
            if (arc.to >= nT) continue;  // Skip boundary nodes in phase 1
            if (!visD[arc.to]) {
                visD[arc.to] = 1;
                inC[arc.eid] = 1;
                qd.push_back(arc.to);
            }
        }
    }
    
    // Count how many faces were reached in phase 1
    int facesReached = 0;
    for (int f = 0; f < nT; ++f) if (visD[f]) facesReached++;
    
    // Phase 2: Add boundary edges to reach boundary nodes
    // For each boundary node, find ONE edge from a visited face
    for (int bc = 0; bc < bCount; ++bc) {
        int bdNode = nT + bc;
        if (visD[bdNode]) continue; // Already reached (shouldn't happen)
        
        // Find a face-to-boundary edge where the face is visited
        for (const auto& arc : dualAdj[bdNode]) {
            if (inT[arc.eid]) continue;
            if (arc.to < nT && visD[arc.to]) {
                // This edge connects a visited face to this boundary node
                inC[arc.eid] = 1;
                visD[bdNode] = 1;
                break;
            }
        }
    }
    
    // Check if all dual nodes were reached
    int reachedD = 0;
    for (int i = 0; i < nD; ++i) if (visD[i]) reachedD++;
    
    if (reachedD != nD) {
        std::cerr << "WARNING: Dual tree incomplete (" << reachedD << "/" << nD 
                  << "). Topology may not be disk.\n";
    }

    // -----------------------
    // 6) Cut graph = T ∪ L where L = edges not in T and not in C*
    // -----------------------
    std::vector<uint8_t> inCut(nE, 0);
    for (int eid = 0; eid < nE; ++eid) {
        bool leftover = (!inT[eid] && !inC[eid]);
        if (inT[eid] || leftover) {
            inCut[eid] = 1;
        }
    }

    // -----------------------
    // 7) Prune dangling branches (protect terminal vertices)
    // -----------------------
    std::vector<int> deg(nV, 0);
    for (int eid = 0; eid < nE; ++eid) {
        if (!inCut[eid]) continue;
        deg[edges[eid].u]++;
        deg[edges[eid].v]++;
    }

    std::deque<int> st;
    for (int v = 0; v < nV; ++v) {
        if (usedV[v] && deg[v] == 1 && !terminal[v]) {
            st.push_back(v);
        }
    }

    while (!st.empty()) {
        int v = st.front(); st.pop_front();
        if (deg[v] != 1) continue;
        if (terminal[v]) continue;

        int found = -1;
        for (int eid : v2e[v]) {
            if (inCut[eid]) { found = eid; break; }
        }
        if (found == -1) continue;

        int a = edges[found].u;
        int b = edges[found].v;
        int other = (v == a) ? b : a;

        inCut[found] = 0;
        deg[a]--;
        deg[b]--;

        if (usedV[other] && deg[other] == 1 && !terminal[other]) {
            st.push_back(other);
        }
    }
    
    // Verify cut graph doesn't disconnect the mesh
    DSU faceDsu(nT);
    for (int eid = 0; eid < nE; ++eid) {
        if (inCut[eid]) continue;
        const InternalEdge& e = edges[eid];
        if (e.f1 == -1) continue;
        faceDsu.unite(e.f0, e.f1);
    }
    std::set<int> faceComponents;
    for (int f = 0; f < nT; ++f) {
        faceComponents.insert(faceDsu.find(f));
    }
    
    if (faceComponents.size() > 1) {
        std::cerr << "WARNING: Cuts disconnect the mesh into " << faceComponents.size() << " pieces!\n";
    }
    
    // Check connectivity of cut graph via DSU on vertices
    DSU dsuCut(nV);
    for (int eid = 0; eid < nE; ++eid) {
        if (!inCut[eid]) continue;
        dsuCut.unite(edges[eid].u, edges[eid].v);
    }
    
    // Check if all boundary components the cut reaches are in the SAME connected component
    std::vector<int> bcCutRoot(bCount, -1);
    for (int v = 0; v < nV; ++v) {
        if (vBC[v] >= 0) {
            for (int eid : v2e[v]) {
                if (inCut[eid]) {
                    if (bcCutRoot[vBC[v]] == -1) {
                        bcCutRoot[vBC[v]] = dsuCut.find(v);
                    }
                    break;
                }
            }
        }
    }
    
    std::set<int> cutComponents;
    for (int bc = 0; bc < bCount; ++bc) {
        if (bcCutRoot[bc] >= 0) {
            cutComponents.insert(dsuCut.find(bcCutRoot[bc]));
        }
    }
    
    // If cut graph components are disconnected, add connecting paths
    if (cutComponents.size() > 1 && bCount > 1) {
        
        // Build adjacency for shortest paths
        std::vector<std::vector<std::pair<int, double>>> wadj(nV);
        for (int eid = 0; eid < nE; ++eid) {
            if (!usedV[edges[eid].u] || !usedV[edges[eid].v]) continue;
            double w = edgeLength(orig, edges[eid].u, edges[eid].v);
            wadj[edges[eid].u].emplace_back(edges[eid].v, w);
            wadj[edges[eid].v].emplace_back(edges[eid].u, w);
        }
        
        // Get representatives of each cut component
        std::vector<int> compReps(cutComponents.begin(), cutComponents.end());
        
        struct PQItem { double d; int v; bool operator>(const PQItem& o) const { return d > o.d; } };
        
        // Connect component 0 to all other components
        for (size_t i = 1; i < compReps.size(); ++i) {
            // Find vertices in each component
            std::vector<int> srcVerts, dstVerts;
            for (int v = 0; v < nV; ++v) {
                if (!usedV[v]) continue;
                // Check if this vertex is on the cut graph
                bool onCut = false;
                for (int eid : v2e[v]) {
                    if (inCut[eid]) { onCut = true; break; }
                }
                if (!onCut) continue;
                
                if (dsuCut.find(v) == dsuCut.find(compReps[0])) srcVerts.push_back(v);
                if (dsuCut.find(v) == dsuCut.find(compReps[i])) dstVerts.push_back(v);
            }
            
            if (srcVerts.empty() || dstVerts.empty()) continue;
            
            // Dijkstra from all srcVerts to find closest dstVert
            std::vector<double> dist(nV, std::numeric_limits<double>::infinity());
            std::vector<int> prev(nV, -1);
            std::priority_queue<PQItem, std::vector<PQItem>, std::greater<PQItem>> pq;
            
            for (int v : srcVerts) {
                dist[v] = 0.0;
                pq.push({0.0, v});
            }
            std::set<int> dstSet(dstVerts.begin(), dstVerts.end());
            
            int target = -1;
            while (!pq.empty()) {
                auto [d, u] = pq.top();
                pq.pop();
                if (d != dist[u]) continue;
                
                if (dstSet.count(u)) {
                    target = u;
                    break;
                }
                
                for (const auto& [to, w] : wadj[u]) {
                    if (dist[to] > d + w) {
                        dist[to] = d + w;
                        prev[to] = u;
                        pq.push({dist[to], to});
                    }
                }
            }
            
            if (target != -1) {
                // Add path edges to cut
                int curr = target;
                while (prev[curr] != -1) {
                    int p = prev[curr];
                    EdgeKey ek(p, curr);
                    auto it = edgeMap.find(ek);
                    if (it != edgeMap.end()) {
                        inCut[it->second] = 1;
                        dsuCut.unite(p, curr);
                    }
                    curr = p;
                }
            }
        }
    }
    
    // -----------------------
    // POST-PROCESSING: Ensure all boundary components will merge
    // -----------------------
    // For boundaries to merge, cut paths must connect them THROUGH the mesh
    // Build a DSU where boundary vertices merge if connected by cut edges OR by boundary edges
    if (bCount > 1) {
        DSU dsuBoundary(nV);
        
        // First, unite all boundary vertices on the same boundary component
        for (int eid = 0; eid < nE; ++eid) {
            if (edges[eid].f1 == -1) {  // Boundary edge
                dsuBoundary.unite(edges[eid].u, edges[eid].v);
            }
        }
        
        // Now unite vertices connected by cut edges
        for (int eid = 0; eid < nE; ++eid) {
            if (inCut[eid]) {
                dsuBoundary.unite(edges[eid].u, edges[eid].v);
            }
        }
        
        // Check how many distinct boundary "super-components" we have
        std::set<int> boundaryRoots;
        for (int bc = 0; bc < bCount; ++bc) {
            if (bcRep[bc] >= 0) {
                boundaryRoots.insert(dsuBoundary.find(bcRep[bc]));
            }
        }
        
        // If we have more than one, we need to add connecting paths
        if (boundaryRoots.size() > 1) {
            
            std::vector<std::vector<std::pair<int, double>>> wadj(nV);
            for (int eid = 0; eid < nE; ++eid) {
                if (!usedV[edges[eid].u] || !usedV[edges[eid].v]) continue;
                double w = edgeLength(orig, edges[eid].u, edges[eid].v);
                wadj[edges[eid].u].emplace_back(edges[eid].v, w);
                wadj[edges[eid].v].emplace_back(edges[eid].u, w);
            }
            
            struct PQItem2 { double d; int v; bool operator>(const PQItem2& o) const { return d > o.d; } };
            
            std::vector<int> breps(boundaryRoots.begin(), boundaryRoots.end());
            
            for (size_t i = 1; i < breps.size(); ++i) {
                // Find boundary vertices in each super-component
                std::vector<int> srcBV, dstBV;
                for (int bc = 0; bc < bCount; ++bc) {
                    if (bcRep[bc] < 0) continue;
                    if (dsuBoundary.find(bcRep[bc]) == dsuBoundary.find(breps[0])) {
                        // Add all boundary vertices from this BC
                        for (int v = 0; v < nV; ++v) {
                            if (vBC[v] == bc) srcBV.push_back(v);
                        }
                    }
                    if (dsuBoundary.find(bcRep[bc]) == dsuBoundary.find(breps[i])) {
                        for (int v = 0; v < nV; ++v) {
                            if (vBC[v] == bc) dstBV.push_back(v);
                        }
                    }
                }
                
                if (srcBV.empty() || dstBV.empty()) continue;
                
                // Dijkstra from srcBV to dstBV
                std::vector<double> dist(nV, std::numeric_limits<double>::infinity());
                std::vector<int> prev(nV, -1);
                std::priority_queue<PQItem2, std::vector<PQItem2>, std::greater<PQItem2>> pq;
                
                for (int v : srcBV) {
                    dist[v] = 0.0;
                    pq.push({0.0, v});
                }
                std::set<int> dstSet(dstBV.begin(), dstBV.end());
                
                int target = -1;
                while (!pq.empty()) {
                    auto [d, u] = pq.top();
                    pq.pop();
                    if (d != dist[u]) continue;
                    
                    if (dstSet.count(u)) {
                        target = u;
                        break;
                    }
                    
                    for (const auto& [to, w] : wadj[u]) {
                        if (dist[to] > d + w) {
                            dist[to] = d + w;
                            prev[to] = u;
                            pq.push({dist[to], to});
                        }
                    }
                }
                
                if (target != -1) {
                    int curr = target;
                    while (prev[curr] != -1) {
                        int p = prev[curr];
                        EdgeKey ek(p, curr);
                        auto it = edgeMap.find(ek);
                        if (it != edgeMap.end()) {
                            inCut[it->second] = 1;
                        }
                        dsuBoundary.unite(p, curr);
                        curr = p;
                    }
                }
            }
        }
    }

    // -----------------------
    // 8) Store cut edges
    // -----------------------
    for (int eid = 0; eid < nE; ++eid) {
        if (inCut[eid] && edges[eid].f1 != -1) { // only stores cuts not on boundary
            cutEdges.insert(EdgeKey(edges[eid].u, edges[eid].v));
        }
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

    // Target markers: boundary vertices OR vertices incident to cutEdges.
    // A singularity is considered "connected" if it lies on the boundary or
    // is an endpoint of an edge in cutEdges.
    std::vector<char> isTarget(nV, false);

    // Mark boundary vertices
    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = orig.triangles[t];
        for (int e = 0; e < 3; ++e) {
            if (orig.triangleAdjacency[t][e] != -1) continue;
            int a = tri[e];
            int b = tri[(e + 1) % 3];
            if (a >= 0 && a < nV) isTarget[a] = true;
            if (b >= 0 && b < nV) isTarget[b] = true;
        }
    }

    // Mark vertices incident to cutEdges
    for (const auto &ek : cutEdges) {
        if (ek.a >= 0 && ek.a < nV) isTarget[ek.a] = true;
        if (ek.b >= 0 && ek.b < nV) isTarget[ek.b] = true;
    }

    // Dijkstra from each singularity to the nearest target (boundary or cut edge).
    using QItem = std::pair<double, int>;
    auto runDijkstra = [&](int source, const std::vector<char> &isT, std::vector<int> &prev) -> int {
        prev.assign(nV, -1);
        std::vector<double> dist(nV, std::numeric_limits<double>::infinity());
        std::priority_queue<QItem, std::vector<QItem>, std::greater<QItem>> pq;
        dist[source] = 0.0;
        pq.push({0.0, source});
        while (!pq.empty()) {
            const auto [d, v] = pq.top();
            pq.pop();
            if (d != dist[v]) continue;
            if (isT[v] && v != source) return v;
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

        // Skip if singularity is already on boundary or connected to a cut edge
        if (isTarget[vSing]) continue;

        const int target = runDijkstra(vSing, isTarget, prev);
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
            singularityPathCutEdges.insert(ek);
            cur = p;
        }
    }

    // add all singularity path cut edges to the main cutEdges set
    for (const auto &ek : singularityPathCutEdges) {
        cutEdges.insert(ek);
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

    // Build edge table with incident faces (same approach as CutMeshDisk.cxx)
    struct InternalEdge {
        int u = -1, v = -1;   // endpoints (u < v)
        int f0 = -1, f1 = -1; // incident faces (f1 == -1 => boundary)
    };
    
    std::unordered_map<EdgeKey, int, EdgeKeyHash> edgeMap;
    edgeMap.reserve(static_cast<size_t>(nT) * 3);
    std::vector<InternalEdge> edges;
    edges.reserve(static_cast<size_t>(nT * 3 / 2));

    auto addEdge = [&](int a, int b, int fidx) {
        EdgeKey k(a, b);
        auto it = edgeMap.find(k);
        if (it == edgeMap.end()) {
            int eid = static_cast<int>(edges.size());
            edgeMap.emplace(k, eid);
            InternalEdge e;
            e.u = k.a; e.v = k.b;
            e.f0 = fidx; e.f1 = -1;
            edges.push_back(e);
        } else {
            InternalEdge& e = edges[it->second];
            if (e.f1 == -1) e.f1 = fidx;
        }
    };

    for (int f = 0; f < nT; ++f) {
        const Triangle& t = orig.triangles[f];
        addEdge(t[0], t[1], f);
        addEdge(t[1], t[2], f);
        addEdge(t[2], t[0], f);
    }

    const int nE = static_cast<int>(edges.size());

    // Helper to find local index of vertex v in triangle t
    auto localIndex = [&](int f, int v) -> int {
        const Triangle& t = orig.triangles[f];
        if (t[0] == v) return 0;
        if (t[1] == v) return 1;
        if (t[2] == v) return 2;
        return -1;
    };

    // DSU over corners (3 per face). Corners are glued across *uncut* interior edges.
    DSU dsu(nT * 3);

    for (int eid = 0; eid < nE; ++eid) {
        const InternalEdge& e = edges[eid];
        if (e.f1 == -1) continue;     // boundary edge: nothing to glue
        
        // Check if this edge is in the cut set
        EdgeKey ek(e.u, e.v);
        if (cutEdges.find(ek) != cutEdges.end()) continue;  // seam edge: do not glue

        int f0 = e.f0, f1 = e.f1;

        int i0u = localIndex(f0, e.u);
        int i1u = localIndex(f1, e.u);
        int i0v = localIndex(f0, e.v);
        int i1v = localIndex(f1, e.v);

        if (i0u < 0 || i1u < 0 || i0v < 0 || i1v < 0) continue;

        dsu.unite(3*f0 + i0u, 3*f1 + i1u);
        dsu.unite(3*f0 + i0v, 3*f1 + i1v);
    }

    // Build new vertex list (one per corner-component)
    std::unordered_map<int, int> rootToNew;
    rootToNew.reserve(static_cast<size_t>(nT) * 3);

    std::vector<int> cornerNew(3 * nT, -1);

    for (int c = 0; c < 3 * nT; ++c) {
        int r = dsu.find(c);
        auto it = rootToNew.find(r);
        if (it == rootToNew.end()) {
            int f = c / 3;
            int i = c % 3;
            int origV = orig.triangles[f][i];
            int newId = static_cast<int>(cut.vertices.size());
            rootToNew.emplace(r, newId);
            cut.vertices.push_back(orig.vertices[origV]);
            cutVertToOrig.push_back(origV);
            cornerNew[c] = newId;
        } else {
            cornerNew[c] = it->second;
        }
    }

    // Build cut faces with new vertex ids
    cut.triangles.resize(nT);
    for (int f = 0; f < nT; ++f) {
        cut.triangles[f] = { cornerNew[3*f + 0], cornerNew[3*f + 1], cornerNew[3*f + 2] };
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

void CutMesh::combFieldDirections() {
    const int nT = static_cast<int>(cut.triangles.size());
    if (nT == 0) return;

    // Perform breadth-first traversal of triangles to align directions.
    std::vector<bool> visT(nT, false); 
    std::deque<int> q;

    // start at triangle 0 and align adjacent triangles accordingly
    visT[0] = true;
    q.push_back(0);

    while (!q.empty()) {
        int curr = q.front();
        q.pop_front();
        double theta_curr = computeAngle(uField[curr]);

        // for each neighbor triangle
        for (int e = 0; e < 3; ++e) {
            int nb = cut.triangleAdjacency[curr][e];
            if (nb == -1 || visT[nb]) continue;
            visT[nb] = true;
            q.push_back(nb);

            // align nb u directions to curr
            double theta_nb = computeAngle(uField[nb]);
            int k = find_rotation_matrix(theta_nb, theta_curr);

            uField[nb] = rotateVector(uField[nb], k);
        }
    }
}

void CutMesh::makeVFieldFromUField() {
    const int nT = static_cast<int>(uField.size());
    vField.resize(nT);
    for (int t = 0; t < nT; ++t) {
        vField[t] = rotateVector(uField[t], 1); // 90 degrees CCW
    }
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
