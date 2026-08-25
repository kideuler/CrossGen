#include "ConeCut.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace {

inline double edgeLength(const Mesh &m, int a, int b) {
    const Point &pa = m.vertices[a];
    const Point &pb = m.vertices[b];
    const double dx = pa[0] - pb[0];
    const double dy = pa[1] - pb[1];
    return std::sqrt(dx * dx + dy * dy);
}

// Disjoint-set union over triangle corners, as in
// HarmonicCut::buildExplicitCutMesh(): corners are glued across every interior
// edge that is not in G, and each surviving group becomes one vertex of Omega.
struct DSU {
    std::vector<int> parent;
    std::vector<int> rank;
    explicit DSU(int n = 0) : parent(n), rank(n, 0) {
        std::iota(parent.begin(), parent.end(), 0);
    }
    int find(int x) {
        int p = x;
        while (parent[p] != p) p = parent[p];
        while (parent[x] != x) { const int nx = parent[x]; parent[x] = p; x = nx; }
        return p;
    }
    void unite(int a, int b) {
        a = find(a); b = find(b);
        if (a == b) return;
        if (rank[a] < rank[b]) std::swap(a, b);
        parent[b] = a;
        if (rank[a] == rank[b]) rank[a]++;
    }
};

} // namespace

ConeCut::ConeCut(std::shared_ptr<Mesh> mesh, const ConeSingularities &cones)
    : orig(std::move(mesh)) {
    if (!orig) throw std::runtime_error("ConeCut: null mesh");
    if (orig->triangles.empty()) throw std::runtime_error("ConeCut: empty mesh");
    if (orig.get() != &cones.getMesh()) {
        throw std::runtime_error("ConeCut: the cone set was measured on a different mesh");
    }

    // --- The voids, by Wang et al. Sec. 4.1 -------------------------------
    harmonic = std::make_unique<HarmonicCut>(orig);
    cutEdges = harmonic->getCutEdges();
    report.voids = harmonic->getReport().voids;
    report.harmonicCuts = harmonic->getReport().cutsMade;
    for (const std::string &m : harmonic->getReport().messages) {
        report.messages.push_back("HarmonicCut: " + m);
    }

    // --- The cones, by Shepherd et al. Sec. 3.2.2 --------------------------
    routeConePaths(cones);

    buildExplicitCutMesh();
    check(cones);
}

// ---------------------------------------------------------------------------
// shortestPathToSet()
// ---------------------------------------------------------------------------
std::vector<int> ConeCut::shortestPathToSet(int source,
                                            const std::vector<char> &targets,
                                            const std::vector<char> &forbidden) const {
    const int nV = static_cast<int>(orig->vertices.size());
    std::vector<double> dist(nV, std::numeric_limits<double>::infinity());
    std::vector<int> prev(nV, -1);

    using QItem = std::pair<double, int>;
    std::priority_queue<QItem, std::vector<QItem>, std::greater<QItem>> pq;
    dist[source] = 0.0;
    pq.push({0.0, source});

    int reached = -1;
    while (!pq.empty()) {
        const auto [d, u] = pq.top();
        pq.pop();
        if (d != dist[u]) continue;

        // The first target reached is where the arc ends -- unless it is
        // another cone, which is neither a place to stop nor a place to pass
        // through. Sec. 3.2.2 wants a small neighbourhood of every singular
        // point to stay one connected component after cutting, and an arc that
        // *ends* on a cone breaks that just as thoroughly as one that crosses
        // it: the cone then has two arcs at it and the graph runs through
        // rather than to it. Earlier cone arcs make all their vertices targets,
        // so the arc still terminates on the same arc, one vertex further
        // along.
        if (u != source && forbidden[u]) continue;
        if (u != source && targets[u]) { reached = u; break; }

        // Neighbours of u, read off the incident triangles.
        const auto &vt = orig->vertexTriangles;
        for (int i = vt.rowPtr[u]; i < vt.rowPtr[u + 1]; ++i) {
            const Triangle &tri = orig->triangles[vt.colIdx[i]];
            for (int c = 0; c < 3; ++c) {
                const int to = tri[c];
                if (to == u) continue;
                const double w = edgeLength(*orig, u, to);
                if (dist[to] > d + w) {
                    dist[to] = d + w;
                    prev[to] = u;
                    pq.push({dist[to], to});
                }
            }
        }
    }

    std::vector<int> path;
    if (reached < 0) return path;
    for (int v = reached; v != -1; v = prev[v]) path.push_back(v);
    std::reverse(path.begin(), path.end()); // source (the cone) first
    return path;
}

// ---------------------------------------------------------------------------
// routeConePaths()  --  Sec. 3.2.2
//
//   targets <- dS
//   sort interior cones by graph distance to dS, ascending
//   for each interior cone c:
//       p <- shortest path from c to (targets union G), weight l_ij,
//            forbidding traversal *through* any other cone
//       G <- G union p
//
// The ordering matters for cut length, not for correctness: taking the cone
// nearest the boundary first lets the ones behind it terminate on the arc it
// just laid rather than running all the way out on their own.
// ---------------------------------------------------------------------------
void ConeCut::routeConePaths(const ConeSingularities &cones) {
    const int nV = static_cast<int>(orig->vertices.size());

    std::vector<char> isTarget(nV, 0);
    for (int v : orig->boundaryVertices) isTarget[v] = 1;
    for (const EdgeKey &e : cutEdges) { isTarget[e.a] = 1; isTarget[e.b] = 1; }

    // Boundary cones count as off-limits too. One already sits in dS and needs
    // no arc, but an arc *ending* on it splits its one-ring exactly as it would
    // an interior cone's, and Sec. 3.2.2 rules that out for every point of P.
    std::vector<char> isCone(nV, 0);
    for (const auto &c : cones.getCones()) isCone[c.vertex] = 1;

    std::vector<ConeSingularities::Cone> interior = cones.interiorCones();
    report.interiorCones = static_cast<int>(interior.size());
    if (interior.empty()) return;

    // Distance from every vertex to the current graph, for the ordering only.
    std::vector<double> toGraph(nV, std::numeric_limits<double>::infinity());
    {
        using QItem = std::pair<double, int>;
        std::priority_queue<QItem, std::vector<QItem>, std::greater<QItem>> pq;
        for (int v = 0; v < nV; ++v) if (isTarget[v]) { toGraph[v] = 0.0; pq.push({0.0, v}); }
        const auto &vt = orig->vertexTriangles;
        while (!pq.empty()) {
            const auto [d, u] = pq.top();
            pq.pop();
            if (d != toGraph[u]) continue;
            for (int i = vt.rowPtr[u]; i < vt.rowPtr[u + 1]; ++i) {
                const Triangle &tri = orig->triangles[vt.colIdx[i]];
                for (int c = 0; c < 3; ++c) {
                    const int to = tri[c];
                    if (to == u) continue;
                    const double w = d + edgeLength(*orig, u, to);
                    if (toGraph[to] > w) { toGraph[to] = w; pq.push({w, to}); }
                }
            }
        }
    }

    std::sort(interior.begin(), interior.end(),
              [&](const ConeSingularities::Cone &a, const ConeSingularities::Cone &b) {
                  if (toGraph[a.vertex] != toGraph[b.vertex])
                      return toGraph[a.vertex] < toGraph[b.vertex];
                  return a.vertex < b.vertex;
              });

    for (const auto &c : interior) {
        if (isTarget[c.vertex]) {
            // Already on the graph -- a cone sitting on a void arc. Sec. 3.2.2
            // forbids this, and check() reports it; nothing to route.
            continue;
        }

        // Every other cone -- interior or boundary -- is off limits both as a
        // waypoint and as a terminus; only the cone being routed is exempt.
        std::vector<char> forbidden = isCone;
        forbidden[c.vertex] = 0;

        std::vector<int> path = shortestPathToSet(c.vertex, isTarget, forbidden);
        if (path.size() < 2) {
            std::ostringstream oss;
            oss << "Cone at vertex " << c.vertex
                << " could not be routed to the cutting graph without passing through "
                << "another cone.";
            report.messages.push_back(oss.str());
            continue;
        }

        ConePath cp;
        cp.cone = c.vertex;
        cp.index = c.index;
        cp.path = path;
        for (size_t k = 0; k + 1 < path.size(); ++k) {
            const EdgeKey ek(path[k], path[k + 1]);
            cp.edges.push_back(ek);
            cp.length += edgeLength(*orig, path[k], path[k + 1]);
            cutEdges.insert(ek);
        }
        cp.reachedBoundary = orig->isBoundaryVertex[path.back()];

        // Everything the arc touches is graph from now on, so a later cone can
        // stop on it. The cone itself is included: a subsequent path must not
        // end on it either, which keeps every cone a leaf of G.
        for (int v : path) isTarget[v] = 1;

        conePaths.push_back(std::move(cp));
        ++report.conesRouted;
    }
}

// ---------------------------------------------------------------------------
// buildExplicitCutMesh()
//
// Mirrors HarmonicCut::buildExplicitCutMesh(); it is repeated here rather than
// reused because the edge set being cut along is the union of both stages'
// arcs, which HarmonicCut never sees.
// ---------------------------------------------------------------------------
void ConeCut::buildExplicitCutMesh() {
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

    cut = Mesh(verts, tris, orig->triangleMatId);

    for (int cv = 0; cv < static_cast<int>(cutVertToOrig.size()); ++cv) {
        const int ov = cutVertToOrig[cv];
        if (ov >= 0 && ov < nV) origToCutVerts[ov].push_back(cv);
    }
}

// ---------------------------------------------------------------------------
// check()
//
// Definition 2.1 asks for two things of the cut: S - G is a topological disk,
// and P is contained in G union dS. Both are measured on Omega itself rather
// than inferred from the number of arcs laid, because an arc that failed to
// find a route leaves the mesh looking almost right.
// ---------------------------------------------------------------------------
void ConeCut::check(const ConeSingularities &cones) {
    const int V = static_cast<int>(cut.vertices.size());
    const int F = static_cast<int>(cut.triangles.size());
    const int E = static_cast<int>(cut.edges.size());
    if (V == 0 || F == 0) {
        report.messages.push_back("Cut mesh is empty.");
        return;
    }
    report.eulerCharacteristic = V - E + F;

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

    std::queue<int> q;
    std::vector<char> seen(V, 0);
    for (int v = 0; v < V; ++v) {
        if (!isBoundaryV[v] || seen[v]) continue;
        ++report.boundaryComponents;
        seen[v] = 1;
        q.push(v);
        while (!q.empty()) {
            const int x = q.front();
            q.pop();
            for (int y : bAdj[x]) if (!seen[y]) { seen[y] = 1; q.push(y); }
        }
    }

    std::vector<char> seenT(F, 0);
    for (int t = 0; t < F; ++t) {
        if (seenT[t]) continue;
        ++report.triangleComponents;
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
    report.trianglesConnected = (report.triangleComponents == 1);
    report.isDisk = report.trianglesConnected && report.boundaryComponents == 1 &&
                    report.eulerCharacteristic == 1;

    // P subset of G union dS: in Omega every child of a cone must be a boundary
    // vertex. A cone still carrying an interior child was not reached by any
    // arc.
    report.allConesOnBoundary = true;
    for (const auto &c : cones.getCones()) {
        for (int cv : origToCutVerts[c.vertex]) {
            if (!cut.isBoundaryVertex[cv]) {
                report.allConesOnBoundary = false;
                std::ostringstream oss;
                oss << "Cone at vertex " << c.vertex << " (index " << c.index
                    << ") is still interior to the cut mesh.";
                report.messages.push_back(oss.str());
                break;
            }
        }
    }

    // Sec. 3.2.2's "to but not through" preference: a cone with more than one
    // child had the graph pass across it, so its one-ring is in pieces. Legal
    // -- Fig. 7 shows a cone in three pieces and Eq. (1) reads the angle sum
    // over all of them -- but the seam bookkeeping in Stage 4 is simpler
    // without it, and a cone arc routed here never does it. What can is a
    // HarmonicCut void arc, which chooses its endpoints by proximity between
    // boundary loops and has no notion of P.
    for (const auto &c : cones.getCones()) {
        if (origToCutVerts[c.vertex].size() > 1) {
            ++report.conesSplitByCut;
            std::ostringstream oss;
            oss << "The cutting graph runs across the " << (c.onBoundary ? "boundary" : "interior")
                << " cone at vertex " << c.vertex << " rather than stopping at it, splitting it "
                << "into " << origToCutVerts[c.vertex].size()
                << " children in Omega (Sec. 3.2.2 prefers, but does not require, otherwise).";
            report.messages.push_back(oss.str());
        }
    }

    if (report.conesRouted != report.interiorCones) {
        std::ostringstream oss;
        oss << "Routed " << report.conesRouted << " of " << report.interiorCones
            << " interior cone(s) to the cutting graph.";
        report.messages.push_back(oss.str());
    }
    if (!report.trianglesConnected) {
        std::ostringstream oss;
        oss << "Cut mesh triangles fall into " << report.triangleComponents << " components.";
        report.messages.push_back(oss.str());
    }
    if (report.boundaryComponents != 1) {
        std::ostringstream oss;
        oss << "Expected 1 boundary component for a disk, got " << report.boundaryComponents << ".";
        report.messages.push_back(oss.str());
    }
    if (report.eulerCharacteristic != 1) {
        std::ostringstream oss;
        oss << "Expected Euler characteristic 1 for a disk, got "
            << report.eulerCharacteristic << ".";
        report.messages.push_back(oss.str());
    }
}

bool ConeCut::writeOBJ(const std::string &filename) const {
    std::ofstream out(filename);
    if (!out) return false;
    for (const auto &p : cut.vertices) out << "v " << p[0] << " " << p[1] << " 0\n";
    for (const auto &t : cut.triangles) {
        out << "f " << (t[0] + 1) << " " << (t[1] + 1) << " " << (t[2] + 1) << "\n";
    }
    return true;
}
