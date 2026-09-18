#include "ATLAS/PlanarDomain.hxx"

#include <algorithm>
#include <cmath>
#include <functional>
#include <set>
#include <sstream>
#include <unordered_map>

namespace {

double angleAt(const Point &p, const Point &a, const Point &b) {
    const Point u = a - p, v = b - p;
    const double nu = normP(u), nv = normP(v);
    if (nu <= 0.0 || nv <= 0.0) return 0.0;
    double c = dotP(u, v) / (nu * nv);
    c = std::max(-1.0, std::min(1.0, c));
    return std::acos(c);
}

} // namespace

PlanarDomain::PlanarDomain(const Mesh &mesh, const Options &opts) : mesh_(mesh), opts_(opts) {
    report_.vertices = static_cast<int>(mesh_.vertices.size());
    report_.edges = static_cast<int>(mesh_.edges.size());
    report_.triangles = static_cast<int>(mesh_.triangles.size());
    report_.boundaryEdges = static_cast<int>(mesh_.boundaryEdges.size());

    const int NV = report_.vertices;
    protectedVertex.assign(NV, 0);
    cornerKind.assign(NV, CornerKind::None);
    interiorAngle.assign(NV, 0.0);
    loopOf.assign(NV, -1);
    interfaceDegree.assign(NV, 0);
    interfaceEdge.assign(mesh_.edges.size(), 0);
    boundaryEdge.assign(mesh_.edges.size(), 0);
    for (int e : mesh_.boundaryEdges) boundaryEdge[e] = 1;

    if (NV == 0 || mesh_.triangles.empty()) {
        report_.messages.push_back("Stage 1: the mesh is empty");
        return;
    }

    Point lo = mesh_.vertices[0], hi = mesh_.vertices[0];
    for (const Point &p : mesh_.vertices) {
        lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
        hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
    }
    report_.scale = normP(hi - lo);
    double edgeSum = 0.0;
    for (const auto &e : mesh_.edges) edgeSum += normP(mesh_.vertices[e[1]] - mesh_.vertices[e[0]]);
    report_.meanEdge = mesh_.edges.empty() ? 0.0 : edgeSum / mesh_.edges.size();

    checkTriangles();
    checkManifold();
    buildLoops();
    buildComponents();
    tagInterfaces();
    tagCorners();

    report_.valid = report_.degenerateTriangles == 0 && report_.invertedTriangles == 0 &&
                    report_.nonManifoldEdges == 0 && report_.nonManifoldVertices == 0 &&
                    report_.componentsWithoutOuterLoop == 0 && report_.eulerMismatches == 0 &&
                    report_.components > 0;
}

// Nondegenerate and counter-clockwise. The .obj loader already turns every face
// counter-clockwise, but a Mesh built from arrays keeps whatever it was given,
// and Sec. 3.1's determinant is only positive for the orientation it assumes.
void PlanarDomain::checkTriangles() {
    const int NT = report_.triangles;
    triangleArea.assign(NT, 0.0);
    const double tiny = opts_.degenerateArea * report_.scale * report_.scale;
    double minArea = 1e300, minAngle = 1e300;
    for (int t = 0; t < NT; ++t) {
        const Triangle &tri = mesh_.triangles[t];
        const Point &a = mesh_.vertices[tri[0]], &b = mesh_.vertices[tri[1]], &c = mesh_.vertices[tri[2]];
        const double A = 0.5 * cross2(b - a, c - a);
        triangleArea[t] = A;
        report_.area += A;
        minArea = std::min(minArea, A);
        if (std::fabs(A) <= tiny) ++report_.degenerateTriangles;
        else if (A < 0.0) ++report_.invertedTriangles;
        for (int i = 0; i < 3; ++i) {
            const int v = tri[i];
            const double ang = angleAt(mesh_.vertices[v], mesh_.vertices[tri[(i + 1) % 3]],
                                       mesh_.vertices[tri[(i + 2) % 3]]);
            interiorAngle[v] += ang;
            minAngle = std::min(minAngle, ang);
        }
    }
    report_.minTriangleArea = minArea;
    report_.minAngleDegrees = minAngle * 180.0 / M_PI;
    if (report_.degenerateTriangles > 0) {
        report_.messages.push_back("Stage 1: " + std::to_string(report_.degenerateTriangles) +
                                   " degenerate triangle(s); Sec. 3.1's bound needs det > 0");
    }
    if (report_.invertedTriangles > 0) {
        report_.messages.push_back("Stage 1: " + std::to_string(report_.invertedTriangles) +
                                   " clockwise triangle(s)");
    }
}

// Mesh's constructor keeps at most two triangles per edge and silently
// overwrites the adjacency of a third, so the edge census is redone here from
// the triangles themselves. The vertex test is the fan: the triangles around a
// manifold vertex are one edge-connected fan (a cycle inside, a path on dS).
void PlanarDomain::checkManifold() {
    std::unordered_map<MeshEdgeKey, int, MeshEdgeKeyHash> uses;
    uses.reserve(mesh_.triangles.size() * 3);
    for (const Triangle &t : mesh_.triangles) {
        for (int i = 0; i < 3; ++i) ++uses[MeshEdgeKey(t[i], t[(i + 1) % 3])];
    }
    for (const auto &kv : uses) {
        if (kv.second > 2) ++report_.nonManifoldEdges;
    }

    const int NV = report_.vertices;
    std::vector<std::vector<int>> vt(NV);
    for (int t = 0; t < report_.triangles; ++t) {
        for (int i = 0; i < 3; ++i) vt[mesh_.triangles[t][i]].push_back(t);
    }
    std::vector<int> boundaryEdgesAt(NV, 0);
    for (int e : mesh_.boundaryEdges) {
        ++boundaryEdgesAt[mesh_.edges[e][0]];
        ++boundaryEdgesAt[mesh_.edges[e][1]];
    }
    for (int v = 0; v < NV; ++v) {
        const std::vector<int> &fan = vt[v];
        if (fan.empty()) continue;   // an unused vertex is harmless
        if (boundaryEdgesAt[v] > 2) { ++report_.nonManifoldVertices; continue; }
        // Edge-connected components of the fan, walking triangleAdjacency
        // across the edges that contain v.
        std::vector<char> seen(fan.size(), 0);
        int pieces = 0;
        for (size_t s = 0; s < fan.size(); ++s) {
            if (seen[s]) continue;
            ++pieces;
            std::vector<size_t> stack{s};
            seen[s] = 1;
            while (!stack.empty()) {
                const size_t k = stack.back();
                stack.pop_back();
                const int t = fan[k];
                for (int le = 0; le < 3; ++le) {
                    const int a = mesh_.triangles[t][le], b = mesh_.triangles[t][(le + 1) % 3];
                    if (a != v && b != v) continue;
                    const int n = mesh_.triangleAdjacency[t][le];
                    if (n < 0) continue;
                    for (size_t j = 0; j < fan.size(); ++j) {
                        if (!seen[j] && fan[j] == n) { seen[j] = 1; stack.push_back(j); }
                    }
                }
            }
        }
        if (pieces != 1) ++report_.nonManifoldVertices;
    }
    if (report_.nonManifoldEdges > 0) {
        report_.messages.push_back("Stage 1: " + std::to_string(report_.nonManifoldEdges) +
                                   " edge(s) shared by three or more triangles");
    }
    if (report_.nonManifoldVertices > 0) {
        report_.messages.push_back("Stage 1: " + std::to_string(report_.nonManifoldVertices) +
                                   " pinched or bow-tie vertex/vertices; Sec. 1.1 requires "
                                   "preprocessing or rejection");
    }
}

// Boundary edges directed as they occur in their (counter-clockwise) triangle,
// which puts the domain on the left, then chained end to start.
void PlanarDomain::buildLoops() {
    const int NV = report_.vertices;
    std::vector<int> outNext(NV, -1), outEdge(NV, -1);
    for (int t = 0; t < report_.triangles; ++t) {
        for (int le = 0; le < 3; ++le) {
            if (mesh_.triangleAdjacency[t][le] >= 0) continue;
            const int a = mesh_.triangles[t][le], b = mesh_.triangles[t][(le + 1) % 3];
            if (outNext[a] >= 0) continue;   // a pinch; already counted above
            outNext[a] = b;
            outEdge[a] = mesh_.triangleEdges[t][le];
        }
    }
    std::vector<char> used(NV, 0);
    for (int s = 0; s < NV; ++s) {
        if (outNext[s] < 0 || used[s]) continue;
        Loop L;
        int v = s;
        while (v >= 0 && !used[v]) {
            used[v] = 1;
            L.vertices.push_back(v);
            L.edges.push_back(outEdge[v]);
            v = outNext[v];
        }
        if (v != s) {
            report_.messages.push_back("Stage 1: a boundary chain does not close into a loop");
            ++report_.nonManifoldVertices;
            continue;
        }
        const int n = static_cast<int>(L.vertices.size());
        for (int i = 0; i < n; ++i) {
            const Point &p = mesh_.vertices[L.vertices[i]];
            const Point &q = mesh_.vertices[L.vertices[(i + 1) % n]];
            L.signedArea += 0.5 * cross2(p, q);
            L.length += normP(q - p);
        }
        L.outer = L.signedArea > 0.0;
        const int li = static_cast<int>(loops.size());
        for (int u : L.vertices) loopOf[u] = li;
        report_.boundaryLength += L.length;
        loops.push_back(std::move(L));
    }
    report_.loops = static_cast<int>(loops.size());
}

// Connected components by triangle adjacency, each with its loops. Sec. 4's
// identity is per component, with h its number of holes, and V - E + F = 1 - h
// is checked here so that h is known from the loops rather than assumed.
void PlanarDomain::buildComponents() {
    const int NT = report_.triangles;
    triangleComponent.assign(NT, -1);
    for (int seed = 0; seed < NT; ++seed) {
        if (triangleComponent[seed] >= 0) continue;
        const int c = static_cast<int>(components.size());
        components.emplace_back();
        std::vector<int> stack{seed};
        triangleComponent[seed] = c;
        while (!stack.empty()) {
            const int t = stack.back();
            stack.pop_back();
            components[c].triangles.push_back(t);
            for (int k = 0; k < 3; ++k) {
                const int n = mesh_.triangleAdjacency[t][k];
                if (n >= 0 && triangleComponent[n] < 0) {
                    triangleComponent[n] = c;
                    stack.push_back(n);
                }
            }
        }
    }
    report_.components = static_cast<int>(components.size());

    for (int li = 0; li < static_cast<int>(loops.size()); ++li) {
        Loop &L = loops[li];
        // The component of the triangle owning the loop's first edge.
        const int e = L.edges.empty() ? -1 : L.edges[0];
        const int t = (e >= 0) ? mesh_.edgeTriangles[e][0] : -1;
        L.component = (t >= 0) ? triangleComponent[t] : -1;
        if (L.component < 0) continue;
        Component &C = components[L.component];
        if (L.outer) {
            if (C.outerLoop >= 0) {
                report_.messages.push_back("Stage 1: a component has two counter-clockwise loops");
                ++report_.componentsWithoutOuterLoop;
            }
            C.outerLoop = li;
        } else {
            C.holes.push_back(li);
            ++report_.holes;
        }
    }

    for (Component &C : components) {
        if (C.outerLoop < 0) ++report_.componentsWithoutOuterLoop;
        std::set<int> vs, es;
        for (int t : C.triangles) {
            for (int k = 0; k < 3; ++k) {
                vs.insert(mesh_.triangles[t][k]);
                es.insert(mesh_.triangleEdges[t][k]);
            }
        }
        C.vertices = static_cast<int>(vs.size());
        C.edges = static_cast<int>(es.size());
        C.eulerCharacteristic = C.vertices - C.edges + static_cast<int>(C.triangles.size());
        if (C.eulerCharacteristic != 1 - static_cast<int>(C.holes.size())) {
            ++report_.eulerMismatches;
            std::ostringstream os;
            os << "Stage 1: a component has V - E + F = " << C.eulerCharacteristic << " but "
               << C.holes.size() << " hole(s)";
            report_.messages.push_back(os.str());
        }
    }
    if (report_.componentsWithoutOuterLoop > 0) {
        report_.messages.push_back("Stage 1: a component has no single outer loop");
    }
}

// The material interface network. Sec. 1.1: interfaces are protected polylines
// the input mesh already embeds (each is a chain of mesh edges whose two
// triangles disagree on the material), and their nodes -- where the network
// branches, lands on dS, kinks or dangles -- are required macrovertices.
void PlanarDomain::tagInterfaces() {
    std::set<int> mats(mesh_.triangleMatId.begin(), mesh_.triangleMatId.end());
    report_.materials = static_cast<int>(mats.size());
    if (!opts_.materialInterfaces) return;
    for (int e = 0; e < report_.edges; ++e) {
        const int t0 = mesh_.edgeTriangles[e][0], t1 = mesh_.edgeTriangles[e][1];
        if (t0 < 0 || t1 < 0) continue;
        if (mesh_.triangleMatId[t0] == mesh_.triangleMatId[t1]) continue;
        interfaceEdge[e] = 1;
        ++report_.interfaceEdges;
        ++interfaceDegree[mesh_.edges[e][0]];
        ++interfaceDegree[mesh_.edges[e][1]];
    }
}

void PlanarDomain::tagCorners() {
    const int NV = report_.vertices;
    const double cornerTol = opts_.cornerAngle * M_PI / 180.0;
    const double kinkTol = opts_.kinkAngle * M_PI / 180.0;

    auto mark = [&](int v, CornerKind k) {
        if (cornerKind[v] == CornerKind::None) cornerKind[v] = k;
        protectedVertex[v] = 1;
    };

    for (int v = 0; v < NV; ++v) {
        const bool onBoundary = loopOf[v] >= 0;
        const int deg = interfaceDegree[v];
        if (deg > 0) {
            if (onBoundary) mark(v, CornerKind::InterfaceLanding);
            else if (deg >= 3) mark(v, CornerKind::InterfaceJunction);
            else if (deg == 1) mark(v, CornerKind::InterfaceDangling);
        }
        if (onBoundary && std::fabs(interiorAngle[v] - M_PI) > cornerTol) {
            mark(v, CornerKind::BoundaryCorner);
        }
    }

    // Kinks: an interior vertex with exactly two interface edges, turning.
    if (opts_.materialInterfaces && report_.interfaceEdges > 0) {
        std::vector<std::vector<int>> nbr(NV);
        for (int e = 0; e < report_.edges; ++e) {
            if (!interfaceEdge[e]) continue;
            nbr[mesh_.edges[e][0]].push_back(mesh_.edges[e][1]);
            nbr[mesh_.edges[e][1]].push_back(mesh_.edges[e][0]);
        }
        for (int v = 0; v < NV; ++v) {
            if (nbr[v].size() != 2 || loopOf[v] >= 0) continue;
            const double a = angleAt(mesh_.vertices[v], mesh_.vertices[nbr[v][0]],
                                     mesh_.vertices[nbr[v][1]]);
            if (std::fabs(a - M_PI) > kinkTol) mark(v, CornerKind::InterfaceKink);
        }
    }

    for (int v : opts_.extraCorners) {
        if (v >= 0 && v < NV) mark(v, CornerKind::UserDesignated);
    }

    for (int v = 0; v < NV; ++v) {
        switch (cornerKind[v]) {
            case CornerKind::BoundaryCorner: ++report_.boundaryCorners; break;
            case CornerKind::InterfaceJunction: ++report_.interfaceJunctions; break;
            case CornerKind::InterfaceLanding: ++report_.interfaceLandings; break;
            case CornerKind::InterfaceKink: ++report_.interfaceKinks; break;
            case CornerKind::InterfaceDangling: ++report_.interfaceDangling; break;
            case CornerKind::UserDesignated: ++report_.userCorners; break;
            default: break;
        }
        if (protectedVertex[v]) {
            ++report_.protectedVertices;
            if (loopOf[v] >= 0) ++loops[loopOf[v]].corners;
        }
    }
    if (report_.interfaceDangling > 0) {
        report_.messages.push_back("Stage 1: " + std::to_string(report_.interfaceDangling) +
                                   " interface end(s) dangle inside a material; kept as corners");
    }
}
