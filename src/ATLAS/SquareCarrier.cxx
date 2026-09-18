#include "ATLAS/SquareCarrier.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <map>
#include <numeric>
#include <sstream>
#include <unordered_map>

namespace {

inline long long edgeKey(int a, int b) {
    const long long lo = std::min(a, b), hi = std::max(a, b);
    return (lo << 32) | hi;
}

double angleBetween(const Point &u, const Point &v) {
    const double nu = normP(u), nv = normP(v);
    if (nu <= 0.0 || nv <= 0.0) return 0.0;
    // atan2 of cross and dot is accurate at every angle, unlike acos near 0 and pi.
    return std::atan2(cross2(u, v), dotP(u, v));
}

} // namespace

// ---------------------------------------------------------------------------
// Construction: Sec. 3, one pass over the triangles.
// ---------------------------------------------------------------------------
SquareCarrier::SquareCarrier(const PlanarDomain &domain, const Options &opts)
    : domain_(&domain), opts_(opts) {
    const Mesh &mesh = domain.getMesh();
    const int NV = static_cast<int>(mesh.vertices.size());
    const int NE = static_cast<int>(mesh.edges.size());
    const int NT = static_cast<int>(mesh.triangles.size());

    // Vertex ids: input vertices keep theirs, then one midpoint per mesh edge,
    // then one centroid per triangle. Shared midpoints are what make the split
    // conforming (Sec. 3: "constructed once per global edge").
    vertices.resize(NV + NE + NT);
    vertexOrigin.resize(NV + NE + NT);
    sourceVertex.assign(NV + NE + NT, -1);
    sourceEdge.assign(NV + NE + NT, -1);
    protectedVertex.assign(NV + NE + NT, 0);
    designatedVertex.assign(NV + NE + NT, 0);
    for (int v = 0; v < NV; ++v) {
        vertices[v] = mesh.vertices[v];
        vertexOrigin[v] = Origin::MeshVertex;
        sourceVertex[v] = v;
        protectedVertex[v] = domain.protectedVertex[v];
    }
    for (int e = 0; e < NE; ++e) {
        vertices[NV + e] = (mesh.vertices[mesh.edges[e][0]] + mesh.vertices[mesh.edges[e][1]]) * 0.5;
        vertexOrigin[NV + e] = Origin::EdgeMidpoint;
        sourceEdge[NV + e] = e;
    }
    for (int t = 0; t < NT; ++t) {
        const Triangle &tri = mesh.triangles[t];
        vertices[NV + NE + t] =
            (mesh.vertices[tri[0]] + mesh.vertices[tri[1]] + mesh.vertices[tri[2]]) / 3.0;
        vertexOrigin[NV + NE + t] = Origin::Centroid;
    }

    cells.reserve(3 * NT);
    for (int t = 0; t < NT; ++t) {
        const Triangle &tri = mesh.triangles[t];
        // triangleEdges[t][i] is the edge (tri[i], tri[i+1]) -- see Mesh's
        // constructor -- so m_{i,i+1} = NV + triangleEdges[t][i].
        const int c = NV + NE + t;
        for (int i = 0; i < 3; ++i) {
            const int vi = tri[i];
            const int mNext = NV + mesh.triangleEdges[t][i];
            const int mPrev = NV + mesh.triangleEdges[t][(i + 2) % 3];
            cells.push_back({vi, mNext, c, mPrev});
            cellMaterial.push_back(mesh.triangleMatId[t]);
            cellTriangle.push_back(t);
            cellOrigin.push_back(CellOrigin::Split);
            cellGroup.push_back(-1);
        }
    }
    initialSplit_ = true;
    rebuild();
}

// ---------------------------------------------------------------------------
// Topology
// ---------------------------------------------------------------------------
void SquareCarrier::rebuild() {
    // Compact: drop every vertex no cell references. A replacement removes the
    // interior vertices of its cavity by removing the cells around them.
    {
        const int NV = numVertices();
        std::vector<int> remap(NV, -1);
        std::vector<char> used(NV, 0);
        for (const auto &q : cells) for (int v : q) used[v] = 1;
        int next = 0;
        for (int v = 0; v < NV; ++v) if (used[v]) remap[v] = next++;
        if (next != NV) {
            auto compact = [&](auto &arr) {
                for (int v = 0; v < NV; ++v) if (remap[v] >= 0) arr[remap[v]] = arr[v];
                arr.resize(next);
            };
            compact(vertices);
            compact(vertexOrigin);
            compact(sourceVertex);
            compact(sourceEdge);
            compact(protectedVertex);
            compact(designatedVertex);
            for (auto &q : cells) for (int &v : q) v = remap[v];
        }
    }

    const int NV = numVertices();
    const int NQ = numCells();

    edges.clear();
    edgeCell.clear();
    edgeSide.clear();
    cellEdges.assign(NQ, {-1, -1, -1, -1});
    neighbor.assign(NQ, {-1, -1, -1, -1});
    neighborSide.assign(NQ, {-1, -1, -1, -1});
    transport.assign(NQ, {});
    report_.nonManifoldEdges = 0;

    std::unordered_map<long long, int> index;
    index.reserve(static_cast<size_t>(NQ) * 2 + 16);
    for (int q = 0; q < NQ; ++q) {
        for (int i = 0; i < 4; ++i) {
            const int a = cells[q][i], b = cells[q][(i + 1) & 3];
            const long long k = edgeKey(a, b);
            auto it = index.find(k);
            if (it == index.end()) {
                const int e = numEdges();
                index.emplace(k, e);
                edges.push_back({std::min(a, b), std::max(a, b)});
                edgeCell.push_back({q, -1});
                edgeSide.push_back({i, -1});
                cellEdges[q][i] = e;
            } else {
                const int e = it->second;
                cellEdges[q][i] = e;
                if (edgeCell[e][1] < 0) {
                    edgeCell[e][1] = q;
                    edgeSide[e][1] = i;
                } else {
                    ++report_.nonManifoldEdges;
                }
            }
        }
    }

    const int NE = numEdges();
    boundaryEdge.assign(NE, 0);
    interfaceEdge.assign(NE, 0);
    boundaryVertex.assign(NV, 0);
    interfaceDegree.assign(NV, 0);
    const bool useInterfaces = domain_->getOptions().materialInterfaces;
    for (int e = 0; e < NE; ++e) {
        const int q0 = edgeCell[e][0], q1 = edgeCell[e][1];
        if (q1 < 0) {
            boundaryEdge[e] = 1;
            boundaryVertex[edges[e][0]] = 1;
            boundaryVertex[edges[e][1]] = 1;
            continue;
        }
        const int i0 = edgeSide[e][0], i1 = edgeSide[e][1];
        neighbor[q0][i0] = q1;
        neighborSide[q0][i0] = i1;
        neighbor[q1][i1] = q0;
        neighborSide[q1][i1] = i0;
        transport[q0][i0] = transportAcross(i0, i1);
        transport[q1][i1] = transportAcross(i1, i0);
        if (useInterfaces && cellMaterial[q0] != cellMaterial[q1]) {
            interfaceEdge[e] = 1;
            ++interfaceDegree[edges[e][0]];
            ++interfaceDegree[edges[e][1]];
        }
    }

    buildRings();
    buildComponents();
}

// The counter-clockwise ring of each vertex. A cell occupies the sector at its
// corner c from the direction of side c (towards corner c+1) round to the
// direction of side c-1 (towards corner c-1), so the next cell
// counter-clockwise is the neighbour across side c-1.
void SquareCarrier::buildRings() {
    const int NV = numVertices();
    const int NQ = numCells();
    valence.assign(NV, 0);
    for (const auto &q : cells) for (int v : q) ++valence[v];
    ringPtr.assign(NV + 1, 0);
    for (int v = 0; v < NV; ++v) ringPtr[v + 1] = ringPtr[v] + valence[v];
    ringCell.assign(ringPtr[NV], -1);
    ringCorner.assign(ringPtr[NV], -1);

    // Unordered first, then walk.
    std::vector<int> fill(ringPtr.begin(), ringPtr.end() - 1);
    for (int q = 0; q < NQ; ++q) {
        for (int c = 0; c < 4; ++c) {
            const int v = cells[q][c];
            ringCell[fill[v]] = q;
            ringCorner[fill[v]] = c;
            ++fill[v];
        }
    }

    manifoldVertex.assign(NV, 1);
    report_.nonManifoldVertices = 0;
    std::vector<int> walkCell, walkCorner;
    for (int v = 0; v < NV; ++v) {
        const int b = ringPtr[v], n = valence[v];
        if (n == 0) continue;
        // Start cell: inside, any; on dS, the one whose side leaving v is a
        // boundary edge (the clockwise-most).
        int start = 0;
        if (boundaryVertex[v]) {
            start = -1;
            for (int k = 0; k < n; ++k) {
                const int q = ringCell[b + k], c = ringCorner[b + k];
                if (neighbor[q][c] < 0) { start = k; break; }
            }
            if (start < 0) { manifoldVertex[v] = 0; ++report_.nonManifoldVertices; continue; }
        }
        walkCell.clear();
        walkCorner.clear();
        int q = ringCell[b + start], c = ringCorner[b + start];
        for (int step = 0; step < n; ++step) {
            walkCell.push_back(q);
            walkCorner.push_back(c);
            const int side = (c + 3) & 3;
            const int r = neighbor[q][side];
            if (r < 0) break;
            const int rc = cornerOf(r, v);
            if (rc < 0) break;
            q = r;
            c = rc;
            if (q == ringCell[b + start] && step + 1 < n) break;
        }
        if (static_cast<int>(walkCell.size()) != n) {
            manifoldVertex[v] = 0;
            ++report_.nonManifoldVertices;
            continue;
        }
        for (int k = 0; k < n; ++k) {
            ringCell[b + k] = walkCell[k];
            ringCorner[b + k] = walkCorner[k];
        }
    }
}

void SquareCarrier::buildComponents() {
    const int NQ = numCells();
    cellComponent.assign(NQ, -1);
    componentCount = 0;
    for (int s = 0; s < NQ; ++s) {
        if (cellComponent[s] >= 0) continue;
        std::vector<int> stack{s};
        cellComponent[s] = componentCount;
        while (!stack.empty()) {
            const int q = stack.back();
            stack.pop_back();
            for (int i = 0; i < 4; ++i) {
                const int r = neighbor[q][i];
                if (r >= 0 && cellComponent[r] < 0) {
                    cellComponent[r] = componentCount;
                    stack.push_back(r);
                }
            }
        }
        ++componentCount;
    }
    // Boundary loops per component: connected pieces of the boundary graph.
    const int NV = numVertices();
    std::vector<int> parent(NV);
    std::iota(parent.begin(), parent.end(), 0);
    std::function<int(int)> find = [&](int x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };
    for (int e = 0; e < numEdges(); ++e) {
        if (!boundaryEdge[e]) continue;
        const int a = find(edges[e][0]), b = find(edges[e][1]);
        if (a != b) parent[a] = b;
    }
    componentHoles.assign(componentCount, 0);
    std::vector<int> loopsOf(componentCount, 0);
    std::vector<char> rootSeen(NV, 0);
    for (int e = 0; e < numEdges(); ++e) {
        if (!boundaryEdge[e]) continue;
        const int r = find(edges[e][0]);
        if (rootSeen[r]) continue;
        rootSeen[r] = 1;
        ++loopsOf[cellComponent[edgeCell[e][0]]];
    }
    for (int c = 0; c < componentCount; ++c) componentHoles[c] = std::max(0, loopsOf[c] - 1);
}

int SquareCarrier::cornerOf(int q, int v) const {
    for (int c = 0; c < 4; ++c) if (cells[q][c] == v) return c;
    return -1;
}

void SquareCarrier::vertexEdges(int v, std::vector<int> &out) const {
    out.clear();
    const int b = ringPtr[v], n = valence[v];
    for (int k = 0; k < n; ++k) out.push_back(cellEdges[ringCell[b + k]][ringCorner[b + k]]);
    if (boundaryVertex[v] && n > 0) {
        const int q = ringCell[b + n - 1], c = ringCorner[b + n - 1];
        out.push_back(cellEdges[q][(c + 3) & 3]);
    }
}

int SquareCarrier::edgeBetween(int a, int b) const {
    for (int k = ringPtr[a]; k < ringPtr[a + 1]; ++k) {
        const int q = ringCell[k], c = ringCorner[k];
        if (cells[q][(c + 1) & 3] == b) return cellEdges[q][c];
        if (cells[q][(c + 3) & 3] == b) return cellEdges[q][(c + 3) & 3];
    }
    return -1;
}

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------
double SquareCarrier::cornerDet(int q, int c) const {
    const Point &p = vertices[cells[q][c]];
    return cross2(vertices[cells[q][(c + 1) & 3]] - p, vertices[cells[q][(c + 3) & 3]] - p);
}

double SquareCarrier::cornerAngle(int q, int c) const {
    const Point &p = vertices[cells[q][c]];
    double a = angleBetween(vertices[cells[q][(c + 1) & 3]] - p, vertices[cells[q][(c + 3) & 3]] - p);
    if (a < 0.0) a += 2.0 * M_PI;
    return a;
}

double SquareCarrier::scaledJacobian(int q, int c) const {
    const Point &p = vertices[cells[q][c]];
    const Point u = vertices[cells[q][(c + 1) & 3]] - p;
    const Point w = vertices[cells[q][(c + 3) & 3]] - p;
    const double d = normP(u) * normP(w);
    return d > 0.0 ? cross2(u, w) / d : -1.0;
}

double SquareCarrier::minScaledJacobian(int q) const {
    double m = 1.0;
    for (int c = 0; c < 4; ++c) m = std::min(m, scaledJacobian(q, c));
    return m;
}

double SquareCarrier::cellArea(int q) const {
    double A = 0.0;
    for (int c = 0; c < 4; ++c) A += cross2(vertices[cells[q][c]], vertices[cells[q][(c + 1) & 3]]);
    return 0.5 * A;
}

Point SquareCarrier::cellCentroid(int q) const {
    Point s{0.0, 0.0};
    for (int c = 0; c < 4; ++c) s = s + vertices[cells[q][c]];
    return s * 0.25;
}

double SquareCarrier::boundaryAngle(int v) const {
    if (sourceVertex[v] >= 0) return domain_->interiorAngle[sourceVertex[v]];
    return M_PI;
}

int SquareCarrier::targetValence(int v) const {
    if (!boundaryVertex[v]) return 4;
    const int t = static_cast<int>(std::lround(boundaryAngle(v) / M_PI_2));
    return std::max(1, std::min(4, t));
}

int SquareCarrier::defect(int v) const {
    return std::abs(valence[v] - targetValence(v));
}

// Sec. 4 and the no-T-junction rule of Sec. 6 together: a vertex that is a
// corner of *some* block must be a corner of *every* block at it, and a block
// can only run through a vertex (as a side-interior or interior point) where
// it has exactly two or four cells arranged straight. So:
//
//   * an interior vertex of valence other than 4 is forced;
//   * a boundary vertex of valence other than 2 is forced (a block side
//     running along dS has two cells at each of its interior vertices);
//   * a protected vertex is forced by Stage 1;
//   * a vertex on an interface is forced unless the interface runs straight
//     through it with two cells on each side.
bool SquareCarrier::forcedMacrovertex(int v) const {
    if (protectedVertex[v]) return true;
    if (boundaryVertex[v]) return valence[v] != 2 || interfaceDegree[v] > 0;
    if (valence[v] != 4) return true;
    if (interfaceDegree[v] == 0) return false;
    if (interfaceDegree[v] != 2) return true;
    // Valence 4 with two interface edges: they must be opposite in the ring.
    std::vector<int> es;
    vertexEdges(v, es);
    if (es.size() != 4) return true;
    return !(interfaceEdge[es[0]] == interfaceEdge[es[2]] &&
             interfaceEdge[es[1]] == interfaceEdge[es[3]] &&
             interfaceEdge[es[0]] != interfaceEdge[es[1]]);
}

double SquareCarrier::meanEdgeLength() const {
    if (edges.empty()) return 0.0;
    double s = 0.0;
    for (const auto &e : edges) s += normP(vertices[e[1]] - vertices[e[0]]);
    return s / edges.size();
}

// ---------------------------------------------------------------------------
// Edits
// ---------------------------------------------------------------------------
void SquareCarrier::apply(const std::vector<Edit> &edits) {
    std::vector<char> removed(cells.size(), 0);
    for (const Edit &ed : edits) for (int q : ed.removeCells) removed[q] = 1;

    for (const Edit &ed : edits) {
        const int base = numVertices();
        const int n = static_cast<int>(ed.newVertices.size());
        for (int k = 0; k < n; ++k) {
            vertices.push_back(ed.newVertices[k]);
            vertexOrigin.push_back(k < static_cast<int>(ed.newOrigin.size()) ? ed.newOrigin[k]
                                                                               : Origin::Template);
            sourceVertex.push_back(-1);
            sourceEdge.push_back(k < static_cast<int>(ed.newSourceEdge.size()) ? ed.newSourceEdge[k] : -1);
            protectedVertex.push_back(0);
            designatedVertex.push_back(0);
        }
        auto decode = [&](int id) { return id >= 0 ? id : base + (-1 - id); };
        for (const auto &q : ed.cells) {
            cells.push_back({decode(q[0]), decode(q[1]), decode(q[2]), decode(q[3])});
            cellMaterial.push_back(ed.material);
            cellTriangle.push_back(-1);
            cellOrigin.push_back(ed.cellOrigin);
            cellGroup.push_back(ed.group);
            removed.push_back(0);
        }
        for (int id : ed.designate) designatedVertex[decode(id)] = 1;
    }

    // Drop the removed cells, keeping order otherwise.
    int w = 0;
    for (int q = 0; q < numCells(); ++q) {
        if (removed[q]) continue;
        cells[w] = cells[q];
        cellMaterial[w] = cellMaterial[q];
        cellTriangle[w] = cellTriangle[q];
        cellOrigin[w] = cellOrigin[q];
        cellGroup[w] = cellGroup[q];
        ++w;
    }
    cells.resize(w);
    cellMaterial.resize(w);
    cellTriangle.resize(w);
    cellOrigin.resize(w);
    cellGroup.resize(w);

    initialSplit_ = false;
    rebuild();
}

// ---------------------------------------------------------------------------
// Validation: Sec. 3.1, Sec. 2.1, Sec. 4, Sec. 9.2
// ---------------------------------------------------------------------------
const SquareCarrier::Report &SquareCarrier::validate() {
    Report &R = report_;
    const int nonManifoldEdges = R.nonManifoldEdges;
    const int nonManifoldVertices = R.nonManifoldVertices;
    R = Report();
    R.nonManifoldEdges = nonManifoldEdges;
    R.nonManifoldVertices = nonManifoldVertices;
    const Mesh &mesh = domain_->getMesh();
    const int NV = numVertices(), NQ = numCells(), NE = numEdges();
    R.vertices = NV;
    R.cells = NQ;
    R.edges = NE;
    R.sourceTriangles = static_cast<int>(mesh.triangles.size());
    R.components = componentCount;
    for (int h : componentHoles) R.holes += h;
    for (int e = 0; e < NE; ++e) {
        if (boundaryEdge[e]) ++R.boundaryEdges;
        if (interfaceEdge[e]) ++R.interfaceEdges;
    }
    for (int q = 0; q < NQ; ++q) {
        switch (cellOrigin[q]) {
            case CellOrigin::Split: ++R.splitCells; break;
            case CellOrigin::Template: ++R.templateCells; break;
            case CellOrigin::Rewrite: ++R.rewriteCells; break;
        }
    }
    for (int v = 0; v < NV; ++v) {
        if (protectedVertex[v]) ++R.protectedVertices;
        if (designatedVertex[v]) ++R.designatedVertices;
    }

    // ---- cells: Sec. 9.1's four-corner certificate, and Sec. 3.1's bound ---
    R.initialSplit = initialSplit_;
    R.minScaledJacobian = 1.0;
    R.minCornerDetRatio = 1e300;
    R.maxCornerDetRatio = -1e300;
    double sjSum = 0.0;
    static const double kRatio[4] = {1.0 / 4.0, 1.0 / 6.0, 1.0 / 12.0, 1.0 / 6.0};
    for (int q = 0; q < NQ; ++q) {
        double minSJ = 1.0;
        bool convex = true;
        for (int c = 0; c < 4; ++c) {
            const double d = cornerDet(q, c);
            if (!(d > 0.0)) convex = false;
            minSJ = std::min(minSJ, scaledJacobian(q, c));
            if (initialSplit_ && cellTriangle[q] >= 0) {
                const double twoA = 2.0 * domain_->triangleArea[cellTriangle[q]];
                const double r = d / twoA;
                R.minCornerDetRatio = std::min(R.minCornerDetRatio, r);
                R.maxCornerDetRatio = std::max(R.maxCornerDetRatio, r);
                if (std::fabs(r - kRatio[c]) > 1e-9 * std::max(1.0, kRatio[c] * 12.0)) {
                    ++R.jacobianBoundViolations;
                }
            }
        }
        if (!convex) ++R.nonConvexCells;
        R.minScaledJacobian = std::min(R.minScaledJacobian, minSJ);
        sjSum += minSJ;
        R.area += cellArea(q);
    }
    R.meanScaledJacobian = NQ > 0 ? sjSum / NQ : 0.0;
    if (!initialSplit_) { R.minCornerDetRatio = 0.0; R.maxCornerDetRatio = 0.0; }
    const double domainArea = domain_->getReport().area;
    R.areaError = domainArea != 0.0 ? std::fabs(R.area - domainArea) / std::fabs(domainArea) : 0.0;

    // ---- transports: g_rq = g_qr^-1, and the shared edge maps onto itself --
    for (int e = 0; e < NE; ++e) {
        const int q = edgeCell[e][0], r = edgeCell[e][1];
        if (r < 0) continue;
        const int i = edgeSide[e][0], j = edgeSide[e][1];
        const SquareTransport &g = transport[q][i];
        const SquareTransport &h = transport[r][j];
        bool ok = g.compose(h).isIdentity() && h.compose(g).isIdentity();
        ok = ok && g.apply(squareCorner(j)) == squareCorner(i + 1) &&
             g.apply(squareCorner(j + 1)) == squareCorner(i);
        ok = ok && cells[r][j] == cells[q][(i + 1) & 3] && cells[r][(j + 1) & 3] == cells[q][i];
        if (!ok) ++R.transportMismatches;
    }

    // ---- Sec. 9.2: angle sums, i.e. the one-rings wind exactly once --------
    std::vector<double> angleSum(NV, 0.0);
    for (int q = 0; q < NQ; ++q) {
        for (int c = 0; c < 4; ++c) angleSum[cells[q][c]] += cornerAngle(q, c);
    }
    for (int v = 0; v < NV; ++v) {
        const double want = boundaryVertex[v] ? boundaryAngle(v) : 2.0 * M_PI;
        if (std::fabs(angleSum[v] - want) > opts_.angleTolerance) ++R.angleSumViolations;
    }

    // ---- Sec. 4: the Euler identity per component, and boundary parity ----
    {
        std::vector<int> vertexComp(NV, -1);
        for (int q = 0; q < NQ; ++q) for (int v : cells[q]) vertexComp[v] = cellComponent[q];
        std::vector<long long> lhs(componentCount, 0);
        for (int v = 0; v < NV; ++v) {
            if (vertexComp[v] < 0) continue;
            lhs[vertexComp[v]] += boundaryVertex[v] ? (2 - valence[v]) : (4 - valence[v]);
        }
        R.eulerHolds = true;
        for (int c = 0; c < componentCount; ++c) {
            const int rhs = 4 * (1 - componentHoles[c]);
            R.eulerLHS += static_cast<int>(lhs[c]);
            R.eulerRHS += rhs;
            if (lhs[c] != rhs) R.eulerHolds = false;
        }
        std::vector<int> bdy(componentCount, 0);
        for (int e = 0; e < NE; ++e) if (boundaryEdge[e]) ++bdy[cellComponent[edgeCell[e][0]]];
        R.boundaryParityEven = true;
        for (int c = 0; c < componentCount; ++c) if (bdy[c] % 2 != 0) R.boundaryParityEven = false;
    }

    // ---- Sec. 1.1: geometric preservation of dS ---------------------------
    //
    // A carrier boundary edge is legitimate iff both its ends lie on one input
    // boundary edge: an input vertex at one of its ends, or a midpoint or
    // inserted point on it. And no input boundary vertex may have vanished.
    {
        auto onSource = [&](int v, int me) {
            if (sourceVertex[v] >= 0) {
                return mesh.edges[me][0] == sourceVertex[v] || mesh.edges[me][1] == sourceVertex[v];
            }
            return sourceEdge[v] == me;
        };
        // Input boundary edges at each input vertex, to avoid the scan above
        // being quadratic on a large dS.
        std::vector<std::vector<int>> bAt(mesh.vertices.size());
        for (int me : mesh.boundaryEdges) {
            bAt[mesh.edges[me][0]].push_back(me);
            bAt[mesh.edges[me][1]].push_back(me);
        }
        std::vector<int> cand;
        for (int e = 0; e < NE; ++e) {
            if (!boundaryEdge[e]) continue;
            const int a = edges[e][0], b = edges[e][1];
            cand.clear();
            if (sourceEdge[a] >= 0) cand.push_back(sourceEdge[a]);
            if (sourceVertex[a] >= 0) for (int me : bAt[sourceVertex[a]]) cand.push_back(me);
            bool ok = false;
            for (int me : cand) {
                if (me >= 0 && mesh.edgeTriangles[me][1] < 0 && onSource(b, me)) { ok = true; break; }
            }
            if (!ok) ++R.boundaryPreservationErrors;
        }
        std::vector<char> present(mesh.vertices.size(), 0);
        for (int v = 0; v < NV; ++v) {
            if (sourceVertex[v] >= 0 && boundaryVertex[v]) present[sourceVertex[v]] = 1;
        }
        for (int mv : mesh.boundaryVertices) if (!present[mv]) ++R.missingBoundaryVertices;
    }

    // ---- irregularity ------------------------------------------------------
    for (int v = 0; v < NV; ++v) {
        const int d = defect(v);
        R.totalDefect += d;
        if (boundaryVertex[v]) { if (d != 0) ++R.irregularBoundary; }
        else if (valence[v] != 4) ++R.irregularInterior;
    }

    R.valid = R.nonManifoldEdges == 0 && R.nonManifoldVertices == 0 && R.nonConvexCells == 0 &&
              R.transportMismatches == 0 && R.angleSumViolations == 0 &&
              R.areaError <= opts_.areaTolerance && R.eulerHolds && R.boundaryParityEven &&
              R.boundaryPreservationErrors == 0 && R.missingBoundaryVertices == 0 &&
              (!initialSplit_ || R.jacobianBoundViolations == 0);

    auto msg = [&](bool bad, const std::string &m) { if (bad) R.messages.push_back("Stage 2: " + m); };
    msg(R.nonManifoldEdges > 0, std::to_string(R.nonManifoldEdges) + " non-manifold edge(s)");
    msg(R.nonManifoldVertices > 0, std::to_string(R.nonManifoldVertices) + " non-manifold vertex/vertices");
    msg(R.nonConvexCells > 0, std::to_string(R.nonConvexCells) + " cell(s) with a non-positive corner");
    msg(R.transportMismatches > 0, std::to_string(R.transportMismatches) + " transport mismatch(es)");
    msg(R.angleSumViolations > 0, std::to_string(R.angleSumViolations) + " vertex/vertices whose one-ring does not wind once");
    msg(R.areaError > opts_.areaTolerance, "the cells do not cover the domain's area");
    msg(!R.eulerHolds, "the Euler identity of Sec. 4 fails");
    msg(!R.boundaryParityEven, "a component has an odd number of boundary edges");
    msg(R.boundaryPreservationErrors > 0, std::to_string(R.boundaryPreservationErrors) + " boundary edge(s) off the input boundary");
    msg(R.missingBoundaryVertices > 0, std::to_string(R.missingBoundaryVertices) + " input boundary vertex/vertices lost");
    msg(initialSplit_ && R.jacobianBoundViolations > 0,
        std::to_string(R.jacobianBoundViolations) + " corner determinant(s) off Sec. 3.1's 2A(1/4,1/6,1/12,1/6)");
    return R;
}

bool SquareCarrier::writeOBJ(const std::string &path) const {
    std::ofstream out(path);
    if (!out) return false;
    out.precision(17);
    out << "# ATLAS square-transport carrier: " << numCells() << " cells\n";
    for (const Point &p : vertices) out << "v " << p[0] << " " << p[1] << " 0\n";
    std::map<int, std::vector<int>> byMat;
    for (int q = 0; q < numCells(); ++q) byMat[cellMaterial[q]].push_back(q);
    for (const auto &kv : byMat) {
        out << "usemtl mat" << kv.first << "\n";
        for (int q : kv.second) {
            out << "f " << cells[q][0] + 1 << " " << cells[q][1] + 1 << " " << cells[q][2] + 1 << " "
                << cells[q][3] + 1 << "\n";
        }
    }
    return static_cast<bool>(out);
}
