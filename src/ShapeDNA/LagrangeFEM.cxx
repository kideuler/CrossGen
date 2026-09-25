// LagrangeFEM.cxx -- see LagrangeFEM.hxx.
#include "LagrangeFEM.hxx"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "Parallel.hxx"

namespace shapedna {

namespace {

// int_{reference triangle} xi^a eta^b = a! b! / (a + b + 2)!
double monomialIntegral(int a, int b) {
    double num = 1.0;
    for (int k = 2; k <= a; ++k) num *= k;
    for (int k = 2; k <= b; ++k) num *= k;
    double den = 1.0;
    for (int k = 2; k <= a + b + 2; ++k) den *= k;
    return num / den;
}

// Solves the n x n system Vd C = I (row-major Vd) by Gaussian elimination with
// partial pivoting; C comes back row-major. n is at most 10 here.
std::vector<double> inverse(int n, std::vector<double> Vd) {
    std::vector<double> C(static_cast<std::size_t>(n) * n, 0.0);
    for (int i = 0; i < n; ++i) C[i * n + i] = 1.0;
    for (int k = 0; k < n; ++k) {
        int piv = k;
        for (int i = k + 1; i < n; ++i)
            if (std::fabs(Vd[i * n + k]) > std::fabs(Vd[piv * n + k])) piv = i;
        if (Vd[piv * n + k] == 0.0) throw std::runtime_error("ReferenceElement: singular Vandermonde matrix");
        if (piv != k)
            for (int j = 0; j < n; ++j) {
                std::swap(Vd[k * n + j], Vd[piv * n + j]);
                std::swap(C[k * n + j], C[piv * n + j]);
            }
        const double inv = 1.0 / Vd[k * n + k];
        for (int j = 0; j < n; ++j) {
            Vd[k * n + j] *= inv;
            C[k * n + j] *= inv;
        }
        for (int i = 0; i < n; ++i) {
            if (i == k) continue;
            const double f = Vd[i * n + k];
            if (f == 0.0) continue;
            for (int j = 0; j < n; ++j) {
                Vd[i * n + j] -= f * Vd[k * n + j];
                C[i * n + j] -= f * C[k * n + j];
            }
        }
    }
    return C;
}

}  // namespace

// ---------------------------------------------------------------------------
// ReferenceElement
// ---------------------------------------------------------------------------

ReferenceElement::ReferenceElement(int deg) : degree(deg) {
    if (degree < 1) throw std::invalid_argument("ReferenceElement: degree must be at least 1");
    const int p = degree;
    for (int j = 0; j <= p; ++j)
        for (int i = 0; i <= p - j; ++i) lattice.push_back({i, j});
    nodes = static_cast<int>(lattice.size());

    // Monomials xi^a eta^b, a + b <= p: as many as there are nodes.
    std::vector<std::array<int, 2>> mono;
    for (int d = 0; d <= p; ++d)
        for (int a = d; a >= 0; --a) mono.push_back({a, d - a});

    // Vd(k, m) = monomial m at node k. Its inverse holds the basis: F_k is
    // sum_m C(m, k) monomial_m, which is 1 at node k and 0 at the others.
    const int n = nodes;
    std::vector<double> Vd(static_cast<std::size_t>(n) * n);
    for (int k = 0; k < n; ++k) {
        const double xi = static_cast<double>(lattice[k][0]) / p;
        const double eta = static_cast<double>(lattice[k][1]) / p;
        for (int m = 0; m < n; ++m) Vd[k * n + m] = std::pow(xi, mono[m][0]) * std::pow(eta, mono[m][1]);
    }
    const std::vector<double> C = inverse(n, Vd);
    auto coef = [&](int m, int k) { return C[m * n + k]; };

    // Each entry computed once, for a <= b, and mirrored: summed separately
    // the two halves would differ in the last bit, and an element matrix that
    // is not exactly symmetric makes the assembled one not exactly symmetric
    // either -- which the Cholesky, reading one triangle, and the products,
    // reading both, would then see differently.
    mass.assign(static_cast<std::size_t>(n) * n, 0.0);
    stiffXX = stiffYY = stiffXY = mass;
    for (int a = 0; a < n; ++a)
        for (int b = a; b < n; ++b) {
            double mab = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0;
            for (int m = 0; m < n; ++m) {
                const double ca = coef(m, a);
                if (ca == 0.0) continue;
                const int am = mono[m][0], bm = mono[m][1];
                for (int q = 0; q < n; ++q) {
                    const double cb = coef(q, b);
                    if (cb == 0.0) continue;
                    const int aq = mono[q][0], bq = mono[q][1];
                    const double cc = ca * cb;
                    mab += cc * monomialIntegral(am + aq, bm + bq);
                    // d/dxi xi^a eta^b = a xi^(a-1) eta^b, and likewise for eta.
                    if (am > 0 && aq > 0) sxx += cc * am * aq * monomialIntegral(am + aq - 2, bm + bq);
                    if (bm > 0 && bq > 0) syy += cc * bm * bq * monomialIntegral(am + aq, bm + bq - 2);
                    if (am > 0 && bq > 0) sxy += cc * am * bq * monomialIntegral(am - 1 + aq, bm + bq - 1);
                    if (bm > 0 && aq > 0) sxy += cc * bm * aq * monomialIntegral(am + aq - 1, bm - 1 + bq);
                }
            }
            mass[a * n + b] = mass[b * n + a] = mab;
            stiffXX[a * n + b] = stiffXX[b * n + a] = sxx;
            stiffYY[a * n + b] = stiffYY[b * n + a] = syy;
            stiffXY[a * n + b] = stiffXY[b * n + a] = sxy;
        }
}

// ---------------------------------------------------------------------------
// LagrangeSpace
// ---------------------------------------------------------------------------

LagrangeSpace::LagrangeSpace(const Mesh &m, int degree) : mesh(m), p(degree) {
    if (p < 1) throw std::invalid_argument("LagrangeSpace: degree must be at least 1");
    nloc = (p + 1) * (p + 2) / 2;
    nV = static_cast<int>(mesh.vertices.size());
    nE = static_cast<int>(mesh.edges.size());
    const int nT = static_cast<int>(mesh.triangles.size());
    perEdge = p - 1;
    perTri = (p - 1) * (p - 2) / 2;

    interiorIndex.assign((p + 1) * (p + 1), -1);
    int q = 0;
    for (int j = 1; j <= p; ++j)
        for (int i = 1; i + j <= p - 1; ++i) interiorIndex[i * (p + 1) + j] = q++;

    const std::size_t total = static_cast<std::size_t>(nV) + static_cast<std::size_t>(nE) * perEdge +
                              static_cast<std::size_t>(nT) * perTri;
    positions.resize(total);
    boundary.assign(total, 0);

    for (int v = 0; v < nV; ++v) {
        positions[v] = mesh.vertices[v];
        boundary[v] = mesh.isBoundaryVertex[v] ? 1 : 0;
    }
    SDNA_OMP(parallel for schedule(static) if(nE > 4096))
    for (int e = 0; e < nE; ++e) {
        const Point &a = mesh.vertices[mesh.edges[e][0]];
        const Point &b = mesh.vertices[mesh.edges[e][1]];
        for (int s = 0; s < perEdge; ++s) {
            const double f = static_cast<double>(s + 1) / p;
            const std::size_t k = static_cast<std::size_t>(nV) + static_cast<std::size_t>(e) * perEdge + s;
            positions[k] = a + (b - a) * f;
            boundary[k] = mesh.isBoundaryEdge[e] ? 1 : 0;
        }
    }

    // The element's nodes in the reference element's local order (j outer, i
    // inner), placing the interior ones as they are met: each interior node
    // belongs to one triangle, so no two threads write the same position.
    elemNodes.resize(static_cast<std::size_t>(nT) * nloc);
    SDNA_OMP(parallel for schedule(static) if(nT > 2048))
    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = mesh.triangles[t];
        const Point &x0 = mesh.vertices[tri[0]];
        const Point e1 = mesh.vertices[tri[1]] - x0;
        const Point e2 = mesh.vertices[tri[2]] - x0;
        int k = 0;
        for (int j = 0; j <= p; ++j)
            for (int i = 0; i <= p - j; ++i, ++k) {
                const int node = latticeNode(t, i, j);
                elemNodes[static_cast<std::size_t>(t) * nloc + k] = node;
                if (node >= nV + nE * perEdge)
                    positions[node] = x0 + e1 * (static_cast<double>(i) / p) + e2 * (static_cast<double>(j) / p);
            }
    }
}

int LagrangeSpace::edgeNode(int t, int localEdge, int steps) const {
    const int e = mesh.triangleEdges[t][localEdge];
    const int start = mesh.triangles[t][localEdge];
    const int sub = (start == mesh.edges[e][0]) ? steps - 1 : p - 1 - steps;
    return nV + e * perEdge + sub;
}

int LagrangeSpace::latticeNode(int t, int i, int j) const {
    const Triangle &tri = mesh.triangles[t];
    const int k = p - i - j;
    if (i == 0 && j == 0) return tri[0];
    if (i == p) return tri[1];
    if (j == p) return tri[2];
    if (j == 0) return edgeNode(t, 0, i);        // v0 -> v1, i steps from v0
    if (k == 0) return edgeNode(t, 1, j);        // v1 -> v2, j steps from v1
    if (i == 0) return edgeNode(t, 2, p - j);    // v2 -> v0, p - j steps from v2
    return nV + nE * perEdge + t * perTri + interiorIndex[i * (p + 1) + j];
}

// ---------------------------------------------------------------------------
// Assembly
// ---------------------------------------------------------------------------
//
// Row by row rather than element by element. An element-by-element scatter has
// every element write into the rows of all its nodes, so two threads on
// neighbouring elements collide; this instead gives each row to one thread,
// which walks the elements around that row's node and adds their entries in
// element order. No atomics, no colouring, and each entry summed in a fixed
// order.
FEMatrices assemble(const LagrangeSpace &space, const ReferenceElement &ref, bool dirichlet) {
    FEMatrices out;
    const Mesh &mesh = space.getMesh();
    const int nT = space.numElements();
    const int nloc = space.nodesPerElement();
    const int nN = space.numNodes();
    const std::vector<int> &en = space.elementNodes();
    if (ref.nodes != nloc) throw std::invalid_argument("assemble: reference element of the wrong degree");

    // Which nodes are unknowns.
    std::vector<char> used(nN, 0);
    for (int x : en) used[x] = 1;
    out.nodeDof.assign(nN, -1);
    for (int v = 0; v < nN; ++v)
        if (used[v] && !(dirichlet && space.onBoundary()[v])) {
            out.nodeDof[v] = static_cast<int>(out.dofNode.size());
            out.dofNode.push_back(v);
        }
    const int nd = static_cast<int>(out.dofNode.size());
    out.dofPositions.resize(nd);
    for (int d = 0; d < nd; ++d) out.dofPositions[d] = space.nodePositions()[out.dofNode[d]];

    // Node -> (element, local slot), in element order.
    std::vector<int> incPtr(nN + 1, 0), inc(static_cast<std::size_t>(nT) * nloc);
    for (int x : en) ++incPtr[x + 1];
    for (int v = 0; v < nN; ++v) incPtr[v + 1] += incPtr[v];
    {
        std::vector<int> fill(incPtr.begin(), incPtr.end() - 1);
        for (std::size_t s = 0; s < en.size(); ++s) inc[fill[en[s]]++] = static_cast<int>(s);
    }

    // Element matrices, all at once.
    const std::size_t nl2 = static_cast<std::size_t>(nloc) * nloc;
    std::vector<double> Ke(static_cast<std::size_t>(nT) * nl2), Me(static_cast<std::size_t>(nT) * nl2);
    int degenerate = 0;
    SDNA_OMP(parallel for schedule(dynamic, 256) if(nT > 1024))
    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = mesh.triangles[t];
        const Point &x0 = mesh.vertices[tri[0]];
        const Point e1 = mesh.vertices[tri[1]] - x0;
        const Point e2 = mesh.vertices[tri[2]] - x0;
        const double d = std::fabs(cross2(e1, e2));
        double *K = Ke.data() + t * nl2;
        double *M = Me.data() + t * nl2;
        const double g11 = dotP(e2, e2), g12 = dotP(e1, e2), g22 = dotP(e1, e1);
        if (!(d > 1e-14 * (g11 + g22))) {
            std::fill(K, K + nl2, 0.0);
            std::fill(M, M + nl2, 0.0);
            SDNA_OMP(atomic)
            ++degenerate;
            continue;
        }
        const double inv = 1.0 / d;
        for (std::size_t k = 0; k < nl2; ++k) {
            K[k] = (g11 * ref.stiffXX[k] - g12 * ref.stiffXY[k] + g22 * ref.stiffYY[k]) * inv;
            M[k] = ref.mass[k] * d;
        }
    }
    out.degenerateElements = degenerate;

    // The pattern: row d holds every unknown sharing an element with node d.
    std::vector<std::vector<int>> rowCols(nd);
    SDNA_OMP(parallel for schedule(dynamic, 256) if(nd > 4096))
    for (int d = 0; d < nd; ++d) {
        const int v = out.dofNode[d];
        std::vector<int> &cols = rowCols[d];
        for (int k = incPtr[v]; k < incPtr[v + 1]; ++k) {
            const int t = inc[k] / nloc;
            for (int b = 0; b < nloc; ++b) {
                const int c = out.nodeDof[en[static_cast<std::size_t>(t) * nloc + b]];
                if (c >= 0) cols.push_back(c);
            }
        }
        std::sort(cols.begin(), cols.end());
        cols.erase(std::unique(cols.begin(), cols.end()), cols.end());
    }

    SparseMatrix &A = out.stiffness;
    A.n = nd;
    A.rowPtr.assign(nd + 1, 0);
    for (int d = 0; d < nd; ++d) A.rowPtr[d + 1] = A.rowPtr[d] + static_cast<int>(rowCols[d].size());
    A.col.resize(A.rowPtr[nd]);
    A.val.assign(A.rowPtr[nd], 0.0);
    SparseMatrix &B = out.mass;
    B.n = nd;
    B.rowPtr = A.rowPtr;
    B.val.assign(A.rowPtr[nd], 0.0);

    SDNA_OMP(parallel for schedule(dynamic, 256) if(nd > 4096))
    for (int d = 0; d < nd; ++d) {
        const int v = out.dofNode[d];
        const int r0 = A.rowPtr[d];
        std::copy(rowCols[d].begin(), rowCols[d].end(), A.col.begin() + r0);
        const int *cb = A.col.data() + r0;
        const int *ce = A.col.data() + A.rowPtr[d + 1];
        for (int k = incPtr[v]; k < incPtr[v + 1]; ++k) {
            const int t = inc[k] / nloc;
            const int a = inc[k] % nloc;
            const double *K = Ke.data() + t * nl2 + static_cast<std::size_t>(a) * nloc;
            const double *M = Me.data() + t * nl2 + static_cast<std::size_t>(a) * nloc;
            for (int b = 0; b < nloc; ++b) {
                const int c = out.nodeDof[en[static_cast<std::size_t>(t) * nloc + b]];
                if (c < 0) continue;
                const int pos = static_cast<int>(std::lower_bound(cb, ce, c) - A.col.data());
                A.val[pos] += K[b];
                B.val[pos] += M[b];
            }
        }
        std::vector<int>().swap(rowCols[d]);
    }
    B.col = A.col;
    return out;
}

// ---------------------------------------------------------------------------
// Refinement
// ---------------------------------------------------------------------------

Mesh refineUniform(const Mesh &mesh, int r) {
    if (r <= 1) return Mesh(mesh.vertices, mesh.triangles, mesh.triangleMatId);
    const LagrangeSpace lattice(mesh, r);
    const int nT = lattice.numElements();
    std::vector<Triangle> tris;
    std::vector<int> mats;
    tris.reserve(static_cast<std::size_t>(nT) * r * r);
    mats.reserve(tris.capacity());
    for (int t = 0; t < nT; ++t) {
        const int mat = mesh.triangleMatId.empty() ? 1 : mesh.triangleMatId[t];
        for (int j = 0; j < r; ++j)
            for (int i = 0; i + j < r; ++i) {
                // The upward triangle at (i, j), and the downward one beside it.
                tris.push_back({lattice.latticeNode(t, i, j), lattice.latticeNode(t, i + 1, j),
                                lattice.latticeNode(t, i, j + 1)});
                mats.push_back(mat);
                if (i + j + 2 <= r) {
                    tris.push_back({lattice.latticeNode(t, i + 1, j), lattice.latticeNode(t, i + 1, j + 1),
                                    lattice.latticeNode(t, i, j + 1)});
                    mats.push_back(mat);
                }
            }
    }
    return Mesh(lattice.nodePositions(), tris, mats);
}

}  // namespace shapedna
