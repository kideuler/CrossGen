#ifndef __SHAPEDNA_LAGRANGEFEM_HXX__
#define __SHAPEDNA_LAGRANGEFEM_HXX__

#include <array>
#include <vector>

#include "mesh/Mesh.hxx"
#include "SparseMatrix.hxx"

// Lagrange finite elements of degree 1-3 on a planar triangle Mesh, and the two
// matrices of docs/shape_dna.md Sec. 4:
//
//     A_lm = int grad F_l . grad F_m  dx      (stiffness)
//     B_lm = int F_l F_m              dx      (mass)
//
// ## Nodes
//
// A degree-p element has its nodes at the lattice points (i/p, j/p) of the
// reference triangle: the three vertices, p-1 on each edge, and (p-1)(p-2)/2
// inside. The global numbering puts the mesh vertices first -- so node v *is*
// vertex v, and an eigenfunction read at the first |V| nodes is read at the
// mesh vertices -- then the edge nodes, p-1 per Mesh::edges entry, ordered from
// edges[e][0] to edges[e][1], then the interior nodes, triangle by triangle.
// Two triangles sharing an edge traverse it in opposite directions, and
// latticeNode() is where that is reconciled.
//
// The same lattice, taken at degree r, is the vertex set of the mesh with every
// edge split into r parts -- Sec. 2's global refinement -- which is how
// refineUniform() gets a conforming refinement without any bookkeeping of its
// own.
//
// ## Element integrals
//
// Sec. 4's point about flat triangles: the element integrals are computed once
// on the reference triangle and carried to each element by its affine map. With
// e1 = x1 - x0, e2 = x2 - x0 and d = det[e1 e2],
//
//     A_e = ( |e2|^2 S_xx - (e1.e2) (S_xy + S_yx) + |e1|^2 S_yy ) / |d|
//     B_e = |d| M_ref
//
// where S and M_ref are the reference integrals of products of the basis
// functions and their derivatives. Those are exact: the basis is expanded in
// monomials (by inverting the Vandermonde matrix at the nodes) and the
// monomials are integrated in closed form, int xi^a eta^b = a! b! / (a+b+2)!.
// No quadrature rule, and nothing to choose.
namespace shapedna {

// The Lagrange basis of degree p on the reference triangle (0,0), (1,0), (0,1),
// and its reference integrals. Local node k is the lattice point lattice[k];
// the local order is j = 0..p, i = 0..p-j within it. Every matrix is
// nodes x nodes, row-major.
struct ReferenceElement {
    int degree = 1;
    int nodes = 3;
    std::vector<std::array<int, 2>> lattice;
    std::vector<double> mass;        // int F_a F_b
    std::vector<double> stiffXX;     // int dF_a/dxi  dF_b/dxi
    std::vector<double> stiffYY;     // int dF_a/deta dF_b/deta
    std::vector<double> stiffXY;     // int dF_a/dxi dF_b/deta + dF_a/deta dF_b/dxi

    explicit ReferenceElement(int degree);
};

// The global nodes of degree-p Lagrange elements on a mesh (any p >= 1; the
// matrices are only assembled for p <= 3).
class LagrangeSpace {
public:
    LagrangeSpace(const Mesh &mesh, int degree);

    int degree() const { return p; }
    int nodesPerElement() const { return nloc; }
    int numNodes() const { return static_cast<int>(positions.size()); }
    int numElements() const { return static_cast<int>(mesh.triangles.size()); }
    const Mesh &getMesh() const { return mesh; }

    const std::vector<Point> &nodePositions() const { return positions; }
    // 1 when the node lies on dS: a boundary vertex, or on a boundary edge.
    const std::vector<char> &onBoundary() const { return boundary; }
    // numElements() x nodesPerElement(), in the reference element's local order.
    const std::vector<int> &elementNodes() const { return elemNodes; }

    // The global node at lattice point (i, j) of triangle t, i + j <= p.
    int latticeNode(int t, int i, int j) const;

private:
    int edgeNode(int t, int localEdge, int steps) const;

    const Mesh &mesh;
    int p = 1, nloc = 3;
    int nV = 0, nE = 0, perEdge = 0, perTri = 0;
    std::vector<int> interiorIndex;       // (p+1) x (p+1) lattice -> interior slot
    std::vector<Point> positions;
    std::vector<char> boundary;
    std::vector<int> elemNodes;
};

// Stiffness and mass on the free nodes. Under Dirichlet conditions (Sec. 5)
// every node on dS is left out of the unknowns; under Neumann they all stay.
// Nodes no triangle uses are never unknowns. The two matrices share one
// sparsity pattern, entry for entry.
struct FEMatrices {
    SparseMatrix stiffness, mass;
    std::vector<int> dofNode;                         // unknown -> node
    std::vector<int> nodeDof;                         // node -> unknown, or -1
    std::vector<std::array<double, 2>> dofPositions;  // for the dissection
    int degenerateElements = 0;                       // zero-area triangles, skipped
};

FEMatrices assemble(const LagrangeSpace &space, const ReferenceElement &ref, bool dirichlet);

// Sec. 2's refinement: every edge split into r equal parts, every triangle into
// r^2 similar ones with its orientation and material. The first |V| vertices of
// the result are the input's, in the same order.
Mesh refineUniform(const Mesh &mesh, int r);

}  // namespace shapedna

#endif  // __SHAPEDNA_LAGRANGEFEM_HXX__
