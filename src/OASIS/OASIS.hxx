#ifndef __OASIS_HXX__
#define __OASIS_HXX__

#include <memory>
#include <vector>

// eigen includes
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "crossfield/CrossField.hxx"
#include "mesh/Mesh.hxx"

// Assembly of the quasi-eigenfunction (QE) KKT system of
//   Ling, Huang, Juttler, Sun, Bao and Wang, "Spectral Quadrangulation with
//   Feature Curve Alignment and Element Size Control", ACM TOG 34(1), 2014.
//   (refs/Huang2014.pdf)
//
// The QE is the minimizer of the least-squares Helmholtz residual (Eq. 12)
//
//     E_Lr(f) = || L_r f - lambda f ||^2
//
// subject to the boundary-alignment conditions (Eq. 6), which make every
// boundary curve an extended minimal integral curve of f:
//
//     df/dn = 0     and     d^2f/dn^2 = xi = 1     on the boundary.
//
// Discretized through the local quadratic fit of Eq. 9, both conditions are
// linear in the vertex values f, giving the constraint (Eq. 11)
//
//     B f = [Y; Z] f = [0; 1] = C.
//
// A quadratic objective with linear constraints is solved directly through its
// KKT system (Eq. 13)
//
//     [ H_Lr  B^T ] [ f  ]   [ 0 ]
//     [ B      0  ] [ nu ] = [ C ],
//
// where nu are the Lagrange multipliers and (Eq. 14)
//
//     H_Lr = L_r^T L_r - lambda (L_r + L_r^T) + lambda^2 I,
//
// i.e. (L_r - lambda I)^T (L_r - lambda I), the normal equations of Eq. 12.
//
// This class only assembles that system. Extraction of the Morse-Smale complex
// and the vibration enhancement pass of Sec. 3.4 are not implemented here.
class OASIS {
public:
    // The quasi-eigenfunction, one value per vertex. Filled by solve().
    Eigen::VectorXd f;

    // Lagrange multipliers nu of Eq. 13, one per constraint row. Nonzero
    // multipliers are expected: they are the price of the alignment
    // conditions, which the Helmholtz residual alone would not satisfy.
    Eigen::VectorXd nu;

    std::shared_ptr<Mesh> mesh;

    // Guiding cross field. Unused by the Eq. 13 assembly, which is the whole
    // point of the paper's boundary conditions: feature alignment needs no
    // direction field. It is kept for the optional orientation term E_Orient of
    // Sec. 5.1, where the combined energy E_Lr(f) + gamma * E_Orient(f) steers
    // mesh orientation away from boundaries and features.
    std::shared_ptr<CrossField> crossField;

    // lambda < 0 is the Helmholtz parameter of Eq. 1. It need not be an
    // eigenvalue; it acts as the global element-size knob (Sec. 5).
    explicit OASIS(std::shared_ptr<CrossField> crossField, double lambda = -1000.0);
    explicit OASIS(std::shared_ptr<Mesh> mesh, double lambda = -1000.0);

    // Element-size field r > 0 of Eq. 8, one entry per vertex. The resulting
    // quad edge length is roughly proportional to 1/sqrt(r). Defaults to the
    // uniform field r == 1.
    void setDensity(const Eigen::VectorXd &density);

    void setLambda(double value) { lambda = value; assembled = false; }
    double getLambda() const { return lambda; }

    // Builds L_r (Eq. 8), B and C (Eq. 11), H_Lr (Eq. 14) and finally the KKT
    // system of Eq. 13 together with its right-hand side [0; C].
    void assemble();

    // Factorizes and solves Eq. 13, filling f and nu. Calls assemble() first if
    // that has not happened yet. Throws if the factorization or solve fails.
    void solve();

    // ||L_r f - lambda f||, the Helmholtz residual of Eq. 12 left over after
    // the alignment constraints. Nonzero by construction: that is precisely
    // what makes f a *quasi*-eigenfunction. It grows as lambda moves away from
    // a true eigenvalue of the constrained problem.
    double residual() const;

    const Eigen::SparseMatrix<double>& getKKT() const { return KKT; }
    const Eigen::VectorXd& getKKTRhs() const { return kktRhs; }
    const Eigen::SparseMatrix<double>& getLaplacian() const { return Lr; }
    const Eigen::SparseMatrix<double>& getHessian() const { return HLr; }
    const Eigen::SparseMatrix<double>& getConstraintMatrix() const { return B; }
    const Eigen::VectorXd& getConstraintRhs() const { return C; }

    // Vertices carrying the Eq. 6 conditions. Each contributes two rows to B:
    // the Y (first derivative) row at index i and the Z (second derivative) row
    // at index constrainedVertices.size() + i. That indexing holds only if no
    // rows were dropped as redundant — compare numConstraints() against
    // 2 * getConstrainedVertices().size() to tell.
    const std::vector<int>& getConstrainedVertices() const { return constrainedVertices; }
    int numConstraints() const { return static_cast<int>(B.rows()); }

    // Access the underlying mesh
    const Mesh& getMesh() const { return *mesh; }
    std::shared_ptr<Mesh> getMeshPtr() const { return mesh; }

private:
    // Indices of the quadratic fit coefficients of Eq. 9, in the order they are
    // produced by fitLocalQuadratic():
    //   f(u,v) = 1/2 (a_uu u^2 + 2 a_uv u v + a_vv v^2) + a_u u + a_v v + a_c
    enum QuadCoeff { AUU = 0, AUV = 1, AVV = 2, AU = 3, AV = 4, AC = 5, NUM_COEFF = 6 };

    void buildDensityLaplacian(); // Eq. 8
    void buildConstraints();      // Eqs. 9-11
    void buildKKT();              // Eqs. 13-14

    // Reduce B (and C) to a full-row-rank system. Eq. 13 is singular otherwise,
    // and SparseLU factorizes the singular matrix without complaining.
    void dropDependentConstraints();

    // Outward unit normal per vertex, length-weighted over the incident
    // boundary edges. Zero for interior vertices.
    std::vector<Point> computeVertexNormals() const;

    // Vertex ids within `rings` rings of v, v itself first.
    std::vector<int> gatherNeighborhood(int v, int rings) const;

    // Least-squares quadratic fit of Eq. 9 around vertex v. On success fills
    // `neighbors` with the stencil vertex ids and `coeffs` with the
    // NUM_COEFF x |neighbors| matrix P satisfying
    //   (a_uu, a_uv, a_vv, a_u, a_v, a_c)^T = P * f[neighbors],
    // so that each a_k is the linear form a_k^T f of the paper. P depends only
    // on the discretization, never on f. `ringsUsed` reports how wide the
    // stencil had to grow to make the fit well posed.
    bool fitLocalQuadratic(int v, std::vector<int> &neighbors, Eigen::MatrixXd &coeffs,
                           int &ringsUsed) const;

    double lambda;     // Helmholtz parameter, negative
    Eigen::VectorXd r; // density field of Eq. 8, one entry per vertex

    std::vector<int> constrainedVertices;   // boundary (and later feature) vertices
    std::vector<Point> constrainedNormals;  // matching outward unit normals

    Eigen::SparseMatrix<double> Lr;  // density-modulated cotangent Laplacian, Eq. 8
    Eigen::SparseMatrix<double> HLr; // Eq. 14
    Eigen::SparseMatrix<double> B;   // Eq. 11, 2*|constrained| x |V|
    Eigen::VectorXd C;               // Eq. 11 right-hand side, [0; 1]
    Eigen::SparseMatrix<double> KKT; // Eq. 13
    Eigen::VectorXd kktRhs;          // [0; C]

    // The KKT matrix is symmetric *indefinite*, so a Cholesky-type solver
    // (SimplicialLDLT) is not applicable; the paper uses UMFPACK, and SparseLU
    // is the direct equivalent available here.
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solverLU;

    bool assembled = false;
};

#endif // __OASIS_HXX__
