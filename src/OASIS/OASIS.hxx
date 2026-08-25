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
// Sec. 5.1 adds one optional term to that objective, taken over unchanged from
//   Huang, Zhang, Ma, Liu, Kobbelt and Bao, "Spectral Quadrangulation with
//   Orientation and Alignment Control", ACM TOG 27(5), 2008, Sec. 3.4:
// a guiding direction field steers the mesh orientation where the boundary
// conditions do not, by minimizing E_Lr(f) + gamma * E_Orient(f) subject to the
// same Eq. 11 constraints. See setOrientationWeight().
//
// This class only assembles that system. Extraction of the Morse-Smale complex
// is not implemented here.
class OASIS {
public:
    // The quasi-eigenfunction, one value per vertex. Filled by solve().
    Eigen::VectorXd f;

    // Lagrange multipliers nu of Eq. 13, one per constraint row. Nonzero
    // multipliers are expected: they are the price of the alignment
    // conditions, which the Helmholtz residual alone would not satisfy.
    Eigen::VectorXd nu;

    std::shared_ptr<Mesh> mesh;

    // Guiding cross field, read per vertex from CrossField::u_k (the 4-symmetry
    // representation vector e^{4 i theta}, so the guiding direction is
    // arg(u_k)/4 modulo pi/2). Nothing in the Eq. 13 assembly needs it, which is
    // the whole point of the paper's boundary conditions: feature alignment
    // needs no direction field. It drives only the optional orientation term
    // E_Orient of Sec. 5.1; see setOrientationWeight().
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

    // Relative weights of Eq. 19. The paper's omega = 0.1 and xi = 100 do NOT
    // carry over, because the three terms are balanced against each other by
    // operator magnitude here rather than in the paper's own normalization.
    // See the defaults below for the measured values and why they were chosen.
    void setVibrationWeight(double omega) { vibrationWeight = omega; }
    void setBoundaryPenaltyWeight(double xi) { boundaryPenaltyWeight = xi; }

    // gamma of the combined energy E_Lr(f) + gamma * E_Orient(f) of Sec. 5.1,
    // with E_Orient the orientation energy of Huang et al. 2008, Eq. 10:
    //
    //     E_Orient(f) = sum_i D_i <Q^i_uv, f>^2 = ||Q_orient f||^2,
    //
    // where Q^i_uv is the mixed second derivative of the Eq. 9 quadratic fit at
    // vertex i, expressed in the frame whose u-axis follows the guiding cross
    // field, and D_i is the vertex area. It vanishes exactly when the Hessian of
    // f is diagonal in that frame, i.e. when the principal directions of f --
    // and with them the arcs of the Morse-Smale complex, which run along them
    // (2008, Sec. 3.4) -- follow the cross field.
    //
    // The cross field enters only through the rotation of that frame, so its
    // 4-symmetry is respected for free: rotating the guiding direction by pi/2
    // flips the sign of Q^i_uv and leaves the energy untouched.
    //
    // As with the vibration weights, the paper's gamma = 100 does not carry
    // over, because the term is normalized here against the magnitude of H_Lr
    // rather than left in the paper's own scaling; see the default below.
    // gamma <= 0, or a null cross field, disables the term entirely and
    // restores the pure 2014 system.
    void setOrientationWeight(double gamma) { orientationWeight = gamma; assembled = false; }
    double getOrientationWeight() const { return orientationWeight; }

    // Per-vertex multiplier on the orientation term, one entry per vertex, all
    // >= 0. The 2014 paper applies its guiding field only "in regions away from
    // boundary curves or features" -- its car hood uses it on the back part
    // alone -- because near a boundary the Eq. 11 constraints already fix the
    // orientation and a conflicting direction field would only pull against
    // them. Defaults to 1 everywhere.
    void setOrientationMask(const Eigen::VectorXd &mask);

    // A ready-made mask for the usual case: 1 at vertices further than
    // `clearance` from the boundary, 0 within it. Straight-line distance, not
    // geodesic, so a thin neck is measured through the gap -- which errs
    // towards guiding less, the safe direction.
    Eigen::VectorXd boundaryClearanceMask(double clearance) const;

    // Mean angle in degrees between the principal directions of the Hessian of
    // f and the guiding cross field, area weighted over the vertices carrying
    // an orientation row. Both frames are 4-symmetric, so this lies in [0, 45]:
    // 0 is perfect alignment, 45 is as misaligned as a cross can get. This is
    // the quantity E_Orient drives down, so comparing it across gamma is the
    // direct measure of whether the term did anything. Returns -1 when there is
    // no cross field or no solution yet.
    //
    // On a solved field it does reach 0 -- a guided solve on data/meshes at a
    // large gamma reports 0.04 degrees. Only when a *continuous* field is
    // sampled onto the mesh and pushed through the Eq. 9 fit does the fit noise
    // put a floor under it, around 4 degrees on these meshes, because the mean
    // is over absolute angles. So compare solves against each other, not a
    // solve against an analytic reference.
    double orientationError() const;

    // Builds L_r (Eq. 8), B and C (Eq. 11), H_Lr (Eq. 14) and finally the KKT
    // system of Eq. 13 together with its right-hand side [0; C].
    void assemble();

    // Factorizes and solves Eq. 13, filling f and nu. Calls assemble() first if
    // that has not happened yet. Throws if the factorization or solve fails.
    void solve();

    // Vibration enhancement, Sec. 3.4. A QE can vibrate strongly in one
    // direction and barely at all in the orthogonal one — the paper's example
    // is a disk-like region with rotational symmetry, where the field varies
    // radially and hardly at all around. Its critical points are then hard to
    // detect or unevenly spaced, and the Morse-Smale complex built on them is
    // poor. This pass pushes the two local vibration amplitudes together by
    // minimizing Eq. 19 with Gauss-Newton, starting from the Eq. 13 solution.
    //
    // The boundary conditions move from hard constraints to a penalty here,
    // deliberately: the paper notes that keeping them as Lagrange multipliers
    // "prevents the scalar field near the boundary from adjusting for desired
    // vibration", which is exactly what this pass needs to do.
    //
    // Requires solve() to have run. Leaves f at the best iterate found.
    void enhanceVibration(int iterations = 10);

    // Mean of the Eq. 18 amplitude-difference energy over all vertices, in
    // [0, 2]. Zero means the two local vibration amplitudes match everywhere;
    // this is the quantity enhanceVibration() drives down, so comparing it
    // before and after is the direct measure of whether the pass did anything.
    double vibrationEnergy() const;

    // ||L_r f - lambda f||, the Helmholtz residual of Eq. 12 left over after
    // the alignment constraints. Nonzero by construction: that is precisely
    // what makes f a *quasi*-eigenfunction. It grows as lambda moves away from
    // a true eigenvalue of the constrained problem.
    double residual() const;

    const Eigen::SparseMatrix<double>& getKKT() const { return KKT; }
    const Eigen::VectorXd& getKKTRhs() const { return kktRhs; }
    const Eigen::SparseMatrix<double>& getLaplacian() const { return Lr; }
    // Eq. 14, plus the weighted orientation block when that term is enabled --
    // it is the (1,1) block of Eq. 13 either way.
    const Eigen::SparseMatrix<double>& getHessian() const { return HLr; }
    const Eigen::SparseMatrix<double>& getOrientationMatrix() const { return Qorient; }
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
    void buildOrientation();      // Eq. 10 of Huang et al. 2008
    void buildKKT();              // Eqs. 13-14

    // Vertex areas D_i = 1/3 * sum of incident triangle areas, i.e. Eq. 2 of
    // Huang et al. 2008 and the same quantity the row scaling of Eq. 8 divides
    // by. Zero for a vertex with no incident area.
    std::vector<double> computeVertexAreas() const;

    // (cos 2t, sin 2t) for the guiding cross direction t at vertex v, which is
    // all the frame rotation the orientation term needs. False when there is no
    // usable direction there (no cross field, or a vanishing representation
    // vector at a singularity).
    bool crossFrame(int v, double &cos2t, double &sin2t) const;

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
    // stencil had to grow to make the fit well posed; `startRings` is where it
    // begins (Sec. 3.3 wants two rings on the boundary, one in the interior).
    bool fitLocalQuadratic(int v, std::vector<int> &neighbors, Eigen::MatrixXd &coeffs,
                           int &ringsUsed, int startRings) const;

    // The Eq. 9 fit at one vertex, stored as the five linear forms the
    // vibration energy needs. W^u and W^v of Eq. 17 are never formed: they are
    // sums of outer products, so f^T W^u f expands to a few squared dot
    // products of these rows with f, which is far cheaper than one sparse
    // matrix per vertex.
    struct VertexFit {
        std::vector<int> stencil;  // global vertex ids
        Eigen::VectorXd au, av;    // first-derivative coefficient rows
        Eigen::VectorXd auu, auv, avv;  // second-derivative coefficient rows
        bool valid = false;
    };
    mutable std::vector<VertexFit> fits;

    // Fit every vertex, caching the coefficient rows for the vibration pass.
    void precomputeFits() const;

    // The two invariants of Eq. 16 that Eq. 18 is built from, evaluated at
    // vertex i: `sum` = A^2_u + A^2_v and `diff` = A^2_u - A^2_v, with their
    // gradients in the stencil values when requested.
    //
    // Sec. 3.4 measures the amplitudes in the *Hessian principal frame* ("using
    // these directions as the local coordinate system"), not the frame the
    // Eq. 9 fit happens to be expressed in. `sum` is a trace and so is the same
    // either way, but `diff` is not: computed naively it reports a rotated but
    // perfectly isotropic field as strongly anisotropic. Rather than
    // eigendecompose the Hessian at every vertex on every iteration, `diff` is
    // evaluated in closed form through the deviatoric identity
    // e1 e1^T - e2 e2^T = H_dev / R, which stays differentiable.
    void vibrationAnisotropy(int i, const Eigen::VectorXd &values,
                             double &sum, double &diff,
                             Eigen::VectorXd *dSum = nullptr,
                             Eigen::VectorXd *dDiff = nullptr) const;

    double lambda;     // Helmholtz parameter, negative
    Eigen::VectorXd r; // density field of Eq. 8, one entry per vertex

    // omega and xi of Eq. 19. The paper's values are 0.1 and 100; neither
    // transfers, because the terms are balanced against operator magnitude here
    // rather than in the paper's own normalization. Both defaults were picked
    // by sweeping against two measured quantities -- the mean E_a the pass is
    // supposed to reduce, and the alignment error it must not destroy:
    //
    //   omega   1 leaves a rotationally symmetric disk stuck at E_a 1.78 (its
    //           unenhanced value); 10 takes it to 0.29.
    //   xi      1e5 lets the boundary drift to 2% of the field amplitude; 1e7
    //           holds it at 5e-4 while costing almost nothing in E_a.
    // omega and xi of Eq. 19, tuned here rather than taken from the paper.
    //
    // These were picked against a measurement that matters more than E_a: on a
    // disk, how much the field varies *around* a circle of constant radius
    // compared with its mean on that circle. A ring-shaped field scores ~0; a
    // properly spotted one scores ~1. Sweeping omega at xi = 1e7 gives
    //
    //   omega     0.1     10      30      100     300
    //   E_a       1.78    0.29    0.16    0.079   0.15
    //   ang/ring  0.05    0.82    0.96    1.08    1.21
    //
    // The paper's omega = 0.1 is a no-op under this normalization: the field
    // stays a set of rings. Note also that E_a alone would have stopped at
    // omega = 10, which still leaves visible ring structure — E_a can be driven
    // partway down without the field ever developing angular variation, so it
    // is a necessary but not sufficient check.
    //
    // xi = 1e7 keeps the alignment conditions to under 1e-3 of the field
    // amplitude; 1e5 lets them drift to 2%, and 1e2 to over half.
    double vibrationWeight = 100.0;
    double boundaryPenaltyWeight = 1e7;

    // gamma of Sec. 5.1. The paper's 100 is stated against its own scaling of
    // E_Orient; here the term is normalized by operator magnitude first (see
    // orientationScale), so this is a *relative* weight against the Helmholtz
    // term. Measured on the 5x5 square of TestOasis carrying a constant cross
    // field 30 degrees off its axes, guided everywhere more than two quads from
    // the boundary -- mean misalignment (orientationError, in degrees) against
    // the Helmholtz residual that keeps the cells square:
    //
    //   gamma      0      0.01   0.1    1      10     100
    //   deg        28.4   28.3   27.1   15.9   9.7    5.6
    //   residual   8.8    10.2   13.4   15.8   17.0   17.5
    //
    // That is the hard case: a constant guiding direction disagrees with the
    // boundary conditions everywhere, so the two can only meet in the middle,
    // and the alignment stalls at 5 degrees no matter what gamma buys. What it
    // does show is where the cost lands -- the residual saturates by gamma = 1,
    // in the band where the guided interior meets the boundary conditions, and
    // further alignment past that is nearly free. (The table keeps the guiding
    // field two quads clear of the boundary, as the paper does. Guiding the
    // whole square instead, the residual never saturates: 21.4 / 33.7 / 69.9 at
    // gamma = 1 / 10 / 100. Mask the boundary region, or expect distorted cells
    // against it -- see setOrientationMask and boundaryClearanceMask.)
    //
    // The normal case is a guiding field that already agrees with the boundary,
    // e.g. one solved by MBO on the same mesh, and there the term does what it
    // says. On data/meshes/singlemat/geom001.obj, guided by its own MBO cross field:
    //
    //   gamma      0      0.1    1      10     100
    //   deg       12.6    5.7    1.4    0.35   0.04
    //   residual   9.0    10.3   10.8   11.0   11.0
    //
    // Alignment keeps improving and the residual stops moving after gamma = 1,
    // so 10 is the default: it lands at a third of a degree here, and it is
    // where the cost has already saturated in the adversarial case too.
    //
    // The mesh-only constructor zeroes this: with no guiding field there is
    // nothing to orient towards.
    double orientationWeight = 10.0;

    // Per-vertex multiplier on E_Orient. Empty means 1 everywhere.
    Eigen::VectorXd orientationMask;

    std::vector<int> constrainedVertices;   // boundary (and later feature) vertices
    std::vector<Point> constrainedNormals;  // matching outward unit normals

    Eigen::SparseMatrix<double> Lr;  // density-modulated cotangent Laplacian, Eq. 8
    Eigen::SparseMatrix<double> HLr; // Eq. 14
    Eigen::SparseMatrix<double> B;   // Eq. 11, 2*|constrained| x |V|
    Eigen::VectorXd C;               // Eq. 11 right-hand side, [0; 1]
    Eigen::SparseMatrix<double> KKT; // Eq. 13
    Eigen::VectorXd kktRhs;          // [0; C]

    // Eq. 10 of Huang et al. 2008, one row per oriented vertex, and the
    // Q^T Q block it contributes to H_Lr. Empty when the term is off.
    Eigen::SparseMatrix<double> Qorient;
    Eigen::SparseMatrix<double> QtQorient;
    std::vector<int> orientedVertices;

    // orientationWeight after normalizing Q^T Q against H_Lr, i.e. the factor
    // actually multiplying the orientation block. Zero when the term is off.
    // The vibration pass reuses it so that both phases minimize the same
    // combined energy.
    double orientationScale = 0.0;

    // The KKT matrix is symmetric *indefinite*, so a Cholesky-type solver
    // (SimplicialLDLT) is not applicable; the paper uses UMFPACK, and SparseLU
    // is the direct equivalent available here.
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solverLU;

    bool assembled = false;
};

#endif // __OASIS_HXX__
