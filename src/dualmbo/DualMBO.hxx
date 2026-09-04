#ifndef __DUALMBO_HXX__
#define __DUALMBO_HXX__

#include <complex>
#include <cmath>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>
#include <limits>

// eigen includes
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "mesh/Mesh.hxx"

class DualMBO {
public:
    Eigen::VectorXcd u_k;      // current solution vector (per triangle)
    Eigen::VectorXcd u_k_prev; // previous solution vector (per triangle)
    double error = std::numeric_limits<double>::max(); // current error for convergence checking
    std::vector<std::pair<int, double>> singularVertices; // (vertex index, cross-field index) pairs
    std::shared_ptr<Mesh> mesh;

    DualMBO(std::shared_ptr<Mesh> mesh, int maxIterations = 100, double gamma = 10.0)
        : mesh(mesh), maxIterations(maxIterations), gamma(gamma) {}

    void initialize();

    // Interior edges the field is to be aligned to, on top of dS.
    //
    // The MBO field is boundary-aligned and nothing else: its only Dirichlet
    // data is the tangent of dS, so on a multi-material domain it reads the
    // material interfaces as ordinary interior edges and runs straight through
    // them. That is wrong for a layout, because an interface is a curve the
    // output has to keep exactly as dS is -- and the cost of ignoring it is not
    // a slightly worse field, it is the wrong *singularities*: each material
    // region carries its own index count, and the cones that count asks for
    // simply are not there in a field that never saw the interface.
    //
    // Passing the interface edges here makes them one-sided Dirichlet edges on
    // both sides, exactly as a boundary edge is on its one side, so the field is
    // tangent to the interface from either material and its holonomy is then
    // read per region. Call before initialize(); an empty set is the old
    // behaviour and is what a single-material mesh gives.
    // The two triangles on such an edge are pinned outright, the way a boundary
    // triangle is. Leaving them in the diffusion instead, so that the alignment
    // only competes with smoothness at the edge penalty weight, was measured
    // and is worse: on a strongly curved interface it drops the +1/-1 pairs
    // that a smoothest aligned field carries (which is what is wanted) but it
    // also collapses the pair a *material* boundary genuinely needs onto the
    // interface itself -- geom001's two cones land 0.14 apart astride the arc
    // instead of 0.72 apart in the middle of their own regions. The pairs that
    // are not wanted come off in Stage 1 instead; see
    // ConeSingularities::cancelDipoles.
    void setAlignedInteriorEdges(const std::vector<int> &edges);

    // Kill the rotational degree of freedom of a disk.
    //
    // A disk is rotationally symmetric, so the boundary-aligned cross field on
    // it is too: rotate the domain and the problem maps to itself, and the four
    // +1/4 cones that Gauss-Bonnet asks for slide around with it. Nothing in
    // the boundary data picks their angular position, so the MBO lands them
    // wherever arithmetic noise puts them, and two runs of the same disk -- or
    // the ten copies of one disk in data/geometry/multimat/bubbles.geo -- need
    // not agree.
    //
    // Pinning the triangle at the center of the disk to u = exp(4i*theta) = 1
    // removes it. Near the center the field is, to leading order, the
    // extension of the boundary data z^4 perturbed off zero: u(z) = z^4 -
    // eps^4 c, whose four simple zeros are the four cones, at radius eps and
    // at angles arg(c)/4 + k*pi/2. The value at the center is u(0) = -eps^4 c,
    // so fixing it fixes arg(c): u(0) = 1 gives c = -1 and puts the cones at
    // pi/4 + k*pi/2, measured from the center. That is the canonical position
    // -- the cones on the diagonals, the cross at the center axis-aligned --
    // and it is reproducible across runs and across disks.
    //
    // On by default. The pin only ever applies to a component of the mesh
    // whose boundary Mesh::computeMaterialCircles accepted as a circle, so a
    // domain with no disk in it is unaffected, and a center triangle that is
    // already Dirichlet from dS or an interface is left alone rather than
    // fought over. Call before initialize().
    void setPinDiskCenters(bool on) { pinDiskCenters = on; }

    // The triangles pinned by the above, one per disk, in materialComponents
    // order. Empty until initialize() has run.
    const std::vector<int>& getDiskCenterTriangles() const { return diskCenterTriangles; }

    // Multiply the tau = D^2/10 heuristic of Step 5 by this factor.
    //
    // 1 is the heuristic and is what every caller in the pipeline uses. It is
    // settable because the heuristic is a heuristic: tau is the diffusion time
    // one MBO step takes, and how far it may be moved before the singularity
    // count starts to depend on it is a measurement, not a derivation. That
    // measurement is paper_tests/E1_Verification, and this is the knob it
    // turns. Call before initialize().
    void setTauScale(double s) { tauScale = s; }
    double getTau() const { return tau; }
    double getGamma() const { return gamma; }

    // Impose the boundary and interface alignment weakly instead of pinning.
    //
    // The assembly of initialize() builds the Nitsche/penalty form for the
    // Dirichlet data whether or not it is pinned: K carries kappa_e on the
    // diagonal of every triangle with an aligned edge and b carries
    // kappa_e * g_e, so A = M + tau*K with that RHS *is* the weak form. Step 6
    // then eliminates those rows, which turns the weak statement into a hard
    // one -- the triangle takes its prescribed cross exactly and the alignment
    // stops competing with smoothness.
    //
    // Hard is the default and is what the pipeline runs. Off leaves the
    // elimination out and solves the weak form as assembled, which is the
    // ablation the paper reports: the two differ only in how much the field is
    // allowed to trade alignment against the edge penalty at the boundary, and
    // the honest way to say which one the results use is to be able to run
    // both. Call before initialize().
    void setHardBoundaryConditions(bool on) { hardBoundary = on; }

    // Weight a triangle's alignment pin by how consistent its edges' data is.
    //
    // The Dirichlet data of a triangle carrying two aligned edges is the
    // kappa-weighted average of exp(4i theta) over them. That average is exact
    // where the two tangents are a quarter turn apart -- both edges then ask for
    // the same cross -- and *cancels* where they are an eighth turn apart, since
    // exp(4i*0) and exp(4i*pi/4) are antipodal. At such a corner the assembled
    // pin has no direction and full strength, which is the worst of both: see
    // the long comment in initialize(). On is the default and scales the pin by
    // the coherence r = |sum kappa g| / sum kappa, which is 1 in the ordinary
    // one-edge case and so changes nothing there. Off restores the old
    // assembly, which is the ablation. Call before initialize().
    void setCornerCoherenceFix(bool on) { cornerCoherenceFix = on; }

    // Coherence below which the averaged direction is not worth pinning to and
    // the triangle is left to diffuse. 45 degrees of disagreement gives r = 0.
    static constexpr double kCornerCoherenceMin = 0.2;

    // How the two elements on an interior edge combine into one penalty weight.
    //
    // At p=0 the volume, consistency and symmetry terms of the interior-penalty form all
    // carry grad u_h and so vanish identically on a piecewise constant. What is
    // left is a(u,v) = sum_e kappa_e [u][v], which means kappa_e is not a free
    // stabilisation parameter here: it *is* the discrete Laplacian, and it has
    // to be chosen to be one.
    //
    // The textbook interior-penalty weight is not one. gamma |e| / min(h_i, h_j)
    // comes from a coercivity bound, not from consistency, and it combines the
    // two half-fluxes in *parallel* where a two-point flux combines them in
    // series. The test is the one any second-order operator has to pass: applied
    // to a linear field sampled at the circumcenters it must give zero. Measured
    // over the 35-model corpus it gives a relative residual of 0.37 -- an O(1)
    // spurious force on a field with nothing in it left to smooth.
    //
    // Orthogonal is the two-point (finite-volume) weight |e| / d_e, with d_e the
    // distance between the two circumcenters. Both circumcenters lie on the
    // perpendicular bisector of e, so c_j - c_i is orthogonal to e for *any*
    // triangle pair and
    //
    //     (phi_j - phi_i) / d_e = grad phi . n_e     exactly for linear phi,
    //
    // so the fluxes telescope and the residual above is 6e-15 instead of 0.37.
    // This is the operator the classical finite-volume scheme uses on the
    // Delaunay/Voronoi pair; the factor 2/3 makes it agree with MinHeight on an
    // equilateral pair, so gamma (and with it the tau ladder) keeps its meaning.
    //
    // HarmonicHeight is the series combination of the same two heights the
    // penalty already uses -- the SWIP choice. It is here because it is the
    // obvious repair and it is worth being able to show that it is not enough:
    // these meshes are near-uniform (area ratios across an edge reach only 2.4),
    // so it moves the residual by 2% and nothing else.
    enum class PenaltyWeight {
        MinHeight,      // gamma |e| / min(h_i, h_j)          -- the textbook one
        HarmonicHeight, // gamma |e| / ((h_i + h_j) / 2)      -- SWIP
        Orthogonal      // (2/3) gamma |e| / d_e              -- two-point/FV
    };
    void setPenaltyWeight(PenaltyWeight w) { penaltyWeight = w; }
    PenaltyWeight getPenaltyWeight() const { return penaltyWeight; }

    // Backward-Euler substeps per MBO step.
    //
    // Threshold dynamics is "run the heat flow for a time tau, then project",
    // and the operator that does the first half is the semigroup exp(-tau L).
    // One backward-Euler solve is not it: it damps a mode of eigenvalue lambda
    // by 1/(1 + tau lambda) where the semigroup damps it by exp(-tau lambda),
    // and the two part company exactly where it matters. At the continuation's
    // floor the stiffest mode has tau*lambda ~ 400, so the semigroup leaves
    // e^-400 = 0 of it and one Euler step leaves 1/401 -- a quarter of a percent
    // of the sharpest content the mesh can carry, surviving every step, right
    // where the projection is most able to turn it into structure.
    //
    // n substeps of tau/n replace 1/(1+tau lambda) by (1 + tau lambda/n)^-n,
    // which is the standard rational approximation of the semigroup and is
    // monotonically closer to it for every mode. The cost is n solves per step
    // against one factorisation, so it is n times the solve time and no more
    // assembly.
    void setDiffusionSubsteps(int n) { diffusionSubsteps = n < 1 ? 1 : n; }

    // The diffusion rate of a typical element: the median over triangles of
    // K_ii / M_ii, in units of 1/length^2. Valid after initialize().
    //
    // This is what a tau-continuation's floor has to be measured against. The
    // floor is the point where one MBO step stops resolving anything the mesh
    // can carry -- the diffusion length ell = h sqrt(tau * lambda) reaching a
    // few edges -- and lambda is a property of the *assembled operator*, not of
    // gamma. Deriving it from gamma alone (ell = sqrt(8 gamma tau), which is
    // K_ii ~ 3 gamma and M_ii ~ h^2/2) is exact only for the MinHeight weight on
    // a uniform mesh, so a run that changes the weight and keeps that formula is
    // silently annealing to a different depth and is not a controlled
    // comparison. Reading the rate off K and M instead makes tauFloorEdges mean
    // the same thing for every weight.
    //
    // The median rather than the maximum: a single sliver would otherwise set
    // the floor for the whole mesh, and the question the floor answers is what
    // the mesh typically resolves. On a uniform mesh with the MinHeight weight
    // it evaluates to 8 gamma / h^2, so the ladder is unchanged there.
    double medianDiffusionRate() const { return medianRate; }

    // Floor on d_e, as a fraction of (h_i + h_j)/2.
    //
    // d_e vanishes when the two triangles are cocircular -- the dual edge has
    // zero length and the transmissibility |e|/d_e is unbounded. It is a
    // measure-zero configuration that a structured or recombined patch sits
    // exactly on, so it has to be handled rather than assumed away: over the
    // corpus's 292069 interior edges, 19 fall below this floor and exactly one
    // is non-Delaunay (d_e < 0). Clamping there under-couples an edge whose two
    // cell centres coincide anyway, which is the harmless direction.
    static constexpr double kMinDualDistance = 0.05;

    // The Dirichlet data as initialize() computed it: triangle -> the unit
    // spin-4 value its incident boundary/interface edges ask for. Non-empty
    // whether or not the pin is hard, because it is the *data* and not the way
    // of imposing it. Read by the alignment-error metrics.
    const std::unordered_map<int, std::complex<double>>& getBoundaryData() const {
        return boundaryTriangleBC;
    }

    // The assembled operator, read-only, for the consistency experiment
    // (paper_tests/E1_Verification, part (e)): K is the edge-penalty stiffness
    // with the Dirichlet diagonal terms in it, M the mass with the pinned rows
    // replaced by identity rows. Neither touches the row of a triangle whose
    // three edges are all interior, which is the only kind that test reads.
    // Valid after initialize().
    const Eigen::SparseMatrix<std::complex<double>>& stiffnessMatrix() const { return K; }
    const Eigen::SparseMatrix<std::complex<double>>& massMatrix() const { return M; }

    void step();

    void computeSingularities();

    void runMBO();

    // Access the underlying mesh
    const Mesh& getMesh() const { return *mesh; }
    std::shared_ptr<Mesh> getMeshPtr() const { return mesh; }

private:
    // sparse complex matrices for the linear system
    Eigen::SparseMatrix<std::complex<double>> M; // Mass matrix (diagonal, real)
    Eigen::SparseMatrix<std::complex<double>> K; // Stiffness matrix (DualMBO edge penalty)
    Eigen::SparseMatrix<std::complex<double>> A; // System matrix (M + tau*K)
    Eigen::VectorXcd b;                          // Boundary right-hand side vector
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> solverLU;      // sparse direct solver
    Eigen::BiCGSTAB<Eigen::SparseMatrix<std::complex<double>>> solverBiCGSTAB; // iterative fallback
    bool useBiCGSTAB = false; // flag to indicate which solver to use

    double tau;          // time step size
    double tauScale = 1.0; // multiplier on the D^2/10 heuristic, see setTauScale
    double gamma;        // edge penalty parameter
    int maxIterations;   // maximum number of iterations
    bool hardBoundary = true; // see setHardBoundaryConditions
    bool cornerCoherenceFix = true; // see setCornerCoherenceFix
    PenaltyWeight penaltyWeight = PenaltyWeight::MinHeight; // see setPenaltyWeight
    int diffusionSubsteps = 1; // see setDiffusionSubsteps
    double medianRate = 0.0;   // see medianDiffusionRate

    // Interior edges promoted to aligned (Dirichlet) edges, per edge of the
    // mesh. Empty when there are none, which is the single-material case.
    std::vector<char> edgeAligned;

    bool pinDiskCenters = true;          // see setPinDiskCenters
    std::vector<int> diskCenterTriangles; // triangles pinned to 1 by it

    // Hard Dirichlet BC per boundary triangle: triangle index -> prescribed exp(4i*theta)
    // Computed from the dominant boundary edge tangent in initialize(); applied after each step.
    // Triangles on an aligned interior edge are pinned the same way.
    std::unordered_map<int, std::complex<double>> boundaryTriangleBC;

    bool isAlignedEdge(int e) const {
        return e >= 0 && e < static_cast<int>(edgeAligned.size()) && edgeAligned[e];
    }
};

#endif // __DUALMBO_HXX__