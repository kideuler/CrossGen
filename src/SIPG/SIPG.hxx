#ifndef __SIPG_HXX__
#define __SIPG_HXX__

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

class SIPG {
public:
    Eigen::VectorXcd u_k;      // current solution vector (per triangle)
    Eigen::VectorXcd u_k_prev; // previous solution vector (per triangle)
    double error = std::numeric_limits<double>::max(); // current error for convergence checking
    std::vector<std::pair<int, double>> singularVertices; // (vertex index, cross-field index) pairs
    std::shared_ptr<Mesh> mesh;

    SIPG(std::shared_ptr<Mesh> mesh, int maxIterations = 100, double gamma = 10.0)
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
    // only competes with smoothness at the SIPG penalty weight, was measured
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
    // allowed to trade alignment against the SIPG penalty at the boundary, and
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

    // The Dirichlet data as initialize() computed it: triangle -> the unit
    // spin-4 value its incident boundary/interface edges ask for. Non-empty
    // whether or not the pin is hard, because it is the *data* and not the way
    // of imposing it. Read by the alignment-error metrics.
    const std::unordered_map<int, std::complex<double>>& getBoundaryData() const {
        return boundaryTriangleBC;
    }

    void step();

    void computeSingularities();

    void runMBO();

    // Access the underlying mesh
    const Mesh& getMesh() const { return *mesh; }
    std::shared_ptr<Mesh> getMeshPtr() const { return mesh; }

private:
    // sparse complex matrices for the linear system
    Eigen::SparseMatrix<std::complex<double>> M; // Mass matrix (diagonal, real)
    Eigen::SparseMatrix<std::complex<double>> K; // Stiffness matrix (SIPG edge penalty)
    Eigen::SparseMatrix<std::complex<double>> A; // System matrix (M + tau*K)
    Eigen::VectorXcd b;                          // Boundary right-hand side vector
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> solverLU;      // sparse direct solver
    Eigen::BiCGSTAB<Eigen::SparseMatrix<std::complex<double>>> solverBiCGSTAB; // iterative fallback
    bool useBiCGSTAB = false; // flag to indicate which solver to use

    double tau;          // time step size
    double tauScale = 1.0; // multiplier on the D^2/10 heuristic, see setTauScale
    double gamma;        // SIPG penalty parameter
    int maxIterations;   // maximum number of iterations
    bool hardBoundary = true; // see setHardBoundaryConditions
    bool cornerCoherenceFix = true; // see setCornerCoherenceFix

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

#endif // __SIPG_HXX__