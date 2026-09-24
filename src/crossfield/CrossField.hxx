#ifndef __CROSSFIELD_HXX__
#define __CROSSFIELD_HXX__

#include <complex>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

// eigen includes
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "mesh/Mesh.hxx"

class CrossField {
public:
    Eigen::VectorXcd u_k; // current solution vector
    Eigen::VectorXcd u_k_prev; // previous solution vector
    double error = std::numeric_limits<double>::max(); // current error for convergence checking
    std::vector<std::pair<int, double>> singularTriangles; // (triangle index, cross-field index) pairs
    std::shared_ptr<Mesh> mesh;

    CrossField(std::shared_ptr<Mesh> mesh, int maxIterations = 100)
    : mesh(mesh), maxIterations(maxIterations) {}; 

    // `seed` 0 takes the seed from the clock, as before; anything else makes
    // the random start of method 1 reproducible, which is what lets two runs of
    // the pipeline be compared against each other rather than against the
    // clock.
    void initialize(int method = 0, unsigned seed = 0);

    // Interior edges the field is to be aligned to, on top of dS: the material
    // interfaces of a multi-material domain. DualMBO::setAlignedInteriorEdges
    // is the same call on the p=0 field, and the reason is the same: a field
    // whose only Dirichlet data is the tangent of dS runs straight through
    // every interface, and each material region then lacks the singularities
    // its own index count asks for.
    //
    // Per vertex rather than per triangle, which is the one thing that
    // changes. DualMBO pins the two triangles either side of an interface edge;
    // here every vertex on one is pinned, to the tangent of the interface, and
    // since every stiffness entry coupling the two materials runs through such
    // a vertex that decouples the regions just as DualMBO's pinned triangles
    // do. The value at a vertex is the length-weighted mean of exp(4i phi)
    // over its incident interface edges -- exact on a smooth run and at a
    // right-angled kink, T or cross, where every edge asks for the same cross.
    // Where they do not (an oblique or 120-degree junction, a 135-degree
    // kink: the ill-posed nodes of data/geometry/multimat) the mean shrinks,
    // and below DualMBO's kCornerCoherenceMin the vertex is left free rather
    // than pinned to a direction no incident curve has; the field then turns
    // there, and the tracing takes that node over (SeparatrixTrace).
    //
    // Boundary vertices keep dS's data. Call before initialize(); an empty set
    // is the old behaviour, which the single-material baseline keeps.
    void setAlignedInteriorEdges(const std::vector<int> &edges);
    // After initialize(): interface vertices pinned, and left free because
    // their incident interface directions disagreed.
    int alignedInterfaceVertices() const { return alignedInterface; }
    int freeInterfaceVertices() const { return freeInterface; }

    // Kill the rotational degree of freedom of a disk:
    // DualMBO::setPinDiskCenters, per vertex. A disk's four +1/4 cones are
    // placed by nothing in the boundary data, so an unpinned field lands them
    // wherever the random start and the arithmetic put them, and the ten
    // disks of data/geometry/multimat/bubbles.geo come out ten different ways.
    // Pinning the vertex nearest the centre of each component
    // Mesh::computeMaterialCircles accepted as a circle to u = 1 puts them on
    // the diagonals, the cross at the centre axis-aligned -- DualMBO's
    // canonical position. Off by default, so the baseline is unchanged; call
    // before initialize().
    void setPinDiskCenters(bool on) { pinDiskCenters = on; }
    int pinnedDiskCenters() const { return pinnedDisks; }

    // Multiply the tau = D^2/10 heuristic by this factor.
    //
    // The same knob DualMBO::setTauScale is, and it exists for the same reason:
    // D^2/10 makes one MBO step diffuse across several domain diameters, so the
    // iteration reaches its fixed point almost at once and the threshold
    // dynamics never runs. It is settable so that the *baseline* can be given
    // the same tau-continuation our method uses and the comparison can say
    // whether the gain is the continuation or the discretisation. 1 is the
    // published heuristic and remains the default. Call before initialize().
    void setTauScale(double s) { tauScale = s; }
    double getTau() const { return tau; }

    void step();

    void computeSingularities();

    void runMBO();

    // Access the underlying mesh
    const Mesh& getMesh() const { return *mesh; }
    std::shared_ptr<Mesh> getMeshPtr() const { return mesh; }

private:

    // sparse complex matrices for the linear system
    Eigen::SparseMatrix<std::complex<double>> M; // Mass matrix
    Eigen::SparseMatrix<std::complex<double>> K; // Stiffness matrix
    Eigen::SparseMatrix<std::complex<double>> A; // System matrix (M + tau*K)
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> solverLU; // sparse direct solver
    Eigen::BiCGSTAB<Eigen::SparseMatrix<std::complex<double>>> solverBiCGSTAB; // sparse iterative solver (fallback)
    bool useBiCGSTAB = false; // flag to indicate which solver to use
    std::vector<char> edgeAligned; // per mesh edge: see setAlignedInteriorEdges
    int alignedInterface = 0, freeInterface = 0;
    bool pinDiskCenters = false;
    int pinnedDisks = 0;
    double tau; // time step size
    double tauScale = 1.0; // multiplier on the D^2/10 heuristic, see setTauScale
    int maxIterations; // maximum number of iterations
};

#endif // __CROSSFIELD_HXX__
