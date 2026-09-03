#ifndef __CROSSFIELD_HXX__
#define __CROSSFIELD_HXX__

#include <complex>
#include <cmath>
#include <memory>
#include <utility>

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

    // Multiply the tau = D^2/10 heuristic by this factor.
    //
    // The same knob SIPG::setTauScale is, and it exists for the same reason:
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
    double tau; // time step size
    double tauScale = 1.0; // multiplier on the D^2/10 heuristic, see setTauScale
    int maxIterations; // maximum number of iterations
};

#endif // __CROSSFIELD_HXX__
