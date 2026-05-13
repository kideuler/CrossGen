#ifndef __SIPG_HXX__
#define __SIPG_HXX__

#include <complex>
#include <cmath>
#include <memory>
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
    double gamma;        // SIPG penalty parameter
    int maxIterations;   // maximum number of iterations
};

#endif // __SIPG_HXX__