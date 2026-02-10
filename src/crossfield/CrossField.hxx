#ifndef __CROSSFIELD_HXX__
#define __CROSSFIELD_HXX__

#include <complex>
#include <cmath>
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

    CrossField(Mesh &mesh, int maxIterations = 100)
    : mesh(mesh), maxIterations(maxIterations) {}; 

    void initialize(int method = 0);

    void step();

    void computeSingularities();

    void runMBO();

private:
    Mesh &mesh;



    // sparse complex matrices for the linear system
    Eigen::SparseMatrix<std::complex<double>> M; // Mass matrix
    Eigen::SparseMatrix<std::complex<double>> K; // Stiffness matrix
    Eigen::SparseMatrix<std::complex<double>> A; // System matrix (M + tau*K)
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> solverLU; // sparse direct solver
    Eigen::BiCGSTAB<Eigen::SparseMatrix<std::complex<double>>> solverBiCGSTAB; // sparse iterative solver (fallback)
    bool useBiCGSTAB = false; // flag to indicate which solver to use
    double tau; // time step size
    int maxIterations; // maximum number of iterations
};

#endif // __CROSSFIELD_HXX__
