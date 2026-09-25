#ifndef __SHAPEDNA_DENSEEIGEN_HXX__
#define __SHAPEDNA_DENSEEIGEN_HXX__

#include <vector>

// Dense symmetric eigensolvers, written for the Shape-DNA module rather than
// taken from Eigen.
//
// The Lanczos iteration reduces the finite-element problem to a small symmetric
// matrix T (a few hundred rows at most: the Krylov basis) and needs its full
// eigen-decomposition after every expansion -- the Rayleigh-Ritz step. That is
// symmetricEigen: Householder reduction to tridiagonal form followed by the
// implicit QL iteration with Wilkinson shifts, the EISPACK tred2/tql2 pair
// (Bowdler, Martin, Reinsch and Wilkinson, Numer. Math. 11, 1968), O(m^3) with
// eigenvectors. m is the basis size, never the mesh size, so this is not where
// the time goes and it runs on one thread.
//
// generalizedSymmetricEigen is the whole problem A x = lambda B x done densely,
// by B = L L^T and the standard problem for L^{-1} A L^{-T}. ShapeDNA falls
// back on it when the mesh is too small for a Krylov basis to make sense, and
// the self-test uses it as a second, independent route to the same spectrum.
//
// All matrices are column-major. Eigenvalues come back ascending, and
// eigenvector i is column i.
namespace shapedna {

// Eigenvalues (and, when `vectors` is not null, orthonormal eigenvectors) of the
// symmetric n x n matrix `a`. Only the lower triangle of `a` is read. Returns
// false only if the QL iteration fails to converge, which for a finite
// symmetric matrix means it was given a NaN.
bool symmetricEigen(int n, const double *a, std::vector<double> &values,
                    std::vector<double> *vectors = nullptr);

// Cholesky factor of the symmetric n x n matrix `a`, in place: on return the
// lower triangle holds L with A = L L^T and the strict upper triangle is zero.
// False (and `a` left part-way) when a pivot is not positive.
bool denseCholesky(int n, double *a);

// A x = lambda B x with A symmetric and B symmetric positive definite, both
// n x n and read from their lower triangles. Eigenvectors are B-orthonormal.
// False when B is not positive definite or the QL iteration fails.
bool generalizedSymmetricEigen(int n, const double *A, const double *B,
                               std::vector<double> &values,
                               std::vector<double> *vectors = nullptr);

}  // namespace shapedna

#endif  // __SHAPEDNA_DENSEEIGEN_HXX__
