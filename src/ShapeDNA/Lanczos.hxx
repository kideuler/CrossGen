#ifndef __SHAPEDNA_LANCZOS_HXX__
#define __SHAPEDNA_LANCZOS_HXX__

#include <string>
#include <vector>

#include "SparseCholesky.hxx"
#include "SparseMatrix.hxx"

// The eigensolver for docs/shape_dna.md Sec. 7: the lowest eigenpairs of the
// generalized symmetric-definite problem
//
//     A x = lambda B x,      A = stiffness (SPSD),  B = mass (SPD),
//
// by shift-invert block Lanczos with thick restarts. Written for this module;
// neither ARPACK nor Spectra nor Eigen's eigensolvers are involved.
//
// ## Shift-invert
//
// The eigenvalues wanted are the smallest, and they are the worst-separated
// part of the spectrum: lambda_k grows like 4 pi k / |Omega| (Weyl) while the
// largest is O(h^-2), so a Krylov method on B^{-1} A finds the top of the
// spectrum and reaches the bottom only after O(1/h) steps. Inverting about a
// shift sigma below lambda_1 turns the problem round:
//
//     Op = (A - sigma B)^{-1} B,     Op x = theta x,    theta = 1 / (lambda - sigma),
//
// so the smallest lambda become the largest, best-separated theta. Op is
// self-adjoint in the B-inner product <x, y>_B = x^T B y, so Lanczos applies
// with every inner product taken in B. One application of Op is a product with
// B and a solve with K = A - sigma B, which SparseCholesky has factored once.
//
// ## Blocks
//
// A symmetric domain has repeated eigenvalues -- every non-radial mode of a
// disk is double, and a square has doubles and worse -- and a single-vector
// Krylov space holds one direction of each eigenspace in exact arithmetic. The
// second copy then appears only through rounding, late, and the Shape-DNA is
// shifted by one from that index on if the iteration stops first. A block of b
// start vectors resolves multiplicities up to b directly, and gives the solve
// b right-hand sides to push through the elimination tree at once.
//
// ## Thick restart
//
// The basis is kept to basisSize vectors. When it is full, the Ritz pairs of
// the projected matrix T = V^T B Op V are computed (symmetricEigen), the best
// `keep` Ritz vectors are kept and the rest discarded, and the Lanczos process
// continues from the residual block -- the symmetric Krylov-Schur restart
// (Wu and Simon, SIAM J. Matrix Anal. Appl. 22, 2000; Stewart 2001). T is kept
// dense and filled from the orthogonalisation coefficients, so after a restart
// it is diagonal plus the arrow coupling the kept vectors to the residual
// block, and nothing special is needed to carry on.
//
// Orthogonality is maintained by full reorthogonalisation, classical Gram-
// Schmidt run twice (CGS2) against the whole basis, which is two block inner
// products and two block updates per step: level-3 work, split over threads.
// A direction that the new block has already lost (the Krylov space is
// invariant, or has run out of room) is replaced by a random vector orthogonal
// to everything, so a breakdown never stops the iteration.
namespace shapedna {

struct LanczosOptions {
    int nev = 10;                // eigenpairs wanted: the smallest lambda above sigma
    int blockSize = 4;
    int basisSize = 0;           // 0: 2 nev + 2 blockSize, rounded to the block
    double tolerance = 1e-10;    // on || Op x - theta x ||_B / theta
    int maxRestarts = 300;
    unsigned long long seed = 1;
    bool vectors = false;        // return the eigenvectors too
};

struct LanczosResult {
    bool converged = false;
    std::vector<double> values;      // lambda, ascending, nev of them
    std::vector<double> residuals;   // || Op x - theta x ||_B / theta for each
    std::vector<double> vectors;     // n x nev, column-major, B-orthonormal (if asked)
    int basisSize = 0, blockSize = 0;
    int restarts = 0;
    int blockSteps = 0;
    int solves = 0;                  // right-hand sides pushed through the factor
    int deflations = 0;              // lost directions replaced by random ones
    std::string message;
};

// K must be the factorisation of A - sigma B. Returns result.converged.
bool shiftInvertLanczos(const SparseMatrix &B, const SparseCholesky &K, double sigma,
                        const LanczosOptions &options, LanczosResult &result);

}  // namespace shapedna

#endif  // __SHAPEDNA_LANCZOS_HXX__
