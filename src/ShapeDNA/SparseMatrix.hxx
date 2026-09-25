#ifndef __SHAPEDNA_SPARSEMATRIX_HXX__
#define __SHAPEDNA_SPARSEMATRIX_HXX__

#include <vector>

// The sparse matrix the Shape-DNA module works in: square, compressed sparse
// row, with *both* triangles of a symmetric matrix stored and the columns of
// every row sorted.
//
// Storing both triangles doubles the memory of the two finite-element matrices
// and buys a product with no scatter in it: row i of y = A x reads x and writes
// y[i] only, so the rows can be split over threads with no atomics and no
// colouring, and each y[i] is summed in the same order whatever the split.
// Both FE matrices come out of the same assembly with the same pattern (see
// LagrangeFEM), so K = A - sigma B is formed entry by entry on `val`.
namespace shapedna {

struct SparseMatrix {
    int n = 0;
    std::vector<int> rowPtr;   // n + 1 offsets into col/val
    std::vector<int> col;      // column of each entry, ascending within a row
    std::vector<double> val;

    long long nnz() const { return rowPtr.empty() ? 0 : static_cast<long long>(rowPtr[n]); }

    // y = A x. Rows in parallel.
    void multiply(const double *x, double *y) const;

    // Y = A X for `nvec` column-major vectors of leading dimension `ld` (>= n):
    // one pass over the matrix for the whole block, which is what the block
    // Lanczos iteration asks for every step.
    void multiply(const double *X, double *Y, int nvec, int ld) const;

    // A(i, j), or 0 when (i, j) is not in the pattern. Binary search.
    double entry(int i, int j) const;
};

}  // namespace shapedna

#endif  // __SHAPEDNA_SPARSEMATRIX_HXX__
