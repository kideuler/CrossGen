#ifndef __SHAPEDNA_SPARSECHOLESKY_HXX__
#define __SHAPEDNA_SPARSECHOLESKY_HXX__

#include <array>
#include <atomic>
#include <vector>

#include "SparseMatrix.hxx"

// A sparse Cholesky factorisation, A = P^T L L^T P, for the shift-invert step of
// the Shape-DNA eigensolver: every Lanczos step solves with K = A - sigma B, and
// this factors K once so each of those solves is two triangular sweeps.
//
// ## Why not Cuthill-McKee
//
// docs/shape_dna.md Sec. 6 renumbers with Cuthill-McKee, which makes the matrix
// banded and lets a banded factorisation stay inside the band. On a planar mesh
// of N nodes the band is ~ sqrt(N) wide, so the factor holds N^1.5 entries and
// costs N^2 flops, and a band factorisation is a sequence of dependent column
// updates that does not split over threads.
//
// Nested dissection does better on both counts. The nodes are split in two by a
// separator, each half is split again, and so on; the separators are numbered
// last. Eliminating one half then never touches the other, so the factor holds
// O(N log N) entries, costs O(N^1.5) flops (George 1973), and the two halves are
// independent work -- all the way down, which is what OpenMP gets to use.
//
// Measured on singlemat/geom001 with the cubic elements the paper uses (82k
// unknowns, 2026-09-25): reverse Cuthill-McKee gives a bandwidth of 920 and an
// envelope of 1.9e7 entries (149 MB) costing ~13 GFlop to factor; this
// dissection gives 1.1e7 entries (86 MB) and 1.8 GFlop, factored in 0.05 s on
// eight threads.
//
// ## The dissection
//
// Geometric, on the node positions the finite-element space already has: the
// node set is cut at the median of its longer bounding-box axis, and the nodes
// on one side with a neighbour on the other become the separator (whichever
// side gives the smaller one). No graph partitioner is needed, and on a
// triangulated planar domain a straight cut is close to the best separator
// there is. Sets of at most Options::leafSize nodes are not cut further.
//
// ## The factorisation
//
// Multifrontal (Duff and Reid 1983; Liu, SIAM Review 34, 1992) on the tree the
// dissection builds: each tree node is a supernode -- a leaf set or a separator,
// its variables consecutive in the new numbering -- whose columns of L are
// eliminated together in one dense frontal matrix. The front of a supernode is
// assembled from its own entries of A and its children's update matrices
// (extend-add), partially factored, and hands the Schur complement on its
// remaining rows to its parent.
//
// Parallelism is on two levels, and neither changes the answer:
//
//   - sibling subtrees are independent, so each is an OpenMP task;
//   - near the root, where there are few fronts but they are large (880 rows
//     for geom001's top separator), the dense partial factorisation of a
//     single front is split by columns with taskloop.
//
// Each entry of every front is updated in a fixed order (children in list
// order, then panels in column order), so the factor is bit-identical for any
// thread count. The solve runs over the same tree: forward substitution bottom
// up, carrying each subtree's contribution to its ancestors in its own buffer
// rather than writing into shared rows, then back substitution top down.
namespace shapedna {

class SparseCholesky {
public:
    struct Options {
        // Node sets at most this large are leaves of the dissection and are
        // eliminated as one dense block. Smaller leaves mean less fill in the
        // leaves and more, smaller fronts.
        int leafSize = 64;
        // A front at least this large factors its own columns in parallel.
        int parallelFront = 256;
    };

    struct Stats {
        int n = 0;
        int supernodes = 0;
        int treeDepth = 0;
        int maxFront = 0;        // rows of the largest frontal matrix
        long long nnzA = 0;      // lower triangle of A, diagonal included
        long long nnzL = 0;      // entries stored for L (the supernodes' dense blocks)
        double flops = 0.0;      // of the numeric factorisation
        double orderSeconds = 0.0, symbolicSeconds = 0.0, numericSeconds = 0.0;
    };

    SparseCholesky() = default;
    explicit SparseCholesky(const Options &opts) : options(opts) {}

    // Orders, analyses and factors A. A must be symmetric positive definite
    // with both triangles stored; coords[i] is where row i sits in the plane,
    // which is what the dissection cuts on. Returns false, with the factor
    // unusable, when a pivot comes out non-positive -- A was not positive
    // definite, or not numerically so.
    bool factor(const SparseMatrix &A, const std::vector<std::array<double, 2>> &coords);

    // X <- A^{-1} X for `nrhs` column-major right-hand sides of leading
    // dimension `ld` (>= n). All of them go through the tree together.
    void solve(double *X, int nrhs, int ld) const;

    int size() const { return n; }
    const Stats &getStats() const { return stats; }
    // perm[new] = old: the dissection ordering.
    const std::vector<int> &permutation() const { return perm; }

private:
    struct Supernode {
        int first = 0, size = 0;            // its pivots, in the new numbering
        int parent = -1;
        int subtreeFirst = 0;               // the subtree's variables are [subtreeFirst, first + size)
        std::vector<int> children;          // ascending, so fixed
        std::vector<int> rows;              // update rows: ascending, all >= first + size
        std::vector<int> relative;          // where each of `rows` sits in the parent's front
        std::vector<double> L;              // (size + rows.size()) x size, column-major
        double work = 0.0;                  // flops of this front alone
        double subtreeWork = 0.0;           // ... and of everything below it
        long long subtreeEntries = 0;       // entries of L in the subtree, for the solve
    };

    void order(const SparseMatrix &A, const std::vector<std::array<double, 2>> &coords);
    void symbolic(const SparseMatrix &A);
    void symbolicSubtree(int s);
    // The recursions take pointers rather than references: they spawn OpenMP
    // tasks, and a reference captured by a task is copied, not shared, unless
    // every construct says otherwise.
    void factorSubtree(int s, std::vector<std::vector<double>> *updates, std::atomic<bool> *bad);
    bool factorFront(int s, std::vector<std::vector<double>> *updates);
    void forwardSubtree(int s, double *Y, int nrhs, std::vector<std::vector<double>> *acc) const;
    void backwardSubtree(int s, double *Y, int nrhs) const;

    Options options;
    Stats stats;
    int n = 0;
    std::vector<int> perm, iperm;
    std::vector<Supernode> nodes;           // sorted by `first`: children before parents
    std::vector<int> roots;

    // The lower triangle of the permuted A, by columns: what the fronts are
    // assembled from.
    std::vector<int> lowPtr, lowRow;
    std::vector<double> lowVal;
};

}  // namespace shapedna

#endif  // __SHAPEDNA_SPARSECHOLESKY_HXX__
