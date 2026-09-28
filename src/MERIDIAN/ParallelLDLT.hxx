#ifndef __PARALLEL_LDLT_HXX__
#define __PARALLEL_LDLT_HXX__

#include <atomic>
#include <memory>
#include <vector>

// eigen includes
#include <Eigen/Sparse>

// ---------------------------------------------------------------------------
// ParallelLDLT
//
// Eigen::SimplicialLDLT<SparseMatrix<double>> -- lower triangle, AMD ordering
// -- with the numeric factorisation spread over OpenMP threads, and a factor
// that is Eigen's own to the last bit. It is a drop-in: analyzePattern() once
// per sparsity pattern, factorize() whenever the values change, and solve(),
// info() and the rest are SimplicialLDLT's, reading the same L and D.
//
// It exists for Stage 6. LayoutEnergy::innerSolve() factorises a proxy Hessian
// of 16k-25k unknowns once per inner iteration, several hundred times a model,
// and on the corpus that is 85% of what MERIDIAN and TORSION spend.
//
// ### Why the answer cannot change
//
// Eigen factorises "up-looking" (LDL, Davis 2005): row k of L is a sparse
// triangular solve against the rows before it,
//
//     y = A(0:k, k)
//     for each i in the pattern of row k, in the order the etree walk gives:
//         l_ki = y_i / d_i
//         y_r -= L(r, i) y_i     for every entry L(r, i), r < k, of column i
//         d_k -= l_ki y_i
//
// and appends l_ki to column i. Every number in the factor is the result of
// that sequence of roundings, so the only freedom there is is *which rows run
// at the same time*, never the order of the operations inside one. Two facts
// about the elimination tree make that freedom real. The pattern of row k lies
// in the subtree of k, so rows in disjoint subtrees touch disjoint data and can
// run at once. And along a chain of the tree -- where the heavy rows are, and
// where subtree parallelism has nothing to offer (measured: it caps out at
// 1.3-1.6x on the Stage 6 Hessians) -- row k+1 needs row k only one entry at a
// time: its step at column i reads the entries of column i above it, the last
// of which row k writes at *its* step at column i. So consecutive rows can
// pipeline, each a step behind the one before, with a counter per column
// saying how far it is filled.
//
// That is the whole scheme. The symbolic pass replays Eigen's etree walk once
// per pattern and records, for every row, its steps (column i, and the slot p
// its entry goes to). The numeric pass hands out the bottom of the tree as
// whole subtrees, which one thread runs in index order without waiting on
// anything, and then the rows above them in index order; a step waits until
// column i holds p - Lp[i] entries (or, for the first entry of a column, until
// row i is done). Any schedule that respects those waits performs exactly the
// operations of the sequential loop on exactly the same values, so the factor
// does not depend on the thread count, and TestMERIDIAN --selftest checks it
// against Eigen's on every thread count.
//
// Inside a row, one more regrouping that leaves every operation where it was:
// a run of steps over consecutive columns of one supernode updates the rows
// below the run in column order, so those updates can wait until the run's
// own rows are done and then be applied row by row from a register instead
// of column by column through memory (factorRow() has the argument). Nearly
// all the work of the Stage 6 Hessians is in such runs, and the grouping alone
// makes one thread 2.2-2.3x faster than Eigen's loop.
//
// The expressions are written in the shapes Eigen writes them in, down to the
// pass-through calls. That is not decoration: clang contracts `d -= l * y`
// into one fused multiply-subtract, but Eigen's `d -= getDiag(l * getSymm(y))`
// rounds the product first, and the factors then differ in the last place.
// Mirroring the shapes makes whatever the compiler does to one happen to the
// other.
//
// ### How many threads
//
// OpenMP's count, capped at the machine's performance cores. A pipeline moves
// at the pace of its slowest stage, so a row of the chain that lands on an
// efficiency core holds up every row behind it: on an M2 (four of each) four
// threads factor the corpus Hessians 4.7-6.1x faster than Eigen, and before
// the cap eight were *slower* than one thread. One thread, or a call from
// inside a parallel region, runs the rows in plain sequence.
//
// And fewer when the machine is not ours. The same property makes the
// pipeline fragile when other processes want the cores: a thread that owes
// the next entry of a column and has been descheduled stalls every row behind
// it. Measured with three runs side by side, a model that factors in 55 s on
// Eigen's single thread took 334 s. So a wait that outlasts a few hundred
// spins gives up its core, and every factorisation measures how much of its
// threads' time they actually got to run -- thread CPU time over wall time,
// in which spinning counts, so it is 1 on a machine to itself. Below 0.85 the
// next factorisation uses half the threads, down to one; after eight clean
// ones it tries twice as many again. The factor is the same on any number of
// threads, so none of this can change an answer, only how fast it comes.
// ---------------------------------------------------------------------------
class ParallelLDLT : public Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> {
public:
    using Matrix = Eigen::SparseMatrix<double>;
    using Base = Eigen::SimplicialLDLT<Matrix>;

    ParallelLDLT() = default;
    ParallelLDLT(const ParallelLDLT &) = delete;
    ParallelLDLT &operator=(const ParallelLDLT &) = delete;

    void analyzePattern(const Matrix &a);
    void factorize(const Matrix &a);
    ParallelLDLT &compute(const Matrix &a) {
        analyzePattern(a);
        factorize(a);
        return *this;
    }

    // Threads for factorize(); 0 (the default) is OpenMP's count capped at the
    // performance cores, as above.
    void setThreads(int threads) { requestedThreads = threads; }
    int lastThreads() const { return usedThreads; }

    // The thread count 0 resolves to.
    static int defaultThreads();

private:
    struct Step {
        int node;   // column i of L
        int pos;    // where L(k, i) goes in L's value array
        int run;    // steps in the run this one starts; 0 inside a run
    };

    template <bool Concurrent>
    void factorRow(int k, double *y, double *ys, bool &failed);

    int n = 0;
    Eigen::Index sourceNonZeros = -1;   // the pattern analyzePattern() saw
    bool scheduled = false;

    // P A P^T, upper triangle, as Eigen builds it; apSource[q] is the index in
    // A's value array of ap's q-th value, so a new factorisation copies rather
    // than re-permutes.
    Matrix ap;
    std::vector<int> apSource;

    // Row k's steps are steps[rowStart[k] .. rowStart[k+1]).
    std::vector<int> rowStart;
    std::vector<Step> steps;

    // Work items: a subtree (its rows ascending) or a single row above them.
    std::vector<int> itemStart;
    std::vector<int> itemRows;
    int scheduledFor = 0;   // the thread count the items were cut for

    // Per column: entries of L written so far. Per row: finished.
    std::unique_ptr<std::atomic<int>[]> filled;
    std::unique_ptr<std::atomic<int>[]> done;

    // One dense accumulator per thread, all zero between rows, and each
    // thread's y_j for the columns of the run it is in.
    std::vector<std::vector<double>> scratch;
    std::vector<std::vector<double>> runScratch;
    int maxRun = 1;

    int requestedThreads = 0;
    int usedThreads = 1;
    // The adaptive count (see the class comment): the threads the next
    // factorisation gets, 0 until the first, and how many in a row have run
    // with their cores to themselves.
    int threadsNow = 0;
    int cleanStreak = 0;
    void adapt(double utilisation, double seconds, int cap);

    void buildItems(int threads);
};

#endif // __PARALLEL_LDLT_HXX__
