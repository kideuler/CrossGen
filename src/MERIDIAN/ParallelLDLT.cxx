#include "MERIDIAN/ParallelLDLT.hxx"

#include <algorithm>
#include <chrono>
#include <ctime>
#include <vector>

#include <sched.h>

#ifdef _OPENMP
#  include <omp.h>
#endif
#if defined(__APPLE__)
#  include <sys/sysctl.h>
#endif

namespace {

// Eigen's getDiag() and getSymm() for a real scalar: the identity, but a call,
// and where the calls are is where the compiler may or may not fuse a multiply
// into the subtraction that follows it. See the class comment.
inline double passThrough(double x) { return x; }

// A spinning thread's hint to the core that it is waiting.
inline void relax() {
#if defined(__aarch64__) || defined(__arm64__)
    __asm__ __volatile__("yield");
#elif defined(__x86_64__) || defined(__i386__)
    __builtin_ia32_pause();
#endif
}

// Waits until ready() holds. A few hundred spins first -- on a machine to
// itself the row being waited on is a step or two ahead and the wait is well
// under a microsecond -- and then the core is given up at every look, because
// a wait longer than that usually means the thread that owes the entry is not
// running: another process has its core, and spinning on only keeps it off.
template <typename Ready>
inline void waitUntil(Ready ready) {
    for (int spin = 0; spin < 512; ++spin) {
        if (ready()) return;
        relax();
    }
    while (!ready()) sched_yield();
}

// This thread's CPU time and the wall clock, in seconds.
inline double threadSeconds() {
    timespec ts;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts);
    return static_cast<double>(ts.tv_sec) + 1e-9 * static_cast<double>(ts.tv_nsec);
}
inline double wallSeconds() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// A factorisation this short says too little about the machine to act on.
constexpr double kMinAdaptSeconds = 2e-3;

// How finely the bottom of the tree is cut: about this many subtrees per
// thread, so that the ones that finish early find more. Four to eight were
// equally good on the corpus Hessians; fewer leaves too much of the tree to
// the pipeline, more makes subtrees too small to be worth a work item.
constexpr int kSubtreesPerThread = 6;

// Below this many unknowns per thread the waits cost more than they save.
constexpr int kRowsPerThread = 256;

} // namespace

int ParallelLDLT::defaultThreads() {
#ifdef _OPENMP
    if (omp_in_parallel()) return 1;
    int threads = omp_get_max_threads();
#  if defined(__APPLE__)
    static const int performanceCores = [] {
        int cores = 0;
        size_t size = sizeof(cores);
        if (sysctlbyname("hw.perflevel0.physicalcpu", &cores, &size, nullptr, 0) != 0) cores = 0;
        return cores;
    }();
    if (performanceCores > 0) threads = std::min(threads, performanceCores);
#  endif
    return std::max(1, threads);
#else
    return 1;
#endif
}

// ---------------------------------------------------------------------------
// analyzePattern()
//
// Eigen's own analysis -- the AMD ordering, the elimination tree, the column
// counts and L's column pointers -- and then, once per pattern, everything
// about the numeric pass that does not depend on the values: where each value
// of P A P^T comes from, and each row's steps, found by replaying the etree
// walk of SimplicialCholeskyBase::factorize_preordered() on the same matrix in
// the same order. The walk also fixes L's row indices, which are therefore
// written here rather than at every factorisation.
// ---------------------------------------------------------------------------
void ParallelLDLT::analyzePattern(const Matrix &a) {
    Base::analyzePattern(a);
    scheduled = false;
    scheduledFor = 0;
    n = static_cast<int>(a.rows());
    sourceNonZeros = a.nonZeros();
    if (!a.isCompressed() || m_info != Eigen::Success || n == 0) return;

    // ap exactly as Base::factorize() builds it, with each value replaced by
    // the index of the value of A it copies.
    Matrix index = a;
    for (Eigen::Index s = 0; s < index.nonZeros(); ++s) index.valuePtr()[s] = static_cast<double>(s);
    ap.resize(n, n);
    Eigen::internal::permute_symm_to_symm<Eigen::Lower, Eigen::Upper, false>(
        index, ap, m_P.size() > 0 ? m_P.indices().data() : nullptr);
    apSource.resize(static_cast<size_t>(ap.nonZeros()));
    for (Eigen::Index q = 0; q < ap.nonZeros(); ++q) apSource[q] = static_cast<int>(ap.valuePtr()[q]);

    const int *Lp = m_matrix.outerIndexPtr();
    int *Li = m_matrix.innerIndexPtr();
    const int *parent = m_parent.data();
    std::vector<int> tags(n, -1), pattern(n), count(n, 0);
    rowStart.assign(n + 1, 0);
    steps.clear();
    steps.reserve(static_cast<size_t>(m_matrix.nonZeros()));
    for (int k = 0; k < n; ++k) {
        int top = n;
        tags[k] = k;
        for (Matrix::InnerIterator it(ap, k); it; ++it) {
            int i = static_cast<int>(it.index());
            if (i > k) continue;
            int len = 0;
            for (; tags[i] != k; i = parent[i]) {
                pattern[len++] = i;
                tags[i] = k;
            }
            while (len > 0) pattern[--top] = pattern[--len];
        }
        for (; top < n; ++top) {
            const int i = pattern[top];
            const int p = Lp[i] + count[i]++;
            Li[p] = k;
            steps.push_back({i, p, 1});
        }
        rowStart[k + 1] = static_cast<int>(steps.size());
    }
    // The replay has to land every entry exactly where Eigen's column counts
    // put it; if it ever did not, the Eigen factorisation is what runs.
    for (int i = 0; i < n; ++i) {
        if (Lp[i] + count[i] != Lp[i + 1]) return;
    }

    // Runs. Column j continues into column j+1 when its rows are j+1 followed
    // by exactly the rows of column j+1 -- consecutive columns of one
    // supernode, which is how AMD lays out the (u, v) pairs of a vertex and
    // the separators above them. A row that steps through such columns one
    // after another is given them as one run; see factorRow().
    std::vector<char> continues(n, 0);
    for (int j = 0; j + 1 < n; ++j) {
        const int len = Lp[j + 1] - Lp[j], next = Lp[j + 2] - Lp[j + 1];
        if (len != next + 1 || Li[Lp[j]] != j + 1) continue;
        bool same = true;
        for (int m = 0; m < next && same; ++m) same = Li[Lp[j] + 1 + m] == Li[Lp[j + 1] + m];
        continues[j] = same ? 1 : 0;
    }
    maxRun = 1;
    for (int k = 0; k < n; ++k) {
        for (int s = rowStart[k]; s < rowStart[k + 1];) {
            int e = s + 1;
            while (e < rowStart[k + 1] && steps[e].node == steps[e - 1].node + 1 && continues[steps[e - 1].node]) ++e;
            steps[s].run = e - s;
            for (int t = s + 1; t < e; ++t) steps[t].run = 0;
            maxRun = std::max(maxRun, e - s);
            s = e;
        }
    }

    filled.reset(new std::atomic<int>[n]);
    done.reset(new std::atomic<int>[n]);
    scheduled = true;
}

// ---------------------------------------------------------------------------
// buildItems()
//
// The work items for a thread count: the subtrees at the bottom of the etree
// whose work is under 1/(threads * kSubtreesPerThread) of the total, each as
// one item run in index order, roots ascending; then every row above them as
// an item of its own, in index order. Handing all the subtrees out first is
// what makes the waits deadlock-free: when a thread takes a row from the top,
// every subtree has an owner that never waits on anyone, and every row below
// it in the top part has an owner too, so the lowest unfinished row can always
// proceed.
// ---------------------------------------------------------------------------
void ParallelLDLT::buildItems(int threads) {
    const int *Lp = m_matrix.outerIndexPtr();
    const int *apOuter = ap.outerIndexPtr();
    const int *parent = m_parent.data();

    // Work of a subtree: per row, its scatter and, per step, the prefix of the
    // column it walks. Children come before parents in index order.
    std::vector<double> subtree(n, 0.0);
    double total = 0.0;
    for (int k = 0; k < n; ++k) {
        double w = 4.0 + (apOuter[k + 1] - apOuter[k]);
        for (int s = rowStart[k]; s < rowStart[k + 1]; ++s) w += 4.0 + (steps[s].pos - Lp[steps[s].node]);
        subtree[k] += w;
        total += w;
        if (parent[k] >= 0) subtree[parent[k]] += subtree[k];
    }
    const double cut = total / static_cast<double>(threads * kSubtreesPerThread);

    // rootOf[k]: the root of the subtree item row k belongs to, or -1 for a
    // row of the top part. Parents come after children, so walking down from
    // n - 1 sees a parent's verdict first.
    std::vector<int> rootOf(n, -1);
    for (int k = n - 1; k >= 0; --k) {
        const int p = parent[k];
        if (p >= 0 && rootOf[p] >= 0) rootOf[k] = rootOf[p];
        else if (subtree[k] <= cut) rootOf[k] = k;
    }

    std::vector<int> itemOfRoot(n, -1);
    int subtrees = 0;
    for (int k = 0; k < n; ++k) {
        if (rootOf[k] == k) itemOfRoot[k] = subtrees++;
    }
    itemStart.assign(subtrees + 1, 0);
    for (int k = 0; k < n; ++k) {
        if (rootOf[k] >= 0) ++itemStart[itemOfRoot[rootOf[k]] + 1];
    }
    for (int t = 0; t < subtrees; ++t) itemStart[t + 1] += itemStart[t];
    itemRows.assign(n, 0);
    std::vector<int> next(itemStart.begin(), itemStart.end() - 1);
    for (int k = 0; k < n; ++k) {
        if (rootOf[k] >= 0) itemRows[next[itemOfRoot[rootOf[k]]]++] = k;
    }
    int at = itemStart[subtrees];
    for (int k = 0; k < n; ++k) {
        if (rootOf[k] < 0) {
            itemRows[at++] = k;
            itemStart.push_back(at);
        }
    }
    scheduledFor = threads;
}

// ---------------------------------------------------------------------------
// factorRow()
//
// Row k of L and d_k: SimplicialCholeskyBase::factorize_preordered()'s loop
// body, with the steps read from the replay instead of re-derived and, when
// Concurrent, a wait before each step and a release after it. `y` is the
// calling thread's accumulator; it is zero on entry and left zero, since every
// entry the row touches is a step of the row, read and cleared in turn.
//
// A single step is Eigen's statement for statement. A run of steps over
// columns j0 .. j1 of one supernode (see analyzePattern()) does the same
// operations in a different grouping. Column j of the run holds rows j+1 ..
// j1 first -- the run's own -- and then the same rows as column j1, so the
// step at column j splits into an update of the run's rows, which the next
// steps read, and an update of the rest, which nothing in the run reads. The
// first part is done step by step, as it was, and contiguously, since those
// rows are consecutive. The second is deferred to the end of the run and done
// row by row: each such y_r then takes the run's multiply-subtracts in column
// order, j0 first, the order the steps would have applied them in, but from a
// register, eight rows at a time, instead of in one pass over memory per
// column. On the Stage 6 Hessians nearly all the work is in runs of 10 to 200
// columns, and the grouping alone makes a thread 2.2x faster than Eigen's.
// ---------------------------------------------------------------------------
template <bool Concurrent>
void ParallelLDLT::factorRow(int k, double *y, double *ys, bool &failed) {
    const int *Lp = m_matrix.outerIndexPtr();
    const int *Li = m_matrix.innerIndexPtr();
    double *Lx = m_matrix.valuePtr();
    double *D = m_diag.data();
    const int *apOuter = ap.outerIndexPtr();
    const int *apInner = ap.innerIndexPtr();
    const double *apValue = ap.valuePtr();

    // Column i's entries up to p are in, and so is d_i.
    auto await = [&](int i, int p) {
        if (!Concurrent) return;
        if (p > Lp[i]) {
            const int need = p - Lp[i];
            waitUntil([&] { return filled[i].load(std::memory_order_acquire) >= need; });
        } else {
            waitUntil([&] { return done[i].load(std::memory_order_acquire) != 0; });
        }
    };

    for (int q = apOuter[k]; q < apOuter[k + 1]; ++q) {
        const int i = apInner[q];
        if (i <= k) y[i] += passThrough(apValue[q]);
    }
    double d = passThrough(y[k]) * m_shiftScale + m_shiftOffset;
    y[k] = 0.0;
    for (int s = rowStart[k]; s < rowStart[k + 1];) {
        const int len = steps[s].run;
        if (len <= 1) {
            const int i = steps[s].node;
            const int p = steps[s].pos;
            await(i, p);
            const double yi = y[i];
            y[i] = 0.0;
            const double lki = yi / passThrough(D[i]);
            for (int q = Lp[i]; q < p; ++q) y[Li[q]] -= passThrough(Lx[q]) * yi;
            d -= passThrough(lki * passThrough(yi));
            Lx[p] = lki;
            if (Concurrent) filled[i].store(p - Lp[i] + 1, std::memory_order_release);
            ++s;
            continue;
        }

        const int j0 = steps[s].node, j1 = j0 + len - 1;
        for (int t = 0; t < len; ++t) {
            const int j = j0 + t;
            const int p = steps[s + t].pos;
            await(j, p);
            const double yj = y[j];
            y[j] = 0.0;
            const double lkj = yj / passThrough(D[j]);
            // The run's own rows, j+1 .. j1: the first j1 - j entries.
            double *yRun = y + j + 1;
            const double *lRun = Lx + Lp[j];
            for (int q = 0; q < j1 - j; ++q) yRun[q] -= passThrough(lRun[q]) * yj;
            d -= passThrough(lkj * passThrough(yj));
            Lx[p] = lkj;
            if (Concurrent) filled[j].store(p - Lp[j] + 1, std::memory_order_release);
            ys[t] = yj;
        }
        // The rest: column j1's rows above k, which every column of the run
        // holds at the same place after its own run rows.
        // The loops over t are each row's chain of multiply-subtracts, and
        // are kept scalar: vectorised as reductions they would round every
        // product before subtracting it, which Eigen's loop does not (see
        // LayoutEnergy::energy() for where that happened).
        const int M = steps[s + len - 1].pos - Lp[j1];
        const int *rows = Li + Lp[j1];
        int m = 0;
        for (; m + 8 <= M; m += 8) {
            double acc[8];
            for (int u = 0; u < 8; ++u) acc[u] = y[rows[m + u]];
#if defined(__clang__)
#pragma clang loop vectorize(disable) interleave(disable)
#endif
            for (int t = 0; t < len; ++t) {
                const double *c = Lx + Lp[j0 + t] + (len - 1 - t) + m;
                const double yt = ys[t];
                for (int u = 0; u < 8; ++u) acc[u] -= passThrough(c[u]) * yt;
            }
            for (int u = 0; u < 8; ++u) y[rows[m + u]] = acc[u];
        }
        for (; m < M; ++m) {
            double acc = y[rows[m]];
#if defined(__clang__)
#pragma clang loop vectorize(disable) interleave(disable)
#endif
            for (int t = 0; t < len; ++t) acc -= passThrough(Lx[Lp[j0 + t] + (len - 1 - t) + m]) * ys[t];
            y[rows[m]] = acc;
        }
        s += len;
    }
    D[k] = d;
    if (d == 0.0) failed = true;
    if (Concurrent) done[k].store(1, std::memory_order_release);
}

// ---------------------------------------------------------------------------
// factorize()
// ---------------------------------------------------------------------------
void ParallelLDLT::factorize(const Matrix &a) {
    if (!scheduled || a.rows() != n || a.nonZeros() != sourceNonZeros || !a.isCompressed()) {
        Base::factorize(a);
        usedThreads = 1;
        return;
    }

    const double *value = a.valuePtr();
    double *apValue = ap.valuePtr();
    const size_t m = apSource.size();
    for (size_t q = 0; q < m; ++q) apValue[q] = value[apSource[q]];
    m_diag.resize(n);

    int cap = requestedThreads > 0 ? requestedThreads : defaultThreads();
#ifdef _OPENMP
    if (omp_in_parallel()) cap = 1;
#else
    cap = 1;
#endif
    cap = std::max(1, std::min(cap, n / kRowsPerThread));
    if (threadsNow <= 0 || threadsNow > cap) threadsNow = cap;
    const int threads = threadsNow;
    usedThreads = threads;
    if (static_cast<int>(scratch.size()) < threads) {
        scratch.resize(threads);
        runScratch.resize(threads);
    }
    for (int t = 0; t < threads; ++t) {
        if (static_cast<int>(scratch[t].size()) != n) scratch[t].assign(n, 0.0);
        if (static_cast<int>(runScratch[t].size()) < maxRun) runScratch[t].assign(maxRun, 0.0);
    }

    bool failed = false;
    const double wall0 = wallSeconds();
    double utilisation = 1.0;
    if (threads == 1) {
        const double cpu0 = threadSeconds();
        double *y = scratch[0].data();
        double *ys = runScratch[0].data();
        for (int k = 0; k < n; ++k) factorRow<false>(k, y, ys, failed);
        const double wall = wallSeconds() - wall0;
        if (wall > 0.0) utilisation = (threadSeconds() - cpu0) / wall;
    } else {
#ifdef _OPENMP
        if (scheduledFor != threads) buildItems(threads);
        for (int i = 0; i < n; ++i) {
            filled[i].store(0, std::memory_order_relaxed);
            done[i].store(0, std::memory_order_relaxed);
        }
        std::atomic<int> nextItem(0);
        std::atomic<bool> anyFailed(false);
        const int items = static_cast<int>(itemStart.size()) - 1;
        // Each thread's share of the time it was busy that it actually ran,
        // measured from the start of the region to its last row.
        std::vector<double> ran(threads, 1.0);
#pragma omp parallel num_threads(threads)
        {
            const int me = omp_get_thread_num();
            const double w0 = wallSeconds(), c0 = threadSeconds();
            double *y = scratch[me].data();
            double *ys = runScratch[me].data();
            bool mine = false;
            for (;;) {
                const int item = nextItem.fetch_add(1, std::memory_order_relaxed);
                if (item >= items) break;
                for (int r = itemStart[item]; r < itemStart[item + 1]; ++r) factorRow<true>(itemRows[r], y, ys, mine);
            }
            const double w = wallSeconds() - w0;
            if (w > 0.0) ran[me] = (threadSeconds() - c0) / w;
            if (mine) anyFailed.store(true, std::memory_order_relaxed);
        }
        failed = anyFailed.load();
        double sum = 0.0;
        for (double r : ran) sum += r;
        utilisation = sum / threads;
#endif
    }
    adapt(utilisation, wallSeconds() - wall0, cap);
    m_info = failed ? Eigen::NumericalIssue : Eigen::Success;
    m_factorizationIsOk = true;
}

void ParallelLDLT::adapt(double utilisation, double seconds, int cap) {
    if (seconds < kMinAdaptSeconds) return;
    if (utilisation < 0.85) {
        threadsNow = std::max(1, threadsNow / 2);
        cleanStreak = 0;
    } else if (utilisation > 0.95 && threadsNow < cap) {
        if (++cleanStreak >= 8) {
            threadsNow = std::min(cap, threadsNow * 2);
            cleanStreak = 0;
        }
    } else {
        cleanStreak = 0;
    }
}
