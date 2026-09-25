#ifndef __SHAPEDNA_PARALLEL_HXX__
#define __SHAPEDNA_PARALLEL_HXX__

// OpenMP plumbing shared by the src/ShapeDNA sources.
//
// Every pragma goes through SDNA_OMP, which compiles to nothing without OpenMP,
// so a build without it is a correct build and only a slower one -- the same
// arrangement as the TMOP smoother (src/mesh/TMOP.cxx).
//
// Every parallel loop in the module is also written so that its answer does not
// depend on the thread count. A loop either writes disjoint outputs with each
// output computed in a fixed order (assembly gathers per row, the dense kernels
// update per column, the dissection numbers disjoint ranges), or it reduces
// over fixed-size chunks whose partial sums are then added in chunk order. None
// of them uses an OpenMP `reduction` clause, whose association follows the
// thread count. So the spectrum is bit-identical on one thread and on eight,
// and ShapeDNA's self-test checks exactly that.
//
// The same property is what lets the heavy loops use dynamic scheduling. On a
// machine with fast and slow cores (an M2 has four of each) a static split
// hands the slow cores as much as the fast ones and everyone waits for them;
// dynamic scheduling moves work, but no work item's arithmetic depends on
// which thread ran it.

#include <cstddef>

#ifdef _OPENMP
#  include <omp.h>
#  define SDNA_PRAGMA(x) _Pragma(#x)
#  define SDNA_OMP(x) SDNA_PRAGMA(omp x)
#else
#  define SDNA_OMP(x)
#endif

namespace shapedna {

inline bool haveOpenMP() {
#ifdef _OPENMP
    return true;
#else
    return false;
#endif
}

inline int maxThreads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

// Sets the OpenMP thread count for its own lifetime and puts the old one back
// afterwards, so a caller asking for one thread does not leave the whole
// process on one thread. `threads` <= 0 leaves the runtime's choice alone.
class ScopedThreads {
public:
    explicit ScopedThreads(int threads) {
#ifdef _OPENMP
        saved = omp_get_max_threads();
        if (threads > 0) omp_set_num_threads(threads);
#else
        (void)threads;
#endif
    }
    ~ScopedThreads() {
#ifdef _OPENMP
        omp_set_num_threads(saved);
#endif
    }
    ScopedThreads(const ScopedThreads &) = delete;
    ScopedThreads &operator=(const ScopedThreads &) = delete;

private:
#ifdef _OPENMP
    int saved = 1;
#endif
};

// Rows per chunk in the chunked loops over long vectors. Fixed, not derived from
// the thread count, because the chunk boundaries are where a reduction's
// association changes: keeping them fixed is what keeps the sums identical.
constexpr std::ptrdiff_t kRowChunk = 2048;

// A plain dot product with four interleaved accumulators, always combined in
// the same order. The four lanes let the compiler use SIMD without being given
// licence to reassociate (-ffast-math); the fixed combination keeps it
// deterministic.
inline double dotSerial(std::ptrdiff_t n, const double *x, const double *y) {
    double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
    std::ptrdiff_t i = 0;
    for (; i + 3 < n; i += 4) {
        s0 += x[i] * y[i];
        s1 += x[i + 1] * y[i + 1];
        s2 += x[i + 2] * y[i + 2];
        s3 += x[i + 3] * y[i + 3];
    }
    for (; i < n; ++i) s0 += x[i] * y[i];
    return (s0 + s1) + (s2 + s3);
}

}  // namespace shapedna

#endif  // __SHAPEDNA_PARALLEL_HXX__
