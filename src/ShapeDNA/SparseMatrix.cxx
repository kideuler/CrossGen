// SparseMatrix.cxx -- see SparseMatrix.hxx.
#include "SparseMatrix.hxx"

#include <algorithm>
#include <cstddef>

#include "Parallel.hxx"

namespace shapedna {

void SparseMatrix::multiply(const double *x, double *y) const {
    SDNA_OMP(parallel for schedule(dynamic, 512) if(n > 4096))
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) s += val[k] * x[col[k]];
        y[i] = s;
    }
}

void SparseMatrix::multiply(const double *X, double *Y, int nvec, int ld) const {
    if (nvec == 1) {
        multiply(X, Y);
        return;
    }
    const std::ptrdiff_t L = ld;
    SDNA_OMP(parallel for schedule(dynamic, 256) if(n > 1024))
    for (int i = 0; i < n; ++i) {
        const int k0 = rowPtr[i], k1 = rowPtr[i + 1];
        for (int v = 0; v < nvec; ++v) {
            const double *x = X + v * L;
            double s = 0.0;
            for (int k = k0; k < k1; ++k) s += val[k] * x[col[k]];
            Y[i + v * L] = s;
        }
    }
}

double SparseMatrix::entry(int i, int j) const {
    const int *b = col.data() + rowPtr[i];
    const int *e = col.data() + rowPtr[i + 1];
    const int *p = std::lower_bound(b, e, j);
    if (p == e || *p != j) return 0.0;
    return val[p - col.data()];
}

}  // namespace shapedna
