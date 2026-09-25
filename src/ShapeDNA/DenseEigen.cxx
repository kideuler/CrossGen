// DenseEigen.cxx -- see DenseEigen.hxx.
#include "DenseEigen.hxx"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace shapedna {

namespace {

// Householder reduction of the symmetric matrix held in V (row-major n x n,
// V[i*n + j]) to tridiagonal form, accumulating the orthogonal transformation
// in V. On return d is the diagonal and e the subdiagonal in e[1..n-1], e[0] = 0.
// The EISPACK tred2 algorithm, scaled rows and all.
void tridiagonalize(int n, std::vector<double> &V, std::vector<double> &d,
                    std::vector<double> &e) {
    auto at = [&](int i, int j) -> double & { return V[static_cast<std::size_t>(i) * n + j]; };

    for (int j = 0; j < n; ++j) d[j] = at(n - 1, j);

    for (int i = n - 1; i > 0; --i) {
        // Scale to avoid under/overflow in the Householder vector.
        double scale = 0.0, h = 0.0;
        for (int k = 0; k < i; ++k) scale += std::fabs(d[k]);
        if (scale == 0.0) {
            e[i] = d[i - 1];
            for (int j = 0; j < i; ++j) {
                d[j] = at(i - 1, j);
                at(i, j) = 0.0;
                at(j, i) = 0.0;
            }
        } else {
            for (int k = 0; k < i; ++k) {
                d[k] /= scale;
                h += d[k] * d[k];
            }
            double f = d[i - 1];
            double g = std::sqrt(h);
            if (f > 0) g = -g;
            e[i] = scale * g;
            h -= f * g;
            d[i - 1] = f - g;
            for (int j = 0; j < i; ++j) e[j] = 0.0;

            // Apply the similarity transformation to the remaining columns.
            for (int j = 0; j < i; ++j) {
                f = d[j];
                at(j, i) = f;
                g = e[j] + at(j, j) * f;
                for (int k = j + 1; k <= i - 1; ++k) {
                    g += at(k, j) * d[k];
                    e[k] += at(k, j) * f;
                }
                e[j] = g;
            }
            f = 0.0;
            for (int j = 0; j < i; ++j) {
                e[j] /= h;
                f += e[j] * d[j];
            }
            const double hh = f / (h + h);
            for (int j = 0; j < i; ++j) e[j] -= hh * d[j];
            for (int j = 0; j < i; ++j) {
                f = d[j];
                g = e[j];
                for (int k = j; k <= i - 1; ++k) at(k, j) -= (f * e[k] + g * d[k]);
                d[j] = at(i - 1, j);
                at(i, j) = 0.0;
            }
        }
        d[i] = h;
    }

    // Accumulate the transformations.
    for (int i = 0; i < n - 1; ++i) {
        at(n - 1, i) = at(i, i);
        at(i, i) = 1.0;
        const double h = d[i + 1];
        if (h != 0.0) {
            for (int k = 0; k <= i; ++k) d[k] = at(k, i + 1) / h;
            for (int j = 0; j <= i; ++j) {
                double g = 0.0;
                for (int k = 0; k <= i; ++k) g += at(k, i + 1) * at(k, j);
                for (int k = 0; k <= i; ++k) at(k, j) -= g * d[k];
            }
        }
        for (int k = 0; k <= i; ++k) at(k, i + 1) = 0.0;
    }
    for (int j = 0; j < n; ++j) {
        d[j] = at(n - 1, j);
        at(n - 1, j) = 0.0;
    }
    at(n - 1, n - 1) = 1.0;
    e[0] = 0.0;
}

// Implicit QL with Wilkinson shifts on the tridiagonal (d, e) left by
// tridiagonalize, rotating the columns of V along. EISPACK tql2. Returns false
// if an eigenvalue has not converged after 60 sweeps.
bool tridiagonalQL(int n, std::vector<double> &V, std::vector<double> &d,
                   std::vector<double> &e) {
    auto at = [&](int i, int j) -> double & { return V[static_cast<std::size_t>(i) * n + j]; };

    for (int i = 1; i < n; ++i) e[i - 1] = e[i];
    e[n - 1] = 0.0;

    double f = 0.0, tst1 = 0.0;
    const double eps = std::numeric_limits<double>::epsilon();
    for (int l = 0; l < n; ++l) {
        // Find a negligible subdiagonal element to split at.
        tst1 = std::max(tst1, std::fabs(d[l]) + std::fabs(e[l]));
        int m = l;
        while (m < n) {
            if (std::fabs(e[m]) <= eps * tst1) break;
            ++m;
        }
        if (m == n) m = n - 1;   // e[n-1] is zero; only a NaN gets here

        if (m > l) {
            int iter = 0;
            do {
                if (++iter > 60) return false;
                // The shift: the eigenvalue of the leading 2x2 nearer d[l].
                double g = d[l];
                double p = (d[l + 1] - g) / (2.0 * e[l]);
                double r = std::hypot(p, 1.0);
                if (p < 0) r = -r;
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                const double dl1 = d[l + 1];
                double h = g - d[l];
                for (int i = l + 2; i < n; ++i) d[i] -= h;
                f += h;

                // The QL sweep, chasing the bulge from m up to l.
                p = d[m];
                double c = 1.0, c2 = c, c3 = c;
                const double el1 = e[l + 1];
                double s = 0.0, s2 = 0.0;
                for (int i = m - 1; i >= l; --i) {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = std::hypot(p, e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);
                    for (int k = 0; k < n; ++k) {
                        h = at(k, i + 1);
                        at(k, i + 1) = s * at(k, i) + c * h;
                        at(k, i) = c * at(k, i) - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
            } while (std::fabs(e[l]) > eps * tst1);
        }
        d[l] += f;
        e[l] = 0.0;
    }
    return true;
}

}  // namespace

bool symmetricEigen(int n, const double *a, std::vector<double> &values,
                    std::vector<double> *vectors) {
    values.assign(n, 0.0);
    if (vectors) vectors->assign(static_cast<std::size_t>(n) * n, 0.0);
    if (n == 0) return true;
    if (n == 1) {
        values[0] = a[0];
        if (vectors) (*vectors)[0] = 1.0;
        return std::isfinite(a[0]);
    }

    // Symmetrise from the lower triangle. Row-major and column-major coincide
    // for a symmetric matrix, so this is also the row-major copy the two
    // EISPACK routines index.
    std::vector<double> V(static_cast<std::size_t>(n) * n);
    for (int j = 0; j < n; ++j)
        for (int i = j; i < n; ++i) {
            const double x = a[i + static_cast<std::size_t>(j) * n];
            V[static_cast<std::size_t>(i) * n + j] = x;
            V[static_cast<std::size_t>(j) * n + i] = x;
        }

    std::vector<double> d(n), e(n);
    tridiagonalize(n, V, d, e);
    if (!tridiagonalQL(n, V, d, e)) return false;

    // Sort ascending, carrying the columns of V along.
    std::vector<int> order(n);
    for (int i = 0; i < n; ++i) order[i] = i;
    std::stable_sort(order.begin(), order.end(), [&](int x, int y) { return d[x] < d[y]; });
    for (int i = 0; i < n; ++i) values[i] = d[order[i]];
    if (vectors) {
        for (int i = 0; i < n; ++i) {
            const int src = order[i];
            double *out = vectors->data() + static_cast<std::size_t>(i) * n;
            for (int k = 0; k < n; ++k) out[k] = V[static_cast<std::size_t>(k) * n + src];
        }
    }
    return true;
}

bool denseCholesky(int n, double *a) {
    for (int k = 0; k < n; ++k) {
        double *ck = a + static_cast<std::size_t>(k) * n;
        const double dkk = ck[k];
        if (!(dkk > 0.0)) return false;
        const double r = std::sqrt(dkk);
        ck[k] = r;
        const double inv = 1.0 / r;
        for (int i = k + 1; i < n; ++i) ck[i] *= inv;
        for (int j = k + 1; j < n; ++j) {
            const double s = ck[j];
            if (s == 0.0) continue;
            double *cj = a + static_cast<std::size_t>(j) * n;
            for (int i = j; i < n; ++i) cj[i] -= ck[i] * s;
        }
    }
    for (int j = 1; j < n; ++j)
        for (int i = 0; i < j; ++i) a[i + static_cast<std::size_t>(j) * n] = 0.0;
    return true;
}

bool generalizedSymmetricEigen(int n, const double *A, const double *B,
                               std::vector<double> &values,
                               std::vector<double> *vectors) {
    const std::size_t N = static_cast<std::size_t>(n);
    std::vector<double> L(B, B + N * N);
    if (!denseCholesky(n, L.data())) return false;

    // C = L^{-1} A L^{-T}, symmetric. First W = L^{-1} A (forward substitution
    // on every column of A, symmetrised from its lower triangle), then
    // C = L^{-1} W^T, using C = C^T.
    std::vector<double> W(N * N);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            W[i + j * N] = (i >= j) ? A[i + j * N] : A[j + i * N];
    auto forward = [&](double *x) {
        for (int k = 0; k < n; ++k) {
            x[k] /= L[k + k * N];
            const double xk = x[k];
            for (int i = k + 1; i < n; ++i) x[i] -= L[i + k * N] * xk;
        }
    };
    for (int j = 0; j < n; ++j) forward(W.data() + j * N);
    std::vector<double> C(N * N);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) C[i + j * N] = W[j + i * N];
    for (int j = 0; j < n; ++j) forward(C.data() + j * N);

    std::vector<double> Y;
    if (!symmetricEigen(n, C.data(), values, vectors ? &Y : nullptr)) return false;
    if (vectors) {
        // x = L^{-T} y, which is B-orthonormal because y is orthonormal.
        for (int j = 0; j < n; ++j) {
            double *y = Y.data() + j * N;
            for (int k = n - 1; k >= 0; --k) {
                double s = y[k];
                for (int i = k + 1; i < n; ++i) s -= L[i + k * N] * y[i];
                y[k] = s / L[k + k * N];
            }
        }
        *vectors = std::move(Y);
    }
    return true;
}

}  // namespace shapedna
