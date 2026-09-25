// Lanczos.cxx -- see Lanczos.hxx. Shift-invert block Lanczos in the B-inner
// product with full reorthogonalisation and Krylov-Schur (thick) restarts.
#include "Lanczos.hxx"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <sstream>

#include "DenseEigen.hxx"
#include "Parallel.hxx"

namespace shapedna {

namespace {

using Index = std::ptrdiff_t;

// A deterministic generator, so a run is reproducible from Options::seed.
struct Rng {
    unsigned long long s;
    explicit Rng(unsigned long long seed) : s(seed ? seed : 88172645463325252ULL) {}
    double next() {   // uniform in [-1, 1)
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        return 2.0 * (static_cast<double>(s >> 11) / 9007199254740992.0) - 1.0;
    }
};

// x . y over fixed chunks, the partial sums added in chunk order: the same
// answer for any thread count.
double dot(Index n, const double *x, const double *y) {
    const Index chunks = (n + kRowChunk - 1) / kRowChunk;
    std::vector<double> part(chunks);
    SDNA_OMP(parallel for schedule(dynamic, 1) if(chunks > 4))
    for (Index ch = 0; ch < chunks; ++ch) {
        const Index r0 = ch * kRowChunk, r1 = std::min(n, r0 + kRowChunk);
        part[ch] = dotSerial(r1 - r0, x + r0, y + r0);
    }
    double s = 0.0;
    for (double p : part) s += p;
    return s;
}

// Rows per block in the kernels below that update tall-skinny blocks: small
// enough that a block of every column involved stays in L1 while it is used.
constexpr Index kSubRows = 256;

// H (cols x b) = V(:, 0:cols)^T W. Each chunk of rows contributes a partial
// cols x b block, reading each basis column once; the partials are then added
// in chunk order, so the sums do not depend on how the chunks were shared out.
void innerProducts(Index n, const double *V, int cols, const double *W, int b, double *H) {
    const Index chunks = (n + kRowChunk - 1) / kRowChunk;
    const std::size_t blk = static_cast<std::size_t>(cols) * b;
    std::vector<double> part(static_cast<std::size_t>(chunks) * blk);
    SDNA_OMP(parallel for schedule(dynamic, 1) if(chunks > 1 && static_cast<double>(n) * blk > 1e5))
    for (Index ch = 0; ch < chunks; ++ch) {
        const Index r0 = ch * kRowChunk, len = std::min(n, r0 + kRowChunk) - r0;
        double *P = part.data() + ch * blk;
        for (int i = 0; i < cols; ++i) {
            const double *v = V + i * n + r0;
            for (int c = 0; c < b; ++c) P[i + c * cols] = dotSerial(len, v, W + c * n + r0);
        }
    }
    for (std::size_t t = 0; t < blk; ++t) {
        double s = 0.0;
        for (Index ch = 0; ch < chunks; ++ch) s += part[ch * blk + t];
        H[t] = s;
    }
}

// W -= V(:, 0:cols) H, a block of rows at a time, each basis column read once
// per block. Every entry of W loses its terms in column order.
void subtractCombination(Index n, const double *V, int cols, const double *H, double *W, int b) {
    const Index blocks = (n + kSubRows - 1) / kSubRows;
    SDNA_OMP(parallel for schedule(dynamic, 4) if(blocks > 8 && static_cast<double>(n) * cols * b > 1e5))
    for (Index blk = 0; blk < blocks; ++blk) {
        const Index r0 = blk * kSubRows, r1 = std::min(n, r0 + kSubRows);
        for (int i = 0; i < cols; ++i) {
            const double *v = V + i * n;
            for (int c = 0; c < b; ++c) {
                const double a = H[i + c * cols];
                if (a == 0.0) continue;
                double *w = W + c * n;
                for (Index r = r0; r < r1; ++r) w[r] -= v[r] * a;
            }
        }
    }
}

// out(:, 0:k) = V(:, 0:m) Y, Y m x k column-major. `out` may be V itself: each
// block of rows is read whole before any of it is written.
void combine(Index n, const double *V, int m, const double *Y, int k, double *out) {
    const Index blocks = (n + kSubRows - 1) / kSubRows;
    SDNA_OMP(parallel for schedule(dynamic, 4) if(blocks > 8))
    for (Index blk = 0; blk < blocks; ++blk) {
        const Index r0 = blk * kSubRows, len = std::min(n, r0 + kSubRows) - r0;
        std::vector<double> tmp(static_cast<std::size_t>(len) * k, 0.0);
        for (int i = 0; i < m; ++i) {
            const double *v = V + i * n + r0;
            for (int o = 0; o < k; ++o) {
                const double a = Y[i + static_cast<std::size_t>(o) * m];
                if (a == 0.0) continue;
                double *t = tmp.data() + o * len;
                for (Index r = 0; r < len; ++r) t[r] += v[r] * a;
            }
        }
        for (int o = 0; o < k; ++o)
            std::copy(tmp.begin() + o * len, tmp.begin() + (o + 1) * len, out + o * n + r0);
    }
}

// The state of one run: the basis and its projected matrix, and the scratch
// the orthogonalisation needs.
struct Krylov {
    const SparseMatrix &B;
    Index n;
    int b;
    std::vector<double> V;       // n x (basis + b)
    std::vector<double> BW;      // n x b
    std::vector<double> H;
    Rng rng;
    int deflations = 0;

    Krylov(const SparseMatrix &B_, int basis, int block, unsigned long long seed)
        : B(B_), n(B_.n), b(block),
          V(static_cast<std::size_t>(B_.n) * (basis + block)),
          BW(static_cast<std::size_t>(B_.n) * block), rng(seed) {}

    double *col(int j) { return V.data() + static_cast<std::size_t>(j) * n; }

    double bnorm(const double *w) {
        B.multiply(w, BW.data());
        return std::sqrt(std::max(0.0, dot(n, w, BW.data())));
    }

    // Classical Gram-Schmidt of the `nw` columns at W against the first `cols`
    // basis columns, in the B-inner product, run up to `passes` times. The
    // coefficients of all passes are summed into Hsum (cols x nw). When
    // `norms` is given it comes back with the B-norm of each column as it was
    // before anything was taken off it.
    //
    // A pass after the first runs only if the one before it lost more than
    // 1 - 1/sqrt(2) of some column's norm (Daniel, Gragg, Kaufman and Stewart,
    // Math. Comp. 30, 1976): below that the first pass left nothing the
    // second would remove above rounding. The test costs nothing, since B W is
    // needed at the start of the next pass anyway.
    void orthogonalize(int cols, double *W, int nw, int passes, std::vector<double> &Hsum,
                       double *norms = nullptr) {
        Hsum.assign(static_cast<std::size_t>(cols) * nw, 0.0);
        std::vector<double> before(nw), now(nw);
        for (int pass = 0; pass < passes; ++pass) {
            B.multiply(W, BW.data(), nw, static_cast<int>(n));
            for (int c = 0; c < nw; ++c)
                now[c] = std::sqrt(std::max(0.0, dot(n, W + c * n, BW.data() + c * n)));
            if (pass == 0 && norms) std::copy(now.begin(), now.end(), norms);
            if (cols == 0) break;
            if (pass > 0) {
                bool again = false;
                for (int c = 0; c < nw; ++c) again = again || now[c] < 0.7071 * before[c];
                if (!again) break;
            }
            before = now;
            H.resize(static_cast<std::size_t>(cols) * nw);
            innerProducts(n, V.data(), cols, BW.data(), nw, H.data());
            subtractCombination(n, V.data(), cols, H.data(), W, nw);
            for (std::size_t t = 0; t < Hsum.size(); ++t) Hsum[t] += H[t];
        }
    }

    // Makes the block of b columns starting at basis column `cols` -- already
    // B-orthogonal to everything before it -- B-orthonormal within itself, and
    // returns its R (b x b, upper triangular) with W = Q R.
    //
    // CholQR2 first (Yamamoto et al., ETNA 44, 2015): G = W^T B W = U^T U,
    // W <- W U^{-1}, twice. That is one block product with B and one b x b
    // Cholesky per pass, where the column-by-column route below costs 3b
    // single products. It is only safe on a well-conditioned block -- its
    // orthogonality error goes as eps cond(W)^2 -- so a block whose columns
    // have cancelled down to near nothing, or are near-dependent, goes the
    // careful way instead.
    void orthonormalizeBlock(int cols, const double *refNorm, std::vector<double> &R) {
        double *W = col(cols);
        std::vector<double> G(static_cast<std::size_t>(b) * b), U(G.size()), acc(G.size(), 0.0);
        for (int c = 0; c < b; ++c) acc[c + c * b] = 1.0;
        for (int pass = 0; pass < 2; ++pass) {
            B.multiply(W, BW.data(), b, static_cast<int>(n));
            innerProducts(n, W, b, BW.data(), b, G.data());
            for (int j = 0; j < b; ++j)
                for (int i = j + 1; i < b; ++i) G[i + j * b] = G[j + i * b] = 0.5 * (G[i + j * b] + G[j + i * b]);
            bool ok = denseCholesky(b, G.data());   // G = L L^T, U = L^T
            double dmin = 0.0, dmax = 0.0;
            if (ok) {
                dmin = dmax = G[0];
                for (int c = 0; c < b; ++c) {
                    const double d = G[c + c * b];
                    dmin = std::min(dmin, d);
                    dmax = std::max(dmax, d);
                    if (pass == 0 && d < 1e-4 * refNorm[c]) ok = false;
                }
                if (dmin < 1e-5 * dmax) ok = false;
            }
            if (!ok) {
                // After a first pass the columns are unit vectors, and that is
                // the scale their losses are measured against.
                std::vector<double> ones(b, 1.0), Rc;
                orthonormalizeColumns(cols, pass == 0 ? refNorm : ones.data(), Rc);
                multiplyUpper(Rc, acc, R);
                return;
            }
            for (int j = 0; j < b; ++j)
                for (int i = 0; i < b; ++i) U[i + j * b] = (i <= j) ? G[j + i * b] : 0.0;
            solveUpperRight(W, U);
            std::vector<double> next;
            multiplyUpper(U, acc, next);
            acc = std::move(next);
        }
        R = std::move(acc);
    }

    // C = A B for b x b upper-triangular A and B.
    void multiplyUpper(const std::vector<double> &A, const std::vector<double> &Bm,
                       std::vector<double> &C) const {
        C.assign(static_cast<std::size_t>(b) * b, 0.0);
        for (int j = 0; j < b; ++j)
            for (int i = 0; i <= j; ++i) {
                double s = 0.0;
                for (int k = i; k <= j; ++k) s += A[i + k * b] * Bm[k + j * b];
                C[i + j * b] = s;
            }
    }

    // W <- W U^{-1} for the b columns at W, U upper triangular: column c
    // becomes (w_c - sum_{c' < c} U(c', c) w_c') / U(c, c), in rows blocks.
    void solveUpperRight(double *W, const std::vector<double> &U) {
        const Index blocks = (n + kSubRows - 1) / kSubRows;
        SDNA_OMP(parallel for schedule(dynamic, 4) if(blocks > 8))
        for (Index blk = 0; blk < blocks; ++blk) {
            const Index r0 = blk * kSubRows, r1 = std::min(n, r0 + kSubRows);
            for (int c = 0; c < b; ++c) {
                double *w = W + c * n;
                for (int cp = 0; cp < c; ++cp) {
                    const double a = U[cp + c * b];
                    const double *q = W + cp * n;
                    for (Index r = r0; r < r1; ++r) w[r] -= q[r] * a;
                }
                const double inv = 1.0 / U[c + c * b];
                for (Index r = r0; r < r1; ++r) w[r] *= inv;
            }
        }
    }

    // The careful route, column by column. Each is orthogonalised against the
    // block's earlier columns twice. A column that has lost almost all of its
    // norm to cancellation -- less than 1e-4 of `refNorm`, what it had before
    // any orthogonalisation -- carries the rounding error of that cancellation
    // at a relative size that matters, so it is taken against the whole basis
    // again for as long as a pass still removes more than 30% of what is left
    // (the Daniel-Gragg-Kaufman-Stewart criterion). What is left at rounding
    // level of refNorm is a direction the Krylov space no longer has; it gets
    // R = 0 and a random replacement, orthogonal to everything.
    void orthonormalizeColumns(int cols, const double *refNorm, std::vector<double> &R) {
        R.assign(static_cast<std::size_t>(b) * b, 0.0);
        std::vector<double> h;
        for (int c = 0; c < b; ++c) {
            double *w = col(cols + c);
            if (c > 0) orthogonalizeWithin(cols, c, w, R);
            double nrm = bnorm(w);
            // Two passes are enough unless nearly everything cancelled ("twice
            // is enough", Kahan and Parlett); only then is it cleaned again.
            for (int extra = 0; extra < 3 && nrm > 0.0 && nrm < 1e-4 * refNorm[c]; ++extra) {
                const double before = nrm;
                orthogonalize(cols + c, w, 1, 1, h);
                for (int i = 0; i < c; ++i) R[i + c * b] += h[cols + i];
                nrm = bnorm(w);
                if (nrm > 0.7 * before) break;
            }
            if (nrm > 1e-13 * refNorm[c] && nrm > 0.0) {
                R[c + c * b] = nrm;
                const double inv = 1.0 / nrm;
                for (Index r = 0; r < n; ++r) w[r] *= inv;
                continue;
            }
            // Lost: replace it.
            ++deflations;
            R[c + c * b] = 0.0;
            for (int attempt = 0; attempt < 5; ++attempt) {
                for (Index r = 0; r < n; ++r) w[r] = rng.next();
                const double start = bnorm(w);
                orthogonalize(cols + c, w, 1, 2, h);
                nrm = bnorm(w);
                if (nrm > 1e-6 * start) break;
            }
            const double inv = nrm > 0.0 ? 1.0 / nrm : 0.0;
            for (Index r = 0; r < n; ++r) w[r] *= inv;
        }
    }

    // Column c of the new block against the block's columns 0..c-1, twice; the
    // coefficients go into R(0:c, c).
    void orthogonalizeWithin(int cols, int c, double *w, std::vector<double> &R) {
        const double *Q = col(cols);
        std::vector<double> hq(c);
        for (int pass = 0; pass < 2; ++pass) {
            B.multiply(w, BW.data());
            innerProducts(n, Q, c, BW.data(), 1, hq.data());
            subtractCombination(n, Q, c, hq.data(), w, 1);
            for (int i = 0; i < c; ++i) R[i + c * b] += hq[i];
        }
    }
};

}  // namespace

bool shiftInvertLanczos(const SparseMatrix &B, const SparseCholesky &K, double sigma,
                        const LanczosOptions &opt, LanczosResult &res) {
    res = LanczosResult();
    const Index n = B.n;
    const int b = std::max(1, opt.blockSize);
    const int nev = std::max(1, opt.nev);
    int M = opt.basisSize > 0 ? opt.basisSize : 2 * nev + 2 * b;
    M = std::max(M, nev + b);
    M = ((M + b - 1) / b) * b;
    res.basisSize = M;
    res.blockSize = b;
    if (static_cast<Index>(M) + b > n) {
        std::ostringstream os;
        os << "the problem has " << n << " unknowns, too few for a basis of " << M
           << " plus a block of " << b << "; solve it densely";
        res.message = os.str();
        return false;
    }

    Krylov kr(B, M, b, opt.seed);
    std::vector<double> T(static_cast<std::size_t>(M) * M, 0.0);
    std::vector<double> Hsum, R, Rlast(static_cast<std::size_t>(b) * b, 0.0);
    std::vector<double> refNorm(b);

    // The start block: random, B-orthonormalised.
    for (Index r = 0; r < n * b; ++r) kr.V[r] = kr.rng.next();
    for (int c = 0; c < b; ++c) refNorm[c] = kr.bnorm(kr.col(c));
    kr.orthonormalizeBlock(0, refNorm.data(), R);

    std::vector<double> w, Yall, theta(M), resid(M), Ysel;
    int kept = 0;
    for (;;) {
        // -- expand the basis from column `kept` to M ----------------------
        for (int j = kept; j < M; j += b) {
            double *W = kr.col(j + b);
            B.multiply(kr.col(j), W, b, static_cast<int>(n));
            K.solve(W, b, static_cast<int>(n));
            res.solves += b;

            const int cols = j + b;
            kr.orthogonalize(cols, W, b, 2, Hsum, refNorm.data());

            // Column block j of T = V^T B Op V, and its mirror. The diagonal
            // block is symmetric in exact arithmetic and is made so.
            for (int c = 0; c < b; ++c)
                for (int i = 0; i < j; ++i) {
                    const double h = Hsum[i + static_cast<std::size_t>(c) * cols];
                    T[i + static_cast<std::size_t>(j + c) * M] = h;
                    T[(j + c) + static_cast<std::size_t>(i) * M] = h;
                }
            for (int c = 0; c < b; ++c)
                for (int c2 = 0; c2 < b; ++c2)
                    T[(j + c2) + static_cast<std::size_t>(j + c) * M] =
                        0.5 * (Hsum[(j + c2) + static_cast<std::size_t>(c) * cols] +
                               Hsum[(j + c) + static_cast<std::size_t>(c2) * cols]);

            kr.orthonormalizeBlock(cols, refNorm.data(), R);
            if (j + b < M) {
                for (int c = 0; c < b; ++c)
                    for (int r = 0; r < b; ++r) {
                        const double x = R[r + static_cast<std::size_t>(c) * b];
                        T[(j + b + r) + static_cast<std::size_t>(j + c) * M] = x;
                        T[(j + c) + static_cast<std::size_t>(j + b + r) * M] = x;
                    }
            } else {
                Rlast = R;
            }
            ++res.blockSteps;
        }

        // -- Rayleigh-Ritz ------------------------------------------------
        if (!symmetricEigen(M, T.data(), w, &Yall)) {
            res.message = "the QL iteration on the projected matrix did not converge";
            return false;
        }
        // Largest theta first, i.e. smallest lambda first.
        for (int i = 0; i < M; ++i) theta[i] = w[M - 1 - i];
        auto ritz = [&](int i) { return Yall.data() + static_cast<std::size_t>(M - 1 - i) * M; };

        // Op V = V T + V_res Rlast E^T, E the last block, so the residual of
        // Ritz pair (theta, V y) is V_res Rlast y(M-b:M), whose B-norm is the
        // 2-norm of Rlast y(M-b:M).
        int good = 0;
        for (int i = 0; i < M; ++i) {
            const double *y = ritz(i);
            double s = 0.0;
            for (int r = 0; r < b; ++r) {
                double t = 0.0;
                for (int c = 0; c < b; ++c) t += Rlast[r + static_cast<std::size_t>(c) * b] * y[M - b + c];
                s += t * t;
            }
            resid[i] = std::sqrt(s);
            if (i < nev && resid[i] <= opt.tolerance * std::fabs(theta[i])) ++good;
        }
        if (good == nev) {
            res.converged = true;
            break;
        }
        if (res.restarts >= opt.maxRestarts) break;

        // -- thick restart ------------------------------------------------
        // Keep the nev wanted Ritz vectors and half of the rest: the unwanted
        // ones nearest the wanted end are what the next expansion converges
        // from, and dropping them throws that progress away.
        int keep = nev + (M - nev) / 2;
        keep = ((keep + b - 1) / b) * b;
        keep = std::min(keep, M - b);
        Ysel.resize(static_cast<std::size_t>(M) * keep);
        for (int i = 0; i < keep; ++i) std::copy(ritz(i), ritz(i) + M, Ysel.begin() + static_cast<std::size_t>(i) * M);
        combine(n, kr.V.data(), M, Ysel.data(), keep, kr.V.data());
        std::copy(kr.col(M), kr.col(M) + n * b, kr.col(keep));

        std::fill(T.begin(), T.end(), 0.0);
        for (int i = 0; i < keep; ++i) {
            T[i + static_cast<std::size_t>(i) * M] = theta[i];
            const double *y = ritz(i);
            for (int r = 0; r < b; ++r) {
                double s = 0.0;
                for (int c = 0; c < b; ++c) s += Rlast[r + static_cast<std::size_t>(c) * b] * y[M - b + c];
                T[(keep + r) + static_cast<std::size_t>(i) * M] = s;
                T[i + static_cast<std::size_t>(keep + r) * M] = s;
            }
        }
        kept = keep;
        ++res.restarts;
    }

    res.deflations = kr.deflations;
    res.values.resize(nev);
    res.residuals.resize(nev);
    for (int i = 0; i < nev; ++i) {
        res.values[i] = sigma + 1.0 / theta[i];
        res.residuals[i] = resid[i] / std::fabs(theta[i]);
    }
    if (opt.vectors) {
        Ysel.resize(static_cast<std::size_t>(M) * nev);
        for (int i = 0; i < nev; ++i)
            std::copy(Yall.data() + static_cast<std::size_t>(M - 1 - i) * M,
                      Yall.data() + static_cast<std::size_t>(M - i) * M,
                      Ysel.begin() + static_cast<std::size_t>(i) * M);
        res.vectors.assign(static_cast<std::size_t>(n) * nev, 0.0);
        combine(n, kr.V.data(), M, Ysel.data(), nev, res.vectors.data());
    }
    if (!res.converged) {
        std::ostringstream os;
        os << "not converged after " << res.restarts << " restarts";
        res.message = os.str();
    }
    return res.converged;
}

}  // namespace shapedna
