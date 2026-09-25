// SparseCholesky.cxx -- see SparseCholesky.hxx. Geometric nested dissection,
// a multifrontal factorisation over the separator tree, and the two triangular
// sweeps over the same tree, each parallel over independent subtrees.
#include "SparseCholesky.hxx"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <utility>

#include "Parallel.hxx"

namespace shapedna {

namespace {

using Clock = std::chrono::steady_clock;

double since(Clock::time_point t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

// Below these a subtree is not worth a task of its own: node sets for the
// dissection, flops for the factorisation, entries of L for a solve. (Unused
// without OpenMP, where the task pragmas they feed compile away.)
[[maybe_unused]] constexpr int kDissectTask = 4096;
[[maybe_unused]] constexpr int kSymbolicTask = 4096;
[[maybe_unused]] constexpr double kFactorTask = 2.0e5;
[[maybe_unused]] constexpr long long kSolveTask = 20000;

// ---------------------------------------------------------------------------
// Nested dissection
// ---------------------------------------------------------------------------
//
// dissect() numbers a node set into [offset, offset + size) -- the two halves
// first, the separator last -- and returns the tree nodes it made for it.
//
// Membership is read off `label`, one tag per half per call, drawn from an
// atomic counter so no two calls ever share one. That is what makes the halves
// safe to dissect concurrently without copying the graph: a node's neighbours
// are either in its own set, whose labels only its own task writes, or in the
// separator of an ancestor, whose labels were written before that task was
// spawned and are never written again.
struct Dissection {
    struct Raw {
        int first, size;
        std::vector<int> children;   // Raw ids
    };

    const SparseMatrix &A;
    const std::vector<std::array<double, 2>> &xy;
    const int leafSize;
    std::vector<int> &perm;
    std::vector<int> label;
    std::atomic<int> nextTag{1};
    std::vector<Raw> raw;

    Dissection(const SparseMatrix &A_, const std::vector<std::array<double, 2>> &xy_,
               int leaf, std::vector<int> &perm_)
        : A(A_), xy(xy_), leafSize(std::max(1, leaf)), perm(perm_), label(A_.n, 0) {}

    int newNode(int first, int size, std::vector<int> children) {
        int id = 0;
        SDNA_OMP(critical(shapedna_dissection))
        {
            id = static_cast<int>(raw.size());
            raw.push_back(Raw{first, size, std::move(children)});
        }
        return id;
    }

    std::vector<int> dissect(std::vector<int> verts, int offset) {
        std::vector<int> roots;
        const int ns = static_cast<int>(verts.size());
        if (ns == 0) return roots;
        if (ns <= leafSize) {
            std::copy(verts.begin(), verts.end(), perm.begin() + offset);
            roots.push_back(newNode(offset, ns, {}));
            return roots;
        }

        // Cut across the longer side of the bounding box, at the median.
        const double inf = std::numeric_limits<double>::infinity();
        double lo[2] = {inf, inf}, hi[2] = {-inf, -inf};
        for (int v : verts)
            for (int d = 0; d < 2; ++d) {
                lo[d] = std::min(lo[d], xy[v][d]);
                hi[d] = std::max(hi[d], xy[v][d]);
            }
        const int axis = (hi[0] - lo[0] >= hi[1] - lo[1]) ? 0 : 1;
        const int half = ns / 2;
        // Ties broken by index: a structured grid puts whole rows of nodes on
        // one coordinate, and the split must still be a fixed function of the set.
        std::nth_element(verts.begin(), verts.begin() + half, verts.end(), [&](int a, int b) {
            const double pa = xy[a][axis], pb = xy[b][axis];
            return pa < pb || (pa == pb && a < b);
        });

        const int tag = nextTag.fetch_add(1);
        const int tagL = 2 * tag, tagR = 2 * tag + 1;
        for (int i = 0; i < ns; ++i) label[verts[i]] = (i < half) ? tagL : tagR;

        // The nodes of each half with a neighbour in the other. Either set
        // separates the halves; take the smaller.
        std::vector<int> sepL, sepR;
        for (int i = 0; i < ns; ++i) {
            const int v = verts[i];
            const int other = (i < half) ? tagR : tagL;
            for (int k = A.rowPtr[v]; k < A.rowPtr[v + 1]; ++k)
                if (label[A.col[k]] == other) {
                    (i < half ? sepL : sepR).push_back(v);
                    break;
                }
        }
        std::vector<int> &sep = (sepL.size() <= sepR.size()) ? sepL : sepR;
        const int sepTag = 2 * nextTag.fetch_add(1);   // carried by no set, ever
        for (int v : sep) label[v] = sepTag;

        std::vector<int> left, right;
        left.reserve(half);
        right.reserve(ns - half);
        for (int v : verts) {
            if (label[v] == tagL) left.push_back(v);
            else if (label[v] == tagR) right.push_back(v);
        }
        const int nl = static_cast<int>(left.size());
        const int nr = static_cast<int>(right.size());
        const int nsep = static_cast<int>(sep.size());

        if (nl == 0 && nr == 0) {
            // Every node touches the other side: a dense little clique. Nothing
            // to gain by cutting it; eliminate it as one block.
            std::copy(verts.begin(), verts.end(), perm.begin() + offset);
            roots.push_back(newNode(offset, ns, {}));
            return roots;
        }

        std::vector<int> rootsL, rootsR;
        [[maybe_unused]] const bool spawn = ns > kDissectTask;
        SDNA_OMP(task shared(left, rootsL) if(spawn))
        rootsL = dissect(std::move(left), offset);
        SDNA_OMP(task shared(right, rootsR) if(spawn))
        rootsR = dissect(std::move(right), offset + nl);
        SDNA_OMP(taskwait)

        std::copy(sep.begin(), sep.end(), perm.begin() + offset + nl + nr);
        std::vector<int> children = std::move(rootsL);
        children.insert(children.end(), rootsR.begin(), rootsR.end());
        // No separator: the halves never touch, and they stay separate trees.
        if (nsep == 0) return children;
        roots.push_back(newNode(offset + nl + nr, nsep, std::move(children)));
        return roots;
    }
};

// ---------------------------------------------------------------------------
// Dense partial Cholesky of one front
// ---------------------------------------------------------------------------
//
// F is m x m, column-major, its lower triangle assembled. Eliminates the first
// `ns` columns: on return F(:, 0:ns) holds those columns of L and the trailing
// (m - ns) x (m - ns) lower triangle holds the Schur complement -- the update
// matrix for the parent.
//
// Right-looking over panels of kPanel columns: factor the panel's diagonal
// block, solve the rows below it against that block, then subtract the panel's
// outer product from every later column. The last two steps are the O(m^2)
// and O(m^3) ones and each splits over independent row blocks and columns, so
// they are what runs as a taskloop on a large front. Every entry is updated
// panel after panel and, within a panel, four columns at a time in a fixed
// grouping, so the result is the same however the loops are split.
constexpr int kPanel = 32;
constexpr int kRowBlock = 64;

bool partialCholesky(double *F, int m, int ns, bool parallel) {
    const std::size_t M = static_cast<std::size_t>(m);
    for (int k0 = 0; k0 < ns; k0 += kPanel) {
        const int k1 = std::min(ns, k0 + kPanel);

        // (1) The panel's diagonal block, unblocked.
        for (int k = k0; k < k1; ++k) {
            double *ck = F + k * M;
            if (!(ck[k] > 0.0)) return false;
            const double r = std::sqrt(ck[k]);
            ck[k] = r;
            const double inv = 1.0 / r;
            for (int i = k + 1; i < k1; ++i) ck[i] *= inv;
            for (int j = k + 1; j < k1; ++j) {
                const double a = ck[j];
                double *cj = F + j * M;
                for (int i = j; i < k1; ++i) cj[i] -= ck[i] * a;
            }
        }
        if (k1 == m) break;

        // (2) The rows below the diagonal block: L21 = F21 L11^{-T}, a block of
        // rows at a time.
        const int nblk = (m - k1 + kRowBlock - 1) / kRowBlock;
        auto panelRows = [=](int b) {
            const int r0 = k1 + b * kRowBlock;
            const int r1 = std::min(m, r0 + kRowBlock);
            for (int k = k0; k < k1; ++k) {
                double *ck = F + k * M;
                for (int j = k0; j < k; ++j) {
                    const double a = F[k + j * M];
                    const double *cj = F + j * M;
                    for (int i = r0; i < r1; ++i) ck[i] -= cj[i] * a;
                }
                const double inv = 1.0 / ck[k];
                for (int i = r0; i < r1; ++i) ck[i] *= inv;
            }
        };

        // (3) Every later column j loses the panel's contribution,
        // F(j:m, j) -= L(j:m, k0:k1) L(j, k0:k1)^T.
        auto trailing = [=](int j) {
            double *cj = F + j * M;
            int k = k0;
            for (; k + 3 < k1; k += 4) {
                const double *c0 = F + k * M, *c1 = c0 + M, *c2 = c1 + M, *c3 = c2 + M;
                const double a0 = c0[j], a1 = c1[j], a2 = c2[j], a3 = c3[j];
                for (int i = j; i < m; ++i)
                    cj[i] -= (c0[i] * a0 + c1[i] * a1) + (c2[i] * a2 + c3[i] * a3);
            }
            for (; k < k1; ++k) {
                const double *ck = F + k * M;
                const double a = ck[j];
                for (int i = j; i < m; ++i) cj[i] -= ck[i] * a;
            }
        };

        if (parallel) {
            SDNA_OMP(taskloop grainsize(1))
            for (int b = 0; b < nblk; ++b) panelRows(b);
            SDNA_OMP(taskloop grainsize(4))
            for (int j = k1; j < m; ++j) trailing(j);
        } else {
            for (int b = 0; b < nblk; ++b) panelRows(b);
            for (int j = k1; j < m; ++j) trailing(j);
        }
    }
    return true;
}

}  // namespace

// ---------------------------------------------------------------------------
// Ordering and analysis
// ---------------------------------------------------------------------------

void SparseCholesky::order(const SparseMatrix &A,
                           const std::vector<std::array<double, 2>> &coords) {
    perm.assign(n, -1);
    Dissection nd(A, coords, options.leafSize, perm);
    std::vector<int> all(n);
    std::iota(all.begin(), all.end(), 0);
    std::vector<int> top;
    SDNA_OMP(parallel)
    SDNA_OMP(single)
    top = nd.dissect(std::move(all), 0);

    // The dissection made its tree nodes in whatever order its tasks ran.
    // Renumber them by where their pivots sit, which is fixed: a post-order,
    // children before parents.
    const int count = static_cast<int>(nd.raw.size());
    std::vector<int> byFirst(count);
    std::iota(byFirst.begin(), byFirst.end(), 0);
    std::sort(byFirst.begin(), byFirst.end(),
              [&](int a, int b) { return nd.raw[a].first < nd.raw[b].first; });
    std::vector<int> id(count);
    for (int k = 0; k < count; ++k) id[byFirst[k]] = k;

    nodes.assign(count, Supernode());
    for (int k = 0; k < count; ++k) {
        const Dissection::Raw &r = nd.raw[byFirst[k]];
        Supernode &S = nodes[k];
        S.first = r.first;
        S.size = r.size;
        for (int c : r.children) S.children.push_back(id[c]);
        std::sort(S.children.begin(), S.children.end());
    }
    for (int k = 0; k < count; ++k) {
        nodes[k].subtreeFirst = nodes[k].first;
        for (int c : nodes[k].children) {
            nodes[c].parent = k;
            nodes[k].subtreeFirst = std::min(nodes[k].subtreeFirst, nodes[c].subtreeFirst);
        }
    }
    roots.clear();
    for (int r : top) roots.push_back(id[r]);
    std::sort(roots.begin(), roots.end());

    iperm.assign(n, -1);
    for (int i = 0; i < n; ++i) iperm[perm[i]] = i;
}

void SparseCholesky::symbolicSubtree(int s) {
    Supernode &S = nodes[s];
    for (int c : S.children) {
        [[maybe_unused]] const int span = nodes[c].first + nodes[c].size - nodes[c].subtreeFirst;
        SDNA_OMP(task if(span > kSymbolicTask))
        symbolicSubtree(c);
    }
    SDNA_OMP(taskwait)

    // The rows of L this supernode's columns reach beyond its own pivots: its
    // own entries of A below the pivots, and whatever its children could not
    // eliminate. Everything in there belongs to an ancestor.
    const int last = S.first + S.size;
    std::vector<int> rows;
    for (int j = S.first; j < last; ++j)
        for (int k = lowPtr[j]; k < lowPtr[j + 1]; ++k)
            if (lowRow[k] >= last) rows.push_back(lowRow[k]);
    for (int c : S.children)
        for (int i : nodes[c].rows)
            if (i >= last) rows.push_back(i);
    std::sort(rows.begin(), rows.end());
    rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
    S.rows = std::move(rows);

    // Where each child's update rows land in this front.
    for (int c : S.children) {
        Supernode &C = nodes[c];
        C.relative.resize(C.rows.size());
        for (std::size_t k = 0; k < C.rows.size(); ++k) {
            const int i = C.rows[k];
            C.relative[k] = (i < last)
                                ? i - S.first
                                : S.size + static_cast<int>(std::lower_bound(S.rows.begin(), S.rows.end(), i) -
                                                            S.rows.begin());
        }
    }

    // Eliminating pivot k of an m-row front updates the (m-k-1)(m-k)/2 entries
    // below and right of it, two flops each.
    const double m = static_cast<double>(S.size + S.rows.size());
    S.work = 0.0;
    for (int k = 0; k < S.size; ++k) S.work += (m - k) * (m - k);
    S.subtreeWork = S.work;
    S.subtreeEntries = static_cast<long long>(m) * S.size;
    for (int c : S.children) {
        S.subtreeWork += nodes[c].subtreeWork;
        S.subtreeEntries += nodes[c].subtreeEntries;
    }
}

void SparseCholesky::symbolic(const SparseMatrix &A) {
    // The permuted matrix's lower triangle, column by column: new column j is
    // old row perm[j], each entry renumbered and kept if it is on or below the
    // diagonal. Counted and filled in parallel, each column in its own slot.
    lowPtr.assign(n + 1, 0);
    SDNA_OMP(parallel for schedule(dynamic, 1024) if(n > 4096))
    for (int j = 0; j < n; ++j) {
        const int v = perm[j];
        int c = 0;
        for (int k = A.rowPtr[v]; k < A.rowPtr[v + 1]; ++k)
            if (iperm[A.col[k]] >= j) ++c;
        lowPtr[j + 1] = c;
    }
    for (int j = 0; j < n; ++j) lowPtr[j + 1] += lowPtr[j];
    lowRow.resize(lowPtr[n]);
    lowVal.resize(lowPtr[n]);
    SDNA_OMP(parallel for schedule(dynamic, 1024) if(n > 4096))
    for (int j = 0; j < n; ++j) {
        const int v = perm[j];
        int p = lowPtr[j];
        for (int k = A.rowPtr[v]; k < A.rowPtr[v + 1]; ++k) {
            const int i = iperm[A.col[k]];
            if (i >= j) {
                lowRow[p] = i;
                lowVal[p] = A.val[k];
                ++p;
            }
        }
    }
    stats.nnzA = lowPtr[n];

    SDNA_OMP(parallel if(nodes.size() > 1))
    SDNA_OMP(single)
    {
        for (int r : roots) {
            SDNA_OMP(task)
            symbolicSubtree(r);
        }
        SDNA_OMP(taskwait)
    }

    stats.supernodes = static_cast<int>(nodes.size());
    stats.nnzL = 0;
    stats.flops = 0.0;
    stats.maxFront = 0;
    for (const Supernode &S : nodes) {
        stats.nnzL += static_cast<long long>(S.size + S.rows.size()) * S.size;
        stats.flops += S.work;
        stats.maxFront = std::max(stats.maxFront, S.size + static_cast<int>(S.rows.size()));
    }
    // Depth, from the roots down; parents come after their children.
    std::vector<int> depth(nodes.size(), 1);
    stats.treeDepth = nodes.empty() ? 0 : 1;
    for (int s = static_cast<int>(nodes.size()) - 1; s >= 0; --s) {
        if (nodes[s].parent >= 0) depth[s] = depth[nodes[s].parent] + 1;
        stats.treeDepth = std::max(stats.treeDepth, depth[s]);
    }
}

// ---------------------------------------------------------------------------
// Numeric factorisation
// ---------------------------------------------------------------------------

bool SparseCholesky::factorFront(int s, std::vector<std::vector<double>> *updates) {
    Supernode &S = nodes[s];
    const int ns = S.size;
    const int nu = static_cast<int>(S.rows.size());
    const int m = ns + nu;
    const std::size_t M = static_cast<std::size_t>(m);
    const int last = S.first + ns;

    std::vector<double> F(M * M, 0.0);

    // This supernode's own columns of A.
    for (int jj = 0; jj < ns; ++jj) {
        const int j = S.first + jj;
        double *cj = F.data() + jj * M;
        for (int k = lowPtr[j]; k < lowPtr[j + 1]; ++k) {
            const int i = lowRow[k];
            const int p = (i < last)
                              ? i - S.first
                              : ns + static_cast<int>(std::lower_bound(S.rows.begin(), S.rows.end(), i) -
                                                      S.rows.begin());
            cj[p] += lowVal[k];
        }
    }

    // Extend-add: the children's update matrices, in child order.
    for (int c : S.children) {
        const Supernode &C = nodes[c];
        const int nc = static_cast<int>(C.rows.size());
        const std::vector<double> &U = (*updates)[c];
        for (int jc = 0; jc < nc; ++jc) {
            double *cj = F.data() + static_cast<std::size_t>(C.relative[jc]) * M;
            const double *uj = U.data() + static_cast<std::size_t>(jc) * nc;
            for (int ic = jc; ic < nc; ++ic) cj[C.relative[ic]] += uj[ic];
        }
        std::vector<double>().swap((*updates)[c]);
    }

    if (!partialCholesky(F.data(), m, ns, m >= options.parallelFront)) return false;

    S.L.assign(F.begin(), F.begin() + M * ns);
    std::vector<double> &U = (*updates)[s];
    U.assign(static_cast<std::size_t>(nu) * nu, 0.0);
    for (int j = 0; j < nu; ++j) {
        const double *src = F.data() + (ns + j) * M + ns;
        double *dst = U.data() + static_cast<std::size_t>(j) * nu;
        for (int i = j; i < nu; ++i) dst[i] = src[i];
    }
    return true;
}

void SparseCholesky::factorSubtree(int s, std::vector<std::vector<double>> *updates,
                                   std::atomic<bool> *bad) {
    for (int c : nodes[s].children) {
        SDNA_OMP(task if(nodes[c].subtreeWork > kFactorTask))
        factorSubtree(c, updates, bad);
    }
    SDNA_OMP(taskwait)
    if (bad->load()) return;
    if (!factorFront(s, updates)) bad->store(true);
}

bool SparseCholesky::factor(const SparseMatrix &A,
                            const std::vector<std::array<double, 2>> &coords) {
    stats = Stats();
    n = A.n;
    stats.n = n;
    nodes.clear();
    roots.clear();
    if (n == 0) return true;

    Clock::time_point t0 = Clock::now();
    order(A, coords);
    stats.orderSeconds = since(t0);

    t0 = Clock::now();
    symbolic(A);
    stats.symbolicSeconds = since(t0);

    t0 = Clock::now();
    std::vector<std::vector<double>> updates(nodes.size());
    std::atomic<bool> bad{false};
    std::vector<std::vector<double>> *up = &updates;
    std::atomic<bool> *flag = &bad;
    SDNA_OMP(parallel)
    SDNA_OMP(single)
    {
        for (int r : roots) {
            SDNA_OMP(task if(nodes[r].subtreeWork > kFactorTask))
            factorSubtree(r, up, flag);
        }
        SDNA_OMP(taskwait)
    }
    stats.numericSeconds = since(t0);

    // The permuted copy of A was only needed to assemble the fronts.
    std::vector<int>().swap(lowPtr);
    std::vector<int>().swap(lowRow);
    std::vector<double>().swap(lowVal);
    return !bad.load();
}

// ---------------------------------------------------------------------------
// Solve
// ---------------------------------------------------------------------------

void SparseCholesky::forwardSubtree(int s, double *Y, int nrhs,
                                    std::vector<std::vector<double>> *acc) const {
    const Supernode &S = nodes[s];
    for (int c : S.children) {
        SDNA_OMP(task if(nodes[c].subtreeEntries > kSolveTask))
        forwardSubtree(c, Y, nrhs, acc);
    }
    SDNA_OMP(taskwait)

    const int ns = S.size;
    const int nu = static_cast<int>(S.rows.size());
    const std::size_t M = static_cast<std::size_t>(ns + nu);
    const std::size_t N = static_cast<std::size_t>(n);

    // What the children could not apply themselves: to this supernode's own
    // entries directly, and to the ancestors' by way of this node's buffer.
    std::vector<double> &a = (*acc)[s];
    a.assign(static_cast<std::size_t>(nu) * nrhs, 0.0);
    for (int c : S.children) {
        const Supernode &C = nodes[c];
        const int nc = static_cast<int>(C.rows.size());
        const std::vector<double> &ac = (*acc)[c];
        for (int r = 0; r < nrhs; ++r) {
            double *y = Y + S.first + r * N;
            double *ar = a.data() + static_cast<std::size_t>(r) * nu;
            const double *cr = ac.data() + static_cast<std::size_t>(r) * nc;
            for (int k = 0; k < nc; ++k) {
                const int p = C.relative[k];
                if (p < ns) y[p] += cr[k];
                else ar[p - ns] += cr[k];
            }
        }
        std::vector<double>().swap((*acc)[c]);
    }

    // Column k of L is read from memory once and applied to every right-hand
    // side while it is still in cache. Each entry still loses its terms in
    // column order, whatever nrhs is.
    for (int k = 0; k < ns; ++k) {
        const double *lk = S.L.data() + k * M;
        const double *lu = lk + ns;
        for (int r = 0; r < nrhs; ++r) {
            double *y = Y + S.first + r * N;
            double *ar = a.data() + static_cast<std::size_t>(r) * nu;
            y[k] /= lk[k];
            const double yk = y[k];
            for (int i = k + 1; i < ns; ++i) y[i] -= lk[i] * yk;
            for (int i = 0; i < nu; ++i) ar[i] -= lu[i] * yk;
        }
    }
}

void SparseCholesky::backwardSubtree(int s, double *Y, int nrhs) const {
    const Supernode &S = nodes[s];
    const int ns = S.size;
    const int nu = static_cast<int>(S.rows.size());
    const std::size_t M = static_cast<std::size_t>(ns + nu);
    const std::size_t N = static_cast<std::size_t>(n);

    // The ancestors' entries are final by now; this supernode's depend on them
    // and on nothing below it. They are gathered once, and each column of L
    // is then read once for all right-hand sides.
    std::vector<double> g(static_cast<std::size_t>(nu) * nrhs);
    for (int r = 0; r < nrhs; ++r)
        for (int i = 0; i < nu; ++i) g[i + static_cast<std::size_t>(r) * nu] = Y[S.rows[i] + r * N];
    for (int k = ns - 1; k >= 0; --k) {
        const double *lk = S.L.data() + k * M;
        for (int r = 0; r < nrhs; ++r) {
            double *y = Y + S.first + r * N;
            double t = y[k] - dotSerial(nu, lk + ns, g.data() + static_cast<std::size_t>(r) * nu);
            for (int i = k + 1; i < ns; ++i) t -= lk[i] * y[i];
            y[k] = t / lk[k];
        }
    }

    for (int c : S.children) {
        SDNA_OMP(task if(nodes[c].subtreeEntries > kSolveTask))
        backwardSubtree(c, Y, nrhs);
    }
    SDNA_OMP(taskwait)
}

void SparseCholesky::solve(double *X, int nrhs, int ld) const {
    if (n == 0 || nrhs <= 0) return;
    const std::size_t N = static_cast<std::size_t>(n);
    const std::size_t LD = static_cast<std::size_t>(ld);
    std::vector<double> Y(N * nrhs);
    double *y = Y.data();

    SDNA_OMP(parallel for schedule(dynamic, 1024) if(n > 4096))
    for (int i = 0; i < n; ++i)
        for (int r = 0; r < nrhs; ++r) y[i + r * N] = X[perm[i] + r * LD];

    std::vector<std::vector<double>> acc(nodes.size());
    std::vector<std::vector<double>> *pacc = &acc;
    SDNA_OMP(parallel if(nodes.size() > 1))
    SDNA_OMP(single)
    {
        for (int r : roots) {
            SDNA_OMP(task if(nodes[r].subtreeEntries > kSolveTask))
            forwardSubtree(r, y, nrhs, pacc);
        }
        SDNA_OMP(taskwait)
        for (int r : roots) {
            SDNA_OMP(task if(nodes[r].subtreeEntries > kSolveTask))
            backwardSubtree(r, y, nrhs);
        }
        SDNA_OMP(taskwait)
    }

    SDNA_OMP(parallel for schedule(dynamic, 1024) if(n > 4096))
    for (int i = 0; i < n; ++i)
        for (int r = 0; r < nrhs; ++r) X[perm[i] + r * LD] = y[i + r * N];
}

}  // namespace shapedna
