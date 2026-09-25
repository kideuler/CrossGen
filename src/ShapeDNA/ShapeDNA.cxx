// ShapeDNA.cxx -- see ShapeDNA.hxx.
#include "ShapeDNA.hxx"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <sstream>

#include "DenseEigen.hxx"
#include "LagrangeFEM.hxx"
#include "Lanczos.hxx"
#include "Parallel.hxx"
#include "SparseCholesky.hxx"

namespace shapedna {

namespace {

using Clock = std::chrono::steady_clock;

double since(Clock::time_point t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

struct Geometry {
    double area = 0.0;
    double boundaryLength = 0.0;
    double cornerTerm = 0.0;
    int components = 0;
    int euler = 0;
};

// Sum of f(i) for i in [0, n), over fixed chunks added in chunk order.
template <class F>
double chunkedSum(int n, F f) {
    const int chunk = static_cast<int>(kRowChunk);
    const int chunks = (n + chunk - 1) / chunk;
    std::vector<double> part(chunks, 0.0);
    SDNA_OMP(parallel for schedule(dynamic, 1) if(chunks > 4))
    for (int c = 0; c < chunks; ++c) {
        double s = 0.0;
        for (int i = c * chunk; i < std::min(n, (c + 1) * chunk); ++i) s += f(i);
        part[c] = s;
    }
    double s = 0.0;
    for (double x : part) s += x;
    return s;
}

// The quantities of the meshed domain the heat trace should give back, and the
// number of pieces (each one a zero mode under Neumann).
Geometry measure(const Mesh &mesh) {
    Geometry g;
    const int nT = static_cast<int>(mesh.triangles.size());
    const int nV = static_cast<int>(mesh.vertices.size());

    g.area = chunkedSum(nT, [&](int t) {
        const Triangle &tri = mesh.triangles[t];
        const Point &x0 = mesh.vertices[tri[0]];
        return 0.5 * std::fabs(cross2(mesh.vertices[tri[1]] - x0, mesh.vertices[tri[2]] - x0));
    });
    const int nB = static_cast<int>(mesh.boundaryEdges.size());
    g.boundaryLength = chunkedSum(nB, [&](int k) {
        const auto &e = mesh.edges[mesh.boundaryEdges[k]];
        return normP(mesh.vertices[e[1]] - mesh.vertices[e[0]]);
    });

    // The interior angle at a boundary vertex is the sum of its triangles'
    // angles there; a straight stretch of dS contributes pi, and nothing.
    const int nBV = static_cast<int>(mesh.boundaryVertices.size());
    g.cornerTerm = chunkedSum(nBV, [&](int k) {
        const int v = mesh.boundaryVertices[k];
        double alpha = 0.0;
        for (int s = mesh.vertexTriangles.rowPtr[v]; s < mesh.vertexTriangles.rowPtr[v + 1]; ++s) {
            const Triangle &tri = mesh.triangles[mesh.vertexTriangles.colIdx[s]];
            const int a = (tri[0] == v) ? 0 : (tri[1] == v) ? 1 : 2;
            const Point u = mesh.vertices[tri[(a + 1) % 3]] - mesh.vertices[v];
            const Point w = mesh.vertices[tri[(a + 2) % 3]] - mesh.vertices[v];
            alpha += std::atan2(std::fabs(cross2(u, w)), dotP(u, w));
        }
        if (!(alpha > 0.0)) return 0.0;
        return (M_PI * M_PI - alpha * alpha) / (24.0 * M_PI * alpha);
    });

    // Pieces: union-find over the triangle adjacency.
    std::vector<int> parent(nT);
    std::iota(parent.begin(), parent.end(), 0);
    auto find = [&](int x) {
        while (parent[x] != x) x = parent[x] = parent[parent[x]];
        return x;
    };
    for (int t = 0; t < nT; ++t)
        for (int k = 0; k < 3; ++k) {
            const int o = mesh.triangleAdjacency[t][k];
            if (o < 0) continue;
            const int a = find(t), b = find(o);
            if (a != b) parent[std::max(a, b)] = std::min(a, b);
        }
    for (int t = 0; t < nT; ++t)
        if (find(t) == t) ++g.components;

    std::vector<char> used(nV, 0);
    for (const Triangle &tri : mesh.triangles)
        for (int v : tri) used[v] = 1;
    const int usedV = static_cast<int>(std::count(used.begin(), used.end(), 1));
    g.euler = usedV - static_cast<int>(mesh.edges.size()) + nT;
    return g;
}

}  // namespace

ShapeDNA::ShapeDNA(const Mesh &m) : mesh(m) {}
ShapeDNA::ShapeDNA(const Mesh &m, const Options &opts) : mesh(m), options(opts) {}

bool ShapeDNA::run() {
    const Clock::time_point t0 = Clock::now();
    report = Report();
    lambda.clear();
    shapeDNA.clear();
    functions.clear();

    ScopedThreads scoped(options.threads);
    report.threads = maxThreads();
    report.openMP = haveOpenMP();

    auto fail = [&](const std::string &why) {
        report.messages.push_back(why);
        report.seconds = since(t0);
        return false;
    };
    if (options.degree < 1 || options.degree > 3) return fail("Options::degree must be 1, 2 or 3");
    if (options.count < 1) return fail("Options::count must be positive");
    if (options.refine < 1) return fail("Options::refine must be at least 1");
    if (mesh.triangles.empty()) return fail("the mesh has no triangles");

    // -- Sec. 2: refine ------------------------------------------------------
    Clock::time_point ts = Clock::now();
    Mesh refined;
    const Mesh *domain = &mesh;
    if (options.refine > 1) {
        refined = refineUniform(mesh, options.refine);
        domain = &refined;
    }
    const Geometry geo = measure(*domain);
    report.area = geo.area;
    report.boundaryLength = geo.boundaryLength;
    report.cornerTerm = geo.cornerTerm;
    report.components = geo.components;
    report.eulerCharacteristic = geo.euler;
    if (!(geo.area > 0.0)) return fail("the mesh has zero area");
    const double weyl = 4.0 * M_PI / geo.area;   // lambda_k ~ weyl * k

    // -- Secs. 3-5: the two matrices on the free nodes ----------------------
    const bool dirichlet = options.boundary == Dirichlet;
    const LagrangeSpace space(*domain, options.degree);
    const ReferenceElement ref(options.degree);
    FEMatrices fe = assemble(space, ref, dirichlet);
    const int nd = fe.stiffness.n;
    report.triangles = space.numElements();
    report.nodes = space.numNodes();
    report.unknowns = nd;
    report.nnz = fe.stiffness.nnz();
    report.degenerateElements = fe.degenerateElements;
    report.secondsAssemble = since(ts);
    if (fe.degenerateElements > 0) {
        std::ostringstream os;
        os << fe.degenerateElements << " zero-area triangle(s) left out of the assembly";
        report.messages.push_back(os.str());
    }
    if (nd == 0) return fail("no unknowns: every node is on the boundary");

    const int zeros = dirichlet ? 0 : geo.components;
    int nev = options.count + zeros;
    if (nev > nd) {
        std::ostringstream os;
        os << "the discretisation has only " << nd << " eigenvalues; asked for " << nev;
        report.messages.push_back(os.str());
        nev = nd;
    }
    report.zeroModes = zeros;

    double sigma = options.shift;
    if (std::isnan(sigma)) sigma = dirichlet ? 0.0 : -0.01 * weyl;
    report.shift = sigma;

    // -- Secs. 6-7: factor K = A - sigma B, then Lanczos ---------------------
    std::vector<double> values, vectors;
    bool dense = nd <= options.denseBelow;
    if (!dense) {
        ts = Clock::now();
        SparseMatrix K = fe.stiffness;
        if (sigma != 0.0) {
            const long long nz = K.nnz();
            SDNA_OMP(parallel for schedule(static) if(nz > 100000))
            for (long long k = 0; k < nz; ++k) K.val[k] -= sigma * fe.mass.val[k];
        }
        SparseCholesky::Options co;
        co.leafSize = options.leafSize;
        SparseCholesky chol(co);
        const bool factored = chol.factor(K, fe.dofPositions);
        const SparseCholesky::Stats &cs = chol.getStats();
        report.supernodes = cs.supernodes;
        report.maxFront = cs.maxFront;
        report.treeDepth = cs.treeDepth;
        report.nnzL = cs.nnzL;
        report.factorFlops = cs.flops;
        report.secondsOrder = cs.orderSeconds;
        report.secondsFactor = since(ts) - cs.orderSeconds;
        if (!factored) {
            std::ostringstream os;
            os << "A - sigma B is not positive definite at sigma = " << sigma
               << " (a pivot came out non-positive); the shift must sit below lambda_1";
            return fail(os.str());
        }

        ts = Clock::now();
        LanczosOptions lo;
        lo.nev = nev;
        lo.blockSize = options.blockSize;
        lo.basisSize = options.basisSize;
        lo.tolerance = options.tolerance;
        lo.maxRestarts = options.maxRestarts;
        lo.seed = options.seed;
        lo.vectors = options.eigenfunctions;
        LanczosResult lr;
        const bool ok = shiftInvertLanczos(fe.mass, chol, sigma, lo, lr);
        report.secondsEigen = since(ts);
        report.basisSize = lr.basisSize;
        report.blockSize = lr.blockSize;
        if (!ok && lr.values.empty()) {
            // Too few unknowns for the basis: the dense solve is the answer.
            report.messages.push_back(lr.message);
            dense = true;
        } else {
            report.converged = ok;
            report.restarts = lr.restarts;
            report.blockSteps = lr.blockSteps;
            report.solves = lr.solves;
            report.deflations = lr.deflations;
            for (double r : lr.residuals) report.maxRitzResidual = std::max(report.maxRitzResidual, r);
            if (!ok) report.messages.push_back(lr.message);
            values = std::move(lr.values);
            vectors = std::move(lr.vectors);
        }
    }
    if (dense) {
        ts = Clock::now();
        const std::size_t N = static_cast<std::size_t>(nd);
        std::vector<double> Ad(N * N, 0.0), Bd(N * N, 0.0);
        for (int i = 0; i < nd; ++i)
            for (int k = fe.stiffness.rowPtr[i]; k < fe.stiffness.rowPtr[i + 1]; ++k) {
                Ad[i + fe.stiffness.col[k] * N] = fe.stiffness.val[k];
                Bd[i + fe.mass.col[k] * N] = fe.mass.val[k];
            }
        std::vector<double> w, Y;
        const bool ok = generalizedSymmetricEigen(nd, Ad.data(), Bd.data(), w,
                                                  options.eigenfunctions ? &Y : nullptr);
        report.secondsEigen += since(ts);
        if (!ok) return fail("the dense generalized eigenproblem failed (mass matrix not positive definite?)");
        report.dense = true;
        report.converged = true;
        values.assign(w.begin(), w.begin() + nev);
        if (options.eigenfunctions) vectors.assign(Y.begin(), Y.begin() + N * nev);
    }

    lambda = values;
    report.eigenvalues = static_cast<int>(lambda.size());

    // -- Sec. 9: drop the zero modes, normalise -----------------------------
    for (int i = 0; i < std::min(zeros, report.eigenvalues); ++i)
        if (std::fabs(lambda[i]) > 1e-6 * weyl) {
            std::ostringstream os;
            os << "Neumann eigenvalue " << i << " should be zero (one per piece of the mesh) and is "
               << lambda[i];
            report.messages.push_back(os.str());
        }
    std::vector<double> raw(lambda.begin() + std::min(zeros, report.eigenvalues), lambda.end());
    double normalizer = 1.0;
    switch (options.normalization) {
        case NoNormalization: normalizer = 1.0; break;
        case AreaNormalization: normalizer = 1.0 / geo.area; break;
        case FirstEigenvalue: normalizer = raw.empty() ? 1.0 : raw[0]; break;
        case WeylSlope: {
            double num = 0.0, den = 0.0;
            for (std::size_t k = 0; k < raw.size(); ++k) {
                num += (k + 1.0) * raw[k];
                den += (k + 1.0) * (k + 1.0);
            }
            normalizer = den > 0.0 ? num / den : 1.0;
            break;
        }
    }
    report.normalizer = normalizer;
    shapeDNA.resize(raw.size());
    for (std::size_t k = 0; k < raw.size(); ++k) shapeDNA[k] = raw[k] / normalizer;

    // -- the eigenfunctions, and their residuals on the assembled matrices --
    if (options.eigenfunctions && !vectors.empty()) {
        const int nV = static_cast<int>(mesh.vertices.size());
        const int nk = report.eigenvalues;
        const std::size_t N = static_cast<std::size_t>(nd);
        functions.assign(static_cast<std::size_t>(nV) * nk, 0.0);
        std::vector<double> Ax(N), Bx(N);
        for (int k = 0; k < nk; ++k) {
            const double *x = vectors.data() + k * N;
            for (int v = 0; v < nV; ++v) {
                const int d = fe.nodeDof[v];
                if (d >= 0) functions[v + static_cast<std::size_t>(k) * nV] = x[d];
            }
            fe.stiffness.multiply(x, Ax.data());
            fe.mass.multiply(x, Bx.data());
            double r2 = 0.0, b2 = 0.0;
            for (std::size_t i = 0; i < N; ++i) {
                const double r = Ax[i] - lambda[k] * Bx[i];
                r2 += r * r;
                b2 += Bx[i] * Bx[i];
            }
            const double scale = std::max(std::fabs(lambda[k]), weyl) * std::sqrt(b2);
            if (scale > 0.0) report.maxTrueResidual = std::max(report.maxTrueResidual, std::sqrt(r2) / scale);
        }
    }

    report.ran = true;
    report.seconds = since(t0);
    return report.converged;
}

double ShapeDNA::distance(const std::vector<double> &a, const std::vector<double> &b) {
    const std::size_t n = std::min(a.size(), b.size());
    double s = 0.0;
    for (std::size_t i = 0; i < n; ++i) s += (a[i] - b[i]) * (a[i] - b[i]);
    return std::sqrt(s);
}

double ShapeDNA::heatTrace(double t) const {
    double z = 0.0;
    for (double l : lambda) z += std::exp(-l * t);
    return z;
}

ShapeDNA::HeatTraceFit ShapeDNA::fitHeatTrace() const {
    HeatTraceFit fit;
    if (lambda.size() < 10 || !(lambda.back() > 0.0)) return fit;
    fit.tMin = 16.0 / lambda.back();
    fit.tMax = 64.0 / lambda.back();

    // Least squares for y(t) = 4 pi t Z(t) = c0 + c1 sqrt(t) + c2 t, in the
    // normalised variable s = sqrt(t / tMax) so the 3 x 3 normal equations
    // are well scaled.
    const int samples = 64;
    double G[3][3] = {{0}}, r[3] = {0};
    std::vector<double> ss(samples), ys(samples);
    for (int k = 0; k < samples; ++k) {
        const double t = fit.tMin * std::pow(fit.tMax / fit.tMin, static_cast<double>(k) / (samples - 1));
        const double s = std::sqrt(t / fit.tMax);
        const double y = 4.0 * M_PI * t * heatTrace(t);
        const double phi[3] = {1.0, s, s * s};
        for (int i = 0; i < 3; ++i) {
            r[i] += phi[i] * y;
            for (int j = 0; j < 3; ++j) G[i][j] += phi[i] * phi[j];
        }
        ss[k] = s;
        ys[k] = y;
    }
    // Gaussian elimination on the 3 x 3 system (symmetric positive definite).
    for (int k = 0; k < 3; ++k)
        for (int i = k + 1; i < 3; ++i) {
            const double f = G[i][k] / G[k][k];
            for (int j = k; j < 3; ++j) G[i][j] -= f * G[k][j];
            r[i] -= f * r[k];
        }
    double c[3];
    for (int i = 2; i >= 0; --i) {
        double s = r[i];
        for (int j = i + 1; j < 3; ++j) s -= G[i][j] * c[j];
        c[i] = s / G[i][i];
    }
    // Back from s to t: c1 sqrt(t) = c1' s means c1 = c1' / sqrt(tMax).
    const double c0 = c[0], c1 = c[1] / std::sqrt(fit.tMax), c2 = c[2] / fit.tMax;
    double res2 = 0.0;
    for (int k = 0; k < samples; ++k) {
        const double e = ys[k] - (c[0] + c[1] * ss[k] + c[2] * ss[k] * ss[k]);
        res2 += e * e;
    }
    fit.area = c0;
    fit.boundaryLength = (options.boundary == Dirichlet ? -2.0 : 2.0) * c1 / std::sqrt(M_PI);
    fit.constant = c2 / (4.0 * M_PI);
    fit.rmsResidual = std::sqrt(res2 / samples) / std::max(std::fabs(c0), 1e-300);
    fit.valid = std::isfinite(c0) && std::isfinite(c1) && std::isfinite(c2);
    return fit;
}

}  // namespace shapedna
