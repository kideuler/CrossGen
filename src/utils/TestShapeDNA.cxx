// Utility to compute the Shape-DNA of planar meshes (src/ShapeDNA), and to
// self-test the module on domains whose spectra are known in closed form.
//
//   TestShapeDNA <mesh.obj> [more.obj ...] [options]
//   TestShapeDNA --selftest
//
// For each mesh it prints the discretisation, what the factorisation and the
// eigensolver did, and the Shape-DNA vector; given several meshes it also
// prints the matrix of pairwise distances between their vectors (Sec. 9 of
// docs/shape_dna.md).
//
//   --count n          length of the DNA (50)
//   --degree p         Lagrange degree, 1-3 (3)
//   --refine r         split every edge into r parts first (1)
//   --neumann          Neumann instead of Dirichlet boundary conditions
//   --norm area|first|weyl|weyl-ratio|none   scale normalisation (area)
//   --threads k        OpenMP threads (runtime default)
//   --block b          Lanczos block size (4)
//   --basis m          Lanczos basis size (2 nev + 2 b)
//   --tol x            Ritz residual tolerance (1e-10)
//   --leaf n           nested-dissection leaf size (64)
//   --dense-below n    solve densely at or below n unknowns (400)
//   --heat             fit the heat trace and compare with the mesh (Sec. 8)
//   --eigenfunctions   also compute eigenvectors and their true residuals
//   --raw              print the raw eigenvalues as well as the DNA
//   --csv file         write one row per mesh: name, then the DNA
//
// Exits non-zero when a mesh fails to converge. --selftest needs no files: it
// builds its meshes in memory.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "TriangleMesher.hpp"
#include "mesh/Mesh.hxx"
#include "ShapeDNA/DenseEigen.hxx"
#include "ShapeDNA/LagrangeFEM.hxx"
#include "ShapeDNA/Parallel.hxx"
#include "ShapeDNA/ShapeDNA.hxx"
#include "ShapeDNA/SparseCholesky.hxx"

using shapedna::ShapeDNA;

namespace {

const char *kPass = "\033[32m[PASS]\033[0m";
const char *kFail = "\033[31m[FAIL]\033[0m";
const char *kInfo = "\033[36m[INFO]\033[0m";

int gFailures = 0;

void verdict(bool ok, const std::string &what) {
    if (!ok) ++gFailures;
    std::cout << "  " << (ok ? kPass : kFail) << " " << what << "\n";
}

void info(const std::string &what) { std::cout << "  " << kInfo << " " << what << "\n"; }

void heading(const std::string &title) {
    std::cout << "\n" << title << "\n" << std::string(title.size(), '-') << "\n";
}

std::string sci(double x, int digits = 2) {
    std::ostringstream os;
    os << std::scientific << std::setprecision(digits) << x;
    return os.str();
}

std::string fix(double x, int digits = 4) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(digits) << x;
    return os.str();
}

double seconds(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// A deterministic generator, so a failure is reproducible without a seed file.
struct Rand {
    unsigned long long s;
    explicit Rand(unsigned long long seed) : s(seed ? seed : 88172645463325252ULL) {}
    double next() {   // uniform in [-1, 1)
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        return 2.0 * (static_cast<double>(s >> 11) / 9007199254740992.0) - 1.0;
    }
};

// ---------------------------------------------------------------------------
// meshes built in memory
// ---------------------------------------------------------------------------

// [0,w] x [0,h] in nx x ny cells, each cut along its (i,j)-(i+1,j+1) diagonal.
Mesh makeRectangle(double w, double h, int nx, int ny) {
    std::vector<Point> verts;
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) verts.push_back(Point{w * i / nx, h * j / ny});
    std::vector<Triangle> tris;
    auto id = [&](int i, int j) { return j * (nx + 1) + i; };
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            tris.push_back({id(i, j), id(i + 1, j), id(i + 1, j + 1)});
            tris.push_back({id(i, j), id(i + 1, j + 1), id(i, j + 1)});
        }
    return Mesh(verts, tris);
}

// The unit square in n x n cells, each cut into four by both diagonals about a
// node at its centre. Unlike makeRectangle's this triangulation has every
// symmetry of the square, 90-degree turns included, and those are what make
// the pairs (m, n), (n, m) one two-dimensional eigenspace: a mesh with only the
// diagonal reflection splits them.
Mesh makeCrissCross(int n) {
    std::vector<Point> verts;
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i) verts.push_back(Point{static_cast<double>(i) / n, static_cast<double>(j) / n});
    auto id = [&](int i, int j) { return j * (n + 1) + i; };
    std::vector<Triangle> tris;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            const int c = static_cast<int>(verts.size());
            verts.push_back(Point{(i + 0.5) / n, (j + 0.5) / n});
            tris.push_back({id(i, j), id(i + 1, j), c});
            tris.push_back({id(i + 1, j), id(i + 1, j + 1), c});
            tris.push_back({id(i + 1, j + 1), id(i, j + 1), c});
            tris.push_back({id(i, j + 1), id(i, j), c});
        }
    return Mesh(verts, tris);
}

void addCircle(double cx, double cy, double r, int n, triangle_wrapper::TriangleMesher2D::MeshInput &in,
               int type) {
    const int base = static_cast<int>(in.vertlist.size());
    std::vector<std::array<int, 2>> seg;
    for (int i = 0; i < n; ++i) {
        const double a = 2.0 * M_PI * i / n;
        in.vertlist.push_back({cx + r * std::cos(a), cy + r * std::sin(a)});
        seg.push_back({base + i, base + (i + 1) % n});
    }
    in.segment_loops.push_back(seg);
    in.type.push_back(type);
}

Mesh triangulate(triangle_wrapper::TriangleMesher2D::MeshInput &in) {
    triangle_wrapper::TriangleMesher2D mesher;
    const auto out = mesher.triangulate(in);
    return Mesh(out.verts, out.triangles);
}

// The unit-radius disk (or an annulus, with `hole` > 0), its circles cut into
// nb chords, meshed to edge length h.
Mesh makeDisk(double radius, int nb, double h, double hole = 0.0) {
    triangle_wrapper::TriangleMesher2D::MeshInput in;
    addCircle(0.0, 0.0, radius, nb, in, 0);
    if (hole > 0.0) {
        addCircle(0.0, 0.0, hole, std::max(16, static_cast<int>(nb * hole / radius)), in, 1);
        in.loop_seed = {{0.0, 0.0}, {0.0, 0.0}};   // the hole's seed: its centre
    }
    in.h = h;
    return triangulate(in);
}

// An L: the unit square with its top-right quarter removed.
Mesh makeL(double h) {
    triangle_wrapper::TriangleMesher2D::MeshInput in;
    in.vertlist = {{0, 0}, {1, 0}, {1, 0.5}, {0.5, 0.5}, {0.5, 1}, {0, 1}};
    in.segment_loops.push_back({{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}});
    in.type.push_back(0);
    in.h = h;
    return triangulate(in);
}

// The same mesh moved by a similarity: rotated by `angle`, scaled by `scale`,
// shifted by (tx, ty), and mirrored in the y axis first when asked. A mirror
// reverses every triangle, so each is turned back to counter-clockwise.
Mesh similar(const Mesh &m, double angle, double scale, double tx, double ty, bool mirror) {
    std::vector<Point> v = m.vertices;
    const double c = std::cos(angle), s = std::sin(angle);
    for (Point &p : v) {
        const double x = mirror ? -p[0] : p[0], y = p[1];
        p = Point{scale * (c * x - s * y) + tx, scale * (s * x + c * y) + ty};
    }
    std::vector<Triangle> t = m.triangles;
    if (mirror)
        for (Triangle &tri : t) std::swap(tri[1], tri[2]);
    return Mesh(v, t, m.triangleMatId);
}

// Two meshes side by side, not touching: one mesh of two pieces.
Mesh disjointUnion(const Mesh &a, const Mesh &b, double dx) {
    std::vector<Point> v = a.vertices;
    for (const Point &p : b.vertices) v.push_back(Point{p[0] + dx, p[1]});
    std::vector<Triangle> t = a.triangles;
    const int off = static_cast<int>(a.vertices.size());
    for (const Triangle &tri : b.triangles) t.push_back({tri[0] + off, tri[1] + off, tri[2] + off});
    return Mesh(v, t);
}

// ---------------------------------------------------------------------------
// spectra known in closed form
// ---------------------------------------------------------------------------

// pi^2 (m^2 / a^2 + n^2 / b^2), m, n >= 1 (Dirichlet) or >= 0 (Neumann).
std::vector<double> rectangleSpectrum(double a, double b, int count, bool neumann) {
    std::vector<double> v;
    const int lo = neumann ? 0 : 1;
    for (int m = lo; m <= 80; ++m)
        for (int n = lo; n <= 80; ++n) v.push_back(M_PI * M_PI * (m * m / (a * a) + n * n / (b * b)));
    std::sort(v.begin(), v.end());
    v.resize(count);
    return v;
}

// Squared zeros of J_m (Dirichlet) and J_m' (Neumann) on the unit disk, each
// with multiplicity 1 for m = 0 and 2 otherwise, ascending. Computed to 13
// digits by bisection on the power series of J_m in 50-digit arithmetic.
std::vector<double> diskSpectrum(bool neumann) {
    struct Z { double z; int m; };
    static const Z dirichlet[] = {
        {2.4048255576958, 0}, {3.8317059702075, 1}, {5.1356223018407, 2}, {5.5200781102863, 0},
        {6.3801618959240, 3}, {7.0155866698156, 1}, {7.5883424345038, 4}, {8.4172441403999, 2},
        {8.6537279129110, 0}, {8.7714838159600, 5}, {9.7610231299817, 3}, {9.9361095242177, 6},
        {10.1734681350627, 1}, {11.0647094885012, 4}, {11.0863700192451, 7}};
    static const Z neumannZ[] = {
        {0.0, 0}, {1.8411837813407, 1}, {3.0542369282271, 2}, {3.8317059702075, 0},
        {4.2011889412105, 3}, {5.3175531260840, 4}, {5.3314427735250, 1}, {6.4156163757002, 5},
        {6.7061331941585, 2}, {7.0155866698156, 0}, {7.5012661446841, 6}, {8.0152365983760, 3}};
    std::vector<double> v;
    if (neumann) {
        for (const Z &z : neumannZ)
            for (int k = 0; k < (z.m == 0 ? 1 : 2); ++k) v.push_back(z.z * z.z);
    } else {
        for (const Z &z : dirichlet)
            for (int k = 0; k < (z.m == 0 ? 1 : 2); ++k) v.push_back(z.z * z.z);
    }
    return v;
}

double maxRelativeError(const std::vector<double> &got, const std::vector<double> &want, int count,
                        int skip = 0) {
    double worst = 0.0;
    for (int i = skip; i < count && i < static_cast<int>(got.size()) && i < static_cast<int>(want.size()); ++i)
        worst = std::max(worst, std::fabs(got[i] - want[i]) / std::fabs(want[i]));
    return worst;
}

ShapeDNA::Options raw(int count, int degree, bool neumann = false) {
    ShapeDNA::Options o;
    o.count = count;
    o.degree = degree;
    o.boundary = neumann ? ShapeDNA::Neumann : ShapeDNA::Dirichlet;
    o.normalization = ShapeDNA::NoNormalization;
    o.denseBelow = 0;   // exercise Lanczos unless it cannot run
    return o;
}

// ---------------------------------------------------------------------------
// self-test
// ---------------------------------------------------------------------------

// The reference element: its integrals against things they must reproduce.
void testReferenceElement() {
    heading("Reference elements (Sec. 4: exact integrals on the reference triangle)");
    for (int p = 1; p <= 3; ++p) {
        const shapedna::ReferenceElement ref(p);
        const int n = ref.nodes;
        double massSum = 0.0, rowSum = 0.0, asym = 0.0;
        for (int a = 0; a < n; ++a) {
            double rx = 0.0, ry = 0.0, rxy = 0.0;
            for (int b = 0; b < n; ++b) {
                massSum += ref.mass[a * n + b];
                rx += ref.stiffXX[a * n + b];
                ry += ref.stiffYY[a * n + b];
                rxy += ref.stiffXY[a * n + b];
                asym = std::max({asym, std::fabs(ref.mass[a * n + b] - ref.mass[b * n + a]),
                                 std::fabs(ref.stiffXX[a * n + b] - ref.stiffXX[b * n + a]),
                                 std::fabs(ref.stiffXY[a * n + b] - ref.stiffXY[b * n + a])});
            }
            rowSum = std::max({rowSum, std::fabs(rx), std::fabs(ry), std::fabs(rxy)});
        }
        // The interpolant of a polynomial of degree <= p is the polynomial, so
        // its energy and mass are integrals of monomials: u = xi^p gives
        // int (du/dxi)^2 = p^2 (2p-2)! / (2p)! and int u^2 = (2p)! / (2p+2)!.
        std::vector<double> u(n);
        for (int k = 0; k < n; ++k) u[k] = std::pow(static_cast<double>(ref.lattice[k][0]) / p, p);
        double e = 0.0, m = 0.0;
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b) {
                e += u[a] * ref.stiffXX[a * n + b] * u[b];
                m += u[a] * ref.mass[a * n + b] * u[b];
            }
        auto fact = [](int k) { double f = 1; for (int i = 2; i <= k; ++i) f *= i; return f; };
        const double eWant = p * p * fact(2 * p - 2) / fact(2 * p);
        const double mWant = fact(2 * p) / fact(2 * p + 2);
        const bool ok = std::fabs(massSum - 0.5) < 1e-14 && rowSum < 1e-13 && asym < 1e-14 &&
                        std::fabs(e - eWant) < 1e-13 && std::fabs(m - mWant) < 1e-14;
        verdict(ok, "P" + std::to_string(p) + ": sum M = 1/2 (" + sci(std::fabs(massSum - 0.5)) +
                        "), stiffness rows sum to 0 (" + sci(rowSum) + "), xi^p energy " + sci(std::fabs(e - eWant)) +
                        " and mass " + sci(std::fabs(m - mWant)) + " off");
    }
    const shapedna::ReferenceElement p1(1);
    const double want[9] = {2, 1, 1, 1, 2, 1, 1, 1, 2};
    double off = 0.0;
    for (int k = 0; k < 9; ++k) off = std::max(off, std::fabs(p1.mass[k] - want[k] / 24.0));
    verdict(off < 1e-16, "P1 mass matrix is the textbook [2 1 1; 1 2 1; 1 1 2] / 24");
}

// The dense solvers against Eigen's.
void testDenseEigen() {
    heading("Dense symmetric eigensolver (tred2 + tql2) against Eigen");
    Rand rng(7);
    for (int n : {1, 2, 17, 120}) {
        Eigen::MatrixXd A(n, n), B(n, n);
        for (int j = 0; j < n; ++j)
            for (int i = j; i < n; ++i) A(i, j) = A(j, i) = rng.next();
        Eigen::MatrixXd G = Eigen::MatrixXd::NullaryExpr(n, n, [&]() { return rng.next(); });
        B = G * G.transpose() + n * Eigen::MatrixXd::Identity(n, n);

        std::vector<double> w, V;
        shapedna::symmetricEigen(n, A.data(), w, &V);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
        double dw = 0.0;
        for (int i = 0; i < n; ++i) dw = std::max(dw, std::fabs(w[i] - es.eigenvalues()(i)));
        Eigen::Map<Eigen::MatrixXd> Vm(V.data(), n, n);
        const double orth = (Vm.transpose() * Vm - Eigen::MatrixXd::Identity(n, n)).norm();
        const double res = (A * Vm - Vm * Eigen::Map<Eigen::VectorXd>(w.data(), n).asDiagonal()).norm();
        const double scale = std::max(1.0, A.norm());

        std::vector<double> gw, GV;
        shapedna::generalizedSymmetricEigen(n, A.data(), B.data(), gw, &GV);
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(A, B);
        double dg = 0.0;
        for (int i = 0; i < n; ++i) dg = std::max(dg, std::fabs(gw[i] - ges.eigenvalues()(i)));
        Eigen::Map<Eigen::MatrixXd> X(GV.data(), n, n);
        const double borth = (X.transpose() * B * X - Eigen::MatrixXd::Identity(n, n)).norm();

        verdict(dw < 1e-13 * scale && orth < 1e-12 && res < 1e-12 * scale && dg < 1e-12 * scale && borth < 1e-11,
                "n = " + std::to_string(n) + ": eigenvalues within " + sci(dw) + " of Eigen, V^T V = I to " +
                    sci(orth) + ", ||AV - VW|| " + sci(res) + "; generalized within " + sci(dg) +
                    ", X^T B X = I to " + sci(borth));
    }
}

// The sparse Cholesky: solves, and solves identically on one thread and many.
void testCholesky() {
    heading("Nested-dissection multifrontal Cholesky");
    const Mesh disk = makeDisk(1.0, 128, 0.05);
    for (int p = 1; p <= 3; ++p) {
        const shapedna::LagrangeSpace space(disk, p);
        const shapedna::ReferenceElement ref(p);
        const shapedna::FEMatrices fe = shapedna::assemble(space, ref, true);
        const shapedna::SparseMatrix &A = fe.stiffness;
        const int n = A.n;
        const int nrhs = 3;
        Rand rng(11);
        std::vector<double> x(static_cast<std::size_t>(n) * nrhs), b(x.size());
        for (double &v : x) v = rng.next();
        A.multiply(x.data(), b.data(), nrhs, n);

        shapedna::SparseCholesky chol;
        const bool ok = chol.factor(A, fe.dofPositions);
        std::vector<double> y = b;
        chol.solve(y.data(), nrhs, n);
        double err = 0.0, nx = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            err = std::max(err, std::fabs(y[i] - x[i]));
            nx = std::max(nx, std::fabs(x[i]));
        }
        // Every row appears in the ordering once.
        std::vector<int> seen(n, 0);
        for (int v : chol.permutation()) ++seen[v];
        const bool perm = std::all_of(seen.begin(), seen.end(), [](int c) { return c == 1; });

        bool same = true;
        if (shapedna::haveOpenMP()) {
            std::vector<double> y1 = b;
            shapedna::ScopedThreads one(1);
            shapedna::SparseCholesky serial;
            serial.factor(A, fe.dofPositions);
            serial.solve(y1.data(), nrhs, n);
            same = std::memcmp(y1.data(), y.data(), y.size() * sizeof(double)) == 0;
        }
        const shapedna::SparseCholesky::Stats &s = chol.getStats();
        const double fill = static_cast<double>(s.nnzL) / static_cast<double>(s.nnzA);
        verdict(ok && perm && err < 1e-10 * nx && same,
                "P" + std::to_string(p) + ", " + std::to_string(n) + " unknowns, " + std::to_string(s.supernodes) +
                    " supernodes, largest front " + std::to_string(s.maxFront) + ", fill " + fix(fill, 1) +
                    "x: 3 solves within " + sci(err / nx) + (shapedna::haveOpenMP() ? ", 1 thread bit-identical" : ""));
    }
}

// Sec. 8's first check: the rectangle.
void testRectangle() {
    heading("Rectangle [0,2] x [0,1] (Sec. 8: lambda = pi^2 (m^2/4 + n^2))");
    const int count = 20;
    const std::vector<double> dir = rectangleSpectrum(2.0, 1.0, count, false);
    const std::vector<double> neu = rectangleSpectrum(2.0, 1.0, count + 1, true);
    // The first 20 at h = 1/20 come out within 4.7e-2, 3.4e-4 and 1.1e-6 for
    // P1-P3 (2026-09-25): discretisation error, as the orders below confirm.
    // The tolerances sit about 2-3x above that.
    const double tol[4] = {0.0, 1e-1, 1e-3, 3e-6};
    for (int p = 1; p <= 3; ++p) {
        const Mesh m = makeRectangle(2.0, 1.0, 40, 20);
        ShapeDNA d(m, raw(count, p));
        const bool ok = d.run();
        const double e = maxRelativeError(d.eigenvalues(), dir, count);
        ShapeDNA nm(m, raw(count, p, true));
        const bool okN = nm.run();
        const double eN = maxRelativeError(nm.eigenvalues(), neu, count + 1, 1);
        verdict(ok && okN && e < tol[p] && eN < tol[p] && std::fabs(nm.eigenvalues()[0]) < 1e-9,
                "P" + std::to_string(p) + ": Dirichlet first 20 within " + sci(e) + ", Neumann within " + sci(eN) +
                    " (lambda_1 = " + sci(nm.eigenvalues()[0]) + "), " + std::to_string(d.getReport().restarts) +
                    " restarts");
    }

    // Convergence: halving h should cut the error in lambda_1 by 2^(2p).
    for (int p = 1; p <= 3; ++p) {
        double err[2];
        for (int k = 0; k < 2; ++k) {
            const int nx = (k == 0) ? 10 : 20;
            const Mesh m = makeRectangle(2.0, 1.0, nx, nx / 2);
            ShapeDNA d(m, raw(4, p));
            d.run();
            err[k] = std::fabs(d.eigenvalues()[0] - dir[0]) / dir[0];
        }
        const double order = std::log2(err[0] / err[1]);
        verdict(order > 2 * p - 0.25,
                "P" + std::to_string(p) + ": lambda_1 error " + sci(err[0]) + " -> " + sci(err[1]) +
                    " as h halves, observed order " + fix(order, 2) + " (expected " + std::to_string(2 * p) + ")");
    }
}

// Exact double eigenvalues: what the block buys.
void testSquareMultiplicity() {
    heading("Unit square, criss-cross mesh: exact multiplicities");
    const int count = 24;
    const std::vector<double> want = rectangleSpectrum(1.0, 1.0, count, false);
    const Mesh m = makeCrissCross(24);
    for (int block : {4, 2, 1}) {
        ShapeDNA::Options o = raw(count, 2);
        o.blockSize = block;
        ShapeDNA d(m, o);
        const bool ok = d.run();
        const double e = maxRelativeError(d.eigenvalues(), want, count);
        // Pairs that are exactly double in the discrete problem must come out
        // equal to the solver's accuracy.
        double split = 0.0;
        for (int i = 0; i + 1 < count; ++i)
            if (std::fabs(want[i + 1] - want[i]) < 1e-12 * want[i])
                split = std::max(split, std::fabs(d.eigenvalues()[i + 1] - d.eigenvalues()[i]) / want[i]);
        const std::string what = "block " + std::to_string(block) + ": first " + std::to_string(count) +
                                 " within " + sci(e) + " of pi^2 (m^2 + n^2), doubles split by " + sci(split) +
                                 ", " + std::to_string(d.getReport().solves) + " solves";
        if (block >= 2) verdict(ok && e < 2e-3 && split < 1e-9, what);
        else info(what + " (single-vector Lanczos: not asserted)");
    }
}

// Sec. 8's second check: the disk.
void testDisk() {
    heading("Unit disk (Sec. 8: lambda = squared zeros of J_m and J_m')");
    // 256 chords: the inscribed polygon has (2 pi / 256)^2 / 6 = 1e-4 less area
    // than the disk, and its eigenvalues are that much higher. That is the
    // geometry the mesh has, and the floor of the comparison.
    const Mesh m = makeDisk(1.0, 256, 0.04);
    const std::vector<double> dir = diskSpectrum(false), neu = diskSpectrum(true);
    ShapeDNA d(m, raw(static_cast<int>(dir.size()), 3));
    const bool ok = d.run();
    const double e = maxRelativeError(d.eigenvalues(), dir, static_cast<int>(dir.size()));
    ShapeDNA n(m, raw(static_cast<int>(neu.size()) - 1, 3, true));
    const bool okN = n.run();
    const double eN = maxRelativeError(n.eigenvalues(), neu, static_cast<int>(neu.size()), 1);
    verdict(ok && e < 5e-4, "Dirichlet: first " + std::to_string(dir.size()) + " within " + sci(e) + " (" +
                                std::to_string(d.getReport().unknowns) + " P3 unknowns)");
    verdict(okN && eN < 5e-4 && std::fabs(n.eigenvalues()[0]) < 1e-9,
            "Neumann: first " + std::to_string(neu.size()) + " within " + sci(eN) + ", lambda_1 = " +
                sci(n.eigenvalues()[0]) + " (left out of the DNA, which starts at " + fix(n.dna()[0], 5) +
                ", j'_11^2 = 3.38996)");
}

// Sec. 8's third check: the heat trace hands back the geometry.
void testHeatTrace() {
    heading("Heat trace of the rectangle (Sec. 8: area, perimeter, corners)");
    const Mesh m = makeRectangle(2.0, 1.0, 40, 20);
    ShapeDNA d(m, raw(200, 3));
    const bool ok = d.run();
    const ShapeDNA::HeatTraceFit fit = d.fitHeatTrace();
    const ShapeDNA::Report &r = d.getReport();
    const double ea = std::fabs(fit.area - 2.0) / 2.0;
    const double el = std::fabs(fit.boundaryLength - 6.0) / 6.0;
    const double ec = std::fabs(fit.constant - 0.25);
    verdict(std::fabs(r.cornerTerm - 0.25) < 1e-12, "the mesh's corner term is 4 x (pi^2 - pi^2/4) / (12 pi^2) = 1/4 (" +
                                                        sci(std::fabs(r.cornerTerm - 0.25)) + ")");
    verdict(ok && fit.valid && ea < 1e-3 && el < 1e-2 && ec < 2e-2,
            "200 eigenvalues, t in [" + sci(fit.tMin) + ", " + sci(fit.tMax) + "]: area " + fix(fit.area, 5) +
                " (2), boundary " + fix(fit.boundaryLength, 4) + " (6), constant " + fix(fit.constant, 4) + " (0.25)");
}

// The DNA is a property of the shape, not of its placement or its mesh.
void testInvariance() {
    heading("Invariance (isometry, scale, remeshing) and discrimination");
    const Mesh a = makeDisk(1.0, 128, 0.06, 0.35);   // an annulus
    ShapeDNA::Options o;
    o.count = 30;
    o.degree = 2;
    ShapeDNA base(a, o);
    base.run();
    const std::vector<double> ref = base.dna();
    double norm = 0.0;
    for (double x : ref) norm += x * x;
    norm = std::sqrt(norm);

    const Mesh moved = similar(a, 0.7, 3.1, 5.0, -2.0, false);
    const Mesh mirrored = similar(a, -1.2, 0.2, -1.0, 4.0, true);
    ShapeDNA sm(moved, o), sr(mirrored, o);
    sm.run();
    sr.run();
    const double dm = ShapeDNA::distance(ref, sm.dna()) / norm;
    const double dr = ShapeDNA::distance(ref, sr.dna()) / norm;
    verdict(dm < 1e-9 && dr < 1e-9, "rotated, scaled x3.1 and moved: " + sci(dm) + "; mirrored and scaled x0.2: " +
                                        sci(dr) + " (relative distance)");

    const Mesh remeshed = makeDisk(1.0, 128, 0.045, 0.35);
    ShapeDNA sq(remeshed, o);
    sq.run();
    const double dq = ShapeDNA::distance(ref, sq.dna()) / norm;
    // A thicker annulus: a different shape, and it should read as one.
    const Mesh other = makeDisk(1.0, 128, 0.06, 0.45);
    ShapeDNA so(other, o);
    so.run();
    const double dother = ShapeDNA::distance(ref, so.dna()) / norm;
    verdict(dq < 2e-3 && dother > 20 * dq, "a different mesh of the same annulus: " + sci(dq) +
                                               "; inner radius 0.35 -> 0.45: " + sci(dother));
}

// Lanczos against the dense solver and Eigen, on one small problem.
void testLanczosVsDense() {
    heading("Block Lanczos against the dense solvers");
    const Mesh m = makeL(0.08);
    for (bool neumann : {false, true}) {
        ShapeDNA::Options o = raw(25, 2, neumann);
        o.eigenfunctions = true;
        ShapeDNA lan(m, o);
        const bool okL = lan.run();
        o.denseBelow = 1 << 30;
        ShapeDNA den(m, o);
        const bool okD = den.run();

        // Eigen on the same assembled matrices, as a third opinion.
        const shapedna::LagrangeSpace space(m, 2);
        const shapedna::ReferenceElement ref(2);
        const shapedna::FEMatrices fe = shapedna::assemble(space, ref, !neumann);
        const int n = fe.stiffness.n;
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n), B = A;
        for (int i = 0; i < n; ++i)
            for (int k = fe.stiffness.rowPtr[i]; k < fe.stiffness.rowPtr[i + 1]; ++k) {
                A(i, fe.stiffness.col[k]) = fe.stiffness.val[k];
                B(i, fe.mass.col[k]) = fe.mass.val[k];
            }
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(A, B);

        // Relative error, except that the Neumann zero mode is measured against
        // the first non-zero eigenvalue: every method gets it to within about
        // eps * lambda_max absolutely (the dense one to 6e-11 here, Lanczos to
        // 1e-12), and relative to zero that means nothing.
        double dd = 0.0, de = 0.0, dde = 0.0;
        const int k = static_cast<int>(lan.eigenvalues().size());
        const double floorScale = den.eigenvalues()[neumann ? 1 : 0];
        for (int i = 0; i < k; ++i) {
            const double s = std::max(floorScale, std::fabs(den.eigenvalues()[i]));
            dd = std::max(dd, std::fabs(lan.eigenvalues()[i] - den.eigenvalues()[i]) / s);
            de = std::max(de, std::fabs(lan.eigenvalues()[i] - ges.eigenvalues()(i)) / s);
            dde = std::max(dde, std::fabs(den.eigenvalues()[i] - ges.eigenvalues()(i)) / s);
        }
        const ShapeDNA::Report &r = lan.getReport();
        verdict(okL && okD && !r.dense && den.getReport().dense && dd < 1e-10 && de < 1e-10 &&
                    r.maxTrueResidual < 1e-8,
                std::string(neumann ? "Neumann" : "Dirichlet") + ", " + std::to_string(n) + " unknowns: Lanczos within " +
                    sci(dd) + " of the dense solver and " + sci(de) + " of Eigen (dense vs Eigen " + sci(dde) +
                    "); ||Ax - lambda Bx|| " + sci(r.maxTrueResidual) + ", " + std::to_string(r.restarts) + " restarts");
    }
}

// The same answer on one thread and on all of them.
void testThreads() {
    heading("Thread count");
    if (!shapedna::haveOpenMP()) {
        info("built without OpenMP; nothing to compare");
        return;
    }
    const Mesh m = makeDisk(1.0, 200, 0.03, 0.3);
    ShapeDNA::Options o;
    o.count = 50;
    o.threads = 1;
    ShapeDNA one(m, o);
    auto t0 = std::chrono::steady_clock::now();
    one.run();
    const double t1 = seconds(t0);
    o.threads = 0;
    ShapeDNA all(m, o);
    t0 = std::chrono::steady_clock::now();
    all.run();
    const double tn = seconds(t0);
    const bool same = one.eigenvalues().size() == all.eigenvalues().size() &&
                      std::memcmp(one.eigenvalues().data(), all.eigenvalues().data(),
                                  one.eigenvalues().size() * sizeof(double)) == 0;
    verdict(same, "P3 annulus, " + std::to_string(all.getReport().unknowns) + " unknowns, 50 eigenvalues: 1 thread and " +
                      std::to_string(all.getReport().threads) + " threads bit-identical (" + fix(t1, 2) + " s vs " +
                      fix(tn, 2) + " s)");
}

// Pieces and holes.
void testTopology() {
    heading("Topology: pieces, holes, zero modes");
    const Mesh ann = makeDisk(1.0, 96, 0.08, 0.4);
    ShapeDNA a(ann, raw(10, 2));
    a.run();
    verdict(a.getReport().eulerCharacteristic == 0 && a.getReport().components == 1,
            "annulus: chi = " + std::to_string(a.getReport().eulerCharacteristic) + ", corner term " +
                fix(a.getReport().cornerTerm, 4) + " (chi / 6 = 0 for a resolved smooth boundary)");

    const Mesh sq = makeRectangle(1.0, 1.0, 12, 12);
    const Mesh two = disjointUnion(sq, makeRectangle(2.0, 1.0, 24, 12), 1.5);
    ShapeDNA d(two, raw(12, 2, true));
    const bool ok = d.run();
    // The union's Neumann spectrum is the two spectra merged.
    std::vector<double> want = rectangleSpectrum(1.0, 1.0, 12, true);
    const std::vector<double> w2 = rectangleSpectrum(2.0, 1.0, 12, true);
    want.insert(want.end(), w2.begin(), w2.end());
    std::sort(want.begin(), want.end());
    const double e = maxRelativeError(d.eigenvalues(), want, 12, 2);
    verdict(ok && d.getReport().components == 2 && d.getReport().zeroModes == 2 && d.dna().size() == 12 &&
                std::fabs(d.eigenvalues()[1]) < 1e-9 && e < 1e-3,
            "two separate rectangles, Neumann: 2 zero modes left out of the DNA, the rest the merged spectra to " +
                sci(e));
}

int selfTest() {
    std::cout << "ShapeDNA self-test (" << (shapedna::haveOpenMP() ? "OpenMP, " : "no OpenMP, ")
              << shapedna::maxThreads() << " threads)\n";
    const auto t0 = std::chrono::steady_clock::now();
    testReferenceElement();
    testDenseEigen();
    testCholesky();
    testRectangle();
    testSquareMultiplicity();
    testDisk();
    testHeatTrace();
    testInvariance();
    testLanczosVsDense();
    testThreads();
    testTopology();
    std::cout << "\n";
    if (gFailures == 0) std::cout << "  " << kPass << " all checks passed";
    else std::cout << "  " << kFail << " " << gFailures << " check(s) failed";
    std::cout << " (" << fix(seconds(t0), 1) << " s)\n";
    return gFailures == 0 ? 0 : 1;
}

// ---------------------------------------------------------------------------
// the driver
// ---------------------------------------------------------------------------

// The file name with its directory, "singlemat/geom001.obj": the corpus uses
// the same names in singlemat/ and multimat/.
std::string baseName(const std::string &path) {
    const std::size_t s = path.find_last_of('/');
    if (s == std::string::npos || s == 0) return path;
    const std::size_t d = path.find_last_of('/', s - 1);
    return d == std::string::npos ? path : path.substr(d + 1);
}

void printVector(const std::vector<double> &v, int digits) {
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (i % 8 == 0) std::cout << (i == 0 ? "   " : "\n   ");
        std::cout << " " << std::setw(digits + 6) << std::fixed << std::setprecision(digits) << v[i];
    }
    std::cout << "\n";
}

}  // namespace

int main(int argc, char **argv) {
    std::vector<std::string> files;
    ShapeDNA::Options opts;
    bool heat = false, printRaw = false;
    std::string csv;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) {
                std::cerr << a << " needs a value\n";
                std::exit(2);
            }
            return argv[++i];
        };
        if (a == "--selftest") return selfTest();
        else if (a == "--count") opts.count = std::stoi(next());
        else if (a == "--degree") opts.degree = std::stoi(next());
        else if (a == "--refine") opts.refine = std::stoi(next());
        else if (a == "--neumann") opts.boundary = ShapeDNA::Neumann;
        else if (a == "--threads") opts.threads = std::stoi(next());
        else if (a == "--block") opts.blockSize = std::stoi(next());
        else if (a == "--basis") opts.basisSize = std::stoi(next());
        else if (a == "--tol") opts.tolerance = std::stod(next());
        else if (a == "--leaf") opts.leafSize = std::stoi(next());
        else if (a == "--dense-below") opts.denseBelow = std::stoi(next());
        else if (a == "--heat") heat = true;
        else if (a == "--eigenfunctions") opts.eigenfunctions = true;
        else if (a == "--raw") printRaw = true;
        else if (a == "--csv") csv = next();
        else if (a == "--norm") {
            const std::string n = next();
            if (n == "area") opts.normalization = ShapeDNA::AreaNormalization;
            else if (n == "first") opts.normalization = ShapeDNA::FirstEigenvalue;
            else if (n == "weyl") opts.normalization = ShapeDNA::WeylSlope;
            else if (n == "weyl-ratio") opts.normalization = ShapeDNA::WeylRatio;
            else if (n == "none") opts.normalization = ShapeDNA::NoNormalization;
            else {
                std::cerr << "unknown normalisation " << n << " (area|first|weyl|weyl-ratio|none)\n";
                return 2;
            }
        } else if (!a.empty() && a[0] == '-') {
            std::cerr << "unknown option " << a << "\n";
            return 2;
        } else files.push_back(a);
    }
    if (files.empty()) {
        std::cerr << "usage: TestShapeDNA <mesh.obj> [more.obj ...] [options] | --selftest\n";
        return 2;
    }

    const char *normName[] = {"raw", "area-normalised", "divided by the first", "divided by the Weyl slope",
                              "area-normalised over 4 pi k"};
    std::vector<std::vector<double>> dnas;
    std::vector<std::string> names;
    int failures = 0;
    for (const std::string &f : files) {
        Mesh mesh;
        try {
            mesh = Mesh(f);
        } catch (const std::exception &e) {
            std::cerr << f << ": " << e.what() << "\n";
            ++failures;
            continue;
        }
        ShapeDNA sd(mesh, opts);
        const bool ok = sd.run();
        const ShapeDNA::Report &r = sd.getReport();
        std::cout << baseName(f) << ": " << r.triangles << " triangles, P" << opts.degree << ", "
                  << (opts.boundary == ShapeDNA::Dirichlet ? "Dirichlet" : "Neumann") << ", " << r.unknowns
                  << " unknowns, " << r.threads << (r.openMP ? " threads" : " thread (no OpenMP)") << "\n";
        std::cout << "  domain   area " << fix(r.area, 6) << ", boundary " << fix(r.boundaryLength, 6) << ", "
                  << r.components << " piece(s), chi " << r.eulerCharacteristic << "\n";
        if (r.dense) {
            std::cout << "  solve    dense, " << fix(r.secondsEigen, 2) << " s\n";
        } else {
            std::cout << "  factor   " << r.supernodes << " supernodes (depth " << r.treeDepth << "), largest front "
                      << r.maxFront << ", nnz(L) " << sci(static_cast<double>(r.nnzL)) << ", "
                      << fix(r.factorFlops * 1e-9, 2) << " GFlop; order " << fix(r.secondsOrder, 3) << " s, factor "
                      << fix(r.secondsFactor, 3) << " s\n";
            std::cout << "  Lanczos  basis " << r.basisSize << " (block " << r.blockSize << "), " << r.restarts
                      << " restarts, " << r.blockSteps << " block steps, " << r.solves << " solves, "
                      << r.deflations << " deflations, max residual " << sci(r.maxRitzResidual) << ", "
                      << fix(r.secondsEigen, 3) << " s\n";
        }
        if (opts.eigenfunctions) std::cout << "  check    max ||Ax - lambda Bx|| / (lambda ||Bx||) = " << sci(r.maxTrueResidual) << "\n";
        std::cout << "  total    assemble " << fix(r.secondsAssemble, 3) << " s, all " << fix(r.seconds, 3) << " s\n";
        for (const std::string &m : r.messages) std::cout << "  note     " << m << "\n";
        if (heat) {
            const ShapeDNA::HeatTraceFit fit = sd.fitHeatTrace();
            if (fit.valid)
                std::cout << "  heat     t in [" << sci(fit.tMin) << ", " << sci(fit.tMax) << "]: area "
                          << fix(fit.area, 5) << ", boundary " << fix(fit.boundaryLength, 4) << ", constant "
                          << fix(fit.constant, 4) << " (corner term of the mesh " << fix(r.cornerTerm, 4)
                          << ", chi/6 = " << fix(r.eulerCharacteristic / 6.0, 4) << ")\n";
            else
                std::cout << "  heat     too few eigenvalues to fit\n";
        }
        if (!ok) {
            ++failures;
            std::cout << "  FAILED\n\n";
            continue;
        }
        if (printRaw) {
            std::cout << "  eigenvalues (" << sd.eigenvalues().size() << "):\n";
            printVector(sd.eigenvalues(), 6);
        }
        std::cout << "  Shape-DNA (" << sd.dna().size() << ", " << normName[opts.normalization] << "):\n";
        printVector(sd.dna(), 4);
        std::cout << "\n";
        dnas.push_back(sd.dna());
        names.push_back(baseName(f));
    }

    if (dnas.size() > 1) {
        std::cout << "Pairwise distances ||DNA_i - DNA_j||:\n";
        std::size_t w = 8;
        for (const std::string &n : names) w = std::max(w, n.size());
        std::cout << std::string(w + 2, ' ');
        for (std::size_t j = 0; j < names.size(); ++j) std::cout << std::setw(10) << j;
        std::cout << "\n";
        for (std::size_t i = 0; i < names.size(); ++i) {
            std::cout << std::setw(3) << i << " " << std::left << std::setw(w - 2) << names[i] << std::right;
            for (std::size_t j = 0; j < names.size(); ++j)
                std::cout << std::setw(10) << std::fixed << std::setprecision(3) << ShapeDNA::distance(dnas[i], dnas[j]);
            std::cout << "\n";
        }
    }

    if (!csv.empty()) {
        std::ofstream out(csv);
        for (std::size_t i = 0; i < dnas.size(); ++i) {
            out << names[i];
            for (double x : dnas[i]) out << "," << std::setprecision(17) << x;
            out << "\n";
        }
    }
    return failures == 0 ? 0 : 1;
}
