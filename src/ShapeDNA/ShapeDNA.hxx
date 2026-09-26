#ifndef __SHAPEDNA_SHAPEDNA_HXX__
#define __SHAPEDNA_SHAPEDNA_HXX__

#include <limits>
#include <string>
#include <vector>

#include "mesh/Mesh.hxx"

// The Shape-DNA of a planar domain (Reuter, Wolter and Peinecke, CAD 38, 2006;
// docs/shape_dna.md): the first n eigenvalues of the Laplacian on the domain,
//
//     -Delta f = lambda f,
//
// normalised for scale. In the plane the Laplace-Beltrami operator is the
// ordinary Laplacian and isometry is congruence, so the vector is invariant
// under rotating, translating and reflecting the domain and, once normalised,
// under scaling it -- and it does not care how the domain was meshed, beyond
// how well the mesh resolves it. Two shapes are compared by the Euclidean
// distance between their vectors.
//
// The steps of the doc, and where each lives:
//
//   1  Pre-process: the Mesh is the representation. Its material ids are
//      ignored -- the spectrum is of the whole domain, interfaces and all --
//      and a mesh of several pieces gives the union of their spectra.
//   2  Refine: Options::refine splits every edge into that many parts
//      (refineUniform).
//   3-4 Weak form and discretisation: Lagrange elements of degree 1-3 with
//      exact reference integrals (LagrangeFEM). Cubic is the default, as in
//      the paper.
//   5  Boundary conditions: Dirichlet (the default) drops the nodes on dS from
//      the unknowns; Neumann keeps them and has lambda_1 = 0, once per
//      connected piece of the mesh, which the DNA leaves out.
//   6  Reorder: nested dissection, not Cuthill-McKee (see SparseCholesky).
//   7  Solve: shift-invert block Lanczos (Lanczos), our own, on our own sparse
//      Cholesky. A problem too small for a Krylov basis is solved densely.
//   8  Verify: fitHeatTrace() reads the area, the boundary length and the
//      constant term back off the spectrum; the TestShapeDNA self-test checks
//      the rectangle and the disk against their closed forms.
//   9  Normalise and truncate: Options::normalization, Options::count.
//
// Every stage that scales with the mesh runs on OpenMP threads, and none of
// them lets the thread count change the answer (Parallel.hxx).
namespace shapedna {

class ShapeDNA {
public:
    enum Boundary : int {
        Dirichlet = 0,   // f = 0 on dS
        Neumann = 1      // df/dn = 0 on dS
    };

    // Sec. 9's three ways of removing scale; d = 2 throughout.
    enum Normalization : int {
        NoNormalization = 0,
        // lambda * |Omega|. Dimensionless, and the default: the one that needs
        // nothing from the spectrum itself.
        AreaNormalization = 1,
        // lambda / lambda_1, lambda_1 the first eigenvalue in the DNA (the
        // first non-zero one under Neumann). The DNA then starts with 1.
        FirstEigenvalue = 2,
        // lambda / c, c the least-squares slope of lambda_k ~ c k over the DNA
        // (Weyl's law: c -> 4 pi / |Omega| as k grows).
        WeylSlope = 3,
        // lambda_k |Omega| / (4 pi k): the area-normalised DNA over Weyl's
        // leading term, k counted from 1 over the DNA (zero modes left out).
        // Scale-free like AreaNormalization, but with the linear trend of
        // Weyl's law divided out, so every entry sits near 1 and what differs
        // between shapes is what the leading term does not know -- the
        // boundary term, ~ |dS| / sqrt(|Omega| k) under Dirichlet, and below
        // it. Without it the n-th entry is ~n times the first on every shape,
        // and a learner reading the DNA as features mostly learns the index.
        // Report::normalizer is 1 / |Omega| here; the 4 pi k is per entry.
        WeylRatio = 4
    };

    struct Options {
        int count = 50;                 // n, the length of the DNA
        int degree = 3;                 // of the Lagrange elements, 1-3
        int refine = 1;                 // split each edge into this many parts first
        Boundary boundary = Dirichlet;
        Normalization normalization = AreaNormalization;

        // The eigensolver. `shift` is sigma in (A - sigma B)^{-1} B; NaN picks
        // 0 under Dirichlet and -0.01 * 4 pi / |Omega| under Neumann, where A
        // is singular and sigma must sit below the zero eigenvalue.
        int blockSize = 4;
        int basisSize = 0;              // 0: 2 (count + zeros) + 2 blockSize
        double tolerance = 1e-10;
        int maxRestarts = 300;
        double shift = std::numeric_limits<double>::quiet_NaN();
        unsigned long long seed = 1;
        // At most this many unknowns and the problem is solved densely instead.
        int denseBelow = 400;

        int leafSize = 64;              // SparseCholesky::Options::leafSize

        // Also compute the eigenfunctions, at the mesh's vertices. Needed for
        // anything but the DNA itself: plotting, and the true residuals in the
        // report.
        bool eigenfunctions = false;

        // Number of OpenMP threads; 0 leaves it to the runtime. The spectrum is
        // bit-identical whatever this is.
        int threads = 0;
    };

    struct Report {
        bool ran = false;
        bool converged = false;
        bool dense = false;             // solved by the dense fallback

        // The discretisation.
        int triangles = 0;              // after refinement
        int nodes = 0;                  // Lagrange nodes
        int unknowns = 0;               // after the boundary condition
        long long nnz = 0;              // of the stiffness matrix
        int degenerateElements = 0;

        // The domain, as meshed. `cornerTerm` is sum over boundary vertices of
        // (pi^2 - alpha^2) / (24 pi alpha), alpha the interior angle there: the
        // constant term of the heat trace of the mesh's polygon (van den Berg
        // and Srisatkunarajah 1988), which is 1/4 for a rectangle and tends to
        // chi / 6 as a smooth boundary is resolved.
        double area = 0.0;
        double boundaryLength = 0.0;
        int components = 0;             // connected pieces
        int eulerCharacteristic = 0;    // V - E + F = pieces - holes
        double cornerTerm = 0.0;

        // The solver.
        double shift = 0.0;
        int zeroModes = 0;              // left out of the DNA (Neumann)
        int eigenvalues = 0;            // computed
        int supernodes = 0, maxFront = 0, treeDepth = 0;
        long long nnzL = 0;
        double factorFlops = 0.0;
        int basisSize = 0, blockSize = 0;
        int restarts = 0, blockSteps = 0, solves = 0, deflations = 0;
        // The largest || Op x - theta x ||_B / theta over the computed pairs:
        // what the Lanczos iteration converged on.
        double maxRitzResidual = 0.0;
        // With Options::eigenfunctions, the largest || A x - lambda B x || /
        // (lambda || B x ||) measured on the assembled matrices directly.
        double maxTrueResidual = 0.0;
        double normalizer = 1.0;        // DNA = raw / normalizer (WeylRatio: also / (4 pi k))

        int threads = 1;
        bool openMP = false;
        double secondsAssemble = 0.0, secondsOrder = 0.0, secondsFactor = 0.0,
               secondsEigen = 0.0, seconds = 0.0;
        std::vector<std::string> messages;
    };

    explicit ShapeDNA(const Mesh &mesh);
    ShapeDNA(const Mesh &mesh, const Options &opts);

    // Assemble, factor, solve, normalise. False when the eigensolver did not
    // converge or the problem could not be set up; the report says which.
    bool run();

    // The spectrum as computed, ascending, zero modes included: count + zeros
    // eigenvalues of the discretised problem.
    const std::vector<double> &eigenvalues() const { return lambda; }
    // The Shape-DNA: `count` normalised eigenvalues, zero modes left out.
    const std::vector<double> &dna() const { return shapeDNA; }
    // Eigenfunction k at mesh vertex v is eigenfunctions()[v + k * |V|],
    // normalised to int f^2 = 1 and zero on dS under Dirichlet. Empty unless
    // Options::eigenfunctions.
    const std::vector<double> &eigenfunctions() const { return functions; }

    const Report &getReport() const { return report; }
    const Options &getOptions() const { return options; }

    // || a - b ||_2 over the first min(|a|, |b|) entries: Sec. 9's comparison.
    static double distance(const std::vector<double> &a, const std::vector<double> &b);

    // Z(t) = sum_i exp(-lambda_i t) over the computed spectrum.
    double heatTrace(double t) const;

    // Sec. 8's check. For small t,
    //
    //     4 pi t Z(t) ~ c0 + c1 sqrt(t) + c2 t,
    //
    // with c0 the area, c1 = -/+ (sqrt(pi) / 2) x boundary length (Dirichlet /
    // Neumann), and c2 / (4 pi) the constant term. Fitted by least squares on
    // log-spaced t in [16, 64] / lambda_max: small enough t for the expansion,
    // large enough that the eigenvalues not computed contribute less than
    // e^-16. A truncated spectrum can only see the domain at the scale
    // sqrt(t), so the recovered numbers get better with more eigenvalues.
    struct HeatTraceFit {
        bool valid = false;
        double area = 0.0;
        double boundaryLength = 0.0;
        double constant = 0.0;
        double tMin = 0.0, tMax = 0.0;
        double rmsResidual = 0.0;       // of the fit, relative to c0
    };
    HeatTraceFit fitHeatTrace() const;

private:
    const Mesh &mesh;
    Options options;
    Report report;
    std::vector<double> lambda, shapeDNA, functions;
};

}  // namespace shapedna

#endif  // __SHAPEDNA_SHAPEDNA_HXX__
