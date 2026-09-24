#ifndef __MESH_TMOP_HXX__
#define __MESH_TMOP_HXX__

#include <cstdint>
#include <string>
#include <vector>
#include "QuadMesh.hxx"

// An in-house Target-Matrix Optimization Paradigm (TMOP) smoother for
// mesh::QuadMesh: a node-local Newton relaxation of the TMOP energy, run as a
// coloured Gauss-Seidel sweep so it parallelises with OpenMP without changing
// its answer.
//
// This is not a binding to MFEM's TMOP. It is the same formulation solved a
// different way -- MFEM assembles a global nonlinear form and hands it to a
// Newton solver; this moves one node at a time, exactly, and mirrors the
// hand-rolled relaxations already in this codebase (TORSION::relaxToKernel,
// DiskTemplate::smooth). What carries over from MFEM is the mathematics and,
// where they coincide, the metric numbering, so `-mid 2` there and
// `Metric::Shape002` here mean the same function.
//
// ## The energy
//
// Each element has a **target Jacobian** W: the 2x2 matrix mapping the
// reference unit square to the element this optimizer would like to see there.
// The mesh carries one per quad in `QuadMesh::targetJacobian`; `prepare()`
// fills it from Options::target when it is empty or when the caller asks for a
// particular one. Against the element's actual Jacobian A the optimizer forms
//
//     T = A W^{-1}
//
// -- the departure from the target, with the target itself divided out, so
// T = I is a perfect element whatever shape or size W asked for -- and
// minimises
//
//     F(x) = sum_q  det(W_q)  sum_k  w_k  mu( T(xi_k, eta_k) )
//
// over the mesh's movable nodes, with (xi_k, w_k) a quadrature rule on the
// reference square and mu one of the metrics below. det(W_q) is the element's
// target area, so the sum is an integral over the *target* configuration and
// elements are weighted by the area they are meant to occupy rather than the
// one they currently do.
//
// ## Which nodes move
//
// Whatever `QuadMesh::nodeType` says, through `QuadMesh::projectStep`: interior
// nodes take their full 2-D Newton step, nodes on a smooth stretch of dS or of
// a material interface slide along that feature, corners, junctions, caller
// pins and one anchor per otherwise-free feature loop do not move at all.
//
// Sliding matters more than it looks. Pinning the whole boundary freezes its
// discretisation density while the interior metric pulls towards its own ideal
// spacing, and the whole mismatch lands in the one ring of elements touching
// the boundary. See QuadMesh::Options::fixAllFeatureNodes.
//
// ## Sliding on a curve
//
// A feature node's degree of freedom is the **parameter** of the curve it sits
// on (QuadMesh::Options::curveSource, the default), not a direction in the
// plane. The local Newton system is still solved in the plane and still
// restricted to the tangent, but what comes out is read as an arc length,
// divided by |C'(u)| to become a step in u, and the node is placed at C(u + du)
// -- so it is on the curve after the step, and after every halving of the step
// the line search tries, since the backtracking halves du and not a chord.
//
// The alternative, and what this did before, is to move along the chord through
// the node's two feature neighbours, refreshing that chord before each node's
// own solve. That keeps the enclosed area exact but leaves the node one
// sagitta off the feature every time, always on the inside of a turn: a
// discretised circle relaxes towards its own inscribed polygon and keeps
// going. It is still available (CurveChord) and is what the area invariant in
// the Report below is written for.
//
// ## Untangling
//
// Every barrier metric is +infinity on an inverted element, so a mesh that
// starts with one cannot be started on directly. Options::untangle runs an
// untangling phase first, on a metric that is finite on an inverted element,
// and hands over to the main metric once nothing is inverted -- the standard
// two-phase arrangement, on by default. Since 2026-09-24 that metric is
// UntangleRegularized (Garanzha et al.'s regularised shape + size energy)
// rather than Untangle022 with a tracked tau_0, which crushes a fold instead
// of opening it; Options::untangler says why and keeps the old one.
//
// Element corners that no admissible move can change and that are flat or
// inverted -- a block corner a layout put on a straight feature -- are
// "frozen" and left out of the energy altogether (findFrozenCorners), so one
// of them neither blocks the smoothing of the rest nor, through the |T|^2 it
// leaves as the only thing to lower, pulls its neighbours onto itself.
namespace mesh {

class TMOP {
public:
    // The metrics, numbered as in MFEM where the two coincide. tau = det T and
    // |T| is the Frobenius norm; each is written so that mu(I) = 0 and mu >= 0.
    //
    // "Barrier" means mu -> +infinity as tau -> 0: such a metric cannot produce
    // an inverted element from a valid one, but also cannot be evaluated on a
    // mesh that already has one, which is what Options::untangle is for.
    enum Metric : int {
        // |T|^2 - 2 tau. Shape, no barrier. Zero exactly when T is a rotation
        // (of any scale: |T|^2 >= 2 tau with equality only there), so it is a
        // shape metric like 002 but finite on an inverted element -- the one
        // here that can be evaluated on a tangled mesh without a tau_0. Nothing
        // in it stops an element flattening on the way, which is why untangling
        // goes through a metric with a (shifted or regularised) barrier.
        Shape004 = 4,
        // 0.5 |T|^2 / tau - 1. **The default.** Pure shape: invariant under
        // scaling T, so it says nothing about how big an element is and
        // everything about whether it is square. Barrier.
        Shape002 = 2,
        // |T|^2 (1 + 1/tau^2) - 4, i.e. |T - T^-T|^2. Shape and size together,
        // barrier. Use when the target size is meaningful -- a per-element size
        // field from the layout, say -- and not just a placeholder.
        ShapeSize007 = 7,
        // 0.5 (|T|^2 - 2 tau) / (tau - tau_0). An untangler: finite wherever
        // tau > tau_0, so putting tau_0 below the worst determinant on the mesh
        // gives every element, inverted or not, a finite energy and a gradient
        // that pushes tau up. Set tau_0 through Options::untangleFloor, or
        // leave it to run() to track. Not scale-invariant once tau_0 < 0: an
        // element shrunk to a point scores zero. See Options::untangler.
        Untangle022 = 22,
        // (tau - 1)^2. Pure size, no barrier: it asks each element to have the
        // target area and does not care what shape it takes to get there.
        Size055 = 55,
        // 0.5 (tau + 1/tau) - 1. Pure size, barrier.
        Size056 = 56,
        // (1 - gamma) * Shape002 + gamma * Size056, the usual way to ask for
        // "mostly shape, a little size". gamma = 0 is Shape002 exactly.
        ShapeSizeCombo = 100,
        // ShapeSizeCombo with tau in both denominators replaced by
        //
        //     chi(tau, eps) = (tau + sqrt(tau^2 + eps^2)) / 2,
        //
        // i.e. (1 - g) (0.5 |T|^2 / chi - 1) + g (0.5 (tau^2 + 1) / chi - 1):
        // the regularisation of Garanzha and Kaporin, in the form of
        // Garanzha et al., "Foldover-free maps in 50 lines of code" (TOG 2021).
        // chi is positive for every tau and tends to tau as eps -> 0, so the
        // metric is finite on an inverted element and is ShapeSizeCombo again
        // once eps is gone. The untangler Options::untangler picks by default;
        // eps rides in the slot Untangle022 uses for tau_0 and g is
        // Options::untangleSizeWeight. Not a shape metric while eps > 0, and
        // that is the point of it -- see Options::untangler.
        UntangleRegularized = 101
    };

    // How the untangling phase gets an inverted element back out.
    enum Untangler : int {
        // Metric 22 with its barrier tau_0 re-placed below the worst
        // determinant every sweep. What this smoother did until 2026-09-24.
        UntangleShifted = 0,
        // Metric UntangleRegularized, eps lowered by the rule of Garanzha et al.
        // as the energy falls. **The default.**
        UntangleRegular = 1
    };

    // Where W comes from.
    enum Target : int {
        // Leave QuadMesh::targetJacobian as the caller set it. Falls back to
        // UniformSquare when it is empty.
        TargetKeep = 0,
        // W = h I for every element: the ideal axis-aligned square of side h,
        // h = Options::targetSize or the mesh's mean edge length. The ordinary
        // choice for shape optimisation.
        TargetUniformSquare = 1,
        // W = the element's own current Jacobian at its centre. Every element
        // starts at T = I at its centre, so the energy starts near zero and the
        // optimizer only removes the *variation* of shape within an element.
        // Size-consistent with wherever the boundary already sits, which makes
        // it the right target for testing the mobility layer in isolation: a
        // shape metric on it cannot introduce a size mismatch of its own.
        TargetCurrentShape = 2,
        // W = h_q I from Options::perQuadSize, for a size field coming from the
        // layout rather than from one global number.
        TargetPerQuadSize = 3
    };

    // Where mu is sampled inside an element.
    enum Quadrature : int {
        // 2x2 Gauss. Exact for the bilinear map's own quadratics and the
        // default: it sees the interior of the element, not only its corners.
        Gauss2x2 = 0,
        // The four corners, weight 1/4 each. The corner Jacobians are the ones
        // whose determinants the scaled-Jacobian report reads, so this drives
        // exactly the number an untangler is judged by.
        Corners = 1
    };

    struct Options {
        Metric metric = Shape002;
        // Weight on the size term of ShapeSizeCombo; ignored by every other
        // metric.
        double gamma = 0.2;

        // Minimise the integral of mu^exponent rather than of mu.
        //
        // At 1 the energy is an average, and an average can be lowered by
        // making most elements better while making the worst one worse -- which
        // is the wrong trade for a mesh, whose usable timestep is set by its
        // worst element and not by its mean. Raising the exponent weights the
        // bad elements by mu^(p-1) and turns the objective towards the maximum;
        // the limit p -> infinity is the worst-element problem itself.
        //
        // 2 by default, because 1 measurably loses. On the four MERIDIAN models
        // the CTest cases use, p = 1 leaves the worst scaled Jacobian at 0.512
        // on singlemat/geom012 -- *below* the 0.625 it started at, the average
        // paying for it with the one element nobody can afford -- while p = 2
        // brings it to 0.731. On the other three, where p = 1 does no harm,
        // p = 2 is still the better of the two on the worst element, and costs
        // about 0.0002 of the mean.
        //
        // Must be 1, or 2 or more: between the two the second derivative of
        // mu^p is singular where mu vanishes, which is exactly where a good
        // element sits. Values in (1, 2) are raised to 2 and reported.
        double exponent = 2.0;

        Target target = TargetUniformSquare;
        double targetSize = 0.0;             // h; <= 0 means the mean edge length
        std::vector<double> perQuadSize;     // for TargetPerQuadSize
        Quadrature quadrature = Gauss2x2;

        // Stop after this many sweeps, or when a sweep moves every node less
        // than moveTolerance times the mean edge length, or when the energy
        // falls by less than energyTolerance of the energy the mesh started
        // with. A negative tolerance turns its own test off, which is how a
        // caller asks for exactly maxSweeps sweeps.
        int maxSweeps = 200;
        double moveTolerance = 1e-6;
        double energyTolerance = 1e-10;

        // Newton safeguards. The local 2x2 Hessian is made positive definite by
        // lifting its smallest eigenvalue to hessianFloor times its largest,
        // and the step is capped at maxStepFraction of the shortest edge at the
        // node so a single Newton step cannot jump a node across its own
        // one-ring. Backtracking halves it up to maxLineSearch times, and a
        // step is accepted only when the node's local energy strictly falls.
        int maxLineSearch = 16;
        double maxStepFraction = 0.4;
        double hessianFloor = 1e-6;

        // Run an untangling phase first when the mesh starts with an inverted
        // element.
        bool untangle = true;
        int untangleMaxSweeps = 100;

        // Which untangler, and why the default changed.
        //
        // Metric 22 is 0.5 (|T|^2 - 2 tau) / (tau - tau_0), and with tau_0 < 0
        // it is not scale-invariant: shrink an element by k and the numerator
        // falls as k^2 while the denominator only falls to -tau_0, so an
        // element crushed to a point costs *nothing*. On an inverted element,
        // whose tau - tau_0 is tiny, shrinking is also by far the steepest way
        // down -- so the phase meant to open a fold closes it instead, dragging
        // the fold's neighbours in after it, and stops the moment the last
        // determinant creeps over zero. The main phase then starts on a mesh
        // with an element of det ~1e-8, whose barrier energy is ~1e20 and
        // whose nodes' Newton steps are capped at a fraction of their now
        // microscopic edges. Traced det_rocket at h = 0.03 is the case:
        // one fold at a singular block corner (scaled Jacobian -0.02) came out
        // as a star of ten nodes within 1e-4 of each other, worst aspect ratio
        // 297, worst scaled Jacobian 0.012.
        //
        // The regularised untangler has a size term, g (tau^2 + 1) / (2 chi),
        // that goes to g / eps as an element collapses, so crushing one is
        // expensive instead of free; and eps is lowered only as fast as the
        // energy falls, so the mesh leaves the phase as a minimiser of a
        // barrier shape + size energy rather than as the first mesh with every
        // determinant positive.
        Untangler untangler = UntangleRegular;
        // g in UntangleRegularized: the weight of the size term.
        // Garanzha et al. use 1/128 in 2-D.
        double untangleSizeWeight = 1.0 / 128.0;
        // UntangleShifted only: tau_0 as a fraction of the worst determinant
        // on the mesh, tau_0 = untangleFloor * min(tau) when that minimum is
        // negative. Above 1 so the barrier stays clear of the worst element.
        double untangleFloor = 1.5;

        // Number of OpenMP threads; 0 leaves it to the runtime.
        //
        // The mesh that comes out is bit-identical whatever this is: the sweep
        // is coloured, so nodes updated together never share an element, and
        // the two quantities that decide when to stop are summed in element
        // order rather than by a floating-point reduction. Setting it to 1 is
        // the quickest way to confirm that on a mesh of your own.
        int threads = 0;

        bool verbose = false;
    };

    struct Report {
        bool ran = false;
        bool converged = false;         // stopped on a tolerance, not on maxSweeps

        int sweeps = 0;                 // of the main phase
        int untangleSweeps = 0;
        int colors = 0;                 // in the Gauss-Seidel colouring
        int threads = 1;                // actually used
        bool openMP = false;            // compiled with OpenMP at all

        int movableNodes = 0;           // free + sliding
        int freeNodes = 0, slidingNodes = 0, fixedNodes = 0;

        // Element corners whose determinant no movable node can change, and
        // which are flat or inverted: never a reason to untangle, and under
        // Quadrature::Corners left out of the energy, their elements sampled
        // at their other corners (see frozenCorner).
        int frozenCorners = 0;

        // Feature curves the mesh was bound to, how many of them are fitted
        // splines rather than the runs' own polylines, and how many sliding
        // nodes ride one. slidingNodes - curveNodes is how many fell back to
        // the chord through their feature neighbours.
        int featureCurves = 0;
        int fittedCurves = 0;
        int curveNodes = 0;
        // The furthest any of those curves bows from the polyline of its own
        // nodes, and the furthest any bound node ended up from its curve. The
        // first is the geometry the smoother was allowed to recover; the second
        // must be zero to rounding, and is the invariant that says the
        // parameter and the position never came apart.
        double curveBow = 0.0;
        double curveDeviation = 0.0;

        // Energy per unit target area, so the number is comparable between
        // meshes and between targets. Both are measured with the main metric,
        // and are NaN when it is a barrier metric and the mesh is still tangled.
        double energyBefore = 0.0, energyAfter = 0.0;

        double minScaledJacobianBefore = 0.0, minScaledJacobianAfter = 0.0;
        double meanScaledJacobianBefore = 0.0, meanScaledJacobianAfter = 0.0;
        int invertedBefore = 0, invertedAfter = 0;
        double worstAspectBefore = 0.0, worstAspectAfter = 0.0;
        double minAreaBefore = 0.0, minAreaAfter = 0.0;

        // Total area of the mesh before and after.
        //
        // Under QuadMesh::Options::CurveChord these agree to rounding even when
        // the whole boundary has redistributed, because a node sliding along
        // the chord through its two feature neighbours moves parallel to the
        // base of the only triangle the polygon's area depends on it through --
        // see QuadMesh::projectStep -- and a difference is then a real defect
        // rather than an accumulation of small ones.
        //
        // On a curve it is the curve, not the chord, that is conserved: a node
        // moving along a discretised circle's interpolant sweeps the area
        // between the polygon and the curve as it goes, so the mesh's area
        // creeps towards the area the curve encloses. That is the correction
        // being asked for, not drift; `curveDeviation` is the invariant to read
        // instead.
        double areaBefore = 0.0, areaAfter = 0.0;

        // The largest distance any single node ended up from where it started,
        // and the total distance walked, both relative to the mean edge length.
        double maxDisplacement = 0.0;
        double meanDisplacement = 0.0;
        // Largest node move in the final sweep, relative to the mean edge
        // length. The convergence measure.
        double lastSweepMove = 0.0;

        double seconds = 0.0;

        std::vector<std::string> messages;
    };

    explicit TMOP(QuadMesh &mesh);
    TMOP(QuadMesh &mesh, const Options &opts);

    // Fill the targets, colour the mesh, untangle if asked, then sweep. Moves
    // the mesh's vertices in place and refreshes its slide tangents and quality
    // when it is done. Returns false only when there is nothing to do (no
    // elements, or no movable node) or when the mesh ends up worse than it
    // started, which it reports rather than reverting.
    bool run();

    const Report &getReport() const { return report; }
    const Options &getOptions() const { return options; }

    // ---- the pieces, exposed for testing and for a caller building its own
    // outer loop --------------------------------------------------------------

    // Fill QuadMesh::targetJacobian per Options::target and precompute W^{-1}
    // and det W. Called by run(); idempotent.
    void prepare();

    // Total TMOP energy of the mesh with the current metric, and of one
    // element. Infinite when a barrier metric meets an inverted element.
    double energy() const;
    double elementEnergy(int q) const;

    // dF/dx_v, the unprojected 2-D gradient of the total energy with respect to
    // node v -- only the elements around v contribute, so this is O(1). The
    // gradient is analytic; `nodeSystem` adds the exact local 2x2 Hessian.
    Point nodeGradient(int v) const;
    void nodeSystem(int v, Point &gradient, double hessian[4]) const;

    // Sum of the energies of the elements around v: what the line search tests.
    // `minTau` (optional) comes back with the smallest det T seen anywhere in
    // those elements, which is what says whether the patch is still valid.
    double nodeEnergy(int v, double *minTau = nullptr) const;

    // The displacement this node would take: the local Newton step, projected
    // onto what its node type allows, capped, and backtracked until the node's
    // local energy strictly falls without inverting anything around it. Zero
    // when nothing was accepted. Leaves the mesh where it found it.
    //
    // `newParam`, when given, comes back with the curve parameter the accepted
    // step lands on, for a node bound to a feature curve -- and with the
    // parameter the node already has for every other node, so a caller can
    // write it back unconditionally. Applying the displacement without it
    // leaves the node's parameter stale and its next step wrong; moveNode()
    // does both.
    Point nodeStep(int v, double *newParam = nullptr);
    // The same, applied. Returns the distance moved.
    double moveNode(int v);

    // Independent sets of the "shares an element" graph. Every node in one
    // bucket can be moved in parallel with every other node in it, because no
    // two of them touch the same quad -- so a coloured sweep and a serial
    // Gauss-Seidel sweep in the same node order compute the same answer.
    const std::vector<std::vector<int>> &colorBuckets() const { return buckets; }

    // One Gauss-Seidel sweep over every movable node, in colour order. Returns
    // the largest displacement any node took, in mesh units.
    double sweep();

    // The metric itself, for tests and for anyone wanting the numbers out:
    // mu(T) and, if `dmu` is non-null, its derivative with respect to T.
    // Returns false when the metric is undefined at this T (a barrier metric on
    // tau <= 0). `exponent` matches Options::exponent.
    static bool evalMetric(Metric m, const Jacobian2 &T, double gamma, double tau0,
                           double *mu, Jacobian2 *dmu = nullptr, double exponent = 1.0);

private:
    // mu written as f(s, t) with s = |T|^2 and t = det T, which every metric
    // here is; the local Hessian needs the five partials of f and nothing else.
    struct MetricPartials {
        bool valid = false;
        double f = 0.0;
        double fs = 0.0, ft = 0.0;
        double fss = 0.0, fst = 0.0, ftt = 0.0;
    };
    static MetricPartials partialsOf(Metric m, double s, double t, double gamma, double tau0);
    // The chain rule for mu -> mu^p, applied to a metric's partials.
    static MetricPartials raiseTo(const MetricPartials &base, double p);

    // Accumulate element q's energy and, when `grad`/`hess` are given, its
    // derivatives with respect to the node at corner `c` of it, over q's
    // quadrature points. `minTau` tracks the smallest det T seen. Returns false
    // when the metric was undefined anywhere in q, in which case the energy is
    // infinite and the derivatives are meaningless.
    bool accumulate(int q, int c, double *energyOut, double *minTau,
                    Point *grad, double *hess) const;

    void buildColoring();
    // Mark the frozen corners: see frozenCorner below.
    void findFrozenCorners();
    double minDeterminant() const;
    int countInverted() const;
    void snapshotQuality(bool before);
    double sweepPhase(Metric m, double tau0, double &maxMove);

    QuadMesh &mesh;
    Options options;
    Report report;

    // Per element, from the targets: W^{-1} and det W.
    std::vector<Jacobian2> targetInverse;
    std::vector<double> targetDet;

    // The metric the current phase is running, which is Untangle022 during
    // untangling and Options::metric afterwards, with its barrier position and
    // its exponent. The untangling phase always runs at exponent 1: its job is
    // to get every determinant positive, and weighting the elements that are
    // already fine would only slow that down.
    Metric activeMetric = Shape002;
    double activeTau0 = 0.0;      // tau_0 for 022, eps for UntangleRegularized
    double activeGamma = 0.0;     // gamma for the combos, the untangler's g
    double activeExponent = 1.0;

    std::vector<int> color;                  // per vertex, -1 for an immovable one
    std::vector<std::vector<int>> buckets;   // color -> its vertices, ascending

    // Per element, a bit per corner that is **frozen**: its determinant is
    // flat or inverted and cannot be changed by any node the smoother may move
    // -- to first order, which for straight features is exactly. The case is a
    // fixed node D whose two neighbours in the element can only slide parallel
    // to D's other edge, and in practice it is a block corner that lies on a
    // straight feature: D a junction or a pinned layout node, and the element's
    // two sides at D both running along the same line, so the corner is 180
    // degrees whatever anyone does. Such a sample is not an obstacle an
    // optimizer can remove, and left in the energy it is worse than useless:
    // its det T is pinned at zero, so every barrier and every regularised
    // untangler reads it as infinite or nearly so, and the only thing left
    // to lower is |T|^2 -- the two neighbours slide down the line onto D and
    // drag the rest of the strip with them. Under Quadrature::Corners a
    // frozen element is sampled at its other corners instead, so it still has
    // a barrier against folding any further; the Gauss points are not corners
    // and are left alone. Either way a mesh whose only non-positive corners
    // are frozen is not sent to the untangler at all.
    std::vector<std::uint8_t> frozenCorner;

    std::vector<Point> startPositions;
    double meanEdge = 1.0;
};

}  // namespace mesh

#endif  // __MESH_TMOP_HXX__
