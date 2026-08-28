#ifndef __LAYOUTENERGY_HXX__
#define __LAYOUTENERGY_HXX__

#include <array>
#include <string>
#include <vector>

// eigen includes
#include <Eigen/Sparse>

#include "MERIDIAN/Immersion.hxx"
#include "MERIDIAN/SubdomainLabels.hxx"
#include "mesh/Mesh.hxx"

// Stage 6 of Shepherd, Gu and Hughes (2022), Section 3.3: minimisation against
// the quadrilateral layout-inducing energies. (docs/shepherd2022.pdf)
//
// This is the stage that actually produces the layout. Everything before it
// exists to hand it a starting map that is already valid where it cannot repair
// itself; everything after it reads what comes out.
//
//     F   = { phi : Omega -> R^2, phi in C^0, phi locally invertible }   Eq. (12)
//     Psi = argmin_{phi in F} sum_{j=1}^{5} lambda_j E_j(phi)            Eq. (13)
//
// The unknowns are (u, v) at every vertex of the *cut* mesh, piecewise linear
// per triangle, so the Jacobian J_t is constant on each triangle. One vertex is
// pinned to kill the translation null space; nothing else is constrained
// directly, which is the point -- Definition 2.1's five properties are asked
// for by energies, not by equations.
//
// ### The five energies
//
//   E1  symmetric Dirichlet, Eq. (14): integral of ||J||^2 + ||J^-1||^2.
//       Preserves Q1. The second term diverges as det J -> 0+, so it is a
//       *barrier*: it stops a triangle inverting and cannot undo an inversion.
//       lambda_1 stays at 1 for the whole continuation.
//
//   E2  boundary alignment, Eq. (15): sum over Gamma_u of l_e^-1 (u_i - u_j)^2,
//       and the same in v over Gamma_v. Buys Q3.
//
//   E3  feature alignment, Eq. (16): identical, on the feature chains. This is
//       what makes a feature of the input a layout edge of the output.
//
//   E4  seam consistency, Eq. (17): for each paired child edge of an arc in
//       Gamma_Hol_k, || (phi+(j) - phi+(i)) - R_k^-1 (phi-(j) - phi-(i)) ||^2.
//       Written on tangents, so the translation part of the transition drops
//       out and only the rotation is held. Buys Q2 and Q4.
//
//   E5  connectivity, Eqs. (18) and (19): for each path of Gamma_topo, the
//       square of the total advance of the constrained coordinate summed over
//       its subcurves. Buys Q5. Remark 3.1 is the reason it is not optional:
//       the same bracket midsurface gives 1898 sliver patches without it and 73
//       usable ones with it.
//
//   E6  interface sectors -- not in the paper, and needed as soon as the
//       features are the interfaces of a multi-material domain rather than the
//       trim curves of one. For each sector of the interface network, with
//       tangents a and b at the node that bound it and q the quarter turns
//       Stage 0b measured between them,
//
//           E6 = sum  (l_a + l_b)/2  || b/l_b  -  R_q (a/l_a) ||^2
//
//       written on tangents in exactly the way E4 is, so that the translation
//       drops out and only the turn is held.
//
//       ### Why E3 is not enough
//
//       E3 holds each interface on a coordinate line, which is a statement
//       about one curve at a time. Everything that makes a network of curves a
//       *layout* lives between them and E3 cannot see it. Three interfaces meet
//       at a triple point: E3 says each is an isoline, and is equally satisfied
//       whether the three leave at 90, 90, 180 degrees in the image -- the
//       sectors the model actually has -- or all three along the same axis,
//       which is degenerate and which the barrier will happily settle for
//       because it is cheaper. An interface reaching dS: E3 says the interface
//       is an isoline and E2 says dS is one, and nothing at all says they are
//       different isolines, so the interface is free to fold onto the boundary.
//       A corner where an interface turns: E3 splits it into two chains and
//       says each is an isoline, which two collinear chains satisfy exactly.
//
//       In every one of those the missing statement is the same, and it is the
//       one the geometry already fixed: how many quarter turns lie between two
//       curves at the point where they meet. That is what E6 asserts, and it is
//       the reason a multi-material layout comes out with its interfaces on it
//       rather than beside it.
//
//       ### Why it is quadratic
//
//       The obvious form of "these two tangents are q quarter turns apart" is a
//       condition on their angles, which is not quadratic in (u, v), and the
//       obvious repair -- cross(R_q a, b) = 0 -- is bilinear, which is worse.
//       Writing it as b = R_q a instead is linear in the unknowns and its
//       square is one more term of exactly the kind E2..E5 already are, at the
//       price of also asserting |b| = |a|: the map stretches the two rays
//       equally at that one vertex. That is a real extra assumption and it is a
//       mild one -- it is local conformality at the node and nowhere else, and
//       E4 has been asserting the same thing across every seam edge since Stage
//       4 -- and it is what makes the term cost nothing to add to the existing
//       proxy.
//
// E5 is real-valued rather than integer-valued, and that is the paper's central
// numerical claim -- the mixed-integer optimisation that [48,49,51,53,60] need
// is avoided entirely, because "the two endpoints lie on a common isoline" is
// an equation in R, not a condition on Z.
//
// ### Why a penalty method and not a Lagrangian one
//
// Sec. 3.3 is explicit that the constraints have to hold *exactly*: a map that
// satisfies them weakly is a map for which Definition 2.1 is false, and a
// layout is then simply not defined. That rules out a Nitsche formulation,
// which the paper says so in as many words, and leaves penalty continuation
// with lambda_2..lambda_5 -> infinity. The schedule is 1e-2 up by a factor of
// ten per outer step, eight to twelve steps.
//
// ### The inner solve
//
// Sec. 3.3 offers composite majorization [92] or a SLIM-style quadratic proxy
// [91]. This is the second. At the current iterate every triangle's J is
// factorised as U Sigma V^T; the rotation part R = U V^T is the nearest
// rotation, and a weight matrix W = U diag(s) U^T is chosen with
//
//     s_i^2 = D_sigma_i / (sigma_i - 1)
//
// so that the proxy (1/2) A ||W (J - R)||_F^2 has *exactly* the gradient of
// A D(sigma) at the current point. That last property is what lets E1's proxy
// be added to E2..E5 -- which are already quadratic and need no proxy at all --
// without biasing the minimum: the quadratic model's gradient is the true
// gradient of the whole objective, so its minimiser is a genuine descent
// direction for the whole objective, and the line search is a real Armijo test
// on the real energy.
//
// The step is then capped so that no triangle inverts:
//
//     s_max = min over triangles of the smallest positive root of det J_t(s) = 0
//     s     = 0.8 min(1, s_max), then backtrack on the energy
//
// det J_t(s) is a quadratic in s for a linear step direction, so the root is
// closed form and the cap costs nothing.
//
// ### The two escapes from a stalled continuation
//
// Both are Sec. 3.3's own. The reference metric of J can be switched between
// the Euclidean geometry the surface came with and the Ricci metric psi_R was
// built from; the two have different local minima and alternating between them
// at successive lambda levels is how the paper gets out of one. And the
// Gamma_u / Gamma_v labelling can be re-read from the current map, since the
// labelling that was right for psi_R need not be right once E2 has pulled the
// boundary around.
//
// Both are applied only when the previous outer step failed to move the worst
// constraint, not at every level. They are escapes, and taking them
// unconditionally does harm: switching the reference at every level keeps
// restarting the distortion term from a different notion of undistorted, and
// the continuation never settles.
//
// ### Why the reference switch is off here, and the relabel is not
//
// The paper's S is a midsurface embedded in R^3, and there "the Euclidean
// geometry the surface came with" means its induced metric -- a curved metric,
// whose curvature the Ricci metric redistributes into the cones but does not
// invent. The two references then really are two nearby notions of undistorted
// and alternating between them is a genuine escape.
//
// This codebase is planar. The domain arrives as a flat triangulation, so its
// Euclidean metric has an angle sum of exactly 2pi at every interior vertex:
// it has **no cones at all**. The Ricci metric has the cone angles Stage 1
// prescribed, 2pi - (pi/2) I. So on a planar input the two references do not
// differ by a choice of local minimum, they differ by the whole cone structure,
// and "J is a rotation" and Q2 are contradictory statements about the same
// vertex. E1 cannot win that argument -- lambda_2..5 reach 1e5 to 1e7 and the
// constraints are met to machine precision -- so it settles it the only way it
// can: the map becomes an isometry of the *model* everywhere except the cones,
// and the entire (pi/2) I of angle defect is crammed into the cone's one ring.
//
// That is not a subtle degradation. On geom003, a disk, psi_R hands over cone
// fans of six triangles at 44.1 to 45.9 degrees apiece -- the conformal answer.
// One switch to the Euclidean reference turns the same fan into
// 26.6, 12.1, 58.6, 148.6, 11.7, 12.5, still summing to the required 270. The
// separatrices leave along the image axes, so two of the three leave inside the
// single fat triangle, the layout's centre patch acquires four reflex corners
// on the model, its Coons blend folds, and Stage 10 inverts elements it cannot
// smooth out. The layout is topologically perfect throughout -- the O-grid a
// disk should get -- and geometrically unusable.
//
// Rescaling the Euclidean reference to the Ricci one's total area does not help
// (measured: identical result), which is the check that this is about the cones
// and not about the 140x global scale between the two.
//
// So `alternateReference` defaults to false. The relabel escape is unaffected
// by any of this and stays on: it re-reads Gamma_u / Gamma_v from the current
// map and says nothing about what "undistorted" means.
//
// ### One departure from a uniform lambda schedule
//
// lambda_2, lambda_3 and lambda_5 start at lambda_init because their
// constraints are badly violated by psi_R and the continuation's job is to
// reach them. lambda_4 does not: psi_R satisfies Q4 exactly, so E4's job is to
// *hold* it while the other three drag the map around, and starting it low lets
// the early levels trade Q4 away for progress elsewhere. Q4 is not free to
// trade. The angular form of its residual is ||d+ - R_k d-|| / l_e, which the
// l_e^-1 weighting of Eq. (17) makes largest on the short seam edges -- and
// those are the ones at the cones, where the same quantity is the cone angle
// sum of Q2. Left at lambda_init, Q2 comes out of the continuation off by 1e-2
// radians on the harder models; started level with lambda_1, by 1e-6.
//
// ### What "converged" is tested on
//
// The continuation stops when *all five* properties hold, not when the four
// energies' residuals are small. The distinction matters because Q2 lags the
// rest by the ratio of the model's extent to its shortest boundary edge at a
// cone: on a model whose cone set has two cones a hundredth of the model apart,
// Q3 arrives at 1e-7 after seven levels and Q2 is still at 1e-2, and a run that
// stopped on the length residuals alone would report a layout it had not
// produced. Where Q2 refuses to come down after the budget is spent, the cause
// is almost always upstream -- clustered cones out of Sec. 3.1's automatic
// placement -- and run() says so by name rather than reporting a number.
class LayoutEnergy {
public:
    // Which metric J is measured against. Sec. 3.3's final paragraph.
    enum class Reference { Ricci = 0, Euclidean = 1 };

    struct Options {
        double lambda1 = 1.0;        // fixed; it is the injectivity barrier
        double lambdaInit = 1e-2;    // lambda_2 .. lambda_5
        double lambdaGrowth = 10.0;
        // Sec. 3.3 quotes 8 to 12 outer steps. That figure belongs to its
        // composite-majorization inner solve; the SLIM-style proxy used here
        // converges more slowly near the barrier, so the schedule is run
        // further rather than the tolerance being loosened. Each extra level
        // is cheap, since an outer step whose inner solve has nothing left to
        // do exits after one or two iterations.
        int outerSteps = 16;
        int innerIterations = 60;

        // Relative to the extent of the image: a step shorter than this ends
        // the inner solve, a constraint residual under this counts as met.
        double innerTolerance = 1e-9;
        double constraintTolerance = 1e-6;
        // Q2, in radians. Tight enough that Sec. 3.4's direction extraction at
        // a cone still lands in the right triangle, loose enough not to trip on
        // the ordinary residual of a converged continuation.
        double angleTolerance = 1e-3;

        double stepSafety = 0.8;     // Sec. 3.3's 0.8 * s_max
        int maxBacktracks = 50;

        Reference reference = Reference::Ricci;

        // Off: on a planar input the Euclidean reference has no cones, so it
        // fights Q2 rather than offering a second local minimum. See "Why the
        // reference switch is off here" above.
        bool alternateReference = false;
        bool relabel = true;

        // Re-read the lengths E6 divides its two tangents by from the map
        // between outer steps, instead of leaving them at the flat metric's.
        //
        // The argument for it is clean and the measurement does not support it,
        // so it is off. Written on flat lengths, E6 asserts the sector *and*
        // that the map stretches the two rays of it equally -- local
        // conformality at the node -- and that second half is a real extra
        // demand where the interfaces meet obliquely. Lagging the lengths
        // removes it and leaves the angle alone, and it does what it says: on
        // geom013, whose groove faces meet dS at 58 and 122 degrees, the worst
        // sector goes from 0.41 rad out to 0.035.
        //
        // The layouts are worse. Over data/meshes/multimat the unmeshed area
        // rises on geom010 (20.2% to 27.1%), geom012 (60.1% to 64.2%) and
        // geom013 (77.2% to 91.1%), and falls on geom007 alone (78.5% to
        // 66.2%). The conformality E6 was accidentally asking for turns out to
        // be regularising the continuation near the nodes, and giving it up
        // costs more than the sector residual it buys.
        bool lagInterfaceScales = false;

        int pinnedVertex = 0;
    };

    struct Report {
        int outerSteps = 0;
        int resumes = 0;            // rounds of Sec. 3.3's repair
        int innerIterations = 0;
        int lineSearchFailures = 0;
        int factorisationFallbacks = 0;
        int relabels = 0;
        int referenceSwitches = 0;

        double lambdaFinal[5] = {0.0, 0.0, 0.0, 0.0, 0.0};   // lambda_2 .. lambda_6
        double energyStart = 0.0, energyEnd = 0.0;
        // The five energies at the end, unweighted.
        double e1 = 0.0, e2 = 0.0, e3 = 0.0, e4 = 0.0, e5 = 0.0;

        double e6 = 0.0;
        // One residual per property of Definition 2.1, each in the units the
        // property is stated in. The three lengths are divided by the extent of
        // the image, so the tolerance is a relative one.
        double maxBoundaryResidual = 0.0;   // Q3
        double maxFeatureResidual = 0.0;    // features
        double maxSeamResidual = 0.0;       // Q4
        double maxTopoResidual = 0.0;       // Q5
        // E6: the worst sector of the interface network, as the length by which
        // the two tangents bounding it miss being a quarter turn apart, over the
        // extent of the image. Reported alongside the count of sectors whose
        // turn came out as a *different* whole number of quarters from the one
        // the model has -- the same distinction Q2 draws between a cone that
        // drifted and a cone that changed valence, and the same reason for
        // drawing it: an interface corner that turned the wrong way is a layout
        // with the interface somewhere else, not a layout with a slack residual.
        double maxInterfaceResidual = 0.0;
        int interfaceCornerChanges = 0;
        int interfaceCorners = 0;
        double minDetJ = 0.0;               // Q1, against the current reference
        double maxConeAngleResidual = 0.0;  // Q2, against the prescribed angle
        double maxRegularAngleResidual = 0.0;
        // Cones whose angle sum settled on a *different* multiple of pi/2 from
        // the one Stage 1 prescribed: the layout still has a cone there, of the
        // wrong valence.
        int coneValenceChanges = 0;
        int worstConeAngleVertex = -1;   // where maxConeAngleResidual was measured

        // The starting residuals, so that the continuation's progress is
        // visible rather than only its endpoint.
        double initialBoundaryResidual = 0.0;
        double initialSeamResidual = 0.0;
        double initialTopoResidual = 0.0;

        int invertedTriangles = 0;

        bool injective = false;       // Q1
        bool constraintsMet = false;  // Q3, Q4, Q5 and the features
        bool anglesHeld = false;      // Q2
        // All five properties of Definition 2.1.
        //
        // Q2 is checked separately from the constraint residuals rather than
        // being taken as a consequence of them, and the difference is not
        // academic. Q3, Q4 and Q5 are measured as the energies control them --
        // a length, relative to the extent of the image. The angle at a vertex
        // is that length divided by the *local* edge, and the l_e^-1 weighting
        // the paper's Eqs. (15) and (17) inherit from their integrals means the
        // shortest edges end up with the largest angular error. So a map can
        // have every constraint residual at 1e-7 and still have a cone whose
        // angle is a tenth of a radian out, if two cones landed a thousandth of
        // the model apart -- and that map is not a layout, whatever the
        // residuals say. Immersion::Report::clusteredConePairs is where the
        // cause shows up.
        bool valid = false;

        std::vector<std::string> messages;
    };

    // `labels` is held by reference and relabelled in place between outer
    // steps, so it has to outlive this object.
    // Two forms rather than a default argument, for the same reason as
    // SubdomainLabels: Options has default member initialisers.
    LayoutEnergy(const Immersion &immersion, SubdomainLabels &labels);
    LayoutEnergy(const Immersion &immersion, SubdomainLabels &labels,
                 const Options &opts);

    // The penalty continuation. Returns true when every constraint came under
    // tolerance with the map still locally injective -- that is, when what came
    // back is a quadrilateral layout in the sense of Definition 2.1.
    bool run();

    // Sec. 3.3's remedy for a layout that came back not being one, and Sec. 4's
    // "on failure" for an arrangement that would not close: continue the
    // continuation from the map as it now stands, at penalties `lambdaBoost`
    // times the ones it reached, after the caller has added the connectivity
    // constraints the separatrices of that map named.
    //
    // From the current phi, not from psi_R. The distinction is the paper's and
    // it is not a matter of speed: psi_R satisfies Q1, Q2 and Q4 and nothing
    // else, so restarting there throws away every level of lambda already paid
    // for and re-enters the same local minimum by the same road. Continuing
    // from phi enters it with the new constraint already switched on, which is
    // the only thing that has changed.
    //
    // Constraint terms are rebuilt first, so any Gamma_topo path added since
    // the last call is picked up. The labelling is *not* re-read: Gamma_u and
    // Gamma_v are re-read only when the continuation stalls, which is a
    // different question from this one.
    bool resume(int outerSteps, double lambdaBoost);

    // Re-read the subdomains into the constraint terms of Eqs. (15) to (19),
    // keeping the map. Called by resume(); public because a caller that has
    // added constraints without wanting to solve again still needs the
    // residuals in Report to mean something.
    void rebuildConstraints();

    // Psi, one planar point per vertex of Omega. Before run() this is psi_R.
    const std::vector<Point>& getUV() const { return uv; }

    // The map, saved and put back. What the repair loop of MERIDIAN::run()
    // uses to make itself monotone: a round of Sec. 3.3's remedy adds
    // constraints and re-minimises, and while that is nearly always an
    // improvement it is a *different* minimisation, not a continuation of the
    // same one, so it can land somewhere worse. Being able to go back is what
    // lets it be tried at all.
    std::vector<Point> saveMap() const { return uv; }
    void loadMap(const std::vector<Point> &m);

    // Per-triangle det J against the current reference, for drawing the
    // distortion and for checking Q1 from outside.
    std::vector<double> determinants() const;

    // Largest relative disagreement between the assembled gradient and a
    // central difference of the energy, over `samples` degrees of freedom.
    // The gradient is the one part of this stage that a wrong sign or a
    // transposed weight matrix leaves *looking* fine -- the line search still
    // finds a decrease, just a worse one -- so this is what the test exercises.
    double checkGradient(int samples = 24, double h = 1e-6) const;

    const Report& getReport() const { return report; }
    const Immersion& getImmersion() const { return *imm; }
    const SubdomainLabels& getLabels() const { return *labels; }

    bool writeOBJ(const std::string &filename) const;

private:
    // A single squared term w (a.x - b)^2, accumulated into the gradient and
    // the quadratic model. Everything except E1 is built out of these, and E1's
    // proxy is too.
    struct Term {
        std::vector<std::pair<int, double>> coeffs;  // (dof, coefficient)
        double gweight = 0.0;    // the geometric weight, l_e^-1 or 1; no lambda
        int which = 0;           // which of E2..E6 this term belongs to
    };

    void buildReference(Reference ref);
    void jacobian(const std::vector<double> &x, int t, double J[4]) const;

    double energy(const std::vector<double> &x, double *e1 = nullptr, double *e2 = nullptr,
                  double *e3 = nullptr, double *e4 = nullptr, double *e5 = nullptr,
                  double *e6 = nullptr) const;

    // The quadratic model at x: gradient of the true objective, and the
    // Hessian of the proxy. Returns false if the map has already inverted.
    //
    // `H` is written *into an existing pattern* -- buildHessianPattern()'s --
    // rather than assembled from a triplet list. Both halves of that matter and
    // neither is a micro-optimisation. The pattern of the proxy Hessian depends
    // on the mesh and on the constraint set and on nothing else: not on x, not
    // on lambda, not on the SLIM weights. Rebuilding it from ~570k triplets at
    // every inner iteration -- sorting them, collapsing the duplicates, and
    // then handing the result to a factorisation that recomputes its fill
    // ordering from scratch -- was a little under half of the whole pipeline's
    // running time on this corpus, all of it spent rediscovering the same
    // answer. Written this way the inner iteration scatters into a fixed array
    // and the ordering is computed once per constraint set.
    bool model(const std::vector<double> &x, std::vector<double> &grad,
               Eigen::SparseMatrix<double> &H) const;

    // The sparsity of that proxy, and the slot in its value array that each
    // local contribution lands in. Rebuilt whenever the constraint set changes
    // -- which is buildConstraintTerms(), and therefore also relabel() and
    // rebuildConstraints() -- and never otherwise.
    //
    // The diagonal is always in the pattern even where nothing writes to it, so
    // that the Tikhonov shift of innerSolve()'s fallback can be added in place
    // without changing the pattern and forcing a re-analysis.
    void buildHessianPattern();

    // The constraint terms. E2 to E5 are exactly quadratic in the unknowns
    // with a zero target, so each is a single squared linear form and the
    // coefficient lists only have to be rebuilt when the labelling changes --
    // never when lambda or x moves.
    void buildConstraintTerms();
    // Rewrite E6's coefficients from the map as it now stands, in place. The
    // dofs each term touches do not change -- only the lengths the two tangents
    // are divided by -- so the proxy's pattern and its symbolic factorisation
    // both survive, which is what makes this affordable once per outer step.
    // See SubdomainLabels::refreshInterfaceScales for why it is done at all.
    void refreshInterfaceTerms();
    // Where E6's terms start in cterms, and how many there are.
    size_t e6Start = 0, e6Count = 0;
    static double residualOf(const Term &t, const std::vector<double> &x) {
        double r = 0.0;
        for (const auto &c : t.coeffs) r += c.second * x[c.first];
        return r;
    }
    bool constraintsUnder(double tol) const;

    // The three verdicts of Definition 2.1 and their conjunction, from the
    // residuals measure() has just taken. Called at the end of a continuation
    // and again whenever the map or the constraint set is changed underneath
    // one -- loadMap() and rebuildConstraints() both do that, and a Report
    // still describing the map that was rolled back from is worse than no
    // Report at all.
    void finaliseVerdict();

    // The outer loop of Sec. 3.3, shared by run() and resume(). The penalties
    // are whatever the caller left in lambda[]; the schedule multiplies them
    // between steps, never before the first, so that the map returned is the
    // one minimised at the largest lambda rather than one level below it.
    bool continuation(int outerSteps);

    double maxStep(const std::vector<double> &x, const std::vector<double> &d) const;
    bool innerSolve(int outer);
    void measure(bool initial);

    // The diagonal of the image *as it currently stands*, which is the length
    // every residual below is divided by and the scale every step length is
    // taken relative to.
    //
    // This has to be re-read rather than taken once from psi_R. E1 measures
    // distortion against a reference triangle, and switching the reference from
    // the Ricci metric to the surface's own Euclidean geometry -- which Sec.
    // 3.3's continuation does, and which run() does here -- rescales the whole
    // map by the ratio of the two metrics' scales. On the meshes in this corpus
    // that ratio is around 40. A residual divided by psi_R's extent after such
    // a switch is understated by exactly that factor, and since constraintsUnder()
    // is what stops the continuation, the penalty loop declares Q3/Q4/Q5 met
    // roughly 40x before they are. The visible symptom is Stage 7: separatrices
    // that pass a few times 1e-6 of the image from a cone instead of landing on
    // it, so nothing snaps and curves that should join two singularities run on
    // to the step cap instead.
    double imageExtent() const;
    void updateExtent() { extent = imageExtent(); }

    static int dofU(int v) { return 2 * v; }
    static int dofV(int v) { return 2 * v + 1; }

    const Immersion *imm = nullptr;
    SubdomainLabels *labels = nullptr;
    Options options;

    const Mesh *cm = nullptr;
    int nV = 0;

    std::vector<Point> uv;
    std::vector<double> x;

    // Per triangle: M = [r1-r0, r2-r0]^-1 stored row-major, and the reference
    // area. The reference layout is always built positively oriented, so
    // det J > 0 is the same statement as "the image triangle is not flipped".
    std::vector<std::array<double, 4>> refM;
    std::vector<double> refArea;
    std::vector<Term> cterms;
    Reference currentReference = Reference::Ricci;

    // --- the proxy Hessian, assembled in place ----------------------------
    // reducedDof maps a dof to its row of the reduced system, or -1 for the two
    // pinned ones. Hproxy carries the pattern; triSlot and termSlot say where
    // each local entry goes in its value array, in the same order model()
    // produces them (a triangle contributes a 6x6 block over
    // (u,v) x (its three vertices); a constraint term a k x k block over its
    // own coefficient list). A -1 slot is an entry against a pinned dof.
    std::vector<int> reducedDof;
    int nReduced = 0;
    Eigen::SparseMatrix<double> Hproxy;
    std::vector<int> triSlot;        // 36 per triangle
    std::vector<int> termSlotStart;  // cterms.size() + 1 prefix sums
    std::vector<int> termSlot;       // k*k per term
    std::vector<int> diagSlot;       // nReduced, for the Tikhonov shift

    // The symbolic half of the factorisation, which depends on the pattern and
    // not on the values, so it survives every inner iteration at every level of
    // lambda. `analysed` is cleared with the pattern.
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
    bool analysed = false;

    double lambda[7] = {1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};   // lambda[1..6]
    // Diagonal of the image, the scale residuals are read against. Kept current
    // by updateExtent() -- see the note there for why a stale one is not a
    // cosmetic error in the reporting but an early stop of the continuation.
    double extent = 1.0;

    Report report;
};

#endif // __LAYOUTENERGY_HXX__
