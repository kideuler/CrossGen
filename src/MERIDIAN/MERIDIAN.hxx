#ifndef __MERIDIAN_HXX__
#define __MERIDIAN_HXX__

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "MERIDIAN/ConeCut.hxx"
#include "MERIDIAN/ConeSingularities.hxx"
#include "MERIDIAN/Immersion.hxx"
#include "MERIDIAN/LayoutEnergy.hxx"
#include "MERIDIAN/RicciFlow.hxx"
#include "MERIDIAN/Separatrices.hxx"
#include "MERIDIAN/SubdomainLabels.hxx"
#include "mesh/Mesh.hxx"

class SIPG;

// Stages 1 to 6 of
//   Shepherd, Gu and Hughes, "Isogeometric model reconstruction of open shells
//   via Ricci flow and quadrilateral layout-inducing energies", Engineering
//   Structures 252 (2022) 113602.   (docs/shepherd2022.pdf)
//
// The paper's contribution is to stop treating a quadrilateral layout as a
// combinatorial local-to-global problem and to characterise it instead as one
// continuous map Psi : S - G -> R^2 satisfying five properties Q1-Q5
// (Definition 2.1). That turns the whole thing into a variational problem, and
// it needs a starting map that is already valid in the properties the
// variational stage cannot repair -- Q1 in particular, since the symmetric
// Dirichlet barrier of Sec. 3.3 can prevent a flipped triangle but never undo
// one. Ricci flow is where that starting map comes from, and getting to it
// takes three stages:
//
//   Stage 1  Sec. 3.1    Cone singularities. Cross field in, integer indices
//                        out, checked against the discrete Gauss-Bonnet
//                        condition sum I(v) = 4 chi(S) of Eq. (4).
//                        -> ConeSingularities
//
//   Stage 2  Sec. 3.2.2  Cutting graph G, and Omega = S - G. Every void opened,
//                        every interior cone dragged to the boundary, the
//                        result a topological disk.
//                        -> ConeCut
//
//   Stage 3  Sec. 3.2.1  Discrete surface Ricci flow. Replace the metric
//                        inherited from the plane by a flat metric whose
//                        curvature is zero everywhere except (pi/2) I(v) at the
//                        cones.
//                        -> RicciFlow
//
//   Stage 4  Sec. 3.2.2  The metric immersion psi_R : Omega -> R^2. Unfold the
//                        flat metric into the plane, one triangle at a time,
//                        and fit the transition of every arc of G, so that Q4
//                        is not merely satisfied but known.
//                        -> Immersion
//
//   Stage 5  Sec. 3.3    Subdomain labelling. Which boundary edges want u held
//                        constant and which want v, the same for the feature
//                        chains, and the connectivity constraints Gamma_topo
//                        that E5 is written against.
//                        -> SubdomainLabels
//
//   Stage 6  Sec. 3.3    Minimise sum_j lambda_j E_j by penalty continuation,
//                        starting from psi_R, until Q1 to Q5 all hold. This is
//                        where the layout is actually produced.
//                        -> LayoutEnergy
//
//   Stage 7  Sec. 4      The separatrices of Psi: the integral curves that
//                        leave each cone along its grid directions, marched
//                        triangle by triangle and continued across the cutting
//                        graph by the transitions of Q4, until each terminates
//                        at a cone or leaves through dS. This is Q5 read back
//                        off the map as curves rather than as a residual, and
//                        it is what Stages 8 to 10 partition and fit.
//                        -> Separatrices
//
// The Gauss-Bonnet check between Stages 1 and 3 is not a formality. Newton's
// system in the flow is Delta du = Kbar - K with Delta a Laplacian, whose
// kernel is the constants; the residual has to be orthogonal to that kernel for
// the system to be consistent at all, and sum(Kbar) = sum(K) = 2 pi chi is
// exactly that condition. An inadmissible cone set does not converge slowly,
// it has no solution. run() therefore stops there rather than handing an
// inconsistent system to the solver.
//
// The other ordering constraint is the one between Stages 4 and 6. E1 of
// Sec. 3.3 is the symmetric Dirichlet energy, whose ||J^-1||_F^2 term diverges
// as det J -> 0 and is therefore a *barrier*: it stops a triangle inverting but
// has no way to undo an inversion that is already there. Stage 6 has to be
// handed a map that is already locally injective, and Stage 4 is where that map
// comes from -- which is why run() reports Stage 4's flipped-face count rather
// than only its own convergence.
//
// What comes out at the end is Psi, one planar point per vertex of the cut
// disk, satisfying Q1 to Q5 to whatever tolerance the continuation reached, and
// the separatrices traced on it. Stages 8 to 10 -- the arrangement of those
// curves into a layout, the spline fit and the analysis handoff -- read them
// and are not implemented here.
class MERIDIAN {
public:
    struct Options {
        // Stage 0: the cross field the cones are read off. The defaults match
        // TestSIPG.
        double sipgGamma = 10.0;
        int sipgMaxSteps = 500;

        // Stage 1
        bool autoRebalance = true;   // restore Eq. (4) by moving boundary cones
        int minBoundaryIndex = -3;
        int maxBoundaryIndex = 1;

        // Stage 3
        double ricciTolerance = 1e-8;
        int ricciMaxIterations = 100;
        bool delaunayFlips = true;

        // Stages 5 and 6. The defaults are the paper's own figures: lambda_1
        // fixed at 1 because it is the injectivity barrier and raising it only
        // fights the constraints, lambda_2..lambda_5 starting at 1e-2 and
        // multiplied by ten per outer step, eight to twelve outer steps.
        bool runLayout = true;
        double lambdaInit = 1e-2;
        double lambdaGrowth = 10.0;
        int outerSteps = 16;
        int innerIterations = 60;
        bool alternateReference = true;
        bool relabelBetweenSteps = true;

        // Gamma_topo. Remark 3.1: without the connectivity constraints a
        // topologically valid layout still exists but is made of slivers --
        // Fig. 13 is 1898 patches against 73 -- so they are seeded
        // automatically from near-misses of the separatrices on psi_R.
        bool seedTopoConstraints = true;
        double topoNearMiss = 0.15;  // as a fraction of the mean cone spacing
        // Q5's "possibly identical" case -- the curve of Fig. 9, out of a cone
        // and back to it -- and the second and later curves joining a pair of
        // cones already joined by one. Both are constraints in their own right
        // and both are easy to drop; see SubdomainLabels::Options.
        bool seedSelfReturns = true;
        bool seedAllConnections = true;

        // Stage 7. The snap tolerance is Sec. 3.4's, as a fraction of the image
        // extent; the step cap is what turns a curve that never terminates into
        // a reported failure of Q5 rather than a hang.
        bool runSeparatrices = true;
        double separatrixSnap = 1e-6;
        int separatrixMaxSteps = 50000;
        int separatrixSnapRings = 2;
        bool separatrixDetectCycles = true;

        // Sec. 3.3's remedy, and the "on failure" of Sec. 4's arrangement, run
        // as a loop rather than left to the reader:
        //
        //     trace the separatrices of Psi
        //     for each that ended at neither a cone nor dS, and for each that
        //         slipped past a cone on its way out through dS,
        //         add the Gamma_topo constraint it names
        //     raise lambda_5 and re-run Stage 6 *from the current phi*
        //     trace again
        //
        // This is what the seeding at Stage 5 cannot do by itself. Gamma_topo
        // is seeded on psi_R, because psi_R is the only map there is when
        // Stage 5 runs, and Stage 6 then moves it a long way -- that is the
        // whole point of Stage 6. A connection that is plain on the finished
        // layout need never have been visible on the map it started from.
        //
        // repairGapLimit is how near a curve has to have passed its cone, as a
        // fraction of the image extent, before the pair counts as one that was
        // meant to be joined. Too large and E5 is asked for a set of integral
        // curves that no map has, which does not converge slowly: the barrier
        // fights a penalty it cannot satisfy and det J is driven to zero.
        int repairPasses = 6;
        double repairGapLimit = 1e-3;
        double repairLambdaBoost = 10.0;
        int repairOuterSteps = 4;
        // How many constraints one pass may add, nearest first; 0 for all of
        // them. Four is what the corpus in data/meshes settles on: it is the
        // largest value at which every model that a repair helps at all is
        // still helped, and small enough that the pass stays satisfiable on the
        // models whose cones Stage 1 left clustered. See
        // SubdomainLabels::adoptCurves for why a cap helps at all.
        int repairMaxPerPass = 4;
    };

    struct Status {
        int mboSteps = 0;
        bool fieldConverged = false;

        bool conesAdmissible = false;
        int interiorCones = 0;
        int boundaryCones = 0;

        bool cutIsDisk = false;
        bool allConesOnBoundary = false;

        bool ricciConverged = false;

        // Every stage did what Definition 2.1 needs of it, so Stage 4 has a
        // flat cone metric on a disk to immerse.
        bool readyForImmersion = false;

        // Stage 4
        bool immersionValid = false;      // everything placed, isometric, no folds
        int immersionFlippedFaces = 0;
        int seamArcs = 0;

        // Stage 5
        int boundaryEdgesU = 0, boundaryEdgesV = 0;
        int featureChains = 0;
        int topoPaths = 0;
        int topoSelfReturns = 0;    // paths back to the cone they left, Fig. 9
        int topoExtraPerPair = 0;   // second and later curves of one cone pair

        // Stage 6
        bool layoutRan = false;
        bool layoutInjective = false;     // det J > 0 on every triangle, Q1
        bool layoutConstrained = false;   // every constraint residual under tolerance
        int outerStepsTaken = 0;

        // Q1 to Q5 all hold, so what came back is a quadrilateral layout in the
        // sense of Definition 2.1.
        bool layoutValid = false;

        // Stage 7
        bool separatricesRan = false;
        int separatrices = 0;
        int separatricesToCone = 0;
        int separatricesToBoundary = 0;
        int separatricesUnresolved = 0;   // capped, cycling, stuck or degenerate
        // Curves that terminated at neither a cone nor dS, or slipped past one
        // on the way out, and were near enough for the pair to be worth a
        // constraint. Zero is what a finished layout looks like.
        int separatricesNearMisses = 0;

        // The repair of Sec. 3.3, as it actually went.
        int repairPasses = 0;
        int repairConstraintsAdded = 0;
        // Every separatrix ended the way Q5 allows and every cone emitted the
        // number its index prescribes. Stage 6 asserts Q5 through the residual
        // of E5, which is a statement about the Gamma_topo paths it was given;
        // this is the same property read off every curve there is.
        bool q5Verified = false;

        std::vector<std::string> messages;
    };

    // Sec. 3.3's repair, as options and as a result. Separate from Options
    // because traceAndRepair() is usable on its own -- the viewer drives Stages
    // 4 to 7 itself, one keypress at a time, and runs this rather than a second
    // copy of it.
    struct RepairOptions {
        int passes = 6;
        int maxPerPass = 4;
        double gapLimit = 1e-3;
        double lambdaBoost = 10.0;
        int outerSteps = 4;
    };

    struct RepairResult {
        // The tracing of the map that came back, which is the map the caller's
        // LayoutEnergy now holds.
        std::unique_ptr<Separatrices> separatrices;
        bool traced = false;
        int passesTaken = 0;
        int constraintsAdded = 0;
        int rolledBack = 0;
    };

    // Stage 7, then Sec. 3.3's remedy applied in a loop until the separatrices
    // all terminate the way Q5 allows or there is nothing further to add. Both
    // stages are taken by reference and left holding the result. Messages are
    // reported through `log`, already prefixed with the stage they came from.
    static RepairResult traceAndRepair(SubdomainLabels &labels, LayoutEnergy &layout,
                                       const Separatrices::Options &trace,
                                       const RepairOptions &opts,
                                       const std::function<void(const std::string &)> &log);

    explicit MERIDIAN(std::shared_ptr<Mesh> mesh);
    MERIDIAN(std::shared_ptr<Mesh> mesh, const Options &opts);
    ~MERIDIAN();

    // Runs the three stages in order and stops at the first one that fails in a
    // way the next stage cannot survive: an inadmissible cone set stops the
    // pipeline, a cut that came out as something other than a disk does not
    // (Ricci flow does not use it, so the metric is still worth computing and
    // the failure is still worth reporting).
    bool run();

    const SIPG& getField() const { return *field; }
    const ConeSingularities& getCones() const { return *cones; }
    const ConeCut& getCut() const { return *cutter; }
    const RicciFlow& getRicci() const { return *ricci; }
    const Immersion& getImmersion() const { return *immersion; }
    const SubdomainLabels& getLabels() const { return *labels; }
    const LayoutEnergy& getLayout() const { return *layout; }
    const Separatrices& getSeparatrices() const { return *separatrices; }

    // Null until the stage that builds them has run.
    bool hasCones() const { return cones != nullptr; }
    bool hasCut() const { return cutter != nullptr; }
    bool hasRicci() const { return ricci != nullptr; }
    bool hasImmersion() const { return immersion != nullptr; }
    bool hasLabels() const { return labels != nullptr; }
    bool hasLayout() const { return layout != nullptr; }
    bool hasSeparatrices() const { return separatrices != nullptr; }

    const Status& getStatus() const { return status; }

    const Mesh& getMesh() const { return *mesh; }
    std::shared_ptr<Mesh> getMeshPtr() const { return mesh; }

private:
    void runField();

    std::shared_ptr<Mesh> mesh;
    Options options;
    Status status;

    std::unique_ptr<SIPG> field;
    std::unique_ptr<ConeSingularities> cones;
    std::unique_ptr<ConeCut> cutter;
    std::unique_ptr<RicciFlow> ricci;
    std::unique_ptr<Immersion> immersion;
    std::unique_ptr<SubdomainLabels> labels;
    std::unique_ptr<LayoutEnergy> layout;
    std::unique_ptr<Separatrices> separatrices;
};

#endif // __MERIDIAN_HXX__
