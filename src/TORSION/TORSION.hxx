#ifndef __TORSION_HXX__
#define __TORSION_HXX__

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "MERIDIAN/Arrangement.hxx"
#include "MERIDIAN/ConeCut.hxx"
#include "MERIDIAN/ConeSingularities.hxx"
#include "MERIDIAN/Immersion.hxx"
#include "MERIDIAN/Interfaces.hxx"
#include "MERIDIAN/LayoutEnergy.hxx"
#include "MERIDIAN/MERIDIAN.hxx"
#include "MERIDIAN/QuadMesh.hxx"
#include "MERIDIAN/RicciFlow.hxx"
#include "MERIDIAN/Separatrices.hxx"
#include "MERIDIAN/SplineFit.hxx"
#include "MERIDIAN/SubdomainLabels.hxx"
#include "TORSION/FieldFrames.hxx"
#include "TORSION/FieldIntegration.hxx"
#include "TORSION/TutteEmbedding.hxx"
#include "mesh/Mesh.hxx"

class SIPG;

// Pipeline B of docs/cf_flow_pipeline.md: the quadrilateral layout of Shepherd,
// Gu and Hughes (2022) with the initial map psi_0 built by **integrating the
// SIPG cross field** instead of by unfolding a Ricci metric.
//
// The paper's Stages 1 to 10 are MERIDIAN's and they are used here verbatim.
// What changes is one thing, and everything in this class is either that thing
// or the plumbing around it:
//
//     MERIDIAN     Stage 3  discrete Ricci flow          -> a flat cone metric
//                  Stage 4  unfold it, triangle by triangle
//
//     TORSION      Stage 3F comb the field, read the matchings, audit the
//                           indices                          -> FieldFrames
//                  Stage 4F one constrained least-squares solve on Omega
//                                                            -> FieldIntegration
//                  Stage 4R untangle, if the solve inverted anything
//                                                     -> TutteEmbedding + Stage 6
//
// Both routes end at the same place -- an `Immersion` -- and every stage after
// it takes that by reference and cannot tell which route filled it. That is the
// whole design: Sec. 1 of the plan calls the Immersion "the hub", and the seam
// between the two pipelines is its constructor and nothing else.
//
// ### What the substitution buys and what it costs
//
//   Bought.  The transitions R_k are *exact from the matchings* rather than
//            least-squares fitted and snapped, so C2 holds by construction and
//            Immersion::Report::maxSnapError becomes a check on the branch
//            selection instead of a rounding step. Boundary and feature
//            alignment come free on the field, so E2 and E3 start small. The
//            initial map costs one sparse solve rather than up to a hundred
//            Newton steps on the flow. And sizing becomes a per-face h_t rather
//            than one global scaling of a Ricci metric.
//
//   Paid.    Local injectivity. A discrete cross field is generically
//            non-integrable, so the least-squares closest map inverts triangles
//            -- generically at the cones and in high-distortion pockets. Stage
//            6's barrier preserves Q1 and cannot restore it, so a map with
//            flips has to be untangled before it is a legal psi_0. Stage 4R
//            does that, out of parts that already existed: a Tutte embedding,
//            which is bijective by theorem, and then Stage 6 itself with E1
//            against the integrated map's own metric, a single mu on E4, and
//            everything else switched off.
//
// ### The one thing that is easy to miss (C4)
//
// E1 of Stage 6 measures J against a reference metric, and the choice of
// reference is the choice of what "undistorted" means. A field-integrated map
// has no flat metric, and if it ships psi_0 and leaves the reference Euclidean
// it inherits, exactly and in full, the failure LayoutEnergy.hxx measures on
// geom003: a planar model's Euclidean metric has no cones at all, so "J is a
// rotation" and Q2 contradict each other at every cone, E1 cannot win, and it
// settles by cramming the whole (pi/2) I of angle defect into one triangle of
// the cone's fan. The layout comes out topologically perfect and geometrically
// unusable.
//
// The plan's answer is LayoutEnergy::Reference::Field: the field's holonomy
// around a cone *is* (pi/2) I, so its frame should carry the cone structure the
// flat metric carries and the Euclidean lengths do not. It was built, and
// Sec. 11's risk register asks for it to be tested standalone before anything
// is built on it. It was, and **it does not work** --
//
//     TestTORSION --ref-test data/meshes/singlemat/geom003.obj
//     det J: Field vs Euclidean differ by at most 9.3e-16 relative
//
// -- because the frame is a scaled rotation, so composing the reference with it
// applies a rotation to the *right* of J, and the symmetric Dirichlet energy
// cannot see one. E1 against Reference::Field is E1 against the Euclidean
// geometry, exactly, and inherits its geom003 failure in full: 656
// quadrilaterals instead of 1001, four of them inverted. LayoutEnergy.hxx
// carries the algebra.
//
// What answers C4 is a reference whose *triangles are a different shape*, and
// the field route has one for free. The integration's seam constraints are
// exact rotations, so psi_0's image angle sum at each cone is already
// 2 pi - (pi/2) I to eight digits before Stage 6 starts: the metric psi_0
// induces on Omega is a flat cone metric with the cones the layout wants, got
// without a flow. That is Reference::Induced and it is this pipeline's default.
// Options::reference selects among the four so that the paragraph above stays a
// measurement rather than a claim.
class TORSION {
public:
    struct Options {
        // --- Stage 0 and 0b, identical to MERIDIAN's -----------------------
        double sipgGamma = 10.0;
        int sipgMaxSteps = 500;

        bool materialInterfaces = true;
        double interfaceKinkAngle = M_PI / 4.0;
        int interfaceLoopSplits = 4;
        bool alignFieldToInterfaces = true;
        bool cancelInterfaceDipoles = true;
        bool prescribeInterfaceCones = true;
        bool interfaceCorners = true;
        bool propagateInterfaceLabels = true;
        bool lagInterfaceScales = false;

        // --- Stage 1 and 2, identical to MERIDIAN's ------------------------
        bool autoRebalance = true;
        int minBoundaryIndex = -3;
        int maxBoundaryIndex = 1;
        bool coneCutsToBoundary = true;
        double coneCutInterfaceAvoidance = 1.0;

        // --- Stage 3F: the frames ------------------------------------------
        // h, the target edge length the frame is scaled by. The image of the
        // integration comes out 1/h times the model, so 1 is the neutral
        // choice; it is Pipeline B's analogue of Pipeline A's global Ricci
        // scaling and it is the natural place for a sizing field.
        double targetEdge = 1.0;
        // Per-face h_t. Empty for uniform.
        std::vector<double> sizing;
        // Comb from a second seed as well and check the two agree. The check
        // that matters is the loop test FieldFrames runs unconditionally; this
        // is the spot check Sec. 13's checklist asks for on top of it.
        bool checkSecondSeed = true;

        // --- Stage 4F: the integration -------------------------------------
        double integrationRegularisation = 1e-10;

        // --- Stage 4R: the untangling --------------------------------------
        // Run it when the integration inverted something. Off leaves psi_0 as
        // the solve returned it, which stops the pipeline at Stage 4's gate --
        // and that is the measurement Sec. 7.1's census is for.
        bool untangle = true;
        // Outer steps of the target-fitting continuation, and the weight on E4
        // relative to E1 at the first of them. mu grows by lambdaGrowth per
        // step exactly as the main continuation's penalties do.
        int untangleOuterSteps = 8;
        int untangleInnerIterations = 60;
        double untangleSeamWeight = 1e-2;

        // --- Stage 6 -------------------------------------------------------
        // C4: what "undistorted" means to E1. See the class comment for why
        // Field is not the answer the plan expected it to be and Induced is.
        //
        //   Induced    the metric psi_0 itself induces on Omega, which has the
        //              cone angles the seam constraints already forced into it.
        //              The default.
        //   Field      Sec. 4's construction, kept so that the measurement
        //              against it can be repeated. Identical to Euclidean.
        //   Euclidean  the model's own geometry, which on a planar input has no
        //              cones at all. The known failure.
        //   Ricci      Pipeline A's flat cone metric. Requires the flow, which
        //              this pipeline otherwise does not run, and is here as the
        //              known-good yardstick rather than as a working default.
        enum class Reference { Induced, Field, Euclidean, Ricci };
        Reference reference = Reference::Induced;

        bool runLayout = true;
        // Sec. 9's closing note: E2 and E3 start small on a field-integrated
        // map because the field is already aligned to dS and to the interfaces,
        // so over-penalising them early only fights E1; E4 starts higher
        // because the untangling introduced seam error. Both are expressed as
        // factors on the common lambda_init rather than as separate schedules,
        // so that turning them off recovers MERIDIAN's behaviour exactly.
        double lambdaInit = 1e-2;
        double lambdaGrowth = 10.0;
        double lambdaAlignmentFactor = 0.1;   // on lambda_2 and lambda_3
        double lambdaSeamFactor = 10.0;       // on lambda_4
        int outerSteps = 16;
        int innerIterations = 60;
        bool relabelBetweenSteps = true;

        bool seedTopoConstraints = true;
        double topoNearMiss = 0.15;
        bool seedSelfReturns = true;
        bool seedAllConnections = true;

        // --- Stages 7 to 10, identical to MERIDIAN's -----------------------
        bool runSeparatrices = true;
        double separatrixSnap = 1e-6;
        int separatrixMaxSteps = 50000;
        int separatrixSnapRings = 2;
        bool separatrixDetectCycles = true;

        bool runArrangement = true;
        double arrangementMerge = 1e-7;
        double arrangementCorner = 0.35;
        double arrangementCollapse = 1e-4;
        bool arrangementTrim = true;

        bool runSplines = true;
        int splineSegments = 3;
        int splineSamples = 8;

        bool runQuadMesh = true;
        double quadTargetEdge = 0.05;
        int quadMinIntervals = 1;
        int quadMaxIntervals = 0;
        bool quadUseSplines = true;
        bool quadInterfacesOnTracedArcs = true;
        int quadSmoothingPasses = 500;
        double quadSmoothingThreshold = 0.0;

        int repairPasses = 6;
        double repairGapLimit = 1e-3;
        double repairLambdaBoost = 10.0;
        int repairOuterSteps = 4;
        int repairMaxPerPass = 4;
        bool repairScoresArrangement = true;
        int repairPatience = 2;
    };

    struct Status {
        // --- Stage 0 / 0b / 1 / 2, as MERIDIAN reports them ----------------
        int mboSteps = 0;
        bool fieldConverged = false;
        bool fieldAlignedToInterfaces = false;
        int coneDipoleUnits = 0;

        int materials = 1;
        int interfaceEdges = 0;
        int interfaceBranches = 0;
        int interfaceNodes = 0;
        int interfaceIllPosedNodes = 0;
        double interfaceWorstSector = 0.0;
        int interfaceConesPrescribed = 0;
        int interfacePrescriptionShift = 0;
        int regions = 0;
        int regionsBalanced = 0;
        int regionQuartersMoved = 0;
        int regionCornersInserted = 0;
        int worstRegionDeficit = 0;

        bool conesAdmissible = false;
        int interiorCones = 0;
        int boundaryCones = 0;

        bool cutIsDisk = false;
        bool allConesOnBoundary = false;
        int cutInteriorJunctions = 0;
        int cutInterfaceEdges = 0;
        int cutInterfaceVertices = 0;
        int cutInterfaceNodes = 0;

        // --- Stage 3F: the frames and the audit ----------------------------
        bool framesRan = false;
        int combedFaces = 0;
        int unreachedFaces = 0;
        int combingDefects = 0;
        double maxFrameJump = 0.0;
        int indexMismatches = 0;
        int clusteredConePairs = 0;
        int highIndexCones = 0;
        int leftHandedFrames = 0;
        double maxMetricDisagreement = 0.0;
        // The spot check of Sec. 13: combing from a second seed differs from
        // the first by one global quarter turn and nothing else.
        bool secondSeedAgrees = true;
        bool framesValid = false;

        // --- Stage 4F: the integration -------------------------------------
        bool integrationRan = false;
        bool integrationSolved = false;
        double integrationSeamResidual = 0.0;
        int integrationFlippedFaces = 0;
        double integrationFlippedAreaFraction = 0.0;
        double integrationMinAreaRatio = 0.0;
        double integrationMaxFitResidual = 0.0;
        double integrationMeanFitResidual = 0.0;
        int integrationFlipsAtCones = 0;
        double integrationNearestFlipToCone = 0.0;

        // --- Stage 4R: the untangling --------------------------------------
        bool untangleRan = false;
        bool tutteValid = false;
        int untangleOuterSteps = 0;
        int untangleFlippedFaces = 0;     // after the pass
        double untangleSeamResidual = 0.0;

        // --- Stage 4: psi_0 ------------------------------------------------
        bool immersionValid = false;
        int immersionFlippedFaces = 0;
        int seamArcs = 0;
        // Q4 read as a check rather than as a rounding: how far the geometry of
        // psi_0 is from the quarter turn the matchings say each arc carries.
        double maxSnapError = 0.0;
        int frameKConflicts = 0;
        // The field metric that psi_0 does not realise -- the non-integrability
        // again, this time per edge of Omega.
        double maxMetricResidual = 0.0;

        // --- Stages 5 to 10, as MERIDIAN reports them ----------------------
        int boundaryEdgesU = 0, boundaryEdgesV = 0;
        int featureChains = 0;
        int topoPaths = 0;
        int topoSelfReturns = 0;
        int topoExtraPerPair = 0;

        int interfaceCorners = 0;
        int interfaceLabelsCorrected = 0;
        double interfaceResidual = 0.0;
        int interfaceCornerChanges = 0;
        bool interfacesAligned = false;

        bool layoutRan = false;
        bool layoutInjective = false;
        bool layoutConstrained = false;
        int outerStepsTaken = 0;
        bool layoutValid = false;
        double layoutConeAngleResidual = 0.0;
        double layoutRegularAngleResidual = 0.0;
        int layoutReferenceLeftHanded = 0;
        // Whether E1's reference actually carried the cones: the angle sum at
        // each cone of the reference metric itself, against 2 pi - (pi/2) I.
        // Reference::Euclidean has no cones at all and this is (pi/2)|I| there;
        // that is the whole of C4 in one number, taken before the continuation
        // rather than inferred from what it produced.
        double referenceConeResidual = 0.0;

        bool separatricesRan = false;
        int separatrices = 0;
        int separatricesToCone = 0;
        int separatricesToBoundary = 0;
        int separatricesUnresolved = 0;
        int separatricesNearMisses = 0;

        bool arrangementRan = false;
        int layoutNodes = 0;
        int layoutArcs = 0;
        int layoutPatches = 0;
        int layoutQuads = 0;
        double layoutCoverage = 0.0;
        bool arrangementValid = false;

        bool splinesRan = false;
        int splinePatches = 0;
        int splineControlPoints = 0;
        double splineMaxDeviation = 0.0;
        bool splinesWatertight = false;
        bool splinesValid = false;

        bool quadMeshRan = false;
        int meshVertices = 0;
        int meshQuads = 0;
        int meshChords = 0;
        int meshUnmeshedPatches = 0;
        double meshMinScaledJacobian = 0.0;
        bool meshConforming = false;
        bool meshValid = false;

        int repairPasses = 0;
        int repairConstraintsAdded = 0;
        bool q5Verified = false;

        std::vector<std::string> messages;
    };

    explicit TORSION(std::shared_ptr<Mesh> mesh);
    TORSION(std::shared_ptr<Mesh> mesh, const Options &opts);
    ~TORSION();

    // Stages 0 to 10 in order, stopping at the first failure the next stage
    // cannot survive. Returns Definition 2.1's verdict on Psi and nothing else.
    //
    // The two hard stops are the same as MERIDIAN's and for the same reasons:
    // an inadmissible cone set has no flat cone metric and no seamless map
    // either, so it stops before the cut is used; and a psi_0 with a flipped
    // face stops before Stage 6, because the barrier of Sec. 3.3 can preserve
    // Q1 and never repair it. The difference is that in this pipeline the
    // second one is *expected* to fire sometimes, which is why Stage 4R sits in
    // front of it.
    bool run();

    const SIPG& getField() const { return *field; }
    const Interfaces& getInterfaces() const { return *interfaces; }
    const ConeSingularities& getCones() const { return *cones; }
    const ConeCut& getCut() const { return *cutter; }
    const FieldFrames& getFrames() const { return *frames; }
    // Only built when Options::reference is Reference::Ricci, and then only to
    // supply E1's reference metric -- psi_0 still comes from the integration.
    const RicciFlow& getReferenceFlow() const { return *referenceFlow; }
    bool hasReferenceFlow() const { return referenceFlow != nullptr; }
    const FieldIntegration& getIntegration() const { return *integration; }
    const TutteEmbedding& getTutte() const { return *tutte; }
    const Immersion& getImmersion() const { return *immersion; }
    const SubdomainLabels& getLabels() const { return *labels; }
    const LayoutEnergy& getLayout() const { return *layout; }
    const Separatrices& getSeparatrices() const { return *separatrices; }
    const Arrangement& getArrangement() const { return *arrangement; }
    const SplineFit& getSplines() const { return *splines; }
    const QuadMesh& getQuadMesh() const { return *quads; }

    bool hasInterfaces() const { return interfaces != nullptr; }
    bool hasCones() const { return cones != nullptr; }
    bool hasCut() const { return cutter != nullptr; }
    bool hasFrames() const { return frames != nullptr; }
    bool hasIntegration() const { return integration != nullptr; }
    bool hasTutte() const { return tutte != nullptr; }
    bool hasImmersion() const { return immersion != nullptr; }
    bool hasLabels() const { return labels != nullptr; }
    bool hasLayout() const { return layout != nullptr; }
    bool hasSeparatrices() const { return separatrices != nullptr; }
    bool hasArrangement() const { return arrangement != nullptr; }
    bool hasSplines() const { return splines != nullptr; }
    bool hasQuadMesh() const { return quads != nullptr; }

    // psi_0 as the integration left it, before any untangling. Kept so that a
    // caller can see what the substitution actually produced rather than only
    // what survived the repair.
    const std::vector<Point>& getIntegratedMap() const { return integratedMap; }
    const Status& getStatus() const { return status; }

    // E1's reference metric for a field-integrated map: the lengths `psi`
    // induces, one per edge of the *input* mesh, which is the indexing
    // LayoutEnergy::Options::referenceLengths is read in.
    //
    // Static and public because it is the answer to C4 and not an
    // implementation detail: any caller that builds psi_0 stage by stage rather
    // than through run() -- the viewer does -- has to hand Stage 6 the same
    // reference, and a second copy of this walk is a second thing to keep in
    // step with it.
    static std::vector<double> inducedLengths(const Mesh &mesh, const ConeCut &cut,
                                              const std::vector<Point> &psi);

    const Mesh& getMesh() const { return *mesh; }
    std::shared_ptr<Mesh> getMeshPtr() const { return mesh; }

private:
    void runField();
    // Stages 0b, 0, 1 and 2, which are MERIDIAN's unchanged. Returns false when
    // the cone set is inadmissible or the cut is not a disk, which are the two
    // things every route downstream depends on.
    bool runFront();
    // Stage 4R. Takes the frames and returns an untangled map, or an empty
    // vector when it could not.
    std::vector<Point> untangle();

    std::shared_ptr<Mesh> mesh;
    Options options;
    Status status;

    std::unique_ptr<SIPG> field;
    std::unique_ptr<Interfaces> interfaces;
    std::unique_ptr<ConeSingularities> cones;
    std::unique_ptr<ConeCut> cutter;
    std::unique_ptr<FieldFrames> frames;
    std::unique_ptr<RicciFlow> referenceFlow;
    std::unique_ptr<Immersion> scaffold;
    std::unique_ptr<FieldIntegration> integration;
    std::unique_ptr<TutteEmbedding> tutte;
    std::unique_ptr<Immersion> immersion;
    std::unique_ptr<SubdomainLabels> labels;
    std::unique_ptr<LayoutEnergy> layout;
    std::unique_ptr<Separatrices> separatrices;
    std::unique_ptr<Arrangement> arrangement;
    std::unique_ptr<SplineFit> splines;
    std::unique_ptr<QuadMesh> quads;

    // The indices the cross field read, before Stage 1 moved any of them; the
    // set Sec. 5.1's audit is run against.
    std::vector<int> fieldIndex;
    std::vector<Point> integratedMap;
};

#endif // __TORSION_HXX__
