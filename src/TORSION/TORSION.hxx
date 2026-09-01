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
#include "TORSION/ConeMetric.hxx"
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
//     TORSION      Stage 3F comb the field, move its singularities onto the
//                           cone set, read the matchings, audit the indices,
//                           and read the axis of every chain of dS - G and of
//                           the interface network      -> FieldFrames
//                  Stage 4F one constrained least-squares solve on Omega, with
//                           the seam *and* those axes as equalities
//                                                            -> FieldIntegration
//                  Stage 4R untangle, if the solve inverted anything: locally
//                           first, and only then through a Tutte embedding
//                                        -> relaxToKernel + TutteEmbedding
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
//            alignment are not merely small on the field but *exact*: the field
//            is pinned tangent to dS and to every interface, so the axis of
//            each chain can be read off the frame and held as an equation in
//            the same solve the seam is held in (Sec. 6.4). The initial map
//            costs one sparse solve rather than up to a hundred Newton steps on
//            the flow -- and the flat cone metric E1 needs costs four or five
//            more, because on a planar model the curvature to be moved is only
//            what sits on the boundary and Newton on the vertex scaling gets it
//            in a handful of Laplacian solves (Sec. 4, ConeMetric). And sizing
//            becomes a per-face h_t rather than one global scaling of a Ricci
//            metric -- which is also where that metric's conformal factor
//            enters the frame.
//
//            Sec. 6.4 is what closes the gap that made this pipeline slow. A
//            Ricci map arrives at Stage 6 with Q3 true to 1e-15, because Stage
//            1 prescribes zero curvature at every non-cone boundary vertex and
//            the flow makes dS geodesic; a least-squares map arrived with a
//            boundary residual the size of the field's non-integrability, and
//            the continuation had to spend its whole budget getting there. On
//            singlemat/geom001 that was 8 outer steps and 236 inner iterations
//            against Pipeline A's 1 and 0 -- 12.1 s against 0.31 s. Held as
//            equalities it is 1 and 0, and 0.24 s.
//
//   Paid.    Local injectivity. A discrete cross field is generically
//            non-integrable, so the least-squares closest map inverts triangles
//            -- generically at the cones and in high-distortion pockets. Stage
//            6's barrier preserves Q1 and cannot restore it, so a map with
//            flips has to be untangled before it is a legal psi_0.
//
//            Rather less of it than it looks, once the frame carries the
//            conformal factor: a good part of what Report::maxFitResidual used
//            to read was the frame asking for an isometry, not the field being
//            non-integrable. See C4 below.
//
//            Stage 4R does that on a ladder. The usual tangle is single figures
//            of faces in one cone's one ring, and for that the map does not
//            need rebuilding: the signed area of each triangle at one vertex is
//            affine in that vertex, so "every triangle here is positive" is an
//            intersection of half planes and its interior is the kernel of the
//            one-ring polygon. Moving the vertex to the centre of it clears the
//            tangle and touches nothing else -- so psi_0 keeps the cone angles,
//            the seam and Sec. 6.4's alignment. Where that needs one of the
//            equalities let go, Sec. 6.5's projection puts it back by fitting
//            the untangled map's own Jacobian under it. Only what none of that
//            clears goes to the Tutte pass, which is bijective by theorem and
//            keeps nothing at all.
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
// What answers C4 is a reference whose *triangles are a different shape*. The
// obvious candidate is the metric psi_0 induces on Omega: the seam constraints
// are exact rotations, so its angle sum at every cone is 2 pi - (pi/2) I to
// eight digits before Stage 6 starts, and Status::referenceConeResidual duly
// reads 1e-10 on it. **That is not enough, and the corpus says so.** The cone
// angles are only the part of the reference that lives at the cones; the rest
// of it is the shape of every other triangle, and E1 against psi_0's own
// lengths says "reproduce psi_0" -- shear, non-integrability and all. So the
// cone *fans* stay however the least-squares fit distributed them, the three
// rays of a +1 cone leave at 200, 60 and 100 degrees instead of 120 each, the
// patch at that cone has a reflex corner, and Stage 10 meshes an element the
// scaled Jacobian reads as -0.39 -- on multimat/geom001, whose layout passes
// Q1 to Q5 with residuals of 1e-7. Q1-Q5 do not see it, which is exactly why
// Reference::Induced survived as a default for as long as it did.
//
// The reference has to be a *flat cone metric*: conformal to the model away
// from the cones, and carrying the prescribed angle defect at them. Pipeline A
// gets one from the flow. The field route does not need the flow to have one,
// because the model is planar: the interior curvature is already zero, the
// whole of the flow's work is moving the boundary's turning into the cones,
// and a conformal change of a triangulation is a vertex scaling whose exact
// Jacobian dK/du is the cotangent Laplacian. So the metric is a Newton solve
// on the model's own triangles -- one linear solve for the first order and
// three or four more for the rest. That is ConeMetric, it is
// Reference::Cone, and it is this pipeline's default. On the corpus it is the
// same metric --ref ricci computes, and it costs a handful of Laplacian solves
// rather than the flow's circle packing, flow triangulation and edge flips.
//
// The same u fixes something upstream of the reference. The frame was
//
//     J*_t = (1/h) R(-theta_t),
//
// which asks the map to be an isometry of the model *everywhere*, and a map
// with cones cannot be one: Sec. 3's own Schwarz-Christoffel integrand has
// modulus prod |z - p_k|^(-I_k/4), which is exp(u) and is not constant. The
// target the integration was being handed was therefore not integrable even
// where the field is, and Report::maxFitResidual was reading a
// non-integrability that was partly the frame's. With the frame scaled --
// Options::conformalSizing, on by default --
//
//     singlemat/geom003     6 inverted faces at Stage 4F  ->  0
//     multimat/geom002      3                             ->  0, and Sec. 6.4's
//                                                            strain 0.27 -> 0.11
//
// Status::referenceConeResidual is the reference's own cone angles, read before
// the continuation is allowed to use it. Options::reference selects among the
// five so that the paragraphs above stay measurements rather than claims.
class TORSION {
public:
    struct Options {
        // --- Stage 0 and 0b, identical to MERIDIAN's -----------------------
        double sipgGamma = 10.0;
        int sipgMaxSteps = 500;

        bool materialInterfaces = true;
        double interfaceKinkAngle = M_PI / 4.0;
        int interfaceLoopSplits = 4;
        // Stages 0c and 11, identical to MERIDIAN's: excise every circular
        // inclusion before the field runs and put its O-grid back from a
        // template after Stage 10. See DiskTemplate for the reasoning and
        // MERIDIAN::Options::diskTemplates for the short version.
        bool diskTemplates = false;
        double diskCoreSquareness = 0.55;
        int diskRingDepth = 0;
        int diskSmoothingPasses = 300;
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

        // Sec. 6.4: read the axis of every chain of dS - G and of every
        // interface branch off the frame, and hold it *exactly* in the
        // integration rather than asking for it with lambda_2 and lambda_3.
        //
        // This is what closes the gap the head of this file calls structural.
        // A Ricci map arrives at Stage 6 with Q3 already true -- Stage 1
        // prescribes zero curvature at every non-cone boundary vertex, so the
        // flow makes dS geodesic and psi_R's boundary is a rectilinear polygon
        // -- and the continuation then exits at its first level on the easy
        // models. A least-squares map arrives with a boundary residual the size
        // of the field's non-integrability, and the continuation has to spend
        // its whole budget driving that down. Holding the axis instead costs
        // one equation per boundary edge in a solve that already has the seam
        // as equalities, and it takes the staircase out at the same time: the
        // axis is decided per chain and the chains end at the cones.
        bool alignInIntegration = true;
        // If the aligned solve inverts more faces than the free one would, keep
        // the free one. One extra factorisation, taken only when the aligned
        // solve inverted something, and it is what makes the alignment a
        // strictly-better-or-equal change on the flip census rather than a
        // trade.
        bool alignmentFallback = true;

        // --- Stage 4F: the integration -------------------------------------
        double integrationRegularisation = 1e-10;

        // --- Stage 4R: the untangling --------------------------------------
        // Run it when the integration inverted something. Off leaves psi_0 as
        // the solve returned it, which stops the pipeline at Stage 4's gate --
        // and that is the measurement Sec. 7.1's census is for.
        bool untangle = true;
        // Sec. 7.2a: try the local repair before Sec. 7.2's global one. It
        // moves the interior vertices around an inverted face to the centre of
        // their one-ring kernels and touches nothing on dOmega, so a map it
        // clears keeps its cone angles, its seam and Sec. 6.4's alignment
        // exactly -- all three of which the Tutte pass throws away and Stage 6
        // then has to rebuild. See relaxToKernel() in TORSION.cxx.
        bool localUntangle = true;
        int localUntangleSweeps = 32;
        // Let a seam vertex move at rungs 0 and 1 of that ladder, with its
        // partner moving under R_k so that Q4 is unchanged. See
        // TORSION::seamPairing().
        //
        // **Off, and measured.** The construction is right and the extra
        // freedom is real, and on this corpus it buys nothing: the tangles that
        // survive rung 0 are not tangles one paired displacement can undo --
        // rungs 0 and 1 clear the same one model with it as without -- and the
        // trajectory it takes at rung 2 leaves the seam far enough out that
        // Sec. 6.5's projection stops being able to put it back, which costs
        // multimat/geom002 its layout. Kept behind a flag because the argument
        // for it does not go away with the measurement: on a model whose tangle
        // straddles the cut rather than sitting beside it, this is the freedom
        // that pass needs.
        bool pairSeamInUntangle = false;
        // Outer steps of the target-fitting continuation, and the weight on E4
        // relative to E1 at the first of them. mu grows by lambdaGrowth per
        // step exactly as the main continuation's penalties do.
        int untangleOuterSteps = 8;
        int untangleInnerIterations = 60;
        double untangleSeamWeight = 1e-2;
        // The weight on E2 and E3 in that fit, relative to E1, on the same
        // schedule as E4's.
        //
        // **Off, and measured.** The argument for it is that the pass starts
        // from a Tutte embedding of a circle, which has none of Sec. 6.4's
        // alignment, and that carrying the alignment across the pass is
        // cheaper than letting Stage 6 rediscover it. The measurement does not
        // support it: on geom004 the cone at the end of the cut comes back
        // 2 pi out -- the fit spends its budget dragging the boundary onto the
        // axes and settles for a cone that has changed valence -- and Q2 then
        // cannot be recovered by anything downstream. Sec. 6.5's re-projection
        // does the same job afterwards, as a projection onto an equality rather
        // than as one more term competing with E1 for the same budget, and it
        // does not have to fight the barrier to do it.
        double untangleAlignWeight = 0.0;
        // Sec. 6.5: after the untangling, fit the untangled map's own Jacobian
        // under the seam and alignment equalities, so that Q3 and Q4 come back
        // exactly. Kept only when it inverts nothing. See
        // FieldIntegration::Options::targetJacobian.
        bool reprojectAfterUntangle = true;

        // --- Stage 6 -------------------------------------------------------
        // C4: what "undistorted" means to E1.
        //
        //   Cone       the flat cone metric of Stage 1's cone set, built from
        //              the model's own triangles by ConeMetric's Newton solve.
        //              The default, and the answer to C4.
        //   Induced    the metric psi_0 itself induces on Omega. It has the
        //              cone *angles* -- the seam equalities forced them in --
        //              and none of the shape: E1 against it says "reproduce
        //              psi_0", including the shear the non-integrability of
        //              the field left in it, so the cone fans stay however the
        //              integration distributed them. Measured on the corpus it
        //              is much nearer Euclidean than to a flat cone metric.
        //   Field      Sec. 4's construction, kept so that the measurement
        //              against it can be repeated. Identical to Euclidean: the
        //              frame is a rotation and E1 is invariant under one.
        //   Euclidean  the model's own geometry, which on a planar input has no
        //              cones at all. The known failure.
        //   Ricci      Pipeline A's flat cone metric, from the flow. The same
        //              object Cone builds, computed the expensive way, and the
        //              yardstick Cone is checked against.
        enum class Reference { Cone, Induced, Field, Euclidean, Ricci };
        Reference reference = Reference::Cone;

        // Scale the frame by the cone metric's conformal factor, so that the
        // integration's target is J*_t = (exp(u_t)/h) R(-theta_t) rather than
        // (1/h) R(-theta_t).
        //
        // The unscaled frame asks the map to be an isometry of the model
        // everywhere at once, and a map with cones cannot be one -- the
        // Schwarz-Christoffel integrand of Sec. 3 has modulus
        // prod |z - p_k|^(-I_k/4), which is exp(u) and is not constant. So the
        // target the integration was being handed was not integrable *even
        // where the field is*, and Report::maxFitResidual was reading a
        // non-integrability that was partly the frame's own. Scaling it costs
        // one multiplication per face; u is already computed for the reference.
        bool conformalSizing = true;

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
        // See SplineFit::Options::fitBoundaryArcs: off, dS and the material
        // interface network are carried as the polylines the mesh has them as
        // rather than approximated.
        bool fitBoundaryArcs = false;
        bool fitInterfaceArcs = false;

        bool runQuadMesh = true;
        double quadTargetEdge = 0.05;
        int quadMinIntervals = 1;
        int quadMaxIntervals = 0;
        bool quadUseSplines = true;
        bool quadFeaturesOnTracedArcs = true;
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

        // --- Stage 3F: the alignment (Sec. 6.4) ----------------------------
        int alignmentChains = 0;
        int alignedBoundaryEdges = 0;
        int alignedInterfaceEdges = 0;
        int alignmentOverrides = 0;
        int alignmentClosedChains = 0;
        double maxAlignmentResidual = 0.0;

        // --- Stage 4F: the integration -------------------------------------
        bool integrationRan = false;
        bool integrationSolved = false;
        int integrationAlignedEdges = 0;
        double integrationAlignResidual = 0.0;
        double integrationAlignStrain = 0.0;
        bool alignmentWasDropped = false;
        double integrationSeamResidual = 0.0;
        int integrationFlippedFaces = 0;
        double integrationFlippedAreaFraction = 0.0;
        double integrationMinAreaRatio = 0.0;
        double integrationMaxFitResidual = 0.0;
        double integrationMeanFitResidual = 0.0;
        int integrationFlipsAtCones = 0;
        double integrationNearestFlipToCone = 0.0;

        // --- Stage 4R: the untangling --------------------------------------
        bool localUntangleRan = false;
        int localUntangleFlippedFaces = 0;   // after Sec. 7.2a, before Sec. 7.2
        // Which rung of Sec. 7.2a's ladder cleared it: 0 held the seam and the
        // alignment throughout, 1 let the alignment go, 2 let the seam go too.
        int localUntangleLevel = -1;
        bool untangleRan = false;
        bool tutteValid = false;
        int untangleOuterSteps = 0;
        int untangleFlippedFaces = 0;     // after the pass
        double untangleSeamResidual = 0.0;
        // Sec. 6.5's re-projection: whether it ran and whether it was kept.
        bool reprojectionRan = false;
        bool reprojectionKept = false;
        int reprojectionFlippedFaces = 0;

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

        // The cone metric of Sec. 4, whether it is E1's reference or not: it is
        // built whenever conformalSizing is on, and its residual is the
        // statement that the flat cone metric of this cone set exists and was
        // found on the model's own triangulation.
        bool coneMetricRan = false;
        bool coneMetricSolved = false;
        bool coneMetricConverged = false;
        int coneMetricNewtonSteps = 0;
        double coneMetricInitialError = 0.0;
        double coneMetricLinearError = 0.0;
        double coneMetricFinalError = 0.0;
        double coneMetricConeResidual = 0.0;
        double coneMetricMinScale = 1.0;
        double coneMetricMaxScale = 1.0;
        int coneMetricNonRealisable = 0;
        bool conformalSizingUsed = false;

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
        // Stages 0c and 11. See MERIDIAN::Status for what each one counts.
        int diskInclusions = 0;
        int diskTrianglesExcised = 0;
        int meshOddLoops = 0;
        int meshParityChordsMoved = 0;
        int meshOddLoopsLeft = 0;
        bool diskTemplatesRan = false;
        int diskTemplatesFilled = 0;
        int diskTemplateBlocks = 0;
        int diskTemplateQuads = 0;
        int diskTemplatesRefused = 0;
        double diskTemplateMinScaledJacobian = 0.0;
        int mergedVertices = 0;
        int mergedQuads = 0;
        double mergedMinScaledJacobian = 0.0;
        bool diskTemplatesValid = false;
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
    const ConeMetric& getConeMetric() const { return *coneMetric; }
    bool hasConeMetric() const { return coneMetric != nullptr; }
    const FieldIntegration& getIntegration() const { return *integration; }
    const TutteEmbedding& getTutte() const { return *tutte; }
    const Immersion& getImmersion() const { return *immersion; }
    const SubdomainLabels& getLabels() const { return *labels; }
    const LayoutEnergy& getLayout() const { return *layout; }
    const Separatrices& getSeparatrices() const { return *separatrices; }
    const Arrangement& getArrangement() const { return *arrangement; }
    const SplineFit& getSplines() const { return *splines; }
    const QuadMesh& getQuadMesh() const { return *quads; }
    const DiskTemplate& getDiskTemplate() const { return *diskFill; }
    bool hasDiskTemplate() const { return diskFill != nullptr; }
    const std::vector<DiskTemplate::Inclusion>& getInclusions() const { return inclusions; }

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

    // The mesh the layout was computed on. With Options::diskTemplates this is
    // the *excised* mesh; getInputMesh() is what came in. See MERIDIAN.
    const Mesh& getMesh() const { return *mesh; }
    const Mesh& getInputMesh() const { return inputMesh ? *inputMesh : *mesh; }
    std::shared_ptr<Mesh> getMeshPtr() const { return mesh; }

private:
    void runField();
    // Stage 0c, MERIDIAN's unchanged: find the circular inclusions and take
    // them out of `mesh`. See DiskTemplate.
    void exciseDisks();
    // Stages 0b, 0, 1 and 2, which are MERIDIAN's unchanged. Returns false when
    // the cone set is inadmissible or the cut is not a disk, which are the two
    // things every route downstream depends on.
    bool runFront();
    // Which coordinates each vertex of Omega may move in without undoing an
    // equality Stage 4F imposed. See TORSION.cxx.
    std::vector<unsigned char> vertexFreedom(int level) const;
    // The two children of each seam vertex of Omega and the quarter turn
    // between their displacements, so that Sec. 7.2a can move a tangle sitting
    // on the seam without letting go of Q4. See TORSION.cxx.
    void seamPairing(std::vector<int> &mate, std::vector<int> &turn) const;
    // Stage 4R. Takes the frames and returns an untangled map, or an empty
    // vector when it could not.
    std::vector<Point> untangle();

    std::shared_ptr<Mesh> mesh;
    // What run() was handed, kept only when Stage 0c replaced it.
    std::shared_ptr<Mesh> inputMesh;
    Options options;
    Status status;

    std::unique_ptr<SIPG> field;
    std::unique_ptr<Interfaces> interfaces;
    std::unique_ptr<ConeSingularities> cones;
    std::unique_ptr<ConeCut> cutter;
    std::unique_ptr<FieldFrames> frames;
    std::unique_ptr<RicciFlow> referenceFlow;
    std::unique_ptr<ConeMetric> coneMetric;
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
    std::unique_ptr<DiskTemplate> diskFill;

    // The circular inclusions Stage 0c took out. Empty unless
    // Options::diskTemplates.
    std::vector<DiskTemplate::Inclusion> inclusions;

    // The indices the cross field read, before Stage 1 moved any of them; the
    // set Sec. 5.1's audit is run against.
    std::vector<int> fieldIndex;
    std::vector<Point> integratedMap;
    // The alignment Stage 4F actually used, after its fallback chose among the
    // three. Sec. 6.5's re-projection has to hold the same one.
    std::vector<int> usedAxis;
};

#endif // __TORSION_HXX__
