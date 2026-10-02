#ifndef __TORSION_HXX__
#define __TORSION_HXX__

#include <array>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

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
// The full definition, not a forward declaration: Options names
// DualMBO::PenaltyWeight, which is a nested type and cannot be forward declared.
#include "dualmbo/DualMBO.hxx"
#include "mesh/Mesh.hxx"

class DualMBO;
class MaterialLayout;

// Pipeline B of docs/cf_flow_pipeline.md: the quadrilateral layout of Shepherd,
// Gu and Hughes (2022) with the initial map psi_0 built by **integrating the
// DualMBO cross field** instead of by unfolding a Ricci metric.
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
//
// ### What a multi-material model asks of Stages 3F to 4R
//
// On a single material the cut touches dS and nothing else. On the multi-
// material corpus a cone inside an enclosed region has no route to dS that does
// not cross the interface around it, concrete's 74 interior cones cross 43
// interface vertices between them, and every place the pipeline had assumed
// one chart, one vertex or one attempt was wrong somewhere. In the order the
// data flows, and each behind an option so the difference stays a measurement:
//
//   Stage 3F  FieldFrames reads an interface branch the cut crosses in the
//             chart each piece of it is combed in (alignAcrossSeams). Voting
//             such a branch as one chart held the far side of every crossing on
//             the perpendicular axis: on concrete, 29 of 58 branches, 109 edges
//             overridden, a strain of 1.83, and the whole alignment dropped.
//
//   Stage 3F  It moves the field's singularities onto the cone set sector by
//             sector and holds every pinned face through the re-smoothing
//             (reconcileSectors). The vertex-level transport let dS absorb
//             anything, walked paths through interface vertices, and rotated
//             the faces the alignment reads its axes off by up to a half turn;
//             every model that needed a transport failed but one.
//
//   Stage 4F  The full alignment is kept whenever Sec. 7.2a can clear what it
//             inverted, and otherwise only the chains the tangle touches are
//             released, a few rounds at most (alignmentReleaseRounds), before
//             anything wholesale. Choosing on the raw flip count dropped the
//             alignment on 20 of the 31 models for tangles of 2 to 30 faces.
//
//   Stage 4R  Sec. 7.2c opens the folds a one-vertex-at-a-time pass cannot: a
//             regularised barrier untangler on every rung of the ladder, kept
//             only if no vertex of S ends a whole turn off the angle sum Q2
//             prescribes (regularisedUntangleRings).
//
//   Stage 4R  Sec. 7.2b pulls a psi_0 that lost any of the alignment back onto
//             it under the barrier before Stage 5 seeds Gamma_topo from it
//             (pullOntoAlignment): an unaligned psi_0 seeds constraints from
//             separatrices the misaligned boundary sent astray, and those
//             cannot be taken back.
//
//   Stage 1   Cancelling the same-region dipoles is retried the other way when
//             it cost the alignment (retryKeepingDipoles): the pair is how the
//             field turns through a strongly curved interface, and without it
//             no map the field integrates to follows the curve.
//
//   Stage 5   The Gamma_topo retry goes one rung tighter and then to no seeding
//             at all when the layout still leaves S partly ungridded
//             (topoNearMissLastRetry, topoRetryUnseeded).
//
// Over data/meshes/multimat at defaults: layouts reaching Definition 2.1 with
// every face four-sided went from 14 of 31 to 23, against 14 for
// Pipeline A on the same corpus; single-material, 33 of 35 to 34.
// docs/cf_flow_pipeline.md Sec. 14 has the model-by-model record.
//
// ### One material region at a time (Options::perMaterial, the default)
//
// Everything in the section above is the price of laying a multi-material model
// out as one map: the interfaces are curves inside Omega, the cut crosses them
// wherever a region is enclosed, and every stage from 3F on has to read a branch
// in two charts. The per-material mode (MaterialLayout) pays none of it. Each
// material region is laid out on its own -- this pipeline, Stages 1 to 8, on
// the region's triangles, with Stage 0's field restricted to them -- so that the
// interfaces are dS of the regions they bound, held on an axis by Q3 like any
// other dS, and no cut crosses one. The field is the same on both sides of an
// interface and tangent to it, so any two regions' layouts can be joined along
// it; what has to be agreed is only where each puts its layout vertices, and
// that is settled in rounds. The vertices one region puts on an interface are
// handed to the region across it as emitters, the neighbour's Stage 5 pairs
// each with a cone it nearly meets or with a separatrix of its own landing
// beside it -- Sec. 3.3's near-miss rule, applied across the interface -- and
// the neighbour's Stages 5 to 8 run again (relayout()) until no region is asked
// for anything new. The regions' arrangements are then glued into one
// Arrangement of S, the vertices matched across each interface made one node,
// and Stages 9 to 11 run on it as on a traced one. What the rounds leave single
// is carried on through the glued faces to dS; where even that does not close
// the layout, the whole model is laid out as well and the one that leaves less
// of S without a grid is kept (perMaterialFallback).
//
// Over data/meshes/multimat at defaults, layouts reaching Definition 2.1 with
// every face four-sided go from the whole-model route's 23 of 31 to 29, 24 of
// them glued, and those with every element's scaled Jacobian positive from 13
// to 18. docs/cf_flow_pipeline.md Sec. 15.
class TORSION {
public:
    struct Options {
        // --- Stage 0 and 0b, identical to MERIDIAN's -----------------------
        double dualMBOGamma = 10.0;
        int dualMBOMaxSteps = 500;

        // How the p=0 edge penalty combines the two elements on an edge. At
        // this order the weight is not a stabilisation parameter, it *is* the
        // discrete Laplacian, so this is the discretisation and not a knob --
        // see DualMBO::PenaltyWeight. Orthogonal is the two-point
        // (finite-volume) weight, the one that is consistent on an affine
        // field, and is what Stage 0 runs.
        DualMBO::PenaltyWeight dualMBOWeight = DualMBO::PenaltyWeight::Orthogonal;

        // Anneal tau down a geometric ladder rather than taking the single
        // tau = D^2/10 the heuristic gives.
        //
        // At the heuristic's tau one MBO step smooths over several domain
        // diameters: tau*K swamps M, the iteration reaches its fixed point in a
        // single solve, and what Stage 0 returns is a harmonic extension
        // followed by one normalisation with the threshold dynamics never
        // running. Annealing tau down, restarting each level from the previous
        // level's field with the pins re-imposed, is the MBO analogue of
        // Ginzburg-Landau's epsilon-continuation and is part of the scheme
        // rather than an optimisation of it.
        //
        // The floor is where one step stops resolving anything the mesh can
        // carry: ell = h sqrt(tau lambda) <= tauFloorEdges * h, with lambda the
        // operator's own median K_ii/M_ii (DualMBO::medianDiffusionRate), so
        // the floor means the same thing whichever weight is assembled.
        bool dualMBOTauContinuation = true;
        double dualMBOTauRatio = 0.25;
        double dualMBOTauFloorEdges = 20.0;
        int dualMBOTauLevelSteps = 2000;

        // Stage 0 override: one unit spin-4 value per triangle of the mesh the
        // layout is computed on, used instead of running the MBO solve.
        //
        // Every stage after Stage 0 reads the field through `DualMBO::u_k_prev`
        // and through nothing else -- ConeSingularities takes it in its
        // constructor and FieldFrames combs it -- so substituting the vector is
        // the whole of substituting the field. That is what makes a controlled
        // comparison possible at all: the experiment holds Stages 0b to 11
        // fixed, swaps this, and attributes the difference in the layout to the
        // field it came from rather than to a second pipeline.
        //
        // Empty runs the solver, which is the ordinary path and what the
        // pipeline itself always does. A vector whose length is not the
        // triangle count is refused and the solver runs instead, with a message
        // -- the length can legitimately change under the caller's feet when
        // Options::diskTemplates excises the inclusions before Stage 0.
        Eigen::VectorXcd externalField;

        // --- The per-material mode (MaterialLayout, Sec. 15) ----------------
        //
        // On a multi-material model, lay out each material region on its own
        // -- Stages 1 to 8 on the region's triangles alone, with this field
        // restricted to them -- match the layout vertices neighbouring regions
        // put on the interfaces they share, extend the ones that find no
        // partner through the neighbour's own map, and glue the result into
        // one arrangement of S for Stages 9 to 11. See MaterialLayout for why
        // and docs/cf_flow_pipeline.md Sec. 15 for the corpus. No effect on a
        // single-material mesh. Off runs Stages 1 to 8 on the whole model at
        // once, which is what this pipeline did until 2026-09-29.
        bool perMaterial = true;
        // How many rounds of matching at most. A round re-runs Stages 5 to 8 of
        // every region whose neighbours asked it for a new layout edge.
        int perMaterialRounds = 10;
        // Regions laid out at once; 0 is one per hardware thread. The regions
        // are independent, and the result does not depend on this.
        int perMaterialThreads = 0;
        // How far apart, in edges of the interface branch, two layout vertices
        // on opposite sides of it may be and still be matched into one vertex
        // of the glued layout (MaterialLayout::Options::matchTolerance).
        double perMaterialMatchTolerance = 3.0;
        // Re-run a region keeping its +1/-1 pairs when cancelling them left it
        // with a face that is not a simple quadrilateral. See MaterialLayout.
        bool perMaterialKeepDipoles = true;
        // When the glued layout is not valid, run the whole model as well and
        // keep whichever leaves less of S without a grid.
        bool perMaterialFallback = true;
        // Weigh where the glued nodes land when matching (no two closer than
        // a quarter of an edge unless a region's own layout already has them
        // so), and bend a matched curve into its glued node rather than move
        // only its last point. See MaterialLayout::Options::spacedMatching and
        // bendMatchedEnds; off, the matching of 2026-09-29.
        bool perMaterialSpacedMatching = true;
        bool perMaterialBendEnds = true;

        // Two hooks the per-material mode sets on each region's own run, and
        // which do nothing at their defaults. extraEmitters: vertices of dS
        // that emit a layout edge although Stage 1 put no cone there, handed
        // to Stage 5's seeding and to Stage 7 alike (see SubdomainLabels::
        // Options::extraEmitters). cancelSingleMaterialDipoles: cancel the
        // +1/-1 pairs Stage 1 finds on a single-material mesh too. Stage 1
        // cancels them only where there are interfaces, and a material region
        // laid out on its own is a single-material mesh: left alone it keeps
        // pairs the whole-model run would have cancelled inside it -- on
        // jelly_roll, 29 patches for a strip that is one.
        std::vector<int> extraEmitters;
        bool cancelSingleMaterialDipoles = false;

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
        // If cancelling cost Stage 4F the alignment, run Stages 1 to 4R again
        // on the same field keeping the pairs, and keep that run if it keeps
        // the alignment. See run().
        bool retryKeepingDipoles = true;
        bool prescribeInterfaceCones = true;
        bool interfaceCorners = true;
        bool propagateInterfaceLabels = true;
        bool seamTurnInterfaceLabels = true;
        bool lagInterfaceScales = false;

        // --- Stage 1 and 2, identical to MERIDIAN's ------------------------
        bool autoRebalance = true;
        int minBoundaryIndex = -3;
        int maxBoundaryIndex = 1;
        bool relocateFlatCones = true;
        double flatConeAngle = 5.0 * M_PI / 6.0;
        // When Stage 1 moved a +1 off a straight stretch of dS and the layout
        // that came out of it is not valid, run Stages 1 to 8 again with it
        // left where it was and keep the better of the two -- on the whole
        // model here, and per region in the per-material mode
        // (MaterialLayout::Options::keepFlatConesOnFailure).
        bool keepFlatConesOnFailure = true;
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
        // If the aligned solve inverts something Sec. 7.2a's ladder cannot clear
        // with the alignment held, try it on fewer chains -- the ones the
        // tangle touches released first, then the unstepped ones, then none --
        // and keep the first that comes out of the ladder; only if none does,
        // the one with the fewest flips. Off keeps the full alignment whatever
        // it inverts. See TORSION::buildPsi0() for why the choice is made on
        // what the ladder leaves rather than on the raw flip count, which is
        // what this did first.
        bool alignmentFallback = true;
        // Decide among those attempts by what the ladder leaves of each, the
        // first it clears winning. Off is the rule this had first, the fewest
        // raw flips, which dropped the alignment on 20 of the 31 multi-material
        // models for tangles of 2 to 30 faces.
        bool alignmentChooseByLadder = true;
        // Before falling back to the alignment on fewer chains or on none, let
        // go of only the chains the tangle touches and solve again, this many
        // times. Zero is the old three-way fallback. See TORSION::run().
        int alignmentReleaseRounds = 4;
        // Sec. 7.2b: when the psi_0 Stage 4R hands on has lost any of the
        // alignment, pull it back on under the barrier before Stage 5 seeds
        // from it. See pullOntoAlignment() in TORSION.cxx for why this is
        // Stage 5's problem and not Stage 6's.
        bool pullOntoAlignment = true;
        int pullOuterSteps = 6;
        // Read an interface branch the cutting graph crosses in the chart each
        // piece of it is combed in. See FieldFrames::Options::alignAcrossSeams:
        // off votes such a branch as one chart and holds the far side of every
        // crossing on the wrong axis.
        bool alignAcrossSeams = true;
        // Move the field's singularities onto Stage 1's cone set sector by
        // sector on dS and on the interface network, by paths that keep off
        // both, and hold the pinned faces through the re-smoothing. See
        // FieldFrames::Options::reconcileSectors.
        bool reconcileSectors = true;

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
        // Sec. 7.2c after Sec. 7.2a on every rung: the regularised untangler,
        // over this many rings round the tangle and at most this many
        // iterations of its eps schedule. Zero rings switches it off.
        int regularisedUntangleRings = 3;
        int regularisedUntangleIterations = 60;
        // The ladder is not climbed for a tangle of more than this fraction of
        // the faces, or this many, whichever is larger: such a tangle is not
        // local, and the attempt goes straight to Sec. 6.4's chain releases.
        double localUntangleFaceFraction = 0.02;
        int localUntangleFaceFloor = 200;
        // Let a seam vertex move at rungs 0 and 1 of that ladder, with its
        // partner moving under R_k so that Q4 is unchanged. See
        // TORSION::seamPairing().
        //
        // Off. It was measured once as buying nothing, and that measurement
        // was of a no-op: vertexFreedom() held every seam child at rungs 0 and
        // 1 whatever the pairing said, so relaxToKernel() never moved a pair.
        // It now frees the paired children the way Sec. 7.2c's regularised
        // pass does (which pairs the seam unconditionally, and is where the
        // freedom mattered: a fold round the tip of a slit). For the kernel pass
        // it is unmeasured since, and left off.
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
        // A second, tighter tolerance to seed from when the first one's layout
        // left a piece of S no grid covers, and zero to switch that off.
        // MERIDIAN::Options::topoNearMissRetry is the argument for it; this
        // pipeline seeds Gamma_topo with the same code and inherits the same
        // asymmetry, that a constraint can be added and never taken back.
        double topoNearMissRetry = 0.04;
        // Two more rungs of the same retry, each run only while the best
        // layout so far still leaves part of S without a grid: a tighter
        // tolerance still, and then no seeding at all, the repair loop alone.
        // See TORSION::run(). Zero, and false, stop at MERIDIAN's two rungs.
        double topoNearMissLastRetry = 0.01;
        bool topoRetryUnseeded = true;
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
        // QuadMesh::Options::collapseSpan: chords every patch of which is
        // thinner than this multiple of the target edge length are contracted
        // and the blocks either side of them merged, so that a layout finer
        // than the elements asked for does not force elements finer than that.
        // Zero keeps every chord.
        double quadCollapseSpan = 0.5;
        bool quadUseSplines = true;
        bool quadFeaturesOnTracedArcs = true;
        // QuadMesh::Options::materialsFromFaces and contractOntoFeatures: an
        // element's material from its layout face, and a contracted run of
        // nodes placed on dS or the interface. Off restores Stage 10 as it was
        // before 2026-10-01.
        bool quadMaterialsFromFaces = true;
        bool quadContractOntoFeatures = true;
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
        // Whether Options::externalField was actually used. It is refused when
        // its length does not match the triangle count, which happens whenever
        // Stage 0c excises circular inclusions before Stage 0 -- so a run with
        // `diskTemplates` on and a disk in the domain solves its own field and
        // is *not* a comparison of the caller's. An experiment that swaps
        // fields has to read this rather than assume.
        bool externalFieldUsed = false;
        bool fieldAlignedToInterfaces = false;
        int coneDipoleUnits = 0;
        int flatConesToCorners = 0;
        int flatConesInside = 0;
        // Whether the moved set's layout was not valid and the one with the +1s
        // left in place is what stands (Options::keepFlatConesOnFailure).
        bool flatConesKept = false;
        // run()'s second attempt without the cancellation: whether it ran, and
        // whether its cone set -- the pairs kept -- is the one that stands.
        bool dipoleRetryRan = false;
        bool dipolesKept = false;

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
        int alignmentSeamCrossings = 0;
        int alignmentClosedChains = 0;
        double maxAlignmentResidual = 0.0;

        // --- Stage 4F: the integration -------------------------------------
        bool integrationRan = false;
        bool integrationSolved = false;
        int integrationAlignedEdges = 0;
        double integrationAlignResidual = 0.0;
        double integrationAlignStrain = 0.0;
        bool alignmentWasDropped = false;
        int alignmentChainsReleased = 0;
        // Whether the map Stage 4R handed on held the whole of Sec. 6.4's
        // alignment exactly, before Sec. 7.2b had anything to pull back.
        bool alignmentHeld = true;
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
        // Sec. 7.2b: whether the pull ran, came back injective, and ended on
        // the exact alignment (projected) or only near it (penalised).
        bool pullRan = false;
        bool pullKept = false;
        bool pullProjected = false;
        double pullBoundaryResidual = 0.0;
        double pullFeatureResidual = 0.0;

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
        // Which seeding tolerance the layout below was built at, and whether
        // Options::topoNearMissRetry had to run a second attempt.
        double topoNearMissUsed = 0.0;
        bool topoNearMissRetried = false;

        int interfaceCorners = 0;
        int interfaceLabelsCorrected = 0;
        int interfaceChainsSeamFlipped = 0;
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

        // --- The per-material mode (MaterialLayout) -------------------------
        // Whether it ran, and whether its glued layout is the one Stages 9 to
        // 11 were run on -- false after a fallback that the whole model won.
        bool perMaterialRan = false;
        bool perMaterialKept = false;
        bool perMaterialFallbackRan = false;
        int perMaterialRegions = 0;
        int perMaterialRegionsValid = 0;
        int perMaterialRounds = 0;
        bool perMaterialConverged = false;
        int perMaterialEmitters = 0;
        int perMaterialMatched = 0;
        int perMaterialUnmatched = 0;
        int perMaterialDipolesKept = 0;
        double perMaterialSeconds = 0.0;

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

    // Stages 5 to 8 again, on the psi_0 run() already built, with `emitters`
    // as Options::extraEmitters. What the per-material mode's matching rounds
    // call on each region: a vertex where a neighbour's layout meets the
    // shared interface is a statement about Psi and not about psi_0 -- it
    // changes which separatrices Stage 5 pairs and which Stage 7 traces, and
    // nothing before -- so a round pays for Stages 5 to 8 only. False if
    // run() never reached Stage 5.
    bool relayout(const std::vector<int> &emitters);

    const DualMBO& getField() const { return *field; }
    // The per-material mode's regions, their layouts and the matching; only
    // after a run() on a multi-material mesh with Options::perMaterial.
    bool hasMaterialLayout() const { return materials != nullptr; }
    const MaterialLayout& getMaterialLayout() const;
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

    // --- Stage 4R's three pieces, public for the same reason -------------
    //
    // Sec. 7.2a is a ladder rather than a single pass, and a caller that drives
    // the stages itself has to climb the same one or it is not running Stage 4R
    // at all -- it is running Sec. 7.2's Tutte pass on every model, which
    // throws away the cone angles, the seam and Sec. 6.4's alignment on the
    // models the local repair would have kept them on. That is a different
    // psi_0, and every stage after it sees the difference.

    // The Jacobian of a piecewise-linear map on Omega, one 2x2 row-major per
    // face: row 0 is grad u, row 1 is grad v, which is the layout
    // FieldIntegration::Options::targetJacobian is read in. A map's own
    // Jacobian is a legal target for the fit that produced it, which is what
    // makes Sec. 6.5's re-projection a swap of one argument.
    static std::vector<std::array<double, 4>> jacobianOf(const Mesh &om,
                                                         const std::vector<Point> &uv);

    // Which coordinates each vertex of Omega may move in without undoing an
    // equality Stage 4F imposed, at rung `level` of Sec. 7.2a's ladder: 0 holds
    // the seam and the alignment, 1 frees the alignment, 2 frees both.
    // `usedAxis` is the axis the integration was actually solved with -- empty
    // when it was solved free. See TORSION.cxx.
    static std::vector<unsigned char> vertexFreedom(const ConeCut &cut,
                                                    const std::vector<int> &usedAxis,
                                                    int level);

    // Sec. 7.2a itself: move the interior vertices around each inverted face to
    // the Chebyshev centre of their one-ring kernels, within `freedom`, for at
    // most `sweeps` sweeps. `partner`/`partnerK` are Options::pairSeamInUntangle's
    // seam pairing and may be empty. Returns the number of faces still
    // inverted, or -1 if the arguments do not match the mesh.
    static int relaxToKernel(const Mesh &om, std::vector<Point> &uv,
                             const std::vector<unsigned char> &freedom,
                             const std::vector<int> &partner,
                             const std::vector<int> &partnerK, int sweeps);

    // Sec. 7.2c: the fold relaxToKernel() cannot open -- a run of triangles
    // turned over together, typically round a cone whose fan the frame asked
    // for the wrong total angle -- opened by minimising Garanzha et al.'s
    // regularised barrier energy over the vertices within `rings` rings of it,
    // within `freedom`. `target` is the frame, read for its scale only. Returns
    // the number of faces still inverted; the map is changed only if that
    // number went down. See TORSION.cxx.
    // `cutToOriginal` is ConeCut::getCutVertexToOriginal() and
    // `prescribedAngle` Q2's angle sum per vertex of S (2 pi - (pi/2) I inside,
    // pi - (pi/2) I on dS), for the one check the result has to pass: that no
    // vertex of S -- a seam vertex included, whose star is split between two
    // children -- is a whole turn off it. Either may be empty, and then only
    // the interior vertices of Omega are held to one turn.
    //
    // `partner`/`partnerK` pair the two children of each seam vertex the way
    // relaxToKernel() reads them: the pair moves together, one child by d and
    // the other by R_k d, which keeps Q4 exact while the fold round a slit's
    // tip opens. Either may be empty.
    static int untangleRegularised(const Mesh &om, std::vector<Point> &uv,
                                   const std::vector<unsigned char> &freedom,
                                   const std::vector<int> &partner,
                                   const std::vector<int> &partnerK,
                                   const std::vector<std::array<double, 4>> &target,
                                   const std::vector<int> &cutToOriginal,
                                   const std::vector<double> &prescribedAngle,
                                   int rings, int maxIterations);

    // --- Stages 4F and 4R as one step, public for the same reason ---------
    //
    // What buildPsi0() reads, and what it hands back. The references are the
    // caller's and must outlive the call; interfaces and coneMetric may be null
    // (a single material; no cone metric was built).
    struct MapStage {
        const Mesh &mesh;                     // S
        const ConeCut &cut;
        const ConeSingularities &cones;
        const FieldFrames &frames;
        const Immersion &scaffold;            // the seam pairing, as run() builds it
        const Interfaces *interfaces;
        const ConeMetric *coneMetric;
    };
    struct Psi0 {
        std::unique_ptr<FieldIntegration> integration;   // the solve that was kept
        std::unique_ptr<TutteEmbedding> tutte;           // only if Sec. 7.2 ran
        std::vector<int> usedAxis;                       // the alignment it holds
        std::vector<Point> integratedMap;                // as the solve left it
        std::vector<Point> psi0;                         // what Stage 4 is handed
    };
    // Sec. 6.4's attempts, Sec. 7.2a-c's ladder, Sec. 7.2's Tutte pass and Sec.
    // 7.2b's pull, in run()'s order and under `options`; the Stage 4F and 4R
    // fields of `status` and its messages are written as run() writes them.
    // False when there is no psi_0 to hand on.
    static bool buildPsi0(const MapStage &in, const Options &options, Status &status,
                          Psi0 &out);

    // Stage 0's field, one level of its tau-continuation at a time. runField()
    // builds and steps every level through these, and so does the viewer's
    // TORSION mode, so that the field it shows -- and lays out, a region at a
    // time or whole -- is the field this class computes. It was not, for a
    // while: the viewer stepped MERIDIAN's field, a single tau at DualMBO's
    // default penalty weight, which on rt_mushroom is two MBO steps from the
    // linear solve and lays out with folds along the interfaces that the
    // pipeline's field does not have.
    //
    // prepareFieldLevel() takes a DualMBO freshly constructed with
    // Options::dualMBOMaxSteps and dualMBOGamma and makes it one level:
    // Options::dualMBOWeight, the interface network as Dirichlet data when
    // Options::alignFieldToInterfaces and the model has more than one material,
    // the level's tau scale, the assembly, and -- when `carried` has one value
    // per triangle -- the previous level's field as the start, with this
    // level's own Dirichlet data re-imposed on it. True when the interfaces
    // were imposed.
    static bool prepareFieldLevel(DualMBO &field, const Interfaces *interfaces,
                                  const Options &options, double tauScale,
                                  const Eigen::VectorXcd &carried);
    // The tau scales the continuation runs through, 1 (the heuristic) first,
    // read off level 0 once it is assembled; {1} when it is off.
    static std::vector<double> fieldTauLadder(const Mesh &mesh, const DualMBO &level0,
                                              const Options &options);
    // The steps any one level of a ladder `levels` long may take, and the
    // test each level stops at before that.
    static int fieldLevelCap(const Options &options, std::size_t levels);
    static bool fieldLevelConverged(const DualMBO &field);

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
    // Stages 5 to 8 at one Gamma_topo seeding tolerance. False when a stage
    // stopped the pipeline, in which case there is nothing to retry.
    // `seed` is whether Gamma_topo is seeded from psi_0 at all; the repair
    // loop runs whenever Options::seedTopoConstraints is on, seeded or not.
    bool runLayoutStages(double nearMiss, bool seed);
    // Stages 5 to 8 with the Gamma_topo retry: runLayoutStages() at
    // Options::topoNearMiss and then down the ladder while the best layout so
    // far leaves part of S without a grid. Starts from statusBeforeLayout.
    bool runLayoutLadder();
    // Stages 1 to 8 on the whole model at once -- this pipeline as it was
    // before the per-material mode. `proceed` is set when there is an
    // arrangement for Stages 9 to 11; otherwise the return value is run()'s.
    // runWholeModelOnce() is one attempt; runWholeModel() makes a second with
    // Stage 1's flat +1s left in place when the first moved some and did not
    // come out valid (Options::keepFlatConesOnFailure).
    bool runWholeModel(bool &proceed);
    bool runWholeModelOnce(bool &proceed);
    // The per-material mode: MaterialLayout on the field Stage 0 solved, its
    // glued arrangement moved into `arrangement`. True when that layout is
    // valid.
    bool runPerMaterial();
    // Stages 9 to 11 on `arrangement`, whichever route built it.
    bool runMeshStages();
    // Stages 0c, 0b and 0, which are MERIDIAN's unchanged and run once.
    void runFieldFront();
    // Stages 1 and 2 on that field, MERIDIAN's unchanged but for whether the
    // same-region dipoles are cancelled. Returns false when the cone set is
    // inadmissible or the cut is not a disk, which are the two things every
    // route downstream depends on.
    bool runConeFront(bool cancelDipoles);
    // Stages 4C, 3F, 4F and 4R, from the cone set to psi0Map. False at the
    // stops no later stage survives.
    bool runMapStages();
    // The two children of each seam vertex of Omega and the quarter turn
    // between their displacements, so that Sec. 7.2a can move a tangle sitting
    // on the seam without letting go of Q4. See TORSION.cxx.
    static void seamPairing(const ConeCut &cut, const Immersion &scaffold,
                            std::vector<int> &mate, std::vector<int> &turn);
    // Sec. 7.2, the Tutte pass. Takes out.integratedMap for its reference and
    // returns an untangled map, or an empty vector when it could not; leaves the
    // embedding in out.tutte.
    static std::vector<Point> untangle(const MapStage &in, const Options &options,
                                       Status &status, Psi0 &out);
    // Sec. 7.2b. Pull an injective map that lost some of Sec. 6.4's alignment
    // back onto it -- labels read off `labelMap`, the solve that held all of
    // it -- under the barrier, and project onto the exact alignment `axis`
    // when that inverts nothing. Empty when there was nothing to pull.
    static std::vector<Point> pullOntoAlignment(const MapStage &in, const Options &options,
                                                Status &status,
                                                const std::vector<Point> &start,
                                                const std::vector<Point> &labelMap,
                                                const std::vector<int> &axis);

    std::shared_ptr<Mesh> mesh;
    // What run() was handed, kept only when Stage 0c replaced it.
    std::shared_ptr<Mesh> inputMesh;
    Options options;
    Status status;

    std::unique_ptr<DualMBO> field;
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
    std::unique_ptr<MaterialLayout> materials;
    // The status Stage 5 starts from, so that relayout() and the retry ladder
    // can run Stages 5 to 8 again without Stage 5's messages piling up.
    Status statusBeforeLayout;

    // The circular inclusions Stage 0c took out. Empty unless
    // Options::diskTemplates.
    std::vector<DiskTemplate::Inclusion> inclusions;

    // The indices the cross field read, before Stage 1 moved any of them; the
    // set Sec. 5.1's audit is run against.
    std::vector<int> fieldIndex;
    std::vector<Point> integratedMap;
    // psi_0 as Stage 4R hands it to the Immersion.
    std::vector<Point> psi0Map;
    // The alignment Stage 4F actually used, after its fallback chose among the
    // three. Sec. 6.5's re-projection has to hold the same one.
    std::vector<int> usedAxis;
};

#endif // __TORSION_HXX__
