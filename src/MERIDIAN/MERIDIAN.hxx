#ifndef __MERIDIAN_HXX__
#define __MERIDIAN_HXX__

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "MERIDIAN/Arrangement.hxx"
#include "MERIDIAN/ConeCut.hxx"
#include "MERIDIAN/ConeSingularities.hxx"
#include "MERIDIAN/DiskTemplate.hxx"
#include "MERIDIAN/Immersion.hxx"
#include "MERIDIAN/Interfaces.hxx"
#include "MERIDIAN/LayoutEnergy.hxx"
#include "MERIDIAN/QuadMesh.hxx"
#include "MERIDIAN/RicciFlow.hxx"
#include "MERIDIAN/Separatrices.hxx"
#include "MERIDIAN/SplineFit.hxx"
#include "MERIDIAN/SubdomainLabels.hxx"
#include "mesh/Mesh.hxx"

class DualMBO;

// Stages 1 to 9 of
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
//   Stage 0b             The material interface network, on a mesh whose
//                        triangles carry more than one material tag. The
//                        paper's features are trim curves and creases and it
//                        takes them as an input to Stage 0; here they are the
//                        boundaries between the physical surfaces of the .geo,
//                        and the graph they form -- its junctions, its landings
//                        on dS, its corners, and the quarter turns the layout
//                        must make in every sector at each of them -- is what
//                        Stages 1, 5, 6, 7 and 8 are then written against.
//                        Nothing on a single-material mesh.
//                        -> Interfaces
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
//                        it is what Stages 8 and 9 partition and fit.
//                        -> Separatrices
//
//   Stage 8  Sec. 4      The arrangement: the curves and dS cut into arcs at
//                        every node they meet, assembled into a half-edge
//                        structure, and its faces enumerated. This is where a
//                        bundle of curves becomes a partition of S, and where
//                        the paper's validation -- four corners to a face, an
//                        arc between two faces, a cone with the arcs its index
//                        prescribes -- is actually run.
//                        -> Arrangement
//
//   Stage 9  Sec. 5      The spline reconstruction: one cubic B-spline per arc,
//                        fitted once and shared by both patches that meet on
//                        it, and a bicubic Coons patch per face built from the
//                        four of them at the control-point level. Fitting each
//                        arc exactly once is what makes the output watertight.
//                        -> SplineFit
//
//   Stage 10             The quadrilateral mesh on those patches. One integer
//                        per *chord* of the layout -- the class of arcs that
//                        must agree because they face each other across a
//                        patch -- chosen against a global target edge length,
//                        then each arc cut once at equal arc length and both
//                        its patches given the same nodes.
//                        -> QuadMesh
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
// ### Multi-material domains
//
// A material interface is a feature in exactly the paper's sense -- a curve of
// the input the output has to keep -- and E3 of Sec. 3.3 is already written for
// features. That is not enough, and the reason is that E3 is a statement about
// one curve at a time while everything that makes a *network* of curves part of
// a layout lives between them. Four things had to be added, and each of them
// fails visibly without the others:
//
//   Stage 0b   the network itself, and the quarter turns the layout must make
//              in each sector at each of its nodes. Two conditions: the local
//              one at a node, which is arithmetic on the incident directions,
//              and the global one on each material region, which is the
//              discrete Gauss-Bonnet count of Eq. (4) restricted to it. The
//              second is not optional and is not a tolerance: a quarter disk
//              inside a square whose interface meets the boundary at right
//              angles at both ends has three quarter turns on its boundary and
//              needs four, and until the fourth is put somewhere E3 and E6 are
//              asking for a map that does not exist.
//
//   Stage 1    the cone index at every node of that network, prescribed from
//              the geometry rather than read off the field. The field is
//              boundary-aligned and has no boundary condition on an interface,
//              so at a junction it reports whatever propagated in from
//              elsewhere.
//
//   Stage 6    E6, the sector condition imposed on Psi. See LayoutEnergy.
//
//   Stages 7   the interfaces as arcs of the layout, and the layout edges that
//   and 8      leave a node whose sector is more than one quadrilateral wide.
//              An interface is a layout edge whether or not a separatrix
//              happens to lie on it, so Stage 8 takes the branches directly;
//              and a node with a three-quarter sector emits two curves into it
//              that no cone emits, so Stage 7 emits from the network's nodes as
//              well as from P, suppressing the rays that would merely lie along
//              the interface a second time.
//
// What comes out is a mesh whose every element lies inside one material and
// whose two sides of every interface share their nodes, which is what an
// analysis code on a multi-material domain needs and what a layout ignoring the
// tags cannot give it.
//
// What comes out at the end is a set of bicubic patches meeting C0 across their
// shared boundary curves and C2 inside, together with everything they were
// derived from: Psi, the separatrices, the arrangement, and a conforming
// all-quadrilateral mesh of the patches at a prescribed edge length. The
// uniform knot-insertion refinement and the export to an analysis solver read
// them and are not implemented here.
//
// run() returns Definition 2.1's verdict and nothing else: whether Psi is a
// quadrilateral layout. Stages 8 and 9 can fail on a map that is a valid layout
// in that sense -- a T-junction is quantisation that has not converged, not a
// broken map -- and they report separately, because the remedy for both is the
// one Sec. 4 gives and Sec. 3.3's repair loop has already applied: more
// lambda_5 and another Gamma_topo constraint.
class MERIDIAN {
public:
    struct Options {
        // Stage 0: the cross field the cones are read off. The defaults match
        // TestDualMBO.
        double dualMBOGamma = 10.0;
        int dualMBOMaxSteps = 500;

        // Stage 0b: the material interface network. No effect on a
        // single-material mesh, where there are no interfaces to find, so this
        // is on by default and the single-material pipeline is unchanged.
        bool materialInterfaces = true;
        // Where an interface turns sharply enough that the layout has to turn
        // with it, and how many arcs a closed interface loop is cut into.
        double interfaceKinkAngle = M_PI / 4.0;
        int interfaceLoopSplits = 4;
        // Leave a disk's rim unsplit: its own cones already cross it four
        // times, and a node anywhere else on it emits a layout edge nothing
        // receives. See Interfaces::Options::splitCircleLoops.
        bool splitCircleLoops = false;
        // Stages 0c and 11: take every circular inclusion out of the layout
        // problem and put its O-grid back from a template. See DiskTemplate,
        // where the reasoning and the parity argument live.
        //
        // The short version. The O-grid is the only block structure a disk has
        // -- a core block and a ring of blocks around it -- and the cross field
        // finds it unaided: an inclusion comes out of Stage 1 with eight cones,
        // four +1 inside at about 0.6 of the radius and four -1 out in the
        // matrix at about 1.4. Asking Stage 6 to place them is what does not
        // work. Eight cones per disk is eighty on
        // data/meshes/multimat/bubbles, the continuation of Sec. 3.3 walls out
        // somewhere near forty interior cones whatever the geometry, and what
        // comes out has 57% of the model unmeshed.
        //
        // So the disks are excised before Stage 0b and the layout is asked for
        // the matrix with holes, which is a different problem and a far easier
        // one: on bubbles it reaches Definition 2.1 with every face a
        // quadrilateral, every element positive and nothing left unmeshed.
        // Stage 11 then fills each hole from the template, which costs no solve.
        //
        // The one thing the excised problem has to be told is that each rim
        // must carry an even number of edges, because a quadrangulation of a
        // disk cannot have an odd boundary. That goes through
        // QuadMesh::Options::evenLoops and is arranged by moving a chord or
        // two by one edge.
        bool diskTemplates = false;
        // How square the core block's boundary is, how many rows of elements
        // the ring has (0 chooses it from the rim spacing), and how hard the
        // template is smoothed. See DiskTemplate::Options.
        double diskCoreSquareness = 0.55;
        int diskRingDepth = 0;
        int diskSmoothingPasses = 300;
        // Align the Stage 0 cross field to the interfaces as well as to dS.
        //
        // Without this the field has no boundary condition on an interface and
        // runs straight through it, so the cones each material region needs in
        // its *interior* are never read off at all -- and Stage 0b then has to
        // put the missing turn somewhere, which means a corner in the middle of
        // an interface that does not turn. See Interfaces::balance and
        // DualMBO::setAlignedInteriorEdges. Off is the behaviour before this,
        // which is what the --no-field-interfaces flag measures against.
        bool alignFieldToInterfaces = true;
        // Annihilate the +1/-1 cone pairs an interface-aligned field puts on the
        // concave side of a strongly curved interface. Both members are in one
        // material region, so the pair is in no region's Gauss-Bonnet count and
        // is a property of the smoothest field rather than of the layout. See
        // ConeSingularities::cancelDipoles; the --no-cancel-dipoles flag is
        // what measures it.
        bool cancelInterfaceDipoles = true;
        // Hand Stage 1 the cone index the interface geometry fixes at every
        // node of that network, overriding what the cross field read there.
        // See ConeSingularities::prescribe for why the field cannot supply it.
        bool prescribeInterfaceCones = true;
        // Stage 6's E6, the sector condition at those nodes. Off leaves the
        // interfaces to E3 alone, which is what the pipeline did before this
        // stage existed and is what the --no-e6 flag is for measuring against.
        bool interfaceCorners = true;
        bool propagateInterfaceLabels = true;
        // See SubdomainLabels::Options::seamTurnInterfaceLabels.
        bool seamTurnInterfaceLabels = true;
        // See LayoutEnergy::Options::lagInterfaceScales; off, and the numbers
        // that say why are there.
        bool lagInterfaceScales = false;

        // Stage 1
        bool autoRebalance = true;   // restore Eq. (4) by moving boundary cones
        int minBoundaryIndex = -3;
        int maxBoundaryIndex = 1;
        // Move every +1 cone off dS where dS runs straight (more than
        // flatConeAngle across the vertex, in radians): to the convex corner at
        // the end of its side when that one has room, otherwise inside. Off,
        // such a cone is a patch corner of pi on the model, which Stage 10
        // meshes as one element with a node in the middle of a side. See
        // ConeSingularities::relocateFlatCones.
        bool relocateFlatCones = true;
        double flatConeAngle = 5.0 * M_PI / 6.0;
        // When the layout that came of that move is not valid, run Stages 1 to
        // 11 again with the cones left where they were and keep the better.
        bool keepFlatConesOnFailure = true;

        // Stage 2. Where a cone's arc of the cutting graph is allowed to
        // stop: dS itself, or (false) the nearest thing already cut, which is
        // the letter of Sec. 3.2.2 and leaves G a tree with junctions in it.
        // See ConeCut::Options::conesToBoundary.
        bool coneCutsToBoundary = true;
        // How dearly the cutting graph avoids the material interface network.
        // See ConeCut::Options::interfaceAvoidance; 0 routes as if the network
        // were not there, which is what Stage 2 did before it knew.
        double coneCutInterfaceAvoidance = 1.0;

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
        // Stage 6's stall escape that swaps E1's reference metric between the
        // Ricci one and the model's own Euclidean one. Off, because a planar
        // model has no cones in its Euclidean metric and the switch then
        // shreds the cone one-rings to satisfy Q2; LayoutEnergy.hxx, "Why the
        // reference switch is off here", has the measurement.
        bool alternateReference = false;
        bool relabelBetweenSteps = true;

        // Gamma_topo. Remark 3.1: without the connectivity constraints a
        // topologically valid layout still exists but is made of slivers --
        // Fig. 13 is 1898 patches against 73 -- so they are seeded
        // automatically from near-misses of the separatrices on psi_R.
        bool seedTopoConstraints = true;
        double topoNearMiss = 0.15;  // as a fraction of the mean cone spacing

        // A second, tighter tolerance to seed from when the first one's layout
        // came out with a piece of S no grid covers. Zero switches the retry
        // off and leaves `topoNearMiss` as the only attempt.
        //
        // The two ends of the tolerance fail differently and only one of them
        // is recoverable. Seed too tightly and connections that exist are
        // missed, which Sec. 3.3's repair then puts back one at a time, nearest
        // miss first, rolling a round back when it lands worse -- so an
        // under-seeded start is a start the loop can walk away from. Seed too
        // widely and E5 is handed a set of integral curves no map has before
        // the continuation begins; nothing downstream can drop a constraint
        // again, and what comes back is a compromise with det J pushed toward
        // zero. The retry therefore runs downward and never upward.
        //
        // It is *conditional* because the wider tolerance is worth having when
        // it works: it joins more cones, so the layout has fewer and larger
        // patches and the elements come out nearer the target edge length. On
        // the nine box-with-disks domains in data/meshes/disks the default
        // tolerance already leaves nothing unmeshed, and re-seeding them at
        // 0.04 would roughly triple the rms log edge ratio for no coverage at
        // all. multimat/bubbles is the case that needs it: ten inclusions leave
        // 12% of S in faces without four corners at 0.15, and none at 0.04.
        //
        // The measure both attempts are scored on is that same fraction -- the
        // area of S in faces that are not quadrilaterals with one arc a side --
        // with the number of unresolved separatrices breaking a tie and the
        // first attempt winning an exact one.
        double topoNearMissRetry = 0.04;
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

        // Stage 8, Sec. 4. mergeTolerance is how near two nodes must be, as a
        // fraction of the diagonal of S, to be one node; cornerTolerance is how
        // far a sector may be from a whole number of quarter turns and still be
        // read as that number.
        bool runArrangement = true;
        double arrangementMerge = 1e-7;
        double arrangementCorner = 0.35;
        double arrangementCollapse = 1e-4;
        // End a separatrix Q5 did not close at the first layout edge it meets,
        // as a T-junction, rather than leaving it with a loose end in the
        // middle of a face. See Arrangement::Options::trimUnresolvedAtCrossings.
        bool arrangementTrim = true;

        // Stage 9, Sec. 5. Three cubic segments per arc is 6x6 control points
        // per patch, which is what the paper uses for the firewall and the
        // speaker frame; four is the 7x7 of the shock house and the chassis.
        bool runSplines = true;
        int splineSegments = 3;
        int splineSamples = 8;
        // Fit the arcs the input gave -- dS and the material interface network
        // -- instead of carrying them as the polylines the mesh has them as.
        // See SplineFit::Options::fitBoundaryArcs.
        bool fitBoundaryArcs = false;
        bool fitInterfaceArcs = false;

        // Stage 10: the quadrilateral mesh on those patches. The target edge
        // length is absolute, in the units of the model; the corpus is
        // normalised into [0,1]^2, so 0.05 is one twentieth of the model and
        // is the default for that reason. See QuadMesh for how one integer per
        // chord is chosen from it.
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
        // Place the nodes along the Stage 9 fits. Off, they go on the traced
        // polylines instead, which tells a meshing artefact apart from a
        // fitting one. See QuadMesh::Options::useSplines.
        bool quadUseSplines = true;
        // See QuadMesh::Options::featuresOnTracedArcs.
        bool quadFeaturesOnTracedArcs = true;
        // QuadMesh::Options::materialsFromFaces and contractOntoFeatures: an
        // element's material from its layout face, and a contracted run of
        // nodes placed on dS or the interface. Off restores Stage 10 as it was
        // before 2026-10-01.
        bool quadMaterialsFromFaces = true;
        bool quadContractOntoFeatures = true;
        // Winslow sweeps over the interior of each block, the boundary held.
        // Zero leaves the transfinite grid alone. See QuadMesh::smooth().
        int quadSmoothingPasses = 500;
        double quadSmoothingThreshold = 0.0;

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
        // Whether a repair round is judged on the arrangement its curves cut S
        // into as well as on the curves themselves. See
        // RepairOptions::scoreArrangement.
        bool repairScoresArrangement = true;
        // Rounds in a row that may fail to improve before the loop gives up.
        int repairPatience = 2;
    };

    struct Status {
        int mboSteps = 0;
        bool fieldConverged = false;
        // Whether the interfaces were Dirichlet data for the field, Stage 0.
        bool fieldAlignedToInterfaces = false;
        // Index units cancelDipoles() annihilated, Stage 1.
        int coneDipoleUnits = 0;
        // Stage 1's +1 cones taken off straight dS (Options::relocateFlatCones):
        // to a corner further along the side, and inside.
        int flatConesToCorners = 0;
        int flatConesInside = 0;
        // Whether that move's layout was not valid and the one with the +1s left
        // in place is what stands (Options::keepFlatConesOnFailure).
        bool flatConesKept = false;
        // Stage 0c: circular inclusions found, and the triangles they took out
        // of the layout problem with them.
        int diskInclusions = 0;
        int diskTrianglesExcised = 0;

        // Stage 0b
        int materials = 1;
        int interfaceEdges = 0;
        int interfaceBranches = 0;
        int interfaceNodes = 0;
        int interfaceIllPosedNodes = 0;
        double interfaceWorstSector = 0.0;   // radians
        int interfaceConesPrescribed = 0;
        int interfacePrescriptionShift = 0;  // index units the field was out by
        int regions = 0;                     // connected pieces of one material
        int regionsBalanced = 0;             // of those, satisfying their own Eq. (4)
        int regionQuartersMoved = 0;
        int regionCornersInserted = 0;
        int worstRegionDeficit = 0;

        bool conesAdmissible = false;
        int interiorCones = 0;
        int boundaryCones = 0;

        bool cutIsDisk = false;
        bool allConesOnBoundary = false;
        // Stage 2 against the interface network: arcs of G lying along an
        // interface, vertices of G on the network at all, and how many of those
        // are nodes. Zero on a single-material mesh, and zero on a
        // multi-material one whose cones all have an interface-free route to dS.
        int cutInterfaceEdges = 0;
        int cutInterfaceVertices = 0;
        int cutInterfaceNodes = 0;

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
        // Which seeding tolerance the layout below was actually built at, and
        // whether Options::topoNearMissRetry had to run a second attempt.
        double topoNearMissUsed = 0.0;
        bool topoNearMissRetried = false;

        // Stage 5/6, the interface terms
        int interfaceCorners = 0;
        int interfaceLabelsCorrected = 0;
        int interfaceChainsSeamFlipped = 0;  // labels a crossing of G turned
        double interfaceResidual = 0.0;      // radians
        int interfaceCornerChanges = 0;
        bool interfacesAligned = false;

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

        // Stage 8
        bool arrangementRan = false;
        int layoutNodes = 0;
        int layoutArcs = 0;
        int layoutPatches = 0;
        int layoutQuads = 0;       // patches with four corners and one arc a side
        double layoutCoverage = 0.0;
        bool arrangementValid = false;

        // Stage 10
        bool quadMeshRan = false;
        int meshVertices = 0;
        int meshQuads = 0;
        int meshChords = 0;
        int meshUnmeshedPatches = 0;
        double meshMinScaledJacobian = 0.0;
        bool meshConforming = false;
        bool meshValid = false;
        // The rim loops Stage 11 needs even, and what the interval assignment
        // had to do about them. See QuadMesh::fixLoopParity.
        int meshOddLoops = 0;
        int meshParityChordsMoved = 0;
        int meshOddLoopsLeft = 0;

        // Stage 11: the O-grid templates.
        bool diskTemplatesRan = false;
        int diskTemplatesFilled = 0;
        int diskTemplateBlocks = 0;
        int diskTemplateQuads = 0;
        int diskTemplatesRefused = 0;
        double diskTemplateMinScaledJacobian = 0.0;
        // Over the matrix mesh and the templates together, which is the mesh a
        // caller actually gets.
        int mergedVertices = 0;
        int mergedQuads = 0;
        double mergedMinScaledJacobian = 0.0;
        bool diskTemplatesValid = false;

        // Stage 9
        bool splinesRan = false;
        int splinePatches = 0;
        int splineControlPoints = 0;   // per arc
        double splineMaxDeviation = 0.0;
        bool splinesWatertight = false;
        bool splinesValid = false;

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

        // Score each round on the arrangement the curves cut S into, not only
        // on the curves.
        //
        // Q5 is a statement about the curves and the loop has to have it, but a
        // caller does not want curves, it wants patches: the failure that costs
        // it a patch is a face with three corners or a cone one arc short, and
        // those are counted in Stage 8 and nowhere earlier. Two rounds can be
        // indistinguishable in unterminated curves and grazes and differ by
        // three broken patches, and without this the loop cannot tell them
        // apart and keeps whichever came second. Stage 8 costs well under a
        // percent of a round, so there is no reason not to look.
        bool scoreArrangement = true;
        Arrangement::Options arrangement;

        // How many rounds in a row may fail to improve before the loop gives
        // up. One is the old behaviour -- stop at the first round that does not
        // pay -- and it is a round too few: a constraint that only helps once
        // its partner is also in costs on the round it is added and pays on the
        // next, which is exactly what a pair of clustered cones does. Whatever
        // the loop ends on, the map and the paths of the best round it saw are
        // what it returns, so patience can only cost time.
        int patience = 2;
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

    // The fraction of S that Stage 10 has no grid for: the area of the faces
    // inside S that are not quadrilaterals with one arc on each side, over the
    // area of S. Zero is Sec. 4's validation passing on the only count that
    // decides whether the model comes out meshed. Static, and taking the
    // arrangement rather than reading a member, so that TORSION scores its own
    // attempts with this function and not a second copy of it.
    static double unmeshableFraction(const Arrangement &arr);

    explicit MERIDIAN(std::shared_ptr<Mesh> mesh);
    MERIDIAN(std::shared_ptr<Mesh> mesh, const Options &opts);
    ~MERIDIAN();

    // Runs the three stages in order and stops at the first one that fails in a
    // way the next stage cannot survive: an inadmissible cone set stops the
    // pipeline, a cut that came out as something other than a disk does not
    // (Ricci flow does not use it, so the metric is still worth computing and
    // the failure is still worth reporting).
    bool run();

    const DualMBO& getField() const { return *field; }
    const Interfaces& getInterfaces() const { return *interfaces; }
    const ConeSingularities& getCones() const { return *cones; }
    const ConeCut& getCut() const { return *cutter; }
    const RicciFlow& getRicci() const { return *ricci; }
    const Immersion& getImmersion() const { return *immersion; }
    const SubdomainLabels& getLabels() const { return *labels; }
    const LayoutEnergy& getLayout() const { return *layout; }
    const Separatrices& getSeparatrices() const { return *separatrices; }
    const Arrangement& getArrangement() const { return *arrangement; }
    const SplineFit& getSplines() const { return *splines; }
    const QuadMesh& getQuadMesh() const { return *quads; }
    const DiskTemplate& getDiskTemplate() const { return *diskFill; }
    const std::vector<DiskTemplate::Inclusion>& getInclusions() const { return inclusions; }

    // Null until the stage that builds them has run.
    bool hasInterfaces() const { return interfaces != nullptr; }
    bool hasCones() const { return cones != nullptr; }
    bool hasCut() const { return cutter != nullptr; }
    bool hasRicci() const { return ricci != nullptr; }
    bool hasImmersion() const { return immersion != nullptr; }
    bool hasLabels() const { return labels != nullptr; }
    bool hasLayout() const { return layout != nullptr; }
    bool hasSeparatrices() const { return separatrices != nullptr; }
    bool hasArrangement() const { return arrangement != nullptr; }
    bool hasSplines() const { return splines != nullptr; }
    bool hasQuadMesh() const { return quads != nullptr; }
    bool hasDiskTemplate() const { return diskFill != nullptr; }

    const Status& getStatus() const { return status; }

    // The mesh the layout was computed on. With Options::diskTemplates this is
    // the *excised* mesh -- the matrix with a hole where each inclusion was --
    // and every index in every stage downstream of Stage 0c refers to it.
    // getInputMesh() is what came in.
    const Mesh& getMesh() const { return *mesh; }
    std::shared_ptr<Mesh> getMeshPtr() const { return mesh; }
    const Mesh& getInputMesh() const { return inputMesh ? *inputMesh : *mesh; }

private:
    void runField();
    // Stage 0c: find the circular inclusions and take them out of `mesh`.
    void exciseDisks();
    // Stages 5 to 8 at one Gamma_topo seeding tolerance. False when a stage
    // stopped the pipeline, in which case there is nothing to retry.
    bool runLayoutStages(double nearMiss);
    // Stages 1 to 11, on the field and the interface network run() built;
    // what run() returns. run() calls it a second time with Stage 1's flat +1s
    // left in place when the first moved some and did not come out valid.
    bool runFromCones(size_t balanceMessagesSeen);

    std::shared_ptr<Mesh> mesh;
    // What run() was handed, kept only when Stage 0c replaced it.
    std::shared_ptr<Mesh> inputMesh;
    Options options;
    Status status;

    std::unique_ptr<DualMBO> field;
    std::unique_ptr<Interfaces> interfaces;
    std::unique_ptr<ConeSingularities> cones;
    std::unique_ptr<ConeCut> cutter;
    std::unique_ptr<RicciFlow> ricci;
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
};

#endif // __MERIDIAN_HXX__
