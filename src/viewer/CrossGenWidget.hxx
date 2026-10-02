#pragma once

#include <QOpenGLWidget>
#include <QTimer>
#include <QPoint>
#include <QString>

#include <memory>
#include <optional>
#include <string>

#include "viewer/ViewerTypes.hxx"
#include "viewer/Render.hxx"
#include "viewer/Geometry.hxx"

#include "mesh/BlockQuadMesh.hxx"
#include "mesh/Mesh.hxx"
#include "mesh/QuadMesh.hxx"
#include "mesh/TMOP.hxx"
#include "mesh/Pillow.hxx"
#include "ATLAS/ATLAS.hxx"
#include "ATLAS/BlockMesh.hxx"
#include "MERIDIAN/Arrangement.hxx"
#include "MERIDIAN/ConeCut.hxx"
#include "MERIDIAN/MERIDIAN.hxx"
#include "MERIDIAN/ConeSingularities.hxx"
#include "MERIDIAN/DiskTemplate.hxx"
#include "MERIDIAN/Immersion.hxx"
#include "MERIDIAN/Interfaces.hxx"
#include "MERIDIAN/LayoutEnergy.hxx"
#include "MERIDIAN/RicciFlow.hxx"
#include "MERIDIAN/Separatrices.hxx"
#include "MERIDIAN/QuadMesh.hxx"
#include "MERIDIAN/SplineFit.hxx"
#include "MERIDIAN/SubdomainLabels.hxx"
#include "Parameterization/CutMesh.hxx"
#include "Parameterization/HarmonicCut.hxx"
#include "Parameterization/UVGParam.hxx"
#include "polyvector/PolyVectors.hxx"
#include "crossfield/CrossField.hxx"
#include "dualMBO/DualMBO.hxx"
#include "ZIPLINE/LayoutBlocks.hxx"
#include "ZIPLINE/PartitionSimplify.hxx"
#include "ZIPLINE/QuadLayout.hxx"
#include "ZIPLINE/SeparatrixTrace.hxx"
#include "ZIPLINE/ZIPLINE.hxx"
#include "medialaxis/MedialAxis.hxx"
#include "medialaxis/MedialAxisMap.hxx"
#include "medialaxis/MedialAxisTMesh.hxx"
#include "quantization/QuantTMeshConvert.hxx"
#include "quantization/TMeshContract.hxx"
#include "quantization/TMeshQuantizer.hxx"
#include "OASIS/OASIS.hxx"
#include "UMBER/BlockLayout.hxx"
#include "UMBER/MotorcycleGraph.hxx"
#include "UMBER/UMBER.hxx"
#include "TORSION/ConeMetric.hxx"
#include "TORSION/FieldFrames.hxx"
#include "TORSION/TORSION.hxx"
#include "TORSION/FieldIntegration.hxx"
#include "TORSION/TutteEmbedding.hxx"
#include "TORSION/MaterialLayout.hxx"

// ── Enumerations mirroring the original viewer state machine ──────────────────

enum class Mode {
    Unselected = 0,
    PolyVector  = 1,
    ZIPLINE     = 2,
    MedialAxis  = 3,
    TORSION     = 4,
    OASIS       = 5,
    UMBER       = 6,
    MERIDIAN    = 7,
    ATLAS       = 8,
};

enum class Phase {
    MeshOnly     = 1,
    CrossField   = 2,
    Singularities = 3,
    CutSeams     = 4,
};

// Mode 2, ZIPLINE (src/ZIPLINE), is Viertel, Osting and Staten (IMR 2019): the
// P1 MBO cross field, its separatrices (Sec. 3), the quad layout with
// T-junctions they cut out, and the partition simplification of Sec. 4.
//
// Quantize and Quantized run the Campen et al. 2015 quantizer on that T-layout
// (a QuadLayout's faces are its blocks, T-junctions and all, so this is the
// same solve the MedialAxisPhase stages below run, on the T-mesh a QuadLayout
// converts to directly). Split the same way, so the console reports the solve
// before the picture that depends on it. They are a detour: the three phases
// after them do not read them.
//
// The last three finish the method the way every other mode is finished:
//
//   Blocks     the simplified layout as the shared BlockDecomposition, on
//              spline geometry (ZIPLINE/LayoutBlocks.hxx): each shared side
//              fitted once as a cubic B-spline, each block its Coons patch.
//              Drawn as every mode draws its blocks -- light-blue sides, green
//              macrovertices -- over the whole layout in grey. A component
//              with a T-junction on a side is not a block, since T-junctions
//              are not meshed yet: the junction is a red disk and the
//              component is shaded, and the shading is the part of the model
//              with no block on it, which the console prints as a coverage.
//   Mesh       mesh/BlockQuadMesh on those blocks, on the dialog and the
//              target edge length every mode shares ('e' re-opens it). The
//              refused components stay shaded and the T-junctions red.
//   Smoothed   mesh::TMOP, the same dialog and solve as the other modes' last
//              phase ('c' re-opens it), with mu at the element corners for the
//              reason ATLAS and UMBER have: a transfinite grid on a block
//              decomposition.
enum class ZIPLINEPhase {
    MeshOnly    = 1,
    CrossField  = 2,
    Stepping    = 3,
    Separatrices = 4,
    Trace       = 5,
    Layout      = 6,
    Simplified  = 7,
    Quantize    = 8,
    Quantized   = 9,
    Blocks      = 10,
    Mesh        = 11,
    Smoothed    = 12,
};

// UMBER borrows the first three DualMBO stages verbatim -- its input *is* a
// converged DualMBO cross field -- and then adds the two solves of the paper:
// the frame field of Sec. 4.2 and the polysquare of Sec. 4.3. Neither is
// animated; both are one blocking L-BFGS run with nothing worth drawing in
// between. The two middle phases split the window, model on the left and
// parameter domain on the right, the way the pipelines' Layout phase shows
// psi_R and Psi.
//
// The two block phases show the same structure twice, and the pair is the
// point. Blocks splits the window and draws the traced iso-lines in both
// domains at once, which is where a ray that went somewhere it should not have
// is visible. Decomposition takes the whole window and draws what survived
// being read as the shared BlockDecomposition
// (BlockLayout::blockDecompositionOf) -- light-blue sides and green
// macrovertices, the routine ATLAS's Blocks phase and the pipelines' Patches
// phase draw theirs with, so that the three methods' block pictures read the
// same way -- over the whole layout in grey underneath. Anywhere grey shows
// through is a face the adapter refused, and that is exactly the part of the
// model the mesh will be missing, so the two colours together are the coverage
// report the console prints in numbers.
//
// The last two phases are the pipelines' Mesh and Smoothed in all but number,
// on the same two classes:
//
//   Mesh       BlockQuadMesh (mesh/BlockQuadMesh.hxx): one count per chord,
//              each macro edge meshed once, transfinite interpolation per
//              block. Opened on a dialog as Stage 10 and ATLAS's Mesh are, and
//              taking the same target edge length, so that meshing a UMBER
//              layout and a MERIDIAN one at one number is a comparison of the
//              layouts. 'e' re-opens it.
//
//   Smoothed   mesh::TMOP, the same dialog, the same solve and the same
//              picture as the other two methods' last phase ('c' re-opens it).
//              UMBER keeps its own settings for the reason ATLAS does: these
//              are transfinite grids on a block decomposition, and mu has to
//              be sampled at the element corners on one of those or the
//              smoother hands back folds the mesh did not go in with.
enum class UMBERPhase {
    MeshOnly   = 1,
    CrossField = 2,
    Stepping   = 3,
    Frames     = 4,
    Polysquare = 5,
    Blocks     = 6,
    Decomposition = 7,
    Mesh       = 8,
    Smoothed   = 9,
};

// The phase sequence shared by the two quadrilateral-layout pipelines, MERIDIAN
// (mode 7) and TORSION (mode 4). They are Stages 1 to 10 of Shepherd, Gu and
// Hughes (2022) either way, they take the same input -- a converged DualMBO cross
// field, whose holonomy Sec. 3.1 reads the cone indices off -- and they differ
// in exactly two of the eleven phases, which is why they share one enum:
//
//                MERIDIAN                       TORSION
//   Flow         Stage 3, discrete Ricci flow   Stage 3F, comb the field and
//                                               read the matchings
//   Metric       the flat cone metric it        Stages 4F and 4R, integrate the
//                produced                       field and untangle what it
//                                               inverted -- psi_0
//
// Everything before them (the interfaces, the field, the cones, the cut) and
// everything after them (the labelling, the continuation, the separatrices, the
// arrangement, the splines, the mesh) is one body of code driven by one phase
// variable, and the mode is only asked about where those two rows differ.
//
// The stages after it are the pipeline proper, and each is chosen to show
// the thing that stage is judged on rather than just what it computed:
//
//   Cones      the cone set and the discrete Gauss-Bonnet condition of Eq. (4).
//              This is the gate: sum I(v) = 4 chi(S) is the solvability
//              condition of the Newton system four stages later, so if it does
//              not hold the flat metric being asked for does not exist and
//              nothing downstream means anything.
//
//   Cut        the cutting graph G and the disk it leaves. The two kinds of arc
//              are drawn apart because they come from different places -- the
//              void arcs from HarmonicCut, the cone arcs from Sec. 3.2.2 -- and
//              behave differently at their ends.
//
//   Flow       In MERIDIAN, the conformal factor u of Eq. (8), the actual
//              unknown of the flow, drawn as a scalar field. On a planar model
//              the interior starts flat and all the curvature sits on the
//              boundary, so what u shows is the transport: the factor swells
//              around the cones that had to absorb it.
//
//              In TORSION, Stage 3F: the field combed to one branch over Omega,
//              drawn as the frame J*_t it defines on every triangle, coloured
//              by the integer a_f the comb assigned. The picture is the thing
//              the stage is judged on -- a_f constant over a patch and stepping
//              only across an arc of G is what "one branch over a disk" looks
//              like, and a step anywhere else is the combing defect the report
//              counts.
//
//   Metric     In MERIDIAN, what the flow produced. Split screen, because the
//              flat cone metric is a set of edge lengths rather than a set of
//              positions and neither half alone says what it is: on the left
//              the model with every edge coloured by how far the flow stretched
//              it, on the right each cone's one-ring unfolded *in that metric*,
//              which is the only place the cone angles themselves can be seen.
//              See viewer::ConeFan.
//
//              In TORSION, Stages 4F and 4R: psi_0 itself. Split screen for a
//              different reason -- the integration produces a map straight away
//              rather than a metric, so the right half is that map, and what it
//              is being judged on is whether it inverted anything. Red faces
//              there are the whole cost of the substitution. 'p' swaps the
//              least-squares map for the untangled one, which is the only way
//              to see what Stage 4R did.
//
//   Layout     Stages 4, 5 and 6 together -- the metric immersion psi_R, the
//              subdomain labelling, and the penalty continuation that turns the
//              first into Psi. They are one phase rather than three because
//              only the first and the last have a picture, and the first is the
//              last's starting point: psi_R and Psi are the same triangulation
//              in the same plane, and what Stage 6 did is the difference
//              between them, which is only visible if the continuation is run
//              before anything is drawn. Split screen -- the model on the left
//              and the parameter domain on the right -- because a map is a
//              thing with two ends. On the field route it is not the first
//              stage to produce one: the Metric phase already did, and this
//              phase is where that map becomes a layout.
//
//   Separatrices  Stage 7, Sec. 4. The integral curves out of the cones,
//              marched over Psi and continued across the cutting graph. Split
//              screen again, and this is the phase the split is really for:
//              Stage 7 stores each curve as barycentric coordinates in a list
//              of triangles, and those are the same numbers in the image and on
//              S, so the left half is the very same curve as the right half and
//              not a second computation of it. What each half answers is
//              different, though. The right is where the curve is straight --
//              every segment axis-parallel, which is what being an integral
//              curve of Psi means -- and where the seam jumps are, since Q4
//              moves the image to the other bank of the cut while the model
//              walks on. The left is where the layout it induces actually is:
//              this is the left half of the paper's Fig. 9 and what Stages 8
//              to 10 partition and fit.
//
//   Patches    Stages 8 and 9 together, Secs. 4 and 5. One phase rather than
//              two because neither has a picture the other does not: Stage 8
//              turns the bundle of curves into a planar subdivision of S --
//              nodes, arcs, faces -- and Stage 9 replaces each arc by the one
//              cubic B-spline both of its faces share, which moves the same
//              lines by less than the width they are drawn at. What is worth
//              seeing is the partition, and it is the same partition either
//              way. So the blocks are drawn as the paper's Fig. 12 draws them,
//              on the model and not in the image: each patch outlined in light
//              blue along its four sides, and every node of the arrangement --
//              the cones, the crossings, the boundary hits -- as a green disk.
//              A block that is missing a corner is a block whose outline runs
//              straight through a green disk without turning, which is the
//              T-junction the validation counts.
//
//              This phase is only entered when Stage 7 finished cleanly. An
//              arrangement of curves that did not close is not a layout, and
//              drawing one as though it were is the one thing this picture
//              must not do.
//
//   Smoothed   Stage 12, TMOP (mesh::TMOP). Not a stage of the paper: the mesh
//              Stage 10 hands over is transfinite interpolation and nothing
//              else, and its element quality is whatever the interval
//              assignment and the patch shapes happened to give. This phase
//              runs the node-local TMOP solve over it.
//
//              It takes the whole window rather than splitting it, and is
//              driven by a dialog like the Stage 10 phase above it and for the
//              same reason -- which metric, which exponent and how many sweeps
//              are judgements about the model that only trying a number
//              settles, so 'c' here re-opens the dialog and smooths the Stage
//              10 mesh again from scratch rather than smoothing the smoothed
//              one further.
//
//              Drawn exactly as the phase before it is, materials and folds and
//              block walls included, because the whole question is what moved:
//              two pictures that differ in how they were drawn cannot answer
//              it. The block walls come from Stage 10's own blocks and still
//              index the smoothed vertex array, so a wall is drawn where the
//              smoothing left it.
enum class PipelinePhase {
    MeshOnly     = 1,
    CrossField   = 2,
    Stepping     = 3,
    Cones        = 4,
    Cut          = 5,
    Flow         = 6,
    Metric       = 7,
    Layout       = 8,
    Separatrices = 9,
    Patches      = 10,
    Mesh         = 11,
    Smoothed     = 12,
};

// ATLAS (mode 8), docs/square_transport_2d_theory_and_implementation.md: a
// block decomposition with no cross field in it, so none of the pipelines'
// first seven phases apply. It has its own six stages up to the blocks, and
// then the pipelines' last two phases, drawn by the same code:
//
//   Domain     Stage 1, PlanarDomain: the input validated and tagged, every
//              protected vertex drawn in the colour of the quarter turns its
//              corner takes out of Sec. 4's identity -- which is the cone
//              index the pipelines would give a boundary cone there.
//
//   Field      Stage 1b, ReferenceField: the DualMBO cross field of either
//              pipeline's Stage 0, solved on the input with TORSION's
//              settings, drawn as they draw theirs -- a cross per triangle and
//              a disk at every cone. Nothing is integrated from it: it is the
//              prior the searches are scored against
//              (docs/atlas_crossfield_guidance.md), so this is the one phase
//              that shows what the layout is being asked for before any of it
//              has been decided. ATLAS::run() solves it again for itself.
//
//   Carrier    Stage 2, SquareCarrier: the three-quad split of the input,
//              valid by construction, with every vertex whose valence is not
//              the regular one drawn as a disk in the cones' colours. On any
//              real mesh that is thousands of them, and the picture is the
//              argument for everything after it (CoarseDomain's header).
//
//   Search     Stages 3 to 6: ATLAS::run(), every search in parallel threads,
//              blocking and announced a frame ahead. What is drawn is the
//              search that won, on the carrier it searched -- for a coarse
//              search the coarse re-triangulation of the domain, with the
//              input's own dS under it -- together with the blocks it found
//              there. 'p' swaps in the carrier the search started from, since
//              what Stages 3 and 6 did is the difference between the two.
//
//   Blocks     Stages 4 and 5 on the input: the winning layout realised on
//              the input domain and certified there, drawn as MERIDIAN's
//              Patches phase draws its layout -- block sides in light blue,
//              macrovertices as green disks, nothing of the model beneath.
//
//   Mesh       BlockMesh: Sec. 11.2's counts and transfinite interpolation per
//              block, opened on a dialog as Stage 10 is ('e' re-opens it).
//
//   Smoothed   mesh::TMOP, the same dialog, the same solve and the same
//              picture as the pipelines' Stage 12 ('c' re-opens it).
enum class ATLASPhase {
    MeshOnly = 1,
    Domain   = 2,
    Field    = 3,
    Carrier  = 4,
    Search   = 5,
    Blocks   = 6,
    Mesh     = 7,
    Smoothed = 8,
};

// OASIS is a one-shot solve driven by a parameter dialog rather than a
// sequence of stages, so it has only "before" and "after".
enum class OASISPhase {
    MeshOnly = 1,
    Field    = 2,
};

// The last two stages are the quantization of Campen et al. 2015: Quantize
// turns the block decomposition into the T-mesh consistency system and
// solves it for integer edge lengths, and Quantized draws the quad grid
// those lengths prescribe. They are split so the console reports the solve
// before the picture that depends on it.
enum class MedialAxisPhase {
    MeshOnly     = 0,
    DelaunayMesh = 1,
    MedialAxis   = 2,
    Map          = 3,
    TMesh        = 4,
    Quantize     = 5,
    Quantized    = 6,
};

// ── Widget ────────────────────────────────────────────────────────────────────

class CrossGenWidget : public QOpenGLWidget {
    Q_OBJECT

public:
    explicit CrossGenWidget(const std::string &meshPath, QWidget *parent = nullptr);
    ~CrossGenWidget() override = default;

protected:
    // QOpenGLWidget interface
    void initializeGL()           override;
    void resizeGL(int w, int h)   override;
    void paintGL()                override;

    // Input events
    void keyPressEvent(QKeyEvent   *event) override;
    void mousePressEvent(QMouseEvent  *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent   *event) override;
    void wheelEvent(QWheelEvent *event)        override;

private slots:
    void onTimer();

private:
    // ── helpers ──────────────────────────────────────────────────────────────
    void doReset();
    void advancePhase();

    // Modal dialog collecting the OASIS parameters. Returns false if cancelled.
    bool promptOASISParameters();

    // Assemble and solve the Eq. 13 KKT system at the current oasisLambda_,
    // with the orientation term of Sec. 5.1 when a guiding field was asked for.
    void runOASIS();

    // Run MBO to convergence (or oasisMBOIterations_) and leave the result in
    // oasisGuide_, to be used as the guiding direction field. Returns false if
    // the solve failed, in which case oasisGuide_ is cleared.
    bool buildOASISGuidingField();

    // Cut the mesh with HarmonicCut (Sec. 4.1) and optimize Eq. (1) on top of
    // the DualMBO field with L-BFGS. Blocking, like runOASIS: there is nothing to
    // draw between the continuation stages.
    void runUMBER();

    // Deform the cut mesh into the parameter domain under the optimized frame
    // field (Sec. 4.3). Blocking for the same reason.
    void runPolysquare();

    // Trace the iso-lines and label the blocks (Sec. 5). Fast enough not to
    // need the announce-a-frame-ahead treatment the two solves get.
    void runBlocks();

    // The blocks as a graph rather than a colouring of the triangles: the
    // nodes, the arcs between them and the faces they bound. Idempotent, and a
    // no-op until runBlocks() has produced something to build from.
    void buildBlockLayout();

    // The traced layout read as the shared BlockDecomposition, with a note of
    // what became of its faces. Called from the one place that can change the
    // structure, so that nothing downstream ever holds a decomposition of a
    // layout that has since been rebuilt.
    void buildUMBERDecomposition();

    // The UMBER mesh dialog: Stage 10's, taking the same target edge length
    // and sharing meshSettings_ with it, so that the three methods' meshes are
    // asked for in one number. 'e' at the Mesh phase re-opens it, and
    // cancelling leaves whatever mesh is already there.
    bool promptUMBERMesh();

    // BlockQuadMesh on the current decomposition at those settings, reported
    // line for line as ATLAS's and Stage 10's meshes are.
    void runUMBERMesh();
    // The mesh dialog and the mesh report both block-decomposition modes
    // (UMBER and ZIPLINE) share.
    bool promptBlockQuadMesh(const BlockDecomposition &decomp, const char *title);
    void reportBlockQuadMesh(const BlockQuadMesh &bqm, double ms);

    // The mesh with the optimized frame, its cuts and its boundary corners --
    // the left half of the split screen, and the whole of the Frames phase.
    void renderUMBERField();

    // Stage 0b: the material interface network, read off the triangle tags.
    // It needs nothing but the mesh -- no field, no cones -- so it runs as soon
    // as MERIDIAN mode is chosen, and its picture is under every phase from the
    // first. On a single-material mesh it finds nothing and says so once.
    void runMERIDIANInterfaces();

    // Stage 1 of Shepherd et al.: read the cone indices off the DualMBO field,
    // check Eq. (4), and rebalance the boundary cones if it does not hold.
    // Cheap; unlike the two below it needs no announcement.
    void runMERIDIANCones();
    // Advance the DualMBO solve by at most `budget` MBO steps, level by level
    // of TORSION's tau-continuation in TORSION mode, setting dualMBOConverged_
    // once the last level stops. The stepping phase calls it a few steps a
    // frame; a stage that needs the finished field calls it with no limit.
    void stepDualMBOField(int budget);

    // Stage 2: HarmonicCut's void arcs plus the Sec. 3.2.2 cone arcs, and the
    // disk they cut S into.
    void runMERIDIANCut();

    // TORSION Stage 3F (docs/cf_flow_pipeline.md Sec. 5), which sits in the
    // phase Pipeline A runs the flow in: comb the DualMBO field to one branch over
    // Omega, read the matchings off it, and audit the indices they imply
    // against the ones the field itself read. Cheap -- one BFS and one pass
    // over the vertices -- so unlike the flow it needs no announcement.
    void runTORSIONFrames();

    // TORSION Stages 4F, 4R and 4 (Secs. 6, 7.2 and 6.2), in the phase Pipeline
    // A shows the flat metric in, and for the same reason the two phases are
    // one on that side: the integration produces a *map*, so what would have
    // been "the metric the flow reached" is here "the map the field integrates
    // to", and the immersion that wraps it is the same object either route
    // ends at. Blocking, and announced a frame ahead.
    void runTORSIONIntegration();

    // TORSION's per-material mode (MaterialLayout, docs/cf_flow_pipeline.md
    // Sec. 15), the default on a multi-material model: every material region
    // laid out by Stages 1 to 8 on its own, matched across the interfaces and
    // glued. It is one blocking step here, run at the cone phase, and the
    // phases after it show what it produced -- each region's cones, cuts and
    // frames on the model, one region's psi_0 and Psi on the right ('[' and ']'
    // step through them), every region's separatrices with the layout vertices
    // the matching joined -- until the patch phase, which fits Stage 9 to the
    // glued arrangement exactly as it fits it to a traced one.
    void runTORSIONPerMaterial();
    void runTORSIONGluedPatches();
    // Whether TORSION is running per material on this model: the mode, the
    // switch ('w' toggles it before the cone phase), and more than one
    // material to run it over.
    bool perMaterialView() const;
    // Fit the right half to the region shown there, at the phases that have
    // one to show; and whether this phase shows one.
    void fitRegionView();
    bool perMaterialRegionPanel() const;

    // Stage 3: the Newton solve on Eq. (10), then the two things drawn from it
    // -- the conformal factor as a scalar field and the unfolded cone fans.
    // Blocking, like runUMBER, and announced a frame ahead for the same reason.
    void runRicciFlow();

    // The model under the flat cone metric: the left half of the Metric phase,
    // and what the Cut and RicciFlow phases draw their cones and graph over.
    void renderMERIDIANModel();

    // Stages 4 to 6 in one go: unfold the flat metric into the plane, label the
    // subdomains, and run the penalty continuation of Eq. (13). Blocking and
    // announced a frame ahead, like the Ricci solve, and by some way the
    // longest of the MERIDIAN stages.
    void runMERIDIANLayout();

    // The connectivity settings of Stages 5 to 7, gathered in one dialog
    // because they are one decision made in three places.
    //
    // What the dialog is for. Q5 asks that every integral curve out of a cone
    // be finite, and E5 is what makes it so -- but only along the paths of
    // Gamma_topo it was given. A direction nothing quantised is a geodesic of a
    // flat cone metric, and it does not end: it winds. So the layout a model
    // comes out with depends on which connections were found, and *that* is a
    // judgement about the model rather than a constant of the method. The paper
    // says as much: it takes Gamma_topo as an input, calls its automatic
    // generation future work, and places the cones of its own reference figure
    // by hand.
    //
    // The defaults are the ones every model in data/meshes settles on, and the
    // dialog shows what each of them comes to in image units for *this* model,
    // which is the only form in which they can be judged.
    //
    // Opened after Stage 4 and before Stage 5, because the numbers it reports
    // -- the mean spacing of the cones, the closest pair of them, the extent of
    // the image -- are all measurements of psi_R and do not exist until the
    // immersion does. Cancelling keeps whatever was last used, so the pipeline
    // runs either way.
    bool promptMERIDIANConnectivity(const Immersion &imm, const std::vector<Point> &uv);

    // Stages 5 to 7 again from the immersion already computed, at whatever the
    // dialog was last left at. What 'c' does at the Separatrices phase: a
    // tolerance is a judgement, and the only way to settle one is to try a
    // number and look.
    void rerunMERIDIANConnectivity();

    // Stage 7, Sec. 4: march the integral curves out of every cone over Psi,
    // continuing them across the cutting graph, until each one terminates the
    // way Q5 allows. Much quicker than the continuation above it -- a fraction
    // of a second on every model in data/meshes -- but a curve that never
    // terminates runs to the step cap, so it is announced a frame ahead like
    // the two blocking stages before it.
    void runMERIDIANSeparatrices();

    // Stages 8 and 9: the arrangement of the traced curves and the bicubic
    // patches fitted to it. Runs only once Stage 7 has come back with Q5
    // verified and a valid layout -- see meridianTraceIsClean().
    void runMERIDIANPatches();
    bool meridianTraceIsClean() const;

    // Stage 10: the quadrilateral mesh itself. A target edge length is a
    // judgement about the model in the same way the connectivity tolerance is,
    // so the phase opens on a dialog and 'c' at it opens the dialog again --
    // one number, tried and looked at.
    //
    // Cancelling leaves whatever mesh is already there, so backing out of the
    // dialog is never destructive.
    bool promptMERIDIANMesh();

    // Builds the mesh at the settings the dialog was left at, with the Winslow
    // smoothing off: this is the grid transfinite interpolation gives, which is
    // the one that answers whether the interval assignment was right. The
    // smoothed mesh is a different question and TestMERIDIAN asks it.
    void runMERIDIANMesh();

    // Stage 12: the TMOP settings, in a dialog for the same reason Stage 10's
    // are. What the numbers on it come to is a question about this mesh -- how
    // many elements there are to sweep over, how many of its nodes are free to
    // move at all -- so the dialog reports those next to them, the way the
    // meshing dialog reports what a target edge length comes to on this model.
    //
    // Cancelling leaves whatever is already on screen, so backing out is never
    // destructive.
    bool promptTMOP();

    // Build a mesh::QuadMesh from Stage 10's mesh (or Stage 11's merged one)
    // and run mesh::TMOP over it at the dialog's settings, leaving the result
    // in smoothMesh_.
    //
    // Always from the Stage 10 mesh and never from the last smoothed one, so
    // that trying a second metric is a fresh attempt and not a further one.
    // This is what makes the before/after numbers in the console mean what
    // they say.
    void runTMOP();

    // Record why a stage produced nothing: to the console in full, to the
    // terminal, and to `pipelineBlocked_` as the short form the overlay keeps
    // on screen for as long as it is true. `what` is one clause, no prefix.
    void blockPipeline(const std::string &stage, const std::string &what);

    // Stage 0c, run the moment either pipeline is chosen and before anything
    // else has been built on the mesh: find the circular inclusions and, if
    // there are any and the dialog is accepted, replace mesh_ by the matrix
    // with a hole where each of them was. Every stage after it is indexed on
    // that mesh, which is why this cannot wait until the phase it belongs to --
    // there is no such phase. See DiskTemplate.
    //
    // On a mesh with no circular inclusion -- every single-material model, and
    // every multi-material one whose regions are not disks -- it finds nothing,
    // says nothing and leaves the pipeline exactly as it was.
    void runDiskExcision();

    // The dialog that decides it, opened only when there is something to
    // decide. Lists what was found, because "excise ten inclusions" is a
    // statement about this model and the radii are how it is judged. Returns
    // false if declined, and then the disks are laid out like any other region.
    bool promptDiskExcision(const std::vector<DiskTemplate::Inclusion> &found);

    // Stage 11: fill each excised rim with its O-grid. Runs straight after
    // Stage 10 and off the same rim arcs, since Stage 10 needed them first to
    // make the rims carry an even number of edges. A no-op when nothing was
    // excised.
    void runDiskTemplates(const std::vector<std::vector<int>> &rims);

    // ATLAS, Stage 1 and Stage 2 on their own: cheap, so each phase builds its
    // own and says what it found. ATLAS::run() rebuilds both internally; they
    // are milliseconds.
    void runATLASDomain();
    void runATLASCarrier();

    // Stage 1b on its own, the same way: the reference cross field, which the
    // search will solve again for itself. Not milliseconds -- 0.3 s median and
    // 1.8 s worst over the corpus -- but the phase is the only place the field
    // can be looked at, and paying for it twice keeps every ATLAS phase built
    // from its own inputs, as the pipelines' are.
    void runATLASField();

    // Stages 1 to 6, ATLAS::run(): every search, in parallel threads. Blocking,
    // and announced a frame ahead like the Ricci solve.
    void runATLAS();

    // The winning cover on the input, reported once when the Blocks phase is
    // first drawn.
    void logATLASBlocks();

    // The TFI mesh on the blocks (BlockMesh). The dialog is Stage 10's -- one
    // target edge length, tried and looked at -- with Stage 10's settings
    // shared, so a comparison with MERIDIAN or TORSION at the same target is
    // the default; the chart switch stands where Stage 10's spline one does.
    bool promptATLASMesh();
    void runATLASMesh();

    // Everything ATLAS mode draws, phase by phase.
    void renderATLAS();

    // The mesh Stage 12 smooths: ATLAS's TFI mesh in ATLAS mode, Stage 11's
    // merged mesh where there is one, Stage 10's otherwise. mesh::QuadMesh
    // takes all three the same way.
    bool haveFinishedMesh() const;
    mesh::QuadMesh finishedMesh(const mesh::QuadMesh::Options &o) const;

    // Whether a parameter domain occupies the right half of the window.
    bool inUVSplitScreen() const;

    // Projection for one half of a split screen, and the line between them.
    void applyHalfOrtho(int x, int vpW, const viewer::ViewState &vs) const;
    void drawSplitDivider(int halfW) const;

    // UMBER mode and both layout pipelines share their first three stages, so
    // the guards that drive the DualMBO solve ask about the stage rather than the
    // mode.
    bool dualMBOStageWantsField() const;
    bool dualMBOStageIsStepping() const;

    // Whether the current mode is one of the two quadrilateral-layout
    // pipelines. Nine of the eleven phases are shared between them and are
    // guarded by this rather than by either mode; the mode itself is only asked
    // about at Flow and Metric, and where a picture belongs to one route alone.
    bool inPipeline() const {
        return mode_ == Mode::MERIDIAN || mode_ == Mode::TORSION;
    }

    // Nodes of the simplified layout that are T-junctions the quantization
    // failed to resolve. A T-junction is found structurally -- three
    // interior darts at a non-singularity; the stored kind can be stale
    // after chord collapses merge nodes -- and one that welds is not
    // returned: with every incident edge quantized >= 1 the grids on both
    // sides place a tick on the junction and it becomes a regular vertex
    // of the result. What is left hanging is a junction with an incident
    // edge forced to zero, or bordering a component the conversion had to
    // skip. Empty until ziplineQuant_ exists.
    std::vector<int> hangingTJunctions() const;
    // The interface network the current mode honours, or null: ZIPLINE's own
    // in mode 2 (Stage 0b with closed loops unsplit), interfaces_ in every
    // other mode.
    const Interfaces *interfaceNetwork() const;

    // ZIPLINE's last three phases: the blocks, the mesh on them, and its dialog.
    void buildZIPLINEBlocks();
    bool promptZIPLINEMesh();
    void runZIPLINEMesh();

    // rendering sub-routines called from paintGL
    void renderMBOAnimation();
    void renderDualMBOAnimation();
    void renderTraceAnimation();
    void renderNormal();
    void renderOverlay(const char *helpText);

    // Everything paintGL does after the clear, factored out so an export can
    // draw the identical scene: the SVG capture has to run the same draw calls
    // a second time, and a raster export runs them into a different
    // framebuffer. A figure that came from a second, export-only drawing path
    // would be a picture of something that was never on screen.
    void drawScene();

    // The screen-space colour keys. Chrome, like the console and the axis, so
    // 'h' hides them too.
    void drawLegends();

    // ── figure export ────────────────────────────────────────────────────────
    // Both write into figureDir_ under a name built from the model, the mode
    // and the phase, so a walk through the pipeline pressing 's' at each stage
    // comes out as a numbered set rather than a pile of overwrites.
    QString nextFigurePath(const char *extension) const;

    // The vector one, and the one to prefer: the scene through the GL feedback
    // buffer, captured as SVG and painted onto a one-page PDF at 1 pt per
    // pixel, since a PDF is what a paper includes. Sizes the feedback buffer by
    // trying and growing, since the primitive count is not known until it has
    // been drawn.
    void exportPdf();

    // The raster fallback, at `scale` times the on-screen framebuffer. For a
    // wireframe figure the SVG is better in every way; this is here for the
    // phases whose picture is a filled field, where a vector file is enormous
    // and gains nothing.
    void exportPng(int scale);

    // per-frame computation (lazy, guarded by has_value / pointer checks)
    void runComputations();

    // convenience. exportScale_ is 1 except while a supersampled raster export
    // is running, and it multiplies here rather than at the call sites because
    // every screen-space thing the viewer draws -- the ortho box, the split
    // divider, the HUD -- is derived from these two.
    int fbw() const { return static_cast<int>(width()  * devicePixelRatio()) * exportScale_; }
    int fbh() const { return static_cast<int>(height() * devicePixelRatio()) * exportScale_; }

    // ── data ─────────────────────────────────────────────────────────────────
    std::shared_ptr<Mesh> mesh_;

    std::optional<PolyField>   field_;
    std::optional<CutMesh>     cutMesh_;
    std::optional<DualMBO>        dualMBOField_;
    // ZIPLINE, every stage of it from the interface network to the mesh on the
    // blocks, the same class and Options TestZIPLINE runs, so that a model
    // traces the same here as on the command line. Mode 2's phases run its
    // stages one at a time (the field and the trace a few steps per frame);
    // the accessors below read them, null until a stage has run. TMOP is not
    // ZIPLINE's here: Stage 12 is the one dialog and solve every mode shares.
    std::unique_ptr<ZIPLINE> zipline_;
    const CrossField *ziplineField() const {
        return zipline_ && zipline_->hasField() ? &zipline_->getField() : nullptr;
    }
    const SeparatrixTrace *ziplineTrace() const {
        return zipline_ && zipline_->hasTrace() ? &zipline_->getTrace() : nullptr;
    }
    const QuadLayout *ziplineLayout() const {
        return zipline_ && zipline_->hasLayout() ? &zipline_->getLayout() : nullptr;
    }
    // Sec. 4's pass; its getLayout() is the simplified layout.
    const PartitionSimplify *ziplineSimplified() const {
        return zipline_ && zipline_->hasSimplification() ? &zipline_->getSimplification() : nullptr;
    }
    const LayoutBlocks *ziplineBlocks() const {
        return zipline_ && zipline_->hasBlocks() ? &zipline_->getBlocks() : nullptr;
    }
    const BlockQuadMesh *ziplineMesh() const {
        return zipline_ && zipline_->hasBlockMesh() ? &zipline_->getBlockMesh() : nullptr;
    }
    // The simplified layout converted to a QuantTMesh and quantized -- the
    // same Sec. 6 solve blockQuant_ below runs, on the block decomposition
    // tracing left rather than the medial axis one. xIdeal is 1 everywhere,
    // for the same reason. A viewer-only detour, not a stage of ZIPLINE, and
    // cleared with it.
    std::optional<QuadLayoutQuant> ziplineQuant_;
    TMeshQuantizer::Report ziplineQuantReport_;
    // Triangles of the model inside a component that is not a block, found
    // once when the blocks are built: the uncovered area, filled.
    std::vector<int> ziplineUncoveredTris_;
    std::shared_ptr<Mesh>      delaunayMesh_;
    std::shared_ptr<MedialAxis> medialAxis_;
    // The Sec. 3 map phi from the boundary to the axis above; holds a
    // shared_ptr to it, so the two are cleared together in reset().
    std::optional<MedialAxisMap> medialAxisMap_;
    // The Sec. 4 coarse block decomposition cut out by the map above.
    std::optional<MedialAxisTMesh> medialAxisTMesh_;
    // The block decomposition welded into a QuantTMesh and quantized. Built
    // from medialAxisTMesh_ and cleared with it. xIdeal is 1 everywhere: the
    // target is the coarsest valid blocking, not a mesh of a given size.
    std::optional<BlockQuant> blockQuant_;
    TMeshQuantizer::Report quantReport_;
    std::optional<OASIS>       oasis_;
    // UMBER runs on the DualMBO field held in dualMBOField_, so it needs no field of
    // its own; the cuts and the optimized frames are all that is added.
    std::optional<HarmonicCut>  umberCut_;
    std::optional<UMBER>        umber_;
    std::optional<Polysquare>   polysquare_;
    std::optional<MotorcycleGraph> blocks_;
    // This points into the one before it, so blocks_ is never rebuilt without
    // clearing it first.
    std::optional<BlockLayout>     blockLayout_;
    // The structure as the shared representation. Rebuilt whenever the layout
    // is, and the one thing the Decomposition phase's picture and the Mesh
    // phase are both read from, so that the picture and the mesh cannot be of
    // two different structures.
    std::optional<BlockDecomposition> umberDecomp_;
    BlockLayout::DecompositionReport  umberDecompReport_;
    std::optional<BlockQuadMesh>      umberMesh_;
    // Cached results of the solve: recomputing them per frame would walk every
    // vertex star for nothing.
    std::vector<std::pair<int, int>>    umberCorners_;   // (vertex, quarter turns)
    std::vector<std::pair<int, double>> umberInternal_;  // what failed to reach the boundary
    // MERIDIAN runs on the DualMBO field in dualMBOField_, like UMBER. Each stage
    // holds a reference to the one before -- ConeCut and RicciFlow both read
    // cones_, and ConeCut checks it was measured on this very mesh -- so cones_
    // is never rebuilt without clearing the two below it first.
    // Stage 0b. SubdomainLabels holds a bare pointer to this, so it is declared
    // before every stage that can be handed one and therefore destroyed after
    // them.
    std::optional<Interfaces>        interfaces_;
    std::optional<ConeSingularities> cones_;
    std::optional<ConeCut>           coneCut_;
    std::optional<RicciFlow>         ricci_;
    // Derived from ricci_ once, because both are a walk over every edge or
    // every cone star and neither belongs in a paint call.
    viewer::FlatMetric               flatMetric_;
    std::vector<viewer::ConeFan>     coneFans_;
    // The conformal factor with its mean removed, ready for the diverging ramp.
    // u is only defined up to an additive constant -- pinning one vertex is
    // what fixes it at all -- so the mean is the honest zero to draw about,
    // not the pinned vertex's value.
    Eigen::VectorXd                  ricciU_;
    double                           ricciUAbsMax_ = 1.0;
    // TORSION's Stages 3F and 4F/4R, filling the same two phases Pipeline A
    // fills with the flow and the metric it produced.
    //
    // scaffold_ is an Immersion over a throwaway map, built only for the arcs
    // of G, their (e+, e-) pairing and their quarter turns -- all three come off
    // ConeCut and the frames rather than off the map -- so it is what the
    // integration's constraint rows are written against. It holds references to
    // coneCut_ and cones_ like every other Immersion here, so it is declared
    // after them and cleared before them.
    // The indices the cross field itself read, snapshotted before Stage 1's
    // prescribe() and rebalance() move any of them on purpose. Sec. 5.1's audit
    // is run against this rather than against the set as it stands, so that it
    // reports a disagreement between the matchings and the field and not the
    // pipeline doing its job.
    std::vector<int>                 fieldIndex_;
    // Sec. 4's flat cone metric: the conformal factor the frame is scaled by
    // and the reference E1 is measured against. Built with the frames, because
    // the frames want its sizing field.
    std::optional<ConeMetric>        coneMetric_;
    std::optional<FieldFrames>       frames_;
    std::optional<Immersion>         scaffold_;
    std::optional<FieldIntegration>  integration_;
    std::optional<TutteEmbedding>    tutte_;
    // psi_0 as the least-squares solve returned it, kept alongside the map that
    // survived Stage 4R so that what the substitution actually cost is on
    // screen rather than only in the report. 'p' swaps the two at the Metric
    // phase, exactly as it swaps psi_R for Psi at the Layout one.
    std::vector<Point>               integratedMap_;
    bool                             showIntegrated_ = false;
    // The alignment axis the integration was actually solved with -- Sec. 6.4's
    // in full, the strict one, or none, whichever Stage 4F's fallback kept.
    // Stage 4R reads it twice: it is what decides which coordinates a vertex
    // may move in at rung 0 of Sec. 7.2a's ladder, and it is what Sec. 6.5's
    // projection puts back afterwards. Empty means the solve was free.
    std::vector<int>                 usedAxis_;
    // Stages 4 to 6. Each holds a reference to the one before it -- Immersion
    // to the cut, the flow and the cones, SubdomainLabels to the immersion,
    // LayoutEnergy to both -- so they are destroyed in the reverse order and
    // never rebuilt without clearing the ones above them.
    std::optional<Immersion>         immersion_;
    std::optional<SubdomainLabels>   meridianLabels_;
    std::optional<LayoutEnergy>      meridianLayout_;
    // psi_R kept alongside Psi: LayoutEnergy moves its copy in place, so the
    // map Stage 6 started from is otherwise gone by the time there is anything
    // to compare it with. 'p' toggles which of the two the right panel shows.
    std::vector<Point>               psiR_;
    bool                             showPsiR_ = false;
    // Stage 7, traced on Psi. Holds a reference to the immersion like the three
    // above it, so it is cleared first and never outlives immersion_.
    std::optional<Separatrices>      separatrices_;
    // The per-material mode's regions, each a TORSION of its own, and the
    // matching between them. The glued arrangement is moved out of it into
    // arrangement_ at the patch phase; nothing else holds on to it.
    std::unique_ptr<MaterialLayout>  materialLayout_;
    bool                             perMaterial_ = TORSION::Options().perMaterial;
    int                              regionShown_ = 0;
    // Stages 8 and 9. Arrangement holds a reference to the separatrices and the
    // labels, SplineFit to the arrangement, so they are cleared before either.
    std::optional<Arrangement>       arrangement_;
    std::optional<SplineFit>         splines_;
    // Stage 10, holding a reference to the fit, so it is cleared before it.
    std::optional<QuadMesh>          quadMesh_;
    // Stages 0c and 11. The inclusions are found before any stage runs and
    // outlive all of them -- Stage 10 wants their rims and Stage 11 wants their
    // circles -- and inputMesh_ is the mesh as it was loaded, kept only so that
    // what was taken out can be said in the report. diskFill_ holds copies of
    // the arrays it was handed rather than a reference to quadMesh_, but it is
    // meaningless without it and is cleared with it.
    std::vector<DiskTemplate::Inclusion> inclusions_;
    std::shared_ptr<Mesh>                inputMesh_;
    std::optional<DiskTemplate>          diskFill_;
    // Stage 12. A copy of whichever of the two above is the finished mesh, with
    // its nodes moved, so the phase before it still has its own picture to
    // draw; cleared whenever either of them is rebuilt.
    std::optional<mesh::QuadMesh>        smoothMesh_;

    // ATLAS. The domain is built on mesh_ and the carrier on the domain, so
    // they are declared in that order and destroyed the other way round.
    // atlas_ owns its own copies of both and everything after them; the mesh
    // copies what it needs from the chosen cover and holds nothing of atlas_,
    // but it is meaningless without it and is cleared with it.
    std::unique_ptr<PlanarDomain>  atlasDomain_;
    // Stage 1b, between the two: it reads the domain's interface edges, and
    // nothing after it in the viewer reads it -- the searches inside atlas_
    // score against atlas_'s own copy, not this one.
    std::unique_ptr<ReferenceField> atlasField_;
    std::unique_ptr<SquareCarrier> atlasCarrier_;
    std::unique_ptr<ATLAS>         atlas_;
    std::optional<BlockMesh>       atlasMesh_;
    // 'p' at the Search phase: the carrier the winning search started from
    // rather than the one it ended on.
    bool atlasShowInitial_ = false;
    // Whether the input carries more than one material: what the 'm' key and
    // the material fill are offered on, since there is no Stage 0b here.
    bool atlasMultiMaterial_ = false;

    // Guiding field for the OASIS orientation term. Held by shared_ptr because
    // OASIS keeps a reference to it for as long as it lives; separate from
    // ZIPLINE's field (zipline_), which follows its own state machine.
    std::shared_ptr<CrossField> oasisGuide_;

    // ── state machine ────────────────────────────────────────────────────────
    Mode           mode_     = Mode::Unselected;
    Phase          phase_    = Phase::MeshOnly;
    ZIPLINEPhase   ziplinePhase_ = ZIPLINEPhase::MeshOnly;
    MedialAxisPhase maPhase_ = MedialAxisPhase::MeshOnly;
    OASISPhase     oasisPhase_ = OASISPhase::MeshOnly;
    UMBERPhase     umberPhase_ = UMBERPhase::MeshOnly;
    PipelinePhase  pipePhase_ = PipelinePhase::MeshOnly;
    ATLASPhase     atlasPhase_ = ATLASPhase::MeshOnly;

    // OASIS parameters and derived display range.
    double oasisLambda_  = 0.0;   // set by the dialog on first use
    double oasisAbsMax_  = 1.0;   // max|f|, the symmetric range for the ramp
    // 0 disables the Sec. 3.4 pass. These two doubles as the dialog's initial
    // state, so a nonzero value is what makes its checkbox start ticked.
    int    oasisVibrationIterations_ = 10;
    double vibrationBefore_ = -1.0;        // mean E_a before the pass, for the log

    // Orientation control (Sec. 5.1). gamma <= 0 disables it, and then no
    // guiding field is computed at all.
    double oasisOrientationWeight_ = 10.0;
    // Iteration cap for the MBO solve that produces the guiding field. MBO
    // stops early on convergence, so this only bounds the wait.
    int    oasisMBOIterations_ = 100;
    // How far the guiding field is kept clear of the boundary, in quad cells.
    // 0 guides everywhere, which the boundary conditions will fight; see the
    // note on setOrientationWeight().
    double oasisGuideClearanceQuads_ = 2.0;

    bool singularitiesLogged_  = false;
    // ZIPLINE's two animations, each announced once when it starts and once
    // when it ends; the counts themselves are zipline_->getStatus()'s.
    bool mboSteppingStarted_   = false;
    bool ziplineTracingStarted_    = false;
    bool ziplineTracingFinished_   = false;
    bool dualMBOSteppingStarted_  = false;
    bool dualMBOConverged_        = false;
    int  dualMBOStepCount_        = 0;
    // TORSION's tau-continuation (TORSION::fieldTauLadder): the tau scales
    // still to run, which of them the field is at, and the steps taken at it.
    // A single level, the heuristic tau, in every other mode.
    std::vector<double> dualMBOLadder_{1.0};
    int  dualMBOLevel_            = 0;
    int  dualMBOLevelSteps_       = 0;
    // The Eq. (1) solve is attempted once per run: a failure leaves umber_
    // empty, and retrying it every frame would only stall the viewer again.
    // It is announced one frame ahead so the notice is on screen while the
    // GUI thread is inside L-BFGS.
    bool umberAnnounced_       = false;
    bool umberAttempted_       = false;
    bool polysquareAnnounced_  = false;
    bool polysquareAttempted_  = false;
    bool blocksAttempted_      = false;
    bool umberMeshAttempted_   = false;
    bool ziplineMeshAttempted_   = false;
    // One-shot discipline for the three MERIDIAN stages. Each is attempted once
    // per run and not retried: a failure leaves its optional empty, and keying
    // off the optional alone would run the whole stage again on every frame --
    // which for the cones means re-running the DualMBO solve sixty times a second.
    // The Ricci solve is additionally announced a frame early, so the notice is
    // on screen while the GUI thread is inside the Newton loop.
    bool interfacesAttempted_  = false;
    bool conesAttempted_       = false;
    bool cutAttempted_         = false;
    bool ricciAnnounced_       = false;
    bool ricciAttempted_       = false;
    // The same discipline for Pipeline B's two. The combing is cheap and needs
    // no announcement; the integration is a sparse saddle solve followed, when
    // it inverted anything, by a whole continuation of its own, so it is
    // announced a frame ahead like the Ricci solve.
    bool framesAttempted_      = false;
    bool integrationAnnounced_ = false;
    bool integrationAttempted_ = false;
    bool materialAnnounced_    = false;
    bool materialAttempted_    = false;
    bool layoutAnnounced_      = false;
    bool layoutAttempted_      = false;
    bool separatricesAnnounced_ = false;
    bool separatricesAttempted_ = false;
    bool patchesAnnounced_      = false;
    bool patchesAttempted_      = false;
    bool meshAttempted_         = false;
    bool tmopAttempted_         = false;
    // MERIDIAN::Options::topoNearMissRetry, in the viewer. Stages 5 to 7 are
    // already re-runnable from the connectivity dialog, so the retry is that
    // path with the number filled in: once per run, and only when Stage 8 came
    // back with a piece of S no grid covers. The two saved numbers are what the
    // first attempt left, so the second can be said against it rather than
    // reported on its own.
    bool   nearMissRetried_   = false;
    double unmeshableBefore_  = -1.0;
    double nearMissBefore_    = 0.0;
    // Stage 0c is attempted once per run, at the moment the mode is chosen.
    bool disksAttempted_        = false;
    // ATLAS's stages, with the same one-shot discipline. run() is announced a
    // frame ahead; the mesh dialog, like Stage 10's, counts as asked once it
    // has opened.
    bool atlasDomainAttempted_  = false;
    // The field solve blocks for up to a couple of seconds, so it is announced
    // a frame ahead as run() is, and built on the frame after.
    bool atlasFieldAnnounced_   = false;
    bool atlasFieldAttempted_   = false;
    bool atlasCarrierAttempted_ = false;
    bool atlasAnnounced_        = false;
    bool atlasAttempted_        = false;
    bool atlasBlocksLogged_     = false;
    bool atlasMeshAttempted_    = false;

    // Why the last phases have nothing to show, in one clause, or empty when
    // nothing is wrong. Kept on screen rather than only in the console: the
    // console holds eight lines and a stage that refuses prints its reason
    // among a dozen others, so by the time the phase that is blank is on
    // screen the reason has scrolled off it. Set by whichever stage refused or
    // failed, cleared by the next one that runs.
    std::string pipelineBlocked_;

    // Stage 0b runs in two halves -- the network before Stage 1, the region
    // balance after it, because the balance needs Stage 1's cones on dS -- and
    // both append to the one message list, so this is how far it has been
    // drained into the console.
    size_t interfaceMessagesSeen_ = 0;

    // What the interface network's two pictures are showing. The network itself
    // is on whenever there is one: it is the input the multi-material path is
    // about, and every stage after it is to be judged against it. The material
    // fill is off, because it competes with the scalar ramps of Stages 3 and 4
    // and because on a single-material model it says nothing at all.
    bool showInterfaces_    = true;
    bool showMaterialFill_  = false;

    // ── figure state ─────────────────────────────────────────────────────────
    // The console and the key-help line, which are the two things on screen
    // that belong to the viewer rather than to the model. On by default and off
    // in a figure: a paper wants the picture, not the legend of shortcuts that
    // produced it. 'h' toggles them, and the exports honour whatever it is set
    // to, so a figure with the console left on is a deliberate one.
    bool showHUD_ = true;

    // Where 's' and 'S' write, and the counter that keeps a walk through the
    // pipeline from overwriting itself.
    QString figureDir_;
    std::string modelName_;
    int figureCounter_ = 0;

    // 1 on screen; the supersampling factor while a raster export is drawing.
    // fbw()/fbh() multiply by it, and viewer::setRenderScale() matches it so
    // line widths and the bitmap font keep their on-screen weight.
    int exportScale_ = 1;

    // Stage 10 settings, surviving a reset the way the connectivity ones do so
    // that the dialog opens on whatever was last tried. The defaults are
    // QuadMesh's own except for the smoothing, which the viewer never runs.
    struct MERIDIANMeshSettings {
        double target      = 0.05;
        int    minEdges    = 1;
        int    maxEdges    = 0;      // 0 = no ceiling
        bool   useSplines  = true;
        // QuadMesh::Options::collapseSpan. Its own default, because a layout
        // finer than the elements asked for is the common case rather than the
        // exception and the contraction is what the target edge length means
        // there; 0 in the dialog turns it off and shows what it was doing.
        double collapseSpan = QuadMesh::Options().collapseSpan;
    };
    MERIDIANMeshSettings meshSettings_;
    // ATLAS's one departure from those: whether the interior nodes come
    // through the block's certified chart (BlockMesh::Options::useChart) or
    // from the Coons blend of its sides -- the counterpart of useSplines.
    bool atlasUseChart_ = true;

    // Stages 0c and 11, surviving a reset like every other judgement in this
    // widget so that the next run opens on whatever was last tried. `excise` is
    // what the Stage 0c dialog was last left at and is only ever asked about on
    // a model that has circular inclusions at all; the other three are
    // DiskTemplate::Options' own defaults and are offered in the Stage 10
    // dialog, because Stage 11 is rerun with the mesh and trying a squareness
    // means meshing again.
    struct DiskTemplateSettings {
        bool   excise      = true;
        double squareness  = DiskTemplate::Options().coreSquareness;
        int    ringDepth   = DiskTemplate::Options().ringDepth;
        int    smoothing   = DiskTemplate::Options().smoothingPasses;
    };
    DiskTemplateSettings diskSettings_;

    // Stage 12's settings, surviving a reset like the two above. Every default
    // is mesh::TMOP::Options' own -- they were measured over the corpus and
    // there is nothing about running in a window that changes what they should
    // be -- except the sweep cap, which is raised because this phase is a
    // picture rather than a test and 200 sweeps stops most models short of
    // where they were going, and the metric, which defaults to 007 (shape and
    // size) here since the viewer is where a per-model size field is actually
    // visible.
    struct TMOPSettings {
        int    metric      = mesh::TMOP::ShapeSize007;
        double gamma       = mesh::TMOP::Options().gamma;
        double exponent    = mesh::TMOP::Options().exponent;
        int    target      = mesh::TMOP::Options().target;
        double targetSize  = mesh::TMOP::Options().targetSize;  // 0 = the mean edge
        bool   corners     = false;   // quadrature: 2x2 Gauss, or the corners
        int    sweeps      = 1000;
        bool   untangle    = mesh::TMOP::Options().untangle;
        int    threads     = mesh::TMOP::Options().threads;     // 0 = the runtime's own
        // QuadMesh::Options, not TMOP::Options: they decide which nodes are
        // allowed to move before the smoother is handed the mesh at all, which
        // is why changing either of them rebuilds it.
        bool   pinFeatures = mesh::QuadMesh::Options().fixAllFeatureNodes;
        double cornerAngle = mesh::QuadMesh::Options().cornerAngle;
        int    curveSource = mesh::QuadMesh::Options().curveSource;
        // Pillow the flat feature corners before smoothing (mesh::Pillow):
        // one layer of quads along each stretch of dS or interface on which an
        // element spans 180 degrees, so that TMOP is handed a corner it can
        // move. Off here and on for the pipelines' Stage 12 below.
        bool   pillow      = false;
    };
    // The pipelines' Stage 12 samples mu at the element corners too, since
    // 2026-10-02. At the 2x2 Gauss points a corner can turn over without the
    // barrier seeing it, and on the pipelines' own meshes it did: over the
    // multimat corpus at 0.05, per material, TMOP left 25 elements inverted
    // that it had made itself or failed to clear (multimat/tooth from a worst
    // of 0.00 to -0.96, artery 0.06 to -0.58), and sampled at the corners 4;
    // 16 -> 5 on singlemat through TORSION and 7 -> 0 through MERIDIAN, with no
    // model's worst element going down by more than 0.03.
    //
    // It pillows first, too, since the same day: a block corner the layout put
    // on a straight feature is an element spanning 180 degrees at a feature
    // node, which no sampling repairs because no admissible move changes the
    // angle. Over the multimat corpus at 0.05, per material, that was basin's
    // worst element going 0.000 -> 0.515 and artery's 0.136 -> 0.256 with the
    // layer in, and no model worse; docs/cf_flow_pipeline.md Sec. 15.10.
    TMOPSettings tmopSettings_ = [] {
        TMOPSettings t;
        t.corners = true;
        t.pillow = true;
        return t;
    }();
    // ATLAS's copy, which samples mu at the element corners for the same
    // reason: those are the Jacobians Sec. 9.1 of the square-transport spec
    // judges an element by -- measured on ATLAS's TFI meshes, where six corpus
    // models came out of TMOP with folds they did not go in with, and none did
    // sampled at the corners.
    TMOPSettings atlasTmopSettings_ = [] { TMOPSettings t; t.corners = true; return t; }();
    // UMBER's copy, which differs the same way and for the same reason: what
    // it smooths is a transfinite grid on a block decomposition, exactly what
    // ATLAS's is, and measured on data/meshes the corner sampling is the
    // difference between the smoother improving the worst element and making
    // it worse -- on singlemat/geom012 the worst scaled Jacobian goes 0.033 ->
    // 0.162 at the corners and 0.033 -> 0.012 at the 2x2 Gauss points.
    TMOPSettings umberTmopSettings_ = [] { TMOPSettings t; t.corners = true; return t; }();
    // ZIPLINE's, for the same reason again: its mesh is BlockQuadMesh's
    // transfinite grid, exactly UMBER's kind of mesh.
    TMOPSettings ziplineTmopSettings_ = [] { TMOPSettings t; t.corners = true; return t; }();
    // The copy above that the current mode's TMOP dialog and solve both use.
    TMOPSettings &tmopSettingsForMode();

    // The connectivity settings of Stages 5 to 7, likewise surviving a reset.
    // Defaults are the library's own, so the dialog opens on the recommended
    // values and this struct only records departures from them.
    struct MERIDIANConnectivity {
        SubdomainLabels::Options labels;
        Separatrices::Options    trace;
        int    repairPasses      = MERIDIAN::Options().repairPasses;
        int    repairMaxPerPass  = MERIDIAN::Options().repairMaxPerPass;
        double repairGapLimit    = MERIDIAN::Options().repairGapLimit;
        double repairLambdaBoost = MERIDIAN::Options().repairLambdaBoost;
        int    repairOuterSteps  = MERIDIAN::Options().repairOuterSteps;
    };
    MERIDIANConnectivity meridianConn_;
    // Whether the dialog has been shown this run. It opens once, on the way
    // into the Layout phase, and after that only when asked for.
    bool meridianConnPrompted_ = false;

    // ── view / camera ────────────────────────────────────────────────────────
    viewer::ViewState view_;     // mesh-space view (left panel)
    viewer::ViewState uvView_;   // UV-space view   (right panel, split screen)
    viewer::Bounds    bounds_;
    double avgEdge_ = 1.0;
    double scale_   = 1.0;
    double pad_     = 0.0;

    // ── mouse drag state ─────────────────────────────────────────────────────
    bool   rightDragging_ = false;
    QPoint lastMousePos_;

    // ── console / log ────────────────────────────────────────────────────────
    viewer::Console console_;

    // ── timer driving animation frames ───────────────────────────────────────
    QTimer *timer_ = nullptr;

    // The MBO convergence test both pipelines use, as `error < 2 N tol` with N
    // the triangle count. 1e-5 is what TORSION::runField() and MERIDIAN's
    // Stage 0 ship with, and it is what the paper_tests layout runs (E4, E5
    // Part 3) measure the pipeline at; a tighter one here would mean the field
    // on screen is not the field either of them lays out.
    static constexpr double DUALMBO_TOL = 1e-5;

    // L-BFGS iteration cap per continuation stage of Eq. (1). The paper stops
    // on the gradient tolerance and so does every mesh in data/meshes at a few
    // hundred to a couple of thousand iterations, so this is headroom rather
    // than a target; the console reports the count so a run that hits it is
    // visible.
    static constexpr int UMBER_LBFGS_ITERATIONS = 3000;
};
