#pragma once

#include <QOpenGLWidget>
#include <QTimer>
#include <QPoint>

#include <memory>
#include <optional>
#include <string>

#include "viewer/ViewerTypes.hxx"
#include "viewer/Render.hxx"
#include "viewer/Geometry.hxx"

#include "mesh/Mesh.hxx"
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
#include "Parameterization/MIQ.hxx"
#include "Parameterization/UVGParam.hxx"
#include "polyvector/PolyVectors.hxx"
#include "crossfield/CrossField.hxx"
#include "sipg/SIPG.hxx"
#include "tracing/PartitionSimplify.hxx"
#include "tracing/QuadLayout.hxx"
#include "tracing/SeparatrixTrace.hxx"
#include "medialaxis/MedialAxis.hxx"
#include "medialaxis/MedialAxisMap.hxx"
#include "medialaxis/MedialAxisTMesh.hxx"
#include "quantization/QuantTMeshConvert.hxx"
#include "quantization/TMeshContract.hxx"
#include "quantization/TMeshQuantizer.hxx"
#include "OASIS/OASIS.hxx"
#include "UMBER/BlockLayout.hxx"
#include "UMBER/ChordCollapse.hxx"
#include "UMBER/MotorcycleGraph.hxx"
#include "UMBER/UMBER.hxx"
#include "TORSION/ConeMetric.hxx"
#include "TORSION/FieldFrames.hxx"
#include "TORSION/TORSION.hxx"
#include "TORSION/FieldIntegration.hxx"
#include "TORSION/TutteEmbedding.hxx"

// ── Enumerations mirroring the original viewer state machine ──────────────────

enum class Mode {
    Unselected = 0,
    PolyVector  = 1,
    MBO         = 2,
    MedialAxis  = 3,
    TORSION     = 4,
    OASIS       = 5,
    UMBER       = 6,
    MERIDIAN    = 7,
};

enum class Phase {
    MeshOnly     = 1,
    CrossField   = 2,
    Singularities = 3,
    CutSeams     = 4,
    UVMesh       = 5,
};

// The last two stages quantize the block decomposition tracing left (Sec. 4
// of Viertel et al. is a block decomposition already -- a QuadLayout's faces
// are its blocks -- so this runs the same Campen et al. 2015 quantizer the
// MedialAxisPhase stages below do, on the T-mesh a QuadLayout converts to
// directly). Split the same way, so the console reports the solve before the
// picture that depends on it.
enum class MBOPhase {
    MeshOnly    = 1,
    CrossField  = 2,
    Stepping    = 3,
    Separatrices = 4,
    Trace       = 5,
    Layout      = 6,
    Simplified  = 7,
    Quantize    = 8,
    Quantized   = 9,
};

// UMBER borrows the first three SIPG stages verbatim -- its input *is* a
// converged SIPG cross field -- and then adds the two solves of the paper:
// the frame field of Sec. 4.2 and the polysquare of Sec. 4.3. Neither is
// animated; both are one blocking L-BFGS run with nothing worth drawing in
// between. The two middle phases split the window, model on the left and
// parameter domain on the right, the way PolyVector and SIPG show their UV
// meshes.
//
// Simplified is the odd one out in two ways. It is driven by a dialog rather
// than by pressing on, the way OASIS mode is, because how thin a chord has to
// be to be worth collapsing is a heuristic and the only way to settle it on a
// given model is to try a number and look -- so 'c' at this phase re-opens the
// dialog and runs the operation again from the structure the tracing left,
// rather than advancing anywhere. And it takes the whole window instead of
// splitting it: the collapse works in model space, and drawing the structure
// before and after over the same mesh at the same scale is what shows which
// blocks it took out.
enum class UMBERPhase {
    MeshOnly   = 1,
    CrossField = 2,
    Stepping   = 3,
    Frames     = 4,
    Polysquare = 5,
    Blocks     = 6,
    Simplified = 7,
};

// The phase sequence shared by the two quadrilateral-layout pipelines, MERIDIAN
// (mode 7) and TORSION (mode 4). They are Stages 1 to 10 of Shepherd, Gu and
// Hughes (2022) either way, they take the same input -- a converged SIPG cross
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
    // the SIPG field with L-BFGS. Blocking, like runOASIS: there is nothing to
    // draw between the continuation stages.
    void runUMBER();

    // Deform the cut mesh into the parameter domain under the optimized frame
    // field (Sec. 4.3). Blocking for the same reason.
    void runPolysquare();

    // Trace the iso-lines and label the blocks (Sec. 5). Fast enough not to
    // need the announce-a-frame-ahead treatment the two solves get.
    void runBlocks();

    // The blocks as a graph rather than a colouring of the triangles, which is
    // what a chord can be walked on. Idempotent, and a no-op until runBlocks()
    // has produced something to build from.
    void buildBlockLayout();

    // Modal dialog collecting the chord collapse settings, with a live count of
    // how many chords they would let through -- the number that decides whether
    // a threshold is the right one. Returns false if cancelled.
    bool promptChordCollapseParameters();

    // Collapse chords greedily at the current settings, always starting from
    // the structure the tracing left rather than from the last result, so that
    // trying a second threshold is a fresh attempt and not a further one.
    void runChordCollapse();

    // The mesh with the optimized frame, its cuts and its boundary corners --
    // the left half of the split screen, and the whole of the Frames phase.
    void renderUMBERField();

    // Stage 0b: the material interface network, read off the triangle tags.
    // It needs nothing but the mesh -- no field, no cones -- so it runs as soon
    // as MERIDIAN mode is chosen, and its picture is under every phase from the
    // first. On a single-material mesh it finds nothing and says so once.
    void runMERIDIANInterfaces();

    // Stage 1 of Shepherd et al.: read the cone indices off the SIPG field,
    // check Eq. (4), and rebalance the boundary cones if it does not hold.
    // Cheap; unlike the two below it needs no announcement.
    void runMERIDIANCones();

    // Stage 2: HarmonicCut's void arcs plus the Sec. 3.2.2 cone arcs, and the
    // disk they cut S into.
    void runMERIDIANCut();

    // TORSION Stage 3F (docs/cf_flow_pipeline.md Sec. 5), which sits in the
    // phase Pipeline A runs the flow in: comb the SIPG field to one branch over
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

    // Stage 4R alone. Returns the untangled map, or an empty vector when there
    // was nothing locally injective to start from -- in which case psi_0 stands
    // as the solve left it and Stage 4's gate is what stops the pipeline.
    std::vector<Point> runTORSIONUntangle();

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
    // dialog was last left at. What 'c' does at the Separatrices phase, in the
    // spirit of UMBER's Simplified phase: a tolerance is a judgement, and the
    // only way to settle one is to try a number and look.
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

    // Whether a parameter domain occupies the right half of the window.
    bool inUVSplitScreen() const;

    // Projection for one half of a split screen, and the line between them.
    void applyHalfOrtho(int x, int vpW, const viewer::ViewState &vs) const;
    void drawSplitDivider(int halfW) const;

    // UMBER mode and both layout pipelines share their first three stages, so
    // the guards that drive the SIPG solve ask about the stage rather than the
    // mode.
    bool sipgStageWantsField() const;
    bool sipgStageIsStepping() const;

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
    // skip. Empty until traceQuant_ exists.
    std::vector<int> hangingTJunctions() const;

    // rendering sub-routines called from paintGL
    void renderMBOAnimation();
    void renderSIPGAnimation();
    void renderTraceAnimation();
    void renderNormal();
    void renderOverlay(const char *helpText);

    // per-frame computation (lazy, guarded by has_value / pointer checks)
    void runComputations();

    // convenience
    int fbw() const { return static_cast<int>(width()  * devicePixelRatio()); }
    int fbh() const { return static_cast<int>(height() * devicePixelRatio()); }

    // ── data ─────────────────────────────────────────────────────────────────
    std::shared_ptr<Mesh> mesh_;

    std::optional<PolyField>   field_;
    std::optional<CutMesh>     cutMesh_;
    std::optional<MIQSolver>   miqSolver_;
    std::optional<CrossField>  crossField_;
    std::optional<SIPG>        sipgField_;
    std::shared_ptr<SeparatrixTrace> separatrixTrace_;
    // Holds a pointer to the trace above, so it must not outlive it: both are
    // cleared together in reset().
    std::optional<QuadLayout>  quadLayout_;
    // Sec. 4 run on the layout above. It keeps a copy of what it was handed, so
    // the layout before simplification survives alongside it and the two can be
    // drawn together.
    std::optional<PartitionSimplify> simplified_;
    // simplified_'s layout converted to a QuantTMesh and quantized -- the
    // same Sec. 6 solve blockQuant_ below runs, on the block decomposition
    // tracing left rather than the medial axis one. xIdeal is 1 everywhere,
    // for the same reason. Built from simplified_ and cleared with it.
    std::optional<QuadLayoutQuant> traceQuant_;
    TMeshQuantizer::Report traceQuantReport_;
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
    // UMBER runs on the SIPG field held in sipgField_, so it needs no field of
    // its own; the cuts and the optimized frames are all that is added.
    std::optional<HarmonicCut>  umberCut_;
    std::optional<UMBER>        umber_;
    std::optional<Polysquare>   polysquare_;
    std::optional<MotorcycleGraph> blocks_;
    // Both point into the one before them -- blockLayout_ into blocks_, and
    // chordCollapse_ into a copy of blockLayout_'s layout -- so blocks_ is
    // never rebuilt without clearing these first.
    std::optional<BlockLayout>     blockLayout_;
    std::optional<ChordCollapse>   chordCollapse_;
    // Cached results of the solve: recomputing them per frame would walk every
    // vertex star for nothing.
    std::vector<std::pair<int, int>>    umberCorners_;   // (vertex, quarter turns)
    std::vector<std::pair<int, double>> umberInternal_;  // what failed to reach the boundary
    // MERIDIAN runs on the SIPG field in sipgField_, like UMBER. Each stage
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

    // Guiding field for the OASIS orientation term. Held by shared_ptr because
    // OASIS keeps a reference to it for as long as it lives; separate from
    // crossField_, which belongs to MBO mode and follows its own state machine.
    std::shared_ptr<CrossField> oasisGuide_;

    // ── state machine ────────────────────────────────────────────────────────
    Mode           mode_     = Mode::Unselected;
    Phase          phase_    = Phase::MeshOnly;
    MBOPhase       mboPhase_ = MBOPhase::MeshOnly;
    MedialAxisPhase maPhase_ = MedialAxisPhase::MeshOnly;
    OASISPhase     oasisPhase_ = OASISPhase::MeshOnly;
    UMBERPhase     umberPhase_ = UMBERPhase::MeshOnly;
    PipelinePhase  pipePhase_ = PipelinePhase::MeshOnly;

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
    bool mboSteppingStarted_   = false;
    bool mboConverged_         = false;
    bool mboTracingStarted_    = false;
    bool mboTracingFinished_   = false;
    int  mboStepCount_         = 0;
    bool sipgSteppingStarted_  = false;
    bool sipgConverged_        = false;
    int  sipgStepCount_        = 0;
    // The Eq. (1) solve is attempted once per run: a failure leaves umber_
    // empty, and retrying it every frame would only stall the viewer again.
    // It is announced one frame ahead so the notice is on screen while the
    // GUI thread is inside L-BFGS.
    bool umberAnnounced_       = false;
    bool umberAttempted_       = false;
    bool polysquareAnnounced_  = false;
    bool polysquareAttempted_  = false;
    bool blocksAttempted_      = false;
    // One-shot discipline for the three MERIDIAN stages. Each is attempted once
    // per run and not retried: a failure leaves its optional empty, and keying
    // off the optional alone would run the whole stage again on every frame --
    // which for the cones means re-running the SIPG solve sixty times a second.
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
    bool layoutAnnounced_      = false;
    bool layoutAttempted_      = false;
    bool separatricesAnnounced_ = false;
    bool separatricesAttempted_ = false;
    bool patchesAnnounced_      = false;
    bool patchesAttempted_      = false;
    bool meshAttempted_         = false;
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

    // Chord collapse settings, surviving a reset the way oasisLambda_ does so
    // that the dialog opens on whatever was tried last.
    ChordCollapse::Settings chordSettings_;

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

    static constexpr int MBO_MAX_STEPS = 500;

    // L-BFGS iteration cap per continuation stage of Eq. (1). The paper stops
    // on the gradient tolerance and so does every mesh in data/meshes at a few
    // hundred to a couple of thousand iterations, so this is headroom rather
    // than a target; the console reports the count so a run that hits it is
    // visible.
    static constexpr int UMBER_LBFGS_ITERATIONS = 3000;
};
