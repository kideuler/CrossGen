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
#include "OASIS/OASIS.hxx"
#include "UMBER/BlockLayout.hxx"
#include "UMBER/ChordCollapse.hxx"
#include "UMBER/MotorcycleGraph.hxx"
#include "UMBER/UMBER.hxx"

// ── Enumerations mirroring the original viewer state machine ──────────────────

enum class Mode {
    Unselected = 0,
    PolyVector  = 1,
    MBO         = 2,
    MedialAxis  = 3,
    SIPG        = 4,
    OASIS       = 5,
    UMBER       = 6,
};

enum class Phase {
    MeshOnly     = 1,
    CrossField   = 2,
    Singularities = 3,
    CutSeams     = 4,
    UVMesh       = 5,
};

enum class MBOPhase {
    MeshOnly    = 1,
    CrossField  = 2,
    Stepping    = 3,
    Separatrices = 4,
    Trace       = 5,
    Layout      = 6,
    Simplified  = 7,
};

enum class SIPGPhase {
    MeshOnly   = 1,
    CrossField = 2,
    Stepping   = 3,
    CutSeams   = 4,
    UVMesh     = 5,
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

// OASIS is a one-shot solve driven by a parameter dialog rather than a
// sequence of stages, so it has only "before" and "after".
enum class OASISPhase {
    MeshOnly = 1,
    Field    = 2,
};

enum class MedialAxisPhase {
    MeshOnly     = 0,
    DelaunayMesh = 1,
    MedialAxis   = 2,
    Map          = 3,
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

    // Whether a parameter domain occupies the right half of the window.
    bool inUVSplitScreen() const;

    // Projection for one half of a split screen, and the line between them.
    void applyHalfOrtho(int x, int vpW, const viewer::ViewState &vs) const;
    void drawSplitDivider(int halfW) const;

    // SIPG mode and UMBER mode share their first three stages, so the guards
    // that drive the SIPG solve ask about the stage rather than the mode.
    bool sipgStageWantsField() const;
    bool sipgStageIsStepping() const;

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
    std::optional<CutMesh>     sipgCutMesh_;
    std::optional<UVGParam>    sipgUVParam_;
    std::shared_ptr<SeparatrixTrace> separatrixTrace_;
    // Holds a pointer to the trace above, so it must not outlive it: both are
    // cleared together in reset().
    std::optional<QuadLayout>  quadLayout_;
    // Sec. 4 run on the layout above. It keeps a copy of what it was handed, so
    // the layout before simplification survives alongside it and the two can be
    // drawn together.
    std::optional<PartitionSimplify> simplified_;
    std::shared_ptr<Mesh>      delaunayMesh_;
    std::shared_ptr<MedialAxis> medialAxis_;
    // The Sec. 3 map phi from the boundary to the axis above; holds a
    // shared_ptr to it, so the two are cleared together in reset().
    std::optional<MedialAxisMap> medialAxisMap_;
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
    // Guiding field for the OASIS orientation term. Held by shared_ptr because
    // OASIS keeps a reference to it for as long as it lives; separate from
    // crossField_, which belongs to MBO mode and follows its own state machine.
    std::shared_ptr<CrossField> oasisGuide_;

    // ── state machine ────────────────────────────────────────────────────────
    Mode           mode_     = Mode::Unselected;
    Phase          phase_    = Phase::MeshOnly;
    MBOPhase       mboPhase_ = MBOPhase::MeshOnly;
    SIPGPhase      sipgPhase_ = SIPGPhase::MeshOnly;
    MedialAxisPhase maPhase_ = MedialAxisPhase::MeshOnly;
    OASISPhase     oasisPhase_ = OASISPhase::MeshOnly;
    UMBERPhase     umberPhase_ = UMBERPhase::MeshOnly;

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

    // Chord collapse settings, surviving a reset the way oasisLambda_ does so
    // that the dialog opens on whatever was tried last.
    ChordCollapse::Settings chordSettings_;

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
