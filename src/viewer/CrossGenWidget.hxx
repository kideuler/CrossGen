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
#include "Parameterization/MIQ.hxx"
#include "Parameterization/UVGParam.hxx"
#include "polyvector/PolyVectors.hxx"
#include "crossfield/CrossField.hxx"
#include "sipg/SIPG.hxx"
#include "tracing/SeparatrixTrace.hxx"
#include "medialaxis/MedialAxis.hxx"
#include "OASIS/OASIS.hxx"

// ── Enumerations mirroring the original viewer state machine ──────────────────

enum class Mode {
    Unselected = 0,
    PolyVector  = 1,
    MBO         = 2,
    MedialAxis  = 3,
    SIPG        = 4,
    OASIS       = 5,
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
};

enum class SIPGPhase {
    MeshOnly   = 1,
    CrossField = 2,
    Stepping   = 3,
    CutSeams   = 4,
    UVMesh     = 5,
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
    Classify     = 3,
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
    std::shared_ptr<Mesh>      delaunayMesh_;
    std::shared_ptr<MedialAxis> medialAxis_;
    std::optional<OASIS>       oasis_;
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
};
