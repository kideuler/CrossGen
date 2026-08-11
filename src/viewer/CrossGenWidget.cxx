// CrossGenWidget.cxx – Qt6 QOpenGLWidget implementation of the CrossGen viewer.
// Replaces the old GLFW-based ViewerMain render loop with Qt event-driven rendering.

#include "viewer/CrossGenWidget.hxx"
#include "viewer/GL.hxx"
#include "viewer/Interaction.hxx"
#include "viewer/Render.hxx"

#include <QCheckBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QKeyEvent>
#include <QLabel>
#include <QSpinBox>
#include <QMouseEvent>
#include <QWheelEvent>
#include <QTimer>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#include "triangle/TriangleMesher.hpp"

// ── helpers ───────────────────────────────────────────────────────────────────

namespace {

std::string formatMs(double ms) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << ms << " ms";
    return oss.str();
}

using Clock = std::chrono::high_resolution_clock;

Phase nextPhase(Phase p) {
    switch (p) {
        case Phase::MeshOnly:      return Phase::CrossField;
        case Phase::CrossField:    return Phase::Singularities;
        case Phase::Singularities: return Phase::CutSeams;
        case Phase::CutSeams:      return Phase::UVMesh;
        case Phase::UVMesh:        return Phase::UVMesh;
    }
    return Phase::UVMesh;
}

MBOPhase nextMBOPhase(MBOPhase p) {
    switch (p) {
        case MBOPhase::MeshOnly:     return MBOPhase::CrossField;
        case MBOPhase::CrossField:   return MBOPhase::Stepping;
        case MBOPhase::Stepping:     return MBOPhase::Separatrices;
        case MBOPhase::Separatrices: return MBOPhase::Trace;
        case MBOPhase::Trace:        return MBOPhase::Layout;
        case MBOPhase::Layout:       return MBOPhase::Simplified;
        case MBOPhase::Simplified:   return MBOPhase::Simplified;
    }
    return MBOPhase::Simplified;
}

SIPGPhase nextSIPGPhase(SIPGPhase p) {
    switch (p) {
        case SIPGPhase::MeshOnly:   return SIPGPhase::CrossField;
        case SIPGPhase::CrossField: return SIPGPhase::Stepping;
        case SIPGPhase::Stepping:   return SIPGPhase::CutSeams;
        case SIPGPhase::CutSeams:   return SIPGPhase::UVMesh;
        case SIPGPhase::UVMesh:     return SIPGPhase::UVMesh;
    }
    return SIPGPhase::UVMesh;
}

UMBERPhase nextUMBERPhase(UMBERPhase p) {
    switch (p) {
        case UMBERPhase::MeshOnly:   return UMBERPhase::CrossField;
        case UMBERPhase::CrossField: return UMBERPhase::Stepping;
        case UMBERPhase::Stepping:   return UMBERPhase::Frames;
        case UMBERPhase::Frames:     return UMBERPhase::Polysquare;
        case UMBERPhase::Polysquare: return UMBERPhase::Blocks;
        case UMBERPhase::Blocks:     return UMBERPhase::Simplified;
        case UMBERPhase::Simplified: return UMBERPhase::Simplified;
    }
    return UMBERPhase::Simplified;
}

MedialAxisPhase nextMedialAxisPhase(MedialAxisPhase p) {
    switch (p) {
        case MedialAxisPhase::MeshOnly:     return MedialAxisPhase::DelaunayMesh;
        case MedialAxisPhase::DelaunayMesh: return MedialAxisPhase::MedialAxis;
        case MedialAxisPhase::MedialAxis:   return MedialAxisPhase::Map;
        case MedialAxisPhase::Map:          return MedialAxisPhase::Map;
    }
    return MedialAxisPhase::Map;
}

const char *phaseName(Phase p) {
    switch (p) {
        case Phase::MeshOnly:      return "1) mesh";
        case Phase::CrossField:    return "2) crossfield";
        case Phase::Singularities: return "3) singularities";
        case Phase::CutSeams:      return "4) cut seams";
        case Phase::UVMesh:        return "5) UV mesh (MIQ)";
    }
    return "?";
}

const char *mboPhaseName(MBOPhase p) {
    switch (p) {
        case MBOPhase::MeshOnly:     return "1) mesh";
        case MBOPhase::CrossField:   return "2) MBO crossfield";
        case MBOPhase::Stepping:     return "3) MBO stepping";
        case MBOPhase::Separatrices: return "4) separatrices";
        case MBOPhase::Trace:        return "5) trace";
        case MBOPhase::Layout:       return "6) quad layout";
        case MBOPhase::Simplified:   return "7) simplified partition";
    }
    return "?";
}

const char *sipgPhaseName(SIPGPhase p) {
    switch (p) {
        case SIPGPhase::MeshOnly:   return "1) mesh";
        case SIPGPhase::CrossField: return "2) SIPG crossfield";
        case SIPGPhase::Stepping:   return "3) SIPG stepping";
        case SIPGPhase::CutSeams:   return "4) cut seams (combed)";
        case SIPGPhase::UVMesh:     return "5) UV mesh (UVGParam)";
    }
    return "?";
}

const char *umberPhaseName(UMBERPhase p) {
    switch (p) {
        case UMBERPhase::MeshOnly:   return "1) mesh";
        case UMBERPhase::CrossField: return "2) SIPG crossfield";
        case UMBERPhase::Stepping:   return "3) SIPG stepping";
        case UMBERPhase::Frames:     return "4) UMBER frame field";
        case UMBERPhase::Polysquare: return "5) polysquare (Sec. 4.3)";
        case UMBERPhase::Blocks:     return "6) block structure (Sec. 5)";
        case UMBERPhase::Simplified: return "7) chord collapse";
    }
    return "?";
}

const char *medialAxisPhaseName(MedialAxisPhase p) {
    switch (p) {
        case MedialAxisPhase::MeshOnly:     return "1) mesh";
        case MedialAxisPhase::DelaunayMesh: return "2) Delaunay re-triangulation";
        case MedialAxisPhase::MedialAxis:   return "3) Medial axis";
        case MedialAxisPhase::Map:          return "4) Boundary map";
    }
    return "?";
}

const char *oasisPhaseName(OASISPhase p) {
    switch (p) {
        case OASISPhase::MeshOnly: return "1) mesh";
        case OASISPhase::Field:    return "2) quasi-eigenfunction";
    }
    return "?";
}

const char *modeName(Mode m) {
    switch (m) {
        case Mode::Unselected: return "unselected";
        case Mode::PolyVector: return "PolyVector";
        case Mode::MBO:        return "MBO";
        case Mode::MedialAxis: return "Medial Axis";
        case Mode::SIPG:       return "SIPG";
        case Mode::OASIS:      return "OASIS";
        case Mode::UMBER:      return "UMBER";
    }
    return "?";
}

} // anonymous namespace

// ── constructor ───────────────────────────────────────────────────────────────

CrossGenWidget::CrossGenWidget(const std::string &meshPath, QWidget *parent)
    : QOpenGLWidget(parent)
{
    mesh_ = std::make_shared<Mesh>(meshPath);

    // Compute view bounds
    bounds_  = viewer::computeBounds(*mesh_);
    avgEdge_ = viewer::averageTriangleEdgeLength(*mesh_);
    scale_   = 0.7 * avgEdge_;

    double dx  = bounds_.maxx - bounds_.minx;
    double dy  = bounds_.maxy - bounds_.miny;
    double ext = std::max(dx, dy);
    if (ext <= 0) ext = 1.0;
    pad_ = 0.1 * ext;

    view_.cx    = 0.5 * (bounds_.minx + bounds_.maxx);
    view_.cy    = 0.5 * (bounds_.miny + bounds_.maxy);
    view_.baseW = dx + 2.0 * pad_;
    view_.baseH = dy + 2.0 * pad_;
    if (view_.baseW <= 0.0) view_.baseW = 1.0;
    if (view_.baseH <= 0.0) view_.baseH = 1.0;
    view_.zoom  = 1.0;

    // Console setup
    console_.setMaxLines(8);
    {
        std::ostringstream oss;
        oss << "Loaded mesh: " << mesh_->triangles.size() << " triangles, "
            << mesh_->vertices.size() << " vertices";
        console_.log(oss.str());
    }

    // Timer drives continuous repaints (animation + general refresh)
    timer_ = new QTimer(this);
    connect(timer_, &QTimer::timeout, this, &CrossGenWidget::onTimer);
    timer_->start(16); // ~60 fps

    setFocusPolicy(Qt::StrongFocus);

    std::cerr << "[Viewer] Phase " << phaseName(phase_)
              << " (press '1' for PolyVector, '2' for MBO, '3' for Medial Axis, '4' for SIPG, '5' for OASIS, '6' for UMBER)\n";
}

// ── QOpenGLWidget overrides ───────────────────────────────────────────────────

void CrossGenWidget::initializeGL() {
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glShadeModel(GL_SMOOTH);
}

void CrossGenWidget::resizeGL(int w, int h) {
    viewer::resizeAndApplyOrtho(view_, fbw(), fbh());
    uvView_.fbw = view_.fbw;
    uvView_.fbh = view_.fbh;
}

void CrossGenWidget::paintGL() {
    // Apply current view/zoom state — must happen here where the GL context is current.
    viewer::applyOrtho(view_);

    // Run any lazy computations that are triggered by the current phase/mode.
    runComputations();

    // Clear
    glClearColor(0.1f, 0.1f, 0.12f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);

    // Choose render path
    bool isMBOStepping = (mode_ == Mode::MBO &&
                          mboPhase_ == MBOPhase::Stepping &&
                          mboSteppingStarted_ &&
                          mboStepCount_ < MBO_MAX_STEPS &&
                          !mboConverged_);

    bool isSIPGStepping = (sipgStageIsStepping() &&
                           sipgSteppingStarted_ &&
                           !sipgConverged_);

    bool isTracing = (mode_ == Mode::MBO &&
                      mboPhase_ == MBOPhase::Trace &&
                      separatrixTrace_ &&
                      !mboTracingFinished_);

    if (isMBOStepping) {
        renderMBOAnimation();
    } else if (isSIPGStepping) {
        renderSIPGAnimation();
    } else if (isTracing) {
        renderTraceAnimation();
    } else {
        renderNormal();
    }
}

// ── slot ──────────────────────────────────────────────────────────────────────

void CrossGenWidget::onTimer() {
    update(); // trigger paintGL
}

// ── key events ────────────────────────────────────────────────────────────────

void CrossGenWidget::keyPressEvent(QKeyEvent *event) {
    if (event->isAutoRepeat()) return;

    switch (event->key()) {
    case Qt::Key_Q:
    case Qt::Key_Escape:
        close();
        break;

    case Qt::Key_R:
        doReset();
        break;

    case Qt::Key_C:
        if (mode_ == Mode::OASIS) {
            // 'c' re-opens the dialog so lambda can be swept without a reset.
            if (promptOASISParameters())
                runOASIS();
        } else if (mode_ != Mode::Unselected) {
            advancePhase();
        }
        break;

    case Qt::Key_1:
        if (mode_ == Mode::Unselected && phase_ == Phase::MeshOnly) {
            mode_ = Mode::PolyVector;
            std::cerr << "[Viewer] Selected mode: " << modeName(mode_) << " (press 'c' to advance)\n";
            console_.log("Selected mode: PolyVector");
        }
        break;

    case Qt::Key_2:
        if (mode_ == Mode::Unselected && phase_ == Phase::MeshOnly) {
            mode_ = Mode::MBO;
            std::cerr << "[Viewer] Selected mode: " << modeName(mode_) << " (press 'c' to advance)\n";
            console_.log("Selected mode: MBO");
        }
        break;

    case Qt::Key_3:
        if (mode_ == Mode::Unselected && phase_ == Phase::MeshOnly) {
            mode_ = Mode::MedialAxis;
            std::cerr << "[Viewer] Selected mode: " << modeName(mode_) << " (press 'c' to advance)\n";
            console_.log("Selected mode: Medial Axis");
        }
        break;

    case Qt::Key_4:
        if (mode_ == Mode::Unselected && phase_ == Phase::MeshOnly) {
            mode_ = Mode::SIPG;
            std::cerr << "[Viewer] Selected mode: " << modeName(mode_) << " (press 'c' to advance)\n";
            console_.log("Selected mode: SIPG");
        }
        break;

    case Qt::Key_5:
        if (mode_ == Mode::Unselected && phase_ == Phase::MeshOnly) {
            // Unlike the other modes, OASIS needs a parameter before it can do
            // anything, so selecting it opens the dialog immediately. Backing
            // out of the dialog leaves the mode unselected.
            if (promptOASISParameters()) {
                mode_ = Mode::OASIS;
                std::cerr << "[Viewer] Selected mode: " << modeName(mode_)
                          << " (press 'c' to change lambda)\n";
                console_.log("Selected mode: OASIS");
                runOASIS();
            }
        }
        break;

    case Qt::Key_6:
        if (mode_ == Mode::Unselected && phase_ == Phase::MeshOnly) {
            mode_ = Mode::UMBER;
            std::cerr << "[Viewer] Selected mode: " << modeName(mode_) << " (press 'c' to advance)\n";
            console_.log("Selected mode: UMBER");
        }
        break;

    default:
        QOpenGLWidget::keyPressEvent(event);
        break;
    }
}

// ── mouse events ─────────────────────────────────────────────────────────────

void CrossGenWidget::mousePressEvent(QMouseEvent *event) {
    if (event->button() == Qt::RightButton) {
        rightDragging_ = true;
        lastMousePos_  = event->pos();
    }
}

void CrossGenWidget::mouseReleaseEvent(QMouseEvent *event) {
    if (event->button() == Qt::RightButton)
        rightDragging_ = false;
}

void CrossGenWidget::mouseMoveEvent(QMouseEvent *event) {
    if (!rightDragging_) return;

    QPoint delta    = event->pos() - lastMousePos_;
    lastMousePos_   = event->pos();

    // Scale from logical pixels to physical pixels for pan calculation
    double dpr = devicePixelRatio();
    bool panRight = inUVSplitScreen() && (lastMousePos_.x() * dpr > fbw() / 2);
    if (panRight)
        viewer::panView(uvView_, delta.x() * dpr, delta.y() * dpr);
    else
        viewer::panView(view_, delta.x() * dpr, delta.y() * dpr);
    update();
}

void CrossGenWidget::wheelEvent(QWheelEvent *event) {
    // On macOS trackpads Qt delivers high-resolution pixel deltas via pixelDelta().
    // angleDelta() is tiny (< 1 degree per gesture tick) on trackpads, so prefer
    // pixelDelta() when available and fall back to angleDelta() for click-wheel mice.
    double scrollSteps = 0.0;
    if (!event->pixelDelta().isNull()) {
        // pixelDelta().y() is in physical pixels; scale to a comfortable zoom rate.
        scrollSteps = event->pixelDelta().y() / 50.0;
    } else {
        double degrees = event->angleDelta().y() / 8.0;
        scrollSteps    = degrees / 15.0;
    }

    if (scrollSteps == 0.0) return;

    // Cursor in physical pixel coordinates
    double dpr = devicePixelRatio();
#if QT_VERSION >= QT_VERSION_CHECK(5, 14, 0)
    double cx = event->position().x() * dpr;
    double cy = event->position().y() * dpr;
#else
    double cx = event->pos().x() * dpr;
    double cy = event->pos().y() * dpr;
#endif
    // zoomView: positive scrollSteps => pow(0.9, positive) < 1 => zoom shrinks => zooms in ✓
    bool zoomRight = inUVSplitScreen() && (cx > fbw() / 2);
    if (zoomRight)
        viewer::zoomView(uvView_, scrollSteps, cx - fbw() / 2, cy);
    else
        viewer::zoomView(view_, scrollSteps, cx, cy);
    update();
}

// ── reset ─────────────────────────────────────────────────────────────────────

void CrossGenWidget::doReset() {
    field_.reset();
    cutMesh_.reset();
    miqSolver_.reset();
    crossField_.reset();
    sipgField_.reset();
    sipgCutMesh_.reset();
    sipgUVParam_.reset();
    separatrixTrace_.reset();
    quadLayout_.reset();
    simplified_.reset();
    delaunayMesh_.reset();
    medialAxisMap_.reset();
    medialAxis_.reset();
    oasis_.reset();
    oasisGuide_.reset();
    umber_.reset();
    umberCut_.reset();
    polysquare_.reset();
    chordCollapse_.reset();
    blockLayout_.reset();
    blocks_.reset();
    umberCorners_.clear();
    umberInternal_.clear();

    mode_     = Mode::Unselected;
    phase_    = Phase::MeshOnly;
    mboPhase_ = MBOPhase::MeshOnly;
    sipgPhase_ = SIPGPhase::MeshOnly;
    maPhase_  = MedialAxisPhase::MeshOnly;
    oasisPhase_ = OASISPhase::MeshOnly;
    umberPhase_ = UMBERPhase::MeshOnly;
    // oasisLambda_ deliberately survives a reset so it can be reused as the
    // dialog's default on the next run.

    singularitiesLogged_  = false;
    mboSteppingStarted_   = false;
    mboConverged_         = false;
    mboTracingStarted_    = false;
    mboTracingFinished_   = false;
    mboStepCount_         = 0;
    sipgSteppingStarted_  = false;
    sipgConverged_        = false;
    sipgStepCount_        = 0;
    umberAnnounced_       = false;
    umberAttempted_       = false;
    polysquareAnnounced_  = false;
    polysquareAttempted_  = false;
    blocksAttempted_      = false;

    view_.cx    = 0.5 * (bounds_.minx + bounds_.maxx);
    view_.cy    = 0.5 * (bounds_.miny + bounds_.maxy);
    view_.baseW = (bounds_.maxx - bounds_.minx) + 2.0 * pad_;
    view_.baseH = (bounds_.maxy - bounds_.miny) + 2.0 * pad_;
    if (view_.baseW <= 0.0) view_.baseW = 1.0;
    if (view_.baseH <= 0.0) view_.baseH = 1.0;
    view_.zoom  = 1.0;
    viewer::applyOrtho(view_);

    console_.clear();
    console_.setMaxLines(8);
    {
        std::ostringstream oss;
        oss << "Loaded mesh: " << mesh_->triangles.size() << " triangles, "
            << mesh_->vertices.size() << " vertices";
        console_.log(oss.str());
    }
    console_.log("[Reset] Restarted viewer.");
    std::cerr << "[Viewer] Reset. Phase " << phaseName(phase_)
              << " (press '1' for PolyVector, '2' for MBO, '3' for Medial Axis, '4' for SIPG, '5' for OASIS, '6' for UMBER)\n";
}

// ── OASIS parameter dialog ───────────────────────────────────────────────────

bool CrossGenWidget::promptOASISParameters() {
    // lambda is the only genuinely required input. The density field r of Eq. 8
    // is left uniform: a *constant* r is redundant with lambda, since the
    // modulated operator is grad^2/r and only the product lambda*r sets the
    // local wavelength. A spatially varying r is the real feature, and needs a
    // field editor rather than a spin box.
    //
    // A quasi-eigenfunction oscillates like cos(u*sqrt(-lambda)), so its
    // critical points — and hence the quad cells of the Morse-Smale complex —
    // sit a distance pi/sqrt(-lambda) apart. That is the number the user
    // actually cares about, so the dialog shows it live alongside the ratio to
    // the background triangle size, which the paper needs at 0.25 or below for
    // the critical points to be resolved at all (Sec. 3.3).
    const double defaultQuad = 8.0 * avgEdge_;
    const double defaultLambda =
        (oasisLambda_ < 0.0) ? oasisLambda_
                             : -(M_PI / defaultQuad) * (M_PI / defaultQuad);

    QDialog dlg(this);
    dlg.setWindowTitle("OASIS quasi-eigenfunction");

    auto *lambdaBox = new QDoubleSpinBox(&dlg);
    lambdaBox->setRange(-1e9, -1e-6);
    lambdaBox->setDecimals(2);
    lambdaBox->setValue(defaultLambda);
    lambdaBox->setSingleStep(10.0);
    lambdaBox->setToolTip("Helmholtz parameter of Eq. 1. Must be negative.\n"
                          "It does not have to be an eigenvalue.");

    auto *derived = new QLabel(&dlg);
    derived->setTextFormat(Qt::PlainText);

    auto updateDerived = [this, lambdaBox, derived]() {
        const double lam = lambdaBox->value();
        const double quad = M_PI / std::sqrt(-lam);
        const double ratio = avgEdge_ / quad;
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(4)
            << "quad edge  pi/sqrt(-lambda) = " << quad << "\n"
            << "mesh edge / quad edge       = " << std::setprecision(3) << ratio;
        if (ratio > 0.25)
            oss << "   <-- too coarse (want <= 0.25)";
        derived->setText(QString::fromStdString(oss.str()));
    };
    QObject::connect(lambdaBox, &QDoubleSpinBox::valueChanged, &dlg, updateDerived);
    updateDerived();

    // Vibration enhancement (Sec. 3.4). On by default. It is a second,
    // nonlinear pass, so it costs time and buys nothing on a field whose
    // amplitudes are already balanced — but where it matters (a rotationally
    // symmetric region, which vibrates radially and shows rings rather than
    // isolated critical points) the Morse-Smale complex is unusable without it.
    auto *vibrationBox = new QCheckBox("run vibration enhancement (Sec. 3.4)", &dlg);
    vibrationBox->setChecked(oasisVibrationIterations_ > 0);
    vibrationBox->setToolTip(
        "Equalizes the two local vibration amplitudes.\n"
        "Worth it when the field vibrates in one direction only —\n"
        "the paper's example is a disk-like, rotationally symmetric region.\n"
        "Trades exact boundary alignment for a penalty, so |Bf-C| grows.");

    auto *iterBox = new QSpinBox(&dlg);
    iterBox->setRange(1, 200);
    iterBox->setValue(oasisVibrationIterations_ > 0 ? oasisVibrationIterations_ : 10);
    iterBox->setEnabled(vibrationBox->isChecked());
    QObject::connect(vibrationBox, &QCheckBox::toggled, iterBox, &QSpinBox::setEnabled);

    // Orientation control (Sec. 5.1). On by default, which means pressing '5'
    // runs an MBO solve before the KKT one — a fraction of the KKT cost on the
    // meshes here. Untick it to get the plain 2014 system back, which is what
    // the boundary conditions alone give.
    auto *orientBox = new QCheckBox("guide orientation with an MBO cross field (Sec. 5.1)", &dlg);
    orientBox->setChecked(oasisOrientationWeight_ > 0.0);
    orientBox->setToolTip(
        "Adds gamma * E_Orient to the objective, so the principal directions\n"
        "of the field — and with them the arcs of the Morse-Smale complex —\n"
        "follow a cross field computed by MBO on this mesh.");

    auto *gammaBox = new QDoubleSpinBox(&dlg);
    gammaBox->setRange(0.01, 1e6);
    gammaBox->setDecimals(2);
    gammaBox->setSingleStep(1.0);
    gammaBox->setValue(oasisOrientationWeight_ > 0.0 ? oasisOrientationWeight_ : 10.0);
    gammaBox->setToolTip(
        "Weight of the orientation energy against the Helmholtz term,\n"
        "after both have been normalized by operator magnitude.\n"
        "10 follows the guiding field closely; below ~1 it barely turns,\n"
        "and far above it the cells stop being square.");

    // MBO stops as soon as the field stops moving, so this is a ceiling on the
    // wait rather than a target. The whole solve blocks the viewer.
    auto *mboIterBox = new QSpinBox(&dlg);
    mboIterBox->setRange(1, 5000);
    mboIterBox->setValue(oasisMBOIterations_);
    mboIterBox->setToolTip("Iteration cap for the MBO cross-field solve.\n"
                           "It exits early once the field converges.");

    auto *clearanceBox = new QDoubleSpinBox(&dlg);
    clearanceBox->setRange(0.0, 50.0);
    clearanceBox->setDecimals(1);
    clearanceBox->setSingleStep(0.5);
    clearanceBox->setValue(oasisGuideClearanceQuads_);
    clearanceBox->setToolTip(
        "Keeps the guiding field this many quad cells clear of the boundary,\n"
        "as the paper does. Within the band the Eq. 11 conditions already fix\n"
        "the orientation, and a guiding direction that disagrees with them is\n"
        "paid for in cell shape. 0 guides everywhere.");

    auto syncOrientation = [orientBox, gammaBox, mboIterBox, clearanceBox]() {
        const bool on = orientBox->isChecked();
        gammaBox->setEnabled(on);
        mboIterBox->setEnabled(on);
        clearanceBox->setEnabled(on);
    };
    QObject::connect(orientBox, &QCheckBox::toggled, &dlg, syncOrientation);
    syncOrientation();

    auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dlg);
    QObject::connect(buttons, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
    QObject::connect(buttons, &QDialogButtonBox::rejected, &dlg, &QDialog::reject);

    auto *form = new QFormLayout(&dlg);
    form->addRow("lambda", lambdaBox);
    form->addRow("implied sizing", derived);
    form->addRow(vibrationBox);
    form->addRow("iterations", iterBox);
    form->addRow(orientBox);
    form->addRow("gamma", gammaBox);
    form->addRow("MBO max iterations", mboIterBox);
    form->addRow("boundary clearance (quads)", clearanceBox);
    form->addRow(new QLabel(QString("mesh: %1 vertices, mean edge %2")
                                .arg(mesh_->vertices.size())
                                .arg(avgEdge_, 0, 'f', 4),
                            &dlg));
    form->addRow(buttons);

    if (dlg.exec() != QDialog::Accepted) return false;

    oasisLambda_ = lambdaBox->value();
    oasisVibrationIterations_ = vibrationBox->isChecked() ? iterBox->value() : 0;
    oasisOrientationWeight_ = orientBox->isChecked() ? gammaBox->value() : 0.0;
    oasisMBOIterations_ = mboIterBox->value();
    oasisGuideClearanceQuads_ = clearanceBox->value();
    return true;
}

// The guiding field is a full MBO run, blocking: the KKT solve that follows
// needs the finished field, and there is nothing useful to draw in between.
bool CrossGenWidget::buildOASISGuidingField() {
    auto t0 = Clock::now();
    try {
        oasisGuide_ = std::make_shared<CrossField>(mesh_, oasisMBOIterations_);
        oasisGuide_->initialize(1);
        oasisGuide_->runMBO();
        oasisGuide_->computeSingularities();
    } catch (const std::exception &e) {
        oasisGuide_.reset();
        console_.log(std::string("[OASIS] guiding field FAILED: ") + e.what());
        std::cerr << "[Viewer] MBO for the OASIS guiding field failed: " << e.what() << "\n";
        return false;
    }
    auto t1 = Clock::now();

    std::ostringstream oss;
    oss << "[OASIS] guiding field: MBO <=" << oasisMBOIterations_ << " iters, error "
        << std::scientific << std::setprecision(2) << oasisGuide_->error << ", "
        << oasisGuide_->singularTriangles.size() << " singularities, "
        << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
    console_.log(oss.str());
    return true;
}

void CrossGenWidget::runOASIS() {
    // Drop the previous solution before rebuilding the field it points into.
    oasis_.reset();
    if (oasisOrientationWeight_ <= 0.0) {
        oasisGuide_.reset();
    } else if (!buildOASISGuidingField()) {
        // Fall through without orientation rather than leaving the user with
        // nothing: the unguided QE is still the paper's result.
        console_.log("[OASIS] continuing without orientation control");
    }

    auto t0 = Clock::now();
    try {
        if (oasisGuide_) {
            oasis_.emplace(oasisGuide_, oasisLambda_);
            oasis_->setOrientationWeight(oasisOrientationWeight_);
            // The clearance is in quad cells, and one quad cell is
            // pi/sqrt(-lambda) across (the QE's half period).
            const double clearance =
                oasisGuideClearanceQuads_ * M_PI / std::sqrt(-oasisLambda_);
            oasis_->setOrientationMask(oasis_->boundaryClearanceMask(clearance));
        } else {
            oasis_.emplace(mesh_, oasisLambda_);
        }
        oasis_->assemble();
        oasis_->solve();
        if (oasisVibrationIterations_ > 0) {
            vibrationBefore_ = oasis_->vibrationEnergy();
            oasis_->enhanceVibration(oasisVibrationIterations_);
        } else {
            vibrationBefore_ = -1.0;
        }
    } catch (const std::exception &e) {
        oasis_.reset();
        oasisPhase_ = OASISPhase::MeshOnly;
        console_.log(std::string("[OASIS] FAILED: ") + e.what());
        std::cerr << "[Viewer] OASIS failed: " << e.what() << "\n";
        return;
    }
    auto t1 = Clock::now();

    const Eigen::VectorXd &f = oasis_->f;
    // The ramp is diverging, so the range must be symmetric about zero.
    oasisAbsMax_ = std::max(std::fabs(f.minCoeff()), std::fabs(f.maxCoeff()));
    if (!(oasisAbsMax_ > 0.0)) oasisAbsMax_ = 1.0;

    oasisPhase_ = OASISPhase::Field;

    {
        std::ostringstream oss;
        oss << "[OASIS] lambda=" << oasisLambda_ << ", "
            << oasis_->numConstraints() << " constraints, solved in "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }
    {
        const double bcErr = (oasis_->getConstraintMatrix() * f -
                              oasis_->getConstraintRhs()).cwiseAbs().maxCoeff();
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(2)
            << "[OASIS] |Bf-C|inf=" << bcErr
            << "  |Lf-lambda f|=" << oasis_->residual()
            << "  quad edge ~" << std::fixed << std::setprecision(4)
            << M_PI / std::sqrt(-oasisLambda_);
        console_.log(oss.str());
    }
    if (oasisGuide_) {
        // The misalignment against the guiding field, which is what the
        // orientation term is there to reduce. A guiding field that agrees with
        // the boundary conditions — an MBO one does — lands well under a
        // degree; several degrees means the two are pulling against each other.
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(2)
            << "[OASIS] orientation gamma=" << oasisOrientationWeight_
            << ", misalignment " << oasis_->orientationError()
            << " deg (0 = aligned, 45 = worst)";
        console_.log(oss.str());
    }
    if (vibrationBefore_ >= 0.0) {
        // Mean E_a is the quantity Sec. 3.4 minimizes, so report it either
        // side of the pass: that is the direct readout of whether it helped.
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(4)
            << "[OASIS] vibration " << oasisVibrationIterations_ << " iters, mean E_a "
            << vibrationBefore_ << " -> " << oasis_->vibrationEnergy() << " (0 = isotropic)";
        console_.log(oss.str());
    }
}

// ── UMBER (Wang et al. 2022, Secs. 4.1-4.2) ──────────────────────────────────

// One blocking call: the cuts, then Eq. (1) over the whole l1 continuation
// schedule. Splitting it across frames would mean stopping L-BFGS mid-stage,
// and an intermediate iterate of a quartic energy is not a field worth drawing.
void CrossGenWidget::runUMBER() {
    umberAttempted_ = true;
    if (!sipgField_.has_value()) return;

    // Eq. (1) starts from a *converged* cross field: the comb in initialize()
    // assumes neighbouring triangles already agree up to a k*90-degree turn,
    // which a half-solved MBO field does not. Advancing out of the stepping
    // phase early therefore finishes the solve here rather than optimizing a
    // field that is still moving.
    if (!sipgConverged_) {
        auto t0 = Clock::now();
        sipgField_->runMBO();
        sipgField_->computeSingularities();
        sipgConverged_ = true;
        auto t1 = Clock::now();
        std::ostringstream oss;
        oss << "[UMBER] finished the SIPG solve first, error " << std::scientific
            << std::setprecision(3) << sipgField_->error << ", "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }

    // --- Cuts, Sec. 4.1 -----------------------------------------------------
    auto tc0 = Clock::now();
    try {
        umberCut_.emplace(mesh_);
    } catch (const std::exception &e) {
        umberCut_.reset();
        console_.log(std::string("[UMBER] cutting FAILED: ") + e.what());
        std::cerr << "[Viewer] HarmonicCut failed: " << e.what() << "\n";
        return;
    }
    auto tc1 = Clock::now();
    {
        const auto &rep = umberCut_->getReport();
        std::ostringstream oss;
        oss << "[UMBER] cuts: " << rep.voids << " void(s), " << rep.cutsMade
            << " cut(s), " << umberCut_->getCutEdges().size() << " edges, "
            << (rep.isDisk ? "disk \033[32m[PASS]\033[0m" : "not a disk \033[31m[FAIL]\033[0m")
            << ", " << formatMs(std::chrono::duration<double, std::milli>(tc1 - tc0).count());
        console_.log(oss.str());
    }

    // --- Frame field, Eq. (1) ----------------------------------------------
    UMBER::EnergyTerms before;
    size_t internalBefore = 0;
    auto t0 = Clock::now();
    try {
        umber_.emplace(*sipgField_, *umberCut_);
        umber_->setMaxIterations(UMBER_LBFGS_ITERATIONS);
        umber_->initialize();
        before = umber_->energy();
        internalBefore = umber_->internalSingularities().size();
        umber_->optimize();
    } catch (const std::exception &e) {
        umber_.reset();
        console_.log(std::string("[UMBER] FAILED: ") + e.what());
        std::cerr << "[Viewer] UMBER failed: " << e.what() << "\n";
        return;
    }
    auto t1 = Clock::now();

    umberCorners_  = umber_->boundarySingularities();
    umberInternal_ = umber_->internalSingularities();

    const UMBER::EnergyTerms after = umber_->energy();
    {
        std::ostringstream oss;
        oss << "[UMBER] E_total " << std::scientific << std::setprecision(3)
            << before.total << " -> " << after.total
            << " (smooth " << before.smooth << " -> " << after.smooth
            << ", align " << after.align << ")";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[UMBER] L-BFGS " << umber_->iterations() << " iters (cap "
            << UMBER_LBFGS_ITERATIONS << "/stage), |grad| " << std::scientific
            << std::setprecision(2) << umber_->gradientNorm() << ", "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }
    {
        // The two numbers Sec. 4.2 is judged on: what left the interior, and
        // what it turned into on the boundary.
        int convex = 0, reflex = 0, other = 0;
        for (const auto &[vid, k] : umberCorners_) {
            if (k == 1) ++convex;
            else if (k == -1) ++reflex;
            else ++other;
        }
        std::ostringstream oss;
        oss << "[UMBER] internal singularities " << internalBefore << " -> "
            << umberInternal_.size() << "; boundary corners " << convex << " convex, "
            << reflex << " reflex";
        if (other > 0) oss << ", " << other << " higher order";
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }
}

// The frame field guided deformation of Sec. 4.3. Blocking like runUMBER, and
// on the same scale: the Poisson solve of Eq. (6) is one factorization and the
// Eq. (9) continuation is a few thousand L-BFGS iterations.
void CrossGenWidget::runPolysquare() {
    polysquareAttempted_ = true;
    if (!umber_.has_value() || !umberCut_.has_value()) return;

    auto t0 = Clock::now();
    try {
        polysquare_.emplace(*umber_, *umberCut_);
        polysquare_->solve();
    } catch (const std::exception &e) {
        polysquare_.reset();
    blocks_.reset();
        console_.log(std::string("[Polysquare] FAILED: ") + e.what());
        std::cerr << "[Viewer] Polysquare failed: " << e.what() << "\n";
        return;
    }
    auto t1 = Clock::now();

    const Polysquare::Report &r = polysquare_->getReport();
    {
        std::ostringstream oss;
        oss << "[Polysquare] " << r.iterations << " L-BFGS iterations, transitions";
        for (int k : polysquare_->getTransitions()) oss << " " << k * 90 << "deg";
        if (polysquare_->getTransitions().empty()) oss << " none";
        oss << ", " << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }
    {
        // The turn count is the readout that says whether this is a polysquare
        // at all: the boundary should turn once per corner of the frame field
        // and nowhere else. Alignment cannot tell -- a staircase is axis
        // aligned too.
        std::ostringstream oss;
        oss << "[Polysquare] boundary turns " << r.turns << " (field asked for "
            << r.expectedTurns << "), alignment " << std::fixed << std::setprecision(2)
            << r.meanAlignDeg << " deg mean / " << r.maxAlignDeg << " worst";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Polysquare] scaled Jacobian min " << std::fixed << std::setprecision(3)
            << r.minScaledJacobian << ", avg " << r.avgScaledJacobian
            << ", flipped triangles " << r.flips;
        console_.log(oss.str());
    }

    viewer::computePolysquareBounds(*polysquare_, uvView_.cx, uvView_.cy,
                                    uvView_.baseW, uvView_.baseH);
    uvView_.zoom = 1.0;
    uvView_.fbw = view_.fbw;
    uvView_.fbh = view_.fbh;
}

// The iso-line tracing of Sec. 5. Milliseconds, so it just runs.
void CrossGenWidget::runBlocks() {
    blocksAttempted_ = true;
    if (!polysquare_.has_value()) return;

    // Both of these point into blocks_, so they go before it is replaced.
    chordCollapse_.reset();
    blockLayout_.reset();

    auto t0 = Clock::now();
    try {
        blocks_.emplace(*polysquare_);
        blocks_->build();
    } catch (const std::exception &e) {
        blocks_.reset();
        console_.log(std::string("[Blocks] FAILED: ") + e.what());
        std::cerr << "[Viewer] MotorcycleGraph failed: " << e.what() << "\n";
        return;
    }
    auto t1 = Clock::now();

    const MotorcycleGraph::Report &r = blocks_->getReport();
    std::ostringstream oss;
    oss << "[Blocks] " << r.blocks << " block(s) from " << r.motorcycles << " iso-line(s), "
        << r.crossings << " crossing(s), " << r.nodes << " node(s), "
        << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
    console_.log(oss.str());
    console_.log("[Blocks] yellow = polysquare corner, green = line leaving the model, "
                 "cyan = crossing");
}

// ── chord collapse ───────────────────────────────────────────────────────────

void CrossGenWidget::buildBlockLayout() {
    if (blockLayout_.has_value() || !blocks_.has_value()) return;

    auto t0 = Clock::now();
    try {
        blockLayout_.emplace(*blocks_);
        blockLayout_->build();
    } catch (const std::exception &e) {
        blockLayout_.reset();
        console_.log(std::string("[BlockLayout] FAILED: ") + e.what());
        std::cerr << "[Viewer] BlockLayout failed: " << e.what() << "\n";
        return;
    }
    auto t1 = Clock::now();

    const QuadLayout::Report &r = blockLayout_->getLayout().getReport();
    std::ostringstream oss;
    oss << "[BlockLayout] " << r.faces << " block(s), " << r.arcs << " side(s), " << r.nodes
        << " node(s), " << r.quadFaces << " four-sided, "
        << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
    console_.log(oss.str());
    // A block that did not come out four-sided is a place the tracing left
    // open, and no chord can be walked through one, so it is worth saying
    // before the dialog reports fewer chords than the model looks like it has.
    if (r.badFaces > 0 || r.arcCrossings > 0) {
        std::ostringstream bad;
        bad << "[BlockLayout] " << r.badFaces << " block(s) not four-sided, " << r.arcCrossings
            << " side(s) crossing: chords through them are blocked";
        console_.log(bad.str());
    }
}

// The dialog reports what the settings would do before they are applied, which
// is the only way to choose the width: the number itself means nothing, and
// "how many chords does it let through, and how much thinner would the next one
// need me to be" is the question actually being asked. Enumerating the chords
// is a walk over the sides of the structure, so it is cheap enough to redo on
// every keystroke.
bool CrossGenWidget::promptChordCollapseParameters() {
    if (!blockLayout_.has_value()) return false;
    const QuadLayout &layout = blockLayout_->getLayout();

    QDialog dlg(this);
    dlg.setWindowTitle("Chord collapse");

    auto *widthBox = new QDoubleSpinBox(&dlg);
    widthBox->setRange(0.0, 10.0);
    widthBox->setDecimals(3);
    widthBox->setSingleStep(0.05);
    widthBox->setValue(chordSettings_.maxWidth);
    widthBox->setToolTip(
        "Rule 1. A chord is collapsed only if every one of its rungs is\n"
        "shorter than this, in units of the mean side length of the block\n"
        "structure. Well under 1: a chord as wide as a block is a partition\n"
        "of the model and not a sliver.");

    auto *aspectBox = new QDoubleSpinBox(&dlg);
    aspectBox->setRange(0.0, 100.0);
    aspectBox->setDecimals(1);
    aspectBox->setSingleStep(0.5);
    aspectBox->setValue(chordSettings_.minAspect);
    aspectBox->setToolTip(
        "Rule 1 from the other side: how many times longer than wide a chord\n"
        "has to be. A short fat chord and a long thin one can have the same\n"
        "rungs and only the second is a sliver. 0 turns it off.");

    auto *maxBox = new QSpinBox(&dlg);
    maxBox->setRange(0, 100000);
    maxBox->setValue(chordSettings_.maxCollapses);
    maxBox->setToolTip("Cap on how many chords the greedy loop takes.\n"
                       "0 leaves the structure alone, for comparison.");

    auto *preview = new QLabel(&dlg);
    preview->setTextFormat(Qt::PlainText);

    auto updatePreview = [&layout, widthBox, aspectBox, preview]() {
        ChordCollapse::Settings s;
        s.maxWidth = widthBox->value();
        s.minAspect = aspectBox->value();
        ChordCollapse probe(layout, s);
        probe.enumerateChords();
        const ChordCollapse::Report &r = probe.getReport();

        // The thinnest chord the width rule is currently turning away, which is
        // exactly how far the threshold would have to move to take one more.
        double nextWidth = -1.0;
        for (const auto &c : probe.getChords()) {
            if (c.block != ChordCollapse::Block::TooThick) continue;
            if (nextWidth < 0.0 || c.width < nextWidth) nextWidth = c.width;
        }

        std::ostringstream oss;
        oss << std::fixed << std::setprecision(3)
            << "mean side = " << probe.getWidthScale() << " in model units\n"
            << r.blocksBefore << " block(s), " << r.chordsSeen << " chord(s), "
            << r.collapsible << " collapsible here";
        for (int i = 1; i < static_cast<int>(ChordCollapse::Block::Count); ++i) {
            if (r.blockCount[i] == 0) continue;
            oss << "\n  " << r.blockCount[i] << " "
                << ChordCollapse::blockName(static_cast<ChordCollapse::Block>(i));
        }
        if (nextWidth >= 0.0)
            oss << "\nthe next one needs a width of " << nextWidth;
        preview->setText(QString::fromStdString(oss.str()));
    };
    QObject::connect(widthBox, &QDoubleSpinBox::valueChanged, &dlg, updatePreview);
    QObject::connect(aspectBox, &QDoubleSpinBox::valueChanged, &dlg, updatePreview);
    updatePreview();

    auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dlg);
    QObject::connect(buttons, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
    QObject::connect(buttons, &QDialogButtonBox::rejected, &dlg, &QDialog::reject);

    auto *form = new QFormLayout(&dlg);
    form->addRow("max rung width (mean sides)", widthBox);
    form->addRow("min length / width", aspectBox);
    form->addRow("max collapses", maxBox);
    form->addRow("at these settings", preview);
    form->addRow(new QLabel("A chord is refused whole if any rung joins two\n"
                            "boundaries or pinches one; a rung with one end on\n"
                            "the boundary always contracts onto that end.",
                            &dlg));
    form->addRow(buttons);

    if (dlg.exec() != QDialog::Accepted) return false;

    chordSettings_.maxWidth = widthBox->value();
    chordSettings_.minAspect = aspectBox->value();
    chordSettings_.maxCollapses = maxBox->value();
    return true;
}

void CrossGenWidget::runChordCollapse() {
    if (!blockLayout_.has_value()) return;

    auto t0 = Clock::now();
    // From the structure the tracing left, every time: the settings are a
    // heuristic being tuned, and collapsing on top of the last result would
    // mean the answer depended on which thresholds had been tried before it.
    chordCollapse_.emplace(blockLayout_->getLayout(), chordSettings_);
    chordCollapse_->run();
    auto t1 = Clock::now();

    const ChordCollapse::Report &r = chordCollapse_->getReport();
    std::ostringstream oss;
    oss << "[Collapse] " << r.blocksBefore << " -> " << r.blocksAfter << " block(s) over "
        << r.collapses << " chord collapse(s), widest " << std::fixed << std::setprecision(2)
        << r.widestCollapsed << " of a mean side, "
        << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
    console_.log(oss.str());

    std::ostringstream why;
    why << "[Collapse] " << r.chordsSeen << " chord(s) left, " << r.collapsible
        << " still collapsible";
    for (int i = 1; i < static_cast<int>(ChordCollapse::Block::Count); ++i) {
        if (r.blockCount[i] == 0) continue;
        why << ", " << r.blockCount[i] << " "
            << ChordCollapse::blockName(static_cast<ChordCollapse::Block>(i));
    }
    console_.log(why.str());
    if (r.rolledBack > 0) {
        std::ostringstream rb;
        rb << "[Collapse] " << r.rolledBack << " collapse(s) undone for leaving a broken layout ("
           << r.rbBlocks << " block count, " << r.rbBad << " not four-sided, " << r.rbCrossings
           << " crossing, " << r.rbArea << " area)";
        console_.log(rb.str());
    }
    console_.log("[Collapse] grey = the structure before, red = after; "
                 "press 'c' to try other settings");
}

// ── phase advancement ────────────────────────────────────────────────────────

void CrossGenWidget::advancePhase() {
    if (mode_ == Mode::PolyVector) {
        Phase old = phase_;
        phase_ = nextPhase(phase_);
        if (phase_ != old)
            std::cerr << "[Viewer] Phase " << phaseName(phase_) << "\n";
    } else if (mode_ == Mode::MBO) {
        MBOPhase old = mboPhase_;
        mboPhase_ = nextMBOPhase(mboPhase_);
        if (mboPhase_ != old)
            std::cerr << "[Viewer] MBO Phase " << mboPhaseName(mboPhase_) << "\n";
    } else if (mode_ == Mode::MedialAxis) {
        MedialAxisPhase old = maPhase_;
        maPhase_ = nextMedialAxisPhase(maPhase_);
        if (maPhase_ != old)
            std::cerr << "[Viewer] Medial Axis Phase " << medialAxisPhaseName(maPhase_) << "\n";
    } else if (mode_ == Mode::SIPG) {
        SIPGPhase old = sipgPhase_;
        sipgPhase_ = nextSIPGPhase(sipgPhase_);
        if (sipgPhase_ != old)
            std::cerr << "[Viewer] SIPG Phase " << sipgPhaseName(sipgPhase_) << "\n";
    } else if (mode_ == Mode::UMBER) {
        UMBERPhase old = umberPhase_;
        umberPhase_ = nextUMBERPhase(umberPhase_);
        if (umberPhase_ != old)
            std::cerr << "[Viewer] UMBER Phase " << umberPhaseName(umberPhase_) << "\n";

        // The last phase is a dialog rather than a picture to press on to, and
        // it stays that way: 'c' at it re-opens the dialog, so a threshold can
        // be tried, looked at, and tried again.
        if (umberPhase_ == UMBERPhase::Simplified) {
            // The blocks are normally built lazily by the frame after the one
            // that asked for them, so make sure they are there rather than
            // assuming a frame has been drawn since.
            if (!blocksAttempted_) runBlocks();
            buildBlockLayout();
            if (blockLayout_.has_value() && promptChordCollapseParameters())
                runChordCollapse();
        }
    }
}

// ── stages shared with SIPG mode ─────────────────────────────────────────────

bool CrossGenWidget::sipgStageWantsField() const {
    return (mode_ == Mode::SIPG  && sipgPhase_  >= SIPGPhase::CrossField) ||
           (mode_ == Mode::UMBER && umberPhase_ >= UMBERPhase::CrossField);
}

bool CrossGenWidget::sipgStageIsStepping() const {
    return (mode_ == Mode::SIPG  && sipgPhase_  == SIPGPhase::Stepping) ||
           (mode_ == Mode::UMBER && umberPhase_ == UMBERPhase::Stepping);
}

// ── lazy computations ────────────────────────────────────────────────────────

void CrossGenWidget::runComputations() {
    // ── MBO: Initialize CrossField ────────────────────────────────────────────
    if (mode_ == Mode::MBO && mboPhase_ >= MBOPhase::CrossField && !crossField_.has_value()) {
        auto t0 = Clock::now();
        crossField_.emplace(mesh_);
        crossField_->initialize(1);
        auto t1 = Clock::now();
        console_.log("[MBO] Initialized CrossField: " +
                     formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count()));
    }

    // ── MBO: Kick off stepping ────────────────────────────────────────────────
    if (mode_ == Mode::MBO && mboPhase_ == MBOPhase::Stepping &&
        crossField_.has_value() && !mboSteppingStarted_) {
        mboSteppingStarted_ = true;
        mboStepCount_       = 0;
        console_.log("[MBO] Starting " + std::to_string(MBO_MAX_STEPS) + " iterations...");
    }

    // ── MBO: Run 2 stepping iterations per frame ──────────────────────────────
    if (mode_ == Mode::MBO && mboPhase_ == MBOPhase::Stepping &&
        mboSteppingStarted_ && mboStepCount_ < MBO_MAX_STEPS && !mboConverged_) {
        double nv = static_cast<double>(crossField_->u_k.size());
        for (int i = 0; i < 2 && mboStepCount_ < MBO_MAX_STEPS; ++i) {
            crossField_->step();
            ++mboStepCount_;
            if (crossField_->error < 2.0 * nv * 1e-7) {
                console_.log("[MBO] Convergence at step " + std::to_string(mboStepCount_) +
                             " error=" + std::to_string(crossField_->error));
                mboConverged_ = true;
                break;
            }
        }
        crossField_->computeSingularities();

        std::ostringstream oss;
        oss << "[MBO] Step " << mboStepCount_ << "/" << MBO_MAX_STEPS;
        console_.log(oss.str());
    }

    // ── MBO: Build SeparatrixTrace ────────────────────────────────────────────
    if (mode_ == Mode::MBO && mboPhase_ >= MBOPhase::Separatrices &&
        crossField_.has_value() && !separatrixTrace_) {
        auto t0 = Clock::now();
        auto cfPtr = std::shared_ptr<CrossField>(&*crossField_, [](CrossField *) {});
        // Solve for where in its triangle the singularity actually sits rather
        // than taking the barycentre. The ports are launched from that point,
        // so a barycentre displaces every separatrix leaving it, and the
        // partition pays for it at the far end: over the sixteen models it
        // costs 40 components and 12 T-junctions, and on geom006 -- a box with
        // a round hole -- the difference is 12 components and none against 14
        // and four.
        separatrixTrace_ = std::make_shared<SeparatrixTrace>(cfPtr, true);
        auto t1 = Clock::now();
        std::ostringstream oss;
        oss << "[Separatrices] Initialized " << separatrixTrace_->separatrices.size()
            << " separatrices from " << separatrixTrace_->singularities.size()
            << " singularities: "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }

    // ── MBO: Step tracing one iteration per frame ─────────────────────────────
    if (mode_ == Mode::MBO && mboPhase_ == MBOPhase::Trace &&
        separatrixTrace_ && !mboTracingFinished_) {
        if (!mboTracingStarted_) {
            mboTracingStarted_ = true;
            console_.log("[Trace] Starting separatrix tracing...");
        }
        separatrixTrace_->stepAndCheck();
        if (separatrixTrace_->finishedTracing) {
            mboTracingFinished_ = true;
            console_.log("[Trace] Tracing complete.");
        }
    }

    // ── MBO: Build the quad layout the separatrices cut out ───────────────────
    if (mode_ == Mode::MBO && mboPhase_ >= MBOPhase::Layout && separatrixTrace_ &&
        !quadLayout_.has_value()) {
        // Skipping ahead past the trace animation is allowed, so finish the
        // tracing here rather than assuming a frame of it has run.
        if (!separatrixTrace_->finishedTracing) {
            separatrixTrace_->run();
            mboTracingFinished_ = true;
        }
        auto t0 = Clock::now();
        quadLayout_.emplace(*separatrixTrace_);
        quadLayout_->build();
        auto t1 = Clock::now();

        const QuadLayout::Report &r = quadLayout_->getReport();
        std::ostringstream oss;
        oss << "[Layout] " << r.faces << " component(s), " << r.quadFaces << " four-sided, "
            << r.nodes << " node(s), " << r.arcs << " arc(s), " << r.tJunctions
            << " T-junction(s), "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
        if (r.arcCrossings || r.danglingEnds) {
            std::ostringstream bad;
            bad << "[Layout] not a valid layout: " << r.arcCrossings
                << " arc crossing(s) with no node, " << r.danglingEnds << " loose end(s)";
            console_.log(bad.str());
        }
    }

    // ── MBO: Sec. 4, collapse the chords the layout does not need ─────────────
    if (mode_ == Mode::MBO && mboPhase_ >= MBOPhase::Simplified && quadLayout_.has_value() &&
        !simplified_.has_value()) {
        auto t0 = Clock::now();
        simplified_.emplace(*quadLayout_);
        simplified_->run();
        auto t1 = Clock::now();

        const PartitionSimplify::Report &r = simplified_->getReport();
        const QuadLayout::Report &sl = simplified_->getLayout().getReport();
        std::ostringstream oss;
        oss << "[Simplify] " << r.componentsBefore << " -> " << sl.faces << " component(s) ("
            << sl.quadFaces << " four-sided), " << r.tJunctionsBefore << " -> " << sl.tJunctions
            << " T-junction(s), over " << r.collapses << " chord collapse(s): "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
        std::ostringstream why;
        why << "[Simplify] " << r.blockedByConditions << " chord(s) blocked by Sec. 4's conditions, "
            << r.blockedByEnergy << " by the energy, " << r.rolledBack
            << " collapse(s) undone for breaking Proposition 2";
        console_.log(why.str());
    }

    // ── SIPG: Initialize ──────────────────────────────────────────────────────
    if (sipgStageWantsField() && !sipgField_.has_value()) {
        auto t0 = Clock::now();
        sipgField_.emplace(mesh_);
        sipgField_->initialize();
        auto t1 = Clock::now();
        console_.log("[SIPG] Initialized: " +
                     formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count()));
    }

    // ── SIPG: Kick off stepping ───────────────────────────────────────────────
    if (sipgStageIsStepping() && sipgField_.has_value() && !sipgSteppingStarted_) {
        sipgSteppingStarted_ = true;
        sipgStepCount_       = 0;
        console_.log("[SIPG] Starting MBO iterations (" +
                     std::to_string(sipgField_->getMesh().triangles.size()) + " tris)...");
    }

    // ── SIPG: Run 2 stepping iterations per frame ─────────────────────────────
    if (sipgStageIsStepping() && sipgSteppingStarted_ && !sipgConverged_) {
        double ntris = static_cast<double>(mesh_->triangles.size());
        for (int i = 0; i < 2 && sipgStepCount_ < 500; ++i) {
            sipgField_->step();
            ++sipgStepCount_;
            if (sipgField_->error < 2.0 * ntris * 1e-8) {
                console_.log("[SIPG] Converged at step " + std::to_string(sipgStepCount_) +
                             " error=" + std::to_string(sipgField_->error));
                sipgConverged_ = true;
                break;
            }
        }
        sipgField_->computeSingularities();

        std::ostringstream stepMsg;
        stepMsg << "[SIPG] Step " << sipgStepCount_ << "  error=" << std::scientific
                << std::setprecision(3) << sipgField_->error;
        console_.log(stepMsg.str());
    }

    // ── UMBER: cuts + Eq. (1) on top of the SIPG field ────────────────────────
    if (mode_ == Mode::UMBER && umberPhase_ >= UMBERPhase::Frames &&
        sipgField_.has_value() && !umberAttempted_) {
        if (!umberAnnounced_) {
            // runComputations() runs at the top of paintGL, so returning here
            // lets this frame draw the notice; the solve, which holds the GUI
            // thread for seconds, starts on the next one.
            console_.log("[UMBER] optimizing the frame field (Eq. 1), this blocks...");
            umberAnnounced_ = true;
        } else {
            runUMBER();
        }
    }

    // ── UMBER: the polysquare of Sec. 4.3 on top of the frame field ──────────
    if (mode_ == Mode::UMBER && umberPhase_ >= UMBERPhase::Polysquare &&
        umber_.has_value() && !polysquareAttempted_) {
        if (!polysquareAnnounced_) {
            console_.log("[Polysquare] deforming the cut mesh (Eq. 6, 9), this blocks...");
            polysquareAnnounced_ = true;
        } else {
            runPolysquare();
        }
    }

    // ── UMBER: the block structure of Sec. 5 ─────────────────────────────────
    if (mode_ == Mode::UMBER && umberPhase_ >= UMBERPhase::Blocks &&
        polysquare_.has_value() && !blocksAttempted_) {
        runBlocks();
    }

    // ── SIPG: Cut seams from converged SIPG field ─────────────────────────────
    if (mode_ == Mode::SIPG && sipgPhase_ >= SIPGPhase::CutSeams &&
        sipgField_.has_value() && !sipgCutMesh_.has_value()) {
        auto t0 = Clock::now();
        sipgCutMesh_.emplace(*sipgField_);
        auto t1 = Clock::now();
        std::ostringstream oss;
        oss << "[SIPG CutSeams] " << sipgCutMesh_->getCutEdges().size() << " cut edges: "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());

        const auto &rep = sipgCutMesh_->sanityCheck();
        console_.log(std::string("[SIPG CutSeams] ") +
                     (rep.looksLikeDisk ? "disk \033[32m[PASS]\033[0m"
                                        : "not a disk \033[31m[FAIL]\033[0m"));
    }

    // ── SIPG: UVGParam parametrization ───────────────────────────────────────
    if (mode_ == Mode::SIPG && sipgPhase_ >= SIPGPhase::UVMesh &&
        sipgCutMesh_.has_value() && !sipgUVParam_.has_value()) {
        auto t0 = Clock::now();
        try {
            sipgUVParam_.emplace(*sipgCutMesh_);
        } catch (const std::exception &e) {
            console_.log(std::string("[UVGParam] ERROR: ") + e.what());
        }
        if (sipgUVParam_.has_value()) {
            auto t1 = Clock::now();
            std::ostringstream oss;
            oss << "[UVGParam] Solved: "
                << sipgUVParam_->getU().size() << " vertices: "
                << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
            console_.log(oss.str());

            viewer::computeUVGParamBounds(*sipgUVParam_, uvView_.cx, uvView_.cy, uvView_.baseW, uvView_.baseH);
            uvView_.zoom = 1.0;
            uvView_.fbw  = view_.fbw;
            uvView_.fbh  = view_.fbh;
        }
    }

    // ── Medial Axis: Delaunay re-triangulation ────────────────────────────────
    if (mode_ == Mode::MedialAxis && maPhase_ >= MedialAxisPhase::DelaunayMesh && !delaunayMesh_) {
        auto t0 = Clock::now();

        // Build boundary loops
        std::unordered_map<int, std::vector<int>> bAdj;
        for (int beIdx : mesh_->boundaryEdges) {
            int a = mesh_->edges[beIdx][0];
            int b = mesh_->edges[beIdx][1];
            bAdj[a].push_back(b);
            bAdj[b].push_back(a);
        }

        std::unordered_set<int> visited;
        std::vector<std::vector<int>> loops;
        for (int bv : mesh_->boundaryVertices) {
            if (visited.count(bv)) continue;
            std::vector<int> loop;
            int prev = -1, curr = bv;
            while (true) {
                visited.insert(curr);
                loop.push_back(curr);
                int next = -1;
                for (int nb : bAdj[curr]) {
                    if (nb != prev && !visited.count(nb)) { next = nb; break; }
                }
                if (next == -1) break;
                prev = curr;
                curr = next;
            }
            if (loop.size() >= 3) loops.push_back(std::move(loop));
        }

        if (loops.empty()) {
            console_.log("[MedialAxis] No boundary loops found!");
        } else {
            std::unordered_map<int, int> vertRemap;
            std::vector<std::array<double, 2>> vertlist;

            for (const auto &loop : loops) {
                for (int vi : loop) {
                    if (!vertRemap.count(vi)) {
                        int newIdx = static_cast<int>(vertlist.size());
                        vertRemap[vi] = newIdx;
                        vertlist.push_back(mesh_->vertices[vi]);
                    }
                }
            }

            std::vector<std::vector<std::array<int, 2>>> segLoops;
            std::vector<int> loopTypes;
            for (size_t li = 0; li < loops.size(); ++li) {
                const auto &loop = loops[li];
                std::vector<std::array<int, 2>> segs;
                for (size_t i = 0; i < loop.size(); ++i) {
                    int a = vertRemap[loop[i]];
                    int b = vertRemap[loop[(i + 1) % loop.size()]];
                    segs.push_back({a, b});
                }
                segLoops.push_back(std::move(segs));
                loopTypes.push_back(li == 0 ? 0 : 1);
            }

            triangle_wrapper::TriangleMesher2D::Options opts;
            opts.just_delaunay = true;
            triangle_wrapper::TriangleMesher2D mesher(opts);

            triangle_wrapper::TriangleMesher2D::MeshInput input;
            input.vertlist      = vertlist;
            input.segment_loops = segLoops;
            input.type          = loopTypes;
            input.h             = 0.0;

            auto output = mesher.triangulate(input);
            delaunayMesh_ = std::make_shared<Mesh>(output.verts, output.triangles);

            auto t1 = Clock::now();
            std::ostringstream oss;
            oss << "[MedialAxis] Delaunay re-triangulation: "
                << delaunayMesh_->triangles.size() << " triangles, "
                << delaunayMesh_->vertices.size() << " vertices: "
                << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
            console_.log(oss.str());
        }
    }

    // ── Medial Axis: compute axis ─────────────────────────────────────────────
    if (mode_ == Mode::MedialAxis && maPhase_ >= MedialAxisPhase::MedialAxis &&
        delaunayMesh_ && !medialAxis_) {
        auto t0 = Clock::now();
        medialAxis_ = std::make_shared<MedialAxis>(delaunayMesh_);
        int rawV = static_cast<int>(medialAxis_->medialVertices.size());
        int rawE = static_cast<int>(medialAxis_->medialEdges.size());
        medialAxis_->deduplicateMedialVertices();
        auto t1 = Clock::now();
        std::ostringstream oss;
        oss << "[MedialAxis] Computed: " << rawV << " -> "
            << medialAxis_->medialVertices.size() << " vertices, "
            << rawE << " -> " << medialAxis_->medialEdges.size() << " edges: "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());

        const auto &st = medialAxis_->stats;
        if (st.degenerateCells || st.centersOutsideDomain || st.edgesCrossingBoundary ||
            st.nonManifoldBoundaryVertices) {
            std::ostringstream warn;
            warn << "[MedialAxis] Sampling warnings: "
                 << st.degenerateCells << " degenerate cells, "
                 << st.centersOutsideDomain << " centers outside, "
                 << st.edgesCrossingBoundary << " edges crossing the boundary, "
                 << st.nonManifoldBoundaryVertices << " pinched boundary vertices";
            console_.log(warn.str());
        }
    }

    // ── Medial Axis: the boundary -> axis map (Sec. 3) ────────────────────────
    if (mode_ == Mode::MedialAxis && maPhase_ >= MedialAxisPhase::Map &&
        medialAxis_ && !medialAxisMap_.has_value()) {
        auto t0 = Clock::now();
        medialAxisMap_.emplace(medialAxis_);
        auto t1 = Clock::now();

        const auto &st = medialAxisMap_->stats();
        std::ostringstream oss;
        oss << "[MedialAxis] Map phi: " << st.fans << " boundary fans ("
            << st.collapsedFans << " collapsed), "
            << medialAxisMap_->spokes().size() << " medial radii: "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
        if (st.skippedVertices || st.fallbackProjections) {
            std::ostringstream warn;
            warn << "[MedialAxis] Map warnings: " << st.skippedVertices
                 << " boundary vertices skipped, " << st.fallbackProjections
                 << " polar rays missed (closest-point fallback)";
            console_.log(warn.str());
        }
        console_.log("[MedialAxis] blue = bijective radii, magenta = polar "
                     "sections (axis endpoints)");
    }

    // ── PolyVector: solve crossfield ──────────────────────────────────────────
    if (mode_ == Mode::PolyVector && phase_ >= Phase::CrossField && !field_.has_value()) {
        auto t0 = Clock::now();
        field_.emplace(mesh_);
        field_->solveForPolyCoeffs();
        auto t1 = Clock::now();
        field_->convertToFieldVectors();
        auto t2 = Clock::now();
        console_.log("[CrossField] Solved poly-coeffs: " +
                     formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count()));
        console_.log("[CrossField] Converted to field vectors: " +
                     formatMs(std::chrono::duration<double, std::milli>(t2 - t1).count()));
    }

    // ── PolyVector: singularities ─────────────────────────────────────────────
    if (mode_ == Mode::PolyVector && phase_ >= Phase::Singularities &&
        field_.has_value() && !singularitiesLogged_) {
        auto t0 = Clock::now();
        field_->computeUSingularities();
        auto t1 = Clock::now();
        std::ostringstream oss;
        oss << "[Singularities] Found " << field_->uSingularities.size() << " singularities: "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
        singularitiesLogged_ = true;
    }

    // ── PolyVector: cut seams ─────────────────────────────────────────────────
    if (mode_ == Mode::PolyVector && phase_ >= Phase::CutSeams &&
        field_.has_value() && !cutMesh_.has_value()) {
        auto t0 = Clock::now();
        cutMesh_.emplace(*field_);
        auto t1 = Clock::now();
        std::ostringstream oss;
        oss << "[CutSeams] Generated " << cutMesh_->getCutEdges().size() << " cut edges: "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());

        std::cerr << "[Viewer] #tri=" << mesh_->triangles.size()
                  << " #vtx=" << mesh_->vertices.size()
                  << " | uSingularities=" << field_->uSingularities.size()
                  << " | cutEdges=" << cutMesh_->getCutEdges().size()
                  << " | singularityPathCutEdges="
                  << cutMesh_->getSingularityPathCutEdges().size() << "\n";
    }

    // ── PolyVector: MIQ parametrization ──────────────────────────────────────
    if (mode_ == Mode::PolyVector && phase_ >= Phase::UVMesh &&
        cutMesh_.has_value() && !miqSolver_.has_value()) {
        auto t0 = Clock::now();
        miqSolver_.emplace(*cutMesh_);
        miqSolver_->solve(100.0, 5.0, false, 10, 5000, true, true, true);
        auto t1 = Clock::now();

        int flips = miqSolver_->numFlips();
        const auto &UV = miqSolver_->getUV();
        std::ostringstream oss;
        oss << "[MIQ] Computed UV mesh: " << UV.rows() << " vertices, "
            << flips << " flips: "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
        std::cerr << "[Viewer] MIQ parametrization: " << UV.rows()
                  << " UV vertices, " << flips << " flipped triangles\n";

        viewer::computeUVMeshBounds(*miqSolver_, uvView_.cx, uvView_.cy, uvView_.baseW, uvView_.baseH);
        uvView_.zoom = 1.0;
        uvView_.fbw  = view_.fbw;
        uvView_.fbh  = view_.fbh;
    }
}

// ── animation render paths ───────────────────────────────────────────────────

void CrossGenWidget::renderMBOAnimation() {
    viewer::drawAxis(view_);
    viewer::drawMesh(*mesh_);
    viewer::drawVertexCrossFieldUK(*mesh_, *crossField_, scale_);

    double ballRadius = 0.5 * avgEdge_;
    for (const auto &sig : crossField_->singularTriangles) {
        int triIdx = sig.first;
        double crossIndex = sig.second;
        if (triIdx < 0 || triIdx >= static_cast<int>(mesh_->triangles.size())) continue;
        const Triangle &tri = mesh_->triangles[triIdx];
        const Point &p0 = mesh_->vertices[tri[0]];
        const Point &p1 = mesh_->vertices[tri[1]];
        const Point &p2 = mesh_->vertices[tri[2]];
        Point centroid = {(p0[0] + p1[0] + p2[0]) / 3.0,
                          (p0[1] + p1[1] + p2[1]) / 3.0};
        if (crossIndex > 0)
            viewer::drawDisk3D(centroid, ballRadius, 0.2f, 0.2f, 0.95f);
        else
            viewer::drawDisk3D(centroid, ballRadius, 0.95f, 0.2f, 0.2f);
    }

    renderOverlay("MBO stepping in progress...\npress 'q' to quit");
}

void CrossGenWidget::renderSIPGAnimation() {
    viewer::drawAxis(view_);
    viewer::drawMesh(*mesh_);
    viewer::drawTriangleCrossField(*mesh_, *sipgField_, scale_);

    double ballRadius = 0.5 * avgEdge_;
    for (const auto &[vertIdx, crossIndex] : sipgField_->singularVertices) {
        if (vertIdx < 0 || vertIdx >= static_cast<int>(mesh_->vertices.size())) continue;
        const Point &c = mesh_->vertices[vertIdx];
        if (crossIndex > 0)
            viewer::drawDisk3D(c, ballRadius, 0.2f, 0.2f, 0.95f);
        else
            viewer::drawDisk3D(c, ballRadius, 0.95f, 0.2f, 0.2f);
    }

    renderOverlay("SIPG stepping...\npress 'q' to quit");
}

void CrossGenWidget::renderTraceAnimation() {
    viewer::drawAxis(view_);
    viewer::drawMesh(*mesh_);

    glLineWidth(3.0f);
    for (const auto &sep : separatrixTrace_->separatrices) {
        if (sep.path.size() < 2) continue;
        if (sep.active)
            glColor3f(0.95f, 0.1f, 0.1f);
        else
            glColor3f(0.1f, 0.9f, 0.2f);
        glBegin(GL_LINE_STRIP);
        for (const auto &tp : sep.path)
            glVertex2d(tp.global_pos[0], tp.global_pos[1]);
        glEnd();
    }
    glLineWidth(1.0f);

    renderOverlay("Tracing separatrices...\npress 'q' to quit");
}

// ── split-screen helpers ─────────────────────────────────────────────────────

// Whether the right half of the window is showing a parameter domain, which is
// what decides where a drag or a scroll lands.
bool CrossGenWidget::inUVSplitScreen() const {
    return (mode_ == Mode::PolyVector && phase_ == Phase::UVMesh && miqSolver_.has_value()) ||
           (mode_ == Mode::SIPG && sipgPhase_ == SIPGPhase::UVMesh && sipgUVParam_.has_value()) ||
           // The chord collapse phase takes the whole window back: what it has
           // to show is the structure before against the structure after, and
           // both of those live in the model.
           (mode_ == Mode::UMBER && umberPhase_ >= UMBERPhase::Polysquare &&
            umberPhase_ <= UMBERPhase::Blocks && polysquare_.has_value());
}

void CrossGenWidget::applyHalfOrtho(int x, int vpW, const viewer::ViewState &vs) const {
    const int h = fbh();
    viewer::ViewState tmp = vs;
    tmp.fbw = vpW;
    tmp.fbh = h;
    double worldW = 1.0, worldH = 1.0;
    viewer::computeWorldBox(tmp, worldW, worldH);
    glViewport(x, 0, vpW, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(tmp.cx - 0.5 * worldW, tmp.cx + 0.5 * worldW,
            tmp.cy - 0.5 * worldH, tmp.cy + 0.5 * worldH, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void CrossGenWidget::drawSplitDivider(int halfW) const {
    const int w = fbw(), h = fbh();
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, static_cast<double>(w), 0.0, static_cast<double>(h), -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glLineWidth(2.0f);
    glColor3f(0.55f, 0.55f, 0.55f);
    glBegin(GL_LINES);
    glVertex2f(static_cast<float>(halfW), 0.0f);
    glVertex2f(static_cast<float>(halfW), static_cast<float>(h));
    glEnd();
    glLineWidth(1.0f);
    glViewport(0, 0, w, h);
}

// The mesh under the optimized frame: both directions of the frame, the cuts
// whose transitions the polysquare will use, and the corners the field put on
// the boundary. Shown on its own in the Frames phase and as the left half of
// the split screen once the parameterization exists, so that a corner can be
// found on the model and in the parameter domain at the same time.
void CrossGenWidget::renderUMBERField() {
    viewer::drawMesh(*mesh_);
    if (!umber_.has_value()) return;

    viewer::drawUField(*mesh_, umber_->getUField(), scale_);
    viewer::drawVField(*mesh_, umber_->getVField(), scale_);

    if (umberCut_.has_value() && !umberCut_->getCutEdges().empty())
        viewer::drawEdgeSetOnMesh(*mesh_, umberCut_->getCutEdges(), 1.0f, 0.2f, 0.9f, 3.5f);
    viewer::drawBoundaryEdges(*mesh_);

    const double ballRadius = 0.5 * avgEdge_;
    for (const auto &[vertIdx, k] : umberCorners_) {
        if (vertIdx < 0 || vertIdx >= static_cast<int>(mesh_->vertices.size())) continue;
        const Point &c = mesh_->vertices[vertIdx];
        if (k == 1)       viewer::drawDisk3D(c, ballRadius, 0.2f, 0.2f, 0.95f);
        else if (k == -1) viewer::drawDisk3D(c, ballRadius, 0.95f, 0.2f, 0.2f);
        else              viewer::drawDisk3D(c, ballRadius, 0.95f, 0.85f, 0.1f);
    }
    for (const auto &[vertIdx, winding] : umberInternal_) {
        if (vertIdx < 0 || vertIdx >= static_cast<int>(mesh_->vertices.size())) continue;
        viewer::drawDisk3D(mesh_->vertices[vertIdx], ballRadius, 0.1f, 0.9f, 0.2f);
    }
}

// ── normal render ─────────────────────────────────────────────────────────────

void CrossGenWidget::renderNormal() {
    if (mode_ == Mode::PolyVector && phase_ == Phase::UVMesh && miqSolver_.has_value()) {
        // ── Split-screen: left = original mesh, right = UV mesh ───────────────
        int w = fbw(), h = fbh();
        int halfW = w / 2;

        // ── Left panel: mesh with cut seams & singularities ───────────────────
        applyHalfOrtho(0, halfW, view_);
        {
            viewer::ViewState leftVs = view_;
            leftVs.fbw = halfW;
            leftVs.fbh = h;
            viewer::drawAxis(leftVs);
        }
        viewer::drawMesh(*mesh_);
        if (field_.has_value() && cutMesh_.has_value()) {
            viewer::drawUField(*mesh_, cutMesh_->getUField(), scale_);
            viewer::drawVField(*mesh_, cutMesh_->getVField(), scale_);
            if (!cutMesh_->getSingularityPathCutEdges().empty())
                viewer::drawEdgeSetOnMesh(*mesh_, cutMesh_->getCutEdges(), 1.0f, 0.75f, 0.1f, 4.0f);
            else
                viewer::drawEdgeSetOnMesh(*mesh_, cutMesh_->getCutEdges(), 1.0f, 0.2f, 0.9f, 3.5f);
            double ballRadius = 0.5 * avgEdge_;
            for (const auto &sig : field_->uSingularities) {
                int vid    = sig.first;
                int index4 = sig.second;
                if (vid < 0 || vid >= static_cast<int>(mesh_->vertices.size())) continue;
                const Point &c = mesh_->vertices[vid];
                if (index4 == 1)
                    viewer::drawDisk3D(c, ballRadius, 0.2f, 0.2f, 0.95f);
                else if (index4 == -1)
                    viewer::drawDisk3D(c, ballRadius, 0.95f, 0.2f, 0.2f);
            }
        }

        // ── Right panel: UV mesh ──────────────────────────────────────────────
        applyHalfOrtho(halfW, w - halfW, uvView_);
        viewer::drawUVMesh(*miqSolver_);
        if (field_.has_value() && cutMesh_.has_value()) {
            double uvRadius = 0.8;
            viewer::drawSingularitiesOnUV(*miqSolver_, *cutMesh_, *field_, uvRadius);
        }

        drawSplitDivider(halfW);
    } else if (mode_ == Mode::SIPG) {
        if (sipgPhase_ == SIPGPhase::UVMesh && sipgUVParam_.has_value()) {
            // ── Split-screen: left = mesh with cut seams, right = UVGParam ────
            int w = fbw(), h = fbh();
            int halfW = w / 2;

            // Left panel: mesh with combed field + cut seams
            applyHalfOrtho(0, halfW, view_);
            {
                viewer::ViewState leftVs = view_;
                leftVs.fbw = halfW;
                leftVs.fbh = h;
                viewer::drawAxis(leftVs);
            }
            viewer::drawMesh(*mesh_);
            if (sipgCutMesh_.has_value()) {
                viewer::drawUField(*mesh_, sipgCutMesh_->getUField(), scale_);
                viewer::drawVField(*mesh_, sipgCutMesh_->getVField(), scale_);
                if (!sipgCutMesh_->getSingularityPathCutEdges().empty())
                    viewer::drawEdgeSetOnMesh(*mesh_, sipgCutMesh_->getCutEdges(), 1.0f, 0.75f, 0.1f, 4.0f);
                else
                    viewer::drawEdgeSetOnMesh(*mesh_, sipgCutMesh_->getCutEdges(), 1.0f, 0.2f, 0.9f, 3.5f);
            }
            if (sipgField_.has_value()) {
                double ballRadius = 0.5 * avgEdge_;
                for (const auto &[vertIdx, crossIndex] : sipgField_->singularVertices) {
                    if (vertIdx < 0 || vertIdx >= static_cast<int>(mesh_->vertices.size())) continue;
                    const Point &c = mesh_->vertices[vertIdx];
                    if (crossIndex > 0)
                        viewer::drawDisk3D(c, ballRadius, 0.2f, 0.2f, 0.95f);
                    else
                        viewer::drawDisk3D(c, ballRadius, 0.95f, 0.2f, 0.2f);
                }
            }

            // Right panel: UVGParam
            applyHalfOrtho(halfW, w - halfW, uvView_);
            viewer::drawUVGParam(*sipgUVParam_);
            viewer::drawFlippedUVTriangles(*sipgUVParam_);
            if (sipgField_.has_value()) {
                double uvRadius = 0.5 * avgEdge_;
                viewer::drawSingularitiesOnUVG(*sipgUVParam_, sipgField_->singularVertices, uvRadius);
            }

            drawSplitDivider(halfW);
        } else {
        viewer::drawAxis(view_);
        viewer::drawMesh(*mesh_);
        if (sipgPhase_ >= SIPGPhase::CrossField && sipgField_.has_value()) {
            if (sipgPhase_ < SIPGPhase::CutSeams || !sipgCutMesh_.has_value()) {
                viewer::drawTriangleCrossField(*mesh_, *sipgField_, scale_);
            }
            double ballRadius = 0.5 * avgEdge_;
            for (const auto &[vertIdx, crossIndex] : sipgField_->singularVertices) {
                if (vertIdx < 0 || vertIdx >= static_cast<int>(mesh_->vertices.size())) continue;
                const Point &c = mesh_->vertices[vertIdx];
                if (crossIndex > 0)
                    viewer::drawDisk3D(c, ballRadius, 0.2f, 0.2f, 0.95f);
                else
                    viewer::drawDisk3D(c, ballRadius, 0.95f, 0.2f, 0.2f);
            }
        }
        if (sipgPhase_ >= SIPGPhase::CutSeams && sipgCutMesh_.has_value()) {
            // Draw combed u and v fields
            viewer::drawUField(*mesh_, sipgCutMesh_->getUField(), scale_);
            viewer::drawVField(*mesh_, sipgCutMesh_->getVField(), scale_);
            // Draw cut edges (singularity paths highlighted differently)
            if (!sipgCutMesh_->getSingularityPathCutEdges().empty())
                viewer::drawEdgeSetOnMesh(*mesh_, sipgCutMesh_->getCutEdges(), 1.0f, 0.75f, 0.1f, 4.0f);
            else
                viewer::drawEdgeSetOnMesh(*mesh_, sipgCutMesh_->getCutEdges(), 1.0f, 0.2f, 0.9f, 3.5f);
            // Draw natural boundary edges in red
            viewer::drawEdgeSetOnMesh(*mesh_, sipgCutMesh_->getNaturalBoundaryEdges(), 1.0f, 0.1f, 0.1f, 3.5f);
        }
        } // end else (non-UVMesh SIPG phases)
    } else if (mode_ == Mode::MBO) {
        viewer::drawAxis(view_);
        viewer::drawMesh(*mesh_);
        if (mboPhase_ >= MBOPhase::CrossField && crossField_.has_value()) {
            if (mboPhase_ >= MBOPhase::Stepping && mboStepCount_ > 0)
                viewer::drawVertexCrossFieldUK(*mesh_, *crossField_, scale_);
            else
                viewer::drawVertexCrossField(*mesh_, *crossField_, scale_);

            if (mboPhase_ < MBOPhase::Separatrices) {
                double ballRadius = 0.5 * avgEdge_;
                for (const auto &sig : crossField_->singularTriangles) {
                    int triIdx = sig.first;
                    double crossIndex = sig.second;
                    if (triIdx < 0 || triIdx >= static_cast<int>(mesh_->triangles.size())) continue;
                    const Triangle &tri = mesh_->triangles[triIdx];
                    const Point &p0 = mesh_->vertices[tri[0]];
                    const Point &p1 = mesh_->vertices[tri[1]];
                    const Point &p2 = mesh_->vertices[tri[2]];
                    Point centroid = {(p0[0] + p1[0] + p2[0]) / 3.0,
                                      (p0[1] + p1[1] + p2[1]) / 3.0};
                    if (crossIndex > 0)
                        viewer::drawDisk3D(centroid, ballRadius, 0.2f, 0.2f, 0.95f);
                    else
                        viewer::drawDisk3D(centroid, ballRadius, 0.95f, 0.2f, 0.2f);
                }
            }
        }

        // Once the layout exists it replaces the raw separatrices: its arcs are
        // the same curves cut at the nodes, plus the pieces of the boundary
        // that close the components, so drawing both would only double the
        // interior lines and still leave the outline out.
        if (mboPhase_ >= MBOPhase::Simplified && simplified_.has_value()) {
            viewer::drawQuadLayoutArcs(simplified_->getLayout(), 3.0f, 0.95f, 0.2f, 0.2f);
            // view_.zoom is 1.0 at fit and shrinks as the view zooms in, so
            // scaling the radius by it makes the markers shrink along with it
            // rather than staying a fixed size in mesh space and so covering
            // more and more of the screen as you zoom in on a component.
            viewer::drawQuadLayoutNodes(simplified_->getLayout(), 0.12 * avgEdge_ * view_.zoom);
        } else if (mboPhase_ >= MBOPhase::Layout && quadLayout_.has_value()) {
            viewer::drawQuadLayoutArcs(*quadLayout_, 3.0f, 0.95f, 0.2f, 0.2f);
        } else if (mboPhase_ >= MBOPhase::Separatrices && separatrixTrace_) {
            glLineWidth(3.0f);
            for (const auto &sep : separatrixTrace_->separatrices) {
                if (sep.path.size() < 2) continue;
                if (sep.active)
                    glColor3f(0.95f, 0.1f, 0.1f);
                else
                    glColor3f(0.1f, 0.9f, 0.2f);
                glBegin(GL_LINE_STRIP);
                for (const auto &tp : sep.path)
                    glVertex2d(tp.global_pos[0], tp.global_pos[1]);
                glEnd();
            }
            glLineWidth(1.0f);
        }
    } else if (mode_ == Mode::OASIS) {
        viewer::drawAxis(view_);
        if (oasisPhase_ == OASISPhase::Field && oasis_.has_value()) {
            // Field first, then the wireframe over it. The wireframe is kept
            // translucent and light: at these mesh densities an opaque one
            // buries the structure it is drawn over, and a light tint stays
            // visible against both the saturated lobes and the near-black
            // neutral around the zero level. The boundary goes on top last,
            // since the alignment claim is about it: the extrema should run
            // square to the boundary.
            viewer::drawScalarField(*mesh_, oasis_->f, oasisAbsMax_);
            viewer::drawMeshOverlay(*mesh_, 0.85f, 0.85f, 0.85f, 0.22f, 1.0f);
            // The guiding crosses go over the field, since the whole claim of
            // the orientation term is that the two line up: the lobes should
            // run along the crosses wherever the guiding field is active.
            if (oasisGuide_)
                viewer::drawVertexCrossFieldUK(*mesh_, *oasisGuide_, scale_);
            viewer::drawBoundaryEdges(*mesh_);
        } else {
            viewer::drawMesh(*mesh_);
        }
    } else if (mode_ == Mode::UMBER && umberPhase_ >= UMBERPhase::Simplified &&
               blockLayout_.has_value()) {
        // ── The structure the tracing left, against what the collapse made ──
        //
        // Both over the one mesh at the one scale, so that a chord that went is
        // a grey line with no red on it and everything else is red over grey.
        // The parameter domain is dropped here: the operation happens in the
        // model, and the picture that answers "which blocks did that remove"
        // is this one.
        viewer::drawAxis(view_);
        viewer::drawMesh(*mesh_);
        if (umberCut_.has_value() && !umberCut_->getCutEdges().empty())
            viewer::drawEdgeSetOnMesh(*mesh_, umberCut_->getCutEdges(), 1.0f, 0.2f, 0.9f, 2.0f);

        viewer::drawQuadLayoutArcs(blockLayout_->getLayout(), 2.0f, 0.45f, 0.45f, 0.5f);
        const QuadLayout &shown = chordCollapse_.has_value() ? chordCollapse_->getLayout()
                                                             : blockLayout_->getLayout();
        viewer::drawQuadLayoutArcs(shown, 4.0f, 0.95f, 0.25f, 0.2f);
        // view_.zoom is 1.0 at fit and shrinks as the view zooms in, so scaling
        // the radius by it keeps the markers the same size on screen instead of
        // swallowing a block once you zoom in on one.
        viewer::drawQuadLayoutNodes(shown, 0.12 * avgEdge_ * view_.zoom);
    } else if (mode_ == Mode::UMBER && umberPhase_ >= UMBERPhase::Polysquare &&
               polysquare_.has_value()) {
        // ── Split-screen: left = mesh and frame, right = the polysquare ─────
        const int w = fbw();
        const int halfW = w / 2;

        // The block structure replaces the frame on the left once it exists:
        // the iso-lines are drawn over the mesh they curve through, and the
        // arrows underneath would only crowd them.
        const bool showBlocks = (umberPhase_ >= UMBERPhase::Blocks && blocks_.has_value());

        applyHalfOrtho(0, halfW, view_);
        {
            viewer::ViewState leftVs = view_;
            leftVs.fbw = halfW;
            leftVs.fbh = fbh();
            viewer::drawAxis(leftVs);
        }
        if (showBlocks) {
            viewer::drawMesh(*mesh_);
            if (umberCut_.has_value()) {
                if (!umberCut_->getCutEdges().empty())
                    viewer::drawEdgeSetOnMesh(*mesh_, umberCut_->getCutEdges(), 1.0f, 0.2f, 0.9f, 2.0f);
                viewer::drawBlockBoundary(*polysquare_, *umberCut_, false, 4.0f);
            }
            viewer::drawBlockEdges(*blocks_, false, 4.0f);
            viewer::drawBlockNodes(*blocks_, false, 0.3 * avgEdge_);
        } else {
            renderUMBERField();
        }

        applyHalfOrtho(halfW, w - halfW, uvView_);
        viewer::drawFlippedPolysquareTriangles(*polysquare_);
        viewer::drawPolysquare(*polysquare_);
        if (umberCut_.has_value()) {
            viewer::drawPolysquareStructure(*polysquare_, *umberCut_);
            // A corner should look the same size on both halves, and the two
            // halves are fit to different world boxes, so the model-space
            // radius is carried over by the ratio between them. Taking it from
            // the live view states rather than the mesh extents keeps the two
            // matched while either panel is zoomed.
            viewer::ViewState left = view_;
            left.fbw = halfW;
            left.fbh = fbh();
            viewer::ViewState right = uvView_;
            right.fbw = w - halfW;
            right.fbh = fbh();
            double leftW = 1.0, leftH = 1.0, rightW = 1.0, rightH = 1.0;
            viewer::computeWorldBox(left, leftW, leftH);
            viewer::computeWorldBox(right, rightW, rightH);

            const double scaleToUV = (leftW > 0.0) ? rightW / leftW : 1.0;
            const double uvRadius = 0.5 * avgEdge_ * scaleToUV;
            if (showBlocks) {
                viewer::drawBlockBoundary(*polysquare_, *umberCut_, true, 4.0f);
                viewer::drawBlockEdges(*blocks_, true, 4.0f);
                viewer::drawBlockNodes(*blocks_, true, 0.3 * avgEdge_ * scaleToUV);
            } else {
                viewer::drawPolysquareCorners(*polysquare_, *umberCut_, umberCorners_, uvRadius);
            }
        }

        drawSplitDivider(halfW);
    } else if (mode_ == Mode::UMBER) {
        viewer::drawAxis(view_);
        if (umberPhase_ < UMBERPhase::Frames || !umber_.has_value()) {
            viewer::drawMesh(*mesh_);
            // The SIPG stages, drawn as SIPG mode draws them: the input field
            // and the singularities Sec. 4.2 is about to move.
            if (umberPhase_ >= UMBERPhase::CrossField && sipgField_.has_value()) {
                viewer::drawTriangleCrossField(*mesh_, *sipgField_, scale_);
                double ballRadius = 0.5 * avgEdge_;
                for (const auto &[vertIdx, crossIndex] : sipgField_->singularVertices) {
                    if (vertIdx < 0 || vertIdx >= static_cast<int>(mesh_->vertices.size())) continue;
                    const Point &c = mesh_->vertices[vertIdx];
                    if (crossIndex > 0)
                        viewer::drawDisk3D(c, ballRadius, 0.2f, 0.2f, 0.95f);
                    else
                        viewer::drawDisk3D(c, ballRadius, 0.95f, 0.2f, 0.2f);
                }
            }
        } else {
            // The optimized frame, its cuts and its corners. Both directions
            // are drawn because the field is non-symmetric here: u and v are
            // distinguishable, unlike the four indistinguishable arms of the
            // cross field it came from. Corners are blue convex (+1), red
            // reflex (-1), yellow higher order; a green disk is a defect that
            // failed to leave the interior.
            renderUMBERField();
        }
    } else if (mode_ == Mode::MedialAxis) {
        viewer::drawAxis(view_);
        // The map phase draws over the boundary alone: the triangulation
        // underneath would bury the radii.
        if (maPhase_ != MedialAxisPhase::Map) {
            if (delaunayMesh_)
                viewer::drawMesh(*delaunayMesh_);
            else
                viewer::drawMesh(*mesh_);
        }

        if (medialAxis_) {
            if (maPhase_ == MedialAxisPhase::Map && medialAxisMap_.has_value()) {
                viewer::drawBoundaryEdges(*medialAxis_->mesh);

                // The map phi drawn as its medial radii (Fig. 7/8 of the
                // paper): each boundary point joined to its image on the
                // axis. Radii first, the axis on top of them. Blue radii are
                // the bijective sections; magenta ones belong to a polar
                // section, where a run of boundary collapses onto a single
                // medial vertex (a discrete axis endpoint).
                glLineWidth(1.0f);
                glBegin(GL_LINES);
                for (const auto &s : medialAxisMap_->spokes()) {
                    if (s.collapsed)
                        glColor4f(0.9f, 0.25f, 0.9f, 0.85f);
                    else
                        glColor4f(0.35f, 0.55f, 0.95f, 0.55f);
                    glVertex2d(s.boundary[0], s.boundary[1]);
                    glVertex2d(s.medial[0], s.medial[1]);
                }
                glEnd();

                viewer::drawMedialAxis(*medialAxis_, avgEdge_ / 6.0);
            } else if (maPhase_ >= MedialAxisPhase::MedialAxis) {
                double ballRadius_ma = avgEdge_ / 5.0;
                viewer::drawMedialAxis(*medialAxis_, ballRadius_ma);
                for (int i = 0; i < static_cast<int>(medialAxis_->interiorAngle.size()); ++i) {
                    if (medialAxis_->isConcaveCorner(i)) {
                        const Point &p = medialAxis_->mesh->vertices[i];
                        viewer::drawDisk3D(p, ballRadius_ma, 0.95f, 0.5f, 0.1f);
                    } else if (medialAxis_->isSharpCorner(i)) {
                        const Point &p = medialAxis_->mesh->vertices[i];
                        viewer::drawDisk3D(p, ballRadius_ma, 0.95f, 0.2f, 0.2f);
                    }
                }
            }
        }
    } else {
        // PolyVector phases 1-4 (not UV)
        viewer::drawAxis(view_);
        if (phase_ >= Phase::MeshOnly)
            viewer::drawMesh(*mesh_);

        if (phase_ >= Phase::CrossField && phase_ < Phase::CutSeams && field_.has_value())
            viewer::drawField(*mesh_, *field_, scale_);

        if (phase_ >= Phase::Singularities && field_.has_value()) {
            double ballRadius = 0.5 * avgEdge_;
            for (const auto &sig : field_->uSingularities) {
                int vid = sig.first;
                int index4 = sig.second;
                if (vid < 0 || vid >= static_cast<int>(mesh_->vertices.size())) continue;
                const Point &c = mesh_->vertices[vid];
                if (index4 == 1)
                    viewer::drawDisk3D(c, ballRadius, 0.2f, 0.2f, 0.95f);
                else if (index4 == -1)
                    viewer::drawDisk3D(c, ballRadius, 0.95f, 0.2f, 0.2f);
            }
        }

        if (phase_ >= Phase::CutSeams && cutMesh_.has_value()) {
            viewer::drawUField(*mesh_, cutMesh_->getUField(), scale_);
            viewer::drawVField(*mesh_, cutMesh_->getVField(), scale_);
            if (!cutMesh_->getSingularityPathCutEdges().empty())
                viewer::drawEdgeSetOnMesh(*mesh_, cutMesh_->getCutEdges(), 1.0f, 0.75f, 0.1f, 4.0f);
            else
                viewer::drawEdgeSetOnMesh(*mesh_, cutMesh_->getCutEdges(), 1.0f, 0.2f, 0.9f, 3.5f);
        }
    }

    // Colour legend for the scalar field, before the text overlay so the
    // console keeps drawing on top.
    if (mode_ == Mode::OASIS && oasisPhase_ == OASISPhase::Field && oasis_.has_value()) {
        viewer::drawScalarFieldLegend(fbw(), fbh(), -oasisAbsMax_, oasisAbsMax_,
                                      "quasi-eigenfunction");
    }

    // Overlay text
    if (mode_ == Mode::Unselected) {
        renderOverlay("press '1' for PolyVector mode\npress '2' for MBO mode\n"
                      "press '3' for Medial Axis mode\npress '4' for SIPG mode\n"
                      "press '5' for OASIS mode\npress '6' for UMBER mode\n"
                      "right-drag to pan, scroll to zoom\n"
                      "press 'r' to restart\npress 'q' to quit");
    } else if (mode_ == Mode::OASIS) {
        renderOverlay("press 'c' to change lambda / orientation\n"
                      "press 'r' to restart\npress 'q' to quit");
    } else if (mode_ == Mode::UMBER && umberPhase_ == UMBERPhase::Simplified) {
        renderOverlay("press 'c' to change the collapse settings\n"
                      "press 'r' to restart\npress 'q' to quit");
    } else {
        renderOverlay("press 'c' to continue\npress 'r' to restart\npress 'q' to quit");
    }
}

// ── overlay helper ────────────────────────────────────────────────────────────

void CrossGenWidget::renderOverlay(const char *helpText) {
    int w = fbw(), h = fbh();
    console_.draw(w, h, 55.0f);
    viewer::drawTextOverlay(w, h, helpText, 10.0f, 20.0f, 0.8f, 0.8f, 0.8f);
}
