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
        case MBOPhase::Trace:        return MBOPhase::Trace;
    }
    return MBOPhase::Trace;
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

MedialAxisPhase nextMedialAxisPhase(MedialAxisPhase p) {
    switch (p) {
        case MedialAxisPhase::MeshOnly:     return MedialAxisPhase::DelaunayMesh;
        case MedialAxisPhase::DelaunayMesh: return MedialAxisPhase::MedialAxis;
        case MedialAxisPhase::MedialAxis:   return MedialAxisPhase::Classify;
        case MedialAxisPhase::Classify:     return MedialAxisPhase::Classify;
    }
    return MedialAxisPhase::Classify;
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

const char *medialAxisPhaseName(MedialAxisPhase p) {
    switch (p) {
        case MedialAxisPhase::MeshOnly:     return "1) mesh";
        case MedialAxisPhase::DelaunayMesh: return "2) Delaunay re-triangulation";
        case MedialAxisPhase::MedialAxis:   return "3) Medial axis";
        case MedialAxisPhase::Classify:     return "4) Classify / polylines";
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
              << " (press '1' for PolyVector, '2' for MBO, '3' for Medial Axis, '4' for SIPG, '5' for OASIS)\n";
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

    bool isSIPGStepping = (mode_ == Mode::SIPG &&
                           sipgPhase_ == SIPGPhase::Stepping &&
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
    bool inUVSplitScreen = (mode_ == Mode::SIPG && sipgPhase_ == SIPGPhase::UVMesh && sipgUVParam_.has_value()) ||
                           (mode_ == Mode::PolyVector && phase_ == Phase::UVMesh && miqSolver_.has_value());
    bool panRight = inUVSplitScreen && (lastMousePos_.x() * dpr > fbw() / 2);
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
    bool inUVSplitScreen = (mode_ == Mode::SIPG && sipgPhase_ == SIPGPhase::UVMesh && sipgUVParam_.has_value()) ||
                           (mode_ == Mode::PolyVector && phase_ == Phase::UVMesh && miqSolver_.has_value());
    bool zoomRight = inUVSplitScreen && (cx > fbw() / 2);
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
    delaunayMesh_.reset();
    medialAxis_.reset();
    oasis_.reset();
    oasisGuide_.reset();

    mode_     = Mode::Unselected;
    phase_    = Phase::MeshOnly;
    mboPhase_ = MBOPhase::MeshOnly;
    sipgPhase_ = SIPGPhase::MeshOnly;
    maPhase_  = MedialAxisPhase::MeshOnly;
    oasisPhase_ = OASISPhase::MeshOnly;
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
              << " (press '1' for PolyVector, '2' for MBO, '3' for Medial Axis, '4' for SIPG, '5' for OASIS)\n";
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
    }
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
        separatrixTrace_ = std::make_shared<SeparatrixTrace>(cfPtr, false);
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

    // ── SIPG: Initialize ──────────────────────────────────────────────────────
    if (mode_ == Mode::SIPG && sipgPhase_ >= SIPGPhase::CrossField && !sipgField_.has_value()) {
        auto t0 = Clock::now();
        sipgField_.emplace(mesh_);
        sipgField_->initialize();
        auto t1 = Clock::now();
        console_.log("[SIPG] Initialized: " +
                     formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count()));
    }

    // ── SIPG: Kick off stepping ───────────────────────────────────────────────
    if (mode_ == Mode::SIPG && sipgPhase_ == SIPGPhase::Stepping &&
        sipgField_.has_value() && !sipgSteppingStarted_) {
        sipgSteppingStarted_ = true;
        sipgStepCount_       = 0;
        console_.log("[SIPG] Starting MBO iterations (" +
                     std::to_string(sipgField_->getMesh().triangles.size()) + " tris)...");
    }

    // ── SIPG: Run 2 stepping iterations per frame ─────────────────────────────
    if (mode_ == Mode::SIPG && sipgPhase_ == SIPGPhase::Stepping &&
        sipgSteppingStarted_ && !sipgConverged_) {
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
    }

    // ── Medial Axis: polylines / classify ─────────────────────────────────────
    if (mode_ == Mode::MedialAxis && maPhase_ >= MedialAxisPhase::Classify &&
        medialAxis_ && medialAxis_->polyLines.empty()) {
        auto t0 = Clock::now();
        medialAxis_->createPolylines();
        medialAxis_->classifyMedialVertices();
        auto t1 = Clock::now();
        std::ostringstream oss;
        oss << "[MedialAxis] Created " << medialAxis_->polyLines.size() << " polylines: "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
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

// ── normal render ─────────────────────────────────────────────────────────────

void CrossGenWidget::renderNormal() {
    if (mode_ == Mode::PolyVector && phase_ == Phase::UVMesh && miqSolver_.has_value()) {
        // ── Split-screen: left = original mesh, right = UV mesh ───────────────
        int w = fbw(), h = fbh();
        int halfW = w / 2;

        // Helper: set viewport + ortho projection for a sub-rectangle
        auto applyHalfOrtho = [&](int x, int vpW, const viewer::ViewState &vs) {
            viewer::ViewState tmp = vs;
            tmp.fbw = vpW;
            tmp.fbh = h;
            double worldW = 1.0, worldH = 1.0;
            viewer::computeWorldBox(tmp, worldW, worldH);
            double left   = tmp.cx - 0.5 * worldW;
            double right  = tmp.cx + 0.5 * worldW;
            double bottom = tmp.cy - 0.5 * worldH;
            double top    = tmp.cy + 0.5 * worldH;
            glViewport(x, 0, vpW, h);
            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();
            glOrtho(left, right, bottom, top, -1, 1);
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
        };

        // ── Left panel: mesh with cut seams & singularities ───────────────────
        applyHalfOrtho(0, halfW, view_);
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

        // ── Dividing line (pixel-space) ───────────────────────────────────────
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

        // Restore full viewport for overlay
        glViewport(0, 0, w, h);
    } else if (mode_ == Mode::SIPG) {
        if (sipgPhase_ == SIPGPhase::UVMesh && sipgUVParam_.has_value()) {
            // ── Split-screen: left = mesh with cut seams, right = UVGParam ────
            int w = fbw(), h = fbh();
            int halfW = w / 2;

            auto applyHalfOrtho = [&](int x, int vpW, const viewer::ViewState &vs) {
                viewer::ViewState tmp = vs;
                tmp.fbw = vpW;
                tmp.fbh = h;
                double worldW = 1.0, worldH = 1.0;
                viewer::computeWorldBox(tmp, worldW, worldH);
                double left   = tmp.cx - 0.5 * worldW;
                double right  = tmp.cx + 0.5 * worldW;
                double bottom = tmp.cy - 0.5 * worldH;
                double top    = tmp.cy + 0.5 * worldH;
                glViewport(x, 0, vpW, h);
                glMatrixMode(GL_PROJECTION);
                glLoadIdentity();
                glOrtho(left, right, bottom, top, -1, 1);
                glMatrixMode(GL_MODELVIEW);
                glLoadIdentity();
            };

            // Left panel: mesh with combed field + cut seams
            applyHalfOrtho(0, halfW, view_);
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

            // Dividing line
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
        } else {
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

        if (mboPhase_ >= MBOPhase::Separatrices && separatrixTrace_) {
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
    } else if (mode_ == Mode::MedialAxis) {
        if (maPhase_ != MedialAxisPhase::Classify) {
            if (delaunayMesh_)
                viewer::drawMesh(*delaunayMesh_);
            else
                viewer::drawMesh(*mesh_);
        }

        if (medialAxis_) {
            if (maPhase_ == MedialAxisPhase::Classify) {
                if (delaunayMesh_)
                    viewer::drawBoundaryEdges(*delaunayMesh_);
                else
                    viewer::drawBoundaryEdges(*mesh_);

                glLineWidth(2.5f);
                glColor3f(0.95f, 0.85f, 0.1f);
                for (const auto &pl : medialAxis_->polyLines) {
                    if (pl.size() < 2) continue;
                    glBegin(GL_LINE_STRIP);
                    for (int vi : pl) {
                        const Point &p = medialAxis_->medialVertices[vi].coord;
                        glVertex2d(p[0], p[1]);
                    }
                    glEnd();
                }
                glLineWidth(1.0f);

                double ballRadius_ma = avgEdge_ / 5.0;
                for (size_t i = 0; i < medialAxis_->medialVertices.size(); ++i) {
                    const auto &mv = medialAxis_->medialVertices[i];
                    if (!mv.active || mv.degree == 2) continue;
                    if (mv.nodeType == TopMakerNodeType::Normal) {
                        viewer::drawDisk3D(mv.coord, ballRadius_ma, 0.2f, 0.2f, 0.95f);
                    } else if (mv.nodeType == TopMakerNodeType::Corner) {
                        if (mv.cornerIndex >= 0 && mv.cornerIndex < static_cast<int>(medialAxis_->mesh->vertices.size())) {
                            const Point &cornerP = medialAxis_->mesh->vertices[mv.cornerIndex];
                            viewer::drawDisk3D(cornerP, ballRadius_ma, 0.95f, 0.2f, 0.2f);
                            glLineWidth(2.0f);
                            glColor3f(0.95f, 0.85f, 0.1f);
                            glBegin(GL_LINES);
                            glVertex2d(cornerP[0], cornerP[1]);
                            glVertex2d(mv.coord[0], mv.coord[1]);
                            glEnd();
                            glLineWidth(1.0f);
                        }
                    } else {
                        viewer::drawDisk3D(mv.coord, ballRadius_ma, 0.1f, 0.9f, 0.2f);
                    }
                }

                double touchRadius = avgEdge_ / 6.0;
                for (size_t i = 0; i < medialAxis_->medialVertices.size(); ++i) {
                    const auto &mv = medialAxis_->medialVertices[i];
                    if (!mv.active || mv.degree == 2) continue;
                    if (mv.nodeType != TopMakerNodeType::Normal) continue;
                    for (int tpIdx : mv.touchPoints) {
                        if (tpIdx >= 0 && tpIdx < static_cast<int>(medialAxis_->mesh->vertices.size())) {
                            const Point &tp = medialAxis_->mesh->vertices[tpIdx];
                            viewer::drawDisk3D(tp, touchRadius, 0.0f, 0.9f, 0.9f);
                        }
                    }
                }
            } else if (maPhase_ >= MedialAxisPhase::MedialAxis) {
                double ballRadius_ma = avgEdge_ / 5.0;
                viewer::drawMedialAxis(*medialAxis_, ballRadius_ma);
                for (int i = 0; i < static_cast<int>(medialAxis_->sharpVertices.size()); ++i) {
                    if (medialAxis_->sharpVertices[i]) {
                        const Point &p = medialAxis_->mesh->vertices[i];
                        viewer::drawDisk3D(p, ballRadius_ma, 0.95f, 0.2f, 0.2f);
                    }
                }
            }
        }
    } else {
        // PolyVector phases 1-4 (not UV)
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
                      "press '5' for OASIS mode\n"
                      "right-drag to pan, scroll to zoom\n"
                      "press 'r' to restart\npress 'q' to quit");
    } else if (mode_ == Mode::OASIS) {
        renderOverlay("press 'c' to change lambda / orientation\n"
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
