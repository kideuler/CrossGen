// CrossGenWidget.cxx – Qt6 QOpenGLWidget implementation of the CrossGen viewer.
// Replaces the old GLFW-based ViewerMain render loop with Qt event-driven rendering.

#include "viewer/CrossGenWidget.hxx"
#include "viewer/Export.hxx"
#include "viewer/GL.hxx"
#include "viewer/Interaction.hxx"
#include "viewer/Render.hxx"

#include <QCheckBox>
#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFileInfo>
#include <QFormLayout>
#include <QImage>
#include <QKeyEvent>
#include <QOpenGLContext>
#include <QOpenGLFramebufferObject>
#include <QOpenGLFunctions>
#include <QPageLayout>
#include <QPageSize>
#include <QPainter>
#include <QPdfWriter>
#include <QSvgRenderer>
#include <QLabel>
#include <QPushButton>
#include <QScreen>
#include <QScrollArea>
#include <QSpinBox>
#include <QVBoxLayout>
#include <QMouseEvent>
#include <QWheelEvent>
#include <QTimer>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
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
        case Phase::CutSeams:      return Phase::CutSeams;
    }
    return Phase::CutSeams;
}

ZIPLINEPhase nextZIPLINEPhase(ZIPLINEPhase p) {
    switch (p) {
        case ZIPLINEPhase::MeshOnly:     return ZIPLINEPhase::CrossField;
        case ZIPLINEPhase::CrossField:   return ZIPLINEPhase::Stepping;
        case ZIPLINEPhase::Stepping:     return ZIPLINEPhase::Separatrices;
        case ZIPLINEPhase::Separatrices: return ZIPLINEPhase::Trace;
        case ZIPLINEPhase::Trace:        return ZIPLINEPhase::Layout;
        case ZIPLINEPhase::Layout:       return ZIPLINEPhase::Simplified;
        case ZIPLINEPhase::Simplified:   return ZIPLINEPhase::Quantize;
        case ZIPLINEPhase::Quantize:     return ZIPLINEPhase::Quantized;
        case ZIPLINEPhase::Quantized:    return ZIPLINEPhase::Blocks;
        case ZIPLINEPhase::Blocks:       return ZIPLINEPhase::Mesh;
        case ZIPLINEPhase::Mesh:         return ZIPLINEPhase::Smoothed;
        case ZIPLINEPhase::Smoothed:     return ZIPLINEPhase::Smoothed;
    }
    return ZIPLINEPhase::Smoothed;
}

UMBERPhase nextUMBERPhase(UMBERPhase p) {
    switch (p) {
        case UMBERPhase::MeshOnly:   return UMBERPhase::CrossField;
        case UMBERPhase::CrossField: return UMBERPhase::Stepping;
        case UMBERPhase::Stepping:   return UMBERPhase::Frames;
        case UMBERPhase::Frames:     return UMBERPhase::Polysquare;
        case UMBERPhase::Polysquare: return UMBERPhase::Blocks;
        case UMBERPhase::Blocks:     return UMBERPhase::Decomposition;
        case UMBERPhase::Decomposition: return UMBERPhase::Mesh;
        case UMBERPhase::Mesh:       return UMBERPhase::Smoothed;
        case UMBERPhase::Smoothed:   return UMBERPhase::Smoothed;
    }
    return UMBERPhase::Smoothed;
}

PipelinePhase nextPipelinePhase(PipelinePhase p) {
    switch (p) {
        case PipelinePhase::MeshOnly:   return PipelinePhase::CrossField;
        case PipelinePhase::CrossField: return PipelinePhase::Stepping;
        case PipelinePhase::Stepping:   return PipelinePhase::Cones;
        case PipelinePhase::Cones:      return PipelinePhase::Cut;
        case PipelinePhase::Cut:        return PipelinePhase::Flow;
        case PipelinePhase::Flow:  return PipelinePhase::Metric;
        case PipelinePhase::Metric:     return PipelinePhase::Layout;
        case PipelinePhase::Layout:     return PipelinePhase::Separatrices;
        case PipelinePhase::Separatrices: return PipelinePhase::Patches;
        case PipelinePhase::Patches:    return PipelinePhase::Mesh;
        case PipelinePhase::Mesh:       return PipelinePhase::Smoothed;
        case PipelinePhase::Smoothed:   return PipelinePhase::Smoothed;
    }
    return PipelinePhase::Smoothed;
}

ATLASPhase nextATLASPhase(ATLASPhase p) {
    switch (p) {
        case ATLASPhase::MeshOnly: return ATLASPhase::Domain;
        case ATLASPhase::Domain:   return ATLASPhase::Field;
        case ATLASPhase::Field:    return ATLASPhase::Carrier;
        case ATLASPhase::Carrier:  return ATLASPhase::Search;
        case ATLASPhase::Search:   return ATLASPhase::Blocks;
        case ATLASPhase::Blocks:   return ATLASPhase::Mesh;
        case ATLASPhase::Mesh:     return ATLASPhase::Smoothed;
        case ATLASPhase::Smoothed: return ATLASPhase::Smoothed;
    }
    return ATLASPhase::Smoothed;
}

MedialAxisPhase nextMedialAxisPhase(MedialAxisPhase p) {
    switch (p) {
        case MedialAxisPhase::MeshOnly:     return MedialAxisPhase::DelaunayMesh;
        case MedialAxisPhase::DelaunayMesh: return MedialAxisPhase::MedialAxis;
        case MedialAxisPhase::MedialAxis:   return MedialAxisPhase::Map;
        case MedialAxisPhase::Map:          return MedialAxisPhase::TMesh;
        case MedialAxisPhase::TMesh:        return MedialAxisPhase::Quantize;
        case MedialAxisPhase::Quantize:     return MedialAxisPhase::Quantized;
        case MedialAxisPhase::Quantized:    return MedialAxisPhase::Quantized;
    }
    return MedialAxisPhase::Quantized;
}

const char *phaseName(Phase p) {
    switch (p) {
        case Phase::MeshOnly:      return "1) mesh";
        case Phase::CrossField:    return "2) crossfield";
        case Phase::Singularities: return "3) singularities";
        case Phase::CutSeams:      return "4) cut seams";
    }
    return "?";
}

const char *ziplinePhaseName(ZIPLINEPhase p) {
    switch (p) {
        case ZIPLINEPhase::MeshOnly:     return "1) mesh";
        case ZIPLINEPhase::CrossField:   return "2) MBO crossfield";
        case ZIPLINEPhase::Stepping:     return "3) MBO stepping";
        case ZIPLINEPhase::Separatrices: return "4) separatrices";
        case ZIPLINEPhase::Trace:        return "5) trace";
        case ZIPLINEPhase::Layout:       return "6) quad layout";
        case ZIPLINEPhase::Simplified:   return "7) simplified partition";
        case ZIPLINEPhase::Quantize:     return "8) Quantization";
        case ZIPLINEPhase::Quantized:    return "9) Quantized block decomposition";
        case ZIPLINEPhase::Blocks:       return "10) block decomposition";
        case ZIPLINEPhase::Mesh:         return "11) quad mesh";
        case ZIPLINEPhase::Smoothed:     return "12) TMOP smoothing";
    }
    return "?";
}

const char *umberPhaseName(UMBERPhase p) {
    switch (p) {
        case UMBERPhase::MeshOnly:   return "1) mesh";
        case UMBERPhase::CrossField: return "2) DualMBO crossfield";
        case UMBERPhase::Stepping:   return "3) DualMBO stepping";
        case UMBERPhase::Frames:     return "4) UMBER frame field";
        case UMBERPhase::Polysquare: return "5) polysquare (Sec. 4.3)";
        case UMBERPhase::Blocks:     return "6) block structure (Sec. 5)";
        case UMBERPhase::Decomposition: return "7) block decomposition";
        case UMBERPhase::Mesh:       return "8) quad mesh";
        case UMBERPhase::Smoothed:   return "9) TMOP smoothing";
    }
    return "?";
}

// The two pipelines differ in exactly two of the eleven rows, so the name is a
// function of the phase and the mode rather than of the phase alone. Every
// other row is the same stage of the same paper reached by the same code.
const char *pipelinePhaseName(PipelinePhase p, Mode m) {
    const bool field = (m == Mode::TORSION);
    switch (p) {
        case PipelinePhase::MeshOnly:   return "1) mesh";
        case PipelinePhase::CrossField: return "2) DualMBO crossfield";
        case PipelinePhase::Stepping:   return "3) DualMBO stepping";
        case PipelinePhase::Cones:      return "4) cone singularities (Sec. 3.1)";
        case PipelinePhase::Cut:        return "5) cutting graph (Sec. 3.2.2)";
        case PipelinePhase::Flow:
            return field ? "6) combed field and matchings (Stage 3F)"
                         : "6) discrete Ricci flow (Sec. 3.2.1)";
        case PipelinePhase::Metric:
            return field ? "7) psi_0 by integration (Stages 4F, 4R)"
                         : "7) flat cone metric";
        case PipelinePhase::Layout:     return "8) layout Psi (Secs. 3.2.2, 3.3)";
        case PipelinePhase::Separatrices: return "9) separatrices (Sec. 4, Q5)";
        case PipelinePhase::Patches:    return "10) arrangement and splines (Secs. 4, 5)";
        case PipelinePhase::Mesh:       return "11) quadrilateral mesh (Sec. 5)";
        case PipelinePhase::Smoothed:   return "12) TMOP smoothing (mesh::TMOP)";
    }
    return "?";
}

// The last two are the pipelines' Mesh and Smoothed phases in all but number.
const char *atlasPhaseName(ATLASPhase p) {
    switch (p) {
        case ATLASPhase::MeshOnly: return "1) mesh";
        case ATLASPhase::Domain:   return "2) planar domain (Stage 1, Sec. 1.1)";
        case ATLASPhase::Field:    return "3) reference cross field (Stage 1b, guidance note Sec. 3)";
        case ATLASPhase::Carrier:  return "4) three-quad carrier (Stage 2, Sec. 3)";
        case ATLASPhase::Search:   return "5) search and rewrite (Stages 3-6, Secs. 7, 8)";
        case ATLASPhase::Blocks:   return "6) blocks on the input (Stages 4-5, Secs. 5, 6)";
        case ATLASPhase::Mesh:     return "7) quadrilateral mesh (Secs. 11.2, 11.3)";
        case ATLASPhase::Smoothed: return "8) TMOP smoothing (mesh::TMOP)";
    }
    return "?";
}

const char *medialAxisPhaseName(MedialAxisPhase p) {
    switch (p) {
        case MedialAxisPhase::MeshOnly:     return "1) mesh";
        case MedialAxisPhase::DelaunayMesh: return "2) Delaunay re-triangulation";
        case MedialAxisPhase::MedialAxis:   return "3) Medial axis";
        case MedialAxisPhase::Map:          return "4) Boundary map";
        case MedialAxisPhase::TMesh:        return "5) Coarse block decomposition";
        case MedialAxisPhase::Quantize:     return "6) Quantization";
        case MedialAxisPhase::Quantized:    return "7) Quantized block decomposition";
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
        case Mode::ZIPLINE:    return "ZIPLINE";
        case Mode::MedialAxis: return "Medial Axis";
        case Mode::TORSION:    return "TORSION";
        case Mode::OASIS:      return "OASIS";
        case Mode::UMBER:      return "UMBER";
        case Mode::MERIDIAN:   return "MERIDIAN";
        case Mode::ATLAS:      return "ATLAS";
    }
    return "?";
}

// The line every mode-selection prompt prints, kept in one place so adding a
// mode does not mean chasing three copies of it.
const char *kModeMenu =
    "press '1' for PolyVector, '2' for ZIPLINE, '3' for Medial Axis, '4' for TORSION, "
    "'5' for OASIS, '6' for UMBER, '7' for MERIDIAN, '8' for ATLAS";

} // anonymous namespace

// ── constructor ───────────────────────────────────────────────────────────────

CrossGenWidget::CrossGenWidget(const std::string &meshPath, QWidget *parent)
    : QOpenGLWidget(parent)
{
    mesh_ = std::make_shared<Mesh>(meshPath);

    // Where the figures go, and what they are called. The directory is only
    // created when the first one is written, so a session that never presses
    // 's' leaves nothing behind.
    modelName_ = QFileInfo(QString::fromStdString(meshPath)).completeBaseName()
                     .toStdString();
    if (modelName_.empty()) modelName_ = "model";
    figureDir_ = qEnvironmentVariableIsSet("CROSSGEN_FIGURE_DIR")
                     ? qEnvironmentVariable("CROSSGEN_FIGURE_DIR")
                     : QStringLiteral("figures");

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

    std::cerr << "[Viewer] Phase " << phaseName(phase_) << " (" << kModeMenu << ")\n";
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
    // Deliberately outside drawScene(): an export re-draws the scene, and a
    // stage must not be advanced or re-solved by having its picture saved.
    runComputations();

    // Clear, in whichever background the theme is on. White is the default,
    // because the figures are what this is for.
    const auto bg = viewer::backgroundColor();
    glClearColor(bg[0], bg[1], bg[2], 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);

    drawScene();
}

void CrossGenWidget::drawScene() {
    // Choose render path
    bool isMBOStepping = (mode_ == Mode::ZIPLINE &&
                          ziplinePhase_ == ZIPLINEPhase::Stepping &&
                          mboSteppingStarted_ &&
                          ziplineField() && !zipline_->fieldFinished());

    bool isDualMBOStepping = (dualMBOStageIsStepping() &&
                           dualMBOSteppingStarted_ &&
                           !dualMBOConverged_);

    bool isTracing = (mode_ == Mode::ZIPLINE &&
                      ziplinePhase_ == ZIPLINEPhase::Trace &&
                      ziplineTrace() &&
                      !ziplineTracingFinished_);

    if (isMBOStepping) {
        renderMBOAnimation();
    } else if (isDualMBOStepping) {
        renderDualMBOAnimation();
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

    case Qt::Key_N:
        // The connectivity dialog of Stages 5 to 7. It used to sit on 'c' at
        // the Patches phase, which was the last one; Stage 10 is now, and 'c'
        // there re-opens the mesh dialog. Which pairs of cones are meant to be
        // joined is still a judgement about the model that only trying a number
        // settles, so it keeps a key of its own.
        if (inPipeline() && pipePhase_ >= PipelinePhase::Patches &&
            immersion_.has_value() && separatricesAttempted_) {
            if (promptMERIDIANConnectivity(*immersion_, psiR_)) {
                rerunMERIDIANConnectivity();
                // Stages 8 and 9 are gone with the trace, and so is the mesh
                // that stood on them; step back to the phase that rebuilds
                // them rather than drawing over a stale one.
                pipePhase_ = PipelinePhase::Patches;
            }
        }
        break;

    case Qt::Key_E:
        // The Stage 10 dialog again. It used to sit on 'c' at the Mesh phase,
        // which was the last one; Stage 12 is now, and 'c' there re-opens the
        // TMOP dialog. A target edge length is still a judgement about the
        // model that only trying a number settles, so it keeps a key of its
        // own -- the same move the connectivity dialog made onto 'n' when this
        // phase displaced it.
        //
        // Offered from Stage 12 as well as from Stage 10, because that is where
        // the reason to re-mesh usually shows up: a target the smoother cannot
        // rescue is a target, and not a smoother, that was wrong.
        if (inPipeline() && pipePhase_ >= PipelinePhase::Mesh && splines_.has_value()) {
            meshAttempted_ = true;
            if (promptMERIDIANMesh()) {
                runMERIDIANMesh();
                // The mesh Stage 12 smoothed is gone with it. Step back to the
                // phase that shows the new one rather than leaving an empty
                // Stage 12 on screen; 'c' from there smooths it.
                pipePhase_ = PipelinePhase::Mesh;
            }
        }
        // ATLAS's mesh the same way, off the same key, for the same reasons.
        if (mode_ == Mode::ATLAS && atlasPhase_ >= ATLASPhase::Mesh && atlas_ &&
            atlas_->hasCover()) {
            atlasMeshAttempted_ = true;
            if (promptATLASMesh()) {
                runATLASMesh();
                atlasPhase_ = ATLASPhase::Mesh;
            }
        }
        // And UMBER's.
        if (mode_ == Mode::UMBER && umberPhase_ >= UMBERPhase::Mesh &&
            umberDecomp_.has_value() && !umberDecomp_->blocks.empty()) {
            umberMeshAttempted_ = true;
            if (promptUMBERMesh()) {
                runUMBERMesh();
                umberPhase_ = UMBERPhase::Mesh;
            }
        }
        // And ZIPLINE's.
        if (mode_ == Mode::ZIPLINE && ziplinePhase_ >= ZIPLINEPhase::Mesh && ziplineBlocks() &&
            !ziplineBlocks()->blocks().empty()) {
            ziplineMeshAttempted_ = true;
            if (promptZIPLINEMesh()) {
                runZIPLINEMesh();
                ziplinePhase_ = ZIPLINEPhase::Mesh;
            }
        }
        break;

    case Qt::Key_I:
        // The interface network, on or off. It is drawn over every phase of
        // either pipeline, which is what makes it useful and also what makes a
        // key to hide it necessary: at the mesh phase it lies on top of the
        // elements it is there to be compared against.
        if ((inPipeline() || mode_ == Mode::UMBER || mode_ == Mode::ZIPLINE) && interfaceNetwork() &&
            interfaceNetwork()->multiMaterial()) {
            showInterfaces_ = !showInterfaces_;
            console_.log(showInterfaces_ ? "[Interfaces] network shown"
                                         : "[Interfaces] network hidden");
            update();
        }
        break;

    case Qt::Key_M:
        // The materials as a fill rather than as the colour of a wireframe
        // edge. Off by default because it competes with the conformal factor
        // and the flat metric for the same triangles.
        if ((inPipeline() && interfaces_.has_value() && interfaces_->multiMaterial()) ||
            (mode_ == Mode::ATLAS && atlasMultiMaterial_) ||
            (mode_ == Mode::UMBER && umberDecompReport_.materials > 1) ||
            (mode_ == Mode::ZIPLINE && ziplineBlocks() &&
             ziplineBlocks()->getReport().materials > 1)) {
            showMaterialFill_ = !showMaterialFill_;
            if (mode_ == Mode::ATLAS)
                console_.log(showMaterialFill_ ? "[ATLAS] cells and elements filled by material"
                                               : "[ATLAS] material fill off");
            else if (mode_ == Mode::UMBER)
                console_.log(showMaterialFill_ ? "[UMBER] elements filled by material"
                                               : "[UMBER] material fill off");
            else if (mode_ == Mode::ZIPLINE)
                console_.log(showMaterialFill_ ? "[Blocks] elements filled by material"
                                               : "[Blocks] material fill off");
            else
                console_.log(showMaterialFill_ ? "[Interfaces] triangles filled by material"
                                               : "[Interfaces] material fill off");
            update();
        }
        break;

    case Qt::Key_P:
        // ATLAS's Search phase is the same kind of pair: the carrier the
        // winning search started from and the one it ended on, and what Stages
        // 3 and 6 did is the difference.
        if (mode_ == Mode::ATLAS && atlasPhase_ == ATLASPhase::Search && atlas_ &&
            atlas_->hasCover()) {
            atlasShowInitial_ = !atlasShowInitial_;
            console_.log(atlasShowInitial_
                             ? "[Search] the carrier the search started from: Stage 2 on its domain"
                             : "[Search] the carrier the search ended on, with the blocks it found");
            update();
            break;
        }
        // Both phases that carry two maps of the same domain in the right half
        // put the swap on this key, and for the same reason: what the stage did
        // is the difference between them, and a difference between two pictures
        // is only visible by putting one where the other was.
        //
        // At TORSION's Metric phase the pair is Stage 4F's least-squares map
        // and what Stage 4R made of it. The whole cost of substituting an
        // integration for a flow is the red faces on the first of them, so
        // hiding that behind the repaired map would hide the one measurement
        // the phase exists to show.
        if (mode_ == Mode::TORSION && pipePhase_ == PipelinePhase::Metric &&
            immersion_.has_value() && !integratedMap_.empty()) {
            showIntegrated_ = !showIntegrated_;
            const std::vector<Point> &uv =
                showIntegrated_ ? integratedMap_ : immersion_->getUV();
            viewer::computeLayoutBounds(uv, uvView_.cx, uvView_.cy,
                                        uvView_.baseW, uvView_.baseH);
            uvView_.zoom = 1.0;
            console_.log(showIntegrated_
                             ? "[psi_0] right panel: the Stage 4F least-squares map, "
                               "red where it inverted"
                             : "[psi_0] right panel: psi_0 as Stage 4 accepted it");
            update();
            break;
        }
        // The right half of the layout phase carries Stage 4's psi_R and Stage
        // 6's Psi, and the whole of what Sec. 3.3 does is the difference
        // between them.
        if (inPipeline() && pipePhase_ == PipelinePhase::Layout &&
            meridianLayout_.has_value() && !psiR_.empty()) {
            showPsiR_ = !showPsiR_;
            const std::vector<Point> &uv = showPsiR_ ? psiR_ : meridianLayout_->getUV();
            viewer::computeLayoutBounds(uv, uvView_.cx, uvView_.cy,
                                        uvView_.baseW, uvView_.baseH);
            uvView_.zoom = 1.0;
            console_.log(showPsiR_ ? "[Layout] right panel: psi_R, the Stage 4 immersion"
                                   : "[Layout] right panel: Psi, the Stage 6 layout");
            update();
        }
        break;

    case Qt::Key_S:
        // The figure keys. 's' is the vector one and the one to reach for: a
        // wireframe, a cut graph, a quad mesh -- everything this viewer draws
        // is lines and flat polygons, which is exactly what a vector page is
        // for, and a paper renders it at whatever resolution the press is. The
        // file is a PDF, which is what a LaTeX \includegraphics wants.
        //
        // 'S' is the raster fallback at 3x the framebuffer, for the phases
        // whose picture is a filled scalar field over every triangle: those
        // are vector files of a hundred thousand paths that no viewer wants to
        // open and that gain nothing over pixels.
        if (event->modifiers() & Qt::ShiftModifier) exportPng(3);
        else                                        exportPdf();
        break;

    case Qt::Key_H:
        // Everything on screen that belongs to the viewer rather than to the
        // model: the console, the key help, the coordinate axis and the colour
        // legends. One key for all four, because a figure wants none of them
        // and a working session wants all of them.
        showHUD_ = !showHUD_;
        viewer::setAxisVisible(showHUD_);
        console_.log(showHUD_ ? "[Figure] console, key help, axis and legends shown"
                              : "[Figure] console, key help, axis and legends hidden");
        update();
        break;

    case Qt::Key_B:
        // White or the dark background the viewer was written against. White is
        // the default; this is here because the dark one is easier to work
        // against for long stretches, and because a figure ought to be
        // checkable against it.
        viewer::setLightBackground(!viewer::lightBackground());
        console_.log(viewer::lightBackground() ? "[Figure] white background"
                                               : "[Figure] dark background");
        update();
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
            mode_ = Mode::ZIPLINE;
            std::cerr << "[Viewer] Selected mode: " << modeName(mode_) << " (press 'c' to advance)\n";
            console_.log("Selected mode: ZIPLINE");
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
            // Stage 0c first, and *before* mode_ is set. It replaces the mesh
            // every stage of the pipeline is built on, and it opens a dialog to
            // decide whether to: a modal dialog runs a nested event loop, that
            // loop paints, and painting is what drives runComputations(). With
            // the mode already set, Stage 0b would be built on the mesh as it
            // stands -- the one about to be thrown away -- while the dialog is
            // still open, and every stage after it would refuse to work with an
            // interface network found on a different mesh.
            runDiskExcision();
            mode_ = Mode::TORSION;
            std::cerr << "[Viewer] Selected mode: " << modeName(mode_) << " (press 'c' to advance)\n";
            console_.log("Selected mode: TORSION (psi_0 by integrating the cross field)");
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

    case Qt::Key_7:
        if (mode_ == Mode::Unselected && phase_ == Phase::MeshOnly) {
            runDiskExcision();   // before mode_, and for the reason given under '4'
            mode_ = Mode::MERIDIAN;
            std::cerr << "[Viewer] Selected mode: " << modeName(mode_) << " (press 'c' to advance)\n";
            console_.log("Selected mode: MERIDIAN (Shepherd, Gu and Hughes 2022, Stages 1-6)");
        }
        break;

    case Qt::Key_8:
        // No Stage 0c: ATLAS lays a circular inclusion out like any other
        // region (Sec. 7.4's annulus is one of its templates), so the mesh is
        // taken as it was loaded.
        if (mode_ == Mode::Unselected && phase_ == Phase::MeshOnly) {
            mode_ = Mode::ATLAS;
            atlasMultiMaterial_ = false;
            for (size_t t = 1; t < mesh_->triangleMatId.size() && !atlasMultiMaterial_; ++t)
                atlasMultiMaterial_ = mesh_->triangleMatId[t] != mesh_->triangleMatId[0];
            std::cerr << "[Viewer] Selected mode: " << modeName(mode_) << " (press 'c' to advance)\n";
            console_.log("Selected mode: ATLAS (square-transport blocking, Stages 1-6, and a TFI mesh)");
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
    dualMBOField_.reset();
    ziplineQuant_.reset();
    ziplineQuantReport_ = TMeshQuantizer::Report{};
    zipline_.reset();
    delaunayMesh_.reset();
    blockQuant_.reset();
    quantReport_ = TMeshQuantizer::Report{};
    medialAxisTMesh_.reset();
    medialAxisMap_.reset();
    medialAxis_.reset();
    oasis_.reset();
    oasisGuide_.reset();
    umber_.reset();
    umberCut_.reset();
    polysquare_.reset();
    umberMesh_.reset();
    umberDecomp_.reset();
    umberDecompReport_ = BlockLayout::DecompositionReport{};
    blockLayout_.reset();
    blocks_.reset();
    umberCorners_.clear();
    umberInternal_.clear();
    // Each stage of either pipeline holds on to the ones before it, so they go
    // in reverse: the layout on the labels and the immersion, the immersion on
    // the cut, the flow and the cones, and both the flow and the cut on the
    // cones.
    diskFill_.reset();
    smoothMesh_.reset();
    pipelineBlocked_.clear();
    // ATLAS's, in reverse of how each stands on the one before it.
    atlasMesh_.reset();
    atlas_.reset();
    atlasCarrier_.reset();
    atlasField_.reset();
    atlasDomain_.reset();
    atlasShowInitial_ = false;
    atlasMultiMaterial_ = false;
    // Stage 0c replaced the mesh every later stage was written on, so the reset
    // has to put the loaded one back before anything is built on it again.
    if (inputMesh_) mesh_ = inputMesh_;
    inputMesh_.reset();
    inclusions_.clear();
    quadMesh_.reset();
    splines_.reset();
    arrangement_.reset();
    separatrices_.reset();
    meridianLayout_.reset();
    meridianLabels_.reset();
    immersion_.reset();
    psiR_.clear();
    showPsiR_ = false;
    ricci_.reset();
    // Pipeline B's three, and the scaffold Immersion that also points at the
    // cut and the cones.
    integration_.reset();
    tutte_.reset();
    scaffold_.reset();
    frames_.reset();
    integratedMap_.clear();
    showIntegrated_ = false;
    fieldIndex_.clear();
    coneCut_.reset();
    cones_.reset();
    // Last of the pipeline stages to go: SubdomainLabels holds a bare pointer
    // to it.
    interfaces_.reset();
    interfaceMessagesSeen_ = 0;
    flatMetric_ = viewer::FlatMetric{};
    coneFans_.clear();
    ricciU_.resize(0);
    ricciUAbsMax_ = 1.0;

    mode_     = Mode::Unselected;
    phase_    = Phase::MeshOnly;
    ziplinePhase_ = ZIPLINEPhase::MeshOnly;
    maPhase_  = MedialAxisPhase::MeshOnly;
    oasisPhase_ = OASISPhase::MeshOnly;
    umberPhase_ = UMBERPhase::MeshOnly;
    pipePhase_ = PipelinePhase::MeshOnly;
    atlasPhase_ = ATLASPhase::MeshOnly;
    // oasisLambda_ deliberately survives a reset so it can be reused as the
    // dialog's default on the next run.

    singularitiesLogged_  = false;
    mboSteppingStarted_   = false;
    ziplineTracingStarted_    = false;
    ziplineTracingFinished_   = false;
    dualMBOSteppingStarted_  = false;
    dualMBOConverged_        = false;
    dualMBOStepCount_        = 0;
    umberAnnounced_       = false;
    umberAttempted_       = false;
    umberMeshAttempted_   = false;
    ziplineMeshAttempted_   = false;
    polysquareAnnounced_  = false;
    polysquareAttempted_  = false;
    blocksAttempted_      = false;
    interfacesAttempted_  = false;
    conesAttempted_       = false;
    cutAttempted_         = false;
    ricciAnnounced_       = false;
    ricciAttempted_       = false;
    framesAttempted_      = false;
    integrationAnnounced_ = false;
    integrationAttempted_ = false;
    layoutAnnounced_      = false;
    layoutAttempted_      = false;
    separatricesAnnounced_ = false;
    separatricesAttempted_ = false;
    patchesAnnounced_      = false;
    patchesAttempted_      = false;
    nearMissRetried_       = false;
    unmeshableBefore_      = -1.0;
    meshAttempted_         = false;
    tmopAttempted_         = false;
    disksAttempted_        = false;
    atlasDomainAttempted_  = false;
    atlasFieldAnnounced_   = false;
    atlasFieldAttempted_   = false;
    atlasCarrierAttempted_ = false;
    atlasAnnounced_        = false;
    atlasAttempted_        = false;
    atlasBlocksLogged_     = false;
    atlasMeshAttempted_    = false;

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
    std::cerr << "[Viewer] Reset. Phase " << phaseName(phase_) << " (" << kModeMenu << ")\n";
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
    if (!dualMBOField_.has_value()) return;

    // Eq. (1) starts from a *converged* cross field: the comb in initialize()
    // assumes neighbouring triangles already agree up to a k*90-degree turn,
    // which a half-solved MBO field does not. Advancing out of the stepping
    // phase early therefore finishes the solve here rather than optimizing a
    // field that is still moving.
    if (!dualMBOConverged_) {
        auto t0 = Clock::now();
        dualMBOField_->runMBO();
        dualMBOField_->computeSingularities();
        dualMBOConverged_ = true;
        auto t1 = Clock::now();
        std::ostringstream oss;
        oss << "[UMBER] finished the DualMBO solve first, error " << std::scientific
            << std::setprecision(3) << dualMBOField_->error << ", "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }

    // --- Cuts, Sec. 4.1 -----------------------------------------------------
    auto tc0 = Clock::now();
    try {
        // No cuts: the holes stay holes and the parameterization is a common
        // polysquare. See the HarmonicCut constructor for why a block
        // decomposition wants that and the closed form of Sec. 4.1 costs it
        // the ring it cannot represent.
        umberCut_.emplace(mesh_, /*openVoids=*/false);
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
        oss << "[UMBER] " << rep.voids << " hole(s), left uncut, chi "
            << rep.eulerCharacteristic << ", "
            << formatMs(std::chrono::duration<double, std::milli>(tc1 - tc0).count());
        console_.log(oss.str());
    }

    // --- Frame field, Eq. (1) ----------------------------------------------
    UMBER::EnergyTerms before;
    size_t internalBefore = 0;
    auto t0 = Clock::now();
    try {
        umber_.emplace(*dualMBOField_, *umberCut_);
        // The interfaces are features of Eq. (4) in exactly the sense dS is;
        // everything after this -- the deformation, the iso-lines, the blocks
        // -- reads them back off the frame field rather than being told again.
        if (interfaces_.has_value() && interfaces_->multiMaterial())
            umber_->setFeatureEdges(interfaces_->interfaceEdges());
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

    // This points into blocks_, so it goes before it is replaced, and
    // everything read off it goes with it.
    blockLayout_.reset();
    umberDecomp_.reset();
    umberMesh_.reset();
    umberMeshAttempted_ = false;

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

// ── the blocks as a graph, and as the shared decomposition ───────────────────

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
    // open: it cannot become a block of the decomposition, so the part of the
    // model it covers is the part the mesh will be missing.
    if (r.badFaces > 0 || r.arcCrossings > 0) {
        std::ostringstream bad;
        bad << "[BlockLayout] " << r.badFaces << " block(s) not four-sided, " << r.arcCrossings
            << " side(s) crossing";
        console_.log(bad.str());
    }
    buildUMBERDecomposition();
}

// ── the structure as the one representation all three methods share ─────────
//
// The traced layout read as a BlockDecomposition. Everything after this point
// reads that and not the layout: the picture, the mesh and the smoothed mesh
// are then all of the same structure by construction rather than by three
// routines agreeing about it.
void CrossGenWidget::buildUMBERDecomposition() {
    umberDecomp_.reset();
    umberDecompReport_ = BlockLayout::DecompositionReport{};
    // The mesh standing on the old structure goes with it.
    umberMesh_.reset();
    smoothMesh_.reset();
    umberMeshAttempted_ = false;
    tmopAttempted_ = false;
    if (!blockLayout_.has_value() || !mesh_) return;

    umberDecomp_.emplace(
        BlockLayout::blockDecompositionOf(blockLayout_->getLayout(), *mesh_,
                                          &umberDecompReport_, "UMBER"));

    const BlockLayout::DecompositionReport &d = umberDecompReport_;
    std::ostringstream oss;
    oss << "[Decomposition] " << d.blocks << " block(s) of " << d.faces << " face(s), "
        << umberDecomp_->edges.size() << " macro edge(s), " << umberDecomp_->vertices.size()
        << " macrovertex/-ices, " << d.materials << " material(s)";
    if (d.notFourSided > 0) oss << ", " << d.notFourSided << " not four-sided";
    if (d.multiArcSides > 0) oss << ", " << d.multiArcSides << " with a split side";
    console_.log(oss.str());

    // The one number that says whether the picture is of the model or of part
    // of it. A refused face is not drawn in another colour, it is a piece of
    // the model with no block on it, and the mesh will have a hole exactly
    // there -- so this is said every time, and said loudly when it is not all
    // of the model.
    const double coverage = BlockLayout::coverageOf(d);
    std::ostringstream cov;
    cov << std::fixed << std::setprecision(1) << "[Decomposition] covering "
        << 100.0 * coverage << "% of the model";
    if (coverage < 0.999)
        cov << " -- the grey outline underneath is the rest, and it gets no elements";
    console_.log(cov.str());
    // A block that straddles an interface is an element no analysis code can
    // integrate, so it is said plainly rather than left in a count: it is the
    // one thing about a multi-material layout that the picture will not show.
    if (d.straddlingBlocks > 0) {
        std::ostringstream bad;
        bad << "[Decomposition] " << d.straddlingBlocks
            << " block(s) sit in more than one material -- the layout does not follow the "
               "interfaces there, and the elements in them will straddle one";
        console_.log(bad.str());
    }
}

// ── UMBER: the quadrilateral mesh on the blocks, and the smoothing ──────────
//
// Stage 10's dialog, on Stage 10's settings object, because the number it asks
// for means the same thing here: a target edge length in the units of the
// model. Sharing it is what makes "mesh this model with UMBER and with
// MERIDIAN at 0.05" a comparison of two layouts rather than of two targets.
bool CrossGenWidget::promptUMBERMesh() {
    if (!umberDecomp_.has_value() || umberDecomp_->blocks.empty()) return false;
    return promptBlockQuadMesh(*umberDecomp_, "UMBER — quadrilateral mesh on the blocks");
}

// The dialog both UMBER and ZIPLINE open on their Mesh phase: the target edge
// length and the chord bounds of mesh/BlockQuadMesh, on the one settings object
// the pipelines' Stage 10 and ATLAS's mesh use too.
bool CrossGenWidget::promptBlockQuadMesh(const BlockDecomposition &decomp, const char *title) {
    if (decomp.blocks.empty() || !mesh_) return false;

    Point lo{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    Point hi{-lo[0], -lo[1]};
    for (const Point &p : mesh_->vertices) {
        lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
        hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
    }
    const double diag = std::hypot(hi[0] - lo[0], hi[1] - lo[1]);
    const double extent = (diag > 0.0) ? diag : 1.0;

    QDialog dlg(this);
    dlg.setWindowTitle(QString::fromUtf8(title));

    auto *targetBox = new QDoubleSpinBox(&dlg);
    targetBox->setRange(1e-4, 10.0);
    targetBox->setDecimals(4);
    targetBox->setSingleStep(0.01);
    targetBox->setValue(meshSettings_.target);
    targetBox->setToolTip(
        "Target length of a mesh edge, in the units of the model.\n"
        "The same number the MERIDIAN, TORSION and ATLAS mesh dialogs take,\n"
        "and shared with them, so the layouts can be meshed alike.");

    auto *derived = new QLabel(&dlg);
    derived->setTextFormat(Qt::PlainText);
    const int blocks = static_cast<int>(decomp.blocks.size());
    const int edges = static_cast<int>(decomp.edges.size());
    auto updateDerived = [targetBox, derived, extent, blocks, edges]() {
        const double h = targetBox->value();
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(4)
            << "diagonal of S                 = " << extent << "\n"
            << "edges across it   diag / h    = " << std::setprecision(1) << (extent / h) << "\n"
            << blocks << " block(s) and " << edges << " macro edge(s) to mesh";
        derived->setText(QString::fromStdString(oss.str()));
    };
    QObject::connect(targetBox, &QDoubleSpinBox::valueChanged, &dlg, updateDerived);
    updateDerived();

    auto *minBox = new QSpinBox(&dlg);
    minBox->setRange(1, 64);
    minBox->setValue(meshSettings_.minEdges);
    minBox->setToolTip("Fewest edges any chord may be given.");

    auto *maxBox = new QSpinBox(&dlg);
    maxBox->setRange(0, 4096);
    maxBox->setValue(meshSettings_.maxEdges);
    maxBox->setSpecialValueText("none");
    maxBox->setToolTip("Most edges any chord may be given. 0 for no ceiling.");

    auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dlg);
    buttons->button(QDialogButtonBox::Ok)->setText("Mesh");
    QObject::connect(buttons, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
    QObject::connect(buttons, &QDialogButtonBox::rejected, &dlg, &QDialog::reject);

    auto *form = new QFormLayout(&dlg);
    form->addRow("target edge length", targetBox);
    form->addRow("implied sizing", derived);
    form->addRow("fewest edges per chord", minBox);
    form->addRow("most edges per chord", maxBox);
    form->addRow(new QLabel("Transfinite interpolation of the four sides per block.\n"
                            "Only a block that folds is smoothed, by Stage 10's\n"
                            "Winslow pass.", &dlg));
    form->addRow(buttons);

    if (dlg.exec() != QDialog::Accepted) return false;

    meshSettings_.target   = targetBox->value();
    meshSettings_.minEdges = minBox->value();
    meshSettings_.maxEdges = maxBox->value();
    return true;
}

// BlockQuadMesh at the dialog's settings, reported line for line as ATLAS's
// mesh and Stage 10's are.
void CrossGenWidget::runUMBERMesh() {
    if (!umberDecomp_.has_value() || umberDecomp_->blocks.empty()) {
        blockPipeline("the mesh not built",
                      "the tracing left no four-sided block to mesh");
        return;
    }
    umberMeshAttempted_ = true;
    umberMesh_.reset();
    // TMOP stood on the mesh that is about to be replaced.
    smoothMesh_.reset();
    tmopAttempted_ = false;
    pipelineBlocked_.clear();

    BlockQuadMesh::Options mo;
    mo.targetEdgeLength = meshSettings_.target;
    mo.minIntervals     = meshSettings_.minEdges;
    mo.maxIntervals     = meshSettings_.maxEdges;
    // Stage 10's selective Winslow pass at the pipelines' setting, so the mesh
    // TMOP is handed and the quality reported here are the same kind of thing
    // for all three methods.
    mo.smoothingPasses    = TORSION::Options().quadSmoothingPasses;
    mo.smoothingThreshold = TORSION::Options().quadSmoothingThreshold;

    auto t0 = Clock::now();
    try {
        umberMesh_.emplace(*umberDecomp_, mo);
    } catch (const std::exception &e) {
        umberMesh_.reset();
        blockPipeline("the mesh failed", e.what());
        return;
    }
    auto t1 = Clock::now();

    reportBlockQuadMesh(*umberMesh_,
                        std::chrono::duration<double, std::milli>(t1 - t0).count());
}

// The console report of a BlockQuadMesh, line for line as ATLAS's mesh and
// Stage 10's are, for whichever mode built it.
void CrossGenWidget::reportBlockQuadMesh(const BlockQuadMesh &bqm, double ms) {
    const BlockQuadMesh::Report &r = bqm.getReport();
    {
        std::ostringstream oss;
        oss << "[Mesh] " << r.quads << " quad(s) on " << r.vertices << " vertices over "
            << r.blocks << " block(s), " << formatMs(ms);
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Mesh] " << r.chords << " chord(s) over " << r.edgesAssigned << " macro edge(s), "
            << r.minIntervals << " to " << r.maxIntervals << " edges each (mean " << std::fixed
            << std::setprecision(2) << r.meanIntervals << ")";
        if (r.clampedChords > 0) oss << ", " << r.clampedChords << " clamped by a bound";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Mesh] edges " << std::fixed << std::setprecision(4) << r.minEdge << " to "
            << r.maxEdge << " against a target of " << r.target << " (worst "
            << std::setprecision(2) << r.worstEdgeRatio << "x, rms log ratio "
            << std::setprecision(3) << r.edgeRatioRms << ")";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Mesh] scaled Jacobian " << std::fixed << std::setprecision(4)
            << r.minScaledJacobian << " worst, " << r.meanScaledJacobian << " mean";
        if (r.smoothedBlocks > 0) {
            oss << " (Winslow on " << r.smoothedBlocks << " folded block(s): " << r.invertedBefore
                << " fold(s), " << r.minScaledJacobianBefore << " worst before)";
        }
        if (r.invertedQuads > 0) oss << " -- " << r.invertedQuads << " element(s) fold";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Mesh] every boundary node on dS, the edges between them within " << std::fixed
            << std::setprecision(4) << r.boundaryDeviation << " of it (" << std::setprecision(3)
            << (r.target > 0.0 ? r.boundaryDeviation / r.target : 0.0) << " of the target)";
        if (r.materials > 1) {
            oss << "; " << r.materials << " materials meeting on " << r.interfaceEdges
                << " element edge(s), within " << std::setprecision(4) << r.interfaceDeviation
                << " of the interfaces";
        }
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Mesh] " << r.interiorEdges << " interior and " << r.boundaryEdges
            << " boundary edge(s), " << r.nonManifoldEdges << " used a third time, " << r.cracks
            << " crack(s): " << (r.conforming ? "conforming [PASS]" : "not conforming [FAIL]");
        console_.log(oss.str());
    }
    for (const std::string &m : r.messages) console_.log("[Mesh] " + m);
    {
        std::ostringstream oss;
        oss << "[Mesh] the mesh " << (r.valid ? "validates [PASS]" : "does not validate [FAIL]");
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }
    console_.log("[Mesh] transfinite grid: grey is a mesh edge, blue a block wall, red a folded "
                 "element. Press 'c' to smooth it, 'e' to mesh again at another target");
}

// ── MERIDIAN: Shepherd, Gu and Hughes (2022), Stages 0b-3 ────────────────────

// Stage 0b, the material interface network. Everything it reads is already on
// the mesh -- which triangle carries which material tag -- so it runs before the
// field and before the cones, and on a single-material model it is one pass over
// the edges that finds nothing.
//
// What it produces is not a picture but a set of demands on the stages after it:
// a cone index at every node of the network for Stage 1, a sector constraint at
// each of them for Stage 6, an emitter for Stage 7 and an arc for Stage 8. The
// console is where those demands are visible before their consequences are, so
// the report is printed in full rather than summarised.
void CrossGenWidget::runMERIDIANInterfaces() {
    interfacesAttempted_ = true;
    interfaceMessagesSeen_ = 0;
    if (!mesh_) return;

    auto t0 = Clock::now();
    try {
        interfaces_.emplace(mesh_);
    } catch (const std::exception &e) {
        interfaces_.reset();
        console_.log(std::string("[Interfaces] FAILED: ") + e.what());
        std::cerr << "[Viewer] Interfaces failed: " << e.what() << "\n";
        return;
    }
    auto t1 = Clock::now();

    const Interfaces::Report &fr = interfaces_->getReport();
    if (!interfaces_->multiMaterial()) {
        // Said once, and not as a warning: the single-material pipeline is
        // exactly what it was, and this is the line that says why nothing else
        // about interfaces appears below.
        console_.log("[Interfaces] Stage 0b: one material, so there are no interfaces "
                     "and the pipeline is the single-material one throughout");
        return;
    }

    {
        std::ostringstream oss;
        oss << "[Interfaces] Stage 0b: " << fr.materials << " material(s), "
            << fr.interfaceEdges << " interface edge(s) in " << fr.branches
            << " branch(es)";
        if (fr.closedLoops > 0) oss << " (" << fr.closedLoops << " closed loop(s))";
        oss << ", " << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Interfaces] " << fr.nodes << " node(s): " << fr.junctions << " junction, "
            << fr.landings << " on dS, " << fr.kinks << " kink, " << fr.loopSplits
            << " loop split";
        if (fr.dangling > 0) oss << ", " << fr.dangling << " dangling";
        console_.log(oss.str());
    }
    {
        // The residual is the whole of what "ill-posed" means here, so it is
        // reported as an angle and not as a flag: a junction 28 degrees off the
        // quarter turns is a property of the domain, and the corpus has several
        // on purpose.
        std::ostringstream oss;
        oss << "[Interfaces] sectors: " << fr.illPosedNodes << " of " << fr.nodes
            << " node(s) more than 10 degrees from whole quarter turns, worst "
            << std::fixed << std::setprecision(1)
            << (fr.worstSectorResidual * 180.0 / M_PI) << " deg";
        if (fr.worstSectorNode >= 0 &&
            fr.worstSectorNode < static_cast<int>(interfaces_->nodes().size())) {
            oss << " at vertex " << interfaces_->nodes()[fr.worstSectorNode].vertex;
        }
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Interfaces] " << fr.prescribedCones
            << " cone(s) prescribed from the geometry, summing to " << std::showpos
            << fr.prescribedIndexSum << std::noshowpos
            << "; Stage 1 will be told them rather than reading them off the field";
        console_.log(oss.str());
    }
    if (fr.dangling > 0 || fr.nonManifoldVertices > 0) {
        std::ostringstream oss;
        oss << "[Interfaces] " << fr.dangling << " dangling branch(es) and "
            << fr.nonManifoldVertices
            << " vertex/-ices where the network is neither a path nor a node -- these are "
            << "tag errors in the input, not something a later stage can repair";
        console_.log(oss.str());
    }
    for (; interfaceMessagesSeen_ < fr.messages.size(); ++interfaceMessagesSeen_)
        console_.log("[Interfaces] " + fr.messages[interfaceMessagesSeen_]);

    console_.log("[Interfaces] white = an interface branch; the disks are its nodes, "
                 "coloured by kind (see the key). Press 'i' to hide them, 'm' to fill "
                 "the triangles by material");
}

// Stage 1, Sec. 3.1. Interior indices are the DualMBO field's winding numbers,
// boundary ones come from the field's rotation across each boundary star
// against the interior angle there; then Eq. (4) is checked, and repaired from
// the boundary cones if the rounding on them cost it.
//
// This is the gate for everything after. sum(Kbar) = 2 pi chi is the
// solvability condition of the Newton system in Stage 3 -- the Laplacian's
// kernel is the constants, so the residual has to be orthogonal to them -- and
// an inadmissible set does not converge slowly, it has no solution at all.
void CrossGenWidget::runMERIDIANCones() {
    if (!dualMBOField_.has_value()) return;   // called before the field: not an attempt
    conesAttempted_ = true;

    // Cone indices are read off a *converged* field: an interior winding number
    // is an exact integer only once the field has stopped moving, and a
    // boundary one is measured against a field that is supposed to be aligned
    // with the boundary. Advancing past the stepping phase early therefore
    // finishes the solve here rather than reading a field still in motion.
    if (!dualMBOConverged_) {
        auto t0 = Clock::now();
        // The pipeline's own loop rather than DualMBO::runMBO(), which stops at
        // a tolerance of its own (1e-7): what this finishes has to be the field
        // TORSION::runField() would have handed Stage 1, or the cone indices
        // read off it below are read off a different field.
        const double ntris = static_cast<double>(mesh_->triangles.size());
        const int cap = TORSION::Options().dualMBOMaxSteps;
        for (int i = dualMBOStepCount_; i < cap; ++i) {
            dualMBOField_->step();
            ++dualMBOStepCount_;
            if (dualMBOField_->error < 2.0 * ntris * DUALMBO_TOL) break;
        }
        dualMBOConverged_ = true;
        auto t1 = Clock::now();
        std::ostringstream oss;
        oss << "[MERIDIAN] finished the DualMBO solve first, error " << std::scientific
            << std::setprecision(3) << dualMBOField_->error << ", "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }

    auto t0 = Clock::now();
    try {
        cones_.emplace(*dualMBOField_);
    } catch (const std::exception &e) {
        cones_.reset();
        console_.log(std::string("[Cones] FAILED: ") + e.what());
        std::cerr << "[Viewer] ConeSingularities failed: " << e.what() << "\n";
        return;
    }
    // The snapshot Sec. 5.1's audit is run against, taken here because this is
    // the last moment at which the indices are the ones the *field* read.
    // Unused by MERIDIAN and cheap, so it is taken in both modes rather than
    // conditionally: what it costs is one vector and what it buys is that the
    // cone stage stays one piece of code.
    fieldIndex_ = cones_->getIndices();

    // The interface network's own indices, before Eq. (4) is checked: they are
    // fixed by the geometry and rebalance() may not move them, so the check and
    // the rebalance below have to see them already in place.
    if (interfaces_.has_value() && interfaces_->multiMaterial()) {
        cones_->prescribe(interfaces_->prescription());
        std::ostringstream oss;
        oss << "[Cones] " << interfaces_->getReport().prescribedCones
            << " index/-ices taken from the interface network, "
            << cones_->prescriptionShift()
            << " unit(s) away from what the field read there -- the field is "
               "boundary-aligned and has no boundary condition on an interface, so at a "
               "junction or a landing it reports whatever propagated in";
        console_.log(oss.str());
    }

    // The pairs the aligned field puts on a curved interface, which no region
    // asked for: a smoothest field hard-aligned to a strongly curved feature
    // carries a +1/-1 dipole on the concave side. Before Eq. (4) is checked
    // (a pair sums to zero and cannot change it) and before Stage 2, which
    // would otherwise cut to each of them. A pair spanning two regions -- the
    // one a region boundary genuinely needs -- is left alone; only a pair that
    // landed inside a single region's own count is cancelled.
    if (interfaces_.has_value() && interfaces_->multiMaterial()) {
        const int nV = static_cast<int>(mesh_->vertices.size());
        std::vector<int> region(nV, -1);
        for (int v = 0; v < nV; ++v) region[v] = interfaces_->regionAt(v);
        const int cancelled = cones_->cancelDipoles(region);
        if (cancelled > 0) {
            std::ostringstream oss;
            oss << "[Cones] cancelled " << cancelled << " +1/-1 cone pair(s) the field put "
                   "inside a single material region";
            console_.log(oss.str());
        }
    }

    auto gb = cones_->gaussBonnet();
    if (!gb.admissible) {
        const int moved = cones_->rebalance();
        gb = cones_->gaussBonnet();
        if (moved > 0) {
            std::ostringstream oss;
            oss << "[Cones] rebalanced " << moved << " index unit(s) onto boundary cones at a "
                << "cost of " << std::fixed << std::setprecision(2) << gb.rebalanceCost
                << " quarter turns";
            console_.log(oss.str());
        } else if (moved < 0) {
            console_.log("[Cones] could not restore Eq. (4): no boundary cone left within "
                         "the allowed index range");
        }
    }
    auto t1 = Clock::now();

    {
        std::ostringstream oss;
        oss << "[Cones] " << cones_->interiorCones().size() << " interior, "
            << cones_->boundaryCones().size() << " boundary; V-E+F = "
            << gb.eulerCharacteristic;
        if (gb.isolatedVertices > 0) oss << " (" << gb.isolatedVertices << " isolated v. excluded)";
        oss << ", " << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Cones] Eq. (4): sum I(v) = " << gb.indexSum << ", 4 chi = " << gb.indexTarget
            << " -- " << (gb.admissible ? "admissible [PASS]" : "NOT admissible [FAIL]");
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }
    {
        std::ostringstream oss;
        oss << "[Cones] Eq. (7): sum K_v - 2 pi chi = " << std::scientific << std::setprecision(2)
            << (gb.curvatureSum - gb.curvatureTarget)
            << (gb.metricConsistent ? " [PASS]" : " [FAIL]");
        console_.log(oss.str());
    }
    if (!gb.admissible) {
        console_.log("[Cones] the flat cone metric asked for does not exist; Ricci flow "
                     "will refuse to start");
    }

    // The second half of Stage 0b: each material region's own Gauss-Bonnet
    // count. It runs here because the cones on dS are half of that count and
    // they are not known until Stage 1 has placed them. Where a region is short
    // a quarter turn, this is what puts one somewhere it can go -- and without
    // it E3 and E6 spend the whole continuation asking for a map that does not
    // exist.
    if (interfaces_.has_value() && interfaces_->multiMaterial()) {
        interfaces_->balance(cones_->getIndices());
        const Interfaces::Report &fr = interfaces_->getReport();
        std::ostringstream oss;
        oss << "[Interfaces] regions: " << fr.regionsBalanced << " of " << fr.regions
            << " satisfy sum I + sum (2 - q) = 4 chi";
        if (fr.quartersMoved > 0 || fr.cornersInserted > 0) {
            oss << ", after moving " << fr.quartersMoved << " quarter turn(s) across a node "
                << "and inserting " << fr.cornersInserted << " corner(s) on a smooth "
                << "interface";
        }
        if (fr.regionsBalanced < fr.regions)
            oss << "; the worst is " << fr.worstRegionDeficit << " quarter(s) out";
        console_.log(oss.str());
        if (fr.cornersInserted > 0) {
            // The counts the picture now shows: inserting a corner splits the
            // branch it went on, so the network drawn from here on is not the
            // one reported above.
            std::ostringstream oss2;
            oss2 << "[Interfaces] the network is now " << fr.branches << " branch(es) and "
                 << fr.nodes << " node(s)";
            console_.log(oss2.str());
        }
        for (; interfaceMessagesSeen_ < fr.messages.size(); ++interfaceMessagesSeen_)
            console_.log("[Interfaces] " + fr.messages[interfaceMessagesSeen_]);
        if (fr.regionsBalanced < fr.regions) {
            console_.log("[Interfaces] a region that cannot be a whole number of "
                         "quadrilaterals is a contradiction and not a slack constraint; "
                         "the remedy is upstream, in where Stage 1 put the boundary cones");
        }
    }
}

// Stage 2, Sec. 3.2.2. HarmonicCut opens the voids and the cone arcs drag every
// interior cone out to the boundary, leaving S - G a disk with P in G union dS.
void CrossGenWidget::runMERIDIANCut() {
    if (!cones_.has_value()) return;       // called before Stage 1: not an attempt
    cutAttempted_ = true;

    auto t0 = Clock::now();
    try {
        // The interface network goes in with the cones: on a multi-material
        // mesh the cut has to be routed around it, or the seam lands on a
        // curve the layout is meant to keep. See ConeCut's class comment.
        coneCut_.emplace(mesh_, *cones_, ConeCut::Options(),
                         interfaces_.has_value() ? &*interfaces_ : nullptr);
    } catch (const std::exception &e) {
        coneCut_.reset();
        // On the banner as well as in the console: everything from Stage 4 on
        // is built on the cut, so this failure is what the four blank phases
        // after it are, and the console will have scrolled past it by then.
        blockPipeline("Stage 2 failed", e.what());
        return;
    }
    auto t1 = Clock::now();

    const ConeCut::Report &r = coneCut_->getReport();
    {
        std::ostringstream oss;
        oss << "[Cut] " << r.harmonicCuts << "/" << r.voids << " void arc(s), "
            << r.conesRouted << "/" << r.interiorCones << " cone arc(s), "
            << coneCut_->getCutEdges().size() << " edges, "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Cut] Omega: chi = " << r.eulerCharacteristic << ", "
            << r.boundaryComponents << " boundary component(s) -- "
            << (r.isDisk ? "disk [PASS]" : "not a disk [FAIL]") << "; P in G u dS "
            << (r.allConesOnBoundary ? "[PASS]" : "[FAIL]");
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }
    if (r.conesSplitByCut > 0) {
        std::ostringstream oss;
        oss << "[Cut] " << r.conesSplitByCut
            << " cone(s) the graph runs across rather than stopping at -- legal, but "
               "Sec. 3.2.2 prefers otherwise";
        console_.log(oss.str());
    }
    console_.log("[Cut] magenta = void arcs (one per hole), amber = cone arcs "
                 "(one per interior cone)");
}

// Stage 3, Sec. 3.2.1. Blocking: a few Newton steps and a flip pass or two,
// milliseconds on these models but a solve all the same, so it is announced a
// frame ahead like the UMBER ones. Everything the last two phases draw is
// derived here rather than per frame -- both are a walk over every edge or
// every cone star.
void CrossGenWidget::runRicciFlow() {
    if (!cones_.has_value()) return;       // called before Stage 1: not an attempt
    ricciAttempted_ = true;

    auto t0 = Clock::now();
    bool ok = false;
    try {
        ricci_.emplace(mesh_, *cones_);
        ok = ricci_->solve();
    } catch (const std::exception &e) {
        ricci_.reset();
        console_.log(std::string("[Ricci] FAILED: ") + e.what());
        std::cerr << "[Viewer] RicciFlow failed: " << e.what() << "\n";
        return;
    }
    auto t1 = Clock::now();

    const RicciFlow::Report &r = ricci_->getReport();
    {
        std::ostringstream oss;
        oss << "[Ricci] ||K - Kbar||_inf " << std::scientific << std::setprecision(2)
            << r.initialError << " -> " << r.finalError << " in " << r.newtonIterations
            << " Newton step(s), " << r.flips << " flip(s), "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }
    {
        std::ostringstream oss;
        oss << "[Ricci] flat cone metric " << (ok ? "reached [PASS]" : "NOT reached [FAIL]")
            << "; cos(phi) in [" << std::fixed << std::setprecision(2)
            << r.minInversiveDistance << ", " << r.maxInversiveDistance << "]";
        if (r.minInversiveDistance >= 1.0) oss << " (all >= 1: the generalised flow of [87])";
        console_.log(oss.str());
    }
    for (const std::string &m : r.messages) console_.log("[Ricci] " + m);
    if (!ok) return;

    // u for the ramp, with its mean removed: only differences in u mean
    // anything, since u -> u + c is a global rescaling of the metric.
    const std::vector<double> &u = ricci_->getU();
    const std::vector<char> &active = cones_->getActiveVertices();
    double sum = 0.0;
    int n = 0;
    for (size_t v = 0; v < u.size(); ++v) if (active[v]) { sum += u[v]; ++n; }
    const double mean = (n > 0) ? sum / n : 0.0;

    ricciU_.resize(static_cast<Eigen::Index>(u.size()));
    ricciUAbsMax_ = 1e-12;
    for (size_t v = 0; v < u.size(); ++v) {
        ricciU_[static_cast<Eigen::Index>(v)] = active[v] ? u[v] - mean : 0.0;
        ricciUAbsMax_ = std::max(ricciUAbsMax_, std::fabs(ricciU_[static_cast<Eigen::Index>(v)]));
    }

    flatMetric_ = viewer::buildFlatMetric(*ricci_);
    coneFans_ = viewer::buildConeFans(*ricci_, *cones_);

    // Fit the right panel to the gallery now, while it is being built, the way
    // the other split-screen modes fit theirs as their parameterization lands.
    if (!coneFans_.empty()) {
        viewer::computeConeFanBounds(coneFans_, uvView_.cx, uvView_.cy,
                                     uvView_.baseW, uvView_.baseH);
        uvView_.zoom = 1.0;
        uvView_.fbw  = view_.fbw;
        uvView_.fbh  = view_.fbh;
    }

    {
        std::ostringstream oss;
        oss << "[Ricci] conformal factor gamma = e^u spans " << std::fixed
            << std::setprecision(2) << std::exp(2.0 * ricciUAbsMax_)
            << "x; edge lengths " << flatMetric_.minRatio << "x to "
            << flatMetric_.maxRatio << "x about the mean";
        console_.log(oss.str());
    }
    if (!flatMetric_.replaced.empty()) {
        std::ostringstream oss;
        oss << "[Ricci] " << flatMetric_.replaced.size()
            << " input edge(s) replaced by flipping (grey underneath, green on top)";
        console_.log(oss.str());
    }
}

// ── TORSION Stage 3F, Sec. 5: the combed field ───────────────────────────────
//
// The phase Pipeline A runs the flow in, doing the same job by other means:
// turning the field into the structure the layout will have. DualMBO stores
// u_k[t] = exp(4 i theta_t) in the global frame, so on a planar model parallel
// transport is identically zero, the matching across an interior edge is the
// integer p_fg the two representatives differ by, and combing is a BFS over the
// faces of Omega carrying an integer a_f.
//
// Cheap -- one BFS and one pass over the vertices -- so unlike the Ricci solve
// it needs no announcement. What it produces that nothing downstream can do
// without is the branch: the quarter turn across each arc of G is read off the
// frames either side of it, which is what makes C2 hold by construction instead
// of by a fit and a rounding.
void CrossGenWidget::runTORSIONFrames() {
    if (!dualMBOField_.has_value() || !coneCut_.has_value() || !cones_.has_value()) return;
    framesAttempted_ = true;

    // Sec. 4, and C4. The frame the plan first wrote down was (1/h) R(-theta),
    // which asks the map to be an isometry of the model everywhere and so is
    // not integrable at a cone whatever the field does; the conformal factor of
    // the flat cone metric is what it is missing, and the same u is the
    // reference metric Stage 6 needs. One Newton solve answers both.
    coneMetric_.reset();
    {
        ConeMetric::Options cmopts;
        try {
            coneMetric_.emplace(*mesh_, *cones_, cmopts);
        } catch (const std::exception &e) {
            console_.log(std::string("[Frames] cone metric FAILED: ") + e.what());
        }
        if (coneMetric_.has_value()) {
            const ConeMetric::Report &cmr = coneMetric_->getReport();
            std::ostringstream oss;
            oss << "[Frames] Sec. 4 cone metric: ||K - Kbar||_inf " << std::scientific
                << std::setprecision(2) << cmr.initialError << " -> " << cmr.finalError
                << " rad in " << cmr.newtonIterations << " Newton step(s), exp(u) in ["
                << std::fixed << std::setprecision(3) << cmr.minScale << ", " << cmr.maxScale
                << "], " << cmr.nonRealisable << " face(s) not realisable "
                << (cmr.solved ? "[PASS]" : "[FAIL]");
            console_.log(oss.str());
            for (const std::string &m : cmr.messages) console_.log("[Frames] " + m);
            if (!cmr.solved) coneMetric_.reset();
        }
    }

    FieldFrames::Options fopts;
    fopts.referenceIndex = fieldIndex_;
    if (coneMetric_.has_value()) {
        const TORSION::Options defaults;
        fopts.sizing = coneMetric_->sizingField(defaults.targetEdge);
        fopts.metricLengths = coneMetric_->edgeLengths();
        if (defaults.targetEdge > 0.0) {
            for (double &l : fopts.metricLengths) l /= defaults.targetEdge;
        }
    }

    auto t0 = Clock::now();
    try {
        // The interfaces go in for Sec. 6.4's alignment: an interface branch is
        // a curve of the layout in exactly the way a chain of dS is, and the
        // field is pinned tangent to it in the same way.
        frames_.emplace(*dualMBOField_, *coneCut_, *cones_, fopts,
                        interfaces_.has_value() ? &*interfaces_ : nullptr);
    } catch (const std::exception &e) {
        frames_.reset();
        console_.log(std::string("[Frames] FAILED: ") + e.what());
        std::cerr << "[Viewer] FieldFrames failed: " << e.what() << "\n";
        return;
    }
    auto t1 = Clock::now();

    const FieldFrames::Report &fr = frames_->getReport();
    {
        std::ostringstream oss;
        oss << "[Frames] combed " << fr.combedFaces << "/" << fr.faces
            << " face(s) of Omega from face " << fr.seedFace << " ("
            << fr.unreachedFaces << " unreached), " << fr.combingDefects
            << " loop defect(s) " << (fr.combingDefects == 0 ? "[PASS]" : "[FAIL]") << ", "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }
    {
        // The roughness of the field, which is what predicts where Stage 4F
        // will invert triangles: delta is the field's own variation across a
        // dual edge once the quarter turn is taken out, and it is bounded by
        // pi/4 by construction, so a value near it is a face where the two
        // representatives were nearly a coin toss apart.
        std::ostringstream oss;
        oss << "[Frames] worst |delta| across a dual edge " << std::fixed
            << std::setprecision(3) << fr.maxFrameJump << " rad of the pi/4 = 0.785 the "
            << "matching leaves; this is where Stage 4F's flips will be";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Frames] Sec. 5.1: I from the matchings agrees with the field at all but "
            << fr.indexMismatches << " interior vertex/-ices "
            << (fr.indexMismatches == 0 ? "[PASS]" : "[FAIL]") << "; sum I = " << fr.indexSum
            << " against 4 chi = " << fr.indexTarget << " "
            << (fr.admissible ? "[PASS]" : "[FAIL]");
        if (fr.highIndexCones > 0) oss << ", " << fr.highIndexCones << " with |I| > 2";
        console_.log(oss.str());
    }
    if (fr.boundaryTurnMismatches > 0) {
        // The one place a field-integrated layout is asked for two different
        // things at one point, and it is not a bug in anything here: Pipeline A
        // gives Stage 6 a geodesic boundary because Stage 1 prescribes zero
        // curvature at every non-cone boundary vertex and the flow drives it
        // there. The field has no such stage -- it follows a curving boundary
        // and quantises it into a staircase -- so every step of that staircase
        // away from a cone is a corner the arrangement will find and Q2 will
        // report. It is why the two routes' arrangements differ.
        std::ostringstream oss;
        oss << "[Frames] the frame turns at " << fr.boundaryTurns << " vertex/-ices of dS, "
            << fr.boundaryTurnMismatches << " of them disagreeing with the cone set (worst "
            << "at vertex " << fr.worstBoundaryTurnVertex << ", rounding residual "
            << std::fixed << std::setprecision(3) << fr.maxBoundaryTurnResidual
            << " rad) -- these are the corners Q2 will find at vertices Stage 1 called "
               "regular";
        console_.log(oss.str());
    }
    if (fr.leftHandedFrames > 0 || fr.sizingFallbacks > 0) {
        std::ostringstream oss;
        oss << "[Frames] " << fr.leftHandedFrames << " frame(s) with det J* <= 0 and "
            << fr.sizingFallbacks << " face(s) whose h_t was not positive -- both make "
               "Stage 6's flip cap meaningless where they occur";
        console_.log(oss.str());
    }
    for (const std::string &m : fr.messages) console_.log("[Frames] " + m);

    console_.log(std::string("[Frames] the branch ") +
                 (fr.valid ? "is single valued over Omega [PASS]"
                           : "is not single valued over Omega [FAIL]") +
                 "; colour = a_f, and a step in a_f anywhere but across an arc of G is a "
                 "leak in the comb");
}

// ── TORSION Stages 4F, 4R and 4, Secs. 6, 7.2 and 6.2: psi_0 ─────────────────
//
// The phase Pipeline A shows the flat metric in, and the reason it is one phase
// there is the reason it is one here read backwards. On that side the metric is
// a set of edge lengths and needs unfolding before there is anything to look
// at; on this side the solve produces a *map* immediately, and what the phase
// is judged on is not what the map looks like but whether it inverted anything.
// The red faces on the left panel's right half are the entire cost of
// substituting an integration for a flow, which is why 'p' keeps the
// least-squares map reachable after Stage 4R has replaced it.
//
// Blocking: one sparse saddle solve, and when it inverted something, a whole
// continuation after it.
void CrossGenWidget::runTORSIONIntegration() {
    if (!frames_.has_value() || !coneCut_.has_value() || !cones_.has_value()) return;
    integrationAttempted_ = true;

    const FieldFrames::Report &ffr = frames_->getReport();
    if (ffr.unreachedFaces > 0 || ffr.combingDefects > 0) {
        console_.log("[psi_0] stopping before the integration: the branch of the field over "
                     "Omega is not single valued, so the seam transitions it would be "
                     "constrained by are not defined");
        return;
    }

    // The scaffold exists for one reason: the constraint rows need the arcs of
    // G, their (e+, e-) pairing and their quarter turns, and all three are
    // Immersion::buildArcs()' work. It is handed Omega's own coordinates as a
    // placeholder map, which nothing read from it depends on.
    try {
        scaffold_.emplace(*coneCut_, *cones_, coneCut_->getCutMesh().vertices,
                          frames_->fieldEdgeLengths(), frames_->combedAngle());
    } catch (const std::exception &e) {
        scaffold_.reset();
        console_.log(std::string("[psi_0] could not build the seam pairing: ") + e.what());
        return;
    }

    const TORSION::Options defaults;

    // Stages 4F and 4R are TORSION::buildPsi0(), the function run() calls:
    // Sec. 6.4's attempts at the alignment, Sec. 7.2a-c's ladder, the Tutte pass
    // behind it and Sec. 7.2b's pull back onto the alignment, in the pipeline's
    // order under its options. This phase used to carry a copy of all of it, and
    // the copy did not keep up with the pipeline.
    auto t0 = Clock::now();
    TORSION::Status st;
    TORSION::Psi0 built;
    bool built0 = false;
    usedAxis_.clear();
    try {
        built0 = TORSION::buildPsi0(
            TORSION::MapStage{*mesh_, *coneCut_, *cones_, *frames_, *scaffold_,
                              interfaces_.has_value() ? &*interfaces_ : nullptr,
                              coneMetric_.has_value() ? &*coneMetric_ : nullptr},
            defaults, st, built);
    } catch (const std::exception &e) {
        integration_.reset();
        console_.log(std::string("[psi_0] FAILED: ") + e.what());
        std::cerr << "[Viewer] TORSION::buildPsi0 failed: " << e.what() << "\n";
        return;
    }
    auto t1 = Clock::now();
    usedAxis_ = built.usedAxis;
    if (built.integration) integration_.emplace(std::move(*built.integration));
    else integration_.reset();
    if (built.tutte) tutte_.emplace(std::move(*built.tutte));
    else tutte_.reset();
    if (!integration_.has_value()) {
        console_.log("[psi_0] the integration could not be built");
        return;
    }

    const FieldIntegration::Report &ir = integration_->getReport();
    {
        std::ostringstream oss;
        oss << "[psi_0] Stages 4F and 4R: " << ir.constraintRows << " constraint row(s) over "
            << ir.seamPairs << " seam pair(s), " << (ir.solvedWithLDLT ? "LDL^T" : "LU fallback")
            << ", seam " << std::scientific << std::setprecision(2) << ir.maxSeamResidual
            << " of the extent, "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }
    if (!ir.solved) {
        console_.log("[psi_0] the field could not be integrated on Omega");
        return;
    }
    {
        // The non-integrability, which is the thing being measured here. A
        // discrete cross field generically has curl X != 0, so the closest map
        // in L^2 inverts triangles; where the per-face fit residual is large is
        // where they are.
        std::ostringstream oss;
        oss << "[psi_0] fit residual mean " << std::fixed << std::setprecision(4)
            << ir.meanFitResidual << ", worst " << ir.maxFitResidual << " at face "
            << ir.worstFitFace << " (0 = the field realised exactly, 1 = a gradient the "
               "size of the target pointing the wrong way)";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[psi_0] the solve kept inverted " << ir.flippedFaces << " face(s), "
            << std::fixed << std::setprecision(2)
            << (ir.totalArea > 0.0 ? 100.0 * ir.flippedArea / ir.totalArea : 0.0)
            << "% of the area, " << ir.flipsAdjacentToCone << " of them in the one ring of "
            << "a cone; nearest flip " << std::scientific << std::setprecision(2)
            << ir.nearestFlipToCone << " of the diagonal from one";
        console_.log(oss.str());
    }
    for (const std::string &m : st.messages) {
        if (m.rfind("Stage 4F: ", 0) == 0) console_.log("[psi_0] " + m.substr(10));
        else if (m.rfind("Stage 4R: ", 0) == 0) console_.log("[Untangle] " + m.substr(10));
    }

    integratedMap_ = std::move(built.integratedMap);
    showIntegrated_ = false;
    if (!built0 || built.psi0.empty()) {
        console_.log("[psi_0] Stage 4R left no map to hand on");
        return;
    }
    std::vector<Point> psi0 = std::move(built.psi0);

    // ── Stage 4: the immersion both routes end at ────────────────────────────
    try {
        immersion_.emplace(*coneCut_, *cones_, psi0, frames_->fieldEdgeLengths(),
                           frames_->combedAngle());
    } catch (const std::exception &e) {
        immersion_.reset();
        console_.log(std::string("[psi_0] Stage 4 could not accept psi_0: ") + e.what());
        return;
    }

    const Immersion::Report &imr = immersion_->getReport();
    {
        std::ostringstream oss;
        oss << "[psi_0] " << imr.arcs << " arc(s), " << imr.seamEdgePairs
            << " seam pair(s), " << imr.frameKConflicts
            << " arc(s) whose frames disagree about their own quarter turn "
            << (imr.frameKConflicts == 0 ? "[PASS]" : "[FAIL]");
        console_.log(oss.str());
    }
    {
        // Q4 read as a check rather than as a rounding, which is what the field
        // route buys: k came off the matchings, so this number is how far the
        // geometry disagrees with them instead of how much it had to be snapped.
        std::ostringstream oss;
        oss << "[psi_0] Q4 the geometry is " << std::scientific << std::setprecision(2)
            << imr.maxSnapError << " rad from the quarter turns the matchings prescribe; "
            << "the field metric is realised to " << imr.maxMetricResidual
            << " relative, which is the non-integrability per edge and not an error";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[psi_0] Q1 " << (imr.flippedFaces == 0 ? "[PASS]" : "[FAIL]") << " ("
            << imr.flippedFaces << " flipped face(s))";
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }
    for (const std::string &m : imr.messages) console_.log("[psi_0] " + m);

    if (imr.flippedFaces > 0) {
        console_.log("[psi_0] Q1 does not hold, and Sec. 3.3's barrier can only preserve it, "
                     "never repair it, so the continuation will not be run. This is the cost "
                     "of the substitution: Sec. 7.3 is the list of what to try");
    }

    viewer::computeLayoutBounds(immersion_->getUV(), uvView_.cx, uvView_.cy,
                                uvView_.baseW, uvView_.baseH);
    uvView_.zoom = 1.0;
    uvView_.fbw  = view_.fbw;
    uvView_.fbh  = view_.fbh;

    console_.log("[psi_0] right panel = psi_0 as Stage 4 accepted it; press 'p' for the "
                 "Stage 4F least-squares map it came from, red where that inverted");
}

// Stages 4, 5 and 6 in one go, Secs. 3.2.2 and 3.3. They are one viewer phase
// because only the ends of the sequence have a picture and the two are the same
// picture moved: psi_R is where Stage 6 starts and Psi is where it stops, drawn
// on the same triangulation in the same plane, so 'p' swaps between them rather
// than 'c' advancing from one to the other.
//
// Blocking, and the longest wait in the mode -- the continuation is sixteen
// outer steps of up to sixty inner Newton-like solves each -- so it is
// announced a frame ahead like the Ricci solve.
//
// The three stages fail differently and are handled differently. A fold out of
// Stage 4 is terminal: E1's barrier keeps local injectivity and cannot restore
// it, so the continuation is not started at all and psi_R is left on screen as
// the thing to look at. A continuation that runs but does not converge is not
// terminal -- what it produced is still a map, still worth drawing, and its
// residuals are what say which property it failed.
// ── MERIDIAN connectivity dialog ─────────────────────────────────────────────

// Sec. 3.3's Gamma_topo, Sec. 3.4's snap tolerance and Sec. 3.3's repair, in
// one place, because they are one decision made three times.
//
// The decision is which pairs of cones are meant to be joined by an integral
// curve. E5 makes the curves that are told to close, close; a direction it was
// not told about is a geodesic of a flat cone metric, which does not end -- it
// winds until the step cap. So which connections were found is what decides
// whether a model comes out with a layout or with a handful of curves running
// forever, and it is a judgement about the model rather than a constant of the
// method. The paper agrees: Gamma_topo is an input there, its automatic
// generation is listed as future work, and the cones of the reference figure
// were placed by hand.
//
// Every field here has a default that the whole of data/meshes settles on, and
// the panel beside them says what each comes to in image units for *this*
// model, which is the only form in which any of them can be judged.
bool CrossGenWidget::promptMERIDIANConnectivity(const Immersion &imm,
                                                const std::vector<Point> &uv) {
    const Immersion::Report &ir = imm.getReport();
    const double extent = std::hypot(ir.uvMax[0] - ir.uvMin[0], ir.uvMax[1] - ir.uvMin[1]);
    const double spacing = SubdomainLabels::meanConeSpacing(imm, uv);

    // The closest two *distinct* cones come in the image. It is the number the
    // widened tolerances have to stay under, and on a model where Stage 1 left a
    // cluster it is very much smaller than the mean spacing.
    double closest = std::numeric_limits<double>::infinity();
    {
        const auto &children = imm.getConeChildren();
        for (size_t i = 0; i < children.size(); ++i) {
            for (size_t j = i + 1; j < children.size(); ++j) {
                for (int a : children[i]) {
                    for (int b : children[j]) closest = std::min(closest, normP(uv[a] - uv[b]));
                }
            }
        }
    }

    QDialog dlg(this);
    dlg.setWindowTitle("MERIDIAN — cone connectivity (Sec. 3.3, Gamma_topo)");

    // ── Stage 5: how the connections are found on psi_R ──────────────────────
    auto *nearMissBox = new QDoubleSpinBox(&dlg);
    nearMissBox->setRange(0.005, 0.60);
    nearMissBox->setDecimals(3);
    nearMissBox->setSingleStep(0.01);
    nearMissBox->setValue(meridianConn_.labels.nearMissTolerance);
    nearMissBox->setToolTip(
        "How near a separatrix of psi_R has to pass a cone for the pair to be\n"
        "seeded as a connectivity constraint, as a fraction of the mean spacing\n"
        "between cones.\n\n"
        "0.15 (recommended) is the largest value at which every model in\n"
        "data/meshes still converges. Too small and Remark 3.1's sliver problem\n"
        "survives, because the constraints that would have pinned the layout\n"
        "together were never seeded. Too large and paths are seeded between cones\n"
        "that were never meant to be joined: E5 then asks for a set of integral\n"
        "curves that no map has, and the continuation does not converge slowly --\n"
        "the barrier fights a penalty it cannot satisfy and det J is driven to 0.");

    auto *selfReturnBox = new QCheckBox("seed curves that return to their own cone", &dlg);
    selfReturnBox->setChecked(meridianConn_.labels.seedSelfReturns);
    selfReturnBox->setToolTip(
        "Q5 allows a separatrix to terminate at a \"possibly identical\"\n"
        "singularity, and the paper's Fig. 9 is exactly that curve: out of a cone\n"
        "along -du, across the cut, back along -dv to where it started.\n\n"
        "Recommended on. With it off, the one direction that most often needs\n"
        "quantising has nothing constraining it, and the curve winds instead of\n"
        "closing.");

    auto *allConnBox = new QCheckBox("keep every distinct curve between a pair of cones", &dlg);
    allConnBox->setChecked(meridianConn_.labels.seedAllConnections);
    allConnBox->setToolTip(
        "Two separatrices can join the same two cones in different directions,\n"
        "and they are two constraints, not one. With this off, only the nearer\n"
        "of the two is kept and the other is free to miss.\n\n"
        "The same curve found from both of its ends still collapses to one\n"
        "either way -- the deduplication is on the pair of ends, not the pair of\n"
        "cones.");

    // ── Stage 0b: what the material interfaces are held to ───────────────────
    //
    // Only shown on a multi-material model, where they are the two switches that
    // decide how much of the interface network the layout is actually made to
    // respect. Both are on, and both are here for the same reason the tolerances
    // above are: what each of them buys is a question about the model, and the
    // only way to settle it is to turn one off and look.
    const bool multiMat = interfaces_.has_value() && interfaces_->multiMaterial();

    auto *e6Box = new QCheckBox("hold the sector at each interface node (E6)", &dlg);
    e6Box->setChecked(meridianConn_.labels.interfaceCorners);
    e6Box->setToolTip(
        "Sec. 3.3's E3 already asks that each interface be a coordinate line, and\n"
        "it is a statement about one curve at a time. What it cannot say is what\n"
        "the map does where two of them meet: a junction whose sectors are one,\n"
        "one and two quarter turns is a condition on a *pair* of tangents, and\n"
        "nothing in E1-E5 asks for it.\n\n"
        "Off is the pipeline as it was before Stage 0b existed, and is what to\n"
        "compare against: the interfaces stay straight and the corners between\n"
        "them turn wherever the rest of the energy finds cheapest.");

    auto *propagateBox = new QCheckBox("label interface chains from the node quantisation", &dlg);
    propagateBox->setChecked(meridianConn_.labels.propagateInterfaceLabels);
    propagateBox->setToolTip(
        "Whether an interface chain is Gamma_u or Gamma_v is decided by walking\n"
        "the network from its nodes, quarter turn by quarter turn, rather than by\n"
        "reading the flux of psi_R along the chain.\n\n"
        "The flux is a measurement of a map that is not yet a layout, and at a\n"
        "triple point it can label all three branches u -- which asks E3 for a map\n"
        "holding u constant on three curves leaving one point in three directions.\n"
        "There is no such map, and the continuation spends every outer step\n"
        "failing to reach it.");

    // ── Stage 7: what counts as arriving ─────────────────────────────────────
    auto *snapBox = new QDoubleSpinBox(&dlg);
    snapBox->setRange(1e-9, 1e-2);
    snapBox->setDecimals(9);
    snapBox->setSingleStep(1e-6);
    snapBox->setValue(meridianConn_.trace.coneSnapTolerance);
    snapBox->setToolTip(
        "Sec. 3.4's snap tolerance, as a fraction of the image extent. The\n"
        "paper's table value is 1e-6 and it is the recommended one.\n\n"
        "It is tempting to raise this when curves miss, and it is the wrong\n"
        "knob. Stage 6 stops when its residuals are under 1e-6 of the extent, so\n"
        "a curve whose connection E5 was actually told about arrives well inside\n"
        "1e-6, while one it was not told about arrives wherever the geometry puts\n"
        "it. Raising the tolerance makes the second kind terminate at the price of\n"
        "a patch corner visibly off its cone. The repair below gives E5 the\n"
        "missing constraint instead, which is Sec. 3.3's own remedy.");

    auto *ringsBox = new QSpinBox(&dlg);
    ringsBox->setRange(0, 32);
    ringsBox->setValue(meridianConn_.trace.coneSnapRings);
    ringsBox->setToolTip(
        "How far from the curve a cone may be, counted in triangles, and still be\n"
        "a candidate to snap to.\n\n"
        "0 is the paper's rule exactly: only the cones at the corners of the\n"
        "triangle being crossed. The restriction is not an optimisation -- Psi\n"
        "overlaps itself, so a cone matched by image distance alone can be one the\n"
        "curve is nowhere near on S -- but it is tighter than it needs to be, and\n"
        "on a fine mesh a one-ring can be narrower than the tolerance. 2 is\n"
        "recommended; the distance test still gates every candidate.");

    auto *stepsBox = new QSpinBox(&dlg);
    stepsBox->setRange(1000, 2000000);
    stepsBox->setSingleStep(10000);
    stepsBox->setValue(meridianConn_.trace.maxSteps);
    stepsBox->setToolTip("Sec. 3.4's cap on triangle crossings. A curve that reaches it is\n"
                         "reported, not truncated silently.");

    auto *cycleBox = new QCheckBox("stop a curve once it is provably winding", &dlg);
    cycleBox->setChecked(meridianConn_.trace.detectCycles);
    cycleBox->setToolTip(
        "Away from its cones Psi is flat and a separatrix is one of its\n"
        "geodesics, so a curve Q5 does not close does not wander off -- it winds,\n"
        "and it winds until the step cap. Stopping it when it re-enters a triangle\n"
        "it has already crossed, in the same direction and on the same isoline to\n"
        "within the snap tolerance, reaches the same conclusion in a few hundred\n"
        "steps instead of fifty thousand, and reports it as a closed orbit rather\n"
        "than as \"ran out of budget\".\n\nRecommended on.");

    // ── Sec. 3.3's repair ────────────────────────────────────────────────────
    auto *repairBox = new QSpinBox(&dlg);
    repairBox->setRange(0, 20);
    repairBox->setValue(meridianConn_.repairPasses);
    repairBox->setToolTip(
        "Rounds of Sec. 3.3's remedy: take the separatrices of Psi that ended\n"
        "nowhere, and the ones that slipped past a cone on their way out through\n"
        "dS, give E5 the constraint each of them names, raise lambda_5, and\n"
        "re-run Stage 6 from the current phi -- never from psi_R.\n\n"
        "This is what the seeding on psi_R cannot do by itself: Stage 6 moves the\n"
        "map a long way, and a connection that is plain on the finished layout\n"
        "need never have been visible on the map it started from.\n\n"
        "0 disables the repair. A round that comes back worse is taken back\n"
        "whole, map and constraints together, so more rounds cannot make the\n"
        "result worse than fewer -- only slower.");

    auto *repairMaxBox = new QSpinBox(&dlg);
    repairMaxBox->setRange(0, 64);
    repairMaxBox->setValue(meridianConn_.repairMaxPerPass);
    repairMaxBox->setToolTip(
        "How many constraints one round may add, nearest miss first; 0 for all\n"
        "of them.\n\n"
        "4 is recommended. On a model whose cones Stage 1 left clustered, a dozen\n"
        "constraints switched on at once move the map somewhere none of them is\n"
        "satisfied, while the nearest few are satisfiable -- and once satisfied\n"
        "they change which of the rest are still wanted.");

    auto *repairGapBox = new QDoubleSpinBox(&dlg);
    repairGapBox->setRange(1e-7, 1e-1);
    repairGapBox->setDecimals(7);
    repairGapBox->setSingleStep(1e-4);
    repairGapBox->setValue(meridianConn_.repairGapLimit);
    repairGapBox->setToolTip(
        "How near a separatrix of Psi has to have passed a cone, as a fraction of\n"
        "the image extent, before the repair counts the pair as one that was meant\n"
        "to be joined. 1e-3 recommended.\n\n"
        "This is also the window the report's near-miss count is taken over, and\n"
        "the threshold at which the viewer colours a curve amber rather than blue.");

    auto *repairBoostBox = new QDoubleSpinBox(&dlg);
    repairBoostBox->setRange(1.0, 1000.0);
    repairBoostBox->setDecimals(1);
    repairBoostBox->setValue(meridianConn_.repairLambdaBoost);
    repairBoostBox->setToolTip("Sec. 3.3's \"raise lambda_5\": the factor the penalties are\n"
                               "multiplied by before each repair round. 10 is the paper's\n"
                               "growth factor and the recommended value.");

    auto *repairOuterBox = new QSpinBox(&dlg);
    repairOuterBox->setRange(1, 40);
    repairOuterBox->setValue(meridianConn_.repairOuterSteps);
    repairOuterBox->setToolTip("Outer penalty steps per repair round.");

    // ── the live panel ───────────────────────────────────────────────────────
    auto *derived = new QLabel(&dlg);
    derived->setTextFormat(Qt::PlainText);
    derived->setStyleSheet("font-family: monospace;");

    auto updateDerived = [&, nearMissBox, snapBox, repairGapBox, derived]() {
        const double seedTol = nearMissBox->value() * spacing;
        const double snapAbs = snapBox->value() * extent;
        const double gapAbs = repairGapBox->value() * extent;

        std::ostringstream oss;
        oss << std::scientific << std::setprecision(3);
        oss << "image extent          " << extent << "\n"
            << "mean cone spacing     " << spacing << "\n"
            << "closest cone pair     ";
        if (std::isfinite(closest)) oss << closest; else oss << "n/a (fewer than two cones)";
        oss << "\n"
            << "seeding tolerance     " << seedTol << "   (Stage 5)\n"
            << "snap tolerance        " << snapAbs << "   (Stage 7)\n"
            << "near-miss window      " << gapAbs << "   (repair)\n";

        // The one arithmetic mistake this dialog exists to prevent: a tolerance
        // wider than the cones are apart, which hands curves to whichever of a
        // clustered pair happens to be a hair nearer. The separation cap in
        // Separatrices stops it doing damage; saying so is better than letting
        // it be silently clamped.
        if (std::isfinite(closest)) {
            if (seedTol > 0.5 * closest) {
                oss << "\nThe seeding tolerance is more than half the distance between the two\n"
                       "closest cones. Constraints will be capped at a quarter of that\n"
                       "distance; consider a smaller fraction of the mean spacing.";
            }
            if (gapAbs > 0.5 * closest) {
                oss << "\nThe near-miss window is more than half the distance between the two\n"
                       "closest cones, so the repair may propose a pair that was never meant\n"
                       "to be joined.";
            }
        }
        if (ir.clusteredConePairs > 0) {
            oss << "\n" << ir.clusteredConePairs
                << " pair(s) of cones are clustered (Sec. 3.1). Q2 -- the cone angles --\n"
                   "is what will limit this model, not the connectivity: Eq. (15)'s 1/l_e\n"
                   "weighting leaves the largest angular error on the shortest edges, which\n"
                   "are the ones between clustered cones. The remedy is Stage 1's, to merge\n"
                   "the cluster into one cone of the summed index.";
        }
        derived->setText(QString::fromStdString(oss.str()));
    };
    QObject::connect(nearMissBox, &QDoubleSpinBox::valueChanged, &dlg, updateDerived);
    QObject::connect(snapBox, &QDoubleSpinBox::valueChanged, &dlg, updateDerived);
    QObject::connect(repairGapBox, &QDoubleSpinBox::valueChanged, &dlg, updateDerived);
    updateDerived();

    auto syncRepair = [repairBox, repairMaxBox, repairGapBox, repairBoostBox, repairOuterBox]() {
        const bool on = repairBox->value() > 0;
        repairMaxBox->setEnabled(on);
        repairGapBox->setEnabled(on);
        repairBoostBox->setEnabled(on);
        repairOuterBox->setEnabled(on);
    };
    QObject::connect(repairBox, &QSpinBox::valueChanged, &dlg, syncRepair);
    syncRepair();

    auto *restore = new QPushButton("Restore recommended", &dlg);
    QObject::connect(restore, &QPushButton::clicked, &dlg,
                     [&, nearMissBox, selfReturnBox, allConnBox, e6Box, propagateBox,
                      snapBox, ringsBox, stepsBox, cycleBox, repairBox, repairMaxBox,
                      repairGapBox, repairBoostBox, repairOuterBox]() {
        const SubdomainLabels::Options l;
        const Separatrices::Options t;
        const MERIDIAN::Options m;
        nearMissBox->setValue(l.nearMissTolerance);
        selfReturnBox->setChecked(l.seedSelfReturns);
        allConnBox->setChecked(l.seedAllConnections);
        e6Box->setChecked(l.interfaceCorners);
        propagateBox->setChecked(l.propagateInterfaceLabels);
        snapBox->setValue(t.coneSnapTolerance);
        ringsBox->setValue(t.coneSnapRings);
        stepsBox->setValue(t.maxSteps);
        cycleBox->setChecked(t.detectCycles);
        repairBox->setValue(m.repairPasses);
        repairMaxBox->setValue(m.repairMaxPerPass);
        repairGapBox->setValue(m.repairGapLimit);
        repairBoostBox->setValue(m.repairLambdaBoost);
        repairOuterBox->setValue(m.repairOuterSteps);
    });

    auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dlg);
    buttons->button(QDialogButtonBox::Cancel)->setText("Keep current");
    buttons->addButton(restore, QDialogButtonBox::ResetRole);
    QObject::connect(buttons, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
    QObject::connect(buttons, &QDialogButtonBox::rejected, &dlg, &QDialog::reject);

    auto section = [&dlg](const char *text) {
        auto *l = new QLabel(text, &dlg);
        l->setStyleSheet("font-weight: bold; margin-top: 8px;");
        return l;
    };

    // The form is tall enough to overflow a laptop screen, so it goes in a
    // scroll area with the buttons pinned below it rather than at the end of
    // the form, where they would be the first thing off the bottom.
    auto *page = new QWidget(&dlg);
    auto *form = new QFormLayout(page);
    form->addRow(section("Stage 5 — seeding Gamma_topo on psi_R"));
    form->addRow("near-miss tolerance (of cone spacing)", nearMissBox);
    form->addRow(selfReturnBox);
    form->addRow(allConnBox);

    if (multiMat) {
        form->addRow(section("Stage 0b — the material interfaces"));
        form->addRow(e6Box);
        form->addRow(propagateBox);
    }

    form->addRow(section("Stage 7 — tracing the separatrices of Psi"));
    form->addRow("snap tolerance (of image extent)", snapBox);
    form->addRow("snap search radius (triangles)", ringsBox);
    form->addRow("step cap", stepsBox);
    form->addRow(cycleBox);

    form->addRow(section("Sec. 3.3 — repair from the traced curves"));
    form->addRow("repair rounds", repairBox);
    form->addRow("constraints per round (0 = all)", repairMaxBox);
    form->addRow("near-miss window (of image extent)", repairGapBox);
    form->addRow("lambda growth per round", repairBoostBox);
    form->addRow("outer steps per round", repairOuterBox);

    form->addRow(section("This model"));
    form->addRow(derived);

    auto *scroll = new QScrollArea(&dlg);
    scroll->setWidget(page);
    scroll->setWidgetResizable(true);
    scroll->setFrameShape(QFrame::NoFrame);

    auto *outer = new QVBoxLayout(&dlg);
    outer->addWidget(scroll, 1);
    outer->addWidget(buttons, 0);

    // Fit the natural height where it fits, and scroll where it does not.
    const QRect avail = screen() ? screen()->availableGeometry() : QRect(0, 0, 1200, 900);
    dlg.resize(std::min(page->sizeHint().width() + 48, avail.width() - 80),
               std::min(page->sizeHint().height() + 80, avail.height() - 80));

    if (dlg.exec() != QDialog::Accepted) return false;

    meridianConn_.labels.nearMissTolerance = nearMissBox->value();
    meridianConn_.labels.seedSelfReturns = selfReturnBox->isChecked();
    meridianConn_.labels.seedAllConnections = allConnBox->isChecked();
    meridianConn_.labels.interfaceCorners = e6Box->isChecked();
    meridianConn_.labels.propagateInterfaceLabels = propagateBox->isChecked();
    meridianConn_.labels.maxTraceSteps = stepsBox->value();
    meridianConn_.trace.coneSnapTolerance = snapBox->value();
    meridianConn_.trace.coneSnapRings = ringsBox->value();
    meridianConn_.trace.maxSteps = stepsBox->value();
    meridianConn_.trace.detectCycles = cycleBox->isChecked();
    meridianConn_.trace.nearMissWindow = repairGapBox->value();
    meridianConn_.repairPasses = repairBox->value();
    meridianConn_.repairMaxPerPass = repairMaxBox->value();
    meridianConn_.repairGapLimit = repairGapBox->value();
    meridianConn_.repairLambdaBoost = repairBoostBox->value();
    meridianConn_.repairOuterSteps = repairOuterBox->value();
    return true;
}

void CrossGenWidget::runMERIDIANLayout() {
    if (!coneCut_.has_value() || !cones_.has_value()) {
        // Said once -- pipelineBlocked_ is only empty the first time through --
        // because this is reached on every frame for as long as it is true.
        if (pipelineBlocked_.empty())
            blockPipeline("Stages 4-11 not run",
                          "there is no cutting graph to unfold: Stage 2 did not produce "
                          "one, and psi_R is an immersion of the disk it cuts");
        return;
    }
    layoutAttempted_ = true;
    pipelineBlocked_.clear();

    // Fit the right panel to whichever map is about to be shown there.
    auto fitPanel = [this](const std::vector<Point> &uv) {
        if (uv.empty()) return;
        viewer::computeLayoutBounds(uv, uvView_.cx, uvView_.cy,
                                    uvView_.baseW, uvView_.baseH);
        uvView_.zoom = 1.0;
        uvView_.fbw  = view_.fbw;
        uvView_.fbh  = view_.fbh;
    };

    // ── Stage 4: the immersion psi_0 ─────────────────────────────────────────
    //
    // The one place the two routes part company, and the only one in this
    // function. Pipeline A unfolds the flat cone metric here, triangle by
    // triangle; Pipeline B has no metric to unfold and built its immersion at
    // the phase before, out of the map the integration produced. Everything
    // below reads immersion_ and cannot tell which of the two filled it, which
    // is the whole point of Sec. 6.2 calling the Immersion the hub.
    auto t0 = Clock::now();
    auto t1 = t0;

    if (mode_ == Mode::TORSION) {
        if (!immersion_.has_value()) return;
        // The same gate the Ricci route applies below, applied to the map the
        // integration produced: the barrier of Sec. 3.3 preserves Q1 and cannot
        // repair it, so a psi_0 with a flipped face is not something to start a
        // continuation from. Stage 4R is what stands in front of this, and when
        // it did not manage it the phase before said so.
        if (immersion_->getReport().flippedFaces > 0) {
            blockPipeline("Stages 5-11 not run",
                          "psi_0 does not satisfy Q1, and Sec. 3.3's barrier can only "
                          "preserve Q1, never repair it; see the Stage 4F and 4R reports "
                          "at the phase before");
            return;
        }
        psiR_ = immersion_->getUV();
        showPsiR_ = false;
        fitPanel(psiR_);
    } else {
    if (!ricci_.has_value()) return;

    t0 = Clock::now();
    try {
        immersion_.emplace(*coneCut_, *ricci_, *cones_);
    } catch (const std::exception &e) {
        immersion_.reset();
        console_.log(std::string("[Immersion] FAILED: ") + e.what());
        std::cerr << "[Viewer] Immersion failed: " << e.what() << "\n";
        return;
    }
    t1 = Clock::now();

    const Immersion::Report &ir = immersion_->getReport();
    psiR_ = immersion_->getUV();
    showPsiR_ = false;
    fitPanel(psiR_);

    {
        std::ostringstream oss;
        oss << "[Immersion] " << ir.placedVertices << "/" << ir.cutVertices
            << " vertices placed, image " << std::fixed << std::setprecision(3)
            << (ir.uvMax[0] - ir.uvMin[0]) << " x " << (ir.uvMax[1] - ir.uvMin[1])
            << (ir.mirrored ? " (reflected)" : "") << ", "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Immersion] metric realised to " << std::scientific << std::setprecision(2)
            << ir.maxMetricResidual << " relative, closure gap " << ir.maxClosureGap
            << "; Q1 " << (ir.flippedFaces == 0 ? "[PASS]" : "[FAIL]") << " ("
            << ir.flippedFaces << " flipped face(s))";
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }
    {
        std::ostringstream oss;
        oss << "[Immersion] " << ir.arcs << " arc(s), " << ir.seamEdgePairs
            << " seam pair(s); Gamma_Hol_0..3 = " << ir.holonomyCount[0] << ", "
            << ir.holonomyCount[1] << ", " << ir.holonomyCount[2] << ", "
            << ir.holonomyCount[3] << "; Q4 snap " << std::scientific
            << std::setprecision(2) << ir.maxSnapError << " rad";
        console_.log(oss.str());
    }
    if (ir.clusteredConePairs > 0) {
        std::ostringstream oss;
        oss << "[Immersion] " << ir.clusteredConePairs << " clustered cone pair(s), closest "
            << std::scientific << std::setprecision(2) << ir.minConeSeparation
            << " of the model apart -- Q2 usually ends up limited by this";
        console_.log(oss.str());
    }
    for (const std::string &m : ir.messages) console_.log("[Immersion] " + m);

    if (ir.flippedFaces > 0 || ir.unplacedVertices > 0) {
        blockPipeline("Stages 5-11 not run",
                      "psi_R does not satisfy Q1 (" + std::to_string(ir.flippedFaces) +
                          " flipped face(s), " + std::to_string(ir.unplacedVertices) +
                          " unplaced vertex/vertices); the continuation can only preserve "
                          "Q1, never repair it");
        return;
    }
    }  // end of the Ricci route's Stage 4

    // The connectivity settings are asked for here rather than on the way into
    // the mode: every number in them is a fraction of something psi_R measures
    // -- the spacing of the cones, the extent of the image -- and none of those
    // exist until Stage 4 has run. Cancelling keeps whatever was last used, so
    // the pipeline runs either way.
    if (!meridianConnPrompted_) {
        meridianConnPrompted_ = true;
        promptMERIDIANConnectivity(*immersion_, psiR_);
    }

    // ── Stage 5: subdomain labelling ─────────────────────────────────────────
    t0 = Clock::now();
    try {
        meridianLabels_.emplace(*immersion_, meridianConn_.labels,
                                interfaces_.has_value() ? &*interfaces_ : nullptr);
    } catch (const std::exception &e) {
        meridianLabels_.reset();
        console_.log(std::string("[Labels] FAILED: ") + e.what());
        std::cerr << "[Viewer] SubdomainLabels failed: " << e.what() << "\n";
        return;
    }
    t1 = Clock::now();

    const SubdomainLabels::Report &sr = meridianLabels_->getReport();
    {
        std::ostringstream oss;
        oss << "[Labels] dS: " << sr.boundaryEdges << " edge(s) -> Gamma_u "
            << sr.boundaryEdgesU << ", Gamma_v " << sr.boundaryEdgesV << " in "
            << sr.boundaryChains << " chain(s)";
        if (sr.ambiguousBoundaryEdges > 0) oss << " (" << sr.ambiguousBoundaryEdges << " near-tied)";
        oss << ", " << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Labels] features: " << sr.featureEdges << " edge(s) in " << sr.featureChains
            << " chain(s); Gamma_topo: " << sr.topoPaths << " path(s) from "
            << sr.separatrices << " separatrix/-ces (" << sr.separatricesToCone
            << " to a cone, " << sr.separatricesCapped << " unresolved)";
        console_.log(oss.str());
    }
    {
        // The two kinds of constraint that are easy to leave out, counted, so
        // that turning either of them off in the dialog shows in the log rather
        // than only in the layout two stages later.
        std::ostringstream oss;
        oss << "[Labels] of those paths, " << sr.topoSelfReturns
            << " return to the cone they left (Q5's \"possibly identical\" case, Fig. 9) and "
            << sr.topoExtraPerPair << " are a second curve between a pair already joined;"
            << " near-miss tolerance " << std::fixed << std::setprecision(3)
            << meridianConn_.labels.nearMissTolerance << " of the mean cone spacing = "
            << std::scientific << std::setprecision(2) << sr.seedSnapTolerance
            << " in image units, searched over " << sr.seedSnapRings << " ring(s)";
        console_.log(oss.str());
    }
    if (sr.interfaceCorners > 0 || sr.featureChainsPropagated > 0) {
        // What Stage 6 will be holding the interfaces to. The corrected count
        // is the one worth watching: a chain the flux of psi_R would have
        // labelled the other way is one where the node quantisation and the map
        // disagree, and E3 alone would have taken the map's word for it.
        std::ostringstream oss;
        oss << "[Labels] interfaces: " << sr.featureChainsPropagated
            << " chain(s) labelled from the node quantisation ("
            << sr.featureLabelsCorrected << " against what their flux said), "
            << sr.interfaceCorners << " sector(s) handed to E6";
        if (sr.interfaceCornersSpanningCut > 0) {
            oss << ", " << sr.interfaceCornersSpanningCut
                << " dropped where the cutting graph runs through the sector";
        }
        console_.log(oss.str());
    }
    for (const std::string &m : sr.messages) console_.log("[Labels] " + m);

    // ── Stage 6: the layout-inducing energies ────────────────────────────────
    //
    // C4, and the second and last place the routes differ. E1 measures J
    // against a reference metric, and a field-integrated map has no flat metric
    // to be measured against; leaving it Euclidean is the failure
    // LayoutEnergy.hxx documents on geom003, where a reference with no cones
    // and Q2 are contradictory statements about the same vertex.
    //
    // psi_0's own induced metric is not the answer either, though it looks like
    // one: the seam constraints are exact rotations, so its angle sum at each
    // cone is already 2 pi - (pi/2) I, and C4 reads 1e-10 on it. The cone
    // angles are only the part of the reference that lives at the cones, and
    // the rest of it is the shape of every other triangle -- E1 against psi_0's
    // lengths says "reproduce psi_0", shear and all, so the cone fans stay
    // however the least-squares fit distributed them and the patch at a cone
    // comes out with a reflex corner. What is used is Sec. 4's flat cone
    // metric, built from the model's own triangles by ConeMetric.
    LayoutEnergy::Options eopts;
    if (mode_ == Mode::TORSION) {
        eopts.reference = LayoutEnergy::Reference::Induced;
        eopts.referenceLengths =
            coneMetric_.has_value()
                ? coneMetric_->edgeLengths()
                : TORSION::inducedLengths(*mesh_, *coneCut_, immersion_->getUV());
        // Sec. 9's closing note: the field is already aligned to dS and to the
        // interfaces, so E2 and E3 start small and over-penalising them early
        // only fights E1; E4 starts higher because the untangling introduced
        // seam error, and Q4 is not free to trade away.
        eopts.lambdaFactor[2] = TORSION::Options().lambdaAlignmentFactor;
        eopts.lambdaFactor[3] = TORSION::Options().lambdaAlignmentFactor;
        eopts.lambdaFactor[4] = TORSION::Options().lambdaSeamFactor;
    }

    t0 = Clock::now();
    bool ok = false;
    try {
        meridianLayout_.emplace(*immersion_, *meridianLabels_, eopts);
        ok = meridianLayout_->run();
    } catch (const std::exception &e) {
        meridianLayout_.reset();
        console_.log(std::string("[Layout] FAILED: ") + e.what());
        std::cerr << "[Viewer] LayoutEnergy failed: " << e.what() << "\n";
        return;
    }
    t1 = Clock::now();

    const LayoutEnergy::Report &er = meridianLayout_->getReport();
    fitPanel(meridianLayout_->getUV());

    {
        std::ostringstream oss;
        oss << "[Layout] " << er.outerSteps << " outer step(s), " << er.innerIterations
            << " inner iteration(s), " << er.referenceSwitches << " reference switch(es), "
            << er.relabels << " relabel(s), "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Layout] Q3 boundary " << std::scientific << std::setprecision(2)
            << er.initialBoundaryResidual << " -> " << er.maxBoundaryResidual
            << ", Q4 seam " << er.initialSeamResidual << " -> " << er.maxSeamResidual
            << ", Q5 topo " << er.initialTopoResidual << " -> " << er.maxTopoResidual;
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Layout] Q1 det J >= " << std::fixed << std::setprecision(4) << er.minDetJ
            << " (" << er.invertedTriangles << " inverted) " << (er.injective ? "[PASS]" : "[FAIL]")
            << "; Q2 cone angles off by " << std::scientific << std::setprecision(2)
            << er.maxConeAngleResidual << " rad " << (er.anglesHeld ? "[PASS]" : "[FAIL]");
        if (er.coneValenceChanges > 0) oss << ", " << er.coneValenceChanges << " changed valence";
        console_.log(oss.str());
    }
    if (er.interfaceCorners > 0) {
        std::ostringstream oss;
        oss << "[Layout] E6: the worst interface sector is " << std::scientific
            << std::setprecision(2) << er.maxInterfaceResidual
            << " rad from the quarter turns the model has there";
        if (er.interfaceCornerChanges > 0) {
            oss << "; " << er.interfaceCornerChanges
                << " turned a different whole number of quarters, which is a layout the "
                << "material regions do not have";
        }
        oss << " " << (er.interfaceCornerChanges == 0 && er.maxInterfaceResidual < 1e-3
                           ? "[PASS]" : "[FAIL]");
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Layout] Psi " << (ok ? "satisfies Q1-Q5: a quadrilateral layout in the sense "
                                        "of Definition 2.1 [PASS]"
                                      : "is not yet a layout; see the residuals above [FAIL]");
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }
    for (const std::string &m : er.messages) console_.log("[Layout] " + m);
    console_.log("[Layout] right panel = Psi; press 'p' for psi_R, the map it started from. "
                 "Blue = Gamma_u, green = Gamma_v, amber/magenta = the two banks of each seam, "
                 "white-cored = a feature chain (on a multi-material model, an interface)");
}

// Stage 7, Sec. 4. The separatrices are what Q5 is actually about -- every
// integral curve of Psi out of a cone is finite -- and Stage 6 can only assert
// that through the residual of E5, which is a statement about the Gamma_topo
// paths Stage 5 happened to seed. This traces every curve there is and reads
// the property off all of them.
//
// It runs on Psi rather than on psi_R even when the continuation fell short,
// and it is worth running then: a curve that reaches the step cap names the
// pair of cones whose connectivity constraint is missing, which is the
// diagnosis Sec. 3.3 asks for when the layout is not yet a layout.
void CrossGenWidget::runMERIDIANSeparatrices() {
    if (!immersion_.has_value() || !meridianLayout_.has_value() ||
        !meridianLabels_.has_value()) return;
    separatricesAttempted_ = true;
    pipelineBlocked_.clear();

    auto t0 = Clock::now();

    MERIDIAN::RepairOptions ropts;
    ropts.passes = meridianConn_.repairPasses;
    ropts.maxPerPass = meridianConn_.repairMaxPerPass;
    ropts.gapLimit = meridianConn_.repairGapLimit;
    ropts.lambdaBoost = meridianConn_.repairLambdaBoost;
    ropts.outerSteps = meridianConn_.repairOuterSteps;

    Separatrices::Options topts = meridianConn_.trace;
    topts.nearMissWindow = meridianConn_.repairGapLimit;
    // The nodes of the interface network emit as well as the cones -- but only
    // the ones with spare sector capacity. A node with a sector three quarters
    // wide needs two layout edges into it that no cone emits, and without them
    // the region outside the corner comes back as one face with five corners
    // that Stage 9 cannot fit and Stage 10 cannot mesh. A node with no spare
    // sector (every quarter already accounted for) has no ray to receive a
    // curve, so it must not be a termination target either -- the same rule
    // Stage 5's own seeding uses (Interfaces::emitterNodes()), shared here so
    // Stage 5 and Stage 7 never disagree about which nodes are live. The rays
    // that would run *along* an interface are suppressed by the tracer itself:
    // the interface is already an arc of the layout.
    if (interfaces_.has_value() && interfaces_->multiMaterial()) {
        const std::vector<int> emitters = interfaces_->emitterNodes();
        topts.extraEmitters.insert(topts.extraEmitters.end(), emitters.begin(), emitters.end());
    }

    MERIDIAN::RepairResult rep;
    try {
        rep = MERIDIAN::traceAndRepair(*meridianLabels_, *meridianLayout_, topts, ropts,
                                       [this](const std::string &m) {
                                           // "Stage 5: ..." and friends, retagged
                                           // with the names the console already
                                           // uses for those stages.
                                           static const char *kTag[3] = {
                                               "[Labels] ", "[Layout] ", "[Separatrices] "};
                                           const size_t colon = m.find(": ");
                                           if (colon == std::string::npos || m.size() < 7) {
                                               console_.log(m);
                                               return;
                                           }
                                           const int stage = m[6] - '5';
                                           console_.log((stage >= 0 && stage < 3 ? kTag[stage]
                                                                                 : "[MERIDIAN] ") +
                                                        m.substr(colon + 2));
                                       });
    } catch (const std::exception &e) {
        separatrices_.reset();
        blockPipeline("Stages 7-11 not run", e.what());
        return;
    }
    auto t1 = Clock::now();

    if (!rep.separatrices) {
        separatrices_.reset();
        blockPipeline("Stages 7-11 not run", "the curves could not be traced");
        return;
    }
    separatrices_.emplace(std::move(*rep.separatrices));

    // The repair moves Psi, so the panel fitted to the map Stage 6 first
    // returned is fitted to the wrong one. Refit before the first frame; and if
    // the Layout phase was left showing psi_R, this phase always shows Psi.
    showPsiR_ = false;
    viewer::computeLayoutBounds(meridianLayout_->getUV(), uvView_.cx, uvView_.cy,
                                uvView_.baseW, uvView_.baseH);
    uvView_.zoom = 1.0;

    const Separatrices::Report &r = separatrices_->getReport();
    {
        std::ostringstream oss;
        oss << "[Separatrices] " << r.emitted << " curve(s) from " << r.cones
            << " cone(s), the indices prescribe " << r.prescribed
            << " (4 - I interior, 1 - I on dS), "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Separatrices] ends: " << r.endedAtCone << " at a cone, " << r.endedAtBoundary
            << " out through dS, " << r.cycled << " on a closed orbit, " << r.capped
            << " at the step cap";
        if (r.stuck > 0 || r.degenerate > 0)
            oss << ", " << r.stuck << " stuck, " << r.degenerate << " degenerate";
        oss << "; " << r.triangleSteps << " triangle crossing(s), " << r.seamCrossings
            << " seam crossing(s)";
        console_.log(oss.str());
    }
    {
        // The two numbers that say the tracing itself is sound, independently of
        // whether the layout is: how far a snap had to reach compared with the
        // tolerance it was allowed, and whether the pullback of each curve is
        // continuous where it crosses the cut, which is the only place a wrong
        // transition would show.
        std::ostringstream oss;
        oss << "[Separatrices] snap tolerance " << std::scientific << std::setprecision(2)
            << r.snapTolerance << " of the image, worst snap taken " << r.maxSnapGap
            << "; pullback continuity " << r.maxPullbackGap << " of the model";
        console_.log(oss.str());
    }
    if (rep.passesTaken > 0 || rep.constraintsAdded > 0 || rep.rolledBack > 0) {
        std::ostringstream oss;
        oss << "[Repair] " << rep.passesTaken << " round(s) of Sec. 3.3's remedy added "
            << rep.constraintsAdded << " connectivity constraint(s) and re-ran Stage 6 from "
            << "the current phi; Gamma_topo now has "
            << meridianLabels_->getReport().topoPaths << " path(s)";
        if (rep.rolledBack > 0) oss << " (" << rep.rolledBack << " round taken back as worse)";
        console_.log(oss.str());

        const LayoutEnergy::Report &er = meridianLayout_->getReport();
        std::ostringstream oss2;
        oss2 << "[Repair] Q3 " << std::scientific << std::setprecision(2)
             << er.maxBoundaryResidual << ", Q4 " << er.maxSeamResidual << ", Q5 "
             << er.maxTopoResidual << ", det J >= " << std::fixed << std::setprecision(4)
             << er.minDetJ << "; Psi "
             << (er.valid ? "satisfies Q1-Q5 [PASS]" : "is not yet a layout [FAIL]");
        console_.log(oss2.str());
    }
    if (r.nearMisses + r.grazes > 0) {
        // The thing the eye cannot pick out of the picture on its own. These
        // are the amber curves.
        std::ostringstream oss;
        oss << "[Separatrices] " << r.nearMisses << " unterminated and " << r.grazes
            << " boundary-bound curve(s) passed within " << std::scientific
            << std::setprecision(2) << r.nearMissWindow
            << " of the image extent of a cone without stopping at it -- each is a "
            << "quadrilateral of poor aspect ratio (Remark 3.1). Press 'c' to re-open the "
            << "connectivity dialog and try a wider near-miss window or more repair rounds";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Separatrices] Q5 " << (r.valid ? "holds on every curve traced [PASS]"
                                                : "not yet verified; see the ends above [FAIL]");
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }
    console_.log("[Separatrices] left = the curves on the model (Fig. 9), right = the same "
                 "curves on Psi. Green ends at a cone, blue leaves through dS, amber grazed "
                 "a cone and went on, red ended nowhere");
}

// Whether Stage 7 came back with something a layout can be read off: every
// curve it traced terminated the way Q5 allows, and Psi is still injective.
// Stage 8's own validation is a different question -- it asks whether the partition is
// four-sided everywhere -- and it is allowed to fail and be reported. What is
// not allowed is to build the partition out of curves that never closed, or on
// a map that folded.
//
// LayoutEnergy::Report::valid is deliberately *not* the test, though it looks
// like the obvious one. It is Q1, Q2 and `constraintsMet`, and that last is the
// continuation's own residual on the Gamma_topo paths Stage 5 seeded -- a
// statement about how far E5 got, not about the curves. Q5 itself is what Stage
// 7 measures, on every separatrix there is rather than on the seeded subset,
// and that is the verdict above. Requiring `valid` on top of it means a model
// whose residual stops just short has its arrangement, its splines and its mesh
// refused over a number that says nothing about whether the curves partition S:
// data/meshes/multimat/bubbles with its inclusions excised is exactly that, Q5
// passing on all 200 traced curves with the E5 residual at 4e-6 against a
// tolerance of 1e-6, and the viewer showed a blank screen for the last two
// phases of a run TestMERIDIAN takes all the way to the templates.
//
// MERIDIAN::run says why in as many words at its own Stage 8: run it even when
// Stage 6 fell short, because the arrangement is where the failure becomes a
// *place* -- this face has three corners, that cone is one arc short -- rather
// than a residual. Q2 likewise: cone angles a little off make some sectors
// something other than quarters, which Stage 8 counts and reports face by face.
// The one thing that cannot be allowed is what this still refuses: curves that
// never closed, and a Psi that folded, on which an arrangement is not a
// partition of anything.
bool CrossGenWidget::meridianTraceIsClean() const {
    if (!separatrices_.has_value() || !meridianLayout_.has_value()) return false;
    const Separatrices::Report &r = separatrices_->getReport();
    if (!r.valid) return false;
    if (r.emitted == 0) return false;
    return meridianLayout_->getReport().injective;
}

// A stage that has nothing to hand the phase after it, said in the three places
// it has to be said: the console, which is where the rest of the run is; the
// terminal, which keeps a transcript the console cannot; and the overlay, which
// is the only one still on screen when the blank phase is.
void CrossGenWidget::blockPipeline(const std::string &stage, const std::string &what) {
    console_.log("[" + stage + "] " + what);
    std::cerr << "[Viewer] " << stage << ": " << what << "\n";
    pipelineBlocked_ = stage + " -- " + what;
}

// Stages 8 and 9: the curves become a planar subdivision of S, and each arc of
// that subdivision becomes the one cubic B-spline both of its faces share.
//
// The two run together because the second is only meaningful on the first and
// neither is worth a phase of its own to look at: Stage 9 moves the lines of
// Stage 8 by less than the width they are drawn at, which is the point -- the
// fit is a reconstruction of the layout, not a change to it. The numbers are
// where they differ, so both reports are printed in full.
void CrossGenWidget::runMERIDIANPatches() {
    // Not "attempted" until the curves it is to be built from are there. This
    // is called from advancePhase() as well as from runComputations(), and the
    // window between them is real: Stage 6 and Stage 7 hold the GUI thread for
    // seconds on a model with forty cones, every keypress made in the meantime
    // is delivered when they return, and two of them walk the phase past this
    // stage before Stage 7 has produced anything. Setting the flag here would
    // then retire Stage 8 for the rest of the run -- the arrangement, the
    // splines, the mesh and the templates all silently absent, on a phase whose
    // whole picture they are.
    if (!separatrices_.has_value() || !meridianLabels_.has_value()) return;
    patchesAttempted_ = true;

    pipelineBlocked_.clear();
    if (!meridianTraceIsClean()) {
        // Refusing rather than drawing a partition that is not one -- and
        // naming which of the reasons it was, because their remedies differ:
        // curves that never closed are what the connectivity dialog can do
        // something about, and a folded Psi is not.
        const Separatrices::Report &sr = separatrices_->getReport();
        std::ostringstream oss;
        if (!meridianLayout_.has_value()) {
            oss << "Stage 6 produced no layout to read";
        } else if (sr.emitted == 0) {
            oss << "Stage 7 emitted no separatrices";
        } else if (!sr.valid) {
            oss << "Q5 failed on the curves Stage 7 traced, so they do not partition S -- "
                   "press 'n' to re-open the connectivity dialog and trace again";
        } else {
            oss << "Psi folded (" << meridianLayout_->getReport().invertedTriangles
                << " inverted triangle(s)), and an arrangement of curves on a map that is "
                   "not injective is not a partition";
        }
        blockPipeline("Stage 8 not built", oss.str());
        return;
    }

    // Built, but on a layout that fell short: said here rather than left to be
    // inferred from Stage 8's own counts, because the picture that follows
    // looks exactly like one from a layout that did not.
    {
        const LayoutEnergy::Report &er = meridianLayout_->getReport();
        if (!er.constraintsMet || !er.anglesHeld) {
            std::ostringstream oss;
            oss << "[Patches] built on a layout Stage 6 did not close out -- ";
            if (!er.constraintsMet)
                oss << "the E5 residual stopped at " << std::scientific << std::setprecision(2)
                    << er.maxTopoResidual << " of the image extent";
            if (!er.constraintsMet && !er.anglesHeld) oss << ", and ";
            if (!er.anglesHeld)
                oss << "a cone is " << std::scientific << std::setprecision(2)
                    << er.maxConeAngleResidual << " rad off the angle Stage 1 prescribed";
            oss << ". Every curve Stage 7 traced still satisfies Q5, so the partition is "
                   "there to look at; where it went wrong is a face with the wrong number "
                   "of corners below, not a residual";
            console_.log(oss.str());
        }
    }

    auto t0 = Clock::now();
    try {
        arrangement_.emplace(*separatrices_, *meridianLabels_);
    } catch (const std::exception &e) {
        arrangement_.reset();
        blockPipeline("Stage 8 failed", e.what());
        return;
    }
    auto t1 = Clock::now();

    const Arrangement::Report &ar = arrangement_->getReport();
    {
        std::ostringstream oss;
        oss << "[Patches] Stage 8: " << ar.nodes << " node(s) (" << ar.coneNodes
            << " cone, " << ar.crossingNodes << " crossing, " << ar.boundaryHitNodes
            << " on dS, " << ar.boundaryCornerNodes << " corner) and " << ar.arcs
            << " arc(s) (" << ar.separatrixArcs << " traced, " << ar.boundaryArcs
            << " of dS), " << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }
    if (ar.interfaceArcs > 0 || ar.interfaceNodes > 0) {
        // mixedPatches is the property the whole multi-material path exists to
        // produce, so it is reported whether it is zero or not: a patch whose
        // interior is not all one material is one every element of which will
        // straddle the interface.
        std::ostringstream oss;
        oss << "[Patches] interfaces: " << ar.interfaceArcs << " arc(s) and "
            << ar.interfaceNodes << " node(s) of the network in the arrangement, "
            << ar.interfaceHitNodes << " separatrix/-ces landed on one; "
            << ar.mixedPatches << " patch(es) of mixed material "
            << (ar.mixedPatches == 0 ? "[PASS]" : "[FAIL]");
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Patches] " << ar.patches << " patch(es), " << ar.quads
            << " with four corners, " << ar.simpleQuads << " with one arc per side; they cover "
            << std::fixed << std::setprecision(6) << ar.areaCoverage << " of S";
        console_.log(oss.str());
    }
    if (ar.tJunctions > 0 || ar.wrongCornerFaces > 0 || ar.coneValenceErrors > 0 ||
        ar.danglingArcs > 0 || ar.unsharedArcs > 0) {
        std::ostringstream oss;
        oss << "[Patches] Sec. 4's validation: " << ar.tJunctions << " T-junction(s), "
            << ar.wrongCornerFaces << " face(s) without four corners, " << ar.coneValenceErrors
            << " cone(s) of the wrong valence, " << ar.danglingArcs << " dangling arc(s), "
            << ar.unsharedArcs << " unshared arc(s) -- the quantisation of Stage 6 has not "
            << "finished on this model";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Patches] Stage 8 " << (ar.valid ? "validates [PASS]"
                                                 : "does not validate [FAIL]");
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }

    // The one number that decides whether the model comes out meshed: the area
    // of S in faces Stage 10 has no grid for. MERIDIAN::run answers it by
    // seeding Gamma_topo again, tighter, and keeping whichever attempt leaves
    // less. The same retry is run here, once, in place -- a tolerance that
    // over-seeds is not recoverable by the repair loop, which can only add
    // constraints, and a tighter one is, which is why it runs downward. What
    // the viewer does not do is keep the first attempt: Stage 8 holds pointers
    // into Stages 5 to 7 and putting those back would mean holding two of each,
    // so the second is reported against the first and 'n' is the way back.
    {
        const double unmeshable = MERIDIAN::unmeshableFraction(*arrangement_);
        const double retryTol = MERIDIAN::Options().topoNearMissRetry;
        if (unmeshableBefore_ >= 0.0) {
            std::ostringstream oss;
            oss << "[Patches] re-seeded: " << std::fixed << std::setprecision(2)
                << 100.0 * unmeshableBefore_ << "% of S had no grid at a near-miss "
                << "tolerance of " << std::setprecision(3) << nearMissBefore_ << ", "
                << std::setprecision(2) << 100.0 * unmeshable << "% at "
                << std::setprecision(3) << meridianConn_.labels.nearMissTolerance;
            if (unmeshable > unmeshableBefore_) {
                oss << " -- the retry left more of S uncovered, so press 'n' and set the "
                    << "tolerance back to " << nearMissBefore_;
            }
            console_.log(oss.str());
            unmeshableBefore_ = -1.0;
        } else if (unmeshable > 0.0 && !nearMissRetried_ && retryTol > 0.0 &&
                   meridianConn_.labels.nearMissTolerance > retryTol) {
            nearMissRetried_ = true;
            unmeshableBefore_ = unmeshable;
            nearMissBefore_ = meridianConn_.labels.nearMissTolerance;
            std::ostringstream oss;
            oss << "[Patches] " << std::fixed << std::setprecision(2) << 100.0 * unmeshable
                << "% of S is in faces Stage 10 has no grid for. Seeding Gamma_topo again at "
                << std::setprecision(3) << retryTol << " and tracing from there; this blocks";
            console_.log(oss.str());
            std::cerr << "[Viewer] " << oss.str() << "\n";
            meridianConn_.labels.nearMissTolerance = retryTol;
            // Everything from Stage 5 down is rebuilt, this arrangement with
            // it, so nothing below may touch it and the frame after this one
            // builds Stage 8 again.
            rerunMERIDIANConnectivity();
            return;
        }
    }

    auto t2 = Clock::now();
    try {
        splines_.emplace(*arrangement_);
    } catch (const std::exception &e) {
        splines_.reset();
        blockPipeline("Stage 9 failed", e.what());
        return;
    }
    auto t3 = Clock::now();

    const SplineFit::Report &sr = splines_->getReport();
    {
        std::ostringstream oss;
        oss << "[Patches] Stage 9: " << sr.curves << " cubic B-spline(s) of "
            << sr.controlPointsPerArc << " control points and " << sr.patches
            << " bicubic patch(es) of " << sr.controlPointsPerArc << "x"
            << sr.controlPointsPerArc << ", "
            << formatMs(std::chrono::duration<double, std::milli>(t3 - t2).count());
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Patches] the fit is off the traced curves by at most " << std::scientific
            << std::setprecision(2) << sr.maxDeviation << " of the model (rms "
            << sr.rmsDeviation << ")";
        console_.log(oss.str());
    }
    {
        // The watertightness rule of Sec. 5, measured: an arc is fitted once
        // and both of its patches are handed the same control points, so the
        // gap between the two copies is zero exactly or the rule was broken.
        std::ostringstream oss;
        oss << "[Patches] shared control points apart by " << std::scientific
            << std::setprecision(3) << sr.maxSeamGap << " over " << sr.sharedArcs
            << " shared arc(s): "
            << (sr.watertight ? "watertight [PASS]" : "not watertight [FAIL]");
        console_.log(oss.str());
    }
    if (sr.foldedPatches > 0) {
        std::ostringstream oss;
        oss << "[Patches] " << sr.foldedPatches
            << " patch(es) fold: the Coons blend of Sec. 5 is bilinear in its four sides and "
            << "a layout face whose sides bow outwards has no such surface. More Bezier "
            << "segments do not fix it";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Patches] Stage 9 " << (sr.valid ? "validates [PASS]"
                                                 : "does not validate [FAIL]");
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }
    console_.log("[Patches] the blocks on the model: light blue is a patch side -- the fitted "
                 "spline where there is one, the traced arc where the fit was skipped -- and "
                 "green is a node of the arrangement");
}

// ── MERIDIAN Stage 10: the mesh dialog ───────────────────────────────────────
//
// One number matters here and the dialog says so: the target edge length, in
// the units the model is in. Everything else on it is either a bound that is
// almost never wanted or a switch for telling a meshing artefact apart from a
// fitting one.
//
// The derived line under the spin box is the point of opening a dialog at all
// rather than hard-coding 0.05. A target is an absolute length, and whether it
// is a sensible one is a question about this model's size -- so the diagonal of
// S and the number of edges that target puts across it are shown next to it,
// and they move as the number does.
bool CrossGenWidget::promptMERIDIANMesh() {
    if (!splines_.has_value() || !arrangement_.has_value() || !mesh_) return false;

    // The same diagonal QuadMesh measures for itself, off the same vertices.
    Point lo{ std::numeric_limits<double>::infinity(),
              std::numeric_limits<double>::infinity() };
    Point hi{ -lo[0], -lo[1] };
    for (const Point &p : mesh_->vertices) {
        lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
        hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
    }
    const double diag = std::hypot(hi[0] - lo[0], hi[1] - lo[1]);
    const double extent = (diag > 0.0) ? diag : 1.0;

    QDialog dlg(this);
    dlg.setWindowTitle("MERIDIAN — Stage 10, quadrilateral mesh");

    auto *targetBox = new QDoubleSpinBox(&dlg);
    targetBox->setRange(1e-4, 10.0);
    targetBox->setDecimals(4);
    targetBox->setSingleStep(0.01);
    targetBox->setValue(meshSettings_.target);
    targetBox->setToolTip(
        "Target length of a mesh edge, in the units of the model.\n"
        "An absolute length, not a fraction: it is the same number\n"
        "whatever the model is, and the line below says what it means here.");

    auto *derived = new QLabel(&dlg);
    derived->setTextFormat(Qt::PlainText);

    const int patches = static_cast<int>(splines_->patches().size());
    auto updateDerived = [targetBox, derived, extent, patches]() {
        const double h = targetBox->value();
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(4)
            << "diagonal of S                 = " << extent << "\n"
            << "edges across it   diag / h    = " << std::setprecision(1)
            << (extent / h) << "\n"
            << patches << " patch(es) to mesh";
        derived->setText(QString::fromStdString(oss.str()));
    };
    QObject::connect(targetBox, &QDoubleSpinBox::valueChanged, &dlg, updateDerived);
    updateDerived();

    // The floor is what keeps a chord across a thin sliver from vanishing; the
    // ceiling is a guard against a target so fine the mesh will not fit in
    // memory. Neither is normally touched, and a chord that hits one is
    // reported as clamped, because the interval assignment did not choose it.
    auto *minBox = new QSpinBox(&dlg);
    minBox->setRange(1, 64);
    minBox->setValue(meshSettings_.minEdges);
    minBox->setToolTip("Fewest edges any chord may be given.");

    auto *maxBox = new QSpinBox(&dlg);
    maxBox->setRange(0, 4096);
    maxBox->setValue(meshSettings_.maxEdges);
    maxBox->setSpecialValueText("none");
    maxBox->setToolTip("Most edges any chord may be given. 0 for no ceiling.");

    auto *collapseBox = new QDoubleSpinBox(&dlg);
    collapseBox->setRange(0.0, 2.0);
    collapseBox->setDecimals(2);
    collapseBox->setSingleStep(0.05);
    collapseBox->setValue(meshSettings_.collapseSpan);
    collapseBox->setSpecialValueText("off");
    collapseBox->setToolTip(
        "Contract a chord every patch of which is thinner than this many target\n"
        "edge lengths, and let the blocks either side of it meet.\n"
        "No integer assignment can make an element bigger than the block it sits\n"
        "in, so where the layout is finer than the target it is the layout, and\n"
        "not the target, that sets the element size — on bubbles a tenth of the\n"
        "faces are under a sixteenth of one element in area. At 0.5 the trade is\n"
        "even: keeping the chord puts every element on it below half the target,\n"
        "contracting it moves the mesh by less than half an element. Off keeps\n"
        "every chord and is what shows the difference.");

    auto *splineBox = new QCheckBox("place the nodes on the Stage 9 splines", &dlg);
    splineBox->setChecked(meshSettings_.useSplines);
    splineBox->setToolTip(
        "On, the sides are sampled along the fitted B-splines and the interiors\n"
        "come from the bicubic patch, so every node lies on the reconstructed\n"
        "surface. Off, the traced polylines are used instead and the interiors\n"
        "are a discrete Coons blend of them — the layout with no fit trusted\n"
        "anywhere, which is what tells a meshing artefact from a fitting one.");

    auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dlg);
    buttons->button(QDialogButtonBox::Ok)->setText("Mesh");
    QObject::connect(buttons, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
    QObject::connect(buttons, &QDialogButtonBox::rejected, &dlg, &QDialog::reject);

    auto *form = new QFormLayout(&dlg);
    form->addRow("target edge length", targetBox);
    form->addRow("implied sizing", derived);
    form->addRow("fewest edges per chord", minBox);
    form->addRow("most edges per chord", maxBox);
    form->addRow("contract chords thinner than", collapseBox);
    form->addRow(splineBox);
    form->addRow(new QLabel("The grid is transfinite interpolation only:\n"
                            "no smoothing is run on it.", &dlg));

    // Stage 11's three settings live in Stage 10's dialog, because Stage 11 is
    // rerun with the mesh and never on its own: the fill is built on the rim
    // the interval assignment produced, so trying another squareness means
    // meshing again. Shown only where something was excised.
    QDoubleSpinBox *squareBox = nullptr;
    QSpinBox *ringBox = nullptr, *smoothBox = nullptr;
    if (!inclusions_.empty()) {
        form->addRow(new QLabel(QString("Stage 11 — the O-grid on %1 excised inclusion(s)")
                                    .arg(static_cast<int>(inclusions_.size())), &dlg));

        squareBox = new QDoubleSpinBox(&dlg);
        squareBox->setRange(0.0, 1.0);
        squareBox->setDecimals(2);
        squareBox->setSingleStep(0.05);
        squareBox->setValue(diskSettings_.squareness);
        squareBox->setToolTip(
            "How square the core block's boundary is.\n"
            "0 is a scaled copy of the rim, which puts a straight angle at each of\n"
            "the four core corners — the very defect meshing the disk as one patch\n"
            "has. 1 is the straight chords between them, which makes the ring 1.6x\n"
            "deeper at the middle of a side than at its ends.");

        ringBox = new QSpinBox(&dlg);
        ringBox->setRange(0, 64);
        ringBox->setValue(diskSettings_.ringDepth);
        ringBox->setSpecialValueText("auto");
        ringBox->setToolTip(
            "Rows of elements between the rim and the core.\n"
            "Auto takes it from the rim spacing, so a ring element is about as\n"
            "deep as it is wide.");

        smoothBox = new QSpinBox(&dlg);
        smoothBox->setRange(0, 5000);
        smoothBox->setValue(diskSettings_.smoothing);
        smoothBox->setToolTip(
            "Laplacian sweeps over the nodes the template placed, the rim held.\n"
            "A move is taken only where it does not lower the worst scaled\n"
            "Jacobian around the node, so this cannot fold an element.\n"
            "Unlike Stage 10's smoothing, which this dialog leaves off, it is\n"
            "what takes the transfinite start to the O-grid picture.");

        form->addRow("core squareness", squareBox);
        form->addRow("ring rows", ringBox);
        form->addRow("template smoothing sweeps", smoothBox);
    }

    form->addRow(buttons);

    if (dlg.exec() != QDialog::Accepted) return false;

    meshSettings_.target     = targetBox->value();
    meshSettings_.minEdges   = minBox->value();
    meshSettings_.maxEdges   = maxBox->value();
    meshSettings_.useSplines = splineBox->isChecked();
    meshSettings_.collapseSpan = collapseBox->value();
    if (squareBox) diskSettings_.squareness = squareBox->value();
    if (ringBox)   diskSettings_.ringDepth  = ringBox->value();
    if (smoothBox) diskSettings_.smoothing  = smoothBox->value();
    return true;
}

// Stage 10, at the dialog's settings and with the smoothing off.
//
// Deliberately the raw transfinite grid: Winslow moves interior nodes off the
// spacing the interval assignment chose for them, which is a good trade on a
// folded block and a bad one everywhere else, and it is the assignment that
// this phase is a picture of. What the smoothing buys is measured in
// TestMERIDIAN, where the before and the after can be put side by side.
void CrossGenWidget::runMERIDIANMesh() {
    if (!splines_.has_value()) {           // Stage 9 not there yet: not an attempt
        blockPipeline("Stage 10 not built", "Stage 9 produced no patches to mesh");
        return;
    }
    meshAttempted_ = true;
    quadMesh_.reset();
    diskFill_.reset();
    // Stage 12 stood on the mesh that is about to be replaced, so it goes with
    // it rather than being drawn over the next one.
    smoothMesh_.reset();
    tmopAttempted_ = false;
    pipelineBlocked_.clear();

    QuadMesh::Options qopts;
    qopts.targetEdgeLength = meshSettings_.target;
    qopts.minIntervals     = meshSettings_.minEdges;
    qopts.maxIntervals     = meshSettings_.maxEdges;
    qopts.useSplines       = meshSettings_.useSplines;
    qopts.collapseSpan     = meshSettings_.collapseSpan;
    // Stage 10's own selective smoothing, at the pipeline's setting. It is not
    // Stage 12: this is the Winslow pass QuadMesh runs on the nodes it is
    // allowed to move, and both pipelines and the paper_tests layout runs (E4,
    // E5 Part 3) read their scaled Jacobian off the mesh *after* it. Zeroing it
    // here would mean the element quality reported below, and the mesh Stage 12
    // is then handed, are not the ones the experiments measure.
    qopts.smoothingPasses    = TORSION::Options().quadSmoothingPasses;
    qopts.smoothingThreshold = TORSION::Options().quadSmoothingThreshold;
    qopts.featuresOnTracedArcs = TORSION::Options().quadFeaturesOnTracedArcs;

    // Which arcs lie on each excised rim. Wanted here rather than at Stage 11
    // because the interval assignment has to know: a rim with an odd number of
    // edges cannot be quadrangulated by any template at all, so the counts are
    // chosen with that parity constraint in them. Empty on a model where
    // nothing was excised, and then this is Stage 10 exactly as it was.
    std::vector<std::vector<int>> rims;
    if (!inclusions_.empty() && arrangement_.has_value()) {
        rims = DiskTemplate::rimArcs(*arrangement_, inclusions_);
        for (const std::vector<int> &r : rims)
            if (!r.empty()) qopts.evenLoops.push_back(r);
    }

    auto t0 = Clock::now();
    try {
        quadMesh_.emplace(*splines_, qopts);
    } catch (const std::exception &e) {
        quadMesh_.reset();
        blockPipeline("Stage 10 failed", e.what());
        return;
    }
    auto t1 = Clock::now();

    const QuadMesh::Report &qr = quadMesh_->getReport();
    {
        std::ostringstream oss;
        oss << "[Mesh] Stage 10: " << qr.quads << " quad(s) on " << qr.vertices
            << " vertices over " << qr.blocks << " block(s), "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }
    {
        // The assignment, and what it bought, read back off the finished mesh
        // rather than asserted: the chords are what was chosen, the edge
        // lengths are what came out.
        std::ostringstream oss;
        oss << "[Mesh] " << qr.chords << " chord(s) over " << qr.arcsAssigned
            << " arc(s), " << qr.minIntervals << " to " << qr.maxIntervals
            << " edges each (mean " << std::fixed << std::setprecision(2)
            << qr.meanIntervals << ")";
        if (qr.clampedChords > 0) oss << ", " << qr.clampedChords << " clamped by a bound";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Mesh] edges " << std::fixed << std::setprecision(4) << qr.minEdge
            << " to " << qr.maxEdge << " against a target of " << qr.target
            << " (worst " << std::setprecision(2) << qr.worstEdgeRatio
            << "x, rms log ratio " << std::setprecision(3) << qr.edgeRatioRms << ")";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Mesh] scaled Jacobian " << std::fixed << std::setprecision(4)
            << qr.minScaledJacobian << " worst, " << qr.meanScaledJacobian << " mean";
        if (qr.invertedQuads > 0) {
            oss << " -- " << qr.invertedQuads << " element(s) fold";
            if (qr.reflexCorners > 0)
                oss << " at " << qr.reflexCorners << " block corner(s) whose angle on the "
                    << "model is more than pi, which no grid can cover unreversed";
        }
        console_.log(oss.str());
    }
    if (qr.unmeshedPatches > 0) {
        std::ostringstream oss;
        oss << "[Mesh] " << qr.unmeshedPatches
            << " patch(es) not meshed -- Stage 8 did not close them as quadrilaterals of "
            << "one arc a side; they are " << std::fixed << std::setprecision(4)
            << qr.unmeshedArea << " of S and are left blank";
        console_.log(oss.str());
    }
    if (qr.collapsedChords > 0) {
        std::ostringstream oss;
        oss << "[Mesh] " << qr.collapsedChords << " chord(s) contracted, merging "
            << qr.collapsedPatches << " face(s) -- " << std::fixed << std::setprecision(4)
            << qr.collapsedArea << " of S -- into their neighbours on " << qr.weldedVertices
            << " welded vertex/vertices. Those faces were thinner than the elements asked "
            << "for, so no interval assignment could have covered them at the target size; "
            << "the dialog's \"contract chords thinner than\" is what to turn off to see it";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Mesh] " << qr.interiorEdges << " interior and " << qr.boundaryEdges
            << " boundary edge(s), " << qr.nonManifoldEdges << " used a third time, "
            << qr.cracks << " crack(s): "
            << (qr.conforming ? "conforming [PASS]" : "not conforming [FAIL]");
        console_.log(oss.str());
    }
    if (qr.materials > 1) {
        std::ostringstream oss;
        oss << "[Mesh] materials: " << qr.materials << " over " << qr.interfaceEdges
            << " element edge(s) on an interface, " << qr.mixedQuads
            << " element(s) straddling one "
            << (qr.materialsPure ? "[PASS]" : "[FAIL]");
        if (qr.unlocatedQuads > 0)
            oss << " (" << qr.unlocatedQuads << " element(s) located off the triangulation)";
        console_.log(oss.str());
    }
    for (const std::string &m : qr.messages) console_.log("[Mesh] " + m);
    {
        std::ostringstream oss;
        oss << "[Mesh] Stage 10 " << (qr.valid ? "validates [PASS]"
                                               : "does not validate [FAIL]");
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }
    if (!qopts.evenLoops.empty()) {
        std::ostringstream oss;
        oss << "[Mesh] rim parity: " << qr.oddLoops << " of " << qopts.evenLoops.size()
            << " rim(s) came out odd, " << qr.parityChordsMoved
            << " chord(s) moved by one edge to fix it, " << qr.oddLoopsLeft << " left odd";
        console_.log(oss.str());
    }
    console_.log("[Mesh] transfinite grid, unsmoothed: grey is a mesh edge, blue a block "
                 "wall, red a folded element. Press 'c' to mesh again at another target");

    runDiskTemplates(rims);
}

// ── Stage 12: TMOP ───────────────────────────────────────────────────────────
//
// The dialog. Three of its rows are the decision and the rest are safety
// catches, so they are ordered that way: the metric, the exponent and the sweep
// cap first, then which nodes are allowed to move, then the target and the
// quadrature, then the two switches -- untangling and the thread count -- that
// only ever answer "is this a defect of the smoother or of the mesh?".
//
// The line under the metric box is the reason this is a dialog and not a
// constant. Every barrier metric is infinite on an inverted element, and Stage
// 10 delivers folds on exactly the models where they are hardest to look at, so
// what a given metric can even be started on is a fact about *this* mesh. The
// counts are read off the mesh that is about to be smoothed and not off the
// last one.
// The TMOP settings the current mode smooths with. ATLAS, UMBER and ZIPLINE
// keep their own (atlasTmopSettings_ and the two after it), which
// differ from the pipelines' in one default for one reason: each smooths a
// transfinite grid on a block decomposition. The dialog and the solve both read
// them through here, so what the dialog sets is what the solve runs -- in UMBER
// mode the dialog used to edit the pipelines' copy while the solve read UMBER's.
CrossGenWidget::TMOPSettings &CrossGenWidget::tmopSettingsForMode() {
    if (mode_ == Mode::ATLAS) return atlasTmopSettings_;
    if (mode_ == Mode::UMBER) return umberTmopSettings_;
    if (mode_ == Mode::ZIPLINE) return ziplineTmopSettings_;
    return tmopSettings_;
}

bool CrossGenWidget::promptTMOP() {
    TMOPSettings &ts = tmopSettingsForMode();
    if (!haveFinishedMesh()) return false;
    const bool atlas = (mode_ == Mode::ATLAS);
    const bool blockMode = atlas || mode_ == Mode::UMBER || mode_ == Mode::ZIPLINE;

    QDialog dlg(this);
    dlg.setWindowTitle(atlas ? "ATLAS — TMOP smoothing"
                             : (mode_ == Mode::UMBER ? "UMBER — TMOP smoothing"
                                                     : (mode_ == Mode::ZIPLINE ? "ZIPLINE — TMOP smoothing"
                                                                               : "Stage 12 — TMOP smoothing")));

    auto *metricBox = new QComboBox(&dlg);
    // Ordered as the header lists them, and paired with the enum value rather
    // than with the row index so that reordering this list cannot silently
    // change what the viewer runs.
    const std::vector<std::pair<int, const char *>> metrics = {
        { mesh::TMOP::Shape002,       "002  shape, barrier  (|T|^2 / 2tau - 1)" },
        { mesh::TMOP::Shape004,       "004  shape, no barrier  (|T|^2 - 2 tau)" },
        { mesh::TMOP::ShapeSize007,   "007  shape and size, barrier" },
        { mesh::TMOP::Size055,        "055  size only, no barrier" },
        { mesh::TMOP::Size056,        "056  size only, barrier" },
        { mesh::TMOP::ShapeSizeCombo, "100  (1-g) 002 + g 056" },
        { mesh::TMOP::Untangle022,    "022  untangler, on its own" },
    };
    for (const auto &m : metrics)
        metricBox->addItem(QString::fromUtf8(m.second), m.first);
    {
        const int at = metricBox->findData(ts.metric);
        metricBox->setCurrentIndex(at >= 0 ? at : 0);
    }
    metricBox->setToolTip(
        "What mu measures. 002 is pure shape: it asks every element to be square\n"
        "and says nothing about how big it is, which is what a mesh coming out of\n"
        "transfinite interpolation at an already-chosen element count wants.\n"
        "A barrier metric is +inf on a folded element, so a mesh with one has to\n"
        "be untangled before it can be started on — see the switch below.");

    auto *gammaBox = new QDoubleSpinBox(&dlg);
    gammaBox->setRange(0.0, 1.0);
    gammaBox->setDecimals(2);
    gammaBox->setSingleStep(0.05);
    gammaBox->setValue(ts.gamma);
    gammaBox->setToolTip("Weight on the size term of metric 100. Ignored by every other metric.");

    auto *powerBox = new QDoubleSpinBox(&dlg);
    powerBox->setRange(1.0, 8.0);
    powerBox->setDecimals(1);
    powerBox->setSingleStep(1.0);
    powerBox->setValue(ts.exponent);
    powerBox->setToolTip(
        "Minimise the integral of mu^p rather than of mu.\n"
        "1 optimises the average, and an average can be lowered by making most\n"
        "elements better and the worst one worse — which is the wrong trade for a\n"
        "mesh, whose usable timestep is set by its worst element. 2 is the default\n"
        "and measurably the better of the two on every model in data/meshes.\n"
        "Values between 1 and 2 are raised to 2: there mu^p has a singular second\n"
        "derivative exactly where a good element sits.");

    auto *sweepBox = new QSpinBox(&dlg);
    sweepBox->setRange(1, 100000);
    sweepBox->setValue(ts.sweeps);
    sweepBox->setToolTip(
        "Most Gauss-Seidel sweeps over the movable nodes. The solve stops before\n"
        "this when a sweep stops moving anything, and says which of the two it was.");

    auto *pinBox = new QCheckBox("pin every feature node instead of sliding it", &dlg);
    pinBox->setChecked(ts.pinFeatures);
    pinBox->setToolTip(
        "On, only interior nodes move and the boundary and the material interfaces\n"
        "stay exactly where Stage 10 put them. Off, a node on a smooth stretch of\n"
        "either slides along the chord through its two feature neighbours, which\n"
        "redistributes the boundary without moving the domain: that chord is\n"
        "parallel to the base of the only triangle the enclosed area depends on\n"
        "that node through, so the area is conserved to rounding either way.\n"
        "Pinning is the isolation knob — it is what tells an interior defect from\n"
        "a boundary one.");

    auto *cornerBox = new QDoubleSpinBox(&dlg);
    cornerBox->setRange(0.0, 90.0);
    cornerBox->setDecimals(1);
    cornerBox->setSingleStep(5.0);
    cornerBox->setValue(ts.cornerAngle);
    cornerBox->setToolTip(
        "A feature node whose two incident feature edges meet at less than\n"
        "(180 - this) degrees is a corner and is fixed; anything straighter slides.\n"
        "Generous on purpose: a node on a discretised circle turns by 360/nSeg,\n"
        "and none of those may be mistaken for corners.");

    auto *curveBox = new QComboBox(&dlg);
    curveBox->addItem("the interpolant of its smooth run", mesh::QuadMesh::Options::CurveSpline);
    curveBox->addItem("the run's own polyline, exactly", mesh::QuadMesh::Options::CurvePolyline);
    curveBox->addItem("the chord through its two neighbours", mesh::QuadMesh::Options::CurveChord);
    {
        const int at = curveBox->findData(ts.curveSource);
        curveBox->setCurrentIndex(at >= 0 ? at : 0);
    }
    curveBox->setToolTip(
        "What a sliding feature node moves along, and so where it is after a step.\n"
        "The chord leaves it one sagitta inside the feature every time, always on\n"
        "the inside of a turn, so a discretised circle relaxes towards its own\n"
        "inscribed polygon. The other two make the node's degree of freedom the\n"
        "parameter of a curve it cannot leave: the polyline itself, which changes\n"
        "no geometry, or the C2 interpolant through the run's nodes, which rounds\n"
        "that circle back out towards the circle.");

    auto *targetBox = new QComboBox(&dlg);
    targetBox->addItem("uniform square, side h", mesh::TMOP::TargetUniformSquare);
    targetBox->addItem("each element's own current shape", mesh::TMOP::TargetCurrentShape);
    {
        const int at = targetBox->findData(ts.target);
        targetBox->setCurrentIndex(at >= 0 ? at : 0);
    }
    targetBox->setToolTip(
        "Where W comes from. The uniform square is the ordinary choice for shape\n"
        "optimisation. The current shape starts every element at T = I, so the\n"
        "solve removes only the variation of shape *within* an element and cannot\n"
        "introduce a size mismatch of its own — which is what to use when the\n"
        "question is whether the node mobility is behaving.");

    auto *sizeBox = new QDoubleSpinBox(&dlg);
    sizeBox->setRange(0.0, 10.0);
    sizeBox->setDecimals(4);
    sizeBox->setSingleStep(0.01);
    sizeBox->setValue(ts.targetSize);
    sizeBox->setSpecialValueText("mean edge");
    sizeBox->setToolTip("h for the uniform-square target. Blank takes the mesh's own mean edge length.");

    auto *cornerQuadBox = new QCheckBox("sample mu at the corners instead of 2x2 Gauss", &dlg);
    cornerQuadBox->setChecked(ts.corners);
    cornerQuadBox->setToolTip(
        "The corner Jacobians are the ones the scaled-Jacobian report reads, so\n"
        "this drives exactly the number the result is judged by. 2x2 Gauss is the\n"
        "default because it sees the interior of the element and not only its rim.");

    auto *untangleBox = new QCheckBox("untangle first where an element is folded", &dlg);
    untangleBox->setChecked(ts.untangle);
    untangleBox->setToolTip(
        "Run metric 022 first, whose barrier sits below the worst determinant on\n"
        "the mesh rather than at zero, until nothing is inverted. Without it a\n"
        "barrier metric cannot be started on a mesh with a fold at all — which is\n"
        "what turning this off is for: it shows which folds Stage 10 left.");

    auto *threadBox = new QSpinBox(&dlg);
    threadBox->setRange(0, 256);
    threadBox->setValue(ts.threads);
    threadBox->setSpecialValueText("auto");
    threadBox->setToolTip(
        "OpenMP threads. The mesh that comes out is bit-identical whatever this\n"
        "is — the sweep is coloured, so two nodes updated together never share an\n"
        "element — so this changes how long the phase takes and nothing else.");

    // What the settings come to on this mesh. Recomputed as the two mobility
    // boxes move, because those are the ones that change the answer: the metric
    // and the exponent decide what is minimised, but the node types decide how
    // much of the mesh is even allowed to move.
    auto *derived = new QLabel(&dlg);
    derived->setTextFormat(Qt::PlainText);
    auto updateDerived = [&, derived, pinBox, cornerBox, curveBox]() {
        mesh::QuadMesh::Options o;
        o.cornerAngle = cornerBox->value();
        o.fixAllFeatureNodes = pinBox->isChecked();
        o.curveSource = static_cast<mesh::QuadMesh::Options::CurveSource>(
            curveBox->currentData().toInt());
        // This runs inside a Qt signal, where an escaping exception is a
        // terminate rather than a message. A mesh the class refuses to build is
        // a real answer to "what would the smoother see here", so it is shown
        // in the label the counts would have gone in.
        mesh::QuadMesh m;
        try {
            m = finishedMesh(o);
        } catch (const std::exception &e) {
            derived->setText(QString("cannot be built: %1").arg(e.what()));
            return;
        }
        const mesh::QuadMesh::Quality &q = m.quality;
        std::ostringstream oss;
        oss << q.quads << " element(s) on " << q.vertices << " node(s)\n"
            << q.freeNodes << " free, " << q.slidingNodes << " sliding, "
            << q.fixedNodes << " fixed\n"
            << "scaled Jacobian " << std::fixed << std::setprecision(4)
            << q.minScaledJacobian << " worst, " << q.meanScaledJacobian << " mean";
        if (q.invertedQuads > 0)
            oss << "\n" << q.invertedQuads << " element(s) folded: a barrier metric "
                << "needs the untangler";
        derived->setText(QString::fromStdString(oss.str()));
    };
    QObject::connect(pinBox, &QCheckBox::toggled, &dlg, updateDerived);
    QObject::connect(cornerBox, &QDoubleSpinBox::valueChanged, &dlg, updateDerived);
    QObject::connect(curveBox, &QComboBox::currentIndexChanged, &dlg, updateDerived);
    updateDerived();

    auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dlg);
    buttons->button(QDialogButtonBox::Ok)->setText("Smooth");
    QObject::connect(buttons, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
    QObject::connect(buttons, &QDialogButtonBox::rejected, &dlg, &QDialog::reject);

    auto *form = new QFormLayout(&dlg);
    form->addRow("metric", metricBox);
    form->addRow("size weight gamma (metric 100)", gammaBox);
    form->addRow("exponent p", powerBox);
    form->addRow("sweeps at most", sweepBox);
    form->addRow(pinBox);
    form->addRow("corner if the feature turns by", cornerBox);
    form->addRow("a feature node slides on", curveBox);
    form->addRow("target", targetBox);
    form->addRow("target size h", sizeBox);
    form->addRow(cornerQuadBox);
    form->addRow(untangleBox);
    form->addRow("threads", threadBox);
    form->addRow("this mesh", derived);
    form->addRow(new QLabel(blockMode ? "Always run on the TFI mesh, never on the last\n"
                                    "smoothed one, so a second setting is a fresh attempt."
                                  : "Always run on the Stage 10 mesh, never on the last\n"
                                    "smoothed one, so a second setting is a fresh attempt.", &dlg));
    form->addRow(buttons);

    if (dlg.exec() != QDialog::Accepted) return false;

    ts.metric      = metricBox->currentData().toInt();
    ts.gamma       = gammaBox->value();
    ts.exponent    = powerBox->value();
    ts.sweeps      = sweepBox->value();
    ts.pinFeatures = pinBox->isChecked();
    ts.cornerAngle = cornerBox->value();
    ts.curveSource = curveBox->currentData().toInt();
    ts.target      = targetBox->currentData().toInt();
    ts.targetSize  = sizeBox->value();
    ts.corners     = cornerQuadBox->isChecked();
    ts.untangle    = untangleBox->isChecked();
    ts.threads     = threadBox->value();
    return true;
}

// Stage 12 itself. Blocking, like the three solves before it, but on a
// different scale: a thousand sweeps over a mesh of this size is a fraction of
// a second, so it is not announced a frame ahead the way the Ricci solve is.
void CrossGenWidget::runTMOP() {
    // ATLAS, UMBER and ZIPLINE keep their own settings; see
    // tmopSettingsForMode().
    const bool atlas = (mode_ == Mode::ATLAS);
    // UMBER and ZIPLINE: the two that smooth a BlockQuadMesh.
    const bool blockMesh = (mode_ == Mode::UMBER) || (mode_ == Mode::ZIPLINE);
    TMOPSettings &ts = tmopSettingsForMode();
    // The pipelines number this Stage 12 after their Stage 10; ATLAS's stages
    // stop at 6, UMBER's at its decomposition and ZIPLINE's at Sec. 4, so
    // there it is named for what it smooths.
    const std::string stage = (atlas || blockMesh) ? "TMOP" : "Stage 12";
    if (!haveFinishedMesh()) {
        blockPipeline(stage + " not run",
                      (atlas || blockMesh) ? "there is no TFI mesh to smooth"
                                       : "Stage 10 produced no mesh to smooth");
        return;
    }
    tmopAttempted_ = true;
    smoothMesh_.reset();
    pipelineBlocked_.clear();

    // Stage 11's merged mesh where there is one, Stage 10's otherwise. The two
    // are the same mesh to this code -- the templates' elements are elements of
    // it -- and DiskTemplate's vertex array extends Stage 10's rather than
    // replacing it, so the blocks the drawing takes its walls from still index
    // what comes back.
    mesh::QuadMesh::Options mopts;
    mopts.cornerAngle = ts.cornerAngle;
    mopts.fixAllFeatureNodes = ts.pinFeatures;
    mopts.curveSource = static_cast<mesh::QuadMesh::Options::CurveSource>(ts.curveSource);

    auto t0 = Clock::now();
    try {
        smoothMesh_.emplace(finishedMesh(mopts));
    } catch (const std::exception &e) {
        smoothMesh_.reset();
        blockPipeline(stage + " failed", e.what());
        return;
    }

    mesh::TMOP::Options topt;
    topt.metric     = static_cast<mesh::TMOP::Metric>(ts.metric);
    topt.gamma      = ts.gamma;
    topt.exponent   = ts.exponent;
    topt.target     = static_cast<mesh::TMOP::Target>(ts.target);
    topt.targetSize = ts.targetSize;
    topt.quadrature = ts.corners ? mesh::TMOP::Corners : mesh::TMOP::Gauss2x2;
    topt.maxSweeps  = ts.sweeps;
    topt.untangle   = ts.untangle;
    topt.threads    = ts.threads;

    // Nothing below the solve is allowed to take the window down with it: a
    // stage that refuses says so on the overlay and leaves the phase before it
    // on screen, which is the discipline every other stage here follows.
    mesh::TMOP smoother(*smoothMesh_, topt);
    bool ok = false;
    try {
        ok = smoother.run();
    } catch (const std::exception &e) {
        smoothMesh_.reset();
        blockPipeline(stage + " failed", e.what());
        return;
    }
    auto t1 = Clock::now();
    const mesh::TMOP::Report &tr = smoother.getReport();

    {
        std::ostringstream oss;
        oss << "[TMOP] " << ((atlas || blockMesh) ? "" : "Stage 12: ") << tr.sweeps << " sweep(s)";
        if (tr.untangleSweeps > 0) oss << " after " << tr.untangleSweeps << " untangling";
        oss << " over " << tr.colors << " colour(s) on " << tr.threads << " thread(s)"
            << (tr.openMP ? "" : " (no OpenMP)") << ", "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count())
            << (tr.converged ? "" : " -- stopped on the sweep cap");
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[TMOP] nodes: " << tr.freeNodes << " free, " << tr.slidingNodes
            << " sliding, " << tr.fixedNodes << " fixed";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        if (tr.featureCurves == 0) {
            oss << "[TMOP] feature nodes slide along the chord through their two neighbours";
        } else {
            oss << "[TMOP] " << tr.curveNodes << " feature node(s) on " << tr.featureCurves
                << " curve(s), " << tr.fittedCurves << " of them interpolants; bow "
                << std::scientific << std::setprecision(2) << tr.curveBow << ", off-curve "
                << tr.curveDeviation;
        }
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[TMOP] scaled Jacobian " << std::fixed << std::setprecision(4)
            << tr.minScaledJacobianBefore << " -> " << tr.minScaledJacobianAfter
            << " worst, " << tr.meanScaledJacobianBefore << " -> "
            << tr.meanScaledJacobianAfter << " mean";
        if (tr.invertedBefore > 0 || tr.invertedAfter > 0)
            oss << ", folds " << tr.invertedBefore << " -> " << tr.invertedAfter;
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[TMOP] worst aspect ratio " << std::fixed << std::setprecision(2)
            << tr.worstAspectBefore << " -> " << tr.worstAspectAfter
            << ", energy " << std::scientific << std::setprecision(3)
            << tr.energyBefore << " -> " << tr.energyAfter;
        console_.log(oss.str());
    }
    {
        // Sliding moves a node along the chord through its two feature
        // neighbours, which leaves the triangle they span -- and so the
        // polygon's area -- unchanged. A difference here is the domain itself
        // moving, which no amount of smoothing is worth, so it is reported
        // whether it happened or not rather than only when it did.
        const double drift = std::fabs(tr.areaAfter - tr.areaBefore);
        const bool held = drift <= 1e-9 * std::fabs(tr.areaBefore);
        std::ostringstream oss;
        oss << "[TMOP] nodes moved " << std::fixed << std::setprecision(3)
            << tr.meanDisplacement << " mean edge(s) on average, " << tr.maxDisplacement
            << " at most; area " << (held ? "unchanged [PASS]" : "changed [FAIL]");
        console_.log(oss.str());
    }
    for (const std::string &m : tr.messages) console_.log("[TMOP] " + m);
    {
        std::ostringstream oss;
        oss << "[TMOP] " << ((atlas || blockMesh) ? "the smoother " : "Stage 12 ")
            << (ok ? "left the mesh better than it found it [PASS]" : "did not improve the mesh [FAIL]");
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }
    console_.log((atlas || blockMesh)
                     ? "[TMOP] the same picture as the TFI mesh, with the nodes where the solve "
                       "left them. Press 'c' to smooth the TFI mesh again at other settings"
                     : "[TMOP] the same picture as Stage 10, with the nodes where the solve "
                       "left them. Press 'c' to smooth the Stage 10 mesh again at other settings");
}

// The mesh Stage 12 is handed, whichever stage finished it. In ATLAS mode it
// is ATLAS's own and nothing else, so a MERIDIAN mesh left over from before a
// reset could not be smoothed in its place.
bool CrossGenWidget::haveFinishedMesh() const {
    if (mode_ == Mode::ATLAS) return atlasMesh_.has_value();
    if (mode_ == Mode::UMBER) return umberMesh_.has_value();
    if (mode_ == Mode::ZIPLINE) return ziplineMesh();
    return quadMesh_.has_value();
}

mesh::QuadMesh CrossGenWidget::finishedMesh(const mesh::QuadMesh::Options &o) const {
    if (mode_ == Mode::ATLAS) return mesh::QuadMesh::from(*atlasMesh_, o);
    if (mode_ == Mode::UMBER) return mesh::QuadMesh::from(*umberMesh_, o);
    if (mode_ == Mode::ZIPLINE) return mesh::QuadMesh::from(*ziplineMesh(), o);
    return diskFill_.has_value() ? mesh::QuadMesh::from(*diskFill_, o)
                                 : mesh::QuadMesh::from(*quadMesh_, o);
}

// ── ATLAS ─────────────────────────────────────────────────────────────────────
//
// docs/square_transport_2d_theory_and_implementation.md, Stages 1 to 6, then the
// pipelines' last two phases on its blocks. Every stage is ATLAS's library code
// with its own defaults; the viewer only decides when each runs and what of it
// is drawn, as it does for the pipelines.

// Stage 1 on its own: validate and tag the domain.
void CrossGenWidget::runATLASDomain() {
    atlasDomainAttempted_ = true;
    pipelineBlocked_.clear();
    auto t0 = Clock::now();
    try {
        atlasDomain_ = std::make_unique<PlanarDomain>(*mesh_, ATLAS::Options().domain);
    } catch (const std::exception &e) {
        atlasDomain_.reset();
        blockPipeline("Stage 1 failed", e.what());
        return;
    }
    auto t1 = Clock::now();
    const PlanarDomain::Report &r = atlasDomain_->getReport();
    {
        std::ostringstream oss;
        oss << "[Domain] Stage 1: " << r.vertices << " vertices, " << r.triangles << " triangles; "
            << r.components << " component(s), " << r.loops << " boundary loop(s), " << r.holes
            << " hole(s), " << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Domain] protected: " << r.boundaryCorners << " boundary corner(s)";
        if (r.materials > 1) {
            oss << ", " << r.interfaceJunctions << " junction(s), " << r.interfaceLandings << " landing(s), "
                << r.interfaceKinks << " kink(s) on " << r.interfaceEdges << " interface edge(s) between "
                << r.materials << " materials";
        }
        console_.log(oss.str());
    }
    for (const std::string &m : r.messages) console_.log("[Domain] " + m);
    {
        std::ostringstream oss;
        oss << "[Domain] Stage 1 " << (r.valid ? "validates: a planar domain as Sec. 1.1 asks [PASS]"
                                              : "does not validate [FAIL]");
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }
    if (!r.valid) {
        blockPipeline("Stage 1 failed", "the mesh is not a domain Sec. 1.1 accepts");
        return;
    }
    console_.log("[Domain] a disk at every protected vertex, coloured by the quarter turns its corner "
                 "takes out of Sec. 4's identity -- a cone's index, in the cones' colours");
}

// Stage 1b on its own: the reference cross field
// (docs/atlas_crossfield_guidance.md, Sec. 3). The solve is either pipeline's
// Stage 0 -- DualMBO with the orthogonal penalty, gamma = 10 and the
// tau-continuation -- on the input mesh with the domain's interfaces as
// aligned edges, and it is the same solve ATLAS::run() will do for itself at
// this stage. Nothing downstream in the viewer reads what is built here: it is
// built so that the phase has something to draw, and so that the console says
// what the search is about to be asked for before it is asked for it.
void CrossGenWidget::runATLASField() {
    atlasFieldAttempted_ = true;
    atlasField_.reset();
    if (!atlasDomain_ || !atlasDomain_->getReport().valid) {
        blockPipeline("Stage 1b not built", "Stage 1 did not validate the domain");
        return;
    }
    const ATLAS::Options o;
    if (o.field.reference != ATLAS::Options::Field::Reference::DualMBO) {
        console_.log("[Field] Stage 1b: ATLAS is configured for the harmonic reference, so there is "
                     "no DualMBO field to show; the searches score their layouts against "
                     "CoarseDomain::crossAngle");
        return;
    }
    pipelineBlocked_.clear();

    std::vector<int> interfaceEdges;
    for (int e = 0; e < static_cast<int>(atlasDomain_->interfaceEdge.size()); ++e) {
        if (atlasDomain_->interfaceEdge[e]) interfaceEdges.push_back(e);
    }
    try {
        atlasField_ = std::make_unique<ReferenceField>(mesh_, interfaceEdges, o.field.solve);
    } catch (const std::exception &e) {
        atlasField_.reset();
        console_.log(std::string("[Field] Stage 1b failed: ") + e.what() +
                     "; the searches would run without a field");
        return;
    }

    const ReferenceField::Report &r = atlasField_->getReport();
    for (const std::string &m : r.messages) console_.log("[Field] " + m);
    if (!r.built) {
        std::ostringstream oss;
        oss << "[Field] Stage 1b: no reference field (" << r.reason
            << "); the searches would run without one";
        console_.log(oss.str());
        atlasField_.reset();
        return;
    }
    {
        std::ostringstream oss;
        oss << "[Field] Stage 1b: DualMBO on " << r.triangles << " triangles, " << r.levels
            << " tau level(s), " << r.steps << " step(s), " << (r.converged ? "converged" : "not converged")
            << ", " << formatMs(1000.0 * r.seconds);
        console_.log(oss.str());
    }
    {
        // The dipole count is worth its line: a +1/-1 pair the smoothest field
        // happens to carry leaves every region's index where it was, and a
        // carrier asked to reproduce it would have to carry a dislocation --
        // the defect v1 of ATLAS stalled on (CavityRewrite's header).
        std::ostringstream oss;
        oss << "[Field] cones: +" << r.rawPlus << " / -" << r.rawMinus << " found";
        if (r.dipoleUnits > 0)
            oss << ", " << r.dipoleUnits << " cancelled as dipoles -> +" << r.conesPlus << " / -" << r.conesMinus;
        console_.log(oss.str());
    }
    {
        // The scale the two terms are read on, said once here rather than at
        // every search: E_dir is area-normalised so that the coarse carriers,
        // their realisations and the fine one are comparable, and w_dir turns
        // it into the annealer's unit, which is blocks.
        std::ostringstream oss;
        oss << "[Field] lookup grid " << r.gridNx << " x " << r.gridNy
            << "; E_dir is area-normalised (0.01 is the whole domain about 3 degrees off) and weighted "
            << "by w_dir = " << std::fixed << std::setprecision(1) << o.field.wDir << " block(s)";
        console_.log(oss.str());
    }
    console_.log("[Field] a cross per triangle, and a disk at every cone the layout is asked to "
                 "reproduce: blue where a valence-3 vertex belongs, red where a valence-5 one does. "
                 "A +1/4 that ends up on the boundary instead is a block corner on a smooth arc, "
                 "which is what E_sing's weights exist to stop. Stage 1's corners are off while "
                 "this is drawn: a convex corner is a blue disk too, and the two would not be "
                 "tellable apart");
    std::cerr << "[Viewer] [Field] Stage 1b: +" << r.conesPlus << " / -" << r.conesMinus << " cone(s)\n";
}

// Stage 2 on its own: the three-quad split, which is the guaranteed incumbent.
void CrossGenWidget::runATLASCarrier() {
    atlasCarrierAttempted_ = true;
    if (!atlasDomain_ || !atlasDomain_->getReport().valid) {
        blockPipeline("Stage 2 not built", "Stage 1 did not validate the domain");
        return;
    }
    pipelineBlocked_.clear();
    auto t0 = Clock::now();
    atlasCarrier_ = std::make_unique<SquareCarrier>(*atlasDomain_, ATLAS::Options().carrier);
    const SquareCarrier::Report &r = atlasCarrier_->validate();
    auto t1 = Clock::now();
    {
        std::ostringstream oss;
        oss << "[Carrier] Stage 2: " << r.cells << " cells from " << r.sourceTriangles
            << " triangles, three each, on " << r.vertices << " vertices, "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Carrier] Sec. 3.1: every corner determinant / 2A_T in [" << std::fixed << std::setprecision(4)
            << r.minCornerDetRatio << ", " << r.maxCornerDetRatio << "], " << r.jacobianBoundViolations
            << " outside (1/12, 1/4); scaled Jacobian " << r.minScaledJacobian << " worst";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Carrier] irregular: " << r.irregularInterior << " inside, " << r.irregularBoundary
            << " on dS, total valence defect " << r.totalDefect << "; Sec. 4: " << r.eulerLHS << " = "
            << r.eulerRHS << (r.eulerHolds ? " [PASS]" : " [FAIL]");
        console_.log(oss.str());
    }
    for (const std::string &m : r.messages) console_.log("[Carrier] " + m);
    {
        std::ostringstream oss;
        oss << "[Carrier] Stage 2 " << (r.valid ? "validates" : "does not validate") << ": the singleton "
            << "blocking, " << r.cells << " blocks, is the incumbent " << (r.valid ? "[PASS]" : "[FAIL]");
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }
    if (!r.valid) {
        blockPipeline("Stage 2 failed", "the three-quad carrier does not validate");
        return;
    }
    console_.log("[Carrier] a disk at every vertex of the wrong valence for where it is: the "
                 "singularities the carrier inherits from the triangulation (Sec. 8.1)");
}

// Stages 1 to 6 in one call: ATLAS::run() and all of its searches.
void CrossGenWidget::runATLAS() {
    atlasAttempted_ = true;
    atlasShowInitial_ = false;
    atlasBlocksLogged_ = false;
    atlasMesh_.reset();
    smoothMesh_.reset();
    atlas_.reset();
    if (!atlasCarrier_ || !atlasCarrier_->getReport().valid) {
        blockPipeline("Stages 3-6 not run", "Stage 2 did not validate the carrier");
        return;
    }
    pipelineBlocked_.clear();

    auto t0 = Clock::now();
    bool ok = false;
    try {
        atlas_ = std::make_unique<ATLAS>(mesh_, ATLAS::Options());
        ok = atlas_->run();
    } catch (const std::exception &e) {
        atlas_.reset();
        blockPipeline("Stages 3-6 failed", e.what());
        return;
    }
    auto t1 = Clock::now();
    const double seconds = std::chrono::duration<double>(t1 - t0).count();

    // One line per search, the winner marked: which kind of carrier the layout
    // came from is the first thing to know about it.
    const ATLAS::Status &st = atlas_->getStatus();
    for (int i = 0; i < atlas_->numSearches(); ++i) {
        const ATLAS::Search &sr = atlas_->getSearch(i);
        std::ostringstream oss;
        oss << "[Search] " << sr.name << ": ";
        if (sr.succeeded()) oss << sr.finalCover()->getReport().blocks << " block(s)";
        else oss << "no valid cover";
        if (sr.coarse && sr.coarseDomain) {
            oss << " from a coarse domain of " << sr.coarseDomain->getReport().triangles << " triangles";
            if (sr.realised) oss << (sr.realisedReport.valid ? ", realised on the input" : ", did not realise");
        }
        if (i == st.chosen) oss << "  <- chosen";
        console_.log(oss.str());
    }
    if (!ok || !atlas_->hasCover()) {
        blockPipeline("Stages 3-6 found no cover", "no search produced a valid conforming blocking");
        return;
    }
    const ATLAS::Search &s = atlas_->getChosen();
    if (s.rewrite && s.annealed) {
        const CavityRewrite::AnnealReport &ar = s.rewrite->getAnnealReport();
        std::ostringstream oss;
        oss << "[Search] Stage 6 annealed: " << ar.moves << " move(s), " << ar.accepted << " accepted; base patches "
            << ar.blocksBefore << " -> " << ar.blocksAfter;
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Search] the " << s.name << " search wins with " << st.bestBlocks << " block(s) (objective "
            << std::fixed << std::setprecision(2) << st.bestObjective << "), " << std::setprecision(1) << seconds
            << " s in all";
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }
    console_.log(s.coarse ? "[Search] drawn: its coarse carrier, the input's dS in grey, block sides in blue. "
                            "Press 'p' for where it started"
                          : "[Search] drawn: its carrier, block sides in blue. Press 'p' for where it started");
}

// The chosen cover on the input, once.
void CrossGenWidget::logATLASBlocks() {
    atlasBlocksLogged_ = true;
    const ATLAS::Search &s = atlas_->getChosen();
    const SquareCarrier &C = atlas_->getCarrier();
    const BlockCover::Report &br = atlas_->getCover().getReport();
    {
        std::ostringstream oss;
        oss << "[Blocks] " << br.blocks << " block(s) over " << C.numCells() << " cell(s) of the "
            << (s.coarse ? "realised" : "fine") << " carrier; " << br.macroVertices << " macrovertices ("
            << br.irregularMacroVertices << " irregular), " << br.macroEdges << " macro edges";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Blocks] the carrier " << (C.getReport().valid ? "passes" : "fails")
            << " Stage 2's validate() on the input; scaled Jacobian " << std::fixed << std::setprecision(4)
            << C.getReport().minScaledJacobian << " worst " << (C.getReport().valid ? "[PASS]" : "[FAIL]");
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Blocks] Sec. 13.2: every cell in one block, every side dS or one neighbour's, no T-junction, "
            << "every protected vertex a corner, Sec. 4 " << br.eulerLHS << " = " << br.eulerRHS << ": "
            << (br.valid ? "a conforming blocking [PASS]" : "not conforming [FAIL]");
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }
    console_.log("[Blocks] block sides in light blue, macrovertices in green. Press 'c' to mesh the blocks");
}

// The mesh dialog: Stage 10's, less what has no counterpart here (the chord
// contraction and the disk templates), with the chart switch where Stage 10's
// spline switch is.
bool CrossGenWidget::promptATLASMesh() {
    if (!atlas_ || !atlas_->hasCover() || !mesh_) return false;

    Point lo{ std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() };
    Point hi{ -lo[0], -lo[1] };
    for (const Point &p : mesh_->vertices) {
        lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
        hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
    }
    const double diag = std::hypot(hi[0] - lo[0], hi[1] - lo[1]);
    const double extent = (diag > 0.0) ? diag : 1.0;

    QDialog dlg(this);
    dlg.setWindowTitle("ATLAS — quadrilateral mesh on the blocks");

    auto *targetBox = new QDoubleSpinBox(&dlg);
    targetBox->setRange(1e-4, 10.0);
    targetBox->setDecimals(4);
    targetBox->setSingleStep(0.01);
    targetBox->setValue(meshSettings_.target);
    targetBox->setToolTip(
        "Target length of a mesh edge, in the units of the model.\n"
        "The same number the MERIDIAN and TORSION mesh dialogs take, and\n"
        "shared with them, so the three layouts can be meshed alike.");

    auto *derived = new QLabel(&dlg);
    derived->setTextFormat(Qt::PlainText);
    const BlockCover::Report &br = atlas_->getCover().getReport();
    const int blocks = br.blocks, edges = br.macroEdges;
    auto updateDerived = [targetBox, derived, extent, blocks, edges]() {
        const double h = targetBox->value();
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(4)
            << "diagonal of S                 = " << extent << "\n"
            << "edges across it   diag / h    = " << std::setprecision(1) << (extent / h) << "\n"
            << blocks << " block(s) and " << edges << " macro edge(s) to mesh";
        derived->setText(QString::fromStdString(oss.str()));
    };
    QObject::connect(targetBox, &QDoubleSpinBox::valueChanged, &dlg, updateDerived);
    updateDerived();

    auto *minBox = new QSpinBox(&dlg);
    minBox->setRange(1, 64);
    minBox->setValue(meshSettings_.minEdges);
    minBox->setToolTip("Fewest edges any chord may be given.");

    auto *maxBox = new QSpinBox(&dlg);
    maxBox->setRange(0, 4096);
    maxBox->setValue(meshSettings_.maxEdges);
    maxBox->setSpecialValueText("none");
    maxBox->setToolTip("Most edges any chord may be given. 0 for no ceiling.");

    auto *chartBox = new QCheckBox("place the interior nodes through the blocks' charts", &dlg);
    chartBox->setChecked(atlasUseChart_);
    chartBox->setToolTip(
        "On, the four sides' chart coordinates are blended by transfinite\n"
        "interpolation and each block's certified chart -- its carrier cells, a\n"
        "valid piecewise-bilinear map (Sec. 5.3) -- is evaluated there, as Stage 10\n"
        "evaluates its spline patches. Off, the interior is the Coons blend of the\n"
        "four sides' points, with no chart trusted inside the block.");

    auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dlg);
    buttons->button(QDialogButtonBox::Ok)->setText("Mesh");
    QObject::connect(buttons, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
    QObject::connect(buttons, &QDialogButtonBox::rejected, &dlg, &QDialog::reject);

    auto *form = new QFormLayout(&dlg);
    form->addRow("target edge length", targetBox);
    form->addRow("implied sizing", derived);
    form->addRow("fewest edges per chord", minBox);
    form->addRow("most edges per chord", maxBox);
    form->addRow(chartBox);
    form->addRow(new QLabel("Transfinite interpolation per block. Only a block that\n"
                            "folds is smoothed, by Stage 10's Winslow pass.", &dlg));
    form->addRow(buttons);

    if (dlg.exec() != QDialog::Accepted) return false;

    meshSettings_.target   = targetBox->value();
    meshSettings_.minEdges = minBox->value();
    meshSettings_.maxEdges = maxBox->value();
    atlasUseChart_         = chartBox->isChecked();
    return true;
}

// BlockMesh at the dialog's settings, reported line for line as Stage 10 is.
void CrossGenWidget::runATLASMesh() {
    if (!atlas_ || !atlas_->hasCover()) {
        blockPipeline("the mesh not built", "Stages 3-6 produced no blocks to mesh");
        return;
    }
    atlasMeshAttempted_ = true;
    atlasMesh_.reset();
    // TMOP stood on the mesh that is about to be replaced.
    smoothMesh_.reset();
    tmopAttempted_ = false;
    pipelineBlocked_.clear();

    BlockMesh::Options mo;
    mo.targetEdgeLength = meshSettings_.target;
    mo.minIntervals     = meshSettings_.minEdges;
    mo.maxIntervals     = meshSettings_.maxEdges;
    mo.useChart         = atlasUseChart_;
    // Stage 10's selective Winslow pass at the pipelines' setting, for the
    // reason runMERIDIANMesh gives: the mesh TMOP is handed, and the quality
    // reported here, are then the same kind of thing for all three methods.
    mo.smoothingPasses    = TORSION::Options().quadSmoothingPasses;
    mo.smoothingThreshold = TORSION::Options().quadSmoothingThreshold;

    auto t0 = Clock::now();
    try {
        atlasMesh_.emplace(atlas_->getCover(), mo);
    } catch (const std::exception &e) {
        atlasMesh_.reset();
        blockPipeline("the mesh failed", e.what());
        return;
    }
    auto t1 = Clock::now();

    const BlockMesh::Report &r = atlasMesh_->getReport();
    {
        std::ostringstream oss;
        oss << "[Mesh] " << r.quads << " quad(s) on " << r.vertices << " vertices over " << r.blocks
            << " block(s), " << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Mesh] " << r.chords << " chord(s) over " << r.arcsAssigned << " macro edge(s), "
            << r.minIntervals << " to " << r.maxIntervals << " edges each (mean " << std::fixed
            << std::setprecision(2) << r.meanIntervals << ")";
        if (r.clampedChords > 0) oss << ", " << r.clampedChords << " clamped by a bound";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Mesh] edges " << std::fixed << std::setprecision(4) << r.minEdge << " to " << r.maxEdge
            << " against a target of " << r.target << " (worst " << std::setprecision(2) << r.worstEdgeRatio
            << "x, rms log ratio " << std::setprecision(3) << r.edgeRatioRms << ")";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Mesh] scaled Jacobian " << std::fixed << std::setprecision(4) << r.minScaledJacobian
            << " worst, " << r.meanScaledJacobian << " mean";
        if (r.smoothedBlocks > 0) {
            oss << " (Winslow on " << r.smoothedBlocks << " folded block(s): " << r.invertedBefore << " fold(s), "
                << r.minScaledJacobianBefore << " worst before)";
        }
        if (r.invertedQuads > 0) oss << " -- " << r.invertedQuads << " element(s) fold";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Mesh] Sec. 11.3: every boundary node on dS, the edges between them within " << std::fixed
            << std::setprecision(4) << r.boundaryDeviation << " of it (" << std::setprecision(3)
            << (r.target > 0.0 ? r.boundaryDeviation / r.target : 0.0) << " of the target)";
        if (r.materials > 1) {
            oss << "; " << r.materials << " materials meeting on " << r.interfaceEdges
                << " element edge(s), within " << std::setprecision(4) << r.interfaceDeviation
                << " of the interfaces";
        }
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Mesh] " << r.interiorEdges << " interior and " << r.boundaryEdges << " boundary edge(s), "
            << r.nonManifoldEdges << " used a third time, " << r.cracks << " crack(s): "
            << (r.conforming ? "conforming [PASS]" : "not conforming [FAIL]");
        console_.log(oss.str());
    }
    for (const std::string &m : r.messages) console_.log("[Mesh] " + m);
    {
        std::ostringstream oss;
        oss << "[Mesh] the mesh " << (r.valid ? "validates [PASS]" : "does not validate [FAIL]");
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }
    console_.log("[Mesh] transfinite grid: grey is a mesh edge, blue a block wall, red a folded element. "
                 "Press 'c' to smooth it, 'e' to mesh again at another target");
}

// What ATLAS mode draws. The later a phase, the more of the model its picture
// replaces -- the blocks and the meshes are drawn alone, as the pipelines draw
// theirs -- but only when that picture exists: a stage that refused leaves the
// last picture that did, so the window is never blank.
void CrossGenWidget::renderATLAS() {
    const bool matFill = showMaterialFill_ && atlasMultiMaterial_;
    const bool haveCover = atlas_ && atlas_->hasCover();

    // The mesh, and after TMOP the same mesh drawn by the same routine off the
    // same block walls: the pipelines' Mesh and Smoothed pictures exactly.
    if (atlasPhase_ >= ATLASPhase::Mesh && atlasMesh_.has_value()) {
        if (atlasPhase_ == ATLASPhase::Smoothed && smoothMesh_.has_value())
            viewer::drawQuadMesh(*smoothMesh_, *atlasMesh_, 1.0f, 2.5f, matFill);
        else
            viewer::drawQuadMesh(*atlasMesh_, 1.0f, 2.5f, matFill);
        return;
    }

    // The blocks on the input, as the shared BlockDecomposition -- the same
    // representation and the same routine MERIDIAN's and TORSION's Patches
    // phase draws its layout with.
    if (atlasPhase_ >= ATLASPhase::Blocks && haveCover) {
        viewer::drawBlockDecomposition(atlas_->getCover().blockDecomposition(), 0.30 * avgEdge_, 3.0f);
        return;
    }

    // The winning search on the carrier it searched. A coarse carrier's
    // boundary is a polygon inscribed in the input's, so the input's own dS
    // goes under it: what Realisation will carry the layout onto.
    if (atlasPhase_ == ATLASPhase::Search && haveCover) {
        const ATLAS::Search &s = atlas_->getChosen();
        const SquareCarrier *C = atlasShowInitial_ ? s.initial.get() : s.best.get();
        if (!C) C = s.best ? s.best.get() : s.initial.get();
        if (C) {
            viewer::drawBoundaryEdges(*mesh_);
            viewer::drawSquareCarrier(*C, 1.25f, 0.0, matFill);
            if (C == s.best.get() && s.bestCover)
                viewer::drawBlockDecomposition(s.bestCover->blockDecomposition(), 0.0, 2.5f);
            viewer::drawCarrierDefects(*C, 0.3 * C->meanEdgeLength());
            return;
        }
    }

    if (atlasPhase_ >= ATLASPhase::Carrier && atlasCarrier_) {
        viewer::drawSquareCarrier(*atlasCarrier_, 1.0f, 0.3 * atlasCarrier_->meanEdgeLength(), matFill);
        return;
    }

    if (matFill) viewer::drawMaterialFill(*mesh_, 0.28f);
    viewer::drawMesh(*mesh_);
    // Stage 1b over the model, as either pipeline draws its Stage 0: the
    // crosses first and the cones over them, so a cone is not hidden by the
    // arms of the crosses around it.
    //
    // Stage 1's corners are left off while it is drawn, and that is deliberate
    // rather than tidiness: a protected convex corner and a +1/4 cone are both
    // drawn as a blue disk -- they are the same index, and the layout has to
    // put a block corner at both -- so with the two layers on top of each
    // other there is no reading which disk came from which stage. The phase
    // before this one is the corners, and nothing in this one is anything but
    // the field.
    if (atlasPhase_ == ATLASPhase::Field && atlasField_) {
        viewer::drawBoundaryEdges(*mesh_);
        viewer::drawReferenceField(*mesh_, *atlasField_, scale_);
        viewer::drawReferenceCones(*atlasField_, 0.5 * avgEdge_);
        return;
    }
    if (atlasPhase_ >= ATLASPhase::Domain && atlasDomain_) {
        viewer::drawBoundaryEdges(*mesh_);
        viewer::drawPlanarDomain(*atlasDomain_, 0.5 * avgEdge_);
    }
}

// ── Stage 0c: the circular inclusions, taken out before anything is built ────
//
// Called from the mode keys and nowhere else. Every stage of either pipeline is
// written on mesh_, so the one moment this can happen is after the mode is
// chosen and before the first of them runs -- there is no phase for it, and the
// phase it would sit in is the one the field is already built by.
//
// What it buys is in DiskTemplate's header: eight cones per inclusion that the
// penalty continuation of Sec. 3.3 does not have to place, which on
// data/geometry/multimat/bubbles is the difference between a layout and 57% of
// the model left unmeshed.
void CrossGenWidget::runDiskExcision() {
    // Once per run, like every other stage here: the mode can only be chosen
    // from Unselected, so this is belt and braces rather than a live guard.
    if (disksAttempted_ || !mesh_) return;
    disksAttempted_ = true;
    inclusions_.clear();
    diskFill_.reset();

    std::vector<std::string> msgs;
    std::vector<DiskTemplate::Inclusion> found = DiskTemplate::detect(*mesh_, msgs);
    // detect() says why it refused a component that looked like a candidate --
    // the model itself being one circle, a disk touching dS -- and those are
    // worth reading even though the answer is "nothing to do".
    for (const std::string &m : msgs) console_.log("[Disks] Stage 0c: " + m);
    if (found.empty()) return;

    if (!promptDiskExcision(found)) {
        std::ostringstream oss;
        oss << "[Disks] Stage 0c declined: the " << found.size()
            << " circular inclusion(s) stay in the layout problem, cones and all";
        console_.log(oss.str());
        return;
    }

    std::shared_ptr<Mesh> excised = DiskTemplate::excise(*mesh_, found);
    if (!excised || excised->triangles.empty()) {
        console_.log("[Disks] Stage 0c: excising the inclusions would leave no mesh "
                     "behind, so they are laid out like any other region");
        return;
    }

    int excisedTris = 0;
    for (const DiskTemplate::Inclusion &inc : found)
        excisedTris += static_cast<int>(inc.triangles.size());

    inputMesh_  = mesh_;
    mesh_       = excised;
    inclusions_ = std::move(found);

    // Anything already standing on the mesh that was just replaced has to go
    // with it. Nothing should be, since this runs before the mode is set and
    // therefore before runComputations() will build any of it -- but "should
    // be" is exactly what a nested event loop breaks, and every stage below
    // checks that its inputs were measured on the mesh it was handed and
    // throws when they were not. Cheap, and it makes the ordering a property
    // of this function rather than of its one caller.
    smoothMesh_.reset();
    quadMesh_.reset();
    splines_.reset();
    arrangement_.reset();
    separatrices_.reset();
    meridianLayout_.reset();
    meridianLabels_.reset();
    immersion_.reset();
    integration_.reset();
    tutte_.reset();
    scaffold_.reset();
    frames_.reset();
    coneMetric_.reset();
    ricci_.reset();
    coneCut_.reset();
    cones_.reset();
    interfaces_.reset();
    dualMBOField_.reset();
    flatMetric_ = viewer::FlatMetric{};
    coneFans_.clear();
    ricciU_.resize(0);
    psiR_.clear();
    integratedMap_.clear();
    fieldIndex_.clear();
    interfaceMessagesSeen_ = 0;
    interfacesAttempted_ = conesAttempted_ = cutAttempted_ = false;
    ricciAttempted_ = ricciAnnounced_ = false;
    framesAttempted_ = integrationAttempted_ = integrationAnnounced_ = false;
    layoutAttempted_ = layoutAnnounced_ = false;
    separatricesAttempted_ = separatricesAnnounced_ = false;
    patchesAttempted_ = patchesAnnounced_ = false;
    meshAttempted_ = tmopAttempted_ = false;
    dualMBOSteppingStarted_ = dualMBOConverged_ = false;
    dualMBOStepCount_ = 0;

    {
        std::ostringstream oss;
        oss << "[Disks] Stage 0c: excised " << inclusions_.size() << " circular "
            << "inclusion(s), " << excisedTris << " of " << inputMesh_->triangles.size()
            << " triangle(s). The layout is asked for the matrix with holes; Stage 11 "
            << "puts each O-grid back from a template.";
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }
}

// The dialog. Opened only when there is something to decide, so a model with no
// circular inclusion never sees it and the pipeline is unchanged for every
// model that had none.
//
// The list is the point of opening a dialog rather than assuming: whether these
// components are the inclusions the model is about, or something the fit
// happened to accept, is a question about this model, and the radii and the
// deviation of each fit are how it is answered.
bool CrossGenWidget::promptDiskExcision(const std::vector<DiskTemplate::Inclusion> &found) {
    QDialog dlg(this);
    dlg.setWindowTitle("Stage 0c — circular inclusions");

    std::ostringstream oss;
    oss << found.size() << " circular inclusion(s) found:\n";
    const std::size_t shown = std::min<std::size_t>(found.size(), 12);
    for (std::size_t i = 0; i < shown; ++i) {
        const DiskTemplate::Inclusion &inc = found[i];
        oss << "  material " << inc.matId << "  r = " << std::fixed << std::setprecision(4)
            << inc.circle.radius << "  at (" << inc.circle.center[0] << ", "
            << inc.circle.center[1] << ")  fit " << std::setprecision(1)
            << (100.0 * inc.circle.maxRelDev) << "%\n" << std::setprecision(4);
    }
    if (shown < found.size()) oss << "  ... and " << (found.size() - shown) << " more\n";

    auto *summary = new QLabel(QString::fromStdString(oss.str()), &dlg);
    summary->setTextFormat(Qt::PlainText);

    auto *why = new QLabel(
        "Excising them asks the layout for the matrix with a hole where each disk\n"
        "was, and fills the holes back in from an O-grid template after Stage 10.\n"
        "A disk carries eight cones — four +1 inside, four −1 in the matrix — and\n"
        "it is placing those, many times over, that the Sec. 3.3 continuation\n"
        "cannot do. Declining lays them out like any other region.", &dlg);
    why->setTextFormat(Qt::PlainText);

    auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dlg);
    buttons->button(QDialogButtonBox::Ok)->setText("Excise");
    buttons->button(QDialogButtonBox::Cancel)->setText("Leave them in");
    QObject::connect(buttons, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
    QObject::connect(buttons, &QDialogButtonBox::rejected, &dlg, &QDialog::reject);

    auto *layout = new QVBoxLayout(&dlg);
    layout->addWidget(summary);
    layout->addWidget(why);
    layout->addWidget(buttons);

    // The dialog remembers the last answer only in so far as it decides which
    // button is the default one, since the decision itself is per model.
    buttons->button(diskSettings_.excise ? QDialogButtonBox::Ok : QDialogButtonBox::Cancel)
        ->setDefault(true);

    const bool accepted = dlg.exec() == QDialog::Accepted;
    diskSettings_.excise = accepted;
    return accepted;
}

// ── Stage 11: the O-grid templates ──────────────────────────────────────────
//
// Straight after Stage 10 and off the same rim arcs it was given, because the
// two are one decision: the fill exists only if the rim carries an even number
// of edges, and that is a property of the interval assignment. A no-op when
// Stage 0c excised nothing.
void CrossGenWidget::runDiskTemplates(const std::vector<std::vector<int>> &rims) {
    diskFill_.reset();
    if (inclusions_.empty() || rims.empty()) return;
    if (!quadMesh_.has_value() || !arrangement_.has_value()) return;

    DiskTemplate::Options topts;
    topts.coreSquareness   = diskSettings_.squareness;
    topts.ringDepth        = diskSettings_.ringDepth;
    topts.smoothingPasses  = diskSettings_.smoothing;

    auto t0 = Clock::now();
    try {
        diskFill_.emplace(quadMesh_->vertices(), quadMesh_->quads(),
                          quadMesh_->quadMaterials(),
                          DiskTemplate::rimVertexLoops(*arrangement_, *quadMesh_, rims),
                          inclusions_, topts);
    } catch (const std::exception &e) {
        diskFill_.reset();
        console_.log(std::string("[Disks] Stage 11 FAILED: ") + e.what());
        return;
    }
    auto t1 = Clock::now();

    const DiskTemplate::Report &dr = diskFill_->getReport();
    {
        std::ostringstream oss;
        oss << "[Disks] Stage 11: " << dr.filled << " of " << dr.inclusions
            << " inclusion(s) templated -- " << dr.blocks << " block(s), " << dr.quads
            << " element(s), " << dr.vertices << " new vertex/vertices, "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }
    if (dr.refusedOdd || dr.refusedShort || dr.refusedOpen) {
        std::ostringstream oss;
        oss << "[Disks] refused: " << dr.refusedOdd << " rim(s) with an odd edge count, "
            << dr.refusedShort << " too short to hold a ring and a core, "
            << dr.refusedOpen << " Stage 10 did not close -- each is left as a hole";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Disks] merged mesh: " << dr.mergedVertices << " vertices, "
            << dr.mergedQuads << " quad(s); scaled Jacobian " << std::fixed
            << std::setprecision(4) << dr.minScaledJacobian << " worst, "
            << dr.meanScaledJacobian << " mean (" << dr.templateMinScaledJacobian
            << " worst on the templates)";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Disks] " << dr.interiorEdges << " interior and " << dr.boundaryEdges
            << " boundary edge(s), " << dr.nonManifoldEdges << " used a third time, "
            << dr.cracks << " crack(s): "
            << (dr.conforming ? "conforming [PASS]" : "not conforming [FAIL]");
        console_.log(oss.str());
    }
    for (const std::string &m : dr.messages) console_.log("[Disks] " + m);
    {
        std::ostringstream oss;
        oss << "[Disks] Stage 11 " << (dr.valid ? "validates [PASS]"
                                                : "does not validate [FAIL]");
        console_.log(oss.str());
        std::cerr << "[Viewer] " << oss.str() << "\n";
    }
}

// Stages 5 to 7 again from the immersion already computed, at whatever the
// dialog was last left at.
//
// Stage 4 is not repeated: psi_R is a property of the flat metric and the cut,
// and nothing in this dialog touches either. What is repeated is the labelling,
// the continuation and the tracing -- which is the expensive part, and the
// point: how near a separatrix has to pass a cone before the pair is said to be
// connected is a judgement about the model, and the only way to settle one is
// to try a number and look.
void CrossGenWidget::rerunMERIDIANConnectivity() {
    if (!immersion_.has_value()) return;

    // Reverse construction order: SplineFit holds a reference to the
    // arrangement, the arrangement to the separatrices and the labels,
    // LayoutEnergy to the immersion and the labels, SubdomainLabels to the
    // immersion. Stage 11 is a copy of Stage 10's arrays and goes with them.
    pipelineBlocked_.clear();
    smoothMesh_.reset();
    diskFill_.reset();
    quadMesh_.reset();
    meshAttempted_ = tmopAttempted_ = false;
    splines_.reset();
    arrangement_.reset();
    patchesAttempted_ = false;
    patchesAnnounced_ = false;
    separatrices_.reset();
    meridianLayout_.reset();
    meridianLabels_.reset();

    auto t0 = Clock::now();
    try {
        meridianLabels_.emplace(*immersion_, meridianConn_.labels,
                                interfaces_.has_value() ? &*interfaces_ : nullptr);
    } catch (const std::exception &e) {
        meridianLabels_.reset();
        console_.log(std::string("[Labels] FAILED: ") + e.what());
        return;
    }
    const SubdomainLabels::Report &sr = meridianLabels_->getReport();
    {
        std::ostringstream oss;
        oss << "[Labels] Gamma_topo re-seeded: " << sr.topoPaths << " path(s) ("
            << sr.topoSelfReturns << " back to their own cone, " << sr.topoExtraPerPair
            << " a second curve of a pair) at a near-miss tolerance of " << std::fixed
            << std::setprecision(3) << meridianConn_.labels.nearMissTolerance
            << " of the mean cone spacing";
        console_.log(oss.str());
    }

    bool ok = false;
    try {
        meridianLayout_.emplace(*immersion_, *meridianLabels_);
        ok = meridianLayout_->run();
    } catch (const std::exception &e) {
        meridianLayout_.reset();
        console_.log(std::string("[Layout] FAILED: ") + e.what());
        return;
    }
    auto t1 = Clock::now();
    {
        const LayoutEnergy::Report &er = meridianLayout_->getReport();
        std::ostringstream oss;
        oss << "[Layout] re-run: " << er.outerSteps << " outer step(s), Q3 "
            << std::scientific << std::setprecision(2) << er.maxBoundaryResidual << ", Q4 "
            << er.maxSeamResidual << ", Q5 " << er.maxTopoResidual << "; Psi "
            << (ok ? "satisfies Q1-Q5 [PASS]" : "is not yet a layout [FAIL]") << ", "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }

    separatricesAttempted_ = false;
    separatricesAnnounced_ = true;   // the notice has been on screen once already
    runMERIDIANSeparatrices();
}

// ── phase advancement ────────────────────────────────────────────────────────

void CrossGenWidget::advancePhase() {
    if (mode_ == Mode::PolyVector) {
        Phase old = phase_;
        phase_ = nextPhase(phase_);
        if (phase_ != old)
            std::cerr << "[Viewer] Phase " << phaseName(phase_) << "\n";
    } else if (mode_ == Mode::ZIPLINE) {
        ZIPLINEPhase old = ziplinePhase_;
        ziplinePhase_ = nextZIPLINEPhase(ziplinePhase_);
        if (ziplinePhase_ != old)
            std::cerr << "[Viewer] ZIPLINE Phase " << ziplinePhaseName(ziplinePhase_) << "\n";

        // The last two phases as every other mode has them: the mesh dialog on
        // the way into Mesh, the TMOP one on the way into Smoothed and on
        // every 'c' there, each marked as asked before it opens so that the
        // catch-up in runComputations() does not open a second one on top.
        if (old == ZIPLINEPhase::Blocks && ziplinePhase_ == ZIPLINEPhase::Mesh) {
            // 'c' pressed again before the frame that builds the blocks.
            if (!ziplineBlocks()) buildZIPLINEBlocks();
            if (ziplineBlocks() && !ziplineBlocks()->blocks().empty()) {
                ziplineMeshAttempted_ = true;
                if (promptZIPLINEMesh()) runZIPLINEMesh();
            }
        }
        if (ziplinePhase_ == ZIPLINEPhase::Smoothed && ziplineMesh()) {
            tmopAttempted_ = true;
            if (promptTMOP()) runTMOP();
        }
    } else if (mode_ == Mode::MedialAxis) {
        MedialAxisPhase old = maPhase_;
        maPhase_ = nextMedialAxisPhase(maPhase_);
        if (maPhase_ != old)
            std::cerr << "[Viewer] Medial Axis Phase " << medialAxisPhaseName(maPhase_) << "\n";
    } else if (mode_ == Mode::UMBER) {
        UMBERPhase old = umberPhase_;
        umberPhase_ = nextUMBERPhase(umberPhase_);
        if (umberPhase_ != old)
            std::cerr << "[Viewer] UMBER Phase " << umberPhaseName(umberPhase_) << "\n";

        // The last two phases as the pipelines and ATLAS have them: the mesh
        // dialog on the way into Mesh, the TMOP one on the way into Smoothed
        // and on every 'c' there, each marked as asked before it opens so that
        // the catch-up in runComputations() does not open a second one on top.
        if (old == UMBERPhase::Decomposition && umberPhase_ == UMBERPhase::Mesh &&
            umberDecomp_.has_value() && !umberDecomp_->blocks.empty()) {
            umberMeshAttempted_ = true;
            if (promptUMBERMesh()) runUMBERMesh();
        }
        if (umberPhase_ == UMBERPhase::Smoothed && umberMesh_.has_value()) {
            tmopAttempted_ = true;
            if (promptTMOP()) runTMOP();
        }
    } else if (inPipeline()) {
        PipelinePhase old = pipePhase_;
        pipePhase_ = nextPipelinePhase(pipePhase_);
        if (pipePhase_ != old)
            std::cerr << "[Viewer] " << modeName(mode_) << " Phase "
                      << pipelinePhaseName(pipePhase_, mode_) << "\n";

        // Mesh is the last phase and it stays that way: 'c' at it re-opens the
        // Stage 10 dialog, so a target edge length can be tried, looked at, and
        // tried again, the way 'c' at UMBER's Simplified phase re-opens the
        // chord dialog. The connectivity dialog, which used to sit on this key
        // at the Patches phase, moved to 'n' -- the two are both judgements
        // about the model and both need re-opening, and one key cannot carry
        // them both.
        if (old == PipelinePhase::Patches && pipePhase_ == PipelinePhase::Mesh) {
            // Stages 8 and 9 are built by the frame after the one that asked
            // for them, so they may not be there yet on the frame that
            // advanced. Build them now rather than opening a dialog on a
            // pipeline that has not run.
            if (!patchesAttempted_) runMERIDIANPatches();
        }
        if (pipePhase_ == PipelinePhase::Mesh && splines_.has_value()) {
            // Marked as asked before the dialog opens, for the same reason the
            // catch-up in runComputations() does it in that order: the dialog
            // paints, painting reaches that catch-up, and a second dialog on
            // top of this one is what it would otherwise open. Cancelling
            // counts as having been asked, and 'c' here is what re-opens it.
            meshAttempted_ = true;
            if (promptMERIDIANMesh()) runMERIDIANMesh();
        }
        // Stage 12, the same way and for the same reasons: marked as asked
        // before the dialog opens, cancelling counts as having been asked, and
        // 'c' at this phase re-opens it. This is now the last phase, so 'c'
        // arrives here rather than at the Mesh one.
        if (pipePhase_ == PipelinePhase::Smoothed && quadMesh_.has_value()) {
            tmopAttempted_ = true;
            if (promptTMOP()) runTMOP();
        }
    } else if (mode_ == Mode::ATLAS) {
        ATLASPhase old = atlasPhase_;
        atlasPhase_ = nextATLASPhase(atlasPhase_);
        if (atlasPhase_ != old)
            std::cerr << "[Viewer] ATLAS Phase " << atlasPhaseName(atlasPhase_) << "\n";

        // The last two phases as the pipelines have them: the mesh dialog on
        // the way into Mesh, the TMOP one on the way into Smoothed and on every
        // 'c' there, each marked as asked before it opens.
        if (old == ATLASPhase::Search && atlasPhase_ == ATLASPhase::Blocks && !atlasAttempted_) {
            // 'c' pressed again before the frame that runs the search: run it
            // now rather than enter a phase whose input does not exist.
            runATLAS();
        }
        if (old == ATLASPhase::Blocks && atlasPhase_ == ATLASPhase::Mesh && atlas_ && atlas_->hasCover()) {
            atlasMeshAttempted_ = true;
            if (promptATLASMesh()) runATLASMesh();
        }
        if (atlasPhase_ == ATLASPhase::Smoothed && atlasMesh_.has_value()) {
            tmopAttempted_ = true;
            if (promptTMOP()) runTMOP();
        }
    }
}

// ── the DualMBO solve, which three modes share ──────────────────────────────────
//
// UMBER's input is a converged DualMBO field, and so is either pipeline's: Sec.
// 3.1 reads the cone indices off a field's holonomy, and TORSION goes on to
// integrate the very same field. So the guards ask about the stage rather than
// the mode.

bool CrossGenWidget::dualMBOStageWantsField() const {
    return (mode_ == Mode::UMBER && umberPhase_ >= UMBERPhase::CrossField) ||
           (inPipeline() && pipePhase_ >= PipelinePhase::CrossField);
}

bool CrossGenWidget::dualMBOStageIsStepping() const {
    return (mode_ == Mode::UMBER && umberPhase_ == UMBERPhase::Stepping) ||
           (inPipeline() && pipePhase_ == PipelinePhase::Stepping);
}

// ZIPLINE keeps the network it traced against (ZIPLINE::findInterfaces, only
// when there is more than one material); every other mode builds interfaces_.
const Interfaces *CrossGenWidget::interfaceNetwork() const {
    if (mode_ == Mode::ZIPLINE)
        return (zipline_ && zipline_->hasInterfaces()) ? &zipline_->getInterfaces() : nullptr;
    return interfaces_.has_value() ? &*interfaces_ : nullptr;
}

std::vector<int> CrossGenWidget::hangingTJunctions() const {
    std::vector<int> out;
    if (!ziplineQuant_.has_value() || !ziplineSimplified()) return out;
    const QuadLayout &sl = ziplineSimplified()->getLayout();
    const auto &nodes = sl.getNodes();
    const auto &arcs = sl.getArcs();
    for (size_t n = 0; n < nodes.size(); ++n) {
        // Structural T-junction: three arcs meeting in the interior at
        // something that is not a singularity, the same test the layout's
        // own report uses (QuadLayout::finish).
        const auto &node = nodes[n];
        if (node.darts.size() != 3 ||
            node.kind == QuadLayout::NodeKind::Singularity) {
            continue;
        }
        bool boundary = false;
        for (const int d : node.darts) {
            if (arcs[QuadLayout::arcOfDart(d)].onBoundary) boundary = true;
        }
        if (boundary) continue;

        // Resolved by the quantization unless an incident edge was forced
        // to zero (no tick lands on the junction, so no grid line leaves
        // it) or an incident arc never became a T-mesh edge at all (it
        // borders a skipped component, where no grid exists to weld with).
        bool hanging = false;
        for (const int d : node.darts) {
            const int e = ziplineQuant_->edgeOfArc[QuadLayout::arcOfDart(d)];
            if (e < 0 || ziplineQuant_->tmesh.edges[e].x <= 0) {
                hanging = true;
                break;
            }
        }
        if (hanging) out.push_back(static_cast<int>(n));
    }
    return out;
}

// ── ZIPLINE's blocks, mesh and dialog ─────────────────────────────────────
//
// Sec. 4's partition read as the shared BlockDecomposition on spline geometry
// (ZIPLINE/LayoutBlocks.hxx), then meshed and smoothed by the classes UMBER
// uses. A component with a T-junction on a side is not a block: T-junctions
// are not meshed yet, so it stays grey with a red disk on the junction, and the
// coverage the console prints is how much of the model is left once those are
// taken out.
void CrossGenWidget::buildZIPLINEBlocks() {
    ziplineUncoveredTris_.clear();
    smoothMesh_.reset();
    ziplineMeshAttempted_ = false;
    tmopAttempted_ = false;
    if (!ziplineSimplified() || !mesh_) return;

    auto t0 = Clock::now();
    try {
        // Discards any mesh on the old blocks with them.
        zipline_->buildBlocks();
    } catch (const std::exception &e) {
        console_.log(std::string("[Blocks] FAILED: ") + e.what());
        std::cerr << "[Viewer] LayoutBlocks failed: " << e.what() << "\n";
        return;
    }
    auto t1 = Clock::now();

    // The uncovered area as triangles of the model, by centroid: what the
    // picture fills so that it reads as a region rather than as an outline.
    {
        const QuadLayout &sl = ziplineSimplified()->getLayout();
        const std::vector<char> &isBlock = ziplineBlocks()->faceIsBlock();
        for (size_t f = 0; f < sl.getFaces().size(); ++f) {
            if (f < isBlock.size() && isBlock[f]) continue;
            std::vector<Point> poly;
            for (const int d : sl.getFaces()[f].darts) {
                std::vector<Point> pts = sl.getArcs()[QuadLayout::arcOfDart(d)].pts;
                if (d & 1) std::reverse(pts.begin(), pts.end());
                poly.insert(poly.end(), pts.begin(), pts.end() - 1);
            }
            if (poly.size() < 3) continue;
            Point lo = poly[0], hi = poly[0];
            for (const Point &p : poly) {
                lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
                hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
            }
            for (size_t t = 0; t < mesh_->triangles.size(); ++t) {
                const Triangle &tri = mesh_->triangles[t];
                const Point c = (mesh_->vertices[tri[0]] + mesh_->vertices[tri[1]] +
                                 mesh_->vertices[tri[2]]) / 3.0;
                if (c[0] < lo[0] || c[0] > hi[0] || c[1] < lo[1] || c[1] > hi[1]) continue;
                int wind = 0;
                for (size_t i = 0; i < poly.size(); ++i) {
                    const Point &a = poly[i], &b = poly[(i + 1) % poly.size()];
                    if (a[1] <= c[1]) {
                        if (b[1] > c[1] && cross2(b - a, c - a) > 0.0) ++wind;
                    } else if (b[1] <= c[1] && cross2(b - a, c - a) < 0.0) {
                        --wind;
                    }
                }
                if (wind != 0) ziplineUncoveredTris_.push_back(static_cast<int>(t));
            }
        }
    }

    const LayoutBlocks::Report &r = ziplineBlocks()->getReport();
    {
        std::ostringstream oss;
        oss << "[Blocks] " << r.blocks << " block(s) of " << r.faces << " component(s), "
            << ziplineBlocks()->decomposition().edges.size() << " macro edge(s), "
            << ziplineBlocks()->decomposition().vertices.size() << " macrovertex/-ices: "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
    }
    if (r.blocks < r.faces) {
        std::ostringstream oss;
        oss << "[Blocks] not blocks: " << r.tJunctionFaces << " with a T-junction on a side ("
            << r.tJunctions << " T-junction(s), in red)";
        if (r.notFourSided > 0) oss << ", " << r.notFourSided << " not four-sided";
        if (r.notDisks > 0) oss << ", " << r.notDisks << " not a disk";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(3) << "[Blocks] splines: " << r.fittedArcs
            << " side(s) fitted as cubic B-splines (up to " << r.maxSegmentsUsed
            << " segments; " << r.meanDeviation << " mean, " << r.maxDeviation
            << " worst deviation, in mean mesh edges), " << r.exactArcs
            << " on the boundary carried exactly";
        if (r.foldedPatches > 0) oss << "; " << r.foldedPatches << " Coons patch(es) fold";
        console_.log(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[Blocks] B-rep: " << r.brepFaces << " face(s), " << r.brepEdges << " edge(s) ("
            << r.brepSharedEdges << " shared, " << r.brepFreeEdges << " free), "
            << (r.brepValid ? "valid" : "NOT valid");
        if (r.brepFailures > 0) oss << ", " << r.brepFailures << " refused by the kernel";
        console_.log(oss.str());
    }
    for (const std::string &m : r.messages) console_.log("[Blocks] " + m);

    // The number to read the method by, said every time and loudly when it is
    // short: what is not covered is not meshed.
    const double coverage = ziplineBlocks()->coverage();
    std::ostringstream cov;
    cov << std::fixed << std::setprecision(1) << "[Blocks] covering " << 100.0 * coverage
        << "% of the model";
    if (coverage < 0.999)
        cov << " -- the shaded components outlined in grey are the rest, and they get no elements";
    console_.log(cov.str());
    std::cerr << "[Viewer] " << cov.str() << "\n";
    if (r.straddlingBlocks > 0) {
        std::ostringstream bad;
        bad << "[Blocks] " << r.straddlingBlocks << " block(s) sit in more than one of the "
            << r.materials << " material(s): the separatrices follow a cross field that is "
               "not aligned to the interfaces, and the elements in those blocks will straddle one";
        console_.log(bad.str());
    }
}

bool CrossGenWidget::promptZIPLINEMesh() {
    if (!ziplineBlocks()) return false;
    return promptBlockQuadMesh(ziplineBlocks()->decomposition(),
                               "ZIPLINE — quadrilateral mesh on the blocks");
}

// BlockQuadMesh at the dialog's settings, as UMBER's Mesh phase runs it.
void CrossGenWidget::runZIPLINEMesh() {
    if (!ziplineBlocks() || ziplineBlocks()->blocks().empty()) {
        blockPipeline("the mesh not built", "the tracing left no conforming block to mesh");
        return;
    }
    ziplineMeshAttempted_ = true;
    smoothMesh_.reset();
    tmopAttempted_ = false;
    pipelineBlocked_.clear();

    BlockQuadMesh::Options mo;
    mo.targetEdgeLength   = meshSettings_.target;
    mo.minIntervals       = meshSettings_.minEdges;
    mo.maxIntervals       = meshSettings_.maxEdges;
    mo.smoothingPasses    = TORSION::Options().quadSmoothingPasses;
    mo.smoothingThreshold = TORSION::Options().quadSmoothingThreshold;

    auto t0 = Clock::now();
    try {
        zipline_->buildMesh(mo);
    } catch (const std::exception &e) {
        blockPipeline("the mesh failed", e.what());
        return;
    }
    auto t1 = Clock::now();
    reportBlockQuadMesh(*ziplineMesh(), std::chrono::duration<double, std::milli>(t1 - t0).count());
    if (ziplineBlocks()->coverage() < 0.999) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(1) << "[Mesh] " << 100.0 * ziplineBlocks()->coverage()
            << "% of the model meshed: the shaded components carry a T-junction or are not "
               "four-sided";
        console_.log(oss.str());
    }
}

// ── lazy computations ────────────────────────────────────────────────────────

void CrossGenWidget::runComputations() {
    // ── ZIPLINE: Initialize CrossField ────────────────────────────────────────
    if (mode_ == Mode::ZIPLINE && ziplinePhase_ >= ZIPLINEPhase::CrossField && !ziplineField()) {
        // Stage 0b first, as TestZIPLINE runs it (ZIPLINE::findInterfaces): the
        // field is aligned to the material interfaces and each disk inclusion's
        // centre is pinned, and the tracing emits from the network's nodes and
        // cuts its branches where separatrices cross them. A network that
        // cannot be built is reported and the model traced as one material.
        if (!zipline_) {
            zipline_ = std::make_unique<ZIPLINE>(mesh_);
            try {
                zipline_->findInterfaces();
            } catch (const std::exception &e) {
                console_.log(std::string("[Interfaces] FAILED: ") + e.what());
                ZIPLINE::Options single;
                single.materialInterfaces = false;
                zipline_ = std::make_unique<ZIPLINE>(mesh_, single);
                zipline_->findInterfaces();
            }
            if (zipline_->hasInterfaces()) {
                const Interfaces::Report &ir = zipline_->getInterfaces().getReport();
                std::ostringstream oss;
                oss << "[Interfaces] " << ir.materials << " material(s), " << ir.interfaceEdges
                    << " interface edge(s) in " << ir.branches << " branch(es) (" << ir.closedLoops
                    << " closed), " << ir.nodes << " node(s): " << ir.junctions << " junction(s), "
                    << ir.landings << " landing(s), " << ir.kinks << " kink(s)";
                if (ir.illPosedNodes > 0) oss << ", " << ir.illPosedNodes << " ill-posed";
                console_.log(oss.str());
            }
        }
        auto t0 = Clock::now();
        zipline_->startField();
        auto t1 = Clock::now();
        console_.log("[MBO] Initialized CrossField: " +
                     formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count()));
    }

    // ── ZIPLINE: Kick off stepping ────────────────────────────────────────────
    if (mode_ == Mode::ZIPLINE && ziplinePhase_ == ZIPLINEPhase::Stepping &&
        ziplineField() && !mboSteppingStarted_) {
        mboSteppingStarted_ = true;
        console_.log("[MBO] Starting " + std::to_string(zipline_->getOptions().fieldMaxSteps) +
                     " iterations...");
    }

    // ── ZIPLINE: Run 2 stepping iterations per frame ──────────────────────────
    if (mode_ == Mode::ZIPLINE && ziplinePhase_ == ZIPLINEPhase::Stepping &&
        mboSteppingStarted_ && ziplineField() && !zipline_->fieldFinished()) {
        zipline_->stepField(2);
        const ZIPLINE::Status &st = zipline_->getStatus();
        if (st.fieldConverged)
            console_.log("[MBO] Convergence at step " + std::to_string(st.fieldSteps) +
                         " error=" + std::to_string(st.fieldError));
        std::ostringstream oss;
        oss << "[MBO] Step " << st.fieldSteps << "/" << zipline_->getOptions().fieldMaxSteps;
        console_.log(oss.str());
    }

    // ── ZIPLINE: Build SeparatrixTrace ────────────────────────────────────────
    if (mode_ == Mode::ZIPLINE && ziplinePhase_ >= ZIPLINEPhase::Separatrices &&
        ziplineField() && !ziplineTrace()) {
        // The field is finished first (ZIPLINE::startTrace), so a 'c' pressed
        // during the stepping animation no longer traces a field that is still
        // moving; say so when that is what happened.
        const int stepsBefore = zipline_->getStatus().fieldSteps;
        const bool finishing = !zipline_->fieldFinished();
        auto t0 = Clock::now();
        zipline_->startTrace();
        auto t1 = Clock::now();
        if (finishing) {
            const ZIPLINE::Status &st = zipline_->getStatus();
            std::ostringstream fin;
            fin << "[MBO] finished the field before tracing it: steps " << stepsBefore << " -> "
                << st.fieldSteps << (st.fieldConverged ? ", converged" : ", out of steps")
                << ", error " << st.fieldError;
            console_.log(fin.str());
        }
        const SeparatrixTrace &trace = *ziplineTrace();
        std::ostringstream oss;
        oss << "[Separatrices] Initialized " << trace.separatrices.size()
            << " separatrices from " << trace.singularities.size()
            << " singularities: "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
        const SeparatrixTrace::Report &tr = trace.getReport();
        if (tr.interfaceNodes > 0 || trace.getInterfaces()) {
            std::ostringstream it;
            it << "[Separatrices] " << tr.interfaceEmitted << " of them from " << tr.interfaceNodes
               << " interface node(s), which absorbed " << tr.absorbedSingularities
               << " singular triangle(s) of the field";
            console_.log(it.str());
        }
        std::ostringstream ph;
        ph << "[Separatrices] Poincare-Hopf: sum of indices "
           << (tr.interiorIndexSum4 + tr.boundaryIndexSum4 + tr.interfaceIndexSum4) << "/4 against "
           << tr.eulerCharacteristic << (tr.poincareHopf ? " [ok]" : " [FAILS]");
        console_.log(ph.str());
    }

    // ── ZIPLINE: Step tracing one iteration per frame ─────────────────────────
    if (mode_ == Mode::ZIPLINE && ziplinePhase_ == ZIPLINEPhase::Trace &&
        ziplineTrace() && !ziplineTracingFinished_) {
        if (!ziplineTracingStarted_) {
            ziplineTracingStarted_ = true;
            console_.log("[Trace] Starting separatrix tracing...");
        }
        if (zipline_->stepTrace()) {
            ziplineTracingFinished_ = true;
            console_.log("[Trace] Tracing complete.");
        }
    }

    // ── ZIPLINE: Build the quad layout the separatrices cut out ───────────────
    if (mode_ == Mode::ZIPLINE && ziplinePhase_ >= ZIPLINEPhase::Layout && ziplineTrace() &&
        !ziplineLayout()) {
        // Skipping ahead past the trace animation is allowed, so finish the
        // tracing here rather than assuming a frame of it has run.
        if (!zipline_->traceFinished()) {
            zipline_->finishTrace();
            ziplineTracingFinished_ = true;
        }
        auto t0 = Clock::now();
        zipline_->buildLayout();
        auto t1 = Clock::now();

        const QuadLayout::Report &r = ziplineLayout()->getReport();
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

    // ── ZIPLINE: Sec. 4, collapse the chords the layout does not need ─────────
    if (mode_ == Mode::ZIPLINE && ziplinePhase_ >= ZIPLINEPhase::Simplified && ziplineLayout() &&
        !ziplineSimplified()) {
        auto t0 = Clock::now();
        zipline_->simplifyLayout();
        auto t1 = Clock::now();

        const PartitionSimplify::Report &r = ziplineSimplified()->getReport();
        const QuadLayout::Report &sl = ziplineSimplified()->getLayout().getReport();
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
        if (r.stems.attempted > 0) {
            std::ostringstream st;
            st << "[Simplify] Sec. 12: " << r.stems.extended << " of " << r.stems.attempted
               << " T-junction stem(s) traced on to the boundary, " << r.tJunctionsBeforeStems
               << " -> " << sl.tJunctions << " T-junction(s)";
            const int refused = r.stems.attempted - r.stems.extended;
            if (refused > 0)
                st << "; the rest " << r.stems.noBoundary << " never reached it, "
                   << r.stems.tangential << " ran alongside a separatrix, " << r.stems.throughNode
                   << " ran into a node, " << r.stems.tooManyCrossings << " crossed too many, "
                   << r.stems.invalid << " left a worse layout";
            console_.log(st.str());
        }
    }

    // ── ZIPLINE: the simplified layout as blocks, then the mesh on them ─────
    //
    // Not after the quantization but beside it: the blocks read the simplified
    // layout, not the quantized grid, so a model the quantizer refuses still
    // has its blocks. The mesh and TMOP dialogs are opened by advancePhase();
    // these two are the catch-ups for a phase reached without passing through
    // it, marked as asked before the dialog opens for the pipelines' reason
    // (the dialog paints, and painting comes back here).
    if (mode_ == Mode::ZIPLINE && ziplinePhase_ >= ZIPLINEPhase::Blocks && ziplineSimplified() &&
        !ziplineBlocks()) {
        buildZIPLINEBlocks();
    }
    if (mode_ == Mode::ZIPLINE && ziplinePhase_ == ZIPLINEPhase::Mesh && ziplineBlocks() &&
        !ziplineBlocks()->blocks().empty() && !ziplineMeshAttempted_) {
        ziplineMeshAttempted_ = true;
        if (promptZIPLINEMesh()) runZIPLINEMesh();
    }
    if (mode_ == Mode::ZIPLINE && ziplinePhase_ == ZIPLINEPhase::Smoothed && ziplineMesh() &&
        !tmopAttempted_) {
        tmopAttempted_ = true;
        if (promptTMOP()) runTMOP();
    }

    // ── ZIPLINE: quantize the block decomposition the layout already is (QGP) ──
    //
    // A simplified QuadLayout's faces are already the blocks -- Sec. 4 of
    // Viertel et al. is a block decomposition in the same sense Sec. 4 of
    // Campen et al. is -- and its arcs are already shared edges, so the
    // conversion is a relabeling rather than the geometric welding the
    // medial axis blocks need. xIdeal is 1 on every edge, same as there:
    // Stage II drives each edge to its minimum and acts as an automatic
    // block-merging operator on top of what chord collapse already did.
    if (mode_ == Mode::ZIPLINE &&
        (ziplinePhase_ == ZIPLINEPhase::Quantize || ziplinePhase_ == ZIPLINEPhase::Quantized) &&
        ziplineSimplified() && !ziplineQuant_.has_value()) {
        auto t0 = Clock::now();
        ziplineQuant_.emplace(makeQuantTMesh(ziplineSimplified()->getLayout()));
        auto t1 = Clock::now();

        const QuadLayoutQuant &lq = *ziplineQuant_;
        std::ostringstream oss;
        oss << "[Trace] T-mesh: " << lq.tmesh.edges.size() << " edges, "
            << lq.tmesh.faces.size() << " faces, " << lq.tmesh.rows.size()
            << " constraints, from " << ziplineSimplified()->getLayout().getFaces().size()
            << " components: "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
        if (lq.skippedFaces > 0) {
            console_.log("[Trace] " + std::to_string(lq.skippedFaces) +
                         " component(s) not four-sided, drawn unquantized "
                         "rather than left as a gap");
        }
        // The report's count is structural (QuadLayout::finish), so it can
        // be trusted after the collapses; the node kinds cannot.
        const int tJunctions = ziplineSimplified()->getLayout().getReport().tJunctions;
        if (tJunctions > 0) {
            console_.log("[Trace] " + std::to_string(tJunctions) +
                         " T-junction(s) in the block structure; the "
                         "quantized grids weld across a T-junction, so "
                         "each becomes a regular vertex of the result");
        }

        if (!lq.ok) {
            console_.log("[Trace] Quantization aborted: " + lq.error);
        } else {
            auto q0 = Clock::now();
            ziplineQuantReport_ = TMeshQuantizer(ziplineQuant_->tmesh).run();
            auto q1 = Clock::now();

            int total = 0, maxLen = 0;
            for (const auto &e : ziplineQuant_->tmesh.edges) {
                total += e.x;
                maxLen = std::max(maxLen, e.x);
            }
            std::ostringstream qss;
            qss << "[Trace] Quantized: " << ziplineQuantReport_.stage1Vectors
                << " stage-I strips, " << ziplineQuantReport_.stage2Moves << "/"
                << ziplineQuantReport_.stage2Tried << " stage-II moves, objective "
                << std::fixed << std::setprecision(3) << ziplineQuantReport_.objective
                << ", " << total << " quads across the boundary (longest edge "
                << maxLen << "): "
                << formatMs(std::chrono::duration<double, std::milli>(q1 - q0).count());
            console_.log(qss.str());
            if (ziplineQuantReport_.forcedZeroEdges > 0) {
                std::ostringstream zss;
                zss << "[Trace] " << ziplineQuantReport_.forcedZeroEdges
                    << " edges are forced to zero by the decomposition itself: "
                       "no consistent assignment can lift them";
                console_.log(zss.str());
            }
            const size_t hanging = hangingTJunctions().size();
            if (hanging > 0) {
                console_.log("[Trace] WARNING: " + std::to_string(hanging) +
                             " T-junction(s) left hanging by zero edges or "
                             "skipped components -- shown in red");
            }
            if (!ziplineQuantReport_.consistent) {
                console_.log("[Trace] WARNING: quantization violates the "
                             "consistency system Ax = 0");
            }
        }
    }

    // ── DualMBO: Initialize ──────────────────────────────────────────────────────
    if (dualMBOStageWantsField() && !dualMBOField_.has_value()) {
        auto t0 = Clock::now();
        // gamma and the step cap are the pipeline's, not DualMBO's constructor
        // defaults: the cap in particular is 100 there and 500 in both
        // pipelines, and it is the cap runMBO() below would otherwise stop at
        // if the stepping phase were walked past before the field converged.
        const TORSION::Options pipe;
        dualMBOField_.emplace(mesh_, pipe.dualMBOMaxSteps, pipe.dualMBOGamma);
        // On a multi-material domain the interfaces are Dirichlet data for the
        // field in exactly the way dS is: a curve the layout has to keep needs
        // the field tangent to it, or neither region either side gets the
        // interior cones its own Gauss-Bonnet count demands. Must be set before
        // initialize(), which does the first assembly.
        // TORSION needs this more than MERIDIAN does rather than less: it
        // integrates the field, so an interface the field ran straight through
        // is an interface the *map* runs straight through.
        if ((inPipeline() || mode_ == Mode::UMBER) && interfaces_.has_value() &&
            interfaces_->multiMaterial()) {
            dualMBOField_->setAlignedInteriorEdges(interfaces_->interfaceEdges());
        }
        dualMBOField_->initialize();
        auto t1 = Clock::now();
        console_.log("[DualMBO] Initialized: " +
                     formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count()));
    }

    // ── DualMBO: Kick off stepping ───────────────────────────────────────────────
    if (dualMBOStageIsStepping() && dualMBOField_.has_value() && !dualMBOSteppingStarted_) {
        dualMBOSteppingStarted_ = true;
        dualMBOStepCount_       = 0;
        console_.log("[DualMBO] Starting MBO iterations (" +
                     std::to_string(dualMBOField_->getMesh().triangles.size()) + " tris)...");
    }

    // ── DualMBO: Run 2 stepping iterations per frame ─────────────────────────────
    if (dualMBOStageIsStepping() && dualMBOSteppingStarted_ && !dualMBOConverged_) {
        double ntris = static_cast<double>(mesh_->triangles.size());
        const int cap = TORSION::Options().dualMBOMaxSteps;
        for (int i = 0; i < 2 && dualMBOStepCount_ < cap; ++i) {
            dualMBOField_->step();
            ++dualMBOStepCount_;
            if (dualMBOField_->error < 2.0 * ntris * DUALMBO_TOL) {
                console_.log("[DualMBO] Converged at step " + std::to_string(dualMBOStepCount_) +
                             " error=" + std::to_string(dualMBOField_->error));
                dualMBOConverged_ = true;
                break;
            }
        }
        dualMBOField_->computeSingularities();

        std::ostringstream stepMsg;
        stepMsg << "[DualMBO] Step " << dualMBOStepCount_ << "  error=" << std::scientific
                << std::setprecision(3) << dualMBOField_->error;
        console_.log(stepMsg.str());
    }

    // ── UMBER: cuts + Eq. (1) on top of the DualMBO field ────────────────────────
    if (mode_ == Mode::UMBER && umberPhase_ >= UMBERPhase::Frames &&
        dualMBOField_.has_value() && !umberAttempted_) {
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

    // ── UMBER: the layout, and the decomposition read off it ─────────────────
    //
    // 'c' pressed twice in quick succession can land a later phase while the
    // structure it stands on has not been built, and then the dialog that
    // belongs to that phase never opens. The three catch-ups below are that
    // case, and each follows the pipelines' discipline: marked as asked
    // *before* a dialog opens, since the dialog runs a nested event loop and
    // that loop paints and painting comes back here. Cancelling counts as
    // having been asked; 'e' and 'c' re-open them.
    if (mode_ == Mode::UMBER && umberPhase_ >= UMBERPhase::Decomposition &&
        blocks_.has_value() && !blockLayout_.has_value()) {
        buildBlockLayout();
    }
    if (mode_ == Mode::UMBER && umberPhase_ == UMBERPhase::Mesh && umberDecomp_.has_value() &&
        !umberDecomp_->blocks.empty() && !umberMeshAttempted_) {
        umberMeshAttempted_ = true;
        if (promptUMBERMesh()) runUMBERMesh();
    }
    if (mode_ == Mode::UMBER && umberPhase_ == UMBERPhase::Smoothed && umberMesh_.has_value() &&
        !tmopAttempted_) {
        tmopAttempted_ = true;
        if (promptTMOP()) runTMOP();
    }

    // ── Both pipelines: Stage 0b, the material interfaces ────────────────────
    //
    // At the first phase, not the cone one: it reads the tags and nothing else,
    // and the picture of what the layout will have to keep is worth having in
    // front of the field rather than after it.
    // UMBER needs it for the same reason and at the same moment: its field,
    // its frame, its deformation and its block structure all have to know the
    // interfaces are there, and every one of them is built from this.
    if ((inPipeline() || mode_ == Mode::UMBER) && !interfacesAttempted_) {
        runMERIDIANInterfaces();
    }

    // ── Both pipelines: Stage 1, the cones and Eq. (4) ───────────────────────
    if (inPipeline() && pipePhase_ >= PipelinePhase::Cones &&
        dualMBOField_.has_value() && !conesAttempted_) {
        runMERIDIANCones();
    }

    // ── Both pipelines: Stage 2, the cutting graph ───────────────────────────
    if (inPipeline() && pipePhase_ >= PipelinePhase::Cut &&
        cones_.has_value() && !cutAttempted_) {
        runMERIDIANCut();
    }

    // ── MERIDIAN: Stage 3, discrete Ricci flow ───────────────────────────────
    if (mode_ == Mode::MERIDIAN && pipePhase_ >= PipelinePhase::Flow &&
        cones_.has_value() && !ricciAttempted_) {
        if (!ricciAnnounced_) {
            // runComputations() runs at the top of paintGL, so returning here
            // lets this frame draw the notice; the solve starts on the next.
            console_.log("[Ricci] solving Eq. (10) for the flat cone metric, this blocks...");
            ricciAnnounced_ = true;
        } else {
            runRicciFlow();
        }
    }

    // ── TORSION: Stage 3F, the combed field and its matchings ────────────────
    //
    // The same phase, and the same place in the argument: this is where the
    // pipeline stops working with the model's own geometry and starts working
    // with the structure the layout is going to have. Pipeline A gets that
    // structure from a flow; this gets it from the field it already has.
    if (mode_ == Mode::TORSION && pipePhase_ >= PipelinePhase::Flow &&
        coneCut_.has_value() && !framesAttempted_) {
        runTORSIONFrames();
    }

    // ── TORSION: Stages 4F, 4R and 4, psi_0 ──────────────────────────────────
    if (mode_ == Mode::TORSION && pipePhase_ >= PipelinePhase::Metric &&
        frames_.has_value() && !integrationAttempted_) {
        if (!integrationAnnounced_) {
            console_.log("[psi_0] integrating the field on Omega and untangling what it "
                         "inverted, this blocks...");
            integrationAnnounced_ = true;
        } else {
            runTORSIONIntegration();
        }
    }

    // ── Both pipelines: Stages 4 to 6, psi_0 through Psi ─────────────────────
    //
    // What has to be there before this runs differs by route -- Pipeline A
    // needs the flow and the metric it produced, because Stage 4 is the
    // unfolding of that metric and happens here; Pipeline B needs the immersion
    // itself, because its Stage 4 already happened at the phase before -- so
    // the guard is the one thing the two have in common at this point, which is
    // that the map psi_0 exists.
    if (inPipeline() && pipePhase_ >= PipelinePhase::Layout && !layoutAttempted_ &&
        ((mode_ == Mode::MERIDIAN && ricci_.has_value() && !flatMetric_.edges.empty()) ||
         (mode_ == Mode::TORSION && immersion_.has_value()))) {
        if (!layoutAnnounced_) {
            console_.log(mode_ == Mode::TORSION
                             ? "[Layout] labelling the subdomains and running the Eq. (13) "
                               "continuation on psi_0, this blocks..."
                             : "[Layout] unfolding psi_R and running the Eq. (13) "
                               "continuation, this blocks...");
            layoutAnnounced_ = true;
        } else {
            runMERIDIANLayout();
        }
    }

    // ── Both pipelines: Stage 7, the separatrices of Psi ─────────────────────
    if (inPipeline() && pipePhase_ >= PipelinePhase::Separatrices &&
        meridianLayout_.has_value() && !separatricesAttempted_) {
        if (!separatricesAnnounced_) {
            console_.log("[Separatrices] tracing the integral curves out of every cone, "
                         "this blocks...");
            separatricesAnnounced_ = true;
        } else {
            runMERIDIANSeparatrices();
        }
    }

    // ── Both pipelines: Stages 8 and 9, the arrangement and the patches ─────
    //
    // Only on a clean trace. Stage 8 builds a planar subdivision out of the
    // curves and Stage 9 fits a surface to its faces; run either on a bundle
    // that did not close and what comes out is a picture of blocks that are
    // not there. What Stage 7 already reported is exactly the precondition --
    // Q5 on every curve, and Psi still a layout after the repair -- so it is
    // read rather than re-derived.
    if (inPipeline() && pipePhase_ >= PipelinePhase::Patches &&
        separatrices_.has_value() && !patchesAttempted_) {
        if (!patchesAnnounced_) {
            console_.log("[Patches] building the arrangement and fitting the splines, "
                         "this blocks...");
            patchesAnnounced_ = true;
        } else {
            runMERIDIANPatches();
        }
    }

    // ── Both pipelines: Stage 10, when the phase arrived before Stage 9 ──────
    //
    // The mesh dialog normally opens on the keypress that enters this phase,
    // which is the right moment: it is a judgement about the model and it wants
    // making with the patches on screen. But that keypress can land while
    // Stages 6 and 7 still have the GUI thread, and then it enters a phase
    // whose input does not exist yet -- the dialog does not open, and nothing
    // afterwards asks for it. So the phase asks again, once, as soon as Stage 9
    // has something to mesh. Cancelling counts as having asked: 'c' re-opens it.
    if (inPipeline() && pipePhase_ == PipelinePhase::Mesh && splines_.has_value() &&
        !meshAttempted_) {
        // Marked as asked *before* the dialog opens, not after: the dialog runs
        // a nested event loop, that loop paints, and painting comes back here.
        // Cancelling counts as having been asked; 'c' re-opens it.
        meshAttempted_ = true;
        if (promptMERIDIANMesh()) runMERIDIANMesh();
    }

    // ── Both pipelines: Stage 12, when the phase arrived before Stage 10 ─────
    //
    // The same catch-up as above and for the same reason: 'c' pressed twice in
    // quick succession can land this phase while Stage 10 has not run, and then
    // the dialog that belongs to it never opens. Cancelling counts as having
    // asked; 'c' re-opens it.
    if (inPipeline() && pipePhase_ == PipelinePhase::Smoothed && quadMesh_.has_value() &&
        !tmopAttempted_) {
        tmopAttempted_ = true;
        if (promptTMOP()) runTMOP();
    }

    // ── ATLAS: Stages 1 and 2, the search, and the catch-ups ─────────────────
    //
    // Each stage once per run, as in the pipelines. The search is announced a
    // frame ahead so the notice is on screen while the GUI thread waits on its
    // threads; the two dialogs are asked for again here in case the phase was
    // entered before its input existed.
    if (mode_ == Mode::ATLAS && atlasPhase_ >= ATLASPhase::Domain && !atlasDomainAttempted_)
        runATLASDomain();
    // Stage 1b is the one ATLAS stage nothing downstream reads -- the searches
    // solve their own -- so it is built only while its own phase is on screen,
    // and skipped outright if 'c' has already carried the viewer past it.
    if (mode_ == Mode::ATLAS && atlasPhase_ == ATLASPhase::Field && !atlasFieldAttempted_) {
        if (!atlasFieldAnnounced_) {
            console_.log("[Field] Stage 1b: solving the reference cross field (DualMBO, TORSION's "
                         "Stage 0 settings); this blocks...");
            atlasFieldAnnounced_ = true;
        } else {
            runATLASField();
        }
    }
    if (mode_ == Mode::ATLAS && atlasPhase_ >= ATLASPhase::Carrier && !atlasCarrierAttempted_)
        runATLASCarrier();
    if (mode_ == Mode::ATLAS && atlasPhase_ >= ATLASPhase::Search && !atlasAttempted_) {
        if (!atlasAnnounced_) {
            std::ostringstream oss;
            const ATLAS::Options o;
            const int coarse = atlasMultiMaterial_ ? 0 : static_cast<int>(o.coarseSpacings.size()) *
                                                             std::max(1, o.coarseSeeds);
            oss << "[Search] running Stages 3-6: the fine search";
            if (coarse > 0) oss << " and " << coarse << " coarse one(s), in parallel threads";
            oss << "; this blocks...";
            console_.log(oss.str());
            atlasAnnounced_ = true;
        } else {
            runATLAS();
        }
    }
    if (mode_ == Mode::ATLAS && atlasPhase_ >= ATLASPhase::Blocks && atlas_ && atlas_->hasCover() &&
        !atlasBlocksLogged_)
        logATLASBlocks();
    if (mode_ == Mode::ATLAS && atlasPhase_ == ATLASPhase::Mesh && atlas_ && atlas_->hasCover() &&
        !atlasMeshAttempted_) {
        atlasMeshAttempted_ = true;
        if (promptATLASMesh()) runATLASMesh();
    }
    if (mode_ == Mode::ATLAS && atlasPhase_ == ATLASPhase::Smoothed && atlasMesh_.has_value() &&
        !tmopAttempted_) {
        tmopAttempted_ = true;
        if (promptTMOP()) runTMOP();
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

    // ── Medial Axis: the coarse block decomposition (Sec. 4) ──────────────────
    if (mode_ == Mode::MedialAxis && maPhase_ >= MedialAxisPhase::TMesh &&
        medialAxisMap_.has_value() && !medialAxisTMesh_.has_value()) {
        auto t0 = Clock::now();
        medialAxisTMesh_.emplace(*medialAxisMap_);
        auto t1 = Clock::now();

        const auto &st = medialAxisTMesh_->stats();
        std::ostringstream oss;
        oss << "[MedialAxis] T-mesh: " << st.blocks << " blocks from " << st.zones
            << " zones (" << st.capZones << " caps) over " << st.coarseEdges
            << " coarse edges, target size " << std::fixed << std::setprecision(3)
            << medialAxisTMesh_->targetSize() << ": "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
        {
            std::ostringstream col;
            col << "[MedialAxis] templates: " << st.colorCounts[0] << " green (2 blocks), "
                << st.colorCounts[1] << " red (3), " << st.colorCounts[2]
                << " blue (1), " << st.colorCounts[3] << " purple (5)";
            console_.log(col.str());
        }
        {
            // A sharp corner left inside a block side is a kink no quad grid
            // can reproduce, so every one of them is forced to be a block
            // corner. Any that could not be is worth seeing.
            std::ostringstream cor;
            cor << "[MedialAxis] sharp corners: " << st.sharpCorners << " found, "
                << st.cornerCuts << " forced into block corners";
            if (st.cornersUnanchored > 0) cor << ", " << st.cornersUnanchored << " UNANCHORED";
            console_.log(cor.str());
        }
        if (st.skippedLoops || st.unpairedEdges) {
            std::ostringstream warn;
            warn << "[MedialAxis] Block warnings: " << st.skippedLoops
                 << " loops skipped, " << st.unpairedEdges << " unpaired edges";
            console_.log(warn.str());
        }
    }

    // ── Medial Axis: quantization of the block decomposition (QGP Sec. 6) ────
    //
    // The blocks carry no shared-edge identifiers, so the converter welds
    // them geometrically first: corners become nodes, block sides are split
    // where another block's corner lands on them, and the T-junctions the
    // red and purple templates leave on the mid-radii are recovered that
    // way. xIdeal is 1 on every edge -- the block-decomposition setting of
    // QGP Sec. 10.1, where Stage II drives each edge to its minimum and so
    // acts as an automatic block-merging operator.
    if (mode_ == Mode::MedialAxis && maPhase_ >= MedialAxisPhase::Quantize &&
        medialAxisTMesh_.has_value() && !blockQuant_.has_value()) {
        auto t0 = Clock::now();
        blockQuant_.emplace(makeQuantTMesh(*medialAxisTMesh_));
        auto t1 = Clock::now();

        const BlockQuant &bq = *blockQuant_;
        std::ostringstream oss;
        oss << "[MedialAxis] T-mesh welded: " << bq.nodes.size() << " nodes, "
            << bq.tmesh.edges.size() << " edges, " << bq.tmesh.faces.size()
            << " faces, " << bq.tmesh.rows.size() << " constraints, from "
            << medialAxisTMesh_->blocks.size() << " blocks: "
            << formatMs(std::chrono::duration<double, std::milli>(t1 - t0).count());
        console_.log(oss.str());
        if (bq.skippedBlocks > 0) {
            std::ostringstream skip;
            skip << "[MedialAxis] " << bq.skippedBlocks
                 << " blocks dropped as unusable cells: " << bq.skippedNonQuad
                 << " not four-sided, " << bq.skippedUnlocatable
                 << " with a corner off their outline, " << bq.skippedDegenerate
                 << " with two sides on one curve, " << bq.skippedOverlapping
                 << " overlapping a pair already in place";
            console_.log(skip.str());
        }

        if (!bq.ok) {
            console_.log("[MedialAxis] Quantization aborted: " + bq.error);
        } else {
            auto q0 = Clock::now();
            quantReport_ = TMeshQuantizer(blockQuant_->tmesh).run();
            auto q1 = Clock::now();

            int total = 0, maxLen = 0;
            for (const auto &e : blockQuant_->tmesh.edges) {
                total += e.x;
                maxLen = std::max(maxLen, e.x);
            }
            std::ostringstream qss;
            qss << "[MedialAxis] Quantized: " << quantReport_.stage1Vectors
                << " stage-I strips, " << quantReport_.stage2Moves << "/"
                << quantReport_.stage2Tried << " stage-II moves, objective "
                << std::fixed << std::setprecision(3) << quantReport_.objective
                << ", " << total << " quads across the boundary (longest edge "
                << maxLen << "): "
                << formatMs(std::chrono::duration<double, std::milli>(q1 - q0).count());
            console_.log(qss.str());
            if (quantReport_.forcedZeroEdges > 0) {
                std::ostringstream zss;
                zss << "[MedialAxis] " << quantReport_.forcedZeroEdges
                    << " edges are forced to zero by the block decomposition "
                       "itself: no consistent assignment can lift them";
                console_.log(zss.str());

                // Clear them out the way QGP Sec. 7.1.1 advises, by
                // deleting the cells they collapse and letting the
                // neighbours meet in the middle.
                const ContractReport cr = contractZeroEdges(*blockQuant_);
                std::ostringstream css;
                css << "[MedialAxis] Zero-cell cleanup: " << cr.mergedCells
                    << " collapsed cells merged into their neighbours, "
                    << cr.pointCells << " point cells dropped, " << cr.splitEdges
                    << " edges split to line the sides up, " << cr.removedEdges
                    << " edges removed";
                console_.log(css.str());
                if (!cr.ok) {
                    console_.log("[MedialAxis] Cleanup rejected, keeping the "
                                 "uncontracted T-mesh: " + cr.error);
                } else if (cr.remainingZero > 0) {
                    std::ostringstream rss;
                    rss << "[MedialAxis] " << cr.remainingZero
                        << " cells could not be merged and stay collapsed";
                    console_.log(rss.str());
                }
            }
            if (!quantReport_.consistent) {
                console_.log("[MedialAxis] WARNING: quantization violates the "
                             "consistency system Ax = 0");
            }
        }
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
}

// ── animation render paths ───────────────────────────────────────────────────

void CrossGenWidget::renderMBOAnimation() {
    viewer::drawAxis(view_);
    viewer::drawMesh(*mesh_);
    viewer::drawVertexCrossFieldUK(*mesh_, *ziplineField(), scale_);

    double ballRadius = 0.5 * avgEdge_;
    for (const auto &sig : ziplineField()->singularTriangles) {
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

void CrossGenWidget::renderDualMBOAnimation() {
    viewer::drawAxis(view_);
    viewer::drawMesh(*mesh_);
    viewer::drawTriangleCrossField(*mesh_, *dualMBOField_, scale_);

    double ballRadius = 0.5 * avgEdge_;
    for (const auto &[vertIdx, crossIndex] : dualMBOField_->singularVertices) {
        if (vertIdx < 0 || vertIdx >= static_cast<int>(mesh_->vertices.size())) continue;
        const Point &c = mesh_->vertices[vertIdx];
        if (crossIndex > 0)
            viewer::drawDisk3D(c, ballRadius, 0.2f, 0.2f, 0.95f);
        else
            viewer::drawDisk3D(c, ballRadius, 0.95f, 0.2f, 0.2f);
    }

    renderOverlay("DualMBO stepping...\npress 'q' to quit");
}

void CrossGenWidget::renderTraceAnimation() {
    viewer::drawAxis(view_);
    viewer::drawMesh(*mesh_);

    viewer::lineWidth(3.0f);
    for (const auto &sep : ziplineTrace()->separatrices) {
        if (sep.path.size() < 2) continue;
        if (sep.active)
            viewer::color3f(0.95f, 0.1f, 0.1f);
        else
            viewer::color3f(0.1f, 0.9f, 0.2f);
        glBegin(GL_LINE_STRIP);
        for (const auto &tp : sep.path)
            glVertex2d(tp.global_pos[0], tp.global_pos[1]);
        glEnd();
    }
    viewer::lineWidth(1.0f);

    renderOverlay("Tracing separatrices...\npress 'q' to quit");
}

// ── split-screen helpers ─────────────────────────────────────────────────────

// Whether the right half of the window is showing a parameter domain, which is
// what decides where a drag or a scroll lands.
bool CrossGenWidget::inUVSplitScreen() const {
    return
           // The chord collapse phase takes the whole window back: what it has
           // to show is the structure before against the structure after, and
           // both of those live in the model.
           (mode_ == Mode::UMBER && umberPhase_ >= UMBERPhase::Polysquare &&
            umberPhase_ <= UMBERPhase::Blocks && polysquare_.has_value()) ||
           // The gallery of unfolded cone fans is not a parameter domain, but
           // it is a second world with its own scale in the right half of the
           // window, which is all this predicate is really asking.
           (mode_ == Mode::MERIDIAN && pipePhase_ == PipelinePhase::Metric &&
            !coneFans_.empty()) ||
           // TORSION's Metric phase is a parameter domain in the ordinary
           // sense, because its Stage 4 already happened: the right half is
           // psi_0 itself rather than a gallery standing in for a metric.
           (mode_ == Mode::TORSION && pipePhase_ == PipelinePhase::Metric &&
            immersion_.has_value()) ||
           // Stages 4 to 6 do have a parameter domain in the ordinary sense:
           // psi_R and Psi are maps of Omega into the plane.
           (inPipeline() && pipePhase_ == PipelinePhase::Layout &&
            immersion_.has_value()) ||
           // Stage 7 draws the same domain again, with the curves on it.
           (inPipeline() && pipePhase_ == PipelinePhase::Separatrices &&
            immersion_.has_value());
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
    viewer::lineWidth(2.0f);
    viewer::color3f(0.55f, 0.55f, 0.55f);
    glBegin(GL_LINES);
    glVertex2f(static_cast<float>(halfW), 0.0f);
    glVertex2f(static_cast<float>(halfW), static_cast<float>(h));
    glEnd();
    viewer::lineWidth(1.0f);
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

// The model at whichever pipeline stage is current: the mesh (or the flat
// metric drawn over it, or the combed frames), the cutting graph once it
// exists, and the cones on top. Shared by the last four phases, by both
// pipelines and by the left half of the split screen, so that the cones stay in
// the same place and the same colours throughout.
void CrossGenWidget::renderMERIDIANModel() {
    // The stretch ramp is dropped again at the separatrix phase: what that one
    // is about is the curves, and a wireframe in the diverging ramp underneath
    // competes with them for exactly the colours they are drawn in.
    const bool showMetric = (pipePhase_ >= PipelinePhase::Metric &&
                             pipePhase_ != PipelinePhase::Separatrices &&
                             pipePhase_ != PipelinePhase::Patches &&
                             !flatMetric_.edges.empty());
    const bool showU = (pipePhase_ == PipelinePhase::Flow && ricciU_.size() > 0);
    // Pipeline B's Flow phase. The combed frames are drawn over a plain
    // wireframe rather than over a filled field, because what is being looked
    // at is where a_f *steps*, and that is a comparison between two adjacent
    // triangles which a fill under them would flood.
    const bool showFrames = (mode_ == Mode::TORSION &&
                             pipePhase_ == PipelinePhase::Flow && frames_.has_value());

    // At the patch and mesh phases nothing of the model underneath is drawn at
    // all -- not the wireframe, not dS, not the cones. The block decomposition
    // is the whole picture there: every patch side is already a curve of S
    // drawn in full, and a triangulation or a boundary loop under it is clutter
    // that was the point of the phase to get past. It is worse at the mesh
    // phase than at the patch one, where the quadrilaterals are near the size
    // of the triangles beneath them and the two grids read as one.
    //
    // Unless there is no decomposition to draw. A stage that refused leaves the
    // phase with nothing of its own, and dropping the model as well turns that
    // into an empty window -- which says "the viewer is broken" when what
    // happened is "Stage 8 refused, and here is why" on the overlay. So the
    // model is kept whenever the picture that was to replace it is not there.
    const bool replaced = quadMesh_.has_value() || arrangement_.has_value();
    if (pipePhase_ >= PipelinePhase::Patches && replaced) return;

    if (showMetric) {
        // The metric replaces the wireframe rather than covering it: every edge
        // is drawn, just coloured by what the flow did to it.
        viewer::drawFlatMetric(*mesh_, flatMetric_, 1.6f);
    } else if (showU) {
        // The conformal factor as a filled field, with the wireframe kept
        // translucent over it so the triangulation stays readable against the
        // saturated ends of the ramp -- the same treatment OASIS mode gives its
        // quasi-eigenfunction.
        viewer::drawScalarField(*mesh_, ricciU_, ricciUAbsMax_);
        viewer::drawMeshOverlay(*mesh_, 0.85f, 0.85f, 0.85f, 0.22f, 1.0f);
    } else {
        // The material fill goes under the wireframe, so the edges the tags
        // separate stay readable over it.
        if (showMaterialFill_ && interfaces_.has_value() && interfaces_->multiMaterial())
            viewer::drawMaterialFill(*mesh_, 0.28f);
        viewer::drawMesh(*mesh_);
    }

    // The cross field stays under the cone phase: the indices were read off its
    // holonomy, and a cone sitting where the field turns is the whole argument
    // for putting one there. It goes once the cutting graph arrives, which
    // would otherwise be lost among the arrows.
    if (pipePhase_ == PipelinePhase::Cones && dualMBOField_.has_value())
        viewer::drawTriangleCrossField(*mesh_, *dualMBOField_, scale_);

    // ... and comes back one phase later on the field route, as one branch of
    // itself rather than as four indistinguishable arms. The cutting graph is
    // drawn over it on purpose here: the arcs of G are exactly where a_f is
    // allowed to step, so the two pictures only mean anything together.
    if (showFrames)
        viewer::drawCombedFrames(*mesh_, *frames_, scale_);

    if (pipePhase_ >= PipelinePhase::Cut && coneCut_.has_value())
        viewer::drawCuttingGraph(*coneCut_, showMetric ? 2.0f : 3.5f);
    viewer::drawBoundaryEdges(*mesh_);

    // The interface network over all of it and under the cones. It is an input
    // rather than a result, so it is drawn at every phase: what each stage is
    // to be judged on is whether its own picture still respects these curves.
    if (showInterfaces_ && interfaces_.has_value())
        viewer::drawInterfaceNetwork(*interfaces_, 0.4 * avgEdge_, showMetric ? 2.0f : 3.0f);

    if (cones_.has_value())
        viewer::drawCones(*mesh_, *cones_, 0.5 * avgEdge_);
}

// ── normal render ─────────────────────────────────────────────────────────────

void CrossGenWidget::renderNormal() {
    if (mode_ == Mode::ZIPLINE && ziplinePhase_ >= ZIPLINEPhase::Blocks && ziplineBlocks() &&
        ziplineSimplified()) {
        // ── The blocks, and the mesh on them ────────────────────────────────
        //
        // What is not a block is shaded and outlined in grey, so that in
        // every one of the three phases the shading is exactly the part of the
        // model with no block on it -- the coverage the console prints, as a
        // picture -- and a red disk sits on every T-junction that put a
        // component there.
        viewer::drawAxis(view_);
        const QuadLayout &sl = ziplineSimplified()->getLayout();
        const auto &faces = sl.getFaces();
        const auto &arcs = sl.getArcs();
        const std::vector<char> &isBlock = ziplineBlocks()->faceIsBlock();
        auto drawRefused = [&](float width) {
            viewer::color3f(0.55f, 0.55f, 0.6f);
            viewer::lineWidth(width);
            for (size_t f = 0; f < faces.size(); ++f) {
                if (f < isBlock.size() && isBlock[f]) continue;
                for (const int d : faces[f].darts) {
                    const auto &pts = arcs[QuadLayout::arcOfDart(d)].pts;
                    glBegin(GL_LINE_STRIP);
                    for (const Point &p : pts) glVertex2d(p[0], p[1]);
                    glEnd();
                }
            }
            viewer::lineWidth(1.0f);
        };
        const double sceneDiag = std::hypot(bounds_.maxx - bounds_.minx, bounds_.maxy - bounds_.miny);
        auto fillUncovered = [&]() {
            if (ziplineUncoveredTris_.empty()) return;
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            viewer::color4f(1.0f, 0.55f, 0.5f, 0.35f);
            glBegin(GL_TRIANGLES);
            for (const int t : ziplineUncoveredTris_) {
                const Triangle &tri = mesh_->triangles[t];
                for (int k = 0; k < 3; ++k)
                    glVertex2d(mesh_->vertices[tri[k]][0], mesh_->vertices[tri[k]][1]);
            }
            glEnd();
        };
        auto drawTJunctions = [&]() {
            const double r = std::max(0.25 * avgEdge_, 0.007 * sceneDiag) * view_.zoom;
            for (const Point &p : ziplineBlocks()->tJunctionPoints())
                viewer::drawDisk3D(p, r, 0.95f, 0.1f, 0.1f);
        };

        if (ziplinePhase_ >= ZIPLINEPhase::Mesh && ziplineMesh()) {
            const bool matFill = showMaterialFill_ && ziplineBlocks()->getReport().materials > 1;
            fillUncovered();
            drawRefused(2.5f);
            if (ziplinePhase_ == ZIPLINEPhase::Smoothed && smoothMesh_.has_value())
                viewer::drawQuadMesh(*smoothMesh_, *ziplineMesh(), 1.0f, 2.5f, matFill);
            else
                viewer::drawQuadMesh(*ziplineMesh(), 1.0f, 2.5f, matFill);
        } else {
            // The same picture UMBER's and ATLAS's block phases draw: the
            // blocks in light blue with green macrovertices, over the layout
            // in grey, so the grey showing through is what has no block.
            // The triangulation light, not in the wireframe colour the
            // earlier phases use: at this phase it is background, and the
            // blocks and the uncovered fill are what the picture is for.
            viewer::drawMeshOverlay(*mesh_, 0.7f, 0.7f, 0.72f, 0.35f, 1.0f);
            viewer::drawBoundaryEdges(*mesh_);
            fillUncovered();
            viewer::drawQuadLayoutArcs(sl, 2.0f, 0.55f, 0.55f, 0.6f);
            viewer::drawBlockDecomposition(ziplineBlocks()->decomposition(),
                                           0.12 * avgEdge_ * view_.zoom, 4.0f);
        }
        if (showInterfaces_ && interfaceNetwork() && ziplinePhase_ == ZIPLINEPhase::Blocks)
            viewer::drawInterfaceNetwork(*interfaceNetwork(), 0.3 * avgEdge_, 2.0f);
        drawTJunctions();
    } else if (mode_ == Mode::ZIPLINE) {
        viewer::drawAxis(view_);
        // The interface network under everything else: the input the field
        // was aligned to and the layout has to keep, so the question every
        // picture after it answers is whether the curves drawn over it follow
        // it. 'i' hides it.
        const bool showNetwork = showInterfaces_ && interfaceNetwork() &&
                                 ziplinePhase_ < ZIPLINEPhase::Quantized;
        // The quantized grid is the payoff of this whole mode, and the
        // triangulation underneath only buries it -- the medial axis mode's
        // Quantized phase drops the mesh the same way.
        if (ziplinePhase_ < ZIPLINEPhase::Quantized) viewer::drawMesh(*mesh_);
        if (showNetwork) viewer::drawInterfaceNetwork(*interfaceNetwork(), 0.3 * avgEdge_, 3.0f);
        // The crossfield answers "why did the separatrices go where they
        // went"; once the quantized grid is up that question is moot and
        // the crosses only clutter the block decomposition it took the
        // rest of the pipeline to produce.
        if (ziplinePhase_ >= ZIPLINEPhase::CrossField && ziplinePhase_ < ZIPLINEPhase::Quantized &&
            ziplineField()) {
            if (ziplinePhase_ >= ZIPLINEPhase::Stepping && zipline_->getStatus().fieldSteps > 0)
                viewer::drawVertexCrossFieldUK(*mesh_, *ziplineField(), scale_);
            else
                viewer::drawVertexCrossField(*mesh_, *ziplineField(), scale_);

            if (ziplinePhase_ < ZIPLINEPhase::Separatrices) {
                double ballRadius = 0.5 * avgEdge_;
                for (const auto &sig : ziplineField()->singularTriangles) {
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
        if (ziplinePhase_ >= ZIPLINEPhase::Quantized && ziplineQuant_.has_value() && ziplineQuant_->ok) {
            // The integer edge lengths as the quad grid they prescribe, same
            // as the medial axis mode's Quantized phase: transfinite curves
            // through each component's four sides -- plus a yellow disk on
            // every grid vertex, interior cell corners included, since after
            // quantization every crossing of two grid lines is a vertex of
            // the finished decomposition. That covers the block-structure
            // T-junctions too: with every edge quantized >= 1 the grids
            // meeting at one weld flush and it is a regular vertex of the
            // result, not a defect to flag.
            //
            // Sized by the model, not by the mesh: vertices are as far
            // apart as cells are, and cells scale with the model, so a
            // radius in mesh edges vanishes on any finely meshed model.
            // view_.zoom is 1.0 at fit and shrinks as the view zooms in, so
            // scaling by it keeps the markers the same size on screen.
            const double sceneDiag = std::hypot(bounds_.maxx - bounds_.minx,
                                                bounds_.maxy - bounds_.miny);
            const double nodeRadius = 0.010 * sceneDiag * view_.zoom;
            viewer::drawQuantizedLayout(*ziplineQuant_, nodeRadius);

            // Red only for the junctions quantization could NOT resolve:
            // an incident edge forced to zero, or a component the
            // conversion skipped. Found structurally, since the stored
            // node kinds go stale once chord collapses merge nodes.
            for (const int n : hangingTJunctions()) {
                viewer::drawDisk3D(ziplineSimplified()->getLayout().getNodes()[n].pos,
                                   1.3 * nodeRadius, 0.95f, 0.1f, 0.1f);
            }
        } else if (ziplinePhase_ >= ZIPLINEPhase::Simplified && ziplineSimplified()) {
            viewer::drawQuadLayoutArcs(ziplineSimplified()->getLayout(), 3.0f, 0.95f, 0.2f, 0.2f);
            // view_.zoom is 1.0 at fit and shrinks as the view zooms in, so
            // scaling the radius by it makes the markers shrink along with it
            // rather than staying a fixed size in mesh space and so covering
            // more and more of the screen as you zoom in on a component.
            viewer::drawQuadLayoutNodes(ziplineSimplified()->getLayout(), 0.12 * avgEdge_ * view_.zoom);
        } else if (ziplinePhase_ >= ZIPLINEPhase::Layout && ziplineLayout()) {
            viewer::drawQuadLayoutArcs(*ziplineLayout(), 3.0f, 0.95f, 0.2f, 0.2f);
        } else if (ziplinePhase_ >= ZIPLINEPhase::Separatrices && ziplineTrace()) {
            viewer::lineWidth(3.0f);
            for (const auto &sep : ziplineTrace()->separatrices) {
                if (sep.path.size() < 2) continue;
                if (sep.active)
                    viewer::color3f(0.95f, 0.1f, 0.1f);
                else
                    viewer::color3f(0.1f, 0.9f, 0.2f);
                glBegin(GL_LINE_STRIP);
                for (const auto &tp : sep.path)
                    glVertex2d(tp.global_pos[0], tp.global_pos[1]);
                glEnd();
            }
            viewer::lineWidth(1.0f);
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
    } else if (mode_ == Mode::UMBER && umberPhase_ >= UMBERPhase::Mesh &&
               umberMesh_.has_value()) {
        // ── The mesh, and after TMOP the same mesh drawn the same way ───────
        //
        // The pipelines' Mesh and Smoothed pictures exactly: one routine, the
        // walls taken off the blocks either way, so that the only difference
        // between the two frames is where the nodes are.
        viewer::drawAxis(view_);
        const bool matFill = showMaterialFill_ && umberDecompReport_.materials > 1;
        if (umberPhase_ == UMBERPhase::Smoothed && smoothMesh_.has_value())
            viewer::drawQuadMesh(*smoothMesh_, *umberMesh_, 1.0f, 2.5f, matFill);
        else
            viewer::drawQuadMesh(*umberMesh_, 1.0f, 2.5f, matFill);
    } else if (mode_ == Mode::UMBER && umberPhase_ >= UMBERPhase::Decomposition &&
               blockLayout_.has_value()) {
        // ── The structure the tracing left, against what the collapse made ──
        //
        // Both over the one mesh at the one scale, so that a chord that went is
        // a grey line with nothing on it and everything else is drawn over
        // grey. The parameter domain is dropped here: the operation happens in
        // the model, and the picture that answers "which blocks did that
        // remove" is this one.
        //
        // What stands afterwards is drawn as the shared BlockDecomposition and
        // not as a QuadLayout, by the routine ATLAS's Blocks phase and the
        // pipelines' Patches phase draw theirs with -- light-blue sides, green
        // macrovertices. The three methods reach a decomposition by three
        // unrelated routes and the picture of one should not say which, and it
        // is also the structure the Mesh phase actually meshes, which the
        // layout underneath it is not: a face the collapse left with a split
        // side is in the layout and is not a block.
        viewer::drawAxis(view_);
        viewer::drawMesh(*mesh_);
        if (umberCut_.has_value() && !umberCut_->getCutEdges().empty())
            viewer::drawEdgeSetOnMesh(*mesh_, umberCut_->getCutEdges(), 1.0f, 0.2f, 0.9f, 2.0f);

        // The interface network over it: an input rather than a result, so the
        // question the picture answers is whether the blue sides still run
        // along these curves.
        if (showInterfaces_ && interfaces_.has_value() && interfaces_->multiMaterial())
            viewer::drawInterfaceNetwork(*interfaces_, 0.4 * avgEdge_, 3.0f);

        viewer::drawQuadLayoutArcs(blockLayout_->getLayout(), 2.0f, 0.45f, 0.45f, 0.5f);
        if (umberDecomp_.has_value()) {
            // view_.zoom is 1.0 at fit and shrinks as the view zooms in, so
            // scaling the radius by it keeps the markers the same size on
            // screen instead of swallowing a block once you zoom in on one.
            viewer::drawBlockDecomposition(*umberDecomp_, 0.12 * avgEdge_ * view_.zoom, 4.0f);
        } else {
            viewer::drawQuadLayoutArcs(blockLayout_->getLayout(), 4.0f, 0.95f, 0.25f, 0.2f);
            viewer::drawQuadLayoutNodes(blockLayout_->getLayout(), 0.12 * avgEdge_ * view_.zoom);
        }
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
            // The three DualMBO stages UMBER shares with both layout pipelines:
            // the input field and the singularities Sec. 4.2 is about to move.
            if (umberPhase_ >= UMBERPhase::CrossField && dualMBOField_.has_value()) {
                viewer::drawTriangleCrossField(*mesh_, *dualMBOField_, scale_);
                double ballRadius = 0.5 * avgEdge_;
                for (const auto &[vertIdx, crossIndex] : dualMBOField_->singularVertices) {
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
    } else if (inPipeline() &&
               (pipePhase_ == PipelinePhase::Layout ||
                pipePhase_ == PipelinePhase::Separatrices) &&
               immersion_.has_value()) {
        // ── Split-screen: left = the model, right = Omega in the plane ───────
        //
        // A map is only readable against the thing it is a map of, so it gets
        // half the window and the model keeps the other. The right half
        // is Psi once the continuation has run and psi_R before it -- or after
        // it, on 'p', which is the only way to see what Stage 6 actually did.
        // At the separatrix phase the map is not a choice: the curves were
        // marched over Psi, and drawing them over psi_R would be drawing them
        // over a map they are not the integral curves of.
        const bool wantPsiR = showPsiR_ && pipePhase_ == PipelinePhase::Layout;
        const std::vector<Point> &uv =
            (meridianLayout_.has_value() && !wantPsiR) ? meridianLayout_->getUV() : psiR_;

        const int w = fbw();
        const int halfW = w / 2;

        applyHalfOrtho(0, halfW, view_);
        {
            viewer::ViewState leftVs = view_;
            leftVs.fbw = halfW;
            leftVs.fbh = fbh();
            viewer::drawAxis(leftVs);
        }
        renderMERIDIANModel();
        // The pullback of the curves, over the model they are a layout of --
        // the left half of Fig. 9. Barycentric coordinates are the same numbers
        // in both worlds, so this is the same curve as the right half and not a
        // second trace of it.
        if (separatrices_.has_value())
            viewer::drawSeparatrices(*separatrices_, Separatrices::Space::Model,
                                     0.35 * avgEdge_, 2.5f);

        applyHalfOrtho(halfW, w - halfW, uvView_);
        {
            // Sized by the image rather than by the mesh, so the cones stay
            // visible on a finely triangulated model, and scaled by the zoom so
            // they keep their size on screen.
            const double diag = std::hypot(uvView_.baseW, uvView_.baseH);
            viewer::drawLayoutUV(*immersion_,
                                 meridianLabels_.has_value() ? &*meridianLabels_ : nullptr,
                                 uv, 0.010 * diag * uvView_.zoom);
            // On Psi every segment of these is axis-parallel -- that is what
            // being an integral curve of the map means -- and every gap in one
            // is a seam crossing, where Q4 moves the image to the far bank of
            // the cut while the pullback on the left walks straight on.
            if (separatrices_.has_value())
                viewer::drawSeparatrices(*separatrices_, Separatrices::Space::Image,
                                         0.007 * diag * uvView_.zoom, 2.5f);
        }

        drawSplitDivider(halfW);
    } else if (mode_ == Mode::TORSION && pipePhase_ == PipelinePhase::Metric &&
               immersion_.has_value()) {
        // ── Split-screen: left = the model, right = psi_0 ────────────────────
        //
        // The same two halves the Layout phase gets, one phase earlier, because
        // on this route Stage 4 has already happened: there is no metric
        // standing between the field and a map. What the right half is being
        // judged on is not what it looks like but whether it inverted anything,
        // which drawLayoutUV fills in red -- and on 'p' it is the Stage 4F
        // least-squares map instead, where those red faces are the entire cost
        // of substituting an integration for a flow.
        //
        // No labels are passed, because Stage 5 has not run: the boundary comes
        // out grey rather than in Gamma_u / Gamma_v, which is honest about
        // there being no labelling yet to colour it by.
        const std::vector<Point> &uv =
            (showIntegrated_ && !integratedMap_.empty()) ? integratedMap_
                                                         : immersion_->getUV();

        const int w = fbw();
        const int halfW = w / 2;

        applyHalfOrtho(0, halfW, view_);
        {
            viewer::ViewState leftVs = view_;
            leftVs.fbw = halfW;
            leftVs.fbh = fbh();
            viewer::drawAxis(leftVs);
        }
        renderMERIDIANModel();

        applyHalfOrtho(halfW, w - halfW, uvView_);
        {
            const double diag = std::hypot(uvView_.baseW, uvView_.baseH);
            viewer::drawLayoutUV(*immersion_, nullptr, uv, 0.010 * diag * uvView_.zoom);
        }

        drawSplitDivider(halfW);
    } else if (mode_ == Mode::MERIDIAN && pipePhase_ == PipelinePhase::Metric &&
               !coneFans_.empty()) {
        // ── Split-screen: left = the flat metric on the model, right = the
        //    cone angles it was driven to ─────────────────────────────────────
        //
        // Neither half says what the metric is on its own. The left shows where
        // it differs from the one the model came with, which is the whole of
        // what the flow changed; the right shows the thing the flow was for,
        // which is invisible on the model because the triangles on screen still
        // carry the angles they were built with.
        const int w = fbw();
        const int halfW = w / 2;

        applyHalfOrtho(0, halfW, view_);
        {
            viewer::ViewState leftVs = view_;
            leftVs.fbw = halfW;
            leftVs.fbh = fbh();
            viewer::drawAxis(leftVs);
        }
        renderMERIDIANModel();

        applyHalfOrtho(halfW, w - halfW, uvView_);
        viewer::drawConeFans(coneFans_);

        drawSplitDivider(halfW);
    } else if (inPipeline()) {
        viewer::drawAxis(view_);
        if (pipePhase_ < PipelinePhase::Cones || !cones_.has_value()) {
            // The three DualMBO stages both pipelines start from: the field whose
            // holonomy Sec. 3.1 is about to turn into cone indices, and the
            // interior singularities it already found.
            if (showMaterialFill_ && interfaces_.has_value() && interfaces_->multiMaterial())
                viewer::drawMaterialFill(*mesh_, 0.28f);
            viewer::drawMesh(*mesh_);
            // Stage 0b needs neither, so its network is already there to see.
            if (showInterfaces_ && interfaces_.has_value())
                viewer::drawInterfaceNetwork(*interfaces_, 0.4 * avgEdge_, 3.0f);
            // And Stage 0c's circles over the holes it left, so that the holes
            // read as something taken out on purpose rather than as a gap in
            // the input.
            viewer::drawInclusionCircles(inclusions_, 1.5f);
            if (pipePhase_ >= PipelinePhase::CrossField && dualMBOField_.has_value()) {
                viewer::drawTriangleCrossField(*mesh_, *dualMBOField_, scale_);
                const double ballRadius = 0.5 * avgEdge_;
                for (const auto &[vertIdx, crossIndex] : dualMBOField_->singularVertices) {
                    if (vertIdx < 0 || vertIdx >= static_cast<int>(mesh_->vertices.size())) continue;
                    const Point &c = mesh_->vertices[vertIdx];
                    if (crossIndex > 0)
                        viewer::drawDisk3D(c, ballRadius, 0.2f, 0.2f, 0.95f);
                    else
                        viewer::drawDisk3D(c, ballRadius, 0.95f, 0.2f, 0.2f);
                }
            }
        } else {
            renderMERIDIANModel();
            // The blocks over the model, once Stages 8 and 9 have produced
            // any: this is the whole of what the phase is for, and it is
            // drawn last so no wireframe edge crosses a patch side.
            if (pipePhase_ >= PipelinePhase::Mesh && quadMesh_.has_value()) {
                // Stage 10 draws its own block walls off the blocks it built,
                // so drawing the BlockDecomposition here too would only lay a
                // second, differently sourced copy of them over the first. A
                // face Stage 10 could not mesh has no wall here, which is the
                // point: the blank is where the mesh is not.
                //
                // Where Stage 11 ran it is the merged mesh that is drawn, not
                // Stage 10's: the templates are elements of that one mesh, and
                // an inclusion the fill refused stays a hole here, which is
                // exactly what the report says of it.
                const bool matFill = showMaterialFill_ && interfaces_.has_value() &&
                                     interfaces_->multiMaterial();
                // Stage 12, once it has run, is the same picture drawn by the
                // same routine off the same block walls, with the nodes where
                // TMOP left them: the phases are meant to be comparable by
                // flipping between them, and a second way of drawing would put
                // its own differences into that comparison. Until it has run --
                // dialog cancelled, or Stage 10 refused -- the unsmoothed mesh
                // stays on screen rather than the window going blank.
                if (pipePhase_ == PipelinePhase::Smoothed && smoothMesh_.has_value())
                    viewer::drawQuadMesh(*smoothMesh_, &*quadMesh_,
                                         diskFill_.has_value() ? &*diskFill_ : nullptr,
                                         1.0f, 2.5f, matFill);
                else if (diskFill_.has_value())
                    viewer::drawQuadMesh(*diskFill_, &*quadMesh_, 1.0f, 2.5f, matFill);
                else
                    viewer::drawQuadMesh(*quadMesh_, 1.0f, 2.5f, matFill);
            } else if (pipePhase_ >= PipelinePhase::Patches && arrangement_.has_value()) {
                const std::string src = (mode_ == Mode::TORSION) ? "TORSION" : "MERIDIAN";
                viewer::drawBlockDecomposition(
                    arrangement_->blockDecomposition(splines_.has_value() ? &*splines_ : nullptr,
                                                     16, src),
                    0.30 * avgEdge_, 3.0f);
            }
            // Over the finished picture, and this is where it earns its place:
            // an element or a patch side that crosses one of these curves
            // rather than running along it is exactly the defect the whole
            // multi-material path exists to prevent, and nothing else in the
            // frame shows it.
            if (showInterfaces_ && pipePhase_ >= PipelinePhase::Patches &&
                interfaces_.has_value())
                viewer::drawInterfaceNetwork(*interfaces_, 0.25 * avgEdge_, 2.0f);
            // The excised circles, until Stage 11 has put the O-grids back:
            // after that the elements are the picture and a circle over them
            // says nothing the rim does not.
            if (!diskFill_.has_value()) viewer::drawInclusionCircles(inclusions_, 1.5f);
        }
    } else if (mode_ == Mode::ATLAS) {
        viewer::drawAxis(view_);
        renderATLAS();
    } else if (mode_ == Mode::MedialAxis) {
        viewer::drawAxis(view_);
        // The map and block phases draw over the boundary alone: the
        // triangulation underneath would bury the radii / zone edges.
        if (maPhase_ < MedialAxisPhase::Map) {
            if (delaunayMesh_)
                viewer::drawMesh(*delaunayMesh_);
            else
                viewer::drawMesh(*mesh_);
        }

        if (medialAxis_) {
            if (maPhase_ >= MedialAxisPhase::Quantized && blockQuant_.has_value() &&
                blockQuant_->ok) {
                // The integer edge lengths as the quad grid they prescribe:
                // transfinite curves through each block's four curved
                // sides, edges only. Grids meeting flush across a block
                // wall are the quantization's consistency made visible.
                viewer::drawQuantizedBlocks(*blockQuant_);
            } else if (maPhase_ >= MedialAxisPhase::TMesh &&
                       medialAxisTMesh_.has_value()) {
                // The coarse block decomposition (Fig. 1-left / Fig. 3c):
                // zones filled with a translucent tint of the class colour
                // that picks their quad template, their walls over the fill,
                // and the downsampled axis in full colour on top.
                viewer::drawMedialTMesh(*medialAxisTMesh_, avgEdge_ / 4.0);
                viewer::drawBoundaryEdges(*medialAxis_->mesh);
            } else if (maPhase_ == MedialAxisPhase::Map && medialAxisMap_.has_value()) {
                viewer::drawBoundaryEdges(*medialAxis_->mesh);

                // The map phi drawn as its medial radii (Fig. 7/8 of the
                // paper): each boundary point joined to its image on the
                // axis. Radii first, the axis on top of them. Blue radii are
                // the bijective sections; magenta ones belong to a polar
                // section, where a run of boundary collapses onto a single
                // medial vertex (a discrete axis endpoint).
                viewer::lineWidth(1.0f);
                glBegin(GL_LINES);
                for (const auto &s : medialAxisMap_->spokes()) {
                    if (s.collapsed)
                        viewer::color4f(0.9f, 0.25f, 0.9f, 0.85f);
                    else
                        viewer::color4f(0.35f, 0.55f, 0.95f, 0.55f);
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
        // PolyVector phases 1-4
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

    // Colour legends. These are chrome too -- a key drawn in screen space in
    // the viewer's own font, not part of the model -- so 'h' takes them with
    // the console, the help and the axis. A figure that wants a key gets one
    // typeset in the caption.
    if (showHUD_) drawLegends();

    // Overlay text. On a multi-material model the two display toggles are
    // appended to whatever the phase's own help says, since they apply at every
    // phase of either pipeline and to none of the other modes.
    const std::string meridianKeys =
        ((inPipeline() || mode_ == Mode::ATLAS)
             ? (pipelineBlocked_.empty() ? std::string() : pipelineBlocked_ + "\n")
             : std::string()) +
        ((inPipeline() && interfaces_.has_value() && interfaces_->multiMaterial())
             ? std::string("press 'i' to show/hide the interface network\n"
                           "press 'm' to fill the triangles by material\n")
             : std::string()) +
        ((mode_ == Mode::ATLAS && atlasMultiMaterial_)
             ? std::string("press 'm' to fill the elements by material\n")
             : std::string());

    // UMBER's own, on the same footing: whatever refused, and the material
    // fill where the blocks landed in more than one material.
    const std::string umberKeys =
        (pipelineBlocked_.empty() ? std::string() : pipelineBlocked_ + "\n") +
        ((mode_ == Mode::UMBER && interfaces_.has_value() && interfaces_->multiMaterial())
             ? std::string("press 'i' to show/hide the interface network\n")
             : std::string()) +
        ((umberDecompReport_.materials > 1)
             ? std::string("press 'm' to fill the elements by material\n")
             : std::string());

    // ZIPLINE's: whatever refused, and the material fill.
    const std::string ziplineKeys =
        (pipelineBlocked_.empty() ? std::string() : pipelineBlocked_ + "\n") +
        ((mode_ == Mode::ZIPLINE && interfaceNetwork())
             ? std::string("press 'i' to show/hide the interface network\n")
             : std::string()) +
        ((mode_ == Mode::ZIPLINE && ziplineBlocks() && ziplineBlocks()->getReport().materials > 1)
             ? std::string("press 'm' to fill the elements by material\n")
             : std::string());

    if (mode_ == Mode::Unselected) {
        renderOverlay("press '1' for PolyVector mode\npress '2' for ZIPLINE mode\n"
                      "press '3' for Medial Axis mode\npress '4' for TORSION mode\n"
                      "press '5' for OASIS mode\npress '6' for UMBER mode\n"
                      "press '7' for MERIDIAN mode\npress '8' for ATLAS mode\n"
                      "right-drag to pan, scroll to zoom\n"
                      "press 'r' to restart\npress 'q' to quit");
    } else if (mode_ == Mode::ATLAS && atlasPhase_ == ATLASPhase::Smoothed) {
        renderOverlay((meridianKeys +
                       "press 'c' to smooth again at other TMOP settings\n"
                       "press 'e' to mesh again at another target edge length\n"
                       "press 'r' to restart\npress 'q' to quit").c_str());
    } else if (mode_ == Mode::ATLAS && atlasPhase_ == ATLASPhase::Mesh) {
        renderOverlay((meridianKeys +
                       "press 'c' to smooth the mesh with TMOP\n"
                       "press 'e' to mesh again at another target edge length\n"
                       "press 'r' to restart\npress 'q' to quit").c_str());
    } else if (mode_ == Mode::ATLAS && atlasPhase_ == ATLASPhase::Blocks) {
        renderOverlay((meridianKeys + "press 'c' to mesh the blocks (TFI)\n"
                       "press 'r' to restart\npress 'q' to quit").c_str());
    } else if (mode_ == Mode::ATLAS && atlasPhase_ == ATLASPhase::Search && atlas_ &&
               atlas_->hasCover()) {
        renderOverlay((meridianKeys +
                       (atlasShowInitial_ ? "press 'p' for the carrier the search ended on\n"
                                          : "press 'p' for the carrier the search started from\n") +
                       "press 'c' to continue\npress 'r' to restart\npress 'q' to quit").c_str());
    } else if (mode_ == Mode::OASIS) {
        renderOverlay("press 'c' to change lambda / orientation\n"
                      "press 'r' to restart\npress 'q' to quit");
    } else if (mode_ == Mode::TORSION && pipePhase_ == PipelinePhase::Metric &&
               immersion_.has_value() && !integratedMap_.empty()) {
        renderOverlay((meridianKeys + "press 'p' to swap psi_0 / the Stage 4F map\n"
                                      "press 'r' to restart\npress 'q' to quit").c_str());
    } else if (inPipeline() && pipePhase_ == PipelinePhase::Layout &&
               meridianLayout_.has_value()) {
        renderOverlay((meridianKeys + "press 'p' to swap psi_R / Psi\n"
                                      "press 'r' to restart\npress 'q' to quit").c_str());
    } else if (inPipeline() && pipePhase_ == PipelinePhase::Smoothed) {
        renderOverlay((meridianKeys +
                       "press 'c' to smooth again at other TMOP settings\n"
                       "press 'e' to mesh again at another target edge length\n"
                       "press 'n' to change the connectivity settings and trace again\n"
                       "press 'r' to restart\npress 'q' to quit").c_str());
    } else if (inPipeline() && pipePhase_ == PipelinePhase::Mesh) {
        renderOverlay((meridianKeys +
                       "press 'c' to smooth the mesh with TMOP (Stage 12)\n"
                       "press 'e' to mesh again at another target edge length\n"
                       "press 'n' to change the connectivity settings and trace again\n"
                       "press 'r' to restart\npress 'q' to quit").c_str());
    } else if (inPipeline() && pipePhase_ == PipelinePhase::Patches) {
        renderOverlay((meridianKeys + "press 'c' to mesh the patches (Stage 10)\n"
                       "press 'n' to change the connectivity settings and trace again\n"
                       "press 'r' to restart\npress 'q' to quit").c_str());
    } else if (mode_ == Mode::ZIPLINE && ziplinePhase_ == ZIPLINEPhase::Blocks) {
        renderOverlay((ziplineKeys + "press 'c' to mesh the blocks\n"
                       "press 'r' to restart\npress 'q' to quit").c_str());
    } else if (mode_ == Mode::ZIPLINE && ziplinePhase_ == ZIPLINEPhase::Mesh) {
        renderOverlay((ziplineKeys + "press 'c' to smooth the mesh with TMOP\n"
                       "press 'e' to mesh again at another target edge length\n"
                       "press 'r' to restart\npress 'q' to quit").c_str());
    } else if (mode_ == Mode::ZIPLINE && ziplinePhase_ == ZIPLINEPhase::Smoothed) {
        renderOverlay((ziplineKeys + "press 'c' to smooth again at other TMOP settings\n"
                       "press 'e' to mesh again at another target edge length\n"
                       "press 'r' to restart\npress 'q' to quit").c_str());
    } else if (mode_ == Mode::UMBER && umberPhase_ == UMBERPhase::Decomposition) {
        renderOverlay((umberKeys + "press 'c' to mesh the blocks\n"
                       "press 'r' to restart\npress 'q' to quit").c_str());
    } else if (mode_ == Mode::UMBER && umberPhase_ == UMBERPhase::Mesh) {
        renderOverlay((umberKeys + "press 'c' to smooth the mesh with TMOP\n"
                       "press 'e' to mesh again at another target edge length\n"
                       "press 'r' to restart\npress 'q' to quit").c_str());
    } else if (mode_ == Mode::UMBER && umberPhase_ == UMBERPhase::Smoothed) {
        renderOverlay((umberKeys + "press 'c' to smooth again at other TMOP settings\n"
                       "press 'e' to mesh again at another target edge length\n"
                       "press 'r' to restart\npress 'q' to quit").c_str());
    } else {
        renderOverlay(((mode_ == Mode::ZIPLINE ? ziplineKeys : meridianKeys) +
                       "press 'c' to continue\npress 'r' to restart\n"
                                      "press 'q' to quit").c_str());
    }
}

// ── legends ─────────────────────────────────────────────────────────────────────────

// The screen-space colour keys, drawn before the text overlay so the console
// keeps drawing on top. Called only while the HUD is shown.
void CrossGenWidget::drawLegends() {
    // ATLAS: the corners' key at Stage 1, the defects' at the two carrier
    // phases; none once the blocks replace the carrier.
    if (mode_ == Mode::ATLAS) {
        if (atlasPhase_ == ATLASPhase::Field && atlasField_)
            viewer::drawReferenceFieldLegend(fbw(), fbh());
        else if (atlasPhase_ == ATLASPhase::Domain && atlasDomain_)
            viewer::drawATLASLegend(fbw(), fbh(), true);
        else if ((atlasPhase_ == ATLASPhase::Carrier && atlasCarrier_) ||
                 (atlasPhase_ == ATLASPhase::Search && atlas_ && atlas_->hasCover()))
            viewer::drawATLASLegend(fbw(), fbh(), false);
        return;
    }
    if (mode_ == Mode::OASIS && oasisPhase_ == OASISPhase::Field && oasis_.has_value()) {
        viewer::drawScalarFieldLegend(fbw(), fbh(), -oasisAbsMax_, oasisAbsMax_,
                                      "quasi-eigenfunction");
    }

    // Either pipeline: the node colours of the interface network, above the
    // cone legend, whenever the network is on screen.
    const bool sepLegendShown = (inPipeline() && cones_.has_value() &&
                                 pipePhase_ == PipelinePhase::Separatrices &&
                                 separatrices_.has_value());
    // ZIPLINE draws the network up to the quantization and at its Blocks
    // phase, and not over the mesh; the key goes with it.
    const bool ziplineNetworkShown =
        mode_ == Mode::ZIPLINE && (ziplinePhase_ < ZIPLINEPhase::Quantized || ziplinePhase_ == ZIPLINEPhase::Blocks);
    if ((inPipeline() || mode_ == Mode::UMBER || ziplineNetworkShown) && showInterfaces_ &&
        interfaceNetwork() && interfaceNetwork()->multiMaterial()) {
        viewer::drawInterfaceLegend(fbw(), fbh(), sepLegendShown);
    }

    // Either pipeline: the cone colours everywhere they are drawn, and
    // whichever ramp -- or, on the field route, whichever palette -- the
    // current phase is using under them.
    if (inPipeline() && cones_.has_value()) {
        viewer::drawConeLegend(fbw(), fbh());
        if (pipePhase_ == PipelinePhase::Separatrices && separatrices_.has_value()) {
            viewer::drawSeparatrixLegend(fbw(), fbh());
        } else if (mode_ == Mode::TORSION && pipePhase_ == PipelinePhase::Flow &&
                   frames_.has_value()) {
            viewer::drawCombedFrameLegend(fbw(), fbh(), *frames_);
        } else if (pipePhase_ == PipelinePhase::Flow && ricciU_.size() > 0) {
            viewer::drawScalarFieldLegend(fbw(), fbh(), -ricciUAbsMax_, ricciUAbsMax_,
                                          "conformal factor u, mean removed");
        } else if (pipePhase_ >= PipelinePhase::Metric &&
                   pipePhase_ != PipelinePhase::Separatrices &&
                   pipePhase_ != PipelinePhase::Patches &&
                   !flatMetric_.edges.empty()) {
            viewer::drawScalarFieldLegend(fbw(), fbh(), -flatMetric_.absMax, flatMetric_.absMax,
                                          "log(l_flat / l_input), mean removed");
        }
    }
}

// ── overlay helper ────────────────────────────────────────────────────────────

void CrossGenWidget::renderOverlay(const char *helpText) {
    // The console and the key help are the viewer talking about itself, so a
    // figure does without them; 'h' is the switch. Everything else on screen
    // belongs to the model and is exported as it stands.
    if (!showHUD_) return;
    int w = fbw(), h = fbh();
    console_.draw(w, h, 55.0f);

    // The figure keys go here rather than into each phase's own help string,
    // for the reason the interface toggles are appended to it: they apply at
    // every phase of every mode, so eleven copies of them would be eleven
    // things to keep in step.
    const std::string help =
        std::string(helpText) +
        "\npress 's' to save a PDF figure, 'S' for a 3x PNG\n"
        "press 'h' to hide this text, 'b' for a " +
        (viewer::lightBackground() ? "dark" : "white") + " background";
    viewer::drawTextOverlay(w, h, help.c_str(), 10.0f, 20.0f, 0.8f, 0.8f, 0.8f);
}

// ── figure export ─────────────────────────────────────────────────────────────

QString CrossGenWidget::nextFigurePath(const char *extension) const {
    // <model>_<mode>_<phase>_<nn>.<ext>, lowercased and with the spaces of the
    // phase names squeezed out, so the files sort into the order the pipeline
    // ran and a caption can be read off the name.
    auto slug = [](QString s, int limit) {
        s = s.toLower();
        QString out;
        for (const QChar c : s) {
            if (c.isLetterOrNumber()) out += c;
            else if (!out.isEmpty() && out.back() != QLatin1Char('-')) out += QLatin1Char('-');
        }
        // The phase names carry their section numbers, which is what makes them
        // worth having in the filename; the prose after them is not, past the
        // first few words. Cut at a word boundary so the result still reads.
        if (out.size() > limit) {
            const int cut = out.lastIndexOf(QLatin1Char('-'), limit);
            out.truncate(cut > 0 ? cut : limit);
        }
        while (out.endsWith(QLatin1Char('-'))) out.chop(1);
        return out;
    };

    const char *phase = "";
    switch (mode_) {
    case Mode::TORSION:
    case Mode::MERIDIAN:   phase = pipelinePhaseName(pipePhase_, mode_); break;
    case Mode::ZIPLINE:        phase = ziplinePhaseName(ziplinePhase_);              break;
    case Mode::MedialAxis: phase = medialAxisPhaseName(maPhase_);        break;
    case Mode::OASIS:      phase = oasisPhaseName(oasisPhase_);          break;
    case Mode::UMBER:      phase = umberPhaseName(umberPhase_);          break;
    case Mode::ATLAS:      phase = atlasPhaseName(atlasPhase_);          break;
    case Mode::PolyVector:
    case Mode::Unselected:  phase = phaseName(phase_);                   break;
    }

    QString stage = slug(QString::fromLatin1(modeName(mode_)), 24);
    const QString ph = slug(QString::fromLatin1(phase), 28);
    if (!ph.isEmpty()) stage += QLatin1Char('_') + ph;

    return QDir(figureDir_).filePath(
        QStringLiteral("%1_%2_%3.%4")
            .arg(QString::fromStdString(modelName_))
            .arg(stage)
            .arg(figureCounter_, 2, 10, QLatin1Char('0'))
            .arg(QString::fromLatin1(extension)));
}

void CrossGenWidget::exportPdf() {
    if (!QDir().mkpath(figureDir_)) {
        console_.log("[Figure] could not create " + figureDir_.toStdString());
        return;
    }
    const QString path = nextFigurePath("pdf");

    makeCurrent();

    // Feedback mode does not rasterise, so by the specification it should not
    // need a draw target -- but on this platform's GL (2.1 over Metal) it does:
    // with no complete framebuffer bound it returns the passthrough tokens and
    // silently drops every primitive, which reads as an empty scene rather than
    // as an error. QOpenGLWidget::makeCurrent() binds the widget's own FBO, so
    // the line above is normally enough; this is here because the failure is
    // invisible if it ever stops being.
    if (QOpenGLContext *ctx = QOpenGLContext::currentContext()) {
        GLint bound = 0;
        glGetIntegerv(GL_FRAMEBUFFER_BINDING, &bound);
        if (bound == 0)
            ctx->functions()->glBindFramebuffer(GL_FRAMEBUFFER, defaultFramebufferObject());
    }

    // The capture has to see exactly the state paintGL leaves before it draws:
    // the same projection, the same viewport, depth testing off.
    viewer::applyOrtho(view_);
    glDisable(GL_DEPTH_TEST);

    // The primitive count is not known until the scene has been drawn, and the
    // feedback buffer cannot be grown once the primitives have gone through it,
    // so on overflow there is nothing to do but draw it again into a bigger
    // one. Four attempts covers 4M to 256M floats, which is about 9M triangles.
    std::string svg;
    viewer::CaptureResult result = viewer::CaptureResult::Overflow;
    for (std::size_t capacity = 4u << 20; capacity <= (256u << 20); capacity *= 4) {
        if (!viewer::beginVectorCapture(capacity)) break;
        drawScene();
        result = viewer::buildVectorCapture(svg, fbw(), fbh());
        if (result != viewer::CaptureResult::Overflow) break;
    }

    doneCurrent();

    // The capture is SVG because that is what the feedback stream writes; the
    // file is PDF because that is what a paper includes. Qt paints one onto the
    // other, and both sides are vector, so nothing is rasterised on the way.
    // One page the size of the framebuffer at 72 dpi, which puts a pixel on a
    // point and leaves the figure at the proportions it had on screen.
    if (result == viewer::CaptureResult::Ok) {
        const int W = fbw(), H = fbh();
        QSvgRenderer renderer(QByteArray::fromStdString(svg));
        if (!renderer.isValid()) {
            console_.log("[Figure] the captured scene could not be re-read for the PDF");
            update();
            return;
        }
        QPdfWriter writer(path);
        writer.setResolution(72);
        writer.setPageSize(QPageSize(QSizeF(W, H), QPageSize::Point,
                                     QStringLiteral("figure"),
                                     QPageSize::ExactMatch));
        writer.setPageMargins(QMarginsF(0, 0, 0, 0), QPageLayout::Point);
        writer.setTitle(QFileInfo(path).completeBaseName());
        QPainter painter;
        if (!painter.begin(&writer)) {
            result = viewer::CaptureResult::WriteFailed;
        } else {
            painter.setRenderHint(QPainter::Antialiasing, true);
            renderer.render(&painter, QRectF(0, 0, W, H));
            painter.end();
        }
    }

    switch (result) {
    case viewer::CaptureResult::Ok: {
        const auto &st = viewer::lastCaptureStats();
        std::ostringstream oss;
        oss << "[Figure] " << path.toStdString() << " -- " << st.polygons
            << " faces, " << st.lines << " segments, " << st.points
            << " points in " << st.elements << " vector elements";
        console_.log(oss.str());
        std::cerr << oss.str() << "\n";
        ++figureCounter_;
        break;
    }
    case viewer::CaptureResult::Overflow:
        console_.log("[Figure] scene too large for the feedback buffer; "
                     "press 'S' for a raster figure instead");
        break;
    case viewer::CaptureResult::Empty:
        console_.log("[Figure] nothing was drawn, so nothing was written");
        break;
    case viewer::CaptureResult::WriteFailed:
        console_.log("[Figure] could not write " + path.toStdString());
        break;
    }
    update();
}

void CrossGenWidget::exportPng(int scale) {
    if (!QDir().mkpath(figureDir_)) {
        console_.log("[Figure] could not create " + figureDir_.toStdString());
        return;
    }
    const QString path = nextFigurePath("png");

    makeCurrent();

    const int w = static_cast<int>(width()  * devicePixelRatio()) * scale;
    const int h = static_cast<int>(height() * devicePixelRatio()) * scale;

    QOpenGLFramebufferObjectFormat fmt;
    fmt.setAttachment(QOpenGLFramebufferObject::CombinedDepthStencil);
    fmt.setSamples(4);
    QOpenGLFramebufferObject fbo(w, h, fmt);
    if (!fbo.isValid()) {
        console_.log("[Figure] could not allocate a framebuffer that large; "
                     "try a smaller scale");
        doneCurrent();
        return;
    }

    // Everything the viewer draws in screen space is derived from fbw()/fbh(),
    // and every line width and glyph is in pixels. Scaling both together is
    // what makes the export a bigger picture of the same thing rather than the
    // same picture with hairlines in it.
    exportScale_ = scale;
    viewer::setRenderScale(static_cast<float>(scale));

    fbo.bind();
    viewer::applyOrtho(view_);
    const auto bg = viewer::backgroundColor();
    glClearColor(bg[0], bg[1], bg[2], 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);
    drawScene();
    const QImage image = fbo.toImage();
    fbo.release();

    exportScale_ = 1;
    viewer::setRenderScale(1.0f);

    // QOpenGLFramebufferObject::release() binds framebuffer 0, which is not
    // the one a QOpenGLWidget paints into; restore it before handing the
    // context back or the next frame goes to the window system's buffer.
    // Through Qt's loader rather than the GL header, which on this platform is
    // the 2.1 one and does not declare the FBO entry points unsuffixed.
    if (QOpenGLContext *ctx = QOpenGLContext::currentContext())
        ctx->functions()->glBindFramebuffer(GL_FRAMEBUFFER, defaultFramebufferObject());
    viewer::resizeAndApplyOrtho(view_, fbw(), fbh());
    uvView_.fbw = view_.fbw;
    uvView_.fbh = view_.fbh;

    doneCurrent();

    if (image.isNull() || !image.save(path)) {
        console_.log("[Figure] could not write " + path.toStdString());
    } else {
        std::ostringstream oss;
        oss << "[Figure] " << path.toStdString() << " -- " << w << "x" << h
            << " (" << scale << "x)";
        console_.log(oss.str());
        std::cerr << oss.str() << "\n";
        ++figureCounter_;
    }
    update();
}
