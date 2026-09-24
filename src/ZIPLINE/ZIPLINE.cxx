#include "ZIPLINE/ZIPLINE.hxx"

#include <limits>
#include <utility>

ZIPLINE::ZIPLINE(std::shared_ptr<Mesh> mesh_) : ZIPLINE(std::move(mesh_), Options()) {}

ZIPLINE::ZIPLINE(std::shared_ptr<Mesh> mesh_, const Options &opts)
    : mesh(std::move(mesh_)), options(opts) {}

ZIPLINE::~ZIPLINE() = default;

bool ZIPLINE::run(Stage last) {
    findInterfaces();
    if (last == Stage::Interfaces) return true;
    startField();
    finishField();
    if (last == Stage::Field) return true;
    startTrace();
    finishTrace();
    if (last == Stage::Trace) return true;
    buildLayout();
    if (last == Stage::Layout) return true;
    simplifyLayout();
    if (last == Stage::Simplified) return true;
    buildBlocks();
    if (last == Stage::Blocks) return true;
    if (!buildMesh()) return false;
    if (last == Stage::Mesh) return true;
    return smooth();
}

void ZIPLINE::discardFrom(Stage from) {
    // Latest first: each stage holds on to the one before it.
    if (from <= Stage::Smoothed) { smoothedMesh.reset(); tmopReport = mesh::TMOP::Report(); }
    if (from <= Stage::Mesh) blockMesh.reset();
    if (from <= Stage::Blocks) blocks.reset();
    if (from <= Stage::Simplified) simplification.reset();
    if (from <= Stage::Layout) layout.reset();
    if (from <= Stage::Trace) trace.reset();
    if (from <= Stage::Field) {
        field.reset();
        fieldDone = false;
        status = Status();
    }
    if (from <= Stage::Interfaces) {
        interfaces.reset();
        interfacesFound = false;
    }
}

// ---------------------------------------------------------------------------
// Stage 0b, as the pipelines run it, with one difference: closed loops are
// not split. A separatrix crossing an inclusion's rim cuts the rim where the
// layout actually turns, and splitting it anywhere else would put a node on
// it that no separatrix ends at (Interfaces::Options::splitCircleLoops has the
// argument). Kept only when there is more than one material, so that every
// stage after it can take "has interfaces" to mean "has a network to honour".
// ---------------------------------------------------------------------------
void ZIPLINE::findInterfaces() {
    discardFrom(Stage::Interfaces);
    interfacesFound = true;
    if (!options.materialInterfaces) return;
    Interfaces::Options io;
    io.splitLoops = false;
    auto found = std::make_unique<Interfaces>(mesh, io);
    if (found->multiMaterial()) interfaces = std::move(found);
}

// ---------------------------------------------------------------------------
// Stage 1. With a network, every interface vertex is Dirichlet at the
// length-weighted mean of exp(4i phi) over its interface edges, which is what
// decouples the materials (CrossField::setAlignedInteriorEdges), and each disk
// inclusion's centre is pinned to u = 1.
// ---------------------------------------------------------------------------
void ZIPLINE::startField() {
    if (!interfacesFound) findInterfaces();
    discardFrom(Stage::Field);
    field = std::make_shared<CrossField>(mesh);
    if (interfaces) {
        field->setAlignedInteriorEdges(interfaces->interfaceEdges());
        field->setPinDiskCenters(options.pinDiskCenters);
    }
    field->initialize(1, options.fieldSeed);
}

void ZIPLINE::advanceField(int steps) {
    if (!field) startField();
    const double tolerance =
        2.0 * static_cast<double>(mesh->vertices.size()) * options.fieldTolerance;
    for (int i = 0; i < steps && !fieldDone; ++i) {
        if (status.fieldSteps >= options.fieldMaxSteps) {
            fieldDone = true;
            break;
        }
        field->step();
        ++status.fieldSteps;
        status.fieldError = field->error;
        if (field->error < tolerance) {
            status.fieldConverged = true;
            fieldDone = true;
        } else if (status.fieldSteps >= options.fieldMaxSteps) {
            fieldDone = true;
        }
    }
}

bool ZIPLINE::stepField(int steps) {
    advanceField(steps);
    field->computeSingularities();
    return fieldDone;
}

void ZIPLINE::finishField() {
    advanceField(std::numeric_limits<int>::max());
    field->computeSingularities();
}

// ---------------------------------------------------------------------------
// Stage 2. The singularity is placed where in its triangle it actually sits,
// not at the barycentre (the `true`): the ports are launched from that point,
// so a barycentre displaces every separatrix leaving it, and the partition
// pays for it at the far end -- over the sixteen models it was measured on it
// cost 40 components and 12 T-junctions, and on geom006, a box with a round
// hole, it was 12 components and no T-junction against 14 and four.
//
// The field is always finished first. Tracing one that is still moving is
// how the viewer used to differ from the driver whenever its MBO animation
// was cut short.
// ---------------------------------------------------------------------------
void ZIPLINE::startTrace() {
    if (!fieldDone) finishField();
    discardFrom(Stage::Trace);
    trace = std::make_unique<SeparatrixTrace>(field, true, options.trace, interfaces.get());
}

bool ZIPLINE::stepTrace() {
    if (!trace) startTrace();
    if (!trace->finishedTracing) trace->stepAndCheck();
    return trace->finishedTracing;
}

void ZIPLINE::finishTrace() {
    if (!trace) startTrace();
    trace->run();
}

void ZIPLINE::buildLayout() {
    if (!traceFinished()) finishTrace();
    discardFrom(Stage::Layout);
    layout = std::make_unique<QuadLayout>(*trace);
    layout->build();
}

void ZIPLINE::simplifyLayout() {
    if (!layout) buildLayout();
    discardFrom(Stage::Simplified);
    simplification = std::make_unique<PartitionSimplify>(*layout, options.simplification);
    if (options.simplify) simplification->run();
}

void ZIPLINE::buildBlocks() {
    if (!simplification) simplifyLayout();
    discardFrom(Stage::Blocks);
    blocks = std::make_unique<LayoutBlocks>(simplification->getLayout(), *mesh);
}

bool ZIPLINE::buildMesh(const BlockQuadMesh::Options &opts) {
    options.mesh = opts;
    return buildMesh();
}

bool ZIPLINE::buildMesh() {
    if (!blocks) buildBlocks();
    discardFrom(Stage::Mesh);
    if (blocks->blocks().empty()) return false;
    blockMesh = std::make_unique<BlockQuadMesh>(blocks->decomposition(), options.mesh);
    return true;
}

// ---------------------------------------------------------------------------
// Stage 7, on a mesh::QuadMesh read off the grid: its boundary and interface
// nodes slide on their curves, and a component with no blocks leaves a rim
// that slides too unless Options::pinOpenSides holds it.
// ---------------------------------------------------------------------------
bool ZIPLINE::smooth() {
    if (!blockMesh && !buildMesh()) return false;
    discardFrom(Stage::Smoothed);
    if (blockMesh->getReport().quads <= 0) return false;
    smoothedMesh = std::make_unique<mesh::QuadMesh>(mesh::QuadMesh::from(*blockMesh));
    if (options.pinOpenSides && !blockMesh->openSideVertices().empty()) {
        for (const int v : blockMesh->openSideVertices()) smoothedMesh->pinVertex(v);
        smoothedMesh->classifyNodes();
        smoothedMesh->buildFeatureCurves();
        smoothedMesh->computeSlideTangents();
    }
    mesh::TMOP smoother(*smoothedMesh, options.tmop);
    smoother.run();
    tmopReport = smoother.getReport();
    return true;
}
