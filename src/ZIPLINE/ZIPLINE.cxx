#include "ZIPLINE/ZIPLINE.hxx"

#include <algorithm>
#include <cmath>
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
    buildConformalField();
    if (last == Stage::Conformal) return true;
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
    if (from <= Stage::Blocks) {
        blocks.reset();
        status.weakChosen = -1;
        status.coverageBeforeRepair = 0.0;
        status.repairRounds = 0;
        status.repairCandidates = 0;
        status.repairEdits.clear();
        status.repairLeftFirstField = false;
    }
    if (from <= Stage::Simplified) simplification.reset();
    if (from <= Stage::Layout) layout.reset();
    if (from <= Stage::Trace) {
        trace.reset();
        status.conformalTraced = false;
    }
    if (from <= Stage::Conformal) {
        conformal.reset();
        pattern.clear();
        absorbed4.clear();
        absorbedTriangles.clear();
        patternRead = false;
        conformalTried = false;
        status.conformalFailure.clear();
    }
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
        field->setWeakInterfaces(options.weakInterfaces);
    }
    field->setWeakBoundary(options.weakBoundary);
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
    if (!conformalTried) buildConformalField();
    discardFrom(Stage::Trace);
    if (hasConformalField()) {
        trace = std::make_unique<SeparatrixTrace>(conformal->crossField(), true, options.trace,
                                                  interfaces.get(), conformal);
        status.conformalTraced = true;
    } else {
        // The first field: Viertel's tracing alone (Options::trace says why).
        SeparatrixTrace::Settings viertel = options.trace;
        viertel.captureAtSingularities = false;
        viertel.stopAtCrossings = true;
        trace = std::make_unique<SeparatrixTrace>(field, true, viertel, interfaces.get());
    }
}

// ---------------------------------------------------------------------------
// Stage 1b. The pattern is read off the first field exactly as a trace of it
// would read it -- FieldTracer with the same settings, singularities in the
// boundary triangles included, each at the point of its triangle where the
// interpolated representation vector vanishes -- so the second field has the
// singularities the first one would have been traced from, and only the field
// round them changes.
//
// Declined, with the reason kept in Status::conformalFailure, on a
// multi-material model (ConformalField's header says why) and wherever
// ConformalField cannot build -- a pattern failing Poincare-Hopf above all,
// which is a singularity the first field dropped or a corner it turned the
// other way, and leaves the Neumann problem without a solution. The trace
// then follows the first field, as before Stage 1b existed.
// ---------------------------------------------------------------------------
void ZIPLINE::buildConformalField() {
    if (!fieldDone) finishField();
    discardFrom(Stage::Conformal);
    conformalTried = true;
    if (!options.secondField) return;
    // With interfaces the pattern is read as the trace of the first field
    // reads it there: a singular triangle at a node of the network is the
    // node's (FieldTracer's absorbAt), and goes to the field as part of that
    // node's corners (ConformalField's header, "Material interfaces").
    const Interfaces *network = (interfaces && options.trace.respectInterfaces) ? interfaces.get() : nullptr;
    std::vector<char> absorbAt;
    if (network) {
        absorbAt.assign(mesh->vertices.size(), 0);
        for (const Interfaces::Node &n : network->nodes())
            if (n.kind != Interfaces::NodeKind::Dangling && n.vertex >= 0) absorbAt[n.vertex] = 1;
        // Held to the interfaces only weakly (Options::weakInterfaces), the
        // first field moves singularities onto them, and a singular triangle
        // with a vertex on an interface is on it as far as the mesh can say:
        // it is taken there, as a node takes one, and that vertex becomes a
        // node below -- a valence-5 point only (FieldTracer's flag 2). A
        // valence-3 point on a smooth interface would be a corner of one
        // quarter on a straight side, and it stays a point of its region
        // (honeycomb: fourteen of them a mesh edge off its junctions, each a
        // flat corner, and the layout unsound).
        if (field->weakInterfaces())
            for (const int e : network->interfaceEdges())
                for (const int v : {mesh->edges[e][0], mesh->edges[e][1]})
                    if (!mesh->isBoundaryVertex[v] && !absorbAt[v]) absorbAt[v] = 2;
    }
    const FieldTracer first(field, true, options.trace.singularitiesAtBoundary, absorbAt);
    for (const Singularity &s : first.getSingularities())
        pattern.push_back(ConformalField::Singular{s.coordinates, s.triangleIndex, s.d});
    if (first.droppedSingularityCount() > 0) {
        pattern.clear();
        status.conformalFailure = "the first field has singular triangles with no ports";
        return;
    }
    absorbed4.assign(mesh->vertices.size(), 0);
    if (network)
        for (const Interfaces::Node &n : network->nodes())
            if (n.kind != Interfaces::NodeKind::Dangling && n.vertex >= 0) absorbed4[n.vertex] = first.absorbedIndexAt(n.vertex);
    absorbedTriangles = first.absorbedTriangles();
    patternRead = true;
    conformal = std::make_shared<ConformalField>(mesh, pattern, conformalOptions(), field.get(), network, absorbed4,
                                                 absorbedTriangles);
    // A weak first field may turn a region's corner where the network has
    // no node (ConformalField::Options::readWeakBoundary): the nodes are
    // made, once -- they stay in the network for every stage after this,
    // Stage 5b's re-solves included -- and the field is built again.
    if (!conformal->valid() && !conformal->offNetworkCorners().empty() && interfaces) {
        status.weakNodesAdded = interfaces->addNodes(conformal->offNetworkCorners());
        if (status.weakNodesAdded > 0)
            conformal = std::make_shared<ConformalField>(mesh, pattern, conformalOptions(), field.get(), network,
                                                         absorbed4, absorbedTriangles);
    }
    if (!conformal->valid()) status.conformalFailure = conformal->getReport().failure;
}

ConformalField::Options ZIPLINE::conformalOptions() const {
    ConformalField::Options o = options.conformal;
    if (field && (field->weakBoundary() || (interfaces && field->weakInterfaces()))) o.readWeakBoundary = true;
    return o;
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
    repairPattern();
    chooseCondition();
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

// ---------------------------------------------------------------------------
// Stage 5b, PatternRepair's search (its header has the edits, and what P2
// does instead). Each candidate edit is taken through Stages 1b to 5 apart
// from the stages above -- an Attempt -- and the best is kept by soundness,
// then coverage, then refused faces, then blocks, once its mesh passes.
//
// Two states are tracked: the one adopted (the stages above, the best seen)
// and the one the search stands at. They differ only after a round that left
// the first field for the second at a loss: on a model whose pattern failed
// Poincare-Hopf the first round's job is to make the corrected field exist at
// all, and what it covers right then is less than what the rounds after it
// make of it (M7_20: 99.98% on the first field, 97.5% once the two valence-3
// points are in, 100% one move later).
// ---------------------------------------------------------------------------
namespace {

struct Score {
    bool sound = false;
    double coverage = 0.0;
    int refused = 0;
    int blocks = 0;
};

// TestZIPLINE's soundness, of a layout and of its simplification: no loose
// end, no separatrix stuck or cycling, no arcs crossing off a node, and the
// faces tiling the model exactly.
Score scoreOf(const SeparatrixTrace &trace, const QuadLayout &layout, const PartitionSimplify &simp,
              const LayoutBlocks &blocks) {
    Score s;
    const QuadLayout::Report &r = layout.getReport();
    const QuadLayout::Report &q = simp.getLayout().getReport();
    int stuck = 0;
    for (const Separatrix &sp : trace.separatrices)
        if (sp.termination_reason == TerminationReason::STUCK ||
            sp.termination_reason == TerminationReason::LIMIT_CYCLE) ++stuck;
    const double area = r.meshArea;
    auto tiles = [&](const QuadLayout::Report &x) {
        return area > 0.0 && std::fabs(x.totalArea / area - 1.0) < 1e-9 && x.faces > 0 &&
               x.danglingEnds == 0 && x.arcCrossings == 0;
    };
    s.sound = stuck == 0 && tiles(r) && tiles(q);
    s.coverage = blocks.coverage();
    const LayoutBlocks::Report &b = blocks.getReport();
    s.refused = b.faces - b.blocks;
    s.blocks = b.blocks;
    return s;
}

// A real gain: more than minGain of the model covered, or as much with fewer
// refused faces.
bool improves(const Score &a, const Score &b, double minGain) {
    if (!a.sound) return false;
    if (!b.sound) return true;
    if (a.coverage > b.coverage + minGain) return true;
    if (a.coverage < b.coverage - minGain) return false;
    return a.refused < b.refused;
}

// The order the candidates of a round are tried in for the mesh check.
bool better(const Score &a, const Score &b) {
    if (a.sound != b.sound) return a.sound;
    if (std::fabs(a.coverage - b.coverage) > 1e-9) return a.coverage > b.coverage;
    if (a.refused != b.refused) return a.refused < b.refused;
    return a.blocks < b.blocks;
}

} // namespace

struct ZIPLINE::Attempt {
    std::vector<ConformalField::Singular> pattern;
    std::shared_ptr<ConformalField> conformal;
    bool second = false;   // traced on the corrected field
    std::unique_ptr<SeparatrixTrace> trace;
    std::unique_ptr<QuadLayout> layout;
    std::unique_ptr<PartitionSimplify> simplification;
    std::unique_ptr<LayoutBlocks> blocks;
    Score score;
    std::vector<std::string> edits;   // every edit that led here
    bool meshed = false;
    double worst = 0.0;    // after checkSweeps of TMOP
    int inverted = 0;
};

ZIPLINE::Attempt ZIPLINE::attempt(std::vector<ConformalField::Singular> p) const {
    Attempt a;
    a.pattern = std::move(p);
    const Interfaces *network = (interfaces && options.trace.respectInterfaces) ? interfaces.get() : nullptr;
    a.conformal = std::make_shared<ConformalField>(mesh, a.pattern, conformalOptions(), field.get(), network,
                                                   absorbed4, absorbedTriangles);
    a.second = a.conformal->valid();
    // An edit is kept only on the corrected field (repairPattern), so one
    // that leaves the pattern without it is not traced at all.
    if (!a.second) return a;
    a.trace = std::make_unique<SeparatrixTrace>(a.conformal->crossField(), true, options.trace, interfaces.get(),
                                                a.conformal);
    a.trace->run();
    a.layout = std::make_unique<QuadLayout>(*a.trace);
    a.layout->build();
    a.simplification = std::make_unique<PartitionSimplify>(*a.layout, options.simplification);
    if (options.simplify) a.simplification->run();
    // Judged on the counts alone (scoreOf); the geometry is built for the
    // few that are meshed or kept (repairPattern).
    LayoutBlocks::Options lb;
    lb.deferGeometry = true;
    a.blocks = std::make_unique<LayoutBlocks>(a.simplification->getLayout(), *mesh, lb);
    a.score = scoreOf(*a.trace, *a.layout, *a.simplification, *a.blocks);
    return a;
}

// The mesh Stages 6 and 7 would make of the attempt's blocks, smoothed for
// PatternRepair::Options::checkSweeps sweeps only: enough to tell an element
// TMOP can open from one it cannot.
void ZIPLINE::checkMesh(const LayoutBlocks &b, double &worst, int &inverted) const {
    worst = 0.0;
    inverted = 0;
    if (b.blocks().empty()) return;
    const BlockQuadMesh grid(b.decomposition(), options.mesh);
    if (grid.getReport().quads <= 0) return;
    mesh::QuadMesh q = mesh::QuadMesh::from(grid);
    mesh::TMOP::Options t = options.tmop;
    t.maxSweeps = std::min(t.maxSweeps, options.repair.checkSweeps);
    mesh::TMOP smoother(q, t);
    smoother.run();
    worst = smoother.getReport().minScaledJacobianAfter;
    inverted = smoother.getReport().invertedAfter;
}

void ZIPLINE::adopt(Attempt &a) {
    // Latest first: each stage holds on to the one before it.
    blocks.reset();
    simplification.reset();
    layout.reset();
    trace.reset();
    conformal.reset();
    pattern = std::move(a.pattern);
    conformal = std::move(a.conformal);
    trace = std::move(a.trace);
    layout = std::move(a.layout);
    simplification = std::move(a.simplification);
    blocks = std::move(a.blocks);
    status.conformalTraced = a.second;
    status.conformalFailure = a.second ? std::string() : conformal->getReport().failure;
    status.repairEdits = a.edits;
}

// ---------------------------------------------------------------------------
// Options::weakChoices: the first field's condition, chosen per model.
//
// The weak conditions (CrossField::setWeakBoundary) are not better or worse
// than the strong ones everywhere -- the ring a strong condition holds round
// an inclusion is a layer of thin blocks on one model and the right O-grid on
// the next (docs/viertel_2019.md, Eighth pass) -- so each pair is taken through
// Stages 1 to 5b on its own, in a ZIPLINE of its own, and the finished runs,
// this one's included, are scored as Stage 5b scores a candidate edit; the
// best replaces this one's stages wholesale, the network it may have added
// nodes to included. Finished, because what a run covers before Stage 5b
// says little about what it covers after: artery's strong run goes from 88%
// to 100% and its weak one from 99.4% to 99.9%, salt_dome's strong one stays
// at 99.3% and its weak one goes from 97.2% to 99.98%, and chosen before the
// repair both came out the wrong way round. Everything a stage holds is
// swapped, not copied, so the pointers each holds into the stage before stay
// good, and the run that lost takes this one's old stages with it.
// ---------------------------------------------------------------------------
void ZIPLINE::chooseCondition() {
    if (options.weakChoices.empty() || !blocks || !trace || !layout || !simplification) return;
    blocks->completeGeometry();   // Stage 5b may have kept a candidate's, built without it
    Score best = scoreOf(*trace, *layout, *simplification, *blocks);
    double bestWorst = 0.0;
    int bestInverted = 0;
    checkMesh(*blocks, bestWorst, bestInverted);
    const double minGain = options.repair.minGain;
    auto beats = [&](const Score &a, double wa, int ia) {
        if (improves(a, best, minGain)) return true;
        if (improves(best, a, minGain) || !a.sound) return false;
        if (ia != bestInverted) return ia < bestInverted;
        return wa > bestWorst + options.chooseMargin;
    };
    std::unique_ptr<ZIPLINE> winner;
    int chosen = -1;
    for (size_t i = 0; i < options.weakChoices.size(); ++i) {
        const auto &c = options.weakChoices[i];
        // A weak interface means nothing without interfaces.
        if (!(c[0] > 0.0) && !(c[1] > 0.0 && interfaces)) continue;
        Options o = options;
        o.weakChoices.clear();
        o.weakBoundary = c[0];
        o.weakInterfaces = c[1];
        auto alt = std::make_unique<ZIPLINE>(mesh, o);
        alt->run(Stage::Blocks);
        if (!alt->blocks || !alt->trace || !alt->layout || !alt->simplification) continue;
        alt->blocks->completeGeometry();
        const Score s = scoreOf(*alt->trace, *alt->layout, *alt->simplification, *alt->blocks);
        double w = 0.0;
        int inv = 0;
        alt->checkMesh(*alt->blocks, w, inv);
        if (!beats(s, w, inv)) continue;
        best = s;
        bestWorst = w;
        bestInverted = inv;
        winner = std::move(alt);
        chosen = static_cast<int>(i);
    }
    if (!winner) return;
    std::swap(interfaces, winner->interfaces);
    std::swap(field, winner->field);
    std::swap(pattern, winner->pattern);
    std::swap(patternRead, winner->patternRead);
    std::swap(absorbed4, winner->absorbed4);
    std::swap(absorbedTriangles, winner->absorbedTriangles);
    std::swap(conformal, winner->conformal);
    std::swap(trace, winner->trace);
    std::swap(layout, winner->layout);
    std::swap(simplification, winner->simplification);
    std::swap(blocks, winner->blocks);
    options.weakBoundary = winner->options.weakBoundary;
    options.weakInterfaces = winner->options.weakInterfaces;
    status.fieldSteps = winner->status.fieldSteps;
    status.fieldConverged = winner->status.fieldConverged;
    status.fieldError = winner->status.fieldError;
    status.conformalTraced = winner->status.conformalTraced;
    status.conformalFailure = winner->status.conformalFailure;
    status.weakNodesAdded = winner->status.weakNodesAdded;
    status.coverageBeforeRepair = winner->status.coverageBeforeRepair;
    status.repairRounds = winner->status.repairRounds;
    status.repairCandidates = winner->status.repairCandidates;
    status.repairEdits = winner->status.repairEdits;
    status.repairLeftFirstField = winner->status.repairLeftFirstField;
    status.weakChosen = chosen;
}

void ZIPLINE::repairPattern() {
    const PatternRepair::Options &R = options.repair;
    if (R.rounds <= 0 || !blocks || !trace || !layout || !simplification) return;
    // Every edit is an edit of the corrected field's pattern.
    if (!options.secondField || !patternRead) return;
    if (blocks->coverage() >= 1.0 - 1e-12) return;
    status.coverageBeforeRepair = blocks->coverage();
    const bool startedOnFirst = !status.conformalTraced;

    Score adopted = scoreOf(*trace, *layout, *simplification, *blocks);
    // The unrepaired blocks' mesh, for the quality floor; made when first needed.
    double floor = -1.0;
    auto passes = [&](Attempt &a) {
        a.blocks->completeGeometry();
        if (!R.checkMesh) return true;
        if (floor < 0.0) {
            double worst = 0.0;
            int inverted = 0;
            checkMesh(*blocks, worst, inverted);
            floor = std::min(R.qualityFloor, 0.5 * std::max(0.0, worst));
        }
        a.meshed = true;
        checkMesh(*a.blocks, a.worst, a.inverted);
        return a.inverted == 0 && a.worst > floor;
    };

    std::unique_ptr<Attempt> at;   // where the search stands, when not at the adopted state
    for (int r = 0; r < R.rounds; ++r) {
        const std::vector<ConformalField::Singular> &P = at ? at->pattern : pattern;
        const bool second = at ? at->second : status.conformalTraced;
        const ConformalField *cf = at ? at->conformal.get() : conformal.get();
        const ConformalField *declined = (!second && cf && !cf->valid()) ? cf : nullptr;
        const Score from = at ? at->score : adopted;
        if (from.coverage >= 1.0 - 1e-12) break;
        std::vector<PatternRepair::Edit> edits =
            PatternRepair::propose(*mesh, P, declined, at ? at->simplification->getLayout() : simplification->getLayout(),
                                   at ? *at->blocks : *blocks, at ? *at->trace : *trace, R);
        // The budget (PatternRepair::Options::maxCandidates*): the most
        // promising first, the order of propose() kept among equals.
        std::stable_sort(edits.begin(), edits.end(), [](const PatternRepair::Edit &x, const PatternRepair::Edit &y) {
            return x.priority > y.priority;
        });
        const int left = std::min(R.maxCandidatesPerRound, R.maxCandidates - status.repairCandidates);
        if (left <= 0) break;
        if (static_cast<int>(edits.size()) > left) edits.resize(left);
        if (edits.empty()) break;
        ++status.repairRounds;
        status.repairCandidates += static_cast<int>(edits.size());

        // Every candidate through Stages 1b-5, independently: each Attempt
        // owns all it builds and reads only the mesh and the first field, so
        // they run in parallel, and the choice below is made serially from
        // results indexed by candidate -- the same for any thread count.
        // Sec. 4's simplification is most of the cost of one (geom026: 133 of
        // 180 ms; geom012: 720 of 810).
        std::vector<Attempt> tried(edits.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (int e = 0; e < static_cast<int>(edits.size()); ++e)
            tried[e] = attempt(PatternRepair::apply(P, edits[e]));
        // On the corrected field an edit must gain; off it, it must reach it.
        std::vector<std::pair<size_t, Attempt>> kept;
        for (size_t e = 0; e < edits.size(); ++e) {
            Attempt &a = tried[e];
            if (!a.second) continue;
            if (second ? !improves(a.score, from, R.minGain) : !a.score.sound) continue;
            kept.emplace_back(e, std::move(a));
        }
        std::stable_sort(kept.begin(), kept.end(),
                         [](const auto &x, const auto &y) { return better(x.second.score, y.second.score); });
        Attempt *chosen = nullptr;
        size_t pick = 0;
        for (auto &k : kept)
            if (passes(k.second)) { chosen = &k.second; pick = k.first; break; }
        if (!chosen) break;
        chosen->edits = at ? at->edits : status.repairEdits;
        chosen->edits.push_back(edits[pick].what);
        if (improves(chosen->score, adopted, R.minGain)) {
            adopted = chosen->score;
            adopt(*chosen);
            at.reset();
        } else {
            at = std::make_unique<Attempt>(std::move(*chosen));
        }
    }
    status.repairLeftFirstField = startedOnFirst && status.conformalTraced;
}
