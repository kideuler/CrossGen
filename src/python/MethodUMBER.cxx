// Mesh.umber(): UMBER (Wang et al. 2022, src/UMBER) taken to the shared
// BlockDecomposition, then BlockQuadMesh and TMOP -- the chain TestUMBER runs.
//
// UMBER is the one method with no class that runs every stage (MERIDIAN,
// TORSION, ZIPLINE and ATLAS each have one), so its stages are chained here, in
// TestUMBER's order and with TestUMBER's settings:
//
//   Stage 0b  Interfaces       the material interface network (multi-material)
//             DualMBO          the initial cross field, aligned to the interfaces
//             HarmonicCut      no cuts by default: the holes stay holes
//   Sec. 4.2  UMBER            the frame field, Eq. (1)
//   Sec. 4.3  Polysquare       the deformation onto a polysquare
//   Sec. 5    MotorcycleGraph  the iso-lines traced from the corners
//             BlockLayout      the blocks as a graph, then blockDecompositionOf
//
// UMBEROptions below collects what TestUMBER lets a caller change, plus the
// setters the stages have that nothing drove before. Its defaults are copied
// from the stages' own member initialisers (cited on each), so a run with no
// keywords is exactly `TestUMBER <model>`. If a stage's default changes, the
// copy here has to follow it.
#include "python/Methods.hxx"
#include "python/PyTables.hxx"
#include "python/PyUtil.hxx"

#include <stdexcept>

#include "MERIDIAN/Interfaces.hxx"
#include "Parameterization/HarmonicCut.hxx"
#include "UMBER/BlockLayout.hxx"
#include "UMBER/MotorcycleGraph.hxx"
#include "UMBER/Polysquare.hxx"
#include "UMBER/UMBER.hxx"
#include "dualmbo/DualMBO.hxx"

#define OPT(path) CG_OPTION(O, path)
#define REP(path) CG_REPORT(O, path)

namespace pycg {
namespace {

struct UMBEROptions {
    // Stage 0b. TestUMBER always looks for interfaces; a single-material
    // model has none and this costs one pass over the edges.
    bool materialInterfaces = true;

    // The DualMBO field UMBER starts from: TestUMBER's MBO_MAX_STEPS and gamma,
    // which are the pipelines' (not DualMBO's constructor defaults of 100/10).
    int dualMBOMaxSteps = 500;
    double dualMBOGamma = 10.0;

    // HarmonicCut(mesh, openVoids): TestUMBER's UMBER_CUT_VOIDS. Off because a
    // block decomposition wants a common polysquare (see HarmonicCut).
    bool cutVoids = false;

    // UMBER, Sec. 4.2 (UMBER.hxx member initialisers).
    struct FrameOptions {
        int maxIterations = 1000;          // maxIterations
        double alignmentWeight = 0.1;      // alignWeight, w_a
        double regularizationWeight = 5.0; // regWeight, w_r
        double tolerance = 1e-6;           // gradTolerance
        bool cancelDipoles = true;         // cancelDipoles
        int maxCornerQuarters = 1;         // maxCornerQuarters
    } frame;

    // Polysquare, Sec. 4.3 (Polysquare.hxx member initialisers).
    struct PolysquareOptions {
        double interfaceWeight = 0.15;     // interfaceWeight
        double cornerWeight = 10.0;        // corWeight
        bool snapBoundary = true;          // snapBoundaryOn
        double cornerSnapTolerance = 0.05; // cornerSnapTol
        int maxIterations = 1000;          // maxIterations
        double tolerance = 1e-6;           // gradTolerance
    } polysquare;

    // MotorcycleGraph::setArrivalTolerance, in mean mesh edges (arrivalTol,
    // and TestUMBER's default).
    double arrivalTolerance = 0.20;
};

const Table<UMBEROptions> &optionsTable() {
    using O = UMBEROptions;
    static const Table<O> t{
        OPT(materialInterfaces),
        OPT(dualMBOMaxSteps),
        OPT(dualMBOGamma),
        OPT(cutVoids),
        OPT(frame.maxIterations),
        OPT(frame.alignmentWeight),
        OPT(frame.regularizationWeight),
        OPT(frame.tolerance),
        OPT(frame.cancelDipoles),
        OPT(frame.maxCornerQuarters),
        OPT(polysquare.interfaceWeight),
        OPT(polysquare.cornerWeight),
        OPT(polysquare.snapBoundary),
        OPT(polysquare.cornerSnapTolerance),
        OPT(polysquare.maxIterations),
        OPT(polysquare.tolerance),
        OPT(arrivalTolerance),
    };
    return t;
}

const Table<BlockLayout::DecompositionReport> &decompositionTable() {
    using O = BlockLayout::DecompositionReport;
    static const Table<O> t{
        REP(faces),     REP(blocks),           REP(notFourSided), REP(multiArcSides),
        REP(unmatchedSides), REP(materials),   REP(straddlingBlocks), REP(coveredArea),
        REP(totalArea),
    };
    return t;
}

const Table<Polysquare::Report> &polysquareTable() {
    using O = Polysquare::Report;
    static const Table<O> t{
        REP(flips),           REP(minScaledJacobian), REP(avgScaledJacobian),
        REP(meanAlignDeg),    REP(maxAlignDeg),       REP(turns),
        REP(expectedTurns),   REP(shortestRun),       REP(worstSegmentDeg),
        REP(suspectSegments), REP(snapWasExact),      REP(cornersSnapped),
        REP(worstCornerSnap), REP(cornerSnapConflicts), REP(runsAligned),
        REP(worstRunAlign),   REP(interfaceEdges),    REP(meanInterfaceAlignDeg),
        REP(maxInterfaceAlignDeg), REP(interfaceTurns), REP(expectedInterfaceTurns),
        REP(transitionDeg),   REP(lengthRatio),       REP(arap),
        REP(l1),              REP(cor),               REP(iterations),
    };
    return t;
}

const Table<MotorcycleGraph::Report> &motorcycleTable() {
    using O = MotorcycleGraph::Report;
    static const Table<O> t{
        REP(motorcycles),   REP(crossings),  REP(reachedBoundary), REP(ranOut),
        REP(skippedCorners), REP(duplicates), REP(arrivals),       REP(continuations),
        REP(continuationsDropped), REP(nodes), REP(blocks),        REP(mergedSlivers),
        REP(smallestBlock), REP(tracedEdges),
    };
    return t;
}

// What the frame field stage has to say for itself, which UMBER keeps as
// accessors rather than as a Report struct.
struct FrameSummary {
    int iterations = 0;
    double gradientNorm = 0.0;
    double energy = 0.0;
    int internalSingularities = 0;
    int boundarySingularities = 0;
    int unwoundHoles = 0;
    double minFrameNorm = 0.0;
};

const Table<FrameSummary> &frameTable() {
    using O = FrameSummary;
    static const Table<O> t{
        REP(iterations), REP(gradientNorm), REP(energy), REP(internalSingularities),
        REP(boundarySingularities), REP(unwoundHoles), REP(minFrameNorm),
    };
    return t;
}

class UmberMethod final : public Method {
public:
    UmberMethod(const Mesh &model, const UMBEROptions &opts)
        : mesh_(std::make_shared<Mesh>(model)), opts_(opts) {}

    void run() {
        try {
            if (opts_.materialInterfaces) interfaces_ = std::make_unique<Interfaces>(mesh_);
        } catch (const std::exception &) {
            // TestUMBER's answer too: carry on as a single-material model.
            interfaces_.reset();
        }
        const bool multi = interfaces_ && interfaces_->multiMaterial();

        field_ = std::make_unique<DualMBO>(mesh_, opts_.dualMBOMaxSteps, opts_.dualMBOGamma);
        if (multi) field_->setAlignedInteriorEdges(interfaces_->interfaceEdges());
        field_->initialize();
        field_->runMBO();
        field_->computeSingularities();

        cut_ = std::make_unique<HarmonicCut>(mesh_, opts_.cutVoids);

        umber_ = std::make_unique<UMBER>(*field_, *cut_);
        if (multi) umber_->setFeatureEdges(interfaces_->interfaceEdges());
        umber_->setMaxIterations(opts_.frame.maxIterations);
        umber_->setAlignmentWeight(opts_.frame.alignmentWeight);
        umber_->setRegularizationWeight(opts_.frame.regularizationWeight);
        umber_->setTolerance(opts_.frame.tolerance);
        umber_->setCancelDipoles(opts_.frame.cancelDipoles);
        umber_->setMaxCornerQuarters(opts_.frame.maxCornerQuarters);
        umber_->initialize();
        umber_->optimize();
        frame_.iterations = umber_->iterations();
        frame_.gradientNorm = umber_->gradientNorm();
        frame_.energy = umber_->energy().total;
        frame_.internalSingularities = static_cast<int>(umber_->internalSingularities().size());
        frame_.boundarySingularities = static_cast<int>(umber_->boundarySingularities().size());
        frame_.unwoundHoles = umber_->unwoundHoles();
        frame_.minFrameNorm = umber_->minFrameNorm();

        poly_ = std::make_unique<Polysquare>(*umber_, *cut_);
        poly_->setInterfaceWeight(opts_.polysquare.interfaceWeight);
        poly_->setCornerWeight(opts_.polysquare.cornerWeight);
        poly_->setSnapBoundary(opts_.polysquare.snapBoundary);
        poly_->setCornerSnapTolerance(opts_.polysquare.cornerSnapTolerance);
        poly_->setMaxIterations(opts_.polysquare.maxIterations);
        poly_->setTolerance(opts_.polysquare.tolerance);
        poly_->solve();

        graph_ = std::make_unique<MotorcycleGraph>(*poly_);
        graph_->setArrivalTolerance(opts_.arrivalTolerance);
        graph_->build();

        layout_ = std::make_unique<BlockLayout>(*graph_);
        layout_->build();
        decomposition_ = BlockLayout::blockDecompositionOf(layout_->getLayout(), *mesh_, &report_);
    }

    const char *name() const override { return "umber"; }
    const BlockDecomposition &decomposition() const override { return decomposition_; }
    double coverage() const override { return BlockLayout::coverageOf(report_); }

    PyObject *report() const override {
        return toDict({source(decompositionTable(), report_),
                       source(frameTable(), frame_, "frame_"),
                       source(polysquareTable(), poly_->getReport(), "polysquare_"),
                       source(motorcycleTable(), graph_->getReport(), "motorcycles_")});
    }
    PyObject *options() const override { return toDict({source(optionsTable(), opts_)}); }

    int mesh(double h, PyObject *kwargs, MeshOutput &out) const override {
        BlockQuadMesh::Options mo;
        mo.targetEdgeLength = h;
        mesh::QuadMesh::Options no;
        if (applyKwargs(kwargs, "mesh", {target(blockQuadMeshOptionsTable(), mo),
                                         target(nodeOptionsTable(), no)}) < 0)
            return -1;
        std::unique_ptr<BlockQuadMesh> grid;
        const bool ok = compute([&] {
            if (decomposition_.blocks.empty()) throw std::runtime_error("umber: no block to mesh");
            grid = std::make_unique<BlockQuadMesh>(decomposition_, mo);
            if (grid->getReport().quads <= 0) throw std::runtime_error("umber: the mesh has no elements");
            out.mesh = std::make_unique<mesh::QuadMesh>(mesh::QuadMesh::from(*grid, no));
        });
        if (!ok) return -1;
        out.report = toDict({source(blockQuadMeshReportTable(), grid->getReport())});
        out.smoothing = blockGridSmoothing();
        return out.report ? 0 : -1;
    }

    double defaultTarget() const override { return BlockQuadMesh::Options().targetEdgeLength; }

private:
    std::shared_ptr<Mesh> mesh_;
    UMBEROptions opts_;
    // In stage order: each later stage holds references into the ones before
    // it, so they are destroyed in reverse.
    std::unique_ptr<Interfaces> interfaces_;
    std::unique_ptr<DualMBO> field_;
    std::unique_ptr<HarmonicCut> cut_;
    std::unique_ptr<UMBER> umber_;
    std::unique_ptr<Polysquare> poly_;
    std::unique_ptr<MotorcycleGraph> graph_;
    std::unique_ptr<BlockLayout> layout_;
    BlockDecomposition decomposition_;
    BlockLayout::DecompositionReport report_;
    FrameSummary frame_;
};

std::shared_ptr<Method> run(const Mesh &model, PyObject *kwargs) {
    UMBEROptions opts;
    if (applyKwargs(kwargs, "umber", {target(optionsTable(), opts)}) < 0) return nullptr;
    std::shared_ptr<UmberMethod> m;
    if (!compute([&] {
            m = std::make_shared<UmberMethod>(model, opts);
            m->run();
        }))
        return nullptr;
    return m;
}

PyObject *defaults(const std::string &stage) {
    if (stage == "method") return toDict({source(optionsTable(), UMBEROptions())});
    if (stage == "mesh") {
        const BlockQuadMesh::Options mo;
        const mesh::QuadMesh::Options nodes;
        return meshingDefaults(mo.targetEdgeLength, {source(blockQuadMeshOptionsTable(), mo),
                                                     source(nodeOptionsTable(), nodes)});
    }
    if (stage == "smooth") return smoothingDefaults(blockGridSmoothing());
    return unknownStage("umber", stage);
}

void buildTables() {
    optionsTable();
    decompositionTable();
    polysquareTable();
    motorcycleTable();
    frameTable();
}

}  // namespace

const MethodSpec umberSpec = {
    "umber",
    "umber(**options) -> BlockDecomposition\n\n"
    "UMBER (Wang et al. 2022): a frame field (Sec. 4.2), its polysquare\n"
    "deformation (Sec. 4.3) and the motorcycle graph of iso-lines traced from\n"
    "the corners (Sec. 5), read off as blocks, as TestUMBER runs it. Faces the\n"
    "graph leaves with a split side are not blocks: read it by `coverage`.\n\n"
    "Keywords: material_interfaces, dual_mbo_max_steps, dual_mbo_gamma,\n"
    "cut_voids, frame_* (UMBER's setters), polysquare_* (Polysquare's) and\n"
    "arrival_tolerance; crossgen.options('umber') lists them with their\n"
    "defaults. Meshing is BlockQuadMesh on the blocks.",
    run,
    defaults,
    buildTables,
};

}  // namespace pycg
