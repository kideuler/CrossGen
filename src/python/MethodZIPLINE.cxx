// Mesh.zipline(): ZIPLINE (src/ZIPLINE/ZIPLINE.hxx), the Viertel et al. IMR 2019
// separatrix partition, run through Stage 5 -- the four-sided components as
// the shared BlockDecomposition. Meshing is Stage 6 (BlockQuadMesh) on those
// blocks, and smoothing is Stage 7's TMOP settings.
//
// The class is run as the driver and viewer run it, with ZIPLINE::Options as it
// stands, so the field seed is ZIPLINE's fixed one and a model gives the same
// layout here as in `TestZIPLINE` and viewer mode 2. Stages 6 and 7 are not
// run through ZIPLINE::buildMesh()/smooth(): those rebuild the object's own
// single mesh, where BlockDecomposition.mesh() may be called any number of
// times at different targets and each result smoothed on its own. What they
// do is repeated here, and it is two lines: a BlockQuadMesh on the blocks, and
// Options::pinOpenSides applied to the mesh::QuadMesh read off it.
#include "python/Methods.hxx"
#include "python/PyTables.hxx"
#include "python/PyUtil.hxx"

#include <stdexcept>

#include "ZIPLINE/ZIPLINE.hxx"

#define OPT(path) CG_OPTION(O, path)
#define REP(path) CG_REPORT(O, path)

namespace pycg {
namespace {

// ZIPLINE::Options, less `mesh` and `tmop`: those are Stages 6 and 7, and are
// BlockDecomposition.mesh()'s and QuadMesh.smooth()'s keywords instead.
const Table<ZIPLINE::Options> &optionsTable() {
    using O = ZIPLINE::Options;
    static const Table<O> t{
        OPT(materialInterfaces),
        OPT(pinDiskCenters),
        OPT(fieldMaxSteps),
        OPT(fieldTolerance),
        OPT(fieldSeed),
        OPT(trace.maxStepsPerSeparatrix),
        OPT(trace.cutAtSingularities),
        OPT(trace.singularityCutRadius),
        OPT(trace.tangentialAngle),
        OPT(trace.boundaryTangentialAngle),
        OPT(trace.boundaryCornerSnap),
        OPT(trace.cornerJoinRadius),
        OPT(trace.evenCornerRays),
        OPT(trace.growByArcLength),
        OPT(trace.stopParallelTangential),
        OPT(trace.resumeAfterTruncation),
        OPT(trace.singularitiesAtBoundary),
        OPT(trace.respectInterfaces),
        OPT(simplify),
        OPT(simplification.zipAngle),
        OPT(simplification.zipLengthOfPatch),
        OPT(simplification.collapseNonZip),
        OPT(simplification.nonZipAngle),
        OPT(simplification.maxDrag),
        OPT(simplification.smoothCollapse),
        OPT(simplification.fixedCorners),
        OPT(simplification.contractRungsWithNodes),
        OPT(simplification.tJunctionsFirst),
        OPT(simplification.sliverSide),
        OPT(simplification.removeSlivers),
        OPT(simplification.extendStems),
        OPT(simplification.stems.minCrossingAngle),
        OPT(simplification.stems.maxCrossings),
        OPT(simplification.stems.maxSteps),
        OPT(simplification.stems.joinRadius),
        OPT(simplification.maxCollapses),
        OPT(pinOpenSides),
    };
    return t;
}

const Table<ZIPLINE::Status> &statusTable() {
    using O = ZIPLINE::Status;
    static const Table<O> t{REP(fieldSteps), REP(fieldConverged), REP(fieldError)};
    return t;
}

const Table<LayoutBlocks::Report> &blocksTable() {
    using O = LayoutBlocks::Report;
    static const Table<O> t{
        REP(faces),           REP(blocks),          REP(notFourSided),
        REP(notDisks),        REP(tJunctionFaces),  REP(tJunctions),
        REP(coveredArea),     REP(totalArea),       REP(tJunctionArea),
        REP(notFourSidedArea), REP(notDiskArea),    REP(macroArcs),
        REP(fittedArcs),      REP(exactArcs),       REP(maxSegmentsUsed),
        REP(unconverged),     REP(maxDeviation),    REP(meanDeviation),
        REP(foldedPatches),   REP(maxCornerGap),    REP(maxSideGap),
        REP(materials),       REP(straddlingBlocks), REP(brepFaces),
        REP(brepEdges),       REP(brepSharedEdges), REP(brepFreeEdges),
        REP(brepFailures),    REP(brepValid),       REP(messages),
    };
    return t;
}

class ZiplineMethod final : public Method {
public:
    ZiplineMethod(const Mesh &model, const ZIPLINE::Options &opts)
        : zipline_(std::make_shared<Mesh>(model), opts) {}

    void run() {
        zipline_.run(ZIPLINE::Stage::Blocks);
        if (!zipline_.hasBlocks()) throw std::runtime_error("zipline: Stage 5 produced no blocks object");
    }

    const char *name() const override { return "zipline"; }
    const BlockDecomposition &decomposition() const override {
        return zipline_.getBlocks().decomposition();
    }
    double coverage() const override { return zipline_.getBlocks().coverage(); }

    PyObject *report() const override {
        return toDict({source(statusTable(), zipline_.getStatus()),
                       source(blocksTable(), zipline_.getBlocks().getReport())});
    }
    PyObject *options() const override {
        return toDict({source(optionsTable(), zipline_.getOptions())});
    }

    int mesh(double h, PyObject *kwargs, MeshOutput &out) const override {
        BlockQuadMesh::Options mo = zipline_.getOptions().mesh;
        mo.targetEdgeLength = h;
        mesh::QuadMesh::Options no;
        if (applyKwargs(kwargs, "mesh", {target(blockQuadMeshOptionsTable(), mo),
                                         target(nodeOptionsTable(), no)}) < 0)
            return -1;
        const bool pin = zipline_.getOptions().pinOpenSides;
        std::unique_ptr<BlockQuadMesh> grid;
        const bool ok = compute([&] {
            if (decomposition().blocks.empty()) throw std::runtime_error("zipline: no block to mesh");
            grid = std::make_unique<BlockQuadMesh>(decomposition(), mo);
            if (grid->getReport().quads <= 0) throw std::runtime_error("zipline: the mesh has no elements");
            out.mesh = std::make_unique<mesh::QuadMesh>(mesh::QuadMesh::from(*grid, no));
            // ZIPLINE::smooth()'s Options::pinOpenSides, on this mesh.
            if (pin && !grid->openSideVertices().empty()) {
                for (const int v : grid->openSideVertices()) out.mesh->pinVertex(v);
                out.mesh->classifyNodes();
                out.mesh->buildFeatureCurves();
                out.mesh->computeSlideTangents();
            }
        });
        if (!ok) return -1;
        out.report = toDict({source(blockQuadMeshReportTable(), grid->getReport())});
        out.smoothing = zipline_.getOptions().tmop;
        return out.report ? 0 : -1;
    }

    double defaultTarget() const override { return zipline_.getOptions().mesh.targetEdgeLength; }

private:
    ZIPLINE zipline_;
};

std::shared_ptr<Method> run(const Mesh &model, PyObject *kwargs) {
    ZIPLINE::Options opts;
    if (applyKwargs(kwargs, "zipline", {target(optionsTable(), opts)}) < 0) return nullptr;
    std::shared_ptr<ZiplineMethod> m;
    if (!compute([&] {
            m = std::make_shared<ZiplineMethod>(model, opts);
            m->run();
        }))
        return nullptr;
    return m;
}

PyObject *defaults(const std::string &stage) {
    const ZIPLINE::Options o;
    if (stage == "method") return toDict({source(optionsTable(), o)});
    if (stage == "mesh") {
        const mesh::QuadMesh::Options nodes;
        return meshingDefaults(o.mesh.targetEdgeLength,
                               {source(blockQuadMeshOptionsTable(), o.mesh),
                                source(nodeOptionsTable(), nodes)});
    }
    if (stage == "smooth") return smoothingDefaults(o.tmop);
    return unknownStage("zipline", stage);
}

void buildTables() {
    optionsTable();
    statusTable();
    blocksTable();
}

}  // namespace

const MethodSpec ziplineSpec = {
    "zipline",
    "zipline(**options) -> BlockDecomposition\n\n"
    "ZIPLINE (Viertel, Osting and Staten, IMR 2019): a cross field, the\n"
    "separatrix partition traced from its singularities, the Sec. 4\n"
    "simplification, and the four-sided components as blocks. Faces with a\n"
    "T-junction on a side are not blocks, so read the result by `coverage`.\n\n"
    "Keywords are the fields of ZIPLINE::Options in snake_case\n"
    "(src/ZIPLINE/ZIPLINE.hxx), e.g. field_seed, simplify,\n"
    "trace_tangential_angle; crossgen.options('zipline') lists them with\n"
    "their defaults. Meshing is BlockQuadMesh on the blocks.",
    run,
    defaults,
    buildTables,
};

}  // namespace pycg
