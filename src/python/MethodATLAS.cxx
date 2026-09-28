// Mesh.atlas(): ATLAS, the square-transport blocking of src/ATLAS, Stages 1-6
// -- the chosen search's certified cover as the shared BlockDecomposition --
// meshed by BlockMesh, which evaluates each block's certified chart rather
// than blending its four sides (TestATLAS --mesh, and viewer mode 8).
#include "python/Methods.hxx"
#include "python/PyTables.hxx"
#include "python/PyUtil.hxx"

#include <stdexcept>

#include "ATLAS/ATLAS.hxx"
#include "ATLAS/BlockMesh.hxx"

#define OPT(path) CG_OPTION(O, path)
#define REP(path) CG_REPORT(O, path)

namespace pycg {

template <>
struct EnumNames<ATLAS::Options::Field::Reference> {
    static constexpr std::pair<ATLAS::Options::Field::Reference, const char *> list[] = {
        {ATLAS::Options::Field::Reference::Harmonic, "harmonic"},
        {ATLAS::Options::Field::Reference::DualMBO, "dual_mbo"},
        {ATLAS::Options::Field::Reference::None, "none"},
    };
};

namespace {

// ATLAS::Options and every scalar of its nested stage options. Left out, as
// no keyword can carry them or none should:
//   * the per-corner overrides of PlanarDomain (vectors indexed by vertex);
//   * the reference field's pointers and the annealer's field callback, which
//     ATLAS::run fills in itself;
//   * templates.wDir/singularity, rewrite.fieldWeight and anneal.wDir/
//     singularity/directedFraction -- ATLAS::search overwrites each from
//     field.* whenever a reference field is on, and none of them does anything
//     without one, so field_w_dir, field_singularity_* and
//     field_directed_fraction are the keywords that reach them;
//   * coarseDomain.maxSpacing, which is set from each of coarseSpacings in
//     turn.
const Table<ATLAS::Options> &optionsTable() {
    using O = ATLAS::Options;
    static const Table<O> t{
        OPT(domain.cornerAngle), OPT(domain.kinkAngle), OPT(domain.materialInterfaces),
        OPT(domain.degenerateArea), OPT(carrier.angleTolerance), OPT(carrier.areaTolerance),
        OPT(runTemplates), OPT(templates.sections), OPT(templates.stars), OPT(templates.ogrids),
        OPT(templates.halfOGrids), OPT(templates.annuli), OPT(templates.halfOGridCorner),
        OPT(templates.reflexAngle), OPT(templates.cornerTolerance), OPT(templates.nodeSnap),
        OPT(templates.coreFraction), OPT(templates.minScaledJacobian),
        OPT(templates.smoothingIterations), OPT(templates.splitBoundary),
        OPT(templates.conePlacement), OPT(templates.chooseByField), OPT(templates.flatNodes),
        OPT(templates.dropBoundary), OPT(templates.generalSections), OPT(templates.flatSections),
        OPT(templates.sectionShift), OPT(templates.propagateNodes), OPT(templates.smallestFirst),
        OPT(rectangles.maxCandidates), OPT(rectangles.maxGrowthSteps), OPT(rectangles.merges),
        OPT(rectangles.growth), OPT(rectangles.maxCutRounds), OPT(rectangles.workFactor),
        OPT(cover.alpha), OPT(cover.beta), OPT(cover.maxSweeps), OPT(cover.singleContact),
        OPT(runRewrite), OPT(rewrite.maxRings), OPT(rewrite.maxCavityCells),
        OPT(rewrite.minScaledJacobian), OPT(rewrite.smoothingIterations), OPT(rewrite.grids),
        OPT(rewrite.stars), OPT(rewrite.maxCornerAngle), OPT(rewrite.realisationsPerWitness),
        OPT(rewrite.cornerPool), OPT(rewrite.maxInsertFactor), OPT(rewrite.lambdaS),
        OPT(rewrite.lambdaQ), OPT(rewrite.boundaryPlusWeight), OPT(rewrite.boundaryMinusWeight),
        OPT(rewrite.coneStars), OPT(rewrite.timeBudget), OPT(rewrite.straddle), OPT(rewrite.straddleReach),
        OPT(rewrite.straddleTemplates), OPT(rewriteRounds),
        OPT(rewritePassesPerRound), OPT(maxRings), OPT(timeBudget), OPT(fineRewrite),
        OPT(coarse), OPT(coarseSpacings), OPT(coarseSeeds), OPT(coarseWithoutTemplates),
        OPT(freeTemplatesMultimat), OPT(multimatOrders), OPT(lockInclusionOGrids), OPT(alignToField), OPT(parallel),
        OPT(coarseDomain.gapFraction), OPT(coarseDomain.curvatureFraction),
        OPT(coarseDomain.grading), OPT(coarseDomain.minAngle), OPT(coarseDomain.minClosedSamples),
        OPT(coarseDomain.minFlipAngle), OPT(coarseDomain.harmoniseCounts),
        OPT(coarseDomain.maxHarmoniseFactor), OPT(coarseDomain.maxHarmoniseAspect),
        OPT(coarseRewriteRounds), OPT(coarsePassesPerRound),
        OPT(coarseMaxRings), OPT(anneal.maxMoves), OPT(anneal.seconds), OPT(anneal.T0),
        OPT(anneal.T1), OPT(anneal.wBlocks), OPT(anneal.wDefect), OPT(anneal.wBoundaryDefect),
        OPT(anneal.wCells), OPT(anneal.wShape), OPT(anneal.shapeFloor),
        OPT(anneal.minScaledJacobian), OPT(anneal.wAlign), OPT(anneal.wBoundaryDefectPlus),
        OPT(anneal.wBoundaryDefectMinus), OPT(anneal.patchTopology), OPT(anneal.maxRings),
        OPT(anneal.realisations), OPT(anneal.straddleFraction), OPT(anneal.freeTemplates), OPT(anneal.seed),
        OPT(realisation.size), OPT(realisation.snap),
        OPT(realisation.smoothingSweeps), OPT(maxRealisations), OPT(maxRealisationsLastResort),
        OPT(field.reference), OPT(field.solve.gamma), OPT(field.solve.maxSteps),
        OPT(field.solve.weight), OPT(field.solve.tauContinuation), OPT(field.solve.tauRatio),
        OPT(field.solve.tauFloorEdges), OPT(field.solve.tauLevelSteps),
        OPT(field.solve.pinDiskCenters), OPT(field.solve.alignToInterfaces),
        OPT(field.solve.cancelDipoles), OPT(field.solve.dipoleRadius), OPT(field.solve.gridStep),
        OPT(field.solve.maxGridNodes), OPT(field.singularity.wPos), OPT(field.singularity.wExtra),
        OPT(field.singularity.wMissing), OPT(field.singularity.wEdgePlus),
        OPT(field.singularity.wEdgeMinus), OPT(field.singularity.r0), OPT(field.wDir),
        OPT(field.templates), OPT(field.greedy), OPT(field.greedyOnFine), OPT(field.anneal),
        OPT(field.directedMoves), OPT(field.directedFraction), OPT(field.arbiter),
        OPT(arbiterShape), OPT(arbiterShapeFloor), OPT(signedDefectsOnFine), OPT(multimatQualityGate), OPT(multimatArbiterWDir),
    };
    return t;
}

const Table<ATLAS::Status> &statusTable() {
    using O = ATLAS::Status;
    static const Table<O> t{
        REP(domainValid),   REP(carrierValid),  REP(coverValid),     REP(triangles),
        REP(initialCells),  REP(chosen),        REP(finalCells),     REP(bestBlocks),
        REP(bestObjective), REP(bestScore),     REP(secondsStage1),  REP(secondsStage2),
        REP(secondsSearch), REP(secondsField),  REP(messages),
    };
    return t;
}

const Table<BlockCover::Report> &coverTable() {
    using O = BlockCover::Report;
    static const Table<O> t{
        REP(candidates),      REP(blocks),          REP(singletonBlocks), REP(largestBlock),
        REP(mapPieces),       REP(macroVertices),   REP(macroEdges),
        REP(irregularMacroVertices), REP(objective), REP(meanDistortion), REP(maxDistortion),
        REP(meanComplexity),  REP(uncovered),       REP(overcovered),     REP(sideMismatches),
        REP(hangingVertices), REP(multipleContacts), REP(protectedNotMacro),
        REP(interfaceInside), REP(uncertified),     REP(eulerLHS),        REP(eulerRHS),
        REP(eulerHolds),      REP(conforming),      REP(valid),           REP(incumbent),
        REP(swaps),           REP(swapTrials),      REP(sweeps),          REP(messages),
    };
    return t;
}

// BlockMesh::Options, less targetEdgeLength (mesh()'s `h`).
const Table<BlockMesh::Options> &meshOptionsTable() {
    using O = BlockMesh::Options;
    static const Table<O> t{
        OPT(minIntervals),    OPT(maxIntervals),       OPT(useChart),
        OPT(smoothingPasses), OPT(smoothingTolerance), OPT(smoothingThreshold),
        OPT(crackTolerance),
    };
    return t;
}

const Table<BlockMesh::Report> &meshReportTable() {
    using O = BlockMesh::Report;
    static const Table<O> t{
        REP(modelExtent),     REP(target),            REP(chords),
        REP(arcsAssigned),    REP(minIntervals),      REP(maxIntervals),
        REP(clampedChords),   REP(meanIntervals),     REP(blocks),
        REP(unmeshedBlocks),  REP(vertices),          REP(quads),
        REP(minEdge),         REP(maxEdge),           REP(meanEdge),
        REP(edgeRatioRms),    REP(worstEdgeRatio),    REP(minScaledJacobian),
        REP(meanScaledJacobian), REP(minScaledJacobianBefore), REP(invertedBefore),
        REP(smoothedBlocks),  REP(smoothingSweeps),   REP(reflexCorners),
        REP(invertedQuads),   REP(minQuadArea),       REP(maxQuadArea),
        REP(meshArea),        REP(domainArea),        REP(boundaryDeviation),
        REP(interfaceDeviation), REP(boundaryNodes),  REP(interiorEdges),
        REP(boundaryEdges),   REP(nonManifoldEdges),  REP(cracks),
        REP(materials),       REP(interfaceEdges),    REP(conforming),
        REP(valid),           REP(messages),
    };
    return t;
}

class AtlasMethod final : public Method {
public:
    AtlasMethod(const Mesh &model, const ATLAS::Options &opts)
        : mesh_(std::make_shared<Mesh>(model)), atlas_(mesh_, opts) {}

    void run() {
        atlas_.run();
        if (atlas_.hasCover()) {
            decomposition_ = atlas_.getCover().blockDecomposition();
            coverage_ = decompositionCoverage(decomposition_, *mesh_);
        }
    }

    const char *name() const override { return "atlas"; }
    const BlockDecomposition &decomposition() const override { return decomposition_; }
    double coverage() const override { return coverage_; }

    PyObject *report() const override {
        if (!atlas_.hasCover()) return toDict({source(statusTable(), atlas_.getStatus())});
        return toDict({source(statusTable(), atlas_.getStatus()),
                       source(coverTable(), atlas_.getCover().getReport(), "cover_")});
    }
    PyObject *options() const override { return toDict({source(optionsTable(), atlas_.getOptions())}); }

    int mesh(double h, PyObject *kwargs, MeshOutput &out) const override {
        BlockMesh::Options mo;
        mo.targetEdgeLength = h;
        mesh::QuadMesh::Options no;
        if (applyKwargs(kwargs, "mesh", {target(meshOptionsTable(), mo),
                                         target(nodeOptionsTable(), no)}) < 0)
            return -1;
        if (!atlas_.hasCover()) {
            PyErr_SetString(PyExc_RuntimeError,
                            "atlas: no cover to mesh -- no search succeeded (see report['messages'])");
            return -1;
        }
        std::unique_ptr<BlockMesh> grid;
        const bool ok = compute([&] {
            grid = std::make_unique<BlockMesh>(atlas_.getCover(), mo);
            if (grid->getReport().quads <= 0) throw std::runtime_error("atlas: the mesh has no elements");
            out.mesh = std::make_unique<mesh::QuadMesh>(mesh::QuadMesh::from(*grid, no));
        });
        if (!ok) return -1;
        out.report = toDict({source(meshReportTable(), grid->getReport())});
        out.smoothing = blockGridSmoothing();
        return out.report ? 0 : -1;
    }

    double defaultTarget() const override { return BlockMesh::Options().targetEdgeLength; }

private:
    std::shared_ptr<Mesh> mesh_;
    ATLAS atlas_;
    BlockDecomposition decomposition_;
    double coverage_ = 0.0;
};

std::shared_ptr<Method> run(const Mesh &model, PyObject *kwargs) {
    ATLAS::Options opts;
    if (applyKwargs(kwargs, "atlas", {target(optionsTable(), opts)}) < 0) return nullptr;
    std::shared_ptr<AtlasMethod> m;
    if (!compute([&] {
            m = std::make_shared<AtlasMethod>(model, opts);
            m->run();
        }))
        return nullptr;
    return m;
}

PyObject *defaults(const std::string &stage) {
    if (stage == "method") return toDict({source(optionsTable(), ATLAS::Options())});
    if (stage == "mesh") {
        const BlockMesh::Options mo;
        const mesh::QuadMesh::Options nodes;
        return meshingDefaults(mo.targetEdgeLength, {source(meshOptionsTable(), mo),
                                                     source(nodeOptionsTable(), nodes)});
    }
    if (stage == "smooth") return smoothingDefaults(blockGridSmoothing());
    return unknownStage("atlas", stage);
}

void buildTables() {
    optionsTable();
    statusTable();
    coverTable();
    meshOptionsTable();
    meshReportTable();
}

}  // namespace

const MethodSpec atlasSpec = {
    "atlas",
    "atlas(**options) -> BlockDecomposition\n\n"
    "ATLAS, the square-transport blocking (src/ATLAS): a certified square\n"
    "carrier, whole-region templates, certified rectangles, an exact cover,\n"
    "and cavity rewrites, searched coarse and fine with a DualMBO reference\n"
    "field. The cover tiles the model, so coverage is 1 when it succeeds.\n\n"
    "Keywords are the fields of ATLAS::Options in snake_case, nested stages\n"
    "prefixed (domain_corner_angle, anneal_seed, field_reference='dual_mbo',\n"
    "coarse_spacings=[0.125, 0.07], ...); crossgen.options('atlas') lists\n"
    "them. Meshing is BlockMesh on the certified charts.",
    run,
    defaults,
    buildTables,
};

}  // namespace pycg
