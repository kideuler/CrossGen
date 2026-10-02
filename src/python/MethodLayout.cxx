// Mesh.meridian() and Mesh.torsion(): Pipelines A and B (src/MERIDIAN,
// src/TORSION) run through Stage 9, and meshed by the Stage 10 they share, with
// Stage 11's O-grid templates where Stage 0c excised a disk.
//
// Both pipelines meet at Immersion and share every stage after it, and their
// classes expose the same accessors, so one binding serves both through Traits.
//
// ### Why Stage 10 is run here and not in run()
//
// A pipeline's run() meshes once, at Options::quadTargetEdge, and has no way to
// mesh again. BlockDecomposition.mesh(h) has to be callable at any h, so the
// pipeline is run with runQuadMesh off and Stages 10 and 11 are built here from
// the pipeline's public stages instead -- the viewer's precedent
// (CrossGenWidget's Stage 10 and runDiskTemplates), and the same calls in the
// same order as the end of MERIDIAN::run() and TORSION::run(): the rim arcs of
// the excised disks become parity constraints on the interval assignment, the
// QuadMesh is built on the spline fit, and where disks were excised the
// DiskTemplate fills them. Stage 10's options default to the pipeline's quad*
// options, which is what run() would have passed. If run() ever changes how it
// builds Stage 10, this has to follow.
//
// The quad* options and the run* switches are therefore not keywords of
// meridian()/torsion(): the former are mesh()'s, and the latter would stop the
// pipeline short of the decomposition it exists to produce.
#include "python/Methods.hxx"
#include "python/PyTables.hxx"
#include "python/PyUtil.hxx"

#include <cmath>
#include <stdexcept>

#include "MERIDIAN/DiskTemplate.hxx"
#include "MERIDIAN/MERIDIAN.hxx"
#include "MERIDIAN/QuadMesh.hxx"
#include "TORSION/TORSION.hxx"

#define OPT(path) CG_OPTION(O, path)
#define REP(path) CG_REPORT(O, path)

namespace pycg {

template <>
struct EnumNames<TORSION::Options::Reference> {
    static constexpr std::pair<TORSION::Options::Reference, const char *> list[] = {
        {TORSION::Options::Reference::Cone, "cone"},
        {TORSION::Options::Reference::Induced, "induced"},
        {TORSION::Options::Reference::Field, "field"},
        {TORSION::Options::Reference::Euclidean, "euclidean"},
        {TORSION::Options::Reference::Ricci, "ricci"},
    };
};

namespace {

// ---- Stage 10 and 11, shared -----------------------------------------------

// ::QuadMesh::Options, less targetEdgeLength (mesh()'s `h`) and evenLoops (the
// rims, which are the pipeline's to supply).
const Table<::QuadMesh::Options> &stage10OptionsTable() {
    using O = ::QuadMesh::Options;
    static const Table<O> t{
        OPT(minIntervals),     OPT(maxIntervals),         OPT(collapseSpan),
        OPT(minLoopEdges),     OPT(arcLengthSamples),     OPT(useSplines),
        OPT(featuresOnTracedArcs), OPT(spanSamples),      OPT(smoothingPasses),
        OPT(smoothingTolerance), OPT(smoothingThreshold), OPT(crackTolerance),
        OPT(materialsFromFaces), OPT(contractOntoFeatures),
    };
    return t;
}

const Table<::QuadMesh::Report> &stage10ReportTable() {
    using O = ::QuadMesh::Report;
    static const Table<O> t{
        REP(modelExtent),     REP(target),          REP(chords),
        REP(arcsAssigned),    REP(minIntervals),    REP(maxIntervals),
        REP(clampedChords),   REP(meanIntervals),   REP(oddLoops),
        REP(parityChordsMoved), REP(oddLoopsLeft),  REP(parityCost),
        REP(blocks),          REP(unmeshedPatches), REP(unmeshedArea),
        REP(collapsedChords), REP(collapsedPatches), REP(collapsedArea),
        REP(weldedVertices),  REP(vertices),        REP(quads),
        REP(minEdge),         REP(maxEdge),         REP(meanEdge),
        REP(edgeRatioRms),    REP(worstEdgeRatio),  REP(minScaledJacobian),
        REP(meanScaledJacobian), REP(minScaledJacobianBefore), REP(invertedBefore),
        REP(smoothedBlocks),  REP(smoothingSweeps), REP(reflexCorners),
        REP(invertedQuads),   REP(minQuadArea),     REP(maxQuadArea),
        REP(meshArea),        REP(patchArea),       REP(interiorEdges),
        REP(boundaryEdges),   REP(nonManifoldEdges), REP(cracks),
        REP(materials),       REP(mixedQuads),      REP(unlocatedQuads),
        REP(relabelledQuads),
        REP(interfaceEdges),  REP(conforming),      REP(materialsPure),
        REP(valid),           REP(messages),
    };
    return t;
}

const Table<DiskTemplate::Report> &stage11ReportTable() {
    using O = DiskTemplate::Report;
    static const Table<O> t{
        REP(inclusions),      REP(filled),          REP(refusedOdd),
        REP(refusedShort),    REP(refusedOpen),     REP(blocks),
        REP(quads),           REP(vertices),        REP(mergedVertices),
        REP(mergedQuads),     REP(minScaledJacobian), REP(meanScaledJacobian),
        REP(templateMinScaledJacobian), REP(invertedQuads), REP(minEdge),
        REP(maxEdge),         REP(smoothingSweeps), REP(interiorEdges),
        REP(boundaryEdges),   REP(nonManifoldEdges), REP(cracks),
        REP(conforming),      REP(valid),           REP(messages),
    };
    return t;
}

// What MERIDIAN::run() and TORSION::run() hand Stage 10, from their options.
template <class PO>
::QuadMesh::Options stage10Options(const PO &o) {
    ::QuadMesh::Options q;
    q.targetEdgeLength = o.quadTargetEdge;
    q.minIntervals = o.quadMinIntervals;
    q.maxIntervals = o.quadMaxIntervals;
    q.collapseSpan = o.quadCollapseSpan;
    q.useSplines = o.quadUseSplines;
    q.featuresOnTracedArcs = o.quadFeaturesOnTracedArcs;
    q.materialsFromFaces = o.quadMaterialsFromFaces;
    q.contractOntoFeatures = o.quadContractOntoFeatures;
    q.smoothingPasses = o.quadSmoothingPasses;
    q.smoothingThreshold = o.quadSmoothingThreshold;
    return q;
}

// ... and Stage 11.
template <class PO>
DiskTemplate::Options stage11Options(const PO &o) {
    DiskTemplate::Options t;
    t.coreSquareness = o.diskCoreSquareness;
    t.ringDepth = o.diskRingDepth;
    t.smoothingPasses = o.diskSmoothingPasses;
    return t;
}

// Stage 12 as TestMERIDIAN and TestTORSION run it (and the viewer's pipeline
// modes): metric 7 sampled at the element corners, exponent 2 -- after the
// flat feature corners have been pillowed, which MeshOutput::pillow asks of
// the first smooth(). Corners since 2026-10-02: at the 2x2 Gauss points TMOP
// turned corners over on these meshes that it could not see (multimat/tooth
// 0.00 -> -0.96).
mesh::TMOP::Options stage12Smoothing() {
    mesh::TMOP::Options t;
    t.metric = mesh::TMOP::ShapeSize007;
    t.quadrature = mesh::TMOP::Corners;
    t.maxSweeps = 1000;
    return t;
}

// ---- the two pipelines ---------------------------------------------------------
template <class P>
struct Traits;

template <>
struct Traits<MERIDIAN> {
    static constexpr const char *name = "meridian";
    static constexpr const char *source = "MERIDIAN";

    static const Table<MERIDIAN::Options> &options() {
        using O = MERIDIAN::Options;
        static const Table<O> t{
            OPT(dualMBOGamma), OPT(dualMBOMaxSteps), OPT(materialInterfaces),
            OPT(interfaceKinkAngle), OPT(interfaceLoopSplits), OPT(splitCircleLoops),
            OPT(diskTemplates), OPT(diskCoreSquareness), OPT(diskRingDepth),
            OPT(diskSmoothingPasses), OPT(alignFieldToInterfaces), OPT(cancelInterfaceDipoles),
            OPT(prescribeInterfaceCones), OPT(interfaceCorners), OPT(propagateInterfaceLabels),
            OPT(seamTurnInterfaceLabels), OPT(lagInterfaceScales), OPT(autoRebalance),
            OPT(minBoundaryIndex), OPT(maxBoundaryIndex), OPT(relocateFlatCones),
            OPT(flatConeAngle), OPT(keepFlatConesOnFailure), OPT(coneCutsToBoundary),
            OPT(coneCutInterfaceAvoidance), OPT(ricciTolerance), OPT(ricciMaxIterations),
            OPT(delaunayFlips), OPT(lambdaInit), OPT(lambdaGrowth), OPT(outerSteps),
            OPT(innerIterations), OPT(alternateReference), OPT(relabelBetweenSteps),
            OPT(seedTopoConstraints), OPT(topoNearMiss), OPT(topoNearMissRetry),
            OPT(seedSelfReturns), OPT(seedAllConnections), OPT(separatrixSnap),
            OPT(separatrixMaxSteps), OPT(separatrixSnapRings), OPT(separatrixDetectCycles),
            OPT(arrangementMerge), OPT(arrangementCorner), OPT(arrangementCollapse),
            OPT(arrangementTrim), OPT(splineSegments), OPT(splineSamples), OPT(fitBoundaryArcs),
            OPT(fitInterfaceArcs), OPT(repairPasses), OPT(repairGapLimit), OPT(repairLambdaBoost),
            OPT(repairOuterSteps), OPT(repairMaxPerPass), OPT(repairScoresArrangement),
            OPT(repairPatience),
        };
        return t;
    }

    // MERIDIAN::Status less the Stage 10 and 11 fields, which run() never
    // reaches here; mesh()'s report carries their numbers instead.
    static const Table<MERIDIAN::Status> &status() {
        using O = MERIDIAN::Status;
        static const Table<O> t{
            REP(mboSteps), REP(fieldConverged), REP(fieldAlignedToInterfaces), REP(coneDipoleUnits),
            REP(diskInclusions), REP(diskTrianglesExcised), REP(materials), REP(interfaceEdges),
            REP(interfaceBranches), REP(interfaceNodes), REP(interfaceIllPosedNodes),
            REP(interfaceWorstSector), REP(interfaceConesPrescribed),
            REP(interfacePrescriptionShift), REP(regions), REP(regionsBalanced),
            REP(regionQuartersMoved), REP(regionCornersInserted), REP(worstRegionDeficit),
            REP(conesAdmissible), REP(interiorCones), REP(boundaryCones), REP(cutIsDisk),
            REP(allConesOnBoundary), REP(cutInterfaceEdges), REP(cutInterfaceVertices),
            REP(cutInterfaceNodes), REP(ricciConverged), REP(readyForImmersion),
            REP(immersionValid), REP(immersionFlippedFaces), REP(seamArcs), REP(boundaryEdgesU),
            REP(boundaryEdgesV), REP(featureChains), REP(topoPaths), REP(topoSelfReturns),
            REP(topoExtraPerPair), REP(topoNearMissUsed), REP(topoNearMissRetried),
            REP(interfaceCorners), REP(interfaceLabelsCorrected), REP(interfaceChainsSeamFlipped),
            REP(interfaceResidual), REP(interfaceCornerChanges), REP(interfacesAligned),
            REP(layoutRan), REP(layoutInjective), REP(layoutConstrained), REP(outerStepsTaken),
            REP(layoutValid), REP(separatricesRan), REP(separatrices), REP(separatricesToCone),
            REP(separatricesToBoundary), REP(separatricesUnresolved), REP(separatricesNearMisses),
            REP(arrangementRan), REP(layoutNodes), REP(layoutArcs), REP(layoutPatches),
            REP(layoutQuads), REP(layoutCoverage), REP(arrangementValid),
            REP(splinesRan), REP(splinePatches), REP(splineControlPoints), REP(splineMaxDeviation),
            REP(splinesWatertight), REP(splinesValid), REP(repairPasses),
            REP(repairConstraintsAdded), REP(q5Verified), REP(messages),
        };
        return t;
    }
};

template <>
struct Traits<TORSION> {
    static constexpr const char *name = "torsion";
    static constexpr const char *source = "TORSION";

    // Less externalField and sizing too: a complex vector per face and a size
    // per vertex, neither of which a keyword can sensibly carry -- and less
    // extraEmitters, a list of vertices the per-material mode sets on each
    // region's own run and nothing else has a use for.
    static const Table<TORSION::Options> &options() {
        using O = TORSION::Options;
        static const Table<O> t{
            OPT(dualMBOGamma), OPT(dualMBOMaxSteps), OPT(dualMBOWeight),
            OPT(dualMBOTauContinuation), OPT(dualMBOTauRatio), OPT(dualMBOTauFloorEdges),
            OPT(dualMBOTauLevelSteps), OPT(perMaterial), OPT(perMaterialRounds),
            OPT(perMaterialThreads), OPT(perMaterialMatchTolerance), OPT(perMaterialKeepDipoles),
            OPT(perMaterialFallback), OPT(perMaterialSpacedMatching), OPT(perMaterialBendEnds),
            OPT(cancelSingleMaterialDipoles),
            OPT(materialInterfaces), OPT(interfaceKinkAngle),
            OPT(interfaceLoopSplits), OPT(diskTemplates), OPT(diskCoreSquareness),
            OPT(diskRingDepth), OPT(diskSmoothingPasses), OPT(alignFieldToInterfaces),
            OPT(cancelInterfaceDipoles), OPT(prescribeInterfaceCones), OPT(interfaceCorners),
            OPT(propagateInterfaceLabels), OPT(seamTurnInterfaceLabels), OPT(lagInterfaceScales),
            OPT(autoRebalance), OPT(minBoundaryIndex), OPT(maxBoundaryIndex),
            OPT(relocateFlatCones), OPT(flatConeAngle), OPT(keepFlatConesOnFailure),
            OPT(coneCutsToBoundary), OPT(coneCutInterfaceAvoidance), OPT(targetEdge),
            OPT(checkSecondSeed), OPT(alignInIntegration), OPT(alignmentFallback),
            OPT(alignAcrossSeams), OPT(alignmentChooseByLadder), OPT(reconcileSectors),
            OPT(alignmentReleaseRounds), OPT(localUntangleFaceFraction),
            OPT(localUntangleFaceFloor),
            OPT(pullOntoAlignment), OPT(pullOuterSteps), OPT(retryKeepingDipoles),
            OPT(integrationRegularisation), OPT(untangle), OPT(localUntangle),
            OPT(localUntangleSweeps), OPT(regularisedUntangleRings),
            OPT(regularisedUntangleIterations), OPT(pairSeamInUntangle), OPT(untangleOuterSteps),
            OPT(untangleInnerIterations), OPT(untangleSeamWeight), OPT(untangleAlignWeight),
            OPT(reprojectAfterUntangle), OPT(reference), OPT(conformalSizing), OPT(lambdaInit),
            OPT(lambdaGrowth), OPT(lambdaAlignmentFactor), OPT(lambdaSeamFactor), OPT(outerSteps),
            OPT(innerIterations), OPT(relabelBetweenSteps), OPT(seedTopoConstraints),
            OPT(topoNearMiss), OPT(topoNearMissRetry), OPT(topoNearMissLastRetry),
            OPT(topoRetryUnseeded), OPT(seedSelfReturns),
            OPT(seedAllConnections), OPT(separatrixSnap), OPT(separatrixMaxSteps),
            OPT(separatrixSnapRings), OPT(separatrixDetectCycles), OPT(arrangementMerge),
            OPT(arrangementCorner), OPT(arrangementCollapse), OPT(arrangementTrim),
            OPT(splineSegments), OPT(splineSamples), OPT(fitBoundaryArcs), OPT(fitInterfaceArcs),
            OPT(repairPasses), OPT(repairGapLimit), OPT(repairLambdaBoost), OPT(repairOuterSteps),
            OPT(repairMaxPerPass), OPT(repairScoresArrangement), OPT(repairPatience),
        };
        return t;
    }

    // TORSION::Status less the Stage 10 and 11 fields, as for MERIDIAN.
    static const Table<TORSION::Status> &status() {
        using O = TORSION::Status;
        static const Table<O> t{
            REP(mboSteps), REP(fieldConverged), REP(externalFieldUsed),
            REP(fieldAlignedToInterfaces), REP(coneDipoleUnits), REP(dipoleRetryRan),
            REP(dipolesKept), REP(materials),
            REP(interfaceEdges), REP(interfaceBranches), REP(interfaceNodes),
            REP(interfaceIllPosedNodes), REP(interfaceWorstSector), REP(interfaceConesPrescribed),
            REP(interfacePrescriptionShift), REP(regions), REP(regionsBalanced),
            REP(regionQuartersMoved), REP(regionCornersInserted), REP(worstRegionDeficit),
            REP(conesAdmissible), REP(interiorCones), REP(boundaryCones), REP(cutIsDisk),
            REP(allConesOnBoundary), REP(cutInteriorJunctions), REP(cutInterfaceEdges),
            REP(cutInterfaceVertices), REP(cutInterfaceNodes), REP(framesRan), REP(combedFaces),
            REP(unreachedFaces), REP(combingDefects), REP(maxFrameJump), REP(indexMismatches),
            REP(clusteredConePairs), REP(highIndexCones), REP(leftHandedFrames),
            REP(maxMetricDisagreement), REP(secondSeedAgrees), REP(framesValid),
            REP(alignmentChains), REP(alignedBoundaryEdges), REP(alignedInterfaceEdges),
            REP(alignmentOverrides), REP(alignmentSeamCrossings), REP(alignmentClosedChains),
            REP(maxAlignmentResidual), REP(integrationRan), REP(integrationSolved),
            REP(integrationAlignedEdges), REP(integrationAlignResidual),
            REP(integrationAlignStrain), REP(alignmentWasDropped), REP(alignmentChainsReleased),
            REP(alignmentHeld),
            REP(integrationSeamResidual), REP(integrationFlippedFaces),
            REP(integrationFlippedAreaFraction), REP(integrationMinAreaRatio),
            REP(integrationMaxFitResidual), REP(integrationMeanFitResidual),
            REP(integrationFlipsAtCones), REP(integrationNearestFlipToCone), REP(localUntangleRan),
            REP(localUntangleFlippedFaces), REP(localUntangleLevel), REP(untangleRan),
            REP(tutteValid), REP(untangleOuterSteps), REP(untangleFlippedFaces),
            REP(untangleSeamResidual), REP(reprojectionRan), REP(reprojectionKept),
            REP(reprojectionFlippedFaces), REP(pullRan), REP(pullKept), REP(pullProjected),
            REP(pullBoundaryResidual), REP(pullFeatureResidual), REP(immersionValid),
            REP(immersionFlippedFaces), REP(seamArcs), REP(maxSnapError), REP(frameKConflicts),
            REP(maxMetricResidual),
            REP(boundaryEdgesU), REP(boundaryEdgesV), REP(featureChains), REP(topoPaths),
            REP(topoSelfReturns), REP(topoExtraPerPair), REP(topoNearMissUsed),
            REP(topoNearMissRetried), REP(interfaceCorners), REP(interfaceLabelsCorrected),
            REP(interfaceChainsSeamFlipped), REP(interfaceResidual), REP(interfaceCornerChanges),
            REP(interfacesAligned), REP(layoutRan), REP(layoutInjective), REP(layoutConstrained),
            REP(outerStepsTaken), REP(layoutValid), REP(layoutConeAngleResidual),
            REP(layoutRegularAngleResidual), REP(layoutReferenceLeftHanded),
            REP(referenceConeResidual), REP(coneMetricRan), REP(coneMetricSolved),
            REP(coneMetricConverged), REP(coneMetricNewtonSteps), REP(coneMetricInitialError),
            REP(coneMetricLinearError), REP(coneMetricFinalError), REP(coneMetricConeResidual),
            REP(coneMetricMinScale), REP(coneMetricMaxScale), REP(coneMetricNonRealisable),
            REP(conformalSizingUsed), REP(separatricesRan), REP(separatrices),
            REP(separatricesToCone), REP(separatricesToBoundary), REP(separatricesUnresolved),
            REP(separatricesNearMisses), REP(arrangementRan), REP(layoutNodes), REP(layoutArcs),
            REP(layoutPatches), REP(layoutQuads), REP(layoutCoverage), REP(arrangementValid),
            REP(splinesRan), REP(splinePatches), REP(splineControlPoints), REP(splineMaxDeviation),
            REP(splinesWatertight), REP(splinesValid), REP(diskInclusions),
            REP(diskTrianglesExcised), REP(repairPasses), REP(repairConstraintsAdded),
            REP(q5Verified), REP(perMaterialRan), REP(perMaterialKept),
            REP(perMaterialFallbackRan), REP(perMaterialRegions), REP(perMaterialRegionsValid),
            REP(perMaterialRounds), REP(perMaterialConverged), REP(perMaterialEmitters),
            REP(perMaterialMatched), REP(perMaterialUnmatched), REP(perMaterialDipolesKept),
            REP(perMaterialSeconds), REP(messages),
        };
        return t;
    }
};

template <class P>
class LayoutMethod final : public Method {
public:
    using Options = typename P::Options;

    LayoutMethod(const Mesh &model, const Options &opts)
        : opts_(opts), pipeline_(std::make_shared<Mesh>(model), opts) {}

    void run() {
        pipeline_.run();
        if (pipeline_.hasArrangement()) {
            decomposition_ = pipeline_.getArrangement().blockDecomposition(
                pipeline_.hasSplines() ? &pipeline_.getSplines() : nullptr, 16, Traits<P>::source);
        }
    }

    const char *name() const override { return Traits<P>::name; }
    const BlockDecomposition &decomposition() const override { return decomposition_; }

    // The area of the arrangement's faces that are blocks -- four corners, one
    // arc a side, the ones Stage 10 meshes -- over the area of S. S is the
    // model after Stage 0c: an excised disk is not part of it, and comes back
    // as an O-grid at mesh time instead.
    double coverage() const override {
        if (!pipeline_.hasArrangement()) return 0.0;
        const Arrangement &arr = pipeline_.getArrangement();
        const double domain = arr.getReport().domainArea;
        if (!(domain > 0.0)) return 0.0;
        const std::vector<Arrangement::Face> &faces = arr.getFaces();
        double blocks = 0.0;
        for (const int f : arr.patchFaces()) {
            if (f < 0 || f >= static_cast<int>(faces.size())) continue;
            if (faces[f].quad && faces[f].simple) blocks += std::fabs(faces[f].area);
        }
        return blocks / domain;
    }

    PyObject *report() const override {
        return toDict({source(Traits<P>::status(), pipeline_.getStatus())});
    }
    PyObject *options() const override { return toDict({source(Traits<P>::options(), opts_)}); }

    int mesh(double h, PyObject *kwargs, MeshOutput &out) const override {
        ::QuadMesh::Options qo = stage10Options(opts_);
        qo.targetEdgeLength = h;
        mesh::QuadMesh::Options no;
        if (applyKwargs(kwargs, "mesh", {target(stage10OptionsTable(), qo),
                                         target(nodeOptionsTable(), no)}) < 0)
            return -1;
        if (!pipeline_.hasSplines()) {
            PyErr_Format(PyExc_RuntimeError,
                         "%s: no spline fit to mesh -- the pipeline stopped before Stage 9 "
                         "(see report['messages'])",
                         Traits<P>::name);
            return -1;
        }
        std::unique_ptr<::QuadMesh> quads;
        std::unique_ptr<DiskTemplate> fill;
        const bool ok = compute([&] {
            const Arrangement &arr = pipeline_.getArrangement();
            const std::vector<DiskTemplate::Inclusion> &inclusions = pipeline_.getInclusions();
            // Stage 10. Each excised rim has to come out with an even number
            // of edges or Stage 11 has nothing it can fill it with.
            const std::vector<std::vector<int>> rims = DiskTemplate::rimArcs(arr, inclusions);
            for (const std::vector<int> &r : rims)
                if (!r.empty()) qo.evenLoops.push_back(r);
            quads = std::make_unique<::QuadMesh>(pipeline_.getSplines(), qo);
            if (quads->getReport().quads <= 0)
                throw std::runtime_error(std::string(Traits<P>::name) + ": the mesh has no elements");
            // Stage 11.
            if (!inclusions.empty()) {
                fill = std::make_unique<DiskTemplate>(
                    quads->vertices(), quads->quads(), quads->quadMaterials(),
                    DiskTemplate::rimVertexLoops(arr, *quads, rims), inclusions,
                    stage11Options(opts_));
                out.mesh = std::make_unique<mesh::QuadMesh>(mesh::QuadMesh::from(*fill, no));
            } else {
                out.mesh = std::make_unique<mesh::QuadMesh>(mesh::QuadMesh::from(*quads, no));
            }
        });
        if (!ok) return -1;
        out.report = fill ? toDict({source(stage10ReportTable(), quads->getReport()),
                                    source(stage11ReportTable(), fill->getReport(),
                                           "disk_templates_")})
                          : toDict({source(stage10ReportTable(), quads->getReport())});
        out.smoothing = stage12Smoothing();
        out.pillow = true;
        return out.report ? 0 : -1;
    }

    double defaultTarget() const override { return opts_.quadTargetEdge; }

private:
    Options opts_;
    P pipeline_;
    BlockDecomposition decomposition_;
};

template <class P>
std::shared_ptr<Method> runLayout(const Mesh &model, PyObject *kwargs) {
    typename P::Options opts;
    if (applyKwargs(kwargs, Traits<P>::name, {target(Traits<P>::options(), opts)}) < 0)
        return nullptr;
    opts.runQuadMesh = false;
    std::shared_ptr<LayoutMethod<P>> m;
    if (!compute([&] {
            m = std::make_shared<LayoutMethod<P>>(model, opts);
            m->run();
        }))
        return nullptr;
    return m;
}

template <class P>
PyObject *layoutDefaults(const std::string &stage) {
    const typename P::Options o;
    if (stage == "method") return toDict({source(Traits<P>::options(), o)});
    if (stage == "mesh") {
        const ::QuadMesh::Options qo = stage10Options(o);
        const mesh::QuadMesh::Options nodes;
        return meshingDefaults(qo.targetEdgeLength, {source(stage10OptionsTable(), qo),
                                                     source(nodeOptionsTable(), nodes)});
    }
    if (stage == "smooth") return smoothingDefaults(stage12Smoothing(), true);
    return unknownStage(Traits<P>::name, stage);
}

template <class P>
void buildLayoutTables() {
    Traits<P>::options();
    Traits<P>::status();
    stage10OptionsTable();
    stage10ReportTable();
    stage11ReportTable();
}

}  // namespace

const MethodSpec meridianSpec = {
    "meridian",
    "meridian(**options) -> BlockDecomposition\n\n"
    "MERIDIAN, Pipeline A (Shepherd, Gu and Hughes 2022): the DualMBO cross\n"
    "field's cones, a cut to the boundary, discrete Ricci flow, the layout\n"
    "energy, separatrices, their arrangement and the spline fit, Stages 0-9.\n"
    "Meshing is Stage 10 on the fitted patches, with Stage 11's O-grids where\n"
    "disk_templates=True excised an inclusion.\n\n"
    "Keywords are the fields of MERIDIAN::Options in snake_case\n"
    "(src/MERIDIAN/MERIDIAN.hxx), less the quad_* ones (BlockDecomposition.mesh\n"
    "takes those) and the run_* switches; crossgen.options('meridian') lists\n"
    "them. For the bubbles model: disk_templates=True, topo_near_miss=0.04.",
    runLayout<MERIDIAN>,
    layoutDefaults<MERIDIAN>,
    buildLayoutTables<MERIDIAN>,
};

const MethodSpec torsionSpec = {
    "torsion",
    "torsion(**options) -> BlockDecomposition\n\n"
    "TORSION, Pipeline B: MERIDIAN's stages with Stages 3-4 replaced by\n"
    "integrating the DualMBO cross field (ConeMetric, FieldFrames,\n"
    "FieldIntegration, TutteEmbedding); every stage from Immersion on is\n"
    "shared, meshing included. On a multi-material model each material region\n"
    "is laid out on its own and the layouts are matched across the interfaces\n"
    "and glued (per_material=True, the default; per_material=False lays the\n"
    "whole model out at once).\n\n"
    "Keywords are the fields of TORSION::Options in snake_case\n"
    "(src/TORSION/TORSION.hxx), less quad_*, run_*, external_field, extra_emitters\n"
    "and sizing;\n"
    "e.g. dual_mbo_weight='orthogonal', reference='cone'.\n"
    "crossgen.options('torsion') lists them.",
    runLayout<TORSION>,
    layoutDefaults<TORSION>,
    buildLayoutTables<TORSION>,
};

}  // namespace pycg
