#include "python/PyTables.hxx"

// Every field of each struct that has a keyword, in the struct's own order.
// OPT is read-write, REP read-only; `O` is the struct the enclosing function
// builds a table for.
#define OPT(path) CG_OPTION(O, path)
#define REP(path) CG_REPORT(O, path)

namespace pycg {

const Table<mesh::QuadMesh::Options> &nodeOptionsTable() {
    using O = mesh::QuadMesh::Options;
    static const Table<O> t{
        OPT(cornerAngle), OPT(fixAllFeatureNodes), OPT(interfacesAreFeatures),
        OPT(curveSource), OPT(curveMaxDeviation),
    };
    return t;
}

const Table<mesh::QuadMesh::Quality> &qualityTable() {
    using O = mesh::QuadMesh::Quality;
    static const Table<O> t{
        REP(vertices),          REP(quads),              REP(freeNodes),
        REP(slidingNodes),      REP(fixedNodes),         REP(minScaledJacobian),
        REP(meanScaledJacobian), REP(invertedQuads),     REP(nonConvexQuads),
        REP(minArea),           REP(maxArea),            REP(totalArea),
        REP(minEdge),           REP(maxEdge),            REP(meanEdge),
        REP(worstAspect),       REP(nonManifoldEdges),   REP(boundaryEdges),
        REP(boundaryLoops),     REP(allCounterClockwise),
    };
    return t;
}

const Table<mesh::TMOP::Options> &tmopOptionsTable() {
    using O = mesh::TMOP::Options;
    static const Table<O> t{
        OPT(metric),          OPT(gamma),           OPT(exponent),
        OPT(target),          OPT(targetSize),      OPT(quadrature),
        OPT(moveTolerance),   OPT(energyTolerance), OPT(maxLineSearch),
        OPT(maxStepFraction), OPT(hessianFloor),    OPT(untangle),
        OPT(untangleMaxSweeps), OPT(untangler),     OPT(untangleSizeWeight),
        OPT(untangleFloor),   OPT(threads),         OPT(verbose),
    };
    return t;
}

const Table<mesh::TMOP::Report> &tmopReportTable() {
    using O = mesh::TMOP::Report;
    static const Table<O> t{
        REP(ran),                     REP(converged),               REP(sweeps),
        REP(untangleSweeps),          REP(colors),                  REP(threads),
        REP(openMP),                  REP(movableNodes),            REP(freeNodes),
        REP(slidingNodes),            REP(fixedNodes),              REP(frozenCorners),
        REP(featureCurves),           REP(fittedCurves),            REP(curveNodes),
        REP(curveBow),                REP(curveDeviation),          REP(energyBefore),
        REP(energyAfter),             REP(minScaledJacobianBefore), REP(minScaledJacobianAfter),
        REP(meanScaledJacobianBefore), REP(meanScaledJacobianAfter), REP(invertedBefore),
        REP(invertedAfter),           REP(worstAspectBefore),       REP(worstAspectAfter),
        REP(minAreaBefore),           REP(minAreaAfter),            REP(areaBefore),
        REP(areaAfter),               REP(maxDisplacement),         REP(meanDisplacement),
        REP(lastSweepMove),           REP(seconds),                 REP(messages),
    };
    return t;
}

const Table<mesh::Pillow::Report> &pillowReportTable() {
    using O = mesh::Pillow::Report;
    static const Table<O> t{
        REP(ran),                     REP(defectsBefore),           REP(defectsAfter),
        REP(flattestBefore),          REP(flattestAfter),           REP(layers),
        REP(closedLayers),            REP(skipped),                 REP(layerQuads),
        REP(splitQuads),              REP(featureSplits),           REP(quadsBefore),
        REP(quadsAfter),              REP(verticesBefore),          REP(verticesAfter),
        REP(minScaledJacobianBefore), REP(minScaledJacobianAfter),  REP(invertedBefore),
        REP(invertedAfter),           REP(messages),
    };
    return t;
}

const Table<BlockQuadMesh::Options> &blockQuadMeshOptionsTable() {
    using O = BlockQuadMesh::Options;
    static const Table<O> t{
        OPT(minIntervals),       OPT(maxIntervals),       OPT(smoothingPasses),
        OPT(smoothingTolerance), OPT(smoothingThreshold), OPT(crackTolerance),
    };
    return t;
}

const Table<BlockQuadMesh::Report> &blockQuadMeshReportTable() {
    using O = BlockQuadMesh::Report;
    static const Table<O> t{
        REP(modelExtent),        REP(target),             REP(chords),
        REP(edgesAssigned),      REP(minIntervals),       REP(maxIntervals),
        REP(clampedChords),      REP(meanIntervals),      REP(blocks),
        REP(unmeshedBlocks),     REP(vertices),           REP(quads),
        REP(minEdge),            REP(maxEdge),            REP(meanEdge),
        REP(edgeRatioRms),       REP(worstEdgeRatio),     REP(minScaledJacobian),
        REP(meanScaledJacobian), REP(minScaledJacobianBefore), REP(invertedBefore),
        REP(smoothedBlocks),     REP(smoothingSweeps),    REP(reflexCorners),
        REP(invertedQuads),      REP(minQuadArea),        REP(maxQuadArea),
        REP(meshArea),           REP(boundaryDeviation),  REP(interfaceDeviation),
        REP(boundaryNodes),      REP(interiorEdges),      REP(boundaryEdges),
        REP(nonManifoldEdges),   REP(cracks),             REP(materials),
        REP(interfaceEdges),     REP(conforming),         REP(valid),
        REP(messages),
    };
    return t;
}

const Table<shapedna::ShapeDNA::Options> &shapeDNAOptionsTable() {
    using O = shapedna::ShapeDNA::Options;
    static const Table<O> t{
        OPT(count),     OPT(degree),      OPT(refine),     OPT(boundary),
        OPT(normalization), OPT(blockSize), OPT(basisSize), OPT(tolerance),
        OPT(maxRestarts), OPT(shift),     OPT(seed),       OPT(denseBelow),
        OPT(leafSize),  OPT(eigenfunctions), OPT(threads),
    };
    return t;
}

const Table<shapedna::ShapeDNA::Report> &shapeDNAReportTable() {
    using O = shapedna::ShapeDNA::Report;
    static const Table<O> t{
        REP(ran),             REP(converged),       REP(dense),
        REP(triangles),       REP(nodes),           REP(unknowns),
        REP(nnz),             REP(degenerateElements), REP(area),
        REP(boundaryLength),  REP(components),      REP(eulerCharacteristic),
        REP(cornerTerm),      REP(shift),           REP(zeroModes),
        REP(eigenvalues),     REP(supernodes),      REP(maxFront),
        REP(treeDepth),       REP(nnzL),            REP(factorFlops),
        REP(basisSize),       REP(blockSize),       REP(restarts),
        REP(blockSteps),      REP(solves),          REP(deflations),
        REP(maxRitzResidual), REP(maxTrueResidual), REP(normalizer),
        REP(threads),         REP(openMP),          REP(secondsAssemble),
        REP(secondsOrder),    REP(secondsFactor),   REP(secondsEigen),
        REP(seconds),         REP(messages),
    };
    return t;
}

const Table<BoundaryFeatures::Options> &boundaryFeaturesOptionsTable() {
    using O = BoundaryFeatures::Options;
    static const Table<O> t{OPT(cornerAngle), OPT(curvedAngle)};
    return t;
}

const Table<BoundaryFeatures::Summary> &boundaryFeaturesSummaryTable() {
    using O = BoundaryFeatures::Summary;
    static const Table<O> t{
        REP(regions),            REP(holes),              REP(euler),
        REP(corners),            REP(cornersOneBlock),    REP(cornersTwoBlocks),
        REP(cornersThreeBlocks), REP(cornersFourBlocks),  REP(acuteCorners),
        REP(ambiguousCorners),   REP(cornerDefect),       REP(singularityBound),
        REP(minimumDefect),      REP(isoperimetricRatio), REP(curvedFraction),
        REP(shortestRun),        REP(interfaceLength),    REP(area),
        REP(perimeter),
    };
    return t;
}

mesh::TMOP::Options blockGridSmoothing() {
    mesh::TMOP::Options t;
    t.metric = mesh::TMOP::ShapeSize007;
    t.quadrature = mesh::TMOP::Corners;
    t.maxSweeps = 1000;
    return t;
}

namespace {

// {key: value} followed by every field of `sources`. Steals `value`.
PyObject *keyedDict(const char *key, PyObject *value, std::initializer_list<Source> sources) {
    PyObject *d = value ? PyDict_New() : nullptr;
    bool ok = d && PyDict_SetItemString(d, key, value) == 0;
    for (const Source &s : sources) ok = ok && addToDict(d, s) == 0;
    Py_XDECREF(value);
    if (ok) return d;
    Py_XDECREF(d);
    return nullptr;
}

}  // namespace

PyObject *meshingDefaults(double h, std::initializer_list<Source> sources) {
    return keyedDict("h", PyFloat_FromDouble(h), sources);
}

PyObject *smoothingDefaults(const mesh::TMOP::Options &defaults, bool pillow) {
    PyObject *d = keyedDict("niters", PyLong_FromLong(defaults.maxSweeps),
                            {source(tmopOptionsTable(), defaults)});
    if (d && PyDict_SetItemString(d, "pillow", pillow ? Py_True : Py_False) < 0) {
        Py_DECREF(d);
        return nullptr;
    }
    return d;
}

void buildSharedTables() {
    nodeOptionsTable();
    qualityTable();
    tmopOptionsTable();
    tmopReportTable();
    pillowReportTable();
    blockQuadMeshOptionsTable();
    blockQuadMeshReportTable();
    shapeDNAOptionsTable();
    shapeDNAReportTable();
    boundaryFeaturesOptionsTable();
    boundaryFeaturesSummaryTable();
}

}  // namespace pycg
