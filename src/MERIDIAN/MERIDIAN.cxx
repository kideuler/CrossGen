#include "MERIDIAN.hxx"

#include <array>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "SIPG/SIPG.hxx"

MERIDIAN::MERIDIAN(std::shared_ptr<Mesh> m) : MERIDIAN(std::move(m), Options()) {}

MERIDIAN::MERIDIAN(std::shared_ptr<Mesh> m, const Options &opts)
    : mesh(std::move(m)), options(opts) {
    if (!mesh) throw std::runtime_error("MERIDIAN: null mesh");
    if (mesh->triangles.empty()) throw std::runtime_error("MERIDIAN: empty mesh");
}

MERIDIAN::~MERIDIAN() = default;

// ---------------------------------------------------------------------------
// runField()
//
// Stage 0 in the paper's terms is the feature-aware triangulation; here the
// mesh arrives already triangulated, so what is left before Sec. 3.1 is the
// cross field the cone indices are read off. Sec. 3.1 uses the frame field of
// [14]; this uses the p=0 DG/SIPG MBO solver, which plays the same role -- a
// boundary-aligned 4-symmetry field whose holonomy is the index.
//
// The same convergence test as TestSIPG: the MBO error is a sum over triangles,
// so the threshold scales with the triangle count.
// ---------------------------------------------------------------------------
void MERIDIAN::runField() {
    field = std::make_unique<SIPG>(mesh, options.sipgMaxSteps, options.sipgGamma);
    // On a multi-material domain the interfaces are Dirichlet data for the field
    // in exactly the way dS is. This is why runField() runs after Stage 0b and
    // not before it: an interface is a curve the layout has to keep, so the
    // field has to be tangent to it, and a field that is not carries none of the
    // cones the regions either side of it need. Empty on a single-material mesh.
    if (interfaces && interfaces->multiMaterial() && options.alignFieldToInterfaces) {
        field->setAlignedInteriorEdges(interfaces->interfaceEdges());
        status.fieldAlignedToInterfaces = true;
    }
    field->initialize();

    const double nTris = static_cast<double>(mesh->triangles.size());
    for (int i = 0; i < options.sipgMaxSteps; ++i) {
        field->step();
        ++status.mboSteps;
        if (field->error < 2.0 * nTris * 1e-5) { status.fieldConverged = true; break; }
    }
    if (!status.fieldConverged) {
        std::ostringstream oss;
        oss << "Cross field did not converge in " << options.sipgMaxSteps
            << " MBO steps (error " << field->error << "); cone placement is unreliable.";
        status.messages.push_back(oss.str());
    }
}

bool MERIDIAN::run() {
    status = Status();
    size_t balanceMessagesSeen = 0;

    // --- Stage 0b: the material interface network -------------------------
    // Before the field, because it says nothing about the field and everything
    // about the domain: which curves the layout has to keep, where they meet,
    // and what the layout has to do where they do. On a single-material mesh it
    // finds nothing and costs one pass over the edges.
    if (options.materialInterfaces) {
        Interfaces::Options iopts;
        iopts.kinkAngle = options.interfaceKinkAngle;
        iopts.loopSplits = options.interfaceLoopSplits;
        interfaces = std::make_unique<Interfaces>(mesh, iopts);
        const Interfaces::Report &fr = interfaces->getReport();
        status.materials = fr.materials;
        status.interfaceEdges = fr.interfaceEdges;
        status.interfaceBranches = fr.branches;
        status.interfaceNodes = fr.nodes;
        status.interfaceIllPosedNodes = fr.illPosedNodes;
        status.interfaceWorstSector = fr.worstSectorResidual;
        for (const std::string &m : fr.messages) status.messages.push_back("Stage 0b: " + m);
        balanceMessagesSeen = fr.messages.size();
    }

    // --- Stage 0: cross field ---------------------------------------------
    runField();

    // --- Stage 1: cone singularities, Sec. 3.1 ----------------------------
    cones = std::make_unique<ConeSingularities>(*field);
    cones->setBoundaryIndexRange(options.minBoundaryIndex, options.maxBoundaryIndex);

    // The interface network's own cone indices, which the field cannot read
    // because it does not know the interfaces are there. Applied before the
    // Gauss-Bonnet check, so that the rebalance below distributes whatever is
    // left over across the boundary cones it is still free to move.
    if (interfaces && interfaces->multiMaterial() && options.prescribeInterfaceCones) {
        cones->prescribe(interfaces->prescription());
        status.interfaceConesPrescribed = interfaces->getReport().prescribedCones;
        status.interfacePrescriptionShift = cones->prescriptionShift();
        if (status.interfacePrescriptionShift > 0) {
            std::ostringstream oss;
            oss << "The interface network moved " << status.interfacePrescriptionShift
                << " index unit(s) from where the cross field put them: the field is "
                << "boundary-aligned and has no boundary condition on an interface, so "
                << "at a junction or a landing it reads whatever propagated in.";
            status.messages.push_back("Stage 1: " + oss.str());
        }
    }

    // The pairs the aligned field puts on a curved interface, which no region
    // asked for. Before the Gauss-Bonnet check because a pair sums to zero and
    // cannot change it, and before Stage 2, which would otherwise cut to each
    // of them.
    if (interfaces && interfaces->multiMaterial() && options.alignFieldToInterfaces &&
        options.cancelInterfaceDipoles) {
        const int nV = static_cast<int>(mesh->vertices.size());
        std::vector<int> region(nV, -1);
        for (int v = 0; v < nV; ++v) region[v] = interfaces->regionAt(v);
        status.coneDipoleUnits = cones->cancelDipoles(region);
        if (status.coneDipoleUnits > 0) {
            std::ostringstream oss;
            oss << "Cancelled " << status.coneDipoleUnits
                << " +1/-1 cone pair(s) the field put inside a single material region: "
                << "an aligned field carries one on the concave side of every strongly "
                << "curved stretch of interface, and a pair inside one region is in no "
                << "region's Gauss-Bonnet count.";
            status.messages.push_back("Stage 1: " + oss.str());
        }
    }

    ConeSingularities::GaussBonnetReport gb = cones->gaussBonnet();
    if (!gb.admissible && options.autoRebalance) {
        const int moved = cones->rebalance();
        if (moved < 0) {
            status.messages.push_back(
                "Could not restore Eq. (4): no boundary cone left within the allowed "
                "index range to take the residual.");
        }
        gb = cones->gaussBonnet();
    }
    status.conesAdmissible = gb.admissible;
    status.interiorCones = static_cast<int>(cones->interiorCones().size());
    status.boundaryCones = static_cast<int>(cones->boundaryCones().size());
    for (const std::string &m : gb.messages) status.messages.push_back("Stage 1: " + m);

    // The discrete Gauss-Bonnet condition is the solvability condition of the
    // Newton system in Stage 3, so an inadmissible set stops the pipeline here
    // rather than producing a metric that means nothing.
    if (!status.conesAdmissible) {
        status.messages.push_back(
            "Stopping before Ricci flow: sum I(v) != 4 chi(S), so the flat cone metric "
            "asked for does not exist.");
        return false;
    }

    // Stage 0b, second pass: the material regions and their own Gauss-Bonnet
    // count. It needs Stage 1's cones on dS -- they are half the count -- which
    // is why it runs here and not with the rest of Stage 0b.
    if (interfaces && interfaces->multiMaterial()) {
        interfaces->balance(cones->getIndices());
        const Interfaces::Report &fr = interfaces->getReport();
        status.interfaceNodes = fr.nodes;
        status.interfaceBranches = fr.branches;
        status.interfaceIllPosedNodes = fr.illPosedNodes;
        status.interfaceWorstSector = fr.worstSectorResidual;
        status.regions = fr.regions;
        status.regionsBalanced = fr.regionsBalanced;
        status.regionQuartersMoved = fr.quartersMoved;
        status.regionCornersInserted = fr.cornersInserted;
        status.worstRegionDeficit = fr.worstRegionDeficit;
        for (size_t i = balanceMessagesSeen; i < fr.messages.size(); ++i) {
            status.messages.push_back("Stage 0b: " + fr.messages[i]);
        }
    }

    // --- Stage 2: cutting graph, Sec. 3.2.2 -------------------------------
    // The interface network goes in with it. Stage 2 is where the layout's two
    // sets of curves first meet: G is free to be routed anywhere and the
    // interfaces are not free at all, so it is G that gives way. See ConeCut's
    // class comment for what a cut that does not costs E3 and E6.
    ConeCut::Options cutOpts;
    cutOpts.conesToBoundary = options.coneCutsToBoundary;
    cutOpts.interfaceAvoidance = options.coneCutInterfaceAvoidance;
    cutter = std::make_unique<ConeCut>(mesh, *cones, cutOpts, interfaces.get());
    const ConeCut::Report &cr = cutter->getReport();
    status.cutIsDisk = cr.isDisk;
    status.allConesOnBoundary = cr.allConesOnBoundary;
    status.cutInterfaceEdges = cr.interfaceEdgesOnCut;
    status.cutInterfaceVertices = cr.interfaceVertsOnCut;
    status.cutInterfaceNodes = cr.interfaceNodesOnCut;
    for (const std::string &m : cr.messages) status.messages.push_back("Stage 2: " + m);

    // --- Stage 3: discrete Ricci flow, Sec. 3.2.1 -------------------------
    // Runs on the *uncut* triangulation: a conformal factor per vertex has
    // nothing to say about seams, and Omega is what Stage 4 lays out, not what
    // the flow acts on.
    ricci = std::make_unique<RicciFlow>(mesh, *cones);
    ricci->setTolerance(options.ricciTolerance);
    ricci->setMaxIterations(options.ricciMaxIterations);
    ricci->setDelaunayFlips(options.delaunayFlips);
    status.ricciConverged = ricci->solve();
    for (const std::string &m : ricci->getReport().messages) {
        status.messages.push_back("Stage 3: " + m);
    }

    status.readyForImmersion = status.conesAdmissible && status.cutIsDisk &&
                               status.allConesOnBoundary && status.ricciConverged &&
                               ricci->getReport().nonRealisableFaces == 0;
    if (!status.readyForImmersion) {
        status.messages.push_back(
            "Stopping before the immersion: the flat cone metric on a disk that Stage 4 "
            "unfolds is not there yet.");
        return false;
    }

    // --- Stage 4: metric immersion, Sec. 3.2.2 ----------------------------
    immersion = std::make_unique<Immersion>(*cutter, *ricci, *cones);
    const Immersion::Report &ir = immersion->getReport();
    status.immersionValid = ir.valid;
    status.immersionFlippedFaces = ir.flippedFaces;
    status.seamArcs = ir.arcs;
    for (const std::string &m : ir.messages) status.messages.push_back("Stage 4: " + m);

    // A fold here is terminal for the rest of the pipeline, not merely
    // untidy: E1's barrier keeps local injectivity, it does not restore it.
    if (ir.flippedFaces > 0 || ir.unplacedVertices > 0) {
        status.messages.push_back(
            "Stopping before the layout energies: psi_R does not satisfy Q1, and Sec. 3.3's "
            "barrier can only preserve Q1, never repair it.");
        return false;
    }
    if (!options.runLayout) return true;

    // --- Stage 5: subdomain labelling, Sec. 3.3 ---------------------------
    SubdomainLabels::Options lopts;
    lopts.seedTopoConstraints = options.seedTopoConstraints;
    lopts.nearMissTolerance = options.topoNearMiss;
    lopts.seedSelfReturns = options.seedSelfReturns;
    lopts.seedAllConnections = options.seedAllConnections;
    lopts.maxTraceSteps = options.separatrixMaxSteps;
    lopts.interfaceCorners = options.interfaceCorners;
    lopts.propagateInterfaceLabels = options.propagateInterfaceLabels;
    labels = std::make_unique<SubdomainLabels>(*immersion, lopts, interfaces.get());
    const SubdomainLabels::Report &lr = labels->getReport();
    status.boundaryEdgesU = lr.boundaryEdgesU;
    status.boundaryEdgesV = lr.boundaryEdgesV;
    status.featureChains = lr.featureChains;
    status.interfaceCorners = lr.interfaceCorners;
    status.interfaceLabelsCorrected = lr.featureLabelsCorrected;
    status.topoPaths = lr.topoPaths;
    status.topoSelfReturns = lr.topoSelfReturns;
    status.topoExtraPerPair = lr.topoExtraPerPair;
    size_t labelMessagesSeen = 0;
    auto drainLabelMessages = [&]() {
        const auto &msgs = labels->getReport().messages;
        for (size_t i = labelMessagesSeen; i < msgs.size(); ++i) {
            status.messages.push_back("Stage 5: " + msgs[i]);
        }
        labelMessagesSeen = msgs.size();
    };
    drainLabelMessages();

    // --- Stage 6: the layout-inducing energies, Sec. 3.3 ------------------
    LayoutEnergy::Options eopts;
    eopts.lambdaInit = options.lambdaInit;
    eopts.lambdaGrowth = options.lambdaGrowth;
    eopts.outerSteps = options.outerSteps;
    eopts.innerIterations = options.innerIterations;
    eopts.alternateReference = options.alternateReference;
    eopts.relabel = options.relabelBetweenSteps;
    eopts.lagInterfaceScales = options.lagInterfaceScales;
    layout = std::make_unique<LayoutEnergy>(*immersion, *labels, eopts);
    status.layoutRan = layout->run();
    const LayoutEnergy::Report &er = layout->getReport();
    status.layoutInjective = er.injective;
    status.layoutConstrained = er.constraintsMet;
    status.outerStepsTaken = er.outerSteps;
    status.layoutValid = er.valid;
    status.interfaceResidual = er.maxInterfaceResidual;
    status.interfaceCornerChanges = er.interfaceCornerChanges;
    status.interfacesAligned = er.interfaceCornerChanges == 0 &&
                               er.maxInterfaceResidual < 1e-3;
    size_t layoutMessagesSeen = 0;
    auto drainLayoutMessages = [&]() {
        const auto &msgs = layout->getReport().messages;
        for (size_t i = layoutMessagesSeen; i < msgs.size(); ++i) {
            status.messages.push_back("Stage 6: " + msgs[i]);
        }
        layoutMessagesSeen = msgs.size();
    };
    drainLayoutMessages();

    if (!options.runSeparatrices) return status.layoutValid;

    // --- Stage 7 and Sec. 3.3's repair ------------------------------------
    Separatrices::Options sopts;
    sopts.coneSnapTolerance = options.separatrixSnap;
    sopts.maxSteps = options.separatrixMaxSteps;
    sopts.coneSnapRings = options.separatrixSnapRings;
    sopts.detectCycles = options.separatrixDetectCycles;
    sopts.nearMissWindow = options.repairGapLimit;
    // The interface network's nodes emit too -- the ones that are singular
    // points of the layout. See Interfaces::emitterNodes.
    if (interfaces && interfaces->multiMaterial()) {
        sopts.extraEmitters = interfaces->emitterNodes();
    }

    RepairOptions ropts;
    // With the seeding switched off, E5 has no constraints at all and the run
    // is Remark 3.1's demonstration -- the bracket at 1898 patches instead of
    // 73. Letting the repair put the constraints back would quietly undo it.
    ropts.passes = options.seedTopoConstraints ? options.repairPasses : 0;
    ropts.maxPerPass = options.repairMaxPerPass;
    ropts.gapLimit = options.repairGapLimit;
    ropts.lambdaBoost = options.repairLambdaBoost;
    ropts.outerSteps = options.repairOuterSteps;
    ropts.scoreArrangement = options.runArrangement && options.repairScoresArrangement;
    ropts.patience = options.repairPatience;
    ropts.arrangement.mergeTolerance = options.arrangementMerge;
    ropts.arrangement.cornerTolerance = options.arrangementCorner;
    ropts.arrangement.collapseTolerance = options.arrangementCollapse;
    ropts.arrangement.trimUnresolvedAtCrossings = options.arrangementTrim;

    RepairResult rep = traceAndRepair(*labels, *layout, sopts, ropts,
                                      [&](const std::string &m) {
                                          status.messages.push_back(m);
                                      });
    separatrices = std::move(rep.separatrices);
    status.repairPasses = rep.passesTaken;
    status.repairConstraintsAdded = rep.constraintsAdded;
    status.topoPaths = labels->getReport().topoPaths;
    status.topoSelfReturns = labels->getReport().topoSelfReturns;
    status.topoExtraPerPair = labels->getReport().topoExtraPerPair;

    const LayoutEnergy::Report &fr = layout->getReport();
    status.layoutInjective = fr.injective;
    status.layoutConstrained = fr.constraintsMet;
    status.outerStepsTaken = fr.outerSteps;
    status.layoutValid = fr.valid;

    if (separatrices) {
        const Separatrices::Report &sr = separatrices->getReport();
        status.separatricesRan = true;
        status.separatrices = sr.emitted;
        status.separatricesToCone = sr.endedAtCone;
        status.separatricesToBoundary = sr.endedAtBoundary;
        status.separatricesUnresolved = sr.capped + sr.cycled + sr.stuck + sr.degenerate;
        status.separatricesNearMisses = sr.nearMisses + sr.grazes;
        // Q5 is the property the paper states -- every integral curve is finite
        // and ends at a cone or at dS -- and a curve that grazed a cone on its
        // way out through the boundary has satisfied it. The near misses are
        // reported beside it rather than folded into it: they are a statement
        // about the *quality* of the layout (Remark 3.1's slivers) and not
        // about its validity, and the repair above is what acts on them.
        status.q5Verified = sr.valid;
    }

    // --- Stage 8: the arrangement, Sec. 4 ---------------------------------
    // Run even when Stage 6 fell short. A layout that does not satisfy Q5 still
    // has an arrangement, and the arrangement is where the failure becomes a
    // *place* -- this face has three corners, that cone is one arc short --
    // rather than a residual. Sec. 4's own "on failure" list is written against
    // exactly this output.
    if (!separatrices || !options.runArrangement) return status.layoutValid;
    Arrangement::Options aopts;
    aopts.mergeTolerance = options.arrangementMerge;
    aopts.cornerTolerance = options.arrangementCorner;
    aopts.collapseTolerance = options.arrangementCollapse;
    aopts.trimUnresolvedAtCrossings = options.arrangementTrim;
    try {
        arrangement = std::make_unique<Arrangement>(*separatrices, *labels, aopts);
    } catch (const std::exception &e) {
        status.messages.push_back(std::string("Stage 8: could not be built: ") + e.what());
        return status.layoutValid;
    }
    const Arrangement::Report &ar = arrangement->getReport();
    status.arrangementRan = true;
    status.layoutNodes = ar.nodes;
    status.layoutArcs = ar.arcs;
    status.layoutPatches = ar.patches;
    status.layoutQuads = ar.simpleQuads;
    status.layoutCoverage = ar.areaCoverage;
    status.arrangementValid = ar.valid;
    for (const std::string &m : ar.messages) status.messages.push_back("Stage 8: " + m);

    // --- Stage 9: the spline reconstruction, Sec. 5 -----------------------
    if (!options.runSplines) return status.layoutValid;
    SplineFit::Options fopts;
    fopts.segments = options.splineSegments;
    fopts.samples = options.splineSamples;
    fopts.fitBoundaryArcs = options.fitBoundaryArcs;
    fopts.fitInterfaceArcs = options.fitInterfaceArcs;
    try {
        splines = std::make_unique<SplineFit>(*arrangement, fopts);
    } catch (const std::exception &e) {
        status.messages.push_back(std::string("Stage 9: could not be built: ") + e.what());
        return status.layoutValid;
    }
    const SplineFit::Report &fr2 = splines->getReport();
    status.splinesRan = true;
    status.splinePatches = fr2.patches;
    status.splineControlPoints = fr2.controlPointsPerArc;
    status.splineMaxDeviation = fr2.maxDeviation;
    status.splinesWatertight = fr2.watertight;
    status.splinesValid = fr2.valid;
    for (const std::string &m : fr2.messages) status.messages.push_back("Stage 9: " + m);

    // --- Stage 10: the quadrilateral mesh -------------------------------
    // Run on whatever Stage 8 closed. A layout with three broken patches still
    // meshes everywhere else, and the fraction of S left uncovered is a more
    // useful statement of the damage than refusing to mesh at all.
    if (!options.runQuadMesh) return status.layoutValid;
    QuadMesh::Options qopts;
    qopts.targetEdgeLength = options.quadTargetEdge;
    qopts.minIntervals = options.quadMinIntervals;
    qopts.maxIntervals = options.quadMaxIntervals;
    qopts.useSplines = options.quadUseSplines;
    qopts.featuresOnTracedArcs = options.quadFeaturesOnTracedArcs;
    qopts.smoothingPasses = options.quadSmoothingPasses;
    qopts.smoothingThreshold = options.quadSmoothingThreshold;
    try {
        quads = std::make_unique<QuadMesh>(*splines, qopts);
    } catch (const std::exception &e) {
        status.messages.push_back(std::string("Stage 10: could not be built: ") + e.what());
        return status.layoutValid;
    }
    const QuadMesh::Report &qr = quads->getReport();
    status.quadMeshRan = true;
    status.meshVertices = qr.vertices;
    status.meshQuads = qr.quads;
    status.meshChords = qr.chords;
    status.meshUnmeshedPatches = qr.unmeshedPatches;
    status.meshMinScaledJacobian = qr.minScaledJacobian;
    status.meshConforming = qr.conforming;
    status.meshValid = qr.valid;
    for (const std::string &m : qr.messages) status.messages.push_back("Stage 10: " + m);

    return status.layoutValid;
}

// ---------------------------------------------------------------------------
// traceAndRepair()  --  Stage 7, then Sec. 3.3's remedy applied in a loop
//
//     trace the separatrices of Psi
//     for each that ended at neither a cone nor dS, and each that slipped past
//         a cone on its way out through dS, add the Gamma_topo constraint it
//         names -- nearest miss first, a few at a time
//     raise lambda_5 and re-run Stage 6 *from the current phi*
//     trace again
//
// From the current phi, not from psi_R. The distinction is the paper's and it
// is not a matter of speed: psi_R satisfies Q1, Q2 and Q4 and nothing else, so
// restarting there throws away every level of lambda already paid for and
// re-enters the same local minimum by the same road. Continuing from phi enters
// it with the new constraint already switched on, which is the only thing that
// has changed.
//
// A round is a *different* minimisation, not a continuation of the same one --
// the constraint set has changed -- so it is not guaranteed to land better than
// where it started, and on a model whose cones Stage 1 left clustered it
// sometimes does not. The map and the paths added are therefore both taken back
// when a round comes out worse, which makes the loop monotone: more rounds
// cannot leave a worse result than fewer, only a slower one. That is what
// allows the repair to be on by default.
//
// Better means, in order: not inverted (a map that folded is never an
// improvement, whatever its residuals); then fewer separatrices that terminated
// nowhere, which is Q5 itself; then fewer that grazed a cone on their way out,
// which is Remark 3.1's sliver count. A tie goes to the newer map, because it
// is the one with more of Gamma_topo switched on.
//
// Static, and taking the two stages by reference rather than reading members,
// so that the viewer -- which drives Stages 4 to 7 itself, one keypress at a
// time -- runs this loop and not a second copy of it.
// ---------------------------------------------------------------------------
MERIDIAN::RepairResult MERIDIAN::traceAndRepair(SubdomainLabels &labels, LayoutEnergy &layout,
                                                const Separatrices::Options &trace,
                                                const RepairOptions &opts,
                                                const std::function<void(const std::string &)> &log) {
    RepairResult out;

    size_t labelMessagesSeen = labels.getReport().messages.size();
    size_t layoutMessagesSeen = layout.getReport().messages.size();
    auto drainLabels = [&]() {
        const auto &m = labels.getReport().messages;
        for (size_t i = labelMessagesSeen; i < m.size(); ++i) log("Stage 5: " + m[i]);
        labelMessagesSeen = m.size();
    };
    auto drainLayout = [&]() {
        const auto &m = layout.getReport().messages;
        for (size_t i = layoutMessagesSeen; i < m.size(); ++i) log("Stage 6: " + m[i]);
        layoutMessagesSeen = m.size();
    };

    // Stage 7 is run even when the continuation fell short. The curves are the
    // only place Q5 is visible as the property it actually is -- every integral
    // curve out of a cone is finite -- and one that does not terminate names
    // the pair of cones whose connectivity constraint is missing, which is the
    // diagnosis Sec. 3.3 asks for when the layout is not yet a layout.
    auto retrace = [&]() -> bool {
        try {
            out.separatrices = std::make_unique<Separatrices>(layout.getImmersion(),
                                                              layout.getUV(), trace);
        } catch (const std::exception &e) {
            out.separatrices.reset();
            log(std::string("Stage 7: could not be run: ") + e.what());
            return false;
        }
        return true;
    };
    if (!retrace()) return out;
    out.traced = true;

    // What "better" means, in order: not inverted -- a map that folded is never
    // an improvement, whatever its residuals; then Q5 itself, the curves that
    // terminated nowhere; then how far the arrangement is from being a set of
    // patches, which is what a caller actually wanted; then Remark 3.1's sliver
    // count, the curves that grazed a cone on their way out.
    //
    // Q5 stays ahead of the arrangement rather than the other way round, and
    // the corpus is what settles it: scored the other way the loop will take a
    // round that mends one broken patch at the price of two more curves that
    // wind, because a curve that winds is only one defect in Stage 8 -- the
    // T-junction trimUnresolvedAtCrossings ends it at -- while being a much
    // worse thing to have. The arrangement is what separates rounds that Q5
    // cannot tell apart, which is the job it was added for.
    auto quality = [&]() {
        const Separatrices::Report &r = out.separatrices->getReport();
        int defects = 0;
        if (opts.scoreArrangement) {
            try {
                Arrangement ar(*out.separatrices, labels, opts.arrangement);
                const Arrangement::Report &a = ar.getReport();
                defects = a.wrongCornerFaces + a.tJunctions + a.coneValenceErrors +
                          a.danglingNodes + a.unsharedArcs;
            } catch (const std::exception &) {
                defects = 1 << 20;   // an arrangement that would not build at all
            }
        }
        return std::array<int, 4>{layout.getReport().invertedTriangles,
                                  r.capped + r.cycled + r.stuck + r.degenerate, defects,
                                  r.nearMisses + r.grazes};
    };

    std::array<int, 4> best = quality();
    std::vector<Point> bestMap = layout.saveMap();
    size_t bestPaths = labels.topoPaths().size();
    bool atBest = true;
    int stale = 0;

    for (int pass = 1; pass <= opts.passes; ++pass) {
        if (best[1] == 0 && best[2] == 0 && best[3] == 0) break;

        const int added = labels.adoptCurves(*out.separatrices, opts.gapLimit, pass,
                                             opts.maxPerPass);
        if (added == 0) {
            // Nothing new to say. Either every unterminated curve already has
            // its constraint -- in which case the continuation is what fell
            // short, not the constraint set -- or none of them came near enough
            // to a cone for the pair to have been meant to join at all.
            break;
        }
        drainLabels();

        layout.resume(opts.outerSteps, opts.lambdaBoost);
        drainLayout();

        const bool traced = retrace();
        const std::array<int, 4> now =
            traced ? quality() : std::array<int, 4>{1 << 20, 1 << 20, 1 << 20, 1 << 20};

        out.passesTaken = pass;
        out.constraintsAdded += added;

        if (now < best) {
            best = now;
            bestMap = layout.saveMap();
            bestPaths = labels.topoPaths().size();
            atBest = true;
            stale = 0;
        } else {
            // Not an improvement. That is not by itself a reason to stop: a
            // round adds the constraint for one pair of cones, and a pair that
            // needs two constraints to be worth anything -- the second cone of
            // a cluster, the far side of a seam -- costs before it pays. What
            // is a reason to stop is `patience` of them in a row, by which
            // point the constraint set is being added to and nothing is coming
            // of it.
            atBest = false;
            ++stale;
            if (!traced || stale >= std::max(1, opts.patience)) break;
        }
    }

    // Back to the best round the loop saw, which is what makes it monotone:
    // more rounds cannot leave a worse result than fewer, only a slower one.
    // That is what allows the repair to be on by default.
    //
    // Order matters. The paths go first so that rebuildConstraints() assembles
    // Eqs. (15) to (19) over the set that is being kept; loadMap() goes last so
    // that the residuals and the verdict in the Report are the ones of the map
    // actually being returned.
    if (!atBest) {
        std::ostringstream oss;
        oss << "Stage 6: the last " << stale << " repair pass(es) did not improve on what "
            << "they found (" << best[1] << " separatrix/ces terminating nowhere, " << best[2]
            << " defect(s) in the arrangement, " << best[3] << " grazing a cone), so their "
            << "constraints and the maps they produced were taken back. What stands is the "
            << "last round that improved.";
        log(oss.str());

        labels.truncateTopoPaths(bestPaths);
        layout.rebuildConstraints();
        layout.loadMap(bestMap);
        ++out.rolledBack;
        retrace();
    }

    // The Stage 7 diagnostics of the map that is actually being returned, and
    // of no other: the intermediate tracings are steps in a search, and logging
    // each of them would bury the one that describes the result.
    if (out.separatrices) {
        for (const std::string &m : out.separatrices->getReport().messages) {
            log("Stage 7: " + m);
        }
    }
    return out;
}
