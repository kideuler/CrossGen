#include "TORSION.hxx"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "SIPG/SIPG.hxx"

TORSION::TORSION(std::shared_ptr<Mesh> m) : TORSION(std::move(m), Options()) {}

TORSION::TORSION(std::shared_ptr<Mesh> m, const Options &opts)
    : mesh(std::move(m)), options(opts) {
    if (!mesh) throw std::runtime_error("TORSION: null mesh");
    if (mesh->triangles.empty()) throw std::runtime_error("TORSION: empty mesh");
}

TORSION::~TORSION() = default;

// ---------------------------------------------------------------------------
// runField()  --  Stage 0, MERIDIAN's unchanged
// ---------------------------------------------------------------------------
void TORSION::runField() {
    field = std::make_unique<SIPG>(mesh, options.sipgMaxSteps, options.sipgGamma);
    // On a multi-material domain the interfaces are Dirichlet data for the
    // field in exactly the way dS is, and that matters more here than in
    // Pipeline A rather than less: this pipeline integrates the field, so an
    // interface the field ran straight through is an interface the *map* runs
    // straight through. See SIPG::setAlignedInteriorEdges for why the hard pin
    // and not a penalty.
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
            << " MBO steps (error " << field->error
            << "); every stage of this pipeline is downstream of it.";
        status.messages.push_back("Stage 0: " + oss.str());
    }
}

// ---------------------------------------------------------------------------
// runFront()  --  Stages 0b, 0, 1 and 2
//
// Identical to MERIDIAN's, deliberately and line for line: the interface
// network, the field, the cone indices with their prescriptions and dipole
// cancellations, the Gauss-Bonnet gate, the per-region balance, and the cut.
// Nothing in any of it depends on how psi_0 will be built, and the plan says so
// -- Stages 0b, 1 and 2 stay as they are -- so the only thing worth doing here
// is to keep it recognisably the same code.
// ---------------------------------------------------------------------------
bool TORSION::runFront() {
    size_t balanceMessagesSeen = 0;

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

    runField();

    cones = std::make_unique<ConeSingularities>(*field);
    cones->setBoundaryIndexRange(options.minBoundaryIndex, options.maxBoundaryIndex);
    // The snapshot Sec. 5.1's audit is run against, taken here because this is
    // the last moment at which the indices are the ones the *field* read.
    // prescribe(), cancelDipoles() and rebalance() all move them on purpose
    // over the next few lines, and an audit run against the moved set would be
    // reporting the pipeline doing its job.
    fieldIndex = cones->getIndices();

    if (interfaces && interfaces->multiMaterial() && options.prescribeInterfaceCones) {
        cones->prescribe(interfaces->prescription());
        status.interfaceConesPrescribed = interfaces->getReport().prescribedCones;
        status.interfacePrescriptionShift = cones->prescriptionShift();
        if (status.interfacePrescriptionShift > 0) {
            std::ostringstream oss;
            oss << "The interface network moved " << status.interfacePrescriptionShift
                << " index unit(s) from where the cross field put them.";
            status.messages.push_back("Stage 1: " + oss.str());
        }
    }

    if (interfaces && interfaces->multiMaterial() && options.alignFieldToInterfaces &&
        options.cancelInterfaceDipoles) {
        const int nV = static_cast<int>(mesh->vertices.size());
        std::vector<int> region(nV, -1);
        for (int v = 0; v < nV; ++v) region[v] = interfaces->regionAt(v);
        status.coneDipoleUnits = cones->cancelDipoles(region);
        if (status.coneDipoleUnits > 0) {
            std::ostringstream oss;
            oss << "Cancelled " << status.coneDipoleUnits
                << " +1/-1 cone pair(s) the field put inside a single material region.";
            status.messages.push_back("Stage 1: " + oss.str());
        }
    }

    ConeSingularities::GaussBonnetReport gb = cones->gaussBonnet();
    if (!gb.admissible && options.autoRebalance) {
        const int moved = cones->rebalance();
        if (moved < 0) {
            status.messages.push_back(
                "Stage 1: could not restore Eq. (4): no boundary cone left within the allowed "
                "index range to take the residual.");
        }
        gb = cones->gaussBonnet();
    }
    status.conesAdmissible = gb.admissible;
    status.interiorCones = static_cast<int>(cones->interiorCones().size());
    status.boundaryCones = static_cast<int>(cones->boundaryCones().size());
    for (const std::string &m : gb.messages) status.messages.push_back("Stage 1: " + m);

    // Eq. (4) is the solvability condition of the Ricci flow's Newton system in
    // Pipeline A. There is no Newton system here -- and it is no less binding:
    // sum I(v) = 4 chi(S) is the condition for a *seamless map with these
    // holonomies to exist at all*, which is the same statement the flow's
    // consistency condition was standing in for.
    if (!status.conesAdmissible) {
        status.messages.push_back(
            "Stopping before the cut: sum I(v) != 4 chi(S), so no seamless map with these "
            "holonomies exists and there is nothing for the integration to find.");
        return false;
    }

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

    ConeCut::Options cutOpts;
    cutOpts.conesToBoundary = options.coneCutsToBoundary;
    cutOpts.interfaceAvoidance = options.coneCutInterfaceAvoidance;
    cutter = std::make_unique<ConeCut>(mesh, *cones, cutOpts, interfaces.get());
    const ConeCut::Report &cr = cutter->getReport();
    status.cutIsDisk = cr.isDisk;
    status.allConesOnBoundary = cr.allConesOnBoundary;
    status.cutInteriorJunctions = cr.interiorJunctions;
    status.cutInterfaceEdges = cr.interfaceEdgesOnCut;
    status.cutInterfaceVertices = cr.interfaceVertsOnCut;
    status.cutInterfaceNodes = cr.interfaceNodesOnCut;
    for (const std::string &m : cr.messages) status.messages.push_back("Stage 2: " + m);

    if (!status.cutIsDisk || !status.allConesOnBoundary) {
        status.messages.push_back(
            "Stopping before the combing: Omega is not a disk with every cone on its boundary, "
            "so the branch of the field over it is not single valued and the BFS that selects "
            "it has no reason to be path-independent.");
        return false;
    }
    return true;
}

// ---------------------------------------------------------------------------
// inducedLengths()
//
// The lengths a map on Omega induces, re-indexed by the edges of S, which is the
// indexing LayoutEnergy::buildReference() reads. A seam edge has two images in
// Omega and they are the same length to the seam residual, because the
// transition holding them together is a rotation; the first found is taken.
//
// This is E1's reference metric in this pipeline, and the reason it is one is
// in the class comment: the field frame cannot be, because a rotation-invariant
// energy cannot see the rotation that carries the cones, and the model's own
// Euclidean geometry has no cones at all. A map that satisfies the seam
// constraints has them exactly, and that is the whole of what is needed.
// ---------------------------------------------------------------------------
std::vector<double> TORSION::inducedLengths(const Mesh &mesh, const ConeCut &cut,
                                            const std::vector<Point> &psi) {
    const Mesh &om = cut.getCutMesh();
    const auto &c2o = cut.getCutVertexToOriginal();
    std::vector<double> len(mesh.edges.size(), 0.0);
    if (psi.size() != om.vertices.size()) return len;

    std::unordered_map<MeshEdgeKey, int, MeshEdgeKeyHash> origEdge;
    origEdge.reserve(mesh.edges.size() * 2);
    for (int e = 0; e < static_cast<int>(mesh.edges.size()); ++e) {
        origEdge.emplace(MeshEdgeKey(mesh.edges[e][0], mesh.edges[e][1]), e);
    }
    for (int e = 0; e < static_cast<int>(om.edges.size()); ++e) {
        auto it = origEdge.find(MeshEdgeKey(c2o[om.edges[e][0]], c2o[om.edges[e][1]]));
        if (it == origEdge.end() || len[it->second] > 0.0) continue;
        len[it->second] = normP(psi[om.edges[e][1]] - psi[om.edges[e][0]]);
    }
    // An edge with no positive image -- a triangle the integration collapsed --
    // falls back to its Euclidean length, which is the reference having no cone
    // information about that one edge rather than having none at all.
    for (size_t e = 0; e < len.size(); ++e) {
        if (!(len[e] > 0.0)) {
            len[e] = normP(mesh.vertices[mesh.edges[e][1]] -
                           mesh.vertices[mesh.edges[e][0]]);
        }
    }
    return len;
}

// ---------------------------------------------------------------------------
// untangle()  --  Stage 4R, Sec. 7.2
//
//   1. Tutte embedding of Omega: bijective by theorem, so zero flips, and
//      useless for anything except being somewhere injective to start from.
//   2. min_phi sum_t A_t [ ||J (J*)^-1||_F^2 + ||J* J^-1||_F^2 ] + mu E4(phi)
//   3. raise mu over a few outer steps.
//
// Step 2 is not a new solver: E1 of Stage 6 is already the bracket above, E4 is
// already the seam term, and the closed-form flip cap in its line search is what
// makes the pass work where the least-squares solve does not. The map it returns
// is locally injective by construction of the step, whether or not the target it
// is fitting is realisable --
//
//     Fitting a map to a per-triangle target Jacobian under a barrier is well
//     posed *even when the target is unrealisable*. Non-integrability makes the
//     target inconsistent, not the fit ill-posed.
//
// -- so the pass is LayoutEnergy with lambda_2, lambda_3, lambda_5 and lambda_6
// switched off, started from the Tutte map instead of from psi_R, and step 3 is
// its ordinary penalty schedule.
//
// **Which reference.** The plan says Reference::Field here, on the argument that
// J* is the target and E1 against J* is the fit. It is not: J* is a scaled
// rotation, so E1 against it is E1 against the model's Euclidean geometry --
// see LayoutEnergy.hxx for the algebra and TestTORSION --ref-test for the
// measurement -- and that reference has no cones, which is the failure this
// whole pipeline has to avoid. What is used instead is the metric the *least-
// squares map* induces. It is the closest thing to the field's own answer that
// a rotation-invariant energy can be given: it carries the cone angles, because
// the seam constraints put them there exactly, and it carries the shape the
// integration settled on everywhere else. The orientation the frame also
// carries is not recoverable by E1 at all, whatever it is composed with, and
// where that matters -- Stage 7's cross-check that a cone's rays leave along
// the field's directions -- it has to be checked rather than energetically
// asked for.
//
// The lengths are read off the map even where it inverted: a flipped triangle
// still has three positive side lengths, they are the ones the integration
// produced, and there are single figures of them.
//
// Returns an empty vector when there was nothing to start from.
// ---------------------------------------------------------------------------
std::vector<Point> TORSION::untangle() {
    const Mesh &omega = cutter->getCutMesh();

    TutteEmbedding::Options topts;
    topts.targetEdge = options.targetEdge;
    tutte = std::make_unique<TutteEmbedding>(omega, topts);
    const TutteEmbedding::Report &tr = tutte->getReport();
    status.tutteValid = tr.valid;
    for (const std::string &m : tr.messages) status.messages.push_back("Stage 4R: " + m);
    if (!tr.valid) {
        status.messages.push_back(
            "Stage 4R: the untangling has nothing locally injective to start from, so it was "
            "not run; psi_0 stands as the integration returned it.");
        return {};
    }

    // A second Immersion over the Tutte map, for the arcs and the seam pairing
    // the energy needs. The combed angles go in again, so the k it holds are
    // the matchings' and not a Procrustes fit of a Tutte map, which would be
    // meaningless.
    std::unique_ptr<Immersion> start;
    try {
        start = std::make_unique<Immersion>(*cutter, *cones, tutte->getUV(),
                                            frames->fieldEdgeLengths(),
                                            frames->combedAngle());
    } catch (const std::exception &e) {
        status.messages.push_back(std::string("Stage 4R: could not wrap the Tutte map: ") + e.what());
        return {};
    }

    // No boundary labelling worth having on a Tutte map and no connectivity to
    // seed from it, so neither is asked for. With lambda_2, lambda_3 and
    // lambda_5 at zero none of it would be read in any case; not building it
    // saves a separatrix trace of a map that means nothing.
    SubdomainLabels::Options lopts;
    lopts.seedTopoConstraints = false;
    lopts.interfaceCorners = false;
    SubdomainLabels startLabels(*start, lopts);

    LayoutEnergy::Options eopts;
    eopts.reference = LayoutEnergy::Reference::Induced;
    eopts.referenceLengths = inducedLengths(*mesh, *cutter, integratedMap);
    // mu is E4's penalty and nothing else is switched on, so it is set through
    // lambdaFactor rather than through lambdaInit: run() floors lambda_4 at
    // lambda_1 (Q4 is not free to trade away in the main continuation), and
    // that floor would swallow a starting mu below 1. The schedule then raises
    // it by lambdaGrowth per outer step, which is step 3 of Sec. 7.2.
    eopts.lambdaInit = 1.0;
    eopts.lambdaFactor[4] = options.untangleSeamWeight;
    eopts.lambdaGrowth = options.lambdaGrowth;
    eopts.outerSteps = options.untangleOuterSteps;
    eopts.innerIterations = options.untangleInnerIterations;
    eopts.alternateReference = false;
    eopts.relabel = false;
    // E1 and E4 only. Everything else is a statement about a layout and this
    // pass is not producing one -- it is producing something injective with the
    // field's directions in it, for Stage 6 to make a layout out of.
    eopts.lambdaFactor[2] = 0.0;
    eopts.lambdaFactor[3] = 0.0;
    eopts.lambdaFactor[5] = 0.0;
    eopts.lambdaFactor[6] = 0.0;


    LayoutEnergy fit(*start, startLabels, eopts);
    fit.run();
    const LayoutEnergy::Report &fr = fit.getReport();
    status.untangleRan = true;
    status.untangleOuterSteps = fr.outerSteps;
    status.untangleFlippedFaces = fr.invertedTriangles;
    status.untangleSeamResidual = fr.maxSeamResidual;
    status.layoutReferenceLeftHanded =
        std::max(status.layoutReferenceLeftHanded, fr.leftHandedFrames);
    for (const std::string &m : fr.messages) status.messages.push_back("Stage 4R: " + m);

    if (fr.invertedTriangles > 0) {
        std::ostringstream oss;
        oss << "The untangling pass came back with " << fr.invertedTriangles
            << " inverted triangle(s), which means the line search's flip cap was defeated -- "
            << "normally by a reference triangle that was already degenerate. Sec. 7.3's "
            << "knobs, in order: field smoothness near the cones, refinement where the "
            << "per-triangle fit residual is largest, a lower h there.";
        status.messages.push_back("Stage 4R: " + oss.str());
    } else {
        std::ostringstream oss;
        oss << "Untangled: from the Tutte embedding, " << fr.outerSteps
            << " outer step(s) of target-Jacobian fitting under the barrier left every "
            << "triangle positively oriented, with the seam " << fr.maxSeamResidual
            << " of the extent from exact. E4 is what Stage 6 raises from here.";
        status.messages.push_back("Stage 4R: " + oss.str());
    }
    return fit.getUV();
}

// ---------------------------------------------------------------------------
bool TORSION::run() {
    status = Status();

    if (!runFront()) return false;

    // --- Stage 3F: the frames, the matchings and the audit (Sec. 5) -------
    FieldFrames::Options fopts;
    fopts.targetEdge = options.targetEdge;
    fopts.sizing = options.sizing;
    fopts.referenceIndex = fieldIndex;
    try {
        frames = std::make_unique<FieldFrames>(*field, *cutter, *cones, fopts);
    } catch (const std::exception &e) {
        status.messages.push_back(std::string("Stage 3F: could not comb the field: ") + e.what());
        return false;
    }
    const FieldFrames::Report &ffr = frames->getReport();
    status.framesRan = true;
    status.combedFaces = ffr.combedFaces;
    status.unreachedFaces = ffr.unreachedFaces;
    status.combingDefects = ffr.combingDefects;
    status.maxFrameJump = ffr.maxFrameJump;
    status.indexMismatches = ffr.indexMismatches;
    status.clusteredConePairs = ffr.clusteredConePairs;
    status.highIndexCones = ffr.highIndexCones;
    status.leftHandedFrames = ffr.leftHandedFrames;
    status.maxMetricDisagreement = ffr.maxMetricDisagreement;
    status.framesValid = ffr.valid;
    for (const std::string &m : ffr.messages) status.messages.push_back("Stage 3F: " + m);

    // Sec. 13's spot check. Two combs of a disk differ by the constant the two
    // seeds differ by and by nothing else, so the *differences* agree exactly;
    // an integer that varies is the BFS having taken a route the other did not
    // and got a different answer for it.
    if (options.checkSecondSeed && ffr.unreachedFaces == 0 && ffr.faces > 1) {
        FieldFrames::Options second = fopts;
        second.seedFace = ffr.faces / 2;
        try {
            FieldFrames other(*field, *cutter, *cones, second);
            const int shift = other.branch()[ffr.seedFace] - frames->branch()[ffr.seedFace];
            for (int t = 0; t < ffr.faces; ++t) {
                if (other.branch()[t] - frames->branch()[t] != shift) {
                    status.secondSeedAgrees = false;
                    break;
                }
            }
        } catch (const std::exception &) {
            status.secondSeedAgrees = false;
        }
        if (!status.secondSeedAgrees) {
            status.messages.push_back(
                "Stage 3F: combing from a second seed gave a different branch. On a disk with "
                "every cone on its boundary the BFS is path-independent, so this is the same "
                "leak across G that the loop check reports, seen from the other side.");
        }
    }

    if (ffr.unreachedFaces > 0 || ffr.combingDefects > 0) {
        status.messages.push_back(
            "Stopping before the integration: the branch of the field over Omega is not single "
            "valued, so the seam transitions it would be constrained by are not defined.");
        return false;
    }

    // --- Stage 4F: the integration (Sec. 6) -------------------------------
    // The scaffold exists for one reason: the constraint rows need the arcs of
    // G, their (e+, e-) pairing and their quarter turns, and all three are
    // Immersion::buildArcs()' work. It is handed Omega's own coordinates as a
    // placeholder map, which it never uses for anything read here -- the arcs
    // and the pairing come off ConeCut, and the k off the frames.
    try {
        scaffold = std::make_unique<Immersion>(*cutter, *cones,
                                               cutter->getCutMesh().vertices,
                                               frames->fieldEdgeLengths(),
                                               frames->combedAngle());
    } catch (const std::exception &e) {
        status.messages.push_back(std::string("Stage 4F: could not build the seam pairing: ") + e.what());
        return false;
    }
    status.seamArcs = scaffold->getReport().arcs;
    status.frameKConflicts = scaffold->getReport().frameKConflicts;

    FieldIntegration::Options iopts;
    iopts.regularisation = options.integrationRegularisation;
    integration = std::make_unique<FieldIntegration>(*cutter, *frames, *scaffold, iopts);
    const FieldIntegration::Report &ir = integration->getReport();
    status.integrationRan = true;
    status.integrationSolved = ir.solved;
    status.integrationSeamResidual = ir.maxSeamResidual;
    status.integrationFlippedFaces = ir.flippedFaces;
    status.integrationFlippedAreaFraction =
        (ir.totalArea > 0.0) ? ir.flippedArea / ir.totalArea : 0.0;
    status.integrationMinAreaRatio = ir.minAreaRatio;
    status.integrationMaxFitResidual = ir.maxFitResidual;
    status.integrationMeanFitResidual = ir.meanFitResidual;
    status.integrationFlipsAtCones = ir.flipsAdjacentToCone;
    status.integrationNearestFlipToCone = ir.nearestFlipToCone;
    for (const std::string &m : ir.messages) status.messages.push_back("Stage 4F: " + m);

    if (!ir.solved) {
        status.messages.push_back("Stopping: the field could not be integrated on Omega.");
        return false;
    }
    integratedMap = integration->getUV();

    // --- Stage 4R: the untangling (Sec. 7.2) ------------------------------
    std::vector<Point> psi0 = integratedMap;
    if (ir.flippedFaces > 0 && options.untangle) {
        std::vector<Point> repaired = untangle();
        if (!repaired.empty()) psi0 = std::move(repaired);
    } else if (ir.flippedFaces == 0) {
        status.messages.push_back(
            "Stage 4F: the integration inverted nothing, so Sec. 7.2's untangling was not "
            "needed. That happens on gently curved, well-aligned models and is not to be "
            "assumed.");
    }

    // --- Stage 4: psi_0 ---------------------------------------------------
    try {
        immersion = std::make_unique<Immersion>(*cutter, *cones, psi0,
                                                frames->fieldEdgeLengths(),
                                                frames->combedAngle());
    } catch (const std::exception &e) {
        status.messages.push_back(std::string("Stage 4: could not accept psi_0: ") + e.what());
        return false;
    }
    const Immersion::Report &imr = immersion->getReport();
    status.immersionValid = imr.valid;
    status.immersionFlippedFaces = imr.flippedFaces;
    status.seamArcs = imr.arcs;
    status.maxSnapError = imr.maxSnapError;
    status.frameKConflicts = imr.frameKConflicts;
    status.maxMetricResidual = imr.maxMetricResidual;
    for (const std::string &m : imr.messages) status.messages.push_back("Stage 4: " + m);

    if (imr.flippedFaces > 0) {
        std::ostringstream oss;
        oss << "Stopping before the layout energies: psi_0 has " << imr.flippedFaces
            << " inverted face(s), so Q1 does not hold, and Sec. 3.3's barrier can only "
            << "preserve Q1 and never repair it. This is the cost of the substitution and "
            << "Sec. 7.3 is the list of what to try: raise the field's smoothness near the "
            << "cones, refine where the per-triangle fit residual is largest, lower h there.";
        status.messages.push_back(oss.str());
        return false;
    }
    if (!options.runLayout) return true;

    // --- Stage 5: the subdomain labelling, MERIDIAN's unchanged -----------
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
    for (const std::string &m : lr.messages) status.messages.push_back("Stage 5: " + m);

    // --- E1's reference metric, C4 ----------------------------------------
    //
    // The lengths psi_0 induces on Omega, re-indexed by the edges of S, which is
    // the indexing buildReference() reads. A seam edge has two images in Omega
    // and they are the same length to the seam residual, because the transition
    // holding them together is a rotation; either one will do and the first is
    // taken.
    std::vector<double> inducedLen;
    if (options.reference == Options::Reference::Induced) {
        inducedLen = inducedLengths(*mesh, *cutter, immersion->getUV());
    } else if (options.reference == Options::Reference::Ricci) {
        // Pipeline A's flat cone metric, computed here for one purpose only:
        // to be E1's notion of undistorted. psi_0 still comes from the
        // integration. This is the known-good yardstick the other two are
        // measured against, and it is the one setting under which this pipeline
        // still needs the flow.
        referenceFlow = std::make_unique<RicciFlow>(mesh, *cones);
        referenceFlow->solve();
        for (const std::string &m : referenceFlow->getReport().messages) {
            status.messages.push_back("Stage 6 (reference flow): " + m);
        }
    }

    // What that reference actually does at the cones, measured before the
    // continuation is allowed to use it. C4 in one number: a reference with no
    // cones reports (pi/2)|I| here, and the plan's whole warning is that E1 and
    // Q2 are then contradictory statements about the same vertex.
    {
        const std::vector<double> *refLen = nullptr;
        std::vector<double> ricciLen;
        if (options.reference == Options::Reference::Induced) {
            refLen = &inducedLen;
        } else if (options.reference == Options::Reference::Ricci && referenceFlow) {
            ricciLen = referenceFlow->originalEdgeLengthsCompleted(nullptr);
            refLen = &ricciLen;
        }
        std::vector<double> euclid;
        if (!refLen) {
            // Field and Euclidean are the same reference; measure the one they
            // both are.
            euclid.assign(mesh->edges.size(), 0.0);
            for (size_t e = 0; e < euclid.size(); ++e) {
                euclid[e] = normP(mesh->vertices[mesh->edges[e][1]] -
                                  mesh->vertices[mesh->edges[e][0]]);
            }
            refLen = &euclid;
        }
        std::vector<double> sum(mesh->vertices.size(), 0.0);
        for (int t = 0; t < static_cast<int>(mesh->triangles.size()); ++t) {
            const double a = (*refLen)[mesh->triangleEdges[t][0]];
            const double b = (*refLen)[mesh->triangleEdges[t][1]];
            const double c = (*refLen)[mesh->triangleEdges[t][2]];
            const double l[3] = {a, b, c};   // l[i] joins tri[i] and tri[i+1]
            for (int i = 0; i < 3; ++i) {
                const double p = l[i], q = l[(i + 2) % 3], r = l[(i + 1) % 3];
                if (!(p > 0.0) || !(q > 0.0)) continue;
                double cosA = (p * p + q * q - r * r) / (2.0 * p * q);
                cosA = std::max(-1.0, std::min(1.0, cosA));
                sum[mesh->triangles[t][i]] += std::acos(cosA);
            }
        }
        const std::vector<int> &idx = cones->getIndices();
        for (const auto &c : cones->getCones()) {
            const double full = mesh->isBoundaryVertex[c.vertex] ? M_PI : 2.0 * M_PI;
            const double want = full - M_PI_2 * idx[c.vertex];
            status.referenceConeResidual =
                std::max(status.referenceConeResidual, std::fabs(sum[c.vertex] - want));
        }
    }

    // --- Stage 6: the layout energies -------------------------------------
    LayoutEnergy::Options eopts;
    eopts.lambdaInit = options.lambdaInit;
    eopts.lambdaGrowth = options.lambdaGrowth;
    eopts.outerSteps = options.outerSteps;
    eopts.innerIterations = options.innerIterations;
    eopts.alternateReference = false;
    eopts.relabel = options.relabelBetweenSteps;
    eopts.lagInterfaceScales = options.lagInterfaceScales;
    switch (options.reference) {
        case Options::Reference::Induced:
            eopts.reference = LayoutEnergy::Reference::Induced;
            eopts.referenceLengths = inducedLen;
            break;
        case Options::Reference::Field:
            eopts.reference = LayoutEnergy::Reference::Field;
            eopts.fieldFrames = frames->frames();
            break;
        case Options::Reference::Euclidean:
            eopts.reference = LayoutEnergy::Reference::Euclidean;
            break;
        case Options::Reference::Ricci:
            eopts.reference = LayoutEnergy::Reference::Induced;
            eopts.referenceLengths = referenceFlow
                ? referenceFlow->originalEdgeLengthsCompleted(nullptr)
                : std::vector<double>();
            break;
    }
    // Sec. 9's closing note. E2 and E3 start an order of magnitude lower
    // because the field is already aligned to dS and to the interfaces, so
    // their residuals start small and over-penalising them early only fights
    // E1; E4 starts higher because the untangling introduced seam error, and
    // because -- exactly as in Pipeline A -- Q4 and Q2 are the same statement
    // in different units and Q4 is not free to trade away.
    eopts.lambdaFactor[2] = options.lambdaAlignmentFactor;
    eopts.lambdaFactor[3] = options.lambdaAlignmentFactor;
    eopts.lambdaFactor[4] = options.lambdaSeamFactor;
    layout = std::make_unique<LayoutEnergy>(*immersion, *labels, eopts);
    status.layoutRan = layout->run();
    {
        const LayoutEnergy::Report &er = layout->getReport();
        status.layoutInjective = er.injective;
        status.layoutConstrained = er.constraintsMet;
        status.outerStepsTaken = er.outerSteps;
        status.layoutValid = er.valid;
        status.layoutConeAngleResidual = er.maxConeAngleResidual;
        status.layoutRegularAngleResidual = er.maxRegularAngleResidual;
        status.layoutReferenceLeftHanded =
            std::max(status.layoutReferenceLeftHanded, er.leftHandedFrames);
        status.interfaceResidual = er.maxInterfaceResidual;
        status.interfaceCornerChanges = er.interfaceCornerChanges;
        status.interfacesAligned = er.interfaceCornerChanges == 0 &&
                                   er.maxInterfaceResidual < 1e-3;
        for (const std::string &m : er.messages) status.messages.push_back("Stage 6: " + m);
    }

    if (!options.runSeparatrices) return status.layoutValid;

    // --- Stage 7 and Sec. 3.3's repair, MERIDIAN's unchanged --------------
    Separatrices::Options sopts;
    sopts.coneSnapTolerance = options.separatrixSnap;
    sopts.maxSteps = options.separatrixMaxSteps;
    sopts.coneSnapRings = options.separatrixSnapRings;
    sopts.detectCycles = options.separatrixDetectCycles;
    sopts.nearMissWindow = options.repairGapLimit;
    if (interfaces && interfaces->multiMaterial()) {
        sopts.extraEmitters = interfaces->emitterNodes();
    }

    MERIDIAN::RepairOptions ropts;
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

    MERIDIAN::RepairResult rep = MERIDIAN::traceAndRepair(
        *labels, *layout, sopts, ropts,
        [&](const std::string &m) { status.messages.push_back(m); });
    separatrices = std::move(rep.separatrices);
    status.repairPasses = rep.passesTaken;
    status.repairConstraintsAdded = rep.constraintsAdded;
    status.topoPaths = labels->getReport().topoPaths;
    status.topoSelfReturns = labels->getReport().topoSelfReturns;
    status.topoExtraPerPair = labels->getReport().topoExtraPerPair;

    {
        const LayoutEnergy::Report &er = layout->getReport();
        status.layoutInjective = er.injective;
        status.layoutConstrained = er.constraintsMet;
        status.outerStepsTaken = er.outerSteps;
        status.layoutValid = er.valid;
        status.layoutConeAngleResidual = er.maxConeAngleResidual;
        status.layoutRegularAngleResidual = er.maxRegularAngleResidual;
    }

    if (separatrices) {
        const Separatrices::Report &sr = separatrices->getReport();
        status.separatricesRan = true;
        status.separatrices = sr.emitted;
        status.separatricesToCone = sr.endedAtCone;
        status.separatricesToBoundary = sr.endedAtBoundary;
        status.separatricesUnresolved = sr.capped + sr.cycled + sr.stuck + sr.degenerate;
        status.separatricesNearMisses = sr.nearMisses + sr.grazes;
        status.q5Verified = sr.valid;
    }

    // --- Stage 8 ----------------------------------------------------------
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

    // --- Stage 9 ----------------------------------------------------------
    if (!options.runSplines) return status.layoutValid;
    SplineFit::Options sfopts;
    sfopts.segments = options.splineSegments;
    sfopts.samples = options.splineSamples;
    sfopts.fitBoundaryArcs = options.fitBoundaryArcs;
    sfopts.fitInterfaceArcs = options.fitInterfaceArcs;
    try {
        splines = std::make_unique<SplineFit>(*arrangement, sfopts);
    } catch (const std::exception &e) {
        status.messages.push_back(std::string("Stage 9: could not be built: ") + e.what());
        return status.layoutValid;
    }
    const SplineFit::Report &spr = splines->getReport();
    status.splinesRan = true;
    status.splinePatches = spr.patches;
    status.splineControlPoints = spr.controlPointsPerArc;
    status.splineMaxDeviation = spr.maxDeviation;
    status.splinesWatertight = spr.watertight;
    status.splinesValid = spr.valid;
    for (const std::string &m : spr.messages) status.messages.push_back("Stage 9: " + m);

    // --- Stage 10 ---------------------------------------------------------
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
