#include "MERIDIAN.hxx"

#include <sstream>
#include <stdexcept>

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

    // --- Stage 0: cross field ---------------------------------------------
    runField();

    // --- Stage 1: cone singularities, Sec. 3.1 ----------------------------
    cones = std::make_unique<ConeSingularities>(*field);
    cones->setBoundaryIndexRange(options.minBoundaryIndex, options.maxBoundaryIndex);

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

    // --- Stage 2: cutting graph, Sec. 3.2.2 -------------------------------
    cutter = std::make_unique<ConeCut>(mesh, *cones);
    const ConeCut::Report &cr = cutter->getReport();
    status.cutIsDisk = cr.isDisk;
    status.allConesOnBoundary = cr.allConesOnBoundary;
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
    labels = std::make_unique<SubdomainLabels>(*immersion, lopts);
    const SubdomainLabels::Report &lr = labels->getReport();
    status.boundaryEdgesU = lr.boundaryEdgesU;
    status.boundaryEdgesV = lr.boundaryEdgesV;
    status.featureChains = lr.featureChains;
    status.topoPaths = lr.topoPaths;
    for (const std::string &m : lr.messages) status.messages.push_back("Stage 5: " + m);

    // --- Stage 6: the layout-inducing energies, Sec. 3.3 ------------------
    LayoutEnergy::Options eopts;
    eopts.lambdaInit = options.lambdaInit;
    eopts.lambdaGrowth = options.lambdaGrowth;
    eopts.outerSteps = options.outerSteps;
    eopts.innerIterations = options.innerIterations;
    eopts.alternateReference = options.alternateReference;
    eopts.relabel = options.relabelBetweenSteps;
    layout = std::make_unique<LayoutEnergy>(*immersion, *labels, eopts);
    status.layoutRan = layout->run();
    const LayoutEnergy::Report &er = layout->getReport();
    status.layoutInjective = er.injective;
    status.layoutConstrained = er.constraintsMet;
    status.outerStepsTaken = er.outerSteps;
    status.layoutValid = er.valid;
    for (const std::string &m : er.messages) status.messages.push_back("Stage 6: " + m);

    if (!options.runSeparatrices) return status.layoutValid;

    // --- Stage 7: separatrix tracing, Sec. 4 ------------------------------
    // Run even when the continuation fell short. The curves are the only place
    // Q5 is visible as the property it actually is -- every integral curve out
    // of a cone is finite -- and a capped one names the pair of cones whose
    // connectivity constraint is missing, which is the diagnosis Sec. 3.3 asks
    // for when the layout is not yet a layout.
    Separatrices::Options sopts;
    sopts.coneSnapTolerance = options.separatrixSnap;
    sopts.maxSteps = options.separatrixMaxSteps;
    try {
        separatrices = std::make_unique<Separatrices>(*immersion, layout->getUV(), sopts);
    } catch (const std::exception &e) {
        status.messages.push_back(std::string("Stage 7: could not be run: ") + e.what());
        return status.layoutValid;
    }
    const Separatrices::Report &sr2 = separatrices->getReport();
    status.separatricesRan = true;
    status.separatrices = sr2.emitted;
    status.separatricesToCone = sr2.endedAtCone;
    status.separatricesToBoundary = sr2.endedAtBoundary;
    status.separatricesUnresolved = sr2.capped + sr2.stuck + sr2.degenerate;
    status.q5Verified = sr2.valid;
    for (const std::string &m : sr2.messages) status.messages.push_back("Stage 7: " + m);

    return status.layoutValid;
}
