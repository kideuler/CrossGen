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
    return status.readyForImmersion;
}
