#include "ATLAS/ATLAS.hxx"

#include <algorithm>
#include <chrono>
#include <sstream>

namespace {

typedef std::chrono::steady_clock Clock;

double since(Clock::time_point t) {
    return std::chrono::duration<double>(Clock::now() - t).count();
}

} // namespace

ATLAS::ATLAS(std::shared_ptr<Mesh> mesh, const Options &opts) : mesh_(std::move(mesh)), opts_(opts) {}

bool ATLAS::run() {
    status_ = Status();
    status_.triangles = static_cast<int>(mesh_->triangles.size());

    // ---- Stage 1 ------------------------------------------------------------
    Clock::time_point t = Clock::now();
    domain_ = std::make_unique<PlanarDomain>(*mesh_, opts_.domain);
    status_.secondsStage1 = since(t);
    status_.domainValid = domain_->getReport().valid;
    for (const std::string &m : domain_->getReport().messages) status_.messages.push_back(m);
    if (!status_.domainValid) {
        status_.messages.push_back("Stopping after Stage 1: the domain does not satisfy Sec. 1.1");
        return false;
    }

    // ---- Stage 2 ------------------------------------------------------------
    t = Clock::now();
    work_ = std::make_unique<SquareCarrier>(*domain_, opts_.carrier);
    initialReport_ = work_->validate();
    initialCarrier_ = std::make_unique<SquareCarrier>(*work_);
    status_.secondsStage2 = since(t);
    status_.carrierValid = initialReport_.valid;
    status_.initialCells = work_->numCells();
    for (const std::string &m : initialReport_.messages) status_.messages.push_back(m);
    if (!status_.carrierValid) {
        status_.messages.push_back("Stopping after Stage 2: the carrier does not validate");
        return false;
    }

    // ---- Stage 3 ------------------------------------------------------------
    t = Clock::now();
    if (opts_.runTemplates) {
        templates_ = std::make_unique<ExplicitTemplates>(*work_, opts_.templates);
        const SquareCarrier::Report &r = work_->validate();
        status_.templatesValid = r.valid;
        if (!r.valid) {
            for (const std::string &m : r.messages) status_.messages.push_back("after Stage 3: " + m);
            status_.messages.push_back("Stopping after Stage 3: a committed template broke the carrier");
            return false;
        }
    }
    templatedCarrier_ = std::make_unique<SquareCarrier>(*work_);
    status_.cellsAfterTemplates = work_->numCells();
    status_.secondsStage3 = since(t);

    // ---- Stages 4 -> 5 -> 6 ---------------------------------------------------
    const Clock::time_point loopStart = Clock::now();
    if (opts_.runRewrite) rewrite_ = std::make_unique<CavityRewrite>(*work_, opts_.rewrite);
    for (int r = 0;; ++r) {
        t = Clock::now();
        auto snap = std::make_shared<SquareCarrier>(*work_);
        auto rects = std::make_shared<RectangleCertifier>(*snap, opts_.rectangles);
        auto cover = std::make_shared<BlockCover>(*rects, opts_.cover);

        Round rd;
        rd.round = r;
        rd.cells = snap->numCells();
        rd.defect = snap->getReport().totalDefect;
        rd.irregular = snap->getReport().irregularInterior + snap->getReport().irregularBoundary;
        rd.candidates = rects->getReport().candidates;
        rd.blocks = cover->getReport().blocks;
        rd.macroVertices = cover->getReport().macroVertices;
        rd.objective = cover->getReport().objective;
        rd.valid = cover->getReport().valid;

        // Stage 4's failure witnesses, as vertices of the working carrier
        // (the snapshot is an exact copy of it).
        std::vector<int> priority;
        for (const RectangleCertifier::Witness &w : rects->getWitnesses()) {
            if (w.vertex >= 0) priority.push_back(w.vertex);
            for (int q : {w.cellA, w.cellB}) {
                if (q < 0) continue;
                for (int v : snap->cells[q]) if (snap->defect(v) > 0) priority.push_back(v);
            }
        }

        if (r == 0) {
            firstCarrier_ = snap;
            firstRects_ = rects;
            firstCover_ = cover;
        }
        const bool better = rd.valid && (!bestCover_ || !bestCover_->getReport().valid ||
                                         rd.objective < status_.bestObjective - 1e-9);
        if (better || !bestCover_) {
            bestCarrier_ = snap;
            bestRects_ = rects;
            bestCover_ = cover;
            status_.bestRound = r;
            status_.bestBlocks = rd.blocks;
            status_.bestObjective = rd.objective;
        }
        rd.seconds = since(t);
        status_.history.push_back(rd);

        if (!rewrite_ || r >= opts_.rewriteRounds) break;
        if (since(loopStart) > opts_.timeBudget) {
            status_.messages.push_back("Stages 4-6: the time budget ran out; the best cover so far is kept");
            break;
        }
        // ---- Stage 6 ----
        t = Clock::now();
        int committed = 0;
        bool broken = false;
        for (int pass = 0; pass < std::max(1, opts_.rewritePassesPerRound); ++pass) {
            // Stage 4's witnesses name vertices of this carrier only until the
            // first pass rewrites it.
            const int c = rewrite_->round(pass == 0 ? priority : std::vector<int>());
            committed += c;
            const SquareCarrier::Report &rep = work_->validate();
            if (!rep.valid) {
                status_.rewritesValid = false;
                for (const std::string &m : rep.messages) status_.messages.push_back("after Stage 6: " + m);
                status_.messages.push_back("Stopping Stage 6: a committed rewrite broke the carrier");
                broken = true;
                break;
            }
            if (c == 0 || since(loopStart) > opts_.timeBudget) break;
        }
        status_.history.back().committed = committed;
        status_.history.back().rings = rewrite_->maxRings();
        status_.history.back().seconds += since(t);
        if (broken) break;
        // A stalled round widens the cavities; a round that commits nothing
        // at the widest radius ends the loop.
        const bool stalled = committed < std::max(3, rd.irregular / 100);
        if (stalled && rewrite_->maxRings() < opts_.maxRings) {
            rewrite_->setMaxRings(rewrite_->maxRings() + 1);
        } else if (committed == 0) {
            break;
        }
    }
    status_.secondsLoop = since(loopStart);

    status_.coverValid = bestCover_ && bestCover_->getReport().valid;
    status_.finalCells = bestCarrier_ ? bestCarrier_->numCells() : 0;
    for (const std::string &m : bestCover_->getReport().messages) status_.messages.push_back(m);
    return status_.domainValid && status_.carrierValid && status_.templatesValid &&
           status_.rewritesValid && status_.coverValid;
}
