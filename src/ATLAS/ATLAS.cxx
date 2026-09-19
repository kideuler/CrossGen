#include "ATLAS/ATLAS.hxx"

#include <algorithm>
#include <chrono>
#include <sstream>
#include <thread>

namespace {

typedef std::chrono::steady_clock Clock;

double since(Clock::time_point t) {
    return std::chrono::duration<double>(Clock::now() - t).count();
}

} // namespace

ATLAS::ATLAS(std::shared_ptr<Mesh> mesh, const Options &opts) : mesh_(std::move(mesh)), opts_(opts) {}

bool ATLAS::run() {
    status_ = Status();
    searches_.clear();
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

    // ---- Stage 2, and the fine search: the guaranteed incumbent ---------------
    t = Clock::now();
    auto fine = std::make_unique<Search>();
    fine->name = "fine";
    fine->work = std::make_unique<SquareCarrier>(*domain_, opts_.carrier);
    fine->initialReport = fine->work->validate();
    fine->initial = std::make_unique<SquareCarrier>(*fine->work);
    fine->carrierValid = fine->initialReport.valid;
    status_.secondsStage2 = since(t);
    status_.carrierValid = fine->carrierValid;
    status_.initialCells = fine->work->numCells();
    for (const std::string &m : fine->initialReport.messages) status_.messages.push_back(m);
    searches_.push_back(std::move(fine));
    if (!status_.carrierValid) {
        status_.messages.push_back("Stopping after Stage 2: the carrier does not validate");
        return false;
    }

    const Clock::time_point searchStart = Clock::now();
    // The coarse searches first: when one succeeds the fine one need not run
    // Stage 6 at all.
    const bool coarseApplies = opts_.coarse && domain_->getReport().interfaceEdges == 0;
    if (opts_.coarse && !coarseApplies) {
        status_.messages.push_back("Coarse searches skipped: CoarseDomain does not yet carry material "
                                   "interfaces, so a multi-material domain is searched on its fine carrier only");
    }
    // The coarse domains are built here, one per spacing and in turn (Triangle
    // is not re-entrant); every search after that only reads its domain, so
    // the searches -- one per spacing and seed, and the fine one -- run in
    // parallel threads, each deterministic on its own.
    if (coarseApplies) {
        for (double h : opts_.coarseSpacings) {
            CoarseDomain::Options co = opts_.coarseDomain;
            co.maxSpacing = h;
            auto domain = std::make_shared<CoarseDomain>(*domain_, co);
            for (int k = 0; k < std::max(1, opts_.coarseSeeds); ++k) {
                auto sp = std::make_unique<Search>();
                sp->coarse = true;
                sp->spacing = h;
                sp->seed = opts_.anneal.seed + 7919u * static_cast<unsigned>(k);
                sp->coarseDomain = domain;
                std::ostringstream os;
                os << "coarse h=" << h;
                if (opts_.coarseSeeds > 1) os << " seed " << k + 1;
                sp->name = os.str();
                searches_.push_back(std::move(sp));
            }
        }
    }
    // Stage 6 on the fine carrier only where there is no coarse search to
    // beat it (or on request): it is the search that stalls (CoarseDomain).
    const bool fineRewrite = opts_.runRewrite && (opts_.fineRewrite || !coarseApplies);
    auto runFine = [&]() {
        Search &f = *searches_[0];
        const Clock::time_point tf = Clock::now();
        search(f, fineRewrite, opts_.rewriteRounds, opts_.rewritePassesPerRound, opts_.maxRings, false);
        f.seconds = since(tf);
    };
    if (opts_.parallel && searches_.size() > 1) {
        std::vector<std::thread> pool;
        pool.emplace_back(runFine);
        for (size_t i = 1; i < searches_.size(); ++i) pool.emplace_back([this, i]() { runCoarse(*searches_[i]); });
        for (std::thread &t : pool) t.join();
    } else {
        for (size_t i = 1; i < searches_.size(); ++i) runCoarse(*searches_[i]);
        runFine();
    }
    status_.secondsSearch = since(searchStart);

    // ---- Pick the answer ------------------------------------------------------
    for (int i = 0; i < numSearches(); ++i) {
        const Search &s = *searches_[i];
        if (!s.succeeded()) continue;
        const double obj = s.finalCover()->getReport().objective;
        if (status_.chosen < 0 || obj < status_.bestObjective - 1e-9) {
            status_.chosen = i;
            status_.bestObjective = obj;
            status_.bestBlocks = s.finalCover()->getReport().blocks;
        }
    }
    if (status_.chosen < 0) {
        // Nothing passed Sec. 13.2: report the fine search's cover anyway.
        if (searches_[0]->bestCover) {
            status_.chosen = 0;
            status_.bestObjective = searches_[0]->bestCover->getReport().objective;
            status_.bestBlocks = searches_[0]->bestCover->getReport().blocks;
        } else {
            return false;
        }
    }
    const Search &c = getChosen();
    status_.coverValid = c.succeeded();
    status_.finalCells = c.finalCarrier()->numCells();
    for (const std::string &m : c.finalCover()->getReport().messages) status_.messages.push_back(m);
    return status_.domainValid && status_.carrierValid && c.templatesValid && c.rewritesValid &&
           status_.coverValid;
}

// ---------------------------------------------------------------------------
// One coarse search: the coarse domain, its carrier, Stages 3-6 there, and the
// realisation of the best round on the input domain.
// ---------------------------------------------------------------------------
void ATLAS::runCoarse(Search &s) {
    const Clock::time_point t0 = Clock::now();
    if (!s.coarseDomain->valid()) {
        s.messages.push_back("CoarseDomain: " + s.coarseDomain->getReport().reason);
        s.seconds = since(t0);
        return;
    }
    s.work = std::make_unique<SquareCarrier>(s.coarseDomain->getDomain(), opts_.carrier);
    s.initialReport = s.work->validate();
    s.initial = std::make_unique<SquareCarrier>(*s.work);
    s.carrierValid = s.initialReport.valid;
    if (!s.carrierValid) {
        for (const std::string &m : s.initialReport.messages) s.messages.push_back("coarse " + m);
        s.seconds = since(t0);
        return;
    }
    search(s, opts_.runRewrite, opts_.coarseRewriteRounds, opts_.coarsePassesPerRound, opts_.coarseMaxRings,
           opts_.anneal.maxMoves > 0);
    if (!s.best) {
        s.seconds = since(t0);
        return;
    }

    // ---- Realise on the input domain: the best round first -------------------
    // A layout can be valid on the coarse proxy and still not realise -- a
    // block that is thin where the true curve bulges into it folds -- so the
    // other rounds' layouts are tried in order of objective until one passes
    // validate() on the input. The report keeps the first failure's reasons.
    const Clock::time_point tr = Clock::now();
    std::vector<Search::Candidate> order = s.candidates;
    std::stable_sort(order.begin(), order.end(),
                     [](const Search::Candidate &a, const Search::Candidate &b) { return a.objective < b.objective; });
    if (order.empty()) order.push_back({0.0, s.bestRound, s.best, s.bestRects, s.bestCover});
    for (const Search::Candidate &cand : order) {
        if (s.realisationAttempts >= std::max(1, opts_.maxRealisations)) break;
        ++s.realisationAttempts;
        auto real = std::make_unique<Realisation>(*s.coarseDomain, *cand.carrier,
                                                  cand.cover && cand.cover->getReport().valid ? cand.cover.get()
                                                                                              : nullptr,
                                                  opts_.realisation);
        if (!real->built()) {
            if (s.realisationAttempts == 1) s.messages.push_back("Realisation: " + real->getReport().reason);
            if (!s.realisation) s.realisation = std::move(real);
            continue;
        }
        auto fineCarrier = std::make_shared<SquareCarrier>(*domain_, opts_.carrier, real->getQuads());
        const SquareCarrier::Report rep = fineCarrier->validate();
        if (!rep.valid) {
            if (s.realisationAttempts == 1) {
                for (const std::string &m : rep.messages) s.messages.push_back("realised " + m);
                s.realisation = std::move(real);
                s.realised = fineCarrier;
                s.realisedReport = rep;
            }
            continue;
        }
        s.realisation = std::move(real);
        s.realised = fineCarrier;
        s.realisedReport = rep;
        s.realisedRound = cand.round;
        s.realisedRects = std::make_shared<RectangleCertifier>(*s.realised, opts_.rectangles);
        s.realisedCover = std::make_shared<BlockCover>(*s.realisedRects, opts_.cover);
        if (s.realisationAttempts > 1) {
            s.messages.push_back("Realisation: the best layout did not realise; round " + std::to_string(cand.round) +
                                 "'s did, on attempt " + std::to_string(s.realisationAttempts));
        }
        break;
    }
    s.secondsRealisation = since(tr);
    s.seconds = since(t0);
}

// ---------------------------------------------------------------------------
// Stages 3-6 on one carrier.
// ---------------------------------------------------------------------------
void ATLAS::search(Search &s, bool rewrite, int rounds, int passes, int maxRings, bool anneal) {
    SquareCarrier &work = *s.work;

    // ---- Stage 3 ------------------------------------------------------------
    if (opts_.runTemplates) {
        s.templates = std::make_unique<ExplicitTemplates>(work, opts_.templates);
        const SquareCarrier::Report &r = work.validate();
        s.templatesValid = r.valid;
        if (!r.valid) {
            for (const std::string &m : r.messages) s.messages.push_back("after Stage 3: " + m);
            s.messages.push_back("Stopping after Stage 3: a committed template broke the carrier");
            return;
        }
    }
    s.templated = std::make_unique<SquareCarrier>(work);

    // ---- Stages 4 -> 5 -> 6 ---------------------------------------------------
    const Clock::time_point loopStart = Clock::now();
    CavityRewrite::Options ro = opts_.rewrite;
    if (rewrite) s.rewrite = std::make_unique<CavityRewrite>(work, ro);
    double bestObjective = 0.0;
    // Stages 4 and 5 on a snapshot of the working carrier, recorded as round
    // r; returns Stage 4's witnesses as vertices of the working carrier.
    auto evaluate = [&](int r) {
        Clock::time_point t = Clock::now();
        auto snap = std::make_shared<SquareCarrier>(work);
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

        std::vector<int> priority;
        for (const RectangleCertifier::Witness &w : rects->getWitnesses()) {
            if (w.vertex >= 0) priority.push_back(w.vertex);
            for (int q : {w.cellA, w.cellB}) {
                if (q < 0) continue;
                for (int v : snap->cells[q]) if (snap->defect(v) > 0) priority.push_back(v);
            }
        }
        if (r == 0) {
            s.first = snap;
            s.firstRects = rects;
            s.firstCover = cover;
        }
        if (s.coarse && rd.valid) s.candidates.push_back({rd.objective, r, snap, rects, cover});
        const bool better = rd.valid && (!s.bestCover || !s.bestCover->getReport().valid ||
                                         rd.objective < bestObjective - 1e-9);
        if (better || !s.bestCover) {
            s.best = snap;
            s.bestRects = rects;
            s.bestCover = cover;
            s.bestRound = r;
            bestObjective = rd.objective;
        }
        rd.seconds = since(t);
        s.history.push_back(rd);
        return priority;
    };

    int r = 0;
    for (;; ++r) {
        const std::vector<int> priority = evaluate(r);
        const Round &rd = s.history.back();
        if (!s.rewrite || r >= rounds) break;
        if (since(loopStart) > opts_.timeBudget) {
            s.messages.push_back("Stages 4-6: the time budget ran out; the best cover so far is kept");
            break;
        }
        // ---- Stage 6 ----
        Clock::time_point t = Clock::now();
        int committed = 0;
        bool broken = false;
        for (int pass = 0; pass < std::max(1, passes); ++pass) {
            // Stage 4's witnesses name vertices of this carrier only until the
            // first pass rewrites it.
            const int c = s.rewrite->round(pass == 0 ? priority : std::vector<int>());
            committed += c;
            const SquareCarrier::Report &rep = work.validate();
            if (!rep.valid) {
                s.rewritesValid = false;
                for (const std::string &m : rep.messages) s.messages.push_back("after Stage 6: " + m);
                s.messages.push_back("Stopping Stage 6: a committed rewrite broke the carrier");
                broken = true;
                break;
            }
            if (c == 0 || since(loopStart) > opts_.timeBudget) break;
        }
        s.history.back().committed = committed;
        s.history.back().rings = s.rewrite->maxRings();
        s.history.back().seconds += since(t);
        if (broken) return;
        // A stalled round widens the cavities; a round that commits nothing
        // at the widest radius ends the loop.
        const bool stalled = committed < std::max(3, rd.irregular / 100);
        if (stalled && s.rewrite->maxRings() < maxRings) {
            s.rewrite->setMaxRings(s.rewrite->maxRings() + 1);
        } else if (committed == 0) {
            break;
        }
    }

    // ---- The annealed search, from the greedy rounds' end state -----------------
    if (anneal && s.rewrite && s.rewritesValid) {
        // Start from the best round, not the last: the greedy rounds only
        // stop when they stall, and the best cover may be an earlier one.
        if (s.best) work = *s.best;
        CavityRewrite::AnnealOptions ao = opts_.anneal;
        ao.seed = s.seed;
        if (s.coarse && s.coarseDomain && opts_.alignToField) {
            const CoarseDomain *cd = s.coarseDomain.get();
            ao.field = [cd](const Point &p, double *w) { return cd->crossAngle(p, w); };
        }
        s.rewrite->anneal(ao);
        const SquareCarrier::Report &rep = work.validate();
        if (!rep.valid) {
            s.rewritesValid = false;
            for (const std::string &m : rep.messages) s.messages.push_back("after annealing: " + m);
            return;
        }
        s.annealed = true;
        // The checkpoints as rounds of their own (realisation fallbacks),
        // then the final state; the working carrier is restored after.
        const SquareCarrier annealedBest = work;
        for (const SquareCarrier &cp : s.rewrite->getCheckpoints()) {
            work = cp;
            evaluate(++r);
        }
        work = annealedBest;
        evaluate(++r);
    }
    s.secondsLoop = since(loopStart);
}
