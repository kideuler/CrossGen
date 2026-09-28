#ifndef __ATLAS_HXX__
#define __ATLAS_HXX__

#include <atomic>
#include <memory>
#include <string>
#include <vector>

#include "ATLAS/BlockCover.hxx"
#include "ATLAS/CavityFill.hxx"
#include "ATLAS/CavityRewrite.hxx"
#include "ATLAS/CoarseDomain.hxx"
#include "ATLAS/ExplicitTemplates.hxx"
#include "ATLAS/PlanarDomain.hxx"
#include "ATLAS/Realisation.hxx"
#include "ATLAS/RectangleCertifier.hxx"
#include "ATLAS/ReferenceField.hxx"
#include "ATLAS/SquareCarrier.hxx"
#include "ATLAS/SquareTransport.hxx"
#include "mesh/Mesh.hxx"

// ATLAS: square-transport blocking of a planar triangle mesh, Stages 1 to 6 of
// docs/square_transport_2d_theory_and_implementation.md.
//
//   Stage 1  PlanarDomain         validate and tag the domain
//   Stage 2  SquareCarrier        the three-quad split and its transports
//   Stage 3  ExplicitTemplates    whole-region replacements (Sec. 7)
//   Stage 4  RectangleCertifier   integer development and certified rectangles
//   Stage 5  BlockCover           a conforming cover of certified blocks
//   Stage 6  CavityRewrite        rewrite the cavities round the obstructions
//
// Stages 4 to 6 repeat within a budget, and the best feasible decomposition
// seen is kept separately from the carrier Stage 6 goes on rewriting (Sec.
// 8.3: "Maintain a separate best feasible decomposition. Exploratory states
// need not be exported."). Each round works on a snapshot of the carrier, so
// the certificates and the cover of the best round stay attached to the exact
// carrier they were issued for.
//
// ### Two carriers
//
// Stages 3-6 run as a *search* on a carrier, and there are two kinds of
// carrier to search on:
//
//   * the fine one, Sec. 3's split of the input itself. It is the guaranteed
//     incumbent -- valid by construction, and its singleton blocking is a
//     valid answer -- and on a domain a Stage 3 template covers whole it is
//     already the best there is. But its topology is the input
//     triangulation's, and on anything else Stage 6 stalls on thousands of
//     isolated singularities (see CoarseDomain for the measurement and why);
//   * coarse ones, the same split of a re-triangulation of the domain at the
//     domain's own resolution (CoarseDomain), one per spacing in
//     Options::coarseSpacings and seed of the annealer (coarseSeeds). The
//     search there is cheap enough for Sec. 8.3's full energy -- greedy
//     rounds, then CavityRewrite::anneal -- and what it finds is carried back
//     onto the input (Realisation) as a new fine carrier, which has to pass
//     the same SquareCarrier::validate() the fine carrier did before Stages 4
//     and 5 are allowed to extract its blocks. If the best layout does not
//     realise, the search's other rounds are tried in order of objective.
//
// The searches share nothing but the read-only domains, so they run in
// parallel, one OpenMP thread each (Options::parallel), and each is
// deterministic.
//
// Every search that produces a valid, conforming cover of the *input* domain
// competes on Stage 5's objective, and the lowest wins; the fine incumbent
// always competes, so a failure of every coarse search is still a valid (fine)
// answer, which the report says it is.
//
// ### A reference cross field, as a prior
//
// With Options::field.reference = DualMBO, Stage 1b solves DualMBO on the
// input (ReferenceField, TORSION's Stage 0 settings) before the searches
// start, and every search reads it: Stage 3 places O-grid, half O-grid and
// star parameters on its cones and picks among families by score, Stage 6's
// annealer adds wDir E_dir + E_sing to its energy and aims some of its moves
// at where the carrier and the field disagree, and the arbiter above ranks
// rounds, realisations and searches by objective + wDir E_dir + E_sing
// (docs/atlas_crossfield_guidance.md). The field only ever scores states and
// sets template parameters; every state is still certified exactly as before,
// so Sec. 14's guarantees do not depend on it.
//
// ### Multi-material domains
//
// An interface cannot be re-subdivided by a cavity on one side of it (Sec.
// 8.4), so on a multi-material domain its sample count is fixed when
// CoarseDomain samples it, and Stage 3 and Stage 6 have to live with it. What
// makes that work (docs/atlas_multimaterial_blocks.md): CoarseDomain
// reconciles the counts before sampling, against the section arrangements
// Stage 3 will cast; Stage 3's templates treat a neighbour's split points as
// flat nodes to meet, not corners, may shed dS midpoints next to an
// interface, and cut regions with holes and awkward corners by general
// sections with star faces; Stage 6's cavities may straddle an interface and
// re-subdivide it on both sides at once, and the annealer may rewrite
// templates. Here, the second seed's searches template smallest region first,
// and the arbiter passes over a coarse search whose worst cell is under half
// the best one's and weighs the field's direction term double.
//
// The contrast with MERIDIAN and TORSION is the spec's central idea: nothing
// here is a cross field to be integrated and quantised afterwards. The carrier
// already is a set of valid charts with integer transitions (Sec. 2.3), every
// later state is one too, and extraction is a finite integer certificate
// (Sec. 5) rather than streamline tracing.
class ATLAS {
public:
    struct Options {
        PlanarDomain::Options domain;
        SquareCarrier::Options carrier;
        bool runTemplates = true;
        ExplicitTemplates::Options templates;
        RectangleCertifier::Options rectangles;
        BlockCover::Options cover;
        bool runRewrite = true;
        CavityRewrite::Options rewrite;
        // Rounds of Stages 4 -> 5 -> 6 after the first Stages 4 and 5.
        int rewriteRounds = 12;
        // Stage 6 passes between two evaluations of Stages 4 and 5. A pass
        // commits only vertex-disjoint cavities, so the next pass finds work
        // the last one had to leave, and it is far cheaper than recertifying.
        int rewritePassesPerRound = 4;
        // Widen Stage 6's cavities by a ring whenever a round commits fewer
        // rewrites than a hundredth of the irregular vertices, up to this many
        // rings (0 = never widen).
        int maxRings = 8;
        // Wall-clock budget for one search's Stage 4-6 loop, seconds.
        double timeBudget = 300.0;

        // ---- the fine search ----
        // Run Stage 6 on the fine carrier as well. Off by default when a
        // coarse search succeeds: it is the slow search that stalls (see
        // CoarseDomain), and the fine carrier's Stages 3-5 already give the
        // incumbent it would have to beat.
        bool fineRewrite = false;

        // ---- the coarse searches ----
        bool coarse = true;
        // CoarseDomain::Options::maxSpacing for each coarse search, as a
        // fraction of the bounding-box diagonal.
        std::vector<double> coarseSpacings{0.125, 0.07};
        // Independent annealed searches per spacing, each with its own seed;
        // the best realised one wins. They and the fine search run in
        // parallel threads.
        int coarseSeeds = 2;
        // On a domain with material interfaces, one more coarse search per
        // spacing with Stage 3 off, on the domain as sampled (without the
        // counts reconciled for templates). It was added when a certified
        // template could still wreck the regions round it -- multimat/
        // geom012's five-block star in one grain turned a 9-block answer into
        // 17 -- and Stage 6 alone did better. Flat nodes, prescribed spokes
        // and general sections now give the templated searches geom012's 9
        // themselves, and over the corpus this search then won only on
        // bubbles, by cutting its disks into D-shaped halves (346 blocks,
        // TMOP worst 0.33, against the templated 369 at 0.36). Off, so a
        // multi-material domain runs the same five searches as any other.
        bool coarseWithoutTemplates = false;
        // On a domain with material interfaces, the annealer may rewrite
        // template cells and give up designated vertices
        // (CavityRewrite::AnnealOptions::freeTemplates).
        bool freeTemplatesMultimat = true;
        // On a domain with material interfaces, the second seed's coarse
        // searches template their regions smallest first
        // (ExplicitTemplates::Options::smallestFirst): which order lets the
        // templates' split points line up depends on the model, and the
        // annealer seeds that differ on a single-material domain barely do
        // on a multi-material one.
        bool multimatOrders = true;
        // On a domain with material interfaces, Stage 6 may not rewrite the
        // O-grid Stage 3 put in a smooth inclusion (a region of no corners),
        // however free AnnealOptions::freeTemplates leaves the other templates
        // (CavityRewrite::lockGroups). Off by default: it is what lets the
        // annealer merge an inclusion's blocks into the layout round it, and
        // on multimat/bubbles that is 369 blocks (TMOP worst 0.36) against
        // 551 with every disk kept a five-block O-grid (0.51, the quality
        // ATLAS had before). On, for the second.
        bool lockInclusionOGrids = false;
        // Score the annealer's layouts against the coarse domain's harmonic
        // reference cross field (CoarseDomain::crossAngle) as well; only with
        // field.reference = Harmonic below.
        bool alignToField = true;
        bool parallel = true;
        CoarseDomain::Options coarseDomain;
        int coarseRewriteRounds = 30;
        int coarsePassesPerRound = 6;
        int coarseMaxRings = 6;
        // After the greedy rounds, the annealed search of Sec. 8.3 on each
        // coarse carrier (CavityRewrite::anneal); 0 moves turns it off.
        CavityRewrite::AnnealOptions anneal;
        Realisation::Options realisation;
        // Candidate coarse layouts tried, best objective first, until one
        // realises on the input; and, only when no coarse search realised in
        // that many, how far each goes on down its list before the fine
        // fallback is taken.
        int maxRealisations = 4;
        int maxRealisationsLastResort = 16;

        // ---- the reference cross field (docs/atlas_crossfield_guidance.md) ----
        struct Field {
            // DualMBO (the default since 2026-09-19): a ReferenceField solved
            // once after Stage 1 and read by every search, with its terms
            // wherever the switches below put them. Harmonic: the annealer's
            // term against CoarseDomain::crossAngle (alignToField above) and
            // nothing else, as before the field existed. None: no field.
            // Measured over the singlemat corpus, 3 seeds, together with the
            // sign-aware defects and patchTopology of CavityRewrite: blocks
            // 678 -> 690, TFI worst SJ median 0.21 -> 0.32, TMOP worst below
            // 0.3 on 5.3 -> 3.0 models, block corners on arcs 43 -> 13
            // (docs/atlas_crossfield_guidance.md, Sec. 11).
            enum class Reference { Harmonic, DualMBO, None };
            Reference reference = Reference::DualMBO;
            ReferenceField::Options solve;
            ReferenceField::SingularityWeights singularity;
            // Blocks per unit of E_dir: 10 makes a layout ~9 degrees off
            // everywhere cost one block.
            double wDir = 10.0;
            // Stage 3: cone-placed O-grid, half O-grid and star parameters,
            // and the family chosen by score (ExplicitTemplates::Options).
            bool templates = true;
            // Stage 6: the greedy rounds' choice among a witness's certified
            // fills, and cone-centred stars (CavityRewrite::Options::field),
            // on the coarse searches' carriers; greedyOnFine on the fine one
            // too. Off there by default: on the multimat fine carriers it
            // lowered E_dir by a fifth on every model (0.41 -> 0.33 on
            // bubbles) and raised the block count by 10-30%.
            bool greedy = true;
            bool greedyOnFine = false;
            // Stage 6: the annealer's energy (CavityRewrite::AnnealOptions).
            bool anneal = true;
            // The field-directed witnesses and cone-centred stars (R1).
            bool directedMoves = true;
            double directedFraction = 0.3;
            // Stage 5's arbiter: where ATLAS compares different carriers --
            // the best round of a search, the order a coarse search's
            // layouts are realised in, and the winning search -- the score is
            // the cover's objective + wDir E_dir + E_sing. Not in the cover's
            // own cost(P): inside one carrier every exact cover has the same
            // per-cell sum (guidance note, Sec. 4.3).
            bool arbiter = true;
        };
        Field field;
        // The arbiter's shape term, independent of the field: wShape (s0 -
        // SJ) / s0 for the carrier's worst cell below s0. The arbiter is
        // otherwise blind to geometry, which is how geom016's 3-block answer
        // with a 0.06 cell wins; the annealer has scored shape all along.
        // 0 = off.
        double arbiterShape = 0.0;
        double arbiterShapeFloor = 0.35;
        // Apply rewrite.boundaryPlusWeight / boundaryMinusWeight in the fine
        // search's greedy rounds too, not only the coarse searches'.
        bool signedDefectsOnFine = false;
        // On a domain with material interfaces, the arbiter passes over a
        // search whose realised carrier's worst cell is below this fraction
        // of the best worst cell any search reached (0 = off). The score is
        // blocks first and nearly blind to geometry (Stage 5's distortion
        // weighs 0.05 a block), and on multi-material domains the searches
        // differ in shape far more than on single-material ones: the search
        // without templates found geom002 in 7 blocks round a 0.04 cell, where
        // a templated one had 12 round a 0.38, and bubbles in 324 blocks whose
        // mesh folded, where one of 346 did not. Measured on the multimat
        // corpus, a worst cell under half the best is always the one whose
        // mesh is visibly worse, and every such search had a sound one within
        // a few blocks of it.
        double multimatQualityGate = 0.5;
        // The arbiter's weight on E_dir on a domain with material interfaces,
        // in place of field.wDir there (0 = field.wDir). Stage 3 now offers a
        // multi-material domain several templated layouts of close block
        // counts that differ in how they sit in the field -- multimat/geom001
        // in 5 blocks with a diagonal splitting a corner of the box (E_dir
        // 0.19), or in TORSION's 8 (0.02) -- and at 10 blocks per unit the
        // arbiter took the 5, whose mesh was worse (TMOP worst 0.67 against
        // 0.83). Replayed over the corpus, 15 or 20 changes that model's
        // choice and no other; 30 starts to trade a worse worst cell on
        // rocket for alignment.
        double multimatArbiterWDir = 20.0;
    };

    struct Round {
        int round = 0;
        int cells = 0;
        int defect = 0;
        int irregular = 0;
        int candidates = 0;
        int blocks = 0;
        int macroVertices = 0;
        double objective = 0.0;
        // What the arbiter ranks rounds by: the objective, plus the field's
        // terms and the shape term when they are on (Options::field.arbiter,
        // arbiterShape); equal to the objective otherwise.
        double score = 0.0;
        double eDir = 0.0, eSing = 0.0, shape = 0.0;
        bool valid = false;
        int committed = 0;      // Stage 6 rewrites after this round's cover
        int rings = 0;          // Stage 6's cavity radius in this round
        double seconds = 0.0;
    };

    // Stages 3-6 on one carrier, and for a coarse one its realisation.
    struct Search {
        std::string name;
        bool coarse = false;
        double spacing = 0.0;                 // coarse: maxSpacing used
        unsigned seed = 0;                    // coarse: the annealer's seed
        bool useTemplates = true;             // Stage 3 runs (Options::runTemplates too)
        bool smallestFirst = false;           // Stage 3's region order
        std::shared_ptr<CoarseDomain> coarseDomain;

        // The searched carrier (fine, or over coarseDomain->getDomain()).
        bool carrierValid = false;
        bool templatesValid = true;
        bool rewritesValid = true;
        std::unique_ptr<SquareCarrier> work;
        std::unique_ptr<SquareCarrier> initial;
        SquareCarrier::Report initialReport;
        std::unique_ptr<SquareCarrier> templated;
        std::unique_ptr<ExplicitTemplates> templates;
        std::unique_ptr<CavityRewrite> rewrite;
        std::vector<Round> history;
        int bestRound = -1;
        bool annealed = false;              // the last round is the annealed carrier
        // A coarse search's valid covers, one per round (the annealer's
        // checkpoints included): what Realisation falls back on, best
        // objective first, when the best one does not realise.
        struct Candidate {
            double objective = 0.0;
            double score = 0.0;
            int round = -1;
            std::shared_ptr<SquareCarrier> carrier;
            std::shared_ptr<RectangleCertifier> rects;
            std::shared_ptr<BlockCover> cover;
        };
        std::vector<Candidate> candidates;
        int realisedRound = -1;             // the round whose carrier was realised
        int realisationAttempts = 0;
        std::shared_ptr<SquareCarrier> best, first;
        std::shared_ptr<RectangleCertifier> bestRects, firstRects;
        std::shared_ptr<BlockCover> bestCover, firstCover;

        // A coarse search's best carrier, realised on the input domain.
        std::unique_ptr<Realisation> realisation;
        std::shared_ptr<SquareCarrier> realised;
        SquareCarrier::Report realisedReport;
        std::shared_ptr<RectangleCertifier> realisedRects;
        std::shared_ptr<BlockCover> realisedCover;

        // The search's answer on the input domain: best* for the fine
        // search, realised* for a coarse one. Null when there is none.
        const SquareCarrier *finalCarrier() const { return coarse ? realised.get() : best.get(); }
        const RectangleCertifier *finalRects() const { return coarse ? realisedRects.get() : bestRects.get(); }
        const BlockCover *finalCover() const { return coarse ? realisedCover.get() : bestCover.get(); }
        bool succeeded() const { return finalCover() && finalCover()->getReport().valid; }
        // The arbiter's score of the final answer, on the input domain.
        double finalScore = 0.0, finalDir = 0.0, finalSing = 0.0, finalShape = 0.0;

        std::vector<std::string> messages;
        double seconds = 0.0, secondsLoop = 0.0, secondsRealisation = 0.0;
    };

    struct Status {
        bool domainValid = false;
        bool carrierValid = false;    // the fine carrier of Stage 2
        bool coverValid = false;      // the chosen cover passed Sec. 13.2's gates
        int triangles = 0;
        int initialCells = 0;
        int chosen = -1;              // index of the winning search
        int finalCells = 0;           // cells of the carrier the chosen cover is on
        int bestBlocks = 0;
        double bestObjective = 0.0;
        double bestScore = 0.0;       // the arbiter's, = bestObjective with no field or shape term
        double secondsStage1 = 0.0, secondsStage2 = 0.0, secondsSearch = 0.0;
        double secondsField = 0.0;
        std::vector<std::string> messages;
    };

    ATLAS(std::shared_ptr<Mesh> mesh, const Options &opts);

    // Runs Stages 1-6. True when a valid conforming blocking came out (the
    // fine fallback counts; Status says which it was).
    bool run();

    const Status &getStatus() const { return status_; }
    const Options &getOptions() const { return opts_; }

    bool hasDomain() const { return static_cast<bool>(domain_); }
    bool hasCarrier() const { return !searches_.empty() && searches_[0]->initial; }
    bool hasCover() const { return status_.chosen >= 0; }

    const PlanarDomain &getDomain() const { return *domain_; }
    // The fine carrier as Stage 2 built it, before anything replaced a cell.
    const SquareCarrier &getInitialCarrier() const { return *searches_[0]->initial; }
    const SquareCarrier::Report &getInitialCarrierReport() const { return searches_[0]->initialReport; }

    int numSearches() const { return static_cast<int>(searches_.size()); }
    const Search &getSearch(int i) const { return *searches_[i]; }
    const Search &getChosen() const { return *searches_[status_.chosen]; }

    // The chosen answer, on the input domain.
    const SquareCarrier &getCarrier() const { return *getChosen().finalCarrier(); }
    const RectangleCertifier &getRectangles() const { return *getChosen().finalRects(); }
    const BlockCover &getCover() const { return *getChosen().finalCover(); }

    // The reference field, when Options::field.reference is DualMBO and it
    // was built; null otherwise.
    const ReferenceField *getField() const { return field_ && field_->built() ? field_.get() : nullptr; }

private:
    // Stages 3-6 on s.work; fills s.history and the best/first snapshots.
    void search(Search &s, bool rewrite, int rounds, int passes, int maxRings, bool anneal);
    void runCoarse(Search &s);
    // Realise s's layouts on the input, resuming after its earlier attempts.
    void realise(Search &s, int maxAttempts);
    // The arbiter's additions to a cover's objective for carrier C: the
    // field's terms and the shape term, each 0 when switched off.
    double arbiterTerms(const SquareCarrier &C, double *eDir, double *eSing, double *shape) const;

    std::shared_ptr<Mesh> mesh_;
    Options opts_;
    Status status_;
    std::unique_ptr<PlanarDomain> domain_;
    std::shared_ptr<ReferenceField> field_;
    std::vector<std::unique_ptr<Search>> searches_;
    // Searches still running while they run side by side; a realisation takes
    // its share of the cores for TMOP from it (see realise()).
    std::atomic<int> activeSearches_{0};
};

#endif // __ATLAS_HXX__
