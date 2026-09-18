#ifndef __ATLAS_HXX__
#define __ATLAS_HXX__

#include <memory>
#include <string>
#include <vector>

#include "ATLAS/BlockCover.hxx"
#include "ATLAS/CavityFill.hxx"
#include "ATLAS/CavityRewrite.hxx"
#include "ATLAS/ExplicitTemplates.hxx"
#include "ATLAS/PlanarDomain.hxx"
#include "ATLAS/RectangleCertifier.hxx"
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
// The contrast with MERIDIAN and TORSION is the spec's central idea: nothing
// here is a cross field to be integrated and quantised afterwards. The carrier
// already is a set of valid charts with integer transitions (Sec. 2.3), every
// later state is one too, and extraction is a finite integer certificate
// (Sec. 5) rather than streamline tracing. When simplification fails, what
// comes out is still a valid conforming blocking -- the fine fallback Sec. 14
// promises -- and the report says which it is.
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
        // Wall-clock budget for the whole Stage 4-6 loop, seconds.
        double timeBudget = 300.0;
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
        bool valid = false;
        int committed = 0;      // Stage 6 rewrites after this round's cover
        int rings = 0;          // Stage 6's cavity radius in this round
        double seconds = 0.0;
    };

    struct Status {
        bool domainValid = false;
        bool carrierValid = false;
        bool templatesValid = true;   // the carrier still validates after Stage 3
        bool rewritesValid = true;    // ... and after every Stage 6 round
        bool coverValid = false;      // the best cover passed Sec. 13.2's gates
        int triangles = 0;
        int initialCells = 0;
        int cellsAfterTemplates = 0;
        int finalCells = 0;           // cells of the carrier the best cover is on
        int bestRound = -1;
        int bestBlocks = 0;
        double bestObjective = 0.0;
        double secondsStage1 = 0.0, secondsStage2 = 0.0, secondsStage3 = 0.0, secondsLoop = 0.0;
        std::vector<Round> history;
        std::vector<std::string> messages;
    };

    ATLAS(std::shared_ptr<Mesh> mesh, const Options &opts);

    // Runs Stages 1-6. True when a valid conforming blocking came out (the
    // fine fallback counts; Status says which it was).
    bool run();

    const Status &getStatus() const { return status_; }
    const Options &getOptions() const { return opts_; }

    bool hasDomain() const { return static_cast<bool>(domain_); }
    bool hasCarrier() const { return static_cast<bool>(initialCarrier_); }
    bool hasTemplates() const { return static_cast<bool>(templates_); }
    bool hasCover() const { return static_cast<bool>(bestCover_); }
    bool hasRewrite() const { return static_cast<bool>(rewrite_); }

    const PlanarDomain &getDomain() const { return *domain_; }
    // The carrier as Stage 2 built it, before anything replaced a cell.
    const SquareCarrier &getInitialCarrier() const { return *initialCarrier_; }
    const SquareCarrier::Report &getInitialCarrierReport() const { return initialReport_; }
    // The carrier right after Stage 3.
    const SquareCarrier &getTemplatedCarrier() const { return *templatedCarrier_; }
    const ExplicitTemplates::Report &getTemplateReport() const { return templates_->getReport(); }
    // The carrier, certificates and cover of the best round.
    const SquareCarrier &getCarrier() const { return *bestCarrier_; }
    const RectangleCertifier &getRectangles() const { return *bestRects_; }
    const BlockCover &getCover() const { return *bestCover_; }
    // The first round's Stages 4 and 5, before any rewrite.
    const RectangleCertifier &getFirstRectangles() const { return *firstRects_; }
    const BlockCover &getFirstCover() const { return *firstCover_; }
    const CavityRewrite::Report &getRewriteReport() const { return rewrite_->getReport(); }

private:
    std::shared_ptr<Mesh> mesh_;
    Options opts_;
    Status status_;

    std::unique_ptr<PlanarDomain> domain_;
    std::unique_ptr<SquareCarrier> work_;             // the carrier Stages 3 and 6 rewrite
    std::unique_ptr<SquareCarrier> initialCarrier_;
    SquareCarrier::Report initialReport_;
    std::unique_ptr<SquareCarrier> templatedCarrier_;
    std::unique_ptr<ExplicitTemplates> templates_;
    std::unique_ptr<CavityRewrite> rewrite_;

    // Shared, because the first round may also be the best one.
    std::shared_ptr<SquareCarrier> bestCarrier_, firstCarrier_;
    std::shared_ptr<RectangleCertifier> bestRects_, firstRects_;
    std::shared_ptr<BlockCover> bestCover_, firstCover_;
};

#endif // __ATLAS_HXX__
