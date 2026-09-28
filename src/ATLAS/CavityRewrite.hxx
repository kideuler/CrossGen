#ifndef __CAVITY_REWRITE_HXX__
#define __CAVITY_REWRITE_HXX__

#include <array>
#include <functional>
#include <map>
#include <random>
#include <string>
#include <vector>

#include <unordered_map>
#include <unordered_set>

#include "ATLAS/CavityFill.hxx"
#include "ATLAS/RectangleCertifier.hxx"
#include "ATLAS/ReferenceField.hxx"
#include "ATLAS/SquareCarrier.hxx"

// Stage 6 of docs/square_transport_2d_theory_and_implementation.md: rewrite
// difficult cavities.
//
// What stops Stages 4 and 5 from finding a coarse blocking is always the same
// thing: vertices whose valence is wrong for where they are. An interior vertex
// of valence 3 or 5 is a corner of every block at it (Sec. 4), and its block
// sides run straight on through regular vertices until they meet dS or another
// such vertex, so every irregular vertex the carrier inherited from the
// triangulation cuts the blocking up (Sec. 8.1). The obstruction witnesses of
// Sec. 8.3 -- a failed development, a conflict in the cover -- all sit on those
// vertices. So they are what this stage attacks, worst first.
//
// ### The move
//
// Around a witness v, a cavity is grown ring by ring: the cells at v, then
// every cell touching those, and so on (Sec. 8.3's "bounded cavities around
// these witnesses"). The cavity must be a disk of one material with no
// protected or designated vertex strictly inside it and no template cell in
// it. Its boundary is kept exactly (CavityFill) -- except that where it runs
// along dS it may be subdivided afresh, which is how a dislocation leaves the
// domain -- and it is refilled with one of Sec. 8.2's families in their
// general form:
//
//   * a structured grid, with four of the boundary vertices as its corners and
//     opposite sides of equal count. For a hexagonal cavity this is the
//     parity-changing 3 -> 2 move, for any cavity it is the move that leaves no
//     irregular vertex inside;
//   * a star, three or five blocks round one new centre -- the 3 <-> 3 flip on
//     a hexagon, and its generalisation.
//
// The corners are chosen by what they do to the boundary. A boundary vertex u
// of the cavity keeps its o_u outside cells and gets one inside cell if it is a
// corner of the fill, two otherwise, so its defect afterwards is
// |o_u + i_u - t_u| with t_u its regular valence (4 inside, from the angle on
// dS). Summed over the loop that is base + sum over corners of delta_u, O(1)
// per corner set, which is what makes enumerating every admissible corner set
// of a grid affordable. The fill with the lowest predicted energy
//
//     E = lambda_S * (total valence defect) + lambda_Q * (cells)
//
// (Sec. 8.3's score with N_S measured as defect and N_Q the cell count) is
// realised, certified by CavityFill::validate, and committed only if E drops.
// Moves that add cells are allowed -- Sec. 8.3 warns that insisting every move
// shrink the carrier traps the search -- but the defect has to pay for them.
//
// Committed cavities in one round are vertex-disjoint, so their predicted
// changes add up exactly; the round is applied as one batch.
//
// Sides that run along dS have a range of counts rather than one: droppable
// points may go and points may be inserted (CavityFill), so a corner set whose
// opposite sides differ can still be a grid, and a star's side can meet its
// spokes. Those corner sets are drawn from a pool of the likeliest corners
// (Options::cornerPool), besides every set whose sides already match.
//
// round() is the greedy half of the stage: fast, local, and blind to how many
// blocks the singularities it leaves will cut the domain into. anneal() is the
// other half, for a coarse carrier, and is described with its options.
//
// ### Across an interface
//
// A cavity of one material may not re-subdivide the interfaces on its loop,
// because the region across has cells on them (Sec. 8.4), and it cannot hold
// an interface inside it. On a multi-material domain that makes every
// interface a wall with its subdivision fixed for good: a singularity can
// never cross one, and a count mismatch between two interfaces -- the two
// sides of a thin layer -- can only be absorbed by singularities inside the
// layer, whose lines then cut every layer round it. That, not the templates,
// was most of the blocks ATLAS put on multimat/det_rocket, icf and rocket.
//
// Sec. 8.4's own remedy is "enlarge the cavity", and a *straddling* cavity
// does: cells of two materials either side of one run of interface, the
// chain, whose ends are on the cavity's loop and whose interior is inside it.
// Side A is filled first with the chain free to be re-subdivided as dS is --
// its input vertices kept, its midpoints droppable, new InterfaceSplit points
// inserted on its input segments -- and side B is then filled on its own loop
// with A's version of the chain fixed, so both sides meet on the same points.
// The two fills are one patch over the whole cavity, certified by
// CavityFill::settle with the chain's kept input vertices as interior anchors,
// and committed as one Edit with a material per cell. Everything validate()
// asks of a single-material fill holds of it, and the interface is the chain
// again by construction; SquareCarrier::validate() audits that it still lies
// on the input's interface segments.
class CavityRewrite {
public:
    struct Options {
        int maxRings = 3;
        int maxCavityCells = 1500;
        // A rewrite may not leave a cell worse than this, or worse than the
        // worst cell it replaces, whichever is lower.
        double minScaledJacobian = 0.2;
        int smoothingIterations = 40;
        bool grids = true;
        bool stars = true;
        // A fill corner must have a cavity angle below this (degrees): a
        // corner cell spanning nearly pi is nearly degenerate.
        double maxCornerAngle = 160.0;
        int realisationsPerWitness = 3;
        // Candidate corners per cavity for the corner sets whose sides are
        // matched by dropping or inserting points on dS.
        int cornerPool = 12;
        // At most this many points inserted on one side, as a multiple of it.
        double maxInsertFactor = 1.0;
        double lambdaS = 1.0;
        double lambdaQ = 0.002;
        // Sign-aware boundary defect in round()'s score (guidance note, R3),
        // relative to an interior unit: a +1/4 on dS (too few cells, a block
        // corner on a smooth run) and a -1/4 there. At 1 and 1 the score is
        // the plain defect count, and a move that carries an interior
        // singularity onto dS changes nothing but the cell count, so any move
        // that also saves cells is taken: that is how the field's interior
        // cones ended up as block corners on arcs before the annealer ever
        // saw the carrier (geom016: E_sing 2.8 -> 6 in two rounds). ATLAS
        // applies them on coarse carriers only (signedDefectsOnFine).
        double boundaryPlusWeight = 1.5;
        double boundaryMinusWeight = 0.5;
        // The reference field in round() (guidance note, Sec. 4.4), or null.
        // Among a witness's certified fills -- the first
        // realisationsPerWitness by score -- the one with the least
        // dE + fieldWeight Delta E_dir is committed, rather than the first:
        // equal-defect fills differ in which way their grid runs, and the
        // greedy rounds, not the annealer, fix most of a coarse carrier's
        // topology. fieldWeight is in defect units per unit of E_dir (the
        // annealer's wDir over its wDefect). coneStars: a star fill's centre
        // goes on a cone of its sign when one is in the cavity's kernel.
        const ReferenceField *field = nullptr;
        double fieldWeight = 5.0;
        bool coneStars = true;
        double timeBudget = 60.0;   // seconds per round
        // Cavities that straddle an interface ("Across an interface" above):
        // round() tries them at witnesses within straddleReach rings of an
        // interface, once the single-material cavities there found nothing.
        bool straddle = true;
        int straddleReach = 2;
        // A straddling cavity may take Stage 3 template cells (never a
        // designated vertex inside it). A single-material cavity never does:
        // there the greedy would pick templates apart for defect alone. Across
        // an interface it is the only way to change the counts a template
        // fixed on its sides, and on det_rocket every interface of the one
        // region no template covers is a template's side.
        bool straddleTemplates = true;
    };

    struct Report {
        int rounds = 0;
        int witnesses = 0;
        int cavities = 0;
        int proposals = 0;
        int realised = 0;
        int certified = 0;
        int committed = 0;
        int gridFills = 0, starFills = 0;
        int fieldChoices = 0;       // commits the field preferred over the first certified fill
        int conePlacedStars = 0;    // star fills centred on a cone
        int straddleCavities = 0, straddleCertified = 0, straddleCommitted = 0;
        int defectBefore = 0, defectAfter = 0;
        int irregularBefore = 0, irregularAfter = 0;
        int cellsBefore = 0, cellsAfter = 0;
        bool timedOut = false;
        std::map<std::string, int> rejections;
    };

    // Sec. 8.3's search proper: one move at a time, scored by the whole
    // carrier's energy
    //
    //     E = wBlocks N_B + sum_v w_v |q_v - t_v| + wCells N_Q
    //         + wShape (cells below shapeFloor) + wAlign (misalignment)
    //         [+ wDir E_dir + E_sing against a ReferenceField],
    //
    // N_B the base complex's patch count (RectangleCertifier::basePatchCount),
    // w_v wDefect inside and wBoundaryDefect on dS, the next two terms the
    // proxies for geometry described with their options, the last the
    // reference field's (docs/atlas_crossfield_guidance.md), and accepted by the
    // Metropolis rule at a temperature falling geometrically from T0 to T1
    // ("include inverse or temporarily enlarging moves; requiring every
    // elementary move to reduce the number of cells can trap the search").
    // The best state seen is what the carrier is left in ("maintain a
    // separate best feasible decomposition"). round() cannot see N_B at all
    // -- its score is local to the cavity -- so it removes singularities but
    // cannot align the ones that must stay; this is what moves them until
    // their separatrices meet.
    //
    // A move is one cavity, round a witness drawn at random (a vertex of
    // wrong valence mostly, any vertex sometimes), of one of three shapes:
    // rings round the witness, rings round one of its edges, or -- at a
    // vertex on dS or a protected one -- rings round part of its fan only.
    // The last is what lets a corner reach any valence: a whole ring leaves
    // a fill at most two cells at the witness, never the three a reflex
    // corner asks for. Among the cavity's proposals one of the locally best
    // few is drawn by a softmax, so moves that only shift a singularity --
    // the ones alignment needs -- are tried as well as the improving ones.
    //
    // Every trial is certified by CavityFill before it is scored, so an
    // accepted uphill move is still a valid carrier. It is meant for a
    // coarse carrier (CoarseDomain), where scoring a trial -- a copy, the
    // edit, and the base complex -- is a fraction of a millisecond.
    struct AnnealOptions {
        // Moves in the schedule; the temperature falls with the move count,
        // so the result depends on the seed alone. `seconds` is a safety cap.
        int maxMoves = 15000;
        double seconds = 60.0;
        double T0 = 1.5, T1 = 0.02;
        double wBlocks = 1.0;
        double wDefect = 2.0;
        double wBoundaryDefect = 3.0;
        double wCells = 0.002;
        // Cell shape: wShape (s0 - SJ)/s0 for every cell whose worst scaled
        // Jacobian SJ is below s0. A layout is only as good as the blocks it
        // can be realised with, and the coarse cells are their proxy.
        double wShape = 1.0;
        double shapeFloor = 0.35;
        // A trial cell below this is refused outright (a fixed floor: the
        // greedy rounds' "no worse than the cavity was" would let thousands
        // of moves ratchet the geometry down).
        double minScaledJacobian = 0.2;
        // Alignment: wAlign w (1 - cos 4(theta_cell - theta_field)) / 2 per
        // cell, against a reference cross field (angle and weight at a point;
        // unset = no term). Blocks, defect and shape are blind to which way a
        // block runs, and a layout whose grid lines cross a rectangle on the
        // diagonal scores as well as one parallel to its walls without it.
        std::function<double(const Point &, double *)> field;
        double wAlign = 0.5;
        // Sign-aware boundary defect (docs/atlas_crossfield_guidance.md, R3):
        // the weight of a unit of boundary defect with too few cells (+1/4:
        // on a smooth run of dS, one cell spanning pi, TFI worst SJ 0.09-0.19
        // on the four models that had them) and with too many (-1/4: three
        // cells of 60 degrees, 0.46-0.70). Negative = wBoundaryDefect, the
        // symmetric weight every layout was scored with before 2026-09-19.
        double wBoundaryDefectPlus = 4.5;
        double wBoundaryDefectMinus = 1.5;
        // The reference cross field's terms (guidance note, Sec. 4.4):
        // wDir E_dir + E_sing, both in blocks like N_B. Null = neither; it
        // replaces `field` above when both are set.
        const ReferenceField *reference = nullptr;
        double wDir = 10.0;
        ReferenceField::SingularityWeights singularity;
        // Field-directed moves (guidance note, R1): this fraction of the
        // witnesses is drawn from where the field disagrees with the carrier
        // -- a boundary defect, an interior singularity no cone accounts for,
        // the carrier vertex nearest a cone nothing sits on -- and a star
        // fill puts its centre on a cone of its own sign when one lies in the
        // cavity's kernel. 0 = the uniform draw above.
        double directedFraction = 0.0;
        // Charge a base patch that is not a disk -- a band round a hole that
        // no line cuts -- as the 1 + 2 (1 - chi) blocks it takes at least
        // (RectangleCertifier::basePatchCount's holeDeficit). Without it N_B
        // counts such a band as one block, and a carrier with no singular
        // vertex round a hole scores as a 3-block answer while Stage 5 can
        // only cover it with singletons (geom021 with the field terms on: N_B
        // 47 -> 3, blocks 176).
        bool patchTopology = true;
        int maxRings = 3;
        // Proposals realised per move, the locally best first.
        int realisations = 2;
        // On a carrier with interfaces, this fraction of the moves is a
        // cavity straddling one ("Across an interface"), grown round the
        // interface vertex nearest the witness.
        double straddleFraction = 0.35;
        // The annealer's cavities may take Stage 3 template cells and hold
        // designated vertices, which the greedy rounds' never may. Stage 3
        // commits one region at a time, and two templates that each
        // certified can still disagree about where their split points sit on
        // the interface between them, each then cutting the other's blocks
        // (multimat/geom004: three templates of 10 blocks, a cover of 31).
        // A designated vertex is a soft macrovertex; only the energy, which
        // counts the cuts it makes, can say which to give up. ATLAS turns
        // this on for domains with interfaces.
        bool freeTemplates = false;
        unsigned seed = 12345;
    };

    struct AnnealReport {
        int moves = 0, cavities = 0, proposals = 0, certified = 0;
        int accepted = 0, uphill = 0, improvements = 0;
        int directed = 0;             // moves whose witness the field chose
        int conePlacedStars = 0;      // star fills centred on a cone
        int straddleMoves = 0, straddleCertified = 0, straddleAccepted = 0;
        bool timedOut = false;
        double energyBefore = 0.0, energyAfter = 0.0;
        int blocksBefore = 0, blocksAfter = 0;
        double seconds = 0.0;
    };

    CavityRewrite(SquareCarrier &carrier, const Options &opts);

    // One pass of Stage 6. `priority` vertices (witnesses Stage 4 recorded) are
    // tried first. Returns the number of rewrites committed.
    int round(const std::vector<int> &priority);

    // The annealed search above; leaves the carrier in its best state.
    const AnnealReport &anneal(const AnnealOptions &ao);
    static double energy(const SquareCarrier &C, const AnnealOptions &ao, int *blocks = nullptr);

    const Report &getReport() const { return report_; }
    const AnnealReport &getAnnealReport() const { return annealReport_; }
    // The annealer's best state at each quarter of its schedule, oldest
    // first (the final best is the carrier itself): fallbacks for a caller
    // that cannot use the final one.
    const std::vector<SquareCarrier> &getCheckpoints() const { return checkpoints_; }

    // Stage 3 template groups (SquareCarrier::cellGroup) no cavity may take,
    // even where AnnealOptions::freeTemplates or Options::straddleTemplates
    // free the rest: ATLAS locks the O-grids of smooth inclusions, which are
    // already what such a region should be (five blocks, no corner lying flat
    // on the circle), and which the annealer otherwise pared down to fewer
    // blocks and worse ones (multimat/bubbles: TMOP worst 0.41 -> 0.36).
    void lockGroups(const std::vector<int> &groups);

    // The cavity radius, in rings round a witness. ATLAS widens it when a
    // round stalls: singular vertices that survive a small radius are the
    // isolated ones, and the only way to cancel an isolated 3-5 pair is a
    // cavity that holds both.
    int maxRings() const { return opts_.maxRings; }
    void setMaxRings(int k) { opts_.maxRings = k; }

private:
    struct Proposal {
        int kind = 0;                 // 0 grid, 1 star
        std::vector<int> corners;     // loop positions
        std::vector<int> counts;      // edges per side after drops and insertions
        std::vector<int> sigma;       // star spoke counts
        int defect = 0;               // predicted defect of the loop + new interior
        double weighted = 0.0;        // the same with Loop::weight
        int cells = 0;
        double geo = 0.0;
    };

    // One cavity loop and what may change on it: a droppable vertex may be
    // left out, and a loop edge on dS may take inserted points.
    struct Loop {
        std::vector<int> ids;
        std::vector<double> angle;
        std::vector<int> outside;
        std::vector<char> droppable;
        std::vector<char> onBoundary;  // loop edge i -> i+1 is on dS
        std::vector<double> weight;    // of each vertex's defect (anneal)
        // Sign-aware weights, when set: of a unit of defect with too few
        // cells, and with too many (AnnealOptions::wBoundaryDefectPlus).
        std::vector<double> weightPlus, weightMinus;
        // Only for the two sides of a straddling cavity, empty otherwise: the
        // target valence at each position (a new point of the chain has none
        // in the carrier yet), and the input feature edge each loop edge may
        // take inserted points on (-1: none), an interface one or dS.
        std::vector<int> target;
        std::vector<int> flexEdge;
        std::vector<char> flexIface;
        // Weighted defect of vertex u at valence `val`.
        double cost(int u, int val, int target) const {
            if (!weightPlus.empty()) {
                return weightPlus[u] * std::max(0, target - val) + weightMinus[u] * std::max(0, val - target);
            }
            return (weight.empty() ? 1.0 : weight[u]) * std::abs(val - target);
        }
    };

    // The cavity of `rings` rings round the seed vertices, never taking a cell
    // in `exclude`, with its one loop; false when it is not a usable disk.
    bool growCavity(const std::vector<int> &seeds, int rings, const std::unordered_set<int> &exclude,
                    std::vector<int> &cav, Loop &L, int &material, std::unordered_set<int> &inCav);

    void propose(const Loop &L, std::vector<Proposal> &out);
    // `loopOut`: the realised loop, corner 0 first, after drops and insertions.
    bool realise(const Loop &L, const Proposal &p, CavityFill::Patch &P, std::vector<int> *loopOut = nullptr) const;

    // A cavity straddling one run of interface: side A (one material) and
    // side B (the other) whose only contact is the chain p ... q, inside the
    // cavity but for its two ends. LA is side A's loop with the chain free to
    // be re-subdivided like dS; side B is filled afterwards against whatever
    // A made of the chain.
    struct Straddle {
        std::vector<int> cav, cavA, cavB;
        std::unordered_set<int> inD, inA, inB;
        int matA = -1, matB = -1;
        int p = -1, q = -1;
        Loop LA;
        CavityFill::Boundary BB;
    };
    // Grow `rings` rings round interface vertex `seed` over both materials
    // there; false when the result is not a straddling cavity as above.
    bool growStraddle(int seed, int rings, Straddle &S);
    // Fill side A with pA and then side B with the best of its proposals
    // against A's chain, in one patch (A's cells first, nA of them), and
    // certify it; `keep` gets the chain's reused interior vertices.
    bool fillStraddle(const Straddle &S, const Proposal &pA, double floor, CavityFill::Patch &P, int &nA,
                      std::vector<int> &keep);
    // An interface vertex eligible to seed a straddling cavity within
    // `reach` rings of v, nearest first, -1 if none.
    int straddleSeed(int v, int reach, std::mt19937 *rng) const;
    // The change in defect a certified straddling fill makes, weighted as
    // round()'s score (dS vertices by boundaryPlus/MinusWeight when signed),
    // over every vertex of the cavity and the patch.
    double straddleDefectDelta(const Straddle &S, const CavityFill::Patch &P, bool signedBoundary) const;
    // Side k of a corner set: its edge count now and the range drops and
    // insertions allow.
    void sideRange(const Loop &L, int from, int to, int &n, int &lo, int &hi) const;
    void reject(const std::string &why) { ++report_.rejections[why]; }
    int totalDefect() const;
    int irregularCount() const;

    // energy(), with E_dir's per-cell terms taken from misalignment_ when the
    // cell was seen before (the annealer's calls), or computed afresh (cache
    // null).
    static double energy(const SquareCarrier &C, const AnnealOptions &ao, int *blocks, CavityRewrite *cache);

    // E_dir's per-cell terms during anneal(), by the bit patterns of the
    // cell's four corners. A move rewrites one cavity and leaves every other
    // cell where it was, so most of the terms the next trial's energy adds up
    // were already computed for the last one. The sum is still taken over
    // every cell in carrier order, so it is the one ReferenceField::
    // directionEnergy() forms, to the bit.
    struct CornerBits {
        std::array<unsigned long long, 8> bits;
        bool operator==(const CornerBits &o) const { return bits == o.bits; }
    };
    struct CornerBitsHash {
        size_t operator()(const CornerBits &k) const;
    };
    double directionEnergy(const SquareCarrier &C, const ReferenceField &F);
    std::unordered_map<CornerBits, double, CornerBitsHash> misalignment_;
    const ReferenceField *misalignmentOf_ = nullptr;

    SquareCarrier &C_;
    Options opts_;
    Report report_;
    AnnealReport annealReport_;
    std::vector<SquareCarrier> checkpoints_;
    // During anneal() with AnnealOptions::freeTemplates.
    bool freeTemplates_ = false;
    std::vector<char> lockedGroup_;
    bool lockedCell(int q) const {
        const int g = C_.cellGroup[q];
        return C_.cellOrigin[q] == SquareCarrier::CellOrigin::Template && g >= 0 &&
               g < static_cast<int>(lockedGroup_.size()) && lockedGroup_[g];
    }
    // During anneal() with directed moves: the cones a star fill may centre
    // on (realise()).
    const std::vector<ReferenceField::Singularity> *starTargets_ = nullptr;
    mutable int conePlacedStars_ = 0;
};

#endif // __CAVITY_REWRITE_HXX__
