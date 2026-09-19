#ifndef __CAVITY_REWRITE_HXX__
#define __CAVITY_REWRITE_HXX__

#include <functional>
#include <map>
#include <string>
#include <vector>

#include <unordered_set>

#include "ATLAS/CavityFill.hxx"
#include "ATLAS/RectangleCertifier.hxx"
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
        double timeBudget = 60.0;   // seconds per round
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
    //         + wShape (cells below shapeFloor) + wAlign (misalignment),
    //
    // N_B the base complex's patch count (RectangleCertifier::basePatchCount),
    // w_v wDefect inside and wBoundaryDefect on dS, the last two terms the
    // proxies for geometry described with their options, and accepted by the
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
        int maxRings = 3;
        // Proposals realised per move, the locally best first.
        int realisations = 2;
        unsigned seed = 12345;
    };

    struct AnnealReport {
        int moves = 0, cavities = 0, proposals = 0, certified = 0;
        int accepted = 0, uphill = 0, improvements = 0;
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
    };

    // The cavity of `rings` rings round the seed vertices, never taking a cell
    // in `exclude`, with its one loop; false when it is not a usable disk.
    bool growCavity(const std::vector<int> &seeds, int rings, const std::unordered_set<int> &exclude,
                    std::vector<int> &cav, Loop &L, int &material, std::unordered_set<int> &inCav);

    void propose(const Loop &L, std::vector<Proposal> &out);
    bool realise(const Loop &L, const Proposal &p, CavityFill::Patch &P) const;
    // Side k of a corner set: its edge count now and the range drops and
    // insertions allow.
    void sideRange(const Loop &L, int from, int to, int &n, int &lo, int &hi) const;
    void reject(const std::string &why) { ++report_.rejections[why]; }
    int totalDefect() const;
    int irregularCount() const;

    SquareCarrier &C_;
    Options opts_;
    Report report_;
    AnnealReport annealReport_;
    std::vector<SquareCarrier> checkpoints_;
};

#endif // __CAVITY_REWRITE_HXX__
