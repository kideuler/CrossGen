#ifndef __CAVITY_REWRITE_HXX__
#define __CAVITY_REWRITE_HXX__

#include <map>
#include <string>
#include <vector>

#include "ATLAS/CavityFill.hxx"
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
// protected, designated or domain-boundary vertex strictly inside it and no
// template cell in it. Its boundary is kept exactly (CavityFill), and it is
// refilled with one of Sec. 8.2's families in their general form:
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

    CavityRewrite(SquareCarrier &carrier, const Options &opts);

    // One pass of Stage 6. `priority` vertices (witnesses Stage 4 recorded) are
    // tried first. Returns the number of rewrites committed.
    int round(const std::vector<int> &priority);

    const Report &getReport() const { return report_; }

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
        std::vector<int> sigma;       // star spoke counts
        int defect = 0;               // predicted defect of the loop + new interior
        double geo = 0.0;
    };

    void propose(const std::vector<int> &loop, const std::vector<double> &angle,
                 const std::vector<int> &outside, std::vector<Proposal> &out);
    bool realise(const std::vector<int> &loop, const Proposal &p, CavityFill::Patch &P) const;
    void reject(const std::string &why) { ++report_.rejections[why]; }
    int totalDefect() const;
    int irregularCount() const;

    SquareCarrier &C_;
    Options opts_;
    Report report_;
};

#endif // __CAVITY_REWRITE_HXX__
