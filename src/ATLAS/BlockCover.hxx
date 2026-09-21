#ifndef __BLOCK_COVER_HXX__
#define __BLOCK_COVER_HXX__

#include <array>
#include <string>
#include <unordered_set>
#include <vector>

#include "ATLAS/RectangleCertifier.hxx"
#include "ATLAS/SquareCarrier.hxx"
#include "mesh/BlockDecomposition.hxx"

// Stage 5 of docs/square_transport_2d_theory_and_implementation.md: select a
// conforming cover.
//
// Sec. 6 states the problem as an exact cover with pairwise exclusions: one
// z_P in {0,1} per candidate, sum_{P ni q} z_P = 1 for every carrier cell, and
// z_P + z_R <= 1 for every incompatible pair, minimising
//
//     sum_P z_P [1 + alpha D(P) + beta C(P)].
//
// Incompatible means the macro complex would not be conforming: a full side of
// one block meeting only part of a side of the other, two blocks touching
// along more than one side, or a corner of one sitting inside a side of the
// other (which also catches blocks that touch at a single vertex). With those
// excluded, every block side is either all domain boundary or exactly one
// full side of exactly one neighbour -- which is Sec. 13.2's "shared-edge
// identity ... no hanging macrovertices" -- and compatible() is written to
// test precisely that, pair by pair.
//
// Sec. 6 is candid that minimising this exactly is hard, and asks for bounded
// enumeration and local cover updates with a feasible incumbent throughout.
// That is what this does:
//
//   * incumbents: all singletons (always feasible -- the carrier is
//     conforming), and each base-complex cover Stage 4 produced, kept if it
//     passes the gates below;
//   * improvement: a candidate P replaces the blocks now covering its cells
//     whenever those blocks lie entirely inside P (so the cover stays exact),
//     P is compatible with every block it touches, and the objective drops.
//     Sweeps repeat until nothing improves.
//
// The result is then audited as a macro complex, independently of how it was
// chosen (Sec. 13.2): coverage, side pairing, hanging vertices, protected
// vertices as macrovertices, interfaces on block sides, and Sec. 4's identity
// on the macrovertices. A search limit never returns a partial partition; the
// worst it returns is the incumbent.
class BlockCover {
public:
    struct Options {
        double alpha = 0.05;       // weight of D(P), the distortion
        double beta = 0.05;        // weight of C(P), the side complexity
        int maxSweeps = 20;
        bool singleContact = true; // Sec. 6: forbid two blocks sharing two sides
    };

    struct Block {
        RectangleCertifier::Certificate cert;
        double cost = 0.0;
    };

    struct MacroEdge {
        std::vector<int> chain;    // carrier vertices
        int blockA = -1, sideA = -1;
        int blockB = -1, sideB = -1;
        bool boundary = false;
        bool interface = false;
    };

    struct Report {
        int candidates = 0;
        int blocks = 0;
        int singletonBlocks = 0;
        int largestBlock = 0;
        int mapPieces = 0;              // carrier cells: the pieces of the block maps
        int macroVertices = 0;
        int macroEdges = 0;
        int irregularMacroVertices = 0;
        double objective = 0.0;
        double meanDistortion = 0.0, maxDistortion = 0.0;
        double meanComplexity = 0.0;
        // The gates.
        int uncovered = 0, overcovered = 0;
        int sideMismatches = 0;
        int hangingVertices = 0;
        int multipleContacts = 0;
        int protectedNotMacro = 0;
        int interfaceInside = 0;
        int uncertified = 0;
        int eulerLHS = 0, eulerRHS = 0;
        bool eulerHolds = false;
        bool conforming = false;
        bool valid = false;
        std::string incumbent;
        int swaps = 0, swapTrials = 0, sweeps = 0;
        std::vector<std::string> messages;
    };

    BlockCover(const RectangleCertifier &rects, const Options &opts);

    const SquareCarrier &getCarrier() const { return C_; }
    const std::vector<Block> &getBlocks() const { return blocks_; }
    const std::vector<MacroEdge> &getMacroEdges() const { return macroEdges_; }
    const std::vector<int> &getMacroVertices() const { return macroVertices_; }
    const Report &getReport() const { return report_; }

    // The macro complex as the shared block-decomposition representation
    // (mesh/BlockDecomposition.hxx): Sec. 13.2's gates already ran in
    // analyze(), so this only reshapes what passed them into the class every
    // pipeline's Blocks/Patches phase draws and can be written out with.
    BlockDecomposition blockDecomposition() const;

    // The macro edges as OBJ polylines over the carrier's vertices (Sec. 13.1).
    bool writeOBJ(const std::string &path) const;

    // Sec. 6's exact-cover program over this cover's candidates -- every
    // block Stage 4 certified, then every singleton cell -- for an external
    // solver: z_P in {0,1}, sum over P containing q of z_P = 1 for every
    // cell, z_P + z_R <= 1 for every incompatible pair, minimise
    // sum z_P cost(P). Only pairs that touch without sharing a cell are
    // listed; two that share a cell are excluded by the equalities already.
    // This class solves the same program by local improvement from the base
    // complex; the model is what an exact (CP-SAT, MIP) backend would take.
    struct Model {
        int cells = 0;
        std::vector<std::vector<int>> candidateCells;
        std::vector<double> cost;
        std::vector<std::pair<int, int>> conflicts;
    };
    Model exactCoverModel() const;
    // As text: "cells candidates", one "cost k c_1 .. c_k" line per candidate,
    // "conflicts", one "i j" line per pair; then the objective and block
    // count this cover reached.
    bool writeModel(const std::string &path) const;

private:
    struct Geo {
        std::array<std::vector<int>, 4> sideEdges;
        std::unordered_set<int> sideInterior;
        std::array<int, 4> corners{{-1, -1, -1, -1}};
    };
    Geo geometry(const RectangleCertifier::Certificate &cert) const;
    RectangleCertifier::Certificate singleton(int q) const;
    double cost(const RectangleCertifier::Certificate &cert) const;
    // Pairwise compatibility; `inR(q)` says whether a carrier cell is in R.
    template <class InP, class InR>
    bool compatible(const Geo &P, InP inP, const Geo &R, InR inR) const;
    // Audit a complete selection (Sec. 13.2) and fill the report's gates.
    bool analyze(const std::vector<RectangleCertifier::Certificate> &sel, Report &rep,
                 std::vector<MacroEdge> *edges, std::vector<int> *mverts) const;

    const RectangleCertifier &R_;
    const SquareCarrier &C_;
    Options opts_;
    std::vector<Block> blocks_;
    std::vector<MacroEdge> macroEdges_;
    std::vector<int> macroVertices_;
    Report report_;
};

#endif // __BLOCK_COVER_HXX__
