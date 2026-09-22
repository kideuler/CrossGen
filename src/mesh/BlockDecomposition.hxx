#ifndef __BLOCK_DECOMPOSITION_HXX__
#define __BLOCK_DECOMPOSITION_HXX__

#include <array>
#include <string>
#include <vector>

#include "mesh/Mesh.hxx"

// The quadrilateral block decomposition -- the object every pipeline's
// "Blocks" or "Patches" phase draws on the model as light-blue sides and
// green macrovertices, once it has one: a set of quadrilaterals, each of
// whose four sides is a macro edge running between two macrovertices, each
// macro edge a polyline on the model shared by at most two blocks.
//
// ATLAS reaches one of these by Sec. 6's exact cover over a square-transport
// carrier (ATLAS::BlockCover, docs/square_transport_2d_theory_and_
// implementation.md Sec. 13.2); MERIDIAN and TORSION by Stage 8's planar
// arrangement of the traced separatrices (MERIDIAN::Arrangement,
// docs/shepherd2022.pdf Sec. 4), the two pipelines sharing that stage and
// everything after it. The three methods -- and the meshes they run the
// search on, a certified square carrier for one, a cut triangulation for the
// other two -- have nothing in common. What they hand back, once found, does:
// a macro complex on the model, every side shared full-length with at most
// one neighbour, every macrovertex the meeting of as many sides as its own
// turning asks for. This class is that answer on its own, with no memory of
// which of the three produced it, so that drawing a decomposition, writing
// it out, or asking any question of "the blocks" rather than of one method's
// own bookkeeping is written once, here, rather than three times over.
//
// Each pipeline exposes an adapter that builds one of these from its own
// result -- BlockCover::blockDecomposition(), Arrangement::
// blockDecomposition() -- and nothing here knows how to build itself. A
// defect in the decomposition (a hanging vertex, a mismatched side, a face
// that did not come out four-sided) is still the source stage's own report
// to give: this class only holds a decomposition once it already is one, and
// a face that failed that test is simply not one of its blocks.
class BlockDecomposition {
public:
    struct MacroVertex {
        Point p{0.0, 0.0};
        bool onBoundary = false;
        bool onInterface = false;
        int valence = 0;   // number of incident macro edges
    };

    // A polyline on the model from macro vertex `from` to macro vertex `to`,
    // endpoints included. `blockA`/`sideA` names it as side `sideA` of block
    // `blockA`, walked corner to corner in the direction `from` -> `to`;
    // `blockB`/`sideB` is the block on its other side, walked the other way,
    // or -1 where there is none -- a side of dS, or a side its own source
    // stage could not match and so left out of the decomposition.
    struct MacroEdge {
        int from = -1, to = -1;
        std::vector<Point> points;

        int blockA = -1, sideA = -1;
        int blockB = -1, sideB = -1;
        bool boundary = false;
        bool interface = false;
        int matLeft = 0, matRight = 0;   // 0 where the source has no materials
    };

    // Four corners and four sides, both in cyclic order: side k runs
    // corners[k] -> corners[(k + 1) % 4], and is edges[k] read forwards
    // (flip[k] false -- this block is that edge's A side) or backwards
    // (flip[k] true -- its B side).
    struct Block {
        std::array<int, 4> corners{{-1, -1, -1, -1}};
        std::array<int, 4> edges{{-1, -1, -1, -1}};
        std::array<bool, 4> flip{{false, false, false, false}};
        int material = 0;
    };

    std::vector<MacroVertex> vertices;
    std::vector<MacroEdge> edges;
    std::vector<Block> blocks;

    // Which pipeline built this -- "ATLAS", "MERIDIAN", "TORSION" or "UMBER"
    // -- kept for messages only; nothing here branches on it.
    std::string source;

    // Side `side` of `block`, corner to corner, in that direction. Empty for
    // an out-of-range block or side.
    std::vector<Point> sidePolyline(int block, int side) const;

    // ── Sampling a block ────────────────────────────────────────────────────
    //
    // The three below are the only things here that compute rather than store,
    // and they are here rather than in a caller because every caller wants the
    // same three and wants them to agree. A block drawn, a block meshed and a
    // block asked which material it sits in have to be sampled identically or
    // the picture is of something other than the mesh; and "the u parameter of
    // side 0" has to mean the same thing to the two blocks that share side 0,
    // or their grids do not meet. Nothing here searches, fits or optimizes:
    // it is the decomposition read at a parameter, no more.

    // `n + 1` points along side `side` of `block`, corner to corner, at equal
    // arc length of the side's own polyline (not at equal parameter, which on
    // a polyline of uneven steps is a different and worse set of points).
    // Empty when the side has no polyline; `n < 1` is taken as 1.
    std::vector<Point> sampleSide(int block, int side, int n) const;

    // One point of side `side` at arc-length fraction `t` in [0, 1], corner to
    // corner. The two blocks sharing a macro edge get the same point for the
    // same place on it, since both read the one polyline and one of them
    // simply reads it backwards -- which is what makes a mesh built on this
    // conforming by construction rather than by tolerance.
    Point pointOnSide(int block, int side, double t) const;

    // The transfinite (Coons) blend of the block's four sides at (u, v) in
    // [0, 1]^2, with u along side 0 (corners 0 -> 1) and v along side 3
    // reversed (corners 0 -> 3). On the four edges of the square it reproduces
    // the sides exactly, so a grid built from it meets its neighbours' grids
    // on the shared sides whatever it does inside.
    Point coonsPoint(int block, double u, double v) const;

    // A point that really is inside block `b`'s outline, for asking a question
    // of the model underneath it -- which material, which triangle. The Coons
    // centre is tried first and is the answer for all but a badly non-convex
    // block; where it falls outside, a vertex of the outline is walked inward
    // instead. Returns false only for a block whose outline is degenerate.
    bool interiorPoint(int block, Point &out) const;

    // The macro edges as OBJ polylines, and the blocks as closed OBJ loops --
    // one pair of writers standing in for what used to be BlockCover::
    // writeOBJ and Arrangement::writeOBJ / writePatchOBJ.
    bool writeEdgesOBJ(const std::string &path) const;
    bool writeBlocksOBJ(const std::string &path) const;
};

#endif // __BLOCK_DECOMPOSITION_HXX__
