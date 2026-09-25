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

    // ── Comparing decompositions ────────────────────────────────────────────
    //
    // Three scores from docs/block_decomposition_metrics.md, each a weighted
    // form of one of its unweighted counts, so that five methods' answers on
    // one model can be put in one column. They read the decomposition as the
    // method drew it: its topology, which smoothing does not change (Sec. 2 of
    // the note), and its own polylines, which smoothing does -- the note's
    // after-smoothing reading needs the meshed block grids, which this class
    // does not hold. A block with a side that has no polyline is left out of
    // all three, as BlockQuadMesh leaves it unmeshed.
    //
    // All three are scores in [0, 1]: 0 is the best a layout can do on that
    // count, and 1 is the worst, which is also what any decomposition that
    // fails covers() gets. Two of the raw counts have no upper bound, so each
    // is bent into [0, 1) by a map named below; the maps are monotone, so they
    // rank exactly as the raw counts do, but a difference near 1 means less
    // than the same difference near 0.
    //
    // Two of them share one notion, the *sector*: at a macrovertex, a run of
    // block corners joined across the macro edges between them, stopping at an
    // edge of dS, at a material interface, and at an edge with no block on its
    // far side. So an interior vertex is one sector of 2 pi, a boundary vertex
    // is one sector of the model's own angle there, and a vertex on an
    // interface is one sector per material, the interfaces acting as boundary
    // exactly as T3 has them. Its angle is the sum of its block corners, each
    // corner measured between the tangents of its two sides, each tangent read
    // from the first tenth of the side (B1) by a least-squares quadratic, so
    // that a side which bends -- dS round a hole, a ring of an O-grid -- does
    // not tilt the tangent by the chord's sagitta. The .cxx has the numbers.

    // T2 and T4 in one number, each irregular vertex weighted by how far it is
    // from regular rather than counted once. A sector of angle theta split into
    // n blocks scores |2 theta / pi - n| -- the number of quarter turns its
    // block count is off by -- and the score is the sum over every sector.
    // An interior vertex of valence d scores |4 - d|, so valence 3 and 5 count
    // one and valence 6 counts two, which is the note's reason for keeping high
    // valence visible (it costs twice the angle). A smooth point of dS scores
    // |3 - d|. A model corner scores its fractional mismatch: a 90-degree
    // corner in one block 0 and in two blocks 1, and a 135-degree corner 1/2
    // either way, which settles T4's "ambiguous corner" by weight instead of
    // by exemption. A boundary sector within 20 degrees of pi is taken as a
    // smooth point (G3's 160-degree rule), since a polygonal dS turns a little
    // at every vertex and that turning is the mesh's, not the model's.
    // Returned as W / (W + N), W that sum and N the number of sectors, so it is
    // the defect per sector and compares across models of different sizes: 0
    // with no irregular vertex, 1/2 at one quarter turn per sector. 0 is not
    // reachable on most models -- the model's own corners force some defect
    // (the quarter disk owes exactly one quarter turn, T3) -- so compare
    // methods on one model rather than reading it against 0.
    double weightedIrregularity() const;

    // Every block has all four sides, and every side is either on dS or shared
    // with a second block: no gap is left anywhere against a block. The three
    // metrics return 1, the worst score, when this fails, so that a method which leaves part of the
    // model out cannot win on what it left out (G2). A decomposition with no
    // blocks does not cover either.
    bool covers() const;

    // B1, area-weighted: the RMS difference between each block
    // corner and the ideal angle of its sector, theta / n. Measuring against
    // theta / n rather than 90 degrees leaves out what the valence forces --
    // weightedIrregularity() has that -- and keeps only how evenly the method
    // spread its blocks round each vertex. Each corner is weighted by its
    // block's area: at one target edge length a block holds elements in
    // proportion to its area, and a corner angle is inherited by the whole of
    // its block's transfinite grid, not by one element, so this is the
    // deviation an average element sits in. Returned over a right angle and
    // capped at 1, 90 degrees being where a corner has folded flat or split a
    // right angle wrongly by all of it: 0 when every vertex splits its sector
    // evenly.
    double weightedCornerDeviation() const;

    // T5, area-weighted: the span ratio S = (longest macro edge) / (shortest)
    // of each chord -- every macro edge of a chord carries one interval count,
    // so S is the ratio of element lengths it forces from one end of the chord
    // to the other -- as a geometric mean over the chords, each weighted by
    // the area of the blocks it runs through (a block twice, if its chord
    // crosses it both ways). A chord's S is paid by every element along it, so
    // a thin channel coupled to a wide region costs in proportion to how much
    // of the model it drags, and a short chord between two nearly equal sides
    // costs almost nothing however large its ratio. Chords are the classes
    // BlockQuadMesh::assignIntervals() builds (side 0 with 2, 1 with 3, over
    // every block), so they are the ones that actually share a count.
    // Returned as 1 - 1/S of that mean, which needs no chosen scale: 0 when
    // every chord's sides are equal, 1/2 when elements along a typical chord
    // differ in size by a factor of two.
    double weightedChordSpan() const;

    // The macro edges as OBJ polylines, and the blocks as closed OBJ loops --
    // one pair of writers standing in for what used to be BlockCover::
    // writeOBJ and Arrangement::writeOBJ / writePatchOBJ.
    bool writeEdgesOBJ(const std::string &path) const;
    bool writeBlocksOBJ(const std::string &path) const;
};

#endif // __BLOCK_DECOMPOSITION_HXX__
