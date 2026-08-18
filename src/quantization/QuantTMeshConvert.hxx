#ifndef _QUANT_TMESH_CONVERT_HXX_
#define _QUANT_TMESH_CONVERT_HXX_

#include "medialaxis/MedialAxisTMesh.hxx"
#include "quantization/QuantTMesh.hxx"
#include "tracing/QuadLayout.hxx"

#include <array>
#include <string>
#include <vector>

// ─── Adapters onto QuantTMesh ────────────────────────────────────────────────
//
// Two ways into the quantizer. A QuadLayout already is a combinatorial
// T-mesh -- its arcs are the edges and its faces list their sides -- so
// that conversion is a relabeling. A block list (MedialAxisTMesh::blocks,
// or any set of four-sided outlines) carries no shared-edge identifiers, so
// those are welded geometrically: block corners become nodes, each block
// side is split where another block's corner lands on it (the
// T-junctions), and coincident pieces are identified into one edge.
//
// In both cases `h` is the target parametric edge length: xIdeal is the
// geometric side length divided by h. Passing h <= 0 sets xIdeal = 1
// everywhere -- the pure block-decomposition setting, where Stage II of the
// quantizer shrinks every edge toward its minimum.
//
// Faces that are not four-sided (non-quad layout faces, cap blocks) are
// skipped and counted; their edges keep a PHANTOM row on that side, i.e.
// they behave like boundary there and impose no constraint.

struct QuadLayoutQuant {
    QuantTMesh tmesh;
    std::vector<int> edgeOfArc;      // layout arc -> tmesh edge, -1 if unused
    std::vector<int> faceOfLayout;   // layout face -> tmesh face, -1 skipped

    // Geometry for rendering, in the same shape BlockQuant carries it: a
    // layout arc already is the finest shared unit between its two faces, so
    // unlike the block welder above this needs no splitting -- one arc is
    // one tmesh edge and its polyline is that edge's geometry directly, in
    // the arc's own a -> b direction.
    std::vector<std::vector<Point>> edgeGeometry;  // tmesh edge -> polyline
    // Whether a side's edge runs against that a -> b direction, per face,
    // per side, per edge of that side: dart parity (see QuadLayout::Arc).
    std::vector<std::array<std::vector<char>, 4>> sideReversed;

    // A face that is not four-sided is not a T-mesh cell and cannot be
    // quantized, but it still occupies real ground in the model -- dropping
    // it silently would leave a hole in the picture, and if it sits on the
    // domain boundary that hole opens straight through the outline. Its true
    // outline (every arc of it, walked in face order and true geometry) is
    // kept here instead, so the renderer can draw it unsubdivided the same
    // way it already draws a block the quantizer collapsed to zero cells.
    std::vector<std::vector<Point>> skippedOutlines;

    int skippedFaces = 0;
    bool ok = false;
    std::string error;
};

QuadLayoutQuant makeQuantTMesh(const QuadLayout &layout, double h = 0.0);

struct BlockQuant {
    QuantTMesh tmesh;
    std::vector<Point> nodes;                     // welded corner points
    std::vector<std::array<int, 2>> edgeNodes;    // endpoints of each edge
    std::vector<std::vector<Point>> edgeGeometry; // polyline of each edge
    std::vector<int> faceOfBlock;                 // block -> face, -1 skipped
    std::vector<int> blockOfFace;                 // face -> block

    // Whether edgeGeometry runs against the walk direction of the side it
    // sits on, per face, per side, per edge of that side. An edge shared by
    // two faces is stored once and so runs backwards along one of them.
    std::vector<std::array<std::vector<char>, 4>> sideReversed;

    int skippedBlocks = 0;      // total of the four below
    int skippedNonQuad = 0;     // not four-sided even after splitting
    int skippedUnlocatable = 0; // a corner not found on its own outline
    int skippedDegenerate = 0;  // two of its own sides are the same curve
    // Blocks that would give some edge a third incident face. An edge of a
    // T-mesh borders at most two cells, so this means the decomposition
    // handed in overlaps itself -- on fan-like models the cap zones can
    // produce a sliver that is exactly the overlap of its two neighbours.
    // The late claimant is dropped so the rest of the model still stands.
    int skippedOverlapping = 0;
    bool ok = false;
    std::string error;
};

// `weldTol` is the distance below which two points are the same place; it
// must be well below the shortest block side. Blocks that are not
// four-sided (caps, and the triangular pieces of the blue and purple
// templates) are skipped and counted, as are blocks so degenerate that two
// of their sides weld onto the same edge.
BlockQuant makeQuantTMesh(const std::vector<TMeshBlock> &blocks,
                          double weldTol, double h = 0.0);

// A default tolerance of 1e-4 of the target block size: on the models in
// data/meshes the shortest genuine block side is a few 1e-3 of that size
// while the near-coincident corners a template occasionally emits are a few
// 1e-6 of it, so this sits with room on both sides.
inline BlockQuant makeQuantTMesh(const MedialAxisTMesh &tm, double h = 0.0) {
    return makeQuantTMesh(tm.blocks, 1e-4 * tm.targetSize(), h);
}

// One side of a face as a curve carrying its quantization.
//
// `pts` is the side's true geometry -- the polylines of its edges walked in
// order -- so nothing of the block's curvature is lost. `tickAt` gives the
// arc length of each quantization tick along it, sideSum + 1 of them: each
// edge is subdivided into its own x pieces, so two faces sharing an edge
// place identical ticks on it and their grids weld.
//
// Together these define a parametrization u in [0, 1] that maps tick i to
// u = i / (tickAt.size() - 1) and follows the real geometry in between --
// what transfinite interpolation needs to fill the block with curves that
// meet its curved sides.
struct SideCurve {
    std::vector<Point> pts;
    std::vector<double> arc;     // cumulative arc length at each pts entry
    std::vector<double> tickAt;  // arc length of each tick
    double length = 0.0;

    int cells() const { return static_cast<int>(tickAt.size()) - 1; }

    // The point at arc length `s`.
    Point atArc(double s) const;

    // The point at parameter u, which passes through the ticks.
    Point at(double u) const;

    // The same curve walked the other way, ticks included.
    SideCurve reversed() const;
};

// Requires a quantized tmesh.
SideCurve sideCurve(const BlockQuant &bq, int face, int side);
SideCurve sideCurve(const QuadLayoutQuant &lq, int face, int side);

#endif  // _QUANT_TMESH_CONVERT_HXX_
