#ifndef __BLOCK_MESH_HXX__
#define __BLOCK_MESH_HXX__

#include <array>
#include <string>
#include <vector>

#include "ATLAS/BlockCover.hxx"
#include "ATLAS/SquareCarrier.hxx"
#include "mesh/Mesh.hxx"

// A quadrilateral mesh on the blocks of an ATLAS cover, at a target edge
// length: Sec. 11.2 of docs/square_transport_2d_theory_and_implementation.md
// ("assigning mesh counts after blocking") and Sec. 11.3 ("distinguish exact
// chart geometry from sampled quads"). It is ATLAS's counterpart of Stage 10 of
// MERIDIAN and TORSION (src/MERIDIAN/QuadMesh.hxx), and deliberately the same
// construction, so that the three methods' meshes differ only in the layout
// they were built on:
//
//   * one integer count per chord, chosen against the target edge length from
//     the isoparametric line lengths of the blocks the chord runs through;
//   * the nodes of every macro edge placed once, at equal arc length, and
//     handed to both blocks that share it;
//   * each block's interior by transfinite interpolation of its four sides;
//   * Winslow smoothing, boundary held, of the blocks that came out folded.
//
// ### What is different: there is nothing to fit
//
// Stage 10 meshes Stage 9's bicubic patches, which were fitted to traced
// curves. A certified ATLAS block needs no fit: it comes with its map (Sec.
// 5.3). Its n_u x n_v carrier cells are the pieces of a piecewise-bilinear
// chart from [0, n_u] x [0, n_v] onto the block, orientation-preserving and
// injective because the carrier passed SquareCarrier::validate(). So where
// Stage 10 blends the four sides' *spline parameters* and evaluates the patch
// there, this blends the four sides' *chart coordinates* and evaluates the
// chart (Options::useChart) -- the same transfinite interpolation, taken in
// the block's own coordinates rather than in the plane. Off, the interior is
// the discrete Coons blend of the four sides' points, exactly Stage 10's mode
// with the splines off.
//
// ### Counts (Sec. 11.2)
//
// Opposite sides of one block carry one count, and a shared side is one
// macro edge, so union-find over the macro edges -- side 0 with side 2 and side
// 1 with side 3 of every block -- gives the classes, the chords. Sec. 11.2
// minimises sum_j w_j (n - r_j)^2 over a class; with r_j = S_j / h and
// w_j = 1 / r_j^2 that is, to first order, the sum of squared *relative*
// errors, and Stage 10 minimises its exact form
//
//     F(N) = sum_k ( log( S_k / (N h) ) )^2,
//
// at whichever of floor(N*) and ceil(N*) is lower, N* = geomean(S_k) / h. The
// same rule here, so a target means the same thing to all three methods. S_k
// is the mean length of the block's isoparametric lines across the chord's
// direction, read off the carrier's own grid lines -- which are exact
// polylines here rather than samples of a surface.
//
// ### Sampling a chart is not reproducing it (Sec. 11.3)
//
// Three things are validated on the result rather than assumed, because Sec.
// 11.3 names them as what sampling can break:
//
//   * the four corner Jacobians of every element: a straight-sided element
//     spanning several pieces of a curved chart can fold even though the
//     chart does not (Report::invertedQuads, minScaledJacobian);
//   * the global incidence: every element edge used once or twice, no two
//     vertices at one point (conforming, cracks);
//   * the boundary error: a boundary node lies on the input's own polyline,
//     since a boundary macro edge *is* that polyline between two corners
//     (Sec. 1.1's preservation, which validate() audited), but the element
//     edge between two nodes is a chord of it and skips the input vertices in
//     between. Report::boundaryDeviation is the furthest any input boundary
//     vertex lies from the mesh boundary, and interfaceDeviation the same for
//     the material interfaces. Protected corners are never skipped: every one
//     is a block corner (Sec. 5.2), so it is a node.
//
// The chart itself remains the authoritative geometry (Sec. 11.3's last
// sentence); this is the sampled mesh a solver or a smoother is handed.
class BlockMesh {
public:
    struct Options {
        // Target edge length in the units of the model, as Stage 10's.
        double targetEdgeLength = 0.05;
        // Floor and (optional) ceiling on the edges per chord; 0 = no ceiling.
        int minIntervals = 1;
        int maxIntervals = 0;
        // Interior nodes through the block's certified chart (on), or the
        // discrete Coons blend of its four sides (off). See above.
        bool useChart = true;
        // Winslow sweeps over the interior nodes of the blocks whose worst
        // element is below smoothingThreshold (0: only folded ones), the
        // boundary held; Stage 10's defaults and Stage 10's reason for them.
        int smoothingPasses = 500;
        double smoothingTolerance = 1e-4;   // of the target edge length
        double smoothingThreshold = 0.0;
        // Two vertices closer than this fraction of the model's diagonal are a
        // crack.
        double crackTolerance = 1e-9;
    };

    // One class of Sec. 11.2's relation: the macro edges one chord of the
    // block complex cuts across, and the one count they all take.
    struct Chord {
        std::vector<int> arcs;          // BlockCover::getMacroEdges() indices
        int intervals = 0;
        double idealIntervals = 0.0;    // N*, before rounding
        double minLength = 0.0, maxLength = 0.0;   // of its macro edges
        double minSpan = 0.0, maxSpan = 0.0;       // of the isolines it was chosen from
        bool clamped = false;           // the floor or the ceiling, not F(N)
    };

    // One meshed block: a structured (ns+1) x (nt+1) grid of vertex indices in
    // the block's own frame, row-major, vert[j * (ns+1) + i] at chart
    // coordinates (i/ns n_u, j/nt n_v) on its sides. The same layout as Stage
    // 10's QuadMesh::Block, so one drawing routine serves both.
    struct Block {
        int block = -1;                 // BlockCover::getBlocks() index
        int ns = 0, nt = 0;
        std::vector<int> vert;
        int material = 0;
    };

    struct Report {
        double modelExtent = 1.0;       // diagonal of the input's bounding box
        double target = 0.0;

        int chords = 0;
        int arcsAssigned = 0;           // macro edges
        int minIntervals = 0, maxIntervals = 0;
        int clampedChords = 0;
        double meanIntervals = 0.0;

        int blocks = 0;                 // blocks meshed
        int unmeshedBlocks = 0;         // a side with no macro edge: the cover was not conforming

        int vertices = 0;
        int quads = 0;

        double minEdge = 0.0, maxEdge = 0.0, meanEdge = 0.0;
        double edgeRatioRms = 0.0;      // rms of log(length / target) over the edges
        double worstEdgeRatio = 1.0;

        double minScaledJacobian = 0.0, meanScaledJacobian = 0.0;
        double minScaledJacobianBefore = 0.0;   // before the Winslow pass
        int invertedBefore = 0;                 // ... and its invertedQuads
        int smoothedBlocks = 0;
        int smoothingSweeps = 0;
        int reflexCorners = 0;          // block corners the mesh turns through more than pi at
        // Elements with a corner determinant that is not positive (Sec. 9.1),
        // which includes every one of non-positive area.
        int invertedQuads = 0;
        double minQuadArea = 0.0, maxQuadArea = 0.0;
        double meshArea = 0.0;
        double domainArea = 0.0;        // the input's, for the chording deficit

        // Sec. 11.3's boundary error, absolute; and at how many nodes.
        double boundaryDeviation = 0.0;
        double interfaceDeviation = 0.0;
        int boundaryNodes = 0;

        int interiorEdges = 0, boundaryEdges = 0, nonManifoldEdges = 0;
        int cracks = 0;
        int materials = 0;
        int interfaceEdges = 0;         // element edges between two materials

        bool conforming = false;        // no third use of an edge, no cracks
        bool valid = false;             // ... every block meshed, nothing folded
        std::vector<std::string> messages;
    };

    // The cover must be of `cover.getCarrier()` and should be valid; a block
    // with a side that is not one macro edge is counted and left out.
    BlockMesh(const BlockCover &cover, const Options &opts);

    const std::vector<Point> &vertices() const { return verts_; }
    const std::vector<std::array<int, 4>> &quads() const { return cells_; }
    const std::vector<int> &quadMaterials() const { return cellMaterial_; }
    const std::vector<Block> &blocks() const { return grids_; }
    const std::vector<Chord> &chords() const { return chords_; }
    // Per macro edge: its count, and its nodes from the chain's first vertex.
    const std::vector<int> &arcIntervals() const { return intervals_; }
    const std::vector<std::vector<int>> &arcVertices() const { return arcNodes_; }
    const Report &getReport() const { return report_; }
    const Options &getOptions() const { return opts_; }

    // An .obj of quadrilateral faces, `usemtl mat<id>` per material as
    // mesh::QuadMesh reads it back.
    bool writeOBJ(const std::string &path) const;

private:
    void buildGrids(const BlockCover &cover);
    void assignIntervals(const BlockCover &cover);
    void meshArcs(const BlockCover &cover);
    void meshBlocks(const BlockCover &cover);
    void smooth();
    void check(const BlockCover &cover);

    const SquareCarrier &C_;
    Options opts_;
    Report report_;

    // Per block: its carrier vertex at each integer chart point, and the
    // macro edge on each side with whether the edge's chain runs the side's
    // way (counter-clockwise round the block).
    std::vector<std::vector<int>> chart_;
    std::vector<std::array<int, 4>> sideArc_;
    std::vector<std::array<char, 4>> sideForward_;

    std::vector<int> chordOf_;                    // per macro edge
    std::vector<Chord> chords_;
    std::vector<int> intervals_;                  // per macro edge
    std::vector<std::vector<int>> arcNodes_;      // per macro edge, intervals + 1
    std::vector<std::vector<double>> arcParams_;  // ... their positions along the chain, in carrier edges

    std::vector<Point> verts_;
    std::vector<std::array<int, 4>> cells_;
    std::vector<int> cellMaterial_;
    std::vector<Block> grids_;
};

#endif // __BLOCK_MESH_HXX__
