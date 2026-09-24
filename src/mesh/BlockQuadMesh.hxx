#ifndef __BLOCK_QUAD_MESH_HXX__
#define __BLOCK_QUAD_MESH_HXX__

#include <array>
#include <string>
#include <vector>

#include "mesh/BlockDecomposition.hxx"
#include "mesh/Mesh.hxx"

// A quadrilateral mesh on a BlockDecomposition, at a target edge length: the
// method-agnostic counterpart of MERIDIAN's Stage 10 (src/MERIDIAN/QuadMesh.hxx)
// and of ATLAS's BlockMesh (src/ATLAS/BlockMesh.hxx), written against the
// shared decomposition and nothing else.
//
// The construction is deliberately the same one all three use, so that two
// methods meshed at one target differ in their layout and in nothing else:
//
//   * one integer count per chord, chosen against the target edge length from
//     the isoparametric line lengths of the blocks the chord runs through;
//   * the nodes of every macro edge placed once, at equal arc length, and
//     handed to both blocks that share it;
//   * each block's interior by transfinite interpolation of its four sides;
//   * Winslow smoothing, boundary held, of the blocks that came out folded.
//
// ### Why this one exists when two already do
//
// Stage 10 and BlockMesh each need something their own pipeline has and the
// decomposition does not carry: Stage 10 blends the four sides' *spline
// parameters* and evaluates Stage 9's fitted bicubic patch there, BlockMesh
// blends the four sides' *chart coordinates* and evaluates the block's
// certified piecewise-bilinear chart. Both are exact geometry inside the
// block, and neither survives the trip through BlockDecomposition, which holds
// the sides and not the surface between them.
//
// A method that has no such interior geometry -- UMBER, whose blocks are the
// faces a motorcycle graph cut out of the model and whose only exact data are
// the traced iso-lines themselves -- has nothing to evaluate and needs none of
// it. The interior is then the discrete Coons blend of the four sides, which
// is exactly BlockMesh's Options::useChart = false and Stage 10's mode with
// the splines off. That is what this class is: their shared construction with
// the one part that differs between them left out, on the one input all three
// can produce.
//
// It is therefore not a replacement for either. Meshing a MERIDIAN layout
// through here rather than through Stage 10 would throw away the bicubic
// patches Stage 9 fitted and mesh the polylines they were fitted to, which is
// a worse mesh of the same layout.
//
// ### What is checked on the result
//
// The same three things BlockMesh checks, for the same reasons (its Sec. 11.3
// note), read back off the mesh rather than assumed:
//
//   * the four corner Jacobians of every element, since a straight-sided
//     element spanning a curved block can fold even where the block does not;
//   * the global incidence -- every element edge used once or twice, no two
//     vertices at one point;
//   * the boundary error: a boundary node lies on the macro edge's own
//     polyline, because that polyline *is* the model's boundary between two
//     macrovertices, but the element edge between two nodes is a chord of it
//     and skips whatever the polyline did in between. Report::boundaryDeviation
//     is the furthest any polyline vertex lies from that chord, and
//     interfaceDeviation the same for the material interfaces.
class BlockQuadMesh {
public:
    struct Options {
        // Target length of a mesh edge, in the units of the model -- the same
        // number Stage 10's and BlockMesh's dialogs take, and meaning the same
        // thing, so that a comparison at one target is a comparison of the
        // layouts.
        double targetEdgeLength = 0.05;
        // Floor and (optional) ceiling on the edges per chord; 0 = no ceiling.
        // Two is what a block needs before its interior has a node at all, but
        // one is enough for the mesh to exist and a layout finer than the
        // target is a real case, so the floor is one.
        int minIntervals = 1;
        int maxIntervals = 0;
        // Winslow sweeps over the interior nodes of the blocks whose worst
        // element is below smoothingThreshold (0: only the folded ones), the
        // boundary held. Stage 10's defaults and Stage 10's reason for them:
        // holding the boundary is what keeps the mesh conforming while the
        // interior is straightened, since the smoothing never moves a node
        // another block can see.
        int smoothingPasses = 500;
        double smoothingTolerance = 1e-4;   // of the target edge length
        double smoothingThreshold = 0.0;
        // Two vertices closer than this fraction of the model's diagonal are a
        // crack.
        double crackTolerance = 1e-9;
    };

    // One class of the "opposite sides carry the same count" relation: the
    // macro edges one chord cuts across, and the single count they all take.
    struct Chord {
        std::vector<int> edges;         // BlockDecomposition::edges indices
        int intervals = 0;
        double idealIntervals = 0.0;    // N*, before rounding
        double minLength = 0.0, maxLength = 0.0;   // of its macro edges
        double minSpan = 0.0, maxSpan = 0.0;       // of the isolines it was chosen from
        bool clamped = false;           // the floor or the ceiling, not F(N)
    };

    // One meshed block: a structured (ns+1) x (nt+1) grid of vertex indices,
    // row-major, vert[j * (ns+1) + i] at parameters (i/ns, j/nt) of the block.
    // The same layout as Stage 10's QuadMesh::Block and BlockMesh::Block, so
    // one drawing routine serves all three.
    struct Block {
        int block = -1;                 // BlockDecomposition::blocks index
        int ns = 0, nt = 0;
        std::vector<int> vert;
        int material = 0;
    };

    struct Report {
        double modelExtent = 1.0;       // diagonal of the decomposition's bounding box
        double target = 0.0;

        int chords = 0;
        int edgesAssigned = 0;          // macro edges
        int minIntervals = 0, maxIntervals = 0;
        int clampedChords = 0;
        double meanIntervals = 0.0;

        int blocks = 0;                 // blocks meshed
        int unmeshedBlocks = 0;         // a side with no macro edge behind it

        int vertices = 0;
        int quads = 0;

        double minEdge = 0.0, maxEdge = 0.0, meanEdge = 0.0;
        double edgeRatioRms = 0.0;      // rms of log(length / target) over the edges
        double worstEdgeRatio = 1.0;

        double minScaledJacobian = 0.0, meanScaledJacobian = 0.0;
        double minScaledJacobianBefore = 0.0;   // before the Winslow pass
        int invertedBefore = 0;                 // ... and its inverted count
        int smoothedBlocks = 0;
        int smoothingSweeps = 0;
        int reflexCorners = 0;          // block corners the mesh turns through more than pi at
        int invertedQuads = 0;          // a corner determinant that is not positive
        double minQuadArea = 0.0, maxQuadArea = 0.0;
        double meshArea = 0.0;

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

    // The decomposition is read and not kept: everything this needs is copied
    // out during construction, so the mesh outlives whatever produced it.
    BlockQuadMesh(const BlockDecomposition &decomposition, const Options &opts);

    // The three mesh::QuadMesh::from() looks for, so that
    //
    //     mesh::QuadMesh out = mesh::QuadMesh::from(blockQuadMesh);
    //
    // hands this to TMOP exactly as the other two methods' meshes are handed
    // to it.
    const std::vector<Point> &vertices() const { return verts_; }
    const std::vector<std::array<int, 4>> &quads() const { return cells_; }
    const std::vector<int> &quadMaterials() const { return cellMaterial_; }

    const std::vector<Block> &blocks() const { return grids_; }
    const std::vector<Chord> &chords() const { return chords_; }
    // Per macro edge: its count, and its nodes from `from` to `to`.
    const std::vector<int> &edgeIntervals() const { return intervals_; }
    const std::vector<std::vector<int>> &edgeVertices() const { return edgeNodes_; }
    const Report &getReport() const { return report_; }
    const Options &getOptions() const { return opts_; }

    // The vertices on macro edges with a block on one side only that are not
    // on dS: the rim of a hole, where the producer left a face out of the
    // decomposition (a T-junction on its side, say) and there are no elements.
    // That rim is not a curve of the model, but it is where the mesh of
    // whatever eventually fills the hole will have to meet this one, so a
    // caller may want to pin it under a smoother (mesh::QuadMesh::pinVertex).
    // Whether that is worth its cost is the caller's to measure; TraceMesh's
    // --pin-open-sides records what it measured.
    const std::vector<int> &openSideVertices() const { return openVerts_; }

    // An .obj of quadrilateral faces, `usemtl mat<id>` per material as
    // mesh::QuadMesh reads it back.
    bool writeOBJ(const std::string &path) const;

private:
    void assignIntervals(const BlockDecomposition &D);
    void meshEdges(const BlockDecomposition &D);
    void meshBlocks(const BlockDecomposition &D);
    void smooth();
    void check(const BlockDecomposition &D);

    Options opts_;
    Report report_;

    std::vector<int> chordOf_;                  // per macro edge
    std::vector<Chord> chords_;
    std::vector<int> intervals_;                // per macro edge
    std::vector<std::vector<int>> edgeNodes_;   // per macro edge, intervals + 1

    std::vector<int> openVerts_;
    std::vector<Point> verts_;
    std::vector<std::array<int, 4>> cells_;
    std::vector<int> cellMaterial_;
    std::vector<Block> grids_;
};

#endif // __BLOCK_QUAD_MESH_HXX__
