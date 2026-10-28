#ifndef __SKELETON_BASE_COMPLEX_HXX__
#define __SKELETON_BASE_COMPLEX_HXX__

#include <array>
#include <vector>

#include "mesh/BlockDecomposition.hxx"
#include "mesh/Mesh.hxx"

// The base complex of a quadrilateral mesh, as the shared BlockDecomposition.
//
// docs/medial_axis.md ends where the paper does, at Step 10: a quad mesh at the
// target size h, with the singularities Cheny et al. guarantee and nothing in
// between. Every other method in this codebase hands back a coarse block
// decomposition instead, and is meshed from it at a target edge length chosen
// separately (mesh/BlockQuadMesh) and smoothed by mesh::TMOP. The two meet at
// the base complex of the Step 10 mesh -- the coarsest block decomposition that
// mesh is a refinement of (Bommes et al., "Quad-mesh generation and
// processing: a survey", Sec. 2.2):
//
//   * its macrovertices are the mesh's singular vertices (valence other than 4
//     inside, other than 3 on the boundary) and the model's corners, which G3
//     asks of every layout;
//   * its sides are the grid lines out of them -- separatrices traced straight
//     through every regular vertex until they reach the boundary or another
//     singular vertex -- with every crossing of two of them a macrovertex too;
//   * its blocks are the patches those lines cut the mesh into, each a
//     structured a x b grid of the mesh's quads, so each is four-sided by
//     construction.
//
// So the blocks carry exactly the singularities Step 10 placed, and meshing
// them at the Step 10 target gives back the Step 10 mesh up to the smoothing.
//
// ### Patches that are not blocks
//
// A patch with no singular vertex or corner on one of its loops -- the ring of
// an annulus, whose Step 10 mesh is perfectly regular -- is not a disk, and
// the separatrices never cut it. Such a patch is cut where its boundary is
// furthest from any macrovertex, by the grid line through that point, and the
// patches are read again, up to Options::maxCutRounds times: the annulus takes
// two cuts and comes out as two blocks. One cut alone would leave a single
// block whose two ends are the same side, which no mesher here handles. A
// patch still not a block after that is left out of the decomposition and
// counted, as ZIPLINE's LayoutBlocks leaves a face with a T-junction out: its
// sides become macro edges with a block on one side only, and Report::coverage
// says how much of the mesh that cost.
//
// ### Geometry
//
// An interior side is the polyline of mesh vertices the separatrix ran
// through. A side on dS is laid on the model's own boundary polyline between
// its two macrovertices instead, as every other adapter lays its boundary
// sides, so that the blocks have exactly the model's outline and not the
// chords between boundary nodes.
class BaseComplex {
public:
    struct Options {
        // Two points closer than this fraction of the mesh's bounding-box
        // diagonal are one place: a model corner on a mesh vertex, a mesh
        // vertex on dS.
        double tolerance = 1e-4;
        // A vertex on the mesh's boundary this close to dS, as a fraction of
        // the diagonal, is moved onto it before anything is read: every
        // vertex of a mesh's boundary is a point of dS, and a block side on dS
        // runs between two of them. SKELETON sets it from its target size h,
        // since a Step 10 boundary vertex can be off dS for two reasons, both
        // bounded by h: Step 4 moves boundary vertices by up to delta_max =
        // 0.05 h and a T-mesh node on a moved vertex keeps the moved position
        // (1.4e-3 of the diagonal on geom031), and a boundary edge quantized
        // to zero welds its two ends, less than about h / 2 apart, into one
        // vertex at their mean, which is off dS by the sagitta of the piece
        // it cut off (3.2e-2 on geom012).
        double snapTolerance = 1e-2;
        // Rounds of cuts through patches that are not blocks.
        int maxCutRounds = 16;
    };

    struct Report {
        int vertices = 0, quads = 0;
        int interiorSingular = 0;     // interior vertices of valence != 4
        int boundarySingular = 0;     // boundary vertices of valence != 3
        int nonManifoldEdges = 0;     // mesh edges with more than two quads
        int corners = 0;              // model corners found on a mesh vertex
        int missedCorners = 0;        // model corners with no mesh vertex on them
        int separatrices = 0;         // grid lines traced
        int cuts = 0;                 // of which: cuts through patches that were not blocks
        int cutRounds = 0;
        int patches = 0;
        int blocks = 0;
        int badPatches = 0;           // patches left out of the decomposition
        int snappedVertices = 0;      // boundary vertices moved onto dS
        double maxSnap = 0.0;         // ... the furthest, as a fraction of the diagonal
        int boundarySides = 0;        // macro edges on dS ...
        int chordSides = 0;           // ... of which not found on dS, left as the mesh's chords
        // Interior sides that cross dS somewhere between their corners: a
        // side cut across a hole or out of the model, which no sound
        // decomposition has. A mesh whose boundary lies on dS and that is not
        // tangled has none but where dS bulges between two of its boundary
        // vertices further than its rows lie from it -- a hole of a few
        // vertices, a concave corner the mesh rounds off.
        int sidesCrossingBoundary = 0;
        double coverage = 0.0;        // area of the blocks / area of the mesh
    };

    // `quads` counter-clockwise on `vertices`. `loops` is dS as closed
    // polylines (first point not repeated), in either orientation; `corners`
    // are the points that must be macrovertices.
    BaseComplex(const std::vector<Point> &vertices, const std::vector<std::array<int, 4>> &quads,
                const std::vector<std::vector<Point>> &loops, const std::vector<Point> &corners,
                const Options &opts);
    BaseComplex(const std::vector<Point> &vertices, const std::vector<std::array<int, 4>> &quads,
                const std::vector<std::vector<Point>> &loops, const std::vector<Point> &corners);

    const BlockDecomposition &decomposition() const { return decomposition_; }
    const Report &getReport() const { return report_; }

    // The mesh the complex was read from, and which of its quads lie in a
    // block: the rest is the part of the model the decomposition leaves out.
    const std::vector<Point> &vertices() const { return verts_; }
    const std::vector<std::array<int, 4>> &quads() const { return quads_; }
    const std::vector<char> &quadInBlock() const { return quadInBlock_; }
    // The singular vertices, for drawing.
    std::vector<Point> singularPoints() const;

    // Each block as the structured grid of mesh vertices it is: (ns + 1) x
    // (nt + 1) of them, row-major, vert[j * (ns + 1) + i] at parameters
    // (i / ns, j / nt) of the block -- BlockQuadMesh::Block's layout, so a
    // block's grid here and its mesh there are indexed alike. A block whose
    // patch did not read back as that grid has none.
    struct Grid {
        int block = -1;               // BlockDecomposition::blocks index
        int ns = 0, nt = 0;
        std::vector<int> vert;
    };
    const std::vector<Grid> &grids() const { return grids_; }

private:
    void buildTopology();
    void snapToBoundary();
    struct Hit { int loop = -1; double s = 0.0, d = 1e300; Point p{0.0, 0.0}; };
    // The closest point of dS to p: its loop, its arc length along it, how far.
    Hit locate(const Point &p) const;
    void findSeeds(const std::vector<Point> &corners);
    // One separatrix from `from` along mesh edge `e`; false if `e` was traced.
    bool trace(int from, int e);
    void traceFrom(int v);
    void updateWalls();
    // The edge at v that continues u -> v straight on: at a regular interior
    // vertex, the one of its four edges that shares no quad with u-v.
    int straightOn(int u, int v) const;
    int edgeBetween(int a, int b) const;
    // Reads the patches; true when every one is a block. Fills cuts_ with one
    // seed per patch that is not.
    bool readPatches();
    void assemble();
    std::vector<Point> boundarySide(const std::vector<int> &seq, bool &exact) const;
    void countCrossings();

    Options opts_;
    Report report_;
    double diagonal_ = 1.0;

    std::vector<Point> verts_;
    std::vector<std::array<int, 4>> quads_;
    std::vector<std::vector<Point>> loops_;
    std::vector<std::vector<double>> loopArc_;   // cumulative length, size + 1 per loop

    std::vector<std::array<int, 2>> edgeVerts_;
    std::vector<std::vector<int>> edgeQuads_;
    std::vector<std::array<int, 4>> quadEdges_;  // edge k runs quads_[q][k] -> [k + 1]
    std::vector<std::vector<int>> vertEdges_;
    std::vector<std::vector<int>> vertQuads_;

    std::vector<char> onBoundary_, singular_, seed_, macro_;
    std::vector<char> marked_, wall_;            // per edge

    struct Patch {
        std::vector<int> quads;
        bool block = false;
        std::vector<int> cycle;    // CCW, starting at a corner, when a block
        std::vector<char> corner;  // parallel to cycle
    };
    std::vector<Patch> patches_;
    std::vector<int> cuts_;                      // seeds chosen by readPatches()

    // The patch's quads row by row from side 0; false if they are not the
    // grid its four sides say.
    bool gridOf(const Patch &P, const std::vector<int> &at, Grid &out) const;

    std::vector<char> quadInBlock_;
    std::vector<Grid> grids_;
    BlockDecomposition decomposition_;
};

#endif // __SKELETON_BASE_COMPLEX_HXX__
