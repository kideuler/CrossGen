#ifndef __MEDIAL_AXIS_HXX__
#define __MEDIAL_AXIS_HXX__

#include "mesh/Mesh.hxx"

#include <array>
#include <cmath>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// ─── Discrete medial axis of a planar domain ────────────────────────────────
//
// Steps 3-5.1 of docs/medial_axis.md (Cheny et al., "Planar quad meshing with
// guarantees", Sec. 3.3-3.5). The axis is the subset of the Voronoi diagram
// of the boundary samples that lies inside the domain. It is represented
// through its dual: a polygonal complex M over the boundary vertices (a
// constrained Delaunay triangulation to begin with). Each active DualCell
// carries one medial vertex at the centre of its stored circle; each complex
// edge shared by two active cells carries one medial edge.
//
// The dual is polygonal rather than strictly triangular on purpose (Sec.
// 3.5): merging two cells is the combinatorial half of the medial-axis
// simplification, and a merged cell is in general a k-gon whose vertices are
// cocircular. Cells therefore store an ordered CCW vertex ring instead of a
// triangle index, and `dualIsTriangle()` distinguishes the two cases. From
// construction on, M is the single source of truth (Step 3.4): nothing is
// ever re-triangulated, and a cell's circle is stored, never recomputed from
// its ring.
//
// simplify() is Step 4. The axis is never edited directly: boundary
// vertices are moved onto a neighbouring circle (4.1) or removed from the
// loop (4.2), within a displacement budget delta_max of their original
// position, and the cells are merged to match, so the axis and the boundary
// map phi built on it stay consistent. That means the boundary this class
// describes after simplify() is not the input boundary any more: the input
// mesh's vertex positions are moved in place, and the original positions and
// loop order are kept in `originalPosition` / `originalNext`, which is what
// later stages draw block sides along (Step 1.6, Step 8.3).

// Boundary vertex connectivity. `next`/`prev` are oriented so that the domain
// interior lies on the LEFT of next: v -> boundaryLinks[v].next. The outer
// loop therefore runs CCW and hole loops run CW.
struct BoundaryVertex {
    int vertexIndex = -1;
    int nextBoundaryVertex = -1;
    int prevBoundaryVertex = -1;
};

// A cell of the dual polygonal complex, i.e. the Delaunay-side dual of one
// medial (Voronoi) vertex.
struct DualCell {
    // Boundary-vertex indices in CCW order. Three of them for a plain
    // Delaunay triangle, more once cells have been merged.
    std::vector<int> verts;

    // The stored circle (c_F, r_F): the position and medial radius of the
    // medial vertex dual to this cell. Set once from the triangle and kept
    // through every merge -- the survivor of a merge keeps its own circle,
    // and the absorbed cell's private vertices are projected onto it (Step
    // 4.1), so the merged ring is cocircular up to rounding.
    // `maxRadiusDeviation` measures that, as a fraction of `radius`; it is 0
    // for an unmerged triangle.
    Point  center{0.0, 0.0};
    double radius = 0.0;
    double maxRadiusDeviation = 0.0;

    // Whether `center` lies inside the domain. Recorded at construction time
    // from the original triangle (whose centroid is trivially interior) and
    // preserved across merges, since a merge never moves the survivor's
    // center.
    bool centerInside = true;

    // Cleared when the cell is absorbed by a merge.
    bool active = true;

    // Index into MedialAxis::medialVertices, or -1 before buildAxis().
    int medialVertex = -1;

    // Number of cells absorbed into this one.
    int timesFused = 0;

    bool dualIsTriangle() const { return verts.size() == 3; }
};

struct MedialVertex {
    Point coord{0.0, 0.0};   // Voronoi vertex = centre of the dual cell's circle
    double radius = 0.0;     // medial radius

    int cell = -1;           // index into MedialAxis::cells
    bool dualIsTriangle = true;

    // False when the circumcenter fell outside the domain. Such a vertex is
    // not on the medial axis at all: every edge incident to it is rejected by
    // the inside filter, leaving it isolated. It is still emitted so that the
    // cell <-> medial vertex correspondence stays a bijection, but consumers
    // should skip it.
    bool insideDomain = true;

    // Medial edge indices, in CCW order around `coord`. Edge k leaves through
    // ring edge (cell.verts[i], cell.verts[i+1]) for increasing i, so the
    // cyclic order matches the dual cell's ring order.
    std::vector<int> incidentEdges;

    // Adjacent medial vertex indices, derived from incidentEdges. Two cells
    // can share two ring edges, so this may hold fewer entries than `degree`;
    // use incidentEdges when the distinction matters.
    std::unordered_set<int> neighbors;

    int degree = 0;          // == incidentEdges.size()
};

// Counters describing how faithful the extracted axis is. Non-zero values in
// the first four fields mean the boundary sampling does not meet the
// eps-sampling assumption the construction relies on.
struct MedialAxisStats {
    int degenerateCells = 0;          // cells whose circumcircle could not be computed
    int centersOutsideDomain = 0;     // active cells whose circumcenter is outside
    int edgesCrossingBoundary = 0;    // candidate medial edges dropped by the inside filter
    int nonManifoldBoundaryVertices = 0; // boundary vertices shared by more than one loop
};

// ─── Step 4 parameters ──────────────────────────────────────────────────────
//
// delta_max and alpha_feat are the guide's [R] defaults (Sec. 1 table): the
// paper names the displacement budget but not its value, and locking feature
// corners is a reconstruction that keeps CAD corners out of the pruning.
struct SimplifyOptions {
    // h, the target element size. <= 0 takes MedialAxis::defaultTargetSize(),
    // the same default the T-mesh downsampling uses.
    double targetSize = 0.0;
    // delta_max: no boundary vertex ever ends up farther than this from its
    // original position. < 0 takes 0.05 * h.
    double maxDisplacement = -1.0;
    // alpha_feat: a boundary vertex whose turning angle exceeds this is a
    // locked feature corner, never moved or removed (Step 1.4).
    double featureAngle = M_PI / 6.0;
    // Two adjacent cells are one circle when every vertex of one lies within
    // this fraction of delta_max of the other's circle: the triangulation of
    // cocircular samples is ambiguous, and whichever choice was made is
    // merged back into one polygon (Step 2.3). Nothing is moved for it, so
    // the merged ring is cocircular to a hundredth of what Step 4 may move a
    // vertex -- enough to take in samples that a CAD circle or a nearly
    // symmetric shape put on one circle up to rounding, which otherwise
    // split one intersection point into two a few 1e-6 apart.
    double cocircularTolerance = 0.01;
    // 4.1 to a fixpoint, then 4.2, repeated while either changes anything
    // (the optional iteration of Step 4), at most this many times.
    int maxRounds = 8;
};

struct SimplifyStats {
    int cocircularMerges = 0;   // Step 2.3: coincident circles merged
    int endpointMerges = 0;     // Step 4.1: endpoint cells absorbed by their neighbour
    int endpointsKept = 0;      // endpoints 4.1 left alone at the fixpoint
    int flattenedChains = 0;    // Step 4.2: polar chains replaced by their chord
    int removedVertices = 0;    // boundary vertices 4.2 deleted from the loop
    int deletedTwoGons = 0;     // endpoint cells 4.2 reduced to a 2-gon and deleted
    int rounds = 0;
    double maxDisplacement = 0.0;  // largest |x - x0| over all boundary vertices
    double budget = 0.0;           // the delta_max that was applied
};

class MedialAxis {
public:
    std::shared_ptr<Mesh> mesh;

    // ── Dual polygonal complex ──
    std::vector<DualCell> cells;

    // Undirected complex edge -> the (up to two) active cells on either side.
    // A key with a single cell is a constrained/boundary edge. This map is the
    // single source of truth for adjacency; medial vertices and edges are
    // derived from it by buildAxis().
    std::unordered_map<MeshEdgeKey, std::array<int, 2>, MeshEdgeKeyHash> edgeCells;

    // ── Medial axis, rebuilt by buildAxis() ──
    std::vector<MedialVertex> medialVertices;
    std::vector<Edge> medialEdges;   // pairs of medialVertices indices
    // The complex edge each medial edge is dual to, i.e. the two boundary
    // samples (p_a, p_b) whose Voronoi cells the medial edge separates. A
    // point inserted on the medial edge maps to the boundary once through
    // each of their cells (Step 5.3, last remark).
    std::vector<Edge> medialEdgeDual;

    // Medial branches: maximal chains of degree-2 vertices delimited by
    // vertices of degree != 2, plus closed cycles of degree-2 vertices (which
    // arise around holes). polyLineIsCycle marks the branches that close on
    // themselves, i.e. whose first and last entry are the same vertex; that
    // covers both a pure degree-2 cycle and a branch that leaves a junction
    // and returns to it. polyLineEdges[b][i] is the medial edge between
    // polyLines[b][i] and polyLines[b][i + 1] -- two cells can share two
    // ring edges, so the vertex pair alone does not name it.
    std::vector<std::vector<int>> polyLines;
    std::vector<std::vector<int>> polyLineEdges;
    std::vector<bool> polyLineIsCycle;

    // ── Boundary ──
    std::vector<BoundaryVertex> boundaryLinks;  // indexed by mesh vertex index
    // Interior angle of the domain at each boundary vertex, in (0, 2*pi), at
    // the current positions. Values > pi mark concave (reflex) corners. NaN
    // for non-boundary vertices.
    std::vector<double> interiorAngle;

    // The boundary as it was before simplify(): x0 per mesh vertex, and the
    // loop successor of every boundary vertex (-1 elsewhere). A vertex 4.2
    // removed keeps both, so originalPolyline() can still walk through it.
    std::vector<Point> originalPosition;
    std::vector<int> originalNext;
    // Locked feature vertices (Step 1.4), from the original turning angles.
    std::vector<char> locked;

    MedialAxisStats stats;
    SimplifyStats simplifyStats;

    explicit MedialAxis(std::shared_ptr<Mesh> mesh);

    // Step 4. Moves and removes boundary vertices (in `mesh`, in place) and
    // merges cells to match, then rebuilds the boundary links, the angles and
    // the axis. Run it once, before building the boundary map.
    void simplify(const SimplifyOptions &opts = SimplifyOptions());

    // (Re)build medialVertices / medialEdges from the active complex. A
    // candidate medial edge is kept only if its Voronoi segment lies inside
    // the domain; rejected candidates are counted in stats.
    void buildAxis();

    // Fill polyLines / polyLineEdges / polyLineIsCycle from the current axis.
    void createPolylines();

    // The original boundary polyline from boundary vertex a forward to b
    // (both included), at the original positions. For two vertices adjacent
    // on the current loop this is the stretch of input boundary the current
    // edge (a, b) stands for, removed vertices and all.
    std::vector<Point> originalPolyline(int a, int b) const;

    // ── Complex queries, for building simplification passes on top ──

    // Ring edge i of `cell` is (verts[i], verts[(i+1) % n]).
    std::array<int, 2> ringEdge(int cell, int edgeIdx) const;

    // Active cell on the far side of ring edge i, or -1 if that edge is
    // constrained (on the domain boundary).
    int neighborAcross(int cell, int edgeIdx) const;

    // True when ring edge i lies on the domain boundary.
    bool isConstrained(int cell, int edgeIdx) const;

    // Number of ring edges of `cell` shared with another active cell. This is
    // the degree of the dual medial vertex before the inside filter is
    // applied; use MedialVertex::degree for the post-filter degree.
    int cellDegree(int cell) const;

    // Number of active cells incident to boundary vertex v. A vertex used by
    // exactly one cell can be moved or removed without disturbing the rest of
    // the complex.
    int vertexCellCount(int v) const;

    // True when the interior angle at boundary vertex v is a sharp convex
    // corner (below `threshold`).
    bool isSharpCorner(int v, double threshold = M_PI * 0.75) const;
    // True when the interior angle at boundary vertex v exceeds pi.
    bool isConcaveCorner(int v) const;

    // Bounding-box diagonal of the domain; the natural length scale for
    // displacement thresholds.
    double boundingBoxDiagonal() const { return bboxDiagonal_; }

    // The h every stage falls back to when given none: a tenth of the
    // bounding-box diagonal.
    double defaultTargetSize() const { return 0.1 * bboxDiagonal_; }

private:
    void buildComplex();
    void rebuildBoundaryLinks();
    void computeInteriorAngles();

    // Constrained edges of the active complex, as position pairs, for the
    // inside/no-crossing test.
    std::vector<std::array<Point, 2>> boundarySegments() const;

    // Merge `absorbed` into `survivor` across their single shared ring edge.
    // Returns false (leaving both untouched) if they do not share exactly one
    // edge or their union would pinch.
    bool mergeCells(int survivor, int absorbed);

    // Whether segment ab properly crosses a current boundary edge, ignoring
    // the edges incident to any vertex in `ignore`.
    bool crossesBoundary(const Point &a, const Point &b,
                         const std::unordered_set<int> &ignore) const;

    // The three passes of simplify(); each returns how much it changed.
    int mergeCocircular(double tolerance, std::vector<int> &regionCells);
    int removeEndpoints(double budget, std::vector<int> &regionCells);
    int flattenArcs(double budget, std::vector<int> &regionCells);

    std::vector<int> vertexCellCount_;
    std::vector<int> regionOf_;   // cell -> connected component of the complex
    double bboxDiagonal_ = 1.0;
};

#endif // __MEDIAL_AXIS_HXX__
