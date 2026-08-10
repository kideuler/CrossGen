#ifndef _MEDIAL_AXIS_HXX_
#define _MEDIAL_AXIS_HXX_

#include "mesh/Mesh.hxx"

#include <array>
#include <cmath>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// ─── Discrete medial axis of a planar domain ────────────────────────────────
//
// The axis is the subset of the Voronoi diagram of the boundary samples that
// lies inside the domain. It is represented through its dual: a polygonal
// complex over the boundary vertices (a constrained Delaunay triangulation to
// begin with). Each active DualCell carries one medial vertex at its
// circumcenter; each complex edge shared by two active cells carries one
// medial edge.
//
// The dual is polygonal rather than strictly triangular on purpose: merging
// two cells is the combinatorial half of the medial-axis simplification, and a
// merged cell is in general a k-gon whose vertices are (exactly or
// approximately) cocircular. Cells therefore store an ordered CCW vertex ring
// instead of a triangle index, and `dualIsTriangle()` distinguishes the two
// cases.

// Merge tolerance for deduplicateMedialVertices, relative to the local medial
// radius.
const double MEDIAL_VERTEX_MERGE_TOLERANCE = 0.01;

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

    // Circumcenter and circumradius: the position and medial radius of the
    // medial vertex dual to this cell.
    //
    // After a merge the ring is only approximately cocircular, because the
    // survivor's circle is kept as-is and the absorbed cell's vertices are not
    // projected onto it. `maxRadiusDeviation` records how far the ring has
    // drifted, as a fraction of `radius`; it is 0 for an unmerged triangle.
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
    Point coord{0.0, 0.0};   // Voronoi vertex = circumcenter of the dual cell
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
    int mergedCells = 0;              // cells absorbed by deduplicateMedialVertices
    int mergesRejectedMultiEdge = 0;  // merges skipped because the pair shared >1 edge
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

    // Medial branches: maximal chains of degree-2 vertices delimited by
    // vertices of degree != 2, plus closed cycles of degree-2 vertices (which
    // arise around holes). polyLineIsCycle marks the branches that close on
    // themselves, i.e. whose first and last entry are the same vertex; that
    // covers both a pure degree-2 cycle and a branch that leaves a junction
    // and returns to it.
    std::vector<std::vector<int>> polyLines;
    std::vector<bool> polyLineIsCycle;

    // ── Boundary ──
    std::vector<BoundaryVertex> boundaryLinks;  // indexed by mesh vertex index
    // Interior angle of the domain at each boundary vertex, in (0, 2*pi).
    // Values > pi mark concave (reflex) corners. NaN for non-boundary vertices.
    std::vector<double> interiorAngle;

    MedialAxisStats stats;

    explicit MedialAxis(std::shared_ptr<Mesh> mesh);

    // Merge medial vertices that sit closer together than
    // `tolerance * localRadius`, leaves first, then chains, then junctions.
    // Merging two medial vertices merges their dual cells, so the complex
    // stays consistent and the surviving cell's ring records every boundary
    // vertex the merged circle touches.
    //
    // This is a placeholder for the boundary-perturbation simplification the
    // axis really wants: it does not move boundary vertices onto a common
    // circle, so `maxRadiusDeviation` grows and the merged cells are only
    // approximately cocircular. Calls buildAxis() when done.
    void deduplicateMedialVertices(double tolerance = MEDIAL_VERTEX_MERGE_TOLERANCE);

    // (Re)build medialVertices / medialEdges from the active complex. A
    // candidate medial edge is kept only if its Voronoi segment lies inside
    // the domain; rejected candidates are counted in stats.
    void buildAxis();

    // Fill polyLines / polyLineIsCycle from the current axis.
    void createPolylines();

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

private:
    void buildComplex();
    void rebuildBoundaryLinks();
    void computeInteriorAngles();

    // Constrained edges of the active complex, as position pairs, for the
    // inside/no-crossing test.
    std::vector<std::array<Point, 2>> boundarySegments() const;

    // Merge `absorbed` into `survivor` across their single shared ring edge.
    // Returns false (leaving both untouched) if they do not share exactly one
    // edge.
    bool mergeCells(int survivor, int absorbed);

    std::vector<int> vertexCellCount_;
    double bboxDiagonal_ = 1.0;
};

#endif // _MEDIAL_AXIS_HXX_
