#pragma once

#include "viewer/ViewerTypes.hxx"
#include "viewer/GL.hxx"

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "Parameterization/CutMesh.hxx"
#include "Parameterization/MIQ.hxx"
#include "Parameterization/UVGParam.hxx"
#include "polyvector/PolyVectors.hxx"
#include "crossfield/CrossField.hxx"
#include "sipg/SIPG.hxx"
#include "medialaxis/MedialAxis.hxx"
#include "medialaxis/MedialAxisTMesh.hxx"
#include "tracing/QuadLayout.hxx"
#include "UMBER/MotorcycleGraph.hxx"
#include "UMBER/Polysquare.hxx"

namespace viewer {

// Simple on-screen console for displaying log messages
class Console {
public:
    // Add a message to the console
    void log(const std::string &msg);

    // Clear all messages
    void clear();

    // Draw the console at the top of the screen.
    // fbw/fbh are the physical framebuffer dimensions (widget size × devicePixelRatio).
    void draw(int fbw, int fbh, float startY = 50.0f) const;

    // Set maximum number of lines to display (default 10)
    void setMaxLines(int n) { maxLines_ = n; }

private:
    std::vector<std::string> lines_;
    int maxLines_ = 10;
};

void drawMesh(const Mesh &m);

// Draw the mesh wireframe as a translucent overlay, for laying edges over
// filled geometry (a scalar field) without burying it. Relies on the blending
// enabled in initializeGL.
void drawMeshOverlay(const Mesh &m, float r, float g, float b, float a, float lineWidth);

void drawEdgeSetOnMesh(const Mesh &m,
                       const std::unordered_set<CutMesh::EdgeKey, CutMesh::EdgeKeyHash> &edges,
                       float r, float g, float b,
                       float lineWidth);

void drawArrow(const Point &p, const Point &dir, double scale, float r, float g, float b);

void drawField(const Mesh &m, const PolyField &field, double scale);

// Draw only the U field (single direction per triangle) from a CutMesh
void drawUField(const Mesh &m, const std::vector<Point> &uField, double scale);

// Draw only the V field (single direction per triangle) from a CutMesh
void drawVField(const Mesh &m, const std::vector<Point> &vField, double scale);

// Draw crossfield on mesh vertices from CrossField u_k_prev (MBO method)
// Each cross direction is (u_k_prev[i])^(1/4)
void drawVertexCrossField(const Mesh &m, const CrossField &cf, double scale);

// Draw crossfield on mesh vertices from CrossField u_k (MBO method, during stepping)
// Each cross direction is (u_k[i])^(1/4)
void drawVertexCrossFieldUK(const Mesh &m, const CrossField &cf, double scale);

// Draw SIPG p=0 cross field: one cross per triangle centroid from SIPG u_k.
void drawTriangleCrossField(const Mesh &m, const SIPG &sipg, double scale);

void drawDisk3D(const Point &center, double radius, float baseR, float baseG, float baseB, int segments = 96);

// Draw simple text overlay in screen coordinates (top-left origin).
// fbw/fbh are the physical framebuffer dimensions.
void drawTextOverlay(int fbw, int fbh, const char *text, float x, float y, float r, float g, float b);

// Compute view bounds for a UV mesh (for initializing the view state once).
void computeUVMeshBounds(const MIQSolver &miq, double &cx, double &cy, double &baseW, double &baseH);

// Draw UV mesh from MIQ parametrization (2D view of UV coordinates).
// Does not modify the view state - call computeUVMeshBounds first to set up the view.
void drawUVMesh(const MIQSolver &miq);

// Draw singularities on the UV mesh using the same coloring as the 3D view.
// Uses the origToCutVerts mapping from CutMesh to find UV coordinates.
void drawSingularitiesOnUV(const MIQSolver &miq, const CutMesh &cutMesh, 
                           const PolyField &field, double radius);

// Draw the medial axis (Voronoi edges + vertices) on top of a mesh.
void drawMedialAxis(const MedialAxis &ma, double vertexRadius);

// The coarse block decomposition of the medial axis: every zone filled with a
// translucent tint of its class colour so the blocking reads at a glance, the
// zone walls (boundary run and the two spokes) over the fill, the downsampled
// axis chains in full class colour on top, and the kept medial vertices as
// block corners. `cornerRadius` sizes the corner disks.
void drawMedialTMesh(const MedialAxisTMesh &tm, double cornerRadius);

// Draw X/Y coordinate axes through the origin, with tick marks spaced at a
// "nice" round interval, so the mesh's scale stays readable no matter what
// mode or phase is on screen. `vs` is the ViewState whose ortho is currently
// bound (the panel's world box sizes the axis extent and tick spacing).
void drawAxis(const ViewState &vs);

// Draw only the boundary edges of a mesh.
void drawBoundaryEdges(const Mesh &m);

// Fill the mesh with a per-vertex scalar field using a diverging blue/red ramp
// with a neutral midpoint. The field is a signed quantity oscillating about
// zero (a quasi-eigenfunction), so the range is taken symmetric: `vmax` should
// be max|f| and zero always lands on the neutral midpoint.
void drawScalarField(const Mesh &m, const Eigen::VectorXd &f, double vmax);

// Screen-space legend for drawScalarField, drawn bottom-left with the field
// extents labelled. fbw/fbh are the physical framebuffer dimensions.
void drawScalarFieldLegend(int fbw, int fbh, double vmin, double vmax, const char *title);

// Compute view bounds for a UVGParam parametrization.
void computeUVGParamBounds(const UVGParam &uvp, double &cx, double &cy, double &baseW, double &baseH);

// Draw UV mesh from UVGParam parametrization (2D view of UV coordinates).
// Does not modify the view state - call computeUVGParamBounds first to set up the view.
void drawUVGParam(const UVGParam &uvp);
void drawFlippedUVTriangles(const UVGParam &uvp);

// Draw SIPG singularities on the UVGParam view.
// singularVertices: (original-mesh vertex index, cross-index) pairs from SIPG.
void drawSingularitiesOnUVG(const UVGParam &uvp,
                             const std::vector<std::pair<int, double>> &singularVertices,
                             double radius);

// ── UMBER polysquare, in the parameter domain ────────────────────────────────

// View bounds for the polysquare (call once, before drawing it).
void computePolysquareBounds(const Polysquare &ps, double &cx, double &cy,
                             double &baseW, double &baseH);

// The deformed mesh as a wireframe at its (u, v) coordinates.
void drawPolysquare(const Polysquare &ps);

// Fill any triangle whose image is inverted. Nothing drawn on a valid map.
void drawFlippedPolysquareTriangles(const Polysquare &ps);

// The two kinds of edge on the boundary of the cut mesh, which look alike in
// the parameter domain and mean opposite things: the boundary of the model,
// which is what the axis alignment applies to, and the banks of the cuts,
// which are interior seams identified by a transition and are free to wander.
void drawPolysquareStructure(const Polysquare &ps, const HarmonicCut &hc);

// The frame field's boundary corners at their images, coloured as on the mesh.
// Every copy of a corner vertex is drawn, so a corner at the end of a cut
// appears once on each bank.
void drawPolysquareCorners(const Polysquare &ps, const HarmonicCut &hc,
                           const std::vector<std::pair<int, int>> &corners, double radius);

// ── UMBER block structure ────────────────────────────────────────────────────

// The edges of the block decomposition: the iso-lines, which are straight in
// the parameter domain and curved on the model. `parameterDomain` picks which
// of the two each segment is drawn in; the segments carry both, so the same
// call serves either half of a split screen.
void drawBlockEdges(const MotorcycleGraph &mg, bool parameterDomain, float lineWidth);

// The boundary of the model is an edge of the block decomposition too -- the
// blocks along it are closed by it -- so it is drawn in the same colour and
// weight as the traced lines. The banks of the cuts are skipped: those are
// interior seams that the blocks run straight through.
void drawBlockBoundary(const Polysquare &ps, const HarmonicCut &hc,
                       bool parameterDomain, float lineWidth);

// The nodes where those edges meet, coloured by what kind of node it is:
// yellow for a corner of the polysquare, green where a line leaves the model,
// cyan where two lines cross.
void drawBlockNodes(const MotorcycleGraph &mg, bool parameterDomain, double radius);

// The arcs of a quad layout: the separatrices cut at the nodes they meet, plus
// the pieces of the boundary that close the components. They are polylines
// following the streamlines of the cross field, so they are drawn as they are
// rather than as chords between their nodes -- the curvature between two nodes
// is the shape of the component's side, and straightening it would show a
// different layout from the one that was built.
void drawQuadLayoutArcs(const QuadLayout &layout, float lineWidth, float r, float g, float b);

// The nodes of a quad layout -- the corners of its blocks -- coloured by what
// kind of node each is: yellow for an irregular node of the field, orange for
// a corner of the model, green where a separatrix ran into the boundary
// squarely, cyan where two separatrices cross, magenta where two were joined
// head-on, and blue for a T-junction.
void drawQuadLayoutNodes(const QuadLayout &layout, double radius);

} // namespace viewer

