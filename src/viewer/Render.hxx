#pragma once

#include "viewer/ViewerTypes.hxx"
#include "viewer/GL.hxx"

#include <array>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "MERIDIAN/Arrangement.hxx"
#include "MERIDIAN/ConeCut.hxx"
#include "MERIDIAN/ConeSingularities.hxx"
#include "MERIDIAN/Immersion.hxx"
#include "MERIDIAN/Interfaces.hxx"
#include "MERIDIAN/RicciFlow.hxx"
#include "MERIDIAN/Separatrices.hxx"
#include "MERIDIAN/QuadMesh.hxx"
#include "MERIDIAN/SplineFit.hxx"
#include "MERIDIAN/SubdomainLabels.hxx"
#include "Parameterization/CutMesh.hxx"
#include "Parameterization/MIQ.hxx"
#include "Parameterization/UVGParam.hxx"
#include "polyvector/PolyVectors.hxx"
#include "crossfield/CrossField.hxx"
#include "sipg/SIPG.hxx"
#include "medialaxis/MedialAxis.hxx"
#include "medialaxis/MedialAxisTMesh.hxx"
#include "quantization/QuantTMeshConvert.hxx"
#include "tracing/QuadLayout.hxx"
#include "UMBER/MotorcycleGraph.hxx"
#include "UMBER/Polysquare.hxx"
#include "TORSION/FieldFrames.hxx"

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

// The quantized block decomposition: the integer edge lengths the QGP
// quantizer settled on, drawn as the quad grid they prescribe -- edges
// only, no fill.
//
// Every grid line is a transfinite (Coons) curve through the four sides of
// its block, evaluated against their real geometry rather than chords
// between the corners, so the cells follow the block's curvature and the
// outermost lines reproduce the curved sides exactly. Cell corners sit on
// the quantization ticks, so the grids of two faces meet flush along the
// edge they share exactly when the quantization is consistent -- a
// violated constraint shows up as a visible mismatch.
void drawQuantizedBlocks(const BlockQuant &bq);

// The same, for the T-mesh a QuadLayout converts to directly (see
// QuantTMeshConvert.hxx): the block decomposition tracing leaves after
// separatrix tracing and chord collapse, rather than the medial axis one.
//
// `vertexRadius` > 0 marks every grid vertex with a disk of that radius --
// block corners and interior cell corners alike, since after quantization
// every crossing of two grid lines is a vertex of the final decomposition.
void drawQuantizedLayout(const QuadLayoutQuant &lq, double vertexRadius = 0.0);

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

// ── MERIDIAN Stage 0b: the material interface network ────────────────────────
//
// On a multi-material model the interfaces are an input, not a result: they are
// the curves the tags already put on the triangulation, and the whole of the
// multi-material path is the pipeline being made to keep them. So they are
// drawn from the first MERIDIAN phase, before there is a field or a cone or a
// layout to draw them against, and they stay under every phase after it.
//
// Two pictures, because they answer two different questions.

// The triangles filled by material, translucently, so the regions read as
// regions rather than as a wireframe that happens to change colour. This is the
// domain the layout has to be compatible with, and on a model whose tags are
// wrong it is the picture that says so before anything else has run. Uniform on
// a single-material mesh, which is why the viewer only turns it on when there
// is more than one material.
void drawMaterialFill(const Mesh &m, float alpha);

// The network itself: every branch as the polyline of mesh edges it is, and
// every node as a disk coloured by what kind of node it is.
//
//   white     a junction -- three or more branches meet inside S
//   green     a landing -- a branch reaches dS
//   orange    a kink -- the interface turns by more than the threshold, and the
//             layout has to turn with it or an element straddles the corner
//   violet    a loop split -- a closed interface carries no node of its own, so
//             it is cut into arcs that the faces inside it can have corners at
//   cyan      a corner Stage 0b's balance() put on a smooth interface to make
//             the region's own Gauss-Bonnet count come out
//   red       a dangling branch: one interface edge ending nowhere, which is a
//             tag error and not a feature
//
// A node whose sectors are more than Options::wellPosedAngle from whole quarter
// turns carries a dark halo. That is not an error either -- an oblique junction
// is a property of the domain, and the corpus has several on purpose -- but it
// is where the layout has the most turning to absorb, so it is worth being able
// to find by eye.
void drawInterfaceNetwork(const Interfaces &itf, double nodeRadius, float lineWidth);

// Screen-space key for those node colours. `separatrixLegendShown` lifts it
// clear of the separatrix legend when both are on screen.
void drawInterfaceLegend(int fbw, int fbh, bool separatrixLegendShown);

// ── MERIDIAN: cones, the cutting graph, and the flat cone metric ─────────────
//
// Shepherd, Gu and Hughes (2022), Secs. 3.1 and 3.2. Everything here is drawn
// over the model at its own coordinates -- the metric Ricci flow produces is a
// set of edge lengths, not a set of positions, so there is no second geometry
// to draw until the immersion of Stage 4 exists. What can be shown is where the
// metric differs from the one the model came with, and what the cone angles it
// was driven to actually came out as; drawConeFans() is the second of those.

// The prescribed cones, one disk each, coloured by index: blue at +1 (a
// valence-three corner), red at -1 (valence five), cyan at -2 (valence six, the
// colour Fig. 11 of the paper uses for them), yellow for anything larger. An
// interior cone carries a pale halo, because it is the kind that needs an arc
// of the cutting graph run out to the boundary; a boundary cone is already
// there and gets none.
void drawCones(const Mesh &m, const ConeSingularities &cones, double radius);

// Screen-space key for those colours, bottom-left, above the scalar legend.
void drawConeLegend(int fbw, int fbh);

// The cutting graph G, with its two kinds of arc kept apart: magenta for the
// void arcs (HarmonicCut's, one per hole) and amber for the cone arcs (Sec.
// 3.2.2's, one per interior cone). The distinction is not cosmetic -- a void
// arc has both ends on the boundary and a cone arc has one end at a
// singularity, which is what makes it a leaf of G.
void drawCuttingGraph(const ConeCut &cut, float lineWidth);

// ── the flat cone metric ─────────────────────────────────────────────────────

// The metric prepared for drawing: one entry per edge of the *flow*
// triangulation (which is the thing that carries the metric, and which the
// weighted-Delaunay flipping may have changed), carrying how far the flow
// stretched that edge.
//
// The stretch is measured as log(l_flat / l_input) with its mean removed. Both
// halves of that matter. The log makes a halving and a doubling symmetric
// about the neutral midpoint of the ramp, and removing the mean throws away the
// global scale -- which is arbitrary, since the Ricci energy is invariant under
// u -> u + c and only one vertex being pinned fixes it at all. What is left is
// the part of the metric change that is real: where the surface had to be
// stretched relative to everywhere else to flatten it.
struct FlatMetric {
    struct Edge {
        int a = -1, b = -1;
        double t = 0.0;             // log(l_flat / l_input), mean removed
        bool newDiagonal = false;   // an edge the flipping introduced
    };
    std::vector<Edge> edges;
    std::vector<std::array<int, 2>> replaced; // input edges the flipping removed
    double absMax = 1.0;                      // max |t|, the symmetric ramp range
    double minRatio = 1.0, maxRatio = 1.0;    // exp of the extremes, for the console
};

FlatMetric buildFlatMetric(const RicciFlow &flow);

// The flow triangulation coloured by that stretch, with the flips called out:
// the input edges the flipping removed in dim slate underneath, the diagonals
// it put in their place in green on top. Everything else is the diverging ramp
// drawScalarField uses, so blue reads "shrunk" and red "stretched".
void drawFlatMetric(const Mesh &m, const FlatMetric &fm, float lineWidth);

// ── the cone angles, unfolded ────────────────────────────────────────────────
//
// The single thing Stage 3 is for is that the angle sum at a cone comes out at
// exactly 2pi - (pi/2) I(v), and nowhere else is it anything but 2pi. That is a
// statement about a metric, so no drawing of the model can show it -- the
// triangles on screen still carry the angles they came with. Laying one cone's
// one-ring out in the plane *under the new metric* can: walk the fan, opening
// each triangle by the angle the flat metric gives it, and the total is the
// cone angle.
//
// A cone of index +1 then closes 90 degrees early and leaves a visible wedge of
// gap; one of index -1 overshoots by 90 and laps itself; one of index -4 wraps
// twice round. The first spoke is drawn white and the last yellow, so the
// discrepancy between them is the cone angle's excess over 2pi, read straight
// off the picture. For an interior cone those two spokes are the *same* mesh
// edge, arrived at from both sides of the fan.
struct ConeFan {
    int vertex = -1;
    int index = 0;
    bool onBoundary = false;
    bool closed = false;                     // interior: the last spoke returns to the first
    double angleSum = 0.0;                   // Theta, the cone angle in the flat metric
    std::array<double, 2> center{{0.0, 0.0}}; // where it sits in the gallery
    std::vector<std::array<double, 2>> ring;  // unfolded one-ring, largest radius 1
};

// One fan per cone, laid out in a grid. Empty if the flow has no cones, or for
// any cone whose star could not be walked (a pinch, or a vertex the flipping
// left with a broken link).
std::vector<ConeFan> buildConeFans(const RicciFlow &flow, const ConeSingularities &cones);

// View bounds for the gallery (call once, before drawing it).
void computeConeFanBounds(const std::vector<ConeFan> &fans,
                          double &cx, double &cy, double &baseW, double &baseH);

void drawConeFans(const std::vector<ConeFan> &fans);

// ── TORSION: the combed field (Stage 3F) ─────────────────────────────────────
//
// The cross field again, but as one branch of it rather than as four
// indistinguishable arms: theta_hat_t = theta_t + (pi/2) a_t, drawn as the
// frame J*_t it defines -- X_t in warm colour, Y_t in cool -- so that the
// picture shows which arm the comb chose and not merely where the cross points.
//
// The colour is the integer a_t, cycled through a fixed palette. That is the
// whole of what the stage is judged on and it is judged by eye: a_t is constant
// over a patch and steps only where an arc of G is crossed, because Omega is a
// disk and the branch is path-independent on it. A step in the middle of a
// patch is the combing defect Report::combingDefects counts, seen at the place
// it happened rather than as a number.
//
// Faces outside Omega -- the ones the cut removed nothing of but the comb was
// never asked about -- carry a_t = 0 like the seed and are drawn as such; the
// count that says whether any such face exists is Report::unreachedFaces, and
// it is zero on every model the cut left a disk.
void drawCombedFrames(const Mesh &m, const FieldFrames &ff, double scale);

// Screen-space key for the branch palette above, under the cone legend.
void drawCombedFrameLegend(int fbw, int fbh, const FieldFrames &ff);

// ── MERIDIAN: the layout Psi on Omega ────────────────────────────────────────
//
// Stages 4 to 6 are the first ones that produce a *map*, so unlike everything
// above they have a second geometry of their own to draw: Omega laid out in the
// plane, at psi_R out of Stage 4 and at Psi out of Stage 6. Both are drawn by
// the same routine, since the difference between them is what the picture is
// for -- the same triangulation, the same seams and the same cones, moved.
//
// What each colour is answering, in the order Definition 2.1 asks it:
//
//   red fill      a triangle with det J < 0. Q1 has failed there, and E1's
//                 barrier cannot repair it, so one is a defect and not a
//                 blemish.
//   blue / green  the boundary edges of Gamma_u and Gamma_v out of Stage 5.
//                 Q3 says each is on a constant-coordinate line, so once the
//                 continuation has converged every blue run is vertical and
//                 every green one horizontal, and a run that is neither is
//                 exactly where the Q3 residual lives.
//   amber/magenta the two banks of each arc of G. Q4 says they are the same
//                 curve up to R_k, so they are drawn apart: an arc whose two
//                 banks are not congruent has not met E4.
//   white-cored   the feature chains, which on a multi-material model are the
//   blue / green  material interfaces. Same rule as dS -- a chain labelled u is
//                 vertical, one labelled v horizontal -- and the angle between
//                 two of them at a node of the network is what E6 holds at a
//                 whole number of right angles.
//   cone disks    the same index colours drawCones() uses on the model, at
//                 every child of the cone in Omega. An interior cone that the
//                 cut opened appears once per child.

// View bounds for the layout panel (call once, when the map lands).
void computeLayoutBounds(const std::vector<Point> &uv,
                         double &cx, double &cy, double &baseW, double &baseH);

// `labels` may be null, which drops the Gamma_u / Gamma_v colouring and leaves
// the boundary in grey -- what psi_R looks like before Stage 5 has run.
void drawLayoutUV(const Immersion &imm, const SubdomainLabels *labels,
                  const std::vector<Point> &uv, double coneRadius);

// ── MERIDIAN: the separatrices of Psi ────────────────────────────────────────
//
// Sec. 4. Stage 7 stores each curve as a list of (triangle, barycentric entry,
// barycentric exit), and barycentric coordinates are the same numbers in the
// image and on S, so the *same* curve can be drawn in either world without
// re-projecting it: Separatrices::polyline() evaluates the same steps against
// Psi or against the model. That is why both halves of this phase are drawn by
// one routine with a Space argument -- they are not two computations, they are
// one curve seen from its two ends, which is the whole point of the picture.
//
// The colour is the end the curve came to, because that is what Q5 asks:
//
//   green    terminated at a cone. The case Q5 wants, and the one that closes
//            a patch of the layout.
//   blue     left transversely through dS. Also allowed, and what happens to
//            the rays of a boundary cone that point out of the model.
//   red      still running when the step cap was reached. Not a failure of the
//            tracer: away from the cones these are geodesics of a flat cone
//            metric, so a direction no Gamma_topo constraint quantised does
//            not close. A red curve names the pair of cones whose connectivity
//            constraint Sec. 3.3 is missing.
//   orange   stuck or degenerate -- a ray that could not be marched at all.
//
// A curve is drawn broken in the image and unbroken on the model, and that
// asymmetry is Q4 rather than a defect: crossing an arc of G jumps the image to
// the other bank of the cut while the pullback walks straight on. `breaks` out
// of polyline() is what puts the gaps in.
//
// The terminus of anything that did not reach a cone gets a disk, since a green
// curve already ends on a cone disk and the other three are the ones worth
// finding.
void drawSeparatrices(const Separatrices &sep, Separatrices::Space space,
                      double endRadius, float lineWidth);

// Screen-space key for those four colours, under the cone legend.
void drawSeparatrixLegend(int fbw, int fbh);

// ── MERIDIAN: the blocks of Stages 8 and 9 ───────────────────────────────────
//
// The layout drawn on the model, which is the left half of the paper's Fig. 12:
// every side of every patch in light blue, every node of the arrangement as a
// green disk.
//
// Both live on S and only on S. Stage 8 is built there on purpose -- Psi
// overlaps itself, so an arrangement computed in the image would invent
// crossings -- and Stage 9's control points are fitted to the arcs' polylines
// on S, so nothing here needs projecting.
//
// `fit` may be null, and where it is non-null a side is drawn as the fitted
// spline sampled rather than as the polyline it was fitted to. The two differ
// by the fit's deviation, which is smaller than the line is wide; what makes
// the distinction worth keeping is that a side with no fit behind it is a side
// of a face Stage 9 skipped, and drawing the raw arc there is honest about it.
//
// Only faces flagged `patch` are outlined. The unbounded face and the holes are
// faces of the subdivision too, and they are not blocks.
void drawLayoutPatches(const Arrangement &arr, const SplineFit *fit,
                       double nodeRadius, float lineWidth, int samples = 16);

// MERIDIAN Stage 10: the quadrilateral mesh the interval assignment and the
// transfinite interpolation produced, drawn on S as edges only.
//
// A quad whose signed area is not positive is filled first, in red, before any
// edge goes down. A folded element is not visible from its wireframe -- the
// four edges of a bowtie look like the four edges of a quad -- and the whole
// reason to look at this picture rather than read the report is to see where
// the folds are, so they are given the one thing the rest of the mesh has not
// got: a fill.
//
// `blockLineWidth` > 0 draws the block walls over the interior edges in a
// second, heavier pass, which is what makes the structure legible: the interior
// of a block is a regular grid, and where two of them meet the rows either
// match or they do not.
// `materialFill` tints every element with the material of the region its
// centroid landed in, which is the one way to see the property the
// multi-material path exists for: an element that straddles an interface is one
// no analysis code can integrate, and with the network drawn over the fill a
// straddling element is a cell the interface runs through rather than along.
void drawQuadMesh(const QuadMesh &qm, float lineWidth, float blockLineWidth,
                  bool materialFill = false);

} // namespace viewer
