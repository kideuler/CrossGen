#pragma once

#include "viewer/ViewerTypes.hxx"
#include "viewer/GL.hxx"

#include <string>
#include <vector>
#include <unordered_set>

#include "IGM/CutMesh.hxx"
#include "IGM/MIQ.hxx"
#include "polyvector/PolyVectors.hxx"
#include "crossfield/CrossField.hxx"
#include "sipg/SIPG.hxx"
#include "medialaxis/MedialAxis.hxx"

namespace viewer {

// Simple on-screen console for displaying log messages
class Console {
public:
    // Add a message to the console
    void log(const std::string &msg);

    // Clear all messages
    void clear();

    // Draw the console at the top of the screen
    void draw(GLFWwindow *window, float startY = 50.0f) const;

    // Set maximum number of lines to display (default 10)
    void setMaxLines(int n) { maxLines_ = n; }

private:
    std::vector<std::string> lines_;
    int maxLines_ = 10;
};

void drawMesh(const Mesh &m);

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
// Must be called with proper orthographic projection set up for screen space.
void drawTextOverlay(GLFWwindow *window, const char *text, float x, float y, float r, float g, float b);

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

// Draw only the boundary edges of a mesh.
void drawBoundaryEdges(const Mesh &m);

} // namespace viewer

