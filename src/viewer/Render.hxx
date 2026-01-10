#pragma once

#include "viewer/ViewerTypes.hxx"
#include "viewer/GL.hxx"

#include <string>
#include <vector>
#include <unordered_set>

#include "polyvec/CutMesh.hxx"
#include "polyvec/IntegerGridMap.hxx"
#include "polyvec/PolyVectors.hxx"

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

void drawDisk3D(const Point &center, double radius, float baseR, float baseG, float baseB, int segments = 96);

// Draw the UV mesh from an IntegerGridMap (clears view and shows UV coordinates).
void drawUVMesh(const IntegerGridMap &igm, float r, float g, float b, float lineWidth = 1.25f);

// Draw simple text overlay in screen coordinates (top-left origin).
// Must be called with proper orthographic projection set up for screen space.
void drawTextOverlay(GLFWwindow *window, const char *text, float x, float y, float r, float g, float b);

} // namespace viewer
