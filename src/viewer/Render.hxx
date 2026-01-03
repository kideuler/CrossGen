#pragma once

#include "viewer/ViewerTypes.hxx"

#include <unordered_set>

#include "polyvec/CutMesh.hxx"
#include "polyvec/PolyVectors.hxx"

namespace viewer {

void drawMesh(const Mesh &m);

void drawEdgeSetOnMesh(const Mesh &m,
                       const std::unordered_set<CutMesh::EdgeKey, CutMesh::EdgeKeyHash> &edges,
                       float r, float g, float b,
                       float lineWidth);

void drawArrow(const Point &p, const Point &dir, double scale, float r, float g, float b);

void drawField(const Mesh &m, const PolyField &field, double scale);

void drawDisk3D(const Point &center, double radius, float baseR, float baseG, float baseB, int segments = 96);

} // namespace viewer
