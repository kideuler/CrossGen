#pragma once

#include "ViewerTypes.hxx"

#include "polyvec/Mesh.hxx"

namespace viewer {

Bounds computeBounds(const Mesh &m);

Point triangleCentroid(const Mesh &m, const Triangle &t);

double averageTriangleEdgeLength(const Mesh &m);

} // namespace viewer
