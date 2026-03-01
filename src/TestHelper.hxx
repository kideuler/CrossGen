#ifndef __TEST_HELPER_HXX__
#define __TEST_HELPER_HXX__

#include "mesh/Mesh.hxx"
#include "triangle/TriangleMesher.hpp"
#include <cmath>
#include <memory>
#include <vector>

namespace TestHelper {

/// Create a triangulated rectangle [xmin,xmax] x [ymin,ymax] with target
/// average edge length h.  Returns a shared_ptr<Mesh> ready for use.
inline std::shared_ptr<Mesh> createBox(
    double xmin, double ymin,
    double xmax, double ymax,
    double h)
{
    using namespace triangle_wrapper;

    // Number of boundary points on each side (at least 2 per side)
    int nx = std::max(2, static_cast<int>(std::round((xmax - xmin) / h)) + 1);
    int ny = std::max(2, static_cast<int>(std::round((ymax - ymin) / h)) + 1);

    double dx = (xmax - xmin) / (nx - 1);
    double dy = (ymax - ymin) / (ny - 1);

    TriangleMesher2D::MeshInput input;

    // Build boundary vertices CCW: bottom -> right -> top -> left
    // Bottom edge (left to right)
    for (int i = 0; i < nx - 1; ++i)
        input.vertlist.push_back({xmin + i * dx, ymin});
    // Right edge (bottom to top)
    for (int i = 0; i < ny - 1; ++i)
        input.vertlist.push_back({xmax, ymin + i * dy});
    // Top edge (right to left)
    for (int i = 0; i < nx - 1; ++i)
        input.vertlist.push_back({xmax - i * dx, ymax});
    // Left edge (top to bottom)
    for (int i = 0; i < ny - 1; ++i)
        input.vertlist.push_back({xmin, ymax - i * dy});

    int npts = static_cast<int>(input.vertlist.size());

    // Build segment loop connecting consecutive boundary vertices
    std::vector<std::array<int, 2>> segments;
    for (int i = 0; i < npts; ++i)
        segments.push_back({i, (i + 1) % npts});

    input.segment_loops.push_back(segments);
    input.type.push_back(0); // exterior
    input.h = h;

    // Triangulate
    TriangleMesher2D::Options opts;
    opts.min_angle_degrees = 20.0;
    TriangleMesher2D mesher(opts);
    auto result = mesher.triangulate(input);

    // Convert to Mesh
    std::vector<Point> verts(result.verts.size());
    for (size_t i = 0; i < result.verts.size(); ++i)
        verts[i] = {result.verts[i][0], result.verts[i][1]};

    std::vector<Triangle> tris(result.triangles.size());
    for (size_t i = 0; i < result.triangles.size(); ++i)
        tris[i] = {result.triangles[i][0], result.triangles[i][1], result.triangles[i][2]};

    return std::make_shared<Mesh>(verts, tris);
}

/// Create a triangulated circle with given center, radius, and target
/// average edge length h.  Returns a shared_ptr<Mesh> ready for use.
inline std::shared_ptr<Mesh> createCircle(
    double cx, double cy,
    double radius,
    double h)
{
    using namespace triangle_wrapper;

    // Number of boundary points (circumference / h, at least 12)
    int npts = std::max(12, static_cast<int>(std::round(2.0 * M_PI * radius / h)));

    TriangleMesher2D::MeshInput input;

    double dtheta = 2.0 * M_PI / npts;
    for (int i = 0; i < npts; ++i) {
        double theta = i * dtheta;
        input.vertlist.push_back({cx + radius * std::cos(theta),
                                  cy + radius * std::sin(theta)});
    }

    // Segment loop
    std::vector<std::array<int, 2>> segments;
    for (int i = 0; i < npts; ++i)
        segments.push_back({i, (i + 1) % npts});

    input.segment_loops.push_back(segments);
    input.type.push_back(0); // exterior
    input.h = h;

    // Triangulate
    TriangleMesher2D::Options opts;
    opts.min_angle_degrees = 20.0;
    TriangleMesher2D mesher(opts);
    auto result = mesher.triangulate(input);

    // Convert to Mesh
    std::vector<Point> verts(result.verts.size());
    for (size_t i = 0; i < result.verts.size(); ++i)
        verts[i] = {result.verts[i][0], result.verts[i][1]};

    std::vector<Triangle> tris(result.triangles.size());
    for (size_t i = 0; i < result.triangles.size(); ++i)
        tris[i] = {result.triangles[i][0], result.triangles[i][1], result.triangles[i][2]};

    return std::make_shared<Mesh>(verts, tris);
}

} // namespace TestHelper

#endif // __TEST_HELPER_HXX__
