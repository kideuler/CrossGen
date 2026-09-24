#ifndef __TRIANGLE_LOCATOR_HXX__
#define __TRIANGLE_LOCATOR_HXX__

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "mesh/Mesh.hxx"

// ---------------------------------------------------------------------------
// Which triangle of a planar mesh a point is in, by a uniform grid of about one
// triangle per cell -- the size at which a uniform grid's build cost and its
// query cost cross over. The tracing needs the answer in two places that have
// nothing else in common: LayoutBlocks, asking which material a block is in,
// and StemExtension, asking which triangle a T-junction sits in before tracing
// on from it.
//
// (UMBER's BlockLayout carries its own copy of the same structure, local to its
// .cxx; it predates this one and is left alone.)
// ---------------------------------------------------------------------------
class TriangleLocator {
public:
    explicit TriangleLocator(const Mesh &m) : mesh(&m) {
        lo = Point{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
        hi = Point{-lo[0], -lo[1]};
        for (const Point &p : m.vertices) {
            lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
            hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
        }
        const int nT = static_cast<int>(m.triangles.size());
        if (nT == 0 || !(hi[0] > lo[0]) || !(hi[1] > lo[1])) return;
        cell = std::sqrt((hi[0] - lo[0]) * (hi[1] - lo[1]) / nT);
        nx = std::max(1, static_cast<int>((hi[0] - lo[0]) / cell) + 1);
        ny = std::max(1, static_cast<int>((hi[1] - lo[1]) / cell) + 1);
        cells.assign(static_cast<size_t>(nx) * ny, {});
        for (int t = 0; t < nT; ++t) {
            const Triangle &tri = m.triangles[t];
            double x0 = hi[0], y0 = hi[1], x1 = lo[0], y1 = lo[1];
            for (int k = 0; k < 3; ++k) {
                const Point &p = m.vertices[tri[k]];
                x0 = std::min(x0, p[0]); y0 = std::min(y0, p[1]);
                x1 = std::max(x1, p[0]); y1 = std::max(y1, p[1]);
            }
            for (int j = clampY(y0); j <= clampY(y1); ++j)
                for (int i = clampX(x0); i <= clampX(x1); ++i)
                    cells[static_cast<size_t>(j) * nx + i].push_back(t);
        }
    }

    // A triangle holding q, closed (a point on an edge is in both), or -1.
    int locate(const Point &q) const {
        if (cells.empty()) return -1;
        const int i = clampX(q[0]), j = clampY(q[1]);
        for (const int t : cells[static_cast<size_t>(j) * nx + i]) {
            const Triangle &tri = mesh->triangles[t];
            const Point &a = mesh->vertices[tri[0]];
            const Point &b = mesh->vertices[tri[1]];
            const Point &c = mesh->vertices[tri[2]];
            const double s0 = cross2(b - a, q - a), s1 = cross2(c - b, q - b), s2 = cross2(a - c, q - c);
            const double eps = -1e-12 * cell * cell;
            if ((s0 >= eps && s1 >= eps && s2 >= eps) || (s0 <= -eps && s1 <= -eps && s2 <= -eps))
                return t;
        }
        return -1;
    }

private:
    int clampX(double x) const {
        return std::min(nx - 1, std::max(0, static_cast<int>((x - lo[0]) / cell)));
    }
    int clampY(double y) const {
        return std::min(ny - 1, std::max(0, static_cast<int>((y - lo[1]) / cell)));
    }

    const Mesh *mesh = nullptr;
    Point lo{0.0, 0.0}, hi{0.0, 0.0};
    double cell = 1.0;
    int nx = 0, ny = 0;
    std::vector<std::vector<int>> cells;
};

#endif // __TRIANGLE_LOCATOR_HXX__
