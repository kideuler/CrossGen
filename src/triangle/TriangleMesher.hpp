#pragma once
/*
  TriangleMesher.hpp
  ------------------
  A small C++ wrapper around Jonathan Richard Shewchuk's "Triangle" library
  (triangle.h / triangle.cpp).

  This wrapper is designed for the following input style:

    std::vector<std::array<double,2>> vertlist;                  // global vertex list
    std::vector<std::vector<std::array<int,2>>> segment_loops;   // each loop is a set of connected segments
    std::vector<int> type;                                       // 0 = exterior, 1 = hole, 2 = region
    double h;                                                    // target average edge length

  and produces:

    std::vector<std::array<int,3>> triangles;       // triangle vertex indices (0-based)
    std::vector<std::array<double,2>> verts;        // output vertex list (includes Steiner points)
    std::vector<int> tri_regions;                   // per-triangle region/material id (Triangle "attribute")

  Notes
  -----
  - You must compile triangle.cpp with TRILIBRARY defined (e.g., -DTRILIBRARY).
  - This wrapper uses Triangle's zero-based indexing by default (switch 'z'),
    so all segment endpoints are expected to be 0-based indices into vertlist.
  - For holes and regions, Triangle needs a "seed point" strictly inside the
    loop. By default, this wrapper computes one automatically; for highly
    concave polygons, automatic seeding can fail. You can provide explicit seed
    points per loop to make it robust.

  Typical switches used by this wrapper:
    p  : PSLG input (segments)
    z  : zero-based indexing
    A  : assign region attributes to triangles
    q  : quality mesh (minimum angle)
    a  : maximum triangle area (derived from h unless you override)
    Q  : quiet
    P  : suppress segment output (smaller out-struct)

  This file is header-only. It depends on triangle.h being available and
  triangle.cpp being linked into your program.
*/

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

// Include Shewchuk Triangle header.
//
// IMPORTANT:
// - The provided triangle.cpp in this repository is a C++-compiled translation unit.
//   Make sure you compile triangle.cpp with -DTRILIBRARY and -DANSI_DECLARATORS.
// - If you instead compile Triangle as C (triangle.c), define TRIANGLE_WRAPPER_C_LINKAGE
//   before including this header.
#if defined(TRIANGLE_WRAPPER_C_LINKAGE) && defined(__cplusplus)
extern "C" {
#endif
#include "triangle.h"
#if defined(TRIANGLE_WRAPPER_C_LINKAGE) && defined(__cplusplus)
} // extern \"C\"
#endif

namespace triangle_wrapper {

// ---------------------------
// Helper: robust-ish geometry
// ---------------------------

namespace detail {

inline bool isFinite(double x) {
  return std::isfinite(x);
}

inline double sqr(double x) { return x * x; }

inline double hypot2(double x, double y) {
  return std::sqrt(x * x + y * y);
}

// Even-odd ray casting test. Returns true for points strictly inside.
// Points on the boundary are treated as inside (practically), by using a tiny epsilon.
inline bool pointInPolygon(
    const std::array<double, 2>& p,
    const std::vector<int>& poly,
    const std::vector<std::array<double, 2>>& verts)
{
  const double px = p[0];
  const double py = p[1];

  bool inside = false;
  const std::size_t n = poly.size();
  if (n < 3) return false;

  for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
    const auto& vi = verts[static_cast<std::size_t>(poly[i])];
    const auto& vj = verts[static_cast<std::size_t>(poly[j])];

    const double xi = vi[0], yi = vi[1];
    const double xj = vj[0], yj = vj[1];

    // Check if edge (vj->vi) crosses horizontal ray to the right of p.
    const bool intersect = ((yi > py) != (yj > py)) &&
                           (px < (xj - xi) * (py - yi) / (yj - yi + 1e-300) + xi);
    if (intersect) inside = !inside;
  }
  return inside;
}

// Compute signed area using shoelace formula (positive if CCW).
inline double signedArea(
    const std::vector<int>& poly,
    const std::vector<std::array<double, 2>>& verts)
{
  const std::size_t n = poly.size();
  if (n < 3) return 0.0;

  long double a = 0.0L;
  for (std::size_t i = 0; i < n; ++i) {
    const auto& p0 = verts[static_cast<std::size_t>(poly[i])];
    const auto& p1 = verts[static_cast<std::size_t>(poly[(i + 1) % n])];
    a += static_cast<long double>(p0[0]) * static_cast<long double>(p1[1]) -
         static_cast<long double>(p1[0]) * static_cast<long double>(p0[1]);
  }
  return static_cast<double>(0.5L * a);
}

// Polygon centroid (area-weighted). If degenerate, returns average of vertices.
inline std::array<double, 2> polygonCentroid(
    const std::vector<int>& poly,
    const std::vector<std::array<double, 2>>& verts)
{
  const std::size_t n = poly.size();
  if (n < 3) {
    // Average
    std::array<double, 2> c{0.0, 0.0};
    for (int idx : poly) {
      c[0] += verts[static_cast<std::size_t>(idx)][0];
      c[1] += verts[static_cast<std::size_t>(idx)][1];
    }
    if (!poly.empty()) {
      c[0] /= static_cast<double>(poly.size());
      c[1] /= static_cast<double>(poly.size());
    }
    return c;
  }

  long double A2 = 0.0L; // twice area
  long double Cx = 0.0L;
  long double Cy = 0.0L;

  for (std::size_t i = 0; i < n; ++i) {
    const auto& p0 = verts[static_cast<std::size_t>(poly[i])];
    const auto& p1 = verts[static_cast<std::size_t>(poly[(i + 1) % n])];

    const long double x0 = p0[0], y0 = p0[1];
    const long double x1 = p1[0], y1 = p1[1];

    const long double cross = x0 * y1 - x1 * y0;
    A2 += cross;
    Cx += (x0 + x1) * cross;
    Cy += (y0 + y1) * cross;
  }

  if (std::fabs(A2) < 1e-300L) {
    // Degenerate: average
    std::array<double, 2> c{0.0, 0.0};
    for (int idx : poly) {
      c[0] += verts[static_cast<std::size_t>(idx)][0];
      c[1] += verts[static_cast<std::size_t>(idx)][1];
    }
    c[0] /= static_cast<double>(poly.size());
    c[1] /= static_cast<double>(poly.size());
    return c;
  }

  const long double inv = 1.0L / (3.0L * A2); // since centroid divides by 6A, and A2=2A
  return { static_cast<double>(Cx * inv), static_cast<double>(Cy * inv) };
}

// Reconstruct a simple loop vertex order from an unordered list of segments.
// Assumes each vertex in the loop has degree 2 (simple cycle).
inline std::vector<int> orderLoopVerticesFromSegments(
    const std::vector<std::array<int, 2>>& segments,
    std::size_t num_global_vertices)
{
  if (segments.size() < 3) {
    throw std::runtime_error("Segment loop must have at least 3 segments.");
  }

  std::vector<std::vector<int>> adj(num_global_vertices);
  std::unordered_set<int> used;
  used.reserve(segments.size() * 2);

  for (const auto& s : segments) {
    const int a = s[0];
    const int b = s[1];
    if (a < 0 || b < 0 ||
        static_cast<std::size_t>(a) >= num_global_vertices ||
        static_cast<std::size_t>(b) >= num_global_vertices)
    {
      throw std::runtime_error("A segment endpoint index is out of range of vertlist.");
    }
    if (a == b) {
      throw std::runtime_error("A segment has identical endpoints (degenerate).");
    }
    adj[static_cast<std::size_t>(a)].push_back(b);
    adj[static_cast<std::size_t>(b)].push_back(a);
    used.insert(a);
    used.insert(b);
  }

  for (int v : used) {
    const auto& nbrs = adj[static_cast<std::size_t>(v)];
    if (nbrs.size() != 2) {
      throw std::runtime_error("Segment loop is not a simple cycle (a vertex does not have degree 2).");
    }
  }

  const int start = segments[0][0];
  int prev = start;
  int curr = segments[0][1];

  std::vector<int> poly;
  poly.reserve(segments.size());
  poly.push_back(start);

  // Walk the cycle.
  for (std::size_t iter = 0; iter < segments.size(); ++iter) {
    poly.push_back(curr);

    const auto& nbrs = adj[static_cast<std::size_t>(curr)];
    const int next = (nbrs[0] == prev) ? nbrs[1] : nbrs[0];

    prev = curr;
    curr = next;

    if (curr == start) {
      break;
    }

    // Safety: prevent infinite loops if the input isn't actually a cycle.
    if (poly.size() > segments.size() + 1) {
      throw std::runtime_error("Failed to reconstruct a loop: detected non-cyclic walk.");
    }
  }

  // poly currently includes start twice? (it shouldn't) - remove if needed.
  if (!poly.empty() && poly.back() == start) {
    poly.pop_back();
  }

  if (poly.size() != segments.size()) {
    // In a simple cycle, #segments == #vertices.
    throw std::runtime_error("Segment loop reconstruction failed: segments do not form a single simple loop.");
  }
  return poly;
}

// Compute a seed point that (attempts to) lie inside the loop.
// Strategy:
//  1) polygon centroid (area-weighted)
//  2) if centroid not inside, try grid sampling in bbox
//  3) if still not inside, try midpoint+epsilon offset along inward normal on edges
inline std::array<double, 2> computeInteriorSeedPoint(
    const std::vector<int>& poly,
    const std::vector<std::array<double, 2>>& verts)
{
  // Bounding box.
  double minx =  std::numeric_limits<double>::infinity();
  double miny =  std::numeric_limits<double>::infinity();
  double maxx = -std::numeric_limits<double>::infinity();
  double maxy = -std::numeric_limits<double>::infinity();

  for (int idx : poly) {
    const auto& v = verts[static_cast<std::size_t>(idx)];
    minx = std::min(minx, v[0]);
    miny = std::min(miny, v[1]);
    maxx = std::max(maxx, v[0]);
    maxy = std::max(maxy, v[1]);
  }

  if (!isFinite(minx) || !isFinite(miny) || !isFinite(maxx) || !isFinite(maxy)) {
    throw std::runtime_error("Non-finite coordinates in loop.");
  }

  // 1) Centroid
  std::array<double, 2> seed = polygonCentroid(poly, verts);
  if (pointInPolygon(seed, poly, verts)) {
    return seed;
  }

  // 2) Grid sampling inside bbox
  const int N = 25;
  const double dx = (maxx - minx);
  const double dy = (maxy - miny);
  if (dx > 0.0 && dy > 0.0) {
    for (int iy = 0; iy < N; ++iy) {
      for (int ix = 0; ix < N; ++ix) {
        const double x = minx + (static_cast<double>(ix) + 0.5) / static_cast<double>(N) * dx;
        const double y = miny + (static_cast<double>(iy) + 0.5) / static_cast<double>(N) * dy;
        const std::array<double, 2> p{x, y};
        if (pointInPolygon(p, poly, verts)) {
          return p;
        }
      }
    }
  }

  // 3) Edge midpoint + inward epsilon
  const double area = signedArea(poly, verts);
  const bool ccw = area > 0.0;
  const double diag = hypot2(dx, dy);
  const double eps = (diag > 0.0) ? (1e-6 * diag) : 1e-6;

  const std::size_t n = poly.size();
  for (std::size_t i = 0; i < n; ++i) {
    const auto& p0 = verts[static_cast<std::size_t>(poly[i])];
    const auto& p1 = verts[static_cast<std::size_t>(poly[(i + 1) % n])];

    const double ex = p1[0] - p0[0];
    const double ey = p1[1] - p0[1];
    const double el = hypot2(ex, ey);
    if (el < 1e-300) continue;

    // Left normal is (-ey, ex). Right normal is (ey, -ex).
    double nx = ccw ? (-ey / el) : (ey / el);
    double ny = ccw ? ( ex / el) : (-ex / el);

    const std::array<double, 2> mid{0.5 * (p0[0] + p1[0]), 0.5 * (p0[1] + p1[1])};
    const std::array<double, 2> p{mid[0] + eps * nx, mid[1] + eps * ny};

    if (pointInPolygon(p, poly, verts)) {
      return p;
    }
  }

  throw std::runtime_error(
      "Failed to automatically compute an interior seed point. "
      "Provide an explicit seed point for this loop.");
}

inline void trifree_if_not_null(void* p) {
  if (!p) return;
  // triangle.h declares trifree(int*). It is safe to cast.
  trifree(reinterpret_cast<int*>(p));
}

inline void free_triangle_output(triangulateio& out) {
  // These are allocated by Triangle (if requested). Free all except holelist/regionlist:
  // holelist & regionlist are only copied from 'in' to 'out' (not allocated).
  trifree_if_not_null(out.pointlist);
  trifree_if_not_null(out.pointattributelist);
  trifree_if_not_null(out.pointmarkerlist);

  trifree_if_not_null(out.trianglelist);
  trifree_if_not_null(out.triangleattributelist);
  trifree_if_not_null(out.trianglearealist);
  trifree_if_not_null(out.neighborlist);

  trifree_if_not_null(out.segmentlist);
  trifree_if_not_null(out.segmentmarkerlist);

  trifree_if_not_null(out.edgelist);
  trifree_if_not_null(out.edgemarkerlist);
  trifree_if_not_null(out.normlist);

  // Leave out.holelist and out.regionlist alone.
  out = triangulateio{}; // reset to nulls/zeros to avoid accidental double-free
}

} // namespace detail

// -------------------------
// Public wrapper structures
// -------------------------

class TriangleMesher2D {
public:
  // Loop types (matches your description).
  enum class LoopType : int {
    Exterior = 0,
    Hole     = 1,
    Region   = 2
  };

  struct Options {
    // Mesh sizing:
    //
    // Triangle controls sizing mainly via area constraints (-a).
    // We map your target edge length h to a max triangle area:
    //   A_max = area_from_h_factor * h^2
    //
    // For an equilateral triangle of side length h, area = sqrt(3)/4 * h^2.
    double area_from_h_factor = 0.43301270189221932338; // sqrt(3)/4

    // If > 0, overrides the derived max area and uses this value directly.
    // If <= 0, max area is derived from h (if h > 0). If both are <= 0, no area constraint is used.
    double max_area_override = 0.0;

    // Quality meshing (minimum angle). If <= 0, quality switch is omitted.
    // Typical values: 20..35. Too high can make meshing fail.
    double min_angle_degrees = 20.0;

    // Triangle verbosity: if false, adds 'Q' (quiet).
    bool verbose = false;

    // If true, adds 'e' and returns edges & edge markers.
    bool output_edges = false;

    // If true, adds 'n' and returns triangle neighbors.
    bool output_neighbors = false;

    // Steiner / segment splitting controls:
    // -Y  : suppresses boundary segment splitting
    // -YY : suppresses all segment splitting (including internal segments)
    bool suppress_boundary_splitting = false;
    bool suppress_all_splitting = false;

    // If >= 0, adds 'S<maxSteiner>' (maximum Steiner points).
    // Set to 0 to forbid Steiner points (may cause failure if constraints can't be met).
    int max_steiner_points = -1;

    // If true, compute only the constrained Delaunay triangulation of the
    // input vertices and segments — no quality refinement, no area
    // constraints, no Steiner points, and no segment splitting.
    // Overrides min_angle_degrees, max_area_override, max_steiner_points,
    // suppress_boundary_splitting, and suppress_all_splitting.
    bool just_delaunay = false;

    // Extra raw switches appended verbatim (advanced use).
    // Example: "D" for conforming Delaunay, "O2" for second-order elements, etc.
    std::string extra_switches;

    // Segment output is usually not needed; add 'P' by default to suppress it.
    bool suppress_segment_output = true;

    // If true, we will require out.numberofcorners == 3.
    // If false, we will take the first 3 corner indices and ignore the rest.
    bool require_linear_triangles = true;
  };

  struct MeshInput {
    std::vector<std::array<double, 2>> vertlist;
    std::vector<std::vector<std::array<int, 2>>> segment_loops;
    std::vector<int> type; // 0 exterior, 1 hole, 2 region
    double h = 0.0;

    // Recommended additions:

    // Optional per-loop region id (only used if type[i] == 2). If empty, ids are assigned 1..K in order.
    std::vector<int> region_id;

    // Optional per-loop interior seed point. If provided, size must equal segment_loops.size().
    // Used when type[i] is Hole or Region. Ignored for Exterior loops.
    std::vector<std::array<double, 2>> loop_seed;

    // Optional per-loop segment marker. If empty, marker = i+1 for all loops.
    // Markers can be used to tag boundary edges if output_edges = true.
    std::vector<int> loop_marker;
  };

  struct MeshOutput {
    std::vector<std::array<double, 2>> verts;
    std::vector<std::array<int, 3>> triangles;
    std::vector<int> tri_regions; // same size as triangles

    // Optional outputs:
    std::vector<std::array<int, 2>> edges;  // if Options::output_edges
    std::vector<int> edge_markers;          // if Options::output_edges
    std::vector<std::array<int, 3>> neighbors; // if Options::output_neighbors (three neighbors per triangle, -1 for boundary)

    // Write mesh to OBJ file format.
    // OBJ uses 1-based indices and supports 2D meshes (z=0).
    void writeOBJ(const std::string& filename) const {
      std::ofstream ofs(filename);
      if (!ofs) {
        throw std::runtime_error("Failed to open file for writing: " + filename);
      }

      ofs << "# OBJ file generated by TriangleMesher2D\n";
      ofs << "# Vertices: " << verts.size() << "\n";
      ofs << "# Triangles: " << triangles.size() << "\n";
      ofs << std::setprecision(16);

      // Write vertices (2D -> 3D with z=0)
      for (const auto& v : verts) {
        ofs << "v " << v[0] << " " << v[1] << " 0\n";
      }

      // Write faces (1-based indexing)
      for (const auto& tri : triangles) {
        ofs << "f " << (tri[0] + 1) << " " << (tri[1] + 1) << " " << (tri[2] + 1) << "\n";
      }

      ofs.close();
      if (!ofs) {
        throw std::runtime_error("Error writing to file: " + filename);
      }
    }
  };

  TriangleMesher2D() = default;
  explicit TriangleMesher2D(Options opt) : opt_(std::move(opt)) {}

  const Options& options() const { return opt_; }
  void setOptions(const Options& opt) { opt_ = opt; }

  // Main entry point using a structured input.
  MeshOutput triangulate(const MeshInput& in) const {
    validateInput_(in);

    // Convert input vertices to TRI_REAL.
    std::vector<TRI_REAL> pointlist;
    pointlist.resize(in.vertlist.size() * 2);

    for (std::size_t i = 0; i < in.vertlist.size(); ++i) {
      const double x = in.vertlist[i][0];
      const double y = in.vertlist[i][1];
      if (!detail::isFinite(x) || !detail::isFinite(y)) {
        throw std::runtime_error("Non-finite coordinate in vertlist.");
      }
      pointlist[2 * i + 0] = static_cast<TRI_REAL>(x);
      pointlist[2 * i + 1] = static_cast<TRI_REAL>(y);
    }

    // Flatten segments and build marker list.
    std::size_t total_segments = 0;
    for (const auto& loop : in.segment_loops) total_segments += loop.size();

    std::vector<int> segmentlist;
    segmentlist.resize(total_segments * 2);

    std::vector<int> segmentmarkerlist;
    segmentmarkerlist.resize(total_segments);

    std::size_t segCursor = 0;
    for (std::size_t li = 0; li < in.segment_loops.size(); ++li) {
      const int marker = loopMarker_(in, li);
      const auto& loop = in.segment_loops[li];
      for (const auto& s : loop) {
        const int a = s[0];
        const int b = s[1];
        if (a < 0 || b < 0 ||
            static_cast<std::size_t>(a) >= in.vertlist.size() ||
            static_cast<std::size_t>(b) >= in.vertlist.size())
        {
          throw std::runtime_error("A segment endpoint index is out of range of vertlist.");
        }

        segmentlist[2 * segCursor + 0] = a;
        segmentlist[2 * segCursor + 1] = b;
        segmentmarkerlist[segCursor] = marker;
        ++segCursor;
      }
    }

    // Build holelist and regionlist.
    std::vector<TRI_REAL> holelist;   // 2 reals per hole
    std::vector<TRI_REAL> regionlist; // 4 reals per region

    // For region id auto assignment.
    int next_auto_region_id = 1;

    for (std::size_t li = 0; li < in.segment_loops.size(); ++li) {
      const int t = in.type[li];
      if (t != static_cast<int>(LoopType::Exterior) &&
          t != static_cast<int>(LoopType::Hole) &&
          t != static_cast<int>(LoopType::Region))
      {
        throw std::runtime_error("type[i] must be 0 (exterior), 1 (hole), or 2 (region).");
      }

      if (t == static_cast<int>(LoopType::Exterior)) continue;

      std::array<double, 2> seed{};
      const bool has_user_seed = (!in.loop_seed.empty() && in.loop_seed.size() == in.segment_loops.size());
      if (has_user_seed) {
        seed = in.loop_seed[li];
        if (!detail::isFinite(seed[0]) || !detail::isFinite(seed[1])) {
          throw std::runtime_error("Non-finite loop_seed coordinate.");
        }
      } else {
        // Reconstruct polygon order from segments and compute an interior point.
        const auto poly = detail::orderLoopVerticesFromSegments(in.segment_loops[li], in.vertlist.size());
        seed = detail::computeInteriorSeedPoint(poly, in.vertlist);
      }

      if (t == static_cast<int>(LoopType::Hole)) {
        holelist.push_back(static_cast<TRI_REAL>(seed[0]));
        holelist.push_back(static_cast<TRI_REAL>(seed[1]));
      } else if (t == static_cast<int>(LoopType::Region)) {
        int rid = regionId_(in, li, next_auto_region_id);
        // Region entry: x, y, attribute, max area
        regionlist.push_back(static_cast<TRI_REAL>(seed[0]));
        regionlist.push_back(static_cast<TRI_REAL>(seed[1]));
        regionlist.push_back(static_cast<TRI_REAL>(rid));
        regionlist.push_back(static_cast<TRI_REAL>(effectiveMaxArea_(in.h)));
      }
    }

    // Prepare Triangle input/output.
    triangulateio tri_in = triangulateio{};
    triangulateio tri_out = triangulateio{};

    tri_in.numberofpoints = static_cast<int>(in.vertlist.size());
    tri_in.numberofpointattributes = 0;
    tri_in.pointlist = pointlist.data();
    tri_in.pointattributelist = nullptr;
    tri_in.pointmarkerlist = nullptr;

    tri_in.numberofsegments = static_cast<int>(total_segments);
    tri_in.segmentlist = segmentlist.data();
    tri_in.segmentmarkerlist = segmentmarkerlist.data();

    tri_in.numberofholes = static_cast<int>(holelist.size() / 2);
    tri_in.holelist = holelist.empty() ? nullptr : holelist.data();

    tri_in.numberofregions = static_cast<int>(regionlist.size() / 4);
    tri_in.regionlist = regionlist.empty() ? nullptr : regionlist.data();

    // Ensure Triangle allocates output arrays we care about by setting to NULL.
    tri_out.pointlist = nullptr;
    tri_out.pointattributelist = nullptr;
    tri_out.pointmarkerlist = nullptr;

    tri_out.trianglelist = nullptr;
    tri_out.triangleattributelist = nullptr;
    tri_out.neighborlist = nullptr;

    tri_out.segmentlist = nullptr;
    tri_out.segmentmarkerlist = nullptr;

    tri_out.edgelist = nullptr;
    tri_out.edgemarkerlist = nullptr;
    tri_out.normlist = nullptr;

    // Build switches.
    const std::string switches = buildSwitches_(in.h, tri_in.numberofregions > 0);
    std::vector<char> swbuf(switches.begin(), switches.end());
    swbuf.push_back('\0');

    // Call Triangle.
    // NOTE: Triangle is not guaranteed to be thread-safe.
    ::triangulate(swbuf.data(), &tri_in, &tri_out, nullptr);

    // Extract output.
    MeshOutput out;
    out.verts.resize(static_cast<std::size_t>(tri_out.numberofpoints));
    for (int i = 0; i < tri_out.numberofpoints; ++i) {
      out.verts[static_cast<std::size_t>(i)][0] = static_cast<double>(tri_out.pointlist[2 * i + 0]);
      out.verts[static_cast<std::size_t>(i)][1] = static_cast<double>(tri_out.pointlist[2 * i + 1]);
    }

    if (tri_out.numberofcorners < 3) {
      detail::free_triangle_output(tri_out);
      throw std::runtime_error("Triangle returned numberofcorners < 3 (unexpected).");
    }
    if (opt_.require_linear_triangles && tri_out.numberofcorners != 3) {
      detail::free_triangle_output(tri_out);
      throw std::runtime_error(
          "Triangle returned higher-order elements (numberofcorners != 3). "
          "Disable Options::require_linear_triangles if you want to accept this and keep only the first 3 corners.");
    }

    const int corners = tri_out.numberofcorners;
    out.triangles.resize(static_cast<std::size_t>(tri_out.numberoftriangles));
    for (int t = 0; t < tri_out.numberoftriangles; ++t) {
      const int base = t * corners;
      out.triangles[static_cast<std::size_t>(t)] = {
        tri_out.trianglelist[base + 0],
        tri_out.trianglelist[base + 1],
        tri_out.trianglelist[base + 2]
      };
    }

    // Triangle attributes -> region ids.
    out.tri_regions.resize(static_cast<std::size_t>(tri_out.numberoftriangles), 0);
    if (tri_out.triangleattributelist != nullptr && tri_out.numberoftriangleattributes > 0) {
      const int aStride = tri_out.numberoftriangleattributes;
      for (int t = 0; t < tri_out.numberoftriangles; ++t) {
        const TRI_REAL a0 = tri_out.triangleattributelist[t * aStride];
        const long long rid_ll = llround(static_cast<double>(a0));
        out.tri_regions[static_cast<std::size_t>(t)] = static_cast<int>(rid_ll);
      }
    }

    // Optional edges.
    if (opt_.output_edges) {
      out.edges.resize(static_cast<std::size_t>(tri_out.numberofedges));
      out.edge_markers.resize(static_cast<std::size_t>(tri_out.numberofedges), 0);
      for (int e = 0; e < tri_out.numberofedges; ++e) {
        out.edges[static_cast<std::size_t>(e)] = {
          tri_out.edgelist[2 * e + 0],
          tri_out.edgelist[2 * e + 1]
        };
        if (tri_out.edgemarkerlist) {
          out.edge_markers[static_cast<std::size_t>(e)] = tri_out.edgemarkerlist[e];
        }
      }
    }

    // Optional neighbors.
    if (opt_.output_neighbors && tri_out.neighborlist) {
      out.neighbors.resize(static_cast<std::size_t>(tri_out.numberoftriangles));
      for (int t = 0; t < tri_out.numberoftriangles; ++t) {
        out.neighbors[static_cast<std::size_t>(t)] = {
          tri_out.neighborlist[3 * t + 0],
          tri_out.neighborlist[3 * t + 1],
          tri_out.neighborlist[3 * t + 2]
        };
      }
    }

    // Free Triangle-owned memory.
    detail::free_triangle_output(tri_out);
    return out;
  }

  // Convenience overload that matches your "legacy" signature:
  void triangulate(
      const std::vector<std::array<double, 2>>& vertlist,
      const std::vector<std::vector<std::array<int, 2>>>& segment_loops,
      const std::vector<int>& type,
      double h,
      std::vector<std::array<int, 3>>& triangles,
      std::vector<std::array<double, 2>>& verts,
      std::vector<int>& tri_regions) const
  {
    MeshInput in;
    in.vertlist = vertlist;
    in.segment_loops = segment_loops;
    in.type = type;
    in.h = h;

    MeshOutput out = triangulate(in);
    triangles = std::move(out.triangles);
    verts = std::move(out.verts);
    tri_regions = std::move(out.tri_regions);
  }

private:
  Options opt_;

  void validateInput_(const MeshInput& in) const {
    if (in.vertlist.empty()) {
      throw std::runtime_error("vertlist is empty.");
    }
    if (in.segment_loops.empty()) {
      throw std::runtime_error("segment_loops is empty (need at least an exterior boundary).");
    }
    if (in.type.size() != in.segment_loops.size()) {
      throw std::runtime_error("type must have the same size as segment_loops.");
    }
    if (!in.loop_seed.empty() && in.loop_seed.size() != in.segment_loops.size()) {
      throw std::runtime_error("If provided, loop_seed must have the same size as segment_loops.");
    }
    if (!in.region_id.empty() && in.region_id.size() != in.segment_loops.size()) {
      throw std::runtime_error("If provided, region_id must have the same size as segment_loops.");
    }
    if (!in.loop_marker.empty() && in.loop_marker.size() != in.segment_loops.size()) {
      throw std::runtime_error("If provided, loop_marker must have the same size as segment_loops.");
    }

    // Basic segment validation.
    for (std::size_t li = 0; li < in.segment_loops.size(); ++li) {
      const auto& loop = in.segment_loops[li];
      if (loop.size() < 3) {
        throw std::runtime_error("Each segment loop must have at least 3 segments.");
      }
      for (const auto& s : loop) {
        const int a = s[0];
        const int b = s[1];
        if (a < 0 || b < 0 ||
            static_cast<std::size_t>(a) >= in.vertlist.size() ||
            static_cast<std::size_t>(b) >= in.vertlist.size())
        {
          throw std::runtime_error("A segment endpoint index is out of range of vertlist.");
        }
      }
    }

    // If h is given, it must be positive.
    if (in.h < 0.0) {
      throw std::runtime_error("h must be >= 0.");
    }
  }

  // Determine max triangle area based on h and options.
  double effectiveMaxArea_(double h) const {
    if (opt_.max_area_override > 0.0) return opt_.max_area_override;
    if (h > 0.0) return opt_.area_from_h_factor * h * h;
    return 0.0;
  }

  // Determine loop marker.
  int loopMarker_(const MeshInput& in, std::size_t li) const {
    if (!in.loop_marker.empty()) return in.loop_marker[li];
    return static_cast<int>(li + 1); // default: 1..N
  }

  // Determine region id.
  int regionId_(const MeshInput& in, std::size_t li, int& next_auto_region_id) const {
    int rid = 0;
    if (!in.region_id.empty()) {
      rid = in.region_id[li];
    } else {
      rid = next_auto_region_id++;
    }
    return rid;
  }

  // Build Triangle switch string.
  std::string buildSwitches_(double h, bool has_regions) const {
    std::ostringstream ss;
    ss.imbue(std::locale::classic());

    // Base:
    ss << "pzA";
    if (!opt_.verbose) ss << "Q";
    if (opt_.suppress_segment_output) ss << "P";

    if (opt_.just_delaunay) {
      // Constrained Delaunay only: no quality, no area constraint,
      // no Steiner points, no segment splitting.
      ss << "YYS0";
    } else {
      // Quality:
      if (opt_.min_angle_degrees > 0.0) {
        ss << "q" << std::setprecision(16) << opt_.min_angle_degrees;
      }

      // Area constraint:
      // NOTE: Triangle's switch parser does not handle scientific notation
      // (e.g. 1.08e-05) because it interprets the 'e' as a switch character.
      // We must use std::fixed to emit a plain decimal representation.
      const double max_area = effectiveMaxArea_(h);
      if (max_area > 0.0) {
        ss << "a" << std::fixed << std::setprecision(20) << max_area;
      }

      // Splitting control:
      if (opt_.suppress_all_splitting) {
        ss << "YY";
      } else if (opt_.suppress_boundary_splitting) {
        ss << "Y";
      }

      // Steiner limit:
      if (opt_.max_steiner_points >= 0) {
        ss << "S" << opt_.max_steiner_points;
      }
    }

    // Optional outputs:
    if (opt_.output_edges) ss << "e";
    if (opt_.output_neighbors) ss << "n";

    // Extra switches:
    if (!opt_.extra_switches.empty()) {
      ss << opt_.extra_switches;
    }

    // If there are no regions, Triangle will still run fine. 'A' just produces 0 attributes.
    (void)has_regions;

    return ss.str();
  }
};

} // namespace triangle_wrapper
