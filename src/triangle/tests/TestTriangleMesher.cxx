/**
 * TestTriangleMesher.cxx
 * ----------------------
 * CTest-based tests for TriangleMesher2D wrapper.
 *
 * Tests:
 *   1. Rectangle [0,2] x [0,1]
 *   2. Circle radius 0.5 centered at (0.5, 0.5)
 *   3. Circle radius 0.5 with an elliptical hole (rx=0.2, ry=0.1)
 *   4. Circle radius 0.5 with an inner elliptical region (two materials)
 *   5. Circle radius 1.0 with two elliptical holes
 *   6. Concentric shells (circles with radii 1.0, 0.8, 0.6, 0.4, 0.2)
 */

#include "../TriangleMesher.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace triangle_wrapper;

// Helper: generate a circle as a list of vertices and segments
void generateCircle(
    double cx, double cy, double radius, int npts,
    std::vector<std::array<double, 2>>& verts,
    std::vector<std::array<int, 2>>& segments,
    int base_index = 0)
{
  const double dtheta = 2.0 * M_PI / static_cast<double>(npts);
  for (int i = 0; i < npts; ++i) {
    double theta = i * dtheta;
    double x = cx + radius * std::cos(theta);
    double y = cy + radius * std::sin(theta);
    verts.push_back({x, y});
  }

  for (int i = 0; i < npts; ++i) {
    int a = base_index + i;
    int b = base_index + (i + 1) % npts;
    segments.push_back({a, b});
  }
}

// Helper: generate an ellipse as a list of vertices and segments
void generateEllipse(
    double cx, double cy, double rx, double ry, int npts,
    std::vector<std::array<double, 2>>& verts,
    std::vector<std::array<int, 2>>& segments,
    int base_index = 0)
{
  const double dtheta = 2.0 * M_PI / static_cast<double>(npts);
  for (int i = 0; i < npts; ++i) {
    double theta = i * dtheta;
    double x = cx + rx * std::cos(theta);
    double y = cy + ry * std::sin(theta);
    verts.push_back({x, y});
  }

  for (int i = 0; i < npts; ++i) {
    int a = base_index + i;
    int b = base_index + (i + 1) % npts;
    segments.push_back({a, b});
  }
}

// Helper: check basic mesh validity and write OBJ file
bool validateMesh(const TriangleMesher2D::MeshOutput& out, const std::string& name, const std::string& obj_filename) {
  if (out.verts.empty()) {
    std::cerr << name << ": FAIL - no vertices produced\n";
    return false;
  }
  if (out.triangles.empty()) {
    std::cerr << name << ": FAIL - no triangles produced\n";
    return false;
  }

  // Check all triangle indices are valid
  for (size_t i = 0; i < out.triangles.size(); ++i) {
    const auto& tri = out.triangles[i];
    for (int j = 0; j < 3; ++j) {
      if (tri[j] < 0 || static_cast<size_t>(tri[j]) >= out.verts.size()) {
        std::cerr << name << ": FAIL - invalid vertex index in triangle " << i << "\n";
        return false;
      }
    }
  }

  // Write OBJ file
  try {
    out.writeOBJ(obj_filename);
    std::cout << name << ": PASS (" << out.verts.size() << " vertices, "
              << out.triangles.size() << " triangles) -> " << obj_filename << "\n";
  } catch (const std::exception& e) {
    std::cerr << name << ": FAIL - could not write OBJ: " << e.what() << "\n";
    return false;
  }

  return true;
}

// Test 1: Rectangle [0,2] x [0,1]
bool testRectangle() {
  TriangleMesher2D::MeshInput input;

  // Four corners
  input.vertlist = {
    {0.0, 0.0},
    {2.0, 0.0},
    {2.0, 1.0},
    {0.0, 1.0}
  };

  // One loop with 4 segments forming the rectangle boundary
  std::vector<std::array<int, 2>> rect_segments = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0}
  };
  input.segment_loops.push_back(rect_segments);
  input.type.push_back(0); // exterior

  input.h = 0.05;

  TriangleMesher2D mesher;
  try {
    auto out = mesher.triangulate(input);
    return validateMesh(out, "Test 1 (Rectangle)", "test1.obj");
  } catch (const std::exception& e) {
    std::cerr << "Test 1 (Rectangle): FAIL - exception: " << e.what() << "\n";
    return false;
  }
}

// Test 2: Circle radius 0.5 centered at (0.5, 0.5)
bool testCircle() {
  TriangleMesher2D::MeshInput input;

  const int npts = 64;
  std::vector<std::array<int, 2>> circle_segments;

  generateCircle(0.5, 0.5, 0.5, npts, input.vertlist, circle_segments, 0);

  input.segment_loops.push_back(circle_segments);
  input.type.push_back(0); // exterior

  input.h = 0.02;

  TriangleMesher2D mesher;
  try {
    auto out = mesher.triangulate(input);
    return validateMesh(out, "Test 2 (Circle)", "test2.obj");
  } catch (const std::exception& e) {
    std::cerr << "Test 2 (Circle): FAIL - exception: " << e.what() << "\n";
    return false;
  }
}

// Test 3: Circle with an elliptical hole
// Outer circle: radius 0.5 centered at (0.5, 0.5)
// Hole: ellipse with rx=0.2, ry=0.1 at the same center
bool testCircleWithHole() {
  TriangleMesher2D::MeshInput input;

  const int npts_outer = 64;
  const int npts_inner = 32;

  // Outer circle
  std::vector<std::array<int, 2>> outer_segments;
  generateCircle(0.5, 0.5, 0.5, npts_outer, input.vertlist, outer_segments, 0);

  // Inner ellipse (hole)
  std::vector<std::array<int, 2>> inner_segments;
  generateEllipse(0.5, 0.5, 0.2, 0.1, npts_inner, input.vertlist, inner_segments, npts_outer);

  input.segment_loops.push_back(outer_segments);
  input.type.push_back(0); // exterior

  input.segment_loops.push_back(inner_segments);
  input.type.push_back(1); // hole

  input.h = 0.02;

  TriangleMesher2D mesher;
  try {
    auto out = mesher.triangulate(input);

    // Additional check: verify the hole by checking no triangles have centroid inside the ellipse
    bool hole_ok = true;
    for (const auto& tri : out.triangles) {
      const auto& v0 = out.verts[tri[0]];
      const auto& v1 = out.verts[tri[1]];
      const auto& v2 = out.verts[tri[2]];
      double cx = (v0[0] + v1[0] + v2[0]) / 3.0;
      double cy = (v0[1] + v1[1] + v2[1]) / 3.0;
      // Check if centroid is inside ellipse: (x-cx)^2/rx^2 + (y-cy)^2/ry^2 < 1
      double dx = (cx - 0.5) / 0.2;
      double dy = (cy - 0.5) / 0.1;
      if (dx * dx + dy * dy < 0.9) { // with tolerance
        hole_ok = false;
        break;
      }
    }

    if (!hole_ok) {
      std::cerr << "Test 3 (Circle with Hole): FAIL - triangles found inside hole\n";
      return false;
    }

    return validateMesh(out, "Test 3 (Circle with Ellipse Hole)", "test3.obj");
  } catch (const std::exception& e) {
    std::cerr << "Test 3 (Circle with Hole): FAIL - exception: " << e.what() << "\n";
    return false;
  }
}

// Test 4: Circle with an inner elliptical region (two materials)
// Outer circle: radius 0.5 centered at (0.5, 0.5)
// Inner region: ellipse with rx=0.2, ry=0.1 at the same center
bool testCircleWithRegion() {
  TriangleMesher2D::MeshInput input;

  const int npts_outer = 64;
  const int npts_inner = 32;

  // Outer circle
  std::vector<std::array<int, 2>> outer_segments;
  generateCircle(0.5, 0.5, 0.5, npts_outer, input.vertlist, outer_segments, 0);

  // Inner ellipse (region boundary)
  std::vector<std::array<int, 2>> inner_segments;
  generateEllipse(0.5, 0.5, 0.2, 0.1, npts_inner, input.vertlist, inner_segments, npts_outer);

  input.segment_loops.push_back(outer_segments);
  input.type.push_back(2); // exterior region (region 1)

  input.segment_loops.push_back(inner_segments);
  input.type.push_back(2); // inner region (region 2)

  // Assign region ids
  input.region_id.resize(2);
  input.region_id[0] = 1; // outer region
  input.region_id[1] = 2; // inner region

  // Provide explicit seed points for both regions
  // Outer region seed: somewhere between inner ellipse and outer circle
  input.loop_seed.push_back({0.5, 0.5 + 0.35}); // outside the ellipse but inside the circle
  // Inner region seed: at the center
  input.loop_seed.push_back({0.5, 0.5});

  input.h = 0.02;

  TriangleMesher2D mesher;
  try {
    auto out = mesher.triangulate(input);

    // Verify that we have triangles with different region ids
    bool has_region_1 = false;
    bool has_region_2 = false;
    for (int r : out.tri_regions) {
      if (r == 1) has_region_1 = true;
      if (r == 2) has_region_2 = true;
    }

    if (!has_region_1 || !has_region_2) {
      std::cerr << "Test 4 (Circle with Region): FAIL - expected triangles in both regions\n";
      std::cerr << "  Has region 1: " << (has_region_1 ? "yes" : "no") << "\n";
      std::cerr << "  Has region 2: " << (has_region_2 ? "yes" : "no") << "\n";
      return false;
    }

    return validateMesh(out, "Test 4 (Circle with Ellipse Region)", "test4.obj");
  } catch (const std::exception& e) {
    std::cerr << "Test 4 (Circle with Region): FAIL - exception: " << e.what() << "\n";
    return false;
  }
}

// Test 5: Circle with two elliptical holes
// Outer circle: radius 1.0 centered at (0.5, 0.5)
// Hole 1: ellipse at (0.5, 0.7) with rx=0.2, ry=0.1
// Hole 2: ellipse at (0.5, 0.2) with rx=0.15, ry=0.075
bool testCircleWithTwoHoles() {
  TriangleMesher2D::MeshInput input;

  const int npts_outer = 128;
  const int npts_hole = 32;

  // Outer circle
  std::vector<std::array<int, 2>> outer_segments;
  generateCircle(0.5, 0.5, 1.0, npts_outer, input.vertlist, outer_segments, 0);

  // First ellipse hole at (0.5, 0.7)
  std::vector<std::array<int, 2>> hole1_segments;
  int base1 = static_cast<int>(input.vertlist.size());
  generateEllipse(0.5, 0.7, 0.2, 0.1, npts_hole, input.vertlist, hole1_segments, base1);

  // Second ellipse hole at (0.5, 0.2)
  std::vector<std::array<int, 2>> hole2_segments;
  int base2 = static_cast<int>(input.vertlist.size());
  generateEllipse(0.5, 0.2, 0.15, 0.075, npts_hole, input.vertlist, hole2_segments, base2);

  input.segment_loops.push_back(outer_segments);
  input.type.push_back(0); // exterior

  input.segment_loops.push_back(hole1_segments);
  input.type.push_back(1); // hole

  input.segment_loops.push_back(hole2_segments);
  input.type.push_back(1); // hole

  input.h = 0.03;

  TriangleMesher2D mesher;
  try {
    auto out = mesher.triangulate(input);

    // Verify holes: no triangle centroids inside either ellipse
    bool holes_ok = true;
    for (const auto& tri : out.triangles) {
      const auto& v0 = out.verts[tri[0]];
      const auto& v1 = out.verts[tri[1]];
      const auto& v2 = out.verts[tri[2]];
      double cx = (v0[0] + v1[0] + v2[0]) / 3.0;
      double cy = (v0[1] + v1[1] + v2[1]) / 3.0;

      // Check hole 1: center (0.5, 0.7), rx=0.2, ry=0.1
      double dx1 = (cx - 0.5) / 0.2;
      double dy1 = (cy - 0.7) / 0.1;
      if (dx1 * dx1 + dy1 * dy1 < 0.9) {
        holes_ok = false;
        break;
      }

      // Check hole 2: center (0.5, 0.2), rx=0.15, ry=0.075
      double dx2 = (cx - 0.5) / 0.15;
      double dy2 = (cy - 0.2) / 0.075;
      if (dx2 * dx2 + dy2 * dy2 < 0.9) {
        holes_ok = false;
        break;
      }
    }

    if (!holes_ok) {
      std::cerr << "Test 5 (Circle with Two Holes): FAIL - triangles found inside a hole\n";
      return false;
    }

    return validateMesh(out, "Test 5 (Circle with Two Ellipse Holes)", "test5.obj");
  } catch (const std::exception& e) {
    std::cerr << "Test 5 (Circle with Two Holes): FAIL - exception: " << e.what() << "\n";
    return false;
  }
}

// Test 6: Concentric shells (multiple regions)
// Circles with radii 1.0, 0.8, 0.6, 0.4, 0.2 centered at (0.5, 0.5)
// Each shell is a separate region (5 regions total)
bool testConcentricShells() {
  TriangleMesher2D::MeshInput input;

  const double radii[] = {1.0, 0.8, 0.6, 0.4, 0.2};
  const int num_shells = 5;
  const int npts = 64;

  // Generate all circles from outermost to innermost
  for (int i = 0; i < num_shells; ++i) {
    std::vector<std::array<int, 2>> segments;
    int base = static_cast<int>(input.vertlist.size());
    generateCircle(0.5, 0.5, radii[i], npts, input.vertlist, segments, base);
    input.segment_loops.push_back(segments);
    input.type.push_back(2); // region
  }

  // Assign region ids (1 = outermost shell, 5 = innermost disk)
  input.region_id.resize(num_shells);
  for (int i = 0; i < num_shells; ++i) {
    input.region_id[i] = i + 1;
  }

  // Provide seed points for each region
  // Each seed should be inside the current circle but outside the next smaller one
  // Region 1: between r=1.0 and r=0.8 -> seed at r=0.9
  // Region 2: between r=0.8 and r=0.6 -> seed at r=0.7
  // Region 3: between r=0.6 and r=0.4 -> seed at r=0.5
  // Region 4: between r=0.4 and r=0.2 -> seed at r=0.3
  // Region 5: inside r=0.2 -> seed at center
  input.loop_seed.push_back({0.5, 0.5 + 0.9});   // region 1
  input.loop_seed.push_back({0.5, 0.5 + 0.7});   // region 2
  input.loop_seed.push_back({0.5, 0.5 + 0.5});   // region 3
  input.loop_seed.push_back({0.5, 0.5 + 0.3});   // region 4
  input.loop_seed.push_back({0.5, 0.5});         // region 5 (center)

  input.h = 0.03;

  TriangleMesher2D mesher;
  try {
    auto out = mesher.triangulate(input);

    // Verify that we have triangles in all 5 regions
    std::vector<bool> has_region(num_shells + 1, false);
    for (int r : out.tri_regions) {
      if (r >= 1 && r <= num_shells) {
        has_region[r] = true;
      }
    }

    bool all_regions_present = true;
    for (int i = 1; i <= num_shells; ++i) {
      if (!has_region[i]) {
        std::cerr << "Test 6 (Concentric Shells): FAIL - missing region " << i << "\n";
        all_regions_present = false;
      }
    }

    if (!all_regions_present) {
      return false;
    }

    return validateMesh(out, "Test 6 (Concentric Shells)", "test6.obj");
  } catch (const std::exception& e) {
    std::cerr << "Test 6 (Concentric Shells): FAIL - exception: " << e.what() << "\n";
    return false;
  }
}

int main(int argc, char* argv[]) {
  // If a test number is provided, run only that test
  int test_num = 0;
  if (argc > 1) {
    test_num = std::atoi(argv[1]);
  }

  bool all_passed = true;

  if (test_num == 0 || test_num == 1) {
    if (!testRectangle()) all_passed = false;
  }
  if (test_num == 0 || test_num == 2) {
    if (!testCircle()) all_passed = false;
  }
  if (test_num == 0 || test_num == 3) {
    if (!testCircleWithHole()) all_passed = false;
  }
  if (test_num == 0 || test_num == 4) {
    if (!testCircleWithRegion()) all_passed = false;
  }
  if (test_num == 0 || test_num == 5) {
    if (!testCircleWithTwoHoles()) all_passed = false;
  }
  if (test_num == 0 || test_num == 6) {
    if (!testConcentricShells()) all_passed = false;
  }

  return all_passed ? 0 : 1;
}
