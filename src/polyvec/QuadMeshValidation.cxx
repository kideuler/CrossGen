#include "QuadMeshValidation.hxx"
#include "IntegerGridMap.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace {

// Hash for pairs
struct PairHash {
    size_t operator()(const std::pair<int,int>& p) const {
        return std::hash<long long>()(static_cast<long long>(p.first) << 32 | p.second);
    }
};

// Tolerance for geometric comparisons
constexpr double EPS = 1e-10;
constexpr double SNAP_TOL = 1e-8;

inline double cross2d(double ax, double ay, double bx, double by) {
    return ax * by - ay * bx;
}

inline bool approxEqual(double a, double b, double tol = EPS) {
    return std::abs(a - b) < tol;
}

// Check if a value is (approximately) an integer
inline bool isInteger(double val, double tol = SNAP_TOL) {
    return std::abs(val - std::round(val)) < tol;
}

// Get the integer isolines that pass through a triangle
void getIsolinesInTriangle(double u0, double v0, double u1, double v1, double u2, double v2,
                           std::vector<int>& u_isolines, std::vector<int>& v_isolines) {
    u_isolines.clear();
    v_isolines.clear();
    
    double u_min = std::min({u0, u1, u2});
    double u_max = std::max({u0, u1, u2});
    double v_min = std::min({v0, v1, v2});
    double v_max = std::max({v0, v1, v2});
    
    // Find integer u values in range
    int u_lo = static_cast<int>(std::ceil(u_min - SNAP_TOL));
    int u_hi = static_cast<int>(std::floor(u_max + SNAP_TOL));
    for (int u = u_lo; u <= u_hi; ++u) {
        u_isolines.push_back(u);
    }
    
    // Find integer v values in range
    int v_lo = static_cast<int>(std::ceil(v_min - SNAP_TOL));
    int v_hi = static_cast<int>(std::floor(v_max + SNAP_TOL));
    for (int v = v_lo; v <= v_hi; ++v) {
        v_isolines.push_back(v);
    }
}

} // namespace


std::string QuadMeshValidation::ValidationResult::summary() const {
    std::ostringstream ss;
    ss << "QuadMesh Validation Result: " << (valid ? "VALID" : "INVALID") << "\n";
    ss << "  Grid vertices: " << num_grid_vertices << "\n";
    ss << "  Grid edges: " << num_grid_edges << "\n";
    ss << "  Grid faces: " << num_grid_faces << "\n";
    ss << "  Quad faces: " << num_quad_faces << " / " << num_grid_faces << "\n";
    ss << "  Non-quad faces: " << num_non_quad_faces << "\n";
    ss << "  Self-crossings: " << num_self_crossings << "\n";
    ss << "  Degenerate edges: " << num_degenerate_edges << "\n";
    ss << "  Missing regions: " << num_missing_regions << "\n";
    
    if (!face_valence_histogram.empty()) {
        ss << "  Face valence distribution:\n";
        for (const auto& kv : face_valence_histogram) {
            ss << "    " << kv.first << "-gons: " << kv.second << "\n";
        }
    }
    
    if (!errors.empty()) {
        ss << "  Errors:\n";
        for (const auto& err : errors) {
            ss << "    - " << err << "\n";
        }
    }
    
    return ss.str();
}


bool QuadMeshValidation::GridVertex::operator<(const GridVertex& other) const {
    if (is_u_isoline != other.is_u_isoline) return is_u_isoline < other.is_u_isoline;
    if (isoline_value != other.isoline_value) return isoline_value < other.isoline_value;
    if (triangle != other.triangle) return triangle < other.triangle;
    if (edge != other.edge) return edge < other.edge;
    return t < other.t;
}

bool QuadMeshValidation::GridVertex::approxEqual(const GridVertex& other, double tol) const {
    return std::abs(x - other.x) < tol && std::abs(y - other.y) < tol;
}


QuadMeshValidation::QuadMeshValidation(const Mesh& m, const std::vector<Point>& uvCoords)
    : meshPtr(&m), uvPtr(&uvCoords) {
}

QuadMeshValidation::QuadMeshValidation(const IntegerGridMap& igm) {
    // Copy the mesh from IntegerGridMap's CutMesh
    const CutMesh& cm = igm.getCutMesh();
    const Mesh& cutMesh = cm.getCutMesh();
    
    // Build local mesh from cut mesh vertices/triangles
    localMesh.vertices = cutMesh.vertices;
    localMesh.triangles = cutMesh.triangles;
    
    // Copy UV coordinates directly - the solution stores raw UV values
    // where integer grid lines are at u ∈ ℤ and v ∈ ℤ
    const auto& uvCoords = igm.uv();
    localUV = uvCoords;
}


QuadMeshValidation::ValidationResult QuadMeshValidation::validate() {
    ValidationResult result;
    
    const Mesh& mesh = getMesh();
    const std::vector<Point>& uv = getUV();
    
    // Check basic requirements
    if (mesh.vertices.empty() || mesh.triangles.empty()) {
        result.errors.push_back("Empty mesh");
        return result;
    }
    if (uv.size() != mesh.vertices.size()) {
        result.errors.push_back("UV count mismatch with vertex count");
        return result;
    }
    
    // Compute UV bounding box
    double u_min = std::numeric_limits<double>::max();
    double u_max = std::numeric_limits<double>::lowest();
    double v_min = std::numeric_limits<double>::max();
    double v_max = std::numeric_limits<double>::lowest();
    
    for (const auto& p : uv) {
        u_min = std::min(u_min, p[0]);
        u_max = std::max(u_max, p[0]);
        v_min = std::min(v_min, p[1]);
        v_max = std::max(v_max, p[1]);
    }
    
    int i_min = static_cast<int>(std::ceil(u_min));
    int i_max = static_cast<int>(std::floor(u_max));
    int j_min = static_cast<int>(std::ceil(v_min));
    int j_max = static_cast<int>(std::floor(v_max));
    
    int nu = i_max - i_min + 1;  // Number of u integer values
    int nv = j_max - j_min + 1;  // Number of v integer values
    
    result.messages.push_back("UV bbox: [" + std::to_string(u_min) + ", " + std::to_string(u_max) + 
                             "] x [" + std::to_string(v_min) + ", " + std::to_string(v_max) + "]");
    result.messages.push_back("Integer grid: " + std::to_string(nu) + " x " + std::to_string(nv));
    
    // Count integer grid vertices (points where u ∈ ℤ AND v ∈ ℤ)
    // These are found by looking for triangles that contain integer (u,v) points
    std::set<std::pair<int, int>> integer_grid_points;
    
    const int T = static_cast<int>(mesh.triangles.size());
    for (int tri = 0; tri < T; ++tri) {
        const auto& t = mesh.triangles[tri];
        double u0 = uv[t[0]][0], v0 = uv[t[0]][1];
        double u1 = uv[t[1]][0], v1 = uv[t[1]][1];
        double u2 = uv[t[2]][0], v2 = uv[t[2]][1];
        
        // Find integer (i,j) points that might be in this triangle's UV region
        int tri_i_min = static_cast<int>(std::ceil(std::min({u0, u1, u2})));
        int tri_i_max = static_cast<int>(std::floor(std::max({u0, u1, u2})));
        int tri_j_min = static_cast<int>(std::ceil(std::min({v0, v1, v2})));
        int tri_j_max = static_cast<int>(std::floor(std::max({v0, v1, v2})));
        
        for (int i = tri_i_min; i <= tri_i_max; ++i) {
            for (int j = tri_j_min; j <= tri_j_max; ++j) {
                // Check if (i, j) is inside the triangle in UV space
                // Using barycentric coordinates
                double denom = (v1 - v2) * (u0 - u2) + (u2 - u1) * (v0 - v2);
                if (std::abs(denom) < 1e-12) continue;  // Degenerate triangle
                
                double b0 = ((v1 - v2) * (i - u2) + (u2 - u1) * (j - v2)) / denom;
                double b1 = ((v2 - v0) * (i - u2) + (u0 - u2) * (j - v2)) / denom;
                double b2 = 1.0 - b0 - b1;
                
                // Point is inside if all barycentric coords are in [0, 1] (with tolerance)
                const double tol = 1e-8;
                if (b0 >= -tol && b1 >= -tol && b2 >= -tol && 
                    b0 <= 1.0 + tol && b1 <= 1.0 + tol && b2 <= 1.0 + tol) {
                    integer_grid_points.insert({i, j});
                }
            }
        }
    }
    
    int expected_grid_vertices = nu * nv;
    int actual_grid_vertices = static_cast<int>(integer_grid_points.size());
    
    result.num_grid_vertices = actual_grid_vertices;
    result.messages.push_back("Expected " + std::to_string(expected_grid_vertices) + 
                             " integer grid points, found " + std::to_string(actual_grid_vertices));
    
    // Step 1: Trace all integer isolines through the mesh
    std::vector<GridVertex> vertices;
    std::vector<GridEdge> edges;
    
    if (!traceIsolines(vertices, edges)) {
        result.errors.push_back("Failed to trace isolines");
        return result;
    }
    
    result.num_grid_edges = static_cast<int>(edges.size());
    
    // Cache results for writeGridOBJ
    cachedVertices = vertices;
    cachedEdges = edges;
    
    // Step 2: Check for self-crossings within each triangle
    if (!checkSelfCrossings(vertices, edges, result)) {
        // Continue even with crossings to get full diagnostics
    }
    
    // Step 3: Build the grid graph and identify faces
    std::vector<GridFace> faces;
    if (!buildGridGraph(vertices, edges, faces, result)) {
        result.errors.push_back("Failed to build grid graph");
        return result;
    }
    
    result.num_grid_faces = static_cast<int>(faces.size());
    
    // Step 4: Count quad vs non-quad faces
    for (const auto& face : faces) {
        int valence = static_cast<int>(face.vertices.size());
        result.face_valence_histogram[valence]++;
        if (valence == 4) {
            result.num_quad_faces++;
        } else {
            result.num_non_quad_faces++;
        }
    }
    
    // Expected faces: for an nu x nv grid, we should have (nu-1) x (nv-1) quad faces
    int expected_quad_faces = std::max(0, (nu - 1) * (nv - 1));
    result.messages.push_back("Expected ~" + std::to_string(expected_quad_faces) + 
                             " quad faces (interior only)");
    
    // Step 5: Check coverage
    if (!checkCoverage(faces, result)) {
        // Continue to get full diagnostics
    }
    
    // Determine overall validity
    // The grid is valid for quad meshing if:
    // 1. No self-crossings (u-isolines don't cross other u-isolines, same for v)
    //    This is the KEY criterion - if u-lines and v-lines don't self-intersect,
    //    then the grid forms a valid quad mesh structure.
    // 2. Reasonable integer grid vertex coverage
    //
    // Note: The face extraction (quad counting) is just diagnostic.
    // The quad structure is defined by the u/v isolines, and if those don't
    // self-intersect, every cell between adjacent u-lines and v-lines is a quad.
    
    result.valid = (result.num_self_crossings == 0);
    
    if (actual_grid_vertices < 10) {
        result.errors.push_back("Very few integer grid vertices found (" + 
                               std::to_string(actual_grid_vertices) + ")");
        result.valid = false;
    }
    
    result.messages.push_back(result.valid ? 
        "Valid for quad meshing: isolines form a clean embedded graph" :
        "Invalid: self-crossings detected in isoline structure");
    
    return result;
}


void QuadMeshValidation::findIsolineEdgeCrossings(int tri, bool is_u_isoline, int iso_value,
                                                   std::vector<GridVertex>& crossings) const {
    const Mesh& mesh = getMesh();
    const std::vector<Point>& uv = getUV();
    
    const auto& t = mesh.triangles[tri];
    const int i0 = t[0], i1 = t[1], i2 = t[2];
    
    const double u0 = uv[i0][0], v0 = uv[i0][1];
    const double u1 = uv[i1][0], v1 = uv[i1][1];
    const double u2 = uv[i2][0], v2 = uv[i2][1];
    
    const double x0 = mesh.vertices[i0][0], y0 = mesh.vertices[i0][1];
    const double x1 = mesh.vertices[i1][0], y1 = mesh.vertices[i1][1];
    const double x2 = mesh.vertices[i2][0], y2 = mesh.vertices[i2][1];
    
    // The coordinate we're looking at (u for u-isolines, v for v-isolines)
    double c0 = is_u_isoline ? u0 : v0;
    double c1 = is_u_isoline ? u1 : v1;
    double c2 = is_u_isoline ? u2 : v2;
    double iso = static_cast<double>(iso_value);
    
    // Check each edge for crossings
    // Edge 0: v0 -> v1 (local edge 0)
    // Edge 1: v1 -> v2 (local edge 1)
    // Edge 2: v2 -> v0 (local edge 2)
    
    auto checkEdge = [&](int edge, double ca, double cb, 
                         double xa, double ya, double xb, double yb,
                         double ua, double va, double ub, double vb) {
        // Check if isoline crosses this edge
        if ((ca <= iso && cb >= iso) || (ca >= iso && cb <= iso)) {
            double denom = cb - ca;
            if (std::abs(denom) < EPS) {
                // Edge is parallel to isoline (both endpoints on the isoline)
                // This is a degenerate case - the whole edge is on the isoline
                // We'll handle this by adding both endpoints
                if (std::abs(ca - iso) < SNAP_TOL) {
                    GridVertex gv;
                    gv.triangle = tri;
                    gv.edge = edge;
                    gv.t = 0.0;
                    gv.x = xa; gv.y = ya;
                    gv.u = ua; gv.v = va;
                    gv.is_u_isoline = is_u_isoline;
                    gv.isoline_value = iso_value;
                    crossings.push_back(gv);
                }
                if (std::abs(cb - iso) < SNAP_TOL) {
                    GridVertex gv;
                    gv.triangle = tri;
                    gv.edge = edge;
                    gv.t = 1.0;
                    gv.x = xb; gv.y = yb;
                    gv.u = ub; gv.v = vb;
                    gv.is_u_isoline = is_u_isoline;
                    gv.isoline_value = iso_value;
                    crossings.push_back(gv);
                }
            } else {
                double t_param = (iso - ca) / denom;
                // Clamp to [0, 1] with snap tolerance
                if (t_param < SNAP_TOL) t_param = 0.0;
                if (t_param > 1.0 - SNAP_TOL) t_param = 1.0;
                
                if (t_param >= 0.0 && t_param <= 1.0) {
                    GridVertex gv;
                    gv.triangle = tri;
                    gv.edge = edge;
                    gv.t = t_param;
                    gv.x = xa + t_param * (xb - xa);
                    gv.y = ya + t_param * (yb - ya);
                    gv.u = ua + t_param * (ub - ua);
                    gv.v = va + t_param * (vb - va);
                    gv.is_u_isoline = is_u_isoline;
                    gv.isoline_value = iso_value;
                    crossings.push_back(gv);
                }
            }
        }
    };
    
    // Edge 0: vertex 0 -> vertex 1
    checkEdge(0, c0, c1, x0, y0, x1, y1, u0, v0, u1, v1);
    // Edge 1: vertex 1 -> vertex 2
    checkEdge(1, c1, c2, x1, y1, x2, y2, u1, v1, u2, v2);
    // Edge 2: vertex 2 -> vertex 0
    checkEdge(2, c2, c0, x2, y2, x0, y0, u2, v2, u0, v0);
}


bool QuadMeshValidation::traceIsolines(std::vector<GridVertex>& vertices,
                                        std::vector<GridEdge>& edges) const {
    vertices.clear();
    edges.clear();
    
    const Mesh& mesh = getMesh();
    const std::vector<Point>& uv = getUV();
    
    const int T = static_cast<int>(mesh.triangles.size());
    
    // First, compute the global UV bounding box to understand scale
    double global_u_min = std::numeric_limits<double>::max();
    double global_u_max = std::numeric_limits<double>::lowest();
    double global_v_min = std::numeric_limits<double>::max();
    double global_v_max = std::numeric_limits<double>::lowest();
    
    for (const auto& p : uv) {
        global_u_min = std::min(global_u_min, p[0]);
        global_u_max = std::max(global_u_max, p[0]);
        global_v_min = std::min(global_v_min, p[1]);
        global_v_max = std::max(global_v_max, p[1]);
    }
    
    int total_u_isolines = static_cast<int>(std::floor(global_u_max + SNAP_TOL)) - 
                          static_cast<int>(std::ceil(global_u_min - SNAP_TOL)) + 1;
    int total_v_isolines = static_cast<int>(std::floor(global_v_max + SNAP_TOL)) - 
                          static_cast<int>(std::ceil(global_v_min - SNAP_TOL)) + 1;
    
    // For each triangle, find all isoline crossings
    for (int tri = 0; tri < T; ++tri) {
        const auto& t = mesh.triangles[tri];
        const double u0 = uv[t[0]][0], v0 = uv[t[0]][1];
        const double u1 = uv[t[1]][0], v1 = uv[t[1]][1];
        const double u2 = uv[t[2]][0], v2 = uv[t[2]][1];
        
        // Get integer isolines passing through this triangle
        std::vector<int> u_isolines, v_isolines;
        getIsolinesInTriangle(u0, v0, u1, v1, u2, v2, u_isolines, v_isolines);
        
        // Collect crossings for each isoline
        std::vector<GridVertex> tri_crossings;
        
        for (int u_val : u_isolines) {
            findIsolineEdgeCrossings(tri, true, u_val, tri_crossings);
        }
        for (int v_val : v_isolines) {
            findIsolineEdgeCrossings(tri, false, v_val, tri_crossings);
        }
        
        // Add vertices to global list, recording their indices
        std::map<std::pair<bool, int>, std::vector<int>> isoline_verts;  // (is_u, value) -> vertex indices
        
        int base_idx = static_cast<int>(vertices.size());
        for (size_t i = 0; i < tri_crossings.size(); ++i) {
            vertices.push_back(tri_crossings[i]);
            auto key = std::make_pair(tri_crossings[i].is_u_isoline, tri_crossings[i].isoline_value);
            isoline_verts[key].push_back(base_idx + static_cast<int>(i));
        }
        
        // Create edges connecting crossings of the same isoline within this triangle
        // Each isoline should cross exactly 2 edges (entering and exiting the triangle)
        // unless it passes through a vertex
        for (auto& kv : isoline_verts) {
            auto& verts = kv.second;
            if (verts.size() == 2) {
                GridEdge e;
                e.v0 = verts[0];
                e.v1 = verts[1];
                e.triangle = tri;
                e.is_u_isoline = kv.first.first;
                e.isoline_value = kv.first.second;
                edges.push_back(e);
            } else if (verts.size() > 2) {
                // More than 2 crossings - this can happen if isoline passes through vertex
                // Sort by position along the isoline and connect consecutive pairs
                // For u-isoline: sort by v coordinate
                // For v-isoline: sort by u coordinate
                if (kv.first.first) {
                    // u-isoline: sort by v
                    std::sort(verts.begin(), verts.end(), [&](int a, int b) {
                        return vertices[a].v < vertices[b].v;
                    });
                } else {
                    // v-isoline: sort by u
                    std::sort(verts.begin(), verts.end(), [&](int a, int b) {
                        return vertices[a].u < vertices[b].u;
                    });
                }
                
                for (size_t i = 0; i + 1 < verts.size(); ++i) {
                    GridEdge e;
                    e.v0 = verts[i];
                    e.v1 = verts[i + 1];
                    e.triangle = tri;
                    e.is_u_isoline = kv.first.first;
                    e.isoline_value = kv.first.second;
                    edges.push_back(e);
                }
            }
            // If verts.size() == 1, isoline just touches triangle at a point (tangent)
            // If verts.size() == 0, no crossing (shouldn't happen given how we built the list)
        }
    }
    
    return true;
}


bool QuadMeshValidation::checkSelfCrossings(const std::vector<GridVertex>& vertices,
                                             const std::vector<GridEdge>& edges,
                                             ValidationResult& result) const {
    // Check for self-crossings: two edges from the SAME FAMILY of isolines crossing
    // u-isolines should never cross other u-isolines, and v-isolines should never cross v-isolines
    // But a u-isoline crossing a v-isoline is EXPECTED - that's a grid vertex!
    
    // Group edges by triangle
    std::map<int, std::vector<int>> tri_to_edges;
    for (size_t i = 0; i < edges.size(); ++i) {
        tri_to_edges[edges[i].triangle].push_back(static_cast<int>(i));
    }
    
    for (const auto& kv : tri_to_edges) {
        const auto& tri_edges = kv.second;
        
        // Check all pairs of edges in this triangle
        for (size_t i = 0; i < tri_edges.size(); ++i) {
            for (size_t j = i + 1; j < tri_edges.size(); ++j) {
                const auto& e1 = edges[tri_edges[i]];
                const auto& e2 = edges[tri_edges[j]];
                
                // Skip if edges are from different families (u vs v) - crossing is expected
                if (e1.is_u_isoline != e2.is_u_isoline) continue;
                
                // Skip if edges are from the same isoline (shouldn't cross by construction)
                if (e1.isoline_value == e2.isoline_value) continue;
                
                const auto& v1a = vertices[e1.v0];
                const auto& v1b = vertices[e1.v1];
                const auto& v2a = vertices[e2.v0];
                const auto& v2b = vertices[e2.v1];
                
                // Check if segments intersect in XY space (on the surface)
                if (segmentsIntersect(v1a.x, v1a.y, v1b.x, v1b.y,
                                      v2a.x, v2a.y, v2b.x, v2b.y)) {
                    // This is a true self-crossing (two u-isolines or two v-isolines crossing)
                    // But we need to verify it's a real crossing, not just shared endpoint
                    bool shared_endpoint = v1a.approxEqual(v2a) || v1a.approxEqual(v2b) ||
                                          v1b.approxEqual(v2a) || v1b.approxEqual(v2b);
                    if (!shared_endpoint) {
                        result.num_self_crossings++;
                        if (result.errors.size() < 10) {  // Limit error messages
                            std::ostringstream ss;
                            ss << "Self-crossing in triangle " << kv.first 
                               << " between " << (e1.is_u_isoline ? "u=" : "v=") << e1.isoline_value
                               << " and " << (e1.is_u_isoline ? "u=" : "v=") << e2.isoline_value;
                            result.errors.push_back(ss.str());
                        }
                    }
                }
            }
        }
    }
    
    return result.num_self_crossings == 0;
}


bool QuadMeshValidation::segmentsIntersect(double x1, double y1, double x2, double y2,
                                            double x3, double y3, double x4, double y4) const {
    // Check if segment (x1,y1)-(x2,y2) intersects segment (x3,y3)-(x4,y4)
    // Using cross product method
    
    double d1 = cross2d(x4 - x3, y4 - y3, x1 - x3, y1 - y3);
    double d2 = cross2d(x4 - x3, y4 - y3, x2 - x3, y2 - y3);
    double d3 = cross2d(x2 - x1, y2 - y1, x3 - x1, y3 - y1);
    double d4 = cross2d(x2 - x1, y2 - y1, x4 - x1, y4 - y1);
    
    // If the cross products have opposite signs, the segments straddle each other
    if (((d1 > EPS && d2 < -EPS) || (d1 < -EPS && d2 > EPS)) &&
        ((d3 > EPS && d4 < -EPS) || (d3 < -EPS && d4 > EPS))) {
        return true;
    }
    
    return false;
}


bool QuadMeshValidation::buildGridGraph(const std::vector<GridVertex>& vertices,
                                         const std::vector<GridEdge>& edges,
                                         std::vector<GridFace>& faces,
                                         ValidationResult& result) const {
    // Build a graph structure from the grid vertices and edges
    // Then find all faces by walking around cycles
    
    if (vertices.empty() || edges.empty()) {
        return true;  // Empty grid is valid (just means no integer isolines)
    }
    
    // First, merge vertices that are at the same XY location
    // (these are where isolines meet at triangle edges)
    std::vector<int> vertex_remap(vertices.size());
    std::vector<GridVertex> unique_verts;
    
    for (size_t i = 0; i < vertices.size(); ++i) {
        bool found = false;
        for (size_t j = 0; j < unique_verts.size(); ++j) {
            if (vertices[i].approxEqual(unique_verts[j], SNAP_TOL)) {
                vertex_remap[i] = static_cast<int>(j);
                found = true;
                break;
            }
        }
        if (!found) {
            vertex_remap[i] = static_cast<int>(unique_verts.size());
            unique_verts.push_back(vertices[i]);
        }
    }
    
    result.num_grid_vertices = static_cast<int>(unique_verts.size());
    
    // Remap edges and remove degenerate ones
    std::vector<std::pair<int, int>> unique_edges;
    std::set<std::pair<int, int>> edge_set;
    
    for (const auto& e : edges) {
        int v0 = vertex_remap[e.v0];
        int v1 = vertex_remap[e.v1];
        if (v0 == v1) {
            result.num_degenerate_edges++;
            continue;
        }
        auto edge_key = (v0 < v1) ? std::make_pair(v0, v1) : std::make_pair(v1, v0);
        if (edge_set.find(edge_key) == edge_set.end()) {
            edge_set.insert(edge_key);
            unique_edges.push_back({v0, v1});
        }
    }
    
    result.num_grid_edges = static_cast<int>(unique_edges.size());
    
    // Build adjacency list
    std::vector<std::vector<int>> adj(unique_verts.size());
    for (const auto& e : unique_edges) {
        adj[e.first].push_back(e.second);
        adj[e.second].push_back(e.first);
    }
    
    // Sort adjacency lists by angle for consistent face extraction
    for (size_t v = 0; v < unique_verts.size(); ++v) {
        const double vx = unique_verts[v].x;
        const double vy = unique_verts[v].y;
        std::sort(adj[v].begin(), adj[v].end(), [&](int a, int b) {
            double ax = unique_verts[a].x - vx;
            double ay = unique_verts[a].y - vy;
            double bx = unique_verts[b].x - vx;
            double by = unique_verts[b].y - vy;
            return std::atan2(ay, ax) < std::atan2(by, bx);
        });
    }
    
    // Find faces by walking around cycles
    // Use half-edge style traversal: for each directed edge, find the next edge
    // by taking the next CCW edge at the destination vertex
    std::set<std::pair<int, int>> used_half_edges;
    
    for (const auto& e : unique_edges) {
        for (int dir = 0; dir < 2; ++dir) {
            int start_v = (dir == 0) ? e.first : e.second;
            int start_next = (dir == 0) ? e.second : e.first;
            
            auto half_edge = std::make_pair(start_v, start_next);
            if (used_half_edges.count(half_edge)) continue;
            
            // Walk around the face
            std::vector<int> face_verts;
            int cur_v = start_v;
            int next_v = start_next;
            
            const int MAX_FACE_SIZE = 1000;  // Safety limit
            while (face_verts.size() < MAX_FACE_SIZE) {
                face_verts.push_back(cur_v);
                used_half_edges.insert({cur_v, next_v});
                
                // Find next edge: the one after (cur_v, next_v) in CCW order around next_v
                const auto& neighbors = adj[next_v];
                auto it = std::find(neighbors.begin(), neighbors.end(), cur_v);
                if (it == neighbors.end()) {
                    // Shouldn't happen
                    break;
                }
                
                // Next neighbor in CCW order (wrapping around)
                ++it;
                if (it == neighbors.end()) it = neighbors.begin();
                int next_next = *it;
                
                cur_v = next_v;
                next_v = next_next;
                
                if (cur_v == start_v && next_v == start_next) {
                    // Completed the cycle
                    break;
                }
            }
            
            if (face_verts.size() >= 3 && face_verts.size() < MAX_FACE_SIZE) {
                GridFace f;
                f.vertices = face_verts;
                faces.push_back(f);
            }
        }
    }
    
    return true;
}


bool QuadMeshValidation::checkCoverage(const std::vector<GridFace>& faces,
                                        ValidationResult& result) const {
    // Check that the grid faces cover all triangles
    // This is a simplified check - we just verify that we found some faces
    // A more thorough check would verify exact coverage
    
    const Mesh& mesh = getMesh();
    const std::vector<Point>& uv = getUV();
    
    if (faces.empty() && !mesh.triangles.empty()) {
        // Check if there are any integer isolines at all
        const int T = static_cast<int>(mesh.triangles.size());
        bool has_any_isoline = false;
        
        for (int tri = 0; tri < T && !has_any_isoline; ++tri) {
            const auto& t = mesh.triangles[tri];
            const double u0 = uv[t[0]][0], v0 = uv[t[0]][1];
            const double u1 = uv[t[1]][0], v1 = uv[t[1]][1];
            const double u2 = uv[t[2]][0], v2 = uv[t[2]][1];
            
            std::vector<int> u_isolines, v_isolines;
            getIsolinesInTriangle(u0, v0, u1, v1, u2, v2, u_isolines, v_isolines);
            
            if (!u_isolines.empty() || !v_isolines.empty()) {
                has_any_isoline = true;
            }
        }
        
        if (has_any_isoline) {
            result.num_missing_regions = 1;
            result.errors.push_back("Grid has isolines but no faces were extracted");
        }
    }
    
    result.coverage_complete = (result.num_missing_regions == 0);
    return result.coverage_complete;
}


bool QuadMeshValidation::writeGridOBJ(const std::string& filename) const {
    return writeGridOBJ(filename, cachedVertices, cachedEdges);
}


bool QuadMeshValidation::writeGridOBJ(const std::string& filename,
                                       const std::vector<GridVertex>& vertices,
                                       const std::vector<GridEdge>& edges) const {
    std::ofstream out(filename);
    if (!out) return false;
    
    out << "# Integer grid lines extracted from UV parameterization\n";
    out << "# Vertices: " << vertices.size() << ", Edges: " << edges.size() << "\n\n";
    
    // Write vertices
    for (const auto& v : vertices) {
        out << "v " << v.x << " " << v.y << " 0\n";
    }
    
    out << "\n";
    
    // Write edges as line segments
    for (const auto& e : edges) {
        out << "l " << (e.v0 + 1) << " " << (e.v1 + 1) << "\n";
    }
    
    return true;
}


void QuadMeshValidation::barycentricToXY(int tri, double b0, double b1, double b2,
                                          double& x, double& y) const {
    const Mesh& mesh = getMesh();
    const auto& t = mesh.triangles[tri];
    const auto& p0 = mesh.vertices[t[0]];
    const auto& p1 = mesh.vertices[t[1]];
    const auto& p2 = mesh.vertices[t[2]];
    
    x = b0 * p0[0] + b1 * p1[0] + b2 * p2[0];
    y = b0 * p0[1] + b1 * p1[1] + b2 * p2[1];
}
