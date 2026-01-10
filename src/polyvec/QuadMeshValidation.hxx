#ifndef __QUAD_MESH_VALIDATION_HXX__
#define __QUAD_MESH_VALIDATION_HXX__

#include <vector>
#include <array>
#include <set>
#include <map>
#include <string>

#include "Mesh.hxx"

// Forward declaration
class IntegerGridMap;

// Validate whether the UV parameterization produces a valid quad mesh
// by tracing integer isolines (u ∈ ℤ and v ∈ ℤ) on the cut mesh.
//
// A valid quad mesh requires:
// 1. No self-crossings of isolines on the surface
// 2. All faces of the embedded grid are quads (4 edges)
// 3. No missing regions (grid covers the entire surface)
//
// This is the key test for whether the parameterization will work for quad meshing.

class QuadMeshValidation {
public:
    struct ValidationResult {
        bool valid = false;
        
        // Grid statistics
        int num_grid_vertices = 0;    // Number of integer grid intersections
        int num_grid_edges = 0;       // Number of grid edges
        int num_grid_faces = 0;       // Number of grid faces
        int num_quad_faces = 0;       // Number of 4-sided faces (should equal num_grid_faces)
        int num_non_quad_faces = 0;   // Number of non-quad faces (should be 0)
        
        // Error counts
        int num_self_crossings = 0;   // Number of places where isolines cross incorrectly
        int num_degenerate_edges = 0; // Number of zero-length or collapsed edges
        int num_missing_regions = 0;  // Number of triangles not covered by the grid
        bool coverage_complete = false;
        
        // Face valence distribution (for debugging)
        std::map<int, int> face_valence_histogram;  // valence -> count
        
        // Detailed error/info messages
        std::vector<std::string> errors;
        std::vector<std::string> messages;
        
        std::string summary() const;
    };
    
    // A point where an integer isoline crosses a triangle edge
    struct GridVertex {
        int triangle = -1;        // Triangle containing this vertex
        int edge = -1;            // Edge of the triangle (0,1,2) or -1 if interior
        double t = 0.0;           // Parameter along edge [0,1], or barycentric coords if interior
        double x = 0.0, y = 0.0;  // XY position on the surface
        double u = 0.0, v = 0.0;  // UV coordinates (one should be integer)
        bool is_u_isoline = true; // true if this is on a u=const line, false for v=const
        int isoline_value = 0;    // The integer value of the constant coordinate
        
        // For identifying unique vertices
        bool operator<(const GridVertex& other) const;
        bool approxEqual(const GridVertex& other, double tol = 1e-10) const;
    };
    
    // An edge of the integer grid
    struct GridEdge {
        int v0 = -1, v1 = -1;     // Indices into gridVertices
        int triangle = -1;        // Triangle this edge segment is in
        bool is_u_isoline = true; // Direction of this edge
        int isoline_value = 0;    // Integer value of the constant coordinate
    };
    
    // A face of the integer grid
    struct GridFace {
        std::vector<int> vertices;  // Vertex indices in order around the face
        int u_min = 0, u_max = 0;   // UV bounds
        int v_min = 0, v_max = 0;
        std::vector<int> triangles; // Triangles that contribute to this face
    };
    
    // Construct validator with mesh and UV coordinates
    QuadMeshValidation(const Mesh& mesh, const std::vector<Point>& uvCoords);
    
    // Construct validator from IntegerGridMap (convenience)
    QuadMeshValidation(const IntegerGridMap& igm);
    
    // Run full validation
    ValidationResult validate();
    
    // Individual validation steps (for debugging)
    bool traceIsolines(std::vector<GridVertex>& vertices,
                       std::vector<GridEdge>& edges) const;
    
    bool buildGridGraph(const std::vector<GridVertex>& vertices,
                        const std::vector<GridEdge>& edges,
                        std::vector<GridFace>& faces,
                        ValidationResult& result) const;
    
    bool checkSelfCrossings(const std::vector<GridVertex>& vertices,
                            const std::vector<GridEdge>& edges,
                            ValidationResult& result) const;
    
    bool checkCoverage(const std::vector<GridFace>& faces,
                       ValidationResult& result) const;

    // Export the grid for visualization (uses cached vertices/edges from validate())
    bool writeGridOBJ(const std::string& filename) const;
    
    // Export the grid for visualization with explicit vertices/edges
    bool writeGridOBJ(const std::string& filename,
                      const std::vector<GridVertex>& vertices,
                      const std::vector<GridEdge>& edges) const;

private:
    const Mesh* meshPtr = nullptr;
    const std::vector<Point>* uvPtr = nullptr;
    
    // Local copies when constructed from IntegerGridMap
    Mesh localMesh;
    std::vector<Point> localUV;
    
    // Cached results from validate()
    mutable std::vector<GridVertex> cachedVertices;
    mutable std::vector<GridEdge> cachedEdges;
    
    // Accessors
    const Mesh& getMesh() const { return meshPtr ? *meshPtr : localMesh; }
    const std::vector<Point>& getUV() const { return uvPtr ? *uvPtr : localUV; }
    
    // Find where isoline u=u_val or v=v_val crosses triangle edges
    void findIsolineEdgeCrossings(int tri, bool is_u_isoline, int iso_value,
                                  std::vector<GridVertex>& crossings) const;
    
    // Compute XY position from barycentric coordinates in a triangle
    void barycentricToXY(int tri, double b0, double b1, double b2,
                         double& x, double& y) const;
    
    // Check if two line segments intersect (for self-crossing detection)
    bool segmentsIntersect(double x1, double y1, double x2, double y2,
                           double x3, double y3, double x4, double y4) const;
};

#endif // __QUAD_MESH_VALIDATION_HXX__
