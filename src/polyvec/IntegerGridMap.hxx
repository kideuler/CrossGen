#ifndef __INTEGER_GRID_MAP_HXX__
#define __INTEGER_GRID_MAP_HXX__

#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "CutMesh.hxx"

// Construct an integer grid map (a seamless parameterization suitable for
// integer-grid extraction) from a CutMesh.
//
// This implements the parameterization stage described in:
//   Bommes et al., "Mixed-Integer Quadrangulation" (MIQ)
//
// Notes for this codebase:
//  * The input mesh is already planar (2D) and CutMesh produces a topological
//    disk by duplicating vertices along cut edges.
//  * We build a per-triangle oriented frame by propagating a consistent
//    orientation of the 4-RoSy field on the cut mesh (BFS).
//  * We solve for (u,v) at cut-mesh vertices by minimizing the MIQ orientation
//    energy, with optional local stiffening to reduce fold-overs.
//  * Seamless gluing across the cut boundary is enforced by MIQ-style
//    transition functions with integer translations and 90-degree rotations.

class IntegerGridMap {
public:
    struct Options {
        Options() = default;
        // Target edge length in parameter space (MIQ uses h).
        // If h <= 0, it will be auto-computed based on mesh size to give
        // approximately `target_grid_lines` integer grid lines across the mesh.
        double h = 0.0;  // 0 means auto-compute
        int target_grid_lines = 100;  // Number of integer grid lines across mesh (when h=0)

        // Local stiffening (MIQ Sec. 4): repeat solve -> update weights.
        int stiffening_iterations = 10;
        int stiffening_smooth_steps = 3;
        double stiffening_c = 0.5;   // multiplier for |4 Δλ|
        double stiffening_d = 5.0;   // clamp

        // Seam constraint mode:
        // - soft_seam_weight > 0: Use soft constraints (penalty-based), allows slight
        //   violation to prevent flipped triangles. Lower values give more flexibility.
        //   Recommended: 0.1 for robustness, 1.0 for tighter gluing.
        // - soft_seam_weight <= 0: Use hard constraints (exact gluing), may cause flips
        //   near singularities.
        double soft_seam_weight = 0.1;  // Default: soft constraints for robustness

        // Barrier-based local injectivity enforcement:
        // Uses a log-barrier term to prevent triangle flips during optimization.
        // The energy becomes: E_total = E_dirichlet + barrier_weight * sum(-log(det(J)))
        // Higher values enforce injectivity more strongly but may increase distortion.
        bool use_barrier = true;
        double barrier_weight = 1.0;       // Weight for barrier term
        int barrier_iterations = 20;       // Max iterations for barrier optimization
        double barrier_eps = 1e-6;         // Minimum allowed triangle area ratio

        // Mixed-integer rounding.
        // For performance, we round in batches (smallest errors first), rather
        // than strictly one-by-one.
        bool do_rounding = true;
        int max_round_iterations = 50;
        int round_batch_size = 200;
        double round_keep_threshold = 1e-4; // fix all vars with error <= threshold

        // Numerical regularization for translation variables (helps conditioning).
        double translation_regularization = 1e-10;
    };

    struct Stats {
        int flipped_triangles = 0;
        double max_lambda = 0.0;
        double mean_lambda = 0.0;
    };

    explicit IntegerGridMap(const CutMesh &cutMesh);

    // Compute the integer grid map. On success, uv() is populated.
    bool solve(const Options &opt);
    // Convenience overload with default options.
    bool solve();

    // Access the cut mesh used for the solve.
    const Mesh& getMesh() const { return cut; }
    
    // Access the underlying CutMesh.
    const CutMesh& getCutMesh() const { return cm; }
    
    // Access the raw solution vector (for validation).
    const Eigen::VectorXd& getSolution() const { return solution; }

    // The effective h value used in the last solve() call.
    double getEffectiveH() const { return effectiveH; }

    // UV coordinates per cut-mesh vertex (same indexing as cut.vertices).
    const std::vector<Point>& uv() const { return uvCoords; }

    // Write the cut mesh with UVs as an OBJ (vt per vertex, faces as v/vt).
    bool writeOBJWithUV(const std::string &filename) const;
    
    // Write the ORIGINAL (uncut) mesh with UVs as texture coordinates.
    // Maps UV coordinates from cut mesh back to original mesh vertices.
    // Each triangle corner gets its own texture coordinate to handle seams correctly.
    bool writeOriginalMeshWithUV(const std::string &filename) const;

    // Compute simple distortion/flip stats for the current uv solution.
    // If h<=0, uses the effective h from last solve.
    Stats computeStats(double h = 0.0) const;

private:
    const CutMesh &cm;
    Mesh cut;

    // Effective h value used in the last solve() call.
    double effectiveH = 1.0;
    
    // Raw solution vector from solve() - used for validation
    Eigen::VectorXd solution;

    // Per-triangle u direction from CutMesh (original triangle ordering).
    std::vector<Point> uField;

    // Oriented frame selection index per triangle (0..3): uAxis = Rot90^k * uField.
    std::vector<int> triK;

    struct BoundaryEdge {
        int tri = -1;
        int localEdge = -1;
        int v0 = -1; // cut vertex id
        int v1 = -1; // cut vertex id
        int o0 = -1; // original vertex id
        int o1 = -1; // original vertex id
    };

    struct Seam {
        int v0a = -1, v0b = -1; // endpoints on "side 0" (cut vertices)
        int v1a = -1, v1b = -1; // endpoints on "side 1" (cut vertices)
        int tri0 = -1, tri1 = -1; // incident triangles
        int rot = 0;              // 0..3, rotation in transition function
        int jVar = -1;            // variable index in the global unknown vector
        int kVar = -1;            // variable index in the global unknown vector
    };

    std::vector<Seam> seams;

    // Chosen cut-mesh vertex id per original singularity (for snapping).
    std::vector<int> singularCutVerts;

    // Final UV per cut-mesh vertex.
    std::vector<Point> uvCoords;

    // Helpers
    static Eigen::Vector2d rot90Vec(const Eigen::Vector2d &v, int k);
    static Eigen::Matrix2d rot90Mat(int k);

    void buildTriangleOrientations();
    void buildSeams();
    void chooseSingularityVertices();

    // Assemble quadratic energy: 0.5 x^T H x - g^T x.
    void assembleEnergy(double h,
                        const std::vector<double> &weights,
                        double translationReg,
                        double softSeamWeight,
                        Eigen::SparseMatrix<double> &H,
                        Eigen::VectorXd &g) const;

    // Build base constraints (seams + anchor), plus any additional fixed-variable
    // constraints of the form x[idx] = value.
    // If softSeamWeight > 0, seam constraints are omitted (handled in energy).
    void assembleConstraints(const std::vector<std::pair<int,double>> &fixed,
                             double softSeamWeight,
                             Eigen::SparseMatrix<double> &C,
                             Eigen::VectorXd &d) const;

    bool solveKKT(const Eigen::SparseMatrix<double> &H,
                  const Eigen::VectorXd &g,
                  const Eigen::SparseMatrix<double> &C,
                  const Eigen::VectorXd &d,
                  Eigen::VectorXd &x) const;

    Stats computeStatsFromX(const Eigen::VectorXd &x, double h, bool verbose = false) const;
    void updateStiffeningWeights(const Eigen::VectorXd &x,
                                 double h,
                                 double c,
                                 double d,
                                 int smoothSteps,
                                 std::vector<double> &weights) const;

    // Untangle flipped triangles using gradient descent with barrier
    bool untangleUV(Eigen::VectorXd &x,
                    const std::vector<std::pair<int,double>> &fixed,
                    double h,
                    int maxIter) const;

    // Barrier-augmented optimization for local injectivity.
    // Uses iterative reweighted least squares with a log-barrier term.
    bool solveWithBarrier(const Options &opt,
                          double h,
                          const std::vector<std::pair<int,double>> &fixed,
                          std::vector<double> &weights,
                          Eigen::VectorXd &x);
    
    // Compute per-triangle Jacobian determinant from current solution.
    void computeTriangleJacobians(const Eigen::VectorXd &x,
                                  std::vector<double> &detJ,
                                  std::vector<double> &areas) const;

    void extractUVFromX(const Eigen::VectorXd &x);
};

#endif // __INTEGER_GRID_MAP_HXX__
