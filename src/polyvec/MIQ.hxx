#ifndef __MIQ_HXX__
#define __MIQ_HXX__

#include <complex>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <CoMISo/Solver/ConstrainedSolver.hh>
#include <CoMISo/Solver/MISolver.hh>

#include "CutMesh.hxx"

// Mixed-Integer Quadrangulation (MIQ)
//
// Implements the algorithm from:
//   "Mixed-Integer Quadrangulation" by D. Bommes, H. Zimmer, L. Kobbelt
//   ACM SIGGRAPH 2009
//
// Given a cut mesh (topological disk) with a combed cross-field, computes
// a seamless global parametrization whose gradient aligns with the field.
// The parametrization is represented as per-vertex UV coordinates on the
// cut mesh.
//
// Usage:
//   CutMesh cutMesh(polyField);
//   MIQSolver miq(cutMesh);
//   miq.solve();
//   const auto& UV = miq.getUV();       // per cut-mesh vertex UVs
//   const auto& FUV = miq.getFUV();     // face UV indices (same as cut mesh faces)
//
class MIQSolver {
public:
    // Seam edge information for constraints.
    // On each seam edge, we have matching vertices across the seam that
    // must satisfy uv0' = R(mismatch) * uv0 + t  (with integer t).
    struct SeamInfo {
        int v0 = -1;       // cut-mesh vertex on one side
        int v0p = -1;      // corresponding cut-mesh vertex on the other side
        int mismatch = 0;  // rotation index (0,1,2,3) for pi/2 increments
        int integerVar = 0; // index of the integer translation variable for this seam
    };

    // System dimensions and seam info.
    struct MeshSystemInfo {
        int numVertVariables = 0;   // number of cut-mesh vertices
        int numIntegerCuts = 0;     // number of integer translation variables
        std::vector<SeamInfo> edgeSeamInfo;
    };

    explicit MIQSolver(const CutMesh& cutMesh);

    // Solve the MIQ parametrization problem.
    // @param gradientSize         Global scaling for the UV gradient (controls quad resolution).
    //                             Larger values -> finer quads.
    // @param stiffness            Weight for stiffness iterations (max delta per iteration).
    // @param directRound          If true, greedily round all integers at once (faster but lower quality).
    // @param iter                 Number of stiffness iterations (0 = no stiffening).
    // @param localIter            Number of local iterations for integer rounding.
    // @param doRound              Enable integer rounding (set false for debugging).
    // @param singularityRound     If true, round singularity coordinates to integers.
    // @param boundaryFeatures     If true, add boundary edges as hard feature constraints.
    //                             For curved 2D boundaries, setting this to false may reduce flips.
    void solve(
        double gradientSize = 30.0,
        double stiffness = 5.0,
        bool directRound = false,
        unsigned int iter = 5,
        unsigned int localIter = 5000,
        bool doRound = true,
        bool singularityRound = true,
        bool boundaryFeatures = true
    );

    // Per cut-mesh vertex UV coordinates (Vcut x 2).
    const Eigen::MatrixXd& getUV() const { return UV; }

    // Face indices into UV (same as cut mesh triangles, nT x 3).
    const Eigen::MatrixXi& getFUV() const { return FUV; }

    // Per-wedge UV coordinates (nT x 6), storing [u0,v0, u1,v1, u2,v2] per triangle.
    const Eigen::MatrixXd& getWUV() const { return WUV; }

    // Number of flipped triangles in the UV parametrization.
    int numFlips() const;

    // Write UV mesh to an OBJ file (z=0).
    bool writeUVMesh(const std::string& filename) const;

    // Write the original (uncut) mesh with UV texture coordinates.
    // This outputs: v (3D positions), vt (UV coords per wedge), f (with v/vt indices).
    bool writeTexturedMesh(const std::string& filename) const;

    // Extract coarse quad mesh by tracing integer isolines.
    // Returns true if successful, false otherwise.
    // The quad mesh is stored in quadVertices_ and quadFaces_.
    bool extractQuadMesh();

    // Write the extracted quad mesh to an OBJ file.
    bool writeQuadMesh(const std::string& filename) const;

    // Get the extracted quad mesh vertices (Nx2, in original mesh coordinates).
    const Eigen::MatrixXd& getQuadVertices() const { return quadVertices_; }

    // Get the extracted quad mesh faces (Mx4, vertex indices).
    const Eigen::MatrixXi& getQuadFaces() const { return quadFaces_; }

private:
    // Input references
    const CutMesh& cutMesh_;

    // Converted Eigen data
    Eigen::MatrixXd V_;       // original mesh vertices (#V x 2)
    Eigen::MatrixXi F_;       // original mesh faces (#F x 3)
    Eigen::MatrixXd Vcut_;    // cut mesh vertices (#Vcut x 2)
    Eigen::MatrixXi Fcut_;    // cut mesh faces (#Fcut x 3) - same connectivity, new indices
    Eigen::MatrixXi TT_;      // triangle-triangle adjacency on cut mesh (#F x 3)
    Eigen::MatrixXi TTi_;     // local edge indices for TT (#F x 3)
    Eigen::MatrixXd PD1_;     // combed u field per face (#F x 2)
    Eigen::MatrixXd PD2_;     // combed v field per face (#F x 2)

    // Mismatch across each edge (in terms of pi/2 rotations, 0-3).
    Eigen::MatrixXi mismatch_;  // #F x 3

    // Per-original-vertex singularity flag.
    Eigen::VectorXi singular_;  // #Vorig x 1 (1 if singular, 0 otherwise)

    // Seam indicator per half-edge. seams_(f, e) != 0 if edge e of face f is a seam.
    Eigen::MatrixXi seams_;     // #F x 3

    // Stiffness per face (for iterative stiffening).
    Eigen::VectorXd stiffnessVector_;

    // Maximum stiffness delta per iteration (the 'd' parameter in stiffening).
    // This is controlled by the 'stiffness' parameter passed to solve().
    double stiffnessWeight_ = 5.0;

    // Whether to add boundary edges as hard feature constraints.
    bool boundaryFeatures_ = true;

    // System info for Poisson solver.
    MeshSystemInfo systemInfo_;

    // Output UV coordinates
    Eigen::MatrixXd UV;       // per cut-mesh vertex, #Vcut x 2
    Eigen::MatrixXi FUV;      // face indices (same as Fcut_)
    Eigen::MatrixXd WUV;      // per-wedge UV, #F x 6

    // ------------------------------------------
    // Setup / preprocessing
    // ------------------------------------------
    void convertMeshToEigen();
    void computeTriangleTriangleAdjacency();
    void computeMismatch();
    void computeSeams();
    void initVertexIndexing();

    // ------------------------------------------
    // Poisson solving
    // ------------------------------------------
    void solvePoisson(
        double gradientSize,
        double gridResolution,
        bool directRound,
        unsigned int localIter,
        bool doRound,
        bool singularityRound
    );

    void buildLaplacianMatrix(double vfscale);
    void buildSeamConstraints();
    void addBoundaryFeatureConstraints();
    void fixBlockedVertex();
    void addSingularityRound();
    void mixedIntegerSolve(double coneGridRes, bool directRound, unsigned int localIter);
    void mapCoords();

    // ------------------------------------------
    // Stiffening for flip reduction
    // ------------------------------------------
    double distortion(int f, double h) const;
    double laplaceDistortion(int f, double h) const;
    bool updateStiffening(double grad_size);
    bool isFlipped(int f) const;
    static bool isFlipped(const Eigen::Vector2d& uv0, const Eigen::Vector2d& uv1, const Eigen::Vector2d& uv2);

    // ------------------------------------------
    // Helper functions
    // ------------------------------------------
    std::complex<double> getRotationComplex(int interval) const;
    int getFirstVertexIndex(int origV) const;

    // Linear system components
    Eigen::SparseMatrix<double> Lhs_;
    Eigen::SparseMatrix<double> Constraints_;
    Eigen::VectorXd rhs_;
    Eigen::VectorXd constraintsRhs_;
    std::vector<double> X_;
    std::vector<int> idsToRound_;
    std::vector<int> hardConstraints_;

    // Sizes
    int nVertVars_ = 0;
    int nIntegerVars_ = 0;
    int numCutConstraint_ = 0;
    int nFixedVars_ = 0;
    int numTotalVars_ = 0;
    int numConstraintEquations_ = 0;
    int numBoundaryConstraints_ = 0;
    bool integerRounding_ = true;

    // Vertex-triangle adjacency
    std::vector<std::vector<int>> VF_;  // vertex -> list of incident faces
    std::vector<std::vector<int>> VFi_; // vertex -> local index in each face

    // Extracted quad mesh
    Eigen::MatrixXd quadVertices_;  // quad mesh vertices in original coordinates
    Eigen::MatrixXi quadFaces_;     // quad mesh faces (Nx4)
};

#endif // __MIQ_HXX__