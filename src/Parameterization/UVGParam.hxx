#ifndef __UVGPARAM_HXX__
#define __UVGPARAM_HXX__

#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "CutMesh.hxx"

// Global UV parameterization via a face-based Poisson system.
//
// Given a CutMesh (topological disk) with combed per-triangle cross-field
// directions (uField, vField), solves:
//
//   L u = b^X,   L v = b^Y
//
// where L is the cotangent Laplacian assembled over the cut mesh and
// b^X, b^Y are the field-alignment load vectors.  One vertex is pinned
// to remove the translational null space.
//
// After construction, per-vertex (u,v) coordinates are available via
// getU() / getV().  A write helper exports a UV-OBJ file.
class UVGParam {
public:
    // Solve the Poisson system on the cut mesh.
    explicit UVGParam(const CutMesh& cutMesh);

    // Per-vertex u and v coordinates (indexed over cut-mesh vertices).
    const Eigen::VectorXd& getU() const { return u_; }
    const Eigen::VectorXd& getV() const { return v_; }

    // Access the underlying cut mesh (for rendering).
    const CutMesh& getCutMesh() const { return cm_; }

    // Returns the number of triangles with negative signed area in UV space
    // (i.e. flipped/inverted triangles). Zero means a valid, flip-free parametrization.
    int numFlippedTriangles() const;

    // Write the parametrized mesh as an OBJ with UV coordinates.
    bool writeOBJ(const std::string& filename) const;

private:
    const CutMesh& cm_;
    Eigen::VectorXd u_;
    Eigen::VectorXd v_;
    Eigen::VectorXd tu_; // translation vector for seam transitions (same size as u_)  
    Eigen::VectorXd tv_; // translation vector for seam transitions (same size as v_)
};

#endif // __UVGPARAM_HXX__