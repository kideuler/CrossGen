#include "UVGParam.hxx"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <Eigen/SparseLU>

// ---------------------------------------------------------------------------
// UVGParam – Global UV parameterization from the .tex derivation
//
// Assembles the weak-form Poisson system
//
//   L u = b^X,   L v = b^Y
//
// where
//   L_{ij} = sum_T  A_T  (grad phi_j · grad phi_i)      (cotangent Laplacian)
//   b^X_i  = sum_T  A_T  (X_T · grad phi_i)
//   b^Y_i  = sum_T  A_T  (Y_T · grad phi_i)
//
// and X_T, Y_T are the combed per-triangle u/v field vectors from CutMesh.
// One vertex (index 0) is pinned (Dirichlet) to remove the translational d.o.f.
// ---------------------------------------------------------------------------

UVGParam::UVGParam(const CutMesh& cutMesh)
    : cm_(cutMesh)
{
    const Mesh& mesh = cutMesh.getCutMesh();
    const int nV = static_cast<int>(mesh.vertices.size());
    const int nT = static_cast<int>(mesh.triangles.size());

    if (nV == 0 || nT == 0) {
        throw std::runtime_error("UVGParam: cut mesh is empty.");
    }

    const std::vector<Point>& xField = cutMesh.getUField(); // X_T  (u direction)
    const std::vector<Point>& yField = cutMesh.getVField(); // Y_T  (v direction)

    if (static_cast<int>(xField.size()) != nT ||
        static_cast<int>(yField.size()) != nT) {
        throw std::runtime_error("UVGParam: field size does not match triangle count.");
    }

    // ------------------------------------------------------------------
    // Assemble  L  and  b^X, b^Y  (Algorithm from the .tex file)
    // ------------------------------------------------------------------
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<size_t>(nT) * 9);

    Eigen::VectorXd bx = Eigen::VectorXd::Zero(nV);
    Eigen::VectorXd by = Eigen::VectorXd::Zero(nV);

    for (int t = 0; t < nT; ++t) {
        const Triangle& tri = mesh.triangles[t];
        const int i = tri[0], j = tri[1], k = tri[2];
        const int idx[3] = {i, j, k};

        const Point& p0 = mesh.vertices[i];
        const Point& p1 = mesh.vertices[j];
        const Point& p2 = mesh.vertices[k];

        // Signed area  (positive for CCW).
        double A2 = (p1[0]-p0[0])*(p2[1]-p0[1]) - (p2[0]-p0[0])*(p1[1]-p0[1]);
        double A  = 0.5 * A2;

        if (std::abs(A) < 1e-15) continue; // degenerate triangle, skip

        // Gradients of the three hat functions (constant per triangle).
        //   grad phi_0 = ( y1-y2, x2-x1 ) / (2A)
        //   grad phi_1 = ( y2-y0, x0-x2 ) / (2A)
        //   grad phi_2 = ( y0-y1, x1-x0 ) / (2A)
        double invA2 = 1.0 / A2;
        Point g[3];
        g[0] = { (p1[1]-p2[1])*invA2,  (p2[0]-p1[0])*invA2 };
        g[1] = { (p2[1]-p0[1])*invA2,  (p0[0]-p2[0])*invA2 };
        g[2] = { (p0[1]-p1[1])*invA2,  (p1[0]-p0[0])*invA2 };

        const Point& X = xField[t];
        const Point& Y = yField[t];

        // Accumulate Laplacian entries and load vector (tex algorithm, lines 9-16).
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                double w = A * (g[a][0]*g[b][0] + g[a][1]*g[b][1]);
                triplets.emplace_back(idx[a], idx[b], w);
            }
            double dotX = X[0]*g[a][0] + X[1]*g[a][1];
            double dotY = Y[0]*g[a][0] + Y[1]*g[a][1];
            bx(idx[a]) += A * dotX;
            by(idx[a]) += A * dotY;
        }
    }

    // ------------------------------------------------------------------
    // Build sparse matrix and pin vertex 0 (Dirichlet BC).
    // ------------------------------------------------------------------
    Eigen::SparseMatrix<double> L(nV, nV);
    L.setFromTriplets(triplets.begin(), triplets.end());

    // Row 0: identity row (L[0,j] = delta_{0j})
    // Zero out row 0 and column 0 contributions, set diagonal to 1.
    for (Eigen::SparseMatrix<double>::InnerIterator it(L, 0); it; ++it)
        it.valueRef() = (it.row() == 0 && it.col() == 0) ? 1.0 : 0.0;
    // Also zero the entries in column 0 (other rows) to keep symmetry.
    for (int col = 0; col < nV; ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(L, col); it; ++it) {
            if (it.row() == 0 && it.col() != 0)
                it.valueRef() = 0.0;
        }
    }
    bx(0) = 0.0;
    by(0) = 0.0;

    L.makeCompressed();

    // ------------------------------------------------------------------
    // Solve with SparseLU (SPD after pinning).
    // ------------------------------------------------------------------
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(L);
    solver.factorize(L);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("UVGParam: factorization of Laplacian failed.");
    }

    u_ = solver.solve(bx);
    v_ = solver.solve(by);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("UVGParam: linear solve failed.");
    }

    std::cout << "[UVGParam] Solved UV parameterization: "
              << nV << " vertices, " << nT << " triangles.\n";
}

// ---------------------------------------------------------------------------
bool UVGParam::writeOBJ(const std::string& filename) const
{
    const Mesh& mesh = cm_.getCutMesh();
    const int nV = static_cast<int>(mesh.vertices.size());
    const int nT = static_cast<int>(mesh.triangles.size());

    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "[UVGParam] Cannot open file: " << filename << "\n";
        return false;
    }

    ofs << "# UVGParam output\n";

    // Geometric vertices (z=0 from the cut mesh 2D positions)
    for (int v = 0; v < nV; ++v) {
        const Point& p = mesh.vertices[v];
        ofs << "v " << p[0] << " " << p[1] << " 0\n";
    }

    // UV texture coordinates
    for (int v = 0; v < nV; ++v) {
        ofs << "vt " << u_(v) << " " << v_(v) << "\n";
    }

    // Faces with UV indices (1-based in OBJ)
    for (int t = 0; t < nT; ++t) {
        const Triangle& tri = mesh.triangles[t];
        ofs << "f "
            << (tri[0]+1) << "/" << (tri[0]+1) << " "
            << (tri[1]+1) << "/" << (tri[1]+1) << " "
            << (tri[2]+1) << "/" << (tri[2]+1) << "\n";
    }

    std::cout << "[UVGParam] Written UV OBJ: " << filename << "\n";
    return true;
}
