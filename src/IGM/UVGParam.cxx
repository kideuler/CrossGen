#include "UVGParam.hxx"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <Eigen/SparseLU>

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
    // Assemble  A  and  [b^X, b^Y]
    // ------------------------------------------------------------------
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<size_t>(2*nT) * 9);

    Eigen::VectorXd b = Eigen::VectorXd::Zero(2*nV);

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
                triplets.emplace_back(idx[a]+nV, idx[b]+nV, w);
            }
            double dotX = X[0]*g[a][0] + X[1]*g[a][1];
            double dotY = Y[0]*g[a][0] + Y[1]*g[a][1];
            b(idx[a]) += A * dotX;
            b(idx[a]+nV) += A * dotY;
        }
    }

    // ------------------------------------------------------------------
    // Build sparse matrix and pin vertex 0 (Dirichlet BC).
    // ------------------------------------------------------------------
    Eigen::SparseMatrix<double> A(2*nV, 2*nV);
    A.setFromTriplets(triplets.begin(), triplets.end());

    // Row 0: identity row (A[0,j] = delta_{0j})
    // Zero out row 0 and column 0 contributions, set diagonal to 1.
    for (Eigen::SparseMatrix<double>::InnerIterator it(A, 0); it; ++it)
        it.valueRef() = (it.row() == 0 && it.col() == 0) ? 1.0 : 0.0;

    // Row nV: identity row (A[nV,j] = delta_{nVj})
    for (Eigen::SparseMatrix<double>::InnerIterator it(A, nV); it; ++it)
        it.valueRef() = (it.row() == nV && it.col() == nV) ? 1.0 : 0.0;

    // Also zero the entries in column 0 (other rows) to keep symmetry.
    for (int col = 0; col < nV; ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, col); it; ++it) {
            if (it.row() == 0 && it.col() != 0)
                it.valueRef() = 0.0;
        }

        for (Eigen::SparseMatrix<double>::InnerIterator it(A, col+nV); it; ++it) {
            if (it.row() == nV && it.col() != nV)
                it.valueRef() = 0.0;
        }
    }
    b(0) = 0.0;
    b(nV) = 0.0;


    // Start building constraint matrix C and constraint RHS d.
    // isoline

    A.makeCompressed();

    // ------------------------------------------------------------------
    // Solve with SparseLU (SPD after pinning).
    // ------------------------------------------------------------------
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("UVGParam: factorization of Laplacian failed.");
    }
    auto x = solver.solve(b);
    u_ = x.segment(0, nV);
    v_ = x.segment(nV, nV);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("UVGParam: linear solve failed.");
    }

    std::cerr << "[UVGParam] Solved UV parameterization: "
              << nV << " vertices, " << nT << " triangles.\n";
}

// ---------------------------------------------------------------------------
int UVGParam::numFlippedTriangles() const
{
    const Mesh& mesh = cm_.getCutMesh();
    const int nV = static_cast<int>(u_.size());
    int count = 0;

    for (const Triangle& tri : mesh.triangles) {
        const int i = tri[0], j = tri[1], k = tri[2];
        if (i < 0 || i >= nV || j < 0 || j >= nV || k < 0 || k >= nV) continue;

        // Signed area = 0.5 * ((uj-ui)*(vk-vi) - (uk-ui)*(vj-vi))
        // Negative means the UV triangle has opposite winding (flipped).
        double signedArea2 = (u_(j) - u_(i)) * (v_(k) - v_(i))
                           - (u_(k) - u_(i)) * (v_(j) - v_(i));
        if (signedArea2 < 0.0) ++count;
    }
    return count;
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

    std::cerr << "[UVGParam] Written UV OBJ: " << filename << "\n";
    return true;
}
