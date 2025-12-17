#include "PolyVectors.hxx"
#include <fstream>
#include <iomanip>

PolyField::PolyField(Mesh &mesh) : mesh(mesh) {
    // Initialize polynomial coefficients
    int nTriangles = static_cast<int>(mesh.triangles.size());
    int nBdyTriangles = static_cast<int>(mesh.boundaryTriangles.size());

    // reserve space in sparse matrix L
    // Assuming L will be an nTriangles x nTriangles operator (e.g., Laplacian-like)
    L.resize(nTriangles, nTriangles);
    L.reserve(3 * nTriangles);

    // fill in L
    for (int t = 0; t < nTriangles; ++t) {
        // initialize all polynomial coefficients to zero
        polyCoeffs.push_back(PolyCoeffs{ std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0) });

        // initialize matrix L
        L.insert(t, t) = 0.0;
        for (int e = 0; e < 3; ++e) {
            int neighbor = mesh.triangleAdjacency[t][e];
            if (neighbor != -1) {
                // interior edge: add connection to neighbor
                L.insert(t, neighbor) = -1.0;
                L.coeffRef(t, t) += 1.0;
            }
        }
    }

    // Loop through boundary triangles to set polynomial coefficients and dirichlet conditions
    for (int i = 0; i < nBdyTriangles; ++i) {
        int t = mesh.boundaryTriangles[i][0];
        int e = mesh.boundaryTriangles[i][1];

        // Get the two vertices of the boundary edge
        const Triangle &tri = mesh.triangles[t];
        int v0 = tri[e];
        int v1 = tri[(e + 1) % 3];  // next vertex in CCW order 
        Point p0 = mesh.vertices[v0];
        Point p1 = mesh.vertices[v1];

        // Compute tangent vector along the edge
        Point tangent = { p1[0] - p0[0], p1[1] - p0[1] };
        // Compute polynomial coefficients from tangent vector
        polyCoeffs[t] = computePolyCoeffsFromTangentVector(tangent);

        // Set Dirichlet condition in L (make row t a unit row)
        // Zero out all existing nonzeros in row t
        for (int col = 0; col < L.outerSize(); ++col) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(L, col); it; ++it) {
                if (it.row() == t) {
                    it.valueRef() = 0.0;
                }
            }
        }
        // Remove zeros from sparsity pattern
        L.prune(0.0);
        // Set the diagonal to 1
        L.coeffRef(t, t) = 1.0;
        
    }

    // factorize L after all modifications
    solver.compute(L);
    if(solver.info() != Eigen::Success) {
        return;
    }
}

void PolyField::solveForPolyCoeffs() {
    // Prepare rhs vector b
    int nTriangles = static_cast<int>(mesh.triangles.size());
    b_re = Eigen::VectorXd::Zero(nTriangles);
    b_im = Eigen::VectorXd::Zero(nTriangles);

    for (int m = 0; m < 2; ++m) {

        // Set rhs for boundary triangles based on polynomial coefficients
        for (int i = 0; i < static_cast<int>(mesh.boundaryTriangles.size()); ++i) {
            int t = mesh.boundaryTriangles[i][m];
            // For simplicity, set rhs to real part of x_0 coefficient
            b_re(t) = polyCoeffs[t][m].real();
            b_im(t) = polyCoeffs[t][m].imag();
        }

        // Solve the linear system L * x = b
        x_re = solver.solve(b_re);
        x_im = solver.solve(b_im);

        // Update polynomial coefficients based on solution x
        for (int t = 0; t < nTriangles; ++t) {
            // Here we simply set x_0 to the solved value and keep x_1 unchanged
            polyCoeffs[t][m] = std::complex<double>(x_re(t), x_im(t));
        }
    }
}

void PolyField::convertToFieldVectors() {
    int nTriangles = static_cast<int>(mesh.triangles.size());
    field.resize(nTriangles);

    for (int t = 0; t < nTriangles; ++t) {
        // Extract polynomial coefficients
        PolyCoeffs coeffs = polyCoeffs[t];

       Eigen::MatrixXcd M(4,4);
       M << std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), -coeffs[0],
           std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
           std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0), -coeffs[1],
           std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0);

       // Compute eigenvalues (no eigenvectors needed)
       Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(M, /* computeEigenvectors = */ false);

       // Put the 4 complex roots into a std::array (double precision)
       Eigen::Matrix<std::complex<double>, 4, 1> ev = es.eigenvalues();
        std::array<std::complex<double>, 4> roots;
        for (int k = 0; k < 4; ++k) {
            roots[k] = ev(k);
        }

        // Compute field vectors from roots
       auto dirs = extract_two_directions(roots);
       field[t].u = dirs.first;
       field[t].v = dirs.second;
    }
}

bool PolyField::writeVTK(const std::string &filename) const {
    // Legacy VTK UnstructuredGrid format, ASCII
    std::ofstream out(filename);
    if (!out) return false;

    const auto &V = mesh.vertices;
    const auto &T = mesh.triangles;
    int nPts = static_cast<int>(V.size());
    int nTri = static_cast<int>(T.size());

    out << "# vtk DataFile Version 3.0\n";
    out << "PolyField export\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    // Points (promote 2D points to 3D with z=0)
    out << "POINTS " << nPts << " float\n";
    out.setf(std::ios::fixed); out << std::setprecision(7);
    for (const auto &p : V) {
        out << static_cast<float>(p[0]) << " "
            << static_cast<float>(p[1]) << " 0\n";
    }

    // Cells: each triangle is written as: 3 i j k
    // Total size = sum over cells of (1 + numIndices) = nTri * (1 + 3)
    out << "\nCELLS " << nTri << " " << nTri * 4 << "\n";
    for (const auto &tri : T) {
        out << 3 << " " << tri[0] << " " << tri[1] << " " << tri[2] << "\n";
    }

    // Cell types: VTK_TRIANGLE == 5
    out << "\nCELL_TYPES " << nTri << "\n";
    for (int i = 0; i < nTri; ++i) out << 5 << "\n";

    // Cell data: per-triangle vectors
    out << "\nCELL_DATA " << nTri << "\n";
    out << "VECTORS field_u float\n";
    for (int i = 0; i < nTri; ++i) {
        float ux = 0.0f, uy = 0.0f;
        if (i < static_cast<int>(field.size())) {
            ux = static_cast<float>(field[i].u[0]);
            uy = static_cast<float>(field[i].u[1]);
        }
        out << ux << " " << uy << " 0\n";
    }
    out << "\nVECTORS field_v float\n";
    for (int i = 0; i < nTri; ++i) {
        float vx = 0.0f, vy = 0.0f;
        if (i < static_cast<int>(field.size())) {
            vx = static_cast<float>(field[i].v[0]);
            vy = static_cast<float>(field[i].v[1]);
        }
        out << vx << " " << vy << " 0\n";
    }

    return true;
}