#ifndef __POLYVECTORS_HXX__
#define __POLYVECTORS_HXX__

#include <complex>
#include <cmath>

// eigen includes
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Mesh.hxx"

typedef std::array<std::complex<double>, 2> PolyCoeffs; // x_0 and x_1 in Polynomial P(z) = z^4 + x_1*z^2 + x_0

struct Field {
    Point u;
    Point v;
};

std::pair<Point, Point> extract_two_directions(const std::array<std::complex<double>, 4>& roots) {
    // Convert roots to 2D unit vectors (Points)
    std::array<Point, 4> vecs;
    for (int k = 0; k < 4; ++k) {
        double rx = roots[k].real();
        double ry = roots[k].imag();
        double len = std::sqrt(rx * rx + ry * ry);
        Point v{0.0, 0.0};
        if (len > 1e-12) {
            v[0] = rx / len;
            v[1] = ry / len;
        }
        vecs[k] = v;
    }

    auto sumSqNorm = [](const Point &a, const Point &b) {
        double sx = a[0] + b[0];
        double sy = a[1] + b[1];
        return sx * sx + sy * sy;
    };

    // Find two roots that are roughly opposite (minimize ||v + w||)
    int i0 = 0, j0 = 1;
    double best = sumSqNorm(vecs[0], vecs[1]);
    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) {
            double score = sumSqNorm(vecs[i], vecs[j]);
            if (score < best) { best = score; i0 = i; j0 = j; }
        }
    }

    // First direction from the best opposite pair
    Point d1 = vecs[i0];

    // Second direction: one of the remaining two (the other should be ~ -d2)
    int k1 = -1, k2 = -1;
    for (int k = 0; k < 4; ++k) {
        if (k != i0 && k != j0) {
            if (k1 < 0) k1 = k; else k2 = k;
        }
    }
    Point d2 = (k1 >= 0) ? vecs[k1] : Point{0.0, 0.0};

    return {d1, d2};
}


PolyCoeffs computePolyCoeffsFromTangentVector(const Point &tangent) {
    // normalize tangent
    double len = std::sqrt(tangent[0]*tangent[0] + tangent[1]*tangent[1]);
    std::complex<double> tx(tangent[0]/len, tangent[1]/len);

    // compute the normal vector (90 degree rotation CCW)
    std::complex<double> tn(-tx.imag(), tx.real());

    // compute coefficients
    PolyCoeffs coeffs;
    coeffs[0] = tx*tx*tn*tn; // x_0
    coeffs[1] = -(tx*tx + tn*tn); // x_1
    return coeffs;
}

class PolyField {
 public:
    std::vector<Field> field; // per-triangle field vectors
    PolyField() = default;

    PolyField(Mesh &mesh); // Initialize field on the given mesh

    void solveForPolyCoeffs(); // Solve for the polynomial coefficients

    PolyCoeffs getPolyCoeffsForTriangle(int triangleIndex) const {
        return polyCoeffs[triangleIndex];
    }

    void convertToFieldVectors(); // Convert polynomial coefficients to field vectors

    // Write mesh and per-triangle field vectors to a VTK legacy unstructured grid file
    // - filename: path to .vtk file
    // Returns true on success, false on failure
    bool writeVTK(const std::string &filename) const;

private:
    Mesh mesh;
    std::vector<PolyCoeffs> polyCoeffs; // per-triangle polynomial coefficients

    // rhs vectors for the linear system
    Eigen::VectorXd b_re;
    Eigen::VectorXd b_im;
    Eigen::VectorXd x_re;
    Eigen::VectorXd x_im;

    // Sparse matrix for the linear system
    Eigen::SparseMatrix<double> L;

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

};

#endif // __POLYVECTORS_HXX__