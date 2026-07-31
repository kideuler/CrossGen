// Meshes a few test domains with TestHelper/TriangleMesher2D, computes the
// OASIS quasi-eigenfunction on each, and writes mesh + field to legacy VTK.
//
// The QE is the scalar field whose Morse-Smale complex becomes the quad mesh in
// Ling et al. 2014 (refs/Huang2014.pdf). Its boundary conditions df/dn = 0 and
// d^2f/dn^2 = 1 make every boundary curve an extended minimal integral curve,
// so the field should show a grid of extrema running square to the boundary
// with no direction field anywhere in the pipeline. The curved domains are the
// interesting cases: the square could be faked with a uniform cross field, the
// half disk and the ellipse could not.
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "OASIS/OASIS.hxx"
#include "TestHelper.hxx"
#include "mesh/Mesh.hxx"

// The quad edge length implied by lambda is pi / sqrt(-lambda), since the QE
// behaves like cos(u * sqrt(-lambda)) along each axis (Eqs. 5 and 7). The
// default targets quads of about 0.25.
static const double DEFAULT_LAMBDA = -158.0;

// Background triangle size. It must be well below the quad size for the
// critical points of the QE to be resolved; the paper uses a factor of 0.25.
static const double DEFAULT_H = 0.05;

// Domain dimensions. All three are on a comparable scale so that one lambda
// gives a sensible element count everywhere.
static const double BOX_SIZE = 5.0;
static const double DISK_RADIUS = 2.5;
static const double ELLIPSE_A = 3.0;
static const double ELLIPSE_B = 2.0;
// The Sec. 3.4 pathological case: rotationally symmetric, so the unenhanced QE
// vibrates radially and shows rings rather than isolated spots.
static const double DISK_RADIUS_FULL = 1.5;

// Gauss-Newton iterations for the vibration pass.
static const int VIBRATION_ITERATIONS = 10;

// Write an unstructured grid of triangles with one scalar per point.
static void writeVTK(const std::string &filename,
                     const Mesh &mesh,
                     const Eigen::VectorXd &field,
                     const std::string &fieldName) {
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Failed to open file for writing: " + filename);
    }

    const size_t numVertices = mesh.vertices.size();
    const size_t numTriangles = mesh.triangles.size();

    out << "# vtk DataFile Version 3.0\n";
    out << "OASIS quasi-eigenfunction\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    out << std::setprecision(16);

    out << "POINTS " << numVertices << " double\n";
    for (const Point &p : mesh.vertices) {
        out << p[0] << " " << p[1] << " 0\n";
    }

    out << "CELLS " << numTriangles << " " << 4 * numTriangles << "\n";
    for (const Triangle &tri : mesh.triangles) {
        out << "3 " << tri[0] << " " << tri[1] << " " << tri[2] << "\n";
    }

    out << "CELL_TYPES " << numTriangles << "\n";
    for (size_t t = 0; t < numTriangles; ++t) {
        out << "5\n";  // VTK_TRIANGLE
    }

    out << "POINT_DATA " << numVertices << "\n";
    out << "SCALARS " << fieldName << " double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (size_t v = 0; v < numVertices; ++v) {
        out << field[static_cast<Eigen::Index>(v)] << "\n";
    }

    out.close();
    if (!out) {
        throw std::runtime_error("Error writing to file: " + filename);
    }
}

// Solve on one domain and write the result. Returns false on failure.
static bool runCase(const std::string &label,
                    const std::shared_ptr<Mesh> &mesh,
                    double lambda,
                    const std::string &outPath,
                    const std::string &vibPath) {
    std::cout << "=== " << label << " ===\n";
    std::cout << "Triangles   : " << mesh->triangles.size() << "\n";
    std::cout << "Vertices    : " << mesh->vertices.size()
              << " (" << mesh->boundaryVertices.size() << " on the boundary)\n";

    OASIS oasis(mesh, lambda);
    try {
        oasis.assemble();
        std::cout << "KKT system  : " << oasis.getKKT().rows() << " x "
                  << oasis.getKKT().cols() << ", " << oasis.getKKT().nonZeros()
                  << " nonzeros, " << oasis.numConstraints() << " constraints\n";
        oasis.solve();
    } catch (const std::exception &e) {
        std::cout << "\033[31m[FAIL]\033[0m OASIS failed: " << e.what() << "\n\n";
        return false;
    }

    const Eigen::VectorXd &f = oasis.f;

    // Sanity: the alignment conditions should hold to solver accuracy, and the
    // Helmholtz residual should be nonzero (this is a *quasi*-eigenfunction).
    const double bcError =
        (oasis.getConstraintMatrix() * f - oasis.getConstraintRhs()).cwiseAbs().maxCoeff();

    std::cout << "QE range    : [" << f.minCoeff() << ", " << f.maxCoeff() << "]\n";
    std::cout << "|B f - C|inf: " << bcError << "\n";
    std::cout << "Helmholtz   : |L_r f - lambda f| = " << oasis.residual() << "\n";

    if (!f.allFinite()) {
        std::cout << "\033[31m[FAIL]\033[0m Solution contains non-finite values\n\n";
        return false;
    }
    if (bcError > 1e-6 * std::max(1.0, f.cwiseAbs().maxCoeff())) {
        std::cout << "\033[31m[FAIL]\033[0m Boundary conditions not satisfied\n\n";
        return false;
    }

    try {
        writeVTK(outPath, *mesh, f, "quasi_eigenfunction");
    } catch (const std::exception &e) {
        std::cout << "\033[31m[FAIL]\033[0m " << e.what() << "\n\n";
        return false;
    }

    // Vibration enhancement (Sec. 3.4), written alongside so the two can be
    // compared directly. On the disk this is the difference between rings and
    // isolated spots.
    const double eaBefore = oasis.vibrationEnergy();
    try {
        oasis.enhanceVibration(VIBRATION_ITERATIONS);
    } catch (const std::exception &e) {
        std::cout << "\033[31m[FAIL]\033[0m vibration enhancement: " << e.what() << "\n\n";
        return false;
    }
    const double eaAfter = oasis.vibrationEnergy();
    const double bcAfter =
        (oasis.getConstraintMatrix() * oasis.f - oasis.getConstraintRhs()).cwiseAbs().maxCoeff();

    std::cout << "mean E_a    : " << eaBefore << " -> " << eaAfter
              << "   (0 = isotropic vibration, 2 = one direction only)\n";
    std::cout << "|B f - C|   : " << bcAfter << "  after enhancement (penalty, not exact)\n";

    try {
        writeVTK(vibPath, *mesh, oasis.f, "quasi_eigenfunction");
    } catch (const std::exception &e) {
        std::cout << "\033[31m[FAIL]\033[0m " << e.what() << "\n\n";
        return false;
    }

    std::cout << "\033[32m[OK]\033[0m Wrote " << outPath << " and " << vibPath << "\n\n";
    return true;
}

int main(int argc, char **argv) {
    if (argc >= 2 && std::string(argv[1]) == "-h") {
        std::cerr << "Usage: " << argv[0] << " [lambda] [h] [output_prefix]\n";
        std::cerr << "  lambda : Helmholtz parameter, negative (default "
                  << DEFAULT_LAMBDA << ")\n";
        std::cerr << "  h      : background triangle edge length (default "
                  << DEFAULT_H << ")\n";
        std::cerr << "  prefix : VTK filename prefix (default oasis)\n";
        return 1;
    }

    const double lambda = (argc >= 2) ? std::stod(argv[1]) : DEFAULT_LAMBDA;
    const double h = (argc >= 3) ? std::stod(argv[2]) : DEFAULT_H;
    const std::string prefix = (argc >= 4) ? argv[3] : "oasis";

    const double quadSize = M_PI / std::sqrt(-lambda);
    std::cout << "lambda      : " << lambda << "  -> quad edge ~ " << quadSize << "\n";
    std::cout << "Triangle h  : " << h << "  (h/quad = " << h / quadSize
              << ", paper wants <= 0.25)\n\n";

    // ------------------------------------------------------------------
    // Build the domains
    // ------------------------------------------------------------------
    std::shared_ptr<Mesh> square, halfDisk, ellipse, disk;
    try {
        square = TestHelper::createBox(0.0, 0.0, BOX_SIZE, BOX_SIZE, h);
        // A 180 degree sweep with equal semi-axes closes the arc with a
        // straight diameter, giving the upper half disk.
        halfDisk = TestHelper::createEllipse(0.0, 0.0, DISK_RADIUS, DISK_RADIUS, 180.0, h);
        ellipse = TestHelper::createEllipse(0.0, 0.0, ELLIPSE_A, ELLIPSE_B, 360.0, h);
        disk = TestHelper::createCircle(0.0, 0.0, DISK_RADIUS_FULL, h);
    } catch (const std::exception &e) {
        std::cout << "\033[31m[FAIL]\033[0m Meshing failed: " << e.what() << "\n";
        return 2;
    }

    // ------------------------------------------------------------------
    // Assemble, solve and write each
    // ------------------------------------------------------------------
    bool ok = true;
    ok &= runCase("Square " + std::to_string(static_cast<int>(BOX_SIZE)) + " x " +
                      std::to_string(static_cast<int>(BOX_SIZE)),
                  square, lambda, prefix + "_square.vtk", prefix + "_square_vib.vtk");
    ok &= runCase("Upper half disk, radius " + std::to_string(DISK_RADIUS),
                  halfDisk, lambda, prefix + "_halfdisk.vtk", prefix + "_halfdisk_vib.vtk");
    ok &= runCase("Ellipse, semi-axes " + std::to_string(ELLIPSE_A) + " and " +
                      std::to_string(ELLIPSE_B),
                  ellipse, lambda, prefix + "_ellipse.vtk", prefix + "_ellipse_vib.vtk");
    ok &= runCase("Full disk, radius " + std::to_string(DISK_RADIUS_FULL) +
                      " (the Sec. 3.4 rings-vs-spots case)",
                  disk, lambda, prefix + "_disk.vtk", prefix + "_disk_vib.vtk");

    if (!ok) {
        std::cout << "\033[31mOne or more cases failed.\033[0m\n";
        return 3;
    }
    std::cout << "\033[32mAll cases passed.\033[0m\n";
    return 0;
}
