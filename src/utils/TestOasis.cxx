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
#include "crossfield/CrossField.hxx"
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

// Orientation control (Sec. 5.1, following Huang et al. 2008). The guiding
// field is a constant cross rotated off the axes of the square, which is the
// 2008 Fig. 5 experiment: nothing about the domain prefers that direction, so
// any alignment in the result comes from the orientation energy alone. 30
// degrees is deliberately not 45, where a cross field is indistinguishable from
// its own diagonal.
static const double GUIDE_ANGLE_DEG = 30.0;
static const double ORIENTATION_GAMMAS[] = {0.0, 0.01, 0.1, 1.0, 10.0, 100.0};

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

// A cross field pointing the same way everywhere, prescribed rather than
// solved: CrossField stores the 4-symmetry representation vector e^{4 i t} per
// vertex, so a constant field is one constant complex number.
static std::shared_ptr<CrossField> constantCrossField(const std::shared_ptr<Mesh> &mesh,
                                                      double degrees) {
    auto field = std::make_shared<CrossField>(mesh);
    const double t = degrees * M_PI / 180.0;
    field->u_k = Eigen::VectorXcd::Constant(static_cast<Eigen::Index>(mesh->vertices.size()),
                                            std::polar(1.0, 4.0 * t));
    field->u_k_prev = field->u_k;
    return field;
}

// Sweep gamma on a square carrying that field. The two numbers that matter pull
// against each other: the misalignment the orientation term is there to remove,
// and the Helmholtz residual that keeps the cells square and evenly spaced.
static bool runOrientationCase(const std::shared_ptr<Mesh> &mesh,
                               double lambda,
                               const std::string &prefix) {
    std::cout << "=== Orientation control, square with a constant cross field at "
              << GUIDE_ANGLE_DEG << " degrees ===\n";

    auto field = constantCrossField(mesh, GUIDE_ANGLE_DEG);

    // Keep the guiding field two quad cells clear of the boundary. The 2014
    // paper guides orientation only "in regions away from boundary curves or
    // features": inside that band the Eq. 11 conditions already dictate which
    // way the field runs, and a guiding direction that disagrees with them can
    // only be bought at the cost of the Helmholtz term, i.e. of cell shape.
    const double band = 2.0 * M_PI / std::sqrt(-lambda);
    const Eigen::VectorXd mask = OASIS(mesh, lambda).boundaryClearanceMask(band);
    std::cout << "Guided      : " << static_cast<int>(mask.sum()) << " of "
              << mesh->vertices.size() << " vertices (everything more than "
              << band << " from the boundary)\n";

    std::cout << "  gamma     misalign(deg)   |L_r f - lambda f|\n";
    bool ok = true;

    for (double gamma : ORIENTATION_GAMMAS) {
        OASIS oasis(field, lambda);
        oasis.setOrientationWeight(gamma);
        oasis.setOrientationMask(mask);
        try {
            oasis.solve();
        } catch (const std::exception &e) {
            std::cout << "\033[31m[FAIL]\033[0m gamma = " << gamma << ": " << e.what() << "\n\n";
            return false;
        }

        std::cout << "  " << std::setw(7) << gamma
                  << "   " << std::setw(11) << std::fixed << std::setprecision(2)
                  << oasis.orientationError()
                  << "   " << std::setw(16) << std::setprecision(3) << oasis.residual()
                  << std::defaultfloat << "\n";

        // Keep the two ends of the sweep to look at: no guidance, and the
        // default weight.
        if (gamma == 0.0 || gamma == oasis.getOrientationWeight()) {
            const std::string path = prefix + (gamma == 0.0 ? "_orient_off.vtk"
                                                            : "_orient_on.vtk");
            try {
                writeVTK(path, *mesh, oasis.f, "quasi_eigenfunction");
            } catch (const std::exception &e) {
                std::cout << "\033[31m[FAIL]\033[0m " << e.what() << "\n\n";
                ok = false;
            }
        }
    }

    // The term is only worth having if it actually turns the field: check the
    // default weight against the unguided solve rather than trusting the table
    // to be read.
    OASIS off(field, lambda), on(field, lambda);
    off.setOrientationWeight(0.0);
    off.setOrientationMask(mask);  // measure over the same region
    on.setOrientationMask(mask);
    try {
        off.solve();
        on.solve();  // default gamma
    } catch (const std::exception &e) {
        std::cout << "\033[31m[FAIL]\033[0m " << e.what() << "\n\n";
        return false;
    }
    const double errOff = off.orientationError();
    const double errOn = on.orientationError();
    if (!(errOn < 0.5 * errOff)) {
        std::cout << "\033[31m[FAIL]\033[0m the default weight barely moved the "
                     "misalignment (" << errOff << " -> " << errOn << " degrees)\n\n";
        return false;
    }

    std::cout << "  misalignment at the default gamma = " << on.getOrientationWeight()
              << ": " << errOff << " -> " << errOn << " degrees\n";

    // The vibration pass minimizes the same combined energy, so it should even
    // out the two local amplitudes without giving the orientation back.
    const double eaBefore = on.vibrationEnergy();
    try {
        on.enhanceVibration(VIBRATION_ITERATIONS);
    } catch (const std::exception &e) {
        std::cout << "\033[31m[FAIL]\033[0m vibration enhancement: " << e.what() << "\n\n";
        return false;
    }
    const double errVib = on.orientationError();
    std::cout << "  after vibration enhancement: mean E_a " << eaBefore << " -> "
              << on.vibrationEnergy() << ", misalignment " << errOn << " -> "
              << errVib << " degrees\n";
    if (errVib > 1.5 * errOn) {
        std::cout << "\033[31m[FAIL]\033[0m the vibration pass undid the orientation\n\n";
        return false;
    }

    if (!ok) return false;
    std::cout << "\033[32m[OK]\033[0m Wrote " << prefix << "_orient_off.vtk and "
              << prefix << "_orient_on.vtk\n\n";
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
    ok &= runOrientationCase(square, lambda, prefix);

    if (!ok) {
        std::cout << "\033[31mOne or more cases failed.\033[0m\n";
        return 3;
    }
    std::cout << "\033[32mAll cases passed.\033[0m\n";
    return 0;
}
