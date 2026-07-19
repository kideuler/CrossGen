// Utility to load a mesh, run the p=0 DG/SIPG MBO cross-field solver,
// and print singularities and per-triangle field values.
#include <iostream>
#include <iomanip>
#include <memory>
#include <string>

#include "sipg/SIPG.hxx"
#include "Parameterization/CutMesh.hxx"
#include "Parameterization/UVGParam.hxx"

static const int MBO_MAX_STEPS = 500;

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <mesh.obj> [gamma] [max_steps]\n";
        std::cerr << "  gamma     : SIPG penalty parameter (default 20.0)\n";
        std::cerr << "  max_steps : max MBO iterations   (default " << MBO_MAX_STEPS << ")\n";
        return 1;
    }

    std::string path = argv[1];
    double gamma    = (argc >= 3) ? std::stod(argv[2]) : 10.0;
    int    maxSteps = (argc >= 4) ? std::stoi(argv[3]) : MBO_MAX_STEPS;

    // ------------------------------------------------------------------
    // Load mesh
    // ------------------------------------------------------------------
    std::shared_ptr<Mesh> mesh;
    try {
        mesh = std::make_shared<Mesh>(path);
    } catch (const std::exception &e) {
        std::cout << "\033[31m[FAIL]\033[0m Failed to load mesh: " << e.what() << "\n";
        return 2;
    }

    // ------------------------------------------------------------------
    // Initialize SIPG solver (assembles M, K, b and factorises A = M + tau*K)
    // ------------------------------------------------------------------
    SIPG sipg(mesh, maxSteps, gamma);
    sipg.initialize();

    // ------------------------------------------------------------------
    // MBO iteration loop
    // ------------------------------------------------------------------
    int  stepCount = 0;
    bool converged = false;
    double ntris   = static_cast<double>(mesh->triangles.size());

    for (int i = 0; i < maxSteps; ++i) {
        sipg.step();
        ++stepCount;
        if (sipg.error < 2.0 * ntris * 1e-5) {
            converged = true;
            break;
        }
    }

    // ------------------------------------------------------------------
    // Detect singularities
    // ------------------------------------------------------------------
    sipg.computeSingularities();

    // ------------------------------------------------------------------
    // Cut mesh from SIPG field
    // ------------------------------------------------------------------
    CutMesh cm(sipg);
    const auto &rep = cm.sanityCheck();

    if (rep.looksLikeDisk && rep.allSingularitiesOnBoundary)
        std::cout << "\033[32m[PASS]\033[0m Cut mesh is a disk with singularities on boundary.\n";
    else
        std::cout << "\033[31m[FAIL]\033[0m Cut mesh sanity check failed.\n";

    // ------------------------------------------------------------------
    // Global UV parameterization
    // ------------------------------------------------------------------
    try {
        UVGParam uvp(cm);

        int flips = uvp.numFlippedTriangles();
        if (flips == 0)
            std::cout << "\033[32m[PASS]\033[0m No flipped triangles in UV space.\n";
        else
            std::cout << "\033[31m[FAIL]\033[0m " << flips << " flipped triangle(s) in UV space.\n";
    } catch (const std::exception &e) {
        std::cout << "\033[31m[FAIL]\033[0m UVGParam failed: " << e.what() << "\n";
    }

    return 0;
}