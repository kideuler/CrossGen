// Utility to load a mesh, run the p=0 DG/SIPG MBO cross-field solver,
// and print singularities and per-triangle field values.
#include <iostream>
#include <iomanip>
#include <memory>
#include <string>

#include "sipg/SIPG.hxx"
#include "IGM/CutMesh.hxx"

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
        std::cerr << "Failed to load mesh: " << e.what() << "\n";
        return 2;
    }

    std::cerr << "Loaded mesh: " << mesh->triangles.size() << " triangles, "
              << mesh->vertices.size() << " vertices, "
              << mesh->boundaryEdges.size() << " boundary edges\n";

    // ------------------------------------------------------------------
    // Initialize SIPG solver (assembles M, K, b and factorises A = M + tau*K)
    // ------------------------------------------------------------------
    SIPG sipg(mesh, maxSteps, gamma);
    sipg.initialize();

    std::cerr << "SIPG initialized (gamma=" << gamma
              << ", tau automatically set)\n";

    // ------------------------------------------------------------------
    // MBO iteration loop
    // ------------------------------------------------------------------
    int  stepCount = 0;
    bool converged = false;
    double ntris   = static_cast<double>(mesh->triangles.size());

    std::cerr << "Running MBO iterations (max " << maxSteps << ")...\n";
    for (int i = 0; i < maxSteps; ++i) {
        sipg.step();
        ++stepCount;

        if (i > 0 && i % 50 == 0) {
            std::cerr << "  step " << stepCount
                      << "  error = " << std::scientific << std::setprecision(4)
                      << sipg.error << "\n";
        }

        if (sipg.error < 2.0 * ntris * 1e-5) {
            converged = true;
            break;
        }
    }

    if (converged) {
        std::cerr << "MBO converged at step " << stepCount
                  << "  (error = " << sipg.error << ")\n";
    } else {
        std::cerr << "MBO reached max steps (" << stepCount
                  << ")  (error = " << sipg.error << ")\n";
    }

    // ------------------------------------------------------------------
    // Detect singularities
    // ------------------------------------------------------------------
    sipg.computeSingularities();

    std::cerr << "Found " << sipg.singularVertices.size() << " singularity/singularities\n";
    for (const auto &[vertIdx, index] : sipg.singularVertices) {
        std::cout << "singularity  vertex=" << vertIdx
                  << "  index=" << std::fixed << std::setprecision(4) << index << "\n";
    }

    // ------------------------------------------------------------------
    // Cut mesh from SIPG field
    // ------------------------------------------------------------------
    std::cerr << "Building cut mesh from SIPG field...\n";
    CutMesh cm(sipg);
    const auto &rep = cm.sanityCheck();

    std::cout << "\n# cut mesh sanity report\n";
    std::cout << "  Triangle components : " << rep.triangleComponents
              << " (connected=" << (rep.trianglesConnected ? "yes" : "no") << ")\n";
    std::cout << "  Boundary components : " << rep.boundaryComponents << "\n";
    std::cout << "  Euler characteristic: " << rep.eulerCharacteristic << "\n";
    std::cout << "  Singularities on boundary: " << (rep.allSingularitiesOnBoundary ? "yes" : "no") << "\n";
    std::cout << "  Looks like disk     : " << (rep.looksLikeDisk ? "yes" : "no") << "\n";
    if (!rep.messages.empty()) {
        for (const auto &msg : rep.messages)
            std::cout << "  - " << msg << "\n";
    }
    if (rep.looksLikeDisk && rep.allSingularitiesOnBoundary)
        std::cout << "\033[32m[PASS]\033[0m\n";
    else
        std::cout << "\033[31m[FAIL]\033[0m\n";
    
    return 0;
}
