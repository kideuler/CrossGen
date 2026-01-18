// Utility to cut a mesh into a topological disk (MIQ-style) and run sanity checks.

#include <iostream>
#include <string>

#include "polyvec/PolyVectors.hxx"
#include "polyvec/CutMesh.hxx"

static void printReport(const CutMesh::SanityReport &rep) {
    std::cout << "Sanity report\n";
    std::cout << "  Triangle components: " << rep.triangleComponents
              << " (connected=" << (rep.trianglesConnected ? "yes" : "no") << ")\n";
    std::cout << "  Boundary components: " << rep.boundaryComponents << "\n";
    std::cout << "  Euler characteristic: " << rep.eulerCharacteristic << "\n";
    std::cout << "  Singularities on boundary: " << (rep.allSingularitiesOnBoundary ? "yes" : "no") << "\n";
    std::cout << "  Looks like disk: " << (rep.looksLikeDisk ? "yes" : "no") << "\n";
    if (!rep.messages.empty()) {
        std::cout << "Messages:\n";
        for (const auto &m : rep.messages) {
            std::cout << "  - " << m << "\n";
        }
    }

    // print overall green pass or red fail
    if (rep.looksLikeDisk && rep.allSingularitiesOnBoundary) {
        std::cout << "\033[32m[PASS]\033[0m\n";
    } else {
        std::cout << "\033[31m[FAIL]\033[0m\n";
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <mesh.obj> [--cut <output.obj>]\n";
        std::cerr << "Options:\n";
        std::cerr << "  --cut <file.obj>     Write cut mesh to file\n";
        return 1;
    }

    const std::string path = argv[1];
    std::string cutOutPath;

    // Parse arguments
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--cut") {
            if (i + 1 < argc) {
                cutOutPath = argv[++i];
            }
        }
    }

    try {
        Mesh m(path);

        // Build the cross field and compute singularities (PolyField stores its own mesh copy).
        PolyField field(m);
        field.solveForPolyCoeffs();
        field.convertToFieldVectors();
        field.computeUSingularities();

        std::cout << "Detected singularities: " << field.uSingularities.size() << "\n";

        CutMesh cm(field);
        auto rep = cm.sanityCheck();
        printReport(rep);

        if (!cutOutPath.empty()) {
            if (cm.writeOBJ(cutOutPath)) {
                std::cout << "Wrote cut mesh to: " << cutOutPath << "\n";
            } else {
                std::cerr << "Failed to write cut mesh OBJ: " << cutOutPath << "\n";
            }
        }

        return rep.looksLikeDisk ? 0 : 2;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 3;
    }
}
