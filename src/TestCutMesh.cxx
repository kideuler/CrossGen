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
        std::cerr << "Usage: " << argv[0] << " <mesh.obj> [cut_mesh_out.obj]\n";
        return 1;
    }

    const std::string path = argv[1];
    const std::string outPath = (argc >= 3) ? argv[2] : std::string();

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

        if (!outPath.empty()) {
            if (cm.writeOBJ(outPath)) {
                std::cout << "Wrote cut mesh to: " << outPath << "\n";
            } else {
                std::cerr << "Failed to write cut mesh OBJ: " << outPath << "\n";
            }
        }

        return rep.looksLikeDisk ? 0 : 2;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 3;
    }
}
