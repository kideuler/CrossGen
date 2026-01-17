// Utility to cut a mesh into a topological disk (MIQ-style) and run sanity checks.
// Optionally computes MIQ parametrization.

#include <iostream>
#include <string>

#include "polyvec/PolyVectors.hxx"
#include "polyvec/CutMesh.hxx"
#include "polyvec/MIQ.hxx"

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
        std::cerr << "Usage: " << argv[0] << " <mesh.obj> [options]\n";
        std::cerr << "Options:\n";
        std::cerr << "  --miq [file.obj]     Compute MIQ parametrization, optionally write UV mesh\n";
        std::cerr << "  --cut <file.obj>     Write cut mesh to file\n";
        std::cerr << "  --quad <file.obj>    Extract and write coarse quad mesh (requires --miq)\n";
        std::cerr << "  --qex <basename>     Output files for QEx: <basename>.obj (textured mesh)\n";
        std::cerr << "                       and <basename>.vval (vertex valences). Requires --miq.\n";
        std::cerr << "  --resolution N       Set quad mesh resolution (default: 200, higher = finer)\n";
        return 1;
    }

    const std::string path = argv[1];
    std::string cutOutPath;
    bool runMIQ = false;
    std::string uvOutPath;
    std::string quadOutPath;
    std::string qexBasename;
    double resolution = 200.0;  // Default resolution (gradientSize)

    // Parse arguments
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--miq") {
            runMIQ = true;
            // Check if next argument is a filename (not another flag)
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                uvOutPath = argv[++i];
            }
        } else if (arg == "--cut") {
            if (i + 1 < argc) {
                cutOutPath = argv[++i];
            }
        } else if (arg == "--quad") {
            if (i + 1 < argc) {
                quadOutPath = argv[++i];
            }
        } else if (arg == "--qex") {
            if (i + 1 < argc) {
                qexBasename = argv[++i];
                runMIQ = true;  // QEx output requires MIQ
            }
        } else if (arg == "--resolution") {
            if (i + 1 < argc) {
                resolution = std::stod(argv[++i]);
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

        // Run MIQ parametrization if requested and mesh is valid
        if (runMIQ && rep.looksLikeDisk) {
            std::cout << "\n--- MIQ Parametrization ---\n";
            std::cout << "Resolution (gradientSize): " << resolution << "\n";
            MIQSolver miq(cm);
            
            // Parameters: gradientSize, stiffness, directRound, iter, localIter, doRound, singularityRound
            miq.solve(resolution, 5.0, false, 50, 5000, true, true);
            
            int flips = miq.numFlips();
            std::cout << "Final flips: " << flips << "\n";
            
            const auto& UV = miq.getUV();
            std::cout << "UV coordinates: " << UV.rows() << " vertices\n";
            
            if (!uvOutPath.empty()) {
                if (miq.writeUVMesh(uvOutPath)) {
                    std::cout << "Wrote UV mesh to: " << uvOutPath << "\n";
                } else {
                    std::cerr << "Failed to write UV mesh OBJ: " << uvOutPath << "\n";
                }
            }
            
            // Output files for QEx if requested
            if (!qexBasename.empty()) {
                std::cout << "\n--- QEx Output ---\n";
                
                std::string texturedPath = qexBasename + ".obj";
                std::string vvalPath = qexBasename + ".vval";
                
                if (miq.writeTexturedMesh(texturedPath)) {
                    std::cout << "Wrote textured mesh to: " << texturedPath << "\n";
                } else {
                    std::cerr << "Failed to write textured mesh: " << texturedPath << "\n";
                }
                
                if (miq.writeVertexValences(vvalPath)) {
                    std::cout << "Wrote vertex valences to: " << vvalPath << "\n";
                } else {
                    std::cerr << "Failed to write vertex valences: " << vvalPath << "\n";
                }
                
                std::cout << "To run QEx: cmdline_tool " << texturedPath << " output_quads.obj " << vvalPath << "\n";
            }
            
            // Extract and write quad mesh if requested
            if (!quadOutPath.empty()) {
                std::cout << "\n--- Quad Mesh Extraction ---\n";
                if (miq.extractQuadMesh()) {
                    const auto& qv = miq.getQuadVertices();
                    const auto& qf = miq.getQuadFaces();
                    std::cout << "Quad mesh: " << qv.rows() << " vertices, " << qf.rows() << " quads\n";
                    
                    if (miq.writeQuadMesh(quadOutPath)) {
                        std::cout << "Wrote quad mesh to: " << quadOutPath << "\n";
                    } else {
                        std::cerr << "Failed to write quad mesh OBJ: " << quadOutPath << "\n";
                    }
                } else {
                    std::cerr << "Failed to extract quad mesh\n";
                }
            }
        } else if (runMIQ && !rep.looksLikeDisk) {
            std::cout << "Skipping MIQ: cut mesh is not a valid topological disk.\n";
        }

        return rep.looksLikeDisk ? 0 : 2;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 3;
    }
}
