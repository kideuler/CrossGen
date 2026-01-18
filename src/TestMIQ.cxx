// Test UV parametrization quality from MIQ solver.
// Checks for flipped triangles (negative area in UV space).

#include <iostream>
#include <string>
#include <cmath>

#include "polyvec/PolyVectors.hxx"
#include "polyvec/CutMesh.hxx"
#include "polyvec/MIQ.hxx"

// Compute signed area of a triangle in UV space using cross product
// Positive = counter-clockwise, Negative = clockwise (flipped)
double signedTriangleArea(const Eigen::Vector2d& p0, 
                          const Eigen::Vector2d& p1, 
                          const Eigen::Vector2d& p2) {
    // Area = 0.5 * ((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0))
    return 0.5 * ((p1.x() - p0.x()) * (p2.y() - p0.y()) - 
                  (p2.x() - p0.x()) * (p1.y() - p0.y()));
}

// Check UV parametrization for flipped triangles
// A triangle is "flipped" if its UV winding is opposite to the original mesh winding
// Returns number of flipped triangles
int checkForFlippedTriangles(const MIQSolver& miq, const CutMesh& cm) {
    const auto& UV = miq.getUV();
    const auto& FUV = miq.getFUV();
    
    // Get original mesh to check winding
    const Mesh& origMesh = cm.getOriginalMesh();
    
    int numFlipped = 0;
    
    for (int f = 0; f < FUV.rows(); ++f) {
        // UV triangle
        Eigen::Vector2d uv0 = UV.row(FUV(f, 0));
        Eigen::Vector2d uv1 = UV.row(FUV(f, 1));
        Eigen::Vector2d uv2 = UV.row(FUV(f, 2));
        double uvArea = signedTriangleArea(uv0, uv1, uv2);
        
        // Original mesh triangle (using cut mesh face indices which match original)
        const auto& tri = origMesh.triangles[f];
        Eigen::Vector2d p0(origMesh.vertices[tri[0]][0], origMesh.vertices[tri[0]][1]);
        Eigen::Vector2d p1(origMesh.vertices[tri[1]][0], origMesh.vertices[tri[1]][1]);
        Eigen::Vector2d p2(origMesh.vertices[tri[2]][0], origMesh.vertices[tri[2]][1]);
        double origArea = signedTriangleArea(p0, p1, p2);
        
        // Flipped if signs differ (or if UV area is zero/degenerate)
        if ((uvArea > 0) != (origArea > 0)) {
            numFlipped++;
        }
    }
    
    return numFlipped;
}

static void printCutMeshReport(const CutMesh::SanityReport &rep) {
    std::cout << "Cut Mesh Sanity Report\n";
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
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <mesh.obj> [options]\n";
        std::cerr << "Options:\n";
        std::cerr << "  --uv <file.obj>      Write UV mesh to file\n";
        std::cerr << "  --qex <basename>     Output files for QEx: <basename>.obj (textured mesh)\n";
        std::cerr << "                       and <basename>.vval (vertex valences)\n";
        std::cerr << "  --quad <file.obj>    Extract and write coarse quad mesh\n";
        std::cerr << "  --resolution N       Set quad mesh resolution (default: 200, higher = finer)\n";
        std::cerr << "  --stiffness S        Set stiffness parameter (default: 5.0)\n";
        std::cerr << "  --iterations N       Set stiffening iterations (default: 50)\n";
        return 1;
    }

    const std::string path = argv[1];
    std::string uvOutPath;
    std::string quadOutPath;
    std::string qexBasename;
    double resolution = 100.0;
    double stiffness = 5.0;
    int iterations = 10;

    // Parse arguments
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--uv") {
            if (i + 1 < argc) {
                uvOutPath = argv[++i];
            }
        } else if (arg == "--quad") {
            if (i + 1 < argc) {
                quadOutPath = argv[++i];
            }
        } else if (arg == "--qex") {
            if (i + 1 < argc) {
                qexBasename = argv[++i];
            }
        } else if (arg == "--resolution") {
            if (i + 1 < argc) {
                resolution = std::stod(argv[++i]);
            }
        } else if (arg == "--stiffness") {
            if (i + 1 < argc) {
                stiffness = std::stod(argv[++i]);
            }
        } else if (arg == "--iterations") {
            if (i + 1 < argc) {
                iterations = std::stoi(argv[++i]);
            }
        }
    }

    try {
        std::cout << "Loading mesh: " << path << "\n";
        Mesh m(path);

        // Build the cross field and compute singularities
        std::cout << "\n--- Cross Field Computation ---\n";
        PolyField field(m);
        field.solveForPolyCoeffs();
        field.convertToFieldVectors();
        field.computeUSingularities();

        std::cout << "Detected singularities: " << field.uSingularities.size() << "\n";

        // Create cut mesh
        std::cout << "\n--- Cut Mesh ---\n";
        CutMesh cm(field);
        auto cutRep = cm.sanityCheck();
        printCutMeshReport(cutRep);

        if (!cutRep.looksLikeDisk) {
            std::cerr << "\033[31m[ERROR]\033[0m Cut mesh is not a valid topological disk. Cannot proceed with MIQ.\n";
            return 2;
        }

        // Run MIQ parametrization
        std::cout << "\n--- MIQ Parametrization ---\n";
        std::cout << "Resolution: " << resolution << "\n";
        std::cout << "Stiffness: " << stiffness << "\n";
        std::cout << "Iterations: " << iterations << "\n";
        
        MIQSolver miq(cm);
        miq.solve(resolution, stiffness, false, iterations, 5000, true, true);
        
        // Check for flipped triangles
        int numFlipped = checkForFlippedTriangles(miq, cm);
        int totalTriangles = miq.getFUV().rows();
        
        std::cout << "\n=== UV Sanity Check ===\n";
        std::cout << "Total triangles: " << totalTriangles << "\n";
        std::cout << "Flipped triangles: " << numFlipped << "\n";
        
        if (numFlipped == 0) {
            std::cout << "\033[32m[PASS]\033[0m No flipped triangles\n";
        } else {
            std::cout << "\033[31m[FAIL]\033[0m " << numFlipped << " flipped triangles ("
                      << (100.0 * numFlipped / totalTriangles) << "%)\n";
        }
        
        // Write UV mesh if requested
        if (!uvOutPath.empty()) {
            if (miq.writeUVMesh(uvOutPath)) {
                std::cout << "\nWrote UV mesh to: " << uvOutPath << "\n";
            } else {
                std::cerr << "Failed to write UV mesh: " << uvOutPath << "\n";
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
                    std::cerr << "Failed to write quad mesh: " << quadOutPath << "\n";
                }
            } else {
                std::cerr << "Failed to extract quad mesh\n";
            }
        }

        // Return based on UV quality
        return (numFlipped == 0) ? 0 : 1;
        
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 3;
    }
}
