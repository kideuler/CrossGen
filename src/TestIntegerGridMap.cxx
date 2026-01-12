// Utility to construct an integer grid map (MIQ-style seamless parameterization)
// on a CutMesh.
//
// This is intentionally similar in spirit to TestCutMesh.cxx:
//   1) Load mesh
//   2) Solve PolyVector field + detect singularities
//   3) Cut mesh into a topological disk (CutMesh)
//   4) Build integer grid map (IntegerGridMap)
//   5) Optionally write an OBJ with vt coordinates
//   6) Test for inverted elements
//   7) Validate quad mesh structure (integer isolines)

#include <iostream>
#include <string>
#include <vector>
#include <filesystem>

#include "polyvec/PolyVectors.hxx"
#include "polyvec/CutMesh.hxx"
#include "polyvec/IntegerGridMap.hxx"
#include "polyvec/QuadMeshValidation.hxx"

#include <cmath>

namespace fs = std::filesystem;

// Check if a value is close to an integer
inline bool isInteger(double val, double tol = 1e-6) {
    return std::abs(val - std::round(val)) < tol;
}

// Test a single mesh file, returns 0 on success (no flipped triangles), non-zero on failure
// If outFlipped is not null, stores the flipped triangle count
// If outValidation is not null, stores the validation result
// If outNonIntegerSingularities is not null, stores the count of non-integer singularities
int testMesh(const std::string &path, const std::string &outPath, bool verbose, 
             int *outFlipped = nullptr, QuadMeshValidation::ValidationResult *outValidation = nullptr,
             int *outNonIntegerSingularities = nullptr) {
    if (outFlipped) *outFlipped = -1;
    if (outNonIntegerSingularities) *outNonIntegerSingularities = -1;
    try {
        Mesh m(path);

        // Build the cross field and compute singularities (PolyField stores its own mesh copy).
        PolyField field(m);
        field.solveForPolyCoeffs();
        field.convertToFieldVectors();
        field.computeUSingularities();

        if (verbose) {
            std::cout << "Detected singularities: " << field.uSingularities.size() << "\n";
        }

        // Cut into a disk (MIQ-style). Cutting implementation is assumed correct.
        CutMesh cm(field);
        auto rep = cm.sanityCheck();
        if (verbose) {
            std::cout << "Cut mesh: looksLikeDisk=" << (rep.looksLikeDisk ? "yes" : "no")
                      << ", singularitiesOnBoundary=" << (rep.allSingularitiesOnBoundary ? "yes" : "no")
                      << ", boundaryComponents=" << rep.boundaryComponents
                      << "\n";
        }
        if (!rep.looksLikeDisk) {
            std::cerr << "[WARN] Cut mesh did not pass the disk sanity check.\n";
        }

        // Build integer grid map.
        IntegerGridMap igm(cm);
        IntegerGridMap::Options opt;
        // opt.h = 0.0 means auto-compute h based on mesh size
        opt.h = 0.0;
        opt.target_grid_lines = 50;      // Coarser grid for faster validation (~2 integer lines across mesh)
        opt.stiffening_iterations = 100;  // Many iterations for continuous phase
        opt.stiffening_c = 2.0;          // Stronger stiffening
        opt.stiffening_d = 20.0;         // Higher clamp
        opt.do_rounding = true;          // Enable rounding by default (with soft seams, should be safe)
        opt.max_round_iterations = 500;  // More rounding iterations

        if (!igm.solve(opt)) {
            std::cerr << "IntegerGridMap::solve failed.\n";
            return 2;
        }

        // Use default h (will use effective h from solve)
        const auto st = igm.computeStats();
        if (outFlipped) *outFlipped = st.flipped_triangles;
        if (verbose) {
            std::cout << "Integer grid map stats: flippedTriangles=" << st.flipped_triangles
                      << ", maxLambda=" << st.max_lambda
                      << ", meanLambda=" << st.mean_lambda
                      << "\n";
        }

        // Check that all singularities are at integer UV positions
        const auto &sings = cm.getSingularities();
        const auto &uvCoords = igm.uv();
        const auto &origToCut = cm.getOriginalToCutVertices();
        int nonIntegerSingularities = 0;
        
        if (verbose && !sings.empty()) {
            std::cout << "Singularity UV positions:\n";
        }
        
        for (size_t i = 0; i < sings.size(); ++i) {
            int origVid = sings[i].first;
            int index4 = sings[i].second;
            const char* typeStr = (index4 == 1) ? "+1" : (index4 == -1) ? "-1" : "??";
            
            if (origVid >= 0 && origVid < static_cast<int>(origToCut.size()) && !origToCut[origVid].empty()) {
                int cutVid = origToCut[origVid][0];
                if (cutVid >= 0 && cutVid < static_cast<int>(uvCoords.size())) {
                    const Point &uv = uvCoords[cutVid];
                    bool uIsInt = isInteger(uv[0]);
                    bool vIsInt = isInteger(uv[1]);
                    
                    if (verbose) {
                        std::cout << "  Singularity " << i << " (index=" << typeStr << "): "
                                  << "UV = (" << uv[0] << ", " << uv[1] << ")";
                        if (uIsInt && vIsInt) {
                            std::cout << " [INTEGER]";
                        } else {
                            std::cout << " [NON-INTEGER]";
                        }
                        std::cout << "\n";
                    }
                    
                    if (!uIsInt || !vIsInt) {
                        nonIntegerSingularities++;
                    }
                }
            }
        }
        
        if (outNonIntegerSingularities) *outNonIntegerSingularities = nonIntegerSingularities;
        
        if (verbose) {
            if (sings.empty()) {
                std::cout << "No singularities in mesh.\n";
            } else {
                std::cout << "Singularity integer check: " << (sings.size() - nonIntegerSingularities) 
                          << "/" << sings.size() << " at integer positions\n";
            }
        }

        // Run quad mesh validation
        QuadMeshValidation validator(igm);
        auto validation = validator.validate();
        if (outValidation) *outValidation = validation;
        
        if (verbose) {
            std::cout << "Quad mesh validation:\n";
            std::cout << "  valid=" << (validation.valid ? "yes" : "no") << "\n";
            std::cout << "  grid_vertices=" << validation.num_grid_vertices << "\n";
            std::cout << "  grid_edges=" << validation.num_grid_edges << "\n";
            std::cout << "  quad_faces=" << validation.num_quad_faces << "\n";
            std::cout << "  non_quad_faces=" << validation.num_non_quad_faces << "\n";
            std::cout << "  self_crossings=" << validation.num_self_crossings << "\n";
            std::cout << "  coverage=" << (validation.coverage_complete ? "complete" : "incomplete") << "\n";
            
            if (!validation.face_valence_histogram.empty()) {
                std::cout << "  face valence histogram: ";
                for (const auto &kv : validation.face_valence_histogram) {
                    std::cout << kv.first << "-gons:" << kv.second << " ";
                }
                std::cout << "\n";
            }
            
            if (!validation.messages.empty()) {
                std::cout << "  messages:\n";
                for (const auto &msg : validation.messages) {
                    std::cout << "    " << msg << "\n";
                }
            }
        }

        if (!outPath.empty()) {
            if (igm.writeOBJWithUV(outPath)) {
                std::cout << "Wrote cut mesh UV OBJ to: " << outPath << "\n";
            } else {
                std::cerr << "Failed to write UV OBJ: " << outPath << "\n";
                return 3;
            }
            
            // Write original (uncut) mesh with UV texture coordinates
            std::string origPath = outPath.substr(0, outPath.rfind('.')) + "_original.obj";
            if (igm.writeOriginalMeshWithUV(origPath)) {
                std::cout << "Wrote original mesh UV OBJ to: " << origPath << "\n";
            }
            
            // Also write the integer grid to an OBJ for visualization
            std::string gridPath = outPath.substr(0, outPath.rfind('.')) + "_grid.obj";
            validator.writeGridOBJ(gridPath);
            std::cout << "Wrote grid OBJ to: " << gridPath << "\n";
        }

        // Return status: 0=perfect, 4=flipped, 6=validation issues, 7=non-integer singularities
        if (st.flipped_triangles > 0) return 4;
        if (nonIntegerSingularities > 0) return 7;
        if (!validation.valid) return 6;
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 5;
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <mesh.obj> [out_uv.obj]\n";
        std::cerr << "       " << argv[0] << " --test-all <meshes_directory>\n";
        return 1;
    }

    const std::string arg1 = argv[1];
    
    // Test all meshes in a directory
    if (arg1 == "--test-all") {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " --test-all <meshes_directory>\n";
            return 1;
        }
        
        const std::string meshDir = argv[2];
        std::vector<std::string> meshFiles;
        
        // Collect all .obj files in the directory
        for (const auto &entry : fs::directory_iterator(meshDir)) {
            if (entry.path().extension() == ".obj") {
                meshFiles.push_back(entry.path().string());
            }
        }
        
        // Sort for consistent ordering
        std::sort(meshFiles.begin(), meshFiles.end());
        
        if (meshFiles.empty()) {
            std::cerr << "No .obj files found in " << meshDir << "\n";
            return 1;
        }
        
        std::cout << "Testing " << meshFiles.size() << " mesh files for inverted elements, singularity positions, and quad mesh validity...\n";
        std::cout << std::string(80, '=') << "\n";
        
        int totalPassed = 0;
        int totalFailed = 0;
        std::vector<std::pair<std::string, int>> failures;
        
        // Totals for validation
        int totalQuadFaces = 0;
        int totalNonQuadFaces = 0;
        int totalSelfCrossings = 0;
        int totalNonIntegerSingularities = 0;
        
        for (const auto &meshPath : meshFiles) {
            std::string filename = fs::path(meshPath).filename().string();
            std::cout << "Testing: " << filename << " ... ";
            std::cout.flush();
            
            int flippedCount = 0;
            int nonIntSingCount = 0;
            QuadMeshValidation::ValidationResult validation;
            int result = testMesh(meshPath, "", false, &flippedCount, &validation, &nonIntSingCount);
            
            if (result == 0) {
                std::cout << "\033[32mPASS\033[0m (quads:" << validation.num_quad_faces 
                          << " non-quads:" << validation.num_non_quad_faces
                          << " crossings:" << validation.num_self_crossings << ")\n";
                totalPassed++;
            } else if (result == 4) {
                std::cout << "\033[33mWARN\033[0m (" << flippedCount << " inverted elements)\n";
                totalFailed++;
                failures.emplace_back(filename, result);
            } else if (result == 7) {
                std::cout << "\033[33mWARN\033[0m (" << nonIntSingCount << " non-integer singularities)\n";
                totalFailed++;
                failures.emplace_back(filename, result);
            } else if (result == 6) {
                std::cout << "\033[33mWARN\033[0m (validation failed: ";
                if (validation.num_self_crossings > 0) std::cout << "crossings:" << validation.num_self_crossings << " ";
                if (validation.num_non_quad_faces > 0) std::cout << "non-quads:" << validation.num_non_quad_faces << " ";
                if (!validation.coverage_complete) std::cout << "incomplete coverage";
                std::cout << ")\n";
                totalFailed++;
                failures.emplace_back(filename, result);
            } else {
                std::cout << "\033[31mFAIL\033[0m (error code " << result << ")\n";
                totalFailed++;
                failures.emplace_back(filename, result);
            }
            
            totalQuadFaces += validation.num_quad_faces;
            totalNonQuadFaces += validation.num_non_quad_faces;
            totalSelfCrossings += validation.num_self_crossings;
            if (nonIntSingCount > 0) totalNonIntegerSingularities += nonIntSingCount;
        }
        
        std::cout << std::string(80, '=') << "\n";
        std::cout << "Summary: " << totalPassed << "/" << meshFiles.size() << " passed";
        if (totalFailed > 0) {
            std::cout << ", " << totalFailed << " failed/warned";
        }
        std::cout << "\n";
        std::cout << "Totals: " << totalQuadFaces << " quads, " << totalNonQuadFaces << " non-quads, " 
                  << totalSelfCrossings << " self-crossings, "
                  << totalNonIntegerSingularities << " non-integer singularities\n";
        
        if (!failures.empty()) {
            std::cout << "\nFiles with issues:\n";
            for (const auto &f : failures) {
                std::cout << "  - " << f.first;
                if (f.second == 4) {
                    std::cout << " (inverted elements)";
                } else if (f.second == 7) {
                    std::cout << " (non-integer singularities)";
                } else if (f.second == 6) {
                    std::cout << " (validation failed)";
                } else {
                    std::cout << " (error " << f.second << ")";
                }
                std::cout << "\n";
            }
        }
        
        return (totalFailed == 0) ? 0 : 1;
    }
    
    // Single mesh test (original behavior)
    const std::string path = argv[1];
    const std::string outPath = (argc >= 3) ? argv[2] : std::string();
    return testMesh(path, outPath, true);
}
