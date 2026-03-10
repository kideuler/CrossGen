// Utility to load a mesh, compute the MBO cross field, and trace separatrices.
#include <iostream>
#include <iomanip>
#include <memory>
#include <string>

#include "crossfield/CrossField.hxx"
#include "tracing/SeparatrixTrace.hxx"

static const int MBO_MAX_STEPS = 500;

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <mesh.obj>\n";
        return 1;
    }

    std::string path = argv[1];
    std::shared_ptr<Mesh> mesh;
    try {
        mesh = std::make_shared<Mesh>(path);
    } catch (const std::exception &e) {
        std::cerr << "Failed to load mesh: " << e.what() << "\n";
        return 2;
    }

    std::cerr << "Loaded mesh: " << mesh->triangles.size() << " triangles, "
              << mesh->vertices.size() << " vertices\n";

    // 1. Initialize and run MBO cross field
    auto crossField = std::make_shared<CrossField>(mesh);
    crossField->initialize(1);

    double nv = static_cast<double>(crossField->u_k.size());
    int stepCount = 0;
    bool converged = false;

    std::cerr << "Running MBO iterations (max " << MBO_MAX_STEPS << ")...\n";
    for (int i = 0; i < MBO_MAX_STEPS; ++i) {
        crossField->step();
        stepCount++;
        if (crossField->error < 2.0 * nv * 1e-7) {
            converged = true;
            break;
        }
    }
    crossField->computeSingularities();

    if (converged) {
        std::cerr << "MBO converged at step " << stepCount
                  << " with error " << crossField->error << "\n";
    } else {
        std::cerr << "MBO reached max steps (" << stepCount
                  << ") with error " << crossField->error << "\n";
    }
    std::cerr << "Found " << crossField->singularTriangles.size() << " singularities\n";

    // 2. Initialize separatrix tracing
    auto separatrixTrace = std::make_shared<SeparatrixTrace>(crossField, true);
    std::cerr << "Initialized " << separatrixTrace->separatrices.size()
              << " separatrices from " << separatrixTrace->singularities.size()
              << " singularities\n";

    // 3. Trace until all separatrices are finished
    int traceStep = 0;
    while (!separatrixTrace->finishedTracing) {
        separatrixTrace->stepAndCheck();
        traceStep++;
        if (traceStep % 1000 == 0) {
            int active = 0;
            for (const auto &sep : separatrixTrace->separatrices) {
                if (sep.active) active++;
            }
            std::cerr << "  Trace step " << traceStep << ": "
                      << active << " active separatrices\n";
        }
    }

    std::cerr << "Tracing complete after " << traceStep << " steps\n";

    // 4. Print summary and detect angle jumps
    for (const auto &sep : separatrixTrace->separatrices) {
        const char *reason = "?";
        switch (sep.termination_reason) {
            case TerminationReason::RUNNING: reason = "RUNNING"; break;
            case TerminationReason::EXIT_BOUNDARY: reason = "EXIT_BOUNDARY"; break;
            case TerminationReason::CONNECT_TANGENTIAL_PRIMARY: reason = "CONNECT_TANGENTIAL_PRIMARY"; break;
            case TerminationReason::CONNECT_TANGENTIAL_SECONDARY: reason = "CONNECT_TANGENTIAL_SECONDARY"; break;
            case TerminationReason::LIMIT_CYCLE: reason = "LIMIT_CYCLE"; break;
            case TerminationReason::ORTHOGONAL_TO_SINGULARITY_SEPARATRIX: reason = "ORTHOGONAL_TO_SINGULARITY_SEPARATRIX"; break;
            case TerminationReason::MAX_STEPS_REACHED: reason = "MAX_STEPS_REACHED"; break;
            case TerminationReason::UNDEFINED: reason = "UNDEFINED"; break;
        }
        std::cerr << "  Separatrix " << sep.id
                  << ": " << sep.path.size() << " points, "
                  << reason << "\n";

        // Detect large direction jumps between consecutive segments
        for (size_t j = 2; j < sep.path.size(); ++j) {
            Point v_prev = sep.path[j-1].global_pos - sep.path[j-2].global_pos;
            Point v_curr = sep.path[j].global_pos - sep.path[j-1].global_pos;
            double dir_prev = std::atan2(v_prev[1], v_prev[0]);
            double dir_curr = std::atan2(v_curr[1], v_curr[0]);
            double jump = std::abs(wrap_pi(dir_curr - dir_prev));
            bool prevSingular = false;
            for (const auto &sing : separatrixTrace->singularities) {
                if (sing.triangleIndex == sep.path[j-1].face_id) { prevSingular = true; break; }
            }
            if (jump > M_PI / 3.0) {
                std::cerr << "    ** JUMP at point " << j << "/" << sep.path.size()
                          << ": angle change = " << (jump * 180.0 / M_PI) << " deg"
                          << "\n      pt[" << j-2 << "] pos=(" << sep.path[j-2].global_pos[0] << "," << sep.path[j-2].global_pos[1] << ")"
                          << " face=" << sep.path[j-2].face_id
                          << " trace_dir=" << std::fixed << std::setprecision(4) << sep.path[j-2].trace_direction
                          << " field_angle=" << sep.path[j-2].field_angle
                          << " edge=" << sep.path[j-2].local_edge_index
                          << " t=" << sep.path[j-2].edge_crossing_t
                          << "\n      pt[" << j-1 << "] pos=(" << sep.path[j-1].global_pos[0] << "," << sep.path[j-1].global_pos[1] << ")"
                          << " face=" << sep.path[j-1].face_id
                          << (prevSingular ? " (SINGULAR)" : "")
                          << " trace_dir=" << sep.path[j-1].trace_direction
                          << " field_angle=" << sep.path[j-1].field_angle
                          << " edge=" << sep.path[j-1].local_edge_index
                          << " t=" << sep.path[j-1].edge_crossing_t
                          << "\n      pt[" << j << "] pos=(" << sep.path[j].global_pos[0] << "," << sep.path[j].global_pos[1] << ")"
                          << " face=" << sep.path[j].face_id
                          << " trace_dir=" << sep.path[j].trace_direction
                          << " field_angle=" << sep.path[j].field_angle
                          << " edge=" << sep.path[j].local_edge_index
                          << " t=" << sep.path[j].edge_crossing_t
                          << "\n";
            }
        }
    }

    return 0;
}
