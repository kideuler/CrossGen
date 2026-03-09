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

    // 4. Print summary
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
    }

    return 0;
}
