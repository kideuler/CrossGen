// Utility to load a mesh, solve a cross field with SIPG, cut it open along the
// Sec. 4.1 cuts, and run the UMBER frame field optimization (Sec. 4.2 of
// Wang et al. 2022) on top, reporting the cut report, the energy split and the
// internal singularities before and after.
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include "Parameterization/HarmonicCut.hxx"
#include "SIPG/SIPG.hxx"
#include "UMBER/UMBER.hxx"

static const int MBO_MAX_STEPS = 500;

static void printEnergy(const char *label, const UMBER::EnergyTerms &e) {
    std::cout << std::scientific << std::setprecision(4)
              << "  " << label
              << " total " << e.total
              << " = smooth " << e.smooth
              << " + w_a * align " << e.align
              << " + w_r * reg " << e.reg << "\n";
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <mesh.obj> [gamma] [max_steps] [lbfgs_iters]\n";
        return 1;
    }

    const std::string path = argv[1];
    const double gamma = (argc >= 3) ? std::stod(argv[2]) : 10.0;
    const int maxSteps = (argc >= 4) ? std::stoi(argv[3]) : MBO_MAX_STEPS;
    const int lbfgsIters = (argc >= 5) ? std::stoi(argv[4]) : 0; // 0 = UMBER default

    std::shared_ptr<Mesh> mesh;
    try {
        mesh = std::make_shared<Mesh>(path);
    } catch (const std::exception &e) {
        std::cout << "\033[31m[FAIL]\033[0m Failed to load mesh: " << e.what() << "\n";
        return 2;
    }

    // --- Initial cross field: p=0 DG/SIPG MBO --------------------------------
    SIPG sipg(mesh, maxSteps, gamma);
    sipg.initialize();
    sipg.runMBO();
    sipg.computeSingularities();

    std::cout << path << ": " << mesh->triangles.size() << " triangles, "
              << sipg.singularVertices.size() << " SIPG singularities\n";

    // --- Cuts, Sec. 4.1 ------------------------------------------------------
    HarmonicCut hc(mesh);
    const auto &cutReport = hc.getReport();
    std::cout << "  boundary loops " << cutReport.boundaryLoops
              << " (" << cutReport.voids << " void(s)), cuts made " << cutReport.cutsMade
              << ", cut edges " << hc.getCutEdges().size() << "\n";
    for (const auto &m : cutReport.messages) std::cout << "    " << m << "\n";
    if (cutReport.isDisk)
        std::cout << "\033[32m[PASS]\033[0m Cut mesh is a disk.\n";
    else
        std::cout << "\033[31m[FAIL]\033[0m Cut mesh is not a disk.\n";

    // --- Frame field optimization, Eq. (1) -----------------------------------
    try {
        UMBER umber(sipg, hc);
        if (lbfgsIters > 0) umber.setMaxIterations(lbfgsIters);
        umber.initialize();

        const UMBER::EnergyTerms before = umber.energy();
        const auto singBefore = umber.internalSingularities();
        std::cout << "  cut edges used as C: " << umber.getCutEdges().size() << "\n";
        printEnergy("initial", before);
        std::cout << "  initial internal singularities: " << singBefore.size() << "\n";

        umber.optimize();

        const UMBER::EnergyTerms after = umber.energy();
        const auto singAfter = umber.internalSingularities();
        printEnergy("final  ", after);
        std::cout << std::scientific << std::setprecision(3)
                  << "  L-BFGS iterations " << umber.iterations()
                  << ", final |grad| " << umber.gradientNorm() << "\n";
        const auto sharp = umber.sharpTurns();
        std::cout << "  final internal singularities: " << singAfter.size()
                  << ", shortest frame vector " << umber.minFrameNorm() << "\n";
        std::cout << "  interior edges turning > 45 deg: " << sharp.size();
        if (!sharp.empty()) std::cout << " (worst " << sharp.front().second << " deg)";
        std::cout << "\n";

        // --- Boundary corners, Sec. 4.2's intended destination ---------------
        const auto corners = umber.boundarySingularities();
        int convex = 0, reflex = 0, other = 0, quarterSum = 0;
        for (const auto &[vid, k] : corners) {
            quarterSum += k;
            if (k == 1) ++convex;
            else if (k == -1) ++reflex;
            else ++other;
        }
        for (const auto &[vid, w] : singAfter) quarterSum += 4 * static_cast<int>(w);

        const int chi = 1 - cutReport.voids;
        std::cout << std::defaultfloat
                  << "  boundary singularities: " << corners.size() << " ("
                  << convex << " convex, " << reflex << " reflex, " << other << " other), "
                  << "index sum " << quarterSum / 4.0 << ", Euler characteristic " << chi << "\n";
        if (quarterSum == 4 * chi)
            std::cout << "\033[32m[PASS]\033[0m Corner indices sum to the Euler characteristic.\n";
        else if (hc.getCutEdges().empty())
            std::cout << "\033[31m[FAIL]\033[0m Corner indices sum to " << quarterSum / 4.0
                      << ", expected " << chi << ".\n";
        else
            std::cout << "  (off by " << (quarterSum - 4 * chi) / 4.0
                      << "; the free transitions on the cuts carry the difference)\n";

        if (after.total <= before.total)
            std::cout << "\033[32m[PASS]\033[0m Energy did not increase.\n";
        else
            std::cout << "\033[31m[FAIL]\033[0m Energy increased.\n";

        if (singAfter.empty())
            std::cout << "\033[32m[PASS]\033[0m Field is free of internal singularities.\n";
        else
            std::cout << "\033[31m[FAIL]\033[0m " << singAfter.size()
                      << " internal singularity(ies) remain.\n";
    } catch (const std::exception &e) {
        std::cout << "\033[31m[FAIL]\033[0m UMBER failed: " << e.what() << "\n";
        return 3;
    }

    return 0;
}
