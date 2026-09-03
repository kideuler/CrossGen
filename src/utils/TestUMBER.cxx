// Utility to load a mesh, solve a cross field with DualMBO, cut it open along the
// Sec. 4.1 cuts, and run the UMBER frame field optimization (Sec. 4.2 of
// Wang et al. 2022) on top, reporting the cut report, the energy split and the
// internal singularities before and after.
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include "Parameterization/HarmonicCut.hxx"
#include "dualmbo/DualMBO.hxx"
#include "UMBER/BlockLayout.hxx"
#include "UMBER/ChordCollapse.hxx"
#include "UMBER/MotorcycleGraph.hxx"
#include "UMBER/Polysquare.hxx"
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
        std::cerr << "Usage: " << argv[0]
                  << " <mesh.obj> [gamma] [max_steps] [lbfgs_iters] [chord_max_width]"
                  << " [arrival_tol]\n";
        return 1;
    }

    const std::string path = argv[1];
    const double gamma = (argc >= 3) ? std::stod(argv[2]) : 10.0;
    const int maxSteps = (argc >= 4) ? std::stoi(argv[3]) : MBO_MAX_STEPS;
    const int lbfgsIters = (argc >= 5) ? std::stoi(argv[4]) : 0; // 0 = UMBER default
    // ChordCollapse::Settings::maxWidth, the heuristic that says how thin a
    // chord has to be to be worth collapsing.
    const double chordMaxWidth =
        (argc >= 6) ? std::stod(argv[5]) : ChordCollapse::Settings().maxWidth;
    // MotorcycleGraph::setArrivalTolerance, in mean mesh edges. 0 turns off the
    // rescue for a ray that has lost its way; see the note on that setter.
    const double arrivalTol = (argc >= 7) ? std::stod(argv[6]) : 0.20;

    std::shared_ptr<Mesh> mesh;
    try {
        mesh = std::make_shared<Mesh>(path);
    } catch (const std::exception &e) {
        std::cout << "\033[31m[FAIL]\033[0m Failed to load mesh: " << e.what() << "\n";
        return 2;
    }

    // --- Initial cross field: p=0 dual-mesh MBO --------------------------------
    DualMBO dualMBO(mesh, maxSteps, gamma);
    dualMBO.initialize();
    dualMBO.runMBO();
    dualMBO.computeSingularities();

    std::cout << path << ": " << mesh->triangles.size() << " triangles, "
              << dualMBO.singularVertices.size() << " DualMBO singularities\n";

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
        UMBER umber(dualMBO, hc);
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

        // --- Polysquare parameterization, Sec. 4.3 --------------------------
        Polysquare poly(umber, hc);
        poly.solve();
        const Polysquare::Report &pr = poly.getReport();

        std::cout << std::defaultfloat << "  polysquare: L-BFGS " << pr.iterations
                  << " iterations, transitions";
        for (int k : poly.getTransitions()) std::cout << " " << k * 90 << "deg";
        if (poly.getTransitions().empty()) std::cout << " none";
        std::cout << " (Eq. 7 residual " << std::fixed << std::setprecision(1)
                  << pr.transitionDeg << " deg)\n";
        std::cout << "    boundary turns: " << pr.turns << " (field asked for "
                  << pr.expectedTurns << "), shortest straight run " << pr.shortestRun
                  << " edges\n";
        std::cout << std::fixed << std::setprecision(2)
                  << "    segments off axis before snapping: " << pr.suspectSegments
                  << " past 10 deg, worst " << pr.worstSegmentDeg << " deg\n";
        std::cout << "    boundary runs put on a shared iso-line: " << pr.runsAligned
                  << ", worst move " << pr.worstRunAlign << " boundary edge(s)\n";
        std::cout << "    corners snapped onto a shared iso-line: " << pr.cornersSnapped
                  << ", worst move " << pr.worstCornerSnap << " boundary edge(s)";
        if (pr.cornerSnapConflicts > 0)
            std::cout << ", " << pr.cornerSnapConflicts << " cluster(s) left alone";
        std::cout << "\n";
        std::cout << std::fixed << std::setprecision(3)
                  << "    boundary alignment: mean " << pr.meanAlignDeg << " deg, worst "
                  << pr.maxAlignDeg << " deg; length ratio " << pr.lengthRatio << "\n"
                  << "    scaled Jacobian: min " << pr.minScaledJacobian << ", avg "
                  << pr.avgScaledJacobian << "; flipped triangles " << pr.flips << "\n"
                  << "    E_arap " << pr.arap << "  E_l1 " << pr.l1 << "  E_cor ";
        if (pr.cor > 0.0) std::cout << std::scientific << std::setprecision(2) << pr.cor << "\n";
        else std::cout << "(off)\n";

        if (pr.flips == 0)
            std::cout << "\033[32m[PASS]\033[0m Parameterization has no flipped triangles.\n";
        else
            std::cout << "\033[31m[FAIL]\033[0m " << pr.flips
                      << " flipped triangle(s) in the parameterization.\n";
        if (pr.meanAlignDeg < 5.0)
            std::cout << "\033[32m[PASS]\033[0m Boundary is axis aligned.\n";
        else
            std::cout << "\033[31m[FAIL]\033[0m Boundary is " << pr.meanAlignDeg
                      << " deg off axis on average.\n";
        // On a holed model the count also picks up the seams, where the cut
        // banks meet the boundary, so it is only exact when there are no cuts.
        if (pr.turns == pr.expectedTurns)
            std::cout << "\033[32m[PASS]\033[0m Boundary turns exactly where the field has a corner.\n";
        else if (hc.getCutEdges().empty())
            std::cout << "\033[31m[FAIL]\033[0m Boundary turns " << pr.turns << " times, but the "
                      << "field asked for " << pr.expectedTurns << ".\n";
        else
            std::cout << "  (" << pr.turns - pr.expectedTurns
                      << " turns beyond the field's corners, at the seams)\n";

        // Strip the directory and the extension for the output names.
        std::string stem = path;
        const size_t slash = stem.find_last_of("/\\");
        if (slash != std::string::npos) stem = stem.substr(slash + 1);
        const size_t dot = stem.find_last_of('.');
        if (dot != std::string::npos) stem = stem.substr(0, dot);

        // --- Meta-block structure, Sec. 5 (motorcycle graph only) -----------
        MotorcycleGraph mg(poly);
        mg.setArrivalTolerance(arrivalTol);
        mg.build();
        const MotorcycleGraph::Report &mr = mg.getReport();
        std::cout << std::defaultfloat
                  << "  blocks: " << mr.blocks << " from " << mr.motorcycles << " motorcycle(s) ("
                  << mr.reachedBoundary << " reached the boundary, "
                  << mr.crossings << " crossing(s)";
        if (mr.ranOut > 0) std::cout << ", " << mr.ranOut << " ran out";
        if (mr.skippedCorners > 0)
            std::cout << ", " << mr.skippedCorners << " direction(s) the polysquare had no room for";
        if (mr.arrivals > 0)
            std::cout << ", " << mr.arrivals << " stopped on a corner they had lost their way past";
        if (mr.duplicates > 0)
            std::cout << ", " << mr.duplicates << " dropped as the same line drawn backwards";
        std::cout << "), " << mr.nodes << " node(s), smallest block " << mr.smallestBlock << " triangles";
        if (mr.mergedSlivers > 0)
            std::cout << ", " << mr.mergedSlivers << " unresolved block(s) folded in";
        std::cout << "\n";

        if (mr.ranOut == 0)
            std::cout << "\033[32m[PASS]\033[0m Every iso-line reached the boundary.\n";
        else if (pr.flips > 0 || mr.skippedCorners > 0)
            std::cout << "  (" << mr.ranOut << " ray(s) stopped early, where the parameterization "
                      << "is already inconsistent)\n";
        else
            std::cout << "\033[31m[FAIL]\033[0m " << mr.ranOut
                      << " iso-line(s) did not reach the boundary.\n";

        // --- The blocks as a graph, and chord collapse -----------------------
        BlockLayout bl(mg);
        bl.build();
        const BlockLayout::Report &br = bl.getReport();
        const QuadLayout::Report lr = bl.getLayout().getReport();
        std::cout << "  block structure: " << lr.faces << " block(s), " << lr.arcs
                  << " side(s), " << lr.nodes << " node(s) on " << br.boundaryLoops
                  << " boundary loop(s)";
        if (br.danglingEnds > 0) std::cout << ", " << br.danglingEnds << " free end(s)";
        if (br.degenerateArcs > 0) std::cout << ", " << br.degenerateArcs << " side(s) dropped";
        if (br.unplacedNodes > 0) std::cout << ", " << br.unplacedNodes << " node(s) unplaced";
        std::cout << "\n";

        // A block that did not come out four-sided is a place the tracing left
        // open, so it is only news where the tracing was clean: a ray the
        // parameterization had no room for, or one that never got out, leaves a
        // block that no number of sides can be right for.
        const bool tracingClean = (pr.flips == 0 && mr.skippedCorners == 0 && mr.ranOut == 0);
        if (lr.badFaces == 0 && lr.arcCrossings == 0)
            std::cout << "\033[32m[PASS]\033[0m Every block came out four-sided ("
                      << lr.faces << " against " << mr.blocks << " from the flood).\n";
        else if (!tracingClean)
            std::cout << "  (" << lr.badFaces << " block(s) not four-sided and "
                      << lr.arcCrossings << " side(s) crossing, where the tracing is already "
                      << "incomplete)\n";
        else
            std::cout << "\033[31m[FAIL]\033[0m " << lr.badFaces << " block(s) not four-sided, "
                      << lr.arcCrossings << " side(s) crossing.\n";

        ChordCollapse::Settings ccs;
        ccs.maxWidth = chordMaxWidth;
        ChordCollapse cc(bl.getLayout(), ccs);
        cc.run();
        const ChordCollapse::Report &cr = cc.getReport();
        std::cout << std::fixed << std::setprecision(2)
                  << "  chord collapse: " << cr.blocksBefore << " -> " << cr.blocksAfter
                  << " block(s) in " << cr.collapses << " collapse(s)";
        if (cr.widestCollapsed > 0.0)
            std::cout << ", widest " << cr.widestCollapsed << " of a mean side";
        if (cr.rolledBack > 0)
            std::cout << ", " << cr.rolledBack << " undone (" << cr.rbBlocks << " block count, "
                      << cr.rbBad << " not four-sided, " << cr.rbCrossings << " crossing, "
                      << cr.rbArea << " area)";
        std::cout << "\n    " << cr.chordsSeen << " chord(s) left, " << cr.collapsible
                  << " still collapsible";
        for (int i = 1; i < static_cast<int>(ChordCollapse::Block::Count); ++i) {
            if (cr.blockCount[i] == 0) continue;
            std::cout << ", " << cr.blockCount[i] << " "
                      << ChordCollapse::blockName(static_cast<ChordCollapse::Block>(i));
        }
        std::cout << "\n";

        // What the operation is answerable for is that it left the structure no
        // worse than it found it, whatever state the tracing handed over.
        const QuadLayout::Report &ar = cc.getLayout().getReport();
        if (ar.badFaces <= lr.badFaces && ar.arcCrossings <= lr.arcCrossings &&
            ar.danglingEnds <= lr.danglingEnds)
            std::cout << "\033[32m[PASS]\033[0m Collapsing left the structure no worse: "
                      << ar.badFaces << " block(s) not four-sided, " << ar.arcCrossings
                      << " side(s) crossing.\n";
        else
            std::cout << "\033[31m[FAIL]\033[0m collapsing made the structure worse: "
                      << lr.badFaces << " -> " << ar.badFaces << " block(s) not four-sided, "
                      << lr.arcCrossings << " -> " << ar.arcCrossings << " side(s) crossing.\n";

        const std::string uvFile = stem + "_polysquare.vtu";
        const std::string srcFile = stem + "_source.vtu";
        const std::string blkFile = stem + "_blocks.vtu";
        const std::string trcFile = stem + "_motorcycles.vtu";
        const std::string layFile = stem + "_blocklayout.vtu";
        const std::string simFile = stem + "_simplified.vtu";
        if (poly.writeVTU(uvFile) && poly.writeSourceVTU(srcFile) &&
            mg.writeBlocksVTU(blkFile) && mg.writeTracesVTU(trcFile) &&
            bl.getLayout().writeFacesVTU(layFile) && cc.getLayout().writeFacesVTU(simFile))
            std::cout << "  wrote " << uvFile << " (parameter domain), " << srcFile
                      << " (input domain, uv as a point field),\n         " << blkFile
                      << " (blocks as cell data), " << trcFile << " (the traces),\n         "
                      << layFile << " (the blocks as polygons) and " << simFile
                      << " (after chord collapse)\n";
        else
            std::cout << "\033[31m[FAIL]\033[0m could not write the VTK output.\n";
    } catch (const std::exception &e) {
        std::cout << "\033[31m[FAIL]\033[0m UMBER failed: " << e.what() << "\n";
        return 3;
    }

    return 0;
}
