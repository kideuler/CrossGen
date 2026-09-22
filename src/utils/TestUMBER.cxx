// Utility to load a mesh, solve a cross field with DualMBO, cut it open along the
// Sec. 4.1 cuts, and run the UMBER frame field optimization (Sec. 4.2 of
// Wang et al. 2022) on top, reporting the cut report, the energy split and the
// internal singularities before and after.
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include "MERIDIAN/Interfaces.hxx"
#include "Parameterization/HarmonicCut.hxx"
#include "dualmbo/DualMBO.hxx"
#include "mesh/BlockQuadMesh.hxx"
#include "mesh/QuadMesh.hxx"
#include "mesh/TMOP.hxx"
#include "UMBER/BlockLayout.hxx"
#include "UMBER/MotorcycleGraph.hxx"
#include "UMBER/Polysquare.hxx"
#include "UMBER/UMBER.hxx"

static const int MBO_MAX_STEPS = 500;
// The pipelines' own default (MERIDIAN::Options::quadTargetEdge). The corpus
// is normalised into [0,1]^2, so this is one twentieth of the model, and it
// means the same thing to all three methods on purpose.
static const double QUAD_TARGET_EDGE = 0.05;

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
                  << " <mesh.obj> [gamma] [max_steps] [lbfgs_iters]"
                  << " [arrival_tol] [target_edge] [interface_weight]\n";
        return 1;
    }

    const std::string path = argv[1];
    const double gamma = (argc >= 3) ? std::stod(argv[2]) : 10.0;
    const int maxSteps = (argc >= 4) ? std::stoi(argv[3]) : MBO_MAX_STEPS;
    const int lbfgsIters = (argc >= 5) ? std::stoi(argv[4]) : 0; // 0 = UMBER default
    // MotorcycleGraph::setArrivalTolerance, in mean mesh edges. 0 turns off the
    // rescue for a ray that has lost its way; see the note on that setter.
    const double arrivalTol = (argc >= 6) ? std::stod(argv[5]) : 0.20;
    // Target mesh edge length, in the units of the model; the number the three
    // methods' mesh dialogs all take.
    const double targetEdge = (argc >= 7) ? std::stod(argv[6]) : QUAD_TARGET_EDGE;
    // Polysquare::setInterfaceWeight: how hard the material interfaces are
    // pulled onto an axis, relative to dS. Negative keeps the default.
    const double ifaceWeight = (argc >= 8) ? std::stod(argv[7]) : -1.0;

    // Checks that must hold whatever the model is. They are the structural
    // invariants of the block decomposition and the mesh on it -- every side
    // resolves, every interface separates two materials, the mesh is
    // conforming -- and nothing that depends on the shape being one the
    // polysquare can flatten. The flips, the internal singularities and the
    // boundary turn count are all printed and none of them is here, because a
    // model whose layout the method cannot produce is a fact about the model.
    int structuralFailures = 0;

    std::shared_ptr<Mesh> mesh;
    try {
        mesh = std::make_shared<Mesh>(path);
    } catch (const std::exception &e) {
        std::cout << "\033[31m[FAIL]\033[0m Failed to load mesh: " << e.what() << "\n";
        return 2;
    }

    // --- Stage 0b: the material interface network ----------------------------
    //
    // An interface is a feature in exactly the sense dS is -- a curve the
    // output has to keep, since an element straddling it carries two materials
    // -- so it is read before anything is built on the mesh and handed to
    // every stage that has to follow it. On a single-material model this finds
    // nothing and costs one pass over the edges.
    std::unique_ptr<Interfaces> interfaces;
    try {
        interfaces = std::make_unique<Interfaces>(mesh);
    } catch (const std::exception &e) {
        std::cout << "  (interface extraction failed: " << e.what()
                  << "; continuing as a single-material model)\n";
    }
    const bool multiMaterial = interfaces && interfaces->multiMaterial();
    if (multiMaterial) {
        const Interfaces::Report &ir = interfaces->getReport();
        std::cout << "  interfaces: " << ir.materials << " material(s), " << ir.interfaceEdges
                  << " interface edge(s) in " << ir.branches << " branch(es), " << ir.nodes
                  << " node(s)\n";
    }

    // --- Initial cross field: p=0 dual-mesh MBO --------------------------------
    DualMBO dualMBO(mesh, maxSteps, gamma);
    if (multiMaterial) dualMBO.setAlignedInteriorEdges(interfaces->interfaceEdges());
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
        if (multiMaterial) umber.setFeatureEdges(interfaces->interfaceEdges());
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

        // The same reading along the interfaces: where the frame turns against
        // a curve the layout has to follow is a corner of the structure just
        // as a corner of dS is.
        if (multiMaterial) {
            const auto icorners = umber.interfaceCorners();
            int conv = 0, refl = 0, high = 0;
            for (const auto &fc : icorners) {
                if (fc.quarters == 1) ++conv;
                else if (fc.quarters == -1) ++refl;
                else ++high;
            }
            std::cout << "  interface corners: " << icorners.size() << " (" << conv
                      << " convex, " << refl << " reflex";
            if (high > 0) std::cout << ", " << high << " higher order";
            std::cout << "), counted on both sides of every branch\n";
        }

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
        if (ifaceWeight >= 0.0) poly.setInterfaceWeight(ifaceWeight);
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
                  << pr.maxAlignDeg << " deg; length ratio " << pr.lengthRatio << "\n";
        if (pr.interfaceEdges > 0) {
            std::cout << "    interface alignment: mean " << pr.meanInterfaceAlignDeg
                      << " deg, worst " << pr.maxInterfaceAlignDeg << " deg over "
                      << pr.interfaceEdges << " edge(s)\n";
        }
        std::cout << std::fixed << std::setprecision(3)
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
        if (mr.continuations > 0)
            std::cout << ", " << mr.continuations << " carrying a line through an interface";
        if (mr.continuationsDropped > 0)
            std::cout << ", " << mr.continuationsDropped << " refused at the generation cap";
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

        // --- The blocks as a graph ------------------------------------------
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

        // --- The structure as the shared block decomposition -----------------
        //
        // The same object ATLAS, MERIDIAN and TORSION hand their blocks back
        // as, so what is checked here is what TestATLAS Case 14 and
        // TestMERIDIAN's Stage 8 check of it check: that the counts agree with
        // the stage's own report, that every side of every block resolves to a
        // macro edge with a polyline on it, and that an interface really does
        // separate two materials.
        BlockLayout::DecompositionReport dr;
        const BlockDecomposition D =
            BlockLayout::blockDecompositionOf(bl.getLayout(), *mesh, &dr);
        std::cout << std::defaultfloat
                  << "  decomposition: " << dr.blocks << " block(s) of " << dr.faces
                  << " face(s), " << D.edges.size() << " macro edge(s), " << D.vertices.size()
                  << " macrovertex/-ices, " << dr.materials << " material(s)";
        if (dr.notFourSided > 0) std::cout << ", " << dr.notFourSided << " not four-sided";
        if (dr.multiArcSides > 0) std::cout << ", " << dr.multiArcSides << " with a split side";
        if (dr.unmatchedSides > 0) std::cout << ", " << dr.unmatchedSides << " unmatched side(s)";
        if (dr.straddlingBlocks > 0)
            std::cout << ", " << dr.straddlingBlocks << " straddling an interface";
        std::cout << "\n";
        // What the counts above add up to. A face that did not qualify is a
        // piece of the model with no block on it, so this is the fraction of
        // the model the mesh below will actually cover.
        const double coverage = BlockLayout::coverageOf(dr);
        std::cout << std::fixed << std::setprecision(1)
                  << "    covering " << 100.0 * coverage << "% of the model";
        if (dr.notFourSided > 0 || dr.multiArcSides > 0) {
            const double lost = dr.totalArea - dr.coveredArea;
            std::cout << "; " << std::setprecision(1) << 100.0 * lost / dr.totalArea
                      << "% left uncovered by " << (dr.notFourSided + dr.multiArcSides)
                      << " refused face(s)";
        }
        std::cout << std::defaultfloat << "\n";

        {
            int interfaceEdges = 0, badSides = 0, badInterfaces = 0;
            for (size_t b = 0; b < D.blocks.size(); ++b) {
                for (int s = 0; s < 4; ++s) {
                    const int e = D.blocks[b].edges[s];
                    if (e < 0 || e >= static_cast<int>(D.edges.size()) ||
                        D.sidePolyline(static_cast<int>(b), s).size() < 2) ++badSides;
                }
            }
            for (const auto &e : D.edges) {
                if (!e.interface) continue;
                ++interfaceEdges;
                if (e.matLeft == e.matRight) ++badInterfaces;
            }
            if (dr.blocks == lr.faces - lr.badFaces - dr.multiArcSides)
                std::cout << "\033[32m[PASS]\033[0m blockDecomposition() carries every "
                          << "four-sided block (" << dr.blocks << ").\n";
            else {
                std::cout << "\033[31m[FAIL]\033[0m blockDecomposition() has " << dr.blocks
                          << " block(s), the layout has " << (lr.faces - lr.badFaces)
                          << " four-sided face(s) of which " << dr.multiArcSides
                          << " have a split side.\n";
                ++structuralFailures;
            }
            if (badSides == 0)
                std::cout << "\033[32m[PASS]\033[0m Every block side resolves to a macro edge "
                          << "with a polyline.\n";
            else {
                std::cout << "\033[31m[FAIL]\033[0m " << badSides
                          << " block side(s) resolve to nothing.\n";
                ++structuralFailures;
            }
            if (badInterfaces == 0)
                std::cout << "\033[32m[PASS]\033[0m Every interface macro edge separates two "
                          << "materials (" << interfaceEdges << ").\n";
            else {
                std::cout << "\033[31m[FAIL]\033[0m " << badInterfaces
                          << " interface macro edge(s) have one material on both sides.\n";
                ++structuralFailures;
            }
            if (dr.straddlingBlocks > 0) {
                std::cout << "\033[31m[FAIL]\033[0m " << dr.straddlingBlocks
                          << " block(s) sit in more than one material.\n";
                ++structuralFailures;
            }
        }

        // --- The quadrilateral mesh on it, and the smoothing ----------------
        //
        // BlockQuadMesh and mesh::TMOP, which is the same pair MERIDIAN's
        // Stage 10 / Stage 12 and ATLAS's BlockMesh / TMOP run, at the same
        // target: what differs between the three is the layout and not the
        // meshing of it.
        BlockQuadMesh::Options qo;
        qo.targetEdgeLength = targetEdge;
        BlockQuadMesh qmesh(D, qo);
        const BlockQuadMesh::Report &qr = qmesh.getReport();
        std::cout << std::fixed << std::setprecision(4)
                  << "  mesh: " << qr.quads << " quad(s) on " << qr.vertices << " vertices over "
                  << qr.blocks << " block(s), " << qr.chords << " chord(s) at "
                  << qr.minIntervals << "-" << qr.maxIntervals << " edge(s) each";
        if (qr.clampedChords > 0) std::cout << " (" << qr.clampedChords << " clamped)";
        if (qr.unmeshedBlocks > 0) std::cout << ", " << qr.unmeshedBlocks << " block(s) unmeshed";
        std::cout << "\n    edges " << qr.minEdge << " to " << qr.maxEdge << " against a target of "
                  << qr.target << " (rms log ratio " << std::setprecision(3) << qr.edgeRatioRms
                  << ")\n    scaled Jacobian " << std::setprecision(4) << qr.minScaledJacobian
                  << " worst, " << qr.meanScaledJacobian << " mean";
        if (qr.smoothedBlocks > 0)
            std::cout << " (Winslow on " << qr.smoothedBlocks << " folded block(s), "
                      << qr.minScaledJacobianBefore << " worst before)";
        if (qr.invertedQuads > 0) std::cout << "; " << qr.invertedQuads << " element(s) fold";
        std::cout << "\n    boundary nodes " << qr.boundaryNodes << ", within "
                  << qr.boundaryDeviation << " of dS";
        if (qr.materials > 1)
            std::cout << "; " << qr.materials << " material(s) meeting on " << qr.interfaceEdges
                      << " element edge(s), within " << qr.interfaceDeviation
                      << " of the interfaces";
        std::cout << "\n";
        for (const std::string &m : qr.messages) std::cout << "    " << m << "\n";

        // A decomposition with no block in it has no mesh to judge, and
        // saying so is not the same as saying the mesh failed: what went
        // wrong is upstream, and the layout checks above have already said it.
        if (qr.quads == 0) {
            std::cout << "  (nothing to mesh: the layout left no four-sided block)\n";
        } else if (qr.conforming) {
            std::cout << "\033[32m[PASS]\033[0m The mesh is conforming (" << qr.interiorEdges
                      << " interior and " << qr.boundaryEdges << " boundary edge(s), "
                      << qr.cracks << " crack(s)).\n";
        } else {
            std::cout << "\033[31m[FAIL]\033[0m The mesh is not conforming: "
                      << qr.nonManifoldEdges << " edge(s) used a third time, " << qr.cracks
                      << " crack(s).\n";
            ++structuralFailures;
        }
        if (qr.quads > 0) {
            if (qr.valid)
                std::cout << "\033[32m[PASS]\033[0m The mesh validates.\n";
            else
                std::cout << "\033[31m[FAIL]\033[0m The mesh does not validate.\n";
        }

        if (qr.quads > 0) {
            mesh::QuadMesh smoothed = mesh::QuadMesh::from(qmesh);
            mesh::TMOP::Options topt;
            topt.metric = mesh::TMOP::ShapeSize007;
            topt.maxSweeps = 1000;
            // mu at the element corners rather than at the 2x2 Gauss points,
            // which is ATLAS's setting and for ATLAS's reason: these are
            // transfinite grids on a block decomposition, and on one of those
            // a corner can turn over without the barrier at the Gauss points
            // ever seeing it -- the smoother then hands back folds the mesh
            // did not go in with.
            topt.quadrature = mesh::TMOP::Corners;
            mesh::TMOP smoother(smoothed, topt);
            const bool improved = smoother.run();
            const mesh::TMOP::Report &tr = smoother.getReport();
            std::cout << std::fixed << std::setprecision(4)
                      << "  TMOP: " << tr.sweeps << " sweep(s), scaled Jacobian "
                      << tr.minScaledJacobianBefore << " -> " << tr.minScaledJacobianAfter
                      << " worst, " << tr.meanScaledJacobianBefore << " -> "
                      << tr.meanScaledJacobianAfter << " mean";
            if (tr.invertedBefore > 0 || tr.invertedAfter > 0)
                std::cout << ", folds " << tr.invertedBefore << " -> " << tr.invertedAfter;
            std::cout << "\n    nodes " << tr.freeNodes << " free, " << tr.slidingNodes
                      << " sliding, " << tr.fixedNodes << " fixed; area "
                      << (std::fabs(tr.areaAfter - tr.areaBefore) <=
                                  1e-9 * std::fabs(tr.areaBefore)
                              ? "unchanged"
                              : "changed")
                      << "\n";
            if (improved)
                std::cout << "\033[32m[PASS]\033[0m Smoothing left the mesh better than it "
                          << "found it.\n";
            else
                std::cout << "\033[31m[FAIL]\033[0m Smoothing did not improve the mesh.\n";
        }

        const std::string uvFile = stem + "_polysquare.vtu";
        const std::string srcFile = stem + "_source.vtu";
        const std::string blkFile = stem + "_blocks.vtu";
        const std::string trcFile = stem + "_motorcycles.vtu";
        const std::string layFile = stem + "_blocklayout.vtu";
        if (poly.writeVTU(uvFile) && poly.writeSourceVTU(srcFile) &&
            mg.writeBlocksVTU(blkFile) && mg.writeTracesVTU(trcFile) &&
            bl.getLayout().writeFacesVTU(layFile))
            std::cout << "  wrote " << uvFile << " (parameter domain), " << srcFile
                      << " (input domain, uv as a point field),\n         " << blkFile
                      << " (blocks as cell data), " << trcFile << " (the traces) and\n         "
                      << layFile << " (the blocks as polygons)\n";
        else
            std::cout << "\033[31m[FAIL]\033[0m could not write the VTK output.\n";
    } catch (const std::exception &e) {
        std::cout << "\033[31m[FAIL]\033[0m UMBER failed: " << e.what() << "\n";
        return 3;
    }

    return structuralFailures > 0 ? 4 : 0;
}
