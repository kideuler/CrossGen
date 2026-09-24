// Viertel, Osting and Staten, IMR 2019, end to end on one model or many: the
// MBO cross field, its separatrices (Sec. 3), the quad layout with T-junctions
// they cut the model into, the partition simplification of Sec. 4, and then
// what this codebase finishes every method with -- the layout as the shared
// block decomposition on spline geometry (tracing/LayoutBlocks.hxx), meshed by
// mesh/BlockQuadMesh and smoothed by mesh::TMOP.
//
//   TraceMesh [options] <mesh.obj> [more.obj ...]
//   TraceMesh --selftest
//
// Per mesh it prints a line of counts, and the number to read the method by:
// what fraction of the model the blocks cover. A component with a T-junction
// on a side is not meshed yet, so everything short of 100% is the area those
// components (and any that did not come out four-sided at all) take up.
//
// Options:
//   -v                      the per-model detail
//   --no-vtu                write no .vtu files (the default writes the layout,
//                           the simplified layout and the separatrices next to
//                           wherever it is run)
//   --mesh <h>              mesh the blocks at target edge length h
//   --tmop <sweeps>         then smooth that mesh (implies --mesh 0.05 if unset)
//   --step <f> --brep <f>   the B-rep of the blocks (one model only)
//   --blocks-obj <f>        the decomposition's macro edges as OBJ polylines
//   --obj <f> --mfem <f>    the mesh (smoothed if --tmop), one model only
//   --cut-radius r, --tangential deg, --square-rays, --no-simplify, --no-stems,
//   --by-arc-length, --no-parallel-stop, --no-resume, --interior-singularities
//                           the tracing and simplification settings; see
//                           SeparatrixTrace::Settings and PartitionSimplify
//   --pin-open-sides        hold the rim of unmeshed components under TMOP
//   --no-interfaces         trace a multi-material model as if it had one
//                           material (the field not aligned to the interfaces)
//   --no-disk-pin           leave disk inclusions' centres free in the field
//   --dump <prefix>         <prefix>_<model>_{layout,blocks,mesh}.obj and
//                           _tj.txt: the simplified layout, the blocks, the mesh
//                           and the T-junctions, for looking at
//
// The exit code is about soundness, not coverage: every layout a plane graph
// whose faces tile the model, every block side resolving to a macro edge that
// runs between the block's own corners, every Coons patch on its own corners,
// and -- with --mesh -- every mesh conforming. How much of the model is covered
// is what the method achieves on the model and is printed, not asserted.
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "MERIDIAN/Interfaces.hxx"
#include "crossfield/CrossField.hxx"
#include "mesh/BlockQuadMesh.hxx"
#include "mesh/QuadMesh.hxx"
#include "mesh/TMOP.hxx"
#include "tracing/LayoutBlocks.hxx"
#include "tracing/PartitionSimplify.hxx"
#include "tracing/QuadLayout.hxx"
#include "tracing/SeparatrixTrace.hxx"

static const int MBO_MAX_STEPS = 500;

static const char *reasonName(TerminationReason r) {
    switch (r) {
        case TerminationReason::RUNNING:            return "RUNNING";
        case TerminationReason::EXIT_BOUNDARY:      return "boundary";
        case TerminationReason::CROSSED_TWICE:      return "crossed-twice";
        case TerminationReason::CUT_AT_SINGULARITY: return "cut-at-singularity";
        case TerminationReason::HETEROCLINIC:       return "heteroclinic";
        case TerminationReason::LIMIT_CYCLE:        return "limit-cycle";
        case TerminationReason::STUCK:              return "STUCK";
        case TerminationReason::TANGENTIAL:         return "tangential";
    }
    return "?";
}

static std::string stemOf(const std::string &path) {
    size_t slash = path.find_last_of("/\\");
    std::string base = (slash == std::string::npos) ? path : path.substr(slash + 1);
    size_t dot = base.find_last_of('.');
    return (dot == std::string::npos) ? base : base.substr(0, dot);
}

struct Flags {
    bool simplify = true;
    PartitionSimplify::Settings simplifySettings;
    SeparatrixTrace::Settings trace;
    double meshTarget = 0.0;   // 0: no mesh
    int tmopSweeps = 0;        // 0: no smoothing
    std::string step, brep, blocksObj, obj, mfem;
    std::string dump;          // prefix for the debugging dumps
    bool interfaces = true;    // honour material interfaces (--no-interfaces)
    bool pinDisks = true;      // ... and pin disk inclusions' centres (--no-disk-pin)
    bool pinOpenSides = false; // hold the rim of unmeshed components under TMOP
};
static Flags gFlags;

struct Outcome {
    std::string name;
    // Two separate questions. `sound` is whether the layout is a layout at all:
    // a plane graph with no arcs crossing off-node, whose faces tile the model
    // exactly, with every separatrix ending somewhere. That has to hold, and it
    // is what the exit code is about. `allQuads` is whether every component
    // came out four-sided, which is the quality of the partition -- what Sec. 4
    // is for, and not something raw tracing is expected to reach on every
    // model.
    bool sound = false;
    bool allQuads = false;
    bool poincareHopf = false;
    int singularities = 0;
    int separatrices = 0;
    int stuck = 0;
    int dangling = 0;
    int faces = 0;
    int quads = 0;
    int tJunctions = 0;
    double areaRatio = 0.0;
    // after Sec. 4
    int sFaces = 0, sQuads = 0, sTJunctions = 0, sCollapses = 0, sRolledBack = 0;
    double sAreaRatio = 0.0;
    bool sSound = false;
    // the blocks
    int blocks = 0;
    int straddling = 0;        // blocks sitting in more than one material
    double coverage = 0.0;
    bool blocksSound = false;
    // the mesh
    bool meshed = false;
    int meshQuads = 0;
    double meshWorst = 0.0;
    bool meshConforming = true;
    bool smoothed = false;
    double tmopWorst = 0.0;
};

// The raw traced curves, before the layout splits them, as one poly-line each.
static void writeSeparatricesVTU(const std::string &file, const SeparatrixTrace &trace) {
    std::ofstream out(file);
    if (!out) return;
    size_t nPts = 0, nCells = 0;
    for (const auto &s : trace.separatrices) { if (s.path.size() < 2) continue; nPts += s.path.size(); ++nCells; }
    out << "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
        << "  <UnstructuredGrid>\n    <Piece NumberOfPoints=\"" << nPts << "\" NumberOfCells=\"" << nCells << "\">\n";
    out << "      <Points>\n        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto &s : trace.separatrices) { if (s.path.size() < 2) continue;
        for (const auto &p : s.path) out << "          " << p.global_pos[0] << " " << p.global_pos[1] << " 0.0\n"; }
    out << "        </DataArray>\n      </Points>\n";
    out << "      <CellData Scalars=\"id\">\n        <DataArray type=\"Int32\" Name=\"id\" format=\"ascii\">\n";
    for (const auto &s : trace.separatrices) if (s.path.size() >= 2) out << "          " << s.id << "\n";
    out << "        </DataArray>\n        <DataArray type=\"Int32\" Name=\"reason\" format=\"ascii\">\n";
    for (const auto &s : trace.separatrices) if (s.path.size() >= 2) out << "          " << (int)s.termination_reason << "\n";
    out << "        </DataArray>\n        <DataArray type=\"Int32\" Name=\"origin\" format=\"ascii\">\n";
    for (const auto &s : trace.separatrices) if (s.path.size() >= 2)
        out << "          " << (s.originKind == SeparatrixOrigin::Singularity ? s.origin_singularity_id : 100 + s.origin_singularity_id) << "\n";
    out << "        </DataArray>\n      </CellData>\n";
    out << "      <Cells>\n        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    size_t base = 0;
    for (const auto &s : trace.separatrices) { if (s.path.size() < 2) continue;
        out << "         "; for (size_t i = 0; i < s.path.size(); ++i) out << " " << base + i; out << "\n"; base += s.path.size(); }
    out << "        </DataArray>\n        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    size_t off = 0;
    for (const auto &s : trace.separatrices) { if (s.path.size() < 2) continue; off += s.path.size(); out << "          " << off << "\n"; }
    out << "        </DataArray>\n        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < nCells; ++i) out << "          4\n";
    out << "        </DataArray>\n      </Cells>\n    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";
}

// ---------------------------------------------------------------------------
// The structural invariants of the decomposition LayoutBlocks hands on: each
// is a statement the construction makes true, so a failure is a bug and not a
// property of the model. Returns the failures, one line each.
// ---------------------------------------------------------------------------
static std::vector<std::string> checkBlocks(const LayoutBlocks &lb) {
    std::vector<std::string> bad;
    const BlockDecomposition &D = lb.decomposition();
    const LayoutBlocks::Report &r = lb.getReport();

    double lo[2] = {1e300, 1e300}, hi[2] = {-1e300, -1e300};
    for (const auto &v : D.vertices)
        for (int k = 0; k < 2; ++k) { lo[k] = std::min(lo[k], v.p[k]); hi[k] = std::max(hi[k], v.p[k]); }
    const double extent = D.vertices.empty() ? 1.0 : std::max(1e-300, std::hypot(hi[0] - lo[0], hi[1] - lo[1]));

    if (static_cast<int>(D.blocks.size()) != r.blocks)
        bad.push_back("decomposition has " + std::to_string(D.blocks.size()) + " blocks, report " +
                      std::to_string(r.blocks));
    for (size_t b = 0; b < D.blocks.size(); ++b) {
        const auto &blk = D.blocks[b];
        for (int k = 0; k < 4; ++k) {
            const std::vector<Point> side = D.sidePolyline(static_cast<int>(b), k);
            if (side.size() < 2) { bad.push_back("block " + std::to_string(b) + " side " + std::to_string(k) + " has no polyline"); continue; }
            const Point &c0 = D.vertices[blk.corners[k]].p;
            const Point &c1 = D.vertices[blk.corners[(k + 1) % 4]].p;
            if (normP(side.front() - c0) > 1e-12 * extent || normP(side.back() - c1) > 1e-12 * extent)
                bad.push_back("block " + std::to_string(b) + " side " + std::to_string(k) +
                              " does not run between its corners");
        }
    }
    for (size_t e = 0; e < D.edges.size(); ++e) {
        const auto &me = D.edges[e];
        if (me.blockA < 0) { bad.push_back("macro edge " + std::to_string(e) + " has no A block"); continue; }
        if (D.blocks[me.blockA].flip[me.sideA]) bad.push_back("macro edge " + std::to_string(e) + ": A side reads it flipped");
        if (me.blockB >= 0 && !D.blocks[me.blockB].flip[me.sideB])
            bad.push_back("macro edge " + std::to_string(e) + ": B side reads it unflipped");
    }
    if (r.maxCornerGap > 1e-10 * extent)
        bad.push_back("a Coons patch misses its corner by " + std::to_string(r.maxCornerGap));
    if (r.maxSideGap > 1e-9 * extent)
        bad.push_back("a Coons patch misses its side curve by " + std::to_string(r.maxSideGap));
    if (r.coveredArea > r.totalArea * (1.0 + 1e-9)) bad.push_back("blocks cover more than the model");
    return bad;
}

static Outcome processMesh(const std::string &path, bool writeVtu, bool verbose, bool single) {
    Outcome out;
    out.name = stemOf(path);

    std::shared_ptr<Mesh> mesh;
    try {
        mesh = std::make_shared<Mesh>(path);
    } catch (const std::exception &e) {
        std::cerr << "  failed to load " << path << ": " << e.what() << "\n";
        return out;
    }

    // Stage 0b, as the pipelines run it: the material interface network, which
    // the field is aligned to and the tracing honours. Closed loops are not
    // split: a separatrix crossing an inclusion's rim cuts it where the layout
    // actually turns (Interfaces::Options::splitCircleLoops has the argument).
    std::unique_ptr<Interfaces> interfaces;
    if (gFlags.interfaces) {
        Interfaces::Options io;
        io.splitLoops = false;
        interfaces = std::make_unique<Interfaces>(mesh, io);
        if (!interfaces->multiMaterial()) interfaces.reset();
    }

    auto crossField = std::make_shared<CrossField>(mesh);
    if (interfaces) {
        crossField->setAlignedInteriorEdges(interfaces->interfaceEdges());
        crossField->setPinDiskCenters(gFlags.pinDisks);
    }
    crossField->initialize(1, 12345);
    const double nv = static_cast<double>(mesh->vertices.size());
    for (int i = 0; i < MBO_MAX_STEPS; ++i) {
        crossField->step();
        if (crossField->error < 2.0 * nv * 1e-9) break;
    }
    crossField->computeSingularities();

    SeparatrixTrace trace(crossField, true, gFlags.trace, interfaces.get());
    trace.run();
    const SeparatrixTrace::Report &tr = trace.getReport();
    out.poincareHopf = tr.poincareHopf;

    QuadLayout layout(trace);
    layout.build();
    const auto &rep = layout.getReport();

    out.singularities = static_cast<int>(trace.singularities.size());
    out.separatrices = static_cast<int>(trace.separatrices.size());
    for (const auto &s : trace.separatrices) {
        if (s.termination_reason == TerminationReason::STUCK ||
            s.termination_reason == TerminationReason::LIMIT_CYCLE) ++out.stuck;
    }
    out.dangling = rep.danglingEnds;
    out.faces = rep.faces;
    out.quads = rep.quadFaces;
    out.tJunctions = rep.tJunctions;
    out.areaRatio = (rep.meshArea > 0.0) ? rep.totalArea / rep.meshArea : 0.0;
    out.sound = (rep.danglingEnds == 0) && (out.stuck == 0) && (rep.arcCrossings == 0) &&
                (std::fabs(out.areaRatio - 1.0) < 1e-9) && rep.faces > 0;
    out.allQuads = out.sound && rep.badFaces == 0;

    PartitionSimplify simp(layout, gFlags.simplifySettings);
    if (gFlags.simplify) simp.run();
    const auto &sr = simp.getReport();
    const auto &sl = simp.getLayout().getReport();
    out.sFaces = sl.faces;
    out.sQuads = sl.quadFaces;
    out.sTJunctions = sl.tJunctions;
    out.sCollapses = sr.collapses;
    out.sRolledBack = sr.rolledBack;
    out.sAreaRatio = (rep.meshArea > 0.0) ? sl.totalArea / rep.meshArea : 0.0;
    out.sSound = sl.arcCrossings == 0 && sl.danglingEnds == 0 && sl.faces > 0 &&
                 std::fabs(out.sAreaRatio - 1.0) < 1e-9;

    const LayoutBlocks blocks(simp.getLayout(), *mesh);
    const LayoutBlocks::Report &br = blocks.getReport();
    out.blocks = br.blocks;
    out.straddling = br.straddlingBlocks;
    out.coverage = blocks.coverage();
    const std::vector<std::string> blockFailures = checkBlocks(blocks);
    out.blocksSound = blockFailures.empty();

    if (verbose) {
        std::cout << "  mesh          " << mesh->triangles.size() << " triangles, "
                  << mesh->vertices.size() << " vertices\n";
        std::cout << "  cross field   " << trace.singularities.size() << " singularities, error "
                  << crossField->error << "\n";
        std::cout << "  Poincare-Hopf sum d " << tr.interiorIndexSum4 << " + sum(2 - q) "
                  << tr.boundaryIndexSum4 << " = " << (tr.interiorIndexSum4 + tr.boundaryIndexSum4)
                  << " against 4 chi = " << 4 * tr.eulerCharacteristic
                  << (tr.poincareHopf ? " [ok]" : " [FAILS]");
        if (tr.interfaceNodes) std::cout << " (" << tr.interfaceIndexSum4 << " of it from "
                                         << tr.interfaceNodes << " interface node(s))";
        if (tr.multipleSingularities) std::cout << "; " << tr.multipleSingularities << " with |d| >= 2";
        if (tr.droppedSingularities) std::cout << "; " << tr.droppedSingularities << " without ports";
        std::cout << "\n";
        int byReason[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        for (const auto &s : trace.separatrices) byReason[static_cast<int>(s.termination_reason)]++;
        std::cout << "  separatrices  " << trace.separatrices.size() << " (";
        bool first = true;
        for (int i = 0; i < 8; ++i) {
            if (!byReason[i]) continue;
            if (!first) std::cout << ", ";
            std::cout << byReason[i] << " " << reasonName(static_cast<TerminationReason>(i));
            first = false;
        }
        std::cout << "); " << tr.heteroclinicJoins << " join(s), " << tr.resumed << " resumed, "
                  << tr.crossingsRetracted << " crossing(s) retracted, "
                  << tr.tangentialBoundaryExits << " tangential landing(s)\n";
        if (interfaces) {
            const Interfaces::Report &ir = interfaces->getReport();
            std::cout << "  interfaces    " << ir.materials << " materials, " << ir.interfaceEdges
                      << " interface edges in " << ir.branches << " branch(es) (" << ir.closedLoops
                      << " closed), " << ir.nodes << " node(s): " << ir.junctions << " junction, "
                      << ir.landings << " landing, " << ir.kinks << " kink, " << ir.illPosedNodes
                      << " ill-posed; field pinned at " << crossField->alignedInterfaceVertices()
                      << " interface vertices, free at " << crossField->freeInterfaceVertices()
                      << ", " << crossField->pinnedDiskCenters() << " disk centre(s) pinned"
                      << "\n                " << tr.interfaceEmitted << " separatrice(s) from the nodes, "
                      << tr.absorbedSingularities << " singular triangle(s) absorbed by them, "
                      << tr.interfaceCrossings << " crossing(s) of an interface (" << tr.interfaceGrazes
                      << " shallower than the tangential angle)\n";
        }
        std::cout << "  corners       " << trace.boundaryCorners.size() << " on the boundary:";
        for (const auto &c : trace.boundaryCorners)
            std::cout << " " << (int)std::lround(c.interiorAngle * 180.0 / M_PI) << "(q" << c.quarters << ")";
        std::cout << "\n";
        for (size_t i = 0; i < trace.singularities.size() && i < 12; ++i) {
            const auto &s = trace.singularities[i];
            std::cout << "    sing " << i << " index " << s.singularityIndex << " at ("
                      << s.coordinates[0] << ", " << s.coordinates[1] << ") residual "
                      << (s.portResidual * 180.0 / M_PI) << " deg, ports";
            for (double a : s.portAngles) std::cout << " " << (int)std::lround(a * 180.0 / M_PI);
            std::cout << "\n";
        }
        std::cout << "  layout        " << rep.nodes << " nodes, " << rep.arcs << " arcs, "
                  << rep.faces << " faces, " << rep.tJunctions << " T-junctions, "
                  << rep.danglingEnds << " loose ends, " << rep.arcCrossings
                  << " un-noded arc crossings, " << rep.unboundedCycles << " outer cycles, "
                  << rep.skewNodes << " skew nodes\n";
        std::cout << "  faces         " << rep.quadFaces << " four-cornered, " << rep.badFaces
                  << " not; area " << std::setprecision(10) << out.areaRatio
                  << " of the model\n" << std::setprecision(6);
        std::cout << "  simplified    " << sl.faces << " components (" << sl.quadFaces
                  << " four-sided), " << sl.tJunctions << " T-junctions, after " << sr.collapses
                  << " chord collapse(s); " << sr.blockedByConditions << " chord(s) blocked by the "
                     "conditions, " << sr.blockedByEnergy << " by the energy, " << sr.rolledBack
                  << " rolled back; area " << std::setprecision(10) << out.sAreaRatio << "\n"
                  << std::setprecision(6);
        if (sr.stems.attempted)
            std::cout << "  stems (Sec 12) " << sr.stems.extended << " of " << sr.stems.attempted
                      << " T-junction stem(s) traced on to the boundary (" << sr.tJunctionsBeforeStems
                      << " -> " << sl.tJunctions << " T-junctions); not: " << sr.stems.noBoundary
                      << " never reached it, " << sr.stems.tangential << " ran alongside a separatrix, "
                      << sr.stems.throughNode << " ran into a node, " << sr.stems.tooManyCrossings
                      << " crossed too many, " << sr.stems.invalid << " left a worse layout\n";
        if (sr.rolledBack)
            std::cout << "  rollbacks     " << sr.rbFailed << " could not be applied, "
                      << sr.rbCrossings << " non-planar, " << sr.rbDangling << " loose ends, "
                      << sr.rbNotFewer << " no fewer components, " << sr.rbMoreT
                      << " more T-junctions, " << sr.rbSing << " lost a singularity, "
                      << sr.rbArea << " changed area, " << sr.rbWorse << " more bad faces\n";
        if (rep.badFaces) {
            int shown = 0;
            for (size_t i = 0; i < layout.getFaces().size() && shown < 8; ++i) {
                const auto &f = layout.getFaces()[i];
                if (f.corners == 4) continue;
                std::cout << "      face " << i << ": " << f.corners << " corners, "
                          << f.darts.size() << " arcs, area " << f.area << " |";
                static const char *kindName[] = {"sing", "corn", "bhit", "xing", "hetc", "tjct", "loose", "inode", "ihit"};
                for (size_t k = 0; k < f.nodes.size(); ++k)
                    std::cout << " " << (f.isCorner[k] ? "*" : "-")
                              << kindName[static_cast<int>(layout.getNodes()[f.nodes[k]].kind)]
                              << (int)std::lround(f.turn[k] * 180.0 / M_PI)
                              << "/" << layout.getNodes()[f.nodes[k]].darts.size();
                std::cout << "\n";
                ++shown;
            }
        }
        std::cout << std::fixed << std::setprecision(1)
                  << "  blocks        " << br.blocks << " of " << br.faces << " components, covering "
                  << 100.0 * out.coverage << "% of the model; refused: " << br.tJunctionFaces
                  << " for a T-junction on a side (" << br.tJunctions << " T-junction node(s), "
                  << 100.0 * br.tJunctionArea / std::max(1e-300, br.totalArea) << "%), "
                  << br.notFourSided << " not four-sided ("
                  << 100.0 * br.notFourSidedArea / std::max(1e-300, br.totalArea) << "%), "
                  << br.notDisks << " not disks\n"
                  << std::setprecision(3)
                  << "  splines       " << br.macroArcs << " macro arc(s): " << br.fittedArcs
                  << " fitted (up to " << br.maxSegmentsUsed << " segments, deviation "
                  << br.meanDeviation << " mean / " << br.maxDeviation << " worst, in mean edges";
        if (br.unconverged) std::cout << ", " << br.unconverged << " above tolerance";
        std::cout << "), " << br.exactArcs << " exact; " << br.foldedPatches
                  << " folded patch(es); corner gap " << std::scientific << std::setprecision(1)
                  << br.maxCornerGap << ", side gap " << br.maxSideGap << std::fixed << "\n"
                  << "  B-rep         " << br.brepFaces << " faces, " << br.brepEdges << " edges ("
                  << br.brepSharedEdges << " shared, " << br.brepFreeEdges << " free), "
                  << (br.brepValid ? "valid" : "NOT valid");
        if (br.brepFailures) std::cout << ", " << br.brepFailures << " refused by the kernel";
        if (br.materials > 1 || br.straddlingBlocks)
            std::cout << "; " << br.materials << " material(s), " << br.straddlingBlocks
                      << " block(s) straddling an interface";
        std::cout << "\n" << std::defaultfloat << std::setprecision(6);
        for (const std::string &m : br.messages) std::cout << "    " << m << "\n";
    }
    for (const std::string &m : blockFailures) std::cout << "  [FAIL] " << m << "\n";

    if (!gFlags.dump.empty()) {
        const std::string pre = gFlags.dump + "_" + out.name;
        std::ofstream lo(pre + "_layout.obj");
        int base = 1;
        for (const auto &a : simp.getLayout().getArcs()) {
            if (a.pts.size() < 2) continue;
            for (const Point &p : a.pts) lo << "v " << p[0] << " " << p[1] << " 0\n";
            lo << "l";
            for (size_t i = 0; i < a.pts.size(); ++i) lo << " " << base + static_cast<int>(i);
            lo << "\n";
            base += static_cast<int>(a.pts.size());
        }
        blocks.decomposition().writeEdgesOBJ(pre + "_blocks.obj");
        std::ofstream tj(pre + "_tj.txt");
        for (const Point &p : blocks.tJunctionPoints()) tj << p[0] << " " << p[1] << "\n";
    }

    if (single && !gFlags.step.empty())
        std::cout << "  STEP " << gFlags.step << (blocks.writeSTEP(gFlags.step) ? "" : " FAILED") << "\n";
    if (single && !gFlags.brep.empty())
        std::cout << "  BREP " << gFlags.brep << (blocks.writeBREP(gFlags.brep) ? "" : " FAILED") << "\n";
    if (single && !gFlags.blocksObj.empty())
        blocks.decomposition().writeEdgesOBJ(gFlags.blocksObj);

    // --- the mesh, and the smoothing ---------------------------------------
    //
    // BlockQuadMesh and mesh::TMOP, the pair UMBER runs, at a target edge
    // length that means what it means to every other method. TMOP samples mu
    // at the element corners: these are transfinite grids on a block
    // decomposition, where a corner can turn over without the Gauss points
    // seeing it (ATLAS's and UMBER's setting, for that reason).
    if (gFlags.meshTarget > 0.0 && br.blocks > 0) {
        BlockQuadMesh::Options qo;
        qo.targetEdgeLength = gFlags.meshTarget;
        const BlockQuadMesh qm(blocks.decomposition(), qo);
        const BlockQuadMesh::Report &qr = qm.getReport();
        out.meshed = true;
        out.meshQuads = qr.quads;
        out.meshWorst = qr.minScaledJacobian;
        out.meshConforming = qr.conforming;
        if (verbose || !qr.conforming) {
            std::cout << std::fixed << std::setprecision(4)
                      << "  mesh          " << qr.quads << " quads on " << qr.vertices << " vertices, "
                      << qr.chords << " chord(s) at " << qr.minIntervals << "-" << qr.maxIntervals
                      << " edges; scaled Jacobian " << qr.minScaledJacobian << " worst, "
                      << qr.meanScaledJacobian << " mean";
            if (qr.invertedQuads) std::cout << ", " << qr.invertedQuads << " folded";
            std::cout << "; " << (qr.conforming ? "conforming" : "NOT conforming") << "\n"
                      << std::defaultfloat << std::setprecision(6);
            for (const std::string &m : qr.messages) std::cout << "    " << m << "\n";
        }
        mesh::QuadMesh qmOut = mesh::QuadMesh::from(qm);
        // --pin-open-sides holds the rim of every unmeshed component, where a
        // mesh of that component would one day have to meet this one. Off by
        // default because it costs the smoother measurably: over the corpus
        // at h = 0.05 and 1000 sweeps the worst scaled Jacobian after TMOP is
        // lower pinned than sliding on every partially covered model where
        // the two differ (geom031 0.32 against 0.74, geom010 0.16 against
        // 0.29, geom007 0.69 against 0.86) -- the ring of elements along a
        // pinned rim cannot redistribute and absorbs the whole mismatch.
        if (gFlags.pinOpenSides && !qm.openSideVertices().empty()) {
            for (const int v : qm.openSideVertices()) qmOut.pinVertex(v);
            qmOut.classifyNodes();
            qmOut.buildFeatureCurves();
            qmOut.computeSlideTangents();
        }
        if (gFlags.tmopSweeps > 0 && qr.quads > 0) {
            mesh::TMOP::Options topt;
            topt.metric = mesh::TMOP::ShapeSize007;
            topt.maxSweeps = gFlags.tmopSweeps;
            topt.quadrature = mesh::TMOP::Corners;
            mesh::TMOP smoother(qmOut, topt);
            smoother.run();
            const mesh::TMOP::Report &tm = smoother.getReport();
            out.smoothed = true;
            out.tmopWorst = tm.minScaledJacobianAfter;
            if (verbose)
                std::cout << std::fixed << std::setprecision(4) << "  TMOP          " << tm.sweeps
                          << " sweep(s), scaled Jacobian " << tm.minScaledJacobianBefore << " -> "
                          << tm.minScaledJacobianAfter << " worst, " << tm.meanScaledJacobianBefore
                          << " -> " << tm.meanScaledJacobianAfter << " mean, folds "
                          << tm.invertedBefore << " -> " << tm.invertedAfter << "\n"
                          << std::defaultfloat << std::setprecision(6);
        }
        if (single && !gFlags.obj.empty()) qmOut.writeOBJ(gFlags.obj);
        if (!gFlags.dump.empty()) qmOut.writeOBJ(gFlags.dump + "_" + out.name + "_mesh.obj");
        if (single && !gFlags.mfem.empty()) qmOut.writeMFEM(gFlags.mfem);
    }

    if (writeVtu) {
        layout.writeArcsVTU(out.name + "_layout_arcs.vtu");
        layout.writeFacesVTU(out.name + "_layout_faces.vtu");
        layout.writeNodesVTU(out.name + "_layout_nodes.vtu");
        simp.getLayout().writeArcsVTU(out.name + "_simple_arcs.vtu");
        simp.getLayout().writeFacesVTU(out.name + "_simple_faces.vtu");
        simp.getLayout().writeNodesVTU(out.name + "_simple_nodes.vtu");
        writeSeparatricesVTU(out.name + "_separatrices.vtu", trace);
    }
    return out;
}

// ===========================================================================
// --selftest: docs/viertel_2019.md Sec. 11, on fields written down in closed
// form over triangulated grids, so that every expected number is known exactly
// rather than read off a previous run.
// ===========================================================================
namespace selftest {

int failures = 0;

void check(bool ok, const std::string &what) {
    std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << what << "\n";
    if (!ok) ++failures;
}

// An n x n grid of unit-square cells scaled to [0, 1]^2, each cell cut along
// alternating diagonals, keeping the cells `keep` accepts (by cell centre),
// each triangle of material `material(cell centre)` when that is given.
std::shared_ptr<Mesh> grid(int n, const std::function<bool(double, double)> &keep,
                           const std::function<int(double, double)> &material = nullptr) {
    std::vector<Point> V;
    std::vector<int> id((n + 1) * (n + 1), -1);
    std::vector<Triangle> T;
    std::vector<int> mat;
    auto vert = [&](int i, int j) {
        int &v = id[j * (n + 1) + i];
        if (v < 0) { v = static_cast<int>(V.size()); V.push_back(Point{double(i) / n, double(j) / n}); }
        return v;
    };
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            if (!keep((i + 0.5) / n, (j + 0.5) / n)) continue;
            const int a = vert(i, j), b = vert(i + 1, j), c = vert(i + 1, j + 1), d = vert(i, j + 1);
            if ((i + j) % 2 == 0) { T.push_back({a, b, c}); T.push_back({a, c, d}); }
            else { T.push_back({a, b, d}); T.push_back({b, c, d}); }
            const int m = material ? material((i + 0.5) / n, (j + 0.5) / n) : 1;
            mat.push_back(m);
            mat.push_back(m);
        }
    }
    if (!material) return std::make_shared<Mesh>(V, T);
    return std::make_shared<Mesh>(V, T, mat);
}

std::shared_ptr<CrossField> field(const std::shared_ptr<Mesh> &m,
                                  const std::function<double(const Point &)> &theta) {
    auto cf = std::make_shared<CrossField>(m);
    cf->u_k.resize(static_cast<Eigen::Index>(m->vertices.size()));
    for (size_t v = 0; v < m->vertices.size(); ++v)
        cf->u_k[static_cast<Eigen::Index>(v)] = std::polar(1.0, 4.0 * theta(m->vertices[v]));
    cf->u_k_prev = cf->u_k;
    return cf;
}

Point barycentre(const Mesh &m, int f) {
    const Triangle &t = m.triangles[f];
    return (m.vertices[t[0]] + m.vertices[t[1]] + m.vertices[t[2]]) / 3.0;
}

int triangleNear(const Mesh &m, const Point &q) {
    int best = 0;
    double bd = 1e300;
    for (size_t f = 0; f < m.triangles.size(); ++f) {
        const double d = normP(barycentre(m, static_cast<int>(f)) - q);
        if (d < bd) { bd = d; best = static_cast<int>(f); }
    }
    return best;
}

double angleGap(double a, double b) { return std::fabs(wrap_pi(a - b)); }

// Sec. 11 tests 1-3: the index, the ports, the hyperbola.
void singularityModel() {
    std::cout << "singularity model: theta = (d/4) atan2(y - y0, x - x0) + theta0\n";
    const auto m = grid(16, [](double, double) { return true; });
    const int f0 = triangleNear(*m, Point{0.53, 0.47});
    const Point a = barycentre(*m, f0);
    const double theta0 = 0.3;
    for (const int d : {1, -1, 2}) {
        auto cf = field(m, [&](const Point &p) {
            return 0.25 * d * std::atan2(p[1] - a[1], p[0] - a[0]) + theta0;
        });
        FieldTracer tracer(cf, false);   // barycentre: the centre the field was written about
        const auto &S = tracer.getSingularities();
        if (std::abs(d) >= 2) {
            // Not the spec's test 2 as written. Under principal matching an
            // index-1/2 singularity cannot sit in one triangle: the three
            // edges subtend 2 pi between them from any point inside, so one
            // subtends at least 2 pi / 3 and the field turns by a third of pi
            // along it -- past the pi/4 the matching can represent. The index
            // then spreads over the triangles round the centre, and what is
            // exact is its total, and that it stays there.
            int sum = 0;
            double far = 0.0;
            for (const Singularity &s : S) {
                sum += s.d;
                far = std::max(far, normP(barycentre(*m, s.triangleIndex) - a));
            }
            check(sum == d && far < 2.0 / 16.0,
                  "d = " + std::to_string(d) + ": the index spreads over " + std::to_string(S.size()) +
                      " triangles round the centre and adds up to d/4");
            continue;
        }
        check(S.size() == 1 && S[0].triangleIndex == f0 && S[0].d == d,
              "d = " + std::to_string(d) + ": exactly the triangle holding the centre is singular, with index d/4");
        if (S.size() != 1) continue;

        // Ports: the cross is radial where (d/4) phi + theta0 = phi mod pi/2.
        const int n = 4 - d;
        double worst = 0.0;
        for (int k = 0; k < n; ++k) {
            const double want = 4.0 * theta0 / n + 2.0 * M_PI * k / n;
            double best = 1e300;
            for (double got : S[0].portAngles) best = std::min(best, angleGap(got, want));
            worst = std::max(worst, best);
        }
        check(S[0].numPorts() == n && worst < 1e-9,
              "d = " + std::to_string(d) + ": the " + std::to_string(n) +
              " ports are the true separatrix directions (worst " + std::to_string(worst) + " rad)");

        // Hyperbola: a streamline entered near a corner of the singular
        // triangle, its samples on r^p cos(p (phi - sigma_k)) = c^p for one k.
        const Triangle &t = m->triangles[f0];
        const Point q = m->vertices[t[0]] * 0.7 + m->vertices[t[1]] * 0.2 + m->vertices[t[2]] * 0.1;
        const double dir = 0.25 * d * std::atan2(q[1] - a[1], q[0] - a[0]) + theta0;
        Walker w = tracer.startAt(f0, q, dir);
        std::vector<TracePoint> path{TracePoint{q, f0, dir}};
        tracer.advance(w, path);
        const double p = n / 4.0;
        double bestSpread = 1e300;
        for (double sigma : S[0].portAngles) {
            double lo = 1e300, hi = -1e300;
            for (size_t i = 0; i + 1 < path.size(); ++i) {   // the last point is a chord's exit
                const Point r = path[i].global_pos - a;
                const double inv = std::pow(normP(r), p) * std::cos(p * wrap_pi(std::atan2(r[1], r[0]) - sigma));
                lo = std::min(lo, inv);
                hi = std::max(hi, inv);
            }
            if (path.size() >= 3 && lo > 0.0) bestSpread = std::min(bestSpread, (hi - lo) / hi);
        }
        check(bestSpread < 1e-9, "d = " + std::to_string(d) + ": the sweep's " +
                                     std::to_string(path.size()) + " samples lie on one hyperbola (spread " +
                                     std::to_string(bestSpread) + ")");
    }
}

// Sec. 11 test 5: a constant field traces straight lines.
void constantField() {
    std::cout << "constant field: streamlines are straight\n";
    const auto m = grid(12, [](double, double) { return true; });
    const double theta0 = 0.2;
    auto cf = field(m, [&](const Point &) { return theta0; });
    FieldTracer tracer(cf, false);
    check(tracer.getSingularities().empty(), "no singular triangle");
    const Point start{0.31, 0.27};
    Walker w = tracer.startAt(triangleNear(*m, start), start, theta0);
    std::vector<TracePoint> path{TracePoint{start, w.tri, theta0}};
    FieldTracer::Status st = FieldTracer::Status::Ok;
    for (int i = 0; i < 1000 && st == FieldTracer::Status::Ok; ++i) st = tracer.advance(w, path);
    double off = 0.0;
    const Point dir{std::cos(theta0), std::sin(theta0)};
    for (const TracePoint &tp : path) off = std::max(off, std::fabs(cross2(dir, tp.global_pos - start)));
    check(st == FieldTracer::Status::Boundary && off < 1e-12,
          "reaches the boundary on a straight line (" + std::to_string(path.size()) + " points, " +
              "worst offset " + std::to_string(off) + ")");
}

// Sec. 11 tests 7 and 8, then the T-junction case this codebase adds: a
// component with a T-junction on a side is not a block.
void layouts() {
    std::cout << "square, constant field: one block\n";
    {
        const auto m = grid(8, [](double, double) { return true; });
        SeparatrixTrace trace(field(m, [](const Point &) { return 0.0; }), false);
        trace.run();
        check(trace.getReport().poincareHopf, "Poincare-Hopf holds");
        QuadLayout layout(trace);
        layout.build();
        const LayoutBlocks lb(layout, *m);
        check(trace.separatrices.empty() && lb.getReport().blocks == 1 && lb.coverage() > 1.0 - 1e-12,
              "no separatrices, one block, all of the model covered");
        check(checkBlocks(lb).empty(), "the decomposition is sound");
    }

    std::cout << "L-shape, constant field: the reflex corner emits two, three blocks\n";
    {
        const auto m = grid(10, [](double x, double y) { return x < 0.5 || y < 0.5; });
        SeparatrixTrace trace(field(m, [](const Point &) { return 0.0; }), false);
        trace.run();
        check(trace.getReport().poincareHopf, "Poincare-Hopf holds");
        QuadLayout layout(trace);
        layout.build();
        const LayoutBlocks lb(layout, *m);
        const auto &r = lb.getReport();
        check(trace.separatrices.size() == 2 && layout.getReport().faces == 3 &&
                  layout.getReport().tJunctions == 0,
              "two separatrices, three components, no T-junction");
        check(r.blocks == 3 && lb.coverage() > 1.0 - 1e-12 && r.tJunctions == 0,
              "three blocks covering the model");
        check(checkBlocks(lb).empty(), "the decomposition is sound");
        BlockQuadMesh::Options qo;
        qo.targetEdgeLength = 0.1;
        const BlockQuadMesh qm(lb.decomposition(), qo);
        check(qm.getReport().conforming && qm.getReport().invertedQuads == 0 && qm.getReport().quads == 75,
              "the mesh is conforming, unfolded, 75 quads (got " + std::to_string(qm.getReport().quads) + ")");
    }

    std::cout << "a T-junction: the component with it on a side is not a block\n";
    {
        using NK = QuadLayout::NodeKind;
        const std::vector<std::pair<Point, NK>> spec = {
            {{0, 0}, NK::BoundaryCorner}, {{1, 0}, NK::BoundaryCorner}, {{1, 1}, NK::BoundaryCorner},
            {{0, 1}, NK::BoundaryCorner}, {{0.5, 0}, NK::BoundaryHit}, {{0.5, 1}, NK::BoundaryHit},
            {{1, 0.5}, NK::BoundaryHit}, {{0.5, 0.5}, NK::TJunction}};
        std::vector<QuadLayout::Node> nodes;
        for (const auto &[p, k] : spec) { QuadLayout::Node n; n.pos = p; n.kind = k; nodes.push_back(n); }
        std::vector<QuadLayout::Arc> arcs;
        auto arc = [&](int a, int b, bool boundary) {
            QuadLayout::Arc c;
            c.pts = {nodes[a].pos, (nodes[a].pos + nodes[b].pos) * 0.5, nodes[b].pos};
            c.a = a; c.b = b; c.onBoundary = boundary; c.separatrix = boundary ? -1 : 0;
            c.length = normP(nodes[b].pos - nodes[a].pos);
            arcs.push_back(c);
        };
        for (const auto &[a, b] : std::vector<std::pair<int, int>>{{0, 4}, {4, 1}, {1, 6}, {6, 2}, {2, 5}, {5, 3}, {3, 0}})
            arc(a, b, true);
        arc(4, 7, false); arc(7, 5, false); arc(6, 7, false);
        QuadLayout layout(1e-9);
        layout.rebuild(nodes, arcs);
        const auto m = grid(8, [](double, double) { return true; });
        const LayoutBlocks lb(layout, *m);
        const auto &r = lb.getReport();
        check(r.faces == 3 && r.blocks == 2 && r.tJunctionFaces == 1,
              "three components, two blocks, one refused for its T-junction");
        check(r.tJunctions == 1 && normP(lb.tJunctionPoints()[0] - Point{0.5, 0.5}) < 1e-12,
              "the T-junction is reported where it is");
        check(std::fabs(lb.coverage() - 0.5) < 1e-12, "half the model covered");
        check(checkBlocks(lb).empty(), "the decomposition is sound");
        BlockQuadMesh::Options qo;
        qo.targetEdgeLength = 0.125;
        const BlockQuadMesh qm(lb.decomposition(), qo);
        check(qm.getReport().conforming && qm.getReport().quads == 32,
              "the two blocks mesh conformingly (" + std::to_string(qm.getReport().quads) + " quads)");
    }
}

// Material interfaces: the P1 MBO field aligned to them, the network's nodes
// emitting by their sectors, and the blocks each in one material.
void interfaces() {
    auto traced = [](const std::shared_ptr<Mesh> &m, int &blocks, double &coverage, int &straddling,
                     int &materials, int &separatrices, int &singularities, bool &ph, bool &sound) {
        Interfaces::Options io;
        io.splitLoops = false;
        Interfaces itf(m, io);
        auto cf = std::make_shared<CrossField>(m);
        cf->setAlignedInteriorEdges(itf.interfaceEdges());
        cf->initialize(1, 7);
        for (int i = 0; i < 200; ++i) { cf->step(); if (cf->error < 1e-10) break; }
        SeparatrixTrace trace(cf, true, SeparatrixTrace::Settings(), &itf);
        trace.run();
        QuadLayout layout(trace);
        layout.build();
        PartitionSimplify simp(layout);
        simp.run();
        const LayoutBlocks lb(simp.getLayout(), *m);
        blocks = lb.getReport().blocks;
        coverage = lb.coverage();
        straddling = lb.getReport().straddlingBlocks;
        materials = lb.getReport().materials;
        separatrices = static_cast<int>(trace.separatrices.size());
        singularities = static_cast<int>(trace.singularities.size());
        ph = trace.getReport().poincareHopf;
        sound = layout.getReport().arcCrossings == 0 && layout.getReport().danglingEnds == 0 &&
                checkBlocks(lb).empty();
    };
    int blocks = 0, straddling = 0, materials = 0, seps = 0, sings = 0;
    double coverage = 0.0;
    bool ph = false, sound = false;

    std::cout << "two materials across a straight interface: two blocks\n";
    traced(grid(8, [](double, double) { return true; }, [](double x, double) { return x < 0.5 ? 1 : 2; }),
           blocks, coverage, straddling, materials, seps, sings, ph, sound);
    check(sings == 0 && seps == 0 && ph, "an aligned field with no singularity, nothing emitted, Poincare-Hopf holds");
    check(blocks == 2 && coverage > 1.0 - 1e-12 && straddling == 0 && materials == 2 && sound,
          "two blocks, one per material, covering the model (" + std::to_string(blocks) + ")");

    std::cout << "an orthogonal T of three materials: its straight sector emits one\n";
    traced(grid(8, [](double, double) { return true; },
                [](double x, double y) { return x < 0.5 ? 1 : (y < 0.5 ? 2 : 3); }),
           blocks, coverage, straddling, materials, seps, sings, ph, sound);
    check(sings == 0 && seps == 1 && ph, "one separatrix, from the junction into its half-turn sector (" +
                                             std::to_string(seps) + ")");
    check(blocks == 4 && coverage > 1.0 - 1e-12 && straddling == 0 && materials == 3 && sound,
          "four blocks, none in two materials, covering the model (" + std::to_string(blocks) + ")");
}

int run() {
    singularityModel();
    constantField();
    layouts();
    interfaces();
    std::cout << (failures ? "\n" + std::to_string(failures) + " check(s) FAILED\n" : "\nall checks passed\n");
    return failures ? 1 : 0;
}

} // namespace selftest

int main(int argc, char **argv) {
    std::vector<std::string> meshes;
    bool writeVtu = true;
    bool verbose = false;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&]() -> const char * { return (i + 1 < argc) ? argv[++i] : ""; };
        if (a == "--selftest") return selftest::run();
        else if (a == "--no-vtu") writeVtu = false;
        else if (a == "--cut-radius") gFlags.trace.singularityCutRadius = std::atof(next());
        else if (a == "--tangential") gFlags.trace.tangentialAngle = std::atof(next()) * M_PI / 180.0;
        else if (a == "--square-rays") gFlags.trace.evenCornerRays = false;
        else if (a == "--no-simplify") gFlags.simplify = false;
        else if (a == "--no-stems") gFlags.simplifySettings.extendStems = false;
        else if (a == "--stem-crossings") gFlags.simplifySettings.stems.maxCrossings = std::atoi(next());
        else if (a == "--by-arc-length") gFlags.trace.growByArcLength = true;
        else if (a == "--no-parallel-stop") gFlags.trace.stopParallelTangential = false;
        else if (a == "--no-resume") gFlags.trace.resumeAfterTruncation = false;
        else if (a == "--interior-singularities") gFlags.trace.singularitiesAtBoundary = false;
        else if (a == "--dump") gFlags.dump = next();
        else if (a == "--pin-open-sides") gFlags.pinOpenSides = true;
        else if (a == "--no-interfaces") gFlags.interfaces = false;
        else if (a == "--no-disk-pin") gFlags.pinDisks = false;
        else if (a == "--mesh") gFlags.meshTarget = std::atof(next());
        else if (a == "--tmop") gFlags.tmopSweeps = std::atoi(next());
        else if (a == "--step") gFlags.step = next();
        else if (a == "--brep") gFlags.brep = next();
        else if (a == "--blocks-obj") gFlags.blocksObj = next();
        else if (a == "--obj") gFlags.obj = next();
        else if (a == "--mfem") gFlags.mfem = next();
        else if (a == "-v") verbose = true;
        else meshes.push_back(a);
    }
    if (gFlags.tmopSweeps > 0 && gFlags.meshTarget <= 0.0) gFlags.meshTarget = 0.05;
    if (meshes.empty()) {
        std::cerr << "Usage: " << argv[0] << " [-v] [--no-vtu] [--mesh h] [--tmop sweeps] "
                     "[--step f] [--brep f] <mesh.obj> [more.obj ...]\n       "
                  << argv[0] << " --selftest\n";
        return 1;
    }
    const bool single = meshes.size() == 1;
    if (!single && !(gFlags.step.empty() && gFlags.brep.empty() && gFlags.blocksObj.empty() &&
                     gFlags.obj.empty() && gFlags.mfem.empty()))
        std::cerr << "note: --step/--brep/--blocks-obj/--obj/--mfem are written for one model only\n";

    std::vector<Outcome> all;
    for (const auto &m : meshes) {
        if (!single || verbose) std::cout << stemOf(m) << "\n";
        all.push_back(processMesh(m, writeVtu, verbose || single, single));
    }

    const bool meshed = gFlags.meshTarget > 0.0;
    std::cout << "\n"
              << std::left << std::setw(12) << "mesh" << std::right << std::setw(5) << "sing"
              << std::setw(3) << "PH" << std::setw(6) << "seps" << std::setw(7) << "faces"
              << std::setw(7) << "quads" << std::setw(6) << "T-jct" << std::setw(6) << "loose"
              << "  |" << std::setw(5) << "coll" << std::setw(6) << "faces" << std::setw(6)
              << "quads" << std::setw(6) << "T-jct" << std::setw(5) << "back"
              << "  |" << std::setw(7) << "blocks" << std::setw(8) << "cover%";
    if (meshed) std::cout << "  |" << std::setw(7) << "quads" << std::setw(8) << "minSJ"
                          << (gFlags.tmopSweeps > 0 ? "    TMOP" : "");
    std::cout << "  sound\n";
    int sound = 0, allQuads = 0, faces = 0, quads = 0, sSound = 0, sFaces = 0, sQuads = 0;
    int collapses = 0, rolled = 0, blocks = 0, tAfter = 0, ph = 0, full = 0;
    int straddling = 0, straddlingModels = 0;
    double coverSum = 0.0;
    for (const auto &o : all) {
        const bool ok = o.sound && o.sSound && o.blocksSound && o.meshConforming;
        std::cout << std::left << std::setw(12) << o.name << std::right << std::setw(5)
                  << o.singularities << std::setw(3) << (o.poincareHopf ? "y" : "N") << std::setw(6)
                  << o.separatrices << std::setw(7) << o.faces << std::setw(7) << o.quads
                  << std::setw(6) << o.tJunctions << std::setw(6) << (o.stuck + o.dangling) << "  |"
                  << std::setw(5) << o.sCollapses << std::setw(6) << o.sFaces << std::setw(6)
                  << o.sQuads << std::setw(6) << o.sTJunctions << std::setw(5) << o.sRolledBack
                  << "  |" << std::setw(7) << o.blocks << std::fixed << std::setprecision(1)
                  << std::setw(8) << 100.0 * o.coverage;
        if (meshed) {
            std::cout << "  |" << std::setw(7) << o.meshQuads << std::setprecision(3) << std::setw(8)
                      << o.meshWorst;
            if (gFlags.tmopSweeps > 0) std::cout << std::setw(8) << o.tmopWorst;
        }
        std::cout << std::defaultfloat << (ok ? "    yes" : "     NO") << "\n";
        if (o.sound) ++sound;
        if (o.allQuads) ++allQuads;
        if (ok) ++sSound;
        if (o.poincareHopf) ++ph;
        if (o.coverage > 1.0 - 1e-9) ++full;
        faces += o.faces;
        quads += o.quads;
        sFaces += o.sFaces;
        sQuads += o.sQuads;
        collapses += o.sCollapses;
        rolled += o.sRolledBack;
        blocks += o.blocks;
        tAfter += o.sTJunctions;
        coverSum += o.coverage;
        straddling += o.straddling;
        if (o.straddling > 0) ++straddlingModels;
    }
    const int n = static_cast<int>(all.size());
    std::cout << std::fixed << std::setprecision(1) << sound << "/" << n
              << " layouts sound before simplification, " << sSound << "/" << n
              << " sound through to the blocks" << (meshed ? " and the mesh" : "") << "; " << allQuads
              << "/" << n << " all four-sided before; Poincare-Hopf holds on " << ph << "/" << n << "\n"
              << "components " << faces << " -> " << sFaces << " over " << collapses
              << " chord collapses (" << rolled << " rolled back); four-sided " << quads << "/"
              << faces << " (" << (faces ? 100.0 * quads / faces : 0.0) << "%) -> " << sQuads << "/"
              << sFaces << " (" << (sFaces ? 100.0 * sQuads / sFaces : 0.0) << "%), " << tAfter
              << " T-junction(s) left\n"
              << blocks << " blocks; mean coverage " << (n ? 100.0 * coverSum / n : 0.0) << "%, "
              << full << "/" << n << " models covered completely\n";
    // A block in two materials is an element no analysis code can integrate.
    // Said here because the table cannot show it and the coverage does not
    // count it; with the interfaces honoured it should be zero.
    if (straddling > 0)
        std::cout << straddling << " block(s) on " << straddlingModels
                  << " model(s) sit in more than one material\n";
    return (sSound == n) ? 0 : 1;
}
