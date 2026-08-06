// Load a mesh, compute the MBO cross field, trace its separatrices and build
// the quad layout with T-junctions they cut the model into -- steps 1 and 2 of
// Algorithm 1 in Viertel, Osting and Staten, IMR 2019.
//
//   TraceMesh <mesh.obj> [more.obj ...]
//
// Per mesh it prints a line of counts and, unless --no-vtu is given, writes
// <stem>_layout_arcs.vtu, <stem>_layout_faces.vtu and <stem>_layout_nodes.vtu
// next to wherever it is run.
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "crossfield/CrossField.hxx"
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
    }
    return "?";
}

static std::string stemOf(const std::string &path) {
    size_t slash = path.find_last_of("/\\");
    std::string base = (slash == std::string::npos) ? path : path.substr(slash + 1);
    size_t dot = base.find_last_of('.');
    return (dot == std::string::npos) ? base : base.substr(0, dot);
}

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
    int singularities = 0;
    int separatrices = 0;
    int stuck = 0;
    int dangling = 0;
    int faces = 0;
    int quads = 0;
    int tJunctions = 0;
    int skew = 0;
    double areaRatio = 0.0;
    // after Sec. 4
    int sFaces = 0, sQuads = 0, sTJunctions = 0, sCollapses = 0, sRolledBack = 0;
    double sAreaRatio = 0.0;
    bool sSound = false;
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

static bool gSimplify = true;
static double gCutRadius = SeparatrixTrace::Settings().singularityCutRadius;
static double gTangential = SeparatrixTrace::Settings().tangentialAngle;
static bool gEvenRays = SeparatrixTrace::Settings().evenCornerRays;

static Outcome processMesh(const std::string &path, bool writeVtu, bool verbose) {
    Outcome out;
    out.name = stemOf(path);

    std::shared_ptr<Mesh> mesh;
    try {
        mesh = std::make_shared<Mesh>(path);
    } catch (const std::exception &e) {
        std::cerr << "  failed to load " << path << ": " << e.what() << "\n";
        return out;
    }

    auto crossField = std::make_shared<CrossField>(mesh);
    crossField->initialize(1, 12345);
    const double nv = static_cast<double>(mesh->vertices.size());
    for (int i = 0; i < MBO_MAX_STEPS; ++i) {
        crossField->step();
        if (crossField->error < 2.0 * nv * 1e-9) break;
    }
    crossField->computeSingularities();

    SeparatrixTrace::Settings settings;
    settings.singularityCutRadius = gCutRadius;
    settings.tangentialAngle = gTangential;
    settings.evenCornerRays = gEvenRays;
    SeparatrixTrace trace(crossField, true, settings);
    trace.run();

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
    out.skew = rep.skewNodes;
    out.areaRatio = (rep.meshArea > 0.0) ? rep.totalArea / rep.meshArea : 0.0;
    out.sound = (rep.danglingEnds == 0) && (out.stuck == 0) && (rep.arcCrossings == 0) &&
                (std::fabs(out.areaRatio - 1.0) < 1e-9) && rep.faces > 0;
    out.allQuads = out.sound && rep.badFaces == 0;

    PartitionSimplify simp(layout);
    if (gSimplify) simp.run();
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

    if (verbose) {
        std::cout << "  mesh          " << mesh->triangles.size() << " triangles, "
                  << mesh->vertices.size() << " vertices\n";
        std::cout << "  cross field   " << trace.singularities.size() << " singularities, error "
                  << crossField->error << "\n";
        int byReason[7] = {0, 0, 0, 0, 0, 0, 0};
        for (const auto &s : trace.separatrices) byReason[static_cast<int>(s.termination_reason)]++;
        std::cout << "  separatrices  " << trace.separatrices.size() << " (";
        bool first = true;
        for (int i = 0; i < 7; ++i) {
            if (!byReason[i]) continue;
            if (!first) std::cout << ", ";
            std::cout << byReason[i] << " " << reasonName(static_cast<TerminationReason>(i));
            first = false;
        }
        std::cout << ")\n";
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
                static const char *kindName[] = {"sing", "corn", "bhit", "xing", "hetc", "tjct", "loose"};
                for (size_t k = 0; k < f.nodes.size(); ++k)
                    std::cout << " " << (f.isCorner[k] ? "*" : "-")
                              << kindName[static_cast<int>(layout.getNodes()[f.nodes[k]].kind)]
                              << (int)std::lround(f.turn[k] * 180.0 / M_PI)
                              << "/" << layout.getNodes()[f.nodes[k]].darts.size();
                std::cout << "\n";
                ++shown;
            }
        }
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

int main(int argc, char **argv) {
    std::vector<std::string> meshes;
    bool writeVtu = true;
    bool verbose = false;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--no-vtu") == 0) writeVtu = false;
        else if (std::strcmp(argv[i], "--cut-radius") == 0 && i + 1 < argc) gCutRadius = std::atof(argv[++i]);
        else if (std::strcmp(argv[i], "--tangential") == 0 && i + 1 < argc) gTangential = std::atof(argv[++i]) * M_PI / 180.0;
        else if (std::strcmp(argv[i], "--square-rays") == 0) gEvenRays = false;
        else if (std::strcmp(argv[i], "--no-simplify") == 0) gSimplify = false;
        else if (std::strcmp(argv[i], "-v") == 0) verbose = true;
        else meshes.push_back(argv[i]);
    }
    if (meshes.empty()) {
        std::cerr << "Usage: " << argv[0] << " [-v] [--no-vtu] <mesh.obj> [more.obj ...]\n";
        return 1;
    }

    std::vector<Outcome> all;
    for (const auto &m : meshes) {
        if (meshes.size() > 1 || verbose) std::cout << stemOf(m) << "\n";
        all.push_back(processMesh(m, writeVtu, verbose || meshes.size() == 1));
    }

    std::cout << "\n"
              << std::left << std::setw(12) << "mesh" << std::right << std::setw(6) << "sing"
              << std::setw(6) << "seps" << std::setw(7) << "faces" << std::setw(7) << "quads"
              << std::setw(7) << "T-jct" << std::setw(7) << "loose" << std::setw(7) << "skew"
              << "  |" << std::setw(6) << "coll" << std::setw(7) << "faces" << std::setw(7)
              << "quads" << std::setw(7) << "T-jct" << std::setw(6) << "back" << "  sound\n";
    int sound = 0, allQuads = 0, faces = 0, quads = 0;
    for (const auto &o : all) {
        std::cout << std::left << std::setw(12) << o.name << std::right << std::setw(6)
                  << o.singularities << std::setw(6) << o.separatrices << std::setw(7) << o.faces
                  << std::setw(7) << o.quads << std::setw(7) << o.tJunctions << std::setw(7)
                  << (o.stuck + o.dangling) << std::setw(7) << o.skew << "  |" << std::setw(6)
                  << o.sCollapses << std::setw(7) << o.sFaces << std::setw(7) << o.sQuads
                  << std::setw(7) << o.sTJunctions << std::setw(6) << o.sRolledBack
                  << ((o.sound && o.sSound) ? "    yes" : "     NO") << "\n";
        if (o.sound) ++sound;
        if (o.allQuads) ++allQuads;
        faces += o.faces;
        quads += o.quads;
    }
    int sSound = 0, sFaces = 0, sQuads = 0, collapses = 0, rolled = 0;
    for (const auto &o : all) {
        if (o.sound && o.sSound) ++sSound;
        sFaces += o.sFaces;
        sQuads += o.sQuads;
        collapses += o.sCollapses;
        rolled += o.sRolledBack;
    }
    std::cout << std::fixed << std::setprecision(1) << sound << "/" << all.size() << " layouts sound before "
              << "simplification, " << sSound << "/" << all.size() << " after; " << allQuads << "/"
              << all.size() << " all four-sided before\n"
              << "components " << faces << " -> " << sFaces << " over " << collapses
              << " chord collapses (" << rolled << " rolled back); four-sided " << quads << "/"
              << faces << " (" << (faces ? 100.0 * quads / faces : 0.0) << "%) -> " << sQuads << "/"
              << sFaces << " (" << (sFaces ? 100.0 * sQuads / sFaces : 0.0) << "%)\n";
    return (sound == static_cast<int>(all.size()) && sSound == static_cast<int>(all.size())) ? 0 : 1;
}
