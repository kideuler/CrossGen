// Utility to run ORACLE (src/ORACLE) on one model or many: the model's
// boundary features, the ZIPLINE probe, the selector's ranking of the five
// methods (py/selector.onnx), and the method the probe rule then keeps -- the
// block decomposition a caller of ORACLE gets, and optionally a mesh on it.
//
//   TestORACLE [options] <mesh.obj> [more.obj ...]
//
// Per model it prints the selector's ranking (expected utility, probability
// of a valid run, predicted metric), every method it ran and why, and the
// decomposition it kept: its method, its block count, its coverage and the
// selector's metric on it.
//
// Options:
//   -v                  the inputs handed to the selector, one per line
//   --selector <f>      the selector to read (default: $CROSSGEN_SELECTOR, else
//                       py/selector.onnx of the source tree)
//   --follow-ranking    follow the selector's ranking even where its answer is
//                       no prediction (ORACLE::Options::guardDomain off)
//   --mesh <h>          mesh the kept blocks with their method's own mesher
//   --tmop <sweeps>     then smooth that mesh with the method's TMOP settings,
//                       pillowing first where the method does (implies --mesh
//                       0.05 if unset)
//   --blocks-obj <f>    the kept decomposition's macro edges as OBJ polylines
//   --obj <f>           the mesh, smoothed if --tmop (one model only)
//
// The exit code is 0 when every model kept a valid decomposition and 1 when
// one did not; 2 on bad arguments; and 77, ctest's skip, when the selector
// cannot be read at all -- a build without ONNX Runtime, or no selector file.
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "ORACLE/ORACLE.hxx"
#include "ORACLE/Selector.hxx"
#include "mesh/Pillow.hxx"
#include "mesh/QuadMesh.hxx"
#include "mesh/TMOP.hxx"

namespace {

struct Flags {
    bool verbose = false;
    ORACLE::Options oracle;
    double meshTarget = 0.0;   // 0: no mesh
    int tmopSweeps = 0;        // 0: no smoothing
    std::string blocksObj, obj;
};

std::string fixed(double v, int digits = 3) {
    if (std::isnan(v)) return "--";
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(digits) << v;
    return oss.str();
}

void usage() {
    std::cerr << "usage: TestORACLE [-v] [--selector f] [--follow-ranking] [--mesh h] [--tmop sweeps] "
                 "[--blocks-obj f] [--obj f] <mesh.obj> [more.obj ...]\n";
}

// The kept blocks meshed, and smoothed, as crossgen's mesh() and smooth()
// would: the method's own mesher and its own TMOP settings. False when the
// mesh could not be built.
bool meshChoice(const ORACLE &o, const Flags &f) {
    const oracle::Candidate &c = o.getChoice();
    oracle::Candidate::MeshSettings s;
    s.target = f.meshTarget;
    std::unique_ptr<oracle::CandidateMesh> cm;
    mesh::QuadMesh qm;
    try {
        cm = c.mesh(s);
        qm = cm->quadMesh(mesh::QuadMesh::Options());
    } catch (const std::exception &e) {
        std::cout << "  mesh: FAILED: " << e.what() << "\n";
        return false;
    }
    const mesh::QuadMesh::Quality &q = qm.quality;
    std::cout << "  mesh (" << c.name() << "'s mesher, h " << f.meshTarget << "): " << q.quads
              << " quads, scaled Jacobian " << fixed(q.minScaledJacobian) << " worst, "
              << fixed(q.meanScaledJacobian) << " mean";
    if (q.invertedQuads > 0) std::cout << ", " << q.invertedQuads << " inverted";
    std::cout << "\n";
    if (f.tmopSweeps > 0) {
        if (c.pillows()) {
            mesh::Pillow pillow(qm);
            pillow.run();
            const mesh::Pillow::Report &pr = pillow.getReport();
            if (pr.ran)
                std::cout << "  pillow: " << pr.defectsBefore << " flat feature corner(s), "
                          << pr.quadsBefore << " -> " << pr.quadsAfter << " quads\n";
        }
        mesh::TMOP::Options t = c.smoothing();
        t.maxSweeps = f.tmopSweeps;
        mesh::TMOP smoother(qm, t);
        smoother.run();
        const mesh::TMOP::Report &r = smoother.getReport();
        std::cout << "  TMOP: " << r.sweeps << " sweep(s), scaled Jacobian "
                  << fixed(r.minScaledJacobianBefore) << " -> " << fixed(r.minScaledJacobianAfter)
                  << " worst, " << fixed(r.meanScaledJacobianBefore) << " -> "
                  << fixed(r.meanScaledJacobianAfter) << " mean\n";
    }
    if (!f.obj.empty()) {
        if (qm.writeOBJ(f.obj)) std::cout << "  mesh -> " << f.obj << "\n";
        else std::cout << "  could not write " << f.obj << "\n";
    }
    return true;
}

}  // namespace

int main(int argc, char **argv) {
    Flags f;
    std::vector<std::string> models;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto value = [&](const char *flag) -> const char * {
            if (i + 1 >= argc) {
                std::cerr << flag << " needs a value\n";
                std::exit(2);
            }
            return argv[++i];
        };
        if (a == "-v") f.verbose = true;
        else if (a == "--selector") f.oracle.selector = value("--selector");
        else if (a == "--follow-ranking") f.oracle.guardDomain = false;
        else if (a == "--mesh") f.meshTarget = std::atof(value("--mesh"));
        else if (a == "--tmop") f.tmopSweeps = std::atoi(value("--tmop"));
        else if (a == "--blocks-obj") f.blocksObj = value("--blocks-obj");
        else if (a == "--obj") f.obj = value("--obj");
        else if (a == "-h" || a == "--help") {
            usage();
            return 0;
        } else if (!a.empty() && a[0] == '-') {
            std::cerr << "unknown option " << a << "\n";
            usage();
            return 2;
        } else {
            models.push_back(a);
        }
    }
    if (models.empty()) {
        usage();
        return 2;
    }
    if (f.tmopSweeps > 0 && !(f.meshTarget > 0.0)) f.meshTarget = 0.05;
    if (models.size() > 1 && (!f.blocksObj.empty() || !f.obj.empty())) {
        std::cerr << "--blocks-obj and --obj take one model\n";
        return 2;
    }

    // Read the selector once up front, so that a build or a checkout that
    // cannot is one clear line and a skip, not the same failure per model.
    const std::string selectorPath =
        f.oracle.selector.empty() ? oracle::Selector::defaultPath() : f.oracle.selector;
    try {
        const oracle::Selector s(selectorPath);
        const oracle::Selector::Description &d = s.description();
        std::cout << "selector " << selectorPath << ": " << d.inputs.size() << " inputs -> ";
        for (size_t i = 0; i < d.methods.size(); ++i) std::cout << (i ? ", " : "") << d.methods[i];
        std::cout << " on " << d.metric << (d.higherBetter ? "" : " (lower better)") << "; probe "
                  << (d.probe.empty() ? "none" : d.probe) << ", fallback "
                  << (d.fallback.empty() ? "none" : d.fallback)
                  << (d.fromFile ? "" : " (the file names no columns: the trainer's defaults assumed)") << "\n";
    } catch (const std::exception &e) {
        std::cerr << "TestORACLE: " << e.what() << "\n";
        return 77;
    }

    int kept = 0, failed = 0;
    for (const std::string &path : models) {
        std::shared_ptr<const Mesh> mesh;
        try {
            mesh = std::make_shared<const Mesh>(path);
        } catch (const std::exception &e) {
            std::cout << path << ": cannot load: " << e.what() << "\n";
            ++failed;
            continue;
        }
        ORACLE o(mesh, f.oracle);
        const bool ok = o.run();
        const ORACLE::Report &r = o.getReport();
        std::cout << path << ":";
        if (r.utility.empty()) {
            std::cout << " no ranking";
            for (const std::string &m : r.messages) std::cout << " -- " << m;
            std::cout << "\n";
            ++failed;
            continue;
        }
        std::cout << "\n";
        if (r.droppedVertices > 0)
            std::cout << "  " << r.droppedVertices << " vertices in no triangle dropped first\n";
        if (f.verbose) {
            for (size_t i = 0; i < r.inputNames.size(); ++i)
                std::cout << "    " << std::left << std::setw(34) << r.inputNames[i] << std::right
                          << std::setprecision(10) << r.inputs[i] << "\n";
        }
        std::cout << "  ranking:";
        for (size_t k = 0; k < r.ranking.size(); ++k) {
            const int i = r.ranking[k];
            std::cout << (k ? " >" : "") << " " << r.methods[i] << " " << fixed(r.utility[i]) << " ("
                      << fixed(100.0 * r.pValid[i], 0) << "% valid, " << r.metric << " "
                      << fixed(r.quality[i]) << ")";
        }
        std::cout << "\n";
        if (r.outOfDomain)
            std::cout << "  the selector is outside its training data on this model: a predicted " << r.metric
                      << " outside [0, 1]" << (f.oracle.guardDomain ? "; ranking not followed" : "; followed anyway")
                      << "\n";
        for (const ORACLE::Run &run : r.runs) {
            std::cout << "  ran " << std::left << std::setw(9) << run.method << std::right << " ("
                      << run.why << "): ";
            if (run.raised) std::cout << "FAILED: " << run.error;
            else
                std::cout << (run.valid ? "valid" : "not valid") << ", " << run.blocks << " block(s), "
                          << fixed(100.0 * run.coverage, 1) << "% covered, " << r.metric << " "
                          << fixed(run.metric);
            std::cout << ", " << fixed(run.seconds, 2) << " s\n";
        }
        std::cout << "  " << r.decision << "\n";
        if (r.chosen >= 0) {
            const ORACLE::Run &c = r.runs[r.chosen];
            std::cout << "  ORACLE: " << c.method << ", " << c.blocks << " block(s), " << r.metric << " "
                      << fixed(c.metric) << (ok ? " [valid]" : " [NOT valid]") << ", "
                      << fixed(r.seconds, 2) << " s in all\n";
        }
        for (const std::string &m : r.messages) std::cout << "  " << m << "\n";
        if (ok) ++kept;
        else ++failed;

        if (o.hasChoice() && !f.blocksObj.empty()) {
            if (o.getDecomposition().writeEdgesOBJ(f.blocksObj))
                std::cout << "  blocks -> " << f.blocksObj << "\n";
            else
                std::cout << "  could not write " << f.blocksObj << "\n";
        }
        if (o.hasChoice() && f.meshTarget > 0.0 && !meshChoice(o, f) && ok) {
            ++failed;
            --kept;
        }
    }
    std::cout << kept << " of " << models.size() << " model(s) kept a valid decomposition\n";
    return failed == 0 ? 0 : 1;
}
