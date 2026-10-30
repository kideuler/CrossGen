// Mesh.oracle(): ORACLE (src/ORACLE/ORACLE.hxx), the selector's choice among
// the five methods -- not a sixth method. It runs the ZIPLINE probe, asks
// py/selector.onnx for its ranking, runs at most two more of the five, and
// keeps one of their runs; the BlockDecomposition this returns is that run's
// blocks, and mesh() meshes them with that method's own mesher (the run's
// oracle::Candidate, kept alive inside this object, as the other bindings keep
// their pipelines).
//
// It is a Mesh method of its own rather than a sixth entry of crossgen.methods,
// because that tuple is the five layout methods -- what ORACLE chooses among --
// and a loop over it that compares methods (test_crossgen.py's, or a corpus
// sweep's) would otherwise run the kept method a second time under another name.
//
// ORACLE's runs are oracle::Candidate's, which repeat this module's
// Mesh.<method>() with no keywords (Candidate.hxx says how, and that they were
// checked equal), so m.oracle() returns the same blocks as m.<chosen>(). The
// selector is the one ORACLE reads when nothing names another:
// $CROSSGEN_SELECTOR, else the source tree's py/selector.onnx. A module built
// without ONNX Runtime has ORACLE but no selector, and m.oracle() raises.
#include "python/Methods.hxx"
#include "python/PyTables.hxx"
#include "python/PyUtil.hxx"

#include <stdexcept>

#include "ORACLE/ORACLE.hxx"
#include "ORACLE/Selector.hxx"

#define OPT(path) CG_OPTION(O, path)
#define REP(path) CG_REPORT(O, path)

namespace pycg {
namespace {

// ORACLE::Options less `selector`, a path, which $CROSSGEN_SELECTOR names
// instead (keywords carry numbers and flags only).
const Table<ORACLE::Options> &optionsTable() {
    using O = ORACLE::Options;
    static const Table<O> t{OPT(guardDomain)};
    return t;
}

// ORACLE::Report less `ranking` (indices into `methods`) and `runs` (a list),
// which report() writes out by hand: the ranking as method names, the runs as
// one dict each.
const Table<ORACLE::Report> &reportTable() {
    using O = ORACLE::Report;
    static const Table<O> t{
        REP(selector),   REP(described),   REP(metric),         REP(higherBetter),
        REP(probe),      REP(fallback),    REP(methods),        REP(inputNames),
        REP(inputs),     REP(utility),     REP(pValid),         REP(quality),
        REP(outOfDomain), REP(chosen),     REP(valid),          REP(decision),
        REP(droppedVertices), REP(seconds), REP(messages),
    };
    return t;
}

const Table<ORACLE::Run> &runTable() {
    using O = ORACLE::Run;
    static const Table<O> t{REP(method), REP(why),  REP(valid),  REP(raised), REP(error),
                            REP(coverage), REP(blocks), REP(metric), REP(seconds)};
    return t;
}

// Candidate::MeshSettings less `target` (BlockDecomposition.mesh's `h`).
const Table<oracle::Candidate::MeshSettings> &meshTable() {
    using O = oracle::Candidate::MeshSettings;
    static const Table<O> t{OPT(minIntervals), OPT(maxIntervals)};
    return t;
}

// Sets d[key] = v and drops the caller's reference to v; false with an
// exception set when v is null or the dict refuses it.
bool put(PyObject *d, const char *key, PyObject *v) {
    if (!v) return false;
    const int rc = PyDict_SetItemString(d, key, v);
    Py_DECREF(v);
    return rc == 0;
}

class OracleMethod final : public Method {
public:
    OracleMethod(const Mesh &model, const ORACLE::Options &opts)
        : oracle_(std::make_shared<const Mesh>(model), opts) {}

    void run() {
        oracle_.run();
        // run() keeps the best run even when none is valid (BlockDecomposition
        // .valid says so then); no choice at all means the selector could not
        // be read or asked, which is an error, not a layout.
        if (!oracle_.hasChoice()) {
            std::string why = "oracle: no run was kept";
            for (const std::string &m : oracle_.getReport().messages) why += "; " + m;
            throw std::runtime_error(why);
        }
    }

    const char *name() const override { return "oracle"; }
    const BlockDecomposition &decomposition() const override { return oracle_.getDecomposition(); }
    double coverage() const override { return oracle_.getChoice().coverage(); }

    // ORACLE's report, plus `method` (the method whose run was kept),
    // `ranking` (the methods by the selector's utility, best first) and `runs`.
    PyObject *report() const override {
        const ORACLE::Report &r = oracle_.getReport();
        PyObject *d = toDict({source(reportTable(), r)});
        if (!d) return nullptr;
        std::vector<std::string> ranking;
        for (const int i : r.ranking)
            if (i >= 0 && i < static_cast<int>(r.methods.size())) ranking.push_back(r.methods[i]);
        PyObject *runs = PyList_New(0);
        bool ok = runs != nullptr;
        for (size_t i = 0; ok && i < r.runs.size(); ++i) {
            PyObject *one = toDict({source(runTable(), r.runs[i])});
            ok = one && PyList_Append(runs, one) == 0;
            Py_XDECREF(one);
        }
        if (!ok) Py_XDECREF(runs);
        // put() takes the reference it is given, whether or not it succeeds.
        ok = ok && put(d, "runs", runs) && put(d, "method", toPy(std::string(oracle_.getChoice().name()))) &&
             put(d, "ranking", toPy(ranking));
        if (!ok) {
            Py_DECREF(d);
            return nullptr;
        }
        return d;
    }
    PyObject *options() const override {
        PyObject *d = toDict({source(optionsTable(), oracle_.getOptions())});
        if (d && !put(d, "selector", toPy(oracle_.getReport().selector))) Py_CLEAR(d);
        return d;
    }

    int mesh(double h, PyObject *kwargs, MeshOutput &out) const override {
        oracle::Candidate::MeshSettings s;
        s.target = h;
        mesh::QuadMesh::Options no;
        if (applyKwargs(kwargs, "mesh", {target(meshTable(), s), target(nodeOptionsTable(), no)}) < 0)
            return -1;
        const oracle::Candidate &c = oracle_.getChoice();
        if (!compute([&] {
                std::unique_ptr<oracle::CandidateMesh> cm = c.mesh(s);
                out.mesh = std::make_unique<mesh::QuadMesh>(cm->quadMesh(no));
            }))
            return -1;
        out.smoothing = c.smoothing();
        out.pillow = c.pillows();
        out.report = PyDict_New();
        if (!out.report || !put(out.report, "method", toPy(std::string(c.name())))) return -1;
        return 0;
    }

    double defaultTarget() const override { return oracle::Candidate::MeshSettings().target; }

private:
    ORACLE oracle_;
};

std::shared_ptr<Method> run(const Mesh &model, PyObject *kwargs) {
    ORACLE::Options opts;
    if (applyKwargs(kwargs, "oracle", {target(optionsTable(), opts)}) < 0) return nullptr;
    std::shared_ptr<OracleMethod> m;
    if (!compute([&] {
            m = std::make_shared<OracleMethod>(model, opts);
            m->run();
        }))
        return nullptr;
    return m;
}

PyObject *defaults(const std::string &stage) {
    if (stage == "method") {
        PyObject *d = toDict({source(optionsTable(), ORACLE::Options())});
        if (d && !put(d, "selector", toPy(oracle::Selector::defaultPath()))) Py_CLEAR(d);
        return d;
    }
    if (stage == "mesh") {
        const oracle::Candidate::MeshSettings s;
        const mesh::QuadMesh::Options nodes;
        return meshingDefaults(s.target, {source(meshTable(), s), source(nodeOptionsTable(), nodes)});
    }
    // All five smooth with the block-grid settings (metric 7, corners, 1000
    // sweeps); whether the first smooth() pillows depends on the method kept,
    // MERIDIAN and TORSION, so it is reported here as the default, off.
    if (stage == "smooth") return smoothingDefaults(blockGridSmoothing());
    return unknownStage("oracle", stage);
}

void buildTables() {
    optionsTable();
    reportTable();
    runTable();
    meshTable();
}

}  // namespace

const MethodSpec oracleSpec = {
    "oracle",
    "oracle(**options) -> BlockDecomposition\n\n"
    "ORACLE (src/ORACLE): the selector's choice among crossgen.methods, not a\n"
    "method of its own. It runs the ZIPLINE probe, asks the selector\n"
    "($CROSSGEN_SELECTOR, else py/selector.onnx) for its ranking, runs at most\n"
    "two more methods and keeps one run: the same blocks as m.<method>() for the\n"
    "method kept, which report['method'] names (report['decision'] says why).\n"
    "mesh() uses that method's mesher. Keyword: guard_domain=True (run the probe\n"
    "and the fallback where the selector's prediction is out of its range).\n"
    "Raises RuntimeError when the selector cannot be read.",
    run,
    defaults,
    buildTables,
};

}  // namespace pycg
