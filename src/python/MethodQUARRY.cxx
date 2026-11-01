// Mesh.quarry(cdf, dcdf, offset, offsetc): QUARRY (src/QUARRY/QUARRY.hxx), the quad blocks SQuadGen's chart
// distance fields describe, cut out of the mesh they were sampled on -- CrossGen's replacement for the release's
// QuadExtraction. The four fields are per triangle of the Mesh, at its centroid, in the Mesh's units (offsets as
// vectors from the centroid), which is what squadgen2d's round trip and SQ-Diffuse inference decode. Meshing is
// BlockQuadMesh on the blocks, and smoothing ZIPLINE's TMOP settings, as for the other methods whose blocks carry
// no interior geometry of their own.
//
// Not one of crossgen.methods: it is not a method of making a layout from a domain, it reads one off fields a
// network made, and ORACLE has nothing to choose it for.
#include "python/Methods.hxx"
#include "python/PyTables.hxx"
#include "python/PyUtil.hxx"

#include <stdexcept>

#include "QUARRY/QUARRY.hxx"
#include "mesh/BlockQuadMesh.hxx"

#define OPT(path) CG_OPTION(O, path)
#define REP(path) CG_REPORT(O, path)

namespace pycg {
namespace {

const Table<QUARRY::Options> &optionsTable() {
    using O = QUARRY::Options;
    static const Table<O> t{
        OPT(charts.noise),   OPT(charts.sameTarget), OPT(charts.core),    OPT(charts.minPeak),
        OPT(duals.noise),    OPT(duals.sameTarget),  OPT(duals.core),     OPT(duals.minPeak),
        OPT(collapseLength), OPT(cornerShare),       OPT(dualReach),      OPT(cornerGap), OPT(cornerSpacing), OPT(dcdfCorner), OPT(cornerSnapEdges), OPT(cornerSnapAngle),
        OPT(unresolvedEdges), OPT(repairRounds),     OPT(sidePasses), OPT(endClearance),
    };
    return t;
}

const Table<QUARRY::Report> &reportTable() {
    using O = QUARRY::Report;
    static const Table<O> t{
        REP(charts.peaks),      REP(charts.noiseMerges), REP(charts.targetMerges), REP(charts.dissolved),
        REP(charts.regions),    REP(charts.splitPieces), REP(charts.unreached),
        REP(duals.peaks),       REP(duals.noiseMerges),  REP(duals.targetMerges),  REP(duals.dissolved),
        REP(duals.regions),     REP(duals.splitPieces),  REP(duals.unreached),
        REP(junctions),         REP(arcs),               REP(collapsedArcs),       REP(addedCorners),
        REP(dissolvedCharts),   REP(snappedCorners),     REP(repairRounds),       REP(blocks),              REP(failedCharts),
        REP(unmatchedSides),    REP(macroVertices),      REP(macroEdges),          REP(boundaryVertices),
        REP(seconds),           REP(covers),             REP(messages),
    };
    return t;
}

// ZIPLINE's Stage 7: metric 7 at the element corners.
mesh::TMOP::Options smoothingOptions() {
    mesh::TMOP::Options t;
    t.metric = mesh::TMOP::ShapeSize007;
    t.quadrature = mesh::TMOP::Corners;
    t.maxSweeps = 1000;
    return t;
}

class QuarryMethod final : public Method {
public:
    QuarryMethod(const Mesh &model, QUARRY::Fields fields, const QUARRY::Options &opts)
        : model_(std::make_shared<const Mesh>(model)), quarry_(model_, std::move(fields), opts) {}

    void run() { quarry_.run(); }

    const char *name() const override { return "quarry"; }
    const BlockDecomposition &decomposition() const override { return quarry_.decomposition(); }
    double coverage() const override { return decompositionCoverage(quarry_.decomposition(), *model_); }

    // The report, and what a picture of the extraction needs: the chart (after Stage 5) and dual region of every
    // triangle, each region's consensus centre, and why each chart that made no block made none.
    PyObject *report() const override {
        PyObject *d = toDict({source(reportTable(), quarry_.getReport())});
        if (!d) return nullptr;
        std::vector<Point> chartCentres, dualCentres;
        std::vector<double> chartPeaks;
        for (const quarry::Basins::Region &g : quarry_.charts().regions()) {
            chartCentres.push_back(g.centre);
            chartPeaks.push_back(g.peakValue);
        }
        for (const quarry::Basins::Region &g : quarry_.duals().regions()) dualCentres.push_back(g.centre);
        std::vector<int> corners, block, junctions, added, dissolved;
        PyObject *problems = PyList_New(0);
        if (!problems) {
            Py_DECREF(d);
            return nullptr;
        }
        for (const QUARRY::ChartStatus &s : quarry_.chartStatus()) {
            corners.push_back(s.corners);
            block.push_back(s.block ? 1 : 0);
            junctions.push_back(s.junctions);
            added.push_back(s.added);
            dissolved.push_back(s.dissolved ? 1 : 0);
            PyObject *p = PyUnicode_FromString(s.problem.c_str());
            if (!p || PyList_Append(problems, p) < 0) {
                Py_XDECREF(p);
                Py_DECREF(problems);
                Py_DECREF(d);
                return nullptr;
            }
            Py_DECREF(p);
        }
        const std::pair<const char *, PyObject *> items[] = {
            {"chart_labels", indexArray(quarry_.chartLabels())},
            {"dual_labels", indexArray(quarry_.duals().labels())},
            {"chart_centres", pointArray(chartCentres)},
            {"chart_peaks", doubleArray(chartPeaks)},
            {"dual_centres", pointArray(dualCentres)},
            {"chart_corners", indexArray(corners)},
            {"chart_junctions", indexArray(junctions)},
            {"chart_added", indexArray(added)},
            {"chart_dissolved", indexArray(dissolved)},
            {"chart_block", indexArray(block)},
            {"chart_problems", problems},
        };
        int rc = 0;
        for (const auto &[key, value] : items) {
            if (!value || PyDict_SetItemString(d, key, value) < 0) rc = -1;
            Py_XDECREF(value);
        }
        if (rc < 0) {
            Py_DECREF(d);
            return nullptr;
        }
        return d;
    }
    PyObject *options() const override { return toDict({source(optionsTable(), quarry_.getOptions())}); }

    int mesh(double h, PyObject *kwargs, MeshOutput &out) const override {
        BlockQuadMesh::Options mo;
        mo.targetEdgeLength = h;
        mesh::QuadMesh::Options no;
        if (applyKwargs(kwargs, "mesh", {target(blockQuadMeshOptionsTable(), mo),
                                         target(nodeOptionsTable(), no)}) < 0)
            return -1;
        std::unique_ptr<BlockQuadMesh> grid;
        const bool ok = compute([&] {
            if (decomposition().blocks.empty()) throw std::runtime_error("quarry: no block to mesh");
            grid = std::make_unique<BlockQuadMesh>(decomposition(), mo);
            if (grid->getReport().quads <= 0) throw std::runtime_error("quarry: the mesh has no elements");
            out.mesh = std::make_unique<mesh::QuadMesh>(mesh::QuadMesh::from(*grid, no));
        });
        if (!ok) return -1;
        out.report = toDict({source(blockQuadMeshReportTable(), grid->getReport())});
        out.smoothing = smoothingOptions();
        return out.report ? 0 : -1;
    }

    double defaultTarget() const override { return BlockQuadMesh::Options().targetEdgeLength; }

private:
    std::shared_ptr<const Mesh> model_;
    QUARRY quarry_;
};

// A per-triangle field keyword, taken out of `kw` (a dict this function owns) so that what is left are options.
bool takeField(PyObject *kw, const char *key, int cols, size_t triangles, std::vector<double> &out) {
    PyObject *o = kw ? PyDict_GetItemString(kw, key) : nullptr;   // borrowed
    if (!o) {
        PyErr_Format(PyExc_TypeError, "quarry() needs the per-triangle field '%s'", key);
        return false;
    }
    Py_INCREF(o);
    PyDict_DelItemString(kw, key);
    PyObject *np = PyImport_ImportModule("numpy");
    PyObject *arr = np ? PyObject_CallMethod(np, "asarray", "Os", o, "float64") : nullptr;
    Py_XDECREF(np);
    Py_DECREF(o);
    if (!arr) return false;
    // cdf and dcdf as (F,) or (F, 1), as the decoder returns them.
    PyObject *shaped = cols == 0 ? PyObject_CallMethod(arr, "reshape", "(i)", -1) : (Py_INCREF(arr), arr);
    Py_DECREF(arr);
    if (!shaped) return false;
    Py_ssize_t rows = 0;
    const bool ok = readDoubles(shaped, cols, key, out, rows);
    Py_DECREF(shaped);
    if (!ok) return false;
    if (static_cast<size_t>(rows) != triangles) {
        PyErr_Format(PyExc_ValueError, "quarry(): '%s' has %zd rows for %zu triangles", key, rows, triangles);
        return false;
    }
    return true;
}

std::shared_ptr<Method> run(const Mesh &model, PyObject *kwargs) {
    PyObject *kw = kwargs ? PyDict_Copy(kwargs) : PyDict_New();
    if (!kw) return nullptr;
    const size_t nt = model.triangles.size();
    QUARRY::Fields f;
    std::vector<double> oc, ov;
    const bool ok = takeField(kw, "cdf", 0, nt, f.cdf) && takeField(kw, "dcdf", 0, nt, f.dcdf) &&
                    takeField(kw, "offsetc", 2, nt, oc) && takeField(kw, "offset", 2, nt, ov);
    QUARRY::Options opts;
    if (!ok || applyKwargs(kw, "quarry", {target(optionsTable(), opts)}) < 0) {
        Py_DECREF(kw);
        return nullptr;
    }
    Py_DECREF(kw);
    f.offsetCentre.resize(nt);
    f.offsetCorner.resize(nt);
    for (size_t t = 0; t < nt; ++t) {
        f.offsetCentre[t] = {oc[2 * t], oc[2 * t + 1]};
        f.offsetCorner[t] = {ov[2 * t], ov[2 * t + 1]};
    }
    std::shared_ptr<QuarryMethod> m;
    if (!compute([&] {
            m = std::make_shared<QuarryMethod>(model, std::move(f), opts);
            m->run();
        }))
        return nullptr;
    return m;
}

PyObject *defaults(const std::string &stage) {
    const QUARRY::Options o;
    if (stage == "method") return toDict({source(optionsTable(), o)});
    if (stage == "mesh") {
        const BlockQuadMesh::Options mo;
        const mesh::QuadMesh::Options nodes;
        return meshingDefaults(mo.targetEdgeLength,
                               {source(blockQuadMeshOptionsTable(), mo), source(nodeOptionsTable(), nodes)});
    }
    if (stage == "smooth") return smoothingDefaults(smoothingOptions());
    return unknownStage("quarry", stage);
}

void buildTables() {
    optionsTable();
    reportTable();
}

}  // namespace

const MethodSpec quarrySpec = {
    "quarry",
    "quarry(cdf, dcdf, offsetc, offset, **options) -> BlockDecomposition\n\n"
    "QUARRY: the quad blocks SQuadGen's chart distance fields describe, cut out of\n"
    "this mesh -- CrossGen's replacement for the release's QuadExtraction. The four\n"
    "fields are per triangle, sampled at its centroid, in the mesh's units: cdf and\n"
    "dcdf as (F,) or (F, 1), offsetc (to the chart centre) and offset (to the\n"
    "subchart's corner) as (F, 2) vectors from the centroid. Charts come from the\n"
    "CDF and offsetc, layout vertices from the DCDF and offset, by the same\n"
    "persistence clustering; a chart's corners are the vertex regions it overlaps,\n"
    "a side on dS is dS, and a side between two charts is their common boundary.\n"
    "A chart that is no quad is left out (b.report['chart_problems']), and the\n"
    "result then does not cover the mesh. Keywords are QUARRY::Options in\n"
    "snake_case (crossgen.options('quarry')). Meshing is BlockQuadMesh.",
    run,
    defaults,
    buildTables,
};

}  // namespace pycg
