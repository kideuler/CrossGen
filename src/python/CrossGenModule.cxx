// crossgen: the Python extension module over CrossGen.
//
//     import crossgen
//     m = crossgen.load("data/meshes/singlemat/geom012.obj")
//     b = m.meridian()                 # or zipline / umber / torsion / atlas
//     q = b.mesh(0.05)                 # quad mesh at a target edge length
//     q.smooth(1000)                   # TMOP, in place
//     e = m.shape_dna(count=50)        # numpy array of normalised eigenvalues
//     f = m.boundary_features()        # corners, holes, T3's bound, as a dict
//
// Three types, each a thin handle on C++ objects the rest of the codebase
// already has:
//
//   Mesh                a ::Mesh (mesh/Mesh.hxx): a triangle mesh with a material
//                       id per triangle, as every method takes it.
//   BlockDecomposition  one method's result: the shared BlockDecomposition
//                       (mesh/BlockDecomposition.hxx) together with the pipeline
//                       that produced it, which is kept because each method
//                       meshes its blocks through its own stage (Methods.hxx).
//   QuadMesh            a mesh::QuadMesh (mesh/QuadMesh.hxx), what TMOP smooths,
//                       with the TMOP settings its method smooths with.
//
// Keyword options are the C++ Options structs' fields in snake_case, with their
// C++ defaults (PyOptions.hxx); crossgen.options() lists them per method and
// stage. Arrays come back as numpy arrays, and every C++ call runs with the GIL
// released, one at a time (PyUtil.hxx).
#include "python/Methods.hxx"
#include "python/PyOptions.hxx"
#include "python/PyTables.hxx"
#include "python/PyUtil.hxx"

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <limits>
#include <new>
#include <set>
#include <stdexcept>
#include <string>

#include "ShapeDNA/ShapeDNA.hxx"
#include "mesh/BoundaryFeatures.hxx"
#include "mesh/Mesh.hxx"

namespace pycg {
namespace {

// A METH_VARARGS | METH_KEYWORDS function as the PyCFunction a PyMethodDef
// holds; through void(*)(void), which is what silences the cast warning.
template <class F>
PyCFunction asCFunction(F f) {
    return reinterpret_cast<PyCFunction>(reinterpret_cast<void (*)(void)>(f));
}

// ---------------------------------------------------------------------------
// Object layouts. The C++ members are constructed with placement new after
// tp_alloc and destroyed by hand in tp_dealloc.
// ---------------------------------------------------------------------------
struct MeshObject {
    PyObject_HEAD
    std::shared_ptr<const Mesh> mesh;
    // BoundaryFeatures at the default options, made on first use: what every
    // BlockDecomposition of this mesh grades its regularity against.
    std::shared_ptr<const BoundaryFeatures> features;
};

struct BlocksObject {
    PyObject_HEAD
    std::shared_ptr<const Method> method;
    // The model's, from the Mesh it was made from (regularity() reads them).
    std::shared_ptr<const BoundaryFeatures> features;
};

struct QuadMeshObject {
    PyObject_HEAD
    std::unique_ptr<mesh::QuadMesh> mesh;
    PyObject *report;                  // the mesher's, a dict
    mesh::TMOP::Options smoothing;     // the method's defaults
    std::string method;
    // Set while smooth() runs with the GIL released. Everything that reads
    // the mesh refuses while it is, rather than read nodes mid-move.
    bool busy;
};

PyTypeObject MeshType = {PyVarObject_HEAD_INIT(nullptr, 0)};
PyTypeObject BlocksType = {PyVarObject_HEAD_INIT(nullptr, 0)};
PyTypeObject QuadMeshType = {PyVarObject_HEAD_INIT(nullptr, 0)};

MeshObject *asMesh(PyObject *o) { return reinterpret_cast<MeshObject *>(o); }
BlocksObject *asBlocks(PyObject *o) { return reinterpret_cast<BlocksObject *>(o); }
QuadMeshObject *asQuad(PyObject *o) { return reinterpret_cast<QuadMeshObject *>(o); }

// A filesystem path argument (str or os.PathLike) as a std::string.
bool pathArg(PyObject *args, const char *fn, std::string &out) {
    PyObject *bytes = nullptr;
    const std::string fmt = std::string("O&:") + fn;
    if (!PyArg_ParseTuple(args, fmt.c_str(), PyUnicode_FSConverter, &bytes)) return false;
    out = PyBytes_AS_STRING(bytes);
    Py_DECREF(bytes);
    return true;
}

PyObject *writeResult(bool ok, const std::string &path) {
    if (!ok) {
        PyErr_Format(PyExc_OSError, "could not write '%s'", path.c_str());
        return nullptr;
    }
    Py_RETURN_NONE;
}

// `kwargs` without `name`, and the value it had (borrowed, or null). The copy
// is a new reference, or null when there are no keywords left; `ok` is false
// with an exception set on failure.
PyObject *popKeyword(PyObject *kwargs, const char *name, PyObject **value, bool &ok) {
    ok = true;
    *value = nullptr;
    if (!kwargs || PyDict_GET_SIZE(kwargs) == 0) return nullptr;
    PyObject *rest = PyDict_Copy(kwargs);
    if (!rest) { ok = false; return nullptr; }
    *value = PyDict_GetItemString(kwargs, name);
    if (*value && PyDict_DelItemString(rest, name) < 0) {
        Py_DECREF(rest);
        ok = false;
        return nullptr;
    }
    return rest;
}

// ===========================================================================
// Mesh
// ===========================================================================
PyObject *newMesh(std::shared_ptr<const Mesh> m) {
    PyObject *o = MeshType.tp_alloc(&MeshType, 0);
    if (!o) return nullptr;
    new (&asMesh(o)->mesh) std::shared_ptr<const Mesh>(std::move(m));
    new (&asMesh(o)->features) std::shared_ptr<const BoundaryFeatures>();
    return o;
}

void Mesh_dealloc(PyObject *self) {
    asMesh(self)->features.~shared_ptr();
    asMesh(self)->mesh.~shared_ptr();
    Py_TYPE(self)->tp_free(self);
}

// The mesh's BoundaryFeatures at the default options, made once. Null with an
// exception set if they could not be.
std::shared_ptr<const BoundaryFeatures> defaultFeatures(PyObject *self) {
    MeshObject *m = asMesh(self);
    if (!m->features) {
        const std::shared_ptr<const Mesh> model = m->mesh;
        std::shared_ptr<const BoundaryFeatures> f;
        if (!compute([&] { f = std::make_shared<const BoundaryFeatures>(*model); })) return nullptr;
        m->features = std::move(f);
    }
    return m->features;
}

// Mesh(vertices, triangles, materials=None)
PyObject *Mesh_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
    static const char *kw[] = {"vertices", "triangles", "materials", nullptr};
    PyObject *V = nullptr, *T = nullptr, *M = Py_None;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|O:Mesh", const_cast<char **>(kw), &V, &T, &M))
        return nullptr;

    std::vector<double> xy;
    std::vector<long long> tri, mat;
    Py_ssize_t nv = 0, nt = 0, nm = 0;
    if (!readDoubles(V, 2, "vertices", xy, nv) || !readIndices(T, 3, "triangles", tri, nt))
        return nullptr;
    if (M != Py_None) {
        if (!readIndices(M, 0, "materials", mat, nm)) return nullptr;
        if (nm != nt) {
            PyErr_Format(PyExc_ValueError, "materials has %zd entries for %zd triangles", nm, nt);
            return nullptr;
        }
    }
    for (const long long i : tri) {
        if (i < 0 || i >= nv) {
            PyErr_Format(PyExc_ValueError, "triangle vertex index %lld is not in [0, %zd)", i, nv);
            return nullptr;
        }
    }

    std::vector<Point> points(static_cast<size_t>(nv));
    for (Py_ssize_t i = 0; i < nv; ++i) points[i] = {xy[2 * i], xy[2 * i + 1]};
    std::vector<Triangle> triangles(static_cast<size_t>(nt));
    for (Py_ssize_t t = 0; t < nt; ++t) {
        Triangle abc = {static_cast<int>(tri[3 * t]), static_cast<int>(tri[3 * t + 1]),
                        static_cast<int>(tri[3 * t + 2])};
        // Counter-clockwise, as the .obj reader leaves every triangle.
        if (cross2(points[abc[1]] - points[abc[0]], points[abc[2]] - points[abc[0]]) < 0.0)
            std::swap(abc[1], abc[2]);
        triangles[t] = abc;
    }
    std::vector<int> matIds(mat.begin(), mat.end());

    std::shared_ptr<const Mesh> m;
    if (!compute([&] { m = std::make_shared<const Mesh>(points, triangles, matIds); })) return nullptr;
    PyObject *o = type->tp_alloc(type, 0);
    if (!o) return nullptr;
    new (&asMesh(o)->mesh) std::shared_ptr<const Mesh>(std::move(m));
    new (&asMesh(o)->features) std::shared_ptr<const BoundaryFeatures>();
    return o;
}

PyObject *Mesh_repr(PyObject *self) {
    const Mesh &m = *asMesh(self)->mesh;
    const std::set<int> materials(m.triangleMatId.begin(), m.triangleMatId.end());
    return PyUnicode_FromFormat("<crossgen.Mesh: %zu vertices, %zu triangles, %zu material(s)>",
                                m.vertices.size(), m.triangles.size(), materials.size());
}

PyObject *Mesh_vertices(PyObject *self, void *) { return pointArray(asMesh(self)->mesh->vertices); }
PyObject *Mesh_triangles(PyObject *self, void *) { return indexArray(asMesh(self)->mesh->triangles); }
PyObject *Mesh_materials(PyObject *self, void *) {
    return indexArray(asMesh(self)->mesh->triangleMatId);
}
PyObject *Mesh_boundaryEdges(PyObject *self, void *) {
    const Mesh &m = *asMesh(self)->mesh;
    std::vector<std::array<int, 2>> e;
    e.reserve(m.boundaryEdges.size());
    for (const int i : m.boundaryEdges) e.push_back(m.edges[i]);
    return indexArray(e);
}
PyObject *Mesh_numVertices(PyObject *self, void *) {
    return PyLong_FromSize_t(asMesh(self)->mesh->vertices.size());
}
PyObject *Mesh_numTriangles(PyObject *self, void *) {
    return PyLong_FromSize_t(asMesh(self)->mesh->triangles.size());
}

PyGetSetDef meshGetSet[] = {
    {"vertices", Mesh_vertices, nullptr, "(n, 2) float64 array of vertex positions.", nullptr},
    {"triangles", Mesh_triangles, nullptr,
     "(m, 3) int64 array of vertex indices, counter-clockwise.", nullptr},
    {"materials", Mesh_materials, nullptr, "(m,) int64 array: the material id of each triangle.",
     nullptr},
    {"boundary_edges", Mesh_boundaryEdges, nullptr,
     "(k, 2) int64 array: the vertex pairs of the edges on the boundary.", nullptr},
    {"num_vertices", Mesh_numVertices, nullptr, "Number of vertices.", nullptr},
    {"num_triangles", Mesh_numTriangles, nullptr, "Number of triangles.", nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr},
};

PyObject *newBlocks(std::shared_ptr<const Method> method, std::shared_ptr<const BoundaryFeatures> features) {
    PyObject *o = BlocksType.tp_alloc(&BlocksType, 0);
    if (!o) return nullptr;
    new (&asBlocks(o)->method) std::shared_ptr<const Method>(std::move(method));
    new (&asBlocks(o)->features) std::shared_ptr<const BoundaryFeatures>(std::move(features));
    return o;
}

// Mesh.zipline(**options), Mesh.umber(...), ...: one instantiation per method.
template <const MethodSpec *S>
PyObject *Mesh_runMethod(PyObject *self, PyObject *args, PyObject *kwargs) {
    if (PyTuple_GET_SIZE(args) != 0) {
        PyErr_Format(PyExc_TypeError, "%s() takes keyword arguments only", S->name);
        return nullptr;
    }
    // Held for the call, so the mesh outlives the run whatever else lets go of it.
    const std::shared_ptr<const Mesh> model = asMesh(self)->mesh;
    std::shared_ptr<const BoundaryFeatures> features = defaultFeatures(self);
    if (!features) return nullptr;
    std::shared_ptr<Method> m = S->run(*model, kwargs);
    if (!m) return nullptr;
    PyObject *b = newBlocks(m, std::move(features));
    if (b && m->decomposition().blocks.empty()) {
        if (PyErr_WarnFormat(PyExc_RuntimeWarning, 1,
                             "%s produced no blocks; BlockDecomposition.report['messages'] says why",
                             S->name) < 0) {
            Py_DECREF(b);
            return nullptr;
        }
    }
    return b;
}

// The five methods, in crossgen.methods order, and the Mesh method each is.
struct MeshMethod {
    const MethodSpec *spec;
    PyObject *(*function)(PyObject *, PyObject *, PyObject *);
};
const MeshMethod kMethods[] = {
    {&ziplineSpec, Mesh_runMethod<&ziplineSpec>},
    {&umberSpec, Mesh_runMethod<&umberSpec>},
    {&meridianSpec, Mesh_runMethod<&meridianSpec>},
    {&torsionSpec, Mesh_runMethod<&torsionSpec>},
    {&atlasSpec, Mesh_runMethod<&atlasSpec>},
};
constexpr size_t kNumMethods = sizeof(kMethods) / sizeof(kMethods[0]);

// Mesh.shape_dna(report=False, **options)
PyObject *Mesh_shapeDNA(PyObject *self, PyObject *args, PyObject *kwargs) {
    if (PyTuple_GET_SIZE(args) != 0) {
        PyErr_SetString(PyExc_TypeError, "shape_dna() takes keyword arguments only");
        return nullptr;
    }
    PyObject *wantReport = nullptr;
    bool ok = true;
    PyObject *rest = popKeyword(kwargs, "report", &wantReport, ok);
    if (!ok) return nullptr;
    bool full = false;
    if (wantReport && fromPy(wantReport, full) < 0) {
        Py_XDECREF(rest);
        return nullptr;
    }
    shapedna::ShapeDNA::Options opts;
    const int rc = applyKwargs(rest, "shape_dna", {target(shapeDNAOptionsTable(), opts)});
    Py_XDECREF(rest);
    if (rc < 0) return nullptr;

    const std::shared_ptr<const Mesh> model = asMesh(self)->mesh;
    std::unique_ptr<shapedna::ShapeDNA> solver;
    bool converged = false;
    if (!compute([&] {
            solver = std::make_unique<shapedna::ShapeDNA>(*model, opts);
            converged = solver->run();
        }))
        return nullptr;

    const shapedna::ShapeDNA::Report &r = solver->getReport();
    if (!converged) {
        const std::string why = r.messages.empty() ? "" : " (" + r.messages.back() + ")";
        if (!r.ran || solver->dna().empty()) {
            PyErr_Format(PyExc_RuntimeError, "shape_dna: the eigenproblem could not be solved%s",
                         why.c_str());
            return nullptr;
        }
        if (PyErr_WarnFormat(PyExc_RuntimeWarning, 1,
                             "shape_dna: the eigensolver did not converge%s", why.c_str()) < 0)
            return nullptr;
    }

    PyObject *dna = doubleArray(solver->dna());
    if (!dna || !full) return dna;

    PyObject *rep = toDict({source(shapeDNAReportTable(), r)});
    PyObject *lambda = rep ? doubleArray(solver->eigenvalues()) : nullptr;
    bool built = lambda && PyDict_SetItemString(rep, "eigenvalues", lambda) == 0;
    Py_XDECREF(lambda);
    // Eigenfunction k at vertex v is functions[v + k |V|]: row k of a (K, |V|)
    // array.
    const std::vector<double> &f = solver->eigenfunctions();
    const Py_ssize_t nv = static_cast<Py_ssize_t>(model->vertices.size());
    if (built && !f.empty() && nv > 0) {
        PyObject *fa = doubleArray(f.data(), static_cast<Py_ssize_t>(f.size()) / nv, nv);
        built = fa && PyDict_SetItemString(rep, "eigenfunctions", fa) == 0;
        Py_XDECREF(fa);
    }
    if (!built) {
        Py_DECREF(dna);
        Py_XDECREF(rep);
        return nullptr;
    }
    PyObject *pair = PyTuple_Pack(2, dna, rep);
    Py_DECREF(dna);
    Py_DECREF(rep);
    return pair;
}

const char *kShapeDNADoc =
    "shape_dna(report=False, **options) -> numpy.ndarray\n\n"
    "The Shape-DNA (Reuter et al. 2006): the first `count` eigenvalues of the\n"
    "Laplacian on the whole domain (materials ignored), normalised for scale.\n"
    "Keywords are ShapeDNA::Options in snake_case: count=50, degree=3 (P1-P3),\n"
    "refine=1, boundary='dirichlet'|'neumann', normalization='area'|'none'|\n"
    "'first_eigenvalue'|'weyl_slope'|'weyl_ratio', tolerance, eigenfunctions=False, threads,\n"
    "...; crossgen.options('shape_dna') lists them. With report=True, returns\n"
    "(dna, report): the solver's report, the raw spectrum under 'eigenvalues',\n"
    "and with eigenfunctions=True a (count, num_vertices) array under\n"
    "'eigenfunctions'.";

// Mesh.boundary_features(report=False, **options)
PyObject *Mesh_boundaryFeatures(PyObject *self, PyObject *args, PyObject *kwargs) {
    if (PyTuple_GET_SIZE(args) != 0) {
        PyErr_SetString(PyExc_TypeError, "boundary_features() takes keyword arguments only");
        return nullptr;
    }
    PyObject *wantReport = nullptr;
    bool ok = true;
    PyObject *rest = popKeyword(kwargs, "report", &wantReport, ok);
    if (!ok) return nullptr;
    bool full = false;
    if (wantReport && fromPy(wantReport, full) < 0) {
        Py_XDECREF(rest);
        return nullptr;
    }
    BoundaryFeatures::Options opts;
    const int rc = applyKwargs(rest, "boundary_features", {target(boundaryFeaturesOptionsTable(), opts)});
    Py_XDECREF(rest);
    if (rc < 0) return nullptr;

    // The default options are the ones regularity() grades against, and are
    // made once per mesh; any others are made for this call.
    const BoundaryFeatures::Options defaults;
    std::shared_ptr<const BoundaryFeatures> f;
    if (opts.cornerAngle == defaults.cornerAngle && opts.curvedAngle == defaults.curvedAngle) {
        f = defaultFeatures(self);
        if (!f) return nullptr;
    } else {
        const std::shared_ptr<const Mesh> model = asMesh(self)->mesh;
        if (!compute([&] { f = std::make_shared<const BoundaryFeatures>(*model, opts); })) return nullptr;
    }

    PyObject *summary = toDict({source(boundaryFeaturesSummaryTable(), f->summary())});
    if (!summary || !full) return summary;

    std::vector<Point> at;
    std::vector<double> angle;
    std::vector<int> blocks, vertex, region;
    for (const BoundaryFeatures::Corner &c : f->corners()) {
        at.push_back(c.p);
        angle.push_back(c.angle * 180.0 / M_PI);
        blocks.push_back(c.blocks);
        vertex.push_back(c.vertex);
        region.push_back(c.region);
    }
    PyObject *rep = PyDict_New();
    PyObject *loops = rep ? PyList_New(static_cast<Py_ssize_t>(f->loops().size())) : nullptr;
    bool built = loops != nullptr;
    for (size_t l = 0; built && l < f->loops().size(); ++l) {
        PyObject *a = indexArray(f->loops()[l].vertices);
        if (!a) built = false;
        else PyList_SET_ITEM(loops, static_cast<Py_ssize_t>(l), a);
    }
    const std::pair<const char *, PyObject *> items[] = {
        {"corners", built ? pointArray(at) : nullptr},
        {"corner_angles", built ? doubleArray(angle) : nullptr},
        {"corner_blocks", built ? indexArray(blocks) : nullptr},
        {"corner_vertices", built ? indexArray(vertex) : nullptr},
        {"corner_regions", built ? indexArray(region) : nullptr},
        {"loops", loops},
    };
    for (const auto &kv : items) {
        built = built && kv.second && PyDict_SetItemString(rep, kv.first, kv.second) == 0;
        Py_XDECREF(kv.second);
    }
    if (!built) {
        Py_DECREF(summary);
        Py_XDECREF(rep);
        return nullptr;
    }
    PyObject *pair = PyTuple_Pack(2, summary, rep);
    Py_DECREF(summary);
    Py_DECREF(rep);
    return pair;
}

const char *kBoundaryFeaturesDoc =
    "boundary_features(report=False, **options) -> dict\n\n"
    "What a quad layout of this mesh owes to its boundary (mesh::BoundaryFeatures):\n"
    "corners (a boundary vertex turning by more than corner_angle=20 degrees), by\n"
    "their ideal block count k = round(2 theta/pi); holes and the Euler\n"
    "characteristic; corner_defect, the quarter turns the corners owe on their own;\n"
    "singularity_bound, T3's interior singularities at the ideal corner split; and\n"
    "minimum_defect, the fewest quarter turns any conforming quad layout can have,\n"
    "which BlockDecomposition.regularity() grades against. Also the isoperimetric\n"
    "ratio, the bent share of the boundary and the shortest corner-to-corner run\n"
    "over sqrt(area). Material regions count separately, interfaces as boundary.\n"
    "With report=True, returns (summary, report): each corner's position, angle in\n"
    "degrees, k, vertex and region, and each boundary loop's vertices.";

// Filled in readyTypes(): the five methods, shape_dna, boundary_features, the sentinel.
PyMethodDef meshMethods[kNumMethods + 3];

// ===========================================================================
// BlockDecomposition
// ===========================================================================
void Blocks_dealloc(PyObject *self) {
    asBlocks(self)->features.~shared_ptr();
    asBlocks(self)->method.~shared_ptr();
    Py_TYPE(self)->tp_free(self);
}

const BlockDecomposition &decompOf(PyObject *self) { return asBlocks(self)->method->decomposition(); }

PyObject *Blocks_repr(PyObject *self) {
    const Method &m = *asBlocks(self)->method;
    const BlockDecomposition &D = m.decomposition();
    char cover[32];
    std::snprintf(cover, sizeof cover, "%.1f%%", 100.0 * m.coverage());
    return PyUnicode_FromFormat(
        "<crossgen.BlockDecomposition from %s: %zu blocks, %zu macro edges, %zu macrovertices, "
        "coverage %s>",
        m.name(), D.blocks.size(), D.edges.size(), D.vertices.size(), cover);
}

PyObject *Blocks_method(PyObject *self, void *) {
    return PyUnicode_FromString(asBlocks(self)->method->name());
}
PyObject *Blocks_numBlocks(PyObject *self, void *) { return PyLong_FromSize_t(decompOf(self).blocks.size()); }
PyObject *Blocks_coverage(PyObject *self, void *) {
    return PyFloat_FromDouble(asBlocks(self)->method->coverage());
}
PyObject *Blocks_vertices(PyObject *self, void *) {
    std::vector<Point> p;
    for (const BlockDecomposition::MacroVertex &v : decompOf(self).vertices) p.push_back(v.p);
    return pointArray(p);
}
PyObject *Blocks_blocks(PyObject *self, void *) {
    std::vector<std::array<int, 4>> c;
    for (const BlockDecomposition::Block &b : decompOf(self).blocks) c.push_back(b.corners);
    return indexArray(c);
}
PyObject *Blocks_blockEdges(PyObject *self, void *) {
    std::vector<std::array<int, 4>> e;
    for (const BlockDecomposition::Block &b : decompOf(self).blocks) e.push_back(b.edges);
    return indexArray(e);
}
PyObject *Blocks_blockMaterials(PyObject *self, void *) {
    std::vector<int> m;
    for (const BlockDecomposition::Block &b : decompOf(self).blocks) m.push_back(b.material);
    return indexArray(m);
}
PyObject *Blocks_edges(PyObject *self, void *) {
    const BlockDecomposition &D = decompOf(self);
    PyObject *list = PyList_New(static_cast<Py_ssize_t>(D.edges.size()));
    if (!list) return nullptr;
    for (size_t i = 0; i < D.edges.size(); ++i) {
        PyObject *a = pointArray(D.edges[i].points);
        if (!a) { Py_DECREF(list); return nullptr; }
        PyList_SET_ITEM(list, static_cast<Py_ssize_t>(i), a);
    }
    return list;
}
PyObject *Blocks_edgeVertices(PyObject *self, void *) {
    std::vector<std::array<int, 2>> e;
    for (const BlockDecomposition::MacroEdge &m : decompOf(self).edges) e.push_back({m.from, m.to});
    return indexArray(e);
}
PyObject *Blocks_report(PyObject *self, void *) { return asBlocks(self)->method->report(); }
PyObject *Blocks_options(PyObject *self, void *) { return asBlocks(self)->method->options(); }

// A decomposition of the model, by both tests: covers() -- every side on dS or
// shared -- and the area check, since a region with no block on any side of
// it leaves no unmatched side behind, but does leave area.
bool isValid(PyObject *self) {
    return asBlocks(self)->method->coverage() >= 0.999 && decompOf(self).covers();
}

PyObject *Blocks_valid(PyObject *self, void *) { return PyBool_FromLong(isValid(self) ? 1 : 0); }

PyGetSetDef blocksGetSet[] = {
    {"method", Blocks_method, nullptr, "The method that made it: 'zipline', 'umber', ...", nullptr},
    {"num_blocks", Blocks_numBlocks, nullptr, "Number of blocks.", nullptr},
    {"coverage", Blocks_coverage, nullptr,
     "Fraction of the model's area inside a block, by the method's own account.", nullptr},
    {"valid", Blocks_valid, nullptr,
     "True when the blocks are a decomposition of the model: every side on the\n"
     "boundary or shared by two blocks, and coverage at least 99.9%. The three\n"
     "qualities are NaN when it is False.",
     nullptr},
    {"vertices", Blocks_vertices, nullptr, "(k, 2) array of macrovertex positions.", nullptr},
    {"blocks", Blocks_blocks, nullptr,
     "(b, 4) array: each block's corners (macrovertex indices), in cyclic order.", nullptr},
    {"block_edges", Blocks_blockEdges, nullptr,
     "(b, 4) array: macro edge k runs from corner k to corner k+1 of its block.", nullptr},
    {"block_materials", Blocks_blockMaterials, nullptr, "(b,) array: each block's material.",
     nullptr},
    {"edges", Blocks_edges, nullptr,
     "List of (n, 2) arrays: each macro edge's polyline on the model, endpoints included.",
     nullptr},
    {"edge_vertices", Blocks_edgeVertices, nullptr,
     "(e, 2) array: the macrovertices each macro edge runs between.", nullptr},
    {"report", Blocks_report, nullptr, "The pipeline's status and reports, as a dict.", nullptr},
    {"options", Blocks_options, nullptr, "The options the method ran with, as a dict.", nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr},
};

PyObject *newQuadMesh(MeshOutput &&out, const char *method) {
    PyObject *o = QuadMeshType.tp_alloc(&QuadMeshType, 0);
    if (!o) {
        Py_XDECREF(out.report);
        return nullptr;
    }
    QuadMeshObject *q = asQuad(o);
    new (&q->mesh) std::unique_ptr<mesh::QuadMesh>(std::move(out.mesh));
    new (&q->smoothing) mesh::TMOP::Options(out.smoothing);
    new (&q->method) std::string(method);
    q->report = out.report;
    q->busy = false;
    out.report = nullptr;
    return o;
}

// BlockDecomposition.mesh(h=0.05, **options)
PyObject *Blocks_mesh(PyObject *self, PyObject *args, PyObject *kwargs) {
    const Method &method = *asBlocks(self)->method;
    PyObject *hArg = nullptr;
    bool ok = true;
    PyObject *rest = popKeyword(kwargs, "h", &hArg, ok);
    if (!ok) return nullptr;
    const Py_ssize_t npos = PyTuple_GET_SIZE(args);
    if (npos > 1 || (npos == 1 && hArg)) {
        Py_XDECREF(rest);
        PyErr_SetString(PyExc_TypeError, "mesh() takes one target edge length, h");
        return nullptr;
    }
    if (npos == 1) hArg = PyTuple_GET_ITEM(args, 0);
    double h = method.defaultTarget();
    if (hArg && fromPy(hArg, h) < 0) {
        Py_XDECREF(rest);
        return nullptr;
    }
    if (!(h > 0.0) || !std::isfinite(h)) {
        Py_XDECREF(rest);
        PyErr_Format(PyExc_ValueError, "mesh(): h must be a positive length, got %R", hArg);
        return nullptr;
    }
    MeshOutput out;
    const int rc = method.mesh(h, rest, out);
    Py_XDECREF(rest);
    if (rc < 0) {
        Py_XDECREF(out.report);
        return nullptr;
    }
    return newQuadMesh(std::move(out), method.name());
}

PyObject *Blocks_writeOBJ(PyObject *self, PyObject *args) {
    std::string path;
    if (!pathArg(args, "write_obj", path)) return nullptr;
    return writeResult(decompOf(self).writeEdgesOBJ(path), path);
}

PyObject *Blocks_writeBlocksOBJ(PyObject *self, PyObject *args) {
    std::string path;
    if (!pathArg(args, "write_blocks_obj", path)) return nullptr;
    return writeResult(decompOf(self).writeBlocksOBJ(path), path);
}

// NaN, not 0, for blocks that are no decomposition of the model: whether a
// method produced one is `valid`'s question, and a score of 0 would answer it
// a second time in every quality column.
constexpr double kNotValid = std::numeric_limits<double>::quiet_NaN();

PyObject *Blocks_regularity(PyObject *self, PyObject *) {
    const BlocksObject *b = asBlocks(self);
    return PyFloat_FromDouble(isValid(self) ? b->method->decomposition().regularity(*b->features) : kNotValid);
}

PyObject *Blocks_angleQuality(PyObject *self, PyObject *) {
    return PyFloat_FromDouble(isValid(self) ? decompOf(self).angleQuality() : kNotValid);
}

PyObject *Blocks_chordQuality(PyObject *self, PyObject *) {
    return PyFloat_FromDouble(isValid(self) ? decompOf(self).chordQuality() : kNotValid);
}

PyMethodDef blocksMethods[] = {
    {"mesh", asCFunction(Blocks_mesh), METH_VARARGS | METH_KEYWORDS,
     "mesh(h=0.05, **options) -> QuadMesh\n\n"
     "A quadrilateral mesh on the blocks at target edge length h, built by the\n"
     "method's own mesher: BlockQuadMesh for zipline and umber, MERIDIAN's\n"
     "Stage 10 (and Stage 11's O-grids) for meridian and torsion, BlockMesh for\n"
     "atlas. Keywords are that mesher's options (min_intervals, max_intervals,\n"
     "smoothing_passes, ...) and mesh::QuadMesh::Options (corner_angle,\n"
     "fix_all_feature_nodes, curve_source, ...), which decide which nodes the\n"
     "smoother may move; crossgen.options(method, 'mesh') lists them."},
    {"write_obj", Blocks_writeOBJ, METH_VARARGS,
     "write_obj(path): the macro edges as OBJ polylines."},
    {"write_blocks_obj", Blocks_writeBlocksOBJ, METH_VARARGS,
     "write_blocks_obj(path): each block as a closed OBJ loop."},
    {"regularity", Blocks_regularity, METH_NOARGS,
     "regularity() -> float in (0, 1], higher better; NaN when not valid\n\n"
     "Irregular macrovertices (T2, T4 of docs/block_decomposition_metrics.md),\n"
     "each weighted by the quarter turns its block count is off by: |4 - d|\n"
     "inside, |3 - d| on smooth boundary, |2 theta/pi - n| at a corner of angle\n"
     "theta in n blocks, and |2 theta/pi - 2| at a model corner no macrovertex\n"
     "sits on; interfaces act as boundary. W is the sum, and the result\n"
     "(1 + W_min) / (1 + W), W_min being the fewest quarter turns any layout of\n"
     "this model can have (Mesh.boundary_features()['minimum_defect']): 1 for a\n"
     "layout as regular as the model allows."},
    {"angle_quality", Blocks_angleQuality, METH_NOARGS,
     "angle_quality() -> float in [0, 1], higher better; NaN when not valid\n\n"
     "Block corner angles (B1): 1 - (area-weighted RMS deviation from the even\n"
     "split of their vertex's sector, 2 pi/d inside, theta/n at a corner) / 90\n"
     "degrees, floored at 0. How evenly the blocks are spread, not what the\n"
     "valence forces; read on the decomposition as drawn, before smoothing,\n"
     "which moves it -- a diagnostic of the method (B4), not a ranking score."},
    {"chord_quality", Blocks_chordQuality, METH_NOARGS,
     "chord_quality() -> float in (0, 1], higher better; NaN when not valid\n\n"
     "Chord span (T5): 1/S, S the geometric mean over chords, weighted by the\n"
     "area each runs through, of longest over shortest macro edge -- the\n"
     "element-size ratio its one interval count forces. 1 when every chord's\n"
     "sides are equal, 1/2 at a factor of two."},
    {nullptr, nullptr, 0, nullptr},
};

// ===========================================================================
// QuadMesh
// ===========================================================================
void Quad_dealloc(PyObject *self) {
    QuadMeshObject *q = asQuad(self);
    Py_XDECREF(q->report);
    q->mesh.~unique_ptr();
    q->smoothing.~Options();
    q->method.~basic_string();
    Py_TYPE(self)->tp_free(self);
}

// The mesh, or null with a RuntimeError while another thread smooths it.
mesh::QuadMesh *quadOf(PyObject *self) {
    QuadMeshObject *q = asQuad(self);
    if (q->busy) {
        PyErr_SetString(PyExc_RuntimeError, "this QuadMesh is being smoothed in another thread");
        return nullptr;
    }
    return q->mesh.get();
}

PyObject *Quad_repr(PyObject *self) {
    mesh::QuadMesh *m = quadOf(self);
    if (!m) return nullptr;
    m->computeQuality();
    char sj[32];
    std::snprintf(sj, sizeof sj, "%.4f", m->quality.minScaledJacobian);
    return PyUnicode_FromFormat(
        "<crossgen.QuadMesh from %s: %zu quads, %zu vertices, min scaled Jacobian %s>",
        asQuad(self)->method.c_str(), m->quads.size(), m->vertices.size(), sj);
}

PyObject *Quad_method(PyObject *self, void *) { return PyUnicode_FromString(asQuad(self)->method.c_str()); }
PyObject *Quad_vertices(PyObject *self, void *) {
    mesh::QuadMesh *m = quadOf(self);
    return m ? pointArray(m->vertices) : nullptr;
}
PyObject *Quad_quads(PyObject *self, void *) {
    mesh::QuadMesh *m = quadOf(self);
    return m ? indexArray(m->quads) : nullptr;
}
PyObject *Quad_materials(PyObject *self, void *) {
    mesh::QuadMesh *m = quadOf(self);
    return m ? indexArray(m->quadMatId) : nullptr;
}
PyObject *Quad_numVertices(PyObject *self, void *) {
    mesh::QuadMesh *m = quadOf(self);
    return m ? PyLong_FromSize_t(m->vertices.size()) : nullptr;
}
PyObject *Quad_numQuads(PyObject *self, void *) {
    mesh::QuadMesh *m = quadOf(self);
    return m ? PyLong_FromSize_t(m->quads.size()) : nullptr;
}
PyObject *Quad_quality(PyObject *self, void *) {
    mesh::QuadMesh *m = quadOf(self);
    if (!m) return nullptr;
    m->computeQuality();
    return toDict({source(qualityTable(), m->quality)});
}
PyObject *Quad_scaledJacobian(PyObject *self, void *) {
    mesh::QuadMesh *m = quadOf(self);
    if (!m) return nullptr;
    std::vector<double> sj(m->quads.size());
    for (size_t i = 0; i < sj.size(); ++i) sj[i] = m->minScaledJacobian(static_cast<int>(i));
    return doubleArray(sj);
}
PyObject *Quad_report(PyObject *self, void *) {
    PyObject *r = asQuad(self)->report;
    Py_INCREF(r);
    return r;
}

PyGetSetDef quadGetSet[] = {
    {"method", Quad_method, nullptr, "The method whose blocks were meshed.", nullptr},
    {"vertices", Quad_vertices, nullptr, "(n, 2) float64 array of node positions.", nullptr},
    {"quads", Quad_quads, nullptr, "(m, 4) int64 array of node indices, counter-clockwise.", nullptr},
    {"materials", Quad_materials, nullptr, "(m,) int64 array: each element's material id.", nullptr},
    {"num_vertices", Quad_numVertices, nullptr, "Number of nodes.", nullptr},
    {"num_quads", Quad_numQuads, nullptr, "Number of elements.", nullptr},
    {"quality", Quad_quality, nullptr,
     "Element quality now (mesh::QuadMesh::Quality): scaled Jacobians, areas, edges, node kinds.",
     nullptr},
    {"scaled_jacobian", Quad_scaledJacobian, nullptr,
     "(m,) array: each element's worst corner scaled Jacobian.", nullptr},
    {"report", Quad_report, nullptr, "The mesher's report, as a dict.", nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr},
};

// QuadMesh.smooth(niters=None, **options)
PyObject *Quad_smooth(PyObject *self, PyObject *args, PyObject *kwargs) {
    QuadMeshObject *q = asQuad(self);
    PyObject *nArg = nullptr;
    bool ok = true;
    PyObject *rest = popKeyword(kwargs, "niters", &nArg, ok);
    if (!ok) return nullptr;
    const Py_ssize_t npos = PyTuple_GET_SIZE(args);
    if (npos > 1 || (npos == 1 && nArg)) {
        Py_XDECREF(rest);
        PyErr_SetString(PyExc_TypeError, "smooth() takes one sweep count, niters");
        return nullptr;
    }
    if (npos == 1) nArg = PyTuple_GET_ITEM(args, 0);
    mesh::TMOP::Options t = q->smoothing;
    if (nArg && nArg != Py_None && fromPy(nArg, t.maxSweeps) < 0) {
        Py_XDECREF(rest);
        return nullptr;
    }
    if (t.maxSweeps < 0) {
        Py_XDECREF(rest);
        PyErr_SetString(PyExc_ValueError, "smooth(): niters must not be negative");
        return nullptr;
    }
    const int rc = applyKwargs(rest, "smooth", {target(tmopOptionsTable(), t)});
    Py_XDECREF(rest);
    if (rc < 0) return nullptr;
    if (!quadOf(self)) return nullptr;

    mesh::TMOP::Report report;
    q->busy = true;
    const bool done = compute([&] {
        mesh::TMOP smoother(*q->mesh, t);
        smoother.run();
        report = smoother.getReport();
    });
    q->busy = false;
    if (!done) return nullptr;
    return toDict({source(tmopReportTable(), report)});
}

PyObject *Quad_writeOBJ(PyObject *self, PyObject *args) {
    std::string path;
    mesh::QuadMesh *m = quadOf(self);
    if (!m || !pathArg(args, "write_obj", path)) return nullptr;
    return writeResult(m->writeOBJ(path), path);
}
PyObject *Quad_writeVTU(PyObject *self, PyObject *args) {
    std::string path;
    mesh::QuadMesh *m = quadOf(self);
    if (!m || !pathArg(args, "write_vtu", path)) return nullptr;
    return writeResult(m->writeVTU(path), path);
}
PyObject *Quad_writeMFEM(PyObject *self, PyObject *args) {
    std::string path;
    mesh::QuadMesh *m = quadOf(self);
    if (!m || !pathArg(args, "write_mfem", path)) return nullptr;
    return writeResult(m->writeMFEM(path), path);
}

PyMethodDef quadMethods[] = {
    {"smooth", asCFunction(Quad_smooth), METH_VARARGS | METH_KEYWORDS,
     "smooth(niters=1000, **options) -> dict\n\n"
     "TMOP (mesh::TMOP) on this mesh, in place, for at most niters sweeps; returns\n"
     "the smoother's report. Starts from the method's own settings (metric 7;\n"
     "mu at the element corners for zipline, umber and atlas, at the 2x2 Gauss\n"
     "points for meridian and torsion), and keywords override any field of\n"
     "mesh::TMOP::Options: metric, exponent, quadrature='corners'|'gauss2x2',\n"
     "untangle, threads, ...; crossgen.options(method, 'smooth') lists them.\n"
     "Calling it again smooths further from where the last call stopped."},
    {"write_obj", Quad_writeOBJ, METH_VARARGS,
     "write_obj(path): quad faces, `usemtl mat<id>` per material."},
    {"write_vtu", Quad_writeVTU, METH_VARARGS, "write_vtu(path): a VTK unstructured grid."},
    {"write_mfem", Quad_writeMFEM, METH_VARARGS,
     "write_mfem(path): an MFEM mesh, the material id as the element attribute."},
    {nullptr, nullptr, 0, nullptr},
};

// ===========================================================================
// Module functions
// ===========================================================================
PyObject *crossgen_load(PyObject *, PyObject *args) {
    std::string path;
    if (!pathArg(args, "load", path)) return nullptr;
    // Mesh(filename) reports an unreadable file as a runtime_error; asking
    // first gives Python its own FileNotFoundError / PermissionError.
    if (FILE *f = std::fopen(path.c_str(), "r")) {
        std::fclose(f);
    } else {
        return PyErr_SetFromErrnoWithFilename(PyExc_OSError, path.c_str());
    }
    std::shared_ptr<const Mesh> m;
    if (!compute([&] { m = std::make_shared<const Mesh>(path); })) return nullptr;
    if (m->triangles.empty()) {
        PyErr_Format(PyExc_ValueError, "'%s' has no triangles", path.c_str());
        return nullptr;
    }
    return newMesh(std::move(m));
}

PyObject *crossgen_options(PyObject *, PyObject *args, PyObject *kwargs) {
    static const char *kw[] = {"method", "stage", nullptr};
    const char *name = nullptr;
    const char *stage = "method";
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s|s:options", const_cast<char **>(kw), &name,
                                     &stage))
        return nullptr;
    const std::string n(name), s(stage);
    if (n == "shape_dna" || n == "boundary_features") {
        if (s != "method") {
            PyErr_Format(PyExc_ValueError, "crossgen.options(): '%s' has only the stage 'method'", name);
            return nullptr;
        }
        if (n == "shape_dna") return toDict({source(shapeDNAOptionsTable(), shapedna::ShapeDNA::Options())});
        return toDict({source(boundaryFeaturesOptionsTable(), BoundaryFeatures::Options())});
    }
    for (const MeshMethod &mm : kMethods)
        if (n == mm.spec->name) return mm.spec->defaults(s);
    PyErr_Format(PyExc_ValueError,
                 "crossgen.options(): unknown method '%s' (expected one of crossgen.methods, "
                 "'shape_dna' or 'boundary_features')",
                 name);
    return nullptr;
}

PyMethodDef moduleMethods[] = {
    {"load", crossgen_load, METH_VARARGS,
     "load(path) -> Mesh\n\n"
     "Read a triangle mesh from an .obj file, with the material id of each\n"
     "triangle from its `usemtl mat<id>` line (as Mesh2Dgmsh writes them)."},
    {"options", asCFunction(crossgen_options), METH_VARARGS | METH_KEYWORDS,
     "options(method, stage='method') -> dict\n\n"
     "The keywords a stage takes, with their defaults. method is one of\n"
     "crossgen.methods, 'shape_dna' or 'boundary_features'; stage is 'method'\n"
     "(Mesh.<method>()), 'mesh' (BlockDecomposition.mesh()) or 'smooth'\n"
     "(QuadMesh.smooth())."},
    {nullptr, nullptr, 0, nullptr},
};

const char *kModuleDoc =
    "CrossGen: cross fields and quad block decompositions of planar\n"
    "multi-material domains.\n\n"
    "    import crossgen\n"
    "    m = crossgen.load('data/meshes/singlemat/geom012.obj')\n"
    "    b = m.meridian()          # zipline, umber, meridian, torsion or atlas\n"
    "    q = b.mesh(0.05)          # quad mesh at a target edge length\n"
    "    q.smooth(1000)            # TMOP, in place\n"
    "    e = m.shape_dna()         # numpy array of the Shape-DNA\n"
    "    f = m.boundary_features() # corners, holes, T3's bound, as a dict\n\n"
    "Every option is a keyword with the C++ default; crossgen.options(method,\n"
    "stage) lists them. Calls release the GIL but run one at a time; use\n"
    "processes to run several models in parallel.";

PyModuleDef moduleDef = {
    PyModuleDef_HEAD_INIT, "crossgen", kModuleDoc, -1, moduleMethods,
    nullptr, nullptr, nullptr, nullptr,
};

int readyTypes() {
    MeshType.tp_name = "crossgen.Mesh";
    MeshType.tp_basicsize = sizeof(MeshObject);
    MeshType.tp_flags = Py_TPFLAGS_DEFAULT;
    MeshType.tp_doc =
        "Mesh(vertices, triangles, materials=None)\n\n"
        "A planar triangle mesh with a material id per triangle -- the input every\n"
        "method takes. crossgen.load(path) reads one from an .obj; the constructor\n"
        "builds one from an (n, 2) array of points, an (m, 3) array of vertex\n"
        "indices and optionally (m,) material ids (default: all 1). Triangles are\n"
        "turned counter-clockwise as the .obj reader turns them.";
    MeshType.tp_new = Mesh_new;
    MeshType.tp_dealloc = Mesh_dealloc;
    MeshType.tp_repr = Mesh_repr;
    MeshType.tp_getset = meshGetSet;
    size_t k = 0;
    for (const MeshMethod &mm : kMethods)
        meshMethods[k++] = {mm.spec->name, asCFunction(mm.function), METH_VARARGS | METH_KEYWORDS,
                            mm.spec->doc};
    meshMethods[k++] = {"shape_dna", asCFunction(Mesh_shapeDNA), METH_VARARGS | METH_KEYWORDS,
                        kShapeDNADoc};
    meshMethods[k++] = {"boundary_features", asCFunction(Mesh_boundaryFeatures),
                        METH_VARARGS | METH_KEYWORDS, kBoundaryFeaturesDoc};
    meshMethods[k] = {nullptr, nullptr, 0, nullptr};
    MeshType.tp_methods = meshMethods;

    BlocksType.tp_name = "crossgen.BlockDecomposition";
    BlocksType.tp_basicsize = sizeof(BlocksObject);
    BlocksType.tp_flags = Py_TPFLAGS_DEFAULT;
    BlocksType.tp_doc =
        "A quad block decomposition of a Mesh, as one method found it: macrovertices,\n"
        "macro edges (polylines on the model) and four-sided blocks. Made by\n"
        "Mesh.zipline(), .umber(), .meridian(), .torsion() and .atlas(); mesh() puts a\n"
        "quad mesh on it with that method's mesher.";
    BlocksType.tp_dealloc = Blocks_dealloc;
    BlocksType.tp_repr = Blocks_repr;
    BlocksType.tp_getset = blocksGetSet;
    BlocksType.tp_methods = blocksMethods;

    QuadMeshType.tp_name = "crossgen.QuadMesh";
    QuadMeshType.tp_basicsize = sizeof(QuadMeshObject);
    QuadMeshType.tp_flags = Py_TPFLAGS_DEFAULT;
    QuadMeshType.tp_doc =
        "A quadrilateral mesh on a BlockDecomposition (mesh::QuadMesh), made by\n"
        "BlockDecomposition.mesh(h). smooth(niters) runs TMOP on it in place.";
    QuadMeshType.tp_dealloc = Quad_dealloc;
    QuadMeshType.tp_repr = Quad_repr;
    QuadMeshType.tp_getset = quadGetSet;
    QuadMeshType.tp_methods = quadMethods;

    if (PyType_Ready(&MeshType) < 0 || PyType_Ready(&BlocksType) < 0 ||
        PyType_Ready(&QuadMeshType) < 0)
        return -1;
    return 0;
}

PyObject *createModule() {
    try {
        buildSharedTables();
        for (const MeshMethod &mm : kMethods) mm.spec->buildTables();
    } catch (const std::exception &e) {
        PyErr_Format(PyExc_ImportError, "crossgen: %s", e.what());
        return nullptr;
    }
    if (readyTypes() < 0) return nullptr;

    PyObject *m = PyModule_Create(&moduleDef);
    if (!m) return nullptr;
    PyObject *names = PyTuple_New(static_cast<Py_ssize_t>(kNumMethods));
    for (size_t i = 0; names && i < kNumMethods; ++i) {
        PyObject *s = PyUnicode_FromString(kMethods[i].spec->name);
        if (!s) { Py_CLEAR(names); break; }
        PyTuple_SET_ITEM(names, static_cast<Py_ssize_t>(i), s);
    }
    const bool ok = names && PyModule_AddObjectRef(m, "Mesh", reinterpret_cast<PyObject *>(&MeshType)) == 0 &&
                    PyModule_AddObjectRef(m, "BlockDecomposition", reinterpret_cast<PyObject *>(&BlocksType)) == 0 &&
                    PyModule_AddObjectRef(m, "QuadMesh", reinterpret_cast<PyObject *>(&QuadMeshType)) == 0 &&
                    PyModule_AddObjectRef(m, "methods", names) == 0 &&
                    PyModule_AddStringConstant(m, "__version__", "0.1.0") == 0;
    Py_XDECREF(names);
    if (!ok) {
        Py_DECREF(m);
        return nullptr;
    }
    return m;
}

}  // namespace
}  // namespace pycg

PyMODINIT_FUNC PyInit_crossgen(void) { return pycg::createModule(); }
