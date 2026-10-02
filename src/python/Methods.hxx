#ifndef __PY_METHODS_HXX__
#define __PY_METHODS_HXX__

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <memory>
#include <string>

#include "mesh/BlockDecomposition.hxx"
#include "mesh/Mesh.hxx"
#include "mesh/QuadMesh.hxx"
#include "mesh/TMOP.hxx"

// The five block-decomposition methods behind one interface, which is all the
// module's BlockDecomposition type knows of them.
//
// Each method reaches the shared BlockDecomposition (mesh/BlockDecomposition.
// hxx) by its own route, and each already has its own way of meshing it, which
// the binding keeps rather than replaces:
//
//   zipline   ZIPLINE, Stages 0b-5 (src/ZIPLINE)        -> BlockQuadMesh
//   umber     the UMBER stages as TestUMBER runs them    -> BlockQuadMesh
//   meridian  MERIDIAN, Stages 0-9 (Pipeline A)          -> MERIDIAN Stage 10 (+ 11)
//   torsion   TORSION, Stages 0-9 (Pipeline B)           -> MERIDIAN Stage 10 (+ 11)
//   atlas     ATLAS, Stages 1-6                          -> ATLAS BlockMesh
//
// That last column is not a detail. BlockQuadMesh meshes the decomposition's
// polylines; MERIDIAN's Stage 10 evaluates the bicubic patches Stage 9 fitted
// and BlockMesh evaluates ATLAS's certified charts, exact interior geometry the
// decomposition does not carry (see BlockQuadMesh.hxx, "Why this one exists").
// Meshing a MERIDIAN layout through BlockQuadMesh would be a worse mesh of the
// same layout, so the binding holds on to each pipeline after it has run and
// meshes through the pipeline's own stage.
//
// The same goes for smoothing: each method's meshes carry the TMOP settings its
// own driver and viewer mode smooth with (metric 7 everywhere, mu at the
// element corners everywhere -- MERIDIAN and TORSION sampled the 2x2 Gauss
// points until 2026-10-02, when their Stage 12 changed), and MERIDIAN's and
// TORSION's also pillow the flat feature corners first, as their Stage 12
// does since the same day (mesh::Pillow).
namespace pycg {

// What BlockDecomposition.mesh() hands to the module: the quad mesh as the
// smoother will see it, the mesher's report, and the method's TMOP settings.
struct MeshOutput {
    std::unique_ptr<mesh::QuadMesh> mesh;
    PyObject *report = nullptr;           // owned: a dict
    mesh::TMOP::Options smoothing;
    // Whether the first smooth() pillows the flat feature corners before it
    // smooths: Stage 12 of MERIDIAN and TORSION.
    bool pillow = false;
};

class Method {
public:
    virtual ~Method() = default;

    // "zipline", "umber", ... -- the Mesh method that made it.
    virtual const char *name() const = 0;

    virtual const BlockDecomposition &decomposition() const = 0;

    // Fraction of the model's area inside a block, by the method's own
    // account (ZIPLINE: LayoutBlocks::coverage; UMBER: BlockLayout::coverageOf;
    // MERIDIAN/TORSION: 1 - MERIDIAN::unmeshableFraction; ATLAS: the blocks'
    // area over the model's).
    virtual double coverage() const = 0;

    // The pipeline's status and reports, and the options it ran with, as
    // dicts. New references, or null with an exception set.
    virtual PyObject *report() const = 0;
    virtual PyObject *options() const = 0;

    // Mesh the decomposition at target edge length `h`, with the method's
    // mesher; `kwargs` (may be null) are the mesher's options and
    // mesh::QuadMesh::Options. Returns 0, or -1 with an exception set.
    virtual int mesh(double h, PyObject *kwargs, MeshOutput &out) const = 0;

    // The `h` mesh() takes when none is given: the mesher's own default.
    virtual double defaultTarget() const = 0;
};

// One Mesh method: its name, its docstring, how to run it, and its defaults.
struct MethodSpec {
    const char *name;
    const char *doc;
    // Run the method on a copy of `mesh` with the options in `kwargs` (may be
    // null). Null with an exception set on a bad keyword or a C++ failure.
    std::shared_ptr<Method> (*run)(const Mesh &mesh, PyObject *kwargs);
    // crossgen.options(name, stage): "method", "mesh" or "smooth". Null with
    // an exception set for an unknown stage.
    PyObject *(*defaults)(const std::string &stage);
    // Build this method's tables (see buildSharedTables).
    void (*buildTables)();
};

extern const MethodSpec ziplineSpec;
extern const MethodSpec umberSpec;
extern const MethodSpec meridianSpec;
extern const MethodSpec torsionSpec;
extern const MethodSpec atlasSpec;

// A ValueError for crossgen.options() given a stage `method` has no table for.
PyObject *unknownStage(const char *method, const std::string &stage);

// Fraction of `mesh`'s area inside D's blocks, each block's outline read as the
// polygon of its four sides.
double decompositionCoverage(const BlockDecomposition &D, const Mesh &mesh);

}  // namespace pycg

#endif // __PY_METHODS_HXX__
