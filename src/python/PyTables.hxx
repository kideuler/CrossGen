#ifndef __PY_TABLES_HXX__
#define __PY_TABLES_HXX__

#include "python/PyOptions.hxx"

#include "ShapeDNA/ShapeDNA.hxx"
#include "dualmbo/DualMBO.hxx"
#include "mesh/BlockQuadMesh.hxx"
#include "mesh/QuadMesh.hxx"
#include "mesh/TMOP.hxx"

// The tables more than one method's binding reads: the mesh a smoother works
// on, the smoother itself, the method-agnostic block mesher, and ShapeDNA.
// Each method's own tables live beside its binding (Method*.cxx).
namespace pycg {

// ---- enumerations shared between methods --------------------------------------
template <>
struct EnumNames<mesh::TMOP::Metric> {
    static constexpr std::pair<mesh::TMOP::Metric, const char *> list[] = {
        {mesh::TMOP::Shape002, "shape002"},
        {mesh::TMOP::Shape004, "shape004"},
        {mesh::TMOP::ShapeSize007, "shape_size007"},
        {mesh::TMOP::Untangle022, "untangle022"},
        {mesh::TMOP::Size055, "size055"},
        {mesh::TMOP::Size056, "size056"},
        {mesh::TMOP::ShapeSizeCombo, "shape_size_combo"},
        {mesh::TMOP::UntangleRegularized, "untangle_regularized"},
    };
};

template <>
struct EnumNames<mesh::TMOP::Target> {
    static constexpr std::pair<mesh::TMOP::Target, const char *> list[] = {
        {mesh::TMOP::TargetKeep, "keep"},
        {mesh::TMOP::TargetUniformSquare, "uniform_square"},
        {mesh::TMOP::TargetCurrentShape, "current_shape"},
        {mesh::TMOP::TargetPerQuadSize, "per_quad_size"},
    };
};

template <>
struct EnumNames<mesh::TMOP::Quadrature> {
    static constexpr std::pair<mesh::TMOP::Quadrature, const char *> list[] = {
        {mesh::TMOP::Gauss2x2, "gauss2x2"},
        {mesh::TMOP::Corners, "corners"},
    };
};

template <>
struct EnumNames<mesh::TMOP::Untangler> {
    static constexpr std::pair<mesh::TMOP::Untangler, const char *> list[] = {
        {mesh::TMOP::UntangleShifted, "shifted"},
        {mesh::TMOP::UntangleRegular, "regular"},
    };
};

template <>
struct EnumNames<mesh::QuadMesh::Options::CurveSource> {
    static constexpr std::pair<mesh::QuadMesh::Options::CurveSource, const char *> list[] = {
        {mesh::QuadMesh::Options::CurveChord, "chord"},
        {mesh::QuadMesh::Options::CurvePolyline, "polyline"},
        {mesh::QuadMesh::Options::CurveSpline, "spline"},
    };
};

template <>
struct EnumNames<DualMBO::PenaltyWeight> {
    static constexpr std::pair<DualMBO::PenaltyWeight, const char *> list[] = {
        {DualMBO::PenaltyWeight::MinHeight, "min_height"},
        {DualMBO::PenaltyWeight::HarmonicHeight, "harmonic_height"},
        {DualMBO::PenaltyWeight::Orthogonal, "orthogonal"},
    };
};

template <>
struct EnumNames<shapedna::ShapeDNA::Boundary> {
    static constexpr std::pair<shapedna::ShapeDNA::Boundary, const char *> list[] = {
        {shapedna::ShapeDNA::Dirichlet, "dirichlet"},
        {shapedna::ShapeDNA::Neumann, "neumann"},
    };
};

template <>
struct EnumNames<shapedna::ShapeDNA::Normalization> {
    static constexpr std::pair<shapedna::ShapeDNA::Normalization, const char *> list[] = {
        {shapedna::ShapeDNA::NoNormalization, "none"},
        {shapedna::ShapeDNA::AreaNormalization, "area"},
        {shapedna::ShapeDNA::FirstEigenvalue, "first_eigenvalue"},
        {shapedna::ShapeDNA::WeylSlope, "weyl_slope"},
        {shapedna::ShapeDNA::WeylRatio, "weyl_ratio"},
    };
};

// ---- tables --------------------------------------------------------------------

// mesh::QuadMesh::Options: which nodes of a finished quad mesh are free,
// sliding or fixed, which decides what the smoother may move. Taken by every
// method's BlockDecomposition.mesh(), since it is fixed when the mesh is built.
const Table<mesh::QuadMesh::Options> &nodeOptionsTable();
const Table<mesh::QuadMesh::Quality> &qualityTable();

// mesh::TMOP::Options, less maxSweeps (QuadMesh.smooth's `niters`) and
// perQuadSize (a per-element vector no keyword can sensibly carry).
const Table<mesh::TMOP::Options> &tmopOptionsTable();
const Table<mesh::TMOP::Report> &tmopReportTable();

// BlockQuadMesh::Options, less targetEdgeLength (BlockDecomposition.mesh's `h`).
const Table<BlockQuadMesh::Options> &blockQuadMeshOptionsTable();
const Table<BlockQuadMesh::Report> &blockQuadMeshReportTable();

const Table<shapedna::ShapeDNA::Options> &shapeDNAOptionsTable();
const Table<shapedna::ShapeDNA::Report> &shapeDNAReportTable();

// The TMOP settings for a transfinite grid on a block decomposition -- ATLAS,
// UMBER and ZIPLINE, and the viewer's copy of each: metric 7, mu sampled at
// the element corners, 1000 sweeps. Corners because at the 2x2 Gauss points a
// corner of such a grid can turn over without the barrier seeing it (measured:
// six ATLAS corpus models, and UMBER geom012 0.033 -> 0.012).
mesh::TMOP::Options blockGridSmoothing();

// crossgen.options(<method>, "mesh"): the mesher's table and the node options'
// over their defaults, with the default target edge length under "h".
PyObject *meshingDefaults(double h, std::initializer_list<Source> sources);

// crossgen.options(<method>, "smooth"): the TMOP table over `defaults`, with
// the default sweep count under "niters".
PyObject *smoothingDefaults(const mesh::TMOP::Options &defaults);

// Build every table once, so that a table defect (two fields one keyword)
// fails `import crossgen` instead of the first call that happens to use it.
void buildSharedTables();

}  // namespace pycg

#endif // __PY_TABLES_HXX__
