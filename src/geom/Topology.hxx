#ifndef __GEOM_TOPOLOGY_HXX__
#define __GEOM_TOPOLOGY_HXX__

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "geom/BSpline.hxx"
#include "geom/Vec.hxx"

// The boundary representation: vertices, edges, faces and the shapes they make,
// as OpenCASCADE topology (TopoDS) behind value types that name none of it.
//
// A curve says where something is; topology says what is connected to what.
// Two faces are joined along an edge not because their boundaries happen to
// coincide but because they hold *the same edge* -- the same kernel object --
// and two edges meet at a vertex because they hold the same vertex. That is the
// form of the watertightness rule SplineFit enforces on its control points, one
// level up: build a Vertex once per node and an Edge once per arc, hand the
// same objects to every edge and face that uses them, and the shape that comes
// out is closed by construction. Shape::freeEdgeCount() then says so rather
// than a tolerance.
//
// Everything is three-dimensional. A 2-D curve or surface given to a
// constructor is placed in the plane z = 0 (lift() in BSpline.hxx), and points
// come back as Vec<3>.
//
// Each class is a handle: copies are the same sub-shape, isSame() says whether
// two handles are, and nothing here modifies a shape after it is built.
namespace geom {

namespace detail {
struct ShapeData;
struct Access;
} // namespace detail

class Vertex {
public:
    Vertex() = default;
    explicit Vertex(const Vec<2> &point);
    explicit Vertex(const Vec<3> &point);

    bool isNull() const { return !data; }
    Vec<3> point() const;
    bool isSame(const Vertex &other) const;

private:
    std::shared_ptr<const detail::ShapeData> data;
    friend struct detail::Access;
};

class Edge {
public:
    Edge() = default;

    // The whole of `curve`, from `start` at the beginning of its domain to `end`
    // at the end of it. The two vertices must be where the curve ends, to the
    // kernel's tolerance (1e-7); they may be the same vertex, for a closed
    // curve. Throws std::invalid_argument otherwise, and on a curve shorter
    // than that tolerance.
    Edge(const BSplineCurve<2> &curve, const Vertex &start, const Vertex &end);
    Edge(const BSplineCurve<3> &curve, const Vertex &start, const Vertex &end);
    // The same with vertices of its own, which no other edge can then share.
    explicit Edge(const BSplineCurve<2> &curve);
    explicit Edge(const BSplineCurve<3> &curve);

    bool isNull() const { return !data; }
    Vertex start() const;
    Vertex end() const;
    // The curve the edge lies on, sharing its kernel object.
    BSplineCurve<3> curve() const;
    double parameterBegin() const;
    double parameterEnd() const;
    Vec<3> evaluate(double u) const;
    double length() const;
    bool isSame(const Edge &other) const;

private:
    std::shared_ptr<const detail::ShapeData> data;
    friend struct detail::Access;
};

class Face {
public:
    Face() = default;

    // The whole of `surface`, bounded by four edges that lie along its four
    // sides. The sides are in the order and orientation of Coons.hxx --
    // bottom (t = t0, s increasing), right (s = s1, t increasing), top (t = t1,
    // s increasing), left (s = s0, t increasing) -- and alongFrame[k] says
    // whether edge k runs that way or the opposite way. Each edge's curve must
    // trace its side at the same parameter, up to that reversal, and the four
    // must meet at shared vertices at the corners.
    //
    // Throws std::invalid_argument if the corners do not share vertices, if an
    // edge is used for two sides, or if an edge is not on its side -- by more
    // than 1e-9 of the surface's size, sampled.
    Face(const BSplineSurface<2> &surface, const std::array<Edge, 4> &sides,
         const std::array<bool, 4> &alongFrame);
    Face(const BSplineSurface<3> &surface, const std::array<Edge, 4> &sides,
         const std::array<bool, 4> &alongFrame);

    bool isNull() const { return !data; }
    BSplineSurface<3> surface() const;
    // Its boundary edges, in the order the boundary visits them.
    std::vector<Edge> edges() const;
    Vec<3> evaluate(double s, double t) const;
    double area() const;
    // The kernel's validity check (BRepCheck): the edges lie on the surface at
    // the parameters the face claims, the boundary closes, and so on.
    bool isValid() const;

private:
    std::shared_ptr<const detail::ShapeData> data;
    friend struct detail::Access;
};

// A collection of faces, and of edges that bound none of them, as one compound.
//
// Faces joined by shared edges are gathered into shells, one per connected
// piece, and turned where needed so that each shared edge is run one way by one
// face and the other way by the other -- the orientation a surface model has to
// have, and what lets a STEP file keep the edges shared (a face left loose in a
// compound is written as a shell of its own, with copies of its edges).
class Shape {
public:
    Shape() = default;
    explicit Shape(const std::vector<Face> &faces, const std::vector<Edge> &looseEdges = {});

    // A shape read back from a file this class wrote, or from any other BREP or
    // STEP source. Throws std::runtime_error if the file cannot be read. Faces
    // on surfaces other than B-splines read in, but Face::surface() is empty
    // for them.
    static Shape readBREP(const std::string &filename);
    static Shape readSTEP(const std::string &filename);

    bool isNull() const { return !data; }
    // Every distinct face and edge, the edges including those on faces.
    std::vector<Face> faces() const;
    std::vector<Edge> edges() const;
    int faceCount() const;
    // Distinct sub-shapes: an edge two faces share counts once.
    int edgeCount() const;
    int vertexCount() const;
    // Edges that bound exactly one face -- the boundary of the surface, on a
    // shape whose faces are joined -- and edges that bound two or more.
    int freeEdgeCount() const;
    int sharedEdgeCount() const;
    double area() const;
    bool isValid() const;

    // The shape in OpenCASCADE's native format, and as a STEP file (AP214, the
    // geometry as it is, no unit conversion). False if the file could not be
    // written.
    bool writeBREP(const std::string &filename) const;
    bool writeSTEP(const std::string &filename) const;

private:
    std::shared_ptr<const detail::ShapeData> data;
    friend struct detail::Access;
};

} // namespace geom

#endif // __GEOM_TOPOLOGY_HXX__
