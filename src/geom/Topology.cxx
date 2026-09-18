#include "geom/Topology.hxx"

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Message.hxx>
#include <Message_Messenger.hxx>
#include <Precision.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <TColgp_Array1OfPnt2d.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>

#include <algorithm>
#include <cmath>
#include <deque>
#include <stdexcept>

#include "geom/detail/Occ.hxx"

namespace geom {

using detail::Access;

namespace {

TopoDS_Vertex vertexOf(const Vertex &v) { return TopoDS::Vertex(Access::shape(v)); }
TopoDS_Edge edgeOf(const Edge &e) { return TopoDS::Edge(Access::shape(e)); }

Edge makeEdge(const Handle(Geom_BSplineCurve) &curve, const Vertex *start, const Vertex *end) {
    if (curve.IsNull()) throw std::invalid_argument("Edge: an empty curve");
    if (start && (start->isNull() || end->isNull())) throw std::invalid_argument("Edge: a null vertex");
    return detail::construct("Edge", [&] {
        const double u0 = curve->FirstParameter(), u1 = curve->LastParameter();
        BRepBuilderAPI_MakeEdge make = start
            ? BRepBuilderAPI_MakeEdge(curve, vertexOf(*start), vertexOf(*end), u0, u1)
            : BRepBuilderAPI_MakeEdge(curve, u0, u1);
        if (!make.IsDone()) {
            throw std::invalid_argument("Edge: the kernel refused the edge (BRepBuilderAPI_EdgeError " +
                                        std::to_string(static_cast<int>(make.Error())) +
                                        "); are the vertices at the ends of the curve?");
        }
        // A curve shorter than the vertex tolerance comes back as a degenerated
        // edge, with no curve at all -- the kind a sphere has at its pole, which
        // is only valid as part of a face.
        if (BRep_Tool::Degenerated(make.Edge())) {
            throw std::invalid_argument("Edge: the curve has no length to be an edge of");
        }
        return Access::wrap<Edge>(make.Edge());
    });
}

// The straight pcurve of one side of a face: from `from` to `to` in the
// surface's (s, t), over the edge's own parameter range.
Handle(Geom2d_BSplineCurve) sideCurve(const gp_Pnt2d &from, const gp_Pnt2d &to, double u0, double u1) {
    TColgp_Array1OfPnt2d poles(1, 2);
    poles(1) = from;
    poles(2) = to;
    TColStd_Array1OfReal knots(1, 2);
    knots(1) = u0;
    knots(2) = u1;
    TColStd_Array1OfInteger mults(1, 2);
    mults(1) = 2;
    mults(2) = 2;
    return new Geom2d_BSplineCurve(poles, knots, mults, 1);
}

Face makeFace(const Handle(Geom_BSplineSurface) &surface, const std::array<Edge, 4> &sides,
              const std::array<bool, 4> &along) {
    if (surface.IsNull()) throw std::invalid_argument("Face: an empty surface");
    for (int k = 0; k < 4; ++k) {
        if (sides[k].isNull()) throw std::invalid_argument("Face: a null edge");
        for (int j = 0; j < k; ++j) {
            if (sides[j].isSame(sides[k])) {
                throw std::invalid_argument("Face: one edge is given for two sides");
            }
        }
    }

    double s0, s1, t0, t1;
    surface->Bounds(s0, s1, t0, t1);
    // Side k's two ends in (s, t), in the direction Coons.hxx runs it.
    const std::array<std::array<gp_Pnt2d, 2>, 4> frame{{
        {{gp_Pnt2d(s0, t0), gp_Pnt2d(s1, t0)}},   // bottom
        {{gp_Pnt2d(s1, t0), gp_Pnt2d(s1, t1)}},   // right
        {{gp_Pnt2d(s0, t1), gp_Pnt2d(s1, t1)}},   // top
        {{gp_Pnt2d(s0, t0), gp_Pnt2d(s0, t1)}},   // left
    }};

    // The corners: where each side begins and ends in the frame has to be one
    // vertex, shared with the neighbouring side.
    auto frameStart = [&](int k) { return along[k] ? sides[k].start() : sides[k].end(); };
    auto frameEnd = [&](int k) { return along[k] ? sides[k].end() : sides[k].start(); };
    if (!frameStart(0).isSame(frameStart(3)) || !frameEnd(0).isSame(frameStart(1)) ||
        !frameEnd(1).isSame(frameEnd(2)) || !frameStart(2).isSame(frameEnd(3))) {
        throw std::invalid_argument("Face: the four edges do not meet at shared vertices");
    }

    // Each edge on its side, at the same parameter. Sampled, against the size
    // of the control net.
    double extent = 0.0;
    {
        const TColgp_Array2OfPnt &poles = surface->Poles();
        gp_Pnt lo = poles(1, 1), hi = poles(1, 1);
        for (int i = poles.LowerRow(); i <= poles.UpperRow(); ++i) {
            for (int j = poles.LowerCol(); j <= poles.UpperCol(); ++j) {
                const gp_Pnt &p = poles(i, j);
                lo.SetCoord(std::min(lo.X(), p.X()), std::min(lo.Y(), p.Y()), std::min(lo.Z(), p.Z()));
                hi.SetCoord(std::max(hi.X(), p.X()), std::max(hi.Y(), p.Y()), std::max(hi.Z(), p.Z()));
            }
        }
        extent = lo.Distance(hi);
    }
    const double tolerance = std::max(1e-9 * extent, 1e-12);
    for (int k = 0; k < 4; ++k) {
        const Handle(Geom_BSplineCurve) curve = Access::curve(sides[k].curve());
        const double u0 = curve->FirstParameter(), u1 = curve->LastParameter();
        for (int i = 0; i <= 8; ++i) {
            const double w = i / 8.0;
            const gp_Pnt2d st(frame[k][0].X() + w * (frame[k][1].X() - frame[k][0].X()),
                              frame[k][0].Y() + w * (frame[k][1].Y() - frame[k][0].Y()));
            const double u = along[k] ? u0 + w * (u1 - u0) : u1 - w * (u1 - u0);
            const double gap = surface->Value(st.X(), st.Y()).Distance(curve->Value(u));
            if (!(gap <= tolerance)) {
                throw std::invalid_argument("Face: edge " + std::to_string(k) + " is " +
                                            std::to_string(gap) + " off its side of the surface");
            }
        }
    }

    return detail::construct("Face", [&] {
        const double tol = Precision::Confusion();
        BRep_Builder builder;
        TopoDS_Face face;
        builder.MakeFace(face, surface, tol);
        TopoDS_Wire wire;
        builder.MakeWire(wire);
        for (int k = 0; k < 4; ++k) {
            TopoDS_Edge e = edgeOf(sides[k]);
            double first = 0.0, last = 0.0;
            BRep_Tool::Range(e, first, last);
            const gp_Pnt2d &a = along[k] ? frame[k][0] : frame[k][1];
            const gp_Pnt2d &b = along[k] ? frame[k][1] : frame[k][0];
            builder.UpdateEdge(e, sideCurve(a, b, first, last), face, tol);
            // The boundary runs anticlockwise in (s, t): bottom and right with
            // the frame, top and left against it.
            const bool forward = (k < 2) == along[k];
            builder.Add(wire, e.Oriented(forward ? TopAbs_FORWARD : TopAbs_REVERSED));
        }
        builder.Add(face, wire);
        return Access::wrap<Face>(face);
    });
}

// OpenCASCADE's translators report to the default messenger, which prints to
// the console. The pipelines' drivers own the console, so the printers are
// taken away for the length of a write and put back after.
class QuietKernel {
public:
    QuietKernel() : saved(Message::DefaultMessenger()->Printers()) {
        Message::DefaultMessenger()->ChangePrinters().Clear();
    }
    ~QuietKernel() { Message::DefaultMessenger()->ChangePrinters() = saved; }

private:
    Message_SequenceOfPrinters saved;
};

// How `edge` is run by `face`, as the face is oriented.
TopAbs_Orientation orientationIn(const TopoDS_Shape &face, const TopoDS_Shape &edge) {
    for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next()) {
        if (it.Current().IsSame(edge)) return it.Current().Orientation();
    }
    return TopAbs_EXTERNAL;
}

int countShapes(const TopoDS_Shape &shape, TopAbs_ShapeEnum type) {
    TopTools_IndexedMapOfShape map;
    TopExp::MapShapes(shape, type, map);
    return map.Extent();
}

} // namespace

// ---------------------------------------------------------------------------
// Vertex
// ---------------------------------------------------------------------------
Vertex::Vertex(const Vec<2> &point) : Vertex(Vec<3>{point[0], point[1], 0.0}) {}

Vertex::Vertex(const Vec<3> &point) {
    *this = detail::construct("Vertex", [&] {
        return Access::wrap<Vertex>(BRepBuilderAPI_MakeVertex(detail::toPnt(point)).Vertex());
    });
}

Vec<3> Vertex::point() const {
    if (!data) return zero<3>();
    return detail::fromPnt<3>(BRep_Tool::Pnt(TopoDS::Vertex(data->shape)));
}

bool Vertex::isSame(const Vertex &other) const {
    return data && other.data && data->shape.IsSame(other.data->shape);
}

// ---------------------------------------------------------------------------
// Edge
// ---------------------------------------------------------------------------
Edge::Edge(const BSplineCurve<2> &curve, const Vertex &start, const Vertex &end)
    : Edge(lift(curve), start, end) {}

Edge::Edge(const BSplineCurve<3> &curve, const Vertex &start, const Vertex &end) {
    *this = makeEdge(Access::curve(curve), &start, &end);
}

Edge::Edge(const BSplineCurve<2> &curve) : Edge(lift(curve)) {}

Edge::Edge(const BSplineCurve<3> &curve) {
    *this = makeEdge(Access::curve(curve), nullptr, nullptr);
}

Vertex Edge::start() const {
    if (!data) return Vertex();
    TopoDS_Vertex first, last;
    TopExp::Vertices(TopoDS::Edge(data->shape), first, last, Standard_False);
    return Access::wrap<Vertex>(first);
}

Vertex Edge::end() const {
    if (!data) return Vertex();
    TopoDS_Vertex first, last;
    TopExp::Vertices(TopoDS::Edge(data->shape), first, last, Standard_False);
    return Access::wrap<Vertex>(last);
}

BSplineCurve<3> Edge::curve() const {
    if (!data) return BSplineCurve<3>();
    double first = 0.0, last = 0.0;
    const Handle(Geom_Curve) c = BRep_Tool::Curve(TopoDS::Edge(data->shape), first, last);
    return Access::curve<3>(Handle(Geom_BSplineCurve)::DownCast(c));
}

double Edge::parameterBegin() const {
    if (!data) return 0.0;
    double first = 0.0, last = 0.0;
    BRep_Tool::Range(TopoDS::Edge(data->shape), first, last);
    return first;
}

double Edge::parameterEnd() const {
    if (!data) return 0.0;
    double first = 0.0, last = 0.0;
    BRep_Tool::Range(TopoDS::Edge(data->shape), first, last);
    return last;
}

Vec<3> Edge::evaluate(double u) const {
    return curve().evaluate(std::max(parameterBegin(), std::min(parameterEnd(), u)));
}

double Edge::length() const {
    if (!data) return 0.0;
    return detail::guard("Edge::length", [&] {
        GProp_GProps props;
        BRepGProp::LinearProperties(data->shape, props);
        return props.Mass();
    });
}

bool Edge::isSame(const Edge &other) const {
    return data && other.data && data->shape.IsSame(other.data->shape);
}

// ---------------------------------------------------------------------------
// Face
// ---------------------------------------------------------------------------
Face::Face(const BSplineSurface<2> &surface, const std::array<Edge, 4> &sides,
           const std::array<bool, 4> &alongFrame)
    : Face(lift(surface), sides, alongFrame) {}

Face::Face(const BSplineSurface<3> &surface, const std::array<Edge, 4> &sides,
           const std::array<bool, 4> &alongFrame) {
    *this = makeFace(Access::surface(surface), sides, alongFrame);
}

BSplineSurface<3> Face::surface() const {
    if (!data) return BSplineSurface<3>();
    return Access::surface<3>(
        Handle(Geom_BSplineSurface)::DownCast(BRep_Tool::Surface(TopoDS::Face(data->shape))));
}

std::vector<Edge> Face::edges() const {
    std::vector<Edge> out;
    if (!data) return out;
    for (TopExp_Explorer it(data->shape, TopAbs_EDGE); it.More(); it.Next()) {
        out.push_back(Access::wrap<Edge>(it.Current().Oriented(TopAbs_FORWARD)));
    }
    return out;
}

Vec<3> Face::evaluate(double s, double t) const { return surface().evaluate(s, t); }

double Face::area() const {
    if (!data) return 0.0;
    return detail::guard("Face::area", [&] {
        GProp_GProps props;
        BRepGProp::SurfaceProperties(data->shape, props);
        return props.Mass();
    });
}

bool Face::isValid() const {
    if (!data) return false;
    return detail::guard("Face::isValid", [&] { return BRepCheck_Analyzer(data->shape).IsValid(); });
}

// ---------------------------------------------------------------------------
// Shape
// ---------------------------------------------------------------------------
// The faces are walked breadth first across their shared edges. The first face
// of each piece keeps its orientation; a neighbour that runs the shared edge the
// same way is reversed. An edge bounding more than two faces joins no two of
// them, and a piece that cannot be oriented consistently -- a Mobius band --
// keeps whatever orientation the walk reached it with.
Shape::Shape(const std::vector<Face> &faces, const std::vector<Edge> &looseEdges) {
    *this = detail::construct("Shape", [&] {
        BRep_Builder builder;
        TopTools_IndexedMapOfShape faceMap;
        TopoDS_Compound all;
        builder.MakeCompound(all);
        for (const Face &f : faces) {
            if (!f.isNull() && faceMap.Add(Access::shape(f)) == faceMap.Extent()) {
                builder.Add(all, Access::shape(f));
            }
        }
        const int n = faceMap.Extent();
        TopTools_IndexedDataMapOfShapeListOfShape edgeFaces;
        TopExp::MapShapesAndUniqueAncestors(all, TopAbs_EDGE, TopAbs_FACE, edgeFaces);

        std::vector<TopoDS_Shape> oriented(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) oriented[i] = faceMap(i + 1);
        std::vector<int> piece(static_cast<std::size_t>(n), -1);
        int pieces = 0;
        for (int seed = 0; seed < n; ++seed) {
            if (piece[seed] >= 0) continue;
            piece[seed] = pieces;
            std::deque<int> queue{seed};
            while (!queue.empty()) {
                const int i = queue.front();
                queue.pop_front();
                for (TopExp_Explorer it(oriented[i], TopAbs_EDGE); it.More(); it.Next()) {
                    const TopTools_ListOfShape &around = edgeFaces.FindFromKey(it.Current());
                    if (around.Extent() != 2) continue;
                    for (const TopoDS_Shape &other : around) {
                        const int j = faceMap.FindIndex(other) - 1;
                        if (j == i || piece[j] >= 0) continue;
                        if (orientationIn(oriented[j], it.Current()) == it.Current().Orientation()) {
                            oriented[j].Reverse();
                        }
                        piece[j] = pieces;
                        queue.push_back(j);
                    }
                }
            }
            ++pieces;
        }

        TopoDS_Compound compound;
        builder.MakeCompound(compound);
        for (int k = 0; k < pieces; ++k) {
            TopoDS_Shell shell;
            builder.MakeShell(shell);
            for (int i = 0; i < n; ++i) {
                if (piece[i] == k) builder.Add(shell, oriented[i]);
            }
            builder.Add(compound, shell);
        }
        for (const Edge &e : looseEdges) {
            if (!e.isNull()) builder.Add(compound, Access::shape(e));
        }
        return Access::wrap<Shape>(compound);
    });
}

Shape Shape::readBREP(const std::string &filename) {
    TopoDS_Shape shape;
    BRep_Builder builder;
    const bool ok = detail::guard("Shape::readBREP", [&] {
        return static_cast<bool>(BRepTools::Read(shape, filename.c_str(), builder));
    });
    if (!ok || shape.IsNull()) throw std::runtime_error("Shape::readBREP: could not read " + filename);
    return Access::wrap<Shape>(shape);
}

Shape Shape::readSTEP(const std::string &filename) {
    TopoDS_Shape shape = detail::guard("Shape::readSTEP", [&] {
        const QuietKernel quiet;
        STEPControl_Reader reader;
        if (reader.ReadFile(filename.c_str()) != IFSelect_RetDone) return TopoDS_Shape();
        reader.TransferRoots();
        return reader.OneShape();
    });
    if (shape.IsNull()) throw std::runtime_error("Shape::readSTEP: could not read " + filename);
    return Access::wrap<Shape>(shape);
}

std::vector<Face> Shape::faces() const {
    std::vector<Face> out;
    if (!data) return out;
    TopTools_IndexedMapOfShape map;
    TopExp::MapShapes(data->shape, TopAbs_FACE, map);
    for (int i = 1; i <= map.Extent(); ++i) out.push_back(Access::wrap<Face>(map(i)));
    return out;
}

std::vector<Edge> Shape::edges() const {
    std::vector<Edge> out;
    if (!data) return out;
    TopTools_IndexedMapOfShape map;
    TopExp::MapShapes(data->shape, TopAbs_EDGE, map);
    for (int i = 1; i <= map.Extent(); ++i) {
        out.push_back(Access::wrap<Edge>(map(i).Oriented(TopAbs_FORWARD)));
    }
    return out;
}

int Shape::faceCount() const { return data ? countShapes(data->shape, TopAbs_FACE) : 0; }
int Shape::edgeCount() const { return data ? countShapes(data->shape, TopAbs_EDGE) : 0; }
int Shape::vertexCount() const { return data ? countShapes(data->shape, TopAbs_VERTEX) : 0; }

int Shape::freeEdgeCount() const {
    if (!data) return 0;
    TopTools_IndexedDataMapOfShapeListOfShape map;
    TopExp::MapShapesAndUniqueAncestors(data->shape, TopAbs_EDGE, TopAbs_FACE, map);
    int count = 0;
    for (int i = 1; i <= map.Extent(); ++i) count += map(i).Extent() == 1;
    return count;
}

int Shape::sharedEdgeCount() const {
    if (!data) return 0;
    TopTools_IndexedDataMapOfShapeListOfShape map;
    TopExp::MapShapesAndUniqueAncestors(data->shape, TopAbs_EDGE, TopAbs_FACE, map);
    int count = 0;
    for (int i = 1; i <= map.Extent(); ++i) count += map(i).Extent() >= 2;
    return count;
}

double Shape::area() const {
    if (!data) return 0.0;
    return detail::guard("Shape::area", [&] {
        GProp_GProps props;
        BRepGProp::SurfaceProperties(data->shape, props);
        return props.Mass();
    });
}

bool Shape::isValid() const {
    if (!data) return false;
    return detail::guard("Shape::isValid", [&] { return BRepCheck_Analyzer(data->shape).IsValid(); });
}

bool Shape::writeBREP(const std::string &filename) const {
    if (!data) return false;
    try {
        return BRepTools::Write(data->shape, filename.c_str());
    } catch (const Standard_Failure &) {
        return false;
    }
}

bool Shape::writeSTEP(const std::string &filename) const {
    if (!data) return false;
    try {
        const QuietKernel quiet;
        STEPControl_Writer writer;
        if (writer.Transfer(data->shape, STEPControl_AsIs) != IFSelect_RetDone) return false;
        return writer.Write(filename.c_str()) == IFSelect_RetDone;
    } catch (const Standard_Failure &) {
        return false;
    }
}

} // namespace geom
