#ifndef __GEOM_DETAIL_OCC_HXX__
#define __GEOM_DETAIL_OCC_HXX__

// The bridge between src/geom's public types and OpenCASCADE.
//
// Only the .cxx files of src/geom include this. The public headers name
// nothing of OpenCASCADE -- each class carries a pointer to one of the structs
// below, declared but not defined there -- and the OpenCASCADE include
// directories are a PRIVATE property of the CrossGenGeom target, so a file
// anywhere else in the repository that tried to include this, or any
// OpenCASCADE header, would not compile. That is the point: the kernel is
// src/geom's business and nobody else's.

#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Standard_Failure.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "geom/BSpline.hxx"
#include "geom/Topology.hxx"
#include "geom/Vec.hxx"

namespace geom {
namespace detail {

// Every curve and surface is held as a 3-D OpenCASCADE object. A 2-D one lies
// in the plane z = 0 and is read back by dropping z, which loses nothing: the
// kernel evaluates each coordinate independently, so x and y come out as the
// same doubles a Geom2d object would have given.
struct CurveData {
    Handle(Geom_BSplineCurve) curve;
};

struct SurfaceData {
    Handle(Geom_BSplineSurface) surface;
};

// A vertex, an edge, a face or a compound. The TopoDS_Shape is a handle to
// shared topology, so two geom objects that hold the same one *are* the same
// sub-shape to the kernel -- which is how edges come to be shared by faces.
struct ShapeData {
    TopoDS_Shape shape;
};

inline gp_Pnt toPnt(const Vec<2> &p) { return gp_Pnt(p[0], p[1], 0.0); }
inline gp_Pnt toPnt(const Vec<3> &p) { return gp_Pnt(p[0], p[1], p[2]); }

template <std::size_t D>
inline Vec<D> fromPnt(const gp_Pnt &p) {
    Vec<D> out = zero<D>();
    out[0] = p.X();
    out[1] = p.Y();
    if (D > 2) out[2] = p.Z();
    return out;
}

template <std::size_t D>
inline Vec<D> fromVec(const gp_Vec &v) {
    Vec<D> out = zero<D>();
    out[0] = v.X();
    out[1] = v.Y();
    if (D > 2) out[2] = v.Z();
    return out;
}

// A flat knot vector (the form every public signature uses) as the distinct
// knots and multiplicities OpenCASCADE's constructors take. Equal means equal
// as doubles. The vector must not be empty.
struct SplitKnots {
    TColStd_Array1OfReal knots;
    TColStd_Array1OfInteger mults;
};

inline SplitKnots splitKnots(const std::vector<double> &flat) {
    std::vector<double> k;
    std::vector<int> m;
    for (const double x : flat) {
        if (!k.empty() && x == k.back()) {
            ++m.back();
        } else {
            k.push_back(x);
            m.push_back(1);
        }
    }
    const int n = static_cast<int>(k.size());
    SplitKnots out{TColStd_Array1OfReal(1, n), TColStd_Array1OfInteger(1, n)};
    for (int i = 0; i < n; ++i) {
        out.knots(i + 1) = k[i];
        out.mults(i + 1) = m[i];
    }
    return out;
}

inline std::vector<double> toVector(const TColStd_Array1OfReal &a) {
    std::vector<double> out;
    out.reserve(static_cast<std::size_t>(a.Length()));
    for (int i = a.Lower(); i <= a.Upper(); ++i) out.push_back(a(i));
    return out;
}

inline TColStd_Array1OfReal toArray(const std::vector<double> &v) {
    TColStd_Array1OfReal a(1, static_cast<int>(v.size()));
    for (std::size_t i = 0; i < v.size(); ++i) a(static_cast<int>(i) + 1) = v[i];
    return a;
}

inline std::string describe(const Standard_Failure &e) {
    const char *msg = e.GetMessageString();
    return std::string("OpenCASCADE: ") + (msg && *msg ? msg : e.DynamicType()->Name());
}

// Runs a kernel call and turns an OpenCASCADE exception into a standard one,
// so that nothing of OpenCASCADE, its exception types included, crosses the
// boundary of src/geom. guard() is for operations on objects already known to
// be valid, where a failure is the kernel's; construct() is for building one
// from caller-supplied data, where it is the caller's.
template <class F>
auto guard(const char *what, F &&f) -> decltype(f()) {
    try {
        return f();
    } catch (const Standard_Failure &e) {
        throw std::runtime_error(std::string(what) + ": " + describe(e));
    }
}

template <class F>
auto construct(const char *what, F &&f) -> decltype(f()) {
    try {
        return f();
    } catch (const Standard_Failure &e) {
        throw std::invalid_argument(std::string(what) + ": " + describe(e));
    }
}

// Reaches into the public classes, which name this struct as a friend.
struct Access {
    template <std::size_t D>
    static Handle(Geom_BSplineCurve) curve(const BSplineCurve<D> &c) {
        return c.data ? c.data->curve : Handle(Geom_BSplineCurve)();
    }
    template <std::size_t D>
    static BSplineCurve<D> curve(Handle(Geom_BSplineCurve) h) {
        BSplineCurve<D> c;
        if (!h.IsNull()) c.data = std::make_shared<const CurveData>(CurveData{std::move(h)});
        return c;
    }

    template <std::size_t D>
    static Handle(Geom_BSplineSurface) surface(const BSplineSurface<D> &s) {
        return s.data ? s.data->surface : Handle(Geom_BSplineSurface)();
    }
    template <std::size_t D>
    static BSplineSurface<D> surface(Handle(Geom_BSplineSurface) h) {
        BSplineSurface<D> s;
        if (!h.IsNull()) s.data = std::make_shared<const SurfaceData>(SurfaceData{std::move(h)});
        return s;
    }

    template <class T>
    static TopoDS_Shape shape(const T &t) {
        return t.data ? t.data->shape : TopoDS_Shape();
    }
    template <class T>
    static T wrap(TopoDS_Shape s) {
        T t;
        if (!s.IsNull()) t.data = std::make_shared<const ShapeData>(ShapeData{std::move(s)});
        return t;
    }
};

// A private copy of a curve or surface, for the operations that modify one in
// place. The public classes share their kernel object between copies and never
// change it, so anything that would has to take its own first.
inline Handle(Geom_BSplineCurve) copyOf(const Handle(Geom_BSplineCurve) &c) {
    return Handle(Geom_BSplineCurve)::DownCast(c->Copy());
}

inline Handle(Geom_BSplineSurface) copyOf(const Handle(Geom_BSplineSurface) &s) {
    return Handle(Geom_BSplineSurface)::DownCast(s->Copy());
}

} // namespace detail
} // namespace geom

#endif // __GEOM_DETAIL_OCC_HXX__
