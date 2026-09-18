#ifndef __GEOM_BSPLINE_HXX__
#define __GEOM_BSPLINE_HXX__

#include <array>
#include <memory>
#include <vector>

#include "geom/Vec.hxx"

// B-spline curves and tensor-product surfaces, of any degree from 1 to
// kMaxDegree, over any knot vector the kernel accepts, in two or three
// dimensions.
//
// ### The kernel
//
// Every curve and surface in the repository is an OpenCASCADE object, and
// src/geom is the only code that knows it. The classes here hold a
// Geom_BSplineCurve or a Geom_BSplineSurface and route evaluation, derivatives,
// knot insertion, reversal, arc length and projection through the kernel; the
// headers name nothing of OpenCASCADE, and the build puts its include
// directories on src/geom's sources alone, so nothing outside this directory
// can call it even by accident. (ctest's Geom_OpenCascadeConfinedToGeom checks
// the same thing from the source side.) A new geometric need is met by adding
// to src/geom, not by reaching past it.
//
// Two-dimensional geometry is carried in three dimensions, in the plane z = 0:
// the kernel's topology (Topology.hxx) is three-dimensional, and a 2-D curve
// that is already a 3-D object can become the edge of a face, be revolved or be
// written to STEP without a conversion. The kernel evaluates each coordinate on
// its own, so reading a 2-D point back loses nothing.
//
// The pieces, and where each one lives:
//
//     geom/BSpline.hxx    knot vectors, the basis, BSplineCurve, BSplineSurface
//     geom/Fitting.hxx    least-squares fitting and interpolation of point runs
//     geom/Polyline.hxx   a polyline as a degree-1 curve on normalised chord length
//     geom/ArcLength.hxx  a tabulated arc length of any parametric curve
//     geom/Coons.hxx      Coons blends, of points, control nets and curves
//     geom/Topology.hxx   Vertex, Edge, Face and Shape: the B-rep, and STEP/BREP out
//
// ### What is still computed here, and why
//
// Only what the kernel does not offer in the form the pipelines need. The
// least-squares fit of Fitting.hxx holds its ends, shares one knot vector across
// every arc and regularises towards the chord, and none of OpenCASCADE's
// approximators does all three -- they choose their own knots, which would give
// opposite sides of a patch different ones. The Coons net of Coons.hxx is taken
// at the Greville abscissae; GeomFill's Coons filling takes it at uniform pole
// indices instead, which is not linearly precise and bows a patch on straight
// sides by 1.5e-2 of its size. Both constructions still evaluate the basis and
// store their results through the kernel.
//
// ### Conventions
//
// A curve with n control points and degree p has n + p + 1 knots, and is
// defined over the domain [knots[p], knots[n]]. evaluate() clamps its argument
// into that domain rather than extrapolating. The kernel requires interior
// knots of multiplicity at most p -- a curve is at least C0 -- and distinct
// knots more than an ulp apart. Every constructor in this directory makes
// clamped knot vectors, and on a clamped one the first and last control points
// are the two ends of the curve exactly, which is what lets two splines that
// share an end point meet at the same double.
//
// A surface's control net is stored row-major, net[j * countS() + i] being the
// control point at (s_i, t_j): i runs along s, j along t. That is the layout
// SplineFit::Patch has always used, and the one Stage 10 samples.
//
// The objects are values. Copies share one kernel object, which is never
// modified in place: the operations that change a curve (insertKnot,
// setControlPoint) give it a fresh one first, so a copy, or an Edge built on
// the curve earlier, keeps the geometry it had.
//
// ### Numerical exactness
//
// The kernel's de Boor evaluation is not the basis-function sum SplineFit used
// before it, and the two differ in the last bit or two (1.2e-15 at most on a
// unit cubic). What the pipelines compare exactly -- control points two patches
// share, corners against nodes -- is copied rather than computed, and stays
// exact; what they sample is to rounding.
namespace geom {

namespace detail {
struct CurveData;
struct SurfaceData;
struct Access;
} // namespace detail

// The kernel's own limit (Geom_BSplineCurve::MaxDegree()).
constexpr int kMaxDegree = 25;

// ---------------------------------------------------------------------------
// Knot vectors and the basis
// ---------------------------------------------------------------------------

// A clamped uniform knot vector on [0, 1] with `segments` polynomial pieces:
// degree + 1 zeros, the interior knots i / segments, and degree + 1 ones. It
// carries segments + degree control points.
std::vector<double> clampedUniformKnots(int degree, int segments);

// The clamped knot vector on [0, 1] that global interpolation of points at the
// given parameters needs: each interior knot is the mean of `degree`
// consecutive parameters (The NURBS Book, Eq. 9.8). `params` must start at 0
// and end at 1; it carries params.size() control points. Throws unless
// degree >= 1 and there are at least degree + 1 parameters.
std::vector<double> averagedKnots(const std::vector<double> &params, int degree);

// The Greville abscissae, (knots[i+1] + ... + knots[i+degree]) / degree, one
// per control point: the parameter each control point "belongs to". A curve
// whose control points sit on a line at their Greville abscissae is that line,
// linearly parameterised, which is why the Coons blend of Coons.hxx is taken at
// them. degree >= 1.
std::vector<double> grevilleAbscissae(const std::vector<double> &knots, int degree);

// The number of control points a knot vector of this degree carries.
inline int controlPointCount(const std::vector<double> &knots, int degree) {
    return static_cast<int>(knots.size()) - degree - 1;
}

// The degree + 1 basis functions that can be non-zero at u, into N, from the
// kernel (BSplCLib). Returns the index of the control point N[0] belongs to, so
// that N[a] is the weight of control point (returned + a). u is clamped into
// the domain; at its end the last basis function is 1. degree 0 is allowed here,
// though no curve can have it.
int basisFunctions(const std::vector<double> &knots, int degree, double u, double *N);

// ---------------------------------------------------------------------------
// BSplineCurve
// ---------------------------------------------------------------------------
template <std::size_t D>
class BSplineCurve {
public:
    BSplineCurve() = default;

    // Throws std::invalid_argument unless 1 <= degree <= kMaxDegree, the knots
    // do not decrease, knots.size() == ctrl.size() + degree + 1, the domain is
    // not empty, and the kernel accepts the knot vector (no interior knot of
    // multiplicity above the degree, no two distinct knots within an ulp).
    BSplineCurve(int degree, const std::vector<double> &knots, const std::vector<Vec<D>> &ctrl);

    bool empty() const { return !data; }
    int size() const;
    int degree() const;
    // The flat knot vector, as the constructor took it.
    std::vector<double> knots() const;
    std::vector<Vec<D>> controlPoints() const;
    Vec<D> controlPoint(int i) const;
    // Moving a control point cannot invalidate the curve; adding or removing
    // one would, so only the former is offered.
    void setControlPoint(int i, const Vec<D> &p);

    double domainBegin() const;
    double domainEnd() const;

    // The point at u, clamped into the domain. The zero vector on an empty curve.
    Vec<D> evaluate(double u) const;
    Vec<D> operator()(double u) const { return evaluate(u); }

    // The order-th derivative with respect to u at u (clamped), order >= 1. Zero
    // above the degree.
    Vec<D> derivativeAt(double u, int order = 1) const;

    // segments + 1 points at uniform steps of the parameter over the domain.
    std::vector<Vec<D>> sample(int segments) const;

    // The derivative as a curve in its own right, one degree lower over the
    // same knots less the first and last. Evaluate that for tangents along a
    // whole curve rather than calling this per point. Throws on a curve of
    // degree 1, whose derivative is piecewise constant, and on one with an
    // interior knot of full multiplicity, whose derivative jumps there -- the
    // kernel can represent neither; derivativeAt() still answers for both.
    BSplineCurve derivative() const;

    // The same point set traversed the other way, over the mirrored knots, so
    // reversed().evaluate(a + b - u) is evaluate(u). Mirroring is exact only if
    // the knots are, which uniform thirds are not: 1 - 2/3 and 1/3 are an ulp
    // apart.
    BSplineCurve reversed() const;

    // Knot insertion of u, `times` times, capped so that no knot exceeds
    // multiplicity `degree`. The curve does not change, geometrically or in
    // parameterisation; it gains control points. u must lie strictly inside the
    // domain. Returns how many copies were inserted.
    int insertKnot(double u, int times = 1);

    // Arc length over the whole domain, or between two parameters, by the
    // kernel's Gauss integration.
    double length() const;
    double length(double u0, double u1) const;
    // The parameter at which arc length s has been travelled from the start of
    // the domain, clamped to it.
    double parameterAtLength(double s) const;

    // The nearest point of the curve to p, over its whole domain -- the ends and
    // any corner (a knot of full multiplicity) included, which the kernel's
    // projection by itself does not look at.
    struct Nearest {
        double parameter = 0.0;
        Vec<D> point{};
        double distance = 0.0;
    };
    Nearest nearest(const Vec<D> &p) const;
    double distance(const Vec<D> &p) const { return nearest(p).distance; }

private:
    std::shared_ptr<const detail::CurveData> data;
    friend struct detail::Access;
};

// ---------------------------------------------------------------------------
// BSplineSurface
// ---------------------------------------------------------------------------
template <std::size_t D>
class BSplineSurface {
public:
    BSplineSurface() = default;

    // `net` is row-major, net[j * countS + i] at (s_i, t_j), with countS and
    // countT fixed by the knot vectors. Throws std::invalid_argument on the
    // same conditions as BSplineCurve, in each direction, or if the net is the
    // wrong size.
    BSplineSurface(int degreeS, const std::vector<double> &knotsS,
                   int degreeT, const std::vector<double> &knotsT,
                   const std::vector<Vec<D>> &net);

    bool empty() const { return !data; }
    int degreeS() const;
    int degreeT() const;
    int countS() const;
    int countT() const;
    std::vector<double> knotsS() const;
    std::vector<double> knotsT() const;
    std::vector<Vec<D>> controlNet() const;

    Vec<D> control(int i, int j) const;
    void setControl(int i, int j, const Vec<D> &p);

    // The point at (s, t), each clamped into its domain.
    Vec<D> evaluate(double s, double t) const;
    Vec<D> operator()(double s, double t) const { return evaluate(s, t); }
    // The two first partial derivatives there, d/ds and d/dt.
    std::array<Vec<D>, 2> partials(double s, double t) const;

    // Row j of the net as a curve along s, and column i as a curve along t. On
    // clamped knots the first and last of each are the surface's four
    // boundary curves exactly.
    BSplineCurve<D> row(int j) const;
    BSplineCurve<D> column(int i) const;

private:
    std::shared_ptr<const detail::SurfaceData> data;
    friend struct detail::Access;
};

extern template class BSplineCurve<2>;
extern template class BSplineCurve<3>;
extern template class BSplineSurface<2>;
extern template class BSplineSurface<3>;

using BSplineCurve2 = BSplineCurve<2>;
using BSplineCurve3 = BSplineCurve<3>;
using BSplineSurface2 = BSplineSurface<2>;
using BSplineSurface3 = BSplineSurface<3>;

// The same curve or surface in the other dimension: a 2-D one placed in the
// plane z = 0, or a 3-D one with z dropped (it is not checked that z was zero).
// The two share one kernel object.
BSplineCurve<3> lift(const BSplineCurve<2> &c);
BSplineCurve<2> flatten(const BSplineCurve<3> &c);
BSplineSurface<3> lift(const BSplineSurface<2> &s);
BSplineSurface<2> flatten(const BSplineSurface<3> &s);

} // namespace geom

#endif // __GEOM_BSPLINE_HXX__
