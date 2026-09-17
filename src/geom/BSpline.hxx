#ifndef __GEOM_BSPLINE_HXX__
#define __GEOM_BSPLINE_HXX__

#include <vector>

#include "geom/Vec.hxx"

// B-spline curves and tensor-product surfaces, of any degree, over any
// non-decreasing knot vector, in two or three dimensions.
//
// This is the spline layer MERIDIAN's Stage 9 (SplineFit) and TORSION's reuse
// of it were built on, taken out of the pipeline and given no opinion about
// arrangements, arcs or patches. The pieces, and where each one lives:
//
//     geom/BSpline.hxx    knot vectors, the basis, BSplineCurve, BSplineSurface
//     geom/Fitting.hxx    least-squares fitting and interpolation of point runs
//     geom/Polyline.hxx   a polyline parameterised by normalised chord length
//     geom/ArcLength.hxx  arc length of any parametric curve, and its inverse
//     geom/Coons.hxx      Coons blends, of points and of control nets
//
// The algorithms are the ones in Piegl and Tiller, *The NURBS Book* (2nd ed.):
// FindSpan and BasisFuns (A2.1, A2.2), the curve point (A3.1), the derivative
// curve (Eq. 3.8), knot insertion (A5.1), global interpolation (A9.1, with the
// averaged knots of Eq. 9.8) and least-squares approximation (Sec. 9.4.1).
//
// ### Conventions
//
// A curve with n control points and degree p has n + p + 1 knots, and is
// defined over the domain [knots[p], knots[n]]. evaluate() clamps its argument
// into that domain rather than extrapolating. Nothing requires the knot vector
// to be clamped, but every constructor in this directory makes clamped ones,
// and on a clamped knot vector the first and last control points are the two
// ends of the curve exactly -- which is what lets two splines that share an end
// point meet at the same double.
//
// A surface's control net is stored row-major, net[j * countS() + i] being the
// control point at (s_i, t_j): i runs along s, j along t. That is the layout
// SplineFit::Patch has always used, and the one Stage 10 samples.
//
// ### Numerical exactness
//
// The arithmetic in the .cxx is written expression for expression as it was in
// SplineFit, not merely to the same formulas: the pipeline's watertightness
// checks compare doubles for exact equality, and on this toolchain a reordered
// sum is a different double. Before rewriting anything in the evaluation or
// fitting paths, rerun the MERIDIAN/TORSION corpus and diff the control points.
namespace geom {

// The basis is evaluated into stack buffers of this size.
constexpr int kMaxDegree = 15;

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
// them.
std::vector<double> grevilleAbscissae(const std::vector<double> &knots, int degree);

// The number of control points a knot vector of this degree carries.
inline int controlPointCount(const std::vector<double> &knots, int degree) {
    return static_cast<int>(knots.size()) - degree - 1;
}

// The knot span containing u: the index k with knots[k] <= u < knots[k+1],
// restricted to [degree, n - 1] so that u at the end of the domain lands in the
// last non-empty span.
int findSpan(const std::vector<double> &knots, int degree, double u);

// The degree + 1 basis functions that are non-zero on `span`, at u, into N.
void basisFunctions(const std::vector<double> &knots, int degree, int span, double u, double *N);

// ---------------------------------------------------------------------------
// BSplineCurve
// ---------------------------------------------------------------------------
template <std::size_t D>
class BSplineCurve {
public:
    BSplineCurve() = default;

    // Throws std::invalid_argument unless 0 <= degree <= kMaxDegree, the knots
    // do not decrease, the domain is not empty, and
    // knots.size() == ctrl.size() + degree + 1.
    BSplineCurve(int degree, std::vector<double> knots, std::vector<Vec<D>> ctrl);

    bool empty() const { return ctrl.empty(); }
    int size() const { return static_cast<int>(ctrl.size()); }
    int degree() const { return p; }
    const std::vector<double> &knots() const { return knot; }
    const std::vector<Vec<D>> &controlPoints() const { return ctrl; }

    // Moving a control point cannot invalidate the curve; adding or removing
    // one would, so only the former is offered.
    Vec<D> &controlPoint(int i) { return ctrl[i]; }
    const Vec<D> &controlPoint(int i) const { return ctrl[i]; }

    double domainBegin() const;
    double domainEnd() const;

    // The point at u, clamped into the domain. The zero vector on an empty curve.
    Vec<D> evaluate(double u) const;
    Vec<D> operator()(double u) const { return evaluate(u); }

    // segments + 1 points at uniform steps of the parameter over the domain.
    std::vector<Vec<D>> sample(int segments) const;

    // The derivative as a curve in its own right, one degree lower over the
    // same knots less the first and last. Evaluate that for tangents rather
    // than calling this per point. Throws on a degree-0 curve.
    BSplineCurve derivative() const;

    // The same point set traversed the other way: control points reversed and
    // knots mirrored across the domain, so reversed().evaluate(a + b - u) is
    // evaluate(u). Mirroring is exact only if the knots are, which uniform
    // thirds are not -- see SplineFit::Report::maxBoundaryGap.
    BSplineCurve reversed() const;

    // Boehm insertion of u, `times` times, capped so that no knot exceeds
    // multiplicity `degree`. The curve does not change, geometrically or in
    // parameterisation; it gains control points. u must lie strictly inside the
    // domain. Returns how many copies were inserted.
    int insertKnot(double u, int times = 1);

private:
    int p = 0;
    std::vector<double> knot;
    std::vector<Vec<D>> ctrl;
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
    BSplineSurface(int degreeS, std::vector<double> knotsS,
                   int degreeT, std::vector<double> knotsT,
                   std::vector<Vec<D>> net);

    bool empty() const { return net.empty(); }
    int degreeS() const { return pS; }
    int degreeT() const { return pT; }
    int countS() const { return nS; }
    int countT() const { return nT; }
    const std::vector<double> &knotsS() const { return knotS; }
    const std::vector<double> &knotsT() const { return knotT; }
    const std::vector<Vec<D>> &controlNet() const { return net; }

    const Vec<D> &control(int i, int j) const { return net[static_cast<std::size_t>(j) * nS + i]; }
    Vec<D> &control(int i, int j) { return net[static_cast<std::size_t>(j) * nS + i]; }

    // The point at (s, t), each clamped into its domain.
    Vec<D> evaluate(double s, double t) const;
    Vec<D> operator()(double s, double t) const { return evaluate(s, t); }

    // Row j of the net as a curve along s, and column i as a curve along t. On
    // clamped knots the first and last of each are the surface's four
    // boundary curves exactly.
    BSplineCurve<D> row(int j) const;
    BSplineCurve<D> column(int i) const;

private:
    int pS = 0, pT = 0;
    int nS = 0, nT = 0;
    std::vector<double> knotS, knotT;
    std::vector<Vec<D>> net;
};

extern template class BSplineCurve<2>;
extern template class BSplineCurve<3>;
extern template class BSplineSurface<2>;
extern template class BSplineSurface<3>;

using BSplineCurve2 = BSplineCurve<2>;
using BSplineCurve3 = BSplineCurve<3>;
using BSplineSurface2 = BSplineSurface<2>;
using BSplineSurface3 = BSplineSurface<3>;

} // namespace geom

#endif // __GEOM_BSPLINE_HXX__
