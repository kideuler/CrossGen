#ifndef __GEOM_POLYLINE_HXX__
#define __GEOM_POLYLINE_HXX__

#include <vector>

#include "geom/BSpline.hxx"
#include "geom/Vec.hxx"

namespace geom {

// The cumulative length along a run of points.
template <std::size_t D>
double polylineLength(const std::vector<Vec<D>> &points);

// A polyline as a curve on [0, 1], parameterised by normalised chord length:
// vertex i sits at parameters()[i], the first at 0 and the last at 1 exactly.
//
// This is how SplineFit carries an arc it must not approximate -- a run of
// edges of dS or of a material interface, which already *is* the geometry --
// and it is the parameterisation fitCurve() fits in by default, so a polyline
// and the spline fitted through it can be compared at the same u. Where the
// polyline has no length at all its vertices are spread uniformly instead.
//
// To the kernel it is a degree-1 B-spline whose knots are those parameters
// (curve()), so it can be the edge of a face and a side of a Coons surface like
// any other curve. The kernel wants its knots more than an ulp apart, so a
// vertex that repeats its predecessor -- or lies within an ulp of it in
// parameter -- is left out of the curve; that changes neither the point set nor
// the parameterisation. points() and parameters() still give every vertex.
template <std::size_t D>
class Polyline {
public:
    Polyline() = default;
    explicit Polyline(std::vector<Vec<D>> points);

    bool empty() const { return pts.empty(); }
    std::size_t size() const { return pts.size(); }
    const std::vector<Vec<D>> &points() const { return pts; }
    const std::vector<double> &parameters() const { return param; }
    double length() const { return total; }

    // The degree-1 curve. Empty unless there are at least two vertices.
    const BSplineCurve<D> &curve() const { return spline; }

    // The point at u, clamped into [0, 1]: the first and last vertex exactly at
    // the ends, the kernel's evaluation of curve() between them. The zero
    // vector on an empty polyline.
    Vec<D> evaluate(double u) const;
    Vec<D> operator()(double u) const { return evaluate(u); }

    // The distance from p to the nearest point of the polyline.
    double distance(const Vec<D> &p) const;

private:
    std::vector<Vec<D>> pts;
    std::vector<double> param;
    double total = 0.0;
    BSplineCurve<D> spline;
};

extern template class Polyline<2>;
extern template class Polyline<3>;
extern template double polylineLength<2>(const std::vector<Vec<2>> &);
extern template double polylineLength<3>(const std::vector<Vec<3>> &);

using Polyline2 = Polyline<2>;
using Polyline3 = Polyline<3>;

} // namespace geom

#endif // __GEOM_POLYLINE_HXX__
