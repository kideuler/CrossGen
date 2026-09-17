#ifndef __GEOM_POLYLINE_HXX__
#define __GEOM_POLYLINE_HXX__

#include <vector>

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

    // The point at u, clamped into [0, 1], interpolated linearly within the
    // segment that contains it. The zero vector on an empty polyline.
    Vec<D> evaluate(double u) const;
    Vec<D> operator()(double u) const { return evaluate(u); }

    // The distance from p to the nearest point of the polyline.
    double distance(const Vec<D> &p) const;

private:
    std::vector<Vec<D>> pts;
    std::vector<double> param;
    double total = 0.0;
};

extern template class Polyline<2>;
extern template class Polyline<3>;
extern template double polylineLength<2>(const std::vector<Vec<2>> &);
extern template double polylineLength<3>(const std::vector<Vec<3>> &);

using Polyline2 = Polyline<2>;
using Polyline3 = Polyline<3>;

} // namespace geom

#endif // __GEOM_POLYLINE_HXX__
