#include "geom/Polyline.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace geom {

template <std::size_t D>
double polylineLength(const std::vector<Vec<D>> &points) {
    double total = 0.0;
    for (std::size_t i = 1; i < points.size(); ++i) total += norm(points[i] - points[i - 1]);
    return total;
}

template <std::size_t D>
Polyline<D>::Polyline(std::vector<Vec<D>> points) : pts(std::move(points)) {
    param.assign(pts.size(), 0.0);
    if (pts.size() < 2) return;
    for (std::size_t i = 1; i < pts.size(); ++i) {
        param[i] = param[i - 1] + norm(pts[i] - pts[i - 1]);
    }
    total = param.back();
    if (total > 0.0) {
        for (double &x : param) x /= total;
    } else {
        for (std::size_t i = 0; i < param.size(); ++i) {
            param[i] = static_cast<double>(i) / (param.size() - 1.0);
        }
    }
    param.back() = 1.0;

    // The vertices the kernel gets: each one further than an ulp in parameter
    // from the last one kept. The final vertex always stays, displacing the
    // one before it if the two are that close.
    std::vector<double> knots{0.0, 0.0};
    std::vector<Vec<D>> ctrl{pts.front()};
    for (std::size_t i = 1; i < pts.size(); ++i) {
        const double last = knots.back();
        const bool apart = param[i] - last > std::nextafter(last, 2.0) - last;
        if (i + 1 == pts.size()) {
            if (!apart && ctrl.size() > 1) {
                knots.pop_back();
                ctrl.pop_back();
            }
            knots.push_back(1.0);
            ctrl.push_back(pts[i]);
        } else if (apart) {
            knots.push_back(param[i]);
            ctrl.push_back(pts[i]);
        }
    }
    knots.push_back(1.0);
    spline = BSplineCurve<D>(1, knots, ctrl);
}

template <std::size_t D>
Vec<D> Polyline<D>::evaluate(double u) const {
    if (pts.empty()) return zero<D>();
    if (pts.size() == 1) return pts.front();
    if (!(u > 0.0)) return pts.front();
    if (u >= 1.0) return pts.back();
    return spline.evaluate(u);
}

template <std::size_t D>
double Polyline<D>::distance(const Vec<D> &p) const {
    if (pts.empty()) return std::numeric_limits<double>::infinity();
    if (pts.size() == 1) return norm(p - pts.front());
    return spline.distance(p);
}

template class Polyline<2>;
template class Polyline<3>;
template double polylineLength<2>(const std::vector<Vec<2>> &);
template double polylineLength<3>(const std::vector<Vec<3>> &);

} // namespace geom
