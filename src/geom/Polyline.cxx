#include "geom/Polyline.hxx"

#include <algorithm>
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
}

template <std::size_t D>
Vec<D> Polyline<D>::evaluate(double u) const {
    if (pts.empty()) return zero<D>();
    if (pts.size() == 1) return pts.front();
    u = std::max(0.0, std::min(1.0, u));
    if (u <= 0.0) return pts.front();
    if (u >= 1.0) return pts.back();
    const std::size_t k = static_cast<std::size_t>(
        std::lower_bound(param.begin(), param.end(), u) - param.begin());
    if (k == 0) return pts.front();
    if (k >= pts.size()) return pts.back();
    const double seg = param[k] - param[k - 1];
    const double w = seg > 0.0 ? (u - param[k - 1]) / seg : 0.0;
    return pts[k - 1] + (pts[k] - pts[k - 1]) * w;
}

template <std::size_t D>
double Polyline<D>::distance(const Vec<D> &p) const {
    if (pts.empty()) return std::numeric_limits<double>::infinity();
    if (pts.size() == 1) return norm(p - pts.front());
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i + 1 < pts.size(); ++i) {
        best = std::min(best, distanceToSegment(p, pts[i], pts[i + 1]));
    }
    return best;
}

template class Polyline<2>;
template class Polyline<3>;
template double polylineLength<2>(const std::vector<Vec<2>> &);
template double polylineLength<3>(const std::vector<Vec<3>> &);

} // namespace geom
