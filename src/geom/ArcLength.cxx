#include "geom/ArcLength.hxx"

#include <algorithm>
#include <cmath>

namespace geom {

double ArcLengthTable::parameterAt(double s) const {
    const int m = static_cast<int>(cum.size()) - 1;
    if (m <= 0 || !(total > 0.0)) return 0.0;
    if (s <= 0.0) return 0.0;
    if (s >= total) return 1.0;
    const int k = static_cast<int>(std::lower_bound(cum.begin(), cum.end(), s) - cum.begin());
    if (k <= 0) return 0.0;
    if (k > m) return 1.0;
    const double seg = cum[k] - cum[k - 1];
    const double w = seg > 0.0 ? (s - cum[k - 1]) / seg : 0.0;
    return (static_cast<double>(k - 1) + w) / m;
}

double ArcLengthTable::lengthAt(double u) const {
    const int m = static_cast<int>(cum.size()) - 1;
    if (m <= 0) return 0.0;
    u = std::max(0.0, std::min(1.0, u));
    const double x = u * m;
    const int k = std::min(m - 1, static_cast<int>(std::floor(x)));
    const double w = x - k;
    return cum[k] + (cum[k + 1] - cum[k]) * w;
}

} // namespace geom
