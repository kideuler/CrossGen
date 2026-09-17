#ifndef __GEOM_ARC_LENGTH_HXX__
#define __GEOM_ARC_LENGTH_HXX__

#include <vector>

#include "geom/Vec.hxx"

namespace geom {

// Arc length along a parametric curve on [0, 1], tabulated at uniform
// parameter, and its inverse.
//
// A spline's parameter is not uniform in length -- chord-length fitting makes
// it close, never equal -- so placing nodes "evenly along a curve" means
// travelling equal lengths and asking which parameter got there. The table
// answers that by linear interpolation between samples, which is what MERIDIAN's
// Stage 10 does for every arc before it puts a node down.
//
// The curve is anything callable as curve(double) returning a Vec<D>: a
// BSplineCurve, a Polyline, or a lambda over something else entirely.
class ArcLengthTable {
public:
    ArcLengthTable() = default;

    // Samples the curve at i / samples for i = 0..samples; at least one step.
    template <class Curve>
    ArcLengthTable(const Curve &curve, int samples) {
        if (samples < 1) samples = 1;
        cum.assign(static_cast<std::size_t>(samples) + 1, 0.0);
        auto prev = curve(0.0);
        for (int i = 1; i <= samples; ++i) {
            const auto p = curve(static_cast<double>(i) / samples);
            cum[i] = cum[i - 1] + geom::norm(p - prev);
            prev = p;
        }
        total = cum.back();
    }

    bool empty() const { return cum.empty(); }
    double length() const { return total; }
    int samples() const { return static_cast<int>(cum.size()) - 1; }
    // Length travelled at each sample, cumulative()[i] at u = i / samples().
    const std::vector<double> &cumulative() const { return cum; }

    // The parameter in [0, 1] at which length s has been travelled, clamped at
    // both ends. Zero on a curve of no length.
    double parameterAt(double s) const;
    // The length travelled by parameter u, clamped into [0, 1].
    double lengthAt(double u) const;

private:
    std::vector<double> cum;
    double total = 0.0;
};

} // namespace geom

#endif // __GEOM_ARC_LENGTH_HXX__
