#ifndef __GEOM_FITTING_HXX__
#define __GEOM_FITTING_HXX__

#include <vector>

#include "geom/BSpline.hxx"
#include "geom/Vec.hxx"

// Turning a run of points into a B-spline curve: by least squares, when the
// points are data the curve should stay near (fitCurve), and by interpolation,
// when they are points it must pass through (interpolateCurve).
//
// fitCurve() is Sec. 5 of Shepherd, Gu and Hughes (2022) -- "chord-length
// parameterise, least-squares fit a cubic B-spline with a fixed uniform knot
// vector" -- with the three things SplineFit found it needed on real traced
// arcs, each of them now an option:
//
//   * the ends held rather than fitted (FitOptions::pinEnds), so that two
//     curves meeting at a point arrive at that point to the last bit instead
//     of each missing it by its own residual;
//   * a short run subdivided before it is fitted (FitOptions::densify), since
//     an arc a few triangles long carries fewer points than the fit has
//     unknowns -- and a point interpolated along a segment of a polyline lies
//     on that polyline exactly, so this invents no data;
//   * a Tikhonov pull towards the straight line between the ends
//     (FitOptions::regularisation), which binds only where the data still
//     cannot determine a control point, and without which such a control point
//     goes anywhere at all.
//
// None of OpenCASCADE's approximators does those three together -- they pick
// their own knots, and a Coons patch needs opposite sides on the same ones --
// so the normal equations are assembled and solved here. The basis values in
// them are the kernel's, and so is the curve that comes out.
namespace geom {

enum class Parameterization {
    Uniform,      // i / (n - 1)
    ChordLength,  // cumulative distance, normalised
    Centripetal,  // cumulative square root of distance, normalised (Lee 1989)
};

// Parameters in [0, 1] for a run of points, one per point, starting at 0 and
// non-decreasing. A run with no length under the chosen measure falls back to
// Uniform.
template <std::size_t D>
std::vector<double> parameterize(const std::vector<Vec<D>> &points,
                                 Parameterization kind = Parameterization::ChordLength);

struct FitOptions {
    int degree = 3;
    // Polynomial pieces of the clamped uniform knot vector, so
    // segments + degree control points. Ignored when `knots` is given.
    int segments = 3;
    // A clamped knot vector on [0, 1] to fit over instead. Curves that must
    // share control-point structure -- opposite sides of a Coons patch, say --
    // should be fitted over the same one.
    std::vector<double> knots;

    Parameterization parameterization = Parameterization::ChordLength;

    // Hold the first and last control points at the first and last data
    // points, exactly.
    bool pinEnds = true;

    // Relative to the largest diagonal of the normal equations. See the header.
    double regularisation = 1e-6;

    // A run with fewer than densify x (control points) points has each of its
    // segments cut evenly until it does not. Zero or less fits the points as
    // given.
    int densify = 4;

    // The curve is sampled at this many uniform steps to measure the deviation;
    // zero or less means max(64, 8 x control points).
    int deviationSamples = 0;
};

// How far a curve strays from the points it was built from: for each point, the
// distance to the nearest place on the curve (sampled), then the largest and
// the root mean square of those.
struct Deviation {
    double max = 0.0;
    double rms = 0.0;
};

template <std::size_t D>
struct FitResult {
    BSplineCurve<D> curve;
    // Measured against the points passed in, not the densified run.
    double maxDeviation = 0.0;
    double rmsDeviation = 0.0;
    double length = 0.0;          // of the run that was fitted
    int samples = 0;              // points in that run, after densification
    // Fewer data points than free control points: the regularisation, not the
    // data, placed some of them. Also set if the normal equations could not be
    // factored, in which case the interior control points lie on the chord.
    bool underdetermined = false;
};

// Least-squares fit. Fewer than two points give a curve collapsed onto the one
// point there is (or the origin). Throws std::invalid_argument on a knot vector
// that is not clamped to [0, 1] or a degree out of range.
template <std::size_t D>
FitResult<D> fitCurve(const std::vector<Vec<D>> &points, const FitOptions &options = FitOptions());

// The curve through every point, point i at parameterize(points)[i], over
// averagedKnots(), solved by the kernel. The degree is lowered to
// points.size() - 1 when there are too few points for it; a single point gives
// a constant curve of degree 1. Throws std::invalid_argument if two consecutive
// points share a parameter, which makes the system singular.
template <std::size_t D>
BSplineCurve<D> interpolateCurve(const std::vector<Vec<D>> &points, int degree = 3,
                                 Parameterization kind = Parameterization::ChordLength);

template <std::size_t D>
Deviation deviation(const BSplineCurve<D> &curve, const std::vector<Vec<D>> &points, int samples);

extern template std::vector<double> parameterize<2>(const std::vector<Vec<2>> &, Parameterization);
extern template std::vector<double> parameterize<3>(const std::vector<Vec<3>> &, Parameterization);
extern template FitResult<2> fitCurve<2>(const std::vector<Vec<2>> &, const FitOptions &);
extern template FitResult<3> fitCurve<3>(const std::vector<Vec<3>> &, const FitOptions &);
extern template BSplineCurve<2> interpolateCurve<2>(const std::vector<Vec<2>> &, int, Parameterization);
extern template BSplineCurve<3> interpolateCurve<3>(const std::vector<Vec<3>> &, int, Parameterization);
extern template Deviation deviation<2>(const BSplineCurve<2> &, const std::vector<Vec<2>> &, int);
extern template Deviation deviation<3>(const BSplineCurve<3> &, const std::vector<Vec<3>> &, int);

} // namespace geom

#endif // __GEOM_FITTING_HXX__
