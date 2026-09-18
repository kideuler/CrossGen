#include "geom/Fitting.hxx"

#include <BSplCLib.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <math_Matrix.hxx>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "geom/Polyline.hxx"
#include "geom/detail/Occ.hxx"

namespace geom {

namespace {

// The non-zero basis functions at a run of parameters, from the kernel, over
// one knot vector converted once rather than per point.
class Basis {
public:
    Basis(const std::vector<double> &knots, int degree)
        : p(degree), flat(detail::toArray(knots)), values(1, 1, 1, degree + 1) {}

    // Fills N[0..p] and returns the index of the control point N[0] weights.
    int at(double u, double *N) {
        int first = 0;
        const int status = detail::guard("fitCurve", [&] {
            return BSplCLib::EvalBsplineBasis(0, p + 1, flat, u, first, values);
        });
        if (status != 0) throw std::runtime_error("fitCurve: the kernel could not evaluate the basis");
        for (int a = 0; a <= p; ++a) N[a] = values(1, a + 1);
        return first - 1;
    }

private:
    int p;
    TColStd_Array1OfReal flat;
    math_Matrix values;
};

// Cholesky on a small dense SPD system, in place, for `rhs` right-hand sides
// laid end to end in b. The systems a fit produces are (control points) square
// -- a handful of rows -- so nothing more is warranted, and a factorisation
// that fails says the fit was rank deficient, which the regularisation is
// there to prevent.
bool choleskySolve(std::vector<double> &A, int n, std::vector<double> &b, int rhs) {
    for (int j = 0; j < n; ++j) {
        double d = A[j * n + j];
        for (int k = 0; k < j; ++k) d -= A[j * n + k] * A[j * n + k];
        if (!(d > 0.0)) return false;
        A[j * n + j] = std::sqrt(d);
        for (int i = j + 1; i < n; ++i) {
            double s = A[i * n + j];
            for (int k = 0; k < j; ++k) s -= A[i * n + k] * A[j * n + k];
            A[i * n + j] = s / A[j * n + j];
        }
    }
    for (int c = 0; c < rhs; ++c) {
        double *x = &b[c * n];
        for (int i = 0; i < n; ++i) {
            double s = x[i];
            for (int k = 0; k < i; ++k) s -= A[i * n + k] * x[k];
            x[i] = s / A[i * n + i];
        }
        for (int i = n - 1; i >= 0; --i) {
            double s = x[i];
            for (int k = i + 1; k < n; ++k) s -= A[k * n + i] * x[k];
            x[i] = s / A[i * n + i];
        }
    }
    return true;
}

void requireClampedUnit(const std::vector<double> &knots, int degree) {
    const int n = controlPointCount(knots, degree);
    if (degree < 1 || degree > kMaxDegree || n < degree + 1) {
        throw std::invalid_argument("fitCurve: degree " + std::to_string(degree) + " with " +
                                    std::to_string(knots.size()) + " knots");
    }
    for (int i = 0; i <= degree; ++i) {
        if (knots[i] != 0.0 || knots[knots.size() - 1 - i] != 1.0) {
            throw std::invalid_argument("fitCurve: the knot vector is not clamped to [0, 1]");
        }
    }
}

} // namespace

// ---------------------------------------------------------------------------
template <std::size_t D>
std::vector<double> parameterize(const std::vector<Vec<D>> &points, Parameterization kind) {
    std::vector<double> t(points.size(), 0.0);
    if (points.size() < 2) return t;
    double total = 0.0;
    if (kind != Parameterization::Uniform) {
        for (std::size_t i = 1; i < points.size(); ++i) {
            const double d = norm(points[i] - points[i - 1]);
            total += (kind == Parameterization::Centripetal) ? std::sqrt(d) : d;
            t[i] = total;
        }
    }
    if (total > 0.0) {
        for (double &x : t) x /= total;
    } else {
        for (std::size_t i = 0; i < t.size(); ++i) t[i] = static_cast<double>(i) / (t.size() - 1.0);
    }
    return t;
}

// ---------------------------------------------------------------------------
template <std::size_t D>
Deviation deviation(const BSplineCurve<D> &curve, const std::vector<Vec<D>> &points, int samples) {
    Deviation out;
    if (points.empty() || curve.empty()) return out;
    const std::vector<Vec<D>> sampled = curve.sample(std::max(1, samples));
    double sum = 0.0;
    for (const Vec<D> &q : points) {
        double best = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i + 1 < sampled.size(); ++i) {
            best = std::min(best, distanceToSegment(q, sampled[i], sampled[i + 1]));
        }
        out.max = std::max(out.max, best);
        sum += best * best;
    }
    out.rms = std::sqrt(sum / points.size());
    return out;
}

// ---------------------------------------------------------------------------
// fitCurve()
//
// The normal equations of min sum_k |C(t_k) - d_k|^2 over the free control
// points, with the pinned ones moved to the right-hand side, plus the Tikhonov
// term. Each data point touches only the degree + 1 basis functions of its
// span -- the kernel's values -- and a basis function that is exactly zero
// there is skipped.
// ---------------------------------------------------------------------------
template <std::size_t D>
FitResult<D> fitCurve(const std::vector<Vec<D>> &points, const FitOptions &options) {
    const int p = options.degree;
    std::vector<double> knots =
        options.knots.empty() ? clampedUniformKnots(p, options.segments) : options.knots;
    requireClampedUnit(knots, p);
    const int n = controlPointCount(knots, p);
    const std::vector<double> greville = grevilleAbscissae(knots, p);

    FitResult<D> result;
    std::vector<Vec<D>> ctrl(static_cast<std::size_t>(n), zero<D>());
    if (points.size() < 2) {
        if (!points.empty()) std::fill(ctrl.begin(), ctrl.end(), points.front());
        result.curve = BSplineCurve<D>(p, knots, ctrl);
        return result;
    }

    // Densify. The last point of each segment is the original point itself and
    // not a + (b - a), which is a rounding away from it: the ends of a run are
    // where other curves meet it, and they have to arrive at the same double.
    std::vector<Vec<D>> data;
    const int want = options.densify * n;
    if (options.densify <= 0 || static_cast<int>(points.size()) >= want) {
        data = points;
    } else {
        const int k = static_cast<int>(std::ceil(static_cast<double>(want) / (points.size() - 1)));
        data.push_back(points.front());
        for (std::size_t i = 1; i < points.size(); ++i) {
            for (int j = 1; j < k; ++j) {
                data.push_back(points[i - 1] + (points[i] - points[i - 1]) *
                                                   (static_cast<double>(j) / k));
            }
            data.push_back(points[i]);
        }
    }

    const std::vector<double> t = parameterize(data, options.parameterization);
    result.length = polylineLength(data);
    result.samples = static_cast<int>(data.size());

    const Vec<D> &p0 = data.front();
    const Vec<D> &pn = data.back();
    const bool pin = options.pinEnds;
    if (pin) {
        ctrl.front() = p0;
        ctrl.back() = pn;
    }

    // Free control points are [first, first + free).
    const int first = pin ? 1 : 0;
    const int free = pin ? n - 2 : n;
    if (free > 0) {
        const int equations = static_cast<int>(data.size()) - (pin ? 2 : 0);
        result.underdetermined = equations < free;

        std::vector<double> A(static_cast<std::size_t>(free) * free, 0.0);
        std::vector<double> rhs(static_cast<std::size_t>(free) * D, 0.0);
        const std::size_t kBegin = pin ? 1 : 0;
        const std::size_t kEnd = pin ? data.size() - 1 : data.size();
        Basis basis(knots, p);
        double N[kMaxDegree + 1];
        int idx[kMaxDegree + 1];
        double val[kMaxDegree + 1];
        for (std::size_t k = kBegin; k < kEnd; ++k) {
            const int firstIndex = basis.at(t[k], N);
            // The residual with the held control points already subtracted.
            Vec<D> b = data[k];
            int nz = 0;
            for (int a = 0; a <= p; ++a) {
                const int i = firstIndex + a;
                if (pin && i == 0) b = b - p0 * N[a];
                else if (pin && i == n - 1) b = b - pn * N[a];
                else if (N[a] != 0.0) { idx[nz] = i - first; val[nz] = N[a]; ++nz; }
            }
            for (int x = 0; x < nz; ++x) {
                for (int y = 0; y < nz; ++y) {
                    A[static_cast<std::size_t>(idx[x]) * free + idx[y]] += val[x] * val[y];
                }
                for (std::size_t d = 0; d < D; ++d) rhs[d * free + idx[x]] += val[x] * b[d];
            }
        }

        // Tikhonov towards the straight line between the ends, at the Greville
        // abscissae. See FitOptions::regularisation.
        double maxDiag = 0.0;
        for (int i = 0; i < free; ++i) maxDiag = std::max(maxDiag, A[static_cast<std::size_t>(i) * free + i]);
        const double reg = std::max(options.regularisation * maxDiag, maxDiag > 0.0 ? 0.0 : 1.0);
        for (int i = 0; i < free; ++i) {
            const double g = greville[i + first];
            const Vec<D> prior = p0 * (1.0 - g) + pn * g;
            A[static_cast<std::size_t>(i) * free + i] += reg;
            for (std::size_t d = 0; d < D; ++d) rhs[d * free + i] += reg * prior[d];
        }

        std::vector<double> chol = A;
        if (!choleskySolve(chol, free, rhs, static_cast<int>(D))) {
            for (int i = 0; i < free; ++i) {
                const double g = greville[i + first];
                ctrl[i + first] = p0 * (1.0 - g) + pn * g;
            }
            result.underdetermined = true;
        } else {
            for (int i = 0; i < free; ++i) {
                for (std::size_t d = 0; d < D; ++d) ctrl[i + first][d] = rhs[d * free + i];
            }
        }
    }

    result.curve = BSplineCurve<D>(p, knots, ctrl);

    // The error the fit actually made, measured against the points it was given
    // and not against its own parameterisation.
    const int ns = options.deviationSamples > 0 ? options.deviationSamples : std::max(64, 8 * n);
    const Deviation dev = deviation(result.curve, points, ns);
    result.maxDeviation = dev.max;
    result.rmsDeviation = dev.rms;
    return result;
}

// ---------------------------------------------------------------------------
// interpolateCurve()
//
// The collocation system N_i(t_k) P_i = Q_k is banded, square and, with
// averaged knots and distinct parameters, non-singular (Schoenberg-Whitney).
// The kernel solves it (BSplCLib::Interpolate, a banded LU).
// ---------------------------------------------------------------------------
template <std::size_t D>
BSplineCurve<D> interpolateCurve(const std::vector<Vec<D>> &points, int degree,
                                 Parameterization kind) {
    const int m = static_cast<int>(points.size());
    if (m == 0) return BSplineCurve<D>();
    if (m == 1) return BSplineCurve<D>(1, {0.0, 0.0, 1.0, 1.0}, {points[0], points[0]});
    const int p = std::max(1, std::min(degree, m - 1));
    if (p > kMaxDegree) throw std::invalid_argument("interpolateCurve: degree out of range");

    const std::vector<double> t = parameterize(points, kind);
    for (int k = 1; k < m; ++k) {
        if (!(t[k] > t[k - 1])) {
            throw std::invalid_argument("interpolateCurve: points " + std::to_string(k - 1) +
                                        " and " + std::to_string(k) + " share a parameter");
        }
    }
    const std::vector<double> knots = averagedKnots(t, p);

    const TColStd_Array1OfReal flat = detail::toArray(knots);
    const TColStd_Array1OfReal params = detail::toArray(t);
    TColStd_Array1OfInteger contact(1, m);
    contact.Init(0);
    TColgp_Array1OfPnt poles(1, m);
    for (int k = 0; k < m; ++k) poles(k + 1) = detail::toPnt(points[k]);
    int problem = 0;
    detail::guard("interpolateCurve", [&] {
        BSplCLib::Interpolate(p, flat, params, contact, poles, problem);
    });
    if (problem != 0) {
        throw std::runtime_error("interpolateCurve: the collocation system is singular");
    }

    std::vector<Vec<D>> ctrl(static_cast<std::size_t>(m));
    for (int i = 0; i < m; ++i) ctrl[i] = detail::fromPnt<D>(poles(i + 1));
    // The ends are interpolated by construction; say so exactly.
    ctrl.front() = points.front();
    ctrl.back() = points.back();
    return BSplineCurve<D>(p, knots, ctrl);
}

template std::vector<double> parameterize<2>(const std::vector<Vec<2>> &, Parameterization);
template std::vector<double> parameterize<3>(const std::vector<Vec<3>> &, Parameterization);
template FitResult<2> fitCurve<2>(const std::vector<Vec<2>> &, const FitOptions &);
template FitResult<3> fitCurve<3>(const std::vector<Vec<3>> &, const FitOptions &);
template BSplineCurve<2> interpolateCurve<2>(const std::vector<Vec<2>> &, int, Parameterization);
template BSplineCurve<3> interpolateCurve<3>(const std::vector<Vec<3>> &, int, Parameterization);
template Deviation deviation<2>(const BSplineCurve<2> &, const std::vector<Vec<2>> &, int);
template Deviation deviation<3>(const BSplineCurve<3> &, const std::vector<Vec<3>> &, int);

} // namespace geom
