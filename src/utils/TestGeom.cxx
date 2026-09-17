// Self-test for the spline library of src/geom.
//
//   TestGeom
//
// Every check is against a closed form, so the program needs no files and
// exits non-zero if any of them fails:
//
//   * the basis is a partition of unity, of any degree, on any knot vector;
//   * a clamped curve interpolates its end control points exactly, and control
//     points on a line at their Greville abscissae give that line, linearly
//     parameterised (linear precision);
//   * derivative(), reversed() and insertKnot() leave the curve they came from
//     where it was;
//   * a tensor-product surface and a Coons net on four straight sides both
//     reproduce the bilinear map, and the Coons net's boundary rows are its
//     sides to the bit;
//   * fitCurve() reproduces what a cubic can represent, holds its ends to the
//     bit, and approximates a circle to the accuracy a cubic should;
//   * interpolateCurve() passes through its points;
//   * a polyline and an arc-length table return the lengths and parameters a
//     ruler would.
//
// MERIDIAN's and TORSION's Stage 9 are built on this code; their own checks
// (TestMERIDIAN, TestTORSION) are what say the pipeline still closes on top of
// it.

#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "geom/ArcLength.hxx"
#include "geom/BSpline.hxx"
#include "geom/Coons.hxx"
#include "geom/Fitting.hxx"
#include "geom/Polyline.hxx"

namespace {

using geom::Vec;
using V2 = Vec<2>;
using V3 = Vec<3>;

const char *kPass = "\033[32m[PASS]\033[0m";
const char *kFail = "\033[31m[FAIL]\033[0m";

const double kPi = 3.14159265358979323846;

int failures = 0;

void check(bool ok, const std::string &what) {
    std::cout << "  " << (ok ? kPass : kFail) << " " << what << "\n";
    if (!ok) ++failures;
}

void heading(const std::string &title) {
    std::cout << "\n" << title << "\n";
    std::cout << std::string(title.size(), '-') << "\n";
}

std::string sci(double x) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(2) << x;
    return oss.str();
}

bool throws(const std::function<void()> &f) {
    try {
        f();
    } catch (const std::invalid_argument &) {
        return true;
    }
    return false;
}

// A clamped, deliberately non-uniform knot vector on [0, 1].
std::vector<double> unevenKnots(int degree) {
    std::vector<double> k(static_cast<size_t>(degree) + 1, 0.0);
    for (double x : {0.1, 0.35, 0.35, 0.6, 0.9}) k.push_back(x);
    for (int i = 0; i <= degree; ++i) k.push_back(1.0);
    return k;
}

// A curve with no special structure: control points on a lopsided spiral.
geom::BSplineCurve<2> wigglyCurve(int degree) {
    std::vector<double> knots = unevenKnots(degree);
    const int n = geom::controlPointCount(knots, degree);
    std::vector<V2> ctrl(n);
    for (int i = 0; i < n; ++i) {
        const double a = 0.9 * i;
        ctrl[i] = V2{(1.0 + 0.3 * i) * std::cos(a), (0.5 + 0.2 * i) * std::sin(a)};
    }
    return geom::BSplineCurve<2>(degree, knots, ctrl);
}

// Control points on the segment a -> b at the Greville abscissae of `knots`.
template <std::size_t D>
std::vector<Vec<D>> lineAtGrevilles(const Vec<D> &a, const Vec<D> &b,
                                    const std::vector<double> &knots, int degree) {
    std::vector<Vec<D>> ctrl;
    for (double g : geom::grevilleAbscissae(knots, degree)) {
        Vec<D> p;
        for (std::size_t d = 0; d < D; ++d) p[d] = a[d] + (b[d] - a[d]) * g;
        ctrl.push_back(p);
    }
    return ctrl;
}

double dist2(const V2 &a, const V2 &b) { return std::hypot(a[0] - b[0], a[1] - b[1]); }

// ---------------------------------------------------------------------------
void testBasis() {
    heading("Knots and the basis");

    const std::vector<double> k = geom::clampedUniformKnots(3, 3);
    const std::vector<double> want = {0, 0, 0, 0, 1.0 / 3, 2.0 / 3, 1, 1, 1, 1};
    check(k == want, "clampedUniformKnots(3, 3) is 0,0,0,0, 1/3, 2/3, 1,1,1,1");
    const std::vector<double> g = geom::grevilleAbscissae(k, 3);
    double gErr = 0.0;
    const std::vector<double> gWant = {0.0, 1.0 / 9, 1.0 / 3, 2.0 / 3, 8.0 / 9, 1.0};
    for (size_t i = 0; i < g.size(); ++i) gErr = std::max(gErr, std::fabs(g[i] - gWant[i]));
    check(g.size() == 6 && gErr < 1e-15, "its Greville abscissae are 0, 1/9, 1/3, 2/3, 8/9, 1");

    double worstSum = 0.0, mostNegative = 0.0;
    for (int p = 0; p <= 5; ++p) {
        const std::vector<double> knots = unevenKnots(p);
        for (int s = 0; s <= 400; ++s) {
            const double u = s / 400.0;
            const int span = geom::findSpan(knots, p, u);
            double N[geom::kMaxDegree + 1];
            geom::basisFunctions(knots, p, span, u, N);
            double sum = 0.0;
            for (int i = 0; i <= p; ++i) {
                sum += N[i];
                mostNegative = std::min(mostNegative, N[i]);
            }
            worstSum = std::max(worstSum, std::fabs(sum - 1.0));
        }
    }
    check(worstSum < 1e-14 && mostNegative >= 0.0,
          "degrees 0-5 on a non-uniform knot vector: the basis is non-negative and sums to 1"
          " (worst " + sci(worstSum) + ")");

    check(throws([] { geom::BSplineCurve<2>(3, {0, 0, 0, 0, 1, 1, 1}, std::vector<V2>(4)); }),
          "a knot vector one short for its control points is refused");
    check(throws([] { geom::BSplineCurve<2>(1, {0, 0, 0.6, 0.4, 1, 1}, std::vector<V2>(4)); }),
          "a decreasing knot vector is refused");
}

// ---------------------------------------------------------------------------
void testCurve() {
    heading("BSplineCurve");

    const geom::BSplineCurve<2> c = wigglyCurve(3);
    check(c.evaluate(0.0) == c.controlPoints().front() && c.evaluate(1.0) == c.controlPoints().back(),
          "a clamped curve starts and ends on its end control points, to the bit");
    check(c.evaluate(-0.5) == c.evaluate(0.0) && c.evaluate(7.0) == c.evaluate(1.0),
          "evaluate() clamps into the domain rather than extrapolating");

    // Linear precision on a non-uniform knot vector.
    for (int p : {1, 2, 3, 5}) {
        const std::vector<double> knots = unevenKnots(p);
        const V2 a{-1.0, 2.0}, b{3.0, -0.5};
        const geom::BSplineCurve<2> line(p, knots, lineAtGrevilles(a, b, knots, p));
        double err = 0.0;
        for (int s = 0; s <= 200; ++s) {
            const double u = s / 200.0;
            err = std::max(err, dist2(line.evaluate(u), V2{a[0] + (b[0] - a[0]) * u,
                                                           a[1] + (b[1] - a[1]) * u}));
        }
        check(err < 1e-14, "degree " + std::to_string(p) + ": control points on a line at their"
              " Grevilles give that line, linearly parameterised (" + sci(err) + ")");
    }

    // derivative() against a central difference.
    const geom::BSplineCurve<2> dc = c.derivative();
    double dErr = 0.0;
    const double h = 1e-6;
    for (int s = 1; s < 100; ++s) {
        const double u = s / 100.0;
        if (std::fabs(u - 0.35) < 2 * h) continue;  // a double knot: only C1 there
        const V2 fd{(c.evaluate(u + h)[0] - c.evaluate(u - h)[0]) / (2 * h),
                    (c.evaluate(u + h)[1] - c.evaluate(u - h)[1]) / (2 * h)};
        dErr = std::max(dErr, dist2(dc.evaluate(u), fd) / (1.0 + std::hypot(fd[0], fd[1])));
    }
    check(dc.degree() == 2 && dErr < 1e-6,
          "derivative() matches a central difference (" + sci(dErr) + " relative)");

    // reversed().
    const geom::BSplineCurve<2> r = c.reversed();
    double rErr = 0.0;
    for (int s = 0; s <= 200; ++s) {
        const double u = s / 200.0;
        rErr = std::max(rErr, dist2(c.evaluate(u), r.evaluate(1.0 - u)));
    }
    check(rErr < 1e-14, "reversed() traces the same points backwards (" + sci(rErr) + ")");

    // insertKnot().
    geom::BSplineCurve<2> k = c;
    const int n0 = k.size();
    const int once = k.insertKnot(0.47);
    const int capped = k.insertKnot(0.35, 5);   // already a double knot
    const int fresh = k.insertKnot(0.8, 7);     // capped at the degree
    double kErr = 0.0;
    for (int s = 0; s <= 256; ++s) {
        const double u = s / 256.0;
        kErr = std::max(kErr, dist2(c.evaluate(u), k.evaluate(u)));
    }
    check(once == 1 && capped == 1 && fresh == 3 && k.size() == n0 + 5,
          "insertKnot() caps multiplicity at the degree: 1 + (3 - 2) + 3 control points gained");
    check(kErr < 1e-14, "insertKnot() leaves the curve where it was (" + sci(kErr) + ")");
    check(throws([&] { k.insertKnot(1.0); }), "insertKnot() refuses a knot outside the domain");
}

// ---------------------------------------------------------------------------
void testSurfaceAndCoons() {
    heading("BSplineSurface and Coons blends");

    const V2 p00{0.0, 0.0}, p10{2.0, 0.0}, p01{-0.5, 1.0}, p11{2.5, 1.5};
    auto bilinear = [&](double s, double t) {
        return V2{(1 - s) * (1 - t) * p00[0] + s * (1 - t) * p10[0] + (1 - s) * t * p01[0] + s * t * p11[0],
                  (1 - s) * (1 - t) * p00[1] + s * (1 - t) * p10[1] + (1 - s) * t * p01[1] + s * t * p11[1]};
    };

    // Different knot vectors, and degrees, in the two directions.
    const std::vector<double> ks = unevenKnots(3);
    const std::vector<double> kt = geom::clampedUniformKnots(2, 4);
    const geom::BSplineCurve<2> bottom(3, ks, lineAtGrevilles(p00, p10, ks, 3));
    const geom::BSplineCurve<2> top(3, ks, lineAtGrevilles(p01, p11, ks, 3));
    const geom::BSplineCurve<2> left(2, kt, lineAtGrevilles(p00, p01, kt, 2));
    const geom::BSplineCurve<2> right(2, kt, lineAtGrevilles(p10, p11, kt, 2));
    const geom::BSplineSurface<2> S = geom::coonsSurface(bottom, top, left, right);

    double err = 0.0, pointErr = 0.0;
    for (int j = 0; j <= 20; ++j) {
        for (int i = 0; i <= 20; ++i) {
            const double s = i / 20.0, t = j / 20.0;
            err = std::max(err, dist2(S.evaluate(s, t), bilinear(s, t)));
            const V2 q = geom::coonsPoint(bottom.evaluate(s), top.evaluate(s), left.evaluate(t),
                                          right.evaluate(t), p00, p10, p01, p11, s, t);
            pointErr = std::max(pointErr, dist2(q, bilinear(s, t)));
        }
    }
    check(S.countS() == bottom.size() && S.countT() == left.size(),
          "coonsSurface() takes its size in s from the bottom and in t from the left side");
    check(err < 1e-14, "the Coons net of four straight sides, at the Grevilles, is the bilinear "
          "map (" + sci(err) + ")");
    check(pointErr < 1e-14, "coonsPoint() on the same four sides agrees (" + sci(pointErr) + ")");
    check(S.row(0).controlPoints() == bottom.controlPoints() &&
              S.row(S.countT() - 1).controlPoints() == top.controlPoints() &&
              S.column(0).controlPoints() == left.controlPoints() &&
              S.column(S.countS() - 1).controlPoints() == right.controlPoints(),
          "its boundary rows are the four side polygons, to the bit");

    double edgeErr = 0.0;
    for (int i = 0; i <= 50; ++i) {
        const double w = i / 50.0;
        edgeErr = std::max(edgeErr, dist2(S.evaluate(w, 0.0), bottom.evaluate(w)));
        edgeErr = std::max(edgeErr, dist2(S.evaluate(1.0, w), right.evaluate(w)));
    }
    check(edgeErr < 1e-14, "the surface runs along its side curves (" + sci(edgeErr) + ")");

    check(throws([&] { geom::coonsSurface(bottom, left, left, right); }),
          "coonsSurface() refuses opposite sides on different knot vectors");
}

// ---------------------------------------------------------------------------
void testFitting() {
    heading("fitCurve and interpolateCurve");

    // A straight run, unevenly spaced: a cubic can be exactly that.
    std::vector<V2> line;
    for (double x : {0.0, 0.05, 0.3, 0.31, 0.7, 0.72, 0.9, 1.4, 2.0}) line.push_back(V2{x, 0.5 * x - 1.0});
    const geom::FitResult<2> fl = geom::fitCurve(line);
    check(fl.curve.size() == 6 && fl.curve.degree() == 3,
          "the defaults are a cubic with three segments, six control points");
    check(fl.curve.controlPoints().front() == line.front() &&
              fl.curve.controlPoints().back() == line.back(),
          "pinEnds holds the first and last control points on the data, to the bit");
    check(fl.maxDeviation < 1e-12, "a straight run is fitted exactly (" + sci(fl.maxDeviation) + ")");
    check(std::fabs(fl.length - geom::polylineLength(line)) < 1e-14,
          "the reported length is the run's");

    geom::FitOptions free;
    free.pinEnds = false;
    const geom::FitResult<2> ff = geom::fitCurve(line, free);
    check(ff.maxDeviation < 1e-9, "with free ends too (" + sci(ff.maxDeviation) + ")");

    // A quarter circle. The knot vector is fixed, so this is approximation in
    // the parameter, not geometric Hermite: the error is O(h^4) with a constant
    // of the size of (pi/2)^4 / 384, about 1e-4 at three segments. The deviation
    // is measured against the curve sampled as a polyline, whose chords sit
    // inside the circle by their sagitta -- 7e-5 at the default 64 samples, which
    // would floor the measurement -- so it is sampled finely here.
    std::vector<V2> arc;
    for (int i = 0; i <= 60; ++i) {
        const double a = 0.5 * kPi * i / 60.0;
        arc.push_back(V2{std::cos(a), std::sin(a)});
    }
    auto radialError = [](const geom::BSplineCurve<2> &c) {
        double e = 0.0;
        for (int i = 0; i <= 2000; ++i) {
            const V2 p = c.evaluate(i / 2000.0);
            e = std::max(e, std::fabs(std::hypot(p[0], p[1]) - 1.0));
        }
        return e;
    };
    geom::FitOptions fine;
    fine.deviationSamples = 8192;
    const geom::FitResult<2> fc = geom::fitCurve(arc, fine);
    const double radial3 = radialError(fc.curve);
    check(fc.maxDeviation < 3e-4 && radial3 < 3e-4,
          "three segments hold a quarter circle to " + sci(fc.maxDeviation) + " at the data and " +
              sci(radial3) + " between");
    geom::FitOptions more = fine;
    more.segments = 6;
    const double radial6 = radialError(geom::fitCurve(arc, more).curve);
    check(radial6 < radial3 / 10,
          "six segments are " + sci(radial3 / radial6) + " times closer: fourth order");

    geom::FitOptions centripetal = fine;
    centripetal.parameterization = geom::Parameterization::Centripetal;
    const double radialC = radialError(geom::fitCurve(arc, centripetal).curve);
    check(std::fabs(radialC - radial3) < 1e-9,
          "on points evenly spaced, centripetal parameterisation is chord length");

    // Three points: densified, the fit is determined; left alone, it is not.
    const std::vector<V2> three = {V2{0, 0}, V2{1, 1}, V2{2, 0}};
    const geom::FitResult<2> dense = geom::fitCurve(three);
    geom::FitOptions sparse;
    sparse.densify = 0;
    const geom::FitResult<2> thin = geom::fitCurve(three, sparse);
    bool finite = true;
    for (const V2 &p : thin.curve.controlPoints()) finite = finite && std::isfinite(p[0]) && std::isfinite(p[1]);
    check(!dense.underdetermined && dense.samples > 20,
          "a three-point run is densified to " + std::to_string(dense.samples) + " points");
    check(thin.underdetermined && finite,
          "without densification it is flagged underdetermined, and the regularisation still "
          "gives finite control points");

    // 3-D, on a helix.
    std::vector<V3> helix;
    for (int i = 0; i <= 80; ++i) {
        const double a = 2.0 * kPi * i / 80.0;
        helix.push_back(V3{std::cos(a), std::sin(a), 0.3 * a});
    }
    geom::FitOptions h;
    h.segments = 16;
    h.deviationSamples = 8192;
    const geom::FitResult<3> fh = geom::fitCurve(helix, h);
    check(fh.maxDeviation < 1e-4, "a turn of a helix, in 3-D, sixteen segments: " + sci(fh.maxDeviation));

    // Interpolation.
    std::vector<V2> wave;
    for (int i = 0; i <= 14; ++i) wave.push_back(V2{0.3 * i, std::sin(0.7 * i) + 0.1 * i * i});
    const geom::BSplineCurve<2> wc = geom::interpolateCurve(wave);
    const std::vector<double> t = geom::parameterize(wave);
    double iErr = 0.0;
    for (size_t i = 0; i < wave.size(); ++i) iErr = std::max(iErr, dist2(wc.evaluate(t[i]), wave[i]));
    check(wc.size() == static_cast<int>(wave.size()) && iErr < 1e-12,
          "interpolateCurve() passes through all 15 points (" + sci(iErr) + ")");
    check(geom::interpolateCurve(three).degree() == 2,
          "three points cannot carry a cubic: the degree drops to 2");
    check(geom::interpolateCurve(std::vector<V2>{V2{4, 5}}).evaluate(0.3) == V2{4, 5},
          "one point interpolates as a constant");
    check(throws([] { geom::interpolateCurve(std::vector<V2>{V2{0, 0}, V2{1, 0}, V2{1, 0}, V2{2, 1}}); }),
          "a repeated point is refused rather than solved as a singular system");
}

// ---------------------------------------------------------------------------
void testPolylineAndArcLength() {
    heading("Polyline and ArcLengthTable");

    const geom::Polyline<2> pl({V2{0, 0}, V2{3, 0}, V2{3, 4}});
    check(pl.length() == 7.0, "(0,0) (3,0) (3,4) has length 7");
    check(pl.parameters().front() == 0.0 && pl.parameters().back() == 1.0 &&
              pl.evaluate(3.0 / 7.0) == V2{3, 0},
          "its parameter is normalised chord length: the corner is at 3/7");
    check(dist2(pl.evaluate(5.0 / 7.0), V2{3, 2}) < 1e-15, "and 5/7 of the way is (3, 2)");
    check(std::fabs(pl.distance(V2{1, 1}) - 1.0) < 1e-15 && pl.distance(V2{3, 2}) == 0.0,
          "distance() is to the nearest segment");

    const geom::Polyline<2> still({V2{1, 1}, V2{1, 1}, V2{1, 1}});
    check(still.length() == 0.0 && still.parameters()[1] == 0.5,
          "a polyline with no length spreads its vertices uniformly instead");

    // A circle at uniform speed: length pi / 2, and half the length is half the parameter.
    const geom::ArcLengthTable quarter(
        [](double u) { return V2{std::cos(0.5 * kPi * u), std::sin(0.5 * kPi * u)}; }, 256);
    check(std::fabs(quarter.length() - 0.5 * kPi) < 1e-5,
          "a quarter circle tabulated at 256 steps: length " + std::to_string(quarter.length()));
    check(std::fabs(quarter.parameterAt(0.5 * quarter.length()) - 0.5) < 1e-12,
          "at uniform speed, half the length is half the parameter");

    // u -> (u^2, 0): the parameter halfway along is 1/sqrt(2), not 1/2.
    const geom::ArcLengthTable accel([](double u) { return V2{u * u, 0.0}; }, 512);
    check(std::fabs(accel.parameterAt(0.5) - std::sqrt(0.5)) < 1e-5,
          "on u -> (u^2, 0) half the length is at u = 1/sqrt(2) (" +
              std::to_string(accel.parameterAt(0.5)) + ")");
    check(std::fabs(accel.lengthAt(accel.parameterAt(0.3)) - 0.3) < 1e-12,
          "lengthAt() inverts parameterAt()");
    check(accel.parameterAt(-1.0) == 0.0 && accel.parameterAt(2.0) == 1.0,
          "parameterAt() clamps at both ends");

    // A fitted spline, measured.
    std::vector<V2> arc;
    for (int i = 0; i <= 40; ++i) {
        const double a = kPi * i / 40.0;
        arc.push_back(V2{2.0 * std::cos(a), 2.0 * std::sin(a)});
    }
    const geom::FitResult<2> f = geom::fitCurve(arc);
    const geom::ArcLengthTable semi(f.curve, 512);
    check(std::fabs(semi.length() - 2.0 * kPi) < 1e-3,
          "a BSplineCurve goes straight into the table: a fitted half circle of radius 2 "
          "measures " + std::to_string(semi.length()));
}

} // namespace

int main() {
    std::cout << "TestGeom -- the spline library of src/geom\n";
    testBasis();
    testCurve();
    testSurfaceAndCoons();
    testFitting();
    testPolylineAndArcLength();

    heading("Self-test result");
    if (failures == 0) {
        std::cout << "  " << kPass << " Every closed-form check held.\n";
    } else {
        std::cout << "  " << kFail << " " << failures << " check(s) failed.\n";
    }
    return failures == 0 ? 0 : 1;
}
