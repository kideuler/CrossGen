// Self-test for the geometry layer of src/geom, and of the OpenCASCADE kernel
// under it.
//
//   TestGeom
//
// Every check is against a closed form, so the program needs no files and
// exits non-zero if any of them fails:
//
//   * the basis is a partition of unity, of any degree, on any knot vector, and
//     the kernel refuses what it cannot represent;
//   * a clamped curve interpolates its end control points exactly, and control
//     points on a line at their Greville abscissae give that line, linearly
//     parameterised (linear precision);
//   * derivative(), reversed() and insertKnot() leave the curve they came from
//     where it was, and arc length and the nearest point are what a ruler says;
//   * a tensor-product surface and a Coons net on four straight sides both
//     reproduce the bilinear map, the Coons net's boundary rows are its sides to
//     the bit, and the Coons surface of a polyline and three cubics is the Coons
//     blend of those curves as evaluated;
//   * fitCurve() reproduces what a cubic can represent, holds its ends to the
//     bit, and approximates a circle to the accuracy a cubic should;
//   * interpolateCurve() passes through its points;
//   * a polyline and an arc-length table return the lengths and parameters a
//     ruler would;
//   * faces that hold one edge make a valid shape that counts it once, and the
//     shape can be written as BREP and as STEP.
//
// MERIDIAN's and TORSION's Stage 9 are built on this code; their own checks
// (TestMERIDIAN, TestTORSION) are what say the pipeline still closes on top of
// it. Nothing here includes OpenCASCADE: it is tested through src/geom only.

#include <array>
#include <cmath>
#include <filesystem>
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
#include "geom/Topology.hxx"

namespace {

using geom::Vec;
using geom::operator+;
using geom::operator-;
using geom::operator*;
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

// A clamped, deliberately non-uniform knot vector on [0, 1], with a double knot
// at 0.35 wherever the degree leaves the curve continuous across one.
std::vector<double> unevenKnots(int degree) {
    std::vector<double> k(static_cast<size_t>(degree) + 1, 0.0);
    for (double x : {0.1, 0.35, 0.35, 0.6, 0.9}) {
        if (degree < 2 && !k.empty() && x == k.back()) continue;
        k.push_back(x);
    }
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
double dist3(const V3 &a, const V3 &b) { return geom::distance(a, b); }

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
    bool indexOk = true;
    for (int p = 0; p <= 5; ++p) {
        const std::vector<double> knots = unevenKnots(p);
        const int n = geom::controlPointCount(knots, p);
        for (int s = 0; s <= 400; ++s) {
            const double u = s / 400.0;
            double N[geom::kMaxDegree + 1];
            const int first = geom::basisFunctions(knots, p, u, N);
            indexOk = indexOk && first >= 0 && first + p < n;
            double sum = 0.0;
            for (int i = 0; i <= p; ++i) {
                sum += N[i];
                mostNegative = std::min(mostNegative, N[i]);
            }
            worstSum = std::max(worstSum, std::fabs(sum - 1.0));
        }
    }
    check(worstSum < 1e-14 && mostNegative >= 0.0 && indexOk,
          "degrees 0-5 on a non-uniform knot vector: the basis is non-negative, sums to 1 and "
          "stays inside the control polygon (worst " + sci(worstSum) + ")");
    double Nend[4];
    const int lastFirst = geom::basisFunctions(k, 3, 1.0, Nend);
    check(lastFirst == 2 && Nend[3] == 1.0, "at the end of a clamped domain the last function is 1");

    check(throws([] { geom::BSplineCurve<2>(3, {0, 0, 0, 0, 1, 1, 1}, std::vector<V2>(4)); }),
          "a knot vector one short for its control points is refused");
    check(throws([] { geom::BSplineCurve<2>(1, {0, 0, 0.6, 0.4, 1, 1}, std::vector<V2>(4)); }),
          "a decreasing knot vector is refused");
    check(throws([] { geom::BSplineCurve<2>(0, {0, 0.5, 1}, std::vector<V2>(2)); }),
          "degree 0 is refused: the kernel's curves are at least linear");
    check(throws([] { geom::BSplineCurve<2>(1, {0, 0, 0.5, 0.5, 1, 1}, std::vector<V2>(4)); }),
          "an interior knot above the degree -- a curve with a gap -- is refused");
    const double half = 0.5, next = std::nextafter(0.5, 1.0);
    check(throws([&] { geom::BSplineCurve<2>(1, {0, 0, half, next, 1, 1}, std::vector<V2>(4)); }),
          "two knots an ulp apart are refused rather than evaluated across a span of nothing");
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
    double daErr = 0.0;
    for (int s = 0; s <= 100; ++s) {
        const double u = s / 100.0;
        daErr = std::max(daErr, dist2(c.derivativeAt(u), dc.evaluate(u)));
    }
    check(daErr < 1e-12, "derivativeAt() is the derivative curve, pointwise (" + sci(daErr) + ")");
    check(throws([] { geom::BSplineCurve<2>(1, {0, 0, 1, 1}, {V2{0, 0}, V2{1, 1}}).derivative(); }),
          "derivative() of a degree-1 curve, which would be piecewise constant, is refused");

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

    // Copies share nothing they can change.
    geom::BSplineCurve<2> moved = c;
    moved.setControlPoint(2, V2{10.0, 10.0});
    check(moved.controlPoint(2) == V2{10.0, 10.0} && c.controlPoint(2) != V2{10.0, 10.0},
          "setControlPoint() on a copy leaves the original alone");

    // Arc length and its inverse, on a line a cubic represents exactly.
    const std::vector<double> uk = unevenKnots(3);
    const geom::BSplineCurve<2> ruler(3, uk, lineAtGrevilles(V2{1, 1}, V2{4, 5}, uk, 3));
    check(std::fabs(ruler.length() - 5.0) < 1e-12 && std::fabs(ruler.length(0.2, 0.7) - 2.5) < 1e-12,
          "length() of (1,1)-(4,5) is 5, and between u = 0.2 and 0.7 it is 2.5 (" +
              sci(std::fabs(ruler.length() - 5.0)) + ")");
    check(std::fabs(ruler.parameterAtLength(1.0) - 0.2) < 1e-10,
          "parameterAtLength(1) on it is 0.2");
    // ... and on a curve with no closed form, against a dense polyline of itself.
    double chord = 0.0;
    std::vector<V2> dense = c.sample(20000);
    for (size_t i = 1; i < dense.size(); ++i) chord += dist2(dense[i], dense[i - 1]);
    check(std::fabs(c.length() - chord) < 1e-6 * chord,
          "length() of a wiggly cubic agrees with 20000 chords (" +
              sci(std::fabs(c.length() - chord) / chord) + " relative)");

    // The nearest point: against brute force, and where it is an end.
    double nErr = 0.0;
    for (const V2 &q : {V2{0.3, 0.4}, V2{-2.0, 1.5}, V2{2.5, -2.0}, V2{0.0, 3.0}}) {
        double brute = 1e300;
        for (const V2 &p : dense) brute = std::min(brute, dist2(p, q));
        const auto near = c.nearest(q);
        nErr = std::max(nErr, std::fabs(near.distance - brute));
        nErr = std::max(nErr, std::fabs(dist2(near.point, q) - near.distance));
        nErr = std::max(nErr, dist2(c.evaluate(near.parameter), near.point));
    }
    check(nErr < 1e-6, "nearest() agrees with a 20000-point search (" + sci(nErr) + ")");
    const auto offEnd = ruler.nearest(V2{-2.0, -3.0});
    check(offEnd.parameter == 0.0 && std::fabs(offEnd.distance - 5.0) < 1e-12,
          "a point beyond the start of a line is nearest to the start, 5 away");
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

    const std::array<V2, 2> dS = S.partials(0.3, 0.6);
    const V2 dsWant = bilinear(1.0, 0.6) - bilinear(0.0, 0.6);
    const V2 dtWant = bilinear(0.3, 1.0) - bilinear(0.3, 0.0);
    check(dist2(dS[0], dsWant) < 1e-12 && dist2(dS[1], dtWant) < 1e-12,
          "partials() of the bilinear map are its two edge vectors");

    // makeCompatible() and a Coons surface on sides that do not share knots:
    // a polyline bottom, a cubic top, quadratic left and right.
    const geom::Polyline<2> zigzag({V2{0, 0}, V2{0.5, -0.2}, V2{1.2, 0.1}, V2{2, 0}});
    const std::vector<double> kc = geom::clampedUniformKnots(3, 3);
    const geom::BSplineCurve<2> arch(3, kc, {V2{-0.5, 1}, V2{0.2, 1.6}, V2{0.8, 1.2},
                                             V2{1.5, 1.9}, V2{2.1, 1.4}, V2{2.5, 1.5}});
    const geom::BSplineCurve<2> lside(2, kt, lineAtGrevilles(V2{0, 0}, V2{-0.5, 1}, kt, 2));
    const std::vector<double> kq = geom::clampedUniformKnots(2, 1);
    const geom::BSplineCurve<2> rside(2, kq, {V2{2, 0}, V2{2.6, 0.7}, V2{2.5, 1.5}});

    const auto [zc, ac] = geom::makeCompatible(zigzag.curve(), arch);
    double mcErr = 0.0;
    for (int i = 0; i <= 400; ++i) {
        const double u = i / 400.0;
        mcErr = std::max(mcErr, dist2(zc.evaluate(u), zigzag.evaluate(u)));
        mcErr = std::max(mcErr, dist2(ac.evaluate(u), arch.evaluate(u)));
    }
    check(zc.degree() == 3 && zc.knots() == ac.knots() && mcErr < 1e-14,
          "makeCompatible() gives a polyline and a cubic one degree and knot vector and moves "
          "neither (" + sci(mcErr) + ")");

    const geom::BSplineSurface<2> mixed = geom::coonsSurface(zigzag.curve(), arch, lside, rside);
    double mixErr = 0.0, sideErr = 0.0;
    for (int j = 0; j <= 30; ++j) {
        for (int i = 0; i <= 30; ++i) {
            const double s = i / 30.0, t = j / 30.0;
            const V2 want = geom::coonsPoint(zigzag.evaluate(s), arch.evaluate(s), lside.evaluate(t),
                                             rside.evaluate(t), V2{0, 0}, V2{2, 0}, V2{-0.5, 1},
                                             V2{2.5, 1.5}, s, t);
            mixErr = std::max(mixErr, dist2(mixed.evaluate(s, t), want));
        }
        const double w = j / 30.0;
        sideErr = std::max(sideErr, dist2(mixed.evaluate(w, 0.0), zigzag.evaluate(w)));
    }
    check(mixErr < 1e-14 && sideErr < 1e-14,
          "the Coons surface of a polyline, a cubic and two quadratics is their Coons blend as "
          "curves (" + sci(mixErr) + "), and runs along the polyline (" + sci(sideErr) + ")");
    check(throws([&] { geom::makeCompatible(arch, geom::BSplineCurve<2>(1, {0, 0, 2, 2}, {V2{0, 0}, V2{1, 1}})); }),
          "makeCompatible() refuses curves over different domains");
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
    const geom::BSplineCurve<2> one = geom::interpolateCurve(std::vector<V2>{V2{4, 5}});
    check(one.degree() == 1 && one.evaluate(0.3) == V2{4, 5},
          "one point interpolates as a constant, of degree 1");
    check(throws([] { geom::interpolateCurve(std::vector<V2>{V2{0, 0}, V2{1, 0}, V2{1, 0}, V2{2, 1}}); }),
          "a repeated point is refused rather than solved as a singular system");
}

// ---------------------------------------------------------------------------
void testPolylineAndArcLength() {
    heading("Polyline and ArcLengthTable");

    const geom::Polyline<2> pl({V2{0, 0}, V2{3, 0}, V2{3, 4}});
    check(pl.length() == 7.0, "(0,0) (3,0) (3,4) has length 7");
    check(pl.parameters().front() == 0.0 && pl.parameters().back() == 1.0 &&
              dist2(pl.evaluate(3.0 / 7.0), V2{3, 0}) < 1e-15,
          "its parameter is normalised chord length: the corner is at 3/7");
    check(dist2(pl.evaluate(5.0 / 7.0), V2{3, 2}) < 1e-15, "and 5/7 of the way is (3, 2)");
    check(pl.curve().degree() == 1 && pl.curve().size() == 3 && pl.curve().knots()[2] == 3.0 / 7.0,
          "to the kernel it is a degree-1 B-spline with its knots at the vertex parameters");
    check(std::fabs(pl.distance(V2{1, 1}) - 1.0) < 1e-15 && pl.distance(V2{3, 2}) < 1e-15 &&
              std::fabs(pl.distance(V2{4, -1}) - std::sqrt(2.0)) < 1e-15 &&
              std::fabs(pl.distance(V2{3, 6}) - 2.0) < 1e-15,
          "distance() is to the nearest segment, the corner or the end");

    const geom::Polyline<2> stutter({V2{0, 0}, V2{1, 0}, V2{1, 0}, V2{1, 1}, V2{1, 1}});
    double stErr = 0.0;
    for (int i = 0; i <= 100; ++i) {
        const double u = i / 100.0;
        const V2 want = u <= 0.5 ? V2{2 * u, 0} : V2{1, 2 * u - 1};
        stErr = std::max(stErr, dist2(stutter.evaluate(u), want));
    }
    check(stutter.size() == 5 && stutter.curve().size() == 3 && stErr < 1e-15,
          "repeated vertices stay in points() but not in the kernel's curve, which is unchanged "
          "by losing them");

    const geom::Polyline<2> still({V2{1, 1}, V2{1, 1}, V2{1, 1}});
    check(still.length() == 0.0 && still.parameters()[1] == 0.5 && still.evaluate(0.7) == V2{1, 1},
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

// ---------------------------------------------------------------------------
// The B-rep: two unit squares side by side, sharing the edge x = 1.
//
//     v01 ---- v11 ---- v21
//      |        |        |
//      |   F0   e   F1   |
//      |        |        |
//     v00 ---- v10 ---- v20
// ---------------------------------------------------------------------------
void testTopology() {
    heading("Vertex, Edge, Face and Shape");

    const std::vector<double> kc = geom::clampedUniformKnots(3, 3);
    auto straight = [&](const V2 &a, const V2 &b) {
        return geom::BSplineCurve<2>(3, kc, lineAtGrevilles(a, b, kc, 3));
    };
    const V2 p00{0, 0}, p10{1, 0}, p20{2, 0}, p01{0, 1}, p11{1, 1}, p21{2, 1};
    const geom::Vertex v00(p00), v10(p10), v20(p20), v01(p01), v11(p11), v21(p21);

    const geom::Edge b0(straight(p00, p10), v00, v10), t0(straight(p01, p11), v01, v11);
    const geom::Edge b1(straight(p10, p20), v10, v20), t1(straight(p11, p21), v11, v21);
    const geom::Edge left(straight(p00, p01), v00, v01), shared(straight(p10, p11), v10, v11);
    // The right side is drawn downwards, so the face has to be told it runs
    // against its frame.
    const geom::Edge right(straight(p21, p20), v21, v20);

    const geom::BSplineSurface<2> S0 = geom::coonsSurface(straight(p00, p10), straight(p01, p11),
                                                          straight(p00, p01), straight(p10, p11));
    const geom::BSplineSurface<2> S1 = geom::coonsSurface(straight(p10, p20), straight(p11, p21),
                                                          straight(p10, p11), straight(p20, p21));
    const geom::Face F0(S0, {b0, shared, t0, left}, {true, true, true, true});
    const geom::Face F1(S1, {b1, right, t1, shared}, {true, false, true, true});

    check(shared.start().isSame(v10) && shared.end().isSame(v11) && !shared.start().isSame(v11),
          "an edge holds the vertices it was given");
    check(std::fabs(shared.length() - 1.0) < 1e-12 && dist3(shared.evaluate(0.5), V3{1, 0.5, 0}) < 1e-15,
          "an edge measures its curve, which lies in z = 0");
    check(F0.isValid() && F1.isValid(), "both faces pass the kernel's check");
    check(std::fabs(F0.area() - 1.0) < 1e-12, "a unit square has area 1");

    int holdsShared = 0;
    for (const geom::Face *F : {&F0, &F1}) {
        for (const geom::Edge &e : F->edges()) holdsShared += e.isSame(shared);
    }
    check(holdsShared == 2 && F0.edges().size() == 4, "each face lists four edges, one of them the shared one");

    const geom::Shape shape({F0, F1});
    check(shape.faceCount() == 2 && shape.edgeCount() == 7 && shape.vertexCount() == 6,
          "the shape has 2 faces, 7 edges and 6 vertices: what is shared is counted once");
    check(shape.sharedEdgeCount() == 1 && shape.freeEdgeCount() == 6,
          "one edge bounds both faces, and the six on the outside bound one each");
    check(shape.isValid() && std::fabs(shape.area() - 2.0) < 1e-12, "it is valid, with area 2");

    // The same two squares with the middle edge built twice: geometrically
    // identical, topologically two edges and a seam.
    const geom::Edge twin(straight(p10, p11), v10, v11);
    const geom::Face F1b(S1, {b1, right, t1, twin}, {true, false, true, true});
    const geom::Shape cracked({F0, F1b});
    check(cracked.sharedEdgeCount() == 0 && cracked.freeEdgeCount() == 8,
          "built from a second copy of the middle edge, nothing is shared: the crack is visible");

    // The same neighbour with its (s, t) frame turned over -- s up the shared
    // edge, t across -- so that its surface faces -z. The shape turns it back.
    const geom::BSplineSurface<2> S1flip = geom::coonsSurface(
        straight(p10, p11), straight(p20, p21), straight(p10, p20), straight(p11, p21));
    const geom::Face F1flip(S1flip, {shared, t1, right, b1}, {true, true, false, true});
    const geom::Shape turned({F0, F1flip});
    check(F1flip.isValid() && turned.isValid() && turned.sharedEdgeCount() == 1,
          "a neighbour whose frame faces the other way is oriented to match: the shape is valid");

    // What the kernel refuses.
    check(throws([&] { geom::Edge(straight(p00, p10), v00, v11); }),
          "an edge whose vertex is not at the end of its curve is refused");
    check(throws([&] { geom::Edge(straight(p00, p00), v00, v00); }),
          "an edge on a curve of no length is refused");
    check(throws([&] { geom::Face(S0, {b0, shared, t0, left}, {true, true, true, false}); }),
          "a face whose edges do not meet at its corners is refused");
    check(throws([&] { geom::Face(S0, {b1, right, t1, shared}, {true, false, true, true}); }),
          "a face whose edges are not on its surface is refused");

    // A patch with a polyline side, as SplineFit builds one: the Coons surface of
    // the sides as they are, bounded by an edge on the polyline itself.
    const geom::Polyline<2> wall({p00, V2{0.3, -0.1}, V2{0.7, 0.08}, p10});
    const geom::Edge wallEdge(wall.curve(), v00, v10);
    const geom::BSplineSurface<2> SW =
        geom::coonsSurface(wall.curve(), straight(p01, p11), straight(p00, p01), straight(p10, p11));
    const geom::Face FW(SW, {wallEdge, shared, t0, left}, {true, true, true, true});
    const geom::Shape walled({FW, F1});
    check(FW.isValid() && walled.isValid() && walled.sharedEdgeCount() == 1,
          "a face on a polyline side is valid and still shares its other edge");

    // A face that is not flat.
    auto lifted = [&](const V3 &a, const V3 &b, double bulge) {
        std::vector<V3> ctrl = lineAtGrevilles(a, b, kc, 3);
        for (size_t i = 1; i + 1 < ctrl.size(); ++i) ctrl[i][2] += bulge;
        return geom::BSplineCurve<3>(3, kc, ctrl);
    };
    const V3 q00{0, 0, 0}, q10{1, 0, 0}, q01{0, 1, 0}, q11{1, 1, 0};
    const geom::BSplineCurve<3> cb = lifted(q00, q10, 0.3), ct = lifted(q01, q11, 0.3);
    const geom::BSplineCurve<3> cl = lifted(q00, q01, -0.3), cr = lifted(q10, q11, -0.3);
    const geom::Vertex w00(q00), w10(q10), w01(q01), w11(q11);
    const geom::Face saddle(geom::coonsSurface(cb, ct, cl, cr),
                            {geom::Edge(cb, w00, w10), geom::Edge(cr, w10, w11),
                             geom::Edge(ct, w01, w11), geom::Edge(cl, w00, w01)},
                            {true, true, true, true});
    check(saddle.isValid() && saddle.area() > 1.0,
          "a saddle in 3-D is a valid face, with more area than its unit shadow (" +
              std::to_string(saddle.area()) + ")");

    // Out to files.
    const std::filesystem::path dir = std::filesystem::temp_directory_path();
    const std::filesystem::path brep = dir / "crossgen_testgeom.brep";
    const std::filesystem::path step = dir / "crossgen_testgeom.step";
    const bool wroteBrep = shape.writeBREP(brep.string());
    const bool wroteStep = shape.writeSTEP(step.string());
    check(wroteBrep && std::filesystem::file_size(brep) > 0 && wroteStep &&
              std::filesystem::file_size(step) > 0,
          "the shape writes as BREP and as STEP");
    for (const auto &[what, back] : {std::make_pair(std::string("BREP"), wroteBrep ? geom::Shape::readBREP(brep.string()) : geom::Shape()),
                                     std::make_pair(std::string("STEP"), wroteStep ? geom::Shape::readSTEP(step.string()) : geom::Shape())}) {
        check(back.faceCount() == 2 && back.edgeCount() == 7 && back.sharedEdgeCount() == 1 &&
                  back.freeEdgeCount() == 6 && std::fabs(back.area() - 2.0) < 1e-9 && back.isValid(),
              "read back from " + what + ", it is still 2 faces sharing 1 of 7 edges, area 2");
    }
    bool missingThrows = false;
    try {
        geom::Shape::readSTEP((dir / "crossgen_no_such_file.step").string());
    } catch (const std::runtime_error &) {
        missingThrows = true;
    }
    check(missingThrows, "reading a file that is not there throws");
    std::error_code ignored;
    std::filesystem::remove(brep, ignored);
    std::filesystem::remove(step, ignored);
}

} // namespace

int main() {
    std::cout << "TestGeom -- the geometry layer of src/geom, on OpenCASCADE\n";
    testBasis();
    testCurve();
    testSurfaceAndCoons();
    testFitting();
    testPolylineAndArcLength();
    testTopology();

    heading("Self-test result");
    if (failures == 0) {
        std::cout << "  " << kPass << " Every closed-form check held.\n";
    } else {
        std::cout << "  " << kFail << " " << failures << " check(s) failed.\n";
    }
    return failures == 0 ? 0 : 1;
}
