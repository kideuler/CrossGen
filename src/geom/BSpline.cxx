#include "geom/BSpline.hxx"

#include <BSplCLib.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <math_Matrix.hxx>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include "geom/detail/Occ.hxx"

namespace geom {

using detail::Access;

namespace {

void requireKnots(const std::vector<double> &knots, int degree, std::size_t count,
                  const char *what) {
    if (degree < 1 || degree > kMaxDegree) {
        throw std::invalid_argument(std::string(what) + ": degree " + std::to_string(degree) +
                                    " is outside [1, " + std::to_string(kMaxDegree) + "]");
    }
    if (knots.size() != count + static_cast<std::size_t>(degree) + 1) {
        throw std::invalid_argument(std::string(what) + ": " + std::to_string(knots.size()) +
                                    " knots for " + std::to_string(count) +
                                    " control points of degree " + std::to_string(degree));
    }
    for (std::size_t i = 1; i < knots.size(); ++i) {
        if (!(knots[i] >= knots[i - 1])) {
            throw std::invalid_argument(std::string(what) + ": knots decrease at " +
                                        std::to_string(i));
        }
    }
    if (!(knots[count] > knots[degree])) {
        throw std::invalid_argument(std::string(what) + ": empty parameter domain");
    }
}

// The kernel's arc-length routines integrate to a tolerance. This one is well
// below anything a mesh can resolve and well above what Gauss quadrature on a
// polynomial span struggles to reach.
constexpr double kLengthTolerance = 1e-12;

} // namespace

// ---------------------------------------------------------------------------
// Knot vectors
// ---------------------------------------------------------------------------
std::vector<double> clampedUniformKnots(int degree, int segments) {
    if (degree < 0) degree = 0;
    if (segments < 1) segments = 1;
    const int n = segments + degree;
    std::vector<double> knot(static_cast<std::size_t>(n + degree + 1), 0.0);
    for (int i = 0; i <= degree; ++i) knot[i] = 0.0;
    for (int i = 1; i < segments; ++i) knot[degree + i] = static_cast<double>(i) / segments;
    for (int i = n; i < n + degree + 1; ++i) knot[i] = 1.0;
    return knot;
}

std::vector<double> averagedKnots(const std::vector<double> &params, int degree) {
    const int n = static_cast<int>(params.size());
    if (degree < 1 || n < degree + 1) {
        throw std::invalid_argument("averagedKnots: needs degree >= 1 and degree + 1 parameters");
    }
    std::vector<double> knot(static_cast<std::size_t>(n + degree + 1), 0.0);
    for (int i = n; i < n + degree + 1; ++i) knot[i] = 1.0;
    for (int j = 1; j < n - degree; ++j) {
        double s = 0.0;
        for (int i = j; i < j + degree; ++i) s += params[i];
        knot[j + degree] = s / degree;
    }
    return knot;
}

std::vector<double> grevilleAbscissae(const std::vector<double> &knots, int degree) {
    const int n = controlPointCount(knots, degree);
    if (degree < 1 || n < 1) {
        throw std::invalid_argument("grevilleAbscissae: needs degree >= 1 and a control point");
    }
    const TColStd_Array1OfReal flat = detail::toArray(knots);
    TColStd_Array1OfReal g(1, n);
    detail::guard("grevilleAbscissae", [&] { BSplCLib::BuildSchoenbergPoints(degree, flat, g); });
    return detail::toVector(g);
}

int basisFunctions(const std::vector<double> &knots, int degree, double u, double *N) {
    const int n = controlPointCount(knots, degree);
    if (degree < 0 || degree > kMaxDegree || n < 1) {
        throw std::invalid_argument("basisFunctions: degree " + std::to_string(degree) + " with " +
                                    std::to_string(knots.size()) + " knots");
    }
    u = std::max(knots[degree], std::min(knots[n], u));
    const TColStd_Array1OfReal flat = detail::toArray(knots);
    math_Matrix basis(1, 1, 1, degree + 1);
    int first = 0;
    const int status = detail::guard("basisFunctions", [&] {
        return BSplCLib::EvalBsplineBasis(0, degree + 1, flat, u, first, basis);
    });
    if (status != 0) throw std::runtime_error("basisFunctions: the kernel could not evaluate the basis");
    for (int a = 0; a <= degree; ++a) N[a] = basis(1, a + 1);
    return first - 1;
}

// ---------------------------------------------------------------------------
// BSplineCurve
// ---------------------------------------------------------------------------
template <std::size_t D>
BSplineCurve<D>::BSplineCurve(int degree, const std::vector<double> &knots,
                              const std::vector<Vec<D>> &ctrl) {
    requireKnots(knots, degree, ctrl.size(), "BSplineCurve");
    TColgp_Array1OfPnt poles(1, static_cast<int>(ctrl.size()));
    for (std::size_t i = 0; i < ctrl.size(); ++i) poles(static_cast<int>(i) + 1) = detail::toPnt(ctrl[i]);
    const detail::SplitKnots k = detail::splitKnots(knots);
    Handle(Geom_BSplineCurve) h = detail::construct("BSplineCurve", [&] {
        return Handle(Geom_BSplineCurve)(new Geom_BSplineCurve(poles, k.knots, k.mults, degree));
    });
    data = std::make_shared<const detail::CurveData>(detail::CurveData{h});
}

template <std::size_t D>
int BSplineCurve<D>::size() const {
    return data ? data->curve->NbPoles() : 0;
}

template <std::size_t D>
int BSplineCurve<D>::degree() const {
    return data ? data->curve->Degree() : 0;
}

template <std::size_t D>
std::vector<double> BSplineCurve<D>::knots() const {
    if (!data) return {};
    const Handle(Geom_BSplineCurve) &c = data->curve;
    TColStd_Array1OfReal flat(1, c->NbPoles() + c->Degree() + 1);
    c->KnotSequence(flat);
    return detail::toVector(flat);
}

template <std::size_t D>
std::vector<Vec<D>> BSplineCurve<D>::controlPoints() const {
    std::vector<Vec<D>> out;
    if (!data) return out;
    const Handle(Geom_BSplineCurve) &c = data->curve;
    out.reserve(static_cast<std::size_t>(c->NbPoles()));
    for (int i = 1; i <= c->NbPoles(); ++i) out.push_back(detail::fromPnt<D>(c->Pole(i)));
    return out;
}

template <std::size_t D>
Vec<D> BSplineCurve<D>::controlPoint(int i) const {
    if (!data || i < 0 || i >= size()) throw std::out_of_range("BSplineCurve::controlPoint");
    return detail::fromPnt<D>(data->curve->Pole(i + 1));
}

template <std::size_t D>
void BSplineCurve<D>::setControlPoint(int i, const Vec<D> &p) {
    if (!data || i < 0 || i >= size()) throw std::out_of_range("BSplineCurve::setControlPoint");
    Handle(Geom_BSplineCurve) h = detail::copyOf(data->curve);
    h->SetPole(i + 1, detail::toPnt(p));
    data = std::make_shared<const detail::CurveData>(detail::CurveData{h});
}

template <std::size_t D>
double BSplineCurve<D>::domainBegin() const {
    return data ? data->curve->FirstParameter() : 0.0;
}

template <std::size_t D>
double BSplineCurve<D>::domainEnd() const {
    return data ? data->curve->LastParameter() : 0.0;
}

template <std::size_t D>
Vec<D> BSplineCurve<D>::evaluate(double u) const {
    if (!data) return zero<D>();
    const Handle(Geom_BSplineCurve) &c = data->curve;
    u = std::max(c->FirstParameter(), std::min(c->LastParameter(), u));
    gp_Pnt p;
    c->D0(u, p);
    return detail::fromPnt<D>(p);
}

template <std::size_t D>
Vec<D> BSplineCurve<D>::derivativeAt(double u, int order) const {
    if (!data) return zero<D>();
    if (order < 1) throw std::invalid_argument("BSplineCurve::derivativeAt: order must be >= 1");
    const Handle(Geom_BSplineCurve) &c = data->curve;
    u = std::max(c->FirstParameter(), std::min(c->LastParameter(), u));
    return detail::fromVec<D>(detail::guard("BSplineCurve::derivativeAt", [&] { return c->DN(u, order); }));
}

template <std::size_t D>
std::vector<Vec<D>> BSplineCurve<D>::sample(int segments) const {
    if (segments < 1) segments = 1;
    std::vector<Vec<D>> out;
    out.reserve(static_cast<std::size_t>(segments) + 1);
    const double a = domainBegin(), b = domainEnd();
    for (int i = 0; i < segments; ++i) {
        out.push_back(evaluate(a + (b - a) * (static_cast<double>(i) / segments)));
    }
    out.push_back(evaluate(b));
    return out;
}

// C'(u) = sum_i p (P_{i+1} - P_i) / (u_{i+p+1} - u_{i+1}) N_{i+1,p-1}(u): the
// hodograph, built from the kernel's control points and knots and handed back
// to it as a curve.
template <std::size_t D>
BSplineCurve<D> BSplineCurve<D>::derivative() const {
    if (!data) return BSplineCurve();
    const int p = degree();
    if (p < 2) {
        throw std::invalid_argument("BSplineCurve::derivative: the derivative of a degree-1 curve "
                                    "is piecewise constant; use derivativeAt()");
    }
    const Handle(Geom_BSplineCurve) &c = data->curve;
    for (int i = 2; i < c->NbKnots(); ++i) {
        if (c->Multiplicity(i) >= p) {
            throw std::invalid_argument("BSplineCurve::derivative: a knot of full multiplicity "
                                        "makes the derivative discontinuous; use derivativeAt()");
        }
    }
    const std::vector<double> knot = knots();
    const std::vector<Vec<D>> ctrl = controlPoints();
    std::vector<Vec<D>> q(ctrl.size() - 1);
    for (std::size_t i = 0; i + 1 < ctrl.size(); ++i) {
        const double h = knot[i + p + 1] - knot[i + 1];
        q[i] = h > 0.0 ? (ctrl[i + 1] - ctrl[i]) * (p / h) : zero<D>();
    }
    return BSplineCurve(p - 1, std::vector<double>(knot.begin() + 1, knot.end() - 1), q);
}

template <std::size_t D>
BSplineCurve<D> BSplineCurve<D>::reversed() const {
    if (!data) return *this;
    Handle(Geom_BSplineCurve) h = detail::copyOf(data->curve);
    detail::guard("BSplineCurve::reversed", [&] { h->Reverse(); });
    return Access::curve<D>(h);
}

template <std::size_t D>
int BSplineCurve<D>::insertKnot(double u, int times) {
    if (!data || times < 1) return 0;
    if (!(u > domainBegin() && u < domainEnd())) {
        throw std::invalid_argument("BSplineCurve::insertKnot: u is not inside the domain");
    }
    const Handle(Geom_BSplineCurve) &c = data->curve;
    int existing = 0;
    for (int i = 1; i <= c->NbKnots(); ++i) {
        if (c->Knot(i) == u) existing = c->Multiplicity(i);
    }
    const int r = std::min(times, c->Degree() - existing);
    if (r <= 0) return 0;
    Handle(Geom_BSplineCurve) h = detail::copyOf(c);
    detail::guard("BSplineCurve::insertKnot", [&] { h->InsertKnot(u, r, 0.0, Standard_True); });
    const int gained = h->NbPoles() - c->NbPoles();
    data = std::make_shared<const detail::CurveData>(detail::CurveData{h});
    return gained;
}

template <std::size_t D>
double BSplineCurve<D>::length() const {
    if (!data) return 0.0;
    return length(domainBegin(), domainEnd());
}

template <std::size_t D>
double BSplineCurve<D>::length(double u0, double u1) const {
    if (!data) return 0.0;
    const double a = domainBegin(), b = domainEnd();
    u0 = std::max(a, std::min(b, u0));
    u1 = std::max(a, std::min(b, u1));
    if (u0 == u1) return 0.0;
    const GeomAdaptor_Curve adaptor(data->curve);
    const double sign = u1 < u0 ? -1.0 : 1.0;
    return sign * detail::guard("BSplineCurve::length", [&] {
        return GCPnts_AbscissaPoint::Length(adaptor, std::min(u0, u1), std::max(u0, u1),
                                            kLengthTolerance);
    });
}

template <std::size_t D>
double BSplineCurve<D>::parameterAtLength(double s) const {
    if (!data) return 0.0;
    const double a = domainBegin(), b = domainEnd();
    const double total = length();
    if (!(s > 0.0) || !(total > 0.0)) return a;
    if (s >= total) return b;
    const GeomAdaptor_Curve adaptor(data->curve);
    return detail::guard("BSplineCurve::parameterAtLength", [&] {
        GCPnts_AbscissaPoint at(kLengthTolerance, adaptor, s, a);
        if (!at.IsDone()) throw std::runtime_error("BSplineCurve::parameterAtLength: no convergence");
        return std::max(a, std::min(b, at.Parameter()));
    });
}

template <std::size_t D>
typename BSplineCurve<D>::Nearest BSplineCurve<D>::nearest(const Vec<D> &p) const {
    Nearest best;
    best.distance = std::numeric_limits<double>::infinity();
    if (!data) return best;
    const Handle(Geom_BSplineCurve) &c = data->curve;
    const gp_Pnt target = detail::toPnt(p);
    auto consider = [&](double u) {
        const gp_Pnt q = c->Value(u);
        const double d = q.Distance(target);
        if (d < best.distance) {
            best.parameter = u;
            best.point = detail::fromPnt<D>(q);
            best.distance = d;
        }
    };
    // The ends and the corners, where the distance can be smallest without its
    // derivative vanishing.
    consider(c->FirstParameter());
    consider(c->LastParameter());
    for (int i = 2; i < c->NbKnots(); ++i) {
        if (c->Multiplicity(i) >= c->Degree()) consider(c->Knot(i));
    }
    // The stationary points in between, from the kernel.
    try {
        GeomAPI_ProjectPointOnCurve projection(target, c, c->FirstParameter(), c->LastParameter());
        for (int i = 1; i <= projection.NbPoints(); ++i) consider(projection.Parameter(i));
    } catch (const Standard_Failure &) {
        // No stationary point was found; the candidates above stand.
    }
    return best;
}

// ---------------------------------------------------------------------------
// BSplineSurface
// ---------------------------------------------------------------------------
template <std::size_t D>
BSplineSurface<D>::BSplineSurface(int degreeS, const std::vector<double> &knotsS,
                                  int degreeT, const std::vector<double> &knotsT,
                                  const std::vector<Vec<D>> &net) {
    const int nS = controlPointCount(knotsS, degreeS);
    const int nT = controlPointCount(knotsT, degreeT);
    if (nS < 1 || nT < 1 ||
        net.size() != static_cast<std::size_t>(nS) * static_cast<std::size_t>(nT)) {
        throw std::invalid_argument("BSplineSurface: the net is not countS x countT");
    }
    requireKnots(knotsS, degreeS, static_cast<std::size_t>(nS), "BSplineSurface (s)");
    requireKnots(knotsT, degreeT, static_cast<std::size_t>(nT), "BSplineSurface (t)");
    TColgp_Array2OfPnt poles(1, nS, 1, nT);
    for (int j = 0; j < nT; ++j) {
        for (int i = 0; i < nS; ++i) {
            poles(i + 1, j + 1) = detail::toPnt(net[static_cast<std::size_t>(j) * nS + i]);
        }
    }
    const detail::SplitKnots ks = detail::splitKnots(knotsS);
    const detail::SplitKnots kt = detail::splitKnots(knotsT);
    Handle(Geom_BSplineSurface) h = detail::construct("BSplineSurface", [&] {
        return Handle(Geom_BSplineSurface)(new Geom_BSplineSurface(
            poles, ks.knots, kt.knots, ks.mults, kt.mults, degreeS, degreeT));
    });
    data = std::make_shared<const detail::SurfaceData>(detail::SurfaceData{h});
}

template <std::size_t D>
int BSplineSurface<D>::degreeS() const { return data ? data->surface->UDegree() : 0; }
template <std::size_t D>
int BSplineSurface<D>::degreeT() const { return data ? data->surface->VDegree() : 0; }
template <std::size_t D>
int BSplineSurface<D>::countS() const { return data ? data->surface->NbUPoles() : 0; }
template <std::size_t D>
int BSplineSurface<D>::countT() const { return data ? data->surface->NbVPoles() : 0; }

template <std::size_t D>
std::vector<double> BSplineSurface<D>::knotsS() const {
    return data ? detail::toVector(data->surface->UKnotSequence()) : std::vector<double>();
}

template <std::size_t D>
std::vector<double> BSplineSurface<D>::knotsT() const {
    return data ? detail::toVector(data->surface->VKnotSequence()) : std::vector<double>();
}

template <std::size_t D>
std::vector<Vec<D>> BSplineSurface<D>::controlNet() const {
    std::vector<Vec<D>> out;
    if (!data) return out;
    const TColgp_Array2OfPnt &poles = data->surface->Poles();
    const int nS = countS(), nT = countT();
    out.reserve(static_cast<std::size_t>(nS) * nT);
    for (int j = 1; j <= nT; ++j) {
        for (int i = 1; i <= nS; ++i) out.push_back(detail::fromPnt<D>(poles(i, j)));
    }
    return out;
}

template <std::size_t D>
Vec<D> BSplineSurface<D>::control(int i, int j) const {
    if (!data || i < 0 || j < 0 || i >= countS() || j >= countT()) {
        throw std::out_of_range("BSplineSurface::control");
    }
    return detail::fromPnt<D>(data->surface->Pole(i + 1, j + 1));
}

template <std::size_t D>
void BSplineSurface<D>::setControl(int i, int j, const Vec<D> &p) {
    if (!data || i < 0 || j < 0 || i >= countS() || j >= countT()) {
        throw std::out_of_range("BSplineSurface::setControl");
    }
    Handle(Geom_BSplineSurface) h = detail::copyOf(data->surface);
    h->SetPole(i + 1, j + 1, detail::toPnt(p));
    data = std::make_shared<const detail::SurfaceData>(detail::SurfaceData{h});
}

template <std::size_t D>
Vec<D> BSplineSurface<D>::evaluate(double s, double t) const {
    if (!data) return zero<D>();
    const Handle(Geom_BSplineSurface) &S = data->surface;
    double s0, s1, t0, t1;
    S->Bounds(s0, s1, t0, t1);
    s = std::max(s0, std::min(s1, s));
    t = std::max(t0, std::min(t1, t));
    gp_Pnt p;
    S->D0(s, t, p);
    return detail::fromPnt<D>(p);
}

template <std::size_t D>
std::array<Vec<D>, 2> BSplineSurface<D>::partials(double s, double t) const {
    if (!data) return {zero<D>(), zero<D>()};
    const Handle(Geom_BSplineSurface) &S = data->surface;
    double s0, s1, t0, t1;
    S->Bounds(s0, s1, t0, t1);
    s = std::max(s0, std::min(s1, s));
    t = std::max(t0, std::min(t1, t));
    gp_Pnt p;
    gp_Vec ds, dt;
    S->D1(s, t, p, ds, dt);
    return {detail::fromVec<D>(ds), detail::fromVec<D>(dt)};
}

template <std::size_t D>
BSplineCurve<D> BSplineSurface<D>::row(int j) const {
    if (!data || j < 0 || j >= countT()) throw std::out_of_range("BSplineSurface::row");
    std::vector<Vec<D>> c(static_cast<std::size_t>(countS()));
    for (int i = 0; i < countS(); ++i) c[i] = control(i, j);
    return BSplineCurve<D>(degreeS(), knotsS(), c);
}

template <std::size_t D>
BSplineCurve<D> BSplineSurface<D>::column(int i) const {
    if (!data || i < 0 || i >= countS()) throw std::out_of_range("BSplineSurface::column");
    std::vector<Vec<D>> c(static_cast<std::size_t>(countT()));
    for (int j = 0; j < countT(); ++j) c[j] = control(i, j);
    return BSplineCurve<D>(degreeT(), knotsT(), c);
}

// ---------------------------------------------------------------------------
BSplineCurve<3> lift(const BSplineCurve<2> &c) { return Access::curve<3>(Access::curve(c)); }
BSplineCurve<2> flatten(const BSplineCurve<3> &c) { return Access::curve<2>(Access::curve(c)); }
BSplineSurface<3> lift(const BSplineSurface<2> &s) { return Access::surface<3>(Access::surface(s)); }
BSplineSurface<2> flatten(const BSplineSurface<3> &s) { return Access::surface<2>(Access::surface(s)); }

template class BSplineCurve<2>;
template class BSplineCurve<3>;
template class BSplineSurface<2>;
template class BSplineSurface<3>;

} // namespace geom
