#include "geom/BSpline.hxx"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>

namespace geom {

namespace {

void requireKnots(const std::vector<double> &knots, int degree, std::size_t count,
                  const char *what) {
    if (degree < 0 || degree > kMaxDegree) {
        throw std::invalid_argument(std::string(what) + ": degree " + std::to_string(degree) +
                                    " is outside [0, " + std::to_string(kMaxDegree) + "]");
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
    if (count > 0 && !(knots[count] > knots[degree])) {
        throw std::invalid_argument(std::string(what) + ": empty parameter domain");
    }
}

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

// Accumulated from zero in order: for degree 3 that is bit for bit the
// (k1 + k2 + k3) / 3.0 SplineFit wrote, since 0 + k1 is k1.
std::vector<double> grevilleAbscissae(const std::vector<double> &knots, int degree) {
    const int n = controlPointCount(knots, degree);
    std::vector<double> g(static_cast<std::size_t>(std::max(n, 0)), 0.0);
    for (int i = 0; i < n; ++i) {
        if (degree == 0) {
            g[i] = 0.5 * (knots[i] + knots[i + 1]);
            continue;
        }
        double s = 0.0;
        for (int j = 1; j <= degree; ++j) s += knots[i + j];
        g[i] = s / degree;
    }
    return g;
}

int findSpan(const std::vector<double> &knots, int degree, double u) {
    const int p = degree;
    const int n = controlPointCount(knots, degree);
    if (u >= knots[n]) return n - 1;
    if (u <= knots[p]) return p;
    int lo = p, hi = n, mid = (lo + hi) / 2;
    while (u < knots[mid] || u >= knots[mid + 1]) {
        if (u < knots[mid]) hi = mid; else lo = mid;
        mid = (lo + hi) / 2;
    }
    return mid;
}

void basisFunctions(const std::vector<double> &U, int degree, int span, double u, double *N) {
    const int p = degree;
    double left[kMaxDegree + 1], right[kMaxDegree + 1];
    N[0] = 1.0;
    for (int j = 1; j <= p; ++j) {
        left[j] = u - U[span + 1 - j];
        right[j] = U[span + j] - u;
        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            const double temp = N[r] / (right[r + 1] + left[j - r]);
            N[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        N[j] = saved;
    }
}

// ---------------------------------------------------------------------------
// BSplineCurve
// ---------------------------------------------------------------------------
template <std::size_t D>
BSplineCurve<D>::BSplineCurve(int degree, std::vector<double> knots, std::vector<Vec<D>> points)
    : p(degree), knot(std::move(knots)), ctrl(std::move(points)) {
    requireKnots(knot, p, ctrl.size(), "BSplineCurve");
}

template <std::size_t D>
double BSplineCurve<D>::domainBegin() const {
    return knot.empty() ? 0.0 : knot[p];
}

template <std::size_t D>
double BSplineCurve<D>::domainEnd() const {
    return knot.empty() ? 0.0 : knot[ctrl.size()];
}

template <std::size_t D>
Vec<D> BSplineCurve<D>::evaluate(double u) const {
    if (ctrl.empty()) return zero<D>();
    u = std::max(knot[p], std::min(knot[ctrl.size()], u));
    const int span = findSpan(knot, p, u);
    double N[kMaxDegree + 1];
    basisFunctions(knot, p, span, u, N);
    Vec<D> out = zero<D>();
    for (int k = 0; k <= p; ++k) out = out + ctrl[span - p + k] * N[k];
    return out;
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

// C'(u) = sum_i p (P_{i+1} - P_i) / (u_{i+p+1} - u_{i+1}) N_{i+1,p-1}(u). A
// zero denominator is a knot of multiplicity p + 1 inside the curve, across
// which it is discontinuous; the derivative is taken as zero there.
template <std::size_t D>
BSplineCurve<D> BSplineCurve<D>::derivative() const {
    if (p < 1) throw std::invalid_argument("BSplineCurve::derivative: degree-0 curve");
    if (ctrl.size() < 2) return BSplineCurve();
    std::vector<Vec<D>> q(ctrl.size() - 1);
    for (std::size_t i = 0; i + 1 < ctrl.size(); ++i) {
        const double h = knot[i + p + 1] - knot[i + 1];
        q[i] = h > 0.0 ? (ctrl[i + 1] - ctrl[i]) * (p / h) : zero<D>();
    }
    return BSplineCurve(p - 1, std::vector<double>(knot.begin() + 1, knot.end() - 1),
                        std::move(q));
}

template <std::size_t D>
BSplineCurve<D> BSplineCurve<D>::reversed() const {
    if (ctrl.empty()) return *this;
    const double a = domainBegin(), b = domainEnd();
    std::vector<double> k(knot.size());
    for (std::size_t i = 0; i < knot.size(); ++i) k[i] = (a + b) - knot[knot.size() - 1 - i];
    std::vector<Vec<D>> c(ctrl.rbegin(), ctrl.rend());
    return BSplineCurve(p, std::move(k), std::move(c));
}

// The NURBS Book, A5.1, with its indices: np is the index of the last control
// point and mp that of the last knot.
template <std::size_t D>
int BSplineCurve<D>::insertKnot(double u, int times) {
    if (ctrl.empty() || times < 1) return 0;
    if (!(u > domainBegin() && u < domainEnd())) {
        throw std::invalid_argument("BSplineCurve::insertKnot: u is not inside the domain");
    }
    const int np = static_cast<int>(ctrl.size()) - 1;
    const int mp = np + p + 1;
    const int k = findSpan(knot, p, u);
    int s = 0;
    for (int i = k; i >= 0 && knot[i] == u; --i) ++s;
    const int r = std::min(times, p - s);
    if (r <= 0) return 0;

    std::vector<double> UQ(knot.size() + r);
    for (int i = 0; i <= k; ++i) UQ[i] = knot[i];
    for (int i = 1; i <= r; ++i) UQ[k + i] = u;
    for (int i = k + 1; i <= mp; ++i) UQ[i + r] = knot[i];

    std::vector<Vec<D>> Q(ctrl.size() + r);
    for (int i = 0; i <= k - p; ++i) Q[i] = ctrl[i];
    for (int i = k - s; i <= np; ++i) Q[i + r] = ctrl[i];

    std::vector<Vec<D>> R(static_cast<std::size_t>(p - s) + 1);
    for (int i = 0; i <= p - s; ++i) R[i] = ctrl[k - p + i];
    int L = 0;
    for (int j = 1; j <= r; ++j) {
        L = k - p + j;
        for (int i = 0; i <= p - j - s; ++i) {
            const double alpha = (u - knot[L + i]) / (knot[i + k + 1] - knot[L + i]);
            R[i] = R[i + 1] * alpha + R[i] * (1.0 - alpha);
        }
        Q[L] = R[0];
        Q[k + r - j - s] = R[p - j - s];
    }
    for (int i = L + 1; i < k - s; ++i) Q[i] = R[i - L];

    knot = std::move(UQ);
    ctrl = std::move(Q);
    return r;
}

// ---------------------------------------------------------------------------
// BSplineSurface
// ---------------------------------------------------------------------------
template <std::size_t D>
BSplineSurface<D>::BSplineSurface(int degreeS, std::vector<double> knotsS,
                                  int degreeT, std::vector<double> knotsT,
                                  std::vector<Vec<D>> points)
    : pS(degreeS), pT(degreeT), knotS(std::move(knotsS)), knotT(std::move(knotsT)),
      net(std::move(points)) {
    nS = controlPointCount(knotS, pS);
    nT = controlPointCount(knotT, pT);
    if (nS < 1 || nT < 1 ||
        net.size() != static_cast<std::size_t>(nS) * static_cast<std::size_t>(nT)) {
        throw std::invalid_argument("BSplineSurface: the net is not countS x countT");
    }
    requireKnots(knotS, pS, static_cast<std::size_t>(nS), "BSplineSurface (s)");
    requireKnots(knotT, pT, static_cast<std::size_t>(nT), "BSplineSurface (t)");
}

template <std::size_t D>
Vec<D> BSplineSurface<D>::evaluate(double s, double t) const {
    if (net.empty()) return zero<D>();
    s = std::max(knotS[pS], std::min(knotS[nS], s));
    t = std::max(knotT[pT], std::min(knotT[nT], t));
    const int si = findSpan(knotS, pS, s), ti = findSpan(knotT, pT, t);
    double Ns[kMaxDegree + 1], Nt[kMaxDegree + 1];
    basisFunctions(knotS, pS, si, s, Ns);
    basisFunctions(knotT, pT, ti, t, Nt);
    Vec<D> out = zero<D>();
    for (int b = 0; b <= pT; ++b) {
        for (int a = 0; a <= pS; ++a) {
            out = out + net[static_cast<std::size_t>(ti - pT + b) * nS + (si - pS + a)] *
                            (Ns[a] * Nt[b]);
        }
    }
    return out;
}

template <std::size_t D>
BSplineCurve<D> BSplineSurface<D>::row(int j) const {
    std::vector<Vec<D>> c(net.begin() + static_cast<std::ptrdiff_t>(j) * nS,
                          net.begin() + static_cast<std::ptrdiff_t>(j + 1) * nS);
    return BSplineCurve<D>(pS, knotS, std::move(c));
}

template <std::size_t D>
BSplineCurve<D> BSplineSurface<D>::column(int i) const {
    std::vector<Vec<D>> c(static_cast<std::size_t>(nT));
    for (int j = 0; j < nT; ++j) c[j] = control(i, j);
    return BSplineCurve<D>(pT, knotT, std::move(c));
}

template class BSplineCurve<2>;
template class BSplineCurve<3>;
template class BSplineSurface<2>;
template class BSplineSurface<3>;

} // namespace geom
