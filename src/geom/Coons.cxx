#include "geom/Coons.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include "geom/detail/Occ.hxx"

namespace geom {

using detail::Access;

template <std::size_t D>
std::vector<Vec<D>> coonsNet(const std::vector<Vec<D>> &bottom, const std::vector<Vec<D>> &top,
                             const std::vector<Vec<D>> &left, const std::vector<Vec<D>> &right,
                             const std::vector<double> &grevilleS,
                             const std::vector<double> &grevilleT) {
    const std::size_t nS = bottom.size(), nT = left.size();
    if (nS < 2 || nT < 2 || top.size() != nS || right.size() != nT ||
        grevilleS.size() != nS || grevilleT.size() != nT) {
        throw std::invalid_argument("coonsNet: opposite sides or their abscissae differ in size");
    }
    const Vec<D> p00 = bottom.front(), p10 = bottom.back();
    const Vec<D> p01 = top.front(), p11 = top.back();

    std::vector<Vec<D>> net(nS * nT, zero<D>());
    for (std::size_t j = 0; j < nT; ++j) {
        const double tt = grevilleT[j];
        for (std::size_t i = 0; i < nS; ++i) {
            const double ss = grevilleS[i];
            Vec<D> q = bottom[i] * (1.0 - tt) + top[i] * tt + left[j] * (1.0 - ss) + right[j] * ss;
            q = q - (p00 * ((1.0 - ss) * (1.0 - tt)) + p10 * (ss * (1.0 - tt)) +
                     p01 * ((1.0 - ss) * tt) + p11 * (ss * tt));
            net[j * nS + i] = q;
        }
    }
    // The boundary rows are the sides, put there literally so that no rounding
    // of the blend can separate two patches that share one by an ulp.
    for (std::size_t i = 0; i < nS; ++i) {
        net[i] = bottom[i];
        net[(nT - 1) * nS + i] = top[i];
    }
    for (std::size_t j = 0; j < nT; ++j) {
        net[j * nS] = left[j];
        net[j * nS + (nS - 1)] = right[j];
    }
    return net;
}

// ---------------------------------------------------------------------------
// makeCompatible()
//
// Degree elevation and knot insertion, both of which the kernel does without
// moving the curve. The one step that is not exact is the snap: two knots a few
// ulps apart would otherwise both be inserted, into a span the kernel refuses
// as too short, so the second curve's knot is moved onto the first curve's.
// That moves the second curve by the same few ulps of its parameter.
// ---------------------------------------------------------------------------
template <std::size_t D>
std::pair<BSplineCurve<D>, BSplineCurve<D>> makeCompatible(const BSplineCurve<D> &a,
                                                           const BSplineCurve<D> &b) {
    if (a.empty() || b.empty()) throw std::invalid_argument("makeCompatible: an empty curve");
    Handle(Geom_BSplineCurve) A = detail::copyOf(Access::curve(a));
    Handle(Geom_BSplineCurve) B = detail::copyOf(Access::curve(b));

    const double a0 = A->FirstParameter(), a1 = A->LastParameter();
    const double snap = 16.0 * std::numeric_limits<double>::epsilon() *
                        std::max({1.0, std::fabs(a0), std::fabs(a1)});
    if (std::fabs(B->FirstParameter() - a0) > snap || std::fabs(B->LastParameter() - a1) > snap) {
        throw std::invalid_argument("makeCompatible: the curves have different parameter domains");
    }

    detail::guard("makeCompatible", [&] {
        const int p = std::max(A->Degree(), B->Degree());
        if (A->Degree() < p) A->IncreaseDegree(p);
        if (B->Degree() < p) B->IncreaseDegree(p);

        // The ends first, so that the domains are the same doubles.
        if (B->Knot(1) != A->Knot(1)) B->SetKnot(1, A->Knot(1));
        if (B->Knot(B->NbKnots()) != A->Knot(A->NbKnots())) B->SetKnot(B->NbKnots(), A->Knot(A->NbKnots()));
        for (int i = 2; i < B->NbKnots(); ++i) {
            for (int k = 2; k < A->NbKnots(); ++k) {
                const double ka = A->Knot(k);
                if (ka != B->Knot(i) && std::fabs(ka - B->Knot(i)) <= snap) {
                    B->SetKnot(i, ka);
                    break;
                }
            }
        }

        auto interior = [](const Handle(Geom_BSplineCurve) &c, TColStd_Array1OfReal &knots,
                           TColStd_Array1OfInteger &mults) {
            const int n = c->NbKnots() - 2;
            if (n <= 0) return false;
            knots.Resize(1, n, false);
            mults.Resize(1, n, false);
            for (int i = 0; i < n; ++i) {
                knots(i + 1) = c->Knot(i + 2);
                mults(i + 1) = c->Multiplicity(i + 2);
            }
            return true;
        };
        TColStd_Array1OfReal ka(1, 1), kb(1, 1);
        TColStd_Array1OfInteger ma(1, 1), mb(1, 1);
        const bool hasA = interior(A, ka, ma);
        const bool hasB = interior(B, kb, mb);
        if (hasB) A->InsertKnots(kb, mb, 0.0, Standard_False);
        if (hasA) B->InsertKnots(ka, ma, 0.0, Standard_False);
    });

    BSplineCurve<D> ca = Access::curve<D>(A), cb = Access::curve<D>(B);
    if (ca.degree() != cb.degree() || ca.knots() != cb.knots()) {
        throw std::runtime_error("makeCompatible: the kernel left the knot vectors different");
    }
    return {std::move(ca), std::move(cb)};
}

template <std::size_t D>
BSplineSurface<D> coonsSurface(const BSplineCurve<D> &bottom, const BSplineCurve<D> &top,
                               const BSplineCurve<D> &left, const BSplineCurve<D> &right) {
    if (bottom.empty() || top.empty() || left.empty() || right.empty()) {
        throw std::invalid_argument("coonsSurface: an empty side");
    }
    auto pair = [](const BSplineCurve<D> &x, const BSplineCurve<D> &y) {
        if (x.degree() == y.degree() && x.knots() == y.knots()) return std::make_pair(x, y);
        return makeCompatible(x, y);
    };
    const auto [B, T] = pair(bottom, top);
    const auto [L, R] = pair(left, right);

    // The Grevilles, carried onto [0, 1] for the blend. On a curve already over
    // [0, 1] that is the identity, to the bit.
    auto unit = [](const BSplineCurve<D> &c) {
        std::vector<double> g = grevilleAbscissae(c.knots(), c.degree());
        const double a = c.domainBegin(), b = c.domainEnd();
        if (a != 0.0 || b != 1.0) {
            for (double &x : g) x = (x - a) / (b - a);
        }
        return g;
    };
    const std::vector<Vec<D>> net = coonsNet(B.controlPoints(), T.controlPoints(),
                                             L.controlPoints(), R.controlPoints(), unit(B), unit(L));
    return BSplineSurface<D>(B.degree(), B.knots(), L.degree(), L.knots(), net);
}

template std::vector<Vec<2>> coonsNet<2>(const std::vector<Vec<2>> &, const std::vector<Vec<2>> &,
                                         const std::vector<Vec<2>> &, const std::vector<Vec<2>> &,
                                         const std::vector<double> &, const std::vector<double> &);
template std::vector<Vec<3>> coonsNet<3>(const std::vector<Vec<3>> &, const std::vector<Vec<3>> &,
                                         const std::vector<Vec<3>> &, const std::vector<Vec<3>> &,
                                         const std::vector<double> &, const std::vector<double> &);
template std::pair<BSplineCurve<2>, BSplineCurve<2>> makeCompatible<2>(const BSplineCurve<2> &,
                                                                      const BSplineCurve<2> &);
template std::pair<BSplineCurve<3>, BSplineCurve<3>> makeCompatible<3>(const BSplineCurve<3> &,
                                                                      const BSplineCurve<3> &);
template BSplineSurface<2> coonsSurface<2>(const BSplineCurve<2> &, const BSplineCurve<2> &,
                                           const BSplineCurve<2> &, const BSplineCurve<2> &);
template BSplineSurface<3> coonsSurface<3>(const BSplineCurve<3> &, const BSplineCurve<3> &,
                                           const BSplineCurve<3> &, const BSplineCurve<3> &);

} // namespace geom
