#include "geom/Coons.hxx"

#include <stdexcept>
#include <utility>

namespace geom {

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

template <std::size_t D>
BSplineSurface<D> coonsSurface(const BSplineCurve<D> &bottom, const BSplineCurve<D> &top,
                               const BSplineCurve<D> &left, const BSplineCurve<D> &right) {
    if (bottom.degree() != top.degree() || bottom.knots() != top.knots() ||
        left.degree() != right.degree() || left.knots() != right.knots()) {
        throw std::invalid_argument("coonsSurface: opposite sides do not share a knot vector");
    }
    std::vector<Vec<D>> net =
        coonsNet(bottom.controlPoints(), top.controlPoints(), left.controlPoints(),
                 right.controlPoints(), grevilleAbscissae(bottom.knots(), bottom.degree()),
                 grevilleAbscissae(left.knots(), left.degree()));
    return BSplineSurface<D>(bottom.degree(), bottom.knots(), left.degree(), left.knots(),
                             std::move(net));
}

template std::vector<Vec<2>> coonsNet<2>(const std::vector<Vec<2>> &, const std::vector<Vec<2>> &,
                                         const std::vector<Vec<2>> &, const std::vector<Vec<2>> &,
                                         const std::vector<double> &, const std::vector<double> &);
template std::vector<Vec<3>> coonsNet<3>(const std::vector<Vec<3>> &, const std::vector<Vec<3>> &,
                                         const std::vector<Vec<3>> &, const std::vector<Vec<3>> &,
                                         const std::vector<double> &, const std::vector<double> &);
template BSplineSurface<2> coonsSurface<2>(const BSplineCurve<2> &, const BSplineCurve<2> &,
                                           const BSplineCurve<2> &, const BSplineCurve<2> &);
template BSplineSurface<3> coonsSurface<3>(const BSplineCurve<3> &, const BSplineCurve<3> &,
                                           const BSplineCurve<3> &, const BSplineCurve<3> &);

} // namespace geom
