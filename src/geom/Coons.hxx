#ifndef __GEOM_COONS_HXX__
#define __GEOM_COONS_HXX__

#include <utility>
#include <vector>

#include "geom/BSpline.hxx"
#include "geom/Vec.hxx"

// Bilinearly blended Coons patches: the surface over the unit square that
// interpolates four boundary curves,
//
//     S(s,t) = (1-t) B(s) + t T(s) + (1-s) L(t) + s R(t)
//              - [ (1-s)(1-t) P00 + s(1-t) P10 + (1-s)t P01 + s t P11 ]
//
// with B the bottom, T the top, L the left and R the right side, and the P the
// four corners. Every function here takes the sides in that order and in the
// same orientation -- *not* cyclic:
//
//     bottom  (0,0) -> (1,0)      top    (0,1) -> (1,1)
//     left    (0,0) -> (0,1)      right  (1,0) -> (1,1)
//
// and takes the corners from the bottom and top sides. The left and right
// sides are expected to meet them; nothing checks that they do.
//
// Two forms. coonsPoint() blends four boundary *points* for one (s, t), which is
// a transfinite grid on sampled sides. coonsNet() blends four control polygons
// at their Greville abscissae, which gives the B-spline surface whose boundary
// rows are the four side curves identically -- not to a tolerance -- so that
// two patches built from one shared curve meet exactly. That is the
// construction of SplineFit (Shepherd, Gu and Hughes 2022, Sec. 5). Taking the
// blend at the Grevilles rather than at i / (n - 1) is what lets the bilinear
// term reproduce a linear function exactly: four straight sides give the
// bilinear patch, not a slightly bowed one. It is also why this is not
// OpenCASCADE's GeomFill Coons filling, which blends at the pole indices and
// misses the bilinear patch on straight B-spline sides by 1.5e-2 of its size.
//
// coonsSurface() is the curve-level form, and by that linear precision it is
// the Coons patch of the four curves *as curves*, not merely of their control
// polygons: sum_j N_j(t) (1 - t_j) = 1 - t with t_j the Grevilles, so the
// tensor-product surface over the blended net is (1 - t) B(s) + t T(s) + ...
// identically. Sides that do not share a knot vector are first given one
// (makeCompatible), which changes neither of them, so any four curves meeting
// at their corners have a Coons surface -- a polyline side and a cubic side
// included.
namespace geom {

template <std::size_t D>
inline Vec<D> coonsPoint(const Vec<D> &bottom, const Vec<D> &top,
                         const Vec<D> &left, const Vec<D> &right,
                         const Vec<D> &p00, const Vec<D> &p10,
                         const Vec<D> &p01, const Vec<D> &p11,
                         double s, double t) {
    return bottom * (1.0 - t) + top * t + left * (1.0 - s) + right * s -
           (p00 * ((1.0 - s) * (1.0 - t)) + p10 * (s * (1.0 - t)) +
            p01 * ((1.0 - s) * t) + p11 * (s * t));
}

// The control net, row-major (net[j * nS + i] at (s_i, t_j)), of the Coons
// blend of four control polygons. bottom and top carry nS points and are blended
// at grevilleS; left and right carry nT and are blended at grevilleT. The
// boundary rows of the result are copies of the four inputs. Throws
// std::invalid_argument on mismatched sizes.
template <std::size_t D>
std::vector<Vec<D>> coonsNet(const std::vector<Vec<D>> &bottom, const std::vector<Vec<D>> &top,
                             const std::vector<Vec<D>> &left, const std::vector<Vec<D>> &right,
                             const std::vector<double> &grevilleS,
                             const std::vector<double> &grevilleT);

// Two curves over one parameter domain, re-expressed over a common degree and
// knot vector without changing either: the lower degree is raised and each
// curve is given the other's knots, by the kernel. Knots of the two that are
// within a few ulps of each other are first made the same double, which is what
// lets a curve and its mirrored neighbour (1 - 2/3 against 1/3) agree. Throws
// std::invalid_argument if either is empty or the domains differ.
template <std::size_t D>
std::pair<BSplineCurve<D>, BSplineCurve<D>> makeCompatible(const BSplineCurve<D> &a,
                                                           const BSplineCurve<D> &b);

// The Coons surface of four curves, in the orientation above. Opposite sides
// that do not already share a degree and a knot vector are made to
// (makeCompatible); where they do, the net is exactly coonsNet() of their
// control points. Each pair must share a parameter domain, which is mapped
// onto [0, 1] for the blend.
template <std::size_t D>
BSplineSurface<D> coonsSurface(const BSplineCurve<D> &bottom, const BSplineCurve<D> &top,
                               const BSplineCurve<D> &left, const BSplineCurve<D> &right);

extern template std::vector<Vec<2>> coonsNet<2>(const std::vector<Vec<2>> &, const std::vector<Vec<2>> &,
                                                const std::vector<Vec<2>> &, const std::vector<Vec<2>> &,
                                                const std::vector<double> &, const std::vector<double> &);
extern template std::vector<Vec<3>> coonsNet<3>(const std::vector<Vec<3>> &, const std::vector<Vec<3>> &,
                                                const std::vector<Vec<3>> &, const std::vector<Vec<3>> &,
                                                const std::vector<double> &, const std::vector<double> &);
extern template std::pair<BSplineCurve<2>, BSplineCurve<2>> makeCompatible<2>(const BSplineCurve<2> &,
                                                                             const BSplineCurve<2> &);
extern template std::pair<BSplineCurve<3>, BSplineCurve<3>> makeCompatible<3>(const BSplineCurve<3> &,
                                                                             const BSplineCurve<3> &);
extern template BSplineSurface<2> coonsSurface<2>(const BSplineCurve<2> &, const BSplineCurve<2> &,
                                                  const BSplineCurve<2> &, const BSplineCurve<2> &);
extern template BSplineSurface<3> coonsSurface<3>(const BSplineCurve<3> &, const BSplineCurve<3> &,
                                                  const BSplineCurve<3> &, const BSplineCurve<3> &);

} // namespace geom

#endif // __GEOM_COONS_HXX__
