#ifndef __GEOM_VEC_HXX__
#define __GEOM_VEC_HXX__

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

// The coordinate type every class in src/geom is written over.
//
// geom::Vec<2> is std::array<double, 2>, which is exactly the global Point of
// mesh/Mesh.hxx, so a std::vector<Point> goes into a spline and comes back out
// without a copy or a conversion, and nothing in this directory has to include
// the mesh library to get it. Vec<3> is the same thing one dimension up, for
// the day a model carries positions in R^3.
//
// The arithmetic lives in namespace geom so that it never meets the global
// Point operators in overload resolution: code outside the namespace keeps
// using Mesh.hxx's, code inside uses these, and the two compute the same
// doubles.
//
// ### Why dot<2> is spelled out
//
// This toolchain contracts a*b + c*d into a fused multiply-add, and which
// product gets rounded first depends on how the expression is written. A loop
// accumulating s += a[i]*b[i] rounds the *first* product and fuses the second;
// Mesh.hxx's dotP writes a[0]*b[0] + a[1]*b[1] and rounds the *second*. The two
// differ in the last bit, and a spline fitted through lengths that differ in
// the last bit is a different spline. The specialisation keeps geom::norm of a
// Point equal to normP of it, bit for bit, which is what lets the MERIDIAN and
// TORSION corpora come out unchanged on top of this code.
namespace geom {

template <std::size_t D>
using Vec = std::array<double, D>;

template <std::size_t D>
inline Vec<D> operator+(const Vec<D> &a, const Vec<D> &b) {
    Vec<D> r;
    for (std::size_t i = 0; i < D; ++i) r[i] = a[i] + b[i];
    return r;
}

template <std::size_t D>
inline Vec<D> operator-(const Vec<D> &a, const Vec<D> &b) {
    Vec<D> r;
    for (std::size_t i = 0; i < D; ++i) r[i] = a[i] - b[i];
    return r;
}

template <std::size_t D>
inline Vec<D> operator*(const Vec<D> &a, double s) {
    Vec<D> r;
    for (std::size_t i = 0; i < D; ++i) r[i] = a[i] * s;
    return r;
}

template <std::size_t D>
inline Vec<D> operator*(double s, const Vec<D> &a) {
    return a * s;
}

template <std::size_t D>
inline Vec<D> zero() {
    Vec<D> r;
    r.fill(0.0);
    return r;
}

template <std::size_t D>
inline double dot(const Vec<D> &a, const Vec<D> &b) {
    double s = 0.0;
    for (std::size_t i = 0; i < D; ++i) s += a[i] * b[i];
    return s;
}

// See the header: the same expression as dotP, not the loop.
template <>
inline double dot<2>(const Vec<2> &a, const Vec<2> &b) {
    return a[0] * b[0] + a[1] * b[1];
}

template <std::size_t D>
inline double norm(const Vec<D> &a) {
    return std::sqrt(dot(a, a));
}

template <std::size_t D>
inline double distance(const Vec<D> &a, const Vec<D> &b) {
    return norm(a - b);
}

// Distance from p to the closed segment [a, b].
template <std::size_t D>
inline double distanceToSegment(const Vec<D> &p, const Vec<D> &a, const Vec<D> &b) {
    const Vec<D> d = b - a;
    const double dd = dot(d, d);
    if (dd <= 0.0) return norm(p - a);
    double s = dot(p - a, d) / dd;
    s = std::max(0.0, std::min(1.0, s));
    return norm(p - (a + d * s));
}

} // namespace geom

#endif // __GEOM_VEC_HXX__
