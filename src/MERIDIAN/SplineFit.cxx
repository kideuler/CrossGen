#include "MERIDIAN/SplineFit.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace {

// Cholesky on a small dense SPD system, in place. The systems here are
// (segments + 1) square -- at most a handful of rows -- so nothing more is
// warranted, and a factorisation that fails says the fit was rank deficient,
// which the regularisation is there to prevent.
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

double pointSegment(const Point &p, const Point &a, const Point &b) {
    const Point d = b - a;
    const double dd = dotP(d, d);
    if (dd <= 0.0) return normP(p - a);
    double s = dotP(p - a, d) / dd;
    s = std::max(0.0, std::min(1.0, s));
    return normP(p - (a + d * s));
}

} // namespace

// ---------------------------------------------------------------------------
SplineFit::SplineFit(const Arrangement &arrangement)
    : SplineFit(arrangement, Options()) {}

SplineFit::SplineFit(const Arrangement &arrangement, const Options &opts)
    : arr(&arrangement), options(opts) {
    if (options.segments < 1) options.segments = 1;
    if (options.samples < 2) options.samples = 2;

    const Mesh &m = arr->getMesh();
    Point lo{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    Point hi{-lo[0], -lo[1]};
    for (const Point &p : m.vertices) {
        lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
        hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
    }
    modelExtent = normP(hi - lo);
    if (!(modelExtent > 0.0)) modelExtent = 1.0;

    buildKnots();
    fitArcs();
    buildPatches();
    check();
}

// ---------------------------------------------------------------------------
// buildKnots()
//
// One knot vector for every arc in the model. Clamped so that the first and
// last control points are the arc's two end nodes, and uniform inside so that
// the interior knots have single multiplicity -- which is what leaves the curve
// C2 across them, and what Stage 10's refinement preserves.
// ---------------------------------------------------------------------------
void SplineFit::buildKnots() {
    const int p = 3;
    const int n = controlPointsPerArc();
    const int s = options.segments;
    knot.assign(n + p + 1, 0.0);
    for (int i = 0; i <= p; ++i) knot[i] = 0.0;
    for (int i = 1; i < s; ++i) knot[p + i] = static_cast<double>(i) / s;
    for (int i = n; i < n + p + 1; ++i) knot[i] = 1.0;

    greville.assign(n, 0.0);
    for (int i = 0; i < n; ++i) {
        greville[i] = (knot[i + 1] + knot[i + 2] + knot[i + 3]) / 3.0;
    }
}

int SplineFit::findSpan(double u) const {
    const int p = 3;
    const int n = controlPointsPerArc();
    if (u >= knot[n]) return n - 1;
    if (u <= knot[p]) return p;
    int lo = p, hi = n, mid = (lo + hi) / 2;
    while (u < knot[mid] || u >= knot[mid + 1]) {
        if (u < knot[mid]) hi = mid; else lo = mid;
        mid = (lo + hi) / 2;
    }
    return mid;
}

void SplineFit::basisFuns(int span, double u, const std::vector<double> &U, double *N) {
    const int p = 3;
    double left[4], right[4];
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

Point SplineFit::evaluateNet(const Curve &c, double u) const {
    if (c.ctrl.empty()) return Point{0.0, 0.0};
    u = std::max(0.0, std::min(1.0, u));
    const int span = findSpan(u);
    double N[4];
    basisFuns(span, u, knot, N);
    Point out{0.0, 0.0};
    for (int k = 0; k <= 3; ++k) out = out + c.ctrl[span - 3 + k] * N[k];
    return out;
}

// An exact arc is walked along its own normalised chord length, which is the
// parameterisation the fit was taken in as well. Using the same one for both is
// what keeps the difference between them -- the correction the Coons blend gets
// -- as small as the two curves actually are apart, rather than inflating it
// with a reparameterisation.
Point SplineFit::evaluate(const Curve &c, double u) const {
    if (!c.exact || c.poly.size() < 2) return evaluateNet(c, u);
    u = std::max(0.0, std::min(1.0, u));
    if (u <= 0.0) return c.poly.front();
    if (u >= 1.0) return c.poly.back();
    const size_t k = static_cast<size_t>(
        std::lower_bound(c.cum.begin(), c.cum.end(), u) - c.cum.begin());
    if (k == 0) return c.poly.front();
    if (k >= c.poly.size()) return c.poly.back();
    const double seg = c.cum[k] - c.cum[k - 1];
    const double w = seg > 0.0 ? (u - c.cum[k - 1]) / seg : 0.0;
    return c.poly[k - 1] + (c.poly[k] - c.poly[k - 1]) * w;
}

// ---------------------------------------------------------------------------
// sideCorrection()
//
// What side `k` of a patch has to be moved by, at patch-frame parameter `w`, to
// put it back on the polyline it was carried from. Zero on a fitted side, and
// zero at both ends of an exact one because fitOne() pinned the end control
// points to the arc's end nodes.
//
// The mapping is the inverse of the one buildPatches() used to orient the four
// control polygons: sides 0 and 1 run with the (s, t) frame and sides 2 and 3
// against it, and either may additionally be traversed against the arc's own
// direction.
// ---------------------------------------------------------------------------
Point SplineFit::sideCorrection(const Patch &p, int k, double w) const {
    if (!p.exactSide[k]) return Point{0.0, 0.0};
    const Curve &c = fitted[p.side[k]];
    double u = (k < 2) ? w : 1.0 - w;
    if (!p.forward[k]) u = 1.0 - u;
    return evaluate(c, u) - evaluateNet(c, u);
}

Point SplineFit::evaluate(const Patch &p, double s, double t) const {
    const int n = controlPointsPerArc();
    if (static_cast<int>(p.net.size()) != n * n) return Point{0.0, 0.0};
    s = std::max(0.0, std::min(1.0, s));
    t = std::max(0.0, std::min(1.0, t));
    const int si = findSpan(s), ti = findSpan(t);
    double Ns[4], Nt[4];
    basisFuns(si, s, knot, Ns);
    basisFuns(ti, t, knot, Nt);
    Point out{0.0, 0.0};
    for (int b = 0; b <= 3; ++b) {
        for (int a = 0; a <= 3; ++a) {
            out = out + p.net[(ti - 3 + b) * n + (si - 3 + a)] * (Ns[a] * Nt[b]);
        }
    }
    // The transfinite correction of the header. No bilinear corner term: every
    // d_k is zero at both of its ends, so the term it would subtract is zero.
    if (p.corrected) {
        out = out + sideCorrection(p, 0, s) * (1.0 - t) + sideCorrection(p, 2, s) * t +
              sideCorrection(p, 3, t) * (1.0 - s) + sideCorrection(p, 1, t) * s;
    }
    return out;
}

// ---------------------------------------------------------------------------
// fitOne()
//
// Sec. 5's "chord-length parameterise, least-squares fit a cubic B-spline with
// a fixed uniform knot vector", with the two ends held rather than fitted. They
// are held because they are the nodes of the arrangement: two arcs meeting at a
// cone have to arrive at the same point or the patches around it do not close,
// and a least-squares fit that is free at the ends misses each of them by its
// own residual. Holding them costs two degrees of freedom and buys exact
// interpolation of every corner of the layout.
// ---------------------------------------------------------------------------
SplineFit::Curve SplineFit::fitOne(const std::vector<Point> &poly) const {
    const int n = controlPointsPerArc();
    Curve c;
    c.ctrl.assign(n, Point{0.0, 0.0});
    if (poly.size() < 2) {
        if (!poly.empty()) for (Point &q : c.ctrl) q = poly.front();
        return c;
    }

    // A short arc -- five or six triangles crossed -- carries fewer points than
    // the fit has control points, and the normal equations are then rank
    // deficient. Subdividing its segments is not inventing data: the arc *is*
    // the polyline, straight between its points, so a point interpolated along
    // one of its segments lies on it exactly. It only fills in equations that
    // the sampling of the triangulation happened not to provide.
    std::vector<Point> data;
    const int want = 4 * n;
    if (static_cast<int>(poly.size()) >= want || poly.size() < 2) {
        data = poly;
    } else {
        const int k = static_cast<int>(std::ceil(static_cast<double>(want) /
                                                 (poly.size() - 1)));
        data.push_back(poly.front());
        for (size_t i = 1; i < poly.size(); ++i) {
            for (int j = 1; j < k; ++j) {
                data.push_back(poly[i - 1] + (poly[i] - poly[i - 1]) *
                                                 (static_cast<double>(j) / k));
            }
            // The last one is the original point itself and not a + (b - a),
            // which is a rounding away from it. The two ends of an arc are the
            // nodes of the layout, and every patch that meets there has to
            // arrive at the same double, not at one an ulp away.
            data.push_back(poly[i]);
        }
    }

    // Chord length.
    std::vector<double> t(data.size(), 0.0);
    double total = 0.0;
    for (size_t i = 1; i < data.size(); ++i) {
        total += normP(data[i] - data[i - 1]);
        t[i] = total;
    }
    c.length = total;
    if (total > 0.0) for (double &x : t) x /= total;
    else for (size_t i = 0; i < t.size(); ++i) t[i] = static_cast<double>(i) / (t.size() - 1.0);
    c.samples = static_cast<int>(data.size());

    const Point &p0 = data.front();
    const Point &pn = data.back();
    c.ctrl.front() = p0;
    c.ctrl.back() = pn;

    const int free = n - 2;
    if (free <= 0) return c;
    c.underdetermined = static_cast<int>(data.size()) - 2 < free;

    std::vector<double> A(static_cast<size_t>(free) * free, 0.0);
    std::vector<double> rhs(static_cast<size_t>(free) * 2, 0.0);

    for (size_t k = 1; k + 1 < data.size(); ++k) {
        const int span = findSpan(t[k]);
        double N[4];
        basisFuns(span, t[k], knot, N);
        // The residual with the two held control points already subtracted.
        Point b = data[k];
        double row[32] = {0.0};
        for (int a = 0; a <= 3; ++a) {
            const int i = span - 3 + a;
            if (i == 0) b = b - p0 * N[a];
            else if (i == n - 1) b = b - pn * N[a];
            else row[i - 1] = N[a];
        }
        for (int i = 0; i < free; ++i) {
            if (row[i] == 0.0) continue;
            for (int j = 0; j < free; ++j) {
                if (row[j] == 0.0) continue;
                A[static_cast<size_t>(i) * free + j] += row[i] * row[j];
            }
            rhs[i] += row[i] * b[0];
            rhs[free + i] += row[i] * b[1];
        }
    }

    // Tikhonov towards the straight line between the ends, at the Greville
    // abscissae. See Options::regularisation.
    double maxDiag = 0.0;
    for (int i = 0; i < free; ++i) maxDiag = std::max(maxDiag, A[static_cast<size_t>(i) * free + i]);
    const double reg = std::max(options.regularisation * maxDiag,
                                maxDiag > 0.0 ? 0.0 : 1.0);
    for (int i = 0; i < free; ++i) {
        const double g = greville[i + 1];
        const Point prior = p0 * (1.0 - g) + pn * g;
        A[static_cast<size_t>(i) * free + i] += reg;
        rhs[i] += reg * prior[0];
        rhs[free + i] += reg * prior[1];
    }

    std::vector<double> chol = A;
    if (!choleskySolve(chol, free, rhs, 2)) {
        for (int i = 0; i < free; ++i) {
            const double g = greville[i + 1];
            c.ctrl[i + 1] = p0 * (1.0 - g) + pn * g;
        }
        c.underdetermined = true;
    } else {
        for (int i = 0; i < free; ++i) c.ctrl[i + 1] = Point{rhs[i], rhs[free + i]};
    }

    // The error the fit actually made, measured against the polyline it was
    // given and not against its own parameterisation: a point of the data is
    // compared with the nearest place on the sampled curve.
    std::vector<Point> sampled;
    const int ns = std::max(64, 8 * n);
    sampled.reserve(ns + 1);
    for (int i = 0; i <= ns; ++i) sampled.push_back(evaluate(c, static_cast<double>(i) / ns));
    double sum = 0.0;
    for (const Point &q : poly) {
        double best = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i + 1 < sampled.size(); ++i) {
            best = std::min(best, pointSegment(q, sampled[i], sampled[i + 1]));
        }
        c.maxDeviation = std::max(c.maxDeviation, best);
        sum += best * best;
    }
    c.rmsDeviation = poly.empty() ? 0.0 : std::sqrt(sum / poly.size());
    return c;
}

// ---------------------------------------------------------------------------
// fitArcs()
//
// Every arc gets its control points, because every arc is a side of some patch
// and the Coons net has to be built from something. What differs is whether
// those control points are also the *curve*: on a separatrix they are, and on
// an arc that came from the input they are not -- the polyline is kept and
// evaluate() returns it. See the header for why the two kinds are not the same
// kind of object.
// ---------------------------------------------------------------------------
void SplineFit::fitArcs() {
    const std::vector<Arrangement::Arc> &list = arr->getArcs();
    report.arcs = static_cast<int>(list.size());
    report.controlPointsPerArc = controlPointsPerArc();
    fitted.assign(list.size(), Curve());
    double sum = 0.0;
    long long count = 0;
    for (size_t a = 0; a < list.size(); ++a) {
        if (list[a].points.size() < 2) continue;
        fitted[a] = fitOne(list[a].points);
        fitted[a].arc = static_cast<int>(a);
        ++report.curves;
        if (fitted[a].underdetermined) ++report.underdetermined;

        const Arrangement::ArcKind kind = list[a].kind;
        const bool exact =
            (kind == Arrangement::ArcKind::Boundary && !options.fitBoundaryArcs) ||
            (kind == Arrangement::ArcKind::Interface && !options.fitInterfaceArcs);
        if (exact) {
            Curve &c = fitted[a];
            c.exact = true;
            c.poly = list[a].points;
            c.cum.assign(c.poly.size(), 0.0);
            for (size_t i = 1; i < c.poly.size(); ++i) {
                c.cum[i] = c.cum[i - 1] + normP(c.poly[i] - c.poly[i - 1]);
            }
            const double total = c.cum.back();
            if (total > 0.0) {
                for (double &x : c.cum) x /= total;
            } else {
                for (size_t i = 0; i < c.cum.size(); ++i) {
                    c.cum[i] = static_cast<double>(i) / (c.cum.size() - 1.0);
                }
            }
            c.cum.back() = 1.0;
            ++report.exactArcs;
            if (c.maxDeviation > report.maxNetDeviation) {
                report.maxNetDeviation = c.maxDeviation;
                report.worstNetArc = static_cast<int>(a);
            }
            continue;
        }

        ++report.fittedArcs;
        if (fitted[a].maxDeviation > report.maxDeviation) {
            report.maxDeviation = fitted[a].maxDeviation;
            report.worstArc = static_cast<int>(a);
        }
        sum += fitted[a].rmsDeviation * fitted[a].rmsDeviation * fitted[a].samples;
        count += fitted[a].samples;
    }
    report.rmsDeviation = count > 0 ? std::sqrt(sum / count) / modelExtent : 0.0;
    report.maxDeviation /= modelExtent;
    report.maxNetDeviation /= modelExtent;
}

// ---------------------------------------------------------------------------
// buildPatches()
//
// The Coons blend of Sec. 5, taken at the control points. The four sides come
// out of Stage 8 in cyclic order, so they are re-oriented once into the (s, t)
// frame of the patch:
//
//     corner 0 -> (0,0)   side 0 = C0, the bottom, s increasing
//     corner 1 -> (1,0)   side 1 = D1, the right,  t increasing
//     corner 2 -> (1,1)   side 2 reversed = C1, the top
//     corner 3 -> (0,1)   side 3 reversed = D0, the left
//
// Reversing a curve over a clamped *uniform* knot vector is exactly reversing
// its control points, because that knot vector is symmetric -- which is the
// second reason, after the shared knot vector, that every arc is given the same
// one. Nothing is refitted, so both patches sharing an arc get the same control
// points to the last bit.
// ---------------------------------------------------------------------------
void SplineFit::buildPatches() {
    const int n = controlPointsPerArc();
    const std::vector<int> &pf = arr->patchFaces();
    report.faces = static_cast<int>(pf.size());

    for (int f : pf) {
        const std::vector<Arrangement::Side> sides = arr->patchSides(f);
        if (sides.size() != 4) { ++report.skipped; continue; }
        bool ok = true;
        for (const Arrangement::Side &sd : sides) {
            if (sd.arc < 0 || static_cast<int>(fitted[sd.arc].ctrl.size()) != n) ok = false;
        }
        if (!ok) { ++report.skipped; continue; }

        Patch p;
        p.face = f;
        p.faceArea = arr->getFaces()[f].area;
        for (int k = 0; k < 4; ++k) {
            p.side[k] = sides[k].arc;
            p.forward[k] = sides[k].forward;
            const Arrangement::Arc &ar = arr->getArcs()[sides[k].arc];
            p.corner[k] = sides[k].forward ? ar.from : ar.to;
            p.exactSide[k] = fitted[sides[k].arc].exact;
            if (p.exactSide[k]) p.corrected = true;
        }
        if (p.corrected) ++report.correctedPatches;

        // The four sides as control polygons, each running the way the (s, t)
        // frame wants it.
        auto oriented = [&](int k, bool flip) {
            std::vector<Point> c = fitted[p.side[k]].ctrl;
            if (p.forward[k] == flip) std::reverse(c.begin(), c.end());
            return c;
        };
        const std::vector<Point> C0 = oriented(0, false);  // (0,0) -> (1,0)
        const std::vector<Point> D1 = oriented(1, false);  // (1,0) -> (1,1)
        const std::vector<Point> C1 = oriented(2, true);   // (0,1) -> (1,1)
        const std::vector<Point> D0 = oriented(3, true);   // (0,0) -> (0,1)

        const Point p00 = C0.front(), p10 = C0.back();
        const Point p01 = C1.front(), p11 = C1.back();

        p.net.assign(static_cast<size_t>(n) * n, Point{0.0, 0.0});
        for (int j = 0; j < n; ++j) {
            const double tt = greville[j];
            for (int i = 0; i < n; ++i) {
                const double ss = greville[i];
                Point q = C0[i] * (1.0 - tt) + C1[i] * tt + D0[j] * (1.0 - ss) + D1[j] * ss;
                q = q - (p00 * ((1.0 - ss) * (1.0 - tt)) + p10 * (ss * (1.0 - tt)) +
                         p01 * ((1.0 - ss) * tt) + p11 * (ss * tt));
                p.net[static_cast<size_t>(j) * n + i] = q;
            }
        }
        // The boundary rows have to come out equal to the side fits, which is
        // the watertightness rule; put them there literally so that no rounding
        // of the blend can separate two patches by an ulp.
        for (int i = 0; i < n; ++i) {
            p.net[i] = C0[i];
            p.net[static_cast<size_t>(n - 1) * n + i] = C1[i];
        }
        for (int j = 0; j < n; ++j) {
            p.net[static_cast<size_t>(j) * n] = D0[j];
            p.net[static_cast<size_t>(j) * n + (n - 1)] = D1[j];
        }

        // Area and the worst sampled cell, which is what says the blend did not
        // fold the patch over.
        const int m = options.samples;
        std::vector<Point> grid(static_cast<size_t>(m + 1) * (m + 1));
        for (int b = 0; b <= m; ++b) {
            for (int a = 0; a <= m; ++a) {
                grid[static_cast<size_t>(b) * (m + 1) + a] =
                    evaluate(p, static_cast<double>(a) / m, static_cast<double>(b) / m);
            }
        }
        double area = 0.0, worst = std::numeric_limits<double>::infinity();
        for (int b = 0; b < m; ++b) {
            for (int a = 0; a < m; ++a) {
                const Point &q00 = grid[static_cast<size_t>(b) * (m + 1) + a];
                const Point &q10 = grid[static_cast<size_t>(b) * (m + 1) + a + 1];
                const Point &q11 = grid[static_cast<size_t>(b + 1) * (m + 1) + a + 1];
                const Point &q01 = grid[static_cast<size_t>(b + 1) * (m + 1) + a];
                const double cell = 0.5 * (cross2(q10 - q00, q11 - q00) +
                                           cross2(q11 - q00, q01 - q00));
                area += cell;
                worst = std::min(worst, cell);
            }
        }
        p.area = area;
        const double mean = area / (m * m);
        p.minCellRatio = (std::fabs(mean) > 0.0) ? worst / mean : 0.0;
        if (p.minCellRatio <= 0.0) ++report.foldedPatches;

        nets.push_back(std::move(p));
    }
    report.patches = static_cast<int>(nets.size());
}

// ---------------------------------------------------------------------------
// check()
//
// The watertightness rule, read back off the result. Two patches that share an
// arc must carry the same control points for it -- not close, the same -- and
// every patch corner must sit on the node of the arrangement it was built from.
// Both are zero when the construction above is right, so anything else is a
// bug rather than a tolerance.
// ---------------------------------------------------------------------------
void SplineFit::check() {
    const int n = controlPointsPerArc();
    std::vector<std::vector<int>> byArc(arr->getArcs().size());
    for (size_t k = 0; k < nets.size(); ++k) {
        for (int s = 0; s < 4; ++s) byArc[nets[k].side[s]].push_back(static_cast<int>(k));
    }

    for (size_t a = 0; a < byArc.size(); ++a) {
        if (byArc[a].size() < 2) continue;
        ++report.sharedArcs;
        // The arc's control points as each patch holds them, read off the
        // boundary row it lives in and put back into the arc's own direction.
        std::vector<std::vector<Point>> seen;
        for (int k : byArc[a]) {
            const Patch &p = nets[k];
            for (int s = 0; s < 4; ++s) {
                if (p.side[s] != static_cast<int>(a)) continue;
                std::vector<Point> row(n);
                for (int i = 0; i < n; ++i) {
                    switch (s) {
                        case 0: row[i] = p.net[i]; break;
                        case 1: row[i] = p.net[static_cast<size_t>(i) * n + (n - 1)]; break;
                        case 2: row[i] = p.net[static_cast<size_t>(n - 1) * n + (n - 1 - i)]; break;
                        default: row[i] = p.net[static_cast<size_t>(n - 1 - i) * n]; break;
                    }
                }
                if (!p.forward[s]) std::reverse(row.begin(), row.end());
                seen.push_back(std::move(row));
            }
        }
        for (size_t i = 1; i < seen.size(); ++i) {
            for (int c = 0; c < n; ++c) {
                report.maxSeamGap = std::max(report.maxSeamGap,
                                             normP(seen[i][c] - seen[0][c]) / modelExtent);
            }
        }
    }

    report.minCellRatio = nets.empty() ? 0.0 : std::numeric_limits<double>::infinity();
    for (const Patch &p : nets) {
        report.patchArea += p.area;
        report.faceArea += p.faceArea;
        report.minCellRatio = std::min(report.minCellRatio, p.minCellRatio);
        for (int k = 0; k < 4; ++k) {
            const Point corner = (k == 0) ? p.net[0]
                               : (k == 1) ? p.net[n - 1]
                               : (k == 2) ? p.net[static_cast<size_t>(n) * n - 1]
                                          : p.net[static_cast<size_t>(n - 1) * n];
            report.maxCornerGap = std::max(report.maxCornerGap,
                                           normP(corner - arr->getNodes()[p.corner[k]].p) /
                                               modelExtent);
        }
    }
    if (!std::isfinite(report.minCellRatio)) report.minCellRatio = 0.0;

    // Does each patch actually lie on its four arcs? For a pure bicubic that is
    // the construction and nothing more can be said; with the transfinite
    // correction in play it is a claim about arithmetic, so it is sampled.
    const int ms = std::max(8, options.samples * 4);
    for (const Patch &p : nets) {
        for (int k = 0; k < 4; ++k) {
            const Curve &c = fitted[p.side[k]];
            if (c.ctrl.empty()) continue;
            for (int i = 0; i <= ms; ++i) {
                const double w = static_cast<double>(i) / ms;
                double u = (k < 2) ? w : 1.0 - w;
                if (!p.forward[k]) u = 1.0 - u;
                const Point on = (k == 0) ? evaluate(p, w, 0.0)
                               : (k == 1) ? evaluate(p, 1.0, w)
                               : (k == 2) ? evaluate(p, w, 1.0)
                                          : evaluate(p, 0.0, w);
                report.maxBoundaryGap = std::max(report.maxBoundaryGap,
                                                 normP(on - evaluate(c, u)) / modelExtent);
            }
        }
    }

    report.watertight = report.maxSeamGap == 0.0 && report.maxCornerGap == 0.0 &&
                        report.maxBoundaryGap <= Report::boundaryTolerance;

    std::ostringstream oss;
    if (report.skipped > 0) {
        oss << report.skipped << " of " << report.faces << " face(s) were not fitted: a Coons "
            << "patch needs four sides of one arc each, and Stage 8 did not give them one. "
            << "Sec. 4's remedy is upstream -- raise lambda_5 and re-run Stage 6 from the "
            << "current phi.";
        report.messages.push_back(oss.str());
        oss.str("");
    }
    if (report.foldedPatches > 0) {
        oss << report.foldedPatches << " patch(es) fold: the bilinear blend of their four "
            << "sides reverses somewhere inside. The layout face is too far from a rectangle "
            << "for a Coons patch to cover it, which is Sec. 5's own limitation and needs "
            << "more patches rather than a better fit.";
        report.messages.push_back(oss.str());
        oss.str("");
    }
    if (report.underdetermined > 0) {
        oss << report.underdetermined << " arc(s) carry fewer sample points than the "
            << report.controlPointsPerArc << " control points asked of them, so the "
            << "regularisation set their interior control points rather than the data.";
        report.messages.push_back(oss.str());
        oss.str("");
    }
    if (!report.watertight) {
        oss << "Two patches disagree about a shared arc by " << report.maxSeamGap
            << " of the model, a corner is " << report.maxCornerGap << " off its node, or a "
            << "patch leaves its own arc by " << report.maxBoundaryGap
            << ". The first two should be exactly zero and the third within "
            << Report::boundaryTolerance << ".";
        report.messages.push_back(oss.str());
        oss.str("");
    }

    report.valid = report.patches > 0 && report.skipped == 0 && report.foldedPatches == 0 &&
                   report.watertight;
}

// ---------------------------------------------------------------------------
bool SplineFit::writeCurvesOBJ(const std::string &filename, int samples) const {
    std::ofstream out(filename);
    if (!out) return false;
    if (samples < 2) samples = 2;
    int base = 1;
    for (const Curve &c : fitted) {
        if (c.ctrl.empty()) continue;
        // An exact arc is written as the polyline it is, not as that polyline
        // resampled: resampling it at `samples` stations would cut its corners
        // for exactly the reason the fit was taken off it.
        if (c.exact && c.poly.size() >= 2) {
            for (const Point &q : c.poly) out << "v " << q[0] << " " << q[1] << " 0\n";
            out << "l";
            for (size_t i = 0; i < c.poly.size(); ++i) out << " " << (base + static_cast<int>(i));
            out << "\n";
            base += static_cast<int>(c.poly.size());
            continue;
        }
        for (int i = 0; i <= samples; ++i) {
            const Point p = evaluate(c, static_cast<double>(i) / samples);
            out << "v " << p[0] << " " << p[1] << " 0\n";
        }
        out << "l";
        for (int i = 0; i <= samples; ++i) out << " " << (base + i);
        out << "\n";
        base += samples + 1;
    }
    return true;
}

bool SplineFit::writeNetOBJ(const std::string &filename) const {
    std::ofstream out(filename);
    if (!out) return false;
    const int n = controlPointsPerArc();
    int base = 1;
    for (size_t k = 0; k < nets.size(); ++k) {
        const Patch &p = nets[k];
        out << "o patch" << k << "_net\n";
        for (const Point &q : p.net) out << "v " << q[0] << " " << q[1] << " 0\n";
        for (int j = 0; j < n; ++j) {
            out << "l";
            for (int i = 0; i < n; ++i) out << " " << (base + j * n + i);
            out << "\n";
        }
        for (int i = 0; i < n; ++i) {
            out << "l";
            for (int j = 0; j < n; ++j) out << " " << (base + j * n + i);
            out << "\n";
        }
        base += n * n;
    }
    return true;
}

bool SplineFit::writeSurfaceOBJ(const std::string &filename, int samples) const {
    std::ofstream out(filename);
    if (!out) return false;
    const int m = samples >= 2 ? samples : options.samples;
    int base = 1;
    for (size_t k = 0; k < nets.size(); ++k) {
        const Patch &p = nets[k];
        out << "o patch" << k << "\n";
        for (int b = 0; b <= m; ++b) {
            for (int a = 0; a <= m; ++a) {
                const Point q = evaluate(p, static_cast<double>(a) / m,
                                         static_cast<double>(b) / m);
                out << "v " << q[0] << " " << q[1] << " 0\n";
            }
        }
        for (int b = 0; b < m; ++b) {
            for (int a = 0; a < m; ++a) {
                const int v00 = base + b * (m + 1) + a;
                out << "f " << v00 << " " << (v00 + 1) << " " << (v00 + m + 2) << " "
                    << (v00 + m + 1) << "\n";
            }
        }
        base += (m + 1) * (m + 1);
    }
    return true;
}
