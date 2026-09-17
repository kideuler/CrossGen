#include "MERIDIAN/SplineFit.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "geom/Coons.hxx"
#include "geom/Fitting.hxx"

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

    // One knot vector for every arc in the model. Clamped so that the first and
    // last control points are the arc's two end nodes, and uniform inside so
    // that the interior knots have single multiplicity -- which is what leaves
    // the curve C2 across them, and what Stage 10's refinement preserves.
    knot = geom::clampedUniformKnots(3, options.segments);
    greville = geom::grevilleAbscissae(knot, 3);

    fitArcs();
    buildPatches();
    check();
}

// An exact arc is walked along its own normalised chord length, which is the
// parameterisation the fit was taken in as well. Using the same one for both is
// what keeps the difference between them -- the correction the Coons blend gets
// -- as small as the two curves actually are apart, rather than inflating it
// with a reparameterisation.
Point SplineFit::evaluate(const Curve &c, double u) const {
    if (!c.exact || c.poly.size() < 2) return c.spline.evaluate(u);
    return c.poly.evaluate(u);
}

// ---------------------------------------------------------------------------
// sideCorrection()
//
// What side `k` of a patch has to be moved by, at patch-frame parameter `w`, to
// put it back on the polyline it was carried from. Zero on a fitted side, and
// zero at both ends of an exact one because the fit pinned the end control
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
    return evaluate(c, u) - c.spline.evaluate(u);
}

Point SplineFit::evaluate(const Patch &p, double s, double t) const {
    const int n = controlPointsPerArc();
    if (p.surface.countS() != n || p.surface.countT() != n) return Point{0.0, 0.0};
    s = std::max(0.0, std::min(1.0, s));
    t = std::max(0.0, std::min(1.0, t));
    Point out = p.surface.evaluate(s, t);
    // The transfinite correction of the header. No bilinear corner term: every
    // d_k is zero at both of its ends, so the term it would subtract is zero.
    if (p.corrected) {
        out = out + sideCorrection(p, 0, s) * (1.0 - t) + sideCorrection(p, 2, s) * t +
              sideCorrection(p, 3, t) * (1.0 - s) + sideCorrection(p, 1, t) * s;
    }
    return out;
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
//
// The fit itself is geom::fitCurve with its defaults: Sec. 5's "chord-length
// parameterise, least-squares fit a cubic B-spline with a fixed uniform knot
// vector", with the two ends held rather than fitted. They are held because
// they are the nodes of the arrangement: two arcs meeting at a cone have to
// arrive at the same point or the patches around it do not close, and a
// least-squares fit that is free at the ends misses each of them by its own
// residual. Holding them costs two degrees of freedom and buys exact
// interpolation of every corner of the layout.
// ---------------------------------------------------------------------------
void SplineFit::fitArcs() {
    const std::vector<Arrangement::Arc> &list = arr->getArcs();
    report.arcs = static_cast<int>(list.size());
    report.controlPointsPerArc = controlPointsPerArc();
    fitted.assign(list.size(), Curve());

    geom::FitOptions fopts;
    fopts.degree = 3;
    fopts.knots = knot;
    fopts.regularisation = options.regularisation;

    double sum = 0.0;
    long long count = 0;
    for (size_t a = 0; a < list.size(); ++a) {
        if (list[a].points.size() < 2) continue;
        geom::FitResult<2> fit = geom::fitCurve(list[a].points, fopts);
        Curve &c = fitted[a];
        c.arc = static_cast<int>(a);
        c.spline = std::move(fit.curve);
        c.maxDeviation = fit.maxDeviation;
        c.rmsDeviation = fit.rmsDeviation;
        c.length = fit.length;
        c.samples = fit.samples;
        c.underdetermined = fit.underdetermined;
        ++report.curves;
        if (c.underdetermined) ++report.underdetermined;

        const Arrangement::ArcKind kind = list[a].kind;
        const bool exact =
            (kind == Arrangement::ArcKind::Boundary && !options.fitBoundaryArcs) ||
            (kind == Arrangement::ArcKind::Interface && !options.fitInterfaceArcs);
        if (exact) {
            c.exact = true;
            c.poly = geom::Polyline<2>(list[a].points);
            ++report.exactArcs;
            if (c.maxDeviation > report.maxNetDeviation) {
                report.maxNetDeviation = c.maxDeviation;
                report.worstNetArc = static_cast<int>(a);
            }
            continue;
        }

        ++report.fittedArcs;
        if (c.maxDeviation > report.maxDeviation) {
            report.maxDeviation = c.maxDeviation;
            report.worstArc = static_cast<int>(a);
        }
        sum += c.rmsDeviation * c.rmsDeviation * c.samples;
        count += c.samples;
    }
    report.rmsDeviation = count > 0 ? std::sqrt(sum / count) / modelExtent : 0.0;
    report.maxDeviation /= modelExtent;
    report.maxNetDeviation /= modelExtent;
}

// ---------------------------------------------------------------------------
// buildPatches()
//
// The Coons blend of Sec. 5, taken at the control points (geom::coonsNet). The
// four sides come out of Stage 8 in cyclic order, so they are re-oriented once
// into the (s, t) frame of the patch:
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
            if (sd.arc < 0 || fitted[sd.arc].spline.size() != n) ok = false;
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
            std::vector<Point> c = fitted[p.side[k]].spline.controlPoints();
            if (p.forward[k] == flip) std::reverse(c.begin(), c.end());
            return c;
        };
        const std::vector<Point> C0 = oriented(0, false);  // (0,0) -> (1,0)
        const std::vector<Point> D1 = oriented(1, false);  // (1,0) -> (1,1)
        const std::vector<Point> C1 = oriented(2, true);   // (0,1) -> (1,1)
        const std::vector<Point> D0 = oriented(3, true);   // (0,0) -> (0,1)

        p.surface = geom::BSplineSurface<2>(3, knot, 3, knot,
                                            geom::coonsNet(C0, C1, D0, D1, greville, greville));

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
                        case 0: row[i] = p.surface.control(i, 0); break;
                        case 1: row[i] = p.surface.control(n - 1, i); break;
                        case 2: row[i] = p.surface.control(n - 1 - i, n - 1); break;
                        default: row[i] = p.surface.control(0, n - 1 - i); break;
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
            const Point corner = (k == 0) ? p.surface.control(0, 0)
                               : (k == 1) ? p.surface.control(n - 1, 0)
                               : (k == 2) ? p.surface.control(n - 1, n - 1)
                                          : p.surface.control(0, n - 1);
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
            if (c.spline.empty()) continue;
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
        if (c.spline.empty()) continue;
        // An exact arc is written as the polyline it is, not as that polyline
        // resampled: resampling it at `samples` stations would cut its corners
        // for exactly the reason the fit was taken off it.
        if (c.exact && c.poly.size() >= 2) {
            for (const Point &q : c.poly.points()) out << "v " << q[0] << " " << q[1] << " 0\n";
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
        for (const Point &q : p.surface.controlNet()) out << "v " << q[0] << " " << q[1] << " 0\n";
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
