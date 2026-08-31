#include "MERIDIAN/QuadMesh.hxx"

#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <fstream>
#include <limits>
#include <unordered_map>

namespace {

// Union-find over the arcs. The classes are the chords; see the header.
int findRoot(std::vector<int> &parent, int x) {
    while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
    return x;
}

void unite(std::vector<int> &parent, int a, int b) {
    a = findRoot(parent, a);
    b = findRoot(parent, b);
    if (a != b) parent[b] = a;
}

// The scaled Jacobian of a planar quadrilateral: the smallest of the four
// corner cross products of the unit vectors along its two incident sides. One
// for a square, zero for a corner that has closed up, negative for a fold.
double scaledJacobian(const Point q[4]) {
    double worst = std::numeric_limits<double>::infinity();
    for (int i = 0; i < 4; ++i) {
        const Point a = q[(i + 1) % 4] - q[i];
        const Point b = q[(i + 3) % 4] - q[i];
        const double na = normP(a), nb = normP(b);
        if (!(na > 0.0) || !(nb > 0.0)) return 0.0;
        worst = std::min(worst, cross2(a, b) / (na * nb));
    }
    return worst;
}

double quadArea(const Point q[4]) {
    return 0.5 * (cross2(q[1] - q[0], q[2] - q[0]) + cross2(q[2] - q[0], q[3] - q[0]));
}

} // namespace

// ---------------------------------------------------------------------------
QuadMesh::QuadMesh(const SplineFit &f) : QuadMesh(f, Options()) {}

QuadMesh::QuadMesh(const SplineFit &f, const Options &opts)
    : fit(&f), arr(&f.getArrangement()), options(opts) {
    const Mesh &m = arr->getMesh();
    Point lo{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    Point hi{-lo[0], -lo[1]};
    for (const Point &p : m.vertices) {
        lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
        hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
    }
    modelExtent = normP(hi - lo);
    if (!(modelExtent > 0.0)) modelExtent = 1.0;
    report.modelExtent = modelExtent;
    report.target = options.targetEdgeLength;

    if (!(options.targetEdgeLength > 0.0)) {
        report.messages.push_back("target edge length must be positive; nothing was meshed");
        return;
    }
    if (options.minIntervals < 1) options.minIntervals = 1;

    facePatch.assign(arr->getFaces().size(), -1);
    for (size_t k = 0; k < fit->patches().size(); ++k) {
        const int f2 = fit->patches()[k].face;
        if (f2 >= 0 && f2 < static_cast<int>(facePatch.size())) {
            facePatch[f2] = static_cast<int>(k);
        }
    }

    assignIntervals();
    meshArcs();
    meshPatches();
    smooth();
    classifyMaterials();
    check();
}

// ---------------------------------------------------------------------------
// tabulate() / paramAtLength() / evaluateArc()
//
// Arc length along one arc, and its inverse. Uniform in the curve parameter,
// which is not uniform in length -- that is the whole reason the table exists.
// The nodes are then placed by inverting it, so they are equally spaced along
// the arc rather than along a parameterisation that only approximates it.
// ---------------------------------------------------------------------------
Point QuadMesh::evaluateArc(int arc, double u) const {
    u = std::max(0.0, std::min(1.0, u));
    const bool onFeature =
        arc >= 0 && arc < static_cast<int>(arr->getArcs().size()) &&
        arr->getArcs()[arc].kind != Arrangement::ArcKind::Separatrix;
    if (options.useSplines && !(onFeature && options.featuresOnTracedArcs) &&
        arc < static_cast<int>(fit->curves().size()) &&
        !fit->curves()[arc].ctrl.empty()) {
        return fit->evaluate(fit->curves()[arc], u);
    }
    // The traced polyline, parameterised by its own normalised chord length.
    const std::vector<Point> &poly = arr->getArcs()[arc].points;
    if (poly.empty()) return Point{0.0, 0.0};
    if (poly.size() == 1) return poly.front();
    std::vector<double> cum(poly.size(), 0.0);
    for (size_t i = 1; i < poly.size(); ++i) cum[i] = cum[i - 1] + normP(poly[i] - poly[i - 1]);
    const double total = cum.back();
    if (!(total > 0.0)) return poly.front();
    const double s = u * total;
    const size_t k = static_cast<size_t>(
        std::lower_bound(cum.begin(), cum.end(), s) - cum.begin());
    if (k == 0) return poly.front();
    if (k >= poly.size()) return poly.back();
    const double seg = cum[k] - cum[k - 1];
    const double w = seg > 0.0 ? (s - cum[k - 1]) / seg : 0.0;
    return poly[k - 1] + (poly[k] - poly[k - 1]) * w;
}

QuadMesh::Table QuadMesh::tabulate(int arc) const {
    Table t;
    const int m = std::max(8, options.arcLengthSamples);
    t.cum.assign(m + 1, 0.0);
    Point prev = evaluateArc(arc, 0.0);
    for (int i = 1; i <= m; ++i) {
        const Point p = evaluateArc(arc, static_cast<double>(i) / m);
        t.cum[i] = t.cum[i - 1] + normP(p - prev);
        prev = p;
    }
    t.length = t.cum.back();
    return t;
}

double QuadMesh::paramAtLength(const Table &t, double s) const {
    const int m = static_cast<int>(t.cum.size()) - 1;
    if (m <= 0 || !(t.length > 0.0)) return 0.0;
    if (s <= 0.0) return 0.0;
    if (s >= t.length) return 1.0;
    const int k = static_cast<int>(
        std::lower_bound(t.cum.begin(), t.cum.end(), s) - t.cum.begin());
    if (k <= 0) return 0.0;
    if (k > m) return 1.0;
    const double seg = t.cum[k] - t.cum[k - 1];
    const double w = seg > 0.0 ? (s - t.cum[k - 1]) / seg : 0.0;
    return (static_cast<double>(k - 1) + w) / m;
}

// ---------------------------------------------------------------------------
// patchSpans()
//
// How far a row of elements has to reach across a patch, in each direction:
// the mean over the transverse parameter of the length of an isoparametric
// line. The two extreme isolines are the patch's own boundary arcs, so on a
// patch that is nearly a rectangle this returns their length and changes
// nothing; on one that tapers it returns what the interior actually spans,
// which is the number the interval assignment needs and the arc lengths do not
// carry. See the header.
// ---------------------------------------------------------------------------
void QuadMesh::patchSpans(int patch, const std::vector<Arrangement::Side> &sides,
                          double &sMean, double &tMean) const {
    sMean = tMean = 0.0;
    const int m = std::max(4, options.spanSamples);
    std::vector<Point> grid(static_cast<size_t>(m + 1) * (m + 1));

    if (options.useSplines && patch >= 0 && patch < static_cast<int>(fit->patches().size())) {
        const SplineFit::Patch &p = fit->patches()[patch];
        for (int j = 0; j <= m; ++j) {
            for (int i = 0; i <= m; ++i) {
                grid[static_cast<size_t>(j) * (m + 1) + i] =
                    fit->evaluate(p, static_cast<double>(i) / m, static_cast<double>(j) / m);
            }
        }
    } else {
        // No surface to sample, so the four sides are blended instead. The
        // stations are equal arc length along each, which is where the nodes
        // will go, so the isoline lengths this measures are the ones the
        // elements will have to span.
        if (sides.size() != 4) return;
        auto side = [&](int k, double c) {
            const bool fwd = sides[k].forward;
            const double frac = (k < 2) ? (fwd ? c : 1.0 - c) : (fwd ? 1.0 - c : c);
            const Table &t = tables[sides[k].arc];
            return evaluateArc(sides[k].arc, paramAtLength(t, t.length * frac));
        };
        std::vector<Point> B(m + 1), T(m + 1), L(m + 1), R(m + 1);
        for (int i = 0; i <= m; ++i) {
            const double c = static_cast<double>(i) / m;
            B[i] = side(0, c); T[i] = side(2, c); L[i] = side(3, c); R[i] = side(1, c);
        }
        for (int j = 0; j <= m; ++j) {
            const double v = static_cast<double>(j) / m;
            for (int i = 0; i <= m; ++i) {
                const double u = static_cast<double>(i) / m;
                grid[static_cast<size_t>(j) * (m + 1) + i] =
                    B[i] * (1.0 - v) + T[i] * v + L[j] * (1.0 - u) + R[j] * u -
                    (B.front() * ((1.0 - u) * (1.0 - v)) + B.back() * (u * (1.0 - v)) +
                     T.front() * ((1.0 - u) * v) + T.back() * (u * v));
            }
        }
    }
    for (int j = 0; j <= m; ++j) {
        double L = 0.0;
        for (int i = 0; i < m; ++i) {
            L += normP(grid[static_cast<size_t>(j) * (m + 1) + i + 1] -
                       grid[static_cast<size_t>(j) * (m + 1) + i]);
        }
        sMean += L;
    }
    for (int i = 0; i <= m; ++i) {
        double L = 0.0;
        for (int j = 0; j < m; ++j) {
            L += normP(grid[static_cast<size_t>(j + 1) * (m + 1) + i] -
                       grid[static_cast<size_t>(j) * (m + 1) + i]);
        }
        tMean += L;
    }
    sMean /= (m + 1);
    tMean /= (m + 1);
}

// ---------------------------------------------------------------------------
// assignIntervals()
//
// The integer problem, and the only part of this stage that is one. Union the
// opposite sides of every meshable patch to get the chords, then give each
// chord the count that minimises the sum of squared log ratios of its arcs'
// edge lengths to the target. See the header for why that objective and not
// another.
// ---------------------------------------------------------------------------
void QuadMesh::assignIntervals() {
    const int nArcs = static_cast<int>(arr->getArcs().size());
    intervals.assign(nArcs, -1);
    chordOf.assign(nArcs, -1);
    if (nArcs == 0) return;

    std::vector<int> parent(nArcs);
    for (int i = 0; i < nArcs; ++i) parent[i] = i;
    std::vector<char> used(nArcs, 0);

    std::vector<std::pair<int, std::vector<Arrangement::Side>>> meshable;
    for (int f : arr->patchFaces()) {
        const std::vector<Arrangement::Side> sides = arr->patchSides(f);
        if (sides.size() != 4) continue;
        bool ok = true;
        for (const Arrangement::Side &sd : sides) if (sd.arc < 0) ok = false;
        if (!ok) continue;
        unite(parent, sides[0].arc, sides[2].arc);
        unite(parent, sides[1].arc, sides[3].arc);
        for (const Arrangement::Side &sd : sides) used[sd.arc] = 1;
        meshable.push_back({f, sides});
    }

    tables.assign(nArcs, Table());
    for (int a = 0; a < nArcs; ++a) if (used[a]) tables[a] = tabulate(a);

    // One (arc, span) pair per direction of every meshable patch: the arc says
    // which class the span belongs to, the span says what that class has to
    // cover. See the header for why this and not the arc's own length.
    std::vector<std::pair<int, double>> spans;
    for (const auto &mf : meshable) {
        double sSpan = 0.0, tSpan = 0.0;
        patchSpans(facePatch[mf.first], mf.second, sSpan, tSpan);
        if (sSpan > 0.0) spans.push_back({mf.second[0].arc, sSpan});
        if (tSpan > 0.0) spans.push_back({mf.second[1].arc, tSpan});
    }

    // Group the used arcs by their root, in first-seen order so that the chord
    // list is stable from run to run.
    std::unordered_map<int, int> slotOf;
    for (int a = 0; a < nArcs; ++a) {
        if (!used[a]) continue;
        const int r = findRoot(parent, a);
        auto it = slotOf.find(r);
        if (it == slotOf.end()) {
            chordOf[a] = static_cast<int>(chordList.size());
            slotOf.emplace(r, chordOf[a]);
            chordList.push_back(Chord());
        } else {
            chordOf[a] = it->second;
        }
        chordList[chordOf[a]].arcs.push_back(a);
    }

    // The spans, gathered by the class of the arc that named them.
    std::vector<std::vector<double>> chordSpans(chordList.size());
    for (const std::pair<int, double> &sp : spans) {
        const int slot = chordOf[sp.first];
        if (slot >= 0) chordSpans[slot].push_back(sp.second);
    }

    const double h = options.targetEdgeLength;
    double sumN = 0.0;
    for (size_t slot = 0; slot < chordList.size(); ++slot) {
        Chord &c = chordList[slot];
        c.minLength = std::numeric_limits<double>::infinity();
        for (int a : c.arcs) {
            const double L = tables[a].length;
            c.minLength = std::min(c.minLength, L);
            c.maxLength = std::max(c.maxLength, L);
        }
        if (!std::isfinite(c.minLength)) c.minLength = 0.0;

        // A chord with no span at all -- every patch it touches was unfitted --
        // falls back to its own arcs, which is what it would have used before.
        std::vector<double> &data = chordSpans[slot];
        if (data.empty()) {
            for (int a : c.arcs) if (tables[a].length > 0.0) data.push_back(tables[a].length);
        }

        double logSum = 0.0;
        int counted = 0;
        c.minSpan = std::numeric_limits<double>::infinity();
        for (double S : data) {
            c.minSpan = std::min(c.minSpan, S);
            c.maxSpan = std::max(c.maxSpan, S);
            if (S > 0.0) { logSum += std::log(S); ++counted; }
        }
        if (!std::isfinite(c.minSpan)) c.minSpan = 0.0;

        if (counted == 0) {
            // Every arc of the chord came out with no length at all. There is
            // nothing to divide, so it takes the floor and is reported.
            c.idealIntervals = 0.0;
            c.intervals = options.minIntervals;
            c.clamped = true;
        } else {
            const double geo = std::exp(logSum / counted);
            c.idealIntervals = geo / h;
            // The minimiser of F over the integers is one of the two integers
            // bracketing N*, and not necessarily the nearer of them.
            const int lo = std::max(1, static_cast<int>(std::floor(c.idealIntervals)));
            const int hi2 = std::max(1, static_cast<int>(std::ceil(c.idealIntervals)));
            auto cost = [&](int N) {
                double s = 0.0;
                for (double S : data) {
                    if (!(S > 0.0)) continue;
                    const double r = std::log(S / (N * h));
                    s += r * r;
                }
                return s;
            };
            c.intervals = (cost(lo) <= cost(hi2)) ? lo : hi2;

            const int floorN = std::max(1, options.minIntervals);
            if (c.intervals < floorN) { c.intervals = floorN; c.clamped = true; }
            if (options.maxIntervals > 0 && c.intervals > options.maxIntervals) {
                c.intervals = options.maxIntervals;
                c.clamped = true;
            }
        }

        c.minEdge = c.minSpan / c.intervals;
        c.maxEdge = c.maxSpan / c.intervals;
        for (int a : c.arcs) intervals[a] = c.intervals;

        if (c.clamped) ++report.clampedChords;
        sumN += c.intervals;
    }

    report.chords = static_cast<int>(chordList.size());
    report.minIntervals = std::numeric_limits<int>::max();
    for (const Chord &c : chordList) {
        report.minIntervals = std::min(report.minIntervals, c.intervals);
        report.maxIntervals = std::max(report.maxIntervals, c.intervals);
        report.arcsAssigned += static_cast<int>(c.arcs.size());
    }
    if (chordList.empty()) report.minIntervals = 0;
    report.meanIntervals = chordList.empty() ? 0.0 : sumN / chordList.size();
}

// ---------------------------------------------------------------------------
// meshArcs()
//
// Every arc that bounds a meshable patch is cut once, here, and the vertices go
// into a list both its patches read. That is the whole of the conformity
// argument: two patches sharing an arc do not each place points on it and hope
// they agree, they use the same indices.
// ---------------------------------------------------------------------------
void QuadMesh::meshArcs() {
    const std::vector<Arrangement::Arc> &list = arr->getArcs();
    arcNodes.assign(list.size(), {});
    arcParams.assign(list.size(), {});
    nodeVert.assign(arr->getNodes().size(), -1);

    auto vertexForNode = [&](int n) {
        if (n < 0 || n >= static_cast<int>(nodeVert.size())) return -1;
        if (nodeVert[n] < 0) {
            nodeVert[n] = static_cast<int>(verts.size());
            verts.push_back(arr->getNodes()[n].p);
        }
        return nodeVert[n];
    };

    for (size_t a = 0; a < list.size(); ++a) {
        const int N = intervals[a];
        if (N <= 0) continue;
        const Arrangement::Arc &ar = list[a];
        arcNodes[a].assign(N + 1, -1);
        arcParams[a].assign(N + 1, 0.0);

        arcNodes[a][0] = vertexForNode(ar.from);
        arcNodes[a][N] = vertexForNode(ar.to);
        arcParams[a][0] = 0.0;
        arcParams[a][N] = 1.0;

        const Table &t = tables[a];
        for (int k = 1; k < N; ++k) {
            const double u = (t.length > 0.0)
                                 ? paramAtLength(t, t.length * k / N)
                                 : static_cast<double>(k) / N;
            arcParams[a][k] = u;
            arcNodes[a][k] = static_cast<int>(verts.size());
            verts.push_back(evaluateArc(static_cast<int>(a), u));
        }
    }
}

// ---------------------------------------------------------------------------
// meshPatches()
//
// The (s, t) frame is Stage 9's: corner 0 at (0,0), side 0 the bottom running
// s from 0 to 1, then round counter-clockwise. Sides 0 and 1 are traversed in
// the direction the frame wants when the arc runs forwards along the face; 2
// and 3 are traversed against it, so their parameters come back reversed. That
// single asymmetry is where every index below comes from.
//
// The interior is not a bilinear blend of the boundary *points*. The four sides
// give a parameter distribution each, those are blended, and the patch is then
// evaluated at the blended parameter -- so an interior node lies on the bicubic
// surface exactly, and the grid follows the arc-length spacing of the sides
// into the middle of the patch instead of drifting off it.
// ---------------------------------------------------------------------------
void QuadMesh::meshPatches() {
    const std::vector<Arrangement::Face> &faces = arr->getFaces();
    const double domain = arr->getReport().domainArea;

    for (int f : arr->patchFaces()) {
        const std::vector<Arrangement::Side> sides = arr->patchSides(f);
        const int pi = (f < static_cast<int>(facePatch.size())) ? facePatch[f] : -1;
        bool ok = sides.size() == 4 && pi >= 0;
        if (ok) {
            for (const Arrangement::Side &sd : sides) {
                if (sd.arc < 0 || intervals[sd.arc] <= 0) ok = false;
            }
        }
        if (ok && (intervals[sides[0].arc] != intervals[sides[2].arc] ||
                   intervals[sides[1].arc] != intervals[sides[3].arc])) {
            // The union-find made this impossible; if it happens the layout is
            // not what patchSides() reported and the patch is left out rather
            // than meshed into a mesh with a hanging node.
            ok = false;
            report.messages.push_back("opposite sides of a patch disagree on their interval "
                                      "count; the patch was left unmeshed");
        }
        if (!ok) {
            ++report.unmeshedPatches;
            if (domain > 0.0) report.unmeshedArea += faces[f].area / domain;
            continue;
        }

        const int ns = intervals[sides[0].arc];
        const int nt = intervals[sides[1].arc];

        Block b;
        b.face = f;
        b.patch = pi;
        b.ns = ns;
        b.nt = nt;
        b.vert.assign(static_cast<size_t>(ns + 1) * (nt + 1), -1);
        auto at = [&](int i, int j) -> int& {
            return b.vert[static_cast<size_t>(j) * (ns + 1) + i];
        };

        // The patch coordinate a point at curve parameter w on side k sits at.
        auto coord = [&](int k, double w) {
            const bool fwd = sides[k].forward;
            if (k < 2) return fwd ? w : 1.0 - w;
            return fwd ? 1.0 - w : w;
        };

        std::vector<double> sB(ns + 1), sT(ns + 1), tL(nt + 1), tR(nt + 1);
        for (int i = 0; i <= ns; ++i) {
            const int mb = sides[0].forward ? i : ns - i;
            const int mt = sides[2].forward ? ns - i : i;
            at(i, 0) = arcNodes[sides[0].arc][mb];
            at(i, nt) = arcNodes[sides[2].arc][mt];
            sB[i] = coord(0, arcParams[sides[0].arc][mb]);
            sT[i] = coord(2, arcParams[sides[2].arc][mt]);
        }
        for (int j = 0; j <= nt; ++j) {
            const int ml = sides[3].forward ? nt - j : j;
            const int mr = sides[1].forward ? j : nt - j;
            at(0, j) = arcNodes[sides[3].arc][ml];
            at(ns, j) = arcNodes[sides[1].arc][mr];
            tL[j] = coord(3, arcParams[sides[3].arc][ml]);
            tR[j] = coord(1, arcParams[sides[1].arc][mr]);
        }

        const SplineFit::Patch &sp = fit->patches()[pi];
        for (int j = 1; j < nt; ++j) {
            const double v = static_cast<double>(j) / nt;
            for (int i = 1; i < ns; ++i) {
                const double u = static_cast<double>(i) / ns;
                const double s = (1.0 - v) * sB[i] + v * sT[i];
                const double t = (1.0 - u) * tL[j] + u * tR[j];
                Point p;
                if (options.useSplines) {
                    p = fit->evaluate(sp, s, t);
                } else {
                    // Discrete Coons on the boundary points, for the same
                    // reason the arcs came off the polylines: no spline is
                    // being trusted anywhere in this mode.
                    const Point &b0 = verts[at(i, 0)], &b1 = verts[at(i, nt)];
                    const Point &l0 = verts[at(0, j)], &l1 = verts[at(ns, j)];
                    const Point &c00 = verts[at(0, 0)], &c10 = verts[at(ns, 0)];
                    const Point &c01 = verts[at(0, nt)], &c11 = verts[at(ns, nt)];
                    p = b0 * (1.0 - v) + b1 * v + l0 * (1.0 - u) + l1 * u -
                        (c00 * ((1.0 - u) * (1.0 - v)) + c10 * (u * (1.0 - v)) +
                         c01 * ((1.0 - u) * v) + c11 * (u * v));
                }
                at(i, j) = static_cast<int>(verts.size());
                verts.push_back(p);
            }
        }

        const int blockId = static_cast<int>(grids.size());
        for (int j = 0; j < nt; ++j) {
            for (int i = 0; i < ns; ++i) {
                cells.push_back({{at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)}});
                cellBlock.push_back(blockId);
            }
        }
        report.patchArea += faces[f].area;
        grids.push_back(std::move(b));
    }

    report.blocks = static_cast<int>(grids.size());
    report.vertices = static_cast<int>(verts.size());
    report.quads = static_cast<int>(cells.size());
}

// ---------------------------------------------------------------------------
// smooth()
//
// Winslow's elliptic smoother, Gauss-Seidel, on the blocks that need it,
// boundary held. See the header for why this is here at all and not a finishing
// touch: the transfinite grid folds on the layout faces a Coons patch itself
// folds on, and this is what takes the fold out.
//
// The coefficients are recomputed at every node from the current grid, which is
// what makes the system elliptic rather than harmonic and what makes it pull a
// crowded row apart instead of merely averaging it. `alpha + gamma` vanishes
// only where the node's four neighbours have collapsed onto it, and there is
// nothing to solve there, so it is left alone.
//
// Because only interior nodes move, no vertex that appears in another block is
// touched, and the conformity established in meshPatches() survives untouched
// -- check() re-measures it afterwards rather than taking that on trust.
// ---------------------------------------------------------------------------
void QuadMesh::smooth() {
    if (options.smoothingPasses <= 0 || cells.empty()) return;

    // What the grid looked like before, so that the report can say what the
    // smoothing was worth.
    double before = std::numeric_limits<double>::infinity();
    for (const std::array<int, 4> &q : cells) {
        const Point c[4] = {verts[q[0]], verts[q[1]], verts[q[2]], verts[q[3]]};
        before = std::min(before, scaledJacobian(c));
        if (!(quadArea(c) > 0.0)) ++report.invertedBefore;
    }
    report.minScaledJacobianBefore = std::isfinite(before) ? before : 0.0;

    const double tol = options.smoothingTolerance * options.targetEdgeLength;
    for (const Block &b : grids) {
        if (b.ns < 2 || b.nt < 2) continue;
        const int w = b.ns + 1;
        auto id = [&](int i, int j) { return b.vert[static_cast<size_t>(j) * w + i]; };

        // Only blocks that came out badly; see Options::smoothingThreshold.
        double worst = std::numeric_limits<double>::infinity();
        for (int j = 0; j < b.nt; ++j) {
            for (int i = 0; i < b.ns; ++i) {
                const Point c[4] = {verts[id(i, j)], verts[id(i + 1, j)],
                                    verts[id(i + 1, j + 1)], verts[id(i, j + 1)]};
                worst = std::min(worst, scaledJacobian(c));
            }
        }
        if (!(worst < options.smoothingThreshold)) continue;
        ++report.smoothedBlocks;

        int sweep = 0;
        for (; sweep < options.smoothingPasses; ++sweep) {
            double moved = 0.0;
            for (int j = 1; j < b.nt; ++j) {
                for (int i = 1; i < b.ns; ++i) {
                    const Point &xe = verts[id(i + 1, j)];
                    const Point &xw = verts[id(i - 1, j)];
                    const Point &xn = verts[id(i, j + 1)];
                    const Point &xs = verts[id(i, j - 1)];
                    const Point ds = (xe - xw) * 0.5;
                    const Point dt = (xn - xs) * 0.5;
                    const double alpha = dotP(dt, dt);
                    const double beta = dotP(ds, dt);
                    const double gamma = dotP(ds, ds);
                    const double denom = 2.0 * (alpha + gamma);
                    if (!(denom > 0.0)) continue;
                    const Point cross = (verts[id(i + 1, j + 1)] - verts[id(i - 1, j + 1)] -
                                         verts[id(i + 1, j - 1)] + verts[id(i - 1, j - 1)]) *
                                        0.25;
                    const Point next = ((xe + xw) * alpha + (xn + xs) * gamma -
                                        cross * (2.0 * beta)) / denom;
                    Point &here = verts[id(i, j)];
                    moved = std::max(moved, normP(next - here));
                    here = next;
                }
            }
            if (moved <= tol) { ++sweep; break; }
        }
        report.smoothingSweeps = std::max(report.smoothingSweeps, sweep);
    }
}

// ---------------------------------------------------------------------------
// check()
//
// What the construction promises, read back off the result: every edge used
// twice or once and never three times, no two vertices at the same point, no
// element folded, and the edge lengths near the target that the chords were
// chosen for. All of it is measured rather than assumed, which is what makes it
// worth having -- an error in the index arithmetic above shows up here as a
// crack or a third use of an edge and nowhere else.
// ---------------------------------------------------------------------------
// classifyMaterials()
//
// The material of every element, and whether any of them straddles an
// interface.
//
// This is the property the whole multi-material path exists to produce, and it
// is measured on the finished elements rather than inferred from the layout,
// because it is the thing that is actually true or false about the output: an
// element straddling an interface carries two materials, and no analysis code
// can integrate one. Five samples per element -- the centroid and the four edge
// midpoints -- located in the input triangulation, which is where the material
// tags live. An element inside one region has all five agree; one lying across
// an interface does not, whatever the layout says about itself.
//
// The edge midpoints are the samples that matter. A centroid alone would call
// an element pure whenever the interface clipped only a corner of it, which is
// exactly the case a nearly-aligned layout produces and exactly the one worth
// catching.
// ---------------------------------------------------------------------------
void QuadMesh::classifyMaterials() {
    const Mesh &m = arr->getMesh();
    cellMaterial.assign(cells.size(), 0);
    report.materials = 0;
    report.mixedQuads = report.unlocatedQuads = report.interfaceEdges = 0;
    if (m.triangleMatId.size() != m.triangles.size()) return;

    std::set<int> seen;
    for (size_t k = 0; k < cells.size(); ++k) {
        const std::array<int, 4> &q = cells[k];
        Point c{0.0, 0.0};
        for (int i = 0; i < 4; ++i) c = c + verts[q[i]];
        c = c / 4.0;

        // Five points of the element's *interior*: the centroid and the four
        // half-way points from it to the corners.
        //
        // The obvious sampling -- the edge midpoints, nudged inwards -- is
        // wrong, and wrong in the direction that reports a defect where there
        // is none. An element edge lying on a curved interface is a chord of
        // it, so its midpoint is off the interface by the sagitta whichever
        // side the element is on; on geom011, whose interface is a cosine with
        // a radius of curvature of about three element lengths, that put 37 of
        // 984 elements on the wrong side of their own edge. Nothing is wrong
        // with those elements: a straight-sided quadrilateral cannot follow a
        // curve exactly and is not meant to. What would be wrong is an element
        // the interface runs *through*, and the interior is where to look for
        // that.
        std::vector<Point> samples{c};
        for (int i = 0; i < 4; ++i) samples.push_back((c + verts[q[i]]) * 0.5);

        int mat = 0;
        bool mixed = false;
        for (const Point &p : samples) {
            const int t = m.findTriangleContainingPoint(p);
            if (t < 0) continue;
            const int id = m.triangleMatId[t];
            if (mat == 0) mat = id;
            else if (id != mat) mixed = true;
        }
        if (mat == 0) {
            // Every sample fell outside the triangulation, which happens where
            // the Coons patch bulges a fraction of an element past a curved
            // piece of dS. The element is still in the model and still in one
            // material; the nearest triangle says which.
            ++report.unlocatedQuads;
            double best = std::numeric_limits<double>::infinity();
            for (size_t t = 0; t < m.triangles.size(); ++t) {
                const Triangle &tri = m.triangles[t];
                const Point g = (m.vertices[tri[0]] + m.vertices[tri[1]] + m.vertices[tri[2]]) / 3.0;
                const double d = normP(g - c);
                if (d < best) { best = d; mat = m.triangleMatId[t]; }
            }
            if (mat == 0) continue;
        }
        cellMaterial[k] = mat;
        seen.insert(mat);
        if (mixed) ++report.mixedQuads;
    }
    report.materials = static_cast<int>(seen.size());

    // Element edges two materials share. Both sides carry the same nodes by
    // construction -- the interface is one arc of the layout, meshed once --
    // so this is a count of what the two regions have in common and not a
    // check that they meet.
    std::map<std::pair<int, int>, std::pair<int, int>> edgeMat;
    for (size_t k = 0; k < cells.size(); ++k) {
        if (cellMaterial[k] == 0) continue;
        for (int i = 0; i < 4; ++i) {
            int a = cells[k][i], b = cells[k][(i + 1) % 4];
            if (a > b) std::swap(a, b);
            auto &slot = edgeMat[{a, b}];
            if (slot.first == 0) slot.first = cellMaterial[k];
            else if (slot.second == 0) slot.second = cellMaterial[k];
        }
    }
    for (const auto &kv : edgeMat) {
        if (kv.second.first != 0 && kv.second.second != 0 &&
            kv.second.first != kv.second.second) {
            ++report.interfaceEdges;
        }
    }
    report.materialsPure = report.mixedQuads == 0;
}

// ---------------------------------------------------------------------------
void QuadMesh::check() {
    report.vertices = static_cast<int>(verts.size());
    report.quads = static_cast<int>(cells.size());
    if (cells.empty()) {
        report.messages.push_back("no patch of the layout could be meshed");
        return;
    }

    const long long n = static_cast<long long>(verts.size());
    std::unordered_map<long long, int> edgeUse;
    edgeUse.reserve(cells.size() * 4);
    double lenSum = 0.0, logSum = 0.0;
    int lenCount = 0;
    report.minEdge = std::numeric_limits<double>::infinity();
    const double h = options.targetEdgeLength;

    for (const std::array<int, 4> &q : cells) {
        for (int k = 0; k < 4; ++k) {
            const int a = q[k], b = q[(k + 1) % 4];
            const long long key = static_cast<long long>(std::min(a, b)) * n + std::max(a, b);
            const int before = edgeUse[key]++;
            if (before == 0) {
                const double L = normP(verts[b] - verts[a]);
                report.minEdge = std::min(report.minEdge, L);
                report.maxEdge = std::max(report.maxEdge, L);
                lenSum += L;
                ++lenCount;
                if (L > 0.0 && h > 0.0) {
                    const double r = std::log(L / h);
                    logSum += r * r;
                    if (std::fabs(r) > std::fabs(std::log(report.worstEdgeRatio))) {
                        report.worstEdgeRatio = L / h;
                    }
                }
            }
        }
    }
    if (!std::isfinite(report.minEdge)) report.minEdge = 0.0;
    report.meanEdge = lenCount > 0 ? lenSum / lenCount : 0.0;
    report.edgeRatioRms = lenCount > 0 ? std::sqrt(logSum / lenCount) : 0.0;

    for (const auto &e : edgeUse) {
        if (e.second == 1) ++report.boundaryEdges;
        else if (e.second == 2) ++report.interiorEdges;
        else ++report.nonManifoldEdges;
    }

    double sjSum = 0.0;
    report.minScaledJacobian = std::numeric_limits<double>::infinity();
    report.minQuadArea = std::numeric_limits<double>::infinity();
    for (const std::array<int, 4> &q : cells) {
        const Point c[4] = {verts[q[0]], verts[q[1]], verts[q[2]], verts[q[3]]};
        const double area = quadArea(c);
        const double sj = scaledJacobian(c);
        report.meshArea += area;
        report.minQuadArea = std::min(report.minQuadArea, area);
        report.maxQuadArea = std::max(report.maxQuadArea, area);
        report.minScaledJacobian = std::min(report.minScaledJacobian, sj);
        sjSum += sj;
        if (!(area > 0.0)) ++report.invertedQuads;
    }
    if (!std::isfinite(report.minScaledJacobian)) report.minScaledJacobian = 0.0;
    if (!std::isfinite(report.minQuadArea)) report.minQuadArea = 0.0;
    report.meanScaledJacobian = sjSum / cells.size();

    // Cracks: two distinct vertices at the same point of the model. A hash grid
    // at four times the tolerance, so a pair either shares a cell or is in one
    // of the eight around it.
    const double tol = options.crackTolerance * modelExtent;
    const double cell = std::max(tol * 4.0, modelExtent * 1e-12);
    std::unordered_map<long long, std::vector<int>> grid;
    auto key = [&](double x, double y) {
        return static_cast<long long>(std::floor(x / cell)) * 73856093LL ^
               static_cast<long long>(std::floor(y / cell)) * 19349663LL;
    };
    for (size_t i = 0; i < verts.size(); ++i) grid[key(verts[i][0], verts[i][1])].push_back(static_cast<int>(i));
    for (size_t i = 0; i < verts.size(); ++i) {
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                auto it = grid.find(key(verts[i][0] + dx * cell, verts[i][1] + dy * cell));
                if (it == grid.end()) continue;
                for (int j : it->second) {
                    if (j <= static_cast<int>(i)) continue;
                    if (normP(verts[j] - verts[i]) <= tol) ++report.cracks;
                }
            }
        }
    }

    // Corners of a block that the model turns through more than a half turn
    // at. The element there is reversed whatever the grid does inside, so this
    // is the part of `invertedQuads` that is the layout's and not the mesh's.
    for (const Block &b : grids) {
        if (b.ns < 1 || b.nt < 1) continue;
        const int w = b.ns + 1;
        auto id = [&](int i, int j) { return b.vert[static_cast<size_t>(j) * w + i]; };
        const int ci[4] = {0, b.ns, b.ns, 0};
        const int cj[4] = {0, 0, b.nt, b.nt};
        const int ai[4] = {1, b.ns, b.ns - 1, 0};
        const int aj[4] = {0, 1, b.nt, b.nt - 1};
        const int bi[4] = {0, b.ns - 1, b.ns, 1};
        const int bj[4] = {1, 0, b.nt - 1, b.nt};
        for (int k = 0; k < 4; ++k) {
            const Point o = verts[id(ci[k], cj[k])];
            const Point u = verts[id(ai[k], aj[k])] - o;
            const Point v = verts[id(bi[k], bj[k])] - o;
            // The interior of the block is on the left of o -> u, so the
            // interior angle is the turn from u round to v that way.
            double a = std::atan2(cross2(u, v), dotP(u, v));
            if (a < 0.0) a += 2.0 * M_PI;
            if (a > M_PI) ++report.reflexCorners;
        }
    }

    report.conforming = report.nonManifoldEdges == 0 && report.cracks == 0;
    report.valid = report.conforming && report.invertedQuads == 0 &&
                   report.unmeshedPatches == 0 && report.mixedQuads == 0;

    if (report.unmeshedPatches > 0) {
        report.messages.push_back(
            std::to_string(report.unmeshedPatches) +
            " face(s) of the layout are not quadrilaterals with one arc a side and were left "
            "unmeshed; the remedy is Sec. 3.3's repair, upstream of this stage");
    }
    if (report.invertedQuads > 0) {
        std::string why = "; smoothing cannot reach them";
        if (report.reflexCorners == 0) why = "; the patch they are in folds under its blend";
        report.messages.push_back(
            std::to_string(report.invertedQuads) + " element(s) came out with non-positive area" +
            (report.reflexCorners > 0
                 ? ", against " + std::to_string(report.reflexCorners) +
                       " corner(s) of the layout that the model turns through more than a "
                       "half turn at" + why
                 : why));
    }
    if (report.nonManifoldEdges > 0) {
        report.messages.push_back(std::to_string(report.nonManifoldEdges) +
                                  " edge(s) are used by more than two elements");
    }
    if (report.cracks > 0) {
        report.messages.push_back(std::to_string(report.cracks) +
                                  " pair(s) of distinct vertices are at the same point");
    }
}

// ---------------------------------------------------------------------------
bool QuadMesh::writeOBJ(const std::string &filename) const {
    std::ofstream out(filename);
    if (!out) return false;
    out << "# MERIDIAN Stage 10: " << report.quads << " quadrilaterals, target edge "
        << options.targetEdgeLength << "\n";
    for (const Point &p : verts) out << "v " << p[0] << " " << p[1] << " 0\n";
    for (size_t k = 0; k < cells.size(); ++k) {
        const std::array<int, 4> &q = cells[k];
        out << "f " << (q[0] + 1) << " " << (q[1] + 1) << " " << (q[2] + 1) << " "
            << (q[3] + 1) << "\n";
    }
    return true;
}

bool QuadMesh::writeVTU(const std::string &filename) const {
    std::ofstream out(filename);
    if (!out) return false;
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << verts.size() << "\" NumberOfCells=\""
        << cells.size() << "\">\n";
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const Point &p : verts) out << "          " << p[0] << " " << p[1] << " 0.0\n";
    out << "        </DataArray>\n      </Points>\n";
    out << "      <Cells>\n";
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const std::array<int, 4> &q : cells) {
        out << "          " << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << "\n";
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (size_t k = 0; k < cells.size(); ++k) out << "          " << 4 * (k + 1) << "\n";
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t k = 0; k < cells.size(); ++k) out << "          9\n";
    out << "        </DataArray>\n      </Cells>\n";
    out << "      <CellData Scalars=\"scaledJacobian\">\n";
    out << "        <DataArray type=\"Float64\" Name=\"scaledJacobian\" format=\"ascii\">\n";
    for (const std::array<int, 4> &q : cells) {
        const Point c[4] = {verts[q[0]], verts[q[1]], verts[q[2]], verts[q[3]]};
        out << "          " << scaledJacobian(c) << "\n";
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"block\" format=\"ascii\">\n";
    for (size_t k = 0; k < cells.size(); ++k) {
        out << "          " << (k < cellBlock.size() ? cellBlock[k] : -1) << "\n";
    }
    out << "        <DataArray type=\"Int32\" Name=\"material\" format=\"ascii\">\n";
    for (size_t k = 0; k < cells.size(); ++k) {
        out << "          " << (k < cellMaterial.size() ? cellMaterial[k] : 0) << "\n";
    }
    out << "        </DataArray>\n";
    out << "        </DataArray>\n      </CellData>\n";
    out << "    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";
    return true;
}
