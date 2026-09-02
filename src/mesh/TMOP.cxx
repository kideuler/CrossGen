// TMOP.cxx -- see TMOP.hxx. A node-local Newton relaxation of the TMOP energy
// over mesh::QuadMesh, swept in colour order so it runs in parallel without
// changing its answer.
#include "TMOP.hxx"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <sstream>

#ifdef _OPENMP
#  include <omp.h>
#  define CG_PRAGMA(x) _Pragma(#x)
#  define CG_OMP(x) CG_PRAGMA(omp x)
#else
#  define CG_OMP(x)
#endif

namespace mesh {

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

// The quadrature rules of TMOP::Quadrature, on the reference square [0,1]^2.
// Both have unit total weight, so summing w * det(W) over an element gives the
// element's target area.
struct QuadPoint { double xi, eta, w; };

const QuadPoint kGauss2x2[4] = {
    // 0.5 +- 1/(2 sqrt 3): the two-point Gauss rule mapped from [-1,1] to [0,1].
    {0.5 - 0.28867513459481288, 0.5 - 0.28867513459481288, 0.25},
    {0.5 + 0.28867513459481288, 0.5 - 0.28867513459481288, 0.25},
    {0.5 + 0.28867513459481288, 0.5 + 0.28867513459481288, 0.25},
    {0.5 - 0.28867513459481288, 0.5 + 0.28867513459481288, 0.25},
};

const QuadPoint kCorners[4] = {
    {0.0, 0.0, 0.25}, {1.0, 0.0, 0.25}, {1.0, 1.0, 0.25}, {0.0, 1.0, 0.25},
};

// The cofactor matrix of a 2x2, which is d(det T)/dT.
inline Jacobian2 cofactor(const Jacobian2 &T) {
    return {T[3], -T[2], -T[1], T[0]};
}

}  // namespace

// ---------------------------------------------------------------------------
// the metrics
// ---------------------------------------------------------------------------

// Every metric here is a function of the two invariants s = |T|^2 and
// t = det T alone, which is what makes the exact local Hessian below short:
// the second derivative of t with respect to T contracts to zero against the
// rank-one perturbations a single node can make, so only the five partials of
// f survive. See the derivation above `accumulate`.
TMOP::MetricPartials TMOP::partialsOf(Metric m, double s, double t,
                                      double gamma, double tau0) {
    MetricPartials p;
    switch (m) {
        case Shape004: {
            p.valid = true;
            p.f = s - 2.0 * t;
            p.fs = 1.0;
            p.ft = -2.0;
            return p;
        }
        case Shape002: {
            if (!(t > 0.0)) return p;
            p.valid = true;
            p.f  = 0.5 * s / t - 1.0;
            p.fs = 0.5 / t;
            p.ft = -0.5 * s / (t * t);
            p.fst = -0.5 / (t * t);
            p.ftt = s / (t * t * t);
            return p;
        }
        case ShapeSize007: {
            if (!(t > 0.0)) return p;
            const double t2 = t * t;
            p.valid = true;
            p.f  = s * (1.0 + 1.0 / t2) - 4.0;
            p.fs = 1.0 + 1.0 / t2;
            p.ft = -2.0 * s / (t2 * t);
            p.fst = -2.0 / (t2 * t);
            p.ftt = 6.0 * s / (t2 * t2);
            return p;
        }
        case Untangle022: {
            const double d = t - tau0;
            if (!(d > 0.0)) return p;
            const double n = s - 2.0 * t;
            p.valid = true;
            p.f  = 0.5 * n / d;
            p.fs = 0.5 / d;
            p.ft = -1.0 / d - 0.5 * n / (d * d);
            p.fst = -0.5 / (d * d);
            p.ftt = 2.0 / (d * d) + n / (d * d * d);
            return p;
        }
        case Size055: {
            p.valid = true;
            p.f = (t - 1.0) * (t - 1.0);
            p.ft = 2.0 * (t - 1.0);
            p.ftt = 2.0;
            return p;
        }
        case Size056: {
            if (!(t > 0.0)) return p;
            p.valid = true;
            p.f  = 0.5 * (t + 1.0 / t) - 1.0;
            p.ft = 0.5 * (1.0 - 1.0 / (t * t));
            p.ftt = 1.0 / (t * t * t);
            return p;
        }
        case ShapeSizeCombo: {
            const double g = std::min(std::max(gamma, 0.0), 1.0);
            const MetricPartials a = partialsOf(Shape002, s, t, 0.0, 0.0);
            const MetricPartials b = partialsOf(Size056, s, t, 0.0, 0.0);
            if (!a.valid || !b.valid) return p;
            p.valid = true;
            p.f   = (1.0 - g) * a.f   + g * b.f;
            p.fs  = (1.0 - g) * a.fs  + g * b.fs;
            p.ft  = (1.0 - g) * a.ft  + g * b.ft;
            p.fss = (1.0 - g) * a.fss + g * b.fss;
            p.fst = (1.0 - g) * a.fst + g * b.fst;
            p.ftt = (1.0 - g) * a.ftt + g * b.ftt;
            return p;
        }
    }
    return p;
}

// The energy is the integral of mu^p, so every partial the local system needs
// picks up the chain rule for that power. Written once here rather than per
// metric: p only ever sees mu and its derivatives, never which metric produced
// them.
//
// At p >= 2 this stays finite where mu vanishes -- f^(p-2) is 1 at p = 2 and 0
// above it -- which matters, because mu = 0 is precisely where a converged
// element sits. Between 1 and 2 it is singular there, which is why
// Options::exponent refuses that range.
TMOP::MetricPartials TMOP::raiseTo(const MetricPartials &b, double p) {
    if (p == 1.0 || !b.valid) return b;
    // Every metric here is non-negative wherever it is defined, so this guard
    // only ever catches rounding at a converged element, where mu is zero and
    // so is every derivative of mu^p for p > 2.
    if (!(b.f > 0.0)) { MetricPartials z; z.valid = true; return z; }
    const double f = b.f;
    const double fm1 = std::pow(f, p - 1.0);
    const double fm2 = std::pow(f, p - 2.0);
    const double c = p * (p - 1.0) * fm2;
    MetricPartials r;
    r.valid = true;
    r.f   = f * fm1;
    r.fs  = p * fm1 * b.fs;
    r.ft  = p * fm1 * b.ft;
    r.fss = c * b.fs * b.fs + p * fm1 * b.fss;
    r.fst = c * b.fs * b.ft + p * fm1 * b.fst;
    r.ftt = c * b.ft * b.ft + p * fm1 * b.ftt;
    return r;
}

bool TMOP::evalMetric(Metric m, const Jacobian2 &T, double gamma, double tau0,
                      double *mu, Jacobian2 *dmu, double exponent) {
    const double s = T[0] * T[0] + T[1] * T[1] + T[2] * T[2] + T[3] * T[3];
    const double t = det2(T);
    const MetricPartials p = raiseTo(partialsOf(m, s, t, gamma, tau0), exponent);
    if (!p.valid) {
        if (mu) *mu = kInf;
        if (dmu) *dmu = Jacobian2{0.0, 0.0, 0.0, 0.0};
        return false;
    }
    if (mu) *mu = p.f;
    if (dmu) {
        // dmu/dT = f_s ds/dT + f_t dt/dT = 2 f_s T + f_t cof(T).
        const Jacobian2 C = cofactor(T);
        for (int k = 0; k < 4; ++k) (*dmu)[k] = 2.0 * p.fs * T[k] + p.ft * C[k];
    }
    return true;
}

// ---------------------------------------------------------------------------
// construction and targets
// ---------------------------------------------------------------------------

TMOP::TMOP(QuadMesh &m) : mesh(m) {}
TMOP::TMOP(QuadMesh &m, const Options &opts) : mesh(m), options(opts) {}

void TMOP::prepare() {
    const int nQ = static_cast<int>(mesh.quads.size());

    mesh.computeQuality();
    meanEdge = mesh.quality.meanEdge > 0.0 ? mesh.quality.meanEdge : 1.0;

    switch (options.target) {
        case TargetKeep:
            if (mesh.targetJacobian.size() != mesh.quads.size())
                mesh.setUniformTargets(options.targetSize);
            break;
        case TargetUniformSquare:
            mesh.setUniformTargets(options.targetSize);
            break;
        case TargetCurrentShape:
            mesh.setTargetsFromCurrentShape();
            break;
        case TargetPerQuadSize:
            mesh.setSizeTargets(options.perQuadSize);
            break;
    }

    targetInverse.assign(nQ, Jacobian2{1.0, 0.0, 0.0, 1.0});
    targetDet.assign(nQ, 1.0);
    int singular = 0;
    for (int q = 0; q < nQ; ++q) {
        const Jacobian2 W = mesh.targetAt(q, 0);
        bool ok = false;
        const Jacobian2 Wi = inv2(W, ok);
        // A target with no inverse is a target that says nothing; fall back to
        // an ideal square of the mesh's mean edge length rather than dropping
        // the element out of the energy altogether.
        if (!ok || !(det2(W) > 0.0)) {
            ++singular;
            const Jacobian2 F{meanEdge, 0.0, 0.0, meanEdge};
            bool ok2 = false;
            targetInverse[q] = inv2(F, ok2);
            targetDet[q] = det2(F);
        } else {
            targetInverse[q] = Wi;
            targetDet[q] = det2(W);
        }
    }
    if (singular > 0) {
        std::ostringstream os;
        os << singular << " element(s) had a singular or inverted target Jacobian; "
           << "an ideal square of the mean edge length was used for them";
        report.messages.push_back(os.str());
    }

    if (options.exponent > 1.0 && options.exponent < 2.0) {
        std::ostringstream os;
        os << "Options::exponent was " << options.exponent
           << "; raised to 2 (the second derivative of mu^p is singular at mu = 0 "
           << "for any p strictly between 1 and 2)";
        report.messages.push_back(os.str());
        options.exponent = 2.0;
    } else if (options.exponent < 1.0) {
        report.messages.push_back("Options::exponent below 1 was raised to 1");
        options.exponent = 1.0;
    }

    activeMetric = options.metric;
    activeTau0 = 0.0;
    activeExponent = options.exponent;
}

// ---------------------------------------------------------------------------
// energy, gradient and the local Hessian
// ---------------------------------------------------------------------------

// The whole of the derivative bookkeeping, in one place.
//
// With W fixed, T = A W^{-1} and A = sum_c x_c (x) grad N_c, so moving node v
// (which sits at corner `c` of this element) changes T by a rank-one matrix:
//
//     dT/dx_v[i] = M_i,   M_i has row i equal to b and its other row zero,
//     b_k = sum_j W^{-1}[j][k] grad N_c[j]      (i.e. b = W^{-T} grad N_c)
//
// Writing mu as f(s, t) with s = |T|^2 and t = det T, that gives
//
//     dE/dx_i   = 2 f_s (T b)_i + f_t (C b)_i,          C = cof(T)
//     d2E/dx_i dx_l = 2 f_s |b|^2 delta_il
//                   + 4 f_ss (Tb)_i (Tb)_l
//                   + 2 f_st [(Tb)_i (Cb)_l + (Cb)_i (Tb)_l]
//                   +   f_tt (Cb)_i (Cb)_l
//
// The second derivative of t with respect to T contributes nothing: its only
// non-zero entries pair an index in row 0 of T with one in row 1, and M_i has
// just one non-zero row, so every such pair either vanishes or cancels against
// its partner. That is why the local Hessian is exact and still this short --
// no finite differencing, and no metric-specific second-derivative code beyond
// the five partials of f.
bool TMOP::accumulate(int q, int c, double *energyOut, double *minTau,
                      Point *grad, double *hess) const {
    const QuadPoint *rule = (options.quadrature == Corners) ? kCorners : kGauss2x2;
    const Jacobian2 &Wi = targetInverse[q];
    const double dW = targetDet[q];

    bool valid = true;
    for (int k = 0; k < 4; ++k) {
        const QuadPoint &qp = rule[k];
        const Jacobian2 A = mesh.jacobianAt(q, qp.xi, qp.eta);
        const Jacobian2 T = mul2(A, Wi);

        const double s = T[0] * T[0] + T[1] * T[1] + T[2] * T[2] + T[3] * T[3];
        const double t = det2(T);
        if (minTau) *minTau = std::min(*minTau, t);

        const MetricPartials p = raiseTo(
            partialsOf(activeMetric, s, t, options.gamma, activeTau0), activeExponent);
        if (!p.valid) {
            valid = false;
            if (energyOut) *energyOut = kInf;
            if (!minTau) return false;   // nothing left to learn from this element
            continue;                     // keep scanning: minTau still wants the rest
        }
        if (energyOut && valid) *energyOut += qp.w * dW * p.f;
        if (!grad && !hess) continue;

        const Point g = QuadMesh::shapeGrad(c, qp.xi, qp.eta);
        // b = W^{-T} grad N_c, with Wi row-major: Wi[j][k] = Wi[2j + k].
        const double b0 = Wi[0] * g[0] + Wi[2] * g[1];
        const double b1 = Wi[1] * g[0] + Wi[3] * g[1];

        const double Tb0 = T[0] * b0 + T[1] * b1;
        const double Tb1 = T[2] * b0 + T[3] * b1;
        const Jacobian2 C = cofactor(T);
        const double Cb0 = C[0] * b0 + C[1] * b1;
        const double Cb1 = C[2] * b0 + C[3] * b1;

        const double scale = qp.w * dW;
        if (grad) {
            (*grad)[0] += scale * (2.0 * p.fs * Tb0 + p.ft * Cb0);
            (*grad)[1] += scale * (2.0 * p.fs * Tb1 + p.ft * Cb1);
        }
        if (hess) {
            const double bb = b0 * b0 + b1 * b1;
            const double diag = 2.0 * p.fs * bb;
            hess[0] += scale * (diag + 4.0 * p.fss * Tb0 * Tb0
                                + 4.0 * p.fst * Tb0 * Cb0 + p.ftt * Cb0 * Cb0);
            hess[3] += scale * (diag + 4.0 * p.fss * Tb1 * Tb1
                                + 4.0 * p.fst * Tb1 * Cb1 + p.ftt * Cb1 * Cb1);
            const double off = scale * (4.0 * p.fss * Tb0 * Tb1
                                        + 2.0 * p.fst * (Tb0 * Cb1 + Cb0 * Tb1)
                                        + p.ftt * Cb0 * Cb1);
            hess[1] += off;
            hess[2] += off;
        }
    }
    return valid;
}

double TMOP::elementEnergy(int q) const {
    if (q < 0 || q >= static_cast<int>(mesh.quads.size())) return 0.0;
    double e = 0.0;
    if (!accumulate(q, 0, &e, nullptr, nullptr, nullptr)) return kInf;
    return e;
}

// Evaluated in parallel, summed in element order.
//
// A `reduction(+:total)` would be shorter and would give a slightly different
// last bit at each thread count, because floating-point addition is not
// associative. That does not matter for the number itself -- but this number
// decides when run() stops, and a mesh that took one more sweep on eight
// threads than on one would make the colouring's whole guarantee untrue where
// anyone would notice it. Per-element results into a buffer, then one ordered
// pass, costs an array and buys exact reproducibility.
double TMOP::energy() const {
    const int nQ = static_cast<int>(mesh.quads.size());
    std::vector<double> per(nQ, 0.0);
    CG_OMP(parallel for schedule(static))
    for (int q = 0; q < nQ; ++q) per[q] = elementEnergy(q);
    double total = 0.0;
    for (int q = 0; q < nQ; ++q) total += per[q];
    return total;
}

double TMOP::nodeEnergy(int v, double *minTau) const {
    if (v < 0 || v + 1 >= static_cast<int>(mesh.vertexQuads.rowPtr.size())) return 0.0;
    double e = 0.0;
    bool valid = true;
    for (int i = mesh.vertexQuads.begin(v); i < mesh.vertexQuads.end(v); ++i) {
        const int q = mesh.vertexQuads.colIdx[i];
        double eq = 0.0;
        if (!accumulate(q, mesh.vertexQuads.corner[i], &eq, minTau, nullptr, nullptr))
            valid = false;
        else if (valid) e += eq;
    }
    return valid ? e : kInf;
}

void TMOP::nodeSystem(int v, Point &gradient, double hessian[4]) const {
    gradient = Point{0.0, 0.0};
    for (int k = 0; k < 4; ++k) hessian[k] = 0.0;
    if (v < 0 || v + 1 >= static_cast<int>(mesh.vertexQuads.rowPtr.size())) return;
    for (int i = mesh.vertexQuads.begin(v); i < mesh.vertexQuads.end(v); ++i)
        accumulate(mesh.vertexQuads.colIdx[i], mesh.vertexQuads.corner[i],
                   nullptr, nullptr, &gradient, hessian);
}

Point TMOP::nodeGradient(int v) const {
    Point g{0.0, 0.0};
    if (v < 0 || v + 1 >= static_cast<int>(mesh.vertexQuads.rowPtr.size())) return g;
    for (int i = mesh.vertexQuads.begin(v); i < mesh.vertexQuads.end(v); ++i)
        accumulate(mesh.vertexQuads.colIdx[i], mesh.vertexQuads.corner[i],
                   nullptr, nullptr, &g, nullptr);
    return g;
}

// ---------------------------------------------------------------------------
// the local solve
// ---------------------------------------------------------------------------

Point TMOP::nodeStep(int v) {
    const Point zero{0.0, 0.0};
    if (v < 0 || v >= static_cast<int>(mesh.vertices.size())) return zero;
    if (v >= static_cast<int>(mesh.nodeType.size())) return zero;
    if (mesh.nodeType[v] == QuadMesh::NodeFixed) return zero;

    // Where the line search starts from, and whether this patch is somewhere
    // the metric is defined at all. Taken before the local system rather than
    // after it because an infinite starting energy makes the system meaningless
    // -- a barrier metric on a patch that still holds an inverted element.
    double tauBefore = kInf;
    const double e0 = nodeEnergy(v, &tauBefore);
    if (!std::isfinite(e0)) return zero;

    // The tangent has to be the one that goes with the *current* positions of
    // this node's feature neighbours: they may have moved earlier in this same
    // sweep, and a stale tangent would walk the node off its own polyline.
    mesh.updateSlideTangent(v);

    Point g{0.0, 0.0};
    double H[4] = {0.0, 0.0, 0.0, 0.0};
    nodeSystem(v, g, H);
    if (!std::isfinite(g[0]) || !std::isfinite(g[1])) return zero;
    if (normP(g) <= 0.0) return zero;

    // Symmetrise; the two off-diagonals are equal in exact arithmetic.
    const double h01 = 0.5 * (H[1] + H[2]);
    double a = H[0], b = h01, d = H[3];

    Point dir{0.0, 0.0};
    const bool sliding = (mesh.nodeType[v] == QuadMesh::NodeSliding);
    if (sliding) {
        const Point &t = mesh.slideTangent[v];
        if (normP(t) <= 0.0) return zero;
        // The Newton system restricted to the one direction the node may take.
        const double gt = dotP(g, t);
        double ht = a * t[0] * t[0] + 2.0 * b * t[0] * t[1] + d * t[1] * t[1];
        const double floorH = options.hessianFloor * std::max(std::fabs(a) + std::fabs(d), 1e-300);
        if (!(ht > floorH)) ht = std::max(floorH, 1e-300);
        dir = t * (-gt / ht);
    } else {
        // Lift the smaller eigenvalue of the symmetric 2x2 to a positive floor,
        // so an indefinite local Hessian -- which a barrier metric does produce
        // near a nearly inverted element -- still gives a descent direction.
        const double tr = a + d;
        const double disc = std::sqrt(std::max((a - d) * (a - d) + 4.0 * b * b, 0.0));
        const double lmax = 0.5 * (tr + disc);
        const double lmin = 0.5 * (tr - disc);
        const double floorH = options.hessianFloor * std::max(std::fabs(lmax), 1e-300);
        if (lmin < floorH) {
            const double lift = floorH - lmin;
            a += lift;
            d += lift;
        }
        const double det = a * d - b * b;
        if (!(std::fabs(det) > 0.0)) {
            // Nothing usable: fall back to steepest descent, scaled below.
            dir = g * -1.0;
        } else {
            dir = Point{-(d * g[0] - b * g[1]) / det, -(a * g[1] - b * g[0]) / det};
        }
    }

    dir = mesh.projectStep(v, dir);
    const double len = normP(dir);
    if (!(len > 0.0) || !std::isfinite(len)) return zero;

    // Cap at a fraction of the shortest edge at this node, so one Newton step
    // can never carry a node across its own one-ring however wild the local
    // curvature is.
    double shortest = kInf;
    for (int n : mesh.vertexNeighbors[v])
        shortest = std::min(shortest, normP(mesh.vertices[n] - mesh.vertices[v]));
    if (!std::isfinite(shortest) || shortest <= 0.0) shortest = meanEdge;
    const double cap = options.maxStepFraction * shortest;
    if (len > cap) dir = dir * (cap / len);

    // Backtracking. The move is accepted only when this node's own patch energy
    // strictly falls; and if the patch was valid before, only when it still is
    // -- a non-barrier metric would otherwise happily fold an element to lower
    // its own number.
    const Point origin = mesh.vertices[v];
    Point accepted = zero;
    double alpha = 1.0;
    for (int k = 0; k < options.maxLineSearch; ++k) {
        const Point trial = dir * alpha;
        mesh.vertices[v] = origin + trial;
        double tauAfter = kInf;
        const double e = nodeEnergy(v, &tauAfter);
        const bool stillValid = !(tauBefore > 0.0) || tauAfter > 0.0;
        if (std::isfinite(e) && e < e0 && stillValid) { accepted = trial; break; }
        alpha *= 0.5;
    }
    mesh.vertices[v] = origin;
    return accepted;
}

double TMOP::moveNode(int v) {
    const Point step = nodeStep(v);
    const double len = normP(step);
    if (len > 0.0) mesh.vertices[v] = mesh.vertices[v] + step;
    return len;
}

// ---------------------------------------------------------------------------
// colouring and sweeps
// ---------------------------------------------------------------------------

// Two nodes may be moved at the same time exactly when no element contains
// both: the local solve reads and writes only the elements around its node, so
// disjoint element stencils means no interference. Colouring the "shares an
// element" graph and running one colour at a time therefore parallelises a
// Gauss-Seidel sweep *without approximating it* -- the answer is the same as
// the serial sweep that visits the nodes in colour order, and the same whatever
// the thread count.
void TMOP::buildColoring() {
    const int nV = static_cast<int>(mesh.vertices.size());
    color.assign(nV, -1);
    buckets.clear();

    std::vector<int> stamp;   // colour -> the vertex that last forbade it
    int maxColor = -1;

    for (int v = 0; v < nV; ++v) {
        if (v >= static_cast<int>(mesh.nodeType.size())) break;
        if (mesh.nodeType[v] == QuadMesh::NodeFixed) continue;

        for (int i = mesh.vertexQuads.begin(v); i < mesh.vertexQuads.end(v); ++i) {
            const Quad &q = mesh.quads[mesh.vertexQuads.colIdx[i]];
            for (int c = 0; c < 4; ++c) {
                const int w = q[c];
                if (w == v || w < 0 || w >= nV) continue;
                const int cw = color[w];
                if (cw < 0) continue;
                if (cw >= static_cast<int>(stamp.size())) stamp.resize(cw + 1, -1);
                stamp[cw] = v;
            }
        }

        int pick = 0;
        while (pick < static_cast<int>(stamp.size()) && stamp[pick] == v) ++pick;
        color[v] = pick;
        maxColor = std::max(maxColor, pick);
        if (pick >= static_cast<int>(stamp.size())) stamp.resize(pick + 1, -1);
    }

    buckets.assign(maxColor + 1, {});
    for (int v = 0; v < nV; ++v)
        if (color[v] >= 0) buckets[color[v]].push_back(v);
    report.colors = static_cast<int>(buckets.size());
}

double TMOP::sweep() {
    double maxMove = 0.0;
    for (std::size_t c = 0; c < buckets.size(); ++c) {
        const std::vector<int> &bucket = buckets[c];
        const int n = static_cast<int>(bucket.size());
        double bucketMax = 0.0;
        CG_OMP(parallel for schedule(static) reduction(max:bucketMax))
        for (int i = 0; i < n; ++i) {
            const double d = moveNode(bucket[i]);
            if (d > bucketMax) bucketMax = d;
        }
        maxMove = std::max(maxMove, bucketMax);
    }
    return maxMove;
}

// ---------------------------------------------------------------------------
// the run
// ---------------------------------------------------------------------------

double TMOP::minDeterminant() const {
    const int nQ = static_cast<int>(mesh.quads.size());
    double worst = kInf;
    CG_OMP(parallel for schedule(static) reduction(min:worst))
    for (int q = 0; q < nQ; ++q) {
        const Jacobian2 &Wi = targetInverse[q];
        // Corners as well as the active rule: a barrier placed below the
        // quadrature points alone can still be above a folded corner.
        for (int k = 0; k < 4; ++k) {
            const QuadPoint &a = kCorners[k];
            const QuadPoint &b = (options.quadrature == Corners) ? kCorners[k] : kGauss2x2[k];
            worst = std::min(worst, det2(mul2(mesh.jacobianAt(q, a.xi, a.eta), Wi)));
            worst = std::min(worst, det2(mul2(mesh.jacobianAt(q, b.xi, b.eta), Wi)));
        }
    }
    return worst;
}

int TMOP::countInverted() const {
    const int nQ = static_cast<int>(mesh.quads.size());
    int bad = 0;
    CG_OMP(parallel for schedule(static) reduction(+:bad))
    for (int q = 0; q < nQ; ++q)
        if (mesh.minScaledJacobian(q) <= 0.0) ++bad;
    return bad;
}

void TMOP::snapshotQuality(bool before) {
    mesh.computeQuality();
    const QuadMesh::Quality &Q = mesh.quality;
    if (before) {
        report.minScaledJacobianBefore = Q.minScaledJacobian;
        report.meanScaledJacobianBefore = Q.meanScaledJacobian;
        report.invertedBefore = Q.invertedQuads;
        report.worstAspectBefore = Q.worstAspect;
        report.minAreaBefore = Q.minArea;
        report.areaBefore = Q.totalArea;
        report.freeNodes = Q.freeNodes;
        report.slidingNodes = Q.slidingNodes;
        report.fixedNodes = Q.fixedNodes;
        report.movableNodes = Q.freeNodes + Q.slidingNodes;
    } else {
        report.minScaledJacobianAfter = Q.minScaledJacobian;
        report.meanScaledJacobianAfter = Q.meanScaledJacobian;
        report.invertedAfter = Q.invertedQuads;
        report.worstAspectAfter = Q.worstAspect;
        report.minAreaAfter = Q.minArea;
        report.areaAfter = Q.totalArea;
    }
}

bool TMOP::run() {
    const std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    report = Report();

#ifdef _OPENMP
    report.openMP = true;
    const int savedThreads = omp_get_max_threads();
    if (options.threads > 0) omp_set_num_threads(options.threads);
    report.threads = omp_get_max_threads();
#else
    report.threads = 1;
    if (options.threads > 1)
        report.messages.push_back("built without OpenMP; Options::threads was ignored");
#endif

    prepare();
    buildColoring();
    snapshotQuality(true);

    if (mesh.quads.empty() || report.movableNodes == 0) {
        report.messages.push_back(mesh.quads.empty()
                                      ? "no elements to smooth"
                                      : "every node is fixed; nothing to smooth");
#ifdef _OPENMP
        omp_set_num_threads(savedThreads);
#endif
        return false;
    }

    startPositions = mesh.vertices;
    double targetArea = 0.0;
    for (double d : targetDet) targetArea += d;
    if (!(targetArea > 0.0)) targetArea = 1.0;

    // -- phase 1: untangling ------------------------------------------------
    // A barrier metric is infinite on an inverted element, so the main phase
    // cannot even be started on a tangled mesh. Metric 22 puts its barrier at
    // tau_0 below the worst determinant instead of at zero, which makes the
    // energy finite everywhere and its gradient push the determinants up; tau_0
    // is re-tracked each sweep, so the barrier follows the mesh as it recovers.
    if (options.untangle && report.invertedBefore > 0) {
        activeMetric = Untangle022;
        activeExponent = 1.0;
        for (int s = 0; s < options.untangleMaxSweeps; ++s) {
            const double tmin = minDeterminant();
            if (tmin > 0.0) break;
            activeTau0 = options.untangleFloor * std::min(tmin, 0.0)
                         - 1e-9 * std::max(1.0, std::fabs(tmin));
            const double move = sweep();
            ++report.untangleSweeps;
            if (move < options.moveTolerance * meanEdge) break;
        }
        const int left = countInverted();
        std::ostringstream os;
        os << "untangling: " << report.invertedBefore << " inverted element(s) at the start, "
           << left << " after " << report.untangleSweeps << " sweep(s)";
        report.messages.push_back(os.str());
    }

    // -- phase 2: the metric the caller asked for ---------------------------
    activeMetric = options.metric;
    activeTau0 = 0.0;
    activeExponent = options.exponent;
    if (activeMetric == Untangle022) {
        // Asked for directly rather than as a phase: place the barrier once,
        // below whatever the mesh currently is, and leave it there.
        const double tmin = minDeterminant();
        activeTau0 = std::min(tmin, 0.0) * options.untangleFloor
                     - 1e-9 * std::max(1.0, std::fabs(tmin));
    }

    double previous = energy();
    const double energyScale = std::max(std::fabs(previous), 1e-300);
    if (!std::isfinite(previous)) {
        report.messages.push_back(
            "the mesh is still tangled and the chosen metric is a barrier metric; "
            "no smoothing was attempted (try Options::untangle, or Metric::Shape004)");
        snapshotQuality(false);
        // Both infinite, and reported as such rather than as zero: the metric
        // genuinely has no value on this mesh, and a zero here would read as a
        // perfect one.
        report.energyBefore = report.energyAfter = previous / targetArea;
        report.ran = false;
        report.seconds = std::chrono::duration<double>(
                             std::chrono::steady_clock::now() - t0).count();
#ifdef _OPENMP
        omp_set_num_threads(savedThreads);
#endif
        return false;
    }
    report.energyBefore = previous / targetArea;

    for (int s = 0; s < options.maxSweeps; ++s) {
        const double move = sweep();
        ++report.sweeps;
        report.lastSweepMove = move / meanEdge;
        const double now = energy();
        if (move < options.moveTolerance * meanEdge) { report.converged = true; break; }
        // Measured against the energy the mesh *started* with, not against what
        // is left of it. A Gauss-Seidel sweep contracts the error by a roughly
        // constant factor, so the relative drop per sweep never gets small and
        // a test against the current energy would only fire at roundoff; what
        // says the run is finished is that a sweep no longer buys a meaningful
        // fraction of what there was to gain.
        if (std::fabs(previous - now) <= options.energyTolerance * energyScale) {
            report.converged = true;
            previous = now;
            break;
        }
        previous = now;
    }
    report.energyAfter = energy() / targetArea;

    // The tangents and the quality are the two position-dependent things the
    // mesh carries; both are stale now and both are cheap.
    mesh.computeSlideTangents();
    snapshotQuality(false);

    double total = 0.0;
    for (std::size_t v = 0; v < mesh.vertices.size(); ++v) {
        const double d = normP(mesh.vertices[v] - startPositions[v]);
        total += d;
        report.maxDisplacement = std::max(report.maxDisplacement, d);
    }
    report.maxDisplacement /= meanEdge;
    report.meanDisplacement = report.movableNodes > 0
                                  ? total / (report.movableNodes * meanEdge)
                                  : 0.0;

    report.ran = true;
    report.seconds = std::chrono::duration<double>(
                         std::chrono::steady_clock::now() - t0).count();

#ifdef _OPENMP
    omp_set_num_threads(savedThreads);
#endif

    // "Worse than it started" is judged on the two numbers a mesh is accepted
    // or rejected on, not on the energy: the energy is what was minimised, so
    // it going down proves nothing about the mesh being usable.
    const bool worse = report.invertedAfter > report.invertedBefore ||
                       (report.invertedAfter == 0 &&
                        report.minScaledJacobianAfter < report.minScaledJacobianBefore - 1e-12);
    if (worse) {
        report.messages.push_back(
            "the mesh came out no better than it went in; it has been left as the "
            "optimizer finished, not reverted");
        return false;
    }
    return true;
}

}  // namespace mesh
