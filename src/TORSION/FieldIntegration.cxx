#include "FieldIntegration.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

namespace {

inline double signedArea(const Point &a, const Point &b, const Point &c) {
    return 0.5 * cross2(b - a, c - a);
}

// R_k as a 2x2 acting on a tangent, row-major.
inline void quarterTurn(int k, double R[4]) {
    switch (((k % 4) + 4) % 4) {
        case 0: R[0] = 1;  R[1] = 0;  R[2] = 0;  R[3] = 1;  break;
        case 1: R[0] = 0;  R[1] = -1; R[2] = 1;  R[3] = 0;  break;
        case 2: R[0] = -1; R[1] = 0;  R[2] = 0;  R[3] = -1; break;
        default: R[0] = 0; R[1] = 1;  R[2] = -1; R[3] = 0;  break;
    }
}

inline int dofU(int v) { return 2 * v; }
inline int dofV(int v) { return 2 * v + 1; }

} // namespace

FieldIntegration::FieldIntegration(const ConeCut &cut, const FieldFrames &frames,
                                   const Immersion &scaffold, const Options &options)
    : opts(options) {
    assemble(cut, frames, scaffold);
    if (report.solved) measure(cut, frames, scaffold);
}

// ---------------------------------------------------------------------------
// assemble()
//
// L is the cotangent stiffness matrix twice over -- once for u, once for v --
// because sum_t A_t grad(phi_i) . grad(phi_j) *is* the cotan Laplacian, and
// nothing here needs it written any other way. b is div X and div Y on the same
// weights. The seam rows go in as they stand and the pin as two more.
// ---------------------------------------------------------------------------
void FieldIntegration::assemble(const ConeCut &cut, const FieldFrames &frames,
                                const Immersion &scaffold) {
    const Mesh &om = cut.getCutMesh();
    const int nV = static_cast<int>(om.vertices.size());
    const int nT = static_cast<int>(om.triangles.size());
    const int n = 2 * nV;
    report.vertices = nV;
    report.faces = nT;

    uv.assign(nV, Point{0.0, 0.0});
    if (nV == 0 || nT == 0) {
        report.messages.push_back("Omega is empty; there is nothing to integrate.");
        return;
    }
    if (frames.frames().size() != static_cast<size_t>(nT)) {
        report.messages.push_back("The frames have one J* per face of S and Omega has a "
                                  "different number of faces; the two do not belong to the "
                                  "same mesh.");
        return;
    }

    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(static_cast<size_t>(nT) * 18 + static_cast<size_t>(nV) * 2);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);

    double diagSum = 0.0;
    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = om.triangles[t];
        const Point &p0 = om.vertices[tri[0]];
        const Point &p1 = om.vertices[tri[1]];
        const Point &p2 = om.vertices[tri[2]];
        const double twoA = cross2(p1 - p0, p2 - p0);
        if (!(std::fabs(twoA) > 0.0)) continue;
        const double A = 0.5 * twoA;

        // grad phi_i = R90(p_{i+2} - p_{i+1}) / 2A, R90(x, y) = (-y, x).
        Point g[3];
        for (int i = 0; i < 3; ++i) {
            const Point e = om.vertices[tri[(i + 2) % 3]] - om.vertices[tri[(i + 1) % 3]];
            g[i] = Point{-e[1] / twoA, e[0] / twoA};
        }

        // The target gradients: the rows of J*_t.
        const std::array<double, 4> &Js = frames.frames()[t];
        const Point X{Js[0], Js[1]};
        const Point Y{Js[2], Js[3]};

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                const double w = A * dotP(g[i], g[j]);
                trip.emplace_back(dofU(tri[i]), dofU(tri[j]), w);
                trip.emplace_back(dofV(tri[i]), dofV(tri[j]), w);
                if (i == j) diagSum += 2.0 * std::fabs(w);
            }
            rhs[dofU(tri[i])] += A * dotP(X, g[i]);
            rhs[dofV(tri[i])] += A * dotP(Y, g[i]);
        }
    }
    const double scaleL = (n > 0) ? std::max(diagSum / n, 1e-30) : 1.0;

    // --- the constraints --------------------------------------------------
    const auto &pairs = scaffold.getSeamPairs();
    const auto &pairArc = scaffold.getSeamPairArc();
    const auto &arcs = scaffold.getArcs();
    report.seamPairs = static_cast<int>(pairs.size());

    std::vector<Eigen::Triplet<double>> ctrip;
    ctrip.reserve(pairs.size() * 16 + 2);
    int row = 0;
    for (size_t p = 0; p < pairs.size(); ++p) {
        const int a = pairArc[p];
        if (a < 0 || a >= static_cast<int>(arcs.size())) continue;
        double R[4];
        quarterTurn(arcs[a].k, R);
        const int ip = pairs[p][0], jp = pairs[p][1];
        const int im = pairs[p][2], jm = pairs[p][3];

        // Row 0:  (u_jp - u_ip) - R00 (u_jm - u_im) - R01 (v_jm - v_im) = 0
        // Row 1:  (v_jp - v_ip) - R10 (u_jm - u_im) - R11 (v_jm - v_im) = 0
        for (int c = 0; c < 2; ++c) {
            const int du = (c == 0) ? dofU(jp) : dofV(jp);
            const int di = (c == 0) ? dofU(ip) : dofV(ip);
            ctrip.emplace_back(row + c, du, 1.0);
            ctrip.emplace_back(row + c, di, -1.0);
            ctrip.emplace_back(row + c, dofU(jm), -R[2 * c + 0]);
            ctrip.emplace_back(row + c, dofU(im),  R[2 * c + 0]);
            ctrip.emplace_back(row + c, dofV(jm), -R[2 * c + 1]);
            ctrip.emplace_back(row + c, dofV(im),  R[2 * c + 1]);
        }
        row += 2;
    }

    // The pin. The energy's null space is a constant added to u and a constant
    // added to v; the tangent constraints are blind to both, so it survives
    // them and has to be removed here.
    int pin = opts.pinnedVertex;
    if (pin < 0 || pin >= nV) pin = 0;
    ctrip.emplace_back(row + 0, dofU(pin), 1.0);
    ctrip.emplace_back(row + 1, dofV(pin), 1.0);
    const int pinRow = row;
    row += 2;
    const int m = row;
    report.constraintRows = m;

    // --- the saddle system ------------------------------------------------
    const double eps = opts.regularisation * scaleL;
    std::vector<Eigen::Triplet<double>> all = trip;
    all.reserve(trip.size() + ctrip.size() * 2 + m);
    for (const auto &c : ctrip) {
        all.emplace_back(n + c.row(), c.col(), c.value());
        all.emplace_back(c.col(), n + c.row(), c.value());
    }
    for (int i = 0; i < m; ++i) all.emplace_back(n + i, n + i, -eps);

    Eigen::SparseMatrix<double> K(n + m, n + m);
    K.setFromTriplets(all.begin(), all.end());
    K.makeCompressed();

    Eigen::VectorXd b = Eigen::VectorXd::Zero(n + m);
    b.head(n) = rhs;
    b[n + pinRow + 0] = opts.pinnedTo[0];
    b[n + pinRow + 1] = opts.pinnedTo[1];

    Eigen::VectorXd sol;
    {
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
        ldlt.compute(K);
        if (ldlt.info() == Eigen::Success) {
            sol = ldlt.solve(b);
            if (ldlt.info() == Eigen::Success && sol.allFinite()) {
                report.factorised = true;
                report.solvedWithLDLT = true;
            }
        }
    }
    if (!report.factorised) {
        // The regularisation is what makes the saddle matrix quasi definite and
        // the LDL^T stable on it; where that is not enough -- a seam graph with
        // junctions, so C with dependent rows -- an LU with pivoting still is.
        Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
        lu.compute(K);
        if (lu.info() != Eigen::Success) {
            report.messages.push_back("The constrained least-squares system would not "
                                      "factorise, by Cholesky or by LU.");
            return;
        }
        sol = lu.solve(b);
        if (lu.info() != Eigen::Success || !sol.allFinite()) {
            report.messages.push_back("The constrained least-squares system would not solve.");
            return;
        }
        report.factorised = true;
        report.messages.push_back(
            "The saddle system needed the LU fallback: the LDL^T found it indefinite past its "
            "regularisation, which is what a seam graph with junctions in it looks like.");
    }

    for (int v = 0; v < nV; ++v) uv[v] = Point{sol[dofU(v)], sol[dofV(v)]};
    report.solved = true;
}

// ---------------------------------------------------------------------------
// measure()
//
// Sec. 7.1's census, and the per-triangle non-integrability behind it.
// ---------------------------------------------------------------------------
void FieldIntegration::measure(const ConeCut &cut, const FieldFrames &frames,
                               const Immersion &scaffold) {
    const Mesh &om = cut.getCutMesh();
    const Mesh &orig = cut.getOriginalMesh();
    const int nT = static_cast<int>(om.triangles.size());

    Point lo = uv.empty() ? Point{0.0, 0.0} : uv[0];
    Point hi = lo;
    for (const Point &p : uv) {
        lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
        hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
    }
    report.imageExtent = std::hypot(hi[0] - lo[0], hi[1] - lo[1]);
    const double extent = (report.imageExtent > 0.0) ? report.imageExtent : 1.0;

    // The seam, as it actually came out.
    const auto &pairs = scaffold.getSeamPairs();
    const auto &pairArc = scaffold.getSeamPairArc();
    const auto &arcs = scaffold.getArcs();
    for (size_t p = 0; p < pairs.size(); ++p) {
        const int a = pairArc[p];
        if (a < 0 || a >= static_cast<int>(arcs.size())) continue;
        const Point dp = uv[pairs[p][1]] - uv[pairs[p][0]];
        const Point dm = uv[pairs[p][3]] - uv[pairs[p][2]];
        const Point r = dp - Immersion::rotateQuarter(dm, arcs[a].k);
        report.maxSeamResidual = std::max(report.maxSeamResidual, normP(r) / extent);
    }

    // Cones, for "how far is each flip from one".
    std::vector<char> coneVertex(orig.vertices.size(), 0);
    std::vector<Point> conePos;
    for (const auto &c : scaffold.getCones().getCones()) {
        coneVertex[c.vertex] = 1;
        conePos.push_back(orig.vertices[c.vertex]);
    }
    const auto &c2o = cut.getCutVertexToOriginal();

    double modelDiag = 0.0;
    {
        Point mlo = orig.vertices.empty() ? Point{0, 0} : orig.vertices[0], mhi = mlo;
        for (const Point &p : orig.vertices) {
            mlo[0] = std::min(mlo[0], p[0]); mlo[1] = std::min(mlo[1], p[1]);
            mhi[0] = std::max(mhi[0], p[0]); mhi[1] = std::max(mhi[1], p[1]);
        }
        modelDiag = std::hypot(mhi[0] - mlo[0], mhi[1] - mlo[1]);
    }
    if (!(modelDiag > 0.0)) modelDiag = 1.0;

    fitResidual.assign(nT, 0.0);
    areaRatio.assign(nT, 0.0);
    report.minAreaRatio = std::numeric_limits<double>::infinity();
    double nearest = std::numeric_limits<double>::infinity();
    double farthest = 0.0;
    double residualSum = 0.0;
    int residualCount = 0;

    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = om.triangles[t];
        const Point &p0 = om.vertices[tri[0]];
        const Point &p1 = om.vertices[tri[1]];
        const Point &p2 = om.vertices[tri[2]];
        const double twoA = cross2(p1 - p0, p2 - p0);
        const double A = 0.5 * twoA;
        report.totalArea += std::fabs(A);
        if (!(std::fabs(twoA) > 0.0)) continue;

        Point g[3];
        for (int i = 0; i < 3; ++i) {
            const Point e = om.vertices[tri[(i + 2) % 3]] - om.vertices[tri[(i + 1) % 3]];
            g[i] = Point{-e[1] / twoA, e[0] / twoA};
        }
        Point gu{0.0, 0.0}, gv{0.0, 0.0};
        for (int i = 0; i < 3; ++i) {
            gu = gu + g[i] * uv[tri[i]][0];
            gv = gv + g[i] * uv[tri[i]][1];
        }

        const std::array<double, 4> &Js = frames.frames()[t];
        const Point X{Js[0], Js[1]}, Y{Js[2], Js[3]};
        const double scale = dotP(X, X) + dotP(Y, Y);
        const Point ru = gu - X, rv = gv - Y;
        const double res = (scale > 0.0) ? std::sqrt((dotP(ru, ru) + dotP(rv, rv)) / scale) : 0.0;
        fitResidual[t] = res;
        residualSum += res;
        ++residualCount;
        if (res > report.maxFitResidual) { report.maxFitResidual = res; report.worstFitFace = t; }

        const double detJs = Js[0] * Js[3] - Js[1] * Js[2];
        const double image = signedArea(uv[tri[0]], uv[tri[1]], uv[tri[2]]);
        const double want = detJs * A;
        const double ratio = (std::fabs(want) > 0.0) ? image / want : 0.0;
        areaRatio[t] = ratio;
        report.minAreaRatio = std::min(report.minAreaRatio, ratio);

        if (ratio <= opts.flipTolerance) {
            ++report.flippedFaces;
            report.flippedArea += std::fabs(A);
            const Point centre = (orig.vertices[orig.triangles[t][0]] +
                                  orig.vertices[orig.triangles[t][1]] +
                                  orig.vertices[orig.triangles[t][2]]) / 3.0;
            double d = std::numeric_limits<double>::infinity();
            for (const Point &c : conePos) d = std::min(d, normP(centre - c));
            if (std::isfinite(d)) {
                nearest = std::min(nearest, d);
                farthest = std::max(farthest, d);
            }
            for (int i = 0; i < 3; ++i) {
                const int ov = c2o[tri[i]];
                if (ov >= 0 && coneVertex[ov]) { ++report.flipsAdjacentToCone; break; }
            }
        }
    }
    if (!std::isfinite(report.minAreaRatio)) report.minAreaRatio = 0.0;
    report.meanFitResidual = residualCount ? residualSum / residualCount : 0.0;
    report.nearestFlipToCone = std::isfinite(nearest) ? nearest / modelDiag : 0.0;
    report.farthestFlipFromCone = farthest / modelDiag;

    if (report.flippedFaces > 0) {
        std::ostringstream oss;
        oss << report.flippedFaces << " of " << nT << " face(s) came out inverted ("
            << (100.0 * report.flippedArea / std::max(report.totalArea, 1e-30))
            << "% of the model by area; " << report.flipsAdjacentToCone
            << " in the one ring of a cone, the rest between "
            << report.nearestFlipToCone << " and " << report.farthestFlipFromCone
            << " of the model's diagonal from the nearest one). A discrete cross field is "
            << "generically non-integrable, so this is the expected cost of the substitution "
            << "and not a defect; Sec. 7.2's untangling is what pays it.";
        report.messages.push_back(oss.str());
    }
    if (report.maxSeamResidual > 1e-9) {
        std::ostringstream oss;
        oss << "The seam came out " << report.maxSeamResidual << " of the extent from exact. "
            << "The constraints are equalities and should be met to rounding; this is the "
            << "multiplier-block regularisation showing through, and lowering it is the fix.";
        report.messages.push_back(oss.str());
    }

    report.valid = report.solved && report.flippedFaces == 0 && report.maxSeamResidual < 1e-9;
}
