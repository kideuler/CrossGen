#include "LayoutEnergy.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>

namespace {

const double kInf = std::numeric_limits<double>::infinity();

// The closed-form SVD of a 2x2 matrix, stored row-major.
//
//     A = U diag(s0, s1) V^T,     U = rot(phi),  V = rot(psi)
//
// Writing E = (a+d)/2, F = (a-d)/2, G = (c+b)/2, H = (c-b)/2 turns the four
// entries into two rotating pairs: (E, H) turns with phi - psi and has modulus
// (s0 + s1)/2, while (F, G) turns with phi + psi and has modulus (s0 - s1)/2.
// So the two moduli give the singular values and the two arguments give the two
// angles, and the nearest rotation U V^T = rot(phi - psi) is just atan2(H, E) --
// no iteration and no eigenvalue solve anywhere.
struct SVD2 {
    double s0 = 1.0, s1 = 1.0;   // s0 >= s1; both positive when det A > 0
    double phi = 0.0;            // the angle of U
    double rot = 0.0;            // the angle of U V^T, the nearest rotation
};

SVD2 svd2(const double J[4]) {
    const double E = 0.5 * (J[0] + J[3]);
    const double F = 0.5 * (J[0] - J[3]);
    const double G = 0.5 * (J[2] + J[1]);
    const double H = 0.5 * (J[2] - J[1]);
    const double Q = std::hypot(E, H);
    const double R = std::hypot(F, G);

    SVD2 o;
    o.s0 = Q + R;
    o.s1 = Q - R;
    const double a1 = std::atan2(G, F);
    const double a2 = std::atan2(H, E);
    o.phi = 0.5 * (a1 + a2);
    o.rot = a2;
    return o;
}

// The symmetric Dirichlet energy of one triangle, as a function of its
// Jacobian. In two dimensions ||J^-1||_F = ||J||_F / |det J|, because the
// adjugate of a 2x2 matrix is its own transpose up to sign, so the whole thing
// is one Frobenius norm and one determinant.
inline double symmetricDirichlet(const double J[4], double det) {
    const double nf = J[0] * J[0] + J[1] * J[1] + J[2] * J[2] + J[3] * J[3];
    return nf * (1.0 + 1.0 / (det * det));
}

// dD/dsigma for that energy, and the SLIM weight built from it.
inline double weightFor(double sigma) {
    const double d = 2.0 * sigma - 2.0 / (sigma * sigma * sigma);
    double w2;
    if (std::fabs(sigma - 1.0) > 1e-8) {
        w2 = d / (sigma - 1.0);
    } else {
        // The limit, which is the second derivative: 2 + 6 sigma^-4 -> 8.
        w2 = 2.0 + 6.0 / (sigma * sigma * sigma * sigma);
    }
    if (!(w2 > 1e-12)) w2 = 1e-12;
    if (w2 > 1e14) w2 = 1e14;
    return std::sqrt(w2);
}

// The rotation by k quarter turns as a 2x2, row-major.
//
// A note on which way round it goes, because getting it backwards is not a
// tolerance issue -- it drives the map to a *different* holonomy and the result
// is not a layout at all. Stage 4 fits each arc's transition as
//
//     psi+ = R_k psi- + t
//
// (Sec. 3.2.2: the rigid motion taking the e- child polyline to the e+ one), so
// the tangents of the two sides satisfy d+ = R_k d- and the seam term is
// ||d+ - R_k d-||^2. Eq. (17) is written with R_k^-1 under the opposite
// ordering convention for the pair; what matters is not which of the two is
// called forward but that the *same* R_k appears in the fit, in this energy,
// and at a seam crossing in the tracer -- which it does. Report::
// initialSeamResidual is the check: psi_R already satisfies Q4, so it starts at
// rounding level with this convention and at O(1) with the other.
inline void quarterTurn(int k, double Q[4]) {
    switch (((k % 4) + 4) % 4) {
        case 0:  Q[0] =  1; Q[1] =  0; Q[2] =  0; Q[3] =  1; break;
        case 1:  Q[0] =  0; Q[1] = -1; Q[2] =  1; Q[3] =  0; break;
        case 2:  Q[0] = -1; Q[1] =  0; Q[2] =  0; Q[3] = -1; break;
        default: Q[0] =  0; Q[1] =  1; Q[2] = -1; Q[3] =  0; break;
    }
}

inline double cornerAngle(const Point &a, const Point &b, const Point &c) {
    const Point u = b - a;
    const Point v = c - a;
    const double nu = normP(u), nv = normP(v);
    if (nu <= 0.0 || nv <= 0.0) return 0.0;
    double t = dotP(u, v) / (nu * nv);
    t = std::max(-1.0, std::min(1.0, t));
    return std::acos(t);
}

} // namespace

// ---------------------------------------------------------------------------
LayoutEnergy::LayoutEnergy(const Immersion &immersion, SubdomainLabels &lab)
    : LayoutEnergy(immersion, lab, Options()) {}

LayoutEnergy::LayoutEnergy(const Immersion &immersion, SubdomainLabels &lab,
                           const Options &opts)
    : imm(&immersion), labels(&lab), options(opts) {
    cm = &imm->getCutMesh();
    nV = static_cast<int>(cm->vertices.size());

    uv = imm->getUV();
    x.assign(2 * nV, 0.0);
    for (int v = 0; v < nV; ++v) { x[dofU(v)] = uv[v][0]; x[dofV(v)] = uv[v][1]; }

    updateExtent();

    if (options.pinnedVertex < 0 || options.pinnedVertex >= nV) options.pinnedVertex = 0;

    currentReference = options.reference;
    buildReference(currentReference);
    buildConstraintTerms();
}

// ---------------------------------------------------------------------------
// buildReference()  --  Sec. 3.3, final paragraph
//
// J is measured against a reference triangle, and the choice of reference is
// the choice of what "undistorted" means. Two are available and the paper
// alternates between them: the Euclidean geometry the surface arrived with, and
// the Ricci metric psi_R was unfolded from.
//
// Either way the reference is rebuilt from its three side lengths rather than
// copied from any coordinates, and always laid down with a positive
// determinant. That makes "det J > 0" and "the image triangle is not inverted"
// the same statement, which is what the line search's flip cap relies on.
// ---------------------------------------------------------------------------
void LayoutEnergy::buildReference(Reference ref) {
    const Mesh &om = imm->getOriginalMesh();
    const std::vector<double> &flat = imm->getFlatEdgeLengths();
    const int nT = static_cast<int>(cm->triangles.size());

    refM.assign(nT, {0.0, 0.0, 0.0, 0.0});
    refArea.assign(nT, 0.0);

    int degenerate = 0;
    for (int t = 0; t < nT; ++t) {
        double l01, l12, l20;
        if (ref == Reference::Ricci) {
            l01 = flat[om.triangleEdges[t][0]];
            l12 = flat[om.triangleEdges[t][1]];
            l20 = flat[om.triangleEdges[t][2]];
        } else {
            const Triangle &tri = cm->triangles[t];
            l01 = normP(cm->vertices[tri[1]] - cm->vertices[tri[0]]);
            l12 = normP(cm->vertices[tri[2]] - cm->vertices[tri[1]]);
            l20 = normP(cm->vertices[tri[0]] - cm->vertices[tri[2]]);
        }
        if (!(l01 > 0.0)) { l01 = 1e-12; ++degenerate; }

        const double rx = (l20 * l20 - l12 * l12 + l01 * l01) / (2.0 * l01);
        double ry2 = l20 * l20 - rx * rx;
        if (!(ry2 > 0.0)) { ry2 = (1e-6 * l01) * (1e-6 * l01); ++degenerate; }
        const double ry = std::sqrt(ry2);

        // [r1 - r0, r2 - r0] = [[l01, rx], [0, ry]], determinant l01 * ry > 0.
        const double det = l01 * ry;
        refM[t][0] = ry / det;          // M00
        refM[t][1] = -rx / det;         // M01
        refM[t][2] = 0.0;               // M10
        refM[t][3] = l01 / det;         // M11
        refArea[t] = 0.5 * det;
    }
    if (degenerate > 0) {
        std::ostringstream oss;
        oss << degenerate << " reference triangle(s) were degenerate and had to be "
            << "thickened; the distortion measured on them means little.";
        report.messages.push_back(oss.str());
    }
}

// ---------------------------------------------------------------------------
void LayoutEnergy::jacobian(const std::vector<double> &xx, int t, double J[4]) const {
    const Triangle &tri = cm->triangles[t];
    const double d00 = xx[dofU(tri[1])] - xx[dofU(tri[0])];
    const double d01 = xx[dofU(tri[2])] - xx[dofU(tri[0])];
    const double d10 = xx[dofV(tri[1])] - xx[dofV(tri[0])];
    const double d11 = xx[dofV(tri[2])] - xx[dofV(tri[0])];
    const auto &M = refM[t];
    J[0] = d00 * M[0] + d01 * M[2];
    J[1] = d00 * M[1] + d01 * M[3];
    J[2] = d10 * M[0] + d11 * M[2];
    J[3] = d10 * M[1] + d11 * M[3];
}

// ---------------------------------------------------------------------------
// buildConstraintTerms()  --  Eqs. (15) to (19)
// ---------------------------------------------------------------------------
void LayoutEnergy::buildConstraintTerms() {
    cterms.clear();

    // E2, Eq. (15). Gamma_u holds u constant, Gamma_v holds v.
    for (const auto &be : labels->boundaryEdges()) {
        if (!(be.length > 0.0) || be.label == SubdomainLabels::Align::None) continue;
        const int c = (be.label == SubdomainLabels::Align::U) ? 0 : 1;
        Term t;
        t.which = 2;
        t.gweight = 1.0 / be.length;
        t.coeffs = {{2 * be.a + c, 1.0}, {2 * be.b + c, -1.0}};
        cterms.push_back(std::move(t));
    }

    // E3, Eq. (16). Same form, one label for the whole chain.
    for (const auto &fc : labels->featureChains()) {
        if (fc.label == SubdomainLabels::Align::None) continue;
        const int c = (fc.label == SubdomainLabels::Align::U) ? 0 : 1;
        for (size_t i = 0; i + 1 < fc.verts.size(); ++i) {
            const double L = fc.lengths[i];
            if (!(L > 0.0)) continue;
            Term t;
            t.which = 3;
            t.gweight = 1.0 / L;
            t.coeffs = {{2 * fc.verts[i] + c, 1.0}, {2 * fc.verts[i + 1] + c, -1.0}};
            cterms.push_back(std::move(t));
        }
    }

    // E4, Eq. (17). Written on tangents, so only the rotation of the
    // transition appears and its translation drops out -- which is why this
    // term can hold Q4 without ever having to know where the seam went.
    const auto &pairs = labels->seamPairs();
    const auto &pairArc = labels->seamPairArc();
    const auto &pairLen = labels->seamPairLength();
    const auto &arcList = labels->arcs();
    for (size_t i = 0; i < pairs.size(); ++i) {
        const double L = pairLen[i];
        if (!(L > 0.0)) continue;
        double Q[4];
        quarterTurn(arcList[pairArc[i]].k, Q);

        const int pa = pairs[i][0], pb = pairs[i][1];
        const int ma = pairs[i][2], mb = pairs[i][3];
        for (int p = 0; p < 2; ++p) {
            Term t;
            t.which = 4;
            t.gweight = 1.0 / L;
            t.coeffs.push_back({2 * pb + p, 1.0});
            t.coeffs.push_back({2 * pa + p, -1.0});
            for (int a = 0; a < 2; ++a) {
                const double q = Q[2 * p + a];
                if (q == 0.0) continue;
                t.coeffs.push_back({2 * mb + a, -q});
                t.coeffs.push_back({2 * ma + a, q});
            }
            cterms.push_back(std::move(t));
        }
    }

    // E5, Eqs. (18) and (19). One term per path of Gamma_topo: the total
    // advance of the constrained coordinate, accumulated across every subcurve
    // and every cut, which is zero exactly when the two endpoints sit on a
    // common isoline. A path can revisit a vertex, so the coefficients are
    // accumulated before being written out.
    for (const auto &tp : labels->topoPaths()) {
        std::map<int, double> acc;
        auto addSite = [&](const SubdomainLabels::Site &s, double sign, const Point &e) {
            for (int c = 0; c < 2; ++c) {
                if (e[c] == 0.0) continue;
                if (s.b < 0) {
                    acc[2 * s.a + c] += sign * e[c];
                } else {
                    acc[2 * s.a + c] += sign * e[c] * (1.0 - s.t);
                    acc[2 * s.b + c] += sign * e[c] * s.t;
                }
            }
        };
        for (const auto &sub : tp.subs) {
            const Point e = SubdomainLabels::axis(sub.dir);
            addSite(sub.to, 1.0, e);
            addSite(sub.from, -1.0, e);
        }
        Term t;
        t.which = 5;
        t.gweight = 1.0;
        for (const auto &kv : acc) {
            if (kv.second != 0.0) t.coeffs.push_back({kv.first, kv.second});
        }
        if (t.coeffs.empty()) continue;
        cterms.push_back(std::move(t));
    }

    // E6, the interface sectors. Written on tangents in the same way E4 is, so
    // that the node's own position drops out and only the turn between the two
    // rays is held:  b / l_b  =  R_q ( a / l_a ).
    //
    // The two rays share the node, so its dofs appear in both halves and the
    // coefficients are accumulated before being written out -- and where q is
    // even they partly cancel, which is exactly right and is why the sum has to
    // be taken before the term is built rather than after.
    if (options.lagInterfaceScales && uv.size() == static_cast<size_t>(nV)) {
        labels->refreshInterfaceScales(uv);
    }
    e6Start = cterms.size();
    for (const auto &ic : labels->interfaceCorners()) {
        if (!(ic.sa > 0.0) || !(ic.sb > 0.0)) continue;
        double Q[4];
        quarterTurn(ic.quarters, Q);
        for (int p = 0; p < 2; ++p) {
            std::map<int, double> acc;
            acc[2 * ic.b + p] += 1.0 / ic.sb;
            acc[2 * ic.vertex + p] -= 1.0 / ic.sb;
            for (int c = 0; c < 2; ++c) {
                const double q = Q[2 * p + c];
                if (q == 0.0) continue;
                acc[2 * ic.a + c] -= q / ic.sa;
                acc[2 * ic.vertex + c] += q / ic.sa;
            }
            Term t;
            t.which = 6;
            t.gweight = 0.5 * (ic.la + ic.lb);
            for (const auto &kv : acc) {
                if (kv.second != 0.0) t.coeffs.push_back({kv.first, kv.second});
            }
            if (t.coeffs.empty()) continue;
            cterms.push_back(std::move(t));
        }
    }

    e6Count = cterms.size() - e6Start;

    // The proxy's pattern is the union of the mesh's blocks and these terms',
    // so it changes exactly when they do and has to be rebuilt with them.
    buildHessianPattern();
}

// ---------------------------------------------------------------------------
// refreshInterfaceTerms()
// ---------------------------------------------------------------------------
void LayoutEnergy::refreshInterfaceTerms() {
    if (!options.lagInterfaceScales) return;
    if (e6Count == 0 || uv.size() != static_cast<size_t>(nV)) return;
    labels->refreshInterfaceScales(uv);

    size_t at = e6Start;
    for (const auto &ic : labels->interfaceCorners()) {
        if (!(ic.sa > 0.0) || !(ic.sb > 0.0)) continue;
        double Q[4];
        quarterTurn(ic.quarters, Q);
        for (int p = 0; p < 2; ++p) {
            std::map<int, double> acc;
            acc[2 * ic.b + p] += 1.0 / ic.sb;
            acc[2 * ic.vertex + p] -= 1.0 / ic.sb;
            for (int c = 0; c < 2; ++c) {
                const double q = Q[2 * p + c];
                if (q == 0.0) continue;
                acc[2 * ic.a + c] -= q / ic.sa;
                acc[2 * ic.vertex + c] += q / ic.sa;
            }
            std::vector<std::pair<int, double>> coeffs;
            for (const auto &kv : acc) {
                if (kv.second != 0.0) coeffs.push_back({kv.first, kv.second});
            }
            if (coeffs.empty()) continue;
            if (at >= e6Start + e6Count) return;
            // The dofs are the same ones as when the pattern was built -- the
            // node and the far end of each ray -- so only the values move. If
            // that ever stopped being true the slots would be wrong, so it is
            // checked rather than assumed.
            if (cterms[at].coeffs.size() != coeffs.size()) { ++at; continue; }
            bool same = true;
            for (size_t k = 0; k < coeffs.size(); ++k) {
                if (cterms[at].coeffs[k].first != coeffs[k].first) same = false;
            }
            if (same) cterms[at].coeffs = std::move(coeffs);
            ++at;
        }
    }
}

// ---------------------------------------------------------------------------
// energy()  --  Eq. (13)
// ---------------------------------------------------------------------------
double LayoutEnergy::energy(const std::vector<double> &xx, double *o1, double *o2,
                            double *o3, double *o4, double *o5, double *o6) const {
    double E1 = 0.0;
    double J[4];
    for (int t = 0; t < static_cast<int>(cm->triangles.size()); ++t) {
        jacobian(xx, t, J);
        const double det = J[0] * J[3] - J[1] * J[2];
        if (!(det > 0.0)) return kInf;   // outside F of Eq. (12): the barrier
        E1 += refArea[t] * symmetricDirichlet(J, det);
    }

    double E[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (const Term &t : cterms) {
        const double r = residualOf(t, xx);
        E[t.which] += t.gweight * r * r;
    }

    if (o1) *o1 = E1;
    if (o2) *o2 = E[2];
    if (o3) *o3 = E[3];
    if (o4) *o4 = E[4];
    if (o5) *o5 = E[5];
    if (o6) *o6 = E[6];

    return lambda[1] * E1 + lambda[2] * E[2] + lambda[3] * E[3] +
           lambda[4] * E[4] + lambda[5] * E[5] + lambda[6] * E[6];
}

// ---------------------------------------------------------------------------
// model()
//
// The gradient of the true objective, and the Hessian of the quadratic proxy.
// They are assembled together because for E2..E5 they come from the same
// squared linear form, and because for E1 the whole point of the SLIM weights
// is that the proxy's gradient *is* the true one -- assembling them apart would
// invite them to disagree.
// ---------------------------------------------------------------------------
bool LayoutEnergy::model(const std::vector<double> &xx, std::vector<double> &grad,
                         Eigen::SparseMatrix<double> &H) const {
    grad.assign(2 * nV, 0.0);
    double *val = H.valuePtr();
    if (val) std::fill(val, val + H.nonZeros(), 0.0);

    // --- E1, by the SLIM proxy -------------------------------------------
    //
    // The 6x6 local block is written down whole rather than assembled from the
    // four rank-one pieces one (p, q) at a time. Collecting them costs nothing
    // and is exact: with
    //
    //     c^{pq}_{ak} = W_{pa} g^q_k,   g^q = (-(M_{0q} + M_{1q}), M_{0q}, M_{1q})
    //
    // the four outer products sum to
    //
    //     sum_{pq} c^{pq}_{ak} c^{pq}_{bl} = (W^T W)_{ab} (sum_q g^q_k g^q_l)
    //
    // and W is symmetric, so the block is 2 w (W^2 (x) G): a Kronecker product
    // of a 2x2 with a 3x3, thirty-six numbers from a dozen multiplications
    // rather than a hundred and forty-four triplets from four heap allocations.
    // The same collection applies to the gradient, where the four residuals
    // W(J - R)_{pq} assemble into the single 2x2 B = W^2 (J - R).
    double J[4];
    const int nT = static_cast<int>(cm->triangles.size());
    for (int t = 0; t < nT; ++t) {
        jacobian(xx, t, J);
        const double det = J[0] * J[3] - J[1] * J[2];
        if (!(det > 0.0)) return false;

        const SVD2 s = svd2(J);
        const double w0 = weightFor(std::max(s.s0, 1e-12));
        const double w1 = weightFor(std::max(s.s1, 1e-12));

        const double c = std::cos(s.phi), sn = std::sin(s.phi);
        // W = U diag(w0, w1) U^T.
        const double W[4] = {w0 * c * c + w1 * sn * sn, (w0 - w1) * c * sn,
                             (w0 - w1) * c * sn,       w0 * sn * sn + w1 * c * c};
        // J - R, R = U V^T being the nearest rotation.
        const double cr = std::cos(s.rot), sr = std::sin(s.rot);
        const double D[4] = {J[0] - cr, J[1] + sr, J[2] - sr, J[3] - cr};

        const double W2[4] = {W[0] * W[0] + W[1] * W[2], W[0] * W[1] + W[1] * W[3],
                              W[2] * W[0] + W[3] * W[2], W[2] * W[1] + W[3] * W[3]};
        const double B[4] = {W2[0] * D[0] + W2[1] * D[2], W2[0] * D[1] + W2[1] * D[3],
                             W2[2] * D[0] + W2[3] * D[2], W2[2] * D[1] + W2[3] * D[3]};

        const auto &M = refM[t];
        const double g[2][3] = {{-(M[0] + M[2]), M[0], M[2]},
                                {-(M[1] + M[3]), M[1], M[3]}};
        double G[3][3];
        for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) G[k][l] = g[0][k] * g[0][l] + g[1][k] * g[1][l];
        }

        const Triangle &tri = cm->triangles[t];
        // The proxy carries a factor of one half so that its gradient is the
        // true gradient of A * D(sigma) and not twice it; without that, E1
        // would silently outweigh E2..E5 by a factor of two.
        const double w = 2.0 * (0.5 * lambda[1] * refArea[t]);

        for (int a = 0; a < 2; ++a) {
            for (int k = 0; k < 3; ++k) {
                grad[2 * tri[k] + a] += w * (B[2 * a] * g[0][k] + B[2 * a + 1] * g[1][k]);
            }
        }

        const int *slot = triSlot.data() + 36 * t;
        for (int a = 0; a < 2; ++a) {
            for (int b = 0; b < 2; ++b) {
                const double wab = w * W2[2 * a + b];
                if (wab == 0.0) continue;
                for (int k = 0; k < 3; ++k) {
                    const int *row = slot + (a * 3 + k) * 6 + b * 3;
                    for (int l = 0; l < 3; ++l) {
                        if (row[l] >= 0) val[row[l]] += wab * G[k][l];
                    }
                }
            }
        }
    }

    // --- E2 to E5, exactly ------------------------------------------------
    for (size_t i = 0; i < cterms.size(); ++i) {
        const Term &t = cterms[i];
        const double w2 = lambda[t.which] * t.gweight;
        const double res = residualOf(t, xx);
        const int k = static_cast<int>(t.coeffs.size());
        const int *slot = termSlot.data() + termSlotStart[i];
        for (int p = 0; p < k; ++p) grad[t.coeffs[p].first] += 2.0 * w2 * res * t.coeffs[p].second;
        for (int p = 0; p < k; ++p) {
            const double cp = 2.0 * w2 * t.coeffs[p].second;
            const int *row = slot + p * k;
            for (int q = 0; q < k; ++q) {
                if (row[q] >= 0) val[row[q]] += cp * t.coeffs[q].second;
            }
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// buildHessianPattern()
//
// Where every local entry of the proxy Hessian lands, worked out once. The
// pattern is the union of a 6x6 block per triangle -- (u, v) at its three
// vertices -- and a k x k block per constraint term over its own coefficient
// list, with the two pinned dofs struck out and the diagonal added whether or
// not anything writes to it. Nothing in that list depends on x, on lambda or on
// the SLIM weights, which is why it can be a fixed array scattered into rather
// than a triplet list rebuilt, sorted and collapsed at every inner iteration.
//
// The diagonal is there for innerSolve()'s Tikhonov fallback: a shift added to
// entries that already exist leaves the pattern alone, and leaving the pattern
// alone is what lets the symbolic factorisation be computed once per constraint
// set instead of once per solve.
// ---------------------------------------------------------------------------
void LayoutEnergy::buildHessianPattern() {
    reducedDof.assign(2 * nV, -1);
    nReduced = 0;
    for (int d = 0; d < 2 * nV; ++d) {
        if (d == dofU(options.pinnedVertex) || d == dofV(options.pinnedVertex)) continue;
        reducedDof[d] = nReduced++;
    }

    const int nT = static_cast<int>(cm->triangles.size());
    triSlot.assign(static_cast<size_t>(36) * nT, -1);

    termSlotStart.assign(cterms.size() + 1, 0);
    size_t total = 0;
    for (size_t i = 0; i < cterms.size(); ++i) {
        termSlotStart[i] = static_cast<int>(total);
        total += cterms[i].coeffs.size() * cterms[i].coeffs.size();
    }
    termSlotStart.back() = static_cast<int>(total);
    termSlot.assign(total, -1);

    analysed = false;
    diagSlot.assign(std::max(0, nReduced), -1);
    if (nReduced == 0) { Hproxy.resize(0, 0); return; }

    std::vector<Eigen::Triplet<double>> pat;
    pat.reserve(static_cast<size_t>(nReduced) + static_cast<size_t>(36) * nT + total);
    for (int i = 0; i < nReduced; ++i) pat.emplace_back(i, i, 0.0);

    std::vector<int> ld(6, -1);
    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = cm->triangles[t];
        for (int a = 0; a < 2; ++a) {
            for (int k = 0; k < 3; ++k) ld[a * 3 + k] = reducedDof[2 * tri[k] + a];
        }
        for (int r = 0; r < 6; ++r) {
            if (ld[r] < 0) continue;
            for (int q = 0; q < 6; ++q) {
                if (ld[q] < 0) continue;
                pat.emplace_back(ld[r], ld[q], 0.0);
            }
        }
    }
    for (const Term &t : cterms) {
        for (const auto &ci : t.coeffs) {
            const int ri = reducedDof[ci.first];
            if (ri < 0) continue;
            for (const auto &cj : t.coeffs) {
                const int rj = reducedDof[cj.first];
                if (rj < 0) continue;
                pat.emplace_back(ri, rj, 0.0);
            }
        }
    }

    Hproxy.resize(nReduced, nReduced);
    Hproxy.setFromTriplets(pat.begin(), pat.end());
    Hproxy.makeCompressed();

    const int *outer = Hproxy.outerIndexPtr();
    const int *inner = Hproxy.innerIndexPtr();
    auto slotOf = [&](int r, int col) -> int {
        const int *b = inner + outer[col], *e = inner + outer[col + 1];
        const int *it = std::lower_bound(b, e, r);
        return (it != e && *it == r) ? static_cast<int>(it - inner) : -1;
    };

    for (int i = 0; i < nReduced; ++i) diagSlot[i] = slotOf(i, i);

    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = cm->triangles[t];
        for (int a = 0; a < 2; ++a) {
            for (int k = 0; k < 3; ++k) ld[a * 3 + k] = reducedDof[2 * tri[k] + a];
        }
        int *slot = triSlot.data() + 36 * t;
        for (int r = 0; r < 6; ++r) {
            for (int q = 0; q < 6; ++q) {
                slot[r * 6 + q] = (ld[r] < 0 || ld[q] < 0) ? -1 : slotOf(ld[r], ld[q]);
            }
        }
    }
    for (size_t i = 0; i < cterms.size(); ++i) {
        const Term &t = cterms[i];
        const int k = static_cast<int>(t.coeffs.size());
        int *slot = termSlot.data() + termSlotStart[i];
        for (int p = 0; p < k; ++p) {
            const int ri = reducedDof[t.coeffs[p].first];
            for (int q = 0; q < k; ++q) {
                const int rj = reducedDof[t.coeffs[q].first];
                slot[p * k + q] = (ri < 0 || rj < 0) ? -1 : slotOf(ri, rj);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// maxStep()  --  Sec. 3.3's flip-free line search
//
// Along x + s d the two edge vectors of every triangle move linearly, so
// det J_t(s) is a quadratic in s and its smallest positive root is where that
// triangle would invert. The cap is the smallest such root over the mesh.
// ---------------------------------------------------------------------------
double LayoutEnergy::maxStep(const std::vector<double> &xx,
                             const std::vector<double> &d) const {
    double smax = kInf;
    for (int t = 0; t < static_cast<int>(cm->triangles.size()); ++t) {
        const Triangle &tri = cm->triangles[t];
        const double D00 = xx[dofU(tri[1])] - xx[dofU(tri[0])];
        const double D01 = xx[dofU(tri[2])] - xx[dofU(tri[0])];
        const double D10 = xx[dofV(tri[1])] - xx[dofV(tri[0])];
        const double D11 = xx[dofV(tri[2])] - xx[dofV(tri[0])];
        const double E00 = d[dofU(tri[1])] - d[dofU(tri[0])];
        const double E01 = d[dofU(tri[2])] - d[dofU(tri[0])];
        const double E10 = d[dofV(tri[1])] - d[dofV(tri[0])];
        const double E11 = d[dofV(tri[2])] - d[dofV(tri[0])];

        const double a = E00 * E11 - E01 * E10;
        const double b = D00 * E11 + E00 * D11 - D01 * E10 - E01 * D10;
        const double c = D00 * D11 - D01 * D10;
        if (!(c > 0.0)) return 0.0;   // already inverted; nothing to step

        double root = kInf;
        if (std::fabs(a) < 1e-300) {
            if (b < 0.0) root = -c / b;
        } else {
            const double disc = b * b - 4.0 * a * c;
            if (disc >= 0.0) {
                const double sq = std::sqrt(disc);
                // The numerically stable pair, then the smallest positive one.
                const double q = -0.5 * (b + (b >= 0.0 ? sq : -sq));
                const double r0 = (q != 0.0) ? q / a : kInf;
                const double r1 = (q != 0.0) ? c / q : kInf;
                for (double r : {r0, r1}) if (r > 0.0) root = std::min(root, r);
            }
        }
        smax = std::min(smax, root);
    }
    return smax;
}

// ---------------------------------------------------------------------------
// innerSolve()
// ---------------------------------------------------------------------------
bool LayoutEnergy::innerSolve(int outer) {
    (void)outer;
    // options.innerTolerance is a step length relative to the image, so the
    // stopping test below needs the extent of the image this solve starts from,
    // not psi_R's.
    updateExtent();
    if (nReduced == 0) return false;

    std::vector<double> grad, dir(2 * nV, 0.0), trial(2 * nV, 0.0);
    Eigen::VectorXd rhs(nReduced), sol;

    bool moved = false;
    for (int it = 0; it < options.innerIterations; ++it) {
        if (!model(x, grad, Hproxy)) {
            report.messages.push_back("The map inverted before the model could be built.");
            return moved;
        }

        for (int d = 0; d < 2 * nV; ++d) if (reducedDof[d] >= 0) rhs[reducedDof[d]] = -grad[d];

        // The ordering and the elimination tree come from the pattern, which
        // has not changed since buildHessianPattern() and will not change until
        // the constraint set does; only the numbers have.
        if (!analysed) { ldlt.analyzePattern(Hproxy); analysed = true; }

        double *val = Hproxy.valuePtr();
        bool solved = false;
        double reg = 0.0, applied = 0.0;
        for (int attempt = 0; attempt < 4 && !solved; ++attempt) {
            if (reg != applied) {
                const double delta = reg - applied;
                for (int i = 0; i < nReduced; ++i) {
                    if (diagSlot[i] >= 0) val[diagSlot[i]] += delta;
                }
                applied = reg;
            }
            ldlt.factorize(Hproxy);
            if (ldlt.info() == Eigen::Success) {
                sol = ldlt.solve(rhs);
                if (ldlt.info() == Eigen::Success && sol.allFinite()) solved = true;
            }
            if (!solved) {
                ++report.factorisationFallbacks;
                // The proxy Hessian is positive semi-definite by construction,
                // so a failure here is conditioning rather than indefiniteness;
                // a Tikhonov shift is the right response and biases the
                // direction towards gradient descent, which the line search
                // then makes safe. Added on the diagonal the pattern already
                // carries, so the symbolic factorisation still holds.
                double biggest = 1.0;
                const int nz = static_cast<int>(Hproxy.nonZeros());
                for (int i = 0; i < nz; ++i) biggest = std::max(biggest, std::fabs(val[i]));
                reg = (reg == 0.0) ? 1e-9 * biggest : reg * 100.0;
            }
        }
        if (!solved) {
            report.messages.push_back("The proxy system could not be factorised.");
            return moved;
        }

        std::fill(dir.begin(), dir.end(), 0.0);
        for (int d = 0; d < 2 * nV; ++d) if (reducedDof[d] >= 0) dir[d] = sol[reducedDof[d]];

        double slope = 0.0, dnorm = 0.0;
        for (int d = 0; d < 2 * nV; ++d) {
            slope += grad[d] * dir[d];
            dnorm = std::max(dnorm, std::fabs(dir[d]));
        }
        if (dnorm < options.innerTolerance * extent) break;   // converged
        if (!(slope < 0.0)) {
            // Cannot happen with a positive semi-definite proxy and a nonzero
            // step, but a regularised solve that stopped early can produce it.
            for (int d = 0; d < 2 * nV; ++d) dir[d] = (reducedDof[d] >= 0) ? -grad[d] : 0.0;
            slope = 0.0;
            for (int d = 0; d < 2 * nV; ++d) slope += grad[d] * dir[d];
            if (!(slope < 0.0)) break;
        }

        const double cap = maxStep(x, dir);
        double s = options.stepSafety * std::min(1.0, cap);
        if (!(s > 0.0)) { ++report.lineSearchFailures; break; }

        const double E0 = energy(x);
        bool stepped = false;
        double accepted = E0;
        for (int back = 0; back < options.maxBacktracks; ++back) {
            for (int d = 0; d < 2 * nV; ++d) trial[d] = x[d] + s * dir[d];
            const double E = energy(trial);
            if (E <= E0 + 1e-4 * s * slope) {
                x.swap(trial);
                accepted = E;
                stepped = true;
                moved = true;
                break;
            }
            s *= 0.5;
        }
        ++report.innerIterations;
        if (!stepped) { ++report.lineSearchFailures; break; }
        // The step-size test above is on the *map*; this one is on the energy.
        // Near the barrier the proxy takes many short steps that each shrink
        // the map by a little and the objective by nothing, and there is no
        // point spending the rest of the budget on them -- the next level of
        // lambda is worth far more than another hundred of these.
        if (E0 - accepted <= 1e-12 * std::max(1.0, std::fabs(E0))) break;
    }

    for (int v = 0; v < nV; ++v) uv[v] = Point{x[dofU(v)], x[dofV(v)]};
    return moved;
}

// ---------------------------------------------------------------------------
double LayoutEnergy::imageExtent() const {
    Point lo{kInf, kInf}, hi{-kInf, -kInf};
    for (const Point &p : uv) {
        lo[0] = std::min(lo[0], p[0]);
        lo[1] = std::min(lo[1], p[1]);
        hi[0] = std::max(hi[0], p[0]);
        hi[1] = std::max(hi[1], p[1]);
    }
    if (!(hi[0] >= lo[0])) return 1.0;   // no vertices
    return std::max(1e-30, std::hypot(hi[0] - lo[0], hi[1] - lo[1]));
}

// ---------------------------------------------------------------------------
// measure()
//
// One residual per property, read off the map rather than off the energies --
// a small E2 with a large lambda_2 says nothing about how far Q3 still is.
// ---------------------------------------------------------------------------
void LayoutEnergy::measure(bool initial) {
    // The residuals below are all lengths divided by the extent, and the extent
    // moves with the map. Re-read it before using it.
    updateExtent();

    report.maxBoundaryResidual = 0.0;
    report.maxFeatureResidual = 0.0;
    report.maxSeamResidual = 0.0;
    report.maxTopoResidual = 0.0;
    report.invertedTriangles = 0;
    report.minDetJ = kInf;

    // Q3.
    for (const auto &be : labels->boundaryEdges()) {
        if (be.label == SubdomainLabels::Align::None) continue;
        const int c = (be.label == SubdomainLabels::Align::U) ? 0 : 1;
        report.maxBoundaryResidual = std::max(report.maxBoundaryResidual,
                                              std::fabs(uv[be.a][c] - uv[be.b][c]) / extent);
    }

    // Features.
    for (const auto &fc : labels->featureChains()) {
        if (fc.label == SubdomainLabels::Align::None) continue;
        const int c = (fc.label == SubdomainLabels::Align::U) ? 0 : 1;
        for (size_t i = 0; i + 1 < fc.verts.size(); ++i) {
            report.maxFeatureResidual =
                std::max(report.maxFeatureResidual,
                         std::fabs(uv[fc.verts[i]][c] - uv[fc.verts[i + 1]][c]) / extent);
        }
    }

    // Q4.
    const auto &pairs = labels->seamPairs();
    const auto &pairArc = labels->seamPairArc();
    const auto &arcList = labels->arcs();
    for (size_t i = 0; i < pairs.size(); ++i) {
        const Point dPlus = uv[pairs[i][1]] - uv[pairs[i][0]];
        const Point dMinus = uv[pairs[i][3]] - uv[pairs[i][2]];
        const int k = arcList[pairArc[i]].k;
        const Point mapped = Immersion::rotateQuarter(dMinus, k);
        report.maxSeamResidual = std::max(report.maxSeamResidual,
                                          normP(dPlus - mapped) / extent);
    }

    // Q5.
    for (const auto &tp : labels->topoPaths()) {
        report.maxTopoResidual = std::max(report.maxTopoResidual,
                                          std::fabs(labels->topoResidual(uv, tp)) / extent);
    }

    // E6: how far each sector of the interface network is from the whole number
    // of quarter turns the model has there.
    //
    // Measured as an *angle* and not as the length E6 minimises. The energy is
    // written on b/l_b - R_q(a/l_a), whose norm also carries the difference in
    // how much the map stretches the two rays; that part is E1's business and
    // is not a defect -- an anisotropic layout is still a layout -- so folding
    // it into the verdict would fail models that are perfectly aligned and
    // merely not conformal. The angle between the two tangents is the whole of
    // what the sector condition says and it is scale-free.
    report.maxInterfaceResidual = 0.0;
    report.interfaceCornerChanges = 0;
    report.interfaceCorners = static_cast<int>(labels->interfaceCorners().size());
    for (const auto &ic : labels->interfaceCorners()) {
        const Point a = uv[ic.a] - uv[ic.vertex];
        const Point b = uv[ic.b] - uv[ic.vertex];
        if (!(normP(a) > 0.0) || !(normP(b) > 0.0)) continue;
        double turn = std::atan2(cross2(a, b), dotP(a, b));   // in (-pi, pi]
        if (turn < 0.0) turn += 2.0 * M_PI;                   // in [0, 2pi)
        const double want = M_PI_2 * ic.quarters;
        double res = turn - want;
        while (res > M_PI) res -= 2.0 * M_PI;
        while (res < -M_PI) res += 2.0 * M_PI;
        report.maxInterfaceResidual = std::max(report.maxInterfaceResidual, std::fabs(res));
        if (std::lround(turn / M_PI_2) != (ic.quarters % 4)) ++report.interfaceCornerChanges;
    }

    // Q1, and the angle sums for Q2.
    std::vector<double> angleAt(nV, 0.0);
    double J[4];
    for (int t = 0; t < static_cast<int>(cm->triangles.size()); ++t) {
        jacobian(x, t, J);
        const double det = J[0] * J[3] - J[1] * J[2];
        report.minDetJ = std::min(report.minDetJ, det);
        if (!(det > 0.0)) ++report.invertedTriangles;

        const Triangle &tri = cm->triangles[t];
        angleAt[tri[0]] += cornerAngle(uv[tri[0]], uv[tri[1]], uv[tri[2]]);
        angleAt[tri[1]] += cornerAngle(uv[tri[1]], uv[tri[2]], uv[tri[0]]);
        angleAt[tri[2]] += cornerAngle(uv[tri[2]], uv[tri[0]], uv[tri[1]]);
    }
    if (!std::isfinite(report.minDetJ)) report.minDetJ = 0.0;

    // Q2: the angle sum of a vertex of S, taken over all its children in Omega
    // -- which is how Eq. (1) is written, so that a vertex the cutting graph
    // ran across is read once rather than once per sector.
    //
    // The comparison is against the angle Stage 1 *prescribed*, 2pi (or pi on
    // dS) less (pi/2) I(v), not against whichever multiple of pi/2 happens to
    // be nearest. The difference matters: a cone that drifted from 3pi/2 to
    // 2pi is still "a multiple of pi/2" and would pass the lax test, while
    // being a cone that has quietly changed valence and a layout that no longer
    // has the patch structure Stage 1 asked for. coneValenceChanges counts
    // exactly those.
    report.maxConeAngleResidual = 0.0;
    report.maxRegularAngleResidual = 0.0;
    report.coneValenceChanges = 0;
    report.worstConeAngleVertex = -1;
    {
        const Mesh &om = imm->getOriginalMesh();
        const auto &o2c = imm->getCut().getOriginalToCutVertices();
        const std::vector<int> &index = imm->getCones().getIndices();
        const std::vector<char> &active = imm->getCones().getActiveVertices();
        for (int v = 0; v < static_cast<int>(om.vertices.size()); ++v) {
            if (!active[v] || o2c[v].empty()) continue;
            double sum = 0.0;
            for (int cv : o2c[v]) sum += angleAt[cv];
            const double full = om.isBoundaryVertex[v] ? M_PI : 2.0 * M_PI;
            const double want = full - M_PI_2 * index[v];
            const double res = std::fabs(sum - want);
            if (index[v] != 0) {
                if (res > report.maxConeAngleResidual) report.worstConeAngleVertex = v;
                report.maxConeAngleResidual = std::max(report.maxConeAngleResidual, res);
                if (std::fabs(sum - M_PI_2 * std::round(sum / M_PI_2)) < 0.25 &&
                    std::lround(sum / M_PI_2) != std::lround(want / M_PI_2)) {
                    ++report.coneValenceChanges;
                }
            } else {
                report.maxRegularAngleResidual = std::max(report.maxRegularAngleResidual, res);
            }
        }
    }

    if (initial) {
        report.initialBoundaryResidual = report.maxBoundaryResidual;
        report.initialSeamResidual = report.maxSeamResidual;
        report.initialTopoResidual = report.maxTopoResidual;
    }
}

bool LayoutEnergy::constraintsUnder(double tol) const {
    return report.maxBoundaryResidual < tol && report.maxFeatureResidual < tol &&
           report.maxSeamResidual < tol && report.maxTopoResidual < tol;
}

// ---------------------------------------------------------------------------
// run()  --  Sec. 3.3's penalty continuation
//
//     lambda_1 <- 1                            fixed
//     lambda_2..5 <- lambda_init
//     for outer = 1 .. N:
//         minimise sum lambda_j E_j
//         lambda_j <- 10 lambda_j,  j = 2..5
//         optionally switch J's reference metric
//         optionally relabel Gamma_u / Gamma_v
//
// The two "optionally" steps are taken only when the previous outer step did
// not move the constraints, because they are the paper's escape from a stall
// and not part of the schedule: switching the reference at every level would
// keep restarting the distortion term from a different notion of undistorted
// and never let the continuation settle.
// ---------------------------------------------------------------------------
bool LayoutEnergy::run() {
    lambda[1] = options.lambda1;
    for (int j = 2; j <= 6; ++j) lambda[j] = options.lambdaInit;
    // E4 is the one constraint psi_R already satisfies exactly, so its penalty
    // has a different job from the others: not to *reach* Q4 over the course of
    // the continuation but to hold it while E2, E3 and E5 drag the map around.
    // Left at lambda_init it cannot do that job -- the first few levels then
    // buy progress elsewhere by trading Q4 away, and Q4 is not free to trade,
    // because the cone angle sums of Q2 are the same statement in different
    // units. The angular form of E4's residual is ||d+ - R_k d-|| / l_e, and
    // the l_e^-1 weighting of Eq. (17) makes that *largest* on the short seam
    // edges, which are the ones at the cones. So E4 starts level with the
    // distortion term and grows from there.
    lambda[4] = std::max(options.lambdaInit, options.lambda1);

    measure(true);
    report.energyStart = energy(x);
    return continuation(options.outerSteps);
}

// ---------------------------------------------------------------------------
void LayoutEnergy::finaliseVerdict() {
    report.injective = report.invertedTriangles == 0;
    report.constraintsMet = constraintsUnder(options.constraintTolerance);
    report.anglesHeld = report.coneValenceChanges == 0 &&
                        report.maxConeAngleResidual < options.angleTolerance;
    report.valid = report.injective && report.constraintsMet && report.anglesHeld;
}

// ---------------------------------------------------------------------------
void LayoutEnergy::loadMap(const std::vector<Point> &m) {
    if (m.size() != uv.size()) return;
    uv = m;
    for (int v = 0; v < nV; ++v) { x[dofU(v)] = uv[v][0]; x[dofV(v)] = uv[v][1]; }
    updateExtent();
    measure(false);
    finaliseVerdict();
}

// ---------------------------------------------------------------------------
void LayoutEnergy::rebuildConstraints() {
    buildConstraintTerms();
    measure(false);
    finaliseVerdict();
}

// ---------------------------------------------------------------------------
// resume()  --  Sec. 3.3's "raise lambda_5 and re-run from the current phi"
// ---------------------------------------------------------------------------
bool LayoutEnergy::resume(int outerSteps, double lambdaBoost) {
    rebuildConstraints();
    if (lambdaBoost > 0.0) {
        for (int j = 2; j <= 6; ++j) lambda[j] *= lambdaBoost;
    }
    ++report.resumes;
    return continuation(outerSteps);
}

// ---------------------------------------------------------------------------
bool LayoutEnergy::continuation(int outerSteps) {
    double previousWorst = kInf;
    double previousAngle = kInf;
    bool stalled = false;
    for (int outer = 0; outer < outerSteps; ++outer) {
        // The schedule is applied *before* the solve, not after, so that the
        // last minimisation is the one at the largest lambda rather than one
        // level below it -- and so that Report::lambdaFinal is the lambda the
        // returned map was actually produced at.
        if (outer > 0) {
            for (int j = 2; j <= 6; ++j) lambda[j] *= options.lambdaGrowth;
            if (stalled && options.alternateReference) {
                currentReference = (currentReference == Reference::Ricci)
                                       ? Reference::Euclidean : Reference::Ricci;
                buildReference(currentReference);
                ++report.referenceSwitches;
            }
            if (stalled && options.relabel) {
                labels->relabel(uv);
                buildConstraintTerms();
                ++report.relabels;
            } else {
                // E6's normalisation is lagged, so it is re-read at every level
                // and not only when the continuation stalls. The relabel branch
                // above has already done it as part of rebuilding the terms.
                refreshInterfaceTerms();
            }
        }

        innerSolve(outer);
        ++report.outerSteps;
        measure(false);

        // Every property, not just the ones the penalties are written on.
        // Stopping at "Q3, Q4 and Q5 are under tolerance" leaves Q2 wherever
        // it happened to be, and Q2 is the one that lags: its residual is the
        // others divided by the *local* edge length, so a model with a short
        // boundary edge at a cone still needs several more levels of lambda
        // after the length residuals have arrived.
        if (constraintsUnder(options.constraintTolerance) &&
            report.coneValenceChanges == 0 &&
            report.maxConeAngleResidual < options.angleTolerance &&
            report.interfaceCornerChanges == 0 &&
            report.maxInterfaceResidual < options.angleTolerance) {
            break;
        }
        if (report.invertedTriangles > 0) {
            report.messages.push_back("The map inverted during the continuation; stopping.");
            break;
        }

        const double worst = std::max(std::max(report.maxBoundaryResidual,
                                               report.maxFeatureResidual),
                                      std::max(report.maxSeamResidual,
                                               report.maxTopoResidual));
        stalled = worst > 0.9 * previousWorst &&
                  report.maxConeAngleResidual > 0.9 * previousAngle;
        previousAngle = report.maxConeAngleResidual;
        previousWorst = worst;
    }

    report.energyEnd = energy(x, &report.e1, &report.e2, &report.e3, &report.e4,
                              &report.e5, &report.e6);
    for (int j = 2; j <= 6; ++j) report.lambdaFinal[j - 2] = lambda[j];

    measure(false);
    finaliseVerdict();

    // Both diagnostics below are about the state this call left the map in, and
    // resume() may be about to change it, so they say which round they came
    // from rather than reading as the last word.
    const std::string round = report.resumes == 0
        ? std::string("The continuation")
        : ("Repair pass " + std::to_string(report.resumes) + " of the continuation");

    if (!report.constraintsMet) {
        std::ostringstream oss;
        oss << round << " stopped with constraints unmet: Q3 " << report.maxBoundaryResidual
            << ", features " << report.maxFeatureResidual
            << ", Q4 " << report.maxSeamResidual
            << ", Q5 " << report.maxTopoResidual
            << " (relative, tolerance " << options.constraintTolerance << ").";
        report.messages.push_back(oss.str());
    }
    if (!report.anglesHeld) {
        std::ostringstream oss;
        oss << "Q2 does not hold: the cone at vertex " << report.worstConeAngleVertex
            << " came out " << report.maxConeAngleResidual
            << " rad away from the angle Stage 1 prescribed for it";
        if (report.coneValenceChanges > 0) {
            oss << ", and " << report.coneValenceChanges
                << " cone(s) settled on a different multiple of pi/2 altogether";
        }
        oss << ".";
        if (imm->getReport().clusteredConePairs > 0) {
            oss << " " << imm->getReport().clusteredConePairs
                << " pair(s) of cones are clustered (closest "
                << imm->getReport().minConeSeparation
                << " of the model apart), which is the usual cause: the boundary edge between "
                << "two adjacent cones is far shorter than the rest, and Eq. (15)'s l_e^-1 "
                << "weighting leaves its angular error correspondingly larger. The remedy is "
                << "Sec. 3.1's, in Stage 1: merge the cluster.";
        }
        report.messages.push_back(oss.str());
    }
    return report.valid;
}

// ---------------------------------------------------------------------------
std::vector<double> LayoutEnergy::determinants() const {
    std::vector<double> out(cm->triangles.size(), 0.0);
    double J[4];
    for (int t = 0; t < static_cast<int>(cm->triangles.size()); ++t) {
        jacobian(x, t, J);
        out[t] = J[0] * J[3] - J[1] * J[2];
    }
    return out;
}

// ---------------------------------------------------------------------------
// checkGradient()
//
// The comparison is deliberately *not* made at the current iterate. psi_R is an
// isometry of the Ricci metric, so at the start of the continuation every J is
// the identity, every singular value is 1, dD/dsigma is 0 and the whole
// gradient vanishes -- and a relative comparison against a gradient of zero
// measures nothing but the roundoff of the difference quotient. So the map is
// first pushed off that stationary point by a deterministic wobble, capped by
// the same flip test the line search uses, and the derivative is checked there,
// where it is genuinely non-zero.
// ---------------------------------------------------------------------------
double LayoutEnergy::checkGradient(int samples, double h) const {
    if (nV == 0 || samples <= 0) return 0.0;

    // A deterministic, high-frequency wobble: nothing about it matters except
    // that it is smooth in nothing and reproducible.
    std::vector<double> w(2 * nV, 0.0);
    double wmax = 0.0;
    for (int d = 0; d < 2 * nV; ++d) {
        if (d == dofU(options.pinnedVertex) || d == dofV(options.pinnedVertex)) continue;
        w[d] = std::sin(1.234 + 2.399963229728653 * d) + 0.5 * std::cos(0.7 * d);
        wmax = std::max(wmax, std::fabs(w[d]));
    }
    if (wmax <= 0.0) return 0.0;

    const double cap = maxStep(x, w);
    double alpha = 0.02 * imageExtent() / wmax;
    if (std::isfinite(cap)) alpha = std::min(alpha, 0.25 * cap);

    std::vector<double> xt(2 * nV);
    for (int back = 0; back < 40; ++back) {
        for (int d = 0; d < 2 * nV; ++d) xt[d] = x[d] + alpha * w[d];
        if (std::isfinite(energy(xt))) break;
        alpha *= 0.5;
    }

    std::vector<double> grad;
    Eigen::SparseMatrix<double> H = Hproxy;
    if (!model(xt, grad, H)) return kInf;

    double scale = 0.0;
    for (double g : grad) scale = std::max(scale, std::fabs(g));
    if (scale < 1e-30) return 0.0;

    const double step = h * alpha * wmax;
    const int stride = std::max(1, (2 * nV) / samples);
    double worst = 0.0;
    std::vector<double> xp = xt, xm = xt;

    for (int d = 0; d < 2 * nV; d += stride) {
        if (d == dofU(options.pinnedVertex) || d == dofV(options.pinnedVertex)) continue;
        xp[d] = xt[d] + step;
        xm[d] = xt[d] - step;
        const double ep = energy(xp);
        const double em = energy(xm);
        xp[d] = xt[d];
        xm[d] = xt[d];
        if (!std::isfinite(ep) || !std::isfinite(em)) continue;
        const double numeric = (ep - em) / (2.0 * step);
        worst = std::max(worst, std::fabs(numeric - grad[d]) / scale);
    }
    return worst;
}

bool LayoutEnergy::writeOBJ(const std::string &filename) const {
    std::ofstream out(filename);
    if (!out) return false;
    for (const Point &p : uv) out << "v " << p[0] << " " << p[1] << " 0\n";
    for (const Triangle &t : cm->triangles) {
        out << "f " << (t[0] + 1) << " " << (t[1] + 1) << " " << (t[2] + 1) << "\n";
    }
    return true;
}
