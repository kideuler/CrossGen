#ifndef __TORSION_CONEMETRIC_HXX__
#define __TORSION_CONEMETRIC_HXX__

#include <string>
#include <vector>

#include "mesh/Mesh.hxx"

class ConeSingularities;

// The flat cone metric of Stage 1's cone set, as a conformal change of the
// model's own metric -- docs/cf_flow_pipeline.md Sec. 4, which is C4.
//
// ### What this is for
//
// E1 of Stage 6 measures J against a reference metric, and Sec. 2.1 of the plan
// records what happens when that reference has no cones: E1 and Q2 become
// contradictory statements about the same vertex, the cone fans come out
// uneven, and the patch at the cone acquires a reflex corner that Stage 10
// meshes as an inverted element. The plan's own answer was Reference::Field --
// compose the Euclidean reference with the field's frame -- and that answer is
// wrong, measurably: the frame is (1/h) R(-theta), E1 is invariant under
// right-multiplication by a rotation, and the reference is therefore the
// Euclidean one exactly. --ref-test reports the two as identical because they
// are.
//
// The cone structure lives in the *scale*, not in the rotation. What carries it
// is the conformal factor of Sec. 3:
//
//     u + i theta  =  - sum_k (I_k / 4) log(z - p_k)  +  holomorphic
//
// u and theta are conjugate parts of one holomorphic function, so a route that
// has combed theta has all the information u needs; and on a *planar* model u
// is available directly and much more cheaply than that, because the interior
// curvature it starts from is already zero.
//
// ### The construction
//
// A conformal change of a triangulation is a vertex scaling (Springborn,
// Schroeder and Pinkall): one u_i per vertex, and
//
//     l_ij(u) = exp( (u_i + u_j) / 2 ) l_ij(0).
//
// Differentiating the angle at i in a triangle ijk with respect to u_j gives
// cot(alpha_k) / 2 -- so the derivative of the angle *defect*
//
//     K_i = (2 pi or pi) - sum_t alpha_i^t
//
// with respect to u is exactly the cotangent Laplacian of the current metric,
// dK = L du, with the usual w_ij = (1/2) sum_t cot(alpha_ij^t). Newton is then
//
//     L du = Kbar - K,     u <- u + du
//
// with Kbar_i = (pi/2) I_i the curvature Stage 1 prescribes. The first step is
// taken at u = 0, where the model is flat in the interior and L is the ordinary
// cotangent Laplacian of the input: that single linear solve is the whole of
// the answer to first order, and the iteration is what takes it to rounding.
//
// L is singular on the constants, which are a global scaling of the metric and
// change nothing; one vertex is pinned to remove them. Solvability needs
// sum(Kbar) = sum(K) = 2 pi chi, which is Eq. (4) -- checked at Stage 1 and
// checked again here rather than assumed.
//
// ### Why not the flow
//
// RicciFlow answers the same question and is used by Pipeline A. It carries a
// circle-packing metric, its own flow triangulation and an edge-flip pass,
// because on a curved input the metric can leave the realisable set. Here the
// input is planar and the curvature being moved is only what sits on the
// boundary, so the same answer comes out of a handful of Laplacian solves on
// the model's own triangulation. Realisability is still not free -- a strong
// enough cone can violate the triangle inequality -- so every step is halved
// until it holds, and Report::nonRealisable says if that ran out.
class ConeMetric {
public:
    struct Options {
        int newtonSteps = 40;
        // ||K - Kbar||_inf, radians. Q2's own tolerance is 1e-6 and C4's is
        // 1e-3, so this is two orders below the tighter of them.
        double tolerance = 1e-9;
        // The smallest fraction of a Newton step that is still tried before the
        // step is abandoned to the triangle inequality.
        double minStepFraction = 1e-3;
        // How far below the triangle inequality a face may sit and still count
        // as realisable, relative to the longest side.
        double realisableMargin = 1e-9;
    };

    struct Report {
        int newtonIterations = 0;
        bool converged = false;
        bool solved = false;
        // ||K - Kbar||_inf before any step -- the curvature the model carries
        // on its boundary and the layout wants in its cones -- and after.
        double initialError = 0.0;
        double linearError = 0.0;   // after the first (linear) solve alone
        double finalError = 0.0;
        // The residual measured the way Stage 6 measures C4: the worst cone
        // whose angle sum in this metric is not 2 pi - (pi/2) I.
        double coneResidual = 0.0;
        double gaussBonnetResidual = 0.0;
        double minScale = 1.0, maxScale = 1.0;    // exp(u) over the vertices
        int nonRealisable = 0;      // faces failing the triangle inequality
        int halvings = 0;           // Newton steps that had to be shortened
        std::vector<std::string> messages;
    };

    ConeMetric(const Mesh &mesh, const ConeSingularities &cones, const Options &opts);
    ConeMetric(const Mesh &mesh, const ConeSingularities &cones)
        : ConeMetric(mesh, cones, Options()) {}

    // The conformal factor, one per vertex of the model.
    const std::vector<double>& conformalFactor() const { return u; }

    // The flat cone metric itself: one length per edge of the model, in the
    // indexing Mesh::edges has them. This is what E1's reference is built from
    // and what Immersion is handed as its flat lengths -- the same role
    // RicciFlow::originalEdgeLengthsCompleted() plays in Pipeline A.
    const std::vector<double>& edgeLengths() const { return len; }

    // exp(u) averaged over a face, which is the factor the frame's target
    // Jacobian is to be scaled by: J*_t = (exp(u_t) / h) R(-theta_t). One over
    // it is the per-face h_t FieldFrames::Options::sizing takes.
    const std::vector<double>& faceScale() const { return fscale; }

    // h_0 / faceScale(), ready for FieldFrames::Options::sizing.
    std::vector<double> sizingField(double targetEdge) const;

    const Report& getReport() const { return report; }

private:
    void solve(const Mesh &mesh, const ConeSingularities &cones);
    void rebuild(const Mesh &mesh);   // len and fscale from u

    std::vector<double> u;
    std::vector<double> len;
    std::vector<double> fscale;
    std::vector<double> baseLen;
    Options opts;
    Report report;
};

#endif
