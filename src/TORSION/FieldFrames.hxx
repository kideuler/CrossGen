#ifndef __TORSION_FIELDFRAMES_HXX__
#define __TORSION_FIELDFRAMES_HXX__

#include <array>
#include <string>
#include <vector>

#include "MERIDIAN/ConeCut.hxx"
#include "MERIDIAN/ConeSingularities.hxx"
#include "mesh/Mesh.hxx"

class SIPG;

// Pipeline B, work item 2 -- docs/cf_flow_pipeline.md Sec. 5: the matchings of
// the SIPG cross field, the branch they select over Omega, and the per-triangle
// target Jacobian that branch defines.
//
// This is the whole of what the field route has to add before an integration
// can be written, and it is small because the setting is planar. SIPG stores
// u_k[t] = exp(4 i theta_t) with theta in the *global* frame -- initialize()
// reads boundary tangents as atan2(dy, dx) and never introduces a per-face
// reference edge -- so parallel transport is identically zero and there is no
// rho_fg to precompute. For each interior edge e = (f, g),
//
//     theta_f  = arg(u_k[f]) / 4                      in (-pi/4, pi/4]
//     p_fg     = round( (2/pi) (theta_g - theta_f) )  in {-1, 0, +1}
//     delta_fg = (theta_g - theta_f) - (pi/2) p_fg    in [-pi/4, pi/4)
//
// p is the quarter turn the two representatives disagree by and delta is what
// is left, which is the field's own smooth variation. Combing is then a BFS
// over the faces of Omega carrying an integer a_f, with
//
//     a_g = a_f - p_fg ,   theta_hat_f = theta_f + (pi/2) a_f
//
// so that theta_hat_g - theta_hat_f = delta_fg across every dual edge the BFS
// crossed: one branch of the direction field, continuous over Omega.
//
// ### Why the BFS is path-independent, and what it is checked against
//
// Omega is a disk and every cone is on its boundary, so no loop in Omega
// encircles a cone and the accumulated a is independent of the route. That is
// exactly what "cut to a disk" buys and it is what makes C2 hold by
// construction, with no fitting, no rounding and no tolerance -- the transition
// across each arc of G is R_k with k read off the frames either side of it.
//
// It is also the one thing that can silently go wrong, so it is measured rather
// than assumed. Every dual edge of Omega the BFS did *not* use closes a loop,
// and on each of them a_g - a_f must equal -p_fg. A single failure means the
// traversal leaked across G somewhere, which is a statement about the arc
// membership test and not about the field. Report::combingDefects counts them.
//
// ### The index audit (Sec. 5.1)
//
// Around an interior vertex the raw angles are single valued, so the ring sum
// telescopes to zero and sum(delta) + (pi/2) sum(p) = 0. The index the field
// carries there is sum(delta) / (pi/2) -- which is exactly what
// ConeSingularities measures -- so from the matchings alone
//
//     I(v) = - sum over the one ring of p
//
// and comparing the two is a check that the two halves of the same
// decomposition agree. Failing here costs one pass over the vertices; finding
// the same thing as a broken arrangement in Stage 8 costs the pipeline.
//
// The comparison is made against the indices the *field* read, not against the
// set Stage 2 was given: prescribe(), rebalance() and cancelDipoles() all move
// indices deliberately and after the fact, and an audit that flagged their work
// would be reporting the pipeline doing its job. Options::referenceIndex is
// where a caller passes the snapshot; left empty, the cone set as it stands is
// used and the deliberate moves show up as mismatches.
//
// ### The frames
//
//     X_t = (1/h_t) ( cos theta_hat_t,  sin theta_hat_t)
//     Y_t = (1/h_t) (-sin theta_hat_t,  cos theta_hat_t)
//     J*_t = [ X_t | Y_t ]^T
//
// J*_t is the Jacobian the map is asked to have on triangle t: the rows are the
// gradients the integration of Sec. 6 fits, and it is the reference E1 of
// Stage 6 measures against once LayoutEnergy::Reference::Field is switched on.
// h_t is the target edge length. It is one number for the whole model by
// default, and a per-face sizing field is the natural generalisation -- it is
// what sets the density of the eventual layout, and it is Pipeline B's analogue
// of the arbitrary global scaling of Pipeline A's Ricci metric. At h = 1 the
// image comes out about the size of the model, which is the neutral choice and
// the one that makes the field metric's edge lengths comparable with the
// model's own.
//
// det J*_t = 1/h_t^2 is positive by construction for any real theta_hat, so a
// left-handed frame cannot arise from a well-formed comb; the check is kept
// anyway, because a NaN in the field or a sizing field with a zero in it both
// arrive here as det J* <= 0 and both are worth having by name before Stage 6's
// flip cap quietly loses its meaning on those triangles.
class FieldFrames {
public:
    struct Options {
        // h, the target quad edge length in model units. The image of the
        // integration scales by 1/h.
        double targetEdge = 1.0;
        // Per-face h_t. Empty for a uniform targetEdge; a zero or negative
        // entry falls back to targetEdge and is counted.
        std::vector<double> sizing;
        // The face the comb starts from; -1 picks the first face of Omega.
        int seedFace = -1;
        // The indices the cross field itself read, before Stage 1 moved any of
        // them. Empty compares against the cone set as it now stands.
        std::vector<int> referenceIndex;
    };

    struct Report {
        int faces = 0;
        int combedFaces = 0;        // reached by the BFS over the dual graph of Omega
        int unreachedFaces = 0;     // > 0 means Omega's dual graph is disconnected
        int seedFace = -1;

        // Dual edges of Omega that close a loop and disagree with the branch:
        // a_g - a_f != -p_fg. Zero on a correctly cut disk, and the first thing
        // to look at when the indices below do not add up.
        int combingDefects = 0;
        // The worst |delta| over the dual edges of Omega, in radians: how far
        // the field turns between two triangles that share an edge. This is the
        // roughness of the field, and where it is large is where Sec. 7's flips
        // will be.
        double maxFrameJump = 0.0;

        // Sec. 5.1.
        int indexMismatches = 0;        // I from the matchings != the field's own
        int worstMismatchVertex = -1;
        int indexSum = 0;               // sum over interior and boundary cones
        int indexTarget = 0;            // 4 chi(S)
        bool admissible = false;        // Eq. (4)
        int clusteredConePairs = 0;     // two cones sharing one triangle
        int highIndexCones = 0;         // |I| > 2: legitimate, rare, reported

        // The boundary half of Sec. 5.1, which the interior ring sum cannot
        // reach: a boundary vertex's fan is not a closed loop, so its index is
        // not a sum of matchings. What is available instead is the thing the
        // layout will actually have to do there. The field is tangent to dS, so
        // the *image* of a boundary edge under the frame -- J*_f applied to the
        // edge -- is axis aligned, and the quarter turn between the two
        // boundary edges at a vertex is a number the layout has no freedom
        // about: Q3 puts both on coordinate lines and Q2 then says the angle
        // sum is pi - (pi/2) I. So the frame prescribes a boundary index at
        // every vertex of dS, and where that disagrees with the cone set the
        // layout is being asked for two different things at one point.
        //
        // This is worth having by name because it is where a field-integrated
        // layout goes wrong in a way a Ricci one does not. Pipeline A gives
        // Stage 6 a *geodesic* boundary -- Stage 1 prescribes zero curvature at
        // every non-cone boundary vertex and the flow drives it there -- so the
        // image boundary is straight except at the cones and E2 has nothing to
        // argue with. The field route has no such stage: the frame follows a
        // curving boundary and quantises it into a staircase, and every step of
        // that staircase that is not at a cone is a corner the arrangement will
        // find and Q2 will report as a non-cone vertex off pi by a quarter turn.
        int boundaryTurns = 0;          // vertices of dS where the frame turns
        int boundaryTurnMismatches = 0; // ... of those, disagreeing with I(v)
        int worstBoundaryTurnVertex = -1;
        double maxBoundaryTurnResidual = 0.0;  // how ambiguous the worst rounding was

        // Sec. 6.3, measured on the field metric g_t = (J*_t)^T J*_t. With one
        // h for the whole model the metric is (1/h^2) I on every triangle and
        // the two answers for a shared edge agree exactly, so this is zero
        // unless a sizing field is in play -- it is the disagreement a *sizing*
        // field introduces, and not, on its own, a measure of integrability.
        // maxFrameJump above is the number that predicts the flips.
        double maxMetricDisagreement = 0.0;

        int leftHandedFrames = 0;   // det J*_t <= 0
        int sizingFallbacks = 0;    // faces whose h_t was not positive

        bool valid = false;         // combed everything, no defects, no bad frames

        std::vector<std::string> messages;
    };

    // `cut` supplies Omega and, through getCutEdges(), the dual edges the comb
    // must not cross. `cones` is read for the audit only.
    FieldFrames(const SIPG &field, const ConeCut &cut, const ConeSingularities &cones,
                const Options &opts);
    FieldFrames(const SIPG &field, const ConeCut &cut, const ConeSingularities &cones)
        : FieldFrames(field, cut, cones, Options()) {}

    // theta_f = arg(u_k[f]) / 4, one per face of S, in (-pi/4, pi/4].
    const std::vector<double>& rawAngle() const { return rawTheta; }
    // theta_hat_f, the combed branch. One per face of S.
    const std::vector<double>& combedAngle() const { return theta; }
    // The integer a_f the comb assigned, for anyone wanting the branch itself.
    const std::vector<int>& branch() const { return period; }

    // J*_t, row-major, one per face of S.
    const std::vector<std::array<double, 4>>& frames() const { return frame; }
    // h_t as actually used, one per face.
    const std::vector<double>& sizes() const { return h; }

    // p_fg per edge of S, oriented from edgeTriangles[e][0] to [1]. Zero on a
    // boundary edge, which has no second face to match to.
    const std::vector<int>& matchings() const { return matching; }
    // delta_fg on the same indexing, in radians.
    const std::vector<double>& residuals() const { return delta; }

    // I(v) recomputed from the matchings alone, per vertex of S. Interior
    // vertices only; a boundary vertex's fan is not a closed loop and its index
    // is ConeSingularities::computeBoundaryIndices' business, so it is left at
    // the value that class read.
    const std::vector<int>& computedIndex() const { return indexFromMatchings; }

    // The boundary index the frame prescribes at each vertex of dS, and zero
    // everywhere else. See Report::boundaryTurns.
    const std::vector<int>& boundaryTurn() const { return frameBoundaryIndex; }

    // One length per edge of S, from the field metric. Sec. 6.3: g_t is
    // isotropic here, so this is the Euclidean length divided by h, averaged
    // over the two incident triangles.
    const std::vector<double>& fieldEdgeLengths() const { return metricLen; }

    const Report& getReport() const { return report; }

private:
    void readField(const SIPG &field);
    void comb(const ConeCut &cut);
    void buildFrames();
    void audit(const ConeCut &cut, const ConeSingularities &cones);

    const Mesh *mesh = nullptr;
    Options opts;

    std::vector<double> rawTheta;
    std::vector<double> theta;
    std::vector<int> period;
    std::vector<int> matching;
    std::vector<double> delta;
    std::vector<double> h;
    std::vector<std::array<double, 4>> frame;
    std::vector<double> metricLen;
    std::vector<int> indexFromMatchings;
    std::vector<int> frameBoundaryIndex;

    Report report;
};

#endif // __TORSION_FIELDFRAMES_HXX__
