#ifndef __TORSION_FIELDFRAMES_HXX__
#define __TORSION_FIELDFRAMES_HXX__

#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "MERIDIAN/ConeCut.hxx"
#include "MERIDIAN/ConeSingularities.hxx"
#include "MERIDIAN/Interfaces.hxx"
#include "mesh/Mesh.hxx"

class DualMBO;

// Pipeline B, work item 2 -- docs/cf_flow_pipeline.md Sec. 5: the matchings of
// the DualMBO cross field, the branch they select over Omega, and the per-triangle
// target Jacobian that branch defines.
//
// This is the whole of what the field route has to add before an integration
// can be written, and it is small because the setting is planar. DualMBO stores
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
        // The image metric, one length per edge of S, for fieldEdgeLengths() to
        // hand on to Immersion. Empty derives it from the sizing as
        // |e| / h averaged over the two incident faces -- which is what a
        // per-face h can offer, and which is not quite a metric: h varies
        // inside a face's own three edges, so a face can come out failing the
        // triangle inequality and Immersion computes its area by Heron. The
        // cone metric is a vertex scaling and has no such problem, so when one
        // is available it is passed here instead.
        std::vector<double> metricLengths;
        // The face the comb starts from; -1 picks the first face of Omega.
        int seedFace = -1;
        // The indices the cross field itself read, before Stage 1 moved any of
        // them. Empty compares against the cone set as it now stands.
        std::vector<int> referenceIndex;

        // Move the field's own singularities onto the cone set before combing,
        // and re-smooth the branch afterwards. See reconcile() in
        // FieldFrames.cxx: without it, a vertex where Stage 1 moved an index --
        // cancelDipoles() and prescribe() both do, on purpose -- is a
        // singularity of the field that Stage 2 did not cut to, and the branch
        // over Omega is then not single valued at all.
        bool reconcileIndices = true;
        // Transport to and from *sectors* rather than vertices wherever the
        // vertex is on dS or on an interface, route every path clear of both,
        // and hold the pinned faces through the re-smoothing. See
        // reconcileSectors() in FieldFrames.cxx for why the vertex-level
        // transport cannot be trusted there. Off is the vertex-level transport,
        // with dS absorbing anything, kept so the difference stays measurable.
        bool reconcileSectors = true;

        // Read the axis of every chain of dS - G and of every interface branch
        // off the frame, so that Stage 4F can hold it exactly. See
        // alignmentAxis() for what that buys and what it costs.
        bool buildAlignment = true;
        // How far a boundary edge's image may be from the axis the frame says
        // it is on, in radians, before its vote is not counted and the chain is
        // decided by the rest. pi/8 is half of the quarter the vote is between,
        // so an edge whose reading is a coin toss abstains rather than being
        // rounded.
        double alignmentVoteTolerance = M_PI / 8.0;
        // Read an interface branch the cutting graph crosses in the chart each
        // of its pieces is combed in, and vote on it in one. Off votes the
        // branch as though it lay in one chart, which is what this class did
        // before and what folds the far side of every crossing through a right
        // angle; kept so the difference stays a measurement. See
        // buildAlignment().
        bool alignAcrossSeams = true;
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

        // reconcile(). The field's singularities against Stage 1's cone set,
        // before anything was done about it, and what it took to agree.
        int reconciledVertices = 0;     // interior vertices where the two differed
        int reconciledUnits = 0;        // index units moved
        int reconciledEdges = 0;        // matchings changed to move them
        int reconcileFailures = 0;      // units with nowhere to go
        // reconcileSectors(): sectors of dS and of the interface network where
        // the frame's quarter count and the layout's differed, how many of
        // those were left because no route could reach a sector wanting the
        // opposite, and how many interior units had to be absorbed by a sector
        // that did not want them (the vertex-level transport's only option).
        int reconciledSectors = 0;
        int unresolvedSectorUnits = 0;
        int absorbedUnits = 0;
        bool branchResmoothed = false;
        double maxBranchCorrection = 0.0;   // the largest |c_f|, radians

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

        // The alignment, Sec. 6.4. One chain of dS - G, or one branch of the
        // interface network, is one curve the layout has to put on one
        // coordinate line, and these count what came out of reading which line
        // that is off the frame.
        int alignmentChains = 0;        // chains and branches given an axis
        int alignmentSteppedChains = 0; // ... of those, with a step overridden
        int alignmentReversedChains = 0;// left free: the field turns them back
        int alignmentAbstained = 0;     // left free: no edge near an axis at all
        int alignedBoundaryEdges = 0;   // edges of Omega on dS that were held
        int alignedInterfaceEdges = 0;  // ... and on the interface network
        // Edges whose own frame reading disagreed with the axis its chain was
        // given. This is the staircase of Report::boundaryTurnMismatches seen
        // per edge rather than per vertex, and it is the number the alignment
        // exists to remove: each of these is a step the field put in the middle
        // of a curve that Stage 1 says runs straight, and holding the chain to
        // one axis is what takes it out.
        int alignmentOverrides = 0;
        // Places where an interface branch crosses the cutting graph and its
        // reading changes chart by a non-zero quarter turn, each taken back
        // before the vote rather than counted as a step of the staircase.
        int alignmentSeamCrossings = 0;
        // Chains left unaligned because they close on themselves with no cone
        // on them. Holding a closed curve to one coordinate line collapses it,
        // so they are reported and left to E2.
        int alignmentClosedChains = 0;
        // The worst |image direction - the nearest quarter turn| over the
        // aligned edges, in radians. Zero to rounding is the field exactly
        // tangent to the curve; a large value is a curve the field is not on,
        // and holding it there will cost the integration.
        double maxAlignmentResidual = 0.0;

        bool valid = false;         // combed everything, no defects, no bad frames

        std::vector<std::string> messages;
    };

    // `cut` supplies Omega and, through getCutEdges(), the dual edges the comb
    // must not cross. `cones` is read for the audit only.
    //
    // `interfaces` is optional and read only by the alignment: an interface
    // branch is a curve of the layout in exactly the way a chain of dS is, and
    // the field is pinned tangent to it in exactly the same way, so it takes
    // the same treatment. Null, or single-material, aligns dS alone.
    FieldFrames(const DualMBO &field, const ConeCut &cut, const ConeSingularities &cones,
                const Options &opts, const Interfaces *interfaces = nullptr);
    FieldFrames(const DualMBO &field, const ConeCut &cut, const ConeSingularities &cones)
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

    // Sec. 6.4 -- which coordinate each edge of *Omega* is to hold constant:
    //
    //     0   hold u   (the image of the edge is vertical)
    //     1   hold v   (the image of the edge is horizontal)
    //    -1   free
    //
    // Q3 says every component of dS - G lies on a constant-coordinate line and
    // E3 says the same of every feature; in Pipeline A both arrive by the time
    // Stage 6 has run, and on a Ricci map the first of them arrives a stage
    // early because the flow makes dS geodesic. The field route has no such
    // stage, and a least-squares map is only as aligned as the field is
    // integrable -- which is what leaves Stage 6 a boundary residual of 1e-1 to
    // work off where Pipeline A hands it 1e-15, and it is the whole of why the
    // continuation runs to convergence here and exits at the first level there.
    //
    // The remedy costs one linear equation per edge, because the statement is
    // linear: u_i = u_j, or v_i = v_j. Which one is not free either -- the
    // field is pinned tangent to dS and to every interface, so J*_t applied to
    // such an edge is axis aligned already and the axis can simply be read.
    // What is read per chain rather than per edge is the point: a field
    // following a curving boundary quantises it into a staircase, and every
    // step of that staircase inside a chain that Stage 1 says runs straight is
    // a corner the arrangement would find and Q2 would report. One axis for the
    // whole chain, by the vote of its edges weighted by their length, is that
    // staircase removed -- and the chains end at the cones, which is where the
    // layout is supposed to turn.
    const std::vector<int>& alignmentAxis() const { return alignAxis; }

    // The same, restricted to the chains whose every edge already agreed with
    // the axis it was given -- no staircase step to override. A step is a right
    // angle the map is being asked to unbend and is where holding the boundary
    // costs the integration its injectivity, so this is the set to fall back to
    // before giving up on the alignment altogether. See FieldFrames.cxx.
    const std::vector<int>& strictAlignmentAxis() const { return strictAxis; }

    // Which chain each held edge of Omega belongs to, numbered in the order
    // they were given an axis, and -1 for an edge held by none. A chain is the
    // unit the alignment is decided in, so it is also the unit it is let go
    // in: Stage 4F releases the chains a tangle touches and keeps the rest.
    const std::vector<int>& alignmentChain() const { return alignChain; }

    // One length per edge of S, from the field metric. Sec. 6.3: g_t is
    // isotropic here, so this is the Euclidean length divided by h, averaged
    // over the two incident triangles.
    const std::vector<double>& fieldEdgeLengths() const { return metricLen; }

    const Report& getReport() const { return report; }

private:
    void readField(const DualMBO &field);
    void comb(const ConeCut &cut);
    void buildFrames();
    // The index each interior vertex's one ring of matchings carries, and the
    // sign with which one edge enters that ring sum.
    void indexRing(std::vector<int> &out) const;
    int ringSign(int v, int e) const;
    void reconcile(const ConeCut &cut, const ConeSingularities &cones);
    void reconcileSectors(const ConeCut &cut, const ConeSingularities &cones,
                          const Interfaces *interfaces);
    void smoothBranch(const ConeCut &cut);
    void audit(const ConeCut &cut, const ConeSingularities &cones);
    void buildAlignment(const ConeCut &cut, const ConeSingularities &cones,
                        const Interfaces *interfaces);

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
    std::vector<char> pinned;      // per face: the field's Dirichlet data
    std::vector<int> alignAxis;    // per edge of Omega
    std::vector<int> strictAxis;   // ... the chains with nothing overridden
    std::vector<int> alignChain;   // per edge of Omega: its chain, or -1

    Report report;
};

#endif // __TORSION_FIELDFRAMES_HXX__
