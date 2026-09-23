#ifndef __UMBER_HXX__
#define __UMBER_HXX__

#include <memory>
#include <unordered_set>
#include <vector>

// eigen includes
#include <Eigen/Dense>

#include "Parameterization/HarmonicCut.hxx"
#include "mesh/Mesh.hxx"

class DualMBO;

// Frame field generation of
//   Wang, Ren, Fang, Lin, Xu, Bao and Huang, "IGA-suitable planar
//   parameterization with patch structure simplification of closed-form
//   polysquare", CMAME 392 (2022) 114678, Sections 4.1 and 4.2.
//   (refs/IGA-suitable-planar-parameterization...pdf)
//
// The polysquare pipeline needs a frame field that is smooth, boundary
// aligned, and -- unlike the field a quadrangulation method would produce --
// free of *internal* singularities, because a polysquare structure has none.
// Methods that carry the 4-symmetry in their unknowns (MIQ, MBO/DualMBO,
// polyvectors) cannot express that requirement: the symmetry is exactly what
// lets a singularity form at no cost. Sec. 4.2 therefore drops the symmetry
// and optimizes a *non-symmetric* field, one plain 2D vector v per triangle
// with its 90-degree rotation v^perp completing the frame. A field that is
// smooth in v has no interior holonomy at all, so smoothing it directly is
// what pushes the singularities out to the boundary.
//
// The symmetry cannot be dropped everywhere, though: a domain with beta holes
// admits beta non-trivial harmonic 1-forms, and it is precisely those extra
// degrees of freedom that promote a common polysquare (an exact form) into a
// closed-form polysquare, the structure that puts no corners on a circular
// ring at all. Sec. 4.1 exposes them as transitions across a set of cuts C
// that break the holes open. Across those edges -- and only there -- the
// smoothness term is measured in the symmetric representation tau(v) of
// Eq. (3), which leaves the k*90-degree transition free; measuring them
// non-symmetrically would pin every transition to the identity and collapse
// the result back to a common polysquare.
//
// So Eq. (1) is minimized over the per-triangle vectors v:
//
//     min_v  E_smooth(v) + w_a E_align(v) + w_r E_reg(v),
//
//   E_smooth  Eq. (2): squared jump of v across every interior edge, in the
//             non-symmetric representation off the cuts and in tau(v) on them,
//             area weighted.
//   E_align   Eq. (4): l1 norm of each boundary triangle's vector rotated into
//             its edge frame, length weighted. The l1 norm (rather than a
//             squared one) is what makes the boundary snap to an axis instead
//             of averaging across a corner; see setSmoothingSchedule() for how
//             it is made differentiable.
//   E_reg     Eq. (5): (||v||^2 - 1)^2, the only thing keeping v off zero,
//             since the frame is encoded by v alone and not by a constrained
//             rotation.
//
// Encoding the 2D frame as a vector rather than as the 2x2 rotation matrix of
// the 3D method it adapts (Fang et al. 2016) is the paper's own speedup: half
// the unknowns and one normalization penalty instead of four orthogonality
// constraints, worth a factor of three to four (Sec. 7.1, Fig. 20). Eq. (25)
// shows the two energies agree up to a factor of 2.
//
// The cuts themselves come from HarmonicCut, which implements Sec. 4.1: one
// short, disjoint cut per void, and nothing else. That is what C has to be for
// the free transitions here to correspond one-to-one with the harmonic forms.
class UMBER {
public:
    using EdgeKey = MeshEdgeKey;
    using EdgeKeyHash = MeshEdgeKeyHash;

    // The three terms of Eq. (1), each already normalized as the paper defines
    // it (by total area or total boundary length) but *not* yet weighted;
    // `total` is E_smooth + w_a E_align + w_r E_reg.
    struct EnergyTerms {
        double smooth = 0.0;
        double align = 0.0;
        double reg = 0.0;
        double total = 0.0;
    };

    // Take the initial field from a converged DualMBO cross field and the cuts
    // from a HarmonicCut built on the same mesh. C is exactly the beta cuts of
    // Sec. 4.1 -- one per void, disjoint, meeting the boundary only at their
    // endpoints -- so the free transitions in E_smooth match the harmonic
    // degrees of freedom one for one.
    //
    // If the cutting failed to open every void (check
    // harmonicCut.getReport()), the field is still optimized, but the voids
    // that stayed closed contribute no harmonic form: their frames are pinned
    // to an exact form, which is the common polysquare behaviour with its
    // forced corners.
    UMBER(const DualMBO &dualMBO, const HarmonicCut &harmonicCut);

    // Direct form: per-triangle initial directions (any representative of the
    // cross is fine, the comb in initialize() fixes the branch) and the cut
    // edge set, keyed by *original* mesh vertex indices.
    UMBER(std::shared_ptr<Mesh> mesh,
          const std::vector<Point> &initialField,
          const std::unordered_set<EdgeKey, EdgeKeyHash> &cuts);

    // Interior edges the frame is to be aligned to on top of dS: the material
    // interfaces of a multi-material domain. Mesh edge indices; call before
    // initialize(), and an empty set is the single-material behaviour.
    //
    // Sec. 4.2 knows only one kind of feature, dS, because the paper's models
    // have only one. An interface is the same kind of thing and for the same
    // reason: it is a curve the output has to keep, since an element
    // straddling it carries two materials and no analysis code can integrate
    // that (MERIDIAN::Interfaces says this at length for Pipeline A). A frame
    // that is not tangent to it cannot produce a polysquare whose blocks
    // follow it, so the alignment has to see it here or nothing downstream
    // can.
    //
    // An interface edge is *two* terms of Eq. (4) and not one: it has a
    // triangle on each side and both are asked for the same axis, where a
    // boundary edge asks its one triangle. That is the same treatment
    // DualMBO::setAlignedInteriorEdges gives the field this one starts from,
    // and it has to be, or the two sides of the interface agree about the
    // curve and disagree about the frame on it.
    //
    // The normalization of Eq. (4) is over the features' total length, so
    // adding interfaces dilutes the weight each boundary edge carries. That is
    // deliberate: E_align is one term about one set of curves the field must
    // follow, and dS is not privileged among them.
    void setFeatureEdges(const std::vector<int> &edges);

    // Whether any were given: the one thing downstream has to branch on.
    bool hasFeatureEdges() const { return !featureEdges.empty(); }
    const std::vector<int>& getFeatureEdges() const { return featureEdges; }

    // w_a of Eq. (1). The paper's 0.1, found experimentally (Sec. 7.1,
    // Fig. 14): too small and the field is smooth but ignores the boundary,
    // too large and the smoothness term is overpowered and the field turns
    // ragged.
    void setAlignmentWeight(double wa) { alignWeight = wa; }
    double getAlignmentWeight() const { return alignWeight; }

    // w_r of Eq. (1). The paper's 5, and it reports insensitivity over
    // [0.5, 500] (Sec. 7.1). Below that the vectors near the boundary take on
    // non-uniform lengths, which makes the 2-norm of Eq. (2) a poor measure of
    // smoothness and lets internal singularities back in; far above it the
    // strongly non-linear energy tends to stall in a local minimum.
    //
    // The low end of that range does not survive here, because E_align is
    // linear in ||v|| and so pays a boundary triangle to shrink: at w_r = 0.5
    // on data/meshes/singlemat/geom003.obj the boundary vectors collapse and E_align
    // reads 2e-4 instead of the ~1.0 a unit aligned field gives, which is a
    // degenerate field that merely scores well. 5, 50 and 500 all hold
    // E_align at 0.89, 0.99, 1.00 with no change in the singularity count, so
    // the paper's default sits at the bottom of the usable band rather than in
    // the middle of it.
    void setRegularizationWeight(double wr) { regWeight = wr; }
    double getRegularizationWeight() const { return regWeight; }

    // The l1 norm of Eq. (4) is not differentiable at 0, which is where the
    // minimizer wants to sit, so |t| is replaced by sqrt(t^2 + eps^2) and eps
    // is walked down, each stage warm-starting from the last. A large eps
    // rounds the corner enough for the quasi-Newton steps to find the basin; a
    // small one is what actually snaps the boundary onto an axis. Passing a
    // single value disables the continuation.
    void setSmoothingSchedule(const std::vector<double> &schedule) { l1Schedule = schedule; }

    // Let a +1/-1 pair of input singularities annihilate each other instead of
    // both leaving through dS.
    //
    // The branch cuts of the comb decide where the defects go, and by default
    // each runs to the nearest boundary, so every quarter singularity of the
    // input becomes a corner of the polysquare. For a pair of *opposite* sign
    // that is the wrong destination and, on a multi-material model, actively
    // harmful: a cross field aligned to a curved interface puts a +1/-1 pair
    // on the concave side of every strongly curved stretch of it, the pair
    // contributes nothing to any Gauss-Bonnet count, and pushing both members
    // out invents two corners the shape does not have -- which E_cor of the
    // deformation then insists on realising. This is the same fact
    // MERIDIAN::ConeSingularities::cancelDipoles is written against, and the
    // same rule applies: a pair in one material region may cancel, a pair that
    // *spans* an interface never may, because each of those two is the whole
    // of its own region's deficit.
    //
    // A cut between the two is topologically sound where a cut between a
    // same-sign pair is not: a loop around a +1/-1 pair has no holonomy and
    // crosses the cut twice in opposite directions, so the two jumps cancel,
    // which is exactly right. A loop around a +1/+1 pair has 180 degrees and
    // needs one crossing per singularity -- see singularitySeams() -- so
    // same-sign pairs are never joined.
    //
    // On by default. Off restores the behaviour where every input singularity
    // is routed to dS.
    void setCancelDipoles(bool on) { cancelDipoles = on; }
    bool getCancelDipoles() const { return cancelDipoles; }

    // Pairs annihilated by the above at the last initialize(), as the two
    // vertices of each. Empty when there were none, which is every
    // single-material model in data/meshes whose field carries no such pair.
    const std::vector<std::pair<int, int>>& cancelledDipoles() const { return dipoles; }

    // Stopping criteria per continuation stage. The paper stops L-BFGS at a
    // gradient magnitude of 1e-6 (Sec. 7.1).
    void setTolerance(double tol) { gradTolerance = tol; }
    void setMaxIterations(int iters) { maxIterations = iters; }

    // Classifies the edges, computes the areas and lengths of Eqs. (2), (4),
    // (5), and builds the initial v: the input cross directions, combed by
    // breadth-first traversal of the dual graph with the cut edges removed, so
    // that neighbouring vectors agree on which of the four cross directions
    // they represent. This is the "global alignment on the cut one M_C"
    // followed by "consistently picking an axis" of Sec. 4.2.
    //
    // If the input field has a quarter singularity, no comb can be consistent
    // everywhere: some edges are left carrying a 90-degree jump, and they line
    // up along a branch cut running from the singularity to the boundary. That
    // is expected -- the jump is a large but local penalty that E_smooth
    // removes by pushing the singularity out along the cut -- but *where* the
    // cut runs is not a detail, because it decides where the singularity
    // leaves and therefore where the optimized field puts its boundary corner.
    // So the cuts are routed on shortest paths to the boundary rather than
    // left to the traversal; see singularitySeams() in the implementation for
    // what happens when they are not.
    //
    // The paper's stated failure mode still applies -- a singularity deep
    // inside the model may not make it out (Sec. 8) -- so check
    // internalSingularities() on the result.
    void initialize();

    // Minimizes Eq. (1) with L-BFGS, over the l1 smoothing schedule. Calls
    // initialize() first if that has not happened yet. The resulting vectors
    // are normalized to unit length before being published, since E_reg only
    // penalizes the deviation rather than constraining it.
    void optimize();

    // The optimized frame: unit v per triangle, and its 90-degree rotation
    // v^perp. Both are in the *non-symmetric* representation, so neighbouring
    // triangles off the cuts carry directions that genuinely match, not
    // matching-up-to-90-degrees.
    const std::vector<Point>& getUField() const { return uField; }
    const std::vector<Point>& getVField() const { return vField; }

    // The same field as exp(4 i theta) per triangle, i.e. in the convention of
    // DualMBO::u_k, for the parts of the pipeline that consume a cross field
    // symmetrically (CutMesh, tracing). Reducing to it discards the branch the
    // optimization just fixed, so prefer getUField() where the branch matters.
    Eigen::VectorXcd crossFieldRepresentation() const;

    // Per-triangle direction angle in (-pi, pi], i.e. atan2 of getUField().
    const std::vector<double>& getAngles() const { return angles; }

    // Internal singularities of the *current* field: interior vertices where
    // the holonomy of v does not close, as (vertex index, winding). This is
    // the direct check on what Sec. 4.2 is for, and it is measured on v rather
    // than on tau(v) -- see the note in the implementation for why that
    // distinction decides whether the answer is right or off by an order of
    // magnitude. The free k*90-degree transitions on the cuts are undone
    // first; they carry the harmonic part of the field, not a defect.
    //
    // Empty is the goal but not a guarantee. Measured on data/meshes, DualMBO
    // field in, Sec. 4.1 cuts, defaults throughout: all sixteen models come
    // out with none, E_smooth falls by roughly an order of magnitude
    // everywhere, and no interior edge anywhere turns past 45 degrees (see
    // sharpTurns()).
    //
    // Getting there depends on the branch cuts of the comb, not on the solver.
    // With the cuts left where a plain breadth-first traversal happened to
    // close, one of the sixteen kept a defect, eleven ended at a higher total
    // energy (up to 14% higher), and several put their boundary corners in
    // places the geometry does not support. Routing each cut to the nearest
    // boundary instead cleared the defect and lowered the energy without
    // touching Eq. (1); singularitySeams() in the implementation says why.
    //
    // The paper's own failure mode still stands (Sec. 8: a singularity deep
    // inside the model may not come out) -- near the boundary the cut a defect
    // drags behind it is short and shrinking it walks the defect out, deep
    // inside there is less of a gradient to follow -- so this is a measurement
    // on sixteen models, not a guarantee.
    //
    // The cut set of Eq. (2) is a separate question from the branch cuts and
    // does not decide this. The argument for the Sec. 4.1 cuts is structural
    // -- beta transitions, one per void, as the closed form requires -- not
    // that they clean up more singularities.
    //
    // The ring of Fig. 10, the case the closed-form structure exists for, does
    // come out clean: one cut, and the optimized frame follows the annulus the
    // whole way round -- on a 23.5k-triangle ring of radii 0.4 and 1, 0.30
    // degrees from circumferential on average and 2.07 at worst, with no
    // internal singularity and none of the four corners an exact-form
    // polysquare would be forced to place.
    std::vector<std::pair<int, double>> internalSingularities() const;

    // The other side of the same coin: the corners the field puts *on* the
    // boundary, as (vertex index, quarter turns). Pushing the defects out here
    // is what Sec. 4.2 is for, and these are exactly the corners of the
    // polysquare the field induces -- +1 at a convex one, -1 at a reflex one,
    // and nothing at all along a run of boundary the frame follows.
    //
    // Measured per boundary vertex as the mismatch between how far the frame
    // turns across the vertex star and how far the boundary itself turns:
    //
    //     k = round( (Theta + pi - Omega) / (pi/2) ),
    //
    // with Omega the interior angle at the vertex and Theta the accumulated
    // rotation of v over the star, walked counter-clockwise from one boundary
    // edge to the other. A straight piece of boundary followed by the frame
    // has Omega = pi and Theta = 0 and so reads 0; a square corner has
    // Omega = pi/2 and gives 1. E_align is what makes the rounding meaningful:
    // it holds the frame on an axis of the boundary edge at both ends of the
    // star, so the quantity being rounded really does sit near a multiple of
    // 90 degrees.
    //
    // The k*90-degree transitions on the cuts are undone first, as in
    // internalSingularities(), and a vertex where the boundary pinches (more
    // than two boundary edges) is skipped, having no single corner to report.
    //
    // A corner is spread over the two or three vertices it takes the frame to
    // swing across, so the quarter turns are not rounded vertex by vertex --
    // that reports a corner split 0.35/0.45/0.20 as no corner at all. They are
    // accumulated along each boundary loop and handed to the vertex that
    // carries the running total past the halfway mark, which lands on the
    // middle of the corner and makes the count sum to 4*chi by construction.
    // Read the reported vertex as the centre of a corner a vertex or two wide,
    // not as an exact location.
    //
    // On a field with no internal singularity the k sum to 4*chi, i.e. 4 on a
    // disk -- the four corners a common polysquare must have. On a holed model
    // the free transitions carry part of the total and the sum need not match.
    std::vector<std::pair<int, int>> boundarySingularities() const;

    // The most quarter turns boundarySingularities() will put on one vertex.
    //
    // 1, and it is a structural limit rather than a tuning knob. A corner of a
    // polysquare turns through one quarter; a vertex carrying two turns the
    // image boundary through 180 degrees, which is an interior angle of zero --
    // a needle of no width. Nothing downstream survives one. The face at a
    // needle's tip has three sides, not four, because 180 + 90 + 90 already
    // closes it, so the tip costs the block decomposition that face whatever
    // the tracing does; on data/meshes/singlemat/geom028, a six-pointed star
    // whose 40.5-degree tips each read 1.55 quarters and round to 2, it cost
    // five faces and 20.9% of the model.
    //
    // The cap does not round the surplus away -- that would break the
    // Gauss-Bonnet count the whole stage rests on -- it carries it to the next
    // vertex along, which spends it as a second ordinary corner. A sharp tip
    // comes out as a blunt end one edge wide, which is what a polysquare's
    // version of a sharp tip is.
    //
    // Raising it restores the old behaviour, and is here to measure against.
    void setMaxCornerQuarters(int q) { maxCornerQuarters = q; }
    int getMaxCornerQuarters() const { return maxCornerQuarters; }

    // ── The same reading, taken along the material interfaces ───────────────
    //
    // boundarySingularities() answers "where does the frame put a corner of
    // the polysquare boundary". On a multi-material domain the interfaces ask
    // the same question of themselves: an interface is a curve the layout has
    // to follow, so where the frame turns against it is a corner of the
    // structure just as surely as a corner of dS is, and Eq. (12) needs the
    // number there for the same reason it needs it on dS.
    //
    // The measurement is the same one. At a vertex, the two feature edges
    // meeting there cut the star into the sector on the left of the walk
    // edgeIn -> v -> edgeOut; over that sector the frame turns by Theta and
    // the model subtends Omega, and
    //
    //     quarters = round( (Theta + pi - Omega) / (pi/2) )
    //
    // is how far the image curve turns there, in quarter turns, positive to
    // the left. A straight run followed by the frame reads 0; a square corner
    // of the sector reads +1. On dS, walked with the material on its left, it
    // is exactly what boundarySingularities() reports.
    //
    // The difference is that an interface has two sides, so it appears twice
    // -- once per direction of the walk -- with the sector measured always on
    // the left. A caller walking a branch with material m on its left takes
    // the entry whose `edgeIn`/`edgeOut` run that way.
    //
    // As on dS, the quarter turns are accumulated along each chain rather than
    // rounded vertex by vertex: a corner is spread over the two or three
    // vertices it takes the frame to swing across, and rounding each on its
    // own reports a corner split 0.35/0.45/0.20 as no corner at all.
    //
    // Empty unless setFeatureEdges() was given something.
    struct FeatureCorner {
        int vertex = -1;
        int edgeIn = -1;    // the feature edge arriving at `vertex`
        int edgeOut = -1;   // the one leaving it
        int quarters = 0;
    };
    std::vector<FeatureCorner> interfaceCorners() const;

    // The interface network cut into chains: each a run of interface edges
    // between two places the network is not simply two edges carrying on --
    // an end, a junction of three or more branches, or a landing on dS, none
    // of which the network itself says how to continue through.
    //
    // Exposed because interfaceCorners() measures along these, and a stage
    // consuming those corners has to walk the same chains in the same order
    // or the corner it looks up is not the corner it is at. Built once, on
    // first use, and empty without setFeatureEdges().
    struct FeatureChain {
        std::vector<int> edges;   // mesh edges, in order along the chain
        std::vector<int> verts;   // edges.size() + 1 of them, in the same order
        bool closed = false;      // verts.front() == verts.back()
    };
    const std::vector<FeatureChain>& interfaceChains() const;

    // Energy at the current field, split into the terms of Eq. (1). Valid
    // after initialize() (the initial field) and after optimize() (the
    // result), so comparing the two shows what the optimization bought.
    EnergyTerms energy() const { return currentEnergy; }

    // 2-norm of the gradient at the last accepted iterate, and the total
    // number of L-BFGS iterations over all continuation stages.
    double gradientNorm() const { return currentGradNorm; }
    int iterations() const { return totalIterations; }

    // Holes whose rotation number the initialization had to unwind -- see
    // unwindHoleHolonomy(). Nonzero means the smoothest field on this model
    // was a cornerless ring and the one with a polysquare had to be asked for.
    int unwoundHoles() const { return holesUnwound; }

    // How many pieces the branch cuts left the comb's dual graph in; 1 is the
    // healthy case. See combInitialField().
    int combPieces() const { return combComponents; }

    // Branch cuts that had to be sent past the hole they were nearest; see
    // singularitySeams().
    int reroutedSeams() const { return seamsRerouted; }

    // Interior edges, cuts excluded, where the frame turns by more than
    // `thresholdDegrees` between the two triangles, as (edge index, degrees),
    // worst first. Not singularities -- the field is single-valued and smooth
    // across them -- but worth knowing about for two reasons.
    //
    // First, they are where the parameterization built on this field will lose
    // its scaled Jacobian: a frame that swings 50 degrees in one triangle step
    // is a distortion hot spot, and usually says the mesh is too coarse to
    // resolve a boundary feature the alignment term is trying to follow.
    //
    // Second, they are landmines for anything downstream that measures the
    // field 4-symmetrically -- DualMBO::computeSingularities, cross-field tracing
    // -- because past 45 degrees the spin-4 winding cannot tell which way the
    // frame turned and will report a singularity that is not there. If such a
    // tool disagrees with internalSingularities(), compare its count against
    // this one first.
    std::vector<std::pair<int, double>> sharpTurns(double thresholdDegrees = 45.0) const;

    // Shortest ||v|| over the triangles, before the output was normalized.
    // E_reg holds it near 1 -- 0.79 to 0.97 across data/meshes at the default
    // weights, the shrinkage sitting on the boundary triangles where E_align,
    // which is linear in ||v||, pulls against it. A much smaller value means
    // w_r has lost, and then the 2-norm of Eq. (2) is measuring length
    // differences rather than smoothness; the paper warns that this is where
    // internal singularities come back (Sec. 7.1). It is not a singularity
    // detector in its own right: the models that keep one show no collapse,
    // because a discrete vortex can turn a full circle over a small vertex
    // star without any vector having to shorten.
    double minFrameNorm() const { return minFrameNorm_; }

    // The cut edges actually used as C in Eq. (2), in original-mesh vertex
    // indices.
    const std::unordered_set<EdgeKey, EdgeKeyHash>& getCutEdges() const { return cutEdges; }

    const Mesh& getMesh() const { return *mesh; }
    std::shared_ptr<Mesh> getMeshPtr() const { return mesh; }

private:
    // One interior edge of Eq. (2). `weight` is A_ab / A_M; `onCut` selects
    // the tau(v) form of the jump.
    struct SmoothEdge {
        int fa = -1;
        int fb = -1;
        double weight = 0.0;
        bool onCut = false;
    };

    // One boundary edge of Eq. (4). The rotation R(-sigma_e) into the edge
    // frame is stored by its cosine and sine; `weight` is l_e / l_dM.
    struct AlignEdge {
        int face = -1;
        double cosSigma = 1.0;
        double sinSigma = 0.0;
        double weight = 0.0;
    };

    // Eq. (1) and its gradient at x (the 2*|T| stacked components of v).
    // `terms` optionally receives the unweighted split.
    double evaluate(const Eigen::VectorXd &x, Eigen::VectorXd &grad,
                    EnergyTerms *terms = nullptr) const;

    // L-BFGS with backtracking Armijo line search on evaluate(), at the
    // current l1 smoothing. Returns the iterations taken; leaves the iterate
    // in `x`. A quasi-Newton method is what the paper uses (L-BFGS from
    // ALGLIB, Sec. 4.2) and the energy is far too non-linear for anything
    // cheaper: tau(v) alone is quartic in the unknowns.
    int runLBFGS(Eigen::VectorXd &x, int maxIter);

    // Quarter singularities of the *input* directions, as (vertex, quarter
    // turns), measured on the 4-symmetric representation. Available before the
    // comb, which is when the branch cuts have to be chosen.
    std::vector<std::pair<int, int>> inputSingularities() const;

    // A short path of interior edges from each of those singularities to the
    // boundary: the branch cuts the comb is allowed to leave its 90-degree
    // jumps on. Where they run decides where the optimized field puts its
    // boundary corners -- see the implementation.
    std::unordered_set<EdgeKey, EdgeKeyHash> singularitySeams() const;

    // Breadth-first comb of the initial directions over the dual graph,
    // stopping at the cut edges and at the branch cuts. See initialize().
    void combInitialField();

    // Put the combed field in the one homotopy class that has a polysquare:
    // zero total rotation around every hole. A no-op on a simply connected
    // model and whenever the cuts of Sec. 4.1 are being used, since a cut mesh
    // has no hole to wind around. See the implementation -- this is a choice
    // Eq. (1) cannot make for itself, because the classes are separated by
    // fields with an interior zero and descent does not cross those.
    bool unwindHoleHolonomy();   // true if it changed anything
    int holesUnwound = 0;
    // Components the comb's dual graph fell into. One is the healthy case:
    // more means the branch cuts separated a piece of the model and the comb
    // had to guess a branch there, which shows up as an interior defect the
    // optimization cannot shed. See combInitialField().
    int combComponents = 0;
    // Branch cuts sent to the outer loop rather than to the hole they were
    // nearest, so that the hole's charge came out a whole number of turns.
    // See singularitySeams().
    mutable int seamsRerouted = 0;
    // The cap of setMaxCornerQuarters(), which see.
    int maxCornerQuarters = 1;

    // The edge two triangles of a vertex star share, -1 if they share none.
    int sharedEdge(int ta, int tb) const;

    // Fills `chains` from `featureEdges`. Cheap, and done once.
    void buildInterfaceChains() const;

    // One label per vertex: which material region it is strictly inside, or -1
    // on an interface. Empty (and the whole mesh one region) without feature
    // edges, which is what makes the dipole rule a no-op there.
    void buildRegions();

    // The triangles of `v`'s star swept from feature edge `eIn` to feature
    // edge `eOut`, on the left of the walk (far end of eIn) -> v -> (far end
    // of eOut). Empty when the sweep does not reach eOut, which is a star the
    // mesh cannot fan -- a pinch, or a non-manifold vertex.
    std::vector<int> sectorFan(int v, int eIn, int eOut) const;

    // Omega and Theta over such a fan: how far the model subtends at `v` and
    // how far the frame turns across it, the k*90-degree transitions on the
    // cuts undone as everywhere else.
    void sectorMeasure(int v, const std::vector<int> &fan, double &omega, double &theta) const;

    // tau of Eq. (3): the spin-4 representation (cos 4t, sin 4t) written as a
    // quartic polynomial in v, so that it stays differentiable in the
    // unknowns rather than going through the angle.
    static void tau(double x, double y, double &t1, double &t2);
    // Jacobian of tau, which has the form [[p, q], [-q, p]].
    static void tauJacobian(double x, double y, double &p, double &q);

    std::shared_ptr<Mesh> mesh;
    std::unordered_set<EdgeKey, EdgeKeyHash> cutEdges;
    // Interior feature edges, as mesh edge indices and as the keyed set the
    // per-vertex walks below test membership against.
    std::vector<int> featureEdges;
    std::unordered_set<EdgeKey, EdgeKeyHash> featureEdgeKeys;
    // Derived from featureEdges on first use; `mutable` because asking for
    // them changes nothing about the field.
    mutable std::vector<FeatureChain> chains;
    mutable bool chainsBuilt = false;

    // Material region per vertex, for the dipole rule: two singularities may
    // only cancel if they carry the same non-negative label. A vertex whose
    // incident triangles disagree sits *on* an interface and gets -1, which
    // never matches anything.
    std::vector<int> vertexRegion;
    bool cancelDipoles = true;
    mutable std::vector<std::pair<int, int>> dipoles;

    std::vector<Point> initialDirections; // per triangle, before combing
    std::vector<Point> uField;            // per triangle, unit v
    std::vector<Point> vField;            // per triangle, v^perp
    std::vector<double> angles;           // per triangle, atan2(v)

    std::vector<double> faceArea;
    double totalArea = 0.0;
    double totalBoundaryLength = 0.0;

    std::vector<SmoothEdge> smoothEdges;
    std::vector<AlignEdge> alignEdges;

    Eigen::VectorXd x; // stacked (v_x, v_y) per triangle, the unknowns

    double alignWeight = 0.1;  // w_a, Sec. 4.2
    double regWeight = 5.0;    // w_r, Sec. 4.2
    std::vector<double> l1Schedule{1e-2, 1e-3, 1e-4};
    double l1Eps = 1e-2;       // current stage of the schedule
    double gradTolerance = 1e-6;
    int maxIterations = 1000;

    EnergyTerms currentEnergy;
    double currentGradNorm = 0.0;
    double minFrameNorm_ = 1.0;
    int totalIterations = 0;
    bool initialized = false;
};

#endif // __UMBER_HXX__
