#ifndef __SUBDOMAINLABELS_HXX__
#define __SUBDOMAINLABELS_HXX__

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "MERIDIAN/Immersion.hxx"
#include "MERIDIAN/Interfaces.hxx"
#include "MERIDIAN/Separatrices.hxx"
#include "mesh/Mesh.hxx"

// Stage 5 of Shepherd, Gu and Hughes (2022), Section 3.3: the subdomain
// labelling that the layout-inducing energies of Eqs. (15) to (19) are summed
// over. (docs/shepherd2022.pdf)
//
// Stage 6 minimises sum_j lambda_j E_j, and four of those five energies are
// integrals over *named subsets* of the boundary of Omega. Nothing so far has
// named them. This stage does, once on psi_R and again between penalty steps
// when the map has moved far enough that a different naming is the right one:
//
//   Gamma_u, Gamma_v          boundary curves of dS along which u (resp. v) is
//                             to be held constant.                    -> E2
//   Gamma_u^feat, Gamma_v^feat  the same for the feature chains.       -> E3
//   Gamma_Hol_k               arcs of the cutting graph, by their quarter-turn
//                             transition. Already fixed by Stage 4.   -> E4
//   Gamma_topo                paths between singular points that ought to be
//                             joined by an integral curve.            -> E5
//
// ### Which coordinate a curve wants held constant
//
// The rule of Sec. 3.3 is the obvious one, and the reason it works is that
// psi_R is already an isometry of a flat metric that has been rotated so the
// boundary sits near the axes: a boundary edge that runs mostly horizontally in
// the image is one whose v barely changes, so v is the coordinate to pin, and
// vice versa. Note the naming convention, which is the paper's and is easy to
// read backwards: Gamma_u is the set on which *u is constant*, du/ds = 0, and
// its edges are the ones whose u already varies least.
//
// Boundary edges are labelled one at a time. Feature chains are not: a chain is
// labelled once, from its total flux, because a per-edge vote on a chain that
// wobbles will flip label halfway along and ask E3 to hold u constant on the
// first half and v on the second, which no map satisfies. That constraint is
// not merely hard, it is inconsistent, and the continuation will spend every
// outer step failing to meet it.
//
// ### Connectivity, Gamma_topo
//
// E5 is what Remark 3.1 says is not optional: without it the bracket of Fig. 13
// comes out as 1898 patches of arbitrary aspect ratio instead of 73 usable
// ones. Its input is a set of paths p -> q between singular points that should
// end up joined by an integral curve, each subdivided by G into subcurves
// bounded by cuts, each subcurve carrying a direction label j_k in
// {+u, +v, -u, -v} so that the pieces compose into one consistently oriented
// curve on the branched cover.
//
// The paper takes those paths as input. The practical recipe it gives for
// producing them -- trace the separatrices on the current map, look for one
// that passes close to a cone without hitting it, and pair them -- is what
// seedTopoPaths() does. The subdivision and the direction labels come out of
// the trace for free: a subcurve is what lies between two seam crossings, and
// crossing an arc of Gamma_Hol_k rotates the direction label by k quarter
// turns, which is the whole content of "consistently oriented on the branched
// covering space".
//
// The tracing is Stage 7's, run early: seedTopoPaths() drives Separatrices
// with a loosened snap tolerance rather than marching on its own. That is not
// only economy. A seeder with its own idea of which rays leave a cone finds a
// different set of curves from the one that later checks Q5, and the
// difference is silent -- it shows up two stages later as a separatrix with no
// constraint behind it, missing its cone by a hair and running on. In
// particular the direction extraction has to be the fan sweep of Sec. 3.4 and
// not a test of the four axis directions against each incident triangle: the
// second miscounts at a valence-five cone, where +u occurs twice in the fan,
// and at a boundary vertex, where after Stage 6 the rays land exactly on fan
// edges. See the Separatrices class comment.
//
// Two kinds of curve are worth a constraint and are easy to leave out:
//
//   * one that comes back to the cone it left. Q5 says "terminate at a
//     (possibly identical) singularity", and Fig. 9 is precisely that curve --
//     out along -du, across the cut, back along -dv to where it started. On
//     this corpus these are the majority of the separatrices that fail to
//     terminate, and a seeder that skips them leaves E5 with nothing to say
//     about the one direction that most needs quantising.
//
//   * a second, distinct curve between a pair of cones that are already joined
//     by one. They are different constraints; keeping only the nearer of the
//     two leaves the other free to miss. Hence the deduplication is by the two
//     *ends* -- cone, child and direction -- and not by the pair of cones, so
//     that the same curve traced from both ends still collapses to one.
//
// A subcurve endpoint is generally *not* a vertex -- it is where the curve met
// a seam, part way along an edge -- so it is stored as a Site, an affine
// combination of two vertices with a fixed weight. That keeps E5 linear in the
// unknowns, which is what makes the paper's claim that this avoids
// mixed-integer optimisation true: the constraint is real-valued and quadratic,
// not an integer grid condition.
class SubdomainLabels {
public:
    using EdgeKey = MeshEdgeKey;
    using EdgeKeyHash = MeshEdgeKeyHash;

    // Which coordinate the curve holds constant. None is for a curve that has
    // not been labelled -- a seam, or a feature chain too short to vote.
    enum class Align { None = 0, U = 1, V = 2 };

    // A point of Omega written as an affine combination of at most two
    // vertices, so that phi(site) stays a linear function of the unknowns:
    //     phi = (1 - t) phi(a) + t phi(b),   or phi(a) when b < 0.
    struct Site {
        int a = -1;
        int b = -1;
        double t = 0.0;
    };

    struct BoundaryEdge {
        int cutEdge = -1;      // index into Omega's edge list
        int a = -1, b = -1;    // its two Omega vertices
        double length = 0.0;   // in the flat metric
        Align label = Align::None;
        int chain = -1;
    };

    // A maximal run of boundary edges carrying the same label, bounded by cones,
    // by the ends of the cutting graph, and by label changes. Not used by E2 --
    // which sums edge by edge -- but it is the object Q3 is a statement about
    // and the one worth printing.
    struct BoundaryChain {
        std::vector<int> edges;   // indices into boundaryEdges()
        Align label = Align::None;
        double length = 0.0;
        double spanU = 0.0, spanV = 0.0;  // extent of the chain in the image
    };

    struct FeatureChain {
        std::vector<int> verts;      // Omega vertices in order
        std::vector<double> lengths; // flat length of each edge of the chain
        Align label = Align::None;
        double fluxU = 0.0, fluxV = 0.0;  // |u(end) - u(start)|, likewise v
        double length = 0.0;
        bool onSeam = false;         // the chain runs along an arc of G

        // Where it came from, when the chains were built from a material
        // interface network rather than from the edges alone. `branch` indexes
        // Interfaces::branches(); `dir` is the direction of travel the node
        // quantisation gave it, in the {+u, +v, -u, -v} numbering, and is what
        // `label` is read off. -1 and 0 when there is no network.
        int branch = -1;
        int dir = 0;
        bool dirKnown = false;
    };

    // One sector of one node of the material interface network, as a constraint
    // on Psi: the two curves that bound it, and the number of quarter turns the
    // layout has to put between them. This is what E6 of Stage 6 is summed over.
    //
    // Both curves are read as a tangent at the *same* child of the node in
    // Omega -- the node is one vertex of S but the cutting graph may have split
    // it into several, and a tangent taken at one child and another taken at a
    // different one are not in the same frame. Where the cut runs through the
    // sector itself there is no common child and the sector is dropped, which
    // is what `spansCut` records; Q4 already holds that pair together across
    // the seam, so nothing is lost but the count is worth having.
    struct InterfaceCorner {
        int node = -1;                 // index into Interfaces::nodes()
        int vertex = -1;               // the node's child in Omega
        int a = -1, b = -1;            // the far end of each of the two rays
        double la = 0.0, lb = 0.0;     // their flat lengths
        // The lengths the two tangents are divided by before they are compared.
        // They start at the flat lengths above and are re-read off the current
        // map between outer steps -- see refreshInterfaceScales(), which is
        // where the reason lives.
        double sa = 0.0, sb = 0.0;
        int quarters = 1;              // q_s: the turn from ray a to ray b
        bool onBoundary = false;
        bool aIsBoundary = false, bIsBoundary = false;  // a ray along dS
        double sector = 0.0;           // the angle it has on the model, radians
    };

    // One path of Gamma_topo, already subdivided by G.
    struct TopoPath {
        struct Sub {
            Site from, to;
            // j_k: the direction whose advance is constrained, as an index into
            // {+u, +v, -u, -v}. It is the normal of the integral curve, and it
            // turns by k quarter turns at every crossing of Gamma_Hol_k.
            int dir = 0;
        };
        std::vector<Sub> subs;
        int fromCone = -1, toCone = -1;   // slots in Immersion::getConeVertices()
        double seedGap = 0.0;             // how near the seeding separatrix passed
        int seamCrossings = 0;
        // Which pass put it here: 0 for the seeding on psi_R, then 1, 2, ...
        // for each round of Sec. 3.3's repair on the map as it then stood. A
        // path added late is one the seeding could not see, because psi_R and
        // Psi are far enough apart that a connection visible on one need not be
        // visible on the other.
        int pass = 0;
    };

    struct Options {
        // Seed Gamma_topo by tracing separatrices and pairing near-misses.
        bool seedTopoConstraints = true;
        // How near a separatrix has to pass a cone to count as a near-miss, as
        // a fraction of the mean spacing between cones in the image.
        //
        // Both ends of the range are real. Too small and Remark 3.1's sliver
        // problem survives, because the constraints that would have pinned the
        // layout together were never seeded. Too large and paths get seeded
        // between cones that were never meant to be joined, and E5 then asks
        // for a set of integral curves that no map has: the continuation does
        // not converge slowly, it converges to a compromise, with det J driven
        // towards zero by the barrier fighting a penalty it cannot satisfy.
        // Both failure modes are visible in Report -- the first as separatrices
        // that reach neither a cone nor dS, the second as a Q5 residual that
        // stops falling while lambda_5 keeps rising.
        //
        // 0.15 is the largest value at which every model in data/meshes still
        // converges; 0.2 seeds 27 more constraints across the corpus and makes
        // one model (singlemat/geom007, nine interior cones on 1840 triangles)
        // infeasible.
        double nearMissTolerance = 0.15;
        int maxTraceSteps = 50000;

        // Seed the curve that comes back to the cone it left -- Q5's "possibly
        // identical" case, and the curve the paper draws in Fig. 9. Off is what
        // this code used to do implicitly, and it is worth having as a switch
        // only to be able to see what it costs: on this corpus turning it off
        // is what leaves separatrices winding until the step cap.
        bool seedSelfReturns = true;

        // Keep every distinct connection rather than one per pair of cones.
        // Two separatrices joining the same two cones in different directions
        // are two constraints; collapsing them to one leaves the other free to
        // miss its cone. The same curve traced from both ends still collapses,
        // because the key is the pair of *ends* and not the pair of cones.
        bool seedAllConnections = true;

        // Feature chains are split where the image direction turns by more than
        // this, so that an L-shaped feature is two chains with two labels rather
        // than one chain with an unsatisfiable one.
        double featureCornerAngle = M_PI / 4.0;

        // Material interfaces count as features. These are 2-D midsurface
        // models, so the paper's creases and trim curves arrive here as the
        // boundaries between the physical surfaces of the .geo -- which is what
        // the viewer already draws as feature edges.
        bool materialFeatures = true;

        // Build E6's sector constraints from the interface network's nodes.
        // Only has an effect when the network was handed in; without one there
        // are no nodes and nothing to constrain.
        bool interfaceCorners = true;

        // Take the label of an interface chain from the node quantisation
        // propagated along the network, rather than from the chain's own flux.
        //
        // The flux rule is right for an isolated feature and wrong at a
        // junction, and the failure is not a near miss. Three interfaces meet
        // at a triple point with sectors of one, one and two quarters; the
        // labels that follow are u, v, u -- the third branch collinear with the
        // first. Read off the flux of psi_R instead, all three can come out u,
        // which asks E3 for a map holding u constant on three curves leaving
        // one point in three directions. There is no such map, so the
        // continuation spends every outer step failing to reach it.
        bool propagateInterfaceLabels = true;
    };

    struct Report {
        int boundaryEdges = 0;
        int boundaryEdgesU = 0, boundaryEdgesV = 0;
        int boundaryChains = 0;
        // Edges whose two labels are almost tied, |du| within 20% of |dv|. A
        // large count means the immersion is not near the axes yet and the
        // labelling is close to arbitrary, which is what relabelling between
        // outer steps is for.
        int ambiguousBoundaryEdges = 0;

        int featureEdges = 0;
        int featureChains = 0;
        int featureChainsU = 0, featureChainsV = 0;
        // Chains whose label came from the network's own quantisation rather
        // than from their flux, and interface nodes the propagation could not
        // reach (a component with no chain long enough to seed it).
        int featureChainsPropagated = 0;
        int interfaceCorners = 0;
        int interfaceCornersSpanningCut = 0;
        // Sectors whose two labels disagree with the flux reading. A large
        // count is not an error -- it is the propagation doing its job -- but
        // it is the number that says how much E3 was being asked for before.
        int featureLabelsCorrected = 0;

        int seamArcs = 0;
        int holonomyCount[4] = {0, 0, 0, 0};

        int separatrices = 0;
        int separatricesToBoundary = 0;
        int separatricesToCone = 0;
        int separatricesCapped = 0;

        int topoPaths = 0;
        int topoSelfReturns = 0;    // of those, curves back to their own cone
        int topoExtraPerPair = 0;   // and second-and-later curves of a pair
        int topoAddedByRepair = 0;  // paths added after the seeding, pass > 0
        double meanConeSpacing = 0.0;
        double seedSnapTolerance = 0.0; // what the seeding traced with, absolute
        int seedSnapRings = 0;
        double maxTopoResidual = 0.0;   // worst |sum of advances| as it stands

        std::vector<std::string> messages;
    };

    // Two forms rather than a default argument: Options carries default member
    // initialisers, which a default argument inside the same class definition
    // cannot see yet.
    explicit SubdomainLabels(const Immersion &immersion);
    SubdomainLabels(const Immersion &immersion, const Options &opts);
    // With the material interface network of Stage 0b. `interfaces` is held by
    // pointer and has to outlive this object; null is the single-material case
    // and gives exactly the two-argument behaviour.
    SubdomainLabels(const Immersion &immersion, const Options &opts,
                    const Interfaces *interfaces);

    // Re-read the lengths E6 divides its two tangents by from the map as it now
    // stands.
    //
    // E6 compares b/s_b against R_q (a/s_a), which asserts two things at once:
    // that the two tangents are q quarter turns apart, which is the sector
    // condition and is what the term is for, and that |b|/s_b = |a|/s_a, which
    // is a statement about how much the map stretches the two rays. With s
    // fixed at the flat lengths the second one says the stretch is the same in
    // both directions -- local conformality at the node -- and on a domain
    // whose interfaces meet obliquely that is a real cost: geom013's groove
    // faces meet dS at 58 and 122 degrees and both have to become 90, so the
    // map there *is* anisotropic, and asking it not to be leaves E1 and E6 at
    // an equilibrium with the sector 0.4 rad out and det J at 1e-5.
    //
    // Re-reading s from the current map between outer steps removes that
    // second assertion: at the point it is taken, |a|/s_a = |b|/s_b = 1 by
    // construction, so what is left of the term is the angle. It is the
    // ordinary lagged normalisation of a scale-invariant constraint, and the
    // problem stays exactly quadratic inside each outer step, which is what the
    // proxy needs.
    void refreshInterfaceScales(const std::vector<Point> &uv);

    // Re-derive the Gamma_u / Gamma_v labels from a map that has moved.
    // Sec. 3.3 lists this as an optional step between penalty levels: the
    // labelling made on psi_R can be the wrong one once E2 has pulled the
    // boundary around, and re-reading it is one of the two escapes from a
    // stalled continuation (the other being the reference metric).
    //
    // The seam labels and Gamma_topo are deliberately *not* rebuilt. Q4's k is
    // a topological quantity fixed at Stage 4, and re-seeding the connectivity
    // constraints mid-continuation would change the problem being solved rather
    // than the way it is being solved.
    void relabel(const std::vector<Point> &uv);

    // Build Gamma_topo from the separatrices of the given map. Called once by
    // the constructor when Options::seedTopoConstraints is set.
    void seedTopoPaths(const std::vector<Point> &uv);

    // Add a connectivity constraint by hand, as a path of Omega vertices, with
    // the direction label of its first subcurve. The path is subdivided by G
    // and the remaining labels follow from the transitions it crosses.
    bool addTopoPath(const std::vector<int> &omegaPath, int firstDir);

    // Add the constraint a traced separatrix names: the path from the cone it
    // left to the cone it reached, or -- when it reached none -- to the one it
    // came nearest. Already subdivided by G and already carrying its direction
    // labels, so nothing is recomputed. Returns false when the curve names no
    // pair, or when the pair of ends is one Gamma_topo already has.
    bool addTopoPath(const Separatrices::Curve &c, int pass = 0);

    // Sec. 3.3's repair, and Sec. 4's "on failure" of the arrangement: take the
    // separatrices of the map as it now stands that terminated at neither a
    // cone nor dS, or slipped past one on their way out, and give E5 the
    // constraint each of them names. `window` is how near a curve has to have
    // passed, as a fraction of the image extent. Returns how many were new.
    //
    // This is what the seeding on psi_R cannot do on its own: psi_R and Psi are
    // far apart -- that is the whole point of Stage 6 -- so a connection that
    // is obvious on the finished layout need not have been visible on the map
    // it started from.
    //
    // `maxAdd` caps how many are taken in one pass, nearest first. The cap
    // earns its place on the models whose cones Stage 1 left clustered: there,
    // a dozen constraints switched on at once move the map somewhere none of
    // them is satisfied, while the nearest two or three are satisfiable and,
    // once satisfied, change which of the rest are still wanted. Zero or
    // negative means all of them.
    int adoptCurves(const Separatrices &sep, double window, int pass, int maxAdd = 0);

    // Undo the tail of Gamma_topo, back to `keep` paths, releasing the pairs
    // they had claimed so a later pass may propose them again. The other half
    // of making the repair monotone: a round that came back worse is taken
    // back whole, constraints and map together.
    void truncateTopoPaths(size_t keep);

    const std::vector<BoundaryEdge>& boundaryEdges() const { return bEdges; }
    const std::vector<BoundaryChain>& boundaryChains() const { return bChains; }
    const std::vector<FeatureChain>& featureChains() const { return fChains; }
    const std::vector<InterfaceCorner>& interfaceCorners() const { return iCorners; }
    const Interfaces* getInterfaces() const { return itf; }
    const std::vector<TopoPath>& topoPaths() const { return tPaths; }

    // The seam terms come straight from Stage 4; they are re-exported here so
    // that Stage 6 has one place to read its subdomains from.
    const std::vector<std::array<int, 4>>& seamPairs() const { return imm->getSeamPairs(); }
    const std::vector<int>& seamPairArc() const { return imm->getSeamPairArc(); }
    const std::vector<double>& seamPairLength() const { return imm->getSeamPairLength(); }
    const std::vector<Immersion::Arc>& arcs() const { return imm->getArcs(); }

    // The mean distance from a cone to its nearest other cone, in the given
    // map: the scale a "near miss" is measured against. Static because it is
    // wanted before this object exists -- a caller choosing the near-miss
    // tolerance for a particular model needs to know what a fraction of the
    // spacing comes to in image units.
    static double meanConeSpacing(const Immersion &imm, const std::vector<Point> &uv);

    // phi at a site, for a given map.
    static Point evaluate(const std::vector<Point> &uv, const Site &s) {
        if (s.b < 0) return uv[s.a];
        return uv[s.a] * (1.0 - s.t) + uv[s.b] * s.t;
    }
    // The unit direction e_d for d in {0,1,2,3} = {+u, +v, -u, -v}.
    static Point axis(int d) {
        switch (((d % 4) + 4) % 4) {
            case 0: return Point{1.0, 0.0};
            case 1: return Point{0.0, 1.0};
            case 2: return Point{-1.0, 0.0};
            default: return Point{0.0, -1.0};
        }
    }

    // The current value of the E5 sum for one path: sum over subcurves of the
    // advance of the constrained coordinate. Zero is what Q5 asks for.
    double topoResidual(const std::vector<Point> &uv, const TopoPath &p) const;

    const Immersion& getImmersion() const { return *imm; }
    const Mesh& getCutMesh() const { return imm->getCutMesh(); }
    const Report& getReport() const { return report; }

private:
    void buildSeamTables();
    void buildBoundary(const std::vector<Point> &uv);
    void buildBoundaryChains(const std::vector<Point> &uv);
    void buildFeatures(const std::vector<Point> &uv);
    void labelFeatures(const std::vector<Point> &uv);
    void buildInterfaceCorners();
    // Directions for every ray of the network, propagated from one seed per
    // connected component through the node quantisation. Returns the per-branch
    // direction of travel out of node0, or an empty vector when there is no
    // network to propagate through.
    std::vector<int> propagateDirections(const std::vector<Point> &uv);

    // The mean distance from a cone to its nearest other cone in `uv`, which is
    // the scale a "near miss" is measured against, and the tracer settings that
    // follow from it.
    double coneSpacing(const std::vector<Point> &uv) const {
        return meanConeSpacing(*imm, uv);
    }
    Separatrices::Options seedTracerOptions(const std::vector<Point> &uv,
                                            double spacing) const;

    // The identity of one end of a connection: the cone, which of its children
    // in Omega, and the direction of travel there. A curve and the same curve
    // traced from its far end give the same unordered pair of these, which is
    // what lets the deduplication keep two genuinely different curves between
    // one pair of cones while still collapsing a curve found twice.
    static long long endToken(int cone, int child, int dir) {
        return ((static_cast<long long>(cone) * 1000003LL + child) * 4) +
               (((dir % 4) + 4) % 4);
    }
    static long long topoKeyOf(long long a, long long b);

    const Immersion *imm = nullptr;
    const Interfaces *itf = nullptr;
    Options options;

    std::vector<BoundaryEdge> bEdges;
    std::vector<BoundaryChain> bChains;
    std::vector<FeatureChain> fChains;
    std::vector<InterfaceCorner> iCorners;
    // Direction of travel out of node0 for each branch of the network, in the
    // {+u, +v, -u, -v} numbering, or -1 where the propagation never reached it.
    std::vector<int> branchDir;
    std::vector<TopoPath> tPaths;
    // The pairs of ends already constrained, so that the repair can add to
    // Gamma_topo across several passes without ever adding the same path twice,
    // and the key each path claimed, so truncateTopoPaths() can release them.
    std::unordered_set<long long> topoKeys;
    std::vector<long long> tPathKeys;

    // Omega edge lookup, and the seam partner of every seam child edge.
    std::unordered_map<EdgeKey, int, EdgeKeyHash> cutEdgeIndex;
    // key -> 2 * pairIndex + (0 for the plus child, 1 for the minus child)
    std::unordered_map<EdgeKey, int, EdgeKeyHash> seamSide;
    // Omega boundary edges whose parent lies in dS rather than in a cut, as a
    // fast test.
    std::vector<char> parentOnBoundary;   // per Omega edge

    Report report;
};

#endif // __SUBDOMAINLABELS_HXX__
