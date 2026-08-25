#ifndef __SUBDOMAINLABELS_HXX__
#define __SUBDOMAINLABELS_HXX__

#include <string>
#include <unordered_map>
#include <vector>

#include "MERIDIAN/Immersion.hxx"
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

        // Feature chains are split where the image direction turns by more than
        // this, so that an L-shaped feature is two chains with two labels rather
        // than one chain with an unsatisfiable one.
        double featureCornerAngle = M_PI / 4.0;

        // Material interfaces count as features. These are 2-D midsurface
        // models, so the paper's creases and trim curves arrive here as the
        // boundaries between the physical surfaces of the .geo -- which is what
        // the viewer already draws as feature edges.
        bool materialFeatures = true;
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

        int seamArcs = 0;
        int holonomyCount[4] = {0, 0, 0, 0};

        int separatrices = 0;
        int separatricesToBoundary = 0;
        int separatricesToCone = 0;
        int separatricesCapped = 0;

        int topoPaths = 0;
        double meanConeSpacing = 0.0;
        double maxTopoResidual = 0.0;   // worst |sum of advances| as it stands

        std::vector<std::string> messages;
    };

    // Two forms rather than a default argument: Options carries default member
    // initialisers, which a default argument inside the same class definition
    // cannot see yet.
    explicit SubdomainLabels(const Immersion &immersion);
    SubdomainLabels(const Immersion &immersion, const Options &opts);

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

    const std::vector<BoundaryEdge>& boundaryEdges() const { return bEdges; }
    const std::vector<BoundaryChain>& boundaryChains() const { return bChains; }
    const std::vector<FeatureChain>& featureChains() const { return fChains; }
    const std::vector<TopoPath>& topoPaths() const { return tPaths; }

    // The seam terms come straight from Stage 4; they are re-exported here so
    // that Stage 6 has one place to read its subdomains from.
    const std::vector<std::array<int, 4>>& seamPairs() const { return imm->getSeamPairs(); }
    const std::vector<int>& seamPairArc() const { return imm->getSeamPairArc(); }
    const std::vector<double>& seamPairLength() const { return imm->getSeamPairLength(); }
    const std::vector<Immersion::Arc>& arcs() const { return imm->getArcs(); }

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

    // One traced separatrix. Terminates at a cone, at dS, or at the step cap.
    struct Trace {
        std::vector<TopoPath::Sub> subs;
        int hitCone = -1;      // cone slot, or -1
        double gap = 0.0;      // how far the snap moved the curve
        bool leftDomain = false;
        bool capped = false;
        bool degenerate = false;   // the ray never entered its starting triangle
        int seamCrossings = 0;
    };
    Trace traceSeparatrix(const std::vector<Point> &uv, int startVertex, int face,
                          int tangentDir, double coneTol) const;

    // Local edge of face f joining local corners m and n.
    static int localEdgeBetween(int m, int n) { return ((m + 1) % 3 == n) ? m : n; }

    const Immersion *imm = nullptr;
    Options options;

    std::vector<BoundaryEdge> bEdges;
    std::vector<BoundaryChain> bChains;
    std::vector<FeatureChain> fChains;
    std::vector<TopoPath> tPaths;

    // Omega edge lookup, and the seam partner of every seam child edge.
    std::unordered_map<EdgeKey, int, EdgeKeyHash> cutEdgeIndex;
    // key -> 2 * pairIndex + (0 for the plus child, 1 for the minus child)
    std::unordered_map<EdgeKey, int, EdgeKeyHash> seamSide;
    // Omega boundary edges whose parent lies in dS, as a fast test, and the
    // same question asked of a vertex -- which the tracer needs when a curve
    // runs into a boundary vertex head-on rather than crossing an edge.
    std::vector<char> parentOnBoundary;   // per Omega edge
    std::vector<char> vertexOnRealBoundary;

    Report report;
};

#endif // __SUBDOMAINLABELS_HXX__
