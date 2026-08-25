#ifndef __IMMERSION_HXX__
#define __IMMERSION_HXX__

#include <memory>
#include <string>
#include <vector>

#include "MERIDIAN/ConeCut.hxx"
#include "MERIDIAN/ConeSingularities.hxx"
#include "MERIDIAN/RicciFlow.hxx"
#include "mesh/Mesh.hxx"

// Stage 4 of Shepherd, Gu and Hughes (2022), the second half of Section 3.2.2:
// the metric immersion psi_R : Omega -> R^2. (docs/shepherd2022.pdf; the
// underlying result that this map is locally injective is Jin et al. [80].)
//
// Stage 3 left a flat cone metric -- one length per edge, curvature zero
// everywhere except (pi/2) I(v) at the cones -- and Stage 2 left the cut disk
// Omega it lives on. Neither is a map. This stage makes one, by the only
// construction a flat metric admits: lay one triangle down in the plane and
// unfold its neighbours around it, each triangle placed from the two vertices
// it shares with the one before and its own three lengths.
//
//     place seed triangle t0 (one edge on the u axis)
//     while the queue is not empty:
//         t <- pop; (i, j) the already-placed edge; k the free vertex
//         p_k = the circle-circle intersection |p_k - p_i| = l_ik,
//                                              |p_k - p_j| = l_jk
//         take the root on the far side of ij from the triangle t came from
//
// Path-independence is what makes this well defined, and it is bought by the
// three preceding stages together: the metric is flat away from the cones
// (Stage 3), Omega is a disk (Stage 2), and the cuts stop at the cones rather
// than running through them (Stage 2 again), so no closed loop in Omega
// encircles a cone and there is no holonomy to accumulate. Drop any one of the
// three and the unfolding depends on the order the queue happened to take.
//
// The map overlaps itself. It is locally injective and globally not -- the
// magenta hatching of the paper's Figs. 10b and 10c -- because the cone angles
// are 2pi - (pi/2) I and a disk with a 3pi cone in it does not fit in the plane
// without folding over. That is the expected output, not a failure; what would
// be a failure is a *flipped triangle*, which is a different thing and is
// counted separately.
//
// ### What comes out, and what the next stages do with it
//
// Sec. 3.2.2's summary: for a genus-zero surface with boundary, psi_R already
// satisfies three of the five properties of Definition 2.1 --
//
//   Q1  local injectivity, from [80] and from the flat metric;
//   Q2  vertex angle sums are 2pi, pi or a multiple of pi/2 at the cones,
//       because that is what Stage 3 drove the metric to;
//   Q4  the transition across each arc of the cutting graph is a rotation by a
//       multiple of pi/2 plus a translation.
//
// -- and does not satisfy Q3 (each component of dS - G on a constant-coordinate
// line) or Q5 (integral curves between cones are finite). Those two are the
// whole job of Stage 6, and E1's barrier can only *keep* Q1, never restore it,
// which is why Stage 6 has to start from a map that already has it.
//
// Q3 is worth a second look, because on these models it usually arrives a stage
// early and it is better to know why than to be surprised by a residual of
// 1e-15 where Sec. 3.3 leads one to expect work. Stage 1 prescribes Kbar = 0 at
// every non-cone vertex, the ones on dS included, so Stage 3 drives the whole
// boundary geodesic and concentrates all of its turning at the cones -- where
// it is pi - (pi/2) I, a multiple of pi/2, because the indices are integers.
// The boundary that comes out is a rectilinear polygon, and alignToAxes() puts
// it on the grid. E2 is then already satisfied and Stage 6's real work is Q5.
//
// It is not automatic in general. A boundary chain that the cutting graph
// interrupts, a cone the field placed in the middle of a straight run, or a
// flow stopped short of its tolerance all leave E2 something to do, and
// LayoutEnergy::Report::initialBoundaryResidual is where that shows up.
//
// Q4 is not merely satisfied here, it is *measured* here: the transition of
// each arc is fitted, its rotation snapped to the nearest quarter turn, and the
// arc filed under Gamma_Hol_k. E4 of Stage 6 then holds the map to that same k
// while everything else moves. Getting k wrong is not recoverable later, so the
// snap error is reported per arc.
//
// ### Arcs
//
// An "arc of G" here is a maximal chain of cutting-graph edges whose interior
// vertices have exactly two G-edges and lie off dS. Junctions of the graph and
// the points where it meets the boundary end an arc; so does a cone, which is a
// leaf and therefore a node by the same rule. That decomposition -- rather than
// ConeCut's list of routed paths -- is what the transition fit needs, because a
// later cone path terminating in the middle of an earlier one splits the
// earlier path into two arcs with two different transitions, and fitting one
// rigid motion across the junction would fit neither.
//
// Each arc carries the ordered pair of child polylines it has in Omega, in
// consistent arc-length orientation: walking the arc from its first node to its
// last, the "+" side is the one to the left of the direction of travel, read
// off the orientation of the incident face. That is the (e+, e-) pairing table
// the seam terms of Stages 5, 6 and 7 are all written against.
class Immersion {
public:
    using EdgeKey = MeshEdgeKey;
    using EdgeKeyHash = MeshEdgeKeyHash;

    // One arc of G, with its seam pairing and its fitted transition.
    struct Arc {
        std::vector<int> parentPath;   // original vertex ids, node to node
        std::vector<int> plusChain;    // Omega vertices on the left of travel
        std::vector<int> minusChain;   // Omega vertices on the right

        // Gamma_Hol_k: the transition taking the minus side to the plus side is
        // psi+ = R_k psi- + t, R_k the rotation by k pi/2.
        int k = 0;
        Point translation{0.0, 0.0};

        double theta = 0.0;       // the fitted rotation before the snap
        double snapError = 0.0;   // |theta - k pi/2|, in radians
        double fitResidual = 0.0; // rms of |R_k psi-(v) + t - psi+(v)| over the arc
        double length = 0.0;      // flat length of the arc

        bool touchesCone = false;     // an endpoint of the arc is a cone
        bool degenerate = false;      // the two sides did not actually separate
    };

    struct Report {
        int cutVertices = 0;
        int placedVertices = 0;
        int unplacedVertices = 0;      // > 0 means Omega was not connected

        // A face whose three flat lengths do not satisfy the triangle
        // inequality: the free vertex has no circle-circle intersection at all
        // and was put at the foot of the perpendicular instead. Non-zero only
        // when the flow flipped an edge away from the input triangulation, and
        // then only on the faces the flip was made for.
        int degenerateFaces = 0;
        int recoveredEdgeLengths = 0;  // input edges refilled after a flip

        // Q1. The immersion is expected to overlap itself; it is not expected
        // to fold. Every face should carry the same orientation sign.
        int flippedFaces = 0;
        double minSignedAreaRatio = 0.0; // min over faces of (image area / flat area)

        // The isometry check: the largest relative disagreement between
        // |psi(i) - psi(j)| and the flat length l_ij, over every edge of Omega.
        // This is what says the unfolding actually realised the metric.
        double maxMetricResidual = 0.0;

        // Path-independence, measured rather than assumed. A vertex reached by
        // the queue a second time is re-derived from a different pair of placed
        // neighbours; on a flat metric over a disk with no cone inside any loop
        // the two answers agree, and this is how far apart they were, relative
        // to the local edge length. A number that is not at rounding level says
        // a cone got enclosed -- which means the cut ran through one.
        double maxClosureGap = 0.0;

        // Q2, read directly off the map: at every cone, the angles of the
        // incident image triangles summed over all its children, compared with
        // the nearest multiple of pi/2. Stage 3 drove the *metric* to this; the
        // number here is what survived the unfolding.
        double maxConeAngleResidual = 0.0;
        double maxRegularAngleResidual = 0.0; // the same at non-cone vertices

        // Q4.
        int arcs = 0;
        int seamEdgePairs = 0;
        int holonomyCount[4] = {0, 0, 0, 0};  // arcs in Gamma_Hol_0 .. Gamma_Hol_3
        double maxSnapError = 0.0;            // worst |theta - k pi/2|
        double maxFitResidual = 0.0;
        int degenerateArcs = 0;
        int boundaryArcEdges = 0;   // G edges that lie in dS and cannot be paired

        // Sec. 3.2.2's closing assumption: the immersion is rotated so that at
        // least one boundary edge is exactly axis aligned.
        double globalRotation = 0.0;
        int alignedBoundaryEdge = -1;

        // Cones that Stage 1 put on top of each other. Sec. 3.1 lists
        // clustering as one of the three ways automatic placement fails, and
        // its remedy is to merge the cluster into one cone of the summed index.
        // Nothing downstream can do that job -- the cone set is an input by the
        // time the metric exists -- but it is worth naming here, because a
        // cluster is what a Q2 residual out of Stage 6 usually turns out to be:
        // the boundary edge between two adjacent cones is a thousandth of the
        // model long, and E2 aligning it to a coordinate axis leaves an angular
        // error of a tenth of a radian at both of them.
        int clusteredConePairs = 0;
        double minConeSeparation = 0.0;   // relative to the image extent

        Point uvMin{0.0, 0.0};
        Point uvMax{0.0, 0.0};
        bool mirrored = false;   // the seed happened to lay the disk down reversed

        bool valid = false;   // placed everything, no folds, metric realised

        std::vector<std::string> messages;
    };

    // The three inputs have to agree with each other: the cut and the flow must
    // have been built on the same mesh, and the cones must be the ones the flow
    // was driven to. The constructor checks and throws rather than producing a
    // layout of a metric that belongs to something else.
    Immersion(const ConeCut &cut, const RicciFlow &ricci, const ConeSingularities &cones);

    // psi_R, one planar point per vertex of Omega.
    const std::vector<Point>& getUV() const { return uv; }

    // Omega itself, and the mesh it was cut from.
    const Mesh& getCutMesh() const { return cutMesh; }
    const Mesh& getOriginalMesh() const { return *orig; }
    const ConeCut& getCut() const { return *cutter; }
    const ConeSingularities& getCones() const { return *cones; }

    // The flat metric in two indexings: per edge of the input mesh, and per
    // edge of Omega. Stage 6 uses the second as the reference metric of the
    // Jacobian, and Stages 5 and 6 both weight their boundary integrals by it.
    const std::vector<double>& getFlatEdgeLengths() const { return flatLen; }
    const std::vector<double>& getCutEdgeLengths() const { return cutLen; }

    // The arcs of G with their seam pairings and transitions.
    const std::vector<Arc>& getArcs() const { return arcs; }

    // Every (e+, e-) child edge pair of every arc, flattened: entry 4i..4i+3 is
    // (i+, j+, i-, j-) as Omega vertices, walking the arc forwards. arcOfPair()
    // says which arc each came from, and so which R_k applies.
    const std::vector<std::array<int, 4>>& getSeamPairs() const { return seamPairs; }
    const std::vector<int>& getSeamPairArc() const { return seamPairArc; }
    const std::vector<double>& getSeamPairLength() const { return seamPairLength; }

    // The rotation by k pi/2 as a 2x2 action, for the seam terms downstream.
    static Point rotateQuarter(const Point &p, int k) { return rotateVector(p, ((k % 4) + 4) % 4); }

    // Children in Omega of each cone of P, and the cone index carried alongside.
    // A cone left as a leaf of G has exactly one child; one the graph ran across
    // has one per sector.
    const std::vector<std::vector<int>>& getConeChildren() const { return coneChildren; }
    const std::vector<int>& getConeIndices() const { return coneIndex; }
    const std::vector<int>& getConeVertices() const { return coneVertex; }
    // Cone number at an Omega vertex, or -1. Indexes the three vectors above.
    const std::vector<int>& getCutVertexCone() const { return cutVertCone; }

    const Report& getReport() const { return report; }

    // Omega laid out at psi_R, as a flat .obj. This is the picture of Figs. 10b
    // and 10c, self-overlaps included.
    bool writeOBJ(const std::string &filename) const;

private:
    void buildLengths(const RicciFlow &ricci);
    void layout();                 // the unfolding
    void normaliseOrientation();   // make every image triangle positively oriented
    void buildArcs();              // maximal chains of G between its nodes
    void fitTransitions();         // Sec. 3.2.2's least-squares rigid fit per arc
    void alignToAxes();            // the global rotation
    void check();

    // Place the free vertex of a triangle from the two placed ones. `away` is a
    // point that p_k must end up on the other side of the line ij from.
    // Returns false when the three lengths do not close, in which case p_k is
    // placed on the line ij and the caller counts it.
    bool placeThird(const Point &pi, const Point &pj, const Point &away,
                    double lik, double ljk, Point &pk) const;

    const ConeCut *cutter = nullptr;
    const ConeSingularities *cones = nullptr;
    std::shared_ptr<Mesh> orig;
    const Mesh &cutMesh;

    std::vector<double> flatLen;   // per edge of the input mesh
    std::vector<double> cutLen;    // per edge of Omega
    std::vector<Point> uv;

    std::vector<Arc> arcs;
    std::vector<std::array<int, 4>> seamPairs;
    std::vector<int> seamPairArc;
    std::vector<double> seamPairLength;

    std::vector<std::vector<int>> coneChildren;
    std::vector<int> coneIndex;
    std::vector<int> coneVertex;
    std::vector<int> cutVertCone;

    Report report;
};

#endif // __IMMERSION_HXX__
