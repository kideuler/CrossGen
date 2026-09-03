#ifndef __PAPER_METRICS_HXX__
#define __PAPER_METRICS_HXX__

#include <complex>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "mesh/Mesh.hxx"

// Sec. 5 of the outline, once, so that every experiment measures the same
// thing. Nothing in here knows which method produced the field it is handed:
// a cross field is one unit spin-4 value per triangle and that is the whole of
// the interface, which is what makes the comparison in E2 to E5 a comparison
// rather than three different measurements.
namespace paper {
namespace metrics {

// --- the DualMBO edge weights, which every energy is measured with -------------

// kappa_e = gamma |e| / h_e with h_e = min_K 2|K| / |e|, exactly as
// DualMBO::initialize assembles it. Zero on boundary edges, which carry no
// coupling term.
struct EdgeWeights {
    std::vector<double> kappa;  // per edge
    std::vector<double> area;   // per triangle
    double gamma = 10.0;
};
EdgeWeights edgeWeights(const Mesh &m, double gamma = 10.0);

// E_face(u) = 1/2 sum_{e interior} kappa_e |u_i - P_e u_j|^2, P_e = 1 in the
// plane. `skipEdge`, when non-empty, drops the edges flagged in it -- that is
// how a multi-material energy is compared fairly, since the interface edges
// carry no coupling in the DualMBO operator and counting them would charge the
// method for a jump it deliberately allows.
double commonEnergy(const Mesh &m, const EdgeWeights &w, const Eigen::VectorXcd &u,
                    const std::vector<char> &skipEdge = {});

// The *whole* discrete energy the DualMBO-MBO scheme descends, Dirichlet data
// included:
//
//   E(u) = 1/2 sum_{e interior, not aligned} kappa_e |u_i - u_j|^2
//        + 1/2 sum_{e on dS}                 kappa_e^1 |u_{K(e)} - g_e|^2
//        + 1/2 sum_{e aligned, each side}    kappa_e^1 |u_K - g_e|^2
//
// with g_e = exp(4 i theta_e) the edge tangent's spin-4 value, kappa_e =
// gamma |e| / min_K h_{K,e} on a coupled edge and kappa_e^1 = gamma |e| / h_{K,e}
// on a one-sided one -- which is exactly what DualMBO::initialize assembles into K
// and b, since |u| = |g_e| = 1 makes 1/2<u,Ku> - Re<u,b> this sum less a
// constant.
//
// `commonEnergy` above is the *comparison* functional of Sec. 5 and deliberately
// omits the boundary terms, because a method that does not impose boundary data
// the same way should not be charged for it. This one is the energy of the
// method itself, and it is what E1(d) watches along the MBO.
double dualMBOEnergy(const Mesh &m, double gamma, const Eigen::VectorXcd &u,
                  const std::vector<char> &alignedEdge = {});

// --- singularities ---------------------------------------------------------

struct Singularity {
    int vertex = -1;
    int index4 = 0;          // 4 x the cross index, so +1 is a +1/4 cone
    bool onBoundary = false;
    double distanceToBoundary = 0.0;  // Euclidean, on the model
    double localEdgeLength = 0.0;     // mean length of the edges at the vertex
};

struct SingularityReport {
    std::vector<Singularity> interior;   // index4 != 0 at an interior vertex
    std::vector<Singularity> boundary;   // index4 != 0 at a boundary vertex

    int interiorIndex4Sum = 0;      // sum of 4 I(v) over interior vertices
    int boundaryIndex4Sum = 0;      // sum of 4 I(v) over boundary vertices, from the field
    int boundaryQuarterSum = 0;     // sum of the corner quarter counts, from the geometry
    int chi = 0;                    // V - E + F of the triangulation

    // Prop. 3, as an integer:
    //
    //     sum_{v interior} 4 I(v) + sum_{v in dS} 4 I(v) - 4 chi(Omega)
    //
    // with the boundary index read off the *field* -- the fan closed with the
    // two boundary rays and the tangent's turn across the exterior wedge. This
    // is the identity, so a non-zero value is the extraction failing (a wrap
    // lost where consecutive faces disagree by more than pi in spin space),
    // not the field being wrong.
    int poincareHopfResidual = 0;

    // The same sum with the boundary read off the *geometry* instead: q(v) is
    // the number of quarter turns the boundary makes at v, which is what a
    // boundary-aligned field carries there when it follows the boundary and
    // nothing else. Non-zero means the field did something else somewhere on
    // dS -- put a cone on it, or rounded a corner the other way -- and
    // `boundaryAnomalies` says at how many vertices.
    int geometricPoincareHopfResidual = 0;

    // Singularities within graph distance 1 of a boundary vertex: the ones the
    // outline names as the ones that typically break tracing.
    int boundaryAdjacent = 0;
    // The nearest singularity to dS, in units of the local edge length there.
    double minBoundaryDistanceH = 0.0;
    // Pairs of singularities closer than twice the local edge length -- the
    // "spurious pair" proxy -- and how many of those are +1/-1.
    int closePairs = 0;
    int closeOppositePairs = 0;

    // Sec. 4.6's extraction condition, checked rather than assumed: the largest
    // spin-4 angle jump between consecutive faces of any vertex star, and how
    // many stars have one at or past pi/4. A star past the condition is one
    // where the winding is a rounding rather than a reading, and the outline
    // asks for it to be reported per run.
    double maxSpinJump = 0.0;
    int starsPastCondition = 0;

    // Boundary vertices at which the field's own index differs from the quarter
    // count the geometry asks for.
    int boundaryAnomalies = 0;
    // Corners whose quarter count is nearly a tie -- |2 tau / pi - round(...)|
    // above 0.4, so a corner within about nine degrees of 45 or 135. At such a
    // corner "what the boundary asks for" is genuinely ambiguous and the two
    // residuals above are entitled to disagree. Sec. 7 lists this as a
    // limitation; this counts it.
    int ambiguousCorners = 0;
    // Boundary vertices the closure could not be applied at (a pinched boundary,
    // or a star with no single gap) and boundary vertices where the closed total
    // was not a multiple of 2 pi. Both break the identity above, and both are
    // counted so that a non-zero residual can be attributed rather than
    // wondered at.
    int boundaryUnclosed = 0;
    int boundaryNonIntegral = 0;

    // The multiset of indices, as a printable "+1x4 -1x2" style string.
    std::string indexMultiset() const;
};

// The interior angle at each vertex: the tip angles of the incident triangles,
// added up. 2 pi at an interior vertex of a flat mesh; less or more at a
// boundary one, which is what makes it the corner angle there.
std::vector<double> interiorAngles(const Mesh &m);

// The number of quarter turns the boundary makes at each boundary vertex, from
// the interior angle alone: round(2 tau / pi) with tau = pi - a the exterior
// turn, which is +1 at a convex right-angle corner, 0 on a straight run, -1 at
// a reflex one and -2 at a spike. The same table CrossField's Dirichlet data
// uses (Viertel, Osting and Staten, Table 1).
std::vector<int> cornerQuarters(const Mesh &m);

SingularityReport singularities(const Mesh &m, const Eigen::VectorXcd &u);

// --- alignment -------------------------------------------------------------

struct AlignmentError {
    int edges = 0;
    double maxAngle = 0.0;      // radians, in [0, pi/4]
    double p95Angle = 0.0;
    double meanAngle = 0.0;
};

// The angle between an edge's tangent and the nearest of the four cross axes of
// each incident face, over the edges flagged in `which`. `oneSided` takes each
// incident face separately, which is what an interface needs -- the two sides
// carry different crosses by design.
AlignmentError alignmentError(const Mesh &m, const Eigen::VectorXcd &u,
                              const std::vector<int> &which);

// The same over dS, where there is one incident face per edge.
AlignmentError boundaryAlignmentError(const Mesh &m, const Eigen::VectorXcd &u);

// --- the B1 conversion, measured -------------------------------------------

struct ConversionReport {
    // Faces whose three corner values summed to nearly zero, so the average has
    // no direction: the conversion is undefined there and takes the corner
    // nearest the centroid instead.
    int degenerateFaces = 0;

    // Sec. 6's key number. The P1 field is continuous, so the transport of its
    // angle from one face centroid to the next through the shared edge is
    // determined; the resampled face field's own jump is that transport reduced
    // to (-pi, pi]. Where the two differ by a full turn, the resampling lost a
    // quarter turn of the cross on that edge and the matching the quantizer
    // reads there is not the P1 field's holonomy.
    int holonomyLostEdges = 0;
    // Edges whose path passed through a near-zero of the interpolated field, so
    // the transport could not be tracked. These are the edges of a P1 singular
    // face and are reported rather than counted above.
    int untrackableEdges = 0;

    // The topological content, before and after. `before` is the P1 field's own
    // singular faces mapped to the vertex of that face nearest its centroid;
    // `after` is the vertex winding of the converted face field.
    int singularitiesBefore = 0;
    int singularitiesAfter = 0;
    int indexSum4Before = 0;
    int indexSum4After = 0;
    // Greedy nearest matching of the two sets at equal index within three local
    // edge lengths; what is left over on each side.
    int created = 0;      // in `after`, unmatched
    int annihilated = 0;  // in `before`, unmatched
    int topologyChanged = 0; // created + annihilated
};

// --- quad-mesh metrics for E4 and E5 ---------------------------------------

struct QuadMetrics {
    int vertices = 0;
    int quads = 0;
    int materials = 0;

    // R4. An interior vertex is regular at valence 4 and a boundary vertex at
    // valence 2 (in quads), so anything else is irregular.
    int irregularInterior = 0;
    int irregularBoundary = 0;
    std::vector<int> irregularPerMaterial;  // by material id - 1, interior only

    // R1. Quads whose four corners do not all carry the same material.
    int mixedQuads = 0;

    double minScaledJacobian = 0.0;
    double meanScaledJacobian = 0.0;
    int invertedQuads = 0;
    double maxAngleDeviation = 0.0;   // degrees from 90
    double meanAngleDeviation = 0.0;

    // R3/R4, the CFL proxy. The characteristic length of a zone as a
    // Lagrangian hydro code takes it: |K| / max diagonal, which is the width of
    // the zone across its long way and is what the sound-speed time step
    // divides. Normalised by the target edge length so the numbers are
    // comparable across models.
    double minZoneDimension = 0.0;
    double p01ZoneDimension = 0.0;   // 1st percentile
    double meanZoneDimension = 0.0;
};

// `quadMaterial` may be empty, in which case every quad is material 1.
QuadMetrics quadMetrics(const std::vector<Point> &verts,
                        const std::vector<std::array<int, 4>> &quads,
                        const std::vector<int> &quadMaterial,
                        double targetEdge);

// --- small shared helpers --------------------------------------------------

// The cross angle of a spin-4 value, in [0, pi/2).
inline double crossAngle(const std::complex<double> &u) {
    double a = std::arg(u) / 4.0;
    while (a < 0.0) a += M_PI_2;
    while (a >= M_PI_2) a -= M_PI_2;
    return a;
}

// The angle between the cross `u` and the direction `theta`, in [0, pi/4].
inline double crossToDirection(const std::complex<double> &u, double theta) {
    const std::complex<double> g = std::exp(std::complex<double>(0.0, 4.0 * theta));
    return std::fabs(std::arg(u * std::conj(g))) / 4.0;
}

// The mean edge length at each vertex.
std::vector<double> localEdgeLength(const Mesh &m);

// V - E + F.
int eulerCharacteristic(const Mesh &m);

} // namespace metrics
} // namespace paper

#endif // __PAPER_METRICS_HXX__
