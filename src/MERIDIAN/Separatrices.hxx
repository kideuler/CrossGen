#ifndef __SEPARATRICES_HXX__
#define __SEPARATRICES_HXX__

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

#include "MERIDIAN/Immersion.hxx"
#include "mesh/Mesh.hxx"

class LayoutEnergy;

// Stage 7 of Shepherd, Gu and Hughes (2022): the separatrices of Psi, traced.
// (docs/shepherd2022.pdf -- Definition 2.1's property Q5, Fig. 9 for the
// picture, and the opening paragraph of Sec. 4 for where they are used.)
//
// Stage 6 returned Psi : Omega -> R^2 satisfying Q1 to Q5. Q5 is the property
// that makes this stage possible and is worth restating in the form it is used
// here, because it is a statement about curves rather than about the map:
//
//     lines emanating from a singularity with constant u or v coordinate,
//     pulled back to S - G, either terminate at a (possibly identical)
//     singularity, terminate transversely to the boundary, or are transverse
//     to the cutting graph G -- in which case they continue inductively across
//     the cut by the transition of Q4 -- and every such line is finite.
//
// Those are the three cases the marching below has to implement and the two it
// is allowed to end on. The third is not an ending: a curve that reaches an arc
// of G leaves Omega through the child edge on one side and re-enters through
// the child edge on the other, at the same point of S, with its direction
// rotated by that arc's k quarter turns. Fig. 9 is one curve doing exactly
// that, leaving along -du and continuing along -dv on the far side of the cut,
// and returning in the end to the singular point it started from.
//
// ### What comes out
//
// A curve is stored the way Stage 9 will want to read it: an ordered list of
// (triangle, barycentric entry, barycentric exit). Not a polyline of points.
// The reason is that phi is affine on each triangle, so barycentric coordinates
// are the *same numbers* in the image and on the surface; the curve is
// therefore pulled back to S exactly, with no projection and no search, and it
// survives any later refinement of the triangulation that keeps the same
// geometry. A polyline in the image would have to be re-projected, and a
// polyline on the surface would lose the quantisation that took an hour of
// penalty continuation to buy.
//
// Within one triangle the curve is a straight segment in both spaces: the
// isoline of a coordinate of an affine map is a line, and its preimage under
// the affine map is a line. So two barycentric triples per triangle are not an
// approximation of the curve, they are the curve.
//
// ### How many curves leave a cone, and in which directions
//
// A cone of index I has cone angle 2pi - (pi/2) I in the interior and
// pi - (pi/2) I on the boundary, and Q2 says Psi realises that angle exactly.
// Dividing by the quarter turn gives the valence: 4 - I interior, 3 - I on the
// boundary. An interior cone emits one separatrix along each of its 4 - I grid
// directions. A boundary cone emits (3 - I) - 2 = 1 - I, because two of its grid
// directions run along dS itself -- those two are sides of the layout already,
// which is why a convex corner (I = +1, valence 2) emits nothing at all.
//
// The directions are found by the sweep of Sec. 3.4's description: walk the
// one-ring fan of the cone in Omega accumulating angle, and emit a ray each
// time the accumulated angle passes a multiple of pi/2. Doing it by
// accumulation rather than by testing the four axis directions against each
// triangle in turn matters for three reasons, and each of them is a real case
// in this corpus:
//
//   * The fan of an interior cone sweeps more than 2pi. At a valence-five cone
//     the direction +u occurs twice, in two different triangles, and both are
//     separatrices. Only an accumulated angle can tell them apart.
//
//   * The fan of a cone in Omega is not one vertex's fan. A cone that the
//     cutting graph merely stops at has one child and a fan bounded by the two
//     copies of the last cut edge; one the graph runs *through* -- Fig. 7 --
//     has one child per sector, and the sweep has to cross from one child to
//     the next through the seam, rotating by that arc's k as it goes. Since k
//     is a whole number of quarter turns, the axis directions survive the
//     crossing, and the accumulated angle runs on through it uninterrupted.
//
//   * The ray count comes out exactly. Testing directions against triangles
//     puts a ray that lands on a shared fan edge in two sectors or in neither,
//     depending on which way a 1e-16 goes; and after Stage 6 the boundary of
//     the image is axis-aligned only to the constraint tolerance, so rays land
//     on fan edges all the time. Placing 4 - I of them at prescribed
//     accumulated angles cannot miscount.
//
// ### Where a curve stops
//
// Two endings are Q5's, and the third is a diagnosis rather than an ending:
//
//   End::Cone      the curve came within the snap tolerance of a cone. This is
//                  what a quantised layout looks like, and the distance it was
//                  snapped by is the residual E5 was driving to zero.
//   End::Boundary  it left transversely through dS.
//   End::Capped    it ran past the step cap without doing either. That is not
//                  a bug in the tracer: it is Q5 failing to hold, and the
//                  paper's remedy -- raise lambda_5, add a Gamma_topo
//                  constraint for the pair that nearly met, re-run Stage 6
//                  from the current phi rather than from psi_R -- is upstream
//                  of here. Curve::nearestConeGap is how near the curve came
//                  to closing, which is what says whether the constraint is
//                  missing or the continuation merely stopped short.
//
// The snap tolerance is the paper's 1e-6 of the image extent, and the cones it
// is tested against are the ones at the corners of the triangle currently being
// crossed -- not every cone in the model. That restriction is not an
// optimisation. Psi is locally injective and globally *not* (that is Q1, and
// the self-overlap is the magenta hatching of Figs. 10b and 10c), so two
// distant pieces of the surface can share a point of the image, and a cone
// matched by image distance alone can be one the curve is nowhere near on S.
// A cone within 1e-6 of the curve on S is a cone whose one-ring the curve is
// inside, so this test finds every cone the tight tolerance admits and no cone
// it does not.
//
// ### Footnote 3
//
// If S is an annulus with no singularities, Q5 is vacuous and there is nothing
// to trace. The paper's footnote 3 says to call an arbitrary regular point the
// surface's only singularity and proceed unchanged; that point has cone angle
// 2pi, so it emits four separatrices, and the construction goes through. It is
// done here for completeness -- every model in data/meshes is a disk, whose
// sum I = 4 chi = 4 cannot be met by an empty cone set.
class Separatrices {
public:
    using EdgeKey = MeshEdgeKey;
    using EdgeKeyHash = MeshEdgeKeyHash;

    // How a curve ended. The first two are Q5's two allowed endings; the rest
    // are ways the walk did not get there.
    enum class End {
        Cone,        // Q5 case 1: terminated at a (possibly identical) cone
        Boundary,    // Q5 case 2: left transversely through dS
        Capped,      // ran past the step cap; Q5 does not hold yet
        Stuck,       // the walk could not be continued at all
        Degenerate   // the ray never entered its starting triangle
    };

    // Which space a traced curve is asked for.
    enum class Space { Model, Image };

    // One triangle crossing. Barycentric coordinates are with respect to the
    // corners of cutMesh.triangles[face], in that order, and are the same
    // numbers in Omega, in the image, and on S -- see the class comment.
    struct Step {
        int face = -1;
        std::array<double, 3> entry{{0.0, 0.0, 0.0}};
        std::array<double, 3> exit{{0.0, 0.0, 0.0}};
    };

    struct Curve {
        std::vector<Step> steps;

        // Where it came from.
        int cone = -1;        // slot into coneVertices(); -1 never happens
        int child = -1;       // the vertex of Omega it left from
        int index = 0;        // I(cone)
        bool boundaryCone = false;
        int ray = 0;          // which of the cone's rays, in fan-sweep order
        int dir = 0;          // direction of travel at emission, 0..3 = +u,+v,-u,-v
        double fanAngle = 0.0; // accumulated fan angle it left at, radians

        // Where it went.
        End end = End::Capped;
        int endDir = 0;       // direction of travel when it stopped
        int toCone = -1;      // slot, when end == Cone
        int toChild = -1;     // the cone's child it was snapped to
        double gap = 0.0;     // how far the snap moved it, in image units

        // The closest any cone came to the curve over the whole trace, whether
        // or not it was close enough to snap to. On a capped curve this is the
        // measure of how near Q5 came to holding. Infinite when the curve never
        // crossed a triangle with a cone at one of its corners.
        double nearestConeGap = 0.0;
        int nearestCone = -1;

        int seamCrossings = 0;
        double imageLength = 0.0;   // length in the image
        double modelLength = 0.0;   // length of the pullback on S
    };

    struct Options {
        // Sec. 3.4's snap tolerance, as a fraction of the image extent. The
        // paper's table value; loosening it turns near-misses into terminations
        // and is a way of asking "how far is Q5 from holding", not a way of
        // making it hold.
        double coneSnapTolerance = 1e-6;

        // The cap of Sec. 3.4. A curve that reaches it is reported, not
        // truncated silently.
        int maxSteps = 50000;

        // How near the end of a boundary cone's fan a ray has to be before it
        // counts as running along dS rather than into the interior, in radians.
        // After Stage 6 the boundary is axis-aligned to the constraint
        // tolerance rather than exactly, so the two directions that Q3 puts on
        // dS miss the fan ends by ~1e-6 and would otherwise be emitted as
        // separatrices that then run along the boundary for the length of the
        // model.
        double fanEndTolerance = 1e-3;
    };

    struct Report {
        int cones = 0;              // cones that emit at least one separatrix
        int designatedCone = -1;    // footnote 3: the regular point pressed into service

        int emitted = 0;
        int prescribed = 0;         // sum of 4 - I (interior), 1 - I (boundary)

        int endedAtCone = 0;
        int endedAtBoundary = 0;
        int capped = 0;
        int stuck = 0;
        int degenerate = 0;

        int seamCrossings = 0;
        int maxSeamCrossings = 0;
        long long triangleSteps = 0;
        int maxTriangleSteps = 0;

        // Worst snap among the curves that reached a cone: the largest distance
        // any of them had to be moved to land on one, in units of the image
        // extent. This is Q5's residual read off the curves rather than off the
        // Gamma_topo sum E5 minimised.
        double maxSnapGap = 0.0;
        // The closest approach of the curve that ended worst -- capped, or
        // stuck -- again relative to the extent. A small number here with a
        // non-zero capped count is a missing Gamma_topo constraint; a large one
        // is a layout whose cones the curve never went near.
        double worstMissGap = 0.0;
        int worstMissCurve = -1;

        // Q2 re-read off the fan sweep: the accumulated angle around each cone
        // against the 2pi - (pi/2) I (interior) or pi - (pi/2) I (boundary)
        // that Stage 1 prescribed. The sweep is what the ray directions are
        // measured in, so this is the number that says whether the ray count
        // can be trusted.
        double maxConeAngleResidual = 0.0;
        int worstConeAngleSlot = -1;
        int fanFailures = 0;        // cones whose fan could not be swept

        // The largest distance, relative to the extent of S, between where one
        // triangle crossing left off and the next took up, over every curve.
        // Within a sheet the two are the same point by construction; across an
        // arc of G they are the same point only if the seam pairing and the
        // parameter carried across it are right. So this is the number that
        // says the transitions of Q4 were applied correctly, measured on the
        // curves rather than assumed from the arcs -- and it is the one thing
        // in this stage that a wrong k or a reversed child polyline would leave
        // looking plausible in the image and wrong on the model.
        double maxPullbackGap = 0.0;

        double extent = 0.0;        // diagonal of the image
        double snapTolerance = 0.0; // the same in absolute image units

        // Every curve ended the way Q5 allows, every cone emitted the number of
        // separatrices its index prescribes, and every fan swept cleanly.
        bool valid = false;

        std::vector<std::string> messages;
    };

    // Trace on any map of Omega -- Psi from Stage 6, or psi_R from Stage 4 when
    // the later stages were not run. The map must have one point per vertex of
    // the immersion's cut mesh.
    Separatrices(const Immersion &immersion, const std::vector<Point> &map);
    Separatrices(const Immersion &immersion, const std::vector<Point> &map,
                 const Options &opts);

    // The usual form: the layout that came out of Stage 6.
    explicit Separatrices(const LayoutEnergy &layout);
    Separatrices(const LayoutEnergy &layout, const Options &opts);

    const std::vector<Curve>& curves() const { return traced; }
    const Report& getReport() const { return report; }

    const Immersion& getImmersion() const { return *imm; }
    const Mesh& getCutMesh() const { return imm->getCutMesh(); }
    const std::vector<Point>& getUV() const { return uv; }

    // The cones the curves were emitted from: Immersion's list, except in the
    // footnote 3 case where it is the one designated regular point.
    const std::vector<int>& coneVertices() const { return coneVertex; }
    const std::vector<int>& coneIndices() const { return coneIndex; }

    // A point of a step, in the image or pulled back onto S. Both are the same
    // barycentric combination, of Psi and of the vertices of S respectively.
    Point point(const Step &s, bool exit, Space space) const;

    // The whole curve as a polyline, with `breaks` holding the indices at which
    // it is interrupted -- which happens in the image at every seam crossing,
    // where the curve jumps to the other side of the cut, and never on the
    // model, where the two sides are the same point of S.
    std::vector<Point> polyline(const Curve &c, Space space,
                                std::vector<int> *breaks = nullptr) const;

    // The traced curves as an .obj of polylines. On the model this is the left
    // half of the paper's Fig. 9 and the input to Stage 9's spline fit; in the
    // image it is the right half.
    bool writeOBJ(const std::string &filename, Space space) const;

private:
    // One triangle's worth of the fan around a cone, in sweep order.
    struct Sector {
        int face = -1;
        int child = -1;
        int lc = 0;            // local index of the child in face
        double startAngle = 0.0; // image angle of the fan-start direction
        double angle = 0.0;    // the image angle at the child in this face
        double acc = 0.0;      // accumulated angle at the start of this sector
    };

    void buildTables();
    void buildEmitters();
    void measureExtent();
    void traceAll();
    void check();

    // The fan of cone `slot` in Omega, in sweep order, crossing seams as it
    // goes. `closed` says the sweep came back to where it started, which is
    // what an interior cone does and a boundary cone does not.
    bool sweepCone(int slot, std::vector<Sector> &out, bool &closed) const;

    Curve trace(int slot, int child, int face, int dirT, int ray, double fanAngle) const;

    // The face in w's fan whose angular sector contains the axis direction
    // `dirT`, or -1. Used to continue a curve that ran into a vertex head-on,
    // either within the fan or across a seam. `exclude` keeps the walk out of
    // the face it is leaving.
    int faceContaining(int w, int dirT, int exclude) const;

    // Barycentric constructors. Every point the tracer produces is either a
    // vertex or a point at a known parameter along a known edge, so these are
    // exact; baryOfPoint() solves and is used only where a curve is snapped to
    // a cone part way across a triangle.
    std::array<double, 3> baryUnit(int face, int v) const;
    std::array<double, 3> baryFromEdge(int face, int x, int y, double t) const;
    std::array<double, 3> baryOfPoint(int face, const Point &p) const;

    static Point axis(int d) {
        switch (((d % 4) + 4) % 4) {
            case 0: return Point{1.0, 0.0};
            case 1: return Point{0.0, 1.0};
            case 2: return Point{-1.0, 0.0};
            default: return Point{0.0, -1.0};
        }
    }
    static int localEdgeBetween(int m, int n) { return ((m + 1) % 3 == n) ? m : n; }

    const Immersion *imm = nullptr;
    Options options;

    std::vector<Point> uv;      // the map being traced, one point per Omega vertex
    std::vector<Curve> traced;

    // The emitters: Immersion's cones, or footnote 3's designated point.
    std::vector<int> coneVertex;                 // original mesh vertex
    std::vector<int> coneIndex;                  // I(v)
    std::vector<std::vector<int>> coneChildren;  // its vertices in Omega
    std::vector<char> coneOnBoundary;            // v lies in dS
    std::vector<int> vertCone;                   // Omega vertex -> slot, or -1

    // Omega's edges by their endpoints; which of its boundary edges came from
    // dS rather than from a cut; and, for the ones that came from a cut, the
    // seam pair they belong to and the side they are on.
    std::unordered_map<EdgeKey, int, EdgeKeyHash> cutEdgeIndex;
    std::unordered_map<EdgeKey, int, EdgeKeyHash> seamSide;  // 2*pair + (0 plus, 1 minus)
    std::vector<char> parentOnBoundary;      // per Omega edge
    std::vector<char> vertexOnRealBoundary;  // per Omega vertex

    // The seam partner of a vertex, for a curve that arrives at a seam head-on
    // rather than crossing it. A vertex where two arcs meet has no single
    // rotation to be carried across by and is marked ambiguous instead.
    struct VertexSeam {
        int partner = -1;
        int k = 0;          // the arc's quarter turns
        int sign = 0;       // +1 leaving the plus side, -1 leaving the minus side
        bool ambiguous = false;
    };
    std::unordered_map<int, VertexSeam> vertexSeam;

    double extent = 1.0;        // diagonal of the image, Psi(Omega)
    double modelExtent = 1.0;   // diagonal of S, for the pulled-back curves
    double snapTol = 0.0;

    Report report;
};

#endif // __SEPARATRICES_HXX__
