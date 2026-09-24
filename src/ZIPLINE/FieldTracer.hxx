#ifndef __FIELD_TRACER_HXX__
#define __FIELD_TRACER_HXX__

#include <array>
#include <memory>
#include <vector>

#include "crossfield/CrossField.hxx"

// ---------------------------------------------------------------------------
// Streamline tracing for a cross field carried as one unit complex number per
// vertex, u_v = e^{4 i theta_v}.
//
// The one thing that makes tracing a cross field harder than tracing a vector
// field is that theta is known only modulo pi/2, so every step has to decide
// which of the four directions of the cross continues the previous one, and a
// single wrong decision turns the streamline through a right angle and
// invalidates everything after it. The decision is never allowed to be a close
// call here:
//
//   * Inside one triangle the three vertex angles are lifted to representatives
//     that agree pairwise to within pi/4 -- which is possible exactly when the
//     triangle is not singular -- and the field is the affine interpolant of
//     that lift. The branch is chosen once, on entry, and held for the whole
//     crossing, so the direction inside a triangle is an ordinary smooth vector
//     field with no ambiguity left in it.
//
//   * Across an edge the two triangles' lifts differ by a constant multiple of
//     pi/2 along the shared edge, and not by anything that varies along it:
//     the relative lift of the two shared vertices is the principal matching of
//     that edge, a property of the edge rather than of either triangle, so it
//     is the same number on both banks. The angle carried out of one triangle
//     and the angle read at the same point in the next therefore differ by an
//     exact multiple of pi/2 up to rounding, and recovering that multiple is a
//     rounding decision with a margin of pi/4 rather than a margin of whatever
//     the field happened to do between two vertices.
//
// That is what stands in place of the usual chain of "if this direction looks
// wrong, try the previous one" fallbacks: there is nothing to fall back from.
// The direction handed forward at an exit point is always the field read on
// the edge, i.e. from the two endpoint angles of the edge alone, so the
// property above holds for every hand-off including the ones out of singular
// triangles, where the interpolant of the triangle does not exist.
//
// Inside a singular triangle the field is not the interpolant of anything --
// the lift is inconsistent around the triangle, which is what makes it
// singular -- and streamlines there are computed from the local model of
// Viertel, Osting and Staten, IMR 2019, Sec. 3.2.2: with the singularity at the
// origin and theta measured from a separatrix, f(z) = e^{i d theta / 4}, and
// under w = z^{(4-d)/8} the streamlines are the hyperbolas x y = A of the first
// quadrant. Points along them are evaluated on a fixed fan of rays out of the
// singularity, shared by every streamline crossing that triangle, which is what
// makes two of them unable to cross tangentially there (Prop. 1 and the
// paragraph after it).
// ---------------------------------------------------------------------------

// A point of a traced streamline. Consecutive points bound a segment lying
// inside a single triangle; `face_id` is that triangle for the segment that
// *ends* at this point, and for the first point of a path it is the triangle
// the path starts out into.
struct TracePoint {
    Point global_pos{0.0, 0.0};
    int face_id = -1;
    double theta = 0.0;  // travel direction on arrival: a genuine angle, not a cross angle
};

// A singularity of the cross field and the directions its separatrices leave
// in. Ports are absolute angles, sorted counter-clockwise and spaced by
// 2 pi / (4 - d).
struct Singularity {
    int triangleIndex = -1;
    double singularityIndex = 0.0;  // d/4: +1/4 for a 3-port, -1/4 for a 5-port
    int d = 0;                      // 4 * index, an integer
    Point coordinates{0.0, 0.0};
    std::array<double, 3> barycentric{{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}};
    std::vector<double> portAngles;
    std::vector<int> portSeparatrixIds;  // filled in by SeparatrixTrace
    // Spread of the per-vertex estimates of the port phase, in radians. Small
    // means the three corners of the singular triangle agree about where the
    // separatrices leave; large means the field is not resolved there and the
    // ports are a guess.
    double portResidual = 0.0;

    int numPorts() const { return static_cast<int>(portAngles.size()); }
};

// Where a streamline is and which way it is going. A walk sits on the boundary
// of the triangle it is about to cross, except at the very first point of a
// path, which may be interior (a singularity, or a launch point).
struct Walker {
    int tri = -1;          // triangle to cross next
    Point pos{0.0, 0.0};
    int entryEdge = -1;    // local edge of `tri` holding pos; -1 when pos is interior to it
    int atVertex = -1;     // global vertex id when pos is a mesh vertex, else -1
    double dir = 0.0;      // the field at pos, on the branch the walk is following
    // The direction that actually got the walk across the last edge. It differs
    // from `dir` by the amount the field turns over one triangle, and where the
    // streamline grazes an edge that difference decides the sign of something:
    // `dir` can point back out of the triangle just entered while this, being
    // what crossed the edge, cannot. It is the one direction always safe to
    // step along, and so what the entry into a triangle falls back on.
    double crossDir = 0.0;
};

class FieldTracer {
public:
    enum class Status {
        Ok,        // crossed a triangle, the walker is on its far side
        Boundary,  // left the mesh
        Cut,       // crossed a separatrix of a singularity orthogonally and stopped there
        Stuck      // no way forward; should not happen on a sane mesh, reported rather than hidden
    };

    // Set when a step ends in Status::Cut.
    struct CutInfo {
        int singularity = -1;  // index into getSingularities()
        int port = -1;         // which of its separatrices was crossed
    };

    // `useExactCenters` solves for the point inside the singular triangle where
    // the interpolated representation vector vanishes, instead of taking the
    // barycentre. It is the better centre when the solve lands inside the
    // triangle and nonsense when it does not, so it falls back on its own.
    //
    // `boundaryTriangles` false skips every triangle with a vertex on the
    // boundary when looking for singularities, which is what
    // CrossField::computeSingularities does; see buildSingularities() for why
    // the default is not to.
    //
    // `absorbAt`, one flag per vertex, names vertices whose incident singular
    // triangles are not singularities of the layout: the nodes of a material
    // interface network, which the layout turns at by their sectors
    // (SeparatrixTrace) and where a field that could not be aligned to every
    // incident interface turns instead. Their index is kept in
    // absorbedIndexSum4() and they emit nothing.
    explicit FieldTracer(std::shared_ptr<CrossField> cf, bool useExactCenters = true,
                         bool boundaryTriangles = true, std::vector<char> absorbAt = {});

    const Mesh &getMesh() const { return *mesh; }
    std::shared_ptr<CrossField> getCrossField() const { return crossField; }

    const std::vector<Singularity> &getSingularities() const { return singularities; }
    std::vector<Singularity> &getSingularities() { return singularities; }

    // -1 when the triangle carries no singularity.
    int singularityOfTriangle(int f) const { return singularityOf[f]; }

    // Four times the sum of the indices of every singular triangle found,
    // i.e. sum d, dropped ones included -- the interior half of the
    // Poincare-Hopf check SeparatrixTrace makes (docs/viertel_2019.md Sec. 3.2).
    int interiorIndexSum4() const { return indexSum4; }
    // Triangles of index |d| >= 2, which the paper does not treat.
    int multipleSingularityCount() const { return multipleSingularities; }
    // Triangles whose index leaves the local model fewer than two sectors or
    // more than eight (d >= 3, d <= -5): no ports, so they emit nothing.
    int droppedSingularityCount() const { return droppedSingularities; }
    // Singular triangles at an `absorbAt` vertex, and four times their index.
    int absorbedSingularityCount() const { return absorbedSingularities; }
    int absorbedIndexSum4() const { return absorbedIndex4; }
    // Per vertex: four times the index absorbed there -- the field's own
    // winding round an `absorbAt` vertex, where it is not a singularity to
    // trace from but a count the node's sectors have to agree with. A
    // triangle touching two such vertices is charged to the first.
    int absorbedIndexAt(int v) const {
        return (v >= 0 && v < static_cast<int>(absorbedAt.size())) ? absorbedAt[v] : 0;
    }

    double averageEdgeLength() const { return avgEdge; }

    // Cross-field angle at a vertex, on the branch nearest `ref`.
    double vertexAngle(int v, double ref) const;

    // Cross-field angle at the point (1-t) * a + t * b of mesh edge e = (a, b),
    // on the branch nearest `ref`. This is the canonical reading of the field
    // on an edge: it uses the two endpoint angles and nothing else, so both
    // triangles sharing the edge agree on it exactly.
    double edgeAngle(int e, double t, double ref) const;

    // Advance across one triangle, appending the points traversed to `path`
    // (one point for a regular triangle, several for a singular one). The
    // walker is left on the far side, ready for the next call.
    Status advance(Walker &w, std::vector<TracePoint> &path, CutInfo *cut = nullptr) const;

    // Place a walker at `p` inside triangle `f` heading in direction `dir`,
    // for launching a streamline from a singularity or a boundary corner.
    Walker startAt(int f, const Point &p, double dir) const;

    // Barycentric coordinates of p in triangle f.
    std::array<double, 3> barycentric(int f, const Point &p) const;

    // Where a ray leaves triangle f, and through which of its local edges.
    // `excludeEdge` is skipped, which is how the edge just entered through is
    // kept from being picked up again at distance ~0.
    bool exitOfRay(int f, const Point &from, const Point &dir, int excludeEdge,
                   Point &hit, int &edge, double &along) const;

    // How many samples per sector the hyperbolic sweep inside a singular
    // triangle uses. The samples sit on a fan of rays shared by every
    // streamline through that triangle, so raising this refines the curve
    // without ever letting two streamlines cross tangentially.
    int sweepSamplesPerSector = 24;

private:
    // The three vertex angles of a non-singular triangle, lifted to agree
    // pairwise to within pi/4 and then shifted as a block so that the
    // interpolant at `at` is the branch nearest `ref`.
    void triangleAngles(int f, const std::array<double, 3> &at, double ref,
                        std::array<double, 3> &theta) const;

    // The step inside a singular triangle: the hyperbola of Prop. 1.
    Status sweepSingular(Walker &w, int sing, std::vector<TracePoint> &path, CutInfo *cut) const;

    // Move a walker sitting exactly on a mesh vertex into whichever incident
    // triangle its direction points into.
    Status leaveVertex(Walker &w) const;

    // Hand a walker that has just reached local edge `edge` of triangle `f` at
    // parameter `along` over to the triangle on the far side, or report the
    // boundary. `outTheta` is the field there, `chordDir` the direction the
    // step actually travelled in.
    Status crossEdge(Walker &w, int f, int edge, double along, const Point &hit,
                     double outTheta, double chordDir) const;

    void buildSingularities(bool useExactCenters, bool boundaryTriangles,
                            const std::vector<char> &absorbAt);

    std::shared_ptr<CrossField> crossField;
    const Mesh *mesh = nullptr;

    std::vector<double> vertexTheta;   // arg(u_v)/4, the principal representative
    std::vector<Singularity> singularities;
    std::vector<int> singularityOf;    // triangle -> index into singularities, or -1
    int indexSum4 = 0;
    int multipleSingularities = 0;
    int droppedSingularities = 0;
    int absorbedSingularities = 0;
    int absorbedIndex4 = 0;
    std::vector<int> absorbedAt;

    double avgEdge = 1.0;
    double vertexSnap = 1e-9;          // distance below which a point is taken to be a vertex
};

#endif // __FIELD_TRACER_HXX__
