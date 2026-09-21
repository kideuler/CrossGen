#ifndef __COARSE_DOMAIN_HXX__
#define __COARSE_DOMAIN_HXX__

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "ATLAS/PlanarDomain.hxx"
#include "ATLAS/SquareCarrier.hxx"
#include "mesh/Mesh.hxx"

// A coarse re-triangulation of a validated domain, on which ATLAS searches for
// the layout's topology before realising it on the input (Realisation).
//
// ### Why the carrier's resolution is the problem
//
// Sec. 8.1 says aggregation cannot remove the singularities the carrier
// inherits from its triangulation; Stage 6 is there to rewrite them. What it
// leaves out is how the difficulty scales. Every triangle of the input puts a
// valence-3 vertex into the carrier and every input vertex one of valence
// about six, so the carrier of an 18k-triangle mesh starts with ~10^4
// singularities, and Stage 6's local cavities remove the dense ones and stall
// on the isolated 3-5 pairs that remain (~500 on singlemat/geom007, whose
// cross field needs nine). An isolated 3-5 pair is a dislocation: the
// development of any loop round it closes in rotation but not in translation,
// so no fill that keeps the loop can remove it. It has to glide to dS or meet
// a partner of opposite Burgers vector, and on a fine carrier both are dozens
// of cells away. The base complex that is left is thousands of blocks.
//
// None of that is about the domain. Sec. 13.3's last regression row asks for
// exactly this separation -- "same boundary, different triangulations: measure
// whether successful coarse replacements remove dependence on interior
// triangle density" -- and the direct way to get it is to search on a carrier
// whose resolution is the domain's, not the mesher's: a few hundred cells,
// where every singularity is a few cells from dS or from its partner, every
// cavity is a sizeable fraction of the domain, and the whole of Stages 3-6 is
// cheap enough to run to convergence. Measured on geom007 before anything
// else changed: the unchanged pipeline on a 96-triangle mesh of the same
// domain ends at 19 irregular vertices and 93 blocks, against 502 and 10956
// on the 7032-triangle input.
//
// ### What is kept exactly, and what is a proxy
//
// The coarse mesh's feature vertices are input feature vertices, sampled along
// each chain between protected corners, so every protected corner is one and
// the coarse feature network is a polygon inscribed in the input's. That
// polygon is only a proxy: the layout found on it is carried back to the
// input's own curves by Realisation, and certified there. Two things are
// transferred so the proxy does not mislead the search:
//
//   * the *layout* angle at a sample is the input's angle there
//     (PlanarDomain::Options::angleOverride, and interfaceAngleOverride along
//     an interface), so a circle sampled eight times is still a smooth curve
//     with valence target 2 everywhere, not an octagon with eight corners;
//   * every coarse feature vertex has a location on the input's curve
//     (locate()): a sample is its own input vertex, and a point Stages 2-6 put
//     on a coarse feature edge sits at the same fraction of the input arc
//     between the edge's two samples.
//
// ### Features: dS and the interface network
//
// Sec. 1.1 asks the same two things of a material interface as of dS -- every
// segment of it survives, and it is a union of macroedges -- so an interface
// is sampled, constrained and carried back exactly as a boundary chain is, and
// the only structural difference is that it has a material on both sides. So
// the curves below ("arcs") are the boundary loops *and* the chains of the
// interface network between its nodes (junctions, landings on dS, kinks,
// dangling ends -- every interface vertex Stage 1 protects), a chain with no
// node being closed like a loop.
//
// Without this, a multi-material domain had no coarse search at all and fell
// back on the fine carrier, which is the case Stage 6 cannot finish: multimat/
// geom001, a quarter disk in a square, ended at 29054 blocks against the 8 its
// layout wants.
//
// The materials of the coarse triangles are *not* read off the input by point
// location, which is wrong in the sliver between a sampled chord and the arc
// it cuts. They are flooded out from the chains themselves: the coarse
// triangle on the left of a sampled interface segment takes the material on
// the left of that chain in the input, the one on its right the other, and the
// flood stops at the constrained segments. Every coarse triangle on one side
// of the network is then of one material by construction, whatever the
// sampling did to the geometry.
//
// ### Sampling
//
// Each chain is sampled at a spacing that is the least of three bounds, graded
// along the loop so it varies by at most Options::grading per unit length:
//
//   * a global one, Options::maxSpacing of the bounding-box diagonal;
//   * the gap: gapFraction of the distance to the nearest feature point that
//     is not a neighbour along the same curve, so a narrow passage between two
//     holes, or between a hole and dS, or between an interface and dS, is at
//     least two coarse triangles wide;
//   * the curvature: curvatureFraction of the local radius, so a chord never
//     strays far from its arc (0.8 puts eight chords on a circle, whose
//     sagitta is then 7.6% of the radius).
//
// Triangle then fills the interior with -YY (no point on any segment, boundary
// or interface: those are all input vertices) and a quality bound, which
// grades the interior to the feature spacing on its own.
//
// ### Corners
//
// In the three-quad split a boundary vertex's valence is its number of
// triangles, and a quality triangulation usually cuts a right angle into two.
// A corner with two cells sends a separatrix across the domain on the
// diagonal, and Stage 6 undoes it badly, since a corner's fan is only ever in
// a cavity whole. So edges are flipped at every boundary sample that has more
// triangles than its layout angle asks for, as long as that lowers the total
// excess on dS and no angle falls below Options::minFlipAngle. Measured on
// singlemat/geom025 (a rounded plate with a V notch): 12 blocks -> 6.
class CoarseDomain {
public:
    struct Options {
        double maxSpacing = 0.125;       // fraction of the bounding-box diagonal
        double gapFraction = 0.6;
        double curvatureFraction = 0.8;
        double grading = 0.4;
        double minAngle = 20.0;          // Triangle's quality bound, degrees
        // Minimum samples on a closed chain (a loop with no protected corner).
        int minClosedSamples = 4;
        // Edge flips that give each boundary sample the number of triangles
        // its layout angle asks for may not make an angle smaller than this.
        double minFlipAngle = 12.0;
    };

    // A point on a fine feature curve: which arc, and the arc-length position
    // along it. `loop` is an index into arcs().
    struct Location {
        int loop = -1;
        double s = 0.0;
    };

    // One fine feature curve, parameterised by arc length: a boundary loop
    // (closed, the domain on the left) or a chain of the interface network.
    // A closed arc has one edge per vertex, the last running back to the
    // first; an open one has n - 1 edges and s.back() == length.
    struct Arc {
        std::vector<int> vertices;       // fine vertex ids, the domain on the left
        std::vector<int> edges;          // fine edge i runs vertices[i] -> vertices[i+1]
        std::vector<double> s;           // arc length at vertices[i]
        double length = 0.0;
        bool closed = true;
        bool onInterface = false;        // a chain of the interface network, not dS
        // The input materials either side of an interface chain, walking it
        // from vertices[0]: what the flood fill seeds the coarse triangles on
        // each side with. Both 0 on a boundary loop.
        int leftMaterial = 0, rightMaterial = 0;
    };

    struct Report {
        bool valid = false;
        std::string reason;
        int chains = 0;
        int samples = 0;                 // feature vertices of the coarse mesh
        int vertices = 0, triangles = 0;
        int flips = 0;                   // boundary-valence edge flips
        int interfaceArcs = 0;           // chains of the interface network
        int interfaceSegments = 0;       // coarse interface edges they became
        int materialRegions = 0;         // regions the flood fill separated
        double spacingMin = 0.0, spacingMax = 0.0;
        double seconds = 0.0;
    };

    CoarseDomain(const PlanarDomain &fine, const Options &opts);

    bool valid() const { return report_.valid; }
    const Report &getReport() const { return report_; }
    const Options &getOptions() const { return opts_; }
    const PlanarDomain &getFine() const { return fine_; }
    const Mesh &getMesh() const { return *mesh_; }
    const PlanarDomain &getDomain() const { return *domain_; }

    // Coarse mesh vertex -> the fine vertex it samples, -1 for an interior
    // Steiner point.
    const std::vector<int> &fineVertexOf() const { return fineVertex_; }

    const std::vector<Arc> &arcs() const { return arcs_; }
    // A fine vertex's index in one arc, -1 if it is not on that arc.
    int arcIndex(int arc, int fineVertex) const;
    // Its arc-length position there, or -1.
    double positionOf(int arc, int fineVertex) const;
    // The arc a coarse feature edge was sampled from, -1 if it is not one.
    int edgeArc(int coarseEdge) const {
        return coarseEdge >= 0 && coarseEdge < static_cast<int>(edgeArc_.size()) ? edgeArc_[coarseEdge] : -1;
    }

    // Where a feature vertex of a carrier over getDomain() lies on the fine
    // curves. A vertex Stages 2-6 put inside a coarse feature edge is on one
    // arc; a sample is on every arc through its input vertex -- a landing is
    // on dS and on its interface chain, a junction on each of its branches --
    // so all of them are returned, and locate() gives the first.
    void locateAll(const SquareCarrier &C, int v, std::vector<Location> &out) const;
    // The first of them; loop < 0 when the vertex is on no feature curve.
    Location locate(const SquareCarrier &C, int v) const;

    // The smoothest boundary-aligned cross field, as a reference direction
    // for scoring layouts (CavityRewrite::AnnealOptions::field), never
    // integrated or traced. It is the harmonic interpolation, over the coarse
    // triangulation, of the representation vector u = exp(4 i theta) of the
    // input boundary's tangent -- Dirichlet on every boundary sample that is
    // not a protected corner, where the tangent jumps. |u| falls towards the
    // field's singularities, which is exactly where its direction means
    // least, so it is returned as the weight. Returns theta in [0, pi/2).
    double crossAngle(const Point &p, double *weight) const;

    // Forward arc length from a to b along an arc. On a closed one it is in
    // (0, length] and a == b gives the whole loop; on an open one it is
    // simply b - a, clamped to the arc.
    double forward(int loop, double a, double b) const;
    // The point at arc position s, and the fine edge it lies on (the edge
    // starting at the vertex at or before s).
    Point pointAt(int loop, double s, int *edge = nullptr, int *index = nullptr) const;

private:
    void buildArcs();
    void buildInterfaceArcs();
    void indexArc(int arc);
    void geodesic(int source, double cap, std::vector<double> &dist, std::vector<int> &touched) const;
    std::vector<double> spacing(int loop) const;
    bool sample(std::vector<std::vector<int>> &samples);
    bool triangulate(const std::vector<std::vector<int>> &samples);
    // Materials of the coarse triangles, flooded out from the sampled chains.
    // Each segment is a consecutive pair of samples in its arc's own
    // direction, so the triangle on its left takes that arc's leftMaterial.
    bool floodMaterials(const std::vector<Point> &V, const std::vector<Triangle> &T,
                        const std::vector<std::array<int, 2>> &segments, const std::vector<int> &segArc,
                        std::vector<int> &matOut);
    void buildField();

    const PlanarDomain &fine_;
    Options opts_;
    Report report_;
    std::vector<Arc> arcs_;
    // (arc, fine vertex) -> index in that arc, as arc * numVertices + vertex.
    std::unordered_map<long long, int> arcIndex_;
    // Fine vertex -> the arcs through it, in arcs_ order.
    std::vector<std::vector<int>> arcsAt_;
    // The feature network as a weighted graph on the fine vertices: every
    // boundary and interface edge, for the geodesic the gap bound measures
    // "far along the curve" with.
    std::vector<std::vector<std::pair<int, double>>> featAdj_;
    int numBoundaryArcs_ = 0;
    std::shared_ptr<Mesh> mesh_;
    std::unique_ptr<PlanarDomain> domain_;
    std::vector<int> fineVertex_;
    // Coarse feature edge -> its start vertex walking its arc forward, and
    // the arc it was sampled from.
    std::vector<int> edgeFrom_;
    std::vector<int> edgeArc_;
    // The field on a regular lookup grid over the bounding box: u per node
    // (0 outside the domain), bilinear in between.
    std::vector<std::array<double, 2>> fieldGrid_;
    int fieldNx_ = 0, fieldNy_ = 0;
    Point fieldLo_{0.0, 0.0};
    double fieldStep_ = 1.0;
    double fieldMax_ = 1.0;
};

#endif // __COARSE_DOMAIN_HXX__
