#ifndef __COARSE_DOMAIN_HXX__
#define __COARSE_DOMAIN_HXX__

#include <memory>
#include <string>
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
// The coarse mesh's boundary vertices are input boundary vertices, sampled
// along each chain of dS between protected corners, so every protected corner
// is one and the coarse boundary is a polygon inscribed in the input's. That
// polygon is only a proxy: the layout found on it is carried back to the
// input's own boundary by Realisation, and certified there. Two things are
// transferred so the proxy does not mislead the search:
//
//   * the *layout* angle at a sample is the input's angle there
//     (PlanarDomain::Options::angleOverride), so a circle sampled eight times
//     is still a smooth curve with valence target 2 everywhere, not an octagon
//     with eight corners;
//   * every coarse boundary vertex has a location on the input boundary
//     (locate()): a sample is its own input vertex, and a point Stages 2-6 put
//     on a coarse boundary edge sits at the same fraction of the input arc
//     between the edge's two samples.
//
// ### Sampling
//
// Each chain is sampled at a spacing that is the least of three bounds, graded
// along the loop so it varies by at most Options::grading per unit length:
//
//   * a global one, Options::maxSpacing of the bounding-box diagonal;
//   * the gap: gapFraction of the distance to the nearest boundary point that
//     is not a neighbour along the same curve, so a narrow passage between two
//     holes, or between a hole and dS, is at least two coarse triangles wide;
//   * the curvature: curvatureFraction of the local radius, so a chord never
//     strays far from its arc (0.8 puts eight chords on a circle, whose
//     sagitta is then 7.6% of the radius).
//
// Triangle then fills the interior with -Y (no point on a boundary segment:
// those are all input vertices) and a quality bound, which grades the
// interior to the boundary spacing on its own.
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

    // A point on the fine boundary: its loop and arc-length position.
    struct Location {
        int loop = -1;
        double s = 0.0;
    };

    // One fine boundary loop, parameterised by arc length.
    struct Arc {
        std::vector<int> vertices;       // fine vertex ids, the domain on the left
        std::vector<int> edges;          // fine edge i runs vertices[i] -> vertices[i+1]
        std::vector<double> s;           // arc length at vertices[i]
        double length = 0.0;
    };

    struct Report {
        bool valid = false;
        std::string reason;
        int chains = 0;
        int samples = 0;                 // boundary vertices of the coarse mesh
        int vertices = 0, triangles = 0;
        int flips = 0;                   // boundary-valence edge flips
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
    // Fine boundary vertex -> its index in its Arc, -1 inside.
    int arcIndex(int fineVertex) const { return arcIndex_[fineVertex]; }
    int arcOf(int fineVertex) const { return fine_.loopOf[fineVertex]; }

    // Where a boundary vertex of a carrier over getDomain() lies on the fine
    // boundary. loop < 0 when it is not on a coarse boundary edge at all.
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

    // Forward arc length from a to b along a loop, in (0, length]; a == b
    // gives the whole loop.
    double forward(int loop, double a, double b) const;
    // The point at arc position s, and the fine edge it lies on (the edge
    // starting at the vertex at or before s).
    Point pointAt(int loop, double s, int *edge = nullptr, int *index = nullptr) const;

private:
    void buildArcs();
    std::vector<double> spacing(int loop) const;
    bool sample(std::vector<std::vector<int>> &samples);
    bool triangulate(const std::vector<std::vector<int>> &samples);
    void buildField();

    const PlanarDomain &fine_;
    Options opts_;
    Report report_;
    std::vector<Arc> arcs_;
    std::vector<int> arcIndex_;
    std::shared_ptr<Mesh> mesh_;
    std::unique_ptr<PlanarDomain> domain_;
    std::vector<int> fineVertex_;
    // Coarse boundary edge -> its start vertex in loop order.
    std::vector<int> edgeFrom_;
    // The field on a regular lookup grid over the bounding box: u per node
    // (0 outside the domain), bilinear in between.
    std::vector<std::array<double, 2>> fieldGrid_;
    int fieldNx_ = 0, fieldNy_ = 0;
    Point fieldLo_{0.0, 0.0};
    double fieldStep_ = 1.0;
    double fieldMax_ = 1.0;
};

#endif // __COARSE_DOMAIN_HXX__
