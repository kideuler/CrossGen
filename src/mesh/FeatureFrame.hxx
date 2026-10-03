#ifndef __FEATURE_FRAME_HXX__
#define __FEATURE_FRAME_HXX__

#include <complex>
#include <string>
#include <vector>

#include "mesh/Mesh.hxx"

// The cross a model's walls ask a quad grid to follow, at every point of the
// model: the harmonic extension of the tangent cross of dS and of the material
// interfaces -- u = exp(4i theta) on them, theta the tangent's angle, and the
// Laplace equation inside. What BlockDecomposition::alignmentQuality() and
// mesh::QuadMesh::alignmentQuality() grade a layout and a mesh against
// (docs/block_decomposition_metrics.md Sec. 7). Read off the triangle mesh
// before any method has run, as BoundaryFeatures is, and for the same reason:
// it is the model's, and every method's answer is graded against one copy.
//
// ### Why a grid should follow it, in a shock code
//
// A front that leaves a wall or an interface -- a reflected shock, a
// transmitted one, the rarefaction behind it, the slip line along an
// interface -- starts out parallel to it and travels along its normal. A
// shock is captured best by cell faces parallel or orthogonal to it: Gnoffo
// (AIAA 2009-599, Sec. II) makes that the condition for good shock capturing
// in the hypersonic regime, and the false diffusion a first-order scheme adds
// across an oblique front grows as sin 2 Delta of the angle Delta between the
// front and the grid lines (de Vahl Davis and Mallinson, Comput. Fluids 4,
// 1976). In a Lagrangian code the same misalignment is where mesh imprinting
// and spurious vorticity begin (Dukowicz and Meltz, J. Comput. Phys. 99, 1992,
// the skewed-grid Saltzman piston; Dobrev, Kolev and Rieben, SIAM J. Sci.
// Comput. 34, 2012, Secs. 1 and 6.3). The frame a grid should hold to is
// therefore the walls' own, carried inward; and where walls of different
// directions are in view at once, there is no one frame to hold to at all.
// Both halves of that are what this class computes.
//
// ### Why the harmonic extension, and what |u| means
//
// u(x) is the average of the walls' crosses weighted by their harmonic
// measure from x -- the share of x's surroundings each stretch of wall takes
// up -- so it is the wall's own frame next to a wall, the frame of the nearer
// and larger walls further in, and a compromise where they meet. |u| <= 1 is
// how far they agree: 1 on a wall and everywhere in a box (whose four walls
// share one cross), r^4 in a disk of radius 1 (a circle's tangents agree with
// one another only near the circle), and 0 at its singular points, round
// which the cross turns as it does round an irregular vertex of a layout.
// That makes |u| the weight a grade
// wants: a grid is held to the walls' direction where there is one, and not
// blamed where there is none.
//
// It is the plain linear solve, not a Ginzburg-Landau field like DualMBO, on
// purpose. It has no parameter, seed or continuation, so it cannot favour the
// method that shares its solver -- TORSION integrates DualMBO, and ATLAS's
// own ReferenceField is DualMBO too. Where the two disagree, the core of a
// round region, DualMBO pins four cones and so prefers one orientation of an
// O-grid's core; the harmonic field gives that core no weight, which for a
// front converging on it or spreading from it is the physics: no orientation
// of a square core is better aligned with a circle than another. ATLAS's
// guidance note (docs/atlas_crossfield_guidance.md Sec. 2.2) measured the two
// fields 2.9 degrees apart on the single-material corpus, median over models.
// Graded against an MBO (Ginzburg-Landau) field of the same walls instead,
// the five methods' blocks on those 35 models score with a correlation of
// 0.94, and ATLAS against the best cross-field method comes out the same way
// on 29 of the 31 models where they differ; the two that flip are the disks
// geom002 and geom003, where the MBO field pins the core.
//
// ### The solve
//
// P1 on the model's own triangles, the cotangent Laplacian, real and
// imaginary parts as two right-hand sides of one sparse Cholesky
// factorisation. A vertex on dS or on an interface is fixed to the mean of
// exp(4i theta) over its feature edges: on a smooth run that is the run's
// cross, at a right-angled corner it is still the corner's (both sides share
// one cross), and at a 135-degree corner the two sides' crosses cancel to 0,
// which is the right answer to "which way should a grid run here". An
// interface is fixed like dS, so each material region is solved against its
// own walls. Cost: a few tens of milliseconds on the corpus meshes.
class FeatureFrame {
public:
    struct Options {
        // Hold the frame to the material interfaces as well as to dS. Off,
        // the interfaces are invisible to it and a region's frame is blended
        // across them with its neighbours'.
        bool interfaces = true;
    };

    struct Report {
        bool solved = false;
        std::string reason;          // why not, when not
        int vertices = 0, triangles = 0;
        int fixedVertices = 0;       // on dS or an interface
        int boundaryEdges = 0, interfaceEdges = 0;
        // Area-weighted mean of |u| over the model: how much of it has a
        // direction to follow at all. 1 for a box, 1/3 for a disk.
        double meanCoherence = 0.0;
    };

    explicit FeatureFrame(const Mesh &model);
    FeatureFrame(const Mesh &model, const Options &opts);

    // u at p: exp(4i theta) of the walls' cross there, scaled by how far they
    // agree, so |u| <= 1. Linear on the model's triangles. A point just
    // outside them -- on the sagitta between a curved wall's chord and its
    // arc -- takes the value of the nearest triangle's nearest point.
    std::complex<double> at(const Point &p) const;

    // u per vertex of the model.
    const std::vector<std::complex<double>> &vertexValues() const { return u_; }

    const Report &getReport() const { return report_; }
    const Options &options() const { return opts_; }

    // exp(4i theta) of a direction: the cross it lies on, as a unit
    // representation vector; 0 for a zero vector.
    static std::complex<double> crossOf(const Point &d);

    // Sums cells one at a time -- the Coons cells of a block, the elements of
    // a mesh -- and grades them together, so that BlockDecomposition's and
    // mesh::QuadMesh's alignmentQuality() are one number read on two objects.
    //
    // A cell is its centre, the two directions its grid lines run in there,
    // and its area. Each direction is off the frame's cross by an angle Delta
    // in [0, 45 degrees]; the cell's misalignment is sin^2 2 Delta averaged
    // over its two directions -- 0 for a cell whose sides are parallel or
    // orthogonal to the walls' cross, 1 for one at 45 degrees to it, and
    // ((1 - cos 4 Delta) / 2, so with no angle to wrap) -- weighted by its area
    // and by |u| at its centre. A sheared cell is off in at least one of its
    // two directions, so skew counts as misalignment, which in a shock code
    // it is.
    class Alignment {
    public:
        explicit Alignment(const FeatureFrame &frame) : frame_(frame) {}
        void add(const Point &centre, const Point &du, const Point &dv, double area);
        // E, the weighted mean of sin^2 2 Delta, as the quality
        // 1 - (2 / pi) asin sqrt(E): 1 for a grid aligned everywhere the walls
        // have a direction, 1/2 for one no better aligned than a grid at a
        // random angle (E = 1/2), 0 for one at 45 degrees everywhere. It is
        // 1 - Delta_e / 45 degrees, Delta_e the single angle off that would
        // give the same E, so a difference reads in degrees. 1 when nothing
        // added had a direction to follow.
        double quality() const;
        double weight() const { return den_; }
        double meanMisalignment() const { return den_ > 0.0 ? num_ / den_ : 0.0; }

    private:
        const FeatureFrame &frame_;
        double num_ = 0.0, den_ = 0.0;
    };

private:
    Options opts_;
    Report report_;
    std::vector<Point> vertices_;
    std::vector<Triangle> triangles_;
    std::vector<std::complex<double>> u_;

    // A uniform grid of buckets over the triangles, about two to a bucket,
    // each triangle in every bucket its bounding box touches: at()'s point
    // location. Mesh::findTriangleContainingPoint walks, which can stall in a
    // model with holes.
    Point lo_{0.0, 0.0};
    double cellW_ = 1.0, cellH_ = 1.0;
    int nx_ = 0, ny_ = 0;
    std::vector<int> bucketStart_, bucketTriangles_;

    void buildBuckets();
};

#endif // __FEATURE_FRAME_HXX__
