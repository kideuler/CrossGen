#ifndef __SPLINE_FIT_HXX__
#define __SPLINE_FIT_HXX__

#include <array>
#include <string>
#include <vector>

#include "MERIDIAN/Arrangement.hxx"
#include "mesh/Mesh.hxx"

// Stage 9 of Shepherd, Gu and Hughes (2022): the spline reconstruction.
// (docs/shepherd2022.pdf, Sec. 5; docs/ricci_flow_pipeline.md Sec. 11.)
//
// Stage 8 left a partition of S into quadrilaterals whose sides are arcs. This
// stage turns each arc into a cubic B-spline and each quadrilateral into a
// bicubic patch, and the whole point of it is one rule:
//
//     **fit each arc exactly once, and give the resulting control points to
//     both patches that share it.**
//
// That is the mechanism -- the only mechanism -- that makes the output
// watertight. Fitting an arc twice, once per patch, would put two curves within
// a fitting tolerance of each other where the model wants one, which is the gap
// problem the entire method exists to remove: it is what a trimmed B-Rep does,
// and Sec. 1 is a list of what it costs downstream. So the arcs are fitted
// first, patch by patch never enters into it, and the boundary rows of a
// patch's control net are *copied* from the arc fits rather than recomputed.
//
// ### Why every arc gets the same number of segments
//
// A Coons patch blends four boundary curves, and blending them at the level of
// the control points -- rather than sampling the surface and refitting --
// requires opposite sides to share a knot vector. Choosing the knot vector per
// arc, from its own length or its own curvature, would break the neighbour's
// fit, because the neighbour needs *its* opposite side to match. The resolution
// is the one the paper's uniform patch counts imply: a globally constant number
// of cubic segments per arc. Sec. 5 uses 7x7 control points (four segments) for
// the shock house and the chassis and 6x6 (three) for the firewall and the
// speaker frame; Options::segments is that number and nothing else varies.
//
// ### The Coons construction, at the control-point level
//
//     S(s,t) = (1-t) C0(s) + t C1(s) + (1-s) D0(t) + s D1(t)
//              - [ (1-s)(1-t) P00 + s(1-t) P10 + (1-s)t P01 + s t P11 ]
//
// applied to the control points at their Greville abscissae rather than to the
// curves at their parameters. The Grevilles are what make the bilinear
// correction term reproduce a linear function exactly, and using i/(n-1)
// instead would leave a small bias along the sides. The four boundary rows come
// out equal to the four side fits identically -- not to a tolerance -- which is
// the watertightness rule above, checked rather than assumed in
// Report::maxSeamGap.
//
// ### The pullback
//
// Sec. 5 begins by pulling each traced arc back to R^3 through its barycentric
// coordinates. Here S is a planar domain -- Mesh carries two coordinates per
// vertex -- so that step is the identity on the polyline Stage 8 already holds,
// and there is nothing to do. It is still the same construction: the arcs keep
// the (face, barycentric entry, barycentric exit) steps they were cut from, so
// a mesh that one day carries positions in R^3 is fitted by the same code with
// the same numbers, and the fit survives a refinement of the triangulation that
// keeps the geometry.
//
// ### What this does not do
//
// Bilinear blending of the four sides cannot represent curvature in the
// interior of a patch that the sides do not already carry; Sec. 5's own caveat.
// On the planar domains this codebase works with that is exact rather than
// merely acceptable, since a Coons patch on four planar curves is planar. On a
// curved shell it is the known limitation and the answer is more patches or a
// manifold spline space.
//
// Continuity is C0 between patches -- they share a boundary curve exactly --
// and C2 inside a patch, which is what the single-multiplicity interior knots
// of a uniform cubic give. Stage 10's refinement preserves both and is not
// implemented here.
class SplineFit {
public:
    // A cubic B-spline over a clamped uniform knot vector with
    // Options::segments spans, so segments + 3 control points, the first and
    // last of which are the arc's two end nodes exactly.
    struct Curve {
        int arc = -1;
        std::vector<Point> ctrl;
        double maxDeviation = 0.0;   // of the traced polyline from the fit
        double rmsDeviation = 0.0;
        double length = 0.0;         // of the traced polyline
        int samples = 0;
        bool underdetermined = false; // fewer data points than control points
    };

    // One bicubic patch, as an (n x n) control net with n = segments + 3, in
    // row-major order: net[j * n + i] is the control point at (s_i, t_j).
    struct Patch {
        int face = -1;                    // the face of the arrangement
        std::array<int, 4> side{{-1, -1, -1, -1}};      // its four arcs
        std::array<bool, 4> forward{{true, true, true, true}};
        std::array<int, 4> corner{{-1, -1, -1, -1}};    // its four nodes
        std::vector<Point> net;

        double area = 0.0;         // of the surface, from a sampled grid
        double faceArea = 0.0;     // of the arrangement face it was fitted to
        double minCellRatio = 0.0; // smallest sampled cell area over the mean;
                                   // <= 0 is a folded patch
    };

    struct Options {
        // Cubic Bezier segments per arc, so segments + 3 control points per
        // arc and (segments + 3)^2 per patch. The paper's two values are 3 and
        // 4 -- 6x6 and 7x7.
        int segments = 3;

        // Tikhonov weight, relative to the largest diagonal of the normal
        // equations, pulling an underdetermined fit towards the straight line
        // between the arc's ends. It binds only where an arc has fewer sample
        // points than it has free control points, which happens on an arc a
        // handful of triangles long; elsewhere it is below the noise of the
        // fit. Without it those arcs come back with control points anywhere at
        // all, and the patch built on them folds.
        double regularisation = 1e-6;

        // Grid resolution used to sample a patch for its area and its worst
        // cell, and by writeSurfaceOBJ() when it is not given one.
        int samples = 8;

        // Fit and report the arcs even where Stage 8 could not close the
        // layout. The curves are useful on their own -- they are the layout
        // drawn on the model -- and refusing to produce them because some face
        // came out with three corners hides more than it protects.
        bool fitAllArcs = true;
    };

    struct Report {
        int arcs = 0;
        int curves = 0;
        int underdetermined = 0;
        int controlPointsPerArc = 0;

        int faces = 0;              // patches the arrangement offered
        int patches = 0;            // ... that were fitted
        int skipped = 0;            // ... that were not, for want of four sides
        int foldedPatches = 0;      // a sampled cell came out reversed

        // How far the fitted curve is from the polyline it was fitted to, over
        // every arc, as a fraction of the diagonal of S. This is the geometric
        // error of the reconstruction and the number Sec. 5's "least-squares
        // fit a cubic B-spline" is answerable for.
        double maxDeviation = 0.0;
        double rmsDeviation = 0.0;
        int worstArc = -1;

        // The watertightness rule, measured rather than assumed: the largest
        // distance between the control points two patches carry for the arc
        // they share. It is zero, exactly, or the rule was broken somewhere.
        double maxSeamGap = 0.0;
        int sharedArcs = 0;

        // Corner interpolation: how far a patch's corner control point is from
        // the node of the arrangement it belongs to.
        double maxCornerGap = 0.0;

        double patchArea = 0.0;
        double faceArea = 0.0;      // the same patches as arrangement faces
        double minCellRatio = 0.0;

        bool watertight = false;
        bool valid = false;

        std::vector<std::string> messages;
    };

    explicit SplineFit(const Arrangement &arrangement);
    SplineFit(const Arrangement &arrangement, const Options &opts);

    const std::vector<Curve>& curves() const { return fitted; }
    const std::vector<Patch>& patches() const { return nets; }
    const Report& getReport() const { return report; }
    const Options& getOptions() const { return options; }
    const Arrangement& getArrangement() const { return *arr; }

    int controlPointsPerArc() const { return options.segments + 3; }

    // A point of one of the fitted arcs, at parameter u in [0, 1].
    Point evaluate(const Curve &c, double u) const;
    // A point of a patch, at (s, t) in [0, 1]^2.
    Point evaluate(const Patch &p, double s, double t) const;

    // The fitted arcs, sampled, as .obj polylines.
    bool writeCurvesOBJ(const std::string &filename, int samples = 32) const;
    // The control nets, as a grid of lines: the picture of Fig. 19.
    bool writeNetOBJ(const std::string &filename) const;
    // The patches themselves, tessellated into a quad mesh -- the reconstructed
    // model, one .obj group per patch.
    bool writeSurfaceOBJ(const std::string &filename, int samples = 0) const;

    // The clamped uniform knot vector these curves are written over, and the
    // Greville abscissae the Coons blend is taken at. Public because Stage 10's
    // knot insertion needs them.
    const std::vector<double>& knots() const { return knot; }
    const std::vector<double>& grevilles() const { return greville; }

private:
    void buildKnots();
    void fitArcs();
    void buildPatches();
    void check();

    Curve fitOne(const std::vector<Point> &poly) const;
    static void basisFuns(int span, double u, const std::vector<double> &U, double *N);
    int findSpan(double u) const;

    const Arrangement *arr = nullptr;
    Options options;

    std::vector<double> knot;
    std::vector<double> greville;
    std::vector<Curve> fitted;     // one per arc of the arrangement, in order
    std::vector<Patch> nets;

    double modelExtent = 1.0;
    Report report;
};

#endif // __SPLINE_FIT_HXX__
