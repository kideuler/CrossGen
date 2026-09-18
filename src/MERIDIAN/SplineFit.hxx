#ifndef __SPLINE_FIT_HXX__
#define __SPLINE_FIT_HXX__

#include <array>
#include <string>
#include <vector>

#include "MERIDIAN/Arrangement.hxx"
#include "geom/BSpline.hxx"
#include "geom/Polyline.hxx"
#include "geom/Topology.hxx"
#include "mesh/Mesh.hxx"

// Stage 9 of Shepherd, Gu and Hughes (2022): the spline reconstruction.
// (docs/shepherd2022.pdf, Sec. 5; docs/ricci_flow_pipeline.md Sec. 11.)
//
// The geometry -- curves, surfaces, the least-squares fit, the Coons blend and
// the B-rep they are assembled into -- lives in src/geom, on the OpenCASCADE
// kernel, and has no idea what an arrangement is. What stays here is the part
// that is Stage 9's own: which arcs are fitted and which carried exactly, one
// knot vector for all of them, the orientation of four arcs into a patch frame,
// which vertex and edge each node and arc becomes, and the checks.
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
// ### Which arcs are fitted, and which are carried exactly
//
// Not every arc of the arrangement is the same kind of object, and only one
// kind is the fit's to approximate.
//
// A separatrix is a curve *the pipeline chose*. Nothing outside this code knows
// where it is, it is smooth by construction, and replacing it with three cubics
// is a modelling decision the fit is entitled to make.
//
// A boundary arc and an interface arc are curves *the input gave*. A boundary
// arc is a run of edges of dS and an interface arc a run of edges of the
// material network, so each one is already an exact polyline through mesh
// vertices -- the true geometry, in the only form the mesh has it. A
// least-squares cubic through such a run does not approximate a smooth curve;
// it rounds off the features the model was drawn to have. On the ICF hohlraum
// the outer wall runs straight for a while and then turns into a fillet inside
// one arc, and three cubics cut the corner by 1.1e-2 of the model -- a fifth of
// an element -- while the separatrices of the same model are fitted to 6.4e-4.
// The interfaces are worse: 2.2e-2, which is what put elements astride a
// material boundary.
//
// So by default only separatrices are fitted. A boundary or interface arc is
// carried as its polyline (Curve::exact), evaluate() walks that polyline, and
// the least-squares control points are still computed -- but they are now only
// the net the Coons blend needs, not the curve. Options::fitBoundaryArcs and
// Options::fitInterfaceArcs put either kind back under the fit, which is what
// tells a fitting artefact from a meshing one.
//
// A patch with an exact side is no longer a pure bicubic. Its net would leave
// the boundary rows on the fitted curves, so the patch is instead the Coons
// surface of its four sides *as they are evaluated* -- the polyline where the
// side is exact, the cubic where it is fitted:
//
//     S(s,t) = Net(s,t) + (1-t) d_0(s) + t d_2(s) + (1-s) d_3(t) + s d_1(t)
//
// with d_k = (exact side k) - (its control-net curve). That is the same thing
// written two ways: the Coons blend is linear in the sides, the net is the blend
// of the fitted curves (geom/Coons.hxx), and the bilinear corner term cancels
// because the fit pins both end control points to the arc's end nodes
// (geom::FitOptions::pinEnds), so every d_k vanishes at both of its ends. It is
// built as one tensor-product B-spline, geom::coonsSurface() over the side
// curves with the polylines raised to degree 3 and every pair of opposite sides
// given one knot vector -- so the patch is a single kernel surface that its own
// edges lie on, rather than a net with a correction added at evaluation time.
// The result is exact on the boundary of the patch, C0 across it -- two patches
// sharing an exact arc both land on the same polyline -- and a smooth blend
// inside. What it costs is that Patch::surface then carries the polyline's
// knots and no longer equals Patch::net. Stage 10 samples through evaluate(),
// so this is invisible to the mesh; a consumer that wants the net alone should
// turn the two options on and accept the deviation.
//
// ### The B-rep
//
// The rule above, fit each arc once and share it, is also what a boundary
// representation is made of, and the output is built as one: a geom::Vertex
// per node of the arrangement, a geom::Edge per arc on the curve evaluate()
// walks, and a geom::Face per patch on its surface, bounded by the edges of its
// four arcs -- the same edge objects, not copies. Two patches are then joined
// along an arc because they hold one edge, and Report::brepSharedEdges counts
// it once. shape() is that model, and writeSTEP() puts it where any CAD tool or
// mesher can read it.
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
        geom::BSplineCurve<2> spline;
        // The arc as an edge of the B-rep: on `spline`, or on `poly`'s curve
        // when the arc is exact, between the vertices of its two nodes. Null if
        // the kernel refused it; Report::messages says why.
        geom::Edge edge;

        // Set on an arc that is carried as its traced polyline rather than
        // approximated: evaluate() walks `poly`, and `spline` is then only the
        // control net the Coons blend is built from. See the header.
        bool exact = false;
        geom::Polyline<2> poly;      // the traced polyline, when exact

        // Of the traced polyline from the control-net curve. On a fitted arc
        // that is the geometric error of the reconstruction. On an exact arc
        // the curve *is* the polyline and this is instead how far a patch on
        // that side departs from its net there.
        double maxDeviation = 0.0;
        double rmsDeviation = 0.0;
        double length = 0.0;         // of the traced polyline
        int samples = 0;
        bool underdetermined = false; // fewer data points than control points
    };

    // One patch: its (n x n) control net with n = segments + 3, in row-major
    // order (net.control(i, j) is the control point at (s_i, t_j)), the surface
    // evaluate() reads, and the face of the B-rep on that surface.
    struct Patch {
        int face = -1;                    // the face of the arrangement
        std::array<int, 4> side{{-1, -1, -1, -1}};      // its four arcs
        std::array<bool, 4> forward{{true, true, true, true}};
        std::array<int, 4> corner{{-1, -1, -1, -1}};    // its four nodes
        // Which sides are carried exactly, and so on which of them the patch
        // departs from its net.
        std::array<bool, 4> exactSide{{false, false, false, false}};
        bool corrected = false;                         // any of them
        // The Coons net of the four fitted curves over the shared knot vector:
        // what two patches sharing an arc must agree on to the bit.
        geom::BSplineSurface<2> net;
        // The patch: the Coons surface of the four sides as evaluated. The net
        // itself unless a side is exact; see the header.
        geom::BSplineSurface<2> surface;
        // `surface` bounded by the edges of the four arcs: the face of the
        // B-rep, as `face` is the face of the arrangement. Null if the kernel
        // refused it; Report::messages says why.
        geom::Face brepFace;

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

        // Put the arcs the *input* gave -- the runs of edges along dS and
        // along the material interface network -- back under the least-squares
        // fit, instead of carrying them exactly. Off by default; the header
        // says why, and the numbers there are what turning either one on costs.
        bool fitBoundaryArcs = false;
        bool fitInterfaceArcs = false;

        // Fit and report the arcs even where Stage 8 could not close the
        // layout. The curves are useful on their own -- they are the layout
        // drawn on the model -- and refusing to produce them because some face
        // came out with three corners hides more than it protects.
        bool fitAllArcs = true;
    };

    struct Report {
        int arcs = 0;
        int curves = 0;              // arcs that came out with a curve at all
        int fittedArcs = 0;          // ... by least squares
        int exactArcs = 0;           // ... by carrying the polyline
        int correctedPatches = 0;    // patches with at least one exact side
        int underdetermined = 0;
        int controlPointsPerArc = 0;

        int faces = 0;              // patches the arrangement offered
        int patches = 0;            // ... that were fitted
        int skipped = 0;            // ... that were not, for want of four sides
        int foldedPatches = 0;      // a sampled cell came out reversed

        // How far the fitted curve is from the polyline it was fitted to, over
        // the arcs that were *fitted*, as a fraction of the diagonal of S. This
        // is the geometric error of the reconstruction and the number Sec. 5's
        // "least-squares fit a cubic B-spline" is answerable for. An exact arc
        // contributes nothing to it because it has no such error.
        double maxDeviation = 0.0;
        double rmsDeviation = 0.0;
        int worstArc = -1;

        // The same measurement on the exact arcs, where it is not an error of
        // the model but how far a patch built on the polyline departs from the
        // net's boundary row. It says how far the control net alone would have
        // been from the input.
        double maxNetDeviation = 0.0;
        int worstNetArc = -1;

        // The watertightness rule, measured rather than assumed: the largest
        // distance between the control points two patches carry for the arc
        // they share. It is zero, exactly, or the rule was broken somewhere.
        double maxSeamGap = 0.0;
        int sharedArcs = 0;

        // Corner interpolation: how far a patch's corner control point is from
        // the node of the arrangement it belongs to.
        double maxCornerGap = 0.0;

        // The same rule read off the surface rather than off the net, which is
        // what a patch with an exact side makes a separate question: how far
        // evaluate(Patch) along a side is from evaluate(Curve) on the arc that
        // side is. It says every patch sits on its own arcs, and therefore that
        // two patches sharing one meet along it -- the sampled statement of
        // maxSeamGap, on the geometry Stage 10 actually reads.
        //
        // Unlike maxSeamGap this one is *not* exactly zero, and cannot be. A
        // side traversed backwards is evaluated by reversing its control points,
        // which reproduces the reversed curve only because the uniform knot
        // vector is symmetric -- and it is not, in binary: 1 - (1/3) and (2/3)
        // are a double apart. So the mirrored basis functions differ in the last
        // bit and the two evaluations of one arc land an ulp apart; on an exact
        // side the polyline's degree elevation adds a few more. It is below
        // 4e-15 of the model on the corpus, and boundaryTolerance is what
        // separates that floor from a real gap.
        double maxBoundaryGap = 0.0;
        static constexpr double boundaryTolerance = 1e-12;

        double patchArea = 0.0;
        double faceArea = 0.0;      // the same patches as arrangement faces
        double minCellRatio = 0.0;

        // The B-rep, counted by the kernel. Every patch should have a face and
        // every arc two patches share should be one edge bounding both of them,
        // so brepFaces == patches and brepSharedEdges == sharedArcs; the free
        // edges are then the boundary of the patched region. brepValid is the
        // kernel's own check of the whole shape (BRepCheck).
        int brepFaces = 0;
        int brepEdges = 0;
        int brepSharedEdges = 0;
        int brepFreeEdges = 0;
        int brepFailures = 0;       // arcs and patches the kernel refused
        bool brepValid = false;

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

    // The B-rep: every face, and the edges of arcs no patch uses, in the plane
    // z = 0.
    const geom::Shape& shape() const { return model; }
    bool writeBREP(const std::string &filename) const { return model.writeBREP(filename); }
    bool writeSTEP(const std::string &filename) const { return model.writeSTEP(filename); }

    // The clamped uniform knot vector these curves are written over, and the
    // Greville abscissae the Coons blend is taken at. Public because Stage 10's
    // knot insertion needs them.
    const std::vector<double>& knots() const { return knot; }
    const std::vector<double>& grevilles() const { return greville; }

private:
    void fitArcs();
    void buildPatches();
    void buildShape();
    void check();

    const Arrangement *arr = nullptr;
    Options options;

    std::vector<double> knot;
    std::vector<double> greville;
    std::vector<Curve> fitted;     // one per arc of the arrangement, in order
    std::vector<Patch> nets;
    geom::Shape model;

    double modelExtent = 1.0;
    Report report;
};

#endif // __SPLINE_FIT_HXX__
