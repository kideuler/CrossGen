#ifndef __CONFORMAL_FIELD_HXX__
#define __CONFORMAL_FIELD_HXX__

#include <memory>
#include <string>
#include <vector>

#include "crossfield/CrossField.hxx"
#include "mesh/Mesh.hxx"

class Interfaces;

// ---------------------------------------------------------------------------
// The corrected cross field of docs/quad_layout_pipeline.md Sec. 6: the second,
// *linear* field solve of Jezdimirovic, Chemin, Reberol, Henrotte and Remacle
// (arXiv:2103.02939, Sec. 5 -- "P2" below), made from the singularity pattern
// of the first one.
//
// ### Why a second field
//
// The first field (Stage 1, the P1 MBO of CrossField) is a Ginzburg-Landau
// minimiser in all but name, and Sec. 4.5 says what it is good for and what it
// is not: its *pattern* -- how many singularities, of which index, roughly
// where -- is reliable, because an index is an integer and a location is wanted
// only approximately; its *directions near a singularity* are not, because a
// singularity is exactly where the penalty failed to keep |u| = 1 and the
// field in its core is whatever the arithmetic left there. FieldTracer reads
// the ports off the three corners of the singular triangle and says how badly
// they disagree (Singularity::portResidual); a separatrix that leaves along a
// port read off a core like that misses the singularity it should have reached
// and runs on, which is the tangential crossing of Sec. 7.8 and, downstream,
// most of what Sec. 4 of Viertel et al. has to collapse.
//
// So the pattern is kept and the field is thrown away. With the singularities
// S = {(X_j, k_j)} given, a field matching them is the log-derivative of a
// conformal map, f = exp(H + i theta) (Sec. 6.1), and H and theta are two
// *harmonic* functions related by Cauchy-Riemann (P2 eq. 10). Nothing is
// nonlinear any more: all of it was pushed into the known positions of the
// singularities.
//
// ### The construction (Secs. 6.2 to 6.5, with one change)
//
//   1. H solves Poisson's equation with a Dirac of weight 2 pi k_j at each X_j
//      and Neumann data the boundary curvature (P2 eqs. 11-14). In P1 weak form
//      the right-hand side at boundary vertex b is its turning angle
//      tau_b = pi - alpha_b, less 2 pi k_b at a corner of Table 1 (Sec. 6.3),
//      and it sums to zero exactly when the pattern satisfies Poincare-Hopf --
//      the solvability condition of a pure Neumann problem *is* Poincare-Hopf,
//      so assembling it checks the pattern for free.
//
//   2. theta is the harmonic conjugate of H. With P1 H and Crouzeix-Raviart
//      theta the discrete Cauchy-Riemann equation is exact and theta is had by
//      integration, not by a solve (Sec. 6.5 and Appendix A): grad H is
//      constant on a triangle, so is w_T = J grad H, and theta at the edge
//      midpoints of a triangle differs by w_T dotted with their separation.
//      Walking the triangles breadth first from one boundary edge, whose
//      midpoint is given the edge's own tangent, fills in every midpoint. The
//      branch cut of P2 Sec. 5.2.1 is implicit in the walk: the edges it never
//      crosses are the cut.
//
//   3. **The one change ("ours").** P2 puts its singularities on mesh vertices
//      and refines the mesh around each "in the form of bicycle spokes"
//      (Sec. 5.5), because H ~ k ln r and theta ~ k phi have 1/r gradients that
//      no coarse P1 or CR space represents. Here the singular part is taken out
//      analytically instead:
//
//          H     = sum_j k_j ln|z - X_j|  + h,
//          theta = sum_j k_j arg(z - X_j) + theta_h,
//
//      and only the regular parts h and theta_h, which have no singularities
//      and nothing steep to resolve, are discretised. The X_j stay exactly where
//      the first field put them, inside their triangles, and the mesh is never
//      touched. h is harmonic everywhere, so its right-hand side has no
//      interior Diracs at all; at a boundary vertex it is
//
//          b_b = tau_b - 2 pi k_b - Delta_b,
//
//      where Delta_b is the exact change of the singular part along the
//      boundary from the midpoint of the edge before b to the midpoint of the
//      edge after it. That is the statement "theta turns with the boundary
//      tangent between consecutive boundary midpoints" (Appendix A item 5)
//      written for the regular part, so tangency on the anchor's loop is exact
//      by construction, as in P2; and the b_b still sum to
//      2 pi (chi - sum k), so Poincare-Hopf is still the solvability test.
//      Near X_j the field is exactly theta_0 + k_j phi, which is the local model
//      FieldTracer's hyperbola sweep (Viertel et al. Sec. 3.2.2) assumes, and
//      the ports are the closed form of Sec. 7.3, phi_m = (theta_0 + m pi/2) /
//      (1 - k): a spoke mesh is how P2 buys that accuracy, and subtracting the
//      singular part is how it is had for nothing.
//
//   4. Holes (Sec. 6.8). theta is fixed by H up to the anchor's constant, so
//      tangency holds exactly on the anchor's loop and only up to a constant
//      offset delta_m on every other loop: b - 1 real conditions on the
//      singularity *positions*, which a pattern read off an MBO field has no
//      reason to meet (the doc's illustration: a square minus a disk with its
//      four valence-5 singularities rotated 20 degrees off the diagonals comes
//      out 17 degrees off tangent on the hole). The offsets are measured and
//      reported. With Options::alignHoles they are then removed by subtracting
//      delta_m times the harmonic measure of hole m (harmonic, 1 on that hole,
//      0 on every other loop). That leaves theta harmonic with the same
//      singularities and tangent on every loop, which is the canonical harmonic
//      map with these singularities -- the epsilon -> 0 limit of the first
//      field itself -- at the price of H no longer being its exact conjugate:
//      the conformal scale exp(H) picks up a factor going round a hole. That
//      factor is a size, and tracing reads only directions.
//
// ### What it hands on
//
// The regular part theta_h per mesh vertex (the CR values of the incident
// triangles, averaged), and the singularities. FieldTracer evaluates the field
// at a point of triangle t as the P1 interpolant of the regular part plus the
// exact singular part (angleAt() here is the same formula), and launches each
// singularity's ports from the closed form. crossField() is the same field as a
// per-vertex CrossField, u_v = exp(4 i theta(x_v)), for drawing and for
// anything that reads the old representation.
//
// ### Material interfaces
//
// With an interface network (MERIDIAN's Stage 0b Interfaces) the field has to
// be tangent to every interface from both sides, which decouples the
// materials exactly as the first field's interface pins do
// (CrossField::setAlignedInteriorEdges): each region -- a component of the
// model cut along the interfaces -- is a domain of its own, bounded by its
// share of dS and by interface curves, and gets its own Neumann problem, its
// own anchor and its own Poincare-Hopf count. That is done literally: the
// mesh is cut along the interface edges (every vertex duplicated once per
// wedge of triangles the interface edges separate at it), and everything
// above runs unchanged on the cut mesh, whose components are the regions.
// A node of the network (a junction, a landing on dS, a kink) is a corner of
// every region it touches, one per sector, and its quarter counts are not
// Table 1's: they are the node's, the fan's own count (4 inside, Table 1's
// reading of the whole dS corner at a landing) less the index of the first
// field's singular triangles the node absorbed, shared out over the sectors
// by largest remainder -- the rule SeparatrixTrace launches from the nodes
// by, so the field turns at a node exactly as the separatrices leave it
// (nodeQuarters()).
//
// The singular part stays global (every singularity, every region), so the
// field is still one value per vertex of the uncut mesh plus one analytic
// term, and FieldTracer crosses an interface without knowing it is there: the
// regions' regular parts are then each harmonic in their own region, the
// solvability condition of each region's problem is still that region's
// Poincare-Hopf (the global term's change round a region's boundary is 2 pi
// times the index inside it), and at an interface vertex the copies the
// regions solved for agree mod pi/2, both being tangent to the interface.
// They are joined by their mean mod pi/2, and the largest disagreement is
// reported (Report::interfaceMismatch).
//
// What it does not do: the mesh adaptation of Sec. 5.5 (made unnecessary by
// step 3), and the interactive pattern edits of Sec. 5.4 beyond two automatic
// ones -- an opposite pair closer than a mesh edge or two is a dipole the MBO
// left behind and cancels (Options::cancelDipoles), and two singularities in
// one triangle are one singularity of their summed index. Stage 5b
// (PatternRepair) makes the rest.
// ---------------------------------------------------------------------------
class ConformalField {
public:
    // One point of the singularity pattern S: where it is, the triangle holding
    // it, and four times its index (+1 for a valence-3 point, -1 for a 5).
    struct Singular {
        Point pos{0.0, 0.0};
        int triangle = -1;
        int d = 0;
    };

    struct Options {
        // Sec. 6.8: remove each hole's tangency offset delta_m with that hole's
        // harmonic measure (see the header). Off is P2 as written: tangent on
        // the anchor's loop only, delta_m reported.
        bool alignHoles = true;
        // Sec. 5.4's edits that need no user. A +/- pair closer than
        // dipoleRadius mean edges is two halves of nothing -- an MBO
        // singularity pair that had not finished annihilating -- and both go.
        bool cancelDipoles = true;
        double dipoleRadius = 1.5;
        // Sec. 12.3 [P2 Sec. 6.3, Fig. 18]: a +1/4 singularity in the fan of
        // a convex corner makes it a valence-2 corner -- the corner's one quad
        // and the singularity's three collapse into a corner with no quad at
        // all, the configuration P2 repairs by moving the singularity into
        // the patch. Here it is moved before the field is solved: along the
        // corner's bisector, cornerRelocation of the way to the centroid of
        // the triangle the corner makes with its two neighbouring corners
        // (each side taken to the nearer of the two, so the triangle is
        // isosceles about the bisector). Within relocationRadius local mean
        // edges of the corner counts as in its fan.
        //
        // A quarter of the way, measured: on geom028, a star whose six tips
        // hold one each, a quarter gives 48 blocks and a half 78 (the
        // singularity then sends separatrices down both sides of its tip),
        // with the same worst element after smoothing, sin 40 degrees -- the
        // tips' own angle, so the most a quad there can have. Not relocating
        // leaves it at 0.15, and on the first field it was a folded mesh.
        bool relocateCornerSingularities = true;
        double cornerRelocation = 0.25;
        double relocationRadius = 1.5;
    };

    struct Report {
        bool built = false;            // false: see `failure`; nothing else is set
        std::string failure;
        int components = 0;
        int loops = 0;
        int holes = 0;
        int singularities = 0;         // after the edits
        int dipolesCancelled = 0;      // pairs
        int merged = 0;                // singularities merged into another in its triangle
        int relocated = 0;             // moved off a convex corner (Options::relocateCornerSingularities)
        int corners = 0;               // Table 1 corners (q != 2)
        int cornersReread = 0;         // re-read off the first field (see the constructor)
        bool poincareHopf = false;
        // When Poincare-Hopf fails with Table 1's corners: by how many quarter
        // turns the pattern falls short of chi, sum over the components of
        // 4 chi - 4 sum k. Positive is the usual case -- a corner the first
        // field turned as valence 2, so a valence-3 point too few (P2 Sec.
        // 6.3) -- and PatternRepair reads it to know how many to add. Zero
        // when the field built.
        int indexDeficit4 = 0;
        // The same per region (component of the model, cut along the
        // interfaces when there are any; regionOf() says which triangle is
        // in which), when Poincare-Hopf fails; empty when the field built.
        std::vector<int> regionDeficit4;
        // Sec. 6.3: the assembled right-hand side of h sums to zero exactly
        // when the pattern is valid. Its sum, relative to the largest entry or
    // to one radian, whichever is larger.
        double loadSum = 0.0;
        // Sec. 6.5's checks. Crouzeix-Raviart consistency: the largest
        // disagreement (mod pi/2) between the two triangles sharing an interior
        // edge at its midpoint, which is the discrete circulation round every
        // vertex and every hole at once -- zero at a regular vertex, a multiple
        // of pi/2 round a hole.
        double maxMidpointMismatch = 0.0;
        // Tangency, |theta - tangent| mod pi/2 at boundary edge midpoints, worst
        // over the anchor loops: exact by construction.
        double anchorTangency = 0.0;
        // Sec. 6.8: delta_m per hole before the correction, in radians, in
        // (-pi/4, pi/4]; how much it varies along its own hole (zero in exact
        // arithmetic); and the worst tangency on any hole after the correction.
        std::vector<double> holeOffsets;
        double holeSpread = 0.0;
        double holeTangency = 0.0;
        // How far the per-vertex average strays from the CR values it was made
        // of, worst over interior vertices: the price of handing a P1 field on.
        double vertexSpread = 0.0;
        // The regular part's variation, the largest |grad theta_h| times the
        // mean edge: a reading of how steep the part left to the mesh is.
        double maxRegularTurnPerEdge = 0.0;
        // With interfaces: the material regions solved for (components of
        // the cut mesh), the interface nodes whose sectors were prescribed,
        // the first field's singular index those nodes absorbed (quarters),
        // and the largest disagreement mod pi/2 between two regions' regular
        // parts at a vertex of an interface that is no node.
        int regions = 0;
        int interfaceNodes = 0;
        int absorbed4 = 0;
        // Quarter turns added to node sectors the first field turned as no
        // quad at all (valence 2), each a valence-3 point the region is then
        // short of (see the constructor); they show in indexDeficit4.
        int splitSectors = 0;
        double interfaceMismatch = 0.0;
    };

    // `pattern` is the first field's singularities (positions inside their
    // triangles). The boundary corners and their indices are Table 1's, read
    // off the mesh the same way CrossField and SeparatrixTrace read them.
    //
    // `first`, when given, is the field the pattern was read off. It is used
    // only if the pattern fails Poincare-Hopf with Table 1's corners: then the
    // corners are re-read from that field (Sec. 2.4, "measuring k_c from a
    // computed field") wherever the two disagree, which makes the pattern
    // consistent by construction -- an index sum read off one field always
    // adds up. That happens at a corner whose rotation sits on Table 1's
    // boundary, a 225 or 135 degree one, where CrossField::initialize's
    // length-weighted corner data falls on one side of the +-pi/4 wrap or the
    // other with the lengths of the two edges (geom033's four 225-degree
    // corners).
    ConformalField(std::shared_ptr<Mesh> mesh, std::vector<Singular> pattern,
                   const Options &opts, const CrossField *first = nullptr);
    // On a multi-material model: `interfaces` is the network the field is to
    // be tangent to (see the header), and absorbed4 the index, per vertex of
    // the mesh in quarters, of the first field's singular triangles absorbed
    // at an interface node (FieldTracer::absorbedIndexAt) -- it is part of
    // the node's corners, since `pattern` no longer holds it. Null interfaces
    // or a single-material network is the constructor above.
    ConformalField(std::shared_ptr<Mesh> mesh, std::vector<Singular> pattern, const Options &opts,
                   const CrossField *first, const Interfaces *interfaces, std::vector<int> absorbed4);
    ConformalField(std::shared_ptr<Mesh> mesh, std::vector<Singular> pattern)
        : ConformalField(std::move(mesh), std::move(pattern), Options()) {}

    // Per mesh vertex: the corner's interior angle in quarter turns as this
    // field turns there -- Table 1's, or the first field's where it had to be
    // re-read -- 2 at a vertex that is no corner and at every interior one.
    // The tracing launches from the corners by these.
    const std::vector<int> &cornerQuarters() const { return quarters_; }
    // With interfaces, per node of the network (Interfaces::nodes() order),
    // the quarter turns of each of its sectors as this field turns there;
    // empty for a node it did not use (a dangling one) and without interfaces.
    const std::vector<std::vector<int>> &nodeQuarters() const { return nodeQuarters_; }
    const Interfaces *getInterfaces() const { return interfaces_; }
    // Per triangle, the region it is in: a component of the mesh, cut along
    // the interfaces when there are any (set whether or not the field built).
    const std::vector<int> &regionOf() const { return regionOf_; }
    // The node sectors given a quad the first field did not turn there (see
    // Report::splitSectors): the node's vertex and a triangle of the sector,
    // whose region (regionOf()) is the one left a valence-3 point short.
    struct SplitCorner { int vertex = -1; int triangle = -1; };
    const std::vector<SplitCorner> &splitCorners() const { return splitCorners_; }

    bool valid() const { return report_.built; }
    const Report &getReport() const { return report_; }
    const Options &getOptions() const { return opts_; }
    const Mesh &getMesh() const { return *mesh_; }
    std::shared_ptr<Mesh> getMeshPtr() const { return mesh_; }

    // The pattern after the edits: what the field has, and all it has.
    const std::vector<Singular> &singularities() const { return sing_; }

    // theta_h at each vertex, in radians, on an arbitrary branch per vertex
    // (it is a cross angle: only its value mod pi/2 means anything, and the
    // tracer lifts the three of a triangle to agree before interpolating).
    const std::vector<double> &regularAngles() const { return regular_; }

    // The singular part sum_j k_j arg(x - X_j), mod pi/2 (each term is
    // multivalued by 2 pi k_j, a multiple of pi/2, so the sum of principal
    // arguments is a fine representative).
    double singularPart(const Point &x) const;

    // The cross angle at x in triangle t, mod pi/2: the P1 interpolant of the
    // regular part plus the singular part. On the branch nearest `ref`.
    double angleAt(int t, const Point &x, double ref = 0.0) const;

    // h, the regular part of log |F'|, per vertex (pinned to 0 at one vertex
    // per component). exp(H) is P2's size map; nothing here uses it.
    const std::vector<double> &regularLogScale() const { return h_; }

    // The same field as a per-vertex CrossField, u_v = exp(4 i theta(x_v)).
    std::shared_ptr<CrossField> crossField() const { return cross_; }

private:
    void build();

    std::shared_ptr<Mesh> mesh_;
    Options opts_;
    Report report_;
    const CrossField *first_ = nullptr;
    // The interface network and the absorbed indices (read during
    // construction), and the per-node sector quarters it settled on.
    const Interfaces *interfaces_ = nullptr;
    std::vector<int> absorbed4_;
    std::vector<std::vector<int>> nodeQuarters_;
    std::vector<int> regionOf_;
    std::vector<SplitCorner> splitCorners_;
    std::vector<int> quarters_;
    std::vector<Singular> sing_;
    std::vector<double> regular_;
    std::vector<double> h_;
    std::shared_ptr<CrossField> cross_;
};

#endif // __CONFORMAL_FIELD_HXX__
