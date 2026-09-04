#ifndef __PAPER_DOMAINS_HXX__
#define __PAPER_DOMAINS_HXX__

#include <memory>
#include <string>
#include <vector>

#include "mesh/Mesh.hxx"

// The domains Sec. 6 of the paper runs on that are not files.
//
// Two kinds live here. The **canonical** ones -- square, disk, annulus,
// L-shape, wedges, polygons -- are the domains of E1 and E2, and they are built
// in memory rather than read from data/meshes for one reason: E1 refines them,
// and a refinement sequence has to be the *same domain* at five mesh sizes, not
// five files that happen to look alike. The exact geometry is also what makes
// the expected answer knowable in advance -- the disk's four +1/4 cones sit on
// the diagonals, the annulus carries none -- and E2 is a comparison against
// that expectation.
//
// The **junction** ones are E5's, and the outline names them: a T-junction, a
// quadruple point, a thin layer, an embedded inclusion, and a junction whose
// sectors are not right angles. data/meshes/multimat has nine models and E5
// runs on all of them, but they are models of things rather than of junctions,
// and the close-ups Fig. 6 wants are easier to read on a domain built to show
// one junction and nothing else.
//
// Every multi-material domain here is built so that **the interfaces are edges
// of the triangulation**, because Sec. 7 says plainly that the pipeline assumes
// that and does not generate it: the interface polylines go into the mesher as
// constrained segments, and the material id of each triangle is Triangle's own
// region attribute, flooded from a seed point per material.
namespace paper {

// A domain, with whatever is known about it before the field is computed.
struct Domain {
    std::string name;
    std::shared_ptr<Mesh> mesh;

    // Four times the sum of the interior singularity indices a correct
    // boundary-aligned cross field must carry on this domain, where that is
    // known in closed form, and `false` in `known` where it is not.
    //
    // It is derivable from the geometry alone: Poincare-Hopf reads
    //
    //     sum_{v interior} 4 I(v)  =  4 chi(Omega) - sum_{v in dS} q(v),
    //
    // with q(v) the number of quarter turns the boundary makes at v (Metrics::
    // cornerQuarters). So this field is not an independent expectation so much
    // as the *statement of the identity*, written down before the code runs so
    // that E2 can report a disagreement rather than a tautology.
    bool knownIndex = false;
    int expectedInteriorIndex4 = 0;

    std::string note;
};

// --- the canonical domains of E1 and E2 -----------------------------------

// [0,1]^2. Four 90-degree corners take the whole of chi, so a correct field has
// no interior singularity at all and is in fact constant.
std::shared_ptr<Mesh> squareDomain(double h);

// The unit disk. No corner carries anything, so the whole of 4 chi = 4 sits in
// the interior: four +1/4 cones, and with DualMBO's disk-centre pin they sit at
// 45 + k*90 degrees (DualMBO::setPinDiskCenters).
std::shared_ptr<Mesh> diskDomain(double h);

// An annulus, r in [rIn, rOut]. chi = 0 and neither rim has a corner, so a
// correct field carries no interior singularity.
std::shared_ptr<Mesh> annulusDomain(double rIn, double rOut, double h);

// The L, [0,1]^2 minus the upper-right quarter: five convex corners and one
// reflex, 5 - 1 = 4 quarter turns, so again no interior singularity is needed.
std::shared_ptr<Mesh> lShapeDomain(double h);

// A circular sector of the given sweep. The apex is the interesting vertex:
// its interior angle is neither a right angle nor straight, so the quarter
// count there is a rounding and the interior has to make up the remainder.
std::shared_ptr<Mesh> wedgeDomain(double sweepDegrees, double radius, double h);

// A regular n-gon, which for n not 4 has no right angle anywhere on it.
std::shared_ptr<Mesh> regularPolygonDomain(int n, double h);

// The unit disk again, but meshed by `rings` concentric rings whose radial
// spacing grows geometrically outward, dr_{i+1} = growth * dr_i, with the first
// spacing chosen so the outermost ring lands exactly on the boundary and the
// angular count of each ring set to keep its triangles about isotropic. The
// triangulation is the constrained Delaunay one of exactly those points -- no
// Steiner points, no quality refinement -- so `growth` is a direct knob on how
// graded the mesh is at a fixed number of rings, and `growth = 1` is the
// uniform concentric mesh.
//
// It exists because the corpus is near-uniform (area ratios across an edge
// reach only 2.2) and the inconsistency of the interior-penalty weight is an
// O(1) *spurious force whose size scales with the mismatch between the two
// elements on an edge*. On a near-uniform mesh a provable defect can still be
// invisible in the field, and the honest way to find out is to grade the mesh
// until it is not near-uniform and measure again -- both the operator's
// residual on an affine field and, on the one domain whose answer is known in
// closed form, where the cones actually land. That is E1(f).
//
// The disk rather than a square because the disk is the domain with a known
// answer that a graded mesh does not trivialise: four +1/4 cones at
// 45 + k*90 degrees from the centre (fixed by the centre pin,
// DualMBO::setPinDiskCenters). On the square and the L every boundary tangent
// is axis-aligned, so exp(4 i theta) = 1 on every boundary edge, the constant
// field is the exact minimiser for *any* weight, and no discretisation can be
// wrong there.
std::shared_ptr<Mesh> gradedDiskDomain(double radius, int rings, double growth);

// The canonical set as E2 runs it, at one mesh size.
std::vector<Domain> canonicalDomains(double h);

// Move every *interior* vertex by up to `frac` of the local edge length, in a
// direction drawn from `seed`. The boundary is left exactly where it was, so
// the domain -- and every expectation above -- is unchanged, and only the
// triangulation of it differs. This is E2's remeshing robustness check.
std::shared_ptr<Mesh> jitter(const Mesh &src, double frac, unsigned seed);

// --- E5's junction domains -------------------------------------------------

// Three materials meeting along a T: the stem runs down from the centre of a
// square and the crossbar runs the width of it.
std::shared_ptr<Mesh> tJunctionDomain(double h);

// Four materials meeting at one point, the two interfaces crossing at right
// angles. The quadruple point is the vertex at which a vertex-based field has
// four mutually exclusive tangency constraints.
std::shared_ptr<Mesh> quadruplePointDomain(double h);

// A thin band of a second material across a block: two nearly parallel
// interfaces a small multiple of the element size apart.
std::shared_ptr<Mesh> thinLayerDomain(double h);

// One circular inclusion in a matrix block.
std::shared_ptr<Mesh> embeddedInclusionDomain(double h);

// Three materials meeting at 120 degrees, so that no sector of the junction is
// a multiple of a right angle and no single cross can be tangent to all three
// interfaces.
//
// `sectorDegrees` is the angle of the first two sectors; the third takes the
// rest. 120 is the paper's domain. Other values exist so that the layout
// stage's behaviour at a junction can be measured as a function of how far the
// sectors are from a quarter turn: 90 is a T-junction in disguise (90/90/180),
// 100 is 100/100/160, and so on.
std::shared_ptr<Mesh> obliqueJunctionDomain(double h, double sectorDegrees = 120.0);

// The five above, at one mesh size. `obliqueSectorDegrees` is passed to the
// oblique domain; 120 is the paper's.
std::vector<Domain> junctionDomains(double h, double obliqueSectorDegrees = 120.0);

// --- the mechanism set -----------------------------------------------------
//
// Domains built to test *where the two discretisations must differ*, rather
// than to sample an application. The distinction matters and the set is kept
// apart from data/meshes/{singlemat,multimat} for it: a corpus you choose after
// seeing the results is not evidence. What makes these legitimate is that each
// family sweeps one parameter, the prediction is stated before the run, and
// each family contains the setting at which the prediction says the two methods
// *agree* -- so the result is a curve with a control at one end and not a win.
//
// The mechanism is the boundary data. A vertex-based field must assign one
// cross to a boundary vertex, and `CrossField` does it by rounding the corner's
// interior angle into one of four quarter-turn classes (Viertel, Osting and
// Staten, Table 1). That rounding is exact at 90, 180 and 270 degrees and wrong
// by up to 22.5 degrees half way between, and the same argument at a material
// junction is Remark 4.1. A face-based field pins each boundary *edge* to its
// own tangent and never rounds. So:
//
//   comb        all interior angles 90 or 270      -> control, no difference
//   star        tip angle swept 90 down to 40      -> error grows, both methods
//                                                     eventually (ours at the
//                                                     corner-face degeneracy)
//   laminate    interfaces meeting dS at an angle  -> junction residual grows
//   grain       Voronoi polycrystal, every triple
//               junction generically oblique       -> the application case

// A rectilinear comb: a block with `teeth` slots cut into it, every interior
// angle 90 or 270 degrees. The control: the corner rounding a vertex field does
// is exact here, so the two methods should agree, and if they do not the
// explanation is not the one this set is testing.
std::shared_ptr<Mesh> combDomain(int teeth, double h);

// A `points`-pointed star whose tip interior angle is `tipDegrees`. The valley
// angle follows, 360 - 360/points - tipDegrees, so one parameter sweeps both
// away from a multiple of 90 degrees. At 90 degrees the tips are right angles
// and only the valleys are odd.
std::shared_ptr<Mesh> starDomain(int points, double tipDegrees, double h);

// `layers` material bands across the unit square, their interfaces at
// `tiltDegrees` to the x axis. Each interface ends on dS, and where it does the
// boundary tangent and the interface tangent are two constraints on one corner:
// a junction of obliquity |tilt - 90| rounded into a right angle. 90 degrees is
// the control.
std::shared_ptr<Mesh> laminateDomain(int layers, double tiltDegrees, double h);

// A Voronoi polycrystal of `grains` cells: the aggregate of the Voronoi cells
// of a jittered set of seeds, one material per cell, meshed so that every cell
// wall is a chain of mesh edges. Almost every interior vertex is a triple
// junction whose three sectors are generically not multiples of a right angle,
// which is the configuration Remark 4.1 is about, and it is what a
// polycrystalline microstructure actually looks like.
std::shared_ptr<Mesh> polycrystalDomain(int grains, unsigned seed, double h);

// The whole set, with the note each Domain carries saying what it tests.
std::vector<Domain> mechanismDomains(double h);

} // namespace paper

#endif // __PAPER_DOMAINS_HXX__
