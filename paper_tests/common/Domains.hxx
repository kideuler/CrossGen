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
std::shared_ptr<Mesh> obliqueJunctionDomain(double h);

// The five above, at one mesh size.
std::vector<Domain> junctionDomains(double h);

} // namespace paper

#endif // __PAPER_DOMAINS_HXX__
