#ifndef __CONESINGULARITIES_HXX__
#define __CONESINGULARITIES_HXX__

#include <complex>
#include <memory>
#include <string>
#include <vector>

// eigen includes
#include <Eigen/Dense>

#include "mesh/Mesh.hxx"

class DualMBO;

// Stage 1 of
//   Shepherd, Gu and Hughes, "Isogeometric model reconstruction of open shells
//   via Ricci flow and quadrilateral layout-inducing energies", Engineering
//   Structures 252 (2022) 113602, Sections 3.1 and 3.2.1.
//   (docs/shepherd2022.pdf)
//
// Ricci flow needs a *prescribed* target curvature, not a measured one: a set
// P of cone singularities with integer indices, and zero curvature everywhere
// else. Eq. (3) defines the index of a vertex from the immersion that does not
// exist yet,
//
//     I(v) = (2/pi) ( 2pi - sum_T angle(v, Psi(T)) )   v not on the boundary
//     I(v) = (2/pi) (  pi - sum_T angle(v, Psi(T)) )   v on the boundary
//
// so in practice it is read off a cross field instead, which is what Sec. 3.1
// does with the frame field of [14]. Here the field is the p=0 dual-mesh MBO
// field of DualMBO, and the two halves of Eq. (3) are recovered as follows.
//
// Interior vertices. The mesh is planar, so sum_T angle(v, T) is exactly 2pi
// and the whole of the index sits in the field's holonomy: the cross turns by
// (pi/2) I(v) over the vertex star. DualMBO::computeSingularities already
// measures that as a winding number of exp(4 i theta), and I(v) is that
// winding number -- integer by construction, with no rounding anywhere.
//
// Boundary vertices. Here the interior angle Omega is whatever the mesh
// happens to have and the field alone cannot say what the index is, so the
// two are combined the way UMBER::boundarySingularities() does it:
//
//     I(v) ~ ( Theta + pi - Omega ) / (pi/2),
//
// Theta being the rotation of the cross across the star walked from one
// boundary edge to the other. A straight run of boundary that the field
// follows has Theta = 0, Omega = pi and reads 0; a square corner has
// Omega = pi/2 and reads +1 (valence two); a reentrant corner has
// Omega = 3pi/2 and reads -1. This one *is* rounded, and it is rounded along
// each boundary loop rather than vertex by vertex, because a discrete field
// spreads a corner over the two or three vertices it takes to swing across --
// see the implementation.
//
// Summed over every vertex the unrounded quantities give 4 chi(S) exactly:
// each interior edge appears in the star walk of both its endpoints with
// opposite sign and cancels, leaving only the total boundary turning 2 pi chi.
// That is the discrete Gauss-Bonnet condition of Eq. (4), and it is what
// admissible() checks. What breaks it is the rounding on the boundary, so
// rebalance() puts it back by moving whole index units between boundary cones
// -- Sec. 3.1's own remedy ("the problematic cones should be repositioned"),
// automated by taking the units from wherever they cost the least.
//
// The target curvature handed to Ricci flow is then Eq. (9)'s Kbar:
//
//     Kbar_i = (pi/2) I(v_i)   for v_i in P,     0 otherwise,
//
// and sum_i Kbar_i = (pi/2) 4 chi = 2 pi chi is the Eq. (7) hypothesis the
// flow needs to have a solution at all.
class ConeSingularities {
public:
    // One prescribed cone. `raw` is the unrounded quarter-turn count the
    // rounding came from -- 1.00 at a clean square corner, 0.35 at a vertex
    // carrying a third of one -- and is worth reading when a cone lands
    // somewhere the geometry does not justify.
    struct Cone {
        int vertex = -1;
        int index = 0;          // I(v), the integer of Eq. (3)
        bool onBoundary = false;
        int valence = 4;        // 4 - I inside, 3 - I on the boundary
        double raw = 0.0;       // I(v) before the integer snap
        bool fromRebalance = false; // index moved by rebalance(), not measured
        bool prescribed = false;    // set by prescribe(), not read off the field
    };

    // Everything Eq. (4) and Eq. (7) have to say about the current cone set.
    // Both are checked: the index form is the admissibility condition the flow
    // needs, the curvature form is a check on the mesh itself.
    struct GaussBonnetReport {
        int V = 0, E = 0, F = 0;        // V counts only vertices carrying a triangle
        int isolatedVertices = 0;       // vertices no triangle references
        int eulerCharacteristic = 0;    // chi = V - E + F
        int boundaryLoops = 0;

        int indexSum = 0;               // sum_v I(v)
        int indexTarget = 0;            // 4 chi, Eq. (4)
        double rawIndexSum = 0.0;       // the same sum before rounding

        double curvatureSum = 0.0;      // sum_v K_v on the input metric, Eq. (6)
        double targetCurvatureSum = 0.0;// sum_v Kbar_v
        double curvatureTarget = 0.0;   // 2 pi chi, Eq. (7)
        double curvatureResidual = 0.0; // |curvatureSum - curvatureTarget|

        int rebalanceUnits = 0;         // index units moved to restore Eq. (4)
        double rebalanceCost = 0.0;     // sum |I - raw| added by that move

        bool admissible = false;        // indexSum == 4 chi
        bool metricConsistent = false;  // curvatureSum == 2 pi chi to tolerance
        std::vector<std::string> messages;
    };

    // From a converged DualMBO field. Run dualMBO.initialize() and step it (or
    // runMBO()) before constructing this: the field is read from
    // DualMBO::u_k_prev, which is what computeSingularities() also uses, and an
    // unconverged field simply gives cones in the wrong places.
    //
    // DualMBO::computeSingularities() is refreshed from the same field on the way
    // through, so anything downstream reading DualMBO::singularVertices sees the
    // same interior windings that became cones here.
    explicit ConeSingularities(DualMBO &dualMBO);

    // Direct form, for a field that did not come from DualMBO: one exp(4 i theta)
    // per triangle, in the same convention as DualMBO::u_k.
    ConeSingularities(std::shared_ptr<Mesh> mesh, const Eigen::VectorXcd &crossField);

    // Overwrite I(v) at vertices whose index the geometry fixes rather than the
    // field -- the nodes of the material interface network of Stage 0b.
    //
    // The field cannot supply these and it is not a matter of resolution. A
    // cross field is a section of a bundle over the *surface*; it knows nothing
    // about the material tags, so at a triple junction it reports whatever the
    // boundary conditions elsewhere happened to propagate, and at a landing it
    // reports the boundary corner it can see and not the interface it cannot.
    // Meanwhile the layout has no freedom there at all: every sector between
    // two incident interfaces has to be a whole number of quadrilaterals, which
    // fixes the cone angle, which fixes I(v). See Interfaces::quantiseNode.
    //
    // A prescribed index is exact in the same sense an interior winding number
    // is, so rebalance() will not move it and reports it as unavailable.
    // Prescribing an index of zero is meaningful and is the common case: it
    // says the geometry has looked and there is no cone here, which a field
    // that read one has to be overruled on.
    void prescribe(const std::vector<std::pair<int, int>> &fixed);

    // Whether prescribe() set the index at each vertex.
    const std::vector<char>& getPrescribed() const { return prescribed; }
    // How much the prescription moved the field's own reading, summed over the
    // vertices it touched: the honest measure of how far the field was from
    // knowing the interfaces were there.
    int prescriptionShift() const { return lastPrescriptionShift; }

    // Index units cancelDipoles() annihilated.
    int cancelledUnits() const { return lastCancelledUnits; }

    // Annihilate interior cones of opposite index that lie in the same group.
    //
    // `group` is one label per vertex; only two cones sharing a non-negative
    // label may cancel. Stage 0b supplies the material region, which is the
    // label that matters: a cross field aligned to a curved interface puts a
    // +1/-1 pair on the concave side of every strongly curved stretch of it --
    // three of them on geom011's corrugation, one at each extremum of the
    // cosine -- and both members of such a pair are in the *same* region. The
    // pair therefore contributes nothing to that region's Gauss-Bonnet count
    // and nothing to Eq. (4), and it is a property of the smoothest field
    // rather than of the layout: removing it takes away two cones, their two
    // cuts and the patches their separatrices carve, and leaves every count
    // exactly where it was.
    //
    // The pair that *spans* an interface is the opposite case and is never
    // touched, because the two labels differ: geom001's +1 sits inside the
    // quarter disk and its -1 outside, and each of them is the whole of its own
    // region's deficit. Cancelling those two is precisely the mistake that
    // leaves the quarter disk one quarter turn short and sends Stage 0b looking
    // for somewhere to put a corner.
    //
    // Closest pair first, so what goes is the tightest cluster -- Sec. 3.1's
    // own remedy for cones the field put somewhere the geometry does not
    // justify, in the case where the two indices sum to zero and repositioning
    // them onto each other annihilates them. Prescribed cones are never moved,
    // and neither are cones on dS: their index is the rounded one and
    // rebalance() is what is written against them.
    //
    // Returns the number of index units cancelled.
    int cancelDipoles(const std::vector<int> &group);

    // Move index units between boundary cones until Eq. (4) holds. Each unit
    // goes to whichever boundary vertex is cheapest, cost being how much
    // further it drags that vertex's index from the unrounded value the field
    // and the geometry actually gave -- so a corner the field read as 0.51 and
    // rounded up to 1 is the first thing to give a unit back.
    //
    // Boundary indices are held in [minBoundaryIndex, maxBoundaryIndex]:
    // I = 1 is already a valence-two corner and there is nothing above it,
    // while the lower end is only there to stop a runaway on a badly broken
    // field. Interior cones are never touched -- they are integers measured
    // with no rounding, so an imbalance is never their fault.
    //
    // Returns the number of units moved (0 if it was already admissible), or
    // -1 if the residual could not be placed at all.
    int rebalance();

    // I(v) for every vertex, zero away from the cones.
    const std::vector<int>& getIndices() const { return index; }

    // Kbar_i = (pi/2) I(v_i): the target curvature of Eq. (9), one per vertex.
    const std::vector<double>& getTargetCurvature() const { return targetCurvature; }

    // K_i of Eq. (6) on the input metric. Zero in the interior for a planar
    // mesh; on the boundary it is pi minus the interior angle, i.e. the
    // turning, which is the curvature the flow has to move into the cones.
    const std::vector<double>& getInputCurvature() const { return inputCurvature; }

    // The cones themselves, interior first then boundary, each group ordered
    // by vertex index.
    const std::vector<Cone>& getCones() const { return cones; }

    // Whether each vertex carries any triangle at all. An .obj can list
    // vertices no face references -- several models in data/meshes do -- and
    // such a vertex is not part of the surface: it adds 1 to V and nothing to
    // E or F, so counting it inflates chi and Eq. (4) then asks for cones that
    // should not exist. Everything here is summed over the active vertices
    // only, and RicciFlow keeps them out of its Newton system for the same
    // reason (an isolated vertex contributes an all-zero Laplacian row).
    const std::vector<char>& getActiveVertices() const { return active; }

    // Just the interior ones, as DualMBO reports them.
    std::vector<Cone> interiorCones() const;
    std::vector<Cone> boundaryCones() const;

    // Recomputed on demand, so it reflects any rebalance() that has happened.
    GaussBonnetReport gaussBonnet() const;

    bool admissible() const { return gaussBonnet().admissible; }

    void setBoundaryIndexRange(int lo, int hi) { minBoundaryIndex = lo; maxBoundaryIndex = hi; }

    const Mesh& getMesh() const { return *mesh; }
    std::shared_ptr<Mesh> getMeshPtr() const { return mesh; }

private:
    void computeInputCurvature();
    void computeInteriorIndices(const Eigen::VectorXcd &field);
    void computeBoundaryIndices(const Eigen::VectorXcd &field);
    void rebuildCones();

    // The star of a boundary vertex, walked from one boundary edge to the
    // other and oriented counter-clockwise. Empty if v is a pinch (more than
    // two boundary edges) or the walk does not close on a boundary edge.
    std::vector<int> boundaryFan(int v) const;

    std::shared_ptr<Mesh> mesh;

    std::vector<char> active;           // vertex carries at least one triangle
    std::vector<int> index;             // I(v) per vertex
    std::vector<double> rawIndex;       // the unrounded value per vertex
    std::vector<char> measured;         // whether rawIndex[v] means anything
    std::vector<char> movedByRebalance; // index units placed by rebalance()
    std::vector<char> prescribed;       // index fixed by prescribe(), not measured
    std::vector<double> inputCurvature; // K_v, Eq. (6)
    std::vector<double> targetCurvature;// Kbar_v, Eq. (9)

    std::vector<Cone> cones;

    int minBoundaryIndex = -3;
    int maxBoundaryIndex = 1;
    int lastPrescriptionShift = 0;
    int lastCancelledUnits = 0;
    int lastRebalanceUnits = 0;
    double lastRebalanceCost = 0.0;
};

#endif // __CONESINGULARITIES_HXX__
