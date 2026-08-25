#ifndef __MERIDIAN_HXX__
#define __MERIDIAN_HXX__

#include <memory>
#include <string>
#include <vector>

#include "MERIDIAN/ConeCut.hxx"
#include "MERIDIAN/ConeSingularities.hxx"
#include "MERIDIAN/RicciFlow.hxx"
#include "mesh/Mesh.hxx"

class SIPG;

// Stages 1 to 3 of
//   Shepherd, Gu and Hughes, "Isogeometric model reconstruction of open shells
//   via Ricci flow and quadrilateral layout-inducing energies", Engineering
//   Structures 252 (2022) 113602.   (docs/shepherd2022.pdf)
//
// The paper's contribution is to stop treating a quadrilateral layout as a
// combinatorial local-to-global problem and to characterise it instead as one
// continuous map Psi : S - G -> R^2 satisfying five properties Q1-Q5
// (Definition 2.1). That turns the whole thing into a variational problem, and
// it needs a starting map that is already valid in the properties the
// variational stage cannot repair -- Q1 in particular, since the symmetric
// Dirichlet barrier of Sec. 3.3 can prevent a flipped triangle but never undo
// one. Ricci flow is where that starting map comes from, and getting to it
// takes three stages:
//
//   Stage 1  Sec. 3.1    Cone singularities. Cross field in, integer indices
//                        out, checked against the discrete Gauss-Bonnet
//                        condition sum I(v) = 4 chi(S) of Eq. (4).
//                        -> ConeSingularities
//
//   Stage 2  Sec. 3.2.2  Cutting graph G, and Omega = S - G. Every void opened,
//                        every interior cone dragged to the boundary, the
//                        result a topological disk.
//                        -> ConeCut
//
//   Stage 3  Sec. 3.2.1  Discrete surface Ricci flow. Replace the metric
//                        inherited from the plane by a flat metric whose
//                        curvature is zero everywhere except (pi/2) I(v) at the
//                        cones.
//                        -> RicciFlow
//
// The Gauss-Bonnet check between Stages 1 and 3 is not a formality. Newton's
// system in the flow is Delta du = Kbar - K with Delta a Laplacian, whose
// kernel is the constants; the residual has to be orthogonal to that kernel for
// the system to be consistent at all, and sum(Kbar) = sum(K) = 2 pi chi is
// exactly that condition. An inadmissible cone set does not converge slowly,
// it has no solution. run() therefore stops there rather than handing an
// inconsistent system to the solver.
//
// What comes out is the flat cone metric, as one conformal factor u_i per
// vertex (Eq. 8) plus the fixed inversive distances, together with the cut disk
// it is to be laid out on. Stage 4 -- isometric propagation of that metric into
// the plane, Sec. 3.2.2's second half -- is the next thing to build, and is not
// implemented here.
class MERIDIAN {
public:
    struct Options {
        // Stage 0: the cross field the cones are read off. The defaults match
        // TestSIPG.
        double sipgGamma = 10.0;
        int sipgMaxSteps = 500;

        // Stage 1
        bool autoRebalance = true;   // restore Eq. (4) by moving boundary cones
        int minBoundaryIndex = -3;
        int maxBoundaryIndex = 1;

        // Stage 3
        double ricciTolerance = 1e-8;
        int ricciMaxIterations = 100;
        bool delaunayFlips = true;
    };

    struct Status {
        int mboSteps = 0;
        bool fieldConverged = false;

        bool conesAdmissible = false;
        int interiorCones = 0;
        int boundaryCones = 0;

        bool cutIsDisk = false;
        bool allConesOnBoundary = false;

        bool ricciConverged = false;

        // Every stage did what Definition 2.1 needs of it, so Stage 4 has a
        // flat cone metric on a disk to immerse.
        bool readyForImmersion = false;

        std::vector<std::string> messages;
    };

    explicit MERIDIAN(std::shared_ptr<Mesh> mesh);
    MERIDIAN(std::shared_ptr<Mesh> mesh, const Options &opts);
    ~MERIDIAN();

    // Runs the three stages in order and stops at the first one that fails in a
    // way the next stage cannot survive: an inadmissible cone set stops the
    // pipeline, a cut that came out as something other than a disk does not
    // (Ricci flow does not use it, so the metric is still worth computing and
    // the failure is still worth reporting).
    bool run();

    const SIPG& getField() const { return *field; }
    const ConeSingularities& getCones() const { return *cones; }
    const ConeCut& getCut() const { return *cutter; }
    const RicciFlow& getRicci() const { return *ricci; }

    // Null until the stage that builds them has run.
    bool hasCones() const { return cones != nullptr; }
    bool hasCut() const { return cutter != nullptr; }
    bool hasRicci() const { return ricci != nullptr; }

    const Status& getStatus() const { return status; }

    const Mesh& getMesh() const { return *mesh; }
    std::shared_ptr<Mesh> getMeshPtr() const { return mesh; }

private:
    void runField();

    std::shared_ptr<Mesh> mesh;
    Options options;
    Status status;

    std::unique_ptr<SIPG> field;
    std::unique_ptr<ConeSingularities> cones;
    std::unique_ptr<ConeCut> cutter;
    std::unique_ptr<RicciFlow> ricci;
};

#endif // __MERIDIAN_HXX__
