#ifndef __PAPER_METHODS_HXX__
#define __PAPER_METHODS_HXX__

#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "Metrics.hxx"
#include "mesh/Mesh.hxx"

// The three methods E2 to E5 compare, behind one interface.
//
// All three produce the same object -- one unit spin-4 value per triangle --
// and that is the fair-comparison protocol of the outline made concrete: the
// energies, the singularities and the alignment errors are then read off the
// same representation by the same code, and the only thing that differs between
// the rows of a table is what put the numbers in the vector.
//
//   SIPG   ours. Face DOFs natively; nothing is converted.
//   B1     P1-MBO (`CrossField`): vertex crosses and face singularities,
//          followed by the documented, best-effort conversion to faces that
//          Sec. 6 requires us to describe and to be fair about. What the
//          conversion costs is measured rather than asserted -- see
//          metrics::ConversionReport and `convertP1ToFaces`.
//   B2     the polyvector field (`PolyField`): face DOFs, a single linear
//          least-squares solve, no MBO. It is in the comparison to separate
//          "face-based" from "MBO" as the source of any improvement.
namespace paper {

struct MethodOptions {
    // SIPG and the common energy alike. The outline picks 10 after E1's sweep.
    double gamma = 10.0;
    int maxSteps = 500;
    // Multiplier on SIPG's tau = D^2/10. E1(b) is the only caller that moves it.
    double tauScale = 1.0;
    bool pinDiskCenters = true;
    // E5(a)'s ablation: hard-pinned Dirichlet rows, or the weak Nitsche form.
    bool hardBoundary = true;
    // Align to material interfaces where the mesh has them. Off is the
    // single-material behaviour and is what a mesh with one material gets
    // regardless.
    bool alignInterfaces = true;
    // B1's random start. B1 is seeded because it has to be: `CrossField`
    // initialises its interior at random and the outline asks for mean and
    // sigma over ten seeds.
    unsigned seed = 1;
    // Keep the per-iteration energies and increment. E1(d) is the only caller
    // that wants them.
    bool recordHistory = false;
    // Take every one of `maxSteps` iterations rather than stopping at the
    // convergence test. E1(d) needs this: at tau = D^2/10 one step diffuses
    // most of the way to the fixed point and the loop exits after two or three,
    // which is a convergence history with nothing in it to look at.
    bool forceSteps = false;
    // The convergence test both MBO methods use, as `error < 2 N tol` with N the
    // number of degrees of freedom -- triangles for SIPG, vertices for B1. 1e-5
    // is what the pipeline ships and is the right number for a layout, which
    // only needs the singularities. It is *not* tight enough for an energy: on
    // the unit disk it stops after three steps at an energy 5% above the fixed
    // point's (E1(d) measures exactly this). The field-quality experiments run
    // at 1e-9 and say so; the pipeline experiments leave it at the shipped value
    // so that what they measure is the pipeline.
    double convergenceTol = 1e-5;

    // --- tau-continuation ---------------------------------------------------
    //
    // The MBO step diffuses for a time tau, and on this operator the length it
    // smooths over is
    //
    //     ell = sqrt(8 gamma tau)
    //
    // (K's diagonal is ~3 kappa ~ 3 gamma and M's is the area ~ h^2/2, so the
    // graph Laplacian M^-1 K has eigenvalues up to ~8 gamma / h^2, and a mode of
    // wavelength ell is damped once 8 gamma tau / ell^2 >~ 1).
    //
    // At the shipped tau = D^2/10 that is ell ~ 2.8 D: one step smooths across
    // several domain diameters, tau*K swamps M, and the iteration reaches its
    // fixed point after a *single* solve -- measured, the energy at step 1 and
    // at step 300 agree to every digit. The scheme is then a harmonic extension
    // followed by a normalisation, which is why it returns the same singularity
    // count as the one-shot polyvector baseline on all 35 models, and the
    // threshold dynamics never actually runs.
    //
    // The fix is the MBO analogue of Ginzburg-Landau's epsilon-continuation:
    // anneal tau down a geometric ladder, re-solving from the previous level's
    // field, until the diffusion length reaches a few mesh edges -- below that
    // the step resolves nothing the mesh can carry. Measured on the 35-model
    // corpus this lowers the common energy by 22% in total and turns a loss to
    // B1 on energy into a win on 27 of 35 models, with boundary alignment still
    // exact and the Poincare-Hopf residual still zero.
    bool tauContinuation = true;
    // Floor, as a multiple of the mean edge length: stop once ell <= this * h.
    double tauFloorEdges = 20.0;
    // Ratio between consecutive levels.
    double tauRatio = 0.25;
    // Per-level iteration cap.
    int tauLevelSteps = 2000;

    // Give B1 the same ladder. Off by default, because the baseline belongs in
    // the table as it is published; on, it answers the reviewer's question --
    // is the gain the continuation, or the face-based discretisation? -- by
    // measuring both methods with the continuation and both without.
    bool b1TauContinuation = false;
};

struct FieldRun {
    std::string method;
    Eigen::VectorXcd u;          // one unit spin-4 value per triangle

    bool ok = false;
    std::string error;

    int iterations = 0;
    bool converged = false;
    double assembleSeconds = 0.0;   // build + factorise
    double solveSeconds = 0.0;      // the iteration loop, or the single solve
    double totalSeconds = 0.0;

    // Filled when MethodOptions::recordHistory, one entry per step: the
    // comparison energy of Sec. 5, the method's own energy with the Dirichlet
    // terms in it, and the MBO increment ||u^{k+1} - u^k||.
    std::vector<double> energyHistory;
    std::vector<double> sipgEnergyHistory;
    std::vector<double> incrementHistory;

    // B1 only.
    bool hasConversion = false;
    metrics::ConversionReport conversion;
    // B1's field before the conversion, one value per *vertex*, kept so that
    // E3 can measure what the resampling did.
    Eigen::VectorXcd p1Vertex;
};

// Interior edges whose two triangles carry different material ids. This is what
// SIPG::setAlignedInteriorEdges is handed and what the interface metrics are
// measured over; it is computed here directly from the material ids rather than
// through `Interfaces` so that the harness depends on the field code and on
// nothing downstream of it.
std::vector<int> interfaceEdges(const Mesh &m);
std::vector<char> interfaceEdgeFlags(const Mesh &m);
bool isMultiMaterial(const Mesh &m);

FieldRun runSIPG(const std::shared_ptr<Mesh> &m, const MethodOptions &o);
FieldRun runP1MBO(const std::shared_ptr<Mesh> &m, const MethodOptions &o);
FieldRun runPolyVector(const std::shared_ptr<Mesh> &m, const MethodOptions &o);

// The three, in the order the tables print them.
std::vector<std::string> methodNames();
FieldRun runMethod(const std::string &name, const std::shared_ptr<Mesh> &m, const MethodOptions &o);

// --- the B1 conversion, in the open ----------------------------------------
//
// Sec. 6 says the strawman has to be a fair one and has to be described in the
// paper, so it is one function and this is it.
//
//   face value  the spin-4 average of the triangle's three corner values,
//               renormalised. A linear interpolant evaluated at the centroid is
//               exactly that average, which is what makes the holonomy check
//               below well posed rather than a comparison of two conventions.
//   matching    per interior edge, the quarter turn between the two face values
//               -- what the quantizer reads.
//   singular
//   faces       mapped to the vertex of the face nearest its centroid.
//
// The report says what that cost. `holonomyLostEdges` is the number the outline
// calls the key one: the P1 field is continuous, so its angle transports
// unambiguously from one centroid to the next across a shared edge; the
// resampled field's jump is that transport reduced to (-pi, pi]. Where they
// differ by a full spin turn, the conversion has lost a quarter turn of the
// cross on that edge and the matching handed to the quantizer is not the P1
// field's holonomy.
FieldRun convertP1ToFaces(const Mesh &m, const Eigen::VectorXcd &uVertex,
                          const std::vector<std::pair<int, double>> &p1SingularFaces);

// --- output for the figures ------------------------------------------------

// A legacy VTK unstructured grid: the triangulation, the two cross axes per
// cell, the material id per cell, and the singularity index per point.
bool writeFieldVTK(const std::string &path, const Mesh &m, const Eigen::VectorXcd &u);

} // namespace paper

#endif // __PAPER_METHODS_HXX__
