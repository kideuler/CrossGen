#ifndef __ZIPLINE_HXX__
#define __ZIPLINE_HXX__

#include <memory>
#include <string>

#include "MERIDIAN/Interfaces.hxx"
#include "ZIPLINE/ConformalField.hxx"
#include "ZIPLINE/LayoutBlocks.hxx"
#include "ZIPLINE/PartitionSimplify.hxx"
#include "ZIPLINE/QuadLayout.hxx"
#include "ZIPLINE/SeparatrixTrace.hxx"
#include "crossfield/CrossField.hxx"
#include "mesh/BlockQuadMesh.hxx"
#include "mesh/Mesh.hxx"
#include "mesh/QuadMesh.hxx"
#include "mesh/TMOP.hxx"

// ZIPLINE: the separatrix partition of
//   Viertel, Osting and Staten, "Coarse quad layouts through robust simplification
//   of cross field separatrix partitions", IMR 2019   (docs/viertel_2019.md)
// taken on to a block decomposition, a mesh and a smoothed mesh the way every
// other method in this codebase is finished.
//
//   Stage 0b  Interfaces          MERIDIAN's material interface network, closed
//                                 loops unsplit (multi-material models only)
//   Stage 1   CrossField          the P1 MBO cross field (baseline B1), aligned
//                                 to the interfaces
//   Stage 1b  ConformalField      the corrected field of docs/
//                                 quad_layout_pipeline.md Sec. 6: Stage 1's
//                                 singularity pattern kept, its field replaced
//                                 by the conformal one with exactly those
//                                 singularities (two linear problems, P2 Sec. 5)
//   Stage 2   SeparatrixTrace     singularities, their ports, the separatrices
//                                 traced out of them (paper Sec. 3); on the
//                                 corrected field also joined to a singular
//                                 point whose critical patch they enter, one
//                                 copy of each such connection kept (docs/
//                                 quad_layout_pipeline.md Secs. 7.5 and 8)
//   Stage 3   QuadLayout          the T-layout the separatrices cut the model into
//   Stage 4   PartitionSimplify   chord collapses (Sec. 4), then StemExtension
//                                 (spec Sec. 12), then Sec. 4 again
//   Stage 5   LayoutBlocks        the four-sided components as the shared
//                                 BlockDecomposition on spline geometry
//   Stage 6   BlockQuadMesh       a transfinite grid on the blocks
//   Stage 7   mesh::TMOP          the grid smoothed, mu sampled at the corners
//
// docs/viertel_2019.md Sec. 14 maps each section of the spec to the stage that
// implements it, and records every place the code departs from the paper.
//
// ### Why one class
//
// Until 2026-09-24 TestZIPLINE and the viewer's mode 2 each ran these stages
// themselves, and they had drifted apart where nobody was looking: the viewer
// started the MBO solve from CrossField::initialize(1) with no seed -- which
// initialize() reads as "seed from the clock", so no two sessions traced the
// same field -- and stopped it at a tolerance a hundred times looser than the
// driver's, and a 'c' pressed during the MBO animation traced a field that had
// not finished converging. A layout the corpus sweep reported was therefore
// not the layout the viewer showed, and a defect seen in the viewer could not
// be reproduced from the command line. Both now run this class with these
// Options, so the same model gives the same field, the same layout and the
// same blocks in both.
//
// ### Running it
//
// run() takes the model as far as the stage it is given. The stages can also
// be run one at a time, for a caller that wants to look in between -- the
// viewer animates Stages 1 and 2 a few steps per frame -- and each one first
// runs whatever it needs that has not run. Running a stage again rebuilds it
// from the stage before and discards every stage after it: re-meshing at
// another target edge length throws away the smoothed mesh, but not the
// blocks.
//
// Nothing here decides whether a result is good. Whether the layout is sound,
// how much of the model the blocks cover and how good the elements are is
// read off the stages' own reports by the caller, the way TestMERIDIAN reads
// MERIDIAN's.
class ZIPLINE {
public:
    enum class Stage { Interfaces, Field, Conformal, Trace, Layout, Simplified, Blocks, Mesh, Smoothed };

    struct Options {
        // ---- Stage 0b ----
        // Find the material interface network and honour it: the field is
        // aligned to it, the tracing emits from its nodes and cuts its
        // branches, and Sec. 4 never moves it. No effect on a single-material
        // model, which has no network.
        bool materialInterfaces = true;
        // Pin each disk inclusion's centre vertex to u = 1 (DualMBO's pin),
        // which puts the inclusion's four cones on its diagonals and makes it
        // an O-grid. Only with interfaces.
        bool pinDiskCenters = true;

        // ---- Stage 1 ----
        // The MBO solve stops at error < 2 N fieldTolerance, N the vertex
        // count, or after fieldMaxSteps steps. 1e-9 is what every corpus
        // number in docs/viertel_2019.md Sec. 14 was measured at.
        int fieldMaxSteps = 500;
        double fieldTolerance = 1e-9;
        // CrossField::initialize(1, seed) starts from a random cross per
        // vertex. Never 0: initialize() takes 0 to mean "seed from the clock".
        unsigned fieldSeed = 12345;

        // ---- Stage 1b ----
        // Trace the corrected field (ConformalField) instead of the MBO one:
        // the MBO's singularity pattern, its field rebuilt as the conformal
        // field with exactly those singularities. Falls back to the MBO field,
        // with the reason in Status::conformalFailure, on a multi-material
        // model (the corrected field knows no interfaces) and wherever the
        // pattern fails Poincare-Hopf (the Neumann problem has no solution).
        bool secondField = true;
        ConformalField::Options conformal;

        // ---- Stages 2 and 4 ----
        // Viertel's tracing, with docs/quad_layout_pipeline.md Sec. 7.5's
        // capture by critical patches and Sec. 8's reduction on top of it
        // (SeparatrixTrace::Settings::captureAtSingularities has the
        // measurements). Applied only when the trace follows Stage 1b's
        // corrected field: the capture trusts a separatrix to arrive along
        // a port when it arrives at all, which the first field's cores do not
        // let it do -- on the multi-material corpus, traced on the first
        // field, it took T-junctions left from 43 to 151 and three layouts
        // unsound.
        SeparatrixTrace::Settings trace = [] {
            SeparatrixTrace::Settings t;
            t.captureAtSingularities = true;
            return t;
        }();
        // Off: the layout goes on to the blocks exactly as it was traced.
        bool simplify = true;
        PartitionSimplify::Settings simplification;

        // ---- Stage 6 ----
        BlockQuadMesh::Options mesh;

        // ---- Stage 7 ----
        // Metric 7 at the element corners. Corners because a transfinite grid
        // on a block decomposition can turn a corner over without the 2x2
        // Gauss points seeing it (ATLAS's and UMBER's setting, for that
        // reason).
        mesh::TMOP::Options tmop = [] {
            mesh::TMOP::Options t;
            t.metric = mesh::TMOP::ShapeSize007;
            t.quadrature = mesh::TMOP::Corners;
            t.maxSweeps = 1000;
            return t;
        }();
        // Hold the rim of every component that got no blocks, where a mesh of
        // it would one day have to meet this one. Off because it costs the
        // smoother on every model where it differs (at h = 0.05, 1000 sweeps:
        // geom031 worst scaled Jacobian 0.32 pinned against 0.74 sliding,
        // geom010 0.16 against 0.29, geom007 0.69 against 0.86): the ring of
        // elements along a pinned rim cannot redistribute and takes the whole
        // mismatch.
        bool pinOpenSides = false;
    };

    struct Status {
        // Stage 1
        int fieldSteps = 0;
        bool fieldConverged = false;   // under the tolerance, not out of steps
        double fieldError = 0.0;
        // Stage 1b: whether the trace follows the corrected field, and if it
        // was asked for and does not, why not.
        bool conformalTraced = false;
        std::string conformalFailure;
    };

    explicit ZIPLINE(std::shared_ptr<Mesh> mesh);
    ZIPLINE(std::shared_ptr<Mesh> mesh, const Options &opts);
    ~ZIPLINE();

    // Every stage up to and including `last`. False only when a stage had
    // nothing to work on: no block to mesh, or no element to smooth.
    bool run(Stage last = Stage::Blocks);

    // ---- one stage at a time ----
    void findInterfaces();                    // Stage 0b
    void startField();                        // Stage 1: the initial field
    // Up to `steps` more MBO steps, then the field's singular triangles for
    // drawing. True once the field is done (converged or out of steps).
    bool stepField(int steps);
    void finishField();
    void buildConformalField();               // Stage 1b
    void startTrace();                        // Stage 2: ports, nothing grown
    // One round of growth. True once every separatrix has stopped.
    bool stepTrace();
    void finishTrace();
    void buildLayout();                       // Stage 3
    void simplifyLayout();                    // Stage 4
    void buildBlocks();                       // Stage 5
    // Stage 6. False when there is no block to mesh.
    bool buildMesh();
    bool buildMesh(const BlockQuadMesh::Options &opts);
    // Stage 7. False when there is no element to smooth.
    bool smooth();

    // ---- what the stages left ----
    // Null until its stage has run. The interfaces stay null on a
    // single-material model, and with Options::materialInterfaces off.
    bool hasInterfaces() const { return interfaces != nullptr; }
    bool hasField() const { return field != nullptr; }
    bool fieldFinished() const { return fieldDone; }
    // Stage 1b ran and built a field (it may have run and declined).
    bool hasConformalField() const { return conformal && conformal->valid(); }
    bool hasTrace() const { return trace != nullptr; }
    bool traceFinished() const { return trace && trace->finishedTracing; }
    bool hasLayout() const { return layout != nullptr; }
    bool hasSimplification() const { return simplification != nullptr; }
    bool hasBlocks() const { return blocks != nullptr; }
    bool hasBlockMesh() const { return blockMesh != nullptr; }
    bool hasSmoothedMesh() const { return smoothedMesh != nullptr; }

    const Interfaces &getInterfaces() const { return *interfaces; }
    const CrossField &getField() const { return *field; }
    // Stage 1b's object, built or not (its report says which, and why).
    // Null until the stage has run, and when it was not asked for.
    const ConformalField *getConformalField() const { return conformal.get(); }
    // The field the separatrices follow: Stage 1b's when it built, Stage 1's
    // otherwise. Only meaningful once the trace exists.
    const CrossField &getTraceField() const {
        return (trace && status.conformalTraced) ? *conformal->crossField() : *field;
    }
    const SeparatrixTrace &getTrace() const { return *trace; }
    // The layout as traced, before Sec. 4.
    const QuadLayout &getLayout() const { return *layout; }
    // Sec. 4's pass: its report, its chords, and the layout it left.
    const PartitionSimplify &getSimplification() const { return *simplification; }
    const QuadLayout &getSimplifiedLayout() const { return simplification->getLayout(); }
    const LayoutBlocks &getBlocks() const { return *blocks; }
    const BlockQuadMesh &getBlockMesh() const { return *blockMesh; }
    // Stage 7's mesh, with the nodes where the smoother left them.
    const mesh::QuadMesh &getSmoothedMesh() const { return *smoothedMesh; }
    const mesh::TMOP::Report &getTMOPReport() const { return tmopReport; }

    const Options &getOptions() const { return options; }
    const Status &getStatus() const { return status; }
    const Mesh &getMesh() const { return *mesh; }
    std::shared_ptr<Mesh> getMeshPtr() const { return mesh; }

private:
    // Drop `from` and every stage after it.
    void discardFrom(Stage from);
    // Stage 1's steps without the singularities: finishField() computes them
    // once, stepField() after every call.
    void advanceField(int steps);

    std::shared_ptr<Mesh> mesh;
    Options options;
    Status status;

    bool interfacesFound = false;
    bool fieldDone = false;
    bool conformalTried = false;
    // Each holds on to the one before it, so they are declared in stage order
    // and destroyed in reverse.
    std::unique_ptr<Interfaces> interfaces;
    std::shared_ptr<CrossField> field;
    std::shared_ptr<ConformalField> conformal;
    std::unique_ptr<SeparatrixTrace> trace;
    std::unique_ptr<QuadLayout> layout;
    std::unique_ptr<PartitionSimplify> simplification;
    std::unique_ptr<LayoutBlocks> blocks;
    std::unique_ptr<BlockQuadMesh> blockMesh;
    std::unique_ptr<mesh::QuadMesh> smoothedMesh;
    mesh::TMOP::Report tmopReport;
};

#endif // __ZIPLINE_HXX__
