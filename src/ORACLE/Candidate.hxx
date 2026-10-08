#ifndef __ORACLE_CANDIDATE_HXX__
#define __ORACLE_CANDIDATE_HXX__

#include <memory>
#include <string>
#include <vector>

#include "ATLAS/BlockMesh.hxx"
#include "MERIDIAN/DiskTemplate.hxx"
#include "MERIDIAN/QuadMesh.hxx"
#include "mesh/BlockDecomposition.hxx"
#include "mesh/BlockQuadMesh.hxx"
#include "mesh/Mesh.hxx"
#include "mesh/QuadMesh.hxx"
#include "mesh/TMOP.hxx"

namespace oracle {

// The five block-decomposition methods ORACLE chooses between, under the names
// py/build_dataset.py writes their columns with.
enum class Method { ZIPLINE, UMBER, MERIDIAN, TORSION, ATLAS };
constexpr int kNumMethods = 5;

// "zipline", "umber", "meridian", "torsion", "atlas".
const char *methodName(Method m);
// The method a column prefix names; false for none of the five.
bool methodNamed(const std::string &name, Method &out);

// A quadrilateral mesh on a Candidate's blocks, made by the method's own
// mesher, so exactly one of the first three is set: BlockQuadMesh for ZIPLINE
// and UMBER, BlockMesh (each block's certified chart) for ATLAS, MERIDIAN's
// Stage 10 on the fitted spline patches for MERIDIAN and TORSION. Each is kept
// as its own class rather than only as the mesh::QuadMesh the smoother takes,
// because the block walls a picture of it draws are that class's blocks. A
// Stage 10 mesh holds a reference into its pipeline, so none of this may
// outlive the Candidate that made it.
struct CandidateMesh {
    std::unique_ptr<BlockQuadMesh> blockQuadMesh;   // zipline, umber
    std::unique_ptr<BlockMesh> blockMesh;           // atlas
    std::unique_ptr<::QuadMesh> stage10;            // meridian, torsion
    // Stage 11's merged mesh where Stage 0c excised a disk -- only with the
    // pipelines' diskTemplates option, which is off by default.
    std::unique_ptr<DiskTemplate> stage11;
    // Nodes the smoother has to hold: ZIPLINE's open sides, under its
    // Options::pinOpenSides.
    std::vector<int> pinned;

    // The mesh as the smoother takes it -- Stage 11's merged one where there
    // is one -- with node options `o` and `pinned` held.
    mesh::QuadMesh quadMesh(const mesh::QuadMesh::Options &o) const;
};

// One of the five, run on one model the way crossgen's Mesh.<method>() runs
// it with no keywords -- the way py/build_dataset.py ran it for the selector's
// training data, so that what the selector predicts is the run this makes:
//
//   zipline   ZIPLINE::run(Stage::Blocks), its LayoutBlocks as the blocks
//   umber     the UMBER stages chained as TestUMBER chains them
//   meridian  MERIDIAN::run() through Stage 9, the arrangement's patches
//   torsion   TORSION::run() through Stage 9, the same
//   atlas     ATLAS::run(), the chosen search's certified cover
//
// with every Options at its default (copy 0 of the dataset's: the default
// seeds), the fitted sides of MERIDIAN and TORSION exported at 64 chords as
// src/python/MethodLayout.cxx exports them, and "valid" and the coverage read
// as crossgen reads them (Methods.hxx's coverage(), BlockDecomposition.valid).
//
// ### Why the runs are copied here and not shared with src/python
//
// The module's Method classes are these same runs, but written against
// Python.h: their options come in as keywords and their reports go out as
// dicts. So the two places that rebuild what no class runs on its own --
// UMBER's chain, which has no driver class, and Stages 10 and 11 of the two
// pipelines, which their run() only builds at one target edge length -- are
// here a second time, in the same calls and order as MethodUMBER.cxx and
// MethodLayout.cxx. If one of those changes, or TestUMBER or MERIDIAN::run()
// does, this has to follow it, or the selector is predicting runs ORACLE no
// longer makes.
class Candidate {
public:
    // Runs `method` on a copy of `model`. Never throws: a method that throws
    // -- and a ZIPLINE that leaves no blocks object, which crossgen raises on
    // -- is a run that failed, as the dataset counted it: raised(), with
    // coverage 0, no blocks, and error() saying why.
    static std::unique_ptr<Candidate> run(Method method, const Mesh &model);

    virtual ~Candidate() = default;

    Method method() const { return method_; }
    const char *name() const { return methodName(method_); }

    // The blocks: empty when the run failed outright.
    const BlockDecomposition &decomposition() const { return decomposition_; }
    // Fraction of the model's area inside a block, by the method's own account.
    double coverage() const { return coverage_; }
    // A decomposition of the model: every side on dS or shared, and coverage
    // of at least 99.9% -- crossgen's BlockDecomposition.valid.
    bool valid() const { return coverage_ >= 0.999 && decomposition_.covers(); }
    bool raised() const { return raised_; }
    const std::string &error() const { return error_; }
    double seconds() const { return seconds_; }

    // What the mesher is asked for; every other setting is the method's own
    // default, as crossgen's BlockDecomposition.mesh(h) leaves them.
    struct MeshSettings {
        double target = 0.05;     // target edge length, the model's units
        int minIntervals = 1;     // fewest edges on a chord
        int maxIntervals = 0;     // most; 0 for no ceiling
    };
    // A mesh on the blocks by the method's own mesher. Throws
    // std::runtime_error when the run left nothing to mesh or the mesher
    // made no element.
    virtual std::unique_ptr<CandidateMesh> mesh(const MeshSettings &s) const = 0;

    // The TMOP settings the method smooths its mesh with, and whether it
    // pillows the flat feature corners first (Stage 12 of MERIDIAN and
    // TORSION) -- crossgen's QuadMesh.smooth() defaults for it.
    virtual mesh::TMOP::Options smoothing() const = 0;
    virtual bool pillows() const { return false; }

protected:
    explicit Candidate(Method m) : method_(m) {}
    // The method itself; run() times it and catches what it throws.
    virtual void execute() = 0;

    Method method_;
    BlockDecomposition decomposition_;
    double coverage_ = 0.0;
    double seconds_ = 0.0;
    bool raised_ = false;
    std::string error_;
};

}  // namespace oracle

#endif // __ORACLE_CANDIDATE_HXX__
