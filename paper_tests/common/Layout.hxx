#ifndef __PAPER_LAYOUT_HXX__
#define __PAPER_LAYOUT_HXX__

#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "Metrics.hxx"
#include "mesh/Mesh.hxx"

// E4 and E5's end of the pipeline: a cross field in, a quadrilateral mesh out,
// through TORSION.
//
// The experiment is a controlled one and this is where the control lives.
// TORSION runs Stages 0b to 11 -- the interface network, the cone set, the cut,
// the frames, the integration, the layout energies, the separatrices, the
// arrangement, the spline fit and the meshing -- and every one of them reads
// the field through `SIPG::u_k_prev` and through nothing else. So the
// experiment hands TORSION a field through Options::externalField and changes
// nothing else at all: same mesh, same target edge length, same seeding
// tolerance, same everything. A difference in the layout is then a difference
// the field made.
//
// Two honest caveats, both of which E4 reports rather than hides:
//
//   1. Stages 1 and 2 do not take the field's singularities as given. Stage 1
//      prescribes cones on the interface network, cancels dipoles inside a
//      region and rebalances to satisfy Gauss-Bonnet, so a field whose indices
//      are wrong is partly repaired before the layout ever sees them. That
//      makes the comparison *conservative*: it understates the cost of a bad
//      field rather than overstating it.
//   2. Not every model reaches a mesh, on either field. The outline says so and
//      the tables report the stage each run stopped at, so "success rate" is a
//      measured column and not an assumption.
namespace paper {

struct LayoutOptions {
    double quadTargetEdge = 0.05;
    bool diskTemplates = false;
    double topoNearMiss = 0.15;
    double topoNearMissRetry = 0.04;
    int sipgMaxSteps = 500;
    double sipgGamma = 10.0;
    // Seconds after which a model is abandoned. Zero is no limit. (Advisory:
    // it is checked between stages, so a stage that runs long overruns it.)
    double timeLimit = 0.0;
};

struct LayoutResult {
    bool ran = false;
    std::string error;

    // Definition 2.1's verdict, which is what TORSION::run() returns.
    bool valid = false;
    // The last stage that completed, as a short name, so a table can say where
    // a model stopped rather than only that it did.
    std::string reachedStage = "none";

    bool framesValid = false;
    bool immersionValid = false;
    bool layoutValid = false;
    bool separatricesRan = false;
    bool arrangementValid = false;
    bool splinesValid = false;
    bool meshValid = false;
    bool meshConforming = false;

    int interiorCones = 0;
    int boundaryCones = 0;
    int integrationFlippedFaces = 0;
    int separatrices = 0;
    int separatricesUnresolved = 0;
    int patches = 0;
    double coverage = 0.0;
    int meshVertices = 0;
    int meshQuads = 0;
    int unmeshedPatches = 0;

    // The pipeline's own reading of the final mesh, kept alongside ours because
    // they are computed differently and a disagreement is worth seeing.
    double pipelineMinScaledJacobian = 0.0;
    int pipelineMixedQuads = 0;
    int pipelineMaterials = 0;

    // Ours, from the element table alone (metrics::quadMetrics).
    metrics::QuadMetrics quality;

    double seconds = 0.0;
    int outerSteps = 0;
    std::vector<std::string> messages;

    // The final mesh, kept when `keepMesh` was asked for, so a figure can be
    // written without running the pipeline twice.
    std::vector<Point> quadVertices;
    std::vector<std::array<int, 4>> quadCells;
    std::vector<int> quadMaterials;
};

// `field` empty runs TORSION's own SIPG solve, which is the pipeline as it
// ships. Non-empty substitutes it. The mesh is deep-copied first, so the caller
// can run several fields against the same domain without one run's Stage 0c
// excision or circle refit reaching the next.
LayoutResult runLayout(const std::shared_ptr<Mesh> &m, const Eigen::VectorXcd &field,
                       const LayoutOptions &o, bool keepMesh = false);

bool writeQuadOBJ(const std::string &path, const LayoutResult &r);

} // namespace paper

#endif // __PAPER_LAYOUT_HXX__
