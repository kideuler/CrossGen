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
// the field through `DualMBO::u_k_prev` and through nothing else. So the
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
    int dualMBOMaxSteps = 500;
    double dualMBOGamma = 10.0;
    // Seconds after which a model is abandoned. Zero is no limit. (Advisory:
    // it is checked between stages, so a stage that runs long overruns it.)
    double timeLimit = 0.0;

    // Named overrides of TORSION::Options fields, applied after the ones above,
    // so that a pipeline knob can be swept from an experiment's command line
    // (`--layout-opt outerSteps=32`) without this header knowing the pipeline.
    // The names are the TORSION::Options member names; an unknown one throws,
    // so a misspelt sweep fails rather than silently runs the default. The
    // recognised set is the list in Layout.cxx.
    std::vector<std::pair<std::string, double>> overrides;
};

struct LayoutResult {
    bool ran = false;
    std::string error;

    // Definition 2.1's verdict, which is what TORSION::run() returns.
    bool valid = false;
    // The last stage that completed, as a short name, so a table can say where
    // a model stopped rather than only that it did.
    std::string reachedStage = "none";

    // Whether the field this call passed in is the field the pipeline used.
    // False means Stage 0c changed the triangle count under it -- disk
    // templates on a domain with a circular inclusion -- and the pipeline
    // solved its own instead, so the row is not a comparison of fields.
    bool externalFieldUsed = false;

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

    // How integrable the field was, as the layout stage found it: the largest
    // and mean relative deviation, edge by edge, of the integrated map from the
    // field's own metric (TORSION Stage 4's fit residual). A field that is a
    // gradient integrates to a map that reproduces it exactly, so this is the
    // one field-quality number that is measured by the consumer rather than by
    // us, and it is zero only for a field with no non-integrable curl in it.
    // The largest frame jump of the combing is kept beside it.
    double integrationFitResidualMax = 0.0;
    double integrationFitResidualMean = 0.0;
    double integrationFlippedAreaFraction = 0.0;
    double maxFrameJump = 0.0;
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

// `field` empty runs TORSION's own DualMBO solve, which is the pipeline as it
// ships. Non-empty substitutes it. The mesh is deep-copied first, so the caller
// can run several fields against the same domain without one run's Stage 0c
// excision or circle refit reaching the next.
LayoutResult runLayout(const std::shared_ptr<Mesh> &m, const Eigen::VectorXcd &field,
                       const LayoutOptions &o, bool keepMesh = false);

bool writeQuadOBJ(const std::string &path, const LayoutResult &r);

// "name=value" -> an entry of LayoutOptions::overrides. Throws on a malformed
// string; the name itself is checked when the layout runs.
std::pair<std::string, double> parseLayoutOverride(const std::string &arg);

} // namespace paper

#endif // __PAPER_LAYOUT_HXX__
