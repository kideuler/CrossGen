#include "Layout.hxx"

#include <chrono>
#include <cmath>
#include <fstream>
#include <stdexcept>

#include "MERIDIAN/DiskTemplate.hxx"
#include "MERIDIAN/QuadMesh.hxx"
#include "TORSION/TORSION.hxx"
#include "mesh/QuadMesh.hxx"
#include "mesh/TMOP.hxx"

namespace paper {
namespace {

using Clock = std::chrono::steady_clock;

// The stage a run got to, read off the Status flags in pipeline order. This is
// the column E4's table calls "reached", and it is what makes a failure
// reportable rather than merely counted.
std::string reachedStage(const TORSION::Status &s) {
    if (s.meshValid)         return "mesh";
    if (s.quadMeshRan)       return "mesh(partial)";
    if (s.splinesValid)      return "splines";
    if (s.arrangementValid)  return "arrangement";
    if (s.separatricesRan)   return "separatrices";
    if (s.layoutValid)       return "layout";
    if (s.layoutRan)         return "layout(partial)";
    if (s.immersionValid)    return "psi_0";
    if (s.framesValid)       return "frames";
    if (s.cutIsDisk)         return "cut";
    if (s.conesAdmissible)   return "cones";
    return "field";
}

// The TORSION::Options fields an experiment may override by name. Booleans take
// 0/1, integers are rounded.
void applyOverrides(TORSION::Options &t, const std::vector<std::pair<std::string, double>> &ov) {
    for (const auto &[key, val] : ov) {
        const int iv = static_cast<int>(std::lround(val));
        const bool bv = val != 0.0;
        if      (key == "targetEdge")               t.targetEdge = val;
        else if (key == "conformalSizing")          t.conformalSizing = bv;
        else if (key == "alignInIntegration")       t.alignInIntegration = bv;
        else if (key == "alignAcrossSeams")         t.alignAcrossSeams = bv;
        else if (key == "alignmentChooseByLadder")  t.alignmentChooseByLadder = bv;
        else if (key == "reconcileSectors")         t.reconcileSectors = bv;
        else if (key == "alignmentReleaseRounds")   t.alignmentReleaseRounds = iv;
        else if (key == "pullOntoAlignment")        t.pullOntoAlignment = bv;
        else if (key == "regularisedUntangleRings") t.regularisedUntangleRings = iv;
        else if (key == "retryKeepingDipoles")      t.retryKeepingDipoles = bv;
        else if (key == "topoNearMissLastRetry")    t.topoNearMissLastRetry = val;
        else if (key == "topoRetryUnseeded")        t.topoRetryUnseeded = bv;
        else if (key == "perMaterial")              t.perMaterial = bv;
        else if (key == "perMaterialRounds")        t.perMaterialRounds = iv;
        else if (key == "perMaterialMatchTolerance") t.perMaterialMatchTolerance = val;
        else if (key == "perMaterialFallback")      t.perMaterialFallback = bv;
        else if (key == "perMaterialSpacedMatching") t.perMaterialSpacedMatching = bv;
        else if (key == "perMaterialBendEnds")      t.perMaterialBendEnds = bv;
        else if (key == "quadMaterialsFromFaces")   t.quadMaterialsFromFaces = bv;
        else if (key == "quadContractOntoFeatures") t.quadContractOntoFeatures = bv;
        else if (key == "untangle")                 t.untangle = bv;
        else if (key == "localUntangle")            t.localUntangle = bv;
        else if (key == "lambdaInit")               t.lambdaInit = val;
        else if (key == "lambdaGrowth")             t.lambdaGrowth = val;
        else if (key == "lambdaAlignmentFactor")    t.lambdaAlignmentFactor = val;
        else if (key == "lambdaSeamFactor")         t.lambdaSeamFactor = val;
        else if (key == "outerSteps")               t.outerSteps = iv;
        else if (key == "innerIterations")          t.innerIterations = iv;
        else if (key == "relabelBetweenSteps")      t.relabelBetweenSteps = bv;
        else if (key == "seedTopoConstraints")      t.seedTopoConstraints = bv;
        else if (key == "topoNearMiss")             t.topoNearMiss = val;
        else if (key == "topoNearMissRetry")        t.topoNearMissRetry = val;
        else if (key == "seedSelfReturns")          t.seedSelfReturns = bv;
        else if (key == "seedAllConnections")       t.seedAllConnections = bv;
        else if (key == "separatrixMaxSteps")       t.separatrixMaxSteps = iv;
        else if (key == "separatrixSnapRings")      t.separatrixSnapRings = iv;
        else if (key == "separatrixDetectCycles")   t.separatrixDetectCycles = bv;
        else if (key == "arrangementCorner")        t.arrangementCorner = val;
        else if (key == "arrangementCollapse")      t.arrangementCollapse = val;
        else if (key == "arrangementTrim")          t.arrangementTrim = bv;
        else if (key == "quadMinIntervals")         t.quadMinIntervals = iv;
        else if (key == "quadMaxIntervals")         t.quadMaxIntervals = iv;
        else if (key == "quadSmoothingPasses")      t.quadSmoothingPasses = iv;
        else if (key == "repairPasses")             t.repairPasses = iv;
        else if (key == "repairGapLimit")           t.repairGapLimit = val;
        else if (key == "repairLambdaBoost")        t.repairLambdaBoost = val;
        else if (key == "repairOuterSteps")         t.repairOuterSteps = iv;
        else if (key == "repairMaxPerPass")         t.repairMaxPerPass = iv;
        else if (key == "repairPatience")           t.repairPatience = iv;
        else if (key == "interfaceCorners")         t.interfaceCorners = bv;
        else if (key == "cancelInterfaceDipoles")   t.cancelInterfaceDipoles = bv;
        else if (key == "prescribeInterfaceCones")  t.prescribeInterfaceCones = bv;
        else if (key == "coneCutsToBoundary")       t.coneCutsToBoundary = bv;
        else throw std::invalid_argument("LayoutOptions::overrides: unknown TORSION option '" + key + "'");
    }
}

} // namespace

LayoutResult runLayout(const std::shared_ptr<Mesh> &m, const Eigen::VectorXcd &field,
                       const LayoutOptions &o, bool keepMesh) {
    LayoutResult r;
    const Clock::time_point t0 = Clock::now();

    // A copy, so that Stage 0c's excision and the circle refit stay inside this
    // run and the caller's mesh is the same object for the next field.
    std::shared_ptr<Mesh> work =
        std::make_shared<Mesh>(m->vertices, m->triangles, m->triangleMatId);

    try {
        TORSION::Options opts;
        opts.dualMBOGamma = o.dualMBOGamma;
        opts.dualMBOMaxSteps = o.dualMBOMaxSteps;
        opts.diskTemplates = o.diskTemplates;
        opts.topoNearMiss = o.topoNearMiss;
        opts.topoNearMissRetry = o.topoNearMissRetry;
        opts.quadTargetEdge = o.quadTargetEdge;
        applyOverrides(opts, o.overrides);
        if (field.size() > 0) opts.externalField = field;

        TORSION pipe(work, opts);
        r.valid = pipe.run();
        r.ran = true;

        const TORSION::Status &s = pipe.getStatus();
        r.reachedStage = reachedStage(s);
        r.externalFieldUsed = s.externalFieldUsed;
        r.framesValid = s.framesValid;
        r.immersionValid = s.immersionValid;
        r.layoutValid = s.layoutValid;
        r.separatricesRan = s.separatricesRan;
        r.arrangementValid = s.arrangementValid;
        r.splinesValid = s.splinesValid;
        r.meshValid = s.meshValid;
        r.meshConforming = s.meshConforming;

        r.interiorCones = s.interiorCones;
        r.boundaryCones = s.boundaryCones;
        r.integrationFlippedFaces = s.integrationFlippedFaces;
        r.integrationFitResidualMax = s.integrationMaxFitResidual;
        r.integrationFitResidualMean = s.integrationMeanFitResidual;
        r.integrationFlippedAreaFraction = s.integrationFlippedAreaFraction;
        r.maxFrameJump = s.maxFrameJump;
        r.separatrices = s.separatrices;
        r.separatricesUnresolved = s.separatricesUnresolved;
        r.patches = s.layoutPatches;
        r.coverage = s.layoutCoverage;
        r.meshVertices = s.meshVertices;
        r.meshQuads = s.meshQuads;
        r.unmeshedPatches = s.meshUnmeshedPatches;
        r.pipelineMinScaledJacobian = s.meshMinScaledJacobian;
        r.outerSteps = s.outerStepsTaken;
        r.messages = s.messages;

        if (pipe.hasQuadMesh()) {
            const ::QuadMesh &qm = pipe.getQuadMesh();
            const ::QuadMesh::Report &qr = qm.getReport();
            r.pipelineMixedQuads = qr.mixedQuads;
            r.pipelineMaterials = qr.materials;

            // Stage 11 puts the disk templates back and produces a mesh of its
            // own; that is the final mesh when it ran.
            const std::vector<Point> *verts = &qm.vertices();
            const std::vector<std::array<int, 4>> *cells = &qm.quads();
            const std::vector<int> *mats = &qm.quadMaterials();
            bool usingDiskTemplate = false;
            if (pipe.hasDiskTemplate()) {
                const DiskTemplate &dt = pipe.getDiskTemplate();
                if (!dt.quads().empty()) {
                    verts = &dt.vertices();
                    cells = &dt.quads();
                    mats = &dt.quadMaterials();
                    usingDiskTemplate = true;
                }
            }
            if (!cells->empty()) {
                r.quality = metrics::quadMetrics(*verts, *cells, *mats, o.quadTargetEdge);

                if (o.tmopSweeps > 0) {
                    // Same source as verts/cells/mats above: the disk-template
                    // merged mesh when Stage 11 produced one, otherwise Stage 10's.
                    mesh::QuadMesh fm = usingDiskTemplate
                                            ? mesh::QuadMesh::from(pipe.getDiskTemplate())
                                            : mesh::QuadMesh::from(qm);
                    mesh::TMOP::Options topt;
                    topt.metric = mesh::TMOP::ShapeSize007;
                    topt.maxSweeps = o.tmopSweeps;
                    topt.exponent = o.tmopPower;
                    mesh::TMOP smoother(fm, topt);
                    r.tmopOk = smoother.run();
                    const mesh::TMOP::Report &tr = smoother.getReport();
                    r.tmopConverged = tr.converged;
                    r.tmopSweepsRun = tr.sweeps;
                    r.tmopSeconds = tr.seconds;
                    r.qualitySmoothed =
                        metrics::quadMetrics(fm.vertices, fm.quads, fm.quadMatId, o.quadTargetEdge);
                    r.smoothed = true;
                    if (keepMesh) {
                        // The smoothed positions are the ones worth keeping for a
                        // figure; the topology is exactly what verts/cells above
                        // already carried.
                        r.quadVertices = fm.vertices;
                        r.quadCells = fm.quads;
                        r.quadMaterials = fm.quadMatId;
                    }
                }
            }
            if (keepMesh && !r.smoothed) {
                r.quadVertices = *verts;
                r.quadCells = *cells;
                r.quadMaterials = *mats;
            }
        }
    } catch (const std::exception &e) {
        r.error = e.what();
    } catch (...) {
        r.error = "unknown exception";
    }

    r.seconds = std::chrono::duration<double>(Clock::now() - t0).count();
    return r;
}

std::pair<std::string, double> parseLayoutOverride(const std::string &arg) {
    const std::size_t eq = arg.find('=');
    if (eq == std::string::npos || eq == 0 || eq + 1 >= arg.size())
        throw std::invalid_argument("--layout-opt expects NAME=VALUE, got '" + arg + "'");
    return {arg.substr(0, eq), std::stod(arg.substr(eq + 1))};
}

bool writeQuadOBJ(const std::string &path, const LayoutResult &r) {
    if (r.quadCells.empty()) return false;
    std::ofstream f(path);
    if (!f) return false;
    f << "# quad mesh from the paper_tests layout run\n";
    for (const Point &p : r.quadVertices) f << "v " << p[0] << " " << p[1] << " 0\n";
    int lastMat = -1;
    for (std::size_t c = 0; c < r.quadCells.size(); ++c) {
        const int mid = c < r.quadMaterials.size() ? r.quadMaterials[c] : 1;
        if (mid != lastMat) { f << "usemtl mat" << mid << "\n"; lastMat = mid; }
        f << "f " << r.quadCells[c][0] + 1 << " " << r.quadCells[c][1] + 1 << " "
          << r.quadCells[c][2] + 1 << " " << r.quadCells[c][3] + 1 << "\n";
    }
    return true;
}

} // namespace paper
