#include "Layout.hxx"

#include <chrono>
#include <fstream>

#include "MERIDIAN/DiskTemplate.hxx"
#include "MERIDIAN/QuadMesh.hxx"
#include "TORSION/TORSION.hxx"

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
        opts.sipgGamma = o.sipgGamma;
        opts.sipgMaxSteps = o.sipgMaxSteps;
        opts.diskTemplates = o.diskTemplates;
        opts.topoNearMiss = o.topoNearMiss;
        opts.topoNearMissRetry = o.topoNearMissRetry;
        opts.quadTargetEdge = o.quadTargetEdge;
        if (field.size() > 0) opts.externalField = field;

        TORSION pipe(work, opts);
        r.valid = pipe.run();
        r.ran = true;

        const TORSION::Status &s = pipe.getStatus();
        r.reachedStage = reachedStage(s);
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
            if (pipe.hasDiskTemplate()) {
                const DiskTemplate &dt = pipe.getDiskTemplate();
                if (!dt.quads().empty()) {
                    verts = &dt.vertices();
                    cells = &dt.quads();
                    mats = &dt.quadMaterials();
                }
            }
            if (!cells->empty())
                r.quality = metrics::quadMetrics(*verts, *cells, *mats, o.quadTargetEdge);
            if (keepMesh) {
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
