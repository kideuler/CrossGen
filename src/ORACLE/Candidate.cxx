#include "ORACLE/Candidate.hxx"

#include <chrono>
#include <cmath>
#include <stdexcept>

#include "ATLAS/ATLAS.hxx"
#include "MERIDIAN/Interfaces.hxx"
#include "MERIDIAN/MERIDIAN.hxx"
#include "Parameterization/HarmonicCut.hxx"
#include "TORSION/TORSION.hxx"
#include "UMBER/BlockLayout.hxx"
#include "UMBER/MotorcycleGraph.hxx"
#include "UMBER/Polysquare.hxx"
#include "UMBER/UMBER.hxx"
#include "ZIPLINE/ZIPLINE.hxx"
#include "dualmbo/DualMBO.hxx"

namespace oracle {

const char *methodName(Method m) {
    switch (m) {
        case Method::ZIPLINE:  return "zipline";
        case Method::UMBER:    return "umber";
        case Method::MERIDIAN: return "meridian";
        case Method::TORSION:  return "torsion";
        case Method::ATLAS:    return "atlas";
    }
    return "?";
}

bool methodNamed(const std::string &name, Method &out) {
    for (int i = 0; i < kNumMethods; ++i) {
        if (name == methodName(static_cast<Method>(i))) {
            out = static_cast<Method>(i);
            return true;
        }
    }
    return false;
}

mesh::QuadMesh CandidateMesh::quadMesh(const mesh::QuadMesh::Options &o) const {
    mesh::QuadMesh m;
    if (stage11) m = mesh::QuadMesh::from(*stage11, o);
    else if (stage10) m = mesh::QuadMesh::from(*stage10, o);
    else if (blockMesh) m = mesh::QuadMesh::from(*blockMesh, o);
    else if (blockQuadMesh) m = mesh::QuadMesh::from(*blockQuadMesh, o);
    else throw std::runtime_error("no mesh");
    // ZIPLINE::smooth()'s Options::pinOpenSides, on this mesh.
    if (!pinned.empty()) {
        for (const int v : pinned) m.pinVertex(v);
        m.classifyNodes();
        m.buildFeatureCurves();
        m.computeSlideTangents();
    }
    return m;
}

namespace {

// crossgen's TMOP settings for a transfinite grid on a block decomposition --
// UMBER's and ATLAS's (src/python/PyTables.cxx, blockGridSmoothing()).
mesh::TMOP::Options blockGridSmoothing() {
    mesh::TMOP::Options t;
    t.metric = mesh::TMOP::ShapeSize007;
    t.quadrature = mesh::TMOP::Corners;
    t.maxSweeps = 1000;
    return t;
}

// Fraction of `model`'s area inside D's blocks, each block's outline read as
// the polygon of its four sides: crossgen's coverage for ATLAS
// (src/python/Methods.cxx, decompositionCoverage()).
double decompositionCoverage(const BlockDecomposition &D, const Mesh &model) {
    double modelArea = 0.0;
    for (const Triangle &t : model.triangles) {
        const Point &a = model.vertices[t[0]], &b = model.vertices[t[1]], &c = model.vertices[t[2]];
        modelArea += 0.5 * std::fabs(cross2(b - a, c - a));
    }
    if (!(modelArea > 0.0)) return 0.0;
    double blockArea = 0.0;
    for (int b = 0; b < static_cast<int>(D.blocks.size()); ++b) {
        double twice = 0.0;
        for (int s = 0; s < 4; ++s) {
            const std::vector<Point> side = D.sidePolyline(b, s);
            for (size_t i = 0; i + 1 < side.size(); ++i) twice += cross2(side[i], side[i + 1]);
        }
        blockArea += 0.5 * std::fabs(twice);
    }
    return blockArea / modelArea;
}

void requireBlocks(const BlockDecomposition &D, const char *name) {
    if (D.blocks.empty()) throw std::runtime_error(std::string(name) + ": no block to mesh");
}

// ---- ZIPLINE (src/python/MethodZIPLINE.cxx) -----------------------------------
class ZiplineCandidate final : public Candidate {
public:
    explicit ZiplineCandidate(const Mesh &model)
        : Candidate(Method::ZIPLINE), zipline_(std::make_shared<Mesh>(model), ZIPLINE::Options()) {}

    std::unique_ptr<CandidateMesh> mesh(const MeshSettings &s) const override {
        requireBlocks(decomposition_, "zipline");
        BlockQuadMesh::Options mo = zipline_.getOptions().mesh;
        mo.targetEdgeLength = s.target;
        mo.minIntervals = s.minIntervals;
        mo.maxIntervals = s.maxIntervals;
        auto out = std::make_unique<CandidateMesh>();
        out->blockQuadMesh = std::make_unique<BlockQuadMesh>(decomposition_, mo);
        if (out->blockQuadMesh->getReport().quads <= 0)
            throw std::runtime_error("zipline: the mesh has no elements");
        if (zipline_.getOptions().pinOpenSides) out->pinned = out->blockQuadMesh->openSideVertices();
        return out;
    }

    mesh::TMOP::Options smoothing() const override { return zipline_.getOptions().tmop; }

private:
    void execute() override {
        zipline_.run(ZIPLINE::Stage::Blocks);
        if (!zipline_.hasBlocks()) throw std::runtime_error("Stage 5 produced no blocks object");
        decomposition_ = zipline_.getBlocks().decomposition();
        coverage_ = zipline_.getBlocks().coverage();
    }

    ZIPLINE zipline_;
};

// ---- UMBER (src/python/MethodUMBER.cxx, TestUMBER's chain) ---------------------
//
// The numbers are MethodUMBER.cxx's UMBEROptions defaults, which are copies of
// the stages' own member initialisers: TestUMBER's DualMBO steps and gamma,
// no harmonic cut of the voids, and the frame field's and the polysquare's
// weights and tolerances.
class UmberCandidate final : public Candidate {
public:
    explicit UmberCandidate(const Mesh &model)
        : Candidate(Method::UMBER), mesh_(std::make_shared<Mesh>(model)) {}

    std::unique_ptr<CandidateMesh> mesh(const MeshSettings &s) const override {
        requireBlocks(decomposition_, "umber");
        BlockQuadMesh::Options mo;
        mo.targetEdgeLength = s.target;
        mo.minIntervals = s.minIntervals;
        mo.maxIntervals = s.maxIntervals;
        auto out = std::make_unique<CandidateMesh>();
        out->blockQuadMesh = std::make_unique<BlockQuadMesh>(decomposition_, mo);
        if (out->blockQuadMesh->getReport().quads <= 0)
            throw std::runtime_error("umber: the mesh has no elements");
        return out;
    }

    mesh::TMOP::Options smoothing() const override { return blockGridSmoothing(); }

private:
    void execute() override {
        try {
            interfaces_ = std::make_unique<Interfaces>(mesh_);
        } catch (const std::exception &) {
            // TestUMBER's answer too: carry on as a single-material model.
            interfaces_.reset();
        }
        const bool multi = interfaces_ && interfaces_->multiMaterial();

        field_ = std::make_unique<DualMBO>(mesh_, 500, 10.0);
        if (multi) field_->setAlignedInteriorEdges(interfaces_->interfaceEdges());
        field_->initialize();
        field_->runMBO();
        field_->computeSingularities();

        cut_ = std::make_unique<HarmonicCut>(mesh_, false);

        umber_ = std::make_unique<UMBER>(*field_, *cut_);
        if (multi) umber_->setFeatureEdges(interfaces_->interfaceEdges());
        umber_->setMaxIterations(1000);
        umber_->setAlignmentWeight(0.1);
        umber_->setRegularizationWeight(5.0);
        umber_->setTolerance(1e-6);
        umber_->setCancelDipoles(true);
        umber_->setMaxCornerQuarters(1);
        umber_->initialize();
        umber_->optimize();

        poly_ = std::make_unique<Polysquare>(*umber_, *cut_);
        poly_->setInterfaceWeight(0.15);
        poly_->setCornerWeight(10.0);
        poly_->setSnapBoundary(true);
        poly_->setCornerSnapTolerance(0.05);
        poly_->setMaxIterations(1000);
        poly_->setTolerance(1e-6);
        poly_->solve();

        graph_ = std::make_unique<MotorcycleGraph>(*poly_);
        graph_->setArrivalTolerance(0.20);
        graph_->build();

        layout_ = std::make_unique<BlockLayout>(*graph_);
        layout_->build();
        BlockLayout::DecompositionReport report;
        decomposition_ = BlockLayout::blockDecompositionOf(layout_->getLayout(), *mesh_, &report);
        coverage_ = BlockLayout::coverageOf(report);
    }

    std::shared_ptr<Mesh> mesh_;
    // In stage order: each later stage holds references into the ones before
    // it, so they are destroyed in reverse.
    std::unique_ptr<Interfaces> interfaces_;
    std::unique_ptr<DualMBO> field_;
    std::unique_ptr<HarmonicCut> cut_;
    std::unique_ptr<UMBER> umber_;
    std::unique_ptr<Polysquare> poly_;
    std::unique_ptr<MotorcycleGraph> graph_;
    std::unique_ptr<BlockLayout> layout_;
};

// ---- MERIDIAN and TORSION (src/python/MethodLayout.cxx) ------------------------

// MethodLayout.cxx's kSideSamples, and its reason: alignment_quality() reads
// the sides between their points, and the 16 chords the viewer draws a long
// curved side with are each off the spline by half their turning.
constexpr int kSideSamples = 64;

// What MERIDIAN::run() and TORSION::run() hand Stage 10, from their options.
template <class PO>
::QuadMesh::Options stage10Options(const PO &o) {
    ::QuadMesh::Options q;
    q.targetEdgeLength = o.quadTargetEdge;
    q.minIntervals = o.quadMinIntervals;
    q.maxIntervals = o.quadMaxIntervals;
    q.collapseSpan = o.quadCollapseSpan;
    q.useSplines = o.quadUseSplines;
    q.featuresOnTracedArcs = o.quadFeaturesOnTracedArcs;
    q.materialsFromFaces = o.quadMaterialsFromFaces;
    q.contractOntoFeatures = o.quadContractOntoFeatures;
    q.smoothingPasses = o.quadSmoothingPasses;
    q.smoothingThreshold = o.quadSmoothingThreshold;
    return q;
}

// ... and Stage 11.
template <class PO>
DiskTemplate::Options stage11Options(const PO &o) {
    DiskTemplate::Options t;
    t.coreSquareness = o.diskCoreSquareness;
    t.ringDepth = o.diskRingDepth;
    t.smoothingPasses = o.diskSmoothingPasses;
    return t;
}

template <class P>
class LayoutCandidate final : public Candidate {
public:
    using Options = typename P::Options;

    LayoutCandidate(Method m, const char *source, const Mesh &model)
        : Candidate(m), source_(source), opts_(pipelineOptions()),
          pipeline_(std::make_shared<Mesh>(model), opts_) {}

    // Stages 10 and 11 the way the end of run() builds them, at `s` instead of
    // Options::quadTargetEdge: the excised rims as parity constraints on the
    // interval assignment, the mesh on the spline fit, and the O-grids.
    std::unique_ptr<CandidateMesh> mesh(const MeshSettings &s) const override {
        if (!pipeline_.hasSplines())
            throw std::runtime_error(std::string(name()) +
                                     ": no spline fit to mesh -- the pipeline stopped before Stage 9");
        ::QuadMesh::Options qo = stage10Options(opts_);
        qo.targetEdgeLength = s.target;
        qo.minIntervals = s.minIntervals;
        qo.maxIntervals = s.maxIntervals;
        const Arrangement &arr = pipeline_.getArrangement();
        const std::vector<DiskTemplate::Inclusion> &inclusions = pipeline_.getInclusions();
        const std::vector<std::vector<int>> rims = DiskTemplate::rimArcs(arr, inclusions);
        for (const std::vector<int> &r : rims)
            if (!r.empty()) qo.evenLoops.push_back(r);
        auto out = std::make_unique<CandidateMesh>();
        out->stage10 = std::make_unique<::QuadMesh>(pipeline_.getSplines(), qo);
        if (out->stage10->getReport().quads <= 0)
            throw std::runtime_error(std::string(name()) + ": the mesh has no elements");
        if (!inclusions.empty()) {
            const ::QuadMesh &q = *out->stage10;
            out->stage11 = std::make_unique<DiskTemplate>(
                q.vertices(), q.quads(), q.quadMaterials(), DiskTemplate::rimVertexLoops(arr, q, rims),
                inclusions, stage11Options(opts_));
        }
        return out;
    }

    // Stage 12 as TestMERIDIAN and TestTORSION run it, after the pillow.
    mesh::TMOP::Options smoothing() const override {
        mesh::TMOP::Options t;
        t.metric = mesh::TMOP::ShapeSize007;
        t.quadrature = mesh::TMOP::Corners;
        t.maxSweeps = 1000;
        return t;
    }
    bool pillows() const override { return true; }

private:
    static Options pipelineOptions() {
        Options o;
        o.runQuadMesh = false;   // Stage 10 is mesh()'s, at the target it is given
        return o;
    }

    void execute() override {
        pipeline_.run();
        if (!pipeline_.hasArrangement()) return;
        decomposition_ = pipeline_.getArrangement().blockDecomposition(
            pipeline_.hasSplines() ? &pipeline_.getSplines() : nullptr, kSideSamples, source_);
        // The area of the faces that are blocks -- four corners, one arc a
        // side, the ones Stage 10 meshes -- over the area of S after Stage 0c.
        const Arrangement &arr = pipeline_.getArrangement();
        const double domain = arr.getReport().domainArea;
        if (!(domain > 0.0)) return;
        const std::vector<Arrangement::Face> &faces = arr.getFaces();
        double blocks = 0.0;
        for (const int f : arr.patchFaces()) {
            if (f < 0 || f >= static_cast<int>(faces.size())) continue;
            if (faces[f].quad && faces[f].simple) blocks += std::fabs(faces[f].area);
        }
        coverage_ = blocks / domain;
    }

    const char *source_;
    Options opts_;
    P pipeline_;
};

// ---- ATLAS (src/python/MethodATLAS.cxx) ---------------------------------------
class AtlasCandidate final : public Candidate {
public:
    explicit AtlasCandidate(const Mesh &model)
        : Candidate(Method::ATLAS), mesh_(std::make_shared<Mesh>(model)), atlas_(mesh_, ATLAS::Options()) {}

    std::unique_ptr<CandidateMesh> mesh(const MeshSettings &s) const override {
        if (!atlas_.hasCover())
            throw std::runtime_error("atlas: no cover to mesh -- no search succeeded");
        BlockMesh::Options mo;
        mo.targetEdgeLength = s.target;
        mo.minIntervals = s.minIntervals;
        mo.maxIntervals = s.maxIntervals;
        auto out = std::make_unique<CandidateMesh>();
        out->blockMesh = std::make_unique<BlockMesh>(atlas_.getCover(), mo);
        if (out->blockMesh->getReport().quads <= 0)
            throw std::runtime_error("atlas: the mesh has no elements");
        return out;
    }

    mesh::TMOP::Options smoothing() const override { return blockGridSmoothing(); }

private:
    void execute() override {
        atlas_.run();
        if (!atlas_.hasCover()) return;
        decomposition_ = atlas_.getCover().blockDecomposition();
        coverage_ = decompositionCoverage(decomposition_, *mesh_);
    }

    std::shared_ptr<Mesh> mesh_;
    ATLAS atlas_;
};

// A method whose pipeline could not even be built: there is nothing in it
// but the reason.
class FailedCandidate final : public Candidate {
public:
    FailedCandidate(Method m, const std::string &why) : Candidate(m) {
        raised_ = true;
        error_ = why;
    }
    std::unique_ptr<CandidateMesh> mesh(const MeshSettings &) const override {
        throw std::runtime_error(std::string(name()) + ": the run failed (" + error_ + ")");
    }
    mesh::TMOP::Options smoothing() const override { return blockGridSmoothing(); }

private:
    void execute() override {}
};

}  // namespace

std::unique_ptr<Candidate> Candidate::run(Method method, const Mesh &model) {
    using Clock = std::chrono::steady_clock;
    const Clock::time_point t0 = Clock::now();
    std::unique_ptr<Candidate> c;
    try {
        switch (method) {
            case Method::ZIPLINE:  c = std::make_unique<ZiplineCandidate>(model); break;
            case Method::UMBER:    c = std::make_unique<UmberCandidate>(model); break;
            case Method::MERIDIAN:
                c = std::make_unique<LayoutCandidate<MERIDIAN>>(method, "MERIDIAN", model);
                break;
            case Method::TORSION:
                c = std::make_unique<LayoutCandidate<TORSION>>(method, "TORSION", model);
                break;
            case Method::ATLAS:    c = std::make_unique<AtlasCandidate>(model); break;
        }
    } catch (const std::exception &e) {
        c = std::make_unique<FailedCandidate>(method, e.what());
    }
    if (!c) c = std::make_unique<FailedCandidate>(method, "unknown method");
    if (!c->raised_) {
        try {
            c->execute();
        } catch (const std::exception &e) {
            c->raised_ = true;
            c->error_ = e.what();
        } catch (...) {
            c->raised_ = true;
            c->error_ = "an exception that is not a std::exception";
        }
        // A run that raised has no result, whatever it had built before it
        // did: coverage 0 and no block count, as the dataset records it.
        if (c->raised_) {
            c->decomposition_ = BlockDecomposition();
            c->coverage_ = 0.0;
        }
    }
    c->seconds_ = std::chrono::duration<double>(Clock::now() - t0).count();
    return c;
}

}  // namespace oracle
