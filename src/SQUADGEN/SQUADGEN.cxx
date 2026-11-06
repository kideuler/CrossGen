#include "SQUADGEN/SQUADGEN.hxx"

#include <chrono>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace {

using Clock = std::chrono::steady_clock;

double since(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

// The fraction of the model's area inside a block, as crossgen's BlockDecomposition.coverage takes it
// (src/python/Methods.cxx): the shoelace over every block's four sides.
double coverageOf(const BlockDecomposition &D, const Mesh &mesh) {
    double modelArea = 0.0;
    for (const Triangle &t : mesh.triangles) {
        const Point &a = mesh.vertices[t[0]], &b = mesh.vertices[t[1]], &c = mesh.vertices[t[2]];
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

}  // namespace

SQUADGEN::SQUADGEN(std::shared_ptr<const Networks> networks, Model model, const Options &opts)
    : networks_(std::move(networks)), model_(std::move(model)), opts_(opts), domain_(model_) {
    if (!networks_) throw std::runtime_error("SQUADGEN: no networks");
}

void SQUADGEN::clearFrom(int stage) {
    if (stage <= 1) {
        tokens_.clear();
        context_.clear();
        tokenSDF_.clear();
        contextSDF_.clear();
        code_.clear();
        pe_.clear();
        report_ = Report();
    }
    if (stage <= 2) z0_.clear();
    if (stage <= 3) {
        latent_.clear();
        report_.sampled = false;
        report_.evaluations = report_.steps = report_.rejected = 0;
        report_.sampleSeconds = 0.0;
    }
    if (stage <= 4) {
        extraction_.reset();
        extractionModel_ = Model();
        fields_ = Networks::Fields();
        report_.levels = report_.triangles = 0;
        report_.decodeSeconds = 0.0;
    }
    if (stage <= 5) {
        quarry_.reset();
        report_.extracted = report_.valid = false;
        report_.coverage = 0.0;
        report_.blocks = 0;
        report_.extractSeconds = 0.0;
    }
}

// ── Stage 1 ──────────────────────────────────────────────────────────────────

void SQUADGEN::condition() {
    const int n = networks_->description().tokens;
    const Clock::time_point t0 = Clock::now();
    // One set of each size from one generator, as domain_condition draws its pair from one numpy generator.
    squadgen::Random rng(opts_.seed);
    std::vector<Point> tokens = domain_.tokenSet(n, rng);
    std::vector<Point> context = domain_.tokenSet(4 * n, rng);
    const double fps = since(t0);
    condition(std::move(tokens), std::move(context));
    report_.conditionSeconds += fps;
}

void SQUADGEN::condition(std::vector<Point> tokens, std::vector<Point> context) {
    const Networks::Description &d = networks_->description();
    if (static_cast<int>(tokens.size()) != d.tokens || static_cast<int>(context.size()) != 4 * d.tokens)
        throw std::runtime_error("SQUADGEN: the networks read " + std::to_string(d.tokens) + " tokens and " +
                                 std::to_string(4 * d.tokens) + " key/value points, given " +
                                 std::to_string(tokens.size()) + " and " + std::to_string(context.size()));
    clearFrom(1);
    const Clock::time_point t0 = Clock::now();
    tokens_ = std::move(tokens);
    context_ = std::move(context);
    tokenSDF_ = domain_.insideSDF(tokens_, d.tau);
    contextSDF_ = domain_.insideSDF(context_, d.tau);
    const std::vector<float> xyz = squadgen::flatten(tokens_);
    code_ = networks_->geometryCode(xyz, tokenSDF_, squadgen::flatten(context_), contextSDF_);
    pe_ = networks_->conditioning(code_, xyz, tokenSDF_);
    report_.tokens = static_cast<int>(tokens_.size());
    report_.contextPoints = static_cast<int>(context_.size());
    report_.segments = static_cast<int>(domain_.segments().size());
    report_.conditionSeconds = since(t0);
}

// ── Stage 2 ──────────────────────────────────────────────────────────────────

void SQUADGEN::noise(unsigned seed) {
    if (!hasCondition()) throw std::runtime_error("SQUADGEN: the noise needs Stage 1's tokens");
    const size_t N = tokens_.size(), C = static_cast<size_t>(networks_->description().latentDim);
    // Offset from the seed Stage 1 drew its tokens with, so that one seed does not give both the same bits.
    squadgen::Random rng(static_cast<uint64_t>(seed) + 0x9E3779B97F4A7C15ULL);
    std::vector<float> z(N * C);
    for (float &v : z) v = static_cast<float>(rng.normal());
    const int blocked = opts_.taubin ? squadgen::taubinSmooth(z, static_cast<int>(C), tokens_, domain_.segments()) : 0;
    setNoise(std::move(z));
    report_.noiseSeed = seed;
    report_.blockedPairs = blocked;
}

void SQUADGEN::setNoise(std::vector<float> z0) {
    if (!hasCondition()) throw std::runtime_error("SQUADGEN: the noise needs Stage 1's tokens");
    if (z0.size() != tokens_.size() * static_cast<size_t>(networks_->description().latentDim))
        throw std::runtime_error("SQUADGEN: the initial latent is not N x C");
    clearFrom(2);
    z0_ = std::move(z0);
    double s = 0.0;
    for (const float v : z0_) s += static_cast<double>(v) * v;
    report_.noiseStd = std::sqrt(s / static_cast<double>(z0_.size()));
}

// ── Stage 3 ──────────────────────────────────────────────────────────────────

bool SQUADGEN::sample(const Sampler::Observer &observe, const std::atomic<bool> *stop) {
    if (!hasNoise()) throw std::runtime_error("SQUADGEN: the ODE needs Stage 2's noise");
    clearFrom(3);
    const Clock::time_point t0 = Clock::now();
    struct Stopped {};
    auto f = [this, stop](const std::vector<float> &x, double t, std::vector<float> &v) {
        if (stop && stop->load()) throw Stopped();
        networks_->velocity(x, t, pe_, v);
    };
    Sampler::Result r;
    try {
        r = Sampler::integrate(f, z0_, opts_.sampler, observe);
    } catch (const Stopped &) {
        report_.sampleSeconds = since(t0);
        return false;
    }
    report_.evaluations = r.evaluations;
    report_.steps = r.steps;
    report_.rejected = r.rejected;
    report_.sampleSeconds = since(t0);
    if (!r.finished) return false;
    latent_ = std::move(r.x);
    report_.sampled = true;
    return true;
}

void SQUADGEN::setLatent(std::vector<float> x1) {
    if (!hasCondition()) throw std::runtime_error("SQUADGEN: a latent needs Stage 1's tokens");
    if (x1.size() != tokens_.size() * static_cast<size_t>(networks_->description().latentDim))
        throw std::runtime_error("SQUADGEN: the latent is not N x C");
    clearFrom(3);
    latent_ = std::move(x1);
    report_.sampled = true;
}

// ── Stage 4 ──────────────────────────────────────────────────────────────────

squadgen::Networks::Fields SQUADGEN::decodeAt(const std::vector<float> &normalised,
                                              const std::vector<float> &framePoints) const {
    return networks_->query(networks_->decoderTokens(networks_->denormalize(normalised)), framePoints);
}

void SQUADGEN::decode() {
    if (!hasLatent()) throw std::runtime_error("SQUADGEN: decoding needs Stage 3's latent");
    clearFrom(4);
    const Clock::time_point t0 = Clock::now();
    int levels = 0;
    extractionModel_ = squadgen::extractionMesh(model_, opts_.hRel, &levels);
    extraction_ = std::make_shared<const Mesh>(extractionModel_.vertices, extractionModel_.triangles);
    fields_ = decodeAt(latent_, squadgen::centroidsInFrame(extractionModel_, frame()));
    report_.levels = levels;
    report_.triangles = static_cast<int>(extractionModel_.triangles.size());
    report_.decodeSeconds = since(t0);
}

// ── Stage 5 ──────────────────────────────────────────────────────────────────

void SQUADGEN::extract() {
    if (!hasFields()) throw std::runtime_error("SQUADGEN: QUARRY needs Stage 4's fields");
    clearFrom(5);
    const Clock::time_point t0 = Clock::now();
    // In the model's units, as write_quarry_input hands them over: the offsets divided by the frame's scale, in
    // float32 as numpy divides a float32 array, then read as doubles.
    const size_t F = fields_.size();
    const float s = static_cast<float>(frame().scale);
    QUARRY::Fields f;
    f.cdf.assign(fields_.cdf.begin(), fields_.cdf.end());
    f.dcdf.assign(fields_.dcdf.begin(), fields_.dcdf.end());
    f.offsetCentre.resize(F);
    f.offsetCorner.resize(F);
    for (size_t t = 0; t < F; ++t) {
        f.offsetCentre[t] = {fields_.offsetc[2 * t] / s, fields_.offsetc[2 * t + 1] / s};
        f.offsetCorner[t] = {fields_.offset[2 * t] / s, fields_.offset[2 * t + 1] / s};
    }
    quarry_ = std::make_unique<QUARRY>(extraction_, std::move(f), opts_.quarry);
    quarry_->run();
    const BlockDecomposition &D = quarry_->decomposition();
    report_.extracted = true;
    report_.coverage = coverageOf(D, *extraction_);
    report_.valid = D.covers() && report_.coverage >= 0.999;
    report_.blocks = static_cast<int>(D.blocks.size());
    report_.extractSeconds = since(t0);
}

void SQUADGEN::run() {
    condition();
    noise(opts_.seed);
    sample();
    decode();
    extract();
}
