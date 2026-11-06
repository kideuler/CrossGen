#ifndef __SQUADGEN_HXX__
#define __SQUADGEN_HXX__

#include <atomic>
#include <memory>
#include <string>
#include <vector>

#include "QUARRY/QUARRY.hxx"
#include "SQUADGEN/Domain.hxx"
#include "SQUADGEN/Networks.hxx"
#include "SQUADGEN/Sampler.hxx"
#include "mesh/BlockDecomposition.hxx"
#include "mesh/Mesh.hxx"

// SQUADGEN: SQuadGen's inference (docs/squadgen_2d.md Sec. 4) in C++, from a domain to the quad blocks a trained
// SQ-Diffuse generates for it -- what `python -m squadgen2d generate` does, without Python. The networks are the ONNX
// graphs `python -m squadgen2d export-onnx` writes (squadgen::Networks); everything around them is ported from the
// squadgen2d module that runs it there, and named after it below.
//
//   Stage 1  condition()  Sec. 4 steps 1-3: the domain in the networks' frame (squadgen::Frame), N tokens and 4N
//                         key/value points by farthest-point sampling of its area, the truncated SDF at each, the
//                         Geom-AE's code at the tokens, and the SiT's per-block conditioning from it, which depends
//                         on the geometry alone and so is made once (sq_diffuse.domain_condition, SiT.cond_embed)
//   Stage 2  noise()      the initial latent: white noise per token and channel, Taubin-smoothed over the tokens
//                         with pairs across dS blocked (Appendix C, sq_diffuse.taubin_smooth)
//   Stage 3  sample()     the ODE from t = 0 to 1, dopri5 at atol 1e-6 and rtol 1e-3 by default (the release's
//                         inference, sq_diffuse.sample_latents), each evaluation one run of the velocity graph; an
//                         observer is shown every accepted step
//   Stage 4  decode()     the extraction mesh (the model split until its median edge is about h_rel x its diagonal,
//                         roundtrip.extraction_mesh) and the latent decoded at every triangle's centroid: cdf, dcdf,
//                         and the offsets to the chart centre and the subchart corner (sq_vae.decode_fields)
//   Stage 5  extract()    QUARRY on those fields, the offsets taken back to the model's units (quarry.py): the
//                         BlockDecomposition, which BlockQuadMesh meshes and TMOP smooths as every method's blocks are
//
// Each stage needs the one before it and clears the ones after it, so a second seed is noise() and on, with the
// condition kept, and a stage can be watched on its own: the viewer's mode 0 runs them one phase at a time, Stage 3
// on a worker thread with the observer drawing the latent as it goes.
//
// ### What differs from squadgen2d
//
// The networks are the same to float32 rounding (export_onnx checks them against torch, TestSQUADGEN --check against
// its reference case). The random draws are not: the token sets and the noise come from squadgen::Random, not numpy
// and torch, so a seed here is not the seed of the same number there, and gives a different sample of the same
// distribution. The reference case carries the draws for a comparison stage by stage. The model is read as
// squadgen2d reads it (squadgen::readOBJ), so the frame and the extraction mesh are the same.
//
// What a sample comes to is not judged here. generate keeps the best of several seeds -- blocks that cover the model,
// a mesh with no inverted element, then the highest alignment quality -- and TestSQUADGEN --samples does the same.
class SQUADGEN {
public:
    using Networks = squadgen::Networks;
    using Model = squadgen::Model;
    using Sampler = squadgen::Sampler;
    using Domain = squadgen::Domain;
    using Frame = squadgen::Frame;

    struct Options {
        unsigned seed = 0;              // the token sets of Stage 1
        bool taubin = true;             // Taubin-smoothed noise (the release's inference), or white (as trained)
        Sampler::Options sampler;       // Stage 3
        double hRel = 0.002;            // Stage 4's extraction mesh, as a fraction of the diagonal (the release's)
        QUARRY::Options quarry;         // Stage 5
    };

    struct Report {
        // Stage 1
        int tokens = 0, contextPoints = 0, segments = 0;
        double conditionSeconds = 0.0;
        // Stage 2
        unsigned noiseSeed = 0;
        int blockedPairs = 0;
        double noiseStd = 0.0;          // over every token and channel, after smoothing
        // Stage 3
        bool sampled = false;
        int evaluations = 0, steps = 0, rejected = 0;
        double sampleSeconds = 0.0;
        // Stage 4
        int levels = 0;                 // subdivisions of the model
        int triangles = 0;              // of the extraction mesh
        double decodeSeconds = 0.0;
        // Stage 5
        bool extracted = false;
        bool valid = false;             // covers() and at least 99.9 % of the area in a block: crossgen's `valid`
        double coverage = 0.0;
        int blocks = 0;
        double extractSeconds = 0.0;
    };

    // The model as squadgen2d reads it (squadgen::readOBJ), or a Mesh's arrays (squadgen::modelOf), and the networks
    // to run on it, which a caller may share between models: loading them takes a few seconds.
    SQUADGEN(std::shared_ptr<const Networks> networks, Model model, const Options &opts);

    void condition();
    // Stage 1 on token sets given in the frame, rather than drawn: the reference case's, or another run's.
    void condition(std::vector<Point> tokens, std::vector<Point> context);
    void noise(unsigned seed);
    // Stage 2 with the initial latent given (normalised, N x C, smoothed or not as the caller wants it).
    void setNoise(std::vector<float> z0);
    // Returns false when the observer stopped it, or `stop` was set -- it is read before every evaluation, so a
    // caller on another thread waits at most one for the solve to end -- and then nothing after it can run.
    bool sample(const Sampler::Observer &observe = {}, const std::atomic<bool> *stop = nullptr);
    // Stage 3 with the latent given (normalised, N x C).
    void setLatent(std::vector<float> x1);
    void decode();
    void extract();
    void run();   // Stages 1 to 5 at the options' seeds

    const Options &options() const { return opts_; }
    Options &options() { return opts_; }
    const Report &getReport() const { return report_; }
    const Networks &networks() const { return *networks_; }
    const Model &model() const { return model_; }
    const Domain &domain() const { return domain_; }
    const Frame &frame() const { return domain_.frame(); }

    bool hasCondition() const { return !pe_.empty(); }
    bool hasNoise() const { return !z0_.empty(); }
    bool hasLatent() const { return !latent_.empty(); }
    bool hasFields() const { return extraction_ != nullptr && fields_.size() > 0; }
    bool hasBlocks() const { return quarry_ != nullptr; }

    // Stage 1, in the frame.
    const std::vector<Point> &tokens() const { return tokens_; }
    const std::vector<Point> &contextPoints() const { return context_; }
    const std::vector<float> &tokenSDF() const { return tokenSDF_; }
    const std::vector<float> &contextSDF() const { return contextSDF_; }
    const std::vector<float> &code() const { return code_; }
    const std::vector<float> &conditioning() const { return pe_; }
    // Stage 2 and 3, normalised latents N x C.
    const std::vector<float> &initialNoise() const { return z0_; }
    const std::vector<float> &latent() const { return latent_; }
    // Stage 4. The extraction mesh in the model's units; the fields per triangle, in the frame's.
    const std::shared_ptr<const Mesh> &extractionMesh() const { return extraction_; }
    const Model &extractionModel() const { return extractionModel_; }
    const Networks::Fields &fields() const { return fields_; }
    // Stage 5.
    const QUARRY &quarry() const { return *quarry_; }
    const BlockDecomposition &decomposition() const { return quarry_->decomposition(); }

    // The fields a normalised latent decodes to at points of the frame (float32 pairs), in frame units: Stage 4 on
    // the extraction mesh's centroids, and what a picture of an unfinished solve shows on a coarser one.
    Networks::Fields decodeAt(const std::vector<float> &normalised, const std::vector<float> &framePoints) const;

private:
    void clearFrom(int stage);

    std::shared_ptr<const Networks> networks_;
    Model model_;
    Options opts_;
    Domain domain_;
    Report report_;
    std::vector<Point> tokens_, context_;
    std::vector<float> tokenSDF_, contextSDF_, code_, pe_;
    std::vector<float> z0_, latent_;
    Model extractionModel_;
    std::shared_ptr<const Mesh> extraction_;
    Networks::Fields fields_;
    std::unique_ptr<QUARRY> quarry_;
};

#endif // __SQUADGEN_HXX__
