#ifndef __SQUADGEN_NETWORKS_HXX__
#define __SQUADGEN_NETWORKS_HXX__

#include <memory>
#include <string>
#include <vector>

namespace squadgen {

// The three networks of SQuadGen's inference (docs/squadgen_2d.md Sec. 4) as the five ONNX graphs
// `python -m squadgen2d export-onnx` writes (py/squadgen/squadgen2d/export_onnx.py, whose docstring is the
// specification of every input and output): the Geom-AE's code at the tokens, the SiT's per-block conditioning and its
// velocity, and the SQ-VAE decoder's trunk and per-point query. Batch 1 throughout; every tensor is float32 and
// row-major, and every point and field is in the domain's normalised frame (squadgen::Frame).
//
//   geometryCode   (N, 2), (N, 1), (4N, 2), (4N, 1) token sets and their truncated SDF  ->  (N, codeDim)
//   conditioning   the code and the N-set again                                        ->  (depth, N, hidden)
//   velocity       a normalised latent (N, C), a time t in [0, 1] and the conditioning ->  (N, C)
//   decoderTokens  a denormalised latent (N, C)                                        ->  (N, vaeDim)
//   query          the decoder tokens and P points (P, 2)                               ->  cdf, dcdf (P), offset,
//                                                                                         offsetc (P, 2)
//
// ### What the files say about themselves
//
// Each graph carries squadgen.* metadata (export_onnx.metadata): N, the latent width C, the SiT's depth and width,
// the code and decoder widths, the truncated SDF's tau, the latent statistics target_mu / target_std the SiT's
// latents are normalised by, and the checkpoint and step it came from. The five carry one export id, and a directory
// whose files disagree on it -- two exports mixed -- is refused, as is one whose graphs' metadata disagrees with
// their own tensor shapes.
//
// ### Why ONNX Runtime is behind this class
//
// For the reason oracle::Selector gives: it is a large dependency only this file needs, so it is linked PRIVATE to
// the SQUADGEN library and its headers reach Networks.cxx alone. A build without it still builds SQUADGEN, and the
// constructor says why it cannot read the networks instead.
//
// Every graph runs on ONNX Runtime's CPU provider. Measured on the 8 GB M2 at N = 2048 (2026-11-02): ~0.55 s a
// velocity, 0.4 s for the decoder tokens, ~10 us a query point, against 0.43 s / 0.31 s / ~7 us in torch on the same
// CPU. Its CoreML provider (the GPU) took 40-50 s a velocity, the graph split into 25 or more partitions, so it is
// not offered. The thread pool does not spin between operators: the viewer samples on a worker thread and its GUI
// thread needs the cores the spinning would hold.
class Networks {
public:
    struct Description {
        std::string directory;
        std::string exportId;       // the same in all five files
        std::string model;          // the SQ-Diffuse checkpoint exported, repository-relative
        int step = 0;               // its training step
        int tokens = 0;             // N
        int latentDim = 0;          // C
        int depth = 0;              // SiT blocks
        int hidden = 0;             // SiT width
        int codeDim = 0;            // Geom-AE code width
        int vaeDim = 0;             // SQ-VAE decoder width
        double tau = 0.1;           // the truncated SDF's, in frame units
        std::vector<double> targetMu, targetStd;   // the latent normalisation, one per channel
    };

    struct Options {
        int threads = 0;            // intra-op threads per graph; 0 is ONNX Runtime's own choice
    };

    // The four fields the decoder's query heads give QUARRY, per point: cdf and dcdf (P), offset and offsetc
    // (P x 2, row-major), all in frame units.
    struct Fields {
        std::vector<float> cdf, dcdf, offset, offsetc;
        size_t size() const { return cdf.size(); }
    };

    // Loads the five graphs in `directory`. Throws std::runtime_error with the reason when it cannot: a graph
    // missing, files from two exports, metadata that disagrees with the graphs, or a build without ONNX Runtime.
    explicit Networks(const std::string &directory, const Options &opts);
    explicit Networks(const std::string &directory) : Networks(directory, Options()) {}
    ~Networks();
    Networks(const Networks &) = delete;
    Networks &operator=(const Networks &) = delete;

    const Description &description() const;

    std::vector<float> geometryCode(const std::vector<float> &xyz, const std::vector<float> &sdf,
                                    const std::vector<float> &xyzKV, const std::vector<float> &sdfKV) const;
    std::vector<float> conditioning(const std::vector<float> &code, const std::vector<float> &xyz,
                                    const std::vector<float> &sdf) const;
    // One function evaluation; v is resized to N x C.
    void velocity(const std::vector<float> &x, double t, const std::vector<float> &pe, std::vector<float> &v) const;
    std::vector<float> decoderTokens(const std::vector<float> &latent) const;
    // The fields at `xy` (P x 2), in chunks of `chunk` points so that the attention's P x N scores stay small.
    Fields query(const std::vector<float> &tokens, const std::vector<float> &xy, size_t chunk = 8192) const;

    // The latent the decoder reads from a normalised one, and back: z = x * target_std + target_mu.
    std::vector<float> denormalize(const std::vector<float> &x) const;

    // Whether this build can read the networks at all.
    static bool available();

    // The directory a caller gets when it names none: $CROSSGEN_SQUADGEN when it is set, best_models/onnx of the
    // source tree this was built from otherwise.
    static std::string defaultDirectory();

    // The five file names, in the order the description above lists the graphs.
    static const std::vector<std::string> &graphNames();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace squadgen

#endif // __SQUADGEN_NETWORKS_HXX__
