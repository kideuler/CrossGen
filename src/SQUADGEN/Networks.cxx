#include "SQUADGEN/Networks.hxx"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

#ifdef CROSSGEN_HAVE_ONNXRUNTIME
#include <onnxruntime_cxx_api.h>
#endif

// CMake sets it to the source tree's best_models/onnx; this is only for a build that does not.
#ifndef CROSSGEN_SQUADGEN_DEFAULT
#define CROSSGEN_SQUADGEN_DEFAULT "best_models/onnx"
#endif

namespace squadgen {
namespace {

enum GraphId { GeomCode = 0, SiTCondition, SiTVelocity, VaeTokens, VaeQuery, GraphCount };

const std::vector<std::string> kGraphs = {"geom_code", "sit_condition", "sit_velocity", "vae_tokens", "vae_query"};

// Each graph's inputs and outputs in the order this file passes them, by the names export_onnx gives them.
const std::vector<std::vector<std::string>> kInputs = {
    {"xyz", "sdf", "xyz_kv", "sdf_kv"}, {"code", "xyz", "sdf"}, {"x", "t", "pe"}, {"z"}, {"tokens", "xy"}};
const std::vector<std::vector<std::string>> kOutputs = {
    {"code"}, {"pe"}, {"v"}, {"tokens"}, {"cdf", "dcdf", "offset", "offsetc"}};

std::string fileOf(const std::string &dir, int g) {
    return dir + (dir.empty() || dir.back() == '/' ? "" : "/") + kGraphs[g] + ".onnx";
}

#ifdef CROSSGEN_HAVE_ONNXRUNTIME
// One per process, as ONNX Runtime asks, made on first use. ORACLE's selector makes its own the same way; the
// runtime keeps a single environment underneath and counts the references to it.
Ort::Env &environment() {
    static Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "SQUADGEN");
    return env;
}

// A loaded graph, its metadata, and the name arrays Run() takes in this file's order.
struct Graph {
    Ort::Session session{nullptr};
    std::map<std::string, std::string> meta;
    std::vector<const char *> in, out;
    std::vector<std::vector<int64_t>> inShapes, outShapes;   // -1 where an axis is dynamic
};

// A comma-separated list of numbers, as export_onnx writes target_mu / target_std.
std::vector<double> numbers(const std::string &s) {
    std::vector<double> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) out.push_back(std::stod(item));
    return out;
}

int integer(const std::map<std::string, std::string> &meta, const std::string &key, const std::string &file) {
    const auto it = meta.find(key);
    if (it == meta.end()) throw std::runtime_error(file + " has no '" + key + "' in its metadata: not an export_onnx graph");
    return std::stoi(it->second);
}

Ort::Value tensor(const Ort::MemoryInfo &mem, const std::vector<float> &data, std::vector<int64_t> shape) {
    // Run() only reads its inputs; the C API takes them as mutable all the same.
    return Ort::Value::CreateTensor<float>(mem, const_cast<float *>(data.data()), data.size(), shape.data(),
                                           shape.size());
}

Ort::Value tensor(const Ort::MemoryInfo &mem, const float *data, size_t n, std::vector<int64_t> shape) {
    return Ort::Value::CreateTensor<float>(mem, const_cast<float *>(data), n, shape.data(), shape.size());
}

void copyOut(const Ort::Value &v, std::vector<float> &into) {
    const size_t n = v.GetTensorTypeAndShapeInfo().GetElementCount();
    const float *p = v.GetTensorData<float>();
    into.assign(p, p + n);
}

void appendOut(const Ort::Value &v, std::vector<float> &into) {
    const size_t n = v.GetTensorTypeAndShapeInfo().GetElementCount();
    const float *p = v.GetTensorData<float>();
    into.insert(into.end(), p, p + n);
}
#endif

}  // namespace

struct Networks::Impl {
    Description d;
#ifdef CROSSGEN_HAVE_ONNXRUNTIME
    // Run() is not const in the C++ API, though a session may run on several threads at once.
    mutable Graph graphs[GraphCount];
    Ort::MemoryInfo mem = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);

    std::vector<Ort::Value> run(int g, std::vector<Ort::Value> &inputs) const {
        Graph &G = graphs[g];
        try {
            return G.session.Run(Ort::RunOptions{nullptr}, G.in.data(), inputs.data(), inputs.size(), G.out.data(),
                                 G.out.size());
        } catch (const Ort::Exception &e) {
            throw std::runtime_error("SQuadGen's " + kGraphs[g] + " failed: " + e.what());
        }
    }
#endif
};

Networks::Networks(const std::string &directory, const Options &opts) : impl_(std::make_unique<Impl>()) {
    Description &d = impl_->d;
    d.directory = directory;
    for (int g = 0; g < GraphCount; ++g)
        if (!std::ifstream(fileOf(directory, g)).good())
            throw std::runtime_error("no " + kGraphs[g] + ".onnx in '" + directory +
                                     "' (python -m squadgen2d export-onnx writes the five graphs; $CROSSGEN_SQUADGEN "
                                     "names another directory)");
#ifndef CROSSGEN_HAVE_ONNXRUNTIME
    (void)opts;
    throw std::runtime_error("this build has no ONNX Runtime to run SQuadGen's networks with: install it (brew install "
                             "onnxruntime) and run cmake again");
#else
    Ort::SessionOptions so;
    so.SetIntraOpNumThreads(opts.threads);
    so.SetInterOpNumThreads(1);
    so.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
    so.AddConfigEntry("session.intra_op.allow_spinning", "0");
    Ort::AllocatorWithDefaultOptions alloc;
    for (int g = 0; g < GraphCount; ++g) {
        const std::string file = fileOf(directory, g);
        Graph &G = impl_->graphs[g];
        try {
            G.session = Ort::Session(environment(), file.c_str(), so);
            const Ort::ModelMetadata md = G.session.GetModelMetadata();
            for (const Ort::AllocatedStringPtr &k : md.GetCustomMetadataMapKeysAllocated(alloc)) {
                const Ort::AllocatedStringPtr v = md.LookupCustomMetadataMapAllocated(k.get(), alloc);
                if (v) G.meta[k.get()] = v.get();
            }
            // The names in this file's order, each checked to be one the graph has.
            auto collect = [&](size_t count, bool input, const std::vector<std::string> &want,
                               std::vector<const char *> &names, std::vector<std::vector<int64_t>> &shapes) {
                std::map<std::string, std::vector<int64_t>> have;
                for (size_t i = 0; i < count; ++i) {
                    const Ort::AllocatedStringPtr name =
                        input ? G.session.GetInputNameAllocated(i, alloc) : G.session.GetOutputNameAllocated(i, alloc);
                    const Ort::TypeInfo info = input ? G.session.GetInputTypeInfo(i) : G.session.GetOutputTypeInfo(i);
                    have[name.get()] = info.GetTensorTypeAndShapeInfo().GetShape();
                }
                if (have.size() != want.size())
                    throw std::runtime_error("it has " + std::to_string(have.size()) + (input ? " inputs" : " outputs") +
                                             ", export_onnx's has " + std::to_string(want.size()));
                for (const std::string &w : want) {
                    const auto it = have.find(w);
                    if (it == have.end())
                        throw std::runtime_error(std::string("it has no ") + (input ? "input" : "output") + " '" + w +
                                                 "'");
                    names.push_back(w.c_str());   // kInputs / kOutputs live as long as the program
                    shapes.push_back(it->second);
                }
            };
            collect(G.session.GetInputCount(), true, kInputs[g], G.in, G.inShapes);
            collect(G.session.GetOutputCount(), false, kOutputs[g], G.out, G.outShapes);
        } catch (const Ort::Exception &e) {
            throw std::runtime_error("cannot read '" + file + "': " + e.what());
        } catch (const std::runtime_error &e) {
            throw std::runtime_error("'" + file + "': " + e.what());
        }
    }

    // The description from the velocity graph, which is the one every sample runs a hundred times; every other
    // graph has to come from the same export.
    const Graph &V = impl_->graphs[SiTVelocity];
    const std::string vfile = fileOf(directory, SiTVelocity);
    auto text = [&](const std::string &key) {
        const auto it = V.meta.find(key);
        if (it == V.meta.end())
            throw std::runtime_error("'" + vfile + "' has no '" + key + "' in its metadata: not an export_onnx graph");
        return it->second;
    };
    d.exportId = text("squadgen.export");
    d.model = text("squadgen.model");
    d.step = integer(V.meta, "squadgen.step", vfile);
    d.tokens = integer(V.meta, "squadgen.tokens", vfile);
    d.latentDim = integer(V.meta, "squadgen.latent_dim", vfile);
    d.depth = integer(V.meta, "squadgen.depth", vfile);
    d.hidden = integer(V.meta, "squadgen.hidden", vfile);
    d.codeDim = integer(V.meta, "squadgen.code_dim", vfile);
    d.vaeDim = integer(V.meta, "squadgen.vae_dim", vfile);
    d.tau = std::stod(text("squadgen.tau"));
    d.targetMu = numbers(text("squadgen.target_mu"));
    d.targetStd = numbers(text("squadgen.target_std"));
    for (int g = 0; g < GraphCount; ++g) {
        const auto it = impl_->graphs[g].meta.find("squadgen.export");
        if (it == impl_->graphs[g].meta.end() || it->second != d.exportId)
            throw std::runtime_error("'" + fileOf(directory, g) + "' is from another export than '" + vfile +
                                     "' (export id " + (it == impl_->graphs[g].meta.end() ? "none" : it->second) +
                                     " vs " + d.exportId + "): export all five again");
    }

    // What the graphs' own shapes can check of that.
    auto require = [&](bool ok, const std::string &what) {
        if (!ok) throw std::runtime_error("the SQuadGen graphs in '" + directory + "': " + what);
    };
    require(d.tokens > 0 && d.latentDim > 0 && d.depth > 0 && d.hidden > 0 && d.codeDim > 0 && d.vaeDim > 0,
            "a width in the metadata is not positive");
    require(static_cast<int>(d.targetMu.size()) == d.latentDim && static_cast<int>(d.targetStd.size()) == d.latentDim,
            "target_mu / target_std do not have one entry per latent channel");
    const std::vector<int64_t> &pe = impl_->graphs[SiTVelocity].inShapes[2];
    require(pe.size() == 4 && pe[0] == d.depth && pe[3] == d.hidden, "the velocity's pe input is not (depth, 1, N, hidden)");
    require(impl_->graphs[SiTVelocity].inShapes[0].back() == d.latentDim, "the velocity's x is not latent_dim wide");
    require(impl_->graphs[GeomCode].outShapes[0].back() == d.codeDim, "the code is not code_dim wide");
    require(impl_->graphs[VaeTokens].outShapes[0].back() == d.vaeDim, "the decoder tokens are not vae_dim wide");
#endif
}

Networks::~Networks() = default;

const Networks::Description &Networks::description() const { return impl_->d; }

std::vector<float> Networks::denormalize(const std::vector<float> &x) const {
    const Description &d = impl_->d;
    const size_t C = static_cast<size_t>(d.latentDim);
    std::vector<float> z(x.size());
    for (size_t i = 0; i < x.size(); ++i)
        z[i] = static_cast<float>(x[i] * d.targetStd[i % C] + d.targetMu[i % C]);
    return z;
}

#ifdef CROSSGEN_HAVE_ONNXRUNTIME
std::vector<float> Networks::geometryCode(const std::vector<float> &xyz, const std::vector<float> &sdf,
                                          const std::vector<float> &xyzKV, const std::vector<float> &sdfKV) const {
    const int64_t n = static_cast<int64_t>(sdf.size()), m = static_cast<int64_t>(sdfKV.size());
    if (xyz.size() != 2 * sdf.size() || xyzKV.size() != 2 * sdfKV.size())
        throw std::runtime_error("geometryCode: a token set and its sdf differ in length");
    const Ort::MemoryInfo &mem = impl_->mem;
    std::vector<Ort::Value> in;
    in.push_back(tensor(mem, xyz, {1, n, 2}));
    in.push_back(tensor(mem, sdf, {1, n, 1}));
    in.push_back(tensor(mem, xyzKV, {1, m, 2}));
    in.push_back(tensor(mem, sdfKV, {1, m, 1}));
    std::vector<float> out;
    copyOut(impl_->run(GeomCode, in)[0], out);
    return out;
}

std::vector<float> Networks::conditioning(const std::vector<float> &code, const std::vector<float> &xyz,
                                          const std::vector<float> &sdf) const {
    const int64_t n = static_cast<int64_t>(sdf.size());
    if (code.size() != static_cast<size_t>(n * impl_->d.codeDim) || xyz.size() != 2 * sdf.size())
        throw std::runtime_error("conditioning: the code, the tokens and their sdf differ in length");
    const Ort::MemoryInfo &mem = impl_->mem;
    std::vector<Ort::Value> in;
    in.push_back(tensor(mem, code, {1, n, impl_->d.codeDim}));
    in.push_back(tensor(mem, xyz, {1, n, 2}));
    in.push_back(tensor(mem, sdf, {1, n, 1}));
    std::vector<float> out;
    copyOut(impl_->run(SiTCondition, in)[0], out);
    return out;
}

void Networks::velocity(const std::vector<float> &x, double t, const std::vector<float> &pe,
                        std::vector<float> &v) const {
    const Description &d = impl_->d;
    const int64_t n = static_cast<int64_t>(x.size()) / d.latentDim;
    if (x.size() != static_cast<size_t>(n * d.latentDim) ||
        pe.size() != static_cast<size_t>(d.depth) * static_cast<size_t>(n) * static_cast<size_t>(d.hidden))
        throw std::runtime_error("velocity: the latent and the conditioning differ in their token count");
    const float tt = static_cast<float>(t);
    const Ort::MemoryInfo &mem = impl_->mem;
    std::vector<Ort::Value> in;
    in.push_back(tensor(mem, x, {1, n, d.latentDim}));
    in.push_back(tensor(mem, &tt, 1, {1}));
    in.push_back(tensor(mem, pe, {d.depth, 1, n, d.hidden}));
    copyOut(impl_->run(SiTVelocity, in)[0], v);
}

std::vector<float> Networks::decoderTokens(const std::vector<float> &latent) const {
    const Description &d = impl_->d;
    const int64_t n = static_cast<int64_t>(latent.size()) / d.latentDim;
    std::vector<Ort::Value> in;
    in.push_back(tensor(impl_->mem, latent, {1, n, d.latentDim}));
    std::vector<float> out;
    copyOut(impl_->run(VaeTokens, in)[0], out);
    return out;
}

Networks::Fields Networks::query(const std::vector<float> &tokens, const std::vector<float> &xy, size_t chunk) const {
    const Description &d = impl_->d;
    const int64_t n = static_cast<int64_t>(tokens.size()) / d.vaeDim;
    const size_t P = xy.size() / 2;
    Fields f;
    f.cdf.reserve(P);
    f.dcdf.reserve(P);
    f.offset.reserve(2 * P);
    f.offsetc.reserve(2 * P);
    if (chunk == 0) chunk = P;
    for (size_t s = 0; s < P; s += chunk) {
        const size_t p = std::min(chunk, P - s);
        std::vector<Ort::Value> in;
        in.push_back(tensor(impl_->mem, tokens, {1, n, d.vaeDim}));
        in.push_back(tensor(impl_->mem, xy.data() + 2 * s, 2 * p, {1, static_cast<int64_t>(p), 2}));
        const std::vector<Ort::Value> out = impl_->run(VaeQuery, in);
        appendOut(out[0], f.cdf);
        appendOut(out[1], f.dcdf);
        appendOut(out[2], f.offset);
        appendOut(out[3], f.offsetc);
    }
    return f;
}
#else
std::vector<float> Networks::geometryCode(const std::vector<float> &, const std::vector<float> &,
                                          const std::vector<float> &, const std::vector<float> &) const {
    throw std::runtime_error("this build has no ONNX Runtime");
}
std::vector<float> Networks::conditioning(const std::vector<float> &, const std::vector<float> &,
                                          const std::vector<float> &) const {
    throw std::runtime_error("this build has no ONNX Runtime");
}
void Networks::velocity(const std::vector<float> &, double, const std::vector<float> &, std::vector<float> &) const {
    throw std::runtime_error("this build has no ONNX Runtime");
}
std::vector<float> Networks::decoderTokens(const std::vector<float> &) const {
    throw std::runtime_error("this build has no ONNX Runtime");
}
Networks::Fields Networks::query(const std::vector<float> &, const std::vector<float> &, size_t) const {
    throw std::runtime_error("this build has no ONNX Runtime");
}
#endif

bool Networks::available() {
#ifdef CROSSGEN_HAVE_ONNXRUNTIME
    return true;
#else
    return false;
#endif
}

std::string Networks::defaultDirectory() {
    const char *env = std::getenv("CROSSGEN_SQUADGEN");
    return (env && *env) ? std::string(env) : std::string(CROSSGEN_SQUADGEN_DEFAULT);
}

const std::vector<std::string> &Networks::graphNames() { return kGraphs; }

}  // namespace squadgen
