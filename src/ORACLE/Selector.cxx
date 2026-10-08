#include "ORACLE/Selector.hxx"

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <stdexcept>

#ifdef CROSSGEN_HAVE_ONNXRUNTIME
#include <onnxruntime_cxx_api.h>
#endif

// CMake sets it to the source tree's py/selector.onnx; this is only for a
// build that does not.
#ifndef CROSSGEN_SELECTOR_DEFAULT
#define CROSSGEN_SELECTOR_DEFAULT "py/selector.onnx"
#endif

namespace oracle {
namespace {

#ifdef CROSSGEN_HAVE_ONNXRUNTIME
// The trainer's defaults, for a file that names nothing (see the header):
// py/build_dataset.py's FEATURES as feat_*, then PROBE = "zipline"'s
// PROBE_COLUMNS -- valid, the four METRICS, coverage, num_blocks -- and its
// METHODS in the order it writes their columns; train_classifier.py's METRIC,
// PROBE and FALLBACK.
const char *const kDefaultInputs[] = {
    "feat_regions",          "feat_holes",
    "feat_euler",            "feat_corners",
    "feat_corners_one_block", "feat_corners_two_blocks",
    "feat_corners_three_blocks", "feat_corners_four_blocks",
    "feat_acute_corners",    "feat_ambiguous_corners",
    "feat_corner_defect",    "feat_singularity_bound",
    "feat_minimum_defect",   "feat_isoperimetric_ratio",
    "feat_curved_fraction",  "feat_shortest_run",
    "feat_interface_length", "probe_zipline_valid",
    "probe_zipline_alignment_quality", "probe_zipline_regularity",
    "probe_zipline_angle_quality", "probe_zipline_chord_quality",
    "probe_zipline_coverage", "probe_zipline_num_blocks",
};
const char *const kDefaultMethods[] = {"zipline", "umber", "meridian", "torsion", "atlas"};
const char *const kDefaultMetric = "alignment_quality";
const char *const kDefaultProbe = "zipline";
const char *const kDefaultFallback = "atlas";

// The output names export_onnx() gives the graph.
const char *const kOutputs[3] = {"metrics", "valid", "utility"};

// A comma-separated list, as export_onnx() writes one, each item trimmed.
std::vector<std::string> splitList(const std::string &s) {
    std::vector<std::string> out;
    size_t start = 0;
    while (start <= s.size()) {
        size_t end = s.find(',', start);
        if (end == std::string::npos) end = s.size();
        std::string item = s.substr(start, end - start);
        const size_t a = item.find_first_not_of(" \t"), b = item.find_last_not_of(" \t");
        item = (a == std::string::npos) ? std::string() : item.substr(a, b - a + 1);
        if (!item.empty()) out.push_back(item);
        start = end + 1;
    }
    return out;
}

std::string joined(const std::vector<std::string> &items) {
    std::string s;
    for (const std::string &i : items) s += (s.empty() ? "" : ", ") + i;
    return s;
}

bool contains(const std::vector<std::string> &items, const std::string &x) {
    for (const std::string &i : items)
        if (i == x) return true;
    return false;
}

// One per process, as ONNX Runtime asks, made on first use.
Ort::Env &environment() {
    static Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "ORACLE");
    return env;
}

// The width of dimension `d`, or -1 where the graph leaves it symbolic.
int64_t widthOf(const std::vector<int64_t> &shape, size_t d) {
    return d < shape.size() ? shape[d] : -1;
}
#endif

}  // namespace

struct Selector::Impl {
    std::string path;
    Description description;
#ifdef CROSSGEN_HAVE_ONNXRUNTIME
    // Run() is not const in the C++ API, though a session may run on several
    // threads at once.
    mutable Ort::Session session{nullptr};
    std::string input;
#endif
};

Selector::Selector(const std::string &path) : impl_(std::make_unique<Impl>()) {
    impl_->path = path;
    if (!std::ifstream(path).good())
        throw std::runtime_error("no selector at '" + path +
                                 "' (py/train_classifier.py --onnx writes one; $CROSSGEN_SELECTOR "
                                 "names another)");
#ifndef CROSSGEN_HAVE_ONNXRUNTIME
    throw std::runtime_error("this build has no ONNX Runtime to read the selector '" + path +
                             "' with: install it (brew install onnxruntime) and run cmake again");
#else
    Description &d = impl_->description;
    int64_t inputWidth = -1, methodRows = -1;
    try {
        Ort::SessionOptions so;
        // One row through a network of a few thousand weights: a thread pool
        // costs more than it saves.
        so.SetIntraOpNumThreads(1);
        so.SetInterOpNumThreads(1);
        impl_->session = Ort::Session(environment(), path.c_str(), so);
        const Ort::Session &s = impl_->session;
        Ort::AllocatorWithDefaultOptions alloc;

        if (s.GetInputCount() != 1)
            throw std::runtime_error("a selector has one input, this graph has " +
                                     std::to_string(s.GetInputCount()));
        impl_->input = s.GetInputNameAllocated(0, alloc).get();
        const Ort::TypeInfo inType = s.GetInputTypeInfo(0);
        const auto inInfo = inType.GetTensorTypeAndShapeInfo();
        if (inInfo.GetElementType() != ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT)
            throw std::runtime_error("its input is not float32");
        inputWidth = widthOf(inInfo.GetShape(), 1);

        bool found[3] = {false, false, false};
        for (size_t i = 0; i < s.GetOutputCount(); ++i) {
            const std::string name = s.GetOutputNameAllocated(i, alloc).get();
            for (int k = 0; k < 3; ++k) {
                if (name != kOutputs[k]) continue;
                found[k] = true;
                if (k == 2)
                    methodRows = widthOf(s.GetOutputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape(), 1);
            }
        }
        for (int k = 0; k < 3; ++k)
            if (!found[k])
                throw std::runtime_error(std::string("it has no '") + kOutputs[k] +
                                         "' output: not a py/train_classifier.py export");

        const Ort::ModelMetadata md = s.GetModelMetadata();
        auto lookup = [&](const char *key, bool &have) {
            const Ort::AllocatedStringPtr v = md.LookupCustomMetadataMapAllocated(key, alloc);
            have = v != nullptr;
            return have ? std::string(v.get()) : std::string();
        };
        bool have = false;
        const std::string inputs = lookup("crossgen.inputs", have);
        d.fromFile = have;
        if (d.fromFile) {
            bool any = false;
            d.inputs = splitList(inputs);
            d.methods = splitList(lookup("crossgen.methods", any));
            const std::vector<std::string> metrics = splitList(lookup("crossgen.metric", any));
            const std::vector<std::string> better = splitList(lookup("crossgen.better", any));
            if (metrics.size() != 1)
                throw std::runtime_error("it ranks on " + std::to_string(metrics.size()) +
                                         " metrics (" + joined(metrics) +
                                         "), and ORACLE keeps the better of two runs on one: "
                                         "retrain it with --metric");
            d.metric = metrics.front();
            d.higherBetter = better.empty() || better.front() != "lower";
            d.probe = lookup("crossgen.probe", any);
            d.fallback = lookup("crossgen.fallback", any);
        } else {
            d.inputs.assign(std::begin(kDefaultInputs), std::end(kDefaultInputs));
            d.methods.assign(std::begin(kDefaultMethods), std::end(kDefaultMethods));
            d.metric = kDefaultMetric;
            d.higherBetter = true;
            d.probe = kDefaultProbe;
            d.fallback = kDefaultFallback;
        }
    } catch (const Ort::Exception &e) {
        throw std::runtime_error("cannot read the selector '" + path + "': " + e.what());
    } catch (const std::runtime_error &e) {
        throw std::runtime_error("the selector '" + path + "': " + e.what());
    }

    // What the graph's own shapes can check. A file that names nothing is
    // only taken for the defaults if it is at least their size.
    const std::string what = d.fromFile ? "its metadata names " : "the trainer's defaults have ";
    if (inputWidth > 0 && inputWidth != static_cast<int64_t>(d.inputs.size()))
        throw std::runtime_error("the selector '" + path + "' takes " + std::to_string(inputWidth) +
                                 " inputs, and " + what + std::to_string(d.inputs.size()) +
                                 (d.fromFile ? "" : " (re-export it with py/train_classifier.py --onnx, "
                                                    "which names its columns)"));
    if (methodRows > 0 && methodRows != static_cast<int64_t>(d.methods.size()))
        throw std::runtime_error("the selector '" + path + "' ranks " + std::to_string(methodRows) +
                                 " methods, and " + what + std::to_string(d.methods.size()));
    if (d.inputs.empty() || d.methods.empty())
        throw std::runtime_error("the selector '" + path + "' names no inputs or no methods");
    if (!d.probe.empty() && !contains(d.methods, d.probe))
        throw std::runtime_error("the selector's probe '" + d.probe + "' is not one of its methods");
    if (!d.fallback.empty() && !contains(d.methods, d.fallback))
        throw std::runtime_error("the selector's fallback '" + d.fallback + "' is not one of its methods");
#endif
}

Selector::~Selector() = default;

const Selector::Description &Selector::description() const { return impl_->description; }
const std::string &Selector::path() const { return impl_->path; }

Selector::Prediction Selector::predict(const std::vector<double> &inputs) const {
    const Description &d = impl_->description;
    if (inputs.size() != d.inputs.size())
        throw std::runtime_error("the selector takes " + std::to_string(d.inputs.size()) +
                                 " inputs, given " + std::to_string(inputs.size()));
#ifndef CROSSGEN_HAVE_ONNXRUNTIME
    throw std::runtime_error("this build has no ONNX Runtime");
#else
    // float32, as the graph was traced; a NaN stays a NaN, which the graph
    // fills where the trainer allowed one (a probe's missing metrics).
    std::vector<float> row(inputs.begin(), inputs.end());
    const int64_t shape[2] = {1, static_cast<int64_t>(row.size())};
    const Ort::MemoryInfo mem = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
    Ort::Value in = Ort::Value::CreateTensor<float>(mem, row.data(), row.size(), shape, 2);
    const char *inName = impl_->input.c_str();
    std::vector<Ort::Value> out;
    try {
        out = impl_->session.Run(Ort::RunOptions{nullptr}, &inName, &in, 1, kOutputs, 3);
    } catch (const Ort::Exception &e) {
        throw std::runtime_error(std::string("the selector failed: ") + e.what());
    }

    const size_t A = d.methods.size();
    // The metrics come out (1, method, metric); one metric is what the
    // description promises, and what a second column would be is unknown.
    auto read = [&](int k, size_t width, std::vector<double> &into) {
        const size_t n = out[k].GetTensorTypeAndShapeInfo().GetElementCount();
        if (n != width)
            throw std::runtime_error(std::string("the selector's '") + kOutputs[k] + "' output has " +
                                     std::to_string(n) + " values for " + std::to_string(A) + " methods");
        const float *v = out[k].GetTensorData<float>();
        into.assign(v, v + n);
    };
    Prediction p;
    read(0, A, p.quality);
    read(1, A, p.valid);
    read(2, A, p.utility);
    return p;
#endif
}

bool Selector::available() {
#ifdef CROSSGEN_HAVE_ONNXRUNTIME
    return true;
#else
    return false;
#endif
}

std::string Selector::defaultPath() {
    const char *env = std::getenv("CROSSGEN_SELECTOR");
    return (env && *env) ? std::string(env) : std::string(CROSSGEN_SELECTOR_DEFAULT);
}

}  // namespace oracle
