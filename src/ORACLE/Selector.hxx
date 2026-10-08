#ifndef __ORACLE_SELECTOR_HXX__
#define __ORACLE_SELECTOR_HXX__

#include <memory>
#include <string>
#include <vector>

namespace oracle {

// The algorithm selector py/train_classifier.py trains and exports with --onnx
// (py/selector.onnx): from one model's inputs -- its boundary features and the
// probe method's own result -- the probability that each method produces a
// valid block decomposition, the quality of the one it would produce on the
// selector's metric, and the expected utility of choosing it. The trainer's
// docstring defines all three; everything the training preprocessed is inside
// the graph, so a row goes in exactly as py/build_dataset.py wrote it, a
// probe metric that does not exist as NaN.
//
// ### What the file says about itself
//
// An ONNX graph knows its input as a width and its outputs as rows, not which
// column or method is which. train_classifier.py's export_onnx() has written
// them into the file's metadata since 2026-10-02 (the crossgen.* keys below),
// so a model retrained on other columns, another metric or another probe is
// read as what it is. A file exported before that carries none, and is taken
// to be the trainer's defaults: py/build_dataset.py's FEATURES and probe
// columns, its five METHODS in their column order, and train_classifier.py's
// METRIC, PROBE and FALLBACK -- checked against the widths the graph declares,
// which is the most that can be checked of a file that does not say. If the
// two scripts' defaults change, kDefault* in the .cxx has to follow them.
//
// ### Why ONNX Runtime is behind this class
//
// It is a large dependency that only this file needs, so it is linked PRIVATE
// to the ORACLE library and its headers reach Selector.cxx alone, the way
// OpenCASCADE's reach src/geom alone: nothing else of ORACLE, the viewer or a
// driver names an ONNX type. A build without it still builds ORACLE, and the
// constructor says why it cannot read a selector instead.
class Selector {
public:
    // What the model was trained on and for.
    struct Description {
        std::vector<std::string> inputs;   // the input columns, in order
        std::vector<std::string> methods;  // the output rows: "zipline", "umber", ...
        std::string metric;                // the one metric it ranks on
        bool higherBetter = true;          // that metric's direction
        std::string probe;                 // run before it is asked; "" for none
        std::string fallback;              // run when nothing else was valid; "" for none
        // Read from the file's crossgen.* metadata, rather than assumed.
        bool fromFile = false;
    };

    // One model's answer, per method in Description::methods order.
    struct Prediction {
        std::vector<double> quality;   // the metric of a valid run; the probe's is its own
        std::vector<double> valid;     // the probability of one; the probe's is its result
        // Higher better. -inf for a probe that was not valid: no longer an option.
        std::vector<double> utility;
    };

    // Loads the selector at `path`. Throws std::runtime_error with the reason
    // when it cannot: no such file, a graph that is not a selector, metadata
    // that disagrees with the graph, or a build without ONNX Runtime.
    explicit Selector(const std::string &path);
    ~Selector();
    Selector(const Selector &) = delete;
    Selector &operator=(const Selector &) = delete;

    const Description &description() const;
    const std::string &path() const;

    // The prediction for one row of inputs, in description().inputs order.
    // Throws std::runtime_error on a row of the wrong length or a failed run.
    Prediction predict(const std::vector<double> &inputs) const;

    // Whether this build can read a selector at all.
    static bool available();

    // The selector a caller gets when it names none: $CROSSGEN_SELECTOR when it
    // is set, py/selector.onnx of the source tree this was built from
    // otherwise.
    static std::string defaultPath();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace oracle

#endif // __ORACLE_SELECTOR_HXX__
