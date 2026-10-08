#ifndef __ORACLE_HXX__
#define __ORACLE_HXX__

#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "ORACLE/Candidate.hxx"
#include "mesh/BlockDecomposition.hxx"
#include "mesh/BoundaryFeatures.hxx"
#include "mesh/FeatureFrame.hxx"
#include "mesh/Mesh.hxx"

// ORACLE: a block decomposition chosen rather than computed. It has no layout
// method of its own. It runs the ones this codebase has -- ZIPLINE, UMBER,
// MERIDIAN, TORSION and ATLAS -- and lets the selector py/train_classifier.py
// trained on their results (py/selector.onnx, read by oracle::Selector) decide
// which:
//
//   1  the inputs    the model's BoundaryFeatures, read as py/build_dataset.py
//                    reads them, vertices no triangle uses dropped first
//   2  the probe     the selector's probe method (ZIPLINE, about 0.1 s) run
//                    first; its validity, metrics, coverage and block count
//                    are inputs as well
//   3  the ranking   the selector's expected utility of each method on its
//                    metric (alignment_quality by default)
//   4  the rule      the probe is kept if it was valid and ranks first;
//                    otherwise the best-ranked other method runs, and the
//                    better of the two valid runs on that metric is kept;
//                    where neither is valid, the fallback (ATLAS) runs
//
// The rule is the one the selector is trained and scored for --
// train_classifier.py's "<probe> first" report row, and the "action" column
// of its --predict -- so the selector is asked once, after one run, and ORACLE
// runs at most three methods and usually one or two. Where the selector's
// answer is no prediction at all (Options::guardDomain), the probe and the
// fallback run instead and the better is kept.
//
// ### What it promises
//
// The runs are the dataset's runs (oracle::Candidate says how they are kept
// so), and a decomposition is kept only if it is valid in crossgen's sense --
// every side on dS or shared, 99.9% of the area covered. What it cannot
// promise is that the method it keeps is the best of the five: the selector
// was trained on 35 single-material models and the faces of the MAMBO parts
// (2026-10-02), and on a model unlike those its ranking is an extrapolation
// -- on any multi-material one, a wild one, which Options::guardDomain
// catches. The report keeps the whole ranking, so a surprising choice can be
// seen for what it was.
//
// ### Running it
//
// run() does all of it, and blocks: a ZIPLINE probe, then at most two more
// methods, which for ATLAS can be minutes. Afterwards the chosen method's run
// is kept -- its blocks, and its pipeline, which getChoice().mesh() meshes
// with the method's own mesher -- and the others are released.
class ORACLE {
public:
    struct Options {
        // The selector to read. Empty: oracle::Selector::defaultPath(), which
        // is $CROSSGEN_SELECTOR or else py/selector.onnx of the source tree.
        std::string selector;
        // Do not follow a ranking the selector cannot have meant. Asked about
        // a model unlike any it was trained on it does not degrade gently: an
        // input that was constant in training is standardised by a deviation
        // of 1e-6, so the 2026-10-02 selector, trained on single-material
        // models alone (feat_regions always 1, feat_interface_length always
        // 0), puts a two-material model 700,000 deviations out and predicts
        // an alignment_quality of 768 for it (multimat/geom001; on its 442
        // training rows every prediction is within [0, 1.003]). A predicted
        // quality outside [-0.5, 1.5] is taken as that, and then ORACLE runs
        // the probe and the fallback and keeps the better -- the trainer's
        // rule with no model -- and says so in the decision.
        bool guardDomain = true;
    };

    // One method run, in the order ORACLE ran them.
    struct Run {
        std::string method;
        std::string why;          // "probe", "ranked first", "fallback"
        bool valid = false;       // a decomposition of the model
        bool raised = false;      // threw, or made no blocks object
        std::string error;
        double coverage = 0.0;
        int blocks = 0;
        // The selector's metric on the run, as crossgen computes it; NaN
        // when the run was not valid.
        double metric = std::numeric_limits<double>::quiet_NaN();
        double seconds = 0.0;
    };

    struct Report {
        // The selector: the file, what it ranks on, and whether the file
        // named its columns or the trainer's defaults had to be assumed.
        std::string selector;
        bool described = false;
        std::string metric;
        bool higherBetter = true;
        std::string probe, fallback;
        std::vector<std::string> methods;     // its rows
        std::vector<std::string> inputNames;  // its columns
        std::vector<double> inputs;           // the row it was given
        // Its answer, per method: the expected utility, the probability of a
        // valid run and the metric of one; and the methods by utility, best
        // first.
        std::vector<double> utility, pValid, quality;
        std::vector<int> ranking;
        // A predicted quality outside [-0.5, 1.5]: the model is outside what
        // the selector was trained on (Options::guardDomain).
        bool outOfDomain = false;

        std::vector<Run> runs;
        int chosen = -1;            // index into runs; -1 when none could be kept
        bool valid = false;         // the kept run is a decomposition of the model
        std::string decision;       // the rule as it went, in one sentence
        int droppedVertices = 0;    // vertices no triangle used, taken out first
        double seconds = 0.0;
        // Why it stopped short of a choice, when it did.
        std::vector<std::string> messages;
    };

    explicit ORACLE(std::shared_ptr<const Mesh> mesh);
    ORACLE(std::shared_ptr<const Mesh> mesh, const Options &opts);
    ~ORACLE();

    // The whole of it. True when the run kept is a valid decomposition. False
    // when none was -- the best of them is still kept, for looking at -- or
    // when the selector could not be read or asked; getReport().messages says
    // which.
    bool run();

    // The run kept, with its pipeline: what to draw and what to mesh.
    bool hasChoice() const { return chosen_ != nullptr; }
    const oracle::Candidate &getChoice() const { return *chosen_; }
    const BlockDecomposition &getDecomposition() const { return chosen_->decomposition(); }

    const Report &getReport() const { return report_; }
    const Options &getOptions() const { return options_; }
    // The model the methods ran on: the one given, less any vertex no
    // triangle uses.
    const Mesh &getMesh() const { return *model_; }
    std::shared_ptr<const Mesh> getMeshPtr() const { return model_; }
    // Its BoundaryFeatures (default options), what the feat_* inputs are read
    // from and regularity() grades against.
    const BoundaryFeatures &getFeatures() const;
    // Its FeatureFrame, what alignment_quality() grades against.
    const FeatureFrame &getFrame() const;

    // Metric `name` -- alignment_quality, regularity, angle_quality or
    // chord_quality -- of `D` on this model, as crossgen's
    // BlockDecomposition methods compute it: NaN for an unknown name, for a
    // decomposition that does not cover, and for a value that is not finite.
    double metric(const std::string &name, const BlockDecomposition &D) const;
    static bool knownMetric(const std::string &name);

private:
    // Runs method `m` for the reason `why`, records it in report_.runs and
    // keeps the run in ran_. Returns its index into report_.runs.
    int runMethod(oracle::Method m, const char *why);
    // Input column `name` of the selector: feat_<summary field> or
    // probe_<method>_<column>, as py/build_dataset.py writes it. False for a
    // column ORACLE cannot compute.
    bool inputValue(const std::string &name, double &out) const;

    std::shared_ptr<const Mesh> model_;
    Options options_;
    Report report_;
    mutable std::unique_ptr<BoundaryFeatures> features_;
    mutable std::unique_ptr<FeatureFrame> frame_;
    // Every run of this call, by method, until the rule has chosen.
    std::vector<std::unique_ptr<oracle::Candidate>> ran_;
    std::unique_ptr<oracle::Candidate> chosen_;
};

#endif // __ORACLE_HXX__
