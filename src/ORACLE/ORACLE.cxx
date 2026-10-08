#include "ORACLE/ORACLE.hxx"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "ORACLE/Selector.hxx"

namespace {

using Clock = std::chrono::steady_clock;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// feat_<key>: Mesh.boundary_features()'s keys, which are BoundaryFeatures::
// Summary's fields in snake_case, and what py/build_dataset.py writes as
// float(value).
bool featureValue(const BoundaryFeatures::Summary &s, const std::string &key, double &out) {
    const std::pair<const char *, double> table[] = {
        {"regions", s.regions},
        {"holes", s.holes},
        {"euler", s.euler},
        {"corners", s.corners},
        {"corners_one_block", s.cornersOneBlock},
        {"corners_two_blocks", s.cornersTwoBlocks},
        {"corners_three_blocks", s.cornersThreeBlocks},
        {"corners_four_blocks", s.cornersFourBlocks},
        {"acute_corners", s.acuteCorners},
        {"ambiguous_corners", s.ambiguousCorners},
        {"corner_defect", s.cornerDefect},
        {"singularity_bound", s.singularityBound},
        {"minimum_defect", s.minimumDefect},
        {"isoperimetric_ratio", s.isoperimetricRatio},
        {"curved_fraction", s.curvedFraction},
        {"shortest_run", s.shortestRun},
        {"interface_length", s.interfaceLength},
        {"area", s.area},
        {"perimeter", s.perimeter},
    };
    for (const auto &kv : table) {
        if (key == kv.first) {
            out = kv.second;
            return true;
        }
    }
    return false;
}

std::string number(double v) {
    if (std::isnan(v)) return "--";
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << v;
    return oss.str();
}

}  // namespace

ORACLE::ORACLE(std::shared_ptr<const Mesh> mesh) : ORACLE(std::move(mesh), Options()) {}

ORACLE::ORACLE(std::shared_ptr<const Mesh> mesh, const Options &opts) : options_(opts) {
    if (!mesh) throw std::invalid_argument("ORACLE: no mesh");
    // py/build_dataset.py drops the vertices no triangle uses before any
    // method runs -- mesh::Mesh keeps every vertex of the file, and ones in no
    // triangle change what UMBER and ATLAS compute -- so the runs the selector
    // learnt from were on that mesh, and these are too. The order of the
    // vertices left is kept.
    std::vector<int> index(mesh->vertices.size(), -1);
    for (const Triangle &t : mesh->triangles)
        for (const int v : t)
            if (v >= 0 && v < static_cast<int>(index.size())) index[v] = 0;
    int used = 0;
    for (int &i : index)
        if (i == 0) i = used++;
    report_.droppedVertices = static_cast<int>(mesh->vertices.size()) - used;
    if (report_.droppedVertices == 0) {
        model_ = std::move(mesh);
        return;
    }
    std::vector<Point> points;
    points.reserve(static_cast<size_t>(used));
    for (size_t v = 0; v < index.size(); ++v)
        if (index[v] >= 0) points.push_back(mesh->vertices[v]);
    std::vector<Triangle> triangles;
    triangles.reserve(mesh->triangles.size());
    for (const Triangle &t : mesh->triangles) {
        Triangle abc = {index[t[0]], index[t[1]], index[t[2]]};
        // Counter-clockwise, as crossgen.Mesh() and the .obj reader leave
        // every triangle.
        if (cross2(points[abc[1]] - points[abc[0]], points[abc[2]] - points[abc[0]]) < 0.0)
            std::swap(abc[1], abc[2]);
        triangles.push_back(abc);
    }
    model_ = std::make_shared<const Mesh>(points, triangles, mesh->triangleMatId);
}

ORACLE::~ORACLE() = default;

const BoundaryFeatures &ORACLE::getFeatures() const {
    if (!features_) features_ = std::make_unique<BoundaryFeatures>(*model_);
    return *features_;
}

const FeatureFrame &ORACLE::getFrame() const {
    if (!frame_) frame_ = std::make_unique<FeatureFrame>(*model_);
    return *frame_;
}

bool ORACLE::knownMetric(const std::string &name) {
    return name == "alignment_quality" || name == "regularity" || name == "angle_quality" ||
           name == "chord_quality";
}

double ORACLE::metric(const std::string &name, const BlockDecomposition &D) const {
    double v = kNaN;
    if (name == "alignment_quality") v = D.alignmentQuality(getFrame());
    else if (name == "regularity") v = D.regularity(getFeatures());
    else if (name == "angle_quality") v = D.angleQuality();
    else if (name == "chord_quality") v = D.chordQuality();
    return std::isfinite(v) ? v : kNaN;
}

int ORACLE::runMethod(oracle::Method m, const char *why) {
    std::unique_ptr<oracle::Candidate> c = oracle::Candidate::run(m, *model_);
    Run r;
    r.method = c->name();
    r.why = why;
    r.valid = c->valid();
    r.raised = c->raised();
    r.error = c->error();
    r.coverage = c->coverage();
    r.blocks = static_cast<int>(c->decomposition().blocks.size());
    r.seconds = c->seconds();
    if (r.valid) r.metric = metric(report_.metric, c->decomposition());
    report_.runs.push_back(r);
    ran_[static_cast<int>(m)] = std::move(c);
    return static_cast<int>(report_.runs.size()) - 1;
}

bool ORACLE::inputValue(const std::string &name, double &out) const {
    const std::string feat = "feat_", probe = "probe_";
    if (name.compare(0, feat.size(), feat) == 0)
        return featureValue(getFeatures().summary(), name.substr(feat.size()), out);
    if (name.compare(0, probe.size(), probe) != 0) return false;

    // probe_<method>_<column>: the probe's one run, as build_dataset.py's
    // row_of() writes it. A run that raised is not valid, has coverage 0 and
    // no block count; a metric exists only for a valid run. Where a value
    // does not exist it goes in as NaN, which the graph fills itself.
    const std::string rest = name.substr(probe.size());
    const std::string prefix = report_.probe + "_";
    if (report_.probe.empty() || rest.compare(0, prefix.size(), prefix) != 0) return false;
    oracle::Method m;
    if (!oracle::methodNamed(report_.probe, m) || !ran_[static_cast<int>(m)]) return false;
    const oracle::Candidate &c = *ran_[static_cast<int>(m)];
    const std::string column = rest.substr(prefix.size());
    if (column == "valid") {
        out = c.valid() ? 1.0 : 0.0;
    } else if (column == "coverage") {
        out = c.raised() ? 0.0 : c.coverage();
    } else if (column == "num_blocks") {
        out = c.raised() ? kNaN : static_cast<double>(c.decomposition().blocks.size());
    } else if (knownMetric(column)) {
        out = c.valid() ? metric(column, c.decomposition()) : kNaN;
    } else {
        return false;
    }
    return true;
}

bool ORACLE::run() {
    const Clock::time_point t0 = Clock::now();
    const int dropped = report_.droppedVertices;
    report_ = Report();
    report_.droppedVertices = dropped;
    chosen_.reset();
    ran_.clear();
    ran_.resize(oracle::kNumMethods);
    // Each method runs at most once, and the rule holds pointers into this
    // while it adds to it.
    report_.runs.reserve(oracle::kNumMethods);
    auto stop = [&](const std::string &why) {
        report_.messages.push_back(why);
        ran_.clear();
        report_.seconds = std::chrono::duration<double>(Clock::now() - t0).count();
        return false;
    };

    // The selector, and what it was trained on.
    report_.selector = options_.selector.empty() ? oracle::Selector::defaultPath() : options_.selector;
    std::unique_ptr<oracle::Selector> selector;
    try {
        selector = std::make_unique<oracle::Selector>(report_.selector);
    } catch (const std::exception &e) {
        return stop(e.what());
    }
    const oracle::Selector::Description &d = selector->description();
    report_.described = d.fromFile;
    report_.metric = d.metric;
    report_.higherBetter = d.higherBetter;
    report_.probe = d.probe;
    report_.fallback = d.fallback;
    report_.methods = d.methods;
    report_.inputNames = d.inputs;
    if (!knownMetric(d.metric))
        return stop("the selector ranks on '" + d.metric + "', which ORACLE cannot compute");
    const int A = static_cast<int>(d.methods.size());
    std::vector<oracle::Method> methods(static_cast<size_t>(A));
    for (int i = 0; i < A; ++i)
        if (!oracle::methodNamed(d.methods[i], methods[i]))
            return stop("the selector ranks '" + d.methods[i] + "', which is not a method ORACLE can run");
    auto indexOf = [&](const std::string &name) {
        for (int i = 0; i < A; ++i)
            if (d.methods[i] == name) return i;
        return -1;
    };
    const int probe = indexOf(d.probe), fallback = indexOf(d.fallback);

    // The inputs: the boundary features and the probe's run.
    const int probeRun = probe >= 0 ? runMethod(methods[probe], "probe") : -1;
    report_.inputs.assign(d.inputs.size(), kNaN);
    for (size_t i = 0; i < d.inputs.size(); ++i)
        if (!inputValue(d.inputs[i], report_.inputs[i]))
            return stop("the selector asks for '" + d.inputs[i] + "', which ORACLE cannot compute");

    // The ranking. A NaN utility -- which the graph does not produce -- would
    // rank last.
    try {
        const oracle::Selector::Prediction p = selector->predict(report_.inputs);
        report_.utility = p.utility;
        report_.pValid = p.valid;
        report_.quality = p.quality;
    } catch (const std::exception &e) {
        return stop(e.what());
    }
    auto utilityOf = [&](int i) {
        const double u = report_.utility[i];
        return std::isnan(u) ? -std::numeric_limits<double>::infinity() : u;
    };
    report_.ranking.resize(static_cast<size_t>(A));
    std::iota(report_.ranking.begin(), report_.ranking.end(), 0);
    std::stable_sort(report_.ranking.begin(), report_.ranking.end(),
                     [&](int a, int b) { return utilityOf(a) > utilityOf(b); });

    // Whether the answer is a prediction at all (Options::guardDomain): the
    // four metrics ORACLE knows are qualities in [0, 1], and on its training
    // rows the selector predicts them within [0, 1.003].
    int wild = -1;
    for (int i = 0; i < A; ++i) {
        const double q = report_.quality[i];
        if (i != probe && !(q >= -0.5 && q <= 1.5)) {
            if (wild < 0 || std::fabs(q) > std::fabs(report_.quality[wild])) wild = i;
        }
    }
    report_.outOfDomain = wild >= 0;
    const bool guarded = report_.outOfDomain && options_.guardDomain;

    // The rule.
    const int top = report_.ranking.front();
    auto better = [&](double a, double b) { return d.higherBetter ? a > b : a < b; };
    std::ostringstream why;
    if (guarded)
        why << "the selector is outside what it was trained on here (it predicts " << d.metric << " "
            << number(report_.quality[wild]) << " for " << d.methods[wild]
            << "), so its ranking is not followed: ";
    int keep = -1;
    if (!guarded && probe >= 0 && report_.runs[probeRun].valid && top == probe) {
        keep = probeRun;
        why << d.probe << " was valid and ranks first: kept it";
    } else {
        // The method run second. Followed, the selector's: the best-ranked
        // method other than the probe, which is the first of the ranking
        // unless the probe was valid and first (kept above), since the graph
        // puts a probe that was not valid last. Not followed, the fallback.
        int next = -1;
        if (guarded) {
            if (fallback >= 0 && fallback != probe) next = fallback;
        } else {
            for (const int i : report_.ranking) {
                if (i != probe) {
                    next = i;
                    break;
                }
            }
        }
        const char *label = guarded ? "fallback" : next == top ? "ranked first" : "ranked next";
        const int nextRun = next >= 0 ? runMethod(methods[next], label) : -1;
        const Run *a = probeRun >= 0 ? &report_.runs[probeRun] : nullptr;
        const Run *b = nextRun >= 0 ? &report_.runs[nextRun] : nullptr;
        const std::string nextName = next >= 0 ? d.methods[next] : std::string("nothing");
        const std::string lead = guarded ? "ran the fallback, " + nextName : nextName + " ranks first: ran it";
        if (a && a->valid && b && b->valid) {
            // A tie keeps the probe: the second run bought nothing.
            keep = better(b->metric, a->metric) ? nextRun : probeRun;
            const Run &k = report_.runs[keep], &other = report_.runs[keep == probeRun ? nextRun : probeRun];
            why << lead << ", and kept " << k.method << ", the better on " << d.metric << " ("
                << number(k.metric) << " against " << number(other.metric) << " for " << other.method << ")";
        } else if (b && b->valid) {
            keep = nextRun;
            why << lead << " and kept it";
            if (a) why << " (" << d.probe << " was not valid)";
        } else if (a && a->valid) {
            keep = probeRun;
            if (b) why << lead << ", which was not valid: kept " << d.probe;
            else why << "kept " << d.probe;
        } else if (guarded || !b) {
            // Nothing left to run: the fallback was the second run, or there
            // is no method but the probe.
            if (a && b) why << "neither " << d.probe << " nor " << nextName << " was valid";
            else if (a) why << d.probe << " was not valid";
            else if (b) why << nextName << " was not valid";
            else why << "nothing ran";
        } else {
            why << (a ? "neither " + d.probe + " nor " + nextName + " was valid"
                      : nextName + " ranks first but was not valid");
            if (fallback >= 0 && fallback != next && fallback != probe) {
                const int fb = runMethod(methods[fallback], "fallback");
                why << ": ran the fallback, " << d.fallback;
                if (report_.runs[fb].valid) {
                    keep = fb;
                    why << ", and kept it";
                } else {
                    why << ", which was not valid either";
                }
            } else if (fallback >= 0) {
                why << ", and " << d.fallback << " is the fallback";
            }
        }
    }

    // Nothing valid: keep whichever run covers the most, so there is still
    // something to look at, and say so.
    report_.valid = keep >= 0;
    if (keep < 0) {
        double best = -1.0;
        for (int r = 0; r < static_cast<int>(report_.runs.size()); ++r) {
            const Run &run = report_.runs[r];
            if (run.blocks > 0 && run.coverage > best) {
                best = run.coverage;
                keep = r;
            }
        }
        if (keep >= 0)
            why << "; no valid decomposition, so " << report_.runs[keep].method << "'s blocks ("
                << std::fixed << std::setprecision(1) << 100.0 * report_.runs[keep].coverage
                << "% covered) are kept to look at";
        else
            why << "; no method produced any blocks";
    }
    report_.chosen = keep;
    if (keep >= 0) {
        oracle::Method m;
        if (oracle::methodNamed(report_.runs[keep].method, m)) chosen_ = std::move(ran_[static_cast<int>(m)]);
    }
    // The others' pipelines go here; only the kept one is meshed.
    ran_.clear();
    report_.decision = why.str();
    report_.seconds = std::chrono::duration<double>(Clock::now() - t0).count();
    return report_.valid;
}
