// Utility to run SQUADGEN (src/SQUADGEN) -- SQuadGen's inference on the ONNX graphs `python -m squadgen2d
// export-onnx` writes -- on models, or to check it against the reference case the export wrote.
//
//   TestSQUADGEN [options] <mesh.obj> [more.obj ...]
//   TestSQUADGEN --check <dir>                       (best_models/onnx/reference)
//
// Per model: Stages 1-5 (tokens and code, noise, the ODE, the decoded fields, QUARRY) for every seed asked for, the
// blocks of each meshed by BlockQuadMesh at 0.01 x the diagonal and smoothed by TMOP as quarry.py does, and the best
// kept the way `python -m squadgen2d generate` keeps it: blocks that cover the model, then a mesh with no inverted
// element, then the highest alignment quality. A sample is `valid` when it passes the first two.
//
// Options:
//   --models <dir>      the exported graphs (default $CROSSGEN_SQUADGEN, else best_models/onnx of the source tree)
//   --samples <k>       seeds per model (default 1); sample j's noise is seed + j
//   --seed <s>          the token sets' seed and the first noise seed (default 0)
//   --solver <m>        dopri5 (default), heun or euler;  --steps <n> for the last two (default 12)
//   --atol <x>, --rtol <x>   dopri5's tolerances (default 1e-6, 1e-3)
//   --white             white initial noise rather than Taubin-smoothed
//   --h-rel <x>         the extraction mesh, as a fraction of the diagonal (default 0.002)
//   --tmop <n>          TMOP sweeps on the mesh of the blocks (default 300; 0 leaves it as meshed)
//   --threads <n>       ONNX Runtime's intra-op threads (default: its own choice)
//   --blocks-obj <f>    the kept sample's macro edges as OBJ polylines (one model only)
//
// --check runs the C++ pipeline on the reference case's inputs and compares every stage with what torch produced
// there: the frame, the truncated SDF, the Geom-AE code, the SiT conditioning and a velocity, Taubin's smoothing of
// the same white noise, the latent dopri5 reaches from the same initial noise, the decoder's tokens and fields, the
// extraction mesh, and QUARRY's blocks on the reference fields and on the ones the C++ chain decoded.
//
// The exit code is 0 when every model kept a valid sample (or every check passed), 1 when one did not, 2 on bad
// arguments, and 77 -- ctest's skip -- when there is no ONNX Runtime or no graphs to run.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "SQUADGEN/SQUADGEN.hxx"
#include "mesh/BlockQuadMesh.hxx"
#include "mesh/BoundaryFeatures.hxx"
#include "mesh/FeatureFrame.hxx"
#include "mesh/QuadMesh.hxx"
#include "mesh/TMOP.hxx"

namespace {

using Clock = std::chrono::steady_clock;
using squadgen::Networks;

double since(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

struct Settings {
    std::string models = Networks::defaultDirectory();
    int samples = 1;
    unsigned seed = 0;
    SQUADGEN::Options opts;
    int tmop = 300;
    int threads = 0;
    std::string blocksObj;
};

// ── one sample, scored as generate scores it ─────────────────────────────────

struct Sample {
    unsigned seed = 0;
    bool covers = false;
    int blocks = 0, failed = 0, quads = 0, inverted = -1;
    double minSJ = std::nan(""), regularity = std::nan(""), angle = std::nan(""), chord = std::nan(""),
           alignment = std::nan("");
    int evaluations = 0;
    double seconds = 0.0;
    bool valid() const { return covers && inverted == 0; }
};

// generate.rank: covers, then no inverted quad, then alignment quality.
bool better(const Sample &a, const Sample &b) {
    auto key = [](const Sample &s) {
        return std::make_tuple(s.covers, s.covers && s.inverted == 0, std::isfinite(s.alignment) ? s.alignment : -1.0);
    };
    return key(a) > key(b);
}

double diagonal(const Mesh &m) {
    Point lo = m.vertices[0], hi = m.vertices[0];
    for (const Point &p : m.vertices) {
        lo = {std::min(lo[0], p[0]), std::min(lo[1], p[1])};
        hi = {std::max(hi[0], p[0]), std::max(hi[1], p[1])};
    }
    return normP(hi - lo);
}

// quarry.py's extract(): the four block qualities when the blocks cover, and the mesh at 0.01 x the diagonal after
// TMOP with ZIPLINE's settings (metric 7 at the element corners).
void score(const SQUADGEN &g, int tmop, Sample &s) {
    const BlockDecomposition &D = g.decomposition();
    const QUARRY::Report &r = g.quarry().getReport();
    s.covers = g.getReport().valid;
    s.blocks = static_cast<int>(D.blocks.size());
    s.failed = r.failedCharts + r.unmatchedSides;
    if (!s.covers || D.blocks.empty()) return;
    const Mesh &model = *g.extractionMesh();
    s.regularity = D.regularity(BoundaryFeatures(model));
    s.angle = D.angleQuality();
    s.chord = D.chordQuality();
    s.alignment = D.alignmentQuality(FeatureFrame(model));
    BlockQuadMesh::Options mo;
    mo.targetEdgeLength = 0.01 * diagonal(model);
    const BlockQuadMesh grid(D, mo);
    if (grid.getReport().quads <= 0) return;
    mesh::QuadMesh q = mesh::QuadMesh::from(grid, mesh::QuadMesh::Options());
    if (tmop > 0) {
        mesh::TMOP::Options t;
        t.metric = mesh::TMOP::ShapeSize007;
        t.quadrature = mesh::TMOP::Corners;
        t.maxSweeps = tmop;
        mesh::TMOP smoother(q, t);
        smoother.run();
    } else {
        q.computeQuality();
    }
    s.quads = q.quality.quads;
    s.inverted = q.quality.invertedQuads;
    s.minSJ = q.quality.minScaledJacobian;
}

std::string fixed(double v, int digits) {
    if (!std::isfinite(v)) return "--";
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(digits) << v;
    return oss.str();
}

bool runModel(const std::shared_ptr<const Networks> &nets, const std::string &path, const Settings &s) {
    const size_t slash = path.find_last_of('/');
    const std::string name = slash == std::string::npos ? path : path.substr(slash + 1);
    const Clock::time_point t0 = Clock::now();
    SQUADGEN g(nets, squadgen::readOBJ(path), s.opts);
    g.condition();
    std::vector<Sample> samples;
    Sample best;
    bool any = false;
    BlockDecomposition kept;
    for (int k = 0; k < s.samples; ++k) {
        Sample smp;
        smp.seed = s.seed + static_cast<unsigned>(k);
        const Clock::time_point t1 = Clock::now();
        g.noise(smp.seed);
        g.sample();
        g.decode();
        g.extract();
        smp.evaluations = g.getReport().evaluations;
        score(g, s.tmop, smp);
        smp.seconds = since(t1);
        std::cout << "  " << name << " seed " << smp.seed << ": "
                  << (smp.covers ? std::to_string(smp.blocks) + " blocks, " +
                                       (smp.valid() ? std::string("valid") : std::to_string(smp.inverted) + " quads inverted") +
                                       ", min SJ " + fixed(smp.minSJ, 3) + ", alignment " + fixed(smp.alignment, 3)
                                 : "not extracted: " + (smp.failed > 0 ? std::to_string(smp.failed) +
                                                                             " charts/sides unmatched"
                                                                       : std::string("blocks not a valid decomposition")))
                  << "  (" << smp.evaluations << " evaluations, " << fixed(g.getReport().sampleSeconds, 1)
                  << " s sampling, " << fixed(smp.seconds, 1) << " s in all)\n";
        if (!any || better(smp, best)) {
            best = smp;
            kept = g.decomposition();
            any = true;
        }
        samples.push_back(smp);
    }
    int valid = 0;
    for (const Sample &x : samples) valid += x.valid();
    std::cout << std::left << std::setw(24) << name << " " << valid << "/" << samples.size() << " valid; "
              << (best.covers ? "kept seed " + std::to_string(best.seed) + ", " + std::to_string(best.blocks) +
                                    " blocks, alignment " + fixed(best.alignment, 3) + ", angle " +
                                    fixed(best.angle, 3) + ", regularity " + fixed(best.regularity, 3) +
                                    (best.valid() ? "" : " (not valid: " + std::to_string(best.inverted) +
                                                             " quads inverted)")
                              : std::string("none extracted"))
              << "  (" << g.getReport().triangles << " extraction triangles, " << fixed(since(t0), 0) << " s)\n";
    if (!s.blocksObj.empty() && best.covers) kept.writeEdgesOBJ(s.blocksObj);
    return best.valid();
}

// ── the reference case ───────────────────────────────────────────────────────

// A .npy array (format 1.0-3.0, C order, little-endian f4, f8, i4 or i8) as doubles, with its shape.
struct Array {
    std::vector<size_t> shape;
    std::vector<double> data;
    size_t rows() const { return shape.empty() ? 0 : shape[0]; }
};

Array readNpy(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("no " + path);
    char magic[6];
    in.read(magic, 6);
    if (std::string(magic, 6) != "\x93NUMPY") throw std::runtime_error(path + " is not a .npy file");
    unsigned char ver[2];
    in.read(reinterpret_cast<char *>(ver), 2);
    uint32_t hlen = 0;
    if (ver[0] == 1) {
        unsigned char b[2];
        in.read(reinterpret_cast<char *>(b), 2);
        hlen = b[0] | (b[1] << 8);
    } else {
        unsigned char b[4];
        in.read(reinterpret_cast<char *>(b), 4);
        hlen = b[0] | (b[1] << 8) | (b[2] << 16) | (static_cast<uint32_t>(b[3]) << 24);
    }
    std::string header(hlen, ' ');
    in.read(header.data(), hlen);
    auto field = [&](const std::string &key) {
        const size_t k = header.find("'" + key + "'");
        if (k == std::string::npos) throw std::runtime_error(path + ": no " + key + " in the header");
        return header.substr(header.find(':', k) + 1);
    };
    const std::string descr = field("descr");
    const std::string dtype = descr.substr(descr.find('\'') + 1, 3);
    if (field("fortran_order").find("False") == std::string::npos)
        throw std::runtime_error(path + " is in Fortran order");
    Array a;
    const std::string shp = field("shape");
    std::stringstream ss(shp.substr(shp.find('(') + 1, shp.find(')') - shp.find('(') - 1));
    std::string item;
    size_t count = 1;
    while (std::getline(ss, item, ','))
        if (item.find_first_not_of(" ") != std::string::npos) {
            a.shape.push_back(std::stoul(item));
            count *= a.shape.back();
        }
    a.data.resize(count);
    auto read = [&](auto tag) {
        using T = decltype(tag);
        std::vector<T> raw(count);
        in.read(reinterpret_cast<char *>(raw.data()), static_cast<std::streamsize>(count * sizeof(T)));
        if (!in) throw std::runtime_error(path + " is shorter than its header says");
        for (size_t i = 0; i < count; ++i) a.data[i] = static_cast<double>(raw[i]);
    };
    if (dtype == "<f4") read(float());
    else if (dtype == "<f8") read(double());
    else if (dtype == "<i4") read(int32_t());
    else if (dtype == "<i8") read(int64_t());
    else throw std::runtime_error(path + ": dtype " + dtype + " is not read here");
    return a;
}

std::vector<float> asFloat(const Array &a) { return std::vector<float>(a.data.begin(), a.data.end()); }

std::vector<Point> asPoints(const Array &a) {
    std::vector<Point> p(a.rows());
    for (size_t i = 0; i < p.size(); ++i)
        p[i] = {static_cast<double>(static_cast<float>(a.data[2 * i])),
                static_cast<double>(static_cast<float>(a.data[2 * i + 1]))};
    return p;
}

// The largest difference, and the largest difference over the largest reference value.
std::pair<double, double> difference(const std::vector<float> &got, const std::vector<double> &want) {
    if (got.size() != want.size()) return {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    double d = 0.0, m = 0.0;
    for (size_t i = 0; i < got.size(); ++i) {
        d = std::max(d, std::fabs(got[i] - want[i]));
        m = std::max(m, std::fabs(want[i]));
    }
    return {d, d / std::max(m, 1e-30)};
}

struct Checker {
    int failures = 0;
    void line(const std::string &what, double value, double tol, const std::string &unit, const std::string &note = "") {
        const bool ok = value <= tol;
        failures += ok ? 0 : 1;
        std::cout << "  " << std::left << std::setw(46) << what << std::right << std::setw(10) << std::scientific
                  << std::setprecision(2) << value << std::defaultfloat << " " << std::left << std::setw(16)
                  << unit + " (<= " + fixed(tol, 6).substr(0, 8) + ")" << (ok ? " PASS" : " FAIL")
                  << (note.empty() ? "" : "  " + note) << "\n";
    }
    void yes(const std::string &what, bool ok, const std::string &note) {
        failures += ok ? 0 : 1;
        std::cout << "  " << std::left << std::setw(46) << what << " " << std::setw(27) << note
                  << (ok ? " PASS" : " FAIL") << "\n";
    }
};

int check(const std::string &dir, const Settings &s) {
    std::map<std::string, std::string> kv;
    {
        std::ifstream in(dir + "/case.txt");
        if (!in) {
            std::cerr << "no reference case in " << dir << " (python -m squadgen2d export-onnx writes one)\n";
            return 77;
        }
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            const size_t sp = line.find(' ');
            kv[line.substr(0, sp)] = sp == std::string::npos ? "" : line.substr(sp + 1);
        }
    }
    const std::string graphs = dir + "/..";
    std::shared_ptr<const Networks> nets;
    try {
        Networks::Options no;
        no.threads = s.threads;
        nets = std::make_shared<const Networks>(graphs, no);
    } catch (const std::exception &e) {
        std::cerr << "cannot load the networks: " << e.what() << "\n";
        return 77;
    }
    const Networks::Description &d = nets->description();
    std::cout << "SQuadGen graphs " << graphs << ": " << d.model << " step " << d.step << ", N = " << d.tokens
              << ", export " << d.exportId << "\nreference " << kv["source"] << ", seed " << kv["seed"] << ", "
              << kv["nfe"] << " evaluations in torch\n";
    auto A = [&](const std::string &k) { return readNpy(dir + "/" + k + ".npy"); };
    Checker c;

    // The frame.
    const squadgen::Model model = squadgen::readOBJ(dir + "/" + kv["obj"]);
    SQUADGEN::Options opts = s.opts;
    opts.hRel = std::stod(kv["h_rel"]);
    opts.sampler.atol = std::stod(kv["atol"]);
    opts.sampler.rtol = std::stod(kv["rtol"]);
    opts.sampler.method = squadgen::Sampler::Method::Dopri5;
    SQUADGEN g(nets, model, opts);
    {
        std::stringstream cs(kv["centre"]);
        double cx = 0.0, cy = 0.0;
        cs >> cx >> cy;
        const squadgen::Frame &f = g.frame();
        const double scale = std::stod(kv["scale"]);
        c.line("frame: centre and scale", std::max({std::fabs(f.centre[0] - cx), std::fabs(f.centre[1] - cy),
                                                    std::fabs(f.scale - scale) / scale}), 1e-12, "abs/rel");
    }

    // Stage 1 on the reference token sets.
    const Array xyz = A("xyz"), xyzKV = A("xyz_kv"), sdf = A("sdf"), sdfKV = A("sdf_kv"), code = A("code");
    const std::vector<Point> tokens = asPoints(xyz), context = asPoints(xyzKV);
    {
        const auto e1 = difference(g.domain().insideSDF(tokens, d.tau), sdf.data);
        const auto e2 = difference(g.domain().insideSDF(context, d.tau), sdfKV.data);
        c.line("truncated SDF at the N and 4N points", std::max(e1.first, e2.first), 1e-5, "abs");
    }
    const std::vector<float> refXYZ = asFloat(xyz), refSDF = asFloat(sdf), refCode = asFloat(code);
    c.line("Geom-AE code, from the reference inputs",
           difference(nets->geometryCode(refXYZ, refSDF, asFloat(xyzKV), asFloat(sdfKV)), code.data).second, 1e-4,
           "of max");
    const std::vector<float> pe = nets->conditioning(refCode, refXYZ, refSDF);
    {
        const Array head = A("pe_head");   // (depth, 64, hidden)
        const size_t n = static_cast<size_t>(d.tokens), H = static_cast<size_t>(d.hidden), K = head.shape[1];
        std::vector<float> mine;
        for (size_t l = 0; l < static_cast<size_t>(d.depth); ++l)
            for (size_t i = 0; i < K; ++i)
                for (size_t h = 0; h < H; ++h) mine.push_back(pe[(l * n + i) * H + h]);
        c.line("SiT conditioning (first 64 tokens)", difference(mine, head.data).second, 1e-4, "of max");
    }
    {
        std::vector<float> v;
        nets->velocity(asFloat(A("z0")), 0.0, pe, v);
        c.line("velocity at the initial noise, t = 0", difference(v, A("v0").data).second, 1e-4, "of max");
    }

    // Stage 2: the same white noise, smoothed here.
    {
        std::vector<float> z = asFloat(A("noise"));
        squadgen::taubinSmooth(z, d.latentDim, tokens, g.domain().segments());
        c.line("Taubin smoothing of the same white noise", difference(z, A("z0").data).first, 1e-3, "abs");
    }

    // Stage 3 from the reference's initial latent, on Stage 1 recomputed here.
    g.condition(tokens, context);
    c.line("Geom-AE code, Stage 1 here", difference(g.code(), code.data).second, 1e-4, "of max");
    g.setNoise(asFloat(A("z0")));
    const Clock::time_point t0 = Clock::now();
    g.sample();
    const double seconds = since(t0);
    const int nfe = std::stoi(kv["nfe"]);
    c.line("latent dopri5 reaches (normalised)", difference(g.latent(), A("x1").data).first, 1e-2, "abs",
           std::to_string(g.getReport().evaluations) + " evaluations, torch " + std::to_string(nfe) + "; " +
               fixed(seconds, 1) + " s");
    c.yes("dopri5's evaluations within one step of torch's", std::abs(g.getReport().evaluations - nfe) <= 6,
          std::to_string(g.getReport().evaluations) + " vs " + std::to_string(nfe));

    // Stage 4: the decoder from the reference latent, at the reference mesh's centroids; then the extraction mesh.
    const Array latent = A("latent"), V = A("V"), T = A("T");
    const std::vector<float> tok = nets->decoderTokens(asFloat(latent));
    c.line("decoder tokens of the reference latent", difference(tok, A("tokens").data).second, 1e-4, "of max");
    squadgen::Model ref;
    for (size_t i = 0; i < V.rows(); ++i) ref.vertices.push_back({V.data[2 * i], V.data[2 * i + 1]});
    for (size_t t = 0; t < T.rows(); ++t)
        ref.triangles.push_back({static_cast<int>(T.data[3 * t]), static_cast<int>(T.data[3 * t + 1]),
                                 static_cast<int>(T.data[3 * t + 2])});
    const Networks::Fields f = nets->query(tok, squadgen::centroidsInFrame(ref, g.frame()));
    {
        double worst = 0.0;
        for (const char *k : {"cdf", "dcdf", "offset", "offsetc"}) {
            const std::vector<float> &mine = std::string(k) == "cdf" ? f.cdf
                                             : std::string(k) == "dcdf" ? f.dcdf
                                             : std::string(k) == "offset" ? f.offset : f.offsetc;
            worst = std::max(worst, difference(mine, A(k).data).first);
        }
        c.line("cdf, dcdf, offset, offsetc of the same", worst, 1e-3, "abs",
               std::to_string(f.size()) + " centroids");
    }
    int levels = 0;
    const squadgen::Model em = squadgen::extractionMesh(model, opts.hRel, &levels);
    c.yes("extraction mesh: the same vertices and triangles",
          em.vertices == ref.vertices && em.triangles == ref.triangles,
          std::to_string(em.triangles.size()) + " triangles, " + std::to_string(levels) + " levels");

    // Stage 5: QUARRY on the reference fields, which must give exactly what quarry.py gave on them, and on the fields
    // the C++ chain decoded from its own latent.
    const int blocks = std::stoi(kv["quarry_blocks"]);
    const bool covers = kv["quarry_covers"] == "1";
    g.setLatent(asFloat(A("x1")));
    g.decode();
    g.extract();
    c.yes("QUARRY on the reference latent's fields", g.getReport().blocks == blocks && g.getReport().valid == covers,
          std::to_string(g.getReport().blocks) + " blocks" + (g.getReport().valid ? ", covers" : "") + " vs " +
              std::to_string(blocks) + (covers ? ", covers" : ""));
    g.setNoise(asFloat(A("z0")));
    g.sample();
    g.decode();
    g.extract();
    c.yes("QUARRY on the C++ chain's fields", g.getReport().blocks == blocks && g.getReport().valid == covers,
          std::to_string(g.getReport().blocks) + " blocks" + (g.getReport().valid ? ", covers" : "") + " vs " +
              std::to_string(blocks) + (covers ? ", covers" : ""));
    std::cout << (c.failures ? "check FAILED: " + std::to_string(c.failures) + " of the comparisons\n"
                             : "check passed\n");
    return c.failures ? 1 : 0;
}

}  // namespace

int main(int argc, char **argv) {
    Settings s;
    std::vector<std::string> files;
    std::string checkDir;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto value = [&]() -> std::string {
            if (i + 1 >= argc) {
                std::cerr << a << " needs a value\n";
                std::exit(2);
            }
            return argv[++i];
        };
        if (a == "--check") checkDir = value();
        else if (a == "--models") s.models = value();
        else if (a == "--samples") s.samples = std::max(1, std::atoi(value().c_str()));
        else if (a == "--seed") s.seed = static_cast<unsigned>(std::atoi(value().c_str()));
        else if (a == "--solver") {
            const std::string m = value();
            if (m == "dopri5") s.opts.sampler.method = squadgen::Sampler::Method::Dopri5;
            else if (m == "heun") s.opts.sampler.method = squadgen::Sampler::Method::Heun;
            else if (m == "euler") s.opts.sampler.method = squadgen::Sampler::Method::Euler;
            else {
                std::cerr << "--solver takes dopri5, heun or euler\n";
                return 2;
            }
        } else if (a == "--steps") s.opts.sampler.steps = std::atoi(value().c_str());
        else if (a == "--atol") s.opts.sampler.atol = std::atof(value().c_str());
        else if (a == "--rtol") s.opts.sampler.rtol = std::atof(value().c_str());
        else if (a == "--white") s.opts.taubin = false;
        else if (a == "--h-rel") s.opts.hRel = std::atof(value().c_str());
        else if (a == "--tmop") s.tmop = std::atoi(value().c_str());
        else if (a == "--threads") s.threads = std::atoi(value().c_str());
        else if (a == "--blocks-obj") s.blocksObj = value();
        else if (!a.empty() && a[0] == '-') {
            std::cerr << "unknown option " << a << "\n";
            return 2;
        } else files.push_back(a);
    }
    s.opts.seed = s.seed;
    if (!checkDir.empty()) {
        try {
            return check(checkDir, s);
        } catch (const std::exception &e) {
            std::cerr << "check failed: " << e.what() << "\n";
            return 1;
        }
    }
    if (files.empty()) {
        std::cerr << "usage: " << argv[0] << " [--models dir] [--samples k] [--seed s] [--solver m] [--steps n] "
                  << "[--atol x] [--rtol x] [--white] [--h-rel x] [--tmop n] [--threads n] [--blocks-obj f] "
                  << "<mesh.obj>...\n       " << argv[0] << " --check <reference dir>\n";
        return 2;
    }
    std::shared_ptr<const Networks> nets;
    try {
        Networks::Options no;
        no.threads = s.threads;
        const Clock::time_point t0 = Clock::now();
        nets = std::make_shared<const Networks>(s.models, no);
        const Networks::Description &d = nets->description();
        std::cout << "SQuadGen " << d.model << " (step " << d.step << ", N = " << d.tokens << ") from " << s.models
                  << ", loaded in " << fixed(since(t0), 1) << " s; " << squadgen::Sampler::name(s.opts.sampler.method)
                  << (s.opts.taubin ? ", Taubin noise" : ", white noise") << ", h_rel " << s.opts.hRel << "\n";
    } catch (const std::exception &e) {
        std::cerr << "cannot load the networks: " << e.what() << "\n";
        return 77;
    }
    int failures = 0;
    for (const std::string &f : files) {
        try {
            failures += runModel(nets, f, s) ? 0 : 1;
        } catch (const std::exception &e) {
            std::cout << f << ": failed: " << e.what() << "\n";
            ++failures;
        }
    }
    return failures ? 1 : 0;
}
