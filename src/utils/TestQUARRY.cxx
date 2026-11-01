// Utility to run QUARRY (src/QUARRY) on the true chart distance fields of a method's own layout: the round trip of
// docs/squadgen_2d.md Sec. 6 step 1, entirely in C++ and with no network in it, so that what it measures is the
// extractor and nothing else.
//
//   TestQUARRY [options] <mesh.obj> [more.obj ...]
//   TestQUARRY --selftest
//
// Per model: ZIPLINE's blocks, their ChartDistanceField (the fields SQ-VAE is trained on, mesh/ChartDistanceField),
// the model mesh refined by midpoint subdivision until its median edge is at most 1.5 h_rel times the bounding-box
// diagonal (what squadgen2d's round trip extracts on), the four fields QUARRY reads at every refined triangle's
// centroid, QUARRY on them, and the comparison: every true chart matched by one found chart with area IoU above
// 1/2, no found chart left over, the decomposition covering the model, and BlockQuadMesh's mesh of it valid with no
// inverted element.
//
// Options:
//   --h-rel <x>        extraction resolution (default 0.004; squadgen2d's round trip uses 0.002)
//   --noise <a>        perturb the fields before extracting, as a decoder would: a smooth field of amplitude a on
//                      CDF and DCDF (clamped to [0, 1]) and of a/10 of the diagonal on both offsets, a few random
//                      waves across the model, plus white noise of a/3 per triangle
//   --seed <n>         the noise's seed (default 1)
//   --blocks-obj <f>   the extracted decomposition's macro edges as OBJ polylines (one model only)
//
// --selftest runs synthetic domains -- a rectangle, an L, a square round a square hole, a disk, a plate with a round
// hole -- each on its true fields and on fields with noise 0.15, and requires every one recovered.
//
// The exit code is 0 when every model was recovered with a valid mesh, 1 when one was not, 2 on bad arguments.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "QUARRY/QUARRY.hxx"
#include "ZIPLINE/TriangleLocator.hxx"
#include "ZIPLINE/ZIPLINE.hxx"
#include "mesh/BlockQuadMesh.hxx"
#include "mesh/ChartDistanceField.hxx"
#include "mesh/Mesh.hxx"
#include "triangle/TriangleMesher.hpp"

namespace {

struct Settings {
    double hRel = 0.004;
    double noise = 0.0;
    unsigned seed = 1;
    std::string blocksObj;
};

struct Outcome {
    int trueCharts = 0, found = 0, matched = 0, blocks = 0;
    bool covers = false, meshValid = false, recovered = false;
    double seconds = 0.0;
    std::string note;
};

double diagonal(const Mesh &m) {
    Point lo = m.vertices[0], hi = m.vertices[0];
    for (const Point &p : m.vertices) {
        lo = {std::min(lo[0], p[0]), std::min(lo[1], p[1])};
        hi = {std::max(hi[0], p[0]), std::max(hi[1], p[1])};
    }
    return normP(hi - lo);
}

// One 1-to-4 midpoint subdivision, every original vertex and every triangle's orientation kept.
Mesh subdivide(const Mesh &m) {
    std::vector<Point> V = m.vertices;
    std::map<std::pair<int, int>, int> mid;
    auto at = [&](int a, int b) {
        const std::pair<int, int> key{std::min(a, b), std::max(a, b)};
        const auto it = mid.find(key);
        if (it != mid.end()) return it->second;
        V.push_back((m.vertices[a] + m.vertices[b]) * 0.5);
        return mid[key] = static_cast<int>(V.size()) - 1;
    };
    std::vector<Triangle> T;
    for (const Triangle &t : m.triangles) {
        const int ab = at(t[0], t[1]), bc = at(t[1], t[2]), ca = at(t[2], t[0]);
        T.push_back({t[0], ab, ca});
        T.push_back({ab, t[1], bc});
        T.push_back({ca, bc, t[2]});
        T.push_back({ab, bc, ca});
    }
    return Mesh(V, T);
}

// The extraction mesh: subdivided until its median edge is at most 1.5 hRel diagonals, four times at most -- the
// rule of squadgen2d.roundtrip.extraction_mesh.
Mesh extractionMesh(const Mesh &model, double hRel) {
    Mesh m(model.vertices, model.triangles);
    const double target = 1.5 * hRel * diagonal(model);
    for (int level = 0; level < 4; ++level) {
        std::vector<double> e;
        for (const auto &ed : m.edges) e.push_back(normP(m.vertices[ed[1]] - m.vertices[ed[0]]));
        std::nth_element(e.begin(), e.begin() + e.size() / 2, e.end());
        if (e[e.size() / 2] <= target) break;
        m = subdivide(m);
    }
    return m;
}

// A smooth scalar field: a few plane waves of random direction and a wavelength between a third and the whole of
// the model's diagonal, scaled to amplitude about 1.
struct Waves {
    std::vector<Point> k;
    std::vector<double> phase;
    Waves(std::mt19937 &rng, double diag) {
        std::uniform_real_distribution<double> u(0.0, 1.0);
        for (int i = 0; i < 4; ++i) {
            const double angle = 2.0 * M_PI * u(rng), wavelength = diag * (0.33 + 0.67 * u(rng));
            k.push_back(Point{std::cos(angle), std::sin(angle)} * (2.0 * M_PI / wavelength));
            phase.push_back(2.0 * M_PI * u(rng));
        }
    }
    double operator()(const Point &p) const {
        double s = 0.0;
        for (size_t i = 0; i < k.size(); ++i) s += std::sin(dotP(k[i], p) + phase[i]);
        return s / std::sqrt(static_cast<double>(k.size()));
    }
};

Outcome roundTrip(const std::shared_ptr<Mesh> &model, const Settings &s, const std::string &name) {
    Outcome out;
    ZIPLINE zipline(model);
    zipline.run(ZIPLINE::Stage::Blocks);
    if (!zipline.hasBlocks() || !zipline.getBlocks().decomposition().covers()) {
        out.note = "ZIPLINE gave no covering layout";
        return out;
    }
    const BlockDecomposition &truth = zipline.getBlocks().decomposition();
    ChartDistanceField::Options co;
    co.samples = 16;
    co.fpsMax = 0;
    const ChartDistanceField cdf(truth, co);
    out.trueCharts = static_cast<int>(truth.blocks.size());

    // The split mesh as a Mesh, for locating points in it.
    std::vector<Triangle> splitTris;
    for (int t = 0; t < cdf.triangleCount(); ++t) {
        const std::array<int, 3> v = cdf.triangle(t);
        splitTris.push_back({v[0], v[1], v[2]});
    }
    const Mesh split(cdf.vertices(), splitTris);
    const TriangleLocator locate(split);
    std::vector<Point> splitCentroid(splitTris.size());
    for (size_t t = 0; t < splitTris.size(); ++t)
        splitCentroid[t] = (split.vertices[splitTris[t][0]] + split.vertices[splitTris[t][1]] +
                            split.vertices[splitTris[t][2]]) / 3.0;

    auto R = std::make_shared<Mesh>(extractionMesh(*model, s.hRel));
    const int nt = static_cast<int>(R->triangles.size());
    const double diag = diagonal(*R);
    std::mt19937 rng(s.seed);
    const Waves w0(rng, diag), w1(rng, diag), w2(rng, diag), w3(rng, diag), w4(rng, diag), w5(rng, diag);
    std::uniform_real_distribution<double> white(-1.0, 1.0);
    QUARRY::Fields f;
    std::vector<int> trueChart(nt, -1);
    std::vector<double> area(nt);
    for (int t = 0; t < nt; ++t) {
        const Point &a = R->vertices[R->triangles[t][0]], &b = R->vertices[R->triangles[t][1]],
                    &c = R->vertices[R->triangles[t][2]];
        const Point p = (a + b + c) / 3.0;
        area[t] = 0.5 * std::fabs(cross2(b - a, c - a));
        int st = locate.locate(p);
        if (st < 0) {   // just outside a chord of a curved dS: the nearest split triangle
            double best = std::numeric_limits<double>::max();
            for (size_t u = 0; u < splitCentroid.size(); ++u) {
                const double d = normP(splitCentroid[u] - p);
                if (d < best) best = d, st = static_cast<int>(u);
            }
        }
        const ChartDistanceField::Fields x = cdf.evaluate(st, p);
        trueChart[t] = cdf.subchartBlock()[cdf.faceSubchart()[st / 2]];
        double cd = x.cdf, dc = x.dcdf;
        Point oc = x.offsets[2], ov = x.offsets[0];
        if (s.noise > 0.0) {
            cd = std::clamp(cd + s.noise * (w0(p) + white(rng) / 3.0), 0.0, 1.0);
            dc = std::clamp(dc + s.noise * (w1(p) + white(rng) / 3.0), 0.0, 1.0);
            oc = oc + Point{w2(p), w3(p)} * (0.1 * s.noise * diag);
            ov = ov + Point{w4(p), w5(p)} * (0.1 * s.noise * diag);
        }
        f.cdf.push_back(cd);
        f.dcdf.push_back(dc);
        f.offsetCentre.push_back(oc);
        f.offsetCorner.push_back(ov);
    }

    QUARRY q(R, f, QUARRY::Options());
    q.run();
    const BlockDecomposition &D = q.decomposition();
    out.seconds = q.getReport().seconds;
    out.blocks = static_cast<int>(D.blocks.size());
    out.covers = D.covers();
    if (!s.blocksObj.empty()) D.writeEdgesOBJ(s.blocksObj);

    // Area IoU of every true chart with every found one.
    const std::vector<int> &L = q.chartLabels();
    int nf = 0;
    for (const int l : L) nf = std::max(nf, l + 1);
    std::vector<double> inter(static_cast<size_t>(out.trueCharts) * (nf + 1), 0.0), aTrue(out.trueCharts, 0.0),
        aFound(nf + 1, 0.0);
    for (int t = 0; t < nt; ++t) {
        const int l = L[t] < 0 ? nf : L[t];
        inter[static_cast<size_t>(trueChart[t]) * (nf + 1) + l] += area[t];
        aTrue[trueChart[t]] += area[t];
        aFound[l] += area[t];
    }
    for (int l = 0; l < nf; ++l) out.found += aFound[l] > 0.0;
    for (int c = 0; c < out.trueCharts; ++c) {
        double best = 0.0;
        for (int l = 0; l < nf; ++l) {
            const double i = inter[static_cast<size_t>(c) * (nf + 1) + l];
            best = std::max(best, i / (aTrue[c] + aFound[l] - i));
        }
        out.matched += best > 0.5;
    }
    if (out.covers) {
        BlockQuadMesh::Options mo;
        mo.targetEdgeLength = 0.01 * diag;
        const BlockQuadMesh mesh(D, mo);
        out.meshValid = mesh.getReport().valid && mesh.getReport().invertedQuads == 0;
    }
    out.recovered = out.covers && out.meshValid && out.found == out.trueCharts && out.matched == out.trueCharts;
    std::cout << std::left << std::setw(26) << name << " charts " << std::setw(3) << out.trueCharts << " found "
              << std::setw(3) << out.found << " matched " << std::setw(3) << out.matched << " blocks "
              << std::setw(3) << out.blocks << (out.covers ? " covers" : " DOES NOT COVER")
              << (out.meshValid ? ", mesh valid" : ", mesh invalid") << "  " << (out.recovered ? "OK" : "FAIL")
              << "  (" << nt << " triangles, " << std::fixed << std::setprecision(3) << out.seconds << " s"
              << std::defaultfloat << (out.note.empty() ? "" : ", " + out.note) << ")\n";
    return out;
}

// ---------------------------------------------------------------------------
// Self-test domains
// ---------------------------------------------------------------------------

std::shared_ptr<Mesh> meshLoops(const std::vector<std::vector<Point>> &loops, const std::vector<int> &types,
                                const std::vector<Point> &seeds, double h) {
    using namespace triangle_wrapper;
    TriangleMesher2D::MeshInput input;
    for (size_t l = 0; l < loops.size(); ++l) {
        const std::vector<Point> &P = loops[l];
        const int first = static_cast<int>(input.vertlist.size());
        for (size_t k = 0; k < P.size(); ++k) {
            const Point a = P[k], b = P[(k + 1) % P.size()];
            const int n = std::max(1, static_cast<int>(std::lround(normP(b - a) / h)));
            for (int i = 0; i < n; ++i) {
                const Point p = a + (b - a) * (static_cast<double>(i) / n);
                input.vertlist.push_back({p[0], p[1]});
            }
        }
        const int last = static_cast<int>(input.vertlist.size());
        std::vector<std::array<int, 2>> segs;
        for (int i = first; i < last; ++i) segs.push_back({i, i + 1 < last ? i + 1 : first});
        input.segment_loops.push_back(segs);
        input.type.push_back(types[l]);
        input.loop_seed.push_back({seeds[l][0], seeds[l][1]});
    }
    input.h = h;
    TriangleMesher2D::Options o;
    o.min_angle_degrees = 25.0;
    const auto result = TriangleMesher2D(o).triangulate(input);
    std::vector<Point> V(result.verts.size());
    for (size_t i = 0; i < V.size(); ++i) V[i] = {result.verts[i][0], result.verts[i][1]};
    std::vector<Triangle> T(result.triangles.size());
    for (size_t i = 0; i < T.size(); ++i) T[i] = {result.triangles[i][0], result.triangles[i][1], result.triangles[i][2]};
    return std::make_shared<Mesh>(V, T);
}

std::vector<Point> circle(double cx, double cy, double r, int n) {
    std::vector<Point> p;
    for (int i = 0; i < n; ++i) p.push_back({cx + r * std::cos(2.0 * M_PI * i / n), cy + r * std::sin(2.0 * M_PI * i / n)});
    return p;
}

int selftest() {
    struct Case {
        std::string name;
        std::shared_ptr<Mesh> mesh;
    };
    const double h = 0.04;
    std::vector<Case> cases = {
        {"rectangle", meshLoops({{{0, 0}, {2, 0}, {2, 1}, {0, 1}}}, {0}, {{1, 0.5}}, h)},
        {"L", meshLoops({{{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}}}, {0}, {{0.5, 0.5}}, h)},
        {"square-hole", meshLoops({{{0, 0}, {2, 0}, {2, 2}, {0, 2}}, {{0.7, 0.7}, {0.7, 1.3}, {1.3, 1.3}, {1.3, 0.7}}},
                                  {0, 1}, {{0.2, 0.2}, {1.0, 1.0}}, h)},
        {"disk", meshLoops({circle(0, 0, 1, 64)}, {0}, {{0, 0}}, h)},
        {"plate-round-hole", meshLoops({{{0, 0}, {3, 0}, {3, 2}, {0, 2}}, circle(1.2, 1.0, 0.45, 48)}, {0, 1},
                                       {{0.1, 0.1}, {1.2, 1.0}}, h)},
    };
    int failures = 0;
    for (const double noise : {0.0, 0.15}) {
        Settings s;
        s.noise = noise;
        std::cout << (noise > 0.0 ? "-- fields with noise 0.15\n" : "-- true fields\n");
        for (const Case &c : cases) failures += roundTrip(c.mesh, s, c.name).recovered ? 0 : 1;
    }
    std::cout << (failures ? "selftest FAILED: " + std::to_string(failures) + " case(s)\n" : "selftest passed\n");
    return failures ? 1 : 0;
}

}  // namespace

int main(int argc, char **argv) {
    Settings s;
    std::vector<std::string> files;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto value = [&]() -> std::string {
            if (i + 1 >= argc) {
                std::cerr << a << " needs a value\n";
                std::exit(2);
            }
            return argv[++i];
        };
        if (a == "--selftest") return selftest();
        if (a == "--h-rel") s.hRel = std::atof(value().c_str());
        else if (a == "--noise") s.noise = std::atof(value().c_str());
        else if (a == "--seed") s.seed = static_cast<unsigned>(std::atoi(value().c_str()));
        else if (a == "--blocks-obj") s.blocksObj = value();
        else if (!a.empty() && a[0] == '-') {
            std::cerr << "unknown option " << a << "\n";
            return 2;
        } else files.push_back(a);
    }
    if (files.empty()) {
        std::cerr << "usage: " << argv[0] << " [--h-rel x] [--noise a] [--seed n] [--blocks-obj f] <mesh.obj>...\n"
                  << "       " << argv[0] << " --selftest\n";
        return 2;
    }
    int failures = 0;
    for (const std::string &f : files) {
        auto model = std::make_shared<Mesh>(f);
        const size_t slash = f.find_last_of('/');
        failures += roundTrip(model, s, slash == std::string::npos ? f : f.substr(slash + 1)).recovered ? 0 : 1;
    }
    return failures ? 1 : 0;
}
