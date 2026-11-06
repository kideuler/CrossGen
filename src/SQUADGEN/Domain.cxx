#include "SQUADGEN/Domain.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace squadgen {
namespace {

double signedArea(const Point &a, const Point &b, const Point &c) {
    return 0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]));
}

double pointSegment(const Point &p, const Point &a, const Point &b) {
    const Point e = b - a;
    const double t = std::clamp(dotP(p - a, e) / std::max(dotP(e, e), 1e-30), 0.0, 1.0);
    const Point d = p - a - e * t;
    return std::sqrt(d[0] * d[0] + d[1] * d[1]);
}

// Every edge as (lower, higher) vertex id, in the order squadgen2d lists them -- all first edges (0-1), then all
// second (1-2), then all third (2-0) -- and the rank of each among the distinct ones, which is np.unique's order.
void uniqueEdges(const std::vector<Triangle> &T, std::vector<std::pair<int, int>> &distinct, std::vector<int> &rank) {
    const size_t F = T.size();
    std::vector<std::pair<int, int>> all(3 * F);
    for (int k = 0; k < 3; ++k)
        for (size_t t = 0; t < F; ++t) {
            const int a = T[t][k], b = T[t][(k + 1) % 3];
            all[k * F + t] = {std::min(a, b), std::max(a, b)};
        }
    std::vector<int> order(all.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = static_cast<int>(i);
    std::sort(order.begin(), order.end(), [&](int i, int j) { return all[i] < all[j]; });
    rank.assign(all.size(), -1);
    distinct.clear();
    for (const int i : order) {
        if (distinct.empty() || distinct.back() != all[i]) distinct.push_back(all[i]);
        rank[i] = static_cast<int>(distinct.size()) - 1;
    }
}

// The median of the lengths of every triangle's first edge, as np.median takes it: the mean of the two middle values
// of an even count.
double medianFirstEdge(const Model &m) {
    std::vector<double> e(m.triangles.size());
    for (size_t t = 0; t < e.size(); ++t) {
        const Point d = m.vertices[m.triangles[t][1]] - m.vertices[m.triangles[t][0]];
        e[t] = std::sqrt(d[0] * d[0] + d[1] * d[1]);
    }
    if (e.empty()) return 0.0;
    std::sort(e.begin(), e.end());
    const size_t n = e.size();
    return n % 2 ? e[n / 2] : 0.5 * (e[n / 2 - 1] + e[n / 2]);
}

// One 1-to-4 midpoint split (roundtrip.subdivide): every vertex and the orientation of every triangle kept.
Model subdivide(const Model &m) {
    std::vector<std::pair<int, int>> edges;
    std::vector<int> rank;
    uniqueEdges(m.triangles, edges, rank);
    Model out;
    out.vertices = m.vertices;
    out.vertices.reserve(m.vertices.size() + edges.size());
    for (const auto &[a, b] : edges) out.vertices.push_back((m.vertices[a] + m.vertices[b]) * 0.5);
    const size_t F = m.triangles.size();
    const int base = static_cast<int>(m.vertices.size());
    out.triangles.resize(4 * F);
    for (size_t t = 0; t < F; ++t) {
        const int a = m.triangles[t][0], b = m.triangles[t][1], c = m.triangles[t][2];
        const int ab = base + rank[t], bc = base + rank[F + t], ca = base + rank[2 * F + t];
        out.triangles[t] = {a, ab, ca};
        out.triangles[F + t] = {ab, b, bc};
        out.triangles[2 * F + t] = {ca, bc, c};
        out.triangles[3 * F + t] = {ab, bc, ca};
    }
    return out;
}

}  // namespace

// ── Random ───────────────────────────────────────────────────────────────────

double Random::uniform() { return static_cast<double>(engine_() >> 11) * (1.0 / 9007199254740992.0); }

double Random::normal() {
    if (haveSpare_) {
        haveSpare_ = false;
        return spare_;
    }
    const double u1 = 1.0 - uniform();   // (0, 1]: the logarithm stays finite
    const double u2 = uniform();
    const double r = std::sqrt(-2.0 * std::log(u1));
    spare_ = r * std::sin(2.0 * M_PI * u2);
    haveSpare_ = true;
    return r * std::cos(2.0 * M_PI * u2);
}

size_t Random::below(size_t n) {
    if (n == 0) return 0;
    return std::min(n - 1, static_cast<size_t>(uniform() * static_cast<double>(n)));
}

// ── the model ────────────────────────────────────────────────────────────────

Model readOBJ(const std::string &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open the OBJ '" + path + "'");
    Model m;
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("v ", 0) == 0) {
            std::istringstream ss(line.substr(2));
            double x = 0.0, y = 0.0;
            ss >> x >> y;
            m.vertices.push_back({x, y});
        } else if (line.rfind("f ", 0) == 0) {
            std::istringstream ss(line.substr(2));
            Triangle t{};
            std::string tok;
            int k = 0;
            while (k < 3 && ss >> tok) t[k++] = std::stoi(tok.substr(0, tok.find('/'))) - 1;
            if (k == 3) m.triangles.push_back(t);
        }
    }
    const int n = static_cast<int>(m.vertices.size());
    for (Triangle &t : m.triangles) {
        for (const int v : t)
            if (v < 0 || v >= n) throw std::runtime_error("'" + path + "' has a face on a vertex it does not list");
        if (signedArea(m.vertices[t[0]], m.vertices[t[1]], m.vertices[t[2]]) < 0.0) std::swap(t[0], t[2]);
    }
    if (m.triangles.empty()) throw std::runtime_error("'" + path + "' has no triangles");
    return m;
}

Model modelOf(const Mesh &mesh) { return Model{mesh.vertices, mesh.triangles}; }

// ── the domain in the frame ──────────────────────────────────────────────────

Domain::Domain(const Model &m) : triangles_(m.triangles) {
    if (m.triangles.empty()) throw std::runtime_error("SQuadGen: the model has no triangles");
    // dS: the edges one triangle uses.
    std::vector<std::pair<int, int>> edges;
    std::vector<int> rank;
    uniqueEdges(m.triangles, edges, rank);
    std::vector<int> uses(edges.size(), 0);
    for (const int r : rank) ++uses[r];
    std::vector<char> onBoundary(m.vertices.size(), 0);
    std::vector<std::pair<int, int>> boundary;
    for (size_t e = 0; e < edges.size(); ++e)
        if (uses[e] == 1) {
            boundary.push_back(edges[e]);
            onBoundary[edges[e].first] = onBoundary[edges[e].second] = 1;
        }

    // The area centroid, and 0.999 over the farthest boundary vertex from it.
    double area = 0.0;
    Point c{0.0, 0.0};
    for (const Triangle &t : m.triangles) {
        const Point &a = m.vertices[t[0]], &b = m.vertices[t[1]], &d = m.vertices[t[2]];
        const double w = std::fabs(signedArea(a, b, d));
        c = c + (a + b + d) / 3.0 * w;
        area += w;
    }
    if (!(area > 0.0)) throw std::runtime_error("SQuadGen: the model has no area");
    frame_.centre = c / area;
    double far = 0.0;
    for (size_t v = 0; v < m.vertices.size(); ++v)
        if (onBoundary[v]) far = std::max(far, normP(m.vertices[v] - frame_.centre));
    if (!(far > 0.0)) throw std::runtime_error("SQuadGen: the model has no boundary");
    frame_.scale = 0.999 / far;

    vertices_.resize(m.vertices.size());
    for (size_t v = 0; v < m.vertices.size(); ++v) vertices_[v] = frame_.toFrame(m.vertices[v]);
    cdf_.resize(triangles_.size());
    double sum = 0.0;
    for (size_t t = 0; t < triangles_.size(); ++t) {
        const Triangle &T = triangles_[t];
        sum += std::fabs(signedArea(vertices_[T[0]], vertices_[T[1]], vertices_[T[2]]));
        cdf_[t] = sum;
    }
    for (const auto &[a, b] : boundary) segments_.push_back({vertices_[a], vertices_[b]});
}

std::vector<Point> Domain::sampleArea(size_t n, Random &rng) const {
    std::vector<Point> out(n);
    for (size_t i = 0; i < n; ++i) {
        const double x = rng.uniform() * cdf_.back();
        const size_t t = std::min(static_cast<size_t>(std::lower_bound(cdf_.begin(), cdf_.end(), x) - cdf_.begin()),
                                  cdf_.size() - 1);
        double u = rng.uniform(), v = rng.uniform();
        if (u + v > 1.0) u = 1.0 - u, v = 1.0 - v;
        const Triangle &T = triangles_[t];
        out[i] = vertices_[T[0]] * (1.0 - u - v) + vertices_[T[1]] * u + vertices_[T[2]] * v;
    }
    return out;
}

std::vector<Point> Domain::tokenSet(int n, Random &rng) const {
    const std::vector<Point> candidates = sampleArea(std::max<size_t>(8 * static_cast<size_t>(n), 8192), rng);
    std::vector<Point> out;
    out.reserve(n);
    for (const int i : farthestPoints(candidates, n, rng)) out.push_back(candidates[i]);
    return out;
}

std::vector<float> Domain::insideSDF(const std::vector<Point> &points, double tau) const {
    std::vector<float> out(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        double d = tau;
        for (const auto &s : segments_) d = std::min(d, pointSegment(points[i], s[0], s[1]));
        out[i] = static_cast<float>(std::max(-1.0, -d / tau));
    }
    return out;
}

// ── farthest points ──────────────────────────────────────────────────────────

std::vector<int> farthestPoints(const std::vector<Point> &candidates, int n, Random &rng) {
    const size_t M = candidates.size();
    if (n <= 0 || M == 0) return {};
    if (static_cast<size_t>(n) > M) throw std::runtime_error("farthestPoints: more points asked for than candidates");
    // Single precision and separate coordinate arrays, which the update loop vectorises over: 8192 picks of 65 536
    // candidates (N = 2048's key/value set) are half a billion distance updates.
    std::vector<float> X(M), Y(M), D(M, std::numeric_limits<float>::infinity());
    for (size_t i = 0; i < M; ++i) X[i] = static_cast<float>(candidates[i][0]), Y[i] = static_cast<float>(candidates[i][1]);
    std::vector<int> chosen(n);
    size_t j = rng.below(M);
    chosen[0] = static_cast<int>(j);
    for (int k = 1; k < n; ++k) {
        const float px = X[j], py = Y[j];
        for (size_t i = 0; i < M; ++i) {
            const float dx = X[i] - px, dy = Y[i] - py;
            D[i] = std::min(D[i], dx * dx + dy * dy);
        }
        j = static_cast<size_t>(std::max_element(D.begin(), D.end()) - D.begin());
        chosen[k] = static_cast<int>(j);
    }
    return chosen;
}

// ── the extraction mesh ──────────────────────────────────────────────────────

Model extractionMesh(const Model &model, double hRel, int *levels) {
    // Without unused vertices, in the order of their ids (np.unique(T, return_inverse=True)).
    std::vector<int> map(model.vertices.size(), -1);
    for (const Triangle &t : model.triangles)
        for (const int v : t) map[v] = 0;
    Model m;
    for (size_t v = 0; v < map.size(); ++v)
        if (map[v] == 0) {
            map[v] = static_cast<int>(m.vertices.size());
            m.vertices.push_back(model.vertices[v]);
        }
    m.triangles.reserve(model.triangles.size());
    for (const Triangle &t : model.triangles) m.triangles.push_back({map[t[0]], map[t[1]], map[t[2]]});

    Point lo = m.vertices.front(), hi = lo;
    for (const Point &p : m.vertices) {
        lo = {std::min(lo[0], p[0]), std::min(lo[1], p[1])};
        hi = {std::max(hi[0], p[0]), std::max(hi[1], p[1])};
    }
    const Point ext = hi - lo;
    const double diag = std::sqrt(ext[0] * ext[0] + ext[1] * ext[1]);
    int level = 0;
    while (level < 4 && medianFirstEdge(m) > 1.5 * hRel * diag) {
        m = subdivide(m);
        ++level;
    }
    if (levels) *levels = level;
    return m;
}

std::vector<float> centroidsInFrame(const Model &m, const Frame &frame) {
    std::vector<float> out(2 * m.triangles.size());
    for (size_t t = 0; t < m.triangles.size(); ++t) {
        const Triangle &T = m.triangles[t];
        const Point c = frame.toFrame((m.vertices[T[0]] + m.vertices[T[1]] + m.vertices[T[2]]) / 3.0);
        out[2 * t] = static_cast<float>(c[0]);
        out[2 * t + 1] = static_cast<float>(c[1]);
    }
    return out;
}

std::vector<float> flatten(const std::vector<Point> &points) {
    std::vector<float> out(2 * points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        out[2 * i] = static_cast<float>(points[i][0]);
        out[2 * i + 1] = static_cast<float>(points[i][1]);
    }
    return out;
}

}  // namespace squadgen
