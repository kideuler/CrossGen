#include "mesh/BoundaryFeatures.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <unordered_map>

namespace {

// One edge of a region's boundary, directed with the region on its left.
struct HalfEdge {
    int a = -1, b = -1;
    bool onBoundary = false;   // dS rather than a material interface
    bool used = false;
};

// Swept counter-clockwise from `out` to `back`, in (0, 2 pi]: the interior
// angle at a loop vertex, `out` leaving it and `back` pointing to where the
// walk came from.
double interiorAngle(const Point &out, const Point &back) {
    double a = std::atan2(cross2(out, back), dotP(out, back));
    if (a <= 0.0) a += 2.0 * M_PI;
    return a;
}

// The header's minimum for one region: sum_c |a_c - n_c| + |4 chi - sum_c
// (2 - n_c)| over n_c >= 1, a_c = 2 theta_c / pi. Each corner starts at the
// cheaper-to-reach lower neighbour of a_c and may move up one, which changes
// the corner's own cost by delta_c and the interior term's argument by one;
// for m corners moved, the m smallest deltas are the ones to move, so trying
// every m over the sorted deltas finds the minimum.
double regionMinimumDefect(const std::vector<double> &a, int euler) {
    double base = 0.0;
    double sum = 0.0;
    std::vector<double> delta;
    for (double x : a) {
        const double lo = std::max(1.0, std::floor(x));
        const double hi = std::max(1.0, std::ceil(x));
        base += std::fabs(x - lo);
        sum += lo;
        if (hi != lo) delta.push_back(std::fabs(x - hi) - std::fabs(x - lo));
    }
    std::sort(delta.begin(), delta.end());
    const double K = 4.0 * euler - 2.0 * static_cast<double>(a.size());
    double best = base + std::fabs(K + sum);
    double moved = base;
    for (size_t m = 0; m < delta.size(); ++m) {
        moved += delta[m];
        best = std::min(best, moved + std::fabs(K + sum + static_cast<double>(m + 1)));
    }
    return best;
}

}  // namespace

BoundaryFeatures::BoundaryFeatures(const Mesh &mesh) : BoundaryFeatures(mesh, Options()) {}

BoundaryFeatures::BoundaryFeatures(const Mesh &mesh, const Options &options) : options_(options) {
    const double cornerTurn = options.cornerAngle * M_PI / 180.0;
    const double curvedTurn = options.curvedAngle * M_PI / 180.0;
    const int nT = static_cast<int>(mesh.triangles.size());
    if (nT == 0) return;

    // Every triangle edge that bounds its region, directed with the triangle
    // on its left. The loader and Mesh's constructors leave triangles counter-
    // clockwise, but a mesh built some other way need not be, and a boundary
    // walked the wrong way round would turn every convex corner reflex.
    auto regionOf = [&](int t) {
        return (t < static_cast<int>(mesh.triangleComponent.size())) ? mesh.triangleComponent[t] : 0;
    };
    std::map<int, std::vector<HalfEdge>> byRegion;
    double area = 0.0;
    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = mesh.triangles[t];
        const double twice = cross2(mesh.vertices[tri[1]] - mesh.vertices[tri[0]],
                                    mesh.vertices[tri[2]] - mesh.vertices[tri[0]]);
        area += 0.5 * std::fabs(twice);
        for (int k = 0; k < 3; ++k) {
            const int nb = mesh.triangleAdjacency[t][k];
            if (nb >= 0 && regionOf(nb) == regionOf(t)) continue;
            HalfEdge h;
            h.a = tri[k];
            h.b = tri[(k + 1) % 3];
            if (twice < 0.0) std::swap(h.a, h.b);
            h.onBoundary = (nb < 0);
            byRegion[regionOf(t)].push_back(h);
        }
    }

    // Chain each region's half-edges into closed loops. Where a region touches
    // itself at a vertex, two ways out leave that vertex; the one that keeps
    // the walk on the wedge it arrived by is the one with the smallest
    // interior angle.
    for (auto &kv : byRegion) {
        const int region = kv.first;
        std::vector<HalfEdge> &H = kv.second;
        std::unordered_map<int, std::vector<int>> leaving;
        for (int i = 0; i < static_cast<int>(H.size()); ++i) leaving[H[i].a].push_back(i);

        for (int s = 0; s < static_cast<int>(H.size()); ++s) {
            if (H[s].used) continue;
            Loop L;
            L.region = region;
            int cur = s;
            while (true) {
                H[cur].used = true;
                L.vertices.push_back(H[cur].a);
                L.onBoundary.push_back(H[cur].onBoundary);
                const int v = H[cur].b;
                const Point back = mesh.vertices[H[cur].a] - mesh.vertices[v];
                int next = -1;
                double best = std::numeric_limits<double>::infinity();
                for (int j : leaving[v]) {
                    if (H[j].used && j != s) continue;
                    const double theta = interiorAngle(mesh.vertices[H[j].b] - mesh.vertices[v], back);
                    if (theta < best) {
                        best = theta;
                        next = j;
                    }
                }
                if (next < 0 || next == s) break;   // closed (or an open chain, which a mesh boundary is not)
                cur = next;
            }
            loops_.push_back(std::move(L));
        }
    }

    // Turning, corners and runs of each loop.
    std::map<int, int> loopsOfRegion;
    std::map<int, std::vector<double>> quarterTurns;   // per region: 2 theta / pi of each corner
    double perimeter = 0.0, curved = 0.0, interfaceLength = 0.0;
    double shortestRun = std::numeric_limits<double>::infinity();
    for (int l = 0; l < static_cast<int>(loops_.size()); ++l) {
        Loop &L = loops_[l];
        const int n = static_cast<int>(L.vertices.size());
        ++loopsOfRegion[L.region];
        L.turn.assign(n, 0.0);
        std::vector<double> edge(n, 0.0);
        for (int i = 0; i < n; ++i) {
            const Point &u = mesh.vertices[L.vertices[(i + n - 1) % n]];
            const Point &v = mesh.vertices[L.vertices[i]];
            const Point &w = mesh.vertices[L.vertices[(i + 1) % n]];
            L.turn[i] = std::atan2(cross2(v - u, w - v), dotP(v - u, w - v));
            edge[i] = normP(w - v);
            L.length += edge[i];
            L.signedArea += 0.5 * cross2(v, w);
        }
        std::vector<bool> isCorner(n, false);
        for (int i = 0; i < n; ++i) {
            if (!(std::fabs(L.turn[i]) > cornerTurn)) continue;
            isCorner[i] = true;
            Corner c;
            c.vertex = L.vertices[i];
            c.region = L.region;
            c.loop = l;
            c.p = mesh.vertices[c.vertex];
            c.angle = M_PI - L.turn[i];
            c.blocks = std::max(1, static_cast<int>(std::lround(2.0 * c.angle / M_PI)));
            c.step = std::min(edge[i], edge[(i + n - 1) % n]);
            L.corners.push_back(static_cast<int>(corners_.size()));
            quarterTurns[L.region].push_back(2.0 * c.angle / M_PI);
            corners_.push_back(c);
        }
        // dS's length and how much of it bends; an interface edge is walked
        // once from each of its two regions.
        for (int i = 0; i < n; ++i) {
            if (!L.onBoundary[i]) {
                interfaceLength += 0.5 * edge[i];
                continue;
            }
            perimeter += edge[i];
            const int j = (i + 1) % n;
            const double bend = std::max(isCorner[i] ? 0.0 : std::fabs(L.turn[i]),
                                         isCorner[j] ? 0.0 : std::fabs(L.turn[j]));
            if (bend > curvedTurn) curved += edge[i];
        }
        // Runs: corner to corner, or the whole loop when it has at most one.
        std::vector<int> at;
        for (int i = 0; i < n; ++i)
            if (isCorner[i]) at.push_back(i);
        if (at.size() <= 1) {
            Run r;
            r.loop = l;
            r.length = L.length;
            for (int i = 0; i < n; ++i)
                if (!isCorner[i]) r.turning += L.turn[i];
            runs_.push_back(r);
        } else {
            for (size_t q = 0; q < at.size(); ++q) {
                const int from = at[q], to = at[(q + 1) % at.size()];
                Run r;
                r.loop = l;
                for (int i = from; i != to; i = (i + 1) % n) {
                    r.length += edge[i];
                    if (i != from) r.turning += L.turn[i];
                }
                runs_.push_back(r);
            }
        }
        for (const Run &r : runs_)
            if (r.loop == l) shortestRun = std::min(shortestRun, r.length);
    }

    Summary &S = summary_;
    S.area = area;
    S.perimeter = perimeter;
    S.regions = static_cast<int>(loopsOfRegion.size());
    for (const auto &kv : loopsOfRegion) {
        const int chi = 2 - kv.second;
        S.holes += kv.second - 1;
        S.euler += chi;
        const std::vector<double> &a = quarterTurns[kv.first];
        double ideal = 0.0;
        for (double x : a) ideal += 2.0 - std::max(1.0, std::round(x));
        S.singularityBound += std::fabs(4.0 * chi - ideal);
        S.minimumDefect += regionMinimumDefect(a, chi);
    }
    for (const Corner &c : corners_) {
        const double deg = c.angle * 180.0 / M_PI;
        ++S.corners;
        if (c.blocks == 1) ++S.cornersOneBlock;
        else if (c.blocks == 2) ++S.cornersTwoBlocks;
        else if (c.blocks == 3) ++S.cornersThreeBlocks;
        else ++S.cornersFourBlocks;
        if (deg < 80.0) ++S.acuteCorners;
        if (std::fabs(deg - 135.0) < 10.0 || std::fabs(deg - 225.0) < 10.0) ++S.ambiguousCorners;
        S.cornerDefect += std::fabs(2.0 * c.angle / M_PI - c.blocks);
    }
    const double rootA = std::sqrt(std::max(area, 0.0));
    if (area > 0.0) {
        S.isoperimetricRatio = perimeter * perimeter / (4.0 * M_PI * area);
        S.shortestRun = std::isfinite(shortestRun) ? shortestRun / rootA : 0.0;
        S.interfaceLength = interfaceLength / rootA;
    }
    if (perimeter > 0.0) S.curvedFraction = curved / perimeter;
}
