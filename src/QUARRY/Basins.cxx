#include "QUARRY/Basins.hxx"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <tuple>

namespace quarry {

Basins::Basins(const Mesh &mesh, const std::vector<double> &value, const std::vector<Point> &target,
               const Options &opts)
    : mesh_(mesh), value_(value), target_(target), opts_(opts) {
    const int nt = static_cast<int>(mesh.triangles.size());
    if (static_cast<int>(value.size()) != nt || static_cast<int>(target.size()) != nt)
        throw std::invalid_argument("quarry::Basins: one value and one target per triangle");
    centroid_.resize(nt);
    area_.resize(nt);
    for (int t = 0; t < nt; ++t) {
        const Point &a = mesh.vertices[mesh.triangles[t][0]];
        const Point &b = mesh.vertices[mesh.triangles[t][1]];
        const Point &c = mesh.vertices[mesh.triangles[t][2]];
        centroid_[t] = (a + b + c) / 3.0;
        area_[t] = 0.5 * std::fabs(cross2(b - a, c - a));
    }

    sweep();
    retest();

    // Stage 2: every surviving region's core is its seed. A triangle merged in from a noise peak is not a seed
    // unless it is also within `core` of the region's own peak, since the offsets there are the decoder's least
    // reliable.
    std::vector<int> regionOf(nt, -1);
    for (int r : roots_) {
        Region g;
        g.peak = sums_[r].peak;
        g.peakValue = sums_[r].peakValue;
        g.centre = sums_[r].wq / sums_[r].w;
        g.position = centroid_[sums_[r].peak];
        g.coreArea = sums_[r].w;
        regionOf[r] = static_cast<int>(regions_.size());
        regions_.push_back(g);
    }
    labels_.assign(nt, -1);
    for (int t = 0; t < nt; ++t) {
        const int k = regionOf[find(t)];
        if (k >= 0 && value_[t] >= regions_[k].peakValue - opts_.core) labels_[t] = k;
    }
    grow(labels_);
    keepMainPieces();

    for (Region &g : regions_) g.area = 0.0, g.triangles = 0;
    for (int t = 0; t < nt; ++t) {
        if (labels_[t] < 0) {
            ++report_.unreached;
            continue;
        }
        regions_[labels_[t]].area += area_[t];
        ++regions_[labels_[t]].triangles;
    }
    report_.regions = static_cast<int>(regions_.size());
}

// The test on the targets: the two cores predict one centre, much closer together than the peaks themselves are.
bool Basins::same(const Sums &a, const Sums &b) const {
    if (a.w <= 0.0 || b.w <= 0.0) return false;
    return normP(a.wq / a.w - b.wq / b.w) < opts_.sameTarget * normP(centroid_[a.peak] - centroid_[b.peak]);
}

int Basins::find(int t) {
    while (parent_[t] != t) {
        parent_[t] = parent_[parent_[t]];
        t = parent_[t];
    }
    return t;
}

// Stage 1, the persistence sweep. Triangles are visited from the highest value down, each attached to the region
// of its highest visited edge-neighbour (gradient ascent on the edge graph), so that a region is always connected
// and its root is its peak. Where a triangle touches a second region, the two are compared at that triangle's
// value, the saddle between them, and merged on either test of the header. A pair that passes neither test may
// meet again lower down and is compared again there: the persistence grows as the saddle falls, and the cores
// stop growing once the sweep is a `core` below both peaks, so the second look decides nothing new unless the
// cores grew in between.
void Basins::sweep() {
    const int nt = static_cast<int>(mesh_.triangles.size());
    std::vector<int> order(nt);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return value_[a] > value_[b]; });
    parent_.assign(nt, -1);
    sums_.assign(nt, Sums());

    auto addCore = [&](int r, int t) {
        Sums &s = sums_[r];
        if (value_[t] < s.peakValue - opts_.core) return;
        s.w += area_[t];
        s.wq = s.wq + target_[t] * area_[t];
    };

    for (const int f : order) {
        int best = -1;
        for (int k = 0; k < 3; ++k) {
            const int g = mesh_.triangleAdjacency[f][k];
            if (g >= 0 && parent_[g] >= 0 && (best < 0 || value_[g] > value_[best])) best = g;
        }
        if (best < 0) {
            parent_[f] = f;
            sums_[f].peak = f;
            sums_[f].peakValue = value_[f];
            addCore(f, f);
            ++report_.peaks;
            continue;
        }
        const int r = find(best);
        parent_[f] = r;
        addCore(r, f);
        for (int k = 0; k < 3; ++k) {
            const int g = mesh_.triangleAdjacency[f][k];
            if (g < 0 || parent_[g] < 0) continue;
            const int a = find(f), b = find(g);
            if (a == b) continue;
            const bool aLower = sums_[a].peakValue < sums_[b].peakValue ||
                                (sums_[a].peakValue == sums_[b].peakValue && a > b);
            const int lo = aLower ? a : b, hi = aLower ? b : a;
            const Sums &L = sums_[lo], &H = sums_[hi];
            bool merge = false;
            if (same(L, H)) {
                merge = true;
                ++report_.targetMerges;
            } else if (L.peakValue - value_[f] < opts_.noise) {
                merge = true;
                ++report_.noiseMerges;
            }
            if (!merge) continue;
            parent_[lo] = hi;
            // The lower core joins the higher one's sums only if it lies in the higher one's core band too: a
            // noise peak's few triangles near a wall would otherwise pull the centre towards the blur.
            if (L.peakValue >= H.peakValue - opts_.core) {
                Sums &h = sums_[hi];
                h.w += L.w;
                h.wq = h.wq + L.wq;
            }
        }
    }
    for (int t = 0; t < nt; ++t) find(t);
}

// After the sweep, every pair of adjacent regions is tested again on the target alone, now that every core is
// complete: a chart the decoder cut in two by a false wall meets its other half at the bottom of that wall, where
// the halves' cores may still have been growing. Regions whose peak never reached minPeak are then dropped, and
// their triangles are grown into by their neighbours in Stage 2.
void Basins::retest() {
    const int nt = static_cast<int>(mesh_.triangles.size());
    for (bool changed = true; changed;) {
        changed = false;
        std::vector<std::pair<int, int>> pairs;
        for (int t = 0; t < nt; ++t)
            for (int k = 0; k < 3; ++k) {
                const int u = mesh_.triangleAdjacency[t][k];
                if (u <= t) continue;
                const int a = find(t), b = find(u);
                if (a != b) pairs.emplace_back(std::min(a, b), std::max(a, b));
            }
        std::sort(pairs.begin(), pairs.end());
        pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
        for (auto [a0, b0] : pairs) {
            const int a = find(a0), b = find(b0);
            if (a == b) continue;
            const Sums &A = sums_[a], &B = sums_[b];
            if (!same(A, B)) continue;
            const bool aLower = A.peakValue < B.peakValue || (A.peakValue == B.peakValue && a > b);
            const int lo = aLower ? a : b, hi = aLower ? b : a;
            parent_[lo] = hi;
            Sums &h = sums_[hi];
            h.w += sums_[lo].w;
            h.wq = h.wq + sums_[lo].wq;
            ++report_.targetMerges;
            changed = true;
        }
        for (int t = 0; t < nt; ++t) find(t);
    }
    for (int t = 0; t < nt; ++t) {
        if (parent_[t] != t) continue;
        if (sums_[t].peakValue < opts_.minPeak || sums_[t].w <= 0.0) {
            ++report_.dissolved;
            continue;
        }
        roots_.push_back(t);
    }
}

// Stage 2's growth: the release's priority growth (face_clustering_via_seeds), every unlabelled triangle taken in
// order of the distance from the centre it predicts to the centre of the region reaching it, and claimed by the
// first region to reach it. A triangle on a blurred wall predicts a point between two centres and goes to the
// nearer, which puts the boundary where the prediction changes hands.
void Basins::grow(std::vector<int> &labels) const {
    using Item = std::tuple<double, int, int>;
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> queue;
    const int nt = static_cast<int>(mesh_.triangles.size());
    auto key = [&](int t, int k) { return normP(target_[t] - regions_[k].centre); };
    for (int t = 0; t < nt; ++t) {
        if (labels[t] < 0) continue;
        for (int e = 0; e < 3; ++e) {
            const int u = mesh_.triangleAdjacency[t][e];
            if (u >= 0 && labels[u] < 0) queue.emplace(key(u, labels[t]), u, labels[t]);
        }
    }
    while (!queue.empty()) {
        const auto [d, t, k] = queue.top();
        queue.pop();
        if (labels[t] >= 0) continue;
        labels[t] = k;
        for (int e = 0; e < 3; ++e) {
            const int u = mesh_.triangleAdjacency[t][e];
            if (u >= 0 && labels[u] < 0) queue.emplace(key(u, k), u, k);
        }
    }
}

// A region the growth left in pieces -- two cores of one merged region that other regions grew between, or a
// region's triangles cut off by a neighbour -- keeps the piece with the most core in it and gives up the rest,
// which the neighbours then grow into as before. A chart is a disk, and a block built from a region in two pieces
// would be no block at all.
void Basins::keepMainPieces() {
    const int nt = static_cast<int>(mesh_.triangles.size());
    const int K = static_cast<int>(regions_.size());
    std::vector<int> piece(nt, -1);
    std::vector<double> pieceCore;
    std::vector<int> pieceRegion;
    std::vector<int> stack;
    for (int t = 0; t < nt; ++t) {
        if (labels_[t] < 0 || piece[t] >= 0) continue;
        const int id = static_cast<int>(pieceCore.size());
        const int k = labels_[t];
        pieceCore.push_back(0.0);
        pieceRegion.push_back(k);
        piece[t] = id;
        stack.assign(1, t);
        while (!stack.empty()) {
            const int s = stack.back();
            stack.pop_back();
            // Core area, with every triangle's own area as a tie-break so that a piece with no core still has a
            // size to be compared by.
            pieceCore[id] += area_[s] * (value_[s] >= regions_[k].peakValue - opts_.core ? 1.0 : 1e-6);
            for (int e = 0; e < 3; ++e) {
                const int u = mesh_.triangleAdjacency[s][e];
                if (u >= 0 && piece[u] < 0 && labels_[u] == k) {
                    piece[u] = id;
                    stack.push_back(u);
                }
            }
        }
    }
    std::vector<int> mainPiece(K, -1);
    for (int p = 0; p < static_cast<int>(pieceCore.size()); ++p) {
        const int k = pieceRegion[p];
        if (mainPiece[k] < 0 || pieceCore[p] > pieceCore[mainPiece[k]]) mainPiece[k] = p;
    }
    bool cut = false;
    for (int t = 0; t < nt; ++t) {
        if (labels_[t] >= 0 && piece[t] != mainPiece[labels_[t]]) {
            labels_[t] = -1;
            cut = true;
        }
    }
    report_.splitPieces = static_cast<int>(pieceCore.size()) - K;
    if (cut) grow(labels_);
}

}  // namespace quarry
