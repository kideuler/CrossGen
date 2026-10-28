#include "SKELETON/BaseComplex.hxx"

#include <algorithm>
#include <climits>
#include <cmath>
#include <deque>
#include <numeric>
#include <tuple>
#include <unordered_map>

namespace {

long long edgeKey(int a, int b) {
    if (a > b) std::swap(a, b);
    return (static_cast<long long>(a) << 32) | static_cast<unsigned>(b);
}

// Signed, by the diagonals: exact for any simple quadrilateral.
double quadArea(const std::vector<Point> &v, const std::array<int, 4> &q) {
    return 0.5 * cross2(v[q[2]] - v[q[0]], v[q[3]] - v[q[1]]);
}

}  // namespace

BaseComplex::BaseComplex(const std::vector<Point> &vertices,
                         const std::vector<std::array<int, 4>> &quads,
                         const std::vector<std::vector<Point>> &loops,
                         const std::vector<Point> &corners)
    : BaseComplex(vertices, quads, loops, corners, Options()) {}

BaseComplex::BaseComplex(const std::vector<Point> &vertices,
                         const std::vector<std::array<int, 4>> &quads,
                         const std::vector<std::vector<Point>> &loops,
                         const std::vector<Point> &corners, const Options &opts)
    : opts_(opts), verts_(vertices), quads_(quads), loops_(loops) {
    report_.vertices = static_cast<int>(verts_.size());
    report_.quads = static_cast<int>(quads_.size());
    if (!verts_.empty()) {
        Point lo = verts_[0], hi = verts_[0];
        for (const Point &p : verts_) {
            lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
            hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
        }
        diagonal_ = std::max(normP(hi - lo), 1e-300);
    }
    for (const auto &L : loops_) {
        std::vector<double> arc(L.size() + 1, 0.0);
        for (size_t i = 0; i < L.size(); ++i)
            arc[i + 1] = arc[i] + normP(L[(i + 1) % L.size()] - L[i]);
        loopArc_.push_back(std::move(arc));
    }

    buildTopology();
    snapToBoundary();
    findSeeds(corners);
    for (int v = 0; v < static_cast<int>(verts_.size()); ++v) {
        if (seed_[v]) traceFrom(v);
    }
    for (int round = 0;; ++round) {
        updateWalls();
        if (readPatches()) break;
        if (round >= opts_.maxCutRounds || cuts_.empty()) break;
        ++report_.cutRounds;
        for (const int v : cuts_) {
            seed_[v] = 1;
            ++report_.cuts;
            traceFrom(v);
        }
    }
    assemble();
}

std::vector<Point> BaseComplex::singularPoints() const {
    std::vector<Point> out;
    for (size_t v = 0; v < verts_.size(); ++v) {
        if (singular_[v]) out.push_back(verts_[v]);
    }
    return out;
}

BaseComplex::Hit BaseComplex::locate(const Point &p) const {
    Hit h;
    for (size_t l = 0; l < loops_.size(); ++l) {
        const auto &L = loops_[l];
        for (size_t i = 0; i < L.size(); ++i) {
            const Point &a = L[i], &b = L[(i + 1) % L.size()];
            const Point ab = b - a;
            const double len2 = dotP(ab, ab);
            const double t = len2 > 0.0 ? std::clamp(dotP(p - a, ab) / len2, 0.0, 1.0) : 0.0;
            const Point q = a + ab * t;
            const double d = normP(p - q);
            if (d < h.d) {
                h.d = d;
                h.loop = static_cast<int>(l);
                h.s = loopArc_[l][i] + t * std::sqrt(len2);
                h.p = q;
            }
        }
    }
    return h;
}

void BaseComplex::snapToBoundary() {
    if (loops_.empty()) return;
    const double tol = opts_.snapTolerance * diagonal_;
    for (size_t v = 0; v < verts_.size(); ++v) {
        if (!onBoundary_[v]) continue;
        const Hit h = locate(verts_[v]);
        if (h.loop < 0 || h.d > tol || h.d == 0.0) continue;
        verts_[v] = h.p;
        ++report_.snappedVertices;
        report_.maxSnap = std::max(report_.maxSnap, h.d / diagonal_);
    }
}

void BaseComplex::buildTopology() {
    const int nv = static_cast<int>(verts_.size());
    vertEdges_.assign(nv, {});
    vertQuads_.assign(nv, {});
    quadEdges_.assign(quads_.size(), {{-1, -1, -1, -1}});
    std::unordered_map<long long, int> index;
    for (size_t q = 0; q < quads_.size(); ++q) {
        for (int k = 0; k < 4; ++k) {
            const int a = quads_[q][k], b = quads_[q][(k + 1) % 4];
            vertQuads_[a].push_back(static_cast<int>(q));
            auto it = index.find(edgeKey(a, b));
            int e;
            if (it == index.end()) {
                e = static_cast<int>(edgeVerts_.size());
                index.emplace(edgeKey(a, b), e);
                edgeVerts_.push_back({a, b});
                edgeQuads_.emplace_back();
                vertEdges_[a].push_back(e);
                vertEdges_[b].push_back(e);
            } else {
                e = it->second;
            }
            edgeQuads_[e].push_back(static_cast<int>(q));
            quadEdges_[q][k] = e;
        }
    }
    onBoundary_.assign(nv, 0);
    singular_.assign(nv, 0);
    for (size_t e = 0; e < edgeVerts_.size(); ++e) {
        if (edgeQuads_[e].size() == 2) continue;
        onBoundary_[edgeVerts_[e][0]] = onBoundary_[edgeVerts_[e][1]] = 1;
        if (edgeQuads_[e].size() > 2) {
            ++report_.nonManifoldEdges;
            singular_[edgeVerts_[e][0]] = singular_[edgeVerts_[e][1]] = 1;
        }
    }
    marked_.assign(edgeVerts_.size(), 0);
    wall_.assign(edgeVerts_.size(), 0);
}

// The valence rule of the header, and the model's corners on top of it. A
// corner the mesh has no vertex on is counted and otherwise ignored: the side
// that runs through it is the method's answer there, and regularity() scores
// it as a corner the layout rounded off.
void BaseComplex::findSeeds(const std::vector<Point> &corners) {
    const int nv = static_cast<int>(verts_.size());
    seed_.assign(nv, 0);
    macro_.assign(nv, 0);
    for (int v = 0; v < nv; ++v) {
        const int valence = static_cast<int>(vertEdges_[v].size());
        if (valence == 0) continue;
        if (onBoundary_[v] ? valence != 3 : valence != 4) singular_[v] = 1;
        if (!singular_[v]) continue;
        seed_[v] = 1;
        if (onBoundary_[v]) ++report_.boundarySingular;
        else ++report_.interiorSingular;
    }
    const double tol = opts_.tolerance * diagonal_;
    for (const Point &c : corners) {
        int best = -1;
        double bestD = tol;
        for (int v = 0; v < nv; ++v) {
            if (!onBoundary_[v]) continue;
            const double d = normP(verts_[v] - c);
            if (d <= bestD) {
                bestD = d;
                best = v;
            }
        }
        if (best < 0) {
            ++report_.missedCorners;
            continue;
        }
        ++report_.corners;
        seed_[best] = 1;
    }
}

int BaseComplex::edgeBetween(int a, int b) const {
    for (const int e : vertEdges_[a]) {
        if (edgeVerts_[e][0] == b || edgeVerts_[e][1] == b) return e;
    }
    return -1;
}

int BaseComplex::straightOn(int u, int v) const {
    if (onBoundary_[v] || vertEdges_[v].size() != 4) return -1;
    const int in = edgeBetween(u, v);
    if (in < 0) return -1;
    std::array<int, 3> skip{{in, -1, -1}};
    int n = 1;
    for (const int q : edgeQuads_[in]) {
        for (int k = 0; k < 4; ++k) {
            if (quads_[q][k] != v) continue;
            for (const int e : {quadEdges_[q][k], quadEdges_[q][(k + 3) % 4]}) {
                if (e != in && n < 3) skip[n++] = e;
            }
        }
    }
    int out = -1;
    for (const int e : vertEdges_[v]) {
        if (std::find(skip.begin(), skip.begin() + n, e) != skip.begin() + n) continue;
        if (out >= 0) return -1;
        out = e;
    }
    return out;
}

// Straight on through every regular vertex, until the line reaches dS, a
// singular vertex, a seed, or a stretch already traced -- whoever traced it
// carried it on from there.
bool BaseComplex::trace(int from, int e) {
    if (marked_[e] || edgeQuads_[e].size() != 2) return false;
    int prev = from;
    for (size_t guard = 0; guard <= edgeVerts_.size(); ++guard) {
        marked_[e] = 1;
        const int cur = edgeVerts_[e][0] == prev ? edgeVerts_[e][1] : edgeVerts_[e][0];
        if (seed_[cur] || singular_[cur] || onBoundary_[cur]) break;
        const int next = straightOn(prev, cur);
        if (next < 0 || marked_[next] || edgeQuads_[next].size() != 2) break;
        prev = cur;
        e = next;
    }
    return true;
}

void BaseComplex::traceFrom(int v) {
    for (const int e : vertEdges_[v]) {
        if (trace(v, e)) ++report_.separatrices;
    }
}

void BaseComplex::updateWalls() {
    for (size_t e = 0; e < edgeVerts_.size(); ++e)
        wall_[e] = (edgeQuads_[e].size() != 2 || marked_[e]) ? 1 : 0;
    for (size_t v = 0; v < verts_.size(); ++v) {
        int walls = 0;
        for (const int e : vertEdges_[v]) walls += wall_[e];
        macro_[v] = (seed_[v] || singular_[v] || walls >= 3) ? 1 : 0;
    }
}

// A patch is a block when its walls close into one simple counter-clockwise
// cycle that no edge appears on twice, turning at exactly four vertices --
// the ones where the patch has a single quad -- and those four are
// macrovertices and no other vertex of the cycle is. With the quads
// counter-clockwise, the walls walked with the patch on their left are that
// cycle.
bool BaseComplex::readPatches() {
    const int nq = static_cast<int>(quads_.size());
    std::vector<int> parent(nq);
    std::iota(parent.begin(), parent.end(), 0);
    auto find = [&](int a) {
        while (parent[a] != a) a = parent[a] = parent[parent[a]];
        return a;
    };
    for (size_t e = 0; e < edgeVerts_.size(); ++e) {
        if (wall_[e]) continue;
        const int a = find(edgeQuads_[e][0]), b = find(edgeQuads_[e][1]);
        if (a != b) parent[std::max(a, b)] = std::min(a, b);
    }
    patches_.clear();
    std::vector<int> patchOf(nq, -1), slot(nq, -1);
    for (int q = 0; q < nq; ++q) {
        const int r = find(q);
        if (slot[r] < 0) {
            slot[r] = static_cast<int>(patches_.size());
            patches_.emplace_back();
        }
        patchOf[q] = slot[r];
        patches_[slot[r]].quads.push_back(q);
    }

    cuts_.clear();
    bool all = true;
    for (size_t pi = 0; pi < patches_.size(); ++pi) {
        Patch &P = patches_[pi];
        std::unordered_map<int, int> count;            // the patch's quads at a vertex
        std::unordered_map<int, std::vector<int>> next;
        std::unordered_map<int, int> in;
        std::unordered_map<int, int> edgeUse;
        bool twice = false;
        for (const int q : P.quads) {
            for (int k = 0; k < 4; ++k) {
                ++count[quads_[q][k]];
                const int e = quadEdges_[q][k];
                if (!wall_[e]) continue;
                const int a = quads_[q][k], b = quads_[q][(k + 1) % 4];
                next[a].push_back(b);
                ++in[b];
                if (++edgeUse[e] > 1) twice = true;
            }
        }
        bool ok = !twice && !next.empty() && in.size() == next.size();
        for (const auto &kv : next) {
            if (kv.second.size() != 1 || in[kv.first] != 1) ok = false;
        }
        std::vector<int> cycle;
        if (ok) {
            int start = INT_MAX;
            for (const auto &kv : next) start = std::min(start, kv.first);
            int v = start;
            do {
                cycle.push_back(v);
                v = next[v][0];
            } while (v != start && cycle.size() <= next.size());
            ok = v == start && cycle.size() == next.size();
        }
        std::vector<char> corner;
        if (ok) {
            int corners = 0;
            for (const int v : cycle) {
                const int c = count[v];
                corner.push_back(c == 1 ? 1 : 0);
                if (c == 1) {
                    ++corners;
                    if (!macro_[v]) ok = false;
                } else if (c >= 3 || macro_[v]) {
                    ok = false;
                }
            }
            ok = ok && corners == 4;
        }
        if (ok) {
            const auto first = std::find(corner.begin(), corner.end(), 1) - corner.begin();
            std::rotate(cycle.begin(), cycle.begin() + first, cycle.end());
            std::rotate(corner.begin(), corner.begin() + first, corner.end());
            P.block = true;
            P.cycle = std::move(cycle);
            P.corner = std::move(corner);
            continue;
        }
        all = false;

        // The cut: the vertex of the patch's walls furthest along them from
        // any macrovertex, on dS where there is one, whose grid line runs
        // into the patch. A breadth-first search over the walls from every
        // macrovertex on them gives the distance; a loop with none on it is
        // infinitely far, which is the ring of the header.
        std::unordered_map<int, std::vector<int>> adj;
        for (const auto &kv : next) {
            for (const int b : kv.second) {
                adj[kv.first].push_back(b);
                adj[b].push_back(kv.first);
            }
        }
        std::unordered_map<int, int> dist;
        std::deque<int> queue;
        for (const auto &kv : adj) {
            if (macro_[kv.first]) {
                dist[kv.first] = 0;
                queue.push_back(kv.first);
            }
        }
        while (!queue.empty()) {
            const int v = queue.front();
            queue.pop_front();
            for (const int w : adj[v]) {
                if (dist.count(w)) continue;
                dist[w] = dist[v] + 1;
                queue.push_back(w);
            }
        }
        int best = -1;
        std::tuple<int, int, int> bestKey{-1, -1, 0};
        for (const auto &kv : adj) {
            const int v = kv.first;
            if (macro_[v] || count[v] != 2) continue;
            auto nx = next.find(v);
            if (nx == next.end() || nx->second.size() != 1 || in[v] != 1) continue;
            // The edge the patch's two quads at v share, which the cut runs
            // along into the patch.
            std::vector<int> mine;
            for (const int q : vertQuads_[v]) {
                if (patchOf[q] == static_cast<int>(pi)) mine.push_back(q);
            }
            if (mine.size() != 2) continue;
            int shared = -1;
            for (int k = 0; k < 4; ++k) {
                for (int l = 0; l < 4; ++l) {
                    const int e = quadEdges_[mine[0]][k];
                    if (e == quadEdges_[mine[1]][l] &&
                        (edgeVerts_[e][0] == v || edgeVerts_[e][1] == v))
                        shared = e;
                }
            }
            if (shared < 0 || wall_[shared]) continue;
            const auto d = dist.find(v);
            const std::tuple<int, int, int> key{onBoundary_[v] ? 1 : 0,
                                                d == dist.end() ? INT_MAX : d->second, -v};
            if (key > bestKey) {
                bestKey = key;
                best = v;
            }
        }
        if (best >= 0) cuts_.push_back(best);
    }
    return all;
}

void BaseComplex::assemble() {
    BlockDecomposition &D = decomposition_;
    D = BlockDecomposition();
    D.source = "SKELETON";
    quadInBlock_.assign(quads_.size(), 0);

    std::vector<int> macroIndex(verts_.size(), -1);
    auto vertexOf = [&](int v) {
        if (macroIndex[v] < 0) {
            macroIndex[v] = static_cast<int>(D.vertices.size());
            BlockDecomposition::MacroVertex mv;
            mv.p = verts_[v];
            mv.onBoundary = onBoundary_[v] != 0;
            D.vertices.push_back(mv);
        }
        return macroIndex[v];
    };
    std::unordered_map<int, int> macroEdgeOf;   // mesh edge -> macro edge

    double total = 0.0, covered = 0.0;
    for (const auto &q : quads_) total += std::abs(quadArea(verts_, q));
    report_.patches = static_cast<int>(patches_.size());

    for (const Patch &P : patches_) {
        if (!P.block) continue;
        const int b = static_cast<int>(D.blocks.size());
        for (const int q : P.quads) {
            quadInBlock_[q] = 1;
            covered += std::abs(quadArea(verts_, quads_[q]));
        }
        std::vector<int> at;
        for (size_t i = 0; i < P.cycle.size(); ++i) {
            if (P.corner[i]) at.push_back(static_cast<int>(i));
        }
        const int n = static_cast<int>(P.cycle.size());
        BlockDecomposition::Block block;
        for (int s = 0; s < 4; ++s) {
            std::vector<int> seq;
            const int end = at[(s + 1) % 4] + (s == 3 ? n : 0);
            for (int i = at[s]; i <= end; ++i) seq.push_back(P.cycle[i % n]);
            block.corners[s] = vertexOf(seq.front());

            const int e0 = edgeBetween(seq[0], seq[1]);
            const auto it = macroEdgeOf.find(e0);
            if (it != macroEdgeOf.end()) {
                // Its other block walked it first, the other way.
                BlockDecomposition::MacroEdge &E = D.edges[it->second];
                if (E.blockB < 0) {
                    E.blockB = b;
                    E.sideB = s;
                }
                block.edges[s] = it->second;
                block.flip[s] = true;
                continue;
            }
            BlockDecomposition::MacroEdge E;
            E.from = vertexOf(seq.front());
            E.to = vertexOf(seq.back());
            E.blockA = b;
            E.sideA = s;
            bool onDS = true;
            for (size_t i = 0; i + 1 < seq.size(); ++i) {
                if (edgeQuads_[edgeBetween(seq[i], seq[i + 1])].size() != 1) onDS = false;
            }
            E.boundary = onDS;
            if (onDS) {
                bool exact = false;
                E.points = boundarySide(seq, exact);
                ++report_.boundarySides;
                if (!exact) ++report_.chordSides;
            } else {
                for (const int v : seq) E.points.push_back(verts_[v]);
            }
            const int id = static_cast<int>(D.edges.size());
            for (size_t i = 0; i + 1 < seq.size(); ++i)
                macroEdgeOf[edgeBetween(seq[i], seq[i + 1])] = id;
            D.edges.push_back(std::move(E));
            block.edges[s] = id;
            block.flip[s] = false;
        }
        D.blocks.push_back(block);
        Grid g;
        g.block = b;
        if (gridOf(P, at, g)) grids_.push_back(std::move(g));
    }
    for (const auto &E : D.edges) {
        ++D.vertices[E.from].valence;
        ++D.vertices[E.to].valence;
    }
    report_.blocks = static_cast<int>(D.blocks.size());
    report_.badPatches = report_.patches - report_.blocks;
    report_.coverage = total > 0.0 ? covered / total : 0.0;
    countCrossings();
}

// Strictly through: each segment's two ends clear of the other's line by more
// than rounding, so a side that starts on dS, or touches it, does not count.
void BaseComplex::countCrossings() {
    report_.sidesCrossingBoundary = 0;
    for (const auto &E : decomposition_.edges) {
        if (E.boundary) continue;
        bool crosses = false;
        for (size_t i = 0; i + 1 < E.points.size() && !crosses; ++i) {
            const Point &a = E.points[i], &b = E.points[i + 1];
            const Point lo{std::min(a[0], b[0]), std::min(a[1], b[1])};
            const Point hi{std::max(a[0], b[0]), std::max(a[1], b[1])};
            for (const auto &L : loops_) {
                for (size_t j = 0; j < L.size() && !crosses; ++j) {
                    const Point &c = L[j], &d = L[(j + 1) % L.size()];
                    if (std::max(c[0], d[0]) < lo[0] || std::min(c[0], d[0]) > hi[0] ||
                        std::max(c[1], d[1]) < lo[1] || std::min(c[1], d[1]) > hi[1])
                        continue;
                    const double eps = 1e-9 * normP(b - a) * normP(d - c);
                    const double o1 = cross2(b - a, c - a), o2 = cross2(b - a, d - a);
                    const double o3 = cross2(d - c, a - c), o4 = cross2(d - c, b - c);
                    crosses = ((o1 > eps && o2 < -eps) || (o1 < -eps && o2 > eps)) &&
                              ((o3 > eps && o4 < -eps) || (o3 < -eps && o4 > eps));
                }
                if (crosses) break;
            }
        }
        if (crosses) ++report_.sidesCrossingBoundary;
    }
}

// Row 0 is side 0, corner 0 to corner 1. With the quads counter-clockwise,
// the quad above a row is the one that has the row's edge r[i] -> r[i + 1]
// as one of its own, and its other two vertices are the next row's i + 1
// and i; the rows stop at side 2.
bool BaseComplex::gridOf(const Patch &P, const std::vector<int> &at, Grid &out) const {
    const int n = static_cast<int>(P.cycle.size());
    const int ns = at[1] - at[0];
    const int nt = at[2] - at[1];
    if (ns < 1 || nt < 1 || at[3] - at[2] != ns || n - at[3] != nt) return false;
    std::unordered_map<long long, int> above;   // directed edge -> quad * 4 + k
    for (const int q : P.quads) {
        for (int k = 0; k < 4; ++k) {
            const long long key = (static_cast<long long>(quads_[q][k]) << 32) |
                                  static_cast<unsigned>(quads_[q][(k + 1) % 4]);
            above[key] = q * 4 + k;
        }
    }
    out.ns = ns;
    out.nt = nt;
    out.vert.assign(static_cast<size_t>(ns + 1) * (nt + 1), -1);
    for (int i = 0; i <= ns; ++i) out.vert[i] = P.cycle[at[0] + i];
    for (int j = 0; j < nt; ++j) {
        for (int i = 0; i < ns; ++i) {
            const int a = out.vert[j * (ns + 1) + i], b = out.vert[j * (ns + 1) + i + 1];
            const auto it = above.find((static_cast<long long>(a) << 32) | static_cast<unsigned>(b));
            if (it == above.end()) return false;
            const int q = it->second / 4, k = it->second % 4;
            const int up0 = quads_[q][(k + 3) % 4], up1 = quads_[q][(k + 2) % 4];
            int &left = out.vert[(j + 1) * (ns + 1) + i];
            int &right = out.vert[(j + 1) * (ns + 1) + i + 1];
            if ((left >= 0 && left != up0) || (right >= 0 && right != up1)) return false;
            left = up0;
            right = up1;
        }
    }
    // The last row is side 2, walked backwards.
    for (int i = 0; i <= ns; ++i) {
        if (out.vert[nt * (ns + 1) + i] != P.cycle[(at[3] - i + n) % n]) return false;
    }
    return static_cast<int>(P.quads.size()) == ns * nt;
}

// The model's own boundary between the side's two ends, walked the way the
// side runs; the side's chords if either end, or its middle, is not on dS.
std::vector<Point> BaseComplex::boundarySide(const std::vector<int> &seq, bool &exact) const {
    exact = false;
    std::vector<Point> chord;
    for (const int v : seq) chord.push_back(verts_[v]);
    if (loops_.empty() || seq.size() < 2) return chord;

    // Which way round the loop the side runs is read off a vertex in its
    // middle. A side of one mesh edge has none -- its midpoint is off dS
    // wherever the edge cuts a corner -- and takes the shorter way, the edge
    // being a single step of the mesh along dS.
    const Hit a = locate(chord.front()), z = locate(chord.back());
    const bool hasMiddle = seq.size() >= 3;
    const Hit m = hasMiddle ? locate(verts_[seq[seq.size() / 2]]) : a;
    const double tol = opts_.tolerance * diagonal_;
    if (a.loop < 0 || a.loop != z.loop || a.loop != m.loop || a.d > tol || z.d > tol || m.d > tol)
        return chord;

    const auto &L = loops_[a.loop];
    const auto &arc = loopArc_[a.loop];
    const int n = static_cast<int>(L.size());
    const double T = arc[n];
    auto forward = [&](double from, double to) {
        double d = to - from;
        if (d < 0.0) d += T;
        return d;
    };
    auto pointAt = [&](double s) {
        const int i = std::clamp(static_cast<int>(std::upper_bound(arc.begin(), arc.end(), s) -
                                                  arc.begin()) - 1, 0, n - 1);
        const double len = arc[i + 1] - arc[i];
        const double t = len > 0.0 ? (s - arc[i]) / len : 0.0;
        return L[i] + (L[(i + 1) % n] - L[i]) * t;
    };
    // Loop vertices strictly between s0 and s1, going forwards.
    auto walk = [&](double s0, double s1) {
        const double span = forward(s0, s1);
        std::vector<Point> out{pointAt(s0)};
        int i0 = static_cast<int>(std::upper_bound(arc.begin(), arc.begin() + n, s0) - arc.begin());
        if (i0 >= n) i0 = 0;
        for (int k = 0; k < n; ++k) {
            const int idx = (i0 + k) % n;
            double off = arc[idx] - s0;
            if (off < 0.0) off += T;
            if (off >= span) break;
            if (normP(L[idx] - out.back()) > 1e-9 * diagonal_) out.push_back(L[idx]);
        }
        out.push_back(pointAt(s1));
        return out;
    };
    std::vector<Point> pts;
    const bool forwards = hasMiddle ? forward(a.s, m.s) <= forward(a.s, z.s)
                                    : forward(a.s, z.s) <= forward(z.s, a.s);
    if (forwards) {
        pts = walk(a.s, z.s);
    } else {
        pts = walk(z.s, a.s);
        std::reverse(pts.begin(), pts.end());
    }
    if (pts.size() < 2) return chord;
    pts.front() = chord.front();
    pts.back() = chord.back();
    exact = true;
    return pts;
}
