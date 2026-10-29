#include "mesh/ChartDistanceField.hxx"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <unordered_map>

#include "mesh/BlockQuadMesh.hxx"
#include "mesh/NpzWriter.hxx"

namespace {

double signedArea(const Point &a, const Point &b, const Point &c) { return 0.5 * cross2(b - a, c - a); }

std::uint64_t directedKey(int tail, int head) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(tail)) << 32) |
           static_cast<std::uint32_t>(head);
}

// The split meshes' half-edges: 4f + k runs from corner k of face f to corner
// k + 1, and its pair is the half-edge the other way, or -1 on dS. Built once
// by buildCoordinates() and read by it only.
struct HalfEdges {
    std::vector<int> pair;
    std::vector<int> edge;      // undirected edge id of each half-edge
    std::unordered_map<std::uint64_t, int> byEnds;
    int edges = 0;

    static int opposite(int h) { return (h & ~3) | ((h + 2) & 3); }

    int find(int tail, int head) const {
        const auto it = byEnds.find(directedKey(tail, head));
        return it == byEnds.end() ? -1 : it->second;
    }
};

}  // namespace

ChartDistanceField::ChartDistanceField(const BlockDecomposition &D, const Options &opts) : opts_(opts) {
    const auto start = std::chrono::steady_clock::now();
    if (D.blocks.empty() || !D.covers())
        throw std::invalid_argument("ChartDistanceField: the blocks do not cover the model "
                                    "(BlockDecomposition::covers() is false)");

    // The target edge length, relative to the decomposition's own extent
    // unless given outright.
    double lo[2] = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    double hi[2] = {-lo[0], -lo[1]};
    for (const BlockDecomposition::MacroEdge &e : D.edges) {
        for (const Point &p : e.points) {
            for (int d = 0; d < 2; ++d) {
                lo[d] = std::min(lo[d], p[d]);
                hi[d] = std::max(hi[d], p[d]);
            }
        }
    }
    const double diagonal = std::hypot(hi[0] - lo[0], hi[1] - lo[1]);
    const double h = opts_.targetEdgeLength > 0.0 ? opts_.targetEdgeLength : opts_.relativeEdgeLength * diagonal;
    if (!(h > 0.0) || !std::isfinite(h))
        throw std::invalid_argument("ChartDistanceField: no positive target edge length");

    BlockQuadMesh::Options mo;
    mo.targetEdgeLength = h;
    mo.minIntervals = std::max(1, opts_.minIntervals);
    mo.smoothingPasses = opts_.smoothingPasses;
    const BlockQuadMesh grid(D, mo);
    const BlockQuadMesh::Report &gr = grid.getReport();
    if (gr.unmeshedBlocks > 0 || grid.blocks().size() != D.blocks.size())
        throw std::invalid_argument("ChartDistanceField: BlockQuadMesh left " + std::to_string(gr.unmeshedBlocks) +
                                    " of " + std::to_string(D.blocks.size()) + " blocks unmeshed");
    if (!gr.conforming)
        throw std::invalid_argument("ChartDistanceField: the block mesh is not conforming");

    report_.target = h;
    report_.blocks = static_cast<int>(D.blocks.size());
    report_.chords = static_cast<int>(grid.chords().size());
    report_.quads = static_cast<int>(grid.quads().size());

    buildSplitMesh(grid, D);
    buildCoordinates();
    sample();

    report_.valid = report_.splitMismatches == 0 && report_.coordinateMismatches == 0 && report_.missingPairs == 0;
    if (report_.invertedFaces > 0)
        report_.messages.push_back(std::to_string(report_.invertedFaces) +
                                   " split-mesh quads are folded: the block mesh inherited a fold "
                                   "the Winslow pass did not undo");
    report_.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

// ---------------------------------------------------------------------------
// buildSplitMesh()
//
// Step 1 of the header: one split per chord, cut through every block of the
// chord, giving each block's grid one extra column and one extra row. The
// block mesh's vertices come first and keep their indices; a vertex on a cut
// grid edge is keyed by that edge, so the two blocks either side of a macro
// edge get the same split vertex -- and get it at the same place, which is
// counted (splitMismatches) rather than assumed, since it holds only if the
// chord's orientation has been carried correctly from block to block.
// ---------------------------------------------------------------------------
void ChartDistanceField::buildSplitMesh(const BlockQuadMesh &grid, const BlockDecomposition &D) {
    const std::vector<Point> &Q = grid.vertices();
    quadVertexCount_ = static_cast<int>(Q.size());
    quadFaces_ = grid.quads();
    verts_ = Q;

    const int nE = static_cast<int>(D.edges.size());
    std::vector<int> chordOf(nE, -1);
    for (size_t c = 0; c < grid.chords().size(); ++c)
        for (int e : grid.chords()[c].edges) chordOf[e] = static_cast<int>(c);

    // Orientation. `along(b, k)` is +1 when block b's grid axis k % 2 increases
    // along macro edge edges[k] from its `from` to its `to` (BlockQuadMesh lays
    // side 0 forwards along i, side 1 forwards along j, sides 2 and 3
    // backwards; a flipped side reverses that). A chord is oriented by giving
    // each of its macro edges a sign -- whether `from` -> `to` runs the chord's
    // way -- and a block axis then runs the chord's way exactly when that sign
    // times along() is +1 on both of its sides. The signs are spread by a walk
    // from one edge per chord, across blocks; a contradiction would be a chord
    // that comes back on itself reversed, which a planar layout cannot have.
    auto along = [&](int b, int k) {
        const int s = k < 2 ? 1 : -1;
        return D.blocks[b].flip[k] ? -s : s;
    };
    std::vector<std::vector<std::pair<int, int>>> sidesOf(nE);   // (block, side) per macro edge
    for (size_t b = 0; b < D.blocks.size(); ++b)
        for (int k = 0; k < 4; ++k) sidesOf[D.blocks[b].edges[k]].push_back({static_cast<int>(b), k});

    std::vector<int> edgeSign(nE, 0);
    std::vector<std::array<int, 2>> axisSign(D.blocks.size(), {{0, 0}});
    for (int seed = 0; seed < nE; ++seed) {
        if (edgeSign[seed] != 0 || sidesOf[seed].empty()) continue;
        edgeSign[seed] = 1;
        std::vector<int> stack{seed};
        while (!stack.empty()) {
            const int e = stack.back();
            stack.pop_back();
            for (const auto &[b, k] : sidesOf[e]) {
                const int o = edgeSign[e] * along(b, k);
                int &axis = axisSign[b][k % 2];
                if (axis != 0 && axis != o)
                    throw std::runtime_error("ChartDistanceField: a chord through block " + std::to_string(b) +
                                             " cannot be oriented");
                axis = o;
                const int k2 = (k + 2) % 4, e2 = D.blocks[b].edges[k2];
                const int want = o * along(b, k2);
                if (edgeSign[e2] == 0) {
                    edgeSign[e2] = want;
                    stack.push_back(e2);
                } else if (edgeSign[e2] != want) {
                    throw std::runtime_error("ChartDistanceField: a chord through block " + std::to_string(b) +
                                             " cannot be oriented");
                }
            }
        }
    }

    // The split of each chord, from its macro edges' node spacing on the
    // block mesh: cell `cell` (counted the chord's way) at fraction `t` of it.
    struct Split {
        int cell = -1;
        double t = 0.5;
    };
    std::vector<Split> splits(grid.chords().size());
    for (size_t c = 0; c < grid.chords().size(); ++c) {
        std::vector<double> profile;
        int arcs = 0;
        for (int e : grid.chords()[c].edges) {
            std::vector<int> nodes = grid.edgeVertices()[e];
            if (nodes.size() < 2 || edgeSign[e] == 0) continue;
            if (edgeSign[e] < 0) std::reverse(nodes.begin(), nodes.end());
            const size_t n = nodes.size() - 1;
            if (profile.empty()) profile.assign(n, 0.0);
            if (profile.size() != n)
                throw std::runtime_error("ChartDistanceField: a chord's macro edges differ in interval count");
            std::vector<double> cum(n);
            double sum = 0.0;
            for (size_t i = 0; i < n; ++i) cum[i] = (sum += normP(Q[nodes[i + 1]] - Q[nodes[i]]));
            if (!(sum > 0.0)) continue;
            for (size_t i = 0; i < n; ++i) profile[i] += cum[i] / sum;
            ++arcs;
        }
        if (arcs == 0) throw std::runtime_error("ChartDistanceField: a chord has no macro edge of positive length");
        // CDFGen's rule (patchsample.cpp 818-833): the first cell whose end the
        // mean profile reaches 1/2 at, and the fraction of that cell before the
        // crossing, kept off both ends of the cell.
        const double half = 0.5 * profile.back();
        for (size_t i = 0; i < profile.size(); ++i) {
            if (profile[i] >= half) {
                const double width = i > 0 ? profile[i] - profile[i - 1] : profile[i];
                const double after = width > 0.0 ? (profile[i] - half) / width : 0.5;
                splits[c].cell = static_cast<int>(i);
                splits[c].t = 1.0 - std::clamp(after, 0.01, 0.99);
                break;
            }
        }
        if (splits[c].cell < 0) splits[c] = {static_cast<int>(profile.size()) - 1, 0.5};
    }

    // Two blocks reach a shared split vertex from the two ends of its edge, so
    // their positions differ by rounding; a wrong cell or fraction is off by a
    // fraction of an edge, many orders more.
    const double tol = 1e-6 * report_.target;
    std::unordered_map<std::uint64_t, int> splitVertexOf;   // block-mesh edge (lo, hi) -> vertex
    auto splitVertex = [&](int a, int b, double t) -> int {
        const Point p = Q[a] + (Q[b] - Q[a]) * t;
        const std::uint64_t key = a < b ? directedKey(a, b) : directedKey(b, a);
        const auto it = splitVertexOf.find(key);
        if (it != splitVertexOf.end()) {
            if (normP(verts_[it->second] - p) > tol) ++report_.splitMismatches;
            return it->second;
        }
        const int v = static_cast<int>(verts_.size());
        verts_.push_back(p);
        splitVertexOf.emplace(key, v);
        return v;
    };

    grids_.clear();
    grids_.reserve(grid.blocks().size());
    std::vector<int> cornerVertices;
    for (const BlockQuadMesh::Block &g : grid.blocks()) {
        const int b = g.block, ns = g.ns, nt = g.nt;
        auto Qid = [&](int i, int j) { return g.vert[static_cast<size_t>(j) * (ns + 1) + i]; };

        Grid G;
        G.block = b;
        G.ns = ns;
        G.nt = nt;
        double area = 0.0;
        for (int j = 0; j < nt; ++j)
            for (int i = 0; i < ns; ++i)
                area += signedArea(Q[Qid(i, j)], Q[Qid(i + 1, j)], Q[Qid(i + 1, j + 1)]) +
                        signedArea(Q[Qid(i + 1, j + 1)], Q[Qid(i, j + 1)], Q[Qid(i, j)]);
        G.ccw = area >= 0.0;

        // The two cuts in this block's own grid: the chord's cell counted from
        // the end the axis starts at.
        double ts = 0.5, tt = 0.5;
        for (int axis = 0; axis < 2; ++axis) {
            const Split &s = splits[chordOf[D.blocks[b].edges[axis]]];
            const int n = axis == 0 ? ns : nt;
            const bool forward = axisSign[b][axis] >= 0;
            const int cell = forward ? s.cell : n - 1 - s.cell;
            const double t = forward ? s.t : 1.0 - s.t;
            if (axis == 0) G.ks = cell, ts = t;
            else G.kt = cell, tt = t;
        }

        G.v.assign(static_cast<size_t>(ns + 2) * (nt + 2), -1);
        for (int q = 0; q <= nt + 1; ++q) {
            for (int r = 0; r <= ns + 1; ++r) {
                const bool fr = r == G.ks + 1, fq = q == G.kt + 1;
                const int i = r <= G.ks ? r : r - 1, j = q <= G.kt ? q : q - 1;
                int v;
                if (!fr && !fq) {
                    v = Qid(i, j);
                } else if (fr && !fq) {
                    v = splitVertex(Qid(G.ks, j), Qid(G.ks + 1, j), ts);
                } else if (!fr) {
                    v = splitVertex(Qid(i, G.kt), Qid(i, G.kt + 1), tt);
                } else {
                    // The chart centre: the bilinear point of the cell both cuts
                    // cross, which is where its two cut lines meet.
                    const Point &p00 = Q[Qid(G.ks, G.kt)], &p10 = Q[Qid(G.ks + 1, G.kt)];
                    const Point &p11 = Q[Qid(G.ks + 1, G.kt + 1)], &p01 = Q[Qid(G.ks, G.kt + 1)];
                    v = static_cast<int>(verts_.size());
                    verts_.push_back(p00 * ((1 - ts) * (1 - tt)) + p10 * (ts * (1 - tt)) + p11 * (ts * tt) +
                                     p01 * ((1 - ts) * tt));
                }
                G.v[static_cast<size_t>(q) * (ns + 2) + r] = v;
            }
        }

        // Faces, each in the subchart of the corner its quadrant holds, with
        // the diagonal through its vertex nearest that corner.
        const int gi = static_cast<int>(grids_.size());
        G.faces.reserve(static_cast<size_t>(ns + 1) * (nt + 1));
        for (int q = 0; q <= nt; ++q) {
            for (int r = 0; r <= ns; ++r) {
                const int a = G.at(r, q), bb = G.at(r + 1, q), c = G.at(r + 1, q + 1), d = G.at(r, q + 1);
                const std::array<int, 4> f = G.ccw ? std::array<int, 4>{a, bb, c, d} : std::array<int, 4>{a, d, c, bb};
                const int qs = r <= G.ks ? 0 : 1, qt = q <= G.kt ? 0 : 1;
                const int corner = qt == 0 ? qs : 3 - qs;   // (0,0) 0, (1,0) 1, (1,1) 2, (0,1) 3
                const int nearest = G.at(r + qs, q + qt);
                const int pos = static_cast<int>(std::find(f.begin(), f.end(), nearest) - f.begin());
                G.faces.push_back(static_cast<int>(faces_.size()));
                faces_.push_back(f);
                split_.push_back(static_cast<signed char>(pos & 1));
                faceSubchart_.push_back(4 * gi + corner);
            }
        }

        // The four subcharts: corner, side mid-points and centre, in the
        // order the face runs (CDFGen walks the dual face the same way), and
        // the chart's other three corners in that order from the first side.
        const std::array<int, 4> cornerV = {G.at(0, 0), G.at(ns + 1, 0), G.at(ns + 1, nt + 1), G.at(0, nt + 1)};
        const std::array<int, 4> midV = {G.at(G.ks + 1, 0), G.at(ns + 1, G.kt + 1), G.at(G.ks + 1, nt + 1),
                                         G.at(0, G.kt + 1)};
        const int centre = G.at(G.ks + 1, G.kt + 1);
        for (int k = 0; k < 4; ++k) {
            const int next = (k + 1) % 4, opp = (k + 2) % 4, prev = (k + 3) % 4;
            if (G.ccw) {
                abcd_.push_back({cornerV[k], midV[k], centre, midV[prev]});
                others_.push_back({cornerV[next], cornerV[opp], cornerV[prev]});
            } else {
                abcd_.push_back({cornerV[k], midV[prev], centre, midV[k]});
                others_.push_back({cornerV[prev], cornerV[opp], cornerV[next]});
            }
            subchartBlock_.push_back(b);
            cornerVertices.push_back(cornerV[k]);
        }
        grids_.push_back(std::move(G));
    }

    // One dual chart (patch) per vertex that is some block's corner, numbered
    // in vertex order as CDFGen numbers its corner vertices.
    patchVertex_ = cornerVertices;
    std::sort(patchVertex_.begin(), patchVertex_.end());
    patchVertex_.erase(std::unique(patchVertex_.begin(), patchVertex_.end()), patchVertex_.end());
    subchartPatch_.resize(cornerVertices.size());
    for (size_t s = 0; s < cornerVertices.size(); ++s)
        subchartPatch_[s] = static_cast<int>(
            std::lower_bound(patchVertex_.begin(), patchVertex_.end(), cornerVertices[s]) - patchVertex_.begin());

    report_.vertices = static_cast<int>(verts_.size());
    report_.faces = static_cast<int>(faces_.size());
    report_.subcharts = static_cast<int>(abcd_.size());
    report_.patches = static_cast<int>(patchVertex_.size());
}

// ---------------------------------------------------------------------------
// buildCoordinates()
//
// Steps 2 and 3 of the header: the ring-averaged length of every edge, the
// coordinate of every vertex on a half side, and every half-edge's reference
// into the table of (start, end) pairs.
// ---------------------------------------------------------------------------
void ChartDistanceField::buildCoordinates() {
    const int nF = static_cast<int>(faces_.size());
    HalfEdges he;
    he.pair.assign(4 * static_cast<size_t>(nF), -1);
    he.edge.assign(4 * static_cast<size_t>(nF), -1);
    he.byEnds.reserve(4 * static_cast<size_t>(nF));
    for (int f = 0; f < nF; ++f) {
        for (int k = 0; k < 4; ++k) {
            const int h = 4 * f + k;
            if (!he.byEnds.emplace(directedKey(faces_[f][k], faces_[f][(k + 1) % 4]), h).second)
                throw std::runtime_error("ChartDistanceField: the split mesh is not an oriented manifold "
                                         "(an edge is used twice in one direction)");
        }
    }
    for (int h = 0; h < 4 * nF; ++h) {
        if (he.edge[h] >= 0) continue;
        const int f = h / 4, k = h % 4;
        const int p = he.find(faces_[f][(k + 1) % 4], faces_[f][k]);
        he.pair[h] = p;
        he.edge[h] = he.edges;
        if (p >= 0) {
            he.pair[p] = h;
            he.edge[p] = he.edges;
        }
        ++he.edges;
    }
    auto tail = [&](int h) { return faces_[h / 4][h % 4]; };
    auto head = [&](int h) { return faces_[h / 4][(h % 4 + 1) % 4]; };

    // Boundary loops: each boundary half-edge is followed by the one leaving
    // its head.
    {
        std::unordered_map<int, int> leaving;
        for (int h = 0; h < 4 * nF; ++h)
            if (he.pair[h] < 0) leaving[tail(h)] = h;
        std::vector<char> seen(4 * static_cast<size_t>(nF), 0);
        for (int h = 0; h < 4 * nF; ++h) {
            if (he.pair[h] >= 0 || seen[h]) continue;
            ++report_.boundaryLoops;
            for (int cur = h; cur >= 0 && !seen[cur];) {
                seen[cur] = 1;
                const auto it = leaving.find(head(cur));
                cur = it == leaving.end() ? -1 : it->second;
            }
        }
    }

    // Step 2: the ring-averaged length of each edge. A ring is grown from an
    // edge both ways, across the face on each side to the opposite edge and on
    // through that edge's other face, until it reaches dS or closes.
    std::vector<double> ringLength(he.edges, 0.0);
    {
        std::vector<int> firstHalf(he.edges, -1);
        for (int h = 0; h < 4 * nF; ++h)
            if (firstHalf[he.edge[h]] < 0) firstHalf[he.edge[h]] = h;
        std::vector<char> done(he.edges, 0);
        std::vector<int> ring;
        for (int e = 0; e < he.edges; ++e) {
            if (done[e]) continue;
            ring.assign(1, e);
            done[e] = 1;
            const int h0 = firstHalf[e];
            for (int start : {h0, he.pair[h0]}) {
                for (int cur = start; cur >= 0;) {
                    const int o = HalfEdges::opposite(cur);
                    if (done[he.edge[o]]) break;
                    done[he.edge[o]] = 1;
                    ring.push_back(he.edge[o]);
                    cur = he.pair[o];
                }
            }
            double sum = 0.0;
            for (int r : ring) sum += normP(verts_[head(firstHalf[r])] - verts_[tail(firstHalf[r])]);
            for (int r : ring) ringLength[r] = sum / static_cast<double>(ring.size());
        }
    }
    auto edgeOf = [&](int a, int b) {
        int h = he.find(a, b);
        if (h < 0) h = he.find(b, a);
        return h < 0 ? -1 : he.edge[h];
    };

    // Step 3a: the coordinate of each vertex on a half side, from its corner
    // (1) to its side's split vertex (0), by remaining ring length.
    std::vector<double> coordinate(verts_.size(), std::numeric_limits<double>::quiet_NaN());
    auto setCoordinate = [&](int v, double c) {
        if (std::isnan(coordinate[v])) coordinate[v] = c;
        else if (std::fabs(coordinate[v] - c) > 1e-9) ++report_.coordinateMismatches;
    };
    // The half sides of a grid, as refined-grid walks from a corner: corner k
    // along its two sides, k's own side (toward corner k + 1) first.
    auto halfSides = [](const Grid &G, int k) {
        const int S = G.ns + 1, T = G.nt + 1, ms = G.ks + 1, mt = G.kt + 1;
        std::array<std::vector<std::array<int, 2>>, 2> w;
        auto walk = [](int r0, int q0, int r1, int q1) {
            std::vector<std::array<int, 2>> out;
            const int dr = (r1 > r0) - (r1 < r0), dq = (q1 > q0) - (q1 < q0);
            for (int r = r0, q = q0;; r += dr, q += dq) {
                out.push_back({r, q});
                if (r == r1 && q == q1) break;
            }
            return out;
        };
        switch (k) {
        case 0: w = {walk(0, 0, ms, 0), walk(0, 0, 0, mt)}; break;
        case 1: w = {walk(S, 0, S, mt), walk(S, 0, ms, 0)}; break;
        case 2: w = {walk(S, T, ms, T), walk(S, T, S, mt)}; break;
        default: w = {walk(0, T, 0, mt), walk(0, T, ms, T)}; break;
        }
        return w;
    };
    for (const Grid &G : grids_) {
        for (int k = 0; k < 4; ++k) {
            for (const auto &walk : halfSides(G, k)) {
                std::vector<double> len(walk.size() - 1);
                double total = 0.0;
                for (size_t i = 0; i + 1 < walk.size(); ++i) {
                    const int e = edgeOf(G.at(walk[i][0], walk[i][1]), G.at(walk[i + 1][0], walk[i + 1][1]));
                    len[i] = e >= 0 ? ringLength[e] : 0.0;
                    total += len[i];
                }
                // The split vertex is exactly 0, as CDFGen leaves it, not the
                // rounding residue of the subtraction.
                double rest = total;
                for (size_t i = 0; i + 1 < walk.size(); ++i) {
                    setCoordinate(G.at(walk[i][0], walk[i][1]), total > 0.0 ? rest / total : 0.0);
                    rest -= len[i];
                }
                setCoordinate(G.at(walk.back()[0], walk.back()[1]), 0.0);
            }
        }
    }

    // Step 3b: one (start, end) pair per half-side edge, read along the
    // edge's first half-edge; then every half-edge of a subchart refers to the
    // half-side edge of its own column (an s-edge) or row (a t-edge) on its
    // subchart's side, signed by whether it runs the same way.
    pairs_.assign(1, {0.0, 0.0});
    std::vector<int> pairOfEdge(he.edges, 0), pairTail(1, -1);
    auto sidePair = [&](int a, int b) -> std::pair<int, int> {   // (pair index, its tail)
        int h = he.find(a, b);
        if (h < 0) h = he.find(b, a);
        if (h < 0) return {0, -1};
        const int e = he.edge[h];
        if (pairOfEdge[e] == 0) {
            const int t = tail(h), hd = head(h);
            pairOfEdge[e] = static_cast<int>(pairs_.size());
            pairs_.push_back({coordinate[t], coordinate[hd]});
            pairTail.push_back(t);
        }
        return {pairOfEdge[e], pairTail[pairOfEdge[e]]};
    };

    halfEdgePair_.assign(4 * static_cast<size_t>(nF), 0);
    for (const Grid &G : grids_) {
        const int ns = G.ns, nt = G.nt;
        for (int q = 0; q <= nt; ++q) {
            for (int r = 0; r <= ns; ++r) {
                const int f = G.faces[static_cast<size_t>(q) * (ns + 1) + r];
                const int qs = r <= G.ks ? 0 : 1, qt = q <= G.kt ? 0 : 1;
                const std::array<std::array<int, 2>, 4> ccw = {{{r, q}, {r + 1, q}, {r + 1, q + 1}, {r, q + 1}}};
                const std::array<std::array<int, 2>, 4> pos =
                    G.ccw ? ccw : std::array<std::array<int, 2>, 4>{ccw[0], ccw[3], ccw[2], ccw[1]};
                for (int k = 0; k < 4; ++k) {
                    const std::array<int, 2> &p0 = pos[k], &p1 = pos[(k + 1) % 4];
                    int a, b;
                    bool forward;
                    if (p0[1] == p1[1]) {        // along s: the side row of this quadrant
                        const int c = std::min(p0[0], p1[0]), row = qt == 0 ? 0 : nt + 1;
                        a = G.at(c, row);
                        b = G.at(c + 1, row);
                        forward = p1[0] > p0[0];
                    } else {                     // along t: the side column
                        const int c = std::min(p0[1], p1[1]), col = qs == 0 ? 0 : ns + 1;
                        a = G.at(col, c);
                        b = G.at(col, c + 1);
                        forward = p1[1] > p0[1];
                    }
                    const auto [index, t] = sidePair(a, b);
                    if (index == 0) {
                        ++report_.missingPairs;
                        continue;
                    }
                    halfEdgePair_[4 * static_cast<size_t>(f) + k] = forward == (t == a) ? index : -index;
                }
            }
        }
    }
    report_.coordinatePairs = static_cast<int>(pairs_.size());
}

std::array<int, 3> ChartDistanceField::triangle(int t) const {
    const std::array<int, 4> &q = faces_[t / 2];
    if (!split_[t / 2]) return t % 2 == 0 ? std::array<int, 3>{q[0], q[1], q[2]} : std::array<int, 3>{q[2], q[3], q[0]};
    return t % 2 == 0 ? std::array<int, 3>{q[3], q[0], q[1]} : std::array<int, 3>{q[1], q[2], q[3]};
}

Point ChartDistanceField::pointAt(int t, const std::array<double, 3> &w) const {
    const std::array<int, 3> v = triangle(t);
    return verts_[v[0]] * w[0] + verts_[v[1]] * w[1] + verts_[v[2]] * w[2];
}

// ---------------------------------------------------------------------------
// sample()
//
// Area-uniform samples, and farthest-point token sets over a larger draw of
// them. One generator, seeded, drawn in a fixed order.
// ---------------------------------------------------------------------------
void ChartDistanceField::sample() {
    const int nT = triangleCount();
    std::vector<double> cum(nT);
    double total = 0.0;
    report_.minTriangleArea = std::numeric_limits<double>::infinity();
    std::vector<char> folded(faces_.size(), 0);
    for (int t = 0; t < nT; ++t) {
        const std::array<int, 3> v = triangle(t);
        const double a = signedArea(verts_[v[0]], verts_[v[1]], verts_[v[2]]);
        report_.minTriangleArea = std::min(report_.minTriangleArea, a);
        if (!(a > 0.0)) folded[t / 2] = 1;
        cum[t] = (total += std::max(0.0, a));
    }
    report_.area = total;
    report_.invertedFaces = static_cast<int>(std::count(folded.begin(), folded.end(), 1));
    if (!(total > 0.0)) throw std::runtime_error("ChartDistanceField: the split mesh has no positive area");

    std::mt19937_64 gen(opts_.seed);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    auto draw = [&](int &tri, std::array<double, 3> &w) {
        const double u = U(gen) * total;
        tri = static_cast<int>(std::min<std::ptrdiff_t>(std::upper_bound(cum.begin(), cum.end(), u) - cum.begin(), nT - 1));
        const double s = std::sqrt(U(gen)), r = U(gen);
        w = {1.0 - s, s * (1.0 - r), s * r};
    };

    const int n = std::max(0, opts_.samples);
    sampleTri_.resize(n);
    sampleBary_.resize(n);
    for (int i = 0; i < n; ++i) draw(sampleTri_[i], sampleBary_[i]);
    report_.samples = n;

    fps_.clear();
    if (opts_.fpsMin <= 0 || opts_.fpsMax < opts_.fpsMin) return;
    std::vector<int> candTri;
    std::vector<std::array<double, 3>> candBary;
    std::vector<Point> candP;
    std::vector<double> dist;
    for (int size = opts_.fpsMin; size <= opts_.fpsMax; size *= 2) {
        TokenSets &set = fps_[size];
        const int m = std::max(size, std::max(opts_.fpsOversample * size, opts_.fpsMinCandidates));
        for (int copy = 0; copy < std::max(1, opts_.fpsCopies); ++copy) {
            candTri.resize(m);
            candBary.resize(m);
            candP.resize(m);
            for (int i = 0; i < m; ++i) {
                draw(candTri[i], candBary[i]);
                candP[i] = pointAt(candTri[i], candBary[i]);
            }
            dist.assign(m, std::numeric_limits<double>::infinity());
            std::vector<int> tri(size);
            std::vector<std::array<double, 3>> bary(size);
            int pick = std::min(m - 1, static_cast<int>(U(gen) * m));
            for (int i = 0; i < size; ++i) {
                tri[i] = candTri[pick];
                bary[i] = candBary[pick];
                const Point p = candP[pick];
                int far = 0;
                double farD = -1.0;
                for (int j = 0; j < m; ++j) {
                    const double dx = candP[j][0] - p[0], dy = candP[j][1] - p[1];
                    const double d = std::min(dist[j], dx * dx + dy * dy);
                    dist[j] = d;
                    if (d > farD) farD = d, far = j;
                }
                pick = far;
            }
            set.triangles.push_back(std::move(tri));
            set.bary.push_back(std::move(bary));
        }
        ++report_.tokenSets;
        if (size > std::numeric_limits<int>::max() / 2) break;
    }
}

// ---------------------------------------------------------------------------
// evaluate()
//
// CDFGen's get_grading_color_all (patchsample.cpp 1536) in the plane: the
// triangle's two quad edges, each with its pair, the coordinate along each by
// the oblique projection parallel to the other, and the gradient of that
// projection; then the release loader's min/max and its tie-breaks.
// ---------------------------------------------------------------------------
ChartDistanceField::Fields ChartDistanceField::evaluate(int t, const Point &p) const {
    Fields out;
    const int f = t / 2, k = t % 2;
    const std::array<int, 3> v = triangle(t);
    const int ea = split_[f] ? (k == 0 ? 3 : 1) : (k == 0 ? 0 : 2);
    const int eb = (ea + 1) % 4;
    auto pairOf = [&](int e) {
        const int id = halfEdgePair_[4 * static_cast<size_t>(f) + e];
        const std::array<double, 2> &c = pairs_[std::abs(id)];
        return id >= 0 ? c : std::array<double, 2>{c[1], c[0]};
    };
    const std::array<double, 2> ca = pairOf(ea), cb = pairOf(eb);
    const Point E0 = verts_[v[1]] - verts_[v[0]], E1 = verts_[v[2]] - verts_[v[1]];
    const double den = cross2(E0, E1);
    const double t0 = std::clamp(cross2(p - verts_[v[0]], E1) / den, 0.0, 1.0);
    const double t1 = std::clamp(cross2(p - verts_[v[1]], E0) / -den, 0.0, 1.0);
    out.c0 = ca[0] + t0 * (ca[1] - ca[0]);
    out.c1 = cb[0] + t1 * (cb[1] - cb[0]);
    auto perp = [](const Point &a) { return Point{a[1], -a[0]}; };
    auto unit = [](const Point &a) { return a / (normP(a) + 1e-10); };
    const Point g0 = unit(perp(E1) * ((ca[1] - ca[0]) / den));
    const Point g1 = unit(perp(E0) * ((cb[1] - cb[0]) / -den));
    out.dcdf = std::clamp(std::min(out.c0, out.c1), 0.0, 1.0);
    out.gdcdf = out.c1 < out.c0 ? g1 : g0;
    out.cdf = std::clamp(1.0 - std::max(out.c0, out.c1), 0.0, 1.0);
    out.gcdf = (1.0 - out.c1 < 1.0 - out.c0) ? g1 * -1.0 : g0 * -1.0;
    const int s = faceSubchart_[f];
    for (int i = 0; i < 4; ++i) out.offsets[i] = verts_[abcd_[s][i]] - p;
    for (int i = 0; i < 3; ++i) out.offsets[4 + i] = verts_[others_[s][i]] - p;
    return out;
}

// ---------------------------------------------------------------------------
// addTo()
//
// The release's keys, shapes and dtypes (patchsample.cpp 223-528): uint32
// indices, float32 coordinates, the (n, 1) columns as (n, 1). edge_faceids and
// edge_info, which only QuadExtraction's debugging reads, are left out.
// ---------------------------------------------------------------------------
void ChartDistanceField::addTo(NpzWriter &out) const {
    using u32 = std::uint32_t;
    auto one = [](size_t n) { return std::vector<size_t>{n, 1}; };
    const size_t nV = verts_.size(), nF = faces_.size(), nS = abcd_.size();

    out.add<std::int8_t>("is_checkerboard", {1, 1}, {2});
    out.add<u32>("num_complex", {1, 1}, {static_cast<u32>(report_.blocks)});
    out.add<u32>("num_boundary", {1, 1}, {static_cast<u32>(report_.boundaryLoops)});
    out.add<float>("invT", {3, 4}, {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0});

    out.add<u32>("quad_mesh_vertex_num", {1, 1}, {static_cast<u32>(quadVertexCount_)});
    std::vector<u32> quadFacet;
    quadFacet.reserve(4 * quadFaces_.size());
    for (const auto &q : quadFaces_)
        for (int v : q) quadFacet.push_back(static_cast<u32>(v));
    out.add<u32>("quad_facet", {quadFaces_.size(), 4}, quadFacet);

    std::vector<float> xyz;
    xyz.reserve(3 * nV);
    for (const Point &p : verts_) xyz.insert(xyz.end(), {static_cast<float>(p[0]), static_cast<float>(p[1]), 0.0f});
    out.add<float>("subdiv_vertex", {nV, 3}, xyz);
    std::vector<u32> facet;
    facet.reserve(4 * nF);
    for (const auto &q : faces_)
        for (int v : q) facet.push_back(static_cast<u32>(v));
    out.add<u32>("subdiv_facet", {nF, 4}, facet);
    out.add<u32>("face2quad", one(nF), std::vector<u32>(faceSubchart_.begin(), faceSubchart_.end()));
    out.add<u32>("quad2patch", one(nS), std::vector<u32>(subchartPatch_.begin(), subchartPatch_.end()));
    out.add<std::int8_t>("quad_split", one(nF), std::vector<std::int8_t>(split_.begin(), split_.end()));

    std::vector<u32> abcd, others;
    for (const auto &a : abcd_)
        for (int v : a) abcd.push_back(static_cast<u32>(v));
    for (const auto &a : others_)
        for (int v : a) others.push_back(static_cast<u32>(v));
    out.add<u32>("offset_abcd_id", {nS, 4}, abcd);
    out.add<u32>("offset_123_id", {nS, 3}, others);

    out.add<std::int32_t>("edge_vertex_color", one(4 * nF),
                          std::vector<std::int32_t>(halfEdgePair_.begin(), halfEdgePair_.end()));
    std::vector<float> store;
    for (const auto &c : pairs_) store.insert(store.end(), {static_cast<float>(c[0]), static_cast<float>(c[1])});
    out.add<float>("edge_color_store", {pairs_.size(), 2}, store);

    std::vector<float> bary;
    for (const auto &w : sampleBary_)
        bary.insert(bary.end(), {static_cast<float>(w[0]), static_cast<float>(w[1]), static_cast<float>(w[2])});
    out.add<float>("sample_bary_coords", {sampleBary_.size(), 3}, bary);
    out.add<u32>("sample_point_tri_id", one(sampleTri_.size()), std::vector<u32>(sampleTri_.begin(), sampleTri_.end()));

    for (const auto &[size, set] : fps_) {
        std::vector<float> w;
        std::vector<u32> tri;
        for (size_t c = 0; c < set.triangles.size(); ++c) {
            for (size_t i = 0; i < set.triangles[c].size(); ++i) {
                const auto &b = set.bary[c][i];
                w.insert(w.end(), {static_cast<float>(b[0]), static_cast<float>(b[1]), static_cast<float>(b[2])});
                tri.push_back(static_cast<u32>(set.triangles[c][i]));
            }
        }
        const size_t copies = set.triangles.size(), n = static_cast<size_t>(size);
        out.add<float>("fps_bary_coords_" + std::to_string(size), {copies, n, 3}, w);
        out.add<u32>("fps_point_tri_id_" + std::to_string(size), {copies, n, 1}, tri);
    }

    out.add<std::int32_t>("subchart_block", {nS}, std::vector<std::int32_t>(subchartBlock_.begin(), subchartBlock_.end()));
    out.add<std::int32_t>("patch_vertex", {patchVertex_.size()},
                          std::vector<std::int32_t>(patchVertex_.begin(), patchVertex_.end()));
    out.addScalar("target_edge_length", report_.target);
}
