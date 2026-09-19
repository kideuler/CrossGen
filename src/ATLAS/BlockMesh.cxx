#include "ATLAS/BlockMesh.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <set>
#include <unordered_map>

#include "geom/Coons.hxx"

namespace {

int findRoot(std::vector<int> &parent, int x) {
    while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
    return x;
}

void unite(std::vector<int> &parent, int a, int b) {
    a = findRoot(parent, a);
    b = findRoot(parent, b);
    if (a != b) parent[b] = a;
}

// The scaled Jacobian of a planar quadrilateral, the worst of its four corners,
// and its signed area: Stage 10's two measures, so the reports compare.
double scaledJacobian(const Point q[4]) {
    double worst = std::numeric_limits<double>::infinity();
    for (int i = 0; i < 4; ++i) {
        const Point a = q[(i + 1) % 4] - q[i];
        const Point b = q[(i + 3) % 4] - q[i];
        const double na = normP(a), nb = normP(b);
        if (!(na > 0.0) || !(nb > 0.0)) return 0.0;
        worst = std::min(worst, cross2(a, b) / (na * nb));
    }
    return worst;
}

double quadArea(const Point q[4]) {
    return 0.5 * (cross2(q[1] - q[0], q[2] - q[0]) + cross2(q[2] - q[0], q[3] - q[0]));
}

double segmentDistance(const Point &p, const Point &a, const Point &b) {
    const Point d = b - a;
    const double L2 = dotP(d, d);
    double t = L2 > 0.0 ? dotP(p - a, d) / L2 : 0.0;
    t = std::min(1.0, std::max(0.0, t));
    return normP(p - (a + d * t));
}

// Cumulative arc length along a chain of carrier vertices.
std::vector<double> chainLengths(const SquareCarrier &C, const std::vector<int> &chain) {
    std::vector<double> d(chain.size(), 0.0);
    for (size_t k = 1; k < chain.size(); ++k) {
        d[k] = d[k - 1] + normP(C.vertices[chain[k]] - C.vertices[chain[k - 1]]);
    }
    return d;
}

} // namespace

// ---------------------------------------------------------------------------
BlockMesh::BlockMesh(const BlockCover &cover, const Options &opts)
    : C_(cover.getCarrier()), opts_(opts) {
    report_.target = opts_.targetEdgeLength;
    const PlanarDomain &D = C_.getDomain();
    report_.modelExtent = D.getReport().scale > 0.0 ? D.getReport().scale : 1.0;
    report_.domainArea = D.getReport().area;
    if (!(opts_.targetEdgeLength > 0.0)) {
        report_.messages.push_back("the target edge length must be positive");
        return;
    }
    if (cover.getBlocks().empty()) {
        report_.messages.push_back("the cover has no blocks");
        return;
    }
    buildGrids(cover);
    assignIntervals(cover);
    meshArcs(cover);
    meshBlocks(cover);
    smooth();
    check(cover);
}

// ---------------------------------------------------------------------------
// buildGrids()
//
// Each block's chart as its carrier vertex at every integer point, read off
// the certificate: cell cells[i + n_u j] sits on the unit square at (i, j), and
// G carries its corner c, at squareCorner(c) in the cell's own coordinates, to
// the block's. And each side's macro edge: BlockCover makes one per block
// side, its chain running counter-clockwise round blockA, so for blockB it
// runs the other way; which way is read off the ends rather than assumed.
// ---------------------------------------------------------------------------
void BlockMesh::buildGrids(const BlockCover &cover) {
    const std::vector<BlockCover::Block> &blocks = cover.getBlocks();
    const std::vector<BlockCover::MacroEdge> &edges = cover.getMacroEdges();
    const int nb = static_cast<int>(blocks.size());
    chart_.assign(nb, {});
    sideArc_.assign(nb, {{-1, -1, -1, -1}});
    sideForward_.assign(nb, {{1, 1, 1, 1}});

    for (int b = 0; b < nb; ++b) {
        const RectangleCertifier::Certificate &c = blocks[b].cert;
        std::vector<int> &g = chart_[b];
        g.assign(static_cast<size_t>(c.nu + 1) * (c.nv + 1), -1);
        for (size_t pos = 0; pos < c.cells.size() && pos < c.G.size(); ++pos) {
            const int q = c.cells[pos];
            if (q < 0 || q >= C_.numCells()) continue;
            for (int k = 0; k < 4; ++k) {
                const IPoint p = c.G[pos].apply(squareCorner(k));
                if (p[0] < 0 || p[0] > c.nu || p[1] < 0 || p[1] > c.nv) continue;
                g[static_cast<size_t>(p[1]) * (c.nu + 1) + p[0]] = C_.cells[q][k];
            }
        }
    }

    for (int a = 0; a < static_cast<int>(edges.size()); ++a) {
        const BlockCover::MacroEdge &me = edges[a];
        const int ends[2][2] = {{me.blockA, me.sideA}, {me.blockB, me.sideB}};
        for (const auto &e : ends) {
            const int b = e[0], s = e[1];
            if (b < 0 || b >= nb || s < 0 || s > 3 || me.chain.empty()) continue;
            const std::vector<int> &side = blocks[b].cert.sides[s];
            sideArc_[b][s] = a;
            sideForward_[b][s] = (!side.empty() && side.front() == me.chain.front()) ? 1 : 0;
        }
    }
}

// ---------------------------------------------------------------------------
// assignIntervals()
//
// Sec. 11.2 with Stage 10's objective; see the header. A block contributes the
// mean length of its u-lines to the chord of its sides 0 and 2, and of its
// v-lines to the chord of its sides 1 and 3.
// ---------------------------------------------------------------------------
void BlockMesh::assignIntervals(const BlockCover &cover) {
    const std::vector<BlockCover::Block> &blocks = cover.getBlocks();
    const std::vector<BlockCover::MacroEdge> &edges = cover.getMacroEdges();
    const int na = static_cast<int>(edges.size());
    const int nb = static_cast<int>(blocks.size());
    const double h = opts_.targetEdgeLength;

    std::vector<int> parent(na);
    std::iota(parent.begin(), parent.end(), 0);
    for (int b = 0; b < nb; ++b) {
        const std::array<int, 4> &s = sideArc_[b];
        if (s[0] >= 0 && s[2] >= 0) unite(parent, s[0], s[2]);
        if (s[1] >= 0 && s[3] >= 0) unite(parent, s[1], s[3]);
    }
    std::unordered_map<int, int> classOf;
    chordOf_.assign(na, -1);
    for (int a = 0; a < na; ++a) {
        const int r = findRoot(parent, a);
        auto it = classOf.find(r);
        if (it == classOf.end()) {
            it = classOf.emplace(r, static_cast<int>(chords_.size())).first;
            chords_.push_back(Chord());
        }
        chordOf_[a] = it->second;
        chords_[it->second].arcs.push_back(a);
    }

    std::vector<std::vector<double>> spans(chords_.size());
    for (int b = 0; b < nb; ++b) {
        const RectangleCertifier::Certificate &c = blocks[b].cert;
        const std::vector<int> &g = chart_[b];
        auto at = [&](int i, int j) { return g[static_cast<size_t>(j) * (c.nu + 1) + i]; };
        bool complete = true;
        for (int v : g) if (v < 0) complete = false;
        if (!complete || c.nu < 1 || c.nv < 1) continue;
        double su = 0.0, sv = 0.0;
        for (int j = 0; j <= c.nv; ++j) {
            for (int i = 0; i < c.nu; ++i) su += normP(C_.vertices[at(i + 1, j)] - C_.vertices[at(i, j)]);
        }
        for (int i = 0; i <= c.nu; ++i) {
            for (int j = 0; j < c.nv; ++j) sv += normP(C_.vertices[at(i, j + 1)] - C_.vertices[at(i, j)]);
        }
        su /= (c.nv + 1);
        sv /= (c.nu + 1);
        if (sideArc_[b][0] >= 0) spans[chordOf_[sideArc_[b][0]]].push_back(su);
        if (sideArc_[b][1] >= 0) spans[chordOf_[sideArc_[b][1]]].push_back(sv);
    }

    intervals_.assign(na, 0);
    long long total = 0;
    report_.minIntervals = std::numeric_limits<int>::max();
    report_.maxIntervals = 0;
    for (size_t k = 0; k < chords_.size(); ++k) {
        Chord &ch = chords_[k];
        ch.minLength = std::numeric_limits<double>::infinity();
        for (int a : ch.arcs) {
            const std::vector<double> d = chainLengths(C_, edges[a].chain);
            const double L = d.empty() ? 0.0 : d.back();
            ch.minLength = std::min(ch.minLength, L);
            ch.maxLength = std::max(ch.maxLength, L);
        }
        std::vector<double> &S = spans[k];
        if (S.empty()) {
            // Only edges of blocks whose chart did not build: take the edges.
            S.push_back(ch.maxLength);
        }
        double logSum = 0.0;
        ch.minSpan = std::numeric_limits<double>::infinity();
        for (double s : S) {
            const double sp = std::max(s, 1e-300);
            logSum += std::log(sp);
            ch.minSpan = std::min(ch.minSpan, sp);
            ch.maxSpan = std::max(ch.maxSpan, sp);
        }
        ch.idealIntervals = std::exp(logSum / S.size()) / h;
        auto F = [&](int N) {
            double f = 0.0;
            for (double s : S) {
                const double r = std::log(std::max(s, 1e-300) / (N * h));
                f += r * r;
            }
            return f;
        };
        const int lo = std::max(1, static_cast<int>(std::floor(ch.idealIntervals)));
        int N = F(lo + 1) < F(lo) ? lo + 1 : lo;
        const int floorN = std::max(1, opts_.minIntervals);
        int clampedN = std::max(N, floorN);
        if (opts_.maxIntervals > 0) clampedN = std::min(clampedN, std::max(floorN, opts_.maxIntervals));
        ch.clamped = clampedN != N;
        ch.intervals = clampedN;
        if (ch.clamped) ++report_.clampedChords;
        for (int a : ch.arcs) intervals_[a] = ch.intervals;
        report_.minIntervals = std::min(report_.minIntervals, ch.intervals);
        report_.maxIntervals = std::max(report_.maxIntervals, ch.intervals);
        total += ch.intervals;
    }
    report_.chords = static_cast<int>(chords_.size());
    report_.arcsAssigned = na;
    if (chords_.empty()) report_.minIntervals = 0;
    report_.meanIntervals = chords_.empty() ? 0.0 : static_cast<double>(total) / chords_.size();
}

// ---------------------------------------------------------------------------
// meshArcs()
//
// The nodes of every macro edge, once, at equal arc length along its chain, so
// the two blocks that share it get the same vertices and the mesh is
// conforming by construction. Corners -- the macrovertices -- are made once
// and shared by every edge that ends there. Each node also keeps where on the
// chain it fell, in carrier edges, which is its chart coordinate along the
// side (a block side's chain steps one unit of the chart per carrier edge).
// ---------------------------------------------------------------------------
void BlockMesh::meshArcs(const BlockCover &cover) {
    const std::vector<BlockCover::MacroEdge> &edges = cover.getMacroEdges();
    const int na = static_cast<int>(edges.size());
    std::unordered_map<int, int> corner;
    auto cornerVertex = [&](int v) {
        auto it = corner.find(v);
        if (it != corner.end()) return it->second;
        const int id = static_cast<int>(verts_.size());
        verts_.push_back(C_.vertices[v]);
        corner.emplace(v, id);
        return id;
    };

    arcNodes_.assign(na, {});
    arcParams_.assign(na, {});
    for (int a = 0; a < na; ++a) {
        const std::vector<int> &chain = edges[a].chain;
        const int N = intervals_[a];
        if (chain.size() < 2 || N < 1) continue;
        const std::vector<double> d = chainLengths(C_, chain);
        const double L = d.back();
        std::vector<int> &nodes = arcNodes_[a];
        std::vector<double> &par = arcParams_[a];
        nodes.assign(N + 1, -1);
        par.assign(N + 1, 0.0);
        nodes[0] = cornerVertex(chain.front());
        nodes[N] = cornerVertex(chain.back());
        par[0] = 0.0;
        par[N] = static_cast<double>(chain.size() - 1);
        size_t m = 0;
        for (int k = 1; k < N; ++k) {
            const double s = L * k / N;
            while (m + 2 < d.size() && d[m + 1] < s) ++m;
            const double seg = d[m + 1] - d[m];
            const double f = seg > 0.0 ? std::min(1.0, std::max(0.0, (s - d[m]) / seg)) : 0.0;
            const Point &p0 = C_.vertices[chain[m]], &p1 = C_.vertices[chain[m + 1]];
            nodes[k] = static_cast<int>(verts_.size());
            verts_.push_back(p0 + (p1 - p0) * f);
            par[k] = static_cast<double>(m) + f;
        }
    }
}

// ---------------------------------------------------------------------------
// meshBlocks()
//
// One (ns+1) x (nt+1) grid per block, its rim looked up from the four macro
// edges and its interior interpolated. A side's node m (from the side's
// start) and its chart coordinate sigma along the side go to
//
//     side 0:  (m, 0),        chart (sigma, 0)
//     side 1:  (ns, m),       chart (n_u, sigma)
//     side 2:  (ns - m, nt),  chart (n_u - sigma, n_v)
//     side 3:  (0, nt - m),   chart (0, n_v - sigma)
//
// and the interior node (i, j) is the chart evaluated at the transfinite
// blend of those coordinates -- s from the bottom and top, t from the left and
// right, exactly as Stage 10 blends its sides' spline parameters.
// ---------------------------------------------------------------------------
void BlockMesh::meshBlocks(const BlockCover &cover) {
    const std::vector<BlockCover::Block> &blocks = cover.getBlocks();
    const int nb = static_cast<int>(blocks.size());
    std::set<int> materials;

    for (int b = 0; b < nb; ++b) {
        const RectangleCertifier::Certificate &c = blocks[b].cert;
        const std::array<int, 4> &sa = sideArc_[b];
        bool ok = c.nu >= 1 && c.nv >= 1;
        for (int s = 0; s < 4 && ok; ++s) {
            if (sa[s] < 0 || arcNodes_[sa[s]].empty()) ok = false;
        }
        if (ok && (intervals_[sa[0]] != intervals_[sa[2]] || intervals_[sa[1]] != intervals_[sa[3]])) {
            // The union-find made this impossible; if it happens the cover is
            // not what its macro edges say, and the block is left out rather
            // than meshed with a hanging node.
            ok = false;
            report_.messages.push_back("opposite sides of a block disagree on their count; the block was "
                                       "left unmeshed");
        }
        const std::vector<int> &g = chart_[b];
        if (ok) for (int v : g) if (v < 0) ok = false;
        if (!ok) {
            ++report_.unmeshedBlocks;
            continue;
        }

        Block B;
        B.block = b;
        B.material = c.material;
        const int ns = intervals_[sa[0]], nt = intervals_[sa[1]];
        B.ns = ns;
        B.nt = nt;
        B.vert.assign(static_cast<size_t>(ns + 1) * (nt + 1), -1);
        auto at = [&](int i, int j) -> int & { return B.vert[static_cast<size_t>(j) * (ns + 1) + i]; };

        // The rim, and each rim node's chart coordinate along its side.
        std::vector<double> xB(ns + 1), xT(ns + 1), yL(nt + 1), yR(nt + 1);
        for (int s = 0; s < 4; ++s) {
            const int a = sa[s];
            const int N = intervals_[a];
            const double len = static_cast<double>(c.sides[s].size() - 1);
            for (int m = 0; m <= N; ++m) {
                const int k = sideForward_[b][s] ? m : N - m;
                const int v = arcNodes_[a][k];
                const double sigma = sideForward_[b][s] ? arcParams_[a][k] : len - arcParams_[a][k];
                switch (s) {
                    case 0: at(m, 0) = v; xB[m] = sigma; break;
                    case 1: at(ns, m) = v; yR[m] = sigma; break;
                    case 2: at(ns - m, nt) = v; xT[ns - m] = c.nu - sigma; break;
                    default: at(0, nt - m) = v; yL[nt - m] = c.nv - sigma; break;
                }
            }
        }

        auto gp = [&](int i, int j) -> const Point & {
            return C_.vertices[g[static_cast<size_t>(j) * (c.nu + 1) + i]];
        };
        auto chart = [&](double s, double t) {
            const int i0 = std::min(std::max(static_cast<int>(std::floor(s)), 0), c.nu - 1);
            const int j0 = std::min(std::max(static_cast<int>(std::floor(t)), 0), c.nv - 1);
            const double a = std::min(1.0, std::max(0.0, s - i0));
            const double e = std::min(1.0, std::max(0.0, t - j0));
            return gp(i0, j0) * ((1.0 - a) * (1.0 - e)) + gp(i0 + 1, j0) * (a * (1.0 - e)) +
                   gp(i0 + 1, j0 + 1) * (a * e) + gp(i0, j0 + 1) * ((1.0 - a) * e);
        };
        for (int j = 1; j < nt; ++j) {
            const double v = static_cast<double>(j) / nt;
            for (int i = 1; i < ns; ++i) {
                const double u = static_cast<double>(i) / ns;
                Point p;
                if (opts_.useChart) {
                    const double s = (1.0 - v) * xB[i] + v * xT[i];
                    const double t = (1.0 - u) * yL[j] + u * yR[j];
                    p = chart(s, t);
                } else {
                    p = geom::coonsPoint(verts_[at(i, 0)], verts_[at(i, nt)], verts_[at(0, j)], verts_[at(ns, j)],
                                         verts_[at(0, 0)], verts_[at(ns, 0)], verts_[at(0, nt)],
                                         verts_[at(ns, nt)], u, v);
                }
                at(i, j) = static_cast<int>(verts_.size());
                verts_.push_back(p);
            }
        }
        for (int j = 0; j < nt; ++j) {
            for (int i = 0; i < ns; ++i) {
                cells_.push_back({{at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)}});
                cellMaterial_.push_back(c.material);
            }
        }
        materials.insert(c.material);
        grids_.push_back(std::move(B));
    }
    report_.blocks = static_cast<int>(grids_.size());
    report_.materials = static_cast<int>(materials.size());
}

// ---------------------------------------------------------------------------
// smooth()
//
// Stage 10's Winslow pass, unchanged: Gauss-Seidel on the interior nodes of
// the blocks whose worst element is below the threshold, the rim held so the
// mesh stays conforming. See src/MERIDIAN/QuadMesh.hxx for why it is there.
// ---------------------------------------------------------------------------
void BlockMesh::smooth() {
    double before = std::numeric_limits<double>::infinity();
    for (const std::array<int, 4> &q : cells_) {
        const Point c[4] = {verts_[q[0]], verts_[q[1]], verts_[q[2]], verts_[q[3]]};
        const double sj = scaledJacobian(c);
        before = std::min(before, sj);
        if (!(sj > 0.0) || !(quadArea(c) > 0.0)) ++report_.invertedBefore;
    }
    report_.minScaledJacobianBefore = std::isfinite(before) ? before : 0.0;
    if (opts_.smoothingPasses <= 0 || cells_.empty()) return;

    const double tol = opts_.smoothingTolerance * opts_.targetEdgeLength;
    for (const Block &b : grids_) {
        if (b.ns < 2 || b.nt < 2) continue;
        const int w = b.ns + 1;
        auto id = [&](int i, int j) { return b.vert[static_cast<size_t>(j) * w + i]; };
        double worst = std::numeric_limits<double>::infinity();
        for (int j = 0; j < b.nt; ++j) {
            for (int i = 0; i < b.ns; ++i) {
                const Point c[4] = {verts_[id(i, j)], verts_[id(i + 1, j)], verts_[id(i + 1, j + 1)],
                                    verts_[id(i, j + 1)]};
                worst = std::min(worst, scaledJacobian(c));
            }
        }
        if (!(worst < opts_.smoothingThreshold)) continue;
        ++report_.smoothedBlocks;

        int sweep = 0;
        for (; sweep < opts_.smoothingPasses; ++sweep) {
            double moved = 0.0;
            for (int j = 1; j < b.nt; ++j) {
                for (int i = 1; i < b.ns; ++i) {
                    const Point &xe = verts_[id(i + 1, j)];
                    const Point &xw = verts_[id(i - 1, j)];
                    const Point &xn = verts_[id(i, j + 1)];
                    const Point &xs = verts_[id(i, j - 1)];
                    const Point ds = (xe - xw) * 0.5;
                    const Point dt = (xn - xs) * 0.5;
                    const double alpha = dotP(dt, dt);
                    const double beta = dotP(ds, dt);
                    const double gamma = dotP(ds, ds);
                    const double denom = 2.0 * (alpha + gamma);
                    if (!(denom > 0.0)) continue;
                    const Point cross = (verts_[id(i + 1, j + 1)] - verts_[id(i - 1, j + 1)] -
                                         verts_[id(i + 1, j - 1)] + verts_[id(i - 1, j - 1)]) *
                                        0.25;
                    const Point next = ((xe + xw) * alpha + (xn + xs) * gamma - cross * (2.0 * beta)) / denom;
                    Point &here = verts_[id(i, j)];
                    moved = std::max(moved, normP(next - here));
                    here = next;
                }
            }
            if (moved <= tol) { ++sweep; break; }
        }
        report_.smoothingSweeps = std::max(report_.smoothingSweeps, sweep);
    }
}

// ---------------------------------------------------------------------------
// check()
//
// Sec. 11.3's three checks on the finished mesh -- the corner Jacobians, the
// incidence, the boundary error -- and Stage 10's edge-length audit, all read
// back off the result.
// ---------------------------------------------------------------------------
void BlockMesh::check(const BlockCover &cover) {
    report_.vertices = static_cast<int>(verts_.size());
    report_.quads = static_cast<int>(cells_.size());
    if (cells_.empty()) {
        report_.messages.push_back("no block of the cover could be meshed");
        return;
    }

    // Edges: use count, lengths against the target, and the materials either
    // side.
    const long long n = static_cast<long long>(verts_.size());
    std::unordered_map<long long, std::array<int, 3>> edgeUse;  // count, first material, mixed
    edgeUse.reserve(cells_.size() * 4);
    double lenSum = 0.0, logSum = 0.0;
    int lenCount = 0;
    report_.minEdge = std::numeric_limits<double>::infinity();
    const double h = opts_.targetEdgeLength;
    for (size_t qi = 0; qi < cells_.size(); ++qi) {
        const std::array<int, 4> &q = cells_[qi];
        for (int k = 0; k < 4; ++k) {
            const int a = q[k], b = q[(k + 1) % 4];
            const long long key = static_cast<long long>(std::min(a, b)) * n + std::max(a, b);
            auto it = edgeUse.find(key);
            if (it == edgeUse.end()) {
                edgeUse.emplace(key, std::array<int, 3>{{1, cellMaterial_[qi], 0}});
                const double L = normP(verts_[b] - verts_[a]);
                report_.minEdge = std::min(report_.minEdge, L);
                report_.maxEdge = std::max(report_.maxEdge, L);
                lenSum += L;
                ++lenCount;
                if (L > 0.0) {
                    const double r = std::log(L / h);
                    logSum += r * r;
                    if (std::fabs(r) > std::fabs(std::log(report_.worstEdgeRatio))) report_.worstEdgeRatio = L / h;
                }
            } else {
                ++it->second[0];
                if (it->second[1] != cellMaterial_[qi]) it->second[2] = 1;
            }
        }
    }
    if (!std::isfinite(report_.minEdge)) report_.minEdge = 0.0;
    report_.meanEdge = lenCount > 0 ? lenSum / lenCount : 0.0;
    report_.edgeRatioRms = lenCount > 0 ? std::sqrt(logSum / lenCount) : 0.0;
    for (const auto &e : edgeUse) {
        if (e.second[0] == 1) ++report_.boundaryEdges;
        else if (e.second[0] == 2) ++report_.interiorEdges;
        else ++report_.nonManifoldEdges;
        if (e.second[2]) ++report_.interfaceEdges;
    }

    // Elements.
    double sjSum = 0.0;
    report_.minScaledJacobian = std::numeric_limits<double>::infinity();
    report_.minQuadArea = std::numeric_limits<double>::infinity();
    for (const std::array<int, 4> &q : cells_) {
        const Point c[4] = {verts_[q[0]], verts_[q[1]], verts_[q[2]], verts_[q[3]]};
        const double area = quadArea(c);
        const double sj = scaledJacobian(c);
        report_.meshArea += area;
        report_.minQuadArea = std::min(report_.minQuadArea, area);
        report_.maxQuadArea = std::max(report_.maxQuadArea, area);
        report_.minScaledJacobian = std::min(report_.minScaledJacobian, sj);
        sjSum += sj;
        // Sec. 9.1: the bilinear map is positive iff all four corner
        // determinants are, so a non-convex element with one corner turned
        // over is invalid even though its area is positive. Stage 10 counts
        // area alone; Sec. 11.3 asks for the corners.
        if (!(sj > 0.0) || !(area > 0.0)) ++report_.invertedQuads;
    }
    if (!std::isfinite(report_.minScaledJacobian)) report_.minScaledJacobian = 0.0;
    if (!std::isfinite(report_.minQuadArea)) report_.minQuadArea = 0.0;
    report_.meanScaledJacobian = sjSum / cells_.size();

    // Cracks: two distinct vertices at one point, on a hash grid.
    const double tol = opts_.crackTolerance * report_.modelExtent;
    const double cell = std::max(tol * 4.0, report_.modelExtent * 1e-12);
    std::unordered_map<long long, std::vector<int>> grid;
    auto key = [&](double x, double y) {
        return static_cast<long long>(std::floor(x / cell)) * 73856093LL ^
               static_cast<long long>(std::floor(y / cell)) * 19349663LL;
    };
    for (size_t i = 0; i < verts_.size(); ++i) grid[key(verts_[i][0], verts_[i][1])].push_back(static_cast<int>(i));
    for (size_t i = 0; i < verts_.size(); ++i) {
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                auto it = grid.find(key(verts_[i][0] + dx * cell, verts_[i][1] + dy * cell));
                if (it == grid.end()) continue;
                for (int j : it->second) {
                    if (j <= static_cast<int>(i)) continue;
                    if (normP(verts_[j] - verts_[i]) <= tol) ++report_.cracks;
                }
            }
        }
    }

    // Block corners the mesh turns through more than a half turn at, as Stage
    // 10 counts them: the element there is reversed whatever the interior does.
    for (const Block &b : grids_) {
        const int w = b.ns + 1;
        auto id = [&](int i, int j) { return b.vert[static_cast<size_t>(j) * w + i]; };
        const int ci[4] = {0, b.ns, b.ns, 0}, cj[4] = {0, 0, b.nt, b.nt};
        const int ai[4] = {1, b.ns, b.ns - 1, 0}, aj[4] = {0, 1, b.nt, b.nt - 1};
        const int bi[4] = {0, b.ns - 1, b.ns, 1}, bj[4] = {1, 0, b.nt - 1, b.nt};
        for (int k = 0; k < 4; ++k) {
            const Point o = verts_[id(ci[k], cj[k])];
            const Point u = verts_[id(ai[k], aj[k])] - o;
            const Point v = verts_[id(bi[k], bj[k])] - o;
            double a = std::atan2(cross2(u, v), dotP(u, v));
            if (a < 0.0) a += 2.0 * M_PI;
            if (a > M_PI) ++report_.reflexCorners;
        }
    }

    // Sec. 11.3's boundary error: how far each carrier vertex of a boundary
    // (or interface) macro edge lies from the chord between the two nodes
    // either side of it. Nodes sit at equal arc length, so the chord a chain
    // vertex falls under is found from its own arc length.
    const std::vector<BlockCover::MacroEdge> &edges = cover.getMacroEdges();
    for (size_t a = 0; a < edges.size(); ++a) {
        const BlockCover::MacroEdge &me = edges[a];
        if (!me.boundary && !me.interface) continue;
        const std::vector<int> &nodes = arcNodes_[a];
        const int N = static_cast<int>(nodes.size()) - 1;
        if (N < 1) continue;
        if (me.boundary) report_.boundaryNodes += N;
        const std::vector<double> d = chainLengths(C_, me.chain);
        const double L = d.back();
        double worst = 0.0;
        for (size_t m = 1; m + 1 < me.chain.size(); ++m) {
            const int k = L > 0.0 ? std::min(N - 1, static_cast<int>(std::floor(d[m] / L * N))) : 0;
            worst = std::max(worst, segmentDistance(C_.vertices[me.chain[m]], verts_[nodes[k]], verts_[nodes[k + 1]]));
        }
        if (me.boundary) report_.boundaryDeviation = std::max(report_.boundaryDeviation, worst);
        else report_.interfaceDeviation = std::max(report_.interfaceDeviation, worst);
    }

    report_.conforming = report_.nonManifoldEdges == 0 && report_.cracks == 0;
    report_.valid = report_.conforming && report_.invertedQuads == 0 && report_.unmeshedBlocks == 0;

    if (report_.unmeshedBlocks > 0) {
        report_.messages.push_back(std::to_string(report_.unmeshedBlocks) +
                                   " block(s) have a side that is not one macro edge and were left unmeshed; "
                                   "the cover was not conforming (Sec. 13.2)");
    }
    if (report_.invertedQuads > 0) {
        report_.messages.push_back(
            std::to_string(report_.invertedQuads) + " element(s) came out with a corner Jacobian that is not positive" +
            (report_.reflexCorners > 0
                 ? ", against " + std::to_string(report_.reflexCorners) +
                       " block corner(s) the mesh turns through more than a half turn at"
                 : std::string("; the chart is valid, so these are straight elements spanning a curved piece "
                               "of it (Sec. 11.3)")));
    }
    if (report_.nonManifoldEdges > 0) {
        report_.messages.push_back(std::to_string(report_.nonManifoldEdges) +
                                   " edge(s) are used by more than two elements");
    }
    if (report_.cracks > 0) {
        report_.messages.push_back(std::to_string(report_.cracks) + " pair(s) of distinct vertices are at the same point");
    }
}

// ---------------------------------------------------------------------------
bool BlockMesh::writeOBJ(const std::string &path) const {
    std::ofstream out(path);
    if (!out) return false;
    out.precision(17);
    out << "# ATLAS BlockMesh: " << cells_.size() << " quads over " << grids_.size() << " blocks\n";
    for (const Point &p : verts_) out << "v " << p[0] << " " << p[1] << " 0\n";
    int current = std::numeric_limits<int>::min();
    for (size_t q = 0; q < cells_.size(); ++q) {
        if (cellMaterial_[q] != current) {
            current = cellMaterial_[q];
            out << "usemtl mat" << current << "\n";
        }
        out << "f";
        for (int v : cells_[q]) out << " " << v + 1;
        out << "\n";
    }
    return static_cast<bool>(out);
}
