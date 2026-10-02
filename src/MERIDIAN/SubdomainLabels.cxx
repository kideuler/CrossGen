#include "SubdomainLabels.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <sstream>
#include <unordered_set>

namespace {

inline int localOf(const Triangle &t, int v) {
    if (t[0] == v) return 0;
    if (t[1] == v) return 1;
    if (t[2] == v) return 2;
    return -1;
}

} // namespace

// ---------------------------------------------------------------------------
SubdomainLabels::SubdomainLabels(const Immersion &immersion)
    : SubdomainLabels(immersion, Options()) {}

SubdomainLabels::SubdomainLabels(const Immersion &immersion, const Options &opts)
    : SubdomainLabels(immersion, opts, nullptr) {}

SubdomainLabels::SubdomainLabels(const Immersion &immersion, const Options &opts,
                                 const Interfaces *interfaces)
    : imm(&immersion), itf(interfaces), options(opts) {
    if (itf && !itf->multiMaterial()) itf = nullptr;
    buildSeamTables();

    report.seamArcs = static_cast<int>(imm->getArcs().size());
    for (const Immersion::Arc &a : imm->getArcs()) {
        if (!a.degenerate) ++report.holonomyCount[((a.k % 4) + 4) % 4];
    }

    const std::vector<Point> &uv = imm->getUV();
    buildBoundary(uv);
    buildBoundaryChains(uv);
    buildFeatures(uv);
    computeSeamTurns();
    buildInterfaceCorners();
    labelFeatures(uv);
    if (options.seedTopoConstraints) seedTopoPaths(uv);
}

// ---------------------------------------------------------------------------
// buildSeamTables()
//
// Three lookups that everything below wants: Omega's edges by their endpoints,
// which Omega boundary edges came from dS rather than from a cut, and, for the
// ones that came from a cut, which seam pair they belong to and on which side.
// ---------------------------------------------------------------------------
void SubdomainLabels::buildSeamTables() {
    const Mesh &cm = imm->getCutMesh();
    const Mesh &om = imm->getOriginalMesh();
    const std::vector<int> &c2o = imm->getCut().getCutVertexToOriginal();

    cutEdgeIndex.reserve(cm.edges.size() * 2);
    for (int e = 0; e < static_cast<int>(cm.edges.size()); ++e) {
        cutEdgeIndex.emplace(EdgeKey(cm.edges[e][0], cm.edges[e][1]), e);
    }

    std::unordered_map<EdgeKey, int, EdgeKeyHash> origEdgeIndex;
    origEdgeIndex.reserve(om.edges.size() * 2);
    for (int e = 0; e < static_cast<int>(om.edges.size()); ++e) {
        origEdgeIndex.emplace(EdgeKey(om.edges[e][0], om.edges[e][1]), e);
    }

    parentOnBoundary.assign(cm.edges.size(), 0);
    for (int e = 0; e < static_cast<int>(cm.edges.size()); ++e) {
        auto it = origEdgeIndex.find(EdgeKey(c2o[cm.edges[e][0]], c2o[cm.edges[e][1]]));
        if (it != origEdgeIndex.end() && om.isBoundaryEdge[it->second]) {
            parentOnBoundary[e] = 1;
        }
    }

    const auto &pairs = imm->getSeamPairs();
    seamSide.reserve(pairs.size() * 4);
    for (int i = 0; i < static_cast<int>(pairs.size()); ++i) {
        seamSide[EdgeKey(pairs[i][0], pairs[i][1])] = 2 * i + 0;   // the plus child
        seamSide[EdgeKey(pairs[i][2], pairs[i][3])] = 2 * i + 1;   // the minus child
    }
}

// ---------------------------------------------------------------------------
// buildBoundary()  --  Sec. 3.3
//
//     if |u(v_i) - u(v_j)| <= |v(v_i) - v(v_j)| :  e in Gamma_u
//     else                                      :  e in Gamma_v
// ---------------------------------------------------------------------------
void SubdomainLabels::buildBoundary(const std::vector<Point> &uv) {
    const Mesh &cm = imm->getCutMesh();
    const std::vector<double> &len = imm->getCutEdgeLengths();

    bEdges.clear();
    report.boundaryEdges = report.boundaryEdgesU = report.boundaryEdgesV = 0;
    report.ambiguousBoundaryEdges = 0;

    for (int e : cm.boundaryEdges) {
        if (!parentOnBoundary[e]) continue;   // a seam child, not a curve of dS

        BoundaryEdge be;
        be.cutEdge = e;
        be.a = cm.edges[e][0];
        be.b = cm.edges[e][1];
        be.length = len[e];

        const double du = std::fabs(uv[be.a][0] - uv[be.b][0]);
        const double dv = std::fabs(uv[be.a][1] - uv[be.b][1]);
        be.label = (du <= dv) ? Align::U : Align::V;

        const double big = std::max(du, dv);
        if (big > 0.0 && std::min(du, dv) > 0.8 * big) ++report.ambiguousBoundaryEdges;

        (be.label == Align::U ? report.boundaryEdgesU : report.boundaryEdgesV) += 1;
        bEdges.push_back(be);
    }
    report.boundaryEdges = static_cast<int>(bEdges.size());
}

// ---------------------------------------------------------------------------
// buildBoundaryChains()
//
// Q3 is a statement about the components of dS - G, so it is worth having them
// as objects even though E2 never needs them: a chain that comes out with two
// labels is a chain whose cone set is wrong, and that is far easier to see here
// than in a residual. A chain ends at a cone (a corner of the layout), where
// the cutting graph meets the boundary, at a junction of more than two boundary
// edges, and wherever the label changes.
// ---------------------------------------------------------------------------
void SubdomainLabels::buildBoundaryChains(const std::vector<Point> &uv) {
    bChains.clear();
    if (bEdges.empty()) { report.boundaryChains = 0; return; }

    const std::vector<int> &vCone = imm->getCutVertexCone();

    std::unordered_map<int, std::vector<int>> at;   // Omega vertex -> bEdge indices
    at.reserve(bEdges.size() * 2);
    for (int i = 0; i < static_cast<int>(bEdges.size()); ++i) {
        at[bEdges[i].a].push_back(i);
        at[bEdges[i].b].push_back(i);
    }

    // A chain cannot run past a cone (that is a corner of the layout) nor past
    // a point where the cutting graph meets dS. The third way a chain ends --
    // the label changing -- is tested in the walk itself, where both edges are
    // to hand.
    auto breaks = [&](int v) -> bool {
        return vCone[v] >= 0 || at[v].size() != 2;
    };

    std::vector<char> used(bEdges.size(), 0);
    for (int seed = 0; seed < static_cast<int>(bEdges.size()); ++seed) {
        if (used[seed]) continue;

        std::vector<int> chain{seed};
        used[seed] = 1;

        // Grow in both directions from the seed edge.
        for (int side = 0; side < 2; ++side) {
            int cur = seed;
            int v = (side == 0) ? bEdges[seed].a : bEdges[seed].b;
            while (true) {
                if (breaks(v)) break;
                int next = -1;
                for (int c : at[v]) if (c != cur && !used[c]) next = c;
                if (next < 0) break;
                if (bEdges[next].label != bEdges[cur].label) break;
                used[next] = 1;
                if (side == 0) chain.insert(chain.begin(), next);
                else chain.push_back(next);
                cur = next;
                v = (bEdges[next].a == v) ? bEdges[next].b : bEdges[next].a;
            }
        }

        BoundaryChain bc;
        bc.label = bEdges[chain.front()].label;
        double lo0 = std::numeric_limits<double>::infinity(), hi0 = -lo0;
        double lo1 = lo0, hi1 = hi0;
        for (int i : chain) {
            bEdges[i].chain = static_cast<int>(bChains.size());
            bc.length += bEdges[i].length;
            for (int v : {bEdges[i].a, bEdges[i].b}) {
                lo0 = std::min(lo0, uv[v][0]); hi0 = std::max(hi0, uv[v][0]);
                lo1 = std::min(lo1, uv[v][1]); hi1 = std::max(hi1, uv[v][1]);
            }
        }
        bc.spanU = hi0 - lo0;
        bc.spanV = hi1 - lo1;
        bc.edges = std::move(chain);
        bChains.push_back(std::move(bc));
    }
    report.boundaryChains = static_cast<int>(bChains.size());
}

// ---------------------------------------------------------------------------
// buildFeatures()
//
// Stage 0 of the paper marks trim curves, dihedral creases and user-designated
// curves as features. These models are planar midsurfaces with material tags,
// so the analogue here is the interface between two physical surfaces -- the
// same edges the viewer already draws as feature edges. E3 then guarantees what
// Sec. 3.3 promises: a feature preserved from the input is a layout edge in the
// output, not something a patch runs across.
//
// The chains are built in Omega rather than in S, so a feature that happens to
// run along an arc of the cutting graph contributes both of its children and
// each gets its own constraint -- which is what E3 needs, since the two sides
// are different unknowns.
// ---------------------------------------------------------------------------
void SubdomainLabels::buildFeatures(const std::vector<Point> &uv) {
    fChains.clear();
    report.featureEdges = 0;
    report.featureChains = report.featureChainsU = report.featureChainsV = 0;
    if (!options.materialFeatures) return;

    const Mesh &om = imm->getOriginalMesh();
    const Mesh &cm = imm->getCutMesh();
    if (om.triangleMatId.size() != om.triangles.size()) return;

    // The feature edges of Omega, and which branch of the network each came
    // from. An interface edge contributes both of its Omega children when it
    // runs along an arc of G, and both inherit the branch.
    std::unordered_set<int> featEdges;
    std::unordered_map<int, int> edgeBranch;   // Omega edge -> branch, or -1
    for (int e = 0; e < static_cast<int>(om.edges.size()); ++e) {
        const int f0 = om.edgeTriangles[e][0];
        const int f1 = om.edgeTriangles[e][1];
        if (f0 < 0 || f1 < 0) continue;                          // dS, handled by E2
        if (om.triangleMatId[f0] == om.triangleMatId[f1]) continue;

        const int a = om.edges[e][0], b = om.edges[e][1];
        const int br = itf ? itf->branchOfEdge(e) : -1;
        for (int f : {f0, f1}) {
            const int la = localOf(om.triangles[f], a);
            const int lb = localOf(om.triangles[f], b);
            if (la < 0 || lb < 0) continue;
            auto it = cutEdgeIndex.find(EdgeKey(cm.triangles[f][la], cm.triangles[f][lb]));
            if (it != cutEdgeIndex.end()) {
                featEdges.insert(it->second);
                edgeBranch[it->second] = br;
            }
        }
    }
    report.featureEdges = static_cast<int>(featEdges.size());
    if (featEdges.empty()) return;

    // Where a chain has to end because the layout turns there: the nodes of the
    // interface network, and the cones, which are the two kinds of point a
    // layout edge is allowed to stop at. Both are read in Omega, so a node the
    // cutting graph split into several children stops every one of them.
    const std::vector<int> &c2o = imm->getCut().getCutVertexToOriginal();
    const std::vector<int> &vertCone = imm->getCutVertexCone();
    std::vector<char> stopVert(cm.vertices.size(), 0);
    for (size_t v = 0; v < cm.vertices.size(); ++v) {
        if (v < vertCone.size() && vertCone[v] >= 0) stopVert[v] = 1;
        if (itf && v < c2o.size() && c2o[v] >= 0 && itf->nodeAt(c2o[v]) >= 0) stopVert[v] = 1;
    }

    // Chain them: nodes where the feature graph is not a simple path, and
    // corners where the image direction turns, both end a chain.
    std::unordered_map<int, std::vector<int>> at;
    for (int e : featEdges) {
        at[cm.edges[e][0]].push_back(e);
        at[cm.edges[e][1]].push_back(e);
    }

    const std::vector<double> &len = imm->getCutEdgeLengths();
    std::unordered_set<int> used;

    // The signed turn of the chain at v, going in through e0 and out through e1.
    auto turnAt = [&](int v, int e0, int e1) -> double {
        auto dir = [&](int e, int from) {
            const int other = (cm.edges[e][0] == from) ? cm.edges[e][1] : cm.edges[e][0];
            return normalizeP(uv[other] - uv[from]);
        };
        const Point in = dir(e0, v) * -1.0;   // the direction of travel arriving at v
        const Point out = dir(e1, v);
        return std::atan2(cross2(in, out), dotP(in, out));
    };

    for (int seed : featEdges) {
        if (used.count(seed)) continue;
        std::vector<int> chainEdges{seed};
        used.insert(seed);

        // A chain is grown until it has turned too far, not merely until it
        // turns sharply at one vertex. Both tests matter and they catch
        // different things. A single sharp vertex is a corner of the feature
        // and obviously ends a chain. But a feature that curves *smoothly*
        // through a right angle -- a filleted material interface, which is
        // most of them -- has no sharp vertex anywhere, and taking it as one
        // chain asks E3 to lay a quarter circle on a single isoline. There is
        // no such map at any finite distortion: the barrier of E1 and the
        // penalty of E3 come to an equilibrium with det J at 1e-4 and the
        // constraint residual stuck. Splitting on the *accumulated* turn
        // instead makes each piece a curve that a coordinate line can plausibly
        // follow, which is what the paper's own feature chains -- creases and
        // trim curves between cone-marked corners -- already are.
        const int seedBranch = edgeBranch.count(seed) ? edgeBranch[seed] : -1;
        for (int side = 0; side < 2; ++side) {
            int cur = seed;
            int v = (side == 0) ? cm.edges[seed][0] : cm.edges[seed][1];
            double turned = 0.0;
            while (at[v].size() == 2) {
                if (stopVert[v]) break;
                int next = -1;
                for (int c : at[v]) if (c != cur) next = c;
                if (next < 0 || used.count(next)) break;
                // With a network in hand a chain is a *branch*, and it runs
                // from node to node however much it curves on the way. The
                // accumulated-turn split below is what has to be done without
                // one: it is a guess at where the corners are, and on an
                // interface it guesses wrong in the expensive direction,
                // cutting a smooth quarter-circle into two chains and asking
                // E3 to hold each of them on a different isoline, which puts a
                // corner of the layout in the middle of a smooth interface.
                if (!itf) {
                    const double t = turnAt(v, cur, next);
                    if (std::fabs(t) > options.featureCornerAngle) break;
                    turned += t;
                    if (std::fabs(turned) > options.featureCornerAngle) break;
                } else {
                    if (edgeBranch.count(next) && edgeBranch[next] != seedBranch) break;
                }
                used.insert(next);
                if (side == 0) chainEdges.insert(chainEdges.begin(), next);
                else chainEdges.push_back(next);
                cur = next;
                v = (cm.edges[next][0] == v) ? cm.edges[next][1] : cm.edges[next][0];
            }
        }

        // Turn the edge run into an ordered vertex run.
        FeatureChain fc;
        {
            int prev = -1;
            const int e0 = chainEdges.front();
            if (chainEdges.size() == 1) {
                prev = cm.edges[e0][0];
            } else {
                const int e1 = chainEdges[1];
                const int s0 = cm.edges[e0][0], s1 = cm.edges[e0][1];
                const bool s1Shared = (cm.edges[e1][0] == s1 || cm.edges[e1][1] == s1);
                prev = s1Shared ? s0 : s1;
            }
            fc.verts.push_back(prev);
            for (int e : chainEdges) {
                const int nxt = (cm.edges[e][0] == prev) ? cm.edges[e][1] : cm.edges[e][0];
                fc.verts.push_back(nxt);
                fc.lengths.push_back(len[e]);
                fc.length += len[e];
                if (seamSide.count(EdgeKey(cm.edges[e][0], cm.edges[e][1]))) fc.onSeam = true;
                prev = nxt;
            }
        }
        fc.branch = seedBranch;
        fChains.push_back(std::move(fc));
    }
    report.featureChains = static_cast<int>(fChains.size());
}

// ---------------------------------------------------------------------------
// buildInterfaceCorners()   --  what E6 is summed over
//
// One entry per sector of the interface network: the two curves that bound it
// and the whole number of quarter turns the layout has to put between them,
// which Stage 0b already worked out from the geometry.
//
// The two tangents are taken at the same child of the node in Omega -- the
// faces either side of each ray say which child, and taking the face on the
// side the sector is on picks the right one whenever the cutting graph does not
// run through the sector itself. Where it does there is no common child and the
// sector is dropped: a tangent at one child and a tangent at another are not in
// the same frame, and forcing a quarter turn between them would be asking for
// the transition Q4 already holds.
//
// A node on dS contributes its two boundary edges as rays as well. That is not
// bookkeeping: an interface meeting the boundary has to meet it at a whole
// number of quarter turns, and without the boundary rays nothing says so --
// E2 puts dS on a coordinate line and E3 puts the interface on one, but which
// one, and which way round, is exactly what the sector count fixes.
// ---------------------------------------------------------------------------
void SubdomainLabels::buildInterfaceCorners() {
    iCorners.clear();
    report.interfaceCorners = report.interfaceCornersSpanningCut = 0;
    if (!itf || !options.interfaceCorners) return;

    const Mesh &om = imm->getOriginalMesh();
    const Mesh &cm = imm->getCutMesh();
    const std::vector<double> &len = imm->getCutEdgeLengths();

    // The child of original vertex `v` seen from face `f`.
    auto childIn = [&](int v, int f) -> int {
        if (f < 0) return -1;
        const int l = localOf(om.triangles[f], v);
        return (l < 0) ? -1 : cm.triangles[f][l];
    };
    auto flatLength = [&](int a, int b) -> double {
        auto it = cutEdgeIndex.find(EdgeKey(a, b));
        return (it == cutEdgeIndex.end()) ? 0.0 : len[it->second];
    };

    const std::vector<Interfaces::Node> &nds = itf->nodes();
    for (size_t n = 0; n < nds.size(); ++n) {
        const Interfaces::Node &nd = nds[n];
        if (nd.rays.size() < 2 || nd.quarters.size() + 1 < nd.rays.size()) continue;

        // Consecutive pairs only, never the wrap-around sector of an interior
        // node. Composing every sector of a node closes the loop with
        // R_(sum q) = R_(4 - I), which is the identity only when the node
        // carries no cone; at one that does, adding the last sector would be
        // asserting that the cone is not there.
        for (size_t k = 0; k + 1 < nd.rays.size(); ++k) {
            const Interfaces::Ray &ra = nd.rays[k];
            const Interfaces::Ray &rb = nd.rays[k + 1];
            const int fa = ra.faceCCW;      // the face on the sector's side of a
            const int fb = rb.faceCW;       // ... and on the sector's side of b
            const int va = childIn(nd.vertex, fa);
            const int vb = childIn(nd.vertex, fb);
            const int wa = childIn(ra.neighbour, fa);
            const int wb = childIn(rb.neighbour, fb);
            if (va < 0 || vb < 0 || wa < 0 || wb < 0) continue;
            if (va != vb) { ++report.interfaceCornersSpanningCut; continue; }

            InterfaceCorner ic;
            ic.node = static_cast<int>(n);
            ic.vertex = va;
            ic.a = wa;
            ic.b = wb;
            ic.la = flatLength(va, wa);
            ic.lb = flatLength(vb, wb);
            ic.sa = ic.la;
            ic.sb = ic.lb;
            ic.quarters = nd.quarters[k];
            ic.onBoundary = nd.onBoundary;
            ic.aIsBoundary = (ra.branch < 0);
            ic.bIsBoundary = (rb.branch < 0);
            ic.sector = nd.sector[k];
            if (!(ic.la > 0.0) || !(ic.lb > 0.0)) continue;
            iCorners.push_back(ic);
        }
    }
    report.interfaceCorners = static_cast<int>(iCorners.size());
}

// ---------------------------------------------------------------------------
// computeSeamTurns()
//
// A branch the cutting graph crosses is not in one chart. The crossing point of
// S has two children in Omega, the edges of the branch before it hang off one
// and the edges after it off the other -- which is why buildFeatures() ends a
// chain there -- and the two sides are related by Stage 4's transition
// psi+ = R_k psi- + t. A direction d on the minus side is d + k on the plus
// side, the same rule addTopoPath() applies to a path of Gamma_topo.
//
// This is not a corner case of multi-material models, it is what every
// *enclosed* region looks like: a cone inside an inclusion has no route to dS
// that does not cross the interface around it (ConeCut can only choose where),
// so every such interface comes in two or more charts. The rocket's
// casing/body interface is a U round two +1 cones, cut twice, each cut a
// quarter turn: its middle chain holds the other coordinate from its two ends.
//
// Branches are walked edge by edge from node0 against Interfaces::Branch::
// edges, so the order is the model's and not a guess from chord directions.
// A branch the walk cannot account for -- one running along a seam, which
// gives two chains over the same edges, or one stopped by a cone, where the
// turn is the cone's sector and not a transition -- keeps turn 0 on every
// chain, which is the behaviour before this existed.
// ---------------------------------------------------------------------------
void SubdomainLabels::computeSeamTurns() {
    report.featureChainsPastSeam = report.featureChainsSeamFlipped = 0;
    report.branchesUnwalked = 0;
    branchTurn.clear();
    for (FeatureChain &fc : fChains) { fc.turn = 0; fc.turnKnown = false; fc.sense = 0; }
    if (!itf) return;
    const std::vector<Interfaces::Branch> &brs = itf->branches();
    branchTurn.assign(brs.size(), 0);
    if (!options.seamTurnInterfaceLabels || brs.empty()) return;

    // (child on one side, child on the other) -> +(arc + 1) plus to minus,
    // -(arc + 1) minus to plus.
    const auto &arcList = imm->getArcs();
    std::unordered_map<long long, int> crossing;
    for (int a = 0; a < static_cast<int>(arcList.size()); ++a) {
        const Immersion::Arc &arc = arcList[a];
        if (arc.degenerate) continue;
        const size_t n = std::min(arc.plusChain.size(), arc.minusChain.size());
        for (size_t i = 0; i < n; ++i) {
            if (arc.plusChain[i] == arc.minusChain[i]) continue;
            crossing.emplace(static_cast<long long>(arc.plusChain[i]) * 1000003LL + arc.minusChain[i],
                             +(a + 1));
            crossing.emplace(static_cast<long long>(arc.minusChain[i]) * 1000003LL + arc.plusChain[i],
                             -(a + 1));
        }
    }

    const std::vector<int> &c2o = imm->getCut().getCutVertexToOriginal();
    std::vector<std::vector<int>> chainsOf(brs.size());
    for (int i = 0; i < static_cast<int>(fChains.size()); ++i) {
        const int b = fChains[i].branch;
        if (b >= 0 && b < static_cast<int>(brs.size())) chainsOf[b].push_back(i);
    }

    struct Span { int chain; int lo, hi; int first, last; };   // Omega ends in branch order
    for (size_t b = 0; b < brs.size(); ++b) {
        const Interfaces::Branch &br = brs[b];
        if (chainsOf[b].empty() || br.verts.size() < 2) continue;
        const int nEdges = static_cast<int>(br.verts.size()) - 1;

        std::unordered_map<EdgeKey, int, EdgeKeyHash> at;
        at.reserve(br.verts.size() * 2);
        for (int i = 0; i < nEdges; ++i) at.emplace(EdgeKey(br.verts[i], br.verts[i + 1]), i);

        std::vector<Span> spans;
        bool ok = true;
        for (int ci : chainsOf[b]) {
            const FeatureChain &fc = fChains[ci];
            if (fc.verts.size() < 2) { ok = false; break; }
            int lo = nEdges, hi = -1, firstIdx = -1, lastIdx = -1;
            for (size_t j = 0; j + 1 < fc.verts.size() && ok; ++j) {
                const int o0 = c2o[fc.verts[j]], o1 = c2o[fc.verts[j + 1]];
                auto it = (o0 < 0 || o1 < 0) ? at.end() : at.find(EdgeKey(o0, o1));
                if (it == at.end()) { ok = false; break; }
                if (j == 0) firstIdx = it->second;
                lastIdx = it->second;
                lo = std::min(lo, it->second);
                hi = std::max(hi, it->second);
            }
            if (!ok || hi - lo + 2 != static_cast<int>(fc.verts.size())) { ok = false; break; }
            // Which end of `verts` sits at br.verts[lo]. With more than one
            // edge the order of the edge indices says; with one, the vertex.
            bool frontFirst = (firstIdx < lastIdx) ||
                              (firstIdx == lastIdx && c2o[fc.verts.front()] == br.verts[lo]);
            spans.push_back({ci, lo, hi, frontFirst ? fc.verts.front() : fc.verts.back(),
                             frontFirst ? fc.verts.back() : fc.verts.front()});
            fChains[ci].sense = frontFirst ? +1 : -1;
        }
        std::sort(spans.begin(), spans.end(),
                  [](const Span &x, const Span &y) { return x.lo < y.lo; });
        if (ok && (spans.front().lo != 0 || spans.back().hi != nEdges - 1)) ok = false;

        std::vector<int> turn(spans.size(), 0);
        for (size_t i = 1; ok && i < spans.size(); ++i) {
            // Contiguous and not overlapping: an overlap is a branch along a
            // seam, whose two children each made a chain over the same edges.
            if (spans[i].lo != spans[i - 1].hi + 1) { ok = false; break; }
            const int a = spans[i - 1].last, c = spans[i].first;
            if (a == c) { ok = false; break; }   // stopped by a cone, not a crossing
            auto it = crossing.find(static_cast<long long>(a) * 1000003LL + c);
            if (it == crossing.end()) { ok = false; break; }
            const int k = arcList[std::abs(it->second) - 1].k;
            turn[i] = (it->second > 0) ? turn[i - 1] - k : turn[i - 1] + k;
        }

        if (!ok) {
            for (int ci : chainsOf[b]) fChains[ci].sense = 0;
            ++report.branchesUnwalked;
            continue;
        }
        for (size_t i = 0; i < spans.size(); ++i) {
            FeatureChain &fc = fChains[spans[i].chain];
            fc.turn = ((turn[i] % 4) + 4) % 4;
            fc.turnKnown = true;
            if (i > 0) ++report.featureChainsPastSeam;
            if (fc.turn % 2) ++report.featureChainsSeamFlipped;
        }
        branchTurn[b] = ((turn.back() % 4) + 4) % 4;
    }
}

// ---------------------------------------------------------------------------
// propagateDirections()
//
// Which of {+u, +v, -u, -v} each branch of the network runs along.
//
// It is one propagation and not a vote, because the node quantisation already
// fixes every branch relative to every other one it can reach: ray k+1 leaves
// q_k quarter turns counter-clockwise of ray k, and the two ends of a branch
// point opposite ways. What is left free is one direction per connected
// component of the network, and that is the only thing read off the map.
//
// Doing it the other way -- a label per chain from its own flux, which is what
// Sec. 3.3 prescribes for an isolated feature -- is not merely less accurate.
// At the triple junction of geom003 the sectors are one, one and two quarters,
// so the three branches must be labelled u, v, u; psi_R's flux reads all three
// as u, and E3 is then asked for a map holding u constant on three curves
// leaving one point in three different directions, which has no solution at any
// penalty. The continuation does not diverge, it converges to a compromise in
// which none of the three is met.
// ---------------------------------------------------------------------------
std::vector<int> SubdomainLabels::propagateDirections(const std::vector<Point> &uv) {
    if (!itf) return {};
    const std::vector<Interfaces::Node> &nds = itf->nodes();
    const std::vector<Interfaces::Branch> &brs = itf->branches();
    if (brs.empty()) return {};

    // Global ray numbering, and the two rays of each branch.
    std::vector<int> rayBase(nds.size() + 1, 0);
    for (size_t n = 0; n < nds.size(); ++n) rayBase[n + 1] = rayBase[n] +
                                            static_cast<int>(nds[n].rays.size());
    const int nRays = rayBase.back();
    std::vector<std::vector<int>> branchRays(brs.size());
    for (size_t n = 0; n < nds.size(); ++n) {
        for (size_t k = 0; k < nds[n].rays.size(); ++k) {
            const int b = nds[n].rays[k].branch;
            if (b >= 0) branchRays[b].push_back(rayBase[n] + static_cast<int>(k));
        }
    }

    // The image displacement of each branch, node0 to node1, summed over the
    // chains that carry it, each rotated back into node0's chart by the seam
    // turns in front of it -- summed unrotated, the pieces of a branch cut
    // twice by a quarter turn cancel and the seed is read off noise. Where
    // computeSeamTurns() walked the branch the chain's orientation is known;
    // elsewhere it is compared against the branch's on the model.
    const Mesh &om = imm->getOriginalMesh();
    const std::vector<int> &c2o = imm->getCut().getCutVertexToOriginal();
    auto rotateQuarters = [](const Point &p, int q) -> Point {
        switch (((q % 4) + 4) % 4) {
            case 0: return p;
            case 1: return Point{-p[1], p[0]};
            case 2: return Point{-p[0], -p[1]};
            default: return Point{p[1], -p[0]};
        }
    };
    std::vector<Point> flux(brs.size(), Point{0.0, 0.0});
    std::vector<double> weight(brs.size(), 0.0);
    for (const FeatureChain &fc : fChains) {
        if (fc.branch < 0 || fc.branch >= static_cast<int>(brs.size()) || fc.verts.size() < 2) {
            continue;
        }
        double sign = fc.sense;
        if (sign == 0.0) {
            const Interfaces::Branch &br = brs[fc.branch];
            const Point mBranch = om.vertices[br.verts.back()] - om.vertices[br.verts.front()];
            const int o0 = c2o[fc.verts.front()], o1 = c2o[fc.verts.back()];
            if (o0 < 0 || o1 < 0) continue;
            const Point mChain = om.vertices[o1] - om.vertices[o0];
            sign = (dotP(mChain, mBranch) >= 0.0) ? 1.0 : -1.0;
        }
        const Point d = (uv[fc.verts.back()] - uv[fc.verts.front()]) * sign;
        flux[fc.branch] = flux[fc.branch] + rotateQuarters(d, -fc.turn);
        weight[fc.branch] += fc.length;
    }
    const auto turnOf = [&](int b) -> int {
        return (b >= 0 && b < static_cast<int>(branchTurn.size())) ? branchTurn[b] : 0;
    };
    // Whether ray k of node n is the node0 end of its branch. A branch that
    // leaves and returns to one node has both rays there, and only the one
    // pointing at the branch's second vertex is its start.
    auto atNode0 = [&](size_t n, size_t k) -> bool {
        const Interfaces::Ray &ray = nds[n].rays[k];
        if (ray.branch < 0) return false;
        const Interfaces::Branch &br = brs[ray.branch];
        if (static_cast<int>(n) != br.node0) return false;
        if (br.node0 != br.node1 || br.verts.size() < 2) return true;
        return ray.neighbour == br.verts[1];
    };

    std::vector<int> rayDir(nRays, -1);
    std::vector<char> seen(nRays, 0);

    // Seed order: the longest branch first, so each component is anchored on
    // the piece whose flux says the most.
    std::vector<int> order(brs.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = static_cast<int>(i);
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return weight[a] > weight[b]; });

    std::vector<int> stack;
    for (int b0 : order) {
        if (branchRays[b0].empty() || seen[branchRays[b0].front()]) continue;

        int best = 0;
        double bestDot = -std::numeric_limits<double>::infinity();
        for (int d = 0; d < 4; ++d) {
            const double v = dotP(flux[b0], axis(d));
            if (v > bestDot) { bestDot = v; best = d; }
        }
        // `best` is the direction of travel out of node0 in node0's chart, so
        // it belongs on the node0 ray. The node1 ray, where that is the only
        // one there is, points back along the branch in node1's chart.
        int seedRay = branchRays[b0].front();
        int seedDir = best + turnOf(b0) + 2;
        for (int r : branchRays[b0]) {
            size_t n = 0;
            while (n + 1 < nds.size() && rayBase[n + 1] <= r) ++n;
            if (atNode0(n, static_cast<size_t>(r - rayBase[n]))) { seedRay = r; seedDir = best; break; }
        }
        rayDir[seedRay] = ((seedDir % 4) + 4) % 4;
        stack.assign(1, seedRay);
        seen[seedRay] = 1;

        while (!stack.empty()) {
            const int r = stack.back();
            stack.pop_back();
            // Which node this ray belongs to.
            size_t n = 0;
            while (n + 1 < nds.size() && rayBase[n + 1] <= r) ++n;
            const int k = r - rayBase[n];
            const Interfaces::Node &nd = nds[n];

            auto push = [&](int to, int d) {
                if (to < 0 || seen[to]) return;
                seen[to] = 1;
                rayDir[to] = ((d % 4) + 4) % 4;
                stack.push_back(to);
            };
            // Along the node's fan: ray k+1 is q_k quarters counter-clockwise.
            if (k + 1 < static_cast<int>(nd.rays.size()) &&
                k < static_cast<int>(nd.quarters.size())) {
                push(r + 1, rayDir[r] + nd.quarters[k]);
            }
            if (k > 0 && k - 1 < static_cast<int>(nd.quarters.size())) {
                push(r - 1, rayDir[r] - nd.quarters[k - 1]);
            }
            // Across a branch: the far end points the other way, in a chart
            // turned by the seams the branch crosses on the way there.
            const int b = nd.rays[k].branch;
            if (b >= 0) {
                const int t = atNode0(n, static_cast<size_t>(k)) ? turnOf(b) : -turnOf(b);
                for (int other : branchRays[b]) if (other != r) push(other, rayDir[r] + t + 2);
            }
        }
    }

    // Back to one direction per branch: the way it leaves node0.
    branchDir.assign(brs.size(), -1);
    for (size_t b = 0; b < brs.size(); ++b) {
        if (branchRays[b].empty()) continue;
        // branchRays holds them in node order, and node0 is the smaller node
        // only by accident, so find the ray that actually sits at node0.
        for (int r : branchRays[b]) {
            size_t n = 0;
            while (n + 1 < nds.size() && rayBase[n + 1] <= r) ++n;
            if (atNode0(n, static_cast<size_t>(r - rayBase[n]))) { branchDir[b] = rayDir[r]; break; }
        }
        if (branchDir[b] < 0 && rayDir[branchRays[b].front()] >= 0) {
            // Only the node1 ray was reached: turn it back into node0's chart.
            branchDir[b] = (((rayDir[branchRays[b].front()] - turnOf(b) + 2) % 4) + 4) % 4;
        }
    }
    return branchDir;
}

// ---------------------------------------------------------------------------
// labelFeatures()
//
// One label per chain, from the *total* flux rather than an edge-by-edge vote.
// A chain that changed its mind halfway would ask E3 to hold u constant on one
// half and v on the other; there is no map that does both, so the continuation
// would grind against a constraint that cannot be met rather than converging.
// ---------------------------------------------------------------------------
void SubdomainLabels::labelFeatures(const std::vector<Point> &uv) {
    report.featureChainsU = report.featureChainsV = 0;
    report.featureChainsPropagated = 0;
    report.featureLabelsCorrected = 0;

    const bool propagate = itf && options.propagateInterfaceLabels;
    if (propagate) propagateDirections(uv);

    for (FeatureChain &fc : fChains) {
        if (fc.verts.size() < 2) { fc.label = Align::None; continue; }
        const Point d = uv[fc.verts.back()] - uv[fc.verts.front()];
        fc.fluxU = std::fabs(d[0]);
        fc.fluxV = std::fabs(d[1]);

        if (fc.fluxU == 0.0 && fc.fluxV == 0.0) {
            // A closed chain, or one that returned to where it started: the net
            // flux says nothing, so fall back to the total variation.
            double tvU = 0.0, tvV = 0.0;
            for (size_t i = 0; i + 1 < fc.verts.size(); ++i) {
                tvU += std::fabs(uv[fc.verts[i + 1]][0] - uv[fc.verts[i]][0]);
                tvV += std::fabs(uv[fc.verts[i + 1]][1] - uv[fc.verts[i]][1]);
            }
            fc.fluxU = tvU;
            fc.fluxV = tvV;
        }
        const Align byFlux = (fc.fluxU <= fc.fluxV) ? Align::U : Align::V;

        // The propagated direction wins where there is one. It is the direction
        // of *travel*, so the coordinate held constant is the other one: a
        // branch running along +-u holds v.
        fc.dirKnown = false;
        if (propagate && fc.branch >= 0 && fc.branch < static_cast<int>(branchDir.size()) &&
            branchDir[fc.branch] >= 0) {
            // branchDir is in node0's chart; the chain may be past a seam.
            fc.dir = (branchDir[fc.branch] + fc.turn) % 4;
            fc.dirKnown = true;
            fc.label = (fc.dir % 2 == 0) ? Align::V : Align::U;
            ++report.featureChainsPropagated;
            if (fc.label != byFlux) ++report.featureLabelsCorrected;
        } else {
            fc.label = byFlux;
        }
        (fc.label == Align::U ? report.featureChainsU : report.featureChainsV) += 1;
    }
}

// ---------------------------------------------------------------------------
// refreshInterfaceScales()
// ---------------------------------------------------------------------------
void SubdomainLabels::refreshInterfaceScales(const std::vector<Point> &uv) {
    for (InterfaceCorner &ic : iCorners) {
        if (ic.vertex < 0 || ic.a < 0 || ic.b < 0) continue;
        if (ic.vertex >= static_cast<int>(uv.size()) || ic.a >= static_cast<int>(uv.size()) ||
            ic.b >= static_cast<int>(uv.size())) {
            continue;
        }
        const double na = normP(uv[ic.a] - uv[ic.vertex]);
        const double nb = normP(uv[ic.b] - uv[ic.vertex]);
        // A tangent the map has collapsed says nothing about a direction, so
        // the flat length stands in rather than a division by nothing.
        ic.sa = (na > 0.0) ? na : ic.la;
        ic.sb = (nb > 0.0) ? nb : ic.lb;
    }
}

// ---------------------------------------------------------------------------
// relabel()
// ---------------------------------------------------------------------------
void SubdomainLabels::relabel(const std::vector<Point> &uv) {
    buildBoundary(uv);
    buildBoundaryChains(uv);
    labelFeatures(uv);
}

// ---------------------------------------------------------------------------
// meanConeSpacing()
//
// The scale a "near miss" is measured against: the mean distance from each cone
// to its nearest other cone, in the image. A tolerance in absolute units means
// nothing here -- the Ricci metric is only defined up to a scale -- and one as
// a fraction of the whole image means nothing either, because what makes two
// cones candidates for being joined is how they stand relative to the other
// cones, not to the model's bounding box.
// ---------------------------------------------------------------------------
double SubdomainLabels::meanConeSpacing(const Immersion &imm, const std::vector<Point> &uv) {
    const auto &children = imm.getConeChildren();
    double acc = 0.0;
    int counted = 0;
    for (size_t i = 0; i < children.size(); ++i) {
        if (children[i].empty()) continue;
        double best = std::numeric_limits<double>::infinity();
        for (size_t j = 0; j < children.size(); ++j) {
            if (i == j || children[j].empty()) continue;
            for (int a : children[i]) {
                for (int b : children[j]) best = std::min(best, normP(uv[a] - uv[b]));
            }
        }
        if (std::isfinite(best)) { acc += best; ++counted; }
    }
    return counted ? acc / counted : 0.0;
}

// ---------------------------------------------------------------------------
// seedTracerOptions()
//
// Stage 7's tracer, set up to look for near misses instead of terminations.
//
// Two of the settings have to move together and the reason is worth stating,
// because getting one of them right and not the other is a silent failure. The
// snap tolerance is widened from Sec. 3.4's 1e-6 of the extent to the near-miss
// tolerance -- that is the point of the seeding. But a cone is only ever
// offered to a curve if it is within Options::coneSnapRings faces of the
// triangle being crossed, and that restriction exists because Psi overlaps
// itself, so it cannot simply be dropped. Widen the tolerance without widening
// the neighbourhood and the seeder quietly stops finding the misses it was
// widened to find.
//
// So the radius is derived from the tolerance: an image distance of `tol` is,
// since psi_R is an isometry of the flat metric, about tol / (mean edge) faces
// on S. Clamped at both ends -- never fewer than two rings, never more than
// thirty-two, which is where the cost of the breadth-first walk starts to show.
// ---------------------------------------------------------------------------
Separatrices::Options SubdomainLabels::seedTracerOptions(const std::vector<Point> &uv,
                                                          double spacing) const {
    const Mesh &cm = imm->getCutMesh();

    Point lo{std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::infinity()};
    Point hi{-lo[0], -lo[1]};
    for (const Point &p : uv) {
        lo[0] = std::min(lo[0], p[0]); hi[0] = std::max(hi[0], p[0]);
        lo[1] = std::min(lo[1], p[1]); hi[1] = std::max(hi[1], p[1]);
    }
    double extent = uv.empty() ? 1.0 : std::hypot(hi[0] - lo[0], hi[1] - lo[1]);
    if (!(extent > 0.0)) extent = 1.0;

    double edgeAcc = 0.0;
    int edgeCount = 0;
    for (const auto &e : cm.edges) {
        edgeAcc += normP(uv[e[0]] - uv[e[1]]);
        ++edgeCount;
    }
    const double meanEdge = edgeCount ? edgeAcc / edgeCount : extent;

    const double tol = options.nearMissTolerance * spacing;

    Separatrices::Options so;
    // The same emitter set Stage 7 will use, so that the curves the seeding
    // sees and the curves Q5 is checked on are the same curves. A seeder with
    // its own idea of which points emit finds a different set, and the
    // difference shows up two stages later as a separatrix with no constraint
    // behind it.
    if (itf) so.extraEmitters = itf->emitterNodes();
    so.extraEmitters.insert(so.extraEmitters.end(), options.extraEmitters.begin(),
                            options.extraEmitters.end());
    so.coneSnapTolerance = tol / extent;
    so.coneSnapRings = (meanEdge > 0.0)
        ? std::max(2, std::min(32, static_cast<int>(std::ceil(tol / meanEdge)) + 1))
        : 2;
    // The separation cap is what stops the widened tolerance handing a curve to
    // the wrong member of a clustered pair. It is not optional here the way it
    // is at Sec. 3.4's tolerance, where nothing is near enough to be confused.
    so.coneSnapSeparationCap = 0.25;
    so.maxSteps = options.maxTraceSteps;
    so.nearMissWindow = tol / extent;
    return so;
}

// ---------------------------------------------------------------------------
long long SubdomainLabels::topoKeyOf(long long a, long long b) {
    const long long lo = std::min(a, b), hi = std::max(a, b);
    // Two 64-bit tokens folded into one. The mix is arbitrary; only collisions
    // between *different* pairs would matter and the token space is far larger
    // than any cone set here.
    return lo * 1000000007LL + hi;
}

// ---------------------------------------------------------------------------
void SubdomainLabels::truncateTopoPaths(size_t keep) {
    if (keep >= tPaths.size()) return;
    for (size_t i = keep; i < tPaths.size() && i < tPathKeys.size(); ++i) {
        topoKeys.erase(tPathKeys[i]);
        if (tPaths[i].pass > 0) --report.topoAddedByRepair;
        if (tPaths[i].fromCone == tPaths[i].toCone) --report.topoSelfReturns;
    }
    tPaths.resize(keep);
    tPathKeys.resize(std::min(tPathKeys.size(), keep));
    report.topoPaths = static_cast<int>(tPaths.size());
}

// ---------------------------------------------------------------------------
// addTopoPath()  --  from a traced separatrix
//
// The curve already is the path: Separatrices subdivides it by G as it marches
// and labels each subcurve with the coordinate E5 is to sum, so there is no
// geometry to recompute and, more to the point, no second opinion about where
// the seams were crossed or which quarter turn was applied there.
// ---------------------------------------------------------------------------
bool SubdomainLabels::addTopoPath(const Separatrices::Curve &c, int pass) {
    const int other = c.otherCone();
    if (c.cone < 0 || other < 0) return false;
    if (!options.seedSelfReturns && other == c.cone) return false;

    const int otherChild = (c.end == Separatrices::End::Cone) ? c.toChild : c.nearestChild;
    if (otherChild < 0) return false;

    const std::vector<Separatrices::Sub> subs = c.connection();
    if (subs.empty()) return false;

    // The two ends, each as (cone, child, direction of travel). The far end is
    // recorded as the direction a trace *from* there would leave in, which is
    // the reverse of the direction the curve arrived in -- so that the same
    // curve found from either end produces the same unordered pair.
    const int arriveDir = (c.end == Separatrices::End::Cone) ? c.endDir : c.nearestDir;
    const long long ka = endToken(c.cone, c.child, c.dir);
    const long long kb = endToken(other, otherChild, (arriveDir + 2) % 4);

    long long dedupA = ka, dedupB = kb;
    if (!options.seedAllConnections) {
        // The old behaviour: one path per pair of cones, whichever was found
        // first. Kept as a switch so the cost of it can be seen rather than
        // argued about.
        dedupA = endToken(std::min(c.cone, other), 0, 0);
        dedupB = endToken(std::max(c.cone, other), 0, 0);
    }
    const long long key = topoKeyOf(dedupA, dedupB);
    if (!topoKeys.insert(key).second) return false;

    TopoPath tp;
    tp.subs.reserve(subs.size());
    for (const Separatrices::Sub &s : subs) {
        TopoPath::Sub o;
        o.from = Site{s.from.a, s.from.b, s.from.t};
        o.to = Site{s.to.a, s.to.b, s.to.t};
        o.dir = s.dir;
        tp.subs.push_back(o);
    }
    tp.fromCone = c.cone;
    tp.toCone = other;
    tp.seedGap = (c.end == Separatrices::End::Cone) ? c.gap : c.nearestConeGap;
    tp.seamCrossings = c.seamCrossings;
    tp.pass = pass;

    if (other == c.cone) ++report.topoSelfReturns;
    if (pass > 0) ++report.topoAddedByRepair;
    tPathKeys.push_back(key);
    tPaths.push_back(std::move(tp));
    report.topoPaths = static_cast<int>(tPaths.size());
    return true;
}

// ---------------------------------------------------------------------------
// seedTopoPaths()  --  Sec. 3.3, "seeding Gamma_topo in practice"
//
// Trace every separatrix of the map with the snap tolerance widened from
// Sec. 3.4's 1e-6 to the near-miss tolerance, and keep as a connectivity
// constraint every curve that reached a cone under it. Widening the tolerance
// is exactly the "look for one that passes close to a cone without hitting it"
// of the recipe: under Sec. 3.4's tolerance it would not have hit, under this
// one it does, and the distance it was moved by is the residual E5 is then
// asked to drive to zero.
// ---------------------------------------------------------------------------
void SubdomainLabels::seedTopoPaths(const std::vector<Point> &uv) {
    tPaths.clear();
    tPathKeys.clear();
    topoKeys.clear();
    report.separatrices = report.separatricesToBoundary = 0;
    report.separatricesToCone = report.separatricesCapped = 0;
    report.topoPaths = report.topoSelfReturns = 0;
    report.topoExtraPerPair = report.topoAddedByRepair = 0;

    const std::vector<int> &coneVert = imm->getConeVertices();
    if (coneVert.empty()) return;

    report.meanConeSpacing = coneSpacing(uv);
    if (!(report.meanConeSpacing > 0.0)) {
        // One cone, or none: there is no pair to join and nothing to seed. Not
        // a failure -- a rectangle with four convex corners is exactly this.
        return;
    }

    const Separatrices::Options so = seedTracerOptions(uv, report.meanConeSpacing);
    report.seedSnapRings = so.coneSnapRings;

    std::unique_ptr<Separatrices> sep;
    try {
        sep = std::make_unique<Separatrices>(*imm, uv, so);
    } catch (const std::exception &e) {
        report.messages.push_back(
            std::string("Gamma_topo could not be seeded: the separatrices of psi_R would "
                        "not trace (") + e.what() + "). Remark 3.1's sliver layout is what "
            "Stage 6 will produce without them.");
        return;
    }

    const Separatrices::Report &sr = sep->getReport();
    report.seedSnapTolerance = sr.snapTolerance;
    report.separatrices = sr.emitted;
    report.separatricesToCone = sr.endedAtCone;
    report.separatricesToBoundary = sr.endedAtBoundary;
    report.separatricesCapped = sr.capped + sr.cycled + sr.stuck;

    // Count the pairs as they go by, so that the two things this used to leave
    // out can be reported rather than merely fixed.
    std::unordered_set<long long> pairsSeen;
    for (const Separatrices::Curve &c : sep->curves()) {
        if (c.end != Separatrices::End::Cone) continue;
        const int lo = std::min(c.cone, c.toCone), hi = std::max(c.cone, c.toCone);
        const long long pk = static_cast<long long>(lo) * 1000003LL + hi;
        const bool firstOfPair = pairsSeen.insert(pk).second;
        if (addTopoPath(c, 0) && !firstOfPair && c.cone != c.toCone) {
            ++report.topoExtraPerPair;
        }
    }

    report.maxTopoResidual = 0.0;
    for (const TopoPath &p : tPaths) {
        report.maxTopoResidual = std::max(report.maxTopoResidual,
                                          std::fabs(topoResidual(uv, p)));
    }

    {
        std::ostringstream oss;
        oss << "Gamma_topo seeded with " << report.topoPaths << " path(s) from "
            << report.separatrices << " separatrix/ces of psi_R, at a near-miss tolerance of "
            << options.nearMissTolerance << " of the mean cone spacing ("
            << report.seedSnapTolerance << " in image units, over " << report.seedSnapRings
            << " ring(s) of faces): " << report.topoSelfReturns << " back to their own cone, "
            << report.topoExtraPerPair << " a second curve between a pair already joined.";
        report.messages.push_back(oss.str());
    }
    if (report.separatricesCapped > 0) {
        std::ostringstream oss;
        oss << report.separatricesCapped << " separatrix/ces of psi_R neither reached a cone "
            << "nor left through dS; Q5 does not hold on psi_R, which is what Stage 6 is for.";
        report.messages.push_back(oss.str());
    }
}

// ---------------------------------------------------------------------------
// adoptCurves()  --  Sec. 3.3's repair, and Sec. 4's "on failure"
//
// The seeding above runs on psi_R because that is the only map there is when
// Stage 5 runs. Stage 6 then moves it a long way -- that is its whole job -- and
// a connection that is plain on the finished layout can have been invisible on
// the map it started from, and the other way round. So after Stage 7 has traced
// the separatrices of Psi, the ones that ended nowhere, and the ones that
// slipped past a cone on their way out through dS, each name a pair E5 was
// never told to join. Adding them and re-running the continuation *from the
// current phi* is the paper's remedy, and re-running it from psi_R would be
// throwing away everything the continuation has already bought.
// ---------------------------------------------------------------------------
int SubdomainLabels::adoptCurves(const Separatrices &sep, double window, int pass,
                                 int maxAdd) {
    // unconstrainedCurves() returns them nearest first, and nearest is the
    // right order to take them in: the curve that came closest to its cone is
    // the one the continuation was closest to getting right, so it is both the
    // most likely to have been meant and the cheapest to satisfy.
    const std::vector<int> wanted = sep.unconstrainedCurves(window);
    int added = 0;
    for (int i : wanted) {
        if (maxAdd > 0 && added >= maxAdd) break;
        if (addTopoPath(sep.curves()[i], pass)) ++added;
    }
    if (added > 0) {
        std::ostringstream oss;
        oss << "Repair pass " << pass << ": " << added << " connectivity constraint(s) added "
            << "for separatrices that terminated at neither a cone nor dS, or passed within "
            << window << " of the image extent of a cone without stopping at it";
        if (maxAdd > 0 && static_cast<int>(wanted.size()) > added) {
            oss << " (the " << added << " nearest of " << wanted.size() << " candidates)";
        }
        oss << ".";
        report.messages.push_back(oss.str());
    }
    report.topoPaths = static_cast<int>(tPaths.size());
    return added;
}
// ---------------------------------------------------------------------------
// addTopoPath()
//
// The manual form. Consecutive vertices either share an edge of Omega -- the
// path stays in one sheet -- or are the two children of the same point of the
// cutting graph, which is a crossing and turns the direction label by that
// arc's k.
// ---------------------------------------------------------------------------
bool SubdomainLabels::addTopoPath(const std::vector<int> &omegaPath, int firstDir) {
    if (omegaPath.size() < 2) return false;

    // (childOnOneSide, childOnTheOther) -> the arc they belong to.
    std::unordered_map<long long, int> crossing;
    const auto &arcList = imm->getArcs();
    for (int a = 0; a < static_cast<int>(arcList.size()); ++a) {
        const Immersion::Arc &arc = arcList[a];
        if (arc.degenerate) continue;
        for (size_t i = 0; i < arc.plusChain.size(); ++i) {
            const long long p = static_cast<long long>(arc.plusChain[i]) * 1000003LL + arc.minusChain[i];
            const long long m = static_cast<long long>(arc.minusChain[i]) * 1000003LL + arc.plusChain[i];
            crossing.emplace(p, +(a + 1));   // plus -> minus
            crossing.emplace(m, -(a + 1));   // minus -> plus
        }
    }

    TopoPath tp;
    int dirT = ((firstDir % 4) + 4) % 4;
    Site subStart{omegaPath.front(), -1, 0.0};

    for (size_t i = 1; i < omegaPath.size(); ++i) {
        const int a = omegaPath[i - 1], b = omegaPath[i];
        if (cutEdgeIndex.count(EdgeKey(a, b))) continue;   // same sheet, keep going

        auto it = crossing.find(static_cast<long long>(a) * 1000003LL + b);
        if (it == crossing.end()) return false;            // not an edge and not a crossing
        const int arcIdx = std::abs(it->second) - 1;
        const int k = arcList[arcIdx].k;
        const bool plusToMinus = it->second > 0;

        tp.subs.push_back({subStart, Site{a, -1, 0.0}, (dirT + 1) % 4});
        dirT = plusToMinus ? (((dirT - k) % 4 + 4) % 4) : (((dirT + k) % 4 + 4) % 4);
        subStart = Site{b, -1, 0.0};
        ++tp.seamCrossings;
    }
    tp.subs.push_back({subStart, Site{omegaPath.back(), -1, 0.0}, (dirT + 1) % 4});

    const std::vector<int> &vCone = imm->getCutVertexCone();
    tp.fromCone = vCone[omegaPath.front()];
    tp.toCone = vCone[omegaPath.back()];

    // A hand-written path claims its pair of ends like a traced one, so that
    // the seeding and the repair do not later propose the same connection --
    // and so that tPathKeys stays parallel to tPaths for truncateTopoPaths().
    const long long key = topoKeyOf(
        endToken(tp.fromCone, omegaPath.front(), firstDir),
        endToken(tp.toCone, omegaPath.back(), (dirT + 2) % 4));
    topoKeys.insert(key);
    tPathKeys.push_back(key);
    tPaths.push_back(std::move(tp));
    report.topoPaths = static_cast<int>(tPaths.size());
    return true;
}

// ---------------------------------------------------------------------------
double SubdomainLabels::topoResidual(const std::vector<Point> &uv, const TopoPath &p) const {
    double sum = 0.0;
    for (const TopoPath::Sub &s : p.subs) {
        sum += dotP(axis(s.dir), evaluate(uv, s.to) - evaluate(uv, s.from));
    }
    return sum;
}
