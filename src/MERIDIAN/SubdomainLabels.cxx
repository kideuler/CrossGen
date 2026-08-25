#include "SubdomainLabels.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <sstream>
#include <unordered_set>

namespace {

// Distance from a point to a segment, and where along the segment that was.
double pointSegmentDistance(const Point &p, const Point &a, const Point &b) {
    const Point ab = b - a;
    const double denom = dotP(ab, ab);
    double t = 0.0;
    if (denom > 0.0) t = std::max(0.0, std::min(1.0, dotP(p - a, ab) / denom));
    return normP(p - (a + ab * t));
}

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
    : imm(&immersion), options(opts) {
    buildSeamTables();

    report.seamArcs = static_cast<int>(imm->getArcs().size());
    for (const Immersion::Arc &a : imm->getArcs()) {
        if (!a.degenerate) ++report.holonomyCount[((a.k % 4) + 4) % 4];
    }

    const std::vector<Point> &uv = imm->getUV();
    buildBoundary(uv);
    buildBoundaryChains(uv);
    buildFeatures(uv);
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
    vertexOnRealBoundary.assign(cm.vertices.size(), 0);
    for (int e = 0; e < static_cast<int>(cm.edges.size()); ++e) {
        auto it = origEdgeIndex.find(EdgeKey(c2o[cm.edges[e][0]], c2o[cm.edges[e][1]]));
        if (it != origEdgeIndex.end() && om.isBoundaryEdge[it->second]) {
            parentOnBoundary[e] = 1;
            vertexOnRealBoundary[cm.edges[e][0]] = 1;
            vertexOnRealBoundary[cm.edges[e][1]] = 1;
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

    // The feature edges of Omega.
    std::unordered_set<int> featEdges;
    for (int e = 0; e < static_cast<int>(om.edges.size()); ++e) {
        const int f0 = om.edgeTriangles[e][0];
        const int f1 = om.edgeTriangles[e][1];
        if (f0 < 0 || f1 < 0) continue;                          // dS, handled by E2
        if (om.triangleMatId[f0] == om.triangleMatId[f1]) continue;

        const int a = om.edges[e][0], b = om.edges[e][1];
        for (int f : {f0, f1}) {
            const int la = localOf(om.triangles[f], a);
            const int lb = localOf(om.triangles[f], b);
            if (la < 0 || lb < 0) continue;
            auto it = cutEdgeIndex.find(EdgeKey(cm.triangles[f][la], cm.triangles[f][lb]));
            if (it != cutEdgeIndex.end()) featEdges.insert(it->second);
        }
    }
    report.featureEdges = static_cast<int>(featEdges.size());
    if (featEdges.empty()) return;

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
        for (int side = 0; side < 2; ++side) {
            int cur = seed;
            int v = (side == 0) ? cm.edges[seed][0] : cm.edges[seed][1];
            double turned = 0.0;
            while (at[v].size() == 2) {
                int next = -1;
                for (int c : at[v]) if (c != cur) next = c;
                if (next < 0 || used.count(next)) break;
                const double t = turnAt(v, cur, next);
                if (std::fabs(t) > options.featureCornerAngle) break;
                turned += t;
                if (std::fabs(turned) > options.featureCornerAngle) break;
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
        fChains.push_back(std::move(fc));
    }
    report.featureChains = static_cast<int>(fChains.size());
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
        fc.label = (fc.fluxU <= fc.fluxV) ? Align::U : Align::V;
        (fc.label == Align::U ? report.featureChainsU : report.featureChainsV) += 1;
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
// traceSeparatrix()  --  the marching of Sec. 3.4, used here only as a seeder
//
// Inside a triangle phi is affine, so an isoline of the constrained coordinate
// is a straight segment and the exit point is a linear interpolation along one
// edge. Three things can end the walk, and they are the three cases Q5 lists:
// the curve reaches a cone, it leaves through dS, or it runs out of steps --
// the last being the diagnosis "Q5 does not hold yet" rather than a bug.
//
// A seam is not an end. The curve leaves Omega through the child edge on one
// side of an arc and re-enters through the child edge on the other, with the
// point carried across by the arc's transition and the direction rotated by its
// k quarter turns. That rotation is exactly what makes the direction label of
// the next subcurve differ from this one's, and it is the whole reason E5 needs
// a label per subcurve instead of one per path.
// ---------------------------------------------------------------------------
SubdomainLabels::Trace SubdomainLabels::traceSeparatrix(const std::vector<Point> &uv,
                                                        int startVertex, int face,
                                                        int tangentDir, double coneTol) const {
    Trace tr;
    const Mesh &cm = imm->getCutMesh();
    const auto &pairs = imm->getSeamPairs();
    const auto &pairArc = imm->getSeamPairArc();
    const auto &arcList = imm->getArcs();
    const std::vector<int> &vCone = imm->getCutVertexCone();
    const std::vector<int> &coneVert = imm->getConeVertices();

    const int startCone = vCone[startVertex];

    const Point &lo = imm->getReport().uvMin;
    const Point &hi = imm->getReport().uvMax;
    const double scale = std::max(hi[0] - lo[0], hi[1] - lo[1]);
    const double tiny = 1e-12 * std::max(scale, 1e-30);
    const double vertexEps = 1e-9 * std::max(scale, 1e-30);

    int f = face;
    int entryEdge = -1;
    int dirT = ((tangentDir % 4) + 4) % 4;
    Site subStart{startVertex, -1, 0.0};
    Point p = uv[startVertex];
    int fanSpins = 0;   // consecutive vertex rotations, so a fan cannot cycle

    for (int step = 0; step < options.maxTraceSteps; ++step) {
        const int tc = dirT % 2;          // the coordinate the curve advances in
        const int cc = 1 - tc;            // the one it holds constant
        const double s = (dirT < 2) ? 1.0 : -1.0;
        const double c0 = p[cc];

        const Triangle &t = cm.triangles[f];
        int best = -1;
        double bestT = 0.0, bestAdv = tiny;
        Point bestP{0.0, 0.0};

        for (int q = 0; q < 3; ++q) {
            if (q == entryEdge) continue;
            const int m = q, n = (q + 1) % 3;
            const double A = uv[t[m]][cc];
            const double B = uv[t[n]][cc];
            const double den = B - A;
            if (std::fabs(den) < 1e-300) continue;   // the edge is on the isoline
            const double tt = (c0 - A) / den;
            if (tt < -1e-9 || tt > 1.0 + 1e-9) continue;
            const Point P = uv[t[m]] * (1.0 - tt) + uv[t[n]] * tt;
            const double adv = (P[tc] - p[tc]) * s;
            if (adv > bestAdv) { bestAdv = adv; best = q; bestT = tt; bestP = P; }
        }
        if (best < 0) {
            // No exit edge, which on a non-degenerate triangle means one thing:
            // the curve arrived exactly at a *vertex*, and continues into a
            // different triangle of that vertex's fan rather than out of this
            // one. It is not a rare case -- Stage 4 puts the boundary on the
            // axes, so an axis-parallel isoline runs into boundary vertices
            // head-on all the time -- and without it the walk simply stops in
            // the middle of the domain and Q5 looks unsatisfiable when it is
            // not. Rotate around the fan to whichever triangle's angular sector
            // contains the direction of travel.
            int w = -1;
            for (int q = 0; q < 3; ++q) {
                if (normP(p - uv[t[q]]) <= vertexEps) { w = t[q]; break; }
            }
            int nf = -1;
            if (w >= 0 && fanSpins < 32) {
                const Point e = axis(dirT);
                const auto &vt = cm.vertexTriangles;
                for (int i = vt.rowPtr[w]; i < vt.rowPtr[w + 1] && nf < 0; ++i) {
                    const int g = vt.colIdx[i];
                    if (g == f) continue;
                    const int lw = localOf(cm.triangles[g], w);
                    if (lw < 0) continue;
                    const Point A = normalizeP(uv[cm.triangles[g][(lw + 1) % 3]] - uv[w]);
                    const Point B = normalizeP(uv[cm.triangles[g][(lw + 2) % 3]] - uv[w]);
                    if (cross2(A, e) > -1e-9 && cross2(e, B) > 1e-9) nf = g;
                }
            }
            if (nf >= 0) {
                ++fanSpins;
                f = nf;
                entryEdge = -1;
                p = uv[w];
                continue;
            }
            // Nothing in the fan takes the direction: either the curve leaves
            // the domain at this vertex, or the fan is a seam and the walk
            // cannot be continued without a transition at a point rather than
            // across an edge, which is beyond what a seeder needs.
            if (w >= 0 && vertexOnRealBoundary[w]) {
                tr.subs.push_back({subStart, Site{w, -1, 0.0}, (dirT + 1) % 4});
                tr.leftDomain = true;
                return tr;
            }
            if (step == 0) tr.degenerate = true;
            tr.capped = true;
            break;
        }
        fanSpins = 0;

        // Q5 case 1: the curve came near a cone. The snap is what turns a
        // near-miss into a connectivity constraint, and how near it came is the
        // residual E5 will be asked to drive to zero.
        int snapped = -1;
        double snappedGap = coneTol;
        for (int slot = 0; slot < static_cast<int>(coneVert.size()); ++slot) {
            if (slot == startCone) continue;
            for (int child : imm->getConeChildren()[slot]) {
                if (child == startVertex) continue;
                const double d = pointSegmentDistance(uv[child], p, bestP);
                if (d < snappedGap) { snappedGap = d; snapped = slot; }
            }
        }
        if (snapped >= 0 && step >= 1) {
            int child = -1;
            double bestD = std::numeric_limits<double>::infinity();
            for (int c : imm->getConeChildren()[snapped]) {
                const double d = pointSegmentDistance(uv[c], p, bestP);
                if (d < bestD) { bestD = d; child = c; }
            }
            tr.subs.push_back({subStart, Site{child, -1, 0.0}, (dirT + 1) % 4});
            tr.hitCone = snapped;
            tr.gap = snappedGap;
            return tr;
        }

        const int e = cm.triangleEdges[f][best];
        const int nb = cm.triangleAdjacency[f][best];
        const int x = t[best], y = t[(best + 1) % 3];

        if (nb >= 0) {
            const Triangle &tn = cm.triangles[nb];
            const int lx = localOf(tn, x), ly = localOf(tn, y);
            entryEdge = (lx < 0 || ly < 0) ? -1 : localEdgeBetween(lx, ly);
            f = nb;
            p = bestP;
            continue;
        }

        // Q5 case 2: out through dS.
        if (parentOnBoundary[e]) {
            tr.subs.push_back({subStart, Site{x, y, bestT}, (dirT + 1) % 4});
            tr.leftDomain = true;
            return tr;
        }

        // Q5 case 3: across a seam.
        auto it = seamSide.find(EdgeKey(x, y));
        if (it == seamSide.end()) { tr.capped = true; break; }
        const int pairIdx = it->second / 2;
        const bool onPlus = (it->second % 2) == 0;
        const auto &pr = pairs[pairIdx];
        const int k = arcList[pairArc[pairIdx]].k;

        int ta = -1, tb = -1, newDir = dirT;
        double tPar = bestT;
        if (onPlus) {
            // Leaving through the plus child, so the transition to apply is the
            // inverse and the direction turns back by k.
            tPar = (x == pr[0]) ? bestT : 1.0 - bestT;
            ta = pr[2]; tb = pr[3];
            newDir = ((dirT - k) % 4 + 4) % 4;
        } else {
            tPar = (x == pr[2]) ? bestT : 1.0 - bestT;
            ta = pr[0]; tb = pr[1];
            newDir = ((dirT + k) % 4 + 4) % 4;
        }

        auto ie = cutEdgeIndex.find(EdgeKey(ta, tb));
        if (ie == cutEdgeIndex.end()) { tr.capped = true; break; }
        const int te = ie->second;
        const int tf = cm.edgeTriangles[te][0] >= 0 ? cm.edgeTriangles[te][0]
                                                    : cm.edgeTriangles[te][1];
        if (tf < 0) { tr.capped = true; break; }

        tr.subs.push_back({subStart, Site{x, y, bestT}, (dirT + 1) % 4});
        ++tr.seamCrossings;

        const Triangle &tt2 = cm.triangles[tf];
        const int la = localOf(tt2, ta), lb = localOf(tt2, tb);
        entryEdge = (la < 0 || lb < 0) ? -1 : localEdgeBetween(la, lb);
        f = tf;
        dirT = newDir;
        subStart = Site{ta, tb, tPar};
        p = evaluate(uv, subStart);
    }

    tr.capped = true;
    return tr;
}

// ---------------------------------------------------------------------------
// seedTopoPaths()  --  Sec. 3.3, "seeding Gamma_topo in practice"
//
// Emit the separatrices of every cone, follow each one, and whenever one passes
// within tolerance of another cone, keep the traced path as a connectivity
// constraint between the two.
//
// The directions are read straight off the map rather than by walking the fan
// and accumulating angle: a ray leaves the cone along +u, +v, -u or -v, and it
// belongs to whichever incident triangle's angular sector contains it. That is
// the same set of rays -- a cone of index I has total angle 2pi - (pi/2) I, so
// four axis directions swept over 5pi/2 give five rays at a valence-five cone
// and three over 3pi/2 give three at a valence-three one, which is exactly the
// n = 4 - I of Sec. 3.4 -- without needing a fan order that a boundary vertex
// may not have.
// ---------------------------------------------------------------------------
void SubdomainLabels::seedTopoPaths(const std::vector<Point> &uv) {
    tPaths.clear();
    report.separatrices = report.separatricesToBoundary = 0;
    report.separatricesToCone = report.separatricesCapped = 0;
    report.topoPaths = 0;

    const Mesh &cm = imm->getCutMesh();
    const std::vector<int> &coneVert = imm->getConeVertices();
    const auto &children = imm->getConeChildren();
    if (coneVert.size() < 2) return;

    // The scale a "near miss" is measured against: the mean distance from each
    // cone to its nearest other cone, in the image.
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
    report.meanConeSpacing = counted ? acc / counted : 0.0;
    if (!(report.meanConeSpacing > 0.0)) return;
    const double tol = options.nearMissTolerance * report.meanConeSpacing;

    // slot pair -> index into tPaths, so a pair of cones found twice keeps only
    // the trace that came closest.
    std::unordered_map<long long, int> seen;

    for (int slot = 0; slot < static_cast<int>(coneVert.size()); ++slot) {
        for (int c : children[slot]) {
            // The directions that run *along* the boundary of Omega at this
            // cone. Sec. 3.4: a cone on the boundary of valence n emits n - 2
            // separatrices into the interior, because two of its n grid
            // directions lie in dS. Since Stage 4 put the boundary on the axes,
            // those two directions are exactly axis directions and would
            // otherwise be emitted and then traced along the boundary edge.
            std::vector<Point> along;
            {
                const auto &vt0 = cm.vertexTriangles;
                for (int i = vt0.rowPtr[c]; i < vt0.rowPtr[c + 1]; ++i) {
                    const int f0 = vt0.colIdx[i];
                    for (int q = 0; q < 3; ++q) {
                        if (cm.triangleAdjacency[f0][q] >= 0) continue;
                        const int a0 = cm.triangles[f0][q];
                        const int b0 = cm.triangles[f0][(q + 1) % 3];
                        if (a0 != c && b0 != c) continue;
                        const int other = (a0 == c) ? b0 : a0;
                        along.push_back(normalizeP(uv[other] - uv[c]));
                    }
                }
            }
            auto runsAlongBoundary = [&](const Point &e) {
                for (const Point &d0 : along) {
                    if (std::fabs(cross2(e, d0)) < 1e-9 && dotP(e, d0) > 0.0) return true;
                }
                return false;
            };

            const auto &vt = cm.vertexTriangles;
            for (int i = vt.rowPtr[c]; i < vt.rowPtr[c + 1]; ++i) {
                const int f = vt.colIdx[i];
                const Triangle &t = cm.triangles[f];
                const int lc = localOf(t, c);
                if (lc < 0) continue;
                const Point A = uv[t[(lc + 1) % 3]] - uv[c];
                const Point B = uv[t[(lc + 2) % 3]] - uv[c];

                const Point An = normalizeP(A);
                const Point Bn = normalizeP(B);

                for (int d = 0; d < 4; ++d) {
                    const Point e = axis(d);
                    // The sector of a positively oriented triangle at c runs
                    // from A to B counter-clockwise. The test is half-open --
                    // closed at A, open at B -- so a direction lying exactly
                    // along a shared edge belongs to one of the two faces and
                    // not to both. That is not a pedantic edge case here:
                    // Stage 4 rotates the whole immersion onto the axes, so a
                    // boundary that came out rectilinear has many edges exactly
                    // parallel to a grid direction.
                    if (!(cross2(An, e) > -1e-9 && cross2(e, Bn) > 1e-9)) continue;
                    if (runsAlongBoundary(e)) continue;

                    Trace tr = traceSeparatrix(uv, c, f, d, tol);
                    // A ray that could not take even one step was not a
                    // separatrix; it left along an edge of the fan.
                    if (tr.degenerate) continue;
                    ++report.separatrices;
                    if (tr.leftDomain) ++report.separatricesToBoundary;
                    else if (tr.hitCone >= 0) ++report.separatricesToCone;
                    else ++report.separatricesCapped;

                    if (tr.hitCone < 0 || tr.hitCone == slot || tr.subs.empty()) continue;

                    const int lo = std::min(slot, tr.hitCone);
                    const int hi = std::max(slot, tr.hitCone);
                    const long long key = static_cast<long long>(lo) * 1000003LL + hi;

                    TopoPath tp;
                    tp.subs = std::move(tr.subs);
                    tp.fromCone = slot;
                    tp.toCone = tr.hitCone;
                    tp.seedGap = tr.gap;
                    tp.seamCrossings = tr.seamCrossings;

                    auto it = seen.find(key);
                    if (it == seen.end()) {
                        seen.emplace(key, static_cast<int>(tPaths.size()));
                        tPaths.push_back(std::move(tp));
                    } else if (tp.seedGap < tPaths[it->second].seedGap) {
                        tPaths[it->second] = std::move(tp);
                    }
                }
            }
        }
    }

    report.topoPaths = static_cast<int>(tPaths.size());
    report.maxTopoResidual = 0.0;
    for (const TopoPath &p : tPaths) {
        report.maxTopoResidual = std::max(report.maxTopoResidual,
                                          std::fabs(topoResidual(uv, p)));
    }
    if (report.separatricesCapped > 0) {
        std::ostringstream oss;
        oss << report.separatricesCapped << " separatrix/ces neither reached a cone nor left "
            << "through dS within the step cap; Q5 does not hold on psi_R, which is what "
            << "Stage 6 is for.";
        report.messages.push_back(oss.str());
    }
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
