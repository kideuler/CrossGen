#include "ATLAS/RectangleCertifier.hxx"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>

namespace {

const long long kCoordLimit = 1LL << 40;
const long long kKeyOffset = 1LL << 31;

inline unsigned long long coordKey(long long x, long long y) {
    return (static_cast<unsigned long long>(x + kKeyOffset) << 32) |
           static_cast<unsigned long long>(y + kKeyOffset);
}

unsigned long long cellSetHash(std::vector<int> cells) {
    std::sort(cells.begin(), cells.end());
    unsigned long long h = 1469598103934665603ULL;
    for (int q : cells) {
        h ^= static_cast<unsigned long long>(q) + 0x9e3779b97f4a7c15ULL;
        h *= 1099511628211ULL;
    }
    h ^= cells.size();
    return h;
}

} // namespace

const char *RectangleCertifier::failureName(Failure f) {
    switch (f) {
        case Failure::None: return "none";
        case Failure::Empty: return "empty";
        case Failure::Disconnected: return "disconnected";
        case Failure::TransportConflict: return "transport conflict (holonomy on a cycle)";
        case Failure::VertexPathDependence: return "vertex path dependence (periodic)";
        case Failure::VertexCollision: return "vertex collision (extra identification)";
        case Failure::DuplicateOccupancy: return "duplicate occupancy";
        case Failure::MissingSquare: return "missing square";
        case Failure::FeatureInside: return "interface inside";
        case Failure::ProtectedSwallowed: return "protected vertex swallowed";
        case Failure::Overflow: return "coordinate overflow";
    }
    return "?";
}

// ---------------------------------------------------------------------------
// Sec. 5.1 and 5.2
// ---------------------------------------------------------------------------
bool RectangleCertifier::certify(const std::vector<int> &cells, Certificate &cert, Witness &w) const {
    ++report_.certifyCalls;
    report_.certifiedCells += static_cast<long long>(cells.size());
    w = Witness();
    auto fail = [&](Failure f, int a = -1, int b = -1, int v = -1) {
        w.kind = f;
        w.cellA = a;
        w.cellB = b;
        w.vertex = v;
        return false;
    };
    const int k = static_cast<int>(cells.size());
    if (k == 0) return fail(Failure::Empty);
    const int NQ = C_.numCells();
    if (static_cast<int>(stamp_.size()) != NQ) {
        stamp_.assign(NQ, 0);
        local_.assign(NQ, -1);
        stampValue_ = 0;
    }
    ++stampValue_;
    for (int i = 0; i < k; ++i) {
        stamp_[cells[i]] = stampValue_;
        local_[cells[i]] = i;
    }
    auto inP = [&](int q) { return q >= 0 && stamp_[q] == stampValue_; };

    // ---- development ----------------------------------------------------
    std::vector<SquareTransport> G(k);
    std::vector<char> seen(k, 0);
    std::vector<int> queue{0};
    seen[0] = 1;
    for (size_t head = 0; head < queue.size(); ++head) {
        const int li = queue[head];
        const int q = cells[li];
        for (int i = 0; i < 4; ++i) {
            const int r = C_.neighbor[q][i];
            if (!inP(r)) continue;
            const int ri = local_[r];
            const SquareTransport Gr = G[li].compose(C_.transport[q][i]);
            if (!seen[ri]) {
                seen[ri] = 1;
                G[ri] = Gr;
                queue.push_back(ri);
            } else if (G[ri] != Gr) {
                return fail(Failure::TransportConflict, q, r);
            }
        }
    }
    if (static_cast<int>(queue.size()) != k) return fail(Failure::Disconnected, cells[0]);

    // ---- vertices: one coordinate each, one vertex per coordinate ---------
    std::unordered_map<int, IPoint> vcoord;
    std::unordered_map<unsigned long long, int> cvert;
    vcoord.reserve(k * 2 + 8);
    cvert.reserve(k * 2 + 8);
    long long x0 = 0, y0 = 0, x1 = 0, y1 = 0;
    bool firstV = true;
    for (int li = 0; li < k; ++li) {
        const int q = cells[li];
        for (int c = 0; c < 4; ++c) {
            const IPoint p = G[li].apply(squareCorner(c));
            if (std::llabs(p[0]) > kCoordLimit || std::llabs(p[1]) > kCoordLimit) return fail(Failure::Overflow, q);
            const int v = C_.cells[q][c];
            auto it = vcoord.find(v);
            if (it != vcoord.end()) {
                if (it->second != p) return fail(Failure::VertexPathDependence, q, -1, v);
            } else {
                vcoord.emplace(v, p);
            }
            const unsigned long long key = coordKey(p[0], p[1]);
            auto jt = cvert.find(key);
            if (jt != cvert.end()) {
                if (jt->second != v) return fail(Failure::VertexCollision, q, -1, v);
            } else {
                cvert.emplace(key, v);
            }
            if (firstV) { x0 = x1 = p[0]; y0 = y1 = p[1]; firstV = false; }
            x0 = std::min(x0, p[0]); x1 = std::max(x1, p[0]);
            y0 = std::min(y0, p[1]); y1 = std::max(y1, p[1]);
        }
    }

    // ---- squares: distinct, and all of the bounding rectangle -------------
    std::unordered_map<unsigned long long, int> occupied;
    occupied.reserve(k * 2 + 8);
    std::vector<IPoint> square(k);
    for (int li = 0; li < k; ++li) {
        IPoint lo{1LL << 50, 1LL << 50};
        for (int c = 0; c < 4; ++c) {
            const IPoint p = G[li].apply(squareCorner(c));
            lo[0] = std::min(lo[0], p[0]);
            lo[1] = std::min(lo[1], p[1]);
        }
        square[li] = lo;
        if (!occupied.emplace(coordKey(lo[0], lo[1]), li).second) {
            return fail(Failure::DuplicateOccupancy, cells[li], cells[occupied[coordKey(lo[0], lo[1])]]);
        }
    }
    const long long nu = x1 - x0, nv = y1 - y0;
    if (nu * nv != k) return fail(Failure::MissingSquare, cells[0]);

    // ---- condition 5: features and protected vertices ---------------------
    for (int li = 0; li < k; ++li) {
        const int q = cells[li];
        for (int i = 0; i < 4; ++i) {
            if (inP(C_.neighbor[q][i]) && C_.interfaceEdge[C_.cellEdges[q][i]]) {
                return fail(Failure::FeatureInside, q, C_.neighbor[q][i]);
            }
        }
    }
    for (const auto &kv : vcoord) {
        if (!C_.protectedVertex[kv.first]) continue;
        const IPoint &p = kv.second;
        const bool corner = (p[0] == x0 || p[0] == x1) && (p[1] == y0 || p[1] == y1);
        if (!corner) return fail(Failure::ProtectedSwallowed, -1, -1, kv.first);
    }

    // ---- the certificate ---------------------------------------------------
    cert = Certificate();
    cert.nu = static_cast<int>(nu);
    cert.nv = static_cast<int>(nv);
    cert.cells.assign(k, -1);
    cert.G.assign(k, SquareTransport());
    const SquareTransport shift = SquareTransport::translation(-x0, -y0);
    for (int li = 0; li < k; ++li) {
        const long long i = square[li][0] - x0, j = square[li][1] - y0;
        const int pos = static_cast<int>(i + nu * j);
        cert.cells[pos] = cells[li];
        cert.G[pos] = shift.compose(G[li]);
    }
    auto vat = [&](long long x, long long y) { return cvert[coordKey(x, y)]; };
    cert.corners = {vat(x0, y0), vat(x1, y0), vat(x1, y1), vat(x0, y1)};
    for (long long x = x0; x <= x1; ++x) cert.sides[0].push_back(vat(x, y0));
    for (long long y = y0; y <= y1; ++y) cert.sides[1].push_back(vat(x1, y));
    for (long long x = x1; x >= x0; --x) cert.sides[2].push_back(vat(x, y1));
    for (long long y = y1; y >= y0; --y) cert.sides[3].push_back(vat(x0, y));
    cert.material = C_.cellMaterial[cells[0]];
    measure(C_, cert);
    return true;
}

// D(P): the conformal distortion |J|_F^2 / (2 det J) - 1 of each piece at its
// corners, averaged -- zero for a square, and blind to size, which Sec. 10's
// "size-normalized" asks for. C(P): the total turning of the four sides in
// quarter turns, the map/side complexity of Sec. 6 in the one form that does
// not depend on how finely the side happens to be sampled.
void RectangleCertifier::measure(const SquareCarrier &C, Certificate &cert) {
    double d = 0.0;
    int n = 0;
    for (int q : cert.cells) {
        for (int c = 0; c < 4; ++c) {
            const Point &p = C.vertices[C.cells[q][c]];
            const Point u = C.vertices[C.cells[q][(c + 1) & 3]] - p;
            const Point w = C.vertices[C.cells[q][(c + 3) & 3]] - p;
            const double det = cross2(u, w);
            if (det > 0.0) d += (dotP(u, u) + dotP(w, w)) / (2.0 * det) - 1.0;
            else d += 1e3;
            ++n;
        }
    }
    cert.distortion = n > 0 ? d / n : 0.0;
    double turn = 0.0;
    for (const auto &side : cert.sides) {
        for (size_t k = 1; k + 1 < side.size(); ++k) {
            const Point a = C.vertices[side[k]] - C.vertices[side[k - 1]];
            const Point b = C.vertices[side[k + 1]] - C.vertices[side[k]];
            turn += std::fabs(std::atan2(cross2(a, b), dotP(a, b)));
        }
    }
    cert.complexity = turn / M_PI_2;
}

// ---------------------------------------------------------------------------
// Candidate generation
// ---------------------------------------------------------------------------
std::vector<std::vector<int>> RectangleCertifier::baseComplex(const std::vector<char> &forced) const {
    const int NE = C_.numEdges(), NV = C_.numVertices(), NQ = C_.numCells();
    std::vector<char> cut(NE, 0);
    for (int e = 0; e < NE; ++e) cut[e] = C_.isFeatureEdge(e) ? 1 : 0;
    std::vector<int> es, es2;
    for (int v = 0; v < NV; ++v) {
        if (!forced[v]) continue;
        C_.vertexEdges(v, es);
        for (int e0 : es) {
            if (C_.isFeatureEdge(e0)) continue;
            int u = v, e = e0;
            for (int guard = 0; guard < NE + 1; ++guard) {
                if (cut[e] == 2) break;
                cut[e] = 2;
                const int x = C_.otherEnd(e, u);
                if (forced[x] || C_.boundaryVertex[x] || C_.valence[x] != 4) break;
                C_.vertexEdges(x, es2);
                if (es2.size() != 4) break;
                int k = 0;
                while (k < 4 && es2[k] != e) ++k;
                if (k == 4) break;
                const int nxt = es2[(k + 2) & 3];
                if (C_.isFeatureEdge(nxt)) break;
                u = x;
                e = nxt;
            }
        }
    }
    std::vector<int> patchOf(NQ, -1);
    std::vector<std::vector<int>> patches;
    for (int s = 0; s < NQ; ++s) {
        if (patchOf[s] >= 0) continue;
        const int id = static_cast<int>(patches.size());
        patches.emplace_back();
        std::vector<int> stack{s};
        patchOf[s] = id;
        while (!stack.empty()) {
            const int q = stack.back();
            stack.pop_back();
            patches[id].push_back(q);
            for (int i = 0; i < 4; ++i) {
                const int r = C_.neighbor[q][i];
                if (r < 0 || patchOf[r] >= 0 || cut[C_.cellEdges[q][i]]) continue;
                patchOf[r] = id;
                stack.push_back(r);
            }
        }
    }
    return patches;
}

int RectangleCertifier::basePatchCount(const SquareCarrier &C) {
    const int NE = C.numEdges(), NV = C.numVertices(), NQ = C.numCells();
    std::vector<char> cut(NE, 0);
    for (int e = 0; e < NE; ++e) cut[e] = C.isFeatureEdge(e) ? 1 : 0;
    std::vector<int> es, es2;
    for (int v = 0; v < NV; ++v) {
        if (!(C.forcedMacrovertex(v) || C.designatedVertex[v])) continue;
        C.vertexEdges(v, es);
        for (int e0 : es) {
            if (C.isFeatureEdge(e0)) continue;
            int u = v, e = e0;
            for (int guard = 0; guard < NE + 1; ++guard) {
                if (cut[e] == 2) break;
                cut[e] = 2;
                const int x = C.otherEnd(e, u);
                if (C.boundaryVertex[x] || C.valence[x] != 4 || C.designatedVertex[x] || C.forcedMacrovertex(x)) break;
                C.vertexEdges(x, es2);
                if (es2.size() != 4) break;
                int k = 0;
                while (k < 4 && es2[k] != e) ++k;
                if (k == 4) break;
                const int nxt = es2[(k + 2) & 3];
                if (C.isFeatureEdge(nxt)) break;
                u = x;
                e = nxt;
            }
        }
    }
    std::vector<int> seen(NQ, 0);
    int patches = 0;
    std::vector<int> stack;
    for (int s = 0; s < NQ; ++s) {
        if (seen[s]) continue;
        ++patches;
        stack.assign(1, s);
        seen[s] = 1;
        while (!stack.empty()) {
            const int q = stack.back();
            stack.pop_back();
            for (int i = 0; i < 4; ++i) {
                const int r = C.neighbor[q][i];
                if (r < 0 || seen[r] || cut[C.cellEdges[q][i]]) continue;
                seen[r] = 1;
                stack.push_back(r);
            }
        }
    }
    return patches;
}

int RectangleCertifier::addCandidate(Certificate &&cert) {
    const unsigned long long h = cellSetHash(cert.cells);
    auto it = seen_.find(h);
    if (it != seen_.end()) return it->second;
    if (static_cast<int>(candidates_.size()) >= opts_.maxCandidates) return -1;
    const int id = static_cast<int>(candidates_.size());
    report_.largestCandidate = std::max(report_.largestCandidate, static_cast<int>(cert.cells.size()));
    candidates_.push_back(std::move(cert));
    seen_.emplace(h, id);
    return id;
}

void RectangleCertifier::certifiedBase(std::vector<char> forced, std::vector<int> &cover,
                                       const char *source, int &patchCount) {
    cover.clear();
    const int NQ = C_.numCells();
    for (int round = 0; round <= opts_.maxCutRounds; ++round) {
        const std::vector<std::vector<int>> patches = baseComplex(forced);
        std::vector<Certificate> certs(patches.size());
        std::vector<int> cuts;
        std::vector<int> patchOf(NQ, -1);
        for (size_t p = 0; p < patches.size(); ++p) for (int q : patches[p]) patchOf[q] = static_cast<int>(p);

        for (size_t p = 0; p < patches.size(); ++p) {
            Witness w;
            if (certify(patches[p], certs[p], w)) continue;
            ++report_.failures;
            ++report_.failuresByKind[static_cast<int>(w.kind)];
            if (round == 0 && witnesses_.size() < 100000) witnesses_.push_back(w);
            // A cut through the failing patch: a vertex on its boundary,
            // preferring dS (one interior edge, so exactly one new line), and
            // among those the farthest from the witness, so a band cut once is
            // cut again opposite rather than beside the first cut.
            std::vector<int> ring;
            std::set<int> vs;
            for (int q : patches[p]) for (int v : C_.cells[q]) vs.insert(v);
            int bestV = -1, bestScore = -1;
            // Graph distance over the patch's vertices from the witness.
            std::map<int, int> dist;
            int src = w.vertex >= 0 ? w.vertex : C_.cells[patches[p][0]][0];
            dist[src] = 0;
            std::vector<int> bfs{src};
            for (size_t h = 0; h < bfs.size(); ++h) {
                const int v = bfs[h];
                for (int r = C_.ringPtr[v]; r < C_.ringPtr[v + 1]; ++r) {
                    const int q = C_.ringCell[r];
                    if (patchOf[q] != static_cast<int>(p)) continue;
                    const int c = C_.ringCorner[r];
                    for (int nb : {C_.cells[q][(c + 1) & 3], C_.cells[q][(c + 3) & 3]}) {
                        if (!dist.count(nb)) { dist[nb] = dist[v] + 1; bfs.push_back(nb); }
                    }
                }
            }
            for (int v : vs) {
                if (forced[v]) continue;
                bool onPatchBoundary = C_.boundaryVertex[v];
                for (int r = C_.ringPtr[v]; r < C_.ringPtr[v + 1] && !onPatchBoundary; ++r) {
                    if (patchOf[C_.ringCell[r]] != static_cast<int>(p)) onPatchBoundary = true;
                }
                if (!onPatchBoundary) continue;
                const int d = dist.count(v) ? dist[v] : 0;
                const int score = d * 2 + (C_.boundaryVertex[v] ? 1 : 0) + (C_.boundaryVertex[v] ? 1000000 : 0);
                if (score > bestScore) { bestScore = score; bestV = v; }
            }
            if (bestV < 0) {
                for (int v : vs) if (!forced[v]) { bestV = v; break; }
            }
            if (bestV >= 0) cuts.push_back(bestV);
        }

        // Sec. 6: two blocks may not touch along two separate sides. A pair of
        // certified patches that do gets a cut at the middle of a shared side,
        // which splits both.
        if (cuts.empty()) {
            std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> contacts;  // (P,R) -> (side of P, mid vertex)
            for (size_t p = 0; p < patches.size(); ++p) {
                for (int s = 0; s < 4; ++s) {
                    const std::vector<int> &side = certs[p].sides[s];
                    if (side.size() < 2) continue;
                    const int e = C_.edgeBetween(side[0], side[1]);
                    if (e < 0 || C_.boundaryEdge[e]) continue;
                    const int q0 = C_.edgeCell[e][0], q1 = C_.edgeCell[e][1];
                    const int other = patchOf[q0] == static_cast<int>(p) ? patchOf[q1] : patchOf[q0];
                    if (other < 0 || other == static_cast<int>(p)) continue;
                    contacts[{static_cast<int>(p), other}].push_back({s, side[side.size() / 2]});
                }
            }
            for (const auto &kv : contacts) {
                if (kv.second.size() < 2 || kv.first.first > kv.first.second) continue;
                // Cut at the middle of the longest shared side, if it has an
                // interior vertex; else at the other.
                int bestMid = -1;
                size_t bestLen = 0;
                for (const auto &sm : kv.second) {
                    const size_t len = certs[kv.first.first].sides[sm.first].size();
                    if (len > 2 && len > bestLen && !forced[sm.second]) { bestLen = len; bestMid = sm.second; }
                }
                if (bestMid >= 0) cuts.push_back(bestMid);
            }
            if (cuts.empty()) {
                for (size_t p = 0; p < patches.size(); ++p) {
                    const int id = addCandidate(std::move(certs[p]));
                    if (id >= 0) candidates_[id].source = source;
                    cover.push_back(id);
                }
                // A candidate budget overrun leaves -1s; the cover is then unusable.
                if (std::find(cover.begin(), cover.end(), -1) != cover.end()) cover.clear();
                patchCount = static_cast<int>(patches.size());
                return;
            }
        }
        for (int v : cuts) {
            if (!forced[v]) { forced[v] = 1; ++report_.cutVertices; }
        }
    }
    report_.cutsConverged = false;
}

RectangleCertifier::RectangleCertifier(const SquareCarrier &carrier, const Options &opts)
    : C_(carrier), opts_(opts) {
    report_.cells = C_.numCells();
    const int NV = C_.numVertices();
    std::vector<char> hard(NV, 0), all(NV, 0);
    for (int v = 0; v < NV; ++v) {
        hard[v] = C_.forcedMacrovertex(v) ? 1 : 0;
        all[v] = hard[v] || C_.designatedVertex[v];
        report_.forcedVertices += hard[v];
        report_.softVertices += (!hard[v] && C_.designatedVertex[v]) ? 1 : 0;
    }

    int nAll = 0, nHard = 0;
    certifiedBase(all, baseAll_, "base (forced + designated)", nAll);
    report_.basePatchesAll = nAll;
    report_.fromBaseAll = static_cast<int>(candidates_.size());
    if (report_.softVertices > 0) {
        const int before = static_cast<int>(candidates_.size());
        certifiedBase(hard, baseHard_, "base (forced)", nHard);
        report_.basePatchesHard = nHard;
        report_.fromBaseHard = static_cast<int>(candidates_.size()) - before;
    } else {
        baseHard_ = baseAll_;
        report_.basePatchesHard = nAll;
    }

    const long long budget = report_.certifiedCells +
                             static_cast<long long>(opts_.workFactor * C_.numCells());

    // Merges: the union of two adjacent base patches, when it certifies.
    if (opts_.merges) {
        const int before = static_cast<int>(candidates_.size());
        const int NQ = C_.numCells();
        std::vector<int> patchOf(NQ, -1);
        const std::vector<int> base = baseAll_;
        for (size_t p = 0; p < base.size(); ++p) for (int q : candidates_[base[p]].cells) patchOf[q] = static_cast<int>(p);
        std::set<std::pair<int, int>> pairs;
        for (size_t p = 0; p < base.size(); ++p) {
            for (int q : candidates_[base[p]].cells) {
                for (int i = 0; i < 4; ++i) {
                    const int r = C_.neighbor[q][i];
                    if (r < 0 || patchOf[r] < 0 || patchOf[r] == static_cast<int>(p)) continue;
                    if (C_.interfaceEdge[C_.cellEdges[q][i]]) continue;
                    pairs.insert({std::min<int>(p, patchOf[r]), std::max<int>(p, patchOf[r])});
                }
            }
        }
        for (const auto &pr : pairs) {
            if (static_cast<int>(candidates_.size()) >= opts_.maxCandidates) break;
            if (report_.certifiedCells > budget) { report_.workBudgetHit = true; break; }
            std::vector<int> cells = candidates_[base[pr.first]].cells;
            const std::vector<int> &other = candidates_[base[pr.second]].cells;
            cells.insert(cells.end(), other.begin(), other.end());
            Certificate cert;
            Witness w;
            if (certify(cells, cert, w)) {
                const int id = addCandidate(std::move(cert));
                if (id >= 0 && candidates_[id].source.empty()) candidates_[id].source = "merge";
            }
        }
        report_.fromMerges = static_cast<int>(candidates_.size()) - before;
    }

    // Growth: boxes extended a full row at a time while they certify. Each
    // step tries the four sides in turn and takes the first that grows, so
    // the box that comes out is maximal -- no side of it can take another row.
    if (opts_.growth) {
        const int before = static_cast<int>(candidates_.size());
        // Largest patches first: they are the ones a grown box can matter for.
        std::vector<int> seeds = baseAll_;
        std::stable_sort(seeds.begin(), seeds.end(), [&](int a, int b) {
            return candidates_[a].cells.size() > candidates_[b].cells.size();
        });
        for (int sid : seeds) {
            if (static_cast<int>(candidates_.size()) >= opts_.maxCandidates) break;
            if (report_.certifiedCells > budget) { report_.workBudgetHit = true; break; }
            Certificate cur = candidates_[sid];
            bool grew = false;
            for (int step = 0; step < opts_.maxGrowthSteps && report_.certifiedCells <= budget; ++step) {
                bool any = false;
                const std::set<int> inCur(cur.cells.begin(), cur.cells.end());
                for (int s = 0; s < 4 && !any; ++s) {
                    const std::vector<int> &side = cur.sides[s];
                    std::vector<int> row;
                    bool ok = true;
                    for (size_t k = 0; k + 1 < side.size() && ok; ++k) {
                        const int e = C_.edgeBetween(side[k], side[k + 1]);
                        if (e < 0 || C_.boundaryEdge[e] || C_.interfaceEdge[e]) { ok = false; break; }
                        const int q0 = C_.edgeCell[e][0], q1 = C_.edgeCell[e][1];
                        const int out = inCur.count(q0) ? q1 : q0;
                        if (inCur.count(out)) ok = false;
                        else row.push_back(out);
                    }
                    if (!ok || row.empty()) continue;
                    std::vector<int> cells = cur.cells;
                    cells.insert(cells.end(), row.begin(), row.end());
                    Certificate next;
                    Witness w;
                    if (!certify(cells, next, w)) continue;
                    cur = std::move(next);
                    any = true;
                }
                if (!any) break;
                grew = true;
            }
            if (grew) {
                const int id = addCandidate(std::move(cur));
                if (id >= 0 && candidates_[id].source.empty()) candidates_[id].source = "growth";
            }
        }
        report_.fromGrowth = static_cast<int>(candidates_.size()) - before;
    }
    report_.candidates = static_cast<int>(candidates_.size());
}
