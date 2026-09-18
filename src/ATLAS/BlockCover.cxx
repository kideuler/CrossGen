#include "ATLAS/BlockCover.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <numeric>
#include <set>

typedef RectangleCertifier::Certificate Certificate;

RectangleCertifier::Certificate BlockCover::singleton(int q) const {
    Certificate c;
    c.cells = {q};
    c.G = {SquareTransport()};
    c.nu = c.nv = 1;
    for (int k = 0; k < 4; ++k) {
        c.corners[k] = C_.cells[q][k];
        c.sides[k] = {C_.cells[q][k], C_.cells[q][(k + 1) & 3]};
    }
    c.material = C_.cellMaterial[q];
    c.source = "singleton";
    RectangleCertifier::measure(C_, c);
    return c;
}

double BlockCover::cost(const Certificate &c) const {
    return 1.0 + opts_.alpha * c.distortion + opts_.beta * c.complexity;
}

BlockCover::Geo BlockCover::geometry(const Certificate &cert) const {
    Geo g;
    g.corners = cert.corners;
    for (int s = 0; s < 4; ++s) {
        const std::vector<int> &chain = cert.sides[s];
        for (size_t k = 0; k + 1 < chain.size(); ++k) g.sideEdges[s].push_back(C_.edgeBetween(chain[k], chain[k + 1]));
        for (size_t k = 1; k + 1 < chain.size(); ++k) g.sideInterior.insert(chain[k]);
    }
    return g;
}

// Sec. 6's incompatibility, as a test of the conforming property itself:
// every contact between P and R is one full side of each, at most one such
// side (singleContact), and no corner of either lies inside a side of the
// other.
template <class InP, class InR>
bool BlockCover::compatible(const Geo &P, InP inP, const Geo &R, InR inR) const {
    int fullSides = 0;
    for (int s = 0; s < 4; ++s) {
        int cnt = 0;
        for (int e : P.sideEdges[s]) {
            if (e < 0 || C_.boundaryEdge[e]) continue;
            const int q0 = C_.edgeCell[e][0], q1 = C_.edgeCell[e][1];
            const int other = inP(q0) ? q1 : q0;
            if (other >= 0 && inR(other)) ++cnt;
        }
        if (cnt == 0) continue;
        if (cnt != static_cast<int>(P.sideEdges[s].size())) return false;
        bool matched = false;
        for (int t = 0; t < 4 && !matched; ++t) {
            if (R.sideEdges[t].size() != P.sideEdges[s].size()) continue;
            std::vector<int> a = P.sideEdges[s], b = R.sideEdges[t];
            std::sort(a.begin(), a.end());
            std::sort(b.begin(), b.end());
            matched = a == b;
        }
        if (!matched) return false;
        ++fullSides;
    }
    if (opts_.singleContact && fullSides > 1) return false;
    for (int c : P.corners) if (R.sideInterior.count(c)) return false;
    for (int c : R.corners) if (P.sideInterior.count(c)) return false;
    return true;
}

// ---------------------------------------------------------------------------
// Sec. 13.2's gates on a complete selection.
// ---------------------------------------------------------------------------
bool BlockCover::analyze(const std::vector<Certificate> &sel, Report &rep,
                         std::vector<MacroEdge> *edgesOut, std::vector<int> *mvertsOut) const {
    const int NQ = C_.numCells(), NV = C_.numVertices();
    const int nb = static_cast<int>(sel.size());
    rep.blocks = nb;
    rep.uncovered = rep.overcovered = rep.sideMismatches = rep.hangingVertices = 0;
    rep.multipleContacts = rep.protectedNotMacro = rep.interfaceInside = 0;
    rep.singletonBlocks = rep.largestBlock = rep.mapPieces = 0;
    rep.objective = rep.meanDistortion = rep.maxDistortion = rep.meanComplexity = 0.0;

    std::vector<int> cellBlock(NQ, -1);
    for (int b = 0; b < nb; ++b) {
        for (int q : sel[b].cells) {
            if (cellBlock[q] >= 0) ++rep.overcovered;
            else cellBlock[q] = b;
        }
        const int k = static_cast<int>(sel[b].cells.size());
        rep.mapPieces += k;
        if (k == 1) ++rep.singletonBlocks;
        rep.largestBlock = std::max(rep.largestBlock, k);
        rep.objective += cost(sel[b]);
        rep.meanDistortion += sel[b].distortion;
        rep.maxDistortion = std::max(rep.maxDistortion, sel[b].distortion);
        rep.meanComplexity += sel[b].complexity;
    }
    if (nb > 0) { rep.meanDistortion /= nb; rep.meanComplexity /= nb; }
    for (int q = 0; q < NQ; ++q) if (cellBlock[q] < 0) ++rep.uncovered;

    std::vector<Geo> geo(nb);
    for (int b = 0; b < nb; ++b) geo[b] = geometry(sel[b]);
    std::vector<int> cornerCount(NV, 0);
    for (int b = 0; b < nb; ++b) for (int c : geo[b].corners) ++cornerCount[c];

    std::map<std::pair<int, int>, int> pairSides;
    std::vector<MacroEdge> edges;
    for (int b = 0; b < nb; ++b) {
        for (int s = 0; s < 4; ++s) {
            const std::vector<int> &se = geo[b].sideEdges[s];
            bool allBoundary = true, anyBoundary = false, allInterface = true;
            std::set<int> others;
            for (int e : se) {
                if (e < 0) { allBoundary = false; others.insert(-2); continue; }
                if (!C_.interfaceEdge[e]) allInterface = false;
                if (C_.boundaryEdge[e]) { anyBoundary = true; continue; }
                allBoundary = false;
                const int q0 = C_.edgeCell[e][0], q1 = C_.edgeCell[e][1];
                const int o = cellBlock[q0] == b ? cellBlock[q1] : cellBlock[q0];
                others.insert(o);
            }
            if (allBoundary) {
                MacroEdge me;
                me.chain = sel[b].sides[s];
                me.blockA = b;
                me.sideA = s;
                me.boundary = true;
                edges.push_back(me);
                continue;
            }
            if (anyBoundary || others.size() != 1 || *others.begin() < 0) { ++rep.sideMismatches; continue; }
            const int r = *others.begin();
            int match = -1;
            for (int t = 0; t < 4 && match < 0; ++t) {
                if (geo[r].sideEdges[t].size() != se.size()) continue;
                std::vector<int> a = se, bb = geo[r].sideEdges[t];
                std::sort(a.begin(), a.end());
                std::sort(bb.begin(), bb.end());
                if (a == bb) match = t;
            }
            if (match < 0) { ++rep.sideMismatches; continue; }
            ++pairSides[{std::min(b, r), std::max(b, r)}];
            if (b < r) {
                MacroEdge me;
                me.chain = sel[b].sides[s];
                me.blockA = b;
                me.sideA = s;
                me.blockB = r;
                me.sideB = match;
                me.interface = allInterface;
                edges.push_back(me);
            }
        }
    }
    for (const auto &kv : pairSides) {
        // Each shared side was counted from both blocks.
        if (kv.second / 2 > 1) ++rep.multipleContacts;
    }
    for (int b = 0; b < nb; ++b) {
        for (int v : geo[b].sideInterior) if (cornerCount[v] > 0) ++rep.hangingVertices;
    }
    for (int v = 0; v < NV; ++v) {
        if (C_.protectedVertex[v] && C_.valence[v] > 0 && cornerCount[v] == 0) ++rep.protectedNotMacro;
    }
    for (int b = 0; b < nb; ++b) {
        for (int q : sel[b].cells) {
            for (int i = 0; i < 4; ++i) {
                const int r = C_.neighbor[q][i];
                if (r >= 0 && cellBlock[r] == b && C_.interfaceEdge[C_.cellEdges[q][i]]) ++rep.interfaceInside;
            }
        }
    }
    rep.interfaceInside /= 2;

    // Sec. 4 on the macro complex: its vertices are the block corners.
    std::vector<int> mverts;
    long long lhs = 0;
    rep.irregularMacroVertices = 0;
    for (int v = 0; v < NV; ++v) {
        if (cornerCount[v] == 0) continue;
        mverts.push_back(v);
        const int q = cornerCount[v];
        if (C_.boundaryVertex[v]) {
            lhs += 2 - q;
            if (q != 2) ++rep.irregularMacroVertices;
        } else {
            lhs += 4 - q;
            if (q != 4) ++rep.irregularMacroVertices;
        }
    }
    long long rhs = 0;
    for (int c = 0; c < C_.componentCount; ++c) rhs += 4 * (1 - C_.componentHoles[c]);
    rep.eulerLHS = static_cast<int>(lhs);
    rep.eulerRHS = static_cast<int>(rhs);
    rep.eulerHolds = lhs == rhs;
    rep.macroVertices = static_cast<int>(mverts.size());
    rep.macroEdges = static_cast<int>(edges.size());

    rep.conforming = rep.uncovered == 0 && rep.overcovered == 0 && rep.sideMismatches == 0 &&
                     rep.hangingVertices == 0 && (!opts_.singleContact || rep.multipleContacts == 0);
    rep.valid = rep.conforming && rep.protectedNotMacro == 0 && rep.interfaceInside == 0 && rep.eulerHolds;
    if (edgesOut) *edgesOut = std::move(edges);
    if (mvertsOut) *mvertsOut = std::move(mverts);
    return rep.valid;
}

// ---------------------------------------------------------------------------
BlockCover::BlockCover(const RectangleCertifier &rects, const Options &opts)
    : R_(rects), C_(rects.getCarrier()), opts_(opts) {
    const std::vector<Certificate> &cands = R_.getCandidates();
    const int NQ = C_.numCells();

    // ---- incumbents ------------------------------------------------------
    std::vector<Certificate> best;
    double bestObj = 1e300;
    std::string bestName;
    auto consider = [&](std::vector<Certificate> sel, const std::string &name) {
        Report rep;
        if (!analyze(sel, rep, nullptr, nullptr)) {
            std::string why;
            if (rep.sideMismatches) why += " side mismatches " + std::to_string(rep.sideMismatches);
            if (rep.hangingVertices) why += " hanging " + std::to_string(rep.hangingVertices);
            if (rep.multipleContacts) why += " double contacts " + std::to_string(rep.multipleContacts);
            if (rep.protectedNotMacro) why += " protected " + std::to_string(rep.protectedNotMacro);
            if (rep.interfaceInside) why += " interface " + std::to_string(rep.interfaceInside);
            if (!rep.eulerHolds) why += " Euler";
            report_.messages.push_back("Stage 5: incumbent '" + name + "' fails the gates:" + why);
            return;
        }
        if (rep.objective < bestObj) {
            bestObj = rep.objective;
            best = std::move(sel);
            bestName = name;
        }
    };
    {
        std::vector<Certificate> singles;
        singles.reserve(NQ);
        for (int q = 0; q < NQ; ++q) singles.push_back(singleton(q));
        consider(std::move(singles), "singletons");
    }
    auto fromCover = [&](const std::vector<int> &ids) {
        std::vector<Certificate> sel;
        for (int id : ids) sel.push_back(cands[id]);
        return sel;
    };
    if (!R_.baseCoverAll().empty()) consider(fromCover(R_.baseCoverAll()), "base complex (forced + designated)");
    if (!R_.baseCoverHard().empty() && R_.baseCoverHard() != R_.baseCoverAll()) {
        consider(fromCover(R_.baseCoverHard()), "base complex (forced)");
    }
    if (best.empty()) {
        // Only possible if the carrier itself is not conforming, which Stage 2
        // would have refused; keep going with the singletons regardless.
        for (int q = 0; q < NQ; ++q) best.push_back(singleton(q));
        bestName = "singletons (unverified)";
    }
    report_.incumbent = bestName;

    // ---- local improvement ------------------------------------------------
    std::vector<Certificate> sel = std::move(best);
    std::vector<double> selCost(sel.size());
    std::vector<char> alive(sel.size(), 1);
    std::vector<Geo> geo(sel.size());
    std::vector<int> cellBlock(NQ, -1);
    for (size_t b = 0; b < sel.size(); ++b) {
        selCost[b] = cost(sel[b]);
        geo[b] = geometry(sel[b]);
        for (int q : sel[b].cells) cellBlock[q] = static_cast<int>(b);
    }
    std::vector<int> order(cands.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
        return cands[a].cells.size() > cands[b].cells.size();
    });
    std::vector<int> stamp(NQ, 0);
    int stampValue = 0;
    for (int sweep = 0; sweep < opts_.maxSweeps; ++sweep) {
        ++report_.sweeps;
        bool changed = false;
        for (int id : order) {
            const Certificate &P = cands[id];
            ++stampValue;
            for (int q : P.cells) stamp[q] = stampValue;
            std::vector<int> S;
            for (int q : P.cells) {
                const int b = cellBlock[q];
                if (std::find(S.begin(), S.end(), b) == S.end()) S.push_back(b);
            }
            if (S.size() == 1 && sel[S[0]].cells.size() == P.cells.size()) continue;   // already there
            bool inside = true;
            double oldCost = 0.0;
            for (int b : S) {
                oldCost += selCost[b];
                for (int q : sel[b].cells) if (stamp[q] != stampValue) { inside = false; break; }
                if (!inside) break;
            }
            if (!inside) continue;
            const double newCost = cost(P);
            if (newCost >= oldCost - 1e-9) continue;
            ++report_.swapTrials;
            const Geo gp = geometry(P);
            std::set<int> nbrs;
            for (int s = 0; s < 4; ++s) {
                for (int v : P.sides[s]) {
                    for (int r = C_.ringPtr[v]; r < C_.ringPtr[v + 1]; ++r) {
                        const int b = cellBlock[C_.ringCell[r]];
                        if (std::find(S.begin(), S.end(), b) == S.end()) nbrs.insert(b);
                    }
                }
            }
            auto inP = [&](int q) { return q >= 0 && stamp[q] == stampValue; };
            bool ok = true;
            for (int b : nbrs) {
                auto inR = [&](int q) { return q >= 0 && cellBlock[q] == b; };
                if (!compatible(gp, inP, geo[b], inR)) { ok = false; break; }
            }
            if (!ok) continue;
            for (int b : S) alive[b] = 0;
            const int nbId = static_cast<int>(sel.size());
            sel.push_back(P);
            selCost.push_back(newCost);
            alive.push_back(1);
            geo.push_back(gp);
            for (int q : P.cells) cellBlock[q] = nbId;
            ++report_.swaps;
            changed = true;
        }
        if (!changed) break;
    }

    std::vector<Certificate> finalSel;
    for (size_t b = 0; b < sel.size(); ++b) if (alive[b]) finalSel.push_back(std::move(sel[b]));
    const int candidates = static_cast<int>(cands.size());
    const std::string incumbent = report_.incumbent;
    const int swaps = report_.swaps, trials = report_.swapTrials, sweeps = report_.sweeps;
    std::vector<std::string> messages = report_.messages;
    analyze(finalSel, report_, &macroEdges_, &macroVertices_);
    report_.candidates = candidates;
    report_.incumbent = incumbent;
    report_.swaps = swaps;
    report_.swapTrials = trials;
    report_.sweeps = sweeps;
    report_.messages = messages;
    if (!report_.valid) report_.messages.push_back("Stage 5: the selected cover fails the Sec. 13.2 gates");
    blocks_.clear();
    for (Certificate &c : finalSel) {
        Block b;
        b.cost = cost(c);
        b.cert = std::move(c);
        blocks_.push_back(std::move(b));
    }
}

bool BlockCover::writeOBJ(const std::string &path) const {
    std::ofstream out(path);
    if (!out) return false;
    out.precision(17);
    out << "# ATLAS macro complex: " << blocks_.size() << " blocks, " << macroEdges_.size()
        << " macro edges, " << macroVertices_.size() << " macro vertices\n";
    for (const Point &p : C_.vertices) out << "v " << p[0] << " " << p[1] << " 0\n";
    for (const MacroEdge &me : macroEdges_) {
        out << "l";
        for (int v : me.chain) out << " " << v + 1;
        out << "\n";
    }
    return static_cast<bool>(out);
}
