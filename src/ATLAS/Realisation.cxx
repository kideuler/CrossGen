#include "ATLAS/Realisation.hxx"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <numeric>

#include <Eigen/Sparse>

#include "mesh/QuadMesh.hxx"
#include "mesh/TMOP.hxx"

namespace {

typedef SquareCarrier::Origin Origin;

std::vector<double> arcFractions(const std::vector<Point> &pts) {
    std::vector<double> f(pts.size(), 0.0);
    for (size_t k = 1; k < pts.size(); ++k) f[k] = f[k - 1] + normP(pts[k] - pts[k - 1]);
    const double L = f.empty() ? 0.0 : f.back();
    for (size_t k = 0; k < f.size(); ++k) {
        f[k] = L > 0.0 ? f[k] / L : static_cast<double>(k) / std::max<size_t>(1, f.size() - 1);
    }
    return f;
}

struct Quality {
    int inverted = 0;               // bad cells: a corner not positive, or at a bad vertex
    double minSJ = 1.0, meanSJ = 0.0;
    std::vector<char> bad;          // per cell
};

// The two local halves of Sec. 9.2 that validate() applies: every corner
// positive, and every one-ring winding exactly once (2 pi inside, the
// domain's own angle on dS). A mesh can pass the first and fail the second --
// a block folded over across dS has every corner positive and its boundary
// vertices at 3 pi -- so both mark cells bad, and the repairs below see both.
Quality measure(const std::vector<Point> &V, const std::vector<std::array<int, 4>> &cells,
                const std::vector<double> &want) {
    Quality Qy;
    Qy.bad.assign(cells.size(), 0);
    std::vector<double> angle(V.size(), 0.0);
    double sum = 0.0;
    for (size_t f = 0; f < cells.size(); ++f) {
        const auto &q = cells[f];
        double m = 1.0;
        for (int c = 0; c < 4; ++c) {
            const Point &p = V[q[c]];
            const Point u = V[q[(c + 1) & 3]] - p, w = V[q[(c + 3) & 3]] - p;
            const double d = normP(u) * normP(w);
            const double det = cross2(u, w);
            if (!(det > 0.0)) Qy.bad[f] = 1;
            m = std::min(m, d > 0.0 ? det / d : -1.0);
            double a = std::atan2(det, dotP(u, w));
            if (a < 0.0) a += 2.0 * M_PI;
            angle[q[c]] += a;
        }
        Qy.minSJ = std::min(Qy.minSJ, m);
        sum += m;
    }
    for (size_t f = 0; f < cells.size(); ++f) {
        if (Qy.bad[f]) continue;
        for (int v : cells[f]) {
            if (std::fabs(angle[v] - want[v]) > 1e-6) { Qy.bad[f] = 1; break; }
        }
    }
    for (char b : Qy.bad) Qy.inverted += b;
    Qy.meanSJ = cells.empty() ? 0.0 : sum / cells.size();
    return Qy;
}

} // namespace

Realisation::Realisation(const CoarseDomain &cd, const SquareCarrier &C, const BlockCover *cover,
                         const Options &opts) {
    const auto t0 = std::chrono::steady_clock::now();
    const PlanarDomain &F = cd.getFine();
    const Mesh &M = F.getMesh();
    const int NVc = C.numVertices(), NQc = C.numCells(), NEc = C.numEdges();
    const double h = opts.size > 0.0 ? opts.size : F.getReport().meanEdge;
    report_.coarseCells = NQc;
    auto fail = [&](const std::string &why) { report_.reason = why; };

    SquareCarrier::Quads &Q = quads_;
    auto addV = [&](const Point &p, Origin o, int sv, int se) {
        Q.vertices.push_back(p);
        Q.origin.push_back(o);
        Q.sourceVertex.push_back(sv);
        Q.sourceEdge.push_back(se);
        Q.designated.push_back(0);
        return static_cast<int>(Q.vertices.size()) - 1;
    };

    // ---- 1. Every coarse boundary vertex on the fine boundary -----------------
    std::vector<CoarseDomain::Location> loc(NVc);
    for (int v = 0; v < NVc; ++v) {
        if (!C.boundaryVertex[v]) continue;
        loc[v] = cd.locate(C, v);
        if (loc[v].loop < 0) { fail("a coarse boundary vertex has no place on the fine boundary"); return; }
    }
    // Snap onto input vertices: the exact samples first (they sit on one), then
    // the rest in loop order. A snap must keep the loop order -- a point
    // pulled past its neighbour would make the edge between them span the
    // whole loop -- and never puts two coarse vertices on one input vertex.
    std::vector<int> snappedTo(NVc, -1);
    std::vector<int> claimed(M.vertices.size(), -1);
    std::vector<double> eff(NVc, 0.0);   // effective arc position
    std::vector<std::vector<int>> onLoop(cd.arcs().size());
    for (int v = 0; v < NVc; ++v) {
        if (!C.boundaryVertex[v] || C.valence[v] == 0) continue;
        eff[v] = loc[v].s;
        onLoop[loc[v].loop].push_back(v);
    }
    for (auto &list : onLoop) {
        std::sort(list.begin(), list.end(), [&](int a, int b) { return loc[a].s < loc[b].s; });
    }
    for (int pass = 0; pass < 2; ++pass) {
        for (int l = 0; l < static_cast<int>(onLoop.size()); ++l) {
            const std::vector<int> &list = onLoop[l];
            const CoarseDomain::Arc &A = cd.arcs()[l];
            const int n = static_cast<int>(A.vertices.size()), m = static_cast<int>(list.size());
            const double Ltot = A.length;
            auto wrapHalf = [&](double d) {
                d = std::fmod(d, Ltot);
                if (d > 0.5 * Ltot) d -= Ltot;
                if (d <= -0.5 * Ltot) d += Ltot;
                return d;
            };
            for (int k = 0; k < m; ++k) {
                const int v = list[k];
                if ((C.sourceVertex[v] >= 0) != (pass == 0)) continue;
                int i = 0;
                cd.pointAt(l, loc[v].s, nullptr, &i);
                const double s0 = A.s[i], s1 = (i + 1 < n) ? A.s[i + 1] : A.length;
                const double seg = s1 - s0;
                const double d0 = loc[v].s - s0, d1 = s1 - loc[v].s;
                const int j = d0 <= d1 ? i : (i + 1) % n;
                const double dist = std::min(d0, d1);
                const int f = A.vertices[j];
                const bool exact = dist <= 1e-12 * std::max(1.0, A.length);
                // Offsets from v's own position: the candidate, and the two
                // neighbours along the loop where they are now.
                const double dc = wrapHalf(A.s[j] - loc[v].s);
                const double dp = m > 1 ? wrapHalf(eff[list[(k + m - 1) % m]] - loc[v].s) : -0.5 * Ltot;
                const double dn = m > 1 ? wrapHalf(eff[list[(k + 1) % m]] - loc[v].s) : 0.5 * Ltot;
                const bool ordered = exact || (dp < dc && dc < dn);
                if ((exact || dist <= opts.snap * seg) && ordered && claimed[f] < 0) {
                    claimed[f] = v;
                    snappedTo[v] = f;
                    eff[v] = A.s[j];
                    ++report_.snapped;
                } else if (exact) {
                    fail("two coarse boundary vertices land on one input vertex");
                    return;
                }
            }
        }
    }

    // ---- 2. Fine ids of the coarse vertices on block sides ---------------------
    // With a cover, only the macro complex is carried over vertex for vertex;
    // the inside of every block is interpolated afresh (step 5).
    std::vector<char> onMacroV(NVc, cover ? 0 : 1), onMacroE(NEc, cover ? 0 : 1);
    if (cover) {
        for (const BlockCover::Block &b : cover->getBlocks()) {
            for (const auto &side : b.cert.sides) {
                for (size_t k = 0; k < side.size(); ++k) {
                    onMacroV[side[k]] = 1;
                    if (k + 1 < side.size()) {
                        const int e = C.edgeBetween(side[k], side[k + 1]);
                        if (e < 0) { fail("a block side is not a chain of coarse edges"); return; }
                        onMacroE[e] = 1;
                    }
                }
            }
        }
        for (int e = 0; e < NEc; ++e) {
            if (C.boundaryEdge[e] && !onMacroE[e]) { fail("a boundary edge is on no block side"); return; }
        }
    }
    vertexMap_.assign(NVc, -1);
    std::vector<int> fineIdOfInput(M.vertices.size(), -1);
    for (int v = 0; v < NVc; ++v) {
        if (C.valence[v] == 0 || !onMacroV[v]) continue;
        if (C.boundaryVertex[v]) {
            if (snappedTo[v] >= 0) {
                const int f = snappedTo[v];
                vertexMap_[v] = addV(M.vertices[f], Origin::MeshVertex, f, -1);
                fineIdOfInput[f] = vertexMap_[v];
            } else {
                int e = -1;
                const Point p = cd.pointAt(loc[v].loop, loc[v].s, &e);
                vertexMap_[v] = addV(p, Origin::BoundarySplit, -1, e);
                ++report_.boundarySplits;
            }
        } else {
            vertexMap_[v] = addV(C.vertices[v], Origin::Realised, -1, -1);
        }
        if (C.designatedVertex[v]) Q.designated[vertexMap_[v]] = 1;
    }
    auto posOf = [&](int v) { return vertexMap_[v] >= 0 ? Q.vertices[vertexMap_[v]] : C.vertices[v]; };

    // ---- 3. Count classes (Sec. 11.2) ------------------------------------------
    std::vector<int> parent(NEc);
    std::iota(parent.begin(), parent.end(), 0);
    std::function<int(int)> find = [&](int x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };
    for (int q = 0; q < NQc; ++q) {
        for (int k = 0; k < 2; ++k) {
            const int a = find(C.cellEdges[q][k]), b = find(C.cellEdges[q][k + 2]);
            if (a != b) parent[a] = b;
        }
    }

    // The input vertices strictly inside each coarse boundary edge, in loop
    // order from the edge's start in loop order.
    std::vector<std::vector<int>> chain(NEc);
    std::vector<int> loopStart(NEc, -1);
    std::vector<double> length(NEc, 0.0);
    for (int e = 0; e < NEc; ++e) {
        const int a = C.edges[e][0], b = C.edges[e][1];
        if (!C.boundaryEdge[e]) {
            length[e] = normP(posOf(b) - posOf(a));
            continue;
        }
        const int q = C.edgeCell[e][0], i = C.edgeSide[e][0];
        const int from = C.cells[q][i], to = C.cells[q][(i + 1) & 3];
        loopStart[e] = from;
        const int l = loc[from].loop;
        if (loc[to].loop != l) { fail("a coarse boundary edge joins two boundary loops"); return; }
        const CoarseDomain::Arc &A = cd.arcs()[l];
        const int n = static_cast<int>(A.vertices.size());
        const double span = cd.forward(l, eff[from], eff[to]);
        length[e] = span;
        const double eps = 1e-12 * std::max(1.0, A.length);
        int i0 = 0;
        cd.pointAt(l, eff[from], nullptr, &i0);
        for (int k = 1; k <= n; ++k) {
            const int j = (i0 + k) % n;
            const double d = cd.forward(l, eff[from], A.s[j]);
            if (d >= span - eps || d >= A.length - eps) break;
            if (d <= eps) continue;
            chain[e].push_back(A.vertices[j]);
        }
    }
    std::vector<int> required(NEc, 1);
    std::vector<double> lenSum(NEc, 0.0);
    std::vector<int> members(NEc, 0);
    for (int e = 0; e < NEc; ++e) {
        const int r = find(e);
        if (C.boundaryEdge[e]) required[r] = std::max(required[r], static_cast<int>(chain[e].size()) + 1);
        lenSum[r] += length[e];
        ++members[r];
    }
    std::vector<int> count(NEc, 0);
    for (int e = 0; e < NEc; ++e) {
        if (find(e) != e) continue;
        ++report_.classes;
        const int want = static_cast<int>(std::lround(lenSum[e] / members[e] / h));
        count[e] = std::max({1, required[e], want});
    }

    // ---- 4. The fine points of every coarse edge on a block side ---------------
    std::vector<std::vector<int>> pts(NEc);
    for (int e = 0; e < NEc; ++e) {
        if (!onMacroE[e]) continue;
        const int n = count[find(e)];
        const int a = C.edges[e][0], b = C.edges[e][1];
        std::vector<int> L;
        if (C.boundaryEdge[e]) {
            const int from = loopStart[e];
            const int to = from == a ? b : a;
            const int l = loc[from].loop;
            const CoarseDomain::Arc &A = cd.arcs()[l];
            L.push_back(vertexMap_[from]);
            for (int f : chain[e]) {
                if (fineIdOfInput[f] < 0) fineIdOfInput[f] = addV(M.vertices[f], Origin::MeshVertex, f, -1);
                L.push_back(fineIdOfInput[f]);
            }
            L.push_back(vertexMap_[to]);
            // The input edge each piece lies on.
            auto pieceEdge = [&](int p, int q) {
                if (Q.origin[p] == Origin::BoundarySplit) return Q.sourceEdge[p];
                if (Q.origin[q] == Origin::BoundarySplit) return Q.sourceEdge[q];
                return A.edges[cd.arcIndex(Q.sourceVertex[p])];
            };
            int extra = n - (static_cast<int>(L.size()) - 1);
            while (extra-- > 0) {
                int best = 0;
                double bestLen = -1.0;
                for (size_t k = 0; k + 1 < L.size(); ++k) {
                    const double len = normP(Q.vertices[L[k + 1]] - Q.vertices[L[k]]);
                    if (len > bestLen) { bestLen = len; best = static_cast<int>(k); }
                }
                const int p = L[best], q = L[best + 1];
                const int id = addV((Q.vertices[p] + Q.vertices[q]) * 0.5, Origin::BoundarySplit, -1,
                                    pieceEdge(p, q));
                ++report_.boundarySplits;
                L.insert(L.begin() + best + 1, id);
            }
            if (from != a) std::reverse(L.begin(), L.end());
        } else {
            const Point pa = Q.vertices[vertexMap_[a]], pb = Q.vertices[vertexMap_[b]];
            L.push_back(vertexMap_[a]);
            for (int k = 1; k < n; ++k) {
                L.push_back(addV(pa + (pb - pa) * (static_cast<double>(k) / n), Origin::Realised, -1, -1));
            }
            L.push_back(vertexMap_[b]);
        }
        pts[e] = std::move(L);
    }
    // The fine points along a chain of coarse vertices, from its first to its last.
    auto chainPoints = [&](const std::vector<int> &cv, std::vector<int> &out) {
        out.clear();
        for (size_t k = 0; k + 1 < cv.size(); ++k) {
            const int e = C.edgeBetween(cv[k], cv[k + 1]);
            std::vector<int> part = pts[e];
            if (C.edges[e][0] != cv[k]) std::reverse(part.begin(), part.end());
            out.insert(out.end(), part.begin() + (k == 0 ? 0 : 1), part.end());
        }
    };

    // ---- 5. Every block (or, without a cover, every coarse cell) as one
    // transfinite grid. Its four sides are fine polylines; its inside,
    // including whatever coarse vertices and edges were there, is new.
    auto fillGrid = [&](const std::array<std::vector<int>, 4> &side, int material,
                        const std::function<int(int, int)> &groupOf) {
        const int n = static_cast<int>(side[0].size()) - 1, m = static_cast<int>(side[1].size()) - 1;
        if (n < 1 || m < 1 || static_cast<int>(side[2].size()) != n + 1 ||
            static_cast<int>(side[3].size()) != m + 1) {
            return false;
        }
        const int W = n + 1;
        std::vector<int> g(static_cast<size_t>(W) * (m + 1), -1);
        auto at = [&](int i, int j) -> int & { return g[i + W * j]; };
        for (int i = 0; i <= n; ++i) { at(i, 0) = side[0][i]; at(i, m) = side[2][n - i]; }
        for (int j = 0; j <= m; ++j) { at(n, j) = side[1][j]; at(0, j) = side[3][m - j]; }
        std::vector<Point> Bp(n + 1), Tp(n + 1), Lp(m + 1), Rp(m + 1);
        for (int i = 0; i <= n; ++i) { Bp[i] = Q.vertices[at(i, 0)]; Tp[i] = Q.vertices[at(i, m)]; }
        for (int j = 0; j <= m; ++j) { Lp[j] = Q.vertices[at(0, j)]; Rp[j] = Q.vertices[at(n, j)]; }
        const std::vector<double> fb = arcFractions(Bp), ft = arcFractions(Tp);
        const std::vector<double> fl = arcFractions(Lp), fr = arcFractions(Rp);
        const Point P00 = Bp[0], P10 = Bp[n], P11 = Tp[n], P01 = Tp[0];
        for (int j = 1; j < m; ++j) {
            for (int i = 1; i < n; ++i) {
                const double den = 1.0 - (ft[i] - fb[i]) * (fr[j] - fl[j]);
                const double u = (fb[i] + fl[j] * (ft[i] - fb[i])) / den;
                const double v = fl[j] + u * (fr[j] - fl[j]);
                const Point X = Bp[i] * (1.0 - v) + Tp[i] * v + Lp[j] * (1.0 - u) + Rp[j] * u -
                                (P00 * ((1.0 - u) * (1.0 - v)) + P10 * (u * (1.0 - v)) + P11 * (u * v) +
                                 P01 * ((1.0 - u) * v));
                at(i, j) = addV(X, Origin::Realised, -1, -1);
            }
        }
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < n; ++i) {
                Q.cells.push_back({at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)});
                Q.material.push_back(material);
                Q.group.push_back(groupOf(i, j));
            }
        }
        return true;
    };
    if (cover) {
        for (const BlockCover::Block &b : cover->getBlocks()) {
            const RectangleCertifier::Certificate &c = b.cert;
            std::array<std::vector<int>, 4> side;
            for (int k = 0; k < 4; ++k) chainPoints(c.sides[k], side[k]);
            // Fine column and row of every coarse column and row, for the
            // coarse cell each fine cell sits in (the local repair's regions).
            std::vector<int> colAt, rowAt;
            for (int i = 0; i < c.nu; ++i) {
                const int n = count[find(C.edgeBetween(c.sides[0][i], c.sides[0][i + 1]))];
                for (int k = 0; k < n; ++k) colAt.push_back(i);
            }
            for (int j = 0; j < c.nv; ++j) {
                const int n = count[find(C.edgeBetween(c.sides[1][j], c.sides[1][j + 1]))];
                for (int k = 0; k < n; ++k) rowAt.push_back(j);
            }
            auto groupOf = [&](int i, int j) {
                const int ci = i < static_cast<int>(colAt.size()) ? colAt[i] : c.nu - 1;
                const int cj = j < static_cast<int>(rowAt.size()) ? rowAt[j] : c.nv - 1;
                return c.cells[ci + c.nu * cj];
            };
            if (!fillGrid(side, c.material, groupOf)) { fail("opposite sides of a block got different counts"); return; }
            ++report_.blocks;
        }
    } else {
        for (int q = 0; q < NQc; ++q) {
            std::array<std::vector<int>, 4> side;
            for (int i = 0; i < 4; ++i) {
                const int e = C.cellEdges[q][i];
                side[i] = pts[e];
                if (C.edges[e][0] != C.cells[q][i]) std::reverse(side[i].begin(), side[i].end());
            }
            if (!fillGrid(side, C.cellMaterial[q], [q](int, int) { return q; })) {
                fail("opposite sides of a coarse cell got different counts");
                return;
            }
        }
    }
    report_.vertices = static_cast<int>(Q.vertices.size());
    report_.cells = static_cast<int>(Q.cells.size());

    // ---- 6. Untangle and smooth the interior -----------------------------------
    std::vector<double> want(Q.vertices.size(), 2.0 * M_PI);
    for (size_t v = 0; v < Q.vertices.size(); ++v) {
        if (Q.origin[v] == Origin::MeshVertex) want[v] = F.interiorAngle[Q.sourceVertex[v]];
        else if (Q.origin[v] == Origin::BoundarySplit) want[v] = M_PI;
    }
    const Quality before = measure(Q.vertices, Q.cells, want);
    report_.invertedInterpolated = before.inverted;
    report_.minScaledJacobianInterpolated = before.minSJ;
    Quality best = before;
    // Where a curve bulges into a thin coarse cell the interpolation folds.
    // The discrete harmonic map with dS held fixed opens such a fold -- every
    // interior node at the mean of its neighbours, so the ones next to the
    // curve are pulled off it -- whereas TMOP's untangler, started on a folded
    // mesh, knows nothing of the domain and can drag nodes across a hole. So
    // the harmonic state is tried first, and the better of it and the
    // interpolation (by inverted cells, then worst scaled Jacobian) is what
    // TMOP starts from. It is solved outright, not by sweeps: a fold a few
    // dozen cells deep is out of reach of any affordable number of Jacobi
    // iterations.
    if (before.inverted > 0 && opts.smoothingSweeps > 0) {
        const int NV = static_cast<int>(Q.vertices.size());
        std::vector<std::vector<int>> nbr(NV);
        for (const auto &q : Q.cells) {
            for (int c = 0; c < 4; ++c) {
                nbr[q[c]].push_back(q[(c + 1) & 3]);
                nbr[q[c]].push_back(q[(c + 3) & 3]);
            }
        }
        for (auto &l : nbr) { std::sort(l.begin(), l.end()); l.erase(std::unique(l.begin(), l.end()), l.end()); }
        std::vector<int> unknown(NV, -1);
        int nu = 0;
        for (int v = 0; v < NV; ++v) if (Q.origin[v] == Origin::Realised && !nbr[v].empty()) unknown[v] = nu++;
        if (nu > 0) {
            std::vector<Eigen::Triplet<double>> T;
            Eigen::VectorXd bx = Eigen::VectorXd::Zero(nu), by = Eigen::VectorXd::Zero(nu);
            for (int v = 0; v < NV; ++v) {
                const int i = unknown[v];
                if (i < 0) continue;
                T.emplace_back(i, i, static_cast<double>(nbr[v].size()));
                for (int w : nbr[v]) {
                    if (unknown[w] >= 0) {
                        T.emplace_back(i, unknown[w], -1.0);
                    } else {
                        bx[i] += Q.vertices[w][0];
                        by[i] += Q.vertices[w][1];
                    }
                }
            }
            Eigen::SparseMatrix<double> L(nu, nu);
            L.setFromTriplets(T.begin(), T.end());
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(L);
            if (solver.info() == Eigen::Success) {
                const Eigen::VectorXd x = solver.solve(bx), y = solver.solve(by);
                std::vector<Point> harmonic = Q.vertices;
                for (int v = 0; v < NV; ++v) {
                    if (unknown[v] >= 0) harmonic[v] = {x[unknown[v]], y[unknown[v]]};
                }
                const Quality qy = measure(harmonic, Q.cells, want);
                report_.invertedHarmonic = qy.inverted;
                if (qy.inverted < best.inverted || (qy.inverted == best.inverted && qy.minSJ > best.minSJ)) {
                    best = qy;
                    Q.vertices = harmonic;
                    report_.harmonic = true;
                }
            }
        }
    }
    // What is left is local: a thin coarse cell along a curve, whose interior
    // side the coarse search put inside the lens between the chord and the
    // arc. A harmonic re-solve over just the coarse cells round each fold --
    // everything outside that region held where it is -- spreads the rows
    // across the region and opens the fold, where the global solve above
    // disturbed the whole mesh to do it. The region grows a ring at a time.
    if (best.inverted > 0) {
        const int NV = static_cast<int>(Q.vertices.size());
        const int NF = static_cast<int>(Q.cells.size());
        std::vector<std::vector<int>> cellsAt(NV);
        for (int f = 0; f < NF; ++f) for (int v : Q.cells[f]) cellsAt[v].push_back(f);
        // Coarse cell adjacency through shared vertices.
        std::vector<std::vector<int>> coarseNbr(NQc);
        for (int v = 0; v < NVc; ++v) {
            for (int r = C.ringPtr[v]; r < C.ringPtr[v + 1]; ++r) {
                for (int s = C.ringPtr[v]; s < C.ringPtr[v + 1]; ++s) {
                    if (r != s) coarseNbr[C.ringCell[r]].push_back(C.ringCell[s]);
                }
            }
        }
        for (int rings = 1; rings <= 4 && best.inverted > 0; ++rings) {
            std::vector<char> region(NQc, 0);
            for (int f = 0; f < NF; ++f) if (best.bad[f]) region[Q.group[f]] = 1;
            for (int g = 0; g < rings; ++g) {
                std::vector<char> next = region;
                for (int q = 0; q < NQc; ++q) if (region[q]) for (int r : coarseNbr[q]) next[r] = 1;
                region.swap(next);
            }
            std::vector<int> unknown(NV, -1);
            int nu = 0;
            for (int v = 0; v < NV; ++v) {
                if (Q.origin[v] != Origin::Realised || cellsAt[v].empty()) continue;
                bool inside = true;
                for (int f : cellsAt[v]) if (!region[Q.group[f]]) { inside = false; break; }
                if (inside) unknown[v] = nu++;
            }
            if (nu == 0) continue;
            std::vector<Eigen::Triplet<double>> T;
            Eigen::VectorXd bx = Eigen::VectorXd::Zero(nu), by = Eigen::VectorXd::Zero(nu);
            for (int v = 0; v < NV; ++v) {
                const int i = unknown[v];
                if (i < 0) continue;
                std::vector<int> nb;
                for (int f : cellsAt[v]) {
                    const auto &q = Q.cells[f];
                    for (int c = 0; c < 4; ++c) {
                        if (q[c] != v) continue;
                        nb.push_back(q[(c + 1) & 3]);
                        nb.push_back(q[(c + 3) & 3]);
                    }
                }
                std::sort(nb.begin(), nb.end());
                nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
                T.emplace_back(i, i, static_cast<double>(nb.size()));
                for (int w : nb) {
                    if (unknown[w] >= 0) {
                        T.emplace_back(i, unknown[w], -1.0);
                    } else {
                        bx[i] += Q.vertices[w][0];
                        by[i] += Q.vertices[w][1];
                    }
                }
            }
            Eigen::SparseMatrix<double> L(nu, nu);
            L.setFromTriplets(T.begin(), T.end());
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(L);
            if (solver.info() != Eigen::Success) continue;
            const Eigen::VectorXd x = solver.solve(bx), y = solver.solve(by);
            std::vector<Point> trial = Q.vertices;
            for (int v = 0; v < NV; ++v) if (unknown[v] >= 0) trial[v] = {x[unknown[v]], y[unknown[v]]};
            const Quality qy = measure(trial, Q.cells, want);
            if (qy.inverted < best.inverted || (qy.inverted == best.inverted && qy.minSJ > best.minSJ)) {
                best = qy;
                Q.vertices = trial;
                ++report_.localRepairs;
            }
        }
    }
    if (opts.smoothingSweeps > 0) {
        mesh::QuadMesh::Options qo;
        qo.fixAllFeatureNodes = true;
        mesh::QuadMesh qm(Q.vertices, Q.cells, Q.material, qo);
        mesh::TMOP::Options to;
        to.metric = mesh::TMOP::Shape002;
        to.maxSweeps = opts.smoothingSweeps;
        to.untangle = true;
        mesh::TMOP smoother(qm, to);
        smoother.run();
        report_.tmopSweeps = smoother.getReport().sweeps;
        report_.untangleSweeps = smoother.getReport().untangleSweeps;
        if (qm.vertices.size() == Q.vertices.size()) {
            const Quality after = measure(qm.vertices, Q.cells, want);
            // An untangling that leaves folds is not trusted (see above); TMOP
            // is kept only when it ends valid and no worse.
            const bool better = after.inverted == 0 && (best.inverted > 0 || after.minSJ > best.minSJ);
            if (better) {
                // Only interior nodes may have moved; put dS back exactly, so
                // that validate() audits the input's own coordinates.
                for (size_t v = 0; v < Q.vertices.size(); ++v) {
                    if (Q.origin[v] == Origin::Realised) Q.vertices[v] = qm.vertices[v];
                }
                best = measure(Q.vertices, Q.cells, want);
                report_.smoothed = true;
            }
        }
    }
    report_.inverted = best.inverted;
    report_.minScaledJacobian = best.minSJ;
    report_.meanScaledJacobian = best.meanSJ;
    report_.built = true;
    report_.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}
