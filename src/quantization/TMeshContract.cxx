#include "quantization/TMeshContract.hxx"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <set>
#include <vector>

namespace {

double polylineLength(const std::vector<Point> &pts) {
    double len = 0.0;
    for (size_t i = 1; i < pts.size(); ++i) len += normP(pts[i] - pts[i - 1]);
    return len;
}

Point pointAt(const std::vector<Point> &pts, double target) {
    double acc = 0.0;
    for (size_t i = 1; i < pts.size(); ++i) {
        const double seg = normP(pts[i] - pts[i - 1]);
        if (acc + seg >= target && seg > 0.0) {
            return pts[i - 1] + (pts[i] - pts[i - 1]) * ((target - acc) / seg);
        }
        acc += seg;
    }
    return pts.back();
}

// The piece of `pts` between two arc lengths, endpoints included.
std::vector<Point> subCurve(const std::vector<Point> &pts, double from, double to) {
    std::vector<Point> out{pointAt(pts, from)};
    double acc = 0.0;
    for (size_t i = 1; i < pts.size(); ++i) {
        const double seg = normP(pts[i] - pts[i - 1]);
        acc += seg;
        if (acc > from && acc < to) out.push_back(pts[i]);
    }
    out.push_back(pointAt(pts, to));
    return out;
}

// The curve halfway between two polylines that run the same way, sampled
// evenly along both. This is the line the two surviving neighbours of a
// removed cell meet on.
std::vector<Point> midCurve(const std::vector<Point> &a,
                            const std::vector<Point> &b) {
    const size_t n = std::max<size_t>(std::max(a.size(), b.size()), 2);
    const double la = polylineLength(a), lb = polylineLength(b);
    std::vector<Point> out;
    out.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        const double t = static_cast<double>(i) / (n - 1);
        out.push_back((pointAt(a, la * t) + pointAt(b, lb * t)) * 0.5);
    }
    return out;
}

// The working T-mesh, mutable while the cleanup runs.
struct Work {
    struct Edge {
        int x = 0;
        double xIdeal = 1.0;
        std::vector<Point> geom;
        bool alive = true;
    };
    struct Face {
        std::array<std::vector<int>, 4> sides;
        // Whether each side entry's geometry runs against the walk.
        std::array<std::vector<char>, 4> rev;
        int block = -1;
        bool alive = true;
    };
    std::vector<Edge> edges;
    std::vector<Face> faces;

    long long sideSum(int f, int s) const {
        long long sum = 0;
        for (int e : faces[f].sides[s]) sum += edges[e].x;
        return sum;
    }

    // Cut edge `e` so that its geometry's first part keeps `keep` quads and
    // a new edge takes the rest. Every face holding `e` is updated, on the
    // correct side of it for the direction that face walks.
    int split(int e, int keep) {
        Edge &E = edges[e];
        const int rest = E.x - keep;
        const double len = polylineLength(E.geom);
        const double at = len * keep / E.x;

        Edge tail;
        tail.x = rest;
        tail.xIdeal = std::max(1e-9, E.xIdeal * rest / E.x);
        tail.geom = subCurve(E.geom, at, len);

        const std::vector<Point> head = subCurve(E.geom, 0.0, at);
        E.xIdeal = std::max(1e-9, E.xIdeal * keep / E.x);
        E.x = keep;
        E.geom = head;

        edges.push_back(std::move(tail));
        const int e2 = static_cast<int>(edges.size()) - 1;

        for (Face &f : faces) {
            if (!f.alive) continue;
            for (int s = 0; s < 4; ++s) {
                for (size_t i = 0; i < f.sides[s].size(); ++i) {
                    if (f.sides[s][i] != e) continue;
                    // Walking with the geometry the head comes first;
                    // against it, the tail does.
                    const bool against = f.rev[s][i] != 0;
                    f.sides[s].insert(f.sides[s].begin() + i + (against ? 0 : 1), e2);
                    f.rev[s].insert(f.rev[s].begin() + i + (against ? 0 : 1),
                                    against ? 1 : 0);
                    ++i;  // step over the entry just inserted
                }
            }
        }
        return e2;
    }
};

// Union-find over edges that also tracks orientation.
//
// Identifying two edges is not enough: they were traced independently and
// may run opposite ways, so the survivor's curve can point against what the
// faces of the other one expect. Carrying a parity bit to the class
// representative lets every face's "runs against the walk" flag be
// corrected when the edges are folded together. Without it the sides of a
// merged cell come out reversed and its patch turns inside out.
struct PDSU {
    std::vector<int> parent;
    std::vector<char> par;  // orientation relative to the parent

    explicit PDSU(size_t n) : parent(n), par(n, 0) {
        std::iota(parent.begin(), parent.end(), 0);
    }
    // Root of a's class, and whether a is flipped with respect to it.
    std::pair<int, int> find(int a) const {
        int p = 0;
        while (parent[a] != a) {
            p ^= par[a];
            a = parent[a];
        }
        return {a, p};
    }
    // `opposite` says a and b run against each other.
    void join(int a, int b, int opposite) {
        const auto [ra, pa] = find(a);
        const auto [rb, pb] = find(b);
        if (ra == rb) return;
        parent[rb] = ra;
        par[rb] = static_cast<char>(pa ^ pb ^ opposite);
    }
    void grow(size_t n) {
        while (parent.size() < n) {
            parent.push_back(static_cast<int>(parent.size()));
            par.push_back(0);
        }
    }
};

// The quad offsets at which a side is cut into edges, interior breaks only.
std::vector<int> breaks(const Work &w, const std::vector<int> &side) {
    std::vector<int> out;
    int acc = 0;
    for (size_t i = 0; i + 1 < side.size(); ++i) {
        acc += w.edges[side[i]].x;
        out.push_back(acc);
    }
    return out;
}

// Give both sides the same interior breaks, splitting edges as needed.
// Both are walked in the direction they face each other across the cell.
int refine(Work &w, std::vector<int> &A, std::vector<char> &Arev,
           std::vector<int> &B, std::vector<char> &Brev, int face, int sideA,
           int sideB) {
    int splits = 0;
    for (int guard = 0; guard < 64; ++guard) {
        std::vector<int> ba = breaks(w, A), bb = breaks(w, B);
        std::vector<int> want;
        std::set_union(ba.begin(), ba.end(), bb.begin(), bb.end(),
                       std::back_inserter(want));
        if (ba == want && bb == want) break;

        bool did = false;
        for (int target : want) {
            for (int pass = 0; pass < 2 && !did; ++pass) {
                std::vector<int> &S = pass == 0 ? A : B;
                const std::vector<int> &have = pass == 0 ? ba : bb;
                if (std::find(have.begin(), have.end(), target) != have.end()) continue;
                // Find the edge straddling `target` and cut it there.
                int acc = 0;
                for (size_t i = 0; i < S.size(); ++i) {
                    const int x = w.edges[S[i]].x;
                    if (acc < target && target < acc + x) {
                        const std::vector<char> &rv = pass == 0 ? Arev : Brev;
                        // `keep` counts along the geometry, the walk may
                        // run the other way.
                        const int fromStart = target - acc;
                        w.split(S[i], rv[i] ? x - fromStart : fromStart);
                        ++splits;
                        did = true;
                        break;
                    }
                    acc += x;
                }
            }
            if (did) break;
        }
        if (!did) break;
        // The split rewrote the face's side lists; pick them up again.
        A = w.faces[face].sides[sideA];
        Arev.assign(w.faces[face].rev[sideA].begin(), w.faces[face].rev[sideA].end());
        B = w.faces[face].sides[sideB];
        Brev.assign(w.faces[face].rev[sideB].begin(), w.faces[face].rev[sideB].end());
        std::reverse(B.begin(), B.end());
        std::reverse(Brev.begin(), Brev.end());
        for (char &c : Brev) c = c ? 0 : 1;
    }
    return splits;
}

}  // namespace

ContractReport contractZeroEdges(BlockQuant &bq) {
    ContractReport report;
    QuantTMesh &tm = bq.tmesh;

    Work w;
    w.edges.resize(tm.edges.size());
    for (size_t e = 0; e < tm.edges.size(); ++e) {
        w.edges[e].x = tm.edges[e].x;
        w.edges[e].xIdeal = tm.edges[e].xIdeal;
        w.edges[e].geom = bq.edgeGeometry[e];
    }
    w.faces.resize(tm.faces.size());
    for (size_t f = 0; f < tm.faces.size(); ++f) {
        for (int s = 0; s < 4; ++s) {
            w.faces[f].sides[s] = tm.faces[f].sides[s];
            w.faces[f].rev[s] = bq.sideReversed[f][s];
            w.faces[f].rev[s].resize(w.faces[f].sides[s].size(), 0);
        }
        w.faces[f].block = bq.blockOfFace[f];
    }

    // Which cells collapsed, and along which axis they still have length.
    std::vector<int> collapsed;   // face ids with exactly one live axis
    for (size_t f = 0; f < w.faces.size(); ++f) {
        const bool zero0 = w.sideSum(static_cast<int>(f), 0) == 0;
        const bool zero1 = w.sideSum(static_cast<int>(f), 1) == 0;
        if (!zero0 && !zero1) continue;
        if (zero0 && zero1) {
            w.faces[f].alive = false;   // the cell is a single point
            ++report.pointCells;
            continue;
        }
        collapsed.push_back(static_cast<int>(f));
    }

    // Reads a collapsed cell's two surviving sides, the second walked
    // backwards so the pair line up across the vanished direction.
    auto liveSides = [&](int f, std::vector<int> &A, std::vector<char> &Arev,
                         std::vector<int> &B, std::vector<char> &Brev) {
        const int live = w.sideSum(f, 0) == 0 ? 1 : 0;
        A.assign(w.faces[f].sides[live].begin(), w.faces[f].sides[live].end());
        Arev.assign(w.faces[f].rev[live].begin(), w.faces[f].rev[live].end());
        B.assign(w.faces[f].sides[live + 2].begin(),
                 w.faces[f].sides[live + 2].end());
        Brev.assign(w.faces[f].rev[live + 2].begin(),
                    w.faces[f].rev[live + 2].end());
        std::reverse(B.begin(), B.end());
        std::reverse(Brev.begin(), Brev.end());
        for (char &c : Brev) c = c ? 0 : 1;
        return live;
    };

    // First bring every collapsed cell's two sides into step. Splitting
    // rewrites side lists, so it all happens before any merging: a split
    // landing on an edge that had already been folded into another cell's
    // class would leave that class describing the wrong curve.
    for (int f : collapsed) {
        std::vector<int> A, B;
        std::vector<char> Arev, Brev;
        const int live = liveSides(f, A, Arev, B, Brev);
        report.splitEdges += refine(w, A, Arev, B, Brev, f, live, live + 2);
    }

    PDSU dsu(w.edges.size());
    for (int f : collapsed) {
        std::vector<int> A, B;
        std::vector<char> Arev, Brev;
        liveSides(f, A, Arev, B, Brev);

        bool pairable = A.size() == B.size();
        for (size_t i = 0; pairable && i < A.size(); ++i) {
            if (w.edges[A[i]].x != w.edges[B[i]].x) pairable = false;
        }
        if (!pairable) {
            // The sides could not be brought into step. Leaving the cell
            // alone keeps the structure sound, at the cost of the zero
            // edges it holds.
            ++report.remainingZero;
            continue;
        }

        for (size_t i = 0; i < A.size(); ++i) {
            // Both curves walked the same way, then averaged: the two
            // neighbours of the vanishing cell meet along the middle of it
            // and between them cover the ground it stood on.
            std::vector<Point> ga = w.edges[A[i]].geom, gb = w.edges[B[i]].geom;
            if (Arev[i]) std::reverse(ga.begin(), ga.end());
            if (Brev[i]) std::reverse(gb.begin(), gb.end());
            std::vector<Point> mid = midCurve(ga, gb);

            dsu.join(A[i], B[i], Arev[i] != Brev[i] ? 1 : 0);

            // `mid` runs along the walk; store it in the representative's
            // own direction so the parity bits describe it correctly.
            const auto [root, flipped] = dsu.find(A[i]);
            if ((Arev[i] != 0) != (flipped != 0)) {
                std::reverse(mid.begin(), mid.end());
            }
            w.edges[root].geom = std::move(mid);
        }
        w.faces[f].alive = false;
        ++report.mergedCells;
    }

    // Rebuild: keep the live faces, map every edge to its class, and drop
    // the zero-length ones from the sides they sat on.
    dsu.grow(w.edges.size());
    std::vector<int> newId(w.edges.size(), -1);
    BlockQuant out;
    out.nodes = bq.nodes;
    out.faceOfBlock.assign(bq.faceOfBlock.size(), -1);
    out.skippedBlocks = bq.skippedBlocks;
    out.skippedNonQuad = bq.skippedNonQuad;
    out.skippedUnlocatable = bq.skippedUnlocatable;
    out.skippedDegenerate = bq.skippedDegenerate;
    out.skippedOverlapping = bq.skippedOverlapping;

    for (size_t f = 0; f < w.faces.size(); ++f) {
        if (!w.faces[f].alive) continue;
        std::array<std::vector<int>, 4> sides;
        std::array<std::vector<char>, 4> rev;
        bool usable = true;
        for (int s = 0; s < 4; ++s) {
            for (size_t i = 0; i < w.faces[f].sides[s].size(); ++i) {
                const int e = w.faces[f].sides[s][i];
                if (w.edges[e].x == 0) continue;   // contracted away
                const auto [rep, flipped] = dsu.find(e);
                if (newId[rep] < 0) {
                    newId[rep] = out.tmesh.addEdge(w.edges[rep].xIdeal);
                    out.tmesh.edges[newId[rep]].x = w.edges[rep].x;
                    out.edgeGeometry.push_back(w.edges[rep].geom);
                    out.edgeNodes.push_back({-1, -1});
                }
                sides[s].push_back(newId[rep]);
                // The face's flag was written against the edge it used to
                // hold; the surviving curve may run the other way, and the
                // parity bit is exactly that difference.
                rev[s].push_back(
                    static_cast<char>((w.faces[f].rev[s][i] != 0) != (flipped != 0)));
            }
            if (sides[s].empty()) usable = false;
        }
        if (!usable) {
            ++report.remainingZero;
            continue;
        }
        const int nf = out.tmesh.addFace(sides[0], sides[1], sides[2], sides[3]);
        out.sideReversed.push_back(rev);
        out.blockOfFace.push_back(w.faces[f].block);
        if (w.faces[f].block >= 0 &&
            w.faces[f].block < static_cast<int>(out.faceOfBlock.size())) {
            out.faceOfBlock[w.faces[f].block] = nf;
        }
    }

    // Rebuild the node list: contraction merged some and the mid-curves
    // moved others, so the old one no longer describes this T-mesh.
    {
        out.nodes.clear();
        double diag = 0.0;
        Point lo{1e300, 1e300}, hi{-1e300, -1e300};
        for (const auto &g : out.edgeGeometry) {
            for (const Point &p : g) {
                lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
                hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
            }
        }
        if (!out.edgeGeometry.empty()) diag = normP(hi - lo);
        const double tol = std::max(1e-12, 1e-9 * diag);
        auto nodeAt = [&](const Point &p) {
            for (size_t i = 0; i < out.nodes.size(); ++i) {
                if (normP(out.nodes[i] - p) < tol) return static_cast<int>(i);
            }
            out.nodes.push_back(p);
            return static_cast<int>(out.nodes.size()) - 1;
        };
        for (size_t e = 0; e < out.edgeGeometry.size(); ++e) {
            out.edgeNodes[e] = {nodeAt(out.edgeGeometry[e].front()),
                                nodeAt(out.edgeGeometry[e].back())};
        }
    }

    report.removedEdges =
        static_cast<int>(tm.edges.size()) - static_cast<int>(out.tmesh.edges.size());
    report.ok = out.tmesh.finalize(&out.error);
    report.error = out.error;
    if (report.ok) {
        out.ok = true;
        bq = std::move(out);
    }
    return report;
}
