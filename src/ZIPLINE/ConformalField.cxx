#include "ZIPLINE/ConformalField.hxx"

#include "MERIDIAN/Interfaces.hxx"

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <algorithm>
#include <functional>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <queue>
#include <sstream>

namespace {

// The representative of a mod pi/2 nearest ref.
inline double liftNear(double ref, double a) { return a + M_PI_2 * std::round((ref - a) / M_PI_2); }
// Into (-pi/4, pi/4]: how far a cross angle is from zero.
inline double wrapQuarter(double a) { return a - M_PI_2 * std::round(a / M_PI_2); }

// Table 1 of Viertel et al. with its own intervals, exactly as CrossField::
// initialize sets the boundary data and SeparatrixTrace::buildBoundary reads the
// corners back: the corners the second field turns at have to be the ones the
// first was told to turn at and the tracing will launch from.
int quartersOf(double interiorAngle) {
    if (interiorAngle < 0.75 * M_PI) return 1;
    if (interiorAngle <= 1.25 * M_PI) return 2;
    if (interiorAngle <= 1.75 * M_PI) return 3;
    return std::max(4, static_cast<int>(std::lround(interiorAngle / M_PI_2)));
}

inline Point midpoint(const Mesh &m, int t, int k) {
    const Triangle &tri = m.triangles[t];
    return (m.vertices[tri[k]] + m.vertices[tri[(k + 1) % 3]]) * 0.5;
}

} // namespace

ConformalField::ConformalField(std::shared_ptr<Mesh> mesh, std::vector<Singular> pattern,
                               const Options &opts, const CrossField *first)
    : mesh_(std::move(mesh)), opts_(opts), first_(first), sing_(std::move(pattern)) {
    build();
    first_ = nullptr;   // read during construction only
}

ConformalField::ConformalField(std::shared_ptr<Mesh> mesh, std::vector<Singular> pattern,
                               const Options &opts, const CrossField *first, const Interfaces *interfaces,
                               std::vector<int> absorbed4, std::vector<std::array<int, 3>> absorbedTriangles)
    : mesh_(std::move(mesh)), opts_(opts), first_(first), sing_(std::move(pattern)) {
    if (interfaces && interfaces->multiMaterial()) {
        interfaces_ = interfaces;
        absorbed4_ = std::move(absorbed4);
        absorbedTriangles_ = std::move(absorbedTriangles);
    }
    build();
    first_ = nullptr;
    absorbed4_.clear();
    absorbedTriangles_.clear();
}

double ConformalField::singularPart(const Point &x) const {
    double s = 0.0;
    for (const Singular &q : sing_) s += 0.25 * q.d * std::atan2(x[1] - q.pos[1], x[0] - q.pos[0]);
    for (const BoundarySingular &q : bsing_) s += 0.5 * q.d * std::atan2(x[1] - q.pos[1], x[0] - q.pos[0]);
    return s;
}

double ConformalField::singularPart(const Point &x, int region) const {
    double s = 0.0;
    for (const Singular &q : sing_) s += 0.25 * q.d * std::atan2(x[1] - q.pos[1], x[0] - q.pos[0]);
    for (const BoundarySingular &q : bsing_)
        if (q.region < 0 || q.region == region) s += 0.5 * q.d * std::atan2(x[1] - q.pos[1], x[0] - q.pos[0]);
    return s;
}

double ConformalField::angleAt(int t, const Point &x, double ref) const {
    const Mesh &m = *mesh_;
    const Triangle &tri = m.triangles[t];
    const Point &a = m.vertices[tri[0]], &b = m.vertices[tri[1]], &c = m.vertices[tri[2]];
    const double den = cross2(b - a, c - a);
    double l1 = 1.0 / 3.0, l2 = 1.0 / 3.0;
    if (std::fabs(den) > 1e-300) {
        l1 = cross2(x - a, c - a) / den;
        l2 = cross2(b - a, x - a) / den;
    }
    if (regionTerms_) {
        const double r0 = regularCorner_[t][0];
        const double r1 = liftNear(r0, regularCorner_[t][1]);
        const double r2 = liftNear(r0, regularCorner_[t][2]);
        return liftNear(ref, (1.0 - l1 - l2) * r0 + l1 * r1 + l2 * r2 + singularPart(x, regionOf_[t]));
    }
    const double r0 = regular_[tri[0]];
    const double r1 = liftNear(r0, regular_[tri[1]]);
    const double r2 = liftNear(r0, regular_[tri[2]]);
    return liftNear(ref, (1.0 - l1 - l2) * r0 + l1 * r1 + l2 * r2 + singularPart(x));
}

// ---------------------------------------------------------------------------
// build()
// ---------------------------------------------------------------------------
void ConformalField::build() {
    report_ = Report{};
    nodeQuarters_.clear();
    splitCorners_.clear();
    offNetwork_.clear();
    auto fail = [&](const std::string &why) { report_.built = false; report_.failure = why; };
    if (mesh_->triangles.empty()) return fail("no triangles");

    // ---- with interfaces: the mesh cut along them (see the header) ---------
    //
    // Every corner of every triangle is joined to the same vertex's corner in
    // the neighbour across each of its two edges at that vertex, unless the
    // edge is an interface edge; each class is one vertex of the cut mesh.
    // Triangles keep their indices, so a singularity's triangle is the same in
    // both. fixedQ holds, per cut vertex, the quarter count an interface node
    // prescribes for its sector, -1 where Table 1 decides.
    std::shared_ptr<Mesh> cutMesh;
    std::vector<int> origOf;
    std::vector<int> fixedQ;
    std::vector<int> landingFan;   // per uncut vertex: Table 1's quarters of a landing node's whole corner
    const int nVo = static_cast<int>(mesh_->vertices.size());
    if (interfaces_) {
        const Mesh &O = *mesh_;
        const int nTo = static_cast<int>(O.triangles.size());
        std::vector<char> onInterface(O.edges.size(), 0);
        for (const int e : interfaces_->interfaceEdges())
            if (e >= 0 && e < static_cast<int>(O.edges.size())) onInterface[e] = 1;
        std::vector<int> parent(3 * static_cast<size_t>(nTo));
        for (size_t i = 0; i < parent.size(); ++i) parent[i] = static_cast<int>(i);
        std::function<int(int)> find = [&](int x) {
            while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
            return x;
        };
        for (int t = 0; t < nTo; ++t)
            for (int k = 0; k < 3; ++k) {
                const int v = O.triangles[t][k];
                // Local edge i runs from corner i to corner i + 1: the two at
                // corner k are k and k + 2.
                for (const int le : {k, (k + 2) % 3}) {
                    const int e = O.triangleEdges[t][le];
                    const int n = O.triangleAdjacency[t][le];
                    if (n < 0 || (e >= 0 && onInterface[e])) continue;
                    for (int kk = 0; kk < 3; ++kk)
                        if (O.triangles[n][kk] == v) {
                            const int a = find(3 * t + k), b = find(3 * n + kk);
                            if (a != b) parent[a] = b;
                        }
                }
            }
        std::vector<int> idOfRoot(parent.size(), -1);
        std::vector<Point> V;
        std::vector<Triangle> T(nTo);
        for (int t = 0; t < nTo; ++t)
            for (int k = 0; k < 3; ++k) {
                int &id = idOfRoot[find(3 * t + k)];
                if (id < 0) {
                    id = static_cast<int>(V.size());
                    V.push_back(O.vertices[O.triangles[t][k]]);
                    origOf.push_back(O.triangles[t][k]);
                }
                T[t][k] = id;
            }
        cutMesh = std::make_shared<Mesh>(V, T, O.triangleMatId);

        // Each node's sector quarters. The node's total is the fan's own count
        // -- 4 inside, Table 1's reading of the whole corner on dS -- less what
        // the node absorbed, as SeparatrixTrace::launchFromInterfaceNodes has
        // it. How the total is shared among the sectors is read off the first
        // field, sector by sector, since each sector is a corner of a
        // different region and each region's Poincare-Hopf count is the first
        // field's: a sector of angle s across which the field turns by Delta,
        // counter-clockwise round the node, is a corner of
        //
        //     q = round((s - Delta) / (pi/2))
        //
        // quarter turns (the corner model z^(q pi / 2s), whose cross turns by
        // s - q pi/2 from one side of the corner to the other). Delta is read
        // on the ring of vertices round the node, never at the node, whose
        // value is the first field's compromise between the branches pinned
        // there. Shared out by the sectors' angles alone, the way the first
        // field's trace does it, a 270-degree corner an interface bisects
        // reads 2 + 1 where the field turned 1 + 2, and the two regions then
        // fail Poincare-Hopf by a quarter each, in opposite directions (ten
        // of the 31 multi-material models). Where the readings do not add up
        // to the total they are adjusted by largest remainder.
        fixedQ.assign(V.size(), -1);
        std::vector<double> fullAngle(nVo, 0.0);
        for (const Triangle &tri : O.triangles)
            for (int k = 0; k < 3; ++k) {
                const Point a = O.vertices[tri[(k + 1) % 3]] - O.vertices[tri[k]];
                const Point b = O.vertices[tri[(k + 2) % 3]] - O.vertices[tri[k]];
                const double na = normP(a), nb = normP(b);
                if (na > 1e-300 && nb > 1e-300)
                    fullAngle[tri[k]] += std::acos(std::max(-1.0, std::min(1.0, dotP(a, b) / (na * nb))));
            }
        // The first field's turn across sector k of node n, counter-clockwise
        // from ray k to ray k + 1 on the ring of the node's neighbours; NaN
        // when the walk does not close the sector.
        auto sectorTurn = [&](const Interfaces::Node &n, size_t k) {
            const double nan = std::numeric_limits<double>::quiet_NaN();
            if (!first_ || static_cast<int>(first_->u_k.size()) != nVo) return nan;
            const size_t R = n.rays.size();
            const int v = n.vertex;
            const int endAt = n.rays[(k + 1) % R].neighbour;
            auto th = [&](int u) { return 0.25 * std::arg(first_->u_k[u]); };
            int t = n.rays[k].faceCCW;
            double turn = 0.0;
            for (size_t guard = 0; guard < 64 && t >= 0; ++guard) {
                int x = -1, y = -1;
                for (int kk = 0; kk < 3; ++kk) {
                    if (O.triangles[t][kk] != v) continue;
                    x = O.triangles[t][(kk + 1) % 3];
                    y = O.triangles[t][(kk + 2) % 3];
                }
                if (x < 0) return nan;
                if (cross2(O.vertices[x] - O.vertices[v], O.vertices[y] - O.vertices[v]) < 0.0) std::swap(x, y);
                turn += wrapQuarter(th(y) - th(x));
                if (y == endAt) return turn;
                // On round the node, across the edge v-y, unless the sector
                // ends there.
                int next = -1;
                for (int le = 0; le < 3; ++le) {
                    const int a = O.triangles[t][le], b = O.triangles[t][(le + 1) % 3];
                    if (!((a == v && b == y) || (a == y && b == v))) continue;
                    const int e = O.triangleEdges[t][le];
                    if (e >= 0 && onInterface[e]) return nan;
                    next = O.triangleAdjacency[t][le];
                }
                t = next;
            }
            return nan;
        };
        const auto &nodes = interfaces_->nodes();
        nodeQuarters_.assign(nodes.size(), {});
        landingFan.assign(nVo, 2);
        for (size_t ni = 0; ni < nodes.size(); ++ni) {
            const Interfaces::Node &n = nodes[ni];
            if (n.kind == Interfaces::NodeKind::Dangling || n.vertex < 0 || n.sector.empty()) continue;
            const int fan = n.onBoundary ? quartersOf(fullAngle[n.vertex]) : 4;
            if (n.onBoundary) landingFan[n.vertex] = fan;
            const int absorbed = n.vertex < static_cast<int>(absorbed4_.size()) ? absorbed4_[n.vertex] : 0;
            report_.absorbed4 += absorbed;
            const int total = std::max(static_cast<int>(n.sector.size()), fan - absorbed);
            // The reading, or the sector's own quarter turns where there is none.
            std::vector<double> want(n.sector.size());
            for (size_t k = 0; k < n.sector.size(); ++k) {
                const double turn = (k < n.rays.size()) ? sectorTurn(n, k) : std::numeric_limits<double>::quiet_NaN();
                want[k] = std::isfinite(turn) ? (n.sector[k] - turn) / M_PI_2 : n.sector[k] / M_PI_2;
            }
            // Rounded, then made to add up to the node's total by largest
            // remainder; then a sector read as no quad at all (a valence-2
            // corner: a sharp sector the first field did not turn, the
            // interface version of a cusp) is given one, P2's split (Sec.
            // 6.3), *without* taking it from the next sector: that region is
            // then a valence-3 point short, which Stage 5b supplies, and the
            // others are left exactly as the first field turned them.
            std::vector<int> r(n.sector.size(), 0);
            int used = 0;
            for (size_t k = 0; k < r.size(); ++k) {
                r[k] = static_cast<int>(std::lround(want[k]));
                used += r[k];
            }
            const int nodeTotal = fan - absorbed;
            while (used < nodeTotal) {
                size_t best = 0;
                double bestGap = -1e300;
                for (size_t k = 0; k < r.size(); ++k)
                    if (want[k] - r[k] > bestGap) { bestGap = want[k] - r[k]; best = k; }
                ++r[best];
                ++used;
            }
            while (used > nodeTotal) {
                size_t best = 0;
                double bestGap = 1e300;
                for (size_t k = 0; k < r.size(); ++k)
                    if (want[k] - r[k] < bestGap) { bestGap = want[k] - r[k]; best = k; }
                --r[best];
                --used;
            }
            std::vector<int> q(r.size(), 1);
            for (size_t k = 0; k < r.size(); ++k) {
                q[k] = std::max(1, r[k]);
                if (r[k] < 1) {
                    report_.splitSectors += 1 - r[k];
                    if (k < n.rays.size()) splitCorners_.push_back(SplitCorner{n.vertex, n.rays[k].faceCCW});
                }
            }
            for (size_t k = 0; k < n.sector.size() && k < n.rays.size(); ++k) {
                const int t = n.rays[k].faceCCW;
                if (t < 0 || t >= nTo) continue;
                for (int kk = 0; kk < 3; ++kk)
                    if (O.triangles[t][kk] == n.vertex) fixedQ[T[t][kk]] = q[k];
            }
            nodeQuarters_[ni] = std::move(q);
            ++report_.interfaceNodes;
        }
    }
    const Mesh &M = cutMesh ? *cutMesh : *mesh_;
    const int nV = static_cast<int>(M.vertices.size());
    const int nT = static_cast<int>(M.triangles.size());
    auto orig = [&](int v) { return cutMesh ? origOf[v] : v; };

    double meanEdge = 0.0;
    for (const auto &e : M.edges) meanEdge += normP(M.vertices[e[1]] - M.vertices[e[0]]);
    meanEdge = M.edges.empty() ? 1.0 : meanEdge / static_cast<double>(M.edges.size());

    std::vector<char> used(nV, 0);
    for (const Triangle &t : M.triangles) used[t[0]] = used[t[1]] = used[t[2]] = 1;

    // ---- the boundary, as loops with the domain on their left --------------
    //
    // Read off the geometry rather than the vertex order of the owning
    // triangle, so nothing below depends on the mesh being counter-clockwise.
    std::vector<int> nextOf(nV, -1), prevOf(nV, -1), edgeAfter(nV, -1);
    for (const int e : M.boundaryEdges) {
        const int t = (M.edgeTriangles[e][0] >= 0) ? M.edgeTriangles[e][0] : M.edgeTriangles[e][1];
        if (t < 0) continue;
        int a = M.edges[e][0], b = M.edges[e][1];
        int c = -1;
        for (int k = 0; k < 3; ++k)
            if (M.triangles[t][k] != a && M.triangles[t][k] != b) c = M.triangles[t][k];
        if (c < 0) continue;
        if (cross2(M.vertices[b] - M.vertices[a], M.vertices[c] - M.vertices[a]) < 0.0) std::swap(a, b);
        if (nextOf[a] >= 0 || prevOf[b] >= 0)
            return fail("the boundary is not a disjoint union of simple loops (a pinched vertex)");
        nextOf[a] = b;
        prevOf[b] = a;
        edgeAfter[a] = e;
    }
    std::vector<std::vector<int>> loops;
    std::vector<int> loopOf(nV, -1);
    for (int v = 0; v < nV; ++v) {
        if (nextOf[v] < 0 || loopOf[v] >= 0) continue;
        std::vector<int> loop;
        int u = v;
        while (u >= 0 && loopOf[u] < 0) {
            loopOf[u] = static_cast<int>(loops.size());
            loop.push_back(u);
            u = nextOf[u];
        }
        if (u != v) return fail("a boundary loop does not close");
        loops.push_back(std::move(loop));
    }
    report_.loops = static_cast<int>(loops.size());

    // Interior angle at each boundary vertex, and from it the turning angle
    // tau_b = pi - alpha_b (Sec. 2.3) and the Table 1 quarters q_b.
    std::vector<double> interior(nV, 0.0);
    for (const Triangle &tri : M.triangles) {
        for (int k = 0; k < 3; ++k) {
            const int v = tri[k];
            if (loopOf[v] < 0) continue;
            const Point a = M.vertices[tri[(k + 1) % 3]] - M.vertices[v];
            const Point b = M.vertices[tri[(k + 2) % 3]] - M.vertices[v];
            const double na = normP(a), nb = normP(b);
            if (na < 1e-300 || nb < 1e-300) continue;
            interior[v] += std::acos(std::max(-1.0, std::min(1.0, dotP(a, b) / (na * nb))));
        }
    }
    std::vector<int> quarters(nV, 2);
    std::vector<double> loopTurning(loops.size(), 0.0);
    for (size_t l = 0; l < loops.size(); ++l)
        for (const int b : loops[l]) {
            quarters[b] = (!fixedQ.empty() && fixedQ[b] >= 0) ? fixedQ[b] : quartersOf(interior[b]);
            if (quarters[b] != 2) ++report_.corners;
            loopTurning[l] += M_PI - interior[b];
        }

    // ---- Options::readWeakBoundary: every boundary quarter off the first field
    //
    // The vertex reading, the quarters a misalignment adds inside an edge,
    // and the absorbed triangles, as the option has them. An interface edge
    // is on two loops, once each way: both see the same misalignment with
    // the opposite sign and give it to the same end (the nearer, the one
    // of the lower vertex number on a tie), so the two sides of a vertex
    // still add up to what the uncut fan turns.
    const bool weak = opts_.readWeakBoundary && first_ && static_cast<int>(first_->u_k.size()) == nVo;
    if (weak) {
        auto th = [&](int v) { return 0.25 * std::arg(first_->u_k[orig(v)]); };
        std::vector<int> wq(nV, 2);
        for (const auto &loop : loops)
            for (const int b : loop) {
                const double phIn = computeAngle(M.vertices[b] - M.vertices[prevOf[b]]);
                const double phOut = computeAngle(M.vertices[nextOf[b]] - M.vertices[b]);
                const double mIn = wrapQuarter(th(b) - phIn), mOut = wrapQuarter(th(b) - phOut);
                wq[b] = static_cast<int>(std::lround((interior[b] + mIn - mOut) / M_PI_2));
            }
        for (const auto &loop : loops)
            for (const int a : loop) {
                const int b = nextOf[a];
                const double ph = computeAngle(M.vertices[b] - M.vertices[a]);
                const double ma = wrapQuarter(th(a) - ph), mb = wrapQuarter(th(b) - ph);
                const double D = wrapQuarter(th(b) - th(a));
                const int j = static_cast<int>(std::lround((ma + D - mb) / M_PI_2));
                if (j == 0) continue;
                // Where along a -> b the misalignment passes +-45 degrees.
                const double t = std::fabs(D) > 1e-12 ? ((j > 0 ? M_PI_4 : -M_PI_4) - ma) / D : 0.5;
                const bool toA = std::fabs(t - 0.5) < 1e-9 ? orig(a) < orig(b) : t < 0.5;
                wq[toA ? a : b] += j;
                ++report_.weakWraps;
            }
        for (const auto &ab : absorbedTriangles_) {
            const int t = ab[0], v = ab[1], d = ab[2];
            if (t < 0 || t >= nT) continue;
            for (int k = 0; k < 3; ++k)
                if (orig(M.triangles[t][k]) == v) {
                    wq[M.triangles[t][k]] -= d;
                    break;
                }
        }
        std::vector<char> isNodeVertex(nVo, 0);
        if (interfaces_)
            for (const auto &n : interfaces_->nodes())
                if (n.kind != Interfaces::NodeKind::Dangling && n.vertex >= 0 && n.vertex < nVo)
                    isNodeVertex[n.vertex] = 1;
        std::vector<int> aTriangle(nV, -1);
        for (int t = 0; t < nT; ++t)
            for (int k = 0; k < 3; ++k) aTriangle[M.triangles[t][k]] = t;
        report_.splitSectors = 0;
        splitCorners_.clear();
        std::vector<char> listed(nVo, 0);
        for (const auto &loop : loops)
            for (const int b : loop) {
                const int o = orig(b);
                const bool onDS = mesh_->isBoundaryVertex[o];
                const int strong = (!fixedQ.empty() && fixedQ[b] >= 0) ? fixedQ[b] : quartersOf(interior[b]);
                if (wq[b] != strong) ++report_.weakCorners;
                // No corner flatter than a quad can have. A corner of q
                // quarters over an angle of more than q times 135 degrees --
                // Table 1's own boundary between one quarter and two -- is a
                // valence-3 point the weak field let out: through a convex
                // corner, read as no quarter at all (P2's valence-2 corner,
                // Sec. 6.3), or through a smooth stretch, read as a corner of
                // one on a straight side (cruciform, dam: a worst element of
                // 0 after smoothing). Each is put back in and the corner given
                // its quarter, which leaves the region's index sum as it was:
                // at a convex corner in a triangle of its fan, which Sec.
                // 12.3's relocation below moves into the patch; elsewhere a
                // slip length (two local edges at least) inside, along the
                // inward normal, where the strong condition would have held
                // it off -- nearer if that is out of the region.
                const int qMin = std::max(1, static_cast<int>(std::ceil(interior[b] / (0.75 * M_PI) - 1e-9)));
                if (wq[b] < qMin) {
                    const int k = qMin - wq[b];
                    int t = aTriangle[b];
                    Point x{0.0, 0.0};
                    if (t >= 0) {
                        const Triangle &tr = M.triangles[t];
                        x = (M.vertices[tr[0]] + M.vertices[tr[1]] + M.vertices[tr[2]]) / 3.0;
                    }
                    if (t >= 0 && qMin > 1) {
                        const Point &X = M.vertices[b];
                        const Point eIn = normalizeP(X - M.vertices[prevOf[b]]);
                        const Point eOut = normalizeP(M.vertices[nextOf[b]] - X);
                        // Left of travel is inside (the loops keep the region on their left).
                        const Point n = normalizeP(Point{-(eIn[1] + eOut[1]), eIn[0] + eOut[0]});
                        const double local = 0.5 * (normP(M.vertices[nextOf[b]] - X) + normP(M.vertices[prevOf[b]] - X));
                        const double slip = onDS ? first_->weakBoundaryLength() : first_->weakInterfaceLength();
                        for (double dist = std::max(slip, 2.0 * local); dist > 0.25 * local; dist *= 0.5) {
                            const Point y = X + n * dist;
                            // A walk from the fan towards y; the cut mesh has no
                            // neighbour across an interface, so it stays in the region.
                            int u = aTriangle[b];
                            for (int guard = 0; guard < 4096 && u >= 0; ++guard) {
                                const Triangle &tu = M.triangles[u];
                                const Point &p0 = M.vertices[tu[0]], &p1 = M.vertices[tu[1]], &p2 = M.vertices[tu[2]];
                                const double den = cross2(p1 - p0, p2 - p0);
                                const double l1 = cross2(y - p0, p2 - p0) / den, l2 = cross2(p1 - p0, y - p0) / den;
                                const double lb[3] = {1.0 - l1 - l2, l1, l2};
                                int worst = -1;
                                double w = -1e-12;
                                for (int kk = 0; kk < 3; ++kk) if (lb[kk] < w) { w = lb[kk]; worst = kk; }
                                if (worst < 0) break;
                                u = M.triangleAdjacency[u][(worst + 1) % 3];
                            }
                            if (u >= 0) { t = u; x = y; break; }
                        }
                    }
                    if (t >= 0) {
                        for (int kk = 0; kk < k; ++kk) sing_.push_back(Singular{x, t, 1});
                        report_.weakReinserted += k;
                    } else {
                        report_.splitSectors += k;
                        splitCorners_.push_back(SplitCorner{o, t});
                    }
                    wq[b] = qMin;
                }
                if (onDS && !isNodeVertex[o] && first_->weakBoundary() && wq[b] != quartersOf(interior[b])) {
                    bsing_.push_back(BoundarySingular{M.vertices[b], b, quartersOf(interior[b]) - wq[b]});
                    ++report_.weakBoundaryPoints;
                }
                if (!onDS && !isNodeVertex[o] && wq[b] != 2 && !listed[o]) {
                    listed[o] = 1;
                    offNetwork_.push_back(o);
                }
                quarters[b] = wq[b];
            }
        report_.offNetwork = static_cast<int>(offNetwork_.size());
        // The node sectors are these readings too (a sector is a corner of
        // its region), and a landing's fan their sum.
        if (interfaces_) {
            const auto &nodes = interfaces_->nodes();
            for (size_t ni = 0; ni < nodes.size() && ni < nodeQuarters_.size(); ++ni) {
                const Interfaces::Node &n = nodes[ni];
                if (nodeQuarters_[ni].empty()) continue;
                int fan = 0;
                for (size_t k = 0; k < nodeQuarters_[ni].size() && k < n.rays.size(); ++k) {
                    const int t = n.rays[k].faceCCW;
                    if (t < 0 || t >= nT) continue;
                    for (int kk = 0; kk < 3; ++kk)
                        if (orig(M.triangles[t][kk]) == n.vertex) nodeQuarters_[ni][k] = wq[M.triangles[t][kk]];
                    fan += nodeQuarters_[ni][k];
                }
                if (n.onBoundary && n.vertex >= 0 && n.vertex < static_cast<int>(landingFan.size()))
                    landingFan[n.vertex] = fan;
            }
        }
        report_.corners = 0;
        for (const auto &loop : loops)
            for (const int b : loop) if (quarters[b] != 2) ++report_.corners;
        if (!offNetwork_.empty()) return fail("interface corners off the network");
    }

    // ---- components, and the outer loop of each -----------------------------
    std::vector<int> compOf(nT, -1);
    int nComp = 0;
    for (int t0 = 0; t0 < nT; ++t0) {
        if (compOf[t0] >= 0) continue;
        std::queue<int> q;
        q.push(t0);
        compOf[t0] = nComp;
        while (!q.empty()) {
            const int t = q.front();
            q.pop();
            for (int k = 0; k < 3; ++k) {
                const int n = M.triangleAdjacency[t][k];
                if (n >= 0 && compOf[n] < 0) { compOf[n] = nComp; q.push(n); }
            }
        }
        ++nComp;
    }
    report_.components = nComp;
    regionOf_ = compOf;
    std::vector<int> compOfVertex(nV, -1);
    for (int t = 0; t < nT; ++t)
        for (int k = 0; k < 3; ++k) compOfVertex[M.triangles[t][k]] = compOf[t];
    // ---- Options::interfaceCornerTerms: a node sector's own term --------------
    //
    // At every sector of an interface node (a boundary vertex of a region
    // whose original is no vertex of dS) turned by a quarter count other than
    // Table 1's reading of its angle, the BoundarySingular term of that
    // difference, in that region only.
    if (interfaces_ && (weak || opts_.interfaceCornerTerms)) {
        for (const auto &loop : loops)
            for (const int b : loop) {
                const int o = orig(b);
                if (mesh_->isBoundaryVertex[o]) continue;
                const int d = quartersOf(interior[b]) - quarters[b];
                if (d == 0) continue;
                bsing_.push_back(BoundarySingular{M.vertices[b], b, d, compOfVertex[b]});
                ++report_.cornerTerms;
            }
        regionTerms_ = report_.cornerTerms > 0;
    }
    std::vector<std::vector<int>> loopsOfComp(nComp);
    for (size_t l = 0; l < loops.size(); ++l) loopsOfComp[compOfVertex[loops[l][0]]].push_back(static_cast<int>(l));
    std::vector<int> outerOf(nComp, -1);
    for (int c = 0; c < nComp; ++c) {
        for (const int l : loopsOfComp[c]) {
            // +2 pi round the outer loop, -2 pi round a hole (Sec. 2.3).
            if (loopTurning[l] > 0.0) {
                if (outerOf[c] >= 0) return fail("a component with two outer loops");
                outerOf[c] = l;
            } else {
                ++report_.holes;
            }
        }
        if (outerOf[c] < 0) return fail("a component with no outer loop");
    }

    // ---- Sec. 5.4's automatic edits -----------------------------------------
    //
    // Locate each singularity's triangle afresh from its position if the one
    // given does not hold it, merge singularities sharing a triangle (their
    // local models cannot both be the one the sweep uses), then cancel dipoles.
    auto bary = [&](int t, const Point &x) {
        const Triangle &tri = M.triangles[t];
        const Point &a = M.vertices[tri[0]], &b = M.vertices[tri[1]], &c = M.vertices[tri[2]];
        const double den = cross2(b - a, c - a);
        const double l1 = cross2(x - a, c - a) / den, l2 = cross2(b - a, x - a) / den;
        return std::array<double, 3>{{1.0 - l1 - l2, l1, l2}};
    };
    for (Singular &s : sing_) {
        if (s.triangle < 0 || s.triangle >= nT) return fail("a singularity with no triangle");
        const auto l = bary(s.triangle, s.pos);
        if (std::min({l[0], l[1], l[2]}) < -1e-9) return fail("a singularity outside its triangle");
    }
    {
        std::vector<Singular> kept;
        std::vector<int> slot(nT, -1);
        for (const Singular &s : sing_) {
            if (s.d == 0) continue;
            int &k = slot[s.triangle];
            if (k < 0) { k = static_cast<int>(kept.size()); kept.push_back(s); continue; }
            // Two in one triangle: one singularity of the summed index, at
            // the index-weighted mean (or simply the mean when they cancel).
            Singular &o = kept[k];
            const int wsum = std::abs(o.d) + std::abs(s.d);
            o.pos = (o.pos * static_cast<double>(std::abs(o.d)) + s.pos * static_cast<double>(std::abs(s.d))) /
                    static_cast<double>(wsum);
            o.d += s.d;
            ++report_.merged;
        }
        sing_.clear();
        for (const Singular &s : kept)
            if (s.d != 0) sing_.push_back(s);
    }
    if (opts_.cancelDipoles && opts_.dipoleRadius > 0.0) {
        const double R = opts_.dipoleRadius * meanEdge;
        struct Pair { double dist; int a, b; };
        std::vector<Pair> pairs;
        for (size_t i = 0; i < sing_.size(); ++i)
            for (size_t j = i + 1; j < sing_.size(); ++j) {
                if (sing_[i].d + sing_[j].d != 0) continue;
                // Across an interface a +- pair is not a dipole: each point is
                // in its own region's index count (concrete: an aggregate's two
                // valence-5 points facing two valence-3 points of the ring
                // round it, cancelled, left both regions failing
                // Poincare-Hopf by two).
                if (compOf[sing_[i].triangle] != compOf[sing_[j].triangle]) continue;
                const double d = normP(sing_[i].pos - sing_[j].pos);
                if (d < R) pairs.push_back({d, static_cast<int>(i), static_cast<int>(j)});
            }
        std::sort(pairs.begin(), pairs.end(), [](const Pair &x, const Pair &y) { return x.dist < y.dist; });
        std::vector<char> gone(sing_.size(), 0);
        for (const Pair &p : pairs) {
            if (gone[p.a] || gone[p.b]) continue;
            gone[p.a] = gone[p.b] = 1;
            ++report_.dipolesCancelled;
        }
        std::vector<Singular> kept;
        for (size_t i = 0; i < sing_.size(); ++i)
            if (!gone[i]) kept.push_back(sing_[i]);
        sing_.swap(kept);
    }
    // Sec. 12.3: a +1/4 singularity hugging a convex corner, moved into the
    // patch (Options::relocateCornerSingularities).
    if (opts_.relocateCornerSingularities && opts_.cornerRelocation > 0.0) {
        // Boundary length from b to the next corner along the loop, either way.
        auto runTo = [&](int b, bool forward) {
            double L = 0.0;
            int v = b;
            for (size_t guard = 0; guard < M.vertices.size(); ++guard) {
                const int n = forward ? nextOf[v] : prevOf[v];
                if (n < 0) break;
                L += normP(M.vertices[n] - M.vertices[v]);
                v = n;
                if (v == b || quarters[v] != 2) break;
            }
            return L;
        };
        auto locate = [&](const Point &x, int hint) {
            // A walk from the hint towards x, then a scan if it fails.
            int t = hint;
            for (int guard = 0; guard < 4096 && t >= 0; ++guard) {
                const auto l = bary(t, x);
                int worst = -1;
                double w = -1e-12;
                for (int k = 0; k < 3; ++k) if (l[k] < w) { w = l[k]; worst = k; }
                if (worst < 0) return t;
                // The edge opposite the most negative coordinate is local edge
                // (worst + 1) % 3.
                t = M.triangleAdjacency[t][(worst + 1) % 3];
            }
            for (int u = 0; u < nT; ++u) {
                const auto l = bary(u, x);
                if (l[0] >= -1e-12 && l[1] >= -1e-12 && l[2] >= -1e-12) return u;
            }
            return -1;
        };
        std::vector<char> holds(nT, 0);
        for (const Singular &q : sing_) holds[q.triangle] = 1;
        for (Singular &q : sing_) {
            if (q.d <= 0) continue;
            int corner = -1;
            double best = 1e300;
            for (const auto &loop : loops)
                for (const int b : loop) {
                    if (quarters[b] != 1) continue;
                    const double local = 0.5 * (normP(M.vertices[nextOf[b]] - M.vertices[b]) +
                                                normP(M.vertices[prevOf[b]] - M.vertices[b]));
                    const double d = normP(q.pos - M.vertices[b]);
                    const Triangle &tr = M.triangles[q.triangle];
                    const bool inFan = tr[0] == b || tr[1] == b || tr[2] == b;
                    if ((inFan || d < opts_.relocationRadius * local) && d < best) { best = d; corner = b; }
                }
            if (corner < 0) continue;
            const Point &c = M.vertices[corner];
            const Point e1 = normalizeP(M.vertices[nextOf[corner]] - c);
            const Point e2 = normalizeP(M.vertices[prevOf[corner]] - c);
            const Point bis = normalizeP(e1 + e2);
            if (normP(bis) < 0.5) continue;
            const double L = std::min(runTo(corner, true), runTo(corner, false));
            const double halfAngle = 0.5 * interior[corner];
            double D = opts_.cornerRelocation * (2.0 / 3.0) * L * std::cos(halfAngle);
            holds[q.triangle] = 0;
            int t = -1;
            Point x{0.0, 0.0};
            for (int tries = 0; tries < 8 && D > 0.25 * best; ++tries, D *= 0.5) {
                x = c + bis * D;
                t = locate(x, q.triangle);
                // Never into another region (one across an interface).
                if (t >= 0 && compOf[t] != compOf[q.triangle]) t = -1;
                if (t >= 0 && !holds[t]) {
                    const auto l = bary(t, x);
                    if (std::min({l[0], l[1], l[2]}) > 0.05) break;
                    // Too near an edge of its triangle: its centroid instead.
                    const Triangle &tr = M.triangles[t];
                    x = (M.vertices[tr[0]] + M.vertices[tr[1]] + M.vertices[tr[2]]) / 3.0;
                    break;
                }
                t = -1;
            }
            if (t < 0) { holds[q.triangle] = 1; continue; }
            q.pos = x;
            q.triangle = t;
            holds[t] = 1;
            ++report_.relocated;
        }
    }
    report_.singularities = static_cast<int>(sing_.size());

    // ---- Poincare-Hopf, per component (Sec. 3.3) -----------------------------
    auto poincareHopf = [&](std::string *why) {
        std::vector<int> sum4(nComp, 0);
        for (const Singular &s : sing_) sum4[compOf[s.triangle]] += s.d;
        for (size_t l = 0; l < loops.size(); ++l)
            for (const int b : loops[l]) sum4[compOfVertex[b]] += 2 - quarters[b];
        // The shortfall over all components, for Report::indexDeficit4; set
        // by the checks that report (why != null), so it is the last such
        // check's -- Table 1's corners, or the re-read ones that made it hold.
        int deficit = 0;
        for (int c = 0; c < nComp; ++c) deficit += 4 * (2 - static_cast<int>(loopsOfComp[c].size())) - sum4[c];
        if (why) {
            report_.indexDeficit4 = 0;
            report_.regionDeficit4.clear();
        }
        for (int c = 0; c < nComp; ++c) {
            const int chi = 2 - static_cast<int>(loopsOfComp[c].size());
            if (sum4[c] != 4 * chi) {
                if (why) {
                    report_.indexDeficit4 = deficit;
                    report_.regionDeficit4.assign(nComp, 0);
                    for (int cc = 0; cc < nComp; ++cc)
                        report_.regionDeficit4[cc] = 4 * (2 - static_cast<int>(loopsOfComp[cc].size())) - sum4[cc];
                    std::ostringstream os;
                    os << "Poincare-Hopf fails: indices add to " << sum4[c] << "/4 against chi = " << chi;
                    *why = os.str();
                }
                return false;
            }
        }
        return true;
    };
    {
        std::string why;
        if (!poincareHopf(&why) && !weak && first_ && static_cast<int>(first_->u_k.size()) == nVo) {
            // Re-read the corners off the first field (see the constructor):
            // the cross turns by the principal matchings along the boundary
            // edges, and across a stretch of boundary the turn Delta and the
            // tangent's turn tau fix the stretch's index, 2 pi k = tau - Delta
            // (Sec. 2.4). Exact over any stretch that begins and ends at
            // vertices carrying the tangent, which every vertex of a straight
            // or smooth run does -- but *not* per corner when two corners are
            // neighbours on the loop (a chamfer: two 225-degree corners a
            // short edge apart, geom033), since each would claim the edge
            // between them. So the stretch is a maximal run of consecutive
            // corners, its index is read whole, and Table 1's split of it is
            // kept and adjusted at the corner nearest its Table 1 threshold
            // until it adds up.
            std::vector<int> reread = quarters;
            int changed = 0;
            auto th = [&](int v) { return 0.25 * std::arg(first_->u_k[orig(v)]); };
            // A node's sector is prescribed, and the first field has one value
            // at the node for all its sectors: no run is read across one.
            auto prescribed = [&](int v) { return !fixedQ.empty() && fixedQ[v] >= 0; };
            auto isCorner = [&](int v) {
                if (quarters[v] != 2) return true;
                const double delta = wrapQuarter(th(v) - th(prevOf[v])) + wrapQuarter(th(nextOf[v]) - th(v));
                return std::lround((M_PI - interior[v] - delta) / M_PI_2) != 0;
            };
            // The nearest Table 1 threshold to an interior angle, in the
            // direction that raises (up) or lowers q.
            auto toThreshold = [&](double a, int q, bool up) {
                const double t = up ? (q == 1 ? 0.75 * M_PI : q == 2 ? 1.25 * M_PI : q == 3 ? 1.75 * M_PI : 1e9)
                                    : (q == 2 ? 0.75 * M_PI : q == 3 ? 1.25 * M_PI : q >= 4 ? 1.75 * M_PI : -1e9);
                return std::fabs(a - t);
            };
            for (const auto &loop : loops) {
                const int m = static_cast<int>(loop.size());
                // Start the walk at a vertex that is no corner, so no run
                // straddles the seam; a loop of nothing but corners is read
                // vertex by vertex.
                int start = -1;
                for (int i = 0; i < m && start < 0; ++i) if (!isCorner(loop[i])) start = i;
                if (start < 0) continue;
                for (int i = 0; i < m;) {
                    const int v0 = loop[(start + i) % m];
                    if (!isCorner(v0)) { ++i; continue; }
                    std::vector<int> run;
                    while (i < m && isCorner(loop[(start + i) % m])) { run.push_back(loop[(start + i) % m]); ++i; }
                    const int before = prevOf[run.front()], after = nextOf[run.back()];
                    bool touches = prescribed(before) || prescribed(after);
                    for (const int v : run) touches = touches || prescribed(v);
                    if (touches) continue;
                    double delta = wrapQuarter(th(run.front()) - th(before));
                    double tau = 0.0;
                    for (size_t r = 0; r < run.size(); ++r) {
                        tau += M_PI - interior[run[r]];
                        const int nx = (r + 1 < run.size()) ? run[r + 1] : after;
                        delta += wrapQuarter(th(nx) - th(run[r]));
                    }
                    const int want = static_cast<int>(std::lround((tau - delta) / M_PI_2));   // sum (2 - q)
                    int have = 0;
                    for (const int v : run) have += 2 - reread[v];
                    while (have != want) {
                        const bool raise = have > want;   // more q, less (2 - q)
                        int best = -1;
                        double bestD = 1e300;
                        for (const int v : run) {
                            if (!raise && reread[v] <= 1) continue;
                            const double d = toThreshold(interior[v], reread[v], raise);
                            if (d < bestD) { bestD = d; best = v; }
                        }
                        if (best < 0) break;
                        reread[best] += raise ? 1 : -1;
                        have += raise ? -1 : 1;
                        ++changed;
                    }
                }
            }
            if (changed > 0) {
                std::swap(quarters, reread);
                if (poincareHopf(nullptr)) {
                    report_.cornersReread = changed;
                    report_.corners = 0;
                    for (const auto &loop : loops)
                        for (const int b : loop) if (quarters[b] != 2) ++report_.corners;
                } else {
                    std::swap(quarters, reread);
                }
            }
        }
        if (!poincareHopf(&why)) return fail(why);
        report_.poincareHopf = true;
    }

    // ---- weak boundary conditions (Options::boundarySnap) ---------------------
    //
    // After the Poincare-Hopf check, so that the corners it may have re-read
    // are settled and a snapped vertex is never re-read: the edit moves a
    // quarter from the pattern to a vertex of the same component and leaves
    // every component's sum as it was.
    if (opts_.boundarySnap > 0.0 && !sing_.empty()) {
        std::vector<char> nodeVertex(nVo, 0);
        if (interfaces_)
            for (const auto &n : interfaces_->nodes())
                if (n.vertex >= 0 && n.vertex < nVo) nodeVertex[n.vertex] = 1;
        // dS proper: not an interface (a cut vertex whose original is no
        // boundary vertex), and not a node, whose sectors are prescribed.
        auto onDS = [&](int v) {
            const int o = orig(v);
            return o >= 0 && o < nVo && mesh_->isBoundaryVertex[o] && !nodeVertex[o] &&
                   (fixedQ.empty() || fixedQ[v] < 0);
        };
        std::vector<char> gone(sing_.size(), 0);
        for (size_t i = 0; i < sing_.size(); ++i) {
            const Singular &s = sing_[i];
            if (s.d >= 0) continue;
            const int comp = compOf[s.triangle];
            double best = 1e300;
            int bv = -1;
            for (const auto &loop : loops) {
                if (compOfVertex[loop[0]] != comp) continue;
                for (const int a : loop) {
                    const int b = nextOf[a];
                    const Point &pa = M.vertices[a], &pb = M.vertices[b];
                    const Point d = pb - pa;
                    const double dd = dotP(d, d);
                    const double t = dd > 0.0 ? std::min(1.0, std::max(0.0, dotP(s.pos - pa, d) / dd)) : 0.0;
                    const double dist = normP(s.pos - (pa + d * t));
                    if (dist < best) { best = dist; bv = (t < 0.5) ? a : b; }
                }
            }
            if (bv < 0 || !onDS(bv)) continue;
            const double local = 0.5 * (normP(M.vertices[nextOf[bv]] - M.vertices[bv]) +
                                        normP(M.vertices[prevOf[bv]] - M.vertices[bv]));
            if (best > opts_.boundarySnap * local) continue;
            const int q = quarters[bv] - s.d;
            if (q <= 0 || interior[bv] / q < opts_.snapMinAngle) continue;
            quarters[bv] = q;
            gone[i] = 1;
            bsing_.push_back(BoundarySingular{M.vertices[bv], bv, s.d});
            ++report_.snapped;
        }
        if (report_.snapped > 0) {
            std::vector<Singular> kept;
            for (size_t i = 0; i < sing_.size(); ++i)
                if (!gone[i]) kept.push_back(sing_[i]);
            sing_.swap(kept);
            report_.singularities = static_cast<int>(sing_.size());
            report_.corners = 0;
            for (const auto &loop : loops)
                for (const int b : loop) if (quarters[b] != 2) ++report_.corners;
            if (!poincareHopf(nullptr)) return fail("Poincare-Hopf fails after the boundary snap");
        }
    }

    // ---- the right-hand side of h (Sec. 6.3, with the singular part out) -----
    //
    // The change of the singular part along a straight segment that misses
    // every X_j, exactly: each term turns by the angle the segment subtends.
    auto singularChange = [&](const Point &p, const Point &q) {
        double s = 0.0;
        for (const Singular &j : sing_) {
            const double a0 = std::atan2(p[1] - j.pos[1], p[0] - j.pos[0]);
            const double a1 = std::atan2(q[1] - j.pos[1], q[0] - j.pos[0]);
            s += 0.25 * j.d * wrap_pi(a1 - a0);
        }
        return s;
    };
    std::vector<double> load(nV, 0.0);
    double loadScale = 0.0;
    for (const auto &loop : loops) {
        for (const int b : loop) {
            const int a = prevOf[b], c = nextOf[b];
            const Point &xb = M.vertices[b];
            const Point mIn = (M.vertices[a] + xb) * 0.5;
            const Point mOut = (xb + M.vertices[c]) * 0.5;
            const double tau = M_PI - interior[b];
            const double cornerTurn = M_PI_2 * (2 - quarters[b]);   // 2 pi k_b
            // The points on dS: from mIn to mOut a term turns by the angle
            // the two half-edges subtend at it, and at its own vertex, where
            // both half-edges end on it, by the sweep through the model from
            // the one to the other, -alpha (the model is on the left).
            double onBoundary = 0.0;
            for (const BoundarySingular &j : bsing_) {
                if (j.region >= 0 && j.region != compOfVertex[b]) continue;
                if (j.vertex == b) { onBoundary += 0.5 * j.d * -interior[b]; continue; }
                const double a0 = std::atan2(mIn[1] - j.pos[1], mIn[0] - j.pos[0]);
                const double a1 = std::atan2(xb[1] - j.pos[1], xb[0] - j.pos[0]);
                const double a2 = std::atan2(mOut[1] - j.pos[1], mOut[0] - j.pos[0]);
                onBoundary += 0.5 * j.d * (wrap_pi(a1 - a0) + wrap_pi(a2 - a1));
            }
            load[b] = tau - cornerTurn - (singularChange(mIn, xb) + singularChange(xb, mOut)) - onBoundary;
            loadScale = std::max(loadScale, std::fabs(load[b]));
        }
    }
    {
        std::vector<double> sum(nComp, 0.0);
        for (int v = 0; v < nV; ++v)
            if (loopOf[v] >= 0) sum[compOfVertex[v]] += load[v];
        for (int c = 0; c < nComp; ++c)
            report_.loadSum = std::max(report_.loadSum, std::fabs(sum[c]) / std::max(1.0, loadScale));
    }

    // ---- the stiffness, and h ------------------------------------------------
    std::vector<Point> gradH(nT, Point{0.0, 0.0});
    std::vector<std::array<Point, 3>> grads(nT);
    std::vector<double> area(nT, 0.0);
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(static_cast<size_t>(nT) * 9);
    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = M.triangles[t];
        const Point &p0 = M.vertices[tri[0]], &p1 = M.vertices[tri[1]], &p2 = M.vertices[tri[2]];
        const double a2 = cross2(p1 - p0, p2 - p0);
        area[t] = 0.5 * std::fabs(a2);
        if (std::fabs(a2) < 1e-300) continue;
        grads[t][0] = Point{(p1[1] - p2[1]) / a2, (p2[0] - p1[0]) / a2};
        grads[t][1] = Point{(p2[1] - p0[1]) / a2, (p0[0] - p2[0]) / a2};
        grads[t][2] = Point{(p0[1] - p1[1]) / a2, (p1[0] - p0[0]) / a2};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                trips.emplace_back(tri[i], tri[j], area[t] * dotP(grads[t][i], grads[t][j]));
    }
    // One vertex per component pinned to 0 (h is defined up to a constant),
    // and every vertex no triangle uses.
    std::vector<char> pinned(nV, 0);
    for (int v = 0; v < nV; ++v) if (!used[v]) pinned[v] = 1;
    for (int c = 0; c < nComp; ++c) pinned[loops[outerOf[c]][0]] = 1;
    std::vector<int> freeIdx(nV, -1);
    int nFree = 0;
    for (int v = 0; v < nV; ++v) if (!pinned[v]) freeIdx[v] = nFree++;
    {
        std::vector<Eigen::Triplet<double>> red;
        red.reserve(trips.size());
        for (const auto &tr : trips)
            if (freeIdx[tr.row()] >= 0 && freeIdx[tr.col()] >= 0)
                red.emplace_back(freeIdx[tr.row()], freeIdx[tr.col()], tr.value());
        Eigen::SparseMatrix<double> K(nFree, nFree);
        K.setFromTriplets(red.begin(), red.end());
        Eigen::VectorXd rhs(nFree);
        for (int v = 0; v < nV; ++v) if (freeIdx[v] >= 0) rhs[freeIdx[v]] = load[v];
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(K);
        if (ldlt.info() != Eigen::Success) return fail("the Neumann problem for h would not factorise");
        const Eigen::VectorXd sol = ldlt.solve(rhs);
        if (ldlt.info() != Eigen::Success) return fail("the Neumann problem for h would not solve");
        h_.assign(nV, 0.0);
        for (int v = 0; v < nV; ++v) if (freeIdx[v] >= 0) h_[v] = sol[freeIdx[v]];
    }
    std::vector<Point> w(nT, Point{0.0, 0.0});   // J grad h: the gradient of theta_h
    for (int t = 0; t < nT; ++t) {
        Point g{0.0, 0.0};
        for (int i = 0; i < 3; ++i) g = g + grads[t][i] * h_[M.triangles[t][i]];
        gradH[t] = g;
        w[t] = Point{-g[1], g[0]};
        report_.maxRegularTurnPerEdge = std::max(report_.maxRegularTurnPerEdge, normP(w[t]) * meanEdge);
    }

    // ---- theta_h, by Crouzeix-Raviart integration (Sec. 6.5) ----------------
    std::vector<std::array<double, 3>> cr(nT, std::array<double, 3>{{0.0, 0.0, 0.0}});
    std::vector<char> seen(nT, 0);
    auto fill = [&](int t, int k, double value) {
        const Point m = midpoint(M, t, k);
        for (int kk = 0; kk < 3; ++kk) cr[t][kk] = value + dotP(w[t], midpoint(M, t, kk) - m);
    };
    auto localEdgeOf = [&](int t, int e) {
        for (int k = 0; k < 3; ++k) if (M.triangleEdges[t][k] == e) return k;
        return -1;
    };
    for (int c = 0; c < nComp; ++c) {
        const auto &outer = loops[outerOf[c]];
        const int v0 = outer[0], v1 = nextOf[v0];
        const int e0 = edgeAfter[v0];
        const int t0 = (M.edgeTriangles[e0][0] >= 0) ? M.edgeTriangles[e0][0] : M.edgeTriangles[e0][1];
        const int k0 = localEdgeOf(t0, e0);
        if (k0 < 0) return fail("mesh edge tables disagree");
        const Point m0 = (M.vertices[v0] + M.vertices[v1]) * 0.5;
        // The anchor: the cross at the midpoint of the first edge of the outer
        // loop has a branch along it (P2 eq. 15's theta(P) = theta_P).
        const double phi0 = computeAngle(M.vertices[v1] - M.vertices[v0]);
        fill(t0, k0, phi0 - singularPart(m0, c));
        seen[t0] = 1;
        std::queue<int> q;
        q.push(t0);
        while (!q.empty()) {
            const int t = q.front();
            q.pop();
            for (int k = 0; k < 3; ++k) {
                const int n = M.triangleAdjacency[t][k];
                if (n < 0 || seen[n]) continue;
                const int kn = localEdgeOf(n, M.triangleEdges[t][k]);
                if (kn < 0) return fail("mesh edge tables disagree");
                fill(n, kn, cr[t][k]);
                seen[n] = 1;
                q.push(n);
            }
        }
    }
    for (int i = 0; i < static_cast<int>(M.edges.size()); ++i) {
        const int t = M.edgeTriangles[i][0], u = M.edgeTriangles[i][1];
        if (t < 0 || u < 0) continue;
        const int kt = localEdgeOf(t, i), ku = localEdgeOf(u, i);
        if (kt < 0 || ku < 0) continue;
        report_.maxMidpointMismatch =
            std::max(report_.maxMidpointMismatch, std::fabs(wrapQuarter(cr[t][kt] - cr[u][ku])));
    }

    // The tangency of the CR field at the midpoints of each loop's edges: the
    // offset delta of Sec. 6.8, which on the anchor's loop is zero by
    // construction.
    auto loopOffsets = [&](const std::vector<int> &loop, double &mean, double &spread) {
        std::complex<double> acc(0.0, 0.0);
        std::vector<double> d;
        for (const int a : loop) {
            const int e = edgeAfter[a], b = nextOf[a];
            const int t = (M.edgeTriangles[e][0] >= 0) ? M.edgeTriangles[e][0] : M.edgeTriangles[e][1];
            const int k = localEdgeOf(t, e);
            const Point m = (M.vertices[a] + M.vertices[b]) * 0.5;
            const double phi = computeAngle(M.vertices[b] - M.vertices[a]);
            const double off = wrapQuarter(cr[t][k] + singularPart(m, compOf[t]) - phi);
            d.push_back(off);
            acc += std::polar(1.0, 4.0 * off);
        }
        mean = 0.25 * std::arg(acc);
        spread = 0.0;
        for (const double x : d) spread = std::max(spread, std::fabs(wrapQuarter(x - mean)));
    };

    // ---- the regular part per vertex -----------------------------------------
    //
    // Each triangle's affine theta_h read at the vertex, lifted round the fan
    // one triangle at a time (a sharp corner's fan turns by up to pi/4 in all,
    // more than a single lift against the first value could take), averaged.
    regular_.assign(nV, 0.0);
    {
        std::vector<std::vector<int>> trisOf(nV);
        for (int t = 0; t < nT; ++t)
            for (int k = 0; k < 3; ++k) trisOf[M.triangles[t][k]].push_back(t);
        auto valueAt = [&](int t, int v) {
            for (int k = 0; k < 3; ++k)
                if (M.triangles[t][k] == v) return cr[t][k] + dotP(w[t], M.vertices[v] - midpoint(M, t, k));
            return 0.0;
        };
        for (int v = 0; v < nV; ++v) {
            const auto &ts = trisOf[v];
            if (ts.empty()) continue;
            // Order the fan by walking across the edges at v: start at a
            // triangle with a boundary edge at v when there is one.
            auto edgesAtV = [&](int t, int out[2]) {
                int n = 0;
                for (int k = 0; k < 3; ++k) {
                    const int a = M.triangles[t][k], b = M.triangles[t][(k + 1) % 3];
                    if (a == v || b == v) out[n++] = k;
                }
            };
            int start = ts[0], startEdge = -1;
            for (const int t : ts) {
                int ek[2];
                edgesAtV(t, ek);
                for (const int k : ek)
                    if (M.triangleAdjacency[t][k] < 0) { start = t; startEdge = k; }
                if (startEdge >= 0) break;
            }
            std::vector<int> fan;
            {
                int t = start, came = startEdge;
                if (came < 0) { int ek[2]; edgesAtV(t, ek); came = ek[0]; }
                for (size_t guard = 0; guard <= ts.size() && t >= 0; ++guard) {
                    fan.push_back(t);
                    int ek[2];
                    edgesAtV(t, ek);
                    const int leave = (ek[0] == came) ? ek[1] : ek[0];
                    const int n = M.triangleAdjacency[t][leave];
                    if (n < 0 || n == start) break;
                    came = localEdgeOf(n, M.triangleEdges[t][leave]);
                    t = n;
                }
                if (fan.size() != ts.size()) fan = ts;   // non-manifold: take them as they come
            }
            double prev = valueAt(fan[0], v), sum = prev, lo = prev, hi = prev;
            for (size_t i = 1; i < fan.size(); ++i) {
                const double x = liftNear(prev, valueAt(fan[i], v));
                sum += x;
                lo = std::min(lo, x);
                hi = std::max(hi, x);
                prev = x;
            }
            regular_[v] = sum / static_cast<double>(fan.size());
            if (loopOf[v] < 0) report_.vertexSpread = std::max(report_.vertexSpread, 0.5 * (hi - lo));
        }
    }

    // ---- holes (Sec. 6.8) ----------------------------------------------------
    std::vector<double> offsetOfLoop(loops.size(), 0.0);
    for (size_t l = 0; l < loops.size(); ++l) {
        double mean = 0.0, spread = 0.0;
        loopOffsets(loops[l], mean, spread);
        offsetOfLoop[l] = mean;
        if (loopTurning[l] > 0.0) {
            report_.anchorTangency = std::max(report_.anchorTangency, std::fabs(mean) + spread);
        } else {
            report_.holeOffsets.push_back(mean);
            report_.holeSpread = std::max(report_.holeSpread, spread);
        }
    }
    if (opts_.alignHoles && report_.holes > 0) {
        // Harmonic measure of each hole: Dirichlet 1 on it and 0 on every
        // other loop of its component, harmonic in between. One factorisation
        // serves them all.
        std::vector<int> idx(nV, -1);
        int n = 0;
        for (int v = 0; v < nV; ++v) if (used[v] && loopOf[v] < 0) idx[v] = n++;
        std::vector<Eigen::Triplet<double>> red;
        red.reserve(trips.size());
        for (const auto &tr : trips)
            if (idx[tr.row()] >= 0 && idx[tr.col()] >= 0) red.emplace_back(idx[tr.row()], idx[tr.col()], tr.value());
        Eigen::SparseMatrix<double> Kii(n, n);
        Kii.setFromTriplets(red.begin(), red.end());
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
        if (n > 0) {
            ldlt.compute(Kii);
            if (ldlt.info() != Eigen::Success) return fail("the harmonic measures would not factorise");
        }
        for (size_t l = 0; l < loops.size(); ++l) {
            if (loopTurning[l] > 0.0 || offsetOfLoop[l] == 0.0) continue;
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);
            for (const auto &tr : trips)
                if (idx[tr.row()] >= 0 && idx[tr.col()] < 0 && loopOf[tr.col()] == static_cast<int>(l))
                    rhs[idx[tr.row()]] -= tr.value();
            Eigen::VectorXd omega = (n > 0) ? Eigen::VectorXd(ldlt.solve(rhs)) : Eigen::VectorXd();
            const int comp = compOfVertex[loops[l][0]];
            for (int v = 0; v < nV; ++v) {
                if (compOfVertex[v] != comp) continue;
                const double om = (idx[v] >= 0) ? omega[idx[v]] : (loopOf[v] == static_cast<int>(l) ? 1.0 : 0.0);
                regular_[v] -= offsetOfLoop[l] * om;
            }
        }
    }
    // Tangency of the per-vertex field the tracer reads, at every boundary
    // edge midpoint, on the holes.
    for (size_t l = 0; l < loops.size(); ++l) {
        if (loopTurning[l] > 0.0) continue;
        for (const int a : loops[l]) {
            const int b = nextOf[a];
            const Point m = (M.vertices[a] + M.vertices[b]) * 0.5;
            const double ra = regular_[a], rb = liftNear(ra, regular_[b]);
            const double th = 0.5 * (ra + rb) + singularPart(m, compOfVertex[a]);
            report_.holeTangency = std::max(
                report_.holeTangency, std::fabs(wrapQuarter(th - computeAngle(M.vertices[b] - M.vertices[a]))));
        }
    }

    // ---- with interfaces: back onto the uncut mesh ----------------------------
    //
    // A vertex of an interface has one copy per region at it; the copies agree
    // mod pi/2 up to the discretisation (both regions are tangent to the
    // interface), and are joined by their mean mod pi/2, lifted to the first.
    // A node's sectors are corners and may honestly differ; only the other
    // interface vertices are measured.
    if (cutMesh) {
        // With a term of one region (BoundarySingular::region) the copies'
        // regular parts no longer agree mod pi/2 at an interface vertex --
        // their fields do -- so they are kept per triangle corner, and what
        // is joined is the field itself, the regular part handed on per
        // vertex being that less singularPart(x).
        auto copyValue = [&](int c) {
            return regionTerms_ ? regular_[c] + singularPart(M.vertices[c], compOfVertex[c]) : regular_[c];
        };
        if (regionTerms_) {
            const Mesh &O = *mesh_;
            regularCorner_.assign(O.triangles.size(), std::array<double, 3>{{0.0, 0.0, 0.0}});
            for (int t = 0; t < static_cast<int>(O.triangles.size()); ++t)
                for (int k = 0; k < 3; ++k)
                    for (int kk = 0; kk < 3; ++kk)
                        if (origOf[M.triangles[t][kk]] == O.triangles[t][k]) {
                            regularCorner_[t][k] = regular_[M.triangles[t][kk]];
                            break;
                        }
        }
        std::vector<std::complex<double>> acc(nVo, std::complex<double>(0.0, 0.0));
        std::vector<int> copies(nVo, 0), firstCopy(nVo, -1);
        for (int c = 0; c < nV; ++c) {
            const int o = origOf[c];
            acc[o] += std::polar(1.0, 4.0 * copyValue(c));
            if (copies[o]++ == 0) firstCopy[o] = c;
        }
        std::vector<char> isNode(nVo, 0);
        for (const auto &n : interfaces_->nodes())
            if (n.vertex >= 0 && n.vertex < nVo) isNode[n.vertex] = 1;
        std::vector<double> regular(nVo, 0.0), h(nVo, 0.0);
        std::vector<int> q(nVo, 2);
        for (int o = 0; o < nVo; ++o) {
            if (firstCopy[o] < 0) continue;
            const double ref = copyValue(firstCopy[o]);
            regular[o] = (copies[o] == 1) ? ref : liftNear(ref, 0.25 * std::arg(acc[o]));
            if (regionTerms_) regular[o] -= singularPart(mesh_->vertices[o]);
            h[o] = h_[firstCopy[o]];
            q[o] = quarters[firstCopy[o]];
        }
        for (int c = 0; c < nV; ++c) {
            const int o = origOf[c];
            if (copies[o] > 1 && !isNode[o]) {
                const double mine = regionTerms_ ? copyValue(c) - singularPart(mesh_->vertices[o]) : regular_[c];
                report_.interfaceMismatch =
                    std::max(report_.interfaceMismatch, std::fabs(wrapQuarter(mine - regular[o])));
            }
        }
        // dS corners as the uncut fan turns: Table 1's reading of the whole
        // corner at a landing (SeparatrixTrace reads it as the node's fan).
        for (const auto &n : interfaces_->nodes())
            if (n.vertex >= 0 && n.vertex < nVo) q[n.vertex] = n.onBoundary ? landingFan[n.vertex] : 2;
        regular_.swap(regular);
        h_.swap(h);
        quarters.swap(q);
        report_.regions = nComp;
    }
    quarters_ = quarters;

    // ---- the field as a CrossField -------------------------------------------
    cross_ = std::make_shared<CrossField>(mesh_);
    cross_->u_k.resize(nVo);
    for (int v = 0; v < nVo; ++v)
        cross_->u_k[v] = std::polar(1.0, 4.0 * (regular_[v] + singularPart(mesh_->vertices[v])));
    cross_->u_k_prev = cross_->u_k;
    cross_->error = 0.0;

    report_.built = true;
}
