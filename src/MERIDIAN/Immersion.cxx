#include "Immersion.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace {

// The edge of a triangle joining two of its local corners. Local edge p joins
// tri[p] and tri[(p+1)%3], the same convention Mesh and RicciFlow use.
inline int localEdgeBetween(int m, int n) { return ((m + 1) % 3 == n) ? m : n; }

inline double safeAcos(double c) {
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    return std::acos(c);
}

// Angle at corner a of a triangle given the three image points.
inline double cornerAngle(const Point &a, const Point &b, const Point &c) {
    const Point u = b - a;
    const Point v = c - a;
    const double nu = normP(u), nv = normP(v);
    if (nu <= 0.0 || nv <= 0.0) return 0.0;
    return safeAcos(dotP(u, v) / (nu * nv));
}

inline double signedArea(const Point &a, const Point &b, const Point &c) {
    return 0.5 * cross2(b - a, c - a);
}

// Area from three side lengths, by the numerically stable form of Heron's
// formula (Kahan): sort descending, then the product below. The textbook form
// loses everything on a sliver, and the flat metric produces slivers.
double heronArea(double a, double b, double c) {
    if (a < b) std::swap(a, b);
    if (a < c) std::swap(a, c);
    if (b < c) std::swap(b, c);
    const double t = (a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c));
    return (t > 0.0) ? 0.25 * std::sqrt(t) : 0.0;
}

} // namespace

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
Immersion::Immersion(const ConeCut &cut, const RicciFlow &ricci, const ConeSingularities &conesIn)
    : cutter(&cut), cones(&conesIn), orig(cut.getOriginalMeshPtr()),
      cutMesh(cut.getCutMesh()) {
    if (!orig) throw std::runtime_error("Immersion: the cut has no source mesh");
    if (&ricci.getMesh() != orig.get()) {
        throw std::runtime_error("Immersion: the flow ran on a different mesh from the cut");
    }
    if (&conesIn.getMesh() != orig.get()) {
        throw std::runtime_error("Immersion: the cones were measured on a different mesh");
    }
    if (cutMesh.triangles.size() != orig->triangles.size()) {
        throw std::runtime_error("Immersion: Omega does not have one face per face of S");
    }

    // The cones, and where each of them went in Omega.
    const auto &o2c = cutter->getOriginalToCutVertices();
    cutVertCone.assign(cutMesh.vertices.size(), -1);
    for (const auto &c : cones->getCones()) {
        const int slot = static_cast<int>(coneVertex.size());
        coneVertex.push_back(c.vertex);
        coneIndex.push_back(c.index);
        coneChildren.push_back(o2c[c.vertex]);
        for (int cv : o2c[c.vertex]) cutVertCone[cv] = slot;
    }

    buildLengths(ricci);
    layout();
    normaliseOrientation();
    buildArcs();
    fitTransitions();
    alignToAxes();
    check();
}

// ---------------------------------------------------------------------------
// The field constructor -- docs/cf_flow_pipeline.md Sec. 6.2
//
// buildLengths() and layout() are the only two things the Ricci route does that
// a field route has already done for itself, so they are the only two skipped.
// The cone bookkeeping above them and the four stages below them are common to
// both, and they are common *by construction*: every one of them reads ConeCut
// and nothing else about the metric.
//
// normaliseOrientation() is skipped as well, and deliberately. It exists
// because the unfolding's seed triangle fixes the handedness of the whole
// immersion to whatever that one face's winding happened to be, which is an
// artefact of the unfolding. An integrated map has no seed: its handedness is
// the field's, det J*_t > 0 everywhere, and a negative total area here would be
// a real fold and not a convention to be reflected away. Reflecting it would
// hide exactly the failure Q1 is here to catch.
// ---------------------------------------------------------------------------
Immersion::Immersion(const ConeCut &cut, const ConeSingularities &conesIn,
                     const std::vector<Point> &psi,
                     const std::vector<double> &flatEdgeLengths,
                     const std::vector<double> &faceAngleIn)
    : cutter(&cut), cones(&conesIn), orig(cut.getOriginalMeshPtr()),
      cutMesh(cut.getCutMesh()) {
    if (!orig) throw std::runtime_error("Immersion: the cut has no source mesh");
    if (&conesIn.getMesh() != orig.get()) {
        throw std::runtime_error("Immersion: the cones were measured on a different mesh");
    }
    if (cutMesh.triangles.size() != orig->triangles.size()) {
        throw std::runtime_error("Immersion: Omega does not have one face per face of S");
    }
    if (psi.size() != cutMesh.vertices.size()) {
        throw std::runtime_error("Immersion: the integrated map has one point per vertex of "
                                 "Omega, and this one does not");
    }
    if (flatEdgeLengths.size() != orig->edges.size()) {
        throw std::runtime_error("Immersion: the field metric has one length per edge of S, "
                                 "and this one does not");
    }
    if (!faceAngleIn.empty() && faceAngleIn.size() != orig->triangles.size()) {
        throw std::runtime_error("Immersion: the combed field angle has one value per face of "
                                 "S, and this one does not");
    }

    fromField = true;
    faceAngle = faceAngleIn;

    const auto &o2c = cutter->getOriginalToCutVertices();
    cutVertCone.assign(cutMesh.vertices.size(), -1);
    for (const auto &c : cones->getCones()) {
        const int slot = static_cast<int>(coneVertex.size());
        coneVertex.push_back(c.vertex);
        coneIndex.push_back(c.index);
        coneChildren.push_back(o2c[c.vertex]);
        for (int cv : o2c[c.vertex]) cutVertCone[cv] = slot;
    }

    flatLen = flatEdgeLengths;
    buildCutLengths();
    uv = psi;
    report.cutVertices = static_cast<int>(cutMesh.vertices.size());
    report.placedVertices = report.cutVertices;
    report.unplacedVertices = 0;

    buildArcs();
    fitTransitions();
    alignToAxes();
    check();
}

// ---------------------------------------------------------------------------
// buildLengths()
//
// The flat metric arrives as one length per edge of the *input* mesh; Omega
// needs it per edge of Omega, and the layout needs it per (face, local edge).
// All three are the same numbers under different indexings, because a face of
// Omega is a face of S with its corners renamed -- ConeCut keeps the face list
// in order, which is what makes the local edge of a cut face and the local edge
// of its parent the same edge.
// ---------------------------------------------------------------------------
void Immersion::buildLengths(const RicciFlow &ricci) {
    int recovered = 0;
    flatLen = ricci.originalEdgeLengthsCompleted(&recovered);
    report.recoveredEdgeLengths = recovered;
    if (recovered > 0) {
        std::ostringstream oss;
        oss << recovered << " input edge(s) the flow flipped away were given their Eq. (11) "
            << "length from the inversive distance they started with.";
        report.messages.push_back(oss.str());
    }

    buildCutLengths();
}

// ---------------------------------------------------------------------------
// buildCutLengths()
//
// flatLen re-indexed by the edges of Omega. Split out of buildLengths() because
// the field constructor is handed flatLen ready-made -- from the field metric
// rather than from the flow -- and still needs this half of it.
// ---------------------------------------------------------------------------
void Immersion::buildCutLengths() {
    std::unordered_map<EdgeKey, int, EdgeKeyHash> origEdgeIndex;
    origEdgeIndex.reserve(orig->edges.size() * 2);
    for (int e = 0; e < static_cast<int>(orig->edges.size()); ++e) {
        origEdgeIndex.emplace(EdgeKey(orig->edges[e][0], orig->edges[e][1]), e);
    }

    const auto &c2o = cutter->getCutVertexToOriginal();
    cutLen.assign(cutMesh.edges.size(), 0.0);
    for (int e = 0; e < static_cast<int>(cutMesh.edges.size()); ++e) {
        const int pa = c2o[cutMesh.edges[e][0]];
        const int pb = c2o[cutMesh.edges[e][1]];
        auto it = origEdgeIndex.find(EdgeKey(pa, pb));
        cutLen[e] = (it == origEdgeIndex.end()) ? 0.0 : flatLen[it->second];
    }
}

// ---------------------------------------------------------------------------
// placeThird()
// ---------------------------------------------------------------------------
bool Immersion::placeThird(const Point &pi, const Point &pj, const Point &away,
                           double lik, double ljk, Point &pk) const {
    const Point d = pj - pi;
    const double dd = normP(d);
    if (!(dd > 0.0)) { pk = pi; return false; }

    const Point e1 = d / dd;
    const Point e2{-e1[1], e1[0]};

    const double x = (lik * lik - ljk * ljk + dd * dd) / (2.0 * dd);
    double y2 = lik * lik - x * x;
    bool ok = true;
    if (!(y2 > 0.0)) { y2 = 0.0; ok = false; }   // the face does not close
    const double y = std::sqrt(y2);

    // cross2(pj - pi, p - pi) is positive on the e2 side of the line, so the
    // sign of the root is fixed by putting p_k opposite the corner the queue
    // arrived from. Doing it this way rather than from a stored per-face
    // orientation means the unfolding never has to agree with the input mesh's
    // winding, only with itself.
    const double sideAway = cross2(d, away - pi);
    const double s = (sideAway > 0.0) ? -1.0 : 1.0;

    pk = pi + e1 * x + e2 * (s * y);
    return ok;
}

// ---------------------------------------------------------------------------
// layout()
// ---------------------------------------------------------------------------
void Immersion::layout() {
    const int nCV = static_cast<int>(cutMesh.vertices.size());
    const int nF = static_cast<int>(cutMesh.triangles.size());

    // Flat lengths per (face, local edge), looked up once.
    std::vector<std::array<double, 3>> faceLen(nF);
    for (int f = 0; f < nF; ++f) {
        for (int p = 0; p < 3; ++p) {
            const int oe = orig->triangleEdges[f][p];
            faceLen[f][p] = (oe >= 0) ? flatLen[oe] : 0.0;
        }
    }

    uv.assign(nCV, Point{0.0, 0.0});
    std::vector<char> placedV(nCV, 0);
    std::vector<char> visitedF(nF, 0);

    // The seed: first edge on the u axis, third vertex above it. Everything
    // after is relative to this, and the global rotation at the end replaces
    // the choice anyway.
    {
        const Triangle &t = cutMesh.triangles[0];
        const double l01 = faceLen[0][0];
        const double l12 = faceLen[0][1];
        const double l20 = faceLen[0][2];
        uv[t[0]] = Point{0.0, 0.0};
        uv[t[1]] = Point{l01, 0.0};
        const double x = (l20 * l20 - l12 * l12 + l01 * l01) / (2.0 * std::max(l01, 1e-300));
        const double y2 = l20 * l20 - x * x;
        if (!(y2 > 0.0)) ++report.degenerateFaces;
        uv[t[2]] = Point{x, (y2 > 0.0) ? std::sqrt(y2) : 0.0};
        placedV[t[0]] = placedV[t[1]] = placedV[t[2]] = 1;
        visitedF[0] = 1;
    }

    std::queue<int> q;
    q.push(0);
    while (!q.empty()) {
        const int f = q.front();
        q.pop();
        const Triangle &tf = cutMesh.triangles[f];

        for (int p = 0; p < 3; ++p) {
            const int g = cutMesh.triangleAdjacency[f][p];
            if (g < 0 || visitedF[g]) continue;

            const int a = tf[p];
            const int b = tf[(p + 1) % 3];
            const int m = tf[(p + 2) % 3];

            const Triangle &tg = cutMesh.triangles[g];
            int qa = -1, qb = -1, qk = -1;
            for (int c = 0; c < 3; ++c) {
                if (tg[c] == a) qa = c;
                else if (tg[c] == b) qb = c;
                else qk = c;
            }
            if (qa < 0 || qb < 0 || qk < 0) continue;   // not the edge we thought

            const int k = tg[qk];
            const double lak = faceLen[g][localEdgeBetween(qa, qk)];
            const double lbk = faceLen[g][localEdgeBetween(qb, qk)];

            Point pk;
            if (!placeThird(uv[a], uv[b], uv[m], lak, lbk, pk)) ++report.degenerateFaces;

            if (placedV[k]) {
                // Reached twice. On a flat metric over a disk the two answers
                // are the same point; keep the first and record the gap.
                const double scale = std::max(lak, 1e-300);
                report.maxClosureGap = std::max(report.maxClosureGap, normP(pk - uv[k]) / scale);
            } else {
                uv[k] = pk;
                placedV[k] = 1;
            }

            visitedF[g] = 1;
            q.push(g);
        }
    }

    report.cutVertices = nCV;
    report.placedVertices = 0;
    for (int v = 0; v < nCV; ++v) if (placedV[v]) ++report.placedVertices;
    report.unplacedVertices = nCV - report.placedVertices;
    if (report.unplacedVertices > 0) {
        std::ostringstream oss;
        oss << report.unplacedVertices << " vertex/vertices of Omega were never reached; "
            << "the cut mesh is not connected.";
        report.messages.push_back(oss.str());
    }
}

// ---------------------------------------------------------------------------
// normaliseOrientation()
//
// The seed was laid down with its third vertex above its first edge, which
// fixes the orientation of the whole immersion to whatever that one face's
// winding happens to be. Stage 6 wants det J > 0 on every triangle against a
// positively oriented reference, so the map is reflected once here if the seed
// came out the wrong way round. A reflection is an isometry, so it changes no
// length and no angle -- only the sign -- and it is applied before the
// transitions are fitted so that the k it hands to Gamma_Hol_k is the k of the
// map that is actually returned.
// ---------------------------------------------------------------------------
void Immersion::normaliseOrientation() {
    double total = 0.0;
    for (const Triangle &t : cutMesh.triangles) {
        total += signedArea(uv[t[0]], uv[t[1]], uv[t[2]]);
    }
    if (total >= 0.0) return;

    for (Point &p : uv) p[0] = -p[0];
    report.mirrored = true;
}

// ---------------------------------------------------------------------------
// buildArcs()
//
// G, as an edge set on the input vertices, decomposed into maximal chains
// between its nodes. A node is a vertex with any number of G-edges other than
// two, or one that lies on dS -- so a cone (a leaf, one edge) is a node, a
// junction where a later cone path landed on an earlier one is a node, and the
// point where an arc reaches the boundary is a node.
// ---------------------------------------------------------------------------
void Immersion::buildArcs() {
    const auto &cutEdges = cutter->getCutEdges();
    if (cutEdges.empty()) return;

    std::unordered_map<int, std::vector<int>> adj;
    adj.reserve(cutEdges.size() * 2);
    for (const EdgeKey &e : cutEdges) {
        adj[e.a].push_back(e.b);
        adj[e.b].push_back(e.a);
    }

    auto isNode = [&](int v) -> bool {
        auto it = adj.find(v);
        if (it == adj.end()) return false;
        return it->second.size() != 2 || orig->isBoundaryVertex[v];
    };

    std::unordered_map<EdgeKey, char, EdgeKeyHash> used;
    used.reserve(cutEdges.size() * 2);

    auto walk = [&](int start, int first) {
        std::vector<int> path{start, first};
        used[EdgeKey(start, first)] = 1;
        int prev = start, cur = first;
        // Guard against a malformed graph rather than trusting the node test.
        for (size_t guard = 0; guard <= cutEdges.size() + 1 && !isNode(cur); ++guard) {
            const std::vector<int> &nb = adj[cur];
            int next = -1;
            for (int w : nb) if (w != prev) { next = w; break; }
            if (next < 0) break;
            used[EdgeKey(cur, next)] = 1;
            path.push_back(next);
            prev = cur;
            cur = next;
        }
        Arc arc;
        arc.parentPath = std::move(path);
        arcs.push_back(std::move(arc));
    };

    for (const auto &kv : adj) {
        const int v = kv.first;
        if (!isNode(v)) continue;
        for (int w : kv.second) {
            if (used.count(EdgeKey(v, w))) continue;
            walk(v, w);
        }
    }

    // A component of G with no node at all is a closed cycle. It cannot arise
    // from ConeCut -- every arc it lays ends on dS or on an earlier arc -- but
    // if it ever did, cutting it open at an arbitrary vertex keeps the pairing
    // well defined rather than losing the component silently.
    for (const auto &kv : adj) {
        const int v = kv.first;
        bool any = false;
        for (int w : kv.second) if (!used.count(EdgeKey(v, w))) any = true;
        if (!any) continue;
        report.messages.push_back(
            "A component of the cutting graph is a closed loop with no junction; it was "
            "opened at an arbitrary vertex to give it a transition.");
        for (int w : kv.second) {
            if (used.count(EdgeKey(v, w))) continue;
            walk(v, w);
        }
    }

    report.arcs = static_cast<int>(arcs.size());
}

// ---------------------------------------------------------------------------
// fitTransitions()  --  Sec. 3.2.2
//
//     theta = arg min over rigid motions of the fitting residual
//     k     = round( theta / (pi/2) )
//     assign the arc to Gamma_Hol_k; store R_k and t
//
// The pairing first. Walking an arc from its first node to its last, each edge
// (a, b) has two incident faces in S, exactly one of which sees the edge as the
// directed a -> b (the mesh is consistently oriented, and ConeCut re-uses that
// orientation for Omega). That one is to the left of the direction of travel
// and is the "+" side; the other is the "-" side. Taking the side from the
// direction of travel rather than from the two faces' indices is what gives the
// consistent arc-length orientation the seam table needs: the pairing at edge
// t and the pairing at edge t+1 then agree at the vertex between them, because
// the faces in the sector between the two edges on one side are glued to each
// other in Omega and share a child.
//
// The fit itself is 2-D Procrustes with the rotation left free: centre both
// polylines, then theta = atan2( sum cross, sum dot ). Only the rotation is
// snapped -- the translation is recovered from the snapped rotation afterwards,
// so it absorbs whatever the snap moved and the residual stays a measure of the
// rotation alone.
// ---------------------------------------------------------------------------
void Immersion::fitTransitions() {
    std::unordered_map<EdgeKey, int, EdgeKeyHash> origEdgeIndex;
    origEdgeIndex.reserve(orig->edges.size() * 2);
    for (int e = 0; e < static_cast<int>(orig->edges.size()); ++e) {
        origEdgeIndex.emplace(EdgeKey(orig->edges[e][0], orig->edges[e][1]), e);
    }

    auto localOf = [&](int f, int v) -> int {
        const Triangle &t = orig->triangles[f];
        if (t[0] == v) return 0;
        if (t[1] == v) return 1;
        if (t[2] == v) return 2;
        return -1;
    };
    // Does face f traverse the edge as a -> b?
    auto directedIn = [&](int f, int a, int b) -> bool {
        const int p = localOf(f, a);
        return p >= 0 && orig->triangles[f][(p + 1) % 3] == b;
    };

    for (size_t ai = 0; ai < arcs.size(); ++ai) {
        Arc &arc = arcs[ai];
        const std::vector<int> &path = arc.parentPath;
        arc.plusChain.clear();
        arc.minusChain.clear();

        std::vector<std::array<int, 4>> pairs;
        std::vector<double> pairLen;
        bool broken = false;

        // The field route's exact transition. Reading it off the frames rather
        // than fitting it is what makes C2 hold by construction: with the same
        // dx on both sides of the seam, dpsi+ = J*_+ dx and dpsi- = J*_- dx, so
        // the rotation taking the minus side to the plus side is
        // J*_+ (J*_-)^-1 = R(theta_hat_- - theta_hat_+), and the combing has
        // already made that a whole number of quarter turns. Which face is
        // which is settled here and nowhere else -- fPlus is the one traversing
        // the arc in the direction of travel -- so this is where the frames are
        // read, and the "+"/"-" of the transition cannot drift away from the
        // "+"/"-" of the pairing.
        int frameK = -1;

        for (size_t t = 0; t + 1 < path.size(); ++t) {
            const int a = path[t], b = path[t + 1];
            auto ie = origEdgeIndex.find(EdgeKey(a, b));
            if (ie == origEdgeIndex.end()) { broken = true; break; }
            const int e = ie->second;
            const int fA = orig->edgeTriangles[e][0];
            const int fB = orig->edgeTriangles[e][1];
            if (fA < 0 || fB < 0) {
                // A cutting-graph edge that lies in dS has one face and no
                // second side to pair with. It is already a boundary of Omega,
                // so it needs no transition; it just cannot be part of a seam.
                ++report.boundaryArcEdges;
                broken = true;
                break;
            }

            int fPlus = -1, fMinus = -1;
            if (directedIn(fA, a, b))      { fPlus = fA; fMinus = fB; }
            else if (directedIn(fB, a, b)) { fPlus = fB; fMinus = fA; }
            else { broken = true; break; }   // inconsistent winding across the edge

            if (!faceAngle.empty()) {
                const double turns = (faceAngle[fMinus] - faceAngle[fPlus]) / M_PI_2;
                const long rounded = std::lround(turns);
                report.maxFrameKResidual = std::max(report.maxFrameKResidual,
                                                    std::fabs(turns - static_cast<double>(rounded)) * M_PI_2);
                const int kk = static_cast<int>(((rounded % 4) + 4) % 4);
                if (frameK < 0) frameK = kk;
                else if (frameK != kk) ++report.frameKConflicts;
            }

            const int pa = cutMesh.triangles[fPlus][localOf(fPlus, a)];
            const int pb = cutMesh.triangles[fPlus][localOf(fPlus, b)];
            const int ma = cutMesh.triangles[fMinus][localOf(fMinus, a)];
            const int mb = cutMesh.triangles[fMinus][localOf(fMinus, b)];

            if (arc.plusChain.empty()) {
                arc.plusChain.push_back(pa);
                arc.minusChain.push_back(ma);
            } else if (arc.plusChain.back() != pa || arc.minusChain.back() != ma) {
                // The sector between this edge and the last one is not glued:
                // an undetected junction. Reported rather than papered over.
                broken = true;
                break;
            }
            arc.plusChain.push_back(pb);
            arc.minusChain.push_back(mb);

            pairs.push_back({pa, pb, ma, mb});
            pairLen.push_back(flatLen[e]);
            arc.length += flatLen[e];
        }

        if (broken || arc.plusChain.size() < 2) {
            arc.degenerate = true;
            ++report.degenerateArcs;
            std::ostringstream oss;
            oss << "Arc " << ai << " of the cutting graph (vertices " << path.front()
                << ".." << path.back() << ") could not be paired into two sides; "
                << "no transition was fitted for it.";
            report.messages.push_back(oss.str());
            continue;
        }

        // The two sides have to be genuinely two sides. If every pair is the
        // same vertex the edge set was cut along an edge that separated nothing.
        bool separated = false;
        for (size_t i = 0; i < arc.plusChain.size(); ++i) {
            if (arc.plusChain[i] != arc.minusChain[i]) { separated = true; break; }
        }
        if (!separated) {
            arc.degenerate = true;
            ++report.degenerateArcs;
            report.messages.push_back("Arc " + std::to_string(ai) +
                                      " has the same vertices on both sides; it did not cut.");
            continue;
        }

        // Procrustes over the corresponding children.
        const size_t n = arc.plusChain.size();
        Point ca{0.0, 0.0}, cb{0.0, 0.0};
        for (size_t i = 0; i < n; ++i) {
            ca = ca + uv[arc.minusChain[i]];
            cb = cb + uv[arc.plusChain[i]];
        }
        ca = ca / static_cast<double>(n);
        cb = cb / static_cast<double>(n);

        double sc = 0.0, sd = 0.0;
        for (size_t i = 0; i < n; ++i) {
            const Point p = uv[arc.minusChain[i]] - ca;
            const Point qq = uv[arc.plusChain[i]] - cb;
            sc += cross2(p, qq);
            sd += dotP(p, qq);
        }
        arc.theta = std::atan2(sc, sd);

        // The fit is always taken. Where the frames supplied a k it is the fit
        // that gives way, and snapError stops being a rounding step and becomes
        // a check on the branch selection: how far the geometry of the map is
        // from the quarter turn the matchings say the arc carries.
        int k = static_cast<int>(std::lround(arc.theta / M_PI_2));
        k = ((k % 4) + 4) % 4;
        if (frameK >= 0) k = frameK;
        arc.k = k;
        arc.snapError = std::fabs(wrap_pi(arc.theta - k * M_PI_2));

        arc.translation = cb - rotateQuarter(ca, k);

        double sum2 = 0.0;
        for (size_t i = 0; i < n; ++i) {
            const Point mapped = rotateQuarter(uv[arc.minusChain[i]], k) + arc.translation;
            sum2 += dotP(mapped - uv[arc.plusChain[i]], mapped - uv[arc.plusChain[i]]);
        }
        arc.fitResidual = std::sqrt(sum2 / static_cast<double>(n));

        for (int cv : {arc.parentPath.front(), arc.parentPath.back()}) {
            if (cutVertCone.empty()) break;
            for (int child : cutter->getOriginalToCutVertices()[cv]) {
                if (cutVertCone[child] >= 0) arc.touchesCone = true;
            }
        }

        ++report.holonomyCount[k];
        report.maxSnapError = std::max(report.maxSnapError, arc.snapError);
        report.maxFitResidual = std::max(report.maxFitResidual, arc.fitResidual);

        for (size_t i = 0; i < pairs.size(); ++i) {
            seamPairs.push_back(pairs[i]);
            seamPairArc.push_back(static_cast<int>(ai));
            seamPairLength.push_back(pairLen[i]);
        }
    }
    report.seamEdgePairs = static_cast<int>(seamPairs.size());
}

// ---------------------------------------------------------------------------
// alignToAxes()
//
// Sec. 3.2.2 closes with the assumption that the immersion has been rotated so
// that at least one boundary edge lies on a line of constant u or constant v.
// Which one is not determined; it matters, because Stage 6's E2 then drags
// every other boundary edge onto the same grid, and starting from the rotation
// that most of the boundary is already near costs far fewer outer steps than
// starting from an arbitrary one.
//
// So the choice is made in two parts. The whole boundary votes, as a
// 4-symmetric field -- sum of l_e exp(4 i theta_e), whose argument over four is
// the rotation that best aligns the boundary as a whole, and which is immune to
// the fact that a boundary edge has no preferred one of its four axis
// directions. Then the single edge whose own direction is nearest that vote is
// aligned *exactly*, which is what the assumption literally asks for.
// ---------------------------------------------------------------------------
void Immersion::alignToAxes() {
    const auto &c2o = cutter->getCutVertexToOriginal();
    std::unordered_map<EdgeKey, int, EdgeKeyHash> origEdgeIndex;
    origEdgeIndex.reserve(orig->edges.size() * 2);
    for (int e = 0; e < static_cast<int>(orig->edges.size()); ++e) {
        origEdgeIndex.emplace(EdgeKey(orig->edges[e][0], orig->edges[e][1]), e);
    }

    double sumRe = 0.0, sumIm = 0.0;
    std::vector<std::pair<int, double>> candidates;   // Omega edge, its direction

    for (int e : cutMesh.boundaryEdges) {
        const int ca = cutMesh.edges[e][0];
        const int cb = cutMesh.edges[e][1];
        auto it = origEdgeIndex.find(EdgeKey(c2o[ca], c2o[cb]));
        if (it == origEdgeIndex.end() || !orig->isBoundaryEdge[it->second]) continue;

        const Point d = uv[cb] - uv[ca];
        const double len = normP(d);
        if (!(len > 0.0)) continue;
        const double th = std::atan2(d[1], d[0]);
        sumRe += len * std::cos(4.0 * th);
        sumIm += len * std::sin(4.0 * th);
        candidates.push_back({e, th});
    }
    if (candidates.empty()) return;

    const double vote = 0.25 * std::atan2(sumIm, sumRe);   // in (-pi/4, pi/4]

    int best = -1;
    double bestGap = std::numeric_limits<double>::infinity();
    double bestRot = 0.0;
    for (const auto &c : candidates) {
        // The rotation that puts this edge exactly on an axis, reduced to the
        // quarter turn it is really asking for.
        double rot = -c.second;
        rot -= M_PI_2 * std::round(rot / M_PI_2);
        const double gap = std::fabs(wrap_pi(4.0 * (rot + vote)) * 0.25);
        if (gap < bestGap) { bestGap = gap; best = c.first; bestRot = rot; }
    }
    if (best < 0) return;

    const double cs = std::cos(bestRot), sn = std::sin(bestRot);
    for (Point &p : uv) {
        const double x = p[0], y = p[1];
        p[0] = cs * x - sn * y;
        p[1] = sn * x + cs * y;
    }
    // A rotation commutes with R_k, so k and the fit residual are untouched and
    // only the translations move with the map.
    for (Arc &a : arcs) {
        const double x = a.translation[0], y = a.translation[1];
        a.translation[0] = cs * x - sn * y;
        a.translation[1] = sn * x + cs * y;
    }

    report.globalRotation = bestRot;
    report.alignedBoundaryEdge = best;
}

// ---------------------------------------------------------------------------
// check()
//
// Everything here is measured on the map that came out, not inferred from the
// steps that made it. Three separate things are being asked:
//
//   * did the unfolding realise the metric      -> maxMetricResidual
//   * is it locally injective, Q1               -> flippedFaces, min area ratio
//   * are the angle sums what Stage 3 asked for -> the two angle residuals, Q2
// ---------------------------------------------------------------------------
void Immersion::check() {
    // The metric, edge by edge.
    for (int e = 0; e < static_cast<int>(cutMesh.edges.size()); ++e) {
        const double want = cutLen[e];
        if (!(want > 0.0)) continue;
        const double got = normP(uv[cutMesh.edges[e][1]] - uv[cutMesh.edges[e][0]]);
        report.maxMetricResidual = std::max(report.maxMetricResidual,
                                            std::fabs(got - want) / want);
    }

    // Q1, and the angle sums for Q2.
    const int nCV = static_cast<int>(cutMesh.vertices.size());
    std::vector<double> angleAt(nCV, 0.0);
    report.minSignedAreaRatio = std::numeric_limits<double>::infinity();

    for (int f = 0; f < static_cast<int>(cutMesh.triangles.size()); ++f) {
        const Triangle &t = cutMesh.triangles[f];
        const Point &p0 = uv[t[0]], &p1 = uv[t[1]], &p2 = uv[t[2]];

        const double area = signedArea(p0, p1, p2);
        const double flat = heronArea(flatLen[orig->triangleEdges[f][0]],
                                      flatLen[orig->triangleEdges[f][1]],
                                      flatLen[orig->triangleEdges[f][2]]);
        if (flat > 0.0) {
            report.minSignedAreaRatio = std::min(report.minSignedAreaRatio, area / flat);
        }
        if (!(area > 0.0)) ++report.flippedFaces;

        angleAt[t[0]] += cornerAngle(p0, p1, p2);
        angleAt[t[1]] += cornerAngle(p1, p2, p0);
        angleAt[t[2]] += cornerAngle(p2, p0, p1);
    }
    if (!std::isfinite(report.minSignedAreaRatio)) report.minSignedAreaRatio = 0.0;

    // Q2. The angle sum of an original vertex is the sum over all its children,
    // which is what Eq. (1) is written to allow; the flat cone metric asks for
    // 2pi (or pi on dS) less the (pi/2) I it was driven to.
    const std::vector<int> &index = cones->getIndices();
    const auto &o2c = cutter->getOriginalToCutVertices();
    const std::vector<char> &active = cones->getActiveVertices();
    for (int v = 0; v < static_cast<int>(orig->vertices.size()); ++v) {
        if (!active[v] || o2c[v].empty()) continue;
        double sum = 0.0;
        for (int cv : o2c[v]) sum += angleAt[cv];
        const double full = orig->isBoundaryVertex[v] ? M_PI : 2.0 * M_PI;
        const double want = full - M_PI_2 * index[v];
        const double res = std::fabs(sum - want);
        if (index[v] != 0) report.maxConeAngleResidual = std::max(report.maxConeAngleResidual, res);
        else report.maxRegularAngleResidual = std::max(report.maxRegularAngleResidual, res);
    }

    if (!uv.empty()) {
        report.uvMin = report.uvMax = uv[0];
        for (const Point &p : uv) {
            report.uvMin[0] = std::min(report.uvMin[0], p[0]);
            report.uvMin[1] = std::min(report.uvMin[1], p[1]);
            report.uvMax[0] = std::max(report.uvMax[0], p[0]);
            report.uvMax[1] = std::max(report.uvMax[1], p[1]);
        }
    }

    // Clustered cones, Sec. 3.1.
    {
        const double extent = std::hypot(report.uvMax[0] - report.uvMin[0],
                                         report.uvMax[1] - report.uvMin[1]);
        double closest = std::numeric_limits<double>::infinity();
        int closestA = -1, closestB = -1;
        for (size_t i = 0; i < coneChildren.size(); ++i) {
            for (size_t j = i + 1; j < coneChildren.size(); ++j) {
                for (int a : coneChildren[i]) {
                    for (int b : coneChildren[j]) {
                        const double d = normP(uv[a] - uv[b]);
                        if (d < closest) { closest = d; closestA = coneVertex[i]; closestB = coneVertex[j]; }
                    }
                }
            }
        }
        if (std::isfinite(closest) && extent > 0.0) {
            report.minConeSeparation = closest / extent;
            for (size_t i = 0; i < coneChildren.size(); ++i) {
                for (size_t j = i + 1; j < coneChildren.size(); ++j) {
                    double d = std::numeric_limits<double>::infinity();
                    for (int a : coneChildren[i]) {
                        for (int b : coneChildren[j]) d = std::min(d, normP(uv[a] - uv[b]));
                    }
                    if (d < 1e-2 * extent) ++report.clusteredConePairs;
                }
            }
            if (report.clusteredConePairs > 0) {
                std::ostringstream oss;
                oss << report.clusteredConePairs << " pair(s) of cones lie within 1% of the "
                    << "model of each other (closest: vertices " << closestA << " and "
                    << closestB << ", " << report.minConeSeparation
                    << " of the extent apart). Sec. 3.1 calls this clustering and merges such a "
                    << "cluster into one cone of the summed index; left alone it shows up as a "
                    << "Q2 residual after Stage 6.";
                report.messages.push_back(oss.str());
            }
        }
    }

    if (report.flippedFaces > 0) {
        std::ostringstream oss;
        oss << report.flippedFaces << " face(s) of Omega are inverted in the immersion; "
            << "Q1 does not hold and Stage 6's barrier cannot repair it.";
        report.messages.push_back(oss.str());
    }
    if (report.maxMetricResidual > 1e-6) {
        std::ostringstream oss;
        if (fromField) {
            // Not a failure and not the same statement. The unfolding realises
            // its metric or it has gone wrong; an integrated map realises the
            // field metric only as far as the field is integrable, and this is
            // that gap measured per edge. It is the number Sec. 7.3 of
            // docs/cf_flow_pipeline.md sends refinement after.
            oss << "The integrated map is off the field metric by up to "
                << report.maxMetricResidual << " relative; that is the "
                << "non-integrability of the field, edge by edge, and not a "
                << "broken layout.";
        } else {
            oss << "The unfolding is off the flat metric by up to "
                << report.maxMetricResidual << " relative.";
        }
        report.messages.push_back(oss.str());
    }
    if (report.frameKConflicts > 0) {
        std::ostringstream oss;
        oss << report.frameKConflicts << " seam edge(s) disagree with the rest of their arc "
            << "about its quarter turn (worst frame residual " << report.maxFrameKResidual
            << " rad). The combing leaked across G somewhere; the arc membership test is "
            << "what to look at, not the field.";
        report.messages.push_back(oss.str());
    }

    report.valid = report.unplacedVertices == 0 && report.flippedFaces == 0 &&
                   report.degenerateFaces == 0 && report.degenerateArcs == 0 &&
                   (fromField ? report.frameKConflicts == 0
                              : report.maxMetricResidual < 1e-6);
}

// ---------------------------------------------------------------------------
bool Immersion::writeOBJ(const std::string &filename) const {
    std::ofstream out(filename);
    if (!out) return false;
    for (const Point &p : uv) out << "v " << p[0] << " " << p[1] << " 0\n";
    for (const Triangle &t : cutMesh.triangles) {
        out << "f " << (t[0] + 1) << " " << (t[1] + 1) << " " << (t[2] + 1) << "\n";
    }
    return true;
}
