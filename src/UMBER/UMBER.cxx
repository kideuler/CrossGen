#include "UMBER.hxx"

#include <algorithm>
#include <cmath>
#include <complex>
#include <deque>
#include <iostream>
#include <limits>
#include <queue>
#include <stdexcept>

#include "dualmbo/DualMBO.hxx"

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
UMBER::UMBER(const DualMBO &dualMBO, const HarmonicCut &harmonicCut) {
    mesh = dualMBO.getMeshPtr();
    if (!mesh) throw std::runtime_error("UMBER: DualMBO carries no mesh");

    if (harmonicCut.getOriginalMeshPtr() != mesh) {
        throw std::runtime_error("UMBER: the HarmonicCut was not built on the DualMBO mesh");
    }

    const int nT = static_cast<int>(mesh->triangles.size());
    if (static_cast<int>(dualMBO.u_k.size()) != nT) {
        throw std::runtime_error("UMBER: the DualMBO field does not match the mesh");
    }

    // u_k[t] = exp(4 i theta_t): any of the four directions is a valid start,
    // combInitialField() picks the branch.
    initialDirections.resize(nT);
    for (int t = 0; t < nT; ++t) {
        const double theta = std::arg(dualMBO.u_k[t]) / 4.0;
        initialDirections[t] = Point{std::cos(theta), std::sin(theta)};
    }

    // C of Eq. (2): the beta cuts of Sec. 4.1, one per void.
    cutEdges = harmonicCut.getCutEdges();
}

UMBER::UMBER(std::shared_ptr<Mesh> meshIn,
             const std::vector<Point> &initialField,
             const std::unordered_set<EdgeKey, EdgeKeyHash> &cuts)
    : mesh(std::move(meshIn)), cutEdges(cuts), initialDirections(initialField) {
    if (!mesh) throw std::runtime_error("UMBER: null mesh");
    if (initialDirections.size() != mesh->triangles.size()) {
        throw std::runtime_error("UMBER: the initial field does not match the mesh");
    }
}

// ---------------------------------------------------------------------------
void UMBER::setFeatureEdges(const std::vector<int> &edges) {
    featureEdges.clear();
    featureEdgeKeys.clear();
    featureEdges.reserve(edges.size());
    for (const int e : edges) {
        if (e < 0 || e >= static_cast<int>(mesh->edges.size())) continue;
        featureEdges.push_back(e);
        featureEdgeKeys.insert(EdgeKey(mesh->edges[e][0], mesh->edges[e][1]));
    }
    initialized = false;   // the quadrature of Eq. (4) has changed
}

// ---------------------------------------------------------------------------
// buildRegions()
//
// The material regions: components of the triangles when the feature edges are
// walls. A vertex strictly inside one gets its label; a vertex whose incident
// triangles disagree sits on an interface and gets -1, which never matches
// anything, so nothing on an interface is ever paired by the dipole rule.
//
// Without feature edges every triangle is in region 0 and every vertex is
// labelled 0, which is the right answer: a single-material model is one
// region, and a +1/-1 pair in it is a pair in one region.
// ---------------------------------------------------------------------------
void UMBER::buildRegions() {
    const int nT = static_cast<int>(mesh->triangles.size());
    const int nV = static_cast<int>(mesh->vertices.size());
    std::vector<int> triRegion(nT, -1);

    std::vector<char> isFeature(mesh->edges.size(), 0);
    for (const int e : featureEdges)
        if (e >= 0 && e < static_cast<int>(isFeature.size())) isFeature[e] = 1;

    int regions = 0;
    std::vector<int> stack;
    for (int t = 0; t < nT; ++t) {
        if (triRegion[t] >= 0) continue;
        const int r = regions++;
        stack.assign(1, t);
        triRegion[t] = r;
        while (!stack.empty()) {
            const int c = stack.back();
            stack.pop_back();
            for (int i = 0; i < 3; ++i) {
                const int nb = mesh->triangleAdjacency[c][i];
                if (nb < 0 || triRegion[nb] >= 0) continue;
                const int e = mesh->triangleEdges[c][i];
                if (e >= 0 && isFeature[e]) continue;
                triRegion[nb] = r;
                stack.push_back(nb);
            }
        }
    }

    vertexRegion.assign(nV, -2);   // -2: not seen yet
    for (int t = 0; t < nT; ++t) {
        for (int i = 0; i < 3; ++i) {
            const int v = mesh->triangles[t][i];
            if (vertexRegion[v] == -2) vertexRegion[v] = triRegion[t];
            else if (vertexRegion[v] != triRegion[t]) vertexRegion[v] = -1;
        }
    }
    for (int v = 0; v < nV; ++v) if (vertexRegion[v] == -2) vertexRegion[v] = -1;
}

// ---------------------------------------------------------------------------
// tau of Eq. (3) and its Jacobian
//
// For a unit v = (cos t, sin t) this is (cos 4t, sin 4t), but it is written as
// a quartic in the components so that the smoothness term on the cuts stays a
// polynomial of the unknowns -- the paper is explicit that the angle is never
// a variable.
// ---------------------------------------------------------------------------
void UMBER::tau(double x, double y, double &t1, double &t2) {
    const double x2 = x * x;
    const double y2 = y * y;
    t1 = x2 * x2 - 6.0 * x2 * y2 + y2 * y2;
    t2 = 4.0 * x2 * x * y - 4.0 * x * y2 * y;
}

void UMBER::tauJacobian(double x, double y, double &p, double &q) {
    const double x2 = x * x;
    const double y2 = y * y;
    // d(tau_1)/dx = d(tau_2)/dy = p,  d(tau_1)/dy = -d(tau_2)/dx = q
    p = 4.0 * x * (x2 - 3.0 * y2);
    q = -4.0 * y * (3.0 * x2 - y2);
}

// ---------------------------------------------------------------------------
// initialize()  --  edge classification, quadrature weights, initial field
// ---------------------------------------------------------------------------
void UMBER::initialize() {
    const int nT = static_cast<int>(mesh->triangles.size());
    if (nT == 0) throw std::runtime_error("UMBER: empty mesh");

    // --- Triangle areas, A_M ------------------------------------------------
    faceArea.assign(nT, 0.0);
    totalArea = 0.0;
    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = mesh->triangles[t];
        const Point &p0 = mesh->vertices[tri[0]];
        const Point &p1 = mesh->vertices[tri[1]];
        const Point &p2 = mesh->vertices[tri[2]];
        const double x10 = p1[0] - p0[0], y10 = p1[1] - p0[1];
        const double x20 = p2[0] - p0[0], y20 = p2[1] - p0[1];
        faceArea[t] = 0.5 * std::fabs(x10 * y20 - x20 * y10);
        totalArea += faceArea[t];
    }
    if (totalArea <= 0.0) throw std::runtime_error("UMBER: mesh has no area");

    // --- Interior edges of Eq. (2), split by membership in C -----------------
    smoothEdges.clear();
    smoothEdges.reserve(mesh->edges.size());
    for (int e = 0; e < static_cast<int>(mesh->edges.size()); ++e) {
        if (mesh->isBoundaryEdge[e]) continue;

        const int fa = mesh->edgeTriangles[e][0];
        const int fb = mesh->edgeTriangles[e][1];
        if (fa < 0 || fb < 0) continue;

        SmoothEdge se;
        se.fa = fa;
        se.fb = fb;
        se.weight = (faceArea[fa] + faceArea[fb]) / totalArea; // A_ab / A_M
        se.onCut = cutEdges.find(EdgeKey(mesh->edges[e][0], mesh->edges[e][1])) != cutEdges.end();
        smoothEdges.push_back(se);
    }

    // --- Feature edges of Eq. (4): dS, and the interfaces where there are any
    //
    // sigma_e is the polar angle of the edge; the l1 norm of R(-sigma_e) v is
    // invariant under sigma_e -> sigma_e + pi/2, so the edge orientation and
    // the choice of which cross direction v represents are both irrelevant
    // here.
    //
    // A boundary edge contributes one term, for its one triangle. An interface
    // edge contributes one per side: both triangles are asked for the same
    // axis, which is what makes the frame agree across the interface rather
    // than merely follow it from one material. See setFeatureEdges.
    alignEdges.clear();
    alignEdges.reserve(mesh->boundaryEdges.size() + 2 * featureEdges.size());
    totalBoundaryLength = 0.0;

    struct Side { int edge; int face; };
    std::vector<Side> sides;
    sides.reserve(mesh->boundaryEdges.size() + 2 * featureEdges.size());
    for (const int e : mesh->boundaryEdges) {
        const int f = mesh->edgeTriangles[e][0] >= 0 ? mesh->edgeTriangles[e][0]
                                                     : mesh->edgeTriangles[e][1];
        if (f >= 0) sides.push_back({e, f});
    }
    for (const int e : featureEdges) {
        if (e < 0 || e >= static_cast<int>(mesh->edges.size())) continue;
        if (mesh->isBoundaryEdge[e]) continue;   // already carried above
        for (int k = 0; k < 2; ++k) {
            const int f = mesh->edgeTriangles[e][k];
            if (f >= 0) sides.push_back({e, f});
        }
    }

    std::vector<double> edgeLen(sides.size(), 0.0);
    for (size_t i = 0; i < sides.size(); ++i) {
        const Point &pa = mesh->vertices[mesh->edges[sides[i].edge][0]];
        const Point &pb = mesh->vertices[mesh->edges[sides[i].edge][1]];
        edgeLen[i] = normP(pb - pa);
        totalBoundaryLength += edgeLen[i];
    }

    for (size_t i = 0; i < sides.size(); ++i) {
        if (edgeLen[i] < 1e-14 || totalBoundaryLength <= 0.0) continue;
        const Point &pa = mesh->vertices[mesh->edges[sides[i].edge][0]];
        const Point &pb = mesh->vertices[mesh->edges[sides[i].edge][1]];
        const double dx = pb[0] - pa[0], dy = pb[1] - pa[1];

        AlignEdge ae;
        ae.face = sides[i].face;
        ae.cosSigma = dx / edgeLen[i];
        ae.sinSigma = dy / edgeLen[i];
        ae.weight = edgeLen[i] / totalBoundaryLength; // l_e / l_dM
        alignEdges.push_back(ae);
    }

    // --- Initial field ------------------------------------------------------
    // The regions first: singularitySeams() needs them to tell a pair that may
    // annihilate from one that may not.
    buildRegions();
    combInitialField();

    x.resize(2 * nT);
    for (int t = 0; t < nT; ++t) {
        x[2 * t + 0] = uField[t][0];
        x[2 * t + 1] = uField[t][1];
    }

    // Report the initial energy at the *final* smoothing, so that it is
    // directly comparable with what optimize() leaves behind: E_align carries
    // an offset of about 2 eps per boundary edge, which would otherwise make
    // the two numbers differ for a reason that has nothing to do with the
    // field.
    l1Eps = l1Schedule.empty() ? 1e-3 : l1Schedule.back();
    Eigen::VectorXd g(2 * nT);
    evaluate(x, g, &currentEnergy);
    currentGradNorm = g.norm();
    totalIterations = 0;
    minFrameNorm_ = 1.0;

    initialized = true;
}

// ---------------------------------------------------------------------------
// inputSingularities()
//
// The spin-4 winding of the *input* directions, before any combing: the
// quarter singularities the cross field carries, as (vertex, quarter turns).
// Measured the way DualMBO::computeSingularities does, on the 4-symmetric
// representation, because that is the only one in which a quarter turn exists.
// ---------------------------------------------------------------------------
std::vector<std::pair<int, int>> UMBER::inputSingularities() const {
    std::vector<std::pair<int, int>> result;
    const int nV = static_cast<int>(mesh->vertices.size());

    for (int v = 0; v < nV; ++v) {
        if (mesh->isBoundaryVertex[v]) continue;

        const auto &vt = mesh->vertexTriangles;
        const int start = vt.rowPtr[v];
        const int end = vt.rowPtr[v + 1];
        const int star = end - start;
        if (star < 2) continue;

        double total = 0.0;
        for (int k = start; k < end; ++k) {
            const int tCur = vt.colIdx[k];
            const int tNext = vt.colIdx[(k - start + 1) % star + start];
            total += wrap_pi(4.0 * (computeAngle(initialDirections[tNext]) -
                                    computeAngle(initialDirections[tCur])));
        }

        const int quarters = static_cast<int>(std::lround(total / (2.0 * M_PI)));
        if (quarters != 0) result.emplace_back(v, quarters);
    }

    return result;
}

// ---------------------------------------------------------------------------
// singularitySeams()
//
// One short path of interior edges from each singularity of the input field to
// the boundary, as an edge set.
//
// The comb below cannot make a quarter turn disappear: around a singularity
// the four cross directions permute, so whatever order the traversal takes,
// some edge is left carrying a 90-degree jump, and the jumps line up along a
// curve from the singularity to somewhere it can end -- another singularity or
// the boundary. That curve is the *branch cut* of the chosen representative,
// and it is what Eq. (1) then has to remove: E_smooth pays for every edge of
// it, so the optimization shortens it, and the singularity comes out of the
// domain at the end where the curve meets the boundary.
//
// Which is why the curve cannot be left to the traversal to place. A plain
// breadth-first comb puts it where its wavefronts happen to collide, which is
// decided by the seed triangle and not by the field. On the half disk in
// data/meshes the DualMBO field carries one quarter singularity either side of
// the centre, symmetrically; the seed lands one branch cut across the whole
// model, from the left singularity through the right one and out of the right
// side of the arc, so both quarter turns leave through that single exit and
// the optimized field ends with its two corners next to each other on the
// lower right instead of one on each side. Nothing downstream can undo that:
// same-sign defects repel, so the pair is a worse minimum than the symmetric
// answer, and L-BFGS has no way down to it.
//
// Routing each cut along a shortest path to the boundary instead makes the
// destination the nearest boundary, which is where a defect would go if it
// could move freely, and makes the cut as short as the mesh allows -- so the
// jump E_smooth has to remove is small to begin with. The paths come from one
// multi-source Dijkstra out of the whole boundary, so two singularities that
// share a watershed share the tail of their cut and the jump there adds, which
// is the correct branch structure for that case rather than a special one.
//
// A cut per singularity is also what the topology requires. Cutting between a
// pair of quarter singularities and stopping there does not work: a loop
// around both has 180 degrees of holonomy, but it crosses such a cut twice, in
// opposite directions, and the two jumps cancel. Only a cut leaving each
// singularity for the boundary gives the loop the one crossing per singularity
// that the holonomy needs.
// ---------------------------------------------------------------------------
std::unordered_set<UMBER::EdgeKey, UMBER::EdgeKeyHash> UMBER::singularitySeams() const {
    std::unordered_set<EdgeKey, EdgeKeyHash> seams;

    const auto sing = inputSingularities();
    if (sing.empty()) return seams;

    const int nV = static_cast<int>(mesh->vertices.size());

    // Primal adjacency, weighted by edge length.
    std::vector<std::vector<std::pair<int, double>>> adj(nV);
    for (const auto &e : mesh->edges) {
        const Point d = mesh->vertices[e[1]] - mesh->vertices[e[0]];
        const double w = normP(d);
        adj[e[0]].emplace_back(e[1], w);
        adj[e[1]].emplace_back(e[0], w);
    }

    // Where a defect is allowed to leave: dS, and -- on a multi-material model
    // -- the interfaces too.
    //
    // An interface is a curve the layout has to follow, so it is a curve the
    // layout is allowed to *turn* on, exactly as dS is, and a defect that
    // reaches one becomes a corner of the interface rather than a corner of
    // the boundary. Routing to dS alone is what forces it past the interface
    // and out to the far edge of the model, and that is not a slightly worse
    // placement: on data/meshes/multimat/geom001 the +1 the field puts inside
    // the quarter disk is the whole of that region's deficit -- its boundary
    // makes three quarter turns and needs four -- so its destination is the
    // middle of the arc and nowhere else. Sent to dS instead it invents a
    // corner on the square and leaves the arc with no corner at all, which is
    // a quarter disk that cannot be one quadrilateral.
    std::vector<char> isExit(nV, 0);
    for (const int v : mesh->boundaryVertices) isExit[v] = 1;
    for (const int e : featureEdges) {
        if (e < 0 || e >= static_cast<int>(mesh->edges.size())) continue;
        isExit[mesh->edges[e][0]] = 1;
        isExit[mesh->edges[e][1]] = 1;
    }

    // Multi-source Dijkstra out of those: dist[v] is the distance to the
    // nearest one and prev[v] the next step towards it.
    std::vector<double> dist(nV, std::numeric_limits<double>::max());
    std::vector<int> prev(nV, -1);
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>,
                        std::greater<std::pair<double, int>>> queue;
    for (int v = 0; v < nV; ++v) {
        if (!isExit[v]) continue;
        dist[v] = 0.0;
        queue.emplace(0.0, v);
    }
    while (!queue.empty()) {
        const auto [d, v] = queue.top();
        queue.pop();
        if (d > dist[v]) continue;
        for (const auto &[nb, w] : adj[v]) {
            if (dist[v] + w < dist[nb]) {
                dist[nb] = dist[v] + w;
                prev[nb] = v;
                queue.emplace(dist[nb], nb);
            }
        }
    }

    // --- The pairs that annihilate rather than leave ------------------------
    //
    // A +1/-1 pair inside one material region contributes nothing to that
    // region's count and is a property of the smoothest field, not of the
    // layout -- so the cut between them is what it wants, and the two defects
    // meet and vanish instead of inventing two corners on dS. See
    // setCancelDipoles for the topology and for the rule that a pair spanning
    // an interface is never one of these.
    //
    // Greedy closest-first over the shortest paths between them, which is
    // MERIDIAN::ConeSingularities::cancelDipoles' own order: what goes is the
    // tightest cluster.
    //
    // Only on a multi-material model, which is MERIDIAN's own gate
    // (MERIDIAN.cxx runs cancelDipoles under multiMaterial()) and for its
    // reason: the pair this removes is the one interface alignment put there.
    // On a single-material model a +1/-1 pair is not that artefact, and
    // cancelling it measurably loses -- over data/meshes/singlemat it takes
    // geom007 from one flipped triangle in the parameterization to three and
    // geom010 from eleven to thirteen, and moves no model the other way.
    dipoles.clear();
    std::vector<char> paired(sing.size(), 0);
    if (cancelDipoles && !featureEdges.empty() && sing.size() >= 2 && !vertexRegion.empty()) {
        // One Dijkstra per singularity, over the same primal graph. There are
        // a handful of these on any model in data/meshes, so the quadratic
        // pairing below costs nothing worth avoiding.
        const int n = static_cast<int>(sing.size());
        std::vector<std::vector<double>> dTo(n);
        std::vector<std::vector<int>> pTo(n);
        for (int i = 0; i < n; ++i) {
            dTo[i].assign(nV, std::numeric_limits<double>::max());
            pTo[i].assign(nV, -1);
            std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>,
                                std::greater<std::pair<double, int>>> q;
            dTo[i][sing[i].first] = 0.0;
            q.emplace(0.0, sing[i].first);
            while (!q.empty()) {
                const auto [d, v] = q.top();
                q.pop();
                if (d > dTo[i][v]) continue;
                for (const auto &[nb, w] : adj[v]) {
                    if (dTo[i][v] + w < dTo[i][nb]) {
                        dTo[i][nb] = dTo[i][v] + w;
                        pTo[i][nb] = v;
                        q.emplace(dTo[i][nb], nb);
                    }
                }
            }
        }

        struct Candidate { double d; int i; int j; };
        std::vector<Candidate> cands;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                // Only a single quarter against a single anti-quarter: a
                // higher-order defect would need as many jumps as its index
                // and one cut cannot carry them.
                if (sing[i].second + sing[j].second != 0) continue;
                if (std::abs(sing[i].second) != 1) continue;
                const int ri = vertexRegion[sing[i].first];
                const int rj = vertexRegion[sing[j].first];
                if (ri < 0 || ri != rj) continue;   // on, or across, an interface
                const double d = dTo[i][sing[j].first];
                if (!(d < std::numeric_limits<double>::max())) continue;
                cands.push_back({d, i, j});
            }
        }
        std::sort(cands.begin(), cands.end(),
                  [](const Candidate &a, const Candidate &b) { return a.d < b.d; });
        for (const Candidate &c : cands) {
            if (paired[c.i] || paired[c.j]) continue;
            int cur = sing[c.j].first;
            int guard = 0;
            bool ok = true;
            std::vector<EdgeKey> path;
            while (cur != sing[c.i].first && guard++ < nV) {
                const int nxt = pTo[c.i][cur];
                if (nxt < 0) { ok = false; break; }
                path.emplace_back(cur, nxt);
                cur = nxt;
            }
            if (!ok || path.empty()) continue;
            for (const EdgeKey &k : path) seams.insert(k);
            paired[c.i] = paired[c.j] = 1;
            dipoles.emplace_back(sing[c.i].first, sing[c.j].first);
        }
    }

    // --- Everything else keeps its route to dS ------------------------------
    for (size_t i = 0; i < sing.size(); ++i) {
        if (paired[i]) continue;
        int cur = sing[i].first;
        int guard = 0;
        while (cur >= 0 && !isExit[cur] && guard++ < nV) {
            const int nxt = prev[cur];
            if (nxt < 0) break; // unreachable: leave this one to the traversal
            seams.insert(EdgeKey(cur, nxt));
            cur = nxt;
        }
    }

    return seams;
}

// ---------------------------------------------------------------------------
// combInitialField()
//
// The input field is 4-symmetric, so each triangle's direction is only defined
// modulo 90 degrees; the non-symmetric energy of Eq. (2) needs one consistent
// choice. Breadth-first over the dual graph, rotating each newly reached
// triangle by the multiple of 90 degrees closest to its parent, does that --
// and stopping at the cut edges is what leaves the transitions there free.
//
// The traversal also stops at the branch cuts of singularitySeams(), so that
// the 90-degree jumps the input's quarter singularities force lie on those
// short paths to the boundary rather than wherever the wavefront closed. See
// singularitySeams() for why that placement decides where the optimized field
// puts its boundary corners.
// ---------------------------------------------------------------------------
void UMBER::combInitialField() {
    const int nT = static_cast<int>(mesh->triangles.size());
    uField = initialDirections;

    const std::unordered_set<EdgeKey, EdgeKeyHash> seams = singularitySeams();

    // Dual adjacency restricted to the non-cut edges.
    std::vector<std::array<int, 3>> combAdj(nT, std::array<int, 3>{-1, -1, -1});
    for (int t = 0; t < nT; ++t) {
        for (int e = 0; e < 3; ++e) {
            const int nb = mesh->triangleAdjacency[t][e];
            if (nb < 0) continue;
            const int edgeIdx = mesh->triangleEdges[t][e];
            if (edgeIdx < 0) continue;
            const EdgeKey ek(mesh->edges[edgeIdx][0], mesh->edges[edgeIdx][1]);
            if (cutEdges.find(ek) != cutEdges.end()) continue; // transition stays free
            if (seams.find(ek) != seams.end()) continue;       // branch cut of the comb
            combAdj[t][e] = nb;
        }
    }

    std::vector<char> visited(nT, 0);
    std::deque<int> queue;
    for (int seed = 0; seed < nT; ++seed) {
        if (visited[seed]) continue;

        // A new seed means the non-cut dual graph is disconnected (the cuts
        // separate a piece of the domain). Each component is combed on its own
        // and the mismatch between them is left to the optimization.
        visited[seed] = 1;
        queue.push_back(seed);

        while (!queue.empty()) {
            const int cur = queue.front();
            queue.pop_front();
            const double thetaCur = computeAngle(uField[cur]);

            for (int e = 0; e < 3; ++e) {
                const int nb = combAdj[cur][e];
                if (nb < 0 || visited[nb]) continue;
                visited[nb] = 1;

                const double thetaNb = computeAngle(uField[nb]);
                const int k = find_rotation_matrix(thetaNb, thetaCur);
                uField[nb] = rotateVector(uField[nb], k);

                queue.push_back(nb);
            }
        }
    }

    vField.resize(nT);
    angles.resize(nT);
    for (int t = 0; t < nT; ++t) {
        vField[t] = rotateVector(uField[t], 1);
        angles[t] = computeAngle(uField[t]);
    }
}

// ---------------------------------------------------------------------------
// evaluate()  --  Eq. (1) and its gradient
// ---------------------------------------------------------------------------
double UMBER::evaluate(const Eigen::VectorXd &v, Eigen::VectorXd &grad,
                       EnergyTerms *terms) const {
    const int nT = static_cast<int>(mesh->triangles.size());
    grad.setZero(2 * nT);

    double eSmooth = 0.0, eAlign = 0.0, eReg = 0.0;

    // --- Smoothness, Eq. (2) ------------------------------------------------
    for (const SmoothEdge &se : smoothEdges) {
        const int ia = 2 * se.fa, ib = 2 * se.fb;
        const double xa = v[ia], ya = v[ia + 1];
        const double xb = v[ib], yb = v[ib + 1];

        if (!se.onCut) {
            // ||v_a - v_b||^2
            const double d1 = xa - xb, d2 = ya - yb;
            eSmooth += se.weight * (d1 * d1 + d2 * d2);

            const double c = 2.0 * se.weight;
            grad[ia]     += c * d1;
            grad[ia + 1] += c * d2;
            grad[ib]     -= c * d1;
            grad[ib + 1] -= c * d2;
        } else {
            // ||tau(v_a) - tau(v_b)||^2: blind to a k*90-degree transition, so
            // the harmonic degrees of freedom of Sec. 4.1 survive here.
            double ta1, ta2, tb1, tb2;
            tau(xa, ya, ta1, ta2);
            tau(xb, yb, tb1, tb2);
            const double d1 = ta1 - tb1, d2 = ta2 - tb2;
            eSmooth += se.weight * (d1 * d1 + d2 * d2);

            double pa, qa, pb, qb;
            tauJacobian(xa, ya, pa, qa);
            tauJacobian(xb, yb, pb, qb);

            // J = [[p, q], [-q, p]], so J^T d = (p*d1 - q*d2, q*d1 + p*d2).
            const double c = 2.0 * se.weight;
            grad[ia]     += c * (pa * d1 - qa * d2);
            grad[ia + 1] += c * (qa * d1 + pa * d2);
            grad[ib]     -= c * (pb * d1 - qb * d2);
            grad[ib + 1] -= c * (qb * d1 + pb * d2);
        }
    }

    // --- Boundary alignment, Eq. (4) ---------------------------------------
    // r = R(-sigma) v; the smoothed |r_1| + |r_2| is minimal exactly when r
    // lies on an axis of the edge frame, i.e. when v is parallel or orthogonal
    // to the edge.
    const double eps2 = l1Eps * l1Eps;
    for (const AlignEdge &ae : alignEdges) {
        const int i = 2 * ae.face;
        const double vx = v[i], vy = v[i + 1];

        const double r1 =  ae.cosSigma * vx + ae.sinSigma * vy;
        const double r2 = -ae.sinSigma * vx + ae.cosSigma * vy;

        const double s1 = std::sqrt(r1 * r1 + eps2);
        const double s2 = std::sqrt(r2 * r2 + eps2);
        eAlign += ae.weight * (s1 + s2);

        // dE/dv = R(-sigma)^T dE/dr
        const double dr1 = ae.weight * r1 / s1;
        const double dr2 = ae.weight * r2 / s2;
        grad[i]     += alignWeight * ( ae.cosSigma * dr1 - ae.sinSigma * dr2);
        grad[i + 1] += alignWeight * ( ae.sinSigma * dr1 + ae.cosSigma * dr2);
    }

    // --- Unit length regularization, Eq. (5) --------------------------------
    for (int t = 0; t < nT; ++t) {
        const int i = 2 * t;
        const double vx = v[i], vy = v[i + 1];
        const double n2 = vx * vx + vy * vy;
        const double d = n2 - 1.0;
        const double w = faceArea[t] / totalArea;

        eReg += w * d * d;

        const double c = regWeight * 4.0 * w * d;
        grad[i]     += c * vx;
        grad[i + 1] += c * vy;
    }

    const double total = eSmooth + alignWeight * eAlign + regWeight * eReg;
    if (terms) {
        terms->smooth = eSmooth;
        terms->align = eAlign;
        terms->reg = eReg;
        terms->total = total;
    }
    return total;
}

// ---------------------------------------------------------------------------
// runLBFGS()  --  limited-memory BFGS with a backtracking Armijo line search
// ---------------------------------------------------------------------------
int UMBER::runLBFGS(Eigen::VectorXd &state, int maxIter) {
    const int m = 10;               // history depth
    const double c1 = 1e-4;         // Armijo constant
    const int maxBacktracks = 60;

    const int n = static_cast<int>(state.size());
    Eigen::VectorXd g(n), gNew(n), d(n), xNew(n);

    double f = evaluate(state, g);
    double gnorm = g.norm();

    std::deque<Eigen::VectorXd> S, Y;
    std::deque<double> rho;

    std::vector<double> alpha(m, 0.0);

    int iter = 0;
    for (; iter < maxIter; ++iter) {
        if (gnorm <= gradTolerance) break;

        // --- Two-loop recursion for the search direction --------------------
        Eigen::VectorXd q = g;
        const int k = static_cast<int>(S.size());
        for (int i = k - 1; i >= 0; --i) {
            alpha[i] = rho[i] * S[i].dot(q);
            q -= alpha[i] * Y[i];
        }
        if (k > 0) {
            const double yy = Y[k - 1].squaredNorm();
            if (yy > 1e-300) q *= S[k - 1].dot(Y[k - 1]) / yy;
        }
        for (int i = 0; i < k; ++i) {
            const double beta = rho[i] * Y[i].dot(q);
            q += S[i] * (alpha[i] - beta);
        }
        d = -q;

        double dg = d.dot(g);
        if (!(dg < 0.0)) { // curvature lost: fall back to steepest descent
            d = -g;
            dg = -g.squaredNorm();
        }

        // --- Line search ----------------------------------------------------
        // The first step is scaled by 1/||g||: the raw gradient of Eq. (1) can
        // be large, and a unit step along it lands far outside the region
        // where the quartic tau term is meaningful.
        double step = (k == 0) ? std::min(1.0, 1.0 / std::max(gnorm, 1e-12)) : 1.0;

        bool progressed = false;
        double fNew = f;
        for (int bt = 0; bt < maxBacktracks; ++bt) {
            xNew = state + step * d;
            fNew = evaluate(xNew, gNew);
            if (std::isfinite(fNew) && fNew <= f + c1 * step * dg) {
                progressed = true;
                break;
            }
            step *= 0.5;
        }
        if (!progressed) break; // no decrease available: converged or stuck

        // --- Curvature pair -------------------------------------------------
        Eigen::VectorXd s = xNew - state;
        Eigen::VectorXd y = gNew - g;
        const double sy = s.dot(y);
        if (sy > 1e-12) {
            if (static_cast<int>(S.size()) == m) {
                S.pop_front();
                Y.pop_front();
                rho.pop_front();
            }
            S.push_back(std::move(s));
            Y.push_back(std::move(y));
            rho.push_back(1.0 / sy);
        }

        state = xNew;
        f = fNew;
        g = gNew;
        gnorm = g.norm();
    }

    currentGradNorm = gnorm;
    return iter;
}

// ---------------------------------------------------------------------------
// optimize()  --  Eq. (1) over the l1 smoothing schedule
// ---------------------------------------------------------------------------
void UMBER::optimize() {
    if (!initialized) initialize();

    const int nT = static_cast<int>(mesh->triangles.size());
    totalIterations = 0;

    const std::vector<double> schedule = l1Schedule.empty()
        ? std::vector<double>{1e-3}
        : l1Schedule;

    for (double eps : schedule) {
        l1Eps = eps;
        totalIterations += runLBFGS(x, maxIterations);
    }

    // Publish the frame. E_reg only penalizes the deviation from unit length,
    // so the result is normalized here. A vector that vanished outright has no
    // direction left to normalize; the initial one stands in so that the
    // output is still a well-formed field, and the substitution is reported,
    // because a field with no direction at a triangle is a failure of Eq. (1),
    // not a detail.
    uField.resize(nT);
    vField.resize(nT);
    angles.resize(nT);
    int degenerate = 0;
    minFrameNorm_ = std::numeric_limits<double>::max();
    for (int t = 0; t < nT; ++t) {
        Point v{x[2 * t + 0], x[2 * t + 1]};
        const double n = normP(v);
        minFrameNorm_ = std::min(minFrameNorm_, n);
        if (n < 1e-8) {
            ++degenerate;
            v = initialDirections[t];
        } else {
            v = v / n;
        }
        uField[t] = v;
        vField[t] = rotateVector(v, 1);
        angles[t] = computeAngle(v);
    }
    if (degenerate > 0) {
        std::cerr << "UMBER: " << degenerate << " triangle(s) ended with a vanishing frame vector; "
                  << "the field is degenerate there.\n";
    }

    Eigen::VectorXd g(2 * nT);
    evaluate(x, g, &currentEnergy); // energy of the result
}

// ---------------------------------------------------------------------------
// sharedEdge()  --  the edge two triangles of a vertex star have in common
// ---------------------------------------------------------------------------
int UMBER::sharedEdge(int ta, int tb) const {
    for (int i = 0; i < 3; ++i) {
        const int e = mesh->triangleEdges[ta][i];
        for (int j = 0; j < 3; ++j) {
            if (mesh->triangleEdges[tb][j] == e) return e;
        }
    }
    return -1;
}

// ---------------------------------------------------------------------------
// internalSingularities()
//
// The holonomy of v itself, not of tau(v).
//
// That distinction is the whole point of Sec. 4.2 and it is easy to get wrong.
// A cross field carries its singularities in the 4-symmetry: the natural test
// is the spin-4 winding, which is what DualMBO::computeSingularities does and
// what the input field needs. But Eq. (1) does not optimize a cross field --
// it optimizes a plain vector field, precisely so that a singularity has
// nowhere to hide. A continuous, non-vanishing v has winding 0 around every
// interior vertex, and then the frame (v, v^perp) it spans has no internal
// singularity, full stop.
//
// Measuring the optimized field with the spin-4 test instead reports
// singularities that are not there. wrap(4 d) cannot tell a 46-degree turn
// from a -44-degree one, so a single interior edge where the frame turns a
// little past 45 degrees registers as a quarter turn. On data/meshes the
// counts from the two tests are 0-vs-4, 0-vs-1, 1-vs-7, 0-vs-1, 2-vs-14: the
// spin-4 number tracks sharpTurns().size() almost exactly, and the field is in
// fact clean. See sharpTurns() for why those edges are still worth knowing
// about.
//
// The k*90-degree transitions on the cuts are undone before measuring -- they
// are free by construction (Sec. 4.1) and carry the harmonic part of the
// field, not a defect. Everywhere else the jump is taken as it is, since
// hiding a large one is exactly the error described above.
//
// The index reported is the winding of v, so 1 is a full turn. In cross terms
// that is four quarter-singularities at one vertex; a non-symmetric field
// cannot express a lone quarter turn at all, which is another way of saying
// what Sec. 4.2 buys.
// ---------------------------------------------------------------------------
std::vector<std::pair<int, double>> UMBER::internalSingularities() const {
    std::vector<std::pair<int, double>> result;
    if (uField.empty()) return result;

    const int nV = static_cast<int>(mesh->vertices.size());
    for (int v = 0; v < nV; ++v) {
        if (mesh->isBoundaryVertex[v]) continue;

        const auto &vt = mesh->vertexTriangles;
        const int start = vt.rowPtr[v];
        const int end = vt.rowPtr[v + 1];
        const int star = end - start;
        if (star < 2) continue;

        double total = 0.0;
        for (int k = start; k < end; ++k) {
            const int tCur = vt.colIdx[k];
            const int tNext = vt.colIdx[(k - start + 1) % star + start];

            double thetaNext = angles[tNext];

            // On a cut, rotate the far side back into the near side's frame:
            // that rotation is Pi_gamma of Eq. (7), a degree of freedom rather
            // than a jump.
            const int e = sharedEdge(tCur, tNext);
            if (e >= 0 &&
                cutEdges.find(EdgeKey(mesh->edges[e][0], mesh->edges[e][1])) != cutEdges.end()) {
                thetaNext += find_rotation_matrix(thetaNext, angles[tCur]) * M_PI_2;
            }

            total += wrap_pi(thetaNext - angles[tCur]);
        }

        const int winding = static_cast<int>(std::lround(total / (2.0 * M_PI)));
        if (winding != 0) result.emplace_back(v, static_cast<double>(winding));
    }

    return result;
}

// ---------------------------------------------------------------------------
// boundarySingularities()
//
// Where internalSingularities() counts what failed to leave, this counts what
// arrived: the corners of the polysquare the optimized frame induces.
//
// The measurement is local and needs no global orientation of the boundary
// loops. Around a boundary vertex the frame turns by Theta and the domain
// itself subtends the interior angle Omega; a frame that follows the boundary
// axis-for-axis has Theta - (Omega - pi) equal to a multiple of 90 degrees,
// and that multiple is the corner. The star is walked counter-clockwise so
// that the sign of Theta and the sign of pi - Omega agree -- the direction is
// read off the signed area of the first triangle in the star rather than from
// the boundary loop, which keeps holes and the outer loop on the same footing.
//
// What that quantity must not be is rounded vertex by vertex. A corner is a
// point defect of the continuum but a discrete field spreads it over the two
// or three vertices it takes the frame to swing across, and each of those
// carries a fraction of the quarter turn: on the circle-like model of
// data/meshes/singlemat/geom016 all four corners come out as triples summing to 1.00 --
// 0.35, 0.45, 0.20 at one of them -- so rounding each vertex on its own
// reports every one of them as nothing and the field appears to have lost the
// four defects it in fact still has.
//
// Summing first is what fixes that, and it is exact rather than a heuristic.
// Over a boundary loop the unrounded quantity telescopes: the Theta cancel in
// pairs and the (pi - Omega) add up to the total turning, so the loop sums to
// 4 * chi whatever the field does -- 3.990 measured on that same model. So the
// loop is walked in order with a running total and a quarter turn is handed to
// whichever vertex carries it past the halfway mark. Isolated corners come out
// exactly as before, smeared ones land on their middle vertex, and the count
// adds up to the Euler characteristic by construction instead of by luck.
// ---------------------------------------------------------------------------
std::vector<std::pair<int, int>> UMBER::boundarySingularities() const {
    std::vector<std::pair<int, int>> result;
    if (uField.empty()) return result;

    const int nV = static_cast<int>(mesh->vertices.size());

    // vertex -> its boundary edges. Two is the manifold case; anything else is
    // a pinch and is left alone below.
    std::vector<std::array<int, 2>> incidentBoundary(nV, std::array<int, 2>{-1, -1});
    std::vector<int> boundaryDegree(nV, 0);
    for (int e : mesh->boundaryEdges) {
        for (int k = 0; k < 2; ++k) {
            const int v = mesh->edges[e][k];
            if (boundaryDegree[v] < 2) incidentBoundary[v][boundaryDegree[v]] = e;
            ++boundaryDegree[v];
        }
    }

    auto otherEnd = [&](int e, int v) {
        return (mesh->edges[e][0] == v) ? mesh->edges[e][1] : mesh->edges[e][0];
    };

    // --- Pass 1: the unrounded quarter turns, vertex by vertex --------------
    std::vector<double> raw(nV, 0.0);
    std::vector<char> measured(nV, 0);

    for (int v : mesh->boundaryVertices) {
        if (boundaryDegree[v] != 2) continue;

        // --- The star, walked from one boundary edge to the other -----------
        // The CSR star cannot be used here: it is sorted by absolute angle, so
        // a fan that straddles the -pi/pi cut comes out in the wrong order,
        // which is harmless for a closed interior loop and fatal for an open
        // boundary one.
        const int firstEdge = incidentBoundary[v][0];
        int entry = firstEdge;
        int cur = (mesh->edgeTriangles[entry][0] >= 0) ? mesh->edgeTriangles[entry][0]
                                                       : mesh->edgeTriangles[entry][1];
        const int guard = mesh->vertexTriangles.vertexDegree(v) + 1;

        std::vector<int> fan;
        int exitEdge = -1;
        while (cur >= 0) {
            fan.push_back(cur);

            int next = -1;
            for (int i = 0; i < 3; ++i) {
                const int e = mesh->triangleEdges[cur][i];
                if (e < 0 || e == entry) continue;
                if (mesh->edges[e][0] == v || mesh->edges[e][1] == v) { next = e; break; }
            }
            if (next < 0) { fan.clear(); break; }
            if (mesh->isBoundaryEdge[next]) { exitEdge = next; break; }

            const int nb = (mesh->edgeTriangles[next][0] == cur) ? mesh->edgeTriangles[next][1]
                                                                 : mesh->edgeTriangles[next][0];
            if (nb < 0 || static_cast<int>(fan.size()) > guard) { fan.clear(); break; }
            entry = next;
            cur = nb;
        }
        if (fan.empty() || exitEdge < 0) continue;

        // --- Interior angle Omega -------------------------------------------
        double omega = 0.0;
        for (int t : fan) {
            const Triangle &tri = mesh->triangles[t];
            int i0 = -1;
            for (int i = 0; i < 3; ++i) if (tri[i] == v) i0 = i;
            if (i0 < 0) continue;
            const Point a = mesh->vertices[tri[(i0 + 1) % 3]] - mesh->vertices[v];
            const Point b = mesh->vertices[tri[(i0 + 2) % 3]] - mesh->vertices[v];
            omega += std::fabs(std::atan2(cross2(a, b), dotP(a, b)));
        }

        // --- Orient the walk counter-clockwise ------------------------------
        // In the first triangle of the star, the counter-clockwise sweep at v
        // runs from the edge to the next vertex to the edge to the previous
        // one when the triangle is positively oriented, and the other way when
        // it is not. The walk started on `firstEdge`, so comparing the two
        // says which way it went.
        {
            const Triangle &tri = mesh->triangles[fan.front()];
            int i0 = -1;
            for (int i = 0; i < 3; ++i) if (tri[i] == v) i0 = i;
            if (i0 < 0) continue;
            const int va = tri[(i0 + 1) % 3];
            const int vb = tri[(i0 + 2) % 3];
            const double s = cross2(mesh->vertices[va] - mesh->vertices[v],
                                    mesh->vertices[vb] - mesh->vertices[v]);
            const int w = otherEnd(firstEdge, v);
            const bool ccw = (w == va) ? (s > 0.0) : (s < 0.0);
            if (!ccw) std::reverse(fan.begin(), fan.end());
        }

        // --- Rotation of the frame across the star --------------------------
        double theta = 0.0;
        for (size_t i = 0; i + 1 < fan.size(); ++i) {
            const int tCur = fan[i];
            const int tNext = fan[i + 1];

            double thetaNext = angles[tNext];
            const int e = sharedEdge(tCur, tNext);
            if (e >= 0 &&
                cutEdges.find(EdgeKey(mesh->edges[e][0], mesh->edges[e][1])) != cutEdges.end()) {
                thetaNext += find_rotation_matrix(thetaNext, angles[tCur]) * M_PI_2;
            }

            theta += wrap_pi(thetaNext - angles[tCur]);
        }

        raw[v] = (theta + M_PI - omega) / M_PI_2;
        measured[v] = 1;
    }

    // --- Pass 2: hand the quarter turns out along each boundary loop --------
    // A running total, and the vertex that carries it past the next halfway
    // mark gets the corner. Equivalent to rounding vertex by vertex wherever a
    // corner is sharp, and correct where it is not.
    std::vector<char> visited(nV, 0);

    for (int seed : mesh->boundaryVertices) {
        if (visited[seed] || boundaryDegree[seed] != 2) continue;

        double running = 0.0;
        long handedOut = 0;
        int v = seed;
        int prevEdge = -1;

        while (v >= 0 && !visited[v]) {
            visited[v] = 1;
            if (measured[v]) running += raw[v];

            const long want = std::lround(running);
            if (want != handedOut) {
                result.emplace_back(v, static_cast<int>(want - handedOut));
                handedOut = want;
            }

            // Step to the next vertex of the loop. A pinch (more than two
            // boundary edges) has no unambiguous next, so the walk stops and
            // whatever is left of `running` is dropped with it.
            if (boundaryDegree[v] != 2) break;
            const int e = (incidentBoundary[v][0] != prevEdge) ? incidentBoundary[v][0]
                                                              : incidentBoundary[v][1];
            prevEdge = e;
            v = otherEnd(e, v);
        }
    }

    return result;
}

// ---------------------------------------------------------------------------
// sectorFan() / sectorMeasure()
//
// The star walk boundarySingularities() does, written once and against a pair
// of feature edges rather than against dS, so that the interfaces are measured
// by the same code and cannot drift from it.
//
// The sweep starts in the triangle on the left of the directed edge
// (far end of eIn) -> v, and steps from triangle to triangle through the
// edges at v until it arrives at eOut. Starting on the left is the whole of
// what picks the side: an interface has a fan on each side of it, and which
// one is measured is which way the caller is walking the chain.
// ---------------------------------------------------------------------------
std::vector<int> UMBER::sectorFan(int v, int eIn, int eOut) const {
    std::vector<int> fan;
    if (eIn < 0 || eOut < 0) return fan;
    const int w = (mesh->edges[eIn][0] == v) ? mesh->edges[eIn][1] : mesh->edges[eIn][0];

    // The triangle with the sector on its left when w -> v is walked.
    int cur = -1;
    for (int k = 0; k < 2; ++k) {
        const int t = mesh->edgeTriangles[eIn][k];
        if (t < 0) continue;
        const Triangle &tri = mesh->triangles[t];
        int c = -1;
        for (int i = 0; i < 3; ++i) if (tri[i] != v && tri[i] != w) c = tri[i];
        if (c < 0) continue;
        if (cross2(mesh->vertices[v] - mesh->vertices[w],
                   mesh->vertices[c] - mesh->vertices[w]) > 0.0) { cur = t; break; }
    }
    if (cur < 0) return fan;

    int entry = eIn;
    const int guard = mesh->vertexTriangles.vertexDegree(v) + 1;
    while (cur >= 0) {
        fan.push_back(cur);
        int next = -1;
        for (int i = 0; i < 3; ++i) {
            const int e = mesh->triangleEdges[cur][i];
            if (e < 0 || e == entry) continue;
            if (mesh->edges[e][0] == v || mesh->edges[e][1] == v) { next = e; break; }
        }
        if (next < 0) { fan.clear(); break; }
        if (next == eOut) return fan;

        const int nb = (mesh->edgeTriangles[next][0] == cur) ? mesh->edgeTriangles[next][1]
                                                             : mesh->edgeTriangles[next][0];
        if (nb < 0 || static_cast<int>(fan.size()) > guard) { fan.clear(); break; }
        entry = next;
        cur = nb;
    }
    fan.clear();   // the sweep ran off the model without reaching eOut
    return fan;
}

void UMBER::sectorMeasure(int v, const std::vector<int> &fan, double &omega,
                          double &theta) const {
    omega = 0.0;
    theta = 0.0;
    for (const int t : fan) {
        const Triangle &tri = mesh->triangles[t];
        int i0 = -1;
        for (int i = 0; i < 3; ++i) if (tri[i] == v) i0 = i;
        if (i0 < 0) continue;
        const Point a = mesh->vertices[tri[(i0 + 1) % 3]] - mesh->vertices[v];
        const Point b = mesh->vertices[tri[(i0 + 2) % 3]] - mesh->vertices[v];
        omega += std::fabs(std::atan2(cross2(a, b), dotP(a, b)));
    }
    for (size_t i = 0; i + 1 < fan.size(); ++i) {
        const int tCur = fan[i], tNext = fan[i + 1];
        double thetaNext = angles[tNext];
        const int e = sharedEdge(tCur, tNext);
        if (e >= 0 &&
            cutEdges.find(EdgeKey(mesh->edges[e][0], mesh->edges[e][1])) != cutEdges.end()) {
            thetaNext += find_rotation_matrix(thetaNext, angles[tCur]) * M_PI_2;
        }
        theta += wrap_pi(thetaNext - angles[tCur]);
    }
}

// ---------------------------------------------------------------------------
// interfaceCorners()
//
// The interface network cut into chains -- a run of interface edges between
// two vertices where the network is not simply two edges meeting -- and each
// chain walked in both directions, since an interface has a layout on each
// side of it and the two need not turn the same way.
//
// The accumulation is boundarySingularities()'s, for its reason: a corner is
// spread over the vertices it takes the frame to swing across, and rounding
// each on its own loses it. A chain is open rather than cyclic, so the
// telescoping identity that makes the boundary sum to 4 chi does not apply
// here and the running total is a smoothing rather than an exact count.
// ---------------------------------------------------------------------------
void UMBER::buildInterfaceChains() const {
    chainsBuilt = true;
    chains.clear();
    if (featureEdges.empty()) return;

    const int nV = static_cast<int>(mesh->vertices.size());
    std::vector<std::vector<int>> atVertex(nV);
    for (const int e : featureEdges) {
        if (e < 0 || e >= static_cast<int>(mesh->edges.size())) continue;
        if (mesh->isBoundaryEdge[e]) continue;
        atVertex[mesh->edges[e][0]].push_back(e);
        atVertex[mesh->edges[e][1]].push_back(e);
    }
    auto otherEnd = [&](int e, int v) {
        return (mesh->edges[e][0] == v) ? mesh->edges[e][1] : mesh->edges[e][0];
    };
    // A vertex the chain may not run through: an end, a junction of three or
    // more branches, or a landing on dS, where which pair of edges continues
    // is not a question the network answers on its own.
    auto isChainNode = [&](int v) {
        return atVertex[v].size() != 2 || mesh->isBoundaryVertex[v];
    };

    std::vector<char> used(mesh->edges.size(), 0);
    auto walkFrom = [&](int v, int e) {
        FeatureChain fc;
        int cur = v, ce = e;
        fc.verts.push_back(v);
        while (ce >= 0 && !used[ce]) {
            used[ce] = 1;
            fc.edges.push_back(ce);
            const int nxt = otherEnd(ce, cur);
            fc.verts.push_back(nxt);
            if (isChainNode(nxt)) break;
            const int e2 = (atVertex[nxt][0] != ce) ? atVertex[nxt][0] : atVertex[nxt][1];
            cur = nxt;
            ce = e2;
        }
        if (fc.edges.empty()) return;
        fc.closed = fc.verts.front() == fc.verts.back();
        chains.push_back(std::move(fc));
    };
    // Open chains first, from their ends, so that a closed loop is only what
    // is left over after every branch has been taken.
    for (int v = 0; v < nV; ++v) {
        if (!isChainNode(v)) continue;
        for (const int e : atVertex[v]) if (!used[e]) walkFrom(v, e);
    }
    for (const int e : featureEdges) {
        if (e < 0 || e >= static_cast<int>(mesh->edges.size()) || used[e]) continue;
        if (mesh->isBoundaryEdge[e]) continue;
        walkFrom(mesh->edges[e][0], e);
    }
}

const std::vector<UMBER::FeatureChain>& UMBER::interfaceChains() const {
    if (!chainsBuilt) buildInterfaceChains();
    return chains;
}

std::vector<UMBER::FeatureCorner> UMBER::interfaceCorners() const {
    std::vector<FeatureCorner> result;
    if (uField.empty() || featureEdges.empty()) return result;

    // --- Each chain, each way round ----------------------------------------
    for (const FeatureChain &fc : interfaceChains()) {
        if (fc.edges.size() < 2) continue;
        const bool closed = fc.closed;

        for (int dir = 0; dir < 2; ++dir) {
            std::vector<int> ce = fc.edges, cv = fc.verts;
            if (dir == 1) {
                std::reverse(ce.begin(), ce.end());
                std::reverse(cv.begin(), cv.end());
            }
            double running = 0.0;
            long handedOut = 0;
            const size_t last = closed ? ce.size() : ce.size() - 1;
            for (size_t i = 0; i < last; ++i) {
                const int v = cv[i + 1];
                const int eIn = ce[i];
                const int eOut = ce[(i + 1) % ce.size()];
                const std::vector<int> fan = sectorFan(v, eIn, eOut);
                if (!fan.empty()) {
                    double omega = 0.0, theta = 0.0;
                    sectorMeasure(v, fan, omega, theta);
                    running += (theta + M_PI - omega) / M_PI_2;
                }
                const long want = std::lround(running);
                if (want != handedOut) {
                    result.push_back(FeatureCorner{v, eIn, eOut,
                                                   static_cast<int>(want - handedOut)});
                    handedOut = want;
                }
            }
        }
    }

    return result;
}

// ---------------------------------------------------------------------------
// sharpTurns()
// ---------------------------------------------------------------------------
std::vector<std::pair<int, double>> UMBER::sharpTurns(double thresholdDegrees) const {
    std::vector<std::pair<int, double>> result;
    if (uField.empty()) return result;

    for (int e = 0; e < static_cast<int>(mesh->edges.size()); ++e) {
        if (mesh->isBoundaryEdge[e]) continue;
        if (cutEdges.find(EdgeKey(mesh->edges[e][0], mesh->edges[e][1])) != cutEdges.end()) continue;

        const int fa = mesh->edgeTriangles[e][0];
        const int fb = mesh->edgeTriangles[e][1];
        if (fa < 0 || fb < 0) continue;

        const double deg = std::fabs(wrap_pi(angles[fb] - angles[fa])) * 180.0 / M_PI;
        if (deg > thresholdDegrees) result.emplace_back(e, deg);
    }

    std::sort(result.begin(), result.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });
    return result;
}

// ---------------------------------------------------------------------------
// crossFieldRepresentation()
// ---------------------------------------------------------------------------
Eigen::VectorXcd UMBER::crossFieldRepresentation() const {
    const int nT = static_cast<int>(uField.size());
    Eigen::VectorXcd u(nT);
    for (int t = 0; t < nT; ++t) {
        double t1, t2;
        tau(uField[t][0], uField[t][1], t1, t2);
        const double mag = std::sqrt(t1 * t1 + t2 * t2);
        u[t] = (mag > 1e-14) ? std::complex<double>(t1 / mag, t2 / mag)
                             : std::complex<double>(1.0, 0.0);
    }
    return u;
}
