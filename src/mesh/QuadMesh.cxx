// QuadMesh.cxx -- see QuadMesh.hxx. Topology, geometry and degree-of-freedom
// bookkeeping only; there is deliberately no smoothing here.
#include "QuadMesh.hxx"
#include "FeatureFrame.hxx"

#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <limits>
#include <map>
#include <iomanip>
#include <stdexcept>

#include "geom/Fitting.hxx"
#include "geom/Polyline.hxx"

namespace mesh {

// ---------------------------------------------------------------------------
// construction
// ---------------------------------------------------------------------------

QuadMesh::QuadMesh(const std::vector<Point> &verts, const std::vector<Quad> &cells)
    : QuadMesh(verts, cells, std::vector<int>{}, Options{}) {}

QuadMesh::QuadMesh(const std::vector<Point> &verts, const std::vector<Quad> &cells,
                   const std::vector<int> &matIds)
    : QuadMesh(verts, cells, matIds, Options{}) {}

QuadMesh::QuadMesh(const std::vector<Point> &verts, const std::vector<Quad> &cells,
                   const std::vector<int> &matIds, const Options &opts)
    : vertices(verts), quads(cells), options(opts) {
    if (matIds.size() == quads.size()) quadMatId = matIds;
    else quadMatId.assign(quads.size(), 1);
    buildTopology();
}

// Parse an OBJ index token "i", "i/j" or "i/j/k" into its (1-based) vertex
// index. The same rule Mesh uses, kept local so the two readers do not have to
// share a private helper.
static bool parseObjIndex(const std::string &tok, int &vertexIndexOut) {
    if (tok.empty()) return false;
    const std::size_t slash = tok.find('/');
    const std::string vi = (slash == std::string::npos) ? tok : tok.substr(0, slash);
    try {
        vertexIndexOut = std::stoi(vi);
        return true;
    } catch (...) {
        return false;
    }
}

static bool parseMaterialId(const std::string &name, int &matIdOut) {
    std::size_t end = name.size();
    while (end > 0 && std::isdigit(static_cast<unsigned char>(name[end - 1]))) --end;
    if (end == name.size()) return false;
    try {
        matIdOut = std::stoi(name.substr(end));
        return true;
    } catch (...) {
        return false;
    }
}

QuadMesh::QuadMesh(const std::string &filename) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("QuadMesh: cannot open " + filename);

    int currentMat = 1;
    std::string line;
    int lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        std::istringstream ss(line);
        std::string tag;
        if (!(ss >> tag)) continue;
        if (tag == "v") {
            double x = 0.0, y = 0.0;
            ss >> x >> y;               // a z, if present, is dropped: this is 2-D
            vertices.push_back(Point{x, y});
        } else if (tag == "usemtl") {
            std::string name;
            ss >> name;
            parseMaterialId(name, currentMat);
        } else if (tag == "f") {
            std::vector<int> idx;
            std::string tok;
            while (ss >> tok) {
                int vi = 0;
                if (parseObjIndex(tok, vi)) idx.push_back(vi);
            }
            // A face that is not a quad is refused rather than triangulated or
            // dropped: either would leave the caller holding a different
            // element count than the file has, which is worse than an error.
            if (idx.size() != 4) {
                throw std::runtime_error("QuadMesh: " + filename + ":" +
                                         std::to_string(lineNo) + " face with " +
                                         std::to_string(idx.size()) +
                                         " vertices; this reader takes quads only");
            }
            Quad q{};
            for (int k = 0; k < 4; ++k) {
                int vi = idx[k];
                if (vi < 0) vi = static_cast<int>(vertices.size()) + vi;  // relative
                else        vi -= 1;                                      // 1-based
                if (vi < 0 || vi >= static_cast<int>(vertices.size())) {
                    throw std::runtime_error("QuadMesh: " + filename + ":" +
                                             std::to_string(lineNo) +
                                             " vertex index out of range");
                }
                q[k] = vi;
            }
            quads.push_back(q);
            quadMatId.push_back(currentMat);
        }
    }
    buildTopology();
}

// ---------------------------------------------------------------------------
// topology
// ---------------------------------------------------------------------------

void QuadMesh::buildTopology() {
    if (quadMatId.size() != quads.size()) quadMatId.assign(quads.size(), 1);
    if (pinned.size() != vertices.size()) pinned.resize(vertices.size(), false);

    orientQuads();
    buildEdges();
    buildBoundary();
    buildVertexQuads();
    buildBoundaryLoops();
    markInterfaceEdges();
    classifyNodes();
    buildFeatureCurves();
    computeSlideTangents();
    computeQuality();
}

int QuadMesh::orientQuads() {
    int flipped = 0;
    for (std::size_t q = 0; q < quads.size(); ++q) {
        if (signedArea(static_cast<int>(q)) < 0.0) {
            std::swap(quads[q][1], quads[q][3]);
            ++flipped;
        }
    }
    return flipped;
}

void QuadMesh::buildEdges() {
    edges.clear();
    quadEdges.assign(quads.size(), std::array<int, 4>{-1, -1, -1, -1});
    edgeQuads.clear();
    nonManifoldEdges.clear();

    std::unordered_map<MeshEdgeKey, int, MeshEdgeKeyHash> edgeMap;
    edgeMap.reserve(quads.size() * 4);
    std::unordered_set<int> nonManifold;

    for (int q = 0; q < static_cast<int>(quads.size()); ++q) {
        for (int s = 0; s < 4; ++s) {
            const MeshEdgeKey key(quads[q][s], quads[q][(s + 1) % 4]);
            auto it = edgeMap.find(key);
            int e;
            if (it == edgeMap.end()) {
                e = static_cast<int>(edges.size());
                edges.push_back(Edge{key.a, key.b});
                edgeQuads.push_back(std::array<int, 2>{q, -1});
                edgeMap.emplace(key, e);
            } else {
                e = it->second;
                if (edgeQuads[e][1] == -1) edgeQuads[e][1] = q;
                // Third and later users are recorded and otherwise dropped.
                // Storing one would mean overwriting a neighbour that is
                // already correct, and every walk downstream -- the vertex fan,
                // the boundary loops -- assumes at most two.
                else nonManifold.insert(e);
            }
            quadEdges[q][s] = e;
        }
    }

    nonManifoldEdges.assign(nonManifold.begin(), nonManifold.end());
    std::sort(nonManifoldEdges.begin(), nonManifoldEdges.end());

    // Side s of quad q faces whichever quad on edge quadEdges[q][s] is not q.
    quadAdjacency.assign(quads.size(), std::array<int, 4>{-1, -1, -1, -1});
    for (int q = 0; q < static_cast<int>(quads.size()); ++q) {
        for (int s = 0; s < 4; ++s) {
            const int e = quadEdges[q][s];
            if (e < 0) continue;
            const int a = edgeQuads[e][0], b = edgeQuads[e][1];
            quadAdjacency[q][s] = (a == q) ? b : (b == q ? a : -1);
        }
    }
}

void QuadMesh::buildBoundary() {
    isBoundaryEdge.assign(edges.size(), false);
    boundaryEdges.clear();
    for (int e = 0; e < static_cast<int>(edges.size()); ++e) {
        if (edgeQuads[e][1] == -1) {
            isBoundaryEdge[e] = true;
            boundaryEdges.push_back(e);
        }
    }

    isBoundaryVertex.assign(vertices.size(), false);
    boundaryVertices.clear();
    for (int e : boundaryEdges) {
        isBoundaryVertex[edges[e][0]] = true;
        isBoundaryVertex[edges[e][1]] = true;
    }
    for (int v = 0; v < static_cast<int>(vertices.size()); ++v)
        if (isBoundaryVertex[v]) boundaryVertices.push_back(v);
}

void QuadMesh::buildVertexQuads() {
    const int nV = static_cast<int>(vertices.size());
    manifoldRing.assign(nV, true);
    vertexNeighbors.assign(nV, {});

    std::vector<std::vector<std::array<int, 2>>> inc(nV);  // (quad, corner)
    for (int q = 0; q < static_cast<int>(quads.size()); ++q)
        for (int c = 0; c < 4; ++c) {
            const int v = quads[q][c];
            if (v >= 0 && v < nV) inc[v].push_back({q, c});
        }

    vertexQuads.rowPtr.assign(nV + 1, 0);
    vertexQuads.colIdx.clear();
    vertexQuads.corner.clear();
    vertexQuads.colIdx.reserve(quads.size() * 4);
    vertexQuads.corner.reserve(quads.size() * 4);

    auto cornerOf = [&](int q, int v) {
        for (int c = 0; c < 4; ++c) if (quads[q][c] == v) return c;
        return -1;
    };

    for (int v = 0; v < nV; ++v) {
        vertexQuads.rowPtr[v] = static_cast<int>(vertexQuads.colIdx.size());
        if (inc[v].empty()) continue;

        // Walking the fan. At corner c of quad q the two sides meeting v are
        // side c (v -> the next corner) and side (c+3)%4 (the previous corner
        // -> v). Crossing side (c+3)%4 turns counter-clockwise about v and
        // crossing side c turns clockwise, so: rewind clockwise until there is
        // no clockwise neighbour (an open fan hits the boundary; a closed one
        // returns to where it started, and there any quad will do), then walk
        // counter-clockwise from there.
        int startQ = inc[v][0][0], startC = inc[v][0][1];
        {
            int q = startQ, c = startC;
            for (std::size_t guard = 0; guard <= inc[v].size(); ++guard) {
                const int nq = quadAdjacency[q][c];
                if (nq < 0 || nq == startQ) break;
                const int nc = cornerOf(nq, v);
                if (nc < 0) break;
                q = nq; c = nc;
            }
            startQ = q; startC = c;
        }

        std::vector<std::array<int, 2>> fan;
        {
            int q = startQ, c = startC;
            std::unordered_set<int> seen;
            while (seen.insert(q).second) {
                fan.push_back({q, c});
                const int nq = quadAdjacency[q][(c + 3) % 4];
                if (nq < 0) break;
                const int nc = cornerOf(nq, v);
                if (nc < 0) break;
                q = nq; c = nc;
            }
        }

        // A fan that did not reach every incident quad is not a single ring:
        // the vertex is a pinch point, or sits on a non-manifold edge. Fall
        // back to index order and say so -- the membership is still right, it
        // is only the ordering that has no meaning there.
        if (fan.size() != inc[v].size()) {
            manifoldRing[v] = false;
            fan = inc[v];
        }

        for (const std::array<int, 2> &qc : fan) {
            vertexQuads.colIdx.push_back(qc[0]);
            vertexQuads.corner.push_back(qc[1]);
        }

        // One-ring vertices in the same order: each quad contributes the corner
        // ahead of v and the one behind it, which on a manifold fan comes out
        // counter-clockwise with each neighbour appearing once. The diagonal
        // corner is not a neighbour here -- these are the vertices joined to v
        // by an element side, which is the stencil a Laplacian or a local TMOP
        // patch solve wants.
        std::vector<int> ring;
        for (const std::array<int, 2> &qc : fan) {
            const int q = qc[0], c = qc[1];
            for (int k : {(c + 1) % 4, (c + 3) % 4}) {
                const int n = quads[q][k];
                if (std::find(ring.begin(), ring.end(), n) == ring.end()) ring.push_back(n);
            }
        }
        vertexNeighbors[v] = ring;
    }
    vertexQuads.rowPtr[nV] = static_cast<int>(vertexQuads.colIdx.size());
}

void QuadMesh::buildBoundaryLoops() {
    boundaryLoops.clear();
    boundaryLoopOf.assign(vertices.size(), std::array<int, 2>{-1, -1});

    // Directed boundary sides taken from the quad that owns them, so the
    // interior is on the left of every one of them. The outer loop therefore
    // comes out counter-clockwise and a hole clockwise.
    std::unordered_map<int, std::vector<int>> outgoing;      // vertex -> next vertices
    std::vector<std::array<int, 2>> sides;
    for (int q = 0; q < static_cast<int>(quads.size()); ++q)
        for (int s = 0; s < 4; ++s)
            if (quadAdjacency[q][s] < 0) {
                const int a = quads[q][s], b = quads[q][(s + 1) % 4];
                outgoing[a].push_back(b);
                sides.push_back({a, b});
            }

    std::unordered_set<long long> used;
    auto key = [](int a, int b) {
        return (static_cast<long long>(a) << 32) ^ static_cast<unsigned>(b);
    };

    for (const std::array<int, 2> &seed : sides) {
        if (used.count(key(seed[0], seed[1]))) continue;

        std::vector<int> loop;
        int a = seed[0], b = seed[1];
        while (used.insert(key(a, b)).second) {
            loop.push_back(a);
            auto it = outgoing.find(b);
            if (it == outgoing.end()) break;
            // The next unwalked side out of b. A pinch vertex has more than
            // one and either continuation closes a valid loop, so the first
            // unused one is as good as any.
            int next = -1;
            for (int cand : it->second)
                if (!used.count(key(b, cand))) { next = cand; break; }
            if (next < 0) break;
            a = b; b = next;
        }

        if (loop.size() >= 3) {
            const int li = static_cast<int>(boundaryLoops.size());
            for (std::size_t k = 0; k < loop.size(); ++k)
                if (boundaryLoopOf[loop[k]][0] < 0)
                    boundaryLoopOf[loop[k]] = {li, static_cast<int>(k)};
            boundaryLoops.push_back(std::move(loop));
        }
    }
}

void QuadMesh::markInterfaceEdges() {
    isInterfaceEdge.assign(edges.size(), false);
    for (int e = 0; e < static_cast<int>(edges.size()); ++e) {
        const int a = edgeQuads[e][0], b = edgeQuads[e][1];
        if (a >= 0 && b >= 0 && quadMatId[a] != quadMatId[b]) isInterfaceEdge[e] = true;
    }

    std::unordered_set<MeshEdgeKey, MeshEdgeKeyHash> marked;
    for (const Edge &m : markedFeatures) marked.insert(MeshEdgeKey(m[0], m[1]));

    isFeatureEdge.assign(edges.size(), false);
    for (int e = 0; e < static_cast<int>(edges.size()); ++e) {
        if (isBoundaryEdge[e]) isFeatureEdge[e] = true;
        else if (options.interfacesAreFeatures && isInterfaceEdge[e]) isFeatureEdge[e] = true;
        else if (marked.count(MeshEdgeKey(edges[e][0], edges[e][1]))) isFeatureEdge[e] = true;
    }

    isFeatureVertex.assign(vertices.size(), false);
    for (int e = 0; e < static_cast<int>(edges.size()); ++e) {
        if (!isFeatureEdge[e]) continue;
        isFeatureVertex[edges[e][0]] = true;
        isFeatureVertex[edges[e][1]] = true;
    }
}

// ---------------------------------------------------------------------------
// degrees of freedom
// ---------------------------------------------------------------------------

// The feature-curve neighbours of every vertex: the other endpoint of each
// incident feature edge. A smooth point of a curve has exactly two; anything
// else -- an end, a junction of three curves, a boundary node where an
// interface lands -- has some other number and cannot slide anywhere.
static std::vector<std::vector<int>> featureNeighbors(const QuadMesh &m) {
    std::vector<std::vector<int>> fn(m.vertices.size());
    for (int e = 0; e < static_cast<int>(m.edges.size()); ++e) {
        if (!m.isFeatureEdge[e]) continue;
        fn[m.edges[e][0]].push_back(m.edges[e][1]);
        fn[m.edges[e][1]].push_back(m.edges[e][0]);
    }
    return fn;
}

void QuadMesh::classifyNodes() {
    const int nV = static_cast<int>(vertices.size());
    nodeType.assign(nV, NodeFree);
    featureNeighborOf.assign(nV, std::array<int, 2>{-1, -1});
    if (pinned.size() != vertices.size()) pinned.resize(nV, false);

    const std::vector<std::vector<int>> fn = featureNeighbors(*this);
    const double cosLimit = std::cos((180.0 - options.cornerAngle) * M_PI / 180.0);

    for (int v = 0; v < nV; ++v) {
        // Recorded whatever the node type turns out to be, so the feature graph
        // stays walkable through corners and pins.
        if (fn[v].size() == 2) featureNeighborOf[v] = {fn[v][0], fn[v][1]};

        if (pinned[v]) { nodeType[v] = NodeFixed; continue; }
        if (!isFeatureVertex[v]) { nodeType[v] = NodeFree; continue; }
        if (options.fixAllFeatureNodes) { nodeType[v] = NodeFixed; continue; }

        // Not exactly two feature edges: an endpoint or a junction, and there
        // is no single curve to slide along.
        if (fn[v].size() != 2) { nodeType[v] = NodeFixed; continue; }

        const Point a = normalizeP(vertices[fn[v][0]] - vertices[v]);
        const Point b = normalizeP(vertices[fn[v][1]] - vertices[v]);
        // Straight through means a . b = -1. The node is a corner when the two
        // directions open up past the tolerance, i.e. their cosine rises above
        // cos(180 - cornerAngle).
        nodeType[v] = (dotP(a, b) > cosLimit) ? NodeFixed : NodeSliding;
    }

    anchorFreeFeatureLoops();
}

// Every node of a closed feature curve with no corner anywhere on it -- a
// finely discretised inclusion rim, a round outer dS -- comes out Sliding, and
// each is individually held to its own tangent. The *loop*, though, is held by
// nothing: its nodes can all creep the same way around it and reparameterise
// the curve, drifting the discretisation without any single node ever leaving
// the polyline. Pinning one node per such loop costs one degree of freedom and
// removes the drift.
//
// This walks the feature graph rather than `boundaryLoops`, so an interior
// material-interface ring is caught exactly as an outer boundary is. It runs
// after the pins and the corner test, so a loop that already has a fixed node
// anywhere on it is left alone: the walk simply reaches that node and stops.
void QuadMesh::anchorFreeFeatureLoops() {
    const int nV = static_cast<int>(vertices.size());
    std::vector<bool> visited(nV, false);

    for (int seed = 0; seed < nV; ++seed) {
        if (visited[seed] || nodeType[seed] != NodeSliding) continue;

        // Walk one way out of the seed for as long as the run stays sliding.
        // It ends either on a non-sliding node -- the run is an open arc
        // between two fixed ends and is already anchored -- or back on the
        // seed, which is the case this pass exists for.
        std::vector<int> run{seed};
        visited[seed] = true;
        int prev = seed, cur = featureNeighborOf[seed][1];
        bool closed = false;
        while (cur >= 0) {
            if (cur == seed) { closed = true; break; }
            if (nodeType[cur] != NodeSliding) break;
            run.push_back(cur);
            visited[cur] = true;
            const std::array<int, 2> &f = featureNeighborOf[cur];
            const int next = (f[0] == prev) ? f[1] : f[0];
            prev = cur;
            cur = next;
        }

        if (closed) {
            // Lowest index rather than the seed itself: the choice has to not
            // depend on which node the outer loop happened to reach first, so
            // a rebuild anchors the same node twice running.
            int anchor = run[0];
            for (int v : run) anchor = std::min(anchor, v);
            nodeType[anchor] = NodeFixed;
            continue;
        }

        // Open run: mark the other half visited too, so the next seed does not
        // walk the same arc again from its middle.
        prev = seed;
        cur = featureNeighborOf[seed][0];
        while (cur >= 0 && cur != seed && nodeType[cur] == NodeSliding) {
            visited[cur] = true;
            const std::array<int, 2> &f = featureNeighborOf[cur];
            const int next = (f[0] == prev) ? f[1] : f[0];
            prev = cur;
            cur = next;
        }
    }
}

// ---------------------------------------------------------------------------
// feature curves
// ---------------------------------------------------------------------------

std::vector<std::vector<int>> QuadMesh::featureChains() const {
    const int nV = static_cast<int>(vertices.size());
    std::vector<std::vector<int>> chains;
    if (static_cast<int>(nodeType.size()) != nV) return chains;

    const std::vector<std::vector<int>> fn = featureNeighbors(*this);
    std::vector<bool> used(nV, false);   // sliding nodes already on a chain

    // Runs start where sliding stops, so every one of them is bounded at both
    // ends by a node that does not move. A closed feature loop has exactly one
    // such node -- anchorFreeFeatureLoops() saw to that -- and comes back as a
    // single run from that anchor round to itself.
    for (int s = 0; s < nV; ++s) {
        if (nodeType[s] == NodeSliding) continue;
        for (int w : fn[s]) {
            if (w < 0 || w >= nV) continue;
            if (nodeType[w] != NodeSliding || used[w]) continue;

            std::vector<int> chain{s};
            int prev = s, cur = w;
            while (cur >= 0 && nodeType[cur] == NodeSliding && !used[cur]) {
                chain.push_back(cur);
                used[cur] = true;
                // A sliding node has exactly two feature neighbours: that is
                // what classifyNodes() requires before it calls one sliding.
                const std::array<int, 2> &f = featureNeighborOf[cur];
                const int next = (f[0] == prev) ? f[1] : f[0];
                prev = cur;
                cur = next;
            }
            // The far end, which closes the run. It is missing only if the
            // feature graph is malformed, and such a run is dropped: a curve
            // with a loose end has no fixed geometry to be pinned to.
            if (cur >= 0 && nodeType[cur] != NodeSliding) {
                chain.push_back(cur);
                chains.push_back(std::move(chain));
            }
        }
    }
    return chains;
}

int QuadMesh::buildFeatureCurves() {
    featureCurves.clear();
    curveOf.assign(vertices.size(), -1);
    curveParam.assign(vertices.size(), 0.0);
    if (options.curveSource == Options::CurveChord) return 0;
    if (nodeType.size() != vertices.size()) return 0;

    for (const std::vector<int> &chain : featureChains()) {
        std::vector<Point> pts;
        pts.reserve(chain.size());
        for (int v : chain) pts.push_back(vertices[v]);

        // The polyline is built whatever the source asks for: it is the
        // fallback, it supplies the parameters, and its length is what the
        // deviation test below is measured against.
        const geom::Polyline<2> poly(pts);
        if (!(poly.length() > 0.0) || poly.curve().empty()) continue;

        FeatureCurve fc;
        fc.chain = chain;
        fc.closed = chain.front() == chain.back();
        fc.length = poly.length();
        fc.curve = poly.curve();
        std::vector<double> param = poly.parameters();

        if (options.curveSource == Options::CurveSpline && pts.size() >= 4) {
            // Interpolation rather than a least-squares fit, so that binding
            // the nodes moves none of them: the curve passes through every one
            // at the parameter it is given here, and only *between* the nodes
            // does it differ from the polyline -- which is the whole of what
            // is wanted, the smooth curve the discretisation came from.
            try {
                const std::vector<double> sparam =
                    geom::parameterize(pts, geom::Parameterization::ChordLength);
                const geom::BSplineCurve<2> spline = geom::interpolateCurve(pts, 3);
                // How far it bows from the polyline, measured where the two are
                // furthest apart: at the middle of a segment, never at a node.
                double bow = 0.0;
                for (std::size_t i = 0; i + 1 < sparam.size(); ++i) {
                    const double um = 0.5 * (sparam[i] + sparam[i + 1]);
                    const Point mid{0.5 * (pts[i][0] + pts[i + 1][0]),
                                    0.5 * (pts[i][1] + pts[i + 1][1])};
                    bow = std::max(bow, normP(spline.evaluate(um) - mid));
                }
                const double meanSeg = poly.length() / static_cast<double>(pts.size() - 1);
                if (bow <= options.curveMaxDeviation * meanSeg) {
                    fc.curve = spline;
                    fc.fitted = true;
                    fc.deviation = bow;
                    param = sparam;
                }
            } catch (const std::invalid_argument &) {
                // Two chain vertices at one parameter: the run has a repeated
                // point and cannot be interpolated. Its polyline still can.
            }
        }

        const int index = static_cast<int>(featureCurves.size());
        for (std::size_t i = 0; i < chain.size(); ++i) {
            const int v = chain[i];
            // Only the sliding nodes are bound. The two ends do not move, and
            // on a closed run they are one vertex with two parameters, which
            // there would be no way to hold.
            if (nodeType[v] != NodeSliding) continue;
            curveOf[v] = index;
            curveParam[v] = param[i];
        }
        featureCurves.push_back(std::move(fc));
    }
    return static_cast<int>(featureCurves.size());
}

double QuadMesh::clampParameter(int v, double u) const {
    if (!isOnCurve(v)) return u;
    const geom::BSplineCurve<2> &c = featureCurves[curveOf[v]].curve;
    return std::min(std::max(u, c.domainBegin()), c.domainEnd());
}

Point QuadMesh::curvePoint(int v, double u) const {
    if (!isOnCurve(v)) return vertices[v];
    return featureCurves[curveOf[v]].curve.evaluate(clampParameter(v, u));
}

Point QuadMesh::curveTangent(int v, double u, double *speed) const {
    if (speed) *speed = 0.0;
    if (!isOnCurve(v)) return Point{0.0, 0.0};
    const Point d = featureCurves[curveOf[v]].curve.derivativeAt(clampParameter(v, u), 1);
    const double s = normP(d);
    if (speed) *speed = s;
    if (!(s > 0.0)) return Point{0.0, 0.0};
    return Point{d[0] / s, d[1] / s};
}

void QuadMesh::setVertexParameter(int v, double u) {
    if (!isOnCurve(v)) return;
    const double uc = clampParameter(v, u);
    curveParam[v] = uc;
    vertices[v] = featureCurves[curveOf[v]].curve.evaluate(uc);
}

double QuadMesh::maxCurveDeviation() const {
    double worst = 0.0;
    for (int v = 0; v < static_cast<int>(curveOf.size()); ++v) {
        if (!isOnCurve(v)) continue;
        worst = std::max(worst, normP(vertices[v] - curvePoint(v, curveParam[v])));
    }
    return worst;
}

void QuadMesh::updateSlideTangent(int v) {
    if (v < 0 || v >= static_cast<int>(vertices.size())) return;
    if (slideTangent.size() != vertices.size())
        slideTangent.assign(vertices.size(), Point{0.0, 0.0});
    if (v >= static_cast<int>(nodeType.size()) ||
        v >= static_cast<int>(featureNeighborOf.size()) ||
        nodeType[v] != NodeSliding) {
        slideTangent[v] = Point{0.0, 0.0};
        return;
    }
    // On a curve the tangent is the curve's own, at the parameter the node
    // currently holds. It needs no refreshing against the neighbours -- the
    // curve does not move when they do -- but it does have to follow the node
    // along the curve, which is why it is taken here rather than once.
    if (isOnCurve(v)) {
        slideTangent[v] = curveTangent(v, curveParam[v]);
        if (normP(slideTangent[v]) > 0.0) return;
        // A cusp, or a curve of no length: fall through to the chord.
    }
    const std::array<int, 2> &fn = featureNeighborOf[v];
    if (fn[0] < 0 || fn[1] < 0) { slideTangent[v] = Point{0.0, 0.0}; return; }
    // The chord through the two feature neighbours. On dS and on a material
    // interface this is not an approximation of some truer curve: those are
    // never splined, so the polyline through the neighbours is the geometry,
    // and a step along this direction keeps the node exactly on the segment it
    // already sits on.
    slideTangent[v] = normalizeP(vertices[fn[1]] - vertices[fn[0]]);
}

void QuadMesh::computeSlideTangents() {
    slideTangent.assign(vertices.size(), Point{0.0, 0.0});
    if (nodeType.size() != vertices.size()) return;
    for (int v = 0; v < static_cast<int>(vertices.size()); ++v) updateSlideTangent(v);
}

Point QuadMesh::projectStep(int v, const Point &displacement, double *newParam) const {
    if (newParam && v >= 0 && v < static_cast<int>(curveParam.size())) *newParam = curveParam[v];
    if (v < 0 || v >= static_cast<int>(nodeType.size())) return Point{0.0, 0.0};
    if (nodeType[v] == NodeFixed) return Point{0.0, 0.0};
    if (nodeType[v] != NodeSliding) return displacement;
    if (v >= static_cast<int>(slideTangent.size())) return Point{0.0, 0.0};

    const Point &t = slideTangent[v];
    if (isOnCurve(v)) {
        // The tangential part of the step, read as an arc length, divided by
        // the speed the parameter runs at to become a parameter step. First
        // order in du -- the curve's own curvature is not corrected for -- but
        // the node still lands exactly on the curve, and a Newton relaxation
        // takes the length it did not travel on its next sweep.
        double speed = 0.0;
        const Point tc = curveTangent(v, curveParam[v], &speed);
        if (!(speed > 0.0)) return Point{0.0, 0.0};
        const double u = clampParameter(v, curveParam[v] + dotP(tc, displacement) / speed);
        if (newParam) *newParam = u;
        return curvePoint(v, u) - vertices[v];
    }
    return t * dotP(t, displacement);
}

bool QuadMesh::markFeatureEdge(int edge) {
    if (edge < 0 || edge >= static_cast<int>(edges.size())) return false;
    markedFeatures.push_back(edges[edge]);
    isFeatureEdge[edge] = true;
    isFeatureVertex[edges[edge][0]] = true;
    isFeatureVertex[edges[edge][1]] = true;
    return true;
}

bool QuadMesh::markFeatureEdge(int va, int vb) {
    for (int e = 0; e < static_cast<int>(edges.size()); ++e)
        if ((edges[e][0] == va && edges[e][1] == vb) ||
            (edges[e][0] == vb && edges[e][1] == va))
            return markFeatureEdge(e);
    return false;
}

void QuadMesh::pinVertex(int v) {
    if (v < 0 || v >= static_cast<int>(vertices.size())) return;
    if (pinned.size() != vertices.size()) pinned.resize(vertices.size(), false);
    pinned[v] = true;
    if (nodeType.size() == vertices.size()) nodeType[v] = NodeFixed;
}

void QuadMesh::unpinAll() {
    pinned.assign(vertices.size(), false);
    classifyNodes();
    buildFeatureCurves();
    computeSlideTangents();
}

// ---------------------------------------------------------------------------
// targets
// ---------------------------------------------------------------------------

void QuadMesh::setUniformTargets(double h) {
    if (!(h > 0.0)) {
        computeQuality();
        h = quality.meanEdge;
    }
    if (!(h > 0.0)) h = 1.0;
    targetJacobian.assign(quads.size(), Jacobian2{h, 0.0, 0.0, h});
}

void QuadMesh::setTargetsFromCurrentShape() {
    targetJacobian.assign(quads.size(), Jacobian2{1.0, 0.0, 0.0, 1.0});
    for (int q = 0; q < static_cast<int>(quads.size()); ++q)
        targetJacobian[q] = jacobianAt(q, 0.5, 0.5);
}

void QuadMesh::setSizeTargets(const std::vector<double> &hPerQuad) {
    targetJacobian.assign(quads.size(), Jacobian2{1.0, 0.0, 0.0, 1.0});
    for (int q = 0; q < static_cast<int>(quads.size()) && q < static_cast<int>(hPerQuad.size()); ++q) {
        const double h = hPerQuad[q] > 0.0 ? hPerQuad[q] : 1.0;
        targetJacobian[q] = Jacobian2{h, 0.0, 0.0, h};
    }
}

// ---------------------------------------------------------------------------
// geometry
// ---------------------------------------------------------------------------

double QuadMesh::shapeFn(int c, double xi, double eta) {
    switch (c) {
        case 0: return (1.0 - xi) * (1.0 - eta);
        case 1: return xi * (1.0 - eta);
        case 2: return xi * eta;
        case 3: return (1.0 - xi) * eta;
        default: return 0.0;
    }
}

Point QuadMesh::shapeGrad(int c, double xi, double eta) {
    switch (c) {
        case 0: return Point{-(1.0 - eta), -(1.0 - xi)};
        case 1: return Point{ (1.0 - eta), -xi};
        case 2: return Point{ eta,          xi};
        case 3: return Point{-eta,         (1.0 - xi)};
        default: return Point{0.0, 0.0};
    }
}

Jacobian2 QuadMesh::jacobianAt(int q, double xi, double eta) const {
    Jacobian2 J{0.0, 0.0, 0.0, 0.0};
    if (q < 0 || q >= static_cast<int>(quads.size())) return J;
    for (int c = 0; c < 4; ++c) {
        const Point &x = vertices[quads[q][c]];
        const Point g = shapeGrad(c, xi, eta);
        J[0] += x[0] * g[0];  // dx/dxi
        J[1] += x[0] * g[1];  // dx/deta
        J[2] += x[1] * g[0];  // dy/dxi
        J[3] += x[1] * g[1];  // dy/deta
    }
    return J;
}

Jacobian2 QuadMesh::cornerJacobian(int q, int c) const {
    static const double XI[4]  = {0.0, 1.0, 1.0, 0.0};
    static const double ETA[4] = {0.0, 0.0, 1.0, 1.0};
    if (c < 0 || c > 3) return Jacobian2{0.0, 0.0, 0.0, 0.0};
    return jacobianAt(q, XI[c], ETA[c]);
}

double QuadMesh::signedArea(int q) const {
    if (q < 0 || q >= static_cast<int>(quads.size())) return 0.0;
    double a = 0.0;
    for (int k = 0; k < 4; ++k) {
        const Point &p = vertices[quads[q][k]];
        const Point &n = vertices[quads[q][(k + 1) % 4]];
        a += cross2(p, n);
    }
    return 0.5 * a;
}

Point QuadMesh::centroid(int q) const {
    if (q < 0 || q >= static_cast<int>(quads.size())) return Point{0.0, 0.0};
    Point c{0.0, 0.0};
    for (int k = 0; k < 4; ++k) c = c + vertices[quads[q][k]];
    return c / 4.0;
}

double QuadMesh::scaledJacobian(int q, int c) const {
    const Jacobian2 J = cornerJacobian(q, c);
    const double la = std::sqrt(J[0] * J[0] + J[2] * J[2]);
    const double lb = std::sqrt(J[1] * J[1] + J[3] * J[3]);
    if (la <= 0.0 || lb <= 0.0) return 0.0;
    return det2(J) / (la * lb);
}

double QuadMesh::minScaledJacobian(int q) const {
    double worst = std::numeric_limits<double>::max();
    for (int c = 0; c < 4; ++c) worst = std::min(worst, scaledJacobian(q, c));
    return worst;
}

void QuadMesh::elementSpans(int q, double &sXi, double &sEta) const {
    sXi = sEta = 0.0;
    if (q < 0 || q >= static_cast<int>(quads.size())) return;
    const Point &p0 = vertices[quads[q][0]], &p1 = vertices[quads[q][1]];
    const Point &p2 = vertices[quads[q][2]], &p3 = vertices[quads[q][3]];
    // The two mid-side spans in each direction, averaged: the size of the
    // element along xi and along eta, which is what a size target compares to.
    sXi  = 0.5 * (normP(p1 - p0) + normP(p2 - p3));
    sEta = 0.5 * (normP(p3 - p0) + normP(p2 - p1));
}

int QuadMesh::findQuadContainingPoint(const Point &p) const {
    auto inTriangle = [&](const Point &a, const Point &b, const Point &c) {
        const double d0 = cross2(b - a, p - a);
        const double d1 = cross2(c - b, p - b);
        const double d2 = cross2(a - c, p - c);
        const bool neg = (d0 < 0) || (d1 < 0) || (d2 < 0);
        const bool pos = (d0 > 0) || (d1 > 0) || (d2 > 0);
        return !(neg && pos);
    };
    for (int q = 0; q < static_cast<int>(quads.size()); ++q) {
        const Point &a = vertices[quads[q][0]], &b = vertices[quads[q][1]];
        const Point &c = vertices[quads[q][2]], &d = vertices[quads[q][3]];
        if (inTriangle(a, b, c) || inTriangle(a, c, d)) return q;
    }
    return -1;
}

double QuadMesh::alignmentQuality(const FeatureFrame &frame) const {
    // The same cell as a Coons cell of BlockDecomposition::alignmentQuality():
    // centre at the mean of the corners, directions the mean of opposite sides.
    FeatureFrame::Alignment grade(frame);
    for (int q = 0; q < static_cast<int>(quads.size()); ++q) {
        const Point &p0 = vertices[quads[q][0]], &p1 = vertices[quads[q][1]];
        const Point &p2 = vertices[quads[q][2]], &p3 = vertices[quads[q][3]];
        grade.add((p0 + p1 + p2 + p3) * 0.25, ((p1 - p0) + (p2 - p3)) * 0.5,
                  ((p3 - p0) + (p2 - p1)) * 0.5, std::fabs(signedArea(q)));
    }
    return grade.quality();
}

void QuadMesh::computeQuality() {
    Quality Q;
    Q.vertices = static_cast<int>(vertices.size());
    Q.quads = static_cast<int>(quads.size());
    Q.nonManifoldEdges = static_cast<int>(nonManifoldEdges.size());
    Q.boundaryEdges = static_cast<int>(boundaryEdges.size());
    Q.boundaryLoops = static_cast<int>(boundaryLoops.size());

    for (std::uint8_t t : nodeType) {
        if (t == NodeFree) ++Q.freeNodes;
        else if (t == NodeSliding) ++Q.slidingNodes;
        else ++Q.fixedNodes;
    }

    if (quads.empty()) { quality = Q; return; }

    Q.minScaledJacobian = std::numeric_limits<double>::max();
    Q.minArea = std::numeric_limits<double>::max();
    Q.allCounterClockwise = true;
    double sjSum = 0.0;

    for (int q = 0; q < static_cast<int>(quads.size()); ++q) {
        double worst = std::numeric_limits<double>::max();
        bool reflex = false;
        for (int c = 0; c < 4; ++c) {
            const double sj = scaledJacobian(q, c);
            worst = std::min(worst, sj);
            if (sj <= 0.0) reflex = true;
        }
        sjSum += worst;
        Q.minScaledJacobian = std::min(Q.minScaledJacobian, worst);
        if (worst <= 0.0) ++Q.invertedQuads;

        const double a = signedArea(q);
        if (a <= 0.0) Q.allCounterClockwise = false;
        else if (reflex) ++Q.nonConvexQuads;
        Q.minArea = std::min(Q.minArea, a);
        Q.maxArea = std::max(Q.maxArea, a);
        Q.totalArea += a;

        double sXi = 0.0, sEta = 0.0;
        elementSpans(q, sXi, sEta);
        const double lo = std::min(sXi, sEta), hi = std::max(sXi, sEta);
        if (lo > 0.0) Q.worstAspect = std::max(Q.worstAspect, hi / lo);
    }
    Q.meanScaledJacobian = sjSum / static_cast<double>(quads.size());

    Q.minEdge = std::numeric_limits<double>::max();
    double edgeSum = 0.0;
    for (const Edge &e : edges) {
        const double L = normP(vertices[e[1]] - vertices[e[0]]);
        Q.minEdge = std::min(Q.minEdge, L);
        Q.maxEdge = std::max(Q.maxEdge, L);
        edgeSum += L;
    }
    if (!edges.empty()) Q.meanEdge = edgeSum / static_cast<double>(edges.size());
    else Q.minEdge = 0.0;

    quality = Q;
}

// ---------------------------------------------------------------------------
// conversion and I/O
// ---------------------------------------------------------------------------

Mesh QuadMesh::toTriangleMesh() const {
    std::vector<Triangle> tris;
    std::vector<int> matIds;
    tris.reserve(quads.size() * 2);
    matIds.reserve(quads.size() * 2);
    for (int q = 0; q < static_cast<int>(quads.size()); ++q) {
        const int a = quads[q][0], b = quads[q][1], c = quads[q][2], d = quads[q][3];
        // Split on the shorter diagonal: on a non-convex quad it is the one
        // that stays inside, and on a convex one it is the one that gives the
        // better-shaped pair.
        const double ac = normP(vertices[c] - vertices[a]);
        const double bd = normP(vertices[d] - vertices[b]);
        if (ac <= bd) {
            tris.push_back(Triangle{a, b, c});
            tris.push_back(Triangle{a, c, d});
        } else {
            tris.push_back(Triangle{a, b, d});
            tris.push_back(Triangle{b, c, d});
        }
        const int m = q < static_cast<int>(quadMatId.size()) ? quadMatId[q] : 1;
        matIds.push_back(m);
        matIds.push_back(m);
    }
    return Mesh(vertices, tris, matIds);
}

bool QuadMesh::writeOBJ(const std::string &filename) const {
    std::ofstream out(filename);
    if (!out) return false;
    out << "# quad mesh: " << vertices.size() << " vertices, "
        << quads.size() << " quads\n";
    for (const Point &p : vertices)
        out << "v " << p[0] << " " << p[1] << " 0\n";

    // Faces grouped by material so `usemtl` reads back the ids this mesh
    // carries, the same convention Mesh's own reader expects.
    std::vector<int> order(quads.size());
    for (std::size_t k = 0; k < order.size(); ++k) order[k] = static_cast<int>(k);
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
        const int ma = a < static_cast<int>(quadMatId.size()) ? quadMatId[a] : 1;
        const int mb = b < static_cast<int>(quadMatId.size()) ? quadMatId[b] : 1;
        return ma < mb;
    });

    int current = std::numeric_limits<int>::min();
    for (int q : order) {
        const int m = q < static_cast<int>(quadMatId.size()) ? quadMatId[q] : 1;
        if (m != current) { out << "usemtl mat" << m << "\n"; current = m; }
        out << "f " << quads[q][0] + 1 << " " << quads[q][1] + 1 << " "
            << quads[q][2] + 1 << " " << quads[q][3] + 1 << "\n";
    }
    return true;
}

bool QuadMesh::writeVTU(const std::string &filename) const {
    std::ofstream out(filename);
    if (!out) return false;
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << vertices.size()
        << "\" NumberOfCells=\"" << quads.size() << "\">\n";

    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const Point &p : vertices) out << "          " << p[0] << " " << p[1] << " 0.0\n";
    out << "        </DataArray>\n      </Points>\n";

    out << "      <Cells>\n";
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const Quad &q : quads)
        out << "          " << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << "\n";
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (std::size_t k = 0; k < quads.size(); ++k) out << "          " << 4 * (k + 1) << "\n";
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::size_t k = 0; k < quads.size(); ++k) out << "          9\n";  // VTK_QUAD
    out << "        </DataArray>\n      </Cells>\n";

    out << "      <CellData Scalars=\"scaledJacobian\">\n";
    out << "        <DataArray type=\"Float64\" Name=\"scaledJacobian\" format=\"ascii\">\n";
    for (int q = 0; q < static_cast<int>(quads.size()); ++q)
        out << "          " << minScaledJacobian(q) << "\n";
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"material\" format=\"ascii\">\n";
    for (std::size_t k = 0; k < quads.size(); ++k)
        out << "          " << (k < quadMatId.size() ? quadMatId[k] : 0) << "\n";
    out << "        </DataArray>\n      </CellData>\n";

    // The node types go out too: a smoother's degrees of freedom are the first
    // thing to look at when it does not move what you expected it to.
    out << "      <PointData Scalars=\"nodeType\">\n";
    out << "        <DataArray type=\"Int32\" Name=\"nodeType\" format=\"ascii\">\n";
    for (std::size_t v = 0; v < vertices.size(); ++v)
        out << "          " << (v < nodeType.size() ? static_cast<int>(nodeType[v]) : 0) << "\n";
    out << "        </DataArray>\n      </PointData>\n";

    out << "    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";
    return true;
}

// ---------------------------------------------------------------------------
// MFEM
// ---------------------------------------------------------------------------

std::vector<std::pair<int, int>> QuadMesh::materialCounts() const {
    std::map<int, int> tally;
    for (std::size_t q = 0; q < quads.size(); ++q)
        ++tally[q < quadMatId.size() ? quadMatId[q] : 1];
    return std::vector<std::pair<int, int>>(tally.begin(), tally.end());
}

bool QuadMesh::writeMFEM(const std::string &filename) const {
    return writeMFEM(filename, MFEMOptions(), nullptr);
}

bool QuadMesh::writeMFEM(const std::string &filename, const MFEMOptions &mopts,
                         MFEMReport *report) const {
    std::ofstream out(filename);
    if (!out) return false;

    MFEMReport rep;

    // -- attributes --------------------------------------------------------
    // MFEM wants them strictly positive. Shift by one constant when they are
    // not, so the partition and the spacing between ids both survive: a caller
    // whose materials are 0 and 1 gets 1 and 2, still two distinct regions in
    // the same order.
    int shift = 0;
    int lowest = std::numeric_limits<int>::max();
    for (std::size_t q = 0; q < quads.size(); ++q)
        lowest = std::min(lowest, q < quadMatId.size() ? quadMatId[q] : 1);
    if (!quads.empty() && lowest < 1) {
        shift = 1 - lowest;
        rep.attributesRemapped = true;
    }
    auto attrOf = [&](int q) {
        return (q < static_cast<int>(quadMatId.size()) ? quadMatId[q] : 1) + shift;
    };

    // -- vertices actually used -------------------------------------------
    // A vertex no element references would be written and then never
    // addressed; MFEM accepts that but every consumer downstream has to carry
    // the gap, so drop them and renumber.
    std::vector<int> newIndex(vertices.size(), -1);
    std::vector<int> used;
    used.reserve(vertices.size());
    for (const Quad &q : quads)
        for (int c = 0; c < 4; ++c) {
            const int v = q[c];
            if (v >= 0 && v < static_cast<int>(vertices.size()) && newIndex[v] < 0) {
                newIndex[v] = static_cast<int>(used.size());
                used.push_back(v);
            }
        }
    rep.vertices = static_cast<int>(used.size());
    rep.unusedVertices = static_cast<int>(vertices.size()) - rep.vertices;

    // -- boundary segments -------------------------------------------------
    // Taken from the owning quad's own side, so each runs with the interior on
    // its left, the orientation MFEM expects of a boundary element. Only dS:
    // an interior material interface is the boundary between two attributes in
    // MFEM's model and is found from them, so emitting it here would be wrong,
    // not merely redundant.
    struct BdrSeg { int a, b, attr; };
    std::vector<BdrSeg> bdr;
    if (mopts.writeBoundary) {
        for (int q = 0; q < static_cast<int>(quads.size()); ++q)
            for (int s = 0; s < 4; ++s) {
                if (quadAdjacency[q][s] >= 0) continue;
                const int a = quads[q][s], b = quads[q][(s + 1) % 4];
                int attr = mopts.boundaryAttribute;
                if (mopts.boundaryByAxis) {
                    const Point d = vertices[b] - vertices[a];
                    if (std::fabs(d[0]) <= mopts.axisTolerance * normP(d))      attr = 1;
                    else if (std::fabs(d[1]) <= mopts.axisTolerance * normP(d)) attr = 2;
                }
                bdr.push_back(BdrSeg{a, b, attr});
            }
    }
    rep.boundaryElements = static_cast<int>(bdr.size());
    rep.elements = static_cast<int>(quads.size());

    // -- census ------------------------------------------------------------
    {
        std::map<int, int> tally;
        for (int q = 0; q < static_cast<int>(quads.size()); ++q) ++tally[attrOf(q)];
        rep.attributeCounts.assign(tally.begin(), tally.end());
        rep.attributes = static_cast<int>(tally.size());
    }

    // -- the file ----------------------------------------------------------
    out << "MFEM mesh v1.0\n\n";
    out << "#\n# Written by CrossGen (mesh::QuadMesh::writeMFEM).\n";
    out << "# MFEM geometry types: SEGMENT = 1, SQUARE = 3.\n";
    out << "# The element attribute is the material id.\n";
    if (mopts.annotate) {
        for (const std::pair<int, int> &kv : rep.attributeCounts)
            out << "#   material " << kv.first << ": " << kv.second << " element(s)\n";
        if (rep.attributesRemapped)
            out << "#   (ids shifted by " << shift << " so the smallest is 1, as MFEM"
                << " requires positive attributes)\n";
    }
    out << "#\n\n";

    out << "dimension\n2\n\n";

    out << "elements\n" << quads.size() << "\n";
    for (int q = 0; q < static_cast<int>(quads.size()); ++q) {
        out << attrOf(q) << " 3";
        for (int c = 0; c < 4; ++c) out << " " << newIndex[quads[q][c]];
        out << "\n";
    }
    out << "\n";

    out << "boundary\n" << bdr.size() << "\n";
    for (const BdrSeg &s : bdr)
        out << s.attr << " 1 " << newIndex[s.a] << " " << newIndex[s.b] << "\n";
    out << "\n";

    out << "vertices\n" << used.size() << "\n2\n";
    out << std::setprecision(mopts.precision);
    for (int v : used) out << vertices[v][0] << " " << vertices[v][1] << "\n";

    if (report) *report = rep;
    return static_cast<bool>(out);
}

}  // namespace mesh
