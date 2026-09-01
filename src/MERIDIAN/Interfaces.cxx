#include "Interfaces.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <stdexcept>

namespace {

// The angle at v in triangle t.
inline double cornerAngle(const Mesh &m, int t, int v) {
    const Triangle &tri = m.triangles[t];
    int i0 = -1;
    for (int i = 0; i < 3; ++i) if (tri[i] == v) i0 = i;
    if (i0 < 0) return 0.0;
    const Point a = m.vertices[tri[(i0 + 1) % 3]] - m.vertices[v];
    const Point b = m.vertices[tri[(i0 + 2) % 3]] - m.vertices[v];
    return std::fabs(std::atan2(cross2(a, b), dotP(a, b)));
}

// Into [0, 2pi).
inline double wrap2pi(double a) {
    a = std::fmod(a, 2.0 * M_PI);
    if (a < 0.0) a += 2.0 * M_PI;
    return a;
}

} // namespace

// ---------------------------------------------------------------------------
Interfaces::Interfaces(std::shared_ptr<Mesh> m) : Interfaces(std::move(m), Options()) {}

Interfaces::Interfaces(std::shared_ptr<Mesh> m, const Options &opts)
    : mesh(std::move(m)), options(opts) {
    if (!mesh) throw std::runtime_error("Interfaces: null mesh");
    if (mesh->triangles.empty()) throw std::runtime_error("Interfaces: empty mesh");

    collectEdges();
    if (edgeList.empty()) return;

    findNodes();
    buildBranches();
    if (options.splitLoops) splitClosedLoops();
    // Before measureNodes(), so that every node knows which region each of its
    // sectors opens into from the start. balance() rebuilds both, and gets the
    // same answer; what this buys is that regionAt() is usable by Stage 1,
    // which runs before balance() does.
    buildRegions();
    measureNodes();

    report.nodes = static_cast<int>(nodeList.size());
    report.branches = static_cast<int>(branchList.size());
    for (const Node &n : nodeList) {
        switch (n.kind) {
            case NodeKind::Junction:  ++report.junctions; break;
            case NodeKind::Landing:   ++report.landings; break;
            case NodeKind::Kink:      ++report.kinks; break;
            case NodeKind::LoopSplit: ++report.loopSplits; break;
            case NodeKind::Dangling:  ++report.dangling; break;
        }
        if (!n.wellPosed) ++report.illPosedNodes;
        if (n.residual > report.worstSectorResidual) {
            report.worstSectorResidual = n.residual;
            report.worstSectorNode = static_cast<int>(&n - nodeList.data());
        }
    }
    for (const Branch &b : branchList) if (b.closed) ++report.closedLoops;

    for (const auto &p : prescription()) {
        if (p.second != 0) ++report.prescribedCones;
        report.prescribedIndexSum += p.second;
    }

    if (report.dangling > 0) {
        std::ostringstream oss;
        oss << report.dangling << " interface branch(es) simply stop inside the model. "
            << "An interface separates two materials and cannot end in the middle of "
            << "one, so this is a tagging error in the input, not a layout problem.";
        report.messages.push_back(oss.str());
    }
    if (report.illPosedNodes > 0) {
        std::ostringstream oss;
        oss << report.illPosedNodes << " of " << report.nodes
            << " interface node(s) have sectors that are not whole quarter turns; the "
            << "worst is off by " << (report.worstSectorResidual * 180.0 / M_PI)
            << " degrees. The layout has to absorb that turn, which it can, but no "
            << "vertex-based cross field is well-posed there.";
        report.messages.push_back(oss.str());
    }
}

// ---------------------------------------------------------------------------
// collectEdges()
//
// An interface edge is an interior edge whose two triangles carry different
// material tags. An edge of dS has one triangle and is not one: the boundary is
// E2's business, and an interface that runs along dS is the boundary.
// ---------------------------------------------------------------------------
void Interfaces::collectEdges() {
    const int nE = static_cast<int>(mesh->edges.size());
    const int nV = static_cast<int>(mesh->vertices.size());

    edgeIsInterface.assign(nE, 0);
    vertDegree.assign(nV, 0);
    vertNode.assign(nV, -1);
    vertBranch.assign(nV, -1);
    edgeBranch.assign(nE, -1);
    vertEdges.assign(nV, {});

    std::set<int> mats;
    if (mesh->triangleMatId.size() == mesh->triangles.size()) {
        for (int id : mesh->triangleMatId) mats.insert(id);
    } else {
        mats.insert(1);
    }
    report.materials = static_cast<int>(mats.size());
    if (report.materials < 2) return;

    for (int e = 0; e < nE; ++e) {
        const int f0 = mesh->edgeTriangles[e][0];
        const int f1 = mesh->edgeTriangles[e][1];
        if (f0 < 0 || f1 < 0) continue;
        if (mesh->triangleMatId[f0] == mesh->triangleMatId[f1]) continue;
        edgeIsInterface[e] = 1;
        edgeList.push_back(e);
        for (int k = 0; k < 2; ++k) {
            const int v = mesh->edges[e][k];
            ++vertDegree[v];
            vertEdges[v].push_back(e);
        }
    }
    report.interfaceEdges = static_cast<int>(edgeList.size());
    for (int v = 0; v < nV; ++v) if (vertDegree[v] > 0) ++report.interfaceVertices;
}

// ---------------------------------------------------------------------------
// findNodes()
//
// The four reasons a vertex of the interface graph is a node, in the order they
// are tested: it is on dS, so the interface lands there; three or more branches
// meet, so it is a junction; only one does, which cannot happen on a
// conformally tagged mesh and is reported rather than repaired; or two do and
// the interface turns through more than Options::kinkAngle between them.
//
// The turn is measured between the two incident *edges* and not between chords
// a few vertices away. At a genuine kink the two edges are the two faces
// meeting there and their directions are the answer; along a smooth arc the
// per-edge turn is the discretisation's own, an order of magnitude under the
// threshold at any usable resolution.
// ---------------------------------------------------------------------------
void Interfaces::findNodes() {
    const int nV = static_cast<int>(mesh->vertices.size());

    auto rayAngle = [&](int v, int e) {
        const int w = (mesh->edges[e][0] == v) ? mesh->edges[e][1] : mesh->edges[e][0];
        const Point d = mesh->vertices[w] - mesh->vertices[v];
        return std::atan2(d[1], d[0]);
    };

    for (int v = 0; v < nV; ++v) {
        if (vertDegree[v] == 0) continue;

        NodeKind kind;
        if (mesh->isBoundaryVertex[v]) {
            // A branch has to end where the model does, so this is a node
            // whatever the options say; prescribeLandings only decides whether
            // its index is handed to Stage 1 or left to the field.
            kind = NodeKind::Landing;
        } else if (vertDegree[v] >= 3) {
            kind = NodeKind::Junction;
        } else if (vertDegree[v] == 1) {
            kind = NodeKind::Dangling;
            ++report.nonManifoldVertices;
        } else {
            const double a0 = rayAngle(v, vertEdges[v][0]);
            const double a1 = rayAngle(v, vertEdges[v][1]);
            // Straight through means the two rays are opposite; the turn of the
            // interface is how far from opposite they are.
            const double turn = M_PI - std::fabs(wrap_pi(a1 - a0));
            if (std::fabs(turn) <= options.kinkAngle) continue;   // not a node
            kind = NodeKind::Kink;
        }

        Node n;
        n.vertex = v;
        n.kind = kind;
        n.onBoundary = mesh->isBoundaryVertex[v];
        vertNode[v] = static_cast<int>(nodeList.size());
        nodeList.push_back(n);
    }
}

// ---------------------------------------------------------------------------
// buildBranches()
//
// Walk out of every node along every incident interface edge until the next
// node, marking edges as used. What is left over when that is done is a closed
// loop carrying no node at all -- an inclusion -- and is emitted as one branch
// with node0 = node1 = -1 for splitClosedLoops() to cut.
// ---------------------------------------------------------------------------
void Interfaces::buildBranches() {
    branchList.clear();
    std::fill(edgeBranch.begin(), edgeBranch.end(), -1);
    std::fill(vertBranch.begin(), vertBranch.end(), -1);
    const int nE = static_cast<int>(mesh->edges.size());
    std::vector<char> used(nE, 0);

    auto other = [&](int e, int v) {
        return (mesh->edges[e][0] == v) ? mesh->edges[e][1] : mesh->edges[e][0];
    };

    auto finish = [&](Branch &b) {
        b.length = 0.0;
        for (size_t i = 0; i + 1 < b.verts.size(); ++i) {
            b.length += normP(mesh->vertices[b.verts[i + 1]] - mesh->vertices[b.verts[i]]);
        }
        b.turning = 0.0;
        for (size_t i = 1; i + 1 < b.verts.size(); ++i) {
            const Point d0 = mesh->vertices[b.verts[i]] - mesh->vertices[b.verts[i - 1]];
            const Point d1 = mesh->vertices[b.verts[i + 1]] - mesh->vertices[b.verts[i]];
            b.turning += std::atan2(cross2(d0, d1), dotP(d0, d1));
        }
        if (!b.edges.empty()) {
            const auto mm = materialsAcross(b.edges.front());
            b.matLeft = mm.first;
            b.matRight = mm.second;
        }
        const int id = static_cast<int>(branchList.size());
        for (size_t i = 0; i < b.verts.size(); ++i) {
            const int v = b.verts[i];
            if (vertNode[v] < 0) vertBranch[v] = id;
        }
        for (int e : b.edges) edgeBranch[e] = id;
        branchList.push_back(std::move(b));
    };

    // Node to node.
    for (size_t ni = 0; ni < nodeList.size(); ++ni) {
        const int start = nodeList[ni].vertex;
        for (int e0 : vertEdges[start]) {
            if (used[e0]) continue;
            Branch b;
            b.node0 = static_cast<int>(ni);
            b.verts.push_back(start);

            int e = e0, v = start;
            while (true) {
                used[e] = 1;
                const int w = other(e, v);
                b.edges.push_back(e);
                b.verts.push_back(w);
                if (vertNode[w] >= 0) { b.node1 = vertNode[w]; break; }
                int nxt = -1;
                for (int c : vertEdges[w]) if (c != e && !used[c]) nxt = c;
                if (nxt < 0) break;          // a loop that came back on itself
                e = nxt;
                v = w;
            }
            finish(b);
        }
    }

    // Whatever is left is a closed loop with no node on it.
    for (int e0 : edgeList) {
        if (used[e0]) continue;
        Branch b;
        b.closed = true;
        const int start = mesh->edges[e0][0];
        b.verts.push_back(start);

        int e = e0, v = start;
        while (true) {
            used[e] = 1;
            const int w = other(e, v);
            b.edges.push_back(e);
            b.verts.push_back(w);
            if (w == start) break;
            int nxt = -1;
            for (int c : vertEdges[w]) if (c != e && !used[c]) nxt = c;
            if (nxt < 0) break;
            e = nxt;
            v = w;
        }
        finish(b);
    }
}

// ---------------------------------------------------------------------------
// circularComponentOfLoop()
//
// Whether a closed interface loop is the rim of a disk, which is what
// Options::splitCircleLoops turns on. The judgement is not made here: the
// loop's two sides name their material components and
// Mesh::computeMaterialCircles has already fitted a circle to each
// component's own boundary and applied its two acceptance tests.
// ---------------------------------------------------------------------------
int Interfaces::circularComponentOfLoop(const Branch &b) const {
    if (b.edges.empty() || mesh->triangleComponent.empty()) return -1;
    int best = -1;
    for (int k = 0; k < 2; ++k) {
        const int t = mesh->edgeTriangles[b.edges.front()][k];
        if (t < 0 || t >= static_cast<int>(mesh->triangleComponent.size())) continue;
        const int c = mesh->triangleComponent[t];
        if (c < 0 || c >= static_cast<int>(mesh->materialComponents.size())) continue;
        if (!mesh->materialComponents[c].circle.isCircle) continue;
        // Both sides can fit a circle only when one of them is an annulus
        // around the other; the inclusion is the smaller of the two and is the
        // one this loop is the whole boundary of.
        if (best < 0 || mesh->materialComponents[c].triangles.size() <
                        mesh->materialComponents[best].triangles.size()) {
            best = c;
        }
    }
    return best;
}

// ---------------------------------------------------------------------------
// splitClosedLoops()
//
// A closed interface loop carries no node, and a closed curve through no
// singular point cannot be a union of layout edges: the face it bounds has no
// corners, so it is not a quadrilateral, so nothing downstream can mesh it. The
// loop has to be cut, and where it is cut is not arbitrary -- an inclusion's
// layout wants its corners where the interface's own direction turns through a
// quarter, which for an ellipse is the four ends of its axes and for a rounded
// square is its four flats.
//
// So the cut points are read off the *unwrapped* tangent angle along the loop,
// which runs monotonically through 2 pi on a convex loop: split at the four
// parameters where it has advanced by a quarter of its total. That is exactly
// "where the tangent points along each of four directions a quarter apart",
// without having to know which four directions the layout will choose, and it
// degrades gracefully on a non-convex loop, where the unwrapped angle still
// totals 2 pi but is not monotone and the crossings are taken in order.
//
// A disk is the one loop this is not asked of, because it is the one loop that
// already has corners before the stage runs: see Options::splitCircleLoops.
// ---------------------------------------------------------------------------
void Interfaces::splitClosedLoops() {
    bool any = false;
    for (const Branch &b : branchList) if (b.closed) { any = true; break; }
    if (!any || options.loopSplits < 1) return;

    for (const Branch &b : branchList) {
        if (!b.closed || b.verts.size() < static_cast<size_t>(options.loopSplits) + 1) continue;

        // A disk's rim is already crossed four times by the separatrices of
        // the cones its own rotational symmetry puts inside it, so it already
        // has the corners splitting was for and a node here is one the layout
        // has nothing to meet. See Options::splitCircleLoops.
        if (!options.splitCircleLoops && circularComponentOfLoop(b) >= 0) continue;

        // Unwrapped tangent angle at each vertex of the loop.
        const size_t n = b.verts.size() - 1;      // verts.back() == verts.front()
        std::vector<double> phi(n + 1, 0.0);
        double acc = 0.0;
        for (size_t i = 0; i < n; ++i) {
            const Point d0 = mesh->vertices[b.verts[(i + 1) % n]] - mesh->vertices[b.verts[i]];
            const Point d1 = mesh->vertices[b.verts[(i + 2) % n]] -
                             mesh->vertices[b.verts[(i + 1) % n]];
            acc += std::atan2(cross2(d0, d1), dotP(d0, d1));
            phi[i + 1] = acc;
        }
        const double total = acc;
        if (std::fabs(total) < M_PI) continue;    // not a turn we can quarter

        for (int k = 0; k < options.loopSplits; ++k) {
            const double target = total * (static_cast<double>(k) / options.loopSplits);
            size_t best = 0;
            double bd = std::numeric_limits<double>::infinity();
            for (size_t i = 0; i < n; ++i) {
                const double d = std::fabs(phi[i] - target);
                if (d < bd) { bd = d; best = i; }
            }
            const int v = b.verts[best];
            if (vertNode[v] >= 0) continue;
            Node nd;
            nd.vertex = v;
            nd.kind = NodeKind::LoopSplit;
            nd.onBoundary = mesh->isBoundaryVertex[v];
            vertNode[v] = static_cast<int>(nodeList.size());
            nodeList.push_back(nd);
        }
    }
    buildBranches();
}

// ---------------------------------------------------------------------------
// materialsAcross()
// ---------------------------------------------------------------------------
std::vector<int> Interfaces::emitterNodes() const {
    std::vector<int> out;
    for (const Node &n : nodeList) {
        if (n.vertex < 0) continue;
        int spare = 0;
        for (int q : n.quarters) spare += q - 1;
        if (spare > 0 || n.index != 0) out.push_back(n.vertex);
    }
    return out;
}

// ---------------------------------------------------------------------------
// regionAt()
//
// The region a vertex is inside, or -1 when it is on the boundary between two
// of them. A vertex on an interface belongs to both regions and to neither for
// the purpose of anything counted per region, so it gets -1 rather than an
// arbitrary one of the two.
// ---------------------------------------------------------------------------
int Interfaces::regionAt(int v) const {
    if (triRegion.empty() || v < 0 || v >= static_cast<int>(mesh->vertices.size())) return -1;
    const auto tris = mesh->vertexTriangles.trianglesForVertex(v);
    int r = -1;
    for (const int *t = tris.first; t != tris.second; ++t) {
        const int rt = triRegion[*t];
        if (rt < 0) return -1;
        if (r < 0) r = rt;
        else if (r != rt) return -1;
    }
    return r;
}

// ---------------------------------------------------------------------------
std::pair<int, int> Interfaces::materialsAcross(int edge) const {
    if (edge < 0 || edge >= static_cast<int>(mesh->edges.size())) return {0, 0};
    const int f0 = mesh->edgeTriangles[edge][0];
    const int f1 = mesh->edgeTriangles[edge][1];
    if (f0 < 0 || f1 < 0) return {0, 0};
    if (mesh->triangleMatId.size() != mesh->triangles.size()) return {0, 0};
    return {mesh->triangleMatId[f0], mesh->triangleMatId[f1]};
}

// ---------------------------------------------------------------------------
// measureNodes()   --  the sectors, and the index they fix
// ---------------------------------------------------------------------------
void Interfaces::measureNodes() {
    // The triangle immediately counter-clockwise of the ray (v, w) at v is the
    // one whose third vertex lies counter-clockwise of w; the other one is
    // clockwise of it. Read from the geometry rather than from the triangle's
    // winding, which the loader does not normalise.
    auto sideFaces = [&](int v, int e, int &ccwFace, int &cwFace) {
        ccwFace = cwFace = -1;
        const int w = (mesh->edges[e][0] == v) ? mesh->edges[e][1] : mesh->edges[e][0];
        const Point d = mesh->vertices[w] - mesh->vertices[v];
        for (int k = 0; k < 2; ++k) {
            const int t = mesh->edgeTriangles[e][k];
            if (t < 0) continue;
            int x = -1;
            for (int i = 0; i < 3; ++i) {
                const int c = mesh->triangles[t][i];
                if (c != v && c != w) x = c;
            }
            if (x < 0) continue;
            const Point dx = mesh->vertices[x] - mesh->vertices[v];
            (cross2(d, dx) > 0.0 ? ccwFace : cwFace) = t;
        }
    };

    for (Node &n : nodeList) {
        const int v = n.vertex;
        n.rays.clear();
        n.branches.clear();
        n.sector.clear();

        // One ray per incident branch end. A branch that leaves v and comes
        // back to it contributes two, which is right: they are two rays.
        std::vector<Ray> rays;
        auto addRay = [&](int branch, int neighbour) {
            Ray r;
            r.branch = branch;
            r.neighbour = neighbour;
            const Point d = mesh->vertices[neighbour] - mesh->vertices[v];
            r.dir = std::atan2(d[1], d[0]);
            for (int e : vertEdges[v]) {
                if (mesh->edges[e][0] == neighbour || mesh->edges[e][1] == neighbour) r.edge = e;
            }
            if (r.edge >= 0) sideFaces(v, r.edge, r.faceCCW, r.faceCW);
            rays.push_back(r);
        };

        for (size_t b = 0; b < branchList.size(); ++b) {
            const Branch &br = branchList[b];
            if (br.verts.size() < 2) continue;
            if (br.verts.front() == v) addRay(static_cast<int>(b), br.verts[1]);
            if (br.verts.back() == v) addRay(static_cast<int>(b), br.verts[br.verts.size() - 2]);
        }
        if (rays.empty()) continue;

        if (!n.onBoundary) {
            // Interior: the rays cut the whole disk, and the sectors -- taken
            // cyclically -- sum to 2 pi.
            std::sort(rays.begin(), rays.end(),
                      [](const Ray &a, const Ray &b) { return a.dir < b.dir; });
            n.rays = std::move(rays);
            const size_t m = n.rays.size();
            for (size_t i = 0; i < m; ++i) {
                n.sector.push_back(wrap2pi(n.rays[(i + 1) % m].dir - n.rays[i].dir));
            }
        } else {
            // On dS: the fan runs from one boundary edge to the other and the
            // sectors sum to the interior angle, not to 2 pi. The two boundary
            // edges are rays of the node in their own right -- the interface
            // has to meet dS at a whole number of quarter turns just as it has
            // to meet another interface at one -- so they bound the list.
            double omega = 0.0;
            const auto tris = mesh->vertexTriangles.trianglesForVertex(v);
            for (const int *t = tris.first; t != tris.second; ++t) {
                omega += cornerAngle(*mesh, *t, v);
            }

            int b0 = -1, b1 = -1;
            for (int e : mesh->boundaryEdges) {
                if (mesh->edges[e][0] != v && mesh->edges[e][1] != v) continue;
                (b0 < 0 ? b0 : b1) = e;
            }
            if (b0 < 0 || b1 < 0) continue;   // a pinch; leave it to the field

            auto rayOf = [&](int e) {
                const int w = (mesh->edges[e][0] == v) ? mesh->edges[e][1] : mesh->edges[e][0];
                const Point d = mesh->vertices[w] - mesh->vertices[v];
                return std::atan2(d[1], d[0]);
            };
            double a0 = rayOf(b0), a1 = rayOf(b1);

            // Which of the two boundary rays the fan starts at -- that is,
            // which way round from it the model lies. Comparing the two sweeps
            // wrap2pi(a1 - a0) and wrap2pi(a0 - a1) against the interior angle
            // is the obvious test and it is exactly wrong on the commonest case
            // there is: a straight run of boundary has an interior angle of pi
            // and both sweeps are pi, so the comparison is a coin toss and half
            // the landings in the corpus came out with their sectors read the
            // wrong way round. The triangle answers it outright. The single
            // triangle on b0 spans the corner at v between b0's ray and its own
            // third vertex, so if the fan runs counter-clockwise out of a0 that
            // third vertex sits at a0 + (that corner angle), and if it runs the
            // other way it sits at a0 - it.
            {
                const int t = (mesh->edgeTriangles[b0][0] >= 0) ? mesh->edgeTriangles[b0][0]
                                                                : mesh->edgeTriangles[b0][1];
                const int nb0 = (mesh->edges[b0][0] == v) ? mesh->edges[b0][1]
                                                          : mesh->edges[b0][0];
                int w = -1;
                for (int i = 0; i < 3; ++i) {
                    const int c = mesh->triangles[t][i];
                    if (c != v && c != nb0) w = c;
                }
                if (w >= 0) {
                    const Point d = mesh->vertices[w] - mesh->vertices[v];
                    const double aw = std::atan2(d[1], d[0]);
                    const double corner = cornerAngle(*mesh, t, v);
                    if (std::fabs(wrap2pi(aw - a0) - corner) >
                        std::fabs(wrap2pi(a0 - aw) - corner)) {
                        std::swap(a0, a1);
                        std::swap(b0, b1);
                    }
                }
            }

            std::sort(rays.begin(), rays.end(), [&](const Ray &x, const Ray &y) {
                return wrap2pi(x.dir - a0) < wrap2pi(y.dir - a0);
            });

            Ray r0, r1;
            r0.branch = r1.branch = -1;
            r0.edge = b0; r1.edge = b1;
            r0.neighbour = (mesh->edges[b0][0] == v) ? mesh->edges[b0][1] : mesh->edges[b0][0];
            r1.neighbour = (mesh->edges[b1][0] == v) ? mesh->edges[b1][1] : mesh->edges[b1][0];
            r0.dir = a0; r1.dir = a1;
            sideFaces(v, b0, r0.faceCCW, r0.faceCW);
            sideFaces(v, b1, r1.faceCCW, r1.faceCW);

            n.rays.push_back(r0);
            for (const Ray &r : rays) n.rays.push_back(r);
            n.rays.push_back(r1);

            double prev = 0.0;
            for (size_t i = 1; i < n.rays.size(); ++i) {
                const double s = (i + 1 == n.rays.size()) ? omega
                                                          : wrap2pi(n.rays[i].dir - a0);
                n.sector.push_back(std::max(0.0, s - prev));
                prev = s;
            }
        }

        for (const Ray &r : n.rays) if (r.branch >= 0) n.branches.push_back(r.branch);

        // Which region each sector is in: the face on the counter-clockwise
        // side of the ray that opens it lies in that sector, by construction.
        n.sectorRegion.assign(n.sector.size(), -1);
        if (!triRegion.empty()) {
            for (size_t k = 0; k < n.sector.size() && k < n.rays.size(); ++k) {
                const int f = n.rays[k].faceCCW;
                if (f >= 0 && f < static_cast<int>(triRegion.size())) {
                    n.sectorRegion[k] = triRegion[f];
                }
            }
        }
        quantiseNode(n);
    }
}

// ---------------------------------------------------------------------------
// quantiseNode()
//
//     q_s = max(1, round(sector_s / (pi/2)))          every sector is a quad
//     I(v) = 4 - sum q_s   (interior)                 cone angle 2pi - (pi/2)I
//     I(v) = 2 - sum q_s   (on dS)                    angle      pi - (pi/2)I
//
// The clamp to at least one is what makes this a layout condition rather than a
// rounding: two interface branches leaving a vertex 20 degrees apart still have
// a quadrilateral between them in the output, because the alternative is an
// element with two of its sides on the interface.
// ---------------------------------------------------------------------------
void Interfaces::quantiseNode(Node &n) const {
    n.quarters.clear();
    for (double s : n.sector) {
        int q = static_cast<int>(std::lround(s / M_PI_2));
        if (q < 1) q = 1;
        n.quarters.push_back(q);
    }

    // Anything balance() moved at this node, put back. The base reading is a
    // local one and balance() is the only thing that knows the region-wide
    // count, so its answer wins wherever it has one.
    auto it = sectorAdjust.find(n.vertex);
    if (it != sectorAdjust.end() && it->second.size() == n.quarters.size()) {
        n.rebalanced = 0;
        for (size_t k = 0; k < n.quarters.size(); ++k) {
            n.quarters[k] = std::max(1, n.quarters[k] + it->second[k]);
            n.rebalanced += std::abs(it->second[k]);
        }
    }

    int total = 0;
    n.residual = 0.0;
    for (size_t k = 0; k < n.quarters.size(); ++k) {
        total += n.quarters[k];
        n.residual = std::max(n.residual,
                              std::fabs(n.sector[k] - n.quarters[k] * M_PI_2));
    }
    n.index = (n.onBoundary ? 2 : 4) - total;
    n.wellPosed = n.residual <= options.wellPosedAngle;
}

// ---------------------------------------------------------------------------
// prescription()
// ---------------------------------------------------------------------------
std::vector<std::pair<int, int>> Interfaces::prescription() const {
    std::vector<std::pair<int, int>> out;
    for (const Node &n : nodeList) {
        if (n.sector.empty()) continue;
        switch (n.kind) {
            case NodeKind::Landing: if (!options.prescribeLandings) continue; break;
            case NodeKind::Kink:    if (!options.prescribeKinks) continue; break;
            case NodeKind::Dangling: continue;
            default: break;
        }
        out.push_back({n.vertex, n.index});
    }
    return out;
}

// ---------------------------------------------------------------------------
// buildRegions()
//
// One region per connected component of the triangles under an adjacency that
// stops at the interfaces. This is the object the balance below is a statement
// about: the piece of one material that has to come out of the pipeline as a
// whole number of quadrilateral patches.
// ---------------------------------------------------------------------------
void Interfaces::buildRegions() {
    const int nT = static_cast<int>(mesh->triangles.size());
    regionList.clear();
    triRegion.assign(nT, -1);

    std::vector<int> stack;
    for (int seed = 0; seed < nT; ++seed) {
        if (triRegion[seed] >= 0) continue;
        const int id = static_cast<int>(regionList.size());
        Region reg;
        reg.material = mesh->triangleMatId.empty() ? 1 : mesh->triangleMatId[seed];

        stack.assign(1, seed);
        triRegion[seed] = id;
        while (!stack.empty()) {
            const int t = stack.back();
            stack.pop_back();
            reg.triangles.push_back(t);
            const Triangle &tri = mesh->triangles[t];
            reg.area += std::fabs(0.5 * cross2(mesh->vertices[tri[1]] - mesh->vertices[tri[0]],
                                               mesh->vertices[tri[2]] - mesh->vertices[tri[0]]));
            for (int k = 0; k < 3; ++k) {
                const int e = mesh->triangleEdges[t][k];
                if (e < 0 || edgeIsInterface[e]) continue;
                const int nb = (mesh->edgeTriangles[e][0] == t) ? mesh->edgeTriangles[e][1]
                                                                : mesh->edgeTriangles[e][0];
                if (nb < 0 || triRegion[nb] >= 0) continue;
                triRegion[nb] = id;
                stack.push_back(nb);
            }
        }

        // chi = V - E + F over the region's own triangles. A vertex or an edge
        // on an interface is used by both regions and is counted once in each,
        // which is what makes each of them a surface with boundary in its own
        // right rather than two halves of one.
        std::unordered_set<int> vs, es;
        for (int t : reg.triangles) {
            for (int k = 0; k < 3; ++k) {
                vs.insert(mesh->triangles[t][k]);
                es.insert(mesh->triangleEdges[t][k]);
            }
        }
        reg.chi = static_cast<int>(vs.size()) - static_cast<int>(es.size()) +
                  static_cast<int>(reg.triangles.size());
        regionList.push_back(std::move(reg));
    }
    report.regions = static_cast<int>(regionList.size());

    for (size_t b = 0; b < branchList.size(); ++b) {
        if (branchList[b].edges.empty()) continue;
        const int e = branchList[b].edges.front();
        for (int k = 0; k < 2; ++k) {
            const int t = mesh->edgeTriangles[e][k];
            if (t < 0) continue;
            std::vector<int> &lst = regionList[triRegion[t]].branches;
            if (std::find(lst.begin(), lst.end(), static_cast<int>(b)) == lst.end()) {
                lst.push_back(static_cast<int>(b));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// measureRegions()
//
//     turning(R)  =  sum_{v inside R} I(v)  +  sum_{v on dR} (2 - q_v)
//
// in quarter turns, against 4 chi(R). The three kinds of boundary vertex are
// read differently and each for a reason:
//
//   * on dS and not on an interface: q = 2 - I(v). Stage 1 already decided what
//     the layout does at a boundary vertex and this is that decision, not a
//     second reading of the same corner.
//   * an interface node: q is the sector of this region, from quantiseNode().
//   * anywhere else on an interface: q = 2. The interface runs straight through
//     and the layout does not turn.
// ---------------------------------------------------------------------------
void Interfaces::measureRegions(const std::vector<int> &coneIndex) {
    const int nV = static_cast<int>(mesh->vertices.size());
    if (regionList.empty()) return;

    // Which sector of which node a (region, vertex) pair is, so that the walk
    // below is a lookup and not a search. A node can open more than one sector
    // into the same region -- a branch that leaves it and comes back does -- so
    // the quarters are summed rather than assigned.
    std::unordered_map<long long, int> sectorOf;
    for (const Node &n : nodeList) {
        for (size_t k = 0; k < n.quarters.size() && k < n.sectorRegion.size(); ++k) {
            if (n.sectorRegion[k] < 0) continue;
            const long long key = static_cast<long long>(n.sectorRegion[k]) * nV + n.vertex;
            sectorOf[key] += n.quarters[k];
        }
    }

    for (size_t r = 0; r < regionList.size(); ++r) {
        Region &reg = regionList[r];

        // Which of its vertices are on its boundary, and which are inside.
        std::unordered_set<int> onBoundary, inside;
        for (int t : reg.triangles) {
            for (int k = 0; k < 3; ++k) {
                const int e = mesh->triangleEdges[t][k];
                if (!(mesh->isBoundaryEdge[e] || edgeIsInterface[e])) continue;
                onBoundary.insert(mesh->edges[e][0]);
                onBoundary.insert(mesh->edges[e][1]);
            }
        }
        for (int t : reg.triangles) {
            for (int k = 0; k < 3; ++k) {
                const int v = mesh->triangles[t][k];
                if (!onBoundary.count(v)) inside.insert(v);
            }
        }

        int turning = 0;
        for (int v : inside) turning += (v < static_cast<int>(coneIndex.size())) ? coneIndex[v] : 0;
        for (int v : onBoundary) {
            int q;
            const long long key = static_cast<long long>(r) * nV + v;
            auto it = sectorOf.find(key);
            if (vertNode[v] >= 0 && it != sectorOf.end()) {
                q = it->second;
            } else if (mesh->isBoundaryVertex[v] && vertDegree[v] == 0) {
                q = 2 - ((v < static_cast<int>(coneIndex.size())) ? coneIndex[v] : 0);
            } else {
                q = 2;
            }
            turning += 2 - q;
        }

        reg.turning = turning;
        reg.deficit = 4 * reg.chi - turning;
    }
}

// ---------------------------------------------------------------------------
// insertCorner()
// ---------------------------------------------------------------------------
bool Interfaces::insertCorner(int branch, int region, const std::vector<int> &coneIndex) {
    Branch &br = branchList[branch];
    if (br.verts.size() < 3) return false;

    // Where the branch has turned half of what it turns in total. On a smooth
    // arc that is the natural place for the layout's corner; on a straight
    // interface every point has turned nothing and the test degrades to the
    // midpoint, which is as good an answer as there is.
    std::vector<double> acc(br.verts.size(), 0.0);
    double total = 0.0;
    for (size_t i = 1; i + 1 < br.verts.size(); ++i) {
        const Point d0 = mesh->vertices[br.verts[i]] - mesh->vertices[br.verts[i - 1]];
        const Point d1 = mesh->vertices[br.verts[i + 1]] - mesh->vertices[br.verts[i]];
        total += std::atan2(cross2(d0, d1), dotP(d0, d1));
        acc[i] = total;
    }
    size_t best = br.verts.size() / 2;
    if (std::fabs(total) > 1e-9) {
        double bd = std::numeric_limits<double>::infinity();
        for (size_t i = 1; i + 1 < br.verts.size(); ++i) {
            const double d = std::fabs(acc[i] - 0.5 * total);
            if (d < bd) { bd = d; best = i; }
        }
    }
    const int v = br.verts[best];
    if (vertNode[v] >= 0) return false;

    Node nd;
    nd.vertex = v;
    nd.kind = NodeKind::Balance;
    nd.onBoundary = mesh->isBoundaryVertex[v];
    vertNode[v] = static_cast<int>(nodeList.size());
    nodeList.push_back(nd);

    buildBranches();
    buildRegions();
    measureNodes();

    // The new node reads as two straight-through sectors, since the interface
    // is smooth there; the layout wants one quarter on the side that was short
    // and three on the side that was over.
    {
        Node &added = nodeList[vertNode[v]];
        if (added.quarters.size() != 2 || added.onBoundary) return false;
        const int want = (added.sectorRegion[0] == region) ? 0 : 1;
        std::vector<int> delta(2, 0);
        delta[want] = 1 - added.quarters[want];
        delta[1 - want] = 3 - added.quarters[1 - want];
        sectorAdjust[v] = delta;
        quantiseNode(added);
    }
    measureRegions(coneIndex);
    return true;
}

// ---------------------------------------------------------------------------
// balance()
// ---------------------------------------------------------------------------
void Interfaces::balance(const std::vector<int> &coneIndex) {
    if (!multiMaterial()) return;
    buildRegions();
    measureNodes();          // sectorRegion needs triRegion, which now exists
    measureRegions(coneIndex);

    // One hop of the flow below: a quarter turn handed from region `give` to
    // region `take` across something the two of them share.
    //
    // A quarter can only be handed between two regions that meet, and the two
    // regions whose counts are wrong need not meet -- in geom013 the ferritic
    // plate has one turn too many and the weld one too few, and the buttering
    // layer sits between them precisely so that they never touch. So the move
    // is a *path* through the region graph and not a swap between neighbours:
    // every region on the way takes a quarter and gives one, which leaves its
    // own count where it was.
    struct Hop {
        int node = -1;       // a quarter moved across an existing node
        int giveSector = -1, takeSector = -1;
        int branch = -1;     // or a new corner put on a branch
        double cost = 0.0;
    };

    // The cheapest way to hand a quarter from `give` to `take`, or a Hop with
    // both fields negative if there is none.
    auto findHop = [&](int give, int take) {
        Hop best;
        best.cost = std::numeric_limits<double>::infinity();
        for (size_t n = 0; n < nodeList.size(); ++n) {
            const Node &nd = nodeList[n];
            for (size_t k = 0; k < nd.quarters.size() && k < nd.sectorRegion.size(); ++k) {
                if (nd.sectorRegion[k] != give || nd.quarters[k] <= 1) continue;
                for (size_t j = 0; j < nd.quarters.size() && j < nd.sectorRegion.size(); ++j) {
                    if (j == k || nd.sectorRegion[j] != take) continue;
                    const double after =
                        std::fabs(nd.sector[k] - (nd.quarters[k] - 1) * M_PI_2) +
                        std::fabs(nd.sector[j] - (nd.quarters[j] + 1) * M_PI_2);
                    if (after < best.cost) {
                        best.cost = after;
                        best.node = static_cast<int>(n);
                        best.giveSector = static_cast<int>(k);
                        best.takeSector = static_cast<int>(j);
                        best.branch = -1;
                    }
                }
            }
        }
        if (best.node >= 0) return best;

        // Nothing to move, so make somewhere to move it: a new corner on a
        // branch the two share, one quarter on the giving side and three on the
        // taking side. Priced above every quarter move, and among branches the
        // one that turns most is preferred -- a corner on an arc that is
        // already bending is one a coordinate line can follow, a corner in the
        // middle of a straight interface is one the layout has to invent.
        for (size_t b = 0; b < branchList.size(); ++b) {
            const Branch &br = branchList[b];
            if (br.edges.empty() || br.verts.size() < 3) continue;
            int ra = -1, rb = -1;
            for (int k = 0; k < 2; ++k) {
                const int t = mesh->edgeTriangles[br.edges.front()][k];
                if (t < 0) continue;
                (ra < 0 ? ra : rb) = triRegion[t];
            }
            if (!((ra == give && rb == take) || (ra == take && rb == give))) continue;
            const double cost = 1e3 + 1.0 / (1.0 + std::fabs(br.turning));
            if (cost < best.cost) {
                best.cost = cost;
                best.node = -1;
                best.branch = static_cast<int>(b);
            }
        }
        return best;
    };

    if (options.balanceRegions) {
        for (int move = 0; move < options.maxBalanceMoves; ++move) {
            int rShort = -1;
            for (size_t r = 0; r < regionList.size(); ++r) {
                if (regionList[r].deficit <= 0) continue;
                if (rShort < 0 || regionList[r].deficit > regionList[rShort].deficit) {
                    rShort = static_cast<int>(r);
                }
            }
            if (rShort < 0) break;

            // Breadth-first for the nearest region that has turned too much,
            // over the hops that exist.
            std::vector<int> prev(regionList.size(), -2);
            std::vector<int> queue{rShort};
            prev[rShort] = -1;
            int target = -1;
            for (size_t qi = 0; qi < queue.size() && target < 0; ++qi) {
                const int a = queue[qi];
                for (size_t b = 0; b < regionList.size(); ++b) {
                    if (prev[b] != -2) continue;
                    if (findHop(a, static_cast<int>(b)).cost ==
                        std::numeric_limits<double>::infinity()) {
                        continue;
                    }
                    prev[b] = a;
                    if (regionList[b].deficit < 0) { target = static_cast<int>(b); break; }
                    queue.push_back(static_cast<int>(b));
                }
            }
            if (target < 0) break;

            std::vector<int> path;
            for (int r = target; r >= 0; r = prev[r]) path.push_back(r);
            std::reverse(path.begin(), path.end());

            bool applied = false;
            for (size_t i = 0; i + 1 < path.size(); ++i) {
                const Hop h = findHop(path[i], path[i + 1]);
                if (h.node >= 0) {
                    Node &nd = nodeList[h.node];
                    std::vector<int> &delta = sectorAdjust[nd.vertex];
                    if (delta.size() != nd.quarters.size()) delta.assign(nd.quarters.size(), 0);
                    delta[h.giveSector] -= 1;
                    delta[h.takeSector] += 1;
                    quantiseNode(nd);
                    ++report.quartersMoved;
                    applied = true;
                } else if (h.branch >= 0) {
                    if (insertCorner(h.branch, path[i], coneIndex)) {
                        ++report.cornersInserted;
                        applied = true;
                    }
                    // insertCorner rebuilt the branches and the nodes, so the
                    // rest of this path is stale; recompute from the top.
                    break;
                } else {
                    break;
                }
            }
            measureRegions(coneIndex);
            if (!applied) break;
        }
    }

    report.regionsBalanced = 0;
    report.worstRegionDeficit = 0;
    for (const Region &reg : regionList) {
        if (reg.deficit == 0) ++report.regionsBalanced;
        if (std::abs(reg.deficit) > std::abs(report.worstRegionDeficit)) {
            report.worstRegionDeficit = reg.deficit;
        }
    }
    if (report.regionsBalanced < report.regions) {
        std::ostringstream oss;
        oss << (report.regions - report.regionsBalanced) << " of " << report.regions
            << " material region(s) do not satisfy their own Gauss-Bonnet count; the "
            << "worst is " << report.worstRegionDeficit << " quarter turn(s) out. Those "
            << "regions cannot be a whole number of quadrilateral patches, and the "
            << "remedy is upstream in Stage 1, which chose where the cones on dS went.";
        report.messages.push_back(oss.str());
    }
    if (report.quartersMoved || report.cornersInserted) {
        std::ostringstream oss;
        oss << "Balanced the material regions by moving " << report.quartersMoved
            << " quarter turn(s) across existing nodes and putting "
            << report.cornersInserted << " new corner(s) on smooth interfaces.";
        report.messages.push_back(oss.str());
    }
}

// ---------------------------------------------------------------------------
bool Interfaces::writeOBJ(const std::string &filename) const {
    std::ofstream out(filename);
    if (!out) return false;
    out << "# MERIDIAN Stage 0b: the material interface network\n";
    out << "# " << branchList.size() << " branch(es), " << nodeList.size() << " node(s)\n";

    int base = 1;
    for (const Branch &b : branchList) {
        for (int v : b.verts) {
            out << "v " << mesh->vertices[v][0] << " " << mesh->vertices[v][1] << " 0\n";
        }
        out << "l";
        for (size_t i = 0; i < b.verts.size(); ++i) out << " " << (base + static_cast<int>(i));
        out << "\n";
        base += static_cast<int>(b.verts.size());
    }
    return true;
}
