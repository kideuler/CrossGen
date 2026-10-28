#include "SKELETON/SKELETON.hxx"

#include <algorithm>
#include <cmath>
#include <map>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "triangle/TriangleMesher.hpp"

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

namespace {

// TMOP on SKELETON's trace mesh. The quads keep the order Step 10 gave them
// (keepOrientation): a quad folded right over has to look folded to the
// untangler. Boundary nodes slide on the trace mesh's own boundary polyline,
// not on a spline through it: some of them are off dS (Step 4's moves, the
// welds of zero-length edges), a spline through those bows further off
// still, and BaseComplex has to find every one of them on dS again
// afterwards (M7_14 left three sides beyond its reach). `pinned` and `held`
// are held; x is replaced, and true returned, only by an untangled result:
// one the smoother gave up on has moved its nodes part of the way, and the
// sides read off it zigzag worse than Step 10's (multimat/geom004, 1.2 ->
// 5.8 rad). run()'s own verdict is not the test: it is false whenever the
// worst element ends up worse than it started, which on a trace mesh is
// usual (geom035: 0.34 -> 0.24) and beside the point -- the grid lines are
// what is read off it, not its elements. `partway` takes the result however
// it came out, for a stage of stageOntoBoundary().
bool smoothTrace(std::vector<Point> &x, const std::vector<std::array<int, 4>> &quads,
                 const std::vector<int> &pinned, const std::vector<char> &held,
                 const mesh::TMOP::Options &opts, mesh::TMOP::Report &report,
                 bool partway = false) {
    try {
        mesh::QuadMesh::Options qo;
        qo.keepOrientation = true;
        qo.curveSource = mesh::QuadMesh::Options::CurvePolyline;
        mesh::QuadMesh trace(x, std::vector<mesh::Quad>(quads.begin(), quads.end()), {}, qo);
        for (const int v : pinned) trace.pinVertex(v);
        for (size_t v = 0; v < held.size(); ++v) {
            if (held[v]) trace.pinVertex(static_cast<int>(v));
        }
        trace.classifyNodes();
        trace.buildFeatureCurves();
        trace.computeSlideTangents();
        mesh::TMOP smoother(trace, opts);
        smoother.run();
        report = smoother.getReport();
        if (partway && trace.vertices.size() == x.size()) {
            x = trace.vertices;
            return report.invertedAfter == 0;
        }
        if (!report.ran || report.invertedAfter != 0 || trace.vertices.size() != x.size()) {
            report.ran = false;
            return false;
        }
        x = trace.vertices;
        return true;
    } catch (const std::exception &) {
        report.ran = false;
        return false;
    }
}

Point closestOn(const std::vector<Point> &L, const Point &p, double &dist) {
    Point best = p;
    dist = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < L.size(); ++i) {
        const Point &a = L[i], &b = L[(i + 1) % L.size()];
        const Point ab = b - a;
        const double len2 = dotP(ab, ab);
        const double t = len2 > 0.0 ? std::clamp(dotP(p - a, ab) / len2, 0.0, 1.0) : 0.0;
        const Point q = a + ab * t;
        const double d = normP(p - q);
        if (d < dist) {
            dist = d;
            best = q;
        }
    }
    return best;
}

// Where each boundary vertex of the trace mesh belongs on dS. Each loop of
// the trace mesh's boundary is laid on the loop of dS most of its vertices
// are closest to -- so that a vertex in a thin channel is not taken across
// it -- in order: the corners of the model it holds (`pinned`) stay where they
// are, and the vertices between two of them are spaced along dS between the
// corners as they are spaced along the trace mesh's own boundary. A loop with
// no corner is anchored at its vertex closest to dS. Not each vertex at its
// own closest point: where Step 10 folded the boundary that scrambles the
// vertices' order along dS (geom003, an exact disc), and a quad whose held
// boundary vertices are out of order is a fold no smoothing undoes. `held`
// marks the vertices on a boundary edge, `adjacent` every vertex's
// neighbours.
void boundaryTargets(const std::vector<Point> &x, const std::vector<std::array<int, 4>> &quads,
                     const std::vector<std::vector<Point>> &loops, const std::vector<int> &pinned,
                     std::vector<Point> &target, std::vector<char> &held,
                     std::vector<std::vector<int>> &adjacent) {
    const int nv = static_cast<int>(x.size());
    std::vector<std::vector<int>> adj(nv);
    adjacent.assign(nv, {});
    {
        std::map<std::pair<int, int>, int> use;
        for (const auto &q : quads) {
            for (int k = 0; k < 4; ++k) ++use[std::minmax(q[k], q[(k + 1) % 4])];
        }
        for (const auto &kv : use) {
            adjacent[kv.first.first].push_back(kv.first.second);
            adjacent[kv.first.second].push_back(kv.first.first);
            if (kv.second != 1) continue;
            adj[kv.first.first].push_back(kv.first.second);
            adj[kv.first.second].push_back(kv.first.first);
        }
    }
    std::vector<char> isPinned(nv, 0);
    for (const int v : pinned) {
        if (v >= 0 && v < nv) isPinned[v] = 1;
    }
    target = x;
    held.assign(nv, 0);
    std::vector<char> seen(nv, 0);
    for (int v0 = 0; v0 < nv; ++v0) {
        if (seen[v0] || adj[v0].empty()) continue;
        // The loop, in order; anything but a simple cycle is laid vertex by
        // vertex at its closest points instead.
        std::vector<int> cycle{v0};
        seen[v0] = 1;
        bool simple = adj[v0].size() == 2;
        for (int prev = -1, cur = v0;;) {
            int next = -1;
            for (const int w : adj[cur]) {
                if (w != prev && !seen[w]) {
                    next = w;
                    break;
                }
            }
            if (next < 0) break;
            if (adj[next].size() != 2) simple = false;
            seen[next] = 1;
            cycle.push_back(next);
            prev = cur;
            cur = next;
        }
        std::vector<int> votes(loops.size(), 0);
        for (const int v : cycle) {
            int best = 0;
            double bestD = std::numeric_limits<double>::infinity();
            for (size_t l = 0; l < loops.size(); ++l) {
                double d = 0.0;
                closestOn(loops[l], x[v], d);
                if (d < bestD) {
                    bestD = d;
                    best = static_cast<int>(l);
                }
            }
            ++votes[best];
        }
        const auto &L = loops[std::max_element(votes.begin(), votes.end()) - votes.begin()];
        for (const int v : cycle) held[v] = 1;
        const int n = static_cast<int>(L.size());
        std::vector<double> arc(n + 1, 0.0);
        for (int i = 0; i < n; ++i) arc[i + 1] = arc[i] + normP(L[(i + 1) % n] - L[i]);
        const double T = arc[n];
        auto locate = [&](const Point &p, double &dist) {
            double best = std::numeric_limits<double>::infinity(), at = 0.0;
            for (int i = 0; i < n; ++i) {
                const Point &a = L[i], ab = L[(i + 1) % n] - L[i];
                const double l2 = dotP(ab, ab);
                const double t = l2 > 0.0 ? std::clamp(dotP(p - a, ab) / l2, 0.0, 1.0) : 0.0;
                const double d = normP(p - (a + ab * t));
                if (d < best) {
                    best = d;
                    at = arc[i] + t * std::sqrt(l2);
                }
            }
            dist = best;
            return at;
        };
        auto pointAt = [&](double sv) {
            sv = std::fmod(sv, T);
            if (sv < 0.0) sv += T;
            const int i = std::clamp(static_cast<int>(std::upper_bound(arc.begin(), arc.end(), sv) - arc.begin()) - 1, 0, n - 1);
            const double len = arc[i + 1] - arc[i];
            return L[i] + (L[(i + 1) % n] - L[i]) * (len > 0.0 ? (sv - arc[i]) / len : 0.0);
        };
        const int c = static_cast<int>(cycle.size());
        std::vector<double> sv(c), dv(c);
        for (int i = 0; i < c; ++i) sv[i] = locate(x[cycle[i]], dv[i]);
        if (!simple || c < 3 || T <= 0.0) {
            for (int i = 0; i < c; ++i) target[cycle[i]] = pointAt(sv[i]);
            continue;
        }
        // Which way round dS the cycle runs: the way most of its steps go.
        int forward = 0;
        for (int i = 0; i < c; ++i) {
            const double step = std::fmod(sv[(i + 1) % c] - sv[i] + 2.0 * T, T);
            forward += step < 0.5 * T ? 1 : -1;
        }
        const double dir = forward >= 0 ? 1.0 : -1.0;
        std::vector<int> anchors;
        for (int i = 0; i < c; ++i) {
            if (isPinned[cycle[i]]) anchors.push_back(i);
        }
        if (anchors.empty())
            anchors.push_back(static_cast<int>(std::min_element(dv.begin(), dv.end()) - dv.begin()));
        for (size_t k = 0; k < anchors.size(); ++k) {
            const int i0 = anchors[k], i1 = anchors[(k + 1) % anchors.size()];
            const int count = ((i1 - i0 + c - 1) % c) + 1;   // steps from i0 to i1; c for one anchor
            double span = std::fmod(dir * (sv[i1] - sv[i0]) + 2.0 * T, T);
            if (anchors.size() == 1 || span <= 0.0) span = T;
            double chain = 0.0;
            std::vector<double> at(count + 1, 0.0);
            for (int j = 1; j <= count; ++j) {
                chain += normP(x[cycle[(i0 + j) % c]] - x[cycle[(i0 + j - 1) % c]]);
                at[j] = chain;
            }
            target[cycle[i0]] = isPinned[cycle[i0]] ? x[cycle[i0]] : pointAt(sv[i0]);
            for (int j = 1; j < count; ++j) {
                const double f = chain > 0.0 ? at[j] / chain : static_cast<double>(j) / count;
                target[cycle[(i0 + j) % c]] = pointAt(sv[i0] + dir * span * f);
            }
        }
    }
}

// Options::boundaryStages: every boundary vertex of the trace mesh moved a
// share of the way to its place on dS (boundaryTargets), the whole boundary
// held and TMOP untangling behind it, `stages` times, and then smoothed with
// the boundary still held. True, with x replaced, when that comes out
// untangled.
bool stageOntoBoundary(std::vector<Point> &x, const std::vector<std::array<int, 4>> &quads,
                       const std::vector<std::vector<Point>> &loops,
                       const std::vector<int> &pinned, int stages,
                       const mesh::TMOP::Options &opts, mesh::TMOP::Report &report) {
    std::vector<Point> target;
    std::vector<char> held;
    std::vector<std::vector<int>> adjacent;
    boundaryTargets(x, quads, loops, pinned, target, held, adjacent);
    const std::vector<Point> start = x;
    mesh::TMOP::Options stage = opts;
    stage.maxSweeps = std::min(opts.maxSweeps, 100);
    for (int s = 1; s <= stages; ++s) {
        const double t = static_cast<double>(s) / stages;
        for (size_t v = 0; v < x.size(); ++v) {
            if (held[v]) x[v] = start[v] + (target[v] - start[v]) * t;
        }
        mesh::TMOP::Report r;
        smoothTrace(x, quads, pinned, held, stage, r, true);   // tangled for now is allowed
    }
    return smoothTrace(x, quads, pinned, held, opts, report);
}

// The trace mesh's interior placed afresh, every interior vertex at the
// mean of its neighbours with the boundary on dS (boundaryTargets) and held
// -- Tutte's embedding, which is fold-free whenever dS is a convex loop --
// and smoothed from there with the boundary still held. For a Step 10 mesh
// too folded to untangle from where it is: an exact disc, whose axis is a
// point and whose caps Step 10 folds over by the hundred (geom003: 114
// quads). True, with x replaced, when that comes out untangled.
bool relaxOntoBoundary(std::vector<Point> &x, const std::vector<std::array<int, 4>> &quads,
                       const std::vector<std::vector<Point>> &loops,
                       const std::vector<int> &pinned, const mesh::TMOP::Options &opts,
                       mesh::TMOP::Report &report) {
    std::vector<Point> target;
    std::vector<char> held;
    std::vector<std::vector<int>> adjacent;
    boundaryTargets(x, quads, loops, pinned, target, held, adjacent);
    const int nv = static_cast<int>(x.size());
    std::vector<int> idx(nv, -1);
    int nf = 0;
    for (int v = 0; v < nv; ++v) {
        if (!held[v] && !adjacent[v].empty()) idx[v] = nf++;
    }
    if (nf > 0) {
        std::vector<Eigen::Triplet<double>> trip;
        Eigen::VectorXd bx = Eigen::VectorXd::Zero(nf), by = Eigen::VectorXd::Zero(nf);
        for (int v = 0; v < nv; ++v) {
            if (idx[v] < 0) continue;
            trip.emplace_back(idx[v], idx[v], static_cast<double>(adjacent[v].size()));
            for (const int w : adjacent[v]) {
                if (idx[w] >= 0) {
                    trip.emplace_back(idx[v], idx[w], -1.0);
                } else {
                    bx[idx[v]] += target[w][0];
                    by[idx[v]] += target[w][1];
                }
            }
        }
        Eigen::SparseMatrix<double> L(nf, nf);
        L.setFromTriplets(trip.begin(), trip.end());
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(L);
        if (ldlt.info() != Eigen::Success) return false;
        const Eigen::VectorXd sx = ldlt.solve(bx), sy = ldlt.solve(by);
        for (int v = 0; v < nv; ++v) {
            if (idx[v] >= 0) x[v] = Point{sx[idx[v]], sy[idx[v]]};
        }
    }
    for (int v = 0; v < nv; ++v) {
        if (held[v]) x[v] = target[v];
    }
    return smoothTrace(x, quads, pinned, held, opts, report);
}

}  // namespace

SKELETON::SKELETON(std::shared_ptr<Mesh> mesh_) : SKELETON(std::move(mesh_), Options()) {}

SKELETON::SKELETON(std::shared_ptr<Mesh> mesh_, const Options &opts)
    : mesh(std::move(mesh_)), options(opts) {}

SKELETON::~SKELETON() = default;

bool SKELETON::run(Stage last) {
    triangulate();
    if (last == Stage::Triangulation) return true;
    buildAxis();
    if (last == Stage::Axis) return true;
    buildMap();
    if (last == Stage::Map) return true;
    buildTMesh();
    if (last == Stage::TMesh) return true;
    if (!quantize()) return false;
    if (last == Stage::Quantized) return true;
    if (!extractQuads()) return false;
    if (last == Stage::QuadMesh) return true;
    if (!buildBlocks()) return false;
    if (last == Stage::Blocks) return true;
    if (!buildMesh()) return false;
    if (last == Stage::Mesh) return true;
    return smooth();
}

void SKELETON::discardFrom(Stage from) {
    // Latest first, as ZIPLINE does.
    if (from <= Stage::Smoothed) { smoothedMesh.reset(); tmopReport = mesh::TMOP::Report(); }
    if (from <= Stage::Mesh) blockMesh.reset();
    if (from <= Stage::Blocks) {
        blocks.reset();
        blockSmoothingReport = mesh::TMOP::Report();
        geometry = BlockGeometry::AsTraced;
    }
    if (from <= Stage::QuadMesh) quads.reset();
    if (from <= Stage::Quantized) { quant.reset(); quantReport = TMeshQuantizer::Report(); }
    if (from <= Stage::TMesh) tmesh.reset();
    if (from <= Stage::Map) map.reset();
    if (from <= Stage::Axis) axis.reset();
    if (from <= Stage::Triangulation) {
        triangulation.reset();
        loops.clear();
        corners.clear();
        status = Status();
    }
}

double SKELETON::h() const {
    if (options.targetSize > 0.0) return options.targetSize;
    return axis ? axis->defaultTargetSize() : 0.0;
}

// ---------------------------------------------------------------------------
// Steps 1-2: dS's own vertices, triangulated with no interior point, which is
// what Step 3.1 builds the dual complex on. The outer loop is the one that
// encloses the most area. Until 2026-10-17 it was taken to be the first loop
// the walk found, and on a model whose first boundary vertex sits on a hole
// that told the triangulator to fill the hole and carve out the domain.
// ---------------------------------------------------------------------------
void SKELETON::triangulate() {
    discardFrom(Stage::Triangulation);

    std::unordered_map<int, std::vector<int>> adjacent;
    for (const int be : mesh->boundaryEdges) {
        const int a = mesh->edges[be][0], b = mesh->edges[be][1];
        adjacent[a].push_back(b);
        adjacent[b].push_back(a);
    }
    std::unordered_set<int> visited;
    std::vector<std::vector<int>> vertexLoops;
    for (const int start : mesh->boundaryVertices) {
        if (visited.count(start)) continue;
        std::vector<int> loop;
        int prev = -1, cur = start;
        while (true) {
            visited.insert(cur);
            loop.push_back(cur);
            int next = -1;
            for (const int nb : adjacent[cur]) {
                if (nb != prev && !visited.count(nb)) {
                    next = nb;
                    break;
                }
            }
            if (next < 0) break;
            prev = cur;
            cur = next;
        }
        if (loop.size() >= 3) vertexLoops.push_back(std::move(loop));
    }
    if (vertexLoops.empty()) throw std::runtime_error("SKELETON: the model has no boundary loop");

    auto area = [&](const std::vector<int> &loop) {
        double a = 0.0;
        for (size_t i = 0; i < loop.size(); ++i)
            a += cross2(mesh->vertices[loop[i]], mesh->vertices[loop[(i + 1) % loop.size()]]);
        return 0.5 * a;
    };
    size_t outer = 0;
    for (size_t l = 1; l < vertexLoops.size(); ++l) {
        if (std::abs(area(vertexLoops[l])) > std::abs(area(vertexLoops[outer]))) outer = l;
    }
    std::swap(vertexLoops[0], vertexLoops[outer]);

    std::unordered_map<int, int> remap;
    std::vector<std::array<double, 2>> vertlist;
    std::vector<std::vector<std::array<int, 2>>> segmentLoops;
    std::vector<int> types;
    for (size_t l = 0; l < vertexLoops.size(); ++l) {
        for (const int v : vertexLoops[l]) {
            if (remap.emplace(v, static_cast<int>(vertlist.size())).second)
                vertlist.push_back(mesh->vertices[v]);
        }
        std::vector<std::array<int, 2>> segments;
        const auto &loop = vertexLoops[l];
        for (size_t i = 0; i < loop.size(); ++i)
            segments.push_back({remap[loop[i]], remap[loop[(i + 1) % loop.size()]]});
        segmentLoops.push_back(std::move(segments));
        types.push_back(l == 0 ? 0 : 1);
    }

    triangle_wrapper::TriangleMesher2D::Options to;
    to.just_delaunay = true;
    triangle_wrapper::TriangleMesher2D mesher(to);
    triangle_wrapper::TriangleMesher2D::MeshInput input;
    input.vertlist = vertlist;
    input.segment_loops = segmentLoops;
    input.type = types;
    input.h = 0.0;
    const auto output = mesher.triangulate(input);
    triangulation = std::make_shared<Mesh>(output.verts, output.triangles);

    // The loops as points, and their corners by the turning rule.
    const double threshold = options.cornerAngle * M_PI / 180.0;
    for (const auto &loop : vertexLoops) {
        std::vector<Point> pts;
        for (const int v : loop) pts.push_back(mesh->vertices[v]);
        const size_t n = pts.size();
        for (size_t i = 0; i < n; ++i) {
            const Point a = pts[i] - pts[(i + n - 1) % n], b = pts[(i + 1) % n] - pts[i];
            if (std::abs(std::atan2(cross2(a, b), dotP(a, b))) > threshold) corners.push_back(pts[i]);
        }
        loops.push_back(std::move(pts));
    }
    status.boundaryLoops = static_cast<int>(loops.size());
    status.corners = static_cast<int>(corners.size());
}

// ---------------------------------------------------------------------------
// Steps 3-4, on a copy of the triangulation: simplify() moves boundary
// vertices of the mesh it was built on, and the triangulation is kept as
// Steps 1-2 left it so that this stage can be run again.
// ---------------------------------------------------------------------------
void SKELETON::buildAxis() {
    if (!triangulation) triangulate();
    discardFrom(Stage::Axis);
    axis = std::make_shared<MedialAxis>(std::make_shared<Mesh>(*triangulation));
    status.rawMedialVertices = static_cast<int>(axis->medialVertices.size());
    status.rawMedialEdges = static_cast<int>(axis->medialEdges.size());
    SimplifyOptions so = options.simplify;
    so.targetSize = h();
    axis->simplify(so);
    status.targetSize = h();
}

void SKELETON::buildMap() {
    if (!axis) buildAxis();
    discardFrom(Stage::Map);
    map = std::make_unique<MedialAxisMap>(axis);
}

void SKELETON::buildTMesh() {
    if (!map) buildMap();
    discardFrom(Stage::TMesh);
    MedialAxisTMesh::Options to = options.tmesh;
    to.targetSize = h();
    tmesh = std::make_unique<MedialAxisTMesh>(*map, to);
}

bool SKELETON::quantize() {
    if (!tmesh) buildTMesh();
    discardFrom(Stage::Quantized);
    quant = std::make_unique<BlockQuant>(makeQuantTMesh(*tmesh, tmesh->targetSize()));
    if (!quant->ok) return false;
    TMeshQuantizer::Options qo = options.quantizer;
    qo.requirements = quant->requirements;
    quantReport = TMeshQuantizer(quant->tmesh, qo).run();
    return quantReport.consistent;
}

bool SKELETON::extractQuads() {
    if (!quant) quantize();
    discardFrom(Stage::QuadMesh);
    if (!quantizationOk()) return false;
    quads = std::make_unique<QuadExtraction>(extractQuadMesh(*quant));
    return !quads->quads.empty();
}

bool SKELETON::buildBlocks() {
    if (!quads && !extractQuads()) return false;
    discardFrom(Stage::Blocks);
    if (quads->quads.empty()) return false;
    QuadExtraction fine = options.blockSubdivisions > 1
                              ? extractQuadMesh(*quant, options.blockSubdivisions)
                              : *quads;
    // The model's corners, as vertices of the trace mesh, held by every step
    // below. mesh::QuadMesh fixes a boundary node turning by more than 45
    // degrees and lets the rest slide; a corner by the 20-degree rule that
    // slid off would no longer be found on a vertex, and the blocks would lose
    // the separatrices it starts (geom034: 115 interior sides -> 99).
    std::vector<int> cornerVertices;
    {
        const double tol = blockOptions().tolerance * axis->boundingBoxDiagonal();
        for (const Point &c : corners) {
            int best = -1;
            double bestD = tol;
            for (size_t v = 0; v < fine.vertices.size(); ++v) {
                if (!fine.onBoundary[v]) continue;
                const double d = normP(fine.vertices[v] - c);
                if (d <= bestD) {
                    bestD = d;
                    best = static_cast<int>(v);
                }
            }
            if (best >= 0) cornerVertices.push_back(best);
        }
    }
    const std::vector<Point> traced = fine.vertices;

    // The geometry the blocks are read off, best first: the trace mesh
    // smoothed with its boundary sliding; failing that -- tangled, or a
    // block side across dS -- with its boundary laid onto dS in stages and
    // held; failing that -- still tangled -- relaxed afresh from dS inwards.
    // Of those that untangle, the one fewest sides cross; as Step 10 left it
    // when none does.
    struct Candidate {
        BlockGeometry kind;
        std::unique_ptr<BaseComplex> blocks;
        mesh::TMOP::Report report;
    };
    std::vector<Candidate> found;
    auto attempt = [&](BlockGeometry kind, std::vector<Point> &x, bool untangled,
                       const mesh::TMOP::Report &report) {
        if (!untangled) return false;
        found.push_back({kind, std::make_unique<BaseComplex>(x, fine.quads, loops, corners,
                                                             blockOptions()),
                         report});
        return true;
    };
    // A pinched vertex -- two fans of quads touching at a concave corner
    // (tooth, M7_20) -- does not stop the smoother; an edge of three quads
    // would.
    const bool smoothing = options.blockSmoothing.maxSweeps > 0 && fine.nonManifoldEdges == 0;
    bool clean = false;
    if (smoothing) {
        std::vector<Point> x = traced;
        mesh::TMOP::Report report;
        if (attempt(BlockGeometry::Smoothed, x,
                    smoothTrace(x, fine.quads, cornerVertices, {}, options.blockSmoothing, report),
                    report))
            clean = found.back().blocks->getReport().sidesCrossingBoundary == 0;
        else
            blockSmoothingReport = report;
    }
    if (smoothing && !clean && options.boundaryStages > 0) {
        std::vector<Point> x = traced;
        mesh::TMOP::Report report;
        attempt(BlockGeometry::Staged, x,
                stageOntoBoundary(x, fine.quads, loops, cornerVertices, options.boundaryStages,
                                  options.blockSmoothing, report),
                report);
        if (found.empty()) {
            std::vector<Point> y = traced;
            mesh::TMOP::Report relaxed;
            attempt(BlockGeometry::Relaxed, y,
                    relaxOntoBoundary(y, fine.quads, loops, cornerVertices, options.blockSmoothing,
                                      relaxed),
                    relaxed);
        }
    }
    int best = -1;
    for (size_t i = 0; i < found.size(); ++i) {
        if (best < 0 || found[i].blocks->getReport().sidesCrossingBoundary <
                            found[best].blocks->getReport().sidesCrossingBoundary)
            best = static_cast<int>(i);
    }
    if (best >= 0) {
        blocks = std::move(found[best].blocks);
        blockSmoothingReport = found[best].report;
        geometry = found[best].kind;
    } else {
        blocks = std::make_unique<BaseComplex>(traced, fine.quads, loops, corners, blockOptions());
        blockSmoothingReport.ran = false;
        geometry = BlockGeometry::AsTraced;
    }
    return !blocks->decomposition().blocks.empty();
}

BaseComplex::Options SKELETON::blockOptions() const {
    BaseComplex::Options bo = options.blocks;
    if (axis && axis->boundingBoxDiagonal() > 0.0)
        bo.snapTolerance = std::max(bo.snapTolerance, 0.5 * h() / axis->boundingBoxDiagonal());
    return bo;
}

bool SKELETON::buildMesh(const BlockQuadMesh::Options &opts) {
    options.mesh = opts;
    return buildMesh();
}

bool SKELETON::buildMesh() {
    if (!blocks && !buildBlocks()) return false;
    discardFrom(Stage::Mesh);
    if (blocks->decomposition().blocks.empty()) return false;
    blockMesh = std::make_unique<BlockQuadMesh>(blocks->decomposition(), options.mesh);
    return true;
}

// As ZIPLINE smooths: dS nodes slide on their curves, nothing is pinned --
// with one difference. BlockQuadMesh's elements are counter-clockwise by
// construction, a transfinite grid on a bent block folds some of them right
// over, and mesh::QuadMesh by default turns any element of negative area
// round: TMOP then sees nothing to untangle, its report counts no fold, and
// the fold stays (geom012: 3 folds reported, 39 there). keepOrientation keeps
// them as they are, folded and visible.
bool SKELETON::smooth() {
    if (!blockMesh && !buildMesh()) return false;
    discardFrom(Stage::Smoothed);
    if (blockMesh->getReport().quads <= 0) return false;
    mesh::QuadMesh::Options qo;
    qo.keepOrientation = true;
    smoothedMesh = std::make_unique<mesh::QuadMesh>(mesh::QuadMesh::from(*blockMesh, qo));
    mesh::TMOP smoother(*smoothedMesh, options.tmop);
    smoother.run();
    tmopReport = smoother.getReport();
    return true;
}
