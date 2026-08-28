#include "MERIDIAN/Arrangement.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace {

// The local indices of face `t`, rotated so that vertex `v` is first and the
// three run counter-clockwise. `o` is the counter-clockwise order of the face.
bool cornerOf(const Triangle &t, const std::array<int, 3> &o, int v,
              int &iv, int &ib, int &ic) {
    for (int k = 0; k < 3; ++k) {
        if (t[o[k]] == v) { iv = o[k]; ib = o[(k + 1) % 3]; ic = o[(k + 2) % 3]; return true; }
    }
    return false;
}

double signedArea(const std::vector<Point> &loop) {
    double a = 0.0;
    for (size_t i = 0; i < loop.size(); ++i) {
        const Point &p = loop[i];
        const Point &q = loop[(i + 1) % loop.size()];
        a += cross2(p, q);
    }
    return 0.5 * a;
}

double pointSegment(const Point &p, const Point &a, const Point &b) {
    const Point d = b - a;
    const double dd = dotP(d, d);
    if (dd <= 0.0) return normP(p - a);
    double s = dotP(p - a, d) / dd;
    s = std::max(0.0, std::min(1.0, s));
    return normP(p - (a + d * s));
}

// Barycentric points are affine along a step, so a curve is clipped by
// interpolating the two triples.
std::array<double, 3> lerpBary(const std::array<double, 3> &a,
                               const std::array<double, 3> &b, double s) {
    return {{a[0] + (b[0] - a[0]) * s, a[1] + (b[1] - a[1]) * s, a[2] + (b[2] - a[2]) * s}};
}

} // namespace

// ---------------------------------------------------------------------------
Arrangement::Arrangement(const Separatrices &separatrices, const SubdomainLabels &labels)
    : Arrangement(separatrices, labels, Options()) {}

Arrangement::Arrangement(const Separatrices &separatrices, const SubdomainLabels &labels,
                         const Options &opts)
    : sep(&separatrices), lab(&labels), imm(&separatrices.getImmersion()),
      orig(&separatrices.getImmersion().getOriginalMesh()),
      cut(&separatrices.getCutMesh()), uv(&separatrices.getUV()), options(opts) {
    if (&labels.getImmersion() != imm) {
        throw std::runtime_error("Arrangement: the labels and the curves belong to "
                                 "different immersions");
    }
    if (uv->size() != cut->vertices.size()) {
        throw std::runtime_error("Arrangement: the map does not cover Omega");
    }

    buildDomain();
    dedupeCurves();
    buildEndNodes();
    buildBoundaryNodes();
    buildInterfaceNodes();
    collectSegments();
    findCrossings();
    findInterfaceCrossings();
    trimUnresolved();
    splitCurves();
    buildBoundaryArcs();
    buildInterfaceArcs();
    collapseShortArcs();
    buildAngles();
    linkHalfEdges();
    extractFaces();
    classifyFaces();
    classifyMaterials();
    if (options.checkFeatures) checkFeatures();
    check();
}

// ---------------------------------------------------------------------------
// buildDomain()
//
// The parts of S the arrangement is built against: how big it is, how its
// triangles are oriented, and dS directed so that the model is always on the
// left. That last one is what separates a patch from a hole later on without
// any point-in-polygon test -- see the class comment.
// ---------------------------------------------------------------------------
void Arrangement::buildDomain() {
    const int nV = static_cast<int>(orig->vertices.size());
    const int nF = static_cast<int>(orig->triangles.size());

    Point lo{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    Point hi{-lo[0], -lo[1]};
    for (const Point &p : orig->vertices) {
        lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
        hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
    }
    modelExtent = normP(hi - lo);
    if (!(modelExtent > 0.0)) modelExtent = 1.0;
    mergeTol = options.mergeTolerance * modelExtent;
    cellSize = std::max(mergeTol * 4.0, modelExtent * 1e-9);

    report.domainArea = 0.0;
    for (int f = 0; f < nF; ++f) {
        const Triangle &t = orig->triangles[f];
        report.domainArea += std::fabs(0.5 * cross2(orig->vertices[t[1]] - orig->vertices[t[0]],
                                                    orig->vertices[t[2]] - orig->vertices[t[0]]));
    }

    vertexNode.assign(nV, -1);
    fanBuilt.assign(nV, 0);
    fanCache.assign(nV, {});
    fanOffset.assign(nV, {});
    fanTotal.assign(nV, 0.0);
    boundaryOut.assign(nV, -1);

    for (int e : orig->boundaryEdges) {
        const int f = orig->edgeTriangles[e][0] >= 0 ? orig->edgeTriangles[e][0]
                                                     : orig->edgeTriangles[e][1];
        if (f < 0) continue;
        const std::array<int, 3> o = ccw(f);
        const Triangle &t = orig->triangles[f];
        const int x = orig->edges[e][0], y = orig->edges[e][1];
        // In the counter-clockwise ordering the model is to the left of every
        // directed side of the triangle, so the boundary edge is directed the
        // way the triangle uses it.
        DirEdge de;
        de.edge = e;
        de.face = f;
        for (int k = 0; k < 3; ++k) {
            const int a = t[o[k]], b = t[o[(k + 1) % 3]];
            if ((a == x && b == y) || (a == y && b == x)) { de.a = a; de.b = b; break; }
        }
        if (de.a < 0) continue;
        boundaryOut[de.a] = static_cast<int>(dirBoundary.size());
        dirBoundary.push_back(de);
    }
    boundaryEvents.assign(dirBoundary.size(), {});
}

std::array<int, 3> Arrangement::ccw(int f) const {
    const Triangle &t = orig->triangles[f];
    const Point &a = orig->vertices[t[0]];
    const Point &b = orig->vertices[t[1]];
    const Point &c = orig->vertices[t[2]];
    if (cross2(b - a, c - a) >= 0.0) return std::array<int, 3>{{0, 1, 2}};
    return std::array<int, 3>{{0, 2, 1}};
}

double Arrangement::imageAngle(int f, int v) const {
    int iv = -1, ib = -1, ic = -1;
    if (!cornerOf(orig->triangles[f], ccw(f), v, iv, ib, ic)) return 0.0;
    const Triangle &ct = cut->triangles[f];
    const Point A = (*uv)[ct[ib]] - (*uv)[ct[iv]];
    const Point B = (*uv)[ct[ic]] - (*uv)[ct[iv]];
    double a = std::atan2(cross2(A, B), dotP(A, B));
    if (a < 0.0) a += 2.0 * M_PI;
    return a;
}

// ---------------------------------------------------------------------------
// vertexFan()
//
// The faces around v, counter-clockwise. For a boundary vertex the fan is a
// run rather than a cycle and has to start at the right end: index 0 is the
// face whose first edge at v is the boundary edge leaving v with the model on
// its left. Every layout angle at a boundary node is then measured from that
// edge, and the sector inside the model is [0, totalAngle] with no wrap in it.
// ---------------------------------------------------------------------------
bool Arrangement::vertexFan(int v, std::vector<int> &fan) const {
    fan.clear();
    const auto &vt = orig->vertexTriangles;
    if (v < 0 || v + 1 >= static_cast<int>(vt.rowPtr.size())) return false;
    const int deg = vt.rowPtr[v + 1] - vt.rowPtr[v];
    if (deg <= 0) return false;

    auto across = [&](int f, int a, int b) {
        for (int k = 0; k < 3; ++k) {
            const int e = orig->triangleEdges[f][k];
            const int x = orig->edges[e][0], y = orig->edges[e][1];
            if ((x == a && y == b) || (x == b && y == a)) {
                const int t0 = orig->edgeTriangles[e][0];
                const int t1 = orig->edgeTriangles[e][1];
                return (t0 == f) ? t1 : t0;
            }
        }
        return -1;
    };

    const int f0 = vt.colIdx[vt.rowPtr[v]];
    int f = f0;
    for (int guard = 0; guard <= deg; ++guard) {   // clockwise, to the fan's start
        int iv, ib, ic;
        if (!cornerOf(orig->triangles[f], ccw(f), v, iv, ib, ic)) return false;
        const int g = across(f, v, orig->triangles[f][ib]);
        if (g < 0 || g == f0) break;
        f = g;
    }
    const int start = f;
    for (int guard = 0; guard <= deg; ++guard) {   // counter-clockwise, collecting
        fan.push_back(f);
        int iv, ib, ic;
        if (!cornerOf(orig->triangles[f], ccw(f), v, iv, ib, ic)) return false;
        const int g = across(f, v, orig->triangles[f][ic]);
        if (g < 0 || g == start) break;
        f = g;
    }
    return static_cast<int>(fan.size()) == deg;
}

void Arrangement::buildFan(int v) {
    if (v < 0 || fanBuilt[v]) return;
    fanBuilt[v] = 1;
    std::vector<int> fan;
    if (!vertexFan(v, fan)) { fanCache[v].clear(); fanTotal[v] = 0.0; return; }
    fanCache[v] = fan;
    fanOffset[v].assign(fan.size(), 0.0);
    double acc = 0.0;
    for (size_t i = 0; i < fan.size(); ++i) {
        fanOffset[v][i] = acc;
        acc += imageAngle(fan[i], v);
    }
    fanTotal[v] = acc;
}

// ---------------------------------------------------------------------------
long long Arrangement::cellOf(const Point &p) const {
    const long long x = static_cast<long long>(std::floor(p[0] / cellSize));
    const long long y = static_cast<long long>(std::floor(p[1] / cellSize));
    return x * 73856093LL ^ (y * 19349663LL);
}

int Arrangement::addNode(const Point &p, NodeKind kind) {
    Node n;
    n.p = p;
    n.kind = kind;
    const int id = static_cast<int>(nodes.size());
    nodes.push_back(n);
    nodeRef.push_back(Point{1.0, 0.0});
    nodeGrid[cellOf(p)].push_back(id);
    return id;
}

int Arrangement::findOrAddCrossing(const Point &p, int face) {
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            const Point q{p[0] + dx * cellSize, p[1] + dy * cellSize};
            auto it = nodeGrid.find(cellOf(q));
            if (it == nodeGrid.end()) continue;
            for (int id : it->second) {
                if (normP(nodes[id].p - p) <= mergeTol) {
                    ++report.mergedNodes;
                    return id;
                }
            }
        }
    }
    const int id = addNode(p, NodeKind::Crossing);
    nodes[id].face = face;
    return id;
}

// ---------------------------------------------------------------------------
// dedupeCurves()
//
// Q5 closes a separatrix at a cone, and the same line is a separatrix of *that*
// cone too, so a finished layout traces every cone-to-cone edge twice, once
// from each end. Both copies are the same layout edge. Left alone they would
// double the arcs at every cone the layout actually got right -- so the models
// that work would be the ones reported as broken -- and would bound a face of
// zero area between them.
//
// They are matched by geometry rather than by their endpoints alone, because
// two cones can genuinely be joined by two different arcs going two different
// ways; what says these two are one arc is that they are the same length and
// lie on top of each other in between. The comparison is taken away from the
// ends on purpose: Stage 7 snapped each copy onto the cone it arrived at, so it
// is exactly at the ends that the two disagree, by the snap gap.
// ---------------------------------------------------------------------------
void Arrangement::dedupeCurves() {
    const std::vector<Separatrices::Curve> &curves = sep->curves();
    curveDropped.assign(curves.size(), 0);

    auto poly = [&](int i) {
        std::vector<Point> out;
        const Separatrices::Curve &c = curves[i];
        if (c.steps.empty()) return out;
        out.push_back(sep->point(c.steps.front(), false, Separatrices::Space::Model));
        for (const Separatrices::Step &st : c.steps) {
            out.push_back(sep->point(st, true, Separatrices::Space::Model));
        }
        return out;
    };
    auto sampleAt = [](const std::vector<Point> &p, const std::vector<double> &acc, double s) {
        const double want = s * acc.back();
        size_t k = 1;
        while (k + 1 < acc.size() && acc[k] < want) ++k;
        const double d = acc[k] - acc[k - 1];
        const double f = d > 0.0 ? (want - acc[k - 1]) / d : 0.0;
        return p[k - 1] + (p[k] - p[k - 1]) * f;
    };

    std::unordered_map<long long, std::vector<int>> bucket;
    for (size_t i = 0; i < curves.size(); ++i) {
        const Separatrices::Curve &c = curves[i];
        if (c.end != Separatrices::End::Cone || c.cone < 0 || c.toCone < 0) continue;
        if (c.steps.empty()) continue;
        const long long lo = std::min(c.cone, c.toCone), hi = std::max(c.cone, c.toCone);
        bucket[lo * 1000003LL + hi].push_back(static_cast<int>(i));
    }

    for (auto &kv : bucket) {
        std::vector<int> &v = kv.second;
        for (size_t a = 0; a < v.size(); ++a) {
            if (curveDropped[v[a]]) continue;
            const std::vector<Point> pa = poly(v[a]);
            std::vector<double> aa(pa.size(), 0.0);
            for (size_t k = 1; k < pa.size(); ++k) aa[k] = aa[k - 1] + normP(pa[k] - pa[k - 1]);
            if (!(aa.back() > 0.0)) continue;

            for (size_t b = a + 1; b < v.size(); ++b) {
                if (curveDropped[v[b]]) continue;
                if (curves[v[b]].cone != curves[v[a]].toCone) continue;
                if (curves[v[b]].toCone != curves[v[a]].cone) continue;

                const std::vector<Point> pb = poly(v[b]);
                std::vector<double> ab(pb.size(), 0.0);
                for (size_t k = 1; k < pb.size(); ++k) ab[k] = ab[k - 1] + normP(pb[k] - pb[k - 1]);
                if (!(ab.back() > 0.0)) continue;

                const double scale = std::max(aa.back(), ab.back());
                const double tol = options.duplicateTolerance * scale;
                if (std::fabs(aa.back() - ab.back()) > tol) continue;
                bool same = true;
                for (double s : {0.25, 0.5, 0.75}) {
                    if (normP(sampleAt(pa, aa, s) - sampleAt(pb, ab, 1.0 - s)) > tol) {
                        same = false;
                        break;
                    }
                }
                if (!same) continue;
                curveDropped[v[b]] = 1;
                ++report.duplicateCurves;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// buildEndNodes()
//
// The cones first, because every other node that lands on one has to find it
// rather than duplicate it; then the two ends of every traced curve. A curve
// that Q5 closed ends on a cone that is already a node. One that left through
// dS puts a node on the boundary edge it crossed, at the parameter it crossed
// it. One that did neither leaves a dangling node, which is not a layout node
// at all and is counted as the defect it is.
// ---------------------------------------------------------------------------
void Arrangement::buildEndNodes() {
    const std::vector<int> &cv = sep->coneVertices();
    const std::vector<int> &ci = sep->coneIndices();
    coneNode.assign(cv.size(), -1);
    for (size_t s = 0; s < cv.size(); ++s) {
        const int v = cv[s];
        if (v < 0 || v >= static_cast<int>(orig->vertices.size())) continue;
        const int id = addNode(orig->vertices[v], NodeKind::Cone);
        Node &n = nodes[id];
        n.cone = static_cast<int>(s);
        n.index = ci[s];
        n.vertex = v;
        n.onBoundary = orig->isBoundaryVertex[v];
        n.valence = n.onBoundary ? (3 - n.index) : (4 - n.index);
        coneNode[s] = id;
        vertexNode[v] = id;
    }

    const std::vector<Separatrices::Curve> &curves = sep->curves();
    curveEvents.assign(curves.size(), {});

    for (size_t c = 0; c < curves.size(); ++c) {
        const Separatrices::Curve &cu = curves[c];
        if (cu.steps.empty() || curveDropped[c]) continue;
        if (!options.keepUnresolved && !Separatrices::resolved(cu.end)) continue;

        const int start = (cu.cone >= 0 && cu.cone < static_cast<int>(coneNode.size()))
                              ? coneNode[cu.cone] : -1;
        if (start < 0) continue;
        curveEvents[c].push_back(Event{0.0, start});

        const double last = static_cast<double>(cu.steps.size());
        int endNode = -1;
        if (cu.end == Separatrices::End::Cone && cu.toCone >= 0 &&
            cu.toCone < static_cast<int>(coneNode.size())) {
            endNode = coneNode[cu.toCone];
            const Point hit = sep->point(cu.steps.back(), true, Separatrices::Space::Model);
            report.maxConeSnapGap = std::max(report.maxConeSnapGap,
                                             normP(hit - nodes[endNode].p) / modelExtent);
        } else if (cu.end == Separatrices::End::Boundary) {
            const Separatrices::Step &st = cu.steps.back();
            const Triangle &t = orig->triangles[st.face];
            int l = 0;
            for (int k = 1; k < 3; ++k) if (st.exit[k] < st.exit[l]) l = k;
            const int m = (l + 1) % 3, n2 = (l + 2) % 3;
            const double den = st.exit[m] + st.exit[n2];
            const double par = den > 0.0 ? st.exit[n2] / den : 0.0;
            const int e = orig->triangleEdges[st.face][m];

            if (par <= 1e-9 || par >= 1.0 - 1e-9) {
                const int v = par <= 1e-9 ? t[m] : t[n2];
                if (vertexNode[v] >= 0) {
                    endNode = vertexNode[v];
                } else {
                    endNode = addNode(orig->vertices[v], NodeKind::BoundaryHit);
                    nodes[endNode].vertex = v;
                    nodes[endNode].onBoundary = true;
                    vertexNode[v] = endNode;
                }
            } else if (e >= 0 && orig->isBoundaryEdge[e]) {
                const Point p = orig->vertices[t[m]] * (1.0 - par) + orig->vertices[t[n2]] * par;
                endNode = findOrAddCrossing(p, st.face);
                if (nodes[endNode].kind == NodeKind::Crossing) {
                    nodes[endNode].kind = NodeKind::BoundaryHit;
                    nodes[endNode].onBoundary = true;
                    nodes[endNode].edge = e;
                    nodes[endNode].face = st.face;
                    const int be = boundaryOut[t[m]] >= 0 && dirBoundary[boundaryOut[t[m]]].edge == e
                                       ? boundaryOut[t[m]] : boundaryOut[t[n2]];
                    if (be >= 0 && dirBoundary[be].edge == e) {
                        nodes[endNode].t = (dirBoundary[be].a == t[m]) ? par : 1.0 - par;
                        nodeRef[endNode] = dirEdgeImage(be);
                        boundaryEvents[be].push_back(Event{nodes[endNode].t, endNode});
                    }
                }
            }
        }
        if (endNode < 0) {
            const Point p = sep->point(cu.steps.back(), true, Separatrices::Space::Model);
            endNode = addNode(p, NodeKind::Dangling);
            nodes[endNode].face = cu.steps.back().face;
        }
        curveEvents[c].push_back(Event{last, endNode});
    }
}

// ---------------------------------------------------------------------------
// buildBoundaryNodes()
//
// Q3 puts every component of dS - G on a line of constant u or constant v, so
// the label of dS changes exactly where the layout has a corner. Stage 1 is
// supposed to have put a cone at each of them; where it did not, the layout
// still has the corner and the arrangement had better carry a node for it, or
// a patch two of whose sides are one arc comes out with three corners and no
// diagnosis. Every such node is reported.
// ---------------------------------------------------------------------------
void Arrangement::buildBoundaryNodes() {
    // dS labels, indexed by the edge of S rather than by the edge of Omega.
    const std::vector<int> &toOrig = imm->getCut().getCutVertexToOriginal();
    std::unordered_map<MeshEdgeKey, int, MeshEdgeKeyHash> label;
    for (const SubdomainLabels::BoundaryEdge &be : lab->boundaryEdges()) {
        if (be.a < 0 || be.b < 0) continue;
        const int a = toOrig[be.a], b = toOrig[be.b];
        if (a < 0 || b < 0) continue;
        label[MeshEdgeKey(a, b)] = static_cast<int>(be.label);
    }
    edgeLabel = label;

    if (!options.boundaryCornerNodes) return;
    for (const DirEdge &de : dirBoundary) {
        const int v = de.b;                       // where this edge ends
        if (vertexNode[v] >= 0) continue;
        buildFan(v);
        // A regular point of dS has image angle pi exactly: Q3 put the boundary
        // on a coordinate line and Stage 3 drove its curvature to zero. Anything
        // else is a quarter turn the cone set does not know about, and it is a
        // corner of the layout whether or not Stage 1 named it one.
        if (std::fabs(fanTotal[v] - M_PI) <= options.cornerTolerance) continue;
        const int id = addNode(orig->vertices[v], NodeKind::BoundaryCorner);
        nodes[id].vertex = v;
        nodes[id].onBoundary = true;
        vertexNode[v] = id;
    }
}

// ---------------------------------------------------------------------------
void Arrangement::collectSegments() {
    faceSegs.assign(orig->triangles.size(), {});
    const std::vector<Separatrices::Curve> &curves = sep->curves();
    const double tiny = 1e-12 * modelExtent;
    for (size_t c = 0; c < curves.size(); ++c) {
        if (curveEvents[c].empty()) continue;
        const Separatrices::Curve &cu = curves[c];
        for (size_t i = 0; i < cu.steps.size(); ++i) {
            Seg s;
            s.curve = static_cast<int>(c);
            s.step = static_cast<int>(i);
            s.a = sep->point(cu.steps[i], false, Separatrices::Space::Model);
            s.b = sep->point(cu.steps[i], true, Separatrices::Space::Model);
            if (normP(s.b - s.a) <= tiny) continue;
            faceSegs[cu.steps[i].face].push_back(static_cast<int>(segs.size()));
            segs.push_back(s);
        }
    }
}

// ---------------------------------------------------------------------------
// findCrossings()
//
// Sec. 4's "crossings are computed triangle-locally". Two separatrices meet
// where two of their steps meet inside one face of S, and that is the only
// place they can meet, because a step is the whole of the curve inside that
// face. Doing it this way is not an optimisation: Psi overlaps itself, so the
// same question asked of the image would answer yes for two curves on opposite
// sides of the model.
//
// A crossing that lands on an edge of the triangulation is found twice, once
// from each side, and the merge tolerance is what makes those one node.
// ---------------------------------------------------------------------------
void Arrangement::findCrossings() {
    for (size_t f = 0; f < faceSegs.size(); ++f) {
        const std::vector<int> &list = faceSegs[f];
        for (size_t i = 0; i + 1 < list.size(); ++i) {
            for (size_t j = i + 1; j < list.size(); ++j) {
                const Seg &p = segs[list[i]];
                const Seg &q = segs[list[j]];
                // Consecutive steps of one curve share an endpoint by
                // construction; that is not a crossing.
                if (p.curve == q.curve && std::abs(p.step - q.step) <= 1) continue;

                const Point r = p.b - p.a;
                const Point s = q.b - q.a;
                const double den = cross2(r, s);
                const double scale = normP(r) * normP(s);
                if (std::fabs(den) <= 1e-12 * scale) { ++report.parallelOverlaps; continue; }
                const Point w = q.a - p.a;
                const double tp = cross2(w, s) / den;
                const double tq = cross2(w, r) / den;
                const double ep = mergeTol / std::max(normP(r), 1e-300);
                const double eq = mergeTol / std::max(normP(s), 1e-300);
                if (tp < -ep || tp > 1.0 + ep) continue;
                if (tq < -eq || tq > 1.0 + eq) continue;

                const Point hit = p.a + r * std::max(0.0, std::min(1.0, tp));
                const int id = findOrAddCrossing(hit, static_cast<int>(f));
                if (p.curve == q.curve) ++report.selfCrossings;
                curveEvents[p.curve].push_back(
                    Event{p.step + std::max(0.0, std::min(1.0, tp)), id});
                curveEvents[q.curve].push_back(
                    Event{q.step + std::max(0.0, std::min(1.0, tq)), id});
            }
        }
    }
}

// ---------------------------------------------------------------------------
// splitCurves()
//
// Every curve, cut at the nodes on it. The clipping is exact: phi is affine on
// a face, so a step is a straight segment in barycentric coordinates and a
// fraction of the way along it is the same fraction of the way between the two
// triples. Nothing is re-projected and nothing is searched for.
// ---------------------------------------------------------------------------
// trimUnresolved()
//
// A curve Q5 did not close, ended at the first layout edge it met.
//
// "Layout edge" means an arc of a curve that *is* closed, or dS. A crossing
// with another unclosed curve will not do: two curves that both wind and are
// both cut at their mutual crossing leave a node of degree two in the middle of
// a face and no more of a subdivision than the loose ends they started with.
// So the nodes worth stopping at are marked from the closed curves and from the
// boundary first, and the winding curve is then cut at the first of them it
// reaches.
//
// A curve that meets none of them is left alone. Nothing is gained by moving
// its loose end somewhere else, and the diagnosis -- danglingNodes, and the
// message that goes with it -- is the one the caller needs.
// ---------------------------------------------------------------------------
void Arrangement::trimUnresolved() {
    if (!options.trimUnresolvedAtCrossings) return;
    const std::vector<Separatrices::Curve> &curves = sep->curves();

    std::vector<char> onLayout(nodes.size(), 0);
    for (size_t c = 0; c < curves.size(); ++c) {
        if (curveDropped[c] || !Separatrices::resolved(curves[c].end)) continue;
        for (const Event &e : curveEvents[c]) {
            if (e.node >= 0 && e.node < static_cast<int>(onLayout.size())) onLayout[e.node] = 1;
        }
    }
    for (size_t n = 0; n < nodes.size(); ++n) {
        if (nodes[n].onBoundary || nodes[n].kind == NodeKind::Cone) onLayout[n] = 1;
    }

    for (size_t c = 0; c < curves.size(); ++c) {
        if (curveDropped[c] || Separatrices::resolved(curves[c].end)) continue;
        std::vector<Event> &ev = curveEvents[c];
        if (ev.size() < 2) continue;
        std::sort(ev.begin(), ev.end(),
                  [](const Event &a, const Event &b) { return a.param < b.param; });

        size_t cut = ev.size();
        for (size_t k = 1; k < ev.size(); ++k) {
            if (ev[k].param <= 1e-9) continue;
            const int n = ev[k].node;
            if (n < 0 || n >= static_cast<int>(onLayout.size())) continue;
            if (nodes[n].kind == NodeKind::Dangling) continue;
            if (!onLayout[n]) continue;
            cut = k;
            break;
        }
        if (cut + 1 >= ev.size()) continue;   // nothing to cut away

        ev.resize(cut + 1);
        ++report.trimmedCurves;
    }
}

// ---------------------------------------------------------------------------
void Arrangement::splitCurves() {
    const std::vector<Separatrices::Curve> &curves = sep->curves();
    const double tiny = 1e-12 * modelExtent;

    for (size_t c = 0; c < curves.size(); ++c) {
        std::vector<Event> &ev = curveEvents[c];
        if (ev.size() < 2) continue;
        const Separatrices::Curve &cu = curves[c];
        std::sort(ev.begin(), ev.end(),
                  [](const Event &a, const Event &b) { return a.param < b.param; });
        // Two events at the same place -- a crossing found from both sides of
        // an edge, or a crossing that fell on a cone -- are one event.
        std::vector<Event> keep;
        for (const Event &e : ev) {
            if (!keep.empty() && (e.node == keep.back().node ||
                                  e.param - keep.back().param < 1e-9)) {
                continue;
            }
            keep.push_back(e);
        }
        if (keep.size() < 2) continue;

        for (size_t k = 0; k + 1 < keep.size(); ++k) {
            const double p0 = keep[k].param, p1 = keep[k + 1].param;
            Arc arc;
            arc.kind = ArcKind::Separatrix;
            arc.curve = static_cast<int>(c);
            arc.from = keep[k].node;
            arc.to = keep[k + 1].node;

            int i0 = static_cast<int>(std::floor(p0)), i1 = static_cast<int>(std::floor(p1));
            double s0 = p0 - i0, s1 = p1 - i1;
            const int last = static_cast<int>(cu.steps.size()) - 1;
            if (i0 > last) { i0 = last; s0 = 1.0; }
            if (i1 > last) { i1 = last; s1 = 1.0; }
            for (int i = i0; i <= i1; ++i) {
                Separatrices::Step st = cu.steps[i];
                if (i == i0) st.entry = lerpBary(cu.steps[i].entry, cu.steps[i].exit, s0);
                if (i == i1) st.exit  = lerpBary(cu.steps[i].entry, cu.steps[i].exit, s1);
                const Point a = sep->point(st, false, Separatrices::Space::Model);
                const Point b = sep->point(st, true, Separatrices::Space::Model);
                if (normP(b - a) <= tiny && arc.steps.size() + 1 < cu.steps.size()) continue;
                arc.steps.push_back(st);
            }
            if (arc.steps.empty()) continue;

            arc.points.push_back(sep->point(arc.steps.front(), false, Separatrices::Space::Model));
            for (const Separatrices::Step &st : arc.steps) {
                const Point p = sep->point(st, true, Separatrices::Space::Model);
                if (normP(p - arc.points.back()) > tiny) arc.points.push_back(p);
            }
            // The two ends are the nodes, exactly. A curve that Stage 7 snapped
            // to a cone stops within the snap tolerance of it, not on it, and
            // the layout wants the corner on the cone.
            arc.points.front() = nodes[arc.from].p;
            arc.points.back() = nodes[arc.to].p;

            arc.faceFrom = arc.steps.front().face;
            arc.faceTo = arc.steps.back().face;
            arc.dirFrom = cu.dir;
            arc.dirTo = cu.endDir;
            arc.family = (cu.dir % 2 == 0) ? Align::V : Align::U;
            arc.familyTurns = (cu.dir % 2) != (cu.endDir % 2);
            arc.dangling = nodes[arc.to].kind == NodeKind::Dangling ||
                           nodes[arc.from].kind == NodeKind::Dangling;
            for (size_t i = 1; i < arc.points.size(); ++i) {
                arc.modelLength += normP(arc.points[i] - arc.points[i - 1]);
            }
            for (const Separatrices::Step &st : arc.steps) {
                arc.imageLength += normP(sep->point(st, true, Separatrices::Space::Image) -
                                         sep->point(st, false, Separatrices::Space::Image));
            }
            arc.degenerate = arc.modelLength <= mergeTol;
            arcs.push_back(std::move(arc));
        }
    }
}

// ---------------------------------------------------------------------------
// buildBoundaryArcs()
//
// dS is part of the layout: the sides of the patches that touch the boundary
// are pieces of it. It is walked in the direction that keeps the model on the
// left and cut at every node that landed on it -- a cone, a separatrix hitting
// it, a change of label -- so that a boundary arc is a side and a boundary
// half-edge running backwards is the outside.
// ---------------------------------------------------------------------------
void Arrangement::buildBoundaryArcs() {
    if (dirBoundary.empty()) return;

    std::vector<char> done(dirBoundary.size(), 0);
    for (size_t seed = 0; seed < dirBoundary.size(); ++seed) {
        if (done[seed]) continue;

        // The loop this edge belongs to, in order.
        std::vector<int> loop;
        int cur = static_cast<int>(seed);
        while (!done[cur]) {
            done[cur] = 1;
            loop.push_back(cur);
            const int nx = boundaryOut[dirBoundary[cur].b];
            if (nx < 0) break;
            cur = nx;
        }
        // Nodes on the loop, in order, as (position in loop, parameter).
        std::vector<std::pair<double, int>> stops;
        for (size_t i = 0; i < loop.size(); ++i) {
            const DirEdge &de = dirBoundary[loop[i]];
            if (vertexNode[de.a] >= 0) {
                stops.push_back({static_cast<double>(i), vertexNode[de.a]});
            }
            std::vector<Event> &ev = boundaryEvents[loop[i]];
            std::sort(ev.begin(), ev.end(),
                      [](const Event &a, const Event &b) { return a.param < b.param; });
            for (const Event &e : ev) stops.push_back({i + e.param, e.node});
        }
        if (stops.empty()) {
            // A boundary loop with nothing on it bounds an annulus, not a
            // quadrilateral. Anchor it so the walk terminates and say so.
            const int v = dirBoundary[loop.front()].a;
            const int id = addNode(orig->vertices[v], NodeKind::BoundaryCorner);
            nodes[id].vertex = v;
            nodes[id].onBoundary = true;
            vertexNode[v] = id;
            stops.push_back({0.0, id});
            ++report.annularFaces;
        }

        const double total = static_cast<double>(loop.size());
        for (size_t k = 0; k < stops.size(); ++k) {
            const double p0 = stops[k].first;
            const double p1 = (k + 1 < stops.size()) ? stops[k + 1].first
                                                     : stops.front().first + total;
            Arc arc;
            arc.kind = ArcKind::Boundary;
            arc.from = stops[k].second;
            arc.to = (k + 1 < stops.size()) ? stops[k + 1].second : stops.front().second;

            arc.points.push_back(nodes[arc.from].p);
            for (double q = std::floor(p0) + 1.0; q < p1 - 1e-12; q += 1.0) {
                const int idx = loop[static_cast<int>(q) % loop.size()];
                arc.verts.push_back(dirBoundary[idx].a);
                arc.points.push_back(orig->vertices[dirBoundary[idx].a]);
            }
            arc.points.push_back(nodes[arc.to].p);
            // Which of the two coordinates dS holds constant here, straight
            // from Stage 5 rather than re-derived: the labels are chain-wide
            // and an edge-by-edge test can disagree with them on a short edge.
            const int mid = loop[static_cast<int>(std::floor((p0 + p1) * 0.5)) % loop.size()];
            auto it = edgeLabel.find(MeshEdgeKey(dirBoundary[mid].a, dirBoundary[mid].b));
            arc.family = (it == edgeLabel.end()) ? Align::None
                                                 : static_cast<Align>(it->second);
            arc.bdFrom = loop[static_cast<int>(std::floor(p0)) % loop.size()];
            arc.bdTo = loop[(static_cast<int>(std::ceil(p1)) - 1 + 2 * static_cast<int>(loop.size()))
                            % loop.size()];
            arc.faceFrom = dirBoundary[arc.bdFrom].face;
            arc.faceTo = dirBoundary[arc.bdTo].face;
            for (size_t i = 1; i < arc.points.size(); ++i) {
                arc.modelLength += normP(arc.points[i] - arc.points[i - 1]);
            }
            arc.degenerate = arc.modelLength <= mergeTol;
            arcs.push_back(std::move(arc));
        }
    }
}

// ---------------------------------------------------------------------------
// buildInterfaceNodes()
//
// The nodes of the material interface network are nodes of the layout. Stage 0b
// already decided which vertices they are and how many quarter turns the layout
// makes in each sector at them, and Stage 6 has just spent its whole
// continuation making Psi agree; all that is left here is to put them in.
//
// A node that already carries one -- a cone the field put on an interface, or a
// landing, which is on dS and may already be a boundary corner -- keeps the one
// it has. The kinds are not interchangeable downstream: a cone's valence is
// checked against its index and an interface node's is not.
// ---------------------------------------------------------------------------
void Arrangement::buildInterfaceNodes() {
    itf = lab->getInterfaces();
    if (!itf || !itf->multiMaterial()) { itf = nullptr; return; }

    for (const Interfaces::Node &nd : itf->nodes()) {
        const int v = nd.vertex;
        if (v < 0 || v >= static_cast<int>(vertexNode.size())) continue;
        if (vertexNode[v] >= 0) continue;   // already a cone, or a corner of dS
        const int id = addNode(orig->vertices[v], NodeKind::InterfaceNode);
        nodes[id].vertex = v;
        nodes[id].onBoundary = orig->isBoundaryVertex[v];
        vertexNode[v] = id;
    }

    // Every interface edge, directed the way its branch runs, so that a branch
    // is a contiguous run of dirInterface and an event on it has a parameter
    // that is an index into that run plus a fraction.
    const std::vector<Interfaces::Branch> &brs = itf->branches();
    branchStart.assign(brs.size(), 0);
    branchCount.assign(brs.size(), 0);
    for (size_t b = 0; b < brs.size(); ++b) {
        branchStart[b] = static_cast<int>(dirInterface.size());
        const Interfaces::Branch &br = brs[b];
        for (size_t i = 0; i + 1 < br.verts.size() && i < br.edges.size(); ++i) {
            DirEdge de;
            de.a = br.verts[i];
            de.b = br.verts[i + 1];
            de.edge = br.edges[i];
            de.face = orig->edgeTriangles[de.edge][0] >= 0 ? orig->edgeTriangles[de.edge][0]
                                                           : orig->edgeTriangles[de.edge][1];
            dirInterface.push_back(de);
        }
        branchCount[b] = static_cast<int>(dirInterface.size()) - branchStart[b];
    }
    interfaceEvents.assign(dirInterface.size(), {});
}

// ---------------------------------------------------------------------------
// findInterfaceCrossings()
//
// Where a separatrix meets an interface. Same question as findCrossings(), and
// asked the same way -- triangle-locally, on the model -- but between a step of
// a curve and an edge of the triangulation rather than between two steps.
//
// A separatrix crosses an interface *at* an edge of the mesh, so the crossing
// is the endpoint of one step and the start of the next and is found twice, once
// from each side. That is exactly the case the node merge tolerance exists for,
// and the event dedupe in splitCurves() takes care of the curve's own copy.
// ---------------------------------------------------------------------------
void Arrangement::findInterfaceCrossings() {
    if (!itf) return;

    // Which entries of dirInterface each face of S carries.
    std::vector<std::vector<int>> faceInterface(orig->triangles.size());
    for (size_t i = 0; i < dirInterface.size(); ++i) {
        for (int k = 0; k < 2; ++k) {
            const int f = orig->edgeTriangles[dirInterface[i].edge][k];
            if (f >= 0) faceInterface[f].push_back(static_cast<int>(i));
        }
    }

    for (size_t f = 0; f < faceSegs.size(); ++f) {
        if (faceInterface[f].empty()) continue;
        for (int si : faceSegs[f]) {
            const Seg &p = segs[si];
            for (int ii : faceInterface[f]) {
                const DirEdge &de = dirInterface[ii];
                const Point qa = orig->vertices[de.a];
                const Point qb = orig->vertices[de.b];

                const Point r = p.b - p.a;
                const Point sdir = qb - qa;
                const double den = cross2(r, sdir);
                const double scale = normP(r) * normP(sdir);
                if (std::fabs(den) <= 1e-12 * scale) continue;   // running along it
                const Point w = qa - p.a;
                const double tp = cross2(w, sdir) / den;
                const double tq = cross2(w, r) / den;
                const double ep = mergeTol / std::max(normP(r), 1e-300);
                const double eq = mergeTol / std::max(normP(sdir), 1e-300);
                if (tp < -ep || tp > 1.0 + ep) continue;
                if (tq < -eq || tq > 1.0 + eq) continue;

                const double cq = std::max(0.0, std::min(1.0, tq));
                const Point hit = qa + sdir * cq;
                const int id = findOrAddCrossing(hit, static_cast<int>(f));
                if (nodes[id].kind == NodeKind::Crossing) nodes[id].kind = NodeKind::InterfaceHit;
                curveEvents[p.curve].push_back(
                    Event{p.step + std::max(0.0, std::min(1.0, tp)), id});
                interfaceEvents[ii].push_back(Event{cq, id});
            }
        }
    }
}

// ---------------------------------------------------------------------------
// buildInterfaceArcs()
//
// Each branch of the network, cut at every node on it. Structurally the same
// walk as buildBoundaryArcs(), and different in one way that matters: a branch
// is a path and not a cycle, so its two ends are always nodes -- Stage 0b made
// them so -- and there is no anchoring case to handle.
// ---------------------------------------------------------------------------
void Arrangement::buildInterfaceArcs() {
    if (!itf || dirInterface.empty()) return;
    const std::vector<Interfaces::Branch> &brs = itf->branches();

    for (size_t b = 0; b < brs.size(); ++b) {
        const int base = branchStart[b], n = branchCount[b];
        if (n <= 0) continue;

        std::vector<std::pair<double, int>> stops;
        for (int i = 0; i < n; ++i) {
            const DirEdge &de = dirInterface[base + i];
            if (vertexNode[de.a] >= 0) stops.push_back({static_cast<double>(i), vertexNode[de.a]});
            std::vector<Event> &ev = interfaceEvents[base + i];
            std::sort(ev.begin(), ev.end(),
                      [](const Event &x, const Event &y) { return x.param < y.param; });
            for (const Event &e : ev) stops.push_back({i + e.param, e.node});
        }
        if (vertexNode[dirInterface[base + n - 1].b] >= 0) {
            stops.push_back({static_cast<double>(n), vertexNode[dirInterface[base + n - 1].b]});
        }
        // Two stops at the same point -- a crossing that landed on a node -- are
        // one stop, and a branch with fewer than two is not an arc.
        std::sort(stops.begin(), stops.end());
        std::vector<std::pair<double, int>> keep;
        for (const auto &st : stops) {
            if (!keep.empty() && (st.second == keep.back().second ||
                                  st.first - keep.back().first < 1e-9)) {
                continue;
            }
            keep.push_back(st);
        }
        if (keep.size() < 2) continue;

        for (size_t k = 0; k + 1 < keep.size(); ++k) {
            const double p0 = keep[k].first, p1 = keep[k + 1].first;
            Arc arc;
            arc.kind = ArcKind::Interface;
            arc.branch = static_cast<int>(b);
            arc.from = keep[k].second;
            arc.to = keep[k + 1].second;
            arc.matLeft = brs[b].matLeft;
            arc.matRight = brs[b].matRight;

            arc.points.push_back(nodes[arc.from].p);
            for (double q = std::floor(p0) + 1.0; q < p1 - 1e-12; q += 1.0) {
                const int idx = base + static_cast<int>(q);
                if (idx < base || idx >= base + n) continue;
                arc.verts.push_back(dirInterface[idx].a);
                arc.points.push_back(orig->vertices[dirInterface[idx].a]);
            }
            arc.points.push_back(nodes[arc.to].p);

            arc.ifFrom = base + std::min(n - 1, static_cast<int>(std::floor(p0)));
            arc.ifTo = base + std::max(0, std::min(n - 1, static_cast<int>(std::ceil(p1)) - 1));
            arc.faceFrom = dirInterface[arc.ifFrom].face;
            arc.faceTo = dirInterface[arc.ifTo].face;

            // Which coordinate it holds constant, straight from Stage 5's
            // propagated direction rather than re-read off the geometry.
            arc.family = Align::None;
            for (const auto &fc : lab->featureChains()) {
                if (fc.branch == static_cast<int>(b)) { arc.family = fc.label; break; }
            }
            for (size_t i = 1; i < arc.points.size(); ++i) {
                arc.modelLength += normP(arc.points[i] - arc.points[i - 1]);
            }
            arc.degenerate = arc.modelLength <= mergeTol;
            arcs.push_back(std::move(arc));
        }
    }
}

// ---------------------------------------------------------------------------
// classifyMaterials()
//
// Which material each patch is in, and whether it is all one.
//
// The point of the whole multi-material extension is that a patch lies inside
// one material, so this is the property worth measuring rather than assuming.
// It is measured by locating points of the patch in the triangulation: the
// midpoint of each of its arcs, stepped a short way into the patch along the
// inward normal, which lands inside the patch for any arc that is not itself
// degenerate. A patch whose samples disagree is a patch an element of the final
// mesh could straddle, and it is counted.
// ---------------------------------------------------------------------------
void Arrangement::classifyMaterials() {
    if (orig->triangleMatId.size() != orig->triangles.size()) return;
    const double step = 1e-3 * modelExtent;

    for (int fi : patches) {
        Face &fc = faces[fi];
        std::vector<int> found;
        for (int h : fc.half) {
            const Arc &ar = arcs[halves[h].arc];
            if (ar.points.size() < 2 || ar.modelLength <= mergeTol) continue;
            const size_t m = ar.points.size() / 2;
            const size_t m0 = (m == 0) ? 0 : m - 1;
            const Point a = ar.points[m0], b = ar.points[m0 + 1];
            const Point d = normalizeP(b - a);
            if (!(normP(d) > 0.0)) continue;
            // The face lies to the left of its own half-edges, so the inward
            // normal is the left normal of the direction of travel.
            const double sgn = halves[h].forward ? 1.0 : -1.0;
            const Point mid = (a + b) * 0.5;
            const Point inward{-d[1] * sgn, d[0] * sgn};
            const int t = orig->findTriangleContainingPoint(mid + inward * step);
            if (t >= 0) found.push_back(orig->triangleMatId[t]);
        }
        if (found.empty()) continue;
        fc.material = found.front();
        for (int m : found) if (m != fc.material) fc.mixed = true;
        if (fc.mixed) ++report.mixedPatches;
    }
}

// ---------------------------------------------------------------------------
// collapseShortArcs()
//
// A sliver arc is not a layout edge and its two ends are one node; see
// Options::collapseTolerance for where they come from. Collapsing is done here
// rather than by widening the node merge tolerance because the thing that is
// too small is an *arc*, and an arc's length is the quantity worth thresholding
// -- two nodes can be far apart in the plane and still be joined by nothing.
//
// The node that survives is the one that means the most: a cone before a
// boundary node before a crossing. Two cones a sliver apart are never merged --
// that is Sec. 3.1's clustering failure, whose remedy is upstream (merge the
// cluster into one cone of the summed index) and cannot be applied once the
// metric has been built around both of them -- so it is counted and left alone.
// ---------------------------------------------------------------------------
void Arrangement::collapseShortArcs() {
    const double tol = options.collapseTolerance * modelExtent;
    if (!(tol > 0.0) || arcs.empty()) return;

    std::vector<int> alias(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i) alias[i] = static_cast<int>(i);
    std::function<int(int)> find = [&](int x) {
        while (alias[x] != x) { alias[x] = alias[alias[x]]; x = alias[x]; }
        return x;
    };
    auto rank = [&](int n) {
        switch (nodes[n].kind) {
            case NodeKind::Cone: return 4;
            case NodeKind::BoundaryCorner: return 3;
            case NodeKind::BoundaryHit: return 2;
            case NodeKind::Crossing: return 1;
            default: return 0;
        }
    };

    std::vector<int> order(arcs.size());
    for (size_t i = 0; i < arcs.size(); ++i) order[i] = static_cast<int>(i);
    std::sort(order.begin(), order.end(),
              [this](int a, int b) { return arcs[a].modelLength < arcs[b].modelLength; });

    std::vector<char> dead(arcs.size(), 0);
    for (int a : order) {
        if (arcs[a].modelLength > tol) break;
        const int u = find(arcs[a].from), v = find(arcs[a].to);
        if (u == v) { dead[a] = 1; ++report.collapsedArcs; continue; }
        if (nodes[u].kind == NodeKind::Cone && nodes[v].kind == NodeKind::Cone) {
            ++report.clusteredCones;
            continue;
        }
        // Only where the merge cannot change what the model is. Two nodes of dS
        // joined by a piece of dS are neighbours along it and merging them
        // shortens the boundary; two interior nodes joined by a sliver are one
        // point inside a face. Anything else pinches: a short separatrix
        // between two points of dS closes the model off across itself, and a
        // sliver from a crossing onto dS drags an interior node onto a boundary
        // whose fan is a half plane its other arcs do not lie in. Both leave a
        // face that is inside the model along part of its boundary and outside
        // along the rest, which is not a face.
        const bool bothOnBoundary = nodes[u].onBoundary && nodes[v].onBoundary;
        const bool bothInside = !nodes[u].onBoundary && !nodes[v].onBoundary;
        const bool safe = (arcs[a].kind == ArcKind::Boundary && bothOnBoundary) ||
                          (arcs[a].kind == ArcKind::Separatrix && bothInside);
        if (!safe) { ++report.sliversKept; continue; }
        const int keep = rank(u) >= rank(v) ? u : v;
        alias[keep == u ? v : u] = keep;
        dead[a] = 1;
        ++report.collapsedArcs;
    }
    if (report.collapsedArcs == 0) return;

    // Compact: the surviving nodes, then the arcs re-anchored onto them.
    std::vector<int> newId(nodes.size(), -1);
    std::vector<Node> keptNodes;
    std::vector<Point> keptRef;
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (find(static_cast<int>(i)) != static_cast<int>(i)) continue;
        newId[i] = static_cast<int>(keptNodes.size());
        keptNodes.push_back(nodes[i]);
        keptRef.push_back(nodeRef[i]);
    }
    for (size_t i = 0; i < nodes.size(); ++i) newId[i] = newId[find(static_cast<int>(i))];

    std::vector<Arc> keptArcs;
    for (size_t a = 0; a < arcs.size(); ++a) {
        if (dead[a]) continue;
        Arc ar = std::move(arcs[a]);
        ar.from = newId[ar.from];
        ar.to = newId[ar.to];
        if (ar.points.size() >= 2) {
            ar.points.front() = keptNodes[ar.from].p;
            ar.points.back() = keptNodes[ar.to].p;
            ar.modelLength = 0.0;
            for (size_t k = 1; k < ar.points.size(); ++k) {
                ar.modelLength += normP(ar.points[k] - ar.points[k - 1]);
            }
        }
        ar.degenerate = ar.modelLength <= mergeTol;
        keptArcs.push_back(std::move(ar));
    }
    nodes.swap(keptNodes);
    nodeRef.swap(keptRef);
    arcs.swap(keptArcs);
}

// ---------------------------------------------------------------------------
Point Arrangement::dirEdgeImage(int be) const {
    if (be < 0 || be >= static_cast<int>(dirBoundary.size())) return Point{1.0, 0.0};
    const DirEdge &de = dirBoundary[be];
    const Triangle &t = orig->triangles[de.face];
    const Triangle &ct = cut->triangles[de.face];
    int la = -1, lb = -1;
    for (int k = 0; k < 3; ++k) {
        if (t[k] == de.a) la = k;
        if (t[k] == de.b) lb = k;
    }
    if (la < 0 || lb < 0) return Point{1.0, 0.0};
    return (*uv)[ct[lb]] - (*uv)[ct[la]];
}

// ---------------------------------------------------------------------------
Point Arrangement::dirInterfaceImage(int ie) const {
    if (ie < 0 || ie >= static_cast<int>(dirInterface.size())) return Point{1.0, 0.0};
    const DirEdge &de = dirInterface[ie];
    if (de.face < 0) return Point{1.0, 0.0};
    const Triangle &t = orig->triangles[de.face];
    const Triangle &ct = cut->triangles[de.face];
    int la = -1, lb = -1;
    for (int k = 0; k < 3; ++k) {
        if (t[k] == de.a) la = k;
        if (t[k] == de.b) lb = k;
    }
    if (la < 0 || lb < 0) return Point{1.0, 0.0};
    return (*uv)[ct[lb]] - (*uv)[ct[la]];
}

// ---------------------------------------------------------------------------
// layoutAngle()
//
// Where an arc leaves a node, measured in Psi rather than on S. Away from the
// cones this is one atan2: Psi is a single global map on Omega, so two arcs in
// two different faces are already in the same frame. At a cone it is not, and
// the reason is the one Stage 7's fan sweep has: the cone angle is
// 2pi - (pi/2) I, which for a negative index is more than a full turn, so two
// arcs can leave along the same image direction and only the angle accumulated
// around the fan separates them. Accumulating also carries the sweep across an
// arc of the cutting graph without any transition being applied, because what
// is summed is the angle *inside* each face and never an absolute direction.
// ---------------------------------------------------------------------------
double Arrangement::layoutAngle(int n, int f, const Point &imageDir, const Point &modelDir) {
    Node &nd = nodes[n];

    if (nd.vertex >= 0) {
        buildFan(nd.vertex);
        const std::vector<int> &fan = fanCache[nd.vertex];
        nd.totalAngle = fanTotal[nd.vertex];
        if (fan.empty() || !(nd.totalAngle > 0.0)) return 0.0;

        // Which face's sector the direction lies in. The arc's own face is
        // where the search starts and is usually the answer, but it is not
        // always: a curve that Stage 7 snapped onto this cone finished in
        // whichever face it happened to be crossing, which can be a face or two
        // round the fan from the one the arc actually leaves the cone through.
        // Trusting it costs a whole face's angle, which is enough to turn a
        // corner into a T-junction, so the sector is *found* rather than
        // assumed -- outwards from the arc's own face, because at a cone of
        // negative index the fan sweeps more than a full turn and the same
        // direction lies in two of its faces.
        const int m = static_cast<int>(fan.size());
        int pos = -1;
        for (int i = 0; i < m; ++i) if (fan[i] == f) { pos = i; break; }
        const bool closed = !orig->isBoundaryVertex[nd.vertex];

        // The angle of `imageDir` inside face fan[i]'s sector, or nothing if it
        // does not lie in it. A direction that comes back a hair under a full
        // turn is a direction along the sector's own start edge, rounded the
        // wrong side of it: with contracted multiply-add even the cross product
        // of a vector with itself is not exactly zero, so this case is the rule
        // and not the exception on the two edges of dS at a boundary cone.
        const double slack = 1e-9;
        auto angleIn = [&](int i, double &a, double &w) {
            const int g = fan[i];
            int iv, ib, ic;
            if (!cornerOf(orig->triangles[g], ccw(g), nd.vertex, iv, ib, ic)) return false;
            const Triangle &ct = cut->triangles[g];
            const Point A = (*uv)[ct[ib]] - (*uv)[ct[iv]];
            a = std::atan2(cross2(A, imageDir), dotP(A, imageDir));
            if (a < 0.0) a += 2.0 * M_PI;
            if (a >= 2.0 * M_PI - slack) a = 0.0;
            w = imageAngle(g, nd.vertex);
            return true;
        };
        auto within = [&](int i, double &acc) {
            double a = 0.0, w = 0.0;
            if (!angleIn(i, a, w)) return false;
            if (a > w + 1e-12) return false;
            acc = fanOffset[nd.vertex][i] + std::min(a, w);
            return true;
        };

        if (pos >= 0) {
            for (int d = 0; d < m; ++d) {
                for (int sgn = 1; sgn >= -1; sgn -= 2) {
                    if (d == 0 && sgn < 0) continue;
                    int i = pos + sgn * d;
                    if (closed) i = ((i % m) + m) % m;
                    else if (i < 0 || i >= m) continue;
                    double acc = 0.0;
                    if (within(i, acc)) return acc;
                }
            }
            // No sector claims it, which takes a direction that is not a
            // direction. Put it at the near end of the arc's own face.
            double a = 0.0, w = 0.0;
            if (angleIn(pos, a, w)) {
                return fanOffset[nd.vertex][pos] + ((a - w <= 2.0 * M_PI - a) ? w : 0.0);
            }
        }

        // The arc is not in this vertex's fan at all, which happens only when
        // Stage 7 snapped a curve to a cone from beyond its one-ring (see
        // Separatrices::Options::coneSnapRings). Place it by the direction it
        // arrives from on S, spread across the fan in proportion to the angle
        // each face carries in the image. Monotone and continuous, which is all
        // the sector test below needs of it.
        ++report.interpolatedAngles;
        double best = 0.0, bestGap = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < fan.size(); ++i) {
            const int g = fan[i];
            int iv, ib, ic;
            if (!cornerOf(orig->triangles[g], ccw(g), nd.vertex, iv, ib, ic)) continue;
            const Triangle &t = orig->triangles[g];
            const Point A = orig->vertices[t[ib]] - orig->vertices[t[iv]];
            const Point B = orig->vertices[t[ic]] - orig->vertices[t[iv]];
            double wm = std::atan2(cross2(A, B), dotP(A, B));
            if (wm < 0.0) wm += 2.0 * M_PI;
            double am = std::atan2(cross2(A, modelDir), dotP(A, modelDir));
            if (am < 0.0) am += 2.0 * M_PI;
            const double gap = (am <= wm) ? 0.0 : std::min(am - wm, 2.0 * M_PI - am);
            if (gap < bestGap) {
                bestGap = gap;
                const double frac = (wm > 0.0) ? std::min(1.0, am / wm) : 0.0;
                best = fanOffset[nd.vertex][i] + imageAngle(g, nd.vertex) * frac;
            }
        }
        return best;
    }

    if (nd.onBoundary) {
        nd.totalAngle = M_PI;
        const Point &r = nodeRef[n];
        double a = std::atan2(cross2(r, imageDir), dotP(r, imageDir));
        if (a < 0.0) a += 2.0 * M_PI;
        // The same question on the half-plane at a boundary node: outside
        // [0, pi] means along dS, and which way along it is nearest decides.
        if (a > M_PI) a = (a - M_PI <= 2.0 * M_PI - a) ? M_PI : 0.0;
        return a;
    }

    nd.totalAngle = 2.0 * M_PI;
    double a = std::atan2(imageDir[1], imageDir[0]);
    if (a < 0.0) a += 2.0 * M_PI;
    if (a > 2.0 * M_PI - 1e-9) a = 0.0;
    return a;
}

// ---------------------------------------------------------------------------
void Arrangement::buildAngles() {
    const double tiny = 1e-14 * modelExtent;
    halves.assign(2 * arcs.size(), HalfEdge());

    for (size_t a = 0; a < arcs.size(); ++a) {
        Arc &ar = arcs[a];
        HalfEdge &h0 = halves[2 * a];
        HalfEdge &h1 = halves[2 * a + 1];
        h0.arc = h1.arc = static_cast<int>(a);
        h0.forward = true;  h0.origin = ar.from; h0.twin = static_cast<int>(2 * a + 1);
        h1.forward = false; h1.origin = ar.to;   h1.twin = static_cast<int>(2 * a);

        Point iFrom{1.0, 0.0}, iTo{1.0, 0.0};
        if (ar.kind == ArcKind::Separatrix) {
            for (const Separatrices::Step &st : ar.steps) {
                const Point d = sep->point(st, true, Separatrices::Space::Image) -
                                sep->point(st, false, Separatrices::Space::Image);
                if (normP(d) > tiny) { iFrom = d; ar.faceFrom = st.face; break; }
            }
            for (size_t i = ar.steps.size(); i-- > 0;) {
                const Separatrices::Step &st = ar.steps[i];
                const Point d = sep->point(st, false, Separatrices::Space::Image) -
                                sep->point(st, true, Separatrices::Space::Image);
                if (normP(d) > tiny) { iTo = d; ar.faceTo = st.face; break; }
            }
        } else if (ar.kind == ArcKind::Interface) {
            iFrom = dirInterfaceImage(ar.ifFrom);
            iTo = dirInterfaceImage(ar.ifTo) * -1.0;
        } else {
            iFrom = dirEdgeImage(ar.bdFrom);
            iTo = dirEdgeImage(ar.bdTo) * -1.0;
        }

        const size_t np = ar.points.size();
        const Point mFrom = np > 1 ? ar.points[1] - ar.points[0] : Point{1.0, 0.0};
        const Point mTo = np > 1 ? ar.points[np - 2] - ar.points[np - 1] : Point{1.0, 0.0};

        h0.angle = layoutAngle(ar.from, ar.faceFrom, iFrom, mFrom);
        h1.angle = layoutAngle(ar.to, ar.faceTo, iTo, mTo);
    }

    for (size_t h = 0; h < halves.size(); ++h) {
        const int o = halves[h].origin;
        if (o >= 0) nodes[o].out.push_back(static_cast<int>(h));
    }
    for (Node &n : nodes) {
        std::sort(n.out.begin(), n.out.end(),
                  [this](int a, int b) { return halves[a].angle < halves[b].angle; });
        n.angle.clear();
        for (size_t i = 0; i < n.out.size(); ++i) {
            halves[n.out[i]].slot = static_cast<int>(i);
            n.angle.push_back(halves[n.out[i]].angle);
        }
        if (!(n.totalAngle > 0.0)) n.totalAngle = n.onBoundary ? M_PI : 2.0 * M_PI;
    }
}

// ---------------------------------------------------------------------------
// linkHalfEdges()
//
// next(h) is the cyclic *predecessor* of twin(h) in the sorted fan at the node
// h ends on. That is the standard rule and it is what puts the interior of the
// face on the left of every half-edge, so that a bounded face comes out with
// positive signed area and the unbounded one with negative.
// ---------------------------------------------------------------------------
void Arrangement::linkHalfEdges() {
    for (HalfEdge &h : halves) {
        const HalfEdge &tw = halves[h.twin];
        const Node &b = nodes[tw.origin];
        const int deg = static_cast<int>(b.out.size());
        if (deg <= 0) { h.next = h.twin; continue; }
        h.next = b.out[(tw.slot - 1 + deg) % deg];
    }
}

void Arrangement::extractFaces() {
    std::vector<int> owner(halves.size(), -1);
    for (size_t h = 0; h < halves.size(); ++h) {
        if (owner[h] >= 0) continue;
        Face fc;
        int cur = static_cast<int>(h);
        for (size_t guard = 0; guard <= halves.size(); ++guard) {
            owner[cur] = static_cast<int>(faces.size());
            halves[cur].face = static_cast<int>(faces.size());
            fc.half.push_back(cur);
            cur = halves[cur].next;
            if (cur == static_cast<int>(h)) break;
            if (owner[cur] >= 0) break;   // should not happen; guards a bad link
        }
        faces.push_back(std::move(fc));
    }
}

// ---------------------------------------------------------------------------
// classifyFaces()
//
// Sec. 4's validation, face by face. The turn at a node is the sector the face
// occupies there, which is a difference of two of the layout angles already
// sorted at that node, so it comes out for free and in the right units: a
// corner is one quarter turn, a T-junction is none, and anything else is a
// separatrix that is missing or went to the wrong cone.
// ---------------------------------------------------------------------------
void Arrangement::classifyFaces() {
    const double quarter = M_PI_2;

    for (Face &fc : faces) {
        const size_t m = fc.half.size();
        if (m == 0) continue;

        std::vector<Point> loop;
        for (size_t i = 0; i < m; ++i) {
            const HalfEdge &h = halves[fc.half[i]];
            const Arc &ar = arcs[h.arc];
            if (h.forward) {
                for (size_t k = 0; k + 1 < ar.points.size(); ++k) loop.push_back(ar.points[k]);
            } else {
                for (size_t k = ar.points.size(); k-- > 1;) loop.push_back(ar.points[k]);
            }
        }
        fc.area = loop.size() >= 3 ? signedArea(loop) : 0.0;

        bool outside = false, inside = false;
        for (size_t i = 0; i < m; ++i) {
            const HalfEdge &h = halves[fc.half[i]];
            if (arcs[h.arc].kind != ArcKind::Boundary) continue;
            if (h.forward) inside = true; else outside = true;
        }
        if (outside && inside) ++report.mixedBoundaryFaces;
        fc.patch = !outside && fc.area > 0.0;

        // The turn at the end of each half-edge.
        fc.turns.assign(m, 0);
        std::vector<char> corner(m, 0);
        for (size_t i = 0; i < m; ++i) {
            const HalfEdge &h = halves[fc.half[i]];
            const HalfEdge &g = halves[h.next];
            const Node &nd = nodes[g.origin];
            double s = halves[h.twin].angle - g.angle;
            if (s <= 1e-9) s += nd.totalAngle;
            s = std::max(0.0, std::min(nd.totalAngle, s));
            const int q = static_cast<int>(std::lround(s / quarter));
            fc.turns[i] = q;
            corner[i] = (q == 1);
            if (fc.patch) {
                const double res = std::fabs(s - q * quarter);
                if (res > report.maxSectorResidual) {
                    report.maxSectorResidual = res;
                    report.worstSectorNode = g.origin;
                }
                if (res > options.cornerTolerance) ++report.ambiguousSectors;
                if (q == 0) ++fc.tJunctions;
            }
        }

        for (size_t i = 0; i < m; ++i) {
            if (corner[i]) fc.corners.push_back(halves[halves[fc.half[i]].next].origin);
        }
        fc.quad = (fc.corners.size() == 4);

        // Sides: the runs of half-edges between one corner and the next.
        int first = -1;
        for (size_t i = 0; i < m; ++i) if (corner[i]) { first = static_cast<int>(i); break; }
        if (first >= 0) {
            std::vector<int> side;
            for (size_t k = 1; k <= m; ++k) {
                const size_t i = (first + k) % m;
                side.push_back(fc.half[i]);
                if (corner[i]) { fc.sides.push_back(side); side.clear(); }
            }
            // The corner list has to start where the sides do.
            fc.corners.clear();
            for (const std::vector<int> &sd : fc.sides) {
                fc.corners.push_back(halves[halves[sd.back()].next].origin);
            }
            std::rotate(fc.corners.begin(), fc.corners.end() - 1, fc.corners.end());
        }
        fc.simple = fc.quad && fc.sides.size() == 4 &&
                    fc.sides[0].size() == 1 && fc.sides[1].size() == 1 &&
                    fc.sides[2].size() == 1 && fc.sides[3].size() == 1;
    }

    report.minPatchArea = std::numeric_limits<double>::infinity();
    for (size_t f = 0; f < faces.size(); ++f) {
        const Face &fc = faces[f];
        if (!fc.patch) continue;
        patches.push_back(static_cast<int>(f));
        report.patchArea += fc.area;
        report.minPatchArea = std::min(report.minPatchArea, fc.area);
        report.maxPatchArea = std::max(report.maxPatchArea, fc.area);
        report.tJunctions += fc.tJunctions;
        if (!fc.quad) ++report.wrongCornerFaces;
    }
    if (!std::isfinite(report.minPatchArea)) report.minPatchArea = 0.0;
}

// ---------------------------------------------------------------------------
// checkFeatures()
//
// "Every feature chain must be a union of arcs." E3 of Stage 6 aligned the
// features to coordinate lines so that the layout could contain them; this is
// where that promise is collected. A chain that runs along dS is covered by the
// boundary arcs and is checked the same way as any other, which is the point --
// the check does not care why an arc is there.
// ---------------------------------------------------------------------------
void Arrangement::checkFeatures() {
    const std::vector<SubdomainLabels::FeatureChain> &chains = lab->featureChains();
    report.featureChains = static_cast<int>(chains.size());
    if (chains.empty()) return;

    const double tol = options.featureTolerance * modelExtent;
    const double cell = std::max(tol, modelExtent / 512.0);
    auto key = [&](double x, double y) {
        return static_cast<long long>(std::floor(x / cell)) * 73856093LL ^
               static_cast<long long>(std::floor(y / cell)) * 19349663LL;
    };
    std::unordered_map<long long, std::vector<std::array<int, 2>>> grid;
    for (size_t a = 0; a < arcs.size(); ++a) {
        const Arc &ar = arcs[a];
        for (size_t k = 0; k + 1 < ar.points.size(); ++k) {
            const Point &p = ar.points[k];
            const Point &q = ar.points[k + 1];
            const int x0 = static_cast<int>(std::floor(std::min(p[0], q[0]) / cell));
            const int x1 = static_cast<int>(std::floor(std::max(p[0], q[0]) / cell));
            const int y0 = static_cast<int>(std::floor(std::min(p[1], q[1]) / cell));
            const int y1 = static_cast<int>(std::floor(std::max(p[1], q[1]) / cell));
            for (int x = x0; x <= x1; ++x) {
                for (int y = y0; y <= y1; ++y) {
                    grid[static_cast<long long>(x) * 73856093LL ^
                         static_cast<long long>(y) * 19349663LL]
                        .push_back({{static_cast<int>(a), static_cast<int>(k)}});
                }
            }
        }
    }

    const std::vector<int> &toOrig = imm->getCut().getCutVertexToOriginal();
    for (const SubdomainLabels::FeatureChain &ch : chains) {
        double worst = 0.0;
        for (size_t i = 0; i + 1 < ch.verts.size(); ++i) {
            const int a = toOrig[ch.verts[i]], b = toOrig[ch.verts[i + 1]];
            if (a < 0 || b < 0) continue;
            const Point mid = (orig->vertices[a] + orig->vertices[b]) * 0.5;
            double best = std::numeric_limits<double>::infinity();
            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    auto it = grid.find(key(mid[0] + dx * cell, mid[1] + dy * cell));
                    if (it == grid.end()) continue;
                    for (const auto &sg : it->second) {
                        const Arc &ar = arcs[sg[0]];
                        best = std::min(best, pointSegment(mid, ar.points[sg[1]],
                                                           ar.points[sg[1] + 1]));
                    }
                }
            }
            if (!std::isfinite(best)) {
                // Nothing in the neighbouring cells: the chain is not merely
                // off an arc, it is nowhere near one. Say how far, rather than
                // leaving the report with an infinity or a misleading zero.
                for (const Arc &ar : arcs) {
                    for (size_t q = 0; q + 1 < ar.points.size(); ++q) {
                        best = std::min(best, pointSegment(mid, ar.points[q], ar.points[q + 1]));
                    }
                }
            }
            worst = std::max(worst, best);
        }
        if (std::isfinite(worst)) {
            report.maxFeatureGap = std::max(report.maxFeatureGap, worst / modelExtent);
        }
        if (worst <= tol) ++report.featureChainsCovered;
    }
}

// ---------------------------------------------------------------------------
// check()
//
// The rest of Sec. 4's list, and the verdict. None of it is repaired here: the
// remedy the paper gives for every one of these is upstream -- raise lambda_5,
// add the Gamma_topo constraint for the pair that nearly met, re-run Stage 6
// from the current phi -- and MERIDIAN::traceAndRepair has already run it to
// exhaustion by the time these curves exist. What is worth doing here is
// saying exactly which face failed and how.
// ---------------------------------------------------------------------------
void Arrangement::check() {
    for (const Node &n : nodes) {
        // A node no half-edge leaves is not part of the layout. The only ones
        // are the loose ends trimUnresolved() cut away; counting them would
        // report a dangling node that is no longer there.
        if (n.out.empty() && n.kind != NodeKind::Cone) continue;
        ++report.nodes;
        switch (n.kind) {
            case NodeKind::Cone: ++report.coneNodes; break;
            case NodeKind::Crossing: ++report.crossingNodes; break;
            case NodeKind::BoundaryHit: ++report.boundaryHitNodes; break;
            case NodeKind::BoundaryCorner: ++report.boundaryCornerNodes; break;
            case NodeKind::InterfaceNode: ++report.interfaceNodes; break;
            case NodeKind::InterfaceHit: ++report.interfaceHitNodes; break;
            default: ++report.danglingNodes; break;
        }
        if (n.kind == NodeKind::Cone) {
            if (n.out.empty()) ++report.isolatedCones;
            else if (static_cast<int>(n.out.size()) != n.valence) ++report.coneValenceErrors;
        }
    }

    report.arcs = static_cast<int>(arcs.size());
    for (const Arc &a : arcs) {
        if (a.kind == ArcKind::Separatrix) ++report.separatrixArcs;
        else if (a.kind == ArcKind::Interface) ++report.interfaceArcs;
        else ++report.boundaryArcs;
        if (a.degenerate) ++report.degenerateArcs;
        if (a.dangling) ++report.danglingArcs;
    }
    for (size_t a = 0; a < arcs.size(); ++a) {
        const int f0 = halves[2 * a].face, f1 = halves[2 * a + 1].face;
        const bool p0 = f0 >= 0 && faces[f0].patch;
        const bool p1 = f1 >= 0 && faces[f1].patch;
        if (arcs[a].kind == ArcKind::Boundary) {
            if (p0 == p1) ++report.unsharedArcs;    // dS must have S on one side only
        } else if (!(p0 && p1)) {
            ++report.unsharedArcs;
        }
    }

    report.faces = static_cast<int>(faces.size());
    report.patches = static_cast<int>(patches.size());
    for (int f : patches) {
        if (faces[f].quad) ++report.quads;
        if (faces[f].simple) ++report.simpleQuads;
    }
    report.areaCoverage = report.domainArea > 0.0 ? report.patchArea / report.domainArea : 0.0;

    std::ostringstream oss;
    if (report.danglingNodes > 0) {
        oss << report.danglingNodes << " separatrix/ces stop in the middle of a face: Q5 did "
            << "not close them, so the layout there is not a partition. Sec. 4's remedy is "
            << "Sec. 3.3's -- a Gamma_topo constraint for the pair that nearly met, and "
            << "Stage 6 re-run from the current phi.";
        report.messages.push_back(oss.str());
        oss.str("");
    }
    if (report.wrongCornerFaces > 0) {
        oss << report.wrongCornerFaces << " of " << report.patches
            << " patch(es) do not have exactly four corners";
        if (report.tJunctions > 0) {
            oss << "; " << report.tJunctions << " node(s) are T-junctions, which is "
                << "quantisation that has not converged (raise lambda_5)";
        }
        oss << ".";
        report.messages.push_back(oss.str());
        oss.str("");
    }
    if (report.quads > report.simpleQuads) {
        oss << (report.quads - report.simpleQuads) << " patch(es) have four corners but more "
            << "than one arc on a side, so a side and its opposite cannot share a knot vector "
            << "and Stage 9 will not fit them as Coons patches.";
        report.messages.push_back(oss.str());
        oss.str("");
    }
    if (report.coneValenceErrors > 0 || report.isolatedCones > 0) {
        oss << report.coneValenceErrors << " cone(s) have a number of incident arcs their "
            << "index does not prescribe";
        if (report.isolatedCones > 0) oss << " and " << report.isolatedCones << " have none at all";
        oss << ".";
        report.messages.push_back(oss.str());
        oss.str("");
    }
    if (report.unsharedArcs > 0) {
        oss << report.unsharedArcs << " arc(s) do not separate two patches (or a patch from "
            << "dS), so the faces they bound are not a partition of S.";
        report.messages.push_back(oss.str());
        oss.str("");
    }
    if (report.mixedBoundaryFaces > 0) {
        oss << report.mixedBoundaryFaces << " face(s) have dS on both sides of them, so they "
            << "are partly inside the model and partly outside. Two arcs were ordered the "
            << "wrong way round at a node they leave along nearly the same direction.";
        report.messages.push_back(oss.str());
        oss.str("");
    }
    if (report.annularFaces > 0) {
        oss << report.annularFaces << " boundary loop(s) carry no node at all: the face inside "
            << "them is an annulus rather than a quadrilateral. Footnote 3's remedy -- "
            << "designate a regular point a singularity -- applies to the layout as well.";
        report.messages.push_back(oss.str());
        oss.str("");
    }
    if (options.checkFeatures && report.featureChainsCovered < report.featureChains) {
        oss << (report.featureChains - report.featureChainsCovered) << " of "
            << report.featureChains << " feature chain(s) are not a union of arcs; the worst "
            << "runs " << report.maxFeatureGap << " of the model from the nearest one.";
        report.messages.push_back(oss.str());
        oss.str("");
    }
    if (report.areaCoverage < 0.999 || report.areaCoverage > 1.001) {
        oss << "The patches cover " << report.areaCoverage << " of S rather than all of it.";
        report.messages.push_back(oss.str());
        oss.str("");
    }

    report.valid = report.patches > 0 && report.danglingNodes == 0 &&
                   report.wrongCornerFaces == 0 && report.simpleQuads == report.patches &&
                   report.coneValenceErrors == 0 && report.isolatedCones == 0 &&
                   report.unsharedArcs == 0 && report.annularFaces == 0 &&
                   report.mixedBoundaryFaces == 0 &&
                   report.areaCoverage > 0.999 && report.areaCoverage < 1.001 &&
                   report.featureChainsCovered == report.featureChains;
}

// ---------------------------------------------------------------------------
std::vector<Arrangement::Side> Arrangement::patchSides(int face) const {
    std::vector<Side> out;
    if (face < 0 || face >= static_cast<int>(faces.size())) return out;
    const Face &fc = faces[face];
    if (!fc.simple) return out;
    for (const std::vector<int> &sd : fc.sides) {
        const HalfEdge &h = halves[sd.front()];
        out.push_back(Side{h.arc, h.forward});
    }
    return out;
}

// ---------------------------------------------------------------------------
bool Arrangement::writeOBJ(const std::string &filename) const {
    std::ofstream out(filename);
    if (!out) return false;
    int base = 1;
    for (const Arc &a : arcs) {
        if (a.points.size() < 2) continue;
        for (const Point &p : a.points) out << "v " << p[0] << " " << p[1] << " 0\n";
        out << "l";
        for (size_t i = 0; i < a.points.size(); ++i) out << " " << (base + static_cast<int>(i));
        out << "\n";
        base += static_cast<int>(a.points.size());
    }
    return true;
}

bool Arrangement::writePatchOBJ(const std::string &filename) const {
    std::ofstream out(filename);
    if (!out) return false;
    int base = 1;
    for (int f : patches) {
        const Face &fc = faces[f];
        std::vector<Point> loop;
        for (int hid : fc.half) {
            const HalfEdge &h = halves[hid];
            const Arc &ar = arcs[h.arc];
            if (h.forward) {
                for (size_t k = 0; k + 1 < ar.points.size(); ++k) loop.push_back(ar.points[k]);
            } else {
                for (size_t k = ar.points.size(); k-- > 1;) loop.push_back(ar.points[k]);
            }
        }
        if (loop.size() < 3) continue;
        for (const Point &p : loop) out << "v " << p[0] << " " << p[1] << " 0\n";
        out << "l";
        for (size_t i = 0; i < loop.size(); ++i) out << " " << (base + static_cast<int>(i));
        out << " " << base << "\n";
        base += static_cast<int>(loop.size());
    }
    return true;
}
