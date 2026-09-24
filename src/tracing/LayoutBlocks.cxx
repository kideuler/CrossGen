#include "tracing/LayoutBlocks.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <unordered_map>
#include <unordered_set>

#include "geom/Coons.hxx"
#include "geom/Fitting.hxx"
#include "tracing/TriangleLocator.hxx"

namespace {

double polylineLength(const std::vector<Point> &p) {
    double L = 0.0;
    for (size_t i = 1; i < p.size(); ++i) L += normP(p[i] - p[i - 1]);
    return L;
}

} // namespace

LayoutBlocks::LayoutBlocks(const QuadLayout &layout, const Mesh &mesh, const Options &opts)
    : opts_(opts) {
    double total = 0.0;
    for (const auto &e : mesh.edges) total += normP(mesh.vertices[e[1]] - mesh.vertices[e[0]]);
    if (!mesh.edges.empty()) meanEdge_ = total / static_cast<double>(mesh.edges.size());

    const auto &nodes = layout.getNodes();
    const auto &larcs = layout.getArcs();
    nodePos_.resize(nodes.size());
    nodeOnBoundary_.assign(nodes.size(), 0);
    for (size_t n = 0; n < nodes.size(); ++n) {
        nodePos_[n] = nodes[n].pos;
        for (const int d : nodes[n].darts)
            if (larcs[QuadLayout::arcOfDart(d)].onBoundary) nodeOnBoundary_[n] = 1;
    }

    buildMacroArcs(layout);
    classifyFaces(layout);
    buildBlocks(layout);
    fitArcs();
    buildPatches();
    buildShape();
    buildDecomposition(mesh);
}

// ---------------------------------------------------------------------------
// buildMacroArcs()  --  chains through pass-through nodes
//
// A pass-through node is one where exactly two arcs meet and the model has no
// corner: two separatrices joined head-on, or what a chord collapse left of a
// node it merged away. QuadLayout already never calls one a corner; here the
// arcs either side of it become one side, so that a block is not refused for
// a node that is not a junction of anything.
// ---------------------------------------------------------------------------
void LayoutBlocks::buildMacroArcs(const QuadLayout &layout) {
    const auto &nodes = layout.getNodes();
    const auto &larcs = layout.getArcs();

    auto passThrough = [&](int n) {
        const auto &nd = nodes[n];
        return nd.darts.size() == 2 && nd.kind != QuadLayout::NodeKind::BoundaryCorner &&
               nd.kind != QuadLayout::NodeKind::InterfaceNode &&
               QuadLayout::arcOfDart(nd.darts[0]) != QuadLayout::arcOfDart(nd.darts[1]);
    };
    auto originOf = [&](int d) {
        const auto &a = larcs[QuadLayout::arcOfDart(d)];
        return (d & 1) ? a.b : a.a;
    };
    auto endOf = [&](int d) { return originOf(d ^ 1); };

    macroOfArc_.assign(larcs.size(), -1);
    auto walk = [&](int d0, bool closed) {
        MacroArc m;
        m.from = originOf(d0);
        int d = d0;
        for (size_t guard = 0; guard <= larcs.size(); ++guard) {
            const int a = QuadLayout::arcOfDart(d);
            if (macroOfArc_[a] >= 0) break;
            macroOfArc_[a] = static_cast<int>(arcs_.size());
            m.darts.push_back(d);
            const int n = endOf(d);
            if (!passThrough(n) || (closed && n == m.from)) break;
            const auto &nd = nodes[n];
            const int twin = d ^ 1;
            d = (nd.darts[0] == twin) ? nd.darts[1] : nd.darts[0];
        }
        m.to = endOf(m.darts.back());
        for (const int dd : m.darts) {
            const auto &a = larcs[QuadLayout::arcOfDart(dd)];
            if (a.onBoundary) m.onBoundary = true;
            if (a.onInterface) m.onInterface = true;
            std::vector<Point> pts = a.pts;
            if (dd & 1) std::reverse(pts.begin(), pts.end());
            if (!m.points.empty()) pts.erase(pts.begin());
            m.points.insert(m.points.end(), pts.begin(), pts.end());
        }
        m.length = polylineLength(m.points);
        arcs_.push_back(std::move(m));
    };

    // From every node that is not a pass-through, along every arc leaving it.
    for (size_t n = 0; n < nodes.size(); ++n) {
        if (passThrough(static_cast<int>(n))) continue;
        for (const int d : nodes[n].darts)
            if (macroOfArc_[QuadLayout::arcOfDart(d)] < 0) walk(d, false);
    }
    // What is left is a cycle of pass-throughs with no other node on it -- a
    // hole nothing landed on. It is a macro arc from a node to itself, and no
    // block can have it for a side.
    for (size_t a = 0; a < larcs.size(); ++a)
        if (macroOfArc_[a] < 0 && larcs[a].a >= 0) walk(static_cast<int>(2 * a), true);
}

// ---------------------------------------------------------------------------
// classifyFaces()  --  which components are conforming quadrilaterals
//
// The T-junctions are found structurally, by what the face walk already
// decided: a slot at which a component does *not* turn a corner, at a node
// where three or more arcs meet, is that component running straight past a
// separatrix that ended on its side (docs/viertel_2019.md Sec. 5's corner
// rule for a TJ). The node kinds the tracing assigned are not asked, because
// a chord collapse merges nodes and the kinds do not survive it.
// ---------------------------------------------------------------------------
void LayoutBlocks::classifyFaces(const QuadLayout &layout) {
    const auto &nodes = layout.getNodes();
    const auto &faces = layout.getFaces();
    report_.faces = static_cast<int>(faces.size());
    faceIsBlock_.assign(faces.size(), 0);

    std::unordered_set<int> tSet;
    for (size_t f = 0; f < faces.size(); ++f) {
        const QuadLayout::Face &fc = faces[f];
        report_.totalArea += fc.area;

        bool tHere = false;
        for (size_t i = 0; i < fc.nodes.size(); ++i) {
            if (fc.isCorner[i]) continue;
            const int n = fc.nodes[i];
            if (nodes[n].darts.size() < 3) continue;
            tHere = true;
            if (tSet.insert(n).second) {
                tNodes_.push_back(n);
                tPoints_.push_back(nodes[n].pos);
            }
        }

        std::unordered_set<int> seenArc, seenNode;
        bool disk = true;
        for (size_t i = 0; i < fc.darts.size(); ++i) {
            if (!seenArc.insert(QuadLayout::arcOfDart(fc.darts[i])).second) disk = false;
            if (!seenNode.insert(fc.nodes[i]).second) disk = false;
        }

        if (!disk) { ++report_.notDisks; report_.notDiskArea += fc.area; continue; }
        if (fc.corners != 4 || fc.sides.size() != 4) {
            ++report_.notFourSided;
            report_.notFourSidedArea += fc.area;
            continue;
        }
        if (tHere) { ++report_.tJunctionFaces; report_.tJunctionArea += fc.area; continue; }
        faceIsBlock_[f] = 1;
    }
    report_.tJunctions = static_cast<int>(tNodes_.size());
}

// ---------------------------------------------------------------------------
// buildBlocks()  --  a block per conforming quadrilateral
//
// Each side of a qualifying face is a run of darts between two corners with
// only pass-throughs between them, which is exactly one macro arc read one way
// or the other. That is checked rather than assumed; a face whose side is not
// is refused and said so, since it would mean the two walks disagree.
// ---------------------------------------------------------------------------
void LayoutBlocks::buildBlocks(const QuadLayout &layout) {
    const auto &faces = layout.getFaces();
    const auto &larcs = layout.getArcs();
    auto originOf = [&](int d) {
        const auto &a = larcs[QuadLayout::arcOfDart(d)];
        return (d & 1) ? a.b : a.a;
    };

    for (size_t f = 0; f < faces.size(); ++f) {
        if (!faceIsBlock_[f]) continue;
        const QuadLayout::Face &fc = faces[f];
        Block b;
        b.face = static_cast<int>(f);
        b.area = fc.area;
        bool ok = true;
        for (int k = 0; k < 4 && ok; ++k) {
            const std::vector<int> &side = fc.sides[k];
            if (side.empty()) { ok = false; break; }
            const int m = macroOfArc_[QuadLayout::arcOfDart(side.front())];
            if (m < 0 || arcs_[m].darts.size() != side.size() || arcs_[m].from == arcs_[m].to) {
                ok = false;
                break;
            }
            const MacroArc &ma = arcs_[m];
            bool fwd = true, bwd = true;
            for (size_t i = 0; i < side.size(); ++i) {
                if (side[i] != ma.darts[i]) fwd = false;
                if (side[i] != (ma.darts[ma.darts.size() - 1 - i] ^ 1)) bwd = false;
            }
            if (!fwd && !bwd) { ok = false; break; }
            b.sides[k] = m;
            b.forward[k] = fwd;
            b.corners[k] = originOf(side.front());
        }
        if (!ok) {
            faceIsBlock_[f] = 0;
            ++report_.notFourSided;
            report_.notFourSidedArea += fc.area;
            report_.messages.push_back("face " + std::to_string(f) +
                                       ": a side is not one macro arc; not a block");
            continue;
        }
        for (int k = 0; k < 4; ++k) arcs_[b.sides[k]].used = true;
        report_.coveredArea += b.area;
        blocks_.push_back(std::move(b));
    }
    report_.blocks = static_cast<int>(blocks_.size());
}

// ---------------------------------------------------------------------------
// fitArcs()  --  every macro arc once
//
// Stage 9's fit (geom::fitCurve: chord-length parameters, clamped uniform
// knots, both ends pinned to the nodes) at the fewest segments that stay
// within Options::fitTolerance mean edges of the traced polyline. The tolerance
// is in edges because the polyline is itself only a discretisation of the
// streamline at that resolution: a separatrix crosses each triangle as one
// chord, so fitting it much closer than a fraction of an edge fits the
// triangulation, not the field. Boundary arcs are carried exactly.
//
// Every macro arc with two distinct ends is fitted, not only the ones a block
// uses: the rest are separatrices the layout never closed a block on, and they
// go into the B-rep as loose edges so that the model written out shows the
// whole partition and not only the part of it that meshes.
// ---------------------------------------------------------------------------
void LayoutBlocks::fitArcs() {
    std::unordered_map<int, geom::Vertex> vertexOf;
    auto vertex = [&](int n) {
        auto it = vertexOf.find(n);
        if (it != vertexOf.end()) return it->second;
        geom::Vertex v(nodePos_[n]);
        vertexOf.emplace(n, v);
        return v;
    };

    const double tol = opts_.fitTolerance * meanEdge_;
    double devSum = 0.0;
    std::string firstRefusal;
    for (size_t i = 0; i < arcs_.size(); ++i) {
        MacroArc &m = arcs_[i];
        if (m.from == m.to || m.points.size() < 2 || !(m.length > 0.0)) continue;
        if (m.used) ++report_.macroArcs;

        const bool exact = (m.onBoundary && !opts_.fitBoundaryArcs) ||
                           (m.onInterface && !opts_.fitInterfaceArcs);
        if (!exact) {
            try {
                geom::FitOptions fo;
                fo.degree = 3;
                fo.pinEnds = true;
                fo.regularisation = opts_.regularisation;
                int segs = std::max(1, opts_.segments);
                for (;;) {
                    fo.segments = segs;
                    geom::FitResult<2> fit = geom::fitCurve(m.points, fo);
                    m.spline = std::move(fit.curve);
                    m.maxDeviation = fit.maxDeviation;
                    m.segments = segs;
                    if (m.maxDeviation <= tol || segs >= opts_.maxSegments) break;
                    segs = std::min(opts_.maxSegments, 2 * segs);
                }
                if (m.used) {
                    ++report_.fittedArcs;
                    if (m.maxDeviation > tol) ++report_.unconverged;
                    report_.maxDeviation = std::max(report_.maxDeviation, m.maxDeviation / meanEdge_);
                    report_.maxSegmentsUsed = std::max(report_.maxSegmentsUsed, m.segments);
                    devSum += m.maxDeviation / meanEdge_;
                }
            } catch (const std::exception &e) {
                // The polyline stands in for a curve the fit could not make.
                m.spline = geom::BSplineCurve<2>();
                report_.messages.push_back("macro arc " + std::to_string(i) +
                                           ": fit failed (" + e.what() + "), carried exactly");
            }
        }
        if (exact || m.spline.empty()) {
            m.exact = true;
            m.poly = geom::Polyline<2>(m.points);
            if (m.used) ++report_.exactArcs;
        }

        try {
            m.edge = geom::Edge(m.exact ? m.poly.curve() : m.spline, vertex(m.from), vertex(m.to));
        } catch (const std::exception &e) {
            ++report_.brepFailures;
            if (firstRefusal.empty()) firstRefusal = "macro arc " + std::to_string(i) + ": " + e.what();
        }
    }
    if (report_.fittedArcs > 0) report_.meanDeviation = devSum / report_.fittedArcs;
    if (!firstRefusal.empty())
        report_.messages.push_back("the kernel refused " + std::to_string(report_.brepFailures) +
                                   " macro arc(s) as edges; the first: " + firstRefusal);
}

// ---------------------------------------------------------------------------
// buildPatches()  --  the Coons patch of each block, and its face
//
// Oriented into the (s, t) frame the way Stage 9 orients an arrangement face:
//
//     corner 0 -> (0,0)   side 0, the bottom, s increasing
//     corner 1 -> (1,0)   side 1, the right,  t increasing
//     corner 2 -> (1,1)   side 2 reversed, the top
//     corner 3 -> (0,1)   side 3 reversed, the left
//
// The layout's faces are counter-clockwise, so the frame is right-handed and
// a patch that folds shows up as a sampled cell of negative area.
// ---------------------------------------------------------------------------
void LayoutBlocks::buildPatches() {
    int refused = 0;
    std::string firstRefusal;
    for (size_t bi = 0; bi < blocks_.size(); ++bi) {
        Block &b = blocks_[bi];
        std::array<bool, 4> alongFrame{};
        std::array<geom::BSplineCurve<2>, 4> side;
        bool edgesOk = true;
        for (int k = 0; k < 4; ++k) {
            const MacroArc &m = arcs_[b.sides[k]];
            alongFrame[k] = (k < 2) == b.forward[k];
            const geom::BSplineCurve<2> &c = m.exact ? m.poly.curve() : m.spline;
            side[k] = alongFrame[k] ? c : c.reversed();
            if (m.edge.isNull()) edgesOk = false;
        }
        try {
            b.surface = geom::coonsSurface(side[0], side[2], side[3], side[1]);
        } catch (const std::exception &e) {
            if (refused++ == 0) firstRefusal = "block " + std::to_string(bi) + ": " + e.what();
            continue;
        }
        if (edgesOk) {
            try {
                b.brepFace = geom::Face(b.surface,
                                        {arcs_[b.sides[0]].edge, arcs_[b.sides[1]].edge,
                                         arcs_[b.sides[2]].edge, arcs_[b.sides[3]].edge},
                                        alongFrame);
            } catch (const std::exception &e) {
                if (refused++ == 0) firstRefusal = "block " + std::to_string(bi) + ": " + e.what();
            }
        }

        // The patch sits on its own corners and its own sides: both follow
        // from the Coons blend, and both are read back.
        const std::array<Point, 4> uv{{Point{0.0, 0.0}, Point{1.0, 0.0}, Point{1.0, 1.0},
                                        Point{0.0, 1.0}}};
        for (int k = 0; k < 4; ++k) {
            const Point c = b.surface.evaluate(uv[k][0], uv[k][1]);
            report_.maxCornerGap = std::max(report_.maxCornerGap, normP(c - nodePos_[b.corners[k]]));
        }
        for (int i = 0; i <= 8; ++i) {
            const double u = i / 8.0;
            const Point onSide[4] = {side[0].evaluate(u), side[1].evaluate(u), side[2].evaluate(u),
                                     side[3].evaluate(u)};
            const Point onPatch[4] = {b.surface.evaluate(u, 0.0), b.surface.evaluate(1.0, u),
                                      b.surface.evaluate(u, 1.0), b.surface.evaluate(0.0, u)};
            for (int k = 0; k < 4; ++k)
                report_.maxSideGap = std::max(report_.maxSideGap, normP(onSide[k] - onPatch[k]));
        }

        // The fold check, on a sampled grid of the patch.
        const int n = std::max(2, opts_.patchSamples);
        std::vector<Point> grid(static_cast<size_t>(n + 1) * (n + 1));
        for (int j = 0; j <= n; ++j)
            for (int i = 0; i <= n; ++i)
                grid[static_cast<size_t>(j) * (n + 1) + i] =
                    b.surface.evaluate(static_cast<double>(i) / n, static_cast<double>(j) / n);
        double area = 0.0, worst = std::numeric_limits<double>::infinity();
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                const Point &q00 = grid[static_cast<size_t>(j) * (n + 1) + i];
                const Point &q10 = grid[static_cast<size_t>(j) * (n + 1) + i + 1];
                const Point &q11 = grid[static_cast<size_t>(j + 1) * (n + 1) + i + 1];
                const Point &q01 = grid[static_cast<size_t>(j + 1) * (n + 1) + i];
                const double cell = 0.5 * (cross2(q10 - q00, q11 - q00) + cross2(q11 - q00, q01 - q00));
                area += cell;
                worst = std::min(worst, cell);
            }
        }
        const double mean = area / (n * n);
        b.minCellRatio = (std::fabs(mean) > 0.0) ? worst / mean : 0.0;
        if (b.minCellRatio <= 0.0) ++report_.foldedPatches;
    }
    report_.brepFailures += refused;
    if (refused > 0)
        report_.messages.push_back("the kernel refused " + std::to_string(refused) +
                                   " block(s) as patches or faces; the first: " + firstRefusal);
}

// ---------------------------------------------------------------------------
// buildShape()  --  the B-rep
// ---------------------------------------------------------------------------
void LayoutBlocks::buildShape() {
    std::vector<geom::Face> faces;
    for (const Block &b : blocks_) if (!b.brepFace.isNull()) faces.push_back(b.brepFace);
    std::vector<geom::Edge> loose;
    for (const MacroArc &m : arcs_) if (!m.used && !m.edge.isNull()) loose.push_back(m.edge);
    if (faces.empty() && loose.empty()) return;
    try {
        shape_ = geom::Shape(faces, loose);
        report_.brepFaces = shape_.faceCount();
        report_.brepEdges = shape_.edgeCount();
        report_.brepSharedEdges = shape_.sharedEdgeCount();
        report_.brepFreeEdges = shape_.freeEdgeCount();
        report_.brepValid = shape_.isValid();
    } catch (const std::exception &e) {
        report_.messages.push_back(std::string("the B-rep could not be assembled: ") + e.what());
    }
}

// ---------------------------------------------------------------------------
// buildDecomposition()  --  the shared representation, on the splines
//
// A macro edge is stored in the direction its first block walks it, which
// makes that block its A side; the second block to claim it walks it the other
// way and reads it flipped. That is BlockDecomposition's convention and the one
// every other producer follows.
// ---------------------------------------------------------------------------
void LayoutBlocks::buildDecomposition(const Mesh &mesh) {
    decomp_ = BlockDecomposition();
    decomp_.source = "Viertel";

    std::unordered_map<int, int> vidx;
    auto macroVertex = [&](int n) {
        auto it = vidx.find(n);
        if (it != vidx.end()) return it->second;
        BlockDecomposition::MacroVertex mv;
        mv.p = nodePos_[n];
        mv.onBoundary = nodeOnBoundary_[n] != 0;
        const int id = static_cast<int>(decomp_.vertices.size());
        decomp_.vertices.push_back(mv);
        vidx.emplace(n, id);
        return id;
    };
    auto sampled = [&](const MacroArc &m) {
        if (m.exact) return m.poly.points();
        const int n = std::max(8, m.segments * std::max(2, opts_.samplesPerSegment));
        std::vector<Point> p = m.spline.sample(n);
        if (!p.empty()) {
            p.front() = nodePos_[m.from];
            p.back() = nodePos_[m.to];
        }
        return p;
    };

    decomp_.blocks.resize(blocks_.size());
    for (size_t bi = 0; bi < blocks_.size(); ++bi) {
        const Block &b = blocks_[bi];
        BlockDecomposition::Block &ob = decomp_.blocks[bi];
        for (int k = 0; k < 4; ++k) {
            ob.corners[k] = macroVertex(b.corners[k]);
            MacroArc &m = arcs_[b.sides[k]];
            if (m.decompositionEdge < 0) {
                BlockDecomposition::MacroEdge oe;
                oe.boundary = m.onBoundary;
                oe.points = sampled(m);
                if (!b.forward[k]) std::reverse(oe.points.begin(), oe.points.end());
                oe.from = macroVertex(b.corners[k]);
                oe.to = macroVertex(b.corners[(k + 1) % 4]);
                oe.blockA = static_cast<int>(bi);
                oe.sideA = k;
                m.decompositionEdge = static_cast<int>(decomp_.edges.size());
                decomp_.edges.push_back(std::move(oe));
            } else {
                BlockDecomposition::MacroEdge &oe = decomp_.edges[m.decompositionEdge];
                oe.blockB = static_cast<int>(bi);
                oe.sideB = k;
            }
            ob.edges[k] = m.decompositionEdge;
            ob.flip[k] = decomp_.edges[m.decompositionEdge].blockB == static_cast<int>(bi);
        }
    }

    // --- materials, read off the model ------------------------------------
    bool multiMaterial = false;
    for (size_t t = 1; t < mesh.triangleMatId.size(); ++t)
        if (mesh.triangleMatId[t] != mesh.triangleMatId[0]) { multiMaterial = true; break; }
    if (!mesh.triangles.empty() && !mesh.triangleMatId.empty()) {
        const TriangleLocator grid(mesh);
        std::unordered_set<int> seen;
        for (size_t b = 0; b < decomp_.blocks.size(); ++b) {
            std::map<int, int> votes;
            Point q;
            if (decomp_.interiorPoint(static_cast<int>(b), q)) {
                const int t = grid.locate(q);
                if (t >= 0) votes[mesh.triangleMatId[t]] += 2;
            }
            if (multiMaterial) {
                for (int j = 1; j <= 3; ++j) {
                    for (int i = 1; i <= 3; ++i) {
                        const int t = grid.locate(decomp_.coonsPoint(static_cast<int>(b), i / 4.0, j / 4.0));
                        if (t >= 0) votes[mesh.triangleMatId[t]] += 1;
                    }
                }
            }
            int best = 0, bestVotes = -1;
            for (const auto &[mat, v] : votes) if (v > bestVotes) { best = mat; bestVotes = v; }
            decomp_.blocks[b].material = best;
            blocks_[b].material = best;
            if (votes.size() > 1) ++report_.straddlingBlocks;
            if (bestVotes >= 0) seen.insert(best);
        }
        report_.materials = static_cast<int>(seen.size());
    }

    for (BlockDecomposition::MacroEdge &oe : decomp_.edges) {
        oe.matLeft = oe.blockA >= 0 ? decomp_.blocks[oe.blockA].material : 0;
        oe.matRight = oe.blockB >= 0 ? decomp_.blocks[oe.blockB].material : 0;
        oe.interface = !oe.boundary && oe.blockA >= 0 && oe.blockB >= 0 && oe.matLeft != oe.matRight;
        for (const int v : {oe.from, oe.to}) {
            if (v < 0) continue;
            ++decomp_.vertices[v].valence;
            if (oe.boundary) decomp_.vertices[v].onBoundary = true;
            if (oe.interface) decomp_.vertices[v].onInterface = true;
        }
    }
}
