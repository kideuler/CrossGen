// Checks on SKELETON (src/SKELETON/SKELETON.hxx), the medial-axis quad mesher
// of docs/medial_axis.md (Cheny et al., "Planar quad meshing with
// guarantees"), and on the blocks, mesh and smoothed mesh it is finished with.
//
//   TestSKELETON <mesh.obj> [--h size] [--mesh h] [--tmop sweeps]
//                [--blocks-obj f] [--obj f]
//   TestSKELETON --selftest
//
// On a model it runs the SKELETON class -- the one the viewer's mode 3 runs --
// and checks, in order: the dual complex of Step 3 (integrity, boundary
// orientation, the inside filter); Step 4's simplify() (displacement budget,
// locked corners left alone, cocircular rings); the branches of Step 5.1; the
// invariants of Sec. 15 on Steps 5-10 -- one contact per medial edge at every
// vertex, every face four-cornered with positive area, the faces tiling the
// domain, a consistent quantization, and the valence identity
// sum_int (4 - deg) + sum_bdry (3 - deg) = 4 chi on the quad mesh it
// prescribes; and the blocks, that mesh's base complex: every side a macro
// edge between its own block's corners, every side on dS laid on the model's
// own polyline, and the mesh's singular vertices the blocks' irregular ones.
// How much of the model the blocks cover is printed, with the
// docs/block_decomposition_metrics.md scores, not asserted: it is what the
// method achieves on the model.
//
//   --h size          the target size h of Steps 4 and 6-9 (default a tenth
//                     of the bounding-box diagonal)
//   --mesh h          mesh the blocks with BlockQuadMesh at target edge length
//                     h, and check that the mesh conforms
//   --tmop sweeps     then smooth it with mesh::TMOP (implies --mesh 0.05)
//   --blocks-obj f    the blocks' macro edges as OBJ polylines
//   --obj f           the mesh (smoothed with --tmop) as an OBJ
//   --no-block-smoothing  read the blocks off the trace mesh as Step 10 left
//                     it (SKELETON::Options::blockSmoothing)
//   --block-subdivisions k  trace the Step 10 mesh at k times its resolution
//                     (SKELETON::Options::blockSubdivisions, default 4)
//   --boundary-stages n  stages that lay the trace mesh's boundary onto dS
//                     where a block side crosses it (Options::boundaryStages,
//                     default 4; 0 never)
//
// TestSKELETON --selftest builds the unit-test domains of Sec. 15 (square,
// 4x2 rectangle, equilateral triangle, annulus, ellipse, rounded rectangle,
// L-shape, noisy square), sampled at s = h / 8 (Step 2's default), checks what
// the guide expects of each, and that each one's blocks cover it and mesh
// conformingly at h / 2.

#include <iostream>
#include <iomanip>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <cmath>
#include <cstdlib>

#include "mesh/BoundaryFeatures.hxx"
#include "mesh/FeatureFrame.hxx"
#include "mesh/Mesh.hxx"
#include "mesh/QuadMesh.hxx"
#include "triangle/TriangleMesher.hpp"
#include "SKELETON/SKELETON.hxx"

// ANSI colour helpers
#define GREEN "\033[32m"
#define RED   "\033[31m"
#define YELLOW "\033[33m"
#define RESET "\033[0m"

// ─── Check harness ──────────────────────────────────────────────────────────
namespace {

int g_fails = 0;

void check(bool ok, const std::string& label, const std::string& detail = "") {
    std::cout << (ok ? GREEN "[PASS]" RESET : RED "[FAIL]" RESET) << " " << label;
    if (!detail.empty()) std::cout << " (" << detail << ")";
    std::cout << "\n";
    if (!ok) ++g_fails;
}

// Every structural invariant the dual complex is supposed to maintain. Run
// before and after simplification, since a merge is the operation most likely
// to break one.
void checkComplexIntegrity(const MedialAxis& ma, const std::string& when) {
    int badRing = 0, badKey = 0, asymmetric = 0, badCount = 0, repeatedVert = 0;
    std::unordered_map<int, int> observedVertexCells;

    for (int ci = 0; ci < static_cast<int>(ma.cells.size()); ++ci) {
        const auto& cell = ma.cells[ci];
        if (!cell.active) continue;

        const int n = static_cast<int>(cell.verts.size());
        if (n < 3) { ++badRing; continue; }

        std::unordered_set<int> seen;
        for (int v : cell.verts) {
            if (!seen.insert(v).second) ++repeatedVert;
            ++observedVertexCells[v];
        }

        for (int i = 0; i < n; ++i) {
            const auto e = ma.ringEdge(ci, i);
            if (!ma.edgeCells.count(MeshEdgeKey(e[0], e[1]))) { ++badKey; continue; }

            const int nb = ma.neighborAcross(ci, i);
            if (nb < 0) continue;
            // The neighbour must name this cell back across the same edge.
            bool mutual = false;
            const int nn = static_cast<int>(ma.cells[nb].verts.size());
            for (int j = 0; j < nn; ++j) {
                if (ma.neighborAcross(nb, j) == ci) { mutual = true; break; }
            }
            if (!mutual) ++asymmetric;
        }
    }

    for (const auto& kv : observedVertexCells) {
        if (ma.vertexCellCount(kv.first) != kv.second) ++badCount;
    }

    check(badRing == 0, when + ": every active cell has >= 3 ring vertices",
          std::to_string(badRing) + " bad");
    check(repeatedVert == 0, when + ": cell rings are simple (no repeated vertex)",
          std::to_string(repeatedVert) + " repeats");
    check(badKey == 0, when + ": every ring edge is present in edgeCells",
          std::to_string(badKey) + " missing");
    check(asymmetric == 0, when + ": cell adjacency is symmetric",
          std::to_string(asymmetric) + " one-sided");
    check(badCount == 0, when + ": vertexCellCount matches the rings",
          std::to_string(badCount) + " stale");
}

void checkAxisConsistency(const MedialAxis& ma, const std::string& when) {
    const int nv = static_cast<int>(ma.medialVertices.size());
    int oob = 0, nonFinite = 0, badBackref = 0, degreeMismatch = 0, selfLoop = 0;

    for (const auto& e : ma.medialEdges) {
        if (e[0] < 0 || e[0] >= nv || e[1] < 0 || e[1] >= nv) ++oob;
        else if (e[0] == e[1]) ++selfLoop;
    }

    for (int i = 0; i < nv; ++i) {
        const auto& mv = ma.medialVertices[i];
        if (!std::isfinite(mv.coord[0]) || !std::isfinite(mv.coord[1])) ++nonFinite;
        if (mv.degree != static_cast<int>(mv.incidentEdges.size())) ++degreeMismatch;

        // Each incident edge must actually name this vertex, and the dual cell
        // must point back at it.
        for (int e : mv.incidentEdges) {
            if (e < 0 || e >= static_cast<int>(ma.medialEdges.size())) { ++badBackref; continue; }
            if (ma.medialEdges[e][0] != i && ma.medialEdges[e][1] != i) ++badBackref;
        }
        if (mv.cell < 0 || mv.cell >= static_cast<int>(ma.cells.size()) ||
            ma.cells[mv.cell].medialVertex != i) {
            ++badBackref;
        }
    }

    check(oob == 0, when + ": medial edge indices in range", std::to_string(oob) + " bad");
    check(selfLoop == 0, when + ": no self-loop medial edges", std::to_string(selfLoop) + " found");
    check(nonFinite == 0, when + ": medial vertex coordinates finite",
          std::to_string(nonFinite) + " bad");
    check(degreeMismatch == 0, when + ": degree == incidentEdges.size()",
          std::to_string(degreeMismatch) + " bad");
    check(badBackref == 0, when + ": incidentEdges / cell back-references agree",
          std::to_string(badBackref) + " bad");
}


// Step 4: what simplify() promises, whatever it chose to do.
void checkSimplification(const MedialAxis& ma, int cellsBefore) {
    const SimplifyStats& ss = ma.simplifyStats;
    std::cout << "  Cocircular merges " << ss.cocircularMerges << ", endpoints pruned "
              << ss.endpointMerges << " (" << ss.endpointsKept << " kept), arcs flattened "
              << ss.flattenedChains << " (" << ss.removedVertices << " vertices, "
              << ss.deletedTwoGons << " 2-gons), rounds " << ss.rounds << "\n";
    std::cout << "  Largest move " << ss.maxDisplacement << " of delta_max " << ss.budget << "\n";

    int active = 0;
    double drift = 0.0;
    for (const auto& cell : ma.cells) {
        if (!cell.active) continue;
        ++active;
        for (int v : cell.verts) {
            drift = std::max(drift, std::abs(normP(ma.mesh->vertices[v] - cell.center) - cell.radius));
        }
    }
    check(cellsBefore - active ==
              ss.cocircularMerges + ss.endpointMerges + ss.deletedTwoGons,
          "deactivated cells == merges + deleted 2-gons",
          std::to_string(cellsBefore - active));

    // Every boundary vertex within delta_max of where the input put it, and
    // every locked feature exactly there, still on the loop.
    const double slack = 1e-12 * ma.boundingBoxDiagonal();
    int overBudget = 0, lockedMoved = 0;
    for (int v = 0; v < static_cast<int>(ma.boundaryLinks.size()); ++v) {
        const bool onLoop = ma.boundaryLinks[v].nextBoundaryVertex >= 0;
        const double moved = normP(ma.mesh->vertices[v] - ma.originalPosition[v]);
        if (onLoop && moved > ss.budget + slack) ++overBudget;
        if (ma.locked[v] && (!onLoop || moved > slack)) ++lockedMoved;
    }
    check(overBudget == 0, "every boundary vertex within delta_max of x0",
          std::to_string(overBudget) + " over");
    check(lockedMoved == 0, "locked feature vertices neither moved nor removed",
          std::to_string(lockedMoved) + " disturbed");
    // Merged rings sit on the survivor's stored circle (Step 3.2): exactly,
    // up to rounding, where 4.1 projected them, and to the cocircular
    // tolerance (a hundredth of delta_max) where Step 2.3's merge found them
    // there.
    check(drift <= 0.01 * ss.budget + slack, "merged cells are cocircular on their stored circle",
          "worst drift " + std::to_string(drift) + " of delta_max " + std::to_string(ss.budget));
}

// Valences of the Step 10 quad mesh: edges per vertex, i.e. quads around an
// interior vertex and one more than that on the boundary.
struct Valences {
    std::map<int, int> interior, boundary;
    int interiorIrregular = 0, boundaryIrregular = 0;
    long long identity = 0;   // sum_int (4 - deg) + sum_bdry (3 - deg)
    std::vector<std::pair<Point, int>> irregularInterior;
};

Valences meshValences(const QuadExtraction& qe) {
    std::vector<int> quads(qe.vertices.size(), 0);
    for (const auto& q : qe.quads) {
        for (int v : q) ++quads[v];
    }
    Valences out;
    for (size_t v = 0; v < qe.vertices.size(); ++v) {
        if (qe.onBoundary[v]) {
            const int deg = quads[v] + 1;
            ++out.boundary[deg];
            out.identity += 3 - deg;
            if (deg != 3) ++out.boundaryIrregular;
        } else {
            const int deg = quads[v];
            ++out.interior[deg];
            out.identity += 4 - deg;
            if (deg != 4) {
                ++out.interiorIrregular;
                out.irregularInterior.push_back({qe.vertices[v], deg});
            }
        }
    }
    return out;
}

// Quads with a corner whose Jacobian is not positive.
int invertedQuads(const QuadExtraction& qe) {
    int inverted = 0;
    for (const auto& q : qe.quads) {
        for (int k = 0; k < 4; ++k) {
            const Point& c = qe.vertices[q[k]];
            if (cross2(qe.vertices[q[(k + 1) % 4]] - c, qe.vertices[q[(k + 3) % 4]] - c) <= 0.0) {
                ++inverted;
                break;
            }
        }
    }
    return inverted;
}

std::string histogram(const std::map<int, int>& h) {
    std::string s;
    for (const auto& kv : h) {
        s += (s.empty() ? "" : " ") + std::to_string(kv.second) + "x" + std::to_string(kv.first);
    }
    return s.empty() ? "none" : s;
}

double signedArea(const std::vector<Point>& pts) {
    double a = 0.0;
    for (size_t i = 0; i < pts.size(); ++i) a += cross2(pts[i], pts[(i + 1) % pts.size()]);
    return 0.5 * a;
}

// The invariants of Sec. 15 on Steps 5-10.
void checkSteps5to10(const SKELETON& sk, const std::string& when) {
    const MedialAxis& ma = sk.getAxis();
    const MedialAxisTMesh& tm = sk.getTMesh();
    const TMeshStats& st = tm.stats();
    std::cout << "  h " << tm.targetSize() << ": " << st.keptVertices << " vertices kept, "
              << st.insertedVertices << " inserted, " << st.zones << " zones (" << st.capZones
              << " caps), " << st.coarseEdges << " coarse edges: " << st.colorCounts[0]
              << " green, " << st.colorCounts[1] << " red, " << st.colorCounts[2] << " blue, "
              << st.colorCounts[3] << " purple; " << tm.blocks.size() << " faces, "
              << st.midpointCuts << " midpoint cuts, " << st.dissolvedVertices
              << " valence-2 dissolved, " << st.refinements << " refinements\n";

    // Step 5.2: a vertex has one contact per incident medial edge, so the
    // walk passes it once per edge.
    int badContacts = 0;
    for (int m = 0; m < tm.nodeCount(); ++m) {
        if (!tm.isKept(m)) continue;
        if (m < static_cast<int>(ma.medialVertices.size()) && ma.medialVertices[m].degree == 0) continue;
        if (tm.contactCount(m) != tm.nodeDegree(m)) ++badContacts;
    }
    check(badContacts == 0, when + ": #contacts == valence at every coarse vertex",
          std::to_string(badContacts) + " off");
    check(st.skippedLoops == 0 && st.unpairedEdges == 0,
          when + ": every loop walked, every coarse edge has two sides",
          std::to_string(st.skippedLoops) + " loops, " + std::to_string(st.unpairedEdges) + " edges");

    // Step 8: four corners, positive area -- and Step 8.6 found nothing it
    // could not repair.
    int notQuad = 0, folded = 0, leftovers = 0;
    double blockArea = 0.0;
    for (const auto& b : tm.blocks) {
        if (b.edge < 0) ++leftovers;
        else if (b.corners.size() != 4) ++notQuad;
        const double a = signedArea(b.outline);
        if (a <= 0.0) ++folded;
        blockArea += a;
    }
    // A cell no template claimed -- the caps of a disc, whose axis is a
    // point; a zone around a tip the sampling was too coarse for -- keeps the
    // corners it has, and the welder gives a three-cornered one a fourth.
    if (leftovers > 0) {
        std::cout << "  " << leftovers << " cells no template claimed\n";
    }
    check(notQuad == 0, when + ": every template face has four corners",
          std::to_string(notQuad) + " not");
    check(folded == 0, when + ": every face is CCW with positive area", std::to_string(folded) + " folded");
    check(st.badFaces == 0, when + ": Step 8.6 geometry check passes",
          std::to_string(st.badFaces) + " bad faces left");

    // The faces tile the domain: their areas add up to the area the input
    // boundary encloses.
    double domainArea = 0.0;
    {
        std::vector<char> seen(ma.originalNext.size(), 0);
        for (int v0 = 0; v0 < static_cast<int>(ma.originalNext.size()); ++v0) {
            if (seen[v0] || ma.originalNext[v0] < 0) continue;
            std::vector<Point> loop;
            for (int v = v0; v >= 0 && !seen[v]; v = ma.originalNext[v]) {
                seen[v] = 1;
                loop.push_back(ma.originalPosition[v]);
            }
            domainArea += signedArea(loop);
        }
    }
    check(std::abs(blockArea - domainArea) <= 1e-9 * std::abs(domainArea),
          when + ": faces tile the domain (areas add up)",
          std::to_string(blockArea) + " vs " + std::to_string(domainArea));

    // Step 9.
    const BlockQuant& bq = sk.getQuantization();
    const TMeshQuantizer::Report& quant = sk.getQuantReport();
    check(bq.ok && bq.skippedBlocks == 0, when + ": every face welded into the T-mesh",
          bq.error + " " + std::to_string(bq.skippedBlocks) + " dropped");
    if (!bq.ok) return;
    check(quant.consistent, when + ": quantization consistent (Ax = 0)");
    check(quant.unmetRequirements == 0,
          when + ": Prop. 1 requirements met (radius chains, features, loops >= 1)",
          std::to_string(quant.unmetRequirements) + " of " +
              std::to_string(bq.requirements.size()) + " unmet");
    int zero = 0, total = 0;
    for (const auto& e : bq.tmesh.edges) {
        zero += e.x == 0;
        total += e.x;
    }
    std::cout << "  Quantized: " << bq.requirements.size() << " requirements ("
              << bq.radiusRequirements << " radius chains, " << bq.featureRequirements
              << " features, " << bq.loopRequirements << " loops), "
              << quant.stage1Vectors << " stage-I strips, " << quant.stage2Moves
              << " stage-II moves, objective " << quant.objective << ", " << zero << " of "
              << bq.tmesh.edges.size() << " T-mesh edges collapsed\n";

    // Step 10: a manifold quad mesh with the domain's loops and Euler
    // characteristic, and the valence identity that goes with it.
    int loops = 0;
    {
        std::vector<char> seen(ma.boundaryLinks.size(), 0);
        for (int v0 = 0; v0 < static_cast<int>(ma.boundaryLinks.size()); ++v0) {
            if (seen[v0] || ma.boundaryLinks[v0].nextBoundaryVertex < 0) continue;
            ++loops;
            for (int u = v0; u >= 0 && !seen[u]; u = ma.boundaryLinks[u].nextBoundaryVertex) seen[u] = 1;
        }
    }
    const int chi = 2 - loops;
    if (!sk.hasQuadMesh()) {
        check(false, when + ": Step 10 produced a quad mesh");
        return;
    }
    const QuadExtraction& qe = sk.getQuadMesh();
    std::cout << "  Quad mesh: " << qe.vertices.size() << " vertices, " << qe.quads.size()
              << " quads, " << qe.collapsedFaces << " faces collapsed, " << qe.droppedQuads
              << " degenerate quads dropped, " << qe.doublets << " doublets merged, "
              << qe.reorientedQuads << " reoriented, " << invertedQuads(qe) << " inverted\n";
    check(qe.ok, when + ": quad mesh is manifold", qe.error);
    check(qe.boundaryLoops == loops, when + ": quad mesh keeps the boundary loops",
          std::to_string(qe.boundaryLoops) + " vs " + std::to_string(loops));
    check(qe.eulerCharacteristic == chi, when + ": Euler characteristic of the domain",
          std::to_string(qe.eulerCharacteristic) + " vs " + std::to_string(chi));
    const Valences v = meshValences(qe);
    std::cout << "  Valences: interior " << histogram(v.interior) << "; boundary "
              << histogram(v.boundary) << "\n";
    check(v.identity == 4 * chi, when + ": valence identity sum(4 - deg) + sum(3 - deg) = 4 chi",
          std::to_string(v.identity) + " vs " + std::to_string(4 * chi));
}

// The blocks: the Step 10 mesh's base complex, as the shared decomposition.
// Sound means every block side resolves to a macro edge that runs between
// that block's own corners and names the block back, every block outline is
// counter-clockwise, every side on dS is the model's own polyline, and --
// where every patch became a block -- the macrovertices that are not
// four-valent inside are exactly the Step 10 mesh's interior singularities.
void checkBlocks(const SKELETON& sk, const std::string& when) {
    if (!sk.hasBlocks()) {
        check(false, when + ": blocks built from the Step 10 mesh");
        return;
    }
    const BaseComplex::Report& r = sk.getBlocks().getReport();
    const BlockDecomposition& D = sk.getDecomposition();
    std::cout << "  Blocks: " << r.blocks << " of " << r.patches << " patches, "
              << D.edges.size() << " macro edges, " << D.vertices.size()
              << " macrovertices; on a " << r.quads << "-quad trace mesh: " << r.interiorSingular
              << " interior and " << r.boundarySingular << " boundary singular vertices, "
              << r.corners << " corners (" << r.missedCorners << " with no vertex), "
              << r.separatrices << " separatrices, " << r.cuts << " cuts in " << r.cutRounds
              << " rounds, " << r.snappedVertices << " vertices snapped onto dS (up to "
              << r.maxSnap << " of the diagonal); covering " << std::fixed << std::setprecision(1) << 100.0 * r.coverage
              << "%" << std::defaultfloat << std::setprecision(6) << "\n";

    const mesh::TMOP::Report& sr = sk.getBlockSmoothingReport();
    if (!sr.ran && sr.sweeps + sr.untangleSweeps > 0)
        std::cout << "  Trace mesh NOT smoothed: " << sr.sweeps << " sweeps after " << sr.untangleSweeps
                  << " untangling left " << sr.invertedAfter << " folds\n";
    if (sr.ran)
        std::cout << std::fixed << std::setprecision(3) << "  Trace mesh "
                  << SKELETON::name(sk.blockGeometry()) << ": " << sr.sweeps
                  << " sweeps, folds " << sr.invertedBefore << " -> " << sr.invertedAfter
                  << ", worst scaled Jacobian " << sr.minScaledJacobianBefore << " -> "
                  << sr.minScaledJacobianAfter << std::defaultfloat << std::setprecision(6) << "\n";
    // How jagged the sides are: the mean total turning of an interior side
    // between its two corners.
    // The turning is all of it, a smooth bend included. A kink is a vertex
    // whose turn differs from the mean of its two neighbours' by more than 15
    // degrees: a side bending steadily round a hole turns by as much at every
    // vertex and has none, a zigzag has one at every vertex it zigzags at.
    {
        double turning = 0.0;
        int interior = 0, kinks = 0;
        for (const auto& E : D.edges) {
            if (E.boundary) continue;
            ++interior;
            std::vector<double> turn(E.points.size(), 0.0);
            for (size_t i = 1; i + 1 < E.points.size(); ++i) {
                const Point a = E.points[i] - E.points[i - 1], c = E.points[i + 1] - E.points[i];
                if (!(normP(a) > 0.0 && normP(c) > 0.0)) continue;
                turn[i] = std::atan2(cross2(a, c), dotP(a, c));
                turning += std::abs(turn[i]);
            }
            for (size_t i = 1; i + 1 < E.points.size(); ++i) {
                if (std::abs(turn[i] - 0.5 * (turn[i - 1] + turn[i + 1])) > M_PI / 12.0) ++kinks;
            }
        }
        std::cout << "  Interior sides: " << interior << ", turning " << std::fixed << std::setprecision(3)
                  << (interior ? turning / interior : 0.0) << " rad on average, " << kinks
                  << " kinks" << std::defaultfloat << std::setprecision(6) << "\n";
    }

    int badSides = 0, folded = 0;
    for (size_t b = 0; b < D.blocks.size(); ++b) {
        const BlockDecomposition::Block& B = D.blocks[b];
        std::vector<Point> outline;
        for (int s = 0; s < 4; ++s) {
            const BlockDecomposition::MacroEdge& E = D.edges[B.edges[s]];
            const int from = B.flip[s] ? E.to : E.from, to = B.flip[s] ? E.from : E.to;
            const bool named = B.flip[s] ? (E.blockB == static_cast<int>(b) && E.sideB == s)
                                         : (E.blockA == static_cast<int>(b) && E.sideA == s);
            if (from != B.corners[s] || to != B.corners[(s + 1) % 4] || !named) ++badSides;
            const std::vector<Point> side = D.sidePolyline(static_cast<int>(b), s);
            outline.insert(outline.end(), side.begin(), side.empty() ? side.end() : side.end() - 1);
        }
        if (signedArea(outline) <= 0.0) ++folded;
    }
    check(badSides == 0, when + ": every block side is a macro edge between its own corners",
          std::to_string(badSides) + " off");
    // A block outline that is not counter-clockwise runs through a fold of
    // the Step 10 mesh (the inverted quads counted above): the
    // decomposition is sound and its geometry is what Step 10 left, for the
    // mesher and TMOP to untangle.
    if (folded > 0)
        std::cout << YELLOW "[WARN]" RESET " " << folded
                  << " block outline(s) not counter-clockwise, through a fold of the Step 10 mesh\n";
    // A side across dS is a defect of the decomposition, and one this stage
    // cannot always mend: what is left after the boundary is laid onto dS is
    // a hole of a few trace vertices or a concave corner the mesh rounds off
    // (Step 11), so it is reported rather than asserted.
    if (r.sidesCrossingBoundary > 0)
        std::cout << YELLOW "[WARN]" RESET " " << r.sidesCrossingBoundary
                  << " interior side(s) cross dS\n";
    check(r.chordSides == 0, when + ": every side on dS laid on the model's own boundary",
          std::to_string(r.chordSides) + " of " + std::to_string(r.boundarySides) +
              " left as chords");
    if (r.badPatches == 0) {
        int irregular = 0;
        for (const auto& mv : D.vertices) irregular += !mv.onBoundary && mv.valence != 4;
        check(irregular == r.interiorSingular,
              when + ": the blocks' irregular interior vertices are Step 10's singularities",
              std::to_string(irregular) + " vs " + std::to_string(r.interiorSingular));
    }
}

// BlockQuadMesh on the blocks: conforming, with what it folded printed.
void checkMesh(const SKELETON& sk, const std::string& when) {
    const BlockQuadMesh::Report& qr = sk.getBlockMesh().getReport();
    std::cout << std::fixed << std::setprecision(4) << "  Mesh: " << qr.quads << " quads on "
              << qr.vertices << " vertices, " << qr.chords << " chords at " << qr.minIntervals
              << "-" << qr.maxIntervals << " edges; scaled Jacobian " << qr.minScaledJacobian
              << " worst, " << qr.meanScaledJacobian << " mean, " << qr.invertedQuads
              << " folded; boundary deviation " << qr.boundaryDeviation
              << std::defaultfloat << std::setprecision(6) << "\n";
    for (const std::string& m : qr.messages) std::cout << "    " << m << "\n";
    check(qr.conforming, when + ": the mesh on the blocks conforms");
}

// ─── Self-test domains (Sec. 15) ────────────────────────────────────────────

// A closed polygon sampled at spacing <= s, every corner kept.
std::vector<Point> samplePolygon(const std::vector<Point>& corners, double s) {
    std::vector<Point> out;
    for (size_t i = 0; i < corners.size(); ++i) {
        const Point& a = corners[i];
        const Point& b = corners[(i + 1) % corners.size()];
        const int n = std::max(1, static_cast<int>(std::ceil(normP(b - a) / s - 1e-9)));
        for (int k = 0; k < n; ++k) out.push_back(a + (b - a) * (static_cast<double>(k) / n));
    }
    return out;
}

// An ellipse sampled at spacing about s.
std::vector<Point> sampleEllipse(const Point& c, double a, double b, double s) {
    const double perimeter = M_PI * (3.0 * (a + b) - std::sqrt((3.0 * a + b) * (a + 3.0 * b)));
    const int n = std::max(8, static_cast<int>(std::ceil(perimeter / s)));
    std::vector<Point> out;
    for (int k = 0; k < n; ++k) {
        const double t = 2.0 * M_PI * k / n;
        out.push_back({c[0] + a * std::cos(t), c[1] + b * std::sin(t)});
    }
    return out;
}

// A rectangle with rounded corners of radius rho.
std::vector<Point> sampleRoundedRect(double w, double hgt, double rho, double s) {
    std::vector<Point> out;
    const Point centres[4] = {{w - rho, rho}, {w - rho, hgt - rho}, {rho, hgt - rho}, {rho, rho}};
    const int arcN = std::max(4, static_cast<int>(std::ceil(0.5 * M_PI * rho / s)));
    for (int c = 0; c < 4; ++c) {
        const double t0 = -0.5 * M_PI + 0.5 * M_PI * c;
        for (int k = 0; k < arcN; ++k) {
            const double t = t0 + 0.5 * M_PI * k / arcN;
            out.push_back({centres[c][0] + rho * std::cos(t), centres[c][1] + rho * std::sin(t)});
        }
        const Point a{centres[c][0] + rho * std::cos(t0 + 0.5 * M_PI),
                      centres[c][1] + rho * std::sin(t0 + 0.5 * M_PI)};
        const Point& nc = centres[(c + 1) % 4];
        const double t1 = t0 + 0.5 * M_PI;
        const Point b{nc[0] + rho * std::cos(t1), nc[1] + rho * std::sin(t1)};
        const int n = std::max(1, static_cast<int>(std::ceil(normP(b - a) / s)));
        for (int k = 0; k < n; ++k) out.push_back(a + (b - a) * (static_cast<double>(k) / n));
    }
    return out;
}

// The constrained Delaunay triangulation of sampled loops, the first one the
// outer boundary and the rest holes, with no interior points (Step 3.1).
std::shared_ptr<Mesh> triangulateLoops(const std::vector<std::vector<Point>>& loops) {
    std::vector<std::array<double, 2>> vertlist;
    std::vector<std::vector<std::array<int, 2>>> segmentLoops;
    std::vector<int> types;
    for (size_t li = 0; li < loops.size(); ++li) {
        const int base = static_cast<int>(vertlist.size());
        const int n = static_cast<int>(loops[li].size());
        std::vector<std::array<int, 2>> segs;
        for (int i = 0; i < n; ++i) {
            vertlist.push_back(loops[li][i]);
            segs.push_back({base + i, base + (i + 1) % n});
        }
        segmentLoops.push_back(std::move(segs));
        types.push_back(li == 0 ? 0 : 1);
    }
    triangle_wrapper::TriangleMesher2D::Options opts;
    opts.just_delaunay = true;
    triangle_wrapper::TriangleMesher2D mesher(opts);
    triangle_wrapper::TriangleMesher2D::MeshInput input;
    input.vertlist = vertlist;
    input.segment_loops = segmentLoops;
    input.type = types;
    input.h = 0.0;
    auto output = mesher.triangulate(input);
    return std::make_shared<Mesh>(output.verts, output.triangles);
}

struct Shape {
    int endpoints = 0;
    int junctions = 0;      // vertices of degree >= 3
    int maxJunctionDegree = 0;
};

Shape axisShape(const MedialAxis& ma) {
    Shape s;
    for (const auto& mv : ma.medialVertices) {
        if (mv.degree == 1) ++s.endpoints;
        if (mv.degree >= 3) {
            ++s.junctions;
            s.maxJunctionDegree = std::max(s.maxJunctionDegree, mv.degree);
        }
    }
    return s;
}

std::unique_ptr<SKELETON> runDomain(const std::string& name,
                                    const std::vector<std::vector<Point>>& loops, double h) {
    std::cout << "\n── " << name << " (h = " << h << ", s = h/8) ──\n";
    SKELETON::Options o;
    o.targetSize = h;
    auto sk = std::make_unique<SKELETON>(triangulateLoops(loops), o);
    sk->triangulate();
    const int cellsBefore = static_cast<int>(sk->getTriangulation().triangles.size());
    sk->buildAxis();
    checkSimplification(sk->getAxis(), cellsBefore);
    sk->buildTMesh();
    if (sk->quantize()) sk->extractQuads();
    checkSteps5to10(*sk, name);
    sk->buildBlocks();
    checkBlocks(*sk, name);
    if (sk->hasBlocks()) {
        const BaseComplex::Report& r = sk->getBlocks().getReport();
        check(r.badPatches == 0 && sk->getDecomposition().covers(),
              name + ": the blocks cover the domain",
              std::to_string(r.badPatches) + " patches not blocks");
        BlockQuadMesh::Options mo;
        mo.targetEdgeLength = 0.5 * h;
        if (sk->buildMesh(mo)) checkMesh(*sk, name);
        else check(false, name + ": the blocks mesh");
    }
    return sk;
}

Valences valencesOf(const SKELETON& sk) {
    return sk.hasQuadMesh() ? meshValences(sk.getQuadMesh()) : Valences();
}

int blockCount(const SKELETON& sk) {
    return sk.hasBlocks() ? static_cast<int>(sk.getDecomposition().blocks.size()) : 0;
}

int runSelfTest() {
    // Square: an X of four red diagonals ending in four blue corners, the
    // centre the outer end of all four, so nothing is singular but the
    // corners themselves.
    {
        const double h = 0.25;
        auto p = runDomain("square", {samplePolygon({{0, 0}, {1, 0}, {1, 1}, {0, 1}}, h / 8)}, h);
        const Shape a = axisShape(p->getAxis());
        const TMeshStats& st = p->getTMesh().stats();
        check(a.endpoints == 4 && a.junctions == 1 && a.maxJunctionDegree == 4,
              "square: X-shaped axis (4 endpoints, one junction of degree 4)",
              std::to_string(a.endpoints) + " endpoints, " + std::to_string(a.junctions) +
                  " junctions of degree up to " + std::to_string(a.maxJunctionDegree));
        check(st.colorCounts[2] == 4 && st.colorCounts[3] == 0 && st.colorCounts[0] == 0,
              "square: blue corners, red diagonals");
        check(valencesOf(*p).interiorIrregular == 0, "square: no interior singular vertex",
              std::to_string(valencesOf(*p).interiorIrregular));
        check(valencesOf(*p).boundaryIrregular == 4 && valencesOf(*p).boundary[2] == 4,
              "square: the four corners are the only boundary irregulars (valence 2)",
              histogram(valencesOf(*p).boundary));
    }
    // 4 x 2 rectangle: two junctions joined by a green edge, the rest as the
    // square -- a regular grid.
    {
        const double h = 0.5;
        auto p = runDomain("4x2 rectangle", {samplePolygon({{0, 0}, {4, 0}, {4, 2}, {0, 2}}, h / 8)}, h);
        const Shape a = axisShape(p->getAxis());
        const TMeshStats& st = p->getTMesh().stats();
        check(a.endpoints == 4 && a.junctions == 2 && a.maxJunctionDegree == 3,
              "rectangle: two junctions of degree 3, four endpoints");
        check(st.colorCounts[0] >= 1 && st.colorCounts[2] == 4,
              "rectangle: green between the junctions, blue corners",
              std::to_string(st.colorCounts[0]) + " green");
        check(valencesOf(*p).interiorIrregular == 0, "rectangle: regular grid inside",
              std::to_string(valencesOf(*p).interiorIrregular));
        check(valencesOf(*p).boundaryIrregular == 4, "rectangle: only the corners irregular on the boundary",
              histogram(valencesOf(*p).boundary));
    }
    // Equilateral triangle: three quads around a valence-3 centre.
    {
        const double h = 0.2;
        const double r3 = std::sqrt(3.0) / 2.0;
        auto p = runDomain("equilateral triangle",
                               {samplePolygon({{0, 0}, {1, 0}, {0.5, r3}}, h / 8)}, h);
        const Shape a = axisShape(p->getAxis());
        check(a.endpoints == 3 && a.junctions == 1 && a.maxJunctionDegree == 3,
              "triangle: three branches from one junction");
        check(p->getTMesh().stats().colorCounts[2] == 3, "triangle: blue corners");
        const bool centre = valencesOf(*p).interiorIrregular == 1 &&
                            valencesOf(*p).irregularInterior[0].second == 3 &&
                            normP(valencesOf(*p).irregularInterior[0].first - Point{0.5, r3 / 3.0}) < 1e-6;
        check(centre, "triangle: the one interior singularity is valence 3 at the centre",
              histogram(valencesOf(*p).interior));
    }
    // Annulus: an all-green cycle, a regular ring.
    {
        const double h = 0.25;
        std::vector<Point> inner = sampleEllipse({0, 0}, 0.5, 0.5, h / 8);
        std::reverse(inner.begin(), inner.end());
        auto p = runDomain("annulus", {sampleEllipse({0, 0}, 1.0, 1.0, h / 8), inner}, h);
        const TMeshStats& st = p->getTMesh().stats();
        check(axisShape(p->getAxis()).endpoints == 0, "annulus: no endpoints");
        check(st.colorCounts[1] == 0 && st.colorCounts[2] == 0 && st.colorCounts[3] == 0 &&
                  st.colorCounts[0] >= 3,
              "annulus: all green, at least three coarse edges");
        check(valencesOf(*p).interiorIrregular == 0 && valencesOf(*p).boundaryIrregular == 0,
              "annulus: regular ring");
    }
    // Ellipse: a short segment with a purple cap at each end; the caps are
    // valence-3 vertices, and nothing else is worse than 3 or 5.
    {
        const double h = 0.2;
        auto p = runDomain("ellipse", {sampleEllipse({0, 0}, 1.0, 0.6, h / 8)}, h);
        const TMeshStats& st = p->getTMesh().stats();
        check(axisShape(p->getAxis()).endpoints == 2, "ellipse: one branch, two endpoints");
        check(st.colorCounts[3] == 2 && st.colorCounts[2] == 0, "ellipse: purple at both ends");
        bool only35 = valencesOf(*p).interiorIrregular > 0;
        for (const auto& kv : valencesOf(*p).interior) {
            if (kv.first != 4 && kv.first != 3 && kv.first != 5) only35 = false;
        }
        check(only35 && valencesOf(*p).boundaryIrregular == 0,
              "ellipse: singularities of valence 3 or 5 only, none on the boundary",
              histogram(valencesOf(*p).interior));
    }
    // Rounded rectangle: a purple cap at each fillet centre.
    {
        const double h = 0.25;
        auto p = runDomain("rounded rectangle", {sampleRoundedRect(3.0, 1.5, 0.3, h / 8)}, h);
        const TMeshStats& st = p->getTMesh().stats();
        check(st.colorCounts[3] == 4 && st.colorCounts[2] == 0, "rounded rectangle: four purple caps",
              std::to_string(st.colorCounts[3]) + " purple, " + std::to_string(st.colorCounts[2]) + " blue");
    }
    // L-shape: a blue cap at each of the five convex corners. The concave
    // corner is left to Step 11, which this implementation does not have.
    {
        const double h = 0.25;
        auto p = runDomain("L-shape", {samplePolygon({{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}}, h / 8)}, h);
        check(p->getTMesh().stats().colorCounts[2] == 5, "L-shape: five blue corners",
              std::to_string(p->getTMesh().stats().colorCounts[2]));
    }
    // Noisy square: Step 4 prunes the branches the noise grows (Fig. 10),
    // leaving the clean square's four endpoints.
    {
        const double h = 0.25;
        std::vector<Point> pts = samplePolygon({{0, 0}, {1, 0}, {1, 1}, {0, 1}}, h / 8);
        std::mt19937 rng(12345);
        std::uniform_real_distribution<double> noise(-0.25, 0.25);
        const double amplitude = 0.05 * h;   // delta_max
        for (Point& q : pts) {
            const bool corner = (q[0] == 0.0 || q[0] == 1.0) && (q[1] == 0.0 || q[1] == 1.0);
            if (corner) continue;
            // Along the side's normal, a quarter of delta_max at most.
            if (q[1] == 0.0 || q[1] == 1.0) q[1] += noise(rng) * amplitude;
            else q[0] += noise(rng) * amplitude;
        }
        const int rawEndpoints = axisShape(MedialAxis(triangulateLoops({pts}))).endpoints;
        auto p = runDomain("noisy square", {pts}, h);
        std::cout << "  Raw axis endpoints: " << rawEndpoints << "\n";
        check(axisShape(p->getAxis()).endpoints == 4, "noisy square: Step 4 leaves the four corner endpoints",
              std::to_string(rawEndpoints) + " -> " + std::to_string(axisShape(p->getAxis()).endpoints));
    }
    return g_fails;
}

} // namespace

// ─── Main ───────────────────────────────────────────────────────────────────
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <mesh.obj> [--h size] [--mesh h] [--tmop sweeps] [--blocks-obj f] [--obj f]"
                     " [--no-block-smoothing] [--block-subdivisions k] [--boundary-stages n]"
                     " | --selftest\n";
        return 1;
    }

    if (std::string(argv[1]) == "--selftest") {
        runSelfTest();
        std::cout << "\n";
        if (g_fails == 0) std::cout << GREEN "[ALL PASSED]" RESET "\n";
        else              std::cout << RED "[" << g_fails << " FAILURE(S)]" RESET "\n";
        return g_fails > 0 ? 1 : 0;
    }

    const std::string path = argv[1];
    double h = 0.0, meshTarget = 0.0;
    int sweeps = 0;
    bool smoothBlocks = true;
    int subdivisions = 0, stages = -1;
    std::string blocksObj, obj;
    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&]() -> const char* {
            if (i + 1 >= argc) {
                std::cerr << a << " needs a value\n";
                std::exit(1);
            }
            return argv[++i];
        };
        if (a == "--h") h = std::atof(next());
        else if (a == "--mesh") meshTarget = std::atof(next());
        else if (a == "--tmop") sweeps = std::atoi(next());
        else if (a == "--blocks-obj") blocksObj = next();
        else if (a == "--obj") obj = next();
        else if (a == "--no-block-smoothing") smoothBlocks = false;
        else if (a == "--block-subdivisions") subdivisions = std::atoi(next());
        else if (a == "--boundary-stages") stages = std::atoi(next());
        else {
            std::cerr << "Unknown option " << a << "\n";
            return 1;
        }
    }
    if (sweeps > 0 && meshTarget <= 0.0) meshTarget = 0.05;
    std::cout << "Loading mesh: " << path << "\n";

    auto originalMesh = std::make_shared<Mesh>(path);
    std::cout << "  Original mesh: " << originalMesh->triangles.size() << " triangles, "
              << originalMesh->vertices.size() << " vertices, "
              << originalMesh->boundaryEdges.size() << " boundary edges\n";

    SKELETON::Options options;
    options.targetSize = h;
    if (sweeps > 0) options.tmop.maxSweeps = sweeps;
    if (!smoothBlocks) options.blockSmoothing.maxSweeps = 0;
    if (subdivisions > 0) options.blockSubdivisions = subdivisions;
    if (stages >= 0) options.boundaryStages = stages;
    SKELETON sk(originalMesh, options);
    try {
        sk.triangulate();
    } catch (const std::exception& e) {
        std::cerr << RED "[FAIL]" RESET " " << e.what() << "\n";
        return 1;
    }
    std::cout << "  Delaunay mesh:  " << sk.getTriangulation().triangles.size() << " triangles, "
              << sk.getTriangulation().vertices.size() << " vertices, "
              << sk.getTriangulation().boundaryEdges.size() << " boundary edges; "
              << sk.getStatus().boundaryLoops << " boundary loops, " << sk.getStatus().corners
              << " corners\n";

    // ── Raw complex (Step 3), on a copy of the triangulation, before Step 4
    // moves anything ──
    const int cellsBefore = static_cast<int>(sk.getTriangulation().triangles.size());
    {
    auto rawAxis = std::make_shared<MedialAxis>(std::make_shared<Mesh>(sk.getTriangulation()));
    const MedialAxis& ma = *rawAxis;
    const std::shared_ptr<Mesh>& delaunayMesh = ma.mesh;
    std::cout << "  Dual cells:      " << ma.cells.size() << "\n";
    std::cout << "  Medial vertices: " << ma.medialVertices.size() << "\n";
    std::cout << "  Medial edges:    " << ma.medialEdges.size() << "\n";
    std::cout << "  Stats: " << ma.stats.degenerateCells << " degenerate, "
              << ma.stats.centersOutsideDomain << " centers outside, "
              << ma.stats.edgesCrossingBoundary << " edges crossing boundary, "
              << ma.stats.nonManifoldBoundaryVertices << " pinched boundary vertices\n";
    std::cout << "\n── Raw complex ──\n";

    check(ma.cells.size() == delaunayMesh->triangles.size(),
          "one dual cell per triangle",
          std::to_string(ma.cells.size()) + " vs " + std::to_string(delaunayMesh->triangles.size()));

    check(ma.medialVertices.size() == ma.cells.size(),
          "one medial vertex per active cell",
          std::to_string(ma.medialVertices.size()) + " vs " + std::to_string(ma.cells.size()));

    // The inside filter can only remove candidates, never invent them.
    {
        int numInternal = 0;
        for (const auto& et : delaunayMesh->edgeTriangles) {
            if (et[0] != -1 && et[1] != -1) ++numInternal;
        }
        const int kept = static_cast<int>(ma.medialEdges.size());
        check(kept + ma.stats.edgesCrossingBoundary == numInternal,
              "kept edges + filtered edges == internal Delaunay edges",
              std::to_string(kept) + " + " + std::to_string(ma.stats.edgesCrossingBoundary) +
              " vs " + std::to_string(numInternal));
    }

    // Cells are normalised to CCW, which is what makes ring order, incident
    // edge order and boundary orientation agree.
    {
        int cw = 0;
        for (const auto& cell : ma.cells) {
            if (!cell.active) continue;
            double area2 = 0.0;
            const int n = static_cast<int>(cell.verts.size());
            for (int i = 0; i < n; ++i) {
                const Point& p = delaunayMesh->vertices[cell.verts[i]];
                const Point& q = delaunayMesh->vertices[cell.verts[(i + 1) % n]];
                area2 += cross2(p, q);
            }
            if (area2 <= 0.0) ++cw;
        }
        check(cw == 0, "every cell ring is CCW", std::to_string(cw) + " reversed");
    }

    checkComplexIntegrity(ma, "raw");
    checkAxisConsistency(ma, "raw");

    // ── Boundary orientation ──
    // next/prev must keep the interior on the left. That makes the outer loop
    // CCW and hole loops CW, and it makes every interior angle land in
    // (0, 2*pi) with sum (k-2)*pi over an outer loop.
    {
        std::vector<bool> seen(delaunayMesh->vertices.size(), false);
        int loops = 0, ccwLoops = 0, badAngleSum = 0, unlinked = 0;
        double totalArea = 0.0;

        for (int bv : delaunayMesh->boundaryVertices) {
            if (ma.boundaryLinks[bv].nextBoundaryVertex == -1 ||
                ma.boundaryLinks[bv].prevBoundaryVertex == -1) {
                ++unlinked;
            }
        }

        for (int bv : delaunayMesh->boundaryVertices) {
            if (seen[bv] || ma.boundaryLinks[bv].nextBoundaryVertex == -1) continue;
            ++loops;

            double area2 = 0.0, angleSum = 0.0;
            int count = 0, cur = bv;
            do {
                seen[cur] = true;
                const int nxt = ma.boundaryLinks[cur].nextBoundaryVertex;
                if (nxt == -1) break;
                area2 += cross2(delaunayMesh->vertices[cur], delaunayMesh->vertices[nxt]);
                angleSum += ma.interiorAngle[cur];
                ++count;
                cur = nxt;
            } while (cur != bv && count < static_cast<int>(delaunayMesh->vertices.size()) + 1);

            totalArea += 0.5 * area2;
            if (area2 > 0.0) ++ccwLoops;
            // The domain-side angles of a CCW outer loop sum to (k-2)*pi. Around
            // a hole the domain is on the outside, so its k reflex angles sum
            // to (k+2)*pi instead.
            const double expected = (area2 > 0.0 ? (count - 2.0) : (count + 2.0)) * M_PI;
            if (std::abs(angleSum - expected) > 1e-6 * std::max(1.0, std::abs(expected))) {
                ++badAngleSum;
            }
        }

        check(unlinked == 0, "every boundary vertex has next and prev",
              std::to_string(unlinked) + " unlinked");
        check(loops > 0, "boundary loops found", std::to_string(loops) + " loops");
        check(ccwLoops == 1, "exactly one CCW loop (the outer one; holes run CW)",
              std::to_string(ccwLoops) + " CCW of " + std::to_string(loops));
        check(totalArea > 0.0, "signed area of all loops is positive (interior on the left)",
              std::to_string(totalArea));
        check(badAngleSum == 0, "interior angles sum to (k-2)*pi per loop",
              std::to_string(badAngleSum) + " loops off");
    }

    // Medial vertices must lie inside the domain to be medial at all.
    {
        int outside = 0;
        for (const auto& cell : ma.cells) {
            if (cell.active && !cell.centerInside) ++outside;
        }
        check(outside == ma.stats.centersOutsideDomain,
              "centersOutsideDomain matches the cell flags",
              std::to_string(outside) + " vs " + std::to_string(ma.stats.centersOutsideDomain));
        std::cout << (outside == 0 ? GREEN "[INFO]" RESET : YELLOW "[WARN]" RESET)
                  << " " << outside << " circumcenters lie outside the domain"
                  << (outside ? "; the boundary sampling is too coarse for a faithful axis" : "")
                  << "\n";
    }
    }

    // ── Step 4 ──
    std::cout << "\n── After simplification (Step 4) ──\n";
    sk.buildAxis();
    {
        const MedialAxis& ma = sk.getAxis();
        std::cout << "  Medial vertices: " << sk.getStatus().rawMedialVertices << " -> "
                  << ma.medialVertices.size() << "\n";
        std::cout << "  Medial edges:    " << ma.medialEdges.size() << "\n";
        checkSimplification(ma, cellsBefore);
        checkComplexIntegrity(ma, "simplified");
        checkAxisConsistency(ma, "simplified");
    }

    // ── Steps 5-10 ──
    sk.buildTMesh();
    if (sk.quantize()) sk.extractQuads();
    const MedialAxis& ma = sk.getAxis();

    // ── Branches ──
    std::cout << "\n── Branches ──\n";
    check(ma.polyLineEdges.size() == ma.polyLines.size(),
          "polyLineEdges is parallel to polyLines");

    int cyclic = 0, pureCycles = 0;
    for (size_t i = 0; i < ma.polyLines.size(); ++i) {
        if (!ma.polyLineIsCycle[i]) continue;
        ++cyclic;
        // A cycle made only of regular vertices has no seed of degree != 2, so
        // it is reachable only by the dedicated cycle pass.
        if (ma.medialVertices[ma.polyLines[i].front()].degree == 2) ++pureCycles;
    }
    std::cout << "  Polylines: " << ma.polyLines.size() << " (" << cyclic << " cyclic, "
              << pureCycles << " of them pure degree-2 cycles)\n";

    check(ma.polyLines.size() == ma.polyLineIsCycle.size(),
          "polyLineIsCycle is parallel to polyLines");

    // Every medial edge must appear in exactly one branch. This is the check
    // that catches both a branch traced twice and a cycle never traced at all.
    {
        std::map<std::pair<int, int>, int> edgeUse;
        for (const auto& e : ma.medialEdges) {
            edgeUse[{std::min(e[0], e[1]), std::max(e[0], e[1])}] = 0;
        }
        int unknown = 0;
        for (const auto& pl : ma.polyLines) {
            for (size_t j = 1; j < pl.size(); ++j) {
                const auto key = std::make_pair(std::min(pl[j - 1], pl[j]),
                                                std::max(pl[j - 1], pl[j]));
                auto it = edgeUse.find(key);
                if (it == edgeUse.end()) ++unknown;
                else ++it->second;
            }
        }
        int never = 0, twice = 0;
        for (const auto& kv : edgeUse) {
            if (kv.second == 0) ++never;
            else if (kv.second > 1) ++twice;
        }
        check(unknown == 0, "branches only use real medial edges",
              std::to_string(unknown) + " phantom");
        check(never == 0, "every medial edge is covered by a branch",
              std::to_string(never) + " uncovered");
        check(twice == 0, "no medial edge is covered twice",
              std::to_string(twice) + " duplicated");
    }

    // Interior vertices of an open branch must be regular; a cycle is closed.
    {
        int badInterior = 0, badCycle = 0;
        for (size_t i = 0; i < ma.polyLines.size(); ++i) {
            const auto& pl = ma.polyLines[i];
            for (size_t j = 1; j + 1 < pl.size(); ++j) {
                if (ma.medialVertices[pl[j]].degree != 2) ++badInterior;
            }
            if (ma.polyLineIsCycle[i]) {
                if (pl.size() < 3 || pl.front() != pl.back()) ++badCycle;
            } else if (pl.size() > 1 && pl.front() == pl.back()) {
                ++badCycle;
            }
        }
        check(badInterior == 0, "branch interiors are all degree-2 vertices",
              std::to_string(badInterior) + " bad");
        check(badCycle == 0, "cycle flag matches closure", std::to_string(badCycle) + " bad");
    }
    // ── Steps 5-10 ──
    std::cout << "\n── Boundary map, T-mesh, quantization, quad mesh (Steps 5-10) ──\n";
    checkSteps5to10(sk, "model");

    // ── Blocks ──
    std::cout << "\n── Blocks (the Step 10 mesh's base complex) ──\n";
    sk.buildBlocks();
    checkBlocks(sk, "model");
    if (sk.hasBlocks()) {
        const BlockDecomposition& D = sk.getDecomposition();
        const BoundaryFeatures features(*originalMesh);
        const FeatureFrame frame(*originalMesh);
        std::cout << std::fixed << std::setprecision(3) << "  Scores: regularity "
                  << D.regularity(features) << ", angle " << D.angleQuality() << ", chord "
                  << D.chordQuality() << ", alignment " << D.alignmentQuality(frame)
                  << (D.covers() ? "" : " (NaN: the blocks do not cover)") << "\n"
                  << std::defaultfloat << std::setprecision(6);
        if (!blocksObj.empty()) D.writeEdgesOBJ(blocksObj);
    }

    // ── Mesh and smoothing ──
    if (meshTarget > 0.0 && sk.hasBlocks()) {
        std::cout << "\n── Mesh (BlockQuadMesh at " << meshTarget << ")"
                  << (sweeps > 0 ? " and TMOP" : "") << " ──\n";
        BlockQuadMesh::Options mo;
        mo.targetEdgeLength = meshTarget;
        if (sk.buildMesh(mo)) {
            checkMesh(sk, "model");
            if (sweeps > 0 && sk.smooth()) {
                const mesh::TMOP::Report& tm = sk.getTMOPReport();
                std::cout << std::fixed << std::setprecision(4) << "  TMOP: " << tm.sweeps
                          << " sweeps, scaled Jacobian " << tm.minScaledJacobianBefore << " -> "
                          << tm.minScaledJacobianAfter << " worst, " << tm.meanScaledJacobianBefore
                          << " -> " << tm.meanScaledJacobianAfter << " mean, folds "
                          << tm.invertedBefore << " -> " << tm.invertedAfter << "\n"
                          << std::defaultfloat << std::setprecision(6);
                std::cout << "        " << tm.untangleSweeps << " untangling sweeps, "
                          << tm.frozenCorners << " frozen corners, nodes " << tm.freeNodes
                          << " free / " << tm.slidingNodes << " sliding / " << tm.fixedNodes
                          << " fixed" << (tm.converged ? "" : ", stopped on the sweep cap") << "\n";
                for (const std::string &msg : tm.messages) std::cout << "        " << msg << "\n";
            }
            if (!obj.empty()) {
                if (sk.hasSmoothedMesh()) sk.getSmoothedMesh().writeOBJ(obj);
                else sk.getBlockMesh().writeOBJ(obj);
            }
        } else {
            check(false, "model: the blocks mesh");
        }
    }

    // ── Summary ──
    std::cout << "\n";
    if (g_fails == 0) std::cout << GREEN "[ALL PASSED]" RESET "\n";
    else              std::cout << RED "[" << g_fails << " FAILURE(S)]" RESET "\n";
    return g_fails > 0 ? 1 : 0;
}
