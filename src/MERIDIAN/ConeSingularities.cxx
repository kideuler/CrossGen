#include "ConeSingularities.hxx"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <sstream>
#include <stdexcept>

#include "dualmbo/DualMBO.hxx"

namespace {

// Rotation of the cross between two triangles, wrapped into (-pi/4, pi/4].
// The field is stored as exp(4 i theta), so the wrapped difference of the
// *representations* is four times the wrapped difference of the directions.
inline double crossTurn(const std::complex<double> &from, const std::complex<double> &to) {
    const std::complex<double> d = to * std::conj(from);
    return std::atan2(d.imag(), d.real()) * 0.25;
}

inline double cornerAngle(const Mesh &m, int t, int v) {
    const Triangle &tri = m.triangles[t];
    int i0 = -1;
    for (int i = 0; i < 3; ++i) if (tri[i] == v) i0 = i;
    if (i0 < 0) return 0.0;
    const Point a = m.vertices[tri[(i0 + 1) % 3]] - m.vertices[v];
    const Point b = m.vertices[tri[(i0 + 2) % 3]] - m.vertices[v];
    return std::fabs(std::atan2(cross2(a, b), dotP(a, b)));
}

} // namespace

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
ConeSingularities::ConeSingularities(DualMBO &dualMBO)
    : ConeSingularities(dualMBO.getMeshPtr(), dualMBO.u_k_prev) {
    // Populate the solver's own singularity list from the same field, so that
    // anything else reading DualMBO::singularVertices later agrees with the cones
    // prescribed here.
    dualMBO.computeSingularities();
}

ConeSingularities::ConeSingularities(std::shared_ptr<Mesh> m, const Eigen::VectorXcd &crossField)
    : mesh(std::move(m)) {
    if (!mesh) throw std::runtime_error("ConeSingularities: null mesh");
    if (mesh->triangles.empty()) throw std::runtime_error("ConeSingularities: empty mesh");
    if (crossField.size() != static_cast<Eigen::Index>(mesh->triangles.size())) {
        throw std::runtime_error("ConeSingularities: cross field has " +
                                 std::to_string(crossField.size()) + " entries for " +
                                 std::to_string(mesh->triangles.size()) +
                                 " triangles (was the solver stepped?)");
    }

    const int nV = static_cast<int>(mesh->vertices.size());
    active.assign(nV, 0);
    for (const Triangle &t : mesh->triangles) {
        active[t[0]] = 1;
        active[t[1]] = 1;
        active[t[2]] = 1;
    }
    index.assign(nV, 0);
    rawIndex.assign(nV, 0.0);
    measured.assign(nV, 0);
    movedByRebalance.assign(nV, 0);
    prescribed.assign(nV, 0);
    targetCurvature.assign(nV, 0.0);

    computeInputCurvature();
    computeInteriorIndices(crossField);
    computeBoundaryIndices(crossField);
    rebuildCones();
}

// ---------------------------------------------------------------------------
// computeInputCurvature()  --  Eq. (6)
//
// On a planar mesh the interior values are zero up to round-off; they are
// computed rather than assumed so that the Eq. (7) check in gaussBonnet() is
// a real check on the triangulation (a fold or a duplicated vertex shows up
// here) and not a restatement of the assumption.
// ---------------------------------------------------------------------------
void ConeSingularities::computeInputCurvature() {
    const int nV = static_cast<int>(mesh->vertices.size());
    std::vector<double> angleSum(nV, 0.0);

    for (int t = 0; t < static_cast<int>(mesh->triangles.size()); ++t) {
        const Triangle &tri = mesh->triangles[t];
        for (int i = 0; i < 3; ++i) angleSum[tri[i]] += cornerAngle(*mesh, t, tri[i]);
    }

    inputCurvature.assign(nV, 0.0);
    for (int v = 0; v < nV; ++v) {
        if (!active[v]) continue; // not part of the surface; carries no curvature
        const double full = mesh->isBoundaryVertex[v] ? M_PI : 2.0 * M_PI;
        inputCurvature[v] = full - angleSum[v];
    }
}

// ---------------------------------------------------------------------------
// computeInteriorIndices()
//
// The same winding number DualMBO::computeSingularities() reports, kept here in
// its integer form: sum of the wrapped cross rotations around a closed star is
// an exact multiple of pi/2, so I(v) = (2/pi) * that sum needs no rounding.
// ---------------------------------------------------------------------------
void ConeSingularities::computeInteriorIndices(const Eigen::VectorXcd &field) {
    const int nV = static_cast<int>(mesh->vertices.size());
    const auto &vt = mesh->vertexTriangles;

    for (int v = 0; v < nV; ++v) {
        if (mesh->isBoundaryVertex[v]) continue;

        const int begin = vt.rowPtr[v];
        const int end = vt.rowPtr[v + 1];
        const int star = end - begin;
        if (star < 3) continue;

        double theta = 0.0;
        for (int k = begin; k < end; ++k) {
            const int tCur = vt.colIdx[k];
            const int tNext = vt.colIdx[(k - begin + 1) % star + begin];
            theta += crossTurn(field[tCur], field[tNext]);
        }

        rawIndex[v] = theta / M_PI_2;
        measured[v] = 1;
        index[v] = static_cast<int>(std::lround(rawIndex[v]));
    }
}

// ---------------------------------------------------------------------------
// boundaryFan()
//
// The CSR star is sorted by absolute angle around the vertex, which orders a
// closed interior loop correctly and an open boundary fan only by accident --
// a fan straddling the -pi/pi cut comes out shuffled. So the fan is walked
// through the triangle-edge tables instead, entering on one boundary edge and
// stopping on the other, then reversed if that walk happened to run clockwise.
// ---------------------------------------------------------------------------
std::vector<int> ConeSingularities::boundaryFan(int v) const {
    std::vector<int> fan;

    int firstEdge = -1;
    int boundaryDegree = 0;
    for (int e : mesh->boundaryEdges) {
        if (mesh->edges[e][0] == v || mesh->edges[e][1] == v) {
            if (boundaryDegree == 0) firstEdge = e;
            ++boundaryDegree;
        }
    }
    if (boundaryDegree != 2) return fan; // a pinch: no single corner to report

    int entry = firstEdge;
    int cur = (mesh->edgeTriangles[entry][0] >= 0) ? mesh->edgeTriangles[entry][0]
                                                   : mesh->edgeTriangles[entry][1];
    const int guard = mesh->vertexTriangles.vertexDegree(v) + 1;

    bool closed = false;
    while (cur >= 0) {
        fan.push_back(cur);

        int next = -1;
        for (int i = 0; i < 3; ++i) {
            const int e = mesh->triangleEdges[cur][i];
            if (e < 0 || e == entry) continue;
            if (mesh->edges[e][0] == v || mesh->edges[e][1] == v) { next = e; break; }
        }
        if (next < 0) break;
        if (mesh->isBoundaryEdge[next]) { closed = true; break; }

        const int nb = (mesh->edgeTriangles[next][0] == cur) ? mesh->edgeTriangles[next][1]
                                                             : mesh->edgeTriangles[next][0];
        if (nb < 0 || static_cast<int>(fan.size()) > guard) break;
        entry = next;
        cur = nb;
    }
    if (!closed) { fan.clear(); return fan; }

    // Orient counter-clockwise. In the first triangle the counter-clockwise
    // sweep at v runs from the edge to the next vertex towards the edge to the
    // previous one when the triangle is positively oriented; comparing that
    // against the edge the walk actually started on says which way it went.
    const Triangle &tri = mesh->triangles[fan.front()];
    int i0 = -1;
    for (int i = 0; i < 3; ++i) if (tri[i] == v) i0 = i;
    if (i0 < 0) { fan.clear(); return fan; }

    const int va = tri[(i0 + 1) % 3];
    const int vb = tri[(i0 + 2) % 3];
    const double s = cross2(mesh->vertices[va] - mesh->vertices[v],
                            mesh->vertices[vb] - mesh->vertices[v]);
    const int w = (mesh->edges[firstEdge][0] == v) ? mesh->edges[firstEdge][1]
                                                   : mesh->edges[firstEdge][0];
    const bool ccw = (w == va) ? (s > 0.0) : (s < 0.0);
    if (!ccw) std::reverse(fan.begin(), fan.end());

    return fan;
}

// ---------------------------------------------------------------------------
// computeBoundaryIndices()
//
// Pass 1 measures ( Theta + pi - Omega ) / (pi/2) at every boundary vertex.
// Pass 2 turns those into integers, and it does so along each boundary loop
// rather than one vertex at a time. A corner is a point defect of the
// continuum but a discrete cross field spreads it over however many vertices
// it takes to swing across, so rounding each vertex on its own can report a
// corner split 0.35 / 0.45 / 0.20 as no corner at all. Walking the loop with a
// running total and handing a unit to whoever carries it past the next halfway
// mark lands the cone on the middle of the corner and makes the loop's cone
// count equal the loop's unrounded total, rounded once. (This is the same
// argument, and the same fix, as UMBER::boundarySingularities().)
//
// The loop total is not an integer on its own -- only the grand total over
// every vertex of the mesh is, by Eq. (4) -- so whatever is left when the walk
// closes is handed to the last vertex rather than dropped.
// ---------------------------------------------------------------------------
void ConeSingularities::computeBoundaryIndices(const Eigen::VectorXcd &field) {
    const int nV = static_cast<int>(mesh->vertices.size());

    // --- Pass 1: the unrounded quarter turns ------------------------------
    for (int v : mesh->boundaryVertices) {
        const std::vector<int> fan = boundaryFan(v);
        if (fan.empty()) continue;

        double omega = 0.0;
        for (int t : fan) omega += cornerAngle(*mesh, t, v);

        double theta = 0.0;
        for (size_t i = 0; i + 1 < fan.size(); ++i) {
            theta += crossTurn(field[fan[i]], field[fan[i + 1]]);
        }

        rawIndex[v] = (theta + M_PI - omega) / M_PI_2;
        measured[v] = 1;
    }

    // --- Pass 2: hand the quarter turns out along each boundary loop ------
    std::vector<std::array<int, 2>> incident(nV, std::array<int, 2>{-1, -1});
    std::vector<int> degree(nV, 0);
    for (int e : mesh->boundaryEdges) {
        for (int k = 0; k < 2; ++k) {
            const int v = mesh->edges[e][k];
            if (degree[v] < 2) incident[v][degree[v]] = e;
            ++degree[v];
        }
    }

    auto otherEnd = [&](int e, int v) {
        return (mesh->edges[e][0] == v) ? mesh->edges[e][1] : mesh->edges[e][0];
    };

    std::vector<char> visited(nV, 0);
    for (int seed : mesh->boundaryVertices) {
        if (visited[seed] || degree[seed] != 2) continue;

        std::vector<int> loop;
        int v = seed;
        int prevEdge = -1;
        while (v >= 0 && !visited[v]) {
            visited[v] = 1;
            loop.push_back(v);

            int nextEdge = -1;
            for (int k = 0; k < 2; ++k) {
                const int e = incident[v][k];
                if (e >= 0 && e != prevEdge) { nextEdge = e; break; }
            }
            if (nextEdge < 0) break;
            prevEdge = nextEdge;
            v = otherEnd(nextEdge, v);
            if (degree[v] != 2) break;
        }

        // A running total, and whoever carries it past the next halfway mark
        // gets the cone. The last vertex takes whatever is left over, so the
        // loop hands out exactly round(sum of raw) however the corners were
        // smeared.
        double running = 0.0;
        long handedOut = 0;
        for (int u : loop) {
            if (measured[u]) running += rawIndex[u];
            const long want = std::lround(running);
            if (want != handedOut) {
                index[u] += static_cast<int>(want - handedOut);
                handedOut = want;
            }
        }
    }

    for (int v : mesh->boundaryVertices) {
        index[v] = std::max(minBoundaryIndex, std::min(maxBoundaryIndex, index[v]));
    }
}

// ---------------------------------------------------------------------------
// rebuildCones()
// ---------------------------------------------------------------------------
void ConeSingularities::rebuildCones() {
    const int nV = static_cast<int>(mesh->vertices.size());
    cones.clear();
    targetCurvature.assign(nV, 0.0);

    for (int v = 0; v < nV; ++v) {
        targetCurvature[v] = M_PI_2 * static_cast<double>(index[v]);
        if (index[v] == 0) continue;
        Cone c;
        c.vertex = v;
        c.index = index[v];
        c.onBoundary = mesh->isBoundaryVertex[v];
        c.valence = (c.onBoundary ? 3 : 4) - c.index;
        c.raw = rawIndex[v];
        c.fromRebalance = movedByRebalance[v] != 0;
        c.prescribed = !prescribed.empty() && prescribed[v] != 0;
        cones.push_back(c);
    }

    std::stable_sort(cones.begin(), cones.end(), [](const Cone &a, const Cone &b) {
        if (a.onBoundary != b.onBoundary) return !a.onBoundary;
        return a.vertex < b.vertex;
    });
}

std::vector<ConeSingularities::Cone> ConeSingularities::interiorCones() const {
    std::vector<Cone> out;
    for (const Cone &c : cones) if (!c.onBoundary) out.push_back(c);
    return out;
}

std::vector<ConeSingularities::Cone> ConeSingularities::boundaryCones() const {
    std::vector<Cone> out;
    for (const Cone &c : cones) if (c.onBoundary) out.push_back(c);
    return out;
}

// ---------------------------------------------------------------------------
// gaussBonnet()  --  Eqs. (4) and (7)
// ---------------------------------------------------------------------------
ConeSingularities::GaussBonnetReport ConeSingularities::gaussBonnet() const {
    GaussBonnetReport rep;
    const int nV = static_cast<int>(mesh->vertices.size());
    for (int v = 0; v < nV; ++v) (active[v] ? rep.V : rep.isolatedVertices) += 1;
    rep.E = static_cast<int>(mesh->edges.size());
    rep.F = static_cast<int>(mesh->triangles.size());
    rep.eulerCharacteristic = rep.V - rep.E + rep.F;

    // Boundary loops, for the report only: chi = 1 - (loops - 1) on a planar
    // mesh, so a mismatch here is a torn or non-manifold triangulation.
    {
        std::vector<int> degree(nV, 0);
        std::vector<std::array<int, 2>> incident(nV, std::array<int, 2>{-1, -1});
        for (int e : mesh->boundaryEdges) {
            for (int k = 0; k < 2; ++k) {
                const int v = mesh->edges[e][k];
                if (degree[v] < 2) incident[v][degree[v]] = e;
                ++degree[v];
            }
        }
        std::vector<char> seen(nV, 0);
        for (int seed : mesh->boundaryVertices) {
            if (seen[seed]) continue;
            ++rep.boundaryLoops;
            int v = seed, prevEdge = -1;
            while (v >= 0 && !seen[v]) {
                seen[v] = 1;
                int nextEdge = -1;
                for (int k = 0; k < 2; ++k) {
                    const int e = incident[v][k];
                    if (e >= 0 && e != prevEdge) { nextEdge = e; break; }
                }
                if (nextEdge < 0) break;
                prevEdge = nextEdge;
                v = (mesh->edges[nextEdge][0] == v) ? mesh->edges[nextEdge][1]
                                                    : mesh->edges[nextEdge][0];
            }
        }
    }

    for (int v = 0; v < nV; ++v) {
        if (!active[v]) continue;
        rep.indexSum += index[v];
        if (measured[v]) rep.rawIndexSum += rawIndex[v];
        rep.curvatureSum += inputCurvature[v];
        rep.targetCurvatureSum += targetCurvature[v];
    }
    rep.indexTarget = 4 * rep.eulerCharacteristic;
    rep.curvatureTarget = 2.0 * M_PI * rep.eulerCharacteristic;
    rep.curvatureResidual = std::fabs(rep.curvatureSum - rep.curvatureTarget);

    rep.rebalanceUnits = lastRebalanceUnits;
    rep.rebalanceCost = lastRebalanceCost;

    rep.admissible = (rep.indexSum == rep.indexTarget);
    rep.metricConsistent = rep.curvatureResidual < 1e-9 * std::max(1.0, std::fabs(rep.curvatureTarget));

    if (!rep.admissible) {
        std::ostringstream oss;
        oss << "Eq. (4) violated: sum I(v) = " << rep.indexSum << ", need 4 chi = "
            << rep.indexTarget << " (off by " << (rep.indexSum - rep.indexTarget) << ").";
        rep.messages.push_back(oss.str());
    }
    if (!rep.metricConsistent) {
        std::ostringstream oss;
        oss << "Eq. (7) violated: sum K_v = " << rep.curvatureSum << ", need 2 pi chi = "
            << rep.curvatureTarget << " (residual " << rep.curvatureResidual << ").";
        rep.messages.push_back(oss.str());
    }
    if (rep.isolatedVertices > 0) {
        std::ostringstream oss;
        oss << rep.isolatedVertices << " vertex/vertices in the file carry no triangle and "
            << "are excluded from chi; counting them would have given chi = "
            << (rep.eulerCharacteristic + rep.isolatedVertices) << " instead of "
            << rep.eulerCharacteristic << ".";
        rep.messages.push_back(oss.str());
    }
    if (rep.boundaryLoops > 0 && rep.eulerCharacteristic != 1 - (rep.boundaryLoops - 1)) {
        std::ostringstream oss;
        oss << "chi = " << rep.eulerCharacteristic << " does not match " << rep.boundaryLoops
            << " boundary loop(s) on a planar mesh; expected "
            << (1 - (rep.boundaryLoops - 1)) << ".";
        rep.messages.push_back(oss.str());
    }
    return rep;
}

// ---------------------------------------------------------------------------
// rebalance()  --  Sec. 3.1, the "manual adjustment" made automatic
//
// The imbalance can only come from the rounding on the boundary, so that is
// where it is paid back. Each unit goes to the boundary vertex where it costs
// the least, cost being the increase in |I(v) - raw(v)|: a vertex the field
// read as 0.51 gives its unit back for 0.02, one read as 1.00 would pay 2.00,
// so the units come off the corners the field was least sure about. That is
// the same currency Sec. 3.1 spends by hand ("adding or removing boundary
// cones") with the choice made by how weak the evidence for each cone was.
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// prescribe()
// ---------------------------------------------------------------------------
void ConeSingularities::prescribe(const std::vector<std::pair<int, int>> &fixed) {
    const int nV = static_cast<int>(mesh->vertices.size());
    if (static_cast<int>(prescribed.size()) != nV) prescribed.assign(nV, 0);
    lastPrescriptionShift = 0;

    for (const auto &p : fixed) {
        const int v = p.first;
        if (v < 0 || v >= nV || !active[v]) continue;
        lastPrescriptionShift += std::abs(p.second - index[v]);
        index[v] = p.second;
        // The unrounded value becomes the prescribed one: it is what rebalance()
        // costs a move against, and a prescribed index is not a rounding of
        // anything, so there is nothing for it to be a rounding of.
        rawIndex[v] = static_cast<double>(p.second);
        measured[v] = 1;
        prescribed[v] = 1;
        movedByRebalance[v] = 0;
    }
    rebuildCones();
}

// ---------------------------------------------------------------------------
// cancelDipoles()
// ---------------------------------------------------------------------------
int ConeSingularities::cancelDipoles(const std::vector<int> &group) {
    lastCancelledUnits = 0;
    const int nV = static_cast<int>(mesh->vertices.size());
    if (static_cast<int>(group.size()) < nV) return 0;

    auto movable = [&](int v) {
        return active[v] && !mesh->isBoundaryVertex[v] && group[v] >= 0 &&
               (prescribed.empty() || prescribed[v] == 0);
    };

    while (true) {
        // The closest pair of opposite sign sharing a group. Interior cones are
        // few -- eight is the corpus's largest -- so the quadratic scan is the
        // whole of the search.
        int bestA = -1, bestB = -1;
        double bestD = std::numeric_limits<double>::infinity();
        for (int a = 0; a < nV; ++a) {
            if (index[a] == 0 || !movable(a)) continue;
            for (int b = a + 1; b < nV; ++b) {
                if (index[b] == 0 || !movable(b)) continue;
                if (group[a] != group[b]) continue;
                if ((index[a] > 0) == (index[b] > 0)) continue;
                const double d = normP(mesh->vertices[a] - mesh->vertices[b]);
                if (d < bestD) { bestD = d; bestA = a; bestB = b; }
            }
        }
        if (bestA < 0) break;

        const int units = std::min(std::abs(index[bestA]), std::abs(index[bestB]));
        index[bestA] -= (index[bestA] > 0 ? units : -units);
        index[bestB] -= (index[bestB] > 0 ? units : -units);
        rawIndex[bestA] = index[bestA];
        rawIndex[bestB] = index[bestB];
        lastCancelledUnits += units;
    }

    if (lastCancelledUnits > 0) rebuildCones();
    return lastCancelledUnits;
}

// ---------------------------------------------------------------------------
int ConeSingularities::rebalance() {
    GaussBonnetReport rep = gaussBonnet();
    if (rep.admissible) { lastRebalanceUnits = 0; lastRebalanceCost = 0.0; return 0; }

    int residual = rep.indexTarget - rep.indexSum; // units still to place
    const int step = (residual > 0) ? 1 : -1;
    int moved = 0;
    double cost = 0.0;

    while (residual != 0) {
        int best = -1;
        double bestDelta = std::numeric_limits<double>::infinity();

        for (int v : mesh->boundaryVertices) {
            if (!prescribed.empty() && prescribed[v]) continue;   // fixed by geometry
            const int next = index[v] + step;
            if (next < minBoundaryIndex || next > maxBoundaryIndex) continue;
            const double before = std::fabs(static_cast<double>(index[v]) - rawIndex[v]);
            const double after = std::fabs(static_cast<double>(next) - rawIndex[v]);
            const double delta = after - before;
            if (delta < bestDelta) { bestDelta = delta; best = v; }
        }

        if (best < 0) {
            lastRebalanceUnits = moved;
            lastRebalanceCost = cost;
            rebuildCones();
            return -1; // nowhere left to put it
        }

        index[best] += step;
        movedByRebalance[best] = 1;
        cost += bestDelta;
        ++moved;
        residual -= step;
    }

    lastRebalanceUnits = moved;
    lastRebalanceCost = cost;
    rebuildCones();
    return moved;
}

// ---------------------------------------------------------------------------
// relocateFlatCones()
//
// See the header for why. Two notes on the how.
//
// "Along the same side" is a walk over dS that stops at the first vertex that
// is a corner (the 20 degree rule of BoundaryFeatures) or carries any index or
// was prescribed. The unit only goes to that vertex when it is a convex corner
// with no cone of its own: never past it, so a quarter is never carried round
// a corner that already turns, and never onto a reflex corner, whose -1 the
// unit would cancel into two elements of a half turn apiece.
//
// "Inside" is the vertex furthest from every feature -- dS, or an edge between
// two vertices no group claims, which is the interface network on a model
// that has one -- among those reachable from the unit's vertex within a
// quarter of the length of its side. In a wedge that distance grows along the
// bisector, in a triangle it peaks at the incentre, and in a layer that
// pinches out it stops growing where the layer reaches its full thickness;
// those are the three places the extra quarter belongs. The quarter is so
// that two units leaving the two ends of one side -- a lens has a pinch-out at
// each -- do not both go to its middle. The search walks through dS as well as
// the inside, because a wedge a few triangles across has interior vertices no
// interior path joins, and it never ends on a cone or within Sec. 3.1's
// clustering distance of one: a cone that close to another is merged with it
// by Stage 6 whatever Stage 1 says (singlemat/geom033, 1e-4 of the model
// apart, came out of Stage 6 a half turn off).
// ---------------------------------------------------------------------------
int ConeSingularities::relocateFlatCones(const std::vector<int> &group, double maxCornerAngle) {
    relocations.clear();
    const int nV = static_cast<int>(mesh->vertices.size());
    const bool grouped = static_cast<int>(group.size()) >= nV;
    auto groupOf = [&](int v) { return grouped ? group[v] : 0; };
    auto isPrescribed = [&](int v) { return !prescribed.empty() && prescribed[v] != 0; };
    // On dS, inputCurvature is the turning pi - Omega.
    auto omega = [&](int v) { return M_PI - inputCurvature[v]; };
    const double cornerTurn = 20.0 * M_PI / 180.0;

    std::vector<std::array<int, 2>> incident(nV, std::array<int, 2>{-1, -1});
    std::vector<int> degree(nV, 0);
    for (int e : mesh->boundaryEdges) {
        for (int k = 0; k < 2; ++k) {
            const int v = mesh->edges[e][k];
            if (degree[v] < 2) incident[v][degree[v]] = e;
            ++degree[v];
        }
    }
    auto otherEnd = [&](int e, int v) {
        return (mesh->edges[e][0] == v) ? mesh->edges[e][1] : mesh->edges[e][0];
    };
    auto stopsWalk = [&](int u) {
        return degree[u] != 2 || index[u] != 0 || isPrescribed(u) ||
               std::fabs(omega(u) - M_PI) > cornerTurn;
    };
    // From v along `first`, to the first vertex that stops a walk; -1 if the
    // loop closes first. `length` is the arc length walked either way.
    auto walk = [&](int v, int first, double &length) -> int {
        length = 0.0;
        int prevEdge = first;
        int u = otherEnd(first, v);
        length += normP(mesh->vertices[u] - mesh->vertices[v]);
        for (int guard = 0; u != v && guard < nV; ++guard) {
            if (stopsWalk(u)) return u;
            const int e = (incident[u][0] == prevEdge) ? incident[u][1] : incident[u][0];
            if (e < 0) return -1;
            const int w = otherEnd(e, u);
            length += normP(mesh->vertices[w] - mesh->vertices[u]);
            prevEdge = e;
            u = w;
        }
        return -1;
    };

    // The features a cone inside keeps its distance from: dS, and on a model
    // with an interface network the edges between two vertices no group
    // claims. Built on first use, with the vertex adjacency.
    std::vector<std::vector<int>> nbr;
    std::vector<std::array<int, 2>> features;
    auto prepare = [&]() {
        nbr.assign(nV, std::vector<int>());
        for (int e = 0; e < static_cast<int>(mesh->edges.size()); ++e) {
            const Edge &E = mesh->edges[e];
            nbr[E[0]].push_back(E[1]);
            nbr[E[1]].push_back(E[0]);
            if (mesh->isBoundaryEdge[e] || (grouped && group[E[0]] < 0 && group[E[1]] < 0))
                features.push_back({E[0], E[1]});
        }
    };
    auto depthOf = [&](int u) {
        const Point &p = mesh->vertices[u];
        double best = std::numeric_limits<double>::infinity();
        for (const std::array<int, 2> &f : features) {
            const Point &a = mesh->vertices[f[0]];
            const Point d = mesh->vertices[f[1]] - a;
            const double L2 = dotP(d, d);
            double t = (L2 > 0.0) ? dotP(p - a, d) / L2 : 0.0;
            t = std::max(0.0, std::min(1.0, t));
            best = std::min(best, normP(p - (a + d * t)));
        }
        return best;
    };
    // How close an inside destination may come to any cone: Sec. 3.1's
    // clustering distance, which Stage 4 warns at and Stage 6 cannot keep two
    // cones apart within, as a fraction of this mesh's own extent.
    double clearance = 0.0;
    {
        Point lo = mesh->vertices.empty() ? Point{0.0, 0.0} : mesh->vertices[0], hi = lo;
        for (int v = 0; v < nV; ++v) {
            if (!active[v]) continue;
            for (int k = 0; k < 2; ++k) {
                lo[k] = std::min(lo[k], mesh->vertices[v][k]);
                hi[k] = std::max(hi[k], mesh->vertices[v][k]);
            }
        }
        clearance = 0.01 * normP(hi - lo);
    }
    auto nearCone = [&](int u) {
        for (int w = 0; w < nV; ++w) {
            if (index[w] != 0 && normP(mesh->vertices[w] - mesh->vertices[u]) < clearance) return true;
        }
        return false;
    };
    // The deepest vertex of v's group within `radius` of it, through the
    // group, that carries no cone and keeps the clearance from every cone; -1
    // if there is none.
    auto deepestInside = [&](int v, double radius) -> int {
        const int g = groupOf(v);
        if (g < 0) return -1;
        std::vector<double> reach;  // sparse would do; the regions are small
        reach.assign(nV, std::numeric_limits<double>::infinity());
        typedef std::pair<double, int> Item;
        std::priority_queue<Item, std::vector<Item>, std::greater<Item>> open;
        reach[v] = 0.0;
        open.push(Item(0.0, v));
        int best = -1;
        double bestDepth = 0.0;
        while (!open.empty()) {
            const Item it = open.top();
            open.pop();
            if (it.first > reach[it.second]) continue;
            for (int w : nbr[it.second]) {
                // Never across a vertex of another group or of none.
                if (!active[w] || groupOf(w) != g) continue;
                const double d = it.first + normP(mesh->vertices[w] - mesh->vertices[it.second]);
                if (d > radius || d >= reach[w]) continue;
                reach[w] = d;
                open.push(Item(d, w));
                if (mesh->isBoundaryVertex[w]) continue;
                if (index[w] != 0 || isPrescribed(w)) continue;
                const double depth = depthOf(w);
                if (!(depth > bestDepth) || nearCone(w)) continue;
                bestDepth = depth;
                best = w;
            }
        }
        return best;
    };
    for (int v : mesh->boundaryVertices) {
        if (!active[v] || index[v] <= 0 || isPrescribed(v) || degree[v] != 2) continue;
        if (omega(v) <= maxCornerAngle) continue;

        // The corner at the far end of the side, either way: the nearer convex
        // one with room.
        double length[2] = {0.0, 0.0};
        int corner = -1;
        double cornerLength = std::numeric_limits<double>::infinity();
        for (int k = 0; k < 2; ++k) {
            const int u = walk(v, incident[v][k], length[k]);
            if (u < 0 || u == v || isPrescribed(u) || index[u] != 0 || degree[u] != 2) continue;
            if (omega(u) > maxCornerAngle || omega(u) > M_PI - cornerTurn) continue;
            if (length[k] < cornerLength) { cornerLength = length[k]; corner = u; }
        }

        // Otherwise inside, no further from it than a quarter of its side.
        int to = corner;
        const bool in = (corner < 0);
        if (in) {
            if (nbr.empty()) prepare();
            index[v] -= 1;   // lifted first, so the search does not see it as a cone
            to = deepestInside(v, 0.25 * (length[0] + length[1]));
            index[v] += 1;
        }
        if (to < 0) continue;

        index[v] -= 1;
        index[to] += 1;
        movedByRebalance[to] = 1;
        Relocation r;
        r.from = v;
        r.to = to;
        r.inside = in;
        r.distance = normP(mesh->vertices[to] - mesh->vertices[v]);
        relocations.push_back(r);
    }

    if (!relocations.empty()) rebuildCones();
    return static_cast<int>(relocations.size());
}
