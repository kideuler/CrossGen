#include "ConeSingularities.hxx"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "SIPG/SIPG.hxx"

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
ConeSingularities::ConeSingularities(SIPG &sipg)
    : ConeSingularities(sipg.getMeshPtr(), sipg.u_k_prev) {
    // Populate the solver's own singularity list from the same field, so that
    // anything else reading SIPG::singularVertices later agrees with the cones
    // prescribed here.
    sipg.computeSingularities();
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
// The same winding number SIPG::computeSingularities() reports, kept here in
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
