#include "Metrics.hxx"

#include "dualmbo/DualMBO.hxx"   // kMinDualDistance, the two-point weight's floor

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <queue>
#include <sstream>
#include <unordered_map>

namespace paper {
namespace metrics {
namespace {

double triangleArea(const Mesh &m, int t) {
    const Triangle &tri = m.triangles[t];
    const Point &p0 = m.vertices[tri[0]];
    const Point a = m.vertices[tri[1]] - p0;
    const Point b = m.vertices[tri[2]] - p0;
    return 0.5 * std::fabs(cross2(a, b));
}

double edgeLength(const Mesh &m, int e) {
    return normP(m.vertices[m.edges[e][1]] - m.vertices[m.edges[e][0]]);
}

double angleDiff(const std::complex<double> &from, const std::complex<double> &to) {
    const std::complex<double> d = to * std::conj(from);
    return std::atan2(d.imag(), d.real());
}

} // namespace

// ---------------------------------------------------------------------------
EdgeWeights edgeWeights(const Mesh &m, double gamma) {
    EdgeWeights w;
    w.gamma = gamma;
    w.area.resize(m.triangles.size());
    for (std::size_t t = 0; t < m.triangles.size(); ++t) w.area[t] = triangleArea(m, static_cast<int>(t));

    w.kappa.assign(m.edges.size(), 0.0);
    for (std::size_t e = 0; e < m.edges.size(); ++e) {
        if (m.isBoundaryEdge[e]) continue;
        const double len = edgeLength(m, static_cast<int>(e));
        if (len < 1e-14) continue;
        const int ti = m.edgeTriangles[e][0];
        const int tj = m.edgeTriangles[e][1];
        if (ti < 0 || tj < 0) continue;
        const double he = std::min(2.0 * w.area[ti] / len, 2.0 * w.area[tj] / len);
        if (he < 1e-16) continue;
        w.kappa[e] = gamma * len / he;
    }
    return w;
}

EdgeWeights edgeWeightsTwoPoint(const Mesh &m, double gamma) {
    EdgeWeights w;
    w.gamma = gamma;
    const int NT = static_cast<int>(m.triangles.size());
    w.area.resize(NT);
    std::vector<Point> circ(NT), cent(NT);
    for (int t = 0; t < NT; ++t) {
        w.area[t] = triangleArea(m, t);
        const Triangle &tri = m.triangles[t];
        const Point &pa = m.vertices[tri[0]];
        const Point &pb = m.vertices[tri[1]];
        const Point &pc = m.vertices[tri[2]];
        const double bx = pb[0] - pa[0], by = pb[1] - pa[1];
        const double cx = pc[0] - pa[0], cy = pc[1] - pa[1];
        const double den = 2.0 * (bx * cy - by * cx);
        const double b2 = bx * bx + by * by, c2 = cx * cx + cy * cy;
        Point cc = pa;
        if (std::fabs(den) > 1e-300) {
            cc[0] = pa[0] + (cy * b2 - by * c2) / den;
            cc[1] = pa[1] + (bx * c2 - cx * b2) / den;
        }
        circ[t] = cc;
        cent[t][0] = (pa[0] + pb[0] + pc[0]) / 3.0;
        cent[t][1] = (pa[1] + pb[1] + pc[1]) / 3.0;
    }

    w.kappa.assign(m.edges.size(), 0.0);
    for (std::size_t e = 0; e < m.edges.size(); ++e) {
        if (m.isBoundaryEdge[e]) continue;
        const Point &ea = m.vertices[m.edges[e][0]];
        const Point &eb = m.vertices[m.edges[e][1]];
        const double len = std::hypot(eb[0] - ea[0], eb[1] - ea[1]);
        if (len < 1e-14) continue;
        const int ti = m.edgeTriangles[e][0];
        const int tj = m.edgeTriangles[e][1];
        if (ti < 0 || tj < 0) continue;
        const double hi = 2.0 * w.area[ti] / len, hj = 2.0 * w.area[tj] / len;
        double nx = -(eb[1] - ea[1]) / len, ny = (eb[0] - ea[0]) / len;
        if ((cent[tj][0] - cent[ti][0]) * nx + (cent[tj][1] - cent[ti][1]) * ny < 0.0) {
            nx = -nx; ny = -ny;
        }
        const double d = (circ[tj][0] - circ[ti][0]) * nx + (circ[tj][1] - circ[ti][1]) * ny;
        const double dFloor = DualMBO::kMinDualDistance * 0.5 * (hi + hj);
        w.kappa[e] = (2.0 / 3.0) * gamma * len / std::max(d, dFloor);
    }
    return w;
}

double boundaryEnergy(const Mesh &m, double gamma, const Eigen::VectorXcd &u,
                      const std::vector<char> &alignedEdge) {
    std::vector<double> area(m.triangles.size());
    for (std::size_t t = 0; t < m.triangles.size(); ++t) area[t] = triangleArea(m, static_cast<int>(t));
    auto aligned = [&](std::size_t e) { return e < alignedEdge.size() && alignedEdge[e]; };

    double E = 0.0;
    for (std::size_t e = 0; e < m.edges.size(); ++e) {
        if (!m.isBoundaryEdge[e] && !aligned(e)) continue;
        const Point d = m.vertices[m.edges[e][1]] - m.vertices[m.edges[e][0]];
        const double len = normP(d);
        if (len < 1e-14) continue;
        const std::complex<double> ge = std::exp(std::complex<double>(0.0, 4.0 * std::atan2(d[1], d[0])));
        for (int side = 0; side < 2; ++side) {
            const int t = m.edgeTriangles[e][side];
            if (t < 0) continue;
            if (m.isBoundaryEdge[e] && side == 1) continue;
            const double h1 = 2.0 * area[t] / len;
            if (h1 < 1e-16) continue;
            E += 0.5 * (gamma * len / h1) * std::norm(u[t] - ge);
        }
    }
    return E;
}

double commonEnergy(const Mesh &m, const EdgeWeights &w, const Eigen::VectorXcd &u,
                    const std::vector<char> &skipEdge) {
    double E = 0.0;
    for (std::size_t e = 0; e < m.edges.size(); ++e) {
        if (w.kappa[e] <= 0.0) continue;
        if (!skipEdge.empty() && skipEdge[e]) continue;
        const int ti = m.edgeTriangles[e][0];
        const int tj = m.edgeTriangles[e][1];
        const std::complex<double> d = u[ti] - u[tj];
        E += w.kappa[e] * std::norm(d);
    }
    return 0.5 * E;
}

double dualMBOEnergy(const Mesh &m, double gamma, const Eigen::VectorXcd &u,
                  const std::vector<char> &alignedEdge) {
    std::vector<double> area(m.triangles.size());
    for (std::size_t t = 0; t < m.triangles.size(); ++t) area[t] = triangleArea(m, static_cast<int>(t));

    auto aligned = [&](std::size_t e) {
        return e < alignedEdge.size() && alignedEdge[e];
    };

    double E = 0.0;
    for (std::size_t e = 0; e < m.edges.size(); ++e) {
        const Point d = m.vertices[m.edges[e][1]] - m.vertices[m.edges[e][0]];
        const double len = normP(d);
        if (len < 1e-14) continue;
        const std::complex<double> ge = std::exp(std::complex<double>(0.0, 4.0 * std::atan2(d[1], d[0])));

        const int ti = m.edgeTriangles[e][0];
        const int tj = m.edgeTriangles[e][1];

        if (!m.isBoundaryEdge[e] && !aligned(e)) {
            if (ti < 0 || tj < 0) continue;
            const double he = std::min(2.0 * area[ti] / len, 2.0 * area[tj] / len);
            if (he < 1e-16) continue;
            E += 0.5 * (gamma * len / he) * std::norm(u[ti] - u[tj]);
            continue;
        }

        // One-sided: dS has one incident triangle, an aligned interior edge two,
        // and each is pinned to the same tangent from its own side.
        for (int side = 0; side < 2; ++side) {
            const int t = m.edgeTriangles[e][side];
            if (t < 0) continue;
            if (m.isBoundaryEdge[e] && side == 1) continue;
            const double h1 = 2.0 * area[t] / len;
            if (h1 < 1e-16) continue;
            E += 0.5 * (gamma * len / h1) * std::norm(u[t] - ge);
        }
    }
    return E;
}

// ---------------------------------------------------------------------------
std::vector<double> localEdgeLength(const Mesh &m) {
    std::vector<double> len(m.vertices.size(), 0.0);
    std::vector<int> cnt(m.vertices.size(), 0);
    for (std::size_t e = 0; e < m.edges.size(); ++e) {
        const double l = edgeLength(m, static_cast<int>(e));
        len[m.edges[e][0]] += l; ++cnt[m.edges[e][0]];
        len[m.edges[e][1]] += l; ++cnt[m.edges[e][1]];
    }
    for (std::size_t v = 0; v < len.size(); ++v) if (cnt[v]) len[v] /= cnt[v];
    return len;
}

int eulerCharacteristic(const Mesh &m) {
    // Only vertices some triangle uses. Several models in data/meshes carry
    // vertices no face references -- geom031 has 309 of them -- and each one
    // adds 1 to V and nothing to E or F, so counting them puts chi out by the
    // same amount and every Poincare-Hopf statement with it.
    std::vector<char> used(m.vertices.size(), 0);
    for (const Triangle &t : m.triangles) { used[t[0]] = 1; used[t[1]] = 1; used[t[2]] = 1; }
    int V = 0;
    for (char c : used) V += c;
    return V - static_cast<int>(m.edges.size()) + static_cast<int>(m.triangles.size());
}

std::vector<double> interiorAngles(const Mesh &m) {
    std::vector<double> a(m.vertices.size(), 0.0);
    for (const Triangle &tri : m.triangles) {
        for (int k = 0; k < 3; ++k) {
            const int v = tri[k];
            const Point &o = m.vertices[v];
            const Point p = m.vertices[tri[(k + 1) % 3]] - o;
            const Point q = m.vertices[tri[(k + 2) % 3]] - o;
            const double np = normP(p), nq = normP(q);
            if (np < 1e-15 || nq < 1e-15) continue;
            a[v] += std::acos(std::max(-1.0, std::min(1.0, dotP(p, q) / (np * nq))));
        }
    }
    return a;
}

std::vector<int> cornerQuarters(const Mesh &m) {
    const std::vector<double> a = interiorAngles(m);
    std::vector<int> q(m.vertices.size(), 0);
    for (std::size_t v = 0; v < m.vertices.size(); ++v) {
        if (!m.isBoundaryVertex[v]) continue;
        q[v] = static_cast<int>(std::lround(2.0 * (M_PI - a[v]) / M_PI));
    }
    return q;
}

// ---------------------------------------------------------------------------
SingularityReport singularities(const Mesh &m, const Eigen::VectorXcd &u) {
    SingularityReport r;
    const int nV = static_cast<int>(m.vertices.size());
    const auto &vt = m.vertexTriangles;
    const std::vector<double> hLocal = localEdgeLength(m);

    // Interior vertices: the winding of the spin-4 value around the star,
    // exactly as DualMBO::computeSingularities reads it.
    for (int v = 0; v < nV; ++v) {
        const int begin = vt.rowPtr[v], end = vt.rowPtr[v + 1];
        const int star = end - begin;
        if (star < 2) continue;

        double total = 0.0, worst = 0.0;
        for (int k = begin; k < end; ++k) {
            const int cur = vt.colIdx[k];
            const int nxt = vt.colIdx[(k - begin + 1) % star + begin];
            const double d = angleDiff(u[cur], u[nxt]);
            total += d;
            worst = std::max(worst, std::fabs(d));
        }

        if (!m.isBoundaryVertex[v]) {
            r.maxSpinJump = std::max(r.maxSpinJump, worst);
            if (worst >= M_PI_4) ++r.starsPastCondition;

            const int winding = static_cast<int>(std::lround(total / (2.0 * M_PI)));
            if (winding != 0) {
                Singularity s;
                s.vertex = v;
                s.index4 = winding;
                s.localEdgeLength = hLocal[v];
                r.interior.push_back(s);
                r.interiorIndex4Sum += winding;
            }
        }
    }

    // --- the boundary index, read off the field ---------------------------
    //
    // A boundary vertex's fan is open, and closing it is the whole of Prop. 3's
    // boundary term. Two pieces:
    //
    //   * the two rays. The field is meant to be tangent to dS, so the jump
    //     from the first boundary ray into the first face -- and out of the last
    //     face into the last ray -- is how far it is from being so. Wrapped,
    //     like every other jump.
    //   * the exterior wedge. There is no face there, and what a
    //     boundary-aligned field does across it is turn with the tangent: by
    //     4 tau in spin space, tau = pi - a the exterior turn at v. This term is
    //     *not* wrapped, and that is exactly the point -- it is where the
    //     quarter turns a corner takes come from, and wrapping it is what would
    //     throw them away.
    //
    // The rounded total is the index. On a field that follows the boundary it
    // reduces to round(2 tau / pi), the geometric quarter count; where it does
    // not, the difference is a cone the field put on dS.
    const std::vector<int> q = cornerQuarters(m);
    const std::vector<double> ang = interiorAngles(m);
    {
        // The spin-4 value of the ray from v to w.
        auto beta = [&](int v, int w) {
            const Point d = m.vertices[w] - m.vertices[v];
            return 4.0 * std::atan2(d[1], d[0]);
        };

        // Where v sits in a triangle, and the two rays out of it there. A CCW
        // triangle listed (v, p, q) sweeps CCW from the ray to p to the ray to
        // q, so p is the fan's entry ray in that triangle and q its exit ray.
        auto cornerOf = [&](int t, int v) {
            for (int k = 0; k < 3; ++k) if (m.triangles[t][k] == v) return k;
            return -1;
        };

        for (int v = 0; v < nV; ++v) {
            if (!m.isBoundaryVertex[v]) continue;
            const int begin = vt.rowPtr[v], end = vt.rowPtr[v + 1];
            const int star = end - begin;
            if (star < 1) continue;

            // The fan is walked by edge adjacency rather than read off the CSR:
            // the CSR sorts the star by centroid angle about a branch cut that
            // has nothing to do with where the boundary is, and on a star of two
            // or three faces there is no cyclic "gap" to find -- the first and
            // last faces of an open fan can perfectly well share an edge with
            // each other's neighbours. Walking is exact and needs no tolerance.
            //
            // Start at the face whose entry ray is on dS, step out through each
            // face's exit ray, and stop when the exit ray is on dS too.
            int t0 = -1, wStart = -1;
            for (int k = begin; k < end && t0 < 0; ++k) {
                const int t = vt.colIdx[k];
                const int c = cornerOf(t, v);
                if (c < 0) continue;
                if (m.isBoundaryEdge[m.triangleEdges[t][c]]) {
                    t0 = t;
                    wStart = m.triangles[t][(c + 1) % 3];
                }
            }
            if (t0 < 0) { ++r.boundaryUnclosed; continue; }

            std::vector<int> chain;
            int wEnd = -1;
            {
                int t = t0;
                for (int guard = 0; guard <= star + 1; ++guard) {
                    chain.push_back(t);
                    const int c = cornerOf(t, v);
                    if (c < 0) { wEnd = -1; break; }
                    const int exitEdge = m.triangleEdges[t][(c + 2) % 3];
                    if (m.isBoundaryEdge[exitEdge]) {
                        wEnd = m.triangles[t][(c + 2) % 3];
                        break;
                    }
                    const int a = m.edgeTriangles[exitEdge][0];
                    const int b = m.edgeTriangles[exitEdge][1];
                    const int nxt = (a == t) ? b : a;
                    if (nxt < 0) { wEnd = -1; break; }
                    t = nxt;
                }
            }
            if (wEnd < 0 || static_cast<int>(chain.size()) != star) {
                // Either the walk left the fan or the fan is not the whole star:
                // a pinched boundary vertex, where dS visits v twice. There is
                // no single index to read there.
                ++r.boundaryUnclosed;
                continue;
            }

            double total = wrap_pi(std::arg(u[chain.front()]) - beta(v, wStart));
            for (int i = 0; i + 1 < star; ++i) total += angleDiff(u[chain[i]], u[chain[i + 1]]);
            total += wrap_pi(beta(v, wEnd) - std::arg(u[chain.back()]));
            total += 4.0 * (M_PI - ang[v]);

            // The total is an exact multiple of 2 pi when the extraction is
            // sound: the wrapped jumps telescope to 4 a_v modulo 2 pi and the
            // wedge term adds 4 tau_v = 4 pi - 4 a_v. A vertex where it is not
            // is a vertex where a wrap was lost, and it is counted rather than
            // rounded away.
            const double turns = total / (2.0 * M_PI);
            if (std::fabs(turns - std::round(turns)) > 0.05) ++r.boundaryNonIntegral;

            const int index4 = static_cast<int>(std::lround(total / (2.0 * M_PI)));
            r.boundaryIndex4Sum += index4;
            if (index4 != 0) {
                Singularity sg;
                sg.vertex = v;
                sg.index4 = index4;
                sg.onBoundary = true;
                sg.localEdgeLength = hLocal[v];
                r.boundary.push_back(sg);
            }
            if (index4 != q[v]) ++r.boundaryAnomalies;

            const double quarters = 2.0 * (M_PI - ang[v]) / M_PI;
            if (std::fabs(quarters - std::round(quarters)) > 0.4) ++r.ambiguousCorners;
        }
    }

    for (int v = 0; v < nV; ++v) if (m.isBoundaryVertex[v]) r.boundaryQuarterSum += q[v];

    r.chi = eulerCharacteristic(m);
    r.poincareHopfResidual = r.interiorIndex4Sum + r.boundaryIndex4Sum - 4 * r.chi;
    r.geometricPoincareHopfResidual =
        r.interiorIndex4Sum + r.boundaryQuarterSum - 4 * r.chi;

    // Distance to dS, and adjacency to it.
    if (!r.interior.empty()) {
        r.minBoundaryDistanceH = std::numeric_limits<double>::max();
        std::vector<std::vector<int>> nbr(nV);
        for (const auto &e : m.edges) {
            nbr[e[0]].push_back(e[1]);
            nbr[e[1]].push_back(e[0]);
        }
        for (Singularity &s : r.interior) {
            double best = std::numeric_limits<double>::max();
            for (int bv : m.boundaryVertices)
                best = std::min(best, normP(m.vertices[bv] - m.vertices[s.vertex]));
            s.distanceToBoundary = best;
            for (int w : nbr[s.vertex]) if (m.isBoundaryVertex[w]) { ++r.boundaryAdjacent; break; }
            if (s.localEdgeLength > 0.0)
                r.minBoundaryDistanceH = std::min(r.minBoundaryDistanceH, best / s.localEdgeLength);
        }
        if (r.minBoundaryDistanceH == std::numeric_limits<double>::max()) r.minBoundaryDistanceH = 0.0;

        for (std::size_t i = 0; i < r.interior.size(); ++i) {
            for (std::size_t j = i + 1; j < r.interior.size(); ++j) {
                const double d = normP(m.vertices[r.interior[i].vertex] - m.vertices[r.interior[j].vertex]);
                const double h = 0.5 * (r.interior[i].localEdgeLength + r.interior[j].localEdgeLength);
                if (h > 0.0 && d < 2.0 * h) {
                    ++r.closePairs;
                    if (r.interior[i].index4 * r.interior[j].index4 < 0) ++r.closeOppositePairs;
                }
            }
        }
    }
    return r;
}

std::string SingularityReport::indexMultiset() const {
    std::map<int, int> count;
    for (const Singularity &s : interior) ++count[s.index4];
    if (count.empty()) return "(none)";
    std::ostringstream oss;
    bool first = true;
    for (auto it = count.rbegin(); it != count.rend(); ++it) {
        if (!first) oss << " ";
        first = false;
        oss << (it->first > 0 ? "+" : "") << it->first << "/4";
        if (it->second > 1) oss << "x" << it->second;
    }
    return oss.str();
}

// ---------------------------------------------------------------------------
namespace {

AlignmentError summarise(std::vector<double> &a) {
    AlignmentError e;
    e.edges = static_cast<int>(a.size());
    if (a.empty()) return e;
    std::sort(a.begin(), a.end());
    e.maxAngle = a.back();
    e.p95Angle = a[std::min(a.size() - 1, static_cast<std::size_t>(0.95 * a.size()))];
    double sum = 0.0;
    for (double x : a) sum += x;
    e.meanAngle = sum / a.size();
    return e;
}

} // namespace

AlignmentError alignmentError(const Mesh &m, const Eigen::VectorXcd &u,
                              const std::vector<int> &which) {
    std::vector<double> a;
    a.reserve(which.size() * 2);
    for (int e : which) {
        if (e < 0 || e >= static_cast<int>(m.edges.size())) continue;
        const Point d = m.vertices[m.edges[e][1]] - m.vertices[m.edges[e][0]];
        if (normP(d) < 1e-14) continue;
        const double theta = std::atan2(d[1], d[0]);
        for (int side = 0; side < 2; ++side) {
            const int t = m.edgeTriangles[e][side];
            if (t < 0) continue;
            a.push_back(crossToDirection(u[t], theta));
        }
    }
    return summarise(a);
}

AlignmentError boundaryAlignmentError(const Mesh &m, const Eigen::VectorXcd &u) {
    return alignmentError(m, u, m.boundaryEdges);
}

// ---------------------------------------------------------------------------
QuadMetrics quadMetrics(const std::vector<Point> &verts,
                        const std::vector<std::array<int, 4>> &quads,
                        const std::vector<int> &quadMaterial,
                        double targetEdge) {
    QuadMetrics q;
    q.vertices = static_cast<int>(verts.size());
    q.quads = static_cast<int>(quads.size());
    if (quads.empty()) return q;

    int maxMat = 1;
    for (int mid : quadMaterial) maxMat = std::max(maxMat, mid);
    q.materials = quadMaterial.empty() ? 1 : maxMat;
    q.irregularPerMaterial.assign(q.materials, 0);

    // Edge use count -> which vertices are on the mesh boundary.
    std::map<std::pair<int, int>, int> useCount;
    std::vector<int> valence(verts.size(), 0);
    std::vector<char> onBoundary(verts.size(), 0);
    for (const auto &c : quads) {
        for (int k = 0; k < 4; ++k) {
            ++valence[c[k]];
            const int a = c[k], b = c[(k + 1) % 4];
            ++useCount[{std::min(a, b), std::max(a, b)}];
        }
    }
    for (const auto &[e, n] : useCount) {
        if (n == 1) { onBoundary[e.first] = 1; onBoundary[e.second] = 1; }
    }

    // The material a vertex belongs to, for the per-region irregular count: the
    // material of its incident quads when they agree, and "on an interface"
    // when they do not. An irregular vertex sitting on the interface is charged
    // to neither region -- it is R1's business, not R4's.
    std::vector<int> vertMat(verts.size(), 0);
    std::vector<char> vertMixed(verts.size(), 0);
    for (std::size_t c = 0; c < quads.size(); ++c) {
        const int mid = quadMaterial.empty() ? 1 : quadMaterial[c];
        for (int k = 0; k < 4; ++k) {
            const int v = quads[c][k];
            if (vertMat[v] == 0) vertMat[v] = mid;
            else if (vertMat[v] != mid) vertMixed[v] = 1;
        }
    }

    for (std::size_t v = 0; v < verts.size(); ++v) {
        if (valence[v] == 0) continue;
        if (onBoundary[v]) {
            if (valence[v] != 2) ++q.irregularBoundary;
        } else if (valence[v] != 4) {
            ++q.irregularInterior;
            if (!vertMixed[v] && vertMat[v] >= 1 && vertMat[v] <= q.materials)
                ++q.irregularPerMaterial[vertMat[v] - 1];
        }
    }

    // R1: a quad whose corners do not all carry the same material would be a
    // mixed zone. On a mesh built from an interface-aligned layout every quad's
    // four corners are in one region, so this counts quads whose *corner*
    // materials disagree -- the only signal available from the element table
    // alone.
    if (!quadMaterial.empty()) {
        for (std::size_t c = 0; c < quads.size(); ++c) {
            bool mixed = false;
            for (int k = 0; k < 4 && !mixed; ++k) {
                const int v = quads[c][k];
                if (vertMixed[v]) continue;             // an interface node, fine
                if (vertMat[v] != quadMaterial[c]) mixed = true;
            }
            if (mixed) ++q.mixedQuads;
        }
    }

    // Quality: the scaled Jacobian at each corner, the angle deviation, and the
    // characteristic zone length.
    double sjSum = 0.0, devSum = 0.0, zoneSum = 0.0;
    q.minScaledJacobian = std::numeric_limits<double>::max();
    q.minZoneDimension = std::numeric_limits<double>::max();
    std::vector<double> zones;
    zones.reserve(quads.size());
    const double scale = (targetEdge > 0.0) ? targetEdge : 1.0;

    for (const auto &c : quads) {
        double worst = std::numeric_limits<double>::max();
        double area = 0.0;
        for (int k = 0; k < 4; ++k) {
            const Point &p = verts[c[k]];
            const Point &pn = verts[c[(k + 1) % 4]];
            area += 0.5 * (p[0] * pn[1] - pn[0] * p[1]);

            const Point a = verts[c[(k + 1) % 4]] - verts[c[k]];
            const Point b = verts[c[(k + 3) % 4]] - verts[c[k]];
            const double na = normP(a), nb = normP(b);
            if (na < 1e-15 || nb < 1e-15) { worst = 0.0; continue; }
            const double sj = cross2(a, b) / (na * nb);
            worst = std::min(worst, sj);
            const double ang = std::acos(std::max(-1.0, std::min(1.0, dotP(a, b) / (na * nb))));
            const double dev = std::fabs(ang - M_PI_2) * 180.0 / M_PI;
            q.maxAngleDeviation = std::max(q.maxAngleDeviation, dev);
            devSum += dev;
        }
        if (worst <= 0.0) ++q.invertedQuads;
        q.minScaledJacobian = std::min(q.minScaledJacobian, worst);
        sjSum += worst;

        const double d0 = normP(verts[c[2]] - verts[c[0]]);
        const double d1 = normP(verts[c[3]] - verts[c[1]]);
        const double diag = std::max(d0, d1);
        const double L = (diag > 1e-15) ? std::fabs(area) / diag : 0.0;
        zones.push_back(L / scale);
        zoneSum += L / scale;
    }

    q.meanScaledJacobian = sjSum / quads.size();
    q.meanAngleDeviation = devSum / (4.0 * quads.size());
    std::sort(zones.begin(), zones.end());
    q.minZoneDimension = zones.front();
    q.p01ZoneDimension = zones[static_cast<std::size_t>(0.01 * (zones.size() - 1))];
    q.meanZoneDimension = zoneSum / zones.size();
    return q;
}

} // namespace metrics
} // namespace paper
