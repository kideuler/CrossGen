#include "TORSION.hxx"

#include "TORSION/ConeMetric.hxx"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "dualmbo/DualMBO.hxx"

// The Jacobian of a piecewise-linear map on Omega, one 2x2 row-major per face,
// in the layout FieldIntegration reads a target in: row 0 is grad u and row 1
// is grad v, which is how J*_t is written too.
//
// This is what makes Sec. 6.5's re-projection a swap of one argument rather
// than a second solver: a map's own Jacobian is a legal target for the fit that
// produced it, and fitting it back under the constraints is a projection.
//
// A member rather than a file-local: the viewer drives Stage 4R itself and its
// re-projection has to fit the same object. See TORSION.hxx.
std::vector<std::array<double, 4>> TORSION::jacobianOf(const Mesh &om,
                                                       const std::vector<Point> &uv) {
    const int nT = static_cast<int>(om.triangles.size());
    std::vector<std::array<double, 4>> J(nT, {1.0, 0.0, 0.0, 1.0});
    if (uv.size() != om.vertices.size()) return J;
    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = om.triangles[t];
        const Point &p0 = om.vertices[tri[0]];
        const Point &p1 = om.vertices[tri[1]];
        const Point &p2 = om.vertices[tri[2]];
        const double twoA = cross2(p1 - p0, p2 - p0);
        if (!(std::fabs(twoA) > 0.0)) continue;
        Point gu{0.0, 0.0}, gv{0.0, 0.0};
        for (int i = 0; i < 3; ++i) {
            const Point e = om.vertices[tri[(i + 2) % 3]] - om.vertices[tri[(i + 1) % 3]];
            const Point g{-e[1] / twoA, e[0] / twoA};
            gu = gu + g * uv[tri[i]][0];
            gv = gv + g * uv[tri[i]][1];
        }
        J[t] = {gu[0], gu[1], gv[0], gv[1]};
    }
    return J;
}

// The signed angle a map on Omega puts at each vertex of S, summed over the
// vertex's children in Omega. A seam vertex's star is split between two
// children, and a repair that wraps one child's share a whole turn further
// round leaves every triangle positive and the vertex 2 pi out -- which is a
// cone changing valence as far as Q2 is concerned, and nothing downstream
// undoes it. Comparing these sums before and after a repair is the check that
// catches it; winding tests on the stars of Omega cannot, since a seam child's
// star is not a loop.
static std::vector<double> angleSumsOnS(const Mesh &om, const std::vector<int> &c2o,
                                        const std::vector<Point> &m, int nOrig) {
    std::vector<double> sum(static_cast<size_t>(std::max(nOrig, 0)), 0.0);
    if (m.size() != om.vertices.size() || c2o.size() != om.vertices.size()) return sum;
    for (const Triangle &tri : om.triangles) {
        for (int i = 0; i < 3; ++i) {
            const int w = tri[i];
            const int o = c2o[w];
            if (o < 0 || o >= nOrig) continue;
            const Point p = m[tri[(i + 1) % 3]] - m[w], q = m[tri[(i + 2) % 3]] - m[w];
            sum[o] += std::atan2(cross2(p, q), dotP(p, q));
        }
    }
    return sum;
}

// What one vertex of Omega may do while Sec. 7.2a untangles around it.
enum : unsigned char { kFixed = 0, kFreeInU = 1, kFreeInV = 2, kFree = 3 };

// What Sec. 6.4's alignment alone lets each vertex of Omega do, the seam left
// out of it: vertexFreedom() at rung 0 holds every seam child outright, and a
// paired move needs to know which of them the alignment would have held too.
static std::vector<unsigned char> alignmentFreedom(const ConeCut &cut,
                                                   const std::vector<int> &usedAxis) {
    const Mesh &om = cut.getCutMesh();
    std::vector<unsigned char> freedom(om.vertices.size(), kFree);
    if (usedAxis.size() != om.edges.size()) return freedom;
    for (size_t e = 0; e < om.edges.size(); ++e) {
        const int hold = usedAxis[e];
        if (hold != 0 && hold != 1) continue;
        const unsigned char keep = (hold == 0) ? kFreeInV : kFreeInU;
        for (int i = 0; i < 2; ++i) {
            unsigned char &f = freedom[om.edges[e][i]];
            f = static_cast<unsigned char>(f & keep);
        }
    }
    return freedom;
}

// ---------------------------------------------------------------------------
// relaxToKernel()  --  Sec. 7.2a, the local untangling
//
// Sec. 7.2's repair throws the map away and starts again from a Tutte embedding
// of a circle. That is the right answer when the integration came back badly
// tangled and it is a wildly disproportionate one when it came back with a
// single inverted face out of eighteen thousand, which is the usual case: the
// circle has none of psi_0's structure -- not the cone angles the seam
// constraints put there, not the boundary Sec. 6.4 held to an axis -- and the
// fit has to rebuild all of it from nothing, at the cost of a continuation and,
// on the harder models, of a cone that comes back 2 pi out.
//
// A handful of inverted faces does not need any of that. Hold every vertex of
// dOmega where it is -- which is every vertex the seam and the alignment have
// anything to say about -- and move the interior ones off the tangle. For one
// interior vertex the statement is exact and convex:
//
//     the signed area of each incident image triangle is an *affine* function
//     of that one vertex's position, so "every triangle at v is positively
//     oriented" is an intersection of half planes, and its interior is the
//     kernel of the one-ring polygon.
//
// Putting v at the Chebyshev centre of that kernel makes every triangle at v
// positive with the largest margin available, and touches nothing else. The
// margin of a triangle is its signed area over half the length of the edge
// opposite v, which is the signed distance from v to that edge, so the centre
// is the solution of
//
//     max r   subject to   u_t . x + e_t >= r   for every incident triangle t
//
// a three-variable linear program whose optimum has three active constraints.
// A one ring has six or so triangles, so the twenty triples are enumerated
// rather than pivoted.
//
// A sweep does this for every interior vertex within two rings of an inverted
// face, and the sweeps stop when nothing inverted is left or nothing moved.
// Where it cannot finish -- a tangle whose only free vertices are on dOmega, or
// one deep enough that no single vertex can be moved out of it -- Sec. 7.2's
// Tutte pass is still there behind it, and now runs on the models that actually
// need it.
//
// `freedom` says what each vertex of Omega is allowed to do, which is the whole
// of what keeps this pass honest about the equalities Stage 4F imposed: a seam
// vertex cannot move at all, because its partner would have to move with it; a
// vertex on a chain Sec. 6.4 held cannot move across the axis it was held to,
// but may slide along it; everything else is free.
//
// Returns the number of faces still inverted.
// ---------------------------------------------------------------------------
int TORSION::relaxToKernel(const Mesh &om, std::vector<Point> &uv,
                           const std::vector<unsigned char> &freedom,
                           const std::vector<int> &partner, const std::vector<int> &partnerK,
                           int sweeps) {
    const int nV = static_cast<int>(om.vertices.size());
    const int nT = static_cast<int>(om.triangles.size());
    if (uv.size() != static_cast<size_t>(nV)) return -1;
    if (freedom.size() != static_cast<size_t>(nV)) return -1;
    const bool paired = partner.size() == static_cast<size_t>(nV) &&
                        partnerK.size() == static_cast<size_t>(nV);

    auto area = [&](int t) {
        const Triangle &tri = om.triangles[t];
        return 0.5 * cross2(uv[tri[1]] - uv[tri[0]], uv[tri[2]] - uv[tri[0]]);
    };
    auto countInverted = [&]() {
        int n = 0;
        for (int t = 0; t < nT; ++t) if (!(area(t) > 0.0)) ++n;
        return n;
    };

    const auto &vt = om.vertexTriangles;
    std::vector<char> onBoundary(nV, 0);
    for (int e : om.boundaryEdges) {
        onBoundary[om.edges[e][0]] = 1;
        onBoundary[om.edges[e][1]] = 1;
    }

    // Every triangle at a vertex being positively oriented is an intersection
    // of half planes on that vertex's position, and each is written here as a
    // margin on the *displacement*: margin(d) = n . d + c, with n a unit normal
    // and c the margin the map has now. A seam vertex brings its partner's half
    // planes in as well, with their normals turned by R_k -- the two move
    // together under phi_+ = R_k phi_- + t, so one displacement moves both and
    // the programme is still two-dimensional.
    std::vector<Point> nrm;
    std::vector<double> cst;
    std::vector<Point> ringA, ringB;   // the v ring only, for the winding test
    size_t nTriV = 0;

    auto rotate = [](const Point &p, int k) { return Immersion::rotateQuarter(p, k); };

    auto push = [&](int w, int k, bool keepRing) {
        const int begin = vt.rowPtr[w], end = vt.rowPtr[w + 1];
        if (end - begin < 3) return false;
        for (int q = begin; q < end; ++q) {
            const Triangle &tri = om.triangles[vt.colIdx[q]];
            int a = -1, b = -1;
            for (int i = 0; i < 3; ++i) {
                if (tri[i] == w) { a = tri[(i + 1) % 3]; b = tri[(i + 2) % 3]; }
            }
            if (a < 0) return false;
            const Point n{0.5 * (uv[a][1] - uv[b][1]), 0.5 * (uv[b][0] - uv[a][0])};
            const double len = normP(n);
            if (!(len > 0.0)) return false;
            const Point un = n / len;
            const double margin = dotP(un, uv[w]) + 0.5 * cross2(uv[a], uv[b]) / len;
            // d moves w by R_k d, so the half plane reads (R_k^T n) . d + margin.
            nrm.push_back(rotate(un, (4 - ((k % 4) + 4) % 4) % 4));
            cst.push_back(margin);
            if (keepRing) { ringA.push_back(uv[a]); ringB.push_back(uv[b]); }
        }
        return true;
    };

    // The one ring of v, seen from where a displacement would put it, has to
    // wind exactly once. Every triangle being positively oriented does not say
    // that: a ring polygon that crosses itself -- which is what a tangle is --
    // has a region where it winds twice, every triangle at a point of it is
    // still positive, and the angle sum there is 4 pi. Q2 reports that as a
    // vertex 2 pi out, and it is the one way this pass can turn a tangle into
    // something worse than a tangle.
    auto windsOnce = [&](const Point &x) {
        double turn = 0.0;
        for (size_t i = 0; i < nTriV; ++i) {
            const Point p = ringA[i] - x, q = ringB[i] - x;
            turn += std::atan2(cross2(p, q), dotP(p, q));
        }
        return turn < 2.0 * M_PI + 0.05;
    };

    // The same test at a vertex other than the one being moved, read off the map
    // as it stands: moving v changes the one-ring polygon of each of its
    // neighbours as well as its own, and a move that unwraps v by wrapping a
    // neighbour is no better than the tangle it started from. An open fan has
    // nothing to wind, so a vertex of dOmega always passes.
    auto ringWrapsOnce = [&](int w) {
        if (onBoundary[w]) return true;   // an open fan; there is nothing to wind
        const int begin = vt.rowPtr[w], end = vt.rowPtr[w + 1];
        double sum = 0.0;
        for (int k = begin; k < end; ++k) {
            const Triangle &tri = om.triangles[vt.colIdx[k]];
            int a = -1, b = -1;
            for (int i = 0; i < 3; ++i) {
                if (tri[i] == w) { a = tri[(i + 1) % 3]; b = tri[(i + 2) % 3]; }
            }
            if (a < 0) return false;
            const Point p = uv[a] - uv[w], q = uv[b] - uv[w];
            sum += std::atan2(std::fabs(cross2(p, q)), dotP(p, q));
        }
        return sum < 2.0 * M_PI + 0.05;
    };
    auto neighbourhoodHolds = [&](int w) {
        if (!ringWrapsOnce(w)) return false;
        const int begin = vt.rowPtr[w], end = vt.rowPtr[w + 1];
        for (int k = begin; k < end; ++k) {
            const Triangle &tri = om.triangles[vt.colIdx[k]];
            for (int i = 0; i < 3; ++i) {
                if (tri[i] != w && !ringWrapsOnce(tri[i])) return false;
            }
        }
        return true;
    };

    if (countInverted() == 0) return 0;

    // The vertices allowed to move: those of an inverted face to start with,
    // and one more ring each time a sweep runs out of things to do. Keeping the
    // set as small as will work is the point -- every vertex this pass moves is
    // one whose position the integration chose and this one is overruling.
    std::vector<char> candidate(nV, 0);
    for (int t = 0; t < nT; ++t) {
        if (area(t) > 0.0) continue;
        for (int i = 0; i < 3; ++i) candidate[om.triangles[t][i]] = 1;
    }
    auto growOneRing = [&]() {
        std::vector<char> grown = candidate;
        for (int t = 0; t < nT; ++t) {
            const Triangle &tri = om.triangles[t];
            if (!candidate[tri[0]] && !candidate[tri[1]] && !candidate[tri[2]]) continue;
            for (int i = 0; i < 3; ++i) grown[tri[i]] = 1;
        }
        candidate.swap(grown);
    };

    std::vector<std::pair<double, Point>> tries;
    int rings = 0;
    for (int sweep = 0; sweep < sweeps; ++sweep) {
        bool moved = false;
        for (int v = 0; v < nV; ++v) {
            if (!candidate[v] || freedom[v] == kFixed) continue;
            int mate = paired ? partner[v] : -1;
            if (mate >= 0 && mate < v) continue;      // handled at the smaller index
            if (mate >= 0 && freedom[mate] == kFixed) continue;
            // The two half-plane sets are written as though each ring moved on
            // its own. A triangle carrying both children moves twice and the
            // model is then wrong about it, so such a pair is left alone -- it
            // happens at the tip of a one-edge slit and nowhere else.
            if (mate >= 0) {
                for (int q = vt.rowPtr[v]; q < vt.rowPtr[v + 1] && mate >= 0; ++q) {
                    const Triangle &tri = om.triangles[vt.colIdx[q]];
                    for (int i = 0; i < 3; ++i) if (tri[i] == mate) mate = -2;
                }
                if (mate == -2) continue;
            }

            const Point wasVpos = uv[v];
            nrm.clear(); cst.clear(); ringA.clear(); ringB.clear();
            if (!push(v, 0, true)) continue;
            nTriV = nrm.size();
            if (mate >= 0 && !push(mate, partnerK[v], false)) continue;
            if (nrm.size() < 3) continue;

            // A box of one mean ring distance, in the same units -- both are
            // distances. It bounds the programme, which a one ring does not
            // when the fan is open, and it bounds the move, which is what makes
            // this a local repair rather than a global one.
            double reach = 0.0;
            for (size_t i = 0; i < nTriV; ++i) reach += normP(ringA[i] - wasVpos);
            reach /= static_cast<double>(nTriV);
            if (!(reach > 0.0)) continue;
            const Point box[4] = {{1.0, 0.0}, {-1.0, 0.0}, {0.0, 1.0}, {0.0, -1.0}};
            for (int i = 0; i < 4; ++i) { nrm.push_back(box[i]); cst.push_back(reach); }

            // One coordinate held is Sec. 6.4's alignment, which this pass may
            // not break: the search is then along the line the vertex may still
            // slide on.
            const bool line = (freedom[v] != kFree);
            const Point dir = (freedom[v] == kFreeInU) ? Point{1.0, 0.0} : Point{0.0, 1.0};

            auto margin = [&](const Point &d) {
                double r = std::numeric_limits<double>::infinity();
                for (size_t i = 0; i < nrm.size(); ++i) r = std::min(r, dotP(nrm[i], d) + cst[i]);
                return r;
            };
            const double current = margin(Point{0.0, 0.0});

            // max r  s.t.  n_i . d + c_i >= r, bounded, so its optimum has
            // three active constraints -- two on a line. A one ring is six or
            // so triangles, so they are enumerated rather than pivoted.
            tries.clear();
            const int m = static_cast<int>(nrm.size());
            if (line) {
                for (int i = 0; i < m; ++i) {
                    for (int j = i + 1; j < m; ++j) {
                        const double ds = dotP(nrm[i], dir) - dotP(nrm[j], dir);
                        if (std::fabs(ds) < 1e-14) continue;
                        const Point d = dir * ((cst[j] - cst[i]) / ds);
                        const double r = margin(d);
                        if (r > current) tries.emplace_back(r, d);
                    }
                }
            } else {
                for (int i = 0; i < m; ++i) {
                    for (int j = i + 1; j < m; ++j) {
                        for (int k = j + 1; k < m; ++k) {
                            const double a11 = nrm[i][0] - nrm[k][0], a12 = nrm[i][1] - nrm[k][1];
                            const double a21 = nrm[j][0] - nrm[k][0], a22 = nrm[j][1] - nrm[k][1];
                            const double det = a11 * a22 - a12 * a21;
                            if (std::fabs(det) < 1e-14) continue;
                            const double b1 = cst[k] - cst[i], b2 = cst[k] - cst[j];
                            const Point d{(b1 * a22 - a12 * b2) / det,
                                          (a11 * b2 - b1 * a21) / det};
                            const double r = margin(d);
                            if (r > current) tries.emplace_back(r, d);
                        }
                    }
                }
            }
            if (tries.empty()) continue;
            std::sort(tries.begin(), tries.end(),
                      [](const std::pair<double, Point> &a, const std::pair<double, Point> &b) {
                          return a.first > b.first;
                      });

            // Best first, and take the best that leaves every ring it touched
            // winding once. Trying rather than filtering is what lets a vertex
            // settle for a smaller margin when the largest one would have
            // wrapped a neighbour.
            // The best displacement whose own ring still winds once, applied,
            // and then kept only if every ring it touched still does -- moving
            // v changes the one-ring polygon of each of its neighbours as well
            // as its own, and a move that unwraps v by wrapping a neighbour is
            // no better than the tangle it started from.
            const Point wasV = wasVpos;
            const Point wasM = (mate >= 0) ? uv[mate] : Point{0.0, 0.0};
            auto place = [&](const Point &d) {
                uv[v] = wasV + d;
                if (mate >= 0) uv[mate] = wasM + rotate(d, partnerK[v]);
            };
            auto restore = [&]() {
                uv[v] = wasV;
                if (mate >= 0) uv[mate] = wasM;
            };

            double bestR = current;
            Point bestD{0.0, 0.0};
            bool have = false;
            for (const auto &t : tries) {
                if (!(t.first > bestR)) continue;
                if (!windsOnce(wasV + t.second)) continue;
                bestR = t.first; bestD = t.second; have = true;
            }
            if (!have) continue;
            place(bestD);
            if (neighbourhoodHolds(v) && (mate < 0 || neighbourhoodHolds(mate))) moved = true;
            else restore();
        }
        if (countInverted() == 0) return 0;
        if (!moved) {
            if (rings >= 3) break;
            growOneRing();
            ++rings;
        }
    }
    return countInverted();
}

// ---------------------------------------------------------------------------
// untangleRegularised()  --  Sec. 7.2c, the fold relaxToKernel() cannot open
//
// relaxToKernel() moves one vertex at a time and only ever to a position where
// its own star is better than it was, which clears a tangle one vertex deep and
// cannot clear a fold: a run of triangles turned over together, where every
// vertex that would have to move is pinned by a neighbour that is itself
// inverted. The fold that matters here is at a cone. A -1 cone opens its fan
// to 5 pi / 2 and a +1 closes it to 3 pi / 2; the frame asks every triangle of
// the fan for the angle it has on the model, the least-squares fit takes up
// the quarter turn it is owed in whichever triangles it can, and on tunnel that
// is fifteen triangles turned over round one -1 cone, of which the kernel pass
// clears two.
//
// What opens a fold is an energy that is finite on an inverted triangle, so
// that the optimisation can pass through one, and that becomes a barrier as it
// converges, so that it cannot stop there. That is Garanzha et al. (2021), and
// the same regularised energy mesh::TMOP's UntangleRegular runs on quads:
//
//     f(J) = [ (1 - g) |J|^2 / 2 + g (det J^2 + 1) / 2 ] / chi(det J, eps)
//     chi(d, eps) = (d + sqrt(eps^2 + d^2)) / 2
//
// with J the map's Jacobian against the frame's own scale, so that |J|^2 / 2
// det J is the conformal distortion and (det^2 + 1) / 2 det the size, g = 1/128
// between them, and eps shrinking with the worst det. Minimised over the
// vertices within a few rings of the tangle, node by node, by Newton with a
// backtracking line search on the node's star; the rest of the map does not
// move, and `freedom` holds Sec. 7.2a's rungs exactly as relaxToKernel() does.
//
// Every triangle positive is not quite enough, for the reason relaxToKernel()
// checks winding: a ring can wind twice with every triangle in it positive. The
// result is kept only if no interior ring of Omega in the region does.
//
// Returns the number of faces still inverted; the map is left as it was if the
// pass did not reduce it.
// ---------------------------------------------------------------------------
int TORSION::untangleRegularised(const Mesh &om, std::vector<Point> &uv,
                                 const std::vector<unsigned char> &freedom,
                                 const std::vector<int> &partner,
                                 const std::vector<int> &partnerK,
                                 const std::vector<std::array<double, 4>> &target,
                                 const std::vector<int> &cutToOriginal,
                                 const std::vector<double> &prescribedAngle, int rings,
                                 int maxIterations) {
    const int nV = static_cast<int>(om.vertices.size());
    const int nT = static_cast<int>(om.triangles.size());
    if (uv.size() != static_cast<size_t>(nV) || freedom.size() != static_cast<size_t>(nV)) return -1;

    auto signedArea = [&](const std::vector<Point> &m, int t) {
        const Triangle &tri = om.triangles[t];
        return 0.5 * cross2(m[tri[1]] - m[tri[0]], m[tri[2]] - m[tri[0]]);
    };
    auto countInverted = [&](const std::vector<Point> &m) {
        int n = 0;
        for (int t = 0; t < nT; ++t) if (!(signedArea(m, t) > 0.0)) ++n;
        return n;
    };
    const int before = countInverted(uv);
    if (before == 0) return 0;

    // Per face: G = M^-1 / s, with M the model's edge matrix and s the frame's
    // scale, so that J = D G for D the image's edge matrix; and the model area.
    std::vector<std::array<double, 4>> G(nT, {0.0, 0.0, 0.0, 0.0});
    std::vector<double> weight(nT, 0.0);
    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = om.triangles[t];
        const Point a = om.vertices[tri[1]] - om.vertices[tri[0]];
        const Point b = om.vertices[tri[2]] - om.vertices[tri[0]];
        const double det = a[0] * b[1] - a[1] * b[0];
        if (!(std::fabs(det) > 0.0)) continue;
        double s = 1.0;
        if (target.size() == static_cast<size_t>(nT)) {
            const std::array<double, 4> &J = target[t];
            const double d = J[0] * J[3] - J[1] * J[2];
            if (d > 0.0 && std::isfinite(d)) s = std::sqrt(d);
        }
        // M = [a b] (columns); M^-1 = [[b1, -b0], [-a1, a0]] / det. Rows of G:
        // g0 = row 0 of M^-1 / s, g1 = row 1.
        G[t] = {b[1] / (det * s), -b[0] / (det * s), -a[1] / (det * s), a[0] / (det * s)};
        weight[t] = 0.5 * std::fabs(det);
    }

    const double g = 1.0 / 128.0;
    const auto &vt = om.vertexTriangles;
    const bool paired = partner.size() == static_cast<size_t>(nV) &&
                        partnerK.size() == static_cast<size_t>(nV);
    auto sharesTriangle = [&](int a, int b) {
        for (int k = vt.rowPtr[a]; k < vt.rowPtr[a + 1]; ++k) {
            const Triangle &tri = om.triangles[vt.colIdx[k]];
            if (tri[0] == b || tri[1] == b || tri[2] == b) return true;
        }
        return false;
    };
    std::vector<char> onBoundary(nV, 0);
    for (int e : om.boundaryEdges) {
        onBoundary[om.edges[e][0]] = 1;
        onBoundary[om.edges[e][1]] = 1;
    }

    // J on face t with the current map.
    auto jacobian = [&](const std::vector<Point> &m, int t, double J[4]) {
        const Triangle &tri = om.triangles[t];
        const Point d1 = m[tri[1]] - m[tri[0]], d2 = m[tri[2]] - m[tri[0]];
        const std::array<double, 4> &Gt = G[t];
        // J = d1 g0^T + d2 g1^T
        J[0] = d1[0] * Gt[0] + d2[0] * Gt[2];
        J[1] = d1[0] * Gt[1] + d2[0] * Gt[3];
        J[2] = d1[1] * Gt[0] + d2[1] * Gt[2];
        J[3] = d1[1] * Gt[1] + d2[1] * Gt[3];
    };
    auto chi = [](double d, double eps) { return 0.5 * (d + std::sqrt(eps * eps + d * d)); };
    auto faceEnergy = [&](const std::vector<Point> &m, int t, double eps) {
        double J[4];
        jacobian(m, t, J);
        const double f2 = J[0] * J[0] + J[1] * J[1] + J[2] * J[2] + J[3] * J[3];
        const double d = J[0] * J[3] - J[1] * J[2];
        return weight[t] * ((1.0 - g) * 0.5 * f2 + g * 0.5 * (d * d + 1.0)) / chi(d, eps);
    };
    auto starEnergy = [&](const std::vector<Point> &m, int v, double eps) {
        double e = 0.0;
        for (int k = vt.rowPtr[v]; k < vt.rowPtr[v + 1]; ++k) e += faceEnergy(m, vt.colIdx[k], eps);
        return e;
    };
    auto minDet = [&](const std::vector<Point> &m, const std::vector<int> &faces) {
        double worst = std::numeric_limits<double>::infinity();
        for (int t : faces) {
            double J[4];
            jacobian(m, t, J);
            worst = std::min(worst, J[0] * J[3] - J[1] * J[2]);
        }
        return worst;
    };
    // The signed angle the map puts at Omega vertex w, over its own star. A
    // vertex of S with several children in Omega -- a seam vertex -- has its
    // angle sum split between them, and a fold can be opened by wrapping one
    // child's share a whole turn further round with every triangle positive.
    // So the check is on the sum over the children, per vertex of S, against
    // the sum Q2 prescribes there -- not against the map the pass started
    // from, whose sums are a turn out wherever a triangle was turned over.
    // Without the prescription, an interior vertex of Omega is held to one
    // turn and nothing else is checked.
    auto starAngle = [&](const std::vector<Point> &m, int w) {
        double sum = 0.0;
        for (int k = vt.rowPtr[w]; k < vt.rowPtr[w + 1]; ++k) {
            const Triangle &tri = om.triangles[vt.colIdx[k]];
            int a = -1, b = -1;
            for (int i = 0; i < 3; ++i) if (tri[i] == w) { a = tri[(i + 1) % 3]; b = tri[(i + 2) % 3]; }
            if (a < 0) continue;
            const Point p = m[a] - m[w], q = m[b] - m[w];
            sum += std::atan2(cross2(p, q), dotP(p, q));
        }
        return sum;
    };
    const bool haveParents = cutToOriginal.size() == static_cast<size_t>(nV) &&
                             !prescribedAngle.empty();
    std::unordered_map<int, std::vector<int>> childrenOf;
    if (haveParents) {
        for (int w = 0; w < nV; ++w) childrenOf[cutToOriginal[w]].push_back(w);
    }
    auto sumsHold = [&](const std::vector<Point> &after, const std::vector<int> &touched) {
        std::vector<int> parents;
        for (int w : touched) {
            if (haveParents) parents.push_back(cutToOriginal[w]);
            for (int k = vt.rowPtr[w]; k < vt.rowPtr[w + 1]; ++k) {
                const Triangle &tri = om.triangles[vt.colIdx[k]];
                for (int i = 0; i < 3; ++i) {
                    if (haveParents) parents.push_back(cutToOriginal[tri[i]]);
                    else if (!onBoundary[tri[i]] &&
                             std::fabs(starAngle(after, tri[i]) - 2.0 * M_PI) > M_PI)
                        return false;
                }
            }
        }
        if (!haveParents) return true;
        std::sort(parents.begin(), parents.end());
        parents.erase(std::unique(parents.begin(), parents.end()), parents.end());
        for (int o : parents) {
            auto it = childrenOf.find(o);
            if (it == childrenOf.end() || o < 0 || o >= static_cast<int>(prescribedAngle.size())) continue;
            double now = 0.0;
            for (int w : it->second) now += starAngle(after, w);
            if (std::fabs(now - prescribedAngle[o]) > M_PI) {
                return false;
            }
        }
        return true;
    };

    std::vector<char> inRegion(nV, 0);
    for (int t = 0; t < nT; ++t) {
        if (signedArea(uv, t) > 0.0) continue;
        for (int i = 0; i < 3; ++i) inRegion[om.triangles[t][i]] = 1;
    }
    auto grow = [&]() {
        std::vector<char> grown = inRegion;
        for (const Triangle &tri : om.triangles) {
            if (!inRegion[tri[0]] && !inRegion[tri[1]] && !inRegion[tri[2]]) continue;
            for (int i = 0; i < 3; ++i) grown[tri[i]] = 1;
        }
        inRegion.swap(grown);
    };
    for (int r = 0; r < rings; ++r) grow();

    std::vector<Point> best = uv;
    int bestLeft = before;
    // Three widening neighbourhoods, and then the whole map: a tangle of
    // hundreds of faces is not local to anything, and started from the
    // integration's own map the energy keeps every star's winding, so a pass
    // over everything is still a repair of that map and not a new one -- which
    // is the difference from Sec. 7.2's Tutte pass, whose circle flattens every
    // cone at the tip of a slit to a half turn and lets the fit settle it a
    // whole turn from its index.
    for (int attempt = 0; attempt < 4 && bestLeft > 0; ++attempt) {
        if (attempt > 0 && attempt < 3) for (int r = 0; r < rings; ++r) grow();
        if (attempt == 3) std::fill(inRegion.begin(), inRegion.end(), 1);
        std::vector<int> movers, faces;
        std::vector<char> faceIn(nT, 0);
        for (int v = 0; v < nV; ++v) {
            if (!inRegion[v] || freedom[v] == kFixed) continue;
            const int mate = paired ? partner[v] : -1;
            // A seam child moves with its partner, handled at the smaller index
            // of the two; a pair whose stars share a triangle -- the tip of a
            // one-edge slit -- is left where it is.
            if (mate >= 0 && (mate < v || freedom[mate] != kFree || freedom[v] != kFree ||
                              sharesTriangle(v, mate))) continue;
            movers.push_back(v);
            for (int k = vt.rowPtr[v]; k < vt.rowPtr[v + 1]; ++k) faceIn[vt.colIdx[k]] = 1;
            if (mate >= 0) {
                movers.push_back(-1 - mate);   // its partner, marked, for the checks
                for (int k = vt.rowPtr[mate]; k < vt.rowPtr[mate + 1]; ++k) faceIn[vt.colIdx[k]] = 1;
            }
        }
        for (int t = 0; t < nT; ++t) if (faceIn[t]) faces.push_back(t);
        if (movers.empty()) continue;

        std::vector<Point> m = best;
        // A region whose inverted count has not come down for a while is not
        // going to clear at this size; the next, wider one is tried instead.
        int leastInverted = std::numeric_limits<int>::max(), sinceProgress = 0;
        for (int it = 0; it < maxIterations; ++it) {
            const double worst = minDet(m, faces);
            const int invertedNow = countInverted(m);
            if (worst > 0.0 && invertedNow == 0) break;
            if (invertedNow < leastInverted) { leastInverted = invertedNow; sinceProgress = 0; }
            else if (++sinceProgress >= 15) break;
            const double eps = std::sqrt(1e-10 + 0.04 * std::min(worst, 0.0) * std::min(worst, 0.0)) +
                               (worst > 0.0 ? 1e-6 : 0.0);
            // Gradient and Hessian of w's star energy in w's own position: on
            // each face J = x w^T + C, affine in x.
            auto starDerivatives = [&](int v, double eps, double &gx, double &gy, double &hxx,
                                       double &hxy, double &hyy) {
                gx = gy = hxx = hxy = hyy = 0.0;
                for (int k = vt.rowPtr[v]; k < vt.rowPtr[v + 1]; ++k) {
                    const int t = vt.colIdx[k];
                    if (!(weight[t] > 0.0)) continue;
                    const Triangle &tri = om.triangles[t];
                    const std::array<double, 4> &Gt = G[t];
                    double w[2];
                    if (tri[1] == v) { w[0] = Gt[0]; w[1] = Gt[1]; }
                    else if (tri[2] == v) { w[0] = Gt[2]; w[1] = Gt[3]; }
                    else { w[0] = -(Gt[0] + Gt[2]); w[1] = -(Gt[1] + Gt[3]); }
                    double J[4];
                    jacobian(m, t, J);
                    const double x0 = m[v][0], x1 = m[v][1];
                    // C = J - x w^T
                    const double c00 = J[0] - x0 * w[0], c01 = J[1] - x0 * w[1];
                    const double c10 = J[2] - x1 * w[0], c11 = J[3] - x1 * w[1];
                    const double ww = w[0] * w[0] + w[1] * w[1];
                    // q = adj(C)^T w, so that det J = det C + q . x
                    const double q0 = c11 * w[0] - c10 * w[1];
                    const double q1 = -c01 * w[0] + c00 * w[1];
                    const double f2 = J[0] * J[0] + J[1] * J[1] + J[2] * J[2] + J[3] * J[3];
                    const double d = J[0] * J[3] - J[1] * J[2];
                    const double root = std::sqrt(eps * eps + d * d);
                    const double c = 0.5 * (d + root);
                    const double c1 = 0.5 * (1.0 + d / root);
                    const double c2 = 0.5 * eps * eps / (root * root * root);
                    const double N = (1.0 - g) * 0.5 * f2 + g * 0.5 * (d * d + 1.0);
                    // grad N = (1-g)(|w|^2 x + C w) + g d q
                    const double Cw0 = c00 * w[0] + c01 * w[1], Cw1 = c10 * w[0] + c11 * w[1];
                    const double n0 = (1.0 - g) * (ww * x0 + Cw0) + g * d * q0;
                    const double n1 = (1.0 - g) * (ww * x1 + Cw1) + g * d * q1;
                    const double A = weight[t];
                    gx += A * (n0 / c - N * c1 * q0 / (c * c));
                    gy += A * (n1 / c - N * c1 * q1 / (c * c));
                    const double k2 = N * (2.0 * c1 * c1 / (c * c * c) - c2 / (c * c));
                    const double hN00 = (1.0 - g) * ww + g * q0 * q0;
                    const double hN01 = g * q0 * q1;
                    const double hN11 = (1.0 - g) * ww + g * q1 * q1;
                    hxx += A * (hN00 / c - 2.0 * c1 * n0 * q0 / (c * c) + k2 * q0 * q0);
                    hxy += A * (hN01 / c - c1 * (n0 * q1 + n1 * q0) / (c * c) + k2 * q0 * q1);
                    hyy += A * (hN11 / c - 2.0 * c1 * n1 * q1 / (c * c) + k2 * q1 * q1);
                }
            };
            for (int sweep = 0; sweep < 4; ++sweep) {
                for (int v : movers) {
                    if (v < 0) continue;   // a partner, moved with its pair
                    const int mate = paired ? partner[v] : -1;
                    double gx, gy, hxx, hxy, hyy;
                    starDerivatives(v, eps, gx, gy, hxx, hxy, hyy);
                    Point e1{1.0, 0.0}, e2{0.0, 1.0};
                    if (mate >= 0) {
                        // phi_mate moves by R_k d, so its star enters through
                        // R^T grad and R^T H R.
                        double mgx, mgy, mhxx, mhxy, mhyy;
                        starDerivatives(mate, eps, mgx, mgy, mhxx, mhxy, mhyy);
                        e1 = Immersion::rotateQuarter(Point{1.0, 0.0}, partnerK[v]);
                        e2 = Immersion::rotateQuarter(Point{0.0, 1.0}, partnerK[v]);
                        auto Hq = [&](const Point &a, const Point &b) {
                            return a[0] * (mhxx * b[0] + mhxy * b[1]) +
                                   a[1] * (mhxy * b[0] + mhyy * b[1]);
                        };
                        gx += e1[0] * mgx + e1[1] * mgy;
                        gy += e2[0] * mgx + e2[1] * mgy;
                        hxx += Hq(e1, e1);
                        hxy += Hq(e1, e2);
                        hyy += Hq(e2, e2);
                    }
                    Point step{0.0, 0.0};
                    if (freedom[v] == kFreeInU) {
                        if (hxx > 1e-300) step = Point{-gx / hxx, 0.0};
                        else step = Point{-gx, 0.0};
                    } else if (freedom[v] == kFreeInV) {
                        if (hyy > 1e-300) step = Point{0.0, -gy / hyy};
                        else step = Point{0.0, -gy};
                    } else {
                        // Make the 2x2 positive definite before solving.
                        const double tr = hxx + hyy;
                        const double disc = std::sqrt(std::max(0.0, 0.25 * (hxx - hyy) * (hxx - hyy) + hxy * hxy));
                        const double lmin = 0.5 * tr - disc;
                        const double shift = (lmin < 1e-8 * std::max(1.0, std::fabs(tr)))
                                                 ? (1e-8 * std::max(1.0, std::fabs(tr)) - lmin)
                                                 : 0.0;
                        const double a = hxx + shift, b = hxy, dd = hyy + shift;
                        const double det = a * dd - b * b;
                        if (!(det > 0.0)) continue;
                        step = Point{-(dd * gx - b * gy) / det, -(a * gy - b * gx) / det};
                    }
                    if (!std::isfinite(step[0]) || !std::isfinite(step[1])) continue;
                    const Point was = m[v];
                    const Point wasM = (mate >= 0) ? m[mate] : Point{0.0, 0.0};
                    auto energyNow = [&]() {
                        return starEnergy(m, v, eps) + (mate >= 0 ? starEnergy(m, mate, eps) : 0.0);
                    };
                    // The stars a move can re-wind: the moved vertices' own and
                    // every neighbour's. A step is refused if any of them turns
                    // by a whole turn -- which happens when a triangle passes
                    // through the degenerate position with its angle there near
                    // a half turn, and is how an energy that may cross inverted
                    // states reopens a -1 cone's fan the short way, a quarter
                    // turn instead of five. With every step winding-preserving,
                    // the pass keeps the angle sums the integration's seam put
                    // at every vertex it started from.
                    std::vector<int> watch;
                    for (int c : {v, mate}) {
                        if (c < 0) continue;
                        for (int k = vt.rowPtr[c]; k < vt.rowPtr[c + 1]; ++k) {
                            const Triangle &tri = om.triangles[vt.colIdx[k]];
                            for (int i = 0; i < 3; ++i) watch.push_back(tri[i]);
                        }
                    }
                    std::sort(watch.begin(), watch.end());
                    watch.erase(std::unique(watch.begin(), watch.end()), watch.end());
                    // With Q2's prescription to hand the test is directional: a
                    // vertex of S may be wound a turn *towards* its prescribed
                    // sum -- a cone the integration left wound the short way,
                    // opened the right way once its tip is free to move -- and
                    // never away from it. Without one, a turn either way is
                    // refused.
                    std::vector<int> parents;
                    if (haveParents) {
                        for (int w : watch) parents.push_back(cutToOriginal[w]);
                        std::sort(parents.begin(), parents.end());
                        parents.erase(std::unique(parents.begin(), parents.end()), parents.end());
                    }
                    auto gapOf = [&](int o) {
                        auto it = childrenOf.find(o);
                        if (it == childrenOf.end() || o < 0 ||
                            o >= static_cast<int>(prescribedAngle.size())) return 0.0;
                        double sum = 0.0;
                        for (int w : it->second) sum += starAngle(m, w);
                        return std::fabs(sum - prescribedAngle[o]);
                    };
                    std::vector<double> gap0(parents.size());
                    for (size_t i = 0; i < parents.size(); ++i) gap0[i] = gapOf(parents[i]);
                    std::vector<double> wound0(watch.size());
                    for (size_t i = 0; i < watch.size(); ++i) wound0[i] = starAngle(m, watch[i]);
                    auto windingKept = [&]() {
                        if (haveParents) {
                            for (size_t i = 0; i < parents.size(); ++i) {
                                if (gapOf(parents[i]) > gap0[i] + 0.5) return false;
                            }
                            return true;
                        }
                        for (size_t i = 0; i < watch.size(); ++i) {
                            if (std::fabs(starAngle(m, watch[i]) - wound0[i]) > M_PI) return false;
                        }
                        return true;
                    };
                    const double e0 = energyNow();
                    double alpha = 1.0;
                    bool accepted = false;
                    for (int bt = 0; bt < 30; ++bt, alpha *= 0.5) {
                        const Point d = step * alpha;
                        m[v] = was + d;
                        if (mate >= 0) m[mate] = wasM + Immersion::rotateQuarter(d, partnerK[v]);
                        if (energyNow() < e0 && windingKept()) { accepted = true; break; }
                    }
                    if (!accepted) {
                        m[v] = was;
                        if (mate >= 0) m[mate] = wasM;
                    }
                }
            }
        }
        const int left = countInverted(m);
        std::vector<int> touched;
        touched.reserve(movers.size());
        for (int v : movers) touched.push_back(v < 0 ? -1 - v : v);
        const bool wound = sumsHold(m, touched);
        if (wound && left < bestLeft) {
            bestLeft = left;
            best.swap(m);
        }
    }
    if (bestLeft < before) uv.swap(best);
    return bestLeft;
}

TORSION::TORSION(std::shared_ptr<Mesh> m) : TORSION(std::move(m), Options()) {}

TORSION::TORSION(std::shared_ptr<Mesh> m, const Options &opts)
    : mesh(std::move(m)), options(opts) {
    if (!mesh) throw std::runtime_error("TORSION: null mesh");
    if (mesh->triangles.empty()) throw std::runtime_error("TORSION: empty mesh");
}

TORSION::~TORSION() = default;

// ---------------------------------------------------------------------------
// runField()  --  Stage 0, MERIDIAN's unchanged
// ---------------------------------------------------------------------------
void TORSION::runField() {
    // One level of the tau-continuation: a fresh operator at this tau, with the
    // same Dirichlet data, started from the field the previous level left.
    auto buildLevel = [&](double tauScale, const Eigen::VectorXcd &carried) {
        auto f = std::make_unique<DualMBO>(mesh, options.dualMBOMaxSteps, options.dualMBOGamma);
        f->setPenaltyWeight(options.dualMBOWeight);
        f->setTauScale(tauScale);
        // On a multi-material domain the interfaces are Dirichlet data for the
        // field in exactly the way dS is, and that matters more here than in
        // Pipeline A rather than less: this pipeline integrates the field, so an
        // interface the field ran straight through is an interface the *map*
        // runs straight through. See DualMBO::setAlignedInteriorEdges for why
        // the hard pin and not a penalty.
        if (interfaces && interfaces->multiMaterial() && options.alignFieldToInterfaces) {
            f->setAlignedInteriorEdges(interfaces->interfaceEdges());
            status.fieldAlignedToInterfaces = true;
        }
        f->initialize();
        // The constraint does not travel with the field: every level re-imposes
        // the Dirichlet data its own assembly computed, and only the free
        // triangles are carried.
        if (carried.size() == f->u_k_prev.size()) {
            Eigen::VectorXcd start = carried;
            for (const auto &[ti, bc] : f->getBoundaryData()) start[ti] = bc;
            f->u_k_prev = start;
            f->u_k = start;
        }
        return f;
    };

    field = buildLevel(1.0, Eigen::VectorXcd());

    // Options::externalField: the caller has a field already and wants Stages
    // 0b to 11 run on it unchanged. initialize() still runs, because the
    // Dirichlet data it computes is read downstream as data even when the
    // solve it factorised is never used.
    if (options.externalField.size() > 0) {
        if (options.externalField.size() == static_cast<Eigen::Index>(mesh->triangles.size())) {
            field->u_k_prev = options.externalField;
            field->u_k = options.externalField;
            field->error = 0.0;
            status.fieldConverged = true;
            status.externalFieldUsed = true;
            status.messages.push_back(
                "Stage 0: using the caller's cross field; the MBO solve was skipped.");
            return;
        }
        std::ostringstream oss;
        oss << "the supplied cross field has " << options.externalField.size()
            << " value(s) and this mesh has " << mesh->triangles.size()
            << " triangle(s), so it was ignored and the MBO solve was run instead.";
        status.messages.push_back("Stage 0: " + oss.str());
    }

    const double nTris = static_cast<double>(mesh->triangles.size());

    // The ladder of tau scales, largest first. The first entry is always 1 --
    // the heuristic's tau -- so the continuation starts from the field the
    // single-tau scheme would have returned and only ever refines it. Its floor
    // is read off the assembled operator, so it cannot be known until level 0
    // has been built.
    std::vector<double> ladder{1.0};
    if (options.dualMBOTauContinuation && options.dualMBOTauRatio > 0.0 &&
        options.dualMBOTauRatio < 1.0) {
        double mnx = 1e300, mxx = -1e300, mny = 1e300, mxy = -1e300;
        for (const Point &p : mesh->vertices) {
            mnx = std::min(mnx, p[0]); mxx = std::max(mxx, p[0]);
            mny = std::min(mny, p[1]); mxy = std::max(mxy, p[1]);
        }
        const double D = std::hypot(mxx - mnx, mxy - mny);
        const double rate = field->medianDiffusionRate();
        if (D > 0.0 && rate > 0.0) {
            const double tau0 = D * D / 10.0;
            const double c = options.dualMBOTauFloorEdges;
            const double tauMin = c * c / rate;
            double tau = tau0;
            while (tau * options.dualMBOTauRatio > tauMin && ladder.size() < 40) {
                tau *= options.dualMBOTauRatio;
                ladder.push_back(ladder.back() * options.dualMBOTauRatio);
            }
        }
    }

    Eigen::VectorXcd carried;
    for (std::size_t level = 0; level < ladder.size(); ++level) {
        if (level > 0) field = buildLevel(ladder[level], carried);
        const int cap = (ladder.size() == 1)
                            ? options.dualMBOMaxSteps
                            : std::min(options.dualMBOMaxSteps, options.dualMBOTauLevelSteps);
        status.fieldConverged = false;
        for (int i = 0; i < cap; ++i) {
            field->step();
            ++status.mboSteps;
            if (field->error < 2.0 * nTris * 1e-5) { status.fieldConverged = true; break; }
        }
        carried = field->u_k_prev;
    }
    if (ladder.size() > 1) {
        std::ostringstream oss;
        oss << "the tau-continuation ran " << ladder.size() << " level(s) at ratio "
            << options.dualMBOTauRatio << " down to tau/tau_0 = " << ladder.back()
            << ", " << status.mboSteps << " MBO step(s) in total.";
        status.messages.push_back("Stage 0: " + oss.str());
    }
    if (!status.fieldConverged) {
        std::ostringstream oss;
        oss << "Cross field did not converge in " << options.dualMBOMaxSteps
            << " MBO steps (error " << field->error
            << "); every stage of this pipeline is downstream of it.";
        status.messages.push_back("Stage 0: " + oss.str());
    }
}

// ---------------------------------------------------------------------------
// runFieldFront() and runConeFront()  --  Stages 0b, 0, 1 and 2
//
// Identical to MERIDIAN's, deliberately and line for line: the interface
// network, the field, the cone indices with their prescriptions and dipole
// cancellations, the Gauss-Bonnet gate, the per-region balance, and the cut.
// Nothing in any of it depends on how psi_0 will be built, and the plan says so
// -- Stages 0b, 1 and 2 stay as they are -- so the only thing worth doing here
// is to keep it recognisably the same code.
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// exciseDisks()  --  Stage 0c, MERIDIAN's unchanged
//
// Before everything, because everything downstream is indexed on the mesh this
// leaves behind. See DiskTemplate.
// ---------------------------------------------------------------------------
void TORSION::exciseDisks() {
    std::vector<std::string> msgs;
    inclusions = DiskTemplate::detect(*mesh, msgs);
    for (const std::string &m : msgs) status.messages.push_back("Stage 0c: " + m);
    status.diskInclusions = static_cast<int>(inclusions.size());
    if (inclusions.empty()) return;

    std::shared_ptr<Mesh> excised = DiskTemplate::excise(*mesh, inclusions);
    if (!excised || excised->triangles.empty()) {
        inclusions.clear();
        status.diskInclusions = 0;
        status.messages.push_back(
            "Stage 0c: excising the circular inclusions would leave no mesh behind, so "
            "they are laid out like any other region.");
        return;
    }

    for (const DiskTemplate::Inclusion &inc : inclusions) {
        status.diskTrianglesExcised += static_cast<int>(inc.triangles.size());
    }
    inputMesh = mesh;
    mesh = excised;

    std::ostringstream oss;
    oss << "Excised " << inclusions.size() << " circular inclusion(s), "
        << status.diskTrianglesExcised << " of " << inputMesh->triangles.size()
        << " triangle(s): the layout is asked for the matrix with holes and Stage 11 "
        << "puts each O-grid back from a template.";
    status.messages.push_back("Stage 0c: " + oss.str());
}

void TORSION::runFieldFront() {
    if (options.diskTemplates) exciseDisks();

    if (options.materialInterfaces) {
        Interfaces::Options iopts;
        iopts.kinkAngle = options.interfaceKinkAngle;
        iopts.loopSplits = options.interfaceLoopSplits;
        interfaces = std::make_unique<Interfaces>(mesh, iopts);
        const Interfaces::Report &fr = interfaces->getReport();
        status.materials = fr.materials;
        status.interfaceEdges = fr.interfaceEdges;
        status.interfaceBranches = fr.branches;
        status.interfaceNodes = fr.nodes;
        status.interfaceIllPosedNodes = fr.illPosedNodes;
        status.interfaceWorstSector = fr.worstSectorResidual;
        for (const std::string &m : fr.messages) status.messages.push_back("Stage 0b: " + m);
    }

    runField();
}

// ---------------------------------------------------------------------------
// runConeFront()  --  Stages 1 and 2, on the field runFieldFront() solved
//
// Separate from Stage 0 because TORSION may run it twice on one field: once as
// MERIDIAN does, cancelling the same-region dipoles, and -- when that cost the
// alignment -- once keeping them. See run().
// ---------------------------------------------------------------------------
bool TORSION::runConeFront(bool cancelDipoles) {
    const size_t balanceMessagesSeen = interfaces ? interfaces->getReport().messages.size() : 0;

    cones = std::make_unique<ConeSingularities>(*field);
    cones->setBoundaryIndexRange(options.minBoundaryIndex, options.maxBoundaryIndex);
    // The snapshot Sec. 5.1's audit is run against, taken here because this is
    // the last moment at which the indices are the ones the *field* read.
    // prescribe(), cancelDipoles() and rebalance() all move them on purpose
    // over the next few lines, and an audit run against the moved set would be
    // reporting the pipeline doing its job.
    fieldIndex = cones->getIndices();

    if (interfaces && interfaces->multiMaterial() && options.prescribeInterfaceCones) {
        cones->prescribe(interfaces->prescription());
        status.interfaceConesPrescribed = interfaces->getReport().prescribedCones;
        status.interfacePrescriptionShift = cones->prescriptionShift();
        if (status.interfacePrescriptionShift > 0) {
            std::ostringstream oss;
            oss << "The interface network moved " << status.interfacePrescriptionShift
                << " index unit(s) from where the cross field put them.";
            status.messages.push_back("Stage 1: " + oss.str());
        }
    }

    if (interfaces && interfaces->multiMaterial() && options.alignFieldToInterfaces &&
        cancelDipoles) {
        const int nV = static_cast<int>(mesh->vertices.size());
        std::vector<int> region(nV, -1);
        for (int v = 0; v < nV; ++v) region[v] = interfaces->regionAt(v);
        status.coneDipoleUnits = cones->cancelDipoles(region);
        if (status.coneDipoleUnits > 0) {
            std::ostringstream oss;
            oss << "Cancelled " << status.coneDipoleUnits
                << " +1/-1 cone pair(s) the field put inside a single material region.";
            status.messages.push_back("Stage 1: " + oss.str());
        }
    }

    ConeSingularities::GaussBonnetReport gb = cones->gaussBonnet();
    if (!gb.admissible && options.autoRebalance) {
        const int moved = cones->rebalance();
        if (moved < 0) {
            status.messages.push_back(
                "Stage 1: could not restore Eq. (4): no boundary cone left within the allowed "
                "index range to take the residual.");
        }
        gb = cones->gaussBonnet();
    }
    status.conesAdmissible = gb.admissible;
    status.interiorCones = static_cast<int>(cones->interiorCones().size());
    status.boundaryCones = static_cast<int>(cones->boundaryCones().size());
    for (const std::string &m : gb.messages) status.messages.push_back("Stage 1: " + m);

    // Eq. (4) is the solvability condition of the Ricci flow's Newton system in
    // Pipeline A. There is no Newton system here -- and it is no less binding:
    // sum I(v) = 4 chi(S) is the condition for a *seamless map with these
    // holonomies to exist at all*, which is the same statement the flow's
    // consistency condition was standing in for.
    if (!status.conesAdmissible) {
        status.messages.push_back(
            "Stopping before the cut: sum I(v) != 4 chi(S), so no seamless map with these "
            "holonomies exists and there is nothing for the integration to find.");
        return false;
    }

    if (interfaces && interfaces->multiMaterial()) {
        interfaces->balance(cones->getIndices());
        const Interfaces::Report &fr = interfaces->getReport();
        status.interfaceNodes = fr.nodes;
        status.interfaceBranches = fr.branches;
        status.interfaceIllPosedNodes = fr.illPosedNodes;
        status.interfaceWorstSector = fr.worstSectorResidual;
        status.regions = fr.regions;
        status.regionsBalanced = fr.regionsBalanced;
        status.regionQuartersMoved = fr.quartersMoved;
        status.regionCornersInserted = fr.cornersInserted;
        status.worstRegionDeficit = fr.worstRegionDeficit;
        for (size_t i = balanceMessagesSeen; i < fr.messages.size(); ++i) {
            status.messages.push_back("Stage 0b: " + fr.messages[i]);
        }
    }

    ConeCut::Options cutOpts;
    cutOpts.conesToBoundary = options.coneCutsToBoundary;
    cutOpts.interfaceAvoidance = options.coneCutInterfaceAvoidance;
    cutter = std::make_unique<ConeCut>(mesh, *cones, cutOpts, interfaces.get());
    const ConeCut::Report &cr = cutter->getReport();
    status.cutIsDisk = cr.isDisk;
    status.allConesOnBoundary = cr.allConesOnBoundary;
    status.cutInteriorJunctions = cr.interiorJunctions;
    status.cutInterfaceEdges = cr.interfaceEdgesOnCut;
    status.cutInterfaceVertices = cr.interfaceVertsOnCut;
    status.cutInterfaceNodes = cr.interfaceNodesOnCut;
    for (const std::string &m : cr.messages) status.messages.push_back("Stage 2: " + m);

    if (!status.cutIsDisk || !status.allConesOnBoundary) {
        status.messages.push_back(
            "Stopping before the combing: Omega is not a disk with every cone on its boundary, "
            "so the branch of the field over it is not single valued and the BFS that selects "
            "it has no reason to be path-independent.");
        return false;
    }
    return true;
}

// ---------------------------------------------------------------------------
// inducedLengths()
//
// The lengths a map on Omega induces, re-indexed by the edges of S, which is the
// indexing LayoutEnergy::buildReference() reads. A seam edge has two images in
// Omega and they are the same length to the seam residual, because the
// transition holding them together is a rotation; the first found is taken.
//
// This is E1's reference metric in this pipeline, and the reason it is one is
// in the class comment: the field frame cannot be, because a rotation-invariant
// energy cannot see the rotation that carries the cones, and the model's own
// Euclidean geometry has no cones at all. A map that satisfies the seam
// constraints has them exactly, and that is the whole of what is needed.
// ---------------------------------------------------------------------------
std::vector<double> TORSION::inducedLengths(const Mesh &mesh, const ConeCut &cut,
                                            const std::vector<Point> &psi) {
    const Mesh &om = cut.getCutMesh();
    const auto &c2o = cut.getCutVertexToOriginal();
    std::vector<double> len(mesh.edges.size(), 0.0);
    if (psi.size() != om.vertices.size()) return len;

    std::unordered_map<MeshEdgeKey, int, MeshEdgeKeyHash> origEdge;
    origEdge.reserve(mesh.edges.size() * 2);
    for (int e = 0; e < static_cast<int>(mesh.edges.size()); ++e) {
        origEdge.emplace(MeshEdgeKey(mesh.edges[e][0], mesh.edges[e][1]), e);
    }
    for (int e = 0; e < static_cast<int>(om.edges.size()); ++e) {
        auto it = origEdge.find(MeshEdgeKey(c2o[om.edges[e][0]], c2o[om.edges[e][1]]));
        if (it == origEdge.end() || len[it->second] > 0.0) continue;
        len[it->second] = normP(psi[om.edges[e][1]] - psi[om.edges[e][0]]);
    }
    // An edge with no positive image -- a triangle the integration collapsed --
    // falls back to its Euclidean length, which is the reference having no cone
    // information about that one edge rather than having none at all.
    for (size_t e = 0; e < len.size(); ++e) {
        if (!(len[e] > 0.0)) {
            len[e] = normP(mesh.vertices[mesh.edges[e][1]] -
                           mesh.vertices[mesh.edges[e][0]]);
        }
    }
    return len;
}

// ---------------------------------------------------------------------------
// vertexFreedom()
//
// What Sec. 7.2a may do to each vertex of Omega without undoing an equality the
// integration imposed.
//
//   the seam        nothing. A seam vertex is one side of a pair held to
//                   phi_+ = R_k phi_- + t_k, and moving it without moving its
//                   partner breaks Q4, which Stage 6's barrier cannot restore
//                   any more than it can restore Q1.
//
//   an aligned      one coordinate. The chain was held to u = const or
//   chain           v = const; the vertex may slide along that line and not
//                   across it. A vertex where two chains of different axes meet
//                   -- which is a cone, and is where the layout turns -- has
//                   both coordinates held and cannot move at all.
//
//   everything      both. An interior vertex, or a boundary chain Sec. 6.4 left
//   else            free, is the integration's choice and nothing downstream
//                   has been told about it yet.
//
// `level` is the rung of Sec. 7.2a's ladder: 0 respects both, 1 lets the
// alignment go and keeps the seam, 2 lets both go. Anything above 0 hands its
// result to Sec. 6.5's projection to put back what it let go, which is what
// makes letting go of an equality a step rather than a loss.
// ---------------------------------------------------------------------------
// `usedAxis` is the axis the integration was solved with, empty when it was
// solved free. Static, and taking the cut and the axis rather than reading them
// off the object, because the viewer climbs this ladder too. See TORSION.hxx.
std::vector<unsigned char> TORSION::vertexFreedom(const ConeCut &cut,
                                                  const std::vector<int> &usedAxis,
                                                  int level) {
    const Mesh &om = cut.getCutMesh();
    const Mesh &orig = cut.getOriginalMesh();
    const auto &c2o = cut.getCutVertexToOriginal();
    const int nV = static_cast<int>(om.vertices.size());
    std::vector<unsigned char> freedom(nV, kFree);

    std::unordered_map<MeshEdgeKey, int, MeshEdgeKeyHash> origEdge;
    origEdge.reserve(orig.edges.size() * 2);
    for (int e = 0; e < static_cast<int>(orig.edges.size()); ++e) {
        origEdge.emplace(MeshEdgeKey(orig.edges[e][0], orig.edges[e][1]), e);
    }
    if (level < 2) {
        // Every seam vertex is held here. A seam vertex *can* move, paired: it
        // moves under phi_+ = R_k phi_- + t together with its partner, which
        // keeps Q4 exact -- but only a caller that moves the pair together may
        // let it, and that caller frees the pairs itself (buildPsi0's
        // pairedFreedom). A caller moving vertices one at a time would break
        // the seam. seamPairing() says which can be paired; a cone tip, where
        // the two chains are one vertex and the rotation has only its fixed
        // point to offer, and a landing, where dS has the vertex as well, never
        // can.
        for (int e : om.boundaryEdges) {
            auto it = origEdge.find(MeshEdgeKey(c2o[om.edges[e][0]], c2o[om.edges[e][1]]));
            const bool fromBoundary = it != origEdge.end() && orig.isBoundaryEdge[it->second];
            if (fromBoundary) continue;   // a curve of dS; the alignment below has it
            freedom[om.edges[e][0]] = kFixed;
            freedom[om.edges[e][1]] = kFixed;
        }
    }
    if (level < 1 && !usedAxis.empty() && usedAxis.size() == om.edges.size()) {
        for (int e = 0; e < static_cast<int>(om.edges.size()); ++e) {
            const int hold = usedAxis[e];
            if (hold != 0 && hold != 1) continue;
            // hold == 0 is u held, which leaves v; hold == 1 is v held.
            const unsigned char keep = (hold == 0) ? kFreeInV : kFreeInU;
            for (int i = 0; i < 2; ++i) {
                unsigned char &f = freedom[om.edges[e][i]];
                f = static_cast<unsigned char>(f & keep);
            }
        }
    }
    return freedom;
}

// ---------------------------------------------------------------------------
// seamPairing()
//
// The two children of each seam vertex of Omega, and the quarter turn between
// their displacements.
//
// The arcs of G hold phi_+ = R_k phi_- + t, so a displacement d of the minus
// child is legal exactly when the plus child takes R_k d -- and then Q4 is not
// merely preserved but unchanged. That is what lets Sec. 7.2a move a tangle
// that sits *on* the seam, which is where the tangle at the tip of a cut sits,
// without letting go of the one property nothing downstream can repair.
//
// Two cases come back unpaired and are held instead:
//
//   the cone tip   where the plus and minus chains are the same vertex. The
//                  relation is then d = R_k d, whose only solution for k != 0
//                  is d = 0.
//
//   a landing      where a child of the seam is also a vertex of dS, so
//                  Sec. 6.4 has an opinion about it too. Two constraints on one
//                  displacement is more than this pass is written to carry.
//
// A vertex that turns up in two arcs with different partners is also left
// unpaired: with conesToBoundary the seam is a forest of disjoint slits and
// that does not arise, but --cut-to-graph makes junctions and it does.
// ---------------------------------------------------------------------------
void TORSION::seamPairing(const ConeCut &cut, const Immersion &scaffold,
                          std::vector<int> &mate, std::vector<int> &turn) {
    const Mesh &om = cut.getCutMesh();
    const int nV = static_cast<int>(om.vertices.size());
    mate.assign(nV, -1);
    turn.assign(nV, 0);

    const auto &pairs = scaffold.getSeamPairs();
    const auto &pairArc = scaffold.getSeamPairArc();
    const auto &arcs = scaffold.getArcs();
    std::vector<char> conflicted(nV, 0);

    auto link = [&](int minus, int plus, int k) {
        if (minus < 0 || plus < 0) return;
        if (minus == plus) { conflicted[minus] = 1; return; }
        for (int i = 0; i < 2; ++i) {
            const int a = i ? plus : minus;
            const int b = i ? minus : plus;
            const int kk = i ? ((4 - ((k % 4) + 4) % 4) % 4) : (((k % 4) + 4) % 4);
            if (mate[a] >= 0 && (mate[a] != b || turn[a] != kk)) conflicted[a] = 1;
            mate[a] = b;
            turn[a] = kk;
        }
    };
    for (size_t p = 0; p < pairs.size(); ++p) {
        const int a = pairArc[p];
        if (a < 0 || a >= static_cast<int>(arcs.size())) continue;
        link(pairs[p][2], pairs[p][0], arcs[a].k);
        link(pairs[p][3], pairs[p][1], arcs[a].k);
    }

    // A child that is also on dS: Sec. 6.4 has it, and one displacement cannot
    // answer to both.
    {
        const Mesh &orig = cut.getOriginalMesh();
        const auto &c2o = cut.getCutVertexToOriginal();
        std::unordered_map<MeshEdgeKey, int, MeshEdgeKeyHash> origEdge;
        origEdge.reserve(orig.edges.size() * 2);
        for (int e = 0; e < static_cast<int>(orig.edges.size()); ++e) {
            origEdge.emplace(MeshEdgeKey(orig.edges[e][0], orig.edges[e][1]), e);
        }
        for (int e : om.boundaryEdges) {
            auto it = origEdge.find(MeshEdgeKey(c2o[om.edges[e][0]], c2o[om.edges[e][1]]));
            if (it == origEdge.end() || !orig.isBoundaryEdge[it->second]) continue;
            conflicted[om.edges[e][0]] = 1;
            conflicted[om.edges[e][1]] = 1;
        }
    }
    for (int v = 0; v < nV; ++v) {
        if (!conflicted[v]) continue;
        if (mate[v] >= 0) { mate[mate[v]] = -1; turn[mate[v]] = 0; }
        mate[v] = -1;
        turn[v] = 0;
    }
    // A pairing has to be mutual after all that striking out.
    for (int v = 0; v < nV; ++v) {
        if (mate[v] >= 0 && mate[mate[v]] != v) { mate[v] = -1; turn[v] = 0; }
    }
}

// ---------------------------------------------------------------------------
// untangle()  --  Stage 4R, Sec. 7.2
//
//   1. Tutte embedding of Omega: bijective by theorem, so zero flips, and
//      useless for anything except being somewhere injective to start from.
//   2. min_phi sum_t A_t [ ||J (J*)^-1||_F^2 + ||J* J^-1||_F^2 ] + mu E4(phi)
//   3. raise mu over a few outer steps.
//
// Step 2 is not a new solver: E1 of Stage 6 is already the bracket above, E4 is
// already the seam term, and the closed-form flip cap in its line search is what
// makes the pass work where the least-squares solve does not. The map it returns
// is locally injective by construction of the step, whether or not the target it
// is fitting is realisable --
//
//     Fitting a map to a per-triangle target Jacobian under a barrier is well
//     posed *even when the target is unrealisable*. Non-integrability makes the
//     target inconsistent, not the fit ill-posed.
//
// -- so the pass is LayoutEnergy with lambda_2, lambda_3, lambda_5 and lambda_6
// switched off, started from the Tutte map instead of from psi_R, and step 3 is
// its ordinary penalty schedule.
//
// **Which reference.** The plan says Reference::Field here, on the argument that
// J* is the target and E1 against J* is the fit. It is not: J* is a scaled
// rotation, so E1 against it is E1 against the model's Euclidean geometry --
// see LayoutEnergy.hxx for the algebra and TestTORSION --ref-test for the
// measurement -- and that reference has no cones, which is the failure this
// whole pipeline has to avoid. What is used instead is the metric the *least-
// squares map* induces. It is the closest thing to the field's own answer that
// a rotation-invariant energy can be given: it carries the cone angles, because
// the seam constraints put them there exactly, and it carries the shape the
// integration settled on everywhere else. The orientation the frame also
// carries is not recoverable by E1 at all, whatever it is composed with, and
// where that matters -- Stage 7's cross-check that a cone's rays leave along
// the field's directions -- it has to be checked rather than energetically
// asked for.
//
// The lengths are read off the map even where it inverted: a flipped triangle
// still has three positive side lengths, they are the ones the integration
// produced, and there are single figures of them.
//
// Returns an empty vector when there was nothing to start from.
// ---------------------------------------------------------------------------
std::vector<Point> TORSION::untangle(const MapStage &in, const Options &options,
                                     Status &status, Psi0 &out) {
    const Mesh &omega = in.cut.getCutMesh();

    TutteEmbedding::Options topts;
    topts.targetEdge = options.targetEdge;
    out.tutte = std::make_unique<TutteEmbedding>(omega, topts);
    const TutteEmbedding::Report &tr = out.tutte->getReport();
    status.tutteValid = tr.valid;
    for (const std::string &m : tr.messages) status.messages.push_back("Stage 4R: " + m);
    if (!tr.valid) {
        status.messages.push_back(
            "Stage 4R: the untangling has nothing locally injective to start from, so it was "
            "not run; psi_0 stands as the integration returned it.");
        return {};
    }

    // A second Immersion over the Tutte map, for the arcs and the seam pairing
    // the energy needs. The combed angles go in again, so the k it holds are
    // the matchings' and not a Procrustes fit of a Tutte map, which would be
    // meaningless.
    std::unique_ptr<Immersion> start;
    try {
        start = std::make_unique<Immersion>(in.cut, in.cones, out.tutte->getUV(),
                                            in.frames.fieldEdgeLengths(),
                                            in.frames.combedAngle());
    } catch (const std::exception &e) {
        status.messages.push_back(std::string("Stage 4R: could not wrap the Tutte map: ") + e.what());
        return {};
    }

    // No connectivity to seed from a Tutte map, so E5 is not asked for; a
    // separatrix trace of a map that means nothing is the most expensive thing
    // this pass could be made to do.
    //
    // The *labels*, on the other hand, are worth having and cannot be read off
    // the Tutte map either -- a circle has no boundary alignment. They are read
    // off psi_0 instead, where Sec. 6.4 has just put them exactly: every chain
    // of dS - G and every branch of the interface network is on a coordinate
    // line there, and those are the lines the untangled map has to come back
    // to. Letting E2 and E3 carry them through this pass rather than leaving
    // Stage 6 to rediscover them is the difference between one inverted face
    // out of eighteen thousand costing a re-solve and costing a continuation.
    SubdomainLabels::Options lopts;
    lopts.seedTopoConstraints = false;
    lopts.interfaceCorners = false;
    lopts.propagateInterfaceLabels = options.propagateInterfaceLabels;
    lopts.seamTurnInterfaceLabels = options.seamTurnInterfaceLabels;
    SubdomainLabels startLabels(*start, lopts, in.interfaces);
    if (options.untangleAlignWeight > 0.0) startLabels.relabel(out.integratedMap);

    LayoutEnergy::Options eopts;
    eopts.reference = LayoutEnergy::Reference::Induced;
    eopts.referenceLengths = inducedLengths(in.mesh, in.cut, out.integratedMap);
    // mu is E4's penalty and nothing else is switched on, so it is set through
    // lambdaFactor rather than through lambdaInit: run() floors lambda_4 at
    // lambda_1 (Q4 is not free to trade away in the main continuation), and
    // that floor would swallow a starting mu below 1. The schedule then raises
    // it by lambdaGrowth per outer step, which is step 3 of Sec. 7.2.
    eopts.lambdaInit = 1.0;
    eopts.lambdaFactor[4] = options.untangleSeamWeight;
    eopts.lambdaGrowth = options.lambdaGrowth;
    eopts.outerSteps = options.untangleOuterSteps;
    eopts.innerIterations = options.untangleInnerIterations;
    eopts.alternateReference = false;
    eopts.relabel = false;
    // E1, E4, and -- since Sec. 6.4 -- E2 and E3, which are what carry the
    // alignment across the pass. E5 and E6 stay off: both are statements about
    // a layout and this pass is not producing one, it is producing something
    // injective with the field's directions in it for Stage 6 to make a layout
    // out of.
    eopts.lambdaFactor[2] = options.untangleAlignWeight;
    eopts.lambdaFactor[3] = options.untangleAlignWeight;
    eopts.lambdaFactor[5] = 0.0;
    eopts.lambdaFactor[6] = 0.0;


    LayoutEnergy fit(*start, startLabels, eopts);
    fit.run();
    const LayoutEnergy::Report &fr = fit.getReport();
    status.untangleRan = true;
    status.untangleOuterSteps = fr.outerSteps;
    status.untangleFlippedFaces = fr.invertedTriangles;
    status.untangleSeamResidual = fr.maxSeamResidual;
    status.layoutReferenceLeftHanded =
        std::max(status.layoutReferenceLeftHanded, fr.leftHandedFrames);
    for (const std::string &m : fr.messages) status.messages.push_back("Stage 4R: " + m);

    if (fr.invertedTriangles > 0) {
        std::ostringstream oss;
        oss << "The untangling pass came back with " << fr.invertedTriangles
            << " inverted triangle(s), which means the line search's flip cap was defeated -- "
            << "normally by a reference triangle that was already degenerate. Sec. 7.3's "
            << "knobs, in order: field smoothness near the cones, refinement where the "
            << "per-triangle fit residual is largest, a lower h there.";
        status.messages.push_back("Stage 4R: " + oss.str());
    } else {
        std::ostringstream oss;
        oss << "Untangled: from the Tutte embedding, " << fr.outerSteps
            << " outer step(s) of target-Jacobian fitting under the barrier left every "
            << "triangle positively oriented, with the seam " << fr.maxSeamResidual
            << " of the extent from exact. E4 is what Stage 6 raises from here.";
        status.messages.push_back("Stage 4R: " + oss.str());
    }
    return fit.getUV();
}

// ---------------------------------------------------------------------------
// pullOntoAlignment()  --  Sec. 7.2b
//
// When Stage 4F could not keep the alignment -- no solve that held it came out
// of Sec. 7.2a's ladder -- what it hands on is an injective map that is only as
// aligned as the field is integrable. Stage 6 would get it aligned in the end:
// Q3 and the features reach 1e-6 on tunnel from exactly such a start. That is
// not where the damage is. Stage 5 seeds Gamma_topo from the separatrices of
// psi_0 *before* Stage 6 runs, pairing the cones psi_0's traced curves nearly
// join, and on an unaligned psi_0 those curves run wherever the misaligned
// boundary sends them. Every constraint seeded from them is one Stage 6 is then
// bound to satisfy, and Gamma_topo can be added to and never taken back. So an
// unaligned psi_0 is not a slower start to the same layout; it is a different
// and worse layout.
//
// The remedy is the one the continuation already is, run for a different
// purpose. The map is injective and the barrier of Sec. 3.3 keeps it so; E2
// and E3, on labels read off the solve that held the whole alignment, pull
// dS and the interfaces onto their axes; E4 keeps the seam; E1 against the flat
// cone metric keeps the rest in shape. No E5 and no E6: this is psi_0, not a
// layout. What comes out is aligned to the continuation's tolerance, and Sec.
// 6.5's projection then makes it exact where that inverts nothing.
// ---------------------------------------------------------------------------
std::vector<Point> TORSION::pullOntoAlignment(const MapStage &in, const Options &options,
                                              Status &status, const std::vector<Point> &start,
                                              const std::vector<Point> &labelMap,
                                              const std::vector<int> &axis) {
    const Mesh &omega = in.cut.getCutMesh();
    if (start.size() != omega.vertices.size() || labelMap.size() != omega.vertices.size()) return {};

    std::unique_ptr<Immersion> from;
    try {
        from = std::make_unique<Immersion>(in.cut, in.cones, start, in.frames.fieldEdgeLengths(),
                                           in.frames.combedAngle());
    } catch (const std::exception &e) {
        status.messages.push_back(std::string("Stage 4R: could not wrap the map to pull: ") + e.what());
        return {};
    }
    if (from->getReport().flippedFaces > 0) return {};

    SubdomainLabels::Options lopts;
    lopts.seedTopoConstraints = false;
    lopts.interfaceCorners = false;
    lopts.propagateInterfaceLabels = options.propagateInterfaceLabels;
    lopts.seamTurnInterfaceLabels = options.seamTurnInterfaceLabels;
    SubdomainLabels labels(*from, lopts, in.interfaces);
    labels.relabel(labelMap);

    LayoutEnergy::Options eopts;
    eopts.reference = LayoutEnergy::Reference::Induced;
    eopts.referenceLengths = in.coneMetric ? in.coneMetric->edgeLengths()
                                        : inducedLengths(in.mesh, in.cut, start);
    eopts.lambdaInit = options.lambdaInit;
    eopts.lambdaGrowth = options.lambdaGrowth;
    eopts.outerSteps = options.pullOuterSteps;
    eopts.innerIterations = options.innerIterations;
    eopts.alternateReference = false;
    eopts.relabel = false;
    eopts.lambdaFactor[2] = 1.0;
    eopts.lambdaFactor[3] = 1.0;
    eopts.lambdaFactor[4] = options.lambdaSeamFactor;
    eopts.lambdaFactor[5] = 0.0;
    eopts.lambdaFactor[6] = 0.0;

    LayoutEnergy fit(*from, labels, eopts);
    fit.run();
    const LayoutEnergy::Report &fr = fit.getReport();
    status.pullRan = true;
    status.pullBoundaryResidual = fr.maxBoundaryResidual;
    status.pullFeatureResidual = fr.maxFeatureResidual;
    if (fr.invertedTriangles > 0) return {};
    status.pullKept = true;
    std::vector<Point> pulled = fit.getUV();

    // Sec. 6.5 onto the exact alignment, kept only if it inverts nothing.
    if (!axis.empty()) {
        FieldIntegration::Options po;
        po.regularisation = options.integrationRegularisation;
        po.alignAxis = axis;
        po.targetJacobian = jacobianOf(omega, pulled);
        FieldIntegration proj(in.cut, in.frames, in.scaffold, po);
        if (proj.getReport().solved && proj.getReport().flippedFaces == 0) {
            status.pullProjected = true;
            return proj.getUV();
        }
    }
    return pulled;
}

// ---------------------------------------------------------------------------
// runLayoutStages()  --  Stages 5 to 8 at one Gamma_topo seeding tolerance
//
// MERIDIAN::runLayoutStages with this pipeline's Stage 6 in the middle of it.
// Separated out for the same reason and driven the same way; see
// Options::topoNearMissRetry.
//
// Returns false when a stage stopped the pipeline, in which case run() is done
// and there is nothing to retry.
// ---------------------------------------------------------------------------
bool TORSION::runLayoutStages(double nearMiss, bool seed) {
    // --- Stage 5: the subdomain labelling, MERIDIAN's unchanged -----------
    SubdomainLabels::Options lopts;
    lopts.seedTopoConstraints = seed;
    lopts.nearMissTolerance = nearMiss;
    lopts.seedSelfReturns = options.seedSelfReturns;
    lopts.seedAllConnections = options.seedAllConnections;
    lopts.maxTraceSteps = options.separatrixMaxSteps;
    lopts.interfaceCorners = options.interfaceCorners;
    lopts.propagateInterfaceLabels = options.propagateInterfaceLabels;
    lopts.seamTurnInterfaceLabels = options.seamTurnInterfaceLabels;
    labels = std::make_unique<SubdomainLabels>(*immersion, lopts, interfaces.get());
    const SubdomainLabels::Report &lr = labels->getReport();
    status.boundaryEdgesU = lr.boundaryEdgesU;
    status.boundaryEdgesV = lr.boundaryEdgesV;
    status.featureChains = lr.featureChains;
    status.interfaceCorners = lr.interfaceCorners;
    status.interfaceLabelsCorrected = lr.featureLabelsCorrected;
    status.interfaceChainsSeamFlipped = lr.featureChainsSeamFlipped;
    status.topoPaths = lr.topoPaths;
    status.topoSelfReturns = lr.topoSelfReturns;
    status.topoExtraPerPair = lr.topoExtraPerPair;
    for (const std::string &m : lr.messages) status.messages.push_back("Stage 5: " + m);

    // --- E1's reference metric, C4 ----------------------------------------
    //
    // The lengths psi_0 induces on Omega, re-indexed by the edges of S, which is
    // the indexing buildReference() reads. A seam edge has two images in Omega
    // and they are the same length to the seam residual, because the transition
    // holding them together is a rotation; either one will do and the first is
    // taken.
    std::vector<double> inducedLen;
    if (options.reference == Options::Reference::Induced) {
        inducedLen = inducedLengths(*mesh, *cutter, immersion->getUV());
    } else if (options.reference == Options::Reference::Cone) {
        // Already built, before Stage 3F. If it failed to come out a metric at
        // all, fall back to the lengths psi_0 induces rather than to the
        // Euclidean ones, which have no cones.
        if (!coneMetric) {
            inducedLen = inducedLengths(*mesh, *cutter, immersion->getUV());
            status.messages.push_back(
                "Stage 6: the cone metric was not realisable, so E1 measures against the "
                "lengths psi_0 induces instead. C4 below is what that is worth here.");
        }
    } else if (options.reference == Options::Reference::Ricci) {
        // Pipeline A's flat cone metric, computed here for one purpose only:
        // to be E1's notion of undistorted. psi_0 still comes from the
        // integration. This is the known-good yardstick the other two are
        // measured against, and it is the one setting under which this pipeline
        // still needs the flow.
        referenceFlow = std::make_unique<RicciFlow>(mesh, *cones);
        referenceFlow->solve();
        for (const std::string &m : referenceFlow->getReport().messages) {
            status.messages.push_back("Stage 6 (reference flow): " + m);
        }
    }

    // What that reference actually does at the cones, measured before the
    // continuation is allowed to use it. C4 in one number: a reference with no
    // cones reports (pi/2)|I| here, and the plan's whole warning is that E1 and
    // Q2 are then contradictory statements about the same vertex.
    {
        const std::vector<double> *refLen = nullptr;
        std::vector<double> ricciLen;
        if (options.reference == Options::Reference::Induced) {
            refLen = &inducedLen;
        } else if (options.reference == Options::Reference::Cone) {
            refLen = coneMetric ? &coneMetric->edgeLengths() : &inducedLen;
        } else if (options.reference == Options::Reference::Ricci && referenceFlow) {
            ricciLen = referenceFlow->originalEdgeLengthsCompleted(nullptr);
            refLen = &ricciLen;
        }
        std::vector<double> euclid;
        if (!refLen) {
            // Field and Euclidean are the same reference; measure the one they
            // both are.
            euclid.assign(mesh->edges.size(), 0.0);
            for (size_t e = 0; e < euclid.size(); ++e) {
                euclid[e] = normP(mesh->vertices[mesh->edges[e][1]] -
                                  mesh->vertices[mesh->edges[e][0]]);
            }
            refLen = &euclid;
        }
        std::vector<double> sum(mesh->vertices.size(), 0.0);
        for (int t = 0; t < static_cast<int>(mesh->triangles.size()); ++t) {
            const double a = (*refLen)[mesh->triangleEdges[t][0]];
            const double b = (*refLen)[mesh->triangleEdges[t][1]];
            const double c = (*refLen)[mesh->triangleEdges[t][2]];
            const double l[3] = {a, b, c};   // l[i] joins tri[i] and tri[i+1]
            for (int i = 0; i < 3; ++i) {
                const double p = l[i], q = l[(i + 2) % 3], r = l[(i + 1) % 3];
                if (!(p > 0.0) || !(q > 0.0)) continue;
                double cosA = (p * p + q * q - r * r) / (2.0 * p * q);
                cosA = std::max(-1.0, std::min(1.0, cosA));
                sum[mesh->triangles[t][i]] += std::acos(cosA);
            }
        }
        const std::vector<int> &idx = cones->getIndices();
        for (const auto &c : cones->getCones()) {
            const double full = mesh->isBoundaryVertex[c.vertex] ? M_PI : 2.0 * M_PI;
            const double want = full - M_PI_2 * idx[c.vertex];
            status.referenceConeResidual =
                std::max(status.referenceConeResidual, std::fabs(sum[c.vertex] - want));
        }
    }

    // --- Stage 6: the layout energies -------------------------------------
    LayoutEnergy::Options eopts;
    eopts.lambdaInit = options.lambdaInit;
    eopts.lambdaGrowth = options.lambdaGrowth;
    eopts.outerSteps = options.outerSteps;
    eopts.innerIterations = options.innerIterations;
    eopts.alternateReference = false;
    eopts.relabel = options.relabelBetweenSteps;
    eopts.lagInterfaceScales = options.lagInterfaceScales;
    switch (options.reference) {
        case Options::Reference::Cone:
            eopts.reference = LayoutEnergy::Reference::Induced;
            eopts.referenceLengths = coneMetric ? coneMetric->edgeLengths() : inducedLen;
            break;
        case Options::Reference::Induced:
            eopts.reference = LayoutEnergy::Reference::Induced;
            eopts.referenceLengths = inducedLen;
            break;
        case Options::Reference::Field:
            eopts.reference = LayoutEnergy::Reference::Field;
            eopts.fieldFrames = frames->frames();
            break;
        case Options::Reference::Euclidean:
            eopts.reference = LayoutEnergy::Reference::Euclidean;
            break;
        case Options::Reference::Ricci:
            eopts.reference = LayoutEnergy::Reference::Induced;
            eopts.referenceLengths = referenceFlow
                ? referenceFlow->originalEdgeLengthsCompleted(nullptr)
                : std::vector<double>();
            break;
    }
    // Sec. 9's closing note. E2 and E3 start an order of magnitude lower
    // because the field is already aligned to dS and to the interfaces, so
    // their residuals start small and over-penalising them early only fights
    // E1; E4 starts higher because the untangling introduced seam error, and
    // because -- exactly as in Pipeline A -- Q4 and Q2 are the same statement
    // in different units and Q4 is not free to trade away.
    eopts.lambdaFactor[2] = options.lambdaAlignmentFactor;
    eopts.lambdaFactor[3] = options.lambdaAlignmentFactor;
    eopts.lambdaFactor[4] = options.lambdaSeamFactor;
    layout = std::make_unique<LayoutEnergy>(*immersion, *labels, eopts);
    status.layoutRan = layout->run();
    {
        const LayoutEnergy::Report &er = layout->getReport();
        status.layoutInjective = er.injective;
        status.layoutConstrained = er.constraintsMet;
        status.outerStepsTaken = er.outerSteps;
        status.layoutValid = er.valid;
        status.layoutConeAngleResidual = er.maxConeAngleResidual;
        status.layoutRegularAngleResidual = er.maxRegularAngleResidual;
        status.layoutReferenceLeftHanded =
            std::max(status.layoutReferenceLeftHanded, er.leftHandedFrames);
        status.interfaceResidual = er.maxInterfaceResidual;
        status.interfaceCornerChanges = er.interfaceCornerChanges;
        status.interfacesAligned = er.interfaceCornerChanges == 0 &&
                                   er.maxInterfaceResidual < 1e-3;
        for (const std::string &m : er.messages) status.messages.push_back("Stage 6: " + m);
    }

    if (!options.runSeparatrices) return false;

    // --- Stage 7 and Sec. 3.3's repair, MERIDIAN's unchanged --------------
    Separatrices::Options sopts;
    sopts.coneSnapTolerance = options.separatrixSnap;
    sopts.maxSteps = options.separatrixMaxSteps;
    sopts.coneSnapRings = options.separatrixSnapRings;
    sopts.detectCycles = options.separatrixDetectCycles;
    sopts.nearMissWindow = options.repairGapLimit;
    if (interfaces && interfaces->multiMaterial()) {
        sopts.extraEmitters = interfaces->emitterNodes();
    }

    MERIDIAN::RepairOptions ropts;
    ropts.passes = options.seedTopoConstraints ? options.repairPasses : 0;
    ropts.maxPerPass = options.repairMaxPerPass;
    ropts.gapLimit = options.repairGapLimit;
    ropts.lambdaBoost = options.repairLambdaBoost;
    ropts.outerSteps = options.repairOuterSteps;
    ropts.scoreArrangement = options.runArrangement && options.repairScoresArrangement;
    ropts.patience = options.repairPatience;
    ropts.arrangement.mergeTolerance = options.arrangementMerge;
    ropts.arrangement.cornerTolerance = options.arrangementCorner;
    ropts.arrangement.collapseTolerance = options.arrangementCollapse;
    ropts.arrangement.trimUnresolvedAtCrossings = options.arrangementTrim;

    MERIDIAN::RepairResult rep = MERIDIAN::traceAndRepair(
        *labels, *layout, sopts, ropts,
        [&](const std::string &m) { status.messages.push_back(m); });
    separatrices = std::move(rep.separatrices);
    status.repairPasses = rep.passesTaken;
    status.repairConstraintsAdded = rep.constraintsAdded;
    status.topoPaths = labels->getReport().topoPaths;
    status.topoSelfReturns = labels->getReport().topoSelfReturns;
    status.topoExtraPerPair = labels->getReport().topoExtraPerPair;

    {
        const LayoutEnergy::Report &er = layout->getReport();
        status.layoutInjective = er.injective;
        status.layoutConstrained = er.constraintsMet;
        status.outerStepsTaken = er.outerSteps;
        status.layoutValid = er.valid;
        status.layoutConeAngleResidual = er.maxConeAngleResidual;
        status.layoutRegularAngleResidual = er.maxRegularAngleResidual;
    }

    if (separatrices) {
        const Separatrices::Report &sr = separatrices->getReport();
        status.separatricesRan = true;
        status.separatrices = sr.emitted;
        status.separatricesToCone = sr.endedAtCone;
        status.separatricesToBoundary = sr.endedAtBoundary;
        status.separatricesUnresolved = sr.capped + sr.cycled + sr.stuck + sr.degenerate;
        status.separatricesNearMisses = sr.nearMisses + sr.grazes;
        status.q5Verified = sr.valid;
    }

    // --- Stage 8 ----------------------------------------------------------
    if (!separatrices || !options.runArrangement) return false;
    Arrangement::Options aopts;
    aopts.mergeTolerance = options.arrangementMerge;
    aopts.cornerTolerance = options.arrangementCorner;
    aopts.collapseTolerance = options.arrangementCollapse;
    aopts.trimUnresolvedAtCrossings = options.arrangementTrim;
    try {
        arrangement = std::make_unique<Arrangement>(*separatrices, *labels, aopts);
    } catch (const std::exception &e) {
        status.messages.push_back(std::string("Stage 8: could not be built: ") + e.what());
        return false;
    }
    const Arrangement::Report &ar = arrangement->getReport();
    status.arrangementRan = true;
    status.layoutNodes = ar.nodes;
    status.layoutArcs = ar.arcs;
    status.layoutPatches = ar.patches;
    status.layoutQuads = ar.simpleQuads;
    status.layoutCoverage = ar.areaCoverage;
    status.arrangementValid = ar.valid;
    for (const std::string &m : ar.messages) status.messages.push_back("Stage 8: " + m);
    return true;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// buildPsi0()  --  Stages 4F and 4R, from the frames to psi_0
//
// Everything between the combed frames and the map the Immersion is built on:
// Sec. 6.4's attempts at the alignment, Sec. 7.2a-c's ladder, Sec. 7.2's Tutte
// pass behind it, and Sec. 7.2b's pull back onto the alignment. Static, and
// taking what it reads as arguments, because the viewer drives Stage 4 itself
// and has to take the same decisions in the same order -- a second copy of this
// is a second thing to keep in step with it, and the copy it had did not keep
// up. Returns false where no psi_0 can be handed on; out.integration is then
// still whatever was solved, for the report.
// ---------------------------------------------------------------------------
bool TORSION::buildPsi0(const MapStage &in, const Options &options, Status &status, Psi0 &out) {
    // Sec. 6.4's three attempts. Holding the boundary and the interfaces
    // exactly is a constraint on the fit, and a constraint on a fit can invert
    // a triangle the free fit would not have -- the place it happens is a chain
    // whose staircase was overridden, because a step is a right angle the map
    // is being asked to unbend. So the alignment is tried in full, then on the
    // chains that needed no overriding, then not at all.
    //
    // **What decides between them is what Sec. 7.2a leaves, not the raw flip
    // count.** A handful of inverted faces in a cone's one ring is what a
    // least-squares fit of a non-integrable field produces with or without the
    // alignment, and the local pass clears it with the alignment held; the free
    // solve may invert nothing and still be the worse psi_0 by far, because
    // everything it did not hold Stage 6 has to find from an unaligned start.
    // Choosing on the raw count threw the alignment away on 20 of the 31
    // multi-material models for tangles of 2 to 30 faces, and on concrete the
    // difference is the whole layout: dropped, Stage 8 came back with 8590
    // patches of which 3203 were not quadrilaterals; kept and cleared at rung
    // 0, 5485 quadrilaterals and nothing else. So each attempt is handed to
    // the ladder in turn, and the first one that comes out of it as a legal
    // psi_0 is kept; only when none does is the one with the fewest flips
    // handed on to Sec. 7.2's Tutte pass, as before.
    auto solveWith = [&](const std::vector<int> &axis) {
        FieldIntegration::Options io;
        io.regularisation = options.integrationRegularisation;
        io.alignAxis = axis;
        return std::make_unique<FieldIntegration>(in.cut, in.frames, in.scaffold, io);
    };
    auto flipsOf = [](const FieldIntegration &fi) {
        return fi.getReport().solved ? fi.getReport().flippedFaces
                                     : std::numeric_limits<int>::max();
    };

    // Sec. 6.5's projection, as a step rather than as a stage: fit a map's own
    // Jacobian back under the seam and alignment equalities. Used after any
    // repair that was allowed to break one of them. Returns an empty map when
    // it did not solve or inverted something, with the count in `flips`.
    auto reproject = [&](const std::vector<Point> &m, const std::vector<int> &axis,
                         int &flips, std::vector<Point> *raw = nullptr) -> std::vector<Point> {
        FieldIntegration::Options po;
        po.regularisation = options.integrationRegularisation;
        po.alignAxis = axis;
        po.targetJacobian = jacobianOf(in.cut.getCutMesh(), m);
        FieldIntegration proj(in.cut, in.frames, in.scaffold, po);
        flips = proj.getReport().flippedFaces;
        if (raw && proj.getReport().solved) *raw = proj.getUV();
        if (!proj.getReport().solved || proj.getReport().flippedFaces > 0) return {};
        return proj.getUV();
    };

    // Sec. 7.2a, on a ladder. Each rung frees one more of the things the
    // integration fixed, and every rung above the first hands what it produced
    // to the projection to put them back:
    //
    //   0   the seam and Sec. 6.4's alignment both held. Nothing to restore, so
    //       nothing can go wrong restoring it.
    //   1   the alignment let go, the seam still held. Q3 comes back from the
    //       projection, or -- if the projection inverts -- Stage 6 has it to
    //       reach, which is where it started.
    //   2   both let go. Q4 is not recoverable by anything downstream, so this
    //       rung is kept only if the projection takes.
    //
    // The ladder exists because the tangles that survive rung 0 are the ones
    // at a cone, and a cone is exactly the vertex the alignment pins in both
    // coordinates. Rung 1 is what unpins it.
    struct Ladder {
        std::vector<Point> map;        // empty: the ladder did not clear it
        int level = -1;
        int leastLeft = 0;             // the fewest inverted faces any rung left
        std::string how;
        bool reprojectionRan = false;
        bool reprojectionKept = false;
        int reprojectionFlips = 0;
        int windingRejects = 0;        // rungs that cleared it by wrapping a star
    };
    // Q2's angle sum at every vertex of S: 2 pi - (pi/2) I inside, seam
    // vertices included, and pi - (pi/2) I on dS. A repair is checked against
    // these and not against the map it started from, whose signed angle sums
    // are a whole turn out wherever it had a triangle turned over.
    const std::vector<double> prescribedAngle = [&]() {
        std::vector<double> want(in.mesh.vertices.size(), 0.0);
        const std::vector<int> &I = in.cones.getIndices();
        for (size_t o = 0; o < want.size(); ++o) {
            const double full = in.mesh.isBoundaryVertex[o] ? M_PI : 2.0 * M_PI;
            want[o] = full - M_PI_2 * (o < I.size() ? I[o] : 0);
        }
        return want;
    }();
    auto angleSumsHold = [&](const std::vector<Point> &m) {
        const std::vector<double> now =
            angleSumsOnS(in.cut.getCutMesh(), in.cut.getCutVertexToOriginal(), m,
                         static_cast<int>(in.mesh.vertices.size()));
        for (size_t o = 0; o < now.size(); ++o) {
            if (std::fabs(now[o] - prescribedAngle[o]) > M_PI) return false;
        }
        return true;
    };
    // The seam pairing, and the freedom it gives: at rungs 0 and 1 each child of
    // a slit that seamPairing() can pair is freed to move with its partner,
    // unless at rung 0 the alignment would have held either of them too. At
    // rung 2 the seam is let go outright and nothing needs pairing.
    auto pairSeam = [&](const std::vector<int> &axis, int level,
                        std::vector<unsigned char> &freedom, std::vector<int> &pmate,
                        std::vector<int> &pturn) {
        pmate.clear();
        pturn.clear();
        if (level >= 2) return;
        seamPairing(in.cut, in.scaffold, pmate, pturn);
        const std::vector<unsigned char> alignOnly =
            (level == 0) ? alignmentFreedom(in.cut, axis)
                         : std::vector<unsigned char>(pmate.size(), kFree);
        for (size_t v = 0; v < pmate.size(); ++v) {
            const int w = pmate[v];
            if (w < 0) continue;
            if (alignOnly[v] == kFree && alignOnly[w] == kFree) {
                freedom[v] = kFree;
            } else {
                pmate[v] = -1;   // held to an axis as well: leave it held
            }
        }
    };
    // What the kernel pass leaves is a fold, and Sec. 7.2c is what opens a fold
    // -- within the same freedom, but with the seam paired rather than held: a
    // child of a slit moves with its partner under R_k, which is Q4 kept
    // exactly and the one freedom a fold round the tip of a slit needs. The tip
    // itself, and a landing, stay put. Returns what is left inverted.
    auto openFold = [&](std::vector<Point> &m, const std::vector<int> &axis, int level,
                        const std::vector<unsigned char> &freedom) -> int {
        if (options.regularisedUntangleRings <= 0) {
            int n = 0;
            const Mesh &cm = in.cut.getCutMesh();
            for (const Triangle &t : cm.triangles) {
                if (!(cross2(m[t[1]] - m[t[0]], m[t[2]] - m[t[0]]) > 0.0)) ++n;
            }
            return n;
        }
        std::vector<unsigned char> pairedFreedom = freedom;
        std::vector<int> pmate, pturn;
        pairSeam(axis, level, pairedFreedom, pmate, pturn);
        return untangleRegularised(in.cut.getCutMesh(), m, pairedFreedom, pmate, pturn,
                                   in.frames.frames(), in.cut.getCutVertexToOriginal(),
                                   prescribedAngle, options.regularisedUntangleRings,
                                   options.regularisedUntangleIterations);
    };
    // A projection that inverted a few faces is not the end of its rung: it put
    // the seam and the alignment back, and what it left is a tangle of exactly
    // the kind rungs 0 and 1 are for -- on a map whose cones are now wound the
    // way the rung above put them, which is what the integration's own map could
    // not offer. Rung 1 hands its result to the projection once more.
    auto settle = [&](const std::vector<Point> &projected,
                      const std::vector<int> &axis) -> std::vector<Point> {
        for (int lv = 0; lv < 2; ++lv) {
            if (lv == 1 && axis.empty()) continue;
            std::vector<Point> trial = projected;
            const std::vector<unsigned char> fz = vertexFreedom(in.cut, axis, lv);
            int left = relaxToKernel(in.cut.getCutMesh(), trial, fz, {}, {},
                                     options.localUntangleSweeps);
            if (left > 0) left = openFold(trial, axis, lv, fz);
            if (left != 0 || !angleSumsHold(trial)) continue;
            if (lv == 0) return trial;
            int f = 0;
            std::vector<Point> again = reproject(trial, axis, f);
            if (!again.empty()) return again;
        }
        return {};
    };
    auto climb = [&](const std::vector<Point> &start, const std::vector<int> &axis,
                     int flips) -> Ladder {
        Ladder out;
        out.leastLeft = flips;
        // A tangle of thousands of faces is not local to anything: the full
        // alignment on a model it contradicts, or a field that is not integrable
        // over whole regions. No rung clears one, and Sec. 7.2c would spend its
        // whole schedule over the whole map three times finding that out -- 2.5
        // minutes on singlemat/geom024's 6749. Such an attempt goes straight to
        // the chain releases, whose smaller tangles are where the ladder works.
        const int nFaces = static_cast<int>(in.cut.getCutMesh().triangles.size());
        const int limit = std::max(options.localUntangleFaceFloor,
                                   static_cast<int>(options.localUntangleFaceFraction * nFaces));
        if (flips > limit) return out;
        const std::vector<int> seamOnly;
        auto project = [&](const std::vector<Point> &m, const std::vector<int> &ax) {
            int f = 0;
            std::vector<Point> raw;
            std::vector<Point> put = reproject(m, ax, f, &raw);
            out.reprojectionRan = true;
            out.reprojectionFlips = f;
            if (put.empty() && !raw.empty() && f > 0) put = settle(raw, ax);
            if (!put.empty()) out.reprojectionKept = true;
            return put;
        };
        for (int level = 0; level < 3; ++level) {
            // Rung 1 frees the alignment, so with no alignment to free it is
            // rung 0 again -- the same freedom, the same pairing and the same
            // answer. Skipping it keeps the rung the message names honest about
            // what was actually tried.
            if (level == 1 && axis.empty()) continue;
            std::vector<Point> local = start;
            const std::vector<unsigned char> freedom = vertexFreedom(in.cut, axis, level);
            std::vector<unsigned char> kernelFreedom = freedom;
            std::vector<int> mate, turn;
            if (options.pairSeamInUntangle) pairSeam(axis, level, kernelFreedom, mate, turn);
            int left = relaxToKernel(in.cut.getCutMesh(), local, kernelFreedom, mate, turn,
                                     options.localUntangleSweeps);
            if (left > 0) left = openFold(local, axis, level, freedom);
            // Whatever cleared it, every vertex of S has to have the angle
            // sum Q2 prescribes, to within less than a turn, or a cone has
            // changed valence under the repair.
            if (left == 0 && !angleSumsHold(local)) {
                ++out.windingRejects;
                continue;
            }
            out.leastLeft = std::min(out.leastLeft, left);
            if (left != 0) continue;

            if (level == 0) {
                out.map = std::move(local);
                out.level = 0;
                out.how = "with the seam and the alignment held throughout";
                return out;
            }
            std::vector<Point> put = project(local, axis);
            if (!put.empty()) {
                out.map = std::move(put);
                out.level = level;
                out.how = "and Sec. 6.5's projection put the seam and the alignment back";
                return out;
            }
            if (level == 1) {
                // The seam was held throughout this rung, so the map is legal
                // as it stands; only Q3 is left for Stage 6.
                out.map = std::move(local);
                out.level = 1;
                out.how = "with Q3 left for Stage 6, the projection onto it having inverted";
                return out;
            }
            if (!axis.empty()) {
                // Both were let go, so Q4 has to come back or the map is not one
                // Stage 6 can be handed. The alignment is what made the
                // projection a large displacement; ask for the seam alone. With
                // no alignment in play the two projections are the same solve
                // and the first one has already failed.
                put = project(local, seamOnly);
                if (!put.empty()) {
                    out.map = std::move(put);
                    out.level = 2;
                    out.how = "and Sec. 6.5's projection put the seam back, Q3 with the "
                              "alignment having been too far to reach";
                    return out;
                }
            }
        }
        return out;
    };

    const std::vector<int> noAxis;
    out.usedAxis = options.alignInIntegration ? in.frames.alignmentAxis() : noAxis;
    out.integration = solveWith(out.usedAxis);

    Ladder ladder;
    bool ladderRan = false;
    // The solve that held the whole alignment, for Sec. 7.2b to read its
    // labels off if the one kept here does not.
    std::vector<Point> fullAlignedMap;
    std::vector<int> fullAxis = out.usedAxis;
    if (out.integration->getReport().solved) fullAlignedMap = out.integration->getUV();
    const bool canClimb = options.untangle && options.localUntangle;

    if (options.alignInIntegration && options.alignmentFallback &&
        out.integration->getReport().alignedEdges > 0 && flipsOf(*out.integration) > 0) {
        const int firstFlips = flipsOf(*out.integration);
        const double strain = out.integration->getReport().alignmentStrain;

        struct Attempt {
            std::unique_ptr<FieldIntegration> fit;
            std::vector<int> axis;
            const char *name;
            Ladder ladder;
            bool climbed = false;
        };
        std::vector<Attempt> attempts;
        int chosen = -1;
        auto tryAttempt = [&](size_t i) {
            Attempt &a = attempts[i];
            if (!a.fit) a.fit = solveWith(a.axis);
            const int flips = flipsOf(*a.fit);
            if (flips == std::numeric_limits<int>::max()) return false;
            if (flips == 0) return true;
            if (!canClimb || !options.alignmentChooseByLadder) return false;
            a.ladder = climb(a.fit->getUV(), a.axis, flips);
            a.climbed = true;
            return !a.ladder.map.empty();
        };
        attempts.push_back({std::move(out.integration), out.usedAxis, "the full alignment", {}, false});
        if (tryAttempt(0)) chosen = 0;

        // Letting go of the whole alignment because one corner of it is out of
        // reach is the answer this used to give, and it is out of all
        // proportion: on tunnel the one thing wrong is a stratum whose right
        // side is 0.04 of the model long, 25 times shorter than its left, and
        // the frame and its neighbours' axes disagree by a factor of eight about
        // how long its image is. So the chains the tangle actually touches are
        // released -- every held edge with an end on an inverted face, whole
        // chain at a time, since a chain is the unit the axis was decided in --
        // and the rest are solved for again. A few rounds, each releasing only
        // what the last solve's tangle reached; where a tangle touches no
        // chain at all, its neighbourhood is widened a ring at a time before
        // concluding that the alignment is not what inverted it.
        const std::vector<int> &chainOf = in.frames.alignmentChain();
        const Mesh &om = in.cut.getCutMesh();
        int releasedChains = 0;
        for (int round = 0; chosen < 0 && round < options.alignmentReleaseRounds &&
                            chainOf.size() == om.edges.size(); ++round) {
            const Attempt &prev = attempts.back();
            if (!prev.fit || !prev.fit->getReport().solved) break;
            const std::vector<int> &axis = prev.axis;
            if (axis.size() != om.edges.size()) break;
            const std::vector<double> &ratio = prev.fit->areaRatios();
            std::vector<char> near(om.vertices.size(), 0);
            for (size_t t = 0; t < om.triangles.size() && t < ratio.size(); ++t) {
                if (ratio[t] > 0.0) continue;
                for (int i = 0; i < 3; ++i) near[om.triangles[t][i]] = 1;
            }
            std::vector<char> release;
            int hit = 0;
            for (int ring = 0; ring < 4 && hit == 0; ++ring) {
                if (ring > 0) {
                    std::vector<char> grown = near;
                    for (const Triangle &tri : om.triangles) {
                        if (!near[tri[0]] && !near[tri[1]] && !near[tri[2]]) continue;
                        for (int i = 0; i < 3; ++i) grown[tri[i]] = 1;
                    }
                    near.swap(grown);
                }
                release.assign(static_cast<size_t>(in.frames.getReport().alignmentChains), 0);
                for (size_t e = 0; e < om.edges.size(); ++e) {
                    if (axis[e] < 0 || chainOf[e] < 0 ||
                        chainOf[e] >= static_cast<int>(release.size())) continue;
                    if (!near[om.edges[e][0]] && !near[om.edges[e][1]]) continue;
                    if (!release[chainOf[e]]) { release[chainOf[e]] = 1; ++hit; }
                }
            }
            if (hit == 0) break;
            std::vector<int> next = axis;
            for (size_t e = 0; e < om.edges.size(); ++e) {
                if (next[e] >= 0 && chainOf[e] >= 0 &&
                    chainOf[e] < static_cast<int>(release.size()) && release[chainOf[e]]) {
                    next[e] = -1;
                }
            }
            releasedChains += hit;
            attempts.push_back({nullptr, std::move(next),
                                "the alignment less the chains the tangle touched", {}, false});
            if (tryAttempt(attempts.size() - 1)) chosen = static_cast<int>(attempts.size()) - 1;
        }
        status.alignmentChainsReleased = releasedChains;

        if (chosen < 0 && in.frames.getReport().alignmentSteppedChains > 0) {
            attempts.push_back({nullptr, in.frames.strictAlignmentAxis(),
                                "the alignment without the overridden chains", {}, false});
            if (tryAttempt(attempts.size() - 1)) chosen = static_cast<int>(attempts.size()) - 1;
        }
        if (chosen < 0) {
            attempts.push_back({nullptr, noAxis, "no alignment at all", {}, false});
            if (tryAttempt(attempts.size() - 1)) chosen = static_cast<int>(attempts.size()) - 1;
        }
        // None of them came out of the ladder: the old rule, fewest flips and
        // the earliest on a tie, and Sec. 7.2's Tutte pass behind it.
        if (chosen < 0) {
            int best = std::numeric_limits<int>::max();
            for (size_t i = 0; i < attempts.size(); ++i) {
                if (attempts[i].fit && flipsOf(*attempts[i].fit) < best) {
                    best = flipsOf(*attempts[i].fit);
                    chosen = static_cast<int>(i);
                }
            }
            if (chosen < 0) chosen = 0;
        }
        ladderRan = attempts[chosen].climbed;
        ladder = std::move(attempts[chosen].ladder);

        if (chosen > 0) {
            std::ostringstream oss;
            oss << "The full alignment inverted " << firstFlips << " face(s) at a strain of "
                << strain << " relative"
                << (canClimb && options.alignmentChooseByLadder
                        ? ", more than Sec. 7.2a could clear with it held,"
                        : ",")
                << " so Sec. 6.4 fell back to " << attempts[chosen].name;
            if (releasedChains > 0) {
                oss << " (" << releasedChains << " of " << in.frames.getReport().alignmentChains
                    << " chain(s) released over " << (attempts.size() - 1) << " round(s))";
            }
            oss << ", which inverts " << flipsOf(*attempts[chosen].fit) << ". The strain is the "
                << "field disagreeing with the cone set about where the boundary turns; where it "
                << "is large the remedy is in Stage 1 and not here.";
            status.messages.push_back("Stage 4F: " + oss.str());
            status.alignmentWasDropped = true;
        } else if (ladderRan && !ladder.map.empty()) {
            std::ostringstream oss;
            oss << "The full alignment inverted " << firstFlips << " face(s) at a strain of "
                << strain << " relative, and was kept: Sec. 7.2a clears them at rung "
                << ladder.level << ", which is a better psi_0 than any solve that lets the "
                << "alignment go and leaves Stage 6 to find it again.";
            status.messages.push_back("Stage 4F: " + oss.str());
        }
        out.integration = std::move(attempts[chosen].fit);
        out.usedAxis = std::move(attempts[chosen].axis);
    }

    const FieldIntegration::Report &ir = out.integration->getReport();
    status.integrationRan = true;
    status.integrationAlignedEdges = ir.alignedEdges;
    status.integrationAlignResidual = ir.maxAlignResidual;
    status.integrationAlignStrain = ir.alignmentStrain;
    status.integrationSolved = ir.solved;
    status.integrationSeamResidual = ir.maxSeamResidual;
    status.integrationFlippedFaces = ir.flippedFaces;
    status.integrationFlippedAreaFraction =
        (ir.totalArea > 0.0) ? ir.flippedArea / ir.totalArea : 0.0;
    status.integrationMinAreaRatio = ir.minAreaRatio;
    status.integrationMaxFitResidual = ir.maxFitResidual;
    status.integrationMeanFitResidual = ir.meanFitResidual;
    status.integrationFlipsAtCones = ir.flipsAdjacentToCone;
    status.integrationNearestFlipToCone = ir.nearestFlipToCone;
    for (const std::string &m : ir.messages) status.messages.push_back("Stage 4F: " + m);

    if (!ir.solved) {
        status.messages.push_back("Stopping: the field could not be integrated on Omega.");
        return false;
    }
    out.integratedMap = out.integration->getUV();

    // --- Stage 4R: the untangling (Sec. 7.2) ------------------------------
    std::vector<Point> psi0 = out.integratedMap;
    if (ir.flippedFaces > 0 && options.untangle) {
        bool untangled = false;
        if (options.localUntangle) {
            // The fallback above may already have climbed it on this solve.
            if (!ladderRan) ladder = climb(out.integratedMap, out.usedAxis, ir.flippedFaces);
            status.localUntangleRan = true;
            status.localUntangleFlippedFaces = ladder.leastLeft;
            status.reprojectionRan = ladder.reprojectionRan;
            status.reprojectionKept = ladder.reprojectionKept;
            status.reprojectionFlippedFaces = ladder.reprojectionFlips;
            if (!ladder.map.empty()) {
                psi0 = ladder.map;
                untangled = true;
                status.localUntangleLevel = ladder.level;
                std::ostringstream oss;
                oss << "Sec. 7.2a: the " << ir.flippedFaces
                    << " inverted face(s) were a local tangle, cleared at rung "
                    << status.localUntangleLevel << " of the ladder " << ladder.how
                    << ". psi_0 keeps the cone angles and the shape the integration gave it, "
                    << "and Sec. 7.2's Tutte pass was not needed.";
                status.messages.push_back("Stage 4R: " + oss.str());
            } else {
                std::ostringstream oss;
                oss << "Sec. 7.2a took the tangle from " << ir.flippedFaces << " face(s) to "
                    << status.localUntangleFlippedFaces << " at best";
                if (status.localUntangleFlippedFaces == 0) {
                    oss << ", and cleared it only on the rung that lets the seam go, where "
                        << "Sec. 6.5's projection could not put it back";
                }
                oss << ", so Sec. 7.2's Tutte pass runs on the whole map -- and psi_0 loses "
                    << "the cone angles, the boundary and the shape the integration gave it.";
                status.messages.push_back("Stage 4R: " + oss.str());
            }
        }

        if (!untangled) {
            std::vector<Point> repaired = untangle(in, options, status, out);
            if (!repaired.empty()) {
                psi0 = std::move(repaired);
                if (options.reprojectAfterUntangle && !out.usedAxis.empty()) {
                    int f = 0;
                    std::vector<Point> put = reproject(psi0, out.usedAxis, f);
                    status.reprojectionRan = true;
                    status.reprojectionFlippedFaces = f;
                    if (!put.empty()) {
                        status.reprojectionKept = true;
                        psi0 = std::move(put);
                        status.messages.push_back(
                            "Stage 4R: Sec. 6.5 put the seam and the alignment back on the "
                            "untangled map as equalities, and nothing inverted doing it.");
                    } else {
                        std::ostringstream oss;
                        oss << "Sec. 6.5's projection onto the alignment inverted "
                            << status.reprojectionFlippedFaces
                            << " face(s), so the untangled map stands as it is and Stage 6 has "
                            << "Q3 to reach rather than to hold. A Tutte map is a long way from "
                            << "the constraint set, which is the argument for making Sec. 7.2a "
                            << "succeed rather than for making the projection cleverer.";
                        status.messages.push_back("Stage 4R: " + oss.str());
                    }
                }
            }
        }
    } else if (ir.flippedFaces == 0) {
        status.messages.push_back(
            "Stage 4F: the integration inverted nothing, so Sec. 7.2's untangling was not "
            "needed. That happens on gently curved, well-aligned models and is not to be "
            "assumed.");
    }

    // --- Stage 4R, Sec. 7.2b: back onto the alignment -----------------------
    // Whether psi_0 holds the whole alignment as it stands: every held edge's
    // held coordinate difference at rounding, relative to the image. It does
    // after a clean solve or Sec. 7.2a at rung 0; it does not after a fallback
    // that dropped chains, the Tutte pass, or rung 1 with its projection
    // refused -- and that, not which of those paths was taken, is what decides
    // whether Stage 5 is handed an unaligned map.
    auto alignmentGap = [&](const std::vector<Point> &m, const std::vector<int> &axis) {
        const Mesh &cm = in.cut.getCutMesh();
        if (m.size() != cm.vertices.size() || axis.size() != cm.edges.size()) return 0.0;
        Point lo = m.front(), hi = m.front();
        for (const Point &p : m) {
            lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
            hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
        }
        const double extent = std::max(std::hypot(hi[0] - lo[0], hi[1] - lo[1]), 1e-300);
        double worst = 0.0;
        for (size_t e = 0; e < cm.edges.size(); ++e) {
            const int hold = axis[e];
            if (hold != 0 && hold != 1) continue;
            worst = std::max(worst, std::fabs(m[cm.edges[e][1]][hold] - m[cm.edges[e][0]][hold]));
        }
        return worst / extent;
    };
    status.alignmentHeld = fullAxis.empty() || fullAlignedMap.size() != psi0.size() ||
                           alignmentGap(psi0, fullAxis) <= 1e-8;
    if (options.pullOntoAlignment && !status.alignmentHeld) {
        std::vector<Point> pulled = pullOntoAlignment(in, options, status, psi0, fullAlignedMap, fullAxis);
        if (!pulled.empty()) {
            psi0 = std::move(pulled);
            if (status.pullProjected) out.usedAxis = fullAxis;
        }
        std::ostringstream oss;
        if (status.pullKept) {
            oss << "Sec. 7.2b: psi_0 had lost some of Sec. 6.4's alignment, and was pulled back "
                << "onto it under the barrier -- Q3 " << status.pullBoundaryResidual
                << ", features " << status.pullFeatureResidual << " -- "
                << (status.pullProjected
                        ? "and Sec. 6.5's projection then made it exact without inverting anything."
                        : "the projection onto the exact alignment inverting, so it stands at "
                          "that.")
                << " Stage 5 seeds Gamma_topo from this map, and an unaligned one seeds it "
                << "from separatrices the misaligned boundary sent astray.";
        } else {
            oss << "Sec. 7.2b could not pull psi_0 back onto the alignment, so Stage 5 seeds "
                << "from it as Stage 4R left it.";
        }
        status.messages.push_back("Stage 4R: " + oss.str());
    }

    out.psi0 = std::move(psi0);
    return true;
}

// ---------------------------------------------------------------------------
// runMapStages()  --  Stages 4C, 3F, 4F and 4R: from the cone set to psi_0
//
// Everything between the cut and the Immersion, on the cone set runConeFront()
// left. Leaves psi_0 in psi0Map; false at the stops no later stage survives.
// ---------------------------------------------------------------------------
bool TORSION::runMapStages() {
    // --- Sec. 4: the flat cone metric of the cone set ----------------------
    //
    // Built here, before the frames, because it is two things at once: the
    // conformal factor the frame's target Jacobian is missing (an unscaled
    // frame asks the map to be an isometry everywhere, which a map with cones
    // cannot be), and the reference metric E1 measures against. Both want the
    // same u and it is one solve.
    coneMetric.reset();
    if (options.reference == Options::Reference::Cone || options.conformalSizing) {
        ConeMetric::Options cmopts;
        try {
            coneMetric = std::make_unique<ConeMetric>(*mesh, *cones, cmopts);
        } catch (const std::exception &e) {
            status.messages.push_back(std::string("Stage 4C: could not build the cone metric: ") +
                                      e.what());
        }
    }
    if (coneMetric) {
        const ConeMetric::Report &cmr = coneMetric->getReport();
        status.coneMetricRan = true;
        status.coneMetricSolved = cmr.solved;
        status.coneMetricConverged = cmr.converged;
        status.coneMetricNewtonSteps = cmr.newtonIterations;
        status.coneMetricInitialError = cmr.initialError;
        status.coneMetricLinearError = cmr.linearError;
        status.coneMetricFinalError = cmr.finalError;
        status.coneMetricConeResidual = cmr.coneResidual;
        status.coneMetricMinScale = cmr.minScale;
        status.coneMetricMaxScale = cmr.maxScale;
        status.coneMetricNonRealisable = cmr.nonRealisable;
        for (const std::string &m : cmr.messages) status.messages.push_back("Stage 4C: " + m);
        if (!cmr.solved) coneMetric.reset();
    }

    // --- Stage 3F: the frames, the matchings and the audit (Sec. 5) -------
    FieldFrames::Options fopts;
    fopts.targetEdge = options.targetEdge;
    fopts.sizing = options.sizing;
    if (options.conformalSizing && coneMetric && options.sizing.empty()) {
        fopts.sizing = coneMetric->sizingField(options.targetEdge);
        // The image metric exactly, rather than the per-face average of |e|/h
        // that a sizing field can offer: this one is a vertex scaling and is
        // realisable on every face, which is what Immersion's Heron areas want.
        fopts.metricLengths = coneMetric->edgeLengths();
        if (options.targetEdge > 0.0) {
            for (double &l : fopts.metricLengths) l /= options.targetEdge;
        }
        status.conformalSizingUsed = true;
    }
    fopts.referenceIndex = fieldIndex;
    fopts.buildAlignment = options.alignInIntegration;
    fopts.alignAcrossSeams = options.alignAcrossSeams;
    fopts.reconcileSectors = options.reconcileSectors;
    try {
        frames = std::make_unique<FieldFrames>(*field, *cutter, *cones, fopts,
                                               interfaces.get());
    } catch (const std::exception &e) {
        status.messages.push_back(std::string("Stage 3F: could not comb the field: ") + e.what());
        return false;
    }
    const FieldFrames::Report &ffr = frames->getReport();
    status.framesRan = true;
    status.combedFaces = ffr.combedFaces;
    status.unreachedFaces = ffr.unreachedFaces;
    status.combingDefects = ffr.combingDefects;
    status.maxFrameJump = ffr.maxFrameJump;
    status.indexMismatches = ffr.indexMismatches;
    status.clusteredConePairs = ffr.clusteredConePairs;
    status.highIndexCones = ffr.highIndexCones;
    status.leftHandedFrames = ffr.leftHandedFrames;
    status.maxMetricDisagreement = ffr.maxMetricDisagreement;
    status.framesValid = ffr.valid;
    status.alignmentChains = ffr.alignmentChains;
    status.alignedBoundaryEdges = ffr.alignedBoundaryEdges;
    status.alignedInterfaceEdges = ffr.alignedInterfaceEdges;
    status.alignmentOverrides = ffr.alignmentOverrides;
    status.alignmentSeamCrossings = ffr.alignmentSeamCrossings;
    status.alignmentClosedChains = ffr.alignmentClosedChains;
    status.maxAlignmentResidual = ffr.maxAlignmentResidual;
    for (const std::string &m : ffr.messages) status.messages.push_back("Stage 3F: " + m);

    // Sec. 13's spot check. Two combs of a disk differ by the constant the two
    // seeds differ by and by nothing else, so the *differences* agree exactly;
    // an integer that varies is the BFS having taken a route the other did not
    // and got a different answer for it.
    if (options.checkSecondSeed && ffr.unreachedFaces == 0 && ffr.faces > 1) {
        FieldFrames::Options second = fopts;
        second.seedFace = ffr.faces / 2;
        // The spot check is about the branch and nothing else, and the
        // alignment is several passes over dS and the interface network.
        second.buildAlignment = false;
        try {
            FieldFrames other(*field, *cutter, *cones, second);
            const int shift = other.branch()[ffr.seedFace] - frames->branch()[ffr.seedFace];
            for (int t = 0; t < ffr.faces; ++t) {
                if (other.branch()[t] - frames->branch()[t] != shift) {
                    status.secondSeedAgrees = false;
                    break;
                }
            }
        } catch (const std::exception &) {
            status.secondSeedAgrees = false;
        }
        if (!status.secondSeedAgrees) {
            status.messages.push_back(
                "Stage 3F: combing from a second seed gave a different branch. On a disk with "
                "every cone on its boundary the BFS is path-independent, so this is the same "
                "leak across G that the loop check reports, seen from the other side.");
        }
    }

    if (ffr.unreachedFaces > 0 || ffr.combingDefects > 0) {
        status.messages.push_back(
            "Stopping before the integration: the branch of the field over Omega is not single "
            "valued, so the seam transitions it would be constrained by are not defined.");
        return false;
    }

    // --- Stage 4F: the integration (Sec. 6) -------------------------------
    // The scaffold exists for one reason: the constraint rows need the arcs of
    // G, their (e+, e-) pairing and their quarter turns, and all three are
    // Immersion::buildArcs()' work. It is handed Omega's own coordinates as a
    // placeholder map, which it never uses for anything read here -- the arcs
    // and the pairing come off ConeCut, and the k off the frames.
    try {
        scaffold = std::make_unique<Immersion>(*cutter, *cones,
                                               cutter->getCutMesh().vertices,
                                               frames->fieldEdgeLengths(),
                                               frames->combedAngle());
    } catch (const std::exception &e) {
        status.messages.push_back(std::string("Stage 4F: could not build the seam pairing: ") + e.what());
        return false;
    }
    status.seamArcs = scaffold->getReport().arcs;
    status.frameKConflicts = scaffold->getReport().frameKConflicts;

    Psi0 built;
    const bool built0 = buildPsi0(MapStage{*mesh, *cutter, *cones, *frames, *scaffold,
                                           interfaces.get(), coneMetric.get()},
                                  options, status, built);
    integration = std::move(built.integration);
    tutte = std::move(built.tutte);
    usedAxis = std::move(built.usedAxis);
    integratedMap = std::move(built.integratedMap);
    if (!built0) return false;
    std::vector<Point> psi0 = std::move(built.psi0);
    psi0Map = std::move(psi0);
    return true;
}

bool TORSION::run() {
    status = Status();

    runFieldFront();

    // --- Stages 1 to 4R, once or twice ------------------------------------
    //
    // Stage 1 cancels every +1/-1 pair it finds inside one material region,
    // MERIDIAN's rule: the pair is in no region's count, so no layout needs it,
    // and a layout without it is simpler. For Pipeline A that is the whole
    // story, because the flow is driven by the cone set and never looks at the
    // field again. This pipeline integrates the field, and a pair the field
    // put there is how the field turned through the curvature between them:
    // cancelled, the smoothest field left has to make the same turn with no
    // singularity to make it at, and on a strongly curved interface that is
    // not a field any map follows. On rt_mushroom, whose interface rolls up
    // through more than a full turn, cancelling its ten units takes the full
    // alignment from 2 inverted faces to 675 and the layout from 210 clean
    // quadrilaterals to a failure; on turbine_blade, 2 to 129. On cruciform,
    // icf and concrete the same cancellation costs nothing and the layout is
    // simpler for it (62 patches against 318 on cruciform).
    //
    // What separates the two is not the distance between the pair -- icf's are
    // forty edges apart and cancel harmlessly, rt_mushroom's include pairs one
    // edge apart that do not -- but whether the alignment survives. So the
    // cancellation is tried first, and if it cancelled anything and Stage 4F
    // could not keep the whole alignment with it, Stages 1 to 4R run again on
    // the same field keeping the pairs, and the second is kept only if it keeps
    // the alignment the first lost. Stage 0 is not repeated; the interface
    // network is restored to what it was before the first run's balance().
    const Status afterField = status;
    std::unique_ptr<Interfaces> pristine =
        interfaces ? std::make_unique<Interfaces>(*interfaces) : nullptr;
    const bool cancel = options.cancelInterfaceDipoles;
    bool ok = runConeFront(cancel) && runMapStages();

    // Lost, here, means what Stage 4R handed on does not hold the whole of Sec.
    // 6.4's alignment exactly -- a fallback that dropped chains, a rung above
    // the first whose projection was refused, the Tutte pass -- whether or not
    // Sec. 7.2b then pulled it most of the way back.
    const bool lostAlignment = !ok || !status.alignmentHeld;
    if (options.retryKeepingDipoles && cancel && status.coneDipoleUnits > 0 && lostAlignment) {
        struct Variant {
            std::unique_ptr<Interfaces> interfaces;
            std::unique_ptr<ConeSingularities> cones;
            std::unique_ptr<ConeCut> cutter;
            std::unique_ptr<FieldFrames> frames;
            std::unique_ptr<ConeMetric> coneMetric;
            std::unique_ptr<Immersion> scaffold;
            std::unique_ptr<FieldIntegration> integration;
            std::unique_ptr<TutteEmbedding> tutte;
            std::vector<int> fieldIndex, usedAxis;
            std::vector<Point> integratedMap, psi0;
            Status status;
            bool ok = false;
        };
        auto stash = [&](bool okNow) {
            Variant v;
            v.interfaces = std::move(interfaces);
            v.cones = std::move(cones);
            v.cutter = std::move(cutter);
            v.frames = std::move(frames);
            v.coneMetric = std::move(coneMetric);
            v.scaffold = std::move(scaffold);
            v.integration = std::move(integration);
            v.tutte = std::move(tutte);
            v.fieldIndex = std::move(fieldIndex);
            v.usedAxis = std::move(usedAxis);
            v.integratedMap = std::move(integratedMap);
            v.psi0 = std::move(psi0Map);
            v.status = status;
            v.ok = okNow;
            return v;
        };
        auto restore = [&](Variant &v) {
            interfaces = std::move(v.interfaces);
            cones = std::move(v.cones);
            cutter = std::move(v.cutter);
            frames = std::move(v.frames);
            coneMetric = std::move(v.coneMetric);
            scaffold = std::move(v.scaffold);
            integration = std::move(v.integration);
            tutte = std::move(v.tutte);
            fieldIndex = std::move(v.fieldIndex);
            usedAxis = std::move(v.usedAxis);
            integratedMap = std::move(v.integratedMap);
            psi0Map = std::move(v.psi0);
            status = v.status;
        };

        Variant first = stash(ok);
        const int cancelled = first.status.coneDipoleUnits;
        status = afterField;
        interfaces = pristine ? std::make_unique<Interfaces>(*pristine) : nullptr;
        const bool okKeep = runConeFront(false) && runMapStages();
        const bool keepSecond = okKeep && status.alignmentHeld;

        std::ostringstream oss;
        oss << "Stage 1 cancelled " << cancelled << " +1/-1 unit(s) inside single regions and "
            << "Stage 4F could not then keep the whole alignment, so Stages 1 to 4R were run "
            << "again keeping them: ";
        if (keepSecond) {
            oss << "with the pairs, the alignment holds, and that run is the one kept -- the "
                << "field turns through its curved interfaces at the pairs, and without them "
                << "no map it integrates to follows the curve.";
            status.dipolesKept = true;
            status.dipoleRetryRan = true;
            status.messages.push_back("Stage 1: " + oss.str());
        } else {
            restore(first);
            ok = first.ok;
            oss << "that lost the alignment too, so the cancelled cone set stands.";
            status.dipoleRetryRan = true;
            status.messages.push_back("Stage 1: " + oss.str());
        }
        if (keepSecond) ok = okKeep;
    }
    if (!ok) return false;

    // --- Stage 4: psi_0 ---------------------------------------------------
    try {
        immersion = std::make_unique<Immersion>(*cutter, *cones, psi0Map,
                                                frames->fieldEdgeLengths(),
                                                frames->combedAngle());
    } catch (const std::exception &e) {
        status.messages.push_back(std::string("Stage 4: could not accept psi_0: ") + e.what());
        return false;
    }
    const Immersion::Report &imr = immersion->getReport();
    status.immersionValid = imr.valid;
    status.immersionFlippedFaces = imr.flippedFaces;
    status.seamArcs = imr.arcs;
    status.maxSnapError = imr.maxSnapError;
    status.frameKConflicts = imr.frameKConflicts;
    status.maxMetricResidual = imr.maxMetricResidual;
    for (const std::string &m : imr.messages) status.messages.push_back("Stage 4: " + m);

    if (imr.flippedFaces > 0) {
        std::ostringstream oss;
        oss << "Stopping before the layout energies: psi_0 has " << imr.flippedFaces
            << " inverted face(s), so Q1 does not hold, and Sec. 3.3's barrier can only "
            << "preserve Q1 and never repair it. This is the cost of the substitution and "
            << "Sec. 7.3 is the list of what to try: raise the field's smoothness near the "
            << "cones, refine where the per-triangle fit residual is largest, lower h there.";
        status.messages.push_back(oss.str());
        return false;
    }
    if (!options.runLayout) return true;

    // --- Stages 5 to 8, at the Gamma_topo tolerance that works ------------
    //
    // MERIDIAN::run's retry, run further down the same ladder. Gamma_topo can
    // be added to and never taken back, so a tolerance that over-seeds is not
    // recoverable by the repair loop and a tighter one is -- that asymmetry is
    // why the retries only ever tighten. MERIDIAN stops at the second rung; this
    // pipeline goes on to a third and then to none at all, because its psi_0
    // is a field's and not a flow's. Where Sec. 7.2b has pulled an unaligned map
    // back onto the alignment, more of its separatrices pass within a given
    // tolerance of a cone than a Ricci map's would, and on tooth that is the
    // difference between 13 seeded paths -- three of them back to their own cone
    // -- and a layout 26 faces short of four-sided, against 187 clean
    // quadrilaterals at a tolerance of 0.01. Each rung runs only while the best
    // layout so far still leaves part of S without a grid, and the best is what
    // is kept.
    const Status statusBeforeLayout = status;
    if (!runLayoutStages(options.topoNearMiss, options.seedTopoConstraints)) {
        return status.layoutValid;
    }
    status.topoNearMissUsed = options.topoNearMiss;
    {
        struct Rung { double nearMiss; bool seed; };
        std::vector<Rung> rungs;
        auto addRung = [&](double nm, bool seed) {
            if (seed && !(nm > 0.0)) return;
            if (seed && nm == options.topoNearMiss) return;
            for (const Rung &r : rungs) if (r.seed == seed && r.nearMiss == nm) return;
            rungs.push_back({nm, seed});
        };
        if (options.seedTopoConstraints) {
            addRung(options.topoNearMissRetry, true);
            addRung(options.topoNearMissLastRetry, true);
            if (options.topoRetryUnseeded) addRung(0.0, false);
        }

        double bestUnmeshable = arrangement ? MERIDIAN::unmeshableFraction(*arrangement) : 1.0;
        int bestUnresolved = status.separatricesUnresolved;
        double bestNearMiss = options.topoNearMiss;
        bool bestSeeded = options.seedTopoConstraints;
        std::ostringstream retryMsg;
        retryMsg << std::fixed << std::setprecision(2);
        int tried = 0;
        for (const Rung &rung : rungs) {
            if (!arrangement || !(bestUnmeshable > 0.0)) break;
            if (tried == 0) {
                retryMsg << "Stage 5: the layout at a near-miss tolerance of " << std::defaultfloat
                         << options.topoNearMiss << std::fixed << " left "
                         << 100.0 * bestUnmeshable << "% of S in faces Stage 10 has no grid for";
            }
            ++tried;

            Status keptStatus = status;
            std::unique_ptr<SubdomainLabels> keptLabels = std::move(labels);
            std::unique_ptr<LayoutEnergy> keptLayout = std::move(layout);
            std::unique_ptr<Separatrices> keptSeparatrices = std::move(separatrices);
            std::unique_ptr<Arrangement> keptArrangement = std::move(arrangement);

            status = statusBeforeLayout;
            const bool reached = runLayoutStages(rung.nearMiss, rung.seed);
            const double unmeshable =
                (reached && arrangement) ? MERIDIAN::unmeshableFraction(*arrangement) : 1.0;
            const int unresolved = status.separatricesUnresolved;
            const bool better = reached && arrangement &&
                                (unmeshable < bestUnmeshable - 1e-12 ||
                                 (unmeshable <= bestUnmeshable + 1e-12 &&
                                  unresolved < bestUnresolved));
            if (rung.seed) {
                retryMsg << "; seeded again at " << std::defaultfloat << rung.nearMiss
                         << std::fixed << ", " << 100.0 * unmeshable << "%";
            } else {
                retryMsg << "; with no seeding at all, the repair loop alone, "
                         << 100.0 * unmeshable << "%";
            }
            if (better) {
                bestUnmeshable = unmeshable;
                bestUnresolved = unresolved;
                bestNearMiss = rung.nearMiss;
                bestSeeded = rung.seed;
            } else {
                status = std::move(keptStatus);
                labels = std::move(keptLabels);
                layout = std::move(keptLayout);
                separatrices = std::move(keptSeparatrices);
                arrangement = std::move(keptArrangement);
            }
        }
        if (tried > 0) {
            retryMsg << ". Kept ";
            if (!bestSeeded) retryMsg << "the unseeded layout";
            else retryMsg << "the one at " << std::defaultfloat << bestNearMiss;
            retryMsg << ". Gamma_topo can only be added to, never taken back, so a tolerance "
                        "that over-seeds is not recoverable by the repair loop and a tighter "
                        "one is; that asymmetry is the whole reason this retry runs in this "
                        "direction and not the other.";
            status.topoNearMissUsed = bestSeeded ? bestNearMiss : 0.0;
            status.topoNearMissRetried = true;
            status.messages.push_back(retryMsg.str());
        }
    }

    // --- Stage 9 ----------------------------------------------------------
    if (!options.runSplines) return status.layoutValid;
    SplineFit::Options sfopts;
    sfopts.segments = options.splineSegments;
    sfopts.samples = options.splineSamples;
    sfopts.fitBoundaryArcs = options.fitBoundaryArcs;
    sfopts.fitInterfaceArcs = options.fitInterfaceArcs;
    try {
        splines = std::make_unique<SplineFit>(*arrangement, sfopts);
    } catch (const std::exception &e) {
        status.messages.push_back(std::string("Stage 9: could not be built: ") + e.what());
        return status.layoutValid;
    }
    const SplineFit::Report &spr = splines->getReport();
    status.splinesRan = true;
    status.splinePatches = spr.patches;
    status.splineControlPoints = spr.controlPointsPerArc;
    status.splineMaxDeviation = spr.maxDeviation;
    status.splinesWatertight = spr.watertight;
    status.splinesValid = spr.valid;
    for (const std::string &m : spr.messages) status.messages.push_back("Stage 9: " + m);

    // --- Stage 10 ---------------------------------------------------------
    if (!options.runQuadMesh) return status.layoutValid;
    QuadMesh::Options qopts;
    qopts.targetEdgeLength = options.quadTargetEdge;
    qopts.minIntervals = options.quadMinIntervals;
    qopts.maxIntervals = options.quadMaxIntervals;
    qopts.collapseSpan = options.quadCollapseSpan;
    qopts.useSplines = options.quadUseSplines;
    qopts.featuresOnTracedArcs = options.quadFeaturesOnTracedArcs;
    qopts.smoothingPasses = options.quadSmoothingPasses;
    qopts.smoothingThreshold = options.quadSmoothingThreshold;
    // Each excised rim has to come out with an even number of edges or Stage 11
    // has nothing it can fill it with. See QuadMesh::fixLoopParity.
    const std::vector<std::vector<int>> rims =
        DiskTemplate::rimArcs(*arrangement, inclusions);
    for (const std::vector<int> &r : rims) if (!r.empty()) qopts.evenLoops.push_back(r);
    try {
        quads = std::make_unique<QuadMesh>(*splines, qopts);
    } catch (const std::exception &e) {
        status.messages.push_back(std::string("Stage 10: could not be built: ") + e.what());
        return status.layoutValid;
    }
    const QuadMesh::Report &qr = quads->getReport();
    status.quadMeshRan = true;
    status.meshVertices = qr.vertices;
    status.meshQuads = qr.quads;
    status.meshChords = qr.chords;
    status.meshUnmeshedPatches = qr.unmeshedPatches;
    status.meshMinScaledJacobian = qr.minScaledJacobian;
    status.meshConforming = qr.conforming;
    status.meshValid = qr.valid;
    status.meshOddLoops = qr.oddLoops;
    status.meshParityChordsMoved = qr.parityChordsMoved;
    status.meshOddLoopsLeft = qr.oddLoopsLeft;
    for (const std::string &m : qr.messages) status.messages.push_back("Stage 10: " + m);

    // --- Stage 11: the O-grid templates -----------------------------------
    if (inclusions.empty()) return status.layoutValid;
    DiskTemplate::Options topts;
    topts.coreSquareness = options.diskCoreSquareness;
    topts.ringDepth = options.diskRingDepth;
    topts.smoothingPasses = options.diskSmoothingPasses;
    diskFill = std::make_unique<DiskTemplate>(
        quads->vertices(), quads->quads(), quads->quadMaterials(),
        DiskTemplate::rimVertexLoops(*arrangement, *quads, rims), inclusions, topts);
    const DiskTemplate::Report &dr = diskFill->getReport();
    status.diskTemplatesRan = true;
    status.diskTemplatesFilled = dr.filled;
    status.diskTemplateBlocks = dr.blocks;
    status.diskTemplateQuads = dr.quads;
    status.diskTemplatesRefused = dr.refusedOdd + dr.refusedShort + dr.refusedOpen;
    status.diskTemplateMinScaledJacobian = dr.templateMinScaledJacobian;
    status.mergedVertices = dr.mergedVertices;
    status.mergedQuads = dr.mergedQuads;
    status.mergedMinScaledJacobian = dr.minScaledJacobian;
    status.diskTemplatesValid = dr.valid;
    for (const std::string &m : dr.messages) status.messages.push_back("Stage 11: " + m);

    return status.layoutValid;
}
