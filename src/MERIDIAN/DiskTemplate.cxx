#include "MERIDIAN/DiskTemplate.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace {

// The scaled Jacobian of a quadrilateral: the smallest of the four corner
// cross products of its unit side vectors. One for a square, zero for a
// degenerate corner, negative for a folded element. Same definition, for the
// same reason, as QuadMesh's.
double scaledJacobian(const Point &p0, const Point &p1, const Point &p2,
                      const Point &p3) {
    const Point p[4] = {p0, p1, p2, p3};
    double worst = std::numeric_limits<double>::infinity();
    for (int k = 0; k < 4; ++k) {
        const Point a = p[(k + 1) % 4] - p[k];
        const Point b = p[(k + 3) % 4] - p[k];
        const double la = normP(a), lb = normP(b);
        if (la <= 0.0 || lb <= 0.0) return 0.0;
        worst = std::min(worst, cross2(a, b) / (la * lb));
    }
    return worst;
}

double signedArea(const std::vector<Point> &poly) {
    double s = 0.0;
    const size_t n = poly.size();
    for (size_t i = 0; i < n; ++i) s += cross2(poly[i], poly[(i + 1) % n]);
    return 0.5 * s;
}

}  // namespace

// ---------------------------------------------------------------------------
// detect()
//
// A circular inclusion is a material component whose own boundary
// Mesh::computeMaterialCircles accepted as a circle *and* which does not touch
// dS. Both halves matter and the second is the one that is easy to leave out:
// a single-material disk model is one component with a perfect circle fit, and
// taking it out would leave nothing behind.
// ---------------------------------------------------------------------------
std::vector<DiskTemplate::Inclusion> DiskTemplate::detect(
        const Mesh &mesh, std::vector<std::string> &messages) {
    std::vector<Inclusion> out;
    const size_t nComp = mesh.materialComponents.size();
    if (nComp < 2) return out;   // one component is the whole model
    if (mesh.triangleComponent.size() != mesh.triangles.size()) return out;

    for (size_t c = 0; c < nComp; ++c) {
        const MaterialComponent &mc = mesh.materialComponents[c];
        if (!mc.circle.isCircle || mc.circle.radius <= 0.0) continue;

        // The rim: every edge with the component on exactly one side. An edge
        // with the component on one side and nothing on the other is dS, and
        // that disqualifies the inclusion outright.
        bool touchesBoundary = false;
        std::unordered_map<int, std::array<int, 2>> nbr;   // rim vertex -> its two rim neighbours
        int rimEdges = 0;
        for (size_t e = 0; e < mesh.edges.size(); ++e) {
            const int t0 = mesh.edgeTriangles[e][0];
            const int t1 = mesh.edgeTriangles[e][1];
            const bool in0 = t0 >= 0 && mesh.triangleComponent[t0] == static_cast<int>(c);
            const bool in1 = t1 >= 0 && mesh.triangleComponent[t1] == static_cast<int>(c);
            if (in0 == in1) continue;
            if (t0 < 0 || t1 < 0) { touchesBoundary = true; break; }
            ++rimEdges;
            const int u = mesh.edges[e][0], v = mesh.edges[e][1];
            auto ins = [&](int p, int qq) {
                auto it = nbr.find(p);
                if (it == nbr.end()) it = nbr.emplace(p, std::array<int, 2>{-1, -1}).first;
                if (it->second[0] < 0) it->second[0] = qq;
                else if (it->second[1] < 0) it->second[1] = qq;
                else it->second[1] = -2;   // three rim edges here: not a simple loop
            };
            ins(u, v);
            ins(v, u);
        }

        std::ostringstream oss;
        if (touchesBoundary) {
            oss << "The circular region at (" << mc.circle.center[0] << ", "
                << mc.circle.center[1] << ") reaches dS, so it is not an inclusion: "
                << "its rim is not a closed curve of the layout and excising it "
                << "would open the model rather than punch a hole in it. Laid out "
                << "like any other region.";
            messages.push_back(oss.str());
            continue;
        }
        if (rimEdges < 3) continue;

        // Chain the rim edges into one closed loop. Anything else -- a vertex
        // with one neighbour or three -- means the component is not a disk and
        // the fit was accidental.
        bool ok = true;
        for (const auto &kv : nbr) {
            if (kv.second[0] < 0 || kv.second[1] < 0) { ok = false; break; }
        }
        std::vector<int> loop;
        if (ok) {
            const int start = nbr.begin()->first;
            int prev = -1, cur = start;
            loop.push_back(cur);
            while (true) {
                const auto &n = nbr.at(cur);
                const int nxt = (n[0] == prev) ? n[1] : n[0];
                if (nxt == start) break;
                if (nxt < 0 || loop.size() > nbr.size()) { ok = false; break; }
                prev = cur;
                cur = nxt;
                loop.push_back(cur);
            }
            if (loop.size() != nbr.size()) ok = false;
        }
        if (!ok) {
            oss << "The circular region at (" << mc.circle.center[0] << ", "
                << mc.circle.center[1] << ") has a rim that is not one simple closed "
                << "loop, so the circle fit does not describe it. Left alone.";
            messages.push_back(oss.str());
            continue;
        }

        // Counter-clockwise about the fitted centre, so that the disk is on the
        // left of the loop and the template's elements come out positive.
        std::vector<Point> poly;
        poly.reserve(loop.size());
        for (int v : loop) poly.push_back(mesh.vertices[v] - mc.circle.center);
        if (signedArea(poly) < 0.0) std::reverse(loop.begin(), loop.end());

        Inclusion inc;
        inc.component = static_cast<int>(c);
        inc.matId = mc.matId;
        inc.circle = mc.circle;
        inc.triangles = mc.triangles;
        inc.rim = loop;
        inc.area = std::fabs(signedArea(poly));
        out.push_back(std::move(inc));
    }
    return out;
}

// ---------------------------------------------------------------------------
// excise()
// ---------------------------------------------------------------------------
std::shared_ptr<Mesh> DiskTemplate::excise(const Mesh &mesh,
                                           std::vector<Inclusion> &inclusions) {
    if (inclusions.empty()) return nullptr;

    std::vector<bool> drop(mesh.triangles.size(), false);
    for (const Inclusion &inc : inclusions) {
        for (int t : inc.triangles) {
            if (t >= 0 && t < static_cast<int>(drop.size())) drop[t] = true;
        }
    }

    std::vector<Triangle> tris;
    std::vector<int> mats;
    tris.reserve(mesh.triangles.size());
    mats.reserve(mesh.triangles.size());
    std::vector<int> renum(mesh.vertices.size(), -1);
    std::vector<Point> pts;
    for (size_t t = 0; t < mesh.triangles.size(); ++t) {
        if (drop[t]) continue;
        Triangle tri = mesh.triangles[t];
        for (int k = 0; k < 3; ++k) {
            int &v = tri[k];
            if (renum[v] < 0) {
                renum[v] = static_cast<int>(pts.size());
                pts.push_back(mesh.vertices[v]);
            }
            v = renum[v];
        }
        tris.push_back(tri);
        mats.push_back(t < mesh.triangleMatId.size() ? mesh.triangleMatId[t] : 1);
    }
    if (tris.empty()) return nullptr;

    for (Inclusion &inc : inclusions) {
        inc.rimExcised.clear();
        inc.rimExcised.reserve(inc.rim.size());
        for (int v : inc.rim) inc.rimExcised.push_back(renum[v]);
    }

    return std::make_shared<Mesh>(pts, tris, mats);
}

// ---------------------------------------------------------------------------
// rimArcs()  --  see the header for why this is geometric
// ---------------------------------------------------------------------------
std::vector<std::vector<int>> DiskTemplate::rimArcs(
        const Arrangement &arr, const std::vector<Inclusion> &inclusions) {
    std::vector<std::vector<int>> out(inclusions.size());
    if (inclusions.empty()) return out;

    const std::vector<Arrangement::Arc> &arcs = arr.getArcs();
    for (size_t a = 0; a < arcs.size(); ++a) {
        if (arcs[a].kind != Arrangement::ArcKind::Boundary || arcs[a].points.empty()) continue;
        for (size_t i = 0; i < inclusions.size(); ++i) {
            const CircleFit &cf = inclusions[i].circle;
            bool on = true;
            for (const Point &p : arcs[a].points) {
                if (std::fabs(normP(p - cf.center) - cf.radius) > 0.05 * cf.radius) {
                    on = false;
                    break;
                }
            }
            if (on) { out[i].push_back(static_cast<int>(a)); break; }
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// rimVertexLoops()
//
// The arcs of a rim are chained through the nodes they share, and each
// contributes the vertices Stage 10 placed along it. A rim with an arc Stage 10
// left unmeshed, or whose arcs do not close into one cycle, comes back empty
// and the fill declines it with a reason.
// ---------------------------------------------------------------------------
std::vector<std::vector<int>> DiskTemplate::rimVertexLoops(
        const Arrangement &arr, const QuadMesh &qm,
        const std::vector<std::vector<int>> &arcsPerRim) {
    std::vector<std::vector<int>> out(arcsPerRim.size());
    const std::vector<Arrangement::Arc> &arcs = arr.getArcs();
    const std::vector<std::vector<int>> &nodesOn = qm.arcVertices();

    for (size_t i = 0; i < arcsPerRim.size(); ++i) {
        const std::vector<int> &mine = arcsPerRim[i];
        if (mine.size() < 2) continue;

        bool complete = true;
        std::unordered_map<int, std::vector<int>> at;   // node -> its rim arcs
        for (int a : mine) {
            if (a >= static_cast<int>(nodesOn.size()) || nodesOn[a].size() < 2) {
                complete = false;
                break;
            }
            at[arcs[a].from].push_back(a);
            at[arcs[a].to].push_back(a);
        }
        if (!complete) continue;
        for (const auto &kv : at) if (kv.second.size() != 2) { complete = false; break; }
        if (!complete) continue;

        std::vector<int> loop;
        std::vector<char> seen(mine.size(), 0);
        std::unordered_map<int, int> slot;
        for (size_t k = 0; k < mine.size(); ++k) slot[mine[k]] = static_cast<int>(k);

        int arc = mine[0];
        int node = arcs[arc].from;
        for (size_t step = 0; step < mine.size(); ++step) {
            seen[slot[arc]] = 1;
            const std::vector<int> &nv = nodesOn[arc];
            const bool forward = arcs[arc].from == node;
            // Every vertex but the last: the next arc contributes it.
            if (forward) for (size_t k = 0; k + 1 < nv.size(); ++k) loop.push_back(nv[k]);
            else         for (size_t k = nv.size() - 1; k > 0; --k) loop.push_back(nv[k]);
            node = forward ? arcs[arc].to : arcs[arc].from;

            const std::vector<int> &here = at[node];
            const int next = here[0] == arc ? here[1] : here[0];
            if (step + 1 < mine.size() && seen[slot[next]]) { complete = false; break; }
            arc = next;
        }
        if (!complete || arc != mine[0]) continue;
        if (loop.size() < 3) continue;
        out[i] = std::move(loop);
    }
    return out;
}

// ---------------------------------------------------------------------------
// chooseSides()
//
// Which four rim nodes are the corners of the template. The core is an a x b
// grid and the ring has a, b, a, b elements along its four sides, so the split
// is one starting node and one integer a with b = N/2 - a; both are searched
// exhaustively, which is N (N/2 - 1) evaluations and costs nothing.
//
// Two things are wanted and they disagree. The four corners want to be a
// quarter turn apart about the centre, because that is what makes the ring even
// in depth; and a and b want to be equal, because the core is an a x b grid on
// a region that is nearly square, so a != b is an aspect ratio of b/a on every
// element of it. The rim's node spacing is whatever the chord assignment left
// -- it was minimising something else -- so both cannot usually be had, and the
// second is worth about as much as fifteen degrees of the first, which is what
// sets the weight below.
// ---------------------------------------------------------------------------
void DiskTemplate::chooseSides(const std::vector<int> &loop, const Point &centre,
                               int &start, int &a, int &b) const {
    const int N = static_cast<int>(loop.size());
    std::vector<double> theta(N);
    for (int k = 0; k < N; ++k) {
        const Point d = verts[loop[k]] - centre;
        theta[k] = std::atan2(d[1], d[0]);
    }

    double best = std::numeric_limits<double>::infinity();
    start = 0; a = N / 4; b = N / 2 - a;
    for (int s = 0; s < N; ++s) {
        for (int aa = 1; aa <= N / 2 - 1; ++aa) {
            const int bb = N / 2 - aa;
            const int idx[5] = {s, s + aa, s + aa + bb, s + 2 * aa + bb, s + N};
            double score = 0.0;
            for (int m = 0; m < 4; ++m) {
                double d = theta[idx[m + 1] % N] - theta[idx[m] % N];
                while (d <= 0.0) d += 2.0 * M_PI;
                while (d > 2.0 * M_PI) d -= 2.0 * M_PI;
                const double r = d - M_PI_2;
                score += r * r;
            }
            const double t = std::log(static_cast<double>(aa) / bb);
            score += 0.5 * t * t;
            if (score < best) { best = score; start = s; a = aa; b = bb; }
        }
    }
}

// ---------------------------------------------------------------------------
// fillOne()
// ---------------------------------------------------------------------------
bool DiskTemplate::fillOne(int inclusion, const std::vector<int> &given,
                           const CircleFit &cf, int matId) {
    // Counter-clockwise about the centre. The loop arrives in whatever
    // direction the rim's arcs happened to chain, and the template's elements
    // are wound off it directly, so getting this wrong turns every one of them
    // inside out rather than producing anything subtler.
    std::vector<int> loop = given;
    {
        std::vector<Point> poly;
        poly.reserve(loop.size());
        for (int v : loop) poly.push_back(verts[v] - cf.center);
        if (signedArea(poly) < 0.0) std::reverse(loop.begin(), loop.end());
    }
    const int N = static_cast<int>(loop.size());
    std::ostringstream oss;
    oss << std::setprecision(4);
    if (N < options.minRimEdges) {
        ++report.refusedShort;
        oss << "The inclusion at (" << cf.center[0] << ", " << cf.center[1]
            << ") has only " << N << " edge(s) on its rim, which is not enough for a "
            << "ring and a core; it is left as a hole. Mesh it at a smaller target "
            << "edge length.";
        report.messages.push_back(oss.str());
        return false;
    }
    if (N % 2 != 0) {
        ++report.refusedOdd;
        oss << "The inclusion at (" << cf.center[0] << ", " << cf.center[1]
            << ") has " << N << " edge(s) on its rim, an odd number, and no "
            << "quadrangulation of a disk has an odd boundary (4F = 2E_i + N). It is "
            << "left as a hole; QuadMesh::Options::evenLoops is what is supposed to "
            << "have made this even.";
        report.messages.push_back(oss.str());
        return false;
    }

    int start = 0, a = 0, b = 0;
    chooseSides(loop, cf.center, start, a, b);

    // How many rows the ring has and how far in the core corners sit.
    //
    // The rim spacing h is the only length scale the template has, and once the
    // rim node count fixes a and b there are exactly two numbers left to choose.
    // They are not independent and there is no rule of thumb that gets both: a
    // deeper ring makes its rows the right depth and shrinks the core into a
    // dense patch in the middle, a shallower one gives a generous core and rows
    // three times too deep at the middle of a side. So the two are searched
    // against what actually matters, which is how far each family of elements
    // is from being h across:
    //
    //     ring row at a corner        (1 - rho) R / c
    //     ring row at the middle      (1 - (1 - 0.293 w) rho) R / c
    //     core cell, the two ways     L(rho) / a  and  L(rho) / b
    //
    // with L the length of one side of the core boundary, which for the rounded
    // square of `coreSquareness` w is between the quarter arc and the chord.
    // The score is the sum of squared logs of the four against h, so a row
    // twice too deep costs what one half as deep does; the search is a few
    // hundred evaluations and runs once per inclusion.
    double perimeter = 0.0;
    for (int k = 0; k < N; ++k) {
        perimeter += normP(verts[loop[(k + 1) % N]] - verts[loop[k]]);
    }
    const double h = perimeter / N;
    const double R = cf.radius;
    const double w0 = std::min(1.0, std::max(0.0, options.coreSquareness));
    // One side of the core boundary, in radii of the core: the quarter arc
    // blended with the chord it subtends.
    const double sideLen = (1.0 - w0) * M_PI_2 + w0 * std::sqrt(2.0);
    // How far in the middle of a side sits, as a fraction of the core radius.
    const double midFrac = 1.0 - w0 * (1.0 - std::sqrt(0.5));

    int c = options.ringDepth;
    double rho = options.coreRadius;
    if (c <= 0 || rho <= 0.0) {
        double best = std::numeric_limits<double>::infinity();
        int bestC = 1;
        double bestRho = 0.65;
        const int cMax = std::max(1, std::min(4, std::min(a, b)));
        for (int cc = (c > 0 ? c : 1); cc <= (c > 0 ? c : cMax); ++cc) {
            for (int step = 25; step <= 90; ++step) {
                const double r = step / 100.0;
                const double dCorner = (1.0 - r) * R / (cc * h);
                const double dMid = (1.0 - midFrac * r) * R / (cc * h);
                const double L = sideLen * r * R;
                const double dA = L / (a * h);
                const double dB = L / (b * h);
                if (dCorner <= 0.0 || dMid <= 0.0) continue;
                double score = 0.0;
                for (double d : {dCorner, dMid, dA, dB}) {
                    const double t = std::log(d);
                    score += t * t;
                }
                if (score < best) { best = score; bestC = cc; bestRho = r; }
            }
        }
        if (c <= 0) c = bestC;
        if (rho <= 0.0) rho = bestRho;
    }
    c = std::max(1, c);
    rho = std::min(0.90, std::max(0.25, rho));

    // The four corners, and the core boundary they and the rim define.
    const int cornerIdx[4] = {start % N, (start + a) % N, (start + a + b) % N,
                              (start + 2 * a + b) % N};
    Point C[4];
    for (int m = 0; m < 4; ++m) {
        C[m] = cf.center + (verts[loop[cornerIdx[m]]] - cf.center) * rho;
    }

    // The core boundary, one node radially in from each rim node.
    //
    // Placed at *uniform* steps between one core corner and the next, not at
    // the rim's own node fractions. The rim's spacing is whatever the chord
    // assignment left, and one side of an inclusion routinely carries twice
    // the nodes per unit length of the side facing it; a core that inherits
    // that gets a Coons blend crowded from both ends and cells a third the
    // size of their neighbours' in the middle. Uniform here puts the
    // unevenness in the ring instead, where one row of elements absorbs it by
    // skewing slightly rather than by changing size.
    const double w = std::min(1.0, std::max(0.0, options.coreSquareness));
    std::vector<Point> coreBoundary(N);
    {
        const int len[4] = {a, b, a, b};
        double ang[5];
        double rad[5];
        for (int m = 0; m < 4; ++m) {
            const Point d = verts[loop[cornerIdx[m]]] - cf.center;
            ang[m] = std::atan2(d[1], d[0]);
            rad[m] = normP(d) * rho;
        }
        ang[4] = ang[0];
        rad[4] = rad[0];
        for (int m = 0; m < 4; ++m) {
            while (ang[m + 1] <= ang[m]) ang[m + 1] += 2.0 * M_PI;
        }

        int k = start;
        for (int m = 0; m < 4; ++m) {
            for (int i = 0; i < len[m]; ++i) {
                const int kk = (k + i) % N;
                const double t = static_cast<double>(i) / len[m];
                const double th = ang[m] + t * (ang[m + 1] - ang[m]);
                const double rr = rad[m] + t * (rad[m + 1] - rad[m]);
                const Point onCircle{cf.center[0] + rr * std::cos(th),
                                     cf.center[1] + rr * std::sin(th)};
                const Point onChord = C[m] + (C[(m + 1) % 4] - C[m]) * t;
                coreBoundary[kk] = onCircle * (1.0 - w) + onChord * w;
            }
            k = (k + len[m]) % N;
        }
    }

    // ---- vertices -------------------------------------------------------
    // ring[k][r], r = 0 on the rim and r = c on the core boundary.
    std::vector<std::vector<int>> ring(N, std::vector<int>(c + 1, -1));
    std::vector<int> fresh;
    for (int k = 0; k < N; ++k) {
        ring[k][0] = loop[k];
        for (int r = 1; r <= c; ++r) {
            const double f = static_cast<double>(r) / c;
            ring[k][r] = static_cast<int>(verts.size());
            verts.push_back(verts[loop[k]] * (1.0 - f) + coreBoundary[k] * f);
            fresh.push_back(ring[k][r]);
        }
    }

    // core[i][j]: its boundary is the inner ring, walked in the same order.
    std::vector<std::vector<int>> core(a + 1, std::vector<int>(b + 1, -1));
    for (int i = 0; i <= a; ++i) core[i][0] = ring[(start + i) % N][c];
    for (int j = 0; j <= b; ++j) core[a][j] = ring[(start + a + j) % N][c];
    for (int i = 0; i <= a; ++i) core[a - i][b] = ring[(start + a + b + i) % N][c];
    for (int j = 0; j <= b; ++j) core[0][b - j] = ring[(start + 2 * a + b + j) % N][c];

    // The interior by transfinite interpolation of the four sides -- a Coons
    // blend, the same one Stage 10 uses on a layout patch, on a boundary that
    // is already nearly a square.
    for (int i = 1; i < a; ++i) {
        const double u = static_cast<double>(i) / a;
        for (int j = 1; j < b; ++j) {
            const double v = static_cast<double>(j) / b;
            const Point s0 = verts[core[i][0]], s2 = verts[core[i][b]];
            const Point d0 = verts[core[0][j]], d1 = verts[core[a][j]];
            const Point p00 = verts[core[0][0]], p10 = verts[core[a][0]];
            const Point p01 = verts[core[0][b]], p11 = verts[core[a][b]];
            Point p = s0 * (1.0 - v) + s2 * v + d0 * (1.0 - u) + d1 * u -
                      (p00 * ((1.0 - u) * (1.0 - v)) + p10 * (u * (1.0 - v)) +
                       p01 * ((1.0 - u) * v) + p11 * (u * v));
            core[i][j] = static_cast<int>(verts.size());
            verts.push_back(p);
            fresh.push_back(core[i][j]);
        }
    }

    // ---- elements -------------------------------------------------------
    const int firstCell = static_cast<int>(cells.size());
    for (int k = 0; k < N; ++k) {
        const int k1 = (k + 1) % N;
        for (int r = 0; r < c; ++r) {
            cells.push_back({ring[k][r], ring[k1][r], ring[k1][r + 1], ring[k][r + 1]});
            cellMaterial.push_back(matId);
        }
    }
    for (int i = 0; i < a; ++i) {
        for (int j = 0; j < b; ++j) {
            cells.push_back({core[i][j], core[i + 1][j], core[i + 1][j + 1], core[i][j + 1]});
            cellMaterial.push_back(matId);
        }
    }

    // ---- blocks ---------------------------------------------------------
    // The core, then the four sides of the ring, each as a structured grid in
    // its own (s, t) frame with t running inward from the rim.
    {
        Block blk;
        blk.inclusion = inclusion;
        blk.core = true;
        blk.ns = a;
        blk.nt = b;
        blk.vert.resize(static_cast<size_t>(a + 1) * (b + 1));
        for (int j = 0; j <= b; ++j)
            for (int i = 0; i <= a; ++i) blk.vert[j * (a + 1) + i] = core[i][j];
        grids.push_back(std::move(blk));
    }
    {
        const int len[4] = {a, b, a, b};
        int k = start;
        for (int m = 0; m < 4; ++m) {
            Block blk;
            blk.inclusion = inclusion;
            blk.ns = len[m];
            blk.nt = c;
            blk.vert.resize(static_cast<size_t>(len[m] + 1) * (c + 1));
            for (int r = 0; r <= c; ++r)
                for (int i = 0; i <= len[m]; ++i)
                    blk.vert[r * (len[m] + 1) + i] = ring[(k + i) % N][r];
            grids.push_back(std::move(blk));
            k = (k + len[m]) % N;
        }
    }

    const int sweeps = smooth(fresh, firstCell, h);
    report.smoothingSweeps = std::max(report.smoothingSweeps, sweeps);
    ++report.filled;
    report.blocks += 5;
    report.quads += static_cast<int>(cells.size()) - firstCell;
    report.vertices += static_cast<int>(fresh.size());
    return true;
}

// ---------------------------------------------------------------------------
// smooth()
//
// Laplacian sweeps over the nodes the template placed, the rim held. A move is
// taken only when the worst scaled Jacobian *around that node* does not fall,
// which makes the whole pass monotone in the quantity the mesh is judged on:
// the transfinite initialisation is already valid, so the worst the smoother
// can do is decline to improve it.
//
// The weaker test -- take the move unless it turns an element over -- is what
// this had first and it is not good enough. Laplacian smoothing equalises
// element *size*, and around the four valence-3 corners of the core that is
// bought by squashing the corner elements: on bubbles it lifted the mean scaled
// Jacobian from 0.9653 to 0.9660 and dropped the worst from 0.727 to 0.552.
// Stage 10's Winslow pass learned the same thing about the same trade
// (QuadMesh::Options::smoothingThreshold) and the answer is the same one.
// ---------------------------------------------------------------------------
int DiskTemplate::smooth(const std::vector<int> &movable, int firstCell, double scale) {
    if (options.smoothingPasses <= 0 || movable.empty()) return 0;

    // Incidence over this inclusion's own elements only: a movable node is
    // interior to the template, so every element using it was just added.
    std::unordered_map<int, std::vector<int>> inc;
    for (int q = firstCell; q < static_cast<int>(cells.size()); ++q) {
        for (int k = 0; k < 4; ++k) inc[cells[q][k]].push_back(q);
    }

    const double tol = options.smoothingTolerance * scale;
    int sweep = 0;
    for (; sweep < options.smoothingPasses; ++sweep) {
        double moved = 0.0;
        for (int v : movable) {
            auto it = inc.find(v);
            if (it == inc.end() || it->second.empty()) continue;

            Point sum{0.0, 0.0};
            int n = 0;
            for (int q : it->second) {
                for (int k = 0; k < 4; ++k) {
                    if (cells[q][k] != v) continue;
                    sum = sum + verts[cells[q][(k + 1) % 4]];
                    sum = sum + verts[cells[q][(k + 3) % 4]];
                    n += 2;
                }
            }
            if (n == 0) continue;

            const Point old = verts[v];
            double before = std::numeric_limits<double>::infinity();
            for (int q : it->second) {
                before = std::min(before, scaledJacobian(verts[cells[q][0]], verts[cells[q][1]],
                                                         verts[cells[q][2]], verts[cells[q][3]]));
            }
            verts[v] = old + (sum / static_cast<double>(n) - old) * 0.6;
            double after = std::numeric_limits<double>::infinity();
            for (int q : it->second) {
                after = std::min(after, scaledJacobian(verts[cells[q][0]], verts[cells[q][1]],
                                                       verts[cells[q][2]], verts[cells[q][3]]));
            }
            if (after < before - 1e-12) { verts[v] = old; continue; }
            moved = std::max(moved, normP(verts[v] - old));
        }
        if (moved < tol) break;
    }
    return sweep + 1;
}

// ---------------------------------------------------------------------------
// check()
// ---------------------------------------------------------------------------
void DiskTemplate::check() {
    report.mergedVertices = static_cast<int>(verts.size());
    report.mergedQuads = static_cast<int>(cells.size());

    std::unordered_map<MeshEdgeKey, int, MeshEdgeKeyHash> use;
    double minSJ = std::numeric_limits<double>::infinity();
    double sumSJ = 0.0;
    double minTemplate = std::numeric_limits<double>::infinity();
    double minE = std::numeric_limits<double>::infinity(), maxE = 0.0;
    const int firstTemplate = report.mergedQuads - report.quads;
    for (int q = 0; q < report.mergedQuads; ++q) {
        const auto &f = cells[q];
        const double sj = scaledJacobian(verts[f[0]], verts[f[1]], verts[f[2]], verts[f[3]]);
        minSJ = std::min(minSJ, sj);
        sumSJ += sj;
        if (q >= firstTemplate) minTemplate = std::min(minTemplate, sj);
        double area = 0.0;
        for (int k = 0; k < 4; ++k) {
            area += cross2(verts[f[k]], verts[f[(k + 1) % 4]]);
            const double L = normP(verts[f[(k + 1) % 4]] - verts[f[k]]);
            minE = std::min(minE, L);
            maxE = std::max(maxE, L);
            ++use[MeshEdgeKey(f[k], f[(k + 1) % 4])];
        }
        if (area <= 0.0) ++report.invertedQuads;
    }
    report.minScaledJacobian = std::isfinite(minSJ) ? minSJ : 0.0;
    report.meanScaledJacobian = report.mergedQuads ? sumSJ / report.mergedQuads : 0.0;
    report.templateMinScaledJacobian = std::isfinite(minTemplate) ? minTemplate : 0.0;
    report.minEdge = std::isfinite(minE) ? minE : 0.0;
    report.maxEdge = maxE;

    for (const auto &kv : use) {
        if (kv.second == 1) ++report.boundaryEdges;
        else if (kv.second == 2) ++report.interiorEdges;
        else ++report.nonManifoldEdges;
    }

    // Cracks: two distinct vertices at the same point, which is what a rim
    // meshed twice rather than shared would leave. Hashed on a grid at the
    // tolerance so this stays linear.
    double diag = 0.0;
    {
        Point lo{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
        Point hi{-lo[0], -lo[1]};
        for (const Point &p : verts) {
            lo[0] = std::min(lo[0], p[0]); lo[1] = std::min(lo[1], p[1]);
            hi[0] = std::max(hi[0], p[0]); hi[1] = std::max(hi[1], p[1]);
        }
        diag = normP(hi - lo);
    }
    const double eps = std::max(1e-12, 1e-9 * diag);
    std::unordered_map<long long, std::vector<int>> grid;
    auto key = [&](const Point &p, int dx, int dy) {
        const long long ix = static_cast<long long>(std::floor(p[0] / eps)) + dx;
        const long long iy = static_cast<long long>(std::floor(p[1] / eps)) + dy;
        return ix * 73856093LL ^ iy * 19349663LL;
    };
    for (int v = 0; v < static_cast<int>(verts.size()); ++v) {
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                auto it = grid.find(key(verts[v], dx, dy));
                if (it == grid.end()) continue;
                for (int u : it->second) {
                    if (normP(verts[v] - verts[u]) <= eps) ++report.cracks;
                }
            }
        }
        grid[key(verts[v], 0, 0)].push_back(v);
    }

    report.conforming = report.nonManifoldEdges == 0 && report.cracks == 0;
    report.valid = report.conforming && report.invertedQuads == 0 &&
                   report.refusedOdd == 0 && report.refusedShort == 0 &&
                   report.refusedOpen == 0 && report.filled == report.inclusions;
}

// ---------------------------------------------------------------------------
DiskTemplate::DiskTemplate(const std::vector<Point> &v,
                           const std::vector<std::array<int, 4>> &q,
                           const std::vector<int> &mat,
                           const std::vector<std::vector<int>> &rimLoops,
                           const std::vector<Inclusion> &inclusions,
                           const Options &opts)
    : options(opts), verts(v), cells(q), cellMaterial(mat) {
    if (cellMaterial.size() != cells.size()) cellMaterial.assign(cells.size(), 1);
    report.inclusions = static_cast<int>(inclusions.size());

    static const std::vector<int> kNoLoop;
    for (size_t i = 0; i < inclusions.size(); ++i) {
        const std::vector<int> &loop = i < rimLoops.size() ? rimLoops[i] : kNoLoop;
        if (loop.empty()) {
            ++report.refusedOpen;
            std::ostringstream oss;
            oss << std::setprecision(4)
                << "The inclusion at (" << inclusions[i].circle.center[0] << ", "
                << inclusions[i].circle.center[1] << ") has no closed rim in the "
                << "quadrilateral mesh -- Stage 10 left one of the patches beside it "
                << "unmeshed -- so there is nothing for the template to attach to. "
                << "Left as a hole.";
            report.messages.push_back(oss.str());
            continue;
        }
        fillOne(static_cast<int>(i), loop, inclusions[i].circle, inclusions[i].matId);
    }

    check();

    if (report.filled > 0) {
        std::ostringstream oss;
        oss << "Templated " << report.filled << " circular inclusion(s) as O-grids: "
            << report.blocks << " block(s), " << report.quads << " element(s) on "
            << report.vertices << " new vertex/vertices, worst scaled Jacobian "
            << std::setprecision(4) << report.templateMinScaledJacobian
            << " over the template elements.";
        report.messages.push_back(oss.str());
    }
}

// ---------------------------------------------------------------------------
bool DiskTemplate::writeOBJ(const std::string &filename) const {
    std::ofstream out(filename);
    if (!out) return false;
    out << "# CrossGen -- matrix layout mesh with templated circular inclusions\n";
    out << std::setprecision(17);
    for (const Point &p : verts) out << "v " << p[0] << " " << p[1] << " 0\n";
    int current = std::numeric_limits<int>::min();
    for (size_t q = 0; q < cells.size(); ++q) {
        const int m = q < cellMaterial.size() ? cellMaterial[q] : 1;
        if (m != current) { out << "usemtl mat" << m << "\n"; current = m; }
        out << "f " << cells[q][0] + 1 << " " << cells[q][1] + 1 << " "
            << cells[q][2] + 1 << " " << cells[q][3] + 1 << "\n";
    }
    return true;
}

bool DiskTemplate::writeVTU(const std::string &filename) const {
    std::ofstream out(filename);
    if (!out) return false;
    out << "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" "
        << "byte_order=\"LittleEndian\">\n  <UnstructuredGrid>\n    <Piece NumberOfPoints=\""
        << verts.size() << "\" NumberOfCells=\"" << cells.size() << "\">\n";
    out << "      <Points>\n        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
        << "format=\"ascii\">\n";
    out << std::setprecision(17);
    for (const Point &p : verts) out << p[0] << " " << p[1] << " 0\n";
    out << "        </DataArray>\n      </Points>\n      <Cells>\n";
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto &c : cells) out << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << "\n";
    out << "        </DataArray>\n        <DataArray type=\"Int32\" Name=\"offsets\" "
        << "format=\"ascii\">\n";
    for (size_t i = 1; i <= cells.size(); ++i) out << i * 4 << "\n";
    out << "        </DataArray>\n        <DataArray type=\"UInt8\" Name=\"types\" "
        << "format=\"ascii\">\n";
    for (size_t i = 0; i < cells.size(); ++i) out << "9\n";
    out << "        </DataArray>\n      </Cells>\n      <CellData Scalars=\"material\">\n";
    out << "        <DataArray type=\"Int32\" Name=\"material\" format=\"ascii\">\n";
    for (size_t q = 0; q < cells.size(); ++q)
        out << (q < cellMaterial.size() ? cellMaterial[q] : 1) << "\n";
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Float64\" Name=\"scaledJacobian\" format=\"ascii\">\n";
    for (const auto &c : cells)
        out << scaledJacobian(verts[c[0]], verts[c[1]], verts[c[2]], verts[c[3]]) << "\n";
    out << "        </DataArray>\n      </CellData>\n    </Piece>\n  </UnstructuredGrid>\n"
        << "</VTKFile>\n";
    return true;
}
