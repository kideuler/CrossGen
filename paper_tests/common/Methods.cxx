#include "Methods.hxx"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>

#include "SIPG/SIPG.hxx"
#include "crossfield/CrossField.hxx"
#include "polyvector/PolyVectors.hxx"

namespace paper {
namespace {

using Clock = std::chrono::steady_clock;
double since(Clock::time_point t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

double angleDiff(const std::complex<double> &from, const std::complex<double> &to) {
    const std::complex<double> d = to * std::conj(from);
    return std::atan2(d.imag(), d.real());
}

// The linear interpolant of the P1 field inside a triangle, at barycentric l.
std::complex<double> interp(const Mesh &m, int t, const Eigen::VectorXcd &uv,
                            const std::array<double, 3> &l) {
    const Triangle &tri = m.triangles[t];
    return l[0] * uv[tri[0]] + l[1] * uv[tri[1]] + l[2] * uv[tri[2]];
}

} // namespace

// ---------------------------------------------------------------------------
std::vector<char> interfaceEdgeFlags(const Mesh &m) {
    std::vector<char> flag(m.edges.size(), 0);
    if (m.triangleMatId.size() != m.triangles.size()) return flag;
    for (std::size_t e = 0; e < m.edges.size(); ++e) {
        if (m.isBoundaryEdge[e]) continue;
        const int a = m.edgeTriangles[e][0], b = m.edgeTriangles[e][1];
        if (a < 0 || b < 0) continue;
        if (m.triangleMatId[a] != m.triangleMatId[b]) flag[e] = 1;
    }
    return flag;
}

std::vector<int> interfaceEdges(const Mesh &m) {
    const std::vector<char> flag = interfaceEdgeFlags(m);
    std::vector<int> out;
    for (std::size_t e = 0; e < flag.size(); ++e) if (flag[e]) out.push_back(static_cast<int>(e));
    return out;
}

bool isMultiMaterial(const Mesh &m) {
    if (m.triangleMatId.empty()) return false;
    const int first = m.triangleMatId.front();
    for (int id : m.triangleMatId) if (id != first) return true;
    return false;
}

namespace {

// The mean edge length, which is what the tau floor is measured in.
double meanEdgeLength(const Mesh &m) {
    double s = 0.0;
    int n = 0;
    for (std::size_t e = 0; e < m.edges.size(); ++e) {
        s += normP(m.vertices[m.edges[e][1]] - m.vertices[m.edges[e][0]]);
        ++n;
    }
    return n ? s / n : 0.0;
}

// The bounding-box diagonal, which is what tau is measured in.
double boundingDiagonal(const Mesh &m) {
    double mnx = 1e300, mxx = -1e300, mny = 1e300, mxy = -1e300;
    for (const Point &p : m.vertices) {
        mnx = std::min(mnx, p[0]); mxx = std::max(mxx, p[0]);
        mny = std::min(mny, p[1]); mxy = std::max(mxy, p[1]);
    }
    return std::hypot(mxx - mnx, mxy - mny);
}

// The ladder of tau scales the continuation walks, largest first. The first
// entry is always 1 -- the shipped tau -- so the continuation starts from the
// field the single-tau scheme would have returned and only ever refines it.
std::vector<double> tauLadder(const Mesh &m, const MethodOptions &o, bool enabled) {
    std::vector<double> ladder{o.tauScale};
    if (!enabled) return ladder;

    const double D = boundingDiagonal(m);
    const double h = meanEdgeLength(m);
    if (D <= 0.0 || h <= 0.0 || o.tauRatio <= 0.0 || o.tauRatio >= 1.0) return ladder;

    // ell = sqrt(8 gamma tau) <= tauFloorEdges * h  =>  tau <= (c h)^2 / (8 gamma)
    const double tau0 = o.tauScale * D * D / 10.0;
    const double tauMin = (o.tauFloorEdges * h) * (o.tauFloorEdges * h) / (8.0 * o.gamma);

    double tau = tau0;
    while (tau * o.tauRatio > tauMin && ladder.size() < 40) {
        tau *= o.tauRatio;
        ladder.push_back(ladder.back() * o.tauRatio);
    }
    return ladder;
}

} // namespace

// ---------------------------------------------------------------------------
FieldRun runSIPG(const std::shared_ptr<Mesh> &m, const MethodOptions &o) {
    FieldRun r;
    r.method = "SIPG";
    const Clock::time_point tAll = Clock::now();
    try {
        const metrics::EdgeWeights w = o.recordHistory ? metrics::edgeWeights(*m, o.gamma)
                                                       : metrics::EdgeWeights{};
        const std::vector<char> skip = o.recordHistory ? interfaceEdgeFlags(*m) : std::vector<char>{};
        const double nT = static_cast<double>(m->triangles.size());
        const std::vector<double> ladder = tauLadder(*m, o, o.tauContinuation);

        Eigen::VectorXcd carried;   // the previous level's field
        for (std::size_t level = 0; level < ladder.size(); ++level) {
            SIPG solver(m, o.maxSteps, o.gamma);
            solver.setPinDiskCenters(o.pinDiskCenters);
            solver.setTauScale(ladder[level]);
            solver.setHardBoundaryConditions(o.hardBoundary);
            if (o.alignInterfaces && isMultiMaterial(*m))
                solver.setAlignedInteriorEdges(interfaceEdges(*m));

            const Clock::time_point tA = Clock::now();
            solver.initialize();
            r.assembleSeconds += since(tA);

            // Every level after the first starts from the one before it. The
            // Dirichlet triangles keep the data this level assembled -- the
            // constraint does not travel with the field, it is re-imposed.
            if (carried.size() == solver.u_k_prev.size()) {
                Eigen::VectorXcd start = carried;
                for (const auto &[ti, bc] : solver.getBoundaryData()) start[ti] = bc;
                solver.u_k_prev = start;
                solver.u_k = start;
            }

            // Only the last level takes the convergence history: the earlier
            // ones are a continuation path, not the answer.
            const bool record = o.recordHistory && level + 1 == ladder.size();
            const int cap = (ladder.size() == 1) ? o.maxSteps
                                                 : std::min(o.maxSteps, o.tauLevelSteps);

            const Clock::time_point tS = Clock::now();
            bool converged = false;
            for (int i = 0; i < cap; ++i) {
                solver.step();
                ++r.iterations;
                if (record) {
                    r.energyHistory.push_back(metrics::commonEnergy(*m, w, solver.u_k_prev, skip));
                    r.sipgEnergyHistory.push_back(metrics::sipgEnergy(*m, o.gamma, solver.u_k_prev, skip));
                    r.incrementHistory.push_back(solver.error);
                }
                if (solver.error < 2.0 * nT * o.convergenceTol) {
                    converged = true;
                    if (!o.forceSteps) break;
                }
            }
            r.solveSeconds += since(tS);
            r.converged = converged;
            carried = solver.u_k_prev;
        }
        r.u = carried;
        r.ok = true;
    } catch (const std::exception &e) {
        r.error = e.what();
    }
    r.totalSeconds = since(tAll);
    return r;
}

// ---------------------------------------------------------------------------
FieldRun runP1MBO(const std::shared_ptr<Mesh> &m, const MethodOptions &o) {
    FieldRun r;
    r.method = "B1";
    const Clock::time_point tAll = Clock::now();
    try {
        // B1 is run as published by default: one tau = D^2/10, its own seeded
        // random start. `b1TauContinuation` gives it the same ladder our method
        // uses, which is what separates "the continuation" from "the
        // discretisation" as the source of any gain -- see MethodOptions.
        const std::vector<double> ladder = tauLadder(*m, o, o.b1TauContinuation);
        const double nV = static_cast<double>(m->vertices.size());
        double assemble = 0.0, solve = 0.0;
        int iters = 0;
        bool conv = false;
        Eigen::VectorXcd carried;
        std::vector<std::pair<int, double>> singular;   // from the last level

        for (std::size_t level = 0; level < ladder.size(); ++level) {
            CrossField cf(m, o.maxSteps);
            cf.setTauScale(ladder[level]);
            const Clock::time_point tA = Clock::now();
            cf.initialize(1, o.seed ? o.seed : 1u);   // method 1: seeded random interior
            assemble += since(tA);

            if (carried.size() == cf.u_k_prev.size()) {
                // Keep this level's Dirichlet data at the boundary vertices and
                // take the previous level's field everywhere else.
                Eigen::VectorXcd start = carried;
                for (int v : m->boundaryVertices) start[v] = cf.u_k_prev[v];
                cf.u_k_prev = start;
                cf.u_k = start;
            }

            const int cap = (ladder.size() == 1) ? o.maxSteps
                                                 : std::min(o.maxSteps, o.tauLevelSteps);
            const Clock::time_point tS = Clock::now();
            conv = false;
            for (int i = 0; i < cap; ++i) {
                cf.step();
                ++iters;
                if (cf.error < 2.0 * nV * o.convergenceTol) { conv = true; break; }
            }
            solve += since(tS);
            carried = cf.u_k_prev;
            if (level + 1 == ladder.size()) {
                cf.computeSingularities();
                singular = cf.singularTriangles;
            }
        }

        // The conversion is a separate, documented step and returns a run of
        // its own; everything the solve knows is copied onto it here.
        r = convertP1ToFaces(*m, carried, singular);
        r.method = "B1";
        r.converged = conv;
        r.iterations = iters;
        r.assembleSeconds = assemble;
        r.solveSeconds = solve;
        r.ok = true;
    } catch (const std::exception &e) {
        r.ok = false;
        r.error = e.what();
    }
    r.totalSeconds = since(tAll);
    return r;
}

// ---------------------------------------------------------------------------
FieldRun runPolyVector(const std::shared_ptr<Mesh> &m, const MethodOptions &o) {
    (void)o;
    FieldRun r;
    r.method = "B2";
    const Clock::time_point tAll = Clock::now();
    try {
        PolyField pf(m);
        const Clock::time_point tS = Clock::now();
        pf.solveForPolyCoeffs();
        pf.convertToFieldVectors();
        r.solveSeconds = since(tS);
        r.iterations = 1;
        r.converged = true;

        const int nT = static_cast<int>(m->triangles.size());
        r.u.resize(nT);
        for (int t = 0; t < nT; ++t) {
            const Point &d = pf.field[t].u;
            const double th = std::atan2(d[1], d[0]);
            r.u[t] = std::exp(std::complex<double>(0.0, 4.0 * th));
            if (!std::isfinite(r.u[t].real()) || !std::isfinite(r.u[t].imag()))
                r.u[t] = std::complex<double>(1.0, 0.0);
        }
        r.ok = true;
    } catch (const std::exception &e) {
        r.error = e.what();
    }
    r.totalSeconds = since(tAll);
    return r;
}

// ---------------------------------------------------------------------------
FieldRun convertP1ToFaces(const Mesh &m, const Eigen::VectorXcd &uVertex,
                          const std::vector<std::pair<int, double>> &p1SingularFaces) {
    FieldRun r;
    r.method = "B1";
    r.hasConversion = true;
    r.p1Vertex = uVertex;

    const int nT = static_cast<int>(m.triangles.size());
    r.u.resize(nT);

    // --- the face values ---------------------------------------------------
    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = m.triangles[t];
        const std::complex<double> avg = (uVertex[tri[0]] + uVertex[tri[1]] + uVertex[tri[2]]) / 3.0;
        const double mag = std::abs(avg);
        if (mag > 1e-9) {
            r.u[t] = avg / mag;
        } else {
            // The interpolant vanishes inside this triangle -- it is a singular
            // face of the P1 field -- so the average carries no direction at all
            // and the best-effort choice is one corner's value. There is no
            // right answer here, and that is precisely the outline's point.
            ++r.conversion.degenerateFaces;
            r.u[t] = uVertex[tri[0]];
        }
    }

    // --- the holonomy check, per interior edge ------------------------------
    const int kSamples = 32;
    for (std::size_t e = 0; e < m.edges.size(); ++e) {
        if (m.isBoundaryEdge[e]) continue;
        const int ti = m.edgeTriangles[e][0], tj = m.edgeTriangles[e][1];
        if (ti < 0 || tj < 0) continue;

        // The barycentric coordinates of the shared edge's midpoint in each
        // triangle: 1/2 on the two endpoints, 0 on the third corner.
        auto midBary = [&](int t) {
            std::array<double, 3> l{0.0, 0.0, 0.0};
            const Triangle &tri = m.triangles[t];
            for (int k = 0; k < 3; ++k)
                if (tri[k] == m.edges[e][0] || tri[k] == m.edges[e][1]) l[k] = 0.5;
            return l;
        };
        const std::array<double, 3> cen{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};

        double transported = 0.0;
        bool trackable = true;
        auto walk = [&](int t, const std::array<double, 3> &from, const std::array<double, 3> &to) {
            std::complex<double> prev = interp(m, t, uVertex, from);
            if (std::abs(prev) < 1e-3) trackable = false;
            for (int s = 1; s <= kSamples && trackable; ++s) {
                const double a = static_cast<double>(s) / kSamples;
                std::array<double, 3> l{};
                for (int k = 0; k < 3; ++k) l[k] = from[k] * (1.0 - a) + to[k] * a;
                const std::complex<double> cur = interp(m, t, uVertex, l);
                if (std::abs(cur) < 1e-3) { trackable = false; break; }
                transported += angleDiff(prev, cur);
                prev = cur;
            }
        };
        walk(ti, cen, midBary(ti));
        walk(tj, midBary(tj), cen);

        if (!trackable) { ++r.conversion.untrackableEdges; continue; }

        const double resampled = angleDiff(r.u[ti], r.u[tj]);
        const int lost = static_cast<int>(std::lround((transported - resampled) / (2.0 * M_PI)));
        if (lost != 0) ++r.conversion.holonomyLostEdges;
    }

    // --- the topological content, before and after --------------------------
    struct Sing { int vertex; int index4; };
    std::vector<Sing> before;
    for (const auto &[t, idx] : p1SingularFaces) {
        if (t < 0 || t >= nT) continue;
        const Triangle &tri = m.triangles[t];
        const Point c = (m.vertices[tri[0]] + m.vertices[tri[1]] + m.vertices[tri[2]]) / 3.0;
        int best = tri[0];
        double bestD = std::numeric_limits<double>::max();
        for (int k = 0; k < 3; ++k) {
            const double d = normP(m.vertices[tri[k]] - c);
            if (d < bestD) { bestD = d; best = tri[k]; }
        }
        before.push_back({best, static_cast<int>(std::lround(idx * 4.0))});
    }

    const metrics::SingularityReport sr = metrics::singularities(m, r.u);
    std::vector<Sing> after;
    for (const metrics::Singularity &s : sr.interior) after.push_back({s.vertex, s.index4});

    r.conversion.singularitiesBefore = static_cast<int>(before.size());
    r.conversion.singularitiesAfter = static_cast<int>(after.size());
    for (const Sing &s : before) r.conversion.indexSum4Before += s.index4;
    for (const Sing &s : after) r.conversion.indexSum4After += s.index4;

    // Greedy nearest matching at equal index within three local edge lengths.
    // A singularity that moved a little is the same singularity; one that has no
    // partner at all was created or destroyed by the resampling.
    const std::vector<double> hLocal = metrics::localEdgeLength(m);
    std::vector<char> takenAfter(after.size(), 0);
    for (const Sing &b : before) {
        int bestJ = -1;
        double bestD = std::numeric_limits<double>::max();
        for (std::size_t j = 0; j < after.size(); ++j) {
            if (takenAfter[j] || after[j].index4 != b.index4) continue;
            const double d = normP(m.vertices[after[j].vertex] - m.vertices[b.vertex]);
            if (d < bestD) { bestD = d; bestJ = static_cast<int>(j); }
        }
        const double reach = 3.0 * std::max(1e-30, hLocal[b.vertex]);
        if (bestJ >= 0 && bestD <= reach) takenAfter[bestJ] = 1;
        else ++r.conversion.annihilated;
    }
    for (char t : takenAfter) if (!t) ++r.conversion.created;
    r.conversion.topologyChanged = r.conversion.created + r.conversion.annihilated;

    r.ok = true;
    return r;
}

// ---------------------------------------------------------------------------
std::vector<std::string> methodNames() { return {"SIPG", "B1", "B2"}; }

FieldRun runMethod(const std::string &name, const std::shared_ptr<Mesh> &m, const MethodOptions &o) {
    if (name == "SIPG") return runSIPG(m, o);
    if (name == "B1")   return runP1MBO(m, o);
    if (name == "B2")   return runPolyVector(m, o);
    FieldRun r;
    r.method = name;
    r.error = "unknown method";
    return r;
}

// ---------------------------------------------------------------------------
bool writeFieldVTK(const std::string &path, const Mesh &m, const Eigen::VectorXcd &u) {
    std::ofstream f(path);
    if (!f) return false;

    f << "# vtk DataFile Version 3.0\ncross field\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    f << "POINTS " << m.vertices.size() << " double\n";
    for (const Point &p : m.vertices) f << p[0] << " " << p[1] << " 0\n";

    f << "CELLS " << m.triangles.size() << " " << 4 * m.triangles.size() << "\n";
    for (const Triangle &t : m.triangles) f << "3 " << t[0] << " " << t[1] << " " << t[2] << "\n";
    f << "CELL_TYPES " << m.triangles.size() << "\n";
    for (std::size_t i = 0; i < m.triangles.size(); ++i) f << "5\n";

    f << "CELL_DATA " << m.triangles.size() << "\n";
    f << "VECTORS crossU double\n";
    for (int t = 0; t < static_cast<int>(m.triangles.size()); ++t) {
        const double a = metrics::crossAngle(u[t]);
        f << std::cos(a) << " " << std::sin(a) << " 0\n";
    }
    f << "VECTORS crossV double\n";
    for (int t = 0; t < static_cast<int>(m.triangles.size()); ++t) {
        const double a = metrics::crossAngle(u[t]) + M_PI_2;
        f << std::cos(a) << " " << std::sin(a) << " 0\n";
    }
    f << "SCALARS material int 1\nLOOKUP_TABLE default\n";
    for (std::size_t t = 0; t < m.triangles.size(); ++t)
        f << (m.triangleMatId.empty() ? 1 : m.triangleMatId[t]) << "\n";

    const metrics::SingularityReport sr = metrics::singularities(m, u);
    std::vector<int> idx(m.vertices.size(), 0);
    for (const metrics::Singularity &s : sr.interior) idx[s.vertex] = s.index4;
    f << "POINT_DATA " << m.vertices.size() << "\n";
    f << "SCALARS index4 int 1\nLOOKUP_TABLE default\n";
    for (int v : idx) f << v << "\n";
    return true;
}

} // namespace paper
