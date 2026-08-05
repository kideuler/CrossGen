// Unit tests for the streamline tracer, on fields whose streamlines are known
// in closed form, so that what is being checked is the tracer and not another
// piece of the pipeline.
//
//   1  a constant field: the streamline is a straight line
//   2  a rigid rotation: the streamline is a circle, and never leaves it
//   3  v = (1, 2x): the streamline through the origin is y = x^2
//   4-7  the model of Sec. 3.2.2 imposed exactly on one singular triangle, at
//        both indices and with the singularity both centred and off-centre:
//        the ports must come out where the model puts them, and every
//        streamline crossing the triangle must stay on its hyperbola
//
// Test 3 is the one that catches a wrong branch: a cross field is only defined
// modulo pi/2, and picking the wrong quarter turn anywhere along the trace
// sends it off at a right angle, which shows up immediately as a large
// deviation from the parabola.
#include "crossfield/CrossField.hxx"
#include "tracing/FieldTracer.hxx"
#include "tracing/SeparatrixTrace.hxx"
#include "TestHelper.hxx"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <vector>

namespace {

void writePathToOBJ(const std::string &filename, const std::vector<TracePoint> &path) {
    std::ofstream out(filename);
    if (!out) return;
    for (const auto &tp : path) out << "v " << tp.global_pos[0] << " " << tp.global_pos[1] << " 0\n";
    out << "l";
    for (size_t i = 1; i <= path.size(); ++i) out << " " << i;
    out << "\n";
}

// Give every vertex the field u = e^{4 i theta(v)}.
void setField(const std::shared_ptr<CrossField> &cf, const std::shared_ptr<Mesh> &mesh,
              const std::function<double(const Point &)> &theta) {
    cf->u_k.resize(mesh->vertices.size());
    for (int i = 0; i < static_cast<int>(mesh->vertices.size()); ++i)
        cf->u_k[i] = std::polar(1.0, 4.0 * theta(mesh->vertices[i]));
}

// Follow the field from `start` in direction `dir` until it leaves the mesh.
std::vector<TracePoint> traceToBoundary(const FieldTracer &tracer, const Mesh &mesh,
                                        const Point &start, double dir, int maxSteps = 100000) {
    std::vector<TracePoint> path;
    const int f = mesh.findTriangleContainingPoint(start);
    if (f < 0) return path;

    TracePoint p0;
    p0.global_pos = start;
    p0.face_id = f;
    p0.theta = dir;
    path.push_back(p0);

    Walker w = tracer.startAt(f, start, dir);
    for (int i = 0; i < maxSteps; ++i)
        if (tracer.advance(w, path) != FieldTracer::Status::Ok) break;
    return path;
}

bool nearBoundary(const std::vector<TracePoint> &path, double radius, double tol) {
    if (path.empty()) return false;
    const Point &p = path.back().global_pos;
    return std::fabs(normP(p) - radius) <= tol;
}

} // namespace

//=============================================================================
// 1: a constant field
//=============================================================================
static bool test1_ConstantFieldTrace() {
    std::cout << "\n=== Test 1: constant field ===\n";
    auto mesh = TestHelper::createEllipse(0.0, 0.0, 1.0, 1.0, 360.0, 0.05);
    auto cf = std::make_shared<CrossField>(mesh);
    setField(cf, mesh, [](const Point &) { return M_PI_4; });

    FieldTracer tracer(cf, false);
    auto path = traceToBoundary(tracer, *mesh, Point{0.0, 0.0}, M_PI_4);
    writePathToOBJ("test1_constant_field.obj", path);
    std::cout << "  " << path.size() << " points, end ("
              << path.back().global_pos[0] << ", " << path.back().global_pos[1] << ")\n";

    if (!nearBoundary(path, 1.0, 0.05)) {
        std::cerr << "FAIL: did not reach the boundary\n";
        return false;
    }
    // A constant field has straight streamlines, so every point must sit on the
    // ray y = x.
    double worst = 0.0;
    for (const auto &tp : path)
        worst = std::max(worst, std::fabs(tp.global_pos[1] - tp.global_pos[0]));
    std::cout << "  max |y - x| = " << worst << "\n";
    if (worst > 1e-9) {
        std::cerr << "FAIL: streamline is not straight\n";
        return false;
    }
    std::cout << "PASS\n";
    return true;
}

//=============================================================================
// 2: a rigid rotation -- the streamline is a circle and closes on itself
//=============================================================================
static bool test2_RigidBodyRotation() {
    std::cout << "\n=== Test 2: rigid rotation ===\n";
    auto mesh = TestHelper::createEllipse(0.0, 0.0, 1.0, 1.0, 360.0, 0.05);
    auto cf = std::make_shared<CrossField>(mesh);
    setField(cf, mesh, [](const Point &p) {
        return (normP(p) < 1e-12) ? M_PI_2 : std::atan2(p[1], p[0]) + M_PI_2;
    });

    FieldTracer tracer(cf, false);
    const double r0 = 0.5;
    // The field is tangent to the circle of radius r0 there, so the streamline
    // is that circle: it never reaches the boundary and the walk runs until the
    // step budget stops it.
    auto path = traceToBoundary(tracer, *mesh, Point{r0, 0.0}, M_PI_2, 4000);
    writePathToOBJ("test2_rotation.obj", path);

    if (path.size() < 100) {
        std::cerr << "FAIL: only " << path.size() << " points; a closed streamline should run on\n";
        return false;
    }
    double worst = 0.0;
    for (const auto &tp : path) worst = std::max(worst, std::fabs(normP(tp.global_pos) - r0));
    std::cout << "  " << path.size() << " points, max |r - " << r0 << "| = " << worst << "\n";
    if (worst > 0.01) {
        std::cerr << "FAIL: the streamline drifted off its circle\n";
        return false;
    }
    std::cout << "PASS\n";
    return true;
}

//=============================================================================
// 3: v = (1, 2x) -- the streamline through the origin is y = x^2
//=============================================================================
static bool test3_ParabolicStreamline() {
    std::cout << "\n=== Test 3: v = (1, 2x), y = x^2 ===\n";
    auto mesh = TestHelper::createEllipse(0.0, 0.0, 1.0, 1.0, 360.0, 0.05);
    auto cf = std::make_shared<CrossField>(mesh);
    setField(cf, mesh, [](const Point &p) { return std::atan2(2.0 * p[0], 1.0); });

    FieldTracer tracer(cf, false);
    auto path = traceToBoundary(tracer, *mesh, Point{0.0, 0.0}, 0.0);
    writePathToOBJ("test3_parabolic.obj", path);

    if (path.size() < 2) {
        std::cerr << "FAIL: empty path\n";
        return false;
    }
    double worst = 0.0;
    for (const auto &tp : path) {
        const double x = tp.global_pos[0], y = tp.global_pos[1];
        worst = std::max(worst, std::fabs(y - x * x));
    }
    const double endR = normP(path.back().global_pos);
    std::cout << "  " << path.size() << " points, end radius " << endR
              << ", max |y - x^2| = " << worst << "\n";

    if (endR < 0.7) {
        std::cerr << "FAIL: did not reach the boundary\n";
        return false;
    }
    if (worst > 0.01) {
        std::cerr << "FAIL: deviation from the exact streamline is too large\n";
        return false;
    }
    std::cout << "PASS\n";
    return true;
}

//=============================================================================
// 4-7: one singular triangle carrying the model of Sec. 3.2.2 exactly
//
// The field is set to f = e^{i d theta / 4} with theta the polar angle about
// the singularity measured from the positive x axis, so the separatrices are at
// 0, 2 pi/(4-d), ... by construction and the ports the tracer fits are checked
// against those. Then a fan of streamlines is sent in through each edge and
// each one is checked against the hyperbola it is supposed to be on: under
// w = z^M with M = (4-d)/8, rho^2 sin(phi) cos(phi) is constant along it.
//=============================================================================
static bool test4_SingularTriangleSweep(double singularityIndex,
                                        std::array<double, 3> singularityBarycenter) {
    std::cout << "\n=== Singular triangle, index " << singularityIndex << ", singularity at ("
              << singularityBarycenter[0] << ", " << singularityBarycenter[1] << ", "
              << singularityBarycenter[2] << ") ===\n";

    const double R = 1.0;
    const Point v0{0.0, R};
    const Point v1{-R * std::sqrt(3.0) / 2.0, -R / 2.0};
    const Point v2{R * std::sqrt(3.0) / 2.0, -R / 2.0};
    const Point v3{-R * std::sqrt(3.0), R};
    const Point v4{0.0, -R * 1.5};
    const Point v5{R * std::sqrt(3.0), R};

    auto mesh = std::make_shared<Mesh>(std::vector<Point>{v0, v1, v2, v3, v4, v5},
                                       std::vector<Triangle>{{0, 1, 2},   // the singular one
                                                             {1, 0, 3},
                                                             {2, 1, 4},
                                                             {0, 2, 5}});

    const Point sc = v0 * singularityBarycenter[0] + v1 * singularityBarycenter[1] +
                     v2 * singularityBarycenter[2];
    const double d = (singularityIndex > 0.0) ? 1.0 : -1.0;

    auto cf = std::make_shared<CrossField>(mesh);
    setField(cf, mesh, [&](const Point &p) {
        return d * std::atan2(p[1] - sc[1], p[0] - sc[0]) / 4.0;
    });
    cf->singularTriangles.emplace_back(0, singularityIndex);

    FieldTracer tracer(cf, true);
    if (tracer.getSingularities().size() != 1) {
        std::cerr << "FAIL: the singularity was not picked up\n";
        return false;
    }
    const Singularity &s = tracer.getSingularities()[0];
    const int nPorts = s.numPorts();
    const double sector = 2.0 * M_PI / nPorts;
    const double M = nPorts / 8.0;

    std::cout << "  centre (" << s.coordinates[0] << ", " << s.coordinates[1] << ") vs exact ("
              << sc[0] << ", " << sc[1] << "), residual "
              << s.portResidual * 180.0 / M_PI << " deg\n  ports";
    for (double a : s.portAngles) std::cout << " " << a * 180.0 / M_PI;
    std::cout << "\n";

    if (nPorts != static_cast<int>(4 - d)) {
        std::cerr << "FAIL: " << nPorts << " ports, expected " << (4 - d) << "\n";
        return false;
    }
    // The centre is where the interpolated representation vector vanishes,
    // which is near the model's singularity but not equal to it: the three
    // vertex values are unit vectors, and where their affine interpolant
    // vanishes is decided by their directions rather than by where the field
    // they were sampled from is singular. So this checks it lands inside the
    // triangle and close, not that it lands exactly.
    const std::array<double, 3> bc = tracer.barycentric(0, s.coordinates);
    if (std::min({bc[0], bc[1], bc[2]}) <= 0.0) {
        std::cerr << "FAIL: the centre came out outside its own triangle\n";
        return false;
    }
    if (normP(s.coordinates - sc) > 0.2 * R) {
        std::cerr << "FAIL: the centre is " << normP(s.coordinates - sc) << " from the model's\n";
        return false;
    }
    // The field was built with a separatrix along the positive x axis, so the
    // ports must land on multiples of the sector angle. This is what pins down
    // the recovery of beta from the corner values: reading the model without
    // its rotation of the direction gives ports that are wrong by a factor of
    // -d/(4-d) and this test says so.
    for (double a : s.portAngles) {
        const double off = std::fabs(wrap_pi(a - sector * std::round(a / sector)));
        if (off > 1.0 * M_PI / 180.0) {
            std::cerr << "FAIL: port at " << a * 180.0 / M_PI << " deg is " << off * 180.0 / M_PI
                      << " deg off a multiple of " << sector * 180.0 / M_PI << "\n";
            return false;
        }
    }

    // A fan of streamlines in through each edge, each entering along the field.
    const std::array<std::array<Point, 2>, 3> edgeEnds{
        {{v0, v1}, {v1, v2}, {v2, v0}}};  // local edge e runs tri[e] -> tri[e+1]
    int traced = 0;
    double worstDrift = 0.0;

    for (int e = 0; e < 3; ++e) {
        const Point &A = edgeEnds[e][0];
        const Point &B = edgeEnds[e][1];
        for (int i = 1; i <= 50; ++i) {
            const double t = i / 51.0;
            const Point P = A * (1.0 - t) + B * t;

            // The branch of the model field at P that points into the triangle.
            const Point rel = P - s.coordinates;
            const double omega = std::atan2(rel[1], rel[0]);
            const Point inward = normalizeP(s.coordinates - P);
            double dir = 0.0;
            double bestDot = -2.0;
            for (int k = 0; k < 4; ++k) {
                const double a = d * omega / 4.0 + k * M_PI_2;
                const double dot = std::cos(a) * inward[0] + std::sin(a) * inward[1];
                if (dot > bestDot) { bestDot = dot; dir = a; }
            }

            Walker w;
            w.tri = 0;
            w.pos = P;
            w.entryEdge = e;
            w.atVertex = -1;
            w.dir = dir;
            w.crossDir = dir;

            std::vector<TracePoint> path;
            path.push_back(TracePoint{P, 0, dir});
            const FieldTracer::Status st = tracer.advance(w, path);

            if (st == FieldTracer::Status::Stuck) {
                std::cerr << "FAIL: streamline from edge " << e << " at t = " << t << " got stuck\n";
                return false;
            }
            if (path.size() < 2) {
                std::cerr << "FAIL: streamline from edge " << e << " at t = " << t
                          << " produced nothing\n";
                return false;
            }

            // Every point should be on one hyperbola of the family: with theta
            // measured from whichever port the sweep referred to, rho^2 sin phi
            // cos phi does not change along the curve. Which port that is comes
            // out of the first point.
            auto invariant = [&](const Point &q, double ref, double sigma) {
                const Point r = q - s.coordinates;
                const double rho = std::pow(normP(r), M);
                double th = sigma * (std::atan2(r[1], r[0]) - ref);
                th = std::fmod(th, 2.0 * M_PI);
                if (th < 0.0) th += 2.0 * M_PI;
                const double phi = M * th;
                return rho * rho * std::sin(phi) * std::cos(phi);
            };
            // A streamline entering along a port runs straight out along it:
            // there the hyperbola has degenerated to its own asymptote and the
            // invariant is zero in every frame, so it is checked for being
            // straight instead.
            {
                const Point d = path.back().global_pos - path.front().global_pos;
                const double len = normP(d);
                double bend = 0.0;
                if (len > 1e-12) {
                    const Point n{-d[1] / len, d[0] / len};
                    for (const auto &tp : path)
                        bend = std::max(bend, std::fabs(dotP(tp.global_pos - path.front().global_pos, n)));
                }
                if (bend < 1e-9 * std::max(len, 1e-12)) { ++traced; continue; }
            }

            // Otherwise keep whichever frame the curve is actually constant in.
            double bestSpread = std::numeric_limits<double>::max();
            for (int k = 0; k < nPorts; ++k) {
                for (double sigma : {1.0, -1.0}) {
                    const double ref = s.portAngles[k];
                    double lo = std::numeric_limits<double>::max(), hi = -lo, mean = 0.0;
                    bool usable = true;
                    for (const auto &tp : path) {
                        const double a = invariant(tp.global_pos, ref, sigma);
                        if (!std::isfinite(a) || a < 0.0) { usable = false; break; }
                        lo = std::min(lo, a);
                        hi = std::max(hi, a);
                        mean += a;
                    }
                    if (!usable || path.empty()) continue;
                    mean /= static_cast<double>(path.size());
                    if (mean < 1e-12) continue;
                    bestSpread = std::min(bestSpread, (hi - lo) / mean);
                }
            }
            if (bestSpread > 0.02) {
                std::cerr << "FAIL: streamline from edge " << e << " at t = " << t
                          << " is not on a hyperbola of the family (spread " << bestSpread << ")\n";
                return false;
            }
            worstDrift = std::max(worstDrift, bestSpread);
            ++traced;
        }
    }

    std::cout << "  " << traced << " streamlines, worst relative drift off the hyperbola "
              << worstDrift << "\n";
    std::cout << "PASS\n";
    return true;
}

//=============================================================================
int main(int argc, char **argv) {
    std::map<int, std::pair<std::string, std::function<bool()>>> tests = {
        {1, {"ConstantFieldTrace", test1_ConstantFieldTrace}},
        {2, {"RigidBodyRotation", test2_RigidBodyRotation}},
        {3, {"ParabolicStreamline", test3_ParabolicStreamline}},
        {4, {"SingularSweep_3port_centred",
             [] { return test4_SingularTriangleSweep(0.25, {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}); }}},
        {5, {"SingularSweep_3port_offcentre",
             [] { return test4_SingularTriangleSweep(0.25, {2.0 / 5.0, 2.0 / 5.0, 1.0 / 5.0}); }}},
        {6, {"SingularSweep_5port_centred",
             [] { return test4_SingularTriangleSweep(-0.25, {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}); }}},
        {7, {"SingularSweep_5port_offcentre",
             [] { return test4_SingularTriangleSweep(-0.25, {2.0 / 5.0, 2.0 / 5.0, 1.0 / 5.0}); }}},
    };

    std::vector<int> toRun;
    for (int i = 1; i < argc; ++i) {
        const int n = std::atoi(argv[i]);
        if (tests.count(n)) toRun.push_back(n);
        else std::cerr << "Unknown test number: " << argv[i] << "\n";
    }
    if (argc == 1) for (const auto &kv : tests) toRun.push_back(kv.first);

    int passed = 0, failed = 0;
    for (const int n : toRun) {
        try {
            tests[n].second() ? ++passed : ++failed;
        } catch (const std::exception &e) {
            std::cerr << "Test " << n << " (" << tests[n].first << ") threw: " << e.what() << "\n";
            ++failed;
        }
    }
    std::cout << "\nPassed " << passed << "/" << (passed + failed) << "\n";
    return failed ? 1 : 0;
}
