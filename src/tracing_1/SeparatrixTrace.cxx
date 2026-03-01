#include "tracing_1/SeparatrixTrace.hxx"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iostream>
// -----------------------------------------------------------------------------
// Separatrix tracing implementation.
//
// This implementation follows the main ideas from:
//   Viertel et al., "Instant Meshes: Interactive Field-Aligned Mesh Generation"
// (IMR 2019), Sec. 3.2/3.3.
//
// Key improvements over the baseline implementation:
//   * Heun (RK2) integration in regular triangles with robust edge clipping.
//   * Hyperbolic tracing through singular triangles via the conformal map
//     g(z)=z^{(4-d)/8}, sampled along fixed rays in the mapped quadrant.
//   * Robust stopping criteria:
//       (1) boundary,
//       (2) repeated orthogonal crossing of the same separatrix,
//       (3) singularity cutoff,
//       (+) self-intersection & tangential merge safety.
//   * Per-triangle segment indexing for scalable intersection tests.
// -----------------------------------------------------------------------------

namespace {

// Feature flags for enabling/disabling stopping criteria
static const bool ENABLE_TANGENTIAL_CROSSING_CHECK = true;  // Tangential merge detection
static const bool ENABLE_ORTHOGONAL_CROSSING_IN_SINGULAR = true;  // Orthogonal crossing check in singular triangles
static const bool ENABLE_SINGULARITY_CUTOFF = true;  // Stop when crossing another separatrix port in singular triangle
static const bool ENABLE_LIMIT_CYCLE_DETECTION = false;  // Stop if separatrix re-enters a previously visited triangle

// Integration configuration.
static const double STEP_SIZE = 0.02; // relative to average edge length
static const int MAX_STEPS = 1000000;

// Robustness tolerances.
static const double EPS_INTERSECT_T = 1e-10;      // param epsilon for segment intersections
static const double EPS_BARY = 1e-9;              // barycentric inside tolerance
static const double MIN_ADVANCE_REL = 1e-7;       // minimum accepted segment length as fraction of avg edge
static const double ANGLE_ORTHO_DOT = 0.3420201433; // |dot| <= cos(70deg) => within 20deg of orthogonal

// Tangential merge parameters.
static const double MERGE_DIR_DOT_THRESHOLD = -0.95; // directions must be strongly anti-parallel (legacy proximity check)

// Robust tangential detection parameters (per the tangential_crossing_robustness.md doc)
static const double TANGENTIAL_ABS_DOT_THRESHOLD = 0.94;  // |dot| > this means nearly parallel/anti-parallel (~20 deg)
static const double TANGENTIAL_OFFSET_ALONG_THRESHOLD = 0.35; // offset must be mostly perpendicular to streamline

// Discrete edge crossing: represents a separatrix crossing two edges of a triangle.
// Stored as (entry_edge, exit_edge) where edges are local indices 0,1,2.
struct EdgeCrossing {
    int sep_id;
    int entry_edge;  // local edge index where separatrix entered
    int exit_edge;   // local edge index where separatrix exited
    Point entry_pos;
    Point exit_pos;
    size_t path_idx; // index in separatrix path of the entry point
};

// Singular tracing sampling: log-uniform in tan(phi) over a wide range.
static const int SINGULAR_RAY_COUNT = 100;
static const double SINGULAR_LOG_TAN_RANGE = 50.0; // tan(phi) in [exp(-R), exp(R)]

// ----- small geometry helpers -----

inline double wrap_2pi(double a) {
    a = std::fmod(a, 2.0 * M_PI);
    if (a < 0.0) a += 2.0 * M_PI;
    return a;
}

inline Point addP(const Point &a, const Point &b) { return Point{a[0] + b[0], a[1] + b[1]}; }
inline Point subP(const Point &a, const Point &b) { return Point{a[0] - b[0], a[1] - b[1]}; }
inline Point mulP(const Point &a, double s) { return Point{a[0] * s, a[1] * s}; }

inline bool isOrthogonal(const Point &a, const Point &b) {
    double na = normP(a);
    double nb = normP(b);
    if (na <= 0.0 || nb <= 0.0) return false;
    double c = std::abs(dotP(a, b)) / (na * nb);
    return c <= ANGLE_ORTHO_DOT;
}

// Clamp negatives to 0 and renormalize to sum=1.
inline std::array<double, 3> sanitizeBarycentric(const std::array<double, 3> &b) {
    std::array<double, 3> out = b;
    for (double &x : out) {
        if (x < 0.0 && x > -1e-8) x = 0.0;
    }
    double s = out[0] + out[1] + out[2];
    if (std::abs(s) < 1e-24) return {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    out[0] /= s;
    out[1] /= s;
    out[2] /= s;
    return out;
}

inline bool baryInside(const std::array<double, 3> &b, double eps = EPS_BARY) {
    return b[0] >= -eps && b[1] >= -eps && b[2] >= -eps;
}

// Segment-segment intersection with parameter outputs.
// Returns true if they intersect properly (not parallel), and outputs
// t in [0,1] on p0->p1 and s in [0,1] on q0->q1.
inline bool segmentSegmentIntersectionParams(const Point &p0, const Point &p1, const Point &q0, const Point &q1,
                                             double &tOut, double &sOut, Point &ip) {
    Point r = subP(p1, p0);
    Point s = subP(q1, q0);
    double denom = cross2(r, s);
    if (std::abs(denom) < 1e-20) return false;

    Point qp = subP(q0, p0);
    double t = cross2(qp, s) / denom;
    double u = cross2(qp, r) / denom;

    if (t < -EPS_INTERSECT_T || t > 1.0 + EPS_INTERSECT_T || u < -EPS_INTERSECT_T || u > 1.0 + EPS_INTERSECT_T)
        return false;

    t = std::max(0.0, std::min(1.0, t));
    u = std::max(0.0, std::min(1.0, u));

    ip = addP(p0, mulP(r, t));
    tOut = t;
    sOut = u;
    return true;
}

// Find the first triangle edge that a segment p0->p1 exits through.
// Assumes p0 lies in (or very near) the triangle.
inline bool segmentTriangleExit(const Mesh &mesh, int triId, const Point &p0, const Point &p1,
                                int &edgeLocalOut, int &neighborOut, Point &ipOut,
                                double tEps = EPS_INTERSECT_T) {
    const Triangle &tri = mesh.triangles[triId];
    Point v[3] = {mesh.vertices[tri[0]], mesh.vertices[tri[1]], mesh.vertices[tri[2]]};

    Point r = subP(p1, p0);
    double bestT = std::numeric_limits<double>::infinity();
    int bestEdge = -1;
    Point bestIp = p0;

    for (int e = 0; e < 3; ++e) {
        int a = e;
        int b = (e + 1) % 3;
        Point q0 = v[a];
        Point q1 = v[b];

        double t, s;
        Point ip;
        if (!segmentSegmentIntersectionParams(p0, p1, q0, q1, t, s, ip)) continue;

        if (t <= tEps || t >= bestT) continue;
        // allow intersection at vertices with small tolerance
        if (s < -1e-12 || s > 1.0 + 1e-12) continue;

        bestT = t;
        bestEdge = e;
        bestIp = ip;
    }

    if (bestEdge < 0 || !std::isfinite(bestT)) return false;

    edgeLocalOut = bestEdge;
    neighborOut = mesh.triangleAdjacency[triId][bestEdge];
    ipOut = bestIp;
    return true;
}

// Average edge length using mesh.edges (unique edges).
inline double computeAverageEdgeLength(const Mesh &mesh) {
    if (mesh.edges.empty()) {
        // fallback: average over triangles
        double total = 0.0;
        size_t cnt = 0;
        for (const auto &tri : mesh.triangles) {
            for (int i = 0; i < 3; ++i) {
                const Point &a = mesh.vertices[tri[i]];
                const Point &b = mesh.vertices[tri[(i + 1) % 3]];
                total += normP(subP(b, a));
                cnt++;
            }
        }
        return cnt ? total / static_cast<double>(cnt) : 1.0;
    }

    double total = 0.0;
    for (const auto &e : mesh.edges) {
        total += normP(subP(mesh.vertices[e[1]], mesh.vertices[e[0]]));
    }
    return total / static_cast<double>(mesh.edges.size());
}

// Point-to-segment distance (and closest point).
inline std::tuple<double, Point, double> pointToSegmentDistance(const Point &p, const Point &a, const Point &b) {
    Point ab = subP(b, a);
    double len2 = dotP(ab, ab);
    if (len2 < 1e-30) {
        return {normP(subP(p, a)), a, 0.0};
    }
    double t = dotP(subP(p, a), ab) / len2;
    t = std::max(0.0, std::min(1.0, t));
    Point cp = addP(a, mulP(ab, t));
    return {normP(subP(p, cp)), cp, t};
}

// Determine which edge of a triangle a point lies on (returns -1 if not on any edge).
inline int pointOnTriangleEdge(const Mesh &mesh, int triId, const Point &p, double tol = 1e-8) {
    const Triangle &tri = mesh.triangles[triId];
    Point v[3] = {mesh.vertices[tri[0]], mesh.vertices[tri[1]], mesh.vertices[tri[2]]};
    
    for (int e = 0; e < 3; ++e) {
        auto [dist, cp, t] = pointToSegmentDistance(p, v[e], v[(e + 1) % 3]);
        if (dist < tol) {
            return e;
        }
    }
    return -1;
}

// Complex power in polar form, with the argument wrapped to [0,2pi).
inline std::complex<double> complexPowPolar(const std::complex<double> &z, double p) {
    double r = std::abs(z);
    if (r <= 0.0) return std::complex<double>(0.0, 0.0);
    double theta = wrap_2pi(std::arg(z));
    return std::polar(std::pow(r, p), theta * p);
}

} // namespace

// -----------------------------------------------------------------------------
// SeparatrixTrace methods
// -----------------------------------------------------------------------------

SeparatrixTrace::SeparatrixTrace(std::shared_ptr<CrossField> cf) : crossField(cf) {
    separatrices.clear();
}

std::array<double, 3> SeparatrixTrace::globalToBarycentric(int triangleIndex, const Point &p) {
    const Mesh &mesh = crossField->getMesh();
    const Triangle &tri = mesh.triangles[triangleIndex];

    const Point &v0 = mesh.vertices[tri[0]];
    const Point &v1 = mesh.vertices[tri[1]];
    const Point &v2 = mesh.vertices[tri[2]];

    // Robust barycentric via area coordinates.
    Point v0v1 = subP(v1, v0);
    Point v0v2 = subP(v2, v0);
    Point v0p = subP(p, v0);

    double denom = cross2(v0v1, v0v2);
    if (std::abs(denom) < 1e-30) {
        return {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    }

    double l1 = cross2(v0p, v0v2) / denom; // weight for v1
    double l2 = cross2(v0v1, v0p) / denom; // weight for v2
    double l0 = 1.0 - l1 - l2;

    return {l0, l1, l2};
}

std::tuple<Point, int, double, int> SeparatrixTrace::rayEdgeIntersection(int triangleIndex, const Point &origin, double direction) {
    const Mesh &mesh = crossField->getMesh();
    const Triangle &tri = mesh.triangles[triangleIndex];

    Point dir = {std::cos(direction), std::sin(direction)};

    Point v[3] = {mesh.vertices[tri[0]], mesh.vertices[tri[1]], mesh.vertices[tri[2]]};

    double bestT = std::numeric_limits<double>::infinity();
    int bestEdge = -1;
    double bestU = 0.0;
    Point bestIp = origin;

    for (int e = 0; e < 3; ++e) {
        int a = e;
        int b = (e + 1) % 3;
        Point e0 = v[a];
        Point e1 = v[b];
        Point edge = subP(e1, e0);

        // Solve origin + t*dir = e0 + u*edge
        double denom = cross2(dir, edge);
        if (std::abs(denom) < 1e-20) continue;

        Point eo = subP(e0, origin);
        double t = cross2(eo, edge) / denom;
        double u = cross2(eo, dir) / denom;

        if (t <= EPS_INTERSECT_T) continue;
        if (u < -1e-12 || u > 1.0 + 1e-12) continue;

        if (t < bestT) {
            bestT = t;
            bestEdge = e;
            bestU = std::max(0.0, std::min(1.0, u));
            bestIp = addP(origin, mulP(dir, t));
        }
    }

    int neighbor = -1;
    if (bestEdge >= 0) neighbor = mesh.triangleAdjacency[triangleIndex][bestEdge];

    return {bestIp, bestEdge, bestU, neighbor};
}

// Interpolate the cross field direction at a barycentric coordinate.
// Returns a representative angle that is continuous w.r.t. prevAngle.
double SeparatrixTrace::interpolateFieldDirection(int triangleIndex, const std::array<double, 3> &bary, double prevAngle) {
    const Mesh &mesh = crossField->getMesh();
    const Triangle &tri = mesh.triangles[triangleIndex];

    // Prefer complex interpolation of the representation vector u^4.
    // This avoids branch-cut issues when vertex arguments straddle +/-pi.
    std::complex<double> u4_interp(0.0, 0.0);
    for (int i = 0; i < 3; ++i) {
        u4_interp += crossField->u_k[tri[i]] * bary[i];
    }

    double baseAngle = 0.0;
    const double mag = std::abs(u4_interp);

    if (mag > 1e-8) {
        baseAngle = std::arg(u4_interp) / 4.0;
    } else {
        // Degenerate interpolation (near cancellation). Fall back to a robust
        // unwrapped angle interpolation on the triangle.
        double raw[3];
        for (int i = 0; i < 3; ++i) raw[i] = std::arg(crossField->u_k[tri[i]]);

        double bestB[3] = {raw[0], raw[1], raw[2]};
        double bestE = std::numeric_limits<double>::infinity();

        // Fix k0 = 0; search small shifts for the other two vertices.
        for (int k1 = -1; k1 <= 1; ++k1) {
            for (int k2 = -1; k2 <= 1; ++k2) {
                double b0 = raw[0];
                double b1 = raw[1] + 2.0 * M_PI * static_cast<double>(k1);
                double b2 = raw[2] + 2.0 * M_PI * static_cast<double>(k2);

                double e01 = b1 - b0;
                double e12 = b2 - b1;
                double e20 = b0 - b2;
                double E = e01 * e01 + e12 * e12 + e20 * e20;

                if (E < bestE) {
                    bestE = E;
                    bestB[0] = b0;
                    bestB[1] = b1;
                    bestB[2] = b2;
                }
            }
        }

        double betaInterp = bestB[0] * bary[0] + bestB[1] * bary[1] + bestB[2] * bary[2];
        baseAngle = betaInterp / 4.0;
    }

    // Choose the closest cross branch to the previous direction to maintain
    // consistent orientation of the streamline.
    double bestAngle = baseAngle;
    double bestDiff = std::abs(wrap_pi(bestAngle - prevAngle));
    for (int k = 1; k < 4; ++k) {
        double cand = baseAngle + static_cast<double>(k) * M_PI_2;
        double diff = std::abs(wrap_pi(cand - prevAngle));
        if (diff < bestDiff) {
            bestDiff = diff;
            bestAngle = cand;
        }
    }

    return bestAngle;
}

// Compute alpha: nearest cross direction to the reference axis from triangle barycenter to a vertex.
double SeparatrixTrace::getFieldAlignmentAngle(int triangleIndex, int vertexLocalIndex) {
    const Mesh &mesh = crossField->getMesh();
    const Triangle &tri = mesh.triangles[triangleIndex];

    int vIdx = tri[vertexLocalIndex];
    std::complex<double> u4 = crossField->u_k[vIdx];
    double fieldAngle = std::arg(u4) / 4.0;

    // Triangle barycenter.
    Point bc{0.0, 0.0};
    for (int i = 0; i < 3; ++i) {
        bc[0] += mesh.vertices[tri[i]][0];
        bc[1] += mesh.vertices[tri[i]][1];
    }
    bc[0] /= 3.0;
    bc[1] /= 3.0;

    Point refVec = subP(mesh.vertices[vIdx], bc);
    double refAngle = std::atan2(refVec[1], refVec[0]);

    // Pick the cross branch closest to the reference axis.
    double best = fieldAngle;
    double bestDiff = std::abs(wrap_pi(best - refAngle));
    for (int k = 1; k < 4; ++k) {
        double cand = fieldAngle + static_cast<double>(k) * M_PI_2;
        double diff = std::abs(wrap_pi(cand - refAngle));
        if (diff < bestDiff) {
            bestDiff = diff;
            best = cand;
        }
    }

    double alpha = wrap_pi(best - refAngle);
    // Normalize to (-pi/4, pi/4]
    while (alpha > M_PI_4) alpha -= M_PI_2;
    while (alpha <= -M_PI_4) alpha += M_PI_2;

    return alpha;
}

std::vector<std::pair<TracePoint, double>> SeparatrixTrace::computeSingularityPorts(int triangleIndex, double crossFieldIndex) {
    std::vector<std::pair<TracePoint, double>> ports;

    const Mesh &mesh = crossField->getMesh();
    const Triangle &tri = mesh.triangles[triangleIndex];

    // Barycenter start.
    Point bc{0.0, 0.0};
    for (int i = 0; i < 3; ++i) {
        bc[0] += mesh.vertices[tri[i]][0];
        bc[1] += mesh.vertices[tri[i]][1];
    }
    bc[0] /= 3.0;
    bc[1] /= 3.0;

    // compute actual singularity center by solving 2x2 linear system
    Eigen::Matrix2d A;
    std::complex<double> u1 = crossField->u_k[tri[0]];
    std::complex<double> u2 = crossField->u_k[tri[1]];
    std::complex<double> u3 = crossField->u_k[tri[2]];

    A << std::real(u1) - std::real(u3), std::real(u2) - std::real(u3),
         std::imag(u1) - std::imag(u3), std::imag(u2) - std::imag(u3);

    Eigen::Vector2d b;
    b << -std::real(u3), -std::imag(u3);
    // uv is solution
    Eigen::Vector2d uv = A.colPivHouseholderQr().solve(b);

    // find the corresponding point in the triangle
    Point p1 = mesh.vertices[tri[0]];
    Point p2 = mesh.vertices[tri[1]];
    Point p3 = mesh.vertices[tri[2]];
    Point intersection = addP(addP(mulP(p1, uv[0]), mulP(p2, uv[1])), mulP(p3, 1.0 - uv[0] - uv[1]));

    TracePoint tp;
    tp.face_id = triangleIndex;
    tp.global_pos = intersection;
    tp.barycentric = {uv[0], uv[1], 1.0 - uv[0] - uv[1]};

    // Align to reference axis through vertex 0.
    double alpha = getFieldAlignmentAngle(triangleIndex, 0);

    // d from index: index = d/4 => d = 4*index.
    int d = static_cast<int>(std::round(4.0 * crossFieldIndex));
    int n = 4 - d; // number of separatrices (signed)
    int numPorts = std::abs(n);
    if (numPorts <= 0) return ports;

    // Reference axis angle (barycenter -> vertex 0).
    Point refVec = subP(mesh.vertices[tri[0]], bc);
    double refAngle = std::atan2(refVec[1], refVec[0]);

    // Distribute ports.
    double step = (n != 0) ? (2.0 * M_PI / static_cast<double>(n)) : 0.0;

    for (int k = 0; k < numPorts; ++k) {
        double theta_port = alpha + static_cast<double>(k) * step;
        double absoluteAngle = wrap_pi(refAngle + theta_port);
        ports.emplace_back(tp, absoluteAngle);
    }

    return ports;
}

void SeparatrixTrace::initializeSeparatrices() {
    separatrices.clear();

    const auto &singularities = crossField->singularTriangles;

    int separatrixId = 0;
    for (size_t singIdx = 0; singIdx < singularities.size(); ++singIdx) {
        int triangleIndex = singularities[singIdx].first;
        double crossFieldIndex = singularities[singIdx].second;

        auto ports = computeSingularityPorts(triangleIndex, crossFieldIndex);

        for (const auto &pr : ports) {
            const TracePoint &startPoint = pr.first;
            double direction = pr.second;

            Separatrix sep;
            sep.id = separatrixId++;
            sep.origin_singularity_id = static_cast<int>(singIdx);
            sep.active = true;
            sep.termination_reason = TerminationReason::RUNNING;
            sep.merged_with_id = -1;

            sep.path.push_back(startPoint);

            auto [intersection, edgeIdx, edgeParam, neighborTri] =
                rayEdgeIntersection(triangleIndex, startPoint.global_pos, direction);

            if (edgeIdx >= 0) {
                TracePoint edgePoint;
                edgePoint.global_pos = intersection;

                if (neighborTri >= 0) {
                    edgePoint.face_id = neighborTri;
                    edgePoint.barycentric = sanitizeBarycentric(globalToBarycentric(neighborTri, intersection));
                } else {
                    edgePoint.face_id = triangleIndex;
                    edgePoint.barycentric = sanitizeBarycentric(globalToBarycentric(triangleIndex, intersection));
                    sep.active = false;
                    sep.termination_reason = TerminationReason::HIT_BOUNDARY;
                }

                sep.path.push_back(edgePoint);
            } else {
                sep.active = false;
                sep.termination_reason = TerminationReason::MAX_STEPS_REACHED;
            }

            separatrices.push_back(sep);
        }
    }
}

// -----------------------------------------------------------------------------
// Main Trace
// -----------------------------------------------------------------------------

void SeparatrixTrace::trace() {
    const Mesh &mesh = crossField->getMesh();

    // Build singularity lookup.
    singularityMap.clear();
    singularityMap.reserve(crossField->singularTriangles.size());
    for (const auto &p : crossField->singularTriangles) {
        singularityMap[p.first] = p.second;
    }

    const double avgEdgeLen = computeAverageEdgeLength(mesh);
    const double baseStep = std::max(1e-12, STEP_SIZE * (avgEdgeLen > 0.0 ? avgEdgeLen : 1.0));
    const double minAdvance = std::max(1e-12, MIN_ADVANCE_REL * (avgEdgeLen > 0.0 ? avgEdgeLen : 1.0));

    // Precompute per-triangle minimum edge length for step clamping.
    std::vector<double> triMinEdge(mesh.triangles.size(), std::numeric_limits<double>::infinity());
    for (size_t t = 0; t < mesh.triangles.size(); ++t) {
        const Triangle &tri = mesh.triangles[t];
        Point v[3] = {mesh.vertices[tri[0]], mesh.vertices[tri[1]], mesh.vertices[tri[2]]};
        double m = std::numeric_limits<double>::infinity();
        for (int e = 0; e < 3; ++e) {
            m = std::min(m, normP(subP(v[(e + 1) % 3], v[e])));
        }
        if (std::isfinite(m) && m > 0.0) triMinEdge[t] = m;
    }

    // Precompute singular triangle port segments and angles.
    struct SingInfo {
        Point center{0.0, 0.0};
        int d = 0;
        int n = 0; // n = 4 - d
        double radius = 0.0; // max distance from center to vertices
        std::vector<double> portAngles; // [0,2pi)
        std::vector<std::pair<Point, Point>> portSegments; // center -> boundary
    };

    std::unordered_map<int, SingInfo> singInfo;
    singInfo.reserve(singularityMap.size());

    for (const auto &kv : singularityMap) {
        int triIdx = kv.first;
        double crossIdx = kv.second;
        SingInfo info;
        info.d = static_cast<int>(std::round(4.0 * crossIdx));
        info.n = 4 - info.d;

        const Triangle &tri = mesh.triangles[triIdx];
        // compute actual singularity center by solving 2x2 linear system
    Eigen::Matrix2d A;
    std::complex<double> u1 = crossField->u_k[tri[0]];
    std::complex<double> u2 = crossField->u_k[tri[1]];
    std::complex<double> u3 = crossField->u_k[tri[2]];

    A << std::real(u1) - std::real(u3), std::real(u2) - std::real(u3),
         std::imag(u1) - std::imag(u3), std::imag(u2) - std::imag(u3);

    Eigen::Vector2d b;
    b << -std::real(u3), -std::imag(u3);
    // uv is solution
    Eigen::Vector2d uv = A.colPivHouseholderQr().solve(b);

    // find the corresponding point in the triangle
    Point p1 = mesh.vertices[tri[0]];
    Point p2 = mesh.vertices[tri[1]];
    Point p3 = mesh.vertices[tri[2]];
    Point intersection = addP(addP(mulP(p1, uv[0]), mulP(p2, uv[1])), mulP(p3, 1.0 - uv[0] - uv[1]));
    info.center = intersection;

    double r = 0.0;
    for (int i = 0; i < 3; ++i) {
        r = std::max(r, normP(subP(mesh.vertices[tri[i]], intersection)));
        }
        info.radius = r;

        singInfo.emplace(triIdx, info);
    }

    for (const auto &sep : separatrices) {
        if (sep.path.size() < 2) continue;
        int originTri = sep.path.front().face_id;
        auto it = singInfo.find(originTri);
        if (it == singInfo.end()) continue;

        Point p0 = sep.path[0].global_pos;
        Point p1 = sep.path[1].global_pos;
        it->second.center = p0;
        it->second.portSegments.emplace_back(p0, p1);
        it->second.portAngles.push_back(wrap_2pi(std::atan2(p1[1] - p0[1], p1[0] - p0[0])));
    }

    for (auto &kv : singInfo) {
        auto &angles = kv.second.portAngles;
        std::sort(angles.begin(), angles.end());
    }

    // Fixed rays in the mapped quadrant (log-uniform tan sampling).
    std::vector<double> singularRays;
    singularRays.reserve(SINGULAR_RAY_COUNT);
    for (int i = 0; i < SINGULAR_RAY_COUNT; ++i) {
        double t = (static_cast<double>(i) + 0.5) / static_cast<double>(SINGULAR_RAY_COUNT);
        double logTan = (2.0 * t - 1.0) * SINGULAR_LOG_TAN_RANGE;
        double tanv = std::exp(logTan);
        singularRays.push_back(std::atan(tanv));
    }
    std::sort(singularRays.begin(), singularRays.end());

    const double mergeDist = std::max(2.0 * baseStep, 0.10 * (avgEdgeLen > 0.0 ? avgEdgeLen : 1.0));

    // Per-triangle segment index for fast intersection tests.
    struct SegmentRef {
        int sep_id;
        size_t start_idx; // index of segment start point in its separatrix path when inserted
        Point a;
        Point b;
    };

    std::vector<std::vector<SegmentRef>> segIndex(mesh.triangles.size());

    // Per-triangle edge crossing index for discrete tangential merge detection.
    // Key: triangle index -> list of edge crossings (entry_edge, exit_edge) by separatrices.
    std::vector<std::vector<EdgeCrossing>> edgeCrossingIndex(mesh.triangles.size());

    // Track the entry edge for each active separatrix in its current triangle
    // Key: separatrix id -> (current triangle, entry edge, entry position)
    std::unordered_map<int, std::tuple<int, int, Point>> sepEntryState;

    // Seed index with initial segments (from initialization).
    const double edgeTolSeed = std::max(1e-8, 1e-4 * avgEdgeLen);
    for (const auto &sep : separatrices) {
        for (size_t i = 0; i + 1 < sep.path.size(); ++i) {
            int triId = sep.path[i].face_id;
            if (triId < 0 || triId >= static_cast<int>(mesh.triangles.size())) continue;
            Point a = sep.path[i].global_pos;
            Point b = sep.path[i + 1].global_pos;
            if (normP(subP(b, a)) < minAdvance) continue;
            segIndex[triId].push_back(SegmentRef{sep.id, i, a, b});
        }
        
        // Initialize entry state for each separatrix from its last point
        if (sep.path.size() >= 2 && sep.active) {
            const TracePoint &lastPt = sep.path.back();
            int tri = lastPt.face_id;
            if (tri >= 0 && tri < static_cast<int>(mesh.triangles.size())) {
                int entryEdge = pointOnTriangleEdge(mesh, tri, lastPt.global_pos, edgeTolSeed);
                if (entryEdge >= 0) {
                    sepEntryState[sep.id] = std::make_tuple(tri, entryEdge, lastPt.global_pos);
                }
            }
        }
    }

    // Trace each separatrix.
    for (auto &sep : separatrices) {
        if (!sep.active) continue;
        if (sep.path.size() < 2) {
            sep.active = false;
            sep.termination_reason = TerminationReason::MAX_STEPS_REACHED;
            std::cout << "Separatrix " << sep.id << " has insufficient initial path points; skipping." << std::endl;
            continue;
        }

        // Track repeated orthogonal crossings.
        std::unordered_map<int, int> crossCount;
        std::unordered_map<int, Point> lastCrossPoint;

        // Track visited triangles for limit cycle detection.
        std::unordered_set<int> visitedTriangles;
        // Initialize with triangles already in the path.
        for (const auto &pt : sep.path) {
            if (pt.face_id >= 0) visitedTriangles.insert(pt.face_id);
        }

        // Previous direction based on last segment.
        Point pPrev0 = sep.path[sep.path.size() - 2].global_pos;
        Point pPrev1 = sep.path[sep.path.size() - 1].global_pos;
        Point initDir = subP(pPrev1, pPrev0);
        if (normP(initDir) < 1e-30) initDir = Point{1.0, 0.0};
        double prevAngle = std::atan2(initDir[1], initDir[0]);

        int stepCount = 0;

        while (sep.active && stepCount < MAX_STEPS) {
            TracePoint currentPt = sep.path.back();
            int currentTri = currentPt.face_id;

            if (currentTri < 0 || currentTri >= static_cast<int>(mesh.triangles.size())) {
                sep.active = false;
                sep.termination_reason = TerminationReason::MAX_STEPS_REACHED;
                std::cout << "Separatrix " << sep.id << " has exited the mesh; stopping trace." << std::endl;
                break;
            }

            struct StepResult {
                std::vector<TracePoint> pts; // excluding currentPt
                bool hitBoundary = false;
            } stepRes;

            auto appendPointInTri = [&](int triId, const Point &pos) {
                TracePoint tp;
                tp.face_id = triId;
                tp.global_pos = pos;
                tp.barycentric = sanitizeBarycentric(globalToBarycentric(triId, pos));
                stepRes.pts.push_back(tp);
            };

            
auto doRegularAdvance = [&]() {
                double h = baseStep;
                if (std::isfinite(triMinEdge[currentTri])) h = std::min(h, 0.25 * triMinEdge[currentTri]);

                // Direction at current point (choose branch closest to prevAngle).
                double theta1 = interpolateFieldDirection(currentTri, currentPt.barycentric, prevAngle);
                Point dir1{std::cos(theta1), std::sin(theta1)};

                // If the point lies on/near an edge, a tiny nudge along the direction
                // can end up in a different adjacent triangle. Use that to choose the
                // correct "side" for the step (prevents spurious t=0 edge exits).
                {
                    double nudge = std::max(1e-6 * avgEdgeLen, 1e-12);
                    Point nudgedPos = addP(currentPt.global_pos, mulP(dir1, nudge));

                    std::array<double, 3> bHere = globalToBarycentric(currentTri, nudgedPos);
                    if (!baryInside(bHere, EPS_BARY)) {
                        for (int e = 0; e < 3; ++e) {
                            int nb = mesh.triangleAdjacency[currentTri][e];
                            if (nb < 0) continue;
                            std::array<double, 3> bNb = globalToBarycentric(nb, nudgedPos);
                            if (baryInside(bNb, EPS_BARY)) {
                                currentTri = nb;
                                currentPt.face_id = nb;
                                currentPt.barycentric = sanitizeBarycentric(bNb);

                                // Recompute direction in the chosen triangle.
                                theta1 = interpolateFieldDirection(currentTri, currentPt.barycentric, prevAngle);
                                dir1 = Point{std::cos(theta1), std::sin(theta1)};
                                break;
                            }
                        }
                    }
                }

                // Heun predictor (Euler step).
                Point predictedPos = addP(currentPt.global_pos, mulP(dir1, h));

                int triPred = currentTri;
                std::array<double, 3> baryPred = globalToBarycentric(currentTri, predictedPos);

                // If the predictor left the current triangle, try evaluating v2 in the
                // adjacent triangle it moved into (1-edge crossing under our step clamp).
                if (!baryInside(baryPred, EPS_BARY)) {
                    int edgeLocal = -1;
                    int neighbor = -1;
                    Point ip;
                    if (segmentTriangleExit(mesh, currentTri, currentPt.global_pos, predictedPos,
                                            edgeLocal, neighbor, ip, EPS_INTERSECT_T) &&
                        neighbor >= 0) {
                        triPred = neighbor;
                        baryPred = globalToBarycentric(triPred, predictedPos);
                    }
                }

                double theta2 = interpolateFieldDirection(triPred, baryPred, theta1);
                Point dir2{std::cos(theta2), std::sin(theta2)};

                // Heun corrector: x_{n+1} = x_n + (h/2) * (v1 + v2).
                Point vAvg = mulP(addP(dir1, dir2), 0.5);
                if (normP(vAvg) < 1e-30) vAvg = dir1;

                Point nextPos = addP(currentPt.global_pos, mulP(vAvg, h));
                std::array<double, 3> baryNext = globalToBarycentric(currentTri, nextPos);

                if (baryInside(baryNext, EPS_BARY)) {
                    appendPointInTri(currentTri, nextPos);
                    return;
                }

                // Clip to boundary and cross edge.
                int edgeLocal = -1;
                int neighbor = -1;
                Point ip;
                if (!segmentTriangleExit(mesh, currentTri, currentPt.global_pos, nextPos, edgeLocal, neighbor, ip, EPS_INTERSECT_T)) {
                    // Fallback to ray intersection.
                    auto [ip2, e2, u2, nb2] = rayEdgeIntersection(currentTri, currentPt.global_pos, std::atan2(vAvg[1], vAvg[0]));
                    if (e2 < 0) {
                        sep.active = false;
                        sep.termination_reason = TerminationReason::MAX_STEPS_REACHED;

                        return;
                    }
                    ip = ip2;
                    neighbor = nb2;
                }

                if (neighbor < 0) {
                    appendPointInTri(currentTri, ip);
                    stepRes.hitBoundary = true;
                    return;
                }

                TracePoint tp;
                tp.face_id = neighbor;
                tp.global_pos = ip;
                tp.barycentric = sanitizeBarycentric(globalToBarycentric(neighbor, ip));
                stepRes.pts.push_back(tp);
            };

            auto doSingularAdvance = [&]() {
                auto it = singInfo.find(currentTri);
                if (it == singInfo.end()) {
                    doRegularAdvance();
                    return;
                }

                const SingInfo &info = it->second;
                const int n = info.n;
                const int d = info.d;
                (void)d;

                if (n == 0 || info.portAngles.empty()) {
                    doRegularAdvance();
                    return;
                }

                const Point a = info.center;

                std::complex<double> zeta(currentPt.global_pos[0] - a[0], currentPt.global_pos[1] - a[1]);
                double zetaAbs = std::abs(zeta);
                if (zetaAbs < 1e-30) {
                    doRegularAdvance();
                    return;
                }

                double theta_q = wrap_2pi(std::arg(zeta));

                // s0: nearest separatrix clockwise from q.
                double s0 = info.portAngles.back();
                auto ub = std::upper_bound(info.portAngles.begin(), info.portAngles.end(), theta_q);
                if (ub == info.portAngles.begin()) {
                    s0 = info.portAngles.back();
                } else {
                    s0 = *(ub - 1);
                }

                std::complex<double> rot = std::polar(1.0, -s0);
                std::complex<double> zeta_r = zeta * rot;

                const double p = static_cast<double>(n) / 8.0;
                if (p == 0.0) {
                    doRegularAdvance();
                    return;
                }

                std::complex<double> wq = complexPowPolar(zeta_r, p);

                // Mapped coordinates in (approx) first quadrant.
                double xq = wq.real();
                double yq = wq.imag();

                const double quadEps = 1e-16;
                xq = std::max(xq, quadEps);
                yq = std::max(yq, quadEps);

                double A = xq * yq;
                if (!(A > 0.0) || !std::isfinite(A)) {
                    doRegularAdvance();
                    return;
                }

                double phi_q = std::atan2(yq, xq);

                // Determine monotone direction along the hyperbola using a small probe.
                Point dirPrev{std::cos(prevAngle), std::sin(prevAngle)};
                double probeLen = std::max(1e-6 * (avgEdgeLen > 0.0 ? avgEdgeLen : 1.0), 1e-12);
                Point probePos = addP(currentPt.global_pos, mulP(dirPrev, probeLen));

                std::complex<double> zetaProbe(probePos[0] - a[0], probePos[1] - a[1]);
                std::complex<double> wProbe = complexPowPolar(zetaProbe * rot, p);

                double xp = std::max(wProbe.real(), quadEps);
                double yp = std::max(wProbe.imag(), quadEps);
                double phi_probe = std::atan2(yp, xp);

                bool incPhi = (phi_probe >= phi_q);

                const double invExp = 8.0 / static_cast<double>(n);
                const double maxR = std::max(1000.0 * info.radius, 1000.0 * (avgEdgeLen > 0.0 ? avgEdgeLen : 1.0));

                auto evalPhiPoint = [&](double phi, Point &outPos, std::array<double, 3> &outBary) -> bool {
                    double t = std::tan(phi);
                    if (!(t > 0.0) || !(A > 0.0)) return false;

                    // Compute (x,y) on hyperbola xy=A along ray at angle phi.
                    double sqrtA = std::sqrt(A);
                    double sqrtT = std::sqrt(t);
                    double x = sqrtA / sqrtT;
                    double y = sqrtA * sqrtT;

                    if (!std::isfinite(x) || !std::isfinite(y)) {
                        // treat as a very far point
                        x = std::min(std::max(x, quadEps), 1e150);
                        y = std::min(std::max(y, quadEps), 1e150);
                    }

                    std::complex<double> w(x, y);

                    // zeta_r = w^{8/n} (cap magnitude for robustness)
                    double r_w = std::abs(w);
                    double theta_w = std::atan2(y, x);
                    double r_z = std::pow(r_w, invExp);
                    if (!std::isfinite(r_z) || r_z > maxR) r_z = maxR;
                    double theta_z = theta_w * invExp;

                    std::complex<double> zeta_r_phi = std::polar(r_z, theta_z);
                    std::complex<double> zeta_phi = zeta_r_phi * std::polar(1.0, s0);

                    outPos = Point{a[0] + zeta_phi.real(), a[1] + zeta_phi.imag()};
                    outBary = globalToBarycentric(currentTri, outPos);
                    return true;
                };

                // Find starting ray index.
                int startIdx = 0;
                if (incPhi) {
                    startIdx = static_cast<int>(std::lower_bound(singularRays.begin(), singularRays.end(), phi_q) - singularRays.begin());
                    startIdx = std::max(0, std::min(startIdx, static_cast<int>(singularRays.size()) - 1));
                } else {
                    startIdx = static_cast<int>(std::upper_bound(singularRays.begin(), singularRays.end(), phi_q) - singularRays.begin()) - 1;
                    startIdx = std::max(0, std::min(startIdx, static_cast<int>(singularRays.size()) - 1));
                }

                auto inRange = [&](int idx) { return idx >= 0 && idx < static_cast<int>(singularRays.size()); };

                Point lastInside = currentPt.global_pos;
                std::array<double, 3> lastBary = currentPt.barycentric;

                int idx = startIdx;
                while (inRange(idx)) {
                    double phi = singularRays[idx];
                    if (std::abs(phi - phi_q) < 1e-15) {
                        idx += (incPhi ? 1 : -1);
                        continue;
                    }

                    Point candPos;
                    std::array<double, 3> candBary;
                    if (!evalPhiPoint(phi, candPos, candBary)) {
                        idx += (incPhi ? 1 : -1);
                        continue;
                    }

                    if (baryInside(candBary, EPS_BARY)) {
                        TracePoint tp;
                        tp.face_id = currentTri;
                        tp.global_pos = candPos;
                        tp.barycentric = sanitizeBarycentric(candBary);
                        stepRes.pts.push_back(tp);
                        lastInside = candPos;
                        lastBary = tp.barycentric;
                        idx += (incPhi ? 1 : -1);
                        continue;
                    }

                    // Leaving the triangle: intersect chord to outside point.
                    int edgeLocal = -1;
                    int neighbor = -1;
                    Point ip;
                    if (!segmentTriangleExit(mesh, currentTri, lastInside, candPos, edgeLocal, neighbor, ip, EPS_INTERSECT_T)) {
                        // If chord intersection failed, try from current point.
                        Point dir = normalizeP(subP(candPos, lastInside));
                        auto [ip2, e2, u2, nb2] = rayEdgeIntersection(currentTri, lastInside, std::atan2(dir[1], dir[0]));
                        if (e2 < 0) {
                            sep.active = false;
                            sep.termination_reason = TerminationReason::MAX_STEPS_REACHED;
                            return;
                        }
                        ip = ip2;
                        neighbor = nb2;
                    }

                    if (neighbor < 0) {
                        appendPointInTri(currentTri, ip);
                        stepRes.hitBoundary = true;
                        return;
                    }

                    TracePoint tp;
                    tp.face_id = neighbor;
                    tp.global_pos = ip;
                    tp.barycentric = sanitizeBarycentric(globalToBarycentric(neighbor, ip));
                    stepRes.pts.push_back(tp);
                    return;
                }

                // Could not find an exit via sampling; fall back.
                doRegularAdvance();
            };

            if (singularityMap.find(currentTri) != singularityMap.end()) {
                doSingularAdvance();
            } else {
                doRegularAdvance();
            }

            if (stepRes.pts.empty()) {
                sep.active = false;
                sep.termination_reason = TerminationReason::MAX_STEPS_REACHED;

                break;
            }

            // Process each generated point as a segment, applying stopping criteria.
            for (size_t k = 0; k < stepRes.pts.size(); ++k) {
                TracePoint nextPt = stepRes.pts[k];

                Point segVec = subP(nextPt.global_pos, currentPt.global_pos);
                double segLen = normP(segVec);
                if (segLen < minAdvance) {
                    sep.active = false;
                    sep.termination_reason = TerminationReason::MAX_STEPS_REACHED;
                    
                    break;
                }

                // --- Stopping criterion 3: singularity cutoff (Sec. 3.3)
                if (ENABLE_SINGULARITY_CUTOFF) {
                auto singIt = singInfo.find(currentPt.face_id);
                if (singIt != singInfo.end()) {
                    const auto &ports = singIt->second.portSegments;
                    for (const auto &ps : ports) {
                        double tI, sI;
                        Point ip;
                        if (!segmentSegmentIntersectionParams(currentPt.global_pos, nextPt.global_pos, ps.first, ps.second, tI, sI, ip))
                            continue;

                        if (tI <= EPS_INTERSECT_T || tI >= 1.0 - EPS_INTERSECT_T) continue;
                        if (sI <= EPS_INTERSECT_T || sI >= 1.0 - EPS_INTERSECT_T) continue;

                        if (!isOrthogonal(segVec, subP(ps.second, ps.first))) continue;

                        // Clip at intersection point.
                        TracePoint ipPt;
                        ipPt.face_id = currentPt.face_id;
                        ipPt.global_pos = ip;
                        ipPt.barycentric = sanitizeBarycentric(globalToBarycentric(ipPt.face_id, ip));

                        // Add clipped segment to path & index.
                        sep.path.push_back(ipPt);
                        if (currentTri >= 0 && currentTri < static_cast<int>(segIndex.size())) {
                            segIndex[currentTri].push_back(SegmentRef{sep.id, sep.path.size() - 2, currentPt.global_pos, ipPt.global_pos});
                        }

                        sep.active = false;
                        sep.termination_reason = TerminationReason::SINGULARITY_CUTOFF;
                        break;
                    }
                    if (!sep.active) break;
                }
                } // end if (ENABLE_SINGULARITY_CUTOFF)

                // --- Stopping criterion 2: repeated orthogonal crossing of the same separatrix.
                // Also stop on any self-intersection.
                bool terminatedByIntersection = false;
                Point intersectionPoint{0.0, 0.0};
                int hitSepId = -1;

                // Check if current triangle is singular (for skipping orthogonal check)
                bool inSingularTri = (singularityMap.find(currentTri) != singularityMap.end());

                // Only segments in the current triangle can intersect interiorly.
                const auto &candidates = segIndex[currentTri];
                for (const auto &sr : candidates) {
                    // Skip the most recent segment of this separatrix to avoid immediate endpoint hits.
                    if (sr.sep_id == sep.id) {
                        if (sep.path.size() >= 2 && sr.start_idx >= sep.path.size() - 2) continue;
                    }

                    double tI, sI;
                    Point ip;
                    if (!segmentSegmentIntersectionParams(currentPt.global_pos, nextPt.global_pos, sr.a, sr.b, tI, sI, ip)) continue;

                    // Ignore intersections at segment endpoints.
                    if (tI <= EPS_INTERSECT_T || tI >= 1.0 - EPS_INTERSECT_T) continue;
                    if (sI <= EPS_INTERSECT_T || sI >= 1.0 - EPS_INTERSECT_T) continue;

                    if (sr.sep_id == sep.id) {
                        // Self-intersection: stop regardless of orthogonality.
                        terminatedByIntersection = true;
                        hitSepId = sep.id;
                        intersectionPoint = ip;
                        break;
                    }

                    // Skip orthogonal crossing check in singular triangles if disabled
                    if (inSingularTri && !ENABLE_ORTHOGONAL_CROSSING_IN_SINGULAR) continue;

                    // Only count (near) orthogonal crossings between separatrices.
                    if (!isOrthogonal(segVec, subP(sr.b, sr.a))) continue;

                    // Debounce repeated hits at nearly the same location.
                    auto itLast = lastCrossPoint.find(sr.sep_id);
                    if (itLast != lastCrossPoint.end()) {
                        if (normP(subP(ip, itLast->second)) < 0.25 * (avgEdgeLen > 0.0 ? avgEdgeLen : 1.0)) {
                            continue;
                        }
                    }
                    lastCrossPoint[sr.sep_id] = ip;

                    int &cnt = crossCount[sr.sep_id];
                    cnt += 1;
                    if (cnt >= 2) {
                        terminatedByIntersection = true;
                        hitSepId = sr.sep_id;
                        intersectionPoint = ip;
                        break;
                    }
                }

                if (terminatedByIntersection) {
                    TracePoint ipPt;
                    ipPt.face_id = currentPt.face_id;
                    ipPt.global_pos = intersectionPoint;
                    ipPt.barycentric = sanitizeBarycentric(globalToBarycentric(ipPt.face_id, intersectionPoint));

                    sep.path.push_back(ipPt);
                    if (currentTri >= 0 && currentTri < static_cast<int>(segIndex.size())) {
                        segIndex[currentTri].push_back(SegmentRef{sep.id, sep.path.size() - 2, currentPt.global_pos, ipPt.global_pos});
                    }

                    sep.active = false;
                    sep.termination_reason = TerminationReason::SELF_INTERSECTION;
                    break;
                }

                // --- Robust tangential merge detection (per tangential_crossing_robustness.md)
                // Checks for two cases:
                //   Case A: True intersection at small angle (nearly parallel/anti-parallel)
                //   Case B: No intersection but nearly coincident (small perpendicular distance)
                //
                // We select the *earliest* tangential event (smallest parameter t along our segment).
                // For anti-parallel: terminate here. For parallel: skip (same direction).

                bool tangentialMerged = false;
                int tangentialMergeTargetId = -1;
                Point tangentialMergePoint{0.0, 0.0};
                int tangentialMergeTri = -1;
                double bestTangentialT = 2.0; // >1 means no event found yet
                size_t tangentialMergeTargetSegIdx = 0; // segment index in target separatrix to clip at

                if (ENABLE_TANGENTIAL_CROSSING_CHECK) {

                // Helper: validate segment reference (ensure the other separatrix has this segment)
                auto segRefValid = [&](const SegmentRef &sr) -> bool {
                    if (sr.sep_id < 0 || sr.sep_id >= static_cast<int>(separatrices.size())) return false;
                    const auto &otherSep = separatrices[sr.sep_id];
                    if (sr.start_idx + 1 >= otherSep.path.size()) return false;
                    // Check that the recorded endpoints match the actual path
                    const Point &pathA = otherSep.path[sr.start_idx].global_pos;
                    const Point &pathB = otherSep.path[sr.start_idx + 1].global_pos;
                    double distA = normP(subP(pathA, sr.a));
                    double distB = normP(subP(pathB, sr.b));
                    double tol = 1e-10;
                    return (distA < tol && distB < tol);
                };

                Point ourDir = normalizeP(segVec);
                double ourSegLen = normP(segVec);
                if (ourSegLen > 0.0) {
                    // Lambda to check tangential merge candidates in a given triangle
                    auto checkTangentialInTri = [&](int triId) {
                        if (triId < 0 || triId >= static_cast<int>(segIndex.size())) return;
                        for (const auto &sr : segIndex[triId]) {
                            if (sr.sep_id == sep.id) continue;
                            
                            // Validate segment reference to prevent ghost intersections
                            if (!segRefValid(sr)) continue;

                            Point otherVec = subP(sr.b, sr.a);
                            double otherLen = normP(otherVec);
                            if (otherLen <= 0.0) continue;
                            Point otherDir = mulP(otherVec, 1.0 / otherLen);

                            double dp = dotP(ourDir, otherDir);
                            double absDp = std::abs(dp);

                            // Only consider if nearly parallel or anti-parallel
                            if (absDp <= TANGENTIAL_ABS_DOT_THRESHOLD) continue;

                            // ----- Case A: Actual intersection -----
                            double tI, sI;
                            Point ip;
                            if (segmentSegmentIntersectionParams(currentPt.global_pos, nextPt.global_pos,
                                                                  sr.a, sr.b, tI, sI, ip)) {
                                // Valid interior intersection
                                if (tI > 0.01 && tI < 0.99 && sI > 0.01 && sI < 0.99) {
                                    if (dp < 0.0 && tI < bestTangentialT) {
                                        // Anti-parallel intersection - record as tangential merge
                                        bestTangentialT = tI;
                                        tangentialMerged = true;
                                        tangentialMergeTargetId = sr.sep_id;
                                        tangentialMergePoint = ip;
                                        tangentialMergeTri = triId;
                                        tangentialMergeTargetSegIdx = sr.start_idx;
                                    }
                                    // If parallel (dp > 0), they're going the same direction - skip
                                }
                            }

                            // ----- Case B: No intersection but nearly coincident -----
                            // Check endpoint distances to the other segment
                            auto [distA, cpA, tA_param] = pointToSegmentDistance(currentPt.global_pos, sr.a, sr.b);
                            auto [distB, cpB, tB_param] = pointToSegmentDistance(nextPt.global_pos, sr.a, sr.b);

                            // Distance threshold based on merge distance (use tighter threshold)
                            double nearCoincidentDist = mergeDist * 0.5;

                            // Check if either endpoint is close to the other segment
                            bool aIsClose = (distA < nearCoincidentDist && tA_param > 0.01 && tA_param < 0.99);
                            bool bIsClose = (distB < nearCoincidentDist && tB_param > 0.01 && tB_param < 0.99);

                            if (aIsClose || bIsClose) {
                                // Check offset perpendicularity: offset from our segment to other segment
                                // should be mostly perpendicular to our direction
                                Point closerPt = aIsClose ? currentPt.global_pos : nextPt.global_pos;
                                Point closerCp = aIsClose ? cpA : cpB;
                                Point offset = subP(closerCp, closerPt);
                                double offsetLen = normP(offset);
                                
                                if (offsetLen > 1e-12) {
                                    Point offsetDir = mulP(offset, 1.0 / offsetLen);
                                    double along = std::abs(dotP(offsetDir, ourDir));
                                    
                                    // If offset is mostly perpendicular (along < threshold), it's tangential
                                    if (along < TANGENTIAL_OFFSET_ALONG_THRESHOLD) {
                                        // Compute parameter t along our segment for this event
                                        double eventT = aIsClose ? 0.0 : 1.0;
                                        
                                        // Only consider if anti-parallel and earlier than best
                                        if (dp < 0.0 && eventT < bestTangentialT) {
                                            bestTangentialT = eventT;
                                            tangentialMerged = true;
                                            tangentialMergeTargetId = sr.sep_id;
                                            tangentialMergePoint = closerCp;
                                            tangentialMergeTri = triId;
                                            tangentialMergeTargetSegIdx = sr.start_idx;
                                        }
                                    }
                                }
                            }

                            // Also check: our segment midpoint distance to other segment
                            Point midPt = mulP(addP(currentPt.global_pos, nextPt.global_pos), 0.5);
                            auto [distMid, cpMid, tMid_param] = pointToSegmentDistance(midPt, sr.a, sr.b);
                            if (distMid < nearCoincidentDist && tMid_param > 0.01 && tMid_param < 0.99) {
                                Point offset = subP(cpMid, midPt);
                                double offsetLen = normP(offset);
                                if (offsetLen > 1e-12) {
                                    Point offsetDir = mulP(offset, 1.0 / offsetLen);
                                    double along = std::abs(dotP(offsetDir, ourDir));
                                    if (along < TANGENTIAL_OFFSET_ALONG_THRESHOLD) {
                                        double eventT = 0.5;
                                        if (dp < 0.0 && eventT < bestTangentialT) {
                                            bestTangentialT = eventT;
                                            tangentialMerged = true;
                                            tangentialMergeTargetId = sr.sep_id;
                                            tangentialMergePoint = cpMid;
                                            tangentialMergeTri = triId;
                                            tangentialMergeTargetSegIdx = sr.start_idx;
                                        }
                                    }
                                }
                            }
                        }
                    };

                    // Check in current triangle and next triangle
                    checkTangentialInTri(currentTri);
                    if (nextPt.face_id != currentTri) {
                        checkTangentialInTri(nextPt.face_id);
                    }
                }
                } // end if (ENABLE_TANGENTIAL_CROSSING_CHECK)

                // Also maintain the edge crossing index for the original discrete check
                double edgeTol = std::max(1e-8, 1e-4 * avgEdgeLen);
                int nextTri = nextPt.face_id;
                if (nextTri != currentTri && nextTri >= 0) {
                    int exitEdge = pointOnTriangleEdge(mesh, currentTri, nextPt.global_pos, edgeTol);
                    auto entryIt = sepEntryState.find(sep.id);
                    if (exitEdge >= 0 && entryIt != sepEntryState.end()) {
                        auto [entryTri, entryEdge, entryPos] = entryIt->second;
                        if (entryTri == currentTri && entryEdge != exitEdge) {
                            edgeCrossingIndex[currentTri].push_back(EdgeCrossing{
                                sep.id, entryEdge, exitEdge, entryPos, nextPt.global_pos, sep.path.size() - 1
                            });
                        }
                    }
                    
                    int entryEdgeInNext = pointOnTriangleEdge(mesh, nextTri, nextPt.global_pos, edgeTol);
                    if (entryEdgeInNext >= 0) {
                        sepEntryState[sep.id] = std::make_tuple(nextTri, entryEdgeInNext, nextPt.global_pos);
                    } else {
                        sepEntryState.erase(sep.id);
                    }
                }

                if (tangentialMerged && tangentialMergeTargetId >= 0) {
                    // Terminate current separatrix at the tangential merge point
                    TracePoint ipPt;
                    ipPt.face_id = tangentialMergeTri;
                    ipPt.global_pos = tangentialMergePoint;
                    ipPt.barycentric = sanitizeBarycentric(globalToBarycentric(tangentialMergeTri, tangentialMergePoint));

                    sep.path.push_back(ipPt);
                    if (currentTri >= 0 && currentTri < static_cast<int>(segIndex.size())) {
                        segIndex[currentTri].push_back(SegmentRef{sep.id, sep.path.size() - 2, currentPt.global_pos, ipPt.global_pos});
                    }

                    sep.active = false;
                    sep.termination_reason = TerminationReason::MERGED_TANGENTIAL;
                    sep.merged_with_id = tangentialMergeTargetId;

                    // Also clip the target separatrix at the merge point
                    Separatrix &targetSep = separatrices[tangentialMergeTargetId];
                    if (targetSep.active || targetSep.termination_reason != TerminationReason::MERGED_TANGENTIAL) {
                        // Clip target's path: keep points up to tangentialMergeTargetSegIdx+1, then add merge point
                        size_t clipIdx = tangentialMergeTargetSegIdx + 1;
                        if (clipIdx < targetSep.path.size()) {
                            std::deque<TracePoint> clippedPath;
                            for (size_t k = 0; k <= clipIdx && k < targetSep.path.size(); ++k) {
                                clippedPath.push_back(targetSep.path[k]);
                            }
                            clippedPath.push_back(ipPt);
                            targetSep.path = clippedPath;
                        } else {
                            // Segment is at the end, just add the merge point
                            targetSep.path.push_back(ipPt);
                        }
                        targetSep.active = false;
                        targetSep.termination_reason = TerminationReason::MERGED_TANGENTIAL;
                        targetSep.merged_with_id = sep.id;
                    }

                    break;
                }

                // --- Stopping criterion: limit cycle detection
                // If the separatrix re-enters a triangle it has already visited (and left), stop.
                // Only check when transitioning to a different triangle.
                if (ENABLE_LIMIT_CYCLE_DETECTION) {
                    if (nextPt.face_id != currentTri) {
                        // We are entering a new triangle - check if we've been there before
                        if (visitedTriangles.find(nextPt.face_id) != visitedTriangles.end()) {
                            sep.active = false;
                            sep.termination_reason = TerminationReason::LIMIT_CYCLE;
                            break;
                        }
                        visitedTriangles.insert(nextPt.face_id);
                    }
                }

                // Accept the step.
                sep.path.push_back(nextPt);
                segIndex[currentTri].push_back(SegmentRef{sep.id, sep.path.size() - 2, currentPt.global_pos, nextPt.global_pos});

                Point lastSeg = subP(nextPt.global_pos, currentPt.global_pos);
                prevAngle = std::atan2(lastSeg[1], lastSeg[0]);

                if (stepRes.hitBoundary) {
                    sep.active = false;
                    sep.termination_reason = TerminationReason::HIT_BOUNDARY;
                    break;
                }

                currentPt = nextPt;
                currentTri = currentPt.face_id;
                stepCount++;

                if (!sep.active) break;
                if (stepCount >= MAX_STEPS) break;
            }

            if (stepCount >= MAX_STEPS && sep.active) {
                sep.active = false;
                sep.termination_reason = TerminationReason::MAX_STEPS_REACHED;
            }
        }
    }

    // --- Phase 2: Retrace separatrices that terminated with LIMIT_CYCLE or MAX_STEPS_REACHED ---
    // These are retraced with the "repeated orthogonal crossing" criterion:
    // Cut off when hitting the same OTHER separatrix twice.
    {
        std::vector<int> retraceIds;
        for (auto &sep : separatrices) {
            if (sep.termination_reason == TerminationReason::LIMIT_CYCLE ||
                sep.termination_reason == TerminationReason::MAX_STEPS_REACHED) {
                retraceIds.push_back(sep.id);
            }
        }

        if (!retraceIds.empty()) {
            std::cout << "Phase 2: Retracing " << retraceIds.size() 
                      << " separatrices with repeated-crossing criterion\n";
        }

        for (int sepId : retraceIds) {
            Separatrix &sep = separatrices[sepId];

            // Store initial path (first 2 points)
            if (sep.path.size() < 2) continue;
            TracePoint initPt0 = sep.path[0];
            TracePoint initPt1 = sep.path[1];

            // Remove old segments from segIndex
            for (size_t triIdx = 0; triIdx < segIndex.size(); ++triIdx) {
                auto &segs = segIndex[triIdx];
                segs.erase(std::remove_if(segs.begin(), segs.end(),
                    [sepId](const SegmentRef &sr) { return sr.sep_id == sepId; }),
                    segs.end());
            }

            // Reset separatrix
            sep.path.clear();
            sep.path.push_back(initPt0);
            sep.path.push_back(initPt1);
            sep.active = true;
            sep.termination_reason = TerminationReason::RUNNING;
            sep.visited_edges.clear();

            // Re-add initial segment to segIndex
            int initTri = initPt0.face_id;
            if (initTri >= 0 && initTri < static_cast<int>(segIndex.size())) {
                segIndex[initTri].push_back(SegmentRef{sep.id, 0, initPt0.global_pos, initPt1.global_pos});
            }

            // Track crossings with other separatrices
            std::unordered_map<int, int> crossCount;
            std::unordered_map<int, Point> lastCrossPoint;

            // Previous direction
            Point initDir = subP(initPt1.global_pos, initPt0.global_pos);
            if (normP(initDir) < 1e-30) initDir = Point{1.0, 0.0};
            double prevAngle = std::atan2(initDir[1], initDir[0]);

            int stepCount = 0;

            while (sep.active && stepCount < MAX_STEPS) {
                TracePoint currentPt = sep.path.back();
                int currentTri = currentPt.face_id;

                if (currentTri < 0 || currentTri >= static_cast<int>(mesh.triangles.size())) {
                    sep.active = false;
                    sep.termination_reason = TerminationReason::MAX_STEPS_REACHED;
                    break;
                }

                struct StepResult {
                    std::vector<TracePoint> pts;
                    bool hitBoundary = false;
                } stepRes;

                auto appendPointInTri = [&](int triId, const Point &pos) {
                    TracePoint tp;
                    tp.face_id = triId;
                    tp.global_pos = pos;
                    tp.barycentric = sanitizeBarycentric(globalToBarycentric(triId, pos));
                    stepRes.pts.push_back(tp);
                };

                // Regular advance (simplified - using Euler step for phase 2)
                {
                    double h = baseStep;
                    if (std::isfinite(triMinEdge[currentTri])) h = std::min(h, 0.25 * triMinEdge[currentTri]);

                    double theta = interpolateFieldDirection(currentTri, currentPt.barycentric, prevAngle);
                    Point dir{std::cos(theta), std::sin(theta)};

                    Point candidate = addP(currentPt.global_pos, mulP(dir, h));
                    std::array<double, 3> bCand = globalToBarycentric(currentTri, candidate);

                    if (baryInside(bCand, EPS_BARY)) {
                        appendPointInTri(currentTri, candidate);
                    } else {
                        // Find edge exit
                        double bestT = 2.0;
                        int bestEdge = -1;
                        Point bestHit{0.0, 0.0};

                        for (int e = 0; e < 3; ++e) {
                            int v0 = mesh.triangles[currentTri][e];
                            int v1 = mesh.triangles[currentTri][(e + 1) % 3];
                            Point p0{mesh.vertices[v0][0], mesh.vertices[v0][1]};
                            Point p1{mesh.vertices[v1][0], mesh.vertices[v1][1]};

                            double t, s;
                            Point ip;
                            if (segmentSegmentIntersectionParams(currentPt.global_pos, candidate, p0, p1, t, s, ip)) {
                                if (t > EPS_INTERSECT_T && t < bestT && s >= -EPS_INTERSECT_T && s <= 1.0 + EPS_INTERSECT_T) {
                                    bestT = t;
                                    bestEdge = e;
                                    bestHit = ip;
                                }
                            }
                        }

                        if (bestEdge >= 0) {
                            int neighbor = mesh.triangleAdjacency[currentTri][bestEdge];
                            if (neighbor < 0) {
                                appendPointInTri(currentTri, bestHit);
                                stepRes.hitBoundary = true;
                            } else {
                                appendPointInTri(neighbor, bestHit);
                            }
                        } else {
                            appendPointInTri(currentTri, candidate);
                        }
                    }
                }

                if (stepRes.pts.empty()) {
                    sep.active = false;
                    sep.termination_reason = TerminationReason::MAX_STEPS_REACHED;
                    break;
                }

                // Process each step point
                for (size_t k = 0; k < stepRes.pts.size(); ++k) {
                    TracePoint nextPt = stepRes.pts[k];

                    Point segVec = subP(nextPt.global_pos, currentPt.global_pos);
                    double segLen = normP(segVec);
                    if (segLen < minAdvance) {
                        sep.active = false;
                        sep.termination_reason = TerminationReason::MAX_STEPS_REACHED;
                        break;
                    }

                    // Check for self-intersection OR repeated crossing of another separatrix
                    bool terminatedByCrossing = false;
                    bool isSelfIntersect = false;
                    Point intersectionPoint{0.0, 0.0};

                    const auto &candidates = segIndex[currentTri];
                    for (const auto &sr : candidates) {
                        double tI, sI;
                        Point ip;
                        if (!segmentSegmentIntersectionParams(currentPt.global_pos, nextPt.global_pos, sr.a, sr.b, tI, sI, ip)) continue;

                        if (tI <= EPS_INTERSECT_T || tI >= 1.0 - EPS_INTERSECT_T) continue;
                        if (sI <= EPS_INTERSECT_T || sI >= 1.0 - EPS_INTERSECT_T) continue;

                        if (sr.sep_id == sep.id) {
                            // Self-intersection: skip most recent segments to avoid false positives
                            if (sep.path.size() >= 2 && sr.start_idx >= sep.path.size() - 2) continue;
                            
                            // Definitive self-intersection - this IS a limit cycle
                            terminatedByCrossing = true;
                            isSelfIntersect = true;
                            intersectionPoint = ip;
                            break;
                        } else {
                            // Phase 2: count ANY crossing (not just orthogonal) since we're looking for limit cycles

                            // Debounce repeated hits at same location
                            auto itLast = lastCrossPoint.find(sr.sep_id);
                            if (itLast != lastCrossPoint.end()) {
                                if (normP(subP(ip, itLast->second)) < 0.25 * (avgEdgeLen > 0.0 ? avgEdgeLen : 1.0)) {
                                    continue;
                                }
                            }
                            lastCrossPoint[sr.sep_id] = ip;

                            int &cnt = crossCount[sr.sep_id];
                            cnt += 1;
                            if (cnt >= 2) {
                                terminatedByCrossing = true;
                                intersectionPoint = ip;
                                break;
                            }
                        }
                    }

                    if (terminatedByCrossing) {
                        TracePoint ipPt;
                        ipPt.face_id = currentTri;
                        ipPt.global_pos = intersectionPoint;
                        ipPt.barycentric = sanitizeBarycentric(globalToBarycentric(currentTri, intersectionPoint));

                        sep.path.push_back(ipPt);
                        segIndex[currentTri].push_back(SegmentRef{sep.id, sep.path.size() - 2, currentPt.global_pos, ipPt.global_pos});

                        sep.active = false;
                        sep.termination_reason = isSelfIntersect ? TerminationReason::SELF_INTERSECTION : TerminationReason::REPEATED_CROSSING;
                        break;
                    }

                    // Accept the step
                    sep.path.push_back(nextPt);
                    segIndex[currentTri].push_back(SegmentRef{sep.id, sep.path.size() - 2, currentPt.global_pos, nextPt.global_pos});

                    Point lastSeg = subP(nextPt.global_pos, currentPt.global_pos);
                    prevAngle = std::atan2(lastSeg[1], lastSeg[0]);

                    if (stepRes.hitBoundary) {
                        sep.active = false;
                        sep.termination_reason = TerminationReason::HIT_BOUNDARY;
                        break;
                    }

                    currentPt = nextPt;
                    currentTri = currentPt.face_id;
                    stepCount++;

                    if (!sep.active) break;
                    if (stepCount >= MAX_STEPS) break;
                }

                if (stepCount >= MAX_STEPS && sep.active) {
                    sep.active = false;
                    sep.termination_reason = TerminationReason::MAX_STEPS_REACHED;
                }
            }
        }

        // Report phase 2 results
        int repeatedCrossingCount = 0;
        for (const auto &sep : separatrices) {
            if (sep.termination_reason == TerminationReason::REPEATED_CROSSING) {
                repeatedCrossingCount++;
            }
        }
        if (!retraceIds.empty()) {
            std::cout << "Phase 2 complete: " << repeatedCrossingCount 
                      << " terminated by repeated crossing\n";
        }
    }

    // --- Post-processing pass: detect and clip any remaining tangential crossings ---
    // This catches cases where both separatrices passed through before either could detect the other.
    // Uses the same robust algorithm as the main loop: Case A (intersection) and Case B (near-coincident).
    
    if (ENABLE_TANGENTIAL_CROSSING_CHECK) {
    // Build a map from separatrix ID to its segments for efficient clipping
    std::unordered_map<int, std::vector<std::pair<size_t, size_t>>> sepToSegments; // sep_id -> [(tri_id, seg_idx_in_tri), ...]
    for (size_t triIdx = 0; triIdx < segIndex.size(); ++triIdx) {
        for (size_t segIdx = 0; segIdx < segIndex[triIdx].size(); ++segIdx) {
            sepToSegments[segIndex[triIdx][segIdx].sep_id].emplace_back(triIdx, segIdx);
        }
    }
    
    // Structure for tangential crossings (both Case A and Case B)
    struct TangentialCrossing {
        int sep_a, sep_b;
        int tri_id;
        Point merge_point;
        double t_a, t_b; // parameters along each segment (t_b may be approximate for Case B)
        size_t seg_idx_a, seg_idx_b; // which segment in the path
        double dot_product;
        bool is_case_b; // true if near-coincident (Case B), false if intersection (Case A)
    };
    std::vector<TangentialCrossing> tangentialCrossings;
    
    // Debug counters
    int caseAcount = 0, caseBcount = 0;
    
    for (size_t triIdx = 0; triIdx < segIndex.size(); ++triIdx) {
        const auto &segs = segIndex[triIdx];
        for (size_t i = 0; i < segs.size(); ++i) {
            for (size_t j = i + 1; j < segs.size(); ++j) {
                const auto &segA = segs[i];
                const auto &segB = segs[j];
                
                if (segA.sep_id == segB.sep_id) continue;
                
                // Skip if either separatrix was already merged with the other
                if (separatrices[segA.sep_id].merged_with_id == segB.sep_id ||
                    separatrices[segB.sep_id].merged_with_id == segA.sep_id) {
                    continue;
                }
                
                // Skip if either was already marked as merged tangential
                if (separatrices[segA.sep_id].termination_reason == TerminationReason::MERGED_TANGENTIAL ||
                    separatrices[segB.sep_id].termination_reason == TerminationReason::MERGED_TANGENTIAL) {
                    continue;
                }
                
                // Compute directions
                Point vecA = subP(segA.b, segA.a);
                Point vecB = subP(segB.b, segB.a);
                double lenA = normP(vecA);
                double lenB = normP(vecB);
                if (lenA <= 0.0 || lenB <= 0.0) continue;
                
                Point dirA = mulP(vecA, 1.0 / lenA);
                Point dirB = mulP(vecB, 1.0 / lenB);
                double dp = dotP(dirA, dirB);
                double absDp = std::abs(dp);
                
                // Only consider if nearly parallel or anti-parallel
                if (absDp <= TANGENTIAL_ABS_DOT_THRESHOLD) continue;
                
                // Only consider if anti-parallel (dp < 0)
                if (dp >= 0.0) continue;
                
                // ----- Case A: Actual intersection -----
                double tI, sI;
                Point ip;
                if (segmentSegmentIntersectionParams(segA.a, segA.b, segB.a, segB.b, tI, sI, ip)) {
                    // Valid interior intersection
                    if (tI > 0.01 && tI < 0.99 && sI > 0.01 && sI < 0.99) {
                        tangentialCrossings.push_back({
                            segA.sep_id, segB.sep_id,
                            static_cast<int>(triIdx), ip,
                            tI, sI,
                            segA.start_idx, segB.start_idx,
                            dp, false
                        });
                        caseAcount++;
                        continue; // Found intersection, no need to check Case B
                    }
                }
                
                // ----- Case B: No intersection but nearly coincident -----
                // Check distances from endpoints of A to segment B and vice versa
                auto [distA0, cpA0, tA0] = pointToSegmentDistance(segA.a, segB.a, segB.b);
                auto [distA1, cpA1, tA1] = pointToSegmentDistance(segA.b, segB.a, segB.b);
                auto [distB0, cpB0, tB0] = pointToSegmentDistance(segB.a, segA.a, segA.b);
                auto [distB1, cpB1, tB1] = pointToSegmentDistance(segB.b, segA.a, segA.b);
                
                // Use merge distance as threshold
                double nearCoincidentDist = mergeDist;
                
                // Find the closest approach
                double minDist = std::min({distA0, distA1, distB0, distB1});
                
                if (minDist < nearCoincidentDist) {
                    // Determine which point is closest and compute merge point
                    Point closerPt, closerCp;
                    double tOnA, tOnB_approx;
                    bool fromA; // true if the closest point is on segment A
                    
                    if (distA0 <= distA1 && distA0 <= distB0 && distA0 <= distB1 && tA0 > 0.01 && tA0 < 0.99) {
                        closerPt = segA.a; closerCp = cpA0; tOnA = 0.0; tOnB_approx = tA0; fromA = true;
                    } else if (distA1 <= distA0 && distA1 <= distB0 && distA1 <= distB1 && tA1 > 0.01 && tA1 < 0.99) {
                        closerPt = segA.b; closerCp = cpA1; tOnA = 1.0; tOnB_approx = tA1; fromA = true;
                    } else if (distB0 <= distA0 && distB0 <= distA1 && distB0 <= distB1 && tB0 > 0.01 && tB0 < 0.99) {
                        closerPt = segB.a; closerCp = cpB0; tOnA = tB0; tOnB_approx = 0.0; fromA = false;
                    } else if (distB1 <= distA0 && distB1 <= distA1 && distB1 <= distB0 && tB1 > 0.01 && tB1 < 0.99) {
                        closerPt = segB.b; closerCp = cpB1; tOnA = tB1; tOnB_approx = 1.0; fromA = false;
                    } else {
                        continue; // All close points are at endpoints, skip
                    }
                    
                    // Check offset perpendicularity
                    Point offset = subP(closerCp, closerPt);
                    double offsetLen = normP(offset);
                    if (offsetLen > 1e-12) {
                        Point offsetDir = mulP(offset, 1.0 / offsetLen);
                        Point refDir = fromA ? dirA : dirB;
                        double along = std::abs(dotP(offsetDir, refDir));
                        
                        if (along < TANGENTIAL_OFFSET_ALONG_THRESHOLD) {
                            // Use the midpoint between closest points as merge point
                            Point mergePt = mulP(addP(closerPt, closerCp), 0.5);
                            tangentialCrossings.push_back({
                                segA.sep_id, segB.sep_id,
                                static_cast<int>(triIdx), mergePt,
                                tOnA, tOnB_approx,
                                segA.start_idx, segB.start_idx,
                                dp, true
                            });
                            caseBcount++;
                        }
                    }
                }
            }
        }
    }
    
    // Print debug info
    if (caseAcount > 0 || caseBcount > 0) {
        std::cout << "Post-processing: found " << caseAcount << " Case A (intersection) and " 
                  << caseBcount << " Case B (near-coincident) tangential crossings\n";
    }
    
    // Process tangential crossings: merge the separatrices at the crossing point
    // For anti-parallel crossings, we connect sepA's path (up to crossing) with sepB's path (reversed from crossing)
    int mergedCount = 0;
    std::unordered_set<int> alreadyMerged; // Track separatrices that have been merged into another
    
    // Sort crossings by earliest occurrence (smallest segment index) to process in order
    std::sort(tangentialCrossings.begin(), tangentialCrossings.end(), 
              [](const TangentialCrossing &a, const TangentialCrossing &b) {
                  return std::min(a.seg_idx_a, a.seg_idx_b) < std::min(b.seg_idx_a, b.seg_idx_b);
              });
    
    for (const auto &cross : tangentialCrossings) {
        // Skip if either separatrix was already merged into another
        if (alreadyMerged.count(cross.sep_a) || alreadyMerged.count(cross.sep_b)) continue;
        
        Separatrix &sepA = separatrices[cross.sep_a];
        Separatrix &sepB = separatrices[cross.sep_b];
        
        // Skip if already processed (merged)
        if (sepA.termination_reason == TerminationReason::MERGED_TANGENTIAL &&
            sepA.merged_with_id == cross.sep_b) continue;
        if (sepB.termination_reason == TerminationReason::MERGED_TANGENTIAL &&
            sepB.merged_with_id == cross.sep_a) continue;
        
        // For anti-parallel tangential crossings:
        // - Keep sepA's path from start up to and including the merge point
        // - Append sepB's path reversed from the merge point back to sepB's start
        // This creates a single continuous path from sepA's origin through the merge to sepB's origin
        
        // Determine clip indices for both separatrices
        size_t clipIdxA = cross.seg_idx_a + 1; // Keep points 0..clipIdxA, then add merge point
        size_t clipIdxB = cross.seg_idx_b + 1; // Will reverse from this point back to 0
        
        // Safety checks
        if (clipIdxA >= sepA.path.size()) clipIdxA = sepA.path.size() - 1;
        if (clipIdxB >= sepB.path.size()) clipIdxB = sepB.path.size() - 1;
        
        // Build the merged path for sepA
        std::deque<TracePoint> mergedPath;
        
        // Add sepA's path up to the crossing segment endpoint
        for (size_t k = 0; k <= clipIdxA && k < sepA.path.size(); ++k) {
            mergedPath.push_back(sepA.path[k]);
        }
        
        // Add the merge point
        TracePoint mergePt;
        mergePt.face_id = cross.tri_id;
        mergePt.global_pos = cross.merge_point;
        mergePt.barycentric = sanitizeBarycentric(globalToBarycentric(cross.tri_id, cross.merge_point));
        mergedPath.push_back(mergePt);
        
        // Add sepB's path reversed from the crossing back to sepB's origin
        // Start from clipIdxB and go backwards to 0
        for (size_t k = clipIdxB; k > 0; --k) {
            mergedPath.push_back(sepB.path[k - 1]);
        }
        
        // Update sepA with the merged path
        sepA.path = mergedPath;
        sepA.termination_reason = TerminationReason::MERGED_TANGENTIAL;
        sepA.merged_with_id = cross.sep_b;
        
        // Mark sepB as absorbed (its path is now part of sepA)
        sepB.active = false;
        sepB.termination_reason = TerminationReason::MERGED_TANGENTIAL;
        sepB.merged_with_id = cross.sep_a;
        alreadyMerged.insert(cross.sep_b);
        
        mergedCount++;
    }
    
    if (mergedCount > 0) {
        std::cout << "Post-processing: merged " << mergedCount << " separatrix pairs at tangential crossings\n";
    }
    } // end if (ENABLE_TANGENTIAL_CROSSING_CHECK)
}
