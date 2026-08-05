#include "tracing/FieldTracer.hxx"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

// Wrap to [0, 2pi).
inline double wrap_2pi(double a) {
    a = std::fmod(a, 2.0 * M_PI);
    if (a < 0.0) a += 2.0 * M_PI;
    return a;
}

// The representative of `candidate` modulo pi/2 nearest `ref`. In every use
// below the two arguments are the same angle read two ways, so the rounding is
// decided by a margin of pi/4 rather than by a tolerance.
inline double liftNear(double ref, double candidate) {
    const double k = std::round((ref - candidate) / M_PI_2);
    return candidate + k * M_PI_2;
}

inline Point dirOf(double angle) { return Point{std::cos(angle), std::sin(angle)}; }

} // namespace

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
FieldTracer::FieldTracer(std::shared_ptr<CrossField> cf, bool useExactCenters)
    : crossField(std::move(cf)) {
    mesh = crossField->mesh.get();

    const int nV = static_cast<int>(mesh->vertices.size());
    vertexTheta.resize(nV);
    for (int v = 0; v < nV; ++v) vertexTheta[v] = std::arg(crossField->u_k[v]) * 0.25;

    double total = 0.0;
    for (const auto &e : mesh->edges) total += normP(mesh->vertices[e[1]] - mesh->vertices[e[0]]);
    if (!mesh->edges.empty()) avgEdge = total / static_cast<double>(mesh->edges.size());
    vertexSnap = 1e-9 * avgEdge;

    buildSingularities(useExactCenters);
}

// ---------------------------------------------------------------------------
// buildSingularities()  --  centres and ports
//
// The centre is the point of the triangle where the interpolated representation
// vector vanishes, which is where the cross has no orientation left; the
// barycentre stands in when that solve puts it outside the triangle, which
// happens when the three vertex values are nearly collinear in the plane and
// the 2x2 system is ill-conditioned.
//
// The ports are fitted rather than read off one corner. Sec. 3.2.2 writes the
// model as f(z) = e^{i d theta / 4} with theta measured from a separatrix, so
// in those coordinates both the position and the direction are measured from
// that separatrix, and putting it back in the mesh's own frame -- where the
// separatrix points at some absolute angle beta -- rotates the direction as
// well as the point:
//
//     phi(omega) = beta + d (omega - beta) / 4    (mod pi/2)
//
// which rearranges to C := phi - d omega / 4 = beta (4 - d) / 4, so each corner
// of the singular triangle gives an estimate of C and
//
//     beta = 4 C / (4 - d).
//
// Dropping the rotation of the direction -- reading the model as
// phi = d(omega - beta)/4 -- happens to give the right answer when beta is
// zero and a wrong one everywhere else, which is a comfortable way to be wrong
// and worth stating explicitly. Two things fall out of the correct form: the
// pi/2 by which phi is ambiguous moves beta by exactly one port spacing, so the
// port set does not depend on which branch each corner is read on; and the
// error in C is multiplied by only 4/(4-d), which is 4/3 at a 3-port
// singularity rather than the 4 it would be otherwise. Averaging the three
// corners rather than trusting one is what makes that error small to begin
// with.
// ---------------------------------------------------------------------------
void FieldTracer::buildSingularities(bool useExactCenters) {
    singularityOf.assign(mesh->triangles.size(), -1);
    singularities.clear();

    for (const auto &[triIdx, cfIndex] : crossField->singularTriangles) {
        if (triIdx < 0 || triIdx >= static_cast<int>(mesh->triangles.size())) continue;

        Singularity s;
        s.triangleIndex = triIdx;
        s.singularityIndex = cfIndex;
        s.d = static_cast<int>(std::lround(cfIndex * 4.0));

        // 4 - d has to be a usable number of sectors for the conformal model:
        // d = 4 collapses it and d > 4 turns the exponent negative. Neither
        // happens for the simple singularities a smooth field produces, but a
        // ruined field should not take the tracer down with it.
        const int nPorts = 4 - s.d;
        if (nPorts < 2 || nPorts > 8) continue;

        const Triangle &tri = mesh->triangles[triIdx];
        const Point &p0 = mesh->vertices[tri[0]];
        const Point &p1 = mesh->vertices[tri[1]];
        const Point &p2 = mesh->vertices[tri[2]];

        bool haveCenter = false;
        if (useExactCenters) {
            const std::complex<double> u0 = crossField->u_k[tri[0]];
            const std::complex<double> u1 = crossField->u_k[tri[1]];
            const std::complex<double> u2 = crossField->u_k[tri[2]];

            Eigen::Matrix2d A;
            A << std::real(u0) - std::real(u2), std::real(u1) - std::real(u2),
                 std::imag(u0) - std::imag(u2), std::imag(u1) - std::imag(u2);
            Eigen::Vector2d b;
            b << -std::real(u2), -std::imag(u2);

            if (std::fabs(A.determinant()) > 1e-14) {
                const Eigen::Vector2d ab = A.colPivHouseholderQr().solve(b);
                const std::array<double, 3> bc{{ab[0], ab[1], 1.0 - ab[0] - ab[1]}};
                const double lo = std::min({bc[0], bc[1], bc[2]});
                if (lo > 0.02) {  // comfortably inside, not just barely
                    s.barycentric = bc;
                    s.coordinates = p0 * bc[0] + p1 * bc[1] + p2 * bc[2];
                    haveCenter = true;
                }
            }
        }
        if (!haveCenter) {
            s.barycentric = {{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}};
            s.coordinates = (p0 + p1 + p2) / 3.0;
        }

        // C = theta - d omega / 4 at each corner, lifted onto a common branch.
        std::array<double, 3> cEst{{0.0, 0.0, 0.0}};
        for (int i = 0; i < 3; ++i) {
            const Point r = mesh->vertices[tri[i]] - s.coordinates;
            const double omega = std::atan2(r[1], r[0]);
            const double raw = vertexTheta[tri[i]] - 0.25 * s.d * omega;
            cEst[i] = (i == 0) ? raw : liftNear(cEst[0], raw);
        }
        const double cBar = (cEst[0] + cEst[1] + cEst[2]) / 3.0;
        s.portResidual = std::max({std::fabs(cEst[0] - cBar), std::fabs(cEst[1] - cBar),
                                   std::fabs(cEst[2] - cBar)});

        const double beta = 4.0 * cBar / static_cast<double>(nPorts);  // nPorts == 4 - d
        const double sector = 2.0 * M_PI / static_cast<double>(nPorts);
        // Sorted counter-clockwise out of [0, sector), which is what the sector
        // lookup during a sweep assumes.
        const double beta0 = wrap_2pi(beta) - sector * std::floor(wrap_2pi(beta) / sector);
        s.portAngles.resize(nPorts);
        for (int k = 0; k < nPorts; ++k) s.portAngles[k] = beta0 + k * sector;
        s.portSeparatrixIds.assign(nPorts, -1);

        singularityOf[triIdx] = static_cast<int>(singularities.size());
        singularities.push_back(std::move(s));
    }
}

// ---------------------------------------------------------------------------
// Field readings
// ---------------------------------------------------------------------------
double FieldTracer::vertexAngle(int v, double ref) const {
    return liftNear(ref, vertexTheta[v]);
}

double FieldTracer::edgeAngle(int e, double t, double ref) const {
    const int a = mesh->edges[e][0];
    const int b = mesh->edges[e][1];
    const double ta = vertexTheta[a];
    const double tb = liftNear(ta, vertexTheta[b]);  // the principal matching of this edge
    return liftNear(ref, ta * (1.0 - t) + tb * t);
}

void FieldTracer::triangleAngles(int f, const std::array<double, 3> &at, double ref,
                                 std::array<double, 3> &theta) const {
    const Triangle &tri = mesh->triangles[f];
    theta[0] = vertexTheta[tri[0]];
    theta[1] = liftNear(theta[0], vertexTheta[tri[1]]);
    theta[2] = liftNear(theta[0], vertexTheta[tri[2]]);

    // Shift the whole triple -- which leaves the interpolant affine -- until its
    // value at the entry point is the branch the walk arrived on.
    const double here = theta[0] * at[0] + theta[1] * at[1] + theta[2] * at[2];
    const double shift = liftNear(ref, here) - here;
    theta[0] += shift;
    theta[1] += shift;
    theta[2] += shift;
}

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------
std::array<double, 3> FieldTracer::barycentric(int f, const Point &p) const {
    const Triangle &tri = mesh->triangles[f];
    const Point &v0 = mesh->vertices[tri[0]];
    const Point &v1 = mesh->vertices[tri[1]];
    const Point &v2 = mesh->vertices[tri[2]];

    const Point e1 = v1 - v0;
    const Point e2 = v2 - v0;
    const Point w = p - v0;
    const double den = cross2(e1, e2);
    if (std::fabs(den) < 1e-30) return {{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}};

    const double l1 = cross2(w, e2) / den;
    const double l2 = cross2(e1, w) / den;
    return {{1.0 - l1 - l2, l1, l2}};
}

bool FieldTracer::exitOfRay(int f, const Point &from, const Point &dir, int excludeEdge,
                            Point &hit, int &edge, double &along) const {
    const Triangle &tri = mesh->triangles[f];
    double best = std::numeric_limits<double>::max();
    edge = -1;

    for (int i = 0; i < 3; ++i) {
        if (i == excludeEdge) continue;
        const Point &A = mesh->vertices[tri[i]];
        const Point &B = mesh->vertices[tri[(i + 1) % 3]];
        const Point e = B - A;

        const double den = cross2(dir, e);
        if (std::fabs(den) < 1e-18) continue;         // parallel to this edge
        const Point w = A - from;
        const double s = cross2(w, e) / den;
        const double u = cross2(w, dir) / den;
        if (s <= 1e-13 * avgEdge) continue;           // behind, or where we started
        if (u < -1e-9 || u > 1.0 + 1e-9) continue;    // misses the segment

        if (s < best) {
            best = s;
            edge = i;
            along = std::min(1.0, std::max(0.0, u));
            hit = from + dir * s;
        }
    }
    return edge >= 0;
}

Walker FieldTracer::startAt(int f, const Point &p, double dir) const {
    Walker w;
    w.tri = f;
    w.pos = p;
    w.entryEdge = -1;
    w.atVertex = -1;
    w.dir = wrap_pi(dir);
    w.crossDir = w.dir;
    return w;
}

// ---------------------------------------------------------------------------
// crossEdge()  --  hand the walk to the triangle on the far side
//
// The parameter along the edge is re-expressed in the neighbour's orientation
// of it rather than the hit point being re-projected, so the two triangles
// agree on the crossing point to the last bit and the polyline stays connected.
// ---------------------------------------------------------------------------
FieldTracer::Status FieldTracer::crossEdge(Walker &w, int f, int edge, double along,
                                           const Point &hit, double outTheta,
                                           double chordDir) const {
    const Triangle &tri = mesh->triangles[f];
    const int va = tri[edge];
    const int vb = tri[(edge + 1) % 3];
    const int ge = mesh->triangleEdges[f][edge];

    // A crossing that lands on a vertex is moved onto the vertex outright.
    // Leaving it a hair to one side is what produces the degenerate slivers of
    // a step: a segment shorter than the rounding error of its own endpoints,
    // whose direction is noise.
    const double edgeLen = normP(mesh->vertices[vb] - mesh->vertices[va]);
    if (along * edgeLen < vertexSnap || (1.0 - along) * edgeLen < vertexSnap) {
        w.atVertex = (along * edgeLen < vertexSnap) ? va : vb;
        w.pos = mesh->vertices[w.atVertex];
        w.dir = wrap_pi(vertexAngle(w.atVertex, outTheta));
        w.crossDir = w.dir;
        w.entryEdge = -1;
        w.tri = f;
        return leaveVertex(w);
    }

    const int nb = mesh->triangleAdjacency[f][edge];
    if (nb < 0) {
        w.pos = hit;
        w.tri = f;
        w.entryEdge = edge;
        w.atVertex = -1;
        w.dir = wrap_pi(outTheta);
        w.crossDir = wrap_pi(chordDir);
        return Status::Boundary;
    }

    // Which local edge of the neighbour this is, and which way round.
    const Triangle &nt = mesh->triangles[nb];
    int ne = -1;
    bool flipped = false;
    for (int i = 0; i < 3; ++i) {
        if (nt[i] == va && nt[(i + 1) % 3] == vb) { ne = i; flipped = false; break; }
        if (nt[i] == vb && nt[(i + 1) % 3] == va) { ne = i; flipped = true;  break; }
    }
    if (ne < 0) return Status::Stuck;  // adjacency and connectivity disagree

    const double t = flipped ? (1.0 - along) : along;
    const Point &A = mesh->vertices[nt[ne]];
    const Point &B = mesh->vertices[nt[(ne + 1) % 3]];

    w.tri = nb;
    w.pos = A * (1.0 - t) + B * t;
    w.entryEdge = ne;
    w.atVertex = -1;
    // The field on the shared edge, read from the edge's own two endpoints, is
    // the same number on both banks; the lift picks up the multiple of pi/2 the
    // two triangles' interpolants differ by, which is what keeps the branch
    // from drifting.
    const double tOnEdge = (mesh->edges[ge][0] == nt[ne]) ? t : 1.0 - t;
    w.dir = wrap_pi(edgeAngle(ge, tOnEdge, outTheta));
    w.crossDir = wrap_pi(chordDir);
    return Status::Ok;
}

// ---------------------------------------------------------------------------
// leaveVertex()  --  push a walk sitting on a vertex into the right triangle
//
// A streamline through a vertex has no edge to be handed across, so the
// triangle it continues in is the one whose wedge at that vertex holds its
// direction. Taking the wedge the direction sits furthest inside, rather than
// the first that contains it, keeps the choice away from the wedge borders,
// where a direction that grazes an edge would otherwise be assigned by
// rounding.
// ---------------------------------------------------------------------------
FieldTracer::Status FieldTracer::leaveVertex(Walker &w) const {
    const int v = w.atVertex;
    if (v < 0) return Status::Stuck;

    const auto range = mesh->vertexTriangles.trianglesForVertex(v);
    const Point d = dirOf(w.dir);

    int best = -1;
    double bestMargin = -std::numeric_limits<double>::max();
    for (const int *it = range.first; it != range.second; ++it) {
        const int f = *it;
        const Triangle &tri = mesh->triangles[f];
        int i0 = -1;
        for (int i = 0; i < 3; ++i) if (tri[i] == v) i0 = i;
        if (i0 < 0) continue;

        const Point a = mesh->vertices[tri[(i0 + 1) % 3]] - mesh->vertices[v];
        const Point b = mesh->vertices[tri[(i0 + 2) % 3]] - mesh->vertices[v];
        const double na = normP(a), nb = normP(b);
        if (na < 1e-18 || nb < 1e-18) continue;

        // Sines of the angles from each side of the wedge to d, signed so that
        // both are positive exactly when d is inside it.
        const double sgn = (cross2(a, b) > 0.0) ? 1.0 : -1.0;
        const double margin = std::min(sgn * cross2(a, d) / na, sgn * cross2(d, b) / nb);
        if (margin > bestMargin) { bestMargin = margin; best = f; }
    }

    if (best < 0) return Status::Stuck;
    w.tri = best;
    if (bestMargin < -1e-9) {
        // Out of every wedge. At an interior vertex that cannot happen, so this
        // is a boundary vertex and the streamline is leaving the mesh; the
        // vertex id is left set for the caller to record where.
        return Status::Boundary;
    }

    w.entryEdge = -1;
    w.pos = mesh->vertices[v];
    w.atVertex = -1;
    w.crossDir = w.dir;
    return Status::Ok;
}

// ---------------------------------------------------------------------------
// advance()  --  one triangle
//
// Regular triangles are integrated with the trapezoid rule iterated to its own
// fixed point. Because the lifted angle is affine over the triangle, its mean
// along the chord equals its value at the chord's midpoint exactly, so what is
// being solved is the midpoint rule and the direction that comes out is the
// average direction along the whole segment rather than the direction at either
// end of it. Two or three iterations settle it.
// ---------------------------------------------------------------------------
FieldTracer::Status FieldTracer::advance(Walker &w, std::vector<TracePoint> &path,
                                         CutInfo *cut) const {
    if (w.tri < 0) return Status::Stuck;

    if (w.atVertex >= 0) {
        const Status s = leaveVertex(w);
        if (s != Status::Ok) return s;
    }

    const int sing = singularityOf[w.tri];
    if (sing >= 0) return sweepSingular(w, sing, path, cut);

    const int f = w.tri;
    const std::array<double, 3> entryBary = barycentric(f, w.pos);

    std::array<double, 3> th;
    triangleAngles(f, entryBary, w.dir, th);
    const double theta0 = th[0] * entryBary[0] + th[1] * entryBary[1] + th[2] * entryBary[2];

    Point hit{0.0, 0.0};
    int edge = -1;
    double along = 0.0;
    double theta = theta0;

    if (!exitOfRay(f, w.pos, dirOf(theta0), w.entryEdge, hit, edge, along)) {
        // The field at the entry point points back out of the edge just
        // crossed: the streamline is grazing that edge, and which side of it
        // the field is on at the crossing point is decided by rounding. The
        // direction that carried the walk over the edge is the one thing known
        // to point inwards, so the grazing stretch is traced along it.
        if (!exitOfRay(f, w.pos, dirOf(w.crossDir), w.entryEdge, hit, edge, along))
            return Status::Stuck;
        theta = w.crossDir;
    }

    for (int iter = 0; iter < 4; ++iter) {
        const std::array<double, 3> hb = barycentric(f, hit);
        const double thetaEnd = th[0] * hb[0] + th[1] * hb[1] + th[2] * hb[2];
        const double thetaMid = 0.5 * (theta0 + thetaEnd);
        if (std::fabs(thetaMid - theta) < 1e-14) break;

        Point h2{0.0, 0.0};
        int e2 = -1;
        double a2 = 0.0;
        if (!exitOfRay(f, w.pos, dirOf(thetaMid), w.entryEdge, h2, e2, a2)) break;
        theta = thetaMid;
        hit = h2;
        edge = e2;
        along = a2;
    }

    const std::array<double, 3> hb = barycentric(f, hit);
    const double thetaExit = th[0] * hb[0] + th[1] * hb[1] + th[2] * hb[2];

    const Status st = crossEdge(w, f, edge, along, hit, thetaExit, theta);

    TracePoint tp;
    // The far side's reading of the crossing point, or the vertex it snapped
    // to; either way the point the next segment starts from, so the polyline
    // has no gap at a triangle boundary.
    tp.global_pos = (st == Status::Stuck) ? hit : w.pos;
    tp.face_id = f;
    tp.theta = wrap_pi(thetaExit);
    path.push_back(tp);

    return st;
}

// ---------------------------------------------------------------------------
// sweepSingular()  --  the hyperbola of Prop. 1, inside one singular triangle
//
// With the centre at the origin, theta measured from a separatrix and
// M = (4 - d)/8, the map w = z^M sends the pair of sectors bounded by the
// separatrices s0, s1, s2 onto the first quadrant and the streamline crossing
// s1 orthogonally onto a branch of x y = A. That branch is asymptotic to s0 and
// to s2, and its closest approach to the singularity is exactly on s1, at
// phi = pi/4.
//
// The other streamline through the same point -- the one crossing s0 -- is the
// same curve seen from the other end: measuring theta clockwise from s1 instead
// of counter-clockwise from s0 turns it into x y = A again. So one code path
// serves both, with a frame (reference separatrix, orientation) picked from the
// incoming direction. In exact arithmetic the branch index is even in exactly
// one of the two frames, which is what selects between them.
//
// Samples sit on a fan of rays out of the centre at fixed multiples of
// sector / sweepSamplesPerSector. Because the fan is shared by every streamline
// crossing the triangle and the hyperbolas are nested, two of them keep their
// order on every ray and so cannot cross tangentially -- the guarantee the
// whole construction exists for.
// ---------------------------------------------------------------------------
FieldTracer::Status FieldTracer::sweepSingular(Walker &w, int sing,
                                               std::vector<TracePoint> &path,
                                               CutInfo *cut) const {
    const Singularity &s = singularities[sing];
    const int f = w.tri;
    const int nPorts = s.numPorts();
    const double sector = 2.0 * M_PI / static_cast<double>(nPorts);
    const double M = static_cast<double>(nPorts) / 8.0;  // (4 - d)/8

    const Point rel = w.pos - s.coordinates;
    const double r0 = normP(rel);

    // The sector holding the entry point.
    int lower = 0;
    double thetaQ = 2.0 * M_PI;
    if (r0 > 1e-12 * avgEdge) {
        const double omega = std::atan2(rel[1], rel[0]);
        for (int k = 0; k < nPorts; ++k) {
            const double t = wrap_2pi(omega - s.portAngles[k]);
            if (t < thetaQ) { thetaQ = t; lower = k; }
        }
    } else {
        thetaQ = 0.0;  // on the singularity: launching a separatrix
    }

    // The frame. Exactly one of the two gives an even branch index; when
    // rounding makes that a close call, the one closer to even wins rather than
    // the walk failing.
    double refAngle = s.portAngles[lower];
    double sigma = 1.0;
    double theta = thetaQ;
    int m = 0;
    {
        double bestErr = std::numeric_limits<double>::max();
        for (int c = 0; c < 2; ++c) {
            const double sg = (c == 0) ? 1.0 : -1.0;
            const double ref = (c == 0) ? s.portAngles[lower] : s.portAngles[lower] + sector;
            const double th = (c == 0) ? thetaQ : sector - thetaQ;
            const double x = (sg * (w.dir - ref) - 0.25 * s.d * th) / M_PI_2;
            const double mEven = 2.0 * std::round(x * 0.5);
            const double err = std::fabs(x - mEven);
            if (err < bestErr) {
                bestErr = err;
                refAngle = ref;
                sigma = sg;
                theta = th;
                m = ((static_cast<int>(mEven) % 4) + 4) % 4;
            }
        }
    }

    const double phiQ = M * theta;
    const int entryEdge = w.entryEdge;

    // On a separatrix the hyperbola degenerates to its own asymptote and the
    // streamline is the straight ray. That covers both a walk sitting on a port
    // ray and a walk launched from the singularity itself, which is on all of
    // them at once and would otherwise be handed A = 0 and every sample back at
    // the centre.
    if (r0 <= 1e-9 * avgEdge || phiQ <= 1e-9 || phiQ >= M_PI_2 - 1e-9) {
        Point hit{0.0, 0.0};
        int edge = -1;
        double along = 0.0;
        if (!exitOfRay(f, w.pos, dirOf(w.dir), entryEdge, hit, edge, along)) return Status::Stuck;
        const double outTheta = w.dir;  // radial, and constant along a separatrix

        const Status st = crossEdge(w, f, edge, along, hit, outTheta, outTheta);
        TracePoint tp;
        tp.global_pos = (st == Status::Stuck) ? hit : w.pos;
        tp.face_id = f;
        tp.theta = wrap_pi(outTheta);
        path.push_back(tp);
        return st;
    }

    // phi decreases when the streamline runs back out along the reference
    // separatrix, increases when it sweeps round towards the next one.
    const double dStep = (m == 2) ? 1.0 : -1.0;
    const double rhoQ = std::pow(r0, M);
    const double A = rhoQ * rhoQ * std::sin(phiQ) * std::cos(phiQ);

    auto pointAt = [&](double th) {
        const double phi = M * th;
        const double denom = std::max(std::sin(phi) * std::cos(phi), 1e-300);
        const double rho = std::sqrt(A / denom);
        const double r = std::pow(rho, 1.0 / M);
        const double abs = refAngle + sigma * th;
        return s.coordinates + Point{r * std::cos(abs), r * std::sin(abs)};
    };
    // Frame angle of an arbitrary point, kept in [0, sector] -- a point that
    // rounds just past the reference separatrix must not come back as 2 pi.
    auto frameThetaOf = [&](const Point &q) {
        const Point rr = q - s.coordinates;
        double t = wrap_2pi(sigma * (std::atan2(rr[1], rr[0]) - refAngle));
        if (t > M_PI) t -= 2.0 * M_PI;
        return std::min(std::max(t, 0.0), sector);
    };

    const double dTheta = sector / static_cast<double>(sweepSamplesPerSector);
    // The point of a shared fan is that samples land on it, so the first step
    // goes to the next ray of the fan and not a fixed distance on.
    double gridIdx = std::floor(theta / dTheta);
    if (dStep > 0.0) gridIdx += 1.0;

    Point prev = w.pos;
    double prevTheta = theta;
    int excludeEdge = entryEdge;
    const int maxSamples = 2 * sweepSamplesPerSector + 8;

    // Running out along the reference separatrix, the curve is asymptotic to it
    // and the fan of rays never reaches it, so the last stretch before the
    // triangle wall is finished as a straight ray in the tangent direction. It
    // is a straight line there to the accuracy of everything else in the step,
    // and without it a streamline entering just beside a port has nowhere to go.
    auto finishStraight = [&](const Point &from, double fromTheta) {
        const double out = refAngle + sigma * (0.25 * s.d * fromTheta + m * M_PI_2);
        Point hit{0.0, 0.0};
        int edge = -1;
        double along = 0.0;
        const int ex = (from == w.pos) ? entryEdge : -1;
        if (!exitOfRay(f, from, dirOf(out), ex, hit, edge, along) &&
            !exitOfRay(f, from, dirOf(w.crossDir), ex, hit, edge, along))
            return Status::Stuck;
        const Status st = crossEdge(w, f, edge, along, hit, out, out);
        TracePoint tp;
        tp.global_pos = (st == Status::Stuck) ? hit : w.pos;
        tp.face_id = f;
        tp.theta = wrap_pi(out);
        path.push_back(tp);
        return st;
    };

    for (int i = 0; i < maxSamples; ++i) {
        double next = gridIdx * dTheta;
        bool atPort = false;
        // The sweep can only ever meet a separatrix at theta = sector: a
        // hyperbola never reaches its own asymptote, so the reference one at
        // theta = 0 is out of reach. Land on it exactly rather than stepping
        // across it.
        if (dStep > 0.0 && next >= sector - 1e-15) { next = sector; atPort = true; }
        if (dStep < 0.0 && next <= 1e-15) return finishStraight(prev, prevTheta);
        if (std::fabs(next - prevTheta) < 1e-15) { gridIdx += dStep; continue; }

        const Point cur = pointAt(next);

        // Does the step leave the triangle? The previous sample is inside it
        // (or on the entry edge), so a ray from it meets exactly one edge.
        const Point seg = cur - prev;
        const double segLen = normP(seg);
        if (segLen > 1e-15 * avgEdge) {
            Point hit{0.0, 0.0};
            int edge = -1;
            double along = 0.0;
            if (exitOfRay(f, prev, seg / segLen, excludeEdge, hit, edge, along) &&
                normP(hit - prev) <= segLen) {
                // The field direction at the exit is the model's, which is
                // smooth and is what the branch on the far side is matched
                // against.
                const double outTheta =
                    refAngle + sigma * (0.25 * s.d * frameThetaOf(hit) + m * M_PI_2);

                const Status st = crossEdge(w, f, edge, along, hit, outTheta,
                                            computeAngle(seg));
                TracePoint tp;
                tp.global_pos = (st == Status::Stuck) ? hit : w.pos;
                tp.face_id = f;
                tp.theta = wrap_pi(outTheta);
                path.push_back(tp);
                return st;
            }
        }

        TracePoint tp;
        tp.global_pos = cur;
        tp.face_id = f;
        tp.theta = wrap_pi(refAngle + sigma * (0.25 * s.d * next + m * M_PI_2));
        path.push_back(tp);

        if (atPort) {
            // Sec. 3.3, third stopping condition: a streamline crossing a
            // separatrix of a singularity orthogonally inside that
            // singularity's own triangle is cut there. That crossing is the
            // closest it ever comes to the singularity, and stopping turns what
            // would otherwise be two streamlines running side by side into a
            // T-junction on the separatrix.
            if (cut) {
                cut->singularity = sing;
                cut->port = (lower + (sigma > 0.0 ? 1 : 0)) % nPorts;
            }
            w.pos = cur;
            w.tri = f;
            w.entryEdge = -1;
            w.atVertex = -1;
            w.dir = tp.theta;
            return Status::Cut;
        }

        prev = cur;
        prevTheta = next;
        gridIdx += dStep;
        excludeEdge = -1;  // only the first step starts on the entry edge
    }

    return finishStraight(prev, prevTheta);
}
