#include "Domains.hxx"

#include <algorithm>
#include <cmath>
#include <map>
#include <random>
#include <stdexcept>

#include "triangle/TriangleMesher.hpp"

namespace paper {
namespace {

using triangle_wrapper::TriangleMesher2D;

// A PSLG under construction, with one property that matters more than anything
// else here: **a point is added once**.
//
// Every multi-material domain in this file has interface chains that end on the
// exterior boundary and meet each other at a junction, so the same coordinate is
// reached from two or three directions. Handing Triangle two vertices at the
// same place is not an error it refuses -- it silently drops the duplicate and
// leaves the segments that referenced it pointing somewhere unintended -- so the
// interface would come out unconstrained on one side of the junction and the
// material ids would flood past it. `point()` keys on the coordinate and returns
// the index that is already there, which makes the junction one vertex by
// construction rather than by tolerance.
struct PSLG {
    TriangleMesher2D::MeshInput in;

    int point(const Point &p) {
        const long long kx = std::llround(p[0] / kQuantum);
        const long long ky = std::llround(p[1] / kQuantum);
        auto it = index.find({kx, ky});
        if (it != index.end()) return it->second;
        const int id = static_cast<int>(in.vertlist.size());
        in.vertlist.push_back({p[0], p[1]});
        index.emplace(std::make_pair(kx, ky), id);
        return id;
    }

    // The segments of the straight run a -> b, cut into pieces of about `h`.
    // The endpoints go through point(), so a chain landing on the boundary or
    // on another chain shares that vertex.
    std::vector<std::array<int, 2>> segmentsAlong(const Point &a, const Point &b, double h) {
        const int n = std::max(1, static_cast<int>(std::ceil(normP(b - a) / h)));
        std::vector<std::array<int, 2>> segs;
        int prev = point(a);
        for (int i = 1; i <= n; ++i) {
            const double t = static_cast<double>(i) / n;
            const int cur = point({a[0] * (1.0 - t) + b[0] * t, a[1] * (1.0 - t) + b[1] * t});
            segs.push_back({prev, cur});
            prev = cur;
        }
        return segs;
    }

    // The closed ring through `corners`, each side cut at `h`.
    std::vector<std::array<int, 2>> ring(const std::vector<Point> &corners, double h) {
        std::vector<std::array<int, 2>> segs;
        for (std::size_t k = 0; k < corners.size(); ++k) {
            const auto side = segmentsAlong(corners[k], corners[(k + 1) % corners.size()], h);
            segs.insert(segs.end(), side.begin(), side.end());
        }
        return segs;
    }

    // The closed ring of `n` points on a circle.
    std::vector<std::array<int, 2>> circle(const Point &c, double r, int n) {
        std::vector<Point> pts;
        for (int i = 0; i < n; ++i) {
            const double a = 2.0 * M_PI * i / n;
            pts.push_back({c[0] + r * std::cos(a), c[1] + r * std::sin(a)});
        }
        std::vector<std::array<int, 2>> segs;
        for (int i = 0; i < n; ++i) segs.push_back({point(pts[i]), point(pts[(i + 1) % n])});
        return segs;
    }

    // An open polyline through the given points, each run cut at `h`.
    std::vector<std::array<int, 2>> polyline(const std::vector<Point> &pts, double h) {
        std::vector<std::array<int, 2>> segs;
        for (std::size_t k = 0; k + 1 < pts.size(); ++k) {
            const auto run = segmentsAlong(pts[k], pts[k + 1], h);
            segs.insert(segs.end(), run.begin(), run.end());
        }
        return segs;
    }

    void exterior(std::vector<std::array<int, 2>> segs) {
        in.segment_loops.push_back(std::move(segs));
        in.type.push_back(0);
        in.loop_seed.push_back({0.0, 0.0});
        in.region_id.push_back(0);
    }
    void hole(std::vector<std::array<int, 2>> segs, const Point &seed) {
        in.segment_loops.push_back(std::move(segs));
        in.type.push_back(1);
        in.loop_seed.push_back({seed[0], seed[1]});
        in.region_id.push_back(0);
    }
    // A material: its seed point and id, plus whatever interface segments are
    // listed under it. Which segments those are carries no meaning -- Triangle
    // floods the id from the seed until it meets *any* segment -- and they are
    // distributed among the materials only because the mesher insists every
    // loop have at least three of them.
    void region(std::vector<std::array<int, 2>> segs, const Point &seed, int id) {
        in.segment_loops.push_back(std::move(segs));
        in.type.push_back(2);
        in.loop_seed.push_back({seed[0], seed[1]});
        in.region_id.push_back(id);
    }

private:
    static constexpr double kQuantum = 1e-10;
    std::map<std::pair<long long, long long>, int> index;
};

std::shared_ptr<Mesh> build(const PSLG &g, double minAngle = 28.0, bool justDelaunay = false) {
    TriangleMesher2D::Options opts;
    opts.min_angle_degrees = minAngle;
    // The constrained Delaunay triangulation of exactly the given points: no
    // quality refinement, no area constraint, no Steiner points. The caller is
    // then in charge of where every vertex is, which is what gradedDiskDomain
    // needs and what no area constraint can give it.
    opts.just_delaunay = justDelaunay;
    TriangleMesher2D mesher(opts);
    TriangleMesher2D::MeshOutput out = mesher.triangulate(g.in);

    std::vector<Point> verts(out.verts.size());
    for (std::size_t i = 0; i < out.verts.size(); ++i) verts[i] = {out.verts[i][0], out.verts[i][1]};

    std::vector<Triangle> tris(out.triangles.size());
    for (std::size_t i = 0; i < out.triangles.size(); ++i)
        tris[i] = {out.triangles[i][0], out.triangles[i][1], out.triangles[i][2]};

    // Triangle's region attribute is 0 wherever no seed reached, and Mesh reads
    // an empty id vector as "one material". Both fold into the same rule: if any
    // region was flooded, shift the ids so the smallest present is 1.
    std::vector<int> mats;
    const bool anyRegion = std::any_of(out.tri_regions.begin(), out.tri_regions.end(),
                                       [](int r) { return r != 0; });
    if (anyRegion) {
        const int lo = *std::min_element(out.tri_regions.begin(), out.tri_regions.end());
        mats.resize(out.tri_regions.size());
        for (std::size_t i = 0; i < out.tri_regions.size(); ++i) mats[i] = out.tri_regions[i] - lo + 1;
    }
    return std::make_shared<Mesh>(verts, tris, mats);
}

int circlePoints(double r, double h) {
    return std::max(24, static_cast<int>(std::round(2.0 * M_PI * r / h)));
}

} // namespace

// ---------------------------------------------------------------------------
// The canonical domains
// ---------------------------------------------------------------------------
std::shared_ptr<Mesh> squareDomain(double h) {
    PSLG g;
    g.in.h = h;
    g.exterior(g.ring({{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}, h));
    return build(g);
}

std::shared_ptr<Mesh> diskDomain(double h) {
    PSLG g;
    g.in.h = h;
    g.exterior(g.circle({0.0, 0.0}, 1.0, circlePoints(1.0, h)));
    return build(g);
}

std::shared_ptr<Mesh> annulusDomain(double rIn, double rOut, double h) {
    PSLG g;
    g.in.h = h;
    g.exterior(g.circle({0.0, 0.0}, rOut, circlePoints(rOut, h)));
    g.hole(g.circle({0.0, 0.0}, rIn, circlePoints(rIn, h)), {0.0, 0.0});
    return build(g);
}

std::shared_ptr<Mesh> lShapeDomain(double h) {
    PSLG g;
    g.in.h = h;
    g.exterior(g.ring({{0.0, 0.0}, {1.0, 0.0}, {1.0, 0.5},
                       {0.5, 0.5}, {0.5, 1.0}, {0.0, 1.0}}, h));
    return build(g);
}

std::shared_ptr<Mesh> wedgeDomain(double sweepDegrees, double radius, double h) {
    const double sweep = sweepDegrees * M_PI / 180.0;
    const int narc = std::max(8, static_cast<int>(std::round(radius * sweep / h)));

    // The ring, in order: the apex, out along the first radius, round the arc,
    // and back down the second. `ring` cuts each run at h and shares endpoints,
    // so the apex is one vertex however many runs end there.
    std::vector<Point> corners;
    corners.push_back({0.0, 0.0});
    for (int i = 0; i <= narc; ++i) {
        const double t = sweep * i / narc;
        corners.push_back({radius * std::cos(t), radius * std::sin(t)});
    }

    PSLG g;
    g.in.h = h;
    g.exterior(g.ring(corners, h));
    // A 45-degree apex is a sharper angle than Triangle's quality bound can
    // honour, and it responds by refusing to terminate rather than by failing.
    return build(g, sweepDegrees < 60.0 ? 20.0 : 25.0);
}

std::shared_ptr<Mesh> regularPolygonDomain(int n, double h) {
    if (n < 3) throw std::runtime_error("regularPolygonDomain: n must be at least 3");
    std::vector<Point> corners;
    for (int i = 0; i < n; ++i) {
        const double a = 2.0 * M_PI * i / n + 0.5 * M_PI / n; // no side axis-aligned
        corners.push_back({std::cos(a), std::sin(a)});
    }
    PSLG g;
    g.in.h = h;
    g.exterior(g.ring(corners, h));
    return build(g);
}

std::shared_ptr<Mesh> gradedDiskDomain(double radius, int rings, double growth) {
    if (radius <= 0.0 || rings < 3 || growth < 1.0)
        throw std::runtime_error("gradedDiskDomain: need radius > 0, rings >= 3, growth >= 1");

    // The ring radii, from a geometric sequence of radial spacings, scaled so
    // the outermost ring is exactly the boundary circle. Fixing the ring count
    // rather than the first spacing is what makes the sweep a comparison: every
    // grading gets the same radial resolution and only the *ratio* between
    // neighbouring elements changes, so a difference between two levels is the
    // grading and not the mesh size.
    std::vector<double> r{0.0}, dr;
    double step = 1.0;
    for (int i = 0; i < rings; ++i) {
        dr.push_back(step);
        r.push_back(r.back() + step);
        step *= growth;
    }
    const double scale = radius / r.back();
    for (double &x : r) x *= scale;
    for (double &x : dr) x *= scale;

    // The points: the centre, then each ring, with the angular count chosen so
    // that a ring's arc spacing is about its radial spacing and the triangles
    // stay near-isotropic. The grading is then radial and the anisotropy stays
    // bounded, so what the sweep varies is the *area ratio across an edge* and
    // not the shape of the elements.
    PSLG g;
    g.in.h = dr.front();
    g.point({0.0, 0.0});
    std::vector<Point> outer;
    for (std::size_t i = 1; i < r.size(); ++i) {
        const int n = std::max(6, static_cast<int>(std::llround(2.0 * M_PI * r[i] / dr[i - 1])));
        for (int k = 0; k < n; ++k) {
            const double a = 2.0 * M_PI * k / n;
            const Point p{r[i] * std::cos(a), r[i] * std::sin(a)};
            g.point(p);
            if (i + 1 == r.size()) outer.push_back(p);
        }
    }

    // Only the outer ring carries segments: it is the boundary of the domain.
    // Every other point is an interior PSLG vertex, which the constrained
    // Delaunay triangulation keeps exactly where it was put.
    std::vector<std::array<int, 2>> ring;
    for (std::size_t i = 0; i < outer.size(); ++i)
        ring.push_back({g.point(outer[i]), g.point(outer[(i + 1) % outer.size()])});
    g.exterior(ring);

    return build(g, 0.0, /*justDelaunay=*/true);
}

std::vector<Domain> canonicalDomains(double h) {
    std::vector<Domain> d;
    auto add = [&](const std::string &name, std::shared_ptr<Mesh> m,
                   bool known, int idx4, const std::string &note) {
        Domain dom;
        dom.name = name;
        dom.mesh = std::move(m);
        dom.knownIndex = known;
        dom.expectedInteriorIndex4 = idx4;
        dom.note = note;
        d.push_back(std::move(dom));
    };

    // The four with a closed-form answer. Each is 4 chi minus the quarter turns
    // the corners take, and on these four every corner angle is a multiple of
    // pi/2 so that count is exact rather than a rounding.
    add("square",  squareDomain(h),            true, 0,
        "4 right-angle corners take 4 chi; the smoothest field is constant");
    add("disk",    diskDomain(h),              true, 4,
        "no corner takes anything: 4 x +1/4 inside, on the diagonals");
    add("annulus", annulusDomain(0.4, 1.0, h), true, 0,
        "chi = 0 and neither rim has a corner");
    add("lshape",  lShapeDomain(h),            true, 0,
        "5 convex corners and 1 reflex: 5 - 1 = 4 quarter turns");

    // The rest have at least one corner whose turn is not a multiple of pi/2,
    // so what the boundary takes is a rounding and the interior total follows
    // from it rather than from the shape. E2 reports the identity instead of
    // asserting a number.
    add("wedge45",  wedgeDomain(45.0, 1.0, h),  false, 0, "apex angle pi/4");
    add("wedge120", wedgeDomain(120.0, 1.0, h), false, 0, "apex angle 2pi/3");
    add("wedge270", wedgeDomain(270.0, 1.0, h), false, 0, "reflex apex, 3pi/2");
    add("pentagon", regularPolygonDomain(5, h), false, 0, "five 108-degree corners");
    return d;
}

// ---------------------------------------------------------------------------
std::shared_ptr<Mesh> jitter(const Mesh &src, double frac, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> unit(-1.0, 1.0);

    // The local edge length at each vertex: the mean length of the edges meeting
    // there, so the displacement is a fraction of the *mesh* and not of the
    // model, and a graded mesh is jittered evenly.
    const int nV = static_cast<int>(src.vertices.size());
    std::vector<double> len(nV, 0.0);
    std::vector<int> cnt(nV, 0);
    for (const auto &e : src.edges) {
        const double l = normP(src.vertices[e[1]] - src.vertices[e[0]]);
        len[e[0]] += l; ++cnt[e[0]];
        len[e[1]] += l; ++cnt[e[1]];
    }

    std::vector<Point> verts = src.vertices;
    for (int v = 0; v < nV; ++v) {
        if (src.isBoundaryVertex[v] || cnt[v] == 0) continue;
        const double hv = len[v] / cnt[v];
        verts[v][0] += frac * hv * unit(rng);
        verts[v][1] += frac * hv * unit(rng);
    }
    // The boundary does not move, so the domain and every expectation about it
    // is unchanged and only the triangulation of it differs. Interior vertices
    // on a material interface do move, which changes the interface's polyline
    // -- E2 jitters single-material domains only, and says so.
    return std::make_shared<Mesh>(verts, src.triangles, src.triangleMatId);
}

// ---------------------------------------------------------------------------
// E5's junction domains
// ---------------------------------------------------------------------------
std::shared_ptr<Mesh> tJunctionDomain(double h) {
    PSLG g;
    g.in.h = h;
    // The interface landing points are corners of the exterior ring, so the
    // chains end on vertices the boundary already has.
    g.exterior(g.ring({{0.0, 0.0}, {0.5, 0.0}, {1.0, 0.0}, {1.0, 0.5},
                       {1.0, 1.0}, {0.0, 1.0}, {0.0, 0.5}}, h));
    g.region(g.polyline({{0.0, 0.5}, {0.5, 0.5}}, h), {0.5, 0.75}, 1);
    g.region(g.polyline({{0.5, 0.5}, {1.0, 0.5}}, h), {0.25, 0.25}, 2);
    g.region(g.polyline({{0.5, 0.5}, {0.5, 0.0}}, h), {0.75, 0.25}, 3);
    return build(g);
}

std::shared_ptr<Mesh> quadruplePointDomain(double h) {
    PSLG g;
    g.in.h = h;
    g.exterior(g.ring({{0.0, 0.0}, {0.5, 0.0}, {1.0, 0.0}, {1.0, 0.5},
                       {1.0, 1.0}, {0.5, 1.0}, {0.0, 1.0}, {0.0, 0.5}}, h));
    g.region(g.polyline({{0.0, 0.5}, {0.5, 0.5}}, h), {0.25, 0.25}, 1);
    g.region(g.polyline({{0.5, 0.5}, {1.0, 0.5}}, h), {0.75, 0.25}, 2);
    g.region(g.polyline({{0.5, 0.0}, {0.5, 0.5}}, h), {0.75, 0.75}, 3);
    g.region(g.polyline({{0.5, 0.5}, {0.5, 1.0}}, h), {0.25, 0.75}, 4);
    return build(g);
}

std::shared_ptr<Mesh> thinLayerDomain(double h) {
    // Four elements thick: thin enough that an irregular vertex inside the layer
    // is a real cost to the time step, wide enough that a layout has somewhere
    // to put one.
    const double half = 2.0 * h;
    const double lo = 0.5 - half, hi = 0.5 + half;

    PSLG g;
    g.in.h = h;
    g.exterior(g.ring({{0.0, 0.0}, {2.0, 0.0}, {2.0, lo}, {2.0, hi}, {2.0, 1.0},
                       {0.0, 1.0}, {0.0, hi}, {0.0, lo}}, h));
    // Three materials, three loops. The lower interface is one, the upper is
    // split in two so that the third material has segments of its own to be
    // listed with; no segment appears twice.
    g.region(g.polyline({{0.0, lo}, {2.0, lo}}, h), {1.0, 0.5 * lo}, 1);
    g.region(g.polyline({{0.0, hi}, {1.0, hi}}, h), {1.0, 0.5}, 2);
    g.region(g.polyline({{1.0, hi}, {2.0, hi}}, h), {1.0, 0.5 * (hi + 1.0)}, 3);
    return build(g);
}

std::shared_ptr<Mesh> embeddedInclusionDomain(double h) {
    const double r = 0.28;
    PSLG g;
    g.in.h = h;
    g.exterior(g.ring({{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}, h));

    // The circle, split into two arcs so that the inclusion and the matrix each
    // have a loop of their own to hang a seed on.
    const int n = circlePoints(r, h);
    auto all = g.circle({0.5, 0.5}, r, n);
    std::vector<std::array<int, 2>> lower(all.begin(), all.begin() + n / 2);
    std::vector<std::array<int, 2>> upper(all.begin() + n / 2, all.end());
    g.region(lower, {0.5, 0.5}, 1);   // the inclusion
    g.region(upper, {0.05, 0.05}, 2); // the matrix
    return build(g);
}

std::shared_ptr<Mesh> obliqueJunctionDomain(double h, double sectorDegrees) {
    PSLG g;
    g.in.h = h;

    // The three rays, at 90 degrees and then each sector further round; the
    // third sector is whatever is left of the turn.
    const double sector[3] = {sectorDegrees, sectorDegrees, 360.0 - 2.0 * sectorDegrees};
    double rayDeg[3];
    rayDeg[0] = 90.0;
    rayDeg[1] = rayDeg[0] + sector[0];
    rayDeg[2] = rayDeg[1] + sector[1];

    // The rays land on the rim, so put a rim vertex exactly where each one does
    // by building the circle out of three arcs that start at those points.
    const int n = circlePoints(1.0, h);
    std::vector<Point> rim;
    for (int k = 0; k < 3; ++k) {
        const int per = std::max(6, static_cast<int>(std::round(n * sector[k] / 360.0)));
        for (int i = 0; i < per; ++i) {
            const double a = (rayDeg[k] + sector[k] * static_cast<double>(i) / per) * M_PI / 180.0;
            rim.push_back({std::cos(a), std::sin(a)});
        }
    }
    std::vector<std::array<int, 2>> ring;
    for (std::size_t i = 0; i < rim.size(); ++i)
        ring.push_back({g.point(rim[i]), g.point(rim[(i + 1) % rim.size()])});
    g.exterior(ring);

    // At 120 degrees: three rays at 90, 210 and 330 degrees, every sector 120
    // degrees, so no sector is a multiple of a right angle and no one cross is
    // tangent to all three interfaces at once.
    for (int k = 0; k < 3; ++k) {
        const double a = rayDeg[k] * M_PI / 180.0;
        const double s = (rayDeg[k] + 0.5 * sector[k]) * M_PI / 180.0;
        g.region(g.polyline({{0.0, 0.0}, {std::cos(a), std::sin(a)}}, h),
                 {0.5 * std::cos(s), 0.5 * std::sin(s)}, k + 1);
    }
    return build(g);
}

// ---------------------------------------------------------------------------
// The mechanism set
// ---------------------------------------------------------------------------
std::shared_ptr<Mesh> combDomain(int teeth, double h) {
    if (teeth < 1) throw std::runtime_error("combDomain: need at least one tooth");
    // A 1 x 1 block with `teeth` slots cut down from the top edge, each one
    // slot wide with a tooth of the same width beside it. Every corner is 90 or
    // 270 degrees, which is the point.
    const int n = 2 * teeth + 1;          // teeth and slots alternating
    const double w = 1.0 / n;             // width of each
    const double d = 0.55;                // slot depth

    std::vector<Point> ring;
    ring.push_back({0.0, 0.0});
    ring.push_back({1.0, 0.0});
    ring.push_back({1.0, 1.0});
    // Walk back along the top, cutting a slot over every even-indexed band.
    for (int i = n - 1; i >= 0; --i) {
        const double xr = (i + 1) * w, xl = i * w;
        if (i % 2 == 1) {                 // a slot
            ring.push_back({xr, 1.0});
            ring.push_back({xr, 1.0 - d});
            ring.push_back({xl, 1.0 - d});
            ring.push_back({xl, 1.0});
        }
    }
    ring.push_back({0.0, 1.0});

    PSLG g;
    g.in.h = h;
    g.exterior(g.ring(ring, h));
    return build(g);
}

std::shared_ptr<Mesh> starDomain(int points, double tipDegrees, double h) {
    if (points < 3) throw std::runtime_error("starDomain: need at least three points");
    const double valley = 360.0 - 360.0 / points - tipDegrees;
    if (tipDegrees <= 5.0 || valley <= 5.0 || valley >= 355.0)
        throw std::runtime_error("starDomain: that tip angle leaves no valley");

    // Outer radius 1; the inner radius that gives the requested tip angle. The
    // tip is at (1,0) and its neighbours at radius r, angle +-pi/points, so
    // tan(tip/2) = r sin(pi/n) / (1 - r cos(pi/n)).
    const double a = M_PI / points;
    const double t = std::tan(0.5 * tipDegrees * M_PI / 180.0);
    const double r = t / (std::sin(a) + t * std::cos(a));
    if (!(r > 0.0) || r >= 1.0) throw std::runtime_error("starDomain: no inner radius for that tip angle");

    std::vector<Point> ring;
    for (int i = 0; i < points; ++i) {
        const double th = 2.0 * a * i;
        ring.push_back({std::cos(th), std::sin(th)});
        ring.push_back({r * std::cos(th + a), r * std::sin(th + a)});
    }

    PSLG g;
    g.in.h = h;
    g.exterior(g.ring(ring, h));
    // A sharp tip is sharper than the quality bound can honour; Triangle
    // responds by not terminating rather than by failing.
    return build(g, tipDegrees < 60.0 ? 20.0 : 25.0);
}

std::shared_ptr<Mesh> laminateDomain(int layers, double tiltDegrees, double h) {
    if (layers < 2) throw std::runtime_error("laminateDomain: need at least two layers");
    if (tiltDegrees < 20.0 || tiltDegrees > 90.0)
        throw std::runtime_error("laminateDomain: tilt outside 20..90 degrees");

    // The domain is the *unit square*, so every corner of dS is a right angle
    // and the corner rounding a vertex field does is exact there. The only
    // oblique thing in the domain is where an interface lands on dS, which is
    // what isolates the junction mechanism from the boundary one: a laminate at
    // 90 degrees is a control in both, and tilting only moves the junctions.
    const double tilt = tiltDegrees * M_PI / 180.0;
    const bool vertical = std::fabs(tiltDegrees - 90.0) < 1e-9;
    const double slope = vertical ? 0.0 : std::tan(tilt);   // dy/dx

    // Interface k leaves the bottom edge at x = k/layers and runs up-right at
    // `tilt`, until it leaves the square through the top or the right side.
    auto foot = [&](int k) { return static_cast<double>(k) / layers; };
    auto exitPoint = [&](double x0) {
        if (vertical) return Point{x0, 1.0};
        const double xTop = x0 + 1.0 / slope;
        if (xTop <= 1.0) return Point{xTop, 1.0};
        return Point{1.0, (1.0 - x0) * slope};
    };

    // The square's own corners, plus every interface foot, so that a foot is a
    // vertex of the boundary rather than a point the mesher may or may not put
    // one at, and plus every exit point for the same reason.
    std::vector<Point> ring;
    ring.push_back({0.0, 0.0});
    for (int k = 1; k < layers; ++k) ring.push_back({foot(k), 0.0});
    ring.push_back({1.0, 0.0});
    {   // up the right side, through any exits on it, bottom to top
        std::vector<double> ys;
        for (int k = 1; k < layers; ++k) {
            const Point e = exitPoint(foot(k));
            if (std::fabs(e[0] - 1.0) < 1e-12) ys.push_back(e[1]);
        }
        std::sort(ys.begin(), ys.end());
        for (double y : ys) ring.push_back({1.0, y});
    }
    ring.push_back({1.0, 1.0});
    {   // back along the top, through any exits on it, right to left
        std::vector<double> xs;
        for (int k = 1; k < layers; ++k) {
            const Point e = exitPoint(foot(k));
            if (std::fabs(e[1] - 1.0) < 1e-12) xs.push_back(e[0]);
        }
        std::sort(xs.begin(), xs.end(), std::greater<double>());
        for (double x : xs) ring.push_back({x, 1.0});
    }
    ring.push_back({0.0, 1.0});

    PSLG g;
    g.in.h = h;
    g.exterior(g.ring(ring, h));

    // `layers` bands need `layers` region loops but there are only layers - 1
    // interfaces, so the first interface is split in two and its halves carry
    // two of them. Which segments a region loop is listed with carries no
    // meaning: Triangle floods the id from the seed until it meets any segment.
    auto bandSeed = [&](int k) {   // band k lies right of interface k
        const double xl = (k == 0) ? 0.0 : foot(k);
        const double xr = (k + 1 <= layers - 1) ? foot(k + 1) : 1.0;
        return Point{0.5 * (xl + xr), 0.02};
    };
    {
        const Point a{foot(1), 0.0}, b = exitPoint(foot(1));
        const Point mid{0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])};
        g.region(g.polyline({a, mid}, h), bandSeed(0), 1);
        g.region(g.polyline({mid, b}, h), bandSeed(1), 2);
    }
    for (int k = 2; k < layers; ++k)
        g.region(g.polyline({{foot(k), 0.0}, exitPoint(foot(k))}, h), bandSeed(k), k + 1);
    return build(g);
}

namespace {

// The circumcenter of three points, or false if they are collinear.
bool circumcenterOf3(const Point &a, const Point &b, const Point &c, Point &out) {
    const double bx = b[0] - a[0], by = b[1] - a[1];
    const double cx = c[0] - a[0], cy = c[1] - a[1];
    const double den = 2.0 * (bx * cy - by * cx);
    if (std::fabs(den) < 1e-14) return false;
    const double b2 = bx * bx + by * by, c2 = cx * cx + cy * cy;
    out[0] = a[0] + (cy * b2 - by * c2) / den;
    out[1] = a[1] + (bx * c2 - cx * b2) / den;
    return true;
}

} // namespace

std::shared_ptr<Mesh> polycrystalDomain(int grains, unsigned seed, double h) {
    if (grains < 3) throw std::runtime_error("polycrystalDomain: need at least three grains");

    // --- the seeds ---------------------------------------------------------
    // A jittered hexagonal packing inside a disk, so the cells come out of
    // comparable size, plus a ring of ghost seeds outside it. The ghosts are
    // what make every *real* cell bounded, which is what lets the aggregate be
    // built without clipping anything: the domain is the union of the real
    // cells, and its boundary is the walls that only one real cell owns.
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> jit(-0.28, 0.28);

    std::vector<Point> site;
    const double R = 1.0;
    const int rings = std::max(1, static_cast<int>(std::ceil((std::sqrt(1.0 + 4.0 * (grains - 1) / 3.0) - 1.0) / 2.0)));
    const double step = R / (rings + 0.5);
    site.push_back({0.0, 0.0});
    for (int i = 1; i <= rings && static_cast<int>(site.size()) < grains; ++i) {
        const int m = 6 * i;
        for (int k = 0; k < m && static_cast<int>(site.size()) < grains; ++k) {
            const double th = 2.0 * M_PI * k / m;
            site.push_back({i * step * std::cos(th), i * step * std::sin(th)});
        }
    }
    for (Point &p : site) { p[0] += jit(rng) * step; p[1] += jit(rng) * step; }
    const int nReal = static_cast<int>(site.size());

    const int nGhost = std::max(12, 3 * nReal / 2);
    const double rg = R + 2.0 * step;
    for (int k = 0; k < nGhost; ++k) {
        const double th = 2.0 * M_PI * k / nGhost;
        site.push_back({rg * std::cos(th), rg * std::sin(th)});
    }
    const int nAll = static_cast<int>(site.size());

    // --- Delaunay of the seeds, by the empty-circumcircle test --------------
    // O(n^4) and entirely adequate at this size; what matters is that every
    // Voronoi vertex is computed exactly once, from its triangle, so the three
    // cells that meet there are given the identical coordinate rather than
    // three roundings of it.
    struct Tri { int a, b, c; Point cc; };
    std::vector<Tri> tri;
    for (int i = 0; i < nAll; ++i)
        for (int j = i + 1; j < nAll; ++j)
            for (int k = j + 1; k < nAll; ++k) {
                Point cc{};
                if (!circumcenterOf3(site[i], site[j], site[k], cc)) continue;
                const double r2 = (site[i][0] - cc[0]) * (site[i][0] - cc[0])
                                + (site[i][1] - cc[1]) * (site[i][1] - cc[1]);
                bool empty = true;
                for (int m = 0; m < nAll && empty; ++m) {
                    if (m == i || m == j || m == k) continue;
                    const double d2 = (site[m][0] - cc[0]) * (site[m][0] - cc[0])
                                    + (site[m][1] - cc[1]) * (site[m][1] - cc[1]);
                    if (d2 < r2 * (1.0 - 1e-12)) empty = false;
                }
                if (empty) tri.push_back({i, j, k, cc});
            }
    if (tri.empty()) throw std::runtime_error("polycrystalDomain: no Delaunay triangles");

    // --- the cells: the circumcenters around each real seed, in order -------
    std::vector<std::vector<int>> incident(nAll);
    for (std::size_t t = 0; t < tri.size(); ++t) {
        incident[tri[t].a].push_back(static_cast<int>(t));
        incident[tri[t].b].push_back(static_cast<int>(t));
        incident[tri[t].c].push_back(static_cast<int>(t));
    }

    std::vector<std::vector<int>> cell(nReal);     // ring of triangle indices
    for (int s = 0; s < nReal; ++s) {
        std::vector<int> ts = incident[s];
        if (ts.size() < 3) throw std::runtime_error("polycrystalDomain: an unbounded cell");
        std::sort(ts.begin(), ts.end(), [&](int x, int y) {
            return std::atan2(tri[x].cc[1] - site[s][1], tri[x].cc[0] - site[s][0])
                 < std::atan2(tri[y].cc[1] - site[s][1], tri[y].cc[0] - site[s][0]);
        });
        cell[s] = ts;
    }

    // --- the walls, each owned once, and the ones only one cell owns --------
    std::map<std::pair<int, int>, std::vector<int>> wallOwners;   // (tri,tri) -> cells
    for (int s = 0; s < nReal; ++s) {
        const std::vector<int> &ts = cell[s];
        for (std::size_t i = 0; i < ts.size(); ++i) {
            const int u = ts[i], v = ts[(i + 1) % ts.size()];
            wallOwners[{std::min(u, v), std::max(u, v)}].push_back(s);
        }
    }

    std::vector<std::pair<int, int>> interior, boundary;
    for (const auto &[w, owners] : wallOwners) {
        if (owners.size() >= 2) interior.push_back(w);
        else                    boundary.push_back(w);
    }
    if (boundary.size() < 3) throw std::runtime_error("polycrystalDomain: no outer boundary");

    // Chain the boundary walls into one loop.
    std::map<int, std::vector<int>> at;
    for (const auto &[u, v] : boundary) { at[u].push_back(v); at[v].push_back(u); }
    for (const auto &[v, nb] : at)
        if (nb.size() != 2) throw std::runtime_error("polycrystalDomain: the outer boundary is not a simple loop");

    std::vector<int> loop;
    {
        int prev = -1, cur = boundary.front().first;
        for (std::size_t guard = 0; guard <= boundary.size(); ++guard) {
            loop.push_back(cur);
            const std::vector<int> &nb = at[cur];
            const int nxt = (nb[0] == prev) ? nb[1] : nb[0];
            prev = cur;
            cur = nxt;
            if (cur == loop.front()) break;
        }
    }
    if (loop.size() != boundary.size())
        throw std::runtime_error("polycrystalDomain: the outer boundary is not one loop");

    // --- the PSLG ----------------------------------------------------------
    PSLG g;
    g.in.h = h;

    std::vector<std::array<int, 2>> ring;
    for (std::size_t i = 0; i < loop.size(); ++i) {
        const auto side = g.segmentsAlong(tri[loop[i]].cc, tri[loop[(i + 1) % loop.size()]].cc, h);
        ring.insert(ring.end(), side.begin(), side.end());
    }
    g.exterior(ring);

    // Every interior wall once, distributed round-robin over the grains so that
    // each region loop has segments of its own to be listed with; which wall
    // goes to which grain carries no meaning, since Triangle floods a region id
    // from its seed until it meets any segment at all.
    std::vector<std::vector<std::array<int, 2>>> perGrain(nReal);
    for (std::size_t i = 0; i < interior.size(); ++i) {
        const auto segs = g.segmentsAlong(tri[interior[i].first].cc, tri[interior[i].second].cc, h);
        auto &dst = perGrain[i % nReal];
        dst.insert(dst.end(), segs.begin(), segs.end());
    }
    for (int s = 0; s < nReal; ++s) {
        if (perGrain[s].size() < 3)
            throw std::runtime_error("polycrystalDomain: too few walls to give every grain a loop");
        g.region(perGrain[s], site[s], s + 1);
    }
    return build(g);
}

std::vector<Domain> mechanismDomains(double h) {
    std::vector<Domain> d;
    auto add = [&](const std::string &name, std::shared_ptr<Mesh> m, const std::string &note) {
        Domain dom;
        dom.name = name;
        dom.mesh = std::move(m);
        dom.note = note;
        d.push_back(std::move(dom));
    };

    // The control first, so a reader meets the null case before the effect.
    add("comb", combDomain(4, h), "control: every interior angle 90 or 270 degrees");

    // The star family sweeps the tip angle. 90 degrees still has odd valleys,
    // so it is not a control; the comb is.
    for (double t : {90.0, 70.0, 55.0, 45.0})
        add("star5t" + std::to_string(static_cast<int>(t)), starDomain(5, t, h),
            "5 points, tip " + std::to_string(static_cast<int>(t)) + " deg, valley "
            + std::to_string(static_cast<int>(360.0 - 72.0 - t)) + " deg");

    // The laminate family sweeps the angle the interfaces meet dS at.
    for (double a : {90.0, 75.0, 60.0, 45.0})
        add("laminate" + std::to_string(static_cast<int>(a)), laminateDomain(5, a, h),
            "5 layers, interfaces at " + std::to_string(static_cast<int>(a))
            + " deg to dS, junction obliquity "
            + std::to_string(static_cast<int>(90.0 - a)) + " deg"
            + (a == 90.0 ? " (control)" : ""));

    // The application case.
    for (int n : {7, 19})
        add("grain" + std::to_string(n), polycrystalDomain(n, 12345u, h),
            std::to_string(n) + "-grain Voronoi polycrystal, generically oblique triple junctions");
    return d;
}

std::vector<Domain> junctionDomains(double h, double obliqueSectorDegrees) {
    std::vector<Domain> d;
    auto add = [&](const std::string &name, std::shared_ptr<Mesh> m, const std::string &note) {
        Domain dom;
        dom.name = name;
        dom.mesh = std::move(m);
        dom.note = note;
        d.push_back(std::move(dom));
    };
    add("tjunction", tJunctionDomain(h),         "3 materials, one triple point, every sector 90 or 180 degrees");
    add("quadruple", quadruplePointDomain(h),    "4 materials meeting at one vertex");
    add("thinlayer", thinLayerDomain(h),         "a band four elements thick across a block");
    add("inclusion", embeddedInclusionDomain(h), "one circular inclusion in a matrix");
    {
        const double s = obliqueSectorDegrees;
        std::string note = "3 materials at 120 degrees: no sector is a right angle";
        if (s != 120.0) {
            const int a = static_cast<int>(std::lround(s)), b = static_cast<int>(std::lround(360.0 - 2.0 * s));
            note = "3 materials, sectors " + std::to_string(a) + "/" + std::to_string(a) + "/" + std::to_string(b) + " degrees";
        }
        add("oblique", obliqueJunctionDomain(h, s), note);
    }
    return d;
}

} // namespace paper
