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

std::shared_ptr<Mesh> build(const PSLG &g, double minAngle = 28.0) {
    TriangleMesher2D::Options opts;
    opts.min_angle_degrees = minAngle;
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

std::shared_ptr<Mesh> obliqueJunctionDomain(double h) {
    PSLG g;
    g.in.h = h;

    // The rays land on the rim, so put a rim vertex exactly where each one does
    // by building the circle out of three arcs that start at those points.
    const int n = circlePoints(1.0, h);
    const int per = std::max(6, n / 3);
    std::vector<Point> rim;
    for (int k = 0; k < 3; ++k) {
        for (int i = 0; i < per; ++i) {
            const double a = (90.0 + 120.0 * (k + static_cast<double>(i) / per)) * M_PI / 180.0;
            rim.push_back({std::cos(a), std::sin(a)});
        }
    }
    std::vector<std::array<int, 2>> ring;
    for (std::size_t i = 0; i < rim.size(); ++i)
        ring.push_back({g.point(rim[i]), g.point(rim[(i + 1) % rim.size()])});
    g.exterior(ring);

    // Three rays at 90, 210 and 330 degrees: every sector is 120 degrees, so no
    // sector is a multiple of a right angle and no one cross is tangent to all
    // three interfaces at once.
    for (int k = 0; k < 3; ++k) {
        const double a = (90.0 + 120.0 * k) * M_PI / 180.0;
        const double s = (90.0 + 120.0 * k + 60.0) * M_PI / 180.0;
        g.region(g.polyline({{0.0, 0.0}, {std::cos(a), std::sin(a)}}, h),
                 {0.5 * std::cos(s), 0.5 * std::sin(s)}, k + 1);
    }
    return build(g);
}

std::vector<Domain> junctionDomains(double h) {
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
    add("oblique",   obliqueJunctionDomain(h),   "3 materials at 120 degrees: no sector is a right angle");
    return d;
}

} // namespace paper
