// CrossGen: Mesh3Dto2Dfacesgmsh
// Main program: load a .step/.stp CAD file, take its planar faces, put each
// in a canonical pose inside [-1,1]^2 on the XY plane, mesh it with h = 2/np,
// and save the faces worth keeping -- the nice ones, each shape once -- as
// separate .obj files.
//
// ### Which faces are kept
//
// py/dataset.csv (2026-09-27) was built from every planar face of the MAMBO
// parts, and its 1148 faces turned out to be 423 shapes. The same face recurs
// within a part (the two caps of a prism) and across parts (67 disks in 32
// parts). A shape repeated n times is weighed n times by a learner, lands on
// both sides of a train/test split whatever the part grouping says, and gives
// its labels n chances to disagree. So a face is written only if it is
//
//   * one connected piece of triangles whose boundary loops touch nowhere;
//   * resolved at the mesh size: its shortest corner-to-corner stretch of
//     boundary is at least --min-run mesh edges long (default 3), and it is
//     at least --min-width mesh edges thick (default 4), thickness read as
//     2 A / P -- the width of a strip, the radius of a disk;
//   * a new shape: no face written before it, by this run or, through
//     --registry, by an earlier one, is the same up to moving, rotating,
//     mirroring and scaling.
//
// "The same" is read off the triangle mesh with mesh::BoundaryFeatures, the
// corner rule every metric and feature of the selection dataset uses: the
// number of boundary loops, the sorted corner angles, the sorted lengths and
// turnings of the stretches between corners, the perimeter, and the two
// principal second moments of area, each made scale-free by the area. Two
// faces agreeing on all of those to within a mesh's worth of rounding are
// taken as one shape; the second is logged as a duplicate of the first. With
// --all every planar face is written, as before this filter existed.
//
// ### The pose
//
// Each face is turned so that its minimum-area bounding rectangle lies along
// the axes, long side on x, and scaled so that the long side spans [-1, 1]. A
// face whose plane leans in the model used to be scaled by the projection of
// its 3D bounding box instead, which is larger than the face, so it came out
// smaller than [-1,1]^2 and meshed coarser relative to its own size than the
// same face upright -- identical shapes at different resolutions. The
// rectangle is found on sampled boundary curves, before meshing; where two
// rectangles tie (a right isosceles triangle has two), the squarer one wins.
//
// ### What is meshed and written
//
// Only the face: the part's other surfaces, and then every curve and point
// not on the face's boundary, are removed before meshing. Those used to stay
// in the model, were meshed as 1D entities, and every one of their nodes was
// written -- a median of ~2,300 vertices per face that no triangle used, some
// inside the face, which changed what UMBER and ATLAS computed on it. The
// .obj now holds exactly the nodes its triangles use.

#include <gmsh.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "mesh/BoundaryFeatures.hxx"
#include "mesh/Mesh.hxx"

namespace fs = std::filesystem;

// --------------------------------------------------------------------------
// Helper: 3-component vector ops
// --------------------------------------------------------------------------
struct Vec3 { double x, y, z; };

static Vec3 cross(const Vec3 &a, const Vec3 &b) {
    return {a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}
static double dot(const Vec3 &a, const Vec3 &b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
static double length(const Vec3 &v) { return std::sqrt(dot(v, v)); }
static Vec3 normalise(const Vec3 &v) {
    double l = length(v);
    return {v.x / l, v.y / l, v.z / l};
}

// --------------------------------------------------------------------------
// Determine whether a Gmsh surface (dim=2) entity is planar.
// We sample a grid of normals via gmsh::model::getNormal and check that
// they are all (nearly) identical.
// Returns true and fills `normal` and `origin` if planar.
// --------------------------------------------------------------------------
static bool isPlanarFace(int surfaceTag, Vec3 &normal, Vec3 &origin) {
    // Get parametric bounds of the surface
    std::vector<double> pmin, pmax;
    gmsh::model::getParametrizationBounds(2, surfaceTag, pmin, pmax);
    double umin = pmin[0], vmin = pmin[1];
    double umax = pmax[0], vmax = pmax[1];

    // Sample normals on a small grid
    const int N = 4;
    std::vector<double> params;
    for (int i = 0; i <= N; ++i)
        for (int j = 0; j <= N; ++j) {
            double u = umin + (umax - umin) * i / N;
            double v = vmin + (vmax - vmin) * j / N;
            params.push_back(u);
            params.push_back(v);
        }

    std::vector<double> normals;
    gmsh::model::getNormal(surfaceTag, params, normals);

    // Reference normal from first sample
    Vec3 n0 = {normals[0], normals[1], normals[2]};
    double tol = 1e-6;
    for (std::size_t i = 3; i < normals.size(); i += 3) {
        Vec3 ni = {normals[i], normals[i + 1], normals[i + 2]};
        Vec3 diff = {ni.x - n0.x, ni.y - n0.y, ni.z - n0.z};
        if (length(diff) > tol)
            return false;
    }

    normal = normalise(n0);

    // Evaluate a point on the surface to serve as the origin of the local frame
    std::vector<double> p0params = {params[0], params[1]};
    std::vector<double> p0coords;
    gmsh::model::getValue(2, surfaceTag, p0params, p0coords);
    origin = {p0coords[0], p0coords[1], p0coords[2]};

    return true;
}

// --------------------------------------------------------------------------
// Determine whether a planar face has a "simple" boundary – i.e. the kind
// of face --prune skips (rectangles, circles, half-disks, etc.).
//
// Heuristic: get the boundary curves and query each curve's geometric type.
// A face is considered simple if ALL boundary edges are elementary analytic
// curves (Line, Circle, Ellipse) AND there are at most `maxEdges` of them.
// --------------------------------------------------------------------------
static bool isSimpleFace(int surfaceTag, int maxEdges = 4) {
    gmsh::vectorpair input = {{2, surfaceTag}};
    gmsh::vectorpair boundary;
    gmsh::model::getBoundary(input, boundary, /*combined=*/false,
                             /*oriented=*/false, /*recursive=*/false);

    int nEdges = 0;
    for (auto &[dim, tag] : boundary) {
        if (dim != 1) continue; // only look at curves
        ++nEdges;
        std::string etype;
        gmsh::model::getType(dim, tag, etype);
        // BSpline, Unknown, or other non-analytic curves → not simple
        if (etype != "Line" && etype != "Circle" && etype != "Ellipse")
            return false;
    }

    return nEdges <= maxEdges;
}

// --------------------------------------------------------------------------
// Build an orthonormal frame (e1, e2) in the plane with given normal.
// --------------------------------------------------------------------------
static void buildLocalFrame(const Vec3 &n, Vec3 &e1, Vec3 &e2) {
    // Pick a vector not parallel to n
    Vec3 up = (std::fabs(n.z) < 0.9) ? Vec3{0, 0, 1} : Vec3{1, 0, 0};
    e1 = normalise(cross(n, up));
    e2 = normalise(cross(n, e1));
}

// --------------------------------------------------------------------------
// The pose: the minimum-area rectangle around points of the plane.
// --------------------------------------------------------------------------
struct Pose {
    Point axis{1.0, 0.0};     // unit vector that becomes +x
    Point centre{0.0, 0.0};   // rectangle centre, in the same coordinates as the points
    double longSide = 0.0;
};

// Andrew's monotone chain: the convex hull, counter-clockwise.
static std::vector<Point> convexHull(std::vector<Point> p) {
    std::sort(p.begin(), p.end());
    p.erase(std::unique(p.begin(), p.end()), p.end());
    if (p.size() < 3) return p;
    std::vector<Point> h(2 * p.size());
    size_t k = 0;
    for (size_t i = 0; i < p.size(); ++i) {
        while (k >= 2 && cross2(h[k - 1] - h[k - 2], p[i] - h[k - 2]) <= 0.0) --k;
        h[k++] = p[i];
    }
    for (size_t i = p.size() - 1, t = k + 1; i-- > 0;) {
        while (k >= t && cross2(h[k - 1] - h[k - 2], p[i] - h[k - 2]) <= 0.0) --k;
        h[k++] = p[i];
    }
    h.resize(k - 1);
    return h;
}

// Some side of the minimum-area enclosing rectangle lies along a hull edge,
// so trying every hull edge's direction finds it. Ties -- a right isosceles
// triangle has two rectangles of one area -- go to the shorter long side, so
// the face fills more of [-1,1]^2 and the choice does not flip on rounding.
static Pose minimumAreaPose(const std::vector<Point> &points) {
    Pose best;
    const std::vector<Point> hull = convexHull(points);
    if (hull.size() < 2) return best;
    double bestArea = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < hull.size(); ++i) {
        const Point d = normalizeP(hull[(i + 1) % hull.size()] - hull[i]);
        if (normP(d) == 0.0) continue;
        const Point n{-d[1], d[0]};
        double lo0 = std::numeric_limits<double>::infinity(), hi0 = -lo0, lo1 = lo0, hi1 = -lo0;
        for (const Point &q : hull) {
            lo0 = std::min(lo0, dotP(q, d));
            hi0 = std::max(hi0, dotP(q, d));
            lo1 = std::min(lo1, dotP(q, n));
            hi1 = std::max(hi1, dotP(q, n));
        }
        const double w = hi0 - lo0, h = hi1 - lo1, area = w * h;
        const double longSide = std::max(w, h);
        const bool better = area < bestArea * (1.0 - 1e-9) ||
                            (area <= bestArea * (1.0 + 1e-9) && longSide < best.longSide);
        if (!better) continue;
        bestArea = area;
        best.longSide = longSide;
        best.axis = (w >= h) ? d : n;
        best.centre = d * (0.5 * (lo0 + hi0)) + n * (0.5 * (lo1 + hi1));
    }
    return best;
}

// Points along the face's boundary curves, in the plane's (e1, e2) frame.
static std::vector<Point> boundarySamples(int surfaceTag, const Vec3 &origin, const Vec3 &e1, const Vec3 &e2) {
    gmsh::vectorpair curves;
    gmsh::model::getBoundary({{2, surfaceTag}}, curves, /*combined=*/false,
                             /*oriented=*/false, /*recursive=*/false);
    std::vector<Point> out;
    const int kSamples = 64;
    for (const auto &c : curves) {
        if (c.first != 1) continue;
        std::vector<double> lo, hi;
        gmsh::model::getParametrizationBounds(1, std::abs(c.second), lo, hi);
        if (lo.empty() || hi.empty()) continue;
        std::vector<double> t(kSamples + 1), xyz;
        for (int i = 0; i <= kSamples; ++i) t[i] = lo[0] + (hi[0] - lo[0]) * i / kSamples;
        gmsh::model::getValue(1, std::abs(c.second), t, xyz);
        for (size_t i = 0; i + 2 < xyz.size(); i += 3) {
            const Vec3 d{xyz[i] - origin.x, xyz[i + 1] - origin.y, xyz[i + 2] - origin.z};
            out.push_back({dot(d, e1), dot(d, e2)});
        }
    }
    return out;
}

// --------------------------------------------------------------------------
// The shape signature, and the test for two faces being one shape.
// --------------------------------------------------------------------------
struct Signature {
    std::string name;                // the .obj it was written as
    int loops = 0;
    double perimeter = 0.0;          // P / sqrt(A)
    double moment[2] = {0.0, 0.0};   // principal second moments of area / A^2, ascending
    std::vector<double> angles;      // corner angles, degrees, ascending
    std::vector<double> runs;        // corner-to-corner lengths / sqrt(A), ascending
    std::vector<double> turnings;    // their turning, degrees, ascending
};

static Signature signatureOf(const Mesh &mesh, const BoundaryFeatures &features) {
    Signature s;
    const double A = features.summary().area;
    const double rootA = std::sqrt(A);
    s.loops = static_cast<int>(features.loops().size());
    s.perimeter = features.summary().perimeter / rootA;
    for (const BoundaryFeatures::Corner &c : features.corners()) s.angles.push_back(c.angle * 180.0 / M_PI);
    for (const BoundaryFeatures::Run &r : features.runs()) {
        s.runs.push_back(r.length / rootA);
        s.turnings.push_back(r.turning * 180.0 / M_PI);
    }
    std::sort(s.angles.begin(), s.angles.end());
    std::sort(s.runs.begin(), s.runs.end());
    std::sort(s.turnings.begin(), s.turnings.end());

    // Second moments about the centroid, exact on each triangle.
    Point c{0.0, 0.0};
    for (const Triangle &t : mesh.triangles) {
        const Point &p0 = mesh.vertices[t[0]], &p1 = mesh.vertices[t[1]], &p2 = mesh.vertices[t[2]];
        const double a = 0.5 * std::fabs(cross2(p1 - p0, p2 - p0));
        c = c + (p0 + p1 + p2) * (a / 3.0);
    }
    c = c / A;
    double sxx = 0.0, syy = 0.0, sxy = 0.0;
    for (const Triangle &t : mesh.triangles) {
        const Point q0 = mesh.vertices[t[0]] - c, q1 = mesh.vertices[t[1]] - c, q2 = mesh.vertices[t[2]] - c;
        const double a = 0.5 * std::fabs(cross2(q1 - q0, q2 - q0));
        sxx += a / 6.0 * (q0[0] * q0[0] + q1[0] * q1[0] + q2[0] * q2[0] + q0[0] * q1[0] + q1[0] * q2[0] + q2[0] * q0[0]);
        syy += a / 6.0 * (q0[1] * q0[1] + q1[1] * q1[1] + q2[1] * q2[1] + q0[1] * q1[1] + q1[1] * q2[1] + q2[1] * q0[1]);
        sxy += a / 12.0 * (2.0 * (q0[0] * q0[1] + q1[0] * q1[1] + q2[0] * q2[1]) + q0[0] * q1[1] + q1[0] * q0[1] +
                           q1[0] * q2[1] + q2[0] * q1[1] + q0[0] * q2[1] + q2[0] * q0[1]);
    }
    const double mean = 0.5 * (sxx + syy), dev = std::sqrt(0.25 * (sxx - syy) * (sxx - syy) + sxy * sxy);
    s.moment[0] = (mean - dev) / (A * A);
    s.moment[1] = (mean + dev) / (A * A);
    return s;
}

// Within a mesh's worth of rounding on every count. The pose gives two copies
// of one shape the same scale, and so the same mesh size, so the tolerances
// only have to absorb the triangulation, not a change of resolution.
static bool sameShape(const Signature &a, const Signature &b) {
    auto near = [](double x, double y, double rel, double abs) {
        return std::fabs(x - y) <= abs + rel * std::max(std::fabs(x), std::fabs(y));
    };
    if (a.loops != b.loops || a.angles.size() != b.angles.size() || a.runs.size() != b.runs.size())
        return false;
    if (!near(a.perimeter, b.perimeter, 3e-3, 0.0)) return false;
    for (int k = 0; k < 2; ++k)
        if (!near(a.moment[k], b.moment[k], 5e-3, 0.0)) return false;
    for (size_t i = 0; i < a.angles.size(); ++i)
        if (!near(a.angles[i], b.angles[i], 0.0, 1.0)) return false;
    for (size_t i = 0; i < a.runs.size(); ++i)
        if (!near(a.runs[i], b.runs[i], 5e-3, 1e-3) || !near(a.turnings[i], b.turnings[i], 0.0, 1.5))
            return false;
    return true;
}

// One line per shape: name loops perimeter moment0 moment1 nAngles angles...
// nRuns runs... turnings...
static std::string toLine(const Signature &s) {
    std::ostringstream o;
    o << std::setprecision(10) << s.name << ' ' << s.loops << ' ' << s.perimeter << ' ' << s.moment[0] << ' '
      << s.moment[1] << ' ' << s.angles.size();
    for (double v : s.angles) o << ' ' << v;
    o << ' ' << s.runs.size();
    for (double v : s.runs) o << ' ' << v;
    for (double v : s.turnings) o << ' ' << v;
    return o.str();
}

static bool fromLine(const std::string &line, Signature &s) {
    std::istringstream in(line);
    size_t nAngles = 0, nRuns = 0;
    if (!(in >> s.name >> s.loops >> s.perimeter >> s.moment[0] >> s.moment[1] >> nAngles)) return false;
    s.angles.resize(nAngles);
    for (double &v : s.angles)
        if (!(in >> v)) return false;
    if (!(in >> nRuns)) return false;
    s.runs.resize(nRuns);
    s.turnings.resize(nRuns);
    for (double &v : s.runs)
        if (!(in >> v)) return false;
    for (double &v : s.turnings)
        if (!(in >> v)) return false;
    return true;
}

// --------------------------------------------------------------------------
// The mesh of the current model's face, with only the nodes its triangles use.
// --------------------------------------------------------------------------
static void faceMesh(int surfaceTag, std::vector<Point> &points, std::vector<Triangle> &triangles) {
    std::vector<std::size_t> nodeTags;
    std::vector<double> coords, params;
    gmsh::model::mesh::getNodes(nodeTags, coords, params, 2, surfaceTag, /*includeBoundary=*/true,
                                /*returnParametricCoord=*/false);
    std::unordered_map<std::size_t, size_t> position;
    for (size_t i = 0; i < nodeTags.size(); ++i) position[nodeTags[i]] = i;

    std::vector<int> types;
    std::vector<std::vector<std::size_t>> elementTags, elementNodeTags;
    gmsh::model::mesh::getElements(types, elementTags, elementNodeTags, 2, surfaceTag);
    std::vector<std::array<std::size_t, 3>> tris;
    for (size_t k = 0; k < types.size(); ++k) {
        const auto &nodes = elementNodeTags[k];
        if (types[k] == 2) {          // 3-node triangles
            for (size_t j = 0; j + 2 < nodes.size(); j += 3) tris.push_back({nodes[j], nodes[j + 1], nodes[j + 2]});
        } else if (types[k] == 9) {   // 6-node triangles: corners only
            for (size_t j = 0; j + 5 < nodes.size(); j += 6) tris.push_back({nodes[j], nodes[j + 1], nodes[j + 2]});
        } else if (types[k] == 3) {   // 4-node quads: split
            for (size_t j = 0; j + 3 < nodes.size(); j += 4) {
                tris.push_back({nodes[j], nodes[j + 1], nodes[j + 2]});
                tris.push_back({nodes[j], nodes[j + 2], nodes[j + 3]});
            }
        }
    }

    // Numbered in gmsh's node order, and counter-clockwise as the .obj reader
    // leaves them.
    std::unordered_map<std::size_t, int> index;
    for (const auto &t : tris)
        for (std::size_t tag : t) index[tag] = -1;
    points.clear();
    for (size_t i = 0; i < nodeTags.size(); ++i) {
        auto it = index.find(nodeTags[i]);
        if (it == index.end() || it->second >= 0) continue;
        it->second = static_cast<int>(points.size());
        points.push_back({coords[3 * i], coords[3 * i + 1]});
    }
    triangles.clear();
    for (const auto &t : tris) {
        if (index[t[0]] < 0 || index[t[1]] < 0 || index[t[2]] < 0) continue;   // a node gmsh did not list
        Triangle abc{index[t[0]], index[t[1]], index[t[2]]};
        if (cross2(points[abc[1]] - points[abc[0]], points[abc[2]] - points[abc[0]]) < 0.0) std::swap(abc[1], abc[2]);
        triangles.push_back(abc);
    }
}

static void writeOBJ(const std::string &path, const std::vector<Point> &points, const std::vector<Triangle> &triangles) {
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Failed to open OBJ file for writing: " + path);
    out.setf(std::ios::fixed, std::ios::floatfield);
    out << std::setprecision(17);
    out << "# OBJ generated by CrossGen Mesh3Dto2Dfacesgmsh\n";
    for (const Point &p : points) out << "v " << p[0] << ' ' << p[1] << ' ' << 0.0 << '\n';
    for (const Triangle &t : triangles) out << "f " << t[0] + 1 << ' ' << t[1] + 1 << ' ' << t[2] + 1 << '\n';
}

// Why a meshed face is not nice, or "" when it is.
static std::string notNice(const Mesh &mesh, const BoundaryFeatures &features, double h, double minRun, double minWidth) {
    if (mesh.triangles.empty()) return "no triangles";
    if (mesh.materialComponents.size() != 1) return "not one connected piece";
    for (const Triangle &t : mesh.triangles) {
        const double a = 0.5 * cross2(mesh.vertices[t[1]] - mesh.vertices[t[0]], mesh.vertices[t[2]] - mesh.vertices[t[0]]);
        if (!(a > 1e-12 * h * h)) return "a degenerate triangle";
    }
    std::set<int> seen;
    size_t onLoops = 0;
    for (const BoundaryFeatures::Loop &L : features.loops()) {
        onLoops += L.vertices.size();
        seen.insert(L.vertices.begin(), L.vertices.end());
    }
    if (seen.size() != onLoops) return "boundary loops touch";
    const BoundaryFeatures::Summary &s = features.summary();
    const double run = s.shortestRun * std::sqrt(s.area);
    if (run < minRun * h) {
        std::ostringstream o;
        o << "a feature " << std::setprecision(3) << run / h << " mesh edges long";
        return o.str();
    }
    const double width = 2.0 * s.area / s.perimeter;
    if (width < minWidth * h) {
        std::ostringstream o;
        o << "only " << std::setprecision(3) << width / h << " mesh edges thick";
        return o.str();
    }
    return "";
}

// --------------------------------------------------------------------------
int main(int argc, char **argv) {
    if (argc < 5) {
        std::cerr << "Usage: Mesh3Dto2Dfacesgmsh <input.step> <output_fmt> <start_id> <np> [options]\n"
                  << "  output_fmt:       printf-style path with %%d, e.g. geom_%%d.obj; %%d is the\n"
                  << "                    face's index among the file's planar faces, from start_id\n"
                  << "  --registry FILE   shapes already written (one line each); a face the same as\n"
                  << "                    one of them is skipped, and each face written is added\n"
                  << "  --min-run K       shortest corner-to-corner stretch, in mesh edges (default 3)\n"
                  << "  --min-width K     thinnest 2 A / P, in mesh edges (default 4)\n"
                  << "  --prune           also skip simple faces (rectangles, circles, half-disks, etc.)\n"
                  << "  --all             write every planar face: no nice or unique test\n";
        return 1;
    }

    const std::string inputStep = argv[1];
    const std::string outputFmt = argv[2]; // e.g. "outdir/geom_%d.obj"
    const int startId = std::stoi(argv[3]);
    const int np = std::stoi(argv[4]);

    bool pruneSimple = false, writeAll = false;
    double minRun = 3.0, minWidth = 4.0;
    std::string registryPath;
    for (int i = 5; i < argc; ++i) {
        const std::string a(argv[i]);
        if (a == "--prune") pruneSimple = true;
        else if (a == "--all") writeAll = true;
        else if (a == "--registry" && i + 1 < argc) registryPath = argv[++i];
        else if (a == "--min-run" && i + 1 < argc) minRun = std::stod(argv[++i]);
        else if (a == "--min-width" && i + 1 < argc) minWidth = std::stod(argv[++i]);
        else {
            std::cerr << "Error: unknown option '" << a << "'\n";
            return 1;
        }
    }
    if (np <= 0) {
        std::cerr << "Error: np must be a positive integer.\n";
        return 1;
    }

    // The shapes written so far, by earlier runs.
    std::vector<Signature> written;
    if (!registryPath.empty()) {
        std::ifstream in(registryPath);
        std::string line;
        while (std::getline(in, line)) {
            Signature s;
            if (fromLine(line, s)) written.push_back(std::move(s));
        }
    }

    try {
        gmsh::initialize();
        gmsh::option::setNumber("General.Terminal", 1);

        // Load STEP file via the OCC kernel
        gmsh::open(inputStep);
        gmsh::model::occ::synchronize();

        // Enumerate all surface (dim=2) entities
        std::vector<std::pair<int, int>> surfaces;
        gmsh::model::getEntities(surfaces, 2);

        std::cout << "Found " << surfaces.size() << " surface(s) in STEP file.\n";

        // Every planar face, in the file's order: its position here is its
        // index in the output name, whichever faces end up skipped.
        struct PlanarFace {
            int tag;
            Vec3 normal;
            Vec3 origin;
            bool simple;
        };
        std::vector<PlanarFace> planarFaces;
        for (auto &[dim, tag] : surfaces) {
            Vec3 n{}, o{};
            if (isPlanarFace(tag, n, o)) planarFaces.push_back({tag, n, o, isSimpleFace(tag)});
        }
        std::cout << "Identified " << planarFaces.size() << " planar face(s).\n";

        if (planarFaces.empty()) {
            std::cerr << "No planar faces found – nothing to mesh.\n";
            gmsh::finalize();
            return 0;
        }

        fs::path fmtPath(outputFmt);
        if (fmtPath.has_parent_path())
            fs::create_directories(fmtPath.parent_path());

        std::ofstream registry;
        if (!registryPath.empty()) registry.open(registryPath, std::ios::app);

        const double h = 2.0 / static_cast<double>(np);
        int nWritten = 0, nSimple = 0, nNotNice = 0, nDuplicate = 0;

        for (size_t f = 0; f < planarFaces.size(); ++f) {
            const PlanarFace &face = planarFaces[f];
            const int faceIndex = startId + static_cast<int>(f);
            char buf[1024];
            std::snprintf(buf, sizeof(buf), outputFmt.c_str(), faceIndex);
            const std::string outPath(buf);
            if (pruneSimple && face.simple) {
                ++nSimple;
                std::cout << "  Face " << faceIndex << ": skipped, a simple face\n";
                continue;
            }

            // ---- Build a fresh model for this face ----
            // Re-open the STEP so that we get a clean copy each time
            std::string modelName = "face_" + std::to_string(faceIndex);
            gmsh::model::add(modelName);
            gmsh::model::setCurrent(modelName);

            gmsh::open(inputStep);
            gmsh::model::occ::synchronize();

            // Remove the volumes, then every other surface, then every curve
            // and point that is not on this face's boundary, so that nothing
            // but the face is meshed.
            {
                std::vector<std::pair<int, int>> vols;
                gmsh::model::getEntities(vols, 3);
                if (!vols.empty())
                    gmsh::model::occ::remove(vols, false);
                gmsh::model::occ::synchronize();
            }
            {
                std::vector<std::pair<int, int>> surfs, toRemove;
                gmsh::model::getEntities(surfs, 2);
                for (auto &[d, t] : surfs)
                    if (t != face.tag) toRemove.push_back({d, t});
                if (!toRemove.empty())
                    gmsh::model::occ::remove(toRemove, false);
                gmsh::model::occ::synchronize();
            }
            try {
                gmsh::vectorpair keepCurves, keepPoints;
                gmsh::model::getBoundary({{2, face.tag}}, keepCurves, false, false, false);
                gmsh::model::getBoundary({{2, face.tag}}, keepPoints, false, false, true);
                std::set<std::pair<int, int>> keep;
                for (const auto &e : keepCurves) keep.insert({e.first, std::abs(e.second)});
                for (const auto &e : keepPoints) keep.insert({e.first, std::abs(e.second)});
                for (int dim : {1, 0}) {
                    std::vector<std::pair<int, int>> ents, toRemove;
                    gmsh::model::getEntities(ents, dim);
                    for (const auto &e : ents)
                        if (!keep.count(e)) toRemove.push_back(e);
                    if (!toRemove.empty()) gmsh::model::occ::remove(toRemove, false);
                    gmsh::model::occ::synchronize();
                }
            } catch (const std::exception &e) {
                // Harmless if it fails: the .obj only takes the face's own nodes.
                std::cerr << "  Face " << faceIndex << ": could not remove stray entities (" << e.what() << ")\n";
            }

            // ---- The pose ----
            Vec3 e1{}, e2{};
            buildLocalFrame(face.normal, e1, e2);
            const Pose pose = minimumAreaPose(boundarySamples(face.tag, face.origin, e1, e2));
            if (!(pose.longSide > 0.0)) {
                std::cout << "  Face " << faceIndex << ": skipped, a degenerate face\n";
                ++nNotNice;
                gmsh::model::remove();
                continue;
            }
            const double s = 2.0 / pose.longSide;
            // (u, v) = ((P - o).e1, (P - o).e2) in the plane; then x along
            // pose.axis, y 90 degrees counter-clockwise from it, about the
            // rectangle's centre:
            //   x' = s ((u - cu) X_u + (v - cv) X_v),  y' = s ((u - cu) Y_u + (v - cv) Y_v).
            const Point X = pose.axis, Y{-X[1], X[0]};
            const Vec3 ex{X[0] * e1.x + X[1] * e2.x, X[0] * e1.y + X[1] * e2.y, X[0] * e1.z + X[1] * e2.z};
            const Vec3 ey{Y[0] * e1.x + Y[1] * e2.x, Y[0] * e1.y + Y[1] * e2.y, Y[0] * e1.z + Y[1] * e2.z};
            const double tx = -s * (dot(ex, face.origin) + dotP(pose.centre, X));
            const double ty = -s * (dot(ey, face.origin) + dotP(pose.centre, Y));
            std::vector<double> A = {
                s * ex.x, s * ex.y, s * ex.z, tx,
                s * ey.x, s * ey.y, s * ey.z, ty,
                0,        0,        0,        0,
                0,        0,        0,        1
            };
            std::vector<std::pair<int, int>> remaining;
            gmsh::model::getEntities(remaining);
            if (!remaining.empty())
                gmsh::model::occ::affineTransform(remaining, A);
            gmsh::model::occ::synchronize();

            // Set uniform mesh size and mesh
            std::vector<std::pair<int, int>> points;
            gmsh::model::getEntities(points, 0);
            if (!points.empty())
                gmsh::model::mesh::setSize(points, h);
            gmsh::model::mesh::generate(2);

            std::vector<Point> xy;
            std::vector<Triangle> tris;
            faceMesh(face.tag, xy, tris);
            gmsh::model::remove();

            if (!writeAll) {
                const Mesh mesh(xy, tris);
                const BoundaryFeatures features(mesh);
                const std::string why = notNice(mesh, features, h, minRun, minWidth);
                if (!why.empty()) {
                    ++nNotNice;
                    std::cout << "  Face " << faceIndex << ": skipped, " << why << '\n';
                    continue;
                }
                Signature sig = signatureOf(mesh, features);
                sig.name = fs::path(outPath).filename().string();
                const auto twin = std::find_if(written.begin(), written.end(),
                                               [&](const Signature &w) { return sameShape(sig, w); });
                if (twin != written.end()) {
                    ++nDuplicate;
                    std::cout << "  Face " << faceIndex << ": skipped, the same shape as " << twin->name << '\n';
                    continue;
                }
                writeOBJ(outPath, xy, tris);
                if (registry.is_open()) registry << toLine(sig) << '\n' << std::flush;
                written.push_back(std::move(sig));
            } else {
                writeOBJ(outPath, xy, tris);
            }
            ++nWritten;
            std::cout << "  Face " << faceIndex << " (surface tag " << face.tag << ") -> " << outPath << '\n';
        }

        gmsh::finalize();
        std::cout << "Done. Wrote " << nWritten << " of " << planarFaces.size() << " planar face mesh(es)";
        if (!writeAll)
            std::cout << "; skipped " << nDuplicate << " duplicate shape(s) and " << nNotNice << " not nice";
        if (pruneSimple) std::cout << ", " << nSimple << " simple";
        std::cout << ".\n";

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        try { gmsh::finalize(); } catch (...) {}
        return 1;
    } catch (...) {
        std::cerr << "Unknown error.\n";
        try { gmsh::finalize(); } catch (...) {}
        return 1;
    }

    return 0;
}
