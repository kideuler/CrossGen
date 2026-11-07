#include "ZIPLINE/PatternRepair.hxx"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <sstream>

#include "ZIPLINE/TriangleLocator.hxx"

namespace {

// A face of the layout as the closed polygon of its arcs' points.
std::vector<Point> facePolygon(const QuadLayout &L, const QuadLayout::Face &f) {
    std::vector<Point> poly;
    for (const int d : f.darts) {
        std::vector<Point> pts = L.getArcs()[QuadLayout::arcOfDart(d)].pts;
        if (d & 1) std::reverse(pts.begin(), pts.end());
        for (size_t i = 0; i + 1 < pts.size(); ++i) poly.push_back(pts[i]);
    }
    return poly;
}

bool inside(const std::vector<Point> &P, const Point &q) {
    bool in = false;
    for (size_t i = 0, j = P.size() - 1; i < P.size(); j = i++)
        if (((P[i][1] > q[1]) != (P[j][1] > q[1])) &&
            (q[0] < (P[j][0] - P[i][0]) * (q[1] - P[i][1]) / (P[j][1] - P[i][1]) + P[i][0]))
            in = !in;
    return in;
}

double distanceToPolygon(const std::vector<Point> &P, const Point &q) {
    double best = 1e300;
    for (size_t i = 0, j = P.size() - 1; i < P.size(); j = i++) {
        const Point a = P[j], ab = P[i] - P[j];
        const double L2 = dotP(ab, ab);
        const double t = L2 > 0.0 ? std::max(0.0, std::min(1.0, dotP(q - a, ab) / L2)) : 0.0;
        best = std::min(best, normP(q - (a + ab * t)));
    }
    return best;
}

// P2's "barycenter of the partition": the area centroid, unless a point of a
// 40 x 40 grid over the face lies more than twice as deep inside it -- a
// curved sliver (the strip beside a cusp is one) can have its centroid
// outside itself or against a side.
Point interiorPoint(const std::vector<Point> &P) {
    double A = 0.0, cx = 0.0, cy = 0.0;
    for (size_t i = 0, j = P.size() - 1; i < P.size(); j = i++) {
        const double c = P[j][0] * P[i][1] - P[i][0] * P[j][1];
        A += c;
        cx += (P[j][0] + P[i][0]) * c;
        cy += (P[j][1] + P[i][1]) * c;
    }
    Point best{cx / (3.0 * A), cy / (3.0 * A)};
    double depth = inside(P, best) ? distanceToPolygon(P, best) : -1.0;
    double lo0 = 1e300, lo1 = 1e300, hi0 = -1e300, hi1 = -1e300;
    for (const Point &p : P) {
        lo0 = std::min(lo0, p[0]); hi0 = std::max(hi0, p[0]);
        lo1 = std::min(lo1, p[1]); hi1 = std::max(hi1, p[1]);
    }
    const int n = 40;
    for (int i = 1; i < n; ++i)
        for (int j = 1; j < n; ++j) {
            const Point q{lo0 + (hi0 - lo0) * i / n, lo1 + (hi1 - lo1) * j / n};
            if (!inside(P, q)) continue;
            const double d = distanceToPolygon(P, q);
            if (depth < 0.0 || (d > 2.0 * depth && d > 0.0)) { depth = d; best = q; }
        }
    return best;
}

// From p along the unit vector u, the distance to the polygon's boundary.
double rayToPolygon(const std::vector<Point> &P, const Point &p, const Point &u) {
    double best = 1e300;
    for (size_t i = 0, j = P.size() - 1; i < P.size(); j = i++) {
        const Point a = P[j], e = P[i] - P[j];
        const double den = cross2(u, e);
        if (std::fabs(den) < 1e-300) continue;
        const Point w = a - p;
        const double t = cross2(w, e) / den, s = cross2(w, u) / den;
        if (t > 1e-9 && s >= 0.0 && s <= 1.0) best = std::min(best, t);
    }
    return best;
}

std::string at(const Point &p) {
    std::ostringstream os;
    os.precision(3);
    os << std::fixed << "(" << p[0] << ", " << p[1] << ")";
    return os.str();
}

} // namespace

const char *PatternRepair::kindName(Kind k) {
    switch (k) {
        case Kind::Deficit:         return "deficit";
        case Kind::DeficitRemoval:  return "deficit removal";
        case Kind::PairAtTJunction: return "3-5 pair at a T-junction";
        case Kind::MoveStem:        return "stem moved";
        case Kind::CancelPair:      return "pair cancelled";
        case Kind::MoveIntoFace:    return "moved into a face";
    }
    return "?";
}

std::vector<ConformalField::Singular> PatternRepair::apply(const std::vector<ConformalField::Singular> &pattern,
                                                           const Edit &edit) {
    std::vector<ConformalField::Singular> out;
    for (size_t j = 0; j < pattern.size(); ++j)
        if (std::find(edit.remove.begin(), edit.remove.end(), static_cast<int>(j)) == edit.remove.end())
            out.push_back(pattern[j]);
    out.insert(out.end(), edit.add.begin(), edit.add.end());
    return out;
}

std::vector<PatternRepair::Edit> PatternRepair::propose(const Mesh &mesh,
                                                        const std::vector<ConformalField::Singular> &pattern,
                                                        const ConformalField *declined, const QuadLayout &L,
                                                        const LayoutBlocks &blocks,
                                                        const SeparatrixTrace &trace, const Options &opts) {
    std::vector<Edit> out;
    const auto &faces = L.getFaces();
    const std::vector<char> &isBlock = blocks.faceIsBlock();
    if (isBlock.size() != faces.size()) return out;
    const TriangleLocator locator(mesh);
    double h = 0.0;
    for (const auto &e : mesh.edges) h += normP(mesh.vertices[e[1]] - mesh.vertices[e[0]]);
    h = mesh.edges.empty() ? 1.0 : h / static_cast<double>(mesh.edges.size());

    auto addAt = [&](const Point &p, int d, Edit &e) {
        const int t = locator.locate(p);
        if (t < 0) return false;
        e.add.push_back(ConformalField::Singular{p, t, d});
        return true;
    };
    // Closer than this a +- pair is cancelled by ConformalField as a dipole
    // (Options::dipoleRadius is 1.5 mean edges), so it would edit nothing.
    const double minPairDistance = 3.0 * h;

    // The refused faces, with P2's barycentre of each.
    std::vector<int> refused, triangles, pentagons;
    std::vector<Point> centre(faces.size(), Point{0.0, 0.0});
    for (size_t i = 0; i < faces.size(); ++i) {
        if (isBlock[i]) continue;
        refused.push_back(static_cast<int>(i));
        centre[i] = interiorPoint(facePolygon(L, faces[i]));
        if (faces[i].corners == 3) triangles.push_back(static_cast<int>(i));
        if (faces[i].corners == 5) pentagons.push_back(static_cast<int>(i));
    }
    if (refused.empty()) return out;

    // ---- 1. Poincare-Hopf deficit, region by region -------------------------
    //
    // Each region (a component of the model cut along its interfaces) has its
    // own count, so the points go where the shortfall is: for a region short
    // of m quarters, m of its own refused faces (its triangles first) get a
    // valence-3 point at their barycentre (P2's split), or the m valence-5
    // points of the region nearest them are removed; a region over, the
    // mirror. The candidates are the combinations of a few choices per region.
    if (opts.deficit && declined && !declined->getReport().regionDeficit4.empty() &&
        declined->regionOf().size() == mesh.triangles.size()) {
        const std::vector<int> &def = declined->getReport().regionDeficit4;
        const std::vector<int> &regionOf = declined->regionOf();
        auto regionAt = [&](const Point &p) {
            const int t = locator.locate(p);
            return t >= 0 ? regionOf[t] : -1;
        };
        struct Need {
            int region = -1, m = 0, addD = 0;
            std::vector<std::vector<Point>> choices;   // m points each, best first
        };
        // m points spread inside one face: along the line between its two
        // farthest-apart points, at its interior point's offset from that
        // line (P2's Fig. 18, bottom: a two-cornered patch gets two valence-3
        // points side by side).
        auto spread = [&](int f, int m) {
            const auto poly = facePolygon(L, faces[f]);
            std::vector<Point> pts;
            if (poly.size() < 3) return pts;
            size_t ia = 0, ib = 0;
            double far = -1.0;
            for (size_t i = 0; i < poly.size(); ++i)
                for (size_t j = i + 1; j < poly.size(); ++j)
                    if (normP(poly[i] - poly[j]) > far) { far = normP(poly[i] - poly[j]); ia = i; ib = j; }
            const Point A = poly[ia], B = poly[ib], u = normalizeP(B - A);
            const Point nrm{-u[1], u[0]};
            const double off = dotP(centre[f] - A, nrm);
            for (int j = 1; j <= m; ++j) {
                Point x = A + (B - A) * (static_cast<double>(j) / (m + 1));
                x = x + nrm * (off - dotP(x - A, nrm));
                if (!inside(poly, x)) return std::vector<Point>();
                pts.push_back(x);
            }
            return pts;
        };
        // The layout's faces at a split sector's corner, in that region.
        auto facesAtCorner = [&](int vertex, int region) {
            std::vector<int> out2;
            if (vertex < 0) return out2;
            const Point &x = mesh.vertices[vertex];
            int node = -1;
            for (size_t k = 0; k < L.getNodes().size(); ++k)
                if (normP(L.getNodes()[k].pos - x) < 1e-9 * (1.0 + normP(x))) { node = static_cast<int>(k); break; }
            if (node < 0) return out2;
            for (size_t i = 0; i < faces.size(); ++i) {
                if (std::find(faces[i].nodes.begin(), faces[i].nodes.end(), node) == faces[i].nodes.end()) continue;
                const Point c = interiorPoint(facePolygon(L, faces[i]));
                const int t = locator.locate(c);
                if (t >= 0 && regionOf[t] == region) out2.push_back(static_cast<int>(i));
            }
            return out2;
        };
        std::vector<Need> needs;
        bool feasible = true;
        for (size_t r = 0; r < def.size() && feasible; ++r) {
            if (def[r] == 0) continue;
            Need nd;
            nd.region = static_cast<int>(r);
            nd.m = std::abs(def[r]);
            nd.addD = def[r] > 0 ? +1 : -1;
            if (nd.m > 3) { feasible = false; break; }
            // Candidate faces, best first: the region's refused faces of the
            // right kind, its other refused faces, then the faces at the
            // corners the field was short at.
            std::vector<int> cand;
            const std::vector<int> &kind = nd.addD > 0 ? triangles : pentagons;
            for (const int i : kind)
                if (regionAt(centre[i]) == nd.region) cand.push_back(i);
            for (const int i : refused)
                if (std::find(cand.begin(), cand.end(), i) == cand.end() && regionAt(centre[i]) == nd.region)
                    cand.push_back(i);
            for (const auto &sc : declined->splitCorners()) {
                if (sc.triangle < 0 || regionOf[sc.triangle] != nd.region) continue;
                for (const int i : facesAtCorner(sc.vertex, nd.region))
                    if (std::find(cand.begin(), cand.end(), i) == cand.end()) {
                        cand.push_back(i);
                        centre[i] = interiorPoint(facePolygon(L, faces[i]));
                    }
            }
            const int n = static_cast<int>(cand.size());
            // One point in each of m faces ...
            std::vector<int> pick(nd.m);
            std::function<void(int, int)> rec = [&](int k, int from) {
                if (k == nd.m) {
                    std::vector<Point> pts;
                    for (const int j : pick) pts.push_back(centre[cand[j]]);
                    nd.choices.push_back(pts);
                    return;
                }
                for (int j = from; j < n && nd.choices.size() < 16; ++j) { pick[k] = j; rec(k + 1, j + 1); }
            };
            if (n >= nd.m) rec(0, 0);
            // ... or all m in one face.
            if (nd.m >= 2)
                for (int j = 0; j < n && nd.choices.size() < 24; ++j) {
                    const auto pts = spread(cand[j], nd.m);
                    if (static_cast<int>(pts.size()) == nd.m) nd.choices.push_back(pts);
                }
            if (nd.choices.empty()) { feasible = false; break; }
            needs.push_back(std::move(nd));
        }
        if (feasible && !needs.empty()) {
            // Odometer over the regions' choices, the first region fastest.
            std::vector<size_t> digit(needs.size(), 0);
            const size_t cap = 96;
            for (size_t made = 0; made < cap; ++made) {
                Edit add, rem;
                add.kind = Kind::Deficit;
                rem.kind = Kind::DeficitRemoval;
                add.what = "a valence-3/5 point at the barycentre of";
                rem.what = "the opposite point removed at";
                bool addOk = true, remOk = true;
                for (size_t g = 0; g < needs.size(); ++g) {
                    const Need &nd = needs[g];
                    for (const Point &x : nd.choices[digit[g]]) {
                        if (!addAt(x, nd.addD, add)) addOk = false;
                        add.what += " " + at(x);
                        // The removal: the nearest opposite point in the region.
                        int q0 = -1;
                        double bd = 1e300;
                        for (size_t q = 0; q < pattern.size(); ++q) {
                            if (pattern[q].d != -nd.addD || pattern[q].triangle < 0 ||
                                regionOf[pattern[q].triangle] != nd.region)
                                continue;
                            const double d = normP(pattern[q].pos - x);
                            if (d < bd && std::find(rem.remove.begin(), rem.remove.end(), static_cast<int>(q)) ==
                                              rem.remove.end()) {
                                bd = d;
                                q0 = static_cast<int>(q);
                            }
                        }
                        if (q0 < 0) remOk = false;
                        else {
                            rem.remove.push_back(q0);
                            rem.what += " " + at(pattern[q0].pos);
                        }
                    }
                }
                add.priority = rem.priority = std::numeric_limits<double>::max();
                if (addOk) out.push_back(add);
                if (remOk) out.push_back(rem);
                size_t g = 0;
                for (; g < needs.size(); ++g) {
                    if (++digit[g] < needs[g].choices.size()) break;
                    digit[g] = 0;
                }
                if (g == needs.size()) break;
            }
        }
    }

    // ---- 4. cancel a pair near a refused face, or move a point into it -----
    if (opts.moveIntoFaces) {
        for (const int i : refused) {
            const auto poly = facePolygon(L, faces[i]);
            double radius = 0.0;
            for (const Point &p : poly) radius = std::max(radius, normP(p - centre[i]));
            std::vector<int> near;
            for (size_t j = 0; j < pattern.size(); ++j)
                if (normP(pattern[j].pos - centre[i]) < 1.5 * radius) near.push_back(static_cast<int>(j));
            for (size_t a = 0; a < near.size(); ++a)
                for (size_t b = a + 1; b < near.size(); ++b)
                    if (pattern[near[a]].d + pattern[near[b]].d == 0) {
                        Edit e;
                        e.kind = Kind::CancelPair;
                        e.priority = std::fabs(faces[i].area);
                        e.remove = {near[a], near[b]};
                        e.what = "the pair at " + at(pattern[near[a]].pos) + " and " + at(pattern[near[b]].pos) +
                                 " removed";
                        out.push_back(e);
                    }
            for (const int j : near) {
                Edit e;
                e.kind = Kind::MoveIntoFace;
                e.priority = std::fabs(faces[i].area);
                e.remove = {j};
                if (!addAt(centre[i], pattern[j].d, e)) continue;
                e.what = "the point at " + at(pattern[j].pos) + " moved to " + at(centre[i]);
                out.push_back(e);
            }
        }
    }

    // ---- 2 and 3: at each T-junction ---------------------------------------
    if (opts.pairsAtTJunctions || opts.moveStems) {
        for (const int n : blocks.tJunctionNodes()) {
            if (n < 0 || n >= static_cast<int>(L.getNodes().size())) continue;
            const QuadLayout::Node &nd = L.getNodes()[n];
            if (nd.darts.size() != 3) continue;
            // The face that runs straight past the T.
            int F = -1;
            for (size_t i = 0; i < faces.size() && F < 0; ++i)
                for (size_t k = 0; k < faces[i].nodes.size(); ++k)
                    if (faces[i].nodes[k] == n && !faces[i].isCorner[k]) { F = static_cast<int>(i); break; }
            if (F < 0) continue;
            // The widest sector at the T is that face's; the two darts bounding
            // it are the bar, the third is the stem.
            int w = 0;
            double widest = -1.0;
            for (int k = 0; k < 3; ++k) {
                double a = nd.angles[(k + 1) % 3] - nd.angles[k];
                while (a < 0.0) a += 2.0 * M_PI;
                while (a >= 2.0 * M_PI) a -= 2.0 * M_PI;
                if (a > widest) { widest = a; w = k; }
            }
            const int stemDart = nd.darts[(w + 2) % 3];
            const double stemAngle = nd.angles[(w + 2) % 3];
            const Point T = nd.pos;
            const Point u{-std::cos(stemAngle), -std::sin(stemAngle)};   // on into F
            const auto poly = facePolygon(L, faces[F]);
            const double D = std::min(rayToPolygon(poly, T, u), 1e3);

            if (opts.pairsAtTJunctions) {
                auto pair = [&](const Point &p3, const Point &p5, const char *how) {
                    if (normP(p3 - p5) < minPairDistance) return;
                    Edit e;
                    e.kind = Kind::PairAtTJunction;
                    e.priority = std::fabs(faces[F].area);
                    if (!addAt(p3, +1, e) || !addAt(p5, -1, e)) return;
                    e.what = std::string("a 3-5 pair at the T-junction at ") + at(T) + " (" + how + ")";
                    out.push_back(e);
                };
                pair(T + u * (0.5 * h), centre[F], "3 on the T, 5 in the face past it");
                pair(centre[F], T + u * (0.5 * h), "5 on the T, 3 in the face past it");
                pair(T + u * (0.3 * D), T + u * (0.7 * D), "along the stem's line, 3 nearer");
                pair(T + u * (0.7 * D), T + u * (0.3 * D), "along the stem's line, 5 nearer");
                pair(T + u * (2.0 * h), T + u * (5.0 * h), "along the stem's line, 3 nearer, close");
                pair(T + u * (5.0 * h), T + u * (2.0 * h), "along the stem's line, 5 nearer, close");
            }
            if (opts.moveStems) {
                // The stem's singular point, as the pattern holds it: the
                // trace's singularity it left from, matched by position.
                const QuadLayout::Arc &arc = L.getArcs()[QuadLayout::arcOfDart(stemDart)];
                if (arc.separatrix < 0 || arc.separatrix >= static_cast<int>(trace.separatrices.size())) continue;
                const Separatrix &sp = trace.separatrices[arc.separatrix];
                if (sp.originKind != SeparatrixOrigin::Singularity || sp.origin_singularity_id < 0 ||
                    sp.origin_singularity_id >= static_cast<int>(trace.singularities.size()))
                    continue;
                const Singularity &sg = trace.singularities[sp.origin_singularity_id];
                int j = -1;
                double bd = 1e300;
                for (size_t q = 0; q < pattern.size(); ++q) {
                    const double d = normP(pattern[q].pos - sg.coordinates);
                    if (pattern[q].d == sg.d && d < bd) { bd = d; j = static_cast<int>(q); }
                }
                if (j < 0 || bd > 3.0 * h) continue;
                for (const int side : {w, (w + 1) % 3}) {
                    const int dart = nd.darts[side];
                    const QuadLayout::Arc &bar = L.getArcs()[QuadLayout::arcOfDart(dart)];
                    const int other = (dart & 1) ? bar.a : bar.b;
                    const Point N = L.getNodes()[other].pos;
                    for (const double scale : {1.0, 0.5, 2.0}) {
                        Edit e;
                        e.kind = Kind::MoveStem;
                        e.priority = std::fabs(faces[F].area);
                        e.remove = {j};
                        if (!addAt(pattern[j].pos + (N - T) * scale, pattern[j].d, e)) continue;
                        std::ostringstream os;
                        os << "the stem's point at " << at(pattern[j].pos) << " moved towards the bar's node at "
                           << at(N) << " (x" << scale << ")";
                        e.what = os.str();
                        out.push_back(e);
                    }
                }
            }
        }
    }
    return out;
}
