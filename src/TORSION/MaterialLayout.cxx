#include "TORSION/MaterialLayout.hxx"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using Clock = std::chrono::steady_clock;

double secondsSince(Clock::time_point t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

long long pairKey(int a, int b) {
    const long long lo = std::min(a, b), hi = std::max(a, b);
    return (lo << 32) + hi;
}

double signedArea(const std::vector<Point> &loop) {
    double a = 0.0;
    for (size_t i = 0; i < loop.size(); ++i) {
        const Point &p = loop[i];
        const Point &q = loop[(i + 1) % loop.size()];
        a += p[0] * q[1] - q[0] * p[1];
    }
    return 0.5 * a;
}

double polylineLength(const std::vector<Point> &pts) {
    double l = 0.0;
    for (size_t i = 1; i < pts.size(); ++i) l += normP(pts[i] - pts[i - 1]);
    return l;
}

// What a region's layout is weighed by: its faces that are not simple
// quadrilaterals -- the ones Stage 10 has no grid for -- and its loose ends.
int defects(const TORSION *t) {
    if (!t || !t->hasArrangement()) return std::numeric_limits<int>::max() / 4;
    const Arrangement::Report &r = t->getArrangement().getReport();
    return (r.patches - r.simpleQuads) + r.danglingNodes + (r.patches == 0 ? 1 : 0);
}

bool acceptable(const TORSION *t) {
    return t && t->hasArrangement() && defects(t) == 0 && t->getStatus().layoutValid;
}

// Whether Stage 1 of this layout moved a +1 off a straight stretch of dS.
bool movedFlatCones(const TORSION *t) {
    return t && t->getStatus().flatConesToCorners + t->getStatus().flatConesInside > 0;
}

// How firmly a station holds its place when it is made one node with a station
// of the other side: a corner of its region's own layout does not move at all,
// an emitter sits on the vertex it was asked for, a node on a vertex of S beats
// one inside an edge. matchBranch() and assemble() place a pair by it alike.
int standing(const MaterialLayout::Station &x) {
    return x.corner ? 3 : x.emitter ? 2 : (x.vertex >= 0 ? 1 : 0);
}

// Move the ends of a separatrix onto the glued nodes it now ends at, bending the
// curve into them rather than moving only its last point. A match moves a curve's
// end along the interface by up to Options::matchTolerance edges, and moved alone
// the last point drags the last traced segment -- often shorter than the move --
// into a segment running along the interface instead of into it: the face between
// it and the interface is a sliver, its Coons patch folds, and Stage 12 smooths the
// fold into the interface. Spread over four times the move the curve turns by about
// twenty degrees at most; an arc shorter than that bends along its whole length, and
// its far end, a node it shares with other arcs, stays where it is either way.
void bendOnto(std::vector<Point> &pts, const Point &head, const Point &tail) {
    if (pts.size() < 2) {
        if (!pts.empty()) pts.front() = head;
        return;
    }
    std::vector<double> c(pts.size(), 0.0);
    for (size_t i = 1; i < pts.size(); ++i) c[i] = c[i - 1] + normP(pts[i] - pts[i - 1]);
    const double total = c.back();
    const Point d0 = head - pts.front(), d1 = tail - pts.back();
    const double l0 = std::min(total, 4.0 * normP(d0)), l1 = std::min(total, 4.0 * normP(d1));
    auto weight = [](double u) { return 1.0 - u * u * (3.0 - 2.0 * u); };
    for (size_t i = 1; i + 1 < pts.size(); ++i) {
        Point shift{0.0, 0.0};
        if (l0 > 0.0 && c[i] < l0) shift = shift + d0 * weight(c[i] / l0);
        if (l1 > 0.0 && total - c[i] < l1) shift = shift + d1 * weight((total - c[i]) / l1);
        pts[i] = pts[i] + shift;
    }
    pts.front() = head;
    pts.back() = tail;
}

// Whether a is the better of two layouts of one region: laid out at all, fewer
// defects, Definition 2.1 reached, fewer patches.
bool better(const TORSION *a, const TORSION *b) {
    const bool la = a && a->hasArrangement(), lb = b && b->hasArrangement();
    if (la != lb) return la;
    if (!la) return false;
    const int da = defects(a), db = defects(b);
    if (da != db) return da < db;
    const bool va = a->getStatus().layoutValid, vb = b->getStatus().layoutValid;
    if (va != vb) return va;
    return a->getArrangement().getReport().patches < b->getArrangement().getReport().patches;
}

} // namespace

// ---------------------------------------------------------------------------
MaterialLayout::Options MaterialLayout::optionsFor(const TORSION::Options &o) {
    Options mo;
    mo.rounds = o.perMaterialRounds;
    mo.threads = o.perMaterialThreads;
    mo.matchTolerance = o.perMaterialMatchTolerance;
    mo.keepDipolesOnFailure = o.perMaterialKeepDipoles;
    mo.keepFlatConesOnFailure = o.keepFlatConesOnFailure;
    mo.kinkAngle = o.interfaceKinkAngle;
    mo.spacedMatching = o.perMaterialSpacedMatching;
    mo.bendMatchedEnds = o.perMaterialBendEnds;
    return mo;
}

// Each region is the pipeline on a single-material mesh, Stages 0 to 8, with
// Stage 0 replaced by the field already solved on the whole model: restricted
// to a region it is aligned to the region's whole boundary, interfaces
// included, and it is the same field on both sides of every interface.
TORSION::Options MaterialLayout::regionOptionsFor(const TORSION::Options &o) {
    TORSION::Options r = o;
    r.perMaterial = false;
    r.diskTemplates = false;
    r.runSplines = false;
    r.runQuadMesh = false;
    r.externalField.resize(0);
    r.extraEmitters.clear();
    r.cancelSingleMaterialDipoles = o.cancelInterfaceDipoles;
    // A region does not make the whole-model run's second attempt with its
    // flat +1s in place: retry() does that here, with acceptable() -- every
    // face a simple quadrilateral as well as Definition 2.1 -- as the test, and
    // again after a round's relayout, which a region's own run never sees.
    r.keepFlatConesOnFailure = false;
    return r;
}

// ---------------------------------------------------------------------------
MaterialLayout::MaterialLayout(std::shared_ptr<Mesh> m, const Eigen::VectorXcd &f,
                               const TORSION::Options &ro, const Options &opts)
    : mesh(std::move(m)), field(f), regionOptions(ro), options(opts) {}

MaterialLayout::~MaterialLayout() = default;

// ---------------------------------------------------------------------------
// run()
// ---------------------------------------------------------------------------
bool MaterialLayout::run() {
    const auto t0 = Clock::now();
    report = Report();
    regions.clear();
    stations.clear();
    arrangement.reset();

    Interfaces::Options io;
    io.kinkAngle = options.kinkAngle;
    io.splitLoops = false;
    io.loopSplits = 0;
    network = std::make_unique<Interfaces>(mesh, io);
    buildRegions();
    buildBranches();
    report.regions = static_cast<int>(regions.size());

    if (field.size() != static_cast<Eigen::Index>(mesh->triangles.size())) {
        report.messages.push_back("the field has " + std::to_string(field.size()) +
                                  " value(s) for " + std::to_string(mesh->triangles.size()) +
                                  " triangle(s); nothing was laid out.");
        report.seconds = secondsSince(t0);
        return false;
    }

    // Largest first, so that the pool drains evenly: concrete's matrix alone is
    // most of the model.
    auto bySize = [&](std::vector<int> &jobs) {
        std::stable_sort(jobs.begin(), jobs.end(), [&](int a, int b) {
            return regions[a].triangles.size() > regions[b].triangles.size();
        });
    };

    // --- Every region on its own -------------------------------------------
    std::vector<int> all(regions.size());
    std::iota(all.begin(), all.end(), 0);
    bySize(all);
    parallel(all, [&](int r) { firstLayout(regions[r]); });

    auto everyRegionLaidOut = [&]() {
        bool ok = true;
        for (size_t r = 0; r < regions.size(); ++r) {
            if (regions[r].laidOut()) continue;
            ok = false;
            std::ostringstream oss;
            oss << "region " << r << " (material " << regions[r].material << ", "
                << regions[r].triangles.size() << " triangle(s)) did not reach Stage 8";
            const TORSION::Status &rs = regions[r].layout ? regions[r].layout->getStatus()
                                                          : TORSION::Status();
            for (const std::string &m : rs.messages) {
                if (m.rfind("Stopping", 0) == 0) { oss << ": " << m; break; }
            }
            report.messages.push_back(oss.str());
        }
        return ok;
    };
    bool laidOut = everyRegionLaidOut();

    // --- The rounds ----------------------------------------------------------
    for (int round = 1; laidOut && round <= options.rounds + 1; ++round) {
        stations = collectStations();
        const std::vector<std::vector<int>> want = wantedEmitters(stations);
        std::vector<int> jobs;
        std::vector<std::vector<int>> next(regions.size());
        int asked = 0;
        for (size_t r = 0; r < regions.size(); ++r) {
            std::set_union(regions[r].emitters.begin(), regions[r].emitters.end(),
                           want[r].begin(), want[r].end(), std::back_inserter(next[r]));
            if (next[r].size() == regions[r].emitters.size()) continue;
            asked += static_cast<int>(next[r].size() - regions[r].emitters.size());
            jobs.push_back(static_cast<int>(r));
        }
        report.askedPerRound.push_back(asked);
        if (jobs.empty()) {
            report.converged = true;
            break;
        }
        // Out of rounds, and not settled: what was asked for last is not laid
        // out, so it is not recorded as emitted either.
        if (round > options.rounds) break;
        for (int r : jobs) regions[r].emitters = std::move(next[r]);
        report.rounds = round;
        bySize(jobs);
        parallel(jobs, [&](int r) { nextLayout(regions[r]); });
        laidOut = everyRegionLaidOut();
    }

    for (const Region &r : regions) {
        if (r.laidOut()) ++report.regionsLaidOut;
        if (acceptable(r.layout.get())) ++report.regionsValid;
        if (r.dipolesKept) ++report.dipolesKept;
        if (r.flatConesKept) ++report.flatConesKept;
        if (r.flatConesKept || movedFlatCones(r.layout.get())) ++report.flatConesMoved;
        report.emitters += static_cast<int>(r.emitters.size());
        report.regionRuns += r.runs;
        report.regionRelayouts += r.relayouts;
    }
    for (size_t r = 0; r < regions.size(); ++r) {
        const Region &R = regions[r];
        if (!R.laidOut() || acceptable(R.layout.get())) continue;
        const Arrangement::Report &ar = R.layout->getArrangement().getReport();
        std::ostringstream oss;
        oss << "region " << r << " (material " << R.material << ", " << R.triangles.size()
            << " triangle(s)): ";
        if (ar.patches > ar.simpleQuads) {
            oss << (ar.patches - ar.simpleQuads) << " of its " << ar.patches
                << " face(s) are not simple quadrilaterals";
        } else if (ar.danglingNodes > 0) {
            oss << ar.danglingNodes << " separatrix/ces of its own end in the middle of a face";
        } else {
            oss << "its Stage 6 did not reach Definition 2.1";
        }
        oss << (R.dipolesKept ? " (with its +1/-1 pairs kept)." : ".");
        report.messages.push_back(oss.str());
    }
    if (report.flatConesKept > 0) {
        std::ostringstream oss;
        oss << report.flatConesKept << " of the " << report.flatConesMoved
            << " region(s) whose Stage 1 moved a +1 off a straight stretch of boundary "
               "could not be laid out with it moved, and keep it where it was: a patch "
               "corner of pi, which Stage 10 meshes as an element with a node in the "
               "middle of a side.";
        report.messages.push_back(oss.str());
    }
    if (laidOut && !report.converged) {
        std::ostringstream oss;
        oss << "the matching was still asking for new layout edges after " << options.rounds
            << " round(s); what it had asked for by then is what is glued.";
        report.messages.push_back(oss.str());
    }

    stations = collectStations();
    report.stations = static_cast<int>(stations.size());
    if (laidOut) report.assembled = assemble(stations);
    if (report.unmatched > 0 || report.extended > 0) {
        std::ostringstream oss;
        oss << report.unmatched << " layout vertex/vertices on the interfaces found no partner "
            << "on the other side; carrying them on through the faces beyond split "
            << report.extended << " face(s)";
        if (report.tJunctionsLeft > 0) {
            oss << ", and " << report.tJunctionsLeft << " node(s) are still in the middle of a "
                << "side, where Stage 10 has no grid for the face.";
        } else {
            oss << ", and every face is now a simple quadrilateral.";
        }
        report.messages.push_back(oss.str());
    }

    bool regionsValid = laidOut;
    for (const Region &r : regions) regionsValid = regionsValid && r.layout && r.layout->getStatus().layoutValid;
    report.valid = report.assembled && arrangement && arrangement->getReport().valid && regionsValid;
    report.seconds = secondsSince(t0);
    return report.valid;
}

// ---------------------------------------------------------------------------
// buildRegions()
//
// One mesh per region, its triangles in the order Interfaces lists them so
// that the field restricts by index. Vertices are renumbered in the order the
// triangles first use them; `local` inverts that.
// ---------------------------------------------------------------------------
void MaterialLayout::buildRegions() {
    const std::vector<Interfaces::Region> &ir = network->regions();
    regions.resize(ir.size());
    for (size_t r = 0; r < ir.size(); ++r) {
        Region &R = regions[r];
        R.material = ir[r].material;
        R.triangles = ir[r].triangles;
        std::vector<Point> V;
        std::vector<Triangle> T;
        T.reserve(R.triangles.size());
        for (int t : R.triangles) {
            Triangle lt{};
            for (int k = 0; k < 3; ++k) {
                const int g = mesh->triangles[t][k];
                auto it = R.local.find(g);
                if (it == R.local.end()) {
                    lt[k] = static_cast<int>(R.vertices.size());
                    R.local.emplace(g, lt[k]);
                    R.vertices.push_back(g);
                    V.push_back(mesh->vertices[g]);
                } else {
                    lt[k] = it->second;
                }
            }
            T.push_back(lt);
        }
        R.mesh = std::make_shared<Mesh>(V, T);
    }

    vertexRegions.assign(mesh->vertices.size(), {});
    const std::vector<int> &tr = network->triangleRegion();
    for (size_t t = 0; t < mesh->triangles.size() && t < tr.size(); ++t) {
        if (tr[t] < 0) continue;
        for (int k = 0; k < 3; ++k) {
            std::vector<int> &list = vertexRegions[mesh->triangles[t][k]];
            if (std::find(list.begin(), list.end(), tr[t]) == list.end()) list.push_back(tr[t]);
        }
    }
    edgeOf.clear();
    for (size_t e = 0; e < mesh->edges.size(); ++e) {
        edgeOf[pairKey(mesh->edges[e][0], mesh->edges[e][1])] = static_cast<int>(e);
    }
}

// ---------------------------------------------------------------------------
// buildBranches()
//
// Arc length along each branch, its interior vertices, and which region is on
// which side. Left is the region whose triangle uses the branch's first edge
// counter-clockwise, and it matters because it fixes the direction every
// region's boundary arc runs along the branch without looking at the arc: a
// region's dS is directed with the region on its left (Arrangement::
// buildDomain), so the region on the branch's left walks it forwards and the
// one on its right walks it backwards.
// ---------------------------------------------------------------------------
void MaterialLayout::buildBranches() {
    const std::vector<Interfaces::Branch> &ib = network->branches();
    const std::vector<int> &tr = network->triangleRegion();
    branches.assign(ib.size(), Branch());
    for (size_t b = 0; b < ib.size(); ++b) {
        Branch &B = branches[b];
        B.verts = ib[b].verts;
        B.closed = ib[b].closed;
        const size_t n = B.verts.size();
        B.cum.assign(n, 0.0);
        for (size_t i = 1; i < n; ++i) {
            const double l = normP(mesh->vertices[B.verts[i]] - mesh->vertices[B.verts[i - 1]]);
            B.cum[i] = B.cum[i - 1] + l;
            B.longestEdge = std::max(B.longestEdge, l);
        }
        B.length = n ? B.cum.back() : 0.0;
        // A closed loop repeats its first vertex at the end; an open branch
        // ends on two nodes of the network, which are not interior to it.
        for (size_t i = 0; i < n; ++i) {
            if (B.closed ? (i + 1 == n) : (i == 0 || i + 1 == n)) continue;
            B.at.emplace(B.verts[i], static_cast<int>(i));
        }
        if (n < 2 || ib[b].edges.empty()) continue;
        const int e = ib[b].edges.front();
        const int a = B.verts[0], c = B.verts[1];
        for (int k = 0; k < 2; ++k) {
            const int t = mesh->edgeTriangles[e][k];
            if (t < 0 || t >= static_cast<int>(tr.size())) continue;
            const Triangle &T = mesh->triangles[t];
            int ia = -1, ic = -1;
            for (int j = 0; j < 3; ++j) {
                if (T[j] == a) ia = j;
                if (T[j] == c) ic = j;
            }
            if (ia < 0 || ic < 0) continue;
            const bool inOrder = (ic == (ia + 1) % 3);
            const bool ccw = cross2(mesh->vertices[T[1]] - mesh->vertices[T[0]],
                                    mesh->vertices[T[2]] - mesh->vertices[T[0]]) > 0.0;
            if (inOrder == ccw) B.left = tr[t];
            else B.right = tr[t];
        }
    }
}

// ---------------------------------------------------------------------------
// Where a point of an edge of branch b is along it. `segmentOf` finds the edge
// {u, w} among the branch's, as the index of its first vertex; an end vertex of
// an open branch is found through the interior vertex beside it.
// ---------------------------------------------------------------------------
int MaterialLayout::segmentOf(int b, int u, int w) const {
    const Branch &B = branches[b];
    const int n = static_cast<int>(B.verts.size());
    const auto iu = B.at.find(u), iw = B.at.find(w);
    if (iu != B.at.end() && iw != B.at.end()) {
        const int i = std::min(iu->second, iw->second), j = std::max(iu->second, iw->second);
        if (j == i + 1) return i;
        if (B.closed && i == 0 && j == n - 2) return j;
        return -1;
    }
    if (iu != B.at.end() || iw != B.at.end()) {
        const int k = (iu != B.at.end()) ? iu->second : iw->second;
        const int end = (iu != B.at.end()) ? w : u;
        if (k == 1 && B.verts[0] == end) return 0;
        if (k == n - 2 && B.verts[n - 1] == end) return n - 2;
        return -1;
    }
    if (n == 2 && ((B.verts[0] == u && B.verts[1] == w) || (B.verts[0] == w && B.verts[1] == u))) {
        return 0;
    }
    return -1;
}

double MaterialLayout::along(int b, int segment, const Point &p) const {
    const Branch &B = branches[b];
    const Point &P = mesh->vertices[B.verts[segment]];
    const Point &Q = mesh->vertices[B.verts[segment + 1]];
    const Point d = Q - P;
    const double l2 = dotP(d, d);
    double t = l2 > 0.0 ? dotP(p - P, d) / l2 : 0.0;
    t = std::max(0.0, std::min(1.0, t));
    return B.cum[segment] + t * (B.cum[segment + 1] - B.cum[segment]);
}

Point MaterialLayout::pointAt(int b, double s) const {
    const Branch &B = branches[b];
    const int n = static_cast<int>(B.verts.size());
    if (n == 0) return Point{0.0, 0.0};
    if (B.closed && B.length > 0.0) {
        s = std::fmod(s, B.length);
        if (s < 0.0) s += B.length;
    }
    s = std::max(0.0, std::min(B.length, s));
    const auto it = std::upper_bound(B.cum.begin(), B.cum.end(), s);
    int i = static_cast<int>(it - B.cum.begin()) - 1;
    i = std::max(0, std::min(n - 2, i));
    const double span = B.cum[i + 1] - B.cum[i];
    const double t = span > 0.0 ? (s - B.cum[i]) / span : 0.0;
    const Point &P = mesh->vertices[B.verts[i]];
    const Point &Q = mesh->vertices[B.verts[i + 1]];
    return P + (Q - P) * t;
}

bool MaterialLayout::isRegionCone(const Region &r, int vertexOfS) const {
    if (!r.layout || !r.layout->hasCones()) return false;
    const auto it = r.local.find(vertexOfS);
    if (it == r.local.end()) return false;
    const std::vector<int> &I = r.layout->getCones().getIndices();
    return it->second < static_cast<int>(I.size()) && I[it->second] != 0;
}

// ---------------------------------------------------------------------------
// layOut() / firstLayout() / nextLayout()
// ---------------------------------------------------------------------------
void MaterialLayout::layOut(Region &r, bool keepDipoles, bool keepFlatCones) {
    TORSION::Options o = regionOptions;
    o.externalField.resize(static_cast<Eigen::Index>(r.triangles.size()));
    for (size_t i = 0; i < r.triangles.size(); ++i) {
        o.externalField[static_cast<Eigen::Index>(i)] = field[r.triangles[i]];
    }
    o.extraEmitters = r.emitters;
    if (keepDipoles) o.cancelSingleMaterialDipoles = false;
    if (keepFlatCones) o.relocateFlatCones = false;
    auto t = std::make_unique<TORSION>(r.mesh, o);
    try {
        t->run();
    } catch (const std::exception &) {
        // Left without an arrangement, which is what everyRegionLaidOut() reads.
    }
    r.layout = std::move(t);
    ++r.runs;
}

// The two retries, in the order they are tried. Keeping the +1/-1 pairs first,
// because it is the older of the two and the one a failure is more often down
// to; then putting back the +1s Stage 1 moved off a straight stretch of the
// region's boundary (TORSION::Options::relocateFlatCones), with the pairs as
// the first retry left them. The second is needed where the region's input mesh
// cannot turn the quarter the move leaves at an acute corner: a corner whose
// fan is a single triangle has one frame, the frame cannot turn inside it, and
// Stage 6 makes the corner a cusp -- multimat/tooth, dam, artery and one
// region of basin, each a layer pinching out between two interfaces. There the
// corner and the flat vertex beside it stay one cluster of two quarters, as
// they were, and the mesh keeps the element with a node in the middle of a
// side; everywhere else the move stands.
void MaterialLayout::retry(Region &r) {
    if (options.keepDipolesOnFailure && !r.dipolesKept && !acceptable(r.layout.get()) &&
        r.layout->getStatus().coneDipoleUnits > 0) {
        std::unique_ptr<TORSION> cancelled = std::move(r.layout);
        layOut(r, true, r.flatConesKept);
        if (better(r.layout.get(), cancelled.get())) r.dipolesKept = true;
        else r.layout = std::move(cancelled);
    }
    if (options.keepFlatConesOnFailure && !r.flatConesKept && !acceptable(r.layout.get()) &&
        movedFlatCones(r.layout.get())) {
        std::unique_ptr<TORSION> moved = std::move(r.layout);
        layOut(r, r.dipolesKept, true);
        if (better(r.layout.get(), moved.get())) r.flatConesKept = true;
        else r.layout = std::move(moved);
    }
}

void MaterialLayout::firstLayout(Region &r) {
    const auto t0 = Clock::now();
    r.dipolesKept = false;
    r.flatConesKept = false;
    layOut(r, false, false);
    retry(r);
    r.seconds += secondsSince(t0);
}

void MaterialLayout::nextLayout(Region &r) {
    const auto t0 = Clock::now();
    r.layout->relayout(r.emitters);
    ++r.relayouts;
    retry(r);
    r.seconds += secondsSince(t0);
}

// ---------------------------------------------------------------------------
// parallel()
//
// A pool of workers over the regions. Each region is a TORSION of its own on a
// mesh of its own, and nothing is shared between two of them but read-only
// inputs, so the order they finish in cannot change what comes out. The linear
// solves inside a region are OpenMP-parallel themselves (ParallelLDLT); each
// worker hands them its share of the threads, so that eight regions at a time
// do not each ask for eight.
// ---------------------------------------------------------------------------
void MaterialLayout::parallel(const std::vector<int> &jobs, const std::function<void(int)> &work) {
    if (jobs.empty()) return;
    int threads = options.threads > 0 ? options.threads
                                      : static_cast<int>(std::thread::hardware_concurrency());
    if (threads < 1) threads = 1;
    const int workers = std::min<int>(threads, static_cast<int>(jobs.size()));
    if (workers <= 1) {
        for (int j : jobs) work(j);
        return;
    }
    const int inner = std::max(1, threads / workers);
    std::atomic<size_t> next{0};
    std::mutex failure;
    std::string firstError;
    std::vector<std::thread> pool;
    pool.reserve(static_cast<size_t>(workers));
    for (int w = 0; w < workers; ++w) {
        pool.emplace_back([&, inner]() {
#ifdef _OPENMP
            omp_set_num_threads(inner);
#else
            (void)inner;
#endif
            for (;;) {
                const size_t k = next.fetch_add(1);
                if (k >= jobs.size()) break;
                try {
                    work(jobs[k]);
                } catch (const std::exception &e) {
                    std::lock_guard<std::mutex> lock(failure);
                    if (firstError.empty()) firstError = e.what();
                }
            }
        });
    }
    for (std::thread &t : pool) t.join();
    if (!firstError.empty()) report.messages.push_back("a region's layout threw: " + firstError);
}

// ---------------------------------------------------------------------------
// collectStations()
//
// Every node a region's own Stage 8 put on its boundary where that boundary is
// an interface: the corners of its layout there, its emitters, and where its
// separatrices land. Nodes on dS are the region's business alone.
// ---------------------------------------------------------------------------
std::vector<MaterialLayout::Station> MaterialLayout::collectStations() const {
    std::vector<Station> out;
    for (size_t r = 0; r < regions.size(); ++r) {
        const Region &R = regions[r];
        if (!R.laidOut()) continue;
        const std::unordered_set<int> emitting(R.emitters.begin(), R.emitters.end());
        const std::vector<Arrangement::Node> &nodes = R.layout->getArrangement().getNodes();
        for (size_t n = 0; n < nodes.size(); ++n) {
            const Arrangement::Node &nd = nodes[n];
            if (!nd.onBoundary || nd.out.empty()) continue;
            const bool cone = nd.kind == Arrangement::NodeKind::Cone;
            if (!cone && nd.kind != Arrangement::NodeKind::BoundaryCorner &&
                nd.kind != Arrangement::NodeKind::BoundaryHit) {
                continue;
            }
            Station st;
            st.region = static_cast<int>(r);
            st.node = static_cast<int>(n);
            st.p = nd.p;
            st.emitter = cone && nd.vertex >= 0 && emitting.count(nd.vertex) > 0;
            st.corner = (cone && !st.emitter) || nd.kind == Arrangement::NodeKind::BoundaryCorner;
            if (nd.vertex >= 0) {
                const int g = R.vertices[nd.vertex];
                st.vertex = g;
                if (network->nodeAt(g) >= 0) {
                    st.branch = -1;
                    out.push_back(st);
                    continue;
                }
                const int b = network->branchAt(g);
                if (b < 0) continue;
                const auto it = branches[b].at.find(g);
                if (it == branches[b].at.end()) continue;
                st.branch = b;
                st.s = branches[b].cum[it->second];
                out.push_back(st);
            } else if (nd.edge >= 0) {
                const int u = R.vertices[R.mesh->edges[nd.edge][0]];
                const int w = R.vertices[R.mesh->edges[nd.edge][1]];
                const auto ie = edgeOf.find(pairKey(u, w));
                if (ie == edgeOf.end()) continue;
                const int b = network->branchOfEdge(ie->second);
                if (b < 0) continue;
                const int seg = segmentOf(b, u, w);
                if (seg < 0) continue;
                st.branch = b;
                st.s = along(b, seg, nd.p);
                out.push_back(st);
            }
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// matchBranch()
//
// The two sides' stations on one branch are two sequences in the order they
// meet along it, and a matching of them is an alignment in the sense of
// sequence alignment: pairs never cross, so that making each pair one node
// keeps both sides' layout vertices in the order their faces have them, and a
// station can stay single. A pair is allowed within Options::matchTolerance
// edges of the branch, and two corners only where they are one vertex -- a
// corner is where a region's own boundary turns, and does not move.
//
// Pairing only neighbours along the branch, which is where this started, is
// not enough, and it fails at exactly the places the per-material mode is
// most needed: a small region whose map opens an acute corner lands three or
// four curves within one edge of it, the far side can only answer those with
// emitters on three or four separate vertices, and no one of them is then the
// neighbour of the curve it was placed for. So the alignment is free to pair
// across a station, and among alignments the one with the most pairs is
// taken, and among those the one that moves the curves' ends least in all.
//
// ### Where the nodes land, and why that is weighed first
//
// Most pairs, on its own, also takes alignments that are pairs only on paper.
// A pair is one node, at its corner if either station is one, else at its
// emitter, else half way (assemble()), and two pairs can each carry a station
// sitting on the same vertex: basin's left region has an emitter and the right
// one a corner on vertex 284 of a branch, the right one a landing half an edge
// before it, and the left an emitter one edge after. Paired in the order that
// makes two pairs -- emitter with landing, corner with the later emitter --
// both nodes belong on vertex 284. assemble() can only push one of them a
// twentieth of an edge on, and the glued layout then has an interface arc a
// twentieth of an edge long with a face either side of it pinched to a point
// there -- elements no smoothing can open. dam does the same one edge from a
// junction of its interfaces, where the push carries a region's corner off the
// vertex its boundary turns at.
//
// So the alignment knows where its nodes will land, and weighs first how far
// short of a quarter of an edge any two consecutive ones are -- or of half
// the room the region they share already had between its own two vertices
// there, where that is less, so that a region whose own layout puts two
// vertices close together is not counted against. A collision is then never
// cheaper than the matching it stands in for, which on basin is the corner
// paired with the emitter on its own vertex and the other two left single for
// the rounds to answer. Where no alignment collides this is the most-pairs
// rule exactly.
//
// On a closed loop the alignment starts after the widest gap between two
// stations, so that nothing that could pair is cut apart.
// ---------------------------------------------------------------------------
MaterialLayout::Pairing MaterialLayout::matchBranch(int b, const std::vector<int> &onBranch,
                                                    const std::vector<Station> &st) const {
    if (!options.spacedMatching) return matchBranchMostPairs(b, onBranch, st);
    Pairing out;
    const Branch &B = branches[b];
    if (onBranch.empty()) return out;

    // Arc length measured from the origin: 0 on an open branch, the middle of
    // the widest gap on a loop.
    std::vector<int> all = onBranch;
    std::sort(all.begin(), all.end(), [&](int x, int y) {
        if (st[x].s != st[y].s) return st[x].s < st[y].s;
        return x < y;
    });
    if (B.closed && B.length > 0.0 && all.size() > 1) {
        double widest = -1.0;
        for (size_t k = 0; k < all.size(); ++k) {
            const double a = st[all[k]].s;
            const double c = (k + 1 < all.size()) ? st[all[k + 1]].s : st[all[0]].s + B.length;
            if (c - a > widest) {
                widest = c - a;
                out.origin = a + 0.5 * (c - a);
            }
        }
        if (out.origin >= B.length) out.origin -= B.length;
    }
    auto rel = [&](int i) {
        double x = st[i].s - out.origin;
        if (B.closed && B.length > 0.0) {
            while (x < 0.0) x += B.length;
            while (x >= B.length) x -= B.length;
        }
        return x;
    };

    std::vector<int> L, R;
    for (int i : onBranch) (st[i].region == B.left ? L : R).push_back(i);
    auto byRel = [&](int x, int y) {
        const double a = rel(x), c = rel(y);
        if (a != c) return a < c;
        return x < y;
    };
    std::sort(L.begin(), L.end(), byRel);
    std::sort(R.begin(), R.end(), byRel);

    const int m = static_cast<int>(L.size()), n = static_cast<int>(R.size());
    const double tol = options.matchTolerance * B.longestEdge;
    auto allowed = [&](int i, int j) {
        const Station &x = st[i], &y = st[j];
        if (x.corner && y.corner && x.vertex != y.vertex) return false;
        return std::fabs(rel(i) - rel(j)) <= tol;
    };

    // A node of the glued layout is a station `x` alone (y < 0) or the pair
    // (x, y); placed as assemble() places it.
    auto place = [&](int x, int y) {
        if (y < 0) return rel(x);
        const int sx = standing(st[x]), sy = standing(st[y]);
        if (sx > 0 || sy > 0) return rel(sy > sx ? y : x);
        return 0.5 * (rel(x) + rel(y));
    };
    // The edge of the branch at arc length r from the origin.
    auto edgeAt = [&](double r) {
        if (B.cum.size() < 2) return B.longestEdge;
        double sAbs = r + out.origin;
        if (B.closed && B.length > 0.0) {
            while (sAbs >= B.length) sAbs -= B.length;
            while (sAbs < 0.0) sAbs += B.length;
        }
        size_t k = static_cast<size_t>(
            std::upper_bound(B.cum.begin(), B.cum.end(), sAbs) - B.cum.begin());
        k = std::min(std::max<size_t>(k, 1), B.cum.size() - 1);
        const double e = B.cum[k] - B.cum[k - 1];
        return e > 0.0 ? e : B.longestEdge;
    };
    // How far short of the room it should have a node at p1 is, following one
    // at p0. `room` is the least distance between two stations of one region,
    // one in each node (infinite when the two nodes share no region).
    auto shortfall = [&](double p0, double p1, double room) {
        const double need = std::min(0.25 * std::min(edgeAt(p0), edgeAt(p1)), 0.5 * room);
        return std::max(0.0, need - (p1 - p0));
    };
    const double inf = std::numeric_limits<double>::infinity();
    auto roomBetween = [&](int x0, int y0, int x1, int y1) {
        double room = inf;
        for (int u : {x0, y0}) {
            if (u < 0) continue;
            for (int w : {x1, y1}) {
                if (w >= 0 && st[u].region == st[w].region) room = std::min(room, std::fabs(rel(w) - rel(u)));
            }
        }
        return room;
    };

    // score[i][j][a]: the best alignment of the first i of L and the first j
    // of R whose last node is the start (a = 0), L[i-1] alone (1), R[j-1]
    // alone (2) or the pair (L[i-1], R[j-1]) (3). Best is least shortfall,
    // then most pairs, then least displacement.
    struct Score {
        double shortfall = 0.0;
        int pairs = 0;
        double moved = 0.0;
        bool reached = false;
    };
    auto better = [](const Score &a, const Score &c) {
        if (!a.reached) return false;
        if (!c.reached) return true;
        const double eps = 1e-12;
        if (a.shortfall < c.shortfall - eps) return true;
        if (c.shortfall < a.shortfall - eps) return false;
        if (a.pairs != c.pairs) return a.pairs > c.pairs;
        return a.moved < c.moved;
    };
    const int W = n + 1;
    auto at = [&](int i, int j, int a) { return (static_cast<size_t>(i) * W + j) * 4 + a; };
    std::vector<Score> score(static_cast<size_t>(m + 1) * W * 4);
    std::vector<int> from(score.size(), -1);
    score[at(0, 0, 0)].reached = true;
    // The last node of a state, as (x, y).
    auto lastNode = [&](int i, int j, int a, int &x, int &y) {
        x = y = -1;
        if (a == 1) x = L[i - 1];
        else if (a == 2) x = R[j - 1];
        else if (a == 3) { x = L[i - 1]; y = R[j - 1]; }
    };
    // What putting node (x, y) after state (i, j, a) falls short by. The
    // start of an open branch is a node of the network, a vertex of every
    // region along the branch.
    auto stepShort = [&](int i, int j, int a, int x, int y) {
        const double p1 = place(x, y);
        if (a == 0) {
            if (B.closed) return 0.0;
            double room = rel(x);
            if (y >= 0) room = std::min(room, rel(y));
            return shortfall(0.0, p1, room);
        }
        int x0, y0;
        lastNode(i, j, a, x0, y0);
        return shortfall(place(x0, y0), p1, roomBetween(x0, y0, x, y));
    };
    auto relax = [&](int i, int j, int a, int i2, int j2, int a2, double extra, bool pair, double moved) {
        Score s = score[at(i, j, a)];
        s.shortfall += extra;
        if (pair) ++s.pairs;
        s.moved += moved;
        if (better(s, score[at(i2, j2, a2)])) {
            score[at(i2, j2, a2)] = s;
            from[at(i2, j2, a2)] = static_cast<int>(at(i, j, a));
        }
    };
    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j <= n; ++j) {
            for (int a = 0; a < 4; ++a) {
                if (!score[at(i, j, a)].reached) continue;
                if (i < m) relax(i, j, a, i + 1, j, 1, stepShort(i, j, a, L[i], -1), false, 0.0);
                if (j < n) relax(i, j, a, i, j + 1, 2, stepShort(i, j, a, R[j], -1), false, 0.0);
                if (i < m && j < n && allowed(L[i], R[j])) {
                    relax(i, j, a, i + 1, j + 1, 3, stepShort(i, j, a, L[i], R[j]), true,
                          std::fabs(rel(L[i]) - rel(R[j])));
                }
            }
        }
    }
    // The end of an open branch is a node of the network too.
    int bestA = -1;
    Score best;
    for (int a = 0; a < 4; ++a) {
        Score s = score[at(m, n, a)];
        if (!s.reached) continue;
        if (!B.closed && a != 0) {
            int x0, y0;
            lastNode(m, n, a, x0, y0);
            double room = B.length - rel(x0);
            if (y0 >= 0) room = std::min(room, B.length - rel(y0));
            s.shortfall += shortfall(place(x0, y0), B.length, room);
        }
        if (bestA < 0 || better(s, best)) { best = s; bestA = a; }
    }

    std::vector<std::pair<int, int>> backwards;
    for (int k = static_cast<int>(at(m, n, bestA)); k >= 0 && from[k] >= 0; k = from[k]) {
        const int a = k % 4;
        const int ij = k / 4;
        const int i = ij / W, j = ij % W;
        if (a == 3) {
            backwards.push_back({L[i - 1], R[j - 1]});
            out.pairs.push_back({L[i - 1], R[j - 1]});
        } else if (a == 1) {
            backwards.push_back({L[i - 1], -1});
            out.single.push_back(L[i - 1]);
        } else if (a == 2) {
            backwards.push_back({R[j - 1], -1});
            out.single.push_back(R[j - 1]);
        }
    }
    out.sequence.assign(backwards.rbegin(), backwards.rend());
    return out;
}

MaterialLayout::Pairing MaterialLayout::matchBranchMostPairs(int b, const std::vector<int> &onBranch,
                                                             const std::vector<Station> &st) const {
    Pairing out;
    const Branch &B = branches[b];
    if (onBranch.empty()) return out;

    std::vector<int> all = onBranch;
    std::sort(all.begin(), all.end(), [&](int x, int y) {
        if (st[x].s != st[y].s) return st[x].s < st[y].s;
        return x < y;
    });
    if (B.closed && B.length > 0.0 && all.size() > 1) {
        double widest = -1.0;
        for (size_t k = 0; k < all.size(); ++k) {
            const double a = st[all[k]].s;
            const double c = (k + 1 < all.size()) ? st[all[k + 1]].s : st[all[0]].s + B.length;
            if (c - a > widest) {
                widest = c - a;
                out.origin = a + 0.5 * (c - a);
            }
        }
        if (out.origin >= B.length) out.origin -= B.length;
    }
    auto rel = [&](int i) {
        double x = st[i].s - out.origin;
        if (B.closed && B.length > 0.0) {
            while (x < 0.0) x += B.length;
            while (x >= B.length) x -= B.length;
        }
        return x;
    };

    std::vector<int> L, R;
    for (int i : onBranch) (st[i].region == B.left ? L : R).push_back(i);
    auto byRel = [&](int x, int y) {
        const double a = rel(x), c = rel(y);
        if (a != c) return a < c;
        return x < y;
    };
    std::sort(L.begin(), L.end(), byRel);
    std::sort(R.begin(), R.end(), byRel);

    const int m = static_cast<int>(L.size()), n = static_cast<int>(R.size());
    const double tol = options.matchTolerance * B.longestEdge;
    auto allowed = [&](int i, int j) {
        const Station &x = st[i], &y = st[j];
        if (x.corner && y.corner && x.vertex != y.vertex) return false;
        return std::fabs(rel(i) - rel(j)) <= tol;
    };
    // score[i][j]: the most pairs among the first i of L and the first j of R,
    // and the least total displacement giving that many.
    struct Score {
        int pairs = 0;
        double moved = 0.0;
        bool operator<(const Score &o) const {
            if (pairs != o.pairs) return pairs < o.pairs;
            return moved > o.moved;
        }
    };
    std::vector<std::vector<Score>> score(m + 1, std::vector<Score>(n + 1));
    std::vector<std::vector<char>> step(m + 1, std::vector<char>(n + 1, 0));
    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j <= n; ++j) {
            if (i == 0 && j == 0) continue;
            Score best;
            char how = 0;
            bool any = false;
            if (i > 0) { best = score[i - 1][j]; how = 1; any = true; }
            if (j > 0 && (!any || best < score[i][j - 1])) { best = score[i][j - 1]; how = 2; any = true; }
            if (i > 0 && j > 0 && allowed(L[i - 1], R[j - 1])) {
                Score s = score[i - 1][j - 1];
                ++s.pairs;
                s.moved += std::fabs(rel(L[i - 1]) - rel(R[j - 1]));
                if (best < s) { best = s; how = 3; }
            }
            score[i][j] = best;
            step[i][j] = how;
        }
    }
    std::vector<std::pair<int, int>> backwards;
    for (int i = m, j = n; i > 0 || j > 0;) {
        const char how = step[i][j];
        if (how == 3) {
            backwards.push_back({L[i - 1], R[j - 1]});
            out.pairs.push_back({L[i - 1], R[j - 1]});
            --i;
            --j;
        } else if (how == 1) {
            backwards.push_back({L[i - 1], -1});
            out.single.push_back(L[i - 1]);
            --i;
        } else {
            backwards.push_back({R[j - 1], -1});
            out.single.push_back(R[j - 1]);
            --j;
        }
    }
    // What the recursion leaves order-free -- a single on one side against a
    // single on the other -- is put in order of arc length.
    out.sequence.assign(backwards.rbegin(), backwards.rend());
    auto key = [&](const std::pair<int, int> &e) {
        return e.second < 0 ? rel(e.first) : 0.5 * (rel(e.first) + rel(e.second));
    };
    for (size_t k = 1; k < out.sequence.size(); ++k) {
        for (size_t q = k; q > 0; --q) {
            const auto &a = out.sequence[q - 1], &c = out.sequence[q];
            const bool bothSingle = a.second < 0 && c.second < 0;
            const bool otherSides = st[a.first].region != st[c.first].region;
            if (bothSingle && otherSides && key(c) < key(a)) std::swap(out.sequence[q - 1], out.sequence[q]);
            else break;
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// wantedEmitters()
//
// What each region is asked for by the stations as they stand. At a node of
// the network, every region around it needs a layout vertex there as soon as
// one of them has one -- a junction is a corner of some region, and the
// branches of the others leaving it are layout edges that have to end on
// something. Along a branch, a station that no neighbour on the far side
// answers is asked of the far side at the nearest vertex of the branch that is
// free: not one of its nodes, not a cone of that region, not a vertex already
// holding a layout vertex of that region, and not already asked for this round.
// ---------------------------------------------------------------------------
std::vector<std::vector<int>> MaterialLayout::wantedEmitters(const std::vector<Station> &st) const {
    std::vector<std::set<int>> want(regions.size());

    // Every node of the network is a vertex of every region around it, as
    // the whole-model route makes it one by emitting from it (Interfaces::
    // emitterNodes): a junction or a kink is a corner of some region, and an
    // interface landing on dS at a shallow angle is a corner of none, but the
    // side of a face that ran on through it from dS onto the interface would
    // be two arcs, and no arc of the glued layout runs past the end of its
    // branch.
    std::map<int, std::set<int>> has;
    for (const Station &x : st) {
        if (x.branch < 0 && x.vertex >= 0) has[x.vertex].insert(x.region);
    }
    for (const Interfaces::Node &nd : network->nodes()) {
        const int v = nd.vertex;
        if (v < 0 || v >= static_cast<int>(vertexRegions.size())) continue;
        const auto it = has.find(v);
        for (int r : vertexRegions[v]) {
            if ((it != has.end() && it->second.count(r)) || !regions[r].laidOut()) continue;
            const auto lv = regions[r].local.find(v);
            if (lv == regions[r].local.end() || isRegionCone(regions[r], v)) continue;
            want[r].insert(lv->second);
        }
    }

    std::vector<std::vector<int>> onBranch(branches.size());
    for (size_t i = 0; i < st.size(); ++i) {
        if (st[i].branch >= 0) onBranch[st[i].branch].push_back(static_cast<int>(i));
    }
    for (size_t b = 0; b < branches.size(); ++b) {
        if (onBranch[b].empty()) continue;
        const Branch &B = branches[b];
        const Pairing pr = matchBranch(static_cast<int>(b), onBranch[b], st);
        for (int i : pr.single) {
            const Station &x = st[i];
            const int other = (x.region == B.left) ? B.right : B.left;
            if (other < 0 || !regions[other].laidOut()) continue;
            const Region &O = regions[other];

            std::vector<std::pair<double, int>> candidates;
            for (const auto &kv : B.at) {
                double d = std::fabs(B.cum[kv.second] - x.s);
                if (B.closed) d = std::min(d, B.length - d);
                if (d <= options.matchTolerance * B.longestEdge) candidates.push_back({d, kv.second});
            }
            std::sort(candidates.begin(), candidates.end());
            for (const auto &c : candidates) {
                const int g = B.verts[c.second];
                const auto lv = O.local.find(g);
                if (lv == O.local.end() || isRegionCone(O, g)) continue;
                if (want[other].count(lv->second)) continue;
                bool occupied = false;
                for (int j : onBranch[b]) {
                    if (st[j].region == other && st[j].vertex == g) { occupied = true; break; }
                }
                if (occupied) continue;
                want[other].insert(lv->second);
                break;
            }
        }
    }

    std::vector<std::vector<int>> out(regions.size());
    for (size_t r = 0; r < regions.size(); ++r) out[r].assign(want[r].begin(), want[r].end());
    return out;
}

// ---------------------------------------------------------------------------
// coveredEdge()
//
// An edge of S that a boundary arc of a region runs along, which says whether
// the arc is on dS or on an interface, and on which branch.
// ---------------------------------------------------------------------------
int MaterialLayout::coveredEdge(const Region &R, const Arrangement::Arc &a) const {
    const std::vector<Arrangement::Node> &nodes = R.layout->getArrangement().getNodes();
    auto edgeBetween = [&](int u, int w) {
        const auto it = edgeOf.find(pairKey(u, w));
        return it == edgeOf.end() ? -1 : it->second;
    };
    auto vertexOf = [&](int n) {
        return (n >= 0 && nodes[n].vertex >= 0) ? R.vertices[nodes[n].vertex] : -1;
    };
    auto edgeOfNode = [&](int n) {
        if (n < 0 || nodes[n].edge < 0) return -1;
        return edgeBetween(R.vertices[R.mesh->edges[nodes[n].edge][0]],
                           R.vertices[R.mesh->edges[nodes[n].edge][1]]);
    };
    if (a.verts.size() >= 2) return edgeBetween(R.vertices[a.verts[0]], R.vertices[a.verts[1]]);
    if (a.verts.size() == 1) {
        const int v = R.vertices[a.verts[0]];
        const int u = vertexOf(a.from);
        if (u >= 0) return edgeBetween(u, v);
        const int e = edgeOfNode(a.from);
        if (e >= 0) return e;
        const int w = vertexOf(a.to);
        if (w >= 0) return edgeBetween(v, w);
        return edgeOfNode(a.to);
    }
    const int e0 = edgeOfNode(a.from);
    if (e0 >= 0) return e0;
    const int e1 = edgeOfNode(a.to);
    if (e1 >= 0) return e1;
    const int u = vertexOf(a.from), w = vertexOf(a.to);
    return (u >= 0 && w >= 0) ? edgeBetween(u, w) : -1;
}

// ---------------------------------------------------------------------------
// assemble()
//
// The regions' arrangements made one. Nodes: the final pairing's pairs are one
// node each, a station left single is a node of its own (a T-junction on the
// far side), every node of the network is one node for all the regions around
// it, and every other node -- a cone, a crossing, a node of dS -- belongs to
// one region and is carried over. Arcs: a region's separatrix and dS arcs are
// carried over, their ends moved onto the node they now are; the interface is
// cut afresh at the merged nodes into arcs shared by the two sides. Faces: each
// region's patches, their half-edges translated, a boundary arc of the region
// on an interface becoming the run of shared interface arcs between its two
// ends. Arrangement's assembling constructor then counts and checks.
// ---------------------------------------------------------------------------
bool MaterialLayout::assemble(const std::vector<Station> &st) {
    using NK = Arrangement::NodeKind;
    using AK = Arrangement::ArcKind;
    Arrangement::Assembly out;
    auto fail = [&](const std::string &why) {
        report.messages.push_back("could not glue the regions: " + why);
        return false;
    };

    std::vector<std::vector<int>> nodeMap(regions.size());
    for (size_t r = 0; r < regions.size(); ++r) {
        nodeMap[r].assign(regions[r].layout->getArrangement().getNodes().size(), -1);
    }
    auto addNode = [&](const Arrangement::Node &nd) {
        out.nodes.push_back(nd);
        return static_cast<int>(out.nodes.size()) - 1;
    };

    // --- Nodes of the network ---------------------------------------------
    std::map<int, int> networkNode;
    auto networkNodeAt = [&](int v) {
        const auto it = networkNode.find(v);
        if (it != networkNode.end()) return it->second;
        Arrangement::Node nd;
        nd.p = mesh->vertices[v];
        nd.kind = NK::InterfaceNode;
        nd.vertex = v;
        nd.onBoundary = v < static_cast<int>(mesh->isBoundaryVertex.size()) &&
                        mesh->isBoundaryVertex[v];
        const int id = addNode(nd);
        networkNode.emplace(v, id);
        return id;
    };
    for (const Station &x : st) {
        if (x.branch < 0 && x.vertex >= 0) nodeMap[x.region][x.node] = networkNodeAt(x.vertex);
    }

    // --- Stations along the branches: the final pairing --------------------
    std::vector<std::vector<int>> onBranch(branches.size());
    for (size_t i = 0; i < st.size(); ++i) {
        if (st[i].branch >= 0) onBranch[st[i].branch].push_back(static_cast<int>(i));
    }
    std::vector<std::vector<std::pair<double, int>>> alongBranch(branches.size());
    report.matched = report.unmatched = 0;
    report.maxMatchOffset = 0.0;
    unmatched.clear();
    for (size_t b = 0; b < branches.size(); ++b) {
        const Branch &B = branches[b];
        std::vector<std::pair<double, int>> inside;
        if (!onBranch[b].empty()) {
            const Pairing pr = matchBranch(static_cast<int>(b), onBranch[b], st);
            auto rel = [&](double sv) {
                double x = sv - pr.origin;
                if (B.closed && B.length > 0.0) {
                    while (x < 0.0) x += B.length;
                    while (x >= B.length) x -= B.length;
                }
                return x;
            };
            // Where each node of the sequence goes, measured from the origin:
            // at the corner if either station is one, else at the emitter,
            // else half way -- and then no nearer to the one before it than a
            // small fraction of an edge, forwards and then backwards, so that
            // two nodes that would have landed on one vertex stay two. They
            // are two layout vertices of each region; made one, the face
            // between them on either side would have a side of no length.
            const size_t count = pr.sequence.size();
            std::vector<double> at(count, 0.0);
            std::vector<int> atVertex(count, -1);
            for (size_t q = 0; q < count; ++q) {
                const auto &e = pr.sequence[q];
                if (e.second >= 0) {
                    const Station &x = st[e.first];
                    const Station &y = st[e.second];
                    const int sx = standing(x), sy = standing(y);
                    if (sx > 0 || sy > 0) {
                        const Station &keep = (sy > sx) ? y : x;
                        at[q] = rel(keep.s);
                        atVertex[q] = keep.vertex;
                    } else {
                        at[q] = 0.5 * (rel(x.s) + rel(y.s));
                    }
                } else {
                    at[q] = rel(st[e.first].s);
                    atVertex[q] = st[e.first].vertex;
                }
            }
            const double gapMin = 0.05 * B.longestEdge;
            const double lo = B.closed ? 0.0 : gapMin;
            const double hi = B.closed ? B.length - gapMin : B.length - gapMin;
            for (size_t q = 0; q < count; ++q) {
                const double floorAt = (q == 0) ? lo : at[q - 1] + gapMin;
                if (at[q] < floorAt) { at[q] = floorAt; atVertex[q] = -1; }
            }
            for (size_t q = count; q-- > 0;) {
                const double ceilAt = (q + 1 == count) ? hi : at[q + 1] - gapMin;
                if (at[q] > ceilAt) { at[q] = ceilAt; atVertex[q] = -1; }
            }
            for (size_t q = 0; q < count; ++q) {
                const auto &e = pr.sequence[q];
                Arrangement::Node nd;
                nd.kind = NK::InterfaceHit;
                if (e.second >= 0) {
                    const double gap = std::fabs(rel(st[e.first].s) - rel(st[e.second].s));
                    ++report.matched;
                    if (B.longestEdge > 0.0) {
                        report.maxMatchOffset = std::max(report.maxMatchOffset, gap / B.longestEdge);
                    }
                } else {
                    unmatched.push_back(e.first);
                    ++report.unmatched;
                }
                double sAt = at[q] + pr.origin;
                if (B.closed && B.length > 0.0 && sAt >= B.length) sAt -= B.length;
                nd.vertex = atVertex[q];
                nd.p = (nd.vertex >= 0) ? mesh->vertices[nd.vertex] : pointAt(static_cast<int>(b), sAt);
                const int id = addNode(nd);
                nodeMap[st[e.first].region][st[e.first].node] = id;
                if (e.second >= 0) nodeMap[st[e.second].region][st[e.second].node] = id;
                inside.push_back({sAt, id});
            }
        }
        if (B.closed) {
            // The sequence already runs round the loop from the origin; the
            // interface arcs are cut in that order.
            alongBranch[b] = std::move(inside);
            continue;
        }
        // An open branch runs from one node of the network to another, and
        // those are its first and last stop.
        if (!B.verts.empty()) {
            alongBranch[b].push_back({0.0, networkNodeAt(B.verts.front())});
            alongBranch[b].insert(alongBranch[b].end(), inside.begin(), inside.end());
            alongBranch[b].push_back({B.length, networkNodeAt(B.verts.back())});
        }
    }

    // --- Every other node, on demand ----------------------------------------
    auto mapNode = [&](int r, int n) {
        int &m = nodeMap[r][n];
        if (m >= 0) return m;
        const Region &R = regions[r];
        const Arrangement::Node &src = R.layout->getArrangement().getNodes()[n];
        if (src.vertex >= 0 && network->nodeAt(R.vertices[src.vertex]) >= 0) {
            m = networkNodeAt(R.vertices[src.vertex]);
            return m;
        }
        Arrangement::Node nd = src;
        nd.out.clear();
        nd.angle.clear();
        nd.cone = -1;
        nd.vertex = src.vertex >= 0 ? R.vertices[src.vertex] : -1;
        nd.face = src.face >= 0 ? R.triangles[src.face] : -1;
        nd.edge = -1;
        if (src.edge >= 0) {
            const int u = R.vertices[R.mesh->edges[src.edge][0]];
            const int w = R.vertices[R.mesh->edges[src.edge][1]];
            const auto ie = edgeOf.find(pairKey(u, w));
            if (ie != edgeOf.end()) {
                nd.edge = ie->second;
                nd.t = (mesh->edges[ie->second][0] == u) ? src.t : 1.0 - src.t;
            }
        }
        m = addNode(nd);
        return m;
    };

    // --- Separatrix and dS arcs, one for one ----------------------------------
    std::vector<std::vector<int>> arcMap(regions.size()), arcBranch(regions.size());
    for (size_t r = 0; r < regions.size(); ++r) {
        const Region &R = regions[r];
        const std::vector<Arrangement::Arc> &arcs = R.layout->getArrangement().getArcs();
        arcMap[r].assign(arcs.size(), -1);
        arcBranch[r].assign(arcs.size(), -1);
        for (size_t a = 0; a < arcs.size(); ++a) {
            const Arrangement::Arc &src = arcs[a];
            if (src.from < 0 || src.to < 0) continue;
            if (src.kind == AK::Boundary) {
                const int e = coveredEdge(R, src);
                const int b = e >= 0 ? network->branchOfEdge(e) : -1;
                if (b >= 0) {
                    // It has to stay on this one branch: running on past one of
                    // the network's nodes means the region has no layout vertex
                    // at a node where its neighbours do, which the rounds exist
                    // to prevent and a face cannot be glued across.
                    for (int v : src.verts) {
                        if (!branches[b].at.count(R.vertices[v])) {
                            std::ostringstream oss;
                            oss << "region " << r << "'s boundary runs past vertex "
                                << R.vertices[v] << " of the interface network without a "
                                << "layout vertex there";
                            return fail(oss.str());
                        }
                    }
                    arcBranch[r][a] = b;
                    continue;
                }
            }
            Arrangement::Arc arc = src;
            arc.from = mapNode(static_cast<int>(r), src.from);
            arc.to = mapNode(static_cast<int>(r), src.to);
            arc.curve = -1;
            for (Separatrices::Step &sp : arc.steps) {
                if (sp.face >= 0) sp.face = R.triangles[sp.face];
            }
            for (int &v : arc.verts) v = R.vertices[v];
            arc.faceFrom = src.faceFrom >= 0 ? R.triangles[src.faceFrom] : -1;
            arc.faceTo = src.faceTo >= 0 ? R.triangles[src.faceTo] : -1;
            arc.bdFrom = arc.bdTo = -1;
            arc.ifFrom = arc.ifTo = -1;
            if (arc.kind == AK::Separatrix && options.bendMatchedEnds) {
                bendOnto(arc.points, out.nodes[arc.from].p, out.nodes[arc.to].p);
            } else if (!arc.points.empty()) {
                arc.points.front() = out.nodes[arc.from].p;
                arc.points.back() = out.nodes[arc.to].p;
            }
            arc.modelLength = polylineLength(arc.points);
            arcMap[r][a] = static_cast<int>(out.arcs.size());
            out.arcs.push_back(std::move(arc));
        }
    }

    // --- The interfaces, cut at the merged nodes -------------------------------
    std::vector<std::vector<int>> branchArcs(branches.size());
    for (size_t b = 0; b < branches.size(); ++b) {
        const Branch &B = branches[b];
        const auto &list = alongBranch[b];
        const int k = static_cast<int>(list.size());
        const int n = static_cast<int>(B.verts.size());
        const int count = B.closed ? k : k - 1;
        for (int j = 0; j < count; ++j) {
            const int jn = (j + 1) % k;
            const double s0 = list[j].first;
            double s1 = list[jn].first;
            // Round a loop the sequence can pass s = 0 anywhere, not only at
            // its last arc; a loop with one node has one arc, all the way round.
            if (B.closed && (s1 <= s0 || k == 1)) s1 += B.length;
            Arrangement::Arc arc;
            arc.kind = AK::Interface;
            arc.from = list[j].second;
            arc.to = list[jn].second;
            arc.branch = static_cast<int>(b);
            arc.matLeft = B.left >= 0 ? regions[B.left].material : 0;
            arc.matRight = B.right >= 0 ? regions[B.right].material : 0;
            arc.points.push_back(out.nodes[arc.from].p);
            const double eps = 1e-12 * std::max(1.0, B.length);
            const int laps = B.closed ? 2 : 1;
            for (int lap = 0; lap < laps; ++lap) {
                for (int i = 0; i < n; ++i) {
                    if (B.closed && i == n - 1) continue;
                    const double si = B.cum[i] + lap * B.length;
                    if (si > s0 + eps && si < s1 - eps) {
                        arc.verts.push_back(B.verts[i]);
                        arc.points.push_back(mesh->vertices[B.verts[i]]);
                    }
                }
            }
            arc.points.push_back(out.nodes[arc.to].p);
            arc.modelLength = polylineLength(arc.points);
            arc.degenerate = !(arc.modelLength > 1e-12);
            branchArcs[b].push_back(static_cast<int>(out.arcs.size()));
            out.arcs.push_back(std::move(arc));
        }
    }

    // --- Half-edges -------------------------------------------------------------
    out.halves.assign(2 * out.arcs.size(), Arrangement::HalfEdge());
    for (size_t a = 0; a < out.arcs.size(); ++a) {
        Arrangement::HalfEdge &f = out.halves[2 * a];
        Arrangement::HalfEdge &g = out.halves[2 * a + 1];
        f.arc = g.arc = static_cast<int>(a);
        f.forward = true;
        g.forward = false;
        f.origin = out.arcs[a].from;
        g.origin = out.arcs[a].to;
        f.twin = static_cast<int>(2 * a + 1);
        g.twin = static_cast<int>(2 * a);
        f.next = f.twin;
        g.next = g.twin;
    }

    // The shared interface arcs a region's own boundary arc on branch b runs
    // along, from the merged node it starts at to the one it ends at.
    auto walk = [&](int b, int startNode, int endNode, bool forwards, std::vector<int> &seq) {
        const Branch &B = branches[b];
        const auto &list = alongBranch[b];
        const std::vector<int> &arcs = branchArcs[b];
        const int k = static_cast<int>(list.size());
        if (k == 0 || arcs.empty()) return false;
        auto findFirst = [&](int node) {
            for (int j = 0; j < k; ++j) if (list[j].second == node) return j;
            return -1;
        };
        auto findLast = [&](int node) {
            for (int j = k - 1; j >= 0; --j) if (list[j].second == node) return j;
            return -1;
        };
        if (!B.closed) {
            if (forwards) {
                const int i0 = findFirst(startNode);
                const int i1 = findLast(endNode);
                if (i0 < 0 || i1 < 0 || i1 <= i0) return false;
                for (int j = i0; j < i1; ++j) seq.push_back(2 * arcs[j]);
            } else {
                const int i0 = findLast(startNode);
                const int i1 = findFirst(endNode);
                if (i0 < 0 || i1 < 0 || i1 >= i0) return false;
                for (int j = i0 - 1; j >= i1; --j) seq.push_back(2 * arcs[j] + 1);
            }
            return true;
        }
        const int i0 = findFirst(startNode), i1 = findFirst(endNode);
        if (i0 < 0 || i1 < 0) return false;
        int j = i0;
        int guard = 0;
        if (forwards) {
            do {
                seq.push_back(2 * arcs[j]);
                j = (j + 1) % k;
            } while (j != i1 && ++guard <= k);
        } else {
            do {
                j = (j - 1 + k) % k;
                seq.push_back(2 * arcs[j] + 1);
            } while (j != i1 && ++guard <= k);
        }
        return guard <= k;
    };

    // --- Faces ------------------------------------------------------------------
    for (size_t r = 0; r < regions.size(); ++r) {
        const Region &R = regions[r];
        const Arrangement &A = R.layout->getArrangement();
        const std::vector<Arrangement::HalfEdge> &halves = A.getHalfEdges();
        const std::vector<Arrangement::Arc> &arcs = A.getArcs();
        const std::vector<Arrangement::Face> &faces = A.getFaces();
        auto expand = [&](int h, std::vector<int> &seq) {
            const Arrangement::HalfEdge &src = halves[h];
            const int a = src.arc;
            if (arcMap[r][a] >= 0) {
                seq.push_back(2 * arcMap[r][a] + (src.forward ? 0 : 1));
                return true;
            }
            const int b = arcBranch[r][a];
            if (b < 0) return false;
            const Arrangement::Arc &ra = arcs[a];
            const int startNode = mapNode(static_cast<int>(r), src.forward ? ra.from : ra.to);
            const int endNode = mapNode(static_cast<int>(r), src.forward ? ra.to : ra.from);
            const bool forwards = ((static_cast<int>(r) == branches[b].left) == src.forward);
            return walk(b, startNode, endNode, forwards, seq);
        };

        for (int f : A.patchFaces()) {
            const Arrangement::Face &src = faces[f];
            Arrangement::Face fc;
            std::map<int, std::vector<int>> expanded;
            for (size_t i = 0; i < src.half.size(); ++i) {
                std::vector<int> seq;
                if (!expand(src.half[i], seq) || seq.empty()) {
                    std::ostringstream oss;
                    oss << "a face of region " << r << " has a side on an interface that does "
                        << "not run between two of the glued nodes";
                    return fail(oss.str());
                }
                // A node the far side put in the middle of this side is one
                // the face passes straight through: two quarter turns of its
                // own, not a corner (Arrangement::classifyFaces()).
                for (size_t m = 0; m < seq.size(); ++m) {
                    fc.half.push_back(seq[m]);
                    const bool last = (m + 1 == seq.size());
                    fc.turns.push_back(last ? (i < src.turns.size() ? src.turns[i] : 2) : 2);
                }
                expanded.emplace(src.half[i], std::move(seq));
            }
            fc.tJunctions += src.tJunctions;
            for (int c : src.corners) fc.corners.push_back(mapNode(static_cast<int>(r), c));
            for (const std::vector<int> &side : src.sides) {
                std::vector<int> s2;
                for (int h : side) {
                    const auto it = expanded.find(h);
                    if (it != expanded.end()) s2.insert(s2.end(), it->second.begin(), it->second.end());
                }
                fc.sides.push_back(std::move(s2));
            }
            fc.quad = src.quad;
            fc.simple = fc.quad && fc.sides.size() == 4 && fc.sides[0].size() == 1 &&
                        fc.sides[1].size() == 1 && fc.sides[2].size() == 1 &&
                        fc.sides[3].size() == 1;
            fc.patch = true;
            fc.material = R.material;
            fc.mixed = false;

            std::vector<Point> loop;
            for (int h : fc.half) {
                const Arrangement::Arc &ar = out.arcs[out.halves[h].arc];
                if (out.halves[h].forward) {
                    for (size_t q = 0; q + 1 < ar.points.size(); ++q) loop.push_back(ar.points[q]);
                } else {
                    for (size_t q = ar.points.size(); q-- > 1;) loop.push_back(ar.points[q]);
                }
            }
            fc.area = loop.size() >= 3 ? signedArea(loop) : 0.0;

            const int id = static_cast<int>(out.faces.size());
            for (size_t i = 0; i < fc.half.size(); ++i) {
                Arrangement::HalfEdge &h = out.halves[fc.half[i]];
                if (h.face >= 0) {
                    std::ostringstream oss;
                    oss << "two faces claim the same side of arc " << h.arc
                        << " (regions overlap along an interface)";
                    return fail(oss.str());
                }
                h.face = id;
                h.next = fc.half[(i + 1) % fc.half.size()];
            }
            out.faces.push_back(std::move(fc));
        }
    }

    // --- What the matching left single, carried on through the faces ---------
    report.extended = extendThroughFaces(out);

    // --- Each node's outgoing half-edges ----------------------------------------
    for (size_t h = 0; h < out.halves.size(); ++h) {
        Arrangement::Node &nd = out.nodes[out.halves[h].origin];
        out.halves[h].slot = static_cast<int>(nd.out.size());
        nd.out.push_back(static_cast<int>(h));
        nd.angle.push_back(0.0);
    }

    Arrangement::Options aopts;
    aopts.mergeTolerance = regionOptions.arrangementMerge;
    aopts.cornerTolerance = regionOptions.arrangementCorner;
    aopts.collapseTolerance = regionOptions.arrangementCollapse;
    aopts.trimUnresolvedAtCrossings = regionOptions.arrangementTrim;
    try {
        arrangement = std::make_unique<Arrangement>(*mesh, std::move(out), aopts);
    } catch (const std::exception &e) {
        return fail(e.what());
    }
    return true;
}

// ---------------------------------------------------------------------------
// extendThroughFaces()
//
// The matching's leftovers, and what a region's own Stage 8 left, are nodes in
// the middle of a face's side: a vertex of the layout on one side of an arc
// that the face on the other side does not have as a corner. The rounds of
// matching answer such a vertex with an emitter, a layout edge traced through
// the far region's own map, and that is the better answer wherever it is
// available. Where it is not -- a region that lands three or four curves within
// one edge of a corner leaves the far side no vertex free to emit them from --
// the curve is carried on through the glued faces instead, each of which is a
// quadrilateral with a Coons parameterisation of its own: the node at arc
// length fraction t along one side is joined to the point at 1 - t along the
// opposite one (the same isoline of the face's (s, t) frame), splitting the face
// in two, and that point is a new node in the middle of a side of the face
// beyond. Carried on until it reaches dS or an existing node, it is the
// extension through the neighbour's parameterisation the per-material mode is
// built on, with the face's parameterisation standing in for the region's.
// ---------------------------------------------------------------------------
namespace {

std::vector<double> cumulativeLength(const std::vector<Point> &p) {
    std::vector<double> c(p.size(), 0.0);
    for (size_t i = 1; i < p.size(); ++i) c[i] = c[i - 1] + normP(p[i] - p[i - 1]);
    return c;
}

Point atLength(const std::vector<Point> &p, const std::vector<double> &c, double l) {
    if (p.empty()) return Point{0.0, 0.0};
    if (l <= 0.0 || p.size() == 1) return p.front();
    if (l >= c.back()) return p.back();
    size_t i = static_cast<size_t>(std::upper_bound(c.begin(), c.end(), l) - c.begin());
    i = std::max<size_t>(1, std::min(i, p.size() - 1));
    const double span = c[i] - c[i - 1];
    const double t = span > 0.0 ? (l - c[i - 1]) / span : 0.0;
    return p[i - 1] + (p[i] - p[i - 1]) * t;
}

} // namespace

int MaterialLayout::extendThroughFaces(Arrangement::Assembly &out) {
    using NK = Arrangement::NodeKind;
    using AK = Arrangement::ArcKind;
    std::vector<Arrangement::Node> &nodes = out.nodes;
    std::vector<Arrangement::Arc> &arcs = out.arcs;
    std::vector<Arrangement::HalfEdge> &halves = out.halves;
    std::vector<Arrangement::Face> &faces = out.faces;

    auto endOf = [&](int h) { return halves[halves[h].twin].origin; };
    auto halfPoints = [&](int h) {
        std::vector<Point> p = arcs[halves[h].arc].points;
        if (!halves[h].forward) std::reverse(p.begin(), p.end());
        return p;
    };
    auto halfLength = [&](int h) { return arcs[halves[h].arc].modelLength; };
    auto sidePolyline = [&](const std::vector<int> &side) {
        std::vector<Point> line;
        for (int h : side) {
            std::vector<Point> p = halfPoints(h);
            if (!line.empty() && !p.empty()) p.erase(p.begin());
            line.insert(line.end(), p.begin(), p.end());
        }
        return line;
    };
    // Corners, sides, flags and area of a face from its cycle and its turns,
    // as Arrangement::classifyFaces() reads them: a corner is a turn of one
    // quarter, and a side is the run of half-edges from one corner to the next.
    auto refresh = [&](Arrangement::Face &fc) {
        const size_t m = fc.half.size();
        fc.corners.clear();
        fc.sides.clear();
        fc.tJunctions = 0;
        int first = -1;
        for (size_t i = 0; i < m; ++i) {
            if (fc.turns[i] == 1 && first < 0) first = static_cast<int>(i);
            if (fc.turns[i] == 0) ++fc.tJunctions;
        }
        if (first >= 0) {
            std::vector<int> side;
            for (size_t k = 1; k <= m; ++k) {
                const size_t i = (first + k) % m;
                side.push_back(fc.half[i]);
                if (fc.turns[i] == 1) {
                    fc.sides.push_back(side);
                    side.clear();
                }
            }
            for (const std::vector<int> &sd : fc.sides) fc.corners.push_back(halves[sd.front()].origin);
        }
        fc.quad = fc.corners.size() == 4;
        fc.simple = fc.quad && fc.sides.size() == 4 && fc.sides[0].size() == 1 &&
                    fc.sides[1].size() == 1 && fc.sides[2].size() == 1 && fc.sides[3].size() == 1;
        std::vector<Point> loop;
        for (int h : fc.half) {
            std::vector<Point> p = halfPoints(h);
            for (size_t q = 0; q + 1 < p.size(); ++q) loop.push_back(p[q]);
        }
        fc.area = loop.size() >= 3 ? signedArea(loop) : 0.0;
        for (size_t i = 0; i < m; ++i) {
            halves[fc.half[i]].next = fc.half[(i + 1) % m];
        }
    };

    // Split arc a at arc length l from its `from` end; returns the new node.
    // The faces on either side get the node as a straight pass-through.
    auto splitArc = [&](int a, double l) {
        const std::vector<Point> pts = arcs[a].points;
        const std::vector<double> c = cumulativeLength(pts);
        l = std::max(0.0, std::min(c.back(), l));
        size_t i = static_cast<size_t>(std::upper_bound(c.begin(), c.end(), l) - c.begin());
        i = std::max<size_t>(1, std::min(i, pts.size() - 1));
        const Point P = atLength(pts, c, l);

        Arrangement::Node nd;
        nd.p = P;
        nd.kind = arcs[a].kind == AK::Interface ? NK::InterfaceHit
                  : arcs[a].kind == AK::Boundary ? NK::BoundaryHit
                                                 : NK::Crossing;
        nd.onBoundary = arcs[a].kind == AK::Boundary;
        const int y = static_cast<int>(nodes.size());
        nodes.push_back(nd);

        Arrangement::Arc second = arcs[a];
        second.from = y;
        second.points.assign(1, P);
        second.points.insert(second.points.end(), pts.begin() + static_cast<long>(i), pts.end());
        std::vector<Point> firstPts(pts.begin(), pts.begin() + static_cast<long>(i));
        firstPts.push_back(P);
        // Interior vertices of a boundary or interface arc sit at points
        // 1 .. n-2; the first i points keep the first i - 1 of them.
        if (!arcs[a].verts.empty()) {
            const std::vector<int> v = arcs[a].verts;
            const size_t keep = std::min(v.size(), i - 1);
            arcs[a].verts.assign(v.begin(), v.begin() + static_cast<long>(keep));
            second.verts.assign(v.begin() + static_cast<long>(keep), v.end());
        }
        arcs[a].points = firstPts;
        arcs[a].to = y;
        arcs[a].steps.clear();
        second.steps.clear();
        arcs[a].modelLength = polylineLength(arcs[a].points);
        second.modelLength = polylineLength(second.points);
        const int b = static_cast<int>(arcs.size());
        arcs.push_back(second);

        halves.resize(2 * arcs.size());
        Arrangement::HalfEdge &f = halves[2 * b];
        Arrangement::HalfEdge &g = halves[2 * b + 1];
        f = Arrangement::HalfEdge();
        g = Arrangement::HalfEdge();
        f.arc = g.arc = b;
        f.forward = true;
        g.forward = false;
        f.origin = y;
        g.origin = arcs[b].to;
        f.twin = 2 * b + 1;
        g.twin = 2 * b;
        halves[2 * a + 1].origin = y;
        f.face = halves[2 * a].face;
        g.face = halves[2 * a + 1].face;
        f.next = f.twin;
        g.next = g.twin;

        // Forwards, a becomes (a, b); backwards, (b', a').
        for (int which = 0; which < 2; ++which) {
            const int h = 2 * a + which;
            const int fid = halves[h].face;
            if (fid < 0) continue;
            Arrangement::Face &fc = faces[fid];
            for (size_t k = 0; k < fc.half.size(); ++k) {
                if (fc.half[k] != h) continue;
                if (which == 0) {
                    fc.half.insert(fc.half.begin() + static_cast<long>(k) + 1, 2 * b);
                    fc.turns.insert(fc.turns.begin() + static_cast<long>(k), 2);
                } else {
                    fc.half.insert(fc.half.begin() + static_cast<long>(k), 2 * b + 1);
                    fc.turns.insert(fc.turns.begin() + static_cast<long>(k), 2);
                }
                break;
            }
            refresh(fc);
        }
        return y;
    };

    // Where the isoline through a point at fraction t along side k of face f
    // meets the opposite side: that side's half-edge and how far along it.
    struct Exit {
        bool ok = false;
        int half = -1;        // the half-edge of the opposite side it lands on
        double into = 0.0;    // arc length along that half-edge
        int node = -1;        // an existing node it lands on instead, or -1
    };
    auto exitOf = [&](int f, int k, double t) {
        Exit e;
        const Arrangement::Face &F = faces[f];
        const std::vector<int> &S2 = F.sides[(k + 2) % 4];
        double total2 = 0.0;
        for (int h : S2) total2 += halfLength(h);
        const double target = (1.0 - t) * total2;
        // Never onto a corner, or so near one that the face it leaves is a
        // sliver: that is a layout this pass cannot improve on.
        const double margin = 0.05 * total2;
        if (!(total2 > 0.0) || target < margin || target > total2 - margin) return e;
        double acc = 0.0;
        for (size_t j = 0; j < S2.size(); ++j) {
            const double l = halfLength(S2[j]);
            if (target <= acc + l || j + 1 == S2.size()) {
                const double into = target - acc;
                if (j + 1 < S2.size() && l - into < margin) { e.node = endOf(S2[j]); e.ok = true; return e; }
                if (j > 0 && into < margin) { e.node = halves[S2[j]].origin; e.ok = true; return e; }
                e.half = S2[j];
                e.into = std::max(0.0, std::min(l, into));
                e.ok = true;
                return e;
            }
            acc += l;
        }
        return e;
    };
    // Where a point at arc length `into` along half-edge h is on the side of
    // the face that owns h: which side, and the fraction along it.
    auto entryOf = [&](int h, double into, int &side, double &t) {
        const int f = halves[h].face;
        if (f < 0) return false;
        const Arrangement::Face &F = faces[f];
        if (!F.patch || !F.quad || F.sides.size() != 4) return false;
        for (int k = 0; k < 4; ++k) {
            double before = 0.0, total = 0.0;
            bool found = false;
            for (int g : F.sides[k]) {
                if (g == h) { found = true; before = total + into; }
                total += halfLength(g);
            }
            if (found && total > 0.0) {
                side = k;
                t = before / total;
                return true;
            }
        }
        return false;
    };

    // Follow the isoline from node x on side k of face f, without changing
    // anything, and say whether it ends -- on dS, or on a node already in the
    // middle of a side -- within a few faces and without coming back to one.
    // A layout whose chord closes on itself round a ring of faces has no end
    // for it, and carrying it round would only lay more curves on the ring.
    auto reaches = [&](int f, int k, double t) {
        std::set<int> seen;
        for (int hop = 0; hop < 24; ++hop) {
            if (!seen.insert(f).second) return false;
            const Exit e = exitOf(f, k, t);
            if (!e.ok) return false;
            if (e.node >= 0) return true;
            const int a = halves[e.half].arc;
            if (arcs[a].kind == AK::Boundary) return true;
            const int twin = halves[e.half].twin;
            const double len = halfLength(e.half);
            int nk = -1;
            double nt = 0.0;
            if (!entryOf(twin, len - e.into, nk, nt)) return false;
            f = halves[twin].face;
            k = nk;
            t = nt;
        }
        return false;
    };

    // One face split: the isoline from node x (on side k of face f) to the
    // opposite side. Returns the node it lands on, and whether that node is
    // new -- a T-node on the face beyond, for the next split to take up.
    auto splitFace = [&](int f, int k, int x, bool &fresh) {
        fresh = false;
        const int k2 = (k + 2) % 4;
        double before = 0.0, total = 0.0;
        bool atX = false;
        for (int h : faces[f].sides[k]) {
            const double l = halfLength(h);
            total += l;
            if (!atX) before += l;
            if (endOf(h) == x) atX = true;
        }
        const double t = total > 0.0 ? before / total : 0.5;
        const Exit e = exitOf(f, k, t);
        if (!e.ok) return -1;
        int y = e.node;
        if (y < 0) {
            const int a = halves[e.half].arc;
            const double fromStart = halves[e.half].forward ? e.into
                                                            : arcs[a].modelLength - e.into;
            y = splitArc(a, fromStart);
            fresh = arcs[a].kind != AK::Boundary;
        }
        if (y == x) return -1;

        // The curve across the face: its Coons isoline at the parameter of x.
        const std::vector<std::vector<int>> sides = faces[f].sides;
        if (sides.size() != 4) return -1;
        const std::vector<Point> B = sidePolyline(sides[k]);
        const std::vector<Point> R = sidePolyline(sides[(k + 1) % 4]);
        const std::vector<Point> T = sidePolyline(sides[k2]);
        const std::vector<Point> L = sidePolyline(sides[(k + 3) % 4]);
        if (B.empty() || R.empty() || T.empty() || L.empty()) return -1;
        const std::vector<double> cB = cumulativeLength(B), cR = cumulativeLength(R),
                                  cT = cumulativeLength(T), cL = cumulativeLength(L);
        const Point c0 = B.front(), c1 = B.back(), c2 = T.front(), c3 = T.back();
        const Point Bt = atLength(B, cB, t * cB.back());
        const Point Tt = atLength(T, cT, (1.0 - t) * cT.back());
        const int samples = std::max<int>(8, static_cast<int>(std::max(R.size(), L.size())));
        Arrangement::Arc across;
        across.kind = AK::Separatrix;
        across.from = x;
        across.to = y;
        for (int q = 0; q <= samples; ++q) {
            const double v = static_cast<double>(q) / samples;
            const Point Rv = atLength(R, cR, v * cR.back());
            const Point Lv = atLength(L, cL, (1.0 - v) * cL.back());
            across.points.push_back(Bt * (1.0 - v) + Tt * v + Lv * (1.0 - t) + Rv * t -
                                    (c0 * ((1.0 - t) * (1.0 - v)) + c1 * (t * (1.0 - v)) +
                                     c2 * (t * v) + c3 * ((1.0 - t) * v)));
        }
        across.points.front() = nodes[x].p;
        across.points.back() = nodes[y].p;
        across.modelLength = polylineLength(across.points);
        const int n = static_cast<int>(arcs.size());
        arcs.push_back(across);
        halves.resize(2 * arcs.size());
        Arrangement::HalfEdge &fw = halves[2 * n];
        Arrangement::HalfEdge &bw = halves[2 * n + 1];
        fw = Arrangement::HalfEdge();
        bw = Arrangement::HalfEdge();
        fw.arc = bw.arc = n;
        fw.forward = true;
        bw.forward = false;
        fw.origin = x;
        bw.origin = y;
        fw.twin = 2 * n + 1;
        bw.twin = 2 * n;

        // Split the face along it: the part through the two corners after x,
        // and the part through the two before it.
        const Arrangement::Face F = faces[f];
        const size_t m = F.half.size();
        size_t ix = m, iy = m;
        for (size_t i = 0; i < m; ++i) {
            if (endOf(F.half[i]) == x) { ix = i; break; }
        }
        for (size_t d = 1; d <= m && ix < m; ++d) {
            const size_t i = (ix + d) % m;
            if (endOf(F.half[i]) == y) { iy = i; break; }
        }
        if (ix == m || iy == m) return -1;
        Arrangement::Face after = F, beforeFace = F;
        after.half.clear();
        after.turns.clear();
        for (size_t d = 1;; ++d) {
            const size_t i = (ix + d) % m;
            after.half.push_back(F.half[i]);
            after.turns.push_back(i == iy ? 1 : F.turns[i]);
            if (i == iy) break;
        }
        after.half.push_back(2 * n + 1);   // y -> x
        after.turns.push_back(1);
        beforeFace.half.clear();
        beforeFace.turns.clear();
        for (size_t d = 1;; ++d) {
            const size_t i = (iy + d) % m;
            beforeFace.half.push_back(F.half[i]);
            beforeFace.turns.push_back(i == ix ? 1 : F.turns[i]);
            if (i == ix) break;
        }
        beforeFace.half.push_back(2 * n);  // x -> y
        beforeFace.turns.push_back(1);
        faces[f] = beforeFace;
        const int g = static_cast<int>(faces.size());
        faces.push_back(after);
        for (int h : faces[f].half) halves[h].face = f;
        for (int h : faces[g].half) halves[h].face = g;
        refresh(faces[f]);
        refresh(faces[g]);
        return y;
    };

    int splits = 0;
    std::set<int> refused;   // T-nodes whose isoline has no end
    for (size_t guard = 0; guard < 4 * faces.size() + 64; ++guard) {
        int f = -1, k = -1, x = -1;
        double t = 0.0;
        for (size_t i = 0; i < faces.size() && f < 0; ++i) {
            const Arrangement::Face &F = faces[i];
            if (!F.patch || !F.quad || F.sides.size() != 4) continue;
            for (int j = 0; j < 4 && f < 0; ++j) {
                const std::vector<int> &sd = F.sides[j];
                double before = 0.0, total = 0.0;
                for (int h : sd) total += halfLength(h);
                for (size_t q = 0; q + 1 < sd.size(); ++q) {
                    before += halfLength(sd[q]);
                    const int node = endOf(sd[q]);
                    if (refused.count(node)) continue;
                    f = static_cast<int>(i);
                    k = j;
                    x = node;
                    t = total > 0.0 ? before / total : 0.5;
                    break;
                }
            }
        }
        if (f < 0) break;
        if (!reaches(f, k, t)) {
            refused.insert(x);
            continue;
        }
        // Carry it on, face by face, to where reaches() saw it end.
        bool fresh = true;
        for (int hop = 0; hop < 24 && fresh; ++hop) {
            const int y = splitFace(f, k, x, fresh);
            if (y < 0) {
                refused.insert(x);
                break;
            }
            ++splits;
            if (!fresh) break;
            // The face beyond: the one that now has y in the middle of a side.
            int next = -1, side = -1;
            for (size_t i = 0; i < faces.size() && next < 0; ++i) {
                const Arrangement::Face &G = faces[i];
                if (!G.patch || !G.quad || G.sides.size() != 4) continue;
                for (int j = 0; j < 4 && next < 0; ++j) {
                    const std::vector<int> &sd = G.sides[j];
                    for (size_t q = 0; q + 1 < sd.size(); ++q) {
                        if (endOf(sd[q]) == y) { next = static_cast<int>(i); side = j; break; }
                    }
                }
            }
            if (next < 0) break;
            f = next;
            k = side;
            x = y;
        }
    }

    report.tJunctionsLeft = 0;
    for (const Arrangement::Face &F : faces) {
        if (!F.patch) continue;
        for (const std::vector<int> &sd : F.sides) {
            if (sd.size() > 1) report.tJunctionsLeft += static_cast<int>(sd.size()) - 1;
        }
    }
    return splits;
}
