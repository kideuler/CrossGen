#include "tracing/PartitionSimplify.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace {

double polylineLength(const std::vector<Point> &p) {
    double L = 0.0;
    for (size_t i = 1; i < p.size(); ++i) L += normP(p[i] - p[i - 1]);
    return L;
}

// The point a fraction `s` of the way along a polyline by arc length.
Point pointAlong(const std::vector<Point> &p, double s) {
    if (p.empty()) return Point{0.0, 0.0};
    if (p.size() == 1) return p[0];
    const double total = polylineLength(p);
    if (total <= 0.0) return p[0];
    double want = std::min(std::max(s, 0.0), 1.0) * total;
    for (size_t i = 1; i < p.size(); ++i) {
        const double seg = normP(p[i] - p[i - 1]);
        if (want <= seg || i + 1 == p.size()) {
            const double t = (seg > 0.0) ? std::min(1.0, want / seg) : 0.0;
            return p[i - 1] * (1.0 - t) + p[i] * t;
        }
        want -= seg;
    }
    return p.back();
}

// A curve of the model -- dS or a material interface -- rather than a curve
// the tracing chose: never deleted, never moved, never blended into anything.
inline bool fixedArc(const QuadLayout::Arc &a) { return a.onBoundary || a.onInterface; }

// The total length of the interface arcs, which a collapse has to leave
// exactly as it found it. dS is protected by the area check; an interface is
// interior, so moving one changes which material is where and not the total
// area, and this is the check that sees it.
double interfaceLength(const QuadLayout &L) {
    double s = 0.0;
    for (const auto &a : L.getArcs()) if (a.onInterface) s += a.length;
    return s;
}

// Resample a polyline at the given fractions of its arc length.
std::vector<Point> sampleAt(const std::vector<Point> &p, const std::vector<double> &s) {
    std::vector<Point> out;
    out.reserve(s.size());
    for (const double f : s) out.push_back(pointAlong(p, f));
    return out;
}

// 3t^2 - 2t^3 on [0, 1]: from 0 to 1 with zero slope at both ends, so a curve
// eased from one position to another by it leaves the first and arrives at the
// second in the direction it had there.
inline double smoothstep(double t) {
    t = std::min(std::max(t, 0.0), 1.0);
    return t * t * (3.0 - 2.0 * t);
}

// Cut a polyline at a fraction of its arc length into the run up to the cut
// and the run from it, both carrying the cut point.
void splitAt(const std::vector<Point> &p, double f, std::vector<Point> &head,
             std::vector<Point> &tail) {
    head.clear();
    tail.clear();
    if (p.size() < 2) { head = p; tail = p; return; }
    const double total = polylineLength(p);
    double want = std::min(std::max(f, 0.0), 1.0) * total;
    head.push_back(p[0]);
    size_t i = 1;
    for (; i < p.size(); ++i) {
        const double seg = normP(p[i] - p[i - 1]);
        if (want <= seg || i + 1 == p.size()) break;
        want -= seg;
        head.push_back(p[i]);
    }
    const double seg = normP(p[i] - p[i - 1]);
    const double t = (seg > 0.0) ? std::min(1.0, want / seg) : 0.0;
    const Point z = p[i - 1] * (1.0 - t) + p[i] * t;
    head.push_back(z);
    tail.push_back(z);
    for (size_t k = i; k < p.size(); ++k) tail.push_back(p[k]);
}

} // namespace

const char *PartitionSimplify::blockName(Block b) {
    switch (b) {
        case Block::None: return "none";
        case Block::RungNotContractible: return "a node sits along a rung";
        case Block::RungJoinsSingularities: return "rung joins two singularities";
        case Block::RungSingularityToBoundary: return "rung joins a singularity to the boundary";
        case Block::StripBetweenBoundaries: return "strip between two boundaries";
        case Block::ZipAgainstBoundary: return "zip against the boundary";
        case Block::WouldDeleteBoundary: return "would delete the boundary";
        case Block::TJunctionHasNowhereToGo: return "T-junction has nowhere to go";
        case Block::NoPatches: return "no patches";
        case Block::Energy: return "energy (Sec. 4.1)";
        case Block::Drag: return "wider than maxDrag";
        case Block::RungAcrossInterface: return "rung between an interface and another curve";
    }
    return "?";
}

PartitionSimplify::PartitionSimplify(const QuadLayout &layout, const Settings &settings)
    : layout_(1e-9), settings_(settings) {
    // A copy is taken and worked on: the caller's layout is the input and stays
    // the input, so the two can be drawn side by side.
    layout_.rebuild(layout.getNodes(), layout.getArcs());

    double total = 0.0;
    int n = 0;
    for (const auto &a : layout_.getArcs()) { total += a.length; ++n; }
    tol_ = (n > 0) ? 1e-6 * total / n : 1e-9;

    // The scale Settings::maxDrag is measured in. Taken from the model, not
    // from the layout: the layout's own arcs get longer as it is simplified,
    // and the question the cap asks -- is this strip discretisation error --
    // is about the mesh.
    if (const Mesh *m = layout.getMesh()) {
        double e = 0.0;
        for (const auto &ed : m->edges) e += normP(m->vertices[ed[1]] - m->vertices[ed[0]]);
        if (!m->edges.empty()) meshEdge_ = e / static_cast<double>(m->edges.size());
    }

    if (const SeparatrixTrace *t = layout.getTrace()) tracer_ = &t->getTracer();

    const auto &r = layout_.getReport();
    report_.componentsBefore = r.faces;
    report_.tJunctionsBefore = r.tJunctions;
    report_.quadsBefore = r.quadFaces;
}

// ---------------------------------------------------------------------------
// Small queries
// ---------------------------------------------------------------------------
bool PartitionSimplify::isSingularity(int node) const {
    return node >= 0 && layout_.getNodes()[node].kind == QuadLayout::NodeKind::Singularity;
}

// docs/viertel_2019.md Sec. 6: "fixed nodes = SING u BSING". A corner of the
// model that is not flat emits separatrices exactly as a singularity does
// (q - 1 of them), and so does a node of the interface network; a strip whose
// one side is such a separatrix is a strip that side has to survive, and a
// rung carrying one is a patch boundary. Counting only interior singularities
// made the side a non-zip keeps an arbitrary one wherever a corner, not a
// singularity, was what made the side special -- harmless while every strip
// collapsed was a sliver, but a wider one drags the corner's own separatrix
// away and leaves the corner a block of its own with a straight angle in it.
bool PartitionSimplify::isFixed(int node) const {
    if (isSingularity(node)) return true;
    if (!settings_.fixedCorners || node < 0) return false;
    const auto k = layout_.getNodes()[node].kind;
    return k == QuadLayout::NodeKind::BoundaryCorner || k == QuadLayout::NodeKind::InterfaceNode;
}

// On dS or on a material interface: both are curves of the model, and a node
// on one may not be pulled off it.
bool PartitionSimplify::isOnBoundary(int node) const {
    if (node < 0) return false;
    for (const int d : layout_.getNodes()[node].darts)
        if (fixedArc(layout_.getArcs()[d >> 1])) return true;
    return false;
}

bool PartitionSimplify::isOnInterface(int node) const {
    if (node < 0) return false;
    for (const int d : layout_.getNodes()[node].darts)
        if (layout_.getArcs()[d >> 1].onInterface) return true;
    return false;
}

// A node of the interface network: a corner of every region it touches, and
// like a singularity never merged into anything.
bool PartitionSimplify::isNetworkNode(int node) const {
    return node >= 0 && layout_.getNodes()[node].kind == QuadLayout::NodeKind::InterfaceNode;
}

// At a T-junction three arcs meet: two carry the side that runs past it and
// one is the separatrix that stopped there. The pair carrying the side are the
// two bounding the widest sector, so the stem is the third.
int PartitionSimplify::stemArcOf(int node) const {
    const auto &nd = layout_.getNodes()[node];
    if (nd.darts.size() != 3) return -1;
    int widest = 0;
    double best = -1.0;
    for (size_t k = 0; k < 3; ++k) {
        double w = nd.angles[(k + 1) % 3] - nd.angles[k];
        w = std::fmod(w, 2.0 * M_PI);
        if (w < 0.0) w += 2.0 * M_PI;
        if (w > best) { best = w; widest = static_cast<int>(k); }
    }
    // The widest sector runs from darts[widest] to darts[widest+1]; the stem is
    // the one that is neither.
    return nd.darts[(widest + 2) % 3] >> 1;
}

void PartitionSimplify::indexDarts() {
    faceOfDart_.assign(2 * layout_.getArcs().size(), -1);
    for (size_t f = 0; f < layout_.getFaces().size(); ++f)
        for (const int d : layout_.getFaces()[f].darts) faceOfDart_[d] = static_cast<int>(f);
}

// ---------------------------------------------------------------------------
// sharedSide()  --  can a chord carry on through this side?
//
// Only if the whole side is one side of the component on the other side of it.
// A T-junction anywhere along it means the two components disagree about where
// the side starts and ends, and the chord stops there (Sec. 4: "no T-junction
// exists between c_i and c_{i+1}").
// ---------------------------------------------------------------------------
bool PartitionSimplify::sharedSide(int f, int s, int &g, int &t) const {
    const auto &face = layout_.getFaces()[f];
    if (s < 0 || s >= static_cast<int>(face.sides.size())) return false;
    const auto &side = face.sides[s];
    if (side.empty()) return false;

    g = faceOfDart_[side[0] ^ 1];
    if (g < 0) return false;   // the far side is the unbounded face
    for (const int d : side)
        if (faceOfDart_[d ^ 1] != g) return false;

    const auto &other = layout_.getFaces()[g];
    if (other.sides.size() != 4) return false;
    for (size_t k = 0; k < other.sides.size(); ++k) {
        if (other.sides[k].size() != side.size()) continue;
        bool same = true;
        for (size_t i = 0; i < side.size(); ++i)
            if (other.sides[k][i] != (side[side.size() - 1 - i] ^ 1)) { same = false; break; }
        if (same) { t = static_cast<int>(k); return true; }
    }
    return false;
}

// ---------------------------------------------------------------------------
// walkChord()
//
// Sides are numbered in the order the darts go round the component, so entering
// by side e means leaving by side e+2 and the two longitudinal sides are e+1
// and e+3. That labelling carries along the chord: the component reached
// through e+2 is entered by the side matching it, and its e+1 continues the
// same longitudinal side. So "left" and "right" stay consistent without
// anything geometric being measured.
// ---------------------------------------------------------------------------
bool PartitionSimplify::walkChord(int seedFace, int seedSide, Chord &out) const {
    const auto &faces = layout_.getFaces();
    if (faces[seedFace].sides.size() != 4) return false;

    out = Chord{};
    out.faces.push_back(seedFace);
    out.entry.push_back(seedSide);

    std::vector<char> onChord(faces.size(), 0);
    onChord[seedFace] = 1;

    // Forward, leaving by the side opposite the one it came in by.
    int cur = seedFace, e = seedSide;
    for (;;) {
        int g = -1, t = -1;
        if (!sharedSide(cur, (e + 2) % 4, g, t)) break;
        if (g == seedFace) { out.cyclic = true; break; }
        if (onChord[g]) break;   // meets itself other than head to tail
        onChord[g] = 1;
        out.faces.push_back(g);
        out.entry.push_back(t);
        cur = g;
        e = t;
    }

    // Backward, in through the side it came in by.
    if (!out.cyclic) {
        cur = seedFace;
        e = seedSide;
        for (;;) {
            int h = -1, u = -1;
            if (!sharedSide(cur, e, h, u)) break;
            if (h == seedFace) { out.cyclic = true; break; }
            if (onChord[h]) break;
            onChord[h] = 1;
            // h is left through side u, so it is entered through u + 2.
            const int hEntry = (u + 2) % 4;
            out.faces.insert(out.faces.begin(), h);
            out.entry.insert(out.entry.begin(), hEntry);
            cur = h;
            e = hEntry;
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// analyse()  --  rungs, longitudinal sides, patches, energy
// ---------------------------------------------------------------------------
void PartitionSimplify::analyse(Chord &c) const {
    const auto &faces = layout_.getFaces();
    const auto &arcs = layout_.getArcs();

    auto startNodeOfSide = [&](int f, int s) {
        const auto &face = faces[f];
        const int d = face.sides[s][0];
        for (size_t i = 0; i < face.darts.size(); ++i)
            if (face.darts[i] == d) return face.nodes[i];
        return -1;
    };
    auto dartsOfSide = [&](int f, int s) { return faces[f].sides[s]; };
    auto makeRung = [&](int f, int s, int endL, int endR) {
        Rung r;
        r.darts = dartsOfSide(f, s);
        r.endL = endL;
        r.endR = endR;
        for (const int d : r.darts) {
            r.length += arcs[d >> 1].length;
            if (fixedArc(arcs[d >> 1])) r.onBoundary = true;
            if (arcs[d >> 1].onInterface) r.onInterface = true;
        }
        // Contracting a rung merges its two ends into one point, so anything
        // sitting between them is squashed too. A plain join is nothing to
        // squash; a T-junction along the rung is a separatrix whose end would
        // be dragged onto the merged node, so those rungs are left alone.
        for (size_t k = 1; k < r.darts.size(); ++k) {
            const int between = (r.darts[k] & 1) ? arcs[r.darts[k] >> 1].b : arcs[r.darts[k] >> 1].a;
            if (layout_.getNodes()[between].darts.size() > 2) r.contractible = false;
        }
        return r;
    };

    const int n = static_cast<int>(c.faces.size());
    c.rungs.clear();
    c.sideL.assign(n, {});
    c.sideR.assign(n, {});

    // The four corners of component i, in the order the darts run:
    //   A = start of the entry side, B of the left side, C of the exit, D of the right.
    // So the entry rung runs A (right end) to B (left end) and the exit rung
    // runs C (left end) to D (right end).
    for (int i = 0; i < n; ++i) {
        const int f = c.faces[i], e = c.entry[i];
        c.sideL[i] = dartsOfSide(f, (e + 1) % 4);
        c.sideR[i] = dartsOfSide(f, (e + 3) % 4);
    }

    if (!c.cyclic) {
        const int f0 = c.faces[0], e0 = c.entry[0];
        c.rungs.push_back(makeRung(f0, e0, startNodeOfSide(f0, (e0 + 1) % 4),
                                   startNodeOfSide(f0, e0)));
    }
    for (int i = 0; i < n; ++i) {
        const int f = c.faces[i], e = c.entry[i];
        c.rungs.push_back(makeRung(f, (e + 2) % 4, startNodeOfSide(f, (e + 2) % 4),
                                   startNodeOfSide(f, (e + 3) % 4)));
    }

    c.minWidth = std::numeric_limits<double>::max();
    c.maxWidth = 0.0;
    for (const auto &r : c.rungs) {
        c.minWidth = std::min(c.minWidth, r.length);
        c.maxWidth = std::max(c.maxWidth, r.length);
    }

    // A T-junction at either end of the chord is one the collapse would put on
    // top of whatever the rung's other end is, which is how a collapse gets rid
    // of one; see Sec. 4 condition 3 and Settings::tJunctionsFirst.
    c.tJunctionEnds = 0;
    if (!c.rungs.empty() && !c.cyclic) {
        const auto &nodes = layout_.getNodes();
        for (const int i : {0, static_cast<int>(c.rungs.size()) - 1})
            for (const int e : {c.rungs[i].endL, c.rungs[i].endR})
                if (e >= 0 && nodes[e].kind == QuadLayout::NodeKind::TJunction) ++c.tJunctionEnds;
    }

    auto patches = patchesOf(c);
    c.collapsible = !patches.empty();
    if (patches.empty()) c.block = Block::NoPatches;
    c.energy = std::numeric_limits<double>::max();
    for (auto &p : patches) {
        if (!patchOk(c, p, c.block)) { c.collapsible = false; break; }
        p.energy = patchEnergy(c, p);
        c.energy = std::min(c.energy, p.energy);
    }
    if (!c.collapsible) c.energy = 0.0;
}

// A patch is a maximal run of components with singularities only on its first
// and last rung, so the cuts are exactly the rungs carrying one.
std::vector<PartitionSimplify::Patch> PartitionSimplify::patchesOf(const Chord &c) const {
    std::vector<Patch> out;
    const int nr = static_cast<int>(c.rungs.size());
    if (nr < 2) return out;

    auto singular = [&](int i) {
        return isFixed(c.rungs[i].endL) || isFixed(c.rungs[i].endR);
    };

    if (c.cyclic) {
        std::vector<int> cuts;
        for (int i = 0; i < nr; ++i) if (singular(i)) cuts.push_back(i);
        if (cuts.empty()) {
            Patch p;
            p.first = 0;
            p.last = nr;   // the whole ring; rung 0 is also rung nr
            out.push_back(p);
            return out;
        }
        for (size_t k = 0; k < cuts.size(); ++k) {
            Patch p;
            p.first = cuts[k];
            p.last = cuts[(k + 1) % cuts.size()];
            if (p.last <= p.first) p.last += nr;
            out.push_back(p);
        }
        return out;
    }

    int start = 0;
    for (int i = 1; i < nr; ++i) {
        if (!singular(i) && i != nr - 1) continue;
        Patch p;
        p.first = start;
        p.last = i;
        out.push_back(p);
        start = i;
    }
    return out;
}

// ---------------------------------------------------------------------------
// patchOk()  --  Sec. 4's three conditions, plus what the boundary allows
// ---------------------------------------------------------------------------
bool PartitionSimplify::patchOk(const Chord &c, Patch &p, Block &why) const {
    const int nr = static_cast<int>(c.rungs.size());
    auto rung = [&](int i) -> const Rung & { return c.rungs[i % nr]; };

    for (int i = p.first; i <= p.last; ++i) {
        const Rung &r = rung(i);
        if (!r.contractible && !settings_.contractRungsWithNodes) {
            why = Block::RungNotContractible;
            return false;
        }
        // 1. Contracting a rung merges its two ends, so two singularities on
        //    one rung would have to become one.
        if (isFixed(r.endL) && isFixed(r.endR)) { why = Block::RungJoinsSingularities; return false; }
        // 2. And a singularity on a rung whose other end is on the boundary
        //    would have to move onto the boundary.
        if ((isSingularity(r.endL) && isOnBoundary(r.endR)) ||
            (isSingularity(r.endR) && isOnBoundary(r.endL))) {
            why = Block::RungSingularityToBoundary;
            return false;
        }
        // A node of the interface network is fixed as a singularity is: it
        // may absorb the other end of a rung but never be absorbed, so a rung
        // with one at each end -- or one and a corner of dS -- cannot go.
        auto fixedKind = [&](int n) {
            return isSingularity(n) || isNetworkNode(n) ||
                   (n >= 0 && layout_.getNodes()[n].kind == QuadLayout::NodeKind::BoundaryCorner);
        };
        if ((isNetworkNode(r.endL) && fixedKind(r.endR)) ||
            (isNetworkNode(r.endR) && fixedKind(r.endL))) {
            why = Block::RungJoinsSingularities;
            return false;
        }
        // A rung that is not itself a piece of an interface, joining a node
        // on an interface to a node on another curve of the model: whichever
        // end went, it would leave its curve. dS alone has the area check to
        // catch that; an interface has only this.
        if (!r.onBoundary && isOnBoundary(r.endL) && isOnBoundary(r.endR) &&
            (isOnInterface(r.endL) || isOnInterface(r.endR))) {
            why = Block::RungAcrossInterface;
            return false;
        }
    }

    // Which side carries the fixed nodes decides zip against non-zip.
    const Rung &rp = rung(p.first);
    const Rung &rq = rung(p.last);
    const bool sLp = isFixed(rp.endL), sRp = isFixed(rp.endR);
    const bool sLq = isFixed(rq.endL), sRq = isFixed(rq.endR);
    p.zip = (sLp && sRq) || (sRp && sLq);
    p.keepL = sLp || sLq;

    // The boundary is never the side that gets deleted, and never gets blended
    // into something else either.
    bool boundaryL = false, boundaryR = false;
    const auto &arcs = layout_.getArcs();
    for (int i = p.first; i < p.last; ++i) {
        const int comp = c.cyclic ? (i % static_cast<int>(c.faces.size()))
                                  : i;   // component i lies between rungs i and i+1
        if (comp < 0 || comp >= static_cast<int>(c.sideL.size())) continue;
        for (const int d : c.sideL[comp]) if (fixedArc(arcs[d >> 1])) boundaryL = true;
        for (const int d : c.sideR[comp]) if (fixedArc(arcs[d >> 1])) boundaryR = true;
    }
    if (boundaryL && boundaryR) { why = Block::StripBetweenBoundaries; return false; }
    if (p.zip && (boundaryL || boundaryR)) { why = Block::ZipAgainstBoundary; return false; }
    if (!p.zip) {
        if (boundaryL && !p.keepL && (sRp || sRq)) { why = Block::WouldDeleteBoundary; return false; }
        if (boundaryR && p.keepL && (sLp || sLq)) { why = Block::WouldDeleteBoundary; return false; }
        if (boundaryL) p.keepL = true;
        if (boundaryR) p.keepL = false;
    }

    // 3. A T-junction at a corner of the patch has to end up somewhere: either
    //    on a singularity, or nose to nose with another T-junction so that the
    //    two become one line through the merged node, or opposite a singularity
    //    across the patch.
    for (const int i : {p.first, p.last}) {
        const Rung &r = rung(i);
        const bool endRung = c.cyclic || (i == 0) || (i == nr - 1);
        if (!endRung) continue;   // an interior cut is a singular rung, handled above
        for (int side = 0; side < 2; ++side) {
            const int T = side ? r.endR : r.endL;
            const int O = side ? r.endL : r.endR;
            if (T < 0 || layout_.getNodes()[T].kind != QuadLayout::NodeKind::TJunction) continue;
            if (isSingularity(O)) continue;                                   // (a)
            const int stem = stemArcOf(T);
            auto rungHas = [&](int arc) {
                for (const int d : r.darts) if ((d >> 1) == arc) return true;
                return false;
            };
            const bool stemIsRung = stem >= 0 && rungHas(stem);
            if (O >= 0 && layout_.getNodes()[O].kind == QuadLayout::NodeKind::TJunction) {
                const int stemO = stemArcOf(O);
                const bool stemOisRung = stemO >= 0 && rungHas(stemO);
                if (stemIsRung && stemOisRung) continue;                       // (b)
            }
            const Rung &far = rung(i == p.first ? p.last : p.first);
            const int diag = side ? far.endL : far.endR;                       // (c)
            if (isSingularity(diag)) continue;
            why = Block::TJunctionHasNowhereToGo;
            return false;
        }
    }

    p.ok = true;
    return true;
}

// ---------------------------------------------------------------------------
// patchEnergy()  --  Sec. 4.1
// ---------------------------------------------------------------------------
double PartitionSimplify::patchEnergy(const Chord &c, const Patch &p) const {
    if (!p.zip && !settings_.collapseNonZip) return -1.0;
    if (!p.zip && !(settings_.nonZipAngle > 0.0)) return 1.0;

    const int nr = static_cast<int>(c.rungs.size());
    // w: the mean length of the patch's rungs.
    double w = 0.0;
    int nw = 0;
    for (int i = p.first; i <= p.last; ++i) { w += c.rungs[i % nr].length; ++nw; }
    if (nw) w /= nw;

    // l: the mean length of the chord's two longitudinal sides. Sec. 4.1 takes
    // this over the whole chord, not the patch, so a short patch on a long
    // chord is judged by how far the chord as a whole has to bend.
    const auto &arcs = layout_.getArcs();
    double lL = 0.0, lR = 0.0;
    for (const auto &side : c.sideL) for (const int d : side) lL += arcs[d >> 1].length;
    for (const auto &side : c.sideR) for (const int d : side) lR += arcs[d >> 1].length;
    const double l = 0.5 * (lL + lR);
    if (l <= 0.0) return -1.0;

    return (p.zip ? settings_.zipAngle : settings_.nonZipAngle) - std::atan(w / l);
}

// ---------------------------------------------------------------------------
// collapse()  --  contract every rung of one chord
//
// The two longitudinal sides become one curve and the strip between them goes.
// Everything else in the layout is untouched except where it was attached to
// one of the two sides: those attachment points move onto the merged curve,
// which is the paper's "the hanging separatrix is simply extended after the
// collapse operation until it crosses the next separatrix" -- the next
// separatrix being the merged curve itself, and the extension being however far
// the strip was wide, which is small or the chord would not have been chosen.
// ---------------------------------------------------------------------------
bool PartitionSimplify::collapse(const Chord &c) {
    const std::vector<QuadLayout::Node> oldNodes = layout_.getNodes();
    const std::vector<QuadLayout::Arc> oldArcs = layout_.getArcs();
    const int nComp = static_cast<int>(c.faces.size());
    const int nRung = static_cast<int>(c.rungs.size());
    if (nComp == 0 || nRung < 2) return false;

    // A run of darts as one oriented polyline, with the nodes along it and
    // where each sits by arc length.
    struct Chain {
        std::vector<int> nodes;
        std::vector<Point> pts;
        std::vector<double> nodeS;  // where each node sits, by arc length
        std::vector<double> ptS;    // and where each point of the polyline does
    };
    auto chainOf = [&](const std::vector<int> &darts) {
        Chain g;
        if (darts.empty()) return g;
        std::vector<double> cum{0.0};
        for (size_t k = 0; k < darts.size(); ++k) {
            const QuadLayout::Arc &a = oldArcs[darts[k] >> 1];
            const bool flip = (darts[k] & 1) != 0;
            const int from = flip ? a.b : a.a;
            const int to = flip ? a.a : a.b;
            if (k == 0) { g.nodes.push_back(from); g.pts.push_back(flip ? a.pts.back() : a.pts.front()); }
            if (flip)
                for (size_t i = a.pts.size() - 1; i-- > 0;) g.pts.push_back(a.pts[i]);
            else
                for (size_t i = 1; i < a.pts.size(); ++i) g.pts.push_back(a.pts[i]);
            g.nodes.push_back(to);
            cum.push_back(polylineLength(g.pts));
        }
        const double total = cum.back();
        g.nodeS.resize(cum.size());
        for (size_t i = 0; i < cum.size(); ++i) g.nodeS[i] = (total > 0.0) ? cum[i] / total : 0.0;
        g.ptS.assign(g.pts.size(), 0.0);
        double run = 0.0;
        for (size_t i = 1; i < g.pts.size(); ++i) {
            run += normP(g.pts[i] - g.pts[i - 1]);
            g.ptS[i] = (total > 0.0) ? run / total : 0.0;
        }
        if (!g.ptS.empty()) g.ptS.back() = 1.0;
        return g;
    };

    std::vector<char> consumed(oldArcs.size(), 0);
    for (const auto &r : c.rungs) for (const int d : r.darts) consumed[d >> 1] = 1;
    for (const auto &s : c.sideL) for (const int d : s) consumed[d >> 1] = 1;
    for (const auto &s : c.sideR) for (const int d : s) consumed[d >> 1] = 1;

    std::vector<int> remap(oldNodes.size());
    for (size_t i = 0; i < remap.size(); ++i) remap[i] = static_cast<int>(i);
    std::vector<Point> pos(oldNodes.size());
    std::vector<QuadLayout::NodeKind> kind(oldNodes.size());
    for (size_t i = 0; i < oldNodes.size(); ++i) { pos[i] = oldNodes[i].pos; kind[i] = oldNodes[i].kind; }

    std::vector<QuadLayout::Arc> newArcs;
    newArcs.reserve(oldArcs.size());

    // Rungs that lie on the boundary keep their geometry: the boundary is a
    // curve of the model, not something the layout may shorten, so the arc
    // carrying on past the dying end swallows the rung instead.
    struct Absorb { int arc; int atNode; int intoNode; std::vector<Point> extra; };
    std::vector<Absorb> absorbs;

    // Which case each patch is -- zip, or which side a non-zip keeps -- is
    // decided by patchOk(), not by patchesOf(), which only finds where the
    // patches start and end. Taking the bare patches here used to collapse
    // every one as a non-zip keeping its right-hand side: a zip was never
    // blended, and wherever the singularity sat on the left the kept curve
    // ran alongside it and was dragged onto it in its last segment.
    auto patches = patchesOf(c);
    for (Patch &p : patches) {
        Block why = Block::None;
        if (!patchOk(c, p, why)) return false;
    }
    for (const Patch &p : patches) {
        // The two sides of this patch, each running from its first rung to its
        // last. sideR is stored against the traversal of its component, which
        // runs the other way, so it is reversed here.
        std::vector<int> dl, dr;
        for (int i = p.first; i < p.last; ++i) {
            const int comp = i % nComp;
            for (const int d : c.sideL[comp]) dl.push_back(d);
        }
        for (int i = p.first; i < p.last; ++i) {
            const int comp = i % nComp;
            for (size_t k = c.sideR[comp].size(); k-- > 0;) dr.push_back(c.sideR[comp][k] ^ 1);
        }
        const Chain L = chainOf(dl);
        const Chain R = chainOf(dr);
        if (L.pts.size() < 2 || R.pts.size() < 2) return false;

        const Rung &rp = c.rungs[p.first % nRung];
        const Rung &rq = c.rungs[p.last % nRung];
        const bool sLp = isFixed(rp.endL), sRq = isFixed(rq.endR);

        // Where the rungs sit along each side, which is what pairs the two.
        std::vector<double> rungSL{0.0}, rungSR{0.0};
        {
            double accL = 0.0, accR = 0.0;
            std::vector<double> lenL, lenR;
            for (int i = p.first; i < p.last; ++i) {
                const int comp = i % nComp;
                double a = 0.0, b = 0.0;
                for (const int d : c.sideL[comp]) a += oldArcs[d >> 1].length;
                for (const int d : c.sideR[comp]) b += oldArcs[d >> 1].length;
                lenL.push_back(a); lenR.push_back(b); accL += a; accR += b;
            }
            double s = 0.0;
            for (double a : lenL) { s += a; rungSL.push_back(accL > 0 ? s / accL : 0.0); }
            s = 0.0;
            for (double b : lenR) { s += b; rungSR.push_back(accR > 0 ? s / accR : 0.0); }
        }

        // The two sides are paired through the rungs, not through their own arc
        // lengths. A common parameter u runs 0, 1, 2, ... one unit per
        // component, so u = k is rung k on both sides at once and the points
        // being averaged are the two ends of a rung rather than two points that
        // merely happen to be the same fraction along curves of different
        // length. Without that the blend of a strip whose two sides differ in
        // length wanders out of the strip and crosses whatever is outside it.
        const int nSeg = p.last - p.first;
        auto atU = [&](const std::vector<double> &rungS, double u) {
            u = std::min(std::max(u, 0.0), static_cast<double>(nSeg));
            const int k = std::min(nSeg - 1, static_cast<int>(std::floor(u)));
            const double t = u - k;
            return rungS[k] * (1.0 - t) + rungS[k + 1] * t;
        };
        auto uOf = [&](const std::vector<double> &rungS, double s) {
            for (int k = 0; k < nSeg; ++k) {
                if (s > rungS[k + 1] + 1e-15) continue;
                const double d = rungS[k + 1] - rungS[k];
                return k + ((d > 1e-15) ? (s - rungS[k]) / d : 0.0);
            }
            return static_cast<double>(nSeg);
        };

        // A zip blends the two sides so that it starts on whichever carries the
        // singularity there and ends on the other; a non-zip is simply the side
        // that survives, point for point.
        //
        // The weight runs with arc length along the patch (docs/viertel_2019.md
        // Sec. 8: "cumulative mean-side length up to r_i"), not with the
        // component count u: the components of a chord are anything but equal
        // in length, and one a fiftieth of the patch long next to a singularity
        // would otherwise take a whole step of the blend, which is a hook in
        // the merged curve right where it leaves the singularity. It is then
        // eased (Settings::smoothCollapse), so that the merged curve leaves
        // each singularity along the separatrix it was traced as and bends only
        // in between: the angles separatrices make at a singularity, which is
        // what Sec. 4.1's energy exists to protect, come out exactly as traced.
        // weightR(u) is how far across from L to R the merged curve is at u.
        auto weightR = [&](double u) {
            double f = (nSeg > 0) ? u / nSeg : 0.0;
            if (settings_.smoothCollapse)
                f = smoothstep(0.5 * (atU(rungSL, u) + atU(rungSR, u)));
            return p.zip ? (sLp && sRq ? f : 1.0 - f) : (p.keepL ? 0.0 : 1.0);
        };
        auto blend = [&](double u) {
            const double alpha = weightR(u);
            if (alpha <= 0.0) return pointAlong(L.pts, atU(rungSL, u));
            if (alpha >= 1.0) return pointAlong(R.pts, atU(rungSR, u));
            return pointAlong(L.pts, atU(rungSL, u)) * (1.0 - alpha) +
                   pointAlong(R.pts, atU(rungSR, u)) * alpha;
        };

        // Everything that has to sit on the merged curve, by parameter: the
        // merged rung nodes, and every node either side carried in between.
        struct Placed { double s; int node; bool isRung; };
        std::vector<Placed> placed;
        for (int i = p.first, k = 0; i <= p.last; ++i, ++k) {
            const Rung &r = c.rungs[i % nRung];
            const double s = k;   // rung k sits at u = k on both sides
            if (r.endL < 0 || r.endR < 0 || r.endL == r.endR) return false;

            // Which of the two survives: nothing may move a singularity, and
            // nothing may pull a boundary node off the boundary.
            int keep = -1;
            if (isFixed(r.endL) || isNetworkNode(r.endL)) keep = r.endL;
            else if (isFixed(r.endR) || isNetworkNode(r.endR)) keep = r.endR;
            else if (isOnBoundary(r.endL) && !isOnBoundary(r.endR)) keep = r.endL;
            else if (isOnBoundary(r.endR) && !isOnBoundary(r.endL)) keep = r.endR;
            else keep = p.keepL ? r.endL : r.endR;
            const int die = (keep == r.endL) ? r.endR : r.endL;

            if (r.onBoundary) {
                // The rung is a piece of the boundary; the merged node stays
                // where the surviving end already is and the boundary arc
                // beyond the dying end takes over the rung's geometry.
                // Of the same kind: an interface piece is carried on by the
                // interface, a boundary piece by dS.
                int other = -1;
                for (const int d : oldNodes[die].darts) {
                    const int arc = d >> 1;
                    if (!fixedArc(oldArcs[arc]) || oldArcs[arc].onInterface != r.onInterface) continue;
                    bool inRung = false;
                    for (const int rd : r.darts) if ((rd >> 1) == arc) inRung = true;
                    if (!inRung) { other = arc; break; }
                }
                if (other < 0) return false;
                Chain rc = chainOf(r.darts);
                if (rc.nodes.front() != die) {
                    std::reverse(rc.pts.begin(), rc.pts.end());
                    std::reverse(rc.nodes.begin(), rc.nodes.end());
                }
                if (rc.nodes.front() != die) return false;

                // Across a zip the merged curve runs between the two sides, so
                // a rung that is a piece of an interface is not contracted onto
                // either end: the merged node goes where the blend crosses the
                // rung, which is a point of the rung and so of the interface,
                // and the interface arcs beyond the two ends take the two
                // halves. The interface is the same curve afterwards, and the
                // merged curve passes through its node instead of hooking out
                // to whichever side of the strip the node would have stayed on.
                int otherKeep = -1;
                if (settings_.smoothCollapse && p.zip && r.onInterface && !isNetworkNode(keep) &&
                    !isSingularity(keep)) {
                    for (const int d : oldNodes[keep].darts) {
                        const int arc = d >> 1;
                        if (arc == other || !oldArcs[arc].onInterface) continue;
                        bool inRung = false;
                        for (const int rd : r.darts) if ((rd >> 1) == arc) inRung = true;
                        if (!inRung) { otherKeep = arc; break; }
                    }
                }
                if (otherKeep >= 0) {
                    // The blend's weight is measured from L; the rung runs from
                    // the dying end to the surviving one.
                    const double a = weightR(s);
                    const double fromDie = (die == r.endL) ? a : 1.0 - a;
                    std::vector<Point> head, tail;
                    splitAt(rc.pts, fromDie, head, tail);   // die..z, z..keep
                    std::reverse(tail.begin(), tail.end()); // keep..z
                    absorbs.push_back(Absorb{other, die, keep, head});
                    absorbs.push_back(Absorb{otherKeep, keep, keep, tail});
                    pos[keep] = head.back();
                } else {
                    absorbs.push_back(Absorb{other, die, keep, rc.pts});
                }
            } else if (!isSingularity(keep) && !isOnBoundary(keep)) {
                pos[keep] = blend(s);
            }

            if (kind[keep] != QuadLayout::NodeKind::InterfaceNode &&
                static_cast<int>(kind[keep]) > static_cast<int>(kind[die]))
                kind[keep] = kind[die];
            remap[die] = keep;
            placed.push_back(Placed{s, keep, true});
        }

        // The nodes carried along either side between the rungs come with it.
        auto carry = [&](const Chain &ch, const std::vector<double> &rungS, bool isL) {
            for (size_t i = 1; i + 1 < ch.nodes.size(); ++i) {
                const int nd = ch.nodes[i];
                bool atRung = false;
                for (const double rs : rungS) if (std::fabs(ch.nodeS[i] - rs) < 1e-12) atRung = true;
                if (atRung) continue;
                const double s = uOf(rungS, ch.nodeS[i]);
                if (!isSingularity(nd) && !isOnBoundary(nd)) pos[nd] = blend(s);
                placed.push_back(Placed{s, nd, false});
                (void)isL;
            }
        };
        carry(L, rungSL, true);
        carry(R, rungSR, false);

        std::sort(placed.begin(), placed.end(),
                  [](const Placed &a, const Placed &b) { return a.s < b.s; });
        placed.erase(std::unique(placed.begin(), placed.end(),
                                 [](const Placed &a, const Placed &b) { return a.node == b.node; }),
                     placed.end());

        // Cut the merged curve at those points and emit one arc between each
        // consecutive pair.
        auto keptSideHas = [&](bool interface) {
            for (int i = p.first; i < p.last; ++i) {
                const int comp = i % nComp;
                const auto &side = p.keepL ? c.sideL[comp] : c.sideR[comp];
                for (const int d : side)
                    if (interface ? oldArcs[d >> 1].onInterface : oldArcs[d >> 1].onBoundary) return true;
            }
            return false;
        };
        const bool mergedIsBoundary = keptSideHas(false);
        const bool mergedIsInterface = keptSideHas(true);
        // The merged curve is followed at the resolution the two sides already
        // have, not resampled at some fixed count. A side kept whole is emitted
        // point for point: resampling it coarsely would cut the corners off a
        // curve the layout is already committed to, and a cut corner is an arc
        // that crosses its neighbour with no node between them.
        std::vector<double> knots;
        for (const double v : L.ptS) knots.push_back(uOf(rungSL, v));
        for (const double v : R.ptS) knots.push_back(uOf(rungSR, v));
        std::sort(knots.begin(), knots.end());
        knots.erase(std::unique(knots.begin(), knots.end(),
                                [](double a, double b) { return std::fabs(a - b) < 1e-12; }),
                    knots.end());

        for (size_t i = 1; i < placed.size(); ++i) {
            const double s0 = placed[i - 1].s, s1 = placed[i].s;
            std::vector<Point> pts;
            pts.push_back(blend(s0));
            for (const double v : knots)
                if (v > s0 + 1e-12 && v < s1 - 1e-12) pts.push_back(blend(v));
            pts.push_back(blend(s1));
            QuadLayout::Arc a;
            a.pts = std::move(pts);
            a.a = placed[i - 1].node;
            a.b = placed[i].node;
            a.onBoundary = mergedIsBoundary;
            a.onInterface = mergedIsInterface;
            a.separatrix = -1;
            newArcs.push_back(std::move(a));
        }
    }

    // Everything the chord did not consume, carried over.
    for (size_t i = 0; i < oldArcs.size(); ++i) {
        if (consumed[i]) continue;
        QuadLayout::Arc a = oldArcs[i];
        for (const Absorb &ab : absorbs) {
            if (ab.arc != static_cast<int>(i)) continue;
            std::vector<Point> extra = ab.extra;
            if (a.a == ab.atNode) {
                std::reverse(extra.begin(), extra.end());
                extra.pop_back();
                extra.insert(extra.end(), a.pts.begin(), a.pts.end());
                a.pts = std::move(extra);
                a.a = ab.intoNode;
            } else if (a.b == ab.atNode) {
                a.pts.pop_back();
                a.pts.insert(a.pts.end(), extra.begin(), extra.end());
                a.b = ab.intoNode;
            }
        }
        newArcs.push_back(std::move(a));
    }

    // Follow the merges, then move every endpoint onto the node it belongs to.
    auto resolve = [&](int n) {
        int guard = 0;
        while (remap[n] != n && guard++ < 64) n = remap[n];
        return n;
    };
    std::vector<QuadLayout::Node> newNodes(oldNodes.size());
    for (size_t i = 0; i < oldNodes.size(); ++i) {
        newNodes[i] = oldNodes[i];
        newNodes[i].pos = pos[i];
        newNodes[i].kind = kind[i];
        newNodes[i].darts.clear();
        newNodes[i].angles.clear();
    }
    std::vector<char> used(newNodes.size(), 0);
    std::vector<QuadLayout::Arc> keptArcs;
    keptArcs.reserve(newArcs.size());
    for (auto &a : newArcs) {
        a.a = resolve(a.a);
        a.b = resolve(a.b);
        if (a.a < 0 || a.b < 0 || a.pts.size() < 2) continue;
        a.pts.front() = newNodes[a.a].pos;
        a.pts.back() = newNodes[a.b].pos;
        if (a.a == a.b && a.pts.size() < 3) continue;   // contracted to nothing
        a.length = polylineLength(a.pts);
        if (a.length <= 0.0) continue;
        used[a.a] = used[a.b] = 1;
        keptArcs.push_back(std::move(a));
    }

    std::vector<int> compact(newNodes.size(), -1);
    std::vector<QuadLayout::Node> finalNodes;
    for (size_t i = 0; i < newNodes.size(); ++i) {
        if (!used[i]) continue;
        compact[i] = static_cast<int>(finalNodes.size());
        finalNodes.push_back(newNodes[i]);
    }
    for (auto &a : keptArcs) { a.a = compact[a.a]; a.b = compact[a.b]; }

    layout_.rebuild(std::move(finalNodes), std::move(keptArcs));

    // A T-junction is a node where a separatrix stopped against a side. Merging
    // can leave one with nothing to stop against -- two of them nose to nose
    // become a single line through a node of valence two -- so the label is
    // brought back into line with what the node now is, or the count of them
    // would never fall.
    {
        std::vector<QuadLayout::Node> ns = layout_.getNodes();
        bool changed = false;
        for (auto &n : ns) {
            if (n.kind != QuadLayout::NodeKind::TJunction) continue;
            if (n.darts.size() == 3) continue;
            n.kind = (n.darts.size() <= 2) ? QuadLayout::NodeKind::Heteroclinic
                                           : QuadLayout::NodeKind::Crossing;
            changed = true;
        }
        if (changed) layout_.rebuild(std::move(ns), layout_.getArcs());
    }
    return true;
}

bool PartitionSimplify::hasSpur() const {
    std::vector<int> faceOfDart(2 * layout_.getArcs().size(), -1);
    const auto &fs = layout_.getFaces();
    for (size_t f = 0; f < fs.size(); ++f)
        for (const int d : fs[f].darts)
            if (d >= 0 && d < static_cast<int>(faceOfDart.size()))
                faceOfDart[d] = static_cast<int>(f);
    for (size_t a = 0; a < layout_.getArcs().size(); ++a)
        if (faceOfDart[2 * a] >= 0 && faceOfDart[2 * a] == faceOfDart[2 * a + 1]) return true;
    return false;
}

// ---------------------------------------------------------------------------
// mergeLenses()  --  a component bounded by two arcs is not a component
//
// Contracting a triangle's short side leaves its other two sides running
// between the same pair of nodes with nothing between them. Deleting one of
// them puts the two components that were either side of the pinched triangle
// against each other, which is what the triangle was standing in the way of.
// ---------------------------------------------------------------------------
int PartitionSimplify::mergeLenses() {
    int merged = 0;
    for (int guard = 0; guard < 64; ++guard) {
        const auto &faces = layout_.getFaces();
        int drop = -1;
        for (const auto &f : faces) {
            if (f.darts.size() != 2) continue;
            // Keep the longer of the two: it is the one the neighbouring
            // components were traced against, and the shorter is the detour.
            const int a0 = f.darts[0] >> 1, a1 = f.darts[1] >> 1;
            if (a0 == a1) continue;
            const auto &arcs = layout_.getArcs();
            if (fixedArc(arcs[a0]) && fixedArc(arcs[a1])) continue;
            if (fixedArc(arcs[a0])) { drop = a1; break; }
            if (fixedArc(arcs[a1])) { drop = a0; break; }
            drop = (arcs[a0].length < arcs[a1].length) ? a0 : a1;
            break;
        }
        if (drop < 0) break;

        std::vector<QuadLayout::Arc> keep;
        keep.reserve(layout_.getArcs().size());
        for (size_t i = 0; i < layout_.getArcs().size(); ++i)
            if (static_cast<int>(i) != drop) keep.push_back(layout_.getArcs()[i]);

        std::vector<QuadLayout::Node> ns = layout_.getNodes();
        for (auto &n : ns) { n.darts.clear(); n.angles.clear(); }
        layout_.rebuild(std::move(ns), std::move(keep));
        ++merged;
    }
    return merged;
}

// ---------------------------------------------------------------------------
// removeOneSliver()  --  contract the shortest side the layout has no use for
// ---------------------------------------------------------------------------
bool PartitionSimplify::removeOneSliver() {
    if (meshEdge_ <= 0.0) return false;
    const double limit = settings_.sliverSide * meshEdge_;

    // Candidates: an arc short enough that its two ends are the same place as
    // far as the mesh can tell, which is a whole side of some component that is
    // not four-sided. Restricting it to those keeps this to repairing the
    // degenerate components rather than quietly re-running Sec. 4 with a
    // different rule -- a short arc between two healthy quads is a real edge.
    const auto &arcs = layout_.getArcs();
    const auto &nodes = layout_.getNodes();
    std::vector<char> wholeSideOfBadFace(arcs.size(), 0);
    for (const auto &f : layout_.getFaces()) {
        if (f.corners == 4) continue;
        for (const auto &side : f.sides)
            if (side.size() == 1) wholeSideOfBadFace[side[0] >> 1] = 1;
    }

    int best = -1;
    double bestLen = limit;
    for (size_t i = 0; i < arcs.size(); ++i) {
        const auto &a = arcs[i];
        if (!wholeSideOfBadFace[i]) continue;
        if (a.a == a.b) continue;
        if (a.length >= bestLen) continue;
        const bool sA = isSingularity(a.a), sB = isSingularity(a.b);
        if (sA && sB) continue;                  // Sec. 4 condition 1
        if ((sA && isOnBoundary(a.b)) || (sB && isOnBoundary(a.a))) continue;  // condition 2
        // A node of the interface network goes nowhere, and takes nothing
        // that is itself fixed.
        const bool nA = isNetworkNode(a.a), nB = isNetworkNode(a.b);
        if ((nA && (nB || sB)) || (nB && sA)) continue;
        if ((nA || nB) && (nodes[a.a].kind == QuadLayout::NodeKind::BoundaryCorner ||
                           nodes[a.b].kind == QuadLayout::NodeKind::BoundaryCorner))
            continue;
        if (fixedArc(a)) {
            // A piece of the boundary may be contracted: the two nodes at its
            // ends become one and the model keeps its shape, because the arc
            // carrying on past the dying end swallows this one's geometry. What
            // may not go is a corner of the model itself, so a piece running
            // between two of them stays.
            if (nodes[a.a].kind == QuadLayout::NodeKind::BoundaryCorner &&
                nodes[a.b].kind == QuadLayout::NodeKind::BoundaryCorner) continue;
        } else if (isOnBoundary(a.a) && isOnBoundary(a.b)) {
            continue;   // a chord across the model: contracting it pinches it shut
        }
        best = static_cast<int>(i);
        bestLen = a.length;
    }
    if (best < 0) return false;

    const auto &a = arcs[best];
    // Whichever end may not move decides where the merged node goes.
    int keep = a.a, die = a.b;
    if (isSingularity(a.b) || (isOnBoundary(a.b) && !isOnBoundary(a.a))) { keep = a.b; die = a.a; }
    if (fixedArc(a) && nodes[a.b].kind == QuadLayout::NodeKind::BoundaryCorner) {
        keep = a.b; die = a.a;    // a corner of the model stays put
    } else if (fixedArc(a) && nodes[a.a].kind == QuadLayout::NodeKind::BoundaryCorner) {
        keep = a.a; die = a.b;
    }
    if (isNetworkNode(a.b)) { keep = a.b; die = a.a; }   // and so does a node of the network
    else if (isNetworkNode(a.a)) { keep = a.a; die = a.b; }
    const Point at = nodes[keep].pos;

    // Contracting a piece of the boundary must not shorten the model, so the
    // boundary arc on the far side of the dying node takes over its geometry:
    // the outline is the same curve afterwards, carried by one arc fewer.
    int absorbInto = -1;
    std::vector<Point> absorbed;
    if (fixedArc(a)) {
        for (const int d : nodes[die].darts) {
            const int arc = d >> 1;
            if (arc == best || !fixedArc(arcs[arc]) || arcs[arc].onInterface != a.onInterface) continue;
            absorbInto = arc;
            break;
        }
        if (absorbInto < 0) return false;
        absorbed = a.pts;                       // oriented from `die` towards `keep`
        if (a.a != die) std::reverse(absorbed.begin(), absorbed.end());
    }

    std::vector<QuadLayout::Node> ns = nodes;
    if (nodes[keep].kind != QuadLayout::NodeKind::InterfaceNode)
        ns[keep].kind = (static_cast<int>(nodes[keep].kind) < static_cast<int>(nodes[die].kind))
                            ? nodes[keep].kind : nodes[die].kind;
    for (auto &n : ns) { n.darts.clear(); n.angles.clear(); }

    std::vector<QuadLayout::Arc> keptArcs;
    keptArcs.reserve(arcs.size());
    for (size_t i = 0; i < arcs.size(); ++i) {
        if (static_cast<int>(i) == best) continue;
        QuadLayout::Arc c = arcs[i];
        if (static_cast<int>(i) == absorbInto) {
            // Splice the dying piece on at whichever end met it, so the arc now
            // runs all the way to the surviving node along the same curve.
            std::vector<Point> extra = absorbed;   // die -> keep
            if (c.a == die) {
                std::reverse(extra.begin(), extra.end());   // keep -> die
                extra.pop_back();
                extra.insert(extra.end(), c.pts.begin(), c.pts.end());
                c.pts = std::move(extra);
                c.a = keep;
            } else if (c.b == die) {
                c.pts.pop_back();
                c.pts.insert(c.pts.end(), extra.begin(), extra.end());
                c.b = keep;
            }
            c.length = polylineLength(c.pts);
            if (c.length > 0.0) keptArcs.push_back(std::move(c));
            continue;
        }
        if (c.a == die) { c.a = keep; c.pts.front() = at; }
        if (c.b == die) { c.b = keep; c.pts.back() = at; }
        if (c.a == c.b && c.pts.size() < 3) continue;
        c.length = polylineLength(c.pts);
        if (c.length <= 0.0) continue;
        keptArcs.push_back(std::move(c));
    }

    // Drop the node that died, and renumber.
    std::vector<char> used(ns.size(), 0);
    for (const auto &c : keptArcs) { used[c.a] = 1; used[c.b] = 1; }
    std::vector<int> compact(ns.size(), -1);
    std::vector<QuadLayout::Node> finalNodes;
    for (size_t i = 0; i < ns.size(); ++i) {
        if (!used[i]) continue;
        compact[i] = static_cast<int>(finalNodes.size());
        finalNodes.push_back(ns[i]);
    }
    for (auto &c : keptArcs) { c.a = compact[c.a]; c.b = compact[c.b]; }

    layout_.rebuild(std::move(finalNodes), std::move(keptArcs));
    mergeLenses();
    return true;
}

// ---------------------------------------------------------------------------
// enumerateChords()  --  every strip, once
//
// A component belongs to two chords, one through each pair of opposite sides,
// and the pair is what identifies which; sides are numbered round the
// component, so the two pairs are the even sides and the odd ones.
// ---------------------------------------------------------------------------
void PartitionSimplify::enumerateChords() {
    indexDarts();
    chords_.clear();
    const auto &faces = layout_.getFaces();
    std::vector<char> seen(faces.size() * 2, 0);

    for (size_t f = 0; f < faces.size(); ++f) {
        if (faces[f].sides.size() != 4) continue;
        for (int s = 0; s < 2; ++s) {
            if (seen[f * 2 + s]) continue;
            Chord c;
            if (!walkChord(static_cast<int>(f), s, c)) continue;
            for (size_t i = 0; i < c.faces.size(); ++i) seen[c.faces[i] * 2 + (c.entry[i] % 2)] = 1;
            analyse(c);
            chords_.push_back(std::move(c));
        }
    }
}

// ---------------------------------------------------------------------------
// run()  --  Algorithm 3
// ---------------------------------------------------------------------------
void PartitionSimplify::run() {
    // The degenerate components go first and are re-checked after every
    // collapse: each one removed is a chord that could not be walked before,
    // so doing them first is what lets Sec. 4 see the layout it was meant to
    // be handed.
    auto sweepSlivers = [&]() {
        while (settings_.removeSlivers) {
            const QuadLayout before = layout_;
            const auto &rb = before.getReport();
            if (!removeOneSliver()) break;
            const auto &ra = layout_.getReport();
            const bool ok = !hasSpur() && ra.arcCrossings == 0 && ra.danglingEnds == 0 && ra.faces > 0 &&
                            std::fabs(interfaceLength(layout_) - interfaceLength(before)) <=
                                1e-9 * std::max(1.0, interfaceLength(before)) &&
                            ra.faces < rb.faces && ra.singularities == rb.singularities &&
                            std::fabs(ra.totalArea - rb.totalArea) <=
                                1e-6 * std::max(1.0, rb.totalArea) &&
                            (rb.faces - ra.faces) >= (rb.quadFaces - ra.quadFaces);
            if (!ok) { layout_ = before; break; }
            ++report_.slivers;
        }
    };
    // A component bounded by two arcs is junk whatever produced it, so clear
    // any the tracing handed over before looking for chords -- not only the
    // ones a contraction creates.
    if (settings_.removeSlivers) {
        const QuadLayout before = layout_;
        const int n = mergeLenses();
        if (n > 0) {
            const auto &rb = before.getReport();
            const auto &ra = layout_.getReport();
            if (hasSpur() || ra.arcCrossings != 0 || ra.danglingEnds != 0 || ra.faces <= 0 ||
                std::fabs(interfaceLength(layout_) - interfaceLength(before)) >
                    1e-9 * std::max(1.0, interfaceLength(before)) ||
                ra.singularities != rb.singularities ||
                std::fabs(ra.totalArea - rb.totalArea) > 1e-6 * std::max(1.0, rb.totalArea))
                layout_ = before;
            else
                report_.slivers += n;
        }
    }
    sweepSlivers();

    auto collapseAll = [&]() {
    for (int iter = 0; iter < settings_.maxCollapses; ++iter) {
        enumerateChords();

        report_.chordsSeen = static_cast<int>(chords_.size());
        report_.blockedByConditions = 0;
        report_.blockedByEnergy = 0;
        report_.blockedByDrag = 0;
        for (int &b : report_.blockCount) b = 0;
        chordSummaries_.clear();
        std::vector<int> candidates;
        for (size_t i = 0; i < chords_.size(); ++i) {
            const Chord &c = chords_[i];
            Block why = c.block;
            if (c.collapsible && c.energy <= 0.0) why = Block::Energy;
            else if (c.collapsible && meshEdge_ > 0.0 &&
                     c.maxWidth > settings_.maxDrag * meshEdge_) why = Block::Drag;
            chordSummaries_.push_back(ChordSummary{c.faces, c.collapsible && c.energy > 0.0,
                                                   c.energy, c.minWidth, why});
            ++report_.blockCount[static_cast<int>(why)];
            if (!c.collapsible) { ++report_.blockedByConditions; continue; }
            if (c.energy <= 0.0) { ++report_.blockedByEnergy; continue; }
            if (meshEdge_ > 0.0 && c.maxWidth > settings_.maxDrag * meshEdge_) {
                ++report_.blockedByDrag; continue;
            }
            candidates.push_back(static_cast<int>(i));
        }
        if (candidates.empty()) break;

        // The thinnest first: collapsing a wide strip distorts more than it
        // simplifies, and taking the thin ones first tends to leave the wide
        // ones no longer collapsible at all.
        std::sort(candidates.begin(), candidates.end(), [&](int a, int b) {
            if (settings_.tJunctionsFirst && (chords_[a].tJunctionEnds > 0) !=
                                             (chords_[b].tJunctionEnds > 0))
                return chords_[a].tJunctionEnds > chords_[b].tJunctionEnds;
            return chords_[a].minWidth < chords_[b].minWidth;
        });

        const QuadLayout before = layout_;
        const auto &rb = before.getReport();
        int lensesBefore = 0;
        for (const auto &f : before.getFaces()) if (f.corners < 3) ++lensesBefore;
        bool done = false;
        for (const int ci : candidates) {
            if (!collapse(chords_[ci])) {
                layout_ = before;
                ++report_.rolledBack;
                ++report_.rbFailed;
                continue;
            }
            const auto &ra = layout_.getReport();
            // Proposition 2, as a test rather than an assumption.
            //
            // Plus one thing Proposition 2 takes for granted: that every arc
            // still has a different component on each side. An arc with the
            // same one on both is a spur hanging into it rather than a wall
            // between two, so that component is not a disc -- its walk goes out
            // along the spur and back, turning 180 degrees at the tip and
            // picking up two more corners each time, which is where the seven-
            // to ten-cornered components came from. It happens where a node the
            // collapse merges sits a fifth of an element from a corner of the
            // model: the arc between them survives with both its ends on the
            // same merged node.
            const bool spur = hasSpur();

            // A component with fewer than three corners is a lens between two
            // arcs that share both ends -- not something a four-sided region
            // can degenerate into by losing a strip, so a collapse that leaves
            // one has merged two sides that were not the two sides of a strip.
            int lenses = 0;
            for (const auto &f : layout_.getFaces()) if (f.corners < 3) ++lenses;

            bool ok = true;
            if (spur) { ++report_.rbSpur; ok = false; }
            else if (lenses > lensesBefore) { ++report_.rbLens; ok = false; }
            else if (ra.arcCrossings != 0) { ++report_.rbCrossings; ok = false; }
            else if (ra.danglingEnds != 0) { ++report_.rbDangling; ok = false; }
            else if (ra.faces <= 0 || ra.faces >= rb.faces) { ++report_.rbNotFewer; ok = false; }
            else if (ra.tJunctions > rb.tJunctions) { ++report_.rbMoreT; ok = false; }
            else if (ra.singularities != rb.singularities) { ++report_.rbSing; ok = false; }
            else if (std::fabs(ra.totalArea - rb.totalArea) >
                     1e-6 * std::max(1.0, rb.totalArea)) { ++report_.rbArea; ok = false; }
            else if (std::fabs(interfaceLength(layout_) - interfaceLength(before)) >
                     1e-9 * std::max(1.0, interfaceLength(before))) { ++report_.rbInterface; ok = false; }
            else if ((rb.faces - ra.faces) < (rb.quadFaces - ra.quadFaces)) {
                ++report_.rbWorse; ok = false;
            }
            if (ok) { ++report_.collapses; done = true; sweepSlivers(); break; }
            ++report_.rolledBack;
            layout_ = before;
        }
        if (!done) break;
    }
    };
    collapseAll();

    // Sec. 12: the T-junctions no collapse could remove, continued to the
    // boundary where the field allows, and the strips the continuations leave
    // collapsed in turn. The collapses cannot bring a T-junction back
    // (Proposition 2, checked above), so what this removes stays removed.
    report_.tJunctionsBeforeStems = layout_.getReport().tJunctions;
    if (settings_.extendStems && tracer_ && layout_.getReport().tJunctions > 0) {
        StemExtension ext(layout_, *tracer_, settings_.stems);
        ext.run();
        report_.stems = ext.getReport();
        if (report_.stems.extended > 0) {
            layout_ = ext.getLayout();
            sweepSlivers();
            collapseAll();
        }
    }

    const auto &r = layout_.getReport();
    report_.componentsAfter = r.faces;
    report_.tJunctionsAfter = r.tJunctions;
    report_.quadsAfter = r.quadFaces;
}
