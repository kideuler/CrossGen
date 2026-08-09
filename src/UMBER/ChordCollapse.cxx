#include "UMBER/ChordCollapse.hxx"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace {

// A polyline resampled at n points evenly spaced along its length, which is
// what lets two curves of different resolutions be averaged point for point.
std::vector<Point> resample(const std::vector<Point> &p, int n) {
    std::vector<Point> out;
    if (p.size() < 2 || n < 2) return p;

    std::vector<double> s(p.size(), 0.0);
    for (size_t i = 1; i < p.size(); ++i) s[i] = s[i - 1] + normP(p[i] - p[i - 1]);
    const double total = s.back();
    if (total <= 0.0) return p;

    out.reserve(n);
    size_t j = 1;
    for (int k = 0; k < n; ++k) {
        const double target = total * static_cast<double>(k) / static_cast<double>(n - 1);
        while (j + 1 < p.size() && s[j] < target) ++j;
        const double d = s[j] - s[j - 1];
        const double t = (d > 0.0) ? (target - s[j - 1]) / d : 0.0;
        out.push_back(p[j - 1] * (1.0 - std::min(1.0, std::max(0.0, t))) +
                      p[j] * std::min(1.0, std::max(0.0, t)));
    }
    return out;
}

double lengthOf(const std::vector<Point> &p) {
    double l = 0.0;
    for (size_t i = 1; i < p.size(); ++i) l += normP(p[i] - p[i - 1]);
    return l;
}

struct DisjointSet {
    std::vector<int> p;
    explicit DisjointSet(size_t n) : p(n) { std::iota(p.begin(), p.end(), 0); }
    int find(int x) { while (p[x] != x) x = p[x] = p[p[x]]; return x; }
    void unite(int a, int b) {
        a = find(a); b = find(b);
        if (a != b) p[b] = a;
    }
};

} // namespace

const char *ChordCollapse::blockName(Block b) {
    switch (b) {
        case Block::None:               return "collapsible";
        case Block::NotFourSided:       return "a block on it is not four-sided";
        case Block::SelfIntersecting:   return "runs through a block twice";
        case Block::TooThick:           return "too thick";
        case Block::JoinsTwoBoundaries: return "would join two boundaries";
        case Block::PinchesBoundary:    return "would pinch the boundary";
        case Block::WouldEmptyDomain:   return "is the whole structure";
        case Block::RolledBack:         return "left a broken layout and was undone";
        default:                        return "?";
    }
}

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
ChordCollapse::ChordCollapse(const QuadLayout &layout, const Settings &settings)
    : layout_(layout), settings_(settings) {
    double total = 0.0;
    for (const auto &a : layout_.getArcs()) total += a.length;
    if (!layout_.getArcs().empty() && total > 0.0)
        widthScale_ = total / static_cast<double>(layout_.getArcs().size());

    report_.blocksBefore = static_cast<int>(layout_.getFaces().size());
    report_.nodesBefore = static_cast<int>(layout_.getNodes().size());
    report_.blocksAfter = report_.blocksBefore;
    report_.nodesAfter = report_.nodesBefore;
}

// ---------------------------------------------------------------------------
// The three lookups the walk is made of
// ---------------------------------------------------------------------------
void ChordCollapse::indexDarts() {
    faceOfDart_.assign(2 * layout_.getArcs().size(), -1);
    const auto &faces = layout_.getFaces();
    for (int f = 0; f < static_cast<int>(faces.size()); ++f)
        for (const int d : faces[f].darts)
            if (d >= 0 && d < static_cast<int>(faceOfDart_.size())) faceOfDart_[d] = f;
}

// Which boundary loop each node is on. Two nodes are on the same one when a
// chain of boundary arcs joins them, which is what rule 3 asks about: the outer
// boundary and the rim of a hole are different loops however close they run.
void ChordCollapse::classifyBoundary() {
    const auto &nodes = layout_.getNodes();
    const auto &arcs = layout_.getArcs();

    DisjointSet ds(nodes.size());
    onBoundaryNode_.assign(nodes.size(), 0);
    for (const auto &a : arcs) {
        if (!a.onBoundary || a.a < 0 || a.b < 0) continue;
        onBoundaryNode_[a.a] = 1;
        onBoundaryNode_[a.b] = 1;
        ds.unite(a.a, a.b);
    }

    loopOfNode_.assign(nodes.size(), -1);
    std::vector<int> idOfRoot(nodes.size(), -1);
    int next = 0;
    for (int n = 0; n < static_cast<int>(nodes.size()); ++n) {
        if (!onBoundaryNode_[n]) continue;
        const int r = ds.find(n);
        if (idOfRoot[r] < 0) idOfRoot[r] = next++;
        loopOfNode_[n] = idOfRoot[r];
    }
}

bool ChordCollapse::usable(int face) const {
    if (face < 0 || face >= static_cast<int>(layout_.getFaces().size())) return false;
    const auto &f = layout_.getFaces()[face];
    return f.darts.size() == 4 && f.corners == 4;
}

int ChordCollapse::sideOf(int face, int arc) const {
    const auto &darts = layout_.getFaces()[face].darts;
    for (size_t i = 0; i < darts.size(); ++i)
        if ((darts[i] >> 1) == arc) return static_cast<int>(i);
    return -1;
}

int ChordCollapse::across(int arc, int f) const {
    const int d0 = faceOfDart_[2 * arc], d1 = faceOfDart_[2 * arc + 1];
    if (d0 == d1) return -1;   // the same block on both sides: not a wall
    return (d0 == f) ? d1 : d0;
}

// ---------------------------------------------------------------------------
// walkChord()  --  the strip through one side
//
// Out of the block on one side of the seed and on through the side opposite,
// then the same the other way, so the chord comes out in order from one end to
// the other. A chord that closes on itself has no ends and is walked once.
// ---------------------------------------------------------------------------
bool ChordCollapse::walkChord(int seed, Chord &c) const {
    const auto &faces = layout_.getFaces();
    const auto &arcs = layout_.getArcs();
    if (seed < 0 || seed >= static_cast<int>(arcs.size())) return false;

    std::vector<char> onFace(faces.size(), 0);
    std::vector<char> onArc(arcs.size(), 0);
    onArc[seed] = 1;
    bool bad = false;

    auto walk = [&](int rung, int f, std::vector<int> &rr, std::vector<int> &ff) {
        while (f >= 0 && !bad && !c.cyclic) {
            if (!usable(f)) { c.block = Block::NotFourSided; bad = true; return; }
            if (onFace[f]) { c.block = Block::SelfIntersecting; bad = true; return; }
            onFace[f] = 1;
            ff.push_back(f);

            const int s = sideOf(f, rung);
            if (s < 0) { c.block = Block::NotFourSided; bad = true; return; }
            const int opp = faces[f].darts[(s + 2) % 4] >> 1;
            if (opp == seed) { c.cyclic = true; return; }
            if (onArc[opp]) { c.block = Block::SelfIntersecting; bad = true; return; }
            onArc[opp] = 1;
            rr.push_back(opp);

            f = across(opp, f);
            rung = opp;
        }
    };

    std::vector<int> rungsF, facesF, rungsB, facesB;
    walk(seed, faceOfDart_[2 * seed], rungsF, facesF);
    if (!bad && !c.cyclic) walk(seed, faceOfDart_[2 * seed + 1], rungsB, facesB);

    c.rungs.clear();
    c.faces.clear();
    c.rungs.insert(c.rungs.end(), rungsB.rbegin(), rungsB.rend());
    c.rungs.push_back(seed);
    c.rungs.insert(c.rungs.end(), rungsF.begin(), rungsF.end());
    c.faces.insert(c.faces.end(), facesB.rbegin(), facesB.rend());
    c.faces.insert(c.faces.end(), facesF.begin(), facesF.end());
    if (bad || c.faces.empty()) return false;

    // Label the two ends of every rung by which long side of the strip they
    // belong to, carrying the labelling from one rung to the next through the
    // long sides that join them. This is what makes "the side that survives"
    // mean the same thing along the whole chord.
    c.endL.assign(c.rungs.size(), -1);
    c.endR.assign(c.rungs.size(), -1);
    c.sideL.assign(c.faces.size(), -1);
    c.sideR.assign(c.faces.size(), -1);

    for (size_t i = 0; i < c.faces.size(); ++i) {
        const auto &f = faces[c.faces[i]];
        const int s = sideOf(c.faces[i], c.rungs[i]);
        if (s < 0) return false;
        const int a = f.nodes[s], b = f.nodes[(s + 1) % 4];
        const int u = f.nodes[(s + 2) % 4], v = f.nodes[(s + 3) % 4];

        if (i == 0) { c.endL[0] = a; c.endR[0] = b; }
        bool flipped;
        if (c.endL[i] == a) flipped = false;
        else if (c.endL[i] == b) flipped = true;
        else return false;   // the labelling did not carry: not a strip

        const size_t next = (i + 1) % c.rungs.size();
        const int nl = flipped ? u : v;
        const int nr = flipped ? v : u;
        if (c.endL[next] < 0) { c.endL[next] = nl; c.endR[next] = nr; }
        else if (c.endL[next] != nl || c.endR[next] != nr) return false;

        c.sideL[i] = (flipped ? f.darts[(s + 1) % 4] : f.darts[(s + 3) % 4]) >> 1;
        c.sideR[i] = (flipped ? f.darts[(s + 3) % 4] : f.darts[(s + 1) % 4]) >> 1;
    }

    c.width = 0.0;
    for (const int r : c.rungs) c.width = std::max(c.width, arcs[r].length);
    c.length = 0.0;
    for (size_t i = 0; i < c.faces.size(); ++i)
        c.length += 0.5 * (arcs[c.sideL[i]].length + arcs[c.sideR[i]].length);
    return true;
}

// ---------------------------------------------------------------------------
// applyRules()  --  the three conditions, over every rung before any of it moves
// ---------------------------------------------------------------------------
void ChordCollapse::applyRules(Chord &c) const {
    c.collapsible = false;
    if (c.block != Block::None || c.faces.empty()) return;

    if (c.faces.size() >= layout_.getFaces().size()) { c.block = Block::WouldEmptyDomain; return; }

    // Rule 3, and the same argument applied to one loop: a rung whose two ends
    // are both on the boundary may only be contracted when it lies on the
    // boundary itself, which is how a chord ends when it runs out of the model.
    const auto &arcs = layout_.getArcs();
    for (size_t i = 0; i < c.rungs.size(); ++i) {
        const int u = c.endL[i], v = c.endR[i];
        if (u < 0 || v < 0 || u == v) { c.block = Block::NotFourSided; return; }
        const int lu = loopOfNode_[u], lv = loopOfNode_[v];
        if (lu < 0 || lv < 0) continue;
        if (lu != lv) { c.block = Block::JoinsTwoBoundaries; return; }
        if (!arcs[c.rungs[i]].onBoundary) { c.block = Block::PinchesBoundary; return; }
    }

    // The same objection from the other direction. A block with boundary down
    // both of its long sides is the model itself where it has thinned to a
    // strip, and merging those two sides is merging two pieces of boundary --
    // the model would be severed there, or a hole closed up, exactly as rule 3
    // forbids at a rung. Where a chord runs into such a block it is refused
    // whole, the same as anywhere else.
    for (size_t i = 0; i < c.faces.size(); ++i) {
        if (c.sideL[i] < 0 || c.sideR[i] < 0 || c.sideL[i] == c.sideR[i]) {
            c.block = Block::NotFourSided;
            return;
        }
        if (!arcs[c.sideL[i]].onBoundary || !arcs[c.sideR[i]].onBoundary) continue;
        const int lu = loopOfNode_[arcs[c.sideL[i]].a], lv = loopOfNode_[arcs[c.sideR[i]].a];
        c.block = (lu != lv) ? Block::JoinsTwoBoundaries : Block::PinchesBoundary;
        return;
    }

    // Rule 1.
    if (c.width > settings_.maxWidth * widthScale_) { c.block = Block::TooThick; return; }
    if (settings_.minAspect > 0.0 && c.length < settings_.minAspect * c.width) {
        c.block = Block::TooThick;
        return;
    }

    c.collapsible = true;
}

// ---------------------------------------------------------------------------
// priorityOf()  --  rule 2, as an order on the nodes
// ---------------------------------------------------------------------------
int ChordCollapse::priorityOf(int node) const {
    if (node < 0 || node >= static_cast<int>(onBoundaryNode_.size())) return 0;
    if (!onBoundaryNode_[node]) return 0;
    return (layout_.getNodes()[node].kind == QuadLayout::NodeKind::BoundaryCorner) ? 2 : 1;
}

// ---------------------------------------------------------------------------
// collapse()  --  contract every rung, merge the two long sides, rebuild
//
// Both halves are the same operation seen twice: the rungs say which nodes are
// one node afterwards and the blocks say which pairs of arcs are one arc, so
// each is a disjoint set and the new layout is what the two of them leave. The
// blocks themselves are not deleted, they simply do not come back: a block of
// the chord has both its rungs contracted to a point and is bounded by one
// curve traversed twice, which the face walk does not produce.
// ---------------------------------------------------------------------------
bool ChordCollapse::collapse(const Chord &c) {
    const QuadLayout saved = layout_;
    const auto &nodesIn = layout_.getNodes();
    const auto &arcsIn = layout_.getArcs();

    const double areaBefore = layout_.getReport().totalArea;
    const int badBefore = layout_.getReport().badFaces;
    const int crossBefore = layout_.getReport().arcCrossings;
    double stripArea = 0.0;
    for (const int f : c.faces) stripArea += layout_.getFaces()[f].area;

    DisjointSet nds(nodesIn.size());
    for (size_t i = 0; i < c.rungs.size(); ++i) nds.unite(c.endL[i], c.endR[i]);

    // Where each merged node goes: the highest-ranking of the nodes that made
    // it, and the middle of them when they rank equally. Rule 2.
    std::vector<std::vector<int>> members(nodesIn.size());
    for (int n = 0; n < static_cast<int>(nodesIn.size()); ++n) members[nds.find(n)].push_back(n);

    std::vector<Point> posOfRoot(nodesIn.size(), Point{0.0, 0.0});
    std::vector<QuadLayout::NodeKind> kindOfRoot(nodesIn.size(), QuadLayout::NodeKind::Crossing);
    std::vector<int> chosen(nodesIn.size(), -1);   // -1: the members met in the middle
    for (int r = 0; r < static_cast<int>(nodesIn.size()); ++r) {
        if (members[r].empty()) continue;
        int best = -1;
        for (const int n : members[r]) best = std::max(best, priorityOf(n));

        Point sum{0.0, 0.0};
        int count = 0, winner = -1;
        QuadLayout::NodeKind kind = QuadLayout::NodeKind::Dangling;
        for (const int n : members[r]) {
            if (priorityOf(n) != best) continue;
            sum = sum + nodesIn[n].pos;
            ++count;
            winner = n;
            if (static_cast<int>(nodesIn[n].kind) < static_cast<int>(kind)) kind = nodesIn[n].kind;
        }
        if (count <= 0) return false;
        posOfRoot[r] = sum / static_cast<double>(count);
        kindOfRoot[r] = kind;
        chosen[r] = (count == 1) ? winner : -1;
    }

    // The two long sides of each block become one arc.
    DisjointSet ads(arcsIn.size());
    for (size_t i = 0; i < c.faces.size(); ++i) ads.unite(c.sideL[i], c.sideR[i]);

    std::vector<std::vector<int>> group(arcsIn.size());
    for (int a = 0; a < static_cast<int>(arcsIn.size()); ++a) {
        if (arcsIn[a].a < 0 || arcsIn[a].b < 0) continue;
        // An arc whose two ends became one node is a rung, and a rung is what
        // the collapse contracts away.
        if (nds.find(arcsIn[a].a) == nds.find(arcsIn[a].b)) continue;
        group[ads.find(a)].push_back(a);
    }

    // How much of an arc's geometry the node rule already endorsed: an arc both
    // of whose ends won their rung is the side the merged nodes are on, so it
    // is the side the merged curve should follow.
    auto scoreOf = [&](int a) {
        double s = 0.0;
        for (const int n : {arcsIn[a].a, arcsIn[a].b}) {
            const int r = nds.find(n);
            s += (chosen[r] < 0) ? 0.5 : (chosen[r] == n ? 1.0 : 0.0);
        }
        return s;
    };

    std::vector<QuadLayout::Node> newNodes;
    std::vector<QuadLayout::Arc> newArcs;
    std::vector<int> newIndex(nodesIn.size(), -1);
    auto nodeIndex = [&](int root) {
        if (newIndex[root] < 0) {
            QuadLayout::Node n;
            n.pos = posOfRoot[root];
            n.kind = kindOfRoot[root];
            n.source = -1;
            newIndex[root] = static_cast<int>(newNodes.size());
            newNodes.push_back(std::move(n));
        }
        return newIndex[root];
    };

    for (int g = 0; g < static_cast<int>(arcsIn.size()); ++g) {
        if (group[g].empty()) continue;

        // The representative: a piece of boundary if there is one, because the
        // boundary is the model and never moves; otherwise the side the node
        // rule already chose.
        int rep = group[g].front();
        double bestScore = -1.0;
        bool anyBoundary = false;
        for (const int a : group[g]) anyBoundary = anyBoundary || arcsIn[a].onBoundary;
        for (const int a : group[g]) {
            if (anyBoundary && !arcsIn[a].onBoundary) continue;
            const double s = scoreOf(a);
            if (s > bestScore) { bestScore = s; rep = a; }
        }

        const int A = nds.find(arcsIn[rep].a), B = nds.find(arcsIn[rep].b);
        std::vector<Point> pts = arcsIn[rep].pts;

        // Two sides that rank equally are one curve the partition split in two,
        // so the merged one runs between them rather than along either.
        if (!anyBoundary && group[g].size() == 2) {
            const int other = (group[g][0] == rep) ? group[g][1] : group[g][0];
            if (std::fabs(scoreOf(other) - bestScore) < 1e-12) {
                std::vector<Point> q = arcsIn[other].pts;
                if (nds.find(arcsIn[other].a) == B && nds.find(arcsIn[other].b) == A)
                    std::reverse(q.begin(), q.end());
                const int n = static_cast<int>(std::max(pts.size(), q.size()));
                const std::vector<Point> ra = resample(pts, n);
                const std::vector<Point> rb = resample(q, n);
                if (ra.size() == rb.size()) {
                    pts.assign(ra.size(), Point{0.0, 0.0});
                    for (size_t i = 0; i < ra.size(); ++i) pts[i] = (ra[i] + rb[i]) * 0.5;
                }
            }
        }

        pts.front() = posOfRoot[A];
        pts.back() = posOfRoot[B];
        std::vector<Point> clean;
        for (const Point &p : pts)
            if (clean.empty() || normP(p - clean.back()) > 1e-14) clean.push_back(p);
        if (clean.size() < 2) continue;
        if (A == B && clean.size() < 3) continue;

        QuadLayout::Arc arc;
        arc.a = nodeIndex(A);
        arc.b = nodeIndex(B);
        arc.separatrix = arcsIn[rep].separatrix;
        arc.onBoundary = anyBoundary;
        arc.length = lengthOf(clean);
        arc.pts = std::move(clean);
        newArcs.push_back(std::move(arc));
    }

    layout_.rebuild(std::move(newNodes), std::move(newArcs));

    // Proposition-style guard: the collapse is supposed to remove the chord's
    // blocks and change nothing else, and all of that is checkable. A collapse
    // that broke something is undone rather than reported, because the rules
    // above are a reading of what is safe and a wrong reading should cost a
    // chord rather than the structure.
    const auto &r = layout_.getReport();
    const bool okBlocks =
        r.faces > 0 &&
        r.faces == static_cast<int>(saved.getFaces().size()) - static_cast<int>(c.faces.size());
    const bool okBad = r.badFaces <= badBefore;
    const bool okCross = r.arcCrossings <= crossBefore;
    // The strip is what changes hands and nothing else should: whichever side
    // the merged curve ends up on, the blocks beside the chord take over as
    // much of it as the structure gives up, so the total moves by at most the
    // strip's own area either way. A collapse that moved more than that has
    // taken a bite out of something it was not supposed to touch.
    const bool okArea = std::fabs(r.totalArea - areaBefore) <=
                        2.0 * stripArea + 1e-9 * std::fabs(areaBefore);
    if (!okBlocks || !okBad || !okCross || !okArea) {
        if (!okBlocks) ++report_.rbBlocks;
        else if (!okBad) ++report_.rbBad;
        else if (!okCross) ++report_.rbCrossings;
        else ++report_.rbArea;
        layout_ = saved;
        return false;
    }
    return true;
}

// ---------------------------------------------------------------------------
// enumerateChords()  --  every side is a rung of exactly one chord
// ---------------------------------------------------------------------------
void ChordCollapse::enumerateChords() {
    indexDarts();
    classifyBoundary();
    chords_.clear();
    summaries_.clear();
    for (int i = 0; i < static_cast<int>(Block::Count); ++i) report_.blockCount[i] = 0;
    report_.chordsSeen = 0;
    report_.collapsible = 0;
    report_.badBlocks = 0;
    for (int f = 0; f < static_cast<int>(layout_.getFaces().size()); ++f)
        if (!usable(f)) ++report_.badBlocks;

    std::vector<char> seen(layout_.getArcs().size(), 0);
    for (int a = 0; a < static_cast<int>(layout_.getArcs().size()); ++a) {
        if (seen[a]) continue;
        Chord c;
        const bool walked = walkChord(a, c);
        seen[a] = 1;
        for (const int r : c.rungs)
            if (r >= 0 && r < static_cast<int>(seen.size())) seen[r] = 1;
        if (walked) applyRules(c);

        ChordSummary s;
        s.blocks = c.faces;
        s.rungs = c.rungs;
        s.cyclic = c.cyclic;
        s.width = c.width / widthScale_;
        s.length = c.length / widthScale_;
        s.collapsible = c.collapsible;
        s.block = c.block;

        ++report_.chordsSeen;
        if (c.collapsible) ++report_.collapsible;
        ++report_.blockCount[static_cast<int>(c.block)];

        chords_.push_back(std::move(c));
        summaries_.push_back(std::move(s));
    }
}

// ---------------------------------------------------------------------------
// run()  --  the thinnest first, and look again after each one
// ---------------------------------------------------------------------------
void ChordCollapse::run() {
    for (int pass = 0; pass < settings_.maxCollapses; ++pass) {
        enumerateChords();

        std::vector<int> order;
        for (int i = 0; i < static_cast<int>(chords_.size()); ++i)
            if (chords_[i].collapsible) order.push_back(i);
        std::sort(order.begin(), order.end(),
                  [&](int x, int y) { return chords_[x].width < chords_[y].width; });

        bool moved = false;
        for (const int i : order) {
            const double w = chords_[i].width / widthScale_;
            if (collapse(chords_[i])) {
                ++report_.collapses;
                report_.widestCollapsed = std::max(report_.widestCollapsed, w);
                moved = true;
                break;
            }
            ++report_.rolledBack;
        }
        if (!moved) break;
    }

    enumerateChords();
    report_.blocksAfter = static_cast<int>(layout_.getFaces().size());
    report_.nodesAfter = static_cast<int>(layout_.getNodes().size());
}

// ---------------------------------------------------------------------------
// collapseChordThrough()  --  one operation, on a side picked by hand
// ---------------------------------------------------------------------------
bool ChordCollapse::collapseChordThrough(int arc) {
    indexDarts();
    classifyBoundary();

    Chord c;
    if (!walkChord(arc, c)) {
        lastBlock_ = (c.block == Block::None) ? Block::NotFourSided : c.block;
        return false;
    }
    applyRules(c);
    if (!c.collapsible) { lastBlock_ = c.block; return false; }

    if (!collapse(c)) {
        ++report_.rolledBack;
        lastBlock_ = Block::RolledBack;
        return false;
    }
    ++report_.collapses;
    report_.widestCollapsed = std::max(report_.widestCollapsed, c.width / widthScale_);
    lastBlock_ = Block::None;

    enumerateChords();
    report_.blocksAfter = static_cast<int>(layout_.getFaces().size());
    report_.nodesAfter = static_cast<int>(layout_.getNodes().size());
    return true;
}
