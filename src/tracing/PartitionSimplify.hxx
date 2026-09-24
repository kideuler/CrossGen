#ifndef __PARTITION_SIMPLIFY_HXX__
#define __PARTITION_SIMPLIFY_HXX__

#include <vector>

#include "tracing/QuadLayout.hxx"
#include "tracing/StemExtension.hxx"

// ---------------------------------------------------------------------------
// Partition simplification: Sec. 4 of Viertel, Osting and Staten, IMR 2019,
// step 3 of their Algorithm 1.
//
// Separatrix tracing gives a valid quad layout that is far finer than it needs
// to be. Singularities on a discrete mesh never line up exactly, so two
// separatrices that ought to be one curve run side by side down the model and
// leave a strip of components between them a fraction of an edge wide. This
// removes those strips.
//
// The operation is one chord collapse, extended from quad meshes to layouts
// with T-junctions. A chord is a strip of components, entered and left through
// opposite sides, running until a T-junction or the boundary stops it. It has
// two longitudinal sides and a ladder of transverse rungs between them.
// Collapsing it contracts every rung to a point, which merges the two
// longitudinal sides into a single curve and removes the strip:
//
//   * where the strip is bounded by singularities on opposite corners -- a zip
//     patch -- the merged curve is a weighted blend of the two sides, running
//     from one singularity to the other, so both stay where the field put them;
//
//   * otherwise -- a non-zip patch -- the side carrying the singularities is
//     kept as it stands and the other is deleted, again so that nothing moves
//     that should not. The boundary is treated the same way: it is never the
//     side that gets deleted.
//
// A chord is divided into patches at every rung carrying a singularity, and
// each patch decides its own case; the collapse is of the whole chord.
//
// Which chords may be collapsed at all is Sec. 4's three conditions, and which
// are worth collapsing is Sec. 4.1's energy. Algorithm 3 is then greedy: take
// the thinnest collapsible chord, collapse it, and look again.
//
// Proposition 2 says a collapse leaves a T-layout with the same irregular nodes,
// strictly fewer components and no more T-junctions than before. Those are
// checkable, so they are checked: a collapse that violates any of them is
// undone and its chord struck off. That is a guard against this reading of
// Sec. 4's conditions being wrong somewhere, not a substitute for them -- the
// count of collapses it rejects is reported, and on the models here it is zero.
// ---------------------------------------------------------------------------

class PartitionSimplify {
public:
    // Why a chord could not be collapsed, for reporting: which of Sec. 4's
    // conditions -- or which of the extra ones the boundary imposes -- stopped
    // it. Knowing this is what tells a partition that is as coarse as the
    // operation allows from one that is merely stuck.
    enum class Block {
        None = 0,
        RungNotContractible,   // a node sits along a rung and would be squashed
        RungJoinsSingularities,// Sec. 4 condition 1
        RungSingularityToBoundary, // condition 2
        StripBetweenBoundaries,
        ZipAgainstBoundary,
        WouldDeleteBoundary,
        TJunctionHasNowhereToGo,   // condition 3
        NoPatches,
        Energy,                // Sec. 4.1
        Drag,                  // Settings::maxDrag
        RungAcrossInterface    // a rung from an interface to another curve of the model
    };
    static const char *blockName(Block b);

    struct Settings {
        // Sec. 4.1: a zip patch is worth collapsing when the angle its
        // diagonal makes with its base is under this, i.e. when the strip is
        // more than about 2.4 times as long as it is wide. Collapsing a fatter
        // one bends the separatrices too far at the singularities it joins.
        double zipAngle = M_PI / 8.0;

        // Sec. 4.1 leaves the energy of a non-zip patch at "a positive
        // constant": deleting a side with nothing on it costs the partition
        // nothing, so those are always worth taking.
        bool collapseNonZip = true;
        // ... unless this is positive, when a non-zip patch is judged by the
        // same aspect ratio a zip is, against this angle instead of zipAngle.
        double nonZipAngle = 0.0;

        // How far a collapse may drag anything, in mean mesh edges.
        //
        // Sec. 4.1's energy is a ratio -- the strip's width against the chord's
        // length -- so it is scale-free, which is right for the thing it is
        // stated to protect: the angles separatrices leave a singularity at. It
        // does not bound how far anything moves. On a chord thirty elements
        // long it admits a strip twelve wide, and contracting that drags every
        // separatrix that ended on the dying side the full width across, which
        // is how a curve ends up nowhere near the streamline it was traced
        // along and the components either side of it end up with a cusp for a
        // corner. On geom007 the greedy loop works through the real slivers in
        // its first four collapses and then spends nine more on strips up to
        // nine elements wide, dragging one node fifteen.
        //
        // The case for collapsing at all is that the two sides are one curve
        // the discretisation split in two, and that is an argument about edge
        // lengths, so the bound is in edge lengths. It is a cap on top of
        // Sec. 4.1, not a replacement: the energy still decides among the
        // strips narrow enough to be discretisation error in the first place.
        //
        // It was 1.5 while collapse() treated every patch as a non-zip keeping
        // its right-hand side (see there) and the model's corners were not
        // fixed nodes (fixedCorners): a wider strip then meant a longer drag in
        // an arbitrary direction. With both fixed it was measured again over
        // the 46 models, meshed and smoothed (TMOP, 50 sweeps):
        //
        //   cap     blocks  T-jcts  covered  mean worst SJ  models < 0.3
        //   1.5      2612     36      38/46       0.734            3
        //   3        2029     23      40/46       0.742            4
        //   4        1852     24      40/46       0.758            3
        //   5        1684     17      40/46       0.750            3
        //   6        1640     16      40/46       0.735            4
        //   none     ~1540    ~18     37/46       0.275           18
        //
        // Four has the best elements of any; five and six trade a tenth fewer
        // blocks for a drop on a handful of models (geom009 0.81 -> 0.64).
        // With no cap at all a non-zip strip -- whose Sec. 4.1 energy is a
        // constant -- of any width is taken, and geom017 and geom019 collapse
        // to a single block with a straight angle for a corner. Judging a
        // non-zip by the zip's aspect ratio instead (nonZipAngle) does not
        // separate those from the good ones, so the cap stays.
        double maxDrag = 4.0;

        // Keep the curve a zip leaves smooth. Two things put a hook in it
        // otherwise:
        //
        //   * the blend stepping per component rather than by arc length, so a
        //     component a fiftieth of the patch long next to a singularity
        //     takes a whole step of the sideways move;
        //   * a rung that is a piece of a material interface staying on one
        //     side of the strip -- it may not leave the interface -- while the
        //     blend runs down the middle, so the merged curve hooks out to it
        //     at every interface it crosses.
        //
        // With this on, the blend runs with arc length and is eased (zero
        // slope at both ends), so the merged curve leaves and meets each
        // singularity along the separatrix traced there; and an interface rung
        // is contracted to the point where the blend crosses it, which is on
        // the interface, the interface arcs either side taking its two halves.
        // Off is the per-component blend, for comparison.
        //
        // What it does not do is bend the arcs a collapse re-attaches. Moving
        // an arc's end drags its last segment across the strip, which is a
        // hook of its own, but easing that move along the arc makes the arc
        // arrive tangent to the curve it was merged onto: the component
        // between them gets a cusp for a corner, the corner rule loses it, and
        // measured over the corpus that cost 49 four-sided components and
        // seven fully covered models. Those hooks survive only where the
        // collapse that would merge the two curves next is itself refused.
        bool smoothCollapse = true;

        // Treat the model's corners and the interface network's nodes as
        // fixed in Sec. 4's patch logic, as docs/viertel_2019.md Sec. 6 does
        // ("fixed nodes = SING u BSING"); see isFixed(). Off counts interior
        // singularities alone, for comparison.
        bool fixedCorners = true;


        // Contract a rung even when a node sits along it rather than only at
        // its ends. That node -- a separatrix that ended on this side -- is
        // squashed onto the merged node, which is the paper's "the hanging
        // separatrix is simply extended after the collapse operation until it
        // crosses the next separatrix" applied to a rung instead of a
        // longitudinal side. Safe only because maxDrag bounds how far the
        // squash can move it; without that it is a licence to fold the layout.
        //
        // Off by default because it buys nothing measurable: over the sixteen
        // models it frees 23 chords and every one of them is then stopped by
        // maxDrag or undone by Proposition 2, leaving the same 825 components
        // and one more T-junction for 9 more rolled-back attempts.
        bool contractRungsWithNodes = false;

        // Sec. 6: "in each case that we observed, all T-junctions could have
        // been removed from the initial partition by collapsing the chords in a
        // different order, which suggests that perhaps a better collapse order
        // would prioritize or even require collapsing chords that end in
        // T-junctions". Algorithm 3 orders on width alone; this takes the
        // chords that would resolve a T-junction first and orders on width
        // within each group.
        bool tJunctionsFirst = true;

        // A component the tracing pinched to nothing: two separatrices that
        // cross twice within a fraction of an element leave a triangle a
        // seventh of an element on its short side and a quarter of a square
        // element in area. Sec. 4 cannot reach one -- a chord is a run of
        // four-sided components, so a three-sided one is not on any chord, and
        // worse, it stops every chord that would otherwise run through it.
        // Contracting its short side removes it and unblocks them.
        //
        // In mean mesh edges. Below the mesh resolution the two nodes are one
        // node as far as the discretisation can tell, which is the same
        // argument maxDrag rests on.
        double sliverSide = 0.5;
        bool removeSlivers = true;

        // docs/viertel_2019.md Sec. 12: once no chord is left to collapse,
        // trace the stem of every T-junction that remains on through the side
        // it stopped on to the boundary (StemExtension), and then collapse
        // again. Needs the field, so it runs only when the layout handed in
        // was built from a trace (QuadLayout::getTrace()).
        bool extendStems = true;
        StemExtension::Settings stems;

        int maxCollapses = 100000;
    };

    struct Report {
        int componentsBefore = 0;
        int componentsAfter = 0;
        int tJunctionsBefore = 0;
        int tJunctionsAfter = 0;
        int quadsBefore = 0;
        int quadsAfter = 0;
        int collapses = 0;
        int chordsSeen = 0;          // at the last pass
        int blockedByConditions = 0; // at the last pass: Sec. 4's three conditions
        int blockedByEnergy = 0;     // at the last pass: Sec. 4.1
        int blockedByDrag = 0;       // at the last pass: wider than Settings::maxDrag
        int slivers = 0;             // degenerate components contracted away
        int rolledBack = 0;          // collapses undone for breaking Proposition 2
        // How many chords each reason accounted for, at the last pass.
        int blockCount[12] = {0};
        // Which part of Proposition 2 the rolled-back ones broke.
        int rbCrossings = 0, rbDangling = 0, rbNotFewer = 0, rbMoreT = 0, rbSing = 0, rbArea = 0,
            rbWorse = 0, rbFailed = 0, rbSpur = 0, rbLens = 0, rbInterface = 0;
        // Settings::extendStems: what the continuation pass did, and the
        // T-junctions there were before it and after the collapses that
        // followed it.
        StemExtension::Report stems;
        int tJunctionsBeforeStems = 0;
    };

    PartitionSimplify(const QuadLayout &layout, const Settings &settings);
    explicit PartitionSimplify(const QuadLayout &layout)
        : PartitionSimplify(layout, Settings()) {}

    // Algorithm 3.
    void run();

    const QuadLayout &getLayout() const { return layout_; }
    const Report &getReport() const { return report_; }

    // The strips the last pass found, for drawing: each is the list of
    // components it runs through, with whether it was collapsible.
    struct ChordSummary {
        std::vector<int> faces;
        bool collapsible = false;
        double energy = 0.0;
        double minWidth = 0.0;
        Block block = Block::None;
    };
    const std::vector<ChordSummary> &getChords() const { return chordSummaries_; }

private:
    // One rung of the ladder: a whole side of a component, with the node at
    // each end labelled by which longitudinal side of the chord it belongs to.
    struct Rung {
        std::vector<int> darts;   // in the traversal order of the component that owns them
        int endL = -1, endR = -1;
        double length = 0.0;
        bool onBoundary = false;   // a piece of dS or of an interface (a curve of the model)
        bool onInterface = false;  // ... of an interface
        bool contractible = true; // no node other than a plain join sits along it
    };

    struct Chord {
        std::vector<int> faces;
        std::vector<int> entry;              // which side of each face it comes in by
        std::vector<Rung> rungs;             // faces + 1 of them, or faces if cyclic
        std::vector<std::vector<int>> sideL; // per component, its darts on each side
        std::vector<std::vector<int>> sideR;
        bool cyclic = false;
        double minWidth = 0.0;
        double maxWidth = 0.0;   // the widest rung: how far a collapse would drag things
        int tJunctionEnds = 0;   // T-junctions at the ends of the chord that it would resolve
        double energy = 0.0;
        bool collapsible = false;
        Block block = Block::None;
    };

    // A run of components of one chord between two rungs carrying
    // singularities, which is what decides zip against non-zip.
    struct Patch {
        int first = 0, last = 0;  // rung indices bounding it
        bool zip = false;
        bool keepL = false;       // non-zip: which side survives
        bool ok = false;
        double energy = 0.0;
    };

    void enumerateChords();
    bool walkChord(int seedFace, int seedSide, Chord &out) const;
    // Whether side `s` of face `f` is shared whole with one side of one other
    // face, which is what lets a chord carry on through it.
    bool sharedSide(int f, int s, int &g, int &t) const;

    void analyse(Chord &c) const;
    std::vector<Patch> patchesOf(const Chord &c) const;
    bool patchOk(const Chord &c, Patch &p, Block &why) const;
    double patchEnergy(const Chord &c, const Patch &p) const;

    bool collapse(const Chord &c);

    bool isSingularity(int node) const;
    // A singularity, or -- with Settings::fixedCorners -- a corner of the
    // model or a node of the interface network: Sec. 6's SING u BSING.
    bool isFixed(int node) const;
    // On dS or on an interface: a curve of the model, which nothing may leave.
    bool isOnBoundary(int node) const;
    bool isOnInterface(int node) const;
    bool isNetworkNode(int node) const;
    // The arc a T-junction's own separatrix ends along, as opposed to the two
    // that carry the side it ended on.
    int stemArcOf(int node) const;

    QuadLayout layout_;
    Settings settings_;
    Report report_;
    std::vector<Chord> chords_;
    std::vector<ChordSummary> chordSummaries_;
    std::vector<int> faceOfDart_;
    const FieldTracer *tracer_ = nullptr;   // from the layout's trace, for extendStems
    double tol_ = 1e-9;
    double meshEdge_ = 0.0;   // mean edge of the model, the scale maxDrag is in

    void indexDarts();

    // Contract the shortest degenerate side in the layout, and merge away any
    // two-sided component that leaves. Returns false when there is none left
    // that may be contracted.
    bool removeOneSliver();

    // Whether the layout has an arc with the same component on both sides: a
    // spur hanging into it rather than a wall between two, so that component is
    // not a disc. Neither operation may leave one.
    bool hasSpur() const;

    // Delete one arc of every two-arc component, which is what a contraction
    // leaves where it pinched a triangle shut. Returns how many it merged.
    int mergeLenses();

    friend class QuadLayout;
};

#endif // __PARTITION_SIMPLIFY_HXX__
