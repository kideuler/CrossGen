#ifndef __PARTITION_SIMPLIFY_HXX__
#define __PARTITION_SIMPLIFY_HXX__

#include <vector>

#include "tracing/QuadLayout.hxx"

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
        double maxDrag = 1.5;


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
        int rolledBack = 0;          // collapses undone for breaking Proposition 2
        // Which part of Proposition 2 the rolled-back ones broke.
        int rbCrossings = 0, rbDangling = 0, rbNotFewer = 0, rbMoreT = 0, rbSing = 0, rbArea = 0,
            rbWorse = 0, rbFailed = 0;
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
    };
    const std::vector<ChordSummary> &getChords() const { return chordSummaries_; }

private:
    // One rung of the ladder: a whole side of a component, with the node at
    // each end labelled by which longitudinal side of the chord it belongs to.
    struct Rung {
        std::vector<int> darts;   // in the traversal order of the component that owns them
        int endL = -1, endR = -1;
        double length = 0.0;
        bool onBoundary = false;
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
        double energy = 0.0;
        bool collapsible = false;
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
    bool patchOk(const Chord &c, Patch &p) const;
    double patchEnergy(const Chord &c, const Patch &p) const;

    bool collapse(const Chord &c);

    bool isSingularity(int node) const;
    bool isOnBoundary(int node) const;
    // The arc a T-junction's own separatrix ends along, as opposed to the two
    // that carry the side it ended on.
    int stemArcOf(int node) const;

    QuadLayout layout_;
    Settings settings_;
    Report report_;
    std::vector<Chord> chords_;
    std::vector<ChordSummary> chordSummaries_;
    std::vector<int> faceOfDart_;
    double tol_ = 1e-9;
    double meshEdge_ = 0.0;   // mean edge of the model, the scale maxDrag is in

    void indexDarts();

    friend class QuadLayout;
};

#endif // __PARTITION_SIMPLIFY_HXX__
