#ifndef __CHORD_COLLAPSE_HXX__
#define __CHORD_COLLAPSE_HXX__

#include <vector>

#include "tracing/QuadLayout.hxx"

// ---------------------------------------------------------------------------
// Chord collapse on the meta-block structure of a polysquare.
//
// The motorcycle graph tessellates the model into four-sided blocks, and it
// makes more of them than the shape needs. Two reflex corners a little way
// apart send their iso-lines down the model side by side and leave a strip of
// blocks between them a fraction of a block wide; the strip is not a feature of
// the shape, it is the distance between two corners, and every block in it is
// paid for in the mesh that follows.
//
// A chord is that strip. In a quad layout it is a dual loop: enter a block
// through one side, leave through the side opposite, and carry on until the
// boundary stops you or you come back to where you started. The sides crossed
// on the way -- the rungs of the ladder -- are what a collapse contracts. Each
// rung becomes a single node, so each block on the chord is squashed flat and
// disappears, and the two long sides of the strip become one curve. The chord's
// blocks are removed and nothing else about the structure changes: the layout
// is still four-sided everywhere, with as many fewer blocks as the chord was
// long.
//
// Which one to collapse is chosen by an edge -- a side of a block -- and the
// rest follows from it, because the chord through a side is determined: at each
// block there is exactly one side opposite the one you came in by. So the walk
// finds every side that has to go along with the one that was picked, and the
// question is only whether the collapse is allowed.
//
// Three rules say when it is not.
//
//  1. Only thin chords. Contracting a rung drags everything on one side of it
//     across to the other, so the rung's length is how far the structure moves
//     and how much of the model the vanishing blocks give up. A strip a
//     fraction of a block wide is discretization -- two corners that ought to
//     have been one -- and a fat one is a real feature of the shape. The chord
//     is judged by its *widest* rung, not its average: one fat rung is enough
//     to make the collapse move the structure somewhere it does not belong,
//     and the fat one is exactly where the strip stops being a sliver. Where
//     the line falls is a heuristic, so it is Settings::maxWidth, and it is
//     stated relative to the mean side length of the structure rather than in
//     model units so that it means the same thing on any model.
//
//  2. A rung with one end on the boundary and one in the interior contracts to
//     the end on the boundary, always. The boundary is the model and the
//     interior node is a consequence of where the iso-lines happened to fall,
//     so moving the interior node costs nothing and moving the boundary changes
//     the shape. The same argument ranks the boundary nodes among themselves:
//     a corner of the polysquare is a corner of the model, while a node where
//     an iso-line ran into the boundary is only a node, so a rung joining the
//     two contracts to the corner. Two nodes of equal standing meet in the
//     middle. The arcs follow the nodes: where one of the two long sides of the
//     strip lies on the boundary it survives unchanged and the other is
//     deleted, and where neither does the two are blended into one.
//
//  3. A rung with both ends on the boundary cannot be contracted unless it lies
//     on the boundary itself. If the two ends are on different loops -- the
//     outer boundary and the rim of a hole -- contracting it joins them, and
//     the model no longer has a hole in it. If they are on the same loop, it
//     pinches the model in two at that point. Either way the collapse would
//     change what the domain is rather than how it is partitioned, so the whole
//     chord is refused; a rung that is itself a piece of boundary is the one
//     case where both ends being on the boundary is ordinary -- that is simply
//     how a chord ends when it runs into the edge of the model.
//
// A collapse either happens whole or not at all: the rules are checked over
// every rung of the chord before any of it moves.
// ---------------------------------------------------------------------------
class ChordCollapse {
public:
    // Why a chord was not collapsed. Knowing this is what tells a structure
    // that is as coarse as the operation allows from one that is merely stuck.
    enum class Block {
        None = 0,
        NotFourSided,       // a block on the chord is not a quad: nothing to walk
        SelfIntersecting,   // the chord runs through one block twice
        TooThick,           // rule 1
        JoinsTwoBoundaries, // rule 3, two loops
        PinchesBoundary,    // rule 3, one loop
        WouldEmptyDomain,   // the chord is every block there is
        RolledBack,         // collapsed, and undone for leaving a broken layout
        Count
    };
    static const char *blockName(Block b);

    struct Settings {
        // Rule 1, in units of the mean side length of the block structure as it
        // was handed over. A chord is collapsible only if every one of its
        // rungs is shorter than this.
        //
        // The useful range is well under 1: at 1 a chord as wide as a typical
        // block is a chord of square blocks, which is a partition of the model
        // and not a sliver. The default takes strips up to a third of a block
        // wide, which on the models here is comfortably above the slivers left
        // by two corners a few elements apart and below anything that reads as
        // a feature.
        double maxWidth = 0.35;

        // Rule 1 again, from the other side: a chord must be at least this many
        // times longer than it is wide. A short fat chord and a long thin one
        // can have the same rungs, and only the second is a sliver. 0 turns it
        // off and leaves maxWidth to decide alone.
        double minAspect = 0.0;

        int maxCollapses = 1000;
    };

    struct Report {
        int blocksBefore = 0;
        int blocksAfter = 0;
        int nodesBefore = 0;
        int nodesAfter = 0;
        int collapses = 0;
        int chordsSeen = 0;      // at the last pass
        int collapsible = 0;     // at the last pass
        int rolledBack = 0;      // collapses undone for leaving a broken layout
        // Which part of the guard the undone ones broke: the block count did
        // not fall by the length of the chord, a block came out not four-sided,
        // two sides crossed, or the structure gave up more of the model than
        // the strip it removed.
        int rbBlocks = 0, rbBad = 0, rbCrossings = 0, rbArea = 0;
        int blockCount[static_cast<int>(Block::Count)] = {0};  // at the last pass
        double widestCollapsed = 0.0;  // in the units of Settings::maxWidth
        // The layout the operation refuses to touch at all: a block that is not
        // four-sided means the motorcycle graph did not close somewhere, and
        // every chord through it is blocked.
        int badBlocks = 0;
    };

    // The layout is copied; the original is left as it was.
    ChordCollapse(const QuadLayout &layout, const Settings &settings);
    explicit ChordCollapse(const QuadLayout &layout)
        : ChordCollapse(layout, Settings()) {}

    // Collapse the thinnest collapsible chord, then look again, until none is
    // left that the rules allow.
    void run();

    // One operation, on the chord through one side of one block -- the entry
    // point for a side picked by hand rather than by the greedy loop. Returns
    // false if the rules refused it, with why in whyBlocked().
    bool collapseChordThrough(int arc);
    Block whyBlocked() const { return lastBlock_; }

    const QuadLayout &getLayout() const { return layout_; }
    const Report &getReport() const { return report_; }

    // The chords of the structure as it stands, for drawing and for choosing
    // one: every side of every block is a rung of exactly one of them.
    struct ChordSummary {
        std::vector<int> blocks;  // faces of the layout, in order along the chord
        std::vector<int> rungs;   // arcs, one more than blocks unless cyclic
        bool cyclic = false;
        double width = 0.0;       // the widest rung, in units of the mean side
        double length = 0.0;      // along the chord, in the same units
        bool collapsible = false;
        Block block = Block::None;
    };
    // Recomputed from the current layout.
    const std::vector<ChordSummary> &getChords() const { return summaries_; }
    void enumerateChords();

    // The scale Settings::maxWidth and ChordSummary::width are in: the mean
    // side length of the block structure as it was handed over.
    double getWidthScale() const { return widthScale_; }

private:
    // A chord with everything the collapse needs: the rungs in order, the two
    // ends of each labelled by which long side of the strip they belong to, and
    // the arcs of those two long sides.
    struct Chord {
        std::vector<int> faces;
        std::vector<int> rungs;             // faces + 1 of them, or faces if cyclic
        std::vector<int> endL, endR;        // per rung, its two ends
        std::vector<int> sideL, sideR;      // per face, the long side arcs
        bool cyclic = false;
        double width = 0.0;
        double length = 0.0;
        bool collapsible = false;
        Block block = Block::None;
    };

    void indexDarts();
    void classifyBoundary();   // which loop each node is on, -1 for the interior
    bool usable(int face) const;
    // Where `arc` sits in face `f`'s cycle, or -1.
    int sideOf(int face, int arc) const;
    // The face on the other side of `arc` from `f`, or -1 at the boundary.
    int across(int arc, int f) const;

    bool walkChord(int seedArc, Chord &out) const;
    void applyRules(Chord &c) const;
    bool collapse(const Chord &c);

    // How much a node's position is worth keeping: 2 a corner of the
    // polysquare, 1 anywhere else on the boundary, 0 in the interior.
    int priorityOf(int node) const;

    QuadLayout layout_;
    Settings settings_;
    Report report_;
    Block lastBlock_ = Block::None;

    std::vector<int> faceOfDart_;
    std::vector<int> loopOfNode_;
    std::vector<char> onBoundaryNode_;
    std::vector<Chord> chords_;
    std::vector<ChordSummary> summaries_;
    double widthScale_ = 1.0;
};

#endif // __CHORD_COLLAPSE_HXX__
