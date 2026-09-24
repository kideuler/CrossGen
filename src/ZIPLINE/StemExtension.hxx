#ifndef __STEM_EXTENSION_HXX__
#define __STEM_EXTENSION_HXX__

#include <memory>
#include <vector>

#include "ZIPLINE/FieldTracer.hxx"
#include "ZIPLINE/QuadLayout.hxx"
#include "ZIPLINE/TriangleLocator.hxx"

// ---------------------------------------------------------------------------
// The optional post-pass of docs/viertel_2019.md Sec. 12 for the T-junctions
// Sec. 4 leaves behind: "for each remaining TJ, resume tracing its stem until
// the boundary; keep the result if the layout stays valid and #TJ drops."
//
// A T-junction is a separatrix that stopped on the side of another -- by
// condition 2, crossing the same separatrix a second time, or condition 3,
// meeting a singularity's own separatrix near it -- and a chord collapse can
// leave one behind when the chord that would have removed it failed the energy
// test (the paper's own account of four of its eight). The component whose
// side it stopped on then has a fifth node that is not a corner of it, which
// this codebase cannot mesh yet, so the component gets no elements
// (LayoutBlocks). Tracing the stem on through that side turns the junction
// into an ordinary crossing: the component is split in two along a streamline
// of the field, each half with the junction for a corner, and every separatrix
// the continuation crosses on its way to the boundary gets an ordinary
// four-valent node where it is crossed.
//
// Two things make a continuation worth keeping, and each is checked on the
// layout it would produce rather than assumed:
//
//   * it reaches the boundary crossing everything it meets squarely. A
//     continuation that meets a separatrix at a shallow angle is running along
//     it, which in the continuum is one curve; one that runs out of steps is
//     winding onto a limit cycle; one that runs into an existing node would
//     make that node something other than a crossing. Each of those is
//     refused, and the T-junction stays.
//   * the layout it gives is still a plane graph tiling the model, with one
//     T-junction fewer and no component that was four-cornered turned into one
//     that is not.
//
// One end other than the boundary is accepted: another T-junction nose to nose
// with this one (Settings::joinRadius). Two stems that stopped a fraction of an
// element apart on either side of a thin strip, pointing at each other, are one
// streamline; the continuation ends on the other junction and both become
// crossings.
//
// It is a pass over the layout after Sec. 4, not a change to the tracing:
// stopping conditions 2 and 3 are what make the collapsible chords Sec. 4
// exists for (the paper's Fig. 15), and tracing every separatrix to the
// boundary from the start would be the motorcycle graph the paper argues
// against. PartitionSimplify runs it once its collapses are done, when the
// layout it was handed came from a trace, and then collapses again: a
// continuation runs alongside existing separatrices as often as any
// separatrix does, and the strips it leaves are what Sec. 4 removes.
// ---------------------------------------------------------------------------
class StemExtension {
public:
    struct Settings {
        // Squarer than this is a crossing; shallower is running alongside.
        double minCrossingAngle = M_PI_4;
        // A continuation crossing more separatrices than this is refused: it
        // would re-cut most of the model to remove one T-junction.
        int maxCrossings = 64;
        // Triangles a continuation may cross before it is taken to be winding
        // onto a limit cycle.
        int maxSteps = 20000;
        // In mean edges: a continuation crossing the side another T-junction
        // stopped on this close to it, heading the way that one's stem leaves
        // it, ends on it instead of carrying on (extend()). Zero turns it off.
        double joinRadius = 1.5;
    };

    struct Report {
        int tJunctionsBefore = 0;
        int tJunctionsAfter = 0;
        int attempted = 0;
        int extended = 0;
        int joined = 0;          // ... of which ended on another T-junction
        // Why the rest were not.
        int noBoundary = 0;      // ran out of steps, or stuck
        int tangential = 0;      // met a separatrix at a shallow angle
        int throughNode = 0;     // ran into an existing node
        int tooManyCrossings = 0;
        int invalid = 0;         // the layout it gave failed the checks
    };

    StemExtension(const QuadLayout &layout, const FieldTracer &tracer, const Settings &settings);
    StemExtension(const QuadLayout &layout, const FieldTracer &tracer)
        : StemExtension(layout, tracer, Settings()) {}

    // Extend every T-junction's stem that can be, one at a time, each against
    // the layout the ones before it left.
    void run();

    const QuadLayout &getLayout() const { return layout_; }
    const Report &getReport() const { return report_; }

private:
    // One T-junction: its node and the dart of its stem, leaving the node.
    struct Junction {
        int node = -1;
        int stemDart = -1;
    };
    std::vector<Junction> findJunctions() const;

    // Try one; on success layout_ is replaced and true returned.
    bool extend(const Junction &j);

    QuadLayout layout_;
    const FieldTracer *tracer_ = nullptr;
    std::unique_ptr<TriangleLocator> locator_;   // built on the first extension
    Settings settings_;
    Report report_;
};

#endif // __STEM_EXTENSION_HXX__
