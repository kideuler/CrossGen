#ifndef __PATTERN_REPAIR_HXX__
#define __PATTERN_REPAIR_HXX__

#include <string>
#include <vector>

#include "ZIPLINE/ConformalField.hxx"
#include "ZIPLINE/LayoutBlocks.hxx"
#include "ZIPLINE/QuadLayout.hxx"
#include "ZIPLINE/SeparatrixTrace.hxx"
#include "mesh/Mesh.hxx"

// ---------------------------------------------------------------------------
// Stage 5b: correcting the layout by editing the singularity pattern and
// solving the corrected field again -- docs/quad_layout_pipeline.md Secs. 9.4,
// 12.2 and 12.3, on the flexibility P2 (arXiv:2103.02939) gives as its reason
// for a second field at all: "the possibility of moving singularities, adding
// or removing some, while still respecting topological constraints. In this
// work, we use that advantage for correcting or adapting the quad layout when
// necessary, e.g., to fix non-quadrilateral partitions or to avoid limit
// cycles" (P2 Sec. 1).
//
// ### What P2 itself does, and what this does instead
//
// P2's Sec. 6.3 corrections are written as edits of the *layout*:
//
//   * a T-junction is removed by merging its stem with the nearest node of
//     the bar that is neither on dS nor another T-junction, "in the parametric
//     space" of the patches either side (Fig. 17) -- a straight segment in
//     (u, v), not a field line, and nothing re-solved or re-traced;
//   * a three-sided patch left by a valence-2 corner is split by "forcing the
//     singularity of valence 2 on the boundary into splitting in two
//     singularities of valence 3": the corner becomes valence 3 and a new
//     valence-3 point goes to the barycentre of the patch, which is split into
//     three quadrilaterals (Fig. 18). Whether the field is solved again with
//     the new point the text does not say; the figure draws the new point with
//     two of its three arms straight through it, which is what tracing the
//     *old* field from a regular point gives (a re-solved valence-3 point has
//     its arms 120 degrees apart).
//
// 3-5 pairs appear in P2 only as a user's edit before the second solve
// (Sec. 4, Fig. 7: a size map), never as a correction.
//
// Here every correction is a pattern edit followed by the whole of Stages 1b
// to 5 again: the corrected field re-solved with the edited pattern (two
// linear problems, cheap), re-traced, laid out, simplified and read as blocks.
// Re-solving is not optional for the case that matters most. On the MAMBO
// corpus 25 of the 32 models short of full coverage were one shape: a cusp, a
// zero-degree corner where an arc meets a line or another arc tangentially,
// which Table 1 calls a valence-3 corner (one quad) and the first field turns
// as valence 2 (none). The pattern then fails Poincare-Hopf by a quarter per
// such corner, the second field cannot be built at all, the trace falls back
// to the first field, and a triangle is left beside the cusp: exactly P2's
// Fig. 18. Adding the valence-3 point to the triangle makes the pattern valid
// with Table 1's corner, so it is the edit that makes the second field exist.
//
// ### The edits proposed (propose())
//
//   1. Poincare-Hopf deficit (ConformalField::Report::regionDeficit4): a
//      region -- a component of the model cut along its interfaces -- whose
//      pattern is m quarters short of chi (m < 0: over). Either a 3 at the
//      barycentre of each of m of the region's refused faces, three-cornered
//      ones first (P2 Sec. 6.3, Fig. 18; 5s and pentagons when m < 0), or --
//      the first field's other way of misreading a corner, a 5 put against a
//      corner it then turned as flat -- the region's m 5s nearest those faces
//      removed. The candidates are combinations of a few m-subsets per region
//      (m is 1 or 2 in practice). On a multi-material model the shortfall is
//      the interface version of a cusp: a sector of a node sharper than 45
//      degrees, which the first field turns as no quad at all and the
//      corrected field is given one (ConformalField's constructor).
//   2. A 3-5 pair at a T-junction (Sec. 9.4: "a 3-5 pair placed on the path of
//      the cycle"). Most T-junctions left after Sec. 4 are cut limit cycles,
//      separatrices that wind round a hole and come back an offset away from
//      where they left; a dislocation on the path absorbs the offset and lets
//      the loop close onto the pair. Six placements, measured to win on
//      different models, so all are tried: the 3 on the T and the 5 inside the
//      face that runs straight past it, or the other way round; and the pair
//      on the line the stem would continue along into that face, the 3 or the
//      5 nearer the T, at 0.3 and 0.7 of the way across the face or at 2 and 5
//      mean edges for a face too small for that.
//   3. P2's Fig. 17, along the field: move the stem's own singular point along
//      the bar by the offset from the T to a neighbouring node of the bar
//      (scaled 1, 1/2 and 2, since the stem does not move rigidly with it), so
//      the stem runs into that node's line instead of stopping short of it.
//   4. Remove a +- pair of the pattern near a refused face (a dipole the first
//      field had not finished annihilating, too far apart for
//      ConformalField::Options::cancelDipoles), or move one singularity near a
//      refused face into it.
//
// ZIPLINE::repairPattern() evaluates every candidate, keeps the best (sound,
// more of the model covered, then fewer refused faces, then fewer blocks),
// checks its mesh, and goes round again; the class here only proposes.
//
// ### Measured (docs/viertel_2019.md Sec. 14, "Fifth pass")
//
// h = 0.05, TMOP 200 sweeps, against the pipeline without this stage:
// singlemat mean coverage 99.67% -> 99.99% (31 -> 33 of 35 models full,
// T-junctions left 10 -> 2), MAMBO 99.71% -> 100.00% (375 -> 406 of 407 full,
// the one inverted mesh gone), every layout sound, the mean worst element
// after smoothing unchanged on singlemat (0.817 -> 0.814) and better on MAMBO
// (0.854 -> 0.862). Each edit alone does a part of it: 1 is the cusps (23 of
// the 25 to full coverage), 2 and 3 the T-junctions, about equally often. On
// the multi-material corpus, with the corrected field solved per region
// (ConformalField's header), 96.11% -> 99.96% (12 -> 23 of 31 exactly full,
// 29 at 99.95% or more), its three inverted meshes gone: there edit 1 is what
// builds the corrected field on nine of the ten models whose regions fail
// Poincare-Hopf at sharp interface sectors.
// ---------------------------------------------------------------------------
class PatternRepair {
public:
    enum class Kind {
        Deficit,          // 1: a 3 per triangle (or 5 per pentagon) for the deficit
        DeficitRemoval,   // 1: the opposite point removed instead
        PairAtTJunction,  // 2
        MoveStem,         // 3
        CancelPair,       // 4
        MoveIntoFace      // 4
    };

    struct Options {
        // Rounds of the search; 0 turns Stage 5b off.
        int rounds = 8;
        bool deficit = true;            // edit 1
        bool pairsAtTJunctions = true;  // edit 2
        bool moveStems = true;          // edit 3
        bool moveIntoFaces = true;      // edit 4 (both kinds)
        // An edit is kept for more than this fraction of the model covered,
        // or as much covered with fewer refused faces.
        double minGain = 1e-4;
        // Before an edit is kept its blocks are meshed (BlockQuadMesh at
        // ZIPLINE::Options::mesh) and smoothed for checkSweeps TMOP sweeps,
        // and the edit refused if an element is left inverted or the worst
        // scaled Jacobian is below min(qualityFloor, half the unrepaired
        // layout's). Without it the search buys coverage with degenerate
        // elements: S0_11 went to full coverage with a worst scaled
        // Jacobian below 0.001 (0.63 unrepaired).
        bool checkMesh = true;
        int checkSweeps = 200;
        double qualityFloor = 0.15;
        // Every candidate costs a run of Stages 1b-5 -- mostly Sec. 4's
        // simplification -- so the search is budgeted: at most this many per
        // round and per model, the deficit edits first (they are what lets
        // the corrected field be built at all) and the rest by the uncovered
        // area of the face they are made for. Unbudgeted, multimat/bubbles
        // (82 singularities, ten inclusions) tried 1089 candidates and took
        // eleven minutes for 93.3% -> 99.7%.
        int maxCandidatesPerRound = 64;
        int maxCandidates = 256;
    };

    struct Edit {
        Kind kind = Kind::Deficit;
        std::string what;                         // for the report
        // Order within the budget (Options::maxCandidatesPerRound): larger
        // first. A deficit edit is first of all; any other, the area of the
        // refused face it is made for.
        double priority = 0.0;
        std::vector<int> remove;                  // indices into the pattern
        std::vector<ConformalField::Singular> add;
    };

    // The candidate edits of `pattern` -- the singularities Stage 1b was
    // given, its own dipole and corner edits not applied -- for the faces of
    // `layout` (the simplified one) that `blocks` refused. `declined` is the
    // corrected field this pattern did not build (its per-region
    // Poincare-Hopf deficits are read off it), null when it built; `trace`
    // is what the layout was built from (for the stems' origins).
    static std::vector<Edit> propose(const Mesh &mesh, const std::vector<ConformalField::Singular> &pattern,
                                     const ConformalField *declined, const QuadLayout &layout,
                                     const LayoutBlocks &blocks, const SeparatrixTrace &trace, const Options &opts);

    // `pattern` with `edit` applied.
    static std::vector<ConformalField::Singular> apply(const std::vector<ConformalField::Singular> &pattern,
                                                       const Edit &edit);

    static const char *kindName(Kind k);
};

#endif // __PATTERN_REPAIR_HXX__
