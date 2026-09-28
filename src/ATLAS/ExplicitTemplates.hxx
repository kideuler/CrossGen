#ifndef __EXPLICIT_TEMPLATES_HXX__
#define __EXPLICIT_TEMPLATES_HXX__

#include <string>
#include <unordered_set>
#include <vector>

#include "ATLAS/CavityFill.hxx"
#include "ATLAS/ReferenceField.hxx"
#include "ATLAS/SquareCarrier.hxx"

// Stage 3 of docs/square_transport_2d_theory_and_implementation.md: attempt
// large explicit replacements.
//
// Sec. 8.1 is the reason this stage exists. Every triangle centroid of the
// carrier is a three-valent interior vertex, and a three-valent vertex cannot
// sit inside a rectangular grid patch, so no amount of grouping the initial
// cells (Stages 4 and 5) recovers a blocking coarser than the triangulation.
// Topology has to be *replaced*, and the cheapest replacement is the biggest:
// a whole material region at once, rebuilt from its boundary alone.
//
// ### The cavity is a whole region
//
// A region is a maximal connected set of cells of one material; its boundary
// is domain boundary and material interface. Replacing all of it means the
// only reconciliation Sec. 8.4 asks for is along that boundary, and there the
// rule of CavityFill applies: retained edges are literally retained, and the
// only change allowed is inserting collinear points on a *domain* boundary
// segment, which has no neighbour to disagree. An interface segment is never
// subdivided, because the region on its other side would be left with a
// hanging node; a template that would need it is rejected.
//
// ### The families (Sec. 7)
//
//   * Sectioned quadrilaterals (Secs. 7.1, 7.2). A region with four protected
//     corners and no reflex one is one block. A region with reflex corners is
//     cut by *sections*: straight chords from each reflex corner, splitting
//     its angle into near-quarter sectors, run to the far boundary (the
//     feature-aware stations of Sec. 7.2 -- a section is inserted at every
//     protected corner that would otherwise sit inside a block side). The
//     faces of the arrangement must each have four corners; then the
//     opposite-side counts are reconciled by the union-find of Sec. 11.2 and
//     each face is filled with its own structured grid. An L is three blocks,
//     a bent corridor three, a T five. A convex region with an even number
//     K >= 6 of corners is cut by chords between corners instead, a fan of
//     K/2 - 1 four-sided faces.
//   * Stars. A convex three- or five-cornered region is K blocks round one
//     centre of valence K -- for K = 3 exactly the three-quad split of Sec. 3
//     one level up, which is Sec. 13.3's "Triangle" row. The side counts fix
//     the spoke counts through n_i = s_{i-1} + s_{i+1}; when that has no
//     positive integer solution, domain-boundary sides are subdivided until
//     it does.
//   * O-grids (Sec. 7.3). A star-shaped region with at most two protected
//     corners: a core quadrilateral about a kernel point c and four shells,
//     every shell node on the ray from c through a boundary vertex. Its
//     certificate is geometric and complete: with c strictly inside the
//     kernel and the core convex, every shell cell is the intersection of a
//     wedge of angle < pi with the far side of one chord and the near side of
//     another, hence convex, and every core cell is the intersection of two
//     strips between non-crossing segments of a convex quadrilateral, hence
//     convex.
//   * Half O-grids. Not one of Sec. 7's families, and needed because of what
//     Sec. 7.3 does to a region whose two protected corners are real corners:
//     it has to make them two of its four split points, and a split point is
//     where two shells meet, so each corner is cut into two cells. On
//     singlemat/geom002, a half disk, that meshed as cells of 27 and 61
//     degrees in each 90-degree corner, and Sec. 4's count then asks the
//     interior for four three-valent vertices, where one cell at each corner
//     would leave it two. A region with two
//     corners joined by a straight side -- every body that crosses the axis
//     of an axisymmetric (r, z) model -- is instead half of an O-grid: mirror
//     it across that side, lay Sec. 7.3's O-grid over the mirror image with its
//     centre on the side and its split points in mirror pairs, and cut it back
//     along the side. What is left is a core resting on the side and three
//     shells round it, four blocks, the two corners each the corner of one
//     side shell, and the side's two outer stretches radial lines of those
//     shells. Sec. 7.3's certificate carries over unchanged: every shell cell
//     is X(s, t) = c + lambda d(s) with lambda increasing in t, so
//     det DX = lambda_t det(q', d) > 0 whether the radial spacing is uniform
//     (the rays in the arc's interior) or the side's own points (the rays
//     along it). The counts close only if the core's side on the straight
//     side has as many edges as the arc opposite it, which a straight dS side
//     meets by Sec. 8.4's inserted points and an interface meets only by
//     widening the core.
//   * Annuli (Sec. 7.4). Two loops, both radial graphs about a point of the
//     hole: four sectors, no centre, the hole kept exactly. When one loop
//     has exactly four protected corners -- a plate with a hole -- the
//     sector rays go through them.
//   * General sections, Sec. 7.5's "graph of corridor patches ... junction
//     cavities" made explicit, for regions that touch an interface. Sections
//     run from the reflex corners of *every* loop, a hole's included, and
//     from each flat node (the macro edge a neighbour's template ends there
//     continues into this region), to the far boundary; the faces they cut
//     are grids when they have four corners and stars when they have three
//     or five, the star's centre being the singularity the region needs.
//     multimat/det_rocket's fuel region, a stadium round a square inclusion,
//     is eight faces this way -- six grids beside the square's sides and
//     corners, two three-cornered stars in the rounded end -- which is the
//     layout TORSION's field finds, singularities where its cones are. The
//     counts are Sec. 11.2's union-find over the grids plus the stars' own
//     equations n_k = s_{k-1} + s_{k+1}, solved by a small search over the
//     counts nothing fixes. Interface sides are fixed on both sides of a
//     face (Sec. 8.4), so where two of them must match and do not, the
//     sections' far ends are moved along the boundary a few vertices until
//     they do.
//
// Each accepted template's corners become *designated* macrovertices of the
// carrier, so the rest of the pipeline can see the layout the template meant
// even where the carrier is regular across it (an annulus is all valence-4).
//
// ### Corners and flat nodes
//
// Every protected or designated vertex on a region's boundary is a *node*: a
// macrovertex the layout must have. Not every node is a *corner* of the
// region's template. Where the boundary runs straight through a node -- a
// star's split point that a template across a smooth interface designated, or
// a T-junction of the interface network seen from its flat side -- the
// region's blocks meet there side by side, and making it a block corner puts
// a corner of nearly pi into a block (multimat/geom001's quarter disk, whose
// neighbour's five-block star had designated a point on their arc, was read as
// a four-cornered region, and multimat/geom012's mat 1 as a five-cornered one
// that no family took). So a node within cornerTolerance of pi, and any
// designated vertex that is not protected (it lies on a smooth feature curve,
// whatever the coarse polygon's own angle there), is not counted as a
// corner. It stays a boundary vertex of the template, where every family
// runs a grid line from it into the region, so Stage 5 still finds the
// macrovertex it needs; a star does better and puts its spoke there
// (tryStar), so that the two sides' layouts meet at one point instead of
// each cutting the other's blocks.
// A rejected template leaves the carrier untouched and records why (Sec. 13.2:
// "the reason each attempted coarse constructor was rejected").
class ExplicitTemplates {
public:
    struct Options {
        bool sections = true;
        bool stars = true;
        bool ogrids = true;
        bool halfOGrids = true;
        bool annuli = true;
        // A half O-grid gives each of its two corners one cell, so it is tried
        // only when both are below this interior angle (degrees). Round it to
        // the nearest quarter turn: a corner past 135 degrees wants two cells,
        // which is what the full O-grid gives it.
        double halfOGridCorner = 135.0;
        // A protected corner with an interior angle past pi by more than this
        // (degrees) is reflex and emits sections.
        double reflexAngle = 20.0;
        // A face angle within this of pi (degrees) is not a block corner.
        double cornerTolerance = 25.0;
        // A section ends on a protected corner if it passes within this
        // fraction of its own length of it.
        double nodeSnap = 0.15;
        // O-grid core corners sit this fraction of the way from c to the
        // boundary split points.
        double coreFraction = 0.5;
        // Scaled-Jacobian floor for an accepted template cell.
        double minScaledJacobian = 0.02;
        // Laplacian sweeps allowed to repair a Coons fill that does not
        // certify as placed. The explicit maps (O-grid, annulus) are never
        // smoothed: they are the certificate.
        int smoothingIterations = 60;
        // Sec. 8.4: may subdivide domain-boundary segments to reconcile counts.
        bool splitBoundary = true;
        // A node the region's boundary runs straight through is not one of
        // its corners (see "Corners and flat nodes" above): a neighbouring
        // template's split point on a smooth interface, or a T-junction seen
        // from its flat side, within cornerTolerance of pi. Off: every
        // protected or designated vertex is a corner, as before 2026-09-28.
        bool flatNodes = true;
        // Next to an interface, a domain-boundary side may also shed its
        // droppable points (the coarse carrier's midpoints) when that is what
        // closes a section's or a star's counts against the interface sides.
        bool dropBoundary = true;
        // "General sections" (see the class comment), for a region that
        // touches an interface: sections from the reflex corners of every
        // loop -- holes included -- and from flat nodes, faces of three and
        // five corners filled as stars, and section ends moved up to
        // sectionShift loop vertices when that is what closes two interface
        // counts that must match.
        bool generalSections = true;
        bool flatSections = true;
        int sectionShift = 3;
        // Vertices of the domain a section passing close should end on, as on
        // a node, indexed by domain vertex: CoarseDomain's samples at the
        // points where the input's own sections landed (see its "Counts
        // across an interface"). Null: none.
        const std::vector<char> *preferredEnds = nullptr;
        // Next to an interface, a grid face designates the point facing each
        // node on its sides: the separatrix from the node runs straight
        // across the grid and Stage 5 cuts there anyway, and marking where it
        // leaves lets the template across the next interface see it as a
        // flat node (multimat/icf: the capsule's O-grid split points, carried
        // out through two one-block shells to the region round them).
        bool propagateNodes = true;
        // The order regions are tried in after those bounded by dS alone:
        // largest first (the default), or smallest first, so that small
        // inner templates -- a capsule's O-grid -- set the split points the
        // large regions round them follow. ATLAS runs both on
        // multi-material domains.
        bool smallestFirst = false;

        // The reference cross field (docs/atlas_crossfield_guidance.md, Sec.
        // 4.2), or null for none. With it every attempt records the region's
        // cones against the template's singular vertices, and:
        //   * conePlacement: where a family has an obvious slot for the
        //     cones, they fill it -- an O-grid's split rays and core corners
        //     on a region with exactly four +1/4 cones, a half O-grid's two
        //     core corners on one with two, a star's centre on the one cone
        //     of its sign (+1/4 for three blocks, -1/4 for five). Geometry
        //     alone leaves these to the principal axes, which on a disk are
        //     degenerate and pick the orientation by rounding;
        //   * chooseByField: every family that applies is built and the one
        //     with the least blocks + wDir E_dir + E_sing over its cells is
        //     committed, instead of the first that validates in the fixed
        //     order below. The half O-grid's dispatch rule is the one case of
        //     this that was written by hand.
        const ReferenceField *field = nullptr;
        bool conePlacement = true;
        bool chooseByField = true;
        double wDir = 10.0;
        ReferenceField::SingularityWeights singularity;
    };

    struct Attempt {
        int region = -1;
        int cells = 0;
        int material = 0;
        int corners = 0;
        int holes = 0;
        std::string family;
        bool accepted = false;
        std::string reason;
        int blocks = 0;
        int newCells = 0;
        int splits = 0;
        double minScaledJacobian = 0.0;
        std::string maps;   // bilinear / ruled / Coons / radial / polar, per block
        // Against the reference field, when there is one: the region's cones,
        // the template's singular vertices (inside, and defects on dS), its
        // E_dir (over the domain's area) and E_sing, and blocks + wDir E_dir
        // + E_sing. A count or sign mismatch here flags a template before
        // anyone looks at a mesh.
        bool scored = false;
        int conesPlus = 0, conesMinus = 0;
        int singPlus = 0, singMinus = 0;
        int edgePlus = 0, edgeMinus = 0;
        double eDir = 0.0, eSing = 0.0, score = 0.0;
        bool conePlaced = false;
        // chooseByField: the other families that validated, with their scores.
        std::string alternatives;
    };

    struct Report {
        int regions = 0;
        int templated = 0;
        int untouched = 0;
        int cellsBefore = 0, cellsAfter = 0;
        int boundarySplits = 0;
        int blocks = 0;               // blocks the accepted templates are made of
        double templatedArea = 0.0;   // fraction of the domain
        std::vector<Attempt> attempts;
    };

    ExplicitTemplates(SquareCarrier &carrier, const Options &opts);

    const Report &getReport() const { return report_; }

private:
    struct Region {
        std::vector<int> cells;
        std::unordered_set<int> set;
        int material = 0;
        int key = -1;           // smallest source triangle: stable across rebuilds
        CavityFill::Boundary B;
        int outer = -1;
        std::vector<int> holes;
        double area = 0.0;
        double h = 0.0;         // mean boundary edge length
        bool touchesInterface = false;
    };

    std::vector<Region> extractRegions() const;
    bool attempt(const Region &R, Attempt &A);
    bool trySections(const Region &R, CavityFill::Patch &P, Attempt &A);
    bool trySectionsGeneral(const Region &R, CavityFill::Patch &P, Attempt &A);
    bool tryStar(const Region &R, CavityFill::Patch &P, Attempt &A);
    bool tryOGrid(const Region &R, CavityFill::Patch &P, Attempt &A);
    bool tryHalfOGrid(const Region &R, CavityFill::Patch &P, Attempt &A);
    bool tryAnnulus(const Region &R, CavityFill::Patch &P, Attempt &A);

    // A built, certified template against the field: fills the Attempt's
    // scored fields on a copy of the carrier with it applied.
    void score(const Region &R, const CavityFill::Patch &P, Attempt &A) const;

    bool isNode(int v) const;
    // Node i of R's loop `loop` is one of the template's corners (see
    // "Corners and flat nodes").
    bool isCorner(const Region &R, int loop, int i) const;
    // Options::propagateNodes, for one grid face of R's template.
    void propagateNodes(const Region &R, CavityFill::Patch &P, const std::array<std::vector<int>, 4> &sides) const;
    // The arc's ids with `extra` collinear points inserted on its domain
    // boundary segments, longest first; empty when that is not allowed.
    std::vector<int> splitArc(CavityFill::Patch &P, const std::vector<int> &ids, int extra) const;
    bool arcSplittable(const std::vector<int> &ids) const;
    std::string classifyGrid(const CavityFill::Patch &P, const std::array<std::vector<int>, 4> &sides) const;

    SquareCarrier &C_;
    Options opts_;
    Report report_;
    int group_ = 0;
    // Chords between protected corners (outer-loop positions) that
    // trySections adds to its reflex sections: how a convex region with an
    // even number K >= 6 of corners is cut into K/2 - 1 four-sided faces.
    std::vector<std::pair<int, int>> chords_;
    // The cones of the region being attempted (empty without a field).
    std::vector<ReferenceField::Singularity> regionCones_;
};

#endif // __EXPLICIT_TEMPLATES_HXX__
