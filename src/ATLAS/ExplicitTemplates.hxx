#ifndef __EXPLICIT_TEMPLATES_HXX__
#define __EXPLICIT_TEMPLATES_HXX__

#include <string>
#include <unordered_set>
#include <vector>

#include "ATLAS/CavityFill.hxx"
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
//   * Annuli (Sec. 7.4). Two loops, both radial graphs about a point of the
//     hole: four sectors, no centre, the hole kept exactly. When one loop
//     has exactly four protected corners -- a plate with a hole -- the
//     sector rays go through them.
//
// Each accepted template's corners become *designated* macrovertices of the
// carrier, so the rest of the pipeline can see the layout the template meant
// even where the carrier is regular across it (an annulus is all valence-4).
// A rejected template leaves the carrier untouched and records why (Sec. 13.2:
// "the reason each attempted coarse constructor was rejected").
class ExplicitTemplates {
public:
    struct Options {
        bool sections = true;
        bool stars = true;
        bool ogrids = true;
        bool annuli = true;
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
    bool tryStar(const Region &R, CavityFill::Patch &P, Attempt &A);
    bool tryOGrid(const Region &R, CavityFill::Patch &P, Attempt &A);
    bool tryAnnulus(const Region &R, CavityFill::Patch &P, Attempt &A);

    bool isNode(int v) const;
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
};

#endif // __EXPLICIT_TEMPLATES_HXX__
