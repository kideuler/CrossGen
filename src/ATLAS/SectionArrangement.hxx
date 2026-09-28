#ifndef __SECTION_ARRANGEMENT_HXX__
#define __SECTION_ARRANGEMENT_HXX__

#include <functional>
#include <string>
#include <vector>

#include "mesh/Mesh.hxx"

// The arrangement of *sections* in a planar region bounded by polygonal loops:
// straight chords from the region's corners (and, optionally, its flat nodes)
// to the far boundary, and the faces they cut the region into. This is Sec.
// 7.2's "feature-aware stations" and Sec. 7.5's "graph of corridor patches ...
// junction cavities" made into one construction, and it is shared by the two
// places that need to agree on it:
//
//   * ExplicitTemplates' general sections, which fill each face on a carrier
//     (a grid for four corners, a star for three or five);
//   * CoarseDomain's count reconciliation, which builds the same arrangement
//     on the input's own curves *before* sampling them, so that the sides the
//     faces will want equal are sampled equally. An interface's subdivision
//     can never change afterwards (Sec. 8.4), so a face whose opposite sides
//     are interfaces sampled at different spacings could never be a grid, and
//     only sampling knows enough, early enough, to prevent it.
//
// ### The construction
//
// A loop vertex is a *node* when the layout must have a macrovertex there and
// a *corner* when, in addition, the region turns there. A corner of interior
// angle a wants q = round(a / quarter turn) cells, so it sends q - 1 sections
// at k a / q from its outgoing side (k = 1 .. q - 1): a reflex corner of 270
// degrees two, one of 150 degrees one, a convex corner none. A flat node -- a
// node the boundary runs straight through -- sends one, at a / 2, when
// Options::flatSections asks. Each section is cast to the first loop edge it
// meets and ends on a node within Options::nodeSnap of the hit, else on the
// nearer end of the hit edge; a section that would leave the region is
// refused (a corner's) or dropped (a flat node's).
//
// The arrangement's graph nodes are the corners, the section ends and the
// sections' crossings; its arcs are the loops between consecutive graph nodes
// and the sections between consecutive hits. Walking the half-edges with the
// region on the left gives the faces, and a face vertex is one of the face's
// corners when it is a region corner or the face turns by more than
// Options::cornerTolerance there. Every face must have three to five corners
// and each side one arc -- a node inside a side would be a T-junction.
class SectionArrangement {
public:
    struct Loop {
        std::vector<Point> X;        // walked with the region on the left
        std::vector<double> angle;   // the region's interior angle at each vertex
        std::vector<char> node;      // a macrovertex the layout must have
        std::vector<char> corner;    // ... at which the region also turns
        // A vertex a section passing close should end on, as on a node (the
        // coarse samples CoarseDomain placed where the input's sections
        // landed). Optional: empty means none.
        std::vector<char> preferred;
    };

    struct Options {
        bool flatSections = true;
        double nodeSnap = 0.15;          // fraction of the section's length
        double cornerTolerance = 0.4363; // radians (25 degrees)
        // A cross field (its angle at a point, and its coherence |u| there),
        // or unset. A section then leaves along the field direction nearest
        // the one its angle rule gives, when the field is coherent next to
        // the corner (coherence >= minCoherence) and that direction is within
        // maxTurn of the rule's and inside the corner: a section is a
        // separatrix of the layout, and the rule's bisector can run where no
        // separatrix would. multimat/geom013's plate has a corner of 150
        // degrees where the weld's root meets the buttering; its bisector
        // grazes the plate's bottom (a face corner of 163 degrees, a
        // T-junction), where the field, and every other method, run
        // straight to the far wall.
        std::function<double(const Point &, double *)> cross;
        double minCoherence = 0.3;
        double maxTurn = 0.6109;         // radians (35 degrees)
    };

    struct End {
        int l = -1, i = -1;              // loop, and position on it
    };
    struct Section {
        End a, b;                        // a is the corner or flat node it leaves
        bool shiftable = false;          // b was not snapped to a node
    };
    struct Node {
        int l = -1, i = -1;              // a loop position, or -1 for a crossing
        Point x;
    };
    struct Arc {
        int a = -1, b = -1;              // graph nodes
        int loop = -1;                   // loop arcs: which loop; -1 for a section arc
        std::vector<int> pos;            // loop arcs: positions a .. b
        int section = -1;                // section arcs: which section
        double length = 0.0;
    };
    struct Half {
        int arc = -1;
        bool rev = false;
        int from = -1, to = -1;
        double outAng = 0.0, inAng = 0.0;
    };
    struct Face {
        std::vector<int> halves;         // counter-clockwise
        std::vector<int> corners;        // indices into halves: the half-edge leaving each corner
    };
    struct Result {
        bool ok = false;
        std::string why;
        std::vector<Node> nodes;
        std::vector<Arc> arcs;
        std::vector<Half> half;
        std::vector<Face> faces;
    };

    // The sections the construction casts. False (with a reason) when a
    // corner's section leaves the region or finds nothing to end on.
    static bool cast(const std::vector<Loop> &L, const Options &o, std::vector<Section> &out, std::string &why);

    // The arrangement of the given sections, faces and all.
    static bool arrange(const std::vector<Loop> &L, const std::vector<Section> &S, const Options &o, Result &r);

    // The chord between two loop positions is inside the region: it leaves
    // each end inside that end's interior sector and meets no loop edge but
    // the ones at its ends.
    static bool clear(const std::vector<Loop> &L, End a, End b);

    // Positions i and j of loop l are neighbours along it (or equal).
    static bool adjacent(const std::vector<Loop> &L, int l, int i, int j);
};

#endif // __SECTION_ARRANGEMENT_HXX__
