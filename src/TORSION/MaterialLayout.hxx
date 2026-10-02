#ifndef __MATERIALLAYOUT_HXX__
#define __MATERIALLAYOUT_HXX__

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

#include "MERIDIAN/Arrangement.hxx"
#include "MERIDIAN/Interfaces.hxx"
#include "TORSION/TORSION.hxx"
#include "mesh/Mesh.hxx"

// TORSION's per-material mode (docs/cf_flow_pipeline.md Sec. 15): the layout of
// a multi-material model as one layout per material region -- Stages 1 to 8 of
// this pipeline run on the region alone -- matched across the interfaces the
// regions share, and glued into one arrangement of S that Stages 9 to 11 take
// as they take any other.
//
// ### Why a region at a time
//
// What the whole-model route got wrong on the multi-material corpus came, one
// way or another, from a single fact: a cone inside an enclosed region has no
// route to dS that does not cross the interface around it, so its slit crosses
// that interface, and every stage from 3F to 6 then has an interface lying in
// two charts to cope with (TORSION.hxx, "What a multi-material model asks of
// Stages 3F to 4R"). A region laid out on its own has none of that. Its
// interfaces are its boundary, held on an axis by Q3 the way dS is, and its
// cuts end on that boundary without crossing anything. The problems are small
// -- concrete's fifteen regions are a 12606-triangle matrix and fourteen stones
// and shells of a few hundred each -- and independent, so they run side by
// side, and a region that fails fails alone.
//
// What makes it possible is that the cross field is the same field on both
// sides of an interface and tangent to it there. Restricted to a region it is
// aligned to the region's whole boundary, and a layout edge that meets an
// interface meets it at right angles from either side. Two regions' layouts can
// therefore always be joined along the interface between them; the only thing
// that can disagree is *where* each puts its layout vertices along it.
//
// ### Matching, and extending what is not matched
//
// A layout vertex on an interface -- a separatrix landing on it, or a corner of
// the region's layout there -- has to be a layout vertex of the region on the
// other side as well, or the face there has a T-junction on its side, is not a
// simple quadrilateral, and Stage 10 cannot mesh it. Two vertices close
// together on opposite sides are one vertex of the glued layout: the isolines
// of the two maps are matched and a separatrix on one side continues as the
// separatrix on the other. A vertex with no partner is extended: the neighbour
// is handed the point as an emitter -- a vertex of its dS from which a layout
// edge must leave (SubdomainLabels::Options::extraEmitters) -- and that edge is
// traced through the neighbour's own map, to wherever on the neighbour's
// boundary the map takes it, where it is a new layout vertex in its turn.
//
// Done after the fact, that cascades. An extension traced through a map that
// did not know it was coming ends wherever the map sends it, usually on the
// next interface, where nothing is waiting for it either; jelly_roll's two
// spiral strips come out as 934 patches that way, against the four its layout
// has. So the emitters are handed to the neighbour's Stage 5 and not only to
// its Stage 7. Gamma_topo is seeded from an emitter's ray as from any cone's
// separatrix (Sec. 3.3): an extension that nearly meets one of the neighbour's
// cones is joined to it and ends there, and a separatrix of the neighbour's own
// that lands beside an emitter is joined to the emitter and lands on it. The
// matching is not a separate step with a tolerance of its own. It is the
// paper's near-miss rule, applied across the interface by the region on the
// far side of it. The neighbour's layout changes, its own layout vertices move,
// and its other neighbours are handed the difference: a round re-runs Stages 5
// to 8 (TORSION::relayout) of every region whose emitters changed, until no
// region's do.
//
// Emitters are only ever added. Recomputed from scratch each round, the
// matching settles into flipping between two states on a third of the corpus;
// accumulated, it settles in two to six rounds, at the price of the odd
// extension a later round would no longer have needed.
//
// ### Gluing
//
// Once no region's emitters change, every layout vertex on an interface has a
// partner on the other side within an edge -- the emitter sits on a vertex of
// the mesh, the landing it answers anywhere along an edge -- and the two are
// one node of the glued arrangement, placed at the corner if either is one,
// else at the emitter, else half way between. The interface between two such
// nodes is one arc, with a patch of each region on either side of it. Each
// region's faces keep the corners and sides its own Stage 8 found for them.
// A face whose far side put a node where its own side has none gets that node
// as a T-junction on its side, is no longer simple, and is counted as such by
// the same check() any arrangement gets.
//
// ### The cone pairs
//
// Stage 1 cancels a +1/-1 pair the field put inside one material region, and
// runs its alignment-driven retry keeping them (TORSION::run()). A region on
// its own is asked to do the same (Options::cancelSingleMaterialDipoles). Where
// the layout that comes out still has a face that is not a simple
// quadrilateral, the region is laid out once more keeping the pairs, and the
// better of the two stands (Options::keepDipolesOnFailure): on dam, cancelling
// the pair in the spillway's apron leaves three faces without four corners,
// and keeping it leaves none.
//
// ### The flat corners
//
// A region's boundary is mostly interface, and a layer that pinches out
// between two interfaces has a corner of 20 to 45 degrees that the field reads
// as two quarter turns. Stage 1 keeps one at the corner and used to put the
// other on the vertex of the interface beside it, where nothing turns: a patch
// corner of pi, which Stage 10 meshes as an element with a node in the middle
// of a side (ConeSingularities::relocateFlatCones). Stage 1 now moves it to
// the corner at the end of that side, or inside the region. Where the
// region's layout cannot be made valid with it there, the region is laid out
// again with it left in place, and the better of the two stands
// (Options::keepFlatConesOnFailure; see retry() in the .cxx).
class MaterialLayout {
public:
    struct Options {
        // Rounds of matching at most, after the first layout of every region.
        int rounds = 10;
        // Regions laid out at once; 0 is one per hardware thread. Each worker
        // gives the linear solves inside a region its share of the rest.
        int threads = 0;
        // How far apart, in edges of the branch, two layout vertices on
        // opposite sides of an interface may be and still be matched into one
        // vertex of the glued layout. The match moves each end of the curves
        // concerned by up to that much along the interface.
        double matchTolerance = 3.0;
        // Re-lay a region out keeping its +1/-1 pairs when cancelling them
        // left a face that is not a simple quadrilateral.
        bool keepDipolesOnFailure = true;
        // Re-lay a region out with its +1 cones where Stage 1's rounding put
        // them, on a straight stretch of its boundary, when moving them off it
        // (TORSION::Options::relocateFlatCones) left the region's layout short
        // of Definition 2.1. The cone set that moved is better for the mesh --
        // the one that stayed is a patch corner of pi -- but a layout that is
        // not valid is worse than either. See firstLayout() in the .cxx.
        bool keepFlatConesOnFailure = true;
        // The interface network's kink angle (Interfaces::Options::kinkAngle).
        double kinkAngle = 0.7853981633974483;
        // Weigh, before the number of pairs, how close together the nodes of
        // the glued layout would land (matchBranch in the .cxx): two pairs that
        // both belong on one vertex are a collision, not two matches. Off is
        // the most-pairs alignment alone, as on 2026-09-29.
        bool spacedMatching = true;
        // Bend a matched separatrix into the glued node it now ends at, over
        // four times the distance its end moved, instead of moving its last
        // point alone (bendOnto in the .cxx).
        bool bendMatchedEnds = true;
    };

    // One material region of S, laid out on its own.
    struct Region {
        int material = 0;
        std::vector<int> triangles;         // triangles of S, in the region mesh's order
        std::vector<int> vertices;          // vertex of the region mesh -> vertex of S
        std::unordered_map<int, int> local; // vertex of S -> vertex of the region mesh
        std::shared_ptr<Mesh> mesh;
        // Stages 0 to 8 on the region, Stage 0 being the field restricted.
        std::unique_ptr<TORSION> layout;
        // Vertices of the region mesh its neighbours asked a layout edge of.
        std::vector<int> emitters;
        bool dipolesKept = false;
        bool flatConesKept = false;   // laid out again with Stage 1's flat +1s in place
        int runs = 0;        // Stages 0 to 8
        int relayouts = 0;   // Stages 5 to 8 again
        double seconds = 0.0;

        bool laidOut() const { return layout && layout->hasArrangement(); }
    };

    // A layout vertex a region put on an interface.
    struct Station {
        int region = -1;
        int node = -1;        // in the region's arrangement
        int branch = -1;      // of getInterfaces(); -1 at a node of the network
        double s = 0.0;       // arc length along the branch
        int vertex = -1;      // the vertex of S it sits on, or -1 inside an edge
        Point p{0.0, 0.0};
        bool corner = false;  // a corner of the region's layout
        bool emitter = false; // one of the region's emitters
    };

    struct Report {
        int regions = 0;
        int regionsLaidOut = 0;   // reached Stage 8
        int regionsValid = 0;     // Stage 6 valid, every face a simple quadrilateral
        int rounds = 0;
        bool converged = false;
        int emitters = 0;
        int stations = 0;         // layout vertices on the interfaces, both sides
        int matched = 0;          // pairs made one node
        int unmatched = 0;        // left single by the matching
        int extended = 0;         // faces split to carry those on to dS
        int tJunctionsLeft = 0;   // what the splitting could not carry on
        double maxMatchOffset = 0.0;   // in edges of the branch
        int dipolesKept = 0;
        // Regions whose Stage 1 moved a +1 off a straight stretch of their
        // boundary, and of those the ones laid out again without the move
        // because the layout with it was not valid.
        int flatConesMoved = 0;
        int flatConesKept = 0;
        // Layout edges asked of the regions in each round, in order: what a
        // settling matching looks like is this falling to zero.
        std::vector<int> askedPerRound;
        int regionRuns = 0;
        int regionRelayouts = 0;
        bool assembled = false;
        bool valid = false;
        double seconds = 0.0;
        std::vector<std::string> messages;
    };

    // What TORSION::run() hands this class from its own options: the mode's
    // settings, and the options each region's own TORSION runs with -- the
    // pipeline's, with the mode itself, Stages 9 to 11 and the disk templates
    // turned off, and Stage 1's cone-pair rule extended to a single material.
    // Static so that a caller driving the mode itself, as the viewer does,
    // runs it exactly as run() does.
    static Options optionsFor(const TORSION::Options &pipeline);
    static TORSION::Options regionOptionsFor(const TORSION::Options &pipeline);

    // `field` is one unit spin-4 value per triangle of `mesh` (DualMBO::
    // u_k_prev); `regionOptions` are what each region's own TORSION runs with,
    // less the field and the emitters, which are set here.
    MaterialLayout(std::shared_ptr<Mesh> mesh, const Eigen::VectorXcd &field,
                   const TORSION::Options &regionOptions, const Options &opts);
    ~MaterialLayout();

    // Lay every region out, match, extend, glue. True when the glued layout
    // is valid and every region's Stage 6 reached Definition 2.1.
    bool run();

    const std::vector<Region>& getRegions() const { return regions; }
    // The interface network the matching is written against: the model's,
    // with closed loops left whole (Interfaces::Options::splitLoops off),
    // because a loop has no node for a split to stand on.
    const Interfaces& getInterfaces() const { return *network; }
    // The layout vertices on the interfaces as the last round left them.
    const std::vector<Station>& getStations() const { return stations; }
    // The stations the final pairing left without a partner: a T-junction on
    // the far side of each, indices into getStations().
    const std::vector<int>& getUnmatched() const { return unmatched; }
    bool hasArrangement() const { return arrangement != nullptr; }
    const Arrangement& getArrangement() const { return *arrangement; }
    std::unique_ptr<Arrangement> takeArrangement() { return std::move(arrangement); }
    const Report& getReport() const { return report; }

private:
    // One branch of the network as the matching reads it: its vertices in
    // order, arc length along them, and the region either side, left being the
    // side the model is on when the branch is walked in the order of `verts`.
    struct Branch {
        std::vector<int> verts;
        std::vector<double> cum;
        std::unordered_map<int, int> at;   // interior vertex of S -> index in verts
        double length = 0.0;
        double longestEdge = 0.0;
        bool closed = false;
        int left = -1, right = -1;
    };
    // Two stations made one, or one on its own. `sequence` is every node of the
    // glued layout on the branch in the order they meet along it -- (i, j) for
    // a pair, (i, -1) for a station left single -- starting, on a closed loop,
    // after the widest gap between two stations.
    struct Pairing {
        std::vector<std::pair<int, int>> pairs;   // indices into the station list
        std::vector<int> single;
        std::vector<std::pair<int, int>> sequence;
        double origin = 0.0;   // where the sequence starts, as arc length
    };

    void buildRegions();
    void buildBranches();
    // Stages 0 to 8 of one region with its current emitters; `keepDipoles`
    // turns Stage 1's cancellation off, `keepFlatCones` its move of the +1s on
    // a straight stretch of boundary.
    void layOut(Region &r, bool keepDipoles, bool keepFlatCones);
    // The first layout of one region, or its relayout with new emitters, with
    // the cone-pair fallback of the class comment and then the flat-cone one
    // (Options::keepFlatConesOnFailure); retry() is the two fallbacks.
    void firstLayout(Region &r);
    void nextLayout(Region &r);
    void retry(Region &r);
    // Run `work` on each of `jobs`, `threads` at a time.
    void parallel(const std::vector<int> &jobs, const std::function<void(int)> &work);

    std::vector<Station> collectStations() const;
    // Which stations on branch b are one node with which: an alignment of the
    // two sides' sequences along the branch. See MaterialLayout.cxx.
    Pairing matchBranch(int b, const std::vector<int> &onBranch,
                        const std::vector<Station> &st) const;
    // The alignment of 2026-09-29: most pairs, then least displacement, with
    // no regard to where the nodes land (Options::spacedMatching off).
    Pairing matchBranchMostPairs(int b, const std::vector<int> &onBranch,
                                 const std::vector<Station> &st) const;
    // The emitters each region is asked for by the stations as they stand.
    std::vector<std::vector<int>> wantedEmitters(const std::vector<Station> &st) const;
    // Glue the regions' arrangements along the final pairing.
    bool assemble(const std::vector<Station> &st);
    // Extend what the matching left single through the glued faces: split a
    // face with a node in the middle of a side by the curve of its own Coons
    // parameterisation through that node, and carry the split on across the
    // face beyond, until dS. Returns the number of splits made.
    int extendThroughFaces(Arrangement::Assembly &out);

    // The edge {u, w} of branch b as the index of its first vertex in
    // `verts`, or -1; how far along the branch a point of edge `segment` is;
    // and the point at arc length s.
    int segmentOf(int b, int u, int w) const;
    double along(int b, int segment, const Point &p) const;
    Point pointAt(int b, double s) const;
    // An edge of S a boundary arc of region R runs along.
    int coveredEdge(const Region &R, const Arrangement::Arc &a) const;
    bool isRegionCone(const Region &r, int vertexOfS) const;

    std::shared_ptr<Mesh> mesh;
    Eigen::VectorXcd field;
    TORSION::Options regionOptions;
    Options options;
    std::unique_ptr<Interfaces> network;
    std::vector<Region> regions;
    std::vector<Branch> branches;
    std::vector<std::vector<int>> vertexRegions;          // vertex of S -> regions around it
    std::unordered_map<long long, int> edgeOf;            // vertex pair of S -> edge
    std::vector<Station> stations;
    std::vector<int> unmatched;
    std::unique_ptr<Arrangement> arrangement;
    Report report;
};

#endif // __MATERIALLAYOUT_HXX__
