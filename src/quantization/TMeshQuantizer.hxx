#ifndef _TMESH_QUANTIZER_HXX_
#define _TMESH_QUANTIZER_HXX_

#include "quantization/QuantTMesh.hxx"

#include <vector>

// ─── QGP quantization, interval-assignment variant ──────────────────────────
//
// Section 6 of Campen, Bommes & Kobbelt, "Quantized Global Parametrization"
// (SIGGRAPH Asia 2015), with the Sec. 10.1 modification: every edge length
// is kept >= minLength (default 1) instead of >= 0. With all lengths
// positive every path in the T-mesh has positive parametric length, so no
// two nodes can collapse in parameter space and the Sec. 5.2 validity
// search is never needed -- Stage I terminates in a provably valid state
// and Stage II preserves validity with a plain scan for x < minLength.
//
// The algorithm navigates the solution monoid of the consistency system
// Ax = 0, x >= 0 by adding and subtracting generating vectors (Hilbert
// basis elements), never materializing the basis. A generating vector is an
// elementary circuit -- or, on a domain with boundary, a source-to-sink
// path -- of the digraph D of Sec. 6.2, whose nodes are the pairs
// (edge, incident row): "x[edge] was just incremented and this row is now
// unbalanced by it". An arc goes to every edge of opposite sign in that
// row, paired with that edge's other row; landing on a PHANTOM row is a
// sink, since the outer face never needs balancing. Geometrically such a
// circuit or path is a quadrilateral strip through the T-mesh whose edge
// lengths can all be incremented together without breaking consistency.
//
// D is never built: the Dijkstra searches enumerate arcs on the fly from
// the row references stored on the edges. Weights live on D-nodes (i.e. on
// T-mesh edges) and use the three-tier scheme of Sec. 6.2, tiers separated
// by factors of eta = #edges, so a strip through an edge already longer
// than ideal is chosen only when nothing cheaper closes up.
//
// Stage I grows the trivial (all-zero) solution until every edge has
// x >= minLength and every requirement set reaches its bound. Stage II then greedily adds or subtracts strips while
// the objective sum (x/xIdeal - 1)^2 improves. With xIdeal = 1 everywhere
// (the block-decomposition setting) Stage II drives every edge toward its
// minimum, acting as an automatic block-merging operator. The loop always
// holds a valid quantization, so a move budget is safe to impose.

class TMeshQuantizer {
public:
    struct Requirement {
        std::vector<int> edges;
        int atLeast = 1;
    };

    struct Options {
        int minLength = 1;       // lower bound kept on every edge length
        int maxStage2Moves = -1; // accepted Stage II moves; < 0 means no cap
        // Sets of edges whose lengths must sum to at least `atLeast`. With
        // minLength = 0 these carry the validity instead of the per-edge
        // bound: an edge may collapse as long as no set does. This is the
        // linear surrogate of Prop. 1 in Cheny et al., "Planar quad meshing
        // with guarantees" (docs/medial_axis.md Step 9) -- every radius chain,
        // boundary feature and boundary loop of the T-mesh gets one.
        std::vector<Requirement> requirements;
        // What Stage II minimizes. Relative is QGP's Sec. 6.3 measure,
        // sum (x / xIdeal - 1)^2, which weighs an edge by the inverse square
        // of its ideal length. Absolute is sum (x - xIdeal)^2, the target fit
        // of docs/medial_axis.md Step 9 (Option B's Delta), which charges
        // every edge in units of the element size alike -- so an edge far
        // shorter than an element that validity keeps at 1 costs about 1,
        // not the inverse square of its length.
        enum class Objective { Relative, Absolute };
        Objective objective = Objective::Relative;
    };

    struct Report {
        int stage1Vectors = 0;  // strips added to leave the trivial solution
        int stage2Moves = 0;    // accepted Stage II improvements
        int stage2Tried = 0;    // tentative Stage II applications
        // Edges no generating vector passes through. Every solution of the
        // consistency system assigns those zero, so the lower bound cannot
        // be met there by any algorithm -- it is a defect of the input
        // T-mesh, not a failure of the search. They are left at zero and
        // the rest is quantized around them.
        int forcedZeroEdges = 0;
        // Requirement sets no generating vector reaches: like forced zeros,
        // a defect of the input, left at zero.
        int unmetRequirements = 0;
        double objective = 0.0;   // of the kind Options::objective names
        bool consistent = false;
    };

    // The T-mesh must be finalized; results land in tmesh.edges[i].x.
    explicit TMeshQuantizer(QuantTMesh &tmesh) : TMeshQuantizer(tmesh, Options()) {}
    TMeshQuantizer(QuantTMesh &tmesh, Options opts);

    Report run();

private:
    // D-node ids: node 2 * e + k means x[e] was incremented and
    // edges[e].row[k] is the row left unbalanced.
    int rowOfNode(int n) const;

    double weight(int e, int s) const;  // s = +1 addition, -1 subtraction

    template <typename F>
    void forEachOut(int n, F &&f) const;
    template <typename F>
    void forEachIn(int n, F &&f) const;

    // Node-weighted Dijkstra from n0 (dir +1 follows arcs, -1 reverses
    // them); dist[n0] starts at weight(e0, s).
    void dijkstra(int n0, int s, int dir, std::vector<double> &dist,
                  std::vector<int> &parent) const;

    // Cheapest elementary circuit or source-to-sink path through edge e0,
    // as per-edge increments v in {0, 1, 2} with v[e0] >= 1. False when no
    // finite-weight strip through e0 exists.
    bool findGeneratingVector(int e0, int s, std::vector<int> &v) const;

    void stageI(Report &report);
    void stageII(Report &report);

    double objective() const;

    // Whether every requirement through an edge of v still holds.
    bool requirementsHold(const std::vector<int> &v) const;

    QuantTMesh &tm_;
    Options opts_;
    std::vector<std::vector<int>> requirementsOf_;  // edge -> requirement ids
    double eta_;                // #edges, the tier separation factor
    std::vector<int> sinks_;    // D-nodes whose unbalanced row is PHANTOM
};

#endif  // _TMESH_QUANTIZER_HXX_
