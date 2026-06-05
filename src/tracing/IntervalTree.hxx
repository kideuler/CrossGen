#ifndef __INTERVAL_TREE_HXX__
#define __INTERVAL_TREE_HXX__

#include <algorithm>
#include <vector>

#include <Eigen/Dense>

#include "Parameterization/UVGParam.hxx"

// ---------------------------------------------------------------------------
// IntervalTree
//
// A centered interval tree built over triangle UV (or VU) ranges.
//
// Construction
// ------------
//   IntervalTree tree(uvParam, IntervalTree::Axis::U);
//
//   For every triangle t the "interval" is
//       [min(u0,u1,u2), max(u0,u1,u2)]
//   where u_i = u_(v_i) from the UVGParam solution.
//
// Query
// -----
//   std::vector<int> hits = tree.query(0.3);   // all triangle indices whose
//                                               // u-range contains 0.3
//
// The implementation is a standard centered interval tree (de Berg et al.)
// with O(n log n) build time and O(log n + k) query time (k = output size).
// ---------------------------------------------------------------------------

class IntervalTree {
public:
    enum class Axis { U, V };

    // Build the tree from a solved UVGParam.
    // axis selects whether to index on U or V coordinates.
    IntervalTree(const UVGParam& param, Axis axis);

    // Return all triangle indices whose [lo, hi] interval contains 'value'.
    std::vector<int> query(double value) const;

private:
    // One entry per triangle.
    struct Interval {
        double lo, hi;
        int triIdx;
    };

    // A node of the centered interval tree.
    struct Node {
        double center;

        // Intervals that straddle 'center', sorted by lo (ascending) for
        // left-side queries and by hi (descending) for right-side queries.
        std::vector<Interval> byLo;  // sorted ascending on lo
        std::vector<Interval> byHi;  // sorted descending on hi

        int left  = -1;  // index into nodes_
        int right = -1;
    };

    std::vector<Node> nodes_;

    // Recursive build; returns node index in nodes_.
    int build(std::vector<Interval>& ivs, int lo, int hi);

    // Recursive query; appends triangle indices to 'out'.
    void query(int nodeIdx, double value, std::vector<int>& out) const;
};

#endif // __INTERVAL_TREE_HXX__