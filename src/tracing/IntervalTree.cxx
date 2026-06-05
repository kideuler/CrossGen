#include "IntervalTree.hxx"

#include <cassert>
#include <numeric>
#include <stdexcept>

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

IntervalTree::IntervalTree(const UVGParam& param, Axis axis)
{
    const Eigen::VectorXd& coordVec = (axis == Axis::U) ? param.getU() : param.getV();
    const Mesh& mesh = param.getCutMesh().getCutMesh();
    const int nT = static_cast<int>(mesh.triangles.size());
    const int nV = static_cast<int>(coordVec.size());

    if (nT == 0) return;

    // Build one Interval per triangle.
    std::vector<Interval> ivs;
    ivs.reserve(nT);

    for (int t = 0; t < nT; ++t) {
        const Triangle& tri = mesh.triangles[t];
        const int i = tri[0], j = tri[1], k = tri[2];
        if (i < 0 || i >= nV || j < 0 || j >= nV || k < 0 || k >= nV) continue;

        double a = coordVec(i), b = coordVec(j), c = coordVec(k);
        Interval iv;
        iv.lo     = std::min({a, b, c});
        iv.hi     = std::max({a, b, c});
        iv.triIdx = t;
        ivs.push_back(iv);
    }

    // Sort by lo to enable median picking.
    std::sort(ivs.begin(), ivs.end(),
              [](const Interval& x, const Interval& y){ return x.lo < y.lo; });

    // Recursively build the tree; index 0 will be the root.
    nodes_.reserve(ivs.size() * 2);
    build(ivs, 0, static_cast<int>(ivs.size()));
}

// ---------------------------------------------------------------------------
// build – recursive helper
// ---------------------------------------------------------------------------

int IntervalTree::build(std::vector<Interval>& ivs, int lo, int hi)
{
    if (lo >= hi) return -1;

    // Choose center as the median of the midpoints of the intervals in [lo,hi).
    // Using the median of lo-values (already sorted) as a proxy is equivalent.
    const int mid = (lo + hi) / 2;
    const double center = 0.5 * (ivs[mid].lo + ivs[mid].hi);

    Node node;
    node.center = center;

    // Partition into left-only, straddling, right-only.
    std::vector<Interval> leftIvs, rightIvs;

    for (int i = lo; i < hi; ++i) {
        const Interval& iv = ivs[i];
        if (iv.hi < center) {
            leftIvs.push_back(iv);
        } else if (iv.lo > center) {
            rightIvs.push_back(iv);
        } else {
            // Straddles center.
            node.byLo.push_back(iv);
            node.byHi.push_back(iv);
        }
    }

    // Sort straddle lists.
    std::sort(node.byLo.begin(), node.byLo.end(),
              [](const Interval& a, const Interval& b){ return a.lo < b.lo; });
    std::sort(node.byHi.begin(), node.byHi.end(),
              [](const Interval& a, const Interval& b){ return a.hi > b.hi; });

    // Allocate this node.
    const int nodeIdx = static_cast<int>(nodes_.size());
    nodes_.push_back(std::move(node));

    // Recurse (must rebuild left/right sub-ranges from our local vectors).
    // Re-use the ivs buffer by sorting leftIvs and rightIvs by lo.
    std::sort(leftIvs.begin(),  leftIvs.end(),
              [](const Interval& a, const Interval& b){ return a.lo < b.lo; });
    std::sort(rightIvs.begin(), rightIvs.end(),
              [](const Interval& a, const Interval& b){ return a.lo < b.lo; });

    if (!leftIvs.empty()) {
        nodes_[nodeIdx].left  = build(leftIvs,  0, static_cast<int>(leftIvs.size()));
    }
    if (!rightIvs.empty()) {
        nodes_[nodeIdx].right = build(rightIvs, 0, static_cast<int>(rightIvs.size()));
    }

    return nodeIdx;
}

// ---------------------------------------------------------------------------
// query
// ---------------------------------------------------------------------------

std::vector<int> IntervalTree::query(double value) const
{
    std::vector<int> result;
    if (nodes_.empty()) return result;
    query(0, value, result);
    return result;
}

void IntervalTree::query(int nodeIdx, double value, std::vector<int>& out) const
{
    if (nodeIdx < 0 || nodeIdx >= static_cast<int>(nodes_.size())) return;

    const Node& node = nodes_[nodeIdx];

    if (value < node.center) {
        // Query is to the left of center.
        // Report all straddling intervals with lo <= value (byLo is ascending).
        for (const Interval& iv : node.byLo) {
            if (iv.lo > value) break;
            out.push_back(iv.triIdx);
        }
        query(node.left, value, out);
    } else {
        // Query is at or to the right of center.
        // Report all straddling intervals with hi >= value (byHi is descending).
        for (const Interval& iv : node.byHi) {
            if (iv.hi < value) break;
            out.push_back(iv.triIdx);
        }
        query(node.right, value, out);
    }
}