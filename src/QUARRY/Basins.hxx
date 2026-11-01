#ifndef __QUARRY_BASINS_HXX__
#define __QUARRY_BASINS_HXX__

#include <string>
#include <vector>

#include "mesh/Mesh.hxx"

namespace quarry {

// Stages 1 and 2 of QUARRY: the triangles of the model grouped into the regions one of SQuadGen's two chart
// distance fields describes -- the charts, read off the CDF and the offset to the chart centre (offsetc), or the
// dual charts, one per layout vertex, read off the DCDF and the offset to the subchart's corner (offset). The
// algorithm is the same for both, because the two fields are the same thing on the two layouts: a value that is 1
// at the region's centre and 0 on its walls, and an offset at every point that names the centre outright.
//
// ### What the release does, and why it over-segments
//
// QuadExtraction (cdfdcdf2quadmesh.cpp, docs/squadgen_2d.md Sec. 5.1) seeds a chart at every triangle whose CDF is
// the strict maximum of `ringsize` vertex rings, grows the seeds through the triangles above the wall threshold in
// order of the distance between predicted centres, and merges two clusters only when the boundary they share is
// mostly high-CDF. Every local maximum the decoder leaves in a chart's plateau is therefore a chart of its own unless
// the boundary test happens to catch it, and a decoded chart whose wall runs through its middle -- a dip of 0.3 the
// network put where the true field has none -- is always two. That is the failure the 2026-10-29 round trip found
// (memory squadgen-roundtrip-2d): true charts matched, extra clusters beside them.
//
// ### What is done here
//
// Two facts about the true fields make a sharper test possible. The CDF of every chart rises to exactly 1 at its
// centre whatever the chart's size or shape, since it is 1 - max(c0, c1) of coordinates normalised per chart; so a
// real chart's peak stands a full unit above the walls round it, where noise stands a few hundredths. And every
// triangle of a chart predicts the same centre, x + offsetc, which a spurious peak cannot change: a bump inside
// chart A, wherever it sits, still points at A's centre. Two peaks are one region when the centres they predict
// are much closer together than the peaks themselves are -- |q_a - q_b| < sameTarget |x_a - x_b| -- and two
// regions when they are not, which is a scale-free test: it needs no bandwidth, no ring count and no wall
// threshold, and it holds for a sliver as for a square.
//
//   Stage 1  Peaks. The persistence clustering of Chazal et al. (ToMATo, J. ACM 2013) on the triangles' edge
//            graph: triangles taken in order of decreasing value, each joining the region of its highest
//            neighbour, and two regions meeting at a saddle merged when the lower peak stands less than `noise`
//            above it (noise) or when their predicted centres agree (one region with a false wall). Each region
//            keeps the area-weighted mean of the predicted centres over its *core*, the triangles within `core` of
//            its peak value: the part of a chart where the decoded offsets are most reliable, far from the walls
//            the network blurs across. Where the region itself is, for the test, is its peak triangle, not its
//            core: the core of a dual chart on dS lies to one side of its vertex, and two vertices a few triangles
//            apart along dS then looked closer in their predictions than in their cores, and were merged. After
//            the sweep, adjacent regions are tested again with their final cores, and regions whose peak never
//            reached `minPeak` are dissolved: a real chart's peak is 1, a small one decoded perhaps 0.4.
//   Stage 2  Regions. The cores grown outwards over every triangle, wall triangles included, each triangle joining
//            the region it is reached from first in order of |x + offset - c|, the distance from the centre it
//            predicts to the region's centre c: the release's growth, but measured against each region's
//            consensus centre, not one seed triangle's own prediction. A region left in more than one piece keeps
//            the piece holding most of its core, and the rest is grown again from its neighbours.
//
// Nothing here knows whether the regions make a quad layout. That is the next stage's question (QUARRY.hxx).
class Basins {
public:
    struct Options {
        // A peak standing less than this above the saddle where it meets a higher one is noise, merged whatever its
        // offsets say. The true fields have no peak below 1 and no saddle above 0, so this only has to clear the
        // decoder's ripple. (Vetoing it where the two predict clearly different centres was tried, 2026-10-31: it
        // kept spurious peaks on blurred walls and saved no chart, since the charts lost to it are two small charts
        // the decoder blurred into one blob, whose offsets disagree by no more than their peaks do.)
        double noise = 0.1;
        // Two peaks are one region when the centres they predict lie closer together than this fraction of the
        // distance between the peaks themselves. Two real charts' peaks sit at their centres, so for them the
        // ratio is 1; a spurious peak predicts its chart's centre from elsewhere, so for it the ratio is near 0.
        double sameTarget = 0.5;
        // A region's core: its triangles within this of its peak value.
        double core = 0.25;
        // A region whose peak value stays below this is no chart and is dissolved into its neighbours. A real chart's
        // peak is 1; a small one's, decoded, can be 0.4.
        double minPeak = 0.25;
    };

    struct Region {
        int peak = -1;                 // triangle of the highest value
        double peakValue = 0.0;
        Point centre{0.0, 0.0};        // area-weighted mean predicted centre over the core
        Point position{0.0, 0.0};      // centroid of the peak triangle
        double coreArea = 0.0;
        double area = 0.0;             // of all its triangles after Stage 2
        int triangles = 0;
    };

    struct Report {
        int peaks = 0;                 // local maxima of the sweep, every one a region at first
        int noiseMerges = 0;           // merged at a saddle less than `noise` below the lower peak
        int targetMerges = 0;          // merged because both predicted the same centre (sweep and re-test)
        int dissolved = 0;             // peak below minPeak
        int regions = 0;
        int splitPieces = 0;           // pieces of a region regrown from their neighbours
        int unreached = 0;             // triangles no region reached (a mesh component with no region)
    };

    // `value` and `target` per triangle of `mesh`: the field and the point it predicts, centroid + offset.
    Basins(const Mesh &mesh, const std::vector<double> &value, const std::vector<Point> &target,
           const Options &opts);

    // Region of every triangle, -1 where none reached it.
    const std::vector<int> &labels() const { return labels_; }
    const std::vector<Region> &regions() const { return regions_; }
    const Report &getReport() const { return report_; }

private:
    void sweep();
    void retest();
    void grow(std::vector<int> &labels) const;
    void keepMainPieces();

    const Mesh &mesh_;
    const std::vector<double> &value_;
    const std::vector<Point> &target_;
    Options opts_;
    Report report_;
    std::vector<Point> centroid_;
    std::vector<double> area_;

    // Union-find over triangles during the sweep, and each root's running sums.
    std::vector<int> parent_;
    struct Sums {
        int peak = -1;
        double peakValue = 0.0;
        double w = 0.0;
        Point wq{0.0, 0.0};
    };
    std::vector<Sums> sums_;          // indexed by root triangle
    int find(int t);
    bool same(const Sums &a, const Sums &b) const;

    std::vector<int> roots_;          // the regions after Stage 1, as root triangles
    std::vector<int> labels_;
    std::vector<Region> regions_;
};

}  // namespace quarry

#endif // __QUARRY_BASINS_HXX__
