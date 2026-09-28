#ifndef __BOUNDARY_FEATURES_HXX__
#define __BOUNDARY_FEATURES_HXX__

#include <vector>

#include "mesh/Mesh.hxx"

// What a quad layout of a model owes to the model's own boundary, read off the
// triangle mesh before any method has run: the corners, the loops, and the
// irregularity the two of them force on every conforming layout.
//
// ### Why a quad layout needs this and the Shape-DNA does not give it
//
// The topology of a quad layout is decided by the model's corners and holes.
// docs/block_decomposition_metrics.md T3 has it as an identity: for a
// conforming layout of a region with Euler characteristic chi,
//
//     sum_interior (4 - d_v)  +  sum_boundary (3 - d_v)  =  4 chi,
//
// so once each corner is given its block count, the interior singularities
// are fixed up to pairs. A Dirichlet spectrum knows the area and the
// perimeter from its first two heat-trace terms, but the corners and the holes
// only from the third, which a hundred eigenvalues do not resolve. Measured on
// py/dataset.csv (2026-09-27), the thirteen numbers of Summary predicted which
// method to pick better than the hundred-eigenvalue Shape-DNA did, on shapes
// the model had seen and on ones it had not, and adding the spectrum to them
// made it worse.
//
// ### Regions, loops, corners
//
// The unit is the region, a connected run of one material (Mesh::
// materialComponents), because that is the unit T3 is stated for: a material
// interface acts as boundary, and the points where interfaces meet dS or each
// other act as corners, exactly as BlockDecomposition's sectors stop at an
// interface. Each region's boundary is walked as closed loops with the region
// on the left, so an outer loop runs counter-clockwise and a hole clockwise.
//
// A loop vertex is a corner where the loop turns by more than cornerAngle
// (20 degrees): G3's "every model corner sharper than about 160 degrees", and
// the same band BlockDecomposition::regularity() takes as a smooth point of
// dS. The corner rule of this class and the smooth band of that metric are the
// one number on purpose -- the two are compared in regularity().
//
// ### The fewest quarter turns (minimumDefect)
//
// regularity() counts W, the quarter turns by which a layout's block counts
// are off: |4 - d| inside, |2 theta/pi - n| at a boundary sector of angle
// theta split into n blocks. For a corner split into n_c blocks, the identity
// above leaves the interior at least |4 chi - sum_c (2 - n_c)| quarter turns
// out (smooth points of dS can take some of it, at the same price), so every
// conforming layout has
//
//     W  >=  min over n_c >= 1 of   sum_c |2 theta_c/pi - n_c|  +  |4 chi - sum_c (2 - n_c)|,
//
// per region. That minimum is minimumDefect(). It is a property of the model
// alone, and the quarter disk (1), the disk (4) and geom006 (4) meet it with
// the layouts every method finds. Only the two integers nearest 2 theta_c/pi
// need trying for each corner -- going further costs a quarter turn at the
// corner for at most one back in the interior -- so it is a sort, not a search.
class BoundaryFeatures {
public:
    struct Options {
        // Turning, in degrees, beyond which a boundary vertex is a corner.
        double cornerAngle = 20.0;
        // Turning, in degrees, beyond which a vertex that is not a corner
        // bends the boundary: a circle at the corpus mesh size turns about a
        // degree per vertex, a straight side of a polygon by rounding only.
        double curvedAngle = 0.5;
    };

    // A corner of one region's boundary.
    struct Corner {
        int vertex = -1;          // mesh vertex
        int region = -1;          // index into Mesh::materialComponents
        int loop = -1;            // index into loops()
        Point p{0.0, 0.0};
        double angle = 0.0;       // interior angle theta, radians, in (0, 2 pi)
        int blocks = 1;           // k = max(1, round(2 theta / pi)): T2/T4's ideal block count
        double step = 0.0;        // the shorter of the two boundary edges at the corner
    };

    // One closed boundary loop of one region, walked with the region on its left.
    struct Loop {
        int region = -1;
        std::vector<int> vertices;
        std::vector<double> turn;       // signed turning at each vertex, radians, + to the left
        std::vector<bool> onBoundary;   // edge i (vertices[i] -> [i + 1]) lies on dS, not an interface
        std::vector<int> corners;       // indices into corners(), in loop order
        double length = 0.0;
        double signedArea = 0.0;        // > 0 for a region's outer loop, < 0 for a hole
    };

    // A stretch of a loop from one corner to the next -- or the whole loop,
    // for one with no corner.
    struct Run {
        int loop = -1;
        double length = 0.0;
        double turning = 0.0;           // radians turned between its two corners, + to the left
    };

    // The numbers a learner reads: all unchanged by moving, rotating,
    // mirroring or scaling the model. Counts and defects are summed over the
    // regions, so a corner where an interface meets dS counts once for each
    // region it bounds, as it does in T3.
    struct Summary {
        int regions = 0;
        int holes = 0;                  // loops beyond the first of each region
        int euler = 0;                  // sum over regions of chi = 2 - loops
        int corners = 0;
        int cornersOneBlock = 0;        // k = 1: theta below 135 degrees
        int cornersTwoBlocks = 0;       // k = 2, but outside the smooth band
        int cornersThreeBlocks = 0;     // k = 3: theta in [225, 315)
        int cornersFourBlocks = 0;      // k >= 4
        int acuteCorners = 0;           // theta below 80 degrees
        int ambiguousCorners = 0;       // within 10 degrees of 135 or 225: T4's either-way corners
        double cornerDefect = 0.0;      // sum_c |2 theta_c/pi - k_c|: what the corners owe on their own
        double singularityBound = 0.0;  // T3's B = sum_regions |4 chi - sum_c (2 - k_c)|
        double minimumDefect = 0.0;     // the fewest quarter turns any conforming layout can have
        double isoperimetricRatio = 0.0;  // P^2 / (4 pi A) of dS: 1 for a disk
        double curvedFraction = 0.0;    // share of dS's length that bends
        double shortestRun = 0.0;       // shortest corner-to-corner stretch of any loop / sqrt(A)
        double interfaceLength = 0.0;   // total interface length / sqrt(A)
        double area = 0.0;              // of the model, in its own units
        double perimeter = 0.0;         // length of dS, in the model's units
    };

    explicit BoundaryFeatures(const Mesh &mesh);
    BoundaryFeatures(const Mesh &mesh, const Options &options);

    const std::vector<Loop> &loops() const { return loops_; }
    const std::vector<Corner> &corners() const { return corners_; }
    const std::vector<Run> &runs() const { return runs_; }
    const Summary &summary() const { return summary_; }
    const Options &options() const { return options_; }

    // The fewest quarter turns W any conforming quad layout of the model can
    // have (the header's derivation); what regularity() measures a layout
    // against.
    double minimumDefect() const { return summary_.minimumDefect; }

private:
    Options options_;
    std::vector<Loop> loops_;
    std::vector<Corner> corners_;
    std::vector<Run> runs_;
    Summary summary_;
};

#endif // __BOUNDARY_FEATURES_HXX__
