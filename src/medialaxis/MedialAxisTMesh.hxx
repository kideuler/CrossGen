#ifndef _MEDIAL_AXIS_TMESH_HXX_
#define _MEDIAL_AXIS_TMESH_HXX_

#include "medialaxis/MedialAxisMap.hxx"

#include <array>
#include <memory>
#include <vector>

// ─── Coarse block decomposition from the medial axis ────────────────────────
//
// Section 4 of "Planar quad meshing with guarantees", up to but not including
// the per-block quad templates: the axis is downsampled to a target element
// size and the map phi cuts the domain into quad-shaped zones, one on each
// side of every coarse medial edge (Fig. 3b/3c, Fig. 1-left). Each coarse
// edge is then classified by its medial angle, which decides the template a
// later meshing stage would apply. The output is deliberately a *block*
// decomposition -- a handful of four-sided zones -- not a refined mesh.
//
// Construction rides entirely on the monotonicity of phi: walking the
// boundary forward moves phi(p) monotonically along the axis, visiting every
// medial vertex once per side. Cutting the walk at the footpoint of every
// *kept* (downsampled) vertex therefore splits the boundary into intervals
// that pair with chains of the axis -- the zones -- with no explicit branch
// bookkeeping. A run of boundary collapsing onto a single kept vertex (the
// discrete axis endpoint of Sec. 3.3) becomes a cap zone instead.

// The four classes of Sec. 4.2. Green/red by the medial angle against 3*pi/4,
// blue/purple for coarse edges incident to an axis endpoint, split by whether
// the endpoint's dual cell is a triangle (a sharp tip) or a merged polygon (a
// round cap).
enum class MedialColor { Green = 0, Red = 1, Blue = 2, Purple = 3 };

const double MEDIAL_COLOR_ANGLE_THRESHOLD = 3.0 * M_PI / 4.0;

// One four-sided zone. Its sides are: the boundary interval, the spoke from
// each end of that interval to the matching end of the medial chain, and the
// chain itself. A cap zone degenerates the chain to a single vertex, so both
// spokes meet there and the zone is the polar region around an endpoint.
struct MedialZone {
    // Boundary side, in walk order. front() and back() are the corners the
    // spokes leave from; the points between them keep the true boundary
    // geometry.
    std::vector<Point> boundarySide;

    // Medial side, oriented with the walk: boundarySide.front() connects to
    // chain.front() and boundarySide.back() to chain.back(). Indices into
    // MedialAxis::medialVertices; a single entry for a cap.
    std::vector<int> chain;

    int coarseEdge = -1;   // owning edge, -1 for caps
    bool cap = false;
    int capVertex = -1;    // the endpoint a cap wraps, -1 otherwise

    MedialColor color = MedialColor::Green;
};

// One edge of the downsampled axis, shared by the two zones flanking it.
// Their union is the "subdomain associated with the medial edge" of Sec. 4.1.
struct CoarseEdge {
    std::vector<int> chain;        // original medial vertices, kept at the ends
    std::array<int, 2> zones{-1, -1};
    MedialColor color = MedialColor::Green;

    // Medial angle statistics along the chain. The mean decides green vs red:
    // a coarse chain absorbs many vertices and its angle legitimately narrows
    // near junctions, so the minimum would flag nearly every edge red. The
    // minimum is still recorded for a later, finer template pass.
    double meanAngle = M_PI;
    double minAngle = M_PI;
};

// One block of the coarse T-mesh: the subdomain of a coarse edge remeshed by
// the quad template of its class (Fig. 17). `outline` is the closed polygon
// of the block, curved sides (boundary runs, axis pieces) at full geometry;
// `corners` are the block's logical corner points, a subset of the outline.
struct TMeshBlock {
    std::vector<Point> outline;
    std::vector<Point> corners;
    MedialColor color = MedialColor::Green;
};

struct TMeshStats {
    int loops = 0;             // boundary loops walked
    int skippedLoops = 0;      // loops with broken links or no kept cut
    int keptVertices = 0;      // medial vertices surviving the downsampling
    int zones = 0;
    int capZones = 0;
    int coarseEdges = 0;
    int unpairedEdges = 0;     // edges that did not find two flanking zones
    int blocks = 0;
    std::array<int, 4> colorCounts{0, 0, 0, 0};  // indexed by MedialColor
};

class MedialAxisTMesh {
public:
    // `targetSize` is the desired coarse edge length in model units; anything
    // <= 0 falls back to a tenth of the bounding-box diagonal. The map is
    // only read during construction, so it need not outlive this object; the
    // axis is shared and kept.
    explicit MedialAxisTMesh(const MedialAxisMap &map, double targetSize = 0.0);

    std::shared_ptr<MedialAxis> axis;

    std::vector<MedialZone> zones;
    std::vector<CoarseEdge> edges;

    // The coarse T-mesh: each edge's subdomain remeshed by its template.
    // Green yields 2 blocks, red and purple 3, and blue 1 (Fig. 17); caps
    // not consumed by a blue or purple template are kept as one block each.
    std::vector<TMeshBlock> blocks;

    const TMeshStats &stats() const { return stats_; }
    double targetSize() const { return targetSize_; }

    bool isKept(int medialVertex) const {
        return medialVertex >= 0 &&
               medialVertex < static_cast<int>(kept_.size()) &&
               kept_[medialVertex] != 0;
    }

private:
    // Junctions and endpoints always survive; branches keep interior samples
    // every ~targetSize of arc length, and cycles keep at least three so a
    // hole never degenerates to fewer than two zones.
    void selectKeptVertices();

    // Walk one boundary loop, cut it at the kept footpoints, and emit zones.
    // `startVertex` is any vertex of the loop; `visited` is shared across
    // loops so each is walked once.
    void walkLoop(const MedialAxisMap &map, int startVertex,
                  std::vector<char> &visited);

    // Pair the zones across each chain into coarse edges, then color them.
    void buildEdges();
    void classify();

    // Apply the Fig. 17 template of each edge's class to its subdomain.
    void buildBlocks();

    std::vector<char> kept_;
    // Footpoints of each passage of every medial vertex, filled by the walk.
    // A regular vertex is passed once per side; the two entries are its two
    // medial radii, whose angle at the vertex is the medial angle of Sec. 3.4.
    std::vector<std::vector<Point>> passages_;
    double targetSize_ = 0.0;
    TMeshStats stats_;
};

#endif // _MEDIAL_AXIS_TMESH_HXX_
