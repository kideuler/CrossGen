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

// A boundary vertex counts as a sharp corner when its interior angle differs
// from a straight boundary by more than this -- convex below 135 degrees, or
// reflex above 225. Such a corner is forced to be a corner of the blocking:
// left in the interior of a block side, it is a kink no quad grid laid over
// that block can reproduce, which is the geometric unrealizability of Sec. 6.
// The threshold sits well clear of the mild reflex angles a discretized arc
// produces, so a circle is not mistaken for a ring of corners.
const double MEDIAL_SHARP_CORNER_TOLERANCE = M_PI / 4.0;

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
    int sharpCorners = 0;      // sharp boundary corners found
    int cornerCuts = 0;        // of those, the ones forced into block corners
    int cornersUnanchored = 0; // corners no incident medial vertex could carry
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
    // One step of the boundary walk: either a plain boundary vertex or the
    // footpoint of one medial vertex's passage.
    struct WalkEvent {
        Point pos{0.0, 0.0};
        int medial = -1;          // the vertex passed here, -1 for a plain step
        int boundaryVertex = -1;  // set on plain steps, for corner lookup
    };

    // A maximal group of consecutive events sharing one medial vertex. A run
    // longer than a point is a stretch of boundary collapsing onto that
    // vertex -- a polar section.
    struct MedialRun {
        int m = -1;
        int first = -1;
        int last = -1;
    };

    // One boundary loop, walked and grouped into runs. Where the cuts fall is
    // decided afterwards, so this is independent of the downsampling.
    struct LoopWalk {
        std::vector<WalkEvent> events;
        std::vector<MedialRun> runs;
        // {run, event} pairs where a sharp corner forces a cut.
        std::vector<std::array<int, 2>> cornerCuts;
    };

    // Walk every boundary loop into its event stream and runs. Runs depend
    // only on the map, not on which vertices survive downsampling, so this
    // comes first and both of the next two steps read it.
    void gatherLoops(const MedialAxisMap &map);

    // Pin every sharp boundary corner to a medial vertex that touches it, so
    // that the corner becomes a cut -- and hence a block corner -- rather
    // than a kink inside a block side. The vertex chosen is the nearest one
    // whose inscribed circle actually contacts the corner, which makes the
    // spoke from corner to vertex a genuine medial radius.
    void anchorCorners(const MedialAxisMap &map);

    // Junctions, endpoints and corner anchors always survive; the arc length
    // *between* consecutive survivors is then filled with samples every
    // ~targetSize, so a forced corner never leaves a sliver beside it. Cycles
    // keep at least three samples so a hole never degenerates to fewer than
    // two zones.
    void selectKeptVertices();

    // Cut each loop at its corner and kept-vertex cuts, and emit the zones
    // between them.
    void emitZones();

    // True when v's interior angle departs from straight by more than
    // MEDIAL_SHARP_CORNER_TOLERANCE.
    bool isSharpBoundaryCorner(int v) const;

    static void appendFanEvents(const MedialFan &fan, const Point &vpos,
                                int vIndex, std::vector<WalkEvent> &events);

    // Pair the zones across each chain into coarse edges, then color them.
    void buildEdges();
    void classify();

    // Apply the Fig. 17 template of each edge's class to its subdomain.
    void buildBlocks();

    std::vector<LoopWalk> loops_;
    // Medial vertices the downsampling is not allowed to drop: the corner
    // anchors. Junction and endpoint survival is handled per branch.
    std::vector<char> mandatory_;
    std::vector<char> kept_;
    // Footpoints of each passage of every medial vertex, filled by the walk.
    // A regular vertex is passed once per side; the two entries are its two
    // medial radii, whose angle at the vertex is the medial angle of Sec. 3.4.
    std::vector<std::vector<Point>> passages_;
    double targetSize_ = 0.0;
    TMeshStats stats_;
};

#endif // _MEDIAL_AXIS_TMESH_HXX_
