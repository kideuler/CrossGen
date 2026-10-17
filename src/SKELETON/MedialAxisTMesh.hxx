#ifndef __MEDIAL_AXIS_TMESH_HXX__
#define __MEDIAL_AXIS_TMESH_HXX__

#include "SKELETON/MedialAxisMap.hxx"

#include <array>
#include <memory>
#include <unordered_map>
#include <vector>

// ─── Coarse block decomposition from the medial axis ────────────────────────
//
// Steps 6-8 of docs/medial_axis.md: Sec. 4 of "Planar quad meshing with
// guarantees". The axis is downsampled to the target size h branch by branch
// (Step 6.1), the map phi cuts the domain into one quad on each side of
// every coarse medial edge (6.2, Fig. 3b), each coarse edge is coloured by
// its medial angle or by the endpoint it ends in (Step 7), and the colour
// picks the Fig. 17 template its subdomain is remeshed with (Step 8). The
// output is the T-mesh of Fig. 3c: four-sided cells with T-junctions on the
// radii, which Step 9 quantizes.
//
// Construction rides entirely on the monotonicity of phi: walking the
// boundary forward moves phi(p) monotonically along the axis, visiting every
// medial vertex once per contact (Step 5.2). Cutting the walk at the
// preimage of every *kept* vertex therefore splits the boundary into the
// arcs b1 -> b2 of the quads, with no explicit branch bookkeeping. The walk
// is laid on the original boundary polyline: each preimage is a parameter on
// a current boundary edge, and the edge stands for its stretch of input
// boundary (Step 1.6, 8.3), so the blocks follow the input boundary even
// where Step 4 moved or removed vertices.
//
// A contact that is a polar section (Step 5.2 (c)) is one stretch of
// boundary collapsing onto one medial vertex. At an axis endpoint it bounds
// the polar region the endpoint's template swallows (Step 6.3). Anywhere
// else its arclength midpoint stands in as the single preimage for both
// subdomains beside it (Step 5.3, [R]).

// The four classes of Sec. 4.2. Green/red by the medial angle against 3*pi/4,
// blue/purple for coarse edges incident to an axis endpoint: blue where the
// endpoint is a sharp tip -- its cell is a triangle, or its polar section
// holds a locked convex corner ([R]) -- purple at a round cap.
enum class MedialColor { Green = 0, Red = 1, Blue = 2, Purple = 3 };

const double MEDIAL_COLOR_ANGLE_THRESHOLD = 3.0 * M_PI / 4.0;

// One side of one coarse edge before the templates: the quad (m1, m2, b2, b1)
// of Step 6.2, or for a cap the polar region of an endpoint (Step 6.3).
struct MedialZone {
    // Boundary side, in walk order, on the original boundary. front() and
    // back() are the preimages the radii run to; the points between keep the
    // true boundary geometry.
    std::vector<Point> boundarySide;
    // Parallel to boundarySide: the boundary vertex that point is, or -1 for
    // a preimage or a point inside an edge. A blue cap reads its corner p
    // from here.
    std::vector<int> boundaryVertexAt;

    // Medial side, oriented with the walk: boundarySide.front() is the
    // preimage of chain.front() and boundarySide.back() that of
    // chain.back(). Node ids (see MedialAxisTMesh::nodeCoord); every fine
    // vertex passed in between is listed, but only the two ends are coarse
    // vertices. A single entry for a cap. Walking with the domain on the
    // left, the zone lies to the right of chain.front() -> chain.back().
    std::vector<int> chain;

    // The radii the zone is bounded by at chain.front() and chain.back(),
    // as ids into MedialAxisTMesh's spoke list: two zones meeting at a
    // preimage share the spoke there. A template point on a radius (Step 8
    // Assembly 1) is keyed by it.
    std::array<int, 2> spokes{-1, -1};

    int coarseEdge = -1;   // owning edge, -1 for caps
    bool cap = false;
    int capVertex = -1;    // the endpoint a cap wraps, -1 otherwise
    // For a blue cap, the index into boundarySide of the convex corner p
    // the polar region is split at (Fig. 17); -1 otherwise.
    int capCorner = -1;

    MedialColor color = MedialColor::Green;
};

// One edge (m1, m2) of the downsampled axis, shared by the two zones
// flanking it. Their union is the subdomain of the edge (Sec. 4.1).
struct CoarseEdge {
    std::vector<int> chain;        // node ids, m1 = front(), m2 = back()
    std::array<int, 2> zones{-1, -1};
    int branch = -1;               // the medial branch it was sampled from
    MedialColor color = MedialColor::Green;

    // theta_{m,E} at m1 and m2 (Step 5.4): the angle between the radii of
    // the two contacts beside the edge, i.e. between the spokes of its two
    // zones. theta_E, the mean of the two ([R]), splits green from red.
    // Endpoint ends measure nothing and keep pi.
    std::array<double, 2> endAngle{M_PI, M_PI};
    double meanAngle = M_PI;
};

// One cell of the T-mesh. `outline` is its closed polygon, CCW, boundary
// sides at full geometry; `corners` are its four logical corners, a subset
// of the outline. A leftover zone that no template claimed keeps whatever
// corners it has.
struct TMeshBlock {
    std::vector<Point> outline;
    std::vector<Point> corners;
    MedialColor color = MedialColor::Green;
    int edge = -1;   // the coarse edge whose template made it, -1 for a leftover
};

struct TMeshStats {
    int loops = 0;             // boundary loops walked
    int skippedLoops = 0;      // loops with broken links or no kept cut
    int keptVertices = 0;      // axis vertices surviving the downsampling
    int insertedVertices = 0;  // samples Step 6.1 had to add on fine edges
    int midpointCuts = 0;      // polar sections at non-endpoints cut at their midpoint
    int zones = 0;
    int capZones = 0;
    int coarseEdges = 0;
    int unpairedEdges = 0;     // edges that did not find two flanking zones
    int blocks = 0;
    int blueByCorner = 0;      // blue only by the locked-convex-corner rule
    int dissolvedVertices = 0; // interior T-mesh vertices of valence 2 removed (Step 8.5)
    // Coarse edges with both ends at one point -- the axis of a disc, which
    // is a point: every sample on one circle. No template applies; their caps
    // stay whole.
    int pointEdges = 0;
    int badFaces = 0;          // faces failing the Step 8.6 check after refinement
    int refinements = 0;       // rounds of local h decrease (Step 8.6)
    std::array<int, 4> colorCounts{0, 0, 0, 0};  // coarse edges by MedialColor
};

class MedialAxisTMesh {
public:
    struct Options {
        // h, the desired coarse edge length in model units; <= 0 falls back
        // to the axis's defaultTargetSize(), a tenth of the bounding-box
        // diagonal.
        double targetSize = 0.0;
        // Step 8.6: a face with no area or a corner above this interior angle
        // is not a usable cell. [R]: the guide says "severely non-convex"
        // without a number.
        double maxCornerAngle = M_PI + M_PI / 9.0;
        // Step 8.6: halve h on the branches that produced bad faces and
        // rebuild, at most this many times.
        int maxRefinements = 4;
    };

    // The map is only read during construction, so it need not outlive this
    // object; the axis is shared and kept.
    explicit MedialAxisTMesh(const MedialAxisMap &map, double targetSize = 0.0);
    MedialAxisTMesh(const MedialAxisMap &map, const Options &opts);

    std::shared_ptr<MedialAxis> axis;

    std::vector<MedialZone> zones;
    std::vector<CoarseEdge> edges;

    // The coarse T-mesh: each edge's subdomain remeshed by its template.
    // Green yields 2 cells, red and purple 3, and blue 4 (Fig. 17).
    std::vector<TMeshBlock> blocks;

    const TMeshStats &stats() const { return stats_; }
    double targetSize() const { return targetSize_; }

    // Nodes are the axis's medial vertices, ids [0, medialVertices.size()),
    // followed by the samples Step 6.1 inserted on fine medial edges.
    int nodeCount() const;
    Point nodeCoord(int node) const;
    double nodeRadius(int node) const;
    int nodeDegree(int node) const;
    bool isKept(int node) const;

    // How many contacts (Step 5.2) the walk found for a node: one run of
    // passages per contact. Equal to the node's degree on a faithful axis.
    int contactCount(int node) const;

    // Every medial radius the zones are cut along, [m, phi^-1(m)]: one per
    // preimage, two at each end of a cap. Step 9 keeps the T-mesh edges
    // along each from collapsing all together (Prop. 1).
    std::vector<std::array<Point, 2>> radii() const;

private:
    // One step of the boundary walk: a boundary vertex, a vertex removed by
    // Step 4.2 (kept for geometry), or the preimage of one node's passage.
    struct WalkEvent {
        Point pos{0.0, 0.0};
        int node = -1;            // the node passed here, -1 for a plain step
        int boundaryVertex = -1;  // set on boundary-vertex steps
    };

    // A maximal group of consecutive passages of one node. A run longer than
    // a point is a polar section.
    struct MedialRun {
        int m = -1;
        int first = -1;
        int last = -1;
    };

    struct LoopWalk {
        std::vector<WalkEvent> events;
        std::vector<MedialRun> runs;
    };

    // A medial radius [m, b]: a node and the preimage of one of its contacts.
    struct Spoke {
        int node = -1;
        Point foot{0.0, 0.0};
    };

    // A sample inserted on a fine medial edge (Step 6.1).
    struct InsertedNode {
        Point coord{0.0, 0.0};
        double radius = 0.0;
    };

    void build(const MedialAxisMap &map);
    void selectNodes();
    void gatherLoops(const MedialAxisMap &map);
    void emitZones();
    void buildEdges();
    void classify();
    void buildBlocks();
    void dissolveValenceTwo(const std::vector<int> &diamondOf);
    // Step 8.6; returns the branches whose faces failed.
    std::vector<int> checkGeometry();

    std::vector<int> fanSequenceNodes(const MedialAxisMap &map, int v,
                                      std::vector<Point> &feet) const;

    double originalTurn(int v) const;

    Options opts_;
    double targetSize_ = 0.0;
    std::vector<double> branchScale_;          // Step 8.6 local h factor per branch
    std::vector<char> kept_;                    // per axis vertex
    std::vector<InsertedNode> inserted_;
    // Inserted node ids by the complex edge their medial edge is dual to.
    std::unordered_map<MeshEdgeKey, std::vector<int>, MeshEdgeKeyHash> insertedOn_;
    // Consecutive coarse vertices of every branch, as an unordered pair key,
    // -> the branch.
    std::unordered_map<long long, int> coarseBranch_;
    std::vector<int> runCount_;
    std::vector<int> originalPrev_;
    std::vector<LoopWalk> loops_;
    std::vector<Spoke> spokes_;
    TMeshStats stats_;
};

#endif // __MEDIAL_AXIS_TMESH_HXX__
