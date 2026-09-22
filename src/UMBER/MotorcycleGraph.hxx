#ifndef __MOTORCYCLEGRAPH_HXX__
#define __MOTORCYCLEGRAPH_HXX__

#include <string>
#include <unordered_set>
#include <vector>

#include "UMBER/Polysquare.hxx"
#include "mesh/Mesh.hxx"

// The meta-block structure of
//   Wang, Ren, Fang, Lin, Xu, Bao and Huang, "IGA-suitable planar
//   parameterization with patch structure simplification of closed-form
//   polysquare", CMAME 392 (2022) 114678, Section 5, first paragraph:
//   "we construct the meta-block of Omega via the motorcycle graph [60]. The
//   resulting patch structure is made up of the tracing iso-lines (iso-curves)
//   on Omega from the corners, which tessellate Omega into a set of
//   non-degenerate quad patches."
//
// Only that: the blocks, not the simplification of Sec. 5.1 onwards.
//
// A polysquare boundary is made of axis-aligned segments, so its image is a
// rectilinear domain, and cutting a rectilinear domain into rectangles is done
// by sending a ray inward from each reflex corner. That is what the
// motorcycles are. A convex corner needs nothing -- its two boundary segments
// are already iso-lines -- and a straight run of boundary needs nothing
// either, so the only sources are the corners the frame field marked with a
// negative index. A reflex corner (270 degrees in the image) has two axis
// directions pointing into the domain and a 360-degree one has three; each
// gets a motorcycle.
//
// The rays are traced through the mesh rather than in the plane, which is what
// makes the result independent of whether the parameterization's image happens
// to overlap itself: a ray is a straight line in the parameter domain, so
// inside one triangle it is a straight line in the model too (phi is affine
// there), and crossing into the next triangle only needs the direction carried
// over. Across a cut the two triangles disagree about which way the parameter
// axes point, and the rotation between them is read off the shared edge --
// its image has the same length on both banks and differs by exactly the
// transition, so no bookkeeping of Pi_gamma is needed here.
//
// The rays are not stopped where they meet each other, which is the one place
// this departs from the motorcycle graph of [60] as usually written. Every
// iso-line is drawn from its corner until it leaves the model, and two of them
// simply cross. The reason is the sentence Sec. 5 opens with: the meta-block
// structure "has neither T-junction nor internal singularity". A ray halted on
// another ray's trail ends in the middle of the domain, and that end is a
// T-junction; only lines carried through to the boundary avoid them. It also
// removes a dependence on arrival order that has nothing to do with the model:
// two rays leaving one corner can be handed the same triangle of its star when
// the mesh is coarse there -- they are 90 degrees apart but the wedge spans
// enough for both -- and they can only leave it through its one free edge, so
// whichever went second would die one step out. That is what cost
// data/meshes/singlemat/geom005 a trace and a block.
class MotorcycleGraph {
public:
    using EdgeKey = MeshEdgeKey;
    using EdgeKeyHash = MeshEdgeKeyHash;

    struct Report {
        int motorcycles = 0;    // rays launched, one per interior axis at a reflex corner
        // Places two iso-lines actually cross, found by intersecting the steps
        // that share a triangle rather than by counting edges twice marked --
        // which overcounts, since two rays out of one corner share the first
        // edge out of the triangle they start in.
        int crossings = 0;
        // Rays that left the model, whether through a boundary edge or by
        // landing exactly on a boundary vertex. The second is what an iso-line
        // running into the corner it was aimed at looks like, and in a
        // polysquare that is where a line out of a reflex corner is supposed to
        // finish -- see the note in run() where exitPoint() comes back empty.
        int reachedBoundary = 0;
        // Rays that neither crashed nor left the model. On a locally
        // injective parameterization this stays zero; a folded triangle
        // reverses the direction a ray reads there, so one can circle instead
        // of going anywhere, and the step cap is what ends it. A boundary loop
        // whose image does not close does the same thing without any triangle
        // being folded: the image is then an unbounded strip and a ray along it
        // has no wall to hit.
        int ranOut = 0;
        // Rays the corner index called for that the parameterization has no
        // room for. Zero unless Sec. 4.3 left something inconsistent there.
        int skippedCorners = 0;
        // Rays dropped for being another ray drawn backwards -- see the end of
        // run(). Two reflex corners that face each other each aim at the other,
        // so the iso-line between them gets traced from both ends, and one of
        // the two is the line.
        int duplicates = 0;
        // Rays ended by arriving at the corner another ray was launched from,
        // rather than by leaving the model -- see setArrivalTolerance().
        int arrivals = 0;
        // Rays launched to carry an iso-line through an interface it ended on,
        // one per landing. Zero on a single-material model. They are counted in
        // `motorcycles` as well, being rays like any other; this says how many
        // of them nobody's corner asked for.
        int continuations = 0;
        // Continuations refused for arriving at the generation cap. Nonzero
        // means a line crossed more interfaces than the cap allows, which on
        // these models means it was going round in circles.
        int continuationsDropped = 0;
        int nodes = 0;         // corners, exits and crossings of the block structure
        int blocks = 0;
        int mergedSlivers = 0;  // blocks the mesh could not resolve, folded into a neighbour
        int smallestBlock = 0;  // triangles in the smallest block
        int tracedEdges = 0;
    };

    // The polysquare must have been solved; its corners and its (u, v) are
    // what the tracing follows. On a multi-material model it also carries the
    // interfaces, and those change the structure in two ways -- see
    // launchFeatureSectors() and floodBlocks().
    explicit MotorcycleGraph(const Polysquare &ps);

    // How close a ray that is already lost has to pass to the corner another
    // ray was launched from for it to count as having arrived there: the trace
    // is snapped onto that corner and stops. In mean mesh edges; 0 turns it off.
    //
    // "Already lost" is doing as much work here as the distance is, and run()
    // says why: rays that are behaving perfectly also pass launch corners in
    // mid flight, closer than the ones that are lost do, so no tolerance
    // separates the two on its own. The length does. Measured over
    // data/meshes at the default, twenty-eight of the thirty models come out
    // bit for bit as they do with this switched off, and the two that change
    // are the two that had a ray going nowhere: geom031 from 7968 blocks to
    // 371, geom035 from 241 to 51.
    //
    // It is not free on those two. Snapping the end of a trace sideways by up
    // to the tolerance can push it across a neighbouring line, and it does:
    // geom035 picks up two places where sides cross with no node between them
    // and geom031 five. Against a trace that never terminated at all that is
    // the better of the two, which is why this is on -- but it is a rescue for
    // a parameterization that is already wrong, not a thing that improves a
    // right one.
    void setArrivalTolerance(double t) { arrivalTol = t; }

    // Launch the rays, then label the triangles.
    void build();

    // Block index per triangle of the *original* mesh, which the cut mesh
    // shares. -1 only if the mesh has a triangle no flood reached, which
    // cannot happen on a connected mesh.
    const std::vector<int>& getBlockOfTriangle() const { return blockOfTriangle; }
    int numBlocks() const { return report_.blocks; }

    // The rays, as polylines in model coordinates.
    const std::vector<std::vector<Point>>& getTraces() const { return traces; }

    // One step of one ray, carried in both domains. A segment lies inside a
    // single triangle, which is what makes its image exact: phi is affine
    // there, so the straight line in the parameter domain and the straight
    // line in the model are the same segment seen twice. Keeping the two ends
    // of a step together also keeps the seams right -- a point on a cut has an
    // image on each bank, and the one that belongs to a step is the one taken
    // in the triangle that step runs through.
    struct Segment {
        Point a{0.0, 0.0}, b{0.0, 0.0};    // model
        Point ua{0.0, 0.0}, ub{0.0, 0.0};  // parameter domain
        int tri = -1;
        int ray = -1;
        // Which step of that ray this is, so a point of the segment can be
        // named as a parameter along the ray's polyline: getTraces()[ray][step]
        // is a and [step + 1] is b.
        int step = -1;
    };
    const std::vector<Segment>& getSegments() const { return segments; }

    // A node of the block structure: where its edges meet.
    //
    // A node carries where it sits as well as where it is, so that the block
    // structure can be cut out of the traces without looking for it a second
    // time -- see BlockLayout, which needs each node's place along the lines
    // through it in order to split them into the sides of the blocks.
    //
    // A place along a ray is given as a parameter into its polyline: the index
    // of the point before it plus the fraction of the way along that step, so
    // that everything on one ray sorts on one number.
    struct Node {
        Point xy{0.0, 0.0};
        Point uv{0.0, 0.0};
        enum Kind { Corner, BoundaryEnd, Crossing } kind = Corner;

        // Corner: the mesh vertex it stands on.
        int vertex = -1;
        // Crossing: the two rays and the parameter of the crossing along each.
        // BoundaryEnd: ray[0] and its last parameter; ray[1] unused.
        int ray[2] = {-1, -1};
        double param[2] = {0.0, 0.0};
        // BoundaryEnd: the boundary edge the ray left through, and how far
        // along it, measured from edges[boundaryEdge][0] to [1].
        int boundaryEdge = -1;
        double alongEdge = 0.0;
    };
    // Every corner of the polysquare, every point an iso-line leaves the model
    // at, and every point two of them cross. Corners include the convex ones,
    // which launch no ray but are corners of a block all the same.
    const std::vector<Node>& getNodes() const { return nodes; }

    // Mesh edges a ray crossed. These are the walls the blocks are flooded
    // between, so a block boundary follows mesh edges even though the ray
    // itself cuts through triangle interiors -- see writeBlocksVTU().
    const std::unordered_set<EdgeKey, EdgeKeyHash>& getTracedEdges() const { return tracedEdges; }

    // Per ray, the boundary edge it left the model through and how far along
    // that edge, measured from edges[e][0] to edges[e][1]. -1 for a ray that
    // never got out -- one of Report::ranOut -- whose trace therefore ends in
    // the middle of the model and leaves the block structure open there.
    const std::vector<int>& getRayExitEdge() const { return exitEdge; }
    const std::vector<double>& getRayExitAlong() const { return exitAlong; }

    const Report& getReport() const { return report_; }

    // The model the blocks were cut out of.
    const Mesh& getMesh() const { return *mesh; }

    // The material interfaces the tracing treated as walls, as mesh edge
    // indices. Empty on a single-material model. The block structure has to
    // carry these as sides of its own -- see BlockLayout -- since a ray ends
    // on one and a block never crosses one.
    const std::vector<int>& getFeatureEdges() const { return poly->getFeatureEdges(); }
    bool isFeatureEdge(int e) const {
        return e >= 0 && e < static_cast<int>(edgeIsFeature.size()) && edgeIsFeature[e];
    }

    // The input mesh with the block index as cell data.
    //
    // A ray crosses triangles rather than following their edges, so a triangle
    // the ray passes through belongs to whichever side reached it first and
    // the block outlines are ragged at the scale of one triangle. The blocks
    // themselves -- how many, which parts of the model they cover, what they
    // are adjacent to -- do not depend on that.
    bool writeBlocksVTU(const std::string &filename) const;

    // The rays themselves, as polylines, to see where the block walls run
    // without the triangle-level raggedness.
    bool writeTracesVTU(const std::string &filename) const;

private:
    // A ray in flight: which triangle it is in, where, the parameter-domain
    // direction it follows, and how far it has come.
    struct Motorcycle {
        int tri = -1;
        int originVertex = -1;  // the corner it was launched from
        Point pos{0.0, 0.0};    // model coordinates
        Point dir{1.0, 0.0};    // parameter-domain direction, unit
        double distance = 0.0;
        int id = -1;
        int steps = 0;
        bool alive = true;
    };

    void launch();
    // The rays a multi-material domain needs on top of those.
    //
    // Sec. 5 sends a ray inward from every reflex corner, because cutting a
    // rectilinear domain into rectangles is done from its reflex corners and
    // nowhere else. On a multi-material domain the domain to be cut up is not
    // S but each material region, and the corners of those are the corners of
    // dS *and* every place an interface turns, meets another interface, or
    // lands on dS. A region whose reflex corners were not fired from is a
    // region that does not come out as quadrilaterals.
    //
    // The sectors are read off the image rather than off the frame field's
    // corner index, which is what the boundary path above uses. Both say the
    // same thing where the polysquare is sound -- the image angle of a sector
    // is q quarter turns and the index is 2 - q -- and reading the image is
    // what makes this independent of how the corner index was accumulated
    // along a chain, which on an interface is a per-side question with two
    // answers rather than dS's one.
    //
    // A vertex with no interface on it is left entirely to launch(), so a
    // single-material model takes this path nowhere and comes out bit for bit
    // as it did.
    void launchFeatureSectors();
    void run();
    void findNodes();
    void floodBlocks();

    // Where a point of triangle f lands in the parameter domain.
    Point imageOfPoint(int f, const Point &p) const;

    // grad phi of one triangle, and its inverse.
    void triangleJacobian(int f, double J[2][2]) const;
    // Where the ray leaves triangle f, and through which edge.
    bool exitPoint(int f, const Point &from, const Point &meshDir,
                   Point &hit, int &edge, double &along) const;

    const Polysquare *poly = nullptr;
    const Mesh *mesh = nullptr;      // the original mesh
    const Mesh *cut = nullptr;       // M_C, same triangles
    const std::vector<Point> *uv = nullptr;

    // The material interfaces, as a flag per mesh edge and per vertex. They
    // are walls for the flood and sector boundaries for the launch, and they
    // are empty on a single-material model.
    std::vector<char> edgeIsFeature;
    std::vector<char> vertexOnFeature;
    bool multiMaterial = false;

    std::vector<Motorcycle> bikes;
    std::vector<std::vector<Point>> traces;
    std::vector<Segment> segments;
    std::vector<Node> nodes;
    std::unordered_set<EdgeKey, EdgeKeyHash> tracedEdges;
    std::vector<int> blockOfTriangle;
    std::vector<char> onWall;   // a ray passed through this triangle
    std::vector<int> exitEdge;      // per ray, where it left the model
    std::vector<double> exitAlong;
    double arrivalTol = 0.20;

    Report report_;
};

#endif // __MOTORCYCLEGRAPH_HXX__
