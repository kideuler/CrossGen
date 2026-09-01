#ifndef __DISK_TEMPLATE_HXX__
#define __DISK_TEMPLATE_HXX__

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "MERIDIAN/QuadMesh.hxx"
#include "mesh/Mesh.hxx"

// Circular inclusions, taken out of the layout problem and put back from a
// template.
//
// ### Why a disk is the one region worth templating
//
// Every other stage of this pipeline solves for a block structure because it
// does not know one. On a circular inclusion it does know one, exactly and in
// advance, and the knowing is not a heuristic: the O-grid -- a core block and
// a ring of blocks around it -- is what the topology of a disk leaves once you
// ask for quadrilaterals, and the cross field finds it unaided on every
// inclusion in the corpus. What does not work is asking Stage 6 to *rediscover*
// it ten times over. A disk carries eight cones (four +1 inside at about 0.6 of
// the radius, four -1 out in the matrix at about 1.4), so ten inclusions are
// eighty cones on top of the four the box already has, and the penalty
// continuation of Sec. 3.3 walls out somewhere around forty interior cones
// whatever the geometry. On data/meshes/multimat/bubbles it leaves 57% of the
// model unmeshed.
//
// So the disks are removed from the problem and the problem that remains is the
// matrix with ten holes. That is not a smaller version of the same difficulty,
// it is a different difficulty and a much smaller one -- measured on bubbles,
// the same pipeline that fails on the two-material model reaches Definition 2.1
// on the excised one with every face a quadrilateral and every element
// positive. The disks are then filled from the template, which costs no solve
// at all.
//
// ### Where the excision is legitimate, and where it is not
//
// Only a component `Mesh::computeMaterialCircles` accepted as a circle, and
// only one that is genuinely an *inclusion*: none of its boundary may be dS.
// The second test is not paperwork. A single-material disk model (geom003) is
// one component whose circle fit is perfect, and excising it would delete the
// model; a disk that touches the boundary is a different topological object
// whose rim is not a closed loop of the layout at all. Both are refused and
// reported.
//
// ### Parity, which is the whole of the difficulty
//
// The rim of an excised disk comes back from Stage 10 as a closed loop of N
// mesh edges, and the fill has to quadrangulate a disk with that boundary.
// Summing |boundary| over the faces of any quadrangulation of a disk gives
// 4F = 2E_interior + N, so **N is even or the fill does not exist** -- not
// "is hard", does not exist, for any template and any number of interior
// singularities. On bubbles at the default target edge length five of the ten
// rims come out odd, so this is the common case and not an edge case.
//
// N is not free but it is *reachable*: it is a sum of Stage 10's per-chord
// interval counts, and moving one chord by one edge changes the parity of every
// rim loop that chord crosses an odd number of times. That is a small linear
// system over GF(2) and it is solved inside QuadMesh, which is why
// QuadMesh::Options::evenLoops exists and why this stage has to hand it the rim
// loops *before* the intervals are chosen. See QuadMesh::fixLoopParity.
//
// The number of layout *arcs* on a rim has the same parity problem and no such
// remedy, which is the reason the fill is built at the level of elements rather
// than at the level of the arrangement: an arc-level fill would need one ring
// block per rim arc and a core quadrangulation with that many boundary arcs,
// and half the rims cannot have one. At the element level the template's own
// corners are free to sit at any rim node, so nothing outside the disk has to
// agree to anything except the number of edges.
//
// ### The template itself
//
//     N rim nodes, split into four sides of a, b, a, b edges
//     a ring of `c` rows of elements inward from the rim
//     a core of a x b elements inside it
//
// The four sides are chosen to be as near a quarter of the rim as the node
// spacing allows, which is a search over the N possible starting nodes and the
// N/2 - 1 possible splits and costs nothing. The ring depth comes from the rim
// spacing, so that a row of ring elements is about as deep as it is wide. The
// core boundary is a rounded square rather than a scaled copy of the rim --
// a scaled copy would put a 180 degree element corner at each of the four core
// corners, which is the same defect meshing the whole disk as one patch has and
// the reason the template exists.
//
// Everything is then Laplacian-smoothed with the rim held, node by node and
// only where the move does not lower the worst scaled Jacobian around it. On a
// convex region with this connectivity that is safe and it is what takes the
// transfinite initialisation to a mesh that looks like the textbook picture.
class DiskTemplate {
public:
    // One circular inclusion: a component of the input mesh, its fitted circle,
    // the triangles that go, and the rim they leave behind.
    struct Inclusion {
        int component = -1;           // index into Mesh::materialComponents
        int matId = 0;
        CircleFit circle;
        std::vector<int> triangles;   // in the full mesh
        std::vector<int> rim;         // rim vertices, CCW about the centre, open
        std::vector<int> rimExcised;  // the same, renumbered into the excised mesh
        double area = 0.0;            // on the model
    };

    struct Options {
        // Where the corners of the core block sit, in radii. Zero picks it from
        // the ring depth so that a ring element is about as deep as it is wide,
        // which is what makes the ring and the core agree in size; a fixed
        // value overrides that.
        double coreRadius = 0.0;
        // How square the core boundary is. 0 is a scaled copy of the rim, which
        // puts a straight angle at every core corner; 1 is the straight chords
        // between the four core corners, which makes the ring 1.6x deeper at
        // the middle of a side than at its ends. Between them is a rounded
        // square and is what the O-grid wants.
        double coreSquareness = 0.55;
        // Rows of elements in the ring. Zero picks it from the rim spacing.
        int ringDepth = 0;
        // Smoothing sweeps over the nodes the template placed, the rim held.
        int smoothingPasses = 300;
        double smoothingTolerance = 1e-4;
        // A rim shorter than this many edges has no room for a ring and a core
        // and is refused rather than meshed into something degenerate.
        int minRimEdges = 8;
    };

    // One block of the template, as a structured (ns+1) x (nt+1) grid of vertex
    // indices into the merged mesh, row-major. Five per inclusion: the core and
    // the four sides of the ring. Kept for the same reason QuadMesh::Block is:
    // the block structure is the useful thing downstream and cannot be
    // recovered from the elements.
    struct Block {
        int inclusion = -1;
        int ns = 0, nt = 0;
        std::vector<int> vert;
        bool core = false;
    };

    struct Report {
        int inclusions = 0;        // circular inclusions offered
        int filled = 0;            // ... and templated
        int refusedOdd = 0;        // rims with an odd number of edges
        int refusedShort = 0;      // rims with too few edges to hold a template
        int refusedOpen = 0;       // rims Stage 10 did not close
        int blocks = 0;
        int quads = 0;             // elements the template added
        int vertices = 0;          // ... and vertices, over the matrix mesh

        // The merged mesh: the matrix mesh Stage 10 produced with the templates
        // in it. Measured on the whole thing, because the point of the stage is
        // that the join is invisible.
        int mergedVertices = 0;
        int mergedQuads = 0;
        double minScaledJacobian = 0.0;
        double meanScaledJacobian = 0.0;
        double templateMinScaledJacobian = 0.0;  // over the template elements alone
        int invertedQuads = 0;
        double minEdge = 0.0, maxEdge = 0.0;
        int smoothingSweeps = 0;   // the most any one inclusion needed

        int interiorEdges = 0;
        int boundaryEdges = 0;
        int nonManifoldEdges = 0;
        int cracks = 0;
        bool conforming = false;
        bool valid = false;        // conforming, nothing inverted, nothing refused

        std::vector<std::string> messages;
    };

    // ---- Stage 0c, before anything else runs -----------------------------

    // The circular inclusions of a mesh. Empty on a single-material mesh, on a
    // mesh whose only circle is the model itself, and on any inclusion whose
    // rim touches dS -- each of those is reported through `messages`.
    static std::vector<Inclusion> detect(const Mesh &mesh,
                                         std::vector<std::string> &messages);

    // The mesh with those inclusions' triangles removed and its vertices
    // renumbered, so that every rim is a boundary loop of what comes back.
    // Fills Inclusion::rimExcised. Returns nullptr if nothing would be left.
    static std::shared_ptr<Mesh> excise(const Mesh &mesh,
                                        std::vector<Inclusion> &inclusions);

    // ---- finding the rims again, between Stages 8 and 11 -----------------

    // Which arcs of the arrangement lie on each rim, one list per inclusion.
    // A rim is dS of the excised mesh, so its arcs are ArcKind::Boundary, and
    // which rim one is on is asked of its geometry rather than of the run of
    // mesh vertices it carries: an arc between two nodes Stage 8 put on the
    // same mesh edge carries no whole vertex at all, and a test on those drops
    // it and leaves the rim in pieces -- measured, six of the ten rims on
    // bubbles. The rim polyline is inscribed in the fitted circle to within one
    // segment's sagitta and the nearest other piece of dS is many radii away,
    // so the geometric test is not close to ambiguous.
    //
    // Stage 10 needs these to make the rims carry an even number of edges
    // (QuadMesh::Options::evenLoops); Stage 11 needs them to find the rim in
    // the finished mesh.
    static std::vector<std::vector<int>> rimArcs(const Arrangement &arr,
                                                 const std::vector<Inclusion> &inclusions);

    // The same rims as closed loops of quadrilateral-mesh vertices, in the
    // order the loop runs and with no repeated first node. Empty for a rim
    // Stage 10 did not close.
    static std::vector<std::vector<int>> rimVertexLoops(
        const Arrangement &arr, const QuadMesh &qm,
        const std::vector<std::vector<int>> &arcsPerRim);

    // ---- Stage 11, after Stage 10 ----------------------------------------

    // `rimLoops[i]` is the closed loop of quadrilateral-mesh vertices around
    // inclusion `i`, in the order the loop runs, with no repeated first node.
    // An empty loop means Stage 10 did not close that rim and the inclusion is
    // reported and left as a hole.
    DiskTemplate(const std::vector<Point> &verts,
                 const std::vector<std::array<int, 4>> &cells,
                 const std::vector<int> &cellMaterial,
                 const std::vector<std::vector<int>> &rimLoops,
                 const std::vector<Inclusion> &inclusions,
                 const Options &opts);

    const std::vector<Point>& vertices() const { return verts; }
    const std::vector<std::array<int, 4>>& quads() const { return cells; }
    const std::vector<int>& quadMaterials() const { return cellMaterial; }
    const std::vector<Block>& blocks() const { return grids; }
    const Report& getReport() const { return report; }

    bool writeOBJ(const std::string &filename) const;
    bool writeVTU(const std::string &filename) const;

private:
    // One inclusion. Returns false, with a reason on the report, when the rim
    // cannot hold a template.
    bool fillOne(int inclusion, const std::vector<int> &loop, const CircleFit &cf,
                 int matId);
    // The four sides the rim is split into: the starting node and the two side
    // lengths, chosen so that the four corners are as near a quarter turn apart
    // as the node spacing allows.
    void chooseSides(const std::vector<int> &loop, const Point &centre,
                     int &start, int &a, int &b) const;
    // Laplacian smoothing over `movable`, the rim held, accepting a move only
    // where it does not lower the worst scaled Jacobian around the node.
    int smooth(const std::vector<int> &movable, int firstCell, double scale);
    void check();

    Options options;

    std::vector<Point> verts;
    std::vector<std::array<int, 4>> cells;
    std::vector<int> cellMaterial;
    std::vector<Block> grids;

    Report report;
};

#endif // __DISK_TEMPLATE_HXX__
