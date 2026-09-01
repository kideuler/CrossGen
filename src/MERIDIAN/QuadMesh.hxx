#ifndef __QUAD_MESH_HXX__
#define __QUAD_MESH_HXX__

#include <array>
#include <string>
#include <vector>

#include "MERIDIAN/SplineFit.hxx"
#include "mesh/Mesh.hxx"

// Stage 10 of the pipeline: the quadrilateral mesh on the layout Stages 8 and 9
// produced. (docs/shepherd2022.pdf Sec. 5 ends at the patches; the meshing of
// them is the "refinement and analysis handoff" it hands on.)
//
// Stage 9 left a set of bicubic patches, each a map from [0,1]^2 to S, glued
// along shared cubic arcs. Meshing that is not a geometry problem -- the
// geometry is already exact, and any point of a patch is one call to
// SplineFit::evaluate. It is an *integer* problem, and it has exactly one
// constraint:
//
//     a patch is meshed as an n_s x n_t grid, so the two sides that face each
//     other across it must be cut into the same number of edges.
//
// Nothing else couples one patch to another. Everything downstream of that
// constraint -- watertightness, conformity, the absence of hanging nodes --
// follows from meshing each arc once and handing the resulting points to both
// patches that share it, which is the same rule, for the same reason, that
// Stage 9 fits each arc once.
//
// ### Chords
//
// "Opposite sides carry the same count" is an equivalence relation on the arcs:
// arc a ~ arc b if they are opposite sides of some patch, closed transitively.
// A class of that relation is a **chord** of the layout -- the sequence of
// patches you sweep through by entering one side and leaving the far one, and
// the set of arcs that sweep cuts across. Assigning intervals is therefore not
// a per-arc decision and cannot be made to be one: every arc of a chord takes
// the same count or the layout does not mesh, and a chord can run the length of
// the model through dozens of patches of quite different sizes.
//
// The classes are found with a union-find over the arcs, which is O(arcs) and
// needs no traversal: for each quadrilateral patch, union side 0 with side 2
// and side 1 with side 3. A chord that closes on itself, and a chord that
// enters the same patch twice in perpendicular directions -- both of which
// happen on the corpus -- come out as the ordinary case, one class, without
// anything special being done about them. That is the reason to build them this
// way rather than by walking.
//
// ### Choosing the count
//
// One integer N per chord, and the patches it runs through generally disagree
// about what it should be. The count is chosen to minimise
//
//     F(N) = sum_k ( log( S_k / (N h) ) )^2
//
// against the target edge length h. The measure is a *ratio*: a row of elements
// twice the target length is as wrong as one half of it, which is what an
// aspect-ratio-driven mesher wants and what a difference of lengths does not
// say. Minimising over the reals gives
//
//     N* = geometricMean(S_k) / h
//
// and the integer answer is whichever of floor(N*) and ceil(N*) has the smaller
// F -- not round(N*), which is the minimiser of the wrong objective and is a
// different integer whenever N* is small, exactly where it matters.
//
// ### What S_k is, and why it is not the arc length
//
// The obvious data are the lengths of the chord's own arcs, and they are the
// wrong ones. An arc is the *edge* of a patch, and N does not only cut the
// edge: it cuts every row of elements across the patch, all the way to the
// opposite side. On a patch that is nearly a rectangle those are the same
// number and the distinction does not arise. On a patch that is not, they are
// not close -- geom024's layout is a single lens-shaped face whose two ends
// taper to 0.026 across and whose middle is 0.83, and judged on its end arcs
// that whole direction takes one element, which then spans the middle in one
// step sixteen times the target.
//
// So S_k is the **mean length of the isoparametric lines across the patch in
// that direction**, one number per patch the chord passes through, read off a
// sampling of the Stage 9 surface. It agrees with the arc length wherever the
// arc length was a fair summary of what the elements have to span, and reports
// what they actually have to span where it was not.
//
// A chord is never given fewer than Options::minIntervals edges. One is enough
// for the mesh to exist; two is what a patch needs before its interior has a
// node at all, and the difference shows up on the models whose layouts have a
// few very short arcs holding a corner together.
//
// ### Where the points go
//
// Along an arc, at equal arc length of the *fitted spline* -- not equal
// parameter, which on a chord-length fit is close but not the same, and not
// along the traced polyline, which is the data the spline was fitted to and is
// rougher than it. Inside a patch, the parameters are blended from the four
// sides by transfinite interpolation and the patch is then evaluated there, so
// interior nodes lie on the reconstructed surface exactly rather than on a
// bilinear guess at it. Boundary nodes are never computed twice: they are
// looked up from the arc, so the two patches sharing it get the same vertex
// index, and conformity is a property of the construction rather than a
// tolerance to be met.
//
// A block that comes out of that with an element turned over is then smoothed,
// and the smoothing is not cosmetic. A transfinite grid inherits whatever the
// four sides do: where one side is far more curved than the side facing it the
// blend crowds its rows together, and on the worst of the corpus it pushes them
// through each other -- geom003's layout has one face that a Coons patch cannot
// cover without folding, and 42 of its 1260 elements come out inverted.
// Winslow's elliptic system,
//
//     alpha x_ss - 2 beta x_st + gamma x_tt = 0,
//     alpha = |x_t|^2,  beta = x_s . x_t,  gamma = |x_s|^2,
//
// solved by Gauss-Seidel on the interior nodes with the boundary held, is the
// standard answer and it takes the corpus from 57 inverted elements to 18. It
// is Laplace's equation with the roles of the parameters and the coordinates
// exchanged, so what it produces is the inverse of a harmonic map onto the
// square, which has no interior extremum to fold at. Holding the boundary is
// what keeps the mesh conforming while it does so -- the smoothing never moves
// a node that another patch can see.
//
// It is run only where it is needed; see Options::smoothingThreshold.
//
// ### The eighteen that are left
//
// They are not a failure of the smoother, and no smoother can remove them.
// Every one sits at a corner of a layout face whose interior angle on the model
// exceeds pi -- geom003's worst face has all four corners between 246 and 302
// degrees, a curved four-pointed region rather than anything like a rectangle.
// A structured grid takes that whole angle in its corner element, so the
// element is reversed as a matter of arithmetic. Report::reflexCorners counts
// them, and the remedy is the one Stage 9 already gives for the same face: the
// layout wants more patches there, which is Sec. 3.3's repair, upstream.
//
// ### What is not meshed
//
// Faces the arrangement could not close as quadrilaterals with one arc a side.
// They are counted, their area is reported as a fraction of S left uncovered,
// and they are left out; a face with three corners has no n_s x n_t grid and
// pretending otherwise would put a hanging node into an otherwise conforming
// mesh. Report::unmeshedPatches is the count and it is the same defect Stage 8
// already reported as wrongCornerFaces -- the remedy is upstream, in Sec. 3.3's
// repair, and not here.
class QuadMesh {
public:
    struct Options {
        // Target edge length, in the units of the model -- an absolute length,
        // not a fraction of anything. The corpus in data/meshes is normalised
        // into [-1,1]^2, so a model is two units across and 0.05 is one
        // fortieth of it. Report::modelExtent is the diagonal that came out, to
        // check this against on a model normalised some other way or not at
        // all.
        double targetEdgeLength = 0.05;

        // Floor and (optional) ceiling on the edges per arc. A chord always
        // gets at least the floor; 0 for maxIntervals means no ceiling.
        int minIntervals = 1;
        int maxIntervals = 0;

        // Closed loops of arcs that must come out with an **even** number of
        // edges in total. Empty unless something downstream needs one.
        //
        // The one thing that does is Stage 11: a rim it is asked to fill with
        // quadrilaterals must have an even number of edges on it, because
        // summing |boundary| over the faces of any quadrangulation of a disk
        // gives 4F = 2E_interior + E_boundary. That is not a property a fill
        // can be clever about -- an odd rim has no quadrangulation at all --
        // and it is not a property the chord assignment produces by accident:
        // on data/meshes/multimat/bubbles five of the ten inclusion rims come
        // out odd at the default target.
        //
        // It is cheap to arrange, though, because the parity is a linear
        // function over GF(2) of the per-chord counts: moving one chord by one
        // edge flips the parity of every loop that chord crosses an odd number
        // of times. See fixLoopParity, which solves that little system for the
        // cheapest set of chords to move and reports what the move cost.
        std::vector<std::vector<int>> evenLoops;

        // Points used to tabulate arc length along each fitted arc before the
        // nodes are placed on it by inverting that table. The fits are cubic
        // over three or four spans, so this is far finer than it needs to be
        // and costs nothing.
        int arcLengthSamples = 256;

        // Sample the Stage 9 spline fits. Off, the arcs are sampled along the
        // traced polylines instead and the interiors are a discrete Coons blend
        // of them -- the layout with no spline trusted anywhere, which is what
        // tells a meshing artefact apart from a fitting one.
        bool useSplines = true;

        // Place the nodes of a *feature* arc -- one on dS, or on the material
        // interface network -- on the traced polyline even when everything else
        // is placed on the Stage 9 fit.
        //
        // The two kinds of arc are not the same kind of object. A separatrix is
        // a curve the pipeline chose, and approximating it with three cubics is
        // a modelling decision the fit is entitled to make. A feature is a curve
        // the input gave: dS is where the model ends, an element that crosses an
        // interface carries two materials, and three cubics are not always
        // enough. geom011's interface is a cosine whose radius of curvature is
        // about three element lengths, the fit misses it by 4.8e-2 of the model
        // -- a whole element -- and 31 of its 984 elements come out with the
        // interface running through them. Meshed on the traced arc instead, none
        // do. The ICF hohlraum makes the same point on dS: a wall that runs
        // straight and then turns into a fillet inside one arc is cut by 1.1e-2.
        //
        // It costs nothing anywhere else: on a feature the fit follows, the two
        // curves agree to the fit's deviation, which is 1e-15 on the
        // straight-sided models in the corpus.
        //
        // Since SplineFit now carries those arcs exactly by default, this is
        // normally the same curve either way; it stays because it is what keeps
        // the node placement on the input when the arcs *are* put back under the
        // fit (SplineFit::Options::fitBoundaryArcs).
        bool featuresOnTracedArcs = true;

        // Stations per direction used to measure a patch's isoparametric line
        // lengths, which are the data the interval assignment is chosen from.
        int spanSamples = 16;

        // Winslow sweeps over the interior nodes of each block, the boundary
        // held. Zero leaves the transfinite grid as it is. A block stops early
        // once no node of it has moved by more than `smoothingTolerance` of the
        // target edge length.
        int smoothingPasses = 500;
        double smoothingTolerance = 1e-4;

        // Only blocks whose worst element is below this scaled Jacobian are
        // smoothed at all. The transfinite grid places its rows at the
        // arc-length spacing the interval assignment was chosen for, and
        // Winslow trades some of that away for interior regularity, so there is
        // no reason to spend it on a block that came out well. Zero means "only
        // where an element is actually inverted"; raise it to catch the nearly
        // folded too, or above 1 to smooth everything.
        double smoothingThreshold = 0.0;

        // Two mesh vertices closer than this fraction of the diagonal of S are
        // a crack: the same point of the model reached by two different
        // vertices. It should never fire, and it is checked rather than assumed
        // for the same reason Stage 9 measures its own watertightness.
        double crackTolerance = 1e-9;
    };

    // One class of the "opposite sides of a patch" relation: the arcs one chord
    // of the layout cuts across, and the single interval count they all take.
    struct Chord {
        std::vector<int> arcs;
        int intervals = 0;
        double idealIntervals = 0.0;  // N*, before rounding
        double minLength = 0.0;       // of its arcs
        double maxLength = 0.0;
        double minSpan = 0.0;         // of the patch isolines it was chosen from
        double maxSpan = 0.0;
        double minEdge = 0.0;         // minSpan divided by `intervals`
        double maxEdge = 0.0;
        bool clamped = false;         // the floor or the ceiling bound, not F(N)
    };

    // One meshed patch: a structured (ns+1) x (nt+1) block of vertex indices in
    // the patch's own (s, t) frame, row-major, so vert[j * (ns+1) + i] is the
    // node at (s_i, t_j). Kept because the block structure is the useful thing
    // downstream -- it is what a multiblock solver, a sweep, or a spline
    // refinement wants, and it cannot be recovered from the quad soup.
    struct Block {
        int face = -1;      // the face of the arrangement
        int patch = -1;     // its index in SplineFit::patches()
        int ns = 0, nt = 0;
        std::vector<int> vert;
    };

    struct Report {
        double modelExtent = 1.0;
        double target = 0.0;          // the requested edge length

        int chords = 0;
        int arcsAssigned = 0;
        int minIntervals = 0, maxIntervals = 0;
        int clampedChords = 0;        // bound by the floor or the ceiling
        double meanIntervals = 0.0;

        // Options::evenLoops: how many came out odd before the parity fix, how
        // many chords it had to move to make them even, and how much that cost
        // in the objective F the counts were chosen against. `oddLoopsLeft` is
        // the ones it could not fix -- a loop whose every chord also crosses
        // another odd loop an odd number of times, or one holding an arc no
        // patch asked for.
        int oddLoops = 0;
        int parityChordsMoved = 0;
        int oddLoopsLeft = 0;
        double parityCost = 0.0;

        int blocks = 0;               // patches meshed
        int unmeshedPatches = 0;      // faces without four sides of one arc
        double unmeshedArea = 0.0;    // of those, as a fraction of S

        int vertices = 0;
        int quads = 0;

        // What the interval assignment actually bought, as edge lengths on the
        // finished mesh. `edgeRatioRms` is the root mean square of
        // log(length / target) over every edge, which is the objective the
        // chords were chosen against, read back off the result.
        double minEdge = 0.0, maxEdge = 0.0, meanEdge = 0.0;
        double edgeRatioRms = 0.0;
        double worstEdgeRatio = 1.0;  // the furthest from 1, above or below

        // Element quality. The scaled Jacobian of a quad is the smallest of the
        // four corner cross products of its unit side vectors: 1 for a square,
        // 0 for a degenerate corner, negative for a folded element.
        double minScaledJacobian = 0.0;
        double meanScaledJacobian = 0.0;
        // The same, before the interior nodes were smoothed, so that what the
        // smoothing bought is visible rather than asserted.
        double minScaledJacobianBefore = 0.0;
        int invertedBefore = 0;
        int smoothedBlocks = 0;
        int smoothingSweeps = 0;      // the most any one block needed

        // Corners of a block whose interior angle on the model is more than pi.
        // A structured grid cannot cover one without reversing the element
        // there, so these are the inversions no smoothing can reach.
        int reflexCorners = 0;
        int invertedQuads = 0;        // non-positive area
        double minQuadArea = 0.0, maxQuadArea = 0.0;
        double meshArea = 0.0;
        double patchArea = 0.0;       // the same faces, as the arrangement has them

        // Conformity, measured rather than assumed. Every edge belongs to two
        // quads or to one; one is the mesh boundary, and there should be no
        // third. `cracks` is pairs of distinct vertices at the same point,
        // which is what a boundary meshed twice would leave.
        int interiorEdges = 0;
        int boundaryEdges = 0;
        int nonManifoldEdges = 0;
        int cracks = 0;

        // Materials, on a multi-material model. Every element carries the
        // material of the region it lies in, and an element that straddles an
        // interface -- centroid in one material, an edge midpoint in another --
        // is one no analysis code can integrate. Zero of those is the property
        // the whole multi-material path exists to produce, so it is measured on
        // the finished elements rather than inferred from the layout.
        int materials = 0;
        int mixedQuads = 0;
        // Elements every sample of which fell outside the triangulation, so the
        // material came from the nearest triangle instead. Not a defect -- a
        // Coons patch may bulge a fraction of an element past a curved piece of
        // dS -- but worth counting, because a large number of them means the
        // fit and the model have parted company.
        int unlocatedQuads = 0;
        // Element edges lying on a material interface. They are the ones two
        // materials share, and both sides carry the same nodes by construction
        // because the interface is a single arc of the layout.
        int interfaceEdges = 0;

        bool conforming = false;      // no third use of an edge, no cracks
        bool materialsPure = false;   // no element straddles an interface
        bool valid = false;           // ... and every patch meshed, none folded

        std::vector<std::string> messages;
    };

    explicit QuadMesh(const SplineFit &fit);
    QuadMesh(const SplineFit &fit, const Options &opts);

    const std::vector<Point>& vertices() const { return verts; }
    const std::vector<std::array<int, 4>>& quads() const { return cells; }
    // The material id of each element, from the region of the input mesh its
    // centroid lies in. All 1s on a single-material model.
    const std::vector<int>& quadMaterials() const { return cellMaterial; }
    const std::vector<Block>& blocks() const { return grids; }
    const std::vector<Chord>& chords() const { return chordList; }
    const Report& getReport() const { return report; }
    const Options& getOptions() const { return options; }

    // Edges assigned to each arc of the arrangement, or -1 where the arc bounds
    // no meshable patch, and the vertices placed along it from `from` to `to`.
    const std::vector<int>& arcIntervals() const { return intervals; }
    const std::vector<std::vector<int>>& arcVertices() const { return arcNodes; }
    // The chord each arc belongs to, or -1.
    const std::vector<int>& arcChord() const { return chordOf; }

    // The mesh as an .obj of quadrilateral faces.
    bool writeOBJ(const std::string &filename) const;
    // As a VTK unstructured grid of VTK_QUAD, carrying the scaled Jacobian and
    // the owning block per cell -- the picture of the block structure.
    bool writeVTU(const std::string &filename) const;

private:
    void assignIntervals();
    // Options::evenLoops, applied to the counts assignIntervals() chose.
    void fixLoopParity(std::vector<std::vector<double>> &chordSpans);
    void meshArcs();
    void meshPatches();
    void smooth();
    void classifyMaterials();
    void check();

    // Arc length along a fitted arc, tabulated at uniform parameter, and its
    // inverse: the parameter at which a given length has been travelled.
    struct Table {
        std::vector<double> cum;   // arcLengthSamples + 1 entries, cum[0] = 0
        double length = 0.0;
    };
    Table tabulate(int arc) const;
    // The mean length of the isoparametric lines of a patch, per direction.
    void patchSpans(int patch, const std::vector<Arrangement::Side> &sides,
                    double &sMean, double &tMean) const;
    double paramAtLength(const Table &t, double s) const;
    Point evaluateArc(int arc, double u) const;

    const SplineFit *fit = nullptr;
    const Arrangement *arr = nullptr;
    Options options;

    std::vector<Point> verts;
    std::vector<std::array<int, 4>> cells;
    std::vector<int> cellMaterial;
    std::vector<int> cellBlock;             // per quad, its Block
    std::vector<Block> grids;
    std::vector<Chord> chordList;

    std::vector<int> intervals;             // per arc
    std::vector<int> chordOf;               // per arc
    std::vector<std::vector<int>> arcNodes; // per arc, intervals + 1 vertices
    std::vector<std::vector<double>> arcParams; // ... and their curve parameters
    std::vector<int> nodeVert;              // arrangement node -> vertex
    std::vector<int> facePatch;             // arrangement face -> SplineFit patch
    std::vector<Table> tables;              // per arc

    double modelExtent = 1.0;
    Report report;
};

#endif // __QUAD_MESH_HXX__
