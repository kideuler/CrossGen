#ifndef __POLYSQUARE_HXX__
#define __POLYSQUARE_HXX__

#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Parameterization/HarmonicCut.hxx"
#include "UMBER/UMBER.hxx"
#include "mesh/Mesh.hxx"

// Frame field guided polysquare generation of
//   Wang, Ren, Fang, Lin, Xu, Bao and Huang, "IGA-suitable planar
//   parameterization with patch structure simplification of closed-form
//   polysquare", CMAME 392 (2022) 114678, Section 4.3.
//
// Takes the frame field UMBER optimized (Sec. 4.2) and the cuts HarmonicCut
// made (Sec. 4.1) and deforms the cut mesh M_C into the parameter domain, so
// that the boundary comes out axis aligned and the two banks of every cut
// agree up to a k*90-degree rotation. The result is the polysquare: one (u, v)
// per vertex of M_C.
//
// Four stages: the three of Sec. 4.3, and the boundary half of the Sec. 5.5
// re-solve, which is what makes the alignment exact rather than close.
//
//   Eq. (6)  A Poisson solve that asks the deformation gradient to be the
//            frame, grad phi = (v v^perp)^T. Linear, and only an initial
//            value: it satisfies neither the seamless condition nor the axis
//            alignment.
//   Eq. (7)  The transition Pi_gamma of each cut, the k*90-degree rotation
//            that best lines the frames up across it. These are the harmonic
//            degrees of freedom of Sec. 4.1 -- the difference between a common
//            polysquare and a closed-form one -- so they are extracted from
//            the field rather than assumed to be the identity.
//   Eq. (9)  min E_arap + w_l1 E_l1 + w_cor E_cor subject to the transitions,
//            a positive Jacobian everywhere, and a fixed boundary length.
//            E_arap (Eq. 10) holds the distortion down, E_l1 (Eq. 11) is what
//            snaps each boundary edge onto an axis, and E_cor (Eq. 12) keeps
//            the corners where the boundary actually turns instead of letting
//            E_l1 place them on a straight run.
//   Eq. (23) Each boundary segment put exactly on its axis and the interior
//            re-solved with it held there. Eq. (9) aligns by a penalty, so it
//            lands near an axis, and Sec. 5 needs the segments to be on one:
//            its Eq. (13) defines a segment's projection h_i only for a
//            segment that really is axis aligned. See snapBoundary().
//
// Four places where this departs from the paper, all deliberate:
//
//  - The transition constraints of Eq. (8) are eliminated rather than imposed.
//    Each cut vertex of M_C has two copies; one bank is a variable and the
//    other is defined as Pi_gamma times it plus a per-cut translation, which
//    is the same condition with two unknowns per cut instead of one
//    constraint per cut edge. It also lets the Poisson stage be seamless from
//    the start, where the paper's Eq. (6) is not.
//
//  - The target angles of Eq. (12) come from the frame field rather than from
//    the polar angles of [57], which is not implemented here. UMBER already
//    reports a corner index per boundary vertex -- the quarter turns the frame
//    makes against the boundary -- and that is the same quantity Eq. (12)
//    wants: 90 degrees where the field puts a corner, 0 along a run it
//    follows. At a geometric corner the two agree by construction, since the
//    field is aligned on both sides of it.
//
//    The naive alternative, taking the target from the boundary's own turning
//    angle, does not work and is worth recording: on a smooth arc every vertex
//    then asks for 0, so the term demands that the whole arc stay straight in
//    the image, which is impossible where the arc has to turn a corner. At the
//    paper's w_cor = 10 that demand overwhelms E_l1 and the boundary never
//    snaps -- measured across data/meshes, 1.2 to 13.6 degrees of mean
//    misalignment with every model's worst edge sitting at 41 to 45 degrees,
//    i.e. no polysquare at all. Whatever [57] computes, it cannot be the raw
//    turning angle; it has to already distribute the turning into quarter
//    turns, which is what the frame field's corner index does.
//
//  - det(grad phi) > 0 is a penalty on the triangles that violate it rather
//    than the log barrier of [11]. A barrier needs a feasible starting point
//    and the Poisson initialization is not always one; the penalty is defined
//    for flipped triangles too, so it can pull them back. Its weight climbs
//    until nothing is flipped, and getReport().flips says whether that worked.
//
//  - The boundary length constraint of Eq. (9) is a gauge rather than a
//    constraint: E_arap already fixes the scale, since the frame is unit
//    length and so grad phi is asked to be a rotation. getReport() carries the
//    length ratio so the drift is visible.
class Polysquare {
public:
    struct Report {
        int flips = 0;                  // triangles with det(grad phi) <= 0
        double minScaledJacobian = 0.0; // over triangles, det / (|col1| |col2|)
        double avgScaledJacobian = 0.0;
        double meanAlignDeg = 0.0;      // boundary edges vs the nearest axis
        double maxAlignDeg = 0.0;
        // Turns of the image boundary, counted over dM only (the cut banks are
        // seams inside the polysquare, not boundary, and are free to wander).
        // This is the number to look at, not the alignment: a staircase is
        // perfectly axis aligned and completely wrong, so only `turns` can
        // tell the two apart. It should come out equal to `expectedTurns`, the
        // number of boundary corners the frame field asked for.
        int turns = 0;
        int expectedTurns = 0;
        int shortestRun = 0;            // edges in the shortest straight run
        // How far a boundary segment ran from its axis at the moment it was
        // snapped onto it. A few degrees is discretization. A large one says
        // Eq. (9) left a diagonal run, and the snap straightened something
        // whose corners really belong somewhere else -- the alignment will
        // read clean afterwards and the structure will still be wrong, so this
        // is the number that catches it.
        double worstSegmentDeg = 0.0;
        int suspectSegments = 0;        // segments past 10 degrees
        // Corner coordinates moved onto a shared iso-line by snapCorners(),
        // the largest such move in mean image boundary edges, and the clusters
        // it had to leave alone because two corners it could not move disagreed.
        int cornersSnapped = 0;
        double worstCornerSnap = 0.0;
        int cornerSnapConflicts = 0;
        // Straight runs of the boundary moved onto a coordinate they share with
        // another run, and the largest such move in mean image boundary edges.
        // A run moves rigidly, so this costs no alignment: it is the corners at
        // its ends that it is for.
        int runsAligned = 0;
        double worstRunAlign = 0.0;
        double transitionDeg = 0.0;     // worst residual of Eq. (7) across a cut
        double lengthRatio = 1.0;       // image boundary length / input length
        double arap = 0.0, l1 = 0.0, cor = 0.0;
        int iterations = 0;
    };

    // Both must be built on the same mesh, and `frames` must have been
    // optimized -- the field is the guide for the whole deformation.
    Polysquare(const UMBER &frames, const HarmonicCut &cut);

    // w_l1 of Eq. (9), one entry per continuation stage. The paper starts at
    // 0.125 and doubles until past 1.0, warm-starting each stage, because
    // starting large lands in a poor local minimum and lets the shape distort
    // (Sec. 7.1, Fig. 14).
    //
    // The stage past 1.0 is the one that matters, and stopping at 1.0 does not
    // give a polysquare: on data/meshes/geom010 the boundary is still 1.9
    // degrees off axis on average and 19 at worst, with a visible diagonal run
    // where a staircase of axis-aligned segments belongs. The 2.0 stage brings
    // that to 0.32 and 1.5. Pushing further keeps buying alignment and pays
    // for it in distortion -- through 8.0 the mean falls to 0.08 degrees while
    // the minimum scaled Jacobian drops from 0.51 to 0.07 -- which is why the
    // schedule stops at the first weight past 1.0 rather than continuing.
    void setL1Schedule(const std::vector<double> &schedule) { l1Schedule = schedule; }

    // w_cor of Eq. (9), the paper's 10.
    //
    // This term is not a refinement, it is what makes the boundary a
    // polysquare boundary. E_l1 asks each edge to lie on an axis and is just
    // as happy with a staircase as with a straight run, since every step of a
    // staircase is axis aligned too -- so E_l1 alone gives a boundary that
    // turns almost every vertex. Measured with w_cor = 0 on data/meshes: the
    // half disk, whose frame field asks for four corners, comes out with 36
    // turns, and the median straight run is zero edges. E_cor is what says
    // "and do not turn here", by asking consecutive boundary edges to be
    // parallel wherever the field reports no corner.
    //
    // At 10 the turn count lands exactly on the number of corners the field
    // asked for, on all sixteen models in data/meshes, holed ones included --
    // 4 turns for 4 corners, 36 for 36 -- with straight runs of tens of edges
    // on the simple shapes. Anything from 1 upwards does that; the paper
    // reports 10 to 100 as the useful band and degeneracy at 1000.
    //
    // 0 disables it, and is worth knowing about only because the two terms do
    // pull against each other: the mean per-edge misalignment roughly triples
    // when it is switched on (0.1 degrees to 0.3), and on one model it leaves
    // a diagonal run at 21 degrees. That trade is worth taking, since a
    // perfectly aligned staircase is not a polysquare and a slightly rounded
    // corner is.
    void setCornerWeight(double w) { corWeight = w; }

    // Whether to finish by putting every boundary segment exactly on its axis
    // and re-solving the interior with it held there -- the boundary half of
    // the paper's Sec. 5.5 re-parameterization, Eq. (23).
    //
    // On by default, because "axis aligned" is a property the rest of the
    // pipeline needs to hold exactly rather than nearly: Eq. (13) defines a
    // segment's projection h_i only for a segment that is genuinely axis
    // aligned, and everything Sec. 5 does rests on that. Eq. (9) gets close --
    // most edges to within a tenth of a degree -- but close is a different
    // thing, and on one model in data/meshes it leaves a fourteen-edge run
    // bulging out to 21 degrees where E_l1 and E_cor deadlock.
    //
    // Turning it off leaves the raw Eq. (9) result, which is what to look at
    // when asking how well the soft alignment did on its own.
    void setSnapBoundary(bool on) { snapBoundaryOn = on; }

    // How close two corners have to be, in mean image boundary edges, for them
    // to be treated as sharing an iso-line -- see snapCorners(). 0 disables it.
    //
    // The discrepancy this is for is small: on data/meshes/geom021 the two
    // corners that should share an iso-line come out 1.4% of a boundary edge
    // apart, and that is the one that costs a whole partition. So the default
    // is a twentieth of an edge -- three times what is needed there and still
    // far under anything the mesh can resolve.
    //
    // Larger is not better, and it is worth saying why, because the instinct is
    // to leave headroom. Past a point the clusters stop being one corner seen
    // twice and start being two corners, and merging those folds the boundary
    // over itself. Measured over data/meshes: at 0.05 no model gains a block
    // that is not four-sided; at 0.10 geom008 and geom009 both do; at 0.25
    // geom008 also gains a pair of sides that cross each other.
    //
    // The trade is not free in the other direction either. On geom012 a wider
    // tolerance leaves the chord collapse more to work with -- 110 blocks
    // against 145 -- because at 0.05 the small shifts it makes are enough for
    // six of that model's collapses to be rolled back for crossing rather than
    // four. That is a heuristic yielding less, though, and the alternative is a
    // structure with blocks in it that are not four-sided at all, so it is the
    // narrower tolerance that is kept.
    void setCornerSnapTolerance(double t) { cornerSnapTol = t; }

    // L-BFGS iterations per continuation stage, and the gradient tolerance.
    void setMaxIterations(int n) { maxIterations = n; }
    void setTolerance(double t) { gradTolerance = t; }

    // Eq. (6), then Eq. (7), then Eq. (9) over the schedule.
    void solve();

    // The polysquare: one (u, v) per vertex of the cut mesh, in the same order
    // as getCutMesh().vertices.
    const std::vector<Point>& getUV() const { return uv; }

    // M_C, whose triangles carry the same indices as the input mesh's.
    const Mesh& getCutMesh() const { return *cutMesh; }

    // The mesh the whole thing was built on.
    const Mesh& getMesh() const { return *orig; }
    std::shared_ptr<Mesh> getMeshPtr() const { return orig; }

    // The frame field's corner index per original vertex, in quarter turns,
    // which is where the boundary of the polysquare turns: +1 convex, -1
    // reflex, 0 everywhere else.
    const std::vector<int>& getBoundaryCorners() const { return boundaryCorner; }

    // k of Pi_gamma = R(k*90 degrees) per cut, in the order HarmonicCut made
    // them. All zero means the field asked for a common polysquare; anything
    // else is the closed form using its harmonic degrees of freedom.
    const std::vector<int>& getTransitions() const { return transitionK; }

    const Report& getReport() const { return report_; }

    // The parameterized mesh as a VTK unstructured grid: points at (u, v, 0),
    // the triangles of M_C, the original position as a point field so the two
    // domains can be compared, and per-triangle distortion as cell data.
    bool writeVTU(const std::string &filename) const;

    // The input mesh with the parameterization carried as a point field, for
    // looking at where the polysquare's corners land on the model.
    bool writeSourceVTU(const std::string &filename) const;

private:
    // Per triangle: the gradients of the three hat functions and the area.
    struct TriGrad {
        Point g[3];
        double area = 0.0;
        double weight = 0.0; // area / total area
    };

    // One boundary edge of Eq. (11), in cut-mesh vertices, walked in loop
    // order so that consecutive entries share a vertex for Eq. (12).
    struct BoundaryEdge {
        int ca = -1, cb = -1;
        double length = 0.0;
        int sharedOrigVertex = -1; // the vertex this edge ends at
        double targetTurn = 0.0;   // theta_i of Eq. (12) at that vertex
        // -theta_gamma where this edge and the next sit on opposite banks of a
        // cut. Part of targetTurn, and kept separately because measuring the
        // structure has to discount it: across a seam the two edges use
        // different copies of the shared vertex, so their images do not even
        // meet, and the angle between them is the transition rather than a
        // corner of the polysquare.
        double seamTurn = 0.0;
        bool cornerPair = false;   // whether Eq. (12) uses this pair
    };

    void buildTopology();
    void extractTransitions();   // Eq. (7)
    void poissonInit();          // Eq. (6)
    void optimize();             // Eq. (9)
    // One straight run of the image boundary: the vertices on it, the axis it
    // is perpendicular to, and the single coordinate they are all given.
    struct Run {
        int axis = 0;
        double h = 0.0;
        double weight = 0.0;      // image length, so a long run outvotes a short one
        std::vector<int> verts;   // cut-mesh vertices
    };

    void snapBoundary();         // Eq. (23), the boundary constraint
    // Give runs that lie on one iso-line one coordinate between them.
    void alignRuns(std::vector<Run> &runs);
    void snapCorners();          // corners that share an iso-line put on one
    int countFlips() const;
    // The mean length of a boundary edge in the image: the scale everything
    // about "the same iso-line" is measured in.
    double meanImageBoundaryEdge() const;

    // The triangle that has the interior on its left when a -> b is walked,
    // which is what makes "the left bank" mean the same thing along a whole
    // cut or boundary loop.
    int leftTriangleOf(int a, int b) const;
    // The triangle on the other side of the edge (a, b) from f.
    int oppositeTriangle(int f, int a, int b) const;
    // The cut-mesh vertex triangle f uses for the original vertex v.
    int cornerOf(int f, int origVertex) const;

    // x (the reduced variables) -> uv (every cut-mesh vertex).
    void expand(const Eigen::VectorXd &x);
    // grad over uv -> grad over x. Consumes gradUV, folding the dependent
    // copies into the banks they follow.
    void fold(std::vector<Point> &gradUV, Eigen::VectorXd &grad) const;
    // Eq. (9) and its gradient at x.
    double evaluate(const Eigen::VectorXd &x, Eigen::VectorXd &grad,
                    double *arapOut = nullptr, double *l1Out = nullptr,
                    double *corOut = nullptr);
    int runLBFGS(Eigen::VectorXd &x, int maxIter);
    void measure();

    std::shared_ptr<Mesh> orig;
    const Mesh *cutMesh = nullptr;
    const HarmonicCut *harmonicCut = nullptr;

    std::vector<Point> frameU, frameV; // per triangle, the rows of Eq. (6)'s target
    std::vector<TriGrad> triGrad;
    double totalArea = 0.0;
    double totalBoundaryLength = 0.0;

    // theta_i of Eq. (12) per original boundary vertex, in quarter turns, as
    // the frame field reported it.
    std::vector<int> boundaryCorner;

    std::vector<int> transitionK;      // per cut
    std::vector<BoundaryEdge> bEdges;  // in loop order, loops back to back
    std::vector<int> loopStart;        // index into bEdges of each loop's first edge

    // Variable layout. varOf[cv] >= 0 indexes the free variable, -1 marks a
    // vertex defined by a transition, -2 the one pinned vertex that removes
    // the global translation.
    std::vector<int> varOf;
    std::vector<int> depSrc;  // cut vertex -> the copy it follows
    std::vector<int> depCut;  // cut vertex -> which cut's transition applies
    int nIndep = 0;
    int nCuts = 0;
    int pinnedVertex = -1;

    Eigen::VectorXd x;        // 2*nIndep free coordinates, then 2*nCuts translations
    std::vector<char> fixedX; // components snapBoundary() nailed down
    std::vector<Point> uv;    // per cut-mesh vertex

    std::vector<double> l1Schedule{0.125, 0.25, 0.5, 1.0, 2.0};
    double l1Weight = 0.125;  // the stage of that schedule currently running
    double l1Eps = 1e-2;
    double corWeight = 10.0;
    bool snapBoundaryOn = true;
    double cornerSnapTol = 0.05;
    double barrierWeight = 1.0;
    double detFloor = 0.05;
    double gradTolerance = 1e-6;
    // Per stage. None of data/meshes reaches the gradient tolerance, so this
    // is what actually ends each stage; 1000 and 3000 land within 0.01 degrees
    // of each other on the alignment, so it is past the point of diminishing
    // returns rather than a cap the result is fighting.
    int maxIterations = 1000;

    Report report_;
};

#endif // __POLYSQUARE_HXX__
