#ifndef __SKELETON_HXX__
#define __SKELETON_HXX__

#include <memory>
#include <vector>

#include "SKELETON/BaseComplex.hxx"
#include "SKELETON/MedialAxis.hxx"
#include "SKELETON/MedialAxisMap.hxx"
#include "SKELETON/MedialAxisTMesh.hxx"
#include "mesh/BlockQuadMesh.hxx"
#include "mesh/Mesh.hxx"
#include "mesh/QuadMesh.hxx"
#include "mesh/TMOP.hxx"
#include "quantization/QuantTMeshConvert.hxx"
#include "quantization/TMeshQuantizer.hxx"

// SKELETON: the medial-axis quad mesher of
//   Cheny et al., "Planar quad meshing with guarantees"   (docs/medial_axis.md)
// taken on to a block decomposition, a mesh and a smoothed mesh the way every
// other method in this codebase is finished.
//
//   Steps 1-2  triangulate()     the constrained Delaunay triangulation of dS's
//                                own vertices, no interior points (Step 3.1)
//   Steps 3-4  buildAxis()       MedialAxis, the polygonal dual complex, and
//                                its simplify(): endpoints pruned, arcs
//                                flattened, within delta_max = 0.05 h
//   Step 5     buildMap()        MedialAxisMap, the boundary map phi
//   Steps 6-8  buildTMesh()      MedialAxisTMesh: downsampling, colouring, the
//                                Fig. 17 templates, the Step 8.6 check
//   Step 9     quantize()        the T-mesh welded (makeQuantTMesh) and
//                                quantized with Prop. 1's requirements
//   Step 10    extractQuads()    the quad mesh those lengths prescribe
//   Blocks     buildBlocks()     the Step 10 mesh refined and smoothed, and its
//                                base complex as the shared BlockDecomposition
//                                (SKELETON/BaseComplex.hxx)
//   Mesh       buildMesh()       mesh/BlockQuadMesh on the blocks
//   Smoothed   smooth()          mesh::TMOP, mu sampled at the corners
//
// The paper stops at Step 10, with a quad mesh at the target size h. The last
// three stages are what makes SKELETON comparable with the other methods: its
// blocks are drawn, graded (BlockDecomposition::regularity() and the rest) and
// meshed at a target edge length the way ZIPLINE's and UMBER's are, so that a
// SKELETON mesh and a ZIPLINE one at one number differ in their layouts and in
// nothing else. h is then the layout's scale and not the elements': it decides
// how finely the axis is sampled and so where the singularities go, and the
// mesh's own size is BlockQuadMesh's target edge length.
//
// SKELETON could mesh through its own geometry instead, as MERIDIAN does
// through Stage 9's splines: Step 10 again with every quad cut k x k, which
// keeps the T-mesh faces' Coons patches inside each block. That was measured
// on 2026-10-17 over singlemat and multimat at a target of 0.05, both meshes
// then given 1000 TMOP sweeps, and it was not better. It left folds on 11
// models, as BlockQuadMesh did, and it was worse on some (geom034 6 against 0,
// rocket 6 against 2). Its worst element was worse on most models (geom011
// 0.25 against 0.62), since Step 10's template slivers come through it. Its
// element count misses the target by up to 4.6x, since k is an integer
// (geom031: 4896 quads against 1053).
//
// ### Why one class
//
// Until 2026-10-17 the viewer's mode 3 and TestSKELETON each ran these steps
// themselves, with their own copy of Steps 1-2, which took the first boundary
// loop found to be the outer one: on geom014 and multimat/turbine_blade that
// loop is a hole, and the triangulation came out of the wrong region. Both now
// run this class, so a model gives the same axis, the same blocks and the same
// mesh in both, as ZIPLINE's driver made it for mode 2.
//
// ### Running it
//
// As ZIPLINE: run() takes the model as far as the stage it is given; each
// stage can be run on its own, and first runs whatever it needs that has not
// run; running a stage again rebuilds it from the stage before and discards
// every stage after it.
//
// Material interfaces are not honoured: Steps 1-2 see dS only, so a
// multi-material model is meshed as one material.
class SKELETON {
public:
    enum class Stage { Triangulation, Axis, Map, TMesh, Quantized, QuadMesh, Blocks, Mesh, Smoothed };

    struct Options {
        // h, the target size of Steps 4 and 6-9. <= 0 takes a tenth of the
        // bounding-box diagonal (MedialAxis::defaultTargetSize()), the guide's
        // default and the one every corpus number in the memory note was
        // measured at. It overrides the targetSize fields of the two below.
        double targetSize = 0.0;
        SimplifyOptions simplify;
        MedialAxisTMesh::Options tmesh;

        // ---- Step 9 ----
        // minLength 0 and the absolute objective, against the quantizer's
        // defaults (ZIPLINE's): with every length held at 1 the templates'
        // slivers cost an element each and 3-18 edges were forced to zero
        // on several models, and the relative objective lets the tiny
        // template edges dominate. Validity comes from Prop. 1's
        // requirements instead, which quantize() fills in.
        TMeshQuantizer::Options quantizer = [] {
            TMeshQuantizer::Options q;
            q.minLength = 0;
            q.objective = TMeshQuantizer::Options::Objective::Absolute;
            return q;
        }();

        // ---- Blocks ----
        // The Step 10 mesh is traced at k times its resolution, so that the
        // block sides follow the T-mesh faces' Coons geometry at h / k
        // rather than as chords h long.
        int blockSubdivisions = 4;
        // That trace mesh smoothed by TMOP before its base complex is read,
        // the model's corners held and the other dS nodes sliding. Step 10
        // places the mesh by one Coons patch per T-mesh face, the patches
        // agree only where two faces meet, and a face quantized to zero is
        // welded away: a grid line kinks at every face it crosses, and the
        // block sides read off it as traced zigzag between their corners.
        // Smoothed, every quad of the trace mesh is pulled towards a square of
        // one size -- the quantized layout's own flat metric, each block a
        // rectangle of x k by y k of them -- and the sides are the straight
        // lines of that metric as the embedding that least distorts it lays
        // them down: smooth curves, the macrovertices moving with them. The
        // topology is the mesh's, so the blocks are the same blocks.
        //
        // Measured 2026-10-17 over singlemat and multimat (66 models), against
        // the sides as traced: kinks (a vertex turning 15 degrees more than
        // the mean of its neighbours) 1578 -> 427, in 41 -> 28 models;
        // angleQuality 0.786 -> 0.852, chordQuality 0.569 -> 0.625,
        // regularity 0.515 -> 0.507; after BlockQuadMesh at 0.05 and 1000
        // TMOP sweeps, the worst element 0.447 -> 0.486 on average and 11 -> 10
        // models still folded. Over MAMBO's 407: angle 0.858 -> 0.923, chord
        // 0.710 -> 0.748. On geom007, the model it was done for, 23 kinks -> 6.
        //
        // Geodesics of a conformal metric were tried first, from the cone
        // angles of the blocks' own corners as TORSION's ConeMetric solves for
        // them: on geom007 13 of the 21 interior sides came out leaving the
        // model or crossing another side. A conformal metric has the layout's
        // angles but not its lengths -- it cannot make a block of 1 x 8 quads
        // a long thin rectangle -- so its geodesic between two corners runs
        // wherever the scale is small, into the holes. A field's metric
        // (DualMBO's) has the same defect and cones the blocks do not have.
        //
        // ShapeSize007 rather than the pure shape metric: Shape002 lets sizes
        // drift to make the quads square, which is a conformal map again, and
        // left geom007's sides more bent (0.57 rad of turning against 0.48).
        // 1000 untangling sweeps rather than 100: the trace mesh comes in with
        // up to a few dozen quads folded over, and at 100 geom034 stayed
        // tangled (427 kinks over the corpus against 480). A trace mesh that
        // still has a fold when TMOP stops is read as Step 10 left it --
        // geom003 and multimat/geom003, the exact discs, and multimat/geom004
        // -- as is one that is not a manifold (tooth). maxSweeps 0 turns it
        // off (TestSKELETON --no-block-smoothing).
        mesh::TMOP::Options blockSmoothing = [] {
            mesh::TMOP::Options t;
            t.metric = mesh::TMOP::ShapeSize007;
            t.quadrature = mesh::TMOP::Corners;
            t.maxSweeps = 500;
            t.untangleMaxSweeps = 1000;
            return t;
        }();
        // When the smoothed trace mesh stays tangled, or a block side read off
        // it crosses dS (BaseComplex::Report::sidesCrossingBoundary), it is
        // smoothed again from Step 10's positions with its boundary laid onto
        // dS in this many stages -- each loop of it spaced along its loop of
        // dS between the model's corners as it is spaced along itself, the
        // whole boundary held and TMOP untangling behind it -- and smoothed
        // once more with the boundary still held. If that does not untangle
        // either, the interior is placed afresh at the mean of its neighbours
        // with the boundary on dS (Tutte's embedding, fold-free when dS is one
        // convex loop: geom003, an exact disc, whose caps Step 10 folds over
        // by the hundred) and smoothed from there. The blocks are read off
        // whichever untangled with the fewest sides across dS.
        //
        // How a side comes to cross: ShapeSize007 asks every quad to be a
        // square of one size, a channel thinner than that size -- geom007's
        // under its right hole -- pushes its rows apart, and with the
        // boundary sliding the nodes at the bottom of the hole slide away
        // along it until one boundary edge spans the circle's cap as a chord
        // and the grown rows fill the cap. No element is inverted, so the
        // smoother sees nothing wrong, and the block side in that row cuts
        // across the hole. Elsewhere it is Step 4's moves and the welds of
        // zero-length edges leaving the boundary off dS by up to h / 2, which
        // BaseComplex's snap afterwards moves while leaving the rows behind.
        // Held on dS throughout, the boundary keeps the spacing Step 10 gave
        // it and the rows have to stay on the model's side of it; putting it
        // on dS in one step instead folds the rows behind it past what the
        // untangler undoes (geom012, geom009).
        //
        // Measured 2026-10-17, with Step 10's doublets merged and its quads
        // oriented by connectivity (extractQuadMesh) and SKELETON::smooth()
        // keeping its elements' orientation. Singlemat and multimat: sides
        // crossing dS 61 in 12 models -> 27 in 4, geom007's gone; folded after
        // BlockQuadMesh and TMOP, counted with the elements' own orientation,
        // 1 model (tooth, one element) where 8 were, the worst element 0.486
        // -> 0.624 on average; the four scores within 0.01. MAMBO: crossing 72
        // in 20 -> 32 in 3, folded 1 (M1_10). What is left is a concave corner
        // the mesh rounds off or pinches (geom014's square holes, rocket,
        // tooth), which is Step 11's to answer, and M1_10, whose Step 10 mesh
        // nothing here untangles.
        //
        // Only where a side crosses, because a held boundary keeps Step 10's
        // spacing, which a box corner's blue template crowds (geom006:
        // angleQuality 0.93 -> 0.62 when every model was staged; the corpus
        // 0.852 -> 0.775). Letting the boundary slide during the stages kept
        // the quality and let the crossings back. Also tried and refused: the trace mesh bent onto the walls'
        // cross field (FeatureFrame) in the least-squares sense first, each
        // quad turned towards the field as far as the walls agree and every
        // boundary vertex put back on dS after each solve, as-rigid-as-
        // possible style. It bends the rows round the holes as asked, but
        // nothing in least squares stops a fold, and it folded what the
        // untangler then could not undo (24 models left tangled, against 5,
        // crossings 45 rather than 30). 0 turns this off.
        int boundaryStages = 4;
        BaseComplex::Options blocks;
        // A model corner is a vertex of dS turning by more than this many
        // degrees, BoundaryFeatures' (G3's) 20-degree rule, measured on dS
        // alone.
        double cornerAngle = 20.0;

        // ---- Mesh ----
        BlockQuadMesh::Options mesh;

        // ---- Smoothed ----
        // ZIPLINE's settings, for ZIPLINE's reason: a transfinite grid on a
        // block decomposition can turn a corner over without the 2x2 Gauss
        // points seeing it.
        mesh::TMOP::Options tmop = [] {
            mesh::TMOP::Options t;
            t.metric = mesh::TMOP::ShapeSize007;
            t.quadrature = mesh::TMOP::Corners;
            t.maxSweeps = 1000;
            t.untangleMaxSweeps = 2000;
            return t;
        }();
    };

    struct Status {
        // Steps 3-4: the dual complex before and after simplify().
        int rawMedialVertices = 0, rawMedialEdges = 0;
        double targetSize = 0.0;      // h as used
        int boundaryLoops = 0;
        int corners = 0;              // of dS, by Options::cornerAngle
    };

    explicit SKELETON(std::shared_ptr<Mesh> mesh);
    SKELETON(std::shared_ptr<Mesh> mesh, const Options &opts);
    ~SKELETON();

    // Every stage up to and including `last`. False when a stage had nothing
    // to work on: a T-mesh that would not weld, an inconsistent
    // quantization, no quad to read blocks off, no block to mesh.
    bool run(Stage last = Stage::Blocks);

    // ---- one stage at a time ----
    // Steps 1-2. Throws when the model has no boundary loop.
    void triangulate();
    void buildAxis();                         // Steps 3-4
    void buildMap();                          // Step 5
    void buildTMesh();                        // Steps 6-8
    // Step 9. False when the blocks would not weld or the lengths found
    // break the consistency system.
    bool quantize();
    // Step 10. False when there is no consistent quantization to read.
    bool extractQuads();
    // False when Step 10 left no quad.
    bool buildBlocks();
    // False when there is no block to mesh.
    bool buildMesh();
    bool buildMesh(const BlockQuadMesh::Options &opts);
    // False when there is no element to smooth.
    bool smooth();

    // ---- what the stages left ----
    bool hasTriangulation() const { return triangulation != nullptr; }
    bool hasAxis() const { return axis != nullptr; }
    bool hasMap() const { return map != nullptr; }
    bool hasTMesh() const { return tmesh != nullptr; }
    bool hasQuantization() const { return quant != nullptr; }
    // A consistent one, which is what Step 10 needs.
    bool quantizationOk() const { return quant && quant->ok && quantReport.consistent; }
    bool hasQuadMesh() const { return quads != nullptr; }
    bool hasBlocks() const { return blocks != nullptr; }
    bool hasBlockMesh() const { return blockMesh != nullptr; }
    bool hasSmoothedMesh() const { return smoothedMesh != nullptr; }

    // The triangulation the axis was built on. Step 4 moves its boundary
    // vertices (within delta_max), so it is not the input's outline any more.
    const Mesh &getTriangulation() const { return *triangulation; }
    const MedialAxis &getAxis() const { return *axis; }
    const MedialAxisMap &getMap() const { return *map; }
    const MedialAxisTMesh &getTMesh() const { return *tmesh; }
    const BlockQuant &getQuantization() const { return *quant; }
    const TMeshQuantizer::Report &getQuantReport() const { return quantReport; }
    // Step 10's mesh, at the target size h.
    const QuadExtraction &getQuadMesh() const { return *quads; }
    const BaseComplex &getBlocks() const { return *blocks; }
    const BlockDecomposition &getDecomposition() const { return blocks->decomposition(); }
    // How the trace mesh the blocks were read off was placed (buildBlocks()
    // tries them in this order), and the smoothing's report; ran is false
    // for AsTraced.
    enum class BlockGeometry {
        AsTraced,    // Step 10's positions: smoothing off, or nothing untangled
        Smoothed,    // TMOP, the boundary sliding
        Staged,      // the boundary laid onto dS in stages and held
        Relaxed      // the interior relaxed afresh from dS, the boundary held
    };
    BlockGeometry blockGeometry() const { return geometry; }
    static const char *name(BlockGeometry g) {
        switch (g) {
            case BlockGeometry::AsTraced: return "as Step 10 left it";
            case BlockGeometry::Smoothed: return "smoothed";
            case BlockGeometry::Staged:   return "smoothed with its boundary laid onto dS in stages";
            case BlockGeometry::Relaxed:  return "relaxed afresh from dS and smoothed";
        }
        return "?";
    }
    const mesh::TMOP::Report &getBlockSmoothingReport() const { return blockSmoothingReport; }
    const BlockQuadMesh &getBlockMesh() const { return *blockMesh; }
    const mesh::QuadMesh &getSmoothedMesh() const { return *smoothedMesh; }
    const mesh::TMOP::Report &getTMOPReport() const { return tmopReport; }

    // dS as closed polylines of the input's vertices (first point not
    // repeated), the outer loop first, and its corners: what Steps 1-2
    // triangulate and what the blocks are laid on.
    const std::vector<std::vector<Point>> &boundaryLoops() const { return loops; }
    const std::vector<Point> &boundaryCorners() const { return corners; }

    const Options &getOptions() const { return options; }
    const Status &getStatus() const { return status; }
    const Mesh &getMesh() const { return *mesh; }
    std::shared_ptr<Mesh> getMeshPtr() const { return mesh; }

private:
    void discardFrom(Stage from);
    double h() const;
    // Options::blocks with the snap tolerance raised to h / 2.
    BaseComplex::Options blockOptions() const;

    std::shared_ptr<Mesh> mesh;
    Options options;
    Status status;

    std::vector<std::vector<Point>> loops;
    std::vector<Point> corners;

    // In stage order, each built from the one before it.
    std::shared_ptr<Mesh> triangulation;
    std::shared_ptr<MedialAxis> axis;
    std::unique_ptr<MedialAxisMap> map;
    std::unique_ptr<MedialAxisTMesh> tmesh;
    std::unique_ptr<BlockQuant> quant;
    TMeshQuantizer::Report quantReport;
    std::unique_ptr<QuadExtraction> quads;
    std::unique_ptr<BaseComplex> blocks;
    mesh::TMOP::Report blockSmoothingReport;
    BlockGeometry geometry = BlockGeometry::AsTraced;
    std::unique_ptr<BlockQuadMesh> blockMesh;
    std::unique_ptr<mesh::QuadMesh> smoothedMesh;
    mesh::TMOP::Report tmopReport;
};

#endif // __SKELETON_HXX__
