#ifndef __QUARRY_HXX__
#define __QUARRY_HXX__

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "QUARRY/Basins.hxx"
#include "mesh/BlockDecomposition.hxx"
#include "mesh/BoundaryFeatures.hxx"
#include "mesh/Mesh.hxx"

// QUARRY: the quad block decomposition SQuadGen's chart distance fields describe, cut out of the model they were
// decoded on. It is the CrossGen counterpart of the release's QuadExtraction (QuadTools, cdfdcdf2quadmesh.cpp;
// docs/squadgen_2d.md Secs. 4 step 8 and 5), the last step of SQuadGen's inference, and like it reads the four
// fields the release hands its extractor -- CDF, DCDF, and the offsets to the chart centre (offsetc) and to the
// subchart's corner (offset) -- at the centroid of every triangle of a fine triangulation of the domain.
//
//   Stage 1-2  Basins (charts)   the charts, from the CDF and offsetc
//   Stage 1-2  Basins (duals)    the dual charts, one per layout vertex, from the DCDF and offset
//   Stage 3    Junctions, arcs   the mesh vertices where three or more of {charts, dS} meet, and the runs of chart
//                                boundary between them; every chart's boundary a cycle of arcs
//   Stage 4    Corners           a chart's corners are the junctions on its boundary; an arc is collapsed where a
//                                vertex was blurred into two, and a chart short of four corners gets the rest on
//                                its stretch of dS, where a dual chart predicts one, else at a corner of the model,
//                                else at a peak of the DCDF
//   Stage 5    Repair            a chart that still cannot be a quad and has no corner left to take is dissolved
//                                into its neighbour, and Stages 3-4 run again
//   Stage 6    Blocks            the shared BlockDecomposition: a side on dS is dS, a side between two charts their
//                                common boundary, smoothed; BlockQuadMesh meshes it
//
// ### What is kept from the release, and what is not
//
// The idea is the release's: a chart is a region of high CDF, grown outwards from its peak in the order of how far
// the centre a triangle predicts lies from the region's own, and the blocks are the regions' outlines. What is not
// kept is everything that made it fragile on decoded fields, which the 2026-10-29 round trip measured (memory
// squadgen-roundtrip-2d): true fields recovered 89 % of held-out layouts and decoded fields 49 %, failing by
// over-segmentation, and 52 of 143 meshes extracted even from the true fields left the domain.
//
//   * Seeds. The release seeds a chart at every strict CDF maximum over 8 vertex rings and never asks whether two
//     seeds predict one chart. Here the peaks go through the persistence clustering of Basins.hxx, where two peaks
//     whose predicted centres agree are one chart however deep the dip between them.
//   * Corners. The release makes polygons from the cluster outlines and then collapses edges in OpenMesh, ordered
//     by edge length and valence, until they come out as quads -- or do not, and a non-quad polygon is left in the
//     mesh. Here the corner count is decided per chart before anything is built: a conforming quad layout has no
//     T-junction, so every junction on a chart's outline is one of its corners, and the only two ways the outline
//     of a correct chart can still be wrong are fixed by name. A vertex of valence four or more that the blur at
//     the walls split into two junctions joined by a short arc gives the two charts along that arc five corners
//     each, and that arc -- and only an arc whose every chart has a corner too many -- is collapsed. A chart whose
//     corner is on dS where no other chart meets it has no junction there, and gets it from the dual charts: the
//     DCDF peaks at every layout vertex, and its offset names the vertex from every triangle of its subcharts.
//   * Geometry. The release takes each patch side as the straight chord between its two corners and pulls the
//     subdivided patch back onto the triangles by projection; on a curved dS the chord leaves the domain. Here a
//     side on dS *is* dS between its two corners, and a side between two charts is the boundary the two regions
//     share, so no block can leave the domain or cross a hole; the mesh is BlockQuadMesh's on the result.
//
// Material interfaces are not read: SQuadGen's fields know none, the charts are the CDF's, and a block takes the
// material of its peak triangle.
//
// What QUARRY cannot do is make a layout out of fields that do not describe one. A chart with more than four
// corners after Stage 4, or with too few and nothing to take them from, is reported (ChartStatus) and left out of
// the blocks, and the decomposition then does not cover the model, which is how every other method here says the
// same thing.
class QUARRY {
public:
    // Per triangle of the model, sampled at its centroid, in the model's units: what infer_sqdiffuse.py writes to
    // extract_mesh.npz, less the release's frame (the offsets here are the model's, not the network's).
    struct Fields {
        std::vector<double> cdf, dcdf;
        std::vector<Point> offsetCentre;   // offsetc: centroid -> the chart's centre
        std::vector<Point> offsetCorner;   // offset:  centroid -> the subchart's chart corner
    };

    struct Options {
        quarry::Basins::Options charts;
        quarry::Basins::Options duals;
        // Stage 4. An arc is collapsed only when it is shorter than this fraction of the square root of the area of
        // the smaller chart beside it -- and only when every chart along it has a corner too many, which is what
        // says it is one vertex blurred into two rather than a short side.
        double collapseLength = 1.5;
        // An arc of at most this many mesh edges is collapsed whatever the corner counts say: two junctions that
        // close are one vertex as far as the extraction mesh can tell, and the true layout's vertices of valence
        // five or more come out of the region growing as exactly that. A chart it leaves one corner short then
        // takes the missing one on dS.
        int unresolvedEdges = 2;
        // A corner taken from a dual chart: the dual chart covers at least this fraction of the chart's area (a
        // true corner's subchart holds about a quarter), and the vertex it predicts lies within dualReach of the
        // chart's stretch of dS, both in units of the square root of the chart's area.
        double cornerShare = 0.03;
        double dualReach = 0.15;
        // A corner is added only clear of the chart's other corners: this many mesh edges away along its outline, or
        // cornerSpacing x sqrt(chart area) away in the plane.
        int cornerGap = 3;
        double cornerSpacing = 0.1;
        // The last place a chart short of corners looks on dS: a local maximum of the DCDF along its stretch of dS
        // at least this high.
        double dcdfCorner = 0.5;
        // Stage 6: a layout vertex on dS this many mesh edges or fewer from a corner of the model turning by
        // cornerSnapAngle or more, along dS with no other layout vertex between, is moved onto it.
        int cornerSnapEdges = 4;
        double cornerSnapAngle = 40.0;   // degrees of turning, the least a corner snapped to may have
        // Stage 5: at most this many charts that cannot be quads are dissolved, one per round.
        int repairRounds = 32;
        // Passes of (1/4, 1/2, 1/4) smoothing on every side between two charts, end points held: the shared
        // boundary of two regions runs along triangle edges, a staircase on the scale of one triangle.
        int sidePasses = 8;
        // A side between two charts that ends on dS is cut back, at that end, to where it is this many mean outline
        // edges from dS, and joined straight to its corner from there.
        double endClearance = 2.0;
    };

    // How each chart (Basins region of the CDF) came out.
    struct ChartStatus {
        int junctions = 0;        // on its outline, as Stage 3 found them
        int corners = 0;          // after Stage 4
        int added = 0;            // of those, taken on dS from the dual charts or the model's corners
        int cycles = 0;           // closed boundary curves: 1 for a disk
        bool dissolved = false;   // given to a neighbour in Stage 5
        bool block = false;       // made a block
        std::string problem;      // empty when it did
    };

    struct Report {
        quarry::Basins::Report charts, duals;
        int junctions = 0;
        int arcs = 0;
        int collapsedArcs = 0;
        int addedCorners = 0;
        int dissolvedCharts = 0;
        int snappedCorners = 0;
        int repairRounds = 0;
        int blocks = 0;
        int failedCharts = 0;     // neither a block nor dissolved
        int unmatchedSides = 0;   // between a block and a chart that is none
        int macroVertices = 0;
        int macroEdges = 0;
        int boundaryVertices = 0;
        double seconds = 0.0;
        bool covers = false;      // BlockDecomposition::covers()
        std::vector<std::string> messages;
    };

    QUARRY(std::shared_ptr<const Mesh> mesh, Fields fields, const Options &opts);

    void run();

    const BlockDecomposition &decomposition() const { return blocks_; }
    const quarry::Basins &charts() const { return *charts_; }
    const quarry::Basins &duals() const { return *duals_; }
    // The chart of every triangle after Stage 5 (Basins' labels with the dissolved charts given away).
    const std::vector<int> &chartLabels() const { return labels_; }
    const std::vector<ChartStatus> &chartStatus() const { return status_; }
    const Report &getReport() const { return report_; }
    const Options &getOptions() const { return opts_; }
    const Mesh &mesh() const { return *mesh_; }

private:
    // One closed boundary curve of a chart, as half-edges of its own triangles with the chart on their left.
    struct Cycle {
        std::vector<int> vertex;      // mesh vertex at the start of half-edge i
        std::vector<int> across;      // chart across half-edge i; -1 on dS, -2 where no chart reached
        std::vector<int> triangle;    // the chart's own triangle of half-edge i
        double signedArea = 0.0;
    };
    // A run of one chart's outline between two consecutive junctions, all of it against one thing: another chart
    // or dS. Cycle vertices begin..end, cyclically; a closed arc (an outline with no junction at all) has
    // begin == end and u == v == -1.
    struct Arc {
        int chart = -1, across = -1;
        int begin = 0, end = 0;
        int u = -1, v = -1;
        int edges = 0;                // mesh edges along it
        double length = 0.0;
        int twin = -1;                // the same arc in the outline of `across`
    };

    std::vector<std::vector<Cycle>> boundaryCycles() const;
    // Stages 3, 4 and 6 on the current labels. Returns the charts Stage 5 should dissolve.
    std::vector<int> assemble();
    int group(int v);

    std::shared_ptr<const Mesh> mesh_;
    Fields fields_;
    Options opts_;
    Report report_;
    std::unique_ptr<BoundaryFeatures> features_;
    // The points each triangle predicts, which Basins reads by reference.
    std::vector<Point> centreTargets_, cornerTargets_;
    std::unique_ptr<quarry::Basins> charts_, duals_;
    std::vector<int> labels_;
    std::vector<int> groupParent_;    // Stage 4's union-find over mesh vertices: junctions one vertex
    std::vector<ChartStatus> status_;
    BlockDecomposition blocks_;
};

#endif // __QUARRY_HXX__
