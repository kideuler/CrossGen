#ifndef __CHART_DISTANCE_FIELD_HXX__
#define __CHART_DISTANCE_FIELD_HXX__

#include <array>
#include <map>
#include <string>
#include <vector>

#include "mesh/BlockDecomposition.hxx"
#include "mesh/Mesh.hxx"

class BlockQuadMesh;
class NpzWriter;

// SQuadGen's Chart Distance Fields (CDF and DCDF, the SQuadGen paper's Sec. 5)
// of a BlockDecomposition, in the data layout the release's CDFGen writes and
// its Python loader reads (docs/squadgen_2d.md Secs. 2.2 and 2.4), so that a
// layout made by any of our methods becomes an SQ-VAE training example.
//
// ### What the fields are
//
// Every chart (here: every block) is cut by two mid-lines into four subcharts,
// one per chart corner. A subchart carries two coordinates c0, c1 in [0, 1],
// each 1 on the chart side through its corner and 0 on the mid-line facing
// that side, so that
//
//     CDF  = 1 - max(c0, c1)   0 on every chart side, 1 at the chart centre,
//     DCDF = min(c0, c1)       0 on every mid-line, 1 at the chart corner,
//
// the second being the first read on the dual layout, whose charts are the
// four subcharts around each corner. Per point they come with unit gradients
// and with the offsets from the point to its subchart's corner (the dual-chart
// centre), its two side mid-points, the chart centre and the chart's other
// three corners: what the release's loader computes as `dcdf, gdcdf, cdf,
// gcdf, offset, offsetb, offsetc, offsetd, offset1-3`.
//
// ### How CDFGen builds them, and what is different here
//
// CDFGen (QuadTools, patchsample.cpp 683-1372) starts from an unstructured
// quad mesh and first has to *find* the charts: singular vertices by valence,
// separatrices traced from them, loops split, charts by flood fill. That search
// is where it fails on planar input -- a domain corner at valence 3 is not a
// singular vertex, and without docs/squadgen_2d.md's Appendix A patch an
// O-grid round a square hole comes out as two charts instead of four. Here
// the charts are given: they are the blocks, and the mesh is the transfinite
// grid BlockQuadMesh lays on each, so a block corner is a chart corner by
// construction and no valence test is made at all.
//
// The rest is CDFGen's algorithm, read in block-grid coordinates where it
// walks half-edges:
//
//   1. Mid-line split (CDFGen 741-993). A chord of BlockQuadMesh -- the macro
//      edges that must share an interval count, CDFGen's "arc group" -- gets
//      one split: the normalised cumulative-length profiles of its macro
//      edges, oriented consistently along the chord, are averaged, and the
//      split sits where that average crosses 1/2, at a fraction of one grid
//      edge kept inside [0.01, 0.99]. Every block the chord crosses is cut
//      along the grid line at that fractional index, which is what CDFGen's
//      quad-by-quad propagation of "the same t on the opposite edge"
//      produces. Two cuts per block; a cell both cross becomes four.
//      [Departure] CDFGen visits each arc from the charts on both sides and
//      skips the second visit unless the arc starts on the boundary, so an
//      interior arc that starts on dS is averaged in twice; here every macro
//      edge counts once.
//   2. Ring-averaged lengths (CDFGen 1017-1070). Every edge of the split mesh
//      is given the mean length of its edge ring -- the column of edges a quad
//      strip crosses, end to end through every block it enters -- the paper's
//      "edge length assignment". It is what makes the two opposite sides of a
//      chart, which share their rings, carry the same coordinate profile, and
//      so makes CDF continuous across a mid-line.
//   3. Coordinates (CDFGen 1072-1185). A chart corner is 1 and a side's split
//      vertex 0; in between, along each half side, the value is the remaining
//      ring length to the split over the whole half side's. Each half-side
//      edge's (start, end) pair is then carried across its quad strip to the
//      mid-line: every half-edge of a subchart refers to the half-side edge of
//      its own column or row, with a sign for direction -- exactly the
//      table CDFGen floods, here indexed rather than flooded.
//   4. Evaluation (CDFGen 1536; Appendix B of the guide). Each quad is split
//      into two triangles along the diagonal pointing at its subchart corner;
//      in a triangle, the coordinate along each of its two quad edges is the
//      oblique projection parallel to the other, mapped linearly onto that
//      edge's pair.
//
// ### Samples
//
// The release stores 50 000 Poisson-disk samples and, per token count n, five
// Poisson-subsampled sets, all as (triangle, barycentric) pairs on the split
// mesh triangulated as in 4. Here the samples are uniform by area and the sets
// are farthest-point samples of max(8n, 8192) area samples, the token sets
// py/squadgen/squadgen2d/geom_ae.py trains the 2D Geom-AE on, so that a token of the
// SQ-VAE sits where the Geom-AE expects one (guide Sec. 3.2, token alignment).
// The draw is seeded and deterministic.
//
// ### Frame
//
// Everything is in the model's own coordinates; `invT` is written as the
// identity. The release's loader re-centres and rescales on load, and does so
// whatever CDFGen did first, so normalising here would only be undone.
class ChartDistanceField {
public:
    struct Options {
        // Target edge length of the block mesh the fields live on. 0 (the
        // default) takes relativeEdgeLength times the diagonal of the
        // decomposition's bounding box, so that models of any size get the same
        // resolution relative to themselves.
        double targetEdgeLength = 0.0;
        double relativeEdgeLength = 0.01;
        // BlockQuadMesh's floor on the edges per chord, and its Winslow passes
        // over the blocks the transfinite blend folds.
        int minIntervals = 1;
        int smoothingPasses = 500;

        // Uniform area samples (sample_bary_coords / sample_point_tri_id).
        int samples = 20000;
        // Token sets of fpsMin, 2 fpsMin, ... up to fpsMax points, fpsCopies of
        // each, by farthest-point sampling of max(fpsOversample * n,
        // fpsMinCandidates) area samples. fpsMax 0 writes none.
        int fpsMin = 512;
        int fpsMax = 4096;
        int fpsCopies = 5;
        int fpsOversample = 8;
        int fpsMinCandidates = 8192;
        unsigned seed = 0;
    };

    struct Report {
        double target = 0.0;            // the block mesh's target edge length
        int blocks = 0;                 // charts (num_complex)
        int chords = 0;
        int quads = 0;                  // the block mesh's
        int vertices = 0;               // the split mesh's
        int faces = 0;
        int subcharts = 0;
        int patches = 0;                // dual charts: one per block corner vertex
        int coordinatePairs = 0;        // edge_color_store, the reserved 0 included
        int boundaryLoops = 0;          // num_boundary
        int invertedFaces = 0;          // split-mesh quads with a triangle of area <= 0
        double minTriangleArea = 0.0;
        double area = 0.0;
        // Consistency checks that must read 0 on any decomposition: a split
        // vertex placed twice at two places, a side vertex given two
        // coordinates, a half-edge left with no coordinate pair.
        int splitMismatches = 0;
        int coordinateMismatches = 0;
        int missingPairs = 0;
        int samples = 0;
        int tokenSets = 0;              // fps sizes written, each fpsCopies times
        double seconds = 0.0;
        bool valid = false;             // built, and every check above at 0
        std::vector<std::string> messages;
    };

    // The fields at one point, as the release's loader computes them.
    // `offsets` are corner (offset), side mid-point (offsetb), chart centre
    // (offsetc), side mid-point (offsetd), then the chart's other three
    // corners (offset1-3), each as the vector from the point.
    struct Fields {
        double c0 = 0.0, c1 = 0.0;
        double dcdf = 0.0, cdf = 0.0;
        Point gdcdf{0.0, 0.0}, gcdf{0.0, 0.0};
        std::array<Point, 7> offsets{};
    };

    // Throws std::invalid_argument for a decomposition that does not cover
    // the model (BlockDecomposition::covers()) or that BlockQuadMesh cannot
    // mesh whole, and std::runtime_error if the split mesh is not a manifold
    // or a chord cannot be oriented. Inverted elements are reported, not
    // thrown: whether a folded example is still worth training on is the
    // caller's choice (Report::invertedFaces).
    ChartDistanceField(const BlockDecomposition &decomposition, const Options &opts);

    // ── The split mesh (subdiv_mesh) ────────────────────────────────────────
    // Its first quadVertexCount() vertices are the block mesh's, in its order,
    // so quadFaces() indexes vertices() too. Faces are counter-clockwise.
    const std::vector<Point> &vertices() const { return verts_; }
    const std::vector<std::array<int, 4>> &faces() const { return faces_; }
    int quadVertexCount() const { return quadVertexCount_; }
    const std::vector<std::array<int, 4>> &quadFaces() const { return quadFaces_; }
    // Per face: 1 when the diagonal toward the subchart corner is 1-3, 0 when
    // it is 0-2 (quad_split). Triangle 2f + k of face f is then (0, 1, 2),
    // (2, 3, 0) or (3, 0, 1), (1, 2, 3).
    const std::vector<signed char> &split() const { return split_; }
    std::array<int, 3> triangle(int t) const;
    int triangleCount() const { return 2 * static_cast<int>(faces_.size()); }

    // ── Charts ──────────────────────────────────────────────────────────────
    const std::vector<int> &faceSubchart() const { return faceSubchart_; }       // face2quad
    const std::vector<int> &subchartPatch() const { return subchartPatch_; }     // quad2patch
    const std::vector<int> &subchartBlock() const { return subchartBlock_; }
    const std::vector<int> &patchVertex() const { return patchVertex_; }         // corner of each patch
    const std::vector<std::array<int, 4>> &subchartABCD() const { return abcd_; }  // offset_abcd_id
    const std::vector<std::array<int, 3>> &subchartOthers() const { return others_; }  // offset_123_id

    // ── Coordinates ─────────────────────────────────────────────────────────
    // Per half-edge 4f + k (face f, from corner k to k + 1), a signed index
    // into coordinatePairs(): +i reads pair i as (tail, head), -i reversed
    // (edge_vertex_color, edge_color_store). Pair 0 is (0, 0) and unused.
    const std::vector<int> &halfEdgePair() const { return halfEdgePair_; }
    const std::vector<std::array<double, 2>> &coordinatePairs() const { return pairs_; }

    // ── Samples ─────────────────────────────────────────────────────────────
    const std::vector<int> &sampleTriangles() const { return sampleTri_; }
    const std::vector<std::array<double, 3>> &sampleBary() const { return sampleBary_; }
    struct TokenSets {
        std::vector<std::vector<int>> triangles;                 // [copy][i]
        std::vector<std::vector<std::array<double, 3>>> bary;    // [copy][i]
    };
    const std::map<int, TokenSets> &tokenSets() const { return fps_; }

    Point pointAt(int triangle, const std::array<double, 3> &bary) const;
    Fields evaluate(int triangle, const Point &p) const;

    const Report &getReport() const { return report_; }
    const Options &getOptions() const { return opts_; }

    // Every array of the release's NPZ (guide Sec. 2.4) under its key and in
    // its dtype, so that data_tools/test_load_new_format.py loads the file
    // unchanged; plus subchart_block, patch_vertex and target_edge_length,
    // which it does not read. Metadata the caller wants beside them goes into
    // `out` before or after; write() is the caller's.
    void addTo(NpzWriter &out) const;

private:
    void buildSplitMesh(const BlockQuadMesh &grid, const BlockDecomposition &D);
    void buildCoordinates();
    void sample();

    Options opts_;
    Report report_;

    int quadVertexCount_ = 0;
    std::vector<std::array<int, 4>> quadFaces_;
    std::vector<Point> verts_;
    std::vector<std::array<int, 4>> faces_;
    std::vector<signed char> split_;
    std::vector<int> faceSubchart_, subchartPatch_, subchartBlock_, patchVertex_;
    std::vector<std::array<int, 4>> abcd_;
    std::vector<std::array<int, 3>> others_;
    std::vector<int> halfEdgePair_;
    std::vector<std::array<double, 2>> pairs_;
    std::vector<int> sampleTri_;
    std::vector<std::array<double, 3>> sampleBary_;
    std::map<int, TokenSets> fps_;

    // Per block of the split mesh, what buildCoordinates() reads: the refined
    // grid of vertex indices, (ns + 2) x (nt + 2) with the two mid-lines in
    // it, the split columns, and whether the block runs counter-clockwise.
    struct Grid {
        int block = -1;
        int ns = 0, nt = 0;             // the block mesh's cells
        int ks = 0, kt = 0;             // the mid-lines sit in cells ks, kt
        bool ccw = true;
        std::vector<int> v;             // (ns + 2) x (nt + 2), row-major in t
        std::vector<int> faces;         // (ns + 1) x (nt + 1) split-mesh faces, row-major in t
        int at(int r, int q) const { return v[static_cast<size_t>(q) * (ns + 2) + r]; }
    };
    std::vector<Grid> grids_;
};

#endif // __CHART_DISTANCE_FIELD_HXX__
