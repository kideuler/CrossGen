#ifndef __SQUADGEN_DOMAIN_HXX__
#define __SQUADGEN_DOMAIN_HXX__

#include <array>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include "mesh/Mesh.hxx"

namespace squadgen {

// A seeded random source that gives the same numbers on every platform: std::mt19937_64 for the bits, and the
// deviates made from them here rather than by <random>'s distributions, whose algorithms the standard leaves to each
// library. Not numpy's or torch's generator, so a seed does not give the token sets or the noise squadgen2d draws for
// it; the reference case of export_onnx carries those for a comparison that needs them.
class Random {
public:
    explicit Random(uint64_t seed) : engine_(seed) {}
    double uniform();          // [0, 1), 53 bits
    double normal();           // Box-Muller; the second deviate of a pair is the next call's
    size_t below(size_t n);    // uniform in [0, n)

private:
    std::mt19937_64 engine_;
    bool haveSpare_ = false;
    double spare_ = 0.0;
};

// The networks' frame (Sec. 4 step 1, squadgen2d geometry.normalisation): the area centroid of Omega at the origin and
// its farthest boundary vertex at radius 0.999. The area centroid rather than the vertex mean, so that the frame does
// not depend on how the domain was meshed. A point of the model maps to (x - centre) * scale.
struct Frame {
    Point centre{0.0, 0.0};
    double scale = 1.0;
    Point toFrame(const Point &p) const { return (p - centre) * scale; }
    Point fromFrame(const Point &q) const { return q / scale + centre; }
};

// A model as squadgen2d reads one (geometry.read_obj): the vertices and triangles its OBJ lists, materials ignored --
// the networks see only the domain. Unlike Mesh's loader, which turns a clockwise face (a, b, c) into (a, c, b), a
// clockwise face is reversed to (c, b, a) here, as roundtrip.extraction_mesh reverses it, so that the extraction mesh
// built on it is squadgen2d's triangle for triangle; 4 % of the corpus' faces are clockwise.
struct Model {
    std::vector<Point> vertices;
    std::vector<Triangle> triangles;
};
Model readOBJ(const std::string &path);
// A Mesh's arrays, as they stand (counter-clockwise already, in the loader's way).
Model modelOf(const Mesh &mesh);

// The domain in the frame, and what Sec. 4 steps 1-3 draw on it.
class Domain {
public:
    explicit Domain(const Model &model);

    const Frame &frame() const { return frame_; }
    // dS as segments in the frame: the edges one triangle uses.
    const std::vector<std::array<Point, 2>> &segments() const { return segments_; }
    double area() const { return cdf_.empty() ? 0.0 : cdf_.back(); }   // in the frame

    // n uniform samples of Omega, in the frame (geometry.sample_area).
    std::vector<Point> sampleArea(size_t n, Random &rng) const;
    // A token set of n points: the farthest-point sampling of max(8n, 8192) area samples (geometry.token_sets).
    std::vector<Point> tokenSet(int n, Random &rng) const;
    // The truncated SDF at points inside Omega, clip(-d / tau, -1, 0), d the distance to dS (geometry.truncated_sdf
    // with known_inside): what the Geom-AE reads at every token.
    std::vector<float> insideSDF(const std::vector<Point> &points, double tau) const;

private:
    Frame frame_;
    std::vector<Point> vertices_;        // in the frame
    std::vector<Triangle> triangles_;
    std::vector<double> cdf_;            // cumulative triangle area, in the frame
    std::vector<std::array<Point, 2>> segments_;
};

// Greedy max-min selection of n of the candidates (geometry.farthest_point_sampling): a random first one, then each
// time the one farthest from those chosen so far. Returns their indices.
std::vector<int> farthestPoints(const std::vector<Point> &candidates, int n, Random &rng);

// The mesh the fields are decoded on and QUARRY extracts from (roundtrip.extraction_mesh): the model without unused
// vertices, split 1-to-4 at its edge midpoints until the median length of its triangles' first edges is at most
// 1.5 hRel times its bounding-box diagonal, four times at most. Vertices and triangles in squadgen2d's order: the
// model's used vertices, then each round's midpoints in the order of their edges' (lower, higher) vertex ids, and the
// four children of every triangle in four runs. `levels`, when given, is how many rounds it took.
Model extractionMesh(const Model &model, double hRel, int *levels = nullptr);

// Every triangle's centroid in the frame, as float32 pairs: the points the decoder is queried at.
std::vector<float> centroidsInFrame(const Model &m, const Frame &frame);

// Points as the float32 pairs the networks read.
std::vector<float> flatten(const std::vector<Point> &points);

}  // namespace squadgen

#endif // __SQUADGEN_DOMAIN_HXX__
