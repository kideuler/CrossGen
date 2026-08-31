#ifndef __TORSION_TUTTEEMBEDDING_HXX__
#define __TORSION_TUTTEEMBEDDING_HXX__

#include <string>
#include <vector>

#include "mesh/Mesh.hxx"

// Pipeline B, work item 4, step 1 -- docs/cf_flow_pipeline.md Sec. 7.2.
//
//     map dOmega to a convex polygon, boundary lengths proportional to arc length
//     solve L phi = 0 with mean-value weights on the interior
//
// Tutte's theorem: a 3-connected planar graph whose outer face is mapped to a
// convex polygon and whose interior vertices are each a convex combination of
// their neighbours has a straight-line embedding with no crossings. ConeCut
// hands over a disk, the mean-value weights of Floater are strictly positive on
// any triangulation, so the barycentric condition holds and the result is
// bijective -- which is to say **zero flipped triangles, by theorem rather than
// by luck**.
//
// That is the whole reason it is here. The integration of Sec. 6 returns the
// closest map in an L^2 sense and "closest" can invert triangles, generically
// near the cones; the symmetric-Dirichlet barrier of Stage 6 preserves local
// injectivity and cannot restore it, so it needs somewhere injective to start
// from. This is that place. It is a bad map in every other respect -- it knows
// nothing of the field, of the seams, or of the cones, and its distortion is
// enormous -- and none of that matters, because the pass that follows it fits
// the field's target Jacobian under the barrier and only needs the barrier to
// have something to hold on to.
//
// The convex polygon is a circle sized so that the image area matches the
// model's, divided by h^2 to match the scale of the field frame. Getting the
// scale roughly right costs nothing and saves the target-fitting pass from
// spending its first outer steps undoing a global factor -- the symmetric
// Dirichlet energy is not scale invariant, and its minimum is at J = rotation.
class TutteEmbedding {
public:
    struct Options {
        // Radius of the target circle. Non-positive picks one giving the same
        // image area as the model's, scaled by 1/h.
        double radius = -1.0;
        double targetEdge = 1.0;   // h, only read when radius is non-positive
    };

    struct Report {
        int boundaryVertices = 0;
        int interiorVertices = 0;
        int boundaryLoops = 0;      // 1 for a disk
        double radius = 0.0;

        bool solved = false;
        int flippedFaces = 0;       // zero, by the theorem; measured anyway
        double minSignedArea = 0.0;
        double totalArea = 0.0;

        bool valid = false;
        std::vector<std::string> messages;
    };

    TutteEmbedding(const Mesh &omega, const Options &opts);
    explicit TutteEmbedding(const Mesh &omega) : TutteEmbedding(omega, Options()) {}

    const std::vector<Point>& getUV() const { return uv; }
    const Report& getReport() const { return report; }

private:
    // The boundary of Omega as one ordered cycle, or empty if it is not one.
    std::vector<int> boundaryLoop(const Mesh &m) const;

    std::vector<Point> uv;
    Report report;
};

#endif // __TORSION_TUTTEEMBEDDING_HXX__
