#ifndef __RECTANGLE_CERTIFIER_HXX__
#define __RECTANGLE_CERTIFIER_HXX__

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

#include "ATLAS/SquareCarrier.hxx"
#include "ATLAS/SquareTransport.hxx"

// Stage 4 of docs/square_transport_2d_theory_and_implementation.md: generate
// and certify rectangles.
//
// ### The certificate (Sec. 5)
//
// certify() takes a connected set of carrier cells P and decides, exactly,
// whether P is a rectangular tensor grid in its own integer coordinates.
// It develops P from a root cell with G_root = I and G_r = G_q o g_qr across
// every shared edge inside P (Sec. 5.1), then checks Sec. 5.2's list:
//
//   1. every assignment is path-independent -- each cell, reached along any
//      edge of P, gets one transform, and each vertex one coordinate;
//   2. distinct cells occupy distinct unit squares and every square of the
//      bounding rectangle is occupied;
//   3. distinct vertices get distinct coordinates, so the incidences are the
//      grid's and there are no extra identifications;
//   4. (follows from 1-3: the boundary is the rectangle's four-sided loop);
//   5. no feature edge (dS, a material interface) lies inside the block, and
//      no protected vertex sits anywhere but at one of its four corners.
//
// Each of those is a real case, not a formality. A band of cells around a hole
// is locally a perfect grid and fails 1; a band cut once fails 3 (the cut's
// two sides are one carrier edge, developed twice); an L-tromino fails 2. That
// is why Sec. 5.2 insists a cell count or a zero holonomy is not enough.
// Failures come back as a Witness naming the offending cells and vertex, which
// Stage 6 reads as an obstruction to attack (Sec. 8.3).
//
// A certified block brings its map with it (Sec. 5.3): on the square of cell
// q it is X_q o G_q^-1, piecewise bilinear, continuous because shared edges
// agree and injective because the carrier is embedded. Nothing is integrated
// and no isoline is traced.
//
// ### Candidate generation
//
// Sec. 12 asks for bounded neighbourhoods and "promising occupied boxes". The
// promising ones are fixed by Sec. 6's conformity rule, which this code takes
// seriously enough to derive its candidates from it: a vertex some block has
// as a corner must be a corner of every block at it, so every *forced*
// macrovertex (SquareCarrier::forcedMacrovertex) extends block sides straight
// through regular vertices until they meet dS or another forced vertex. The
// cells between those lines -- the base complex of the forced set -- are the
// finest decomposition with no T-junction, and each of its patches is a
// candidate. From there:
//
//   * the same with the designated vertices of Stage 3 added to the forced
//     set (the layout the templates meant);
//   * a patch that fails its certificate (a periodic band round a hole) gets
//     a cut at the witness vertex and the lines are retraced, until every
//     patch certifies;
//   * unions of two adjacent patches that certify, and boxes grown row by
//     row from each patch while the certificate holds.
//
// Singletons are always candidates (Sec. 6) and are implicit.
class RectangleCertifier {
public:
    struct Options {
        int maxCandidates = 200000;
        int maxGrowthSteps = 64;
        bool merges = true;
        bool growth = true;
        // Cut rounds for patches that do not certify.
        int maxCutRounds = 200;
        // Work budget for the merges and the growth, in cells certified, as a
        // multiple of the carrier's cell count. On a fine carrier nearly every
        // merged or grown box is incompatible with the singletons round it,
        // and certifying them all costs far more than it can buy.
        double workFactor = 8.0;
    };

    enum class Failure {
        None,
        Empty,
        Disconnected,
        TransportConflict,   // a cell reached two ways with two transforms: a cycle with holonomy
        VertexPathDependence,// a vertex gets two coordinates: periodic or self-touching
        VertexCollision,     // two vertices at one coordinate: an extra identification
        DuplicateOccupancy,  // two cells on one square
        MissingSquare,       // the occupied squares do not fill the bounding rectangle
        FeatureInside,       // an interface edge inside the block
        ProtectedSwallowed,  // a protected vertex inside the block or a side
        Overflow
    };

    struct Witness {
        Failure kind = Failure::None;
        int cellA = -1, cellB = -1;
        int vertex = -1;
    };

    struct Certificate {
        std::vector<int> cells;              // row-major: square (i, j) at i + nu * j
        std::vector<SquareTransport> G;      // per cell: its coordinates -> block's
        int nu = 0, nv = 0;
        std::array<int, 4> corners{{-1, -1, -1, -1}};   // (0,0) (nu,0) (nu,nv) (0,nv)
        std::array<std::vector<int>, 4> sides;          // vertex chains, counter-clockwise
        int material = 0;
        double distortion = 0.0;             // D(P): mean conformal distortion of its pieces
        double complexity = 0.0;             // C(P): side turning, in quarter turns
        std::string source;
    };

    struct Report {
        int cells = 0;
        int forcedVertices = 0, softVertices = 0, cutVertices = 0;
        int basePatchesHard = 0, basePatchesAll = 0;
        int candidates = 0;                  // not counting singletons
        int fromBaseHard = 0, fromBaseAll = 0, fromMerges = 0, fromGrowth = 0;
        int certifyCalls = 0;
        long long certifiedCells = 0;
        bool workBudgetHit = false;
        int failures = 0;
        std::array<int, 11> failuresByKind{};
        int largestCandidate = 0;
        bool cutsConverged = true;
    };

    RectangleCertifier(const SquareCarrier &carrier, const Options &opts);

    // The certificate itself, usable on any cell set.
    bool certify(const std::vector<int> &cells, Certificate &cert, Witness &w) const;

    const std::vector<Certificate> &getCandidates() const { return candidates_; }
    // Base-complex covers: exact partitions into certified patches, when they
    // exist (with the designated vertices forced, and without).
    const std::vector<int> &baseCoverAll() const { return baseAll_; }
    const std::vector<int> &baseCoverHard() const { return baseHard_; }
    const std::vector<Witness> &getWitnesses() const { return witnesses_; }
    const Report &getReport() const { return report_; }
    const SquareCarrier &getCarrier() const { return C_; }

    static const char *failureName(Failure f);
    // The number of patches of C's base complex -- lines traced straight
    // through regular vertices from every forced macrovertex and every
    // designated one -- without certifying any of them. It is the block
    // count Stages 4 and 5 start from, cheap enough for Stage 6 to evaluate
    // after every trial move on a coarse carrier (N_B of Sec. 8.3's score).
    //
    // A patch need not be a disk: on a domain with holes, a carrier with no
    // singular vertex round a hole has a band there that no line cuts, which
    // counts as one patch and is no block at all -- single contact needs three
    // to go round a hole (Sec. 6). When `holeDeficit` is given it receives
    // sum over patches of (1 - chi), chi = V - E + F of the patch, so that a
    // caller can charge for it (CavityRewrite::AnnealOptions::patchTopology).
    static int basePatchCount(const SquareCarrier &C, int *holeDeficit = nullptr);
    // Distortion and complexity of a certified block.
    static void measure(const SquareCarrier &C, Certificate &cert);

private:
    // Trace the base complex of `forced`; returns the patches as cell lists.
    std::vector<std::vector<int>> baseComplex(const std::vector<char> &forced) const;
    // Base complex with cuts until every patch certifies. Appends the
    // certified patches' candidate indices to `cover` (empty on failure).
    void certifiedBase(std::vector<char> forced, std::vector<int> &cover, const char *source,
                       int &patchCount);
    int addCandidate(Certificate &&cert);

    const SquareCarrier &C_;
    Options opts_;
    std::vector<Certificate> candidates_;
    std::unordered_map<unsigned long long, int> seen_;
    std::vector<int> baseAll_, baseHard_;
    std::vector<Witness> witnesses_;
    mutable Report report_;
    // Scratch for certify(): cell -> local index, stamped.
    mutable std::vector<int> stamp_, local_;
    mutable int stampValue_ = 0;
};

#endif // __RECTANGLE_CERTIFIER_HXX__
