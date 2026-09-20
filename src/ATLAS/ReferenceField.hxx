#ifndef __REFERENCE_FIELD_HXX__
#define __REFERENCE_FIELD_HXX__

#include <array>
#include <complex>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "ATLAS/SquareCarrier.hxx"
#include "dualmbo/DualMBO.hxx"
#include "mesh/Mesh.hxx"

// A reference cross field for ATLAS: DualMBO on the input mesh, as a prior the
// search's layouts are scored against, never as a field to integrate
// (docs/atlas_crossfield_guidance.md, Sec. 3).
//
// ### Why a prior, and why this field
//
// Nothing in Stages 3-6 is blind to topology -- every state is a certified
// carrier -- but everything is blind to *direction*: the base-complex count
// N_B, the valence defect and the cell shapes score a layout whose blocks run
// across a rectangle on the diagonal exactly as they score one that runs along
// its walls. The annealer has carried a term for that since v2, against a
// harmonic interpolation of the boundary tangent on the coarse triangulation
// (CoarseDomain::crossAngle), and measured over the corpus it does nothing
// distinguishable from search noise (guidance note, Sec. 2.3). Part of the
// reason is the reference itself: on a disk the harmonic extension of
// exp(4i theta) = z^4 is z^4, one index-1 singularity at the centre with
// weight r^4, so over half the disk -- exactly where an O-grid's core and
// orientation are decided -- it has no direction at all. DualMBO resolves the
// same disk into four +1/4 cones on the diagonals (DualMBO::setPinDiskCenters
// fixes their angle) and has a direction everywhere but at four points.
//
// The solve is TORSION's Stage 0, line for line (TORSION::runField): the
// orthogonal two-point penalty weight, gamma = 10, and the tau-continuation at
// ratio 0.25 down to a 20-edge floor, the interfaces as aligned edges. That
// loop is copied here rather than factored out of TORSION so the paper's
// pipelines stay byte-identical; the guidance note's Sec. 9 records it as the
// thing to factor into src/dualmbo later. It costs 0.2 s median, 1.7 s worst
// over the corpus, against ATLAS's 3-13 s.
//
// ### What is read off it
//
//   * the representation vector u = exp(4i phi) as a lookup grid, so a cell
//     of any carrier -- the fine one, a coarse one, a template's -- can be
//     averaged over: u_bar_q = (1/|q|) int_q u dA, with the coherence
//     rho_q = |u_bar_q| in [0, 1] saying how much direction the field has at
//     the cell's own scale (near 1 where it is smooth across q, near 0 on a
//     cell that straddles a cone, where u turns by 2 pi);
//   * its cones, the interior vertices of index +-1/4, with opposite-signed
//     pairs in one material region cancelled first as MERIDIAN's Stage 1 does
//     (ConeSingularities::cancelDipoles): a +1/-1 pair of the smoothest field
//     leaves every region's index count where it was, and asking the carrier
//     to reproduce it would be asking for a dislocation.
//
// ### The two terms (guidance note, Secs. 3.2 and 3.3)
//
// E_dir, the direction term, is sum_q m_q / |Omega| with
//
//     m_q = |q| (rho_q - Re(conj(c_q) u_bar_q)) / 2 = |q| rho_q (1 - cos 4(theta_q - phi_q)) / 2,
//
// c_q = exp(4i theta_q) the cell's own cross (CavityRewrite::energy's
// definition: its two mean edge directions, the second turned back a quarter).
// It is area-normalised, so the coarse carrier, its realisation and searches
// at different spacings are on one scale: 0.01 is the whole domain ~3 degrees
// off, 0.1 ~9 degrees RMS.
//
// E_sing, the singularity term, is a minimum-cost matching of the carrier's
// singular vertices against the cones, same sign only. A direction score is a
// bulk measure and misses the failure the field exposes most clearly: ATLAS
// pushes the field's interior cones onto dS, and a +1/4 there -- a block
// corner on a smooth arc -- leaves one cell spanning pi (geom010/016/022/025,
// TFI worst SJ 0.09-0.19), while E_dir of geom016's mesh is lower than that
// of geom003's textbook O-grid. So every boundary defect pays a flat cost by
// sign whatever it is matched to, and only an interior singularity near its
// cone is cheap. Weights are in blocks, the annealer's unit.
class ReferenceField {
public:
    struct Options {
        // DualMBO as TORSION's Stage 0 runs it (TORSION::Options' defaults).
        double gamma = 10.0;
        int maxSteps = 500;
        DualMBO::PenaltyWeight weight = DualMBO::PenaltyWeight::Orthogonal;
        bool tauContinuation = true;
        double tauRatio = 0.25;
        double tauFloorEdges = 20.0;
        int tauLevelSteps = 2000;
        bool pinDiskCenters = true;
        bool alignToInterfaces = true;
        // Cancel opposite-signed cone pairs sharing a material region,
        // closest first, as MERIDIAN's Stage 1 does. dipoleRadius > 0 limits
        // it to pairs closer than that fraction of the diagonal.
        bool cancelDipoles = true;
        double dipoleRadius = 0.0;
        // Lookup grid step, as a fraction of the input's mean edge length,
        // and a cap on its size.
        double gridStep = 0.5;
        int maxGridNodes = 4000000;
    };

    // Weights of E_sing, in blocks (guidance note, Sec. 3.4).
    struct SingularityWeights {
        double wPos = 0.5;        // interior singularity at d from its cone: wPos min(1, d / r0)
        double wExtra = 1.0;      // interior singularity with no cone
        double wMissing = 1.0;    // cone with no singularity at all
        double wEdgePlus = 1.5;   // +1/4 on dS: a block corner on a smooth run, one cell spanning pi
        double wEdgeMinus = 0.25; // -1/4 on dS: three cells where two belong, 60 degrees each
        double r0 = 0.1;          // the layout's length scale, fraction of the diagonal
    };

    // One index unit: a cone of the field, or a singular vertex of a carrier
    // (a vertex of defect k is k units at one point).
    struct Singularity {
        Point x{0.0, 0.0};
        int sign = 0;             // +1: index +1/4 (too few cells), -1: -1/4
        bool boundary = false;    // on dS (a carrier's boundary defect)
        int vertex = -1;          // mesh vertex (cone) or carrier vertex
        int region = -1;          // cones: material component, -1 astride an interface
    };

    struct SingularityReport {
        double energy = 0.0;
        int cones = 0;
        int interiorPlus = 0, interiorMinus = 0;
        int boundaryPlus = 0, boundaryMinus = 0;
        int matched = 0;          // interior singularities matched to a cone
        int absorbed = 0;         // cones taken by a boundary defect
        int missing = 0;          // cones left over
        int extra = 0;            // interior singularities left over
        int dipoles = 0;          // opposite pairs of the carrier cancelled within r0
        double meanDistance = 0.0;  // matched pairs, in units of r0
        // Where the carrier and the field disagree: the carrier vertices of
        // boundary defects, of unmatched interior singularities and of
        // matched ones more than r0/2 from their cone; and the cones no
        // interior singularity sits within r0/2 of. What directed moves aim at.
        std::vector<int> looseVertices;
        std::vector<Point> looseCones;
    };

    struct Report {
        bool built = false;
        std::string reason;
        int triangles = 0;
        int levels = 0, steps = 0;
        bool converged = false;
        int rawPlus = 0, rawMinus = 0;       // cones before dipole cancellation
        int conesPlus = 0, conesMinus = 0;   // after
        int dipoleUnits = 0;
        int gridNx = 0, gridNy = 0;
        double seconds = 0.0;
        std::vector<std::string> messages;
    };

    ReferenceField(std::shared_ptr<Mesh> mesh, const std::vector<int> &interfaceEdges, const Options &opts);

    bool built() const { return report_.built; }
    const Report &getReport() const { return report_; }
    const Options &getOptions() const { return opts_; }
    double diagonal() const { return diag_; }
    double area() const { return area_; }
    // exp(4i phi) per input triangle, unit.
    const Eigen::VectorXcd &triangleValues() const { return u_; }

    // ---- the field at a point and over a cell ------------------------------

    // u at p, bilinear on the lookup grid; |u| <= 1, smaller near cones.
    std::complex<double> valueAt(const Point &p) const;
    // The cross angle at p in [0, pi/2), and |u| there as its coherence.
    double crossAt(const Point &p, double *coherence) const;
    // u_bar over a quadrilateral (counter-clockwise): nine points of its
    // bilinear map, weighted by the map's Jacobian.
    std::complex<double> average(const std::array<Point, 4> &quad) const;
    // c_q = exp(4i theta_q), the quadrilateral's own cross.
    static std::complex<double> cellCross(const std::array<Point, 4> &quad);
    static double quadArea(const std::array<Point, 4> &quad);

    // ---- E_dir -----------------------------------------------------------

    double misalignment(const std::array<Point, 4> &quad) const;       // m_q, area units
    double misalignment(const SquareCarrier &C, int q) const;
    double directionEnergy(const SquareCarrier &C) const;              // sum m_q / |Omega|

    // ---- E_sing -----------------------------------------------------------

    const std::vector<Singularity> &cones() const { return cones_; }
    // The cones inside a region: its outer loop, less its holes.
    std::vector<Singularity> conesIn(const std::vector<Point> &outer,
                                     const std::vector<std::vector<Point>> &holes = {}) const;
    // A carrier's singular vertices as index units. `only`, when given,
    // restricts them to the vertices it marks.
    static std::vector<Singularity> carrierSingularities(const SquareCarrier &C,
                                                         const std::vector<char> *only = nullptr);
    double singularityEnergy(const std::vector<Singularity> &S, const std::vector<Singularity> &F,
                             const SingularityWeights &w, SingularityReport *rep = nullptr) const;
    double singularityEnergy(const SquareCarrier &C, const SingularityWeights &w,
                             SingularityReport *rep = nullptr) const;

private:
    void solve(const std::vector<int> &interfaceEdges);
    void buildGrid();
    void findCones(const DualMBO &f);

    std::shared_ptr<Mesh> mesh_;
    Options opts_;
    Report report_;
    double diag_ = 1.0, area_ = 0.0;
    Eigen::VectorXcd u_;
    std::vector<Singularity> cones_;
    // u on a regular grid over the bounding box (bilinear in between); nodes
    // outside the domain are filled from their neighbours, so a coarse cell
    // whose chord cuts inside a concave arc still reads the arc's direction.
    std::vector<std::complex<double>> grid_;
    int nx_ = 0, ny_ = 0;
    Point lo_{0.0, 0.0};
    double step_ = 1.0;
};

#endif // __REFERENCE_FIELD_HXX__
