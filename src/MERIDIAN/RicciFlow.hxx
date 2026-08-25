#ifndef __RICCIFLOW_HXX__
#define __RICCIFLOW_HXX__

#include <array>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

// eigen includes
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "mesh/Mesh.hxx"

class ConeSingularities;

// Stage 3 of Shepherd, Gu and Hughes (2022), Section 3.2.1: discrete surface
// Ricci flow. (docs/shepherd2022.pdf; the underlying flow is Gu et al. [86] and
// its generalisation to irregular meshes is Yang et al. [87].)
//
// The input carries the metric it inherits from the plane. What the layout
// needs instead is a *flat cone metric*: curvature zero everywhere except at
// the prescribed cones of Stage 1, where it is exactly (pi/2) I(v). On a
// planar model the interior curvature already is zero, so the flow's whole job
// is to move the curvature that lives on the boundary -- corners, hole rims --
// into the interior cones where the layout wants it. The resulting metric is
// locally but not globally flat, and the immersion it induces overlaps itself
// (the magenta hatching of the paper's Fig. 10); that is expected, not a bug.
//
// The metric is parameterised by one circle radius gamma_i per vertex and a
// fixed angle phi_ij per edge, glued by the law of cosines of Eq. (11):
//
//     l_ij^2 = gamma_i^2 + gamma_j^2 + 2 gamma_i gamma_j cos(phi_ij)
//
// with cos(phi_ij) computed once from the input lengths and never touched
// again. The unknowns are the conformal factors u_i = log(gamma_i) of Eq. (8),
// and the flow of Eq. (9) is
//
//     du_i/dt = Kbar_i - K_i.
//
// Rather than integrate that, Eq. (10) recasts it as the critical point of a
// convex energy, and only the critical point is wanted, so the integral is
// never evaluated -- gradient and Hessian in closed form are enough:
//
//     grad = K - Kbar,     Hess_ij = -w_ij,  Hess_ii = sum_j w_ij,
//     w_ij = ( h_ij^k + h_ij^l ) / l_ij
//
// where h_ij^k is the signed distance from the power centre of triangle ijk to
// the edge ij. That Hessian is the discrete Laplace-Beltrami operator of the
// current metric; with all gamma equal the power centre is the circumcentre and
// it reduces to the familiar cotangent Laplacian, which is the sanity check
// worth keeping in mind for the formulas in the implementation.
//
// Two details decide whether this converges on a real structural mesh.
//
// The inversive distance is not clamped. Thurston's circle packing wants
// cos(phi) in [0, 1], but initialising gamma from the input metric gives
// cos(phi) >= 1 on every edge -- the circles overlap in the inversive sense --
// and clipping that back to 1 throws away the metric it was measured from.
// Yang et al. [87] generalise the flow to admit it, and Sec. 3.2.1 relies on
// that generalisation.
//
// The triangulation is kept weighted-Delaunay. As u moves, a triangle can stop
// satisfying the triangle inequality and the metric stops being realisable at
// all. The fix is to flip to the power-Delaunay triangulation of the current
// metric, the criterion being exactly w_ij >= 0 -- so the same quantity that
// assembles the Hessian decides the flips, and keeping it non-negative is also
// what makes the Hessian positive semi-definite and Newton well behaved. The
// flipped mesh is a different triangulation of the same surface: vertices,
// curvatures and the flat metric are unchanged (a flip re-cuts one planar quad
// along its other diagonal, which preserves every vertex's angle sum), only the
// combinatorics move. So this class keeps its own "flow triangulation" and
// leaves the caller's Mesh alone; see getFlowFaces() and originalEdgeLengths()
// for what comes back out.
class RicciFlow {
public:
    using EdgeKey = MeshEdgeKey;
    using EdgeKeyHash = MeshEdgeKeyHash;

    struct Report {
        int newtonIterations = 0;
        int flips = 0;                  // total edge flips over the whole solve
        int flipPasses = 0;
        bool converged = false;
        double initialError = 0.0;      // ||K - Kbar||_inf before the first step
        double finalError = 0.0;        // ||K - Kbar||_inf at the end
        double gaussBonnetResidual = 0.0; // |sum Kbar - 2 pi chi|; must be ~0 to start
        int nonRealisableFaces = 0;     // faces failing the triangle inequality at the end
        double minInversiveDistance = 0.0; // smallest cos(phi_ij) in play
        double maxInversiveDistance = 0.0;
        // Smallest w_ij at the end, split by where the edge lives. The
        // weighted-Delaunay criterion h^k + h^l >= 0 is a statement about
        // interior edges, which have two triangles to compare; a boundary edge
        // has only h^k and is negative whenever that triangle's opposite angle
        // is obtuse -- normal, unfixable by flipping (there is no quad to
        // re-cut), and the exact analogue of a lone negative cotangent on the
        // boundary of the cotangent Laplacian. So only the interior figure is
        // a health check.
        double minInteriorEdgeWeight = 0.0;
        double minBoundaryEdgeWeight = 0.0;
        int flippedAwayOriginalEdges = 0; // input edges the flow triangulation no longer has
        std::vector<std::string> messages;
    };

    // The target curvature comes straight from Stage 1: Kbar_i = (pi/2) I(v_i).
    // Construct with an admissible cone set -- sum Kbar_i must equal 2 pi chi,
    // or the linear system in every Newton step is inconsistent, since the
    // Laplacian's kernel is the constants and the residual has to be orthogonal
    // to it. solve() refuses to start otherwise.
    RicciFlow(std::shared_ptr<Mesh> mesh, const ConeSingularities &cones);

    // Direct form: one target curvature per vertex.
    RicciFlow(std::shared_ptr<Mesh> mesh, const std::vector<double> &targetCurvature);

    // ||K - Kbar||_inf at which the Newton iteration stops. The paper's own
    // figure is 1e-8, tighter than a graphics application would bother with,
    // because the cone angles have to be exact multiples of pi/2 for Q2 to hold
    // and not merely close.
    void setTolerance(double tol) { tolerance = tol; }
    void setMaxIterations(int n) { maxIterations = n; }

    // The vertex held at u = 0. The energy is invariant under u -> u + c (a
    // global scaling of the metric), so one value has to be pinned or the
    // Hessian is singular. Any vertex will do; the default is 0.
    void setPinnedVertex(int v) { pinnedVertex = v; }

    // Whether to maintain the weighted-Delaunay triangulation by flipping.
    // On by default and it should stay on -- see the class comment. Turning it
    // off is for isolating a convergence failure, not for production.
    void setDelaunayFlips(bool on) { delaunayFlips = on; }

    // Newton on Eq. (10) with a flip-and-realisability line search. Returns
    // true when ||K - Kbar||_inf came under the tolerance.
    bool solve();

    const Report& getReport() const { return report; }

    // The conformal factors of Eq. (8). This is the actual output of the stage:
    // together with the fixed cos(phi_ij) it *is* the flat cone metric.
    const std::vector<double>& getU() const { return u; }
    // gamma_i = exp(u_i), the circle radii.
    std::vector<double> getRadii() const;

    // Curvature of the current metric, Eq. (6), and the target it was driven
    // to. Comparing them is the verification of this stage.
    const std::vector<double>& getCurvature() const { return curvature; }
    const std::vector<double>& getTargetCurvature() const { return target; }
    double maxCurvatureError() const;

    // The flow triangulation: same vertices as the input mesh, possibly
    // different faces after flips, with the metric that Stage 4 lays out.
    const std::vector<Triangle>& getFlowFaces() const { return faces; }
    // Length of every edge of the flow triangulation under the flat metric.
    std::unordered_map<EdgeKey, double, EdgeKeyHash> flowEdgeLengths() const;

    // Flat length of each edge of the *input* mesh, indexed as Mesh::edges.
    // An edge the flow flipped away has no length in this metric without an
    // unfolding, and comes back as -1; Report::flippedAwayOriginalEdges counts
    // them. With no flips this is the whole metric in the caller's indexing.
    std::vector<double> originalEdgeLengths() const;

    // The same, with the flipped-away edges filled in rather than left at -1.
    //
    // Stage 4 lays out the *geometric* triangulation -- Omega is built on the
    // input faces, and the feature and boundary bookkeeping of Stages 5 to 8 is
    // keyed on them -- so it needs a length for every input edge, including the
    // ones the flow no longer carries. The practical note of Sec. 3.2.1 is to
    // keep the two triangulations apart and "map back at the end"; this is that
    // mapping.
    //
    // The fill is not an approximation. cos(phi_ij) is measured once from the
    // input metric and is a property of the *edge*, not of the triangulation it
    // sits in: the flow only ever moves gamma. So an input edge that was
    // flipped away still has the length Eq. (11) gives it from its own
    // cos(phi_ij) and the converged u, which is exactly the length the flow
    // would have produced for it had the flip never happened. What a flip does
    // change is the metric's *realisability* on those faces -- the flip
    // happened because a face stopped closing -- so `recovered` comes back with
    // the count and Immersion reports any face that still fails to close.
    std::vector<double> originalEdgeLengthsCompleted(int *recovered = nullptr) const;

    // Largest relative disagreement between the assembled Hessian and a central
    // finite difference of K with respect to u, over `samples` vertices. The
    // Hessian is the one thing here that cannot be checked by looking at the
    // answer -- a wrong sign or a wrong power-centre formula still converges,
    // just slowly and to the same place -- so this is what the test exercises.
    // Cheap: one perturbation pair per sampled vertex.
    double checkHessian(int samples = 20, double h = 1e-6) const;

    const Mesh& getMesh() const { return *mesh; }
    std::shared_ptr<Mesh> getMeshPtr() const { return mesh; }

private:
    // One edge of the flow triangulation, rebuilt whenever the faces change.
    // (i, j) is the edge as face f0 sees it, so that a flip can name the two
    // apexes k and l without re-deriving the orientation.
    struct FlowEdge {
        int i = -1, j = -1;
        int f0 = -1, f1 = -1;   // f1 < 0 on the boundary
        int k = -1, l = -1;     // apex of f0, apex of f1
        double cosPhi = 1.0;
        double length = 0.0;
        double weight = 0.0;    // w_ij
    };

    void initialiseMetric();
    void rebuildEdges();                 // faces -> edge table, cos(phi) carried over
    void rebuildLengths();               // u, cos(phi) -> l_ij, Eq. (11)
    void rebuildAnglesAndCurvature();    // Eqs. (5) and (6)
    void rebuildWeights();               // h_ij^k and w_ij

    // The same three steps without touching any state, so that the line search
    // and checkHessian() can evaluate a trial u without having to undo it.
    std::vector<double> lengthsAt(const std::vector<double> &uTrial) const;
    void curvatureFrom(const std::vector<double> &len,
                       std::vector<double> &K,
                       std::vector<double> &angles) const;

    // True when every face of the flow triangulation satisfies the triangle
    // inequality strictly at the given lengths.
    bool realisable(const std::vector<double> &len, int *badFace = nullptr) const;

    // Change in the Ricci energy of Eq. (10) between u and u + lambda*du.
    // The energy itself is never needed and never evaluated -- only a critical
    // point is wanted -- but its *difference* along a segment is a 1-D integral
    // of (K - Kbar) . du, which Gauss-Legendre gets to machine precision in a
    // handful of curvature evaluations. That is what makes the line search a
    // real Armijo test rather than a guess.
    double energyDelta(const std::vector<double> &du, double lambda) const;

    // One pass of power-Delaunay restoration. Returns the number of flips made.
    int flipPass();
    bool tryFlip(int edgeIndex);

    // Signed distance from the power centre of triangle (a, b, c) to edge ab,
    // positive on the side of c. Needs the lengths and radii of that triangle.
    static double powerHeight(double lab, double lbc, double lca,
                              double ga, double gb, double gc);

    // Assemble the Laplacian, drop the pinned row and column, solve for the
    // Newton direction. Returns false if the factorisation fails.
    bool newtonDirection(std::vector<double> &du);

    std::shared_ptr<Mesh> mesh;
    // Vertices that carry at least one triangle. An .obj may list vertices no
    // face references, and such a vertex is not part of the surface: it has no
    // circle, no curvature and an all-zero Laplacian row, which would make the
    // Newton system singular however many vertices are pinned. It is kept out
    // of chi, out of the system and out of the error norm.
    std::vector<char> active;
    int activeCount = 0;
    std::vector<double> target;      // Kbar
    std::vector<double> u;           // conformal factors
    std::vector<double> curvature;   // K
    std::vector<double> angleSum;    // sum of interior angles per vertex

    std::vector<Triangle> faces;                 // the flow triangulation
    std::vector<FlowEdge> edges;                 // rebuilt from faces
    // Local edge p of face f joins tri[p] and tri[(p+1)%3], so the edge
    // opposite local vertex m is local edge (m+1)%3.
    std::vector<std::array<int, 3>> faceEdges;
    std::unordered_map<EdgeKey, int, EdgeKeyHash> edgeIndex;
    std::unordered_map<EdgeKey, double, EdgeKeyHash> cosPhi; // persistent across flips
    // cos(phi_ij) as first measured, keyed on the *input* edges. `cosPhi` loses
    // an entry every time a flip retires an edge; this one never changes, which
    // is what originalEdgeLengthsCompleted() reads.
    std::unordered_map<EdgeKey, double, EdgeKeyHash> initialCosPhi;

    double tolerance = 1e-8;
    int maxIterations = 100;
    int pinnedVertex = 0;
    bool delaunayFlips = true;

    Report report;
};

#endif // __RICCIFLOW_HXX__
