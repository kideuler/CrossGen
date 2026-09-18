#ifndef __CAVITY_FILL_HXX__
#define __CAVITY_FILL_HXX__

#include <array>
#include <string>
#include <unordered_set>
#include <vector>

#include "ATLAS/SquareCarrier.hxx"

// The machinery Stages 3 and 6 share: describing a cavity of the carrier,
// building a structured replacement for it, and certifying that replacement
// before it is allowed anywhere near the carrier.
//
// A cavity is a set of carrier cells. Its replacement (a Patch) is a set of new
// cells over the cavity's own boundary vertices plus new vertices. The boundary
// is never moved: the one change allowed to it is inserting collinear points on
// a *domain* boundary segment (Sec. 8.4 -- "Changing the number or positions of
// boundary subdivisions may require replacing or refining neighboring cells";
// on dS there is no neighbour, so it requires nothing). That single rule is
// what makes Sec. 8.4's three agreements hold by construction rather than by
// checking: the geometric boundary curve, the incidences and the edge
// parameterisations of every retained edge are all literally unchanged,
// because the retained edges are the same pairs of the same vertices.
//
// ### The certificate (Sec. 9.2)
//
// validate() refuses a patch unless all of the following hold:
//
//   1. every new cell is strictly convex, i.e. its bilinear map has a positive
//      Jacobian at all four corners -- Sec. 9.1's algebraic certificate, since
//      det DX is affine in (u, v) and so is positive on [0,1]^2 iff it is at
//      the corners;
//   2. the patch is a consistently oriented, edge-manifold complex whose
//      boundary is exactly the cavity's boundary, each domain-boundary edge
//      possibly subdivided by collinear inserted points;
//   3. every new interior vertex has a one-ring whose angles sum to 2 pi, and
//      every retained boundary vertex has the same interior angle inside the
//      patch as it had inside the cavity (pi at an inserted point);
//   4. the areas agree.
//
// (1) and (3) make the patch an orientation-preserving local homeomorphism,
// and (2) says it restricts to a homeomorphism of the boundary onto the
// cavity's own boundary curves; a local homeomorphism of a compact region that
// is injective on the boundary is injective, so the new cells tile the cavity
// exactly, with no overlap and no gap. That is the whole of Sec. 9.2's
// "simple oriented outer boundary, noncrossing internal edges, disjoint cell
// interiors, exact coverage", obtained from local tests, and it is why (4) is
// a redundant cross-check rather than the test: Sec. 9.2 warns that a gap and
// an overlap can cancel in the area, and they cannot get past (3).
class CavityFill {
public:
    typedef SquareCarrier::Origin Origin;

    // A proposed replacement. Vertex ids are carrier ids when >= 0 and new
    // vertices when < 0 (-1 - k is newVertices[k]), matching SquareCarrier::Edit.
    struct Patch {
        std::vector<Point> newVertices;
        std::vector<Origin> newOrigin;
        std::vector<int> newSourceEdge;
        std::vector<std::array<int, 4>> cells;
        std::vector<int> designated;
        int blocks = 0;          // blocks the pattern is made of
        std::string kind;

        int addVertex(const Point &p, Origin o, int sourceEdge = -1) {
            newVertices.push_back(p);
            newOrigin.push_back(o);
            newSourceEdge.push_back(sourceEdge);
            return -static_cast<int>(newVertices.size());
        }
        Point at(const SquareCarrier &C, int id) const {
            return id >= 0 ? C.vertices[id] : newVertices[-1 - id];
        }
        void set(int id, const Point &p) { newVertices[-1 - id] = p; }
    };

    // The boundary of a set of cells: its loops, walked with the cavity on the
    // left, and what surrounds each loop vertex.
    struct Boundary {
        std::vector<std::vector<int>> loops;       // vertex ids
        std::vector<std::vector<int>> loopEdges;   // carrier edge loops[k][i] -> loops[k][i+1]
        std::vector<std::vector<double>> angle;    // the cavity's interior angle there
        std::vector<std::vector<int>> inside;      // cavity cells at the vertex
        std::vector<std::vector<int>> outside;     // other cells at the vertex
        std::vector<double> signedArea;            // per loop
        std::vector<int> interior;                 // cavity vertices on no loop
        bool manifold = true;
        std::string reason;
    };

    struct Verdict {
        bool valid = false;
        std::string reason;
        double minScaledJacobian = 0.0;
        double meanScaledJacobian = 0.0;
    };

    static Boundary boundaryOf(const SquareCarrier &C, const std::vector<int> &cells,
                               const std::unordered_set<int> &inCavity);

    // A structured n x m grid over a four-sided region. sides[0..3] are the
    // bottom (corner 0 -> 1, n edges), right (1 -> 2, m edges), top (2 -> 3,
    // n edges) and left (3 -> 0, m edges) boundary node ids, consecutive sides
    // sharing their end node. Interior nodes are placed by discrete
    // transfinite interpolation with arc-length blending (the algebraic Coons
    // construction: exact on the four sides, the bilinear map when the sides
    // are straight). Cell (i, j) has corner 0 at node (i, j), so every
    // transport inside the grid is the identity. Returns false when the sides
    // do not close up or opposite counts differ.
    static bool grid(const SquareCarrier &C, Patch &P, const std::array<std::vector<int>, 4> &sides,
                     Origin interiorOrigin, std::vector<int> *nodesOut = nullptr);

    // n edges on the straight segment a -> b: the two ends and n - 1 new points.
    static std::vector<int> segment(const SquareCarrier &C, Patch &P, int a, int b, int n,
                                    Origin o);

    // K blocks round one new centre (K = 3 or 5): arcs[k] runs from corner k to
    // corner k+1, and n_k = arcs[k].size() - 1 must equal s_{k-1} + s_{k+1}.
    // Block k sits at corner k and is bounded by the tail of arc k-1, the
    // head of arc k, and the spokes from the split points m_{k-1}, m_k to the
    // centre -- the three-quad split of Sec. 3 one level up when K = 3.
    static bool star(const SquareCarrier &C, Patch &P, const std::vector<std::vector<int>> &arcs,
                     const std::vector<int> &sigma, const Point &centre, Origin interiorOrigin,
                     std::vector<int> *splitPoints = nullptr, int *centreId = nullptr);

    // n_i = s_{i-1} + s_{i+1} in positive integers. Solvable only for odd K
    // (for K = 4 the circulant is singular -- that is one block, not a star).
    static bool solveStar(const std::vector<int> &n, std::vector<int> &s);

    // The least subdivision of the sides that makes a star solvable: positive
    // integers s with s_{i-1} + s_{i+1} >= n_i, equality where the side may not
    // be subdivided (an interface), minimising the points added and then the
    // distance from the real solution. Exhaustive over the free spoke counts
    // (two of them for K = 3, three in a window round the real solution for
    // K = 5), so a long side opposite two short ones -- the case a fixed
    // budget of increments misses -- is still found. False when no such s.
    static bool repairStar(const std::vector<int> &n, const std::vector<char> &splittable,
                           std::vector<int> &s);

    // Jacobi-Laplacian smoothing of the patch's free vertices: every new vertex
    // except a point inserted on dS. The cavity boundary does not move.
    static void smooth(const SquareCarrier &C, Patch &P, int iterations);

    // Place the patch, then smooth it only if the placement is not already
    // certified, keeping the best certified state seen.
    static Verdict settle(const SquareCarrier &C, const std::vector<int> &cavity, Patch &P,
                          double minScaledJacobian, int smoothingIterations);

    static Verdict validate(const SquareCarrier &C, const std::vector<int> &cavity, const Patch &P,
                            double minScaledJacobian);

    // Min over the loop's edges of the signed distance from c to the edge's
    // line, positive on the left. Positive iff c is strictly inside the
    // loop's kernel, i.e. the loop is star-shaped about c (Sec. 7.3).
    static double kernelDepth(const std::vector<Point> &loop, const Point &c);
    // Approximately the deepest kernel point, by ascent from `start`.
    static Point kernelCenter(const std::vector<Point> &loop, const Point &start);
    static Point areaCentroid(const std::vector<Point> &loop);
    static double signedArea(const std::vector<Point> &loop);
};

#endif // __CAVITY_FILL_HXX__
