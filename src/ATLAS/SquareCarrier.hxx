#ifndef __SQUARE_CARRIER_HXX__
#define __SQUARE_CARRIER_HXX__

#include <array>
#include <string>
#include <vector>

#include "ATLAS/PlanarDomain.hxx"
#include "ATLAS/SquareTransport.hxx"
#include "mesh/Mesh.hxx"

// Stage 2 of docs/square_transport_2d_theory_and_implementation.md: the
// guaranteed carrier, and the square-transport complex every later stage reads
// and rewrites.
//
// ### The construction (Sec. 3)
//
// Each counter-clockwise triangle T = [v0, v1, v2] becomes three quadrilaterals,
//
//     Q_i = [v_i, m_{i,i+1}, c_T, m_{i-1,i}],
//
// the barycentric region where lambda_i is the largest coordinate, with the
// midpoint of each mesh edge built once and shared by both its triangles, so
// the carrier is conforming by construction. Corner 0 of Q_i is v_i and the
// corners run counter-clockwise, which fixes the cell's coordinates: its
// bilinear map is, on the reference triangle,
//
//     X(u, v) = (u (3 - v) / 6, v (3 - u) / 6),    det DX = (3 - u - v) / 12,
//
// so the four corner determinants of every initial cell are exactly
//
//     2 A_T (1/4, 1/6, 1/12, 1/6)
//
// and validate() checks that identity cell by cell rather than just the sign.
// N_Q = 3 N_T, every cell is convex, and the singleton blocking -- every cell
// its own block -- is already a valid, conforming answer. That is the incumbent
// Sec. 6 needs before any search is allowed to run.
//
// ### The field is the charts (Sec. 2)
//
// Nothing here stores a direction. A cell's corner order is its coordinate
// system; the coframe A_q = (DX_q)^-1 is a property of the four positions; and
// the transition to each neighbour is the integer transform of
// SquareTransport.hxx, fixed by the two side indices alone. The stored
// `transport` table is that transform per side, and validate() checks both
// halves of Sec. 2.1's requirement on it: g_rq = g_qr^-1, and g_qr carries
// r's copy of the shared edge onto q's, endpoint by endpoint.
//
// ### Edits (Secs. 8.4 and 12, "Replacement transaction")
//
// Stages 3 and 6 replace cavities: a set of cells goes, new vertices and cells
// come. An Edit names both; apply() performs a batch of them and rebuilds all
// derived topology. It does no geometric checking itself -- that is
// CavityFill::validate(), which every caller runs *before* building an Edit,
// so an invalid state is never committed and no rollback is needed. The one
// invariant apply() does rely on is that the domain boundary is only ever
// re-subdivided: collinear points inserted on an input boundary segment, or
// ones strictly inside it dropped, never an input vertex (Sec. 1.1's
// geometric preservation). validate() audits that from the provenance each
// vertex carries.
//
// A carrier may also be given outright as quadrilaterals (Quads): that is how
// Realisation hands over a layout found on a CoarseDomain, and validate() is
// then the whole verdict on it -- the same audit, with nothing assumed.
class SquareCarrier {
public:
    enum class Origin : unsigned char {
        MeshVertex,     // an input vertex, same position
        EdgeMidpoint,   // m_ij of Sec. 3
        Centroid,       // c_T of Sec. 3
        Template,       // an interior vertex of a Stage 3 replacement
        BoundarySplit,  // a point inserted on a domain boundary segment (Sec. 8.4)
        Rewrite,        // an interior vertex of a Stage 6 rewrite
        Realised,       // an interior vertex of a coarse layout realised here (Realisation)
        // A point inserted along a material interface. It is an interior
        // vertex of the domain, so nothing about dS applies to it, but it
        // lies on a curve the layout must follow and so is held fixed by
        // every repair, exactly as a BoundarySplit is.
        InterfaceSplit
    };

    enum class CellOrigin : unsigned char { Split, Template, Rewrite, Realised };

    // A carrier given outright as quadrilaterals, rather than split from the
    // domain's triangles: what Realisation builds when it carries a coarse
    // layout onto the fine domain. Every array is per vertex or per cell, as
    // in the authoritative state below; a boundary vertex must be an input
    // vertex (sourceVertex) or lie on an input boundary edge (sourceEdge), and
    // validate() audits exactly that.
    struct Quads {
        std::vector<Point> vertices;
        std::vector<Origin> origin;
        std::vector<int> sourceVertex;
        std::vector<int> sourceEdge;
        std::vector<char> designated;
        std::vector<std::array<int, 4>> cells;    // counter-clockwise
        std::vector<int> material;
        std::vector<int> group;
    };

    struct Options {
        // Tolerances of the validator. The angle sums of Sec. 9.2 are sums of
        // a handful of acos values, good to ~1e-12; 1e-7 separates rounding
        // from a genuinely doubled or missing sector by six orders.
        double angleTolerance = 1e-7;
        double areaTolerance = 1e-9;   // relative
    };

    // A replacement transaction. Vertex ids in `cells` and `designate` are
    // carrier ids when >= 0 and new vertices when < 0: -1 - k is
    // newVertices[k]. Every edit in a batch has its own new-vertex numbering.
    struct Edit {
        std::vector<int> removeCells;
        std::vector<Point> newVertices;
        std::vector<Origin> newOrigin;
        std::vector<int> newSourceEdge;           // mesh edge a BoundarySplit lies on
        std::vector<std::array<int, 4>> cells;    // counter-clockwise
        int material = 1;
        // Per cell when not empty, in place of `material`: a Stage 6 cavity
        // that straddles an interface refills both materials at once.
        std::vector<int> cellMaterials;
        std::vector<int> designate;               // macrovertices the replacement made
        CellOrigin cellOrigin = CellOrigin::Template;
        int group = -1;                           // which replacement, for reporting
    };

    struct Report {
        int vertices = 0, cells = 0, edges = 0;
        int boundaryEdges = 0, interfaceEdges = 0;
        int sourceTriangles = 0;
        int splitCells = 0, templateCells = 0, rewriteCells = 0, realisedCells = 0;
        int components = 0, holes = 0;
        int nonManifoldEdges = 0;
        int nonManifoldVertices = 0;
        int nonConvexCells = 0;              // some corner determinant <= 0
        double minScaledJacobian = 0.0;
        double meanScaledJacobian = 0.0;
        // Sec. 3.1: only while the carrier is still the initial split.
        bool initialSplit = false;
        int jacobianBoundViolations = 0;
        double minCornerDetRatio = 0.0;      // corner det / 2 A_T, expected >= 1/12
        double maxCornerDetRatio = 0.0;      // expected <= 1/4
        int transportMismatches = 0;
        int angleSumViolations = 0;          // Sec. 9.2 local embedding
        double area = 0.0;
        double areaError = 0.0;              // relative, against the domain
        // Sec. 4: sum_int (4 - q_v) + sum_bdy (2 - q_v) and 4 sum_c (1 - h_c).
        int eulerLHS = 0, eulerRHS = 0;
        bool eulerHolds = false;
        bool boundaryParityEven = false;
        int boundaryPreservationErrors = 0;  // a boundary edge off every source segment
        int missingBoundaryVertices = 0;     // an input boundary vertex that vanished
        int interfacePreservationErrors = 0; // the same two, for the interfaces
        int missingInterfaceVertices = 0;
        int irregularInterior = 0;           // valence != 4
        int irregularBoundary = 0;           // valence != target on dS
        int totalDefect = 0;                 // sum |q_v - target_v|
        int protectedVertices = 0, designatedVertices = 0;
        bool valid = false;
        std::vector<std::string> messages;
    };

    SquareCarrier(const PlanarDomain &domain, const Options &opts);
    explicit SquareCarrier(const PlanarDomain &domain) : SquareCarrier(domain, Options()) {}
    // A carrier from explicit quadrilaterals over the same domain. Nothing is
    // assumed about them: validate() is the whole verdict.
    SquareCarrier(const PlanarDomain &domain, const Options &opts, const Quads &quads);

    const PlanarDomain &getDomain() const { return *domain_; }
    const Options &getOptions() const { return opts_; }
    const Report &getReport() const { return report_; }

    // Recompute the report. Returns it; report.valid is the verdict.
    const Report &validate();

    // Apply a batch of cavity-disjoint edits and rebuild the topology.
    void apply(const std::vector<Edit> &edits);

    // ---- authoritative state ---------------------------------------------

    std::vector<Point> vertices;
    std::vector<Origin> vertexOrigin;
    std::vector<int> sourceVertex;       // input vertex, or -1
    std::vector<int> sourceEdge;         // input edge a midpoint or split lies on, or -1
    std::vector<char> protectedVertex;   // required macrovertex (Stage 1)
    std::vector<char> designatedVertex;  // macrovertex a replacement created (soft)

    std::vector<std::array<int, 4>> cells;
    std::vector<int> cellMaterial;
    std::vector<int> cellTriangle;       // source triangle, -1 after a replacement
    std::vector<CellOrigin> cellOrigin;
    std::vector<int> cellGroup;          // replacement id, -1 for split cells

    // ---- derived by rebuild() --------------------------------------------

    std::vector<std::array<int, 2>> edges;             // (min, max)
    std::vector<std::array<int, 4>> cellEdges;         // side i: corner i -> i+1
    std::vector<std::array<int, 4>> neighbor;          // cell across side i, -1 on dS
    std::vector<std::array<int, 4>> neighborSide;      // its side index
    std::vector<std::array<SquareTransport, 4>> transport;  // g_{q, neighbor}
    std::vector<std::array<int, 2>> edgeCell;          // the one or two cells
    std::vector<std::array<int, 2>> edgeSide;
    std::vector<char> boundaryEdge;
    std::vector<char> interfaceEdge;
    std::vector<char> boundaryVertex;
    std::vector<int> interfaceDegree;
    std::vector<int> valence;                          // incident cells
    std::vector<char> manifoldVertex;
    std::vector<int> cellComponent;
    int componentCount = 0;
    std::vector<int> componentHoles;

    // Incident cells of each vertex in counter-clockwise order, with the
    // corner the vertex occupies; around a vertex on dS the walk starts at the
    // cell whose side leaving the vertex is a boundary edge.
    std::vector<int> ringPtr;
    std::vector<int> ringCell;
    std::vector<int> ringCorner;

    // ---- queries ---------------------------------------------------------

    int numVertices() const { return static_cast<int>(vertices.size()); }
    int numCells() const { return static_cast<int>(cells.size()); }
    int numEdges() const { return static_cast<int>(edges.size()); }

    bool isFeatureEdge(int e) const { return boundaryEdge[e] || interfaceEdge[e]; }
    int otherEnd(int e, int v) const { return edges[e][0] == v ? edges[e][1] : edges[e][0]; }
    // The edges at v in counter-clockwise order: valence of them inside,
    // valence + 1 on dS (from one boundary edge round to the other).
    void vertexEdges(int v, std::vector<int> &out) const;
    int edgeBetween(int a, int b) const;
    int cornerOf(int q, int v) const;

    // Corner c's determinant of the bilinear map, the cross product of the two
    // sides leaving it (Sec. 9.1: the map is positive iff all four are).
    double cornerDet(int q, int c) const;
    double cornerAngle(int q, int c) const;
    double scaledJacobian(int q, int c) const;
    double minScaledJacobian(int q) const;
    double cellArea(int q) const;
    Point cellCentroid(int q) const;

    // The domain's interior angle at a vertex on dS: the input's own at an
    // input vertex, pi at a midpoint or an inserted point. This is geometry,
    // what the cells' angles must sum to there.
    double boundaryAngle(int v) const;
    // The angle the layout sees there (PlanarDomain::targetAngle): the same
    // unless the domain is a coarse proxy for a curved boundary.
    double layoutAngle(int v) const;
    // Regular valence: 4 inside, round(layoutAngle / (pi/2)) clamped to
    // [1, 4] on dS.
    int targetValence(int v) const;
    int defect(int v) const;
    // A vertex no block may run through (Sec. 4 and 5.2): irregular, a
    // protected vertex, or a feature vertex whose sectors are not 2 + 2.
    bool forcedMacrovertex(int v) const;

    double meanEdgeLength() const;

    bool writeOBJ(const std::string &path) const;

private:
    void rebuild();
    void buildRings();
    void buildComponents();

    const PlanarDomain *domain_;
    Options opts_;
    Report report_;
    bool initialSplit_ = true;
};

#endif // __SQUARE_CARRIER_HXX__
