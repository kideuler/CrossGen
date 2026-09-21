#ifndef __PLANAR_DOMAIN_HXX__
#define __PLANAR_DOMAIN_HXX__

#include <string>
#include <vector>

#include "mesh/Mesh.hxx"

// Stage 1 of docs/square_transport_2d_theory_and_implementation.md: validate
// and tag the planar domain.
//
// The triangle mesh is the authoritative geometry (Sec. 1.1), and everything
// after this stage trusts it: the carrier of Stage 2 inherits its embedding
// cell by cell, and Sec. 3.1's positive-Jacobian proof is only a proof for a
// nondegenerate, consistently oriented, manifold input. So this stage is where
// an input that does not satisfy those hypotheses is refused, rather than
// discovered three stages later as an inverted element nobody can explain.
//
// ### What is checked, and why each one is fatal
//
//   * Nondegenerate, counter-clockwise triangles. Sec. 3.1 multiplies the
//     reference determinant (3 - u - v)/12 by det[v1 - v0, v2 - v0]; a zero or
//     negative factor is a zero or negative cell.
//   * Manifold edges and vertices. An edge with three triangles has no
//     well-defined neighbour, so the transport g_qr across it is not a function;
//     a pinched vertex (two boundary loops touching, a bow-tie) is excluded by
//     Sec. 1.1 in so many words.
//   * One outer loop per connected component, and V - E + F = 1 - h there.
//     Sec. 4's identity is written per component with h holes, and it is only
//     a check if h is known independently of the thing being checked.
//
// ### What is tagged
//
// Sec. 1.1 separates two requirements: geometric preservation (every boundary
// and interface segment survives, possibly several to a block side) and layout
// preservation (designated corners are macrovertices, designated curves are
// unions of macroedges). This stage supplies the second list:
//
//   * boundary corners, from the turning angle -- a proposal with an explicit
//     threshold, as Stage 1 of Sec. 12 insists, not a classifier;
//   * the material interface network's nodes: junctions (three or more
//     branches), landings on dS, kinks, and dangling ends;
//   * any vertex the caller names.
//
// A polyline vertex is *not* protected merely for being one: a circle cut into
// 64 segments has 64 vertices and no corners, and protecting them would force a
// 64-sided blocking on it (the point Sec. 1.1 makes about geometric versus
// layout preservation).
class PlanarDomain {
public:
    struct Options {
        // A boundary vertex whose interior angle is further than this from pi
        // (degrees) is a protected corner. 30 degrees keeps a circle cut into
        // twelve or more segments free of corners and catches every drawn
        // corner in the corpus.
        double cornerAngle = 30.0;
        // The same test along a material interface: two interface edges meeting
        // further than this from straight make a kink.
        double kinkAngle = 30.0;
        // Treat the edges between different material ids as features.
        bool materialInterfaces = true;
        // Vertices the caller designates as macrovertices, whatever their angle.
        std::vector<int> extraCorners;
        // The interior angle the *layout* is to see at a vertex, where the
        // mesh's own angle is only a proxy for it (NaN, or a short vector,
        // means "the mesh's own"). CoarseDomain needs this: a coarse
        // re-triangulation of a circle is an inscribed polygon whose vertices
        // turn by 45 degrees each, and a layout that believed those angles
        // would put a block corner at every one of them. The curve itself has
        // no corner there, and the fine mesh's angle says so.
        std::vector<double> angleOverride;
        // The same, for the angle the layout sees *along an interface* at a
        // vertex with two interface edges: the kink test reads this where it
        // is finite and the polyline's own turn where it is not. A coarse
        // re-triangulation needs it for the same reason as angleOverride --
        // an interface arc sampled eight times is a polygon that turns at
        // every sample, and believing those turns would protect all eight as
        // kinks and force a block corner at each.
        std::vector<double> interfaceAngleOverride;
        // A triangle whose area is below this fraction of the squared bounding
        // box diagonal is degenerate.
        double degenerateArea = 1e-14;
    };

    enum class CornerKind : unsigned char {
        None,
        BoundaryCorner,     // turning angle over Options::cornerAngle
        InterfaceJunction,  // three or more interface edges meet inside
        InterfaceLanding,   // an interface reaches dS
        InterfaceKink,      // an interface turns by more than kinkAngle
        InterfaceDangling,  // a single interface edge ends here (a tagging error)
        UserDesignated
    };

    // A boundary loop, walked with the domain on the left: counter-clockwise
    // for the outer loop of a component, clockwise for a hole.
    struct Loop {
        std::vector<int> vertices;
        std::vector<int> edges;       // mesh edge from vertices[i] to vertices[i+1]
        int component = -1;
        double signedArea = 0.0;
        bool outer = false;
        int corners = 0;              // protected vertices on it
        double length = 0.0;
    };

    struct Component {
        std::vector<int> triangles;
        int outerLoop = -1;
        std::vector<int> holes;       // loop indices
        int vertices = 0, edges = 0;
        int eulerCharacteristic = 0;  // V - E + F, which must be 1 - h
    };

    struct Report {
        int vertices = 0, edges = 0, triangles = 0, boundaryEdges = 0;
        int degenerateTriangles = 0;
        int invertedTriangles = 0;
        int nonManifoldEdges = 0;
        int nonManifoldVertices = 0;   // pinches and bow-ties
        int components = 0;
        int loops = 0;
        int holes = 0;
        int componentsWithoutOuterLoop = 0;
        int eulerMismatches = 0;       // components with V - E + F != 1 - h
        int materials = 0;
        int interfaceEdges = 0;
        int boundaryCorners = 0;
        int interfaceJunctions = 0, interfaceLandings = 0, interfaceKinks = 0;
        int interfaceDangling = 0;
        int userCorners = 0;
        int protectedVertices = 0;
        double area = 0.0;
        double minTriangleArea = 0.0;
        double minAngleDegrees = 0.0;
        double scale = 0.0;            // bounding-box diagonal
        double meanEdge = 0.0;
        double boundaryLength = 0.0;
        bool valid = false;
        std::vector<std::string> messages;
    };

    PlanarDomain(const Mesh &mesh, const Options &opts);
    explicit PlanarDomain(const Mesh &mesh) : PlanarDomain(mesh, Options()) {}

    const Mesh &getMesh() const { return mesh_; }
    const Options &getOptions() const { return opts_; }
    const Report &getReport() const { return report_; }

    // ---- per vertex ------------------------------------------------------
    std::vector<char> protectedVertex;     // a required macrovertex
    std::vector<CornerKind> cornerKind;
    std::vector<double> interiorAngle;     // sum of incident triangle angles
    // The angle corners and valence targets are read from: interiorAngle,
    // except where Options::angleOverride says otherwise. Geometry (the angle
    // sums of Sec. 9.2) always uses interiorAngle.
    std::vector<double> targetAngle;
    std::vector<int> loopOf;               // boundary loop, -1 if interior
    std::vector<int> interfaceDegree;      // interface edges at the vertex

    // ---- per mesh edge ---------------------------------------------------
    std::vector<char> interfaceEdge;       // two triangles of different material
    std::vector<char> boundaryEdge;

    // ---- per triangle ----------------------------------------------------
    std::vector<int> triangleComponent;
    std::vector<double> triangleArea;      // signed

    std::vector<Loop> loops;
    std::vector<Component> components;

private:
    void checkTriangles();
    void checkManifold();
    void buildLoops();
    void buildComponents();
    void tagInterfaces();
    void tagCorners();

    const Mesh &mesh_;
    Options opts_;
    Report report_;
};

#endif // __PLANAR_DOMAIN_HXX__
