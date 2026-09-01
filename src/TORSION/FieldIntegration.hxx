#ifndef __TORSION_FIELDINTEGRATION_HXX__
#define __TORSION_FIELDINTEGRATION_HXX__

#include <array>
#include <string>
#include <vector>

#include "MERIDIAN/ConeCut.hxx"
#include "MERIDIAN/Immersion.hxx"
#include "TORSION/FieldFrames.hxx"
#include "mesh/Mesh.hxx"

// Pipeline B, work item 3 -- docs/cf_flow_pipeline.md Sec. 6: the seamless
// integration of the cross field on Omega.
//
//     min_{u,v}   sum_t A_t ( ||grad u|_t - X_t||^2 + ||grad v|_t - Y_t||^2 )
//     s.t.        phi_+ = R_k phi_- + t_k   on every seam pair of an arc of Gamma_Hol_k
//
// Unconstrained this is a plain cotan-Laplacian Poisson system, L u = div X and
// L v = div Y, one per coordinate and no coupling between them. The seam is
// what couples them, because R_k for odd k exchanges u and v, so the system is
// assembled over all 2n unknowns at once.
//
// ### The constraints, and why they are written on tangents
//
// t_k is an unknown of the problem and not data. Written vertex by vertex the
// constraint set on one arc is phi_+(i) = R_k phi_-(i) + t_k for every paired
// vertex i, which is 2(n+1) equations in the map and 2 more unknowns; subtract
// consecutive pairs and t_k drops out, leaving exactly
//
//     (phi_+(j) - phi_+(i)) = R_k (phi_-(j) - phi_-(i))
//
// per seam edge -- and the equation that was subtracted away is precisely the
// one that *defines* t_k, so nothing has been lost. The two forms have the same
// solution set in phi. This one has no extra unknowns, is homogeneous, and is
// the same expression E4 of Stage 6 is written on, which means the map that
// comes out of here and the map Stage 6 goes on to hold are being held to the
// identical statement. It also makes the cone tip -- where the plus and minus
// chains are the same vertex, and where the vertex form degenerates into
// (I - R_k) phi = t_k -- an ordinary case needing no special handling.
//
// **t_k is left real.** Rounding it is what turns a seamless map into an
// integer-grid map, and buying Q5 with a mixed-integer program is exactly what
// the paper this pipeline feeds is built to avoid: E5 is real-valued and that
// is Sec. 3.3's central numerical claim. Leaving t_k free is the compatibility
// point between the field route and the rest of MERIDIAN.
//
// ### The solve
//
// [[L, C^T], [C, -eps I]] phi/lambda = [b, 0], one factorisation. The energy's
// null space is the constants in u and in v, which the homogeneous tangent
// constraints do not touch, so one vertex is pinned -- two more rows of C -- to
// remove it. The -eps I on the multiplier block makes the saddle matrix quasi
// definite, so an LDL^T without pivoting is stable on it; it also absorbs the
// redundant rows a seam graph with junctions produces, which is what lets
// --cut-to-graph run at all. eps is scaled off the diagonal of L and the seam
// residual it costs is reported, so the price is visible rather than assumed
// small.
//
// ### What comes back, and what it is not
//
// A map with the field's directions in it and no guarantee whatever of local
// injectivity. A discrete cross field is generically non-integrable -- curl X
// is not zero -- so the least-squares closest map inverts triangles, near the
// cones and in high-distortion pockets. That is the entire cost of substituting
// this stage for Ricci flow, and Report::flippedFaces is where the bill arrives.
// Sec. 7.2's untangling is what pays it; zero flips means it need not run, and
// that does happen on gently curved, well-aligned models.
//
// Report::maxFitResidual is the non-integrability itself, per triangle and
// after the fact: how far the gradient the solve settled for is from the
// gradient the field asked for, relative to |X| = 1/h. Where it is large is
// where the flips are, and it is the number Sec. 7.3's refinement should be
// pointed at.
class FieldIntegration {
public:
    struct Options {
        // The vertex of Omega pinned to remove the translation null space, and
        // where it is pinned to.
        int pinnedVertex = 0;
        Point pinnedTo{0.0, 0.0};
        // Multiplier-block regularisation, relative to the mean diagonal of L.
        double regularisation = 1e-10;

        // Sec. 6.4: which coordinate each edge of Omega holds constant --
        // 0 for u, 1 for v, -1 for free -- as FieldFrames::alignmentAxis()
        // reads it off the frame. Empty leaves Q3 and E3 entirely to Stage 6,
        // which is what this stage did before.
        //
        // These rows are of exactly the kind the seam rows already are:
        // homogeneous, on a difference of two vertex values, and an equality
        // rather than a penalty. What they buy is that psi_0 satisfies Q3 --
        // and, at every boundary vertex, Q2 -- the moment the solve returns,
        // instead of after a penalty continuation has driven a residual of 1e-1
        // down to 1e-6. Pipeline A gets the same thing from a different
        // direction: Stage 1 prescribes zero curvature at every non-cone
        // boundary vertex, so the flow makes dS geodesic and psi_R's boundary
        // is a rectilinear polygon before Stage 6 sees it.
        //
        // What they cost is that the least-squares fit no longer gets to choose
        // the boundary. Where the field is genuinely off the axis its chain was
        // given, the map is being asked for something the field did not offer,
        // and Report::alignmentStrain is that disagreement measured after the
        // fact.
        std::vector<int> alignAxis;

        // The Jacobian the fit is asked for, one 2x2 row-major per face of
        // Omega. Empty means the frames', which is Sec. 6 as written.
        //
        // It is here for one caller: the re-projection of Sec. 6.5. Stage 4R
        // returns an injective map that has lost the alignment -- it was built
        // from a Tutte embedding of a circle and E2 only pulled it back part of
        // the way -- and the cheapest way to put the alignment back is to fit
        // *that map's own Jacobian* under the same constraints. Its target is
        // integrable by construction, being a gradient, so the unconstrained
        // answer is the map itself and the constrained one is the nearest map
        // to it that satisfies Q3 and Q4 exactly.
        std::vector<std::array<double, 4>> targetJacobian;

        // det J_t at or below this counts as inverted, relative to the model
        // triangle's area.
        double flipTolerance = 0.0;
    };

    struct Report {
        int vertices = 0;
        int faces = 0;
        int seamPairs = 0;
        int alignedEdges = 0;
        int constraintRows = 0;

        bool factorised = false;
        bool solvedWithLDLT = false;   // false means the LU fallback ran
        bool solved = false;

        // How exactly the seam came out, as the worst
        // ||(phi+(j) - phi+(i)) - R_k (phi-(j) - phi-(i))|| over the seam pairs,
        // divided by the extent of the image. This is the price of the
        // regularisation above and it should be at rounding level.
        double maxSeamResidual = 0.0;

        // The alignment as it came out: the worst held coordinate difference
        // over the aligned edges, divided by the extent -- which should be at
        // rounding, the rows being equalities -- and how far the aligned edges'
        // images ended up from the direction the *field* wanted for them,
        // relative to their length. The second is the price of Sec. 6.4 and is
        // not an error: it is the field's own disagreement with the cone set,
        // moved from Stage 8, where it would have been a spurious node, to
        // here, where it is a strained triangle.
        double maxAlignResidual = 0.0;
        double alignmentStrain = 0.0;

        // Sec. 7.1's census.
        int flippedFaces = 0;
        double flippedArea = 0.0;        // model area of the inverted faces
        double totalArea = 0.0;
        double minAreaRatio = 0.0;       // min over faces of image area / (det J* * model area)
        // Of the inverted faces, how near the nearest cone is, in model units
        // and relative to the model's diagonal. Sec. 7.1 asks for it because
        // "the flips are at the cones" and "the flips are all over" call for
        // different remedies from Sec. 7.3.
        double nearestFlipToCone = 0.0;
        double farthestFlipFromCone = 0.0;
        int flipsAdjacentToCone = 0;     // in the one ring of a cone

        // The non-integrability, per triangle, relative to |X| = 1/h.
        double maxFitResidual = 0.0;
        double meanFitResidual = 0.0;
        int worstFitFace = -1;

        double imageExtent = 0.0;
        bool valid = false;              // solved, seam exact, nothing inverted

        std::vector<std::string> messages;
    };

    // `scaffold` supplies the arcs of G, their (e+, e-) pairing and their
    // quarter turns. Building it costs one Immersion over a throwaway map; all
    // three of those are derived from ConeCut and the frames alone, so what is
    // read here does not depend on the map the scaffold was handed.
    FieldIntegration(const ConeCut &cut, const FieldFrames &frames,
                     const Immersion &scaffold, const Options &opts);
    FieldIntegration(const ConeCut &cut, const FieldFrames &frames,
                     const Immersion &scaffold)
        : FieldIntegration(cut, frames, scaffold, Options()) {}

    const std::vector<Point>& getUV() const { return uv; }
    // Per face of Omega: ||grad u - X||^2 + ||grad v - Y||^2, normalised by
    // |X|^2, so 0 is the field realised exactly and 1 is a gradient the size of
    // the target pointing the wrong way. This is the map of where the field
    // could not be integrated.
    const std::vector<double>& fitResiduals() const { return fitResidual; }
    // Per face: image area over (det J* * model area). Negative is inverted.
    const std::vector<double>& areaRatios() const { return areaRatio; }

    const Report& getReport() const { return report; }

private:
    void assemble(const ConeCut &cut, const FieldFrames &frames, const Immersion &scaffold);
    void measure(const ConeCut &cut, const FieldFrames &frames, const Immersion &scaffold);

    Options opts;
    std::vector<std::array<double, 4>> target;   // the frames', or Options'
    std::vector<Point> uv;
    std::vector<double> fitResidual;
    std::vector<double> areaRatio;
    Report report;
};

#endif // __TORSION_FIELDINTEGRATION_HXX__
