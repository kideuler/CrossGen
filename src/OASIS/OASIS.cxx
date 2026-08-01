#include "OASIS.hxx"

#include <algorithm>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>

namespace {

// Number of rings used for the local quadratic fit of Eq. 9. The paper uses a
// "diameter-2" neighborhood: two rings at boundary vertices, one ring in the
// interior. Only boundary vertices are fitted here, so two rings it is.
constexpr int kBoundaryFitRings = 2;

// ...but a two-ring stencil is only 8-10 points for 6 unknowns, and if those
// points happen to lie on a conic (a ring of vertices following a circular
// boundary, for instance) the design matrix of Eq. 9 is *exactly* rank
// deficient. Widen the stencil in that case rather than giving up on the
// vertex: a boundary vertex whose alignment constraints are dropped is left
// unpinned, and the KKT solve then diverges there.
constexpr int kMaxFitRings = 4;

// xi of Sec. 3.1. Any xi > 0 only rescales f globally and leaves the critical
// points, and therefore the Morse-Smale complex, unchanged.
constexpr double kXi = 1.0;

constexpr double kEps = 1e-14;

double maxAbsEntry(const Eigen::SparseMatrix<double> &A) {
    double m = 0.0;
    for (int k = 0; k < A.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
            m = std::max(m, std::fabs(it.value()));
    return m;
}

} // namespace

OASIS::OASIS(std::shared_ptr<CrossField> crossField, double lambda)
    : mesh(crossField ? crossField->getMeshPtr() : nullptr),
      crossField(crossField),
      lambda(lambda) {
    if (!mesh) {
        throw std::invalid_argument("OASIS: cross field has no mesh");
    }
    r = Eigen::VectorXd::Ones(static_cast<int>(mesh->vertices.size()));
}

OASIS::OASIS(std::shared_ptr<Mesh> mesh, double lambda)
    : mesh(mesh), crossField(nullptr), lambda(lambda) {
    if (!mesh) {
        throw std::invalid_argument("OASIS: null mesh");
    }
    r = Eigen::VectorXd::Ones(static_cast<int>(mesh->vertices.size()));
    // No guiding field, so no orientation term: this is the pure 2014 system.
    orientationWeight = 0.0;
}

void OASIS::setDensity(const Eigen::VectorXd &density) {
    if (density.size() != static_cast<Eigen::Index>(mesh->vertices.size())) {
        throw std::invalid_argument("OASIS::setDensity: expected one entry per vertex");
    }
    if (density.minCoeff() <= 0.0) {
        throw std::invalid_argument("OASIS::setDensity: density must be strictly positive");
    }
    r = density;
    assembled = false;
}

void OASIS::setOrientationMask(const Eigen::VectorXd &mask) {
    if (mask.size() != static_cast<Eigen::Index>(mesh->vertices.size())) {
        throw std::invalid_argument("OASIS::setOrientationMask: expected one entry per vertex");
    }
    if (mask.minCoeff() < 0.0) {
        throw std::invalid_argument("OASIS::setOrientationMask: mask must be non-negative");
    }
    orientationMask = mask;
    assembled = false;
}

Eigen::VectorXd OASIS::boundaryClearanceMask(double clearance) const {
    const int numVertices = static_cast<int>(mesh->vertices.size());
    Eigen::VectorXd mask = Eigen::VectorXd::Ones(numVertices);
    if (clearance <= 0.0) return mask;

    // Brute force against every boundary vertex. |V| * |dV| is a few million
    // distance evaluations on the meshes this runs on, once per assembly.
    for (int v = 0; v < numVertices; ++v) {
        double best = std::numeric_limits<double>::max();
        for (int b : mesh->boundaryVertices) {
            best = std::min(best, normP(mesh->vertices[v] - mesh->vertices[b]));
            if (best <= clearance) break;
        }
        mask[v] = (best > clearance) ? 1.0 : 0.0;
    }
    return mask;
}

void OASIS::assemble() {
    if (lambda >= 0.0) {
        std::cerr << "OASIS: lambda should be negative (Eq. 1), got " << lambda << std::endl;
    }
    buildDensityLaplacian();
    buildConstraints();
    buildOrientation();
    buildKKT();
    assembled = true;
}

void OASIS::solve() {
    if (!assembled) assemble();

    const int numVertices = static_cast<int>(mesh->vertices.size());
    const int nc2 = static_cast<int>(B.rows());

    solverLU.compute(KKT);
    if (solverLU.info() != Eigen::Success) {
        throw std::runtime_error("OASIS::solve: SparseLU factorization of the KKT matrix failed");
    }

    const Eigen::VectorXd z = solverLU.solve(kktRhs);
    if (solverLU.info() != Eigen::Success) {
        throw std::runtime_error("OASIS::solve: SparseLU solve of the KKT system failed");
    }

    f = z.head(numVertices);
    nu = z.tail(nc2);

    // SparseLU reports Success even when it has factorized a singular matrix,
    // in which case f comes back at ~1e17 and every constraint is violated. The
    // constraints are the one thing the solution is guaranteed to satisfy
    // exactly, so check them rather than trusting the solver's status.
    if (nc2 > 0) {
        const double bcError = (B * f - C).cwiseAbs().maxCoeff();
        const double amplitude = std::max(1.0, f.cwiseAbs().maxCoeff());
        if (!f.allFinite() || bcError > 1e-6 * amplitude) {
            throw std::runtime_error(
                "OASIS::solve: the KKT system is singular or ill-conditioned "
                "(boundary conditions violated by " + std::to_string(bcError) + ")");
        }
    }
}

double OASIS::residual() const {
    if (f.size() != static_cast<Eigen::Index>(mesh->vertices.size())) return -1.0;
    return (Lr * f - lambda * f).norm();
}

std::vector<double> OASIS::computeVertexAreas() const {
    const int numVertices = static_cast<int>(mesh->vertices.size());
    const int numTriangles = static_cast<int>(mesh->triangles.size());

    std::vector<double> areas(numVertices, 0.0);
    for (int t = 0; t < numTriangles; ++t) {
        const Triangle &tri = mesh->triangles[t];
        const Point &p0 = mesh->vertices[tri[0]];
        const Point &p1 = mesh->vertices[tri[1]];
        const Point &p2 = mesh->vertices[tri[2]];

        const double area = 0.5 * std::fabs(cross2(p1 - p0, p2 - p0));
        if (area < kEps) continue;

        for (int k = 0; k < 3; ++k) {
            areas[tri[k]] += area / 3.0;
        }
    }
    return areas;
}

// Eq. 8: the cotangent formula modulated by the density r,
//
//   (grad^2_r f)_i = 3 / (2 * Area_i * r_i) * sum_{j in N(i)} (cot a_ij + cot b_ij) (f_j - f_i)
//
// where Area_i is the total area of the triangles incident to vertex i and
// a_ij, b_ij are the angles opposite edge ij. The row scaling makes L_r
// nonsymmetric, which is why Eq. 14 carries both L_r and L_r^T.
void OASIS::buildDensityLaplacian() {
    const int numVertices = static_cast<int>(mesh->vertices.size());
    const int numTriangles = static_cast<int>(mesh->triangles.size());

    const std::vector<double> vertexArea = computeVertexAreas();

    // Unscaled cotangent weights.
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(numTriangles * 12);

    for (int t = 0; t < numTriangles; ++t) {
        const Triangle &tri = mesh->triangles[t];
        const Point &p0 = mesh->vertices[tri[0]];
        const Point &p1 = mesh->vertices[tri[1]];
        const Point &p2 = mesh->vertices[tri[2]];

        const double area2 = std::fabs(cross2(p1 - p0, p2 - p0));
        if (0.5 * area2 < kEps) continue;

        // The angle at corner k is opposite edge (i, j).
        for (int k = 0; k < 3; ++k) {
            const int i = tri[(k + 1) % 3];
            const int j = tri[(k + 2) % 3];

            const Point e1 = mesh->vertices[i] - mesh->vertices[tri[k]];
            const Point e2 = mesh->vertices[j] - mesh->vertices[tri[k]];

            // cot = cos / sin = dot / |cross|, and |cross| is 2 * area for
            // every corner of the triangle.
            const double cot = dotP(e1, e2) / area2;

            trips.emplace_back(i, j, cot);
            trips.emplace_back(i, i, -cot);
            trips.emplace_back(j, i, cot);
            trips.emplace_back(j, j, -cot);
        }
    }

    Lr.resize(numVertices, numVertices);
    Lr.setFromTriplets(trips.begin(), trips.end());

    // Row scaling 3 / (2 * Area_i * r_i) of Eq. 8, written as 1 / (2 D_i r_i)
    // since D_i is a third of the total incident area Area_i.
    Eigen::VectorXd rowScale(numVertices);
    for (int v = 0; v < numVertices; ++v) {
        const double denom = 2.0 * vertexArea[v] * r[v];
        rowScale[v] = (denom > kEps) ? (1.0 / denom) : 0.0;
    }
    Lr = rowScale.asDiagonal() * Lr;
    Lr.makeCompressed();
}

std::vector<Point> OASIS::computeVertexNormals() const {
    const int numVertices = static_cast<int>(mesh->vertices.size());

    std::vector<Point> normals(numVertices, Point{0.0, 0.0});
    std::vector<double> weights(numVertices, 0.0);

    for (int edgeIdx : mesh->boundaryEdges) {
        // A boundary edge belongs to exactly one triangle; read the edge off
        // that triangle so the CCW orientation is the triangle's.
        const int triIdx = mesh->edgeTriangles[edgeIdx][0];
        if (triIdx < 0) continue;
        const Triangle &tri = mesh->triangles[triIdx];

        int localEdge = -1;
        for (int e = 0; e < 3; ++e) {
            if (mesh->triangleEdges[triIdx][e] == edgeIdx) {
                localEdge = e;
                break;
            }
        }
        if (localEdge < 0) continue;

        const int v0 = tri[localEdge];
        const int v1 = tri[(localEdge + 1) % 3];

        const Point e = mesh->vertices[v1] - mesh->vertices[v0];
        const double len = normP(e);
        if (len < kEps) continue;

        // Outward normal of a CCW triangle: the edge rotated clockwise.
        const Point n = Point{e[1] / len, -e[0] / len};

        for (int v : {v0, v1}) {
            normals[v] = normals[v] + n * len;
            weights[v] += len;
        }
    }

    for (int v = 0; v < numVertices; ++v) {
        if (weights[v] > kEps) {
            normals[v] = normalizeP(normals[v] / weights[v]);
        }
    }
    return normals;
}

std::vector<int> OASIS::gatherNeighborhood(int v, int rings) const {
    std::vector<int> collected;
    std::unordered_set<int> seen;

    collected.push_back(v);
    seen.insert(v);

    std::vector<int> frontier{v};
    std::vector<int> next;

    for (int ring = 0; ring < rings; ++ring) {
        next.clear();
        for (int u : frontier) {
            auto range = mesh->vertexTriangles.trianglesForVertex(u);
            for (const int *it = range.first; it != range.second; ++it) {
                const Triangle &tri = mesh->triangles[*it];
                for (int k = 0; k < 3; ++k) {
                    if (seen.insert(tri[k]).second) {
                        collected.push_back(tri[k]);
                        next.push_back(tri[k]);
                    }
                }
            }
        }
        frontier = next;
    }
    return collected;
}

// Eq. 9: fit a quadratic to the neighborhood in a local frame and read the
// coefficients back as linear forms in f.
bool OASIS::fitLocalQuadratic(int v, std::vector<int> &neighbors, Eigen::MatrixXd &coeffs,
                              int &ringsUsed, int startRings) const {
    for (int rings = startRings; rings <= kMaxFitRings; ++rings) {
        neighbors = gatherNeighborhood(v, rings);
        ringsUsed = rings;

        const int m = static_cast<int>(neighbors.size());
        if (m < NUM_COEFF) continue;

        // The paper parameterizes the neighborhood into the tangent plane with
        // the discrete exponential map of Schmidt et al. The domain here is
        // already planar, so that map is the identity and the local coordinates
        // are just offsets from p_v. Curved surfaces would replace this block.
        const Point &pv = mesh->vertices[v];

        std::vector<Point> local(m);
        double scale = 0.0;
        for (int j = 0; j < m; ++j) {
            local[j] = mesh->vertices[neighbors[j]] - pv;
            scale += normP(local[j]);
        }
        // Fit in units of the mean stencil radius; a stencil of size ~1e-3
        // would otherwise put ~1e-12 next to 1 in the design matrix.
        scale /= (m - 1);
        if (scale < kEps) continue;

        Eigen::MatrixXd A(m, static_cast<int>(NUM_COEFF));
        for (int j = 0; j < m; ++j) {
            const double u = local[j][0] / scale;
            const double w = local[j][1] / scale;
            A(j, AUU) = 0.5 * u * u;
            A(j, AUV) = u * w;
            A(j, AVV) = 0.5 * w * w;
            A(j, AU) = u;
            A(j, AV) = w;
            A(j, AC) = 1.0;
        }

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
        const bool full = svd.rank() >= static_cast<Eigen::Index>(NUM_COEFF);
        if (!full && rings < kMaxFitRings) continue;  // widen and retry

        // Solving against the identity gives the least-squares pseudo-inverse,
        // one column at a time: P = A^+, so a = P f[neighbors]. If the stencil
        // is still degenerate at the widest ring, this is the minimum-norm fit
        // rather than a unique one — an approximate constraint, but the vertex
        // stays pinned, which is what matters for the solve.
        coeffs = svd.solve(Eigen::MatrixXd::Identity(m, m));

        // Undo the coordinate scaling: a quadratic term picked up scale^2, a
        // linear term scale, the constant nothing.
        coeffs.row(AUU) /= scale * scale;
        coeffs.row(AUV) /= scale * scale;
        coeffs.row(AVV) /= scale * scale;
        coeffs.row(AU) /= scale;
        coeffs.row(AV) /= scale;

        return true;
    }
    return false;
}

// Eqs. 10 and 11. With the fit coefficients in hand the directional derivatives
// at a boundary vertex with normal (nu, nv) are
//
//   df/dn    ~  (nu a_u + nv a_v)^T f
//   d^2f/dn^2 ~ (nu^2 a_uu + 2 nu nv a_uv + nv^2 a_vv)^T f,
//
// stacked into B f = [Y; Z] f = [0; xi] = C.
void OASIS::buildConstraints() {
    const int numVertices = static_cast<int>(mesh->vertices.size());
    const std::vector<Point> normals = computeVertexNormals();

    constrainedVertices.clear();
    constrainedNormals.clear();

    // Boundary curves only. Interior feature curves obey the same conditions
    // (Sec. 3.3) but need the one-sided stencil of Fig. 5(b), which splits the
    // neighborhood along the feature line; not handled yet.
    std::vector<Eigen::Triplet<double>> yTrips;
    std::vector<Eigen::Triplet<double>> zTrips;
    std::vector<int> neighbors;
    Eigen::MatrixXd coeffs;
    int widened = 0;
    int dropped = 0;

    for (int v : mesh->boundaryVertices) {
        if (normP(normals[v]) < kEps) continue;
        int ringsUsed = kBoundaryFitRings;
        if (!fitLocalQuadratic(v, neighbors, coeffs, ringsUsed, kBoundaryFitRings)) {
            ++dropped;
            continue;
        }
        if (ringsUsed > kBoundaryFitRings) ++widened;

        const int row = static_cast<int>(constrainedVertices.size());
        const double nu = normals[v][0];
        const double nw = normals[v][1];

        for (int j = 0; j < static_cast<int>(neighbors.size()); ++j) {
            const double y = nu * coeffs(AU, j) + nw * coeffs(AV, j);
            const double z = nu * nu * coeffs(AUU, j)
                           + 2.0 * nu * nw * coeffs(AUV, j)
                           + nw * nw * coeffs(AVV, j);
            yTrips.emplace_back(row, neighbors[j], y);
            zTrips.emplace_back(row, neighbors[j], z);
        }

        constrainedVertices.push_back(v);
        constrainedNormals.push_back(normals[v]);
    }

    const int nc = static_cast<int>(constrainedVertices.size());
    if (widened > 0) {
        std::cerr << "OASIS: widened the fit stencil at " << widened << " of " << nc
                  << " boundary vertices (degenerate two-ring neighborhood)" << std::endl;
    }
    if (dropped > 0) {
        // Every dropped vertex is one the solve leaves unpinned, which shows up
        // as a blow-up in that region rather than a graceful degradation.
        std::cerr << "OASIS: WARNING - dropped alignment constraints at " << dropped
                  << " boundary vertices; the solution there is not trustworthy"
                  << std::endl;
    }
    if (nc == 0) {
        std::cerr << "OASIS: no alignment constraints were built; the Eq. 13 "
                     "system is singular and admits only f = 0" << std::endl;
    }

    // Y occupies rows [0, nc), Z rows [nc, 2 nc).
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(yTrips.size() + zTrips.size());
    trips.insert(trips.end(), yTrips.begin(), yTrips.end());
    for (const auto &t : zTrips) {
        trips.emplace_back(t.row() + nc, t.col(), t.value());
    }

    B.resize(2 * nc, numVertices);
    B.setFromTriplets(trips.begin(), trips.end());
    B.makeCompressed();

    C.resize(2 * nc);
    C.head(nc).setZero();  // df/dn = 0
    C.tail(nc).setConstant(kXi);  // d^2f/dn^2 = xi

    dropDependentConstraints();
}

// The KKT matrix of Eq. 13 is nonsingular only if B has full row rank. That
// holds for a well-shaped mesh, but not always: where the domain is only a
// couple of triangles thick, the fit stencils of nearby boundary vertices
// overlap so heavily that their constraint rows become linearly dependent, and
// B loses rank. A singular KKT matrix is worse than it sounds, because
// SparseLU factorizes it and reports success, returning a solution of order
// 1e17 instead of failing.
//
// Dependent rows carry no information the kept rows do not, so drop them and
// keep an equivalent full-rank system. The one caveat: if a dropped row was
// *inconsistent* with the rows implying it, the original constraints were
// infeasible and no formulation would have satisfied them all; here that
// inconsistency is resolved silently in favour of the kept rows.
void OASIS::dropDependentConstraints() {
    if (B.rows() == 0) return;

    // Column-pivoted QR of B^T: the first `rank` pivots select independent
    // columns of B^T, i.e. independent rows of B.
    Eigen::SparseMatrix<double> Bt = B.transpose();
    Bt.makeCompressed();

    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> qr;
    qr.compute(Bt);
    if (qr.info() != Eigen::Success) {
        std::cerr << "OASIS: rank check of the constraint matrix failed; "
                     "proceeding with all rows" << std::endl;
        return;
    }

    const Eigen::Index rank = qr.rank();
    if (rank >= B.rows()) return;  // already full row rank, nothing to do

    std::vector<int> keep;
    keep.reserve(static_cast<size_t>(rank));
    for (Eigen::Index k = 0; k < rank; ++k) {
        keep.push_back(static_cast<int>(qr.colsPermutation().indices()(k)));
    }
    std::sort(keep.begin(), keep.end());

    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(static_cast<size_t>(B.nonZeros()));
    std::vector<int> newRow(static_cast<size_t>(B.rows()), -1);
    for (size_t i = 0; i < keep.size(); ++i) newRow[keep[i]] = static_cast<int>(i);

    for (int col = 0; col < B.outerSize(); ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(B, col); it; ++it) {
            const int nr = newRow[static_cast<size_t>(it.row())];
            if (nr >= 0) trips.emplace_back(nr, static_cast<int>(it.col()), it.value());
        }
    }

    Eigen::VectorXd newC(static_cast<Eigen::Index>(keep.size()));
    for (size_t i = 0; i < keep.size(); ++i) newC[static_cast<Eigen::Index>(i)] = C[keep[i]];

    std::cerr << "OASIS: constraint matrix was rank deficient (" << rank << " of "
              << B.rows() << " rows independent); dropped the redundant rows"
              << std::endl;

    B.resize(rank, static_cast<Eigen::Index>(mesh->vertices.size()));
    B.setFromTriplets(trips.begin(), trips.end());
    B.makeCompressed();
    C = newC;
}

// ─── Orientation control, Huang et al. 2008 Sec. 3.4 / 2014 Sec. 5.1 ────────

bool OASIS::crossFrame(int v, double &cos2t, double &sin2t) const {
    if (!crossField) return false;

    const Eigen::VectorXcd &u =
        (crossField->u_k.size() == static_cast<Eigen::Index>(mesh->vertices.size()))
            ? crossField->u_k
            : crossField->u_k_prev;
    if (u.size() != static_cast<Eigen::Index>(mesh->vertices.size())) return false;

    const std::complex<double> z = u[v];
    if (std::abs(z) < 1e-12) return false;  // singularity: no direction here

    // u = e^{4 i t}, so the guiding direction is t = arg(u)/4 modulo pi/2 and
    // the frame rotation the Hessian needs is 2t = arg(u)/2. Which of the four
    // branches arg() lands on is irrelevant: they differ by pi/2 in t, hence by
    // pi in 2t, which only flips the sign of the row built below.
    const double twoT = 0.5 * std::arg(z);
    cos2t = std::cos(twoT);
    sin2t = std::sin(twoT);
    return true;
}

// Eq. 10 of Huang et al. 2008. The paper re-parameterizes the neighborhood of
// every vertex with its u-axis along the guiding direction and then penalizes
// the mixed coefficient c_uv of the Eq. 9 fit. Refitting per vertex is not
// necessary here: the Hessian is a tensor, so its entries in the guiding frame
// are a fixed linear combination of the entries in the frame the fit already
// uses. With the fit frame axes (x, y), the guiding direction at angle t and
// H = [[h_xx, h_xy], [h_xy, h_yy]],
//
//   c_uv = e_u^T H e_v = cos(2t) h_xy - 1/2 sin(2t) (h_xx - h_yy),
//
// which stays a linear form in f, exactly as the paper needs. Each row is
// scaled by sqrt(D_i) so that ||Q_orient f||^2 is the area-weighted sum.
void OASIS::buildOrientation() {
    const int numVertices = static_cast<int>(mesh->vertices.size());

    Qorient.resize(0, numVertices);
    QtQorient.resize(numVertices, numVertices);
    orientedVertices.clear();
    orientationScale = 0.0;

    if (orientationWeight <= 0.0) return;
    if (!crossField) {
        std::cerr << "OASIS: orientation weight is set but there is no cross "
                     "field; the orientation term is off" << std::endl;
        return;
    }

    precomputeFits();
    const std::vector<double> areas = computeVertexAreas();

    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(static_cast<size_t>(numVertices) * 8);
    int noDirection = 0;

    for (int v = 0; v < numVertices; ++v) {
        if (!fits[v].valid) continue;

        const double mask = (orientationMask.size() == static_cast<Eigen::Index>(numVertices)) ? orientationMask[v] : 1.0;
        const double weight = areas[v] * mask;
        if (weight <= kEps) continue;

        double cos2t, sin2t;
        if (!crossFrame(v, cos2t, sin2t)) {
            ++noDirection;
            continue;
        }

        const VertexFit &fit = fits[v];
        const int row = static_cast<int>(orientedVertices.size());
        const double rowScale = std::sqrt(weight);

        for (size_t k = 0; k < fit.stencil.size(); ++k) {
            const Eigen::Index ki = static_cast<Eigen::Index>(k);
            const double q = cos2t * fit.auv[ki] - 0.5 * sin2t * (fit.auu[ki] - fit.avv[ki]);
            if (q != 0.0) trips.emplace_back(row, fit.stencil[k], rowScale * q);
        }
        orientedVertices.push_back(v);
    }

    if (orientedVertices.empty()) {
        std::cerr << "OASIS: the cross field yields no usable guiding direction; "
                     "the orientation term is off" << std::endl;
        return;
    }
    if (noDirection > 0) {
        // Singularities of the cross field, where the representation vector
        // vanishes and there is no direction to follow. Leaving those vertices
        // out is the right thing: it lets the singularity sit wherever the
        // Helmholtz term prefers.
        std::cerr << "OASIS: no guiding direction at " << noDirection
                  << " vertices (cross-field singularities); they carry no "
                     "orientation penalty" << std::endl;
    }

    Qorient.resize(static_cast<Eigen::Index>(orientedVertices.size()), numVertices);
    Qorient.setFromTriplets(trips.begin(), trips.end());
    Qorient.makeCompressed();

    QtQorient = Eigen::SparseMatrix<double>(Qorient.transpose()) * Qorient;
    QtQorient.makeCompressed();
}

double OASIS::orientationError() const {
    if (f.size() != static_cast<Eigen::Index>(mesh->vertices.size())) return -1.0;
    if (!crossField) return -1.0;
    precomputeFits();

    const int numVertices = static_cast<int>(mesh->vertices.size());
    const std::vector<double> areas = computeVertexAreas();

    // Per vertex: the deviatoric part of the Hessian in the guiding frame,
    //   q = 1/2 (h_uu - h_vv),  s = h_uv,
    // whose polar angle is twice the angle from the guiding frame to the
    // Hessian principal frame. R = sqrt(q^2 + s^2) is the anisotropy, and where
    // it vanishes the Hessian is umbilic and has no principal directions to
    // compare against.
    std::vector<double> deviation, weight, anisotropy;
    double maxR = 0.0;
    Eigen::VectorXd values;

    for (int v = 0; v < numVertices; ++v) {
        if (!fits[v].valid) continue;
        const double mask = (orientationMask.size() == static_cast<Eigen::Index>(numVertices)) ? orientationMask[v] : 1.0;
        if (areas[v] * mask <= kEps) continue;

        double cos2t, sin2t;
        if (!crossFrame(v, cos2t, sin2t)) continue;

        const VertexFit &fit = fits[v];
        values.resize(static_cast<Eigen::Index>(fit.stencil.size()));
        for (size_t k = 0; k < fit.stencil.size(); ++k)
            values[static_cast<Eigen::Index>(k)] = f[fit.stencil[k]];

        const double hxx = fit.auu.dot(values);
        const double hxy = fit.auv.dot(values);
        const double hyy = fit.avv.dot(values);

        const double qg = 0.5 * (hxx - hyy);
        const double q = cos2t * qg + sin2t * hxy;
        const double s = cos2t * hxy - sin2t * qg;

        // Half the polar angle, folded into [-45, 45] degrees: the two
        // principal directions are interchangeable, so a 90 degree offset is
        // the same alignment (and gives the same s = 0).
        double angle = 0.5 * std::atan2(s, q) * 180.0 / M_PI;
        angle -= 90.0 * std::round(angle / 90.0);

        deviation.push_back(std::fabs(angle));
        weight.push_back(areas[v] * mask);
        anisotropy.push_back(std::sqrt(q * q + s * s));
        maxR = std::max(maxR, anisotropy.back());
    }

    double total = 0.0, totalWeight = 0.0;
    for (size_t k = 0; k < deviation.size(); ++k) {
        if (anisotropy[k] <= 1e-6 * maxR) continue;  // umbilic, no frame to compare
        total += weight[k] * deviation[k];
        totalWeight += weight[k];
    }
    return (totalWeight > kEps) ? total / totalWeight : 0.0;
}

// Eq. 14 then Eq. 13.
void OASIS::buildKKT() {
    const int numVertices = static_cast<int>(mesh->vertices.size());
    const int nc2 = static_cast<int>(B.rows());

    // H_Lr = L_r^T L_r - lambda (L_r + L_r^T) + lambda^2 I
    //      = (L_r - lambda I)^T (L_r - lambda I)
    const Eigen::SparseMatrix<double> LrT = Eigen::SparseMatrix<double>(Lr.transpose());
    Eigen::SparseMatrix<double> identity(numVertices, numVertices);
    identity.setIdentity();

    HLr = LrT * Lr;
    HLr -= lambda * (Lr + LrT);
    HLr += (lambda * lambda) * identity;
    HLr.makeCompressed();

    // Sec. 5.1: minimizing E_Lr(f) + gamma E_Orient(f) under the same
    // constraints only changes the (1,1) block, since the orientation energy is
    // ||Q_orient f||^2 and its normal-equation block is Q^T Q.
    //
    // gamma is applied to a *normalized* Q^T Q rather than to the raw one. The
    // two terms are not commensurable as written: H_Lr carries entries of order
    // 1/h^4 while Q^T Q, being an area-weighted second derivative, carries
    // h^2/h^4 = 1/h^2, so the same gamma would mean something different on
    // every mesh and at every lambda. Matching their largest entries first
    // makes gamma a scale-free knob, at the price of the paper's number not
    // transferring.
    orientationScale = 0.0;
    if (Qorient.rows() > 0) {
        const double nH = maxAbsEntry(HLr);
        const double nQ = maxAbsEntry(QtQorient);
        if (nH > kEps && nQ > kEps) {
            orientationScale = orientationWeight * nH / nQ;
            HLr += orientationScale * QtQorient;
            HLr.makeCompressed();
        }
    }

    // [ H_Lr  B^T ] [ f  ]   [ 0 ]
    // [ B      0  ] [ nu ] = [ C ]
    const int n = numVertices + nc2;

    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(static_cast<size_t>(HLr.nonZeros()) + 2 * static_cast<size_t>(B.nonZeros()));

    for (int col = 0; col < HLr.outerSize(); ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(HLr, col); it; ++it) {
            trips.emplace_back(static_cast<int>(it.row()), static_cast<int>(it.col()), it.value());
        }
    }
    for (int col = 0; col < B.outerSize(); ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(B, col); it; ++it) {
            const int i = static_cast<int>(it.row());
            const int j = static_cast<int>(it.col());
            trips.emplace_back(numVertices + i, j, it.value());  // B
            trips.emplace_back(j, numVertices + i, it.value());  // B^T
        }
    }

    KKT.resize(n, n);
    KKT.setFromTriplets(trips.begin(), trips.end());
    KKT.makeCompressed();

    kktRhs.resize(n);
    kktRhs.head(numVertices).setZero();
    kktRhs.tail(nc2) = C;
}

// ─── Vibration enhancement, Sec. 3.4 ────────────────────────────────────────

// Cache the Eq. 9 fit at every vertex. Sec. 3.3 asks for a two-ring stencil on
// the boundary and a one-ring stencil in the interior; both widen on demand.
void OASIS::precomputeFits() const {
    const int numVertices = static_cast<int>(mesh->vertices.size());
    if (static_cast<int>(fits.size()) == numVertices) return;  // already cached

    fits.assign(numVertices, VertexFit{});

    std::vector<int> neighbors;
    Eigen::MatrixXd coeffs;
    int failed = 0;

    for (int v = 0; v < numVertices; ++v) {
        const int start = mesh->isBoundaryVertex[v] ? kBoundaryFitRings : 1;
        int ringsUsed = start;
        if (!fitLocalQuadratic(v, neighbors, coeffs, ringsUsed, start)) {
            ++failed;
            continue;
        }
        VertexFit &fit = fits[v];
        fit.stencil = neighbors;
        fit.au  = coeffs.row(AU).transpose();
        fit.av  = coeffs.row(AV).transpose();
        fit.auu = coeffs.row(AUU).transpose();
        fit.auv = coeffs.row(AUV).transpose();
        fit.avv = coeffs.row(AVV).transpose();
        fit.valid = true;
    }

    if (failed > 0) {
        std::cerr << "OASIS: no usable quadratic fit at " << failed << " of "
                  << numVertices << " vertices; they carry no vibration penalty"
                  << std::endl;
    }
}

// Eq. 16 with Eq. 17 expanded. W^u = a_u a_u^T/(-lambda) + (a_uu a_uu^T +
// a_uv a_uv^T)/lambda^2 is a sum of outer products, so f^T W^u f is a handful
// of squared dot products of the fit rows with f -- no per-vertex matrix.
//
// Writing g = grad f and h_** for the Hessian entries in the fit frame, and
// p = (h_uu + h_vv)/2, q = (h_uu - h_vv)/2, s = h_uv, R = sqrt(q^2 + s^2):
//
//   A_u^2 + A_v^2 = |g|^2/(-lambda) + (h_uu^2 + h_vv^2 + 2 h_uv^2)/lambda^2
//
// which is frame independent (both terms are traces). In the Hessian principal
// frame the Hessian eigenvalues are p +/- R and the gradient splits through
// e1 e1^T - e2 e2^T = H_dev/R, giving
//
//   A_u^2 - A_v^2 = [q(g_u^2 - g_v^2) + 2 s g_u g_v] / (R * -lambda)
//                 + 4 p R / lambda^2.
//
// With s = 0 and q > 0 this collapses to the frame-naive expression, as it
// must: that is the case where the fit frame already is the principal frame.
void OASIS::vibrationAnisotropy(int i, const Eigen::VectorXd &values,
                                double &sum, double &diff,
                                Eigen::VectorXd *dSum, Eigen::VectorXd *dDiff) const {
    const VertexFit &fit = fits[i];

    const double gu  = fit.au.dot(values);
    const double gv  = fit.av.dot(values);
    const double huu = fit.auu.dot(values);
    const double huv = fit.auv.dot(values);
    const double hvv = fit.avv.dot(values);

    const double invNegLambda = 1.0 / (-lambda);   // lambda < 0, so this is > 0
    const double invLambdaSq = 1.0 / (lambda * lambda);

    const double p = 0.5 * (huu + hvv);
    const double q = 0.5 * (huu - hvv);
    const double sdev = huv;

    sum = (gu * gu + gv * gv) * invNegLambda
        + (huu * huu + hvv * hvv + 2.0 * huv * huv) * invLambdaSq;

    // R = 0 is an umbilic Hessian: the principal directions are undefined and
    // the paper's local model degenerates. The deviatoric term vanishes with R
    // there, so floor R and let that term go to zero rather than blow up.
    const double R2 = q * q + sdev * sdev;
    const double R = std::sqrt(R2);
    const double scale = std::max(std::fabs(p), 1.0) * 1e-8;

    const double N = q * (gu * gu - gv * gv) + 2.0 * sdev * gu * gv;

    if (R <= scale) {
        diff = 0.0;
        if (dSum) {
            *dSum = 2.0 * invNegLambda * (gu * fit.au + gv * fit.av)
                  + 2.0 * invLambdaSq * (huu * fit.auu + hvv * fit.avv
                                         + 2.0 * huv * fit.auv);
        }
        if (dDiff) dDiff->setZero(values.size());
        return;
    }

    diff = N * invNegLambda / R + 4.0 * p * R * invLambdaSq;

    if (dSum) {
        *dSum = 2.0 * invNegLambda * (gu * fit.au + gv * fit.av)
              + 2.0 * invLambdaSq * (huu * fit.auu + hvv * fit.avv
                                     + 2.0 * huv * fit.auv);
    }
    if (dDiff) {
        const Eigen::VectorXd dp = 0.5 * (fit.auu + fit.avv);
        const Eigen::VectorXd dq = 0.5 * (fit.auu - fit.avv);
        const Eigen::VectorXd &ds = fit.auv;

        const Eigen::VectorXd dN = (gu * gu - gv * gv) * dq
                                 + q * 2.0 * (gu * fit.au - gv * fit.av)
                                 + 2.0 * gu * gv * ds
                                 + 2.0 * sdev * (gv * fit.au + gu * fit.av);
        const Eigen::VectorXd dR = (q * dq + sdev * ds) / R;

        *dDiff = invNegLambda * (dN * R - N * dR) / R2
               + 4.0 * invLambdaSq * (dp * R + p * dR);
    }
}

double OASIS::vibrationEnergy() const {
    if (f.size() != static_cast<Eigen::Index>(mesh->vertices.size())) return -1.0;
    precomputeFits();

    double total = 0.0;
    int counted = 0;
    Eigen::VectorXd values;

    for (int i = 0; i < static_cast<int>(fits.size()); ++i) {
        if (!fits[i].valid) continue;
        const auto &stencil = fits[i].stencil;
        values.resize(static_cast<Eigen::Index>(stencil.size()));
        for (size_t k = 0; k < stencil.size(); ++k) values[static_cast<Eigen::Index>(k)] = f[stencil[k]];

        double sum, diff;
        vibrationAnisotropy(i, values, sum, diff);
        if (sum <= kEps) continue;

        // Eq. 18 collapses: with Abar^2 = (Au^2 + Av^2)/2 the two bracketed
        // terms are equal and opposite, so E_a = 2 * ((Au^2 - Av^2)/(Au^2 + Av^2))^2.
        const double ratio = diff / sum;
        total += 2.0 * ratio * ratio;
        ++counted;
    }
    return (counted > 0) ? total / counted : 0.0;
}

void OASIS::enhanceVibration(int iterations) {
    if (f.size() != static_cast<Eigen::Index>(mesh->vertices.size())) {
        throw std::runtime_error("OASIS::enhanceVibration: call solve() first");
    }
    if (iterations <= 0) return;

    precomputeFits();

    const int numVertices = static_cast<int>(mesh->vertices.size());
    const Eigen::Index numRows = B.rows();

    // Vertices that actually contribute a vibration residual.
    std::vector<int> active;
    active.reserve(static_cast<size_t>(numVertices));
    for (int i = 0; i < numVertices; ++i)
        if (fits[i].valid) active.push_back(i);
    if (active.empty()) {
        std::cerr << "OASIS: no vertices carry a vibration penalty; nothing to do"
                  << std::endl;
        return;
    }

    Eigen::SparseMatrix<double> identity(numVertices, numVertices);
    identity.setIdentity();
    const Eigen::SparseMatrix<double> M = Lr - lambda * identity;  // L_r - lambda I
    const Eigen::SparseMatrix<double> Mt = Eigen::SparseMatrix<double>(M.transpose());
    const Eigen::SparseMatrix<double> Bt = Eigen::SparseMatrix<double>(B.transpose());
    const Eigen::SparseMatrix<double> BtB = Bt * B;

    // Builds the vibration residual vector and its Jacobian at the given f.
    Eigen::VectorXd values;
    Eigen::VectorXd dSum, dDiff;
    auto evaluate = [&](const Eigen::VectorXd &x, Eigen::VectorXd &res,
                        Eigen::SparseMatrix<double> *jac) {
        res.setZero(static_cast<Eigen::Index>(active.size()));
        std::vector<Eigen::Triplet<double>> trips;
        if (jac) trips.reserve(active.size() * 12);

        for (size_t a = 0; a < active.size(); ++a) {
            const int i = active[a];
            const auto &stencil = fits[i].stencil;
            values.resize(static_cast<Eigen::Index>(stencil.size()));
            for (size_t k = 0; k < stencil.size(); ++k)
                values[static_cast<Eigen::Index>(k)] = x[stencil[k]];

            double sum, diff;
            vibrationAnisotropy(i, values, sum, diff, jac ? &dSum : nullptr,
                                jac ? &dDiff : nullptr);
            if (sum <= kEps) continue;
            // E_a = r^2 with r = sqrt(2) * (Au^2 - Av^2) / (Au^2 + Av^2).
            const double kRoot2 = std::sqrt(2.0);
            res[static_cast<Eigen::Index>(a)] = kRoot2 * diff / sum;

            if (jac) {
                // d/df [ diff/sum ] = (sum * d(diff) - diff * d(sum)) / sum^2
                const double invSumSq = 1.0 / (sum * sum);
                for (size_t k = 0; k < stencil.size(); ++k) {
                    const double g = kRoot2 * invSumSq *
                                     (sum * dDiff[static_cast<Eigen::Index>(k)] -
                                      diff * dSum[static_cast<Eigen::Index>(k)]);
                    if (g != 0.0)
                        trips.emplace_back(static_cast<int>(a), stencil[k], g);
                }
            }
        }
        if (jac) {
            jac->resize(static_cast<Eigen::Index>(active.size()), numVertices);
            jac->setFromTriplets(trips.begin(), trips.end());
            jac->makeCompressed();
        }
    };

    // Eq. 19's three terms live on wildly different scales here (see the
    // header note on omega and xi). Balance them by *operator* magnitude, not
    // by residual norm: what matters in the normal equations is how big
    // w * G^T G and w * B^T B are next to H_Lr, and H_Lr carries entries of
    // order 1/h^4. Scaling by residual norms instead leaves the boundary term
    // orders of magnitude too weak, and the alignment conditions then drift
    // away entirely during the pass.
    Eigen::VectorXd vibRes;
    Eigen::SparseMatrix<double> G0;
    evaluate(f, vibRes, &G0);

    const Eigen::SparseMatrix<double> H0 = Mt * M;
    const Eigen::SparseMatrix<double> G0t = Eigen::SparseMatrix<double>(G0.transpose());
    const Eigen::SparseMatrix<double> GtG0 = G0t * G0;

    const double nH = maxAbsEntry(H0);
    const double nG = maxAbsEntry(GtG0);
    const double nB = (numRows > 0) ? maxAbsEntry(BtB) : 0.0;

    const double wVib = (nG > kEps) ? vibrationWeight * nH / nG : 0.0;
    const double wBnd = (nB > kEps) ? boundaryPenaltyWeight * nH / nB : 0.0;

    // The orientation term of Sec. 5.1 rides along at the weight assemble()
    // already normalized, so this pass minimizes the same combined energy the
    // KKT solve did. Left out, it would be free for the vibration term to trade
    // the orientation away; included, evening out the two local amplitudes has
    // to pay for any misalignment it introduces.
    const bool useOrient = orientationScale > 0.0 && Qorient.rows() > 0;

    auto totalEnergy = [&](const Eigen::VectorXd &x, const Eigen::VectorXd &res) {
        const double eHelm = (M * x).squaredNorm();
        const double eVib = wVib * res.squaredNorm();
        const double eBnd = (numRows > 0) ? wBnd * (B * x - C).squaredNorm() : 0.0;
        const double eOri = useOrient ? orientationScale * (Qorient * x).squaredNorm() : 0.0;
        return eHelm + eVib + eBnd + eOri;
    };

    double energy = totalEnergy(f, vibRes);
    const double startEnergy = energy;
    const double startVibration = vibrationEnergy();

    // Gauss-Newton with Levenberg damping. The paper uses plain Gauss-Newton
    // and CHOLMOD; the damping is here because J^T J picks up the squared
    // conditioning of L_r and a plain step occasionally overshoots into a
    // worse energy.
    double damping = 1e-8;
    Eigen::SparseMatrix<double> G;
    Eigen::VectorXd best = f;
    int accepted = 0;

    for (int iter = 0; iter < iterations; ++iter) {
        evaluate(f, vibRes, &G);

        const Eigen::SparseMatrix<double> Gt = Eigen::SparseMatrix<double>(G.transpose());
        Eigen::SparseMatrix<double> normal = Mt * M;
        normal += wVib * (Gt * G);
        if (numRows > 0) normal += wBnd * BtB;
        if (useOrient) normal += orientationScale * QtQorient;

        Eigen::VectorXd grad = Mt * (M * f) + wVib * (Gt * vibRes);
        if (numRows > 0) grad += wBnd * (Bt * (B * f - C));
        if (useOrient) grad += orientationScale * (QtQorient * f);

        bool stepTaken = false;
        for (int attempt = 0; attempt < 8; ++attempt) {
            Eigen::SparseMatrix<double> damped = normal;
            for (int d = 0; d < numVertices; ++d) damped.coeffRef(d, d) += damping;
            damped.makeCompressed();

            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
            ldlt.compute(damped);
            if (ldlt.info() != Eigen::Success) {
                damping *= 10.0;
                continue;
            }
            const Eigen::VectorXd delta = ldlt.solve(-grad);
            if (ldlt.info() != Eigen::Success || !delta.allFinite()) {
                damping *= 10.0;
                continue;
            }

            const Eigen::VectorXd candidate = f + delta;
            Eigen::VectorXd candRes;
            evaluate(candidate, candRes, nullptr);
            const double candEnergy = totalEnergy(candidate, candRes);

            if (candEnergy < energy) {
                f = candidate;
                energy = candEnergy;
                vibRes = candRes;
                damping = std::max(damping * 0.1, 1e-12);
                stepTaken = true;
                ++accepted;
                break;
            }
            damping *= 10.0;
        }
        if (!stepTaken) break;  // damping could not find a descent step
        best = f;
    }

    f = best;

    std::cerr << "OASIS: vibration enhancement took " << accepted << " of "
              << iterations << " steps; energy " << startEnergy << " -> " << energy
              << ", mean E_a " << startVibration << " -> " << vibrationEnergy()
              << std::endl;

    // The boundary conditions are a penalty in this phase, not hard
    // constraints, so they end up approximately rather than exactly satisfied.
    if (numRows > 0) {
        const double bcError = (B * f - C).cwiseAbs().maxCoeff();
        const double amplitude = std::max(1.0, f.cwiseAbs().maxCoeff());
        if (bcError > 1e-2 * amplitude) {
            std::cerr << "OASIS: WARNING - boundary conditions drifted to "
                      << bcError << " (relative " << bcError / amplitude
                      << ") during vibration enhancement" << std::endl;
        }
    }
}
