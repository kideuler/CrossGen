#include "OASIS.hxx"

#include <algorithm>
#include <iostream>
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

void OASIS::assemble() {
    if (lambda >= 0.0) {
        std::cerr << "OASIS: lambda should be negative (Eq. 1), got " << lambda << std::endl;
    }
    buildDensityLaplacian();
    buildConstraints();
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

    // Total incident area per vertex.
    std::vector<double> vertexArea(numVertices, 0.0);
    for (int t = 0; t < numTriangles; ++t) {
        const Triangle &tri = mesh->triangles[t];
        const Point &p0 = mesh->vertices[tri[0]];
        const Point &p1 = mesh->vertices[tri[1]];
        const Point &p2 = mesh->vertices[tri[2]];

        const double area = 0.5 * std::fabs(cross2(p1 - p0, p2 - p0));
        if (area < kEps) continue;

        for (int k = 0; k < 3; ++k) {
            vertexArea[tri[k]] += area;
        }
    }

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

    // Row scaling 3 / (2 * Area_i * r_i).
    Eigen::VectorXd rowScale(numVertices);
    for (int v = 0; v < numVertices; ++v) {
        const double denom = 2.0 * vertexArea[v] * r[v];
        rowScale[v] = (denom > kEps) ? (3.0 / denom) : 0.0;
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
                              int &ringsUsed) const {
    for (int rings = kBoundaryFitRings; rings <= kMaxFitRings; ++rings) {
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
        if (!fitLocalQuadratic(v, neighbors, coeffs, ringsUsed)) {
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
