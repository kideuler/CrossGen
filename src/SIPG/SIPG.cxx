#include "SIPG.hxx"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_set>

// ---------------------------------------------------------------------------
void SIPG::setAlignedInteriorEdges(const std::vector<int> &edges) {
    edgeAligned.assign(mesh->edges.size(), 0);
    for (int e : edges) {
        if (e >= 0 && e < static_cast<int>(edgeAligned.size())) edgeAligned[e] = 1;
    }
}

// ---------------------------------------------------------------------------
// initialize()  --  Algorithm 1 from Section 10 of the paper
//
// Assembles the p=0 DG/SIPG matrices M, K and boundary RHS b, then
// factorises the system matrix A = M + tau*K and sets the initial field.
//
// DOFs are one complex number per triangle (not per vertex).
// ---------------------------------------------------------------------------
void SIPG::initialize() {
    const int NT = static_cast<int>(mesh->triangles.size());

    // -----------------------------------------------------------------------
    // Step 1 – Pre-compute triangle areas
    // -----------------------------------------------------------------------
    std::vector<double> area(NT, 0.0);
    for (int t = 0; t < NT; ++t) {
        const Triangle &tri = mesh->triangles[t];
        const Point &p0 = mesh->vertices[tri[0]];
        const Point &p1 = mesh->vertices[tri[1]];
        const Point &p2 = mesh->vertices[tri[2]];
        double x10 = p1[0] - p0[0], y10 = p1[1] - p0[1];
        double x20 = p2[0] - p0[0], y20 = p2[1] - p0[1];
        area[t] = 0.5 * std::fabs(x10 * y20 - x20 * y10);
    }

    // -----------------------------------------------------------------------
    // Step 2 – Assemble M (diagonal: M_ii = A_i)
    // -----------------------------------------------------------------------
    std::vector<Eigen::Triplet<std::complex<double>>> massTrips;
    massTrips.reserve(NT);
    for (int t = 0; t < NT; ++t) {
        massTrips.emplace_back(t, t, std::complex<double>(area[t], 0.0));
    }

    // -----------------------------------------------------------------------
    // Step 3 – Assemble K and b via edge loop
    //
    // For p=0 the volume gradient terms vanish; only edge penalty terms remain.
    // -----------------------------------------------------------------------
    std::vector<Eigen::Triplet<std::complex<double>>> stiffTrips;
    stiffTrips.reserve(mesh->edges.size() * 4);

    b.resize(NT);
    b.setZero();

    // --- Interior edges ---
    for (int edgeIdx = 0; edgeIdx < static_cast<int>(mesh->edges.size()); ++edgeIdx) {
        if (mesh->isBoundaryEdge[edgeIdx]) continue; // handled separately below
        if (isAlignedEdge(edgeIdx)) continue;        // ditto, one side at a time

        int ti = mesh->edgeTriangles[edgeIdx][0]; // K_i
        int tj = mesh->edgeTriangles[edgeIdx][1]; // K_j

        // Compute edge length
        const Point &ea = mesh->vertices[mesh->edges[edgeIdx][0]];
        const Point &eb = mesh->vertices[mesh->edges[edgeIdx][1]];
        double dx = eb[0] - ea[0], dy = eb[1] - ea[1];
        double edgeLen = std::sqrt(dx * dx + dy * dy);
        if (edgeLen < 1e-14) continue;

        // Normal heights: h_i^e = 2*A_i / |e|
        double hi = 2.0 * area[ti] / edgeLen;
        double hj = 2.0 * area[tj] / edgeLen;
        double he = std::min(hi, hj);

        // Penalty weight: kappa_e = gamma * |e| / h_e
        double kappa = gamma * edgeLen / he;

        // Planar connection factor P_e = 1 (one global frame)
        // K_ii += kappa,  K_ij -= kappa,  K_ji -= kappa,  K_jj += kappa
        stiffTrips.emplace_back(ti, ti, std::complex<double>( kappa, 0.0));
        stiffTrips.emplace_back(ti, tj, std::complex<double>(-kappa, 0.0));
        stiffTrips.emplace_back(tj, ti, std::complex<double>(-kappa, 0.0));
        stiffTrips.emplace_back(tj, tj, std::complex<double>( kappa, 0.0));
    }

    // --- Boundary edges, and the interior edges the caller asked to align to ---
    //
    // Both are one-sided Dirichlet data of the same kind: a curve the field has
    // to be tangent to, seen from the triangle on one side of it. A boundary
    // edge has one such triangle and an aligned interior edge has two, and the
    // only difference in the assembly is that the second one is also visited.
    // Also accumulate weighted-average BC value per triangle for hard pinning.
    std::unordered_map<int, std::complex<double>> bcWeightedSum; // numerator: sum(kappa * ge)
    std::unordered_map<int, double>               bcWeightTotal; // denominator: sum(kappa)

    // Pin triangle `ti` to the tangent of `edgeIdx`. A triangle can be pinned by
    // more than one edge -- a corner of dS, a triangle straddling a junction --
    // and the pins are then averaged, which is what makes exp(4i theta) the
    // right thing to average: two edges that meet at a right angle carry the
    // same value and reinforce, and only edges genuinely off the same cross
    // pull against each other.
    //
    // The averaging is only *accumulated* here; what it becomes in K and b is
    // decided below, once every edge incident on the triangle has been seen,
    // because the strength of the pin has to depend on whether the edges agree.
    auto pinToEdge = [&](int edgeIdx, int ti) {
        if (ti < 0) return;

        // Edge endpoints (stored as (min, max) vertex indices, but we need oriented order)
        // Use the triangle's local ordering to get CCW orientation
        int localEdge = -1;
        for (int e = 0; e < 3; ++e) {
            if (mesh->triangleEdges[ti][e] == edgeIdx) {
                localEdge = e;
                break;
            }
        }
        if (localEdge < 0) return;
        const Triangle &tri = mesh->triangles[ti];
        int va = tri[localEdge];
        int vb = tri[(localEdge + 1) % 3];

        const Point &pa = mesh->vertices[va];
        const Point &pb = mesh->vertices[vb];
        double dx = pb[0] - pa[0], dy = pb[1] - pa[1];
        double edgeLen = std::sqrt(dx * dx + dy * dy);
        if (edgeLen < 1e-14) return;

        // Normal height and penalty
        double hi = 2.0 * area[ti] / edgeLen;
        if (hi < 1e-14) return;
        double kappa = gamma * edgeLen / hi;

        // Tangent angle theta_e and spin-4 value g_e = exp(4i * theta_e). The
        // orientation of the edge does not enter: exp(4i theta) is what a cross
        // is, and a cross has no head.
        double theta = std::atan2(dy, dx);
        std::complex<double> ge = std::exp(std::complex<double>(0.0, 4.0 * theta));

        bcWeightedSum[ti] += kappa * ge;
        bcWeightTotal[ti] += kappa;
    };

    for (int edgeIdx : mesh->boundaryEdges) {
        pinToEdge(edgeIdx, mesh->edgeTriangles[edgeIdx][0]); // the only triangle on it
    }
    for (int edgeIdx = 0; edgeIdx < static_cast<int>(edgeAligned.size()); ++edgeIdx) {
        if (!edgeAligned[edgeIdx] || mesh->isBoundaryEdge[edgeIdx]) continue;
        pinToEdge(edgeIdx, mesh->edgeTriangles[edgeIdx][0]);
        pinToEdge(edgeIdx, mesh->edgeTriangles[edgeIdx][1]);
    }

    // --- Disk centers ------------------------------------------------------
    //
    // One more Dirichlet triangle per disk, pinned to u = 1, which is what
    // fixes the rotation the disk's symmetry otherwise leaves free -- see
    // SIPG::setPinDiskCenters for why that value and where the cones then go.
    // It is the same kind of pin as a boundary one, just prescribed by a point
    // rather than read off an edge tangent, so it goes through the same
    // bcWeightedSum/K/b path and gets eliminated with the rest in Step 6.
    diskCenterTriangles.clear();
    if (pinDiskCenters) {
        for (const auto &comp : mesh->materialComponents) {
            if (!comp.circle.isCircle) continue;
            int tc = comp.centerTriangle;
            if (tc < 0) continue;

            // Already Dirichlet from dS or an interface. That only happens on a
            // disk so coarse its center triangle touches its own boundary, and
            // there the tangent data is the better constraint of the two.
            if (bcWeightTotal.count(tc)) continue;

            // Weight to match an edge pin's: gamma * |e| / h with |e| the
            // triangle's longest edge and h = 2A/|e| its height on that edge.
            const Triangle &tri = mesh->triangles[tc];
            double eMax = 0.0;
            for (int k = 0; k < 3; ++k) {
                eMax = std::max(eMax, normP(mesh->vertices[tri[(k + 1) % 3]] - mesh->vertices[tri[k]]));
            }
            if (area[tc] < 1e-30 || eMax < 1e-14) continue;
            double kappa = gamma * eMax * eMax / (2.0 * area[tc]);

            const std::complex<double> gc(1.0, 0.0); // exp(4i*0): cross on the axes
            bcWeightedSum[tc] += kappa * gc;
            bcWeightTotal[tc] += kappa;
            diskCenterTriangles.push_back(tc);
        }
    }

    // --- Emit the Dirichlet data into K and b, and decide what is pinnable ---
    //
    // A triangle with a single aligned edge is the ordinary case: one tangent,
    // nothing to reconcile, and the pin is that tangent at full strength.
    //
    // A triangle carrying two of them is not, and the naive average is wrong in
    // a way that is worth spelling out. The data is exp(4i theta) per edge, so
    // two edges a multiple of 90 degrees apart carry the *same* value and
    // reinforce -- that is the whole reason the spin-4 embedding is the right
    // thing to average in. Two edges 45 degrees apart carry exp(4i*0) = +1 and
    // exp(4i*pi/4) = -1 and cancel outright. The sum is then zero, its direction
    // is whatever the arithmetic noise left, and the assembled system asks the
    // triangle to be equal to nothing at all at full penalty -- it is pulled to
    // the origin and normalisation turns that into noise. Measured at the apex
    // of a 45-degree wedge: 41.8 degrees of boundary misalignment, on the one
    // element where alignment is most visible.
    //
    // The honest reading is that a cross cannot be tangent to both edges, that
    // the constraint set is inconsistent, and that the element is too coarse to
    // resolve the turn. So weight the pin by how consistent its data actually is,
    //
    //     r = |sum_e kappa_e g_e| / sum_e kappa_e  in [0, 1],
    //
    // and impose  r * (sum kappa_e) * ghat  with ghat the unit average. Since
    // r * (sum kappa) * ghat = sum kappa_e g_e, the right-hand side is exactly
    // what it always was and only the diagonal changes: the penalty's minimiser
    // becomes the *unit* vector ghat instead of the short vector r*ghat, so the
    // pull toward the origin is gone. r = 1 at a consistent corner reproduces
    // the old assembly exactly; r -> 0 at a 45-degree one withdraws the
    // constraint and lets the triangle take the value its neighbours imply,
    // which is the smoothest boundary-compatible answer available on this mesh.
    //
    // Below kCornerCoherenceMin the averaged direction carries no information
    // worth pinning to, so such a triangle is left out of the hard-pin set and
    // diffuses instead. Everything else is pinned exactly as before.
    boundaryTriangleBC.clear();
    for (const auto &[ti, wsum] : bcWeightedSum) {
        const double wtot = bcWeightTotal[ti];
        if (wtot < 1e-30) continue;
        const double mag = std::abs(wsum);
        const double r = cornerCoherenceFix ? mag / wtot : 1.0;

        stiffTrips.emplace_back(ti, ti, std::complex<double>(r * wtot, 0.0));
        b[ti] += wsum;

        if (mag > 1e-14 && r >= (cornerCoherenceFix ? kCornerCoherenceMin : 0.0))
            boundaryTriangleBC[ti] = wsum / mag;
    }

    // -----------------------------------------------------------------------
    // Step 4 – Build sparse matrices
    // -----------------------------------------------------------------------
    M.resize(NT, NT);
    K.resize(NT, NT);
    M.setFromTriplets(massTrips.begin(), massTrips.end());
    K.setFromTriplets(stiffTrips.begin(), stiffTrips.end());

    // -----------------------------------------------------------------------
    // Step 5 – Choose tau = D^2 / 10 (same heuristic as CrossField)
    // -----------------------------------------------------------------------
    double minX = std::numeric_limits<double>::max(),  maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max(),  maxY = std::numeric_limits<double>::lowest();
    for (const auto &p : mesh->vertices) {
        minX = std::min(minX, p[0]); maxX = std::max(maxX, p[0]);
        minY = std::min(minY, p[1]); maxY = std::max(maxY, p[1]);
    }
    double D = std::sqrt((maxX - minX) * (maxX - minX) + (maxY - minY) * (maxY - minY));
    tau = tauScale * D * D / 10.0;

    // -----------------------------------------------------------------------
    // Step 6 – Form A = M + tau*K, then eliminate boundary triangle DOFs
    //
    // Following CrossField's approach: for each boundary triangle row, zero out
    // the entire row in both M and A, then put 1 on the diagonal.  This turns
    // the boundary triangle equation into "u_tilde[ti] = u_k_prev[ti] = BC",
    // while interior rows retain full coupling to the (fixed) BC values through
    // the off-diagonal K columns – exactly as CrossField handles vertex BCs.
    // -----------------------------------------------------------------------
    A = M + tau * K;

    const std::complex<double> czero(0.0, 0.0);
    const std::complex<double> cone (1.0, 0.0);

    // Build a fast lookup for pinned triangles (on dS or on an aligned edge).
    //
    // Empty when the alignment is imposed weakly: A and b are then left as the
    // Nitsche form assembled them and no row is eliminated, which is exactly
    // the ablation setHardBoundaryConditions describes.
    std::unordered_set<int> bndTriSet;
    if (hardBoundary) {
        for (const auto &[ti, _] : boundaryTriangleBC) bndTriSet.insert(ti);
    }

    // Zero boundary rows in M (column-major iteration)
    for (int col = 0; col < M.outerSize(); ++col) {
        for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(M, col); it; ++it) {
            if (bndTriSet.count(it.row())) {
                it.valueRef() = (it.row() == col) ? cone : czero;
            }
        }
    }
    M.prune(czero);

    // Zero boundary rows in A (column-major iteration)
    for (int col = 0; col < A.outerSize(); ++col) {
        for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(A, col); it; ++it) {
            if (bndTriSet.count(it.row())) {
                it.valueRef() = (it.row() == col) ? cone : czero;
            }
        }
    }
    A.prune(czero);

    // -----------------------------------------------------------------------
    // Step 7 – Factorise the modified A
    // -----------------------------------------------------------------------
    solverLU.compute(A);
    if (solverLU.info() != Eigen::Success) {
        std::cerr << "SIPG: SparseLU factorization failed, falling back to BiCGSTAB" << std::endl;
        useBiCGSTAB = true;
        solverBiCGSTAB.setTolerance(1e-10);
        solverBiCGSTAB.setMaxIterations(1000);
        solverBiCGSTAB.compute(A);
        if (solverBiCGSTAB.info() != Eigen::Success) {
            throw std::runtime_error("SIPG: Both SparseLU and BiCGSTAB factorization failed");
        }
    } else {
        useBiCGSTAB = false;
    }

    // -----------------------------------------------------------------------
    // Step 8 – Initialise field: interior = 1, boundary = BC value
    // -----------------------------------------------------------------------
    u_k_prev.resize(NT);
    u_k_prev.setOnes();

    for (const auto &[ti, bc] : boundaryTriangleBC) {
        u_k_prev[ti] = bc;
    }

    u_k.resize(NT);
    u_k = u_k_prev;
}

// ---------------------------------------------------------------------------
// step()  --  one MBO iteration (Algorithm 2 from Section 10)
// ---------------------------------------------------------------------------
void SIPG::step() {
    // Form RHS: r^k = M * u^k + tau * b
    Eigen::VectorXcd rhs = M * u_k_prev + tau * b;

    // Solve A * u_tilde = r^k
    Eigen::VectorXcd u_tilde;
    if (useBiCGSTAB) {
        u_tilde = solverBiCGSTAB.solve(rhs);
        if (solverBiCGSTAB.info() != Eigen::Success) {
            throw std::runtime_error("SIPG: BiCGSTAB solve failed");
        }
    } else {
        u_tilde = solverLU.solve(rhs);
        if (solverLU.info() != Eigen::Success) {
            throw std::runtime_error("SIPG: SparseLU solve failed");
        }
    }

    // Project each triangle's value back onto the unit circle
    u_k.resize(u_tilde.size());
    for (int i = 0; i < u_tilde.size(); ++i) {
        double mag = std::abs(u_tilde[i]);
        if (mag > 1e-14) {
            u_k[i] = u_tilde[i] / mag;
        } else {
            u_k[i] = u_k_prev[i]; // keep previous if magnitude is zero
        }
    }

    // Update error for convergence checking
    error = (u_k - u_k_prev).norm();

    // Advance
    u_k_prev = u_k;
}

// ---------------------------------------------------------------------------
// computeSingularities()
//
// For p=0 DG each triangle carries a single complex value, so singularities
// are detected by looking at the holonomy around each interior vertex:
// sum the angle differences of the four adjacent triangle values and check
// for non-zero winding number.
// ---------------------------------------------------------------------------
void SIPG::computeSingularities() {
    singularVertices.clear();

    // Helper: smallest-angle difference between two unit complex numbers in (-pi, pi]
    auto angleDiff = [](std::complex<double> z_start, std::complex<double> z_end) -> double {
        std::complex<double> d = z_end * std::conj(z_start);
        return std::atan2(d.imag(), d.real());
    };

    // For each interior vertex, sum angle differences around the triangle star.
    // A non-zero winding number indicates a singularity; we report the vertex index.
    const int NV = static_cast<int>(mesh->vertices.size());
    for (int v = 0; v < NV; ++v) {
        if (mesh->isBoundaryVertex[v]) continue;

        const auto &vt = mesh->vertexTriangles;
        int start = vt.rowPtr[v];
        int end   = vt.rowPtr[v + 1];
        int star  = end - start;
        if (star < 2) continue;

        double totalAngle = 0.0;
        for (int k = start; k < end; ++k) {
            int tCur  = vt.colIdx[k];
            int tNext = vt.colIdx[(k - start + 1) % star + start];
            totalAngle += angleDiff(u_k_prev[tCur], u_k_prev[tNext]);
        }

        int winding = static_cast<int>(std::round(totalAngle / (2.0 * M_PI)));
        if (winding != 0) {
            double crossIndex = winding / 4.0;
            singularVertices.emplace_back(v, crossIndex);
        }
    }
}

// ---------------------------------------------------------------------------
// runMBO()
// ---------------------------------------------------------------------------
void SIPG::runMBO() {
    int iteration = 0;
    double ntris = static_cast<double>(mesh->triangles.size());
    while (iteration < maxIterations && error > 2.0 * ntris * 1e-7) {
        step();
        ++iteration;
    }
}

