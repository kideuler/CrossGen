#include "SIPG.hxx"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_set>

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

    // --- Boundary edges ---
    for (int edgeIdx : mesh->boundaryEdges) {
        int ti = mesh->edgeTriangles[edgeIdx][0]; // K_i (only triangle on this edge)

        // Edge endpoints (stored as (min, max) vertex indices, but we need oriented order)
        // Use the triangle's local ordering to get CCW orientation
        int localEdge = -1;
        for (int e = 0; e < 3; ++e) {
            if (mesh->triangleEdges[ti][e] == edgeIdx) {
                localEdge = e;
                break;
            }
        }
        const Triangle &tri = mesh->triangles[ti];
        int va = tri[localEdge];
        int vb = tri[(localEdge + 1) % 3];

        const Point &pa = mesh->vertices[va];
        const Point &pb = mesh->vertices[vb];
        double dx = pb[0] - pa[0], dy = pb[1] - pa[1];
        double edgeLen = std::sqrt(dx * dx + dy * dy);
        if (edgeLen < 1e-14) continue;

        // Normal height and penalty
        double hi = 2.0 * area[ti] / edgeLen;
        double kappa = gamma * edgeLen / hi;

        // Boundary tangent angle theta_e and spin-4 value g_e = exp(4i * theta_e)
        double theta = std::atan2(dy, dx);
        std::complex<double> ge = std::exp(std::complex<double>(0.0, 4.0 * theta));

        // K_ii += kappa,  b_i += kappa * g_e
        stiffTrips.emplace_back(ti, ti, std::complex<double>(kappa, 0.0));
        b[ti] += kappa * ge;
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
    tau = D * D / 10.0;

    // -----------------------------------------------------------------------
    // Step 6 – Form and factorize A = M + tau*K
    // -----------------------------------------------------------------------
    A = M + tau * K;

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
    // Step 7 – Initialize field u^0 to the boundary-driven values where
    //          boundary edges exist, and ones elsewhere
    // -----------------------------------------------------------------------
    u_k_prev.resize(NT);
    u_k_prev.setOnes(); // start from constant unit field

    // Warm-start: one diffusion step to propagate boundary data into interior
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

