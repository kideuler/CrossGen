#include "CrossField.hxx"
#include <cstdlib>
#include <ctime>
#include <limits>
#include <iostream>

void CrossField::setAlignedInteriorEdges(const std::vector<int> &edges) {
    edgeAligned.assign(mesh->edges.size(), 0);
    for (int e : edges)
        if (e >= 0 && e < static_cast<int>(edgeAligned.size()) && !mesh->isBoundaryEdge[e])
            edgeAligned[e] = 1;
}

void CrossField::initialize(int method, unsigned seed) {
    // For all boundary vertices, set dirichlet boundary conditions to be (nx + i*ny)^4 where (nx, ny) is the outward normal of the boundary edge

    int numVertices = static_cast<int>(mesh->vertices.size());
    u_k_prev.resize(numVertices);
    u_k_prev.setZero();

    // For each boundary vertex, compute weighted average normal from adjacent boundary edges
    // The normal at a node is the weighted average of adjacent boundary edge normals, weighted by edge length

    // First, collect all boundary edges (vertex pairs) from boundaryTriangles and cornerTriangles
    // Each boundary edge contributes to its two endpoint vertices

    // Map: vertex -> accumulated (weighted normal, total weight)
    std::vector<double> normalX(numVertices, 0.0);
    std::vector<double> normalY(numVertices, 0.0);
    std::vector<double> weights(numVertices, 0.0);

    // Process all boundary edges using the new edge data structures
    // For boundary edges, we need to get the correct orientation from the triangle
    for (int edgeIdx : mesh->boundaryEdges) {
        // Get the triangle that owns this boundary edge (first one, since second is -1)
        int triIdx = mesh->edgeTriangles[edgeIdx][0];
        const Triangle &tri = mesh->triangles[triIdx];
        
        // Find which local edge (0, 1, or 2) this is in the triangle
        int localEdge = -1;
        for (int e = 0; e < 3; ++e) {
            if (mesh->triangleEdges[triIdx][e] == edgeIdx) {
                localEdge = e;
                break;
            }
        }
        
        // Get vertices in CCW order from the triangle
        int v0 = tri[localEdge];
        int v1 = tri[(localEdge + 1) % 3];

        const Point &p0 = mesh->vertices[v0];
        const Point &p1 = mesh->vertices[v1];

        // Edge vector from v0 to v1 (in CCW order)
        double ex = p1[0] - p0[0];
        double ey = p1[1] - p0[1];

        // Edge length
        double len = std::sqrt(ex * ex + ey * ey);
        if (len < 1e-14) continue;

        // Outward normal: rotate edge 90 degrees clockwise (for CCW-oriented triangles)
        // For CCW triangle with edge going v0->v1, outward normal is (ey, -ex) normalized
        double nx = ey / len;
        double ny = -ex / len;

        // Accumulate weighted normal for both vertices
        normalX[v0] += nx * len;
        normalY[v0] += ny * len;
        weights[v0] += len;

        normalX[v1] += nx * len;
        normalY[v1] += ny * len;
        weights[v1] += len;
    }

    // The interior angle at each boundary vertex: the tip angles of the
    // triangles meeting there, added up. This is what Table 1 of Viertel,
    // Osting and Staten (IMR 2019) classifies the vertex by.
    std::vector<double> interiorAngle(numVertices, 0.0);
    for (const Triangle &tri : mesh->triangles) {
        for (int k = 0; k < 3; ++k) {
            const int v = tri[k];
            if (!mesh->isBoundaryVertex[v]) continue;
            const Point &a = mesh->vertices[v];
            const Point u = mesh->vertices[tri[(k + 1) % 3]] - a;
            const Point w = mesh->vertices[tri[(k + 2) % 3]] - a;
            const double nu = normP(u), nw = normP(w);
            if (nu < 1e-15 || nw < 1e-15) continue;
            interiorAngle[v] += std::acos(std::max(-1.0, std::min(1.0, dotP(u, w) / (nu * nw))));
        }
    }

    // Now set Dirichlet BCs for boundary vertices, Sec. 2.2 and Table 1.
    //
    // d is the bisector of the outward normals of the boundary edges meeting
    // at the vertex. On a smooth stretch of boundary the two normals agree, so
    // d is the normal and a cross aligned to it has an arm along the boundary,
    // which is what boundary alignment means. At a corner they do not: at a
    // ninety degree one they are ninety degrees apart, so d sits forty-five
    // degrees from each, and a cross aligned to d has its arms forty-five
    // degrees off *both* boundary edges -- the worst possible alignment,
    // exactly where alignment matters most. Turning d by a further π/4 at a
    // vertex whose index is ±1/4 puts the arms back along the two edges.
    //
    // Since u = d^4, that quarter turn is just a change of sign.
    for (int v : mesh->boundaryVertices) {
        if (weights[v] > 1e-14) {
            // Compute normalized weighted average normal
            double nx = normalX[v] / weights[v];
            double ny = normalY[v] / weights[v];

            // Normalize the averaged normal
            double nlen = std::sqrt(nx * nx + ny * ny);
            if (nlen > 1e-14) {
                nx /= nlen;
                ny /= nlen;
            }

            const double a = interiorAngle[v];
            int quarters;                                   // the index, times four
            if (a < 0.75 * M_PI)       quarters =  1;       // a convex corner
            else if (a <= 1.25 * M_PI) quarters =  0;       // effectively straight
            else if (a <= 1.75 * M_PI) quarters = -1;       // reflex
            else                       quarters = -2;

            std::complex<double> n(nx, ny);
            // normalize
            n /= std::abs(n);
            std::complex<double> u = std::pow(n, 4);
            if (quarters == 1 || quarters == -1) u = -u;
            u_k_prev[v] = u;
        }
    }

    // initialize stiffness and mass matrices for MBO method
    // We will solve (M + tau*K) u^{k+1} = M u^k where M is the mass matrix and K is the stiffness matrix.

    int numTriangles = static_cast<int>(mesh->triangles.size());

    // Build triplets for sparse matrix construction
    std::vector<Eigen::Triplet<std::complex<double>>> massTrips;
    std::vector<Eigen::Triplet<std::complex<double>>> stiffTrips;
    massTrips.reserve(numTriangles * 9);  // 3x3 per triangle
    stiffTrips.reserve(numTriangles * 9);

    for (int t = 0; t < numTriangles; ++t) {
        const Triangle &tri = mesh->triangles[t];
        int v0 = tri[0], v1 = tri[1], v2 = tri[2];

        // Get vertex coordinates
        const Point &p0 = mesh->vertices[v0];
        const Point &p1 = mesh->vertices[v1];
        const Point &p2 = mesh->vertices[v2];

        // Compute edge vectors
        double x10 = p1[0] - p0[0], y10 = p1[1] - p0[1];
        double x20 = p2[0] - p0[0], y20 = p2[1] - p0[1];

        // Triangle area (2 * area = |cross product|)
        double area2 = x10 * y20 - x20 * y10;
        double area = 0.5 * std::fabs(area2);

        if (area < 1e-14) continue;

        // --- Mass matrix (consistent mass) ---
        // For P1 elements: M_ij = (area/12) * (1 + delta_ij)
        // Diagonal: area/6, Off-diagonal: area/12
        std::complex<double> diagMass(area / 6.0, 0.0);
        std::complex<double> offMass(area / 12.0, 0.0);

        int verts[3] = {v0, v1, v2};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (i == j) {
                    massTrips.emplace_back(verts[i], verts[j], diagMass);
                } else {
                    massTrips.emplace_back(verts[i], verts[j], offMass);
                }
            }
        }

        // --- Stiffness matrix ---
        // Gradient of basis functions for P1 elements:
        // grad(phi_i) = (1 / (2*area)) * [y_{i+1} - y_{i+2}, x_{i+2} - x_{i+1}]
        // where indices are cyclic (0,1,2)

        double invArea2 = 1.0 / area2;

        // Compute gradients of basis functions (scaled by 2*area, then divide)
        // grad(phi_0) = [y1 - y2, x2 - x1] / (2*area)
        // grad(phi_1) = [y2 - y0, x0 - x2] / (2*area)
        // grad(phi_2) = [y0 - y1, x1 - x0] / (2*area)

        double gradPhi[3][2];
        gradPhi[0][0] = (p1[1] - p2[1]) * invArea2;
        gradPhi[0][1] = (p2[0] - p1[0]) * invArea2;
        gradPhi[1][0] = (p2[1] - p0[1]) * invArea2;
        gradPhi[1][1] = (p0[0] - p2[0]) * invArea2;
        gradPhi[2][0] = (p0[1] - p1[1]) * invArea2;
        gradPhi[2][1] = (p1[0] - p0[0]) * invArea2;

        // K_ij = area * (grad(phi_i) . grad(phi_j))
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double dot = gradPhi[i][0] * gradPhi[j][0] + gradPhi[i][1] * gradPhi[j][1];
                std::complex<double> val(area * dot, 0.0);
                stiffTrips.emplace_back(verts[i], verts[j], val);
            }
        }
    }

    // A vertex no triangle references has an empty row in M and K, so A is
    // structurally singular and the direct factorisation is refused -- nine of
    // the 24 corpus models carry such vertices (geom031 has 309). Give each one
    // a unit mass entry so the row exists, then pin it below like a boundary
    // vertex: identity row, value 1, never read by anything since no element
    // touches it. This changes nothing about the field and lets the direct
    // solver do its job.
    std::vector<char> usedVertex(numVertices, 0);
    for (const Triangle &tri : mesh->triangles) { usedVertex[tri[0]] = 1; usedVertex[tri[1]] = 1; usedVertex[tri[2]] = 1; }
    for (int v = 0; v < numVertices; ++v)
        if (!usedVertex[v]) massTrips.emplace_back(v, v, std::complex<double>(1.0, 0.0));

    // Assemble sparse matrices
    M.resize(numVertices, numVertices);
    K.resize(numVertices, numVertices);
    M.setFromTriplets(massTrips.begin(), massTrips.end());
    K.setFromTriplets(stiffTrips.begin(), stiffTrips.end());

    // Compute tau = D^2 / 10 where D is the diameter of the domain (diagonal of bounding box)
    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();
    for (const auto &p : mesh->vertices) {
        minX = std::min(minX, p[0]);
        maxX = std::max(maxX, p[0]);
        minY = std::min(minY, p[1]);
        maxY = std::max(maxY, p[1]);
    }
    double D = std::sqrt((maxX - minX) * (maxX - minX) + (maxY - minY) * (maxY - minY));
    tau = tauScale * D * D / 10.0;

    // Compute system matrix A = M + tau * K
    A = M + tau * K;

    // Apply Dirichlet boundary conditions to M and A:
    // For each boundary vertex, set diagonal to 1 and zero out the rest of the row
    // Note: Eigen SparseMatrix is column-major, so we need to iterate carefully
    std::complex<double> one(1.0, 0.0);
    std::complex<double> zero(0.0, 0.0);
    
    // Create a set for fast boundary lookup
    std::unordered_set<int> boundarySet(mesh->boundaryVertices.begin(), mesh->boundaryVertices.end());

    // The unreferenced vertices, pinned (see the mass assembly above).
    for (int v = 0; v < numVertices; ++v) {
        if (usedVertex[v]) continue;
        boundarySet.insert(v);
        u_k_prev[v] = one;
    }

    // The interface vertices, pinned to the interface's own tangent (see
    // setAlignedInteriorEdges). exp(4i phi) is the same for a direction and
    // its reverse and for a direction and its quarter turn, so an edge's
    // orientation does not matter and a right-angled corner averages to
    // itself.
    alignedInterface = freeInterface = 0;
    if (!edgeAligned.empty()) {
        std::vector<std::complex<double>> sum(numVertices, zero);
        std::vector<double> weight(numVertices, 0.0);
        for (int e = 0; e < static_cast<int>(edgeAligned.size()); ++e) {
            if (!edgeAligned[e]) continue;
            const int a = mesh->edges[e][0], b = mesh->edges[e][1];
            const Point d = mesh->vertices[b] - mesh->vertices[a];
            const double len = normP(d);
            if (len < 1e-14) continue;
            const std::complex<double> g = std::polar(len, 4.0 * std::atan2(d[1], d[0]));
            for (const int v : {a, b}) { sum[v] += g; weight[v] += len; }
        }
        // DualMBO::kCornerCoherenceMin, for the same reason: below it the
        // mean is a direction no incident interface has.
        constexpr double kCoherenceMin = 0.2;
        for (int v = 0; v < numVertices; ++v) {
            if (!(weight[v] > 0.0) || mesh->isBoundaryVertex[v] || boundarySet.count(v)) continue;
            if (std::abs(sum[v]) < kCoherenceMin * weight[v]) { ++freeInterface; continue; }
            boundarySet.insert(v);
            u_k_prev[v] = sum[v] / std::abs(sum[v]);
            ++alignedInterface;
        }
    }

    // The disk centres (see setPinDiskCenters): one vertex each, to the cross
    // on the axes, unless it is already data.
    pinnedDisks = 0;
    if (pinDiskCenters) {
        for (const auto &comp : mesh->materialComponents) {
            if (!comp.circle.isCircle || comp.centerTriangle < 0) continue;
            const Triangle &t = mesh->triangles[comp.centerTriangle];
            int best = -1;
            double bestD = std::numeric_limits<double>::max();
            for (int k = 0; k < 3; ++k) {
                const double d = normP(mesh->vertices[t[k]] - comp.circle.center);
                if (d < bestD) { bestD = d; best = t[k]; }
            }
            if (best < 0 || boundarySet.count(best)) continue;
            boundarySet.insert(best);
            u_k_prev[best] = one;
            ++pinnedDisks;
        }
    }

    // For column-major matrices, we iterate over all columns and check each entry's row
    // Zero out rows for boundary vertices in M
    for (int col = 0; col < M.outerSize(); ++col) {
        for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(M, col); it; ++it) {
            int row = it.row();
            if (boundarySet.count(row)) {
                it.valueRef() = (row == col) ? one : zero;
            }
        }
    }
    
    // Zero out rows for boundary vertices in A
    for (int col = 0; col < A.outerSize(); ++col) {
        for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(A, col); it; ++it) {
            int row = it.row();
            if (boundarySet.count(row)) {
                it.valueRef() = (row == col) ? one : zero;
            }
        }
    }

    // Prune zeros from the matrices
    M.prune(zero);
    A.prune(zero);

    // Factorize the system matrix for efficient solves
    // Try SparseLU first (direct solver), fall back to BiCGSTAB (iterative) if it fails
    solverLU.compute(A);
    if (solverLU.info() != Eigen::Success) {
        std::cerr << "SparseLU factorization failed, falling back to BiCGSTAB iterative solver" << std::endl;
        useBiCGSTAB = true;
        solverBiCGSTAB.setTolerance(1e-10);
        solverBiCGSTAB.setMaxIterations(1000);
        solverBiCGSTAB.compute(A);
        if (solverBiCGSTAB.info() != Eigen::Success) {
            throw std::runtime_error("Both SparseLU and BiCGSTAB factorization failed");
        }
    } else {
        useBiCGSTAB = false;
    }

    // Initialize non-boundary vertices based on method
    if (method == 0) {
        // TODO: implement method 0 initialization
    } else if (method == 1) {
        // Set non-boundary vertices to random normalized complex numbers
        std::srand(seed ? seed : static_cast<unsigned>(std::time(nullptr)));
        for (int v = 0; v < numVertices; ++v) {
            if (boundarySet.find(v) == boundarySet.end()) {
                // Generate random angle and create unit complex number
                double theta = 2.0 * M_PI * (static_cast<double>(std::rand()) / RAND_MAX);
                u_k_prev[v] = std::complex<double>(std::cos(theta), std::sin(theta));
            }
        }
    }

}

void CrossField::step() {
    // Solve (M + tau*K) u^{k+1} = M u^k
    // Compute RHS = M * u_k_prev
    Eigen::VectorXcd rhs = M * u_k_prev;
    
    // Solve A * u_k = rhs using the appropriate solver
    if (useBiCGSTAB) {
        u_k = solverBiCGSTAB.solve(rhs);
        if (solverBiCGSTAB.info() != Eigen::Success) {
            std::cerr << "BiCGSTAB solve failed at iteration " << solverBiCGSTAB.iterations() << std::endl;
            throw std::runtime_error("BiCGSTAB solve failed");
        }
    } else {
        u_k = solverLU.solve(rhs);
        if (solverLU.info() != Eigen::Success) {
            throw std::runtime_error("SparseLU solve failed");
        }
    }

    // reproject onto unit circle
    for (int i = 0; i < u_k.size(); ++i) {
        if (std::abs(u_k[i]) > 1e-14) {
            u_k[i] /= std::abs(u_k[i]);
        }
    }

    // update error for convergence checking
    error = (u_k - u_k_prev).norm();

    // Update for next iteration
    u_k_prev = u_k;
}

void CrossField::computeSingularities() {
    // Clear previous results
    singularTriangles.clear();

    // Helper lambda to compute the smallest angle difference between two complex numbers
    // Returns value in (-pi, pi]
    auto angleDiff = [](std::complex<double> z_start, std::complex<double> z_end) -> double {
        // diff_complex = z_end * conj(z_start), phase gives the angle difference
        std::complex<double> diff_complex = z_end * std::conj(z_start);
        return std::atan2(diff_complex.imag(), diff_complex.real());
    };

    int numTriangles = static_cast<int>(mesh->triangles.size());

    // Iterate over every triangle in the mesh
    for (int t = 0; t < numTriangles; ++t) {
        const Triangle &tri = mesh->triangles[t];

        // Get the indices of the three vertices (assumed CCW ordered)
        int i = tri[0];
        int j = tri[1];
        int k = tri[2];

        // Check if any of the vertices are boundary vertices; if so, skip this triangle
        if (mesh->isBoundaryVertex[i] || mesh->isBoundaryVertex[j] || mesh->isBoundaryVertex[k]) {
            continue;
        }

        // Retrieve the complex values of the field at these vertices
        std::complex<double> val_i = u_k_prev[i];
        std::complex<double> val_j = u_k_prev[j];
        std::complex<double> val_k = u_k_prev[k];

        // Compute the change in angle along each edge of the triangle
        double d_theta_1 = angleDiff(val_i, val_j); // Edge i -> j
        double d_theta_2 = angleDiff(val_j, val_k); // Edge j -> k
        double d_theta_3 = angleDiff(val_k, val_i); // Edge k -> i

        // Sum the angle changes around the loop
        double total_angle_change = d_theta_1 + d_theta_2 + d_theta_3;

        // The winding number is the total change divided by 2*pi
        // It should be an integer (or very close to one due to float error)
        int winding_number = static_cast<int>(std::round(total_angle_change / (2.0 * M_PI)));

        // Check if a singularity exists (non-zero winding number)
        if (winding_number != 0) {
            // The cross field index is 1/N * winding number of representation field
            // Here N=4
            double cross_index = winding_number / 4.0;

            // Record the singularity (triangle index, cross-field index)
            singularTriangles.emplace_back(t, cross_index);
        }
    }
}

void CrossField::runMBO() {
    int iteration = 0;
    double nverts = static_cast<double>(mesh->vertices.size());
    while (iteration < maxIterations && error > 2 * nverts * 1e-5) {
        step();
        iteration++;
    }
}