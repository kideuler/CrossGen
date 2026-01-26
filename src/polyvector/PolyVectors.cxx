#include "PolyVectors.hxx"
#include <fstream>
#include <iomanip>
#include <unordered_set>

static constexpr double kPi = 3.141592653589793238462643383279502884;

static inline double wrapToPi(double a) {
    // Map angle to (-pi, pi]
    constexpr double TWO_PI = 2.0 * kPi;
    a = std::fmod(a + kPi, TWO_PI);
    if (a < 0.0) a += TWO_PI;
    a -= kPi;
    // Move -pi to +pi for consistency
    if (a <= -kPi + 1e-15) a += TWO_PI;
    return a;
}

static inline double betaFromDir4(const Point &d) {
    const double a = std::atan2(d[1], d[0]);
    return wrapToPi(4.0 * a);
}

static inline Point subPoint(const Point &a, const Point &b) {
    return Point{a[0] - b[0], a[1] - b[1]};
}

static inline Point triCentroid(const Mesh &m, const Triangle &t) {
    const Point &a = m.vertices[t[0]];
    const Point &b = m.vertices[t[1]];
    const Point &c = m.vertices[t[2]];
    return Point{ (a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0 };
}

static inline int otherSharedVertexBesides(const Triangle &ta, const Triangle &tb, int v) {
    // Returns w if triangles share edge (v,w), else -1.
    int w = -1;
    for (int i = 0; i < 3; ++i) {
        int a = ta[i];
        if (a == v) continue;
        if (tb[0] == a || tb[1] == a || tb[2] == a) {
            w = a;
            break;
        }
    }
    return w;
}

static inline bool angleBetweenCCW(double a, double b, double x) {
    // All angles in [-pi,pi] or any; works by mapping to [0,2pi).
    constexpr double TWO_PI = 2.0 * kPi;
    auto norm = [](double t) {
        double r = std::fmod(t, TWO_PI);
        if (r < 0.0) r += TWO_PI;
        return r;
    };
    a = norm(a);
    b = norm(b);
    x = norm(x);
    const double ab = std::fmod(b - a + TWO_PI, TWO_PI);
    const double ax = std::fmod(x - a + TWO_PI, TWO_PI);
    return ax <= ab;
}

static inline std::vector<int> boundaryNeighborsInTriangle(const Mesh &m, int triIdx, int v) {
    std::vector<int> out;
    const Triangle &tri = m.triangles[triIdx];
    for (int e = 0; e < 3; ++e) {
        if (m.triangleAdjacency[triIdx][e] != -1) continue;
        int a = tri[e];
        int b = tri[(e + 1) % 3];
        if (a == v) out.push_back(b);
        else if (b == v) out.push_back(a);
    }
    // Dedup (can happen in degenerate meshes)
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

PolyField::PolyField(Mesh &mesh) : mesh(mesh) {
    // Initialize polynomial coefficients
    int nTriangles = static_cast<int>(mesh.triangles.size());
    int nBdyTriangles = static_cast<int>(mesh.boundaryTriangles.size());

    // reserve space in sparse matrix L
    // Assuming L will be an nTriangles x nTriangles operator (e.g., Laplacian-like)
    L.resize(nTriangles, nTriangles);
    L.reserve(3 * nTriangles);

    // fill in L
    for (int t = 0; t < nTriangles; ++t) {
        // initialize all polynomial coefficients to zero
        polyCoeffs.push_back(PolyCoeffs{ std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0) });

        // initialize matrix L
        L.insert(t, t) = 0.0;
        for (int e = 0; e < 3; ++e) {
            int neighbor = mesh.triangleAdjacency[t][e];
            if (neighbor != -1) {
                // interior edge: add connection to neighbor
                L.insert(t, neighbor) = -1.0;
                L.coeffRef(t, t) += 1.0;
            }
        }
    }

    // Loop through boundary triangles to set polynomial coefficients and dirichlet conditions
    for (int i = 0; i < nBdyTriangles; ++i) {
        int t = mesh.boundaryTriangles[i][0];
        int e = mesh.boundaryTriangles[i][1];

        // Get the two vertices of the boundary edge
        const Triangle &tri = mesh.triangles[t];
        int v0 = tri[e];
        int v1 = tri[(e + 1) % 3];  // next vertex in CCW order 
        Point p0 = mesh.vertices[v0];
        Point p1 = mesh.vertices[v1];

        // Compute tangent vector along the edge
        Point tangent = { p1[0] - p0[0], p1[1] - p0[1] };
        // Compute polynomial coefficients from tangent vector
        polyCoeffs[t] = computePolyCoeffsFromTangentVector(tangent);

        // Set Dirichlet condition in L (make row t a unit row)
        // Zero out all existing nonzeros in row t
        for (int col = 0; col < L.outerSize(); ++col) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(L, col); it; ++it) {
                if (it.row() == t) {
                    it.valueRef() = 0.0;
                }
            }
        }
        // Remove zeros from sparsity pattern
        L.prune(0.0);
        // Set the diagonal to 1
        L.coeffRef(t, t) = 1.0;
        
    }

    // Also handle corner triangles (triangles with 2 boundary edges)
    // Use the first boundary edge's tangent as the constraint
    for (const auto &ct : mesh.cornerTriangles) {
        int t = ct[0];
        int e = ct[1]; // first boundary edge

        const Triangle &tri = mesh.triangles[t];
        int v0 = tri[e];
        int v1 = tri[(e + 1) % 3];
        Point p0 = mesh.vertices[v0];
        Point p1 = mesh.vertices[v1];

        Point tangent = { p1[0] - p0[0], p1[1] - p0[1] };
        polyCoeffs[t] = computePolyCoeffsFromTangentVector(tangent);

        // Set Dirichlet condition
        for (int col = 0; col < L.outerSize(); ++col) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(L, col); it; ++it) {
                if (it.row() == t) {
                    it.valueRef() = 0.0;
                }
            }
        }
        L.prune(0.0);
        L.coeffRef(t, t) = 1.0;
    }

    // factorize L after all modifications
    solver.compute(L);
    if(solver.info() != Eigen::Success) {
        return;
    }
}

void PolyField::solveForPolyCoeffs() {
    // Prepare rhs vector b
    int nTriangles = static_cast<int>(mesh.triangles.size());
    b_re = Eigen::VectorXd::Zero(nTriangles);
    b_im = Eigen::VectorXd::Zero(nTriangles);

    for (int m = 0; m < 2; ++m) {

        // Set rhs for boundary triangles based on polynomial coefficients
        for (int i = 0; i < static_cast<int>(mesh.boundaryTriangles.size()); ++i) {
            int t = mesh.boundaryTriangles[i][0]; // triangle index
            b_re(t) = polyCoeffs[t][m].real();
            b_im(t) = polyCoeffs[t][m].imag();
        }

        // Also set rhs for corner triangles
        for (const auto &ct : mesh.cornerTriangles) {
            int t = ct[0]; // triangle index
            b_re(t) = polyCoeffs[t][m].real();
            b_im(t) = polyCoeffs[t][m].imag();
        }

        // Solve the linear system L * x = b
        x_re = solver.solve(b_re);
        x_im = solver.solve(b_im);

        // Update polynomial coefficients based on solution x
        for (int t = 0; t < nTriangles; ++t) {
            // Here we simply set x_0 to the solved value and keep x_1 unchanged
            polyCoeffs[t][m] = std::complex<double>(x_re(t), x_im(t));
        }
    }
}

void PolyField::convertToFieldVectors() {
    int nTriangles = static_cast<int>(mesh.triangles.size());
    field.resize(nTriangles);

    for (int t = 0; t < nTriangles; ++t) {
        // Extract polynomial coefficients
        PolyCoeffs coeffs = polyCoeffs[t];

       Eigen::MatrixXcd M(4,4);
       M << std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), -coeffs[0],
           std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
           std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0), -coeffs[1],
           std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0);

       // Compute eigenvalues (no eigenvectors needed)
       Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(M, /* computeEigenvectors = */ false);

       // Put the 4 complex roots into a std::array (double precision)
       Eigen::Matrix<std::complex<double>, 4, 1> ev = es.eigenvalues();
        std::array<std::complex<double>, 4> roots;
        for (int k = 0; k < 4; ++k) {
            roots[k] = ev(k);
        }

        // Compute field vectors from roots
       auto dirs = extract_two_directions(roots);
       field[t].u = dirs.first;
       field[t].v = dirs.second;
    }

    // Keep `uSingularities` in sync with the newly computed field.
    computeUSingularities();

    // Also compute triangle rotations for the new field.
    computeTriangleRotations();
}

void PolyField::computeUSingularities() {
    uSingularities.clear();

    const int nV = static_cast<int>(mesh.vertices.size());
    const int nT = static_cast<int>(mesh.triangles.size());
    if (nV == 0 || nT == 0 || static_cast<int>(field.size()) != nT) return;

    // Build a lightweight boundary-vertex marker and boundary adjacency (from boundary edges).
    std::vector<std::vector<int>> bNbr(nV);
    bNbr.reserve(nV);
    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = mesh.triangles[t];
        for (int e = 0; e < 3; ++e) {
            if (mesh.triangleAdjacency[t][e] != -1) continue;
            int a = tri[e];
            int b = tri[(e + 1) % 3];
            if (a >= 0 && a < nV && b >= 0 && b < nV) {
                bNbr[a].push_back(b);
                bNbr[b].push_back(a);
            }
        }
    }
    for (int v = 0; v < nV; ++v) {
        auto &nbrs = bNbr[v];
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    }

    // Traverse vertices using CSR one-rings.
    for (int v = 0; v < nV; ++v) {
        const int deg = mesh.vertexTriangles.vertexDegree(v);
        if (deg <= 0) continue;

        auto [tBegin, tEnd] = mesh.vertexTriangles.trianglesForVertex(v);
        std::vector<int> tris(tBegin, tEnd);
        if (tris.empty()) continue;

        const bool isBoundary = !bNbr[v].empty();

        // Precompute beta (=4*angle) for each incident triangle.
        std::vector<double> betaTri(tris.size(), 0.0);
        for (size_t i = 0; i < tris.size(); ++i) {
            int t = tris[i];
            if (t < 0 || t >= nT) continue;
            betaTri[i] = betaFromDir4(field[t].u);
        }

        // Interior vertex: close the CCW cycle.
        if (!isBoundary) {
            if (tris.size() < 2) continue;
            double sum = 0.0;
            for (size_t i = 0; i < tris.size(); ++i) {
                const double a = betaTri[i];
                const double b = betaTri[(i + 1) % tris.size()];
                sum += wrapToPi(b - a);
            }
            const double k = sum / (2.0 * kPi);
            const int index4 = static_cast<int>(std::lround(k));
            if (index4 != 0) {
                uSingularities.emplace_back(v, index4);
            }
            continue;
        }

        // Boundary vertex: open fan + close with the two boundary edges (tangents).
        // Identify the "gap" in the one-ring (pair of consecutive incident triangles that do not share an edge).
        int gap = -1;
        std::vector<double> triAngles(tris.size(), 0.0);
        {
            const Point &pv = mesh.vertices[v];
            for (size_t i = 0; i < tris.size(); ++i) {
                const Point c = triCentroid(mesh, mesh.triangles[tris[i]]);
                triAngles[i] = std::atan2(c[1] - pv[1], c[0] - pv[0]);
            }
        }

        // First attempt: use edge-adjacency test.
        if (tris.size() >= 2) {
            int bestGap = -1;
            double bestGapAngle = -1.0;
            for (size_t i = 0; i < tris.size(); ++i) {
                const int ta = tris[i];
                const int tb = tris[(i + 1) % tris.size()];
                const int shared = otherSharedVertexBesides(mesh.triangles[ta], mesh.triangles[tb], v);
                const bool shareEdge = (shared != -1);
                // Angular gap between consecutive centroid directions (positive CCW).
                double da = triAngles[(i + 1) % tris.size()] - triAngles[i];
                if (da < 0.0) da += 2.0 * kPi;
                if (!shareEdge) {
                    if (da > bestGapAngle) {
                        bestGapAngle = da;
                        bestGap = static_cast<int>(i);
                    }
                }
            }
            gap = bestGap;
        }

        // Fallback: pick the largest angular gap (works even for deg==1 or if adjacency test fails).
        if (gap < 0) {
            int bestGap = 0;
            double bestGapAngle = -1.0;
            for (size_t i = 0; i < tris.size(); ++i) {
                double da = triAngles[(i + 1) % tris.size()] - triAngles[i];
                if (da < 0.0) da += 2.0 * kPi;
                if (da > bestGapAngle) {
                    bestGapAngle = da;
                    bestGap = static_cast<int>(i);
                }
            }
            gap = bestGap;
        }

        const int m = static_cast<int>(tris.size());
        const int start = (gap + 1) % m;
        std::vector<int> chainTris(m);
        std::vector<double> chainBeta(m);
        std::vector<double> chainAng(m);
        for (int i = 0; i < m; ++i) {
            const int idx = (start + i) % m;
            chainTris[i] = tris[idx];
            chainBeta[i] = betaTri[idx];
            chainAng[i] = triAngles[idx];
        }

        // Determine the two boundary edges (neighbors) to close the fan.
        int wStart = -1;
        int wEnd = -1;

        // Start side: boundary edge incident to v in the first triangle of the chain.
        {
            const auto bn = boundaryNeighborsInTriangle(mesh, chainTris.front(), v);
            if (!bn.empty()) wStart = bn[0];
        }
        // End side: boundary edge incident to v in the last triangle of the chain.
        {
            const auto bn = boundaryNeighborsInTriangle(mesh, chainTris.back(), v);
            if (!bn.empty()) wEnd = bn[0];
        }

        // If we only saw one boundary neighbor in triangles (e.g., deg==1 corner triangle),
        // use the boundary-neighbor list to recover both edges and order them around the fan.
        if ((wStart < 0 || wEnd < 0 || wStart == wEnd) && bNbr[v].size() >= 2) {
            int a = bNbr[v][0];
            int b = bNbr[v][1];

            const Point &pv = mesh.vertices[v];
            const Point &pa = mesh.vertices[a];
            const Point &pb = mesh.vertices[b];
            const double angA = std::atan2(pa[1] - pv[1], pa[0] - pv[0]);
            const double angB = std::atan2(pb[1] - pv[1], pb[0] - pv[0]);

            // If we have at least one triangle, its centroid direction should lie inside the fan.
            const double angC = chainAng.empty() ? angA : chainAng[m / 2];
            if (angleBetweenCCW(angA, angB, angC)) {
                wStart = a;
                wEnd = b;
            } else {
                wStart = b;
                wEnd = a;
            }
        }

        if (wStart < 0 || wEnd < 0) {
            // Can't robustly close the boundary fan; fall back to interior-style cycle.
            if (tris.size() < 2) continue;
            double sum = 0.0;
            for (size_t i = 0; i < tris.size(); ++i) {
                const double a = betaTri[i];
                const double b = betaTri[(i + 1) % tris.size()];
                sum += wrapToPi(b - a);
            }
            const double k = sum / (2.0 * kPi);
            const int index4 = static_cast<int>(std::lround(k));
            if (index4 != 0) uSingularities.emplace_back(v, index4);
            continue;
        }

        const Point &pv = mesh.vertices[v];
        const Point &pS = mesh.vertices[wStart];
        const Point &pE = mesh.vertices[wEnd];
        const double betaStart = betaFromDir4(subPoint(pS, pv));
        const double betaEnd = betaFromDir4(subPoint(pE, pv));

        double sum = 0.0;
        // boundary start -> first triangle
        sum += wrapToPi(chainBeta[0] - betaStart);
        // triangles along the fan
        for (int i = 0; i < m - 1; ++i) {
            sum += wrapToPi(chainBeta[i + 1] - chainBeta[i]);
        }
        // last triangle -> boundary end
        sum += wrapToPi(betaEnd - chainBeta[m - 1]);
        // boundary corner turn (end -> start)
        sum += wrapToPi(betaStart - betaEnd);

        const double k = sum / (2.0 * kPi);
        const int index4 = static_cast<int>(std::lround(k));
        if (index4 != 0) {
            uSingularities.emplace_back(v, index4);
        }
    }
}

void PolyField::computeTriangleRotations() {
    const int nT = static_cast<int>(mesh.triangles.size());
    fieldTriangleRotation.resize(nT, {-1,-1,-1});

    for (int t = 0; t < nT; ++t) {
        const Triangle &tri = mesh.triangles[t];
        for (int e = 0; e < 3; ++e) {
            if (mesh.triangleAdjacency[t][e] == -1) {
                // Boundary edge: no rotation
                continue;
            }
            int neighbor = mesh.triangleAdjacency[t][e];

            double theta_curr = computeAngle(field[t].u);
            double theta_neigh = computeAngle(field[neighbor].u);
            int k = find_rotation_matrix( theta_curr, theta_neigh);
            fieldTriangleRotation[t][e] = k;
        }
    }
}

bool PolyField::writeVTK(const std::string &filename) const {
    // Legacy VTK UnstructuredGrid format, ASCII
    std::ofstream out(filename);
    if (!out) return false;

    const auto &V = mesh.vertices;
    const auto &T = mesh.triangles;
    int nPts = static_cast<int>(V.size());
    int nTri = static_cast<int>(T.size());

    out << "# vtk DataFile Version 3.0\n";
    out << "PolyField export\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    // Points (promote 2D points to 3D with z=0)
    out << "POINTS " << nPts << " float\n";
    out.setf(std::ios::fixed); out << std::setprecision(7);
    for (const auto &p : V) {
        out << static_cast<float>(p[0]) << " "
            << static_cast<float>(p[1]) << " 0\n";
    }

    // Cells: each triangle is written as: 3 i j k
    // Total size = sum over cells of (1 + numIndices) = nTri * (1 + 3)
    out << "\nCELLS " << nTri << " " << nTri * 4 << "\n";
    for (const auto &tri : T) {
        out << 3 << " " << tri[0] << " " << tri[1] << " " << tri[2] << "\n";
    }

    // Cell types: VTK_TRIANGLE == 5
    out << "\nCELL_TYPES " << nTri << "\n";
    for (int i = 0; i < nTri; ++i) out << 5 << "\n";

    // Cell data: per-triangle vectors
    out << "\nCELL_DATA " << nTri << "\n";
    out << "VECTORS field_u float\n";
    for (int i = 0; i < nTri; ++i) {
        float ux = 0.0f, uy = 0.0f;
        if (i < static_cast<int>(field.size())) {
            ux = static_cast<float>(field[i].u[0]);
            uy = static_cast<float>(field[i].u[1]);
        }
        out << ux << " " << uy << " 0\n";
    }
    out << "\nVECTORS field_v float\n";
    for (int i = 0; i < nTri; ++i) {
        float vx = 0.0f, vy = 0.0f;
        if (i < static_cast<int>(field.size())) {
            vx = static_cast<float>(field[i].v[0]);
            vy = static_cast<float>(field[i].v[1]);
        }
        out << vx << " " << vy << " 0\n";
    }

    return true;
}
