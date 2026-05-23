#include "UVGParam.hxx"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <queue>

#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>

UVGParam::UVGParam(const CutMesh& cutMesh)
    : cm_(cutMesh)
{
    const Mesh& mesh = cutMesh.getCutMesh();
    const int nV = static_cast<int>(mesh.vertices.size());
    const int nT = static_cast<int>(mesh.triangles.size());

    if (nV == 0 || nT == 0) {
        throw std::runtime_error("UVGParam: cut mesh is empty.");
    }

    const std::vector<Point>& xField = cutMesh.getUField(); // X_T  (u direction)
    const std::vector<Point>& yField = cutMesh.getVField(); // Y_T  (v direction)

    if (static_cast<int>(xField.size()) != nT ||
        static_cast<int>(yField.size()) != nT) {
        throw std::runtime_error("UVGParam: field size does not match triangle count.");
    }

    // ------------------------------------------------------------------
    // Assemble  A  and  [b^X, b^Y]
    // ------------------------------------------------------------------
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<size_t>(2*nT) * 9);

    Eigen::VectorXd b = Eigen::VectorXd::Zero(2*nV);

    for (int t = 0; t < nT; ++t) {
        const Triangle& tri = mesh.triangles[t];
        const int i = tri[0], j = tri[1], k = tri[2];
        const int idx[3] = {i, j, k};

        const Point& p0 = mesh.vertices[i];
        const Point& p1 = mesh.vertices[j];
        const Point& p2 = mesh.vertices[k];

        // Signed area  (positive for CCW).
        double A2 = (p1[0]-p0[0])*(p2[1]-p0[1]) - (p2[0]-p0[0])*(p1[1]-p0[1]);
        double A  = 0.5 * A2;

        // Gradients of the three hat functions (constant per triangle).
        //   grad phi_0 = ( y1-y2, x2-x1 ) / (2A)
        //   grad phi_1 = ( y2-y0, x0-x2 ) / (2A)
        //   grad phi_2 = ( y0-y1, x1-x0 ) / (2A)
        double invA2 = 1.0 / A2;
        Point g[3];
        g[0] = { (p1[1]-p2[1])*invA2,  (p2[0]-p1[0])*invA2 };
        g[1] = { (p2[1]-p0[1])*invA2,  (p0[0]-p2[0])*invA2 };
        g[2] = { (p0[1]-p1[1])*invA2,  (p1[0]-p0[0])*invA2 };

        const Point& X = xField[t];
        const Point& Y = yField[t];

        // Accumulate Laplacian entries and load vector.
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                double w = A * (g[a][0]*g[b][0] + g[a][1]*g[b][1]);
                triplets.emplace_back(idx[a], idx[b], w);
                triplets.emplace_back(idx[a]+nV, idx[b]+nV, w);
            }
            double dotX = X[0]*g[a][0] + X[1]*g[a][1];
            double dotY = Y[0]*g[a][0] + Y[1]*g[a][1];
            b(idx[a]) += A * dotX;
            b(idx[a]+nV) += A * dotY;
        }
    }

    // ------------------------------------------------------------------
    // Build sparse matrix.  Dirichlet pinning (vertex 0) is handled as
    // explicit constraint rows in C below — do NOT corrupt A here.
    // ------------------------------------------------------------------
    Eigen::SparseMatrix<double> A(2*nV, 2*nV);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();


    // Start building constraint matrix C and constraint RHS d.
    // Rows of C: 2 Dirichlet pins + isoline constraints + seam jump constraints.
    std::vector<Eigen::Triplet<double>> triplets_C;
    triplets_C.reserve(static_cast<size_t>(cm_.getNaturalBoundaryEdges().size()) * 2 + 4);
    // Over-allocate d: 2 Dirichlet + natural boundary + up to 2 per cut vertex (seam).
    Eigen::VectorXd d = Eigen::VectorXd::Zero(
        2 + static_cast<int>(cm_.getNaturalBoundaryEdges().size()) + 2 * nV + 4);

    // Natural boundary edge keys store ORIGINAL mesh vertex indices.
    // Map to cut mesh vertex indices via origToCutVerts.
    // Boundary vertices not on any cut have exactly one copy in the cut mesh.
    const auto& origToCut = cm_.getOriginalToCutVertices();

    // loop through natural boundary edges and add isoline constraints per Section 4.
    int constraintIdx = 0;

    // Dirichlet: pin u_0 = 0 and v_0 = 0 via constraint rows (vertex 0 of cut mesh).
    triplets_C.emplace_back(0,  constraintIdx, 1.0);  d(constraintIdx) = 0.0;  constraintIdx++;
    triplets_C.emplace_back(nV, constraintIdx, 1.0);  d(constraintIdx) = 0.0;  constraintIdx++;

    const auto &naturalBoundaryEdgeTriangle = cm_.getNaturalBoundaryEdgeTriangle();
    for (const auto& ek : cm_.getNaturalBoundaryEdges()) {
        // ek.a / ek.b are original mesh vertex indices — map to cut mesh vertices.
        if (ek.a < 0 || ek.a >= static_cast<int>(origToCut.size())) continue;
        if (ek.b < 0 || ek.b >= static_cast<int>(origToCut.size())) continue;
        if (origToCut[ek.a].empty() || origToCut[ek.b].empty()) continue;

        // Triangle index for crossfield lookup.
        int t = naturalBoundaryEdgeTriangle.at(ek);
        const Triangle& tri_t = mesh.triangles[t];

        // Find the cut-mesh copy of an original vertex that belongs to triangle t.
        // This is essential at seam/boundary intersections where a vertex has >1 copy.
        auto findCutVert = [&](int origV) -> int {
            for (int cv : origToCut[origV]) {
                if (cv == tri_t[0] || cv == tri_t[1] || cv == tri_t[2]) return cv;
            }
            return origToCut[origV][0]; // fallback (should not happen for valid mesh)
        };

        const int i = findCutVert(ek.a);
        const int j = findCutVert(ek.b);

        // Edge direction in 2D using the correctly identified cut vertices.
        Point edgeVec = mesh.vertices[j] - mesh.vertices[i];

        Point ut = cm_.getUField()[t];
        Point vt = cm_.getVField()[t];

        // Classify the edge per Section 4:
        //   |<t, e^u>| > |<t, e^v>|  =>  edge tangent to u-field  =>  v = constant
        //   otherwise                 =>  edge tangent to v-field  =>  u = constant
        double dotU = std::abs(dotP(edgeVec, ut));
        double dotV = std::abs(dotP(edgeVec, vt));

        if (dotU > dotV) {
            // Edge aligned with u-field → v must be constant: v_j - v_i = 0
            // col(v_i) = nV + i,  col(v_j) = nV + j
            triplets_C.emplace_back(nV + i, constraintIdx, -1.0);
            triplets_C.emplace_back(nV + j, constraintIdx,  1.0);
        } else {
            // Edge aligned with v-field → u must be constant: u_j - u_i = 0
            // col(u_i) = i,  col(u_j) = j
            triplets_C.emplace_back(i, constraintIdx, -1.0);
            triplets_C.emplace_back(j, constraintIdx,  1.0);
        }
        d(constraintIdx) = 0.0;
        constraintIdx++;
    }

    // ------------------------------------------------------------------
    // Section 5: Seam jump constraints with unknown per-component translation.
    // The field is combed, so all seam transitions are identity (R_s = I).
    // For each connected seam component s with tau unknowns (tau_s^u, tau_s^v):
    //   u_{i_R} - u_{i_L} - tau_s^u = 0
    //   v_{i_R} - v_{i_L} - tau_s^v = 0
    // for every original mesh vertex i on the seam that has exactly 2 cut copies.
    // ------------------------------------------------------------------

    // Build adjacency list of original mesh vertices along cut edges.
    const auto& cutEdges = cm_.getCutEdges();
    std::unordered_map<int, std::vector<int>> cutAdj;
    cutAdj.reserve(cutEdges.size() * 2);
    for (const auto& ek : cutEdges) {
        cutAdj[ek.a].push_back(ek.b);
        cutAdj[ek.b].push_back(ek.a);
    }

    // BFS to find connected components of the cut graph.
    std::unordered_map<int, int> vertToSeamComp; // original vertex -> component index
    vertToSeamComp.reserve(cutAdj.size());
    int numSeamComponents = 0;
    for (const auto& [startV, nbrs] : cutAdj) {
        if (vertToSeamComp.count(startV)) continue;
        std::queue<int> q;
        q.push(startV);
        vertToSeamComp[startV] = numSeamComponents;
        while (!q.empty()) {
            int v = q.front(); q.pop();
            for (int nb : cutAdj[v]) {
                if (!vertToSeamComp.count(nb)) {
                    vertToSeamComp[nb] = numSeamComponents;
                    q.push(nb);
                }
            }
        }
        numSeamComponents++;
    }

    // Tau unknowns: tau_s^u at primal DOF index 2*nV + 2*s,
    //               tau_s^v at primal DOF index 2*nV + 2*s + 1.
    // Full primal DOF count:
    const int nDofs = 2 * nV + 2 * numSeamComponents;

    // For each original vertex on the cut with exactly 2 cut copies, add 2 constraint rows.
    for (const auto& [origV, comp] : vertToSeamComp) {
        if (origV < 0 || origV >= static_cast<int>(origToCut.size())) continue;
        const auto& copies = origToCut[origV];
        if (copies.size() != 2) continue;  // skip junction vertices (>2 copies)

        const int iL = copies[0];
        const int iR = copies[1];
        if (iL < 0 || iL >= nV || iR < 0 || iR >= nV) continue;

        const int tauUdof = 2 * nV + 2 * comp;      // primal DOF col for tau_s^u
        const int tauVdof = 2 * nV + 2 * comp + 1;  // primal DOF col for tau_s^v

        // u_{iR} - u_{iL} - tau_s^u = 0
        triplets_C.emplace_back(iR,      constraintIdx,  1.0);
        triplets_C.emplace_back(iL,      constraintIdx, -1.0);
        triplets_C.emplace_back(tauUdof, constraintIdx, -1.0);
        d(constraintIdx) = 0.0;
        constraintIdx++;

        // v_{iR} - v_{iL} - tau_s^v = 0
        triplets_C.emplace_back(nV + iR, constraintIdx,  1.0);
        triplets_C.emplace_back(nV + iL, constraintIdx, -1.0);
        triplets_C.emplace_back(tauVdof, constraintIdx, -1.0);
        d(constraintIdx) = 0.0;
        constraintIdx++;
    }

    std::cerr << "[UVGParam] Seam components: " << numSeamComponents
              << ", total constraints: " << constraintIdx << "\n";

    // Form global KKT matrix  [A    C^T]   size (nDofs + nConstraints)^2
    //                         [C    0  ]
    // where nDofs = 2*nV + 2*numSeamComponents.
    // A occupies the top-left 2*nV x 2*nV block (tau DOFs have zero energy rows/cols).
    // Triplets_C rows are primal DOF indices in [0, nDofs); cols are constraint indices.

    std::vector<Eigen::Triplet<double>> kktTrips;
    kktTrips.reserve(A.nonZeros() + 2 * static_cast<int>(triplets_C.size()));

    // Top-left: A  (lives in [0, 2*nV) x [0, 2*nV))
    for (int col = 0; col < A.outerSize(); ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, col); it; ++it) {
            kktTrips.emplace_back(static_cast<int>(it.row()), static_cast<int>(it.col()), it.value());
        }
    }

    // Top-right: C^T  and  Bottom-left: C
    // triplets_C: (row = primal DOF index, col = constraint index, value)
    for (const auto& tc : triplets_C) {
        int primalRow = tc.row();           // primal DOF in [0, nDofs)
        int constrCol = tc.col();           // constraint index in [0, constraintIdx)
        double val    = tc.value();
        // C  block: row = nDofs + constrCol,  col = primalRow
        kktTrips.emplace_back(nDofs + constrCol, primalRow, val);
        // C^T block: row = primalRow,          col = nDofs + constrCol
        kktTrips.emplace_back(primalRow, nDofs + constrCol, val);
    }

    const int kktSize = nDofs + constraintIdx;
    Eigen::SparseMatrix<double> globalMat(kktSize, kktSize);
    globalMat.setFromTriplets(kktTrips.begin(), kktTrips.end());
    globalMat.makeCompressed();

    // ------------------------------------------------------------------
    // Solve with SparseLU on the indefinite KKT (saddle-point) system.
    // SimplicialLDLT/LLT require positive definiteness, which KKT systems
    // do not have (the zero bottom-right block makes them indefinite).
    // ------------------------------------------------------------------
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(kktSize);
    rhs.segment(0, 2 * nV) = b;
    // rhs[2*nV .. nDofs) = 0 — tau unknowns have no objective contribution.
    rhs.segment(nDofs, constraintIdx) = d.head(constraintIdx);

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(globalMat);
    solver.factorize(globalMat);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("UVGParam: factorization of augmented KKT system failed.");
    }
    auto x = solver.solve(rhs);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("UVGParam: linear solve failed.");
    }

    u_ = x.segment(0, nV);
    v_ = x.segment(nV, nV);

    // Extract per-vertex translation vectors from the tau DOFs.
    // tau_s^u = x[2*nV + 2*s],  tau_s^v = x[2*nV + 2*s + 1]
    // Scatter to per-cut-vertex arrays (zero for vertices not on any seam).
    const auto& cutToOrig = cm_.getCutVertexToOriginal();
    tu_ = Eigen::VectorXd::Zero(nV);
    tv_ = Eigen::VectorXd::Zero(nV);
    for (int cv = 0; cv < nV; ++cv) {
        const int origV = (cv < static_cast<int>(cutToOrig.size())) ? cutToOrig[cv] : -1;
        if (origV < 0) continue;
        auto it = vertToSeamComp.find(origV);
        if (it == vertToSeamComp.end()) continue;
        const int s = it->second;
        tu_(cv) = x(2 * nV + 2 * s);
        tv_(cv) = x(2 * nV + 2 * s + 1);
    }

    std::cerr << "[UVGParam] Solved UV parameterization: "
              << nV << " vertices, " << nT << " triangles.\n";
}

// ---------------------------------------------------------------------------
int UVGParam::numFlippedTriangles() const
{
    const Mesh& mesh = cm_.getCutMesh();
    const int nV = static_cast<int>(u_.size());
    int count = 0;

    for (const Triangle& tri : mesh.triangles) {
        const int i = tri[0], j = tri[1], k = tri[2];
        if (i < 0 || i >= nV || j < 0 || j >= nV || k < 0 || k >= nV) continue;

        // Signed area = 0.5 * ((uj-ui)*(vk-vi) - (uk-ui)*(vj-vi))
        // Negative means the UV triangle has opposite winding (flipped).
        double signedArea2 = (u_(j) - u_(i)) * (v_(k) - v_(i))
                           - (u_(k) - u_(i)) * (v_(j) - v_(i));
        if (signedArea2 < 0.0) ++count;
    }
    return count;
}

// ---------------------------------------------------------------------------
bool UVGParam::writeOBJ(const std::string& filename) const
{
    const Mesh& mesh = cm_.getCutMesh();
    const int nV = static_cast<int>(mesh.vertices.size());
    const int nT = static_cast<int>(mesh.triangles.size());

    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "[UVGParam] Cannot open file: " << filename << "\n";
        return false;
    }

    ofs << "# UVGParam output\n";

    // Geometric vertices (z=0 from the cut mesh 2D positions)
    for (int v = 0; v < nV; ++v) {
        const Point& p = mesh.vertices[v];
        ofs << "v " << p[0] << " " << p[1] << " 0\n";
    }

    // UV texture coordinates
    for (int v = 0; v < nV; ++v) {
        ofs << "vt " << u_(v) << " " << v_(v) << "\n";
    }

    // Faces with UV indices (1-based in OBJ)
    for (int t = 0; t < nT; ++t) {
        const Triangle& tri = mesh.triangles[t];
        ofs << "f "
            << (tri[0]+1) << "/" << (tri[0]+1) << " "
            << (tri[1]+1) << "/" << (tri[1]+1) << " "
            << (tri[2]+1) << "/" << (tri[2]+1) << "\n";
    }

    std::cerr << "[UVGParam] Written UV OBJ: " << filename << "\n";
    return true;
}
