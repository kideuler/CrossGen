#include "MIQ.hxx"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <tuple>

// ============================================================================
// Construction
// ============================================================================

MIQSolver::MIQSolver(const CutMesh& cutMesh)
    : cutMesh_(cutMesh)
{
    convertMeshToEigen();
    computeTriangleTriangleAdjacency();
    computeMismatch();
    computeSeams();
    initVertexIndexing();
}

// ============================================================================
// Convert mesh data to Eigen format
// ============================================================================

void MIQSolver::convertMeshToEigen() {
    const Mesh& orig = cutMesh_.getOriginalMesh();
    const Mesh& cut  = cutMesh_.getCutMesh();

    // Original mesh vertices (2D, but we store as Nx2)
    const int nVorig = static_cast<int>(orig.vertices.size());
    V_.resize(nVorig, 2);
    for (int i = 0; i < nVorig; ++i) {
        V_(i, 0) = orig.vertices[i][0];
        V_(i, 1) = orig.vertices[i][1];
    }

    // Original mesh faces
    const int nF = static_cast<int>(orig.triangles.size());
    F_.resize(nF, 3);
    for (int i = 0; i < nF; ++i) {
        F_(i, 0) = orig.triangles[i][0];
        F_(i, 1) = orig.triangles[i][1];
        F_(i, 2) = orig.triangles[i][2];
    }

    // Cut mesh vertices
    const int nVcut = static_cast<int>(cut.vertices.size());
    Vcut_.resize(nVcut, 2);
    for (int i = 0; i < nVcut; ++i) {
        Vcut_(i, 0) = cut.vertices[i][0];
        Vcut_(i, 1) = cut.vertices[i][1];
    }

    // Cut mesh faces (use same topology as original, but indices refer to cut vertices)
    // The cut mesh has the same number of faces as the original.
    const int nFcut = static_cast<int>(cut.triangles.size());
    Fcut_.resize(nFcut, 3);
    for (int i = 0; i < nFcut; ++i) {
        Fcut_(i, 0) = cut.triangles[i][0];
        Fcut_(i, 1) = cut.triangles[i][1];
        Fcut_(i, 2) = cut.triangles[i][2];
    }

    // Combed field directions (u and v)
    const auto& uField = cutMesh_.getUField();
    const auto& vField = cutMesh_.getVField();
    PD1_.resize(nF, 2);
    PD2_.resize(nF, 2);
    for (int i = 0; i < nF; ++i) {
        PD1_(i, 0) = uField[i][0];
        PD1_(i, 1) = uField[i][1];
        PD2_(i, 0) = vField[i][0];
        PD2_(i, 1) = vField[i][1];
    }

    // Singularity flags
    singular_.setZero(nVorig);
    const auto& sings = cutMesh_.getSingularities();
    for (const auto& s : sings) {
        if (s.first >= 0 && s.first < nVorig) {
            singular_(s.first) = 1;
        }
    }

    // Build vertex-triangle adjacency (on original mesh)
    VF_.clear();
    VFi_.clear();
    VF_.resize(nVorig);
    VFi_.resize(nVorig);
    for (int f = 0; f < nF; ++f) {
        for (int k = 0; k < 3; ++k) {
            int v = F_(f, k);
            VF_[v].push_back(f);
            VFi_[v].push_back(k);
        }
    }
}

// ============================================================================
// Triangle-triangle adjacency
// ============================================================================

void MIQSolver::computeTriangleTriangleAdjacency() {
    // We use the cut mesh's triangleAdjacency which was computed in buildExplicitCutMesh.
    // But we need to store the local edge index as well.
    const Mesh& cut = cutMesh_.getCutMesh();
    const int nF = static_cast<int>(cut.triangles.size());

    TT_.resize(nF, 3);
    TTi_.resize(nF, 3);
    TT_.setConstant(-1);
    TTi_.setConstant(-1);

    // Build edge map: EdgeKey -> (face, local edge)
    struct HalfEdge { int f; int e; };
    std::unordered_map<CutMesh::EdgeKey, HalfEdge, CutMesh::EdgeKeyHash> edgeMap;

    for (int f = 0; f < nF; ++f) {
        for (int e = 0; e < 3; ++e) {
            int v0 = Fcut_(f, e);
            int v1 = Fcut_(f, (e + 1) % 3);
            CutMesh::EdgeKey key(v0, v1);

            auto it = edgeMap.find(key);
            if (it == edgeMap.end()) {
                edgeMap[key] = {f, e};
            } else {
                // Found the neighboring half-edge
                int f2 = it->second.f;
                int e2 = it->second.e;
                TT_(f, e) = f2;
                TTi_(f, e) = e2;
                TT_(f2, e2) = f;
                TTi_(f2, e2) = e;
            }
        }
    }
}

// ============================================================================
// Compute mismatch across edges
// ============================================================================

void MIQSolver::computeMismatch() {
    // Mismatch quantifies the rotation (in multiples of pi/2) needed to align
    // the cross-field direction from face f0 to face f1 across edge e.
    // For a combed field (already aligned), mismatch should be 0 across non-seam edges.
    //
    // mismatch(f, e) = rotation index k in {0,1,2,3} such that:
    //   u_f1 ≈ R(k * pi/2) * u_f0

    const int nF = static_cast<int>(F_.rows());
    mismatch_.resize(nF, 3);
    mismatch_.setZero();

    auto angle = [](double x, double y) { return std::atan2(y, x); };
    auto wrapPi = [](double a) {
        a = std::fmod(a + M_PI, 2.0 * M_PI);
        if (a < 0) a += 2.0 * M_PI;
        return a - M_PI;
    };

    for (int f0 = 0; f0 < nF; ++f0) {
        double theta0 = angle(PD1_(f0, 0), PD1_(f0, 1));

        for (int e = 0; e < 3; ++e) {
            int f1 = TT_(f0, e);
            if (f1 == -1) continue;  // boundary edge

            double theta1 = angle(PD1_(f1, 0), PD1_(f1, 1));

            // Find k that minimizes |wrapPi(theta1 - theta0 - k*pi/2)|
            int bestK = 0;
            double bestDiff = std::abs(wrapPi(theta1 - theta0));
            for (int k = 1; k < 4; ++k) {
                double diff = std::abs(wrapPi(theta1 - theta0 - k * M_PI_2));
                if (diff < bestDiff) {
                    bestDiff = diff;
                    bestK = k;
                }
            }
            mismatch_(f0, e) = bestK;
        }
    }
}

// ============================================================================
// Compute seams (cut edges)
// ============================================================================

void MIQSolver::computeSeams() {
    // seams_(f, e) is non-zero if edge e of face f is a seam (cut edge).
    // A seam edge is one where the original mesh vertices are the same,
    // but cut mesh vertices are different (i.e., the edge was duplicated).
    //
    // We identify seams by checking if adjacent faces share the same cut-mesh
    // vertex indices. If not, it's a seam.

    const int nF = static_cast<int>(F_.rows());
    seams_.resize(nF, 3);
    seams_.setZero();

    const auto& cutEdges = cutMesh_.getCutEdges();
    const Mesh& orig = cutMesh_.getOriginalMesh();

    for (int f = 0; f < nF; ++f) {
        for (int e = 0; e < 3; ++e) {
            int f1 = TT_(f, e);
            if (f1 == -1) continue;  // boundary, not a seam edge

            // Get original vertex indices for this edge
            int origV0 = F_(f, e);
            int origV1 = F_(f, (e + 1) % 3);

            // Check if this edge is in the cut set
            CutMesh::EdgeKey key(origV0, origV1);
            if (cutEdges.find(key) != cutEdges.end()) {
                seams_(f, e) = 1;
            }
        }
    }
}

// ============================================================================
// Initialize vertex indexing and seam info
// ============================================================================

void MIQSolver::initVertexIndexing() {
    // systemInfo holds:
    //  - numVertVariables: number of cut-mesh vertices
    //  - numIntegerCuts: number of integer translation variables (one per seam)
    //  - edgeSeamInfo: per-seam-vertex constraints

    systemInfo_.numVertVariables = static_cast<int>(Vcut_.rows());
    systemInfo_.edgeSeamInfo.clear();

    const int nF = static_cast<int>(F_.rows());

    // Build seam structure similar to libigl's VertexIndexing
    // For each seam edge, we need to track:
    //  - v0, v0p: corresponding cut-mesh vertices on each side
    //  - mismatch: rotation across the seam
    //  - integerVar: which integer variable controls this seam

    // First, find all seam edges and group them into connected seams
    // Each connected seam gets one integer variable

    // Track which edges we've visited
    Eigen::MatrixXi visited = Eigen::MatrixXi::Zero(nF, 3);

    // For vertex-vertex adjacency on seams
    std::vector<std::list<std::tuple<int, int, int, int, int>>> VVSeam(V_.rows());
    // Each entry: (neighbor_v, f0, k0, f1, k1)

    for (int f0 = 0; f0 < nF; ++f0) {
        for (int k0 = 0; k0 < 3; ++k0) {
            if (seams_(f0, k0) == 0) continue;
            if (visited(f0, k0)) continue;

            int f1 = TT_(f0, k0);
            if (f1 == -1) continue;

            int k1 = TTi_(f0, k0);

            // Mark both directions as visited
            visited(f0, k0) = 1;
            visited(f1, k1) = 1;

            // Original vertex indices
            int v0 = F_(f0, k0);
            int v1 = F_(f0, (k0 + 1) % 3);

            VVSeam[v0].push_back({v1, f0, k0, f1, k1});
            VVSeam[v1].push_back({v0, f0, k0, f1, k1});
        }
    }

    // Find start vertices (vertices that start/end a seam branch)
    std::vector<int> startVertices;
    std::vector<bool> isStartVertex(V_.rows(), false);
    for (int v = 0; v < static_cast<int>(V_.rows()); ++v) {
        if (!VVSeam[v].empty() && (VVSeam[v].size() != 2 || singular_(v) != 0)) {
            startVertices.push_back(v);
            isStartVertex[v] = true;
        }
    }

    // Walk along seams from each start vertex
    struct VertexInfo {
        int v;
        int f0, k0;
        int f1, k1;
    };

    std::vector<std::vector<VertexInfo>> verticesPerSeam;

    for (int startV : startVertices) {
        auto& neighbors = VVSeam[startV];
        size_t neighborSize = neighbors.size();

        for (size_t j = 0; j < neighborSize; ++j) {
            std::vector<VertexInfo> thisSeam;

            VertexInfo startVertex;
            startVertex.v = startV;
            startVertex.f0 = startVertex.k0 = startVertex.f1 = startVertex.k1 = -1;

            VertexInfo currentVertex = startVertex;
            thisSeam.push_back(currentVertex);

            auto nextIt = neighbors.begin();
            std::advance(nextIt, j > 0 ? 0 : 0);  // get first remaining
            if (nextIt == neighbors.end()) break;

            auto [nextV, nf0, nk0, nf1, nk1] = *nextIt;
            neighbors.erase(nextIt);

            VertexInfo prevVertex = startVertex;

            while (true) {
                prevVertex = currentVertex;
                currentVertex.v = nextV;
                currentVertex.f0 = nf0;
                currentVertex.k0 = nk0;
                currentVertex.f1 = nf1;
                currentVertex.k1 = nk1;

                thisSeam.push_back(currentVertex);

                auto& currNeighbors = VVSeam[nextV];

                // Remove previous vertex from neighbors
                for (auto it = currNeighbors.begin(); it != currNeighbors.end(); ++it) {
                    if (std::get<0>(*it) == prevVertex.v) {
                        currNeighbors.erase(it);
                        break;
                    }
                }

                if (currNeighbors.size() == 1 && !isStartVertex[currentVertex.v]) {
                    auto it = currNeighbors.begin();
                    std::tie(nextV, nf0, nk0, nf1, nk1) = *it;
                    currNeighbors.erase(it);
                } else {
                    break;
                }
            }

            if (thisSeam.size() > 1) {
                verticesPerSeam.push_back(thisSeam);
            }
        }
    }

    // Build edgeSeamInfo
    int integerVar = 0;
    for (const auto& seam : verticesPerSeam) {
        if (seam.size() < 2) continue;

        // Determine which side of the seam to use
        int priorVertexIdx = -1;
        if (seam.size() > 2) {
            const auto& v1 = seam[1];
            const auto& v2 = seam[2];
            int cv1_f0_k1 = Fcut_(v1.f0, (v1.k0 + 1) % 3);
            int cv2_f0_k0 = Fcut_(v2.f0, v2.k0);
            int cv2_f1_k1 = Fcut_(v2.f1, v2.k1);

            if (cv1_f0_k1 == cv2_f0_k0 || cv1_f0_k1 == cv2_f1_k1) {
                priorVertexIdx = Fcut_(v1.f0, v1.k0);
            } else {
                priorVertexIdx = Fcut_(v1.f1, v1.k1);
            }
        } else {
            const auto& v1 = seam[1];
            priorVertexIdx = Fcut_(v1.f0, v1.k0);
        }

        // Process each vertex on the seam
        for (size_t i = 1; i < seam.size(); ++i) {
            const auto& vertex = seam[i];

            int f, k, ff;
            if (priorVertexIdx == Fcut_(vertex.f0, vertex.k0)) {
                f = vertex.f0;
                ff = vertex.f1;
                k = vertex.k0;
            } else {
                f = vertex.f1;
                ff = vertex.f0;
                k = vertex.k1;
            }

            // Get seam vertex info
            int vtx0 = Fcut_(f, k);
            int vtx1 = Fcut_(f, (k + 1) % 3);
            int kff = TTi_(f, k);
            int vtx0p = Fcut_(ff, (kff + 1) % 3);
            int vtx1p = Fcut_(ff, kff);
            int mm = mismatch_(f, k);

            SeamInfo si;
            si.v0 = vtx0;
            si.v0p = vtx0p;
            si.mismatch = mm;
            si.integerVar = integerVar;
            systemInfo_.edgeSeamInfo.push_back(si);

            if (i == seam.size() - 1) {
                SeamInfo si2;
                si2.v0 = vtx1;
                si2.v0p = vtx1p;
                si2.mismatch = mm;
                si2.integerVar = integerVar;
                systemInfo_.edgeSeamInfo.push_back(si2);
            }

            priorVertexIdx = vtx1;
        }

        integerVar++;
    }

    systemInfo_.numIntegerCuts = integerVar;
}

// ============================================================================
// Main solve routine
// ============================================================================

void MIQSolver::solve(
    double gradientSize,
    double stiffness,
    bool directRound,
    unsigned int iter,
    unsigned int localIter,
    bool doRound,
    bool singularityRound,
    bool boundaryFeatures)
{
    // Store parameters for use in solvePoisson and updateStiffening
    stiffnessWeight_ = stiffness;
    boundaryFeatures_ = boundaryFeatures;

    // Normalize gradient size by mesh diameter
    Eigen::Vector2d minCorner = Vcut_.colwise().minCoeff();
    Eigen::Vector2d maxCorner = Vcut_.colwise().maxCoeff();
    double diameter = (maxCorner - minCorner).norm();
    if (diameter > 1e-10) {
        gradientSize = gradientSize / diameter;
    }

    // Initialize stiffness vector to 1.0 per face (as in libigl miq.impl)
    stiffnessVector_ = Eigen::VectorXd::Constant(F_.rows(), 1.0);

    // Initialize output matrices
    const int nF = static_cast<int>(F_.rows());
    WUV = Eigen::MatrixXd::Zero(nF, 6);
    UV = Eigen::MatrixXd::Zero(Vcut_.rows(), 2);
    FUV = Fcut_;

    if (iter > 0) {
        // Iterative stiffening to reduce flips
        for (unsigned int i = 0; i < iter; ++i) {
            solvePoisson(gradientSize, 1.0, directRound, localIter, doRound, singularityRound);
            int nflips = numFlips();
            std::cout << "ITERATION " << i << " FLIPS " << nflips << std::endl;

            bool folded = updateStiffening(gradientSize);
            if (!folded) break;
        }
    } else {
        solvePoisson(gradientSize, 1.0, directRound, localIter, doRound, singularityRound);
    }

    int nflips = numFlips();
    std::cout << "**** END OPTIMIZING #FLIPS " << nflips << " ****" << std::endl;
}

// ============================================================================
// Poisson solver
// ============================================================================

void MIQSolver::solvePoisson(
    double gradientSize,
    double gridResolution,
    bool directRound,
    unsigned int localIter,
    bool doRound,
    bool singularityRound)
{
    integerRounding_ = doRound;
    idsToRound_.clear();
    hardConstraints_.clear();

    // Find sizes
    nVertVars_ = systemInfo_.numVertVariables;
    nIntegerVars_ = systemInfo_.numIntegerCuts;
    numCutConstraint_ = static_cast<int>(systemInfo_.edgeSeamInfo.size());

    // Count boundary edges (true boundary, not seam cuts)
    // These will get hard feature constraints if enabled
    numBoundaryConstraints_ = 0;
    if (boundaryFeatures_) {
        const int nF = static_cast<int>(Fcut_.rows());
        for (int f = 0; f < nF; ++f) {
            for (int e = 0; e < 3; ++e) {
                if (TT_(f, e) == -1) {  // True boundary edge
                    numBoundaryConstraints_++;
                }
            }
        }
        std::cout << "Found " << numBoundaryConstraints_ << " boundary edges for hard feature constraints" << std::endl;
    } else {
        std::cout << "Boundary feature constraints disabled" << std::endl;
    }

    // Find a vertex to fix (preferably a singularity)
    nFixedVars_ = 0;
    for (int v = 0; v < static_cast<int>(V_.rows()); ++v) {
        if (singular_(v)) {
            hardConstraints_.push_back(v);
            nFixedVars_ = 1;
            break;
        }
    }
    if (nFixedVars_ == 0 && V_.rows() > 0) {
        hardConstraints_.push_back(0);
        nFixedVars_ = 1;
    }

    // Total constraints: seam constraints (2 per seam vertex) + fixed vertex (2) + boundary features (1 per edge)
    numConstraintEquations_ = numCutConstraint_ * 2 + nFixedVars_ * 2 + numBoundaryConstraints_;
    numTotalVars_ = (nVertVars_ + nIntegerVars_) * 2;

    // Allocate matrices
    Lhs_.resize(nVertVars_ * 2, nVertVars_ * 2);
    Constraints_.resize(numConstraintEquations_, numTotalVars_);
    rhs_.resize(nVertVars_ * 2);
    rhs_.setZero();
    constraintsRhs_.resize(numConstraintEquations_);
    constraintsRhs_.setZero();

    // Build Laplacian matrix and RHS
    buildLaplacianMatrix(gradientSize);

    // Build seam constraints
    buildSeamConstraints();

    // Add boundary feature constraints (hard edges snap to grid lines)
    if (boundaryFeatures_) {
        addBoundaryFeatureConstraints();
    }

    // Fix blocked vertices
    fixBlockedVertex();

    // Add singularity rounding
    if (singularityRound) {
        addSingularityRound();
    }

    // Solve mixed-integer system
    mixedIntegerSolve(gridResolution, directRound, localIter);

    // Map coordinates back
    mapCoords();
}

// ============================================================================
// Build Laplacian matrix
// ============================================================================

void MIQSolver::buildLaplacianMatrix(double vfscale) {
    // Build the Poisson energy:
    //   E = sum_f ||grad(u) - PD1||^2 + ||grad(v) - PD2||^2
    //
    // This leads to a Laplacian system.

    const int nF = static_cast<int>(Fcut_.rows());
    const int nV = static_cast<int>(Vcut_.rows());

    // Compute gradient operator and cotangent weights
    // For a 2D mesh, the gradient is computed per-triangle.

    // Build cotangent Laplacian directly for the 2D case
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(nF * 9);

    // Also build RHS from field alignment
    Eigen::VectorXd rhsU = Eigen::VectorXd::Zero(nV);
    Eigen::VectorXd rhsV = Eigen::VectorXd::Zero(nV);

    for (int f = 0; f < nF; ++f) {
        int i0 = Fcut_(f, 0);
        int i1 = Fcut_(f, 1);
        int i2 = Fcut_(f, 2);

        Eigen::Vector2d p0(Vcut_(i0, 0), Vcut_(i0, 1));
        Eigen::Vector2d p1(Vcut_(i1, 0), Vcut_(i1, 1));
        Eigen::Vector2d p2(Vcut_(i2, 0), Vcut_(i2, 1));

        Eigen::Vector2d e0 = p2 - p1;
        Eigen::Vector2d e1 = p0 - p2;
        Eigen::Vector2d e2 = p1 - p0;

        double area2 = std::abs(e2[0] * e1[1] - e2[1] * e1[0]);
        if (area2 < 1e-12) continue;

        double area = 0.5 * area2;
        double stiff = stiffnessVector_(f);

        // Cotangent weights for each edge
        // cot(angle at vertex i) = (e_j . e_k) / |e_j x e_k|
        auto cotWeight = [&](const Eigen::Vector2d& ea, const Eigen::Vector2d& eb) {
            double cross = ea[0] * eb[1] - ea[1] * eb[0];
            double dot = ea.dot(eb);
            if (std::abs(cross) < 1e-12) return 0.0;
            return dot / std::abs(cross);
        };

        double c0 = cotWeight(-e1, e2);   // angle at vertex 0
        double c1 = cotWeight(-e2, e0);   // angle at vertex 1
        double c2 = cotWeight(-e0, e1);   // angle at vertex 2

        // Scale by stiffness
        c0 *= stiff;
        c1 *= stiff;
        c2 *= stiff;

        // Add to Laplacian (for both u and v components)
        // L[i,j] for u-component is at (2*i, 2*j)
        // L[i,j] for v-component is at (2*i+1, 2*j+1)

        // Edge (1,2) opposite to vertex 0
        triplets.push_back({2*i1, 2*i1, c0});
        triplets.push_back({2*i2, 2*i2, c0});
        triplets.push_back({2*i1, 2*i2, -c0});
        triplets.push_back({2*i2, 2*i1, -c0});

        triplets.push_back({2*i1+1, 2*i1+1, c0});
        triplets.push_back({2*i2+1, 2*i2+1, c0});
        triplets.push_back({2*i1+1, 2*i2+1, -c0});
        triplets.push_back({2*i2+1, 2*i1+1, -c0});

        // Edge (2,0) opposite to vertex 1
        triplets.push_back({2*i2, 2*i2, c1});
        triplets.push_back({2*i0, 2*i0, c1});
        triplets.push_back({2*i2, 2*i0, -c1});
        triplets.push_back({2*i0, 2*i2, -c1});

        triplets.push_back({2*i2+1, 2*i2+1, c1});
        triplets.push_back({2*i0+1, 2*i0+1, c1});
        triplets.push_back({2*i2+1, 2*i0+1, -c1});
        triplets.push_back({2*i0+1, 2*i2+1, -c1});

        // Edge (0,1) opposite to vertex 2
        triplets.push_back({2*i0, 2*i0, c2});
        triplets.push_back({2*i1, 2*i1, c2});
        triplets.push_back({2*i0, 2*i1, -c2});
        triplets.push_back({2*i1, 2*i0, -c2});

        triplets.push_back({2*i0+1, 2*i0+1, c2});
        triplets.push_back({2*i1+1, 2*i1+1, c2});
        triplets.push_back({2*i0+1, 2*i1+1, -c2});
        triplets.push_back({2*i1+1, 2*i0+1, -c2});

        // RHS from field alignment
        // grad(u) should align with PD1, grad(v) should align with -PD2
        // Using gradient formula: grad = (1/2A) * sum_i (n x e_i) * u_i
        // where n is the normal (0,0,1) in 2D

        Eigen::Vector2d pd1(PD1_(f, 0), PD1_(f, 1));
        Eigen::Vector2d pd2(PD2_(f, 0), PD2_(f, 1));

        // The gradient operator G applied to u gives grad(u)
        // We want grad(u) = PD1 * vfscale
        // So RHS contribution is G^T * area * stiff * PD1 * vfscale

        // Gradient basis vectors for this triangle
        Eigen::Vector2d g0 = Eigen::Vector2d(-e0[1], e0[0]) / area2;  // gradient of basis function 0
        Eigen::Vector2d g1 = Eigen::Vector2d(-e1[1], e1[0]) / area2;
        Eigen::Vector2d g2 = Eigen::Vector2d(-e2[1], e2[0]) / area2;

        double scale = area * stiff * vfscale * 0.5;

        // u-component RHS (align with PD1)
        rhsU(i0) += scale * g0.dot(pd1);
        rhsU(i1) += scale * g1.dot(pd1);
        rhsU(i2) += scale * g2.dot(pd1);

        // v-component RHS (align with -PD2 for proper orientation)
        rhsV(i0) -= scale * g0.dot(pd2);
        rhsV(i1) -= scale * g1.dot(pd2);
        rhsV(i2) -= scale * g2.dot(pd2);
    }

    Lhs_.setFromTriplets(triplets.begin(), triplets.end());

    // Interleave RHS: [u0, v0, u1, v1, ...]
    for (int i = 0; i < nV; ++i) {
        rhs_(2*i) = rhsU(i);
        rhs_(2*i+1) = rhsV(i);
    }
}

// ============================================================================
// Build seam constraints
// ============================================================================

void MIQSolver::buildSeamConstraints() {
    // For each seam vertex pair (v0, v0p), we have:
    //   R(mismatch) * uv0 - uv0p + t = 0
    // where t is an integer translation variable.

    int constrRow = 0;

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(numCutConstraint_ * 8);

    for (int i = 0; i < numCutConstraint_; ++i) {
        const SeamInfo& si = systemInfo_.edgeSeamInfo[i];

        // Adjust mismatch (libigl convention)
        int interval = si.mismatch;
        if (interval == 1) interval = 3;
        else if (interval == 3) interval = 1;

        int p0 = si.v0;
        int p0p = si.v0p;

        std::complex<double> rot = getRotationComplex(interval);

        int integerVar = nVertVars_ + si.integerVar;

        if (integerRounding_) {
            idsToRound_.push_back(integerVar * 2);
            idsToRound_.push_back(integerVar * 2 + 1);
        }

        // Constraint: rot * (u0 + i*v0) - (u0p + i*v0p) + (tu + i*tv) = 0
        // Real part:  rot.real*u0 - rot.imag*v0 - u0p + tu = 0
        // Imag part:  rot.imag*u0 + rot.real*v0 - v0p + tv = 0

        triplets.push_back({constrRow, 2*p0, rot.real()});
        triplets.push_back({constrRow, 2*p0+1, -rot.imag()});
        triplets.push_back({constrRow+1, 2*p0, rot.imag()});
        triplets.push_back({constrRow+1, 2*p0+1, rot.real()});

        triplets.push_back({constrRow, 2*p0p, -1.0});
        triplets.push_back({constrRow+1, 2*p0p+1, -1.0});

        triplets.push_back({constrRow, 2*integerVar, 1.0});
        triplets.push_back({constrRow+1, 2*integerVar+1, 1.0});

        constraintsRhs_(constrRow) = 0;
        constraintsRhs_(constrRow+1) = 0;

        constrRow += 2;
    }

    Constraints_.setFromTriplets(triplets.begin(), triplets.end());
}

// ============================================================================
// Add boundary feature constraints (hard edges)
// ============================================================================

void MIQSolver::addBoundaryFeatureConstraints() {
    // For each boundary edge (true boundary, not seam), we add a constraint
    // that both vertices of the edge have the same U or V coordinate,
    // effectively snapping the edge to an integer grid line.
    //
    // Which coordinate (U or V) depends on which direction the edge is
    // more aligned with in the cross-field:
    //   - If edge aligns with PD1 (u-direction), constrain V coordinates equal
    //   - If edge aligns with PD2 (v-direction), constrain U coordinates equal
    //
    // The constraint is: v1_coord - v2_coord = 0
    // Both vertices are also added to idsToRound_ to force integer rounding.

    const int nF = static_cast<int>(Fcut_.rows());
    int constrRow = numCutConstraint_ * 2;  // Start after seam constraints

    for (int f = 0; f < nF; ++f) {
        for (int e = 0; e < 3; ++e) {
            if (TT_(f, e) != -1) continue;  // Not a boundary edge

            // Get vertex indices in cut mesh
            int v0 = Fcut_(f, e);
            int v1 = Fcut_(f, (e + 1) % 3);

            // Get vertex positions
            Eigen::Vector2d p0(Vcut_(v0, 0), Vcut_(v0, 1));
            Eigen::Vector2d p1(Vcut_(v1, 0), Vcut_(v1, 1));

            // Edge direction
            Eigen::Vector2d edgeDir = (p1 - p0).normalized();

            // Get field directions for this face
            Eigen::Vector2d pd1(PD1_(f, 0), PD1_(f, 1));
            Eigen::Vector2d pd2(PD2_(f, 0), PD2_(f, 1));
            pd1.normalize();
            pd2.normalize();

            // Determine which direction the edge is more aligned with
            // using absolute dot product (direction doesn't matter)
            double alignWithU = std::abs(edgeDir.dot(pd1));  // alignment with u-direction
            double alignWithV = std::abs(edgeDir.dot(pd2));  // alignment with v-direction

            // offset = 0 means constrain U coordinates, offset = 1 means constrain V
            // If edge aligns with pd1 (u-dir), we want v1-v2 same V coord (offset=1)
            // If edge aligns with pd2 (v-dir), we want v1-v2 same U coord (offset=0)
            int offset = (alignWithU > alignWithV) ? 1 : 0;

            // Add both vertices to rounding list for this coordinate
            if (integerRounding_) {
                idsToRound_.push_back(v0 * 2 + offset);
                idsToRound_.push_back(v1 * 2 + offset);
            }

            // Add constraint: v0_coord - v1_coord = 0
            Constraints_.coeffRef(constrRow, v0 * 2 + offset) = 1.0;
            Constraints_.coeffRef(constrRow, v1 * 2 + offset) = -1.0;
            constraintsRhs_(constrRow) = 0.0;

            constrRow++;
        }
    }
}

// ============================================================================
// Fix blocked vertices
// ============================================================================

void MIQSolver::fixBlockedVertex() {
    int offsetRow = numCutConstraint_ * 2 + numBoundaryConstraints_;

    for (int i = 0; i < static_cast<int>(hardConstraints_.size()); ++i) {
        int v = hardConstraints_[i];

        // Get first cut-mesh vertex corresponding to this original vertex
        int index = getFirstVertexIndex(v);
        int indexVert = index * 2;

        int indexRow = offsetRow + i * 2;

        // Add fixing constraint
        Constraints_.coeffRef(indexRow, indexVert) += 1;
        Constraints_.coeffRef(indexRow+1, indexVert+1) += 1;

        // Fix to origin
        constraintsRhs_(indexRow) = 0;
        constraintsRhs_(indexRow+1) = 0;
    }
}

// ============================================================================
// Add singularity rounding
// ============================================================================

void MIQSolver::addSingularityRound() {
    for (int v = 0; v < static_cast<int>(V_.rows()); ++v) {
        if (singular_(v)) {
            int index0 = getFirstVertexIndex(v);
            idsToRound_.push_back(index0 * 2);
            idsToRound_.push_back(index0 * 2 + 1);
        }
    }
}

// ============================================================================
// Mixed-integer solve
// ============================================================================

void MIQSolver::mixedIntegerSolve(double coneGridRes, bool directRound, unsigned int localIter) {
    const double PENALIZATION = 0.000001;

    const int sizeMatrix = (nVertVars_ + nIntegerVars_) * 2;
    const int scalarSize = nVertVars_ * 2;

    // Allocate solution vector
    X_.resize(sizeMatrix, 0.0);

    // Build CoMISo matrices
    COMISO::ConstrainedSolver::ColMatrix A(sizeMatrix, sizeMatrix);
    COMISO::ConstrainedSolver::RowMatrix C(numConstraintEquations_, sizeMatrix + 1);
    COMISO::ConstrainedSolver::Vector B(sizeMatrix);
    B.setZero();

    // Copy LHS matrix
    for (int k = 0; k < Lhs_.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Lhs_, k); it; ++it) {
            A.coeffRef(it.row(), it.col()) += it.value();
        }
    }

    // Copy constraints matrix
    for (int k = 0; k < Constraints_.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Constraints_, k); it; ++it) {
            C.coeffRef(it.row(), it.col()) += it.value();
        }
    }

    // Add penalization for integer variables
    for (int i = 0; i < nIntegerVars_ * 2; ++i) {
        int index = scalarSize + i;
        A.coeffRef(index, index) = PENALIZATION;
    }

    // Copy RHS
    for (int i = 0; i < scalarSize; ++i) {
        B(i) = rhs_(i) * coneGridRes;
    }

    // Copy constraint RHS (goes in last column of C with negative sign)
    for (int i = 0; i < numConstraintEquations_; ++i) {
        C.coeffRef(i, sizeMatrix) = -constraintsRhs_(i) * coneGridRes;
    }

    // Sort and unique the rounding indices
    std::sort(idsToRound_.begin(), idsToRound_.end());
    auto newEnd = std::unique(idsToRound_.begin(), idsToRound_.end());
    idsToRound_.resize(std::distance(idsToRound_.begin(), newEnd));

    // Solve
    COMISO::ConstrainedSolver solver;
    solver.misolver().set_local_iters(localIter);

    if (directRound) {
        solver.misolver().set_direct_rounding();
    }

    COMISO::ConstrainedSolver::Vector _X(sizeMatrix);
    solver.solve(C, A, _X, B, idsToRound_, 0.0, false);

    // Copy solution back
    for (int i = 0; i < sizeMatrix; ++i) {
        X_[i] = _X(i);
    }
}

// ============================================================================
// Map coordinates
// ============================================================================

void MIQSolver::mapCoords() {
    const int nF = static_cast<int>(Fcut_.rows());

    // Per-wedge UV
    for (int f = 0; f < nF; ++f) {
        for (int k = 0; k < 3; ++k) {
            int indexUV = Fcut_(f, k);
            WUV(f, k*2 + 0) = X_[indexUV * 2];
            WUV(f, k*2 + 1) = X_[indexUV * 2 + 1];
        }
    }

    // Per-vertex UV
    for (int i = 0; i < static_cast<int>(Vcut_.rows()); ++i) {
        UV(i, 0) = X_[i * 2];
        UV(i, 1) = X_[i * 2 + 1];
    }
}

// ============================================================================
// Helper functions
// ============================================================================

std::complex<double> MIQSolver::getRotationComplex(int interval) const {
    assert(interval >= 0 && interval < 4);
    switch (interval) {
        case 0: return {1, 0};
        case 1: return {0, 1};
        case 2: return {-1, 0};
        default: return {0, -1};
    }
}

int MIQSolver::getFirstVertexIndex(int origV) const {
    // Return the first cut-mesh vertex corresponding to original vertex origV
    const auto& mapping = cutMesh_.getOriginalToCutVertices();
    if (origV >= 0 && origV < static_cast<int>(mapping.size()) && !mapping[origV].empty()) {
        return mapping[origV][0];
    }
    // Fallback: use VF adjacency
    if (origV >= 0 && origV < static_cast<int>(VF_.size()) && !VF_[origV].empty()) {
        int f = VF_[origV][0];
        int i = VFi_[origV][0];
        return Fcut_(f, i);
    }
    return 0;
}

// ============================================================================
// Flip detection and distortion
// ============================================================================

bool MIQSolver::isFlipped(const Eigen::Vector2d& uv0, const Eigen::Vector2d& uv1, const Eigen::Vector2d& uv2) {
    Eigen::Vector2d e0 = uv1 - uv0;
    Eigen::Vector2d e1 = uv2 - uv0;
    double area = e0(0) * e1(1) - e0(1) * e1(0);
    return area <= 0;
}

bool MIQSolver::isFlipped(int f) const {
    Eigen::Vector2d uv0(WUV(f, 0), WUV(f, 1));
    Eigen::Vector2d uv1(WUV(f, 2), WUV(f, 3));
    Eigen::Vector2d uv2(WUV(f, 4), WUV(f, 5));
    return isFlipped(uv0, uv1, uv2);
}

int MIQSolver::numFlips() const {
    int count = 0;
    for (int f = 0; f < static_cast<int>(F_.rows()); ++f) {
        if (isFlipped(f)) count++;
    }
    return count;
}

double MIQSolver::distortion(int f, double h) const {
    if (h <= 0) return 10.0;

    Eigen::Vector2d uv0(WUV(f, 0), WUV(f, 1));
    Eigen::Vector2d uv1(WUV(f, 2), WUV(f, 3));
    Eigen::Vector2d uv2(WUV(f, 4), WUV(f, 5));

    int i0 = Fcut_(f, 0);
    int i1 = Fcut_(f, 1);
    int i2 = Fcut_(f, 2);

    Eigen::Vector2d p0(Vcut_(i0, 0), Vcut_(i0, 1));
    Eigen::Vector2d p1(Vcut_(i1, 0), Vcut_(i1, 1));
    Eigen::Vector2d p2(Vcut_(i2, 0), Vcut_(i2, 1));

    // Compute area
    Eigen::Vector2d e1 = p1 - p0;
    Eigen::Vector2d e2 = p2 - p0;
    double area2 = std::abs(e1(0) * e2(1) - e1(1) * e2(0));

    if (area2 < 1e-12) return 10.0;

    double area2_inv = 1.0 / area2;

    // Edge vectors rotated 90 degrees (for gradient computation)
    Eigen::Vector2d neg_t0(-e1(1) + e2(1), e1(0) - e2(0));  // perpendicular to edge opposite vertex 0
    Eigen::Vector2d neg_t1(e2(1), -e2(0));                   // perpendicular to edge opposite vertex 1
    Eigen::Vector2d neg_t2(-e1(1), e1(0));                   // perpendicular to edge opposite vertex 2

    Eigen::Vector2d diffu = (neg_t0 * uv0(0) + neg_t1 * uv1(0) + neg_t2 * uv2(0)) * area2_inv;
    Eigen::Vector2d diffv = (neg_t0 * uv0(1) + neg_t1 * uv1(1) + neg_t2 * uv2(1)) * area2_inv;

    // First fundamental form
    double I00 = diffu.dot(diffu);
    double I01 = diffu.dot(diffv);
    double I11 = diffv.dot(diffv);

    // Eigenvalues
    double trI = I00 + I11;
    double diffDiag = I00 - I11;
    double sqrtDet = std::sqrt(std::max(0.0, diffDiag * diffDiag + 4 * I01 * I01));
    double sig1 = 0.5 * (trI + sqrtDet);
    double sig2 = 0.5 * (trI - sqrtDet);

    if (std::abs(sig2) < 1e-8) sig2 = 0;

    sig1 = std::sqrt(std::max(0.0, sig1));
    sig2 = std::sqrt(std::max(0.0, sig2));

    double tao = isFlipped(f) ? -1.0 : 1.0;
    double factor = tao / h;
    double lam = std::abs(factor * sig1 - 1) + std::abs(factor * sig2 - 1);

    return lam;
}

double MIQSolver::laplaceDistortion(int f, double h) const {
    double mydist = distortion(f, h);
    double lapl = 0;
    for (int i = 0; i < 3; ++i) {
        if (TT_(f, i) != -1) {
            lapl += mydist - distortion(TT_(f, i), h);
        }
    }
    return lapl;
}

bool MIQSolver::updateStiffening(double grad_size) {
    bool flipped = numFlips() > 0;
    if (!flipped) return false;

    double maxL = 0;
    double maxD = 0;

    // Parameters as in libigl miq.impl:
    // c = 1.0: multiplier for Laplacian distortion
    // d = stiffnessWeight_ (default 5.0): maximum stiffness delta per iteration
    const double c = 1.0;
    const double d = stiffnessWeight_;

    for (int i = 0; i < static_cast<int>(Fcut_.rows()); ++i) {
        double dist = distortion(i, grad_size);
        if (dist > maxD) maxD = dist;

        double absLap = std::abs(laplaceDistortion(i, grad_size));
        if (absLap > maxL) maxL = absLap;

        double stiffDelta = std::min(c * absLap, d);
        stiffnessVector_(i) += stiffDelta;
    }

    std::cout << "Maximum Distortion " << maxD << std::endl;
    std::cout << "Maximum Laplacian " << maxL << std::endl;

    return flipped;
}

// ============================================================================
// Output
// ============================================================================

bool MIQSolver::writeUVMesh(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out) return false;

    // Write UV vertices (as 3D with z=0)
    for (int i = 0; i < static_cast<int>(UV.rows()); ++i) {
        out << "v " << UV(i, 0) << " " << UV(i, 1) << " 0\n";
    }

    // Write faces
    for (int f = 0; f < static_cast<int>(FUV.rows()); ++f) {
        out << "f " << (FUV(f, 0) + 1) << " " << (FUV(f, 1) + 1) << " " << (FUV(f, 2) + 1) << "\n";
    }

    return true;
}

bool MIQSolver::writeTexturedMesh(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out) return false;

    // Get original mesh data from cutMesh
    const Mesh& origMesh = cutMesh_.getOriginalMesh();
    const auto& origVerts = origMesh.vertices;      // Original vertex positions
    const auto& origTris = origMesh.triangles;      // Original triangle connectivity
    
    // Write original mesh vertices (v lines)
    // These are the original 3D positions (though z=0 for 2D meshes)
    for (size_t i = 0; i < origVerts.size(); ++i) {
        out << "v " << origVerts[i][0] << " " << origVerts[i][1] << " 0\n";
    }

    // Write texture coordinates (vt lines)
    // We use per-wedge UV coordinates (WUV) which has one UV per face corner
    // WUV is nF x 6, storing [u0,v0, u1,v1, u2,v2] per triangle
    const int nF = static_cast<int>(WUV.rows());
    for (int f = 0; f < nF; ++f) {
        out << "vt " << WUV(f, 0) << " " << WUV(f, 1) << "\n";  // corner 0
        out << "vt " << WUV(f, 2) << " " << WUV(f, 3) << "\n";  // corner 1
        out << "vt " << WUV(f, 4) << " " << WUV(f, 5) << "\n";  // corner 2
    }

    // Write faces with vertex/texture indices (f v/vt lines)
    // Each face references original vertices and per-wedge texture coords
    for (int f = 0; f < nF; ++f) {
        int v0 = origTris[f][0] + 1;  // OBJ is 1-indexed
        int v1 = origTris[f][1] + 1;
        int v2 = origTris[f][2] + 1;
        
        // Texture coordinate indices (3 per face, sequential)
        int vt0 = f * 3 + 1;
        int vt1 = f * 3 + 2;
        int vt2 = f * 3 + 3;
        
        out << "f " << v0 << "/" << vt0 << " " 
                    << v1 << "/" << vt1 << " " 
                    << v2 << "/" << vt2 << "\n";
    }

    return true;
}

// ============================================================================
// Integer Isoline Tracing for Quad Mesh Extraction
// ============================================================================

// Helper: compute the position on an edge where a coordinate (u or v) equals a given integer value
// Returns the parametric t in [0,1] where the crossing occurs, or -1 if no crossing
static double edgeCrossingParam(double val0, double val1, int intVal) {
    double d = static_cast<double>(intVal);
    if ((val0 <= d && val1 >= d) || (val0 >= d && val1 <= d)) {
        double denom = val1 - val0;
        if (std::abs(denom) < 1e-12) return -1.0;
        double t = (d - val0) / denom;
        if (t >= 0.0 && t <= 1.0) return t;
    }
    return -1.0;
}

bool MIQSolver::extractQuadMesh() {
    // Integer isoline tracing algorithm:
    // 1. Find all integer isoline crossings on triangle edges
    // 2. Build quad corners at intersection points of U and V isolines
    // 3. Connect quad corners into quad faces
    
    if (UV.rows() == 0 || Fcut_.rows() == 0) {
        std::cerr << "extractQuadMesh: No UV data available\n";
        return false;
    }

    const int nF = static_cast<int>(Fcut_.rows());
    
    // Find the range of integer U and V values
    double uMin = UV.col(0).minCoeff();
    double uMax = UV.col(0).maxCoeff();
    double vMin = UV.col(1).minCoeff();
    double vMax = UV.col(1).maxCoeff();
    
    int uIntMin = static_cast<int>(std::floor(uMin));
    int uIntMax = static_cast<int>(std::ceil(uMax));
    int vIntMin = static_cast<int>(std::floor(vMin));
    int vIntMax = static_cast<int>(std::ceil(vMax));
    
    std::cout << "UV range: U=[" << uMin << "," << uMax << "], V=[" << vMin << "," << vMax << "]\n";
    std::cout << "Integer grid: U=[" << uIntMin << "," << uIntMax << "], V=[" << vIntMin << "," << vIntMax << "]\n";
    
    // Map from (intU, intV) -> quad vertex index
    // A quad corner exists at every integer (u,v) point that lies inside the domain
    std::map<std::pair<int,int>, int> gridToVertex;
    std::vector<Eigen::Vector2d> quadVerts;  // positions in original mesh coords
    
    // For each triangle, find which integer grid points lie inside it
    // Also find edge crossings for the isolines
    
    // Structure to hold isoline segment info
    struct IsoSegment {
        int tri;           // triangle index
        int isoVal;        // integer value of the isoline
        bool isU;          // true for U-isoline, false for V-isoline
        Eigen::Vector2d p0, p1;  // endpoints in original coordinates
        double u0, v0, u1, v1;   // UV coordinates of endpoints
    };
    
    std::vector<IsoSegment> uIsolines, vIsolines;
    
    // For each triangle, trace through it
    for (int f = 0; f < nF; ++f) {
        // Get UV coordinates
        Eigen::Vector2d uv0(UV(Fcut_(f, 0), 0), UV(Fcut_(f, 0), 1));
        Eigen::Vector2d uv1(UV(Fcut_(f, 1), 0), UV(Fcut_(f, 1), 1));
        Eigen::Vector2d uv2(UV(Fcut_(f, 2), 0), UV(Fcut_(f, 2), 1));
        
        // Get original coordinates
        Eigen::Vector2d p0(Vcut_(Fcut_(f, 0), 0), Vcut_(Fcut_(f, 0), 1));
        Eigen::Vector2d p1(Vcut_(Fcut_(f, 1), 0), Vcut_(Fcut_(f, 1), 1));
        Eigen::Vector2d p2(Vcut_(Fcut_(f, 2), 0), Vcut_(Fcut_(f, 2), 1));
        
        // Range of U values in this triangle
        double triUMin = std::min({uv0(0), uv1(0), uv2(0)});
        double triUMax = std::max({uv0(0), uv1(0), uv2(0)});
        double triVMin = std::min({uv0(1), uv1(1), uv2(1)});
        double triVMax = std::max({uv0(1), uv1(1), uv2(1)});
        
        // For each U-isoline that might cross this triangle
        for (int iu = static_cast<int>(std::floor(triUMin)); iu <= static_cast<int>(std::ceil(triUMax)); ++iu) {
            double uVal = static_cast<double>(iu);
            if (uVal < triUMin || uVal > triUMax) continue;
            
            // Find crossings on each edge
            std::vector<std::pair<double, int>> crossings;  // (t-param, edge index)
            
            // Edge 0: vertex 0 -> vertex 1
            double t01 = edgeCrossingParam(uv0(0), uv1(0), iu);
            if (t01 >= 0) crossings.push_back({t01, 0});
            
            // Edge 1: vertex 1 -> vertex 2
            double t12 = edgeCrossingParam(uv1(0), uv2(0), iu);
            if (t12 >= 0) crossings.push_back({t12, 1});
            
            // Edge 2: vertex 2 -> vertex 0
            double t20 = edgeCrossingParam(uv2(0), uv0(0), iu);
            if (t20 >= 0) crossings.push_back({t20, 2});
            
            if (crossings.size() >= 2) {
                // Compute the positions of the two crossings
                auto getPos = [&](double t, int e) -> Eigen::Vector2d {
                    if (e == 0) return p0 + t * (p1 - p0);
                    if (e == 1) return p1 + t * (p2 - p1);
                    return p2 + t * (p0 - p2);
                };
                auto getUV = [&](double t, int e) -> Eigen::Vector2d {
                    if (e == 0) return uv0 + t * (uv1 - uv0);
                    if (e == 1) return uv1 + t * (uv2 - uv1);
                    return uv2 + t * (uv0 - uv2);
                };
                
                IsoSegment seg;
                seg.tri = f;
                seg.isoVal = iu;
                seg.isU = true;
                seg.p0 = getPos(crossings[0].first, crossings[0].second);
                seg.p1 = getPos(crossings[1].first, crossings[1].second);
                Eigen::Vector2d uvA = getUV(crossings[0].first, crossings[0].second);
                Eigen::Vector2d uvB = getUV(crossings[1].first, crossings[1].second);
                seg.u0 = uvA(0); seg.v0 = uvA(1);
                seg.u1 = uvB(0); seg.v1 = uvB(1);
                uIsolines.push_back(seg);
            }
        }
        
        // For each V-isoline that might cross this triangle
        for (int iv = static_cast<int>(std::floor(triVMin)); iv <= static_cast<int>(std::ceil(triVMax)); ++iv) {
            double vVal = static_cast<double>(iv);
            if (vVal < triVMin || vVal > triVMax) continue;
            
            std::vector<std::pair<double, int>> crossings;
            
            double t01 = edgeCrossingParam(uv0(1), uv1(1), iv);
            if (t01 >= 0) crossings.push_back({t01, 0});
            
            double t12 = edgeCrossingParam(uv1(1), uv2(1), iv);
            if (t12 >= 0) crossings.push_back({t12, 1});
            
            double t20 = edgeCrossingParam(uv2(1), uv0(1), iv);
            if (t20 >= 0) crossings.push_back({t20, 2});
            
            if (crossings.size() >= 2) {
                auto getPos = [&](double t, int e) -> Eigen::Vector2d {
                    if (e == 0) return p0 + t * (p1 - p0);
                    if (e == 1) return p1 + t * (p2 - p1);
                    return p2 + t * (p0 - p2);
                };
                auto getUV = [&](double t, int e) -> Eigen::Vector2d {
                    if (e == 0) return uv0 + t * (uv1 - uv0);
                    if (e == 1) return uv1 + t * (uv2 - uv1);
                    return uv2 + t * (uv0 - uv2);
                };
                
                IsoSegment seg;
                seg.tri = f;
                seg.isoVal = iv;
                seg.isU = false;
                seg.p0 = getPos(crossings[0].first, crossings[0].second);
                seg.p1 = getPos(crossings[1].first, crossings[1].second);
                Eigen::Vector2d uvA = getUV(crossings[0].first, crossings[0].second);
                Eigen::Vector2d uvB = getUV(crossings[1].first, crossings[1].second);
                seg.u0 = uvA(0); seg.v0 = uvA(1);
                seg.u1 = uvB(0); seg.v1 = uvB(1);
                vIsolines.push_back(seg);
            }
        }
    }
    
    std::cout << "Found " << uIsolines.size() << " U-isoline segments, " 
              << vIsolines.size() << " V-isoline segments\n";
    
    // Now find intersection points between U and V isolines within each triangle
    // These are the quad vertices
    
    // Map from (triangle, intU, intV) to vertex index
    std::map<std::tuple<int, int, int>, int> triGridToVertex;
    
    for (int f = 0; f < nF; ++f) {
        // Get UV coordinates
        Eigen::Vector2d uv0(UV(Fcut_(f, 0), 0), UV(Fcut_(f, 0), 1));
        Eigen::Vector2d uv1(UV(Fcut_(f, 1), 0), UV(Fcut_(f, 1), 1));
        Eigen::Vector2d uv2(UV(Fcut_(f, 2), 0), UV(Fcut_(f, 2), 1));
        
        // Get original coordinates
        Eigen::Vector2d p0(Vcut_(Fcut_(f, 0), 0), Vcut_(Fcut_(f, 0), 1));
        Eigen::Vector2d p1(Vcut_(Fcut_(f, 1), 0), Vcut_(Fcut_(f, 1), 1));
        Eigen::Vector2d p2(Vcut_(Fcut_(f, 2), 0), Vcut_(Fcut_(f, 2), 1));
        
        // Range of integer values in this triangle
        double triUMin = std::min({uv0(0), uv1(0), uv2(0)});
        double triUMax = std::max({uv0(0), uv1(0), uv2(0)});
        double triVMin = std::min({uv0(1), uv1(1), uv2(1)});
        double triVMax = std::max({uv0(1), uv1(1), uv2(1)});
        
        // For each integer (u,v) point that might be in this triangle
        for (int iu = static_cast<int>(std::floor(triUMin)); iu <= static_cast<int>(std::ceil(triUMax)); ++iu) {
            for (int iv = static_cast<int>(std::floor(triVMin)); iv <= static_cast<int>(std::ceil(triVMax)); ++iv) {
                Eigen::Vector2d uvPt(static_cast<double>(iu), static_cast<double>(iv));
                
                // Check if this point is inside the triangle using barycentric coords
                // Compute barycentric coordinates
                Eigen::Vector2d v0 = uv1 - uv0;
                Eigen::Vector2d v1 = uv2 - uv0;
                Eigen::Vector2d v2 = uvPt - uv0;
                
                double d00 = v0.dot(v0);
                double d01 = v0.dot(v1);
                double d11 = v1.dot(v1);
                double d20 = v2.dot(v0);
                double d21 = v2.dot(v1);
                
                double denom = d00 * d11 - d01 * d01;
                if (std::abs(denom) < 1e-12) continue;
                
                double bary1 = (d11 * d20 - d01 * d21) / denom;
                double bary2 = (d00 * d21 - d01 * d20) / denom;
                double bary0 = 1.0 - bary1 - bary2;
                
                const double eps = -1e-6;
                if (bary0 >= eps && bary1 >= eps && bary2 >= eps) {
                    // Point is inside triangle - compute position in original coords
                    Eigen::Vector2d origPos = bary0 * p0 + bary1 * p1 + bary2 * p2;
                    
                    // Check if we already have this grid point
                    auto key = std::make_pair(iu, iv);
                    if (gridToVertex.find(key) == gridToVertex.end()) {
                        int idx = static_cast<int>(quadVerts.size());
                        gridToVertex[key] = idx;
                        quadVerts.push_back(origPos);
                    }
                    
                    // Also store triangle-local mapping for building quads
                    triGridToVertex[std::make_tuple(f, iu, iv)] = gridToVertex[key];
                }
            }
        }
    }
    
    std::cout << "Found " << quadVerts.size() << " quad vertices (integer grid points)\n";
    
    if (quadVerts.empty()) {
        std::cerr << "extractQuadMesh: No quad vertices found\n";
        return false;
    }
    
    // Build quad faces: for each integer grid cell (iu, iv) -> (iu+1, iv) -> (iu+1, iv+1) -> (iu, iv+1)
    // A quad exists if all four corners are present
    std::vector<Eigen::Vector4i> quads;
    
    for (const auto& kv : gridToVertex) {
        int iu = kv.first.first;
        int iv = kv.first.second;
        
        // Check if all four corners of the cell exist
        auto c00 = gridToVertex.find({iu, iv});
        auto c10 = gridToVertex.find({iu + 1, iv});
        auto c11 = gridToVertex.find({iu + 1, iv + 1});
        auto c01 = gridToVertex.find({iu, iv + 1});
        
        if (c00 != gridToVertex.end() && c10 != gridToVertex.end() && 
            c11 != gridToVertex.end() && c01 != gridToVertex.end()) {
            // Only add this quad if we haven't seen it yet
            // (we iterate from (iu, iv), so only add when this is the bottom-left corner)
            Eigen::Vector4i quad(c00->second, c10->second, c11->second, c01->second);
            
            // Avoid duplicates by only adding from the bottom-left corner
            bool alreadyAdded = false;
            for (const auto& q : quads) {
                if ((q(0) == quad(0) && q(1) == quad(1) && q(2) == quad(2) && q(3) == quad(3)) ||
                    (q(0) == quad(1) && q(1) == quad(2) && q(2) == quad(3) && q(3) == quad(0)) ||
                    (q(0) == quad(2) && q(1) == quad(3) && q(2) == quad(0) && q(3) == quad(1)) ||
                    (q(0) == quad(3) && q(1) == quad(0) && q(2) == quad(1) && q(3) == quad(2))) {
                    alreadyAdded = true;
                    break;
                }
            }
            if (!alreadyAdded) {
                quads.push_back(quad);
            }
        }
    }
    
    std::cout << "Built " << quads.size() << " quad faces\n";
    
    // Convert to Eigen matrices
    quadVertices_.resize(quadVerts.size(), 2);
    for (size_t i = 0; i < quadVerts.size(); ++i) {
        quadVertices_(i, 0) = quadVerts[i](0);
        quadVertices_(i, 1) = quadVerts[i](1);
    }
    
    quadFaces_.resize(quads.size(), 4);
    for (size_t i = 0; i < quads.size(); ++i) {
        quadFaces_(i, 0) = quads[i](0);
        quadFaces_(i, 1) = quads[i](1);
        quadFaces_(i, 2) = quads[i](2);
        quadFaces_(i, 3) = quads[i](3);
    }
    
    return true;
}

bool MIQSolver::writeQuadMesh(const std::string& filename) const {
    if (quadVertices_.rows() == 0 || quadFaces_.rows() == 0) {
        std::cerr << "writeQuadMesh: No quad mesh data. Call extractQuadMesh() first.\n";
        return false;
    }
    
    std::ofstream out(filename);
    if (!out) return false;
    
    // Write vertices (as 3D with z=0)
    for (int i = 0; i < static_cast<int>(quadVertices_.rows()); ++i) {
        out << "v " << quadVertices_(i, 0) << " " << quadVertices_(i, 1) << " 0\n";
    }
    
    // Write quad faces
    for (int f = 0; f < static_cast<int>(quadFaces_.rows()); ++f) {
        out << "f " << (quadFaces_(f, 0) + 1) << " " << (quadFaces_(f, 1) + 1) 
            << " " << (quadFaces_(f, 2) + 1) << " " << (quadFaces_(f, 3) + 1) << "\n";
    }
    
    return true;
}