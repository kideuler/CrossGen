#include "IntegerGridMap.hxx"

#include <algorithm>
#include <cmath>
#include <deque>
#include <fstream>
#include <iostream>
#include <limits>
#include <unordered_map>

#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>

namespace {

inline double sqr(double x) { return x * x; }

inline double signed_area2(const Point &a, const Point &b, const Point &c) {
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}

// Compute gradients of barycentric basis functions for triangle (p0,p1,p2).
// Returns 2x3 matrix G such that grad(f) = G * [f0 f1 f2]^T.
inline bool tri_grad_matrix(const Point &p0, const Point &p1, const Point &p2, 
                           Eigen::Matrix<double,2,3> &G, double &area) {
    const double det = signed_area2(p0, p1, p2); // 2 * signed area
    const double eps = 1e-14;
    if (std::abs(det) < eps) {
        area = 0.0;
        G.setZero();
        return false;
    }
    area = 0.5 * std::abs(det);
    const double invDet = 1.0 / det;

    // grad phi0
    G(0,0) = (p1[1] - p2[1]) * invDet;
    G(1,0) = (p2[0] - p1[0]) * invDet;
    // grad phi1
    G(0,1) = (p2[1] - p0[1]) * invDet;
    G(1,1) = (p0[0] - p2[0]) * invDet;
    // grad phi2
    G(0,2) = (p0[1] - p1[1]) * invDet;
    G(1,2) = (p1[0] - p0[0]) * invDet;

    return true;
}

inline Eigen::Vector2d to_vec2(const Point &p) {
    return Eigen::Vector2d(p[0], p[1]);
}

inline Point to_point(const Eigen::Vector2d &v) {
    return Point{v.x(), v.y()};
}

} // namespace


Eigen::Vector2d IntegerGridMap::rot90Vec(const Eigen::Vector2d &v, int k) {
    const int kk = ((k % 4) + 4) % 4;
    switch (kk) {
        case 0: return v;
        case 1: return Eigen::Vector2d(-v.y(), v.x());
        case 2: return Eigen::Vector2d(-v.x(), -v.y());
        case 3: return Eigen::Vector2d(v.y(), -v.x());
        default: return v;
    }
}

Eigen::Matrix2d IntegerGridMap::rot90Mat(int k) {
    const int kk = ((k % 4) + 4) % 4;
    Eigen::Matrix2d R;
    switch (kk) {
        case 0: R << 1, 0, 0, 1; break;
        case 1: R << 0, -1, 1, 0; break;
        case 2: R << -1, 0, 0, -1; break;
        case 3: R << 0, 1, -1, 0; break;
        default: R << 1, 0, 0, 1; break;
    }
    return R;
}


IntegerGridMap::IntegerGridMap(const CutMesh &cutMesh)
    : cm(cutMesh), cut(cutMesh.getCutMesh()), uField(cutMesh.getUField()) {
    triK.assign(cut.triangles.size(), -1);
    buildTriangleOrientations();
    buildSeams();
    chooseSingularityVertices();
    uvCoords.assign(cut.vertices.size(), Point{0.0, 0.0});
}


void IntegerGridMap::buildTriangleOrientations() {
    const int T = static_cast<int>(cut.triangles.size());
    if (T == 0) return;

    // Normalize raw u direction per triangle.
    std::vector<Eigen::Vector2d> d(T);
    for (int t = 0; t < T; ++t) {
        Eigen::Vector2d v = to_vec2(uField[t]);
        const double n = v.norm();
        if (n > 1e-14) v /= n;
        else v = Eigen::Vector2d(1.0, 0.0);
        d[t] = v;
    }

    triK.assign(T, -1);
    std::deque<int> q;
    triK[0] = 0;
    q.push_back(0);

    while (!q.empty()) {
        const int t = q.front();
        q.pop_front();
        const Eigen::Vector2d uAxis_t = rot90Vec(d[t], triK[t]);

        for (int e = 0; e < 3; ++e) {
            const int nt = cut.triangleAdjacency[t][e];
            if (nt < 0) continue;
            if (triK[nt] != -1) continue;

            // Choose k for neighbor that best aligns its u-axis with uAxis_t.
            int bestK = 0;
            double bestDot = -std::numeric_limits<double>::infinity();
            for (int k = 0; k < 4; ++k) {
                const Eigen::Vector2d cand = rot90Vec(d[nt], k);
                const double dot = uAxis_t.dot(cand);
                if (dot > bestDot) {
                    bestDot = dot;
                    bestK = k;
                }
            }
            triK[nt] = bestK;
            q.push_back(nt);
        }
    }

    // Any unvisited triangles (should not happen) get a default orientation.
    for (int t = 0; t < T; ++t) {
        if (triK[t] == -1) triK[t] = 0;
    }
}


void IntegerGridMap::buildSeams() {
    seams.clear();

    const auto &cutEdges = cm.getCutEdges();
    const auto &cutToOrig = cm.getCutVertexToOriginal();

    std::unordered_map<CutMesh::EdgeKey, std::vector<BoundaryEdge>, CutMesh::EdgeKeyHash> groups;
    groups.reserve(cutEdges.size() * 2 + 8);

    const int T = static_cast<int>(cut.triangles.size());
    for (int t = 0; t < T; ++t) {
        for (int e = 0; e < 3; ++e) {
            if (cut.triangleAdjacency[t][e] != -1) continue; // not boundary
            const int v0 = cut.triangles[t][e];
            const int v1 = cut.triangles[t][(e + 1) % 3];
            const int o0 = cutToOrig[v0];
            const int o1 = cutToOrig[v1];
            const CutMesh::EdgeKey key(o0, o1);
            if (cutEdges.find(key) == cutEdges.end()) continue; // original boundary edge
            BoundaryEdge be;
            be.tri = t;
            be.localEdge = e;
            be.v0 = v0;
            be.v1 = v1;
            be.o0 = o0;
            be.o1 = o1;
            groups[key].push_back(be);
        }
    }

    // Build paired seams.
    seams.reserve(groups.size());
    const int V = static_cast<int>(cut.vertices.size());
    int seamIndex = 0;
    for (auto &kv : groups) {
        auto &vec = kv.second;
        if (vec.size() != 2) {
            // In a valid explicit cut mesh, each cut edge should appear exactly
            // twice on the boundary. If not, skip defensively.
            continue;
        }
        const BoundaryEdge &e0 = vec[0];
        const BoundaryEdge &e1 = vec[1];

        Seam s;
        s.tri0 = e0.tri;
        s.tri1 = e1.tri;

        // Match endpoints by original vertex ids.
        if (e0.o0 == e1.o0 && e0.o1 == e1.o1) {
            s.v0a = e0.v0; s.v0b = e0.v1;
            s.v1a = e1.v0; s.v1b = e1.v1;
        } else if (e0.o0 == e1.o1 && e0.o1 == e1.o0) {
            s.v0a = e0.v0; s.v0b = e0.v1;
            s.v1a = e1.v1; s.v1b = e1.v0;
        } else {
            // Should not happen.
            continue;
        }

        // Rotation in the transition function.
        // Our triK encodes which cross arm is chosen as the local u-axis.
        // Frame of tri0 = R^{k0}(base), Frame of tri1 = R^{k1}(base).
        // To transform coordinates from side 0's frame to side 1's frame,
        // we need R^{k1 - k0} = R^{-(k0 - k1)}.
        // The constraint is: uv(v1) = R * uv(v0) + t
        // where R rotates from tri0's frame to tri1's frame.
        const int k0 = triK[s.tri0];
        const int k1 = triK[s.tri1];
        s.rot = ((k1 - k0) % 4 + 4) % 4;

        // Assign variable indices for (j,k) translations.
        const int base = 2 * V + 2 * seamIndex;
        s.jVar = base;
        s.kVar = base + 1;
        seams.push_back(s);
        ++seamIndex;
    }
}


void IntegerGridMap::chooseSingularityVertices() {
    singularCutVerts.clear();
    const auto &sing = cm.getSingularities();
    const auto &origToCut = cm.getOriginalToCutVertices();
    singularCutVerts.reserve(sing.size());

    for (const auto &s : sing) {
        const int ov = s.first;
        if (ov < 0 || ov >= static_cast<int>(origToCut.size())) continue;
        const auto &copies = origToCut[ov];
        if (copies.empty()) continue;
        // Pick the first copy; CutMesh ensures singularities lie on boundary.
        singularCutVerts.push_back(copies.front());
    }
}


void IntegerGridMap::assembleEnergy(double h,
                                   const std::vector<double> &weights,
                                   double translationReg,
                                   double softSeamWeight,
                                   Eigen::SparseMatrix<double> &H,
                                   Eigen::VectorXd &g) const {
    const int V = static_cast<int>(cut.vertices.size());
    const int T = static_cast<int>(cut.triangles.size());
    const int S = static_cast<int>(seams.size());
    const int N = 2 * V + 2 * S;

    H.resize(N, N);
    g.setZero(N);

    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(T * 18 * 2 + 2 * S + 8 * S);

    // Normalize raw u direction per triangle.
    std::vector<Eigen::Vector2d> d(T);
    for (int t = 0; t < T; ++t) {
        Eigen::Vector2d v = to_vec2(uField[t]);
        const double n = v.norm();
        if (n > 1e-14) v /= n;
        else v = Eigen::Vector2d(1.0, 0.0);
        d[t] = v;
    }

    const double hh = h;
    for (int t = 0; t < T; ++t) {
        const auto tri = cut.triangles[t];
        const int i0 = tri[0];
        const int i1 = tri[1];
        const int i2 = tri[2];
        const Point &p0 = cut.vertices[i0];
        const Point &p1 = cut.vertices[i1];
        const Point &p2 = cut.vertices[i2];

        Eigen::Matrix<double,2,3> G;
        double area;
        
        if (!tri_grad_matrix(p0, p1, p2, G, area) || area <= 0.0) continue;

        const double w = (t < static_cast<int>(weights.size()) ? weights[t] : 1.0);
        const double s = w * area;

        // Oriented target vectors.
        // vAxis is always CCW 90° rotation of uAxis to ensure positive Jacobian determinant.
        const Eigen::Vector2d uAxis = rot90Vec(d[t], triK[t]);
        const Eigen::Vector2d vAxis = rot90Vec(uAxis, 1);  // CCW rotation

        // Local stiffness K = (h^2 * s) * (G^T G)
        Eigen::Matrix3d K = Eigen::Matrix3d::Zero();
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                const double dot = G.col(a).dot(G.col(b));
                K(a,b) = (hh * hh) * s * dot;
            }
        }

        const int idx[3] = { i0, i1, i2 };
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                if (K(a,b) == 0.0) continue;
                // u block
                trip.emplace_back(idx[a], idx[b], K(a,b));
                // v block
                trip.emplace_back(idx[a] + V, idx[b] + V, K(a,b));
            }
        }

        // RHS: (h * s) * (G^T target)
        Eigen::Vector3d Fu, Fv;
        for (int a = 0; a < 3; ++a) {
            Fu[a] = hh * s * G.col(a).dot(uAxis);
            Fv[a] = hh * s * G.col(a).dot(vAxis);
        }

        g[idx[0]] += Fu[0];
        g[idx[1]] += Fu[1];
        g[idx[2]] += Fu[2];
        g[idx[0] + V] += Fv[0];
        g[idx[1] + V] += Fv[1];
        g[idx[2] + V] += Fv[2];
    }

    // Soft seam constraints: instead of hard constraints, add penalty terms
    // This allows slight violation to prevent flips, while encouraging gluing
    // The constraint is: u(v1) = R*uv(v0) + (j,k)
    // Energy: w_seam * |u(v1) - R*uv(v0) - j|^2 + w_seam * |v(v1) - R*uv(v0) - k|^2
    if (softSeamWeight > 0.0) {
        for (const Seam &s : seams) {
            const Eigen::Matrix2d R = rot90Mat(s.rot);
            
            auto addSoftConstraint = [&](int v0, int v1) {
                // Constraint: u(v1) - R00*u(v0) - R01*v(v0) - j = 0
                // Energy: w * (u(v1) - R00*u(v0) - R01*v(v0) - j)^2
                // Expanded: w * [u(v1)^2 + R00^2*u(v0)^2 + R01^2*v(v0)^2 + j^2
                //              - 2*u(v1)*R00*u(v0) - 2*u(v1)*R01*v(v0) - 2*u(v1)*j
                //              + 2*R00*u(v0)*R01*v(v0) + 2*R00*u(v0)*j + 2*R01*v(v0)*j]
                const double w = softSeamWeight;
                const double R00 = R(0,0), R01 = R(0,1), R10 = R(1,0), R11 = R(1,1);
                
                // u equation: u(v1) - R00*u(v0) - R01*v(v0) - j = 0
                trip.emplace_back(v1, v1, w);                        // u(v1)^2
                trip.emplace_back(v0, v0, w * R00 * R00);            // R00^2 * u(v0)^2
                trip.emplace_back(v0 + V, v0 + V, w * R01 * R01);    // R01^2 * v(v0)^2
                trip.emplace_back(s.jVar, s.jVar, w);                // j^2
                
                trip.emplace_back(v1, v0, -w * R00);                 // -2*u(v1)*R00*u(v0) / 2
                trip.emplace_back(v0, v1, -w * R00);
                trip.emplace_back(v1, v0 + V, -w * R01);             // -2*u(v1)*R01*v(v0) / 2
                trip.emplace_back(v0 + V, v1, -w * R01);
                trip.emplace_back(v1, s.jVar, -w);                   // -2*u(v1)*j / 2
                trip.emplace_back(s.jVar, v1, -w);
                
                trip.emplace_back(v0, v0 + V, w * R00 * R01);        // 2*R00*u(v0)*R01*v(v0) / 2
                trip.emplace_back(v0 + V, v0, w * R00 * R01);
                trip.emplace_back(v0, s.jVar, w * R00);              // 2*R00*u(v0)*j / 2
                trip.emplace_back(s.jVar, v0, w * R00);
                trip.emplace_back(v0 + V, s.jVar, w * R01);          // 2*R01*v(v0)*j / 2
                trip.emplace_back(s.jVar, v0 + V, w * R01);
                
                // v equation: v(v1) - R10*u(v0) - R11*v(v0) - k = 0
                trip.emplace_back(v1 + V, v1 + V, w);                // v(v1)^2
                trip.emplace_back(v0, v0, w * R10 * R10);            // R10^2 * u(v0)^2  
                trip.emplace_back(v0 + V, v0 + V, w * R11 * R11);    // R11^2 * v(v0)^2
                trip.emplace_back(s.kVar, s.kVar, w);                // k^2
                
                trip.emplace_back(v1 + V, v0, -w * R10);
                trip.emplace_back(v0, v1 + V, -w * R10);
                trip.emplace_back(v1 + V, v0 + V, -w * R11);
                trip.emplace_back(v0 + V, v1 + V, -w * R11);
                trip.emplace_back(v1 + V, s.kVar, -w);
                trip.emplace_back(s.kVar, v1 + V, -w);
                
                trip.emplace_back(v0, v0 + V, w * R10 * R11);
                trip.emplace_back(v0 + V, v0, w * R10 * R11);
                trip.emplace_back(v0, s.kVar, w * R10);
                trip.emplace_back(s.kVar, v0, w * R10);
                trip.emplace_back(v0 + V, s.kVar, w * R11);
                trip.emplace_back(s.kVar, v0 + V, w * R11);
            };
            
            addSoftConstraint(s.v0a, s.v1a);
            addSoftConstraint(s.v0b, s.v1b);
        }
    }

    // Small regularization on translation variables helps conditioning.
    if (translationReg > 0.0) {
        for (int s = 0; s < S; ++s) {
            trip.emplace_back(seams[s].jVar, seams[s].jVar, translationReg);
            trip.emplace_back(seams[s].kVar, seams[s].kVar, translationReg);
        }
    }

    H.setFromTriplets(trip.begin(), trip.end());
    H.makeCompressed();
}


void IntegerGridMap::assembleConstraints(const std::vector<std::pair<int,double>> &fixed,
                                        double softSeamWeight,
                                        Eigen::SparseMatrix<double> &C,
                                        Eigen::VectorXd &d) const {
    const int V = static_cast<int>(cut.vertices.size());
    const int S = static_cast<int>(seams.size());
    const int N = 2 * V + 2 * S;

    // If soft seam constraints are used, skip hard seam constraints
    const bool useHardSeams = (softSeamWeight <= 0.0);
    const int seamCount = useHardSeams ? static_cast<int>(seams.size()) : 0;
    const int baseRows = 4 * seamCount + 2; // seams + anchor
    const int rows = baseRows + static_cast<int>(fixed.size());

    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve( (4 * seamCount) * 4 + 2 + fixed.size() );
    d.setZero(rows);

    int r = 0;
    if (useHardSeams) {
    for (const Seam &s : seams) {
        const Eigen::Matrix2d R = rot90Mat(s.rot);

        auto add_point_constraints = [&](int v0, int v1) {
            // u(v1) - R00*u(v0) - R01*v(v0) - j = 0
            trip.emplace_back(r, v1, 1.0);
            trip.emplace_back(r, v0, -R(0,0));
            trip.emplace_back(r, v0 + V, -R(0,1));
            trip.emplace_back(r, s.jVar, -1.0);
            ++r;
            // v(v1) - R10*u(v0) - R11*v(v0) - k = 0
            trip.emplace_back(r, v1 + V, 1.0);
            trip.emplace_back(r, v0, -R(1,0));
            trip.emplace_back(r, v0 + V, -R(1,1));
            trip.emplace_back(r, s.kVar, -1.0);
            ++r;
        };

        add_point_constraints(s.v0a, s.v1a);
        add_point_constraints(s.v0b, s.v1b);
    }
    }

    // Anchor to remove the translational nullspace.
    const int anchor = 0;
    trip.emplace_back(r, anchor, 1.0);
    d[r] = 0.0;
    ++r;
    trip.emplace_back(r, anchor + V, 1.0);
    d[r] = 0.0;
    ++r;

    // Fixed-variable constraints (for rounding): x[idx] = value.
    for (const auto &fv : fixed) {
        const int idx = fv.first;
        if (idx < 0 || idx >= N) continue;
        trip.emplace_back(r, idx, 1.0);
        d[r] = fv.second;
        ++r;
    }

    C.resize(r, N);
    d.conservativeResize(r);
    C.setFromTriplets(trip.begin(), trip.end());
    C.makeCompressed();
}


bool IntegerGridMap::solveKKT(const Eigen::SparseMatrix<double> &H,
                             const Eigen::VectorXd &g,
                             const Eigen::SparseMatrix<double> &C,
                             const Eigen::VectorXd &d,
                             Eigen::VectorXd &x) const {
    const int N = static_cast<int>(H.rows());
    const int P = static_cast<int>(C.rows());
    const int K = N + P;

    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(H.nonZeros() + 2 * C.nonZeros() + K);

    // H block with small diagonal regularization to ensure positive definiteness
    const double reg = 1e-8;
    for (int k = 0; k < H.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(H, k); it; ++it) {
            trip.emplace_back(it.row(), it.col(), it.value());
        }
    }
    // Add regularization to H diagonal
    for (int i = 0; i < N; ++i) {
        trip.emplace_back(i, i, reg);
    }

    // C^T and C blocks
    for (int k = 0; k < C.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(C, k); it; ++it) {
            const int r = it.row();
            const int c = it.col();
            const double v = it.value();
            // top-right: C^T
            trip.emplace_back(c, N + r, v);
            // bottom-left: C
            trip.emplace_back(N + r, c, v);
        }
    }

    // Add small negative regularization to the (2,2) block to stabilize the saddle-point system.
    for (int i = 0; i < P; ++i) {
        trip.emplace_back(N + i, N + i, -reg);
    }

    Eigen::SparseMatrix<double> KKT(K, K);
    KKT.setFromTriplets(trip.begin(), trip.end());
    KKT.makeCompressed();

    Eigen::VectorXd rhs(K);
    rhs.head(N) = g;
    rhs.tail(P) = d;

    // Use SparseLU with the regularized system
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(KKT);
    solver.factorize(KKT);
    if (solver.info() != Eigen::Success) {
        std::cerr << "[IntegerGridMap] KKT factorization failed." << std::endl;
        return false;
    }

    Eigen::VectorXd sol = solver.solve(rhs);
    if (solver.info() != Eigen::Success) {
        std::cerr << "[IntegerGridMap] KKT solve failed." << std::endl;
        return false;
    }

    x = sol.head(N);
    return true;
}


IntegerGridMap::Stats IntegerGridMap::computeStatsFromX(const Eigen::VectorXd &x, double h, bool /*verbose*/) const {
    Stats st;
    const int V = static_cast<int>(cut.vertices.size());
    const int T = static_cast<int>(cut.triangles.size());

    double sumLambda = 0.0;
    int count = 0;
    int flipped = 0;
    double maxLambda = 0.0;

    for (int t = 0; t < T; ++t) {
        const auto tri = cut.triangles[t];
        const int i0 = tri[0], i1 = tri[1], i2 = tri[2];
        const Point &p0 = cut.vertices[i0];
        const Point &p1 = cut.vertices[i1];
        const Point &p2 = cut.vertices[i2];

        Eigen::Matrix<double,2,3> G;
        double area;
        
        if (!tri_grad_matrix(p0, p1, p2, G, area) || area <= 0.0) continue;

        const double u0 = x[i0];
        const double u1 = x[i1];
        const double u2 = x[i2];
        const double v0 = x[i0 + V];
        const double v1 = x[i1 + V];
        const double v2 = x[i2 + V];

        const Eigen::Vector2d gradu = G.col(0) * u0 + G.col(1) * u1 + G.col(2) * u2;
        const Eigen::Vector2d gradv = G.col(0) * v0 + G.col(1) * v1 + G.col(2) * v2;

        Eigen::Matrix2d J;
        J.row(0) = gradu.transpose();
        J.row(1) = gradv.transpose();

        // The Jacobian determinant tells us if the mapping preserves orientation.
        // det(J) = gradu_x * gradv_y - gradu_y * gradv_x
        // For MIQ, we want det(J) > 0 (orientation-preserving parameterization).
        const double detJ = J.determinant();
        if (detJ <= 0.0) {
            ++flipped;
        }

        const Eigen::JacobiSVD<Eigen::Matrix2d> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);
        const Eigen::Vector2d sig = svd.singularValues();
        
        const double tau = (detJ >= 0.0) ? 1.0 : -1.0;

        const double lambda = std::abs(tau * sig[0] / h - 1.0) + std::abs(tau * sig[1] / h - 1.0);
        sumLambda += lambda;
        maxLambda = std::max(maxLambda, lambda);
        ++count;
    }

    st.flipped_triangles = flipped;
    st.max_lambda = maxLambda;
    st.mean_lambda = (count > 0) ? (sumLambda / static_cast<double>(count)) : 0.0;
    return st;
}


void IntegerGridMap::updateStiffeningWeights(const Eigen::VectorXd &x,
                                            double h,
                                            double c,
                                            double dClamp,
                                            int smoothSteps,
                                            std::vector<double> &weights) const {
    const int T = static_cast<int>(cut.triangles.size());
    const int V = static_cast<int>(cut.vertices.size());

    std::vector<double> lambda(T, 0.0);
    // Compute per-triangle lambda.
    for (int t = 0; t < T; ++t) {
        const auto tri = cut.triangles[t];
        const int i0 = tri[0], i1 = tri[1], i2 = tri[2];
        const Point &p0 = cut.vertices[i0];
        const Point &p1 = cut.vertices[i1];
        const Point &p2 = cut.vertices[i2];

        Eigen::Matrix<double,2,3> G;
        double area;
        
        if (!tri_grad_matrix(p0, p1, p2, G, area) || area <= 0.0) {
            lambda[t] = 0.0;
            continue;
        }

        const double u0 = x[i0];
        const double u1 = x[i1];
        const double u2 = x[i2];
        const double v0 = x[i0 + V];
        const double v1 = x[i1 + V];
        const double v2 = x[i2 + V];

        const Eigen::Vector2d gradu = G.col(0) * u0 + G.col(1) * u1 + G.col(2) * u2;
        const Eigen::Vector2d gradv = G.col(0) * v0 + G.col(1) * v1 + G.col(2) * v2;

        Eigen::Matrix2d J;
        J.row(0) = gradu.transpose();
        J.row(1) = gradv.transpose();

        const double detJ = J.determinant();
        const Eigen::JacobiSVD<Eigen::Matrix2d> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);
        const Eigen::Vector2d sig = svd.singularValues();

        const double tau = (detJ >= 0.0) ? 1.0 : -1.0;
        
        // For flipped triangles (detJ < 0), use a large penalty to encourage fixing them
        if (tau < 0.0) {
            lambda[t] = 10.0;  // Large distortion value for flipped triangles
        } else {
            lambda[t] = std::abs(tau * sig[0] / h - 1.0) + std::abs(tau * sig[1] / h - 1.0);
        }
    }

    // Uniform Laplacian on dual mesh.
    std::vector<double> delta(T, 0.0);
    for (int t = 0; t < T; ++t) {
        double acc = 0.0;
        int deg = 0;
        for (int e = 0; e < 3; ++e) {
            const int nt = cut.triangleAdjacency[t][e];
            if (nt < 0) continue;
            acc += (lambda[nt] - lambda[t]);
            ++deg;
        }
        delta[t] = acc; // sum(nei - self)
    }

    if (weights.size() != static_cast<size_t>(T)) weights.assign(T, 1.0);

    for (int t = 0; t < T; ++t) {
        // Increase weight more aggressively for triangles with high distortion
        // and their neighbors
        const double upd = std::min(c * std::abs(4.0 * delta[t]), dClamp);
        weights[t] += upd;
        
        // Also directly increase weights on flipped triangles
        if (lambda[t] >= 10.0) {
            weights[t] += dClamp;
        }
    }

    // Smooth weights a few steps.
    std::vector<double> tmp(T, 1.0);
    for (int it = 0; it < smoothSteps; ++it) {
        for (int t = 0; t < T; ++t) {
            double sum = weights[t];
            int cnt = 1;
            for (int e = 0; e < 3; ++e) {
                const int nt = cut.triangleAdjacency[t][e];
                if (nt < 0) continue;
                sum += weights[nt];
                ++cnt;
            }
            tmp[t] = sum / static_cast<double>(cnt);
        }
        weights.swap(tmp);
    }
}


// Untangle flipped triangles using gradient descent with a barrier function.
// This moves free vertices to eliminate flipped triangles while respecting fixed constraints.
bool IntegerGridMap::untangleUV(Eigen::VectorXd &x,
                               const std::vector<std::pair<int,double>> &fixed,
                               double h,
                               int maxIter) const {
    const int V = static_cast<int>(cut.vertices.size());
    const int T = static_cast<int>(cut.triangles.size());
    
    // Build set of fixed variable indices
    std::vector<char> isFixed(2 * V, 0);
    for (const auto &fv : fixed) {
        if (fv.first >= 0 && fv.first < 2 * V) {
            isFixed[fv.first] = 1;
        }
    }
    
    // Count initial flips
    auto countFlips = [&]() {
        int flips = 0;
        for (int t = 0; t < T; ++t) {
            const auto tri = cut.triangles[t];
            const int i0 = tri[0], i1 = tri[1], i2 = tri[2];
            const double u0 = x[i0], u1 = x[i1], u2 = x[i2];
            const double v0 = x[i0+V], v1 = x[i1+V], v2 = x[i2+V];
            // Signed area in UV space
            const double area = (u1-u0)*(v2-v0) - (u2-u0)*(v1-v0);
            if (area <= 0.0) ++flips;
        }
        return flips;
    };
    
    int initialFlips = countFlips();
    if (initialFlips == 0) return true;
    
    // Gradient descent with barrier
    const double eps = 1e-8;
    double stepSize = 0.1;
    
    for (int iter = 0; iter < maxIter; ++iter) {
        Eigen::VectorXd grad = Eigen::VectorXd::Zero(2 * V);
        
        // Compute gradient of barrier function
        for (int t = 0; t < T; ++t) {
            const auto tri = cut.triangles[t];
            const int i0 = tri[0], i1 = tri[1], i2 = tri[2];
            const double u0 = x[i0], u1 = x[i1], u2 = x[i2];
            const double v0 = x[i0+V], v1 = x[i1+V], v2 = x[i2+V];
            
            // Signed area in UV space: A = (u1-u0)*(v2-v0) - (u2-u0)*(v1-v0)
            const double area = (u1-u0)*(v2-v0) - (u2-u0)*(v1-v0);
            
            // Barrier: -log(area) for area > eps, quadratic penalty otherwise
            double weight;
            if (area > eps) {
                weight = -1.0 / area;  // derivative of -log(area)
            } else {
                // Strong repulsion for nearly-degenerate or flipped triangles
                weight = -1000.0 * (eps - area + 1.0);
            }
            
            // Gradient of area w.r.t. each vertex
            // dA/du0 = -(v2-v1), dA/dv0 = (u2-u1)
            // dA/du1 = (v2-v0),  dA/dv1 = -(u2-u0)
            // dA/du2 = -(v1-v0), dA/dv2 = (u1-u0)
            grad[i0]   += weight * (-(v2-v1));
            grad[i0+V] += weight * (u2-u1);
            grad[i1]   += weight * (v2-v0);
            grad[i1+V] += weight * (-(u2-u0));
            grad[i2]   += weight * (-(v1-v0));
            grad[i2+V] += weight * (u1-u0);
        }
        
        // Zero out gradient for fixed variables
        for (int i = 0; i < 2*V; ++i) {
            if (isFixed[i]) grad[i] = 0.0;
        }
        
        // Line search
        double gradNorm = grad.norm();
        if (gradNorm < 1e-10) break;
        
        Eigen::VectorXd dir = -grad / gradNorm;
        
        // Backtracking line search
        double alpha = stepSize;
        Eigen::VectorXd xNew = x + alpha * dir;
        int newFlips = 0;
        for (int t = 0; t < T; ++t) {
            const auto tri = cut.triangles[t];
            const int i0 = tri[0], i1 = tri[1], i2 = tri[2];
            const double u0 = xNew[i0], u1 = xNew[i1], u2 = xNew[i2];
            const double v0 = xNew[i0+V], v1 = xNew[i1+V], v2 = xNew[i2+V];
            const double area = (u1-u0)*(v2-v0) - (u2-u0)*(v1-v0);
            if (area <= 0.0) ++newFlips;
        }
        
        // Reduce step if we didn't improve
        int currentFlips = countFlips();
        while (newFlips >= currentFlips && alpha > 1e-8) {
            alpha *= 0.5;
            xNew = x + alpha * dir;
            newFlips = 0;
            for (int t = 0; t < T; ++t) {
                const auto tri = cut.triangles[t];
                const int i0 = tri[0], i1 = tri[1], i2 = tri[2];
                const double u0 = xNew[i0], u1 = xNew[i1], u2 = xNew[i2];
                const double v0 = xNew[i0+V], v1 = xNew[i1+V], v2 = xNew[i2+V];
                const double area = (u1-u0)*(v2-v0) - (u2-u0)*(v1-v0);
                if (area <= 0.0) ++newFlips;
            }
        }
        
        if (newFlips < currentFlips) {
            x = xNew;
            if (newFlips == 0) return true;
        } else {
            stepSize *= 0.5;
            if (stepSize < 1e-10) break;
        }
    }
    
    return countFlips() == 0;
}


void IntegerGridMap::extractUVFromX(const Eigen::VectorXd &x) {
    const int V = static_cast<int>(cut.vertices.size());
    uvCoords.resize(V);
    for (int i = 0; i < V; ++i) {
        uvCoords[i][0] = x[i];
        uvCoords[i][1] = x[i + V];
    }
}


bool IntegerGridMap::solve(const Options &opt) {
    const int V = static_cast<int>(cut.vertices.size());
    const int S = static_cast<int>(seams.size());
    const int N = 2 * V + 2 * S;

    if (V == 0 || cut.triangles.empty()) {
        std::cerr << "[IntegerGridMap] Empty mesh." << std::endl;
        return false;
    }

    // Compute h based on mesh bounding box if not specified.
    double h = opt.h;
    if (h <= 0.0) {
        // Find bounding box of the mesh
        double minX = std::numeric_limits<double>::max();
        double maxX = std::numeric_limits<double>::lowest();
        double minY = std::numeric_limits<double>::max();
        double maxY = std::numeric_limits<double>::lowest();
        for (const auto& v : cut.vertices) {
            minX = std::min(minX, v[0]);
            maxX = std::max(maxX, v[0]);
            minY = std::min(minY, v[1]);
            maxY = std::max(maxY, v[1]);
        }
        double rangeX = maxX - minX;
        double rangeY = maxY - minY;
        double maxRange = std::max(rangeX, rangeY);
        // h is chosen so that target_grid_lines integer grid lines span the mesh
        h = maxRange / static_cast<double>(opt.target_grid_lines);
        std::cout << "[IntegerGridMap] Mesh bbox: [" << minX << "," << maxX << "] x [" << minY << "," << maxY << "]"
                  << ", range = " << maxRange << ", auto h = " << h << std::endl;
    }

    // Store the effective h for later use by computeStats
    effectiveH = h;

    std::vector<double> weights(cut.triangles.size(), 1.0);

    Eigen::SparseMatrix<double> H;
    Eigen::VectorXd g;
    Eigen::VectorXd x;

    // --- Local stiffening loop (continuous variables) ---
    for (int it = 0; it < std::max(1, opt.stiffening_iterations); ++it) {
        assembleEnergy(h, weights, opt.translation_regularization, opt.soft_seam_weight, H, g);
        Eigen::SparseMatrix<double> C;
        Eigen::VectorXd d;
        assembleConstraints({}, opt.soft_seam_weight, C, d);

        if (!solveKKT(H, g, C, d, x)) return false;

        Stats st = computeStatsFromX(x, h);
        // Update weights except at last iteration.
        if (it + 1 < opt.stiffening_iterations) {
            updateStiffeningWeights(x, h, opt.stiffening_c, opt.stiffening_d, opt.stiffening_smooth_steps, weights);
        }
        // Early out if already orientation-preserving.
        if (st.flipped_triangles == 0 && it >= 1) break;
    }

    // --- Mixed-integer rounding (greedy with quality check) ---
    // Following MIQ paper: round one variable at a time, choosing the one
    // closest to an integer, and re-solve after each rounding.
    std::vector<std::pair<int,double>> fixed;
    fixed.reserve(2 * seams.size());

    std::vector<int> integerVars;
    integerVars.reserve(2 * seams.size());

    // Seam translations are integers.
    for (const Seam &s : seams) {
        integerVars.push_back(s.jVar);
        integerVars.push_back(s.kVar);
    }

    // NOTE: We intentionally do NOT round singularity positions to integers.
    // While the MIQ paper suggests this for clean quad meshing, it severely
    // constrains the solution and often leads to flipped triangles.
    // The seam translations are sufficient for seamless parameterization.

    if (opt.do_rounding && !integerVars.empty()) {
        std::vector<char> isFixed(N, 0);
        for (const auto &fv : fixed) {
            if (fv.first >= 0 && fv.first < N) isFixed[fv.first] = 1;
        }

        // Extra pre-rounding stiffening to get variables closer to integers
        for (int preIt = 0; preIt < 10; ++preIt) {
            assembleEnergy(h, weights, opt.translation_regularization, opt.soft_seam_weight, H, g);
            Eigen::SparseMatrix<double> C;
            Eigen::VectorXd d;
            assembleConstraints(fixed, opt.soft_seam_weight, C, d);
            if (!solveKKT(H, g, C, d, x)) return false;
            
            Stats st = computeStatsFromX(x, h);
            if (st.flipped_triangles == 0) break;
            
            updateStiffeningWeights(x, h, opt.stiffening_c * 2.0, opt.stiffening_d * 2.0, opt.stiffening_smooth_steps, weights);
        }

        std::vector<int> remaining;
        remaining.reserve(integerVars.size());
        for (int idx : integerVars) {
            if (idx >= 0 && idx < N && !isFixed[idx]) remaining.push_back(idx);
        }

        // Get baseline quality
        Stats baseStats = computeStatsFromX(x, h);

        int roundIter = 0;
        while (!remaining.empty() && roundIter < opt.max_round_iterations) {
            // Find variable closest to an integer (smallest rounding error).
            int bestIdx = -1;
            double bestErr = std::numeric_limits<double>::infinity();
            
            for (int idx : remaining) {
                const double val = x[idx];
                const double r = std::round(val);
                const double err = std::abs(val - r);
                if (err < bestErr) {
                    bestErr = err;
                    bestIdx = idx;
                }
            }
            
            if (bestIdx < 0) break;
            
            // Try BOTH floor and ceil, pick the one with fewer flipped triangles
            const double val = x[bestIdx];
            const double floorVal = std::floor(val);
            const double ceilVal = std::ceil(val);
            
            struct RoundingOption {
                double roundedVal;
                Eigen::VectorXd resultX;
                Stats stats;
                bool valid;
            };
            
            std::vector<RoundingOption> options(2);
            options[0].roundedVal = floorVal;
            options[1].roundedVal = ceilVal;
            
            for (auto &roundOpt : options) {
                std::vector<std::pair<int,double>> testFixed = fixed;
                testFixed.emplace_back(bestIdx, roundOpt.roundedVal);
                
                Eigen::SparseMatrix<double> C;
                Eigen::VectorXd d;
                assembleConstraints(testFixed, opt.soft_seam_weight, C, d);
                
                roundOpt.valid = solveKKT(H, g, C, d, roundOpt.resultX);
                if (roundOpt.valid) {
                    roundOpt.stats = computeStatsFromX(roundOpt.resultX, h);
                }
            }
            
            // Find best valid option (fewest flipped triangles, then lowest distortion)
            int bestOption = -1;
            for (int i = 0; i < 2; ++i) {
                if (!options[i].valid) continue;
                if (bestOption < 0 ||
                    options[i].stats.flipped_triangles < options[bestOption].stats.flipped_triangles ||
                    (options[i].stats.flipped_triangles == options[bestOption].stats.flipped_triangles &&
                     options[i].stats.mean_lambda < options[bestOption].stats.mean_lambda)) {
                    bestOption = i;
                }
            }
            
            if (bestOption < 0) {
                // Neither option worked, skip this variable
                remaining.erase(std::remove(remaining.begin(), remaining.end(), bestIdx), remaining.end());
                ++roundIter;
                continue;
            }
            
            Stats testStats = options[bestOption].stats;
            double bestRounded = options[bestOption].roundedVal;
            Eigen::VectorXd testX = options[bestOption].resultX;
            
            // Accept the rounding if it doesn't increase flipped triangles much
            // Allow just 1 new flip
            const int maxNewFlips = baseStats.flipped_triangles + 1;
            
            if (testStats.flipped_triangles <= maxNewFlips || bestErr < 1e-6) {
                // Accept this rounding
                fixed.emplace_back(bestIdx, bestRounded);
                isFixed[bestIdx] = 1;
                x = testX;
                
                // Apply local stiffening after each snap to help recover from distortion
                if (testStats.flipped_triangles > 0) {
                    for (int stiffIt = 0; stiffIt < 5; ++stiffIt) {
                        updateStiffeningWeights(x, h, opt.stiffening_c * 2.0, opt.stiffening_d * 2.0, opt.stiffening_smooth_steps, weights);
                        assembleEnergy(h, weights, opt.translation_regularization, opt.soft_seam_weight, H, g);
                        Eigen::SparseMatrix<double> Cs;
                        Eigen::VectorXd ds;
                        assembleConstraints(fixed, opt.soft_seam_weight, Cs, ds);
                        if (!solveKKT(H, g, Cs, ds, x)) return false;
                        
                        Stats stiffStats = computeStatsFromX(x, h);
                        if (stiffStats.flipped_triangles == 0) break;
                    }
                }
                
                baseStats = computeStatsFromX(x, h);
                
                // Remove from remaining
                remaining.erase(std::remove(remaining.begin(), remaining.end(), bestIdx), remaining.end());
                
                // If error was very small, try to batch more variables with similar small errors
                if (bestErr < opt.round_keep_threshold) {
                    std::vector<int> toFixBatch;
                    for (int idx : remaining) {
                        const double bval = x[idx];
                        const double r = std::round(bval);
                        if (std::abs(bval - r) < opt.round_keep_threshold) {
                            toFixBatch.push_back(idx);
                        }
                    }
                    
                    // Fix all nearly-integer variables at once
                    for (int idx : toFixBatch) {
                        const double bval = x[idx];
                        const double r = std::round(bval);
                        fixed.emplace_back(idx, r);
                        isFixed[idx] = 1;
                    }
                    
                    // Remove from remaining
                    for (int idx : toFixBatch) {
                        remaining.erase(std::remove(remaining.begin(), remaining.end(), idx), remaining.end());
                    }
                    
                    if (!toFixBatch.empty()) {
                        // Re-solve with all newly fixed variables
                        Eigen::SparseMatrix<double> Cb;
                        Eigen::VectorXd db;
                        assembleConstraints(fixed, opt.soft_seam_weight, Cb, db);
                        if (!solveKKT(H, g, Cb, db, x)) return false;
                        
                        // Apply stiffening after batch fix too
                        for (int stiffIt = 0; stiffIt < 5; ++stiffIt) {
                            Stats batchStats = computeStatsFromX(x, h);
                            if (batchStats.flipped_triangles == 0) break;
                            updateStiffeningWeights(x, h, opt.stiffening_c * 2.0, opt.stiffening_d * 2.0, opt.stiffening_smooth_steps, weights);
                            assembleEnergy(h, weights, opt.translation_regularization, opt.soft_seam_weight, H, g);
                            assembleConstraints(fixed, opt.soft_seam_weight, Cb, db);
                            if (!solveKKT(H, g, Cb, db, x)) return false;
                        }
                        
                        baseStats = computeStatsFromX(x, h);
                    }
                }
            } else {
                // Both roundings increased flips too much, skip this variable
                remaining.erase(std::remove(remaining.begin(), remaining.end(), bestIdx), remaining.end());
            }

            ++roundIter;
        }

        // If some remain (due to max_round_iterations), round them all at once.
        // Use the rounding direction that's closest to current value.
        if (!remaining.empty()) {
            for (int idx : remaining) {
                fixed.emplace_back(idx, std::round(x[idx]));
            }
            Eigen::SparseMatrix<double> C;
            Eigen::VectorXd d;
            assembleConstraints(fixed, opt.soft_seam_weight, C, d);
            if (!solveKKT(H, g, C, d, x)) return false;
        }

        // More aggressive stiffening iterations with fixed integers to eliminate flipped triangles.
        // Increase stiffening weights more aggressively on flipped triangles.
        const int maxPostRoundIterations = 100;
        for (int it = 0; it < maxPostRoundIterations; ++it) {
            Stats st = computeStatsFromX(x, h);
            if (st.flipped_triangles == 0) break;
            
            // Use stronger stiffening for post-rounding, increasing with iteration
            const double stiffMult = 5.0 + it * 1.0;  // Stronger multiplier
            updateStiffeningWeights(x, h, opt.stiffening_c * stiffMult, opt.stiffening_d * stiffMult, opt.stiffening_smooth_steps, weights);
            assembleEnergy(h, weights, opt.translation_regularization, opt.soft_seam_weight, H, g);
            Eigen::SparseMatrix<double> C;
            Eigen::VectorXd d;
            assembleConstraints(fixed, opt.soft_seam_weight, C, d);
            if (!solveKKT(H, g, C, d, x)) return false;
        }
        
        // If still have flips after stiffening, try untangling with gradient descent
        Stats finalStats = computeStatsFromX(x, h);
        if (finalStats.flipped_triangles > 0) {
            // Build fixed list for UV variables only (exclude seam translations)
            std::vector<std::pair<int,double>> uvFixed;
            const int V = static_cast<int>(cut.vertices.size());
            for (const auto &fv : fixed) {
                if (fv.first < 2 * V) {  // Only UV variables, not seam translations
                    uvFixed.push_back(fv);
                }
            }
            untangleUV(x, uvFixed, h, 1000);
        }
    }

    // Store solution for validation
    solution = x;
    
    extractUVFromX(x);
    return true;
}

bool IntegerGridMap::solve() {
    return solve(Options());
}


IntegerGridMap::Stats IntegerGridMap::computeStats(double h) const {
    const int V = static_cast<int>(cut.vertices.size());
    if (uvCoords.size() != static_cast<size_t>(V)) return Stats{};

    // Use effectiveH if h is not specified (<=0)
    const double actualH = (h > 0.0) ? h : effectiveH;

    Eigen::VectorXd x(2 * V);
    for (int i = 0; i < V; ++i) {
        x[i] = uvCoords[i][0];
        x[i + V] = uvCoords[i][1];
    }
    return computeStatsFromX(x, actualH);
}


bool IntegerGridMap::writeOBJWithUV(const std::string &filename) const {
    std::ofstream out(filename);
    if (!out.is_open()) return false;

    const int V = static_cast<int>(cut.vertices.size());
    const int T = static_cast<int>(cut.triangles.size());

    for (int i = 0; i < V; ++i) {
        out << "v " << cut.vertices[i][0] << " " << cut.vertices[i][1] << " 0\n";
    }
    for (int i = 0; i < V; ++i) {
        out << "vt " << uvCoords[i][0] << " " << uvCoords[i][1] << "\n";
    }

    for (int t = 0; t < T; ++t) {
        const auto tri = cut.triangles[t];
        // OBJ is 1-based.
        out << "f "
            << tri[0] + 1 << "/" << tri[0] + 1 << " "
            << tri[1] + 1 << "/" << tri[1] + 1 << " "
            << tri[2] + 1 << "/" << tri[2] + 1 << "\n";
    }

    return true;
}


bool IntegerGridMap::writeOriginalMeshWithUV(const std::string &filename) const {
    std::ofstream out(filename);
    if (!out.is_open()) return false;

    const Mesh &orig = cm.getOriginalMesh();
    const auto &cutToOrig = cm.getCutVertexToOriginal();
    
    const int origV = static_cast<int>(orig.vertices.size());
    const int origT = static_cast<int>(orig.triangles.size());
    const int cutV = static_cast<int>(cut.vertices.size());
    const int cutT = static_cast<int>(cut.triangles.size());
    
    if (uvCoords.size() != static_cast<size_t>(cutV)) {
        std::cerr << "[IntegerGridMap] UV coordinates not computed yet." << std::endl;
        return false;
    }
    
    // Write original mesh vertices (v)
    for (int i = 0; i < origV; ++i) {
        out << "v " << orig.vertices[i][0] << " " << orig.vertices[i][1] << " 0\n";
    }
    
    // For UV coordinates, we need one vt per triangle corner because
    // vertices on the cut seam have different UVs on each side.
    // The cut mesh has the same triangles as the original, but with
    // duplicated vertices along cuts.
    
    // Write texture coordinates: one per cut-mesh vertex
    // These correspond to the cut mesh vertices, which map to original vertices
    for (int i = 0; i < cutV; ++i) {
        out << "vt " << uvCoords[i][0] << " " << uvCoords[i][1] << "\n";
    }
    
    // Write faces using original vertex indices but cut-mesh UV indices
    // The cut mesh triangles have the same connectivity structure, just with
    // potentially different vertex indices where cuts were made.
    // 
    // For each triangle t in the cut mesh:
    //   - cut.triangles[t] gives cut-mesh vertex indices (for UV lookup)
    //   - cutToOrig maps these to original vertex indices (for position lookup)
    
    for (int t = 0; t < cutT; ++t) {
        const auto &cutTri = cut.triangles[t];
        
        // Get original vertex indices
        const int ov0 = cutToOrig[cutTri[0]];
        const int ov1 = cutToOrig[cutTri[1]];
        const int ov2 = cutToOrig[cutTri[2]];
        
        // Cut-mesh vertex indices are used for UV (vt) lookup
        const int uv0 = cutTri[0];
        const int uv1 = cutTri[1];
        const int uv2 = cutTri[2];
        
        // OBJ format: f v/vt v/vt v/vt (1-based indices)
        out << "f "
            << ov0 + 1 << "/" << uv0 + 1 << " "
            << ov1 + 1 << "/" << uv1 + 1 << " "
            << ov2 + 1 << "/" << uv2 + 1 << "\n";
    }
    
    return true;
}

