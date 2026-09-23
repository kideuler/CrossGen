#include "Polysquare.hxx"

#include <algorithm>
#include <cmath>
#include <deque>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace {

// The closest rotation to a 2x2 matrix, i.e. the rotation part of its polar
// decomposition: the angle that maximizes tr(R^T J).
inline void polarRotation(const double J[2][2], double R[2][2]) {
    const double c = J[0][0] + J[1][1];
    const double s = J[1][0] - J[0][1];
    const double n = std::sqrt(c * c + s * s);
    if (n < 1e-14) { R[0][0] = 1.0; R[0][1] = 0.0; R[1][0] = 0.0; R[1][1] = 1.0; return; }
    const double cc = c / n, ss = s / n;
    R[0][0] = cc; R[0][1] = -ss;
    R[1][0] = ss; R[1][1] = cc;
}

inline double det2(const double J[2][2]) {
    return J[0][0] * J[1][1] - J[0][1] * J[1][0];
}

// Angle to the nearest axis, in (-45, 45] degrees.
inline double axisDeviation(double dx, double dy) {
    double a = std::atan2(dy, dx);
    a = std::fmod(a, M_PI_2);
    if (a > M_PI_4) a -= M_PI_2;
    if (a <= -M_PI_4) a += M_PI_2;
    return a;
}

} // namespace

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
Polysquare::Polysquare(const UMBER &frames, const HarmonicCut &cut) {
    orig = frames.getMeshPtr();
    if (!orig) throw std::runtime_error("Polysquare: the frame field carries no mesh");
    if (cut.getOriginalMeshPtr() != orig) {
        throw std::runtime_error("Polysquare: the HarmonicCut was not built on the frame field's mesh");
    }
    if (frames.getUField().size() != orig->triangles.size()) {
        throw std::runtime_error("Polysquare: the frame field has not been optimized yet");
    }

    harmonicCut = &cut;
    cutMesh = &cut.getCutMesh();
    frameU = frames.getUField();
    frameV = frames.getVField();

    // theta_i of Eq. (12) at each boundary vertex, in quarter turns: the corner
    // index the frame field already decided on.
    boundaryCorner.assign(orig->vertices.size(), 0);
    for (const auto &[v, k] : frames.boundarySingularities()) {
        if (v >= 0 && v < static_cast<int>(boundaryCorner.size())) boundaryCorner[v] = k;
    }

    // The interfaces the frame was aligned to, carried through so that the
    // measurement below asks of them exactly what it asks of dS. They are read
    // off the frame field and never passed separately: a deformation that
    // followed a curve the field was not aligned to would be asking the two
    // terms of Eq. (9) to pull against each other by construction.
    featureEdges = frames.getFeatureEdges();
    edgeIsFeature.assign(orig->edges.size(), 0);
    for (const int e : featureEdges)
        if (e >= 0 && e < static_cast<int>(edgeIsFeature.size())) edgeIsFeature[e] = 1;

    // theta_i of Eq. (12) where the curve is an interface: the frame field's
    // own corner index, keyed by the directed pair of edges it was measured
    // across. An interface has a layout on each side and the two need not turn
    // the same way, which is why the key is directed and not a vertex.
    featureChains = frames.interfaceChains();
    const long long nE = static_cast<long long>(orig->edges.size());
    for (const UMBER::FeatureCorner &fc : frames.interfaceCorners()) {
        interfaceCorner[static_cast<long long>(fc.edgeIn) * nE + fc.edgeOut] = fc.quarters;
    }
}

// ---------------------------------------------------------------------------
// Local topology helpers
// ---------------------------------------------------------------------------
int Polysquare::leftTriangleOf(int a, int b) const {
    // The triangle that lists a -> b in its own winding order.
    int found = -1;
    const auto &vt = orig->vertexTriangles;
    for (int k = vt.rowPtr[a]; k < vt.rowPtr[a + 1] && found < 0; ++k) {
        const int t = vt.colIdx[k];
        const Triangle &tri = orig->triangles[t];
        for (int i = 0; i < 3; ++i) {
            if (tri[i] == a && tri[(i + 1) % 3] == b) { found = t; break; }
        }
    }
    if (found < 0) return -1;

    // If that triangle is wound the other way the interior is on its right
    // instead, and the left one is the neighbour across the edge.
    const Triangle &tri = orig->triangles[found];
    const double s = cross2(orig->vertices[tri[1]] - orig->vertices[tri[0]],
                            orig->vertices[tri[2]] - orig->vertices[tri[0]]);
    if (s > 0.0) return found;
    return oppositeTriangle(found, a, b);
}

int Polysquare::oppositeTriangle(int f, int a, int b) const {
    for (int k = 0; k < 3; ++k) {
        const int e = orig->triangleEdges[f][k];
        if (e < 0) continue;
        if ((orig->edges[e][0] == a && orig->edges[e][1] == b) ||
            (orig->edges[e][0] == b && orig->edges[e][1] == a)) {
            const int o0 = orig->edgeTriangles[e][0];
            const int o1 = orig->edgeTriangles[e][1];
            return (o0 == f) ? o1 : o0;
        }
    }
    return -1;
}

int Polysquare::cornerOf(int f, int origVertex) const {
    const Triangle &t = orig->triangles[f];
    for (int i = 0; i < 3; ++i) if (t[i] == origVertex) return cutMesh->triangles[f][i];
    return -1;
}

// ---------------------------------------------------------------------------
// buildTopology()
//
// Hat gradients per triangle, the variable layout that eliminates Eq. (8), and
// the boundary edges of Eqs. (11) and (12) in loop order.
// ---------------------------------------------------------------------------
void Polysquare::buildTopology() {
    const int nT = static_cast<int>(orig->triangles.size());
    const int nVc = static_cast<int>(cutMesh->vertices.size());

    // --- grad lambda_i and areas -------------------------------------------
    // For a linear function on a triangle, grad lambda_i is the inward normal
    // of the opposite edge over twice the area.
    triGrad.assign(nT, TriGrad());
    totalArea = 0.0;
    for (int f = 0; f < nT; ++f) {
        const Triangle &t = orig->triangles[f];
        const Point &p0 = orig->vertices[t[0]];
        const Point &p1 = orig->vertices[t[1]];
        const Point &p2 = orig->vertices[t[2]];
        const double area2 = cross2(p1 - p0, p2 - p0); // signed
        TriGrad &g = triGrad[f];
        g.area = 0.5 * std::fabs(area2);
        if (std::fabs(area2) < 1e-18) {
            g.g[0] = g.g[1] = g.g[2] = Point{0.0, 0.0};
            continue;
        }
        const Point p[3] = {p0, p1, p2};
        for (int i = 0; i < 3; ++i) {
            const Point e = p[(i + 2) % 3] - p[(i + 1) % 3];
            g.g[i] = Point{-e[1] / area2, e[0] / area2}; // perp(e) / (2A)
        }
        totalArea += g.area;
    }
    if (totalArea <= 0.0) throw std::runtime_error("Polysquare: mesh has no area");
    for (int f = 0; f < nT; ++f) triGrad[f].weight = triGrad[f].area / totalArea;

    // --- Variable layout ----------------------------------------------------
    // One bank of each cut is free, the other follows it through Eq. (8).
    const auto &cuts = harmonicCut->getCuts();
    nCuts = static_cast<int>(cuts.size());
    varOf.assign(nVc, 0);
    depSrc.assign(nVc, -1);
    depCut.assign(nVc, -1);

    const auto &origToCut = harmonicCut->getOriginalToCutVertices();
    int unresolved = 0;

    for (int c = 0; c < nCuts; ++c) {
        const auto &path = cuts[c].path;
        for (size_t i = 0; i < path.size(); ++i) {
            const int p = path[i];
            const auto &copies = origToCut[p];
            if (copies.size() != 2) { ++unresolved; continue; }

            // A path edge at p, taken in path order, so the left bank is the
            // same side for every vertex of the cut.
            const int a = (i + 1 < path.size()) ? path[i] : path[i - 1];
            const int b = (i + 1 < path.size()) ? path[i + 1] : path[i];
            const int fL = leftTriangleOf(a, b);
            if (fL < 0) { ++unresolved; continue; }

            const int cvL = cornerOf(fL, p);
            if (cvL < 0) { ++unresolved; continue; }
            const int cvR = (copies[0] == cvL) ? copies[1] : copies[0];
            if (cvR == cvL) { ++unresolved; continue; }

            depSrc[cvR] = cvL;
            depCut[cvR] = c;
        }
    }
    if (unresolved > 0) {
        std::cerr << "Polysquare: " << unresolved << " cut vertex(es) could not be paired across "
                  << "their cut; those seams are left free.\n";
    }

    // One vertex is held at the origin to take out the global translation, and
    // it is taken from the interior: snapBoundary() later fixes a coordinate of
    // every boundary vertex, and a vertex already nailed to the origin could
    // not accept one.
    pinnedVertex = -1;
    for (int cv = 0; cv < nVc && pinnedVertex < 0; ++cv) {
        if (depSrc[cv] < 0 && !cutMesh->isBoundaryVertex[cv]) pinnedVertex = cv;
    }
    for (int cv = 0; cv < nVc && pinnedVertex < 0; ++cv) {
        if (depSrc[cv] < 0) pinnedVertex = cv;
    }

    nIndep = 0;
    for (int cv = 0; cv < nVc; ++cv) {
        if (depSrc[cv] >= 0) { varOf[cv] = -1; continue; }
        if (cv == pinnedVertex) { varOf[cv] = -2; continue; }
        varOf[cv] = nIndep++;
    }

    // --- Boundary edges of Eq. (11), in loop order -------------------------
    // The boundary of M_C minus the cuts is the boundary of M, so the loops of
    // the input mesh are what is walked; each edge is mapped to the copies its
    // one adjacent triangle uses.
    bEdges.clear();
    loopStart.clear();

    std::vector<char> isCutEndpoint(orig->vertices.size(), 0);
    for (const auto &cut : cuts) {
        if (!cut.path.empty()) {
            isCutEndpoint[cut.path.front()] = 1;
            isCutEndpoint[cut.path.back()] = 1;
        }
    }

    const auto &loops = harmonicCut->getBoundaryLoops();
    const int outerLoop = harmonicCut->getOuterLoop();

    for (size_t li = 0; li < loops.size(); ++li) {
        std::vector<int> ring = loops[li];
        if (ring.size() < 3) continue;

        // Walk with the interior on the left: counter-clockwise around the
        // outer loop, clockwise around a void. That is the convention the
        // corner index of Eq. (12) is measured in.
        double area2 = 0.0;
        for (size_t i = 0; i < ring.size(); ++i) {
            const Point &a = orig->vertices[ring[i]];
            const Point &b = orig->vertices[ring[(i + 1) % ring.size()]];
            area2 += a[0] * b[1] - b[0] * a[1];
        }
        const bool wantPositive = (static_cast<int>(li) == outerLoop);
        if ((area2 > 0.0) != wantPositive) std::reverse(ring.begin(), ring.end());

        loopStart.push_back(static_cast<int>(bEdges.size()));

        for (size_t i = 0; i < ring.size(); ++i) {
            const int a = ring[i];
            const int b = ring[(i + 1) % ring.size()];

            const int f = leftTriangleOf(a, b);
            if (f < 0) continue;

            BoundaryEdge be;
            be.ca = cornerOf(f, a);
            be.cb = cornerOf(f, b);
            if (be.ca < 0 || be.cb < 0) continue;
            be.length = normP(orig->vertices[b] - orig->vertices[a]);

            // Eq. (12) acts on this edge and the next, which meet at b.
            be.sharedOrigVertex = b;
            be.targetTurn = boundaryCorner[b] * M_PI_2;
            be.cornerPair = true;

            bEdges.push_back(be);
        }
    }
    loopStart.push_back(static_cast<int>(bEdges.size()));

    // Where a cut lands on the boundary the two edges meeting there use
    // different copies of that vertex, so their images are in frames a
    // transition apart. Comparing them as they stand is meaningless, and
    // dropping the pair is what leaves the boundary free to kink at exactly
    // the point a cut arrives -- worth two extra turns per cut, measurably.
    // The paper's answer is to carry the transition into the target angle,
    // theta_i -> theta_i -/+ theta_gamma, which is what this does; the sign is
    // whichever bank the second edge is on.
    for (size_t l = 0; l + 1 < loopStart.size(); ++l) {
        const int b0 = loopStart[l], b1 = loopStart[l + 1];
        const int n = b1 - b0;
        if (n < 2) continue;
        for (int i = 0; i < n; ++i) {
            BoundaryEdge &ei = bEdges[b0 + i];
            const BoundaryEdge &ej = bEdges[b0 + (i + 1) % n];
            if (ei.cb == ej.ca) continue; // same copy: no seam here

            // One of the two is the dependent copy, and its image carries
            // Pi_gamma^-1 relative to the other.
            if (varOf[ej.ca] == -1) {
                ei.seamTurn = -transitionK[depCut[ej.ca]] * M_PI_2;
            } else if (varOf[ei.cb] == -1) {
                ei.seamTurn = transitionK[depCut[ei.cb]] * M_PI_2;
            } else {
                ei.cornerPair = false; // unpaired seam: nothing to compare
                continue;
            }
            ei.targetTurn += ei.seamTurn;
        }
    }

    totalBoundaryLength = 0.0;
    for (const auto &be : bEdges) totalBoundaryLength += be.length;
    if (totalBoundaryLength <= 0.0) throw std::runtime_error("Polysquare: mesh has no boundary");

    buildFeatureChains();

    totalFeatureLength = totalBoundaryLength;
    for (const auto &fe : fEdges) totalFeatureLength += fe.length;
}

// ---------------------------------------------------------------------------
// buildFeatureChains()
//
// The material interfaces as edges of Eq. (11) and pairs of Eq. (12), exactly
// as dS is. An interface is a curve the output has to keep, so the image has
// to put it on an axis and turn it only where the frame says to -- and if the
// deformation is not asked for that, nothing downstream can follow it: the
// iso-lines of the image are what cut the model into blocks, and an interface
// that is not one of them is an interface every block across it straddles.
//
// Three differences from the boundary, all forced by what an interface is:
//
//  * It is interior, so it has one image rather than two banks. The copies are
//    taken from the triangle on the left of the walk, which is the bank a walk
//    with that material on its left would use, and is the side the frame
//    field measured its corner on.
//
//  * A chain is open unless the network closes on itself, so the pairs of
//    Eq. (12) do not wrap: the first and last edges of an open chain have a
//    neighbour on one side only, and where the chain ends -- a junction, a
//    landing on dS -- what the layout does is the junction's business and not
//    this chain's.
//
//  * theta_i comes from UMBER::interfaceCorners() rather than from
//    boundaryCorner, since the corner an interface turns is a property of the
//    side being walked.
// ---------------------------------------------------------------------------
void Polysquare::buildFeatureChains() {
    fEdges.clear();
    chainStart.clear();
    chainClosed.clear();
    if (featureChains.empty()) return;

    const long long nE = static_cast<long long>(orig->edges.size());
    for (const UMBER::FeatureChain &fc : featureChains) {
        if (fc.edges.size() < 2 || fc.verts.size() != fc.edges.size() + 1) continue;

        const int base = static_cast<int>(fEdges.size());
        bool ok = true;
        std::vector<BoundaryEdge> built;
        built.reserve(fc.edges.size());
        for (size_t i = 0; i < fc.edges.size(); ++i) {
            const int a = fc.verts[i], b = fc.verts[i + 1];
            const int f = leftTriangleOf(a, b);
            if (f < 0) { ok = false; break; }
            BoundaryEdge be;
            be.ca = cornerOf(f, a);
            be.cb = cornerOf(f, b);
            if (be.ca < 0 || be.cb < 0) { ok = false; break; }
            be.length = normP(orig->vertices[b] - orig->vertices[a]);
            be.sharedOrigVertex = b;
            // The pair (this edge, the next) meets at b; the target is what
            // the frame turns there, walked this way.
            const int eIn = fc.edges[i];
            const int eOut = fc.edges[(i + 1) % fc.edges.size()];
            const auto it = interfaceCorner.find(static_cast<long long>(eIn) * nE + eOut);
            be.targetTurn = (it != interfaceCorner.end()) ? it->second * M_PI_2 : 0.0;
            be.cornerPair = true;
            built.push_back(be);
        }
        if (!ok || built.size() < 2) continue;

        chainStart.push_back(base);
        chainClosed.push_back(fc.closed ? 1 : 0);
        for (auto &be : built) fEdges.push_back(be);
    }
    chainStart.push_back(static_cast<int>(fEdges.size()));

    // A cut crossing an interface puts the two edges either side of the
    // crossing on different copies of the shared vertex, exactly as it does on
    // dS, and the answer is the same: carry the transition into the target
    // angle, or drop the pair where neither copy is the dependent one and
    // there is nothing to compare.
    for (size_t c = 0; c + 1 < chainStart.size(); ++c) {
        const int f0 = chainStart[c], f1 = chainStart[c + 1];
        const int n = f1 - f0;
        if (n < 2) continue;
        const bool closed = chainClosed[c] != 0;
        for (int i = 0; i < n; ++i) {
            if (!closed && i == n - 1) { fEdges[f0 + i].cornerPair = false; continue; }
            BoundaryEdge &ei = fEdges[f0 + i];
            const BoundaryEdge &ej = fEdges[f0 + (i + 1) % n];
            if (ei.cb == ej.ca) continue;
            if (varOf[ej.ca] == -1) {
                ei.seamTurn = -transitionK[depCut[ej.ca]] * M_PI_2;
            } else if (varOf[ei.cb] == -1) {
                ei.seamTurn = transitionK[depCut[ei.cb]] * M_PI_2;
            } else {
                ei.cornerPair = false;
                continue;
            }
            ei.targetTurn += ei.seamTurn;
        }
    }
}

// ---------------------------------------------------------------------------
// extractTransitions()  --  Eq. (7)
//
// One rotation per cut, the multiple of 90 degrees that best carries the frame
// on the left bank onto the frame on the right. Trying all four and keeping
// the smallest error is what the paper does; there are only four.
// ---------------------------------------------------------------------------
void Polysquare::extractTransitions() {
    const auto &cuts = harmonicCut->getCuts();
    transitionK.assign(cuts.size(), 0);

    for (size_t c = 0; c < cuts.size(); ++c) {
        double best = std::numeric_limits<double>::max();
        int bestK = 0;

        for (int k = 0; k < 4; ++k) {
            double err = 0.0;
            const auto &path = cuts[c].path;
            for (size_t i = 0; i + 1 < path.size(); ++i) {
                const int a = path[i], b = path[i + 1];
                const int fL = leftTriangleOf(a, b);
                if (fL < 0) continue;
                const int fR = oppositeTriangle(fL, a, b);
                if (fR < 0) continue;

                const double len = normP(orig->vertices[b] - orig->vertices[a]);
                const Point rotated = rotateVector(frameU[fL], k);
                const Point d = rotated - frameU[fR];
                err += len * dotP(d, d);
            }
            if (err < best) { best = err; bestK = k; }
        }
        transitionK[c] = bestK;
    }
}

// ---------------------------------------------------------------------------
// expand() / fold()  --  the elimination of Eq. (8)
//
// phi on the right bank is Pi_gamma^-1 phi on the left plus a translation. The
// inverse is what Eq. (8) asks for: Eq. (7) aligns the *frames*, and the frame
// is the gradient of phi, so the map transitions the other way (the paper
// writes it as Pi_ab = Pi_gamma^-1).
// ---------------------------------------------------------------------------
void Polysquare::expand(const Eigen::VectorXd &state) {
    const int nVc = static_cast<int>(cutMesh->vertices.size());
    uv.assign(nVc, Point{0.0, 0.0});

    for (int cv = 0; cv < nVc; ++cv) {
        const int j = varOf[cv];
        if (j >= 0) uv[cv] = Point{state[2 * j], state[2 * j + 1]};
        else if (j == -2) uv[cv] = Point{0.0, 0.0};
    }
    for (int cv = 0; cv < nVc; ++cv) {
        if (varOf[cv] != -1) continue;
        const int c = depCut[cv];
        const int inv = (4 - transitionK[c]) % 4;
        const int base = 2 * nIndep + 2 * c;
        uv[cv] = rotateVector(uv[depSrc[cv]], inv) + Point{state[base], state[base + 1]};
    }
}

void Polysquare::fold(std::vector<Point> &gradUV, Eigen::VectorXd &grad) const {
    grad.setZero(2 * nIndep + 2 * nCuts);

    // The dependent copies first, so their share reaches the bank they follow.
    const int nVc = static_cast<int>(cutMesh->vertices.size());
    for (int cv = 0; cv < nVc; ++cv) {
        if (varOf[cv] != -1) continue;
        const int c = depCut[cv];
        const Point g = gradUV[cv];
        // d uv[cv] / d uv[src] = R(-k*90), whose transpose is R(+k*90).
        gradUV[depSrc[cv]] = gradUV[depSrc[cv]] + rotateVector(g, transitionK[c]);
        const int base = 2 * nIndep + 2 * c;
        grad[base] += g[0];
        grad[base + 1] += g[1];
    }
    for (int cv = 0; cv < nVc; ++cv) {
        const int j = varOf[cv];
        if (j < 0) continue;
        grad[2 * j] += gradUV[cv][0];
        grad[2 * j + 1] += gradUV[cv][1];
    }
}

// ---------------------------------------------------------------------------
// poissonInit()  --  Eq. (6)
//
// min_phi \int ||grad phi - (v v^perp)^T||_F^2, as a linear least squares in
// the reduced variables. The paper solves it without the seamless condition
// and takes the result as an initial value; imposing Eq. (8) here as well
// costs nothing (the elimination is already built) and starts the non-linear
// stage on a seamless map rather than one that has to be dragged onto the
// constraint.
// ---------------------------------------------------------------------------
void Polysquare::poissonInit() {
    const int nT = static_cast<int>(orig->triangles.size());
    const int nX = 2 * nIndep + 2 * nCuts;

    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(static_cast<size_t>(nT) * 36);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(4 * nT);

    // Row (f, r, c) is grad phi[r][c] - T[r][c], weighted by sqrt(area).
    for (int f = 0; f < nT; ++f) {
        const TriGrad &tg = triGrad[f];
        const double sw = std::sqrt(tg.area);
        if (sw <= 0.0) continue;

        const Point target[2] = {frameU[f], frameV[f]};

        for (int r = 0; r < 2; ++r) {
            for (int cc = 0; cc < 2; ++cc) {
                const int row = 4 * f + 2 * r + cc;
                rhs[row] = sw * target[r][cc];

                for (int i = 0; i < 3; ++i) {
                    const int cv = cutMesh->triangles[f][i];
                    const double g = sw * tg.g[i][cc];

                    // uv[cv] = M * x[src] + t, so row picks up M's r-th row.
                    int src = cv;
                    double m[2] = {0.0, 0.0};
                    m[r] = 1.0;
                    if (varOf[cv] == -1) {
                        src = depSrc[cv];
                        const int inv = (4 - transitionK[depCut[cv]]) % 4;
                        // r-th row of R(inv*90)
                        const Point e0 = rotateVector(Point{1.0, 0.0}, inv);
                        const Point e1 = rotateVector(Point{0.0, 1.0}, inv);
                        m[0] = e0[r];
                        m[1] = e1[r];
                        const int base = 2 * nIndep + 2 * depCut[cv];
                        trips.emplace_back(row, base + r, g);
                    }
                    const int j = varOf[src];
                    if (j < 0) continue; // pinned, or an unpaired seam vertex
                    if (m[0] != 0.0) trips.emplace_back(row, 2 * j + 0, g * m[0]);
                    if (m[1] != 0.0) trips.emplace_back(row, 2 * j + 1, g * m[1]);
                }
            }
        }
    }

    Eigen::SparseMatrix<double> C(4 * nT, nX);
    C.setFromTriplets(trips.begin(), trips.end());

    Eigen::SparseMatrix<double> A = C.transpose() * C;
    for (int i = 0; i < nX; ++i) A.coeffRef(i, i) += 1e-10; // the map is affine, not linear
    const Eigen::VectorXd b = C.transpose() * rhs;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Polysquare: the Eq. (6) system could not be factorized");
    }
    x = solver.solve(b);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Polysquare: the Eq. (6) system could not be solved");
    }
    expand(x);
}

// ---------------------------------------------------------------------------
// evaluate()  --  Eq. (9) and its gradient
// ---------------------------------------------------------------------------
double Polysquare::evaluate(const Eigen::VectorXd &state, Eigen::VectorXd &grad,
                            double *arapOut, double *l1Out, double *corOut) {
    expand(state);

    const int nT = static_cast<int>(orig->triangles.size());
    std::vector<Point> gradUV(cutMesh->vertices.size(), Point{0.0, 0.0});

    double eArap = 0.0, eBar = 0.0;

    // --- E_arap, Eq. (10), and the positive-Jacobian penalty ---------------
    for (int f = 0; f < nT; ++f) {
        const TriGrad &tg = triGrad[f];
        if (tg.area <= 0.0) continue;

        double J[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
        for (int i = 0; i < 3; ++i) {
            const Point &p = uv[cutMesh->triangles[f][i]];
            for (int r = 0; r < 2; ++r)
                for (int c = 0; c < 2; ++c) J[r][c] += p[r] * tg.g[i][c];
        }

        double R[2][2];
        polarRotation(J, R);

        double D[2][2];
        for (int r = 0; r < 2; ++r)
            for (int c = 0; c < 2; ++c) D[r][c] = J[r][c] - R[r][c];

        eArap += tg.weight * (D[0][0] * D[0][0] + D[0][1] * D[0][1] +
                              D[1][0] * D[1][0] + D[1][1] * D[1][1]);

        double dEdJ[2][2];
        for (int r = 0; r < 2; ++r)
            for (int c = 0; c < 2; ++c) dEdJ[r][c] = 2.0 * tg.weight * D[r][c];

        // det(grad phi) > 0, as a penalty on whatever falls below the floor.
        const double d = det2(J);
        if (d < detFloor) {
            const double s = (detFloor - d) / detFloor;
            eBar += barrierWeight * tg.weight * s * s;
            const double dEdDet = -2.0 * barrierWeight * tg.weight * s / detFloor;
            // d det / dJ is the cofactor matrix
            dEdJ[0][0] += dEdDet * J[1][1];
            dEdJ[0][1] += dEdDet * -J[1][0];
            dEdJ[1][0] += dEdDet * -J[0][1];
            dEdJ[1][1] += dEdDet * J[0][0];
        }

        for (int i = 0; i < 3; ++i) {
            const int cv = cutMesh->triangles[f][i];
            for (int r = 0; r < 2; ++r)
                gradUV[cv][r] += dEdJ[r][0] * tg.g[i][0] + dEdJ[r][1] * tg.g[i][1];
        }
    }

    // --- E_l1, Eq. (11): the boundary onto an axis -------------------------
    // |t| is replaced by sqrt(t^2 + eps^2) as in Sec. 4.2, and the direction
    // is normalized so the term measures the angle and not the length.
    double eL1 = 0.0;
    const double eps2 = l1Eps * l1Eps;
    // dS and the interfaces alike: the term asks a curve the output has to
    // keep to lie on an axis, and which kind of curve it is does not enter.
    auto axisTerm = [&](const BoundaryEdge &be, double scale) {
        if (scale <= 0.0) return;
        const Point d = uv[be.cb] - uv[be.ca];
        const double n2 = dotP(d, d);
        if (n2 < 1e-24) return;
        const double n = std::sqrt(n2);
        const double w = scale * be.length / totalFeatureLength;

        const double sx = std::sqrt(d[0] * d[0] + eps2);
        const double sy = std::sqrt(d[1] * d[1] + eps2);
        const double S = sx + sy;
        eL1 += w * S / n;

        const Point g{w * ((d[0] / sx) / n - S * d[0] / (n2 * n)),
                      w * ((d[1] / sy) / n - S * d[1] / (n2 * n))};
        gradUV[be.cb] = gradUV[be.cb] + g * l1Weight;
        gradUV[be.ca] = gradUV[be.ca] - g * l1Weight;
    };
    for (const auto &be : bEdges) axisTerm(be, 1.0);
    for (const auto &fe : fEdges) axisTerm(fe, interfaceWeight);

    // --- E_cor, Eq. (12): where the boundary is allowed to turn ------------
    double eCor = 0.0;
    if (corWeight > 0.0) {
        // One pair of consecutive edges, wherever the two curves meet. Written
        // once and applied to dS's loops and the interfaces' chains alike, so
        // that "do not turn here" means the same thing on both.
        auto cornerTerm = [&](const BoundaryEdge &ei, const BoundaryEdge &ej, double scale) {
                if (!ei.cornerPair || scale <= 0.0) return;

                const Point di = uv[ei.cb] - uv[ei.ca];
                const Point dj = uv[ej.cb] - uv[ej.ca];
                const double ni2 = dotP(di, di), nj2 = dotP(dj, dj);
                if (ni2 < 1e-24 || nj2 < 1e-24) return;
                const double ni = std::sqrt(ni2), nj = std::sqrt(nj2);
                const Point ui = di / ni, uj = dj / nj;

                const double ct = std::cos(ei.targetTurn), st = std::sin(ei.targetTurn);
                const Point rui{ct * ui[0] - st * ui[1], st * ui[0] + ct * ui[1]};

                const double wgeom = scale * (ei.length + ej.length) / totalFeatureLength;
                const double w = corWeight * wgeom;
                const Point diff = rui - uj;
                eCor += wgeom * dotP(diff, diff);

                // dE/d ui = -2 R^T uj, dE/d uj = -2 R ui, then through the
                // normalization: d u / d d = (I - u u^T)/n.
                const Point dEdui{-2.0 * (ct * uj[0] + st * uj[1]),
                                  -2.0 * (-st * uj[0] + ct * uj[1])};
                const Point dEduj{-2.0 * rui[0], -2.0 * rui[1]};

                const Point gi = (dEdui - ui * dotP(ui, dEdui)) / ni;
                const Point gj = (dEduj - uj * dotP(uj, dEduj)) / nj;

                gradUV[ei.cb] = gradUV[ei.cb] + gi * w;
                gradUV[ei.ca] = gradUV[ei.ca] - gi * w;
                gradUV[ej.cb] = gradUV[ej.cb] + gj * w;
                gradUV[ej.ca] = gradUV[ej.ca] - gj * w;
        };

        for (size_t l = 0; l + 1 < loopStart.size(); ++l) {
            const int b0 = loopStart[l], b1 = loopStart[l + 1];
            const int n = b1 - b0;
            if (n < 2) continue;
            for (int i = 0; i < n; ++i)
                cornerTerm(bEdges[b0 + i], bEdges[b0 + (i + 1) % n], 1.0);
        }
        // An open chain's last edge has no successor on it -- buildFeatureChains
        // clears its cornerPair -- so the wrap below only ever fires on a
        // closed one.
        for (size_t c = 0; c + 1 < chainStart.size(); ++c) {
            const int f0 = chainStart[c], f1 = chainStart[c + 1];
            const int n = f1 - f0;
            if (n < 2) continue;
            for (int i = 0; i < n; ++i)
                cornerTerm(fEdges[f0 + i], fEdges[f0 + (i + 1) % n], interfaceWeight);
        }
    }

    fold(gradUV, grad);

    // --- The boundary alignment of Eq. (23), as a pull ----------------------
    //
    // Each coordinate snapBoundary() gave a target is drawn towards it, and
    // nothing is held. The paper imposes this as a constraint (M psi = 0) and
    // so did this, by zeroing those components of the gradient, and that is
    // the single change that costs the most on a model with holes: a
    // constraint that cannot be met while the map stays injective is met
    // anyway, and the fold it needs to get there is permanent, because the
    // vertices that would have to move to undo it are the ones being held.
    // Measured over data/meshes/singlemat with the holes kept, the hard form
    // flips 670 triangles on geom014, 196 on geom032 and 40 on geom012, and
    // turns the image boundary 34, 32 and 108 times where the frame field
    // asked for 26, 24 and 46; with the same alignment asked for as a pull,
    // every one of those is zero flips and the turn count lands exactly on
    // what was asked, on all three.
    //
    // It is a pull rather than a demand for the same reason
    // setInterfaceWeight's term is: where the two can both be had, the weight
    // climbs high enough that the difference is far below anything the tracing
    // can resolve, and where they cannot, the one that matters wins. A
    // boundary a thousandth of an edge off its axis costs the motorcycle graph
    // nothing -- it reads angles and rounds them to quarter turns -- and a
    // folded triangle costs it the block.
    double eSnap = 0.0;
    if (snapHard) {
        // The constraint form, for the attempt described in optimize(): the
        // values are already on their targets and these components of the
        // gradient are zeroed, so every search direction L-BFGS builds out of
        // them is zero there too and they never move.
        for (int i = 0; i < static_cast<int>(hasTarget.size()) && i < grad.size(); ++i) {
            if (hasTarget[i]) grad[i] = 0.0;
        }
    } else if (snapWeight > 0.0) {
        for (int i = 0; i < static_cast<int>(hasTarget.size()) && i < grad.size(); ++i) {
            if (!hasTarget[i]) continue;
            const double d = state[i] - targetX[i];
            eSnap += targetW[i] * d * d;
            grad[i] += snapWeight * 2.0 * targetW[i] * d;
        }
    }

    if (arapOut) *arapOut = eArap;
    if (l1Out) *l1Out = eL1;
    if (corOut) *corOut = eCor;

    return eArap + l1Weight * eL1 + corWeight * eCor + eBar + snapWeight * eSnap;
}

// ---------------------------------------------------------------------------
// runLBFGS()  --  the same quasi-Newton loop UMBER uses on Eq. (1)
// ---------------------------------------------------------------------------
int Polysquare::runLBFGS(Eigen::VectorXd &state, int maxIter) {
    const int m = 10;
    const double c1 = 1e-4;
    const int maxBacktracks = 60;

    const int n = static_cast<int>(state.size());
    Eigen::VectorXd g(n), gNew(n), d(n), xNew(n);

    double f = evaluate(state, g);
    double gnorm = g.norm();

    std::deque<Eigen::VectorXd> S, Y;
    std::deque<double> rho;
    std::vector<double> alpha(m, 0.0);

    int iter = 0;
    for (; iter < maxIter; ++iter) {
        if (gnorm <= gradTolerance) break;

        Eigen::VectorXd q = g;
        const int k = static_cast<int>(S.size());
        for (int i = k - 1; i >= 0; --i) {
            alpha[i] = rho[i] * S[i].dot(q);
            q -= alpha[i] * Y[i];
        }
        if (k > 0) {
            const double yy = Y[k - 1].squaredNorm();
            if (yy > 1e-300) q *= S[k - 1].dot(Y[k - 1]) / yy;
        }
        for (int i = 0; i < k; ++i) {
            const double beta = rho[i] * Y[i].dot(q);
            q += S[i] * (alpha[i] - beta);
        }
        d = -q;

        double dg = d.dot(g);
        if (!(dg < 0.0)) { d = -g; dg = -g.squaredNorm(); }

        double step = (k == 0) ? std::min(1.0, 1.0 / std::max(gnorm, 1e-12)) : 1.0;
        bool progressed = false;
        double fNew = f;
        for (int bt = 0; bt < maxBacktracks; ++bt) {
            xNew = state + step * d;
            fNew = evaluate(xNew, gNew);
            if (std::isfinite(fNew) && fNew <= f + c1 * step * dg) { progressed = true; break; }
            step *= 0.5;
        }
        if (!progressed) break;

        Eigen::VectorXd s = xNew - state;
        Eigen::VectorXd y = gNew - g;
        const double sy = s.dot(y);
        if (sy > 1e-12) {
            if (static_cast<int>(S.size()) == m) { S.pop_front(); Y.pop_front(); rho.pop_front(); }
            S.push_back(std::move(s));
            Y.push_back(std::move(y));
            rho.push_back(1.0 / sy);
        }

        state = xNew;
        f = fNew;
        g = gNew;
        gnorm = g.norm();
    }

    return iter;
}

// ---------------------------------------------------------------------------
// optimize()  --  Eq. (9) over the w_l1 continuation
// ---------------------------------------------------------------------------
void Polysquare::optimize() {
    report_.iterations = 0;

    const std::vector<double> schedule = l1Schedule.empty() ? std::vector<double>{1.0} : l1Schedule;
    for (size_t s = 0; s < schedule.size(); ++s) {
        l1Weight = schedule[s];
        // The smoothing follows the weight down: a rounded corner is what lets
        // the early stages move, a sharp one is what finally snaps the edges
        // onto an axis.
        l1Eps = std::max(1e-3, 1e-2 * std::pow(0.5, static_cast<double>(s)));
        report_.iterations += runLBFGS(x, maxIterations);
    }

    // Anything still folded gets the penalty raised on it. The paper keeps the
    // map feasible throughout with a log barrier; this cannot, so it leans on
    // the penalty afterwards instead, and says so when that is not enough.
    for (int round = 0; round < 4; ++round) {
        expand(x);
        if (countFlips() == 0) break;
        barrierWeight *= 10.0;
        report_.iterations += runLBFGS(x, maxIterations);
    }

    // Eq. (23): the boundary onto its axes, and the interior with it.
    //
    // E_l1 and E_cor have done their work by now -- the corners are where they
    // are going to be -- so what is left is one coordinate per boundary
    // segment, and E_arap and the barrier against it. The weight climbs
    // through a continuation for the reason the l1 schedule does: the first
    // stage moves the boundary while the interior can still follow it, and
    // starting at the last one lands in a fold it cannot leave.
    if (snapBoundaryOn) {
        expand(x);
        snapBoundary();
        // And then the corners that face each other across an iso-line onto
        // one, which snapBoundary() cannot do: it works a segment at a time and
        // this is a statement about two of them.
        snapCorners();

        // E_l1 and E_cor stay switched on through this, where the paper's
        // Eq. (23) drops them. They had nothing to act on while the boundary
        // was being *held*, and against a pull they are allies of it rather
        // than rivals: E_l1 wants each edge on an axis, which is what the pull
        // wants, and E_cor wants no turn where the field marked none, which is
        // what the segmentation above is. Switching them off costs geom011
        // 100% coverage down to 39%.
        // --- Exactly if exactly can be had, and a pull if it cannot --------
        //
        // The constraint form first, which is the paper's: every boundary
        // segment put on its axis and the interior re-solved with it held
        // there. Where the image is already close to a polysquare that costs
        // nothing to impose and is worth a great deal, because "on the axis"
        // is not a matter of degree downstream -- a ray traced from one corner
        // towards another passes it on the wrong side if the two are a
        // rounding apart, carries on into the model, and costs a block at each
        // end. On data/meshes/multimat it is worth 14.5 points of coverage on
        // rocket and 28.7 on geom013, neither of which folds a triangle under
        // it.
        //
        // It is kept only if the map survives it. Where it is infeasible --
        // which is where Sec. 4.3 has left a diagonal run, and it is the same
        // models every time -- imposing it folds the interior permanently,
        // because the vertices that would have to move to undo the fold are
        // the ones being held: 670 triangles on singlemat/geom014, 196 on
        // geom032, and with them an image boundary that turns 34 and 32 times
        // where the frame field asked for 26 and 24. Those models fall back to
        // asking for the same alignment as a pull, through a continuation on
        // its weight, and finish with no folds at all and the turn count
        // exactly as asked.
        //
        // The order matters and cost a measurement to find. Running the pull
        // first and then trying to impose the result is the obvious
        // arrangement and it is worse: the pull settles the map somewhere the
        // constraint can no longer be reached from, so a model that would have
        // taken the exact form is refused it -- multimat/rocket falls back and
        // keeps 65% of the model where imposing it straight away keeps 79.5%.
        const Eigen::VectorXd relaxed = x;
        const double relaxedBarrier = barrierWeight;

        for (int i = 0; i < static_cast<int>(hasTarget.size()) && i < x.size(); ++i)
            if (hasTarget[i]) x[i] = targetX[i];
        snapHard = true;
        report_.iterations += runLBFGS(x, maxIterations);
        for (int round = 0; round < 4; ++round) {
            expand(x);
            if (countFlips() == 0) break;
            barrierWeight *= 10.0;
            report_.iterations += runLBFGS(x, maxIterations);
        }
        expand(x);
        measure();
        const int hardMiss = std::abs(report_.turns - report_.expectedTurns);
        const Eigen::VectorXd imposed = x;

        // And the same alignment as a pull, from the same starting point.
        snapHard = false;
        barrierWeight = relaxedBarrier;
        x = relaxed;
        for (const double w : snapSchedule) {
            snapWeight = w;
            report_.iterations += runLBFGS(x, maxIterations);
            for (int round = 0; round < 3; ++round) {
                expand(x);
                if (countFlips() == 0) break;
                barrierWeight *= 10.0;
                report_.iterations += runLBFGS(x, maxIterations);
            }
        }
        expand(x);
        measure();
        const int softMiss = std::abs(report_.turns - report_.expectedTurns);

        // --- Which of the two to keep ---------------------------------------
        //
        // On how far the image boundary's turn count lands from the number the
        // frame field asked for, and the imposed form on a tie.
        //
        // Folds are the obvious thing to compare and they are the wrong one.
        // A couple of folded triangles out of a few thousand cost nothing
        // measurable, so refusing the imposed form over them gives up exact
        // iso-lines to buy nothing; and the folds that do matter never come
        // alone. The run that folds 670 triangles on singlemat/geom014 is the
        // same run that turns the image boundary 34 times where the field
        // asked for 26, and it is the extra turns that cost the blocks,
        // because each one is a corner of a face that nothing traced a line
        // from. Measured over both corpora, comparing folds instead costs the
        // single-material mean 94.7% down to 93.1%.
        //
        // The turn count is the closest thing reachable from here to the
        // number the stage is really judged by -- how much of the model comes
        // out as four-sided blocks -- which Polysquare cannot see.
        // Folds do not enter it, and a veto on them was tried and refused.
        // Refusing the imposed form wherever it folds more than a thousandth
        // of the mesh takes the single-material mean from 94.7% to 93.6% and
        // -- the reason it was tried -- does not even buy the mesh quality it
        // was meant to: on geom012 the pull it falls back to hands the mesher
        // 27 folded elements where the imposed form hands it 9. The folds in
        // the output mesh come from the shape of the blocks, not from the
        // flips in the map that found them.
        report_.snapWasExact = (hardMiss <= softMiss);
        if (report_.snapWasExact) {
            snapHard = true;
            x = imposed;
            expand(x);
        }
    }

    expand(x);
    measure();
}

// ---------------------------------------------------------------------------
// snapBoundary()  --  Eq. (23), the boundary half of it
//
// Eq. (9) aligns the boundary by a penalty, so it lands near an axis rather
// than on one: across data/meshes most edges finish inside a tenth of a
// degree, but one model keeps a fourteen-edge run bulging out to 21 degrees
// where E_l1 and E_cor deadlock. That residual is not something the rest of
// the pipeline can absorb. Sec. 5 needs every boundary segment to have one
// normal and one projection h_i -- Eq. (13) asks for <c_i+ - c_i-, N_i> = 0
// outright -- and a curved run has neither, so its anti-degeneracy constraints
// cannot even be written down. Tracing iso-lines does not help: the motorcycle
// graph partitions the interior from the corners it is given and leaves the
// boundary exactly as it found it.
//
// The paper does not ask Eq. (9) to be exact either. Its Sec. 5.5 re-solve,
// Eq. (23), minimizes E_arap alone subject to det > 0 and M psi = 0, where the
// second constraint pins every boundary vertex of a segment to that segment's
// projection. So the alignment is imposed, not optimized for, and Eq. (9)'s
// job is only to settle where the corners go. This is that step without the
// simplification in front of it: the corners are the ones already found, each
// segment between them gets a single coordinate, and the interior is re-solved
// with those held fixed.
//
// Segments break at corners and also at the seams, since the two sides of a
// cut are different copies of a vertex and their images do not meet -- they
// are separate segments of the boundary of the polysquare, however continuous
// they look on the model.
// ---------------------------------------------------------------------------
void Polysquare::snapBoundary() {
    hasTarget.assign(2 * nIndep + 2 * nCuts, 0);
    targetX.assign(hasTarget.size(), 0.0);
    targetW.assign(hasTarget.size(), 0.0);
    report_.worstSegmentDeg = 0.0;
    report_.suspectSegments = 0;
    if (bEdges.empty()) return;

    // The deviations are counted in mean image boundary edges, so that w_s
    // means something comparable from model to model.
    const double scale = meanImageBoundaryEdge();
    const double norm = (scale > 0.0) ? 1.0 / (scale * scale) : 1.0;

    int conflicts = 0;
    // The segments are collected first and pinned afterwards, because what a
    // segment's coordinate should be is not a question about that segment
    // alone -- see alignRuns().
    std::vector<Run> runs;

    auto pin = [&](int cv, int axis, double value, double weight) {
        const int j = varOf[cv];
        if (j < 0) return; // the gauge vertex, which is interior by construction
        const int idx = 2 * j + axis;
        if (hasTarget[idx] && std::fabs(targetX[idx] - value) > 1e-9) {
            ++conflicts; // two segments want the same coordinate of one vertex
            return;
        }
        hasTarget[idx] = 1;
        targetX[idx] = value;
        targetW[idx] = weight;
    };

    for (size_t l = 0; l + 1 < loopStart.size(); ++l) {
        const int b0 = loopStart[l], b1 = loopStart[l + 1];
        const int n = b1 - b0;
        if (n < 1) continue;

        // Where this loop breaks into segments: at the seams, and wherever
        // the image actually turns.
        //
        // Reading the breaks off theta_i instead -- break where the frame
        // field marked a corner, and nowhere else -- is the tidier rule and it
        // does not survive contact with the constraint form below. A place the
        // image turns that the field did not mark is then inside a segment,
        // and putting that segment on one axis means straightening a real
        // corner out of the boundary, which is exactly the move that folds the
        // interior: multimat/rocket goes from an exact alignment with no folds
        // to two folded triangles and falls back to the pull, keeping 65% of
        // the model where the threshold keeps 79.5%. The threshold respects
        // the image it has; theta_i describes the image it was asked for.
        std::vector<char> breakAfter(n, 0);
        int breaks = 0;
        for (int i = 0; i < n; ++i) {
            const BoundaryEdge &ei = bEdges[b0 + i];
            const BoundaryEdge &ej = bEdges[b0 + (i + 1) % n];

            bool brk = (ei.cb != ej.ca); // a seam
            if (!brk) {
                const Point di = uv[ei.cb] - uv[ei.ca];
                const Point dj = uv[ej.cb] - uv[ej.ca];
                if (normP(di) > 1e-18 && normP(dj) > 1e-18) {
                    const double turn = std::fabs(
                        wrap_pi(computeAngle(dj) - computeAngle(di) - ei.seamTurn));
                    brk = turn > 30.0 * M_PI / 180.0;
                }
            }
            if (brk) { breakAfter[i] = 1; ++breaks; }
        }

        // Start on an edge that follows a break, so that a segment is never
        // split across the seam of the walk. A loop with no break at all -- the
        // two rings of an annulus, which a closed-form polysquare leaves
        // cornerless -- is one segment, and correctly so: its image is a single
        // straight line.
        int first = 0;
        if (breaks > 0) {
            for (int i = 0; i < n; ++i) {
                if (breakAfter[i]) { first = (i + 1) % n; break; }
            }
        }

        // The loop's segments, in order round it.
        std::vector<std::vector<int>> segs;
        {
            int taken = 0;
            while (taken < n) {
                std::vector<int> seg;
                while (taken < n) {
                    const int idx = (first + taken) % n;
                    seg.push_back(b0 + idx);
                    ++taken;
                    if (breaks > 0 && breakAfter[idx]) break;
                }
                if (!seg.empty()) segs.push_back(std::move(seg));
            }
        }

        // --- The axes, by alternating round the loop ------------------------
        //
        // Letting each segment pick its own axis from its own span is the
        // obvious rule and it is the one that has to go. Consecutive segments
        // are separated by exactly one corner, every corner turns exactly one
        // quarter (see UMBER::setMaxCornerQuarters), and a quarter turn swaps
        // the axis -- so which axis a segment is on is not a measurement, it is
        // a consequence of the one before it. Measuring it lets two neighbours
        // come out on the *same* axis, and then they ask the vertex between
        // them for two different values of one coordinate. `pin` sees the
        // conflict, keeps the first and drops the second, and the boundary
        // reverses direction there: an image that goes out along a line and
        // straight back down it, which is a needle of no width and costs its
        // face exactly as a needle from a two-quarter corner does. On
        // data/meshes/singlemat/geom012 that was sixteen vertices, and sixteen
        // of its twenty-five unasked-for turns measured 180.0 degrees.
        //
        // The parity always works out, so this can be imposed rather than
        // attempted: a loop sums to +4 or -4 quarters with every corner at
        // plus or minus one, so #(+1) - #(-1) = +-4 and the number of corners
        // is even. Alternating round an even cycle closes up.
        //
        // Which of the two phases is a genuine choice, and it is made for the
        // loop as a whole by how well each fits the image -- a segment costs
        // its span across the axis it would be put on -- rather than per
        // segment, which is what let them disagree in the first place.
        std::vector<int> segAxis(segs.size(), 1);
        {
            std::vector<double> spanX(segs.size(), 0.0), spanY(segs.size(), 0.0);
            for (size_t s = 0; s < segs.size(); ++s) {
                for (const int e : segs[s]) {
                    const Point d = uv[bEdges[e].cb] - uv[bEdges[e].ca];
                    spanX[s] += std::fabs(d[0]);
                    spanY[s] += std::fabs(d[1]);
                }
            }
            // Phase p puts segment s on axis (p + s) % 2. The cost of putting a
            // segment on axis 1 (a constant y, so it runs in x) is the span it
            // has in y, which is what would have to be flattened away.
            double cost[2] = {0.0, 0.0};
            for (int p = 0; p < 2; ++p) {
                for (size_t s = 0; s < segs.size(); ++s)
                    cost[p] += ((p + s) % 2 == 1) ? spanY[s] : spanX[s];
            }
            // An odd number of segments cannot alternate round the loop at all.
            // It means a corner went missing upstream, and there is nothing
            // structural left to impose, so each segment keeps its own reading.
            const bool alternates = (segs.size() % 2 == 0) && segs.size() >= 2;
            for (size_t s = 0; s < segs.size(); ++s) {
                if (alternates) {
                    const int p = (cost[0] <= cost[1]) ? 0 : 1;
                    segAxis[s] = static_cast<int>((p + s) % 2);
                } else {
                    segAxis[s] = (spanX[s] >= spanY[s]) ? 1 : 0;
                }
            }
        }

        for (size_t si = 0; si < segs.size(); ++si) {
            const std::vector<int> &seg = segs[si];
            const int axis = segAxis[si];

            {   // How far the segment as a whole is from the axis it is about
                // to be snapped onto. A few degrees is discretization; a large
                // value means Eq. (9) left a genuinely diagonal run here and
                // the snap is about to straighten something that should have
                // been a staircase or a corner somewhere else.
                //
                // Against the axis it is actually given, not against whichever
                // is nearer: now that the axes alternate round the loop a
                // segment can be handed the further one, and that is precisely
                // the case worth hearing about.
                Point chord{0.0, 0.0};
                for (int e : seg) chord = chord + (uv[bEdges[e].cb] - uv[bEdges[e].ca]);
                if (normP(chord) > 1e-12) {
                    const double across = std::fabs(chord[axis]);
                    const double along = std::fabs(chord[1 - axis]);
                    const double dev = std::atan2(across, along) * 180.0 / M_PI;
                    report_.worstSegmentDeg = std::max(report_.worstSegmentDeg, dev);
                    if (dev > 10.0) ++report_.suspectSegments;
                }
            }

            // h_i: the least-squares constant over the segment, weighted by
            // edge length so a fan of short edges cannot outvote a long one.
            double num = 0.0, den = 0.0;
            for (int e : seg) {
                const BoundaryEdge &be = bEdges[e];
                const double w = be.length;
                num += w * 0.5 * (uv[be.ca][axis] + uv[be.cb][axis]);
                den += w;
            }
            const double h = (den > 0.0) ? num / den : 0.0;

            // A segment that reaches the far bank of a cut is left alone,
            // which is the paper's own rule (Sec. 5.4: "we do not change the
            // projection values of a cut and the segments adjacent to the
            // cut"). The two banks are one rigid transition apart, so a cut's
            // two ends cannot both be moved to wherever their segments want:
            // that is two positions asked of two degrees of freedom, and it is
            // over-determined as soon as both ends are pinned on the same
            // axis. Snapping them anyway does not fail quietly -- it leaves
            // one vertex stranded 37 degrees off its segment on geom011, with
            // five turns that the frame field never asked for, and raising the
            // weight from 1e3 to 1e6 does not move it, because the constraint
            // is infeasible rather than weak. Left free, and with E_l1 and
            // E_cor still acting on them through the re-solve, the same
            // segments finish 0.65 degrees off.
            bool touchesFarBank = false;
            for (int e : seg) if (varOf[bEdges[e].ca] == -1) touchesFarBank = true;
            if (varOf[bEdges[seg.back()].cb] == -1) touchesFarBank = true;
            if (touchesFarBank) continue;

            Run run;
            run.axis = axis;
            run.h = h;
            run.weight = den;
            for (int e : seg) run.verts.push_back(bEdges[e].ca);
            run.verts.push_back(bEdges[seg.back()].cb);
            runs.push_back(std::move(run));
        }
    }

    alignRuns(runs);
    for (const Run &r : runs) {
        // A run's image length, spread over the vertices on it, so that a long
        // run outweighs a short one vertex for vertex too. It is a length and
        // not a share of the boundary: dividing through by the total makes w_s
        // scale-free, which reads better and measures worse -- the pull on a
        // model with a long boundary is then weak enough to leave it halfway,
        // and halfway is the one place it must not stop. Over
        // data/meshes/singlemat that costs geom018 and geom032 a fifth of the
        // model each and the corpus mean 93.3% down to 92.3%.
        const double w = (r.verts.empty()) ? 0.0
                                           : norm * r.weight / static_cast<double>(r.verts.size());
        for (const int cv : r.verts) pin(cv, r.axis, r.h, w);
    }

    if (conflicts > 0) {
        std::cerr << "Polysquare: " << conflicts << " boundary vertex(es) wanted two different "
                  << "values for the same coordinate; the first was kept.\n";
    }
}

// ---------------------------------------------------------------------------
// meanImageBoundaryEdge()
// ---------------------------------------------------------------------------
double Polysquare::meanImageBoundaryEdge() const {
    if (bEdges.empty()) return 0.0;
    double total = 0.0;
    for (const auto &be : bEdges) total += normP(uv[be.cb] - uv[be.ca]);
    return total / static_cast<double>(bEdges.size());
}

// ---------------------------------------------------------------------------
// alignRuns()  --  runs that lie on one iso-line given one coordinate
//
// snapBoundary() puts each straight run of the boundary on a coordinate of its
// own, computed from that run and nothing else. Two runs that are meant to be
// the same iso-line therefore come out at two numbers, and a motorcycle sent
// down one of them misses whatever sits on the other.
//
// data/meshes/singlemat/geom031 is what that costs. It is a comb: thirty reflex corners
// whose teeth all sit at one u in the shape they are meant to be. They come out
// spread over nine different u, 2.6e-3 apart end to end -- a tenth of a boundary
// edge -- so every iso-line launched from a tooth misses every other tooth, and
// the sixty rays cross each other 4336 times instead of ending on one another.
// The result is 4381 blocks for a shape that wants a few dozen.
//
// Doing this a run at a time rather than a corner at a time is the whole point.
// A corner has its coordinate from the run it lies on, so moving the corner
// alone would tilt that run off its axis -- and axis alignment is the property
// Sec. 5 cannot do without. A run moves rigidly: every vertex on it shifts by
// the same amount, it stays exactly as straight as it was, and only its
// position changes. So this buys the corners without spending the alignment.
//
// What it must not do is merge two runs that really are at different
// coordinates, and the guard is that nothing moves further than the tolerance:
// runs join a cluster only while consecutive values are within it and the whole
// cluster spans no more than twice it.
// ---------------------------------------------------------------------------
void Polysquare::alignRuns(std::vector<Run> &runs) {
    if (cornerSnapTol <= 0.0 || runs.empty()) return;
    const double meanEdge = meanImageBoundaryEdge();
    if (!(meanEdge > 0.0)) return;
    const double tol = cornerSnapTol * meanEdge;

    for (int axis = 0; axis < 2; ++axis) {
        std::vector<int> idx;
        for (int i = 0; i < static_cast<int>(runs.size()); ++i)
            if (runs[i].axis == axis) idx.push_back(i);
        if (idx.size() < 2) continue;

        std::sort(idx.begin(), idx.end(),
                  [&](int a, int b) { return runs[a].h < runs[b].h; });

        // Single linkage: runs join while each is within the tolerance of the
        // one before it. It has to be the gap between neighbours rather than
        // the width of the whole cluster, because the thing this is for is a
        // chain -- geom031's nine values are 1% of a boundary edge apart in
        // sequence and 10% apart end to end, and a cluster capped at twice the
        // tolerance splits it down the middle and aligns neither half to the
        // other. The cap is kept, ten times wider, only so that a genuine
        // staircase of steps each under the tolerance cannot be linked into one
        // line: at the default that is half a boundary edge, against the 0.08
        // that is the widest any model here actually moves a run.
        const double maxSpan = 10.0 * tol;
        size_t i = 0;
        while (i < idx.size()) {
            size_t j = i + 1;
            while (j < idx.size() && runs[idx[j]].h - runs[idx[j - 1]].h <= tol &&
                   runs[idx[j]].h - runs[idx[i]].h <= maxSpan) ++j;
            if (j - i >= 2) {
                double num = 0.0, den = 0.0;
                for (size_t k = i; k < j; ++k) {
                    num += runs[idx[k]].weight * runs[idx[k]].h;
                    den += runs[idx[k]].weight;
                }
                const double target = (den > 0.0) ? num / den : runs[idx[i]].h;
                for (size_t k = i; k < j; ++k) {
                    const double move = std::fabs(runs[idx[k]].h - target);
                    if (move <= 1e-12) continue;
                    runs[idx[k]].h = target;
                    ++report_.runsAligned;
                    report_.worstRunAlign = std::max(report_.worstRunAlign, move / meanEdge);
                }
            }
            i = j;
        }
    }
}

// ---------------------------------------------------------------------------
// snapCorners()  --  corners that share an iso-line, put on one
//
// The motorcycle graph traces an axis-aligned line inward from every reflex
// corner, and in a polysquare that line usually runs into another corner: two
// corners of one notch face each other across it, and the u they sit at is the
// same u. That is what makes the line stop where the structure says it should.
//
// Eq. (9) does not know it. E_l1 asks each boundary *edge* to lie on an axis
// and E_cor asks consecutive edges on a run to be parallel; snapBoundary() then
// puts each *segment* on a single coordinate of its own. Every one of those is
// a statement about one segment, and none of them says that two segments which
// ought to be collinear carry the same number. So two corners that face each
// other come out a little apart, and "a little" is enough: the ray leaves one
// of them, passes the other on the wrong side, and carries on into the model.
//
// data/meshes/singlemat/geom021 is the case that shows what it costs. Its two reflex
// corners sit at u = 0.125905 and u = 0.126426 -- 5.2e-4 apart, 1.4% of a mesh
// edge -- because the segment they are on is the one the cut lands on, so
// snapBoundary() leaves one half of it alone (see the note on the far bank
// there) and the two halves keep different coordinates. The ray from the first
// passes the second at 6.4e-4 and misses. On that model it then has nowhere
// else to go and circles 71 times before the step cap stops it, but the miss
// itself is the ordinary failure and the circling is a second defect on top.
//
// So: cluster the corners by each coordinate separately and give a cluster one
// value. A corner whose coordinate snapBoundary() already fixed is the one that
// says what the value is -- moving it would pull it off a segment that is
// genuinely aligned -- and only the corners it left free are moved. That bounds
// what this can do: no corner moves further than the tolerance, and no aligned
// segment is touched at all.
// ---------------------------------------------------------------------------
void Polysquare::snapCorners() {
    report_.cornersSnapped = 0;
    report_.worstCornerSnap = 0.0;
    report_.cornerSnapConflicts = 0;
    if (cornerSnapTol <= 0.0 || bEdges.empty()) return;

    // The scale to measure "the same iso-line" in: a boundary edge of the
    // image. Two corners closer than a fraction of one are not two places as
    // far as the mesh is concerned.
    const double meanEdge = meanImageBoundaryEdge();
    if (!(meanEdge > 0.0)) return;
    const double tol = cornerSnapTol * meanEdge;

    const std::vector<int> &toOrig = harmonicCut->getCutVertexToOriginal();

    // A corner's coordinate, and whether anything may move it. Fixed means
    // either snapBoundary() pinned it to a segment or the vertex is the far
    // bank of a cut and follows its transition; both are already answerable to
    // something else.
    struct Corner {
        double value = 0.0;
        int cv = -1;
        bool fixed = false;
    };

    for (int axis = 0; axis < 2; ++axis) {
        std::vector<Corner> cs;
        for (int cv = 0; cv < static_cast<int>(uv.size()); ++cv) {
            if (cv >= static_cast<int>(toOrig.size())) continue;
            const int ov = toOrig[cv];
            if (ov < 0 || ov >= static_cast<int>(boundaryCorner.size())) continue;
            if (boundaryCorner[ov] == 0) continue;
            const int j = varOf[cv];
            Corner c;
            c.cv = cv;
            c.fixed = (j < 0) || hasTarget[2 * j + axis];
            // Where the corner is going, not where it currently sits: a corner
            // snapBoundary() has already aimed at a segment will end up on
            // that segment, and clustering on anything else would pair it with
            // a neighbour it is about to leave.
            c.value = (j >= 0 && hasTarget[2 * j + axis]) ? targetX[2 * j + axis]
                                                          : uv[cv][axis];
            cs.push_back(c);
        }
        if (cs.size() < 2) continue;

        std::sort(cs.begin(), cs.end(),
                  [](const Corner &a, const Corner &b) { return a.value < b.value; });

        // Single linkage at the tolerance, with the span of a cluster capped so
        // that a chain of corners a tolerance apart cannot drag one of them
        // across the model. Nothing moves further than the tolerance.
        size_t i = 0;
        while (i < cs.size()) {
            size_t j = i + 1;
            while (j < cs.size() && cs[j].value - cs[j - 1].value <= tol &&
                   cs[j].value - cs[i].value <= 2.0 * tol) ++j;

            // What the cluster's coordinate is: whatever the corners that
            // cannot move already say, and the average otherwise.
            double target = 0.0;
            int fixedCount = 0;
            double fixedLo = 0.0, fixedHi = 0.0;
            for (size_t k = i; k < j; ++k) {
                if (!cs[k].fixed) continue;
                if (fixedCount == 0) { fixedLo = fixedHi = cs[k].value; }
                fixedLo = std::min(fixedLo, cs[k].value);
                fixedHi = std::max(fixedHi, cs[k].value);
                ++fixedCount;
            }
            if (fixedCount > 0) {
                // Two corners that cannot move and do not agree: there is no
                // value that puts both on one line, so the cluster is left as
                // it is rather than moved onto neither of them.
                if (fixedHi - fixedLo > 1e-12) {
                    ++report_.cornerSnapConflicts;
                    i = j;
                    continue;
                }
                target = fixedLo;
            } else {
                for (size_t k = i; k < j; ++k) target += cs[k].value;
                target /= static_cast<double>(j - i);
            }

            for (size_t k = i; k < j; ++k) {
                if (cs[k].fixed) continue;
                const double move = std::fabs(cs[k].value - target);
                if (move <= 1e-12) continue;
                const int idx = 2 * varOf[cs[k].cv] + axis;
                hasTarget[idx] = 1;
                targetX[idx] = target;
                // A corner is weighted as a whole where a run vertex carries
                // only its share of a length, which puts it about two orders
                // higher. That is what the structure is made of: the
                // motorcycle graph traces from corner to corner, and a ray
                // that passes the corner it was aimed at on the wrong side
                // carries on into the model and costs a block, where a run
                // bulging a degree off its axis between two corners that are
                // both in place costs nothing at all.
                targetW[idx] = 1.0 / (meanEdge * meanEdge);
                ++report_.cornersSnapped;
                report_.worstCornerSnap = std::max(report_.worstCornerSnap, move / meanEdge);
            }
            i = j;
        }
    }
}

int Polysquare::countFlips() const {
    const int nT = static_cast<int>(orig->triangles.size());
    int flips = 0;
    for (int f = 0; f < nT; ++f) {
        if (triGrad[f].area <= 0.0) continue;
        double J[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
        for (int i = 0; i < 3; ++i) {
            const Point &p = uv[cutMesh->triangles[f][i]];
            for (int r = 0; r < 2; ++r)
                for (int c = 0; c < 2; ++c) J[r][c] += p[r] * triGrad[f].g[i][c];
        }
        if (det2(J) <= 0.0) ++flips;
    }
    return flips;
}

// ---------------------------------------------------------------------------
// solve()
// ---------------------------------------------------------------------------
void Polysquare::solve() {
    // The transitions come first: buildTopology() needs them, both to define
    // the dependent bank of each cut and to carry theta_gamma into the corner
    // targets where a cut meets the boundary.
    extractTransitions();
    buildTopology();
    poissonInit();
    optimize();
}

// ---------------------------------------------------------------------------
// measure()
// ---------------------------------------------------------------------------
void Polysquare::measure() {
    const int nT = static_cast<int>(orig->triangles.size());

    report_.flips = 0;
    report_.minScaledJacobian = std::numeric_limits<double>::max();
    double sjSum = 0.0;
    int sjCount = 0;

    for (int f = 0; f < nT; ++f) {
        if (triGrad[f].area <= 0.0) continue;
        double J[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
        for (int i = 0; i < 3; ++i) {
            const Point &p = uv[cutMesh->triangles[f][i]];
            for (int r = 0; r < 2; ++r)
                for (int c = 0; c < 2; ++c) J[r][c] += p[r] * triGrad[f].g[i][c];
        }
        const double d = det2(J);
        if (d <= 0.0) ++report_.flips;

        const double c0 = std::sqrt(J[0][0] * J[0][0] + J[1][0] * J[1][0]);
        const double c1 = std::sqrt(J[0][1] * J[0][1] + J[1][1] * J[1][1]);
        const double sj = (c0 * c1 > 1e-18) ? d / (c0 * c1) : 0.0;
        report_.minScaledJacobian = std::min(report_.minScaledJacobian, sj);
        sjSum += sj;
        ++sjCount;
    }
    report_.avgScaledJacobian = sjCount ? sjSum / sjCount : 0.0;
    if (sjCount == 0) report_.minScaledJacobian = 0.0;

    // Turns of the image boundary. A quarter turn is 90 degrees and the noise
    // is a fraction of one, so anything past 30 counts.
    report_.turns = 0;
    report_.expectedTurns = 0;
    report_.shortestRun = std::numeric_limits<int>::max();
    for (size_t l = 0; l + 1 < loopStart.size(); ++l) {
        const int b0 = loopStart[l], b1 = loopStart[l + 1];
        const int n = b1 - b0;
        if (n < 2) continue;
        int run = 0;
        for (int i = 0; i < n; ++i) {
            const BoundaryEdge &ei = bEdges[b0 + i];
            const BoundaryEdge &ej = bEdges[b0 + (i + 1) % n];
            if (boundaryCorner[ei.sharedOrigVertex] != 0) ++report_.expectedTurns;

            const Point di = uv[ei.cb] - uv[ei.ca];
            const Point dj = uv[ej.cb] - uv[ej.ca];
            if (normP(di) < 1e-18 || normP(dj) < 1e-18) continue;
            const double turn =
                std::fabs(wrap_pi(computeAngle(dj) - computeAngle(di) - ei.seamTurn));
            if (turn > 30.0 * M_PI / 180.0) {
                ++report_.turns;
                report_.shortestRun = std::min(report_.shortestRun, run);
                run = 0;
            } else {
                ++run;
            }
        }
    }
    if (report_.turns == 0) report_.shortestRun = static_cast<int>(bEdges.size());

    // Boundary alignment and length.
    double devSum = 0.0, devMax = 0.0, imageLen = 0.0;
    for (const auto &be : bEdges) {
        const Point d = uv[be.cb] - uv[be.ca];
        const double n = normP(d);
        imageLen += n;
        if (n < 1e-18) continue;
        const double dev = std::fabs(axisDeviation(d[0], d[1])) * 180.0 / M_PI;
        devSum += dev;
        devMax = std::max(devMax, dev);
    }
    report_.meanAlignDeg = bEdges.empty() ? 0.0 : devSum / bEdges.size();
    report_.maxAlignDeg = devMax;
    report_.lengthRatio = imageLen / totalBoundaryLength;

    // The same reading on the interfaces. An interface edge is interior, so it
    // has one image and not two banks -- unless a cut happens to run along it,
    // and then the two triangles disagree and the edge is read from the one on
    // the left, which is the bank a walk with the material on its left would
    // take.
    report_.interfaceEdges = 0;
    report_.meanInterfaceAlignDeg = 0.0;
    report_.maxInterfaceAlignDeg = 0.0;
    {
        double sum = 0.0, worst = 0.0;
        int n = 0;
        for (const int e : featureEdges) {
            if (e < 0 || e >= static_cast<int>(orig->edges.size())) continue;
            const int a = orig->edges[e][0], b = orig->edges[e][1];
            const int f = leftTriangleOf(a, b);
            if (f < 0) continue;
            const int ca = cornerOf(f, a), cb = cornerOf(f, b);
            if (ca < 0 || cb < 0) continue;
            const Point d = uv[cb] - uv[ca];
            if (normP(d) < 1e-18) continue;
            const double dev = std::fabs(axisDeviation(d[0], d[1])) * 180.0 / M_PI;
            sum += dev;
            worst = std::max(worst, dev);
            ++n;
        }
        report_.interfaceEdges = n;
        report_.meanInterfaceAlignDeg = n ? sum / n : 0.0;
        report_.maxInterfaceAlignDeg = worst;
    }

    // Eq. (8) holds by construction, so there is nothing to measure there.
    // What is worth measuring is the assumption underneath it: that a single
    // k*90-degree rotation really does line the frames up across each cut. A
    // large residual here means the cut runs through a part of the field that
    // is not smooth, and the parameterization inherits that.
    report_.transitionDeg = 0.0;
    for (size_t c = 0; c < harmonicCut->getCuts().size(); ++c) {
        const auto &path = harmonicCut->getCuts()[c].path;
        for (size_t i = 0; i + 1 < path.size(); ++i) {
            const int fL = leftTriangleOf(path[i], path[i + 1]);
            if (fL < 0) continue;
            const int fR = oppositeTriangle(fL, path[i], path[i + 1]);
            if (fR < 0) continue;
            const Point rotated = rotateVector(frameU[fL], transitionK[c]);
            const double deg = std::fabs(wrap_pi(computeAngle(rotated) - computeAngle(frameU[fR]))) *
                               180.0 / M_PI;
            report_.transitionDeg = std::max(report_.transitionDeg, deg);
        }
    }

    Eigen::VectorXd g(x.size());
    evaluate(x, g, &report_.arap, &report_.l1, &report_.cor);
}

// ---------------------------------------------------------------------------
// VTK output
// ---------------------------------------------------------------------------
namespace {

void writeHeader(std::ofstream &out, size_t nPoints, size_t nCells) {
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << nCells << "\">\n";
}

void writeCells(std::ofstream &out, const std::vector<Triangle> &tris) {
    out << "      <Cells>\n";
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto &t : tris) out << "          " << t[0] << " " << t[1] << " " << t[2] << "\n";
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (size_t i = 1; i <= tris.size(); ++i) out << "          " << i * 3 << "\n";
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < tris.size(); ++i) out << "          5\n";
    out << "        </DataArray>\n";
    out << "      </Cells>\n";
}

void writeFooter(std::ofstream &out) {
    out << "    </Piece>\n";
    out << "  </UnstructuredGrid>\n";
    out << "</VTKFile>\n";
}

} // namespace

bool Polysquare::writeVTU(const std::string &filename) const {
    if (uv.empty()) return false;
    std::ofstream out(filename);
    if (!out) return false;

    const int nT = static_cast<int>(orig->triangles.size());
    writeHeader(out, uv.size(), cutMesh->triangles.size());

    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto &p : uv) out << "          " << p[0] << " " << p[1] << " 0.0\n";
    out << "        </DataArray>\n";
    out << "      </Points>\n";

    out << "      <PointData Vectors=\"source\">\n";
    out << "        <DataArray type=\"Float64\" Name=\"source\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (size_t i = 0; i < uv.size(); ++i) {
        const Point &s = cutMesh->vertices[i];
        out << "          " << s[0] << " " << s[1] << " 0.0\n";
    }
    out << "        </DataArray>\n";
    out << "      </PointData>\n";

    out << "      <CellData Scalars=\"scaledJacobian\">\n";
    out << "        <DataArray type=\"Float64\" Name=\"scaledJacobian\" format=\"ascii\">\n";
    for (int f = 0; f < nT; ++f) {
        double J[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
        for (int i = 0; i < 3; ++i) {
            const Point &p = uv[cutMesh->triangles[f][i]];
            for (int r = 0; r < 2; ++r)
                for (int c = 0; c < 2; ++c) J[r][c] += p[r] * triGrad[f].g[i][c];
        }
        const double d = det2(J);
        const double c0 = std::sqrt(J[0][0] * J[0][0] + J[1][0] * J[1][0]);
        const double c1 = std::sqrt(J[0][1] * J[0][1] + J[1][1] * J[1][1]);
        out << "          " << ((c0 * c1 > 1e-18) ? d / (c0 * c1) : 0.0) << "\n";
    }
    out << "        </DataArray>\n";
    out << "      </CellData>\n";

    writeCells(out, cutMesh->triangles);
    writeFooter(out);
    return true;
}

bool Polysquare::writeSourceVTU(const std::string &filename) const {
    if (uv.empty()) return false;
    std::ofstream out(filename);
    if (!out) return false;

    writeHeader(out, cutMesh->vertices.size(), cutMesh->triangles.size());

    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto &p : cutMesh->vertices) out << "          " << p[0] << " " << p[1] << " 0.0\n";
    out << "        </DataArray>\n";
    out << "      </Points>\n";

    out << "      <PointData Vectors=\"uv\">\n";
    out << "        <DataArray type=\"Float64\" Name=\"uv\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto &p : uv) out << "          " << p[0] << " " << p[1] << " 0.0\n";
    out << "        </DataArray>\n";
    out << "      </PointData>\n";

    writeCells(out, cutMesh->triangles);
    writeFooter(out);
    return true;
}
