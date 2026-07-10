#include "UVIsoTrace.hxx"

UVIsoTrace::UVIsoTrace(std::shared_ptr<UVGParam> uvParam, int nU, int nV)
    : uvParam_(uvParam)
    , uIntervalTree_(*uvParam, IntervalTree::Axis::U)
    , vIntervalTree_(*uvParam, IntervalTree::Axis::V)
{
    nU_ = nU;
    nV_ = nV;
    const Eigen::VectorXd &u = uvParam_->getU();
    const Eigen::VectorXd &v = uvParam_->getV();
    if (u.size() == 0 || v.size() == 0) {
        uvMin_ = {0.0, 0.0};
        uvMax_ = {1.0, 1.0};
    } else {
        uvMin_ = {u.minCoeff(), v.minCoeff()};
        uvMax_ = {u.maxCoeff(), v.maxCoeff()};
    }

    // Compute step sizes for uniform sampling
    deltaU_ = (uvMax_[0] - uvMin_[0]) / static_cast<double>(nU_-1);
    deltaV_ = (uvMax_[1] - uvMin_[1]) / static_cast<double>(nV_-1);

    // Compute which triangles are valid for tracing (not near singularities or flipped)
    const Mesh& mesh = uvParam_->getCutMesh().getCutMesh();
    int nT = static_cast<int>(mesh.triangles.size());
    isValidTriangle_.resize(nT, true);
    for (int t = 0; t < nT; ++t) {
        const Triangle& tri = mesh.triangles[t];
        int i = tri[0], j = tri[1], k = tri[2];

        // Check if triangle is flipped in UV space
        double signedArea2 = (u(j) - u(i)) * (v(k) - v(i)) - (u(k) - u(i)) * (v(j) - v(i));
        if (signedArea2 <= 0.0) {
            isValidTriangle_[t] = false;
            continue;
        }

        // check if any vertex within this triangle is singular
        const std::vector<bool>& isSingularVertex = uvParam_->getCutMesh().getIsSingularVertex();
        if (isSingularVertex[i] || isSingularVertex[j] || isSingularVertex[k]) {
            isValidTriangle_[t] = false;
            continue;
        }
    }
}


void UVIsoTrace::traceIsolines() {
    const Mesh& mesh = uvParam_->getCutMesh().getCutMesh();
    const Eigen::VectorXd& uCoords = uvParam_->getU();
    const Eigen::VectorXd& vCoords = uvParam_->getV();

    // CutMesh references needed for seam crossing
    const CutMesh& cutMesh  = uvParam_->getCutMesh();
    const auto& cutToOrig   = cutMesh.getCutVertexToOriginal();
    const auto& origToCut   = cutMesh.getOriginalToCutVertices();
    const auto& cutEdgeSet  = cutMesh.getCutEdges();

    // ---------------------------------------------------------------
    // Helpers
    // ---------------------------------------------------------------

    // Coordinate of a cut-mesh vertex along U or V
    auto getCoord = [&](int vert, bool isU) -> double {
        return isU ? uCoords(vert) : vCoords(vert);
    };

    // Local edge index (0-2) in triangle t whose endpoints are va and vb.
    // Returns -1 if not found.
    auto findLocalEdge = [&](int t, int va, int vb) -> int {
        const Triangle& tri = mesh.triangles[t];
        for (int e = 0; e < 3; ++e) {
            int a = tri[e], b = tri[(e + 1) % 3];
            if ((a == va && b == vb) || (a == vb && b == va)) return e;
        }
        return -1;
    };

    // Find the two local edge indices crossed by isoline coord=c in triangle t.
    // An edge is "crossed" when the two endpoint coordinates strictly straddle c
    // (or one equals c while the other is on the other side, with non-zero span).
    // Returns {-1,-1} if fewer than 2 distinct crossings found (degenerate).
    auto findCrossedEdges = [&](int t, double c, bool isU) -> std::pair<int, int> {
        const Triangle& tri = mesh.triangles[t];
        int crossed[3], nc = 0;
        for (int e = 0; e < 3; ++e) {
            double ca = getCoord(tri[e], isU);
            double cb = getCoord(tri[(e + 1) % 3], isU);
            if (std::abs(ca - cb) < 1e-12) continue;   // skip degenerate edges
            double lo = std::min(ca, cb), hi = std::max(ca, cb);
            if (lo <= c && c <= hi) crossed[nc++] = e;
        }
        if (nc == 2) return {crossed[0], crossed[1]};
        return {-1, -1};
    };

    // Compute the spatial point and UV coordinates where isoline coord=c
    // crosses edge e of triangle t.
    auto computeCrossing = [&](int t, int e, double c, bool isU) -> std::pair<Point, Point> {
        const Triangle& tri = mesh.triangles[t];
        int va = tri[e], vb = tri[(e + 1) % 3];
        double ca = getCoord(va, isU), cb = getCoord(vb, isU);
        double param = (c - ca) / (cb - ca);
        param = std::max(0.0, std::min(1.0, param));

        const Point& pa = mesh.vertices[va];
        const Point& pb = mesh.vertices[vb];
        Point sp  = pa + (pb - pa) * param;
        Point uv  = { uCoords(va) + param * (uCoords(vb) - uCoords(va)),
                      vCoords(va) + param * (vCoords(vb) - vCoords(va)) };
        return {sp, uv};
    };

    // ---------------------------------------------------------------
    // Seam crossing: when the tracer hits a mesh-boundary edge that is
    // actually a cut/seam edge (not a true boundary), jump to the
    // matching triangle on the other side of the seam.
    //
    // After the BFS-based combing in CutMesh::combFieldDirections(), the
    // two triangles on opposite sides of a seam can have u-field directions
    // that differ by a multiple of 90°.  For a 90° / 270° mismatch the
    // isoline TYPE must switch (u-isoline ↔ v-isoline) so that the curve
    // continues in the same physical direction in 3-D space.  The isoline
    // VALUE is then taken from the other coordinate at the crossing point
    // on the far side of the seam.  For a 0° / 180° mismatch the type is
    // preserved and only the value is shifted by the seam translation τ.
    //
    // A second fix: the original code picked the first non-va copy of each
    // seam vertex, which is wrong at junction vertices (original vertices
    // with more than two cut-mesh copies).  We now search for the unique
    // copy-pair (cva, cvb) that co-occurs in the same triangle on the
    // other side of the seam.
    // ---------------------------------------------------------------
    auto crossSeam = [&](int t, int exitEdge, double c, bool isU,
                         int& outNextTri, int& outNewEntry, double& outC, bool& outIsU) -> bool
    {
        const Triangle& tri = mesh.triangles[t];
        const int va = tri[exitEdge], vb = tri[(exitEdge + 1) % 3];

        if (va >= static_cast<int>(cutToOrig.size()) ||
            vb >= static_cast<int>(cutToOrig.size())) return false;

        const int origVa = cutToOrig[va];
        const int origVb = cutToOrig[vb];
        if (origVa < 0 || origVb < 0) return false;

        // Only continue if this is a seam (cut) edge, not a true boundary
        CutMesh::EdgeKey ek(origVa, origVb);
        if (cutEdgeSet.find(ek) == cutEdgeSet.end()) return false;

        if (origVa >= static_cast<int>(origToCut.size()) ||
            origVb >= static_cast<int>(origToCut.size())) return false;

        // Find the triangle on the OTHER side of the seam.
        // Search all cut-mesh copies of origVa (≠ va) and origVb (≠ vb)
        // for a pair that appears together in the same triangle.  This
        // handles junction vertices (original vertices with more than two
        // cut-mesh copies) by always finding the face on the other side of
        // THIS specific seam edge rather than some other seam edge.
        int va2 = -1, vb2 = -1, nextTri = -1;
        for (int cva : origToCut[origVa]) {
            if (cva == va || nextTri >= 0) continue;
            auto [beg, end] = mesh.vertexTriangles.trianglesForVertex(cva);
            for (const int* it = beg; it != end && nextTri < 0; ++it) {
                const Triangle& cand = mesh.triangles[*it];
                for (int cvb : origToCut[origVb]) {
                    if (cvb == vb) continue;
                    if (cand[0] == cvb || cand[1] == cvb || cand[2] == cvb) {
                        va2 = cva; vb2 = cvb; nextTri = *it;
                        break;
                    }
                }
            }
        }
        if (nextTri < 0 || va2 < 0 || vb2 < 0) return false;

        const int newEntry = findLocalEdge(nextTri, va2, vb2);
        if (newEntry < 0) return false;

        // ---------------------------------------------------------------
        // Compute the interpolation parameter for the crossing point on
        // the exit edge (va→vb) using the current tracking coordinate.
        // The same parameter locates the identical physical point on the
        // paired edge (va2→vb2) on the other side of the seam.
        // ---------------------------------------------------------------
        const double ca_exit = isU ? uCoords(va) : vCoords(va);
        const double cb_exit = isU ? uCoords(vb) : vCoords(vb);
        const double exitSpan = cb_exit - ca_exit;
        const double param = (std::abs(exitSpan) > 1e-12)
                              ? std::max(0.0, std::min(1.0, (c - ca_exit) / exitSpan))
                              : 0.5;

        // Evaluate both u and v on the R-side (other side of seam) at the crossing.
        const double uR = uCoords(va2) + param * (uCoords(vb2) - uCoords(va2));
        const double vR = vCoords(va2) + param * (vCoords(vb2) - vCoords(va2));

        // ---------------------------------------------------------------
        // Determine the rotation index k_LR (the rotation from the L frame
        // to the R frame) using the combed per-triangle u-field directions.
        //
        // find_rotation_matrix(theta_R, theta_L) returns k_RL such that
        //   rotateVector(uField_R, k_RL) ≈ uField_L   (R→L)
        // The inverse (L→R) is k_LR = (4 - k_RL) % 4.
        //
        // The full transition formula is:
        //   [u_R, v_R]^T = R_{k_LR} * [u_L, v_L]^T + [τ^u, τ^v]^T
        //
        // rotateVector convention (Mesh.hxx):
        //   k=0: (u,v)→(u,v)       identity
        //   k=1: (u,v)→(-v,u)      90° CCW  → u-iso becomes v-iso: v_R = u_L + τ^v
        //   k=2: (u,v)→(-u,-v)     180°     → u-iso stays u-iso:   u_R = -u_L + τ^u
        //   k=3: (u,v)→(v,-u)      270° CCW → u-iso becomes v-iso: v_R = -u_L + τ^v
        //
        // For all k, the new tracking value equals the R-side coordinate
        // at the crossing point (since τ is already baked into uR/vR):
        //   k even (0,2): no axis swap → outIsU unchanged, outC = same coord on R side
        //   k odd  (1,3): axis swap   → outIsU flipped,   outC = other coord on R side
        // ---------------------------------------------------------------
        const std::vector<Point>& uFieldAll = cutMesh.getUField();
        const double theta_L = computeAngle(uFieldAll[t]);
        const double theta_R = computeAngle(uFieldAll[nextTri]);
        const int k_RL  = find_rotation_matrix(theta_R, theta_L);
        const int k_LR  = (4 - k_RL) % 4;
        const bool axisSwap = (k_LR % 2 != 0);

        outNextTri  = nextTri;
        outNewEntry = newEntry;
        outIsU      = axisSwap ? !isU : isU;
        // (isU != axisSwap) is true when the new tracked coordinate is u, false for v.
        outC        = (isU != axisSwap) ? uR : vR;

        return true;
    };

    // ---------------------------------------------------------------
    // Trace one side of an isoline starting from startTri.
    //
    //   fromEdge  – the local edge index in startTri through which the
    //               isoline *enters* (the exit edge is the other crossing).
    //   pushBack  – true  → append to back of deques (forward direction)
    //               false → prepend to front        (backward direction)
    //
    // curC may change when crossing seams (seam translation is applied).
    // ---------------------------------------------------------------
    auto traceOneSide = [&](IntegralCurve& curve,
                            int startTri, int fromEdge,
                            double c, bool isU,
                            bool pushBack) -> TraceTerminationReason
    {
        int    curTri    = startTri;
        int    entryEdge = fromEdge;
        double curC      = c;     // updated at each seam crossing
        bool   curIsU    = isU;   // may switch at seams with 90°/270° rotation

        for (int step = 0; step < MAX_TRACE_STEPS; ++step) {
            if (!isValidTriangle_[curTri])
                return TraceTerminationReason::HIT_INVALID_TRIANGLE;

            auto [eA, eB] = findCrossedEdges(curTri, curC, curIsU);
            if (eA < 0)
                return TraceTerminationReason::HIT_INVALID_TRIANGLE;  // degenerate

            // Exit through the crossing edge that is NOT the entry
            int exitEdge = (entryEdge == eA) ? eB : eA;

            auto [sp, uv] = computeCrossing(curTri, exitEdge, curC, curIsU);

            if (pushBack) {
                curve.points.push_back(sp);
                curve.uv_points.push_back(uv);
                curve.traversedFaces.push_back(curTri);
            } else {
                curve.points.push_front(sp);
                curve.uv_points.push_front(uv);
                curve.traversedFaces.push_front(curTri);
            }

            int nextTri = mesh.triangleAdjacency[curTri][exitEdge];
            if (nextTri < 0) {
                // Boundary hit – attempt a seam transition before declaring exit
                int    seamNext  = -1, seamEntry = -1;
                double seamC     = curC;
                bool   seamIsU   = curIsU;
                if (crossSeam(curTri, exitEdge, curC, curIsU, seamNext, seamEntry, seamC, seamIsU)
                        && isValidTriangle_[seamNext]) {
                    curTri    = seamNext;
                    entryEdge = seamEntry;
                    curC      = seamC;
                    curIsU    = seamIsU;   // update isoline type after rotation seam
                    continue;
                }
                return TraceTerminationReason::EXIT_BOUNDARY;
            }

            // Find the entry edge in the next triangle (the shared edge)
            const Triangle& ct = mesh.triangles[curTri];
            int va = ct[exitEdge], vb = ct[(exitEdge + 1) % 3];
            int newEntry = findLocalEdge(nextTri, va, vb);
            if (newEntry < 0)
                return TraceTerminationReason::HIT_INVALID_TRIANGLE;  // mesh inconsistency

            curTri    = nextTri;
            entryEdge = newEntry;
        }
        return TraceTerminationReason::MAX_STEPS_REACHED;
    };

    // Pick the most severe of two termination reasons
    auto combineReasons = [](TraceTerminationReason a, TraceTerminationReason b) {
        if (a == TraceTerminationReason::MAX_STEPS_REACHED ||
            b == TraceTerminationReason::MAX_STEPS_REACHED)
            return TraceTerminationReason::MAX_STEPS_REACHED;
        if (a == TraceTerminationReason::HIT_INVALID_TRIANGLE ||
            b == TraceTerminationReason::HIT_INVALID_TRIANGLE)
            return TraceTerminationReason::HIT_INVALID_TRIANGLE;
        return TraceTerminationReason::EXIT_BOUNDARY;
    };

    // ---------------------------------------------------------------
    // Trace U isolines
    // ---------------------------------------------------------------
    TrianglesInterstedByUIsolines.clear();

    for (int iu = 0; iu < nU_; ++iu) {
        const double uVal = uvMin_[0] + iu * deltaU_;

        std::vector<int> candidates = uIntervalTree_.query(uVal);

        // Build a per-isoline pending set of valid candidate triangles.
        // Also accumulate into the global membership set.
        std::unordered_set<int> pending;
        for (int t : candidates) {
            if (!isValidTriangle_[t]) continue;
            pending.insert(t);
            TrianglesInterstedByUIsolines.insert(t);
        }

        // Start new curve segments until every candidate triangle is covered.
        // This handles non-injective mappings where one u-value crosses the mesh
        // in multiple disconnected strips.
        while (!pending.empty()) {
            const int seedTri = *pending.begin();

            auto [eA, eB] = findCrossedEdges(seedTri, uVal, true);
            if (eA < 0) { pending.erase(seedTri); continue; }  // degenerate seed

            IntegralCurve curve;
            curve.traceType        = IsoTraceType::U_ISOLINE;
            curve.coordinateValue  = uVal;

            // Forward  – enter seedTri from eB side, exit through eA, push_back
            auto r1 = traceOneSide(curve, seedTri, eB, uVal, true, true);
            // Backward – enter seedTri from eA side, exit through eB, push_front
            auto r2 = traceOneSide(curve, seedTri, eA, uVal, true, false);

            curve.terminationReason = combineReasons(r1, r2);

            // Remove all traversed triangles from the pending set
            for (int t : curve.traversedFaces) pending.erase(t);

            if (curve.terminationReason == TraceTerminationReason::EXIT_BOUNDARY)
                integralCurves_.push_back(std::move(curve));
        }
    }

    // ---------------------------------------------------------------
    // Trace V isolines
    // ---------------------------------------------------------------
    for (int iv = 0; iv < nV_; ++iv) {
        const double vVal = uvMin_[1] + iv * deltaV_;

        std::vector<int> candidates = vIntervalTree_.query(vVal);

        std::unordered_set<int> pending;
        for (int t : candidates) {
            if (isValidTriangle_[t]) pending.insert(t);
        }

        while (!pending.empty()) {
            const int seedTri = *pending.begin();

            auto [eA, eB] = findCrossedEdges(seedTri, vVal, false);
            if (eA < 0) { pending.erase(seedTri); continue; }

            IntegralCurve curve;
            curve.traceType        = IsoTraceType::V_ISOLINE;
            curve.coordinateValue  = vVal;

            auto r1 = traceOneSide(curve, seedTri, eB, vVal, false, true);
            auto r2 = traceOneSide(curve, seedTri, eA, vVal, false, false);

            curve.terminationReason = combineReasons(r1, r2);

            for (int t : curve.traversedFaces) pending.erase(t);

            if (curve.terminationReason == TraceTerminationReason::EXIT_BOUNDARY)
                integralCurves_.push_back(std::move(curve));
        }
    }
}


bool UVIsoTrace::writeVTK(const std::string& filename) const
{
    const Mesh& mesh        = uvParam_->getCutMesh().getCutMesh();
    const Eigen::VectorXd& uCoords = uvParam_->getU();
    const Eigen::VectorXd& vCoords = uvParam_->getV();
    const int nMeshVerts = static_cast<int>(mesh.vertices.size());
    const int nTris      = static_cast<int>(mesh.triangles.size());

    // Collect curves with at least 2 points (degenerate single-point curves
    // cannot be represented as VTK polylines).
    std::vector<const IntegralCurve*> curves;
    curves.reserve(integralCurves_.size());
    for (const auto& c : integralCurves_) {
        if (c.points.size() >= 2) curves.push_back(&c);
    }

    // Count total curve vertices and connectivity size.
    int nCurveVerts = 0;
    int totalConn   = 4 * nTris;  // each triangle: "3 i j k" = 4 tokens
    for (const auto* c : curves) {
        const int n = static_cast<int>(c->points.size());
        nCurveVerts += n;
        totalConn   += 1 + n;     // polyline: "n p0 p1 ... p_{n-1}"
    }

    const int nPts   = nMeshVerts + nCurveVerts;
    const int nCells = nTris + static_cast<int>(curves.size());

    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "[writeVTK] Cannot open: " << filename << "\n";
        return false;
    }
    ofs << std::scientific << std::setprecision(10);

    // ------------------------------------------------------------------
    // Header
    // ------------------------------------------------------------------
    ofs << "# vtk DataFile Version 2.0\n"
        << "CrossGen UV Isolines\n"
        << "ASCII\n"
        << "DATASET UNSTRUCTURED_GRID\n\n";

    // ------------------------------------------------------------------
    // Points: mesh vertices first, then all curve sample points
    // ------------------------------------------------------------------
    ofs << "POINTS " << nPts << " double\n";
    for (int i = 0; i < nMeshVerts; ++i) {
        const Point& p = mesh.vertices[i];
        ofs << p[0] << " " << p[1] << " 0.0\n";
    }
    for (const auto* c : curves)
        for (const Point& p : c->points)
            ofs << p[0] << " " << p[1] << " 0.0\n";
    ofs << "\n";

    // ------------------------------------------------------------------
    // Cells
    // ------------------------------------------------------------------
    ofs << "CELLS " << nCells << " " << totalConn << "\n";

    // Mesh triangles
    for (const Triangle& tri : mesh.triangles)
        ofs << "3 " << tri[0] << " " << tri[1] << " " << tri[2] << "\n";

    // Curve polylines – point indices are offset by nMeshVerts
    {
        int offset = nMeshVerts;
        for (const auto* c : curves) {
            const int n = static_cast<int>(c->points.size());
            ofs << n;
            for (int k = 0; k < n; ++k) ofs << " " << (offset + k);
            ofs << "\n";
            offset += n;
        }
    }
    ofs << "\n";

    // ------------------------------------------------------------------
    // Cell types: 5 = VTK_TRIANGLE, 4 = VTK_POLY_LINE
    // ------------------------------------------------------------------
    ofs << "CELL_TYPES " << nCells << "\n";
    for (int i = 0; i < nTris; ++i)              ofs << "5\n";
    for (size_t i = 0; i < curves.size(); ++i)   ofs << "4\n";
    ofs << "\n";

    // ------------------------------------------------------------------
    // Cell data
    // ------------------------------------------------------------------
    ofs << "CELL_DATA " << nCells << "\n";

    // cell_type: 0 = mesh triangle, 1 = u-isoline, 2 = v-isoline
    ofs << "SCALARS cell_type int 1\nLOOKUP_TABLE default\n";
    for (int i = 0; i < nTris; ++i) ofs << "0\n";
    for (const auto* c : curves)
        ofs << (c->traceType == IsoTraceType::U_ISOLINE ? "1" : "2") << "\n";
    ofs << "\n";

    // coord_value: constant u/v isoline value (0.0 for mesh triangles)
    ofs << "SCALARS coord_value double 1\nLOOKUP_TABLE default\n";
    for (int i = 0; i < nTris; ++i) ofs << "0.0\n";
    for (const auto* c : curves) ofs << c->coordinateValue << "\n";
    ofs << "\n";

    // ------------------------------------------------------------------
    // Point data: UV parameterization values at every point
    // ------------------------------------------------------------------
    ofs << "POINT_DATA " << nPts << "\n";

    ofs << "SCALARS u_coord double 1\nLOOKUP_TABLE default\n";
    for (int i = 0; i < nMeshVerts; ++i) ofs << uCoords(i) << "\n";
    for (const auto* c : curves)
        for (const Point& uv : c->uv_points) ofs << uv[0] << "\n";
    ofs << "\n";

    ofs << "SCALARS v_coord double 1\nLOOKUP_TABLE default\n";
    for (int i = 0; i < nMeshVerts; ++i) ofs << vCoords(i) << "\n";
    for (const auto* c : curves)
        for (const Point& uv : c->uv_points) ofs << uv[1] << "\n";
    ofs << "\n";

    std::cerr << "[writeVTK] Written " << filename
              << "  (" << nTris << " triangles, "
              << curves.size() << " isocurves, "
              << nPts << " points)\n";
    return true;
}


void UVIsoTrace::printQueries() const {

    //print u queries
    for (double uVal = uvMin_[0]; uVal <= uvMax_[0]; uVal += deltaU_) {
        std::vector<int> triIndices = uIntervalTree_.query(uVal);
        std::cout << "u = " << uVal << ": " << triIndices.size() << " intersecting triangles\n";
        // print the triangle indices for debugging
        for (int idx : triIndices) {
            std::cout << "  Triangle " << idx << (isValidTriangle_[idx] ? " (valid)" : " (invalid)") << "\n";
        }   
    }

    // print v queries
    for (double vVal = uvMin_[1]; vVal <= uvMax_[1]; vVal += deltaV_) {
        std::vector<int> triIndices = vIntervalTree_.query(vVal);
        std::cout << "v = " << vVal << ": " << triIndices.size() << " intersecting triangles\n";
        // print the triangle indices for debugging
        for (int idx : triIndices) {
            std::cout << "  Triangle " << idx << (isValidTriangle_[idx] ? " (valid)" : " (invalid)") << "\n";
        }      
    }
}