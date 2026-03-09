#include "tracing/SeparatrixTrace.hxx"

SeparatrixTrace::SeparatrixTrace(std::shared_ptr<CrossField> cf, bool useActualSingularityCoordinates) : crossField(cf) {
    separatrices.clear();
    singularities.clear();

    // Initialize isSingularTriangle lookup vector
    isSingularTriangle.assign(crossField->mesh->triangles.size(), false);
    for (const auto& [triIdx, cfIndex] : crossField->singularTriangles) {
        if (triIdx >= 0 && triIdx < static_cast<int>(isSingularTriangle.size())) {
            isSingularTriangle[triIdx] = true;
        }
    }

    // Fill singularity vector and map from cross field data
    for (const auto& [triIdx, cfIndex] : crossField->singularTriangles) {
        Singularity s;
        s.triangleIndex = triIdx;
        s.singularityIndex = cfIndex; // Map cross field index to singularity index

        s.numPorts = (cfIndex < 0.0) ? 5 : 3;

        // compute coordinate as either the triangle centroid or the actual singularity center by solving the linear system
        const Triangle &tri = crossField->mesh->triangles[triIdx];
        if (useActualSingularityCoordinates) {
            Eigen::Matrix2d A;
            std::complex<double> u1 = crossField->u_k[tri[0]];
            std::complex<double> u2 = crossField->u_k[tri[1]];
            std::complex<double> u3 = crossField->u_k[tri[2]];

            A << std::real(u1) - std::real(u3), std::real(u2) - std::real(u3),
                std::imag(u1) - std::imag(u3), std::imag(u2) - std::imag(u3);

            Eigen::Vector2d b;
            b << -std::real(u3), -std::imag(u3);
            // uv is solution
            Eigen::Vector2d uv = A.colPivHouseholderQr().solve(b);
            
            s.barycentric = {uv[0], uv[1], 1.0 - uv[0] - uv[1]};

            // find the corresponding point in the triangle
            Point p1 = crossField->mesh->vertices[tri[0]];
            Point p2 = crossField->mesh->vertices[tri[1]];
            Point p3 = crossField->mesh->vertices[tri[2]];
            
            s.coordinates = p1 * s.barycentric[0] + p2 * s.barycentric[1] + p3 * s.barycentric[2];
        } else {
            // Use triangle centroid as singularity coordinate
            Point bc{0.0, 0.0};
            for (int i = 0; i < 3; ++i) {
                bc[0] += crossField->mesh->vertices[tri[i]][0];
                bc[1] += crossField->mesh->vertices[tri[i]][1];
            }
            bc[0] /= 3.0;
            bc[1] /= 3.0;
            s.coordinates = bc;
            s.barycentric = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        }

        // get fields angle
        int vIdx = tri[0]; // can pick any vertex since it's a singularity
        std::complex<double> u4 = crossField->u_k[vIdx];
        double fieldAngle = std::arg(u4) / 4.0;
        Point refVec = crossField->mesh->vertices[vIdx] - s.coordinates;
        double refAngle = std::atan2(refVec[1], refVec[0]);
        s.refAngle = refAngle;

        double alpha = makeAngleSamePhase(refAngle, fieldAngle) - refAngle;
        s.alpha = wrap_pi(alpha);
        
        singularities.push_back(s);
        singularityMap[triIdx] = singularities.size() - 1;

        // create a separatrix for each port emanating from this singularity
        double dalpha = (s.numPorts == 5) ? 2.0 * M_PI / 5.0 : 2.0 * M_PI / 3.0;

        for (int port = 0; port < s.numPorts; ++port) {
            // Port direction: reference angle + field alignment offset + port spacing
            double portAngle = wrap_pi(s.refAngle + s.alpha + port * dalpha);

            // Create the origin trace point at the singularity center
            TracePoint tp0;
            tp0.face_id = triIdx;
            tp0.global_pos = s.coordinates;
            tp0.barycentric = s.barycentric;
            tp0.field_angle = portAngle;
            tp0.trace_direction = portAngle;

            // Trace the ray from the singularity center to the triangle boundary
            auto [exitPos, exitEdge, exitT, exitNeighbor] = rayEdgeIntersection(triIdx, s.coordinates, portAngle);

            // Create the separatrix
            Separatrix sep;
            sep.id = static_cast<int>(separatrices.size());
            sep.origin_singularity_id = static_cast<int>(singularities.size()) - 1;
            sep.origin_singularity_port = port;
            sep.active = true;
            sep.termination_reason = TerminationReason::RUNNING;

            // Add the singularity center as the first point
            sep.path.push_back(tp0);

            if (exitEdge >= 0) {
                // Build the second trace point at the exit edge
                TracePoint tp1;
                tp1.global_pos = exitPos;
                tp1.face_id = triIdx;
                tp1.barycentric = globalToBarycentric(triIdx, exitPos);
                tp1.edge_id = crossField->mesh->triangleEdges[triIdx][exitEdge];
                tp1.local_edge_index = exitEdge;
                tp1.edge_crossing_t = exitT;

                // Compute the field angle at the exit point by interpolating vertex angles
                const Triangle &exitTri = crossField->mesh->triangles[triIdx];
                double theta0 = std::arg(crossField->u_k[exitTri[0]]) / 4.0;
                double theta1 = std::arg(crossField->u_k[exitTri[1]]) / 4.0;
                double theta2 = std::arg(crossField->u_k[exitTri[2]]) / 4.0;
                theta1 = makeAngleSamePhase(theta0, theta1);
                theta2 = makeAngleSamePhase(theta0, theta2);
                double thetaInterp = theta0 * tp1.barycentric[0] + theta1 * tp1.barycentric[1] + theta2 * tp1.barycentric[2];
                tp1.field_angle = makeAngleSamePhase(portAngle, thetaInterp);

                // Actual trace direction from center to exit
                Point traceVec = exitPos - s.coordinates;
                tp1.trace_direction = std::atan2(traceVec[1], traceVec[0]);

                sep.path.push_back(tp1);

                // If the exit edge is on the boundary, deactivate immediately
                if (exitNeighbor < 0) {
                    sep.active = false;
                    sep.termination_reason = TerminationReason::EXIT_BOUNDARY;
                }
            } else {
                // Could not find an exit edge (shouldn't happen for interior singularities)
                sep.active = false;
                sep.termination_reason = TerminationReason::UNDEFINED;
            }

            // Record the separatrix ID in the singularity's port list
            singularities.back().portSeparatrixIds[port] = sep.id;

            separatrices.push_back(std::move(sep));
        }
    }
}

int SeparatrixTrace::findPhaseDifference(double referenceAngle, double angleCandidate) {
    // Compute the phase difference in terms of multiples of 90 degrees (pi/2 radians).
    // Returns k such that angleCandidate + k * (pi/2) is closest to referenceAngle.
    // k can be -2, -1, 0, 1, or 2 (we only need to search ±2 since that covers ±π)

    double diff = wrap_pi(referenceAngle - angleCandidate);
    int k = static_cast<int>(std::round(diff / (M_PI / 2.0)));
    
    // Clamp k to [-2, 2] range - beyond this we're more than π away which shouldn't happen
    // with properly wrapped angles
    if (k > 2) k = 2;
    if (k < -2) k = -2;

    return k;
}

double SeparatrixTrace::makeAngleSamePhase(double referenceAngle, double angleCandidate){
    // We want to find the angle that is equivalent to angleCandidate mod pi/2 that is closest to referenceAngle.
    // This ensures that we follow the same "branch" of the cross field direction at this point for continuity.

    int k = findPhaseDifference(referenceAngle, angleCandidate);
    return angleCandidate + k * (M_PI / 2.0);
}

std::array<double, 3> SeparatrixTrace::globalToBarycentric(int triangleIndex, const Point &p) {
    const Triangle &tri = crossField->mesh->triangles[triangleIndex];

    const Point &v0 = crossField->mesh->vertices[tri[0]];
    const Point &v1 = crossField->mesh->vertices[tri[1]];
    const Point &v2 = crossField->mesh->vertices[tri[2]];

    // Robust barycentric via area coordinates.
    Point v0v1 = v1 - v0;
    Point v0v2 = v2 - v0;
    Point v0p = p - v0;

    double denom = cross2(v0v1, v0v2);
    if (std::abs(denom) < 1e-30) {
        return {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    }

    double l1 = cross2(v0p, v0v2) / denom; // weight for v1
    double l2 = cross2(v0v1, v0p) / denom; // weight for v2
    double l0 = 1.0 - l1 - l2;

    std::array<double, 3> bary = {l0, l1, l2};

    // Detect and snap coordinates that are very close to 0 (point on edge)
    // Also clamp coordinates that are slightly outside [0, 1] due to numerical error
    for (int i = 0; i < 3; ++i) {
        if (std::abs(bary[i]) < EPS_BARY) {
            bary[i] = 0.0;
        } else if (bary[i] < 0.0 && bary[i] > -EPS_BARY * 10) {
            bary[i] = 0.0;
        } else if (bary[i] > 1.0 && bary[i] < 1.0 + EPS_BARY * 10) {
            bary[i] = 1.0;
        }
    }

    // Renormalize to ensure sum equals 1 after snapping
    double sum = bary[0] + bary[1] + bary[2];
    if (std::abs(sum) > 1e-30 && std::abs(sum - 1.0) > 1e-14) {
        bary[0] /= sum;
        bary[1] /= sum;
        bary[2] /= sum;
    }

    return bary;
}


std::tuple<Point, int, double, int> SeparatrixTrace::rayEdgeIntersection(int triangleIndex, const Point &origin, double direction, int excludeEdge) {
    const Triangle &tri = crossField->mesh->triangles[triangleIndex];

    Point dir = {std::cos(direction), std::sin(direction)};

    Point v[3] = {crossField->mesh->vertices[tri[0]], crossField->mesh->vertices[tri[1]], crossField->mesh->vertices[tri[2]]};

    double bestT = std::numeric_limits<double>::infinity();
    int bestEdge = -1;
    double bestU = 0.0;
    Point bestIp = origin;

    for (int e = 0; e < 3; ++e) {
        // Skip the excluded edge (e.g., the entry edge)
        if (e == excludeEdge) continue;
        
        int a = e;
        int b = (e + 1) % 3;
        Point e0 = v[a];
        Point e1 = v[b];
        Point edge = e1 - e0;

        // Solve origin + t*dir = e0 + u*edge
        double denom = cross2(dir, edge);
        if (std::abs(denom) < 1e-20) continue;

        Point eo = e0 - origin;
        double t = cross2(eo, edge) / denom;
        double u = cross2(eo, dir) / denom;

        if (t <= EPS_INTERSECT_T) continue;
        if (u < -1e-12 || u > 1.0 + 1e-12) continue;

        if (t < bestT) {
            bestT = t;
            bestEdge = e;
            bestU = std::max(0.0, std::min(1.0, u));
            bestIp = origin + (dir * t);
        }
    }

    int neighbor = -1;
    if (bestEdge >= 0) neighbor = crossField->mesh->triangleAdjacency[triangleIndex][bestEdge];

    return {bestIp, bestEdge, bestU, neighbor};
}

std::tuple<Point, int, double, int> SeparatrixTrace::rayEdgeIntersection(int triangleIndex, const Point &origin, const Point &direction, int excludeEdge) {
    const Triangle &tri = crossField->mesh->triangles[triangleIndex];

    // Normalize the direction vector
    Point dir = normalizeP(direction);

    Point v[3] = {crossField->mesh->vertices[tri[0]], crossField->mesh->vertices[tri[1]], crossField->mesh->vertices[tri[2]]};

    double bestT = std::numeric_limits<double>::infinity();
    int bestEdge = -1;
    double bestU = 0.0;
    Point bestIp = origin;

    for (int e = 0; e < 3; ++e) {
        // Skip the excluded edge (e.g., the entry edge)
        if (e == excludeEdge) continue;
        
        int a = e;
        int b = (e + 1) % 3;
        Point e0 = v[a];
        Point e1 = v[b];
        Point edge = e1 - e0;

        // Solve origin + t*dir = e0 + u*edge
        double denom = cross2(dir, edge);
        if (std::abs(denom) < 1e-20) continue;

        Point eo = e0 - origin;
        double t = cross2(eo, edge) / denom;
        double u = cross2(eo, dir) / denom;

        if (t <= EPS_INTERSECT_T) continue;
        if (u < -1e-12 || u > 1.0 + 1e-12) continue;

        if (t < bestT) {
            bestT = t;
            bestEdge = e;
            bestU = std::max(0.0, std::min(1.0, u));
            bestIp = origin + (dir * t);
        }
    }

    int neighbor = -1;
    if (bestEdge >= 0) neighbor = crossField->mesh->triangleAdjacency[triangleIndex][bestEdge];

    return {bestIp, bestEdge, bestU, neighbor};
}

int SeparatrixTrace::findCrossedEdge(const std::array<double, 3>& b_prev, const std::array<double, 3>& b_next, double& t_exit) {
    int exit_edge = -1;
    double min_t = 2.0; // Initialize to a value strictly > 1.0

    for (int i = 0; i < 3; ++i) {
        if (b_next[i] < 0.0) {
            // Calculate the fraction of the step where this coordinate hit exactly 0.0
            // Denominator is guaranteed > 0 because b_prev[i] >= 0 and b_next[i] < 0
            double t = b_prev[i] / (b_prev[i] - b_next[i]);
            
            if (t < min_t) {
                min_t = t;
                // If coordinate i goes negative, we crossed the edge opposite to vertex i.
                // In your convention, the edge opposite vertex i is (i + 1) % 3.
                exit_edge = (i + 1) % 3;
            }
        }
    }
    
    t_exit = (exit_edge != -1) ? min_t : 1.0;
    return exit_edge;
}

std::tuple<bool, Point> SeparatrixTrace::edgeEdgeIntersection(const Point& p0, const Point& p1, const Point& p2, const Point& p3) {
    // Computes the intersection of segment p0->p1 with segment p2->p3.
    // Returns (true, intersection_point) if the segments intersect, (false, {0,0}) otherwise.
    //
    // Parameterize: p0 + t*(p1 - p0) = p2 + u*(p3 - p2), solve for t and u.
    // Both t and u must be in [0, 1] for the segments to intersect.

    Point d1 = {p1[0] - p0[0], p1[1] - p0[1]}; // direction of segment 1
    Point d2 = {p3[0] - p2[0], p3[1] - p2[1]}; // direction of segment 2

    double denom = cross2(d1, d2);
    if (std::abs(denom) < 1e-20) {
        // Segments are parallel (or degenerate)
        return {false, {0.0, 0.0}};
    }

    Point d02 = {p2[0] - p0[0], p2[1] - p0[1]};
    double t = cross2(d02, d2) / denom;
    double u = cross2(d02, d1) / denom;

    if (t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0) {
        Point ip = {p0[0] + t * d1[0], p0[1] + t * d1[1]};
        return {true, ip};
    }

    return {false, {0.0, 0.0}};
}


void SeparatrixTrace::stepHeuns(Separatrix& sep) {
    // Heun's method step: https://en.wikipedia.org/wiki/Heun%27s_method
    // 1. Euler step to get a candidate point
    // 2. Interpolate field direction at candidate point
    // 3. Average the original and candidate directions for better accuracy

    if (!sep.active) return;

    if (sep.path.empty()) {
        std::cerr << "Warning: stepHeuns called on separatrix with empty path." << std::endl;
        sep.active = false;
        sep.termination_reason = TerminationReason::UNDEFINED;
        return;
    }

    const TracePoint& current = sep.path.back();

    if (current.local_edge_index == -1) { 
        std::cerr << "Warning: stepHeuns called with TracePoint that has no local_edge_index set. This may lead to incorrect behavior." << std::endl; 
        sep.active = false;
        sep.termination_reason = TerminationReason::UNDEFINED;
        return; 
    }

    // find next triangle and barycentric coordinates using the current local_edge_index and direction
    int nextTri = crossField->mesh->triangleAdjacency[current.face_id][current.local_edge_index];
    if (nextTri < 0) {
        // Hit boundary, deactivate the separatrix
        sep.active = false;
        sep.termination_reason = TerminationReason::EXIT_BOUNDARY;
        return;
    }

    // use edge_crossing_t to find barycentric coordinates in next triangle
    // The edge in the current triangle is defined by local vertices local_edge_index and (local_edge_index+1)%3
    // We need to find the corresponding edge in the neighbor triangle
    
    const Triangle &currentTri = crossField->mesh->triangles[current.face_id];
    const Triangle &nextTriVerts = crossField->mesh->triangles[nextTri];
    
    // Get the two vertex indices of the shared edge from the current triangle
    int v_a = currentTri[current.local_edge_index];
    int v_b = currentTri[(current.local_edge_index + 1) % 3];
    
    // Find which local edge in the next triangle corresponds to this shared edge
    // The shared edge will have the same two vertex indices but possibly in reverse order
    int neighborEdge = -1;
    bool reversed = false;
    for (int e = 0; e < 3; ++e) {
        int nv_a = nextTriVerts[e];
        int nv_b = nextTriVerts[(e + 1) % 3];
        if (nv_a == v_a && nv_b == v_b) {
            neighborEdge = e;
            reversed = false;
            break;
        } else if (nv_a == v_b && nv_b == v_a) {
            neighborEdge = e;
            reversed = true;
            break;
        }
    }
    
    if (neighborEdge < 0) {
        std::cerr << "Error: Could not find shared edge in neighbor triangle." << std::endl;
        sep.active = false;
        sep.termination_reason = TerminationReason::UNDEFINED;
        return;
    }
    
    // Compute barycentric coordinates in the new triangle
    // The crossing point lies on edge (neighborEdge, (neighborEdge+1)%3) of the new triangle
    // If reversed, the parameter t is flipped: t_new = 1 - t_old
    double t = reversed ? (1.0 - current.edge_crossing_t) : current.edge_crossing_t;
    
    // Barycentric coordinates: the point is on the edge from vertex neighborEdge to vertex (neighborEdge+1)%3
    // So the third vertex (neighborEdge+2)%3 has weight 0
    // Vertex neighborEdge has weight (1-t), vertex (neighborEdge+1)%3 has weight t
    std::array<double, 3> bary = {0.0, 0.0, 0.0};
    bary[neighborEdge] = 1.0 - t;
    bary[(neighborEdge + 1) % 3] = t;
    bary[(neighborEdge + 2) % 3] = 0.0;

    // find angles of the three vertices of the next triangle
    std::complex<double> u0 = crossField->u_k[crossField->mesh->triangles[nextTri][0]];
    std::complex<double> u1 = crossField->u_k[crossField->mesh->triangles[nextTri][1]];
    std::complex<double> u2 = crossField->u_k[crossField->mesh->triangles[nextTri][2]];

    double theta0 = std::arg(u0)/4.0;
    double theta1 = std::arg(u1)/4.0;
    double theta2 = std::arg(u2)/4.0;

    theta1 = makeAngleSamePhase(theta0, theta1);
    theta2 = makeAngleSamePhase(theta0, theta2);
    
    double theta_entry = theta0 * bary[0] + theta1 * bary[1] + theta2 * bary[2];

    // Match the interpolated field angle to the field angle from the previous triangle's exit point.
    // current.field_angle is the field at the exit of the previous triangle (same geometric point as our entry).
    // This ensures continuity of the field across edges, which is crucial when the cross-field
    // representation (u^4) loses branch information (e.g., when u^4 ≈ 1 everywhere).
    int k = findPhaseDifference(current.field_angle, theta_entry);

    // Debug phase matching
    double theta_entry_raw = theta_entry;
    double adjusted = theta_entry + k * (M_PI / 2.0);
    if (std::abs(wrap_pi(adjusted - current.field_angle)) > 0.5) {
        std::cerr << "Debug phase: raw_entry=" << theta_entry_raw 
                  << " current.field_angle=" << current.field_angle 
                  << " k=" << k << " adjusted=" << adjusted << std::endl;
    }

    // Debug: check if the phase matching introduces a large jump
    double adjusted_theta = theta_entry + k * (M_PI / 2.0);
    double angle_diff = std::abs(wrap_pi(adjusted_theta - current.field_angle));
    if (angle_diff > M_PI / 4.0 + 0.1) {
        std::cerr << "Warning: Large phase jump detected! k=" << k 
                  << ", theta_entry=" << theta_entry 
                  << ", current.field_angle=" << current.field_angle 
                  << ", adjusted=" << adjusted_theta 
                  << ", diff=" << angle_diff << std::endl;
    }

    // make all angles same branch
    theta_entry += k * (M_PI / 2.0);
    theta0 += k * (M_PI / 2.0);
    theta1 += k * (M_PI / 2.0);
    theta2 += k * (M_PI / 2.0);

    // compute the global position of the entry point in the new triangle
    const Point &p0 = crossField->mesh->vertices[nextTriVerts[0]];
    const Point &p1 = crossField->mesh->vertices[nextTriVerts[1]];
    const Point &p2 = crossField->mesh->vertices[nextTriVerts[2]];
    Point entryPos = p0 * bary[0] + p1 * bary[1] + p2 * bary[2];

    // =========== HEUN'S METHOD ===========
    // Step 1: Euler step using direction at entry point
    // Exclude the entry edge (neighborEdge) to avoid backtracking
    Point v1 = {std::cos(theta_entry), std::sin(theta_entry)};
    auto [exitPos1, exitEdge1, exitT1, exitNeighbor1] = rayEdgeIntersection(nextTri, entryPos, v1, neighborEdge);

    if (exitEdge1 < 0) {
        // The computed direction doesn't lead to a valid exit.
        // This can happen due to numerical issues or if the direction points back toward the entry edge.
        // Try using the field angle directly without phase adjustment as a fallback.
        std::cerr << "Warning: stepHeuns Euler step failed in triangle " << nextTri 
                  << ", theta_entry=" << theta_entry << ", field_angle=" << current.field_angle << std::endl;
        
        // Fallback: try using current.field_angle directly (the direction we were traveling)
        Point v_fallback = {std::cos(current.field_angle), std::sin(current.field_angle)};
        auto [fbExitPos, fbExitEdge, fbExitT, fbExitNeighbor] = rayEdgeIntersection(nextTri, entryPos, v_fallback, neighborEdge);
        
        if (fbExitEdge < 0) {
            std::cerr << "Warning: stepHeuns fallback also failed in triangle " << nextTri << std::endl;
            sep.active = false;
            sep.termination_reason = TerminationReason::UNDEFINED;
            return;
        }
        
        // Use fallback result
        TracePoint next;
        next.face_id = nextTri;
        next.barycentric = globalToBarycentric(nextTri, fbExitPos);
        next.global_pos = fbExitPos;
        next.field_angle = current.field_angle; // Keep the same field angle
        next.trace_direction = std::atan2(fbExitPos[1] - entryPos[1], fbExitPos[0] - entryPos[0]);
        next.edge_id = crossField->mesh->triangleEdges[nextTri][fbExitEdge];
        next.local_edge_index = fbExitEdge;
        next.edge_crossing_t = fbExitT;
        sep.path.push_back(next);
        sep.visited_edges.insert(next.edge_id);

        // Check if the exit edge is on the boundary
        if (fbExitNeighbor < 0) {
            sep.active = false;
            sep.termination_reason = TerminationReason::EXIT_BOUNDARY;
        }
        return;
    }

    // Step 2: Compute field direction at the trial exit point
    std::array<double, 3> exitBary1 = globalToBarycentric(nextTri, exitPos1);
    double theta_exit1 = theta0 * exitBary1[0] + theta1 * exitBary1[1] + theta2 * exitBary1[2];

    // Step 3: Average the entry and exit directions (Heun's corrector)
    double theta_avg = 0.5 * (theta_entry + theta_exit1);
    Point v_avg = {std::cos(theta_avg), std::sin(theta_avg)};

    // Step 4: Retrace with the averaged direction (also exclude entry edge)
    auto [exitPos, exitEdge, exitT, exitNeighbor] = rayEdgeIntersection(nextTri, entryPos, v_avg, neighborEdge);

    if (exitEdge < 0) {
        // Fall back to Euler step if Heun's fails
        exitPos = exitPos1;
        exitEdge = exitEdge1;
        exitT = exitT1;
        exitNeighbor = exitNeighbor1;
    }

    // Compute barycentric coordinates of the final exit point
    std::array<double, 3> exitBary = globalToBarycentric(nextTri, exitPos);

    // Compute the field angle at the EXIT point for phase continuity in next step
    double exitTheta = theta0 * exitBary[0] + theta1 * exitBary[1] + theta2 * exitBary[2];

    // Compute the actual traced direction (from entry to exit)
    // This is the true direction we moved, used for phase continuity in the next step
    Point traceVec = exitPos - entryPos;
    double actualTraceDir = std::atan2(traceVec[1], traceVec[0]);

    // Get the global edge ID for the exit edge
    int exitEdgeId = crossField->mesh->triangleEdges[nextTri][exitEdge];

    // Build the result TracePoint and append to the separatrix path
    TracePoint next;
    next.face_id = nextTri;
    next.barycentric = exitBary;
    next.global_pos = exitPos;
    next.field_angle = exitTheta;  // Use field angle at EXIT point
    next.trace_direction = actualTraceDir; // Store actual traced direction for phase continuity
    next.edge_id = exitEdgeId;
    next.local_edge_index = exitEdge;
    next.edge_crossing_t = exitT;

    sep.path.push_back(next);
    sep.visited_edges.insert(exitEdgeId);

    // Check if the exit edge is on the boundary
    if (crossField->mesh->triangleAdjacency[nextTri][exitEdge] < 0) {
        sep.active = false;
        sep.termination_reason = TerminationReason::EXIT_BOUNDARY;
    }

    //add seperatrix id to triangleSeparatrixMap for next triangle
    triangleSeparatrixMap[nextTri].first.push_back(sep.id);
    if (triangleSeparatrixMap[nextTri].first.size() > 1) {
        triangleSeparatrixMap[nextTri].second = true; // mark as multiple separatrices passing through this triangle
        Intersections.push(nextTri); // add next triangle to Intersections queue for intersection checking
    }

    // check whether the next triangle is singular
    if (sep.active) {
        nextTri = crossField->mesh->triangleAdjacency[next.face_id][next.local_edge_index];
        sep.in_singularity_zone = isSingularTriangle[nextTri];
    }
    
}


void SeparatrixTrace::stepViertel(Separatrix& sep, bool stopAtOrthogonal) {
    if (!sep.active) return;
    if (sep.path.empty()) {
        sep.active = false;
        sep.termination_reason = TerminationReason::UNDEFINED;
        return;
    }

    const TracePoint& current = sep.path.back();
    if (current.local_edge_index == -1) { 
        sep.active = false;
        sep.termination_reason = TerminationReason::UNDEFINED;
        return; 
    }

    // 1. Identify upcoming triangle
    int nextTri = crossField->mesh->triangleAdjacency[current.face_id][current.local_edge_index];
    if (nextTri < 0) {
        sep.active = false;
        sep.termination_reason = TerminationReason::EXIT_BOUNDARY;
        return;
    }

    auto sigIt = singularityMap.find(nextTri);

    // If the next triangle is not actually singular, fall back to Heun's method
    if (sigIt == singularityMap.end()) {
        sep.in_singularity_zone = false;
        stepHeuns(sep);
        return;
    }

    // 3. Find neighbor edge and entry coordinates
    const Triangle &currentTri = crossField->mesh->triangles[current.face_id];
    const Triangle &nextTriVerts = crossField->mesh->triangles[nextTri];
    int v_a = currentTri[current.local_edge_index];
    int v_b = currentTri[(current.local_edge_index + 1) % 3];
    
    int neighborEdge = -1;
    bool reversed = false;
    for (int e = 0; e < 3; ++e) {
        int nv_a = nextTriVerts[e];
        int nv_b = nextTriVerts[(e + 1) % 3];
        if (nv_a == v_a && nv_b == v_b) {
            neighborEdge = e; reversed = false; break;
        } else if (nv_a == v_b && nv_b == v_a) {
            neighborEdge = e; reversed = true; break;
        }
    }
    
    if (neighborEdge < 0) {
        sep.active = false;
        sep.termination_reason = TerminationReason::UNDEFINED;
        return;
    }
    
    double t_entry = reversed ? (1.0 - current.edge_crossing_t) : current.edge_crossing_t;
    std::array<double, 3> bary = {0.0, 0.0, 0.0};
    bary[neighborEdge] = 1.0 - t_entry;
    bary[(neighborEdge + 1) % 3] = t_entry;
    
    Point p0 = crossField->mesh->vertices[nextTriVerts[0]];
    Point p1 = crossField->mesh->vertices[nextTriVerts[1]];
    Point p2 = crossField->mesh->vertices[nextTriVerts[2]];
    Point q = p0 * bary[0] + p1 * bary[1] + p2 * bary[2];

    const Singularity& s = singularities[sigIt->second];

    // singularity center in global coordinates
    Point sc = s.coordinates;
    double d_type = s.singularityIndex*4;

    double r_q = normP(q - sc); // radius 
    double theta_raw = wrap_pi(std::atan2(q[1] - sc[1], q[0] - sc[0])); // angle of entry point around singularity

    // find the port whose angle is closest to theta_q in the clockwise direction
    // i.e., the port with the smallest positive circular distance wrap_pi(theta_q - portAngle)
    double bestDist = std::numeric_limits<double>::infinity();
    double dalpha = (s.numPorts == 5) ? 2.0 * M_PI / 5.0 : 2.0 * M_PI / 3.0;
    int bestPort = -1;
    double bestTheta = 0.0;
    for (int port = 0; port < s.numPorts; ++port) {
        double portAngle = wrap_pi(s.refAngle + s.alpha + port * dalpha);
        // Circular distance: how far theta_raw is ahead of portAngle (going counterclockwise)
        double dist = wrap_pi(theta_raw - portAngle);
        if (dist < 0) dist += 2.0 * M_PI; // map to [0, 2*pi) so "behind" ports get large distance
        if (dist > 1e-10 && dist < bestDist) { // exclude dist≈0 (exactly on the port)
            bestDist = dist;
            bestPort = port;
            bestTheta = portAngle;
        }
    }

    double theta_q = theta_raw - bestTheta; // angle of entry point relative to the chosen port direction, must lie in range [0, 2*pi]
    if (theta_q < 0) theta_q += 2.0 * M_PI;

    double M = (4.0 - d_type) / 8.0;
    double M_inv = 1.0 / M;
    double M_qr = std::pow(r_q, M); // conformally mapped radius
    double phi_q = theta_q * M; // conformally mapped angle

    // 1. Map the incoming FIELD ANGLE to the conformal tangent angle (psi)
    // We strictly use current.field_angle instead of current.trace_direction to shed numerical noise
    double gamma_local = current.trace_direction - bestTheta;
    double psi = (M - 1.0) * theta_q + gamma_local;

    // 2. Determine which hyperbolic family this streamline belongs to.
    double delta_psi = wrap_pi(psi + phi_q);
    bool isSweepFamily = (std::abs(delta_psi) < M_PI / 4.0 || std::abs(delta_psi) > 3.0 * M_PI / 4.0);

    // 3. Determine the direction to step phi.
    // In polar coordinates, the sign of dphi ALWAYS matches the sign of sin(tangent_angle - position_angle).
    double dphi_sign = (std::sin(psi - phi_q) >= 0) ? 1.0 : -1.0;
    double dphi = dphi_sign * dphi_singularity_zone;

    // 4. Calculate the correct hyperbola constant based on the chosen family
    double A_q;
    if (isSweepFamily) {
        A_q = M_qr * M_qr * std::sin(phi_q) * std::cos(phi_q);
    } else {
        A_q = M_qr * M_qr * std::cos(2.0 * phi_q);
    }

    double phi = phi_q;
    double rho = M_qr;
    Point q_current = q;
    double phi_next = phi;
    double rho_next = rho;
    Point q_next = q_current;
    double r_next, theta_next;
    bool exitedTriangle = false;

    // Validate the initial stepping direction with a trial step.
    // If the first step immediately exits through the entry edge (same edge we came
    // from), the dphi sign is wrong — flip it. This guards against incorrect conformal
    // tangent angle causing the hyperbola to be traversed outward.
    {
        double phi_trial = phi_q + dphi;
        double rho_trial;
        if (isSweepFamily) {
            rho_trial = std::sqrt(std::abs(A_q) / std::abs(std::sin(phi_trial) * std::cos(phi_trial)));
        } else {
            double cos2 = std::cos(2.0 * phi_trial);
            rho_trial = std::sqrt(std::abs(A_q) / std::max(std::abs(cos2), 1e-12));
        }
        double r_trial = std::pow(rho_trial, M_inv);
        double theta_trial = phi_trial * M_inv;
        Point q_trial = sc + Point{r_trial * std::cos(theta_trial + bestTheta),
                                   r_trial * std::sin(theta_trial + bestTheta)};

        std::array<double, 3> bary_trial = globalToBarycentric(nextTri, q_trial);
        // Check if trial point exits through the entry edge (neighborEdge)
        // The entry edge is 'neighborEdge'; the barycentric coordinate opposite to that edge's
        // starting vertex goes negative when we cross back out.
        // Edge e is between vertices e and (e+1)%3, opposite vertex is (e+2)%3.
        int oppositeVertex = (neighborEdge + 2) % 3;
        if (bary_trial[oppositeVertex] < -EPS_BARY) {
            // Trial step went backward through the entry edge — flip direction
            dphi_sign = -dphi_sign;
            dphi = -dphi;
        }
    }

    // get angle for intersecting orthogonal separatrix family if stopAtOrthogonal is true
    int bestPortOrthogonal = bestPort + int(dphi_sign > 0);
    if (bestPortOrthogonal >= s.numPorts) { bestPortOrthogonal -= s.numPorts; }
    if (bestPortOrthogonal < 0) { bestPortOrthogonal += s.numPorts; }
    int orthoSepid = s.portSeparatrixIds[bestPortOrthogonal];
    // get points 1 and 2 from the separatrix
    Point ortho_p0 = separatrices[orthoSepid].path[0].global_pos;
    Point ortho_p1 = separatrices[orthoSepid].path[1].global_pos;

    int i = 0;
    while (i < maxStepsInSingularityZone && !exitedTriangle) {
        // compute next angle on hyperbola
        phi_next = phi + dphi;

        // Safe radius computation depending on the family
        if (isSweepFamily) {
            rho_next = std::sqrt(std::abs(A_q) / (std::sin(phi_next) * std::cos(phi_next)));
        } else {
            double cos2 = std::cos(2.0 * phi_next);
            rho_next = std::sqrt(std::abs(A_q) / std::abs(cos2));
        }

        // apply inverse mapping to get back to original domain
        r_next = std::pow(rho_next, M_inv);
        theta_next = phi_next * M_inv;

        Point q_next = sc + Point{r_next * std::cos(theta_next + bestTheta), r_next * std::sin(theta_next + bestTheta)};

        // check if we have crossed a separatrix orthogonal to the current one, if stopAtOrthogonal is true
        if (stopAtOrthogonal) {
            auto [intersects, ip] = edgeEdgeIntersection(q_current, q_next, ortho_p0, ortho_p1);
            if (intersects) {
                // construct a TracePoint at the intersection and add to path, then terminate
                TracePoint tp;
                tp.face_id = nextTri;
                tp.global_pos = ip;
                tp.barycentric = globalToBarycentric(nextTri, ip);
                
                sep.path.push_back(tp);
                sep.active = false;
                sep.termination_reason = TerminationReason::ORTHOGONAL_TO_SINGULARITY_SEPARATRIX;
                return;
            }
        }

        // determine if we have exited the triangle by checking the barycentric coordinates of q_next in nextTri
        std::array<double, 3> bary_next = globalToBarycentric(nextTri, q_next);
        if (bary_next[0] < 0.0 || bary_next[1] < 0.0 || bary_next[2] < 0.0) {
            exitedTriangle = true;
        }

        // create trace point and add to path
        if (!exitedTriangle) { // if we didnt exit the triangle, we can safely add the point
            TracePoint tp;
            tp.face_id = nextTri;
            tp.global_pos = q_next;
            tp.barycentric = globalToBarycentric(nextTri, q_next);
            tp.field_angle = bestTheta; // use the port direction as the field angle for this point since we're following that direction
            tp.trace_direction = std::atan2(q_next[1] - sc[1], q_next[0] - sc[0]);
            sep.path.push_back(tp);

            // update phi and rho for next iteration
            phi = phi_next;
            rho = rho_next;
            q_current = q_next;
        } else { // if we exited the triangle we need to find the inters
            // use bary_next to find which edge we exited and the t parameter along that edge
            double t_exit;
            std::array<double, 3> bary_prev = globalToBarycentric(nextTri, q_current); // barycentrics of the entry point
            int exit_edge = findCrossedEdge(bary_prev, bary_next, t_exit);
            if (exit_edge < 0) {
                std::cerr << "Error: Could not find exit edge in stepViertel. bary_prev=" << bary_prev[0] << " " << bary_prev[1] << " " << bary_prev[2] << " bary_next=" << bary_next[0] << " " << bary_next[1] << " " << bary_next[2] << " q_current=" << q_current[0] << " " << q_current[1] << std::endl;
                sep.active = false;
                sep.termination_reason = TerminationReason::UNDEFINED;
                return;
            }

            // evaluate the exit point on the edge using t_exit
            Point exitPos = q_current * (1.0 - t_exit) + q_next * t_exit;

            // add the exit point as the final trace point
            TracePoint tp;
            tp.face_id = nextTri;
            tp.global_pos = exitPos;
            tp.barycentric = globalToBarycentric(nextTri, exitPos);
            tp.trace_direction = std::atan2(exitPos[1] - q_current[1], exitPos[0] - q_current[0]);
            tp.edge_id = crossField->mesh->triangleEdges[nextTri][exit_edge];
            tp.local_edge_index = exit_edge;
            tp.edge_crossing_t = t_exit;

            // field angle at the exit point by interpolating vertex angles
            const Triangle &exitTri = crossField->mesh->triangles[nextTri];
            double theta0 = std::arg(crossField->u_k[exitTri[0]]) / 4.0;
            double theta1 = std::arg(crossField->u_k[exitTri[1]]) / 4.0;
            double theta2 = std::arg(crossField->u_k[exitTri[2]]) / 4.0;
            theta1 = makeAngleSamePhase(theta0, theta1);
            theta2 = makeAngleSamePhase(theta0, theta2);
            double thetaInterp = theta0 * tp.barycentric[0] + theta1 * tp.barycentric[1] + theta2 * tp.barycentric[2];
            tp.field_angle = thetaInterp;

            sep.path.push_back(tp);

            // Update in_singularity_zone: check if the triangle beyond the exit edge is also singular
            int beyondTri = crossField->mesh->triangleAdjacency[nextTri][exit_edge];
            if (beyondTri < 0) {
                sep.active = false;
                sep.termination_reason = TerminationReason::EXIT_BOUNDARY;
                sep.in_singularity_zone = false;
            } else {
                sep.in_singularity_zone = isSingularTriangle[beyondTri];
            }
        }

        i++;
    }
        
}


void SeparatrixTrace::stepAndCheck() {
    // Step all active separatrices and check for intersections
    finishedTracing = true;
    for (Separatrix& sep : separatrices) {
        if (!sep.active) continue;

        // Step the separatrix 
        if (!sep.in_singularity_zone) {
            stepHeuns(sep);
        } else {
            stepViertel(sep, true); // stop at orthogonal crossings to singularity separatrices
        }
        finishedTracing = false; // if any separatrix is still active, we're not finished;
    }

}