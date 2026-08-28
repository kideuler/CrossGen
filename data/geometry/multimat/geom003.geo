// Orthogonal T-junction: the CONTROL of the junction-angle family
// (geom003-geom007), which holds the outer disk fixed and varies only the
// angles of the radial cuts.
//
// Three materials meet at the centre. The interfaces are a full diameter
// (direction 0 deg) plus one stem at 90 deg, so every interface direction is
// congruent mod 90 -- a single axis-aligned cross satisfies all three tangency
// constraints at once. The vertex-based field is therefore WELL-POSED here.
// That is the point of this domain: it is the case where P1 gets away with it,
// and the baseline for reading geom004 and geom006, which differ from it only
// in the stem angle.
//
// Each cut meets the circle at a right angle, so the three boundary junctions
// are benign too, and the centre is the only place anything can go wrong.

Point(1) = {0, 0, 0, 1.0};  // junction: all materials meet here
Point(2) = {1, 0, 0, 1.0};  // cut at 0 deg
Point(3) = {6.12323399574e-17, 1, 0, 1.0};  // cut at 90 deg
Point(4) = {-1, 1.22464679915e-16, 0, 1.0};  // cut at 180 deg
Point(5) = {-1.83697019872e-16, -1, 0, 1.0};

// Radial material interfaces, shared by the two sectors on either side.
Line(1) = {1, 2};
Line(2) = {1, 3};
Line(3) = {1, 4};

// Outer boundary arcs.
Circle(4) = {2, 1, 3};
Circle(5) = {3, 1, 4};
Circle(6) = {4, 1, 5};
Circle(7) = {5, 1, 2};

Curve Loop(1) = {1, 4, -2};
Plane Surface(1) = {1};
Curve Loop(2) = {2, 5, -3};
Plane Surface(2) = {2};
Curve Loop(3) = {3, 6, 7, -1};
Plane Surface(3) = {3};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
