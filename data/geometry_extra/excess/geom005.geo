// Orthogonal quadruple point: two perpendicular diameters, four materials.
//
// The second control of the family. Both interface directions (0 and 90 deg)
// are congruent mod 90, so even with four materials meeting at one vertex an
// axis-aligned cross satisfies every tangency constraint. Four materials at a
// vertex is not by itself what breaks a vertex-based field -- the angles are.
// Pairs with geom006 exactly as geom003 pairs with geom004.

Point(1) = {0, 0, 0, 1.0};  // junction: all materials meet here
Point(2) = {1, 0, 0, 1.0};  // cut at 0 deg
Point(3) = {6.12323399574e-17, 1, 0, 1.0};  // cut at 90 deg
Point(4) = {-1, 1.22464679915e-16, 0, 1.0};  // cut at 180 deg
Point(5) = {-1.83697019872e-16, -1, 0, 1.0};  // cut at 270 deg

// Radial material interfaces, shared by the two sectors on either side.
Line(1) = {1, 2};
Line(2) = {1, 3};
Line(3) = {1, 4};
Line(4) = {1, 5};

// Outer boundary arcs.
Circle(5) = {2, 1, 3};
Circle(6) = {3, 1, 4};
Circle(7) = {4, 1, 5};
Circle(8) = {5, 1, 2};

Curve Loop(1) = {1, 5, -2};
Plane Surface(1) = {1};
Curve Loop(2) = {2, 6, -3};
Plane Surface(2) = {2};
Curve Loop(3) = {3, 7, -4};
Plane Surface(3) = {3};
Curve Loop(4) = {4, 8, -1};
Plane Surface(4) = {4};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
