// Symmetric Y-junction at 120 degrees: three equal sectors.
//
// The physically canonical triple point -- three materials at mutual 120 deg
// is the equilibrium configuration when the three interfacial tensions are
// equal, so this is the junction that actually turns up in a multi-material
// simulation rather than one constructed to be awkward.
//
// Directions 90, 210, 330 deg reduce mod 90 to 0, 30, 60: maximally far from
// congruent, and symmetric, so no averaging scheme can favour one material
// without breaking the symmetry the geometry has. Unlike geom004 and geom006,
// where a reviewer might argue the angles were chosen adversarially, this one
// is chosen by physics.

Point(1) = {0, 0, 0, 1.0};  // junction: all materials meet here
Point(2) = {6.12323399574e-17, 1, 0, 1.0};  // cut at 90 deg
Point(3) = {-0.866025403784, 0.5, 0, 1.0};
Point(4) = {-0.866025403784, -0.5, 0, 1.0};  // cut at 210 deg
Point(5) = {-1.83697019872e-16, -1, 0, 1.0};
Point(6) = {0.866025403784, -0.5, 0, 1.0};  // cut at 330 deg
Point(7) = {0.866025403784, 0.5, 0, 1.0};

// Radial material interfaces, shared by the two sectors on either side.
Line(1) = {1, 2};
Line(2) = {1, 4};
Line(3) = {1, 6};

// Outer boundary arcs.
Circle(4) = {2, 1, 3};
Circle(5) = {3, 1, 4};
Circle(6) = {4, 1, 5};
Circle(7) = {5, 1, 6};
Circle(8) = {6, 1, 7};
Circle(9) = {7, 1, 2};

Curve Loop(1) = {1, 4, 5, -2};
Plane Surface(1) = {1};
Curve Loop(2) = {2, 6, 7, -3};
Plane Surface(2) = {2};
Curve Loop(3) = {3, 8, 9, -1};
Plane Surface(3) = {3};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
