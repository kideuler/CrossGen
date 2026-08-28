// Skewed quadruple point: geom005 with the second diameter rotated to 60 deg,
// giving sector angles 60/120/60/120.
//
// Interface directions 0 and 60 deg are not congruent mod 90, so the four
// crosses meeting at the centre cannot be reconciled into one. This is the
// worst interior junction in the set: four materials, no consistent vertex
// value, and two of the four sectors narrow enough that a misplaced
// singularity has nowhere to go.

Point(1) = {0, 0, 0, 1.0};  // junction: all materials meet here
Point(2) = {1, 0, 0, 1.0};  // cut at 0 deg
Point(3) = {0.5, 0.866025403784, 0, 1.0};  // cut at 60 deg
Point(4) = {-0.5, 0.866025403784, 0, 1.0};
Point(5) = {-1, 1.22464679915e-16, 0, 1.0};  // cut at 180 deg
Point(6) = {-0.5, -0.866025403784, 0, 1.0};  // cut at 240 deg
Point(7) = {0.5, -0.866025403784, 0, 1.0};

// Radial material interfaces, shared by the two sectors on either side.
Line(1) = {1, 2};
Line(2) = {1, 3};
Line(3) = {1, 5};
Line(4) = {1, 6};

// Outer boundary arcs.
Circle(5) = {2, 1, 3};
Circle(6) = {3, 1, 4};
Circle(7) = {4, 1, 5};
Circle(8) = {5, 1, 6};
Circle(9) = {6, 1, 7};
Circle(10) = {7, 1, 2};

Curve Loop(1) = {1, 5, -2};
Plane Surface(1) = {1};
Curve Loop(2) = {2, 6, 7, -3};
Plane Surface(2) = {2};
Curve Loop(3) = {3, 8, -4};
Plane Surface(3) = {3};
Curve Loop(4) = {4, 9, 10, -1};
Plane Surface(4) = {4};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
