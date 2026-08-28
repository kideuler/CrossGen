// Oblique T-junction: geom003 with the stem swung from 90 deg to 60 deg, and
// nothing else changed.
//
// Interface directions are now 0 deg (the diameter) and 60 deg (the stem).
// Mod 90 those are 0 and 60, so no single cross is tangent to all three
// interfaces: the vertex-based field is ILL-POSED at the centre, and any P1
// treatment has to average conflicting constraints or duplicate the vertex.
// Face-based DOFs give each of the three incident sectors its own cross and
// the question does not arise.
//
// The minimal-pair comparison against geom003 is the cleanest figure in the
// multi-material set: one parameter changes, well-posedness flips.

Point(1) = {0, 0, 0, 1.0};  // junction: all materials meet here
Point(2) = {1, 0, 0, 1.0};  // cut at 0 deg
Point(3) = {0.5, 0.866025403784, 0, 1.0};  // cut at 60 deg
Point(4) = {-0.5, 0.866025403784, 0, 1.0};
Point(5) = {-1, 1.22464679915e-16, 0, 1.0};  // cut at 180 deg
Point(6) = {-1.83697019872e-16, -1, 0, 1.0};

// Radial material interfaces, shared by the two sectors on either side.
Line(1) = {1, 2};
Line(2) = {1, 3};
Line(3) = {1, 5};

// Outer boundary arcs.
Circle(4) = {2, 1, 3};
Circle(5) = {3, 1, 4};
Circle(6) = {4, 1, 5};
Circle(7) = {5, 1, 6};
Circle(8) = {6, 1, 2};

Curve Loop(1) = {1, 4, -2};
Plane Surface(1) = {1};
Curve Loop(2) = {2, 5, 6, -3};
Plane Surface(2) = {2};
Curve Loop(3) = {3, 7, 8, -1};
Plane Surface(3) = {3};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
