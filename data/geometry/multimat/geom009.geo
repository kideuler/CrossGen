// Two inclusions embedded in a matrix: a horizontal ellipse and a square
// tilted 30 degrees.
//
// Closed interfaces, no junction vertices anywhere -- this domain isolates the
// other half of the multi-material argument, the claim that each material
// diffuses independently. Neither inclusion may perturb the field in the
// matrix beyond what its own boundary tangency imposes, and the matrix may not
// reach through an inclusion to the far side.
//
// The two inclusions ask different questions. The ellipse is a closed smooth
// interface bounding a disk-topology material, so its interior field has to
// carry index +1 on its own and should show four +1/4 singularities inside;
// the matrix meanwhile has to wrap around it. A circle would ask the same
// question but is a bad test: its rotational symmetry leaves the angular
// positions of those four singularities undetermined, so they land wherever
// arithmetic noise puts them and cannot be compared between runs or between
// discretizations. Stretching it to a 0.45 x 0.26 horizontal ellipse breaks
// the symmetry and pins them to the two ends of the major axis, which is a
// prediction the method can actually be scored against. The tilted square asks
// whether the matrix field can hold a 30 degree rotation against the
// axis-aligned outer wall across a short distance, and puts four sharp corners
// of one material in the interior of another.
Point(1) = {0, 0, 0, 1.0};
Point(2) = {2, 0, 0, 1.0};
Point(3) = {2, 2, 0, 1.0};
Point(4) = {0, 2, 0, 1.0};

Point(5) = {0.65, 0.70, 0, 1.0};   // ellipse centre
Point(6) = {1.10, 0.70, 0, 1.0};   // +major (a = 0.45)
Point(7) = {0.65, 0.96, 0, 1.0};   // +minor (b = 0.26)
Point(8) = {0.20, 0.70, 0, 1.0};   // -major
Point(9) = {0.65, 0.44, 0, 1.0};   // -minor

Point(10) = {1.517114, 1.787110, 0, 1.0};   // tilted square, corner at 75 deg
Point(11) = {0.962890, 1.467114, 0, 1.0};   // 165 deg
Point(12) = {1.282886, 0.912890, 0, 1.0};   // 255 deg
Point(13) = {1.837110, 1.232886, 0, 1.0};   // 345 deg

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Quarter ellipse arcs. GEO-kernel Ellipse takes {start, centre, point on the
// major axis, end} and, like Circle, may not span 180 degrees or more, so the
// interface is built from four quarters. Point 6 serves as both the first
// start point and the major-axis direction marker for all four.
Ellipse(5) = {6, 5, 6, 7};
Ellipse(6) = {7, 5, 6, 8};
Ellipse(7) = {8, 5, 6, 9};
Ellipse(8) = {9, 5, 6, 6};

Line(9)  = {10, 11};
Line(10) = {11, 12};
Line(11) = {12, 13};
Line(12) = {13, 10};

Curve Loop(1) = {1, 2, 3, 4};      // outer wall
Curve Loop(2) = {5, 6, 7, 8};      // ellipse interface
Curve Loop(3) = {9, 10, 11, 12};   // tilted square interface

Plane Surface(1) = {1, 2, 3};      // matrix, with both inclusions as holes
Plane Surface(2) = {2};
Plane Surface(3) = {3};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
