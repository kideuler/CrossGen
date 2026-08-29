// Two inclusions embedded in a matrix: a horizontal stadium (oblong) and a
// square
// tilted 30 degrees.
//
// Closed interfaces, no junction vertices anywhere -- this domain isolates the
// other half of the multi-material argument, the claim that each material
// diffuses independently. Neither inclusion may perturb the field in the
// matrix beyond what its own boundary tangency imposes, and the matrix may not
// reach through an inclusion to the far side.
//
// The two inclusions ask different questions. The stadium is a closed
// interface bounding a disk-topology material, so its interior field has to
// carry index +1 on its own and should show four +1/4 singularities inside;
// the matrix meanwhile has to wrap around it. A circle would ask the same
// question but is a bad test: its rotational symmetry leaves the angular
// positions of those four singularities undetermined, so they land wherever
// arithmetic noise puts them and cannot be compared between runs or between
// discretizations. Stretching it into an oblong -- two semicircular caps of
// radius 0.28 joined by straight sides of length 0.24, so a 0.80 x 0.56
// footprint -- breaks the symmetry and pins them to the two ends of the long
// axis, which is a prediction the method can actually be scored against. The
// straight sides also give the interface two flat stretches where the field
// must run parallel to a line rather than to a continuously turning normal,
// which the pure ellipse never tests. The tilted square asks
// whether the matrix field can hold a 30 degree rotation against the
// axis-aligned outer wall across a short distance, and puts four sharp corners
// of one material in the interior of another.
Point(1) = {0, 0, 0, 1.0};
Point(2) = {2, 0, 0, 1.0};
Point(3) = {2, 2, 0, 1.0};
Point(4) = {0, 2, 0, 1.0};

// Stadium: cap centres at x = 0.65 +/- 0.12, cap radius r = 0.28.
Point(5) = {0.77, 0.70, 0, 1.0};   // right cap centre
Point(6) = {0.53, 0.70, 0, 1.0};   // left cap centre
Point(7) = {1.05, 0.70, 0, 1.0};   // right extreme
Point(8) = {0.77, 0.98, 0, 1.0};   // top right, start of upper straight
Point(9) = {0.53, 0.98, 0, 1.0};   // top left, end of upper straight
Point(14) = {0.25, 0.70, 0, 1.0};  // left extreme
Point(15) = {0.53, 0.42, 0, 1.0};  // bottom left, start of lower straight
Point(16) = {0.77, 0.42, 0, 1.0};  // bottom right, end of lower straight

Point(10) = {1.517114, 1.787110, 0, 1.0};   // tilted square, corner at 75 deg
Point(11) = {0.962890, 1.467114, 0, 1.0};   // 165 deg
Point(12) = {1.282886, 0.912890, 0, 1.0};   // 255 deg
Point(13) = {1.837110, 1.232886, 0, 1.0};   // 345 deg

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// The two caps, each split into quarter circles because a GEO-kernel Circle
// may not span 180 degrees or more, alternating with the straight sides.
Circle(5) = {7, 5, 8};     // right cap, upper quarter
Line(6)   = {8, 9};        // upper straight side
Circle(7) = {9, 6, 14};    // left cap, upper quarter
Circle(8) = {14, 6, 15};   // left cap, lower quarter
Line(13)  = {15, 16};      // lower straight side
Circle(14) = {16, 5, 7};   // right cap, lower quarter

Line(9)  = {10, 11};
Line(10) = {11, 12};
Line(11) = {12, 13};
Line(12) = {13, 10};

Curve Loop(1) = {1, 2, 3, 4};      // outer wall
Curve Loop(2) = {5, 6, 7, 8, 13, 14};   // stadium interface
Curve Loop(3) = {9, 10, 11, 12};   // tilted square interface

Plane Surface(1) = {1, 2, 3};      // matrix, with both inclusions as holes
Plane Surface(2) = {2};
Plane Surface(3) = {3};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
