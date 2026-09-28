// A lined tunnel in dipping rock. The opening (a void) is D-shaped -- a
// semicircular crown of radius 0.5 on vertical walls, with a flat invert --
// and is lined with concrete (mat 1), 0.12 thick round the crown and walls
// and 0.16 under the invert. Round the lining is a ring of grouted rock
// (mat 2), 0.25 thick, smooth all the way round. The host is three rock
// strata (mats 3, 4 and 5 from the top) whose bedding dips 25 degrees to the
// right; the upper bedding plane runs straight through where the tunnel was
// driven, so it now ends on the grouted ring at two points, and the lower one
// passes 0.16 beneath the ring.
//
// What this domain is for. A hole inside nested material rings. The void is
// bounded by one material and nested inside two rings (lining, grout), each
// an annulus whose Euler characteristic is zero,
// and the rings' inner boundaries carry corners (the invert corners of the
// opening and of the lining) that their outer boundaries do not. A field has
// to follow three nested closed curves of different smoothness around a
// hole while the host around them is layered at an angle, and the only
// junctions are where the upper bedding plane lands on the ring: at the
// crown, at 149 degrees round from the springline, and on the right wall.
// Interface directions there are {59, 155} and {90, 155} -- mod 90, {59, 65}
// and {0, 65} -- so both are ILL-POSED, the first only by six degrees, the
// second by twenty-five. The bedding meets the side walls at 65/115 degrees.
//
// Resolution at NP = 100 (h = 0.042 here): the lining is 2.9 elements thick
// (3.8 under the invert), the grouted ring 6, and the gap between the ring
// and the lower bedding plane 3.8. All arc centres lie inside the domain.

Point(1) = {-2.1, -2.2, 0, 1.0};
Point(2) = { 2.1, -2.2, 0, 1.0};
Point(3) = { 2.1,  2.0, 0, 1.0};
Point(4) = {-2.1,  2.0, 0, 1.0};

// Bedding planes y = 0.35 - x tan 25 (upper) and y = -1.15 - x tan 25
// (lower), where they meet the walls.
Point(5) = {-2.1,  1.329246, 0, 1.0};
Point(6) = {-2.1, -0.170754, 0, 1.0};
Point(7) = { 2.1, -0.629246, 0, 1.0};
Point(8) = { 2.1, -2.129246, 0, 1.0};

Point(10) = {0, 0.25, 0, 1.0};   // centre of every crown arc

// The opening.
Point(11) = { 0.5, -0.25, 0, 1.0};
Point(12) = { 0.5,  0.25, 0, 1.0};
Point(13) = { 0.0,  0.75, 0, 1.0};
Point(14) = {-0.5,  0.25, 0, 1.0};
Point(15) = {-0.5, -0.25, 0, 1.0};

// Outer face of the lining.
Point(21) = { 0.62, -0.41, 0, 1.0};
Point(22) = { 0.62,  0.25, 0, 1.0};
Point(23) = { 0.00,  0.87, 0, 1.0};
Point(24) = {-0.62,  0.25, 0, 1.0};
Point(25) = {-0.62, -0.41, 0, 1.0};

// Outer face of the grouted ring: the lining's face offset by 0.25, its two
// lower corners rounded about the lining's corners. Points 31 and 34 are
// where the upper bedding plane meets it.
Point(31) = { 0.87,        -0.0556876626, 0, 1.0};
Point(32) = { 0.87,         0.25,         0, 1.0};
Point(33) = { 0.00,         1.12,         0, 1.0};
Point(34) = {-0.7458955163, 0.6978167915, 0, 1.0};
Point(35) = {-0.87,         0.25,         0, 1.0};
Point(36) = {-0.87,        -0.41,         0, 1.0};
Point(37) = {-0.62,        -0.66,         0, 1.0};
Point(38) = { 0.62,        -0.66,         0, 1.0};
Point(39) = { 0.87,        -0.41,         0, 1.0};

// Outer boundary, counter-clockwise.
Line(1) = {1, 2};
Line(2) = {2, 8};
Line(3) = {8, 7};
Line(4) = {7, 3};
Line(5) = {3, 4};
Line(6) = {4, 5};
Line(7) = {5, 6};
Line(8) = {6, 1};

// Bedding planes.
Line(9)  = {5, 34};   // upper, left of the tunnel
Line(10) = {31, 7};   // upper, right of the tunnel
Line(11) = {6, 8};    // lower

// Grouted ring's outer face, counter-clockwise.
Line(21)   = {31, 32};
Circle(22) = {32, 10, 33};
Circle(23) = {33, 10, 34};
Circle(24) = {34, 10, 35};
Line(25)   = {35, 36};
Circle(26) = {36, 25, 37};
Line(27)   = {37, 38};
Circle(28) = {38, 21, 39};
Line(29)   = {39, 31};

// Lining's outer face.
Line(31)   = {21, 22};
Circle(32) = {22, 10, 23};
Circle(33) = {23, 10, 24};
Line(34)   = {24, 25};
Line(35)   = {25, 21};

// The opening.
Line(41)   = {11, 12};
Circle(42) = {12, 10, 13};
Circle(43) = {13, 10, 14};
Line(44)   = {14, 15};
Line(45)   = {15, 11};

Curve Loop(11) = {41, 42, 43, 44, 45};                   // opening
Curve Loop(12) = {31, 32, 33, 34, 35};                   // lining's face
Curve Loop(13) = {21, 22, 23, 24, 25, 26, 27, 28, 29};   // ring's face

Plane Surface(1) = {12, 11};    // lining
Plane Surface(2) = {13, 12};    // grouted ring
Curve Loop(3) = {5, 6, 9, -23, -22, -21, 10, 4};
Plane Surface(3) = {3};         // upper stratum
Curve Loop(4) = {-10, -29, -28, -27, -26, -25, -24, -9, 7, 11, 3};
Plane Surface(4) = {4};         // middle stratum
Curve Loop(5) = {8, 1, 2, -11};
Plane Surface(5) = {5};         // lower stratum

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
