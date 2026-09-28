// One pole of a surface-magnet synchronous machine: a 90-degree sector of a
// four-pole, twelve-slot motor, cut on two radial lines where the geometry
// repeats. From the centre out: the steel shaft (mat 6, r < 0.2), the rotor
// iron (mat 1, to r = 0.55), one arc-shaped permanent magnet (mat 2, 0.06
// thick, spanning 70 of the pole's 90 degrees), the air gap (mat 3, 0.03),
// and the stator iron (mat 4, from r = 0.64 to 1) with three semi-closed
// slots, each holding a coil (mat 5) under an air pocket that opens into the
// gap through a 0.05-wide neck between tooth-tip shoes. The teeth are
// parallel-sided (0.12 wide), so the slots widen outward. This is the
// geometry every electrical-machine finite-element code (FEMM, Maxwell, ...)
// is first run on.
//
// What this domain is for. Two things. First, a THIN LAYER
// THAT IS ALSO THE CONNECTING TISSUE: the air is one region -- the gap, the
// two half pole-gaps beside the magnet, and three slot pockets reached
// through necks -- that winds between every other material and is 3
// elements thick in the gap and 5 across a neck (at NP = 100; h = 0.01
// here). Second, a DENSITY OF SMALL FEATURES IN POLAR ARRANGEMENT: every
// arc is centred on the one point, so radial and circumferential directions
// are the natural cross everywhere, and a field that is not globally polar
// shows it.
//
// Junctions: the magnet's corners on the rotor, and the places the shaft,
// rotor, gap and stator reach the radial cut lines, are right angles in the
// polar sense (directions congruent mod 90, WELL-POSED). The slot flanks are
// parallel to the teeth, not radial, so they meet the coil tops at 75 and 105
// degrees: the six junctions of stator, coil and air are ILL-POSED by 15
// degrees, and the tooth-tip shoes turn through 73 degrees onto the necks.
// Copper is one material on three surfaces. The sector's apex is a 90-degree
// corner, so nothing here is a full circle.

// ---- generated geometry ----
Point(1) = {0.0000000000, 0.0000000000, 0, 1.0};   // centre of every arc
Point(2) = {0.2000000000, 0.0000000000, 0, 1.0};   // shaft
Point(3) = {0.0000000000, 0.2000000000, 0, 1.0};
Point(4) = {0.5500000000, 0.0000000000, 0, 1.0};   // rotor surface
Point(5) = {0.5416442642, 0.0955064977, 0, 1.0};   // magnet
Point(6) = {0.0955064977, 0.5416442642, 0, 1.0};   // magnet
Point(7) = {0.0000000000, 0.5500000000, 0, 1.0};
Point(8) = {0.6007327293, 0.1059253884, 0, 1.0};   // magnet top
Point(9) = {0.1059253884, 0.6007327293, 0, 1.0};   // magnet top
Point(10) = {0.6400000000, 0.0000000000, 0, 1.0};   // stator bore
Point(11) = {0.0000000000, 0.6400000000, 0, 1.0};
Point(12) = {1.0000000000, 0.0000000000, 0, 1.0};   // stator outside
Point(13) = {0.0000000000, 1.0000000000, 0, 1.0};
Point(14) = {0.6241911814, 0.1413696185, 0, 1.0};   // slot at 15 deg: neck on the bore
Point(15) = {0.6536407797, 0.1492606146, 0, 1.0};   // neck meets the tooth-tip shoe
Point(16) = {0.6406998275, 0.1975569059, 0, 1.0};
Point(17) = {0.6112502292, 0.1896659098, 0, 1.0};
Point(18) = {0.7086163747, 0.0600000000, 0, 1.0};   // shoe meets the tooth flank
Point(19) = {0.6436797821, 0.3023466631, 0, 1.0};
Point(20) = {0.7396746602, 0.0600000000, 0, 1.0};   // coil top (slot wedge line)
Point(21) = {0.6705770462, 0.3178758058, 0, 1.0};
Point(22) = {0.8949660872, 0.0600000000, 0, 1.0};   // slot bottom
Point(23) = {0.8050633671, 0.3955215194, 0, 1.0};
Point(24) = {0.4698806107, 0.4345252716, 0, 1.0};   // slot at 45 deg: neck on the bore
Point(25) = {0.4914392129, 0.4560838739, 0, 1.0};   // neck meets the tooth-tip shoe
Point(26) = {0.4560838739, 0.4914392129, 0, 1.0};
Point(27) = {0.4345252716, 0.4698806107, 0, 1.0};
Point(28) = {0.5836797821, 0.4062697116, 0, 1.0};   // shoe meets the tooth flank
Point(29) = {0.4062697116, 0.5836797821, 0, 1.0};
Point(30) = {0.6105770462, 0.4217988543, 0, 1.0};   // coil top (slot wedge line)
Point(31) = {0.4217988543, 0.6105770462, 0, 1.0};
Point(32) = {0.7450633671, 0.4994445678, 0, 1.0};   // slot bottom
Point(33) = {0.4994445678, 0.7450633671, 0, 1.0};
Point(34) = {0.1896659098, 0.6112502292, 0, 1.0};   // slot at 75 deg: neck on the bore
Point(35) = {0.1975569059, 0.6406998275, 0, 1.0};   // neck meets the tooth-tip shoe
Point(36) = {0.1492606146, 0.6536407797, 0, 1.0};
Point(37) = {0.1413696185, 0.6241911814, 0, 1.0};
Point(38) = {0.3023466631, 0.6436797821, 0, 1.0};   // shoe meets the tooth flank
Point(39) = {0.0600000000, 0.7086163747, 0, 1.0};
Point(40) = {0.3178758058, 0.6705770462, 0, 1.0};   // coil top (slot wedge line)
Point(41) = {0.0600000000, 0.7396746602, 0, 1.0};
Point(42) = {0.3955215194, 0.8050633671, 0, 1.0};   // slot bottom
Point(43) = {0.0600000000, 0.8949660872, 0, 1.0};

Line(1) = {1, 2};
Circle(2) = {2, 1, 3};
Line(3) = {3, 1};
Line(4) = {2, 4};
Circle(5) = {4, 1, 5};
Circle(6) = {5, 1, 6};
Circle(7) = {6, 1, 7};
Line(8) = {7, 3};
Line(9) = {5, 8};
Circle(10) = {8, 1, 9};
Line(11) = {9, 6};
Line(12) = {4, 10};
Circle(13) = {10, 1, 14};
Line(14) = {14, 15};
Line(15) = {15, 18};
Line(16) = {18, 20};
Line(17) = {20, 21};
Line(18) = {21, 19};
Line(19) = {19, 16};
Line(20) = {16, 17};
Circle(21) = {17, 1, 24};
Line(22) = {24, 25};
Line(23) = {25, 28};
Line(24) = {28, 30};
Line(25) = {30, 31};
Line(26) = {31, 29};
Line(27) = {29, 26};
Line(28) = {26, 27};
Circle(29) = {27, 1, 34};
Line(30) = {34, 35};
Line(31) = {35, 38};
Line(32) = {38, 40};
Line(33) = {40, 41};
Line(34) = {41, 39};
Line(35) = {39, 36};
Line(36) = {36, 37};
Circle(37) = {37, 1, 11};
Line(38) = {11, 7};
Line(39) = {10, 12};
Circle(40) = {12, 1, 13};
Line(41) = {13, 11};
Line(42) = {41, 43};
Line(43) = {43, 42};
Line(44) = {42, 40};
Line(45) = {31, 33};
Line(46) = {33, 32};
Line(47) = {32, 30};
Line(48) = {21, 23};
Line(49) = {23, 22};
Line(50) = {22, 20};

Curve Loop(1) = {4, 5, 6, 7, 8, -2};
Plane Surface(1) = {1};   // rotor iron
Curve Loop(2) = {9, 10, 11, -6};
Plane Surface(2) = {2};   // magnet
Curve Loop(3) = {12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, -7, -11, -10, -9, -5};
Plane Surface(3) = {3};   // air: gap, pole gaps, slot openings
Curve Loop(4) = {39, 40, 41, -37, -36, -35, -34, 42, 43, 44, -32, -31, -30, -29, -28, -27, -26, 45, 46, 47, -24, -23, -22, -21, -20, -19, -18, 48, 49, 50, -16, -15, -14, -13};
Plane Surface(4) = {4};   // stator iron
Curve Loop(5) = {-50, -49, -48, -17};
Plane Surface(5) = {5};   // coil, slot at 15 deg
Curve Loop(6) = {-47, -46, -45, -25};
Plane Surface(6) = {6};   // coil, slot at 45 deg
Curve Loop(7) = {-44, -43, -42, -33};
Plane Surface(7) = {7};   // coil, slot at 75 deg
Curve Loop(8) = {1, 2, 3};
Plane Surface(8) = {8};   // shaft

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5, 6, 7};
Physical Surface(6) = {8};
