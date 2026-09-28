// A honeycomb of three materials: regular hexagons of edge 0.25, coloured so
// that the three cells round every vertex all differ, cut by a 1.73 x 1.5
// box. The box's side walls run through the centres of alternate rows of
// cells (and along the vertical edges of the rows between), and its top and
// bottom through the centres of the first and last rows, so every cut cell
// is exactly half a hexagon and nothing is a sliver.
//
// What this domain is for. The 120-degree triple junction, many times over.
// When three phases meet with equal interfacial tension they meet at 120
// degrees (Plateau's rules for a foam, Herring's for an annealed
// polycrystal), so this is the equilibrium microstructure, not an adversarial
// one. Every interior vertex here is such a junction: directions 0, 120 and
// 240 up to rotation, which reduce mod 90 to three distinct values, so every
// one is ILL-POSED, identically. Where bubbles asks one inclusion question
// ten times, this asks one junction
// question twenty-eight times, and whether a method answers it the same way
// at every copy is directly readable.
//
// It also has no cheap answer. A hexagon is three quadrilaterals round a
// valence-3 centre, which puts two quads at three of its corners and one at
// the other three; round a junction the three cells then contribute 3 to 6
// quads, but each cell has three "two" corners for every two vertices it
// owns, so the average junction gets 4.5 and some must be irregular. The
// best layout is a global compromise, not a local template.
//
// Materials 1, 2 and 3 are one material each on many surfaces, and no two
// cells of one material touch, even at a point. The walls meet the cell
// edges at 60/120 degrees (side walls) and 90 degrees (top and bottom). At NP
// = 100 (h = 0.0173 here) a cell edge is 14 elements long.
//
// Coordinates are generated (exact to 1e-10) from the construction above;
// the cells are listed by centre.

// ---- generated geometry: R = 0.25, domain [0, 4 sqrt(3) R] x [0, 6 R] ----
Point(1) = {0.4330127019, 0.0000000000, 0, 1.0};
Point(2) = {0.4330127019, 0.1250000000, 0, 1.0};
Point(3) = {0.2165063509, 0.2500000000, 0, 1.0};
Point(4) = {-0.0000000000, 0.1250000000, 0, 1.0};
Point(5) = {-0.0000000000, 0.0000000000, 0, 1.0};
Point(6) = {0.8660254038, 0.0000000000, 0, 1.0};
Point(7) = {0.8660254038, 0.1250000000, 0, 1.0};
Point(8) = {0.6495190528, 0.2500000000, 0, 1.0};
Point(9) = {1.2990381057, 0.0000000000, 0, 1.0};
Point(10) = {1.2990381057, 0.1250000000, 0, 1.0};
Point(11) = {1.0825317547, 0.2500000000, 0, 1.0};
Point(12) = {1.7320508076, 0.0000000000, 0, 1.0};
Point(13) = {1.7320508076, 0.1250000000, 0, 1.0};
Point(14) = {1.5155444566, 0.2500000000, 0, 1.0};
Point(15) = {0.2165063509, 0.5000000000, 0, 1.0};
Point(16) = {0.0000000000, 0.6250000000, 0, 1.0};
Point(17) = {0.6495190528, 0.5000000000, 0, 1.0};
Point(18) = {0.4330127019, 0.6250000000, 0, 1.0};
Point(19) = {1.0825317547, 0.5000000000, 0, 1.0};
Point(20) = {0.8660254038, 0.6250000000, 0, 1.0};
Point(21) = {1.5155444566, 0.5000000000, 0, 1.0};
Point(22) = {1.2990381057, 0.6250000000, 0, 1.0};
Point(23) = {1.7320508076, 0.6250000000, 0, 1.0};
Point(24) = {0.4330127019, 0.8750000000, 0, 1.0};
Point(25) = {0.2165063509, 1.0000000000, 0, 1.0};
Point(26) = {-0.0000000000, 0.8750000000, 0, 1.0};
Point(27) = {0.8660254038, 0.8750000000, 0, 1.0};
Point(28) = {0.6495190528, 1.0000000000, 0, 1.0};
Point(29) = {1.2990381057, 0.8750000000, 0, 1.0};
Point(30) = {1.0825317547, 1.0000000000, 0, 1.0};
Point(31) = {1.7320508076, 0.8750000000, 0, 1.0};
Point(32) = {1.5155444566, 1.0000000000, 0, 1.0};
Point(33) = {0.2165063509, 1.2500000000, 0, 1.0};
Point(34) = {0.0000000000, 1.3750000000, 0, 1.0};
Point(35) = {0.6495190528, 1.2500000000, 0, 1.0};
Point(36) = {0.4330127019, 1.3750000000, 0, 1.0};
Point(37) = {1.0825317547, 1.2500000000, 0, 1.0};
Point(38) = {0.8660254038, 1.3750000000, 0, 1.0};
Point(39) = {1.5155444566, 1.2500000000, 0, 1.0};
Point(40) = {1.2990381057, 1.3750000000, 0, 1.0};
Point(41) = {1.7320508076, 1.3750000000, 0, 1.0};
Point(42) = {0.4330127019, 1.5000000000, 0, 1.0};
Point(43) = {-0.0000000000, 1.5000000000, 0, 1.0};
Point(44) = {0.8660254038, 1.5000000000, 0, 1.0};
Point(45) = {1.2990381057, 1.5000000000, 0, 1.0};
Point(46) = {1.7320508076, 1.5000000000, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 1};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 2};
Line(9) = {1, 6};
Line(10) = {9, 10};
Line(11) = {10, 11};
Line(12) = {11, 7};
Line(13) = {6, 9};
Line(14) = {12, 13};
Line(15) = {13, 14};
Line(16) = {14, 10};
Line(17) = {9, 12};
Line(18) = {15, 16};
Line(19) = {16, 4};
Line(20) = {3, 15};
Line(21) = {17, 18};
Line(22) = {18, 15};
Line(23) = {8, 17};
Line(24) = {19, 20};
Line(25) = {20, 17};
Line(26) = {11, 19};
Line(27) = {21, 22};
Line(28) = {22, 19};
Line(29) = {14, 21};
Line(30) = {23, 21};
Line(31) = {13, 23};
Line(32) = {24, 25};
Line(33) = {25, 26};
Line(34) = {26, 16};
Line(35) = {18, 24};
Line(36) = {27, 28};
Line(37) = {28, 24};
Line(38) = {20, 27};
Line(39) = {29, 30};
Line(40) = {30, 27};
Line(41) = {22, 29};
Line(42) = {31, 32};
Line(43) = {32, 29};
Line(44) = {23, 31};
Line(45) = {33, 34};
Line(46) = {34, 26};
Line(47) = {25, 33};
Line(48) = {35, 36};
Line(49) = {36, 33};
Line(50) = {28, 35};
Line(51) = {37, 38};
Line(52) = {38, 35};
Line(53) = {30, 37};
Line(54) = {39, 40};
Line(55) = {40, 37};
Line(56) = {32, 39};
Line(57) = {41, 39};
Line(58) = {31, 41};
Line(59) = {42, 43};
Line(60) = {43, 34};
Line(61) = {36, 42};
Line(62) = {44, 42};
Line(63) = {38, 44};
Line(64) = {45, 44};
Line(65) = {40, 45};
Line(66) = {46, 45};
Line(67) = {41, 46};

Curve Loop(1) = {1, 2, 3, 4, 5};   // cell at (0.217, 0.000)
Plane Surface(1) = {1};
Curve Loop(2) = {6, 7, 8, -1, 9};   // cell at (0.650, 0.000)
Plane Surface(2) = {2};
Curve Loop(3) = {10, 11, 12, -6, 13};   // cell at (1.083, 0.000)
Plane Surface(3) = {3};
Curve Loop(4) = {14, 15, 16, -10, 17};   // cell at (1.516, 0.000)
Plane Surface(4) = {4};
Curve Loop(5) = {18, 19, -3, 20};   // cell at (0.000, 0.375)
Plane Surface(5) = {5};
Curve Loop(6) = {21, 22, -20, -2, -8, 23};   // cell at (0.433, 0.375)
Plane Surface(6) = {6};
Curve Loop(7) = {24, 25, -23, -7, -12, 26};   // cell at (0.866, 0.375)
Plane Surface(7) = {7};
Curve Loop(8) = {27, 28, -26, -11, -16, 29};   // cell at (1.299, 0.375)
Plane Surface(8) = {8};
Curve Loop(9) = {30, -29, -15, 31};   // cell at (1.732, 0.375)
Plane Surface(9) = {9};
Curve Loop(10) = {32, 33, 34, -18, -22, 35};   // cell at (0.217, 0.750)
Plane Surface(10) = {10};
Curve Loop(11) = {36, 37, -35, -21, -25, 38};   // cell at (0.650, 0.750)
Plane Surface(11) = {11};
Curve Loop(12) = {39, 40, -38, -24, -28, 41};   // cell at (1.083, 0.750)
Plane Surface(12) = {12};
Curve Loop(13) = {42, 43, -41, -27, -30, 44};   // cell at (1.516, 0.750)
Plane Surface(13) = {13};
Curve Loop(14) = {45, 46, -33, 47};   // cell at (0.000, 1.125)
Plane Surface(14) = {14};
Curve Loop(15) = {48, 49, -47, -32, -37, 50};   // cell at (0.433, 1.125)
Plane Surface(15) = {15};
Curve Loop(16) = {51, 52, -50, -36, -40, 53};   // cell at (0.866, 1.125)
Plane Surface(16) = {16};
Curve Loop(17) = {54, 55, -53, -39, -43, 56};   // cell at (1.299, 1.125)
Plane Surface(17) = {17};
Curve Loop(18) = {57, -56, -42, 58};   // cell at (1.732, 1.125)
Plane Surface(18) = {18};
Curve Loop(19) = {59, 60, -45, -49, 61};   // cell at (0.217, 1.500)
Plane Surface(19) = {19};
Curve Loop(20) = {62, -61, -48, -52, 63};   // cell at (0.650, 1.500)
Plane Surface(20) = {20};
Curve Loop(21) = {64, -63, -51, -55, 65};   // cell at (1.083, 1.500)
Plane Surface(21) = {21};
Curve Loop(22) = {66, -65, -54, -57, 67};   // cell at (1.516, 1.500)
Plane Surface(22) = {22};

Physical Surface(1) = {3, 6, 9, 12, 15, 18, 21};
Physical Surface(2) = {1, 4, 7, 10, 13, 16, 19, 22};
Physical Surface(3) = {2, 5, 8, 11, 14, 17, 20};
