// A planar MOSFET in cross-section, as a TCAD device simulator meshes it.
// Eight materials: the silicon substrate (mat 1); n+ source and drain wells
// (mat 2) whose junction with the substrate curves round under the gate
// edge; shallow-trench-isolation oxide (mat 3) with tapered walls at both
// ends; the gate oxide (mat 4); the polysilicon gate (mat 5); nitride
// sidewall spacers (mat 6) with the usual rounded outer face; the pre-metal
// dielectric (mat 7) over everything; and two tungsten contact plugs (mat 8)
// through it onto the source and drain.
//
// What this domain is for. Density of small features against a large one,
// which is the defining property of device meshes. The substrate is 2 wide
// and 0.8 deep; the gate oxide under it
// is 0.05 thick (2.5 elements at NP = 100, h = 0.02 here), and the strip of
// substrate between a well's edge and the gate's edge is as narrow. Around
// that there are eighteen interior triple junctions in 0.6 x 0.5, and eight
// more where interfaces meet the walls and the top. Almost all are
// orthogonal -- the spacers' faces and the wells' junction curves meet the
// silicon surface at right angles, as real spacers and diffusion fronts do --
// so the domain is mostly a well-posed control; the exceptions are the four
// where a well and the dielectric meet the trench walls, which lean 8.13
// degrees off vertical (directions {0, 81.87} mod 90, ILL-POSED by 8
// degrees). An orthogonal junction a layout gets wrong here costs a block
// that a region fifty times larger has to absorb.
//
// The pre-metal dielectric is one material on three surfaces, cut apart by
// the contact plugs, which run from the wells to the top of the domain. The
// device is mirror-symmetric about x = 0, so a method's layout should be too.

Point(1) = {-1, -0.8,  0, 1.0};
Point(2) = { 1, -0.8,  0, 1.0};
Point(3) = { 1,  0.45, 0, 1.0};
Point(4) = {-1,  0.45, 0, 1.0};
Point(5) = {-1, -0.35, 0, 1.0};   // trench bottom on the walls
Point(6) = {-1,  0,    0, 1.0};   // silicon surface on the walls
Point(7) = { 1, -0.35, 0, 1.0};
Point(8) = { 1,  0,    0, 1.0};
Point(9)  = {-0.55, 0.45, 0, 1.0};   // contact plugs at the top
Point(10) = {-0.40, 0.45, 0, 1.0};
Point(11) = { 0.40, 0.45, 0, 1.0};
Point(12) = { 0.55, 0.45, 0, 1.0};

// The silicon surface, y = 0, left to right.
Point(21) = {-0.75, 0, 0, 1.0};   // trench edge
Point(22) = {-0.55, 0, 0, 1.0};   // contact
Point(23) = {-0.40, 0, 0, 1.0};   // contact
Point(24) = {-0.24, 0, 0, 1.0};   // spacer foot
Point(25) = {-0.17, 0, 0, 1.0};   // edge of the source well
Point(26) = {-0.12, 0, 0, 1.0};   // gate edge
Point(27) = { 0.12, 0, 0, 1.0};
Point(28) = { 0.17, 0, 0, 1.0};
Point(29) = { 0.24, 0, 0, 1.0};
Point(30) = { 0.40, 0, 0, 1.0};
Point(31) = { 0.55, 0, 0, 1.0};
Point(32) = { 0.75, 0, 0, 1.0};

// Trench walls, x = +-(0.75 - y / 7), and where the wells' floors meet them.
Point(41) = {-0.7714285714, -0.15, 0, 1.0};
Point(42) = {-0.80,         -0.35, 0, 1.0};
Point(43) = { 0.7714285714, -0.15, 0, 1.0};
Point(44) = { 0.80,         -0.35, 0, 1.0};

// Wells: floor at y = -0.15, and a quarter ellipse (semi-axes 0.10 across,
// 0.15 down) from the floor up to the surface, meeting it at a right angle.
Point(45) = {-0.27, -0.15, 0, 1.0};
Point(46) = { 0.27, -0.15, 0, 1.0};
Point(47) = {-0.27,  0,    0, 1.0};   // ellipse centres
Point(48) = { 0.27,  0,    0, 1.0};

// Gate stack: oxide to y = 0.05, polysilicon to 0.35.
Point(51) = {-0.12, 0.05, 0, 1.0};
Point(52) = {-0.12, 0.35, 0, 1.0};
Point(53) = { 0.12, 0.05, 0, 1.0};
Point(54) = { 0.12, 0.35, 0, 1.0};

// The silicon surface.
Line(101) = {6, 21};
Line(102) = {21, 22};
Line(103) = {22, 23};
Line(104) = {23, 24};
Line(105) = {24, 25};
Line(106) = {25, 26};
Line(107) = {26, 27};   // the channel
Line(108) = {27, 28};
Line(109) = {28, 29};
Line(110) = {29, 30};
Line(111) = {30, 31};
Line(112) = {31, 32};
Line(113) = {32, 8};

// Trench walls and floors.
Line(121) = {21, 41};
Line(122) = {41, 42};
Line(123) = {42, 5};
Line(124) = {32, 43};
Line(125) = {43, 44};
Line(126) = {44, 7};

// Wells: floors and junction curves.
Line(131) = {41, 45};
Line(132) = {46, 43};
Ellipse(133) = {25, 47, 45, 45};
Ellipse(134) = {28, 48, 46, 46};

// Gate stack.
Line(141) = {26, 51};
Line(142) = {51, 52};
Line(143) = {51, 53};
Line(144) = {27, 53};
Line(145) = {53, 54};
Line(146) = {52, 54};

// Spacers' outer faces: quarter ellipses about the gate's lower corners,
// vertical at the foot and horizontal at the top of the gate.
Ellipse(151) = {24, 26, 52, 52};
Ellipse(152) = {29, 27, 54, 54};

// Contact plugs' sides.
Line(161) = {22, 9};
Line(162) = {23, 10};
Line(163) = {30, 11};
Line(164) = {31, 12};

// Top and walls.
Line(171) = {3, 12};
Line(172) = {12, 11};
Line(173) = {11, 10};
Line(174) = {10, 9};
Line(175) = {9, 4};
Line(181) = {4, 6};
Line(182) = {6, 5};
Line(183) = {5, 1};
Line(184) = {1, 2};
Line(185) = {2, 7};
Line(186) = {7, 8};
Line(187) = {8, 3};

Curve Loop(1) = {184, 185, -126, -125, -132, -134, -108, -107, -106, 133, -131, 122, 123, 183};
Plane Surface(1) = {1};                                     // substrate
Curve Loop(2) = {131, -133, -105, -104, -103, -102, 121};
Plane Surface(2) = {2};                                     // source
Curve Loop(3) = {132, -124, -112, -111, -110, -109, 134};
Plane Surface(3) = {3};                                     // drain
Curve Loop(4) = {-123, -122, -121, -101, 182};
Plane Surface(4) = {4};                                     // trench, left
Curve Loop(5) = {126, 186, -113, 124, 125};
Plane Surface(5) = {5};                                     // trench, right
Curve Loop(6) = {107, 144, -143, -141};
Plane Surface(6) = {6};                                     // gate oxide
Curve Loop(7) = {143, 145, -146, -142};
Plane Surface(7) = {7};                                     // gate
Curve Loop(8) = {105, 106, 141, 142, -151};
Plane Surface(8) = {8};                                     // spacer, left
Curve Loop(9) = {108, 109, 152, -145, -144};
Plane Surface(9) = {9};                                     // spacer, right
Curve Loop(10) = {103, 162, 174, -161};
Plane Surface(10) = {10};                                   // contact, source
Curve Loop(11) = {111, 164, 172, -163};
Plane Surface(11) = {11};                                   // contact, drain
Curve Loop(12) = {101, 102, 161, 175, 181};
Plane Surface(12) = {12};                                   // dielectric, left
Curve Loop(13) = {112, 113, 187, 171, -164};
Plane Surface(13) = {13};                                   // dielectric, right
Curve Loop(14) = {104, 151, 146, -152, 110, 163, 173, -162};
Plane Surface(14) = {14};                                   // dielectric over the gate

Physical Surface(1) = {1};
Physical Surface(2) = {2, 3};
Physical Surface(3) = {4, 5};
Physical Surface(4) = {6};
Physical Surface(5) = {7};
Physical Surface(6) = {8, 9};
Physical Surface(7) = {12, 13, 14};
Physical Surface(8) = {10, 11};
