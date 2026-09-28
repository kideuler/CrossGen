// An internal ply drop-off in a composite laminate. Six plies alternate
// between 0-degree (mat 1) and 90-degree (mat 2) fibre orientation; the
// middle pair (plies 3 and 4) is cut off at x = 1.2 to thin the laminate from
// six plies to four, the two plies above step down over the cut ends on a
// 25-degree ramp, and the triangular gap the ramp leaves behind the cut ends
// fills with neat resin (mat 3). The tool side (y = 0) stays flat. This is
// how every tapered composite part -- a wing skin, a blade root -- changes
// thickness, and the resin pocket is where it delaminates, which is why it is
// modelled in exactly this form.
//
// What this domain is for. Many long, parallel, thin layers that must all be
// followed at once, meeting a small acute feature. Every ply is 0.15 thick
// (5 elements at NP = 100, h = 0.03 here) and runs the length of the part,
// so a layout that is not ply-by-ply puts element edges across a ply, which
// is the one thing a laminate model cannot tolerate. Then, at the drop:
//
//  * The resin pocket is a triangle with corners of 90, 65 and 25 degrees:
//    the 25-degree tip, where the ramp lands on ply 2, is a sliver in a
//    material that exists nowhere else.
//  * Four junctions surround it. The cut ends give two orthogonal T's (ply
//    2 | ply 3 | resin and ply 3 | ply 4 | resin: directions 0, 90, 180 and
//    90, 180, 270, WELL-POSED); the top of the cut end and the ramp's foot
//    are ILL-POSED (directions {180, 270, 335} and {0, 155, 180}, i.e.
//    {0, 65} mod 90), and both carry an acute resin corner.
//  * The ramp puts 25-degree kinks into the ply 5 | ply 6 interface and into
//    the top surface, just above the 20-degree threshold at which
//    BoundaryFeatures calls a turn a corner.
//
// Same-orientation plies never touch: 0-degree plies 1, 3, 5 are one material
// on three surfaces and 90-degree plies 2, 4, 6 another, and dropping a 0/90
// pair keeps the thin end alternating (ply 5 lands on ply 2). The covering
// plies keep their thickness on the ramp (their interfaces are offset along
// the ramp's normal), so the kinks in successive interfaces are staggered.

// Tool side and the walls.
Point(1)  = {0, 0,    0, 1.0};
Point(2)  = {3, 0,    0, 1.0};
Point(3)  = {3, 0.15, 0, 1.0};
Point(4)  = {3, 0.30, 0, 1.0};
Point(5)  = {3, 0.45, 0, 1.0};
Point(6)  = {3, 0.60, 0, 1.0};
Point(7)  = {0, 0.90, 0, 1.0};
Point(8)  = {0, 0.75, 0, 1.0};
Point(9)  = {0, 0.60, 0, 1.0};
Point(10) = {0, 0.45, 0, 1.0};
Point(11) = {0, 0.30, 0, 1.0};
Point(12) = {0, 0.15, 0, 1.0};

// The drop: cut ends of plies 3 and 4 at x = 1.2, and the pocket's tip where
// the 25-degree ramp (run 0.3 / tan 25 = 0.643352) meets ply 2.
Point(21) = {1.2,      0.30, 0, 1.0};
Point(22) = {1.2,      0.45, 0, 1.0};
Point(23) = {1.2,      0.60, 0, 1.0};
Point(24) = {1.843352, 0.30, 0, 1.0};
// Kinks of the ply 5 | ply 6 interface and of the top surface: the ramp line
// offset by one and two ply thicknesses along its normal.
Point(25) = {1.233254, 0.75, 0, 1.0};
Point(26) = {1.876606, 0.45, 0, 1.0};
Point(27) = {1.266508, 0.90, 0, 1.0};
Point(28) = {1.909860, 0.60, 0, 1.0};

// Outer boundary, counter-clockwise.
Line(1)  = {1, 2};
Line(2)  = {2, 3};
Line(3)  = {3, 4};
Line(4)  = {4, 5};
Line(5)  = {5, 6};
Line(6)  = {6, 28};     // top surface, thin section
Line(7)  = {28, 27};    // top surface, ramp
Line(8)  = {27, 7};     // top surface, thick section
Line(9)  = {7, 8};
Line(10) = {8, 9};
Line(11) = {9, 10};
Line(12) = {10, 11};
Line(13) = {11, 12};
Line(14) = {12, 1};

// Ply interfaces, and the resin pocket's three sides.
Line(21) = {12, 3};     // ply 1 | ply 2
Line(22) = {11, 21};    // ply 2 | ply 3
Line(23) = {21, 24};    // ply 2 | resin
Line(24) = {24, 4};     // ply 2 | ply 5
Line(25) = {10, 22};    // ply 3 | ply 4
Line(26) = {9, 23};     // ply 4 | ply 5
Line(27) = {21, 22};    // cut end of ply 3 | resin
Line(28) = {22, 23};    // cut end of ply 4 | resin
Line(29) = {23, 24};    // ramp: ply 5 | resin
Line(30) = {8, 25};     // ply 5 | ply 6
Line(31) = {25, 26};
Line(32) = {26, 5};

Curve Loop(1) = {1, 2, -21, 14};
Plane Surface(1) = {1};                          // ply 1, 0 deg
Curve Loop(2) = {21, 3, -24, -23, -22, 13};
Plane Surface(2) = {2};                          // ply 2, 90 deg
Curve Loop(3) = {22, 27, -25, 12};
Plane Surface(3) = {3};                          // ply 3, 0 deg (dropped)
Curve Loop(4) = {25, 28, -26, 11};
Plane Surface(4) = {4};                          // ply 4, 90 deg (dropped)
Curve Loop(5) = {26, 29, 24, 4, -32, -31, -30, 10};
Plane Surface(5) = {5};                          // ply 5, 0 deg
Curve Loop(6) = {30, 31, 32, 5, 6, 7, 8, 9};
Plane Surface(6) = {6};                          // ply 6, 90 deg
Curve Loop(7) = {23, -29, -28, -27};
Plane Surface(7) = {7};                          // resin pocket

Physical Surface(1) = {1, 3, 5};
Physical Surface(2) = {2, 4, 6};
Physical Surface(3) = {7};
