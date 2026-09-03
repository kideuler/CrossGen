// Flip-chip package cross-section: a laminate substrate (mat 1), an array of
// five solder bumps (mat 3) standing in a bed of underfill (mat 2), carrying a
// silicon die (mat 4). Four materials, twenty interior triple junctions and a
// stepped free boundary, all of it orthogonal.
//
// This is the benign counterpart to the weld (geom013), and it is the reason to
// have it: every junction here is a place where a horizontal interface meets a
// vertical one, so the incident directions are 0 and 90 and are congruent mod
// 90. A single cross is tangent to both, the junction has a well-posed
// vertex value, and converting the face-based field to vertices costs nothing.
// The synthetic controls geom003 and geom005 make the same point with one
// junction each; this one makes it twenty times over on a layout an engineer
// would actually hand you, which is what a quantitative junction table needs if
// the benign column is not to rest on two rows.
//
// It is also the thin-feature case at NP=70: the bump pitch leaves the
// underfill 3.2 elements wide at its narrowest, the bumps are 3.5 elements
// across and the die-to-substrate gap is 5.8 elements deep, so the block
// structure R3 asks for has to survive in regions only a few cells thick.
//
// The underfill is one material carried by six disconnected surfaces -- in a
// section plane the bumps cut it into pieces, though it is connected in the
// part -- so Physical Surface(2) names all six. Two modelling simplifications:
// the bumps are drawn as rectangles rather than the usual barrel, and the
// underfill fillet that would bulge past the die edge is omitted. Both would
// replace an exactly orthogonal junction with a nearly orthogonal one and cost
// the domain the property it is here to supply.

Point(1) = {-1.2000, -0.4500, 0, 1.0};
Point(2) = {1.2000, -0.4500, 0, 1.0};
Point(3) = {1.2000, -0.1500, 0, 1.0};
Point(4) = {0.8500, -0.1500, 0, 1.0};
Point(5) = {-1.2000, -0.1500, 0, 1.0};
Point(6) = {-0.8500, -0.1500, 0, 1.0};
Point(7) = {0.8500, 0.0500, 0, 1.0};
Point(8) = {0.8500, 0.4000, 0, 1.0};
Point(9) = {-0.8500, 0.4000, 0, 1.0};
Point(10) = {-0.8500, 0.0500, 0, 1.0};
Point(11) = {-0.7400, -0.1500, 0, 1.0};
Point(12) = {-0.7400, 0.0500, 0, 1.0};
Point(13) = {-0.6200, -0.1500, 0, 1.0};
Point(14) = {-0.6200, 0.0500, 0, 1.0};
Point(15) = {-0.4000, -0.1500, 0, 1.0};
Point(16) = {-0.4000, 0.0500, 0, 1.0};
Point(17) = {-0.2800, -0.1500, 0, 1.0};
Point(18) = {-0.2800, 0.0500, 0, 1.0};
Point(19) = {-0.0600, -0.1500, 0, 1.0};
Point(20) = {-0.0600, 0.0500, 0, 1.0};
Point(21) = {0.0600, -0.1500, 0, 1.0};
Point(22) = {0.0600, 0.0500, 0, 1.0};
Point(23) = {0.2800, -0.1500, 0, 1.0};
Point(24) = {0.2800, 0.0500, 0, 1.0};
Point(25) = {0.4000, -0.1500, 0, 1.0};
Point(26) = {0.4000, 0.0500, 0, 1.0};
Point(27) = {0.6200, -0.1500, 0, 1.0};
Point(28) = {0.6200, 0.0500, 0, 1.0};
Point(29) = {0.7400, -0.1500, 0, 1.0};
Point(30) = {0.7400, 0.0500, 0, 1.0};

// Outer boundary.
Line(1) = {1, 2};            // package bottom
Line(2) = {2, 3};            // substrate right wall
Line(3) = {4, 3};            // substrate top, exposed right of the die
Line(4) = {5, 6};            // substrate top, exposed left of the die
Line(5) = {1, 5};            // substrate left wall
Line(6) = {4, 7};            // underfill free surface, right
Line(7) = {7, 8};            // die side wall, right
Line(8) = {9, 8};            // die top
Line(9) = {6, 10};           // underfill free surface, left
Line(10) = {10, 9};           // die side wall, left

// Material interfaces. Each is used forward in one Curve Loop and reversed in
// the other, so the two surfaces share its nodes and the interface is conformal.
Line(11) = {6, 11};           // substrate | underfill
Line(12) = {10, 12};          // underfill | die
Line(13) = {11, 13};          // substrate | bump
Line(14) = {12, 14};          // bump | die
Line(15) = {13, 15};          // substrate | underfill
Line(16) = {14, 16};          // underfill | die
Line(17) = {15, 17};          // substrate | bump
Line(18) = {16, 18};          // bump | die
Line(19) = {17, 19};          // substrate | underfill
Line(20) = {18, 20};          // underfill | die
Line(21) = {19, 21};          // substrate | bump
Line(22) = {20, 22};          // bump | die
Line(23) = {21, 23};          // substrate | underfill
Line(24) = {22, 24};          // underfill | die
Line(25) = {23, 25};          // substrate | bump
Line(26) = {24, 26};          // bump | die
Line(27) = {25, 27};          // substrate | underfill
Line(28) = {26, 28};          // underfill | die
Line(29) = {27, 29};          // substrate | bump
Line(30) = {28, 30};          // bump | die
Line(31) = {29, 4};           // substrate | underfill
Line(32) = {30, 7};           // underfill | die
Line(33) = {11, 12};          // bump 1, left wall
Line(34) = {13, 14};          // bump 1, right wall
Line(35) = {15, 16};          // bump 2, left wall
Line(36) = {17, 18};          // bump 2, right wall
Line(37) = {19, 20};          // bump 3, left wall
Line(38) = {21, 22};          // bump 3, right wall
Line(39) = {23, 24};          // bump 4, left wall
Line(40) = {25, 26};          // bump 4, right wall
Line(41) = {27, 28};          // bump 5, left wall
Line(42) = {29, 30};          // bump 5, right wall

Curve Loop(1) = {1, 2, -3, -31, -29, -27, -25, -23, -21, -19, -17, -15, -13, -11, -4, -5};
Plane Surface(1) = {1};   // laminate substrate
Curve Loop(2) = {12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 7, -8, -10};
Plane Surface(2) = {2};   // silicon die
Curve Loop(3) = {13, 34, -14, -33};
Plane Surface(3) = {3};   // solder bump 1
Curve Loop(4) = {17, 36, -18, -35};
Plane Surface(4) = {4};   // solder bump 2
Curve Loop(5) = {21, 38, -22, -37};
Plane Surface(5) = {5};   // solder bump 3
Curve Loop(6) = {25, 40, -26, -39};
Plane Surface(6) = {6};   // solder bump 4
Curve Loop(7) = {29, 42, -30, -41};
Plane Surface(7) = {7};   // solder bump 5
Curve Loop(8) = {11, 33, -12, -9};
Plane Surface(8) = {8};   // underfill 1
Curve Loop(9) = {15, 35, -16, -34};
Plane Surface(9) = {9};   // underfill 2
Curve Loop(10) = {19, 37, -20, -36};
Plane Surface(10) = {10};   // underfill 3
Curve Loop(11) = {23, 39, -24, -38};
Plane Surface(11) = {11};   // underfill 4
Curve Loop(12) = {27, 41, -28, -40};
Plane Surface(12) = {12};   // underfill 5
Curve Loop(13) = {31, 6, -32, -42};
Plane Surface(13) = {13};   // underfill 6

Physical Surface(1) = {1};                             // laminate substrate
Physical Surface(2) = {8, 9, 10, 11, 12, 13};          // underfill
Physical Surface(3) = {3, 4, 5, 6, 7};                 // solder bumps
Physical Surface(4) = {2};                             // silicon die
