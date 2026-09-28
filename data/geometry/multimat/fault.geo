// A normal fault through layered sediment, with a throw larger than a bed.
// Five horizontal units (mat 1 at the bottom to mat 5 at the top, each 0.40
// thick) are cut by a fault dipping 60 degrees to the right, and the hanging
// wall on the right has dropped 0.60 -- one and a half beds -- so that across
// the fault each unit faces the one or two below it rather than itself. This
// is the geometry of a fault-seal ("juxtaposition" or Allan diagram) study,
// where what matters is exactly which unit sits against which.
//
// What this domain is for. The fault is a straight interface carrying a
// ladder of SIX interior triple junctions, alternately on its two sides:
// every bed boundary that reaches it ends there, and because the throw is
// not a multiple of the bed thickness, no bed boundary on one side lines up
// with one on the other. Each junction has interface directions 0 (the bed
// boundary), 60 and 240 (the fault): reduced mod 90 that is {0, 60}, so all
// six are ILL-POSED, the same way, at an angle geology rather than an
// adversary chose, and close enough along one straight line that the
// singularities they force have to be placed in relation to each other.
//
// Two further things fall out of the physics and are deliberate:
//
//  * A fault is an interface only where it juxtaposes DIFFERENT units. Above
//    y = 1.65 unit 5 faces unit 5, so the mesh cannot see the fault there at
//    all, and the fault's visible trace ends at (1.447, 1.65) in a kink of
//    the 4|5 interface rather than in a junction. The top unit is one surface
//    that wraps over the end of the fault.
//
//  * Units 2, 3 and 4 are each one material on two disconnected surfaces
//    (one per side of the fault), since the throw separates them completely.
//
// The fault meets the bottom wall at 60/120 degrees; every bed boundary
// meets a side wall at a right angle. At NP = 100 (h = 0.04 here) a bed is 10
// elements deep, and the shortest stretch of fault between two junctions
// (0.23 long) is about 6 elements.

Point(1) = {0, 0, 0, 1.0};
Point(2) = {4, 0, 0, 1.0};
Point(3) = {4, 2, 0, 1.0};
Point(4) = {0, 2, 0, 1.0};

// The fault, x = 2.4 - y / tan(60), at every height where something meets it.
Point(10) = {2.400000, 0.00, 0, 1.0};   // bottom wall
Point(11) = {2.255662, 0.25, 0, 1.0};   // hanging wall 2|3
Point(12) = {2.140192, 0.45, 0, 1.0};   // footwall 1|2
Point(13) = {2.024722, 0.65, 0, 1.0};   // hanging wall 3|4
Point(14) = {1.909252, 0.85, 0, 1.0};   // footwall 2|3
Point(15) = {1.793782, 1.05, 0, 1.0};   // hanging wall 4|5
Point(16) = {1.678312, 1.25, 0, 1.0};   // footwall 3|4
Point(17) = {1.447372, 1.65, 0, 1.0};   // footwall 4|5: the fault's visible tip

// Bed boundaries on the walls: footwall on the left, hanging wall (dropped
// 0.60) on the right.
Point(21) = {0, 0.45, 0, 1.0};
Point(22) = {0, 0.85, 0, 1.0};
Point(23) = {0, 1.25, 0, 1.0};
Point(24) = {0, 1.65, 0, 1.0};
Point(31) = {4, 0.25, 0, 1.0};
Point(32) = {4, 0.65, 0, 1.0};
Point(33) = {4, 1.05, 0, 1.0};

// Outer boundary, counter-clockwise.
Line(1)  = {1, 10};
Line(2)  = {10, 2};
Line(3)  = {2, 31};
Line(4)  = {31, 32};
Line(5)  = {32, 33};
Line(6)  = {33, 3};
Line(7)  = {3, 4};
Line(8)  = {4, 24};
Line(9)  = {24, 23};
Line(10) = {23, 22};
Line(11) = {22, 21};
Line(12) = {21, 1};

// The fault, bottom to top.
Line(20) = {10, 11};
Line(21) = {11, 12};
Line(22) = {12, 13};
Line(23) = {13, 14};
Line(24) = {14, 15};
Line(25) = {15, 16};
Line(26) = {16, 17};

// Bed boundaries: footwall (wall to fault) and hanging wall (fault to wall).
Line(30) = {21, 12};
Line(31) = {22, 14};
Line(32) = {23, 16};
Line(33) = {24, 17};
Line(40) = {11, 31};
Line(41) = {13, 32};
Line(42) = {15, 33};

Curve Loop(1) = {1, 20, 21, -30, 12};          // unit 1, footwall
Plane Surface(1) = {1};
Curve Loop(2) = {30, 22, 23, -31, 11};         // unit 2, footwall
Plane Surface(2) = {2};
Curve Loop(3) = {31, 24, 25, -32, 10};         // unit 3, footwall
Plane Surface(3) = {3};
Curve Loop(4) = {32, 26, -33, 9};              // unit 4, footwall
Plane Surface(4) = {4};
Curve Loop(12) = {2, 3, -40, -20};             // unit 2, hanging wall
Plane Surface(12) = {12};
Curve Loop(13) = {40, 4, -41, -22, -21};       // unit 3, hanging wall
Plane Surface(13) = {13};
Curve Loop(14) = {41, 5, -42, -24, -23};       // unit 4, hanging wall
Plane Surface(14) = {14};
Curve Loop(5) = {42, 6, 7, 8, 33, -26, -25};   // unit 5, both walls
Plane Surface(5) = {5};

Physical Surface(1) = {1};
Physical Surface(2) = {2, 12};
Physical Surface(3) = {3, 13};
Physical Surface(4) = {4, 14};
Physical Surface(5) = {5};
