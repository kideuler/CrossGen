// A welded cruciform joint: the standard fatigue specimen of welded-steel
// design (the IIW recommendations test weld details on exactly this shape). A
// through plate (mat 1) carries two attachment plates (mat 2) on its faces,
// each held by two fillet welds (mat 3) with 0.3 legs and flat 45-degree
// faces. Plates are 0.3 thick; the through plate is 3.0 long and the
// attachments stand 1.05 proud of it.
//
// What this domain is for. The angle a rounding rule finds hardest, in a
// geometry nobody would call contrived. At each of the eight weld toes the
// weld's own corner is 45 degrees -- half a quarter turn, the exact tie in
// round(angle / 90) -- and the domain's boundary turns there by 45 degrees
// the other way, so the domain's corner is 225 degrees, 2.5 quarter turns,
// the same tie. Whatever a method decides at one toe it has to decide at all
// eight, and by symmetry the eight decisions should agree; if they do not,
// the method is breaking a tie with noise. This is the rounding question with
// no right answer to fall back on.
//
// The four weld roots are orthogonal triple junctions (through plate,
// attachment, weld; interface directions 0, 180 and 90 or 270), congruent mod
// 90 and so WELL-POSED: the difficulty is entirely in the toes, which are
// boundary junctions.
//
// At NP = 100 (h = 0.03 here) the plates are 10 elements thick and the weld
// legs 10 elements long.

// Outer boundary, counter-clockwise from the through plate's lower-left
// corner; the Qs on a weld face are its two toes.
Point(1)  = {-1.50, -0.15, 0, 1.0};
Point(2)  = {-0.45, -0.15, 0, 1.0};   // toe
Point(3)  = {-0.15, -0.45, 0, 1.0};   // toe
Point(4)  = {-0.15, -1.20, 0, 1.0};
Point(5)  = { 0.15, -1.20, 0, 1.0};
Point(6)  = { 0.15, -0.45, 0, 1.0};   // toe
Point(7)  = { 0.45, -0.15, 0, 1.0};   // toe
Point(8)  = { 1.50, -0.15, 0, 1.0};
Point(9)  = { 1.50,  0.15, 0, 1.0};
Point(10) = { 0.45,  0.15, 0, 1.0};   // toe
Point(11) = { 0.15,  0.45, 0, 1.0};   // toe
Point(12) = { 0.15,  1.20, 0, 1.0};
Point(13) = {-0.15,  1.20, 0, 1.0};
Point(14) = {-0.15,  0.45, 0, 1.0};   // toe
Point(15) = {-0.45,  0.15, 0, 1.0};   // toe
Point(16) = {-1.50,  0.15, 0, 1.0};

// Weld roots: where an attachment's corner sits on the through plate.
Point(21) = {-0.15, -0.15, 0, 1.0};
Point(22) = { 0.15, -0.15, 0, 1.0};
Point(23) = { 0.15,  0.15, 0, 1.0};
Point(24) = {-0.15,  0.15, 0, 1.0};

// Outer boundary; Lines 2, 6, 10 and 14 are the weld faces.
For k In {1:15}
  Line(k) = {k, k + 1};
EndFor
Line(16) = {16, 1};

// Interfaces. Through plate | weld, through plate | attachment, and
// attachment | weld, each written once.
Line(21) = {2, 21};
Line(22) = {21, 22};
Line(23) = {22, 7};
Line(24) = {10, 23};
Line(25) = {23, 24};
Line(26) = {24, 15};
Line(27) = {3, 21};
Line(28) = {22, 6};
Line(29) = {23, 11};
Line(30) = {14, 24};

Curve Loop(1) = {1, 21, 22, 23, 7, 8, 9, 24, 25, 26, 15, 16};
Plane Surface(1) = {1};                 // through plate
Curve Loop(2) = {-27, 3, 4, 5, -28, -22};
Plane Surface(2) = {2};                 // lower attachment
Curve Loop(3) = {29, 11, 12, 13, 30, -25};
Plane Surface(3) = {3};                 // upper attachment
Curve Loop(4) = {2, 27, -21};
Plane Surface(4) = {4};                 // weld, lower left
Curve Loop(5) = {6, -23, 28};
Plane Surface(5) = {5};                 // weld, lower right
Curve Loop(6) = {10, -29, -24};
Plane Surface(6) = {6};                 // weld, upper right
Curve Loop(7) = {14, -26, -30};
Plane Surface(7) = {7};                 // weld, upper left

Physical Surface(1) = {1};
Physical Surface(2) = {2, 3};
Physical Surface(3) = {4, 5, 6, 7};
