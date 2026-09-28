// A zoned embankment dam on a rock foundation, in cross-section: the model
// seepage and slope-stability analyses are run on. A clay core (mat 1) with
// 1:0.25 faces rises from a cutoff trench keyed into the foundation (mat 4),
// and a grout curtain (mat 5) continues the seal down to the base of the
// section. Sand filters (mat 2) line both faces of the core; the downstream
// one is a chimney drain that turns at its foot into a blanket drain running
// along the foundation to the downstream toe. Rockfill shells (mat 3) with
// 1:1.5 faces make up the rest of the embankment and join over the top of the
// core, so the crest is shell.
//
// What this domain is for. Thin bands at design slopes, meeting each other
// and the foundation. The filters are 3.4 elements thick at NP = 100 (h =
// 0.046 here), inclined at 76 degrees; the downstream filter is one L-shaped
// region -- chimney plus blanket -- with a re-entrant corner where the two
// meet; the blanket, the crest cap and the curtain are 3.3, 3.3 and 2.6
// elements thick. Seven interior junctions. Where the upstream filter's outer
// face stands on the foundation, and at the two corners of the core's crest,
// the directions are {0, 76} or {0, 14} mod 90; at the core's two feet the
// trench walls add a third, {0, 33, 76} upstream and {0, 14, 57}
// downstream. Those five are ILL-POSED, at angles a designer set. The
// curtain's two T's on the trench floor are orthogonal and WELL-POSED. The
// foundation is one
// material on two surfaces, split by the trench and the curtain, and the
// shell is one surface that wraps over the core.
//
// The upstream toe is a 33.7-degree wedge of shell on dS, and the blanket
// drain ends on the downstream face, so both toes are boundary junctions.

// Foundation block and the crest.
Point(1)  = {-2.3,   -0.5,  0, 1.0};
Point(2)  = { 2.3,   -0.5,  0, 1.0};
Point(3)  = { 2.3,    0,    0, 1.0};
Point(4)  = {-2.3,    0,    0, 1.0};
Point(5)  = {-1.95,   0,    0, 1.0};   // upstream toe
Point(6)  = { 1.95,   0,    0, 1.0};   // downstream toe
Point(7)  = {-0.15,   1.2,  0, 1.0};   // crest
Point(8)  = { 0.15,   1.2,  0, 1.0};

// Core: faces x = +-(0.1 + 0.25 (1.05 - y)); cutoff trench to y = -0.25.
Point(11) = {-0.3625, 0,    0, 1.0};
Point(12) = { 0.3625, 0,    0, 1.0};
Point(13) = {-0.1,    1.05, 0, 1.0};
Point(14) = { 0.1,    1.05, 0, 1.0};
Point(15) = {-0.2,   -0.25, 0, 1.0};
Point(16) = { 0.2,   -0.25, 0, 1.0};

// Grout curtain under the trench, 0.12 wide, to the base of the section.
Point(21) = {-0.06,  -0.25, 0, 1.0};
Point(22) = { 0.06,  -0.25, 0, 1.0};
Point(23) = {-0.06,  -0.5,  0, 1.0};
Point(24) = { 0.06,  -0.5,  0, 1.0};

// Filters: 0.16 wide measured horizontally, parallel to the core's faces.
// The downstream chimney stands on a blanket drain 0.15 thick that runs out
// to the downstream face (x = 1.95 - 1.5 y there).
Point(31) = {-0.5225, 0,    0, 1.0};
Point(32) = {-0.26,   1.05, 0, 1.0};
Point(33) = { 0.485,  0.15, 0, 1.0};
Point(34) = { 0.26,   1.05, 0, 1.0};
Point(35) = { 1.725,  0.15, 0, 1.0};

// Outer boundary, counter-clockwise.
Line(1)  = {1, 23};
Line(2)  = {23, 24};
Line(3)  = {24, 2};
Line(4)  = {2, 3};
Line(5)  = {3, 6};
Line(6)  = {6, 35};    // end of the blanket drain
Line(7)  = {35, 8};    // downstream face
Line(8)  = {8, 7};     // crest
Line(9)  = {7, 5};     // upstream face
Line(10) = {5, 4};
Line(11) = {4, 1};

// Interfaces.
Line(21) = {5, 31};    // shell | foundation
Line(22) = {31, 11};   // upstream filter | foundation
Line(23) = {11, 15};   // trench wall
Line(24) = {15, 21};   // trench floor
Line(25) = {21, 22};   // core | curtain
Line(26) = {22, 16};   // trench floor
Line(27) = {16, 12};   // trench wall
Line(28) = {12, 6};    // blanket drain | foundation
Line(29) = {21, 23};   // curtain | foundation
Line(30) = {22, 24};
Line(31) = {11, 13};   // core | upstream filter
Line(32) = {12, 14};   // core | downstream filter
Line(33) = {13, 14};   // core | shell
Line(34) = {31, 32};   // upstream filter | shell
Line(35) = {32, 13};
Line(36) = {14, 34};   // downstream filter | shell
Line(37) = {34, 33};   // chimney's outer face
Line(38) = {33, 35};   // top of the blanket drain

Curve Loop(1) = {1, -29, -24, -23, -22, -21, 10, 11};
Plane Surface(1) = {1};                          // foundation, upstream
Curve Loop(2) = {3, 4, 5, -28, -27, -26, 30};
Plane Surface(2) = {2};                          // foundation, downstream
Curve Loop(3) = {2, -30, -25, 29};
Plane Surface(3) = {3};                          // grout curtain
Curve Loop(4) = {23, 24, 25, 26, 27, 32, -33, -31};
Plane Surface(4) = {4};                          // core and cutoff
Curve Loop(5) = {22, 31, -35, -34};
Plane Surface(5) = {5};                          // upstream filter
Curve Loop(6) = {28, 6, -38, -37, -36, -32};
Plane Surface(6) = {6};                          // chimney and blanket drain
Curve Loop(7) = {21, 34, 35, 33, 36, 37, 38, 7, 8, 9};
Plane Surface(7) = {7};                          // rockfill shells

Physical Surface(1) = {4};
Physical Surface(2) = {5, 6};
Physical Surface(3) = {7};
Physical Surface(4) = {1, 2};
Physical Surface(5) = {3};
