// A sedimentary basin: three flat-lying fill units (mat 2 lowest, mat 3,
// mat 4 at the surface) resting in a trough cut into basement rock (mat 1),
// with a lenticular channel sand (mat 5) buried in the middle unit. The
// trough is asymmetric -- a 35-degree flank on the left, 50 on the right,
// joined by a flat floor -- as a half-graben or a glacial valley would be.
// Basins like this are what site-response and groundwater models are built
// on, and their characteristic geometry is the ONLAP PINCH-OUT: every fill
// unit thins to nothing where it laps onto the basement.
//
// What this domain is for. Acute material wedges. Each bed boundary meets a
// flank at the flank's angle, so each fill unit ends in a sliver of 35
// degrees on the left and 50 on the right: six wedges in the three fill
// units, two of them where the surface unit laps onto the ground surface at
// the margins and four inside the domain. Those four are also triple
// junctions (basement, unit k, unit k+1) with directions {0, 145} or
// {0, 50}: mod 90 that is {0, 55} and {0, 50}, so all four are ILL-POSED.
// A 35-degree corner rounds to one quarter turn with a 55-degree error, and a
// 50-degree corner sits five degrees past the half-way point of the rounding,
// so the two flanks ask the rounding question from opposite sides.
//
// The channel sand is a biconvex lens -- two circular arcs meeting at
// 40-degree tips -- inside unit 3, touching nothing. It is an inclusion with
// two acute corners and no junction, the multi-material cousin of the
// singlemat crescent (geom024) that no whole-region template fits. It is not
// a full circle, so its interior singularities are pinned by its tips.
//
// Resolution at NP = 100 (h = 0.05 here): the fill units are 6, 9 and 7
// elements deep; the lens is 3.2 elements thick at its middle and leaves 2.9
// elements of unit 3 above and below it. The lens arcs' centres lie outside
// the domain (their radius is 1.31), but inside the box's width, so they do
// not change Mesh2Dgmsh's normalisation.

Point(1) = {0,   0,   0, 1.0};
Point(2) = {5,   0,   0, 1.0};
Point(3) = {5,   1.6, 0, 1.0};
Point(4) = {4.4, 1.6, 0, 1.0};   // right margin: basement meets the surface
Point(5) = {0.6, 1.6, 0, 1.0};   // left margin
Point(6) = {0,   1.6, 0, 1.0};

// Left flank, 35 degrees: x = 0.6 + (1.6 - y) / tan(35).
Point(11) = {1.099852, 1.25, 0, 1.0};   // 3|4 pinches out
Point(12) = {1.742518, 0.80, 0, 1.0};   // 2|3 pinches out
Point(13) = {2.170963, 0.50, 0, 1.0};   // foot of the flank
// Right flank, 50 degrees: x = 4.4 - (1.6 - y) / tan(50).
Point(14) = {3.476990, 0.50, 0, 1.0};
Point(15) = {3.728720, 0.80, 0, 1.0};
Point(16) = {4.106315, 1.25, 0, 1.0};

// Channel lens: tips 0.9 apart on y = 1.025, each arc 0.08 deep.
Point(21) = {2.45, 1.025,    0, 1.0};
Point(22) = {3.35, 1.025,    0, 1.0};
Point(23) = {2.90, -0.200625, 0, 1.0};  // centre of the upper arc
Point(24) = {2.90, 2.250625,  0, 1.0};  // centre of the lower arc

// Outer boundary.
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

// The unconformity (basement | fill), left margin to right margin.
Line(11) = {5, 11};
Line(12) = {11, 12};
Line(13) = {12, 13};
Line(14) = {13, 14};
Line(15) = {14, 15};
Line(16) = {15, 16};
Line(17) = {16, 4};

// Bed boundaries, flank to flank.
Line(21) = {12, 15};   // 2|3
Line(22) = {11, 16};   // 3|4

// The lens.
Circle(31) = {22, 23, 21};   // top
Circle(32) = {21, 24, 22};   // bottom

Curve Loop(1) = {1, 2, 3, -17, -16, -15, -14, -13, -12, -11, 5, 6};
Plane Surface(1) = {1};                 // basement
Curve Loop(2) = {14, 15, -21, 13};
Plane Surface(2) = {2};                 // lower fill
Curve Loop(5) = {31, 32};
Curve Loop(3) = {21, 16, -22, 12};
Plane Surface(3) = {3, 5};              // middle fill, with the lens cut out
Curve Loop(4) = {22, 17, 4, 11};
Plane Surface(4) = {4};                 // upper fill
Plane Surface(5) = {5};                 // channel sand

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
