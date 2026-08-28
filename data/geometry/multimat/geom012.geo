// Polycrystalline grain structure: five grains meeting at four interior triple
// junctions, all at generic angles.
//
// The other multi-material domains isolate one junction so it can be studied;
// this one supplies volume. A junction-vertex table needs more than one row per
// model to say anything, and here a single mesh contributes four interior
// triple junctions and four more where an interface terminates on the wall,
// none of them axis-aligned and none of them congruent mod 90. At J1, for
// instance, the three incident interface directions reduce to about 73, 12 and
// 79 degrees, so no cross is tangent to more than one of them.
//
// Grain structures are also a real target: a polycrystal is meshed for
// micromechanics much the way a multi-material hydro domain is, and grain
// boundaries have to be mesh lines for the same reason material interfaces do.
Point(1) = {0, 0, 0, 1.0};
Point(2) = {0.55, 0, 0, 1.0};      // interface meets the bottom wall
Point(3) = {1.70, 0, 0, 1.0};      // interface meets the bottom wall
Point(4) = {2, 0, 0, 1.0};
Point(5) = {2, 1.60, 0, 1.0};      // interface meets the right wall
Point(6) = {2, 2, 0, 1.0};
Point(7) = {0.75, 2, 0, 1.0};      // interface meets the top wall
Point(8) = {0, 2, 0, 1.0};

Point(9)  = {0.70, 0.75, 0, 1.0};  // J1
Point(10) = {1.35, 0.55, 0, 1.0};  // J2
Point(11) = {1.25, 1.35, 0, 1.0};  // J3
Point(12) = {0.55, 1.45, 0, 1.0};  // J4

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};

Line(9)  = {2, 9};    // wall to J1
Line(10) = {3, 10};   // wall to J2
Line(11) = {5, 11};   // wall to J3
Line(12) = {7, 12};   // wall to J4
Line(13) = {9, 10};   // J1-J2
Line(14) = {10, 11};  // J2-J3
Line(15) = {11, 12};  // J3-J4
Line(16) = {12, 9};   // J4-J1

Curve Loop(1) = {13, 14, 15, 16};          // central grain
Plane Surface(1) = {1};
Curve Loop(2) = {2, 10, -13, -9};          // bottom grain
Plane Surface(2) = {2};
Curve Loop(3) = {3, 4, 11, -14, -10};      // right grain
Plane Surface(3) = {3};
Curve Loop(4) = {5, 6, 12, -15, -11};      // top-right grain
Plane Surface(4) = {4};
Curve Loop(5) = {7, 8, 1, 9, -16, -12};    // left grain
Plane Surface(5) = {5};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
