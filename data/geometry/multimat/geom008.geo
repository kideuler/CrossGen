// Tapered thin layer: three stacked materials where the middle one is a few
// elements thick and gets thinner from left to right.
//
// This is the R4/CFL domain rather than a junction domain -- both interfaces
// run wall to wall, so there is no interior junction at all. What it tests is
// what happens when a material region is too thin to hold a transition. At
// NP=70 the layer is about seven elements deep on the left and under three on
// the right, so an irregular vertex landing inside it has nowhere to go, and
// the smallest zone dimension there sets the hydro time step for the whole
// mesh. Reported per material region, not globally: an extraordinary point in
// the thin layer costs far more than the same point in either bulk.
//
// The two interfaces are deliberately not parallel and not quite axis-aligned
// (about 2.9 and -0.4 degrees), so the field has to tilt slightly to stay
// tangent to both while crossing a gap with no room to turn in.
Point(1) = {0, 0, 0, 1.0};
Point(2) = {3, 0, 0, 1.0};
Point(3) = {3, 0.60, 0, 1.0};   // lower interface, right end
Point(4) = {3, 0.73, 0, 1.0};   // upper interface, right end (layer 0.13 deep)
Point(5) = {3, 1.2, 0, 1.0};
Point(6) = {0, 1.2, 0, 1.0};
Point(7) = {0, 0.75, 0, 1.0};   // upper interface, left end
Point(8) = {0, 0.45, 0, 1.0};   // lower interface, left end (layer 0.30 deep)

Line(1) = {1, 2};   // bottom wall
Line(2) = {2, 3};   // right wall, bulk 1
Line(3) = {3, 4};   // right wall, layer
Line(4) = {4, 5};   // right wall, bulk 3
Line(5) = {5, 6};   // top wall
Line(6) = {6, 7};   // left wall, bulk 3
Line(7) = {7, 8};   // left wall, layer
Line(8) = {8, 1};   // left wall, bulk 1

Line(9)  = {8, 3};  // lower interface, shared by materials 1 and 2
Line(10) = {7, 4};  // upper interface, shared by materials 2 and 3

Curve Loop(1) = {1, 2, -9, 8};
Plane Surface(1) = {1};
Curve Loop(2) = {9, 3, -10, 7};
Plane Surface(2) = {2};
Curve Loop(3) = {10, 4, 5, 6};
Plane Surface(3) = {3};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
