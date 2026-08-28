// Concentric shells on a half disk: a core plus three shells, cut on the
// symmetry axis.
//
// This is the target application drawn literally. A two-dimensional
// axisymmetric (r-z) multi-material hydrodynamics problem -- a layered capsule,
// a lined shell -- is exactly a set of nested shells meshed in the half plane
// r >= 0, and the flat side here is that axis.
//
// It stresses interfaces the junction domains do not: every interface is
// curved, they are nested rather than crossing, and their radii are chosen so
// the shells thin outward (0.45 core, then 0.30, 0.15 and 0.10 deep), so the
// outermost is about three elements at NP=70. A field that is not tangent to
// each arc puts mesh lines across the shells, which is the R1 failure that
// forces a mixed-zone closure model right where the density jump is.
//
// The eight vertices where a shell meets the axis are boundary junctions at
// right angles -- each arc crosses the axis perpendicular -- so they are
// benign, and the difficulty is entirely in following four nested curved
// interfaces at once.
Point(1) = {0, 0, 0, 1.0};  // centre, on the axis

Point(2) = {0, -0.45, 0, 1.0};
Point(3) = {0.45, 0, 0, 1.0};
Point(4) = {0, 0.45, 0, 1.0};
Point(5) = {0, -0.75, 0, 1.0};
Point(6) = {0.75, 0, 0, 1.0};
Point(7) = {0, 0.75, 0, 1.0};
Point(8) = {0, -0.9, 0, 1.0};
Point(9) = {0.9, 0, 0, 1.0};
Point(10) = {0, 0.9, 0, 1.0};
Point(11) = {0, -1, 0, 1.0};
Point(12) = {1, 0, 0, 1.0};
Point(13) = {0, 1, 0, 1.0};

// Interface arcs, bottom-to-right and right-to-top for each radius.
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {5, 1, 6};
Circle(4) = {6, 1, 7};
Circle(5) = {8, 1, 9};
Circle(6) = {9, 1, 10};
Circle(7) = {11, 1, 12};
Circle(8) = {12, 1, 13};

// Axis segments, one pair per shell; the core's pair meets at the centre.
Line(9) = {4, 1};
Line(10) = {1, 2};
Line(11) = {7, 4};
Line(12) = {2, 5};
Line(13) = {10, 7};
Line(14) = {5, 8};
Line(15) = {13, 10};
Line(16) = {8, 11};

Curve Loop(1) = {1, 2, 9, 10};
Plane Surface(1) = {1};
Curve Loop(2) = {3, 4, 11, -2, -1, 12};
Plane Surface(2) = {2};
Curve Loop(3) = {5, 6, 13, -4, -3, 14};
Plane Surface(3) = {3};
Curve Loop(4) = {7, 8, 15, -6, -5, 16};
Plane Surface(4) = {4};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
