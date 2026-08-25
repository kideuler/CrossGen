// Quarter disk (mat 1, from singlemat/geom001) sharing its arc, conformally,
// with the surrounding box (mat 2, from singlemat/geom004) that completes it
// to a full [0,2]x[0,2] square. The shared Circle(3) is the material
// interface: used forward in the inner loop, reversed in the outer loop, so
// Gmsh meshes it once and both surfaces share its nodes.
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1.5, 0, 0, 1.0};
Point(3) = {0, 1.5, 0, 1.0};
Point(4) = {2, 0, 0, 1.0};
Point(5) = {2, 2, 0, 1.0};
Point(6) = {0, 2, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {1, 3};
Circle(3) = {2, 1, 3};
Line(4) = {2, 4};
Line(5) = {4, 5};
Line(6) = {5, 6};
Line(7) = {6, 3};

Curve Loop(1) = {2, -3, -1};
Plane Surface(1) = {1};

Curve Loop(2) = {4, 5, 6, 7, -3};
Plane Surface(2) = {2};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
