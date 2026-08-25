// Gmsh project created on Tue Dec 16 15:01:00 2025
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {2, 0, 0, 1.0};
//+
Point(3) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {1, 3};
//+
Ellipse(3) = {2, 1, 2, 3};
//+
Circle(4) = {1.1, 0.3, 0, 0.25, 0, 2*Pi};
//+
Circle(5) = {0.4, 0.6, 0, 0.25, 0, 2*Pi};
//+
Curve Loop(1) = {3, -2, 1};
//+
Curve Loop(2) = {5};
//+
Curve Loop(3) = {4};
//+
Plane Surface(1) = {1, 2, 3};
