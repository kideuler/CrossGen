// Gmsh project created on Tue Dec 16 15:58:17 2025
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {1, 0, 0, 1.0};
//+
Circle(1) = {0.5, 0.5, 0, 0.25, 0, 2*Pi};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Point(6) = {0, 1, 0, 1.0};
//+
Line(4) = {3, 6};
//+
Line(5) = {6, 1};
//+
Curve Loop(1) = {5, 2, 3, 4};
//+
Curve Loop(2) = {1};
//+
Plane Surface(1) = {1, 2};
