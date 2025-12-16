// Gmsh project created on Tue Dec 16 15:53:25 2025
SetFactory("OpenCASCADE");
//+
Point(1) = {1, 0, 0, 1.0};
//+
Point(2) = {2, 0, 0, 1.0};
//+
Point(3) = {2, 2, 0, 1.0};
//+
Point(4) = {0, 2, 0, 1.0};
//+
Point(5) = {0, 1, 0, 1.0};
//+
Point(6) = {0, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Circle(5) = {1, 6, 5};
//+
Curve Loop(1) = {4, -5, 1, 2, 3};
//+
Plane Surface(1) = {1};
