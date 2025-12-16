// Gmsh project created on Tue Dec 16 15:55:45 2025
SetFactory("OpenCASCADE");
//+
Point(1) = {-2, 0, 0, 1.0};
//+
Point(2) = {-1, 0, 0, 1.0};
//+
Point(3) = {0, 0, 0, 1.0};
//+
Point(4) = {1, 0, 0, 1.0};
//+
Point(5) = {2, 0, 0, 1.0};
//+
Point(6) = {0, 1, 0, 1.0};
//+
Point(7) = {2, 2, 0, 1.0};
//+
Point(8) = {-2, 2, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {4, 5};
//+
Line(3) = {5, 7};
//+
Line(4) = {7, 8};
//+
Line(5) = {8, 1};
//+
Circle(6) = {4, 3, 6};
//+
Circle(7) = {6, 3, 2};
//+
Curve Loop(1) = {7, -1, -5, -4, -3, -2, 6};
//+
Plane Surface(1) = {1};
