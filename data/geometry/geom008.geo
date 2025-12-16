// Gmsh project created on Tue Dec 16 15:03:47 2025
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {2, 0, 0, 1.0};
//+
Point(3) = {-2, 0, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Ellipse(1) = {2, 1, 2, 4};
//+
Ellipse(2) = {4, 1, 3, 3};
//+
Line(3) = {3, 1};
//+
Line(4) = {1, 2};
//+
Circle(5) = {-1.3, 0.3, 0, 0.2, 0, 2*Pi};
//+
Circle(6) = {0, 0.6, 0, 0.2, 0, 2*Pi};
//+
Circle(7) = {1.2, 0.3, 0, 0.2, 0, 2*Pi};
//+
Curve Loop(1) = {2, 3, 4, 1};
//+
Curve Loop(2) = {5};
//+
Curve Loop(3) = {6};
//+
Curve Loop(4) = {7};
//+
Plane Surface(1) = {1, 2, 3, 4};
