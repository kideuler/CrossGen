// Gmsh project created on Tue Dec 16 15:14:45 2025
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 01, 0, 1.0};
//+
Point(4) = {0, 01, 0, 1.0};
//+
Circle(1) = {0.2, 0.2, 0, 0.1, 0, 2*Pi};
//+
Circle(2) = {0.4, 0.4, 0, 0.1, 0, 2*Pi};
//+
Circle(3) = {0.2, 0.6, 0, 0.1, 0, 2*Pi};
//+
Circle(4) = {0.3, 0.8, 0, 0.1, 0, 2*Pi};
//+
Circle(5) = {0.8, 0.8, 0, 0.1, 0, 2*Pi};
//+
Circle(6) = {0.6, 0.6, 0, 0.1, 0, 2*Pi};
//+
Circle(7) = {0.7, 0.3, 0, 0.1, 0, 2*Pi};
//+
Circle(8) = {0.5, 0.2, 0, 0.1, 0, 2*Pi};
//+
Line(9) = {1, 2};
//+
Line(10) = {2, 3};
//+
Line(11) = {3, 4};
//+
Line(12) = {4, 1};
//+
Curve Loop(1) = {11, 12, 9, 10};
//+
Curve Loop(2) = {4};
//+
Curve Loop(3) = {3};
//+
Curve Loop(4) = {5};
//+
Curve Loop(5) = {6};
//+
Curve Loop(6) = {2};
//+
Curve Loop(7) = {7};
//+
Curve Loop(8) = {8};
//+
Curve Loop(9) = {1};
//+
Plane Surface(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9};
