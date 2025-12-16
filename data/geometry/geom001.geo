// Gmsh project created on Tue Dec 16 14:30:38 2025
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {1, 3};
//+
Circle(3) = {2, 1, 3};
//+
Curve Loop(1) = {2, -3, -1};
//+
Plane Surface(1) = {1};
