// Gmsh project created on Tue Dec 30 14:54:57 2025
SetFactory("OpenCASCADE");
//+
Ellipse(1) = {0, -0, 0, 0.5, 0.25, 0, 2*Pi};
//+
Circle(2) = {-0, 0.2, 0, 0.3, 0, 2*Pi};
//+
Circle(3) = {0, -0.5, 0, 0.3, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {3};
//+
Curve Loop(4) = {3};
//+
Plane Surface(3) = {4};
//+
BooleanUnion{ Surface{1}; Delete; }{ Surface{2}; Delete; }
//+
BooleanUnion{ Surface{1}; Delete; }{ Surface{3}; Delete; }
