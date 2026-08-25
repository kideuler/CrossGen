// Gmsh project created on Sat Jan 17 16:42:31 2026
SetFactory("OpenCASCADE");
//+
Ellipse(1) = {0, 0.5, 0, 0.5, 0.25, 0, 2*Pi};
//+
Rectangle(1) = {-0.8, -0.5, 0, 1.5, 1, 0};
//+
Curve Loop(2) = {1};
//+
Plane Surface(2) = {2};
//+
BooleanUnion{ Surface{1}; Delete; }{ Surface{2}; Delete; }
//+
Disk(2) = {-0.01, 0.48, 0, 0.4, 0.25};
//+
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }
//+
Rectangle(2) = {0.38, -0.48, -0, 0.3, 0.3, 0};
//+
Rectangle(3) = {-0.78, 0.18, -0, 0.3, 0.3, 0};
//+
BooleanDifference{ Surface{1}; Delete; }{ Surface{3}; Surface{2}; Delete; }
