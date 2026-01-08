// Gmsh project created on Wed Jan  7 19:16:28 2026
SetFactory("OpenCASCADE");
//+
Ellipse(1) = {-0, 0, 0, 1, 0.75, 0, 2*Pi};
//+
Ellipse(2) = {-0, -1, 0, 1, 0.75, 0, 2*Pi};
//+
BooleanUnion{ Curve{2}; Delete; }{ Curve{1}; Delete; }
//+
Delete {
  Curve{2}; Curve{5}; 
}
//+
Circle(7) = {-0.624, 0.253, 0, 0.2, 0, 2*Pi};
//+
Circle(8) = {-0.076, 0.402, 0, 0.2, 0, 2*Pi};
//+
Circle(9) = {-0.248, -0.178, 0, 0.2, 0, 2*Pi};
//+
Circle(10) = {0.586, 0.119, 0, 0.2, 0, 2*Pi};
//+
Circle(11) = {0.248, -0.396, 0, 0.2, 0, 2*Pi};
//+
Circle(12) = {-0.633, -1.008, 0, 0.2, 0, 2*Pi};
//+
Circle(13) = {-0.253, -0.655, 0, 0.2, 0, 2*Pi};
//+
Rectangle(1) = {-0.359, -1.624, 0, 0.5, 0.5, 0};
//+
Rectangle(2) = {0.262, -1.225, 0, 0.5, 0.5, 0};
//+
Delete {
  Surface{1}; 
}
//+
Delete {
  Surface{2}; 
}
//+
Curve Loop(3) = {4, 3, 1, 6};
//+
Curve Loop(4) = {7};
//+
Curve Loop(5) = {8};
//+
Curve Loop(6) = {10};
//+
Curve Loop(7) = {11};
//+
Curve Loop(8) = {9};
//+
Curve Loop(9) = {13};
//+
Curve Loop(10) = {12};
//+
Curve Loop(11) = {16, 17, 14, 15};
//+
Curve Loop(12) = {20, 21, 18, 19};
//+
Plane Surface(1) = {3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
