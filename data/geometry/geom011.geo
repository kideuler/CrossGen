// Gmsh project created on Wed Jan  7 15:41:40 2026
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 1.5, 1, 0};
//+
Point(5) = {0.75, 1.75, 0, 1.0};
//+
Delete {
  Surface{1}; 
}
//+
Delete {
  Curve{3}; 
}
//+
Point(6) = {0.75, 1, 0, 1.0};
//+
Circle(5) = {3, 6, 5};
//+
Circle(6) = {5, 6, 4};
//+
Circle(7) = {0.3, 1, 0, 0.2, 0, 2*Pi};
//+
Circle(8) = {1.2, 1, 0, 0.2, 0, 2*Pi};
//+
Circle(9) = {0.5, 0.4, 0, 0.2, 0, 2*Pi};
//+
Circle(10) = {1, 0.4, 0, 0.2, 0, 2*Pi};
//+
Ellipse(11) = {0.75, 1, 0, 0.5, 0.25, 0, 2*Pi};
//+
BooleanUnion{ Curve{11}; Delete; }{ Curve{7}; Curve{8}; Delete; }
//+
Delete {
  Curve{13}; Curve{11}; Curve{15}; 
}
//+
Delete {
  Curve{17}; Curve{18}; 
}
//+
Curve Loop(2) = {4, 1, 2, 5, 6};
//+
Curve Loop(3) = {12, 16, 14, 19};
//+
Curve Loop(4) = {9};
//+
Curve Loop(5) = {10};
//+
Plane Surface(1) = {2, 3, 4, 5};
