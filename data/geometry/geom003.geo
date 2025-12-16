// Gmsh project created on Tue Dec 16 14:40:45 2025
//+
SetFactory("OpenCASCADE");
Disk(1) = {0, 0, 0, 0.5, 0.5};
//+
Curve Loop(3) = {1};
//+
Plane Surface(3) = {3};
