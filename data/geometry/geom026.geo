// Free-form blob with two holes of different size.
//
// A closed spline boundary: curvature varies continuously and there is not a
// single corner to pin a singularity to, so every singularity in the field is
// placed by the boundary curvature alone. The two holes then add index that has
// to be reconciled with it.
SetFactory("OpenCASCADE");
Point(1) = { 1.60,  0.00, 0, 1.0};
Point(2) = { 1.10,  0.95, 0, 1.0};
Point(3) = { 0.10,  1.25, 0, 1.0};
Point(4) = {-0.85,  0.90, 0, 1.0};
Point(5) = {-1.45,  0.10, 0, 1.0};
Point(6) = {-1.05, -0.85, 0, 1.0};
Point(7) = { 0.05, -1.20, 0, 1.0};
Point(8) = { 1.05, -0.80, 0, 1.0};
BSpline(1) = {1, 2, 3, 4, 5, 6, 7, 8, 1};
Curve Loop(1) = {1};
Plane Surface(1) = {1};
Disk(2) = {-0.55, 0.10, 0, 0.30, 0.30};
Disk(3) = { 0.60, -0.20, 0, 0.18, 0.18};
BooleanDifference(4) = { Surface{1}; Delete; }{ Surface{2}; Surface{3}; Delete; };
