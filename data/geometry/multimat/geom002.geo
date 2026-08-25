// Half disk (mat 1, from singlemat/geom002) filling the semicircular notch
// of an outer block (mat 2) that now wraps all the way around it -- legs on
// both sides down to the ground line y=0, joined across the top -- forming
// an arch. The notch's two quarter-arcs are the disk's own Circle(1) and
// Circle(2), reused directly rather than redeclared, so the interface is
// conformal on both sides.
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {0, 1, 0, 1.0};
Point(4) = {-1, 0, 0, 1.0};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Line(3) = {4, 1};
Line(4) = {1, 2};
Curve Loop(1) = {2, 3, 4, 1};
Plane Surface(1) = {1};

Point(5) = {2, 0, 0, 1.0};
Point(6) = {2, 2, 0, 1.0};
Point(7) = {-2, 2, 0, 1.0};
Point(8) = {-2, 0, 0, 1.0};
Line(5) = {2, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {8, 4};
Curve Loop(2) = {5, 6, 7, 8, 9, -2, -1};
Plane Surface(2) = {2};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
