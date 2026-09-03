// One off-centre disk: the same problem with none of the symmetry, so that
// a layout that only works when the four cones land on the axes shows up.
//
// Generated for the box-with-disks robustness suite: a 1 x 1 matrix
// block (material 1) with 1 circular inclusion(s) (material 2).
// Written in the same idiom as multimat/bubbles.geo -- four quarter arcs
// per circle, the hole loops reused in the matrix surface so the
// interfaces are conformal on both sides, and the mesh size left to
// Mesh2Dgmsh's uniform field (2/NP after normalisation).

Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};

cx[] = {0.3800};
cy[] = {0.5700};
cr[] = {0.1800};

holes[] = {};
For k In {0:0}
  p = 100 + 5*k;
  c = 200 + 4*k;
  Point(p)     = {cx[k],         cy[k],         0, 1.0};
  Point(p + 1) = {cx[k] + cr[k], cy[k],         0, 1.0};
  Point(p + 2) = {cx[k],         cy[k] + cr[k], 0, 1.0};
  Point(p + 3) = {cx[k] - cr[k], cy[k],         0, 1.0};
  Point(p + 4) = {cx[k],         cy[k] - cr[k], 0, 1.0};
  Circle(c)     = {p + 1, p, p + 2};
  Circle(c + 1) = {p + 2, p, p + 3};
  Circle(c + 2) = {p + 3, p, p + 4};
  Circle(c + 3) = {p + 4, p, p + 1};
  Curve Loop(300 + k) = {c, c + 1, c + 2, c + 3};
  Plane Surface(300 + k) = {300 + k};
  holes[] += 300 + k;
EndFor

Plane Surface(1) = {1, holes[]};
Physical Surface(1) = {1};
Physical Surface(2) = {300:300};
