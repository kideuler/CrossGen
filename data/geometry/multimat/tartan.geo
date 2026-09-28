// A tartan: a 5 x 4 Cartesian grid of rectangles whose columns and rows have
// unequal widths, including four thin stripes, coloured with four materials
// so that the four cells round every crossing all differ. Colour of cell
// (i, j) is 1 + (i mod 2) + 2 (j mod 2), so a material never touches itself,
// not even diagonally at a point.
//
// What this domain is for. It is the CONTROL at scale, the way triple_point
// is the control for one junction: twelve interior quadruple points, every one
// an orthogonal crossing (directions 0, 90, 180, 270 -- congruent mod 90, so
// WELL-POSED), and fourteen places an interface lands on the wall, all at
// right angles. The field is constant, no singularity is needed anywhere, and
// the minimal block decomposition is the grid itself: exactly 20 blocks, one
// per cell. Nothing about the answer is in doubt, so a method's count minus 20
// is pure overhead and its worst element is a measure of its own damage.
//
// The unequal widths are the point of not using a checkerboard. The thin
// stripes (0.20 and 0.22 against a 2.69-wide box) are about seven elements
// across at NP = 100, and the wide cells beside them are four times that, so
// a method that splits blocks to even out aspect ratios, or that samples every
// interface at one spacing, gives itself away here without any junction being
// involved. Four materials on twenty surfaces also exercise the
// one-material-many-surfaces path.

xs[] = {0, 0.80, 1.02, 1.57, 1.79, 2.69};   // column edges: 0.80 0.22 0.55 0.22 0.90
ys[] = {0, 0.60, 0.80, 1.55, 1.80};         // row edges:    0.60 0.20 0.75 0.25
nx = #xs[] - 1;
ny = #ys[] - 1;

// Grid point (i, j) is Point(1 + i + (nx + 1) j).
For j In {0:ny}
  For i In {0:nx}
    Point(1 + i + (nx + 1) * j) = {xs[i], ys[j], 0, 1.0};
  EndFor
EndFor

// Horizontal segment from (i, j) to (i + 1, j) is Line(1000 + i + nx j); the
// vertical one from (i, j) to (i, j + 1) is Line(2000 + i + (nx + 1) j). Each
// is shared by the two cells either side of it, so the grid is conformal.
For j In {0:ny}
  For i In {0:nx - 1}
    Line(1000 + i + nx * j) = {1 + i + (nx + 1) * j, 2 + i + (nx + 1) * j};
  EndFor
EndFor
For j In {0:ny - 1}
  For i In {0:nx}
    Line(2000 + i + (nx + 1) * j) = {1 + i + (nx + 1) * j, 1 + i + (nx + 1) * (j + 1)};
  EndFor
EndFor

m1[] = {};
m2[] = {};
m3[] = {};
m4[] = {};
For j In {0:ny - 1}
  For i In {0:nx - 1}
    s = 3000 + i + nx * j;
    Curve Loop(s) = {1000 + i + nx * j, 2000 + (i + 1) + (nx + 1) * j,
                     -(1000 + i + nx * (j + 1)), -(2000 + i + (nx + 1) * j)};
    Plane Surface(s) = {s};
    c = Fmod(i, 2) + 2 * Fmod(j, 2);
    If (c == 0)
      m1[] += s;
    ElseIf (c == 1)
      m2[] += s;
    ElseIf (c == 2)
      m3[] += s;
    Else
      m4[] += s;
    EndIf
  EndFor
EndFor

Physical Surface(1) = {m1[]};
Physical Surface(2) = {m2[]};
Physical Surface(3) = {m3[]};
Physical Surface(4) = {m4[]};
