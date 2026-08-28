// Corrugated interface between two materials: y = 0.15*cos(2*pi*x) across a
// rectangle, two full wavelengths.
//
// The perturbed planar interface is the standard multi-material hydrodynamics
// test configuration (a Richtmyer-Meshkov initial condition), so this is the
// domain in the set an application reader will recognise on sight.
//
// It is built so the interface is the only thing under test. The cosine has
// zero slope where it meets the side walls, so the interface tangent is
// horizontal and the wall tangent vertical there -- congruent mod 90, so both
// boundary junctions are benign and cannot be blamed for anything. What is
// left is a long, strongly curved interface separating two materials whose
// bulk fields both want to be axis-aligned: if diffusion leaks across it, the
// corrugation prints through into the other material, and that print-through
// is directly measurable.
nwave = 120;   // interface sample count
amp = 0.15;

For i In {0:nwave}
  xw = -1 + 2*i/nwave;
  Point(100 + i) = {xw, amp*Cos(2*Pi*xw), 0, 1.0};
EndFor

Point(1) = {-1, -0.6, 0, 1.0};
Point(2) = { 1, -0.6, 0, 1.0};
Point(3) = { 1,  0.6, 0, 1.0};
Point(4) = {-1,  0.6, 0, 1.0};

// Spline endpoints are Point(100) on the left wall and Point(100+nwave) on the
// right wall, so the walls are split by the interface and both surfaces share
// exactly these curves.
Spline(1) = {100:100+nwave};

Line(2) = {100, 1};            // left wall, below the interface
Line(3) = {1, 2};              // bottom wall
Line(4) = {2, 100 + nwave};    // right wall, below the interface
Line(5) = {100 + nwave, 3};    // right wall, above the interface
Line(6) = {3, 4};              // top wall
Line(7) = {4, 100};            // left wall, above the interface

Curve Loop(1) = {2, 3, 4, -1};
Plane Surface(1) = {1};
Curve Loop(2) = {1, 5, 6, 7};
Plane Surface(2) = {2};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
