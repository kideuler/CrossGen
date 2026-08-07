// Six-pointed star: convex and reflex corners alternating all the way round,
// with a round hole in the middle.
//
// Twelve corners, six of them reflex, laid out with exact rotational symmetry,
// and a hole whose index has to be shared out among them. The symmetry means
// every singularity has five others that ought to be its mirror image, so any
// asymmetry in the partition is numerical.
SetFactory("OpenCASCADE");
R_out = 1.0;
R_in  = 0.45;
For k In {0:5}
  a_out = k * Pi/3;
  a_in  = a_out + Pi/6;
  Point(1 + 2*k) = {R_out*Cos(a_out), R_out*Sin(a_out), 0, 1.0};
  Point(2 + 2*k) = {R_in *Cos(a_in),  R_in *Sin(a_in),  0, 1.0};
EndFor
For k In {1:11}
  Line(k) = {k, k+1};
EndFor
Line(12) = {12, 1};
Curve Loop(1) = {1:12};
Plane Surface(1) = {1};
Disk(2) = {0, 0, 0, 0.18, 0.18};
BooleanDifference(3) = { Surface{1}; Delete; }{ Surface{2}; Delete; };
