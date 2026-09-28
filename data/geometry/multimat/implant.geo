// A threaded dental implant in the jaw, drawn in the half plane r >= 0 with
// the implant's axis on x = 0, as axisymmetric implant models are (and as
// geom010 and icf are). A titanium fixture (mat 1) with five V-threads --
// 60-degree included angle, 0.07 deep, 0.16 pitch, flat crests and roots --
// is seated in bone: a cortical shell (mat 2) 0.21 thick at the top over
// cancellous bone (mat 3). Its smooth collar stands 0.15 proud of the bone
// crest, and its apex is a flat end with a chamfer.
//
// What this domain is for. A periodic sawtooth interface. The titanium |
// bone boundary has 24 corners in a row, alternately convex and reflex for
// each side, and the bone fills every thread root, so the two materials
// interlock tooth into tooth: each bone tooth is 3.3 elements wide at its tip
// and each thread 2 at its crest (NP = 100, h = 0.015 here), so the corners
// of one material are the teeth of the other.
//
// The cortical boundary does not stop at a convenient place: it lands half
// way up a thread flank, so the one interior triple junction (titanium,
// cortical, cancellous) has directions {0, 30} mod 90 and is ILL-POSED, at
// the angle the thread standard sets. The bone crest meets the collar at a
// right angle, and every other landing on the axis and the walls is a right
// angle too. Scale is exaggerated about twofold against a real implant
// (pitch ~0.8 mm on a 4 mm fixture) so that the threads can be resolved.

// ---- generated geometry ----
Point(1) = {0.000000, -0.980000, 0, 1.0};   // apex, on the axis
Point(2) = {0.120000, -0.980000, 0, 1.0};
Point(3) = {0.160000, -0.919585, 0, 1.0};   // first thread root
Point(4) = {0.160000, -0.870415, 0, 1.0};
Point(5) = {0.230000, -0.830000, 0, 1.0};
Point(6) = {0.230000, -0.800000, 0, 1.0};
Point(7) = {0.160000, -0.759585, 0, 1.0};
Point(8) = {0.160000, -0.710415, 0, 1.0};
Point(9) = {0.230000, -0.670000, 0, 1.0};
Point(10) = {0.230000, -0.640000, 0, 1.0};
Point(11) = {0.160000, -0.599585, 0, 1.0};
Point(12) = {0.160000, -0.550415, 0, 1.0};
Point(13) = {0.230000, -0.510000, 0, 1.0};
Point(14) = {0.230000, -0.480000, 0, 1.0};
Point(15) = {0.160000, -0.439585, 0, 1.0};
Point(16) = {0.160000, -0.390415, 0, 1.0};
Point(17) = {0.230000, -0.350000, 0, 1.0};
Point(18) = {0.230000, -0.320000, 0, 1.0};
Point(19) = {0.160000, -0.279585, 0, 1.0};
Point(20) = {0.160000, -0.230415, 0, 1.0};
Point(21) = {0.195359, -0.210000, 0, 1.0};   // cortical boundary meets a thread flank
Point(22) = {0.230000, -0.190000, 0, 1.0};
Point(23) = {0.230000, -0.160000, 0, 1.0};
Point(24) = {0.160000, -0.119585, 0, 1.0};
Point(25) = {0.160000, -0.070415, 0, 1.0};
Point(26) = {0.230000, -0.030000, 0, 1.0};
Point(27) = {0.230000, 0.000000, 0, 1.0};   // collar meets the bone crest
Point(28) = {0.000000, 0.150000, 0, 1.0};   // top of the implant
Point(29) = {0.230000, 0.150000, 0, 1.0};
Point(30) = {0.000000, -1.350000, 0, 1.0};
Point(31) = {0.800000, -1.350000, 0, 1.0};
Point(32) = {0.800000, -0.210000, 0, 1.0};
Point(33) = {0.800000, 0.000000, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 17};
Line(17) = {17, 18};
Line(18) = {18, 19};
Line(19) = {19, 20};
Line(20) = {20, 21};
Line(21) = {21, 22};
Line(22) = {22, 23};
Line(23) = {23, 24};
Line(24) = {24, 25};
Line(25) = {25, 26};
Line(26) = {26, 27};
Line(27) = {27, 29};
Line(28) = {29, 28};
Line(29) = {28, 1};
Line(30) = {21, 32};
Line(31) = {32, 33};
Line(32) = {33, 27};
Line(33) = {30, 31};
Line(34) = {31, 32};
Line(35) = {1, 30};

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
Plane Surface(1) = {1};   // titanium
Curve Loop(2) = {30, 31, 32, -26, -25, -24, -23, -22, -21};
Plane Surface(2) = {2};   // cortical bone
Curve Loop(3) = {33, 34, -30, -20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 35};
Plane Surface(3) = {3};   // cancellous bone

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
