// Rounded plate with a vee notch: smooth boundary meeting a sharp one.
//
// The rounded corners carry no index at all -- they are a straight boundary as
// far as Table 1 is concerned -- while the notch has a sharp convex tip and two
// reflex corners where it meets the edge. Mixing the two in one model is what
// separates a corner the geometry really has from one the discretisation
// invented, which is what the interior-angle bands are deciding.
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 2.4, 1.6, 0.35};
Point(100) = {0.9, 1.6, 0, 1.0};
Point(101) = {1.2, 0.7, 0, 1.0};
Point(102) = {1.5, 1.6, 0, 1.0};
Point(103) = {1.5, 1.9, 0, 1.0};
Point(104) = {0.9, 1.9, 0, 1.0};
Line(100) = {100, 101};
Line(101) = {101, 102};
Line(102) = {102, 103};
Line(103) = {103, 104};
Line(104) = {104, 100};
Curve Loop(100) = {100, 101, 102, 103, 104};
Plane Surface(100) = {100};
BooleanDifference(200) = { Surface{1}; Delete; }{ Surface{100}; Delete; };
