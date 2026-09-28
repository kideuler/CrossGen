// A wound ("jelly roll") cell in cross-section, as in a cylindrical battery
// or capacitor. Two electrode strips (mat 1 and mat 2) are wound together
// round a hollow mandrel -- a void, r < 0.2 -- and fill the steel can (mat
// 3, 0.08 thick) out to r = 1. The two strips are separated by two spiral
// interfaces, the second the first turned through half a turn, each winding
// twice round the mandrel. A real roll has dozens of turns; two is what the
// resolution holds.
//
// What this domain is for. Its answer is tiny and nothing local finds it.
// Each strip is a topological disk with exactly four corners -- the two
// places its spirals leave the mandrel and the two where they reach the can
// -- so each is ONE block, a long curved rectangle wound twice round the
// hole, and the two share both spirals whole. With the can split where the
// spirals land, the minimal block decomposition has four blocks. A layout
// built from local decisions (a singularity here, a separatrix there) cannot
// see that a strip closes up two turns later, and the count it produces
// instead measures how badly it fails to. The strips are 0.10 thick at
// mid-radius (4.6 elements at NP = 100, h = 0.0216 here) and much thicker
// near the ends, which is where the pitch is.
//
// Each spiral leaves the mandrel and reaches the can RADIALLY, so the four
// landings are right angles and no strip has an acute corner: the spirals are
// r = 0.2 + 0.8 t, theta = 4 pi (t - sin(2 pi t) / 2 pi), t in [0, 1], whose
// winding rate vanishes at both ends. Two of the four landings are on the
// rim of the mandrel, which is a void rather than a disk, so nothing here is
// a full-circle material region.

r_in  = 0.2;       // mandrel
r_can = 1.0;       // inside of the can
r_out = 1.08;      // outside of the can
Theta = 4 * Pi;    // two turns
N     = 160;       // spline points per spiral

Point(1) = {0, 0, 0, 1.0};   // centre of every circle

// Quarter points of the three circles (a GEO Circle spans less than 180
// degrees). The spirals start at the mandrel's 0 and 180-degree points and
// end at the same angles on the can, two turns later.
Point(2)  = { r_in,  0,     0, 1.0};
Point(3)  = { 0,     r_in,  0, 1.0};
Point(4)  = {-r_in,  0,     0, 1.0};
Point(5)  = { 0,    -r_in,  0, 1.0};
Point(6)  = { r_can, 0,     0, 1.0};
Point(7)  = { 0,     r_can, 0, 1.0};
Point(8)  = {-r_can, 0,     0, 1.0};
Point(9)  = { 0,    -r_can, 0, 1.0};
Point(10) = { r_out, 0,     0, 1.0};
Point(11) = { 0,     r_out, 0, 1.0};
Point(12) = {-r_out, 0,     0, 1.0};
Point(13) = { 0,    -r_out, 0, 1.0};

// Interior points of the two spirals.
For i In {1:N - 1}
  t  = i / N;
  r  = r_in + (r_can - r_in) * t;
  th = Theta * (t - Sin(2 * Pi * t) / (2 * Pi));
  Point(1000 + i) = {r * Cos(th),      r * Sin(th),      0, 1.0};
  Point(3000 + i) = {r * Cos(th + Pi), r * Sin(th + Pi), 0, 1.0};
EndFor
Spline(1) = {2, 1001:1000 + N - 1, 6};   // first spiral, mandrel to can
Spline(2) = {4, 3001:3000 + N - 1, 8};   // second, half a turn behind

Circle(11) = {2, 1, 3};    // mandrel
Circle(12) = {3, 1, 4};
Circle(13) = {4, 1, 5};
Circle(14) = {5, 1, 2};
Circle(21) = {6, 1, 7};    // inside of the can
Circle(22) = {7, 1, 8};
Circle(23) = {8, 1, 9};
Circle(24) = {9, 1, 6};
Circle(31) = {10, 1, 11};  // outside of the can
Circle(32) = {11, 1, 12};
Circle(33) = {12, 1, 13};
Circle(34) = {13, 1, 10};

// Each strip: out along one spiral, half way round the can, back in along
// the other spiral, half way round the mandrel.
Curve Loop(1) = {1, 21, 22, -2, -12, -11};
Plane Surface(1) = {1};              // first electrode
Curve Loop(2) = {2, 23, 24, -1, -14, -13};
Plane Surface(2) = {2};              // second electrode
Curve Loop(3) = {31, 32, 33, 34};
Curve Loop(4) = {21, 22, 23, 24};
Plane Surface(3) = {3, 4};           // can

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
