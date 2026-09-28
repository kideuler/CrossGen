// The triple-point problem: the standard multi-material ALE benchmark (used by
// Galera, Maire & Breil, JCP 2010, and in the ReALE paper of Loubere et al.,
// JCP 2010, among many others), drawn at its usual size. A 7 x 3 box holds
// three materials: a high-pressure driver on the left (mat 1, x < 1) and two
// low-pressure materials of different density stacked on the right (mat 2
// below y = 1.5, mat 3 above). The shock from the driver runs faster through
// one than the other, and the vortex that rolls up at the triple point is why
// the problem exists. It is the geometry every multi-material hydro code has
// been run on, which is reason enough to have it here.
//
// What this domain is for. It is the CONTROL with a known answer. The single
// interior junction, at (1, 1.5), is an orthogonal T -- interface directions
// 0, 90 and 270, all congruent mod 90, so WELL-POSED -- and the three places
// an interface lands on the wall are right angles. The minimal conforming
// block decomposition is therefore four rectangles: the left strip has to be
// cut at y = 1.5 so that the T is a vertex of both of its blocks, and nothing
// else has to be cut. Every block beyond four is overhead a method added, and
// that can be read off without a reference solution. It is geom003's question
// (is an orthogonal T handled for free?) asked on a rectangle, where there is
// no curved boundary to muddy the answer.
//
// At NP = 100 (h = 2/NP after Mesh2Dgmsh's normalisation to a box of max
// dimension 2, so 0.07 here) the left strip is 14 elements wide.

Point(1) = {0, 0,   0, 1.0};
Point(2) = {1, 0,   0, 1.0};   // interface 1|2 lands on the bottom wall
Point(3) = {7, 0,   0, 1.0};
Point(4) = {7, 1.5, 0, 1.0};   // interface 2|3 lands on the right wall
Point(5) = {7, 3,   0, 1.0};
Point(6) = {1, 3,   0, 1.0};   // interface 1|3 lands on the top wall
Point(7) = {0, 3,   0, 1.0};
Point(8) = {1, 1.5, 0, 1.0};   // the triple point

// Outer boundary, counter-clockwise, split wherever an interface lands.
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 1};

// Material interfaces, each written once and used by both of its surfaces.
Line(8)  = {2, 8};   // driver | bottom
Line(9)  = {8, 6};   // driver | top
Line(10) = {8, 4};   // bottom | top

Curve Loop(1) = {1, 8, 9, 6, 7};      // driver
Plane Surface(1) = {1};
Curve Loop(2) = {2, 3, -10, -8};      // dense low-pressure material
Plane Surface(2) = {2};
Curve Loop(3) = {10, 4, 5, -9};       // light low-pressure material
Plane Surface(3) = {3};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
