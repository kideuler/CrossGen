// Pillow block: a bearing housing on a flanged base, with slotted mounting
// feet and a boss for the grease fitting.
//
// The shape is a union of four primitives -- base, boss, web, fitting pad --
// and the interesting part is what a union does to a boundary. Three of the
// joins are not designed corners at all; they are the incidental curves where
// one primitive breaks the surface of another, and they are where the field is
// hardest to guess:
//
//   * the web sides cut the top face of the base at 72 degrees, giving two
//     corners the mesh measures at 251.7 -- reflex, -1/4, two rays each;
//   * the same sides cross the boss circle further up, and those two measure
//     219.5, five degrees under the 225 threshold, so the mesh rounds them to
//     pi: no node, no index, nothing launched, on a kink that is perfectly
//     visible in the outline;
//   * the fitting pad crosses the top of the boss nearly square on, and those
//     two measure 261.0, so they count.
//
// Three pairs of joins of the same kind on one casting, two pairs classified
// and one not, and the difference is only the angle at which a straight line
// happens to cut a circle. That pair is the closest thing here to a coin toss,
// and it does not flip: it measures 218.6, 219.5 and 219.8 at np 75, 150 and
// 250, walking towards the threshold and not reaching it.
//
// Below that, the base is a strip 0.32 thick with two mounting slots cut in it,
// each 0.35 long and 0.15 wide, leaving 0.085 of material above and below. The
// slots are holes and so is the bore, so the Euler characteristic is -2. The
// corners give 6 convex and 4 reflex, +1/2, and the field puts the remaining
// -5/2 on the boundary as 20 boundary singularities, 6 convex and 14 reflex,
// with no internal singularities left at the end.
//
//     np    triangles   motorcycles   raw layout   after collapse
//     50        1841        28          241          223, 1 ray short
//     75        3893        28          120          106, 2 not four-sided
//    150       15076        28           61           32, all four-sided
//
// This one wants a fine mesh, which is the opposite of what geom029 reports.
// The field and the corner inventory are the same at all three sizes; what
// changes is the tracing, which crosses itself 216 times at np=50 and 35 times
// at np=150. The 0.075 slot ends, with 0.085 of material above and below them,
// are what the coarse meshes cannot carry. data/meshes carries np=50, so the
// shipped mesh is the one that misbehaves; np=150 is the one to quote.
SetFactory("OpenCASCADE");

// Base flange, boss, web, and the pad the grease nipple screws into.
Rectangle(1) = {0, 0, 0, 2.4, 0.32};
Disk(2) = {1.2, 1.0, 0, 0.56};

Point(100) = {0.42, 0.0, 0, 1.0};
Point(101) = {1.98, 0.0, 0, 1.0};
Point(102) = {1.65, 1.0, 0, 1.0};
Point(103) = {0.75, 1.0, 0, 1.0};
Line(100) = {100, 101};
Line(101) = {101, 102};
Line(102) = {102, 103};
Line(103) = {103, 100};
Curve Loop(100) = {100, 101, 102, 103};
Plane Surface(3) = {100};

Rectangle(4) = {1.12, 1.40, 0, 0.16, 0.28};

BooleanUnion(10) = { Surface{1}; Delete; }
                  { Surface{2}; Surface{3}; Surface{4}; Delete; };

// The bore, and the two mounting slots: milled 0.15 wide with round ends, so
// the housing can be shuffled on its bolts before it is tightened down.
Disk(20) = {1.2, 1.0, 0, 0.34};

Disk(21) = {0.14, 0.16, 0, 0.075};
Disk(22) = {0.34, 0.16, 0, 0.075};
Rectangle(23) = {0.14, 0.085, 0, 0.20, 0.15};

Disk(24) = {2.06, 0.16, 0, 0.075};
Disk(25) = {2.26, 0.16, 0, 0.075};
Rectangle(26) = {2.06, 0.085, 0, 0.20, 0.15};

BooleanUnion(30) = { Surface{21}; Delete; }{ Surface{22}; Surface{23}; Delete; };
BooleanUnion(31) = { Surface{24}; Delete; }{ Surface{25}; Surface{26}; Delete; };

BooleanDifference(40) = { Surface{10}; Delete; }
                       { Surface{20}; Surface{30}; Surface{31}; Delete; };
