// Connecting rod: two bearing eyes joined by a waisted shank, with the big-end
// cap bolts drilled through.
//
// Four holes, so the Euler characteristic is -3 and this is the most connected
// model in data/geometry after geom023's perforated plate -- but unlike that
// plate the holes are all different sizes and none of them is where the shape
// is widest. The two bolt holes sit in the big-end ring with 0.055 of material
// inside them and 0.055 outside, ligaments thinner than anything except the
// tongues of geom020, and the field has to get round both sides of each.
//
// The shank is the other half of the model. It is a channel about 1.5 long and
// 0.17 wide at its narrowest, tapering the whole way, so the field runs almost
// parallel to both flanks for its full length: the same configuration as
// geom022's dumbbell neck, but tapered rather than uniform, and the chords the
// layout finds along it are correspondingly not parallel.
//
// The two ends of the shank are what the model is really for, because they are
// the same physical feature and the mesh reads them differently. Measured:
// where the flank leaves the big end the boundary turns through 8 degrees, an
// interior angle of 187.6, which rounds to pi -- not a corner at all, no node,
// nothing launched. Where it meets the small end the turn is 66 degrees, an
// interior angle of 246.0, which rounds to 270: a reflex corner at -1/4
// launching two rays. Both readings are stable, 186.6 to 187.9 and 243.8 to
// 246.9 across np 75 to 250. A real rod is filleted at both ends; drawn this
// way, the classification is visibly a property of the angle and not of what
// the feature is called, and both sides of the 225-degree threshold appear on
// one part.
//
// Above np=75 those two joins are the only corners on the whole boundary, so it
// carries -1/2 and the remaining -5/2 goes somewhere else. It goes onto the
// boundary rather than into the interior: the field puts 12 reflex boundary
// singularities on the four smooth hole loops, summing to -3, and UMBER ends
// with no internal singularities at all.
//
//     np    triangles   corners        layout            note
//     50        1064     0 / 16        33 -> 27 blocks    2 not four-sided
//     75        2308     0 / 2         24 -> 18 blocks    2 not four-sided
//    150        8852     0 / 2         27 -> 21 blocks    5 not four-sided
//
// np=50, which is what data/meshes carries, is the one to be careful with: at
// that size the 0.055 bolt holes are 8 or 9 elements round, so the mesh reads
// each of them as a ring of 230-degree corners rather than as a circle, and the
// corner count goes from 2 to 16. It is the same effect as geom029's leading
// edge and it is the reason this model is worth meshing finer than the default.
// At np=150 the parameterization also produces four flipped triangles, which
// the coarser meshes do not.
SetFactory("OpenCASCADE");

Rb  = 0.50;   // big end, outside
rb  = 0.28;   // big end, bore
Rs  = 0.20;   // small end, outside
rs  = 0.115;  // small end, bore
xs  = 1.9;    // small end centre, on the x axis
ab  = 62*Pi/180;   // where the shank leaves the big end
as_ = 25*Pi/180;   // and where it meets the small end, measured from -x

Point(1) = {0, 0, 0, 1.0};    // big end centre
Point(2) = {xs, 0, 0, 1.0};   // small end centre

// Attachments and the two intermediate points that let each 236- and 310-degree
// arc be built from arcs of less than half a turn.
Point(10) = { Rb*Cos(ab),  Rb*Sin(ab), 0, 1.0};
Point(11) = { Rb*Cos(-ab), Rb*Sin(-ab), 0, 1.0};
Point(12) = {-Rb, 0, 0, 1.0};
Point(13) = { Rb*Cos(2*Pi/3), Rb*Sin(2*Pi/3), 0, 1.0};
Point(14) = { Rb*Cos(-2*Pi/3), Rb*Sin(-2*Pi/3), 0, 1.0};

Point(20) = {xs + Rs*Cos(Pi - as_), Rs*Sin(Pi - as_), 0, 1.0};        // upper join
Point(21) = {xs + Rs*Cos(Pi + as_),  Rs*Sin(Pi + as_), 0, 1.0};       // lower join
Point(22) = {xs + Rs, 0, 0, 1.0};
Point(23) = {xs + Rs*Cos(Pi/3),  Rs*Sin(Pi/3), 0, 1.0};
Point(24) = {xs + Rs*Cos(-Pi/3), Rs*Sin(-Pi/3), 0, 1.0};

// The shank flanks. Four intermediate points each, tapering from 0.44 of half
// width at the big end to 0.085 at the small one.
Point(30) = {0.55, -0.33,  0, 1.0};
Point(31) = {0.90, -0.22,  0, 1.0};
Point(32) = {1.25, -0.145, 0, 1.0};
Point(33) = {1.50, -0.105, 0, 1.0};
Point(40) = {0.55,  0.33,  0, 1.0};
Point(41) = {0.90,  0.22,  0, 1.0};
Point(42) = {1.25,  0.145, 0, 1.0};
Point(43) = {1.50,  0.105, 0, 1.0};

// Outer boundary, counter-clockwise from the top of the big end round the back.
Circle(1) = {10, 1, 13};
Circle(2) = {13, 1, 12};
Circle(3) = {12, 1, 14};
Circle(4) = {14, 1, 11};
Spline(5) = {11, 30, 31, 32, 33, 21};
Circle(6) = {21, 2, 24};
Circle(7) = {24, 2, 22};
Circle(8) = {22, 2, 23};
Circle(9) = {23, 2, 20};
Spline(10) = {20, 43, 42, 41, 40, 10};

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
Plane Surface(1) = {1};

// The bores and the two cap bolts, at 105 and 255 degrees on a 0.39 bolt
// circle: centred in the big-end ring, which is what makes both ligaments the
// same thickness.
Disk(2) = {0, 0, 0, rb};
Disk(3) = {xs, 0, 0, rs};
Disk(4) = {0.39*Cos(105*Pi/180), 0.39*Sin(105*Pi/180), 0, 0.055};
Disk(5) = {0.39*Cos(255*Pi/180), 0.39*Sin(255*Pi/180), 0, 0.055};

BooleanDifference(6) = { Surface{1}; Delete; }
                      { Surface{2}; Surface{3}; Surface{4}; Surface{5}; Delete; };
