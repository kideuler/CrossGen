// Hydraulic manifold block: a rectangle of aluminium with the passages drilled
// and milled into it, sectioned through the galleries.
//
// This is the model for internal features. Everything in data/geometry up to
// here either has its complications on the outside (geom020's slots, geom028's
// star) or has holes that are round (geom023, geom032). Here the outline is a
// plain rectangle -- four corners, +1 of index, and nothing else -- and every
// difficult thing is inside it:
//
//   * a blind main gallery with rounded ends, whose left end stops 0.17 short
//     of the face, so the material closes round it through a neck that is the
//     only path between the top half of the block and the bottom half;
//   * the port that breaks through the top face, with a counterbored face
//     round its mouth: a notch, not a hole, and the step at the counterbore
//     puts four more corners in a space 0.4 wide;
//   * an L-shaped milled cross gallery, entirely internal, six corners of which
//     five are reflex;
//   * four bolt holes, the lower two with 0.07 of metal between them and the
//     main gallery.
//
// The topology is five holes -- the L, and the four bolt holes -- against a
// notch that is not one, so the Euler characteristic is -4. The corner
// inventory, measured, is 11 convex and 7 reflex: four from the outline, six
// convex and two reflex from the port and its counterbore, and five reflex
// against one convex from the L. That sums to +1, so the boundary corners
// account for 1 of the -4 with the wrong sign and the field has to find -5
// elsewhere. It finds it entirely on the boundary -- 38 boundary singularities,
// 11 convex and 27 reflex, summing to -4 with the corners -- and ends with no
// internal singularities at all.
//
// The other thing to note is that every one of those corners is at exactly 90
// or 270 degrees, and that fourteen of the eighteen are on curves the
// silhouette gives no hint of. Whatever the layout ends up looking like, it is
// driven from the inside out.
//
//     np    triangles   corners        boundary sing.   raw -> collapsed
//     50        2813     11 / 23       38, sum -4        100 -> 65
//     75        5969     11 / 7        38, sum -4        100 -> 61
//    150       22198     11 / 7        38, sum -4         97 -> 68
//
// Every block comes out four-sided at all three sizes, every iso-line reaches
// the boundary, and the parameterization has no flipped triangles: of the five
// engineering models this is the one that goes through the whole pipeline clean
// and it is the natural regression case. The extra reflex corners at np=50 are
// the bolt-hole circles and the rounded gallery ends being read as polygons;
// they do not disturb the field, which produces the same 38 singularities at
// every size. data/meshes carries np=50.
SetFactory("OpenCASCADE");

Rectangle(1) = {0, 0, 0, 3.0, 1.6};

// Main gallery: drilled 0.26 across on the centreline y = 0.45, blind at both
// ends. The left end stops at x = 0.17, which is the neck the whole block
// hangs together by.
Disk(2) = {0.30, 0.45, 0, 0.13};
Disk(3) = {2.55, 0.45, 0, 0.13};
Rectangle(4) = {0.30, 0.32, 0, 2.25, 0.26};

// The port: 0.20 wide, up from the gallery and out through the top face, with
// a 0.40 wide counterbore for the fitting.
Rectangle(5) = {0.65, 0.45, 0, 0.20, 1.20};
Rectangle(6) = {0.55, 1.44, 0, 0.40, 0.21};

// L-shaped cross gallery, milled, so its corners are sharp and it does not
// reach any face.
Rectangle(7) = {2.00, 0.80, 0, 0.20, 0.45};
Rectangle(8) = {1.55, 1.05, 0, 0.65, 0.20};

// Bolt holes, 0.14 across, on a 0.20 inset from each corner.
Disk(9)  = {0.20, 0.16, 0, 0.07};
Disk(10) = {2.80, 0.16, 0, 0.07};
Disk(11) = {0.20, 1.40, 0, 0.07};
Disk(12) = {2.80, 1.40, 0, 0.07};

BooleanUnion(20) = { Surface{2}; Delete; }
                  { Surface{3}; Surface{4}; Surface{5}; Surface{6}; Delete; };
BooleanUnion(21) = { Surface{7}; Delete; }{ Surface{8}; Delete; };

BooleanDifference(30) = { Surface{1}; Delete; }
                       { Surface{20}; Surface{21}; Surface{9}; Surface{10};
                         Surface{11}; Surface{12}; Delete; };
