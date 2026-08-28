// Dissimilar-metal butt weld: a single-V groove joining a ferritic steel plate
// (mat 1, left) to an austenitic stainless plate (mat 4, right), with a
// buttering layer (mat 2) deposited on the ferritic groove face and the weld
// filler (mat 3) filling the groove. This is the standard remedy for welding
// two alloys that must not be fused directly, and it is the reason the domain
// has four materials rather than three: the buttering exists precisely to keep
// mat 1 and mat 4 from ever sharing an interface.
//
// What this domain is for. The root point (curve endpoints at 0.000, -0.140) is a
// triple junction where the buttering, the weld and the austenitic plate meet,
// at interface directions 118, 58 and 90 degrees. Reduced mod 90 those are
// 0, 28 and 58 -- not congruent, so no single cross can be tangent to all
// three, and a vertex-based field has no well-posed value there. Unlike the
// synthetic junction family (geom004, geom006, geom007) the angles here are not
// chosen to make that happen: they are the groove geometry a welding procedure
// specification would call for, so the ill-posedness cannot be dismissed as
// adversarial.
//
// The buttering is also the thin feature. At 0.11 thick against a 0.80 plate it
// is roughly three elements across at NP=70, and its inner boundary carries a
// kink at (-0.110, -0.167) where the layer turns off the root face and onto the
// groove face -- an interface that is continuous but not smooth, which R1
// alignment has to follow around the corner.
//
// The boundary junctions are deliberately mixed. Where the two vertical root
// interfaces meet the bottom edge the incident directions are 0 and 90, so
// those two are benign; the three on the top edge meet it at 118, 118 and 58
// degrees and are not.

Point(1) = {-1.100000, -0.400000, 0, 1.0};
Point(2) = {-0.110000, -0.400000, 0, 1.0};
Point(3) = {0.000000, -0.400000, 0, 1.0};
Point(4) = {1.100000, -0.400000, 0, 1.0};
Point(5) = {1.100000, 0.400000, 0, 1.0};
Point(6) = {0.337429, 0.400000, 0, 1.0};
Point(7) = {-0.287123, 0.400000, 0, 1.0};
Point(8) = {-0.411706, 0.400000, 0, 1.0};
Point(9) = {-1.100000, 0.400000, 0, 1.0};
Point(10) = {0.000000, -0.140000, 0, 1.0};
Point(11) = {-0.110000, -0.167426, 0, 1.0};

// Outer boundary, counter-clockwise from the bottom-left corner. It is split at
// every point where a material interface lands on it so that each piece belongs
// to exactly one surface.
Line(1) = {1, 2};  // bottom, ferritic
Line(2) = {2, 3};  // bottom, buttering
Line(3) = {3, 4};  // bottom, austenitic
Line(4) = {4, 5};  // right wall, austenitic
Line(5) = {5, 6};  // top, austenitic
Line(6) = {6, 7};  // top, weld
Line(7) = {7, 8};  // top, buttering
Line(8) = {8, 9};  // top, ferritic
Line(9) = {9, 1};  // left wall, ferritic

// Material interfaces. Each is written once and used twice, forward in one
// Curve Loop and reversed in the other, so Gmsh meshes it a single time and the
// two surfaces share its nodes.
Line(10) = {3, 10};  // root face, buttering | austenitic
Line(11) = {10, 6};  // groove face, weld | austenitic
Line(12) = {10, 7};  // fusion line, buttering | weld
Line(13) = {2, 11};  // root offset, ferritic | buttering
Line(14) = {11, 8};  // groove offset, ferritic | buttering

Curve Loop(1) = {1, 13, 14, 8, 9};           // ferritic plate
Plane Surface(1) = {1};
Curve Loop(2) = {2, 10, 12, 7, -14, -13};    // buttering layer
Plane Surface(2) = {2};
Curve Loop(3) = {11, 6, -12};                // weld metal
Plane Surface(3) = {3};
Curve Loop(4) = {3, 4, 5, -11, -10};         // austenitic plate
Plane Surface(4) = {4};

Physical Surface(1) = {1};   // ferritic plate
Physical Surface(2) = {2};   // buttering layer
Physical Surface(3) = {3};   // weld metal
Physical Surface(4) = {4};   // austenitic plate
