// A diseased artery in cross-section: the geometry of plaque-rupture stress
// analysis. The wall is three rings -- adventitia (mat 2, outermost, 0.15
// thick), media (mat 1, 0.13 thick, between the external and internal
// elastic laminae at r = 0.85 and 0.72) and a thickened intima of fibrous
// tissue (mat 3) -- round the lumen, which is a void. Inside the intima sits
// a lipid-rich necrotic core (mat 4): a crescent lying against the internal
// elastic lamina over 130 degrees of it, 0.20 deep, separated from the lumen
// only by the fibrous cap, 0.08 thick at its thinnest. A small plate of
// calcification (mat 5) sits in the lower shoulder of the plaque. The cap's
// thickness is the clinical variable (thin caps rupture), which is why these
// models exist.
//
// What this domain is for. The domain is an annulus: three smooth closed
// interfaces nested round a hole, each ring of Euler characteristic zero, so
// the field
// needs no singularity in them at all and a layout of the rings alone needs
// none either. Everything interesting is then at the plaque:
//
//  * The lipid core is a crescent whose two cusps lie ON an interface (the
//    internal elastic lamina), so each cusp is a triple junction -- intima,
//    lipid, media -- where two circles cross at 28.4 degrees: interface
//    directions 28 degrees apart, ILL-POSED, and an acute corner that is not
//    on dS.
//  * The fibrous cap is a thin layer (4 elements at NP = 100, h = 0.02 here)
//    between a material and a void.
//  * The intima is a region with two holes, one a void (the lumen) and one a
//    material (the calcification), and a boundary that is part circle, part
//    crescent.
//
// No full-circle material region: the rings are annuli, the lumen is an
// empty ellipse, and the calcification is an ellipse with its long axis at
// -20 degrees, which pins the singularities inside it.

// ---- generated geometry ----
Point(1) = {0.0000000000, 0.0000000000, 0, 1.0};   // centre of the three rings
Point(2) = {1.0000000000, 0.0000000000, 0, 1.0};
Point(3) = {0.0000000000, 1.0000000000, 0, 1.0};
Point(4) = {-1.0000000000, 0.0000000000, 0, 1.0};
Point(5) = {-0.0000000000, -1.0000000000, 0, 1.0};
Point(6) = {0.8500000000, 0.0000000000, 0, 1.0};
Point(7) = {0.0000000000, 0.8500000000, 0, 1.0};
Point(8) = {-0.8500000000, 0.0000000000, 0, 1.0};
Point(9) = {-0.0000000000, -0.8500000000, 0, 1.0};
Point(10) = {0.3042851485, -0.6525416067, 0, 1.0};   // lower cusp of the lipid core, on the IEL
Point(11) = {0.3042851485, 0.6525416067, 0, 1.0};   // upper cusp
Point(12) = {-0.6954665949, 0.1863497125, 0, 1.0};
Point(13) = {-0.0627521348, -0.7172601826, 0, 1.0};
Point(14) = {-0.5748329293, 0.0000000000, 0, 1.0};   // centre of the lipid core's inner arc
Point(15) = {0.0200000000, 0.0000000000, 0, 1.0};   // lumen centre
Point(16) = {0.4400000000, 0.0000000000, 0, 1.0};
Point(17) = {0.0200000000, 0.3300000000, 0, 1.0};
Point(18) = {-0.4000000000, 0.0000000000, 0, 1.0};
Point(19) = {0.0200000000, -0.3300000000, 0, 1.0};
Point(20) = {0.1000000000, -0.5200000000, 0, 1.0};   // calcification centre
Point(21) = {0.1845723359, -0.5507818129, 0, 1.0};
Point(22) = {0.1119707050, -0.4871107583, 0, 1.0};
Point(23) = {0.0154276641, -0.4892181871, 0, 1.0};
Point(24) = {0.0880292950, -0.5528892417, 0, 1.0};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};
Circle(9) = {10, 1, 11};
Circle(10) = {11, 1, 12};
Circle(11) = {12, 1, 13};
Circle(12) = {13, 1, 10};
Circle(13) = {10, 14, 11};
Ellipse(14) = {16, 15, 16, 17};
Ellipse(15) = {17, 15, 16, 18};
Ellipse(16) = {18, 15, 16, 19};
Ellipse(17) = {19, 15, 16, 16};
Ellipse(18) = {21, 20, 21, 22};
Ellipse(19) = {22, 20, 21, 23};
Ellipse(20) = {23, 20, 21, 24};
Ellipse(21) = {24, 20, 21, 21};

Curve Loop(11) = {1, 2, 3, 4};
Curve Loop(12) = {5, 6, 7, 8};
Curve Loop(13) = {9, 10, 11, 12};
Curve Loop(14) = {14, 15, 16, 17};
Curve Loop(15) = {18, 19, 20, 21};
Curve Loop(16) = {10, 11, 12, 13};
Curve Loop(17) = {9, -13};
Plane Surface(1) = {12, 13};       // media
Plane Surface(2) = {11, 12};       // adventitia
Plane Surface(3) = {16, 14, 15};   // intima (fibrous tissue), lumen and calcification cut out
Plane Surface(4) = {17};           // lipid core
Plane Surface(5) = {15};           // calcification

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
