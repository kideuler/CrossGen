// Particulate dispersion: ten disks (mat 2) scattered through a 1 x 0.5
// matrix block (mat 1). The disks are the inclusions of a two-phase
// composite -- bubbles in a melt, voids in a casting, second-phase particles
// in an alloy. Radii run 0.032 to 0.057 against a domain 1 wide, and the
// inclusions take up 13% of the area.
//
// The block is written in unit coordinates, x in [0,1] and y in [0,0.5], so
// the numbers in this file are already the fractions of the domain that
// Mesh2Dgmsh's normalisation would report anyway: it rescales to a bounding
// box of max dimension 2, which here is exactly a factor of 2 on everything.
//
// What this domain is for. Two things nothing else in the corpus does:
//
//  * Count. Mat 2 is ONE material on TEN disconnected surfaces
//    (Physical Surface(2) = {300:309}), against geom014's six and icf's two,
//    and the matrix is one surface with ten holes punched in it, so its
//    Euler characteristic is 1 - 10 = -9 and it needs forty -1/4 cones to
//    the disks' forty +1/4.
//
//  * Ten independent copies of the same local problem. Every inclusion poses
//    the identical question to the field and the layout, at ten places and
//    at radii spanning a factor of 1.8, so whether the answers agree is
//    directly readable.
//
// There are no junctions of any kind: no interface touches dS and no two
// interfaces touch, so every interface is a smooth closed curve and the
// junction table has nothing to say about this model. It isolates the
// many-inclusions question from the well-posedness question entirely.
//
// CAVEAT, stated because the rest of data/geometry/multimat obeys the
// opposite rule. These ARE full circles, which the corpus otherwise forbids
// (see geom009's comment): a disk is rotationally symmetric, so the four
// +1/4 interior singularities its Gauss-Bonnet deficit demands have no
// preferred angular position and land wherever arithmetic noise puts them.
// On this model that is deliberate rather than an oversight -- ten
// independent copies of exactly that indeterminacy is the point, and whether
// the ten disks pick consistent angles is itself the measurement. Do not
// copy the pattern into a domain meant to be a control.

// Matrix block, 1 x 0.5.
Point(1) = {0,   0,   0, 1.0};
Point(2) = {1,   0,   0, 1.0};
Point(3) = {1,   0.5, 0, 1.0};
Point(4) = {0,   0.5, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};

// Ten inclusions, interspersed over x in (0.125, 0.875) and staggered in y.
// Measured clearances: the closest two disks are 0.031 apart rim to rim and
// the closest wall approach is 0.043, both comfortably above the mesh size
// the corpus meshes at (2/NP before normalisation, so 0.01 here at NP = 100),
// which is what keeps the matrix between them more than one element thick.
cx[] = {0.180, 0.244, 0.325, 0.406, 0.475, 0.556, 0.638, 0.706, 0.788, 0.828};
cy[] = {0.130, 0.356, 0.212, 0.388, 0.106, 0.275, 0.419, 0.169, 0.331, 0.119};
cr[] = {0.052, 0.035, 0.044, 0.057, 0.032, 0.049, 0.038, 0.055, 0.040, 0.046};

// The GEO kernel's Circle may not span 180 degrees or more, so each disk is
// four quarter arcs through its axis points.
//
// Arc discretisation is left to Mesh2Dgmsh's uniform size field (2/NP after
// normalisation) rather than pinned with Transfinite Curve, so the disk
// boundaries refine along with the rest of the mesh as NP increases. This
// trades away NP-independent interface discretisation (arc segment count
// used to be fixed at 24 per disk regardless of NP) for disks that actually
// get finer when asked to.
holes[] = {};
For k In {0:9}
  p = 100 + 5*k;
  c = 200 + 4*k;
  Point(p)     = {cx[k],         cy[k],         0, 1.0};  // centre
  Point(p + 1) = {cx[k] + cr[k], cy[k],         0, 1.0};
  Point(p + 2) = {cx[k],         cy[k] + cr[k], 0, 1.0};
  Point(p + 3) = {cx[k] - cr[k], cy[k],         0, 1.0};
  Point(p + 4) = {cx[k],         cy[k] - cr[k], 0, 1.0};
  Circle(c)     = {p + 1, p, p + 2};
  Circle(c + 1) = {p + 2, p, p + 3};
  Circle(c + 2) = {p + 3, p, p + 4};
  Circle(c + 3) = {p + 4, p, p + 1};
  Curve Loop(300 + k) = {c, c + 1, c + 2, c + 3};
  Plane Surface(300 + k) = {300 + k};
  holes[] += 300 + k;
EndFor

// The matrix is the block with the ten loops punched out of it; each hole
// loop is the disk's own Curve Loop, reused rather than redeclared, so the
// interfaces are conformal on both sides.
Plane Surface(1) = {1, holes[]};

Physical Surface(1) = {1};
Physical Surface(2) = {300:309};