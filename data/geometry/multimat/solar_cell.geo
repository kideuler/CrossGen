// A textured crystalline-silicon solar cell in cross-section. The front of
// the silicon wafer (mat 1) is etched into pyramids, drawn here as a regular
// row of four; a silicon-nitride anti-reflection coating (mat 2) of uniform
// thickness follows the texture; a screen-printed silver finger (mat 5)
// sits across one valley, filling it; the module's encapsulant (mat 3) fills
// the rest of the texture up to a flat face, under the cover glass (mat 4);
// an aluminium back contact (mat 6) covers the rear.
//
// What this domain is for. Kinks at angles set by crystallography. An
// alkaline etch exposes the {111} planes, so every facet stands at
// arctan(sqrt 2) = 54.74 degrees to the wafer and every pyramid apex is 70.53
// degrees -- nobody chose these numbers, and mod 90 they are 54.74 and 35.26,
// about as far from a quarter turn as angles get. The coating is a thin band
// (3 elements at NP = 100, h = 0.02 here) with SEVEN kinks in each of its two
// interfaces, alternately convex and reflex, and because a conformal offset
// of a zigzag by a uniform normal thickness is the zigzag translated
// straight up (by 0.06 sqrt 3), the band's two faces are exactly parallel.
// A layout either follows the band through every kink or cuts across it.
//
// Around that: the finger's two feet are triple junctions (coating, finger,
// encapsulant) where its rounded top meets a facet, directions about
// {71, 125} -- mod 90 {71, 35}, ILL-POSED; at both side walls the texture is
// in a valley, so the coating ends in 35-degree wedges against the walls; and
// the flat layers above and below are there to be easy.
//
// Proportions are a caricature (pyramids are ~5 microns, the wafer 180, the
// coating 75 nanometres): the coating is thickened and the wafer thinned so
// that one domain holds both, but every angle is exact. Nothing is a circle.

H  = 0.25 * Sqrt(2);        // pyramid height: half-base 0.25 times tan 54.74
dv = 0.06 * Sqrt(3);        // coating's vertical offset for 0.06 normal thickness

Point(1) = {0, -0.72, 0, 1.0};
Point(2) = {2, -0.72, 0, 1.0};
Point(3) = {2, -0.57, 0, 1.0};
Point(4) = {0, -0.57, 0, 1.0};

// The silicon surface: valleys at x = 0, 0.5, 1, 1.5, 2 and apexes between.
For k In {0:8}
  Point(11 + k) = {0.25 * k, H * Fmod(k, 2), 0, 1.0};
EndFor

// The coating's outer face: the same zigzag, dv higher, with the finger's
// feet inserted where they sit on the facets either side of the x = 1 valley.
Point(21) = {0.00, dv,     0, 1.0};
Point(22) = {0.25, H + dv, 0, 1.0};
Point(23) = {0.50, dv,     0, 1.0};
Point(24) = {0.75, H + dv, 0, 1.0};
Point(25) = {0.82, 0.72 * H + dv, 0, 1.0};   // finger's left foot
Point(26) = {1.00, dv,     0, 1.0};
Point(27) = {1.18, 0.72 * H + dv, 0, 1.0};   // finger's right foot
Point(28) = {1.25, H + dv, 0, 1.0};
Point(29) = {1.50, dv,     0, 1.0};
Point(30) = {1.75, H + dv, 0, 1.0};
Point(31) = {2.00, dv,     0, 1.0};

// The finger's rounded top.
Point(41) = {0.855, 0.460, 0, 1.0};
Point(42) = {0.920, 0.525, 0, 1.0};
Point(43) = {1.000, 0.545, 0, 1.0};
Point(44) = {1.080, 0.525, 0, 1.0};
Point(45) = {1.145, 0.460, 0, 1.0};

// Encapsulant | glass, and the top.
Point(51) = {2, 0.62, 0, 1.0};
Point(52) = {0, 0.62, 0, 1.0};
Point(53) = {2, 0.87, 0, 1.0};
Point(54) = {0, 0.87, 0, 1.0};

// Back contact, and the walls and top.
Line(1)  = {1, 2};
Line(2)  = {2, 3};
Line(3)  = {3, 4};       // aluminium | silicon
Line(4)  = {4, 1};
Line(5)  = {3, 19};
Line(6)  = {19, 31};
Line(7)  = {31, 51};
Line(8)  = {51, 53};
Line(9)  = {53, 54};
Line(10) = {54, 52};
Line(11) = {52, 21};
Line(12) = {21, 11};
Line(13) = {11, 4};

// Silicon | coating: the texture, facet by facet.
For k In {0:7}
  Line(21 + k) = {11 + k, 12 + k};
EndFor

// Coating's outer face: encapsulant above, except under the finger.
For k In {0:9}
  Line(31 + k) = {21 + k, 22 + k};
EndFor

Spline(41) = {25, 41, 42, 43, 44, 45, 27};   // finger | encapsulant
Line(42)   = {52, 51};                       // encapsulant | glass

Curve Loop(1) = {-3, 5, -28, -27, -26, -25, -24, -23, -22, -21, 13};
Plane Surface(1) = {1};                      // silicon
Curve Loop(2) = {21, 22, 23, 24, 25, 26, 27, 28, 6,
                 -40, -39, -38, -37, -36, -35, -34, -33, -32, -31, 12};
Plane Surface(2) = {2};                      // anti-reflection coating
Curve Loop(3) = {31, 32, 33, 34, 41, 37, 38, 39, 40, 7, -42, 11};
Plane Surface(3) = {3};                      // encapsulant
Curve Loop(4) = {42, 8, 9, 10};
Plane Surface(4) = {4};                      // glass
Curve Loop(5) = {35, 36, -41};
Plane Surface(5) = {5};                      // silver finger
Curve Loop(6) = {1, 2, 3, 4};
Plane Surface(6) = {6};                      // aluminium back contact

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};
