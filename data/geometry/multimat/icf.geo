// Indirect-drive ICF target: a layered capsule inside a cylindrical gold
// hohlraum, drawn in the half plane r >= 0 with the symmetry axis on x = 0.
//
// This is geom010's concentric-shell capsule put back into the machine that
// drives it. The capsule (mat 1 DT vapour, mat 2 DT ice, mat 3 ablator) hangs
// at the centre of a helium-filled can (mat 4) whose wall is gold (mat 5); the
// laser entrance hole at each end of the can is spanned by a thin polyimide
// window (mat 6). Six materials, and the length scales span two decades: the
// 0.05-deep ice shell and the 0.06-thick wall against a 2.32-tall domain.
//
// What this domain is for, beyond geom010. Three things it adds:
//
//  * Two interior triple junctions, at (0.280, +/-1.100), where the window, the
//    gold lip and the fill gas meet. The incident interface directions are 0
//    (window | gas), 90 (gold | window) and 0 (gold | gas); reduced mod 90 they
//    are all 0, so both are WELL-POSED. They are junctions the hardware
//    dictates, not junctions chosen to be benign, which is what makes them
//    worth having in the benign column next to geom003/geom005/geom014.
//
//  * The window is one material on two disconnected surfaces
//    (Physical Surface(6) = {5,6}), as geom014's underfill is.
//
//  * Curved and straight interfaces in one domain with nothing shared between
//    them: R1 alignment has to follow three nested arcs near the axis and, four
//    capsule radii away, an axis-parallel wall and two axis-normal end caps.
//    A single global cross field has to do both at once.
//
// Every boundary junction is benign by construction. The three capsule arcs
// cross the axis perpendicular, and the window faces meet the axis at (0,
// +/-1.100) and (0, +/-1.160) at right angles.
//
// No full circles: the capsule shells are half disks cut on the axis, so the
// interior singularities are pinned by the cut rather than free to rotate.

Point(1) = {0, 0, 0, 1.0};  // capsule centre, on the axis

// Capsule shells: DT vapour 0.12, DT ice to 0.17, ablator to 0.24. The radius
// ratio to the can (0.24 : 0.60) and its aspect (1.10 : 0.60) are the NIF
// proportions; the shells are thin for the same reason the real ones are.
Point(2) = {0, -0.12, 0, 1.0};
Point(3) = {0.12, 0, 0, 1.0};
Point(4) = {0, 0.12, 0, 1.0};
Point(5) = {0, -0.17, 0, 1.0};
Point(6) = {0.17, 0, 0, 1.0};
Point(7) = {0, 0.17, 0, 1.0};
Point(8) = {0, -0.24, 0, 1.0};
Point(9) = {0.24, 0, 0, 1.0};
Point(10) = {0, 0.24, 0, 1.0};

// Hohlraum cavity: inner radius 0.60, inner half-length 1.10, laser entrance
// hole of radius 0.28 at each end.
Point(11) = {0, 1.10, 0, 1.0};
Point(12) = {0.28, 1.10, 0, 1.0};   // triple junction: gas | gold | window
Point(13) = {0.60, 1.10, 0, 1.0};
Point(14) = {0.60, -1.10, 0, 1.0};
Point(15) = {0.28, -1.10, 0, 1.0};  // triple junction: gas | gold | window
Point(16) = {0, -1.10, 0, 1.0};

// Outer surface of the can: wall and end caps both 0.06 thick, and the windows
// sit flush in the holes.
Point(17) = {0, 1.16, 0, 1.0};
Point(18) = {0.28, 1.16, 0, 1.0};
Point(19) = {0.66, 1.16, 0, 1.0};
Point(20) = {0.66, -1.16, 0, 1.0};
Point(21) = {0.28, -1.16, 0, 1.0};
Point(22) = {0, -1.16, 0, 1.0};

// Capsule interfaces, two quarter arcs per radius (a GEO Circle may not span
// 180 degrees).
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {5, 1, 6};
Circle(4) = {6, 1, 7};
Circle(5) = {8, 1, 9};
Circle(6) = {9, 1, 10};

// Axis segments of the capsule, one pair per shell.
Line(7) = {4, 1};
Line(8) = {1, 2};
Line(9) = {7, 4};
Line(10) = {2, 5};
Line(11) = {10, 7};
Line(12) = {5, 8};

// Cavity boundary, walking up the axis, around the can and back down. Each
// interface is written once and reused reversed, so the surfaces share nodes.
Line(13) = {10, 11};  // axis, gas above the capsule
Line(14) = {11, 12};  // top window inner face, window | gas
Line(15) = {12, 13};  // top lip, gold | gas
Line(16) = {13, 14};  // cavity wall, gold | gas
Line(17) = {14, 15};  // bottom lip, gold | gas
Line(18) = {15, 16};  // bottom window inner face, window | gas
Line(19) = {16, 8};   // axis, gas below the capsule

// Top window.
Line(20) = {12, 18};  // hole edge, gold | window
Line(21) = {18, 17};  // outer face
Line(22) = {17, 11};  // axis

// Bottom window.
Line(23) = {15, 21};  // hole edge, gold | window
Line(24) = {21, 22};  // outer face
Line(25) = {22, 16};  // axis

// Outer surface of the gold.
Line(26) = {18, 19};  // top end cap
Line(27) = {19, 20};  // wall
Line(28) = {20, 21};  // bottom end cap

Curve Loop(1) = {1, 2, 7, 8};                        // DT vapour
Plane Surface(1) = {1};
Curve Loop(2) = {3, 4, 9, -2, -1, 10};               // DT ice
Plane Surface(2) = {2};
Curve Loop(3) = {5, 6, 11, -4, -3, 12};              // ablator
Plane Surface(3) = {3};
Curve Loop(4) = {-6, -5, -19, -18, -17, -16, -15, -14, -13};  // helium fill gas
Plane Surface(4) = {4};
Curve Loop(5) = {14, 20, 21, 22};                    // top window
Plane Surface(5) = {5};
Curve Loop(6) = {18, -25, -24, -23};                   // bottom window
Plane Surface(6) = {6};
Curve Loop(7) = {-20, 15, 16, 17, 23, -28, -27, -26};  // gold wall
Plane Surface(7) = {7};

//+
Delete {
  Surface{7}; 
}
//+
Delete {
  Surface{4}; 
}
//+
Delete {
  Curve{16}; 
}
//+
Delete {
  Curve{27}; 
}
//+
Delete {
  Curve{15}; 
}
//+
Delete {
  Curve{26}; 
}
//+
Delete {
  Curve{17}; 
}
//+
Delete {
  Curve{28}; 
}
//+
Point(23) = {0.5, 1.16, 0, 1.0};
//+
Point(24) = {0.5, 1.05, 0, 1.0};
//+
Point(25) = {0.66, 1.05, 0, 1.0};
//+
Point(26) = {0.5, -1.16, 0, 1.0};
//+
Point(27) = {0.5, -1.04, 0, 1.0};
//+
Point(28) = {0.66, -1.04, 0, 1.0};
//+
Ellipse(26) = {23, 24, 25, 25};
//+
Ellipse(27) = {26, 27, 28, 28};
//+
Line(28) = {18, 23};
//+
Line(29) = {21, 26};
//+
Line(30) = {28, 25};
//+
Point(29) = {0.5, 1.1, 0, 1.0};
//+
Point(30) = {0.59, 1.05, 0, 1.0};
//+
Point(31) = {0.5, -1.1, 0, 1.0};
//+
Point(32) = {0.59, -1.04, 0, 1.0};
//+
Ellipse(31) = {29, 24, 30, 30};
//+
Ellipse(32) = {31, 27, 32, 32};
//+
Line(33) = {15, 31};
//+
Line(34) = {12, 29};
//+
Line(35) = {32, 30};
//+
Curve Loop(8) = {-20, 34, 31, -35, -32, -33, 23, 29, 27, 30, -26, -28};  // gold wall, CCW
//+
Plane Surface(7) = {8};
//+
Curve Loop(9) = {-6, -5, -19, -18, 33, 32, 35, -31, -34, -14, -13};  // helium fill gas, CCW
//+
Plane Surface(8) = {9};

// Physical groups. Declared after the rounding edits above, which deleted the
// original square-cornered gas (Surface 4) and gold (Surface 7) surfaces and
// rebuilt them as Surface 8 and Surface 7 respectively.
Physical Surface(1) = {1};      // DT vapour
Physical Surface(2) = {2};      // DT ice
Physical Surface(3) = {3};      // ablator
Physical Surface(4) = {8};      // helium fill gas
Physical Surface(5) = {7};      // gold wall
Physical Surface(6) = {5, 6};   // LEH windows
