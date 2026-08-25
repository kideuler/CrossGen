// Internally cooled turbine blade, in section: the metal, not the gas.
//
// geom029 and geom030 put an aerofoil in a flow domain, where the section is a
// hole and the trailing edge is a 344-degree reflex corner. This is the same
// object from the other side. The domain is the casting itself, so the trailing
// edge is now a 16-degree convex corner -- which rounds to a quarter turn and
// sits at +1/4, launching nothing -- and the cooling passages that a solver
// would never see are holes.
//
// The section is a NACA 4-digit shape with 6% camber and 20% thickness, which
// is about right for a cooled high-pressure blade and, more to the point here,
// thick enough to put four passages inside: two round ones near the leading
// edge, a mid-chord racetrack, and the thin slot that in the real part feeds
// trailing-edge ejection. The ribs between them are 0.046 and 0.048 of chord,
// and the land between the slot and the suction surface is 0.022.
//
// What makes it a hard model is that the boundary carries almost no index while
// the topology demands a great deal. Four holes put the Euler characteristic at
// -3. Against that the boundary offers, at np=150 and above, exactly one
// corner: the trailing edge. Its wedge is 27 degrees on the geometry and the
// mesh reads it as 48 at np=150 and 59 at np=250 -- both round to a quarter
// turn, so it is +1/4 whatever the resolution. The leading edge is not a corner
// here at all: its radius is 4.4% of chord, twenty times geom029's, and a
// uniform mesh resolves it as the curve it is. The passages are bounded by
// circles and tangent lines and contribute nothing.
//
// So one corner, +1/4, against a Euler characteristic of -3. The field puts the
// missing -13/4 on the boundary rather than in the interior: 20 boundary
// singularities, 4 convex and 16 reflex, summing to -3 against the corners, and
// UMBER ends with no internal singularities. Where they land is the thing to
// look at, because the passages are placed the way a casting is placed and not
// the way a test case is -- not on a lattice, not the same size, not symmetric
// about anything -- so nothing about the outline predicts them.
//
//     np    triangles   corners       boundary sing.   layout
//     50         742     1 / 4        22, sum -3        55 -> 44 blocks
//     75        1575     1 / 3        20, sum -3        49 -> 25 blocks
//    150        5994     1 / 0        20, sum -3        52 -> 33 blocks
//    250       16267     1 / 0
//
// Four blocks come out not four-sided at np=75 and 150, two at np=50, and chord
// collapse does not make it worse. The corner column is the interesting one:
// below np=150 the ends of the racetrack and of the trailing-edge slot are so
// few elements round that the mesh reads each as a corner of 225 to 270
// degrees -- three of them at np=75, four at np=50 -- which is index the
// geometry does not have. They disappear at np=150. data/meshes carries np=50,
// so the shipped mesh is on the wrong side of that.
SetFactory("OpenCASCADE");

N = 70;        // points per surface
mc = 0.06;     // maximum camber, fraction of chord
pc = 0.40;     // its position
tt = 0.20;     // thickness

// NACA 4-digit: thickness distribution laid off normal to the camber line. The
// last thickness coefficient is -0.1036, which closes the trailing edge to a
// point rather than leaving the 0.0021 gap the original -0.1015 gives.
Macro Section
  yt = 5*tt*( 0.2969*Sqrt(xc) - 0.1260*xc - 0.3516*xc*xc
            + 0.2843*xc*xc*xc - 0.1036*xc*xc*xc*xc );
  If (xc < pc)
    yc = mc/(pc*pc) * (2*pc*xc - xc*xc);
    dy = 2*mc/(pc*pc) * (pc - xc);
  Else
    yc = mc/((1-pc)*(1-pc)) * ((1 - 2*pc) + 2*pc*xc - xc*xc);
    dy = 2*mc/((1-pc)*(1-pc)) * (pc - xc);
  EndIf
  th = Atan(dy);
Return

// Suction side, leading edge to trailing edge, cosine spaced.
For i In {0:N}
  xc = 0.5*(1 - Cos(Pi*i/N));
  Call Section;
  Point(1000 + i) = {xc - yt*Sin(th), yc + yt*Cos(th), 0, 1.0};
EndFor

// Pressure side, back again, sharing the two tips.
For i In {1:N-1}
  xc = 0.5*(1 - Cos(Pi*(N-i)/N));
  Call Section;
  Point(2000 + i) = {xc + yt*Sin(th), yc - yt*Cos(th), 0, 1.0};
EndFor

Spline(1) = {1000:1000+N, 2001:2000+N-1, 1000};
Curve Loop(1) = {1};
Plane Surface(1) = {1};

// Leading-edge passages, on the camber line.
Disk(2) = {0.22, 0.0479, 0, 0.055};
Disk(3) = {0.42, 0.0599, 0, 0.050};

// Mid-chord racetrack: two ends and the straight between them.
Disk(4) = {0.55, 0.054, 0, 0.032};
Disk(5) = {0.63, 0.054, 0, 0.032};
Rectangle(6) = {0.55, 0.022, 0, 0.08, 0.064};

// Trailing-edge slot, the same construction an eighth the width.
Disk(7) = {0.72, 0.035, 0, 0.012};
Disk(8) = {0.82, 0.035, 0, 0.012};
Rectangle(9) = {0.72, 0.023, 0, 0.10, 0.024};

BooleanUnion(10) = { Surface{4}; Delete; }{ Surface{5}; Surface{6}; Delete; };
BooleanUnion(11) = { Surface{7}; Delete; }{ Surface{8}; Surface{9}; Delete; };
BooleanDifference(12) = { Surface{1}; Delete; }
                       { Surface{2}; Surface{3}; Surface{10}; Surface{11}; Delete; };
