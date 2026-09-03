// Involute spur gear, 14 teeth, on a keyed hub: a machined part, drawn the way
// it is actually defined rather than as a decorative cog.
//
// The flanks are true involutes of the base circle, the tips and roots are arcs
// of the addendum and dedendum circles, and the bore carries a keyway. What
// that buys is a boundary that alternates convex and reflex corners fourteen
// times over on a curved carrier -- geom028's star does the same alternation
// but on straight edges and only six times -- and 60 corners is the most in
// data/geometry by a factor of three.
//
// The index cancels exactly, and that is the point of the model. Measured, the
// mesh finds 30 corners it rounds to a quarter turn: 28 tooth tips at 123
// degrees, where the involute meets the tip land, and the two at the mouth of
// the keyway at 103. It finds 30 it rounds to three quarters: 28 roots at 269.5
// and the two at the bottom of the keyway at 270. Thirty at +1/4 against thirty
// at -1/4 is zero, the bore makes the model an annulus so the Euler
// characteristic is zero too, and the field therefore has no index to place
// anywhere. It does not place any: DualMBO finds no singularities to begin with
// and UMBER ends with none, at every mesh size tried. Nothing else in the set
// combines that many corners with that little to do.
//
// What breaks is the layout, and it breaks badly enough that this is the
// failing case of geom031-035 the way geom024 is the failing case of the
// earlier group. The 60 corners launch 60 motorcycles into a field that is very
// nearly a rotation about the centre, so nothing stops them: they run round the
// annulus and cross each other thousands of times, and the flood then has to
// call almost every crossing a block.
//
//     np    triangles   corners   motorcycles   raw layout      after collapse
//     50        3873     30 / 30      60         7968 blocks     7878, 4 rays short
//     75        8129     30 / 30      60        19793 blocks     3992, 2 not 4-sided
//    150       28609     30 / 30      60        31209 blocks    31126, 6 rays short
//
// The corner inventory is identical at all three and the field is clean at all
// three; it is only the tracing that is hopeless. Refining makes it worse, not
// better, which is the useful signal: the problem is not resolution, it is that
// a nearly rotational field on an annulus gives a motorcycle nowhere to stop.
//
// The teeth are also the thinnest repeated feature here that is not a slot:
// about 0.065 of tip land against a 0.1 module, so np has to reach 150 before a
// tooth is more than two elements wide. data/meshes carries np=50.
SetFactory("OpenCASCADE");

z  = 14;             // teeth
m  = 0.1;            // module: pitch diameter / number of teeth
a  = 20*Pi/180;      // pressure angle, the standard 20 degrees
Nf = 12;             // points per involute flank

rp = m*z/2;          // pitch radius
rb = rp*Cos(a);      // base radius: the involute is unwound from this circle
ra = rp + m;         // addendum (tip) radius
rf = rp - 1.25*m;    // dedendum (root) radius

inv_a = Tan(a) - a;              // involute function at the pressure angle
psi   = Pi/(2*z) + inv_a;        // half tooth angle, measured at the base circle
tmax  = Sqrt((ra/rb)*(ra/rb) - 1);  // involute parameter at the tip

Point(1) = {0, 0, 0, 1.0};       // centre, shared by every arc

// One tooth per pass: right flank base-to-tip, tip arc, left flank tip-to-base,
// and the root arc that runs on to the next tooth. The involute is
// r = rb*Sqrt(1+t^2) at polar angle t - Atan(t), so the flank feet sit on the
// base circle at +-psi and the tooth is symmetric about its own centre line.
For k In {0:z-1}
  phi = 2*Pi*k/z;
  b = 1000 + 100*k;

  For i In {0:Nf}
    t  = tmax*i/Nf;
    r  = rb*Sqrt(1 + t*t);
    th = phi - psi + (t - Atan(t));
    Point(b + i) = {r*Cos(th), r*Sin(th), 0, 1.0};
  EndFor

  For i In {0:Nf}
    t  = tmax*(Nf - i)/Nf;
    r  = rb*Sqrt(1 + t*t);
    th = phi + psi - (t - Atan(t));
    Point(b + 20 + i) = {r*Cos(th), r*Sin(th), 0, 1.0};
  EndFor

  // Feet of the two undercuts, on the root circle directly below the flanks.
  th = phi - psi;
  Point(b + 50) = {rf*Cos(th), rf*Sin(th), 0, 1.0};
  th = phi + psi;
  Point(b + 51) = {rf*Cos(th), rf*Sin(th), 0, 1.0};
EndFor

cl[] = {};
For k In {0:z-1}
  b = 1000 + 100*k;
  c = 2000 + 10*k;
  kn = k + 1;
  If (k == z-1)
    kn = 0;
  EndIf

  Line(c + 1)   = {b + 50, b};                    // undercut, root circle to base circle
  Spline(c + 2) = {b : b + Nf};                   // right flank
  Circle(c + 3) = {b + Nf, 1, b + 20};            // tip land
  Spline(c + 4) = {b + 20 : b + 20 + Nf};         // left flank
  Line(c + 5)   = {b + 20 + Nf, b + 51};          // undercut back down
  Circle(c + 6) = {b + 51, 1, 1000 + 100*kn + 50};// root, on to the next tooth

  cl[] += {c + 1, c + 2, c + 3, c + 4, c + 5, c + 6};
EndFor

Curve Loop(1) = cl[];
Plane Surface(1) = {1};

// Bore and keyway, cut as one: a 0.5 diameter bore with a 0.1 wide key seat
// 0.045 deep, which is roughly what a standard key for this size of shaft asks
// for and, more usefully here, the single feature that breaks the symmetry.
Disk(2) = {0, 0, 0, 0.25};
Rectangle(3) = {-0.05, 0, 0, 0.10, 0.295};
BooleanUnion(4) = { Surface{2}; Delete; }{ Surface{3}; Delete; };
BooleanDifference(5) = { Surface{1}; Delete; }{ Surface{4}; Delete; };
