// NACA 0012 in a wind tunnel: the external-flow domain of a fluid simulation.
//
// Here for the trailing edge. Closed to a point, the two surfaces meet at about
// 16 degrees, so measured in the fluid -- which is the side the domain is on --
// the interior angle is around 344 degrees, and Table 1 puts that in its last
// band: index -1/2, the only one that launches three rays. Nothing else in
// data/geometry reaches that band; the reflex corners of geom017 to geom028 are
// all 270-degree ones at -1/4, and a sharp trailing edge in a flow domain is
// the natural way a real model produces the other.
//
// The rest of it is the usual CFD arrangement and tests other things: the model
// is an annulus, so the field carries index it cannot put on the boundary; the
// leading edge is the tightest curvature in the set; and the run downstream is
// a long, almost parallel channel of the kind that produces very long chords in
// the layout.
//
// Mesh size. Mesh2Dgmsh sizes elements uniformly, so every chord of empty
// tunnel is paid for in elements across the section. Measured:
//
//     np    triangles   section    leading edge   layout
//     75        6666     2.5 el     277 deg        13 components, all four-sided
//    150       26172     5.0 el     248 deg        13 components, all four-sided
//    250       71952     8.3 el     227 deg        13 components, all four-sided
//
// The layout is the same at all three, which is the point: the coarse mesh is
// not lying, it just does not look like an aerofoil. data/meshes carries the
// np=150 one, and make_meshes.sh's default of 75 reproduces the same layout.
//
// The leading edge is the one thing a uniform mesh cannot do. Its radius is
// 1.1% of the chord, well under an element at any size worth meshing, so the
// nose is read as a corner rather than as smooth boundary and its angle keeps
// moving as the mesh is refined -- it is heading for 180 degrees and will not
// get there. The trailing edge, being a corner the geometry really has, sits in
// the -1/2 band at every resolution.

SetFactory("OpenCASCADE");

N  = 60;      // points per surface
tt = 0.12;    // NACA 00tt thickness, as a fraction of chord

// y_t(x) with the last coefficient closing the trailing edge to a point
// (-0.1036 rather than the -0.1015 that leaves it open), because an open one
// would be two more corners and a blunt base rather than the -1/2 corner this
// model is for.
Macro NacaThickness
  yt = 5*tt*( 0.2969*Sqrt(xc) - 0.1260*xc - 0.3516*xc*xc
            + 0.2843*xc*xc*xc - 0.1036*xc*xc*xc*xc );
Return

// Upper surface, leading edge to trailing edge, cosine spaced so the points
// bunch where the curvature is.
For i In {0:N}
  xc = 0.5*(1 - Cos(Pi*i/N));
  Call NacaThickness;
  Point(1000 + i) = {xc, yt, 0, 1.0};
EndFor

// Lower surface, back from the trailing edge but stopping short of both ends:
// the two tips are shared with the upper surface.
For i In {1:N-1}
  xc = 0.5*(1 - Cos(Pi*(N-i)/N));
  Call NacaThickness;
  Point(2000 + i) = {xc, -yt, 0, 1.0};
EndFor

Spline(1) = {1000:1000+N, 2001:2000+N-1, 1000};
Curve Loop(1) = {1};
Plane Surface(1) = {1};

// The tunnel: one chord ahead, 1.7 behind for the wake, 0.9 either side.
// Closer in than a solver would want, but the mesher here sizes elements
// uniformly, so every chord of empty tunnel is paid for in elements across the
// section -- and it is the section this model exists to test.
Rectangle(2) = {-1.0, -0.9, 0, 3.6, 1.8};
BooleanDifference(3) = { Surface{2}; Delete; }{ Surface{1}; Delete; };
