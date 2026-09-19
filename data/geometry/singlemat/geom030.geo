// Three NACA 0012 sections stacked one above another: a linear cascade, the
// arrangement used for a turbine or compressor blade row.
//
// What geom029 has one of, this has three of, and the interesting part is what
// sits between them. Each section contributes a -1/2 trailing edge, and the
// gaps are channels about half a chord across bounded by curved walls on both
// sides -- so a streamline entering one runs the length of the blade with the
// field nearly parallel to it, which is where separatrices travel side by side
// for a long way before anything separates them.
//
// It is also the most index this set asks a field to carry. The domain has
// three holes, so the three trailing edges at -1/2, the three leading edges the
// mesh reads as corners, and the four square corners of the tunnel all have to
// reconcile against a Euler characteristic of -2, and whatever is left over the
// interior singularities have to supply.
//
// The index works out exactly, which is worth writing down because it is the
// one model here where every term is non-zero: four tunnel corners at +1/4,
// three trailing edges at -1/2, three leading edges the mesh reads at -1/4,
// which is -5/4 on the boundary; the domain has three holes so its Euler
// characteristic is -2; and the interior therefore owes -3/4. The field
// produces exactly three -1/4 singularities, and they sit at x = -0.54, one
// ahead of each leading edge, where the stagnation structure would be.
//
// Mesh size: same section-to-tunnel ratio as geom029, so np=150 puts about five
// elements across each section. Measured:
//
//     np    triangles   raw layout                simplified
//     75        8780     53 components, all 4-sided   51, one T-junction
//    150       34520     67 components, 4 not         31, all four-sided
//    250       94824     64 components, 4 not         42, all four-sided
//
// Unlike geom029 this one does not settle: the raw partition and the number of
// chords Sec. 4 finds to collapse both move with the mesh. That is the model
// earning its place rather than a fault in it -- the channels between the
// sections are the configuration the greedy collapse order is least stable on,
// and nothing else in the set shows it this clearly. data/meshes carries np=150.
SetFactory("OpenCASCADE");

N  = 60;      // points per surface
tt = 0.12;    // NACA 00tt thickness, as a fraction of chord
pitch = 0.6;  // spacing between sections, in chords

// y_t(x), closed to a point at the trailing edge (-0.1036 rather than -0.1015)
// so that each section ends in the -1/2 corner this model is for.
Macro NacaThickness
  yt = 5*tt*( 0.2969*Sqrt(xc) - 0.1260*xc - 0.3516*xc*xc
            + 0.2843*xc*xc*xc - 0.1036*xc*xc*xc*xc );
Return

For j In {0:2}
  y0 = (j - 1) * pitch;
  base = 1000 + 1000*j;

  // Upper surface, leading edge to trailing edge, cosine spaced.
  For i In {0:N}
    xc = 0.5*(1 - Cos(Pi*i/N));
    Call NacaThickness;
    Point(base + i) = {xc, y0 + yt, 0, 1.0};
  EndFor

  // Lower surface back again, stopping short of both tips: they are shared.
  For i In {1:N-1}
    xc = 0.5*(1 - Cos(Pi*(N-i)/N));
    Call NacaThickness;
    Point(base + N + i) = {xc, y0 - yt, 0, 1.0};
  EndFor

  Spline(100 + j) = {base : base+N, base+N+1 : base+N+N-1, base};
  Curve Loop(100 + j) = {100 + j};
  Plane Surface(100 + j) = {100 + j};
EndFor

// The tunnel, sized as geom029's so the sections resolve the same way.
Rectangle(1) = {-0.9, -1.2, 0, 3.6, 2.4};
BooleanDifference(200) = { Surface{1}; Delete; }
                        { Surface{100}; Surface{101}; Surface{102}; Delete; };
