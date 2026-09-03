// Cross: four reflex corners, arranged symmetrically.
//
// Every one of the four 270-degree corners launches two rays, and the shape's
// symmetry means the field's singularities come out symmetrically placed too --
// which is the case Sec. 4 is written for. On a discrete mesh they never line
// up exactly, so the separatrices that ought to coincide run alongside each
// other instead and leave the thin strips the chord collapse exists to remove.
SetFactory("OpenCASCADE");
Rectangle(1) = {0.7, 0.0, 0, 0.6, 2.0};
Rectangle(2) = {0.0, 0.7, 0, 2.0, 0.6};
BooleanUnion(3) = { Surface{1}; Delete; }{ Surface{2}; Delete; };
