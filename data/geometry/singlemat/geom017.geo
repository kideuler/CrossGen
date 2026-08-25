// L-shape: the canonical reflex corner.
//
// One 270-degree corner, which Table 1 gives index -1/4 and which therefore
// launches two rays into the domain instead of none. The rest of the boundary
// is straight, so whatever the field does is the corner's doing.
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 2, 2};
Rectangle(2) = {1, 1, 0, 1.1, 1.1};
BooleanDifference(3) = { Surface{1}; Delete; }{ Surface{2}; Delete; };
