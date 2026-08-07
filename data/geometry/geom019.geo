// Tee: two reflex corners, and no symmetry to help.
//
// The cross field has to reconcile a wide flange with a narrow stem. The two
// 270-degree corners sit close enough that their rays interact, which is where
// separatrices end up crossing each other rather than reaching the boundary.
SetFactory("OpenCASCADE");
Rectangle(1) = {0.0, 1.4, 0, 2.0, 0.6};
Rectangle(2) = {0.7, 0.0, 0, 0.6, 1.5};
BooleanUnion(3) = { Surface{1}; Delete; }{ Surface{2}; Delete; };
