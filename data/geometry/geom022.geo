// Dumbbell: two round lobes joined by a neck a tenth of their width.
//
// The neck is a long channel whose field is essentially parallel everywhere, so
// separatrices entering it run its whole length side by side before they can
// diverge. That is the configuration that produces the very long, very thin
// chords -- the ones whose collapse would drag a curve across the model if the
// width of the strip were not bounded.
SetFactory("OpenCASCADE");
Disk(1) = {-1.0, 0, 0, 0.65, 0.65};
Disk(2) = { 1.0, 0, 0, 0.65, 0.65};
Rectangle(3) = {-1.0, -0.13, 0, 2.0, 0.26};
BooleanUnion(4) = { Surface{1}; Delete; }{ Surface{2}; Surface{3}; Delete; };
