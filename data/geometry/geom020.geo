// Slotted plate: three deep, narrow slots cut in from one edge.
//
// Six reflex corners in a row, and between each pair of slots a tongue only a
// few elements wide. Thin features are where separatrix tracing produces
// components a fraction of an element across, so this is a direct test of
// whether the layout stays four-sided when the geometry itself is slender.
SetFactory("OpenCASCADE");
Rectangle(1) = {0.0, 0.0, 0, 3.0, 1.6};
Rectangle(2) = {0.5, 0.6, 0, 0.22, 1.1};
Rectangle(3) = {1.4, 0.6, 0, 0.22, 1.1};
Rectangle(4) = {2.3, 0.6, 0, 0.22, 1.1};
BooleanDifference(5) = { Surface{1}; Delete; }{ Surface{2}; Surface{3}; Surface{4}; Delete; };
