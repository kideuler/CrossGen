// Crescent: a curved channel with only two corners on the whole boundary.
//
// This one is here because it fails, and it is the smallest thing that fails
// that way. The two tips measure 39 degrees, so Table 1 gives each index 1/4
// and the boundary supplies 1/2 in total; the model is a disc, so the field
// owes another 1/2 from two interior singularities. It produces none, and the
// partition is a single component with two sides.
//
// So either the diffusion-generated field is not meeting the index constraint
// here, or the singularities are there and are not being seen -- note that
// CrossField::computeSingularities() skips every triangle touching a boundary
// vertex, which on a shape this curved is most of the ones near the tips.
// Nothing else in data/geometry separates those two possibilities, because
// every other model has enough corners to satisfy the constraint without the
// interior having to.
SetFactory("OpenCASCADE");
Disk(1) = {0.00, 0, 0, 1.00, 1.00};
Disk(2) = {0.62, 0, 0, 0.78, 0.78};
BooleanDifference(3) = { Surface{1}; Delete; }{ Surface{2}; Delete; };
