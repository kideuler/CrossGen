// Perforated plate: a three-by-two grid of holes.
//
// Each hole forces index into the field around it, so this is the many-
// singularity end of the range with a regular structure the layout ought to
// reproduce. If the partition comes out irregular on a geometry this regular,
// the irregularity is the algorithm's and not the model's.
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 3.0, 2.0};
Disk(2) = {0.6, 0.6, 0, 0.25, 0.25};
Disk(3) = {1.5, 0.6, 0, 0.25, 0.25};
Disk(4) = {2.4, 0.6, 0, 0.25, 0.25};
Disk(5) = {0.6, 1.4, 0, 0.25, 0.25};
Disk(6) = {1.5, 1.4, 0, 0.25, 0.25};
Disk(7) = {2.4, 1.4, 0, 0.25, 0.25};
BooleanDifference(8) = { Surface{1}; Delete; }
                      { Surface{2}; Surface{3}; Surface{4};
                        Surface{5}; Surface{6}; Surface{7}; Delete; };
