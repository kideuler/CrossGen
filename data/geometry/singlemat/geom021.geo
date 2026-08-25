// Annulus with a blind notch: the limit cycle of Fig. 16.
//
// The notch stops short of the hole, so the model is still an annulus and a
// streamline can still circle the hole forever -- Poincare-Bendixson case 3 of
// Sec. 3.3. That is the point: cutting all the way through would leave a disc
// with four corners, which is one quad and no test of anything. Here the notch
// gives the separatrices somewhere to start while the ring of material behind
// it still admits a closed streamline, so the second stopping condition (cut a
// separatrix that meets the same one twice) has to be what ends the trace.
//
// Fig. 16 is the case the paper cannot finish: the T-junction the limit cycle
// leaves cannot be removed without inserting a 3- and 5-valent pair.
SetFactory("OpenCASCADE");
Disk(1) = {0, 0, 0, 1.0, 1.0};
Disk(2) = {0, 0, 0, 0.45, 0.45};
BooleanDifference(3) = { Surface{1}; Delete; }{ Surface{2}; Delete; };
Rectangle(4) = {0.62, -0.09, 0, 0.6, 0.18};
BooleanDifference(5) = { Surface{3}; Delete; }{ Surface{4}; Delete; };
