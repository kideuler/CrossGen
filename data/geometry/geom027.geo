// Staircase: five treads, so five reflex corners and five convex ones
// alternating along a single run of boundary.
//
// Every corner wants to place a singularity nearby and they are close enough
// together to compete. It is the case where the boundary index assignment and
// the field's interior singularities have to agree over a long stretch, and
// where they disagree the partition shows it as components that are not
// four-sided.
SetFactory("OpenCASCADE");
Point(1) = {0.0, 0.0, 0, 1.0};
Point(2) = {2.5, 0.0, 0, 1.0};
Point(3) = {2.5, 0.5, 0, 1.0};
Point(4) = {2.0, 0.5, 0, 1.0};
Point(5) = {2.0, 1.0, 0, 1.0};
Point(6) = {1.5, 1.0, 0, 1.0};
Point(7) = {1.5, 1.5, 0, 1.0};
Point(8) = {1.0, 1.5, 0, 1.0};
Point(9) = {1.0, 2.0, 0, 1.0};
Point(10) = {0.5, 2.0, 0, 1.0};
Point(11) = {0.5, 2.5, 0, 1.0};
Point(12) = {0.0, 2.5, 0, 1.0};
Line(1) = {1, 2};    Line(2) = {2, 3};    Line(3) = {3, 4};
Line(4) = {4, 5};    Line(5) = {5, 6};    Line(6) = {6, 7};
Line(7) = {7, 8};    Line(8) = {8, 9};    Line(9) = {9, 10};
Line(10) = {10, 11}; Line(11) = {11, 12}; Line(12) = {12, 1};
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
Plane Surface(1) = {1};
