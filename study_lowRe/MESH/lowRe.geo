// gmsh lowRe.geo -3 -o lowRe.msh
// gmsh lowRe.geo

// ----------------------
// Parameters
// ----------------------

R = 0.5; // Radius (Diameter = 1.0)
L_far = 30.0;
L_wake = 50.0;
Ev = 1.0;

lc_wall = 0.02; // Fine mesh near wall
lc_far = 1.0; // Coarse mesh at farfield

// ---------------------
// Point
// ---------------------

Point(10) = {0, 0, 0, lc_wall};
Point(1) = {R, 0, 0, lc_wall};
Point(2) = {0, R, 0, lc_wall};
Point(3) = {-R, 0, 0, lc_wall};
Point(4) = {0, -R, 0, lc_wall};

Point(5) = {-L_far, L_far, 0, lc_far};
Point(6) = {-L_far, -L_far, 0, lc_far};
Point(7) = {L_wake, -L_far, 0, lc_far};
Point(8) = {L_wake, L_far, 0, lc_far};

// ---------------------
// Line
// ---------------------

Circle(1) = {1, 10, 2};
Circle(2) = {2, 10, 3};
Circle(3) = {3, 10, 4};
Circle(4) = {4, 10, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(1) = {2, 1};

// ---------------------
// Volume
// ---------------------

// Make the entire surface a structured grid
// Transfinite Surface {1};
Recombine Surface {1};

// Extrude with thickness Ev
layers[] = Extrude {0, 0, Ev} {
    Surface{1}; Layers{1}; Recombine;
};

// ----------------------
// Physical Groups 
// ----------------------
Physical Surface("canal_top") = {33};
Physical Surface("canal_sidewalls") = {1, 50};
Physical Surface("canal_inlet") = {21};
Physical Surface("canal_outlet") = {29};
Physical Surface("canal_bottom") = {25};

Physical Surface("cylinder_walls") = {37, 41, 45, 49};

Physical Volume("Fluid") = {1};
