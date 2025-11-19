
// Schäfer & Turek (1996) 3D circular cylinder case (3D-Z)
// gmsh cylinder3d_circle.geo -3 -o cylinder3d.msh
// gmsh cylinder3d.msh &

// pillar
Dp = 0.1;
e = 0.1;
r = e/Dp;

// domain
H = 0.41; // channel height
L = 2.5; // channel length
p = 0.2; // pillar depth
Ev = 0.41;

xc = 0.45 + Dp/2;
yc = 0.15 + Dp/2;
R  = Dp/2;

// Target mesh size
// lc = 0.1;
// N_bord = 50;
// N_bas = (N_bord - 1)*Dp/e + 1;
// N_walls = 100;
// N_inlet = (1.0*N_walls - 1)*L/H + 1;
// N_outlet = (N_walls - 1)*L/H + 1; 

lc = 0.005;
N_bord = 0;
N_bas = 0;
N_walls = 0;
N_inlet = 0;
N_outlet = 0; 

// ---------------------
// Point
// ---------------------

Point(1) = {xc+R, yc, 0, lc};
Point(2) = {xc, yc+R, 0, lc};
Point(3) = {xc-R, yc, 0, lc};
Point(4) = {xc, yc-R, 0, lc};
Point(10) = {xc, yc, 0, lc};

Point(5) = {0, 0, 0, lc};
Point(6) = {L, 0, 0, lc};
Point(7) = {L, H, 0, lc};
Point(8) = {0, H, 0, lc};

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

// Line 5 & 7 devided by N_walls
// "Using Progression 1" means equally spaced

// Transfinite Curve {5, 7} = N_walls Using Progression 1;
// Transfinite Curve {8} = N_inlet  Using Progression 1;
// Transfinite Curve {6} = N_outlet Using Progression 1;
// Transfinite Curve {1, 3} = N_bord Using Progression 1;
// Transfinite Curve {2, 4} = N_bas  Using Progression 1;

// ---------------------
// Volume
// ---------------------

// Make the entire surface a structured grid
// Transfinite Surface {1};
Recombine Surface {1};

// Extrude with thickness Ev
layers[] = Extrude {0, 0, Ev} {
    Surface{1}; Layers{20}; Recombine;
};

Physical Surface("inlet") = {33};
Physical Surface("outlet") = {25};
Physical Surface("top") = {29};
Physical Surface("bottom") = {21};
Physical Surface("walls") = {1, 50};
Physical Surface("pillar_walls") = {37, 41, 45, 49};

Physical Volume("vol") = {1};
