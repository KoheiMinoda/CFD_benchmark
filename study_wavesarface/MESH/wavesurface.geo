// wavesurface.geo
// gmsh wavesurface.geo -3 -o wavesurface.msh
// gmsh wavesurface.msh &

// ---------- Parameters ----------
d = 1.0; // water depth [m]
H_air = 0.5; // height of air region above SWL [m]
Lw = 15.0; // length from wavemaker to start of beach [m]
Lb = 5.0; // length of beach (bottom slope) [m]
Ltot = Lw + Lb;// total tank length [m]

// 2D Mesh resolution
Nx_w = 150; // number of cells along x in flat-bottom region
Nx_b = 80; // number of cells along x along beach
Ny_tot = 60; // number of cells in vertical direction (water + air)

// 3D thickness (quasi-2D)
Ev = 0.05; // thickness in z [m]
Nz = 1; // number of layers in z

// (Not essential when using Transfinite, but kept for safety)
lc = 0.5 * d; // characteristic length scale (fallback)

// ---------- 2D Points ----------
Point(1) = {0.0, -d, 0.0, lc}; // bottom left
Point(2) = {0.0, H_air, 0.0, lc}; // top left
Point(3) = {Lw, -d, 0.0, lc}; // junction flat bottom / beach (bottom)
Point(4) = {Lw, H_air, 0.0, lc}; // junction flat/top over beach
Point(5) = {Ltot, 0.0, 0.0, lc}; // end of beach (bottom meets SWL)
Point(6) = {Ltot, H_air, 0.0, lc}; // top right

// ---------- 2D Lines ----------
Line(1) = {1, 2}; // left wall
Line(2) = {2, 4}; // top from x=0 to x=Lw
Line(3) = {4, 6}; // top from x=Lw to x=Ltot
Line(4) = {6, 5}; // right wall downward
Line(5) = {5, 3}; // sloping beach bottom
Line(6) = {3, 1}; // flat bottom

Line(7) = {4, 3}; // internal interface

// ---------- 2D Surfaces ----------
// Surface 1: rectangle
Line Loop(10) = {1, 2, 7, 6};
Plane Surface(11) = {10};

// Surface 2: quadrilateral with sloping bottom
Line Loop(12) = {-7, 3, 4, 5};
Plane Surface(13) = {12};

// ---------- Structured 2D mesh ----------
Transfinite Curve{1, 4, 7} = Ny_tot Using Progression 1; // vertical lines
Transfinite Curve{2, 6} = Nx_w Using Progression 1;  // flat-bottom region
Transfinite Curve{3, 5} = Nx_b Using Progression 1;  // beach region

Transfinite Surface{11};
Transfinite Surface{13};
Recombine Surface{11, 13}; // quadrilateral elements in 2D

// --------------------------------------------------------------------
// 3D extrusion in z (thin layer)
// --------------------------------------------------------------------
layers1[] = Extrude {0, 0, Ev} {
    Surface{11}; Layers{Nz}; Recombine;
};

layers2[] = Extrude {0, 0, Ev} {
    Surface{13}; Layers{Nz}; Recombine;
};

// layers1[0], layers2[0] : "back" surfaces (z = Ev)
// layers1[1], layers2[1] : volumes created by extrusion

// ---------------- Physical groups for Code_Saturne ----------------
Physical Volume("fluid") = {layers1[1], layers2[1]};

Physical Surface("front") = {11, 13};
Physical Surface("back") = {35, 57};

Physical Surface("top") = {26, 48};
Physical Surface("bottom") = {34};
Physical Surface("beach") = {56};

Physical Surface("wavemaker") = {22};
Physical Surface("outlet") = {52};
