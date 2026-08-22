/*----------  Parameters  ----------*/

H = 0.04;
lc_cyl = 0.0010;
lc_far = 0.0100;

Lup   =  20*H; // front
Lwake =  50*H; // back
Lside =   7*H; // side
Lz    = 0.1*H; // Z

N_cyl   = (H/lc_cyl) +1;
N_limit = ((2*Lside)/lc_far) +1;
N_side  = ((Lup+Lwake)/lc_far) +1;
N_z     = 1; 

/*----------  Point  ----------*/

Point(1) = { H/2,  H/2, 0, lc_cyl};
Point(2) = { H/2, -H/2, 0, lc_cyl};
Point(3) = {-H/2, -H/2, 0, lc_cyl};
Point(4) = {-H/2,  H/2, 0, lc_cyl};

Point(5) = { -Lup,  Lside, 0, lc_far};
Point(6) = { -Lup, -Lside, 0, lc_far};
Point(7) = {Lwake, -Lside, 0, lc_far};
Point(8) = {Lwake,  Lside, 0, lc_far};

/*----------  Lines  ----------*/

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(1) = {2, 1}; // Cut out by writing in the order (outside, inside)

Transfinite Curve {1, 2, 3, 4} = N_cyl Using Progression 1;

Transfinite Curve {5, 7} = N_limit Using Progression 1;
Transfinite Curve {6, 8} = N_side  Using Progression 1;

/*----------  Volume  ----------*/

Recombine Surface {1};

layers[] = Extrude {0, 0, Lz} { 
    Surface{1}; Layers{N_z}; Recombine;
};

Physical Surface("canal_top")       = {33};
Physical Surface("canal_bottom")    = {25};
Physical Surface("canal_inlet")     = {21};
Physical Surface("canal_outlet")    = {29};
Physical Surface("canal_sidewalls") = {1, 50};

Physical Surface("cyl_top")         = {37};
Physical Surface("cyl_bottom")      = {45};
Physical Surface("cyl_limits_xmin") = {41};
Physical Surface("cyl_limits_xmax") = {49};

Physical Volume("vol") = {1};
