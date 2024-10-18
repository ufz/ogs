lc = 3.0;
nh = 8;
//+
Point(1) = {0, 0, 0, lc};
//+
Point(2) = {18, 0, 0, lc};
//+
Point(3) = {18, 2, 0, lc};
//+
Point(4) = {15, 2, 0, lc};
//+
Point(5) = {6, 6.5, 0, lc};
//+
Point(6) = {0, 6.5, 0, lc};
//+
Point(7) = {15, 0, 0, lc};
//+
Point(8) = {6, 0, 0, lc};

//+
Line(1) = {1, 8};
//+
Line(2) = {8, 7};
//+
Line(3) = {7, 2};
//+
Line(4) = {2, 3};
//+
Line(5) = {3, 4};
//+
Line(6) = {4, 5};
//+
Line(7) = {5, 6};
//+
Line(8) = {6, 1};
//+
Line(9) = {5, 8};
//+
Line(10) = {4, 7};
//+
Curve Loop(1) = {8, 1, -9, 7};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, 2, -10, 6};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {10, 3, 4, 5};
//+
Plane Surface(3) = {3};

Extrude {0, 0, 20} {
	Surface{1,2,3};
}

Transfinite Line "*" = nh;
Transfinite Line {17,18,22,26,44,48,66,70} = 3*nh;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
//+
//+
Physical Surface("cote", 77) = {32, 54, 76, 1, 2, 3};
//+
Physical Surface("bas", 78) = {23, 45, 67};
//+
Physical Surface("derriere", 79) = {19};
//+
Physical Surface("devant", 80) = {71};
//+
Physical Volume("bloc", 81) = {2, 1, 3};
