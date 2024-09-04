// GMSH input file created by OpenGeoSys 6.5.2-325-g593bae04

Point(0) = {0, 0, 0, 0.20000000000000001};
Point(1) = {0, 1, 0, 0.20000000000000001};
Point(2) = {1, 0, 0, 0.20000000000000001};
Point(3) = {1, 1, 0, 0.20000000000000001};
Line(0) = {0,1};
Line(1) = {1,3};
Line(2) = {3,2};
Line(3) = {2,0};
Line Loop(4) = {0,1,2,3};
Plane Surface(0) = {4};
