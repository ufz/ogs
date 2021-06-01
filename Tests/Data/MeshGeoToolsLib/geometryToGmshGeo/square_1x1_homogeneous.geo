// GMSH input file created by OpenGeoSys 6.4.0-502-gddafba0cf.dirty

Point(0) = {0, 0, 0, 0.01};
Point(1) = {0, 1, 0, 0.01};
Point(2) = {1, 0, 0, 0.01};
Point(3) = {1, 1, 0, 0.01};
Line(0) = {0,1};
Line(1) = {1,3};
Line(2) = {3,2};
Line(3) = {2,0};
Line Loop(4) = {0,1,2,3};
Plane Surface(0) = {4};
