// Gmsh project created on Sun Oct 24 17:56:08 2010

Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {0, 1, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};
Line Loop(4) = {2, 3, 1};
Plane Surface(5) = {4};
