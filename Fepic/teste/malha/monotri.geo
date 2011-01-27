// Gmsh project created on Sun Jan 16 14:44:28 2011
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {0, 1, 0, 1.0};
Line(1) = {2, 1};
Line(2) = {1, 3};
Line(3) = {3, 2};
Line Loop(4) = {3, 1, 2};
Plane Surface(5) = {4};
