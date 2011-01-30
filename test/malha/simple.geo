// Gmsh project created on Sun Jan 16 16:57:36 2011
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};
Line(1) = {2, 1};
Line(2) = {1, 4};
Line(3) = {4, 3};
Line(4) = {3, 2};
Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};
