// Gmsh project created on Mon Oct 25 15:56:52 2010
Point(1) = {0, 0, 0, 1.0};
Point(2) = {0.5, 0, 0, 1.0};
Point(3) = {0, 0.5, 0, 1.0};
Point(4) = {-0.5, 0, 0, 1.0};
Point(5) = {0, -0.5, 0, 1.0};
Circle(1) = {5, 1, 2};
Circle(2) = {2, 1, 3};
Circle(3) = {3, 1, 4};
Circle(4) = {4, 1, 5};
Line Loop(5) = {2, 3, 4, 1};
Plane Surface(6) = {5};
Physical Line(1) = {2, 3, 4, 1};
Physical Surface(2) = {6};