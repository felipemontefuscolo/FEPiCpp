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

Plane Surface(10) = {5};

Physical Line(1) = {1, 4};
Physical Line(2) = {3, 2};
Physical Point(1) = {1, 2};
Physical Point(2) = {3, 4};


Physical Surface(3) = {10};
