#include <iostream>
#include "Fepic/src/mesh/mesh.hpp"
#include "Fepic/src/util/misc.hpp"
#include "Fepic/src/misc/ivec.hpp"
#include "Fepic/src/shapefunctions/parametric_pts.hpp"
#include "Fepic/src/shapefunctions/shapefunctions.hpp"
#include "Fepic/src/quadrature/quadrature.hpp"
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <limits>
#include <cstdlib>
#include <sstream>
#include <algorithm>


int const DIM = 3;

using namespace std;
using namespace Eigen;


typedef DefaultTraits<DIM, Simplex<DIM> > MyT;
typedef iMesh<MyT> 						  MyMesh;
typedef MyMesh::CellT 					  MyCell;
typedef MyMesh::PointT 					  MyPoint;
typedef MyMesh::EdgeT 					  MyEdge;

int main(int argc, char *argv[]) {

	stringstream ss;
	iMesh<MyT> 	malha;
    MyCell 		c;
    MyPoint 	p;
    MyEdge 		e;
    
	if(DIM == 3)
	{
		malha.readFileMSH("tet1.msh");
	}
	else
	{
		malha.readFileMSH("malha/tritri.msh");
	}

	return 0.;
}

