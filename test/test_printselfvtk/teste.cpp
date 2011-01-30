#include "prefix_header.hpp"
//#include <iostream>
//#include "Fepic/src/mesh/mesh.hpp"
//#include "Fepic/src/util/misc.hpp"
//#include "Fepic/src/shapefunctions/parametric_pts.hpp"
//#include "Fepic/src/shapefunctions/shapefunctions.hpp"
//#include "Fepic/src/quadrature/quadrature.hpp"
//#include <cmath>
//#include <Eigen/Dense>
//#include <Eigen/Geometry>
//#include <vector>
//#include <cstring>
//#include <iomanip>
//#include <sstream>
//#include <limits>
//#include <cstdlib>
//#include <sstream>
//#include <algorithm>


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
    MyCell 		c({0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26});
    MyPoint 	p;
    MyEdge 		e;
    
	//if(DIM == 3)
	//{
		//malha.readFileMSH("malha/tettet.msh");
	//}
	//else
	//{
		//malha.readFileMSH("malha/tritri.msh");
	//}
	 
    //malha.printInfo(std::cout);

    MyCell::iCellData a;
    
    int n;
    
    if (argc>1)
        n=atoi(argv[1]);
    else
        n=1;
    
    auto temp = a.getMinimesh(n);
    auto coord= genTetParametricPtsINT(n);
    //auto coord= genTriParametricPtsINT(n);
    
    for (uint i = 0; i < coord.size(); i++)
    {
        malha.addPoint(MyPoint{coord[i]});
    }
    
    for (uint i = 0; i < temp.size(); i++)
    {
        Fepic::vectorui v(temp[i].begin(), temp[i].end());
        malha.addCell(MyCell(v));
    }
    
    malha.writeVTK();
    
    cout << temp.size() << endl;
    
    c.printSelfVTK(cout,n); cout << endl;

	return 0.;
}
//////sed -i '/__/d' log
