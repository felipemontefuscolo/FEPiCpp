#include "prefix_header.hpp"
//#include <iostream>
//#include "Fepic/mesh/mesh.hpp"
//#include "Fepic/misc/misc.hpp"
//#include "Fepic/shapefunctions/parametric_pts.hpp"
//#include "Fepic/shapefunctions/shapefunctions.hpp"
//#include "Fepic/quadrature/quadrature.hpp"
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
typedef iMesh<MyT> 					      MyMesh;
typedef MyMesh::CellT 					  MyCell;
typedef MyMesh::PointT 					  MyPoint;
typedef MyMesh::EdgeT 					  MyEdge;

int main(int argc, char *argv[]) {

	//stringstream ss;
	iMesh<MyT> 	malha;
    //MyCell 		c({0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30});
    //MyPoint 	p;
    //MyEdge 		e;
    
	if(DIM == 3)
	{
		//malha.readFileMSH("malha/monotet.msh");
        malha.readFileMSH("malha/tetbig.msh");
	}
	else
	{
		malha.readFileMSH("malha/simple.msh");
	}
	 
    //malha.printInfo(std::cout);
    
    //malha.setOrder(2);
    
    //malha.writeVTK();


    malha.writeFileState();
    
    //auto v = MyCell::getOppFLN(malha.getOrder());
    
    //auto A = malha.getCell(0)->getBorderNodes(0, malha);
    //auto B = malha.getCell(1)->getBorderNodes(1, malha);
    
    //auto temp(B);
    //for (uint i = 0; i < B.size(); ++i)
    //{
        //B[i] = temp[v[0][i]];
    //}
    
    
    //printArray(A, cout, A.size()); cout << endl;
    //printArray(B, cout, B.size()); cout << endl;
    
    return 0;
}
//////sed -i '/__/d' log
