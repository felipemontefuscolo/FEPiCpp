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


int const DIM = 2;

using namespace std;
using namespace Eigen;


typedef DefaultTraits<DIM, Simplex<DIM> > MyT;
typedef iMesh<MyT> MyMesh;
typedef MyMesh::CellT MyCell;
typedef MyMesh::PointT MyPoint;
typedef MyMesh::EdgeT MyEdge;

int main(int argc, char *argv[]) {

    stringstream ss;
    iMesh<MyT> malha;
    MyCell c;
    MyPoint p;
    MyEdge e;
	
    Quadrature<Simplex<DIM> > Q;
    
    Q.setOrder(1);
    
	if(DIM == 3)
	{
		malha.readFileMSH("malha/tettet.msh");
	}
	else
	{
		malha.readFileMSH("malha/tritri.msh");
	}
	 
    
	int n = 2;
	if (argc==2) ss << argv[1];
	ss >> n;
	
	MyPoint *ppt;
	
	ShapeFunction<Simplex<DIM> > Phi(n,0);
	//for (int i = 0; i < Phi.getNumDof(); i++)
	//{
		//cout << Phi({1./3,1./3}, i) << endl;		
	//}
	
	vector<VectorXd> func(Phi.getNumDof(), VectorXd(malha.getNumNodes()));

    double sum = 0;
	for (uint i = 0; i < malha.getNumNodes(); i++)
	{
		ppt = malha.getNode(i);
		
        
		for (int k = 0; k < Phi.getNumDof(); k++)
		{
			//func[k][i] = Phi(ppt->getCoord(), k);
			double x = ppt->getCoord(0);
			double y = ppt->getCoord(1);
            func[k][i] = Phi.gradL(x, y, k)[0];
            //func[k][i] = -((384.*x*x-192.*x+16.)*y+512.*x*x*x-672.*x*x+224.*x-16.)/3.;
            sum += func[k][i];
		}
		
	}
    
    //cout << sum/func[0].size() << endl;
	
	//for (uint i = 0; i < malha.getNumNodes(); i++)
	//{
		//ppt = malha.getNode(i);
		
		//for (int k = 0; k < n_nodes; k++)
		//{
			////func[k][i] = Phi(k, ppt->getCoord());
			//func[k][i] = Phi.gradL(k, ppt->getCoord())[2];
		//}
		
	//}
	
	malha.writeVTK();
	for (int k = 0; k < Phi.getNumDof(); k++)
	{	
		stringstream ss;
		ss << "func" << std::setfill('0') << std::setw(5) << k;
		//cout << ss.str().data() << endl;
		malha.addScalarVTK(ss.str().data(), func[k], func[k].size());
		ss.str(std::string());
        
	}	

	//malha.getCell(5)->printSelfVTK(std::cout);

	//int IAJSHDKJAD[]={0,0,0};
	//vector<int> a(3);
	//int BBB[] = {1,1}
	//uint nodes[] = {21,5,13};
	//uint halfID;
	//cout << iMesh<MyT>::theseVerticesFormAMHalf(nodes, halfID) << endl;
	
	//cell = malha.getCell(0);
	
	//malha.writeFileState();
    //malha.writeVTK();
    //malha.addNodeLabelVTK("Node_Labes");
    //malha.addNodeHalfVTK("Node_Half");

	//cout << "sucesso\n";
	return 0.;
}

