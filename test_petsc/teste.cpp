#include <iostream>
#include "../mesh/mesh.hpp"
#include "../misc/misc.hpp"
#include "../misc/ivec.hpp"
#include "../mypetsc/mypetsc.hpp"
#include <cmath>
#include "Fepic/src/custom_eigen/custom_eigen.hpp"
#include <vector>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <limits>
#include <bitset>
#include "petscksp.h"

using namespace std;
//using namespace Eigen;

const int DIM=2;
typedef unsigned unsigned;
typedef Default_Traits<DIM, Simplex<DIM> > MyT;
typedef iMesh<MyT> MyMesh;

int main(int argc, char *argv[]) {

	PetscInitialize(&argc,&argv,0,0);
	
	//PetscErrorCode  ierr;
	//KSP    			ksp;
	//PC	   			pc;
	
	MyMesh malha;
	
	malha.readFileMsh("mesh/circleA.msh");
	
	malha.writeVtk();
	
	vector<double> value(malha.numValidNodes());
	vector<double>::iterator vit = value.begin();
	
	unsigned counter=0;
	for (; vit != value.end() ; ++vit)
	{
		*vit = counter++;
	}
	
	malha.addNodeScalarVtk("value", value);
	counter=0;
	vit = value.begin();
	for (; vit != value.end() ; ++vit)
	{
		*vit = malha.getNode(counter++)->getCoord(0);
	}	
	malha.addNodeScalarVtk("other", value);
	
    PetscFinalize();
    
	return 0.;
}
