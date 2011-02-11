#include <Fepic/Mesh>


int const DIM = 2;

using namespace std;
using namespace Eigen;


typedef DefaultTraits<DIM, Simplex<DIM> >  MyT;
typedef iMesh<MyT>                         MyMesh;
typedef MyMesh::CellT                      MyCell;
typedef MyMesh::PointT                     MyPoint;

int main(int argc, char *argv[]) {

  //stringstream ss;
  MyMesh malha;
  //MyCell c({0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30});
  //MyPoint p;
  //MyEdge e;

  //if (argc>1)
  //{
    //malha.readFileMsh(argv[1]);
  //}
  //else
  //{
    //std::cout << "file path(malha/simple.msh) : ";
    //std::string name;
    //cin >> name;
    //malha.readFileMsh(name.data()); 
  //}

  if(DIM==3)
  {
    //malha.readFileMsh("malha/monotet.msh");
    malha.readFileMsh("malha/tetbig.msh");
  }
  else
  {
    malha.readFileMsh("malha/simple.msh");
  }
  
  //malha.printInfo(std::cout);
  
  //malha.setOrder(2);
  
  malha.writeVtk();
  malha.writeFsf();
  
  //malha.writeFileState();
  
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
