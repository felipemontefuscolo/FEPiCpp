#include <Fepic/Mesh>


int const DIM = 3;

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
    //malha.readFileFsf("malha/monotet.fsf");
    malha.readFileMsh("malha/tetbig.msh"); // 6.85 e 1.20
  }
  else
  {
    malha.readFileFsf("malha/simple.fsf");
    malha.readFileMsh("malha/simple.msh");
  }
  
  //malha.printInfo(std::cout);
  
  //malha.setOrder(2);
  
  malha.writeVtk();
  //malha.writeFsf("malha/monotet-2.fsf");
  
  //cout << MyCell::getMinimesh(1) << endl;

  
  return 0;
}
//////sed -i '/__/d' log
