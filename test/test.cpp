//#include <Fepic/DofHandler>
#include <Fepic/Mesh>
//#include <tr1/memory>

using namespace std;
//using namespace Eigen;

void printMesh(Mesh *mesh);


//int main(int argc, char *argv[])
int main()
{

  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  SMesh<Triangle3,3> *mesh;
  int dim = 3;

  mesh = (SMesh<Triangle3,3>*)Mesh::create(TRIANGLE3,dim);
  
  mesh->qBuildAdjacency(true);
  
  msh_reader.readFileMsh("tests/meshes/singular_tri3a.msh", mesh);
  //vtk_printer.attachMesh(mesh);
  //vtk_printer.writeVtk("meshes/out/simpedge2.vtk");

  
  delete mesh;
  
  
  
  
  
  
  
  
  
  
  //VarDofs vvv;
  
  //const char* mesh_name = "malha/tri1.msh";
  //Hexahedron8   hex1;
  //Hexahedron20  hex2;
  //Hexahedron27  hex3;
  //Tetrahedron4  tet1;
  //Tetrahedron10 tet2;
  //Triangle3     tri1;
  //Triangle6     tri2;
  //Quadrangle4   qua1;
  //Quadrangle8   qua2;
  //Quadrangle9   qua3;
  //Facet         f[DIM];

  

  ////std::tr1::shared_ptr<Mesh> mesh(new SMesh<Triangle3>);
  //std::tr1::shared_ptr<Mesh> mesh(Mesh::create(TRIANGLE3));
  ////Mesh * mesh = Mesh::create(TRIANGLE3,-1);
  
  //MeshIoMsh Writer;
  //MeshIoVtk vtk_printer(mesh.get());

  //Writer.readFileMsh(mesh_name, mesh.get());
  
  //mesh->printInfo();
  
  //vtk_printer.setOutputFileName(mesh_name);
  //vtk_printer.writeVtk();
  
  //printMesh(mesh.get());
  

  return 0;
}
//////sed -i '/__/d' log



//void printMesh(Mesh *mesh)
//{
  //printf("NODES \n");
  //for (int i = 0; i < mesh->numNodes(); ++i)
  //{
    //double x[2];
    //mesh->getNode(i)->getCoord(x);
    //printf("%lf %lf\n", *x, x[1]);
  //}
  //printf("\n");
  
  //printf("CELLS\n");
  //for (int i = 0; i < mesh->numCells(); ++i)
  //{
    //int nds[3];
    //mesh->getCell(i)->getNodesId(nds);
    //printf("%d %d %d\n", *nds, nds[1], nds[2]);
  //}
  //printf("\n");
  
//}

