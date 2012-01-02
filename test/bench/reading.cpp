#include <iostream>
#include <Fepic/Mesh>
#include <tr1/memory>
#include <omp.h>

using namespace std;
//using namespace Eigen;


void readTet(std::tr1::shared_ptr<Mesh> & mesh)
{
  
  mesh = std::tr1::shared_ptr<Mesh>(Mesh::create(TETRAHEDRON4));
  
  const char* mesh_name = "meshes/big_tet.msh";
  
  MeshIoMsh reader;
  
  reader.readFileMsh(mesh_name, mesh.get());
  
  
}

void f(void *) {};

void test(std::tr1::shared_ptr<Mesh> & mesh)
{
  int iCs0[FEPIC_MAX_ICELLS];
  int iCs1[FEPIC_MAX_ICELLS];
  int viCs0[FEPIC_MAX_ICELLS];
  int viCs1[FEPIC_MAX_ICELLS];

  int n_nodes = mesh->numNodes();
  mesh->timer.restart();
  for (int i = 0; i < n_nodes; ++i)
  {
    Point * p = mesh->getNode(i);
    int C = p->getIncidCell();
    int vC= p->getPosition();
    mesh->vertexStar(C, vC, iCs0, viCs0);
  }
  mesh->timer.elapsed("vertexStar() all mesh");
  
  mesh->timer.restart();
  for (int i = 0; i < n_nodes; ++i)
  {
    Point * p = mesh->getNode(i);
    int C = p->getIncidCell();
    int vC= p->getPosition();
    mesh->vertexStarWithoutColor(C, vC, iCs1, viCs1);
  }  
  mesh->timer.elapsed("vertexStarWithoutColor() all mesh");

  
}


//int main(int argc, char *argv[])
int main()
{
  std::tr1::shared_ptr<Mesh> mesh;
  
  readTet(mesh);

  test(mesh);
  
  mesh->printInfo();
  mesh->printStatistics();
  printf("===============================\n");
  mesh->timer.printTimes();
  mesh->timer.printMethod();    
  
  return 0;
}

