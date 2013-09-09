// Copyright 2005, Google Inc.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <gtest/gtest.h>
#include <Fepic/Mesh>
#include <Fepic/DofHandler>
#include <vector>
#include <tr1/tuple>
#include <tr1/memory>
#include <algorithm>
#include <functional>
#include <iostream>


using std::tr1::shared_ptr;
using namespace std;
using namespace Eigen;

//TEST(RemoveCellTest, WithTri3)
//{
//  MeshIoMsh msh_reader;
//  MeshIoVtk vtk_printer;
//  Mesh *mesh = NULL;  
//
//  ECellType cell_t     = TRIANGLE3;
//  const char* mesh_in  = "meshes/sing_tri3.msh";
//  const char* mesh_out = "meshes/outtest/mod_tri3-03.vtk";
//  
//  mesh = Mesh::create(cell_t);
//  msh_reader.readFileMsh(mesh_in, mesh);
//  
//  MeshTools::removeCell(mesh->getCellPtr(0), mesh);
//  MeshTools::removeCell(mesh->getCellPtr(2), mesh);
//  
//  //cell_iterator cell = mesh->cellBegin();
//  //cell_iterator cell_end = mesh->cellEnd();
// 
//  vtk_printer.attachMesh(mesh);
//  vtk_printer.writeVtk(mesh_out);
//  vtk_printer.printPointIcellVtk();
//  vtk_printer.printPointPositionVtk();
//  
//  delete mesh;
//}
//

void getAllDofs(std::vector<int> & dat, DofHandler const& DofH, Mesh const* mesh)
{
  int dofs[20];
  
  dat.clear();
  
  // points
  for (int k = 0; k < mesh->numNodesTotal(); ++k)
  {
    Point const* point = mesh->getNodePtr(k);
    //if (!mesh->isVertex(point) || point->isDisabled())
    if ( point->isDisabled())
      continue;

    for (int i = 0; i < DofH.numVars(); ++i)
    {
      DofH.getVariable(i).getVertexDofs(dofs, point);
      for (int j = 0; j < DofH.getVariable(i).numDofsPerVertex(); ++j)
        dat.push_back(dofs[j]);
    }
  }
  // cells
  for (int k = 0; k < mesh->numCellsTotal(); ++k)
  {
    Cell const* cell = mesh->getCellPtr(k);
    if (cell->isDisabled())
      continue;

    for (int i = 0; i < DofH.numVars(); ++i)
    {
      DofH.getVariable(i).getCellDofs(dofs, cell);
      for (int j = 0; j < DofH.getVariable(i).numDofsPerCell(); ++j)
        dat.push_back(dofs[j]);
    }
  }
  // facets
  for (int k = 0; k < mesh->numFacetsTotal(); ++k)
  {
    Facet const* facet = mesh->getFacetPtr(k);
    if (facet->isDisabled())
      continue;

    for (int i = 0; i < DofH.numVars(); ++i)
    {
      DofH.getVariable(i).getFacetDofs(dofs, facet);
      for (int j = 0; j < DofH.getVariable(i).numDofsPerFacet(); ++j)
        dat.push_back(dofs[j]);
    }
  }
  // facets
  for (int k = 0; k < mesh->numFacetsTotal(); ++k)
  {
    Facet const* facet = mesh->getFacetPtr(k);
    if (facet->isDisabled())
      continue;

    for (int i = 0; i < DofH.numVars(); ++i)
    {
      DofH.getVariable(i).getFacetDofs(dofs, facet);
      for (int j = 0; j < DofH.getVariable(i).numDofsPerFacet(); ++j)
        dat.push_back(dofs[j]);
    }
  }
  // corners
  if (mesh->cellDim() == 3)
    for (int k = 0; k < mesh->numCornersTotal(); ++k)
    {
      Corner const* corner = mesh->getCornerPtr(k);
      if (corner->isDisabled())
        continue;

      for (int i = 0; i < DofH.numVars(); ++i)
      {
        DofH.getVariable(i).getCornerDofs(dofs, corner);
        for (int j = 0; j < DofH.getVariable(i).numDofsPerCorner(); ++j)
          dat.push_back(dofs[j]);
      }
    }  
}



TEST(DofHandlerTest, AssignsDofsTri3)
{
  MeshIoMsh msh_reader;
  //MeshIoVtk vtk_printer;
  Mesh *mesh = NULL;  

  ECellType cell_t     = TRIANGLE3;
  const char* mesh_in  = "meshes/simptri3.msh";  
  //const char* mesh_out = "meshes/outtest/dof_tri3.vtk";


  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);

  
  DofHandler DofH(mesh);
  //                         ndpv,  ndpr,  ndpf,  ndpc
  DofH.addVariable("altura",    1,     0,     0,     0); // 12
  DofH.addVariable("vetor",     2,     0,     0,     0); // 24
  DofH.addVariable("foo",       0,     1,     0,     0); // 0
  DofH.addVariable("bar",       0,     0,     1,     0); // 25
  DofH.addVariable("moo",       0,     0,     0,     1); // 14
  
  DofH.SetUp();
  
  EXPECT_EQ(75, DofH.numDofs());
  
  
  // checking variable 1
  int dofs[] = {-2,-2,-2,-2,-2,-2,-2,-2};
  int * dofs_end = dofs+sizeof(dofs)/sizeof(int);
  int nn;
  //Cell *cell;
  CellElement *l;
  for (int k = 0; k < mesh->numNodesTotal(); ++k)
  {
    l = mesh->getNodePtr(k);
    DofH.getVariable(1).getVertexDofs(dofs, l);
    nn = count_if(dofs, dofs_end, bind2nd(greater_equal<int>(),0) ); // non-negative
    EXPECT_EQ(2, nn) << "node id = " << mesh->getCellPtr(l->getIncidCell())->getNodeId(l->getPosition());;
    fill(dofs, dofs_end, -1);
  }
  for (int k = 0; k < mesh->numNodesTotal(); ++k)
  {
    l = mesh->getNodePtr(k);
    DofH.getVariable(1).getCornerDofs(dofs, l);
    nn = count_if(dofs, dofs_end, bind2nd(greater_equal<int>(),0) ); // non-negative
    EXPECT_EQ(2, nn) << "node id = " << mesh->getCellPtr(l->getIncidCell())->getNodeId(l->getPosition());;
    fill(dofs, dofs_end, -1);
  }
  
  MeshTools::removeCell(mesh->getCellPtr(2), mesh);
  MeshTools::removeCell(mesh->getCellPtr(3), mesh);
  
  DofH.SetUp();
  
  //11x3 + 22 + 12
  EXPECT_EQ(11, DofH.getVariable(0).numPositiveDofs());
  EXPECT_EQ(22, DofH.getVariable(1).numPositiveDofs());
  EXPECT_EQ( 0, DofH.getVariable(2).numPositiveDofs());
  EXPECT_EQ(22, DofH.getVariable(3).numPositiveDofs());
  EXPECT_EQ(12, DofH.getVariable(4).numPositiveDofs());
  
  EXPECT_EQ(67, DofH.numDofs());
  
  
  // ckeck
  std::vector<int> dat;
  
  EXPECT_EQ(5, DofH.numVars());
  
  getAllDofs(dat,DofH,mesh);

  std::sort(dat.begin(), dat.end());

  EXPECT_TRUE(dat[0] >= 0); // grau de liberdade definidos em todo lugar

  dat.erase(std::remove(dat.begin(), dat.end(), -1), dat.end()); // remove -1
  dat.erase(std::unique(dat.begin(), dat.end()), dat.end());     // remove duplicated
  
  EXPECT_EQ(DofH.numDofs(), (int)dat.size());
  
  int counter = 0;
  for (unsigned i = 0; i < dat.size(); ++i)
  {
    EXPECT_EQ(counter, dat[i]);
    ++counter;
  }
  
  
  
  delete mesh;
}

class MyGetDataVtkTri6Test : public DefaultGetDataVtk
{
public:
  MyGetDataVtkTri6Test(int * i = NULL, Mesh* m=NULL, DofHandler* df=NULL) :
                      DefaultGetDataVtk(NULL,i), mesh_ptr(m), dofh_ptr(df){}
  
  int get_data_i(int nodeid) const
  {
    int dof;
    Point const* p = mesh_ptr->getNodePtr(nodeid);
    if (!mesh_ptr->isVertex(p))
      return -1;
    dofh_ptr->getVariable(0).getVertexDofs(&dof, p);
    return data_i[dof];
  }
  
  Mesh *mesh_ptr;
  DofHandler *dofh_ptr;
  
  virtual ~MyGetDataVtkTri6Test() {}
  
  // int  * data_i; from MyGetDataVtk
};

TEST(DofHandlerTest, AssignsDofsTri6)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh = NULL;  

  ECellType cell_t     = TRIANGLE6;
  const char* mesh_in  = "meshes/simptri6.msh";  
  const char* mesh_out = "meshes/outtest/dof_tri6.vtk";

  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);
  
  DofHandler DofH(mesh);
  //                         ndpv,  ndpr,  ndpf,  ndpc
  DofH.addVariable("numero",    1,     0,     0,     0); // 12
  DofH.addVariable("vetor",     2,     0,     0,     0); // 24
  DofH.addVariable("foo",       0,     1,     0,     0); //  0
  DofH.addVariable("bar",       0,     0,     1,     0); // 25
  DofH.addVariable("moo",       0,     0,     0,     1); // 14
  
  DofH.SetUp();
  
  EXPECT_EQ(75, DofH.numDofs());
  
  MeshTools::removeCell(mesh->getCellPtr(2), mesh);
  MeshTools::removeCell(mesh->getCellPtr(3), mesh);
  
  DofH.SetUp();
  
  //11x3 + 22 + 12 + 12
  EXPECT_EQ(11, DofH.getVariable(0).numPositiveDofs());
  EXPECT_EQ(22, DofH.getVariable(1).numPositiveDofs());
  EXPECT_EQ( 0, DofH.getVariable(2).numPositiveDofs());
  EXPECT_EQ(22, DofH.getVariable(3).numPositiveDofs());
  EXPECT_EQ(12, DofH.getVariable(4).numPositiveDofs());

  EXPECT_EQ(67, DofH.numDofs());
  
  
  //// .getVariable(0)
  //int *dat = DofH.data();
  //
  //int counter = 0;
  //for (int i = 0; i < DofH.totalSize(); ++i)
  //{
  //  if ((*dat) >= 0)
  //    EXPECT_EQ(counter++, (*dat));
  //    
  //}
  
  vtk_printer.attachMesh(mesh);
  vtk_printer.writeVtk(mesh_out);
  
  std::vector<int> dados(DofH.numDofs());
  
  point_iterator pt = mesh->pointBegin();
  point_iterator pt_end = mesh->pointEnd();
  
  // for variable "0" ("numero")
  for (; pt != pt_end; ++pt)
  {
    if (!mesh->isVertex(&*pt))
      continue;
      
    int dof(0); // "numero"
    DofH.getVariable(0).getVertexDofs(&dof, &*pt);
    dados[dof] = mesh->getNodeContigId(mesh->getPointId(&*pt));
    //dados[dof] = 0;
  }
  
  MyGetDataVtkTri6Test get_data(dados.data(), mesh, &DofH);
  
  vtk_printer.addNodeIntVtk(DofH.getVariable(0).getName(), get_data);
  
  // checking dofs
  std::vector<int> dat;
  
  EXPECT_EQ(5, DofH.numVars());
  
  getAllDofs(dat,DofH,mesh);

  std::sort(dat.begin(), dat.end());

  EXPECT_TRUE(dat[0] == -1); // -1 nos nohs de controle

  dat.erase(std::remove(dat.begin(), dat.end(), -1), dat.end()); // remove -1
  dat.erase(std::unique(dat.begin(), dat.end()), dat.end());     // remove duplicated
  
  EXPECT_EQ(DofH.numDofs(), (int)dat.size());
  
  int counter = 0;
  for (unsigned i = 0; i < dat.size(); ++i)
  {
    EXPECT_EQ(counter, dat[i]);
    ++counter;
  }
  
  
  delete mesh;
}

TEST(DofHandlerTest, AssignsDofsTet10)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh = NULL;  

  ECellType cell_t     = TETRAHEDRON10;
  const char* mesh_in  = "meshes/simptet10.msh";  
  const char* mesh_out = "meshes/outtest/dof_tet10.vtk";

  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);
  
  DofHandler DofH(mesh);
  //                         ndpv,  ndpr,  ndpf,  ndpc
  DofH.addVariable("numero",    1,     0,     0,     0); //  9
  DofH.addVariable("foo",       0,     1,     0,     0); // 26
  DofH.addVariable("bar",       0,     0,     1,     0); // 30
  DofH.addVariable("moo",       0,     0,     0,     2); // 24
  
  DofH.SetUp();
  
  EXPECT_EQ( 9, DofH.getVariable(0).numPositiveDofs());
  EXPECT_EQ(26, DofH.getVariable(1).numPositiveDofs());
  EXPECT_EQ(30, DofH.getVariable(2).numPositiveDofs());
  EXPECT_EQ(24, DofH.getVariable(3).numPositiveDofs());
  
  EXPECT_EQ(89, DofH.numDofs());
  
  //DofH.SetUp();
  //
  MeshTools::removeCell(mesh->getCellPtr(2), mesh);
  MeshTools::removeCell(mesh->getCellPtr(3), mesh);
  //
  DofH.SetUp();
  //
  
  EXPECT_EQ( 9, DofH.getVariable(0).numPositiveDofs());
  EXPECT_EQ(26, DofH.getVariable(1).numPositiveDofs());
  EXPECT_EQ(28, DofH.getVariable(2).numPositiveDofs());
  EXPECT_EQ(20, DofH.getVariable(3).numPositiveDofs());
  ////
  EXPECT_EQ(83, DofH.numDofs());
  //
  //
  //// .getVariable(0)
  //int *dat = DofH.data();
  //
  //int counter = 0;
  //for (int i = 0; i < DofH.totalSize(); ++i)
  //{
  //  if ((*dat) >= 0)
  //    EXPECT_EQ(counter++, (*dat));
  //    
  //}
  //
  vtk_printer.attachMesh(mesh);
  vtk_printer.writeVtk(mesh_out);
  //
  std::vector<int> dados(DofH.numDofs());
  
  point_iterator pt = mesh->pointBegin();
  point_iterator pt_end = mesh->pointEnd();
  
  // for variable "0" ("numero")
  for (; pt != pt_end; ++pt)
  {
    if (!mesh->isVertex(&*pt))
      continue;
      
    int dof(0); // "numero"
    DofH.getVariable(0).getVertexDofs(&dof, &*pt);
    dados[dof] = mesh->getNodeContigId(mesh->getPointId(&*pt));
    //dados[dof] = 0;
  }
  
  MyGetDataVtkTri6Test get_data(dados.data(), mesh, &DofH);
  
  vtk_printer.addNodeIntVtk(DofH.getVariable(0).getName(), get_data);
  
  // checking dofs
  std::vector<int> dat;
  
  EXPECT_EQ(4, DofH.numVars());
  
  getAllDofs(dat,DofH,mesh);

  std::sort(dat.begin(), dat.end());

  EXPECT_TRUE(dat[0] == -1); // -1 nos nohs de controle

  dat.erase(std::remove(dat.begin(), dat.end(), -1), dat.end()); // remove -1
  dat.erase(std::unique(dat.begin(), dat.end()), dat.end());     // remove duplicated
  
  EXPECT_EQ(DofH.numDofs(), (int)dat.size());
  
  int counter = 0;
  for (unsigned i = 0; i < dat.size(); ++i)
  {
    EXPECT_EQ(counter, dat[i]);
    ++counter;
  }  
  
  
  delete mesh;
}


TEST(DofHandlerTest, BubbleTri3)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh = NULL;  

  ECellType cell_t     = TRIANGLE3;
  const char* mesh_in  = "meshes/simptri3.msh";  
  const char* mesh_out = "meshes/outtest/metis_tri3.vtk";

  ShapeFunction *phi, *psi;
  phi = ShapeFunction::create(TRIANGLE3,    P2);
  psi = ShapeFunction::create(TRIANGLE3,    P0);

  //std::cout << phi->numDofsAssociatedToVertice() << " " <<
  //             phi->numDofsAssociatedToCorner()  << " " <<
  //             phi->numDofsAssociatedToFacet()   << " " <<
  //             phi->numDofsAssociatedToCell() << std::endl;

  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);
  
  DofHandler DofH(mesh);
  
  DofH.addVariable("velo",   phi, 2);
  DofH.addVariable("press",  psi, 1);
  
  //MeshTools::removeCell(mesh->getCellPtr(0), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(1), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(10), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(6), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(7), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(9), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(2), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(3), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(12), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(4), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(5), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(8), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(11), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(13), mesh);
  
  DofH.SetUp();
  //std::cout << "num dofs = " << DofH.numDofs() << std::endl;
  
  
  //EXPECT_EQ(12, DofH.getVariable(0).numDofs());
  //EXPECT_EQ( 4, DofH.getVariable(1).numDofs());
  //////
  //EXPECT_EQ(16, DofH.numDofs());  
  
  //bool relations[2][2] = {{1,1}, {1,1}};
  //DofH.setVariablesRelationship(*relations);
  
  //DofH.metisRenumber();
  //DofH.boostMinimumDegreeRenumber();
  //DofH.boostCuthillMcKeeRenumber();
  DofH.CuthillMcKeeRenumber();
  
  //ArrayXi var_cell_dofs(8);
  //ArrayXi var_cell_dofs2(3);
  //
  //DofH.getVariable(0).getCellDofs(var_cell_dofs.data(), mesh->getCellPtr(13));
  //DofH.getVariable(1).getCellDofs(var_cell_dofs2.data(), mesh->getCellPtr(13));
  //
  //std::cout << var_cell_dofs.transpose() << std::endl;
  //std::cout << var_cell_dofs2.transpose() << std::endl;

  
  DofH.printSparsityMatlab("Jac");
  vtk_printer.attachMesh(mesh);
  vtk_printer.writeVtk(mesh_out);
  
  // checking dofs
  std::vector<int> dat;
  
  EXPECT_EQ(2, DofH.numVars());
  
  getAllDofs(dat,DofH,mesh);

  std::sort(dat.begin(), dat.end());

  EXPECT_TRUE(dat[0] >= 0); // grau de liberdade definidos em todo lugar

  dat.erase(std::remove(dat.begin(), dat.end(), -1), dat.end()); // remove -1
  dat.erase(std::unique(dat.begin(), dat.end()), dat.end());     // remove duplicated
  
  EXPECT_EQ(DofH.numDofs(), (int)dat.size());
  
  int counter = 0;
  for (unsigned i = 0; i < dat.size(); ++i)
  {
    EXPECT_EQ(counter, dat[i]);
    ++counter;
  }  
  
  delete phi;
  delete psi;
  delete mesh;
}

TEST(DofHandlerTest, BubbleTet4)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh = NULL;  

  ECellType cell_t     = TRIANGLE3;
  const char* mesh_in  = "meshes/simptri3.msh";  
  const char* mesh_out = "meshes/outtest/metis_tet4.vtk";

  ShapeFunction *phi, *psi;
  phi = ShapeFunction::create(TETRAHEDRON4,    P1ph);
  psi = ShapeFunction::create(TETRAHEDRON4,    P1);

  //std::cout << phi->numDofsAssociatedToVertice() << " " <<
  //             phi->numDofsAssociatedToCorner()  << " " <<
  //             phi->numDofsAssociatedToFacet()   << " " <<
  //             phi->numDofsAssociatedToCell() << std::endl;

  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);
  
  DofHandler DofH(mesh);
  
  DofH.addVariable("velo",   phi, 2);
  DofH.addVariable("press",  psi, 1);
  
  //MeshTools::removeCell(mesh->getCellPtr(0), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(1), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(10), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(6), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(7), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(9), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(2), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(3), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(12), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(4), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(5), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(8), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(11), mesh);
  //MeshTools::removeCell(mesh->getCellPtr(13), mesh);
  
  DofH.SetUp();
  //std::cout << "num dofs = " << DofH.numDofs() << std::endl;
  
  
  //EXPECT_EQ(12, DofH.getVariable(0).numDofs());
  //EXPECT_EQ( 4, DofH.getVariable(1).numDofs());
  //////
  //EXPECT_EQ(16, DofH.numDofs());  
  
  //bool relations[2][2] = {{1,1}, {1,0}};
  //DofH.setVariablesRelationship(*relations);
  
  //DofH.metisRenumber();
  //DofH.boostMinimumDegreeRenumber();
  //DofH.boostCuthillMcKeeRenumber();
  DofH.CuthillMcKeeRenumber();
  
  //ArrayXi var_cell_dofs(8);
  //ArrayXi var_cell_dofs2(3);
  //
  //DofH.getVariable(0).getCellDofs(var_cell_dofs.data(), mesh->getCellPtr(13));
  //DofH.getVariable(1).getCellDofs(var_cell_dofs2.data(), mesh->getCellPtr(13));
  //
  //std::cout << var_cell_dofs.transpose() << std::endl;
  //std::cout << var_cell_dofs2.transpose() << std::endl;

  
  DofH.printSparsityMatlab("Jac");
  vtk_printer.attachMesh(mesh);
  vtk_printer.writeVtk(mesh_out);
  
  // checking dofs
  std::vector<int> dat;
  
  EXPECT_EQ(2, DofH.numVars());
  
  getAllDofs(dat,DofH,mesh);

  std::sort(dat.begin(), dat.end());

  EXPECT_TRUE(dat[0] >= 0); // grau de liberdade definidos em todo lugar

  dat.erase(std::remove(dat.begin(), dat.end(), -1), dat.end()); // remove -1
  dat.erase(std::unique(dat.begin(), dat.end()), dat.end());     // remove duplicated
  
  EXPECT_EQ(DofH.numDofs(), (int)dat.size());
  
  int counter = 0;
  for (unsigned i = 0; i < dat.size(); ++i)
  {
    EXPECT_EQ(counter, dat[i]);
    ++counter;
  }  
  
  delete phi;
  delete psi;
  delete mesh;
}



//
//
//
// TAGS TEST
//
//
// 



TEST(DofHandlerTest, TagsDofsTri3)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh = NULL;  

  ECellType cell_t     = TRIANGLE3;
  const char* mesh_in  = "meshes/simptri3.msh";  
  //const char* mesh_out = "meshes/outtest/dof_tri3.vtk";

  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);
  
  int ntags=1;
  int tags=1;
  
  DofHandler DofH(mesh);
  //                         ndpv,  ndpr,  ndpf,  ndpc
  DofH.addVariable("altura",    1,     0,     0,     0, ntags, &tags); // 4
  DofH.addVariable("vetor",     2,     0,     0,     0, ntags, &tags); // 8
  DofH.addVariable("foo",       0,     1,     0,     0, ntags, &tags); // 0
  DofH.addVariable("bar",       0,     0,     1,     0, ntags, &tags); // 4
  DofH.addVariable("moo",       0,     0,     0,     1, ntags, &tags); // 0
  
  
  DofH.SetUp();
  
  EXPECT_EQ(16, DofH.numDofs());
  
  
  MeshTools::removeCell(mesh->getCellPtr(2), mesh);
  MeshTools::removeCell(mesh->getCellPtr(3), mesh);
  
  DofH.SetUp();
  
  // 15
  EXPECT_EQ( 3, DofH.getVariable(0).numPositiveDofs());
  EXPECT_EQ( 6, DofH.getVariable(1).numPositiveDofs());
  EXPECT_EQ( 0, DofH.getVariable(2).numPositiveDofs());
  EXPECT_EQ( 3, DofH.getVariable(3).numPositiveDofs());
  EXPECT_EQ( 0, DofH.getVariable(4).numPositiveDofs());

  EXPECT_EQ(12, DofH.numDofs());
  
  
  // checking dofs
  std::vector<int> dat;
  
  EXPECT_EQ(5, DofH.numVars());
  
  getAllDofs(dat,DofH,mesh);

  std::sort(dat.begin(), dat.end());

  EXPECT_EQ(-1,dat[0]); // tags

  dat.erase(std::remove(dat.begin(), dat.end(), -1), dat.end()); // remove -1
  dat.erase(std::unique(dat.begin(), dat.end()), dat.end());     // remove duplicated
  
  EXPECT_EQ(DofH.numDofs(), (int)dat.size());
  
  int counter = 0;
  for (unsigned i = 0; i < dat.size(); ++i)
  {
    EXPECT_EQ(counter, dat[i]);
    ++counter;
  }  
  
  
  delete mesh;
}


TEST(DofHandlerTest, TagsDofsTet10)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh = NULL;  

  ECellType cell_t     = TETRAHEDRON10;
  const char* mesh_in  = "meshes/simptet10.msh";  
  const char* mesh_out = "meshes/outtest/dof_tet10.vtk";

  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);

  
  int ntags=0;
  int tags=0;
  
  DofHandler DofH(mesh);
  //                         ndpv,  ndpr,  ndpf,  ndpc
  DofH.addVariable("numero",    1,     0,     0,     0, ntags, &tags); //  9
  DofH.addVariable("foo",       0,     1,     0,     0, ntags, &tags); // 26
  DofH.addVariable("bar",       0,     0,     1,     0, ntags, &tags); // 30
  DofH.addVariable("moo",       0,     0,     0,     2, ntags, &tags); // 24

  DofH.SetUp();
  
  EXPECT_EQ( 9, DofH.getVariable(0).numPositiveDofs());
  EXPECT_EQ(26, DofH.getVariable(1).numPositiveDofs());
  EXPECT_EQ(30, DofH.getVariable(2).numPositiveDofs());
  EXPECT_EQ(24, DofH.getVariable(3).numPositiveDofs());
  
  EXPECT_EQ(89, DofH.numDofs());
  
  //DofH.SetUp();
  //
  MeshTools::removeCell(mesh->getCellPtr(2), mesh);
  MeshTools::removeCell(mesh->getCellPtr(3), mesh);
  //
  DofH.SetUp();
  //
  
  EXPECT_EQ( 9, DofH.getVariable(0).numPositiveDofs());
  EXPECT_EQ(26, DofH.getVariable(1).numPositiveDofs());
  EXPECT_EQ(28, DofH.getVariable(2).numPositiveDofs());
  EXPECT_EQ(20, DofH.getVariable(3).numPositiveDofs());
  ////
  EXPECT_EQ(83, DofH.numDofs());
  //
  //
  //// .getVariable(0)
  //int *dat = DofH.data();
  //
  //int counter = 0;
  //for (int i = 0; i < DofH.totalSize(); ++i)
  //{
  //  if ((*dat) >= 0)
  //    EXPECT_EQ(counter++, (*dat));
  //    
  //}
  //
  vtk_printer.attachMesh(mesh);
  vtk_printer.writeVtk(mesh_out);
  //
  std::vector<int> dados(DofH.numDofs());
  
  point_iterator pt = mesh->pointBegin();
  point_iterator pt_end = mesh->pointEnd();
  
  // for variable "0" ("numero")
  for (; pt != pt_end; ++pt)
  {
    if (!mesh->isVertex(&*pt))
      continue;
      
    int dof(0); // "numero"
    DofH.getVariable(0).getVertexDofs(&dof, &*pt);
    dados[dof] = mesh->getNodeContigId(mesh->getPointId(&*pt));
    //dados[dof] = 0;
  }
  
  MyGetDataVtkTri6Test get_data(dados.data(), mesh, &DofH);
  
  vtk_printer.addNodeIntVtk(DofH.getVariable(0).getName(), get_data);


  // checking dofs
  std::vector<int> dat;
  
  EXPECT_EQ(4, DofH.numVars());
  
  getAllDofs(dat,DofH,mesh);

  std::sort(dat.begin(), dat.end());

  EXPECT_TRUE(dat[0] == -1); // -1 nos nohs de controle

  dat.erase(std::remove(dat.begin(), dat.end(), -1), dat.end()); // remove -1
  dat.erase(std::unique(dat.begin(), dat.end()), dat.end());     // remove duplicated
  
  EXPECT_EQ(DofH.numDofs(), (int)dat.size());
  
  int counter = 0;
  for (unsigned i = 0; i < dat.size(); ++i)
  {
    EXPECT_EQ(counter, dat[i]);
    ++counter;
  }  

  
  delete mesh;
}



TEST(DofHandlerTest, TagsLinkTri3)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh = NULL;  

  ECellType cell_t     = TRIANGLE3;
  const char* mesh_in  = "meshes/simptri3.msh";  
  //const char* mesh_out = "meshes/outtest/dof_tri3.vtk";

  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);
  
  int ntags=1;
  int tags=1;
  
  DofHandler DofH(mesh);
  //                         ndpv,  ndpr,  ndpf,  ndpc
  DofH.addVariable("altura",    1,     0,     0,     0, ntags, &tags);
  DofH.addVariable("vetor",     2,     0,     0,     0, ntags, &tags);
  
  
  DofH.SetUp();
  
  EXPECT_EQ(12, DofH.numDofs());

  // .getVariable(0)
  //int *dat = DofH.data();

#define PRINT_DAT                            \
  for (int i = 0; i < DofH.totalSize(); ++i) \
  {                                          \
    std::cout.width (3);                     \
    std::cout << dat[i] << "\t";             \
    if (i==11)                               \
      std::cout <<"|\t";                     \
  }                                          \
  std::cout << std::endl;


  //PRINT_DAT
  
  int dofs1[3] = {-1,-1,-1};
  int dofs2[3] = {-1,-1,-1};
  
  DofH.getVariable(0).getVertexDofs(dofs1,   mesh->getNodePtr(7));
  DofH.getVariable(1).getVertexDofs(dofs1+1, mesh->getNodePtr(7));

  DofH.getVariable(0).getVertexDofs(dofs2,   mesh->getNodePtr(0));
  DofH.getVariable(1).getVertexDofs(dofs2+1, mesh->getNodePtr(0));

  DofH.linkDofs(3, dofs1, dofs2);

  EXPECT_EQ(9, DofH.numDofs()) << dofs1[0] << " "<< dofs1[1] << " "<< dofs1[2] << "\n"<< dofs2[0] << " "<< dofs2[1] << " "<< dofs2[2] << " ";


  // ckeck
  std::vector<int> dat;
  
  EXPECT_EQ(2, DofH.numVars());
  
  getAllDofs(dat,DofH,mesh);

  std::sort(dat.begin(), dat.end());

  EXPECT_TRUE(dat[0] == -1); // grau de liberdade nao definidos em todo lugar

  dat.erase(std::remove(dat.begin(), dat.end(), -1), dat.end()); // remove -1
  dat.erase(std::unique(dat.begin(), dat.end()), dat.end());     // remove duplicated
  
  EXPECT_EQ(DofH.numDofs(), (int)dat.size());
  
  int counter = 0;
  for (unsigned i = 0; i < dat.size(); ++i)
  {
    EXPECT_EQ(counter, dat[i]);
    ++counter;
  }
  
  //PRINT_DAT
  
  delete mesh;
}


TEST(DofHandlerTest, TagsLink2Tri3)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh = NULL;  

  ECellType cell_t     = TRIANGLE3;
  const char* mesh_in  = "meshes/simptri3.msh";  
  //const char* mesh_out = "meshes/outtest/dof_tri3.vtk";

  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);
  
  int ntags=1;
  int tags=1;
  
  DofHandler DofH(mesh);
  //                         ndpv,  ndpr,  ndpf,  ndpc
  DofH.addVariable("altura",    1,     0,     0,     0, ntags, &tags);
  DofH.addVariable("vetor",     2,     0,     0,     0, ntags, &tags);
  
  
  DofH.SetUp();
  
  EXPECT_EQ(12, DofH.numDofs());



  //PRINT_DAT
  int dofs1[3] = {-1,-1,-1};
  int dofs2[3] = {-1,-1,-1};
  
  //DofH.getVariable(0).getVertexDofs(dofs1,   mesh->getNodePtr(7));
  //DofH.getVariable(1).getVertexDofs(dofs1+1, mesh->getNodePtr(7));
  //
  //DofH.getVariable(0).getVertexDofs(dofs2,   mesh->getNodePtr(0));
  //DofH.getVariable(1).getVertexDofs(dofs2+1, mesh->getNodePtr(0));

  DofH.linkDofs(3, dofs1, dofs2); // do nothing

  EXPECT_EQ(12, DofH.numDofs()) << dofs1[0] << " "<< dofs1[1] << " "<< dofs1[2] << "\n"<< dofs2[0] << " "<< dofs2[1] << " "<< dofs2[2] << " ";


  // ckeck
  std::vector<int> dat;
  
  EXPECT_EQ(2, DofH.numVars());
  
  getAllDofs(dat,DofH,mesh);

  std::sort(dat.begin(), dat.end());

  EXPECT_TRUE(dat[0] == -1); // grau de liberdade nao definidos em todo lugar

  dat.erase(std::remove(dat.begin(), dat.end(), -1), dat.end()); // remove -1
  dat.erase(std::unique(dat.begin(), dat.end()), dat.end());     // remove duplicated
  
  EXPECT_EQ(DofH.numDofs(), (int)dat.size());
  
  int counter = 0;
  for (unsigned i = 0; i < dat.size(); ++i)
  {
    EXPECT_EQ(counter, dat[i]);
    ++counter;
  }
  
  //PRINT_DAT
  
  delete mesh;
}


TEST(DofHandlerTest, TagsLink3Tri3)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh = NULL;  

  ECellType cell_t     = TRIANGLE6;
  const char* mesh_in  = "meshes/simptri6.msh";  
  const char* mesh_out = "meshes/outtest/dof_tri6.vtk";

  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);
  
  DofHandler DofH(mesh);
  //                         ndpv,  ndpr,  ndpf,  ndpc
  DofH.addVariable("numero",    1,     0,     0,     0); // 12
  DofH.addVariable("vetor",     2,     0,     0,     0); // 24
  DofH.addVariable("foo",       0,     1,     0,     0); //  0
  DofH.addVariable("bar",       0,     0,     1,     0); // 25
  DofH.addVariable("moo",       0,     0,     0,     1); // 14
  
  DofH.SetUp();
  
  EXPECT_EQ(75, DofH.numDofs());
  
  MeshTools::removeCell(mesh->getCellPtr(2), mesh);
  MeshTools::removeCell(mesh->getCellPtr(3), mesh);
  
  DofH.SetUp();
  
  //11x3 + 22 + 12 + 12
  EXPECT_EQ(11, DofH.getVariable(0).numPositiveDofs());
  EXPECT_EQ(22, DofH.getVariable(1).numPositiveDofs());
  EXPECT_EQ( 0, DofH.getVariable(2).numPositiveDofs());
  EXPECT_EQ(22, DofH.getVariable(3).numPositiveDofs());
  EXPECT_EQ(12, DofH.getVariable(4).numPositiveDofs());

  EXPECT_EQ(67, DofH.numDofs());
  

  // .getVariable(0)
  //int *dat = DofH.data();

  int dofs1[] = {-1,-1,-1,-1};
  int dofs2[] = {-1,-1,-1,-1};
  
  DofH.getVariable(0).getVertexDofs(dofs1,   mesh->getNodePtr(4));
  DofH.getVariable(1).getVertexDofs(dofs1+1, mesh->getNodePtr(4));

  DofH.getVariable(0).getVertexDofs(dofs2,   mesh->getNodePtr(13));
  DofH.getVariable(1).getVertexDofs(dofs2+1, mesh->getNodePtr(13));

  DofH.getVariable(4).getCellAssociatedDofs(dofs1+3, mesh->getCellPtr(10));
  DofH.getVariable(4).getCellAssociatedDofs(dofs2+3, mesh->getCellPtr(11));

  DofH.linkDofs(4, dofs1, dofs2);

  EXPECT_EQ(63, DofH.numDofs()) << dofs1[0] << " "<< dofs1[1] << " "<< dofs1[2] << " "<< dofs1[3] << "\n"<< dofs2[0] << " "<< dofs2[1] << " "<< dofs2[2] << " "<< dofs2[3] << " ";


  // ckeck
  std::vector<int> dat;
  
  EXPECT_EQ(5, DofH.numVars());
  
  getAllDofs(dat,DofH,mesh);

  std::sort(dat.begin(), dat.end());

  EXPECT_TRUE(dat[0] == -1);

  dat.erase(std::remove(dat.begin(), dat.end(), -1), dat.end()); // remove -1
  dat.erase(std::unique(dat.begin(), dat.end()), dat.end());     // remove duplicated
  
  EXPECT_EQ(DofH.numDofs(), (int)dat.size());
  
  int counter = 0;
  for (unsigned i = 0; i < dat.size(); ++i)
  {
    EXPECT_EQ(counter, dat[i]);
    ++counter;
  }
  
  //PRINT_DAT
  
  
  delete mesh;
}


TEST(DofHandlerTest, CopyFunctionTri3)
{
  

  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh = NULL;  

  ECellType cell_t     = TRIANGLE3;
  const char* mesh_in  = "meshes/simptri3.msh";  
  //const char* mesh_out = "meshes/outtest/dof_tri3.vtk";

  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);
  
  DofHandler DofH;
  
  {
    DofHandler DofH_temp(mesh);
    //                               ndpv,  ndpr,  ndpf,  ndpc
    DofH_temp.addVariable("altura",    1,     0,     0,     0);
    DofH_temp.addVariable("vetor",     2,     0,     0,     0);
    
    DofH_temp.SetUp();
    
    DofH.copy(DofH_temp);
  }
  
  EXPECT_EQ(36, DofH.numDofs());

  // .getVariable(0)
  //int *dat = DofH.data();


#define PRINT_DAT                            \
  for (int i = 0; i < DofH.totalSize(); ++i) \
  {                                          \
    std::cout.width (3);                     \
    std::cout << dat[i] << "\t";             \
    if (i==11)                               \
      std::cout <<"|\t";                     \
  }                                          \
  std::cout << std::endl;


  //PRINT_DAT
  
  int dofs1[3] = {-1,-1,-1};
  int dofs2[3] = {-1,-1,-1};
  
  DofH.getVariable(0).getVertexDofs(dofs1,   mesh->getNodePtr(2));
  DofH.getVariable(1).getVertexDofs(dofs1+1, mesh->getNodePtr(2));

  DofH.getVariable(0).getVertexDofs(dofs2,   mesh->getNodePtr(6));
  DofH.getVariable(1).getVertexDofs(dofs2+1, mesh->getNodePtr(6));

  DofH.linkDofs(3, dofs1, dofs2);

  EXPECT_EQ(33, DofH.numDofs()) << dofs1[0] << " "<< dofs1[1] << " "<< dofs1[2] << "\n"<< dofs2[0] << " "<< dofs2[1] << " "<< dofs2[2] << " ";


  // ckeck
  std::vector<int> dat;
  
  EXPECT_EQ(2, DofH.numVars());
  
  getAllDofs(dat,DofH,mesh);

  std::sort(dat.begin(), dat.end());

  EXPECT_TRUE(dat[0] >= 0); // grau de liberdade definidos em todo lugar

  dat.erase(std::remove(dat.begin(), dat.end(), -1), dat.end()); // remove -1
  dat.erase(std::unique(dat.begin(), dat.end()), dat.end());     // remove duplicated
  
  EXPECT_EQ(DofH.numDofs(), (int)dat.size());
  
  int counter = 0;
  for (unsigned i = 0; i < dat.size(); ++i)
  {
    EXPECT_EQ(counter, dat[i]);
    ++counter;
  }

  
  
  delete mesh;  
}



//
//
// SPLITTING AT INTERSECTION
//

TEST(DofHandlerTest, SplitAtIntersectionVolTri3)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh = NULL;  

  ECellType cell_t     = TRIANGLE3;
  const char* mesh_in  = "meshes/1level_tri3.msh";  
  const char* mesh_out = "meshes/outtest/1level_tri3.vtk";

  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);
  
  vtk_printer.attachMesh(mesh);
  vtk_printer.writeVtk(mesh_out);
//  vtk_printer.printPointIcellVtk();
//  vtk_printer.printPointPositionVtk();  
  
  DofHandler DofH(mesh);
  //                         ndpv,  ndpr,  ndpf,  ndpc
  DofH.addVariable("altura",    1,     0,     0,     0); // 16
  DofH.addVariable("vetor",     0,     0,     1,     0); // 18
  DofH.addVariable("coisa",     0,     0,     0,     1); // 33
  
  DofH.SetUp();
  
  EXPECT_EQ(67, DofH.numDofs());


  // split
  int regs[] = {71,70};
  int nr = 2;
  DofH.getVariable(0).setType(SPLITTED_BY_REGION_CELL, 0, 0);
  DofH.getVariable(1).setType(SPLITTED_BY_REGION_CELL, nr, regs);
  DofH.getVariable(2).setType(SPLITTED_BY_REGION_CELL, nr, regs);
  
  DofH.SetUp();

  EXPECT_EQ(74, DofH.numDofs());

  // .getVariable(0)
  //int *dat = DofH.data();

  //PRINT_DAT


  // ckeck
  std::vector<int> dat;
  
  EXPECT_EQ(3, DofH.numVars());
  
  getAllDofs(dat,DofH,mesh);
  
  std::sort(dat.begin(), dat.end());
  
  EXPECT_TRUE(dat[0] == 0);
  
  dat.erase(std::remove(dat.begin(), dat.end(), -1), dat.end()); // remove -1
  dat.erase(std::unique(dat.begin(), dat.end()), dat.end());     // remove duplicated
  
  EXPECT_EQ(DofH.numDofs(), (int)dat.size());
  
  int counter = 0;
  for (unsigned i = 0; i < dat.size(); ++i)
  {
    EXPECT_EQ(counter, dat[i]);
    ++counter;
  }
  
  //PRINT_DAT
  
  delete mesh;
}


