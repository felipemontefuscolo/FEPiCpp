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
//  MeshTools::removeCell(mesh->getCell(0), mesh);
//  MeshTools::removeCell(mesh->getCell(2), mesh);
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
TEST(AssignsDofsTest, WithTri3)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
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
    l = mesh->getNode(k);
    DofH.getVariable(1).getVertexDofs(dofs, l);
    nn = count_if(dofs, dofs_end, bind2nd(greater_equal<int>(),0) ); // non-negative
    EXPECT_EQ(2, nn) << "node id = " << mesh->getCell(l->getIncidCell())->getNodeId(l->getPosition());;
    fill(dofs, dofs_end, -1);
  }
  for (int k = 0; k < mesh->numNodesTotal(); ++k)
  {
    l = mesh->getNode(k);
    DofH.getVariable(1).getCornerDofs(dofs, l);
    nn = count_if(dofs, dofs_end, bind2nd(greater_equal<int>(),0) ); // non-negative
    EXPECT_EQ(2, nn) << "node id = " << mesh->getCell(l->getIncidCell())->getNodeId(l->getPosition());;
    fill(dofs, dofs_end, -1);
  }
  
  
  MeshTools::removeCell(mesh->getCell(2), mesh);
  MeshTools::removeCell(mesh->getCell(3), mesh);
  
  DofH.SetUp();
  
  //11x3 + 22 + 12
  EXPECT_EQ(11, DofH.getVariable(0).numDofs());
  EXPECT_EQ(22, DofH.getVariable(1).numDofs());
  EXPECT_EQ( 0, DofH.getVariable(2).numDofs());
  EXPECT_EQ(22, DofH.getVariable(3).numDofs());
  EXPECT_EQ(12, DofH.getVariable(4).numDofs());

  EXPECT_EQ(67, DofH.numDofs());
  
  
  // .getVariable(0)
  int *dat = DofH.data();

  int counter = 0;
  for (int i = 0; i < DofH.totalSize(); ++i)
  {
    if ((*dat) >= 0)
      EXPECT_EQ(counter++, (*dat));
      
    //std::cout << (*dat++) << std::endl;
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
    Point const* p = mesh_ptr->getNode(nodeid);
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

TEST(AssignsDofsTest, WithTri6)
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
  
  MeshTools::removeCell(mesh->getCell(2), mesh);
  MeshTools::removeCell(mesh->getCell(3), mesh);
  
  DofH.SetUp();
  
  //11x3 + 22 + 12 + 12
  EXPECT_EQ(11, DofH.getVariable(0).numDofs());
  EXPECT_EQ(22, DofH.getVariable(1).numDofs());
  EXPECT_EQ( 0, DofH.getVariable(2).numDofs());
  EXPECT_EQ(22, DofH.getVariable(3).numDofs());
  EXPECT_EQ(12, DofH.getVariable(4).numDofs());

  EXPECT_EQ(67, DofH.numDofs());
  
  
  // .getVariable(0)
  int *dat = DofH.data();

  int counter = 0;
  for (int i = 0; i < DofH.totalSize(); ++i)
  {
    if ((*dat) >= 0)
      EXPECT_EQ(counter++, (*dat));
      
  }
  
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
  
  delete mesh;
}

TEST(AssignsDofsTest, WithTet10)
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
  
  EXPECT_EQ( 9, DofH.getVariable(0).numDofs());
  EXPECT_EQ(26, DofH.getVariable(1).numDofs());
  EXPECT_EQ(30, DofH.getVariable(2).numDofs());
  EXPECT_EQ(24, DofH.getVariable(3).numDofs());
  
  EXPECT_EQ(89, DofH.numDofs());
  
  //DofH.SetUp();
  //
  MeshTools::removeCell(mesh->getCell(2), mesh);
  MeshTools::removeCell(mesh->getCell(3), mesh);
  //
  DofH.SetUp();
  //
  
  EXPECT_EQ( 9, DofH.getVariable(0).numDofs());
  EXPECT_EQ(26, DofH.getVariable(1).numDofs());
  EXPECT_EQ(28, DofH.getVariable(2).numDofs());
  EXPECT_EQ(20, DofH.getVariable(3).numDofs());
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
  
  delete mesh;
}

TEST(BubbleTri3Test, WithTri3)
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
  
  //MeshTools::removeCell(mesh->getCell(0), mesh);
  //MeshTools::removeCell(mesh->getCell(1), mesh);
  //MeshTools::removeCell(mesh->getCell(10), mesh);
  //MeshTools::removeCell(mesh->getCell(6), mesh);
  //MeshTools::removeCell(mesh->getCell(7), mesh);
  //MeshTools::removeCell(mesh->getCell(9), mesh);
  //MeshTools::removeCell(mesh->getCell(2), mesh);
  //MeshTools::removeCell(mesh->getCell(3), mesh);
  //MeshTools::removeCell(mesh->getCell(12), mesh);
  //MeshTools::removeCell(mesh->getCell(4), mesh);
  //MeshTools::removeCell(mesh->getCell(5), mesh);
  //MeshTools::removeCell(mesh->getCell(8), mesh);
  //MeshTools::removeCell(mesh->getCell(11), mesh);
  //MeshTools::removeCell(mesh->getCell(13), mesh);
  
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
  //DofH.getVariable(0).getCellDofs(var_cell_dofs.data(), mesh->getCell(13));
  //DofH.getVariable(1).getCellDofs(var_cell_dofs2.data(), mesh->getCell(13));
  //
  //std::cout << var_cell_dofs.transpose() << std::endl;
  //std::cout << var_cell_dofs2.transpose() << std::endl;

  
  DofH.printSparsityMatlab("Jac");
  vtk_printer.attachMesh(mesh);
  vtk_printer.writeVtk(mesh_out);
  
  std::sort(DofH.data(), DofH.data()+DofH.totalSize());
  int counter = 0;
  for (int i = 0; i < DofH.totalSize(); ++i)
  {
    if (DofH.data() < 0) continue;
    EXPECT_EQ(counter, DofH.data()[i]);
    counter++;
  }
  
  delete phi;
  delete psi;
  delete mesh;
}

TEST(BubbleTet4Test, WithTet4)
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
  
  //MeshTools::removeCell(mesh->getCell(0), mesh);
  //MeshTools::removeCell(mesh->getCell(1), mesh);
  //MeshTools::removeCell(mesh->getCell(10), mesh);
  //MeshTools::removeCell(mesh->getCell(6), mesh);
  //MeshTools::removeCell(mesh->getCell(7), mesh);
  //MeshTools::removeCell(mesh->getCell(9), mesh);
  //MeshTools::removeCell(mesh->getCell(2), mesh);
  //MeshTools::removeCell(mesh->getCell(3), mesh);
  //MeshTools::removeCell(mesh->getCell(12), mesh);
  //MeshTools::removeCell(mesh->getCell(4), mesh);
  //MeshTools::removeCell(mesh->getCell(5), mesh);
  //MeshTools::removeCell(mesh->getCell(8), mesh);
  //MeshTools::removeCell(mesh->getCell(11), mesh);
  //MeshTools::removeCell(mesh->getCell(13), mesh);
  
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
  //DofH.getVariable(0).getCellDofs(var_cell_dofs.data(), mesh->getCell(13));
  //DofH.getVariable(1).getCellDofs(var_cell_dofs2.data(), mesh->getCell(13));
  //
  //std::cout << var_cell_dofs.transpose() << std::endl;
  //std::cout << var_cell_dofs2.transpose() << std::endl;

  
  DofH.printSparsityMatlab("Jac");
  vtk_printer.attachMesh(mesh);
  vtk_printer.writeVtk(mesh_out);
  
  std::sort(DofH.data(), DofH.data()+DofH.totalSize());
  int counter = 0;
  for (int i = 0; i < DofH.totalSize(); ++i)
  {
    if (DofH.data() < 0) continue;
    EXPECT_EQ(counter, DofH.data()[i]);
    counter++;
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



TEST(TagsDofsTest, WithTri3)
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
  
  
  MeshTools::removeCell(mesh->getCell(2), mesh);
  MeshTools::removeCell(mesh->getCell(3), mesh);
  
  DofH.SetUp();
  
  // 15
  EXPECT_EQ( 3, DofH.getVariable(0).numDofs());
  EXPECT_EQ( 6, DofH.getVariable(1).numDofs());
  EXPECT_EQ( 0, DofH.getVariable(2).numDofs());
  EXPECT_EQ( 3, DofH.getVariable(3).numDofs());
  EXPECT_EQ( 0, DofH.getVariable(4).numDofs());

  EXPECT_EQ(12, DofH.numDofs());
  
  
  // .getVariable(0)
  int *dat = DofH.data();

  int counter = 0;
  for (int i = 0; i < DofH.totalSize(); ++i)
  {
    if ((*dat) >= 0)
      EXPECT_EQ(counter++, (*dat));
      
    //std::cout << (*dat++) << std::endl;
  }
  
  
  delete mesh;
}

TEST(TagsDofsTest, WithTet10)
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
  
  EXPECT_EQ( 9, DofH.getVariable(0).numDofs());
  EXPECT_EQ(26, DofH.getVariable(1).numDofs());
  EXPECT_EQ(30, DofH.getVariable(2).numDofs());
  EXPECT_EQ(24, DofH.getVariable(3).numDofs());
  
  EXPECT_EQ(89, DofH.numDofs());
  
  //DofH.SetUp();
  //
  MeshTools::removeCell(mesh->getCell(2), mesh);
  MeshTools::removeCell(mesh->getCell(3), mesh);
  //
  DofH.SetUp();
  //
  
  EXPECT_EQ( 9, DofH.getVariable(0).numDofs());
  EXPECT_EQ(26, DofH.getVariable(1).numDofs());
  EXPECT_EQ(28, DofH.getVariable(2).numDofs());
  EXPECT_EQ(20, DofH.getVariable(3).numDofs());
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
  
  delete mesh;
}

















