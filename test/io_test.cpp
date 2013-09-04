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

// --gtest_catch_exceptions=0



#include <gtest/gtest.h>
#include <Fepic/Mesh>
#include <tr1/array>
#include <tr1/tuple>
#include <algorithm>
#include <tr1/memory>


bool is_binary = false;


class GetDataVtk : public DefaultGetDataVtk
{
public:

  /**
   */  

  GetDataVtk() : DefaultGetDataVtk() {}
  //DefaultGetDataVtk(int * i = NULL, Real * r = NULL) : data_r(r), data_i(i) {}
  
  virtual Real get_data_r(int id) const
  {
    return 0.1234567890123456789012345678901234567890;
  }
  virtual int get_data_i(int id) const
  {
    return 1234567890;
  }
  
  virtual ~GetDataVtk() {}
  
};


TEST(IOTest, identifiesMshMeshTypeTest)
{
  MeshIoMsh msh_reader;
  int sd;
  ECellType ct;

  ct = msh_reader.identifiesMshMeshType("meshes/simpedge2.msh",sd);
  EXPECT_EQ(EDGE2, ct);
  EXPECT_EQ(2, sd);

  ct = msh_reader.identifiesMshMeshType("meshes/simpedge3.msh",sd);
  EXPECT_EQ(EDGE3, ct);
  EXPECT_EQ(2, sd);

  ct = msh_reader.identifiesMshMeshType("meshes/simptri3.msh",sd);
  EXPECT_EQ(TRIANGLE3, ct);
  EXPECT_EQ(2, sd);

  ct = msh_reader.identifiesMshMeshType("meshes/simptri6.msh",sd);
  EXPECT_EQ(TRIANGLE6, ct);
  EXPECT_EQ(2, sd);

  ct = msh_reader.identifiesMshMeshType("meshes/simpquad4.msh",sd);
  EXPECT_EQ(QUADRANGLE4, ct);
  EXPECT_EQ(2, sd);

  ct = msh_reader.identifiesMshMeshType("meshes/simpquad8.msh",sd);
  EXPECT_EQ(QUADRANGLE8, ct);
  EXPECT_EQ(2, sd);

  ct = msh_reader.identifiesMshMeshType("meshes/simpquad9.msh",sd);
  EXPECT_EQ(QUADRANGLE9, ct);
  EXPECT_EQ(2, sd);

  ct = msh_reader.identifiesMshMeshType("meshes/simptet4.msh",sd);
  EXPECT_EQ(TETRAHEDRON4, ct);
  EXPECT_EQ(3, sd);

  ct = msh_reader.identifiesMshMeshType("meshes/simptet4.msh",sd);
  EXPECT_EQ(TETRAHEDRON4, ct);
  EXPECT_EQ(3, sd);

  ct = msh_reader.identifiesMshMeshType("meshes/simptet10.msh",sd);
  EXPECT_EQ(TETRAHEDRON10, ct);
  EXPECT_EQ(3, sd);

  ct = msh_reader.identifiesMshMeshType("meshes/simphex8.msh",sd);
  EXPECT_EQ(HEXAHEDRON8, ct);
  EXPECT_EQ(3, sd);

  ct = msh_reader.identifiesMshMeshType("meshes/simphex20.msh",sd);
  EXPECT_EQ(HEXAHEDRON20, ct);
  EXPECT_EQ(3, sd);

  ct = msh_reader.identifiesMshMeshType("meshes/simphex27.msh",sd);
  EXPECT_EQ(HEXAHEDRON27, ct);
  EXPECT_EQ(3, sd);

}

TEST(IOTest, WriteVtkEdge2Test)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int dim = 2;

  mesh = Mesh::create(EDGE2,dim);
  msh_reader.readFileMsh("meshes/simpedge2.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.setBinaryOutput( is_binary);
  vtk_printer.writeVtk("meshes/outtest/simpedge2.vtk");
  vtk_printer.printPointTagVtk();
  vtk_printer.printPointIcellVtk();
  vtk_printer.printPointPositionVtk();  
  vtk_printer.addNodeScalarVtk("foo",GetDataVtk());
  vtk_printer.printCellIdVtk();
  
  delete mesh;
}

TEST(IOTest, WriteVtkEdge3Test)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int dim = 2;

  mesh = Mesh::create(EDGE3,dim);
  msh_reader.readFileMsh("meshes/simpedge3.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.setBinaryOutput( is_binary);
  vtk_printer.writeVtk("meshes/outtest/simpedge3.vtk");
  vtk_printer.printPointTagVtk();
  vtk_printer.printPointIcellVtk();
  vtk_printer.printPointPositionVtk();
  vtk_printer.addNodeScalarVtk("foo",GetDataVtk());
  vtk_printer.printCellIdVtk();

  delete mesh;
}

TEST(IOTest, WriteVtkTri3Test)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int dim = 2;

  mesh = Mesh::create(TRIANGLE3,dim);
  msh_reader.readFileMsh("meshes/simptri3.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.setBinaryOutput( is_binary);
  vtk_printer.writeVtk("meshes/outtest/simptri3.vtk");
  vtk_printer.printPointTagVtk();
  vtk_printer.printPointIcellVtk();
  vtk_printer.printPointPositionVtk();
  vtk_printer.addNodeScalarVtk("foo",GetDataVtk());
  vtk_printer.printCellIdVtk();

  delete mesh;
}

TEST(IOTest, WriteVtkTri6Test)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int dim = 2;

  mesh = Mesh::create(TRIANGLE6,dim);
  msh_reader.readFileMsh("meshes/simptri6.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.setBinaryOutput( is_binary);
  vtk_printer.writeVtk("meshes/outtest/simptri6.vtk");
  vtk_printer.printPointTagVtk();
  vtk_printer.printPointIcellVtk();
  vtk_printer.printPointPositionVtk();
  vtk_printer.addNodeScalarVtk("foo",GetDataVtk());
  vtk_printer.printCellIdVtk();

  delete mesh;
}

TEST(IOTest, WriteVtkQuad4Test)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int dim = 2;

  mesh = Mesh::create(QUADRANGLE4,dim);
  msh_reader.readFileMsh("meshes/simpquad4.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.setBinaryOutput( is_binary);
  vtk_printer.writeVtk("meshes/outtest/simpquad4.vtk");
  vtk_printer.printPointTagVtk();
  vtk_printer.printPointIcellVtk();
  vtk_printer.printPointPositionVtk();
  vtk_printer.addNodeScalarVtk("foo",GetDataVtk());
  vtk_printer.printCellIdVtk();
  
  //mesh->printInfo();
  
  delete mesh;
}

TEST(IOTest, WriteVtkQuad8Test)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int dim = 2;

  mesh = Mesh::create(QUADRANGLE8,dim);
  msh_reader.readFileMsh("meshes/simpquad8.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.setBinaryOutput( is_binary);
  vtk_printer.writeVtk("meshes/outtest/simpquad8.vtk");
  vtk_printer.printPointTagVtk();
  vtk_printer.printPointIcellVtk();
  vtk_printer.printPointPositionVtk();
  vtk_printer.addNodeScalarVtk("foo",GetDataVtk());
  vtk_printer.printCellIdVtk();

  delete mesh;
}

TEST(IOTest, WriteVtkQuad9Test)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int dim = 2;

  mesh = Mesh::create(QUADRANGLE9,dim);
  msh_reader.readFileMsh("meshes/simpquad9.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.setBinaryOutput( is_binary);
  vtk_printer.writeVtk("meshes/outtest/simpquad9.vtk");
  vtk_printer.printPointTagVtk();
  vtk_printer.printPointIcellVtk();
  vtk_printer.printPointPositionVtk();
  vtk_printer.addNodeScalarVtk("foo",GetDataVtk());
  vtk_printer.printCellIdVtk();
  
  delete mesh;
}

TEST(IOTest, WriteVtkTet4Test)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int dim = 3;

  mesh = Mesh::create(TETRAHEDRON4,dim);
  msh_reader.readFileMsh("meshes/simptet4.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.setBinaryOutput( is_binary);
  vtk_printer.writeVtk("meshes/outtest/simptet4.vtk");
  vtk_printer.printPointTagVtk();
  vtk_printer.printPointIcellVtk();
  vtk_printer.printPointPositionVtk();
  vtk_printer.addNodeScalarVtk("foo",GetDataVtk());
  vtk_printer.printCellIdVtk();
  
  delete mesh;



}

TEST(IOTest, WriteVtkTet10Test)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int dim = 3;

  mesh = Mesh::create(TETRAHEDRON10,dim);
  msh_reader.readFileMsh("meshes/simptet10.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.setBinaryOutput( is_binary);
  vtk_printer.writeVtk("meshes/outtest/simptet10.vtk");
  vtk_printer.printPointTagVtk();
  vtk_printer.printPointIcellVtk();
  vtk_printer.printPointPositionVtk();
  vtk_printer.addNodeScalarVtk("foo",GetDataVtk());
  vtk_printer.printCellIdVtk();

  delete mesh;
}

TEST(IOTest, WriteVtkHex8Test)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int dim = 3;

  mesh = Mesh::create(HEXAHEDRON8,dim);
  msh_reader.readFileMsh("meshes/simphex8.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.setBinaryOutput( is_binary);
  vtk_printer.writeVtk("meshes/outtest/simphex8.vtk");
  vtk_printer.printPointTagVtk();
  vtk_printer.printPointIcellVtk();
  vtk_printer.printPointPositionVtk();
  vtk_printer.addNodeScalarVtk("foo",GetDataVtk());
  vtk_printer.printCellIdVtk();
  
  delete mesh;
}

TEST(IOTest, WriteVtkHex20Test)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int dim = 3;

  mesh = Mesh::create(HEXAHEDRON20,dim);
  msh_reader.readFileMsh("meshes/simphex20.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.setBinaryOutput( is_binary);
  vtk_printer.writeVtk("meshes/outtest/simphex20.vtk");
  vtk_printer.printPointTagVtk();
  vtk_printer.printPointIcellVtk();
  vtk_printer.printPointPositionVtk();
  vtk_printer.addNodeScalarVtk("foo",GetDataVtk());
  vtk_printer.printCellIdVtk();
  
  delete mesh;
}

TEST(IOTest, WriteVtkHex27Test)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int dim = 3;

  mesh = Mesh::create(HEXAHEDRON27,dim);
  msh_reader.readFileMsh("meshes/simphex27.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.setBinaryOutput( is_binary);
  vtk_printer.writeVtk("meshes/outtest/simphex27.vtk");
  vtk_printer.printPointTagVtk();
  vtk_printer.printPointIcellVtk();
  vtk_printer.printPointPositionVtk();
  vtk_printer.addNodeScalarVtk("foo",GetDataVtk());
  vtk_printer.printCellIdVtk();

  delete mesh;
}


//TEST(IOTest, CheckMeshSizeTet4)
//{
//  MeshIoMsh msh_reader;
//  MeshIoVtk vtk_printer;
//  Mesh *mesh;
//  int dim = 3;
//
//  mesh = Mesh::create(TETRAHEDRON4,dim);
//  msh_reader.readFileMsh("meshes/bigtet4.msh", mesh);
//  vtk_printer.attachMesh(mesh);
//  vtk_printer.setBinaryOutput( is_binary);
//  vtk_printer.writeVtk("meshes/outtest/bigtet4.vtk");
//  vtk_printer.printPointTagVtk();
//  vtk_printer.printPointIcellVtk();
//  vtk_printer.printPointPositionVtk();
//  vtk_printer.addNodeScalarVtk("foo",GetDataVtk());
//  vtk_printer.printCellIdVtk();
//
//  delete mesh;
//}










