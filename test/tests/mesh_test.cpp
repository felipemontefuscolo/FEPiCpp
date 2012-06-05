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
#include <tr1/array>
#include <tr1/tuple>
#include <algorithm>


TEST(Edge2ReadVtkTest, ReadVtk)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int dim = 2;

  mesh = Mesh::create(EDGE2,dim);
  msh_reader.readFileMsh("meshes/simpedge2.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.writeVtk("meshes/out/simpedge2.vtk");  
  
  int icells[][2] = {{1,7},{2,0},{3,1},{4,2},{5,3},{6,4},{7,5},{0,6}};
  int ans[2];
  
  Edge2 *cell;
  
  for (int i = 0; i < mesh->numCells(); ++i)
  {
    cell = (Edge2*)mesh->getCell(i);
    ans[0] = cell->Edge2::getIncidCell(0);
    ans[1] = cell->Edge2::getIncidCell(1);
    
    EXPECT_TRUE(sameElements(icells[i], icells[i]+2, ans, ans+2))
      << "icells[0] = " << icells[i][0] << "; icells[1] = " << icells[i][1] << "\n"
      << "ans[0]    = " << ans[0]       << "; ans[1]    = " << ans[1];
  }
  
  delete mesh;
}


//TEST(EdgeRead3VtkTest, ReadVtk)
//{
  //MeshIoMsh msh_reader;
  //MeshIoVtk vtk_printer;
  //Mesh *mesh;
  //int dim = 2;

  //mesh = Mesh::create(EDGE3,dim);
  //msh_reader.readFileMsh("meshes/simpedge3.msh", mesh);
  //vtk_printer.attachMesh(mesh);
  //vtk_printer.writeVtk("meshes/out/simpedge3.vtk");  
    
  //delete mesh;
//}


class Tri3VertexStarTest : public testing::TestWithParam<int> { // vertex id
public:

  virtual void SetUp() {
    mesh = Mesh::create(TRIANGLE3);
    msh_reader.readFileMsh("meshes/simptri3.msh", mesh);
    vtk_printer.attachMesh(mesh);
    vtk_printer.writeVtk("meshes/out/simptri3.vtk");

    id = GetParam();

    p = mesh->getNode(id);

    n = (int)(mesh->vertexStar(p, iCs, viCs)-iCs);

  };

  virtual void TearDown()
  {
    delete mesh;
  }

  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], viCs[64];
  int n, id;
  Point *p;
};

int Tri3VerticesId[] = {0,1,2,3,4,5,6,7,8,9,10,11};

TEST_P(Tri3VertexStarTest, vertexStarTest) {


  // vertex star
  int s[12][6] = {{2,3,-1,-1,-1,-1},
                  {7,6,-1,-1,-1,-1},
                  {0,1,-1,-1,-1,-1},
                  {4,5,-1,-1,-1,-1},
                  {2,9,7,-1,-1,-1},
                  {3,4,12,-1,-1,-1},
                  {0,5,13,-1,-1,-1},
                  {1,10,6,-1,-1,-1},
                  {6,7,8,9,10,-1},
                  {0,1,8,10,11,13},
                  {4,5,11,12,13,-1},
                  {2,3,9,8,11,12}};

  int q[12] = {2,2,2,2,3,3,3,3,5,6,5,6};


  EXPECT_EQ(q[id], n);
  EXPECT_TRUE(sameElements(iCs,iCs+n,s[id],s[id]+n));

}
INSTANTIATE_TEST_CASE_P(VertexStar, Tri3VertexStarTest, ::testing::ValuesIn(Tri3VerticesId));




class Tri6NodeStarTest : public testing::TestWithParam<int>
{ // node id
public:

  virtual void SetUp() {
    mesh = Mesh::create(TRIANGLE6);
    msh_reader.readFileMsh("meshes/simptri6.msh", mesh);
    vtk_printer.attachMesh(mesh);
    vtk_printer.writeVtk("meshes/out/simptri6.vtk");

    id = GetParam();

    p = mesh->getNode(id);

    n = (int)(mesh->nodeStar(p, iCs, viCs)-iCs);

  };

  virtual void TearDown()
  {
    delete mesh;
  }

  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], viCs[64];
  int n, id;
  Point *p;
};

int Tri6NodesId[] = {2, 5, 23, 16};

TEST_P(Tri6NodeStarTest, nodeStarTest)
{

  // node star
  int s[24][6];
  s[ 2][0]=0; s[2][1]=1;
  s[ 5][0]=7;
  s[23][0]=2; s[23][1]=9;
  s[16][0]=6; s[16][1]=7; s[16][2]=9; s[16][3]=8; s[16][4]=10;

  int q[24];
  q[ 2] = 2;
  q[ 5] = 1;
  q[23] = 2;
  q[16] = 5;

  EXPECT_EQ(q[id], n);
  EXPECT_TRUE(sameElements(iCs,iCs+n,s[id],s[id]+n));

}
INSTANTIATE_TEST_CASE_P(NodeStar, Tri6NodeStarTest, ::testing::ValuesIn(Tri6NodesId));




class Quad4VertexStarTest : public testing::TestWithParam<int> { // vertex id
public:

  virtual void SetUp() {
    mesh = Mesh::create(QUADRANGLE4);
    msh_reader.readFileMsh("meshes/simpquad4.msh", mesh);
    vtk_printer.attachMesh(mesh);
    vtk_printer.writeVtk("meshes/out/simpquad4.vtk");

    id = GetParam();

    p = mesh->getNode(id);

    n = (int)(mesh->vertexStar(p, iCs, viCs)-iCs);

  };

  virtual void TearDown()
  {
    delete mesh;
  }

  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], viCs[64];
  int n, id;
  Point *p;
};

int Quad4VerticesId[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
//int Quad4NodesId[] = {0};

TEST_P(Quad4VertexStarTest, vertexStarTest) {


  // vertex star
  int s[16][4] = {{ 8,-1,-1,-1},
                  { 2,-1,-1,-1},
                  { 0,-1,-1,-1},
                  { 6,-1,-1,-1},
                  { 8, 5,-1,-1},
                  { 2, 5,-1,-1},
                  { 2, 1,-1,-1},
                  { 1, 0,-1,-1},
                  { 0, 3,-1,-1},
                  { 3, 6,-1,-1},
                  { 6, 7,-1,-1},
                  { 7, 8,-1,-1},
                  { 0, 1, 3, 4},
                  { 1, 2, 4, 5},
                  { 3, 4, 6, 7},
                  { 4, 5, 7, 8}};

  int q[16] = {1,1,1,1,2,2,2,2,2,2,2,2,4,4,4,4};


  EXPECT_EQ(q[id], n);
  EXPECT_TRUE(sameElements(iCs,iCs+n,s[id],s[id]+n));

}
INSTANTIATE_TEST_CASE_P(VertexStar, Quad4VertexStarTest, ::testing::ValuesIn(Quad4VerticesId));



class Quad8NodeStarTest : public testing::TestWithParam<int>
{ // node id
public:

  virtual void SetUp() {
    mesh = Mesh::create(QUADRANGLE8);
    msh_reader.readFileMsh("meshes/simpquad8.msh", mesh);
    vtk_printer.attachMesh(mesh);
    vtk_printer.writeVtk("meshes/out/simpquad8.vtk");
    //
    id = GetParam();
    //
    p = mesh->getNode(id);
    //
    n = (int)(mesh->nodeStar(p, iCs, viCs)-iCs);


  };

  virtual void TearDown()
  {
    delete mesh;
  }

  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], viCs[64];
  int n, id;
  Point *p;
};

int Quad8NodesId[] = {3, 6, 29, 24};

TEST_P(Quad8NodeStarTest, nodeStarTest)
{

  // node star
  int s[30][4];
  s[ 3][0]=6;
  s[24][0]=0; s[24][1]=1; s[24][2]=3; s[24][3]=4;
  s[29][0]=0; s[29][1]=3;
  s[ 6][0]=8;

  int q[30];
  q[ 3] = 1;
  q[24] = 4;
  q[29] = 2;
  q[ 6] = 1;

  EXPECT_EQ(q[id], n);
  EXPECT_TRUE(sameElements(iCs,iCs+n,s[id],s[id]+n));

}
INSTANTIATE_TEST_CASE_P(NodeStar, Quad8NodeStarTest, ::testing::ValuesIn(Quad8NodesId));




class Quad9NodeStarTest : public testing::TestWithParam<int>
{ // node id
public:

  virtual void SetUp() {
    mesh = Mesh::create(QUADRANGLE9);
    msh_reader.readFileMsh("meshes/simpquad9.msh", mesh);
    vtk_printer.attachMesh(mesh);
    vtk_printer.writeVtk("meshes/out/simpquad9.vtk");

    id = GetParam();

    p = mesh->getNode(id);

    n = (int)(mesh->nodeStar(p, iCs, viCs)-iCs);

  };

  virtual void TearDown()
  {
    delete mesh;
  }

  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], viCs[64];
  int n, id;
  Point *p;
};

int Quad9NodesId[] = {0, 27, 4, 8, 30};

TEST_P(Quad9NodeStarTest, nodeStarTest)
{

  // node star
  int s[31][4];
  s[ 0][0]=8;
  s[27][0]=4; s[27][1]=5; s[27][2]=7; s[27][3]=8;
  s[ 4][0]=8; s[ 4][1]=5;
  s[ 8][0]=2;
  s[30][0]=0;

  int q[31];
  q[ 0] = 1;
  q[27] = 4;
  q[ 4] = 2;
  q[ 8] = 1;
  q[30] = 1;

  EXPECT_EQ(q[id], n);
  EXPECT_TRUE(sameElements(iCs,iCs+n,s[id],s[id]+n));

}
INSTANTIATE_TEST_CASE_P(NodeStar, Quad9NodeStarTest, ::testing::ValuesIn(Quad9NodesId));




class Tet4VertexStarTest : public testing::TestWithParam<int> { // vertex id
public:

  virtual void SetUp() {
    mesh = Mesh::create(TETRAHEDRON4);
    mesh->qBuildAdjacency(true);
    msh_reader.readFileMsh("meshes/simptet4.msh", mesh);
    vtk_printer.attachMesh(mesh);
    vtk_printer.writeVtk("meshes/out/simptet4.vtk");

    id = GetParam();

    p = mesh->getNode(id);

    n = (int)(mesh->vertexStar(p, iCs, viCs)-iCs);

  };

  virtual void TearDown()
  {
    delete mesh;
  }

  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], viCs[64];
  int n, id;
  Point *p;
};

int Tet4VerticesId[] = {0,1,2,3,4,5,6,7,8};
//int Tet4NodesId[] = {0};

TEST_P(Tet4VertexStarTest, vertexStarTest) {


  // vertex star
  int s[9][12] = {{ 6, 7, 8, 9,-1,-1,-1,-1,-1,-1,-1,-1},
                  { 0, 5, 6, 9,11,10,-1,-1,-1,-1,-1,-1},
                  { 2, 5,11,-1,-1,-1,-1,-1,-1,-1,-1,-1},
                  { 2, 4, 5, 8, 9,-1,-1,-1,-1,-1,-1,-1},
                  { 0, 1,10,-1,-1,-1,-1,-1,-1,-1,-1,-1},
                  { 1, 2, 3, 4,10,11,-1,-1,-1,-1,-1,-1},
                  { 3, 4, 7, 8,-1,-1,-1,-1,-1,-1,-1,-1},
                  { 0, 1, 3, 6, 7,-1,-1,-1,-1,-1,-1,-1},
                  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11}};

  int q[9] = {4,6,3,5,3,6,4,5,12};

  EXPECT_EQ(q[id], n);
  EXPECT_TRUE(sameElements(iCs,iCs+n,s[id],s[id]+n));

}
INSTANTIATE_TEST_CASE_P(VertexStar, Tet4VertexStarTest, ::testing::ValuesIn(Tet4VerticesId));




class Tet10NodeStarTest : public testing::TestWithParam<int>
{ // node id
public:

  virtual void SetUp() {
    mesh = Mesh::create(TETRAHEDRON10);
    msh_reader.readFileMsh("meshes/simptet10.msh", mesh);
    vtk_printer.attachMesh(mesh);
    vtk_printer.writeVtk("meshes/out/simptet10.vtk");

    id = GetParam();

    p = mesh->getNode(id);

    n = (int)(mesh->nodeStar(p, iCs, viCs)-iCs);

  };

  virtual void TearDown()
  {
    delete mesh;
  }

  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], viCs[64];
  int n, id;
  Point *p;
};

int Tet10NodesId[] = {22,26,30};
//int Tri3NodesId[] = {0,1,2,3,4,5,6,7,8,9,10,11};

TEST_P(Tet10NodeStarTest, nodeStarTest)
{

  // node star
  int s[31][11];
  s[22][0]=5; s[22][1]=1;
  for (int i = 0; i < 12; ++i)
  {
    s[26][i] = i;
  }
  s[30][0]=1;s[30][1]=2;s[30][2]=3;s[30][3]=7;s[30][4]=8;

  int q[31];
  q[22] = 2;
  q[26] = 12;
  q[30] = 5;

  EXPECT_EQ(q[id], n);
  EXPECT_TRUE(sameElements(iCs,iCs+n,s[id],s[id]+n));

}
INSTANTIATE_TEST_CASE_P(NodeStar, Tet10NodeStarTest, ::testing::ValuesIn(Tet10NodesId));






class Hex8VertexStarTest : public testing::TestWithParam<int> { // vertex id
public:

  virtual void SetUp() {
    mesh = Mesh::create(HEXAHEDRON8);
    mesh->qBuildAdjacency(true);
    msh_reader.readFileMsh("meshes/simphex8.msh", mesh);
    vtk_printer.attachMesh(mesh);
    vtk_printer.writeVtk("meshes/out/simphex8.vtk");

    id = GetParam();

    p = mesh->getNode(id);

    n = (int)(mesh->vertexStar(p, iCs, viCs)-iCs);

  };

  virtual void TearDown()
  {
    delete mesh;
  }

  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], viCs[64];
  int n, id;
  Point *p;
};

int Hex8VerticesId[] = {0,25,26};
//int Hex8NodesId[] = {0};

TEST_P(Hex8VertexStarTest, vertexStarTest) {


  // vertex star
  int s[27][8];
  s[0][0] = 1;
  s[25][0] = 0;
  s[25][1] = 2;
  s[25][2] = 4;
  s[25][3] = 6;
  s[26][0] = 0;
  s[26][1] = 2;
  s[26][2] = 4;
  s[26][3] = 6;
  s[26][4] = 1;
  s[26][5] = 3;
  s[26][6] = 5;
  s[26][7] = 7;

  int q[27] = {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,8};

  EXPECT_EQ(q[id], n);
  EXPECT_TRUE(sameElements(iCs,iCs+n,s[id],s[id]+n));
}
INSTANTIATE_TEST_CASE_P(VertexStar, Hex8VertexStarTest, ::testing::ValuesIn(Hex8VerticesId));



class Hex20NodeStarTest : public testing::TestWithParam<int>
{ // node id
public:

  virtual void SetUp() {
    mesh = Mesh::create(HEXAHEDRON20);
    msh_reader.readFileMsh("meshes/simphex20.msh", mesh);
    vtk_printer.attachMesh(mesh);
    vtk_printer.writeVtk("meshes/out/simphex20.vtk");

    id = GetParam();

    p = mesh->getNode(id);

    n = (int)(mesh->nodeStar(p, iCs, viCs)-iCs);

  };

  virtual void TearDown()
  {
    delete mesh;
  }

  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], viCs[64];
  int n, id;
  Point *p;
};

int Hex20NodesId[] = {74,79,58};

TEST_P(Hex20NodeStarTest, nodeStarTest)
{

  // node star
  int s[80][8];
  s[79][0]=2; s[79][1]=3; s[79][2]=6; s[79][3]=7;
  for (int i = 0; i < 8; ++i)
  {
    s[74][i] = i;
  }
  s[58][0]=6;s[58][1]=7;

  int q[80];
  q[79] = 4;
  q[74] = 8;
  q[58] = 2;

  EXPECT_EQ(q[id], n);
  EXPECT_TRUE(sameElements(iCs,iCs+n,s[id],s[id]+n));

}
INSTANTIATE_TEST_CASE_P(NodeStar, Hex20NodeStarTest, ::testing::ValuesIn(Hex20NodesId));



class Hex27NodeStarTest : public testing::TestWithParam<int>
{ // node id
public:

  virtual void SetUp() {
    mesh = Mesh::create(HEXAHEDRON27);
    msh_reader.readFileMsh("meshes/simphex27.msh", mesh);
    vtk_printer.attachMesh(mesh);
    vtk_printer.writeVtk("meshes/out/simphex27.vtk");

    id = GetParam();

    p = mesh->getNode(id);

    n = (int)(mesh->nodeStar(p, iCs, viCs)-iCs);

  };

  virtual void TearDown()
  {
    delete mesh;
  }

  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], viCs[64];
  int n, id;
  Point *p;
};

int Hex27NodesId[] = {98,101,112,105,68};

TEST_P(Hex27NodeStarTest, nodeStarTest)
{

  // node star
  int s[113][8];
  s[112][0]=2; s[112][1]=3;
  for (int i = 0; i < 8; ++i)
  {
    s[98][i] = i;
  }
  for (int i = 0; i < 4; ++i)
  {
    s[101][i] = i;
  }
  s[105][0]=0;
  s[68][0]=6; s[68][1]=7;


  int q[113];
  q[112] = 2;
  q[98] = 8;
  q[101] = 4;
  q[105] = 1;
  q[68] = 2;

  EXPECT_EQ(q[id], n);
  EXPECT_TRUE(sameElements(iCs,iCs+n,s[id],s[id]+n));

}
INSTANTIATE_TEST_CASE_P(NodeStar, Hex27NodeStarTest, ::testing::ValuesIn(Hex27NodesId));



class Tet4EdgeStarTest : public testing::TestWithParam<int> { // vertex id
public:

  virtual void SetUp() {
    mesh = Mesh::create(TETRAHEDRON4);
    msh_reader.readFileMsh("meshes/simptet4.msh", mesh);

    id = GetParam();

    f = mesh->getCorner(id);

    n = (int)(mesh->edgeStar(f, iCs, eiCs)-iCs);

    mesh->getCornerNodesId(f, nds.data());

    if (nds[0]>nds[1])
      std::swap(nds[0], nds[1]);

  };

  virtual void TearDown()
  {
    delete mesh;
  }

  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], eiCs[64];
  int n, id;
  std::tr1::array<int, 2> nds;
  Corner *f;
};

int Tet4EdgesId[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};

std::tr1::array<int, 2> _retArray(int i, int j) {std::tr1::array<int, 2> v;v[0]=i;v[1]=j; return v;}

TEST_P(Tet4EdgeStarTest, vertexStarTest) {


  // nodes of the edges
  std::tr1::array<int, 2> nodes[]= {_retArray(0,1), _retArray(0,3),
                                    _retArray(0,6), _retArray(0,7),
                                    _retArray(0,8), _retArray(1,2),
                                    _retArray(1,3), _retArray(1,4),
                                    _retArray(1,5), _retArray(1,7),
                                    _retArray(1,8), _retArray(2,3),
                                    _retArray(2,5), _retArray(2,8),
                                    _retArray(3,5), _retArray(3,6),
                                    _retArray(3,8), _retArray(4,5),
                                    _retArray(4,7), _retArray(4,8),
                                    _retArray(5,6), _retArray(5,7),
                                    _retArray(5,8), _retArray(6,7),
                                    _retArray(6,8), _retArray(7,8)};

  int id2 = std::distance(nodes, std::find(nodes, nodes+26, nds));

  // incident cells
  int s[26][6] = {{ 9, 6,-1,-1,-1,-1},{ 8, 9,-1,-1,-1,-1},
                  { 7, 8,-1,-1,-1,-1},{ 6, 7,-1,-1,-1,-1},
                  { 6, 7, 8, 9,-1,-1},{ 5,11,-1,-1,-1,-1},
                  { 5, 9,-1,-1,-1,-1},{ 0,10,-1,-1,-1,-1},
                  {10,11,-1,-1,-1,-1},{ 0, 6,-1,-1,-1,-1},
                  { 0, 6,10,11, 5, 9},{ 5, 2,-1,-1,-1,-1},
                  { 2,11,-1,-1,-1,-1},{11, 2, 5,-1,-1,-1},
                  { 2, 4,-1,-1,-1,-1},{ 4, 8,-1,-1,-1,-1},
                  { 2, 4, 8, 9, 5,-1},{ 1,10,-1,-1,-1,-1},
                  { 0, 1,-1,-1,-1,-1},{ 0, 1,10,-1,-1,-1},
                  { 3, 4,-1,-1,-1,-1},{ 1, 3,-1,-1,-1,-1},
                  { 1, 3, 4, 2,11,10},{ 3, 7,-1,-1,-1,-1},
                  { 8, 7, 3, 4,-1,-1},{ 3, 7, 6, 0, 1,-1}};

  // num of incident cells
  int q[26] = {2,2,2,2,4,2,2,2,2,2,6,2,2,3,2,2,5,2,2,3,2,2,6,2,4,5};

  EXPECT_EQ(q[id2], n);
  EXPECT_TRUE(sameElements(iCs,iCs+n,s[id2],s[id2]+n));

}
INSTANTIATE_TEST_CASE_P(EdgeStar, Tet4EdgeStarTest, ::testing::ValuesIn(Tet4EdgesId));





TEST(Tet4ConnectedVtcsTest, ConnectedVtcs)
{
  MeshIoMsh msh_reader;
  Mesh *mesh;
  int iNs[256];
  int n, id;
  Point *p;

  mesh = Mesh::create(TETRAHEDRON4);
  msh_reader.readFileMsh("meshes/simptet4.msh", mesh);

  id = 8;

  p = mesh->getNode(id);

  n = (int)(mesh->connectedVtcs(p, iNs)-iNs);  
    
  // connected nodes
  int s[9][8];

  s[8][ 0]=0;
  s[8][ 1]=1;
  s[8][ 2]=2;
  s[8][ 3]=3;
  s[8][ 4]=4;
  s[8][ 5]=5;
  s[8][ 6]=6;
  s[8][ 7]=7;

  int q[9];
  q[8] = 8;

  if (!sameElements(iNs,iNs+n,s[id],s[id]+n))
  {
    //std::sort(iNs,iNs+n);
    //std::sort(s[id], s[id]+n);
    //printf("  actual:");
    //for (int i = 0; i < n; ++i)
      //printf("%d ", iNs[i]);
    //printf("\nexpected:");
    //for (int i = 0; i < n; ++i)
      //printf("%d ", s[id][i]);
    //printf("\n");
  }

  EXPECT_EQ(q[id], n);
  EXPECT_TRUE(sameElements(iNs,iNs+n,s[id],s[id]+n));

  delete mesh;  
}



class Tet10ConnectedNodesTest : public testing::TestWithParam<int>
{ // node id
public:

  virtual void SetUp() {
    mesh = Mesh::create(TETRAHEDRON10);
    msh_reader.readFileMsh("meshes/simptet10.msh", mesh);

    id = GetParam();

    p = mesh->getNode(id);

    n = (int)(mesh->connectedNodes(p, iNs)-iNs);

  };

  virtual void TearDown()
  {
    delete mesh;
  }

  MeshIoMsh msh_reader;
  Mesh *mesh;
  int iNs[256];
  int n, id;
  Point *p;
};

//    int Tet10NodesId[] = {22,26,30};

TEST_P(Tet10ConnectedNodesTest, ConnectedNodesTest)
{

  // connected nodes
  int s[31][34];

  s[22][ 0]=2;
  s[22][ 1]=10;
  s[22][ 2]=17;
  s[22][ 3]=30;
  s[22][ 4]=5;
  s[22][ 5]=26;
  s[22][ 6]=3;
  s[22][ 7]=29;
  s[22][ 8]=31;
  s[22][ 9]=18;
  s[22][10]=13;
  s[22][11]=34;
  s[22][12]=6;

  int k = 0;
  for (int i = 0; i < 35; ++i)
  {
    if (i==26) continue;
    s[26][k++] = i;
  }

  s[30][ 0]=2;
  s[30][ 1]=3;
  s[30][ 2]=5;
  s[30][ 3]=26;
  s[30][ 4]=17;
  s[30][ 5]=10;
  s[30][ 6]=22;
  s[30][ 7]=29;
  s[30][ 8]=31;
  s[30][ 9]=12;
  s[30][10]=32;
  s[30][11]=21;
  s[30][12]=4;
  s[30][13]=1;
  s[30][14]=16;
  s[30][15]=9;
  s[30][16]=33;
  s[30][17]=8;
  s[30][18]=20;
  s[30][19]=0;
  s[30][20]=11;
  s[30][21]=27;

  int q[31];
  q[22] = 13;
  q[26] = 34;
  q[30] = 22;

  if (!sameElements(iNs,iNs+n,s[id],s[id]+n))
  {
    //std::sort(iNs,iNs+n);
    //std::sort(s[id], s[id]+n);
    //printf("  actual:");
    //for (int i = 0; i < n; ++i)
      //printf("%d ", iNs[i]);
    //printf("\nexpected:");
    //for (int i = 0; i < n; ++i)
      //printf("%d ", s[id][i]);
    //printf("\n");
  }

  EXPECT_EQ(q[id], n);
  EXPECT_TRUE(sameElements(iNs,iNs+n,s[id],s[id]+n));



}
INSTANTIATE_TEST_CASE_P(ConnectedNodes, Tet10ConnectedNodesTest, ::testing::ValuesIn(Tet10NodesId));



// iterators
TEST(PointIteratorsTest, TetIteratorsTest)
{
  MeshIoMsh msh_reader;
  Mesh *mesh;
  mesh = Mesh::create(TETRAHEDRON4);
  msh_reader.readFileMsh("meshes/uni_tet.msh", mesh);
  
  for (int i = 0; i < mesh->numNodes(); ++i)
  {
    mesh->getNode(i)->setTag(0);
  }
  
  // ====================== traversing the mesh ====================================
  
  {
    point_iterator point = mesh->pointBegin();
    point_iterator point_end = mesh->pointEnd();
    
    for (; point != point_end; ++point)
      point->setTag(point->getTag()+1);
    // check
    for (int i = 0; i < mesh->numNodes(); ++i)
    {
      EXPECT_EQ( 1, mesh->getNode(i)->getTag());
      mesh->getNode(i)->setTag(0); // reseting
    }
  }
  
  
};







TEST(SingleCellTestTri3, SingleSimpMesh)
{

  MeshIoMsh msh_reader;
  Mesh *mesh;  

  ECellType cell_t     = TRIANGLE3;
  const char* mesh_in  = "meshes/simptri3.msh";
  
  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);

  bool singles[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    
  Point const* p;
  for (int i = 0; i < 12; ++i)
  {
    p = mesh->getNode(i);
    EXPECT_EQ(singles[i], mesh->inSingleCell(p));
  }

  Facet const* f;
  int n_sing_facets = 0;
  for (int i = 0; i < mesh->numFacets(); ++i)
  {
    f = mesh->getFacet(i);
    if (mesh->inSingleCell(f))
      ++n_sing_facets;
  }
  EXPECT_EQ(8, n_sing_facets);
  delete mesh;  
  
}
             
                          
TEST(SingleCellTestTri3, SingleSingMesh)
{

  MeshIoMsh msh_reader;
  Mesh *mesh;  

  ECellType cell_t     = TRIANGLE3;
  const char* mesh_in  = "meshes/sing_tri3.msh";
  
  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);

  bool singles[] = {1,0,0,0,1,0,0};
    
  Point const*p;
  for (int i = 0; i < 7; ++i)
  {
    p = mesh->getNode(i);
    EXPECT_EQ(singles[i], mesh->inSingleCell(p));
  }

  //std::cout << mesh->getCell(0)->getIncidCell(0) << std::endl;
  //std::cout << mesh->getCell(0)->getIncidCell(1) << std::endl;
  //std::cout << mesh->getCell(0)->getIncidCell(2) << std::endl;

  delete mesh;  
  
}

TEST(SingleCellTestTri6, SingleSimpMesh)
{

  MeshIoMsh msh_reader;
  Mesh *mesh;  

  ECellType cell_t     = TRIANGLE6;
  const char* mesh_in  = "meshes/simptri6.msh";
  
  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);

  bool sing_pts[] = {0,0,0,0,0,1,1,0,1,1,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
  Point const*p;
  for (int i = 0; i < 37; ++i)
  {
    p = mesh->getNode(i);
    EXPECT_EQ(sing_pts[i], mesh->inSingleCell(p));
  }

  delete mesh;  
  
}
             
                          
TEST(SingleCellTestTri6, SingleSingMesh)
{

  MeshIoMsh msh_reader;
  Mesh *mesh;  
  
  ECellType cell_t     = TRIANGLE6;
  const char* mesh_in  = "meshes/sing_tri6.msh";
  
  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);
  
  bool sing_pts[] = {1,0,0,0,1,0,1,1,1,1,0,1,1,1,0,0,0,0};

  Point const*p;
  for (int i = 0; i < 18; ++i)
  {
    p = mesh->getNode(i);
    EXPECT_EQ(sing_pts[i], mesh->inSingleCell(p));
  }
  
  //std::cout << mesh->getCell(0)->getIncidCell(0) << std::endl;
  //std::cout << mesh->getCell(0)->getIncidCell(1) << std::endl;
  //std::cout << mesh->getCell(0)->getIncidCell(2) << std::endl;
  
  delete mesh;  
  
}



TEST(SingleCellTestTet4, SingleSimpMesh)
{

  MeshIoMsh msh_reader;
  Mesh *mesh;  

  ECellType cell_t     = TETRAHEDRON4;
  const char* mesh_in  = "meshes/simptet4.msh";
  
  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);

  bool singles[] = {0,0,0,0,0,0,0,0,0};
    
  Point const* p;
  for (int i = 0; i < 8; ++i)
  {
    p = mesh->getNode(i);
    EXPECT_EQ(singles[i], mesh->inSingleCell(p));
  }

  // check facets
  Facet const* f;
  int n_sing_facets = 0;
  for (int i = 0; i < mesh->numFacets(); ++i)
  {
    f = mesh->getFacet(i);
    if (mesh->inSingleCell(f))
      ++n_sing_facets;
  }
  EXPECT_EQ(12, n_sing_facets);
  
  
  // check corners
  Corner const* r;
  int n_sing_corners = 0;
  for (int i = 0; i < mesh->numCorners(); ++i)
  {
    r = mesh->getCorner(i);
    if (mesh->inSingleCell(r))
      ++n_sing_corners;
  }
  EXPECT_EQ(0, n_sing_corners);
  
  
  delete mesh;
}


// PROVAVELMENTE ESTA MALHA É INVÁLIDA
TEST(SingleCellTestTet4, SingleSingMesh)
{

  MeshIoMsh msh_reader;
  Mesh *mesh;  

  ECellType cell_t     = TETRAHEDRON4;
  const char* mesh_in  = "meshes/sing_tet4.msh";
  
  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);

  bool singles[] = {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0};
    
  Point const* p;
  for (int i = 0; i < 15; ++i)
  {
    p = mesh->getNode(i);
    EXPECT_EQ(singles[i], mesh->inSingleCell(p));
  }
  
  // check facets
  Facet const* f;
  int n_sing_facets = 0;
  for (int i = 0; i < mesh->numFacets(); ++i)
  {
    f = mesh->getFacet(i);
    if (mesh->inSingleCell(f))
      ++n_sing_facets;
  }
  EXPECT_EQ(24, n_sing_facets);
  
  // check corners
  Corner const* r;
  int n_sing_corners = 0;
  for (int i = 0; i < mesh->numCorners(); ++i)
  {
    r = mesh->getCorner(i);
    if (mesh->inSingleCell(r))
      ++n_sing_corners;
  }
  EXPECT_EQ(3, n_sing_corners);
  
  
  delete mesh;
}





TEST(SingleCellTestTet10, SingleSingMesh)
{

  MeshIoMsh msh_reader;
  Mesh *mesh;  

  ECellType cell_t     = TETRAHEDRON10;
  const char* mesh_in  = "meshes/sing_tet10.msh";
  
  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);

    
  Point const* p;
  for (int i = 0; i < mesh->numNodes(); ++i)
  {
    p = mesh->getNode(i);
    EXPECT_EQ( i==24 || i==26, mesh->inSingleCell(p)) << " i = " << i;
  }
  
  // check facets
  Facet const* f;
  int n_sing_facets = 0;
  for (int i = 0; i < mesh->numFacets(); ++i)
  {
    f = mesh->getFacet(i);
    if (mesh->inSingleCell(f))
    {
      ++n_sing_facets;
    }
  }
  EXPECT_EQ(28, n_sing_facets);
  
  // check corners
  Corner const* r;
  int n_sing_corners = 0;
  for (int i = 0; i < mesh->numCorners(); ++i)
  {
    r = mesh->getCorner(i);
    if (mesh->inSingleCell(r))
      ++n_sing_corners;
  }
  EXPECT_EQ(2, n_sing_corners);
  
  
  
  delete mesh;
}



TEST(IncidentFacetsTestTri6, IncidentFacets)
{

  MeshIoMsh msh_reader;
  Mesh *mesh;  

  ECellType cell_t     = TRIANGLE6;
  const char* mesh_in  = "meshes/simptri6.msh";
  
  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);
  
  int iFs[32], viFs[32];
  int *iFs_end;
  int fnodes[32];
  
  // num incident facets
  int const nif[] = {3,3,3,3,4,0,0,4,0,0,4,0,0,4,0,0,5,6,5,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  for (int i = 0; i < 37; ++i)
  {
    if (!mesh->isVertex(mesh->getNode(i)))
      continue;
    iFs_end = mesh->incidentFacets(mesh->getNode(i), iFs, viFs);
    
    EXPECT_EQ(nif[i], iFs_end - iFs) << "at node i = " << i << std::endl;
    
    for (int *it = iFs; it != iFs_end; ++it)
    {
      EXPECT_TRUE(*it>=0 && *it < 25) << "at node i = " << i << std::endl;
     
      for (int *at = iFs; at != iFs_end ; ++at)
      {
        if (at != it)
          EXPECT_TRUE(*at != *it) << "at node i = " << i << std::endl;
      }
      
      mesh->getFacetNodesId(*it,fnodes);
      
      EXPECT_TRUE(fnodes[viFs[it-iFs]] == i) << "at node i = " << i << std::endl;
    }
  }

  

  delete mesh;  
  
}
             
  



