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

TEST(MeshTest, SizeOfElementsTest)
{
  std::cout << "sizeof(Edge2)         " << sizeof(Edge2)          << std::endl;
  std::cout << "sizeof(Edge3)         " << sizeof(Edge3)          << std::endl;
  std::cout << "sizeof(Triangle3)     " << sizeof(Triangle3)     << std::endl;
  std::cout << "sizeof(Triangle6)     " << sizeof(Triangle6)     << std::endl;
  std::cout << "sizeof(Tetrahedron4)  " << sizeof(Tetrahedron4)  << std::endl;
  std::cout << "sizeof(Tetrahedron10) " << sizeof(Tetrahedron10) << std::endl;
  std::cout << "sizeof(Hexahedron8)   " << sizeof(Hexahedron8)   << std::endl;
  std::cout << "sizeof(Hexahedron20)  " << sizeof(Hexahedron20)  << std::endl;
  std::cout << "sizeof(Hexahedron27)  " << sizeof(Hexahedron27)  << std::endl;
                
  std::cout << "sizeof(Facet)         " << sizeof(Facet)          << std::endl;
  std::cout << "sizeof(Corner)        " << sizeof(Corner)         << std::endl;
  std::cout << "sizeof(Point1d)       " << sizeof(Point1d)        << std::endl;
  std::cout << "sizeof(Point2d)       " << sizeof(Point2d)        << std::endl;
  std::cout << "sizeof(Point3d)       " << sizeof(Point3d)        << std::endl;
  
  EXPECT_TRUE(1);
}

TEST(IOTest, identifiesMshMeshTypeTest)
{
  MeshIoMsh msh_reader;
  int sd;
  ECellType ct;

  ct = msh_reader.identifiesMshMeshType("meshes/simpedge2.msh",sd);
  EXPECT_EQ(EDGE2, ct);
  EXPECT_EQ(sd, 2);

  ct = msh_reader.identifiesMshMeshType("meshes/simpedge3.msh",sd);
  EXPECT_EQ(EDGE3, ct);
  EXPECT_EQ(sd, 2);

  ct = msh_reader.identifiesMshMeshType("meshes/simptri3.msh",sd);
  EXPECT_EQ(TRIANGLE3, ct);
  EXPECT_EQ(sd, 2);

  ct = msh_reader.identifiesMshMeshType("meshes/simptri6.msh",sd);
  EXPECT_EQ(TRIANGLE6, ct);
  EXPECT_EQ(sd, 2);

  ct = msh_reader.identifiesMshMeshType("meshes/simpquad4.msh",sd);
  EXPECT_EQ(QUADRANGLE4, ct);
  EXPECT_EQ(sd, 2);

  ct = msh_reader.identifiesMshMeshType("meshes/simpquad8.msh",sd);
  EXPECT_EQ(QUADRANGLE8, ct);
  EXPECT_EQ(sd, 2);

  ct = msh_reader.identifiesMshMeshType("meshes/simpquad9.msh",sd);
  EXPECT_EQ(QUADRANGLE9, ct);
  EXPECT_EQ(sd, 2);

  ct = msh_reader.identifiesMshMeshType("meshes/simptet4.msh",sd);
  EXPECT_EQ(TETRAHEDRON4, ct);
  EXPECT_EQ(sd, 3);

  ct = msh_reader.identifiesMshMeshType("meshes/simptet4.msh",sd);
  EXPECT_EQ(TETRAHEDRON4, ct);
  EXPECT_EQ(sd, 3);

  ct = msh_reader.identifiesMshMeshType("meshes/simptet10.msh",sd);
  EXPECT_EQ(TETRAHEDRON10, ct);
  EXPECT_EQ(sd, 3);

  ct = msh_reader.identifiesMshMeshType("meshes/simphex8.msh",sd);
  EXPECT_EQ(HEXAHEDRON8, ct);
  EXPECT_EQ(sd, 3);

  ct = msh_reader.identifiesMshMeshType("meshes/simphex20.msh",sd);
  EXPECT_EQ(HEXAHEDRON20, ct);
  EXPECT_EQ(sd, 3);

  ct = msh_reader.identifiesMshMeshType("meshes/simphex27.msh",sd);
  EXPECT_EQ(HEXAHEDRON27, ct);
  EXPECT_EQ(sd, 3);

}

TEST(ReadVtktest, WithEdge2)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int dim = 2;

  mesh = Mesh::create(EDGE2,dim);
  msh_reader.readFileMsh("meshes/simpedge2.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.writeVtk("meshes/outtest/simpedge2.vtk");

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

TEST(ReadVtkTest, withEdge3)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int dim = 2;

  mesh = Mesh::create(EDGE3,dim);
  msh_reader.readFileMsh("meshes/simpedge3.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.writeVtk("meshes/outtest/simpedge3.vtk");

  delete mesh;
}

TEST(VertexStarTest, WithTri3)
{

  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], viCs[64];
  int n, id;
  Point *p;

  mesh = Mesh::create(TRIANGLE3);
  msh_reader.readFileMsh("meshes/simptri3.msh", mesh);
  vtk_printer.attachMesh(mesh);
  //vtk_printer.writeVtk("meshes/outtest/simptri3.vtk");

  for (int k = 0; k < 12; ++k)
  {
    id = k;

    p = mesh->getNode(id);

    n = (int)(mesh->vertexStar(p, iCs, viCs)-iCs);

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

  delete mesh;
}

TEST(VertexStarTest, WithQuad4)
{

  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], viCs[64];
  int n, id;
  Point *p;

  mesh = Mesh::create(QUADRANGLE4);
  msh_reader.readFileMsh("meshes/simpquad4.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.writeVtk("meshes/outtest/simpquad4.vtk");

  for (int k = 0; k < 16; ++k)
  {
    id = k;

    p = mesh->getNode(id);

    n = (int)(mesh->vertexStar(p, iCs, viCs)-iCs);

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

  delete mesh;
}

TEST(VertexStarTest, WithTet4)
{

  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], viCs[64];
  int n, id;
  Point *p;

  mesh = Mesh::create(TETRAHEDRON4);
  mesh->qBuildAdjacency(true);
  msh_reader.readFileMsh("meshes/simptet4.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.writeVtk("meshes/outtest/simptet4.vtk");

  int Tet4VerticesId[] = {0,1,2,3,4,5,6,7,8};

  for (int *it = Tet4VerticesId; it != Tet4VerticesId+sizeof(Tet4VerticesId)/sizeof(int); ++it)
  {
    id = *it;

    p = mesh->getNode(id);

    n = (int)(mesh->vertexStar(p, iCs, viCs)-iCs);

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
  delete mesh;
}

std::tr1::array<int, 2> _retArray(int i, int j) {std::tr1::array<int, 2> v;v[0]=i;v[1]=j; return v;}

TEST(VertexStarTest, WithTet4Case2)
{


  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], eiCs[64];
  int n, id;
  std::tr1::array<int, 2> nds;
  Corner *f;

  mesh = Mesh::create(TETRAHEDRON4);
  msh_reader.readFileMsh("meshes/simptet4.msh", mesh);

  int Tet4EdgesId[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};

  for (int *it=Tet4EdgesId; it!=Tet4EdgesId+sizeof(Tet4EdgesId)/sizeof(int); ++it)
  {
    id = *it;

    f = mesh->getCorner(id);

    n = (int)(mesh->edgeStar(f, iCs, eiCs)-iCs);

    mesh->getCornerNodesId(f, nds.data());

    if (nds[0]>nds[1])
      std::swap(nds[0], nds[1]);

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
  delete mesh;
}

TEST(VertexStarTest, WithHex8)
{

  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], viCs[64];
  int n, id;
  Point *p;

  mesh = Mesh::create(HEXAHEDRON8);
  mesh->qBuildAdjacency(true);
  msh_reader.readFileMsh("meshes/simphex8.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.writeVtk("meshes/outtest/simphex8.vtk");

  int Hex8VerticesId[] = {0,25,26};

  for (int *it=Hex8VerticesId; it!=Hex8VerticesId+sizeof(Hex8VerticesId)/sizeof(int);++it)
  {
    id = *it;

    p = mesh->getNode(id);

    n = (int)(mesh->vertexStar(p, iCs, viCs)-iCs);

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

  delete mesh;
}

TEST(NodeStarTest, WithTri6)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], viCs[64];
  int n, id;
  Point *p;

  mesh = Mesh::create(TRIANGLE6);
  msh_reader.readFileMsh("meshes/simptri6.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.writeVtk("meshes/outtest/simptri6.vtk");

  for (int k = 0; k < 24; ++k)
  {
    if (!(k==2 || k == 5 || k ==23 || k == 16))
      continue;
    id = k;

    p = mesh->getNode(id);

    n = (int)(mesh->nodeStar(p, iCs, viCs)-iCs);

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

  delete mesh;
}

TEST(NodeStarTest, WithQuad8)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], viCs[64];
  int n, id;
  Point *p;

  mesh = Mesh::create(QUADRANGLE8);
  msh_reader.readFileMsh("meshes/simpquad8.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.writeVtk("meshes/outtest/simpquad8.vtk");
  //
  int Quad8NodesId[] = {3, 6, 29, 24};
  int *it = Quad8NodesId;

  for (; it != Quad8NodesId+sizeof(Quad8NodesId)/sizeof(int); ++it)
  {
    id = *it;
    //
    p = mesh->getNode(id);
    //
    n = (int)(mesh->nodeStar(p, iCs, viCs)-iCs);

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

  delete mesh;
}

TEST(NodeStarTest, WithQuad9)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], viCs[64];
  int n, id;
  Point *p;

  mesh = Mesh::create(QUADRANGLE9);
  msh_reader.readFileMsh("meshes/simpquad9.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.writeVtk("meshes/outtest/simpquad9.vtk");

  int Quad9NodesId[] = {0, 27, 4, 8, 30};
  int *it = Quad9NodesId;

  for (; it != Quad9NodesId+sizeof(Quad9NodesId)/sizeof(int); ++it)
  {
    id = *it;

    p = mesh->getNode(id);

    n = (int)(mesh->nodeStar(p, iCs, viCs)-iCs);

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
  delete mesh;
}

TEST(NodeStarTest, WithTet10)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], viCs[64];
  int n, id;
  Point *p;

  mesh = Mesh::create(TETRAHEDRON10);
  msh_reader.readFileMsh("meshes/simptet10.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.writeVtk("meshes/outtest/simptet10.vtk");

  int Tet10NodesId[] = {22,26,30};

  for (int *it = Tet10NodesId; it != Tet10NodesId+sizeof(Tet10NodesId)/sizeof(int); ++it)
  {

    id = *it;

    p = mesh->getNode(id);

    n = (int)(mesh->nodeStar(p, iCs, viCs)-iCs);

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

  delete mesh;
}

TEST(NodeStarTest, WithHex20)
{

  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], viCs[64];
  int n, id;
  Point *p;

  mesh = Mesh::create(HEXAHEDRON20);
  msh_reader.readFileMsh("meshes/simphex20.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.writeVtk("meshes/outtest/simphex20.vtk");

  int Hex20NodesId[] = {74,79,58};

  for (int *it=Hex20NodesId; it!= Hex20NodesId+sizeof(Hex20NodesId)/sizeof(int);++it)
  {
    id = *it;

    p = mesh->getNode(id);

    n = (int)(mesh->nodeStar(p, iCs, viCs)-iCs);

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
  delete mesh;
}

TEST(NodeStarTest, WithHex27)
{

  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int iCs[64], viCs[64];
  int n, id;
  Point *p;

  mesh = Mesh::create(HEXAHEDRON27);
  msh_reader.readFileMsh("meshes/simphex27.msh", mesh);
  vtk_printer.attachMesh(mesh);
  vtk_printer.writeVtk("meshes/outtest/simphex27.vtk");

  int Hex27NodesId[] = {98,101,112,105,68};

  for (int *it=Hex27NodesId; it!= Hex27NodesId+sizeof(Hex27NodesId)/sizeof(int);++it)
  {
    id = *it;

    p = mesh->getNode(id);

    n = (int)(mesh->nodeStar(p, iCs, viCs)-iCs);

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
  delete mesh;
}

TEST(ConnectedNodesTest, WithTet10)
{
  MeshIoMsh msh_reader;
  Mesh *mesh;
  int iNs[256];
  int n, id;
  Point *p;

  mesh = Mesh::create(TETRAHEDRON10);
  msh_reader.readFileMsh("meshes/simptet10.msh", mesh);

  int Tet10NodesId[] = {22,26,30};

  for (int *it = Tet10NodesId; it!=Tet10NodesId+sizeof(Tet10NodesId)/sizeof(int);++it)
  {
    id = *it;

    p = mesh->getNode(id);

    n = (int)(mesh->connectedNodes(p, iNs)-iNs);

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
  delete mesh;
}

TEST(ConnectedVtcs, WithTet4)
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

// iterators
TEST(TetIteratorsTest, PointIteratorsTest)
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

  delete mesh;
};

TEST(SingleCell, WithTri3Case1)
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

TEST(SingleCell, WithTri3Case2)
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

TEST(SingleCell, WithTri6Case1)
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

TEST(SingleCell, WithTri6Case2)
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

TEST(SingleCell, WithTet4Case1)
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

TEST(SingleCell, WithTet4Case2)
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
  EXPECT_EQ(26, n_sing_facets);

  // check corners
  Corner const* r;
  int n_sing_corners = 0;
  for (int i = 0; i < mesh->numCorners(); ++i)
  {
    r = mesh->getCorner(i);
    if (mesh->inSingleCell(r))
      ++n_sing_corners;
  }
  EXPECT_EQ(14, n_sing_corners);


  delete mesh;
}

TEST(SingleCell, WithTet10)
{

  MeshIoMsh msh_reader;
  Mesh *mesh;

  ECellType cell_t     = TETRAHEDRON10;
  const char* mesh_in  = "meshes/sing_tet10.msh";

  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);

  int iCs[512];
  int viCs[512];
  int *cit;

  Point const* p;
  for (int i = 0; i < mesh->numNodes(); ++i)
  {
    p = mesh->getNode(i);
    cit = mesh->nodeStar(p, iCs, viCs);
    EXPECT_EQ( (cit-iCs)==1, mesh->inSingleCell(p)) << " i = " << i;
    //EXPECT_EQ( i==24 || i==26, mesh->inSingleCell(p)) << " i = " << i;
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
  EXPECT_EQ(26, n_sing_facets);

  // check corners
  Corner const* r;
  int n_sing_corners = 0;
  for (int i = 0; i < mesh->numCorners(); ++i)
  {
    r = mesh->getCorner(i);
    if (mesh->inSingleCell(r))
      ++n_sing_corners;
  }
  EXPECT_EQ(14, n_sing_corners);



  delete mesh;
}

TEST(IncidentFacets, WithTri6)
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

TEST(AuxSetConnectedComponentIdTest, WithTri3)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  SMesh<Triangle3,3> *mesh;
  int dim = 3;

  mesh = (SMesh<Triangle3,3>*)Mesh::create(TRIANGLE3,dim);
  msh_reader.readFileMsh("meshes/singular_tri3a.msh", mesh);
  vtk_printer.attachMesh(mesh);
  //vtk_printer.writeVtk("meshes/outtest/simpedge2.vtk");

  mesh->_setConnectedComponentsId(mesh->getCell(0), 3);
  mesh->_setConnectedComponentsId(mesh->getCell(4), 6);
  mesh->_setConnectedComponentsId(mesh->getCell(8), 10);
  mesh->_setConnectedComponentsId(mesh->getCell(6), 99);
  mesh->_setConnectedComponentsId(mesh->getCell(12), 13);

  EXPECT_EQ(3 , mesh->getCell( 0)->getConnectedComponentId());

  EXPECT_EQ(6 , mesh->getCell( 1)->getConnectedComponentId());
  EXPECT_EQ(6 , mesh->getCell( 2)->getConnectedComponentId());
  EXPECT_EQ(6 , mesh->getCell( 3)->getConnectedComponentId());
  EXPECT_EQ(6 , mesh->getCell( 4)->getConnectedComponentId());
  EXPECT_EQ(6 , mesh->getCell( 5)->getConnectedComponentId());

  EXPECT_EQ(10, mesh->getCell( 7)->getConnectedComponentId());
  EXPECT_EQ(10, mesh->getCell( 8)->getConnectedComponentId());
  EXPECT_EQ(10, mesh->getCell( 9)->getConnectedComponentId());
  EXPECT_EQ(10, mesh->getCell(10)->getConnectedComponentId());
  EXPECT_EQ(10, mesh->getCell(11)->getConnectedComponentId());

  EXPECT_EQ(99, mesh->getCell( 6)->getConnectedComponentId());

  EXPECT_EQ(13, mesh->getCell(12)->getConnectedComponentId());

  delete mesh;
}

TEST(SetConnectedComponentIdTest, WithTri3)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int dim = 3;

  mesh = Mesh::create(TRIANGLE3,dim);
  msh_reader.readFileMsh("meshes/singular_tri3a.msh", mesh);
  vtk_printer.attachMesh(mesh);
  //vtk_printer.writeVtk("meshes/outtest/simpedge2.vtk");

  //mesh->setUpConnectedComponentsId();

  int n_connected_components = mesh->numConnectedComponents();
  cell_iterator cell = mesh->cellBegin();
  cell_iterator cell_end = mesh->cellEnd();

  for (; cell != cell_end; ++cell)
  {
    EXPECT_TRUE(cell->getConnectedComponentId() > 0);
  }

  EXPECT_EQ(5, n_connected_components);

  int c_ini[] = {-1,-1,-1,-1,-1};
  int cc_id[] = {-1,-1,-1,-1,-1};

  mesh->getConnectedComponentsPicks(cc_id, c_ini);

  for (int i = 0; i < n_connected_components; ++i)
  {
    EXPECT_TRUE(cc_id[i] >= 0);
    EXPECT_TRUE(c_ini[i] >= 0);
  }

  for (int i = 0; i < n_connected_components; ++i)
    for (int j = i+1; j < n_connected_components; ++j)
      EXPECT_TRUE(mesh->getCell(c_ini[i])->getConnectedComponentId() != mesh->getCell(c_ini[j])->getConnectedComponentId());


  delete mesh;
}

TEST(NextBoundaryFacetTest, WithTri3)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
  int dim = 3;

  mesh = Mesh::create(TRIANGLE3,dim);
  msh_reader.readFileMsh("meshes/singular_tri3a.msh", mesh);
  vtk_printer.attachMesh(mesh);
  //vtk_printer.writeVtk("meshes/outtest/simpedge2.vtk");

  Facet *f0, *fit;

  fit = f0 = mesh->getFacet( mesh->getCell(0)->getFacetId(0) );
  fit = mesh->nextBoundaryFacet(fit);
  EXPECT_EQ(0, fit->getIncidCell());
  EXPECT_EQ(1, fit->getPosition());
  fit = mesh->nextBoundaryFacet(fit);
  EXPECT_EQ(0, fit->getIncidCell());
  EXPECT_EQ(2, fit->getPosition());
  fit = mesh->nextBoundaryFacet(fit);
  EXPECT_EQ(0, fit->getIncidCell());
  EXPECT_EQ(0, fit->getPosition());
  EXPECT_EQ(f0, fit);
  
  fit = f0 = mesh->getFacet( mesh->getCell(1)->getFacetId(0) );
  fit = mesh->nextBoundaryFacet(fit);
  EXPECT_EQ(4, fit->getIncidCell());
  EXPECT_EQ(2, fit->getPosition());
  fit = mesh->nextBoundaryFacet(fit);
  EXPECT_EQ(3, fit->getIncidCell());
  EXPECT_EQ(0, fit->getPosition());
  fit = mesh->nextBoundaryFacet(fit);
  EXPECT_EQ(5, fit->getIncidCell());
  EXPECT_EQ(2, fit->getPosition());
  fit = mesh->nextBoundaryFacet(fit);
  EXPECT_EQ(5, fit->getIncidCell());
  EXPECT_EQ(0, fit->getPosition());
  fit = mesh->nextBoundaryFacet(fit);
  EXPECT_EQ(1, fit->getIncidCell());
  EXPECT_EQ(0, fit->getPosition());
  EXPECT_EQ(f0, fit);


  delete mesh;
}

TEST(SetBoundaryComponentIdTest, WithTri3)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  SMesh<Triangle3,3> *mesh;
  int dim = 3;

  mesh = (SMesh<Triangle3,3>*)Mesh::create(TRIANGLE3,dim);
  msh_reader.readFileMsh("meshes/singular_tri3a.msh", mesh);
  //vtk_printer.attachMesh(mesh);
  //vtk_printer.writeVtk("meshes/outtest/simpedge2.vtk");

  int const n_bound_comps = mesh->numBoundaryComponents();
  
  EXPECT_EQ(5, n_bound_comps);

  int f_inis[] = {-1,-1,-1,-1,-1};
  int bc_ids[] = {-1,-1,-1,-1,-1};

  mesh->getBoundaryComponentsPicks(bc_ids, f_inis);

  for (int i = 0; i < n_bound_comps; ++i)
  {
    Facet *f = mesh->nextBoundaryFacet( mesh->getFacet(f_inis[i]) );
    
    do
    {
      EXPECT_EQ(bc_ids[i], f->getBoundaryComponentId());
      f = mesh->nextBoundaryFacet(f);
    }
    while (f != mesh->getFacet(f_inis[i]));
    
  }
  

  delete mesh;
}

TEST(PushIncidCell2Point, WithSingularVertex)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  SMesh<Triangle3,3> *mesh;
  int dim = 3;

  mesh = (SMesh<Triangle3,3>*)Mesh::create(TRIANGLE3,dim);
  
  mesh->qBuildAdjacency(false);
  
  msh_reader.readFileMsh("meshes/singular_tri3a.msh", mesh);
  //vtk_printer.attachMesh(mesh);
  //vtk_printer.writeVtk("meshes/outtest/simpedge2.vtk");

  Point* p = mesh->getNode(0);

  EXPECT_EQ(0, p->numConnectedComps());
  
  EXPECT_EQ(-1, p->getIncidCell());

  // initial setting
  p->setIncidCell(4);
  p->setPosition(0);
  EXPECT_EQ(4, p->getIncidCell());
  EXPECT_EQ(0, p->getPosition());
  
  mesh->getCell(4)->setConnectedComponentId(1);
  mesh->getCell(3)->setConnectedComponentId(1);
  mesh->getCell(9)->setConnectedComponentId(2);
  mesh->getCell(10)->setConnectedComponentId(2);
  mesh->getCell(12)->setConnectedComponentId(3);

  mesh->pushIncidCell2Point(p, 4,0);
  mesh->pushIncidCell2Point(p, 3,0);
  mesh->pushIncidCell2Point(p, 9,0);
  mesh->pushIncidCell2Point(p, 10,0);
  
  p->replacesIncidCell(12,3,1); // do nothing

  EXPECT_EQ(2, p->numConnectedComps());
  
  p->replacesIncidCell(3,-1,-1);

  EXPECT_EQ(1, p->numConnectedComps());

  mesh->pushIncidCell2Point(p,10,0);
  
  EXPECT_EQ(1, p->numConnectedComps());
  
  mesh->pushIncidCell2Point(p,3,0);
  
  EXPECT_EQ(2, p->numConnectedComps());
  
  delete mesh;
}

TEST(NodeSingularityTest, WithTri3)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  //std::tr1::shared_ptr< SMesh<Triangle3,3> > mesh;
  Mesh* mesh;
  int dim = 2;

  //mesh.reset( (SMesh<Triangle3,3>*)Mesh::create(TRIANGLE3,dim) );
  mesh =Mesh::create(TRIANGLE3,dim) ;
  
  //mesh->qBuildAdjacency(true);
  
  //msh_reader.readFileMsh("meshes/singular_tri3a.msh", mesh.get());
  msh_reader.readFileMsh("meshes/singular_tri3a.msh", mesh);
  //vtk_printer.attachMesh(mesh);
  //vtk_printer.writeVtk("meshes/outtest/simpedge2.vtk");

  bool is_sing[] = {1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0};

  for (int i = 0; i < static_cast<int>( sizeof(is_sing)/sizeof(bool) ); ++i)
  {
    EXPECT_TRUE(mesh->getNode(i)->isSingular()==is_sing[i])  << "i= " << i << std::endl;
  }
  
  
  delete mesh;
  
}

TEST(InBoundaryTest, WithTri3)
{
  MeshIoMsh msh_reader;
  Mesh* mesh;
  int dim = 2;

  mesh =Mesh::create(TRIANGLE3,dim) ;
  
  msh_reader.readFileMsh("meshes/moderate_tri3.msh", mesh);

  int const n_nodes = mesh->numNodesTotal();
  Point *p;
  
  for (int i = 0; i < n_nodes; ++i)
  {
    p = mesh->getNode(i);
    EXPECT_TRUE(p->getIncidCell() >= 0);
  }
  

  int const in_boundary_nds[]  = {1,7,8,9,0,6,5,4,2,21,20,19,18,17,16,3,15,14,13,12,11,10};
  int const interior_nds[]     = {41,47,26,33,32,35,23,44,25,28,37,38,22,43,31,30,39,45,27,36,24,42,34,29,40,46};

  for (int i = 0; (unsigned)i < sizeof(in_boundary_nds)/sizeof(int); ++i)
  {
    p = mesh->getNode(in_boundary_nds[i]);
    EXPECT_TRUE(mesh->inBoundary(p)) << "node " << in_boundary_nds[i] << "; ic "<< p->getIncidCell() << std::endl;
  }
  
  for (int i = 0; (unsigned)i < sizeof(interior_nds)/sizeof(int); ++i)
  {
    p = mesh->getNode(interior_nds[i]);
    EXPECT_FALSE(mesh->inBoundary(p)) << "node " << interior_nds[i] << "; ic "<< p->getIncidCell() << std::endl;
  }

  delete mesh;
  
}


















