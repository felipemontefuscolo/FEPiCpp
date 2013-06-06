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
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>


using std::tr1::shared_ptr;
using std::tr1::tuple;
using std::tr1::get;
using std::vector;
using std::cout;
using std::endl;
using namespace Eigen;

void checkConsistencyTri(Mesh *mesh)
{
  cell_iterator cell = mesh->cellBegin();
  cell_iterator cell_end = mesh->cellEnd();
  
  int const nfpc = mesh->numFacetsPerCell();
  //  int const nnpc = mesh->numNodesPerCell();
  int const nnpf = mesh->numNodesPerFacet();
  
  int f_nds[nnpf];
  
  //  Point *p;
  Facet *f;
  Cell  *c;
  
  for (; cell != cell_end; ++cell)
  {
    int myid = mesh->getCellId(&*cell);
    for (int i = 0; i < nfpc; ++i)
    {
      if (cell->getIncidCell(i) >= 0)
      {
        // verifica o vizinho
        c = mesh->getCellPtr(cell->getIncidCell(i));
        int pos = cell->getIncidCellPos(i);
        EXPECT_EQ(myid, c->getIncidCell(pos))  << "myid=" << myid<<"; c->getIncidCell(pos)="<<c->getIncidCell(pos)<<"; i="<<i<<"; pos="<<pos;
        
        // verifica a face
        f = mesh->getFacetPtr(cell->getFacetId(i));
        int icf = f->getIncidCell();
        EXPECT_TRUE(icf==myid || icf==cell->getIncidCell(i))<<"myid="<<myid<<"; icf="<<icf<<"; i="<<i<<"; cell->getIncidCell(i)="<<cell->getIncidCell(i)<<"\n";
        if (icf==myid)
          EXPECT_TRUE(f->getPosition() == i) << "myid=" << myid<<"; f->getPosition()="<<f->getPosition()<<"; i="<<i<<"\n";
        else
        {
          EXPECT_TRUE(f->getPosition() == pos) << "myid=" << myid<<"; f->getPosition()="<<f->getPosition()<<"; pos="<<pos<<"; i="<<i<<"\n";
        }
      }
      else // bordo
      {
        // verifica a face
        f = mesh->getFacetPtr(cell->getFacetId(i));
        int icf = f->getIncidCell();
        // só pode ser o myid, pq do outro lado não tem ngm
        EXPECT_TRUE(icf==myid);
        
        // verifica se os nós da face estão no contorno
        mesh->getFacetNodesId(f, f_nds);
        for (int j = 0; j < nnpf; ++j)
        {
          EXPECT_TRUE(mesh->inBoundary(mesh->getNodePtr(f_nds[j])))
            << "node="<<f_nds[j]<<"; icell="<<mesh->getNodePtr(f_nds[j])->getIncidCell()
            << "; pos="<< mesh->getNodePtr(f_nds[j])->getPosition();
        }
        
      }
    }
    
  }
  
  std::vector<int> ics;
  std::vector<int> cc_ids;
  std::vector<int>::iterator it;
  for (point_iterator point = mesh->pointBegin(); point != mesh->pointEnd(); ++point)
  {
    point->getAllIncidences(ics);
    
    cc_ids.clear();

    int myid = point.index();
    
    for (int i = 0; i < (int)ics.size()/2; ++i)
    {
      int ic = ics.at(2*i);
      int pos = ics.at(2*i+1);
      Cell *c = mesh->getCellPtr(ic);
      ASSERT_TRUE(c);
      EXPECT_EQ(myid, c->getNodeId(pos)) << "ic = " << ic << "; pos = " << pos;
      cc_ids.push_back(c->getConnectedComponentId());
    }
    
    std::sort(cc_ids.begin(), cc_ids.end());
    it = unique (cc_ids.begin(), cc_ids.end());
    
    // checks if all incident cells have distinct connected component id
    ASSERT_EQ( std::distance(cc_ids.begin(), it) , (int)cc_ids.size());
    
  }
  
  
}

TEST(MtoolsFlipTest, WithTri3)
{
  
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh = NULL;  

  ECellType cell_t      = TRIANGLE3;
  const char* mesh_in  = "meshes/circle_tri3.msh";
  const char* mesh_out = "meshes/outtest/circle.vtk";
  
  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);
  
  vtk_printer.attachMesh(mesh);
  vtk_printer.isFamily(true);
  vtk_printer.setOutputFileName(mesh_out);
  
  checkConsistencyTri(mesh);
  
  vtk_printer.writeVtk();
  
  MeshToolsTri::flipEdge(mesh->getCellPtr( 3), 2, mesh);
  MeshToolsTri::flipEdge(mesh->getCellPtr( 6), 0, mesh);
  MeshToolsTri::flipEdge(mesh->getCellPtr(13), 1, mesh);
  
  checkConsistencyTri(mesh);
  
  vtk_printer.writeVtk();
  
  delete mesh;
  
  EXPECT_TRUE(true);
}

TEST(MtoolsFlipTest, WithTri6)
{
  
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh = NULL;  

  ECellType cell_t      = TRIANGLE6;
  const char* mesh_in  = "meshes/circle_tri6.msh";
  //const char* mesh_out = "meshes/circle_tri6.vtk";
  
  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);
  
  vtk_printer.attachMesh(mesh);
  //vtk_printer.isFamily(true);
  //vtk_printer.setOutputFileName(mesh_out);
  
  checkConsistencyTri(mesh);
  
  //vtk_printer.writeVtk();
  //
  MeshToolsTri::flipEdge(mesh->getCellPtr( 3), 2, mesh);
  MeshToolsTri::flipEdge(mesh->getCellPtr( 6), 0, mesh);
  MeshToolsTri::flipEdge(mesh->getCellPtr(13), 1, mesh);
  //
  checkConsistencyTri(mesh);
  
  //vtk_printer.writeVtk();
  
  delete mesh;
  
  EXPECT_TRUE(true);
}

TEST(MtoolsinCircle2dTest, WithTri3)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh = NULL;  

  ECellType cell_t      = TRIANGLE3;
  const char* mesh_in  = "meshes/simptri3.msh";
  //const char* mesh_out = "meshes/outtest/circle.vtk";
  
  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);

  Facet *f;
  int const n_facets_total = mesh->numFacetsTotal();
  
  for (int i = 0; i < n_facets_total; ++i)
  {
    f = mesh->getFacetPtr(i);
    EXPECT_TRUE(MeshToolsTri::inCircle2d(f, mesh));
  }

  Real coord[2] = {0.45, 0.18};
  Real coord_old[2];
  
  mesh->getNodePtr(8)->getCoord(coord_old,2);
  mesh->getNodePtr(8)->setCoord(coord,2);
  
  EXPECT_FALSE(MeshToolsTri::inCircle2d(mesh->getCellPtr(9), 1, mesh));
  EXPECT_FALSE(MeshToolsTri::inCircle2d(mesh->getCellPtr(2), 1, mesh));
  
  delete mesh;
  
  EXPECT_TRUE(true);
}

TEST(MtoolsFlippingMovingPointsTest, WithTri3)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh = NULL;  

  ECellType cell_t      = TRIANGLE3;
  const char* mesh_in  = "meshes/circle_tri3.msh";
  const char* mesh_out = "meshes/outtest/circle.vtk";
  
  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);

  vtk_printer.attachMesh(mesh);
  vtk_printer.isFamily(true);
  vtk_printer.setOutputFileName(mesh_out);
  
  checkConsistencyTri(mesh);
  
  Facet *f;
  Point *p;
  int const n_facets_total = mesh->numFacetsTotal();
  int const n_nodes_total  = mesh->numNodesTotal();
  double const dt = 0.1;
  
  // time loop
  for (double t=0; t<10; t+=dt)
  {
    //vtk_printer.writeVtk();
    checkConsistencyTri(mesh);
    
    // move points
    for (int i = 0; i < n_nodes_total; ++i)
    {
      double X[2];
      p = mesh->getNodePtr(i);
      p->getCoord(X,2);
      double a = 1. - sqrt(X[0]*X[0] + X[1]*X[1]);
      double Xnew[2] = {X[0] + dt*a*(-X[1]), X[1] + dt*a*X[0]};
      p->setCoord(Xnew,2);
    }
    
    // Delaunay
    for (int i = 0; i < n_facets_total; ++i)
    {
      f = mesh->getFacetPtr(i);
      if (!MeshToolsTri::inCircle2d(f, mesh))
      {
        MeshToolsTri::flipEdge(f, mesh);
      }
    }
    
  }
  //vtk_printer.writeVtk();

  delete mesh;
  
  EXPECT_TRUE(true);  
  
}

TEST(MtoolsFlippingMovingPointsTest, WithTri6)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh = NULL;  

  ECellType cell_t      = TRIANGLE6;
  const char* mesh_in  = "meshes/circle_tri6.msh";
  const char* mesh_out = "meshes/outtest/circle.vtk";
  
  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);

  vtk_printer.attachMesh(mesh);
  vtk_printer.isFamily(true);
  vtk_printer.setOutputFileName(mesh_out);
  
  checkConsistencyTri(mesh);
  
  Facet *f;
  Point *p;
  int const n_facets_total = mesh->numFacetsTotal();
  int const n_nodes_total  = mesh->numNodesTotal();
  double const dt = 0.1;
  
  // time loop
  for (double t=0; t<10; t+=dt)
  {
    //vtk_printer.writeVtk();
    checkConsistencyTri(mesh);
    
    // move points
    for (int i = 0; i < n_nodes_total; ++i)
    {
      double X[2];
      p = mesh->getNodePtr(i);
      p->getCoord(X,2);
      double a = 1. - sqrt(X[0]*X[0] + X[1]*X[1]);
      double Xnew[2] = {X[0] + dt*a*(-X[1]), X[1] + dt*a*X[0]};
      p->setCoord(Xnew,2);
    }
    
    // Delaunay
    for (int i = 0; i < n_facets_total; ++i)
    {
      f = mesh->getFacetPtr(i);
      if (!MeshToolsTri::inCircle2d(f, mesh))
      {
        MeshToolsTri::flipEdge(f, mesh, true);
      }
    }
    
  }
  //vtk_printer.writeVtk();

  delete mesh;
  
  EXPECT_TRUE(true);  
  
}

TEST(MtoolsReadMeshTest, WithTri3)
{
  Mesh *mesh = NULL;
  mesh = Mesh::create(TRIANGLE3);
  
  MeshTools mtools;
  MeshIoVtk vtk_printer;
  
  vtk_printer.attachMesh(mesh);
  
  int const nodes[] = {0,1,3,
                        1,4,3,
                        1,2,4,
                        0,3,5,
                        4,2,7,
                        4,7,6,
                        3,4,6,
                        3,6,5};
                        
  Real const xyz[] = { 0.0, 0.0,
                        0.5, 0.0,
                        1.0, 0.0,
                        0.3, 0.5,
                        0.6, 0.5,
                        0.0, 1.0,
                        0.5, 1.0,
                        1.0, 1.0 };

  int const nnpc = 3;
  int const n_cells = sizeof(nodes)/sizeof(int)/nnpc;
  int const n_nodes = 1 + *std::max_element(nodes, nodes + n_cells*nnpc);
  
  mtools.readMesh(n_nodes, n_cells, nodes, xyz, mesh);
  
  checkConsistencyTri(mesh);
  
  vtk_printer.setOutputFileName("meshes/outtest/array2mesh.vtk");
  vtk_printer.writeVtk();
  
  delete mesh;
}

TEST(MtoolInsertVertexOnEdgeTest, WithTri3)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh = NULL;  
  MeshToolsTri mtools;

  ECellType cell_t      = TRIANGLE3;
  const char* mesh_in  = "meshes/complex_tri3.msh";
  const char* mesh_out = "meshes/outtest/complex.vtk";
  
  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);

  int const n_facets = mesh->numFacetsTotal();

  vtk_printer.isFamily(true);
  vtk_printer.attachMesh(mesh);
  vtk_printer.setOutputFileName(mesh_out);

  int tag;


  for (int i = 0; i < n_facets; ++i)
  //for (int i = 3; i < 5; ++i)
  {
    Facet* edge = mesh->getFacetPtr(i);
    
    if (edge->isDisabled()) return;

    //checkConsistencyTri(mesh);

    // old tag
    tag = edge->getTag();

    EXPECT_TRUE(tag > 0);

    //cout << "cell: "<< mesh->getCellId(mesh->getCellPtr(edge->getIncidCell()))<<"; node_id:"<< mesh->getCellPtr(edge->getIncidCell())->getNodeId(edge->getPosition())
    //     << "; edge position:" << edge->getPosition() <<"; numm_nodes: " <<mesh->numNodes()<< endl;
    //mtools.insertVertexOnEdge(edge, 0.5, mesh);
    //vtk_printer.writeVtk();
    int nd_id = mtools.insertVertexOnEdge(edge->getIncidCell(), edge->getPosition(), 0.5, mesh);
    EXPECT_TRUE(nd_id > 0);
    if (mesh->inBoundary(edge))
      EXPECT_TRUE(mesh->inBoundary( mesh->getNodePtr(nd_id) ));
    
    EXPECT_EQ(tag, mesh->getNodePtr(nd_id)->getTag());
    
    checkConsistencyTri(mesh);
  }
  vtk_printer.writeVtk();
  vtk_printer.printPointTagVtk("node_tags");
  
  
  delete mesh;
}

TEST(MtoolCollapseEdge2dTest, WithTri3)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh = NULL;  
  MeshToolsTri mtools;

  ECellType cell_t      = TRIANGLE3;
  const char* mesh_in  = "meshes/complex_tri3.msh";
  const char* mesh_out = "meshes/outtest/collapse.vtk";
  
  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);

  int const n_facets = mesh->numFacetsTotal();

  vtk_printer.isFamily(true);
  vtk_printer.attachMesh(mesh);
  vtk_printer.setOutputFileName(mesh_out);

  int nd_id = mtools.collapseEdge2d(96,0,.5,mesh);
  checkConsistencyTri(mesh);
  nd_id = mtools.collapseEdge2d(109,0,.5,mesh);
  checkConsistencyTri(mesh);
  nd_id = mtools.collapseEdge2d(102,2,.5,mesh);
  checkConsistencyTri(mesh);
  nd_id = mtools.collapseEdge2d(114,0,.5,mesh);
  checkConsistencyTri(mesh);
  

  vtk_printer.writeVtk();
  
  
  delete mesh;
}


TEST(ImprimeTest, TestTest)
{
  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh = NULL;  
  //MeshToolsTri mtools;

  ECellType cell_t      = TRIANGLE6;
  const char* mesh_in  = "meshes/1level_tri6.msh";
  const char* mesh_out = "meshes/1level_tri6.vtk";
  
  mesh = Mesh::create(cell_t);
  msh_reader.readFileMsh(mesh_in, mesh);

  //int const n_facets = mesh->numFacetsTotal();

  //vtk_printer.isFamily(true);
  vtk_printer.attachMesh(mesh);
  vtk_printer.setOutputFileName(mesh_out);
  vtk_printer.writeVtk();
  
  delete mesh;
}

#if 0

//
//TEST(searchConvexPointTest, Tri3Test)
//{
//  MeshIoMsh msh_reader;
//  MeshIoVtk vtk_printer;
//  Mesh *mesh = NULL;  
//
//  ECellType cell_t      = TRIANGLE3;
//  const char* mesh_in  = "meshes/uni_tri.msh";
//  //const char* mesh_out = "meshes/uni_tri.vtk";
//  
//  mesh = Mesh::create(cell_t);
//  msh_reader.readFileMsh(mesh_in, mesh);
//  
//  //                       cell     x        y
//  typedef std::tr1::tuple<int, double, double> location;
//  
//  location loc[] = {location( 66, 0.08, 0.55),
//                    location( 28, 0.97, 0.01),
//                    location( 27, 0.01, 0.98),
//                    location(375, 0.04, 0.01),
//                    location(376, 0.01, 0.03),
//                    location( 80, 0.30, 0.54),
//                    location( -1, 0.50, 0.53),
//                    location( -1, 0.00, 1.00000001),
//                    location( -1, 1.00000001, 0.0),
//                    location( 28, 1.00, 0.00)};
//  
//  std::pair<bool, Cell *> bc;
//  double X[2];
//  int cid;
//  
//  /* initialize random seed: */
//  srand ( 180387 );  
//  
//  
//  for (int jj=0; jj<10; ++jj)
//  {
//    for (uint i = 0; i < sizeof(loc)/sizeof(location); ++i)
//    //for (uint i = 0; i < 1; ++i)
//    {
//      //cout << get<0>(loc[i]) << " " << get<1>(loc[i]) << " " << get<2>(loc[i]) << endl;
//      X[0] = get<1>(loc[i]);
//      X[1] = get<2>(loc[i]);
//      bc = MeshToolsTri::searchConvexPoint(X, mesh->getCellPtr(rand()%mesh->numCells()), mesh);
//      if (bc.first)
//        cid = mesh->getCellId(bc.second);
//      else
//        cid = -1;
//      
//      EXPECT_EQ(get<0>(loc[i]), cid);
//    }
//  }
//  
//  checkConsistencyTri(mesh);
//  
//  delete mesh;
//}
//
//class SomeFunc
//{public:
//  Real operator() (Real const*x) const
//  {
//    return x[1] - SomeFunc::Y();
//  }
//  static double Y()
//  {
//    return 0.03;
//  }
//};
//
//
//TEST(CreatePathTest, Tri3Test)
//{
//  MeshIoMsh msh_reader;
//  MeshIoVtk vtk_printer;
//  Mesh *mesh = NULL;  
//
//  ECellType cell_t      = TRIANGLE3;
//  const char* mesh_in  = "meshes/uni_tri.msh";
//  //const char* mesh_out = "meshes/uni_tri.vtk";
//  
//  mesh = Mesh::create(cell_t);
//  msh_reader.readFileMsh(mesh_in, mesh);
//  
//  Real x0[] = {0.015, 0.005};
//  
//  std::vector<int> cells_id;
//  std::vector<Real> slices;
//  
//  Cell *c0 = mesh->getCellPtr(313);
//  
//  MeshToolsTri::createPath(x0, c0, SomeFunc(), cells_id, slices,  mesh);
//  
//  EXPECT_EQ(38u, cells_id.size());
//  EXPECT_EQ(114u, slices.size());
//
//  EXPECT_EQ(28, cells_id.back());
//  
//  
//  // conferi os vizinhos
//  for (int c = 0; c < (int)cells_id.size(); ++c)
//  {
//    c0 = mesh->getCellPtr(cells_id[c]);
//    int myid = mesh->getCellId(c0);
//    
//    int f_free = -1;
//    for (int i = 0; i < 3; ++i) {if (slices[3*c + i]<0) {f_free = i; break;}}
//    
//    EXPECT_TRUE(f_free>=0) << cout<<" SLICE: "<<slices[3*c]<<" "<<slices[3*c + 1]<<" "<<slices[3*c + 2]<<std::endl ;
//    
//    for (int kk = 0; kk < 2; ++kk) // sides
//    {
//      int f = (f_free+kk+1)%3;
//      
//      int oc_id = c0->getIncidCell(f);
//      if (oc_id < 0)
//        break;
//      int pos = c0->getIncidCellPos(f);
//      
//      std::vector<int>::iterator it = find(cells_id.begin(), cells_id.end(), oc_id);
//      EXPECT_FALSE(it == cells_id.end()) << cout << "ID: " << oc_id << ", MY: " << myid;
//      
//      EXPECT_EQ(myid, mesh->getCellPtr(oc_id)->getIncidCell(pos));
//      
//      int octh = std::distance(cells_id.begin(),it);
//      EXPECT_DOUBLE_EQ(1.0, slices[3*octh+pos]+slices[3*c+f])
//        << cout << "\n\tOslices= "<<slices[3*octh+0]<<" "<<slices[3*octh+1]<<" "<<slices[3*octh+2]<<"\n"
//                << "\t slices= "<<slices[3*c+0]<<" "<<slices[3*c+1]<<" "<<slices[3*c+2]<<"\n"
//                << "\t pos= " << pos << " " << " f=" << f << " my=" << myid << " oid=" << oc_id << endl;
//      
//    
//          
//    }
//    
//  }
//  
//  
//  //for (int i = 0; i < slices.size(); ++i)
//  //{
//  //  std::cout << slices[i] << std::endl;
//  //}
//  
//  
//  checkConsistencyTri(mesh);
//  
//  delete mesh;
//}
//
//class SomeFunc2
//{public:
//  Real operator() (Real const*x) const
//  {
//    return x[1] - SomeFunc2::Y();
//  }
//  static double Y()
//  {
//    return 0.48;
//  }
//};
//
//
//TEST(CreatePathTest2, Tri3Test)
//{
//  MeshIoMsh msh_reader;
//  MeshIoVtk vtk_printer;
//  Mesh *mesh = NULL;  
//
//  ECellType cell_t      = TRIANGLE3;
//  const char* mesh_in  = "meshes/uni_tri.msh";
//  //const char* mesh_out = "meshes/uni_tri.vtk";
//  
//  mesh = Mesh::create(cell_t);
//  msh_reader.readFileMsh(mesh_in, mesh);
//  
//  Real x0[] = {SomeFunc2::Y(), SomeFunc2::Y()};
//  
//  std::vector<int> cells_id;
//  std::vector<Real> slices;
//  
//  Cell *c0 = mesh->getCellPtr(246);
//  
//  MeshToolsTri::createPath(x0, c0, SomeFunc2(), cells_id, slices,  mesh);
//
//
//
//
//  
//  // conferi os vizinhos
//  for (int c = 0; c < (int)cells_id.size(); ++c)
//  {
//    c0 = mesh->getCellPtr(cells_id[c]);
//    int myid = mesh->getCellId(c0);
//    
//    int f_free = -1;
//    for (int i = 0; i < 3; ++i) {if (slices[3*c + i]<0) {f_free = i; break;}}
//    
//    EXPECT_TRUE(f_free>=0) << cout<<" SLICE: "<<slices[3*c]<<" "<<slices[3*c + 1]<<" "<<slices[3*c + 2]<<std::endl ;
//    
//    for (int kk = 0; kk < 2; ++kk) // sides
//    {
//      int f = (f_free+kk+1)%3;
//      
//      int oc_id = c0->getIncidCell(f);
//      if (oc_id < 0)
//        break;
//      int pos = c0->getIncidCellPos(f);
//      
//      std::vector<int>::iterator it = find(cells_id.begin(), cells_id.end(), oc_id);
//      EXPECT_FALSE(it == cells_id.end()) << cout << "ID: " << oc_id << ", MY: " << myid;
//      
//      EXPECT_EQ(myid, mesh->getCellPtr(oc_id)->getIncidCell(pos));
//      
//      int octh = std::distance(cells_id.begin(),it);
//      EXPECT_DOUBLE_EQ(1.0, slices[3*octh+pos]+slices[3*c+f])
//        << cout << "\n\tOslices= "<<slices[3*octh+0]<<" "<<slices[3*octh+1]<<" "<<slices[3*octh+2]<<"\n"
//                << "\t slices= "<<slices[3*c+0]<<" "<<slices[3*c+1]<<" "<<slices[3*c+2]<<"\n"
//                << "\t pos= " << pos << " " << " f=" << f << " my=" << myid << " oid=" << oc_id << endl;
//
//    }
//    
//  }
//  
//  
//  delete mesh;
//}
//
//TEST(cutConvexPartTest, Tri3Test)
//{
//  MeshIoMsh msh_reader;
//  MeshIoVtk vtk_printer;
//  Mesh *mesh = NULL;  
//
//  ECellType cell_t      = TRIANGLE3;
//  const char* mesh_in  = "meshes/uni_tri.msh";
//  //const char* mesh_out = "meshes/uni_tri.vtk";
//  
//  mesh = Mesh::create(cell_t);
//  msh_reader.readFileMsh(mesh_in, mesh);
//  
//  Real x0[] = {SomeFunc2::Y(), SomeFunc2::Y()};
//  
//  Cell *c0 = mesh->getCellPtr(246);
//  
//  MeshToolsTri::cutConvexPart(x0, c0, SomeFunc2(), mesh);
//
//  
//  delete mesh;
//}
//
//




//
//TEST(AssignsDofsTest, Tri3Test)
//{
//  MeshIoMsh msh_reader;
//  MeshIoVtk vtk_printer;
//  Mesh *mesh = NULL;  
//
//  ECellType cell_t     = TRIANGLE3;
//  const char* mesh_in  = "meshes/simptri3.msh";  
//  //const char* mesh_out = "meshes/outtest/dof_tri3.vtk";
//
//  mesh = Mesh::create(cell_t);
//  msh_reader.readFileMsh(mesh_in, mesh);
//  
//  DofHandler DofH(mesh);
//  //                         ndpv,  ndpr,  ndpf,  ndpc
//  DofH.addVariable("altura",    1,     0,     0,     0); // 12
//  DofH.addVariable("vetor",     2,     0,     0,     0); // 24
//  DofH.addVariable("foo",       0,     1,     0,     0); // 12
//  DofH.addVariable("bar",       0,     0,     1,     0); // 25
//  DofH.addVariable("moo",       0,     0,     0,     1); // 14
//  
//  DofH.SetUp();
//  
//  EXPECT_EQ(87, DofH.numDofs());
//  
//  
//  MeshTools::removeCell(mesh->getCellPtr(2), mesh);
//  MeshTools::removeCell(mesh->getCellPtr(3), mesh);
//  
//  DofH.SetUp();
//  
//  //11x3 + 22 + 12
//  EXPECT_EQ(11, DofH.getVariable(0).numDofs());
//  EXPECT_EQ(22, DofH.getVariable(1).numDofs());
//  EXPECT_EQ(11, DofH.getVariable(2).numDofs());
//  EXPECT_EQ(22, DofH.getVariable(3).numDofs());
//  EXPECT_EQ(12, DofH.getVariable(4).numDofs());
//
//  EXPECT_EQ(78, DofH.numDofs());
//  
//  
//  // .getVariable(0)
//  int *dat = DofH.data();
//
//  int counter = 0;
//  for (int i = 0; i < DofH.totalSize(); ++i)
//  {
//    if ((*dat) >= 0)
//      EXPECT_EQ(counter++, (*dat));
//      
//    //std::cout << (*dat++) << std::endl;
//  }
//  
//  
//  delete mesh;
//}
//
//class MyGetDataVtkTri6Test : public DefaultGetDataVtk
//{
//public:
//  MyGetDataVtkTri6Test(int * i = NULL, Mesh* m=NULL, DofHandler* df=NULL) :
//                      DefaultGetDataVtk(NULL,i), mesh_ptr(m), dofh_ptr(df){}
//  
//  int get_data_i(int nodeid) const
//  {
//    int dof;
//    Point const* p = mesh_ptr->getNodePtr(nodeid);
//    if (!mesh_ptr->isVertex(p))
//      return -1;
//    dofh_ptr->getVariable(0).getVertexDofs(&dof, p);
//    return data_i[dof];
//  }
//  
//  Mesh *mesh_ptr;
//  DofHandler *dofh_ptr;
//  
//  virtual ~MyGetDataVtkTri6Test() {}
//  
//  // int  * data_i; from MyGetDataVtk
//};
//
//TEST(AssignsDofsTest, Tri6Test)
//{
//  MeshIoMsh msh_reader;
//  MeshIoVtk vtk_printer;
//  Mesh *mesh = NULL;  
//
//  ECellType cell_t     = TRIANGLE6;
//  const char* mesh_in  = "meshes/simptri6.msh";  
//  const char* mesh_out = "meshes/outtest/dof_tri6.vtk";
//
//  mesh = Mesh::create(cell_t);
//  msh_reader.readFileMsh(mesh_in, mesh);
//  
//  DofHandler DofH(mesh);
//  //                         ndpv,  ndpr,  ndpf,  ndpc
//  DofH.addVariable("numero",    1,     0,     0,     0); // 12
//  DofH.addVariable("vetor",     2,     0,     0,     0); // 24
//  DofH.addVariable("foo",       0,     1,     0,     0); // 12
//  DofH.addVariable("bar",       0,     0,     1,     0); // 25
//  DofH.addVariable("moo",       0,     0,     0,     1); // 14
//  
//  DofH.SetUp();
//  
//  EXPECT_EQ(87, DofH.numDofs());
//  
//  MeshTools::removeCell(mesh->getCellPtr(2), mesh);
//  MeshTools::removeCell(mesh->getCellPtr(3), mesh);
//  
//  DofH.SetUp();
//  
//  //11x3 + 22 + 12 + 12
//  EXPECT_EQ(11, DofH.getVariable(0).numDofs());
//  EXPECT_EQ(22, DofH.getVariable(1).numDofs());
//  EXPECT_EQ(11, DofH.getVariable(2).numDofs());
//  EXPECT_EQ(22, DofH.getVariable(3).numDofs());
//  EXPECT_EQ(12, DofH.getVariable(4).numDofs());
//
//  EXPECT_EQ(78, DofH.numDofs());
//  
//  
//  // .getVariable(0)
//  int *dat = DofH.data();
//
//  int counter = 0;
//  for (int i = 0; i < DofH.totalSize(); ++i)
//  {
//    if ((*dat) >= 0)
//      EXPECT_EQ(counter++, (*dat));
//      
//  }
//  
//  vtk_printer.attachMesh(mesh);
//  vtk_printer.writeVtk(mesh_out);
//  
//  std::vector<int> dados(DofH.numDofs());
//  
//  point_iterator pt = mesh->pointBegin();
//  point_iterator pt_end = mesh->pointEnd();
//  
//  // for variable "0" ("numero")
//  for (; pt != pt_end; ++pt)
//  {
//    if (!mesh->isVertex(&*pt))
//      continue;
//      
//    int dof(0); // "numero"
//    DofH.getVariable(0).getVertexDofs(&dof, &*pt);
//    dados[dof] = mesh->getNodeContigId(mesh->getPointId(&*pt));
//    //dados[dof] = 0;
//  }
//  
//  MyGetDataVtkTri6Test get_data(dados.data(), mesh, &DofH);
//  
//  vtk_printer.addNodeIntVtk(DofH.getVariable(0).getName(), get_data);
//  
//  delete mesh;
//}
//
//
//
//
//
//TEST(AssignsDofsTest, Tet10Test)
//{
//  MeshIoMsh msh_reader;
//  MeshIoVtk vtk_printer;
//  Mesh *mesh = NULL;  
//
//  ECellType cell_t     = TETRAHEDRON10;
//  const char* mesh_in  = "meshes/simptet10.msh";  
//  const char* mesh_out = "meshes/outtest/dof_tet10.vtk";
//
//  mesh = Mesh::create(cell_t);
//  msh_reader.readFileMsh(mesh_in, mesh);
//  
//  DofHandler DofH(mesh);
//  //                         ndpv,  ndpr,  ndpf,  ndpc
//  DofH.addVariable("numero",    1,     0,     0,     0); //  9
//  DofH.addVariable("foo",       0,     1,     0,     0); // 26
//  DofH.addVariable("bar",       0,     0,     1,     0); // 30
//  DofH.addVariable("moo",       0,     0,     0,     2); // 24
//  
//  DofH.SetUp();
//  
//  EXPECT_EQ( 9, DofH.getVariable(0).numDofs());
//  EXPECT_EQ(26, DofH.getVariable(1).numDofs());
//  EXPECT_EQ(30, DofH.getVariable(2).numDofs());
//  EXPECT_EQ(24, DofH.getVariable(3).numDofs());
//  
//  EXPECT_EQ(89, DofH.numDofs());
//  
//  //DofH.SetUp();
//  //
//  MeshTools::removeCell(mesh->getCellPtr(2), mesh);
//  MeshTools::removeCell(mesh->getCellPtr(3), mesh);
//  //
//  DofH.SetUp();
//  //
//  
//  EXPECT_EQ( 9, DofH.getVariable(0).numDofs());
//  EXPECT_EQ(26, DofH.getVariable(1).numDofs());
//  EXPECT_EQ(28, DofH.getVariable(2).numDofs());
//  EXPECT_EQ(20, DofH.getVariable(3).numDofs());
//  ////
//  EXPECT_EQ(83, DofH.numDofs());
//  //
//  //
//  //// .getVariable(0)
//  //int *dat = DofH.data();
//  //
//  //int counter = 0;
//  //for (int i = 0; i < DofH.totalSize(); ++i)
//  //{
//  //  if ((*dat) >= 0)
//  //    EXPECT_EQ(counter++, (*dat));
//  //    
//  //}
//  //
//  vtk_printer.attachMesh(mesh);
//  vtk_printer.writeVtk(mesh_out);
//  //
//  std::vector<int> dados(DofH.numDofs());
//  
//  point_iterator pt = mesh->pointBegin();
//  point_iterator pt_end = mesh->pointEnd();
//  
//  // for variable "0" ("numero")
//  for (; pt != pt_end; ++pt)
//  {
//    if (!mesh->isVertex(&*pt))
//      continue;
//      
//    int dof(0); // "numero"
//    DofH.getVariable(0).getVertexDofs(&dof, &*pt);
//    dados[dof] = mesh->getNodeContigId(mesh->getPointId(&*pt));
//    //dados[dof] = 0;
//  }
//  
//  MyGetDataVtkTri6Test get_data(dados.data(), mesh, &DofH);
//  
//  vtk_printer.addNodeIntVtk(DofH.getVariable(0).getName(), get_data);
//  
//  delete mesh;
//}
//
//
//TEST(BubbleTri3Test, Tri3Test)
//{
//  MeshIoMsh msh_reader;
//  MeshIoVtk vtk_printer;
//  Mesh *mesh = NULL;  
//
//  ECellType cell_t     = TRIANGLE3;
//  const char* mesh_in  = "meshes/simptri3.msh";  
//  const char* mesh_out = "meshes/outtest/metis_tri3.vtk";
//
//  ShapeFunction *phi, *psi;
//  phi = ShapeFunction::create(TRIANGLE3,    P2);
//  psi = ShapeFunction::create(TRIANGLE3,    P0);
//
//  //std::cout << phi->numDofsAssociatedToVertice() << " " <<
//  //             phi->numDofsAssociatedToCorner()  << " " <<
//  //             phi->numDofsAssociatedToFacet()   << " " <<
//  //             phi->numDofsAssociatedToCell() << std::endl;
//
//  mesh = Mesh::create(cell_t);
//  msh_reader.readFileMsh(mesh_in, mesh);
//  
//  DofHandler DofH(mesh);
//  
//  DofH.addVariable("velo",   phi, 2);
//  DofH.addVariable("press",  psi, 1);
//  
//  //MeshTools::removeCell(mesh->getCellPtr(0), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(1), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(10), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(6), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(7), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(9), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(2), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(3), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(12), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(4), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(5), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(8), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(11), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(13), mesh);
//  
//  DofH.SetUp();
//  std::cout << "num dofs = " << DofH.numDofs() << std::endl;
//  
//  
//  //EXPECT_EQ(12, DofH.getVariable(0).numDofs());
//  //EXPECT_EQ( 4, DofH.getVariable(1).numDofs());
//  //////
//  //EXPECT_EQ(16, DofH.numDofs());  
//  
//  //bool relations[2][2] = {{1,1}, {1,1}};
//  //DofH.setVariablesRelationship(*relations);
//  
//  //DofH.metisRenumber();
//  //DofH.boostMinimumDegreeRenumber();
//  //DofH.boostCuthillMcKeeRenumber();
//  DofH.CuthillMcKeeRenumber();
//  
//  //ArrayXi var_cell_dofs(8);
//  //ArrayXi var_cell_dofs2(3);
//  //
//  //DofH.getVariable(0).getCellDofs(var_cell_dofs.data(), mesh->getCellPtr(13));
//  //DofH.getVariable(1).getCellDofs(var_cell_dofs2.data(), mesh->getCellPtr(13));
//  //
//  //std::cout << var_cell_dofs.transpose() << std::endl;
//  //std::cout << var_cell_dofs2.transpose() << std::endl;
//
//  
//  DofH.printSparsityMatlab("Jac");
//  vtk_printer.attachMesh(mesh);
//  vtk_printer.writeVtk(mesh_out);
//  
//  std::sort(DofH.data(), DofH.data()+DofH.totalSize());
//  int counter = 0;
//  for (int i = 0; i < DofH.totalSize(); ++i)
//  {
//    if (DofH.data() < 0) continue;
//    EXPECT_EQ(counter, DofH.data()[i]);
//    counter++;
//  }
//  
//  delete phi;
//  delete psi;
//  delete mesh;
//}
//
//
//
//TEST(BubbleTet4Test, Tet4Test)
//{
//  MeshIoMsh msh_reader;
//  MeshIoVtk vtk_printer;
//  Mesh *mesh = NULL;  
//
//  ECellType cell_t     = TRIANGLE3;
//  const char* mesh_in  = "meshes/simptri3.msh";  
//  const char* mesh_out = "meshes/outtest/metis_tet4.vtk";
//
//  ShapeFunction *phi, *psi;
//  phi = ShapeFunction::create(TETRAHEDRON4,    P1ph);
//  psi = ShapeFunction::create(TETRAHEDRON4,    P1);
//
//  //std::cout << phi->numDofsAssociatedToVertice() << " " <<
//  //             phi->numDofsAssociatedToCorner()  << " " <<
//  //             phi->numDofsAssociatedToFacet()   << " " <<
//  //             phi->numDofsAssociatedToCell() << std::endl;
//
//  mesh = Mesh::create(cell_t);
//  msh_reader.readFileMsh(mesh_in, mesh);
//  
//  DofHandler DofH(mesh);
//  
//  DofH.addVariable("velo",   phi, 2);
//  DofH.addVariable("press",  psi, 1);
//  
//  //MeshTools::removeCell(mesh->getCellPtr(0), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(1), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(10), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(6), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(7), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(9), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(2), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(3), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(12), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(4), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(5), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(8), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(11), mesh);
//  //MeshTools::removeCell(mesh->getCellPtr(13), mesh);
//  
//  DofH.SetUp();
//  std::cout << "num dofs = " << DofH.numDofs() << std::endl;
//  
//  
//  //EXPECT_EQ(12, DofH.getVariable(0).numDofs());
//  //EXPECT_EQ( 4, DofH.getVariable(1).numDofs());
//  //////
//  //EXPECT_EQ(16, DofH.numDofs());  
//  
//  //bool relations[2][2] = {{1,1}, {1,0}};
//  //DofH.setVariablesRelationship(*relations);
//  
//  //DofH.metisRenumber();
//  //DofH.boostMinimumDegreeRenumber();
//  //DofH.boostCuthillMcKeeRenumber();
//  DofH.CuthillMcKeeRenumber();
//  
//  //ArrayXi var_cell_dofs(8);
//  //ArrayXi var_cell_dofs2(3);
//  //
//  //DofH.getVariable(0).getCellDofs(var_cell_dofs.data(), mesh->getCellPtr(13));
//  //DofH.getVariable(1).getCellDofs(var_cell_dofs2.data(), mesh->getCellPtr(13));
//  //
//  //std::cout << var_cell_dofs.transpose() << std::endl;
//  //std::cout << var_cell_dofs2.transpose() << std::endl;
//
//  
//  DofH.printSparsityMatlab("Jac");
//  vtk_printer.attachMesh(mesh);
//  vtk_printer.writeVtk(mesh_out);
//  
//  std::sort(DofH.data(), DofH.data()+DofH.totalSize());
//  int counter = 0;
//  for (int i = 0; i < DofH.totalSize(); ++i)
//  {
//    if (DofH.data() < 0) continue;
//    EXPECT_EQ(counter, DofH.data()[i]);
//    counter++;
//  }
//  
//  delete phi;
//  delete psi;
//  delete mesh;
//}
//
//
///*
// * 
// * 
// *  TAGS TEST
// * 
// * 
// */ 
//
//
//
//TEST(TagsDofsTest, Tri3Test)
//{
//  MeshIoMsh msh_reader;
//  MeshIoVtk vtk_printer;
//  Mesh *mesh = NULL;  
//
//  ECellType cell_t     = TRIANGLE3;
//  const char* mesh_in  = "meshes/simptri3.msh";  
//  //const char* mesh_out = "meshes/outtest/dof_tri3.vtk";
//
//  mesh = Mesh::create(cell_t);
//  msh_reader.readFileMsh(mesh_in, mesh);
//  
//  int ntags=1;
//  int tags=1;
//  
//  DofHandler DofH(mesh);
//  //                         ndpv,  ndpr,  ndpf,  ndpc
//  DofH.addVariable("altura",    1,     0,     0,     0, ntags, &tags); // 4
//  DofH.addVariable("vetor",     2,     0,     0,     0, ntags, &tags); // 8
//  DofH.addVariable("foo",       0,     1,     0,     0, ntags, &tags); // 4
//  DofH.addVariable("bar",       0,     0,     1,     0, ntags, &tags); // 4
//  DofH.addVariable("moo",       0,     0,     0,     1, ntags, &tags); // 0
//  
//  
//  DofH.SetUp();
//  
//  EXPECT_EQ(20, DofH.numDofs());
//  
//  
//  MeshTools::removeCell(mesh->getCellPtr(2), mesh);
//  MeshTools::removeCell(mesh->getCellPtr(3), mesh);
//  
//  DofH.SetUp();
//  
//  // 15
//  EXPECT_EQ( 3, DofH.getVariable(0).numDofs());
//  EXPECT_EQ( 6, DofH.getVariable(1).numDofs());
//  EXPECT_EQ( 3, DofH.getVariable(2).numDofs());
//  EXPECT_EQ( 3, DofH.getVariable(3).numDofs());
//  EXPECT_EQ( 0, DofH.getVariable(4).numDofs());
//
//  EXPECT_EQ(15, DofH.numDofs());
//  
//  
//  // .getVariable(0)
//  int *dat = DofH.data();
//
//  int counter = 0;
//  for (int i = 0; i < DofH.totalSize(); ++i)
//  {
//    if ((*dat) >= 0)
//      EXPECT_EQ(counter++, (*dat));
//      
//    //std::cout << (*dat++) << std::endl;
//  }
//  
//  
//  delete mesh;
//}
//
//
//
//
//TEST(TagsDofsTest, Tet10Test)
//{
//  MeshIoMsh msh_reader;
//  MeshIoVtk vtk_printer;
//  Mesh *mesh = NULL;  
//
//  ECellType cell_t     = TETRAHEDRON10;
//  const char* mesh_in  = "meshes/simptet10.msh";  
//  const char* mesh_out = "meshes/outtest/dof_tet10.vtk";
//
//  mesh = Mesh::create(cell_t);
//  msh_reader.readFileMsh(mesh_in, mesh);
//
//  
//  int ntags=1;
//  int tags=1;
//  
//  DofHandler DofH(mesh);
//  //                         ndpv,  ndpr,  ndpf,  ndpc
//  DofH.addVariable("numero",    1,     0,     0,     0, ntags, &tags); //  9
//  DofH.addVariable("foo",       0,     1,     0,     0, ntags, &tags); // 26
//  DofH.addVariable("bar",       0,     0,     1,     0, ntags, &tags); // 30
//  DofH.addVariable("moo",       0,     0,     0,     2, ntags, &tags); // 24
//
//  DofH.SetUp();
//  
//  EXPECT_EQ( 9, DofH.getVariable(0).numDofs());
//  EXPECT_EQ(26, DofH.getVariable(1).numDofs());
//  EXPECT_EQ(30, DofH.getVariable(2).numDofs());
//  EXPECT_EQ(24, DofH.getVariable(3).numDofs());
//  
//  EXPECT_EQ(89, DofH.numDofs());
//  
//  //DofH.SetUp();
//  //
//  MeshTools::removeCell(mesh->getCellPtr(2), mesh);
//  MeshTools::removeCell(mesh->getCellPtr(3), mesh);
//  //
//  DofH.SetUp();
//  //
//  
//  EXPECT_EQ( 9, DofH.getVariable(0).numDofs());
//  EXPECT_EQ(26, DofH.getVariable(1).numDofs());
//  EXPECT_EQ(28, DofH.getVariable(2).numDofs());
//  EXPECT_EQ(20, DofH.getVariable(3).numDofs());
//  ////
//  EXPECT_EQ(83, DofH.numDofs());
//  //
//  //
//  //// .getVariable(0)
//  //int *dat = DofH.data();
//  //
//  //int counter = 0;
//  //for (int i = 0; i < DofH.totalSize(); ++i)
//  //{
//  //  if ((*dat) >= 0)
//  //    EXPECT_EQ(counter++, (*dat));
//  //    
//  //}
//  //
//  vtk_printer.attachMesh(mesh);
//  vtk_printer.writeVtk(mesh_out);
//  //
//  std::vector<int> dados(DofH.numDofs());
//  
//  point_iterator pt = mesh->pointBegin();
//  point_iterator pt_end = mesh->pointEnd();
//  
//  // for variable "0" ("numero")
//  for (; pt != pt_end; ++pt)
//  {
//    if (!mesh->isVertex(&*pt))
//      continue;
//      
//    int dof(0); // "numero"
//    DofH.getVariable(0).getVertexDofs(&dof, &*pt);
//    dados[dof] = mesh->getNodeContigId(mesh->getPointId(&*pt));
//    //dados[dof] = 0;
//  }
//  
//  MyGetDataVtkTri6Test get_data(dados.data(), mesh, &DofH);
//  
//  vtk_printer.addNodeIntVtk(DofH.getVariable(0).getName(), get_data);
//  
//  delete mesh;
//}
//
//








#endif



