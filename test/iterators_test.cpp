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
#include <vector>
#include <tr1/tuple>
#include <algorithm>

typedef std::tr1::tuple<ECellType, const char*> TupleT;

using std::cout;
using std::endl;

class IteratorTest : public testing::TestWithParam<TupleT> { // vertex id
protected:

  //virtual ~IteratorTest() {delete mesh;}
  
  virtual void SetUp() {
    
    MeshIoMsh msh_reader; // TEM QUE SER AQUI ... ESTRANHO¹²³
    
    //Mesh *mesh;
    
    ECellType cell_t     = std::tr1::get<0>(GetParam());
    const char* mesh_in  = std::tr1::get<1>(GetParam());
    
    mesh = Mesh::create(cell_t);
    //mesh = new SMesh<Triangle3,2>;
    mesh->qBuildAdjacency(false);
    msh_reader.readFileMsh(mesh_in, mesh);

  };

  virtual void TearDown()
  {
    delete mesh;
  }

  //MeshIoMsh msh_reader;  // MT ESTRANHO, NAO PODE SER AQUI ... ME PARECE BUG DO GTEST.
  //MeshIoVtk vtk_printer;
  Mesh *mesh;
  
};

//const TupleT meshes_t[] = {TupleT(TRIANGLE3,     "meshes/iter_tri3.msh" )};
const TupleT meshes_t[] = {TupleT(TRIANGLE3,     "meshes/iter_tri3.msh" ),
                           TupleT(TRIANGLE6,     "meshes/iter_tri6.msh" ),
                           TupleT(QUADRANGLE4,   "meshes/iter_qua4.msh" ),
                           TupleT(QUADRANGLE8,   "meshes/iter_qua8.msh" ),
                           TupleT(QUADRANGLE9,   "meshes/iter_qua9.msh" ),
                           TupleT(TETRAHEDRON4,  "meshes/iter_tet4.msh" ),
                           TupleT(TETRAHEDRON10, "meshes/iter_tet10.msh"),
                           TupleT(HEXAHEDRON8,   "meshes/iter_hex8.msh" ),
                           TupleT(HEXAHEDRON20,  "meshes/iter_hex20.msh"),
                           TupleT(HEXAHEDRON27,  "meshes/iter_hex27.msh")};


TEST_P(IteratorTest, PointIteratorsTest)
{
  int const n_nodes = mesh->numNodes();
  
  for (int i = 0; i < n_nodes; ++i)
  {
    mesh->getNodePtr(i)->setTag(0);
  }
  
  // ====================== traversing the mesh ====================================
  
  {
    point_iterator point = mesh->pointBegin();
    point_iterator point_end = mesh->pointEnd();
    
    for (; point != point_end; ++point)
      point->setTag(point->getTag()+1);
    // check
    for (int i = 0; i < n_nodes; ++i)
    {
      EXPECT_EQ( 1, mesh->getNodePtr(i)->getTag());
      mesh->getNodePtr(i)->setTag(0); // reseting
    }
  }

  //EXPECT_EQ(1,1);
}



TEST_P(IteratorTest, CellIteratorsTest)
{
  
  for (int i = 0; i < mesh->numCells(); ++i)
  {
    mesh->getCellPtr(i)->setTag(0);
  }
  
  // ====================== traversing the mesh ====================================
  
  {
    cell_iterator cell = mesh->cellBegin();
    cell_iterator cell_end = mesh->cellEnd();
    
    for (; cell != cell_end; ++cell)
      cell->setTag(cell->getTag()+1);
    // check
    for (int i = 0; i < mesh->numCells(); ++i)
    {
      EXPECT_EQ( 1, mesh->getCellPtr(i)->getTag());
      mesh->getCellPtr(i)->setTag(0); // reseting
    }
  }
  
 
  
}


TEST_P(IteratorTest, FacetIteratorsTest)
{
  
  for (int i = 0; i < mesh->numFacets(); ++i)
  {
    mesh->getFacetPtr(i)->setTag(0);
  }
  
  // ====================== traversing the mesh ====================================
  
  {
    facet_iterator facet = mesh->facetBegin();
    facet_iterator facet_end = mesh->facetEnd();
    
    for (; facet != facet_end; ++facet)
      facet->setTag(facet->getTag()+1);
    // check
    for (int i = 0; i < mesh->numFacets(); ++i)
    {
      EXPECT_EQ( 1, mesh->getFacetPtr(i)->getTag());
      mesh->getFacetPtr(i)->setTag(0); // reseting
    }
  }
  
  
  
}


TEST_P(IteratorTest, CornerIteratorsTest)
{
  
  if (mesh->cellDim()>2)
  {
    for (int i = 0; i < mesh->numCorners(); ++i)
    {
      mesh->getCornerPtr(i)->setTag(0);
    }
    
    // ====================== traversing the mesh ====================================
    
    {
      corner_iterator corner = mesh->cornerBegin();
      corner_iterator corner_end = mesh->cornerEnd();
      
      for (; corner != corner_end; ++corner)
        corner->setTag(corner->getTag()+1);
      // check
      for (int i = 0; i < mesh->numCorners(); ++i)
      {
        EXPECT_EQ( 1, mesh->getCornerPtr(i)->getTag());
        mesh->getCornerPtr(i)->setTag(0); // reseting
      }
    }
    
  }
  
  
}


INSTANTIATE_TEST_CASE_P(IteratorTestCase, IteratorTest, ::testing::ValuesIn(meshes_t));



