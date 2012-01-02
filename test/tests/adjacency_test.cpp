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


class AdjacencyTest : public testing::TestWithParam<TupleT> { // vertex id
public:

  virtual void SetUp() {
    
    ECellType cell_t     = std::tr1::get<0>(GetParam());
    const char* mesh_in  = std::tr1::get<1>(GetParam());
    
    mesh = Mesh::create(cell_t);
    
    msh_reader.readFileMsh(mesh_in, mesh);
    

  };

  virtual void TearDown()
  {
    delete mesh;
  }

  MeshIoMsh msh_reader;
  MeshIoVtk vtk_printer;
  Mesh *mesh;
};

const TupleT ct_in_2d3d[] =  {TupleT(TRIANGLE3,     "meshes/simptri3.msh" ),
                            TupleT(TRIANGLE6,     "meshes/simptri6.msh" ),
                            TupleT(QUADRANGLE4,   "meshes/simpquad4.msh" ),
                            TupleT(QUADRANGLE8,   "meshes/simpquad8.msh" ),
                            TupleT(QUADRANGLE9,   "meshes/simpquad9.msh" ),
                            TupleT(TETRAHEDRON4,  "meshes/simptet4.msh" ),
                            TupleT(TETRAHEDRON10, "meshes/simptet10.msh"),
                            TupleT(HEXAHEDRON8,   "meshes/simphex8.msh" ),
                            TupleT(HEXAHEDRON20,  "meshes/simphex20.msh"),
                            TupleT(HEXAHEDRON27,  "meshes/simphex27.msh")};

const TupleT ct_in_3d[] =  {TupleT(TETRAHEDRON4,  "meshes/simptet4.msh" ),
                            TupleT(TETRAHEDRON10, "meshes/simptet10.msh"),
                            TupleT(HEXAHEDRON8,   "meshes/simphex8.msh" ),
                            TupleT(HEXAHEDRON20,  "meshes/simphex20.msh"),
                            TupleT(HEXAHEDRON27,  "meshes/simphex27.msh")};

TEST_P(AdjacencyTest, FacetTest) {

  //const int n_cells = mesh->numCells();
  //const int n_nodes = mesh->numNodes();

  Cell const* icell;
  int oth;
  
  cell_iterator cell = mesh->cellBegin();
  cell_iterator cell_end = mesh->cellEnd();
  int cell_id =  0;
  int icell_id;
  for (; cell != cell_end ; ++cell)
  {
    for (int i = 0; i < mesh->numFacetsPerCell(); ++i)
    {
      icell_id = cell->getIncidCell(i);
      if (icell_id < 0)
        break;
      icell = mesh->getCell(icell_id);
      oth = cell->getIncidCellPos(i);
      
      EXPECT_EQ(cell->getFacetId(i), icell->getFacetId(oth))
        << "cell id ="<< cell_id<<"; icell id="<<icell_id ;
    }
  }

}

TEST_P(AdjacencyTest, RemoveTest) {

  //const int n_cells = mesh->numCells();
  //const int n_nodes = mesh->numNodes();

  //if (mesh->cellType() != TRIANGLE3)
  //  return;

  cell_iterator cell = mesh->cellBegin();
  cell_iterator cell_end = mesh->cellEnd();

  for (; cell != cell_end ; ++cell)
  {
    for (int i = 0; i < mesh->numFacetsPerCell(); ++i)
    {
      int facet_id = cell->getFacetId(i);
      EXPECT_TRUE(facet_id >= 0);
      EXPECT_FALSE(mesh->getFacet(facet_id)->disabled());
    }
  }
  
  MeshTools::removeCell(mesh->getCell(1), mesh);
  MeshTools::removeCell(mesh->getCell(2), mesh);
  MeshTools::removeCell(mesh->getCell(3), mesh);
  
  cell = mesh->cellBegin();
  cell_end = mesh->cellEnd();
  
  for (; cell != cell_end ; ++cell)
  {
    for (int i = 0; i < mesh->numFacetsPerCell(); ++i)
    {
      int facet_id = cell->getFacetId(i);
      EXPECT_TRUE(facet_id >= 0);
      EXPECT_FALSE(mesh->getFacet(facet_id)->disabled())
        << mesh->getCellId(&*cell);
    }
  }

  
  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  
  for (; point != point_end; ++point)
  {
    EXPECT_TRUE(point->getIncidCell() >= 0);
    EXPECT_TRUE(point->getPosition() >= 0);
    cell = cell_iterator(mesh,mesh->getCell(point->getIncidCell()));
    EXPECT_EQ(mesh->getPointId(&*point), cell->getNodeId(point->getPosition()));
  }
  
  facet_iterator facet = mesh->facetBegin();
  facet_iterator facet_end = mesh->facetEnd();
  
  for (; facet != facet_end; ++facet)
  {
    EXPECT_TRUE(facet->getIncidCell() >= 0);
    EXPECT_TRUE(facet->getPosition() >= 0);
    cell = cell_iterator(mesh,mesh->getCell(facet->getIncidCell()));
    EXPECT_EQ(mesh->getFacetId(&*facet), cell->getFacetId(facet->getPosition()));
  }
  
  corner_iterator corner = mesh->cornerBegin();
  corner_iterator corner_end = mesh->cornerEnd();
  
  for (; corner != corner_end; ++corner)
  {
    EXPECT_TRUE(corner->getIncidCell() >= 0);
    EXPECT_TRUE(corner->getPosition() >= 0);
    cell = cell_iterator(mesh,mesh->getCell(corner->getIncidCell()));
    EXPECT_EQ(mesh->getCornerId(&*corner), cell->getCornerId(corner->getPosition()));
  }  
}


INSTANTIATE_TEST_CASE_P(FacetCase, AdjacencyTest, ::testing::ValuesIn(ct_in_2d3d));


