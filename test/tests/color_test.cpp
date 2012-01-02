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

typedef std::tr1::tuple<ECellType, const char*, const char*> TupleT;


class ColoringTest : public testing::TestWithParam<TupleT> { // vertex id
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

const TupleT ct_in_out[] = {//TupleT(TRIANGLE3,     "meshes/cortri3.msh" , "meshes/out/cortri3.vtk" ),
                            //TupleT(TRIANGLE6,     "meshes/cortri6.msh" , "meshes/out/cortri6.vtk" ),
                            //TupleT(QUADRANGLE4,   "meshes/corqua4.msh" , "meshes/out/corqua4.vtk" ),
                            //TupleT(QUADRANGLE8,   "meshes/corqua8.msh" , "meshes/out/corqua8.vtk" ),
                            //TupleT(QUADRANGLE9,   "meshes/corqua9.msh" , "meshes/out/corqua9.vtk" ),
                            //TupleT(TETRAHEDRON4,  "meshes/cortet4.msh" , "meshes/out/cortet4.vtk" ),
                            //TupleT(TETRAHEDRON10, "meshes/cortet10.msh", "meshes/out/cortet10.vtk"),
                            //TupleT(HEXAHEDRON8,   "meshes/corhex8.msh" , "meshes/out/corhex8.vtk" ),
                            //TupleT(HEXAHEDRON20,  "meshes/corhex20.msh", "meshes/out/corhex20.vtk"),
                            TupleT(HEXAHEDRON27,  "meshes/corhex27.msh", "meshes/out/corhex27.vtk")};
                            

TEST_P(ColoringTest, ResultTest) {

  const int n_cells = mesh->numCells();
  const int n_nodes = mesh->numNodes();
  
  std::vector<int> c_colors(n_cells);
  std::vector<int> n_colors(n_nodes);
  
  for (int j = 0; j < n_cells; ++j)
    c_colors[j] = static_cast<int>( mesh->getCell(j)->getColor() );
    
  for (int j = 0; j < n_nodes; ++j)
    n_colors[j] = static_cast<int>( mesh->getNode(j)->getColor() );
  
  // checking
  const int n_c_colors = mesh->numCellColors();
  const int n_n_colors = mesh->numNodeColors();
  const int nodes_p_c= mesh->nodesPerCell();
  Cell const*cell1;
  Cell const*cell2;
  int nodes1[nodes_p_c];
  int nodes2[nodes_p_c];
  
  EXPECT_TRUE(n_c_colors > 0);
  EXPECT_TRUE(n_n_colors > 0);
  
  for (int i = 0; i < n_c_colors; ++i)
  {
    for (int c = 0; c < n_cells; ++c)
    {
      cell1 = mesh->getCell(c);
      if (cell1->getColor() != (char)i)
        continue;
      cell1->getNodesId(nodes1);
      for (int k = 0; k < n_cells; ++k)
      {
        cell2 = mesh->getCell(k);
        if (cell2->getColor() != (char)i || (c==k))
          continue;
        cell2->getNodesId(nodes2);
        
        EXPECT_TRUE( nodes1+nodes_p_c == std::find_first_of(nodes1,nodes1+nodes_p_c,nodes2,nodes2+nodes_p_c) )
          << "Cell1 = " << c << ", Cell2 = " << k;
      }
    }
  }
  

  // printing
  const char* mesh_out = std::tr1::get<2>(GetParam());
  vtk_printer.attachMesh(mesh);
  vtk_printer.writeVtk(mesh_out);
  vtk_printer.addCellIntVtk("cell_color", DefaultGetDataVtk(NULL,c_colors.data()));
  vtk_printer.addNodeIntVtk("node_color", DefaultGetDataVtk(NULL,n_colors.data()));

}
INSTANTIATE_TEST_CASE_P(ColoringCase, ColoringTest, ::testing::ValuesIn(ct_in_out));

