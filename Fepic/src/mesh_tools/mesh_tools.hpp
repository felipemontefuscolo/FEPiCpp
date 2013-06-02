// This file is part of FEPiC++, a toolbox for finite element codes.
//
// FEPiC++ is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// FEPiC++ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// FEPiC++. If not, see <http://www.gnu.org/licenses/>.


#ifndef FEPIC_MESH_TOOLS_HPP
#define FEPIC_MESH_TOOLS_HPP

#include "../mesh/mesh.hpp"
#include <utility>
#include <vector>

class Mesh;
class Cell;

class MeshTools
{
public:
  
  //friend class Mesh;
  
  static int checkConsistency(Mesh *mesh);
  static void removeCell(Cell * cell, Mesh *mesh);  
  static void removeTriCell(Cell * cell, Mesh *mesh);  
  static void readMesh(int n_nodes, int n_cells, int const* nodes, Real const* xyz, Mesh *mesh);
  
};

/// specific for triangular meshes
class MeshToolsTri
{
public:
  // =====================================================================================================
  // Métodos de adaptação ================================================================================
  // =====================================================================================================
  static bool flipEdge(Cell * cell, int fid, Mesh *mesh, bool move_edge_nds = true);
  static bool flipEdge(Facet* edge, Mesh* mesh, bool move_edge_nds = true)
  {
    return flipEdge(mesh->getCellPtr(edge->getIncidCell()), edge->getPosition(), mesh, move_edge_nds);
  }    
  static int insertVertexOnEdge(int cell_A_id, int fidA, Real t, Mesh *mesh);
  static int insertVertexOnEdge(Facet * edge, Real t, Mesh *mesh)
  {
    return insertVertexOnEdge(edge->getIncidCell(), edge->getPosition(), t, mesh);
  }   
  static int collapseEdge2d(int cell, int fidA, Real t, Mesh *mesh);
  static int collapseEdge2d(Facet const* edge, Real t, Mesh *mesh)
  {
    return collapseEdge2d(edge->getIncidCell(), edge->getPosition(), t,mesh);
  }
  
  // =====================================================================================================
  // Métodos de testes geométricos =======================================================================
  // =====================================================================================================
  static bool inCircle2d(Cell const* cell, const int fid, Mesh const* mesh);
  static bool inCircle2d(Facet const* edge, Mesh const* mesh)
  {
     return inCircle2d(mesh->getCellPtr(edge->getIncidCell()), edge->getPosition(),mesh);
  }   
  /** Search a point int the mesh.
   * @param x the point coordinates.
   * @param c0 a first cell to start the search.
   * @param mesh the mesh.
   * @return a "pair" type, where pair.first is a bool that is true if
   * the point was found, and pair.second is a pointer to the cell where
   * the point is. If the point was not found, pair.second returns closest
   * cell.
   * @note the mesh has to be convex, otherwise the result is not guaranteed.
   */ 
  static std::pair<bool, Cell *> searchConvexPoint(Real const* x, Cell const* c0, Mesh const* mesh);
  static Real area(Cell const* cell, Mesh const* mesh); 
  static void centroid(Cell const* cell, Real *coord_ctrd, Mesh const* mesh);
  static Real edgeSize(Cell const* cell, const int fid, Mesh const* mesh);  
  static Real meanRatio(Cell const* cell, Mesh const* mesh);
  
};

/// specific for tetrahedral meshes
class MeshToolsTetr
{
public:
  // =====================================================================================================
  // Métodos de adaptação ================================================================================
  // =====================================================================================================
   static int insertVertexOnEdge3d(int cell_A_id, int face_Am_id, Real t, Mesh *mesh);
  // static int colapseEdge3d(int cell_A_id, int corner_Am_pos, Real t, Mesh *mesh); // não passou por um teste
  
  // =====================================================================================================
  // Métodos de testes geométricos =======================================================================
  // =====================================================================================================
  static Real volume(Cell const* cell, Mesh const* mesh); 
  
  // =====================================================================================================
  // Métodos de uso para as adaptações ===================================================================
  // =====================================================================================================
  static int calcAnchors(int *idsA, int *idsB, Mesh *mesh);
  static int calcAnchors(int cellA_id, int facetA_id, int cellB_id, int facetB_id, Mesh *mesh);
  static int getOpositFacetPos(int vtx_pos);
  static void getOpositVertexsPos(int vtx_t_pos, int vtx_b_pos, int* vtx_r_pos, int* vtx_l_pos);
  static void getOpositCornersPos(int vtx_t_pos, int vtx_b_pos, int* corner_tr_pos, int* corner_tl_pos,
																int* corner_br_pos, int* corner_bl_pos, 
																int* corner_op_pos);
																  
  // =====================================================================================================
  // Método de teste 3D ==================================================================================
  // =====================================================================================================
  static bool checkMesh(Mesh *mesh);
};


#endif // FEPIC_IMESH_IMPL_HPP
