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

#ifndef FEPIC_CELLCORE_HPP
#define FEPIC_CELLCORE_HPP

#include "cell_base.hpp"

#include "Fepic/src/mesh/enums.hpp"
#include "Fepic/src/util/assert.hpp"
#include "Fepic/src/util/forward_declarations.hpp"
#include "Fepic/src/util/typedefs.hpp"


// 1D Cells members
#define FEP_DEF_1D_CELLS_MEMBERS                \
                                                   \
  int m_icells[n_facets]; /* incident cells id*/  \
  union                                          \
  {                                               \
    int m_nodes[n_nodes]; /* nodes id; N=order */   \
    int m_facets[n_vertices]; /* alias to m_nodes */ \
    int m_corners[n_vertices]; /* dummy */                \
  };                                               \
  union                                        \
  {                                             \
    char m_icells_pos[n_facets];                 \
    char m_icells_anchors[n_facets];             \
  }

// 2D Cells members
#define FEP_DEF_2D_CELLS_MEMBERS                           \
                                                    \
  int m_facets[n_facets];  /* facets id  */                  \
  int m_icells[n_facets]; /* incident cells id */            \
  union                                                     \
  {                                                          \
    int m_nodes[n_nodes];   /* nodes id */                    \
    int m_corners[n_vertices]; /* alias for m_nodes */         \
  };                                                         \
  union                                                    \
  {                                                         \
    char m_icells_pos[n_facets]; /* positions of icells */   \
    char m_icells_anchors[n_facets]; /* dummy */             \
  }                                                         \
  
  
// 3D Cells members
#define FEP_DEF_3D_CELLS_MEMBERS                                    \
  char m_icells_pos[n_facets];  /* positions of icells */            \
  char m_icells_anchors[n_facets];/* anchors of icells*/ \
  int  m_facets[n_facets];   /* facets id */                         \
  int  m_icells[n_facets];    /* incident cells id */                \
  int  m_corners[n_corners];  /* edges id */                         \
  /* nodes id */                                                     \
  int  m_nodes[n_nodes]


template<typename CellT>
class iCellCore  : public Cell
{
#if !defined(THIS) && !defined(CONST_THIS)
  #define THIS static_cast<CellT*>(this)
  #define CONST_THIS static_cast<const CellT*>(this)
#endif

protected:
  iCellCore()
  {
    for (int i = 0; i < CellT::n_facets; ++i)
      THIS->m_icells[i] = -1;
    if (CellT::dim > 1)
      for (int i = 0; i < CellT::n_facets; ++i)
        THIS->m_facets[i] = -1;
    else
    {
      THIS->m_nodes[0] = -1;
      THIS->m_nodes[1] = -1;
    }
    if (CellT::dim==3)
      for (int i = 0; i < CellT::n_corners; ++i)
        THIS->m_corners[i] = -1;
    for (int i = 0; i < CellT::n_nodes; ++i)
      THIS->m_nodes[i] = -1;
  };
  //iCellCore(iCellCore const&) {};

public:

  void copy(Cell const& c)
  {
    (*THIS) = *static_cast<const CellT*>(&c);
  }

  Cell* clone() const
  {
    return new CellT(*CONST_THIS);
  }

  int dim() const
  {
    return CellT::dim;
  }

  EMshTag getMshTag() const
  {
    return CellT::msh_tag;
  }

  int getCornerId(int corner) const
  {
    if (CellT::dim==3)
      return CONST_THIS->m_corners[corner];
    else
      return CONST_THIS->m_nodes[corner];
  }

  /** Get nodes of a corner
   * @param[in] f the corner local id.
   * @param[out] corner_nds where to put the nodes.
   * @note the vector corner_nds must have enough space allocated (num nodes per corner).
   */
  void getCornerNodesId(int f, int *corner_nds) const;

  /** Get vertices of a corner
   * @param[in] f the corner local id.
   * @param[out] vtcs where to put the vertices.
   * @note the vector vtcs must have enough space allocated (num vertices per corner).
   */
  void getCornerVerticesId(int f, int *corner_vtcs) const;

  /** Returns the dimension of the entity where a node lives.
   * @param local_id local id of the node.
   * @return dimension of the entity where node lives.
   */
  static int getDimFromNodeLives(int local_id)
  {
    if (local_id < CellT::n_vertices)
      return 0;
    else if (local_id < CellT::n_vertices + CellT::n_corners)
      return 1;
    else if (local_id < CellT::n_vertices + CellT::n_corners + CellT::n_facets)
      return 2;
    else
      return 3;
  }

  /** @brief Returns the id of a node of an edge.
   *  @param eC which edge of the cell.
   *  @param nCe which node of the edge.
   */
  int getEdgeNodeId(int eC, int nCe) const
  {
    if (CellT::dim==1)
    {
      return CONST_THIS->m_nodes[nCe];
      ++eC; // avoid compiler warnings
    }
    else if (CellT::dim==2)
    {
      return CONST_THIS->m_nodes[CellT::m_table_fC_x_nC[eC][nCe]];
    }
    else
      return CONST_THIS->m_nodes[CellT::m_table_bC_x_nC[eC][nCe]];
  }

  int getFacetId(int facet) const
  {
    if (CellT::dim > 1)
      return CONST_THIS->m_facets[facet];
    else
      return CONST_THIS->m_nodes[facet];
  }

  /** Get nodes of a facet
   * @param[in] f the facet local id.
   * @param[out] facet_nds where to put the nodes.
   * @note the vector facet_nds must have enough space allocated (num nodes per facet).
   */
  void getFacetNodesId(int f, int * facet_nds) const;

  /** Get vertices of a facet
   * @param[in] f the facet local id.
   * @param[out] vtcs where to put the vertices.
   * @note the vector vtcs must have enough space allocated (num vertices per facet).
   */
  void getFacetVerticesId(int f, int * facet_vtcs) const;

  void getFacetCornersId(int f, int * corns) const;

  int getIncidCell(int facet) const
  {
    return CONST_THIS->m_icells[facet];
  }

  char getIncidCellAnch(int facet) const
  {
    if (CellT::dim==3)
      return CONST_THIS->m_icells_anchors[facet];
    else
      return -1;
  }

  char getIncidCellPos(int facet) const
  {
    return CONST_THIS->m_icells_pos[facet];
  }

  int getNodeId(int const ith) const
  {
    return CONST_THIS->m_nodes[ith];
  }

  void getNodesId(int * begin) const
  {
    std::copy(CONST_THIS->m_nodes, // from
              CONST_THIS->m_nodes+CellT::n_nodes,
              begin);              // to
  }

  void getVerticesId(int * begin) const
  {
    std::copy(CONST_THIS->m_nodes, // from
              CONST_THIS->m_nodes+CellT::n_vertices,
              begin);              // to
  }

  void getFacetsId(int * begin) const
  {
    if (CellT::dim < 2)
    {
      begin[0] = CONST_THIS->m_nodes[0];
      begin[1] = CONST_THIS->m_nodes[1];
    }
    else
      std::copy(CONST_THIS->m_facets, // from
                CONST_THIS->m_facets+CellT::n_facets,
                begin);              // to
  }

  void getCornersId(int * begin) const
  {
    if (CellT::dim == 3)
      std::copy(CONST_THIS->m_corners, // from
                CONST_THIS->m_corners+CellT::n_corners,
                begin);              // to
    else
      std::copy(CONST_THIS->m_nodes, // from
                CONST_THIS->m_nodes+CellT::n_vertices,
                begin);
  }

  bool hasEdgeNodes() const
  {
    return CellT::has_edge_nodes;
  }
  bool hasFaceNodes() const
  {
    return CellT::has_face_nodes;
  }
  bool hasVolumeNodes() const
  {
    return CellT::has_volume_nodes;
  }

  /** Check if the vertices form a corner of this cell, if so returns corner's id.
   * @param[in] vtcs vector with the ids of the vertices.
   * @param[out] lrid local id of the corner that has those vertices. If the vertices form a corner
   *  in a anticyclically manner, then the negative of the corner's id is returned.
   * @return true if a corner were found, false otherwise.
   */
  bool isCorner(int const* vtcs, int &lrid) const;

  bool inBoundary() const;

  /** Check if the vertices form a facet of this cell, if so returns facet's id.
   * @param[in] vtcs vector with the ids of the vertices.
   * @param[out] lfid local id of the facet that has those vertices. If the vertices form a facet
   *  in a anticyclically manner, then the negative of the facet's id is returned.
   * @return true if a facet were found, false otherwise.
   */
  bool isFacet(int const* vtcs, int &lfid) const;

  bool isParametric() const
  {
    return CellT::n_nodes > CellT::n_vertices;
  }

  int  numCorners() const
  {
    return CellT::n_corners;
  }

  int  numCornersPerFacet() const
  {
    return CellT::Derived::n_facets;
  }

  int numFacets() const
  {
    return CellT::n_facets;
  }

  int numNodes() const
  {
    return CellT::n_nodes;
  }

  int numNodesPerCorner() const
  {
    return CellT::n_nodes_per_corner;
  }

  int numNodesPerFacet() const
  {
    return CellT::n_nodes_per_facet;
  }

  int  numVertices() const
  {
    return CellT::n_vertices;
  }

  int numVerticesPerFacet() const
  {
    return CellT::n_vertices_per_facet;
  }

  int  numVerticesPerCorner() const
  {
    return CellT::n_vertices_per_corner;
  }

  void reset()
  {
    resetIncidCells();
    resetFacets();
    resetCorners();
    resetNodes();
  }

  void resetFacets()
  {
    if (CellT::dim > 1)
    {
      for (int i = 0; i < CellT::n_facets; ++i)
      {
        THIS->m_facets[i] = -1;
      }
    }
    else
    {
      THIS->m_nodes[0] = -1;
      THIS->m_nodes[1] = -1;
    }
  }

  void resetCorners()
  {
    if (CellT::dim==3)
      for (int i = 0; i < CellT::n_corners; ++i)
        THIS->m_corners[i] = -1;
    //else if (CellT::dim==2)
    //  for (int i = 0; i < CellT::n_vertices; ++i)
    //    THIS->m_nodes[i] = -1;
  }

  void resetIncidCells()
  {
    for (int i = 0; i < CellT::n_facets; ++i)
    {
      THIS->m_icells[i] = -1;
    }
  }

  void resetNodes()
  {
    for (int i = 0; i < CellT::n_nodes; ++i)
    {
      THIS->m_nodes[i] = -1;
    }    
  }

  void setCornerId(int corner, int cornerid)
  {
    if (CellT::dim ==3)
      THIS->m_corners[corner] = cornerid;

  }

  void setFacetId(int facet, int facetid)
  {
    if (CellT::dim > 1)
      THIS->m_facets[facet] = facetid;
    else
      THIS->m_nodes[facet] = facetid;
  }

  void setIncidCell(int facet, int icellid)
  {
    THIS->m_icells[facet] = icellid;
  }

  void setIncidCellAnch(int facet, int anch)
  {
    if (CellT::dim == 3)
      THIS->m_icells_anchors[facet] = static_cast<char>(  anch  );
  }

  void setIncidCellPos(int facet, int pos)
  {
    THIS->m_icells_pos[facet] = static_cast<char>(  pos );
  }

  /** is same to:
   *  setIncidCell(facet, icellid);
   *  setIncidCellPos(facet, pos);
   *  setIncidCellAnch(facet, anch);
   */
  void setIncidence(int facet, int icellid, char pos, char anch=0)
  {
    THIS->m_icells[facet] = icellid;
    THIS->m_icells_pos[facet] = pos;
    if (CellT::dim == 3)
      THIS->m_icells_anchors[facet] = anch;
  }

  void setNodeId(int const ith, int const nodeid)
  {
    THIS->m_nodes[ith] = nodeid;
  }


  /*  Set all members of the cell. Pass NULL to ignore an members.
   *  
   * 
   */ 
  void setAllMembers(int const* nodes, int const* corners, int const* facets, int const* icells, int const* icells_pos,
                      int const* icells_ancs, int const* conn_comp_id, int const* tag, int const* flags)
  {
    // nodes
    if (nodes != NULL)
      for (int i = 0; i < CellT::n_nodes; ++i)
        THIS->m_nodes[i] = nodes[i];
    
    if (CellT::dim>2 && corners!= NULL)
      for (int i = 0; i < CellT::n_corners; ++i)
        THIS->m_corners[i] = corners[i];
    
    if (CellT::dim>1 && facets!= NULL)
      for (int i = 0; i < CellT::n_facets; ++i)
        THIS->m_facets[i] = facets[i];
    
    if (icells!=NULL)
      for (int i = 0; i < CellT::n_facets; ++i)
        THIS->m_icells[i] = icells[i];
    
    if (icells_pos!=NULL)
      for (int i = 0; i < CellT::n_facets; ++i)
        THIS->m_icells_pos[i] = icells_pos[i];
    
    if (CellT::dim > 2 && icells_ancs!=NULL)
      for (int i = 0; i < CellT::n_facets; ++i)
        THIS->m_icells_anchors[i] = icells_ancs[i];
    
    if (conn_comp_id != NULL)
      THIS->setConnectedComponentId(*conn_comp_id);
      
    if (tag != NULL)
      THIS->setTag(*tag);
    
    if (flags != NULL)
      THIS->setFlags(*flags);
    
  }


  int table_fC_x_vC(int i, int j) const {return CellT::m_table_fC_x_vC[i][j];}
  int table_fC_x_nC(int i, int j) const {return CellT::m_table_fC_x_nC[i][j];}
  int table_vC_x_fC(int i, int j) const {return CellT::m_table_vC_x_fC[i][j];}
  int table_fC_x_bC(int i, int j) const {return CellT::m_table_fC_x_bC[i][j];}
  int table_bC_x_vC(int i, int j) const {return CellT::m_table_bC_x_vC[i][j];}
  int table_bC_x_nC(int i, int j) const {return CellT::m_table_bC_x_nC[i][j];}
  int table_bC_x_fC(int i, int j) const {return CellT::m_table_bC_x_fC[i][j];}

  ~iCellCore() {};


#undef THIS
#undef CONST_THIS

};




#endif




