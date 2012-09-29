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
  int _icells[n_facets]; /* incident cells id*/  \
  union                                          \
  {                                               \
    int _nodes[n_nodes]; /* nodes id; N=order */   \
    int _facets[n_vertices]; /* alias to _nodes */ \
    int _corners[n_vertices]; /* dummy */                \
  };                                               \
  union                                        \
  {                                             \
    char _icells_pos[n_facets];                 \
    char _icells_anchors[n_facets];             \
  }

// 2D Cells members
#define FEP_DEF_2D_CELLS_MEMBERS                           \
                                                    \
  int _facets[n_facets];  /* facets id  */                  \
  int _icells[n_facets]; /* incident cells id */            \
  union                                                     \
  {                                                          \
    int _nodes[n_nodes];   /* nodes id */                    \
    int _corners[n_vertices]; /* alias for _nodes */         \
  };                                                         \
  union                                                    \
  {                                                         \
    char _icells_pos[n_facets]; /* positions of icells */   \
    char _icells_anchors[n_facets]; /* dummy */             \
  }                                                         \
  
  
// 3D Cells members
#define FEP_DEF_3D_CELLS_MEMBERS                                    \
  char _icells_pos[n_facets];  /* positions of icells */            \
  char _icells_anchors[n_vertices_per_facet];/* anchors of icells*/ \
  int  _facets[n_facets];   /* facets id */                         \
  int  _icells[n_facets];    /* incident cells id */                \
  int  _corners[n_corners];  /* edges id */                         \
  /* nodes id */                                                     \
  int  _nodes[n_nodes]


template<typename CellT>
class _CellCore  : public Cell
{
#if !defined(THIS) && !defined(CONST_THIS)
  #define THIS static_cast<CellT*>(this)
  #define CONST_THIS static_cast<const CellT*>(this)
#endif

protected:
  _CellCore()
  {
    for (int i = 0; i < CellT::n_facets; ++i)
      THIS->_icells[i] = -1;
    if (CellT::dim > 1)
      for (int i = 0; i < CellT::n_facets; ++i)
        THIS->_facets[i] = -1;
    else
    {
      THIS->_nodes[0] = -1;
      THIS->_nodes[1] = -1;
    }
    if (CellT::dim==3)
      for (int i = 0; i < CellT::n_corners; ++i)
        THIS->_corners[i] = -1;
    for (int i = 0; i < CellT::n_nodes; ++i)
      THIS->_nodes[i] = -1;
  };
  //_CellCore(_CellCore const&) {};

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
      return CONST_THIS->_corners[corner];
    else
      return CONST_THIS->_nodes[corner];
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
      return CONST_THIS->_nodes[nCe];
      ++eC; // avoid compiler warnings
    }
    else if (CellT::dim==2)
    {
      return CONST_THIS->_nodes[CellT::_table_fC_x_nC[eC][nCe]];
    }
    else
      return CONST_THIS->_nodes[CellT::_table_bC_x_nC[eC][nCe]];
  }

  int getFacetId(int facet) const
  {
    if (CellT::dim > 1)
      return CONST_THIS->_facets[facet];
    else
      return CONST_THIS->_nodes[facet];
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
    return CONST_THIS->_icells[facet];
  }

  char getIncidCellAnch(int facet) const
  {
    if (CellT::dim==3)
      return CONST_THIS->_icells_anchors[facet];
    else
      return -1;
  }

  char getIncidCellPos(int facet) const
  {
    return CONST_THIS->_icells_pos[facet];
  }

  int getNodeId(int const ith) const
  {
    return CONST_THIS->_nodes[ith];
  }

  void getNodesId(int * begin) const
  {
    std::copy(CONST_THIS->_nodes, // from
              CONST_THIS->_nodes+CellT::n_nodes,
              begin);              // to
  }

  void getVerticesId(int * begin) const
  {
    std::copy(CONST_THIS->_nodes, // from
              CONST_THIS->_nodes+CellT::n_vertices,
              begin);              // to
  }

  void getFacetsId(int * begin) const
  {
    if (CellT::dim < 2)
    {
      begin[0] = CONST_THIS->_nodes[0];
      begin[1] = CONST_THIS->_nodes[1];
    }
    else
      std::copy(CONST_THIS->_facets, // from
                CONST_THIS->_facets+CellT::n_facets,
                begin);              // to
  }

  void getCornersId(int * begin) const
  {
    if (CellT::dim == 3)
      std::copy(CONST_THIS->_corners, // from
                CONST_THIS->_corners+CellT::n_corners,
                begin);              // to
    else
      std::copy(CONST_THIS->_nodes, // from
                CONST_THIS->_nodes+CellT::n_vertices,
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
        THIS->_facets[i] = -1;
      }
    }
    else
    {
      THIS->_nodes[0] = -1;
      THIS->_nodes[1] = -1;
    }
  }

  void resetCorners()
  {
    if (CellT::dim==3)
      for (int i = 0; i < CellT::n_corners; ++i)
        THIS->_corners[i] = -1;
    //else if (CellT::dim==2)
    //  for (int i = 0; i < CellT::n_vertices; ++i)
    //    THIS->_nodes[i] = -1;
  }

  void resetIncidCells()
  {
    for (int i = 0; i < CellT::n_facets; ++i)
    {
      THIS->_icells[i] = -1;
    }
  }

  void resetNodes()
  {
    for (int i = 0; i < CellT::n_nodes; ++i)
    {
      THIS->_nodes[i] = -1;
    }    
  }

  void setCornerId(int corner, int cornerid)
  {
    if (CellT::dim ==3)
      THIS->_corners[corner] = cornerid;

  }

  void setFacetId(int facet, int facetid)
  {
    if (CellT::dim > 1)
      THIS->_facets[facet] = facetid;
    else
      THIS->_nodes[facet] = facetid;
  }

  void setIncidCell(int facet, int icellid)
  {
    THIS->_icells[facet] = icellid;
  }

  void setIncidCellAnch(int facet, int anch)
  {
    if (CellT::dim == 3)
      THIS->_icells_anchors[facet] = static_cast<char>(  anch  );
  }

  void setIncidCellPos(int facet, int pos)
  {
    THIS->_icells_pos[facet] = static_cast<char>(  pos );
  }

  /** is same to:
   *  setIncidCell(facet, icellid);
   *  setIncidCellPos(facet, pos);
   *  setIncidCellAnch(facet, anch);
   */
  void setIncidence(int facet, int icellid, char pos, char anch=0)
  {
    THIS->_icells[facet] = icellid;
    THIS->_icells_pos[facet] = pos;
    if (CellT::dim == 3)
      THIS->_icells_anchors[facet] = anch;
  }

  void setNodeId(int const ith, int const nodeid)
  {
    THIS->_nodes[ith] = nodeid;
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
        THIS->_nodes[i] = nodes[i];
    
    if (CellT::dim>2 && corners!= NULL)
      for (int i = 0; i < CellT::n_corners; ++i)
        THIS->_corners[i] = corners[i];
    
    if (CellT::dim>1 && facets!= NULL)
      for (int i = 0; i < CellT::n_facets; ++i)
        THIS->_facets[i] = facets[i];
    
    if (icells!=NULL)
      for (int i = 0; i < CellT::n_facets; ++i)
        THIS->_icells[i] = icells[i];
    
    if (icells_pos!=NULL)
      for (int i = 0; i < CellT::n_facets; ++i)
        THIS->_icells_pos[i] = icells_pos[i];
    
    if (CellT::dim > 2 && icells_ancs!=NULL)
      for (int i = 0; i < CellT::n_facets; ++i)
        THIS->_icells_anchors[i] = icells_ancs[i];
    
    if (conn_comp_id != NULL)
      THIS->setConnectedComponentId(*conn_comp_id);
      
    if (tag != NULL)
      THIS->setTag(*tag);
    
    if (flags != NULL)
      THIS->setFlags(*flags);
    
  }


  int table_fC_x_vC(int i, int j) const {return CellT::_table_fC_x_vC[i][j];}
  int table_fC_x_nC(int i, int j) const {return CellT::_table_fC_x_nC[i][j];}
  int table_vC_x_fC(int i, int j) const {return CellT::_table_vC_x_fC[i][j];}
  int table_fC_x_bC(int i, int j) const {return CellT::_table_fC_x_bC[i][j];}
  int table_bC_x_vC(int i, int j) const {return CellT::_table_bC_x_vC[i][j];}
  int table_bC_x_nC(int i, int j) const {return CellT::_table_bC_x_nC[i][j];}
  int table_bC_x_fC(int i, int j) const {return CellT::_table_bC_x_fC[i][j];}

  ~_CellCore() {};


#undef THIS
#undef CONST_THIS

};




#endif




