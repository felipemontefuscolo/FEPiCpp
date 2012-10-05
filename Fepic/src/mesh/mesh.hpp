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

#ifndef FEPIC_MESH_HPP
#define FEPIC_MESH_HPP

#include <vector>
#include <deque>
#include "../util/macros.hpp"
#include "../util/timer.hpp"
#include "elements/cellcore.hpp"
#include "elements/elements.hpp"
#include "../util/list_type.hpp"
#include "mesh_iterators.hpp"
#include "entityhandler.hpp"
#include "boost/ptr_container/ptr_vector.hpp"
#include "boost/ptr_container/ptr_deque.hpp"


class Mesh;

namespace fep_internal
{
  typedef SeqList<std::deque<Point>     , SetVector<int> > Mesh_PointContainer ;
  typedef SeqList<std::deque<Corner>    , SetVector<int> > Mesh_CornerContainer;
  typedef SeqList<std::deque<Facet>     , SetVector<int> > Mesh_FacetContainer ;
  typedef SeqList<boost::ptr_deque<Cell>, SetVector<int> > Mesh_CellContainer  ;
  
  typedef typename Mesh_PointContainer::iterator  Mesh_PointContainer_iterator;
  typedef typename Mesh_CornerContainer::iterator Mesh_CornerContainer_iterator;
  typedef typename Mesh_FacetContainer::iterator  Mesh_FacetContainer_iterator;
  typedef typename Mesh_CellContainer::iterator   Mesh_CellContainer_iterator;

  typedef typename Mesh_PointContainer::const_iterator  Mesh_PointContainer_const_iterator;  
  typedef typename Mesh_CornerContainer::const_iterator Mesh_CornerContainer_const_iterator;
  typedef typename Mesh_FacetContainer::const_iterator  Mesh_FacetContainer_const_iterator;
  typedef typename Mesh_CellContainer::const_iterator   Mesh_CellContainer_const_iterator;

  typedef MeshIterator<Mesh_PointContainer_iterator , Mesh_PointContainer , Mesh>  point_iterator ;
  typedef MeshIterator<Mesh_CornerContainer_iterator, Mesh_CornerContainer, Mesh>  corner_iterator;
  typedef MeshIterator<Mesh_FacetContainer_iterator , Mesh_FacetContainer , Mesh>  facet_iterator ;
  typedef MeshIterator<Mesh_CellContainer_iterator  , Mesh_CellContainer  , Mesh>  cell_iterator  ;
  
  typedef MeshIterator<Mesh_PointContainer_const_iterator , Mesh_PointContainer , Mesh>  point_const_iterator ;
  typedef MeshIterator<Mesh_CornerContainer_const_iterator, Mesh_CornerContainer, Mesh>  corner_const_iterator;
  typedef MeshIterator<Mesh_FacetContainer_const_iterator , Mesh_FacetContainer , Mesh>  facet_const_iterator ;
  typedef MeshIterator<Mesh_CellContainer_const_iterator  , Mesh_CellContainer  , Mesh>  cell_const_iterator  ;  
  
  //typedef MeshIterator<CellElement> abstract_iterator; // TODO
  typedef EntityHandler<Cell  , Mesh>      cell_handler  ;
  typedef EntityHandler<Point , Mesh>      point_handler ;
  typedef EntityHandler<Facet , Mesh>      facet_handler ;
  typedef EntityHandler<Corner, Mesh>      corner_handler;

  typedef EntityHandler<Cell   const, Mesh const> cell_const_handler;
  typedef EntityHandler<Point  const, Mesh const> point_const_handler;
  typedef EntityHandler<Facet  const, Mesh const> facet_const_handler;
  typedef EntityHandler<Corner const, Mesh const> corner_const_handler;
  
}


//
// Public
//

typedef fep_internal::point_iterator  point_iterator ;
typedef fep_internal::corner_iterator corner_iterator;
typedef fep_internal::facet_iterator  facet_iterator ;
typedef fep_internal::cell_iterator   cell_iterator  ;

typedef fep_internal::point_const_iterator  point_const_iterator ;
typedef fep_internal::corner_const_iterator corner_const_iterator;
typedef fep_internal::facet_const_iterator  facet_const_iterator ;
typedef fep_internal::cell_const_iterator   cell_const_iterator  ;

typedef fep_internal::cell_handler   cell_handler  ;
typedef fep_internal::point_handler  point_handler ;
typedef fep_internal::facet_handler  facet_handler ;
typedef fep_internal::corner_handler corner_handler;

typedef fep_internal::cell_const_handler   cell_const_handler  ;
typedef fep_internal::point_const_handler  point_const_handler ;
typedef fep_internal::facet_const_handler  facet_const_handler ;
typedef fep_internal::corner_const_handler corner_const_handler;


class Mesh
{
  // iterators
  template<class,class,class> friend class MeshIterator;
  template<class,class>       friend class EntityHandler;
  
public:
  typedef fep_internal::Mesh_PointContainer  PointList ;
  typedef fep_internal::Mesh_CornerContainer CornerList;
  typedef fep_internal::Mesh_FacetContainer  FacetList ;
  typedef fep_internal::Mesh_CellContainer   CellList  ;

  typedef  fep_internal::Mesh_PointContainer_iterator  CellIteratorT;
  typedef  fep_internal::Mesh_CornerContainer_iterator PointIteratorT;
  typedef  fep_internal::Mesh_FacetContainer_iterator  FacetIteratorT;
  typedef  fep_internal::Mesh_CellContainer_iterator   CornerIteratorT;

  typedef  fep_internal::Mesh_PointContainer_const_iterator  CellConstIteratorT;
  typedef  fep_internal::Mesh_CornerContainer_const_iterator PointConstIteratorT;
  typedef  fep_internal::Mesh_FacetContainer_const_iterator  FacetConstIteratorT;
  typedef  fep_internal::Mesh_CellContainer_const_iterator   CornerConstIteratorT;
  
protected:

  friend class MeshIoMsh;

  // ==========================================================
  // BEGIN
  //                      Mesh Attributes
  //
  //  
  CellList      m_cells_list;
  PointList     m_points_list;
  FacetList     m_facets_list;
  CornerList    m_corners_list;

  std::map<int, int>  m_connected_comp_l; // connected component vs initial cell id list
  std::map<int, int>  m_boundary_comp_l;  // boundary component vs initialfacet id list  

  int m_spacedim;

  // Cell attributes
  bool m_is_parametric_cell;
  int m_cell_dim;
  int m_n_nodes_per_cell;
  int m_n_nodes_per_facet;
  int m_n_nodes_per_corner;
  int m_n_vertices_per_cell;
  int m_n_vertices_per_facet;
  int m_n_vertices_per_corner;
  int m_n_facets_per_cell;
  int m_n_corners_per_cell;
  int m_n_corners_per_facet;
  bool m_cell_has_edge_nodes;
  bool m_cell_has_face_nodes;
  bool m_cell_has_volume_nodes;
  
public:  
  Timer timer; // time measure
protected:
  ECellType m_cellm_fep_tag;
  EMshTag   m_cell_msh_tag;
  bool      m_dont_build_adjacency;  
  
  //  
  //
  //                      Mesh Attributes 
  //  END
  // ==========================================================

  /// constructor
  Mesh(ECellType fept=UNDEFINED_CELLT, int spacedim = -1);
public:
  ~Mesh() {};


public:

  static Mesh* create(ECellType type, int spacedim = -1);

  int cellDim() const
  {
    return this->m_cell_dim;
  }
  ECellType cellType() const
  {
    return this->m_cellm_fep_tag;
  }
  EMshTag cellMshTag() const
  {
    return this->m_cell_msh_tag;
  }

  bool isVertex(CellElement const* p) const
  {
    return p->getPosition() < this->numVerticesPerCell();
  }  
  
  bool inSingleCell(Point const* p) const;
  bool inSingleCell(Corner const* p) const;
  bool inSingleCell(Facet const* p) const;

  bool cellHasEdgeNodes() const
  {
    return this->m_cell_has_edge_nodes;
  }
  bool cellHasFaceNodes() const
  {
    return this->m_cell_has_face_nodes;
  }
  bool cellHasVolumeNodes() const
  {
    return this->m_cell_has_volume_nodes;
  };

  int nodesPerCell() const
  {
    return this->m_n_nodes_per_cell;
  }

  int nodesPerFacet() const
  {
    return this->m_n_nodes_per_facet;
  }

  int nodesPerCorner() const
  {
    return this->m_n_nodes_per_corner;
  }

  int verticesPerCell() const
  {
    return this->m_n_vertices_per_cell;
  }

  int verticesPerFacet() const
  {
    return this->m_n_vertices_per_facet;
  }

  int verticesPerCorner() const
  {
    return this->m_n_vertices_per_corner;
  }

  void printInfo() const;
  void printStatistics() const;

  /** Retorna a n-ésima celula (adivinha o tipo de célula pelo tipo da malha)
  */
  Cell* getCellPtr(int nth)
  {
    if (unsigned(nth)<this->m_cells_list.totalSize())
      return &m_cells_list[nth];
    else
      return NULL;
  }
  /** Retorna a n-ésima facet
  */
  Facet* getFacetPtr(int nth)
  {
    if (unsigned(nth)<this->m_facets_list.totalSize())
      return &m_facets_list[nth];
    else
      return NULL;
  }
  /** Retorna a n-ésima corner
  */
  Corner* getCornerPtr(int nth)
  {
    if (unsigned(nth)<this->m_corners_list.totalSize())
      return &m_corners_list[nth];
    else
      return NULL;
  }
  /** Retorna o n-ésimo nó da malha.
  */
  Point* getNodePtr(int nth)
  {
    if (unsigned(nth)<this->m_points_list.totalSize())
      return &m_points_list[nth];
    else
      return NULL;
  }

  /** Retorna a n-ésima celula (adivinha o tipo de célula pelo tipo da malha)
  */
  Cell const* getCellPtr(int nth) const
  {
    if (unsigned(nth)<this->m_cells_list.totalSize())
      return &m_cells_list[nth];
    else
      return NULL;
  }
  /** Retorna a n-ésima facet
  */
  Facet const* getFacetPtr(int nth) const
  {
    if (unsigned(nth)<this->m_facets_list.totalSize())
      return &m_facets_list[nth];
    else
      return NULL;
  }
  /** Retorna a n-ésima corner
  */
  Corner const* getCornerPtr(int nth) const
  {
    if (unsigned(nth)<this->m_corners_list.totalSize())
      return &m_corners_list[nth];
    else
      return NULL;
  }
  /** Retorna o n-ésimo nó da malha.
  */
  Point const* getNodePtr(int nth) const
  {
    if (unsigned(nth)<this->m_points_list.totalSize())
      return &m_points_list[nth];
    else
      return NULL;
  }

  cell_handler handler(cell_iterator & it)
  { return cell_handler(this, &*it, it.index()); }

  facet_handler handler(facet_iterator & it)
  { return facet_handler(this, &*it, it.index()); }
  
  corner_handler handler(corner_iterator & it)
  { return corner_handler(this, &*it, it.index()); }  

  point_handler handler(point_iterator & it)
  { return point_handler(this, &*it, it.index()); }

  cell_const_handler handler(cell_const_iterator const& it)
  { return cell_const_handler(this, &*it, it.index()); }

  facet_const_handler handler(facet_const_iterator const& it)
  { return facet_const_handler(this, &*it, it.index()); }
  
  corner_const_handler handler(corner_const_iterator const& it)
  { return corner_const_handler(this, &*it, it.index()); }  

  point_const_handler handler(point_const_iterator const& it)
  { return point_const_handler(this, &*it, it.index()); }
  

private:
  template<class EntityType>
  EntityType* entityPtr(int ith);

public:

  cell_handler getCell(int nth)
  {
    if (unsigned(nth)<this->m_cells_list.totalSize())
      return cell_handler(this, &m_cells_list[nth], nth);
    else
      return cell_handler(this, NULL, -1);
  }
  facet_handler getFacet(int nth)
  {
    if (unsigned(nth)<this->m_facets_list.totalSize())
      return facet_handler(this, &m_facets_list[nth], nth);
    else
      return facet_handler(this, NULL, -1);
  }
  corner_handler getCorner(int nth)
  {
    if (unsigned(nth)<this->m_corners_list.totalSize())
      return corner_handler(this, &m_corners_list[nth], nth);
    else
      return corner_handler(this, NULL, -1);
  }
  point_handler getNode(int nth)
  {
    if (unsigned(nth)<this->m_points_list.totalSize())
      return point_handler(this, &m_points_list[nth], nth);
    else
      return point_handler(this, NULL, -1);
  }


  cell_const_handler getCell(int nth) const
  {
    if (unsigned(nth)<this->m_cells_list.totalSize())
      return cell_const_handler(this, &m_cells_list[nth], nth);
    else
      return cell_const_handler(this, NULL, -1);
  }
  facet_const_handler getFacet(int nth) const
  {
    if (unsigned(nth)<this->m_facets_list.totalSize())
      return facet_const_handler(this, &m_facets_list[nth], nth);
    else
      return facet_const_handler(this, NULL, -1);
  }
  corner_const_handler getCorner(int nth) const
  {
    if (unsigned(nth)<this->m_corners_list.totalSize())
      return corner_const_handler(this, &m_corners_list[nth], nth);
    else
      return corner_const_handler(this, NULL, -1);
  }
  point_const_handler getNode(int nth) const
  {
    if (unsigned(nth)<this->m_points_list.totalSize())
      return point_const_handler(this, &m_points_list[nth], nth);
    else
      return point_const_handler(this, NULL, -1);
  }


  void disablePoint(int id)
  { m_points_list.disable(id); }
  
  void disableCorner(int id)
  { m_corners_list.disable(id); }
  
  void disableFacet(int id)
  { m_facets_list.disable(id); }
  
  void disableCell(int id)
  { m_cells_list.disable(id); }


  void getCellNodesId(Cell const* cell, int *result) const
  {
    cell->getNodesId(result);
  }
  void getCellVerticesId(Cell const* cell, int *result) const
  {
    cell->getVerticesId(result);
  }
  void getCellFacetsId(Cell const* cell, int *result) const
  {
    cell->getFacetsId(result);
  }
  void getCellCornersId(Cell const* cell, int *result) const
  {
    cell->getCornersId(result);
  }
  void getFacetNodesId(Facet const* facet, int *result) const
  {
    int icell = facet->getIncidCell();
    int pos   = facet->getPosition();
    this->getCellPtr(icell)->getFacetNodesId(pos, result);
  }
  void getCornerNodesId(CellElement const* corner, int *result) const
  {
    int icell = corner->getIncidCell();
    int pos   = corner->getPosition();
    this->getCellPtr(icell)->getCornerNodesId(pos, result);
  }


  void getCellNodesId(int id, int *result) const
  {
    this->getCellNodesId(this->getCellPtr(id), result);
  }
  void getFacetNodesId(int id, int *result) const
  {
    this->getFacetNodesId(this->getFacetPtr(id), result);
  }
  void getCornerNodesId(int id, int *result) const
  {
    this->getCornerNodesId(this->getCornerPtr(id), result);
  }


  int getCellContigId(int id) const
  {
    return m_cells_list.contiguousId(id);
  }
  int getFacetContigId(int id) const
  {
    return m_facets_list.contiguousId(id);
  }
  int getCornerContigId(int id) const
  {
    return m_corners_list.contiguousId(id);
  }
  int getNodeContigId(int id) const
  {
    return m_points_list.contiguousId(id);
  }

  void getCellsContigId(int* first, int const* last, int* result) const
  {
    m_cells_list.contiguousIds(first, last, result);
  }
  void getFacetsContigId(int* first, int const* last, int* result)const
  {
    m_facets_list.contiguousIds(first, last, result);
  }
  void getCornersContigId(int* first, int const* last, int* result) const
  {
    m_corners_list.contiguousIds(first, last, result);
  }
  void getNodesContigId(int* first, int const* last, int* result) const
  {
    m_points_list.contiguousIds(first, last, result);
  }


  void getCellNodesContigId(Cell const* cell, int* result) const
  {
    this->getCellNodesId(cell, result);
    for (int i = 0; i < this->numNodesPerCell(); ++i)
    {
      *result = m_points_list.contiguousId(*result);
      ++result;
    }
  }
  void getFacetNodesContigId(Facet const* facet, int* result) const
  {
    this->getFacetNodesId(facet, result);
    for (int i = 0; i < this->numNodesPerFacet(); ++i)
    {
      *result = m_points_list.contiguousId(*result);
      ++result;
    }
  }
  void getCornerNodesContigId(CellElement const* corner, int* result) const
  {
    this->getCornerNodesId(corner, result);
    for (int i = 0; i < this->numNodesPerCorner(); ++i)
    {
      *result = m_points_list.contiguousId(*result);
      ++result;
    }
  }
  
  /** Retorna na matriz X as coordenadas dos nós passados em map.
  *  As colunas de X correspondem a dimensão enquanto as linhas
  *  correspondem aos graus de liberdade.
  * @param[in] first an iterator at the beginning of a vector that stores the indices.
  * @param[in] last an iterator at the end of a vector that stores the indices.
  * @param[out] X an iterator at the beginning of a matrix that stores the coordinates.
  * */
  template<class Iterator1, class Iterator2>
  void getNodesCoords(Iterator1 first, Iterator1 last, Iterator2 X)
  {
    int const sdim = this->spaceDim();
    Real const* c;
    while(first != last)
    {
      c = this->getNodePtr(*first++)->getCoord();
      for (int d = 0; d < sdim; ++d)
        *X++ = *c++;
    }
  }

  // TODO: linear only ...
  void getCenterCoord(Cell const* cell, Real* x) const;
  void getCenterCoord(Facet const* facet, Real* x) const;
  void getCenterCoord(Corner const* corner, Real* x) const;

  /** Check if the vertices form a facet of this mesh, if so returns facet's id.
   * @param[in] vtcs vector with the ids of the vertices.
   * @param[out] fid id of the facet that has those vertices. If the vertices form a facet
   *  in a anticyclically manner, then the negative of the facet's id is returned.
   * @return true if a facet were found, false otherwise.
   */ 
  bool getFacetIdFromVertices(int const* vtcs, int &fid);

  /** Check if the vertices form a corner of this mesh, if so returns corner's id.
   * @param[in] vtcs vector with the ids of the vertices.
   * @param[out] fid id of the corner that has those vertices. If the vertices form a corner
   *  in a anticyclically manner, then the negative of the corner's id is returned.
   * @return true if a corner were found, false otherwise.
   */ 
  bool  getCornerIdFromVertices(int const* vtcs, int &rid);

  // ---------------------------------------------------- EDGE STAR ---------------------------------------------------

private:
  // TODO: implementar versão para 1d
  int* vertexStar_1D(int C, int vC, int *iCs, int *viCs) const;

  /** INTERNAL USE ONLY\n
   *  Returns all incident cells of a vertex.
   *  @param[in] C an incident cell of a vertex to start the search (initial cell).
   *  @param[in] vC vertex's local id in the initial cell.
   *  @param[out] iCs vector to put the incident cells.
   *  @param[out] viCs vertex's local ids in the incident cells.
   *  @return a pointer to the element following the end of the sequence iCs.
   *  @note the vectors iCs and viCs must have enough space to store the data.
   *  @note does not support high order nodes, only vertices.
   */
  int* vertexStar_2D(int C, int vC, int *iCs, int *viCs) const;

  /** INTERNAL USE ONLY\n
   *  Returns all incident cells of a vertex.
   *  @param[in] C an incident cell of a vertex to start the search (initial cell).
   *  @param[in] vC vertex's local id in the initial cell.
   *  @param[out] iCs vector to put the incident cells.
   *  @param[out] viCs vertex's local ids in the incident cells.
   *  @return a pointer to the element following the end of the sequence iCs.
   *  @note the vectors iCs and viCs must have enough space to store the data.
   * */
  int* vertexStar_3D(int C, int vC, int *iCs, int *viCs) const;

public:

  //virtual int* edgeStar(Corner const*, int *iCs, int *eiCs) const = 0;
  /** Returns all incident cells of a vertex.
   *  @param[in] C an incident cell of a vertex to start the search (initial cell).
   *  @param[in] vC vertex's local id in the initial cell.
   *  @param[out] iCs vector to put the incident cells.
   *  @param[out] vCs vertex's local ids in the initial cells.
   *  @return a pointer to the element following the end of the sequence iCs.
   *  @note the vectors iCs and viCs must have enough space to store the data.
   *  @note does not support high order nodes, only vertices.
   */
  int* vertexStar(int C, int vC, int *iCs, int *viCs) const
  {
    switch (this->cellDim())
    {
      case 1:
        return vertexStar_1D(C,vC,iCs,viCs);
        break;
      case 2:
        return vertexStar_2D(C,vC,iCs,viCs);
        break;
      case 3:
        return vertexStar_3D(C,vC,iCs,viCs);
        break;
      default:
        return NULL;
    }
  }
  /** Returns all incident cells of a vertex.
   *  @param[in] p a pointer to the vertex.
   *  @param[out] iCs vector to put the incident cells.
   *  @param[out] vCs vertex's local ids in the initial cells.
   *  @return a pointer to the element following the end of the sequence iCs.
   *  @note the vectors iCs and viCs must have enough space to store the data.
   *  @note does not support high order nodes, only vertices.
   */
  int* vertexStar(Point const* p, int *iCs, int *viCs) const
  {
    return this->vertexStar(p->getIncidCell(),
                            p->getPosition(),
                            iCs, viCs);
  }
  
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  /** Returns all incident cells of a node.
   *  @param[in] C an incident cell of a node to start the search (initial cell).
   *  @param[in] nC node's local id in the initial cell.
   *  @param[out] iCs vector to put the incident cells.
   *  @param[out] niCs node's local ids in the initial cells.
   *  @return a pointer to the element following the end of the sequence iCs.
   *  @note the vectors iCs and niCs must have enough space to store the data.
   *  @note this function assumes that each edge, face or volume has only one
   *        node in its interior, i.e., cells of third order are not supported.
   */
  int* nodeStar(int C, int nC, int *iCs, int *niCs) const
  {
    int const nvpc = this->numVerticesPerCell();
    int const nrpc = this->numCornersPerCell();
    int const nfpc = this->numFacetsPerCell();
    
    if (!this->cellHasEdgeNodes() || nC<nvpc)
    {
      return this->vertexStar(C, nC, iCs, niCs);
    }

    else if (!this->cellHasFaceNodes() || nC<nvpc + nrpc)
    {
      // number of nodes within the cell.
      int const eC = nC - nvpc;
      int * iCs_end = this->edgeStar(C, eC, iCs, niCs);
      *niCs++ = nC;
      ++iCs;
      for (; iCs!=iCs_end; ++iCs,++niCs)
        *niCs += nvpc;
      return iCs_end;
    }

    else if (!this->cellHasVolumeNodes() || nC<nvpc + nrpc + nfpc)
    {
      int const fC = nC - nvpc - nrpc;
      int * iCs_end = this->faceStar(C, fC, iCs, niCs);
      *niCs++ = nC;
      ++iCs;
      for (; iCs!=iCs_end; ++iCs,++niCs)
        *niCs += nvpc + nrpc;
      return iCs_end;
    }
    else
    {
      *iCs++ = C;
      *niCs  = nC;
      return iCs;
    }
  }
  
  /** Returns all incident cells of a node.
   *  @param[in] p a pointer to the node.
   *  @param[out] iCs vector to put the incident cells.
   *  @param[out] niCs node's local ids in the initial cells.
   *  @return a pointer to the end of the sequence iCs.
   *  @note the vectors iCs and niCs must have enough space to store the data.
   *  @note this function assumes that each edge, face or volume has only one
   *        node in its interior, i.e., cells of third order are not supported.
   */
  int* nodeStar(Point const* point, int *iCs, int *niCs) const;

private:
  // TODO: implementar uma versão para 1D
  int* edgeStar_1D(int C, int eC, int *iCs, int *eiCs) const;

  /** INTERNA USE ONLY\n
   * Returns all incidents cells of a edge.
   * @param[in] C an incident cell of a edge to start the search (initial cell).
   * @param[in] eC edge's local id in the initial cell.
   * @param[out] iCs vector to put the incident cells.
   * @param[out] eiCs edges's local ids in the incident cells.
   * @return a pointer to the element following the end of the sequence iCs.
   */
  int* edgeStar_2D(int C, int eC, int *iCs, int *eiCs) const;

  /** INTERNA USE ONLY\n
   * Returns all incidents cells of a edge.
   * @param[in] C an incident cell of a edge to start the search (initial cell).
   * @param[in] eC edge's local id in the initial cell.
   * @param[out] iCs vector to put the incident cells.
   * @param[out] eiCs edges's local ids in the incident cells.
   * @return the a pointer to the element following the end of the sequence iCs.
   */
  int* edgeStar_3D(int C, int eC, int *iCs, int *eiCs) const;

public:

  /** Returns all incidents cells of a edge.
   * @param[in] C an incident cell of a edge to start the search (initial cell).
   * @param[in] eC edge's local id in the initial cell.
   * @param[out] iCs vector to put the incident cells.
   * @param[out] eiCs edges's local ids in the incident cells.
   * @return a pointer to the element following the end of the sequence iCs.
   */
  int* edgeStar(int C, int eC, int *iCs, int *eiCs) const
  {
    switch (this->cellDim())
    {
      case 3:
        return this->edgeStar_3D(C, eC, iCs, eiCs);
        break;
      case 2:
        return this->edgeStar_2D(C, eC, iCs, eiCs);
        break;
      case 1:
        return this->edgeStar_1D(C, eC, iCs, eiCs);
        break;
      default:
        return NULL;
    }
  }

  /** Returns all incidents cells of a edge.
   * @param[in] e a pointer to the edge.
   * @param[out] iCs vector to put the incident cells.
   * @param[out] eiCs edges's local ids in the incident cells.
   * @return a pointer to the element following the end of the sequence iCs.
   * @note if an edge is not a Facet, use the Corner version function.
   */
  int* edgeStar(CellElement const* e, int *iCs, int *eiCs) const
  {
    return this->edgeStar(e->getIncidCell(),
                          e->getPosition(),
                          iCs, eiCs);
  }

  // ---------------------------------------------------- FACE STAR ---------------------------------------------------
private:
  int* faceStar_3D(int C, int fC, int *iCs, int *fiCs) const;

  int* faceStar_2D(int C, int fC, int *iCs, int *fiCs) const;

  int* faceStar_1D(int C, int fC, int *iCs, int *fiCs) const;
public:

  int* faceStar(int C, int fC, int *iCs, int *fiCs) const
  {
    switch (this->cellDim())
    {
      case 3:
        return faceStar_3D(C, fC, iCs, fiCs);
        break;
      case 2:
        return faceStar_2D(C, fC, iCs, fiCs);
        break;
      case 1:
        return faceStar_1D(C, fC, iCs, fiCs);
        break;
      default:
        return NULL;
    }
  }

  /** @brief Returns all vertices that are connected to a vertex.
   *  @param[in] p a pointer to the vertex.
   *  @param[out] iVs vector with the connected vertices.
   *  @return a pointer to the element following the end of the sequence iVs.
   */
  int* connectedVtcs(Point const* p, int *iVs) const;

  /** @brief Returns all vertices that are connected to a vertex, as well the incident cells.
   *  @param[in] p a pointer to the vertex.
   *  @param[out] iVs vector with the connected vertices.
   *  @param[out] iCs vector with the incident cells.
   *  @param[out] viCs vector with the p local-ids in each cell.
   *  @return a pointer to the element following the end of the sequence iVs.
   *  @note iVs[k] = getCellPtr(iCd[k])->getNodeId(viCs[k]);
   */
  int* connectedVtcs(Point const* p, int *iVs, int *iCs, int *viCs) const;

  /** @brief Returns all nodes that are connected to a node.
   *  @param[in] p a pointer to the node.
   *  @param[out] iNs vector with the connected nodes.
   *  @return a pointer to the element following the end of the sequence iNs.
   */
  int* connectedNodes(Point const* p, int *iNs) const;

  /** INTERNAL USE ONLY\n
   *  @brief Returns all incident facets to a vertex.
   *  @param[in] p a pointer to the vertex.
   *  @param[out] iFs vector to put the incident facets.
   *  @param[out] viFs vertex's local ids in the incident facets.
   *  @return a pointer to the element following the end of the sequence iFs.
   *  @note the vectors iCs and viCs must have enough space to store the data.
   *  @note does not support high order nodes, only vertices.
   */
  int* incidentFacets_1D(Point const* p, int *iFs, int *viFs) const;

  /** INTERNAL USE ONLY\n
   *  @brief Returns all incident facets to a vertex.
   *  @param[in] p a pointer to the vertex.
   *  @param[out] iFs vector to put the incident facets.
   *  @param[out] viFs vertex's local ids in the incident facets.
   *  @return a pointer to the element following the end of the sequence iFs.
   *  @note the vectors iCs and viCs must have enough space to store the data.
   *  @note does not support high order nodes, only vertices.
   */
  int* incidentFacets_2D(Point const* p, int *iFs, int *viFs) const;

  /** INTERNAL USE ONLY\n
   *  @brief Returns all incident facets to a vertex.
   *  @param[in] p a pointer to the vertex.
   *  @param[out] iFs vector to put the incident facets.
   *  @param[out] viFs vertex's local ids in the incident facets.
   *  @return a pointer to the element following the end of the sequence iFs.
   *  @note the vectors iCs and viCs must have enough space to store the data.
   *  @note does not support high order nodes, only vertices.
   */
  int* incidentFacets_3D(Point const* p, int *iFs, int *viFs) const;

  /** @brief Returns all incident facets to a vertex.
   *  @param[in] p a pointer to the vertex.
   *  @param[out] iFs vector to put the incident facets.
   *  @param[out] viFs vertex's local ids in the incident facets.
   *  @return a pointer to the element following the end of the sequence iFs.
   *  @note the vectors iCs and viCs must have enough space to store the data.
   *  @note does not support high order nodes, only vertices.
   */
  int* incidentFacets(Point const* p, int *iFs, int *viFs) const
  {
    FEPIC_CHECK(this->isVertex(p), "invalid C or vC", std::invalid_argument);
    
    switch (this->cellDim())
    {
      case 1:
        return this->incidentFacets_1D(p, iFs, viFs);
        break;
      case 2:
        return this->incidentFacets_2D(p, iFs, viFs);
        break;
      case 3:
        return this->incidentFacets_3D(p, iFs, viFs);
        break;
      default:
        return NULL;
    }
  }

  /** @brief Returns all incident facets to a node.
   *  @param[in] nodeid id oh the node.
   *  @param[out] iFs vector to put the incident facets.
   *  @param[out] viFs vertex's local ids in the incident facets.
   *  @return a pointer to the element following the end of the sequence iFs.
   *  @note the vectors iCs and viCs must have enough space to store the data.
   *  @note does not support high order nodes, only vertices.
   */
  int* incidentFacets(int nodeid, int *iFs, int *viFs) const
  {
    FEPIC_CHECK(this->isVertex(this->getNodePtr(nodeid)), "invalid C or vC", std::invalid_argument);
    return this->incidentFacets(this->getNodePtr(nodeid), iFs, viFs);
  }
  
  Facet* nextBoundaryFacet_1D(Facet const*f);

  Facet* nextBoundaryFacet_2D(Facet const*f);
  
  Facet* nextBoundaryFacet_3D(Facet const*f);
  
  /** given a boundary facet, return the next boundary facet.
   */ 
  Facet* nextBoundaryFacet(Facet const*f)
  {
    switch (this->cellDim())
    {
      case 2:
        return this->nextBoundaryFacet_2D(f);
        break;
      case 3:
        return this->nextBoundaryFacet_3D(f);
        break;
      case 1:
        return this->nextBoundaryFacet_1D(f);
        break;                
      default:
        return NULL;
    }
  }
  
  void pushIncidCell2Point(Point *pt, int iC, int pos);


  void qBuildAdjacency(bool b)
  {
    m_dont_build_adjacency = b;
  };
  bool qBuildAdjacency()
  {
    return m_dont_build_adjacency;
  };


  void buildCorners_1D();
  void buildCorners_2D();
  void buildCorners_3D();
  void buildCorners()
  {
    switch (this->cellDim())
    {
      case 3:
        this->buildCorners_3D();
        break;
      case 1:
        this->buildCorners_1D();
        break;
      case 2:
        this->buildCorners_2D();
        break;
      default:
        FEPIC_ASSERT(false, "invalid cell dimension", std::runtime_error);
    }
  }

  
  /// Constrói as informações topológicas da malha
  void buildAdjacency()
  {
    this->timer.restart();
    this->buildCellsAdjacency();
    this->timer.elapsed("buildCellsAdjacency()");

    // only for 3D
    this->timer.restart();
    this->buildCorners();
    this->timer.elapsed("buildCorners_Template()");

    // its need to be here, because these func use getCellId() wich needs
    // adjacency information.
    this->timer.restart();
    this->setUpConnectedComponentsId();
    this->setUpBoundaryComponentsId();
    this->timer.elapsed("setUpConnected/BoundaryComponentsId()");

    // its need etUpConnectedComponentsId() function because of singular vertices
    this->timer.restart();
    this->buildNodesAdjacency();
    this->timer.elapsed("buildNodesAdjacency()");

  }
  
  void buildCellsAdjacency();
  void buildNodesAdjacency();
  /// Constrói as informações topológicas da malha

  bool inBoundary(Point const* p) const
  {
    return p->inBoundary();
  }
  bool inBoundary(Facet const* f) const
  {
    Cell const* icell = this->getCellPtr(f->getIncidCell());
    if (icell->getIncidCell(f->getPosition()) < 0)
      return true;
    else
      return false;
  }
  bool inBoundary(Corner const* f) const
  {
    Cell const* icell = this->getCellPtr(f->getIncidCell());
    if (!icell->inBoundary())
      return false;
    if (this->cellDim() == 2)
    {
      Point const* point = this->getNodePtr( icell->getNodeId(  f->getPosition() ) );
      return this->inBoundary(point);
    }
    else
    {
      int const m = f->getPosition();
      if (icell->getIncidCell(icell->table_bC_x_fC(m,0)) < 0) return true;
      if (icell->getIncidCell(icell->table_bC_x_fC(m,1)) < 0) return true;
      return false;
    }
  }


  int pushCell(Cell const* C);
  int pushPoint(Point const* P);
  int pushFacet(Facet const* h);
  int pushCorner(Corner const* b);

  Cell*   pushCell(int *id);
  Point*  pushPoint(int *id);
  Facet*  pushFacet(int *id);
  Corner* pushCorner(int *id);

  Cell*   createCell() const;
  Point*  createPoint() const;
  Facet*  createFacet() const;
  Corner* createCorner() const;

  /** Retorna o número células
  *  @note não conta com o/a(s) marcado/a(s) como killed.
  */
  int numCells() const
  {
    return static_cast<int>( m_cells_list.size() );
  }

  /** Retorna o número de células.
  * @note incluindo o/a(s) marcado/a(s) como killed.
  */
  int numCellsTotal() const
  {
    return static_cast<int>( m_cells_list.totalSize() );
  }

  /** Retorna o número de nós.
  *  @note não conta com o/a(s) marcado/a(s) como killed.
  */
  int numNodes() const
  {
    return static_cast<int>( m_points_list.size() );
  }

  /** Retorna no número de nós.
  *  @note incluindo o/a(s) marcado/a(s) como killed.
  */
  int numNodesTotal() const
  {
    return static_cast<int>( m_points_list.totalSize() );
  }

  int numVertices() const;

  /** Retorna número de facets.
  * @note não conta com o/a(s) marcado/a(s) como killed.
  */
  int numFacets() const
  {
    return static_cast<int>( m_facets_list.size() );
  }

  /** Retorna o número de facets.
  *  @note incluindo o/a(s) marcado/a(s) como killed.
  */
  int numFacetsTotal() const
  {
    return static_cast<int>( m_facets_list.totalSize() );
  }

  /** Retorna número de corners.
  * @note não conta com o/a(s) marcado/a(s) como killed.
  */
  int numCorners() const
  {
    if (this->m_cell_dim == 3)
      return static_cast<int>( m_corners_list.size() );
    else
    // FIXME: high orders nodes are not corners
      return static_cast<int>( m_points_list.totalSize() );
  }

  /** Retorna o número de corners.
  *  @note incluindo o/a(s) marcado/a(s) como killed.
  */
  int numCornersTotal() const
  {
    if (this->m_cell_dim == 3)
      return static_cast<int>( m_corners_list.totalSize() );
    else
      // FIXME: high orders nodes are not corners
      return static_cast<int>( m_points_list.totalSize() );
  }


  int numNodesPerCell() const
  {
    return this->m_n_nodes_per_cell;
  }
  int numNodesPerFacet() const
  {
    return this->m_n_nodes_per_facet;
  }
  int numNodesPerCorner() const
  {
    return this->m_n_nodes_per_corner;
  }

  int numVerticesPerCell() const
  {
    return this->m_n_vertices_per_cell;
  }
  int numVerticesPerFacet() const
  {
    return this->m_n_vertices_per_facet;
  }
  int numVerticesPerCorner() const
  {
    return this->m_n_vertices_per_corner;
  }

  int numFacetsPerCell() const
  {
    return this->m_n_facets_per_cell;
  }
  int numCornersPerCell() const
  {
    return this->m_n_corners_per_cell;
  }
  int numCornersPerFacet() const
  {
    return this->m_n_corners_per_facet;
  }


  static unsigned estimateNumFacets(unsigned nc, ECellType t);
  static unsigned estimateNumCorners(unsigned nc, ECellType t);

  int spaceDim() const
  {
    return m_spacedim;
  }

  int getCellId(Cell const* a) const
  {
    FEPIC_CHECK(m_dont_build_adjacency, "this function can not be called without adjacency", std::runtime_error);
    
    int ic_id = a->getIncidCell(0);
    
    if (ic_id >= 0)
      return this->getCellPtr(ic_id)->getIncidCell(a->getIncidCellPos(0));
    else
      return this->getFacetPtr(a->getFacetId(0))->getIncidCell();
  }
  
  int getPointId(Point const* a) const
  {
    FEPIC_CHECK(m_dont_build_adjacency, "this function can not be called without adjacency", std::runtime_error);
    return this->getCellPtr(a->getIncidCell())->getNodeId(a->getPosition());
  }
  
  int getFacetId(Facet const* a) const
  {
    FEPIC_CHECK(m_dont_build_adjacency, "this function can not be called without adjacency", std::runtime_error);
    return this->getCellPtr(a->getIncidCell())->getFacetId(a->getPosition());
  }
  
  int getCornerId(Corner const* a) const
  {
    FEPIC_CHECK(m_dont_build_adjacency, "this function can not be called without adjacency", std::runtime_error);
    return this->getCellPtr(a->getIncidCell())->getCornerId(a->getPosition());
  }

  // --------------------------------------------------- ADJACENCY -------------------------------------------------------

  /** NOT FOR USERS
   * auxiliary function 
   * set the connected component id starting from c_ini with id cc_id.
   * @warning isVisited() flags are used.
   */ 
  void fi_setConnectedComponentsId(cell_handler c_ini, int cc_id);
  void fi_setBoundaryComponentsId(facet_handler f_ini, int bc_id);

  void setUpConnectedComponentsId();
  int numConnectedComponents() const
  {
    return static_cast<int>( m_connected_comp_l.size() );
  }
  
  void setUpBoundaryComponentsId();
  int numBoundaryComponents() const
  {
    return static_cast<int>( m_boundary_comp_l.size() );
  }

  void getConnectedComponentsPicks(int *comps, int *cells) const
  {
    int k = 0;
    for (std::map<int,int>::const_iterator it = m_connected_comp_l.begin(); it != m_connected_comp_l.end(); ++it)
    {
      comps[k] = (*it).first;
      cells[k] = (*it).second;
      ++k;
    }
    
  }
  void getBoundaryComponentsPicks(int *comps, int *facets) const
  {
    int k = 0;
    for (std::map<int,int>::const_iterator it = m_boundary_comp_l.begin(); it != m_boundary_comp_l.end(); ++it)
    {
      comps[k] = (*it).first;
      facets[k] = (*it).second;
      ++k;
    }
    
  }

protected:
  Cell*   incEnabledCell(int &id);
  Point*  incEnabledPoint(int &id);
  Facet*  incEnabledFacet(int &id);
  Corner* incEnabledCorner(int &id);

  Cell*   decEnabledCell(int &id);
  Point*  decEnabledPoint(int &id);
  Facet*  decEnabledFacet(int &id);
  Corner* decEnabledCorner(int &id);

public:

  cell_iterator   cellBegin();
  cell_iterator   cellEnd();
  point_iterator  pointBegin();
  point_iterator  pointEnd();
  facet_iterator  facetBegin();
  facet_iterator  facetEnd();
  corner_iterator cornerBegin();
  corner_iterator cornerEnd();

  cell_iterator   cellBegin(int tid, int nthreads);
  cell_iterator   cellEnd(int tid, int nthreads);
  point_iterator  pointBegin(int tid, int nthreads);
  point_iterator  pointEnd(int tid, int nthreads);
  facet_iterator  facetBegin(int tid, int nthreads);
  facet_iterator  facetEnd(int tid, int nthreads);
  corner_iterator cornerBegin(int tid, int nthreads);
  corner_iterator cornerEnd(int tid, int nthreads);


};




// =====================================================================================
// =====================================================================================
//                                  ITERATORS
// =====================================================================================
// =====================================================================================

cell_iterator   FEP_STRONG_INLINE Mesh::cellBegin()  { return cell_iterator  (this, m_cells_list.begin()  );}
cell_iterator   FEP_STRONG_INLINE Mesh::cellEnd()    { return cell_iterator  (this, m_cells_list.end()    );}
point_iterator  FEP_STRONG_INLINE Mesh::pointBegin() { return point_iterator (this, m_points_list.begin() );}
point_iterator  FEP_STRONG_INLINE Mesh::pointEnd()   { return point_iterator (this, m_points_list.end()   );}
facet_iterator  FEP_STRONG_INLINE Mesh::facetBegin() { return facet_iterator (this, m_facets_list.begin() );}
facet_iterator  FEP_STRONG_INLINE Mesh::facetEnd()   { return facet_iterator (this, m_facets_list.end()   );}
corner_iterator FEP_STRONG_INLINE Mesh::cornerBegin(){ return corner_iterator(this, m_corners_list.begin());}
corner_iterator FEP_STRONG_INLINE Mesh::cornerEnd()  { return corner_iterator(this, m_corners_list.end()  );}

cell_iterator   FEP_STRONG_INLINE Mesh::cellBegin  (int tid, int nthreads) { return cell_iterator  (this, m_cells_list.begin  (tid, nthreads)  );}
cell_iterator   FEP_STRONG_INLINE Mesh::cellEnd    (int tid, int nthreads) { return cell_iterator  (this, m_cells_list.end    (tid, nthreads)  );}
point_iterator  FEP_STRONG_INLINE Mesh::pointBegin (int tid, int nthreads) { return point_iterator (this, m_points_list.begin (tid, nthreads)  );}
point_iterator  FEP_STRONG_INLINE Mesh::pointEnd   (int tid, int nthreads) { return point_iterator (this, m_points_list.end   (tid, nthreads)  );}
facet_iterator  FEP_STRONG_INLINE Mesh::facetBegin (int tid, int nthreads) { return facet_iterator (this, m_facets_list.begin (tid, nthreads)  );}
facet_iterator  FEP_STRONG_INLINE Mesh::facetEnd   (int tid, int nthreads) { return facet_iterator (this, m_facets_list.end   (tid, nthreads)  );}
corner_iterator FEP_STRONG_INLINE Mesh::cornerBegin(int tid, int nthreads) { return corner_iterator(this, m_corners_list.begin(tid, nthreads)  );}
corner_iterator FEP_STRONG_INLINE Mesh::cornerEnd  (int tid, int nthreads) { return corner_iterator(this, m_corners_list.end  (tid, nthreads)  );}

// =====================================================================================
// =====================================================================================


template<class CT>
struct ArgConverter // CT = SeqList value type
{
  template<class U> // U = push function argument
  CT const& operator() (U const* u)
  {
    return *static_cast<Cell const*>(u);
  }
};

template<class CT>
struct ArgConverter<CT*>
{
  template<class U>
  CT* operator() (U const* u)
  {
    return new CT( *static_cast<Cell const*>(u) );
  }
};


/** Add a cell in the mesh's list.
 *  @param cell a pointer to the cell.
 *  @return the id of the new cell.
*/
FEP_STRONG_INLINE
int Mesh::pushCell(Cell const* cell)
{ return m_cells_list.insert(cell->clone()); }

/** Add a node in the mesh's list.
 *  @param node a pointer to the node.
 *  @return the id of the new node.
*/
FEP_STRONG_INLINE
int Mesh::pushPoint(Point const* node)
{ return m_points_list.insert(*node); }

/** Add a facet in the mesh's list.
 *  @param facet a pointer to the facet.
 *  @return the id of the new facet.
*/
FEP_STRONG_INLINE
int Mesh::pushFacet(Facet const* facet)
{ return m_facets_list.insert(*facet); }

/** Add a corner in the mesh's list.
 *  @param corner a pointer to the corner.
 *  @return the id of the new corner.
*/
FEP_STRONG_INLINE
int Mesh::pushCorner(Corner const* corner)
{ return m_corners_list.insert(*corner); }


/** Create a cell in the mesh's list.
 *  @param[out] cell_id a pointer to store the id of the new cell, can
 *              be (int*)NULL pointer.
 *  @return a pointer to the new cell.
*/
FEP_STRONG_INLINE
Cell* Mesh::pushCell(int *cell_id)
{
  int const tmp = m_cells_list.insert(Cell::create(this->cellType()));
  if (cell_id)
    *cell_id = tmp;
  return this->getCellPtr(tmp);
}

/** Create a node in the mesh's list.
 *  @param[out] node_id a pointer to store the id of the new node, can
 *              be (int*)NULL pointer.
 *  @return a pointer to the new node.
*/
FEP_STRONG_INLINE
Point* Mesh::pushPoint(int *node_id)
{
  int const tmp = m_points_list.insert(Point());
  if (node_id)
    *node_id = tmp;
  return this->getNodePtr(tmp);
}

/** Create a facet in the mesh's list.
 *  @param[out] facet_id a pointer to store the id of the new facet, can
 *              be (int*)NULL pointer.
 *  @return a pointer to the new facet.
*/
FEP_STRONG_INLINE
Facet* Mesh::pushFacet(int *facet_id)
{
  int const tmp = m_facets_list.insert(Facet());
  if (facet_id)
    *facet_id = tmp;
  return this->getFacetPtr(tmp);
}

/** Create a corner in the mesh's list.
 *  @param[out] corner_id a pointer to store the id of the new corner, can
 *              be (int*)NULL pointer.
 *  @return a pointer to the new corner.
*/
FEP_STRONG_INLINE
Corner* Mesh::pushCorner(int *corner_id)
{
  int const tmp = m_corners_list.insert(Corner());
  if (corner_id)
    *corner_id = tmp;  
  return this->getCornerPtr(tmp);
}


FEP_STRONG_INLINE
Cell* Mesh::createCell() const
{ return Cell::create(this->cellType()); }

FEP_STRONG_INLINE
Point* Mesh::createPoint() const
{ return Point::create(); }

FEP_STRONG_INLINE
Facet* Mesh::createFacet() const
{ return Facet::create(); }

FEP_STRONG_INLINE
Corner* Mesh::createCorner() const
{ return Corner::create(); }



#endif


