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
#include "boost/ptr_container/ptr_vector.hpp"
#include "boost/ptr_container/ptr_deque.hpp"


typedef _MeshIterator<Cell>        cell_iterator;
typedef _MeshIterator<Point>       point_iterator;
typedef _MeshIterator<Facet>       facet_iterator;
typedef _MeshIterator<Corner>      corner_iterator;
typedef _MeshIterator<CellElement> abstract_iterator; // TODO

typedef _MeshIterator<Cell>        const_cell_iterator;
typedef _MeshIterator<Point>       const_point_iterator;
typedef _MeshIterator<Facet>       const_facet_iterator;
typedef _MeshIterator<Corner>      const_corner_iterator;
typedef _MeshIterator<CellElement> const_abstract_iterator; // TODO

// TODO: criar uma classe para cell_handler, etc.
// atribuir funções especiais para os handlers
typedef BaseHandler<Cell>    cell_handler;
typedef BaseHandler<Point>   point_handler;
typedef BaseHandler<Facet>   facet_handler;
typedef BaseHandler<Corner>  corner_handler;


class Mesh
{
public:
  typedef SeqList<boost::ptr_deque<Cell>, SetVector<int> >  CellList;
  typedef SeqList<std::deque<Point>,     SetVector<int> >  PointList;
  typedef SeqList<std::deque<Facet>,     SetVector<int> >  FacetList;
  typedef SeqList<std::deque<Corner>,    SetVector<int> >  CornerList;

  typedef typename CellList  ::iterator CellIteratorT;
  typedef typename PointList ::iterator PointIteratorT;
  typedef typename FacetList ::iterator FacetIteratorT;
  typedef typename CornerList::iterator CornerIteratorT;

  typedef typename CellList  ::const_iterator CellConstIteratorT;
  typedef typename PointList ::const_iterator PointConstIteratorT;
  typedef typename FacetList ::const_iterator FacetConstIteratorT;
  typedef typename CornerList::const_iterator CornerConstIteratorT;
protected:

  // iterators
  template<class> friend class _MeshIterator;
  friend class MeshIoMsh;

  // ==========================================================
  // BEGIN
  //                      Mesh Attributes
  //
  //  
  CellList      _cellL;
  PointList     _pointL;
  FacetList     _facetL;
  CornerList    _cornerL;

  std::map<int, int>  _connected_compL; // connected component vs initial cell id list
  std::map<int, int>  _boundary_compL;  // boundary component vs initialfacet id list  

  int _spacedim;

  // Cell attributes
  bool _is_parametric_cell;
  int _cell_dim;
  int _n_nodes_per_cell;
  int _n_nodes_per_facet;
  int _n_nodes_per_corner;
  int _n_vertices_per_cell;
  int _n_vertices_per_facet;
  int _n_vertices_per_corner;
  int _n_facets_per_cell;
  int _n_corners_per_cell;
  int _n_corners_per_facet;
  bool _cell_has_edge_nodes;
  bool _cell_has_face_nodes;
  bool _cell_has_volume_nodes;
  
public:  
  Timer timer; // time measure
protected:
  ECellType _cell_fep_tag;
  EMshTag   _cell_msh_tag;
  bool      _dont_build_adjacency;  
  
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
    return this->_cell_dim;
  }
  ECellType cellType() const
  {
    return this->_cell_fep_tag;
  }
  EMshTag cellMshTag() const
  {
    return this->_cell_msh_tag;
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
    return this->_cell_has_edge_nodes;
  }
  bool cellHasFaceNodes() const
  {
    return this->_cell_has_face_nodes;
  }
  bool cellHasVolumeNodes() const
  {
    return this->_cell_has_volume_nodes;
  };

  int nodesPerCell() const
  {
    return this->_n_nodes_per_cell;
  }

  int nodesPerFacet() const
  {
    return this->_n_nodes_per_facet;
  }

  int nodesPerCorner() const
  {
    return this->_n_nodes_per_corner;
  }

  int verticesPerCell() const
  {
    return this->_n_vertices_per_cell;
  }

  int verticesPerFacet() const
  {
    return this->_n_vertices_per_facet;
  }

  int verticesPerCorner() const
  {
    return this->_n_vertices_per_corner;
  }

  void printInfo() const;
  void printStatistics() const;

  /** Retorna a n-ésima celula (adivinha o tipo de célula pelo tipo da malha)
  */
  Cell* getCellPtr(int nth)
  {
    if (unsigned(nth)<this->_cellL.totalSize())
      return &_cellL[nth];
    else
      return NULL;
  }
  /** Retorna a n-ésima facet
  */
  Facet* getFacetPtr(int nth)
  {
    if (unsigned(nth)<this->_facetL.totalSize())
      return &_facetL[nth];
    else
      return NULL;
  }
  /** Retorna a n-ésima corner
  */
  Corner* getCornerPtr(int nth)
  {
    if (unsigned(nth)<this->_cornerL.totalSize())
      return &_cornerL[nth];
    else
      return NULL;
  }
  /** Retorna o n-ésimo nó da malha.
  */
  Point* getNodePtr(int nth)
  {
    if (unsigned(nth)<this->_pointL.totalSize())
      return &_pointL[nth];
    else
      return NULL;
  }

  /** Retorna a n-ésima celula (adivinha o tipo de célula pelo tipo da malha)
  */
  Cell const* getCellPtr(int nth) const
  {
    if (unsigned(nth)<this->_cellL.totalSize())
      return &_cellL[nth];
    else
      return NULL;
  }
  /** Retorna a n-ésima facet
  */
  Facet const* getFacetPtr(int nth) const
  {
    if (unsigned(nth)<this->_facetL.totalSize())
      return &_facetL[nth];
    else
      return NULL;
  }
  /** Retorna a n-ésima corner
  */
  Corner const* getCornerPtr(int nth) const
  {
    if (unsigned(nth)<this->_cornerL.totalSize())
      return &_cornerL[nth];
    else
      return NULL;
  }
  /** Retorna o n-ésimo nó da malha.
  */
  Point const* getNodePtr(int nth) const
  {
    if (unsigned(nth)<this->_pointL.totalSize())
      return &_pointL[nth];
    else
      return NULL;
  }


  cell_handler getCell(int nth)
  {
    if (unsigned(nth)<this->_cellL.totalSize())
      return cell_handler(this, &_cellL[nth], nth);
    else
      return cell_handler(this, NULL, -1);
  }
  facet_handler getFacet(int nth)
  {
    if (unsigned(nth)<this->_facetL.totalSize())
      return facet_handler(this, &_facetL[nth], nth);
    else
      return facet_handler(this, NULL, -1);
  }
  corner_handler getCorner(int nth)
  {
    if (unsigned(nth)<this->_cornerL.totalSize())
      return corner_handler(this, &_cornerL[nth], nth);
    else
      return corner_handler(this, NULL, -1);
  }
  point_handler getNode(int nth)
  {
    if (unsigned(nth)<this->_pointL.totalSize())
      return point_handler(this, &_pointL[nth], nth);
    else
      return point_handler(this, NULL, -1);
  }

  void disablePoint(int id)
  {
    _pointL.disable(id);
  }
  void disableCorner(int id)
  {
    _cornerL.disable(id);
  }
  void disableFacet(int id)
  {
    _facetL.disable(id);
  }
  void disableCell(int id)
  {
    _cellL.disable(id);
  }


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
    return _cellL.contiguousId(id);
  }
  int getFacetContigId(int id) const
  {
    return _facetL.contiguousId(id);
  }
  int getCornerContigId(int id) const
  {
    return _cornerL.contiguousId(id);
  }
  int getNodeContigId(int id) const
  {
    return _pointL.contiguousId(id);
  }

  void getCellsContigId(int* first, int const* last, int* result) const
  {
    _cellL.contiguousIds(first, last, result);
  }
  void getFacetsContigId(int* first, int const* last, int* result)const
  {
    _facetL.contiguousIds(first, last, result);
  }
  void getCornersContigId(int* first, int const* last, int* result) const
  {
    _cornerL.contiguousIds(first, last, result);
  }
  void getNodesContigId(int* first, int const* last, int* result) const
  {
    _pointL.contiguousIds(first, last, result);
  }


  void getCellNodesContigId(Cell const* cell, int* result) const
  {
    this->getCellNodesId(cell, result);
    for (int i = 0; i < this->numNodesPerCell(); ++i)
    {
      *result = _pointL.contiguousId(*result);
      ++result;
    }
  }
  void getFacetNodesContigId(Facet const* facet, int* result) const
  {
    this->getFacetNodesId(facet, result);
    for (int i = 0; i < this->numNodesPerFacet(); ++i)
    {
      *result = _pointL.contiguousId(*result);
      ++result;
    }
  }
  void getCornerNodesContigId(CellElement const* corner, int* result) const
  {
    this->getCornerNodesId(corner, result);
    for (int i = 0; i < this->numNodesPerCorner(); ++i)
    {
      *result = _pointL.contiguousId(*result);
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
    _dont_build_adjacency = b;
  };
  bool qBuildAdjacency()
  {
    return _dont_build_adjacency;
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
    return static_cast<int>( _cellL.size() );
  }

  /** Retorna o número de células.
  * @note incluindo o/a(s) marcado/a(s) como killed.
  */
  int numCellsTotal() const
  {
    return static_cast<int>( _cellL.totalSize() );
  }

  /** Retorna o número de nós.
  *  @note não conta com o/a(s) marcado/a(s) como killed.
  */
  int numNodes() const
  {
    return static_cast<int>( _pointL.size() );
  }

  /** Retorna no número de nós.
  *  @note incluindo o/a(s) marcado/a(s) como killed.
  */
  int numNodesTotal() const
  {
    return static_cast<int>( _pointL.totalSize() );
  }

  int numVertices() const;

  /** Retorna número de facets.
  * @note não conta com o/a(s) marcado/a(s) como killed.
  */
  int numFacets() const
  {
    return static_cast<int>( _facetL.size() );
  }

  /** Retorna o número de facets.
  *  @note incluindo o/a(s) marcado/a(s) como killed.
  */
  int numFacetsTotal() const
  {
    return static_cast<int>( _facetL.totalSize() );
  }

  /** Retorna número de corners.
  * @note não conta com o/a(s) marcado/a(s) como killed.
  */
  int numCorners() const
  {
    if (this->_cell_dim == 3)
      return static_cast<int>( _cornerL.size() );
    else
    // FIXME: high orders nodes are not corners
      return static_cast<int>( _pointL.totalSize() );
  }

  /** Retorna o número de corners.
  *  @note incluindo o/a(s) marcado/a(s) como killed.
  */
  int numCornersTotal() const
  {
    if (this->_cell_dim == 3)
      return static_cast<int>( _cornerL.totalSize() );
    else
      // FIXME: high orders nodes are not corners
      return static_cast<int>( _pointL.totalSize() );
  }


  int numNodesPerCell() const
  {
    return this->_n_nodes_per_cell;
  }
  int numNodesPerFacet() const
  {
    return this->_n_nodes_per_facet;
  }
  int numNodesPerCorner() const
  {
    return this->_n_nodes_per_corner;
  }

  int numVerticesPerCell() const
  {
    return this->_n_vertices_per_cell;
  }
  int numVerticesPerFacet() const
  {
    return this->_n_vertices_per_facet;
  }
  int numVerticesPerCorner() const
  {
    return this->_n_vertices_per_corner;
  }

  int numFacetsPerCell() const
  {
    return this->_n_facets_per_cell;
  }
  int numCornersPerCell() const
  {
    return this->_n_corners_per_cell;
  }
  int numCornersPerFacet() const
  {
    return this->_n_corners_per_facet;
  }


  static unsigned estimateNumFacets(unsigned nc, ECellType t);
  static unsigned estimateNumCorners(unsigned nc, ECellType t);

  int spaceDim() const
  {
    return _spacedim;
  }

  int getCellId(Cell const* a) const
  {
    FEPIC_CHECK(_dont_build_adjacency, "this function can not be called without adjacency", std::runtime_error);
    
    int ic_id = a->getIncidCell(0);
    
    if (ic_id >= 0)
      return this->getCellPtr(ic_id)->getIncidCell(a->getIncidCellPos(0));
    else
      return this->getFacetPtr(a->getFacetId(0))->getIncidCell();
  }
  
  int getPointId(Point const* a) const
  {
    FEPIC_CHECK(_dont_build_adjacency, "this function can not be called without adjacency", std::runtime_error);
    return this->getCellPtr(a->getIncidCell())->getNodeId(a->getPosition());
  }
  
  int getFacetId(Facet const* a) const
  {
    FEPIC_CHECK(_dont_build_adjacency, "this function can not be called without adjacency", std::runtime_error);
    return this->getCellPtr(a->getIncidCell())->getFacetId(a->getPosition());
  }
  
  int getCornerId(Corner const* a) const
  {
    FEPIC_CHECK(_dont_build_adjacency, "this function can not be called without adjacency", std::runtime_error);
    return this->getCellPtr(a->getIncidCell())->getCornerId(a->getPosition());
  }

  // --------------------------------------------------- ADJACENCY -------------------------------------------------------

  /** NOT FOR USERS
   * auxiliary function 
   * set the connected component id starting from c_ini with id cc_id.
   * @warning isVisited() flags are used.
   */ 
  void _setConnectedComponentsId(cell_handler c_ini, int cc_id);
  void _setBoundaryComponentsId(facet_handler f_ini, int bc_id);

  void setUpConnectedComponentsId();
  int numConnectedComponents() const
  {
    return static_cast<int>( _connected_compL.size() );
  }
  
  void setUpBoundaryComponentsId();
  int numBoundaryComponents() const
  {
    return static_cast<int>( _boundary_compL.size() );
  }

  void getConnectedComponentsPicks(int *comps, int *cells) const
  {
    int k = 0;
    for (std::map<int,int>::const_iterator it = _connected_compL.begin(); it != _connected_compL.end(); ++it)
    {
      comps[k] = (*it).first;
      cells[k] = (*it).second;
      ++k;
    }
    
  }
  void getBoundaryComponentsPicks(int *comps, int *facets) const
  {
    int k = 0;
    for (std::map<int,int>::const_iterator it = _boundary_compL.begin(); it != _boundary_compL.end(); ++it)
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




#endif


