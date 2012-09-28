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

template<class CT> class SMesh;

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

  //
  // Mesh Attributes
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
  //
  // End Mesh Attributes
  //

  /// constructor
  Mesh(ECellType fept=UNDEFINED_CELLT, int spacedim = -1);

  virtual Cell*   incEnabledCell(int &id) = 0;
  virtual Point*  incEnabledPoint(int &id) = 0;
  virtual Facet*  incEnabledFacet(int &id) = 0;
  virtual Corner* incEnabledCorner(int &id) = 0;

  virtual Cell*   decEnabledCell(int &id) = 0;
  virtual Point*  decEnabledPoint(int &id) = 0;
  virtual Facet*  decEnabledFacet(int &id) = 0;
  virtual Corner* decEnabledCorner(int &id) = 0;

public:

  // time measure
  Timer timer;

//protected:
  ECellType _cell_fep_tag;
  EMshTag   _cell_msh_tag;
  bool      _dont_build_adjacency;

  virtual cell_iterator   cellBegin() = 0;
  virtual cell_iterator   cellEnd() = 0;
  virtual point_iterator  pointBegin() = 0;
  virtual point_iterator  pointEnd() = 0;
  virtual facet_iterator  facetBegin() = 0;
  virtual facet_iterator  facetEnd() = 0;
  virtual corner_iterator cornerBegin() = 0;
  virtual corner_iterator cornerEnd() = 0;

  virtual cell_iterator   cellBegin(int tid, int nthreads) = 0;
  virtual cell_iterator   cellEnd(int tid, int nthreads) = 0;
  virtual point_iterator  pointBegin(int tid, int nthreads) = 0;
  virtual point_iterator  pointEnd(int tid, int nthreads) = 0;
  virtual facet_iterator  facetBegin(int tid, int nthreads) = 0;
  virtual facet_iterator  facetEnd(int tid, int nthreads) = 0;
  virtual corner_iterator cornerBegin(int tid, int nthreads) = 0;
  virtual corner_iterator cornerEnd(int tid, int nthreads) = 0;

  static Mesh* create(ECellType type, int spacedim = -1);

  virtual int cellDim() const = 0;
  virtual ECellType cellType() const = 0;
  virtual EMshTag cellMshTag() const = 0;

  virtual bool isVertex(CellElement const* p) const = 0;
  virtual bool inSingleCell(Point const* p) const = 0;
  virtual bool inSingleCell(Corner const* p) const = 0;
  virtual bool inSingleCell(Facet const* p) const = 0;

  virtual int nodesPerCell() const = 0;
  virtual int nodesPerFacet() const = 0;
  virtual int nodesPerCorner() const = 0;

  virtual int verticesPerCell() const = 0;
  virtual int verticesPerFacet() const = 0;
  virtual int verticesPerCorner() const = 0;

  virtual void printInfo() const = 0;
  virtual void printStatistics() const = 0;

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

  Cell const* getCellPtr(int nth) const
  {
    FEPIC_CHECK(unsigned(nth)<this->_cellL.totalSize(), "invalid index", std::out_of_range);
    return &_cellL[nth];
  }
  Facet const* getFacetPtr(int nth) const
  {
    FEPIC_CHECK(unsigned(nth)<this->_facetL.totalSize(), "invalid index", std::out_of_range);
    return &_facetL[nth];
  }
  Corner const* getCornerPtr(int nth) const
  {
    FEPIC_CHECK(unsigned(nth)<this->_cornerL.totalSize(), "invalid index", std::out_of_range);
    return &_cornerL[nth];
  }
  Point const* getNodePtr(int nth) const
  {
    FEPIC_CHECK(unsigned(nth)<this->_pointL.totalSize(), "invalid index", std::out_of_range);
    return &_pointL[nth];
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

  virtual bool getFacetIdFromVertices(int const* vtcs, int &fid) =0;
  virtual bool getCornerIdFromVertices(int const* vtcs, int &fid) =0;

  virtual int* edgeStar(int C, int eC, int *iCs, int *eiCs) const = 0;
  virtual int* edgeStar(CellElement const*, int *iCs, int *eiCs) const = 0;
  //virtual int* edgeStar(Corner const*, int *iCs, int *eiCs) const = 0;
  virtual int* vertexStar(int C, int vC, int *iCs, int *viCs) const = 0;
  virtual int* vertexStar(Point const* point, int *iCs, int *viCs) const = 0;
  virtual int* nodeStar(int C, int nC, int *iCs, int *niCs) const = 0;
  virtual int* nodeStar(Point const* point, int *iCs, int *niCs) const = 0;
  virtual int* connectedVtcs(Point const* p, int *iVs) const = 0;
  virtual int* connectedVtcs(Point const* p, int *iVs, int *iCs, int *viCs) const = 0;
  virtual int* connectedNodes(Point const* p, int *iVs) const = 0;
  virtual int* incidentFacets(Point const* p, int *iFs, int *viFs) const = 0;
  virtual int* incidentFacets(int nodeid, int *iFs, int *viFs) const = 0;
  virtual Facet* nextBoundaryFacet(Facet const*f) const = 0;
  virtual void pushIncidCell2Point(Point *pt, int iC, int pos) = 0;


  void qBuildAdjacency(bool b)
  {
    _dont_build_adjacency = b;
  };
  bool qBuildAdjacency()
  {
    return _dont_build_adjacency;
  };
  virtual void buildAdjacency() = 0;
  virtual void buildNodesAdjacency() = 0;

  virtual void setUpConnectedComponentsId() = 0;
  virtual int  numConnectedComponents() const = 0;
  virtual void setUpBoundaryComponentsId() = 0;
  virtual int  numBoundaryComponents() const = 0;
  
  virtual void getConnectedComponentsPicks(int *comps, int *cells) const = 0;
  virtual void getBoundaryComponentsPicks(int *comps, int *facets) const = 0;

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

  virtual Cell*   createCell() const = 0;
  virtual Point*  createPoint() const = 0;
  virtual Facet*  createFacet() const = 0;
  virtual Corner* createCorner() const = 0;

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


  virtual ~Mesh() = 0;

};


/* ====================================================

          ____  __  __ _____ ____  _   _
         / ___||  \/  | ____/ ___|| | | |
         \___ \| |\/| |  _| \___ \| |_| |
          ___) | |  | | |___ ___) |  _  |
         |____/|_|  |_|_____|____/|_| |_|



*/ // =================================================

template<class Cell_Type>
class SMesh : public Mesh
{
  template<class> friend class _MeshIterator;
  friend class Mesh;
  friend class MeshIoMsh;

public:
  typedef Cell_Type               CellT;
  typedef Point                   PointT;
  typedef Facet                   FacetT;
  typedef Corner                  CornerT;
  typedef SMesh<CellT>            MeshT;

  explicit SMesh(int spacedim) : Mesh(CellT::fep_tag, spacedim)
  {
  };

  ~SMesh(){};

private:
  SMesh(SMesh const&) : Mesh() {};

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

  // --------------------------------------------------- VERTEX STAR ---------------------------------------------------

  // TODO: implementar versão para 1d
  template<int celldim>
  int* vertexStar_Template(int C, int vC, int *iCs, int *viCs, typename EnableIf<(celldim==1)>::type* = NULL) const;

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
  template<int celldim>
  int* vertexStar_Template(int C, int vC, int *iCs, int *viCs, typename EnableIf<(celldim==2)>::type* = NULL) const;

  /** INTERNAL USE ONLY\n
   *  Returns all incident cells of a vertex.
   *  @param[in] C an incident cell of a vertex to start the search (initial cell).
   *  @param[in] vC vertex's local id in the initial cell.
   *  @param[out] iCs vector to put the incident cells.
   *  @param[out] viCs vertex's local ids in the incident cells.
   *  @return a pointer to the element following the end of the sequence iCs.
   *  @note the vectors iCs and viCs must have enough space to store the data.
   * */
  template<int celldim>
  int* vertexStar_Template(int C, int vC, int *iCs, int *viCs, typename EnableIf<(celldim==3)>::type* = NULL) const;


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
    return this->MeshT::vertexStar_Template<CellT::dim>(C, vC, iCs, viCs);
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
    return this->MeshT::vertexStar_Template<CellT::dim>(static_cast<PointT const*>(p)->PointT::getIncidCell(),
                                                        static_cast<PointT const*>(p)->PointT::getPosition(),
                                                        iCs, viCs);
  }


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
  int* nodeStar(int C, int nC, int *iCs, int *niCs) const;

  /** Returns all incident cells of a node.
   *  @param[in] p a pointer to the node.
   *  @param[out] iCs vector to put the incident cells.
   *  @param[out] niCs node's local ids in the initial cells.
   *  @return a pointer to the end of the sequence iCs.
   *  @note the vectors iCs and niCs must have enough space to store the data.
   *  @note this function assumes that each edge, face or volume has only one
   *        node in its interior, i.e., cells of third order are not supported.
   */
  int* nodeStar(Point const* p, int *iCs, int *niCs) const;

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
  template<int celldim>
  int* incidentFacets_Template(Point const* p, int *iFs, int *viFs, typename EnableIf<(celldim==1)>::type* = NULL) const;

  /** INTERNAL USE ONLY\n
   *  @brief Returns all incident facets to a vertex.
   *  @param[in] p a pointer to the vertex.
   *  @param[out] iFs vector to put the incident facets.
   *  @param[out] viFs vertex's local ids in the incident facets.
   *  @return a pointer to the element following the end of the sequence iFs.
   *  @note the vectors iCs and viCs must have enough space to store the data.
   *  @note does not support high order nodes, only vertices.
   */
  template<int celldim>
  int* incidentFacets_Template(Point const* p, int *iFs, int *viFs, typename EnableIf<(celldim==2)>::type* = NULL) const;

  /** INTERNAL USE ONLY\n
   *  @brief Returns all incident facets to a vertex.
   *  @param[in] p a pointer to the vertex.
   *  @param[out] iFs vector to put the incident facets.
   *  @param[out] viFs vertex's local ids in the incident facets.
   *  @return a pointer to the element following the end of the sequence iFs.
   *  @note the vectors iCs and viCs must have enough space to store the data.
   *  @note does not support high order nodes, only vertices.
   */
  template<int celldim>
  int* incidentFacets_Template(Point const* p, int *iFs, int *viFs, typename EnableIf<(celldim==3)>::type* = NULL) const;

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
    FEPIC_CHECK(this->MeshT::isVertex(p), "invalid C or vC", std::invalid_argument);
    return this->MeshT::incidentFacets_Template<CellT::dim>(p, iFs, viFs);
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
    FEPIC_CHECK(this->MeshT::isVertex(this->MeshT::getNodePtr(nodeid)), "invalid C or vC", std::invalid_argument);
    return this->MeshT::incidentFacets_Template<CellT::dim>(this->MeshT::getNodePtr(nodeid), iFs, viFs);
  }

  template<int celldim>
  Facet* nextBoundaryFacet_template(Facet const*f, typename EnableIf<(celldim==1)>::type* = NULL) const;

  template<int celldim>
  Facet* nextBoundaryFacet_template(Facet const*f, typename EnableIf<(celldim==2)>::type* = NULL) const;
  
  template<int celldim>
  Facet* nextBoundaryFacet_template(Facet const*f, typename EnableIf<(celldim==3)>::type* = NULL) const;
  
  /** given a boundary facet, return the next boundary facet.
   */ 
  Facet* nextBoundaryFacet(Facet const*f) const
  {
    return this->MeshT::nextBoundaryFacet_template<CellT::dim>(f);
  }


  void pushIncidCell2Point(Point *pt, int iC, int pos);

  // ---------------------------------------------------- EDGE STAR ---------------------------------------------------

  /** INTERNA USE ONLY\n
   * Returns all incidents cells of a edge.
   * @param[in] C an incident cell of a edge to start the search (initial cell).
   * @param[in] eC edge's local id in the initial cell.
   * @param[out] iCs vector to put the incident cells.
   * @param[out] eiCs edges's local ids in the incident cells.
   * @return the a pointer to the element following the end of the sequence iCs.
   */
  template<int celldim>
  int* edgeStar_Template(int C, int eC, int *iCs, int *eiCs, typename EnableIf<(celldim==3)>::type* = NULL) const;

  /** INTERNA USE ONLY\n
   * Returns all incidents cells of a edge.
   * @param[in] C an incident cell of a edge to start the search (initial cell).
   * @param[in] eC edge's local id in the initial cell.
   * @param[out] iCs vector to put the incident cells.
   * @param[out] eiCs edges's local ids in the incident cells.
   * @return a pointer to the element following the end of the sequence iCs.
   */
  template<int celldim>
  int* edgeStar_Template(int C, int eC, int *iCs, int *eiCs, typename EnableIf<(celldim==2)>::type* = NULL) const;

  // TODO: implementar uma versão para 1D
  template<int celldim>
  int* edgeStar_Template(int C, int eC, int *iCs, int *eiCs, typename EnableIf<(celldim==1)>::type* = NULL) const;

  /** Returns all incidents cells of a edge.
   * @param[in] C an incident cell of a edge to start the search (initial cell).
   * @param[in] eC edge's local id in the initial cell.
   * @param[out] iCs vector to put the incident cells.
   * @param[out] eiCs edges's local ids in the incident cells.
   * @return a pointer to the element following the end of the sequence iCs.
   */
  int* edgeStar(int C, int eC, int *iCs, int *eiCs) const
  {
    return this->edgeStar_Template<CellT::dim>(C, eC, iCs, eiCs);
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
    return this->edgeStar_Template<CellT::dim>(e->getIncidCell(),
                                               e->getPosition(),
                                               iCs, eiCs);
  }

  ///** Returns all incidents cells of a edge.
  // * @param[in] e a pointer to the edge.
  // * @param[out] iCs vector to put the incident cells.
  // * @param[out] eiCs edges's local ids in the incident cells.
  // * @return a pointer to the element following the end of the sequence iCs.
  // * @note if an edge is not a Corner, use the Facet version function.
  // */
  //int* edgeStar(Corner const* e, int *iCs, int *eiCs) const
  //{
  //  return this->edgeStar_Template<CellT::dim>(static_cast<CornerT const*>(e)->CornerT::getIncidCell(),
  //                                             static_cast<CornerT const*>(e)->CornerT::getPosition(),
  //                                             iCs, eiCs);
  //}

  // ---------------------------------------------------- FACE STAR ---------------------------------------------------

  template<int celldim>
  int* faceStar_Template(int C, int fC, int *iCs, int *fiCs, typename EnableIf<(celldim==3)>::type* = NULL) const;

  template<int celldim>
  int* faceStar_Template(int C, int fC, int *iCs, int *fiCs, typename EnableIf<(celldim==2)>::type* = NULL) const;

  template<int celldim>
  int* faceStar_Template(int C, int fC, int *iCs, int *fiCs, typename EnableIf<(celldim==1)>::type* = NULL) const;

  int* faceStar(int C, int fC, int *iCs, int *fiCs) const
  {
    return faceStar_Template<CellT::dim>(C, fC, iCs, fiCs);
  }

  // --------------------------------------------------- ADJACENCY -------------------------------------------------------

  void buildCellsAdjacency();

  template<int celldim>
  void buildCorners_Template(typename EnableIf<(celldim==1)>::type* = NULL);

  template<int celldim>
  void buildCorners_Template(typename EnableIf<(celldim==2)>::type* = NULL);

  template<int celldim>
  void buildCorners_Template(typename EnableIf<(celldim==3)>::type* = NULL);

  void buildNodesAdjacency();

  /// Constrói as informações topológicas da malha
  void buildAdjacency();

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

  // ----------------------------------------------------------------------------------------------------------------------

  // TODO: linear only ...
  void getCenterCoord(Cell const* cell, Real* x) const;
  void getCenterCoord(Facet const* facet, Real* x) const;
  void getCenterCoord(Corner const* corner, Real* x) const;

  bool isVertex(CellElement const* p) const
  {
    return p->getPosition() < CellT::n_vertices;
  }

  bool inSingleCell(Point const* p) const;
  bool inSingleCell(Corner const* p) const;
  bool inSingleCell(Facet const* p) const;

  /** @brief estimate of how the containers will grow.
   *  @param factor the size
   */


  int cellDim() const
  {
    return CellT::dim;
  }

  ECellType cellType() const
  {
    return CellT::fep_tag;
  }

  EMshTag cellMshTag() const
  {
    return CellT::msh_tag;
  }

  int nodesPerCell() const
  {
    return CellT::n_nodes;
  }

  int nodesPerFacet() const
  {
    return CellT::n_nodes_per_facet;
  }

  int nodesPerCorner() const
  {
    return CellT::Derived::n_nodes_per_facet;
  }

  int verticesPerCell() const
  {
    return CellT::n_vertices;
  }

  int verticesPerFacet() const
  {
    return CellT::n_vertices_per_facet;
  }

  int verticesPerCorner() const
  {
    return CellT::Derived::n_vertices_per_facet;
  }

  void printInfo() const;

  void printStatistics() const;

  /// create a cell (but not put in the mesh)
  Cell*   createCell() const;
  /// create a point (but not put in the mesh)
  Point*  createPoint() const;
  /// create a facet (but not put in the mesh)
  Facet*  createFacet() const;
  /// create a corner (but not put in the mesh)
  Corner* createCorner() const;
          
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

};


#endif // SMesh


