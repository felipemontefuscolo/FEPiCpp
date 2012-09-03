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
#include "../util/macros.hpp"
#include "../util/timer.hpp"
#include "elements/cellcore.hpp"
#include "elements/elements.hpp"
#include "../util/list_type.hpp"
#include "mesh_iterators.hpp"

template<class CT, int SD> class SMesh;

typedef _MeshIterator<Cell>   cell_iterator;
typedef _MeshIterator<Point>  point_iterator;
typedef _MeshIterator<Facet>  facet_iterator;
typedef _MeshIterator<Corner> corner_iterator;


class Mesh
{
protected:
  explicit Mesh(int sd=0, ECellType fept=UNDEFINED_CELLT);

  // iterators
  template<class> friend class _MeshIterator;
  friend class MeshIoMsh;

  virtual Cell*   incCell(Cell*) = 0;
  virtual Point*  incPoint(Point*) = 0;
  virtual Facet*  incFacet(Facet*) = 0;
  virtual Corner* incCorner(Corner*) = 0;

  virtual Cell*   decCell(Cell*) = 0;
  virtual Point*  decPoint(Point*) = 0;
  virtual Facet*  decFacet(Facet*) = 0;
  virtual Corner* decCorner(Corner*) = 0;


public:

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

  typedef Mesh* (*CreatorMemFunPtr)();

  template<class CT, int SD>
  static Mesh* create()
  {
    return new SMesh<CT,SD>;
  }

  static Mesh* create(ECellType type, int spacedim = -1);



  virtual int cellDim() const = 0;
  virtual int spaceDim() const = 0;
  virtual ECellType cellType() const = 0;
  virtual EMshTag cellMshTag() const = 0;

  virtual bool isVertex(Point const* p) const = 0;
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

  virtual Cell* getCell(int nth) = 0;
  virtual Facet* getFacet(int nth) = 0;
  virtual Corner* getCorner(int nth) = 0;
  virtual Point* getNode(int nth) = 0;

  virtual const Cell* getCell(int nth) const = 0;
  virtual const Facet* getFacet(int nth) const = 0;
  virtual const Corner* getCorner(int nth) const = 0;
  virtual const Point* getNode(int nth) const = 0;

  virtual void disablePoint(int id) = 0;
  virtual void disableCorner(int id) = 0;
  virtual void disableFacet(int id) = 0;
  virtual void disableCell(int id) = 0;

  virtual int getCellId(Cell const*) const = 0;
  virtual int getPointId(Point const*) const = 0;
  virtual int getFacetId(Facet const*) const = 0;
  virtual int getCornerId(Corner const*) const = 0;

  virtual void getCellNodesId(Cell const* cell, int* result) const = 0;
  virtual void getCellVerticesId(Cell const* cell, int* result) const = 0;
  virtual void getCellFacetsId(Cell const* cell, int* result) const = 0;
  virtual void getCellCornersId(Cell const* cell, int* result) const = 0;
  virtual void getFacetNodesId(Facet const* facet, int* result) const = 0;
  virtual void getCornerNodesId(Corner const* corner, int* result) const = 0;

  virtual void getCellNodesId(int id, int *result) const = 0;
  virtual void getFacetNodesId(int id, int *result) const = 0;
  virtual void getCornerNodesId(int id, int *result) const = 0;

  virtual int getCellContigId(int id) const = 0;
  virtual int getFacetContigId(int id) const = 0;
  virtual int getCornerContigId(int id) const = 0;
  virtual int getNodeContigId(int id) const = 0;
  virtual void getCellsContigId(int* first, int const* last, int* result) const = 0;
  virtual void getFacetsContigId(int* first, int const* last, int* result)const = 0;
  virtual void getCornersContigId(int* first, int const* last, int* result) const = 0;
  virtual void getNodesContigId(int* first, int const* last, int* result) const = 0;

  virtual void getCellNodesContigId(Cell const* cell, int* result) const = 0;
  virtual void getFacetNodesContigId(Facet const* facet, int* result) const = 0;
  virtual void getCornerNodesContigId(Corner const* corner, int* result) const = 0;

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
    Real const* c;
    while(first != last)
    {
      c = this->getNode(*first++)->getCoord();
      for (int d = 0; d < this->spaceDim(); ++d)
        *X++ = *c++;
    }
  }

  // TODO: linear only ...
  virtual void getCenterCoord(Cell const* cell, Real* x) const = 0;
  virtual void getCenterCoord(Facet const* facet, Real* x) const = 0;
  virtual void getCenterCoord(Corner const* corner, Real* x) const = 0;

  virtual bool getFacetIdFromVertices(int const* vtcs, int &fid) =0;
  virtual bool getCornerIdFromVertices(int const* vtcs, int &fid) =0;

  virtual int* edgeStar(int C, int eC, int *iCs, int *eiCs) const = 0;
  virtual int* edgeStar(Facet  const*, int *iCs, int *eiCs) const = 0;
  virtual int* edgeStar(Corner const*, int *iCs, int *eiCs) const = 0;
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


  virtual void qBuildAdjacency(bool b)
  {
    _build_adjacency = b;
  };
  virtual bool qBuildAdjacency()
  {
    return _build_adjacency;
  };
  virtual void buildAdjacency() = 0;
  virtual void buildNodesAdjacency() = 0;

  virtual void setUpConnectedComponentsId() = 0;
  virtual int  numConnectedComponents() const = 0;
  virtual void setUpBoundaryComponentsId() = 0;
  virtual int  numBoundaryComponents() const = 0;
  
  virtual void getConnectedComponentsPicks(int *comps, int *cells) const = 0;
  virtual void getBoundaryComponentsPicks(int *comps, int *facets) const = 0;

  virtual bool inBoundary(Point const* p) const = 0;
  virtual bool inBoundary(Facet const* p) const = 0;
  virtual bool inBoundary(Corner const* p) const = 0;

  virtual int pushCell(Cell const* C) = 0;
  virtual int pushPoint(Point const* P) = 0;
  virtual int pushFacet(Facet const* h) = 0;
  virtual int pushCorner(Corner const* b) = 0;

  virtual Cell*   pushCell() = 0;
  virtual Point*  pushPoint() = 0;
  virtual Facet*  pushFacet() = 0;
  virtual Corner* pushCorner() = 0;

  virtual Cell*   createCell() const = 0;
  virtual Point*  createPoint() const = 0;
  virtual Facet*  createFacet() const = 0;
  virtual Corner* createCorner() const = 0;

  virtual int numCells() const = 0;
  virtual int numCellsTotal() const = 0;
  virtual int numNodes() const = 0;
  virtual int numNodesTotal() const = 0;
  virtual int numVertices() const = 0;
  virtual int numFacets() const = 0;
  virtual int numFacetsTotal() const = 0;
  virtual int numCorners() const = 0;
  virtual int numCornersTotal() const = 0;

  virtual int numNodesPerCell() const = 0;
  virtual int numNodesPerFacet() const = 0;
  virtual int numNodesPerCorner() const = 0;

  virtual int numVerticesPerCell() const = 0;
  virtual int numVerticesPerFacet() const = 0;
  virtual int numVerticesPerCorner() const = 0;

  virtual int numFacetsPerCell() const = 0;
  virtual int numCornersPerCell()  const = 0;
  virtual int numCornersPerFacet() const = 0;

  virtual void resizePointL(unsigned size) = 0;
  virtual void resizeCellL(unsigned size) = 0;
  virtual void resizeFacetL(unsigned size) = 0;
  virtual void resizeCornerL(unsigned size) = 0;
  virtual void reservePointL(unsigned size) = 0;
  virtual void reserveCellL(unsigned size) = 0;
  virtual void reserveFacetL(unsigned size) = 0;
  virtual void reserveCornerL(unsigned size) = 0;

  static unsigned estimateNumFacets(unsigned nc, ECellType t);
  static unsigned estimateNumCorners(unsigned nc, ECellType t);

  /** @brief estimate of how the containers will grow.
   *  @param factor the size
   */
  virtual void setGrowFactor(float factor) = 0;

  virtual ~Mesh() {};


  // time measure
  Timer timer;

//protected:
  ECellType _cell_fep_tag;
  EMshTag   _cell_msh_tag;
  bool      _build_adjacency;

};


/* ====================================================

          ____  __  __ _____ ____  _   _
         / ___||  \/  | ____/ ___|| | | |
         \___ \| |\/| |  _| \___ \| |_| |
          ___) | |  | | |___ ___) |  _  |
         |____/|_|  |_|_____|____/|_| |_|



*/ // =================================================

template<class Cell_Type, int Spacedim = Cell_Type::dim>
class SMesh : public Mesh
{
  template<class> friend class _MeshIterator;
  friend class Mesh;
  friend class MeshIoMsh;

public:

  typedef Cell_Type               CellT;
  typedef PointX<Spacedim>        PointT;
  typedef Facet                   FacetT;
  typedef Corner                  CornerT;
  typedef SMesh<CellT, Spacedim>  MeshT;

  //typedef std::vector<CellT>    CellList;
  //typedef std::vector<PointT>   PointList;
  //typedef std::vector<Facet>    FacetList;
  //typedef std::vector<Corner>   CornerList;

  typedef SeqList<CellT, std::vector<CellT>, SetVector<int> >     CellList;
  typedef SeqList<PointT, std::vector<PointT>, SetVector<int> >   PointList;
  typedef SeqList<FacetT, std::vector<FacetT>, SetVector<int> >   FacetList;
  typedef SeqList<CornerT, std::vector<CornerT>, SetVector<int> > CornerList;

  typedef typename CellList  ::iterator CellIteratorT;
  typedef typename PointList ::iterator PointIteratorT;
  typedef typename FacetList ::iterator FacetIteratorT;
  typedef typename CornerList::iterator CornerIteratorT;

  typedef typename CellList  ::const_iterator CellConstIteratorT;
  typedef typename PointList ::const_iterator PointConstIteratorT;
  typedef typename FacetList ::const_iterator FacetConstIteratorT;
  typedef typename CornerList::const_iterator CornerConstIteratorT;


  explicit SMesh() : Mesh(-1, CellT::fep_tag)
  {
  };

  ~SMesh(){};

protected:
  SMesh(SMesh const&) : Mesh() {};

  Cell*   incCell(Cell*);
  Point*  incPoint(Point*);
  Facet*  incFacet(Facet*);
  Corner* incCorner(Corner*);

  Cell*   decCell(Cell*);
  Point*  decPoint(Point*);
  Facet*  decFacet(Facet*);
  Corner* decCorner(Corner*);


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
   *  @note iVs[k] = getCell(iCd[k])->getNodeId(viCs[k]);
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
    FEPIC_CHECK(this->MeshT::isVertex(this->MeshT::getNode(nodeid)), "invalid C or vC", std::invalid_argument);
    return this->MeshT::incidentFacets_Template<CellT::dim>(this->MeshT::getNode(nodeid), iFs, viFs);
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
  int* edgeStar(Facet const* e, int *iCs, int *eiCs) const
  {
    return this->edgeStar_Template<CellT::dim>(static_cast<FacetT const*>(e)->FacetT::getIncidCell(),
                                               static_cast<FacetT const*>(e)->FacetT::getPosition(),
                                               iCs, eiCs);
  }

  /** Returns all incidents cells of a edge.
   * @param[in] e a pointer to the edge.
   * @param[out] iCs vector to put the incident cells.
   * @param[out] eiCs edges's local ids in the incident cells.
   * @return a pointer to the element following the end of the sequence iCs.
   * @note if an edge is not a Corner, use the Facet version function.
   */
  int* edgeStar(Corner const* e, int *iCs, int *eiCs) const
  {
    return this->edgeStar_Template<CellT::dim>(static_cast<CornerT const*>(e)->CornerT::getIncidCell(),
                                               static_cast<CornerT const*>(e)->CornerT::getPosition(),
                                               iCs, eiCs);
  }

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
   * @warning visited() flags are used.
   */ 
  void _setConnectedComponentsId(Cell * c_ini, int cc_id);
  void _setBoundaryComponentsId(Facet * f_ini, int bc_id);

  void setUpConnectedComponentsId();
  int numConnectedComponents() const
  {
    return _connected_compL.size();
  }
  
  void setUpBoundaryComponentsId();
  int numBoundaryComponents() const
  {
    return _boundary_compL.size();
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

  bool isVertex(Point const* p) const
  {
    return static_cast<PointT const*>(p)->PointT::getPosition() < CellT::n_vertices;
  }

  bool inSingleCell(Point const* p) const;
  bool inSingleCell(Corner const* p) const;
  bool inSingleCell(Facet const* p) const;

  /** @brief estimate of how the containers will grow.
   *  @param factor the size
   */
  void setGrowFactor(float factor)
  {
    _cellL.setGrowFactor(factor);
    _pointL.setGrowFactor(factor);
    _facetL.setGrowFactor(factor);
    _cornerL.setGrowFactor(factor);
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

  int cellDim() const
  {
    return CellT::dim;
  }

  int spaceDim() const
  {
    return Spacedim;
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

  /** Retorna a n-ésima celula (adivinha o tipo de célula pelo tipo da malha)
  */
  CellT* getCell(int nth)
  {
    FEPIC_CHECK(unsigned(nth)<this->_cellL.total_size(), "invalid index", std::out_of_range);
    return &_cellL[nth];
  }
  /** Retorna a n-ésima facet
  */
  FacetT* getFacet(int nth)
  {
    //printf("nth=%d, face.size()=%d\n", nth, (int)_facetL.size());
    FEPIC_CHECK(unsigned(nth)<this->_facetL.total_size(), "invalid index", std::out_of_range);
    return &_facetL[nth];
  }
  /** Retorna a n-ésima corner
  */
  CornerT* getCorner(int nth)
  {
    FEPIC_CHECK(unsigned(nth)<this->_cornerL.total_size(), "invalid index", std::out_of_range);
    return &_cornerL[nth];
  }
  /** Retorna o n-ésimo nó da malha.
  */
  PointT* getNode(int nth)
  {
    FEPIC_CHECK(unsigned(nth)<this->_pointL.total_size(), "invalid index", std::out_of_range);
    return &_pointL[nth];
  }


  const CellT* getCell(int nth) const
  {
    FEPIC_CHECK(unsigned(nth)<this->_cellL.total_size(), "invalid index", std::out_of_range);
    return &_cellL[nth];
  }
  const FacetT* getFacet(int nth) const
  {
    FEPIC_CHECK(unsigned(nth)<this->_facetL.total_size(), "invalid index", std::out_of_range);
    return &_facetL[nth];
  }
  const CornerT* getCorner(int nth) const
  {
    FEPIC_CHECK(unsigned(nth)<this->_cornerL.total_size(), "invalid index", std::out_of_range);
    return &_cornerL[nth];
  }
  const PointT* getNode(int nth) const
  {
    FEPIC_CHECK(unsigned(nth)<this->_pointL.total_size(), "invalid index", std::out_of_range);
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

  int getCellId(Cell const* a) const
  {
    return _cellL.getDataId(static_cast<CellT const*>(a));
  }
  int getPointId(Point const* a) const
  {
    return _pointL.getDataId(static_cast<PointT const*>(a));
  }
  int getFacetId(Facet const* a) const
  {
    return _facetL.getDataId(static_cast<FacetT const*>(a));
  }
  int getCornerId(Corner const* a) const
  {
    return _cornerL.getDataId(static_cast<CornerT const*>(a));
  }

  void getCellNodesId(Cell const* cell, int *result) const
  {
    static_cast<CellT const*>(cell)->CellT::getNodesId(result);
  }
  void getCellVerticesId(Cell const* cell, int *result) const
  {
    static_cast<CellT const*>(cell)->CellT::getVerticesId(result);
  }
  void getCellFacetsId(Cell const* cell, int *result) const
  {
    static_cast<CellT const*>(cell)->CellT::getFacetsId(result);
  }
  void getCellCornersId(Cell const* cell, int *result) const
  {
    static_cast<CellT const*>(cell)->CellT::getCornersId(result);
  }
  void getFacetNodesId(Facet const* facet, int *result) const
  {
    int icell = static_cast<FacetT const*>(facet)->FacetT::getIncidCell();
    int pos   = static_cast<FacetT const*>(facet)->FacetT::getPosition();
    this->MeshT::getCell(icell)->CellT::getFacetNodesId(pos, result);
  }
  void getCornerNodesId(Corner const* corner, int *result) const
  {
    int icell = static_cast<CornerT const*>(corner)->CornerT::getIncidCell();
    int pos   = static_cast<CornerT const*>(corner)->CornerT::getPosition();
    this->MeshT::getCell(icell)->CellT::getCornerNodesId(pos, result);
  }

  void getCellNodesId(int id, int *result) const
  {
    this->getCellNodesId(this->MeshT::getCell(id), result);
  }
  void getFacetNodesId(int id, int *result) const
  {
    this->MeshT::getFacetNodesId(this->MeshT::getFacet(id), result);
  }
  void getCornerNodesId(int id, int *result) const
  {
    this->MeshT::getCornerNodesId(this->MeshT::getCorner(id), result);
  }

  void getCellNodesContigId(Cell const* cell, int* result) const
  {
    this->MeshT::getCellNodesId(cell, result);
    for (int i = 0; i < CellT::n_nodes; ++i)
    {
      *result = _pointL.contiguousId(*result);
      ++result;
    }
  }
  void getFacetNodesContigId(Facet const* facet, int* result) const
  {
    this->MeshT::getFacetNodesId(facet, result);
    for (int i = 0; i < CellT::n_nodes_per_facet; ++i)
    {
      *result = _pointL.contiguousId(*result);
      ++result;
    }
  }
  void getCornerNodesContigId(Corner const* corner, int* result) const
  {
    this->MeshT::getCornerNodesId(corner, result);
    for (int i = 0; i < CellT::n_nodes_per_corner; ++i)
    {
      *result = _pointL.contiguousId(*result);
      ++result;
    }
  }

  bool inBoundary(Point const* p) const
  {
    
    return static_cast<PointT const*>(p)->PointT::inBoundary();
    
    // OLD
    ////return this->MeshT::getCell(static_cast<PointT const*>(p)->PointT::getIncidCell())->CellT::inBoundary();
    //CellT const* icell = this->MeshT::getCell(static_cast<PointT const*>(p)->PointT::getIncidCell());
    //if (!icell->inBoundary())
    //  return false;
    //
    //const int m = static_cast<PointT const*>(p)->PointT::getPosition();
    //if (CellT::dim==1)
    //{
    //  if (m>1)
    //    return false;
    //  else
    //    return true;
    //}
    //if (this->MeshT::isVertex(p))
    //{
    //  if (icell->CellT::getIncidCell(CellT::table_vC_x_fC[m][0]) < 0) return true;
    //  if (icell->CellT::getIncidCell(CellT::table_vC_x_fC[m][1]) < 0) return true;
    //  if (CellT::dim>2)
    //    if (icell->CellT::getIncidCell(CellT::table_vC_x_fC[m][2]) < 0) return true;
    //  return false;
    //}
    //else
    //{
    //  const int q = m - CellT::n_vertices;
    //  if (CellT::dim==2)
    //  {
    //    if (q>=CellT::n_facets)
    //      return false;
    //    return icell->CellT::getIncidCell(q)<0;
    //  }
    //  else
    //  {
    //    const int r = q - CellT::n_corners;
    //    const int s = r - CellT::n_facets;
    //
    //    if (s>=0) return false;
    //
    //    if (r>=0) return icell->CellT::getIncidCell(r)<0;
    //
    //    if (icell->CellT::getIncidCell(CellT::table_bC_x_fC[q][0]) < 0) return true;
    //    if (icell->CellT::getIncidCell(CellT::table_bC_x_fC[q][1]) < 0) return true;
    //    return false;
    //  }
    //}


  }
  bool inBoundary(Facet const* f) const
  {
    CellT const* icell = this->MeshT::getCell(static_cast<FacetT const*>(f)->FacetT::getIncidCell());
    if (icell->CellT::getIncidCell(static_cast<FacetT const*>(f)->FacetT::getPosition()) < 0)
      return true;
    else
      return false;
  }
  bool inBoundary(Corner const* f) const
  {
    CellT const* icell = this->MeshT::getCell(static_cast<CornerT const*>(f)->CornerT::getIncidCell());
    if (!icell->CellT::inBoundary())
      return false;
    if (CellT::dim==2)
    {
      PointT const* point = this->MeshT::getNode( icell->CellT::getNodeId(  static_cast<CornerT const*>(f)->CornerT::getPosition() ) );
      return this->MeshT::inBoundary(point);
    }
    else
    {
      int const m = static_cast<CornerT const*>(f)->CornerT::getPosition();
      if (icell->CellT::getIncidCell(CellT::table_bC_x_fC[m][0]) < 0) return true;
      if (icell->CellT::getIncidCell(CellT::table_bC_x_fC[m][1]) < 0) return true;
      return false;
    }
  }

  /** Adiciona uma célula e retorna seu id.
  */
  int pushCell(Cell const* C);

  /** Adiciona um ponto e retorna seu id.
  */
  int pushPoint(Point const* P);

  /** Adiciona uma facet
  *  @param h A facet-xxxx a ser adicionada.
  *  @return A posição da facet-xxxx na lista
  */
  int pushFacet(Facet const* h);

  /** Adiciona uma corner-xxxx
  *  @param h A corner-xxxx a ser adicionada.
  *  @return A posição da corner-xxxx na lista
  */
  int pushCorner(Corner const* h);

  /// adds an element in the mesh and returns a pointer to it.
  Cell*   pushCell();
  /// adds a point in the mesh and returns a pointer to it.
  Point*  pushPoint();
  /// adds a facet in the mesh and returns a pointer to it.
  Facet*  pushFacet();
  /// adds a corner in the mesh and returns a pointer to it.
  Corner* pushCorner();
         
  /// create a cell (but not put in the mesh)
  Cell*   createCell() const;
  /// create a point (but not put in the mesh)
  Point*  createPoint() const;
  /// create a facet (but not put in the mesh)
  Facet*  createFacet() const;
  /// create a corner (but not put in the mesh)
  Corner* createCorner() const;
          
  /** Retorna o número células
  *  @note não conta com o/a(s) marcado/a(s) como killed.
  */
  int numCells() const
  {
    return _cellL.size();
  }

  /** Retorna o número de células.
  * @note incluindo o/a(s) marcado/a(s) como killed.
  */
  int numCellsTotal() const
  {
    return _cellL.total_size();
  }

  /** Retorna o número de nós.
  *  @note não conta com o/a(s) marcado/a(s) como killed.
  */
  int numNodes() const
  {
    return _pointL.size();
  }

  /** Retorna no número de nós.
  *  @note incluindo o/a(s) marcado/a(s) como killed.
  */
  int numNodesTotal() const
  {
    return _pointL.total_size();
  }

  int numVertices() const;

  /** Retorna número de facets.
  * @note não conta com o/a(s) marcado/a(s) como killed.
  */
  int numFacets() const
  {
    return _facetL.size();
  }

  /** Retorna o número de facets.
  *  @note incluindo o/a(s) marcado/a(s) como killed.
  */
  int numFacetsTotal() const
  {
    return _facetL.total_size();
  }

  /** Retorna número de corners.
  * @note não conta com o/a(s) marcado/a(s) como killed.
  */
  int numCorners() const
  {
    return _cornerL.size();
  }

  /** Retorna o número de corners.
  *  @note incluindo o/a(s) marcado/a(s) como killed.
  */
  int numCornersTotal() const
  {
    return _cornerL.total_size();
  }

  int numNodesPerCell() const
  {
    return CellT::n_nodes;
  }
  int numNodesPerFacet() const
  {
    return CellT::n_nodes_per_facet;
  }
  int numNodesPerCorner() const
  {
    return CellT::n_nodes_per_corner;
  }

  int numVerticesPerCell() const
  {
    return CellT::n_vertices;
  }
  int numVerticesPerFacet() const
  {
    return CellT::n_vertices_per_facet;
  }
  int numVerticesPerCorner() const
  {
    return CellT::n_vertices_per_corner;
  }

  int numFacetsPerCell() const
  {
    return CellT::n_facets;
  }
  int numCornersPerCell() const
  {
    return CellT::n_corners;
  }
  int numCornersPerFacet() const
  {
    return CellT::Derived::n_facets;
  }

  void resizePointL(unsigned size);
  void resizeCellL(unsigned size);
  void resizeFacetL(unsigned size);
  void resizeCornerL(unsigned size);

  void reservePointL(unsigned size);
  void reserveCellL(unsigned size);
  void reserveFacetL(unsigned size);
  void reserveCornerL(unsigned size);

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

private:

  // entities
  CellList      _cellL;
  PointList     _pointL;
  FacetList     _facetL;
  CornerList    _cornerL;

  std::map<int, int>  _connected_compL; // connected component vs initial cell id
  std::map<int, int>  _boundary_compL;  // boundary component vs initialfacet id
};


#endif // SMesh


