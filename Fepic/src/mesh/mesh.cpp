#include "mesh.hpp"
#include "boost/scoped_ptr.hpp"
#include <vector>
#include <set>
#include "../util/misc2.hpp"
#include "../util/sorted_vec.hpp"
#include <map>
#include <deque>
#include <cmath>
#ifdef FEP_HAS_OPENMP
  #include <omp.h>
#endif

#ifdef FEP_HAS_OMPTL
  #include <omptl/omptl_algorithm>
#else
  //using namespace omptl = std;
  namespace omptl = std;
#endif
#ifdef __GLIBCXX__
#  include <tr1/array>
#  include <tr1/memory>
#  include <tr1/cstdint>
#  include <tr1/type_traits>
#else
#  ifdef __IBMCPP__
#    define __IBMCPP_TR1__
#  endif
#  include <array>
#  include <memory>
#  include <cstdint>
#  include <type_traits>
#endif
#include "boost/scoped_ptr.hpp"

Mesh::Mesh(ECellType const fept, int spacedim)
{
  FEPIC_CHECK(spacedim>=1 && spacedim<4, "invalid space dimension", std::runtime_error);
  
  Cell * cell = Cell::create(fept);

  FEPIC_CHECK(cell != NULL, "invalid cell type", std::runtime_error);

  m_spacedim = spacedim; // BROKEN
  m_cellm_fep_tag = fept;
  m_cell_msh_tag = ctype2mshTag(fept);
  //m_cell_msh_tag = cell->getMshTag();
  m_build_adjacency = true;
  
  m_is_parametric_cell = cell->isParametric();
  m_cell_dim = cell->dim();
  m_n_nodes_per_cell = cell->numNodes();
  m_n_nodes_per_facet = cell->numNodesPerFacet();
  m_n_nodes_per_corner = cell->numNodesPerCorner();
  m_n_vertices_per_cell = cell->numVertices();
  m_n_vertices_per_facet = cell->numVerticesPerFacet();
  m_n_vertices_per_corner = cell->numVerticesPerCorner();
  m_n_facets_per_cell = cell->numFacets();
  m_n_corners_per_cell = cell->numCorners();
  m_n_corners_per_facet = cell->numCornersPerFacet();
  m_cell_has_edge_nodes = cell->hasEdgeNodes();
  m_cell_has_face_nodes = cell->hasFaceNodes();
  m_cell_has_volume_nodes = cell->hasVolumeNodes();

  delete cell;

  timer = Timer();
}

/** @param nc_ number of cells.
 *  @param type mesh cell type.
 */
unsigned Mesh::estimateNumFacets(unsigned nc_, ECellType type)
{
  ECellClass cc = ctype2cclass(type);
  double nc = static_cast<float>(nc_);  // nc = num cells

#define _ROUND_2_INT_(x) static_cast<unsigned>(floor((x) + 0.5))

  switch (cc)
  {
    case TRIANGLE:
      return _ROUND_2_INT_(1.5*nc + sqrt(2.*nc));
    case QUADRANGLE:
      return _ROUND_2_INT_(2.*(nc + sqrt(nc)));
    case TETRAHEDRON:
      return _ROUND_2_INT_(2.*nc+pow(6.,1./3)*pow(nc,2./3));
    case HEXAHEDRON:
      return _ROUND_2_INT_(3.*nc+3.*pow(nc,2./3));
    default:
      return 0;
  }

#undef _ROUND_2_INT_

}
unsigned Mesh::estimateNumCorners(unsigned nc_, ECellType type)
{
  ECellClass cc = ctype2cclass(type);
  double nc = static_cast<float>(nc_); // nc = num cells

#define _ROUND_2_INT_(x) static_cast<unsigned>(floor((x) + 0.5))

  switch (cc)
  {
    case TRIANGLE:
      return _ROUND_2_INT_(0.5*nc);
    case QUADRANGLE:
      return _ROUND_2_INT_(nc);
    case TETRAHEDRON:
      return _ROUND_2_INT_((7.*nc+9.*pow(6.,1./3)*pow(nc,2./3)+3*pow(6.,2./3)*pow(nc,1./3))/6.);
    case HEXAHEDRON:
      return _ROUND_2_INT_(3*nc+6*pow(nc,2./3)+3*pow(nc,1./3));
    default:
      return 0;
  }

#undef _ROUND_2_INT_
}



// =====================================================================================
// =====================================================================================

template<> FEP_STRONG_INLINE Cell  * Mesh::entityPtr<Cell  >(int ith) {return this->getCellPtr(ith);}
template<> FEP_STRONG_INLINE Facet * Mesh::entityPtr<Facet >(int ith) {return this->getFacetPtr(ith);}
template<> FEP_STRONG_INLINE Point * Mesh::entityPtr<Point >(int ith) {return this->getNodePtr(ith);}
template<> FEP_STRONG_INLINE Corner* Mesh::entityPtr<Corner>(int ith) {return this->getCornerPtr(ith);}

// -----------------------------------------------------------------------------------

int Mesh::numVertices() const
{
  int num_vtcs = 0;
  int const numm_nodes_total = this->numNodesTotal();
  FEP_PRAGMA_OMP(parallel shared(num_vtcs))
  {
    int num_vtcs_local = 0;
    Point const* p;

    FEP_PRAGMA_OMP(for)
    for (int i=0; i<numm_nodes_total; ++i)
    {
      p = this->getNodePtr(i);
      if (p->isDisabled())
        continue;
      if (this->isVertex(p))
        ++num_vtcs_local;
    }

    FEP_PRAGMA_OMP(critical)
    num_vtcs += num_vtcs_local;
  }
  return num_vtcs;
}


void Mesh::printInfo() const
{
  printf("elem type:  %s\n",       ctypeName(this->cellType()));
  printf("space dim:  %d\n", this->spaceDim()                  );
  printf("# vertices: %d\n", this->numVertices()               );
  printf("# nodes:    %d\n", this->numNodes()                  );
  printf("# cells:    %d\n", this->numCells()                  );
  printf("# facets:   %d\n", this->numFacets()                 );
  printf("# corners:  %d\n", this->numCorners()                );
}

void Mesh::printStatistics() const
{
  printf("\n"
         "Cell list: \n"
         "  -size:           %d\n"
         "Node list: \n"
         "  -size:           %d\n"
         "Facet list: \n"
         "  -size:           %d\n"
         "Corner list: \n"
         "  -size:           %d\n",
         (int)m_cells_list.totalSize(),
         (int)m_points_list.totalSize(),
         (int)m_facets_list.totalSize(),
         (int)m_corners_list.totalSize()
         );
}



// -------------------------------------------------- VERTEX STAR ---------------------------------------------------

int* Mesh::vertexStar_1D(int C, int vC, int *iCs, int *viCs) const
{
  // do nothing
  // annoying compiler
  C++; vC++; iCs++; viCs++;
  return iCs;
}

int* Mesh::vertexStar_2D(int C, int vC, int *iCs, int *viCs) const
{
  int const nvpc = this->numVerticesPerCell();
  int const nfpc = this->numFacetsPerCell();
  
  FEPIC_CHECK(vC<nvpc && C>=0, "invalid C or vC", std::invalid_argument);

  Point const* pt = this->getNodePtr(this->getCellPtr(C)->getNodeId(vC));
  int const n_connected_comps = pt->numConnectedComps();

  for (int cc = 0; cc < n_connected_comps; ++cc)
  {
    //C = pt->getIncidCell();
    //vC = pt->getPosition();
    pt->getIthIncidCell(cc,C,vC);
    
    Cell const* cell;
    int g, D=C, vD=vC, q=0;

    *iCs++  = C;
    *viCs++ = vC;

    for (;;)
    {
      cell = this->getCellPtr(D);
      g = (nfpc + vD -q) % nfpc;
      D = cell->getIncidCell(g);
      if (D<0)
      {
        if(q) break;
        D = C;
        vD = vC;
        q = 1;
        continue;
      }
      if (D==C)
        break;
      g = cell->getIncidCellPos(g);
      vD = (g + 1 - q) % nfpc;
      *iCs++  = D;
      *viCs++ = vD;
      //printf("DEBUG     %d, %d de %d, %d\n",D, vD, C, vC);
    }
    
  }

  *iCs = -1;
  *viCs = -1;
  return iCs;
}

int* Mesh::vertexStar_3D(int C, int vC, int *iCs, int *viCs) const
{
  const int nvpc = this->numVerticesPerCell();
  const int nfpc = this->numFacetsPerCell();
  const int nvpf = this->numVerticesPerFacet();
  
  FEPIC_CHECK(vC<nvpc && C>=0, "invalid C or vC", std::invalid_argument);

  Cell const*cell = this->getCellPtr(C);
  int f, D, vD, g, anc, vf, vg, gv;
  int const* iCs_beg = iCs;
  int* iCs_end  = iCs;
  int* viCs_end = viCs;
  int* iter;
  int pos;

  // fv = f-th incident facet of the vertex vC
  // 0 means fv = {}
  // 1<<0 means fv = {0}
  // 1<<1 means fv = {1}
  // 1<<2 means fv = {2}
  // 1<<0 | 1<<1 means fv = {0,1}
  // and so on ...
  // this vector store the facets of each cell that already been covered.
  int iCs_fv[FEPIC_MAX_ICELLS];
  int *fv_end = iCs_fv;

  *iCs_end++  = C;
  *viCs_end++ = vC;
  *fv_end++   = 0;

  // for notation, see file "ini_static_mem.cpp"
  for (;;)
  {
    pos = (int)(iCs - iCs_beg);
    // put the 3 neighbors in the stack
    for (int fv = 0; fv < 3; ++fv) // each vertice has three faces for both tetrahedron and hexahedron
    {
      // (C, f, vf, fv) <===> (D, g, vg, gv)

      FEPIC_CHECK(pos < FEPIC_MAX_ICELLS, "FEPIC_MAX_ICELLS exceeded", std::runtime_error);

      // ask if this facet has already been covered
      if ((fv+1+(fv==2)) & iCs_fv[pos])
        continue;

      f = cell->table_vC_x_fC( vC , fv );
      D = cell->getIncidCell(f);

      if (D<0)
        continue;

      if ((iter = std::find((int*)iCs_beg, iCs_end, D)) != iCs_end)
      {
        if (iter > iCs)
        {
          g = cell->getIncidCellPos(f);
          anc = cell->getIncidCellAnch(f);
          vf = cell->table_vC_x_fC( vC , fv+3 );
          vg = (nfpc+1 - vf - anc) % nvpf;
          gv = cell->table_fC_x_vC( g , vg + nvpf);
          iCs_fv[(int)(iter-iCs_beg)] |= (1 << gv);
        }
        continue;
      }
      g = cell->getIncidCellPos(f);
      anc = cell->getIncidCellAnch(f);

      vf = cell->table_vC_x_fC( vC , fv+3 );
      vg = (nfpc+1 - vf - anc) % nvpf;
      vD  = cell->table_fC_x_vC( g , vg );

      *iCs_end++  = D;
      *viCs_end++ = vD;

      gv = cell->table_fC_x_vC( g , vg + nvpf);
      *fv_end++ = (1<<gv) ;

    }
    if (++iCs != iCs_end)
    {
      cell = this->getCellPtr(*iCs);
      vC = *++viCs;
    }
    else
      break;
  }

  *iCs_end=-1;
  *viCs_end=-1;
  return iCs_end;
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
int* Mesh::nodeStar(Point const* p, int *iCs, int *niCs) const
{
  int const nrpc = this->numCornersPerCell();
  int const nfpc = this->numFacetsPerCell();
  int const nvpc = this->numVerticesPerCell();
  
  int C = p->getIncidCell();
  int nC= p->getPosition();

  FEPIC_CHECK(C>=0 && nC>=0, "invalid argument", std::runtime_error);

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


/** @brief Returns all vertices that are connected to a vertex.
 *  @param[in] p a pointer to the vertex.
 *  @param[out] iVs vector with the connected vertices.
 *  @return a pointer to the element following the end of the sequence iVs.
 */
int* Mesh::connectedVtcs(Point const* p, int *iVs) const
{
  int const nvpc = this->numVerticesPerCell();
  Cell const* cell;
  int id;
  int const* iVs_beg = iVs;
  int iCs[FEPIC_MAX_ICELLS];
  int viCs[FEPIC_MAX_ICELLS];

  const int n = static_cast<int>(this->vertexStar(p, iCs, viCs) - iCs);

  for (int ic = 0; ic < n; ++ic)
  {
    cell = this->getCellPtr(iCs[ic]);
    for (int j = 0; j < nvpc; ++j)
      if (j != viCs[ic])
      {
        id = cell->getNodeId(j);
        if (!checkValue(iVs_beg, static_cast<int const*>(iVs), id))
          *iVs++ = id;
      }
  }

  *iVs = -1;
  return iVs;
}


/** @brief Returns all vertices that are connected to a vertex, as well the incident cells.
 *  @param[in] p a pointer to the vertex.
 *  @param[out] iVs vector with the connected vertices.
 *  @param[out] iCs vector with the incident cells.
 *  @param[out] viCs vector with the p local-ids in each cell.
 *  @return a pointer to the element following the end of the sequence iVs.
 *  @note iVs[k] = getCellPtr(iCd[k])->getNodeId(viCs[k]);
 */
int* Mesh::connectedVtcs(Point const* p, int *iVs, int *iCs, int *viCs) const
{
  int const nvpc = this->numVerticesPerCell();
  Cell const* cell;
  int id;
  int const* iVs_beg = iVs;
  //int iCs[FEPIC_MAX_ICELLS];
  //int viCs[FEPIC_MAX_ICELLS];

  const int n = static_cast<int>(this->vertexStar(p, iCs, viCs) - iCs);

  for (int ic = 0; ic < n; ++ic)
  {
    cell = this->getCellPtr(iCs[ic]);
    for (int j = 0; j < nvpc; ++j)
      if (j != viCs[ic])
      {
        id = cell->getNodeId(j);
        if (!checkValue(iVs_beg, static_cast<int const*>(iVs), id))
          *iVs++ = id;
      }
  }

  *iVs = -1;
  return iVs;
}



/** @brief Returns all nodes that are connected to a node.
 *  @param[in] p a pointer to the node.
 *  @param[out] iNs vector with the connected nodes.
 *  @return a pointer to the element following the end of the sequence iNs.
 */
int* Mesh::connectedNodes(Point const* p, int *iNs) const
{
  int const nnpc = this->numNodesPerCell();
  Cell const* cell;
  int id;
  int const* iNs_beg = iNs;
  int iCs[FEPIC_MAX_ICELLS];
  int viCs[FEPIC_MAX_ICELLS];

  const int n = static_cast<int>(this->nodeStar(p, iCs, viCs) - iCs);

  for (int ic = 0; ic < n; ++ic)
  {
    cell = this->getCellPtr(iCs[ic]);
    for (int j = 0; j < nnpc; ++j)
      if (j != viCs[ic])
      {
        id = cell->getNodeId(j);
        if (!checkValue(iNs_beg, static_cast<int const*>(iNs), id))
          *iNs++ = id;
      }
  }

  *iNs = -1;
  return iNs;
}

// ----------------------------------------------------- INCID FACETS -----------------------------------------------

// TODO: implementar versão para 1d
int* Mesh::incidentFacets_1D(Point const* p, int *iFs, int *viFs) const
{
  // do nothing
  // annoying compiler
  p++; iFs++; viFs++;
  printf("not implemented yet\n");
  throw;
  return iFs;
}

int* Mesh::incidentFacets_2D(Point const* p, int *iFs, int *viFs) const
{
  //FEPIC_CHECK(unsigned(vC)<CT::n_vertices && C>=0, "invalid C or vC", std::invalid_argument);
  //vertexStar_Template(int C, int vC, int *iFs, int *viFs

  int const nfpc = this->numFacetsPerCell();
  Cell const* cell;
  int C = p->getIncidCell();
  int vC= p->getPosition();
  int g, D=C, vD=vC, q=0;
  int fnodes[32];
  int const nodeid = this->getPointId(p);

  if (!this->isVertex(p))
  {
    cell = this->getCellPtr(D);
    g = (nfpc + vD -q) % nfpc;
    *iFs++ = cell->getFacetId(g);
    this->getFacetNodesId(cell->getFacetId(g),fnodes);
    if (nodeid == fnodes[0])
      *viFs++ = 0;
    else if (nodeid == fnodes[1])
      *viFs++ = 1;
    else
    {
      FEPIC_CHECK(false, "implementation error, please contact any developer", std::runtime_error);
    }
  }
  else
  {
    for (;;)
    {
      cell = this->getCellPtr(D);
      g = (nfpc + vD -q)%nfpc;
      D = cell->getIncidCell(g);
      *iFs++ = cell->getFacetId(g);
      this->getFacetNodesId(cell->getFacetId(g),fnodes);
      if (nodeid == fnodes[0])
        *viFs++ = 0;
      else if (nodeid == fnodes[1])
        *viFs++ = 1;
      else
      {
        FEPIC_CHECK(false, "implementation error, please contact any developer", std::runtime_error);
      }
      if (D<0)
      {
        if(q) break;
        D = C;
        vD = vC;
        q = 1;
        continue;
      }
      if (D==C)
        break;
      g = cell->getIncidCellPos(g);
      vD = (g + 1 - q) % nfpc;
      //printf("DEBUG     %d, %d de %d, %d\n",D, vD, C, vC);
    }
  }
  *iFs = -1;
  *viFs = -1;
  return iFs;
}


int* Mesh::incidentFacets_3D(Point const* p, int *iFs, int *viFs) const
{
  // do nothing
  // annoying compiler
  p++; iFs++; viFs++;
  printf("not implemented yet\n");
  throw;
  return iFs;
}


// ---------------------------------------------------- EDGE STAR ---------------------------------------------------


int* Mesh::edgeStar_3D(int C, int eC, int *iCs, int *eiCs) const
{
  FEPIC_CHECK(C>=0 && eC>=0, "invalid argument", std::invalid_argument);

  int const nvpf = this->numVerticesPerFacet();
  int const nfpc = this->numFacetsPerCell();

  int const C_ini = C, eC_ini = eC;
  int D=C, eD=eC, q=0, s=0, f, g, anch, eCf, eDg;
  Cell const* cell = this->getCellPtr(C);

  *iCs++  = C_ini;
  *eiCs++ = eC_ini;

  for(;;)
  {
    f = cell->table_bC_x_fC(eC,q);
    D = cell->getIncidCell(f);
    //this->edgeStarNextCell_Template<cdim>(D, eD, q, D, eD, q);
    if (D==C_ini)
    {
      break;
    }
    if (D<0)
    {
      if(s==0)
      {
        cell = this->getCellPtr(C_ini);
        eC=eC_ini; q=1; s=1; continue;
      }
      else
        break;
    }
    g = cell->getIncidCellPos(f);
    anch = cell->getIncidCellAnch(f);
    eCf = cell->table_bC_x_fC(eC,q+2);
    eDg = (nfpc - eCf - anch) % nvpf;
    eD = cell->table_fC_x_bC(g,eDg);
    if(cell->table_bC_x_fC(eD,0)==g)
      q = 1;
    else
      q = 0;

    *iCs++  = D;
    *eiCs++ = eD;

    cell = this->getCellPtr(D);
    eC = eD;

  }

  *iCs = -1;
  *eiCs = -1;
  return iCs;
}


int* Mesh::edgeStar_2D(int C, int eC, int *iCs, int *eiCs) const
{
  FEPIC_CHECK(C>=0 && eC>=0, "invalid argument", std::invalid_argument);

  *iCs++  = C;
  *eiCs++ = eC;
  *iCs    = this->getCellPtr(C)->getIncidCell(eC);
  if (*iCs>=0)
  {
    ++iCs;
    *eiCs++ = this->getCellPtr(C)->getIncidCellPos(eC);
  }

  *iCs = -1;
  *eiCs = -1;
  return iCs;
}


// TODO: implementar uma versão para 1D
int* Mesh::edgeStar_1D(int C, int eC, int *iCs, int *eiCs) const
{
  FEPIC_ASSERT(false, "not implemented yes", std::runtime_error);

  // apenas para o compilador não encher o saco
  C++; eC++; iCs++; eiCs++;
  return 0;
}

  // ---------------------------------------------------- FACE STAR ---------------------------------------------------

int* Mesh::faceStar_3D(int C, int fC, int *iCs, int *fiCs) const
{
  *iCs++ = C;
  *fiCs++= fC;
  Cell const*const cell = this->getCellPtr(C);
  *iCs = cell->getIncidCell(fC);
  if (*iCs>=0)
  {
    ++iCs;
    *fiCs++ = cell->getIncidCellPos(fC);
  }
  return iCs;
}

int* Mesh::faceStar_2D(int C, int fC, int *iCs, int *fiCs) const
{
  *iCs++ = C;
  *fiCs = fC;
  return iCs;
}

int* Mesh::faceStar_1D(int C, int fC, int *iCs, int *fiCs) const
{
  FEPIC_ASSERT(true, "this function should not be called", std::invalid_argument);
  iCs++; fiCs++; fC++; C++;
  return iCs;
}


  // ------------------------------------------------- NEXT BOUNDARY FACET --------------------------------------------
Facet* Mesh::nextBoundaryFacet_1D(Facet const*)
{
  printf("nextBoundaryFacet not implemented for 1d cell yet.\n");
  throw;
  return NULL;
}

Facet* Mesh::nextBoundaryFacet_2D(Facet const* fct)
{
  int const nvpc = this->numVerticesPerCell();
  
  int iC = fct->getIncidCell();
  int iC_pos = fct->getPosition();
  int neighbor = this->getCellPtr(iC)->getIncidCell((iC_pos+1)%nvpc);
  
  FEPIC_CHECK(this->getCellPtr(iC)->getIncidCell(iC_pos) < 0, "nextBoundaryFacet: must be a boundary facet", std::invalid_argument);
  
  while (neighbor >= 0)
  {
    iC_pos = this->getCellPtr(iC)->getIncidCellPos((iC_pos+1)%nvpc);
    iC = neighbor;
    neighbor = this->getCellPtr(iC)->getIncidCell((iC_pos+1)%nvpc);
  }
  return this->getFacetPtr(this->getCellPtr(iC)->getFacetId((iC_pos+1)%nvpc));
}

Facet* Mesh::nextBoundaryFacet_3D(Facet const*)
{
  printf("nextBoundaryFacet not implemented for 3d cell yet.\n");
  throw;
  return NULL;
}


/** set iC as the incident cell to the point pt = (iC,pos). If an incident cell has the same
  *  connected component of the iC, then it is replaced.
  */ 
void Mesh::pushIncidCell2Point(Point *pt, int iC, int pos)
{
  FEPIC_CHECK(iC<this->numCellsTotal(),"invalid index iC", std::invalid_argument);
  
  Cell* icell = this->getCellPtr(iC);
  Point* p = this->getNodePtr(icell->getNodeId(pos));
  
  FEPIC_CHECK(p == pt, "incompatible argument:pt must match the (iC,pos)", std::invalid_argument);
  
  int const iC_CCid = icell->getConnectedComponentId();
  
  FEPIC_CHECK(iC_CCid>=0, "input iC has no connected component id", std::invalid_argument);
  
  const int nm_icells = p->numConnectedComps();
  
  int oic, opos;
  
  for (int i = 0; i < nm_icells; ++i)
  {
    p->getIthIncidCell(i,oic,opos);
    if (oic<0)
    {
      printf("i=%d iC_CCid=%d nm_icells=%d \n",i,iC_CCid, nm_icells);
    }
    if ((this->getCellPtr(oic))->getConnectedComponentId() == iC_CCid)
    {
      p->replacesIncidCell(oic, iC, pos);
      return;
    }
  }
  //else
  p->pushIncidCell(iC, pos);
  
}

// --------------------------------------------------- ADJACENCY -------------------------------------------------------



template <class Iter1, class Iter2, class T>
bool common1 (Iter1 first1, Iter1 last1,
              Iter2 first2, Iter2 last2,
              T & result)
{
  while (first1!=last1 && first2!=last2)
  {
    if (*first1<*first2) ++first1;
    else if (*first2<*first1) ++first2;
    else {
      result = *first1;
      return true;
    }
  }
  return false;
}

template <class Iter1, class Iter2, class Iter3, class T>
bool common2 (Iter1 first1, Iter1 last1,
              Iter2 first2, Iter2 last2,
              Iter3 first3, Iter3 last3,
              T & result)
{

  while (first1!=last1 && first2!=last2)
  {
    if (*first1<*first2) ++first1;
    else if (*first2<*first1) ++first2;
    else {
       if ( last3 != std::find(first3,last3,*first1) )
       {
         result = *first1;
         return true;
       }
       else
         ++first1;
    }
  }
  return false;

}

/// @note constroi as facets também
void Mesh::buildCellsAdjacency()
{
  const int cdim            = this->cellDim();
  const int n_vtx_per_facet = this->numVerticesPerFacet();
  const int n_facets        = this->numFacetsPerCell();
  //const int n_anch          = cdim ? n_vtx_per_facet : 1;

  //int const n_cells = this->numCells();
  int const n_cells_total = this->numCellsTotal();

  // vertices star list
  std::vector<SetVector<int> >  star( this->numNodesTotal() );

  // reserving memory
  int n_neibs = cdim<3 ? 10 : 40;
  for (unsigned i = 0; i < star.size(); ++i)
    star[i].reserve(n_neibs);

  Cell * cell;
  
  // constructing the vertices star list
  for (int i = 0; i < n_cells_total; ++i)
  {
    cell = getCellPtr(i);
    if (cell->isDisabled())
      continue;
    for (int j = 0; j < cell->numVertices(); ++j)
      star.at( cell->getNodeId(j) ).insert(i);
      
    for (int f = 0; f < n_facets; ++f)
      //cell->setFacetId(f, -1);
      cell->setIncidCell(f, -1);
  }

  //for (int i = 0; i < (int)star.size(); ++i)
  //{
  //  for (int j = 0; j < (int)star[i].size(); ++j)
  //  {
  //    std::cout << (*(star[i].data() + j)) << " ";
  //  }
  //  if (star[i].empty())
  //    std::cout << "--";
  //  std::cout << std::endl;
  //}

  

  std::vector<int> f_vtcs(n_vtx_per_facet);  
  std::tr1::shared_ptr<Facet> facet(this->createFacet());
  int facet_id;
  
  // building adjacency
  // also it creates the facets
  // the strategy: we use the information that for each facet there are
  //               at most two incident cells. The mate cell is given by
  //               the intersection of the vertices star lists.
  for (int i = 0; i < n_cells_total; ++i)
  {
    cell = getCellPtr(i);
    if (cell->isDisabled())
      continue;
    int neighbor_id;
    bool found;
    
    // do not count itself
    for (int k = 0; k < cell->numVertices(); ++k)
      star[ cell->getNodeId(k) ].erase(i);
    
    for (int f = 0; f < n_facets; ++f)
    {
      //if (cell->getFacetId(f) >= 0)
      if (cell->getIncidCell(f) >= 0)
        continue;
      cell->getFacetVerticesId(f, f_vtcs.data());

      if (cdim==1)
      {
        found = ! star[f_vtcs[0]].empty();
        if (found)
          neighbor_id = star[f_vtcs[0]].back();
      }
      else
      if (cdim == 2)
      { 
        int *beg1 = star[f_vtcs[0]].data();
        int *end1 = beg1 + star[f_vtcs[0]].size();
        int *beg2 = star[f_vtcs[1]].data();
        int *end2 = beg2 + star[f_vtcs[1]].size();
        
        found = common1(beg1, end1, beg2, end2, neighbor_id);
      }
      else
      if (cdim == 3)
      { 
        int *beg1 = star[f_vtcs[0]].data();
        int *end1 = beg1 + star[f_vtcs[0]].size();
        int *beg2 = star[f_vtcs[1]].data();
        int *end2 = beg2 + star[f_vtcs[1]].size();
        int *beg3 = star[f_vtcs[2]].data();
        int *end3 = beg3 + star[f_vtcs[2]].size();
        
        found = common2(beg1, end1, beg2, end2, beg3, end3, neighbor_id);
      }
      
      Cell * cellv;
      int anchor;
      int fv;
      if (found)
      {
        cellv = this->getCellPtr(neighbor_id);
        
        if ( cellv->isFacet(f_vtcs.data(), &fv, &anchor) )
        {
          fv = abs(fv);
          anchor = abs(anchor);
          
          cell->setIncidCell(f, neighbor_id);
          cell->setIncidCellPos(f, fv);

          cellv->setIncidCell(fv, i);
          cellv->setIncidCellPos(fv, f);

          if (cdim>2)
          {
            cell->setIncidCellAnch(f, anchor);
            cellv->setIncidCellAnch(fv, anchor);
          }          
        }
        else
        {
          FEPIC_CHECK(false, "something is wrong ...", std::runtime_error);
        }
        
      }
      
      if (cdim > 1)
      {
        // create the facet
        facet->setIncidCell(i);
        facet->setPosition(f);
        facet_id = this->pushFacet(facet.get());
        cell->setFacetId(f, facet_id);
        if (found)
          cellv->setFacetId(fv, facet_id);
      }
      
      //std::cout << "Cell: " << i << "; nodes: "<<f_vtcs[0]<<" "<<f_vtcs[1]<<" "<<f_vtcs[2]<< "; found = " << found;
      //if (found)
      // std::cout << ", who: " << neighbor_id;
      //std::cout << std::endl;
      
    } // end facets


  }



}


void Mesh::buildCorners_1D()
{
  m_corners_list.clear();
}

// assume that the vertices already have incident cells.
void Mesh::buildCorners_2D()
{
  m_corners_list.clear();
}

void Mesh::buildCorners_3D()
{

  int const num_cells  = this->numCellsTotal();
  int const nrpc       = this->numCornersPerCell();

  m_corners_list.clear();

  //// the CellList iterator must be "random access" for the algorithms that follows
  //// otherwise, the algorithms must be reimplementeds.
  //typedef typename CellIteratorT::iterator_category _category;
  //FEP_STATIC_ASSERT_ITERATOR((std::tr1::is_same<_category,std::random_access_iterator_tag>::value));

  FEP_PRAGMA_OMP(parallel for)
  for (unsigned c = 0; c<m_cells_list.totalSize(); ++c)
    m_cells_list[c].resetCorners();

  // TODO: think a way to parallelize this ...
  {

    Cell const* cell;
    //std::tr1::shared_ptr<CornerT> corner(new CornerT); // Corner
    boost::scoped_ptr<Corner> corner(this->createCorner());
    int corner_id;

    int iCs[FEPIC_MAX_ICELLS], eiCs[FEPIC_MAX_ICELLS]; // MAX_ICELLS
    int const* iCs_it, *eiCs_it;
    int const* iCs_end;

    for (int C = 0; C < num_cells; ++C)
    {
      cell = this->getCellPtr(C);
      if (cell->isDisabled())
        continue;

      for (int eC = 0; eC < nrpc; ++eC)
      {
        // se já foi vistitado pula
        if (cell->getCornerId(eC)>=0)
          continue;

        iCs_end = this->edgeStar(C, eC, iCs, eiCs);

        // create a edge
        corner->setIncidCell(C);
        corner->setPosition(eC);
        //corner->CornerT::setAnchor(0);
        corner_id = this->pushCorner(corner.get());

        for (iCs_it = iCs, eiCs_it = eiCs; iCs_it != iCs_end; ++iCs_it)
        {
          this->getCellPtr(*iCs_it)->setCornerId(*eiCs_it, corner_id);
          ++eiCs_it;
        }
      }
    }

  } // end parallel

}

/** atribui a cada nó, uma célula e sua posição
 *  @warning buildCellsAdjacency() must be called before this function.
 */
void Mesh::buildNodesAdjacency()
{
  int const num_cells    = this->numCellsTotal();
  int const nodes_p_cell = this->numNodesPerCell();
  int const numm_nodes    = this->numNodesTotal();
  int const nnpf         = this->numNodesPerFacet();
  int const cdim         = this->cellDim();
  int const nfpc         = this->numFacetsPerCell();

  FEP_PRAGMA_OMP(parallel for)
  for (int i=0; i<numm_nodes; ++i)
  {
    //m_points_list[i].setIncidCell(-1);
    //m_points_list[i].setPosition(-1);
    m_points_list[i].clearIncidences();
  }

  //FEP_PRAGMA_OMP(parallel default(none))
  {
    Cell const* cell;
    Point * point;
    //
    //// first the cells that is not in boundary
    ////FEP_PRAGMA_OMP(for)
    //for (int C = 0; C < num_cells; ++C)
    //{
    //  cell = this->getCellPtr(C);
    //  if (cell->isDisabled())
    //    continue;
    //
    //  for (int n = 0; n < nodes_p_cell; ++n)
    //  {
    //    point = this->getNodePtr(cell->getNodeId(n));
    //    //FEP_PRAGMA_OMP(critical)
    //    {
    //      point->PointT::setIncidCell(C);
    //      point->PointT::setPosition(n);
    //    }
    //  }
    //}

    std::vector<int> fnodes(nnpf);
    // then the cells that is in boundary
    //FEP_PRAGMA_OMP(for)
    for (int C = 0; C < num_cells; ++C)
    {
      cell = this->getCellPtr(C);
      if (cell->isDisabled())
        continue;

      for (int j = 0; j < nfpc; ++j)
      {
        cell->getFacetNodesId(j,fnodes.data());
        //if (cell->getIncidCell(j) >= 0)
        //  continue;
        for (int n = 0; n < nnpf; ++n)
        {
          point = this->getNodePtr(fnodes[n]);
          this->pushIncidCell2Point(point,C,cell->table_fC_x_nC(j,n));
          if (cell->getIncidCell(j) < 0)
            point->setAsBoundary(true);
          //if (cell->getIncidCell(j) < 0)
          //{
          //  this->pushIncidCell2Point(point,C,CT::m_table_fC_x_nC[j][n]);
          //  point->PointT::setAsBoundary(true);
          //}
          //else // if the point is not in boundary, it cant be singular, so
          //{
          //  point->PointT::setIncidCell(C);
          //  point->PointT::setPosition(CT::m_table_fC_x_nC[j][n]);
          //}
        }
      }
      
      // interior node, if exists
      if ((cdim==2 && this->cellHasFaceNodes()) || (cdim==3 && this->cellHasVolumeNodes()) || (cdim==1 && this->cellHasEdgeNodes()))
      {
        point = this->getNodePtr(cell->getNodeId(nodes_p_cell-1));
        point->setIncidCell(C);
        point->setPosition(nodes_p_cell-1);
      }
    }
    
  } // end parallel


}

void Mesh::fi_setConnectedComponentsId(cell_handler c_ini, int cc_id)
{
  //FEPIC_ASSERT((unsigned)this->getCellId(c_ini) < this->numCellsTotal(), " invalid pointer",std::invalid_argument );
  int const nfpc = this->numFacetsPerCell();
  std::deque<int> cells2setup; // TODO: usar deque<> talvez seja melhor
  int const n_cells_total = this->numCellsTotal();
  //std::list<int>::iterator current_id;
  Cell *oc, *current;

  cells2setup.push_back(c_ini.index());
  this->getCellPtr(cells2setup.front())->setVisitedTo(true);
  
  while (!cells2setup.empty())
  {
    FEPIC_ASSERT((int)cells2setup.size()<=n_cells_total, "Infinite loop at fi_setConnectedComponentsId", std::runtime_error);
    
    current = this->getCellPtr(cells2setup.front());
    
    for (int i = 0; i < nfpc; ++i)
    {
      if (current->getIncidCell(i) < 0)
        continue;
      oc = this->getCellPtr(current->getIncidCell(i));
      if (!oc->isVisited())
      {
        cells2setup.push_back(current->getIncidCell(i));
        oc->setVisitedTo(true);
      }
    }
    
    current->setConnectedComponentId(cc_id);
    cells2setup.pop_front();

  }
  
  //const int n_cells_total = m_cells_list.size();
  
  //clear flags
  FEP_PRAGMA_OMP(parallel for)
  for (int i = 0; i < n_cells_total; ++i)
  {
    this->getCellPtr(i)->setVisitedTo(false);
  }  
  
}

void Mesh::setUpConnectedComponentsId()
{
  const int n_cells_total = static_cast<int> (m_cells_list.size() );
  m_connected_comp_l.clear();
  
  //clear conn comp ids
  FEP_PRAGMA_OMP(parallel for)
  for (int i = 0; i < n_cells_total; ++i)
  {
    this->getCellPtr(i)->setConnectedComponentId(-1);
  }
  
  // start from 1
  int id = 1;
  
  //Cell *cellt;
  //cell_iterator cell     = this->cellBegin();
  //cell_iterator cell_end = this->cellEnd();
  //for (; cell != cell_end; ++cell)
  //{
  //  cellt = static_cast<CT*>(cell.getPtr());
  //  if ( cellt->getConnectedComponentId() >= 0)
  //    continue;
  //  this->fi_setConnectedComponentsId(cell,id);
  //   m_connected_comp_l.insert(std::pair<int,int>(id, cell.index()));
  //  ++id;
  //}
  
  cell_handler cell;
  for (int i = 0; i < n_cells_total; ++i)
  {
    cell = this->getCell(i);
    if ( !cell.isValid() || cell->getConnectedComponentId() >= 0)
      continue;
    this->fi_setConnectedComponentsId(cell,id);
     m_connected_comp_l.insert(std::pair<int,int>(id, cell.index()));
    ++id;
  }
  
  
}

void Mesh::fi_setBoundaryComponentsId(facet_handler f_ini, int bc_id)
{
  int const cdim = this->cellDim();
  if (cdim==1 || cdim==3)
  {
    //printf("warning: setBoundaryComponentsId not implemented for 1d and 3d cells.\n");
    return;
  }

  Facet* fit = nextBoundaryFacet(f_ini.ptr());
  f_ini->setBoundaryComponentId(bc_id);
  while (fit != f_ini.ptr() )
  {
    fit->setBoundaryComponentId(bc_id);
    fit = nextBoundaryFacet(fit);
    //printf("this->getFacetId(fit) = %d    this->getFacetId(f_ini) = %d \n", (int)this->getFacetId(fit), (int)this->getFacetId(f_ini));
  }
  
}


void Mesh::setUpBoundaryComponentsId()
{
  int const n_facets_total = this->numFacetsTotal();
  m_boundary_comp_l.clear();
    //clear conn comp ids
  FEP_PRAGMA_OMP(parallel for)
  for (int i = 0; i < n_facets_total; ++i)
  {
    this->getFacetPtr(i)->setBoundaryComponentId(-1);
  }
  // start from 1
  int id = 1;
  
  Facet* facett;
  facet_iterator facet = this->facetBegin();
  facet_iterator facet_end = this->facetEnd();
  for ( ;facet != facet_end; ++facet)
  {
    facett = &*facet;
    if (!this->inBoundary(facett) || facett->getBoundaryComponentId() >= 0 || facett->isDisabled())
      continue;
    this->fi_setBoundaryComponentsId(this->handler(facet),id);
     m_boundary_comp_l.insert(std::pair<int,int>(id, facet.index()));
    ++id;
  }  
  
  
}


/** Check if the vertices form a facet of this mesh, if so returns facet's id.
* @param vtcs the vertices.
* @return the id of the facet. If vertices do not form a facet, then function returns -1.
* @note vertices form a facet when they are cyclically equal to the facet's vertices.
*/
bool Mesh::getFacetIdFromVertices(int const* vtcs, int &fid)
{
  /* ideia: pega todas as células incidentes ao primeiro vértice
   * de vtcs, e procura pela facet em cada célula
   */

  Point const*const pt = this->getNodePtr(*vtcs);
  Cell const*      cell;
  int trash[FEPIC_MAX_ICELLS], iCs[FEPIC_MAX_ICELLS]; // MAX_ICELLS
  int fC;
  int anc;
  bool found;

  // células incidentes ao primeiro vértice da facet
  int* iCs_end = this->vertexStar(pt->getIncidCell(), pt->getPosition(), iCs, trash);

  for(int* iC = iCs; iC!=iCs_end; ++iC)
  {
    cell = this->getCellPtr(*iC);
    found = cell->isFacet(vtcs, &fC, &anc);

    if (!found) continue;

    if (fC<0) fid = -cell->getFacetId(abs(fC));
    else      fid = +cell->getFacetId(fC);
    return true;
  }
  return false;
}

/** Check if the vertices form a corner of this mesh, if so returns corner's id.
* @param vtcs the vertices.
* @return the id of the corner. If vertices do not form a corner, then function returns -1.
* @note vertices form a corner when they are cyclically equal to the corner's vertices.
*/
bool Mesh::getCornerIdFromVertices(int const* vtcs, int &rid)
{
  /* ideia: pega todas as células incidentes ao primeiro vértice
   * de vtcs, e procura pela corner em cada célula
   */
  
  int const cdim = this->cellDim();
  Point const*const pt = this->getNodePtr(*vtcs);
  Cell const*      cell;
  int trash[FEPIC_MAX_ICELLS], iCs[FEPIC_MAX_ICELLS]; // MAX_ICELLS
  int bfC;
  bool found;

  if (cdim==2) // corner is a point
  {
    if (!this->isVertex(pt))
      return -1;
    rid = this->getCellPtr(pt->getIncidCell())->getCornerId(pt->getPosition());
    return true;
  }

  // células incidentes ao primeiro vértice da corner
  int* iCs_end = this->vertexStar(pt->getIncidCell(), pt->getPosition(), iCs, trash);

  FEPIC_CHECK(iCs_end-iCs<=FEPIC_MAX_ICELLS,"INCREASE FEPIC_MAX_ICELLS ", std::runtime_error);

  for(int* iC = iCs; iC!=iCs_end; ++iC)
  {
    cell = this->getCellPtr(*iC);
    found = cell->isCorner(vtcs, bfC);

    if (!found) continue;

    if (bfC<0) rid = -cell->getCornerId(abs(bfC));
    else       rid = +cell->getCornerId(bfC);
    return true;
  }
  return false;
}


bool Mesh::inSingleCell(Point const* p) const
{
  int const nvpc = this->numVerticesPerCell();
  int const nrpc = this->numCornersPerCell();
  int const nfpc = this->numFacetsPerCell();
  int const cdim = this->cellDim();
  
  int iC  = p->getIncidCell();
  int viC = p->getPosition();

  Cell const* cell = this->getCellPtr(iC);

  if (viC<nvpc) // is a vertex
  {
    int counter = 0;
    for (int fv = 0; fv<cdim; ++fv) // for each incident facet
    {
      int f = cell->table_vC_x_fC( viC , fv );
      if (cell->getIncidCell(f) < 0)
        ++counter;
    }
    if (counter == cdim)
      return true;
    else
      return false;
  }
  else
  if (viC<nvpc + nrpc) // is an edge vertex
  {
    if (cdim == 1)
      return true;
    else
    if (cdim == 2)
      return cell->getIncidCell(viC - nvpc) < 0;
    else
    //if (cdim == 3)
    {
      int counter = 0;
      for (int fe = 0; fe < 2; ++fe) // for each incident facet
      {
        const int edge_lid = viC - nvpc;
        int f = cell->table_bC_x_fC( edge_lid , fe );
        if (cell->getIncidCell(f) < 0)
          ++counter;
      }
      if (counter == 2)
        return true;
      else
        return false;
    }
  
  }
  else
  if (viC<nvpc + nrpc + nfpc) // is an face vertex
  {
    if (cdim < 3)
      return true;

    int f = viC - nvpc - nrpc;
    return cell->getIncidCell(f) < 0;
  }
  else // is a volume point
  {
    return true;
  }

}

// only dim = 3
bool Mesh::inSingleCell(Corner const* r) const
{
  int const cdim = this->cellDim();
  if (cdim < 2)
    return false;

  int iC  = r->getIncidCell();
  int eiC = r->getPosition();

  Cell const* cell = this->getCellPtr(iC);

  if (cdim==2)
  {
    Point const* p = this->getNodePtr(cell->getNodeId(eiC));
    return this->inSingleCell(p);
  }

  int counter = 0;
  for (int fe = 0; fe < 2; ++fe) // for each incident facet
  {
    int f = cell->table_bC_x_fC( eiC , fe );
    if (cell->getIncidCell(f) < 0)
      ++counter;
  }
  if (counter == 2)
    return true;
  else
    return false;

}


bool Mesh::inSingleCell(Facet const* fa) const
{
  int iC  = fa->getIncidCell();
  int fiC = fa->getPosition();

  if (this->getCellPtr(iC)->getIncidCell(fiC)<0)
    return true;
  else
    return false;

}

// LINEAR CASE ONLY
void Mesh::getCenterCoord(Cell const* cell, Real* Xc) const
{
  int const sdim = this->spaceDim();
  int const nvpc = this->numVerticesPerCell();
  int i;
  for (i=0; i<sdim; ++i) Xc[i]=0;

  for (i = 0; i < nvpc; ++i)
    for (int j = 0; j < sdim; ++j)
      Xc[j] += getNodePtr(cell->getNodeId(i))->getCoord(j);

  for (i=0; i<sdim; ++i)
    Xc[i] /= nvpc;

}

void Mesh::getCenterCoord(Facet const* facet, Real* Xc) const
{
  int const nnpc = this->numNodesPerCell();
  int const nvpf = this->numVerticesPerFacet();
  int const sdim = spaceDim();
  
  std::vector<int> ids(nnpc);
  
  for (int i=0; i<sdim; ++i)
    Xc[i]=0;

  this->getFacetNodesId(facet,ids.data());

  for (int i = 0; i < nvpf; ++i)
    for (int j = 0; j < sdim; ++j)
      Xc[j] += this->getNodePtr( ids[i] )->getCoord(j);

  for (int i=0; i<sdim; ++i)
    Xc[i] /= nvpf;

}
void Mesh::getCenterCoord(Corner const* corner, Real* Xc) const
{
  int const cdim = this->cellDim();
  if (cdim<3)
  {
    Xc = NULL;
    return;
  }
  int const sdim = this->spaceDim();
  int const nnpr = this->numNodesPerCorner();
  int const nvpr = this->numVerticesPerCorner();
  
  std::vector<int> ids(nnpr + (cdim<3)); // (cdim<3) is to avoid compiler annoying

  for (int i=0; i<sdim; ++i)
    Xc[i]=0;

  this->getCornerNodesId(corner,ids.data());

  for (int i = 0; i < nvpr; ++i)
    for (int j = 0; j < sdim; ++j)
      Xc[j] += this->getNodePtr( ids[i] )->getCoord(j);

  for (int i=0; i<sdim; ++i)
    Xc[i] /= nvpr + (cdim==1?1:0);

}

Mesh* Mesh::create(ECellType type, int spacedim)
{

  if (static_cast<unsigned>(spacedim-1)>2)
    spacedim = ctypeDim(type);

  switch (type)
  {
    case EDGE2        : {return new Mesh(type, spacedim);}
    case EDGE3        : {return new Mesh(type, spacedim);}
    case TRIANGLE3    : {return new Mesh(type, spacedim);}
    case TRIANGLE6    : {return new Mesh(type, spacedim);}
    case QUADRANGLE4  : {return new Mesh(type, spacedim);}
    case QUADRANGLE8  : {return new Mesh(type, spacedim);}
    case QUADRANGLE9  : {return new Mesh(type, spacedim);}
    case TETRAHEDRON4 : {return new Mesh(type, spacedim);}
    case TETRAHEDRON10: {return new Mesh(type, spacedim);}
    case HEXAHEDRON8  : {return new Mesh(type, spacedim);}
    case HEXAHEDRON20 : {return new Mesh(type, spacedim);}
    case HEXAHEDRON27 : {return new Mesh(type, spacedim);}
    default:
    {
      FEPIC_CHECK(false, "invalid or not supported mesh type", std::invalid_argument);
      return NULL;
    }
  }

  //~ switch (type)
  //~ {
    //~ case EDGE2        : {return static_cast<Mesh*>(  new SMesh<Edge2>          (spacedim) );}
    //~ case EDGE3        : {return static_cast<Mesh*>(  new SMesh<Edge3>          (spacedim) );}
    //~ case TRIANGLE3    : {return static_cast<Mesh*>(  new SMesh<Triangle3>     (spacedim) );}
    //~ case TRIANGLE6    : {return static_cast<Mesh*>(  new SMesh<Triangle6>     (spacedim) );}
    //~ case QUADRANGLE4  : {return static_cast<Mesh*>(  new SMesh<Quadrangle4>   (spacedim) );}
    //~ case QUADRANGLE8  : {return static_cast<Mesh*>(  new SMesh<Quadrangle8>   (spacedim) );}
    //~ case QUADRANGLE9  : {return static_cast<Mesh*>(  new SMesh<Quadrangle9>   (spacedim) );}
    //~ case TETRAHEDRON4 : {return static_cast<Mesh*>(  new SMesh<Tetrahedron4>  (spacedim) );}
    //~ case TETRAHEDRON10: {return static_cast<Mesh*>(  new SMesh<Tetrahedron10> (spacedim) );}
    //~ case HEXAHEDRON8  : {return static_cast<Mesh*>(  new SMesh<Hexahedron8>   (spacedim) );}
    //~ case HEXAHEDRON20 : {return static_cast<Mesh*>(  new SMesh<Hexahedron20>  (spacedim) );}
    //~ case HEXAHEDRON27 : {return static_cast<Mesh*>(  new SMesh<Hexahedron27>  (spacedim) );}
    //~ default:
    //~ {
      //~ FEPIC_CHECK(false, "invalid or not supported mesh type", std::invalid_argument);
      //~ return NULL;
    //~ }
  //~ }
}












/* ========================================================
                                                          |
    ____       _       ____   _  __  _   _   ____         |
   | __ )     / \     / ___| | |/ / | | | | |  _ \        |
   |  _ \    / _ \   | |     | ' /  | | | | | |_) |       |
   | |_) |  / ___ \  | |___  | . \  | |_| | |  __/        |
   |____/  /_/   \_\  \____| |_|\_\  \___/  |_|           |
                                                          |
                                                          |
*/ //=====================================================



//   backup: 2013/09/03
//   /// @note constroi as facets também
//   void Mesh::buildCellsAdjacency()
//   {
//     const int cdim            = this->cellDim();
//     const int n_vtx_per_facet = this->numVerticesPerFacet();
//     const int n_facets        = this->numFacetsPerCell();
//     const int n_anch          = cdim ? n_vtx_per_facet : 1;
//   
//     //typedef std::tr1::array<int, n_vtx_per_facet> VecxT;
//     typedef std::vector<int>             VecxT;
//     typedef std::tr1::array<int, 2>      Vec2T;
//   
//     typedef std::pair<VecxT, Vec2T> PairT; // < (facet vertices) , (cell, ith) >
//     typedef std::vector<PairT>      MapT;
//     typedef MapT::const_iterator ConstIterT;
//   
//   
//     int const n_cells = this->numCells();
//     int const n_cells_total = this->numCellsTotal();
//   
//     MapT table(n_facets*n_cells);
//   
//     VecxT    facet_vtcs(n_vtx_per_facet);
//     Vec2T    cell_ith;
//   
//     m_facets_list.clear();
//   
//     // constroi uma tabela com as células e seus vizinhos
//     FEP_PRAGMA_OMP(parallel private(cell_ith,facet_vtcs) shared(table) default(none))
//     {
//       Cell const* cell;
//       int ii;
//       unsigned t;
//   
//       FEP_PRAGMA_OMP(for schedule (static) nowait)
//       for (int k = 0; k < n_cells_total; ++k)
//       {
//         cell = this->getCellPtr(k);
//         if (cell->isDisabled())
//           continue;
//   
//         ii = this->getCellContigId(k);
//         //ii = k;
//   
//         for (int j = 0; j < n_facets; ++j)
//         {
//           cell->getFacetVerticesId(j, facet_vtcs.data());
//           cell_ith[0] = k;
//           cell_ith[1] = j;
//           t = n_facets*ii + j;
//           //table[n_facets*ii + j] = std::make_pair(facet_vtcs, cell_ith);
//           table[t].first  = facet_vtcs;
//           table[t].second = cell_ith;
//   
//         }
//   
//       }
//   
//     }
//   
//     omptl::sort(table.begin(), table.end(), pair_less<PairT>());
//   
//     //// the CellList iterator must be "random access" for the algorithms that follows
//     //// otherwise, the algorithms must be reimplemented.
//     //typedef typename CellIteratorT::iterator_category _category;
//     //FEP_STATIC_ASSERT_ITERATOR((std::tr1::is_same<_category,std::random_access_iterator_tag>::value));
//   
//     // reseting cells neighbors
//     FEP_PRAGMA_OMP(parallel for) // WARNING: VALID ONLY FOR std::vector<> ....
//     for (int i=0; i<n_cells_total; ++i)
//     {
//       m_cells_list[i].resetIncidCells();
//       if (cdim>1)
//         m_cells_list[i].resetFacets();
//     }
//   
//     // build adjacency and create facets
//     FEP_PRAGMA_OMP(parallel private(facet_vtcs) shared(table) default(none))
//     {
//       int otherC, otherith, thisC, thisith;
//       int a;
//       std::tr1::shared_ptr<Facet> facet(this->createFacet());
//       int facet_id;
//   
//       bool found;
//       ConstIterT mid, table_end = table.end(), table_beg = table.begin();
//   
//       FEP_PRAGMA_OMP(for schedule (guided) nowait)
//       for (ConstIterT kit = table.begin(); kit < table_end; ++kit) // table loop
//       {
//   
//         std::reverse_copy(kit->first.begin(), kit->first.end(), facet_vtcs.begin());
//   
//         found = false;
//   
//         for (a = 0; a != n_anch; ++a) // ancora
//         {
//           mid = binary_find(table_beg, kit, facet_vtcs, pair_less<PairT>(), pair_eq<PairT>());
//   
//           if (mid != kit) // se econtrou uma face em comum
//           {
//   
//             otherC = mid->second[0];
//             otherith = mid->second[1];
//             thisC = kit->second[0];
//             thisith =  kit->second[1];
//   
//             (this->getCellPtr(thisC))->setIncidCell(thisith, otherC);
//             (this->getCellPtr(thisC))->setIncidCellPos(thisith, otherith);
//   
//             (this->getCellPtr(otherC))->setIncidCell(otherith, thisC);
//             (this->getCellPtr(otherC))->setIncidCellPos(otherith, thisith);
//   
//             if (cdim==3)
//             {
//               (this->getCellPtr(thisC))->setIncidCellAnch(thisith, a);
//               (this->getCellPtr(otherC))->setIncidCellAnch(otherith, a);
//             }
//   
//             found = true;
//             break;
//           }
//           std::rotate(facet_vtcs.begin(), facet_vtcs.begin()+1, facet_vtcs.end());
//         }
//         if (!found)
//         {
//           thisC = kit->second[0];
//           thisith =  kit->second[1];
//           if( cdim > 1) // border facet
//           {
//             // create a facet
//             facet->setIncidCell(thisC);
//             facet->setPosition(thisith);
//             //facet->FacetT::setAnchor(-1);
//             FEP_PRAGMA_OMP(critical)
//             facet_id = this->pushFacet(facet.get());
//             this->getCellPtr(thisC)->setFacetId(thisith, facet_id);
//           }
//         }
//   
//       }
//   
//   
//     } // end parallel
//   
//     // assigns facets to cells that remained
//     if (cdim > 1)
//     {
//       FEP_PRAGMA_OMP(parallel)
//       {
//         Cell * cell;
//         Cell const* icell;
//         int oth;
//         int facet_id;
//         FEP_PRAGMA_OMP(for)
//         for (int i=0; i<n_cells_total; ++i)
//         {
//           cell = getCellPtr(i);
//           if (cell->isDisabled())
//             continue;
//           for (int j = 0; j < n_facets; ++j)
//           {
//             if (cell->getFacetId(j) < 0)
//             {
//               icell = getCellPtr(cell->getIncidCell(j));
//               oth = cell->getIncidCellPos(j);
//               facet_id = icell->getFacetId(oth);
//               cell->setFacetId(j, facet_id);
//             }
//           }
//         }
//     
//       }
//     }
//   
//   
//   }




