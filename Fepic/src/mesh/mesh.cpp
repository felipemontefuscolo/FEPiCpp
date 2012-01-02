#include "mesh.hpp"

#include <vector>
#include <set>
#include "../util/misc2.hpp"
#include "../util/sorted_vec.hpp"
#include <map>
#include <cmath>
#if (FEP_HAS_OPENMP)
  #include <omp.h>
#endif

#if (FEP_HAS_OMPTL)
  #include <omptl/omptl_algorithm>
#else
  using namespace omptl = std;
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



Mesh::Mesh(int sd, ECellType fept)
{
  //_spacedim = sd; // BROKEN
  _cell_fep_tag = fept;
  _cell_msh_tag = ctype2mshTag(fept);
  _build_adjacency = true;

  timer = Timer();
};

/** @param nc_ number of cells.
 *  @param type mesh cell type.
 */
unsigned Mesh::estimateNumFacets(unsigned nc_, ECellType type)
{
  ECellClass cc = ctype2cclass(type);
  double nc = static_cast<float>(nc_);

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

};
unsigned Mesh::estimateNumCorners(unsigned nc_, ECellType type)
{
  ECellClass cc = ctype2cclass(type);
  double nc = static_cast<float>(nc_);

#define _ROUND_2_INT_(x) static_cast<unsigned>(floor((x) + 0.5))

  switch (cc)
  {
    case TETRAHEDRON:
      return _ROUND_2_INT_((7.*nc+9.*pow(6.,1./3)*pow(nc,2./3)+3*pow(6.,2./3)*pow(nc,1./3))/6.);
    case HEXAHEDRON:
      return _ROUND_2_INT_(3*nc+6*pow(nc,2./3)+3*pow(nc,1./3));
    default:
      return 0;
  }

#undef _ROUND_2_INT_
};



template<class CT, int SD> Cell*   SMesh<CT,SD>::incCell(Cell* a)         { return static_cast<Cell*>( &(*++(CellIteratorT(&_cellL, static_cast<CellT*>(a)))) );}
template<class CT, int SD> Point*  SMesh<CT,SD>::incPoint(Point* a)       { return static_cast<Point*>( &(*++(PointIteratorT(&_pointL, static_cast<PointT*>(a)))) );}
template<class CT, int SD> Facet*  SMesh<CT,SD>::incFacet(Facet* a)       { return static_cast<Facet*>( &(*++(FacetIteratorT(&_facetL, static_cast<FacetT*>(a)))) );}
template<class CT, int SD> Corner* SMesh<CT,SD>::incCorner(Corner* a)     { return static_cast<Corner*>( &(*++(CornerIteratorT(&_cornerL, static_cast<CornerT*>(a)))) );}
// -------------------------------------------- dec ------------------------
template<class CT, int SD> Cell*   SMesh<CT,SD>::decCell(Cell* a)         { return static_cast<Cell*>( &(*--(CellIteratorT(&_cellL, static_cast<CellT*>(a)))) );}
template<class CT, int SD> Point*  SMesh<CT,SD>::decPoint(Point* a)       { return static_cast<Point*>( &(*--(PointIteratorT(&_pointL, static_cast<PointT*>(a)))) );}
template<class CT, int SD> Facet*  SMesh<CT,SD>::decFacet(Facet* a)       { return static_cast<Facet*>( &(*--(FacetIteratorT(&_facetL, static_cast<FacetT*>(a)))) );}
template<class CT, int SD> Corner* SMesh<CT,SD>::decCorner(Corner* a)     { return static_cast<Corner*>( &(*--(CornerIteratorT(&_cornerL, static_cast<CornerT*>(a)))) );}
// -------------------------------------------- colors ------------------------
template<class CT, int SD> Cell*   SMesh<CT,SD>::incColorCell(Cell* a)    { return static_cast<Cell*>( &(*++(CellColorIteratorT(&_cellL, static_cast<CellT*>(a)))) );}
template<class CT, int SD> Point*  SMesh<CT,SD>::incColorPoint(Point* a)  { return static_cast<Point*>( &(*++(PointColorIteratorT(&_pointL, static_cast<PointT*>(a)))) );}
template<class CT, int SD> Facet*  SMesh<CT,SD>::incColorFacet(Facet* a)  { return static_cast<Facet*>( &(*++(FacetColorIteratorT(&_facetL, static_cast<FacetT*>(a)))) );}
template<class CT, int SD> Corner* SMesh<CT,SD>::incColorCorner(Corner* a){ return static_cast<Corner*>( &(*++(CornerColorIteratorT(&_cornerL, static_cast<CornerT*>(a)))) );}

// =====================================================================================
// =====================================================================================
//                                  ITERATORS
// =====================================================================================
// =====================================================================================



template<class CT,int SD> cell_iterator   SMesh<CT,SD>::cellBegin()  { return cell_iterator  (this, &(*_cellL.begin()));};
template<class CT,int SD> cell_iterator   SMesh<CT,SD>::cellEnd()    { return cell_iterator  (this, &(*_cellL.end()));};
template<class CT,int SD> point_iterator  SMesh<CT,SD>::pointBegin() { return point_iterator (this, &(*_pointL.begin()));};
template<class CT,int SD> point_iterator  SMesh<CT,SD>::pointEnd()   { return point_iterator (this, &(*_pointL.end()));};
template<class CT,int SD> facet_iterator  SMesh<CT,SD>::facetBegin() { return facet_iterator (this, &(*_facetL.begin()));};
template<class CT,int SD> facet_iterator  SMesh<CT,SD>::facetEnd()   { return facet_iterator (this, &(*_facetL.end()));};
template<class CT,int SD> corner_iterator SMesh<CT,SD>::cornerBegin(){ return corner_iterator(this, &(*_cornerL.begin()));};
template<class CT,int SD> corner_iterator SMesh<CT,SD>::cornerEnd()  { return corner_iterator(this, &(*_cornerL.end()));};

template<class CT,int SD> cell_iterator   SMesh<CT,SD>::cellBegin  (int tid, int nthreads) { return cell_iterator  (this, &(*_cellL.begin  (tid,nthreads)));};
template<class CT,int SD> cell_iterator   SMesh<CT,SD>::cellEnd    (int tid, int nthreads) { return cell_iterator  (this, &(*_cellL.end    (tid,nthreads)));};
template<class CT,int SD> point_iterator  SMesh<CT,SD>::pointBegin (int tid, int nthreads) { return point_iterator (this, &(*_pointL.begin (tid,nthreads)));};
template<class CT,int SD> point_iterator  SMesh<CT,SD>::pointEnd   (int tid, int nthreads) { return point_iterator (this, &(*_pointL.end   (tid,nthreads)));};
template<class CT,int SD> facet_iterator  SMesh<CT,SD>::facetBegin (int tid, int nthreads) { return facet_iterator (this, &(*_facetL.begin (tid,nthreads)));};
template<class CT,int SD> facet_iterator  SMesh<CT,SD>::facetEnd   (int tid, int nthreads) { return facet_iterator (this, &(*_facetL.end   (tid,nthreads)));};
template<class CT,int SD> corner_iterator SMesh<CT,SD>::cornerBegin(int tid, int nthreads) { return corner_iterator(this, &(*_cornerL.begin(tid,nthreads)));};
template<class CT,int SD> corner_iterator SMesh<CT,SD>::cornerEnd  (int tid, int nthreads) { return corner_iterator(this, &(*_cornerL.end  (tid,nthreads)));};

template<class CT,int SD> cell_color_iterator   SMesh<CT,SD>::cellBegin  (EColor k) { return cell_color_iterator  (this, &(*_cellL.begin  (k)));};
template<class CT,int SD> cell_color_iterator   SMesh<CT,SD>::cellEnd    (EColor k) { return cell_color_iterator  (this, &(*_cellL.end    (k)));};
template<class CT,int SD> point_color_iterator  SMesh<CT,SD>::pointBegin (EColor k) { return point_color_iterator (this, &(*_pointL.begin (k)));};
template<class CT,int SD> point_color_iterator  SMesh<CT,SD>::pointEnd   (EColor k) { return point_color_iterator (this, &(*_pointL.end   (k)));};
template<class CT,int SD> facet_color_iterator  SMesh<CT,SD>::facetBegin (EColor k) { return facet_color_iterator (this, &(*_facetL.begin (k)));};
template<class CT,int SD> facet_color_iterator  SMesh<CT,SD>::facetEnd   (EColor k) { return facet_color_iterator (this, &(*_facetL.end   (k)));};
template<class CT,int SD> corner_color_iterator SMesh<CT,SD>::cornerBegin(EColor k) { return corner_color_iterator(this, &(*_cornerL.begin(k)));};
template<class CT,int SD> corner_color_iterator SMesh<CT,SD>::cornerEnd  (EColor k) { return corner_color_iterator(this, &(*_cornerL.end  (k)));};

template<class CT,int SD> cell_color_iterator   SMesh<CT,SD>::cellBegin  (EColor k, int tid, int nthreads) { return cell_color_iterator  (this, &(*_cellL.begin  (k,tid,nthreads)));};
template<class CT,int SD> cell_color_iterator   SMesh<CT,SD>::cellEnd    (EColor k, int tid, int nthreads) { return cell_color_iterator  (this, &(*_cellL.end    (k,tid,nthreads)));};
template<class CT,int SD> point_color_iterator  SMesh<CT,SD>::pointBegin (EColor k, int tid, int nthreads) { return point_color_iterator (this, &(*_pointL.begin (k,tid,nthreads)));};
template<class CT,int SD> point_color_iterator  SMesh<CT,SD>::pointEnd   (EColor k, int tid, int nthreads) { return point_color_iterator (this, &(*_pointL.end   (k,tid,nthreads)));};
template<class CT,int SD> facet_color_iterator  SMesh<CT,SD>::facetBegin (EColor k, int tid, int nthreads) { return facet_color_iterator (this, &(*_facetL.begin (k,tid,nthreads)));};
template<class CT,int SD> facet_color_iterator  SMesh<CT,SD>::facetEnd   (EColor k, int tid, int nthreads) { return facet_color_iterator (this, &(*_facetL.end   (k,tid,nthreads)));};
template<class CT,int SD> corner_color_iterator SMesh<CT,SD>::cornerBegin(EColor k, int tid, int nthreads) { return corner_color_iterator(this, &(*_cornerL.begin(k,tid,nthreads)));};
template<class CT,int SD> corner_color_iterator SMesh<CT,SD>::cornerEnd  (EColor k, int tid, int nthreads) { return corner_color_iterator(this, &(*_cornerL.end  (k,tid,nthreads)));};





// =====================================================================================
// =====================================================================================
// =====================================================================================
// =====================================================================================

template<class CT, int SD>
int SMesh<CT,SD>::numVertices() const
{
  int num_vtcs = 0;
  int const num_nodes_total = MeshT::numNodesTotal();
  #pragma omp parallel shared(num_vtcs)
  {
    int num_vtcs_local = 0;
    PointT const* p;
    
    #pragma omp for
    for (int i=0; i<num_nodes_total; ++i)
    {
      p = MeshT::getNode(i);
      if (p->disabled())
        continue;
      if (MeshT::isVertex(p))
        ++num_vtcs_local;
    }
    
    #pragma omp critical
    num_vtcs += num_vtcs_local;
  }
  return num_vtcs;
}

/** Adiciona uma célula e retorna seu id.
*/
template<class CT, int SD>
int SMesh<CT,SD>::pushCell(Cell const* C)
{
  return _cellL.insert(*static_cast<CT const*>(C));
}

/** Adiciona um ponto e retorna seu id.
*/
template<class CT, int SD>
int SMesh<CT,SD>::pushPoint(Point const* P)
{
  return _pointL.insert(*static_cast<PointT const*>(P));
}

/** Adiciona uma facet
*  @param h A facet-xxxx a ser adicionada.
*  @return A posição da facet-xxxx na lista
*/
template<class CT, int SD>
int SMesh<CT,SD>::pushFacet(Facet const* F)
{
  return _facetL.insert(*static_cast<FacetT const*>(F));
}

/** Adiciona uma corner-xxxx
*  @param h A corner-xxxx a ser adicionada.
*  @return A posição da corner-xxxx na lista
*/
template<class CT, int SD>
int SMesh<CT,SD>::pushCorner(Corner const* B)
{
  return _cornerL.insert(*static_cast<CornerT const*>(B));
}


template<class CT, int SD>
void SMesh<CT,SD>::resizePointL(unsigned size)
{
  _pointL.resize(size);
}
template<class CT, int SD>
void SMesh<CT,SD>::resizeCellL(unsigned size)
{
  _cellL.resize(size);
}
template<class CT, int SD>
void SMesh<CT,SD>::resizeFacetL(unsigned size)
{
  _facetL.resize(size);
}
template<class CT, int SD>
void SMesh<CT,SD>::resizeCornerL(unsigned size)
{
  _cornerL.resize(size);
}


template<class CT, int SD>
void SMesh<CT,SD>::reservePointL(unsigned size)
{
  _pointL.resize(size);
}
template<class CT, int SD>
void SMesh<CT,SD>::reserveCellL(unsigned size)
{
  _cellL.resize(size);
}
template<class CT, int SD>
void SMesh<CT,SD>::reserveFacetL(unsigned size)
{
  _facetL.resize(size);
}
template<class CT, int SD>
void SMesh<CT,SD>::reserveCornerL(unsigned size)
{
  _cornerL.resize(size);
}

template<class CT, int SD>
void SMesh<CT,SD>::printInfo() const
{
  printf("elem type: %s\n",      nameForCtype(CellT::fep_tag));
  printf("space dim: %d\n",      SD                          );
  printf("# nodes:   %d\n",      numNodes()                  );
  printf("# cells:   %d\n",      numCells()                  );
  printf("# facets:  %d\n",      numFacets()                 );
  printf("# corners: %d\n",      numCorners()                );
  printf("# cell colors:  %d\n", numCellColors()             );
  printf("# node colors:  %d\n", numNodeColors()             );
}

template<class CT, int SD>
void SMesh<CT,SD>::printStatistics() const
{
  int cc = _cellL.capacity(),
      pc = _pointL.capacity(),
      fc = _facetL.capacity(),
      bc = _cornerL.capacity();
  printf("\n"
         "Cell list: \n"
         "  -allocated:      %d\n"
         "  -used:           %d (%d%%)\n"
         "  -# reservations: %d\n"
         "Node list: \n"
         "  -allocated:      %d\n"
         "  -used:           %d (%d%%)\n"
         "  -# reservations: %d\n"
         "Facet list: \n"
         "  -allocated:      %d\n"
         "  -used:           %d (%d%%)\n"
         "  -# reservations: %d\n"
         "Corner list: \n"
         "  -allocated:      %d\n"
         "  -used:           %d (%d%%)\n"
         "  -# reservations: %d\n",
         (int)_cellL.capacity(),
         (int)_cellL.total_size(), (int)(100*_cellL.total_size()/(cc?cc:1)),
         (int)_cellL.numReserveCalls(),
         (int)_pointL.capacity(),
         (int)_pointL.total_size(), (int)(100*_pointL.total_size()/(pc?pc:1)),
         (int)_pointL.numReserveCalls(),
         (int)_facetL.capacity(),
         (int)_facetL.total_size(), (int)(100*_facetL.total_size()/(fc?fc:1)),
         (int)_facetL.numReserveCalls(),
         (int)_cornerL.capacity(),
         (int)_cornerL.total_size(), (int)(100*_cornerL.total_size()/(bc?bc:1)),
         (int)_cornerL.numReserveCalls()
         );
}



// -------------------------------------------------- VERTEX STAR ---------------------------------------------------

// TODO: implementar versão para 1d
template<class CT, int SD>
template<int celldim>
int* SMesh<CT,SD>::vertexStar_Template(int C, int vC, int *iCs, int *viCs, typename EnableIf<(celldim==1)>::type*) const
{
  // do nothing
  // annoying compiler
  C++; vC++; iCs++; viCs++;
  return iCs;
}

template<class CT, int SD>
template<int celldim>
int* SMesh<CT,SD>::vertexStar_Template(int C, int vC, int *iCs, int *viCs, typename EnableIf<(celldim==2)>::type*) const
{
  FEPIC_CHECK(unsigned(vC)<CT::n_vertices && C>=0, "invalid C or vC", std::invalid_argument);

  CT const* cell;
  int g, D=C, vD=vC, q=0;

  *iCs++  = C;
  *viCs++ = vC;

  for (;;)
  {
    cell = this->MeshT::getCell(D);
    g = (CT::n_facets + vD -q)%CT::n_facets;
    D = cell->CT::getIncidCell(g);
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
    g = cell->CT::getIncidCellPos(g);
    vD = (g + 1 - q) % CT::n_facets;
    *iCs++  = D;
    *viCs++ = vD;
    //printf("DEBUG     %d, %d de %d, %d\n",D, vD, C, vC);
  }

  return iCs;
}

template<class CT, int SD>
template<int celldim>
int* SMesh<CT,SD>::vertexStar_Template(int C, int vC, int *iCs, int *viCs, typename EnableIf<(celldim==3)>::type*) const
{
  FEPIC_CHECK(unsigned(vC)<CT::n_vertices && C>=0, "invalid C or vC", std::invalid_argument);

  CT const*cell = this->MeshT::getCell(C);
  int f, D, vD, g, anc, vf, vg;
  int* iCs_end  = iCs;
  int* viCs_end = viCs;
  Uint_color_t visited_colors(0), one(1), color;
  //cell->CT::visited(true);

  *iCs_end++ = C;
  *viCs_end++ = vC;

  visited_colors |= one << cell->CT::getColor();

  for (;;)
  {
    // put the 3 neighbors in the stack
    for (int fv = 0; fv < 3; ++fv) // each vertice has three faces for both tetrahedron and hexahedron
    {
      // (C, f, vf, k) <===> (D, g, vg, h)

      f = CT::table_vC_x_fC[vC][fv];
      D = cell->CT::getIncidCell(f);
      //if (D<0 || this->MeshT::getCell(D)->CT::visited())

      if (D<0)
        continue;

      // if already in the stack
      if ((color = one << MeshT::getCell(D)->CT::getColor()) & visited_colors)
        continue;

      visited_colors |= color;

      g = cell->CT::getIncidCellPos(f);
      anc = cell->CT::getIncidCellAnch(f);

      vf = CT::table_vC_x_fC[vC][fv+3];
      vg = (CT::n_facets+1 - vf - anc)%CT::n_vertices_per_facet;
      vD  = CT::table_fC_x_vC[g][vg];

      *iCs_end++ = D;
      *viCs_end++ = vD;

      //this->MeshT::getCell(D)->CT::visited(true);

    }
    if (++iCs != iCs_end)
    {
      cell = MeshT::getCell(*iCs);
      vC = *++viCs;
    }
    else
      break;
  }

  //for (;iCs_beg != iCs_end; ++iCs_beg)
    //this->MeshT::getCell(*iCs_beg)->CellT::visited(false);

  return iCs_end;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
template<class CT, int SD>
template<int celldim>
int* SMesh<CT,SD>::vertexStar_WithoutColor(int C, int vC, int *iCs, int *viCs,
                                           typename EnableIf<(celldim<3)>::type*) const
{
  return vertexStar_Template<CT::dim>(C,vC,iCs,viCs);
}

template<class CT, int SD>
template<int celldim>
int* SMesh<CT,SD>::vertexStar_WithoutColor(int C, int vC, int *iCs, int *viCs,
                                           typename EnableIf<(celldim==3)>::type*) const
{
  FEPIC_CHECK(unsigned(vC)<CT::n_vertices && C>=0, "invalid C or vC", std::invalid_argument);

  CT const*cell = this->MeshT::getCell(C);
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

      f = CT::table_vC_x_fC[ vC ][ fv ];
      D = cell->CT::getIncidCell(f);

      if (D<0)
        continue;

      if ((iter = std::find((int*)iCs_beg, iCs_end, D)) != iCs_end)
      {
        if (iter > iCs)
        {
          g = cell->CT::getIncidCellPos(f);
          anc = cell->CT::getIncidCellAnch(f);
          vf = CT::table_vC_x_fC[ vC ][ fv+3 ];
          vg = (CT::n_facets+1 - vf - anc)%CT::n_vertices_per_facet;
          gv = CT::table_fC_x_vC[ g ][ vg + CT::n_vertices_per_facet];
          iCs_fv[(int)(iter-iCs_beg)] |= (1 << gv);
        }
        continue;
      }
      g = cell->CT::getIncidCellPos(f);
      anc = cell->CT::getIncidCellAnch(f);

      vf = CT::table_vC_x_fC[ vC ][ fv+3 ];
      vg = (CT::n_facets+1 - vf - anc)%CT::n_vertices_per_facet;
      vD  = CT::table_fC_x_vC[ g ][ vg ];

      *iCs_end++  = D;
      *viCs_end++ = vD;

      gv = CT::table_fC_x_vC[ g ][ vg + CT::n_vertices_per_facet];
      *fv_end++ = (1<<gv) ;

    }
    if (++iCs != iCs_end)
    {
      cell = this->MeshT::getCell(*iCs);
      vC = *++viCs;
    }
    else
      break;
  }

  //for (;iCs_beg != iCs_end; ++iCs_beg)
    //this->MeshT::getCell(*iCs_beg)->CellT::visited(false);

  return iCs_end;
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
template<class CT, int SD>
int* SMesh<CT,SD>::nodeStar(int C, int nC, int *iCs, int *niCs) const
{
  if (!CellT::has_edge_nodes || nC<CellT::n_vertices)
  {
    return this->MeshT::vertexStar_Template<CellT::dim>(C, nC, iCs, niCs);
  }

  else if (!CellT::has_face_nodes || nC<CellT::n_vertices + CellT::n_corners)
  {
    // number of nodes within the cell.
    int const eC = nC - CellT::n_vertices;
    int * iCs_end = this->edgeStar_Template<CellT::dim>(C, eC, iCs, niCs);
    *niCs++ = nC;
    ++iCs;
    for (; iCs!=iCs_end; ++iCs,++niCs)
      *niCs += CellT::n_vertices;
    return iCs_end;
  }

  else if (!CellT::has_volume_nodes || nC<CellT::n_vertices + CellT::n_corners + CellT::n_facets)
  {
    int const fC = nC - CellT::n_vertices - CellT::n_corners;
    int * iCs_end = this->MeshT::faceStar_Template<CellT::dim>(C, fC, iCs, niCs);
    *niCs++ = nC;
    ++iCs;
    for (; iCs!=iCs_end; ++iCs,++niCs)
      *niCs += CellT::n_vertices + CellT::n_corners;
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
template<class CT, int SD>
int* SMesh<CT,SD>::nodeStar(Point const* p, int *iCs, int *niCs) const
{
  int C = static_cast<PointT const*>(p)->getIncidCell();
  int nC= static_cast<PointT const*>(p)->getPosition();

  if (!CellT::has_edge_nodes || nC<CellT::n_vertices)
  {
    return this->MeshT::vertexStar_Template<CellT::dim>(C, nC, iCs, niCs);
  }

  else if (!CellT::has_face_nodes || nC<CellT::n_vertices + CellT::n_corners)
  {
    // number of nodes within the cell.
    int const eC = nC - CellT::n_vertices;
    int * iCs_end = this->edgeStar_Template<CellT::dim>(C, eC, iCs, niCs);
    *niCs++ = nC;
    ++iCs;
    for (; iCs!=iCs_end; ++iCs,++niCs)
      *niCs += CellT::n_vertices;
    return iCs_end;
  }

  else if (!CellT::has_volume_nodes || nC<CellT::n_vertices + CellT::n_corners + CellT::n_facets)
  {
    int const fC = nC - CellT::n_vertices - CellT::n_corners;
    int * iCs_end = this->MeshT::faceStar_Template<CellT::dim>(C, fC, iCs, niCs);
    *niCs++ = nC;
    ++iCs;
    for (; iCs!=iCs_end; ++iCs,++niCs)
      *niCs += CellT::n_vertices + CellT::n_corners;
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
template<class CT, int SD>
int* SMesh<CT,SD>::connectedVtcs(Point const* p, int *iVs) const
{
  CellT const* cell;
  int id;
  int const* iVs_beg = iVs;
  int iCs[FEPIC_MAX_ICELLS];
  int viCs[FEPIC_MAX_ICELLS];

  const int n = static_cast<int>(this->MeshT::vertexStar(p, iCs, viCs) - iCs);

  for (int ic = 0; ic < n; ++ic)
  {
    cell = this->MeshT::getCell(iCs[ic]);
    for (int j = 0; j < CellT::n_vertices; ++j)
      if (j != viCs[ic])
      {
        id = cell->CellT::getNodeId(j);
        if (!checkValue(iVs_beg, static_cast<int const*>(iVs), id))
          *iVs++ = id;
      }
  }
  return iVs;
}


/** @brief Returns all nodes that are connected to a node.
 *  @param[in] p a pointer to the node.
 *  @param[out] iNs vector with the connected nodes.
 *  @return a pointer to the element following the end of the sequence iNs.
 */
template<class CT, int SD>
int* SMesh<CT,SD>::connectedNodes(Point const* p, int *iNs) const
{
  CellT const* cell;
  int id;
  int const* iNs_beg = iNs;
  int iCs[FEPIC_MAX_ICELLS];
  int viCs[FEPIC_MAX_ICELLS];

  const int n = static_cast<int>(this->MeshT::nodeStar(p, iCs, viCs) - iCs);

  for (int ic = 0; ic < n; ++ic)
  {
    cell = this->MeshT::getCell(iCs[ic]);
    for (int j = 0; j < CellT::n_nodes; ++j)
      if (j != viCs[ic])
      {
        id = cell->CellT::getNodeId(j);
        if (!checkValue(iNs_beg, static_cast<int const*>(iNs), id))
          *iNs++ = id;
      }
  }
  return iNs;
}

// ---------------------------------------------------- EDGE STAR ---------------------------------------------------


template<class CT, int SD>
template<int celldim>
int* SMesh<CT,SD>::edgeStar_Template(int C, int eC, int *iCs, int *eiCs, typename EnableIf<(celldim==3)>::type*) const
{
  FEPIC_CHECK(C>=0 && eC>=0, "invalid argument", std::invalid_argument);

  int const C_ini = C, eC_ini = eC;
  int D=C, eD=eC, q=0, s=0, f, g, anch, eCf, eDg;
  CT const* cell = this->MeshT::getCell(C);

  *iCs++  = C_ini;
  *eiCs++ = eC_ini;

  for(;;)
  {
    f = CT::table_bC_x_fC[eC][q];
    D = cell->CT::getIncidCell(f);
    //this->MeshT::edgeStarNextCell_Template<CT::dim>(D, eD, q, D, eD, q);
    if (D==C_ini)
    {
      break;
    }
    if (D<0)
    {
      if(s==0)
      {
        cell = this->MeshT::getCell(C_ini);
        eC=eC_ini; q=1; s=1; continue;
      }
      else
        break;
    }
    g = cell->CT::getIncidCellPos(f);
    anch = cell->CT::getIncidCellAnch(f);
    eCf = CT::table_bC_x_fC[eC][q+2];
    eDg = (CT::n_facets - eCf - anch)%CT::n_vertices_per_facet;
    eD = CT::table_fC_x_bC[g][eDg];
    if(CT::table_bC_x_fC[eD][0]==g)
      q = 1;
    else
      q = 0;

    *iCs++  = D;
    *eiCs++ = eD;

    cell = this->MeshT::getCell(D);
    eC = eD;

  }

  return iCs;
}


template<class CT, int SD>
template<int celldim>
int* SMesh<CT,SD>::edgeStar_Template(int C, int eC, int *iCs, int *eiCs, typename EnableIf<(celldim==2)>::type*) const
{
  FEPIC_CHECK(C>=0 && eC>=0, "invalid argument", std::invalid_argument);

  *iCs++  = C;
  *eiCs++ = eC;
  *iCs    = this->MeshT::getCell(C)->CT::getIncidCell(eC);
  if (*iCs>=0)
  {
    ++iCs;
    *eiCs++ = this->MeshT::getCell(C)->CT::getIncidCellPos(eC);
  }

  return iCs;
}


// TODO: implementar uma versão para 1D
template<class CT, int SD>
template<int celldim>
int* SMesh<CT,SD>::edgeStar_Template(int C, int eC, int *iCs, int *eiCs, typename EnableIf<(celldim==1)>::type*) const
{
  FEPIC_ASSERT(false, "not implemented yes", std::runtime_error);

  // apenas para o compilador não encher o saco
  C++; eC++; iCs++; eiCs++;
  return 0;
}

  // ---------------------------------------------------- FACE STAR ---------------------------------------------------

template<class CT, int SD>
template<int celldim>
int* SMesh<CT,SD>::faceStar_Template(int C, int fC, int *iCs, int *fiCs, typename EnableIf<(celldim==3)>::type*) const
{
  *iCs++ = C;
  *fiCs++= fC;
  CellT const*const cell = this->MeshT::getCell(C);
  *iCs = cell->CellT::getIncidCell(fC);
  if (*iCs>=0)
  {
    ++iCs;
    *fiCs++ = cell->CT::getIncidCellPos(fC);
  }  
  return iCs;
}

template<class CT, int SD>
template<int celldim>
int* SMesh<CT,SD>::faceStar_Template(int C, int fC, int *iCs, int *fiCs, typename EnableIf<(celldim==2)>::type*) const
{
  *iCs++ = C;
  *fiCs = fC;
  return iCs;
}

template<class CT, int SD>
template<int celldim>
int* SMesh<CT,SD>::faceStar_Template(int C, int fC, int *iCs, int *fiCs, typename EnableIf<(celldim==1)>::type*) const
{
  FEPIC_ASSERT(true, "this function should not be called", std::invalid_argument);
  iCs++; fiCs++; fC++; C++;
  return iCs;
}

// --------------------------------------------------- ADJACENCY -------------------------------------------------------


/// Constrói as informações topológicas da malha
template<class CT, int SD>
void SMesh<CT,SD>::buildAdjacency()
{
  this->timer.restart();
  this->buildCellsAdjacency();
  this->timer.elapsed("buildCellsAdjacency()");

  this->timer.restart();
  this->buildNodesAdjacency();
  this->timer.elapsed("buildNodesAdjacency()");

  this->timer.restart();
  this->setColorAllCells();
  this->timer.elapsed("setColorAllCells()");
  
  this->timer.restart();
  this->setColorAllPoints();
  this->timer.elapsed("setColorAllPoints()");

  this->timer.restart();
  this->setColorAllFacets();
  this->timer.elapsed("setColorAllFacets()");
  
  // only for 3D
  this->timer.restart();
  this->buildCorners_Template<CellT::dim>();
  this->timer.elapsed("buildCorners_Template()");
  
  this->timer.restart();
  this->setColorAllCorners();
  this->timer.elapsed("setColorAllCorners()");  
}

/// @note constroi as facets também
template<class CT, int SD>
void SMesh<CT,SD>::buildCellsAdjacency()
{
  const int n_vtx_per_facet = CT::n_vertices_per_facet;
  const int n_facets        = CT::n_facets;
  const int n_anch          = CT::dim==3 ? n_vtx_per_facet : 1;

  typedef std::tr1::array<int, n_vtx_per_facet> VecxT;
  typedef std::tr1::array<int, 2>               Vec2T;

  typedef std::pair<VecxT, Vec2T> PairT; // < (facet vertices) , (cell, ith) >
  typedef std::vector<PairT>      MapT;
  typedef typename MapT::const_iterator ConstIterT;


  int const n_cells = this->MeshT::numCells();
  int const n_cells_total = this->MeshT::numCellsTotal();

  MapT table(n_facets*n_cells);

  VecxT    facet_vtcs;
  Vec2T    cell_ith;

  _facetL.clear();

  #pragma omp parallel private(cell_ith,facet_vtcs) shared(table) default(none)
  {
    CT const *cell;
    int ii;
    unsigned t;

    #pragma omp for schedule (static) nowait
    for (int k = 0; k < n_cells_total; ++k)
    {
      cell = static_cast<CT const*>(this->MeshT::getCell(k));
      if (cell->CT::disabled())
        continue;

      ii = this->MeshT::getCellContigId(k);
      //ii = k;

      for (int j = 0; j < n_facets; ++j)
      {
        cell->CT::getFacetVerticesId(j, facet_vtcs.data());
        cell_ith[0] = k;
        cell_ith[1] = j;
        t = n_facets*ii + j;
        //table[n_facets*ii + j] = std::make_pair(facet_vtcs, cell_ith);
        table[t].first  = facet_vtcs;
        table[t].second = cell_ith;
      }
    }

  }

  omptl::sort(table.begin(), table.end(), pair_less<PairT>());

  //// the CellList iterator must be "random access" for the algorithms that follows
  //// otherwise, the algorithms must be reimplementeds.
  //typedef typename CellIteratorT::iterator_category _category;
  //FEP_STATIC_ASSERT_ITERATOR((std::tr1::is_same<_category,std::random_access_iterator_tag>::value));

  // reseting cells neighbors
  #pragma omp parallel for     // WARNING: VALID ONLY FOR std::vector<> ....
  for (int i=0; i<n_cells_total; ++i)
  {
    _cellL[i].resetIncidCells();
    _cellL[i].resetFacetsId();
  }

  //// reseting cells neighbors
  //#pragma omp parallel for     // WARNING: VALID ONLY FOR std::vector<> ....
  //for (CellIteratorT cit = _cellL.begin(); cit < _cellL.end(); ++cit)
  //  cit->resetIncidCells();


  //
  // OBS: this version is a litte slower
  //

  //// build adjacency and create facets
  //#pragma omp parallel private(facet_vtcs) shared(table) default(none)
  //{
    //int otherC, otherith, thisC, thisith;
    //int a;
    //std::tr1::shared_ptr<FacetT> facet(new FacetT);
    //int facet_id;

    //bool found;
    //ConstIterT mid, table_end = table.end();

    //#pragma omp for schedule (dynamic) nowait
    //for (ConstIterT kit = table.begin(); kit < table_end; ++kit) // table loop
    //{

      //std::reverse_copy(kit->first.begin(), kit->first.end(), facet_vtcs.begin());

      //found = false;

      //for (a = 0; a != n_anch; ++a) // ancora
      //{
        //mid = binary_find(kit+1, table_end, facet_vtcs, pair_less<PairT>(), pair_eq<PairT>());

        //if (mid != table.end()) // se econtrou uma face em comum
        //{

          //otherC = mid->second[0];
          //otherith = mid->second[1];
          //thisC = kit->second[0];
          //thisith =  kit->second[1];


          //this->MeshT::getCell(thisC)->CT::setIncidCell(thisith, otherC);
          //this->MeshT::getCell(thisC)->CT::setIncidCellPos(thisith, otherith);

          //this->MeshT::getCell(otherC)->CT::setIncidCell(otherith, thisC);
          //this->MeshT::getCell(otherC)->CT::setIncidCellPos(otherith, thisith);

          //if (CT::dim==3)
          //{
            //this->MeshT::getCell(thisC)->CT::setIncidCellAnch(thisith, a);
            //this->MeshT::getCell(otherC)->CT::setIncidCellAnch(otherith, a);
          //}

          //found = true;
          //break;
        //}
        //std::rotate(facet_vtcs.begin(), facet_vtcs.begin()+1, facet_vtcs.end());
      //}
      //if(!found && CT::dim > 1) // border facet
      //{

        //thisC = kit->second[0];
        //thisith =  kit->second[1];

        //// create a facet
        //facet->FacetT::setIncidCell(thisC);
        //facet->FacetT::setPosition(thisith);
        //facet->FacetT::setAnchor(-1);
        //#pragma omp critical
        //facet_id = this->MeshT::pushFacet(facet.get());
        //this->MeshT::getCell(thisC)->CT::setFacetId(thisith, facet_id);

      //}

    //}

  //} // end parallel



  //
  // OBS: this version is a litte faster
  //

  // build adjacency and create facets
  #pragma omp parallel private(facet_vtcs) shared(table) default(none)
  {
    int otherC, otherith, thisC, thisith;
    int a;
    std::tr1::shared_ptr<FacetT> facet(new FacetT);
    int facet_id;

    bool found;
    ConstIterT mid, table_end = table.end(), table_beg = table.begin();

    #pragma omp for schedule (guided) nowait
    for (ConstIterT kit = table.begin(); kit < table_end; ++kit) // table loop
    {

      std::reverse_copy(kit->first.begin(), kit->first.end(), facet_vtcs.begin());

      found = false;

      for (a = 0; a != n_anch; ++a) // ancora
      {
        mid = binary_find(table_beg, kit, facet_vtcs, pair_less<PairT>(), pair_eq<PairT>());

        if (mid != kit) // se econtrou uma face em comum
        {

          otherC = mid->second[0];
          otherith = mid->second[1];
          thisC = kit->second[0];
          thisith =  kit->second[1];


          this->MeshT::getCell(thisC)->CT::setIncidCell(thisith, otherC);
          this->MeshT::getCell(thisC)->CT::setIncidCellPos(thisith, otherith);

          this->MeshT::getCell(otherC)->CT::setIncidCell(otherith, thisC);
          this->MeshT::getCell(otherC)->CT::setIncidCellPos(otherith, thisith);

          if (CT::dim==3)
          {
            this->MeshT::getCell(thisC)->CT::setIncidCellAnch(thisith, a);
            this->MeshT::getCell(otherC)->CT::setIncidCellAnch(otherith, a);
          }

          found = true;
          break;
        }
        std::rotate(facet_vtcs.begin(), facet_vtcs.begin()+1, facet_vtcs.end());
      }
      if(!found && CT::dim > 1) // border facet
      {

        thisC = kit->second[0];
        thisith =  kit->second[1];

        // create a facet
        facet->FacetT::setIncidCell(thisC);
        facet->FacetT::setPosition(thisith);
        facet->FacetT::setAnchor(-1);
        #pragma omp critical
        facet_id = this->MeshT::pushFacet(facet.get());
        this->MeshT::getCell(thisC)->CT::setFacetId(thisith, facet_id);
      }

    }

  } // end parallel

  // assigns facets to cells that remained
  #pragma omp parallel
  {
    CT * cell;
    CT const* icell;
    int oth;
    int facet_id;
    #pragma omp for
    for (int i=0; i<n_cells_total; ++i)
    {
      cell = MeshT::getCell(i);
      if (cell->disabled())
        continue;
      for (int j = 0; j < CT::n_facets; ++j)
      {
        if (cell->CT::getFacetId(j) < 0)
        {
          icell = MeshT::getCell(cell->CT::getIncidCell(j));
          oth = cell->CT::getIncidCellPos(j);
          facet_id = icell->CT::getFacetId(oth);
          cell->CT::setFacetId(j, facet_id);
        }
      }
    }
    
  }

}


// para dim = 1 e 2
template<class CT, int SD>
template<int celldim>
void SMesh<CT,SD>::buildCorners_Template(typename EnableIf<(celldim<3)>::type*)
{
  // nada
}

template<class CT, int SD>
template<int celldim>
void SMesh<CT,SD>::buildCorners_Template(typename EnableIf<(celldim==3)>::type*)
{

  int const num_cells  = this->MeshT::numCellsTotal();
  int const n_corners  = CT::n_corners;

  _cornerL.clear();

  //// the CellList iterator must be "random access" for the algorithms that follows
  //// otherwise, the algorithms must be reimplementeds.
  //typedef typename CellIteratorT::iterator_category _category;
  //FEP_STATIC_ASSERT_ITERATOR((std::tr1::is_same<_category,std::random_access_iterator_tag>::value));

  #pragma omp parallel for
  for (unsigned c = 0; c<_cellL.total_size(); ++c)
    _cellL[c].CT::resetCornersId();

  //#pragma omp parallel for
  //for (CellIteratorT cit = _cellL.begin(); cit < _cellL.end(); ++cit)
  //  cit->CT::resetCornersId();


  // think a way to parallelize this ...
  {

    CT const* cell;
    //std::tr1::shared_ptr<CornerT> corner(new CornerT); // Corner
    CornerT *corner = new CornerT;
    int corner_id;

    int iCs[FEPIC_MAX_ICELLS], eiCs[FEPIC_MAX_ICELLS]; // MAX_ICELLS
    int const* iCs_it, *eiCs_it;
    int const* iCs_end;

    for (int C = 0; C < num_cells; ++C)
    {
      cell = this->MeshT::getCell(C);
      if (cell->CT::disabled())
        continue;


      for (int eC = 0; eC < n_corners; ++eC)
      {
        // se já foi vistitado pula
        if (cell->CT::getCornerId(eC)>=0)
          continue;

        iCs_end = this->MeshT::edgeStar(C, eC, iCs, eiCs);

        // create a edge
        corner->CornerT::setIncidCell(C);
        corner->CornerT::setPosition(eC);
        corner->CornerT::setAnchor(0);
        corner_id = this->MeshT::pushCorner(corner);

        for (iCs_it = iCs, eiCs_it = eiCs; iCs_it != iCs_end; ++iCs_it)
        {
          this->MeshT::getCell(*iCs_it)->CT::setCornerId(*eiCs_it, corner_id);
          ++eiCs_it;
        }

      }


    }

    delete corner;

  } // end parallel

}

/** atribui a cada nó, uma célula e sua posição
 */
template<class CT, int SD>
void SMesh<CT,SD>::buildNodesAdjacency()
{
  int const num_cells    = this->numCellsTotal();
  int const nodes_p_cell = CT::n_nodes;

  //#pragma omp parallel default(none)
  {
    CT const* cell;
    PointT * point;

    // first the cells that is not in boundary
    //#pragma omp for
    for (int C = 0; C < num_cells; ++C)
    {
      cell = this->MeshT::getCell(C);
      if (cell->CT::disabled())
        continue;

      if(!cell->CT::inBoundary())
        for (int n = 0; n < nodes_p_cell; ++n)
        {
          point = this->MeshT::getNode(cell->CT::getNodeId(n));
          //#pragma omp critical
          {
          point->PointT::setIncidCell(C);
          point->PointT::setPosition(n);
          }
        }
    }


    // then the cells that is in boundary
    //#pragma omp for
    for (int C = 0; C < num_cells; ++C)
    {
      cell = this->MeshT::getCell(C);
      if (cell->CT::disabled())
        continue;

      if(cell->CT::inBoundary())
        for (int n = 0; n < nodes_p_cell; ++n)
        {
          point = this->MeshT::getNode(cell->CT::getNodeId(n));
          //#pragma omp critical
          {
          point->PointT::setIncidCell(C);
          point->PointT::setPosition(n);
          }
        }
    }

  } // end parallel

}



/** Check if the vertices form a facet of this mesh, if so returns facet's id.
* @param vtcs the vertices.
* @return the id of the facet. If vertices do not form a facet, then function returns -1.
* @note vertices form a facet when they are cyclically equal to the facet's vertices.
*/
template<class CT, int SD>
int SMesh<CT,SD>::getFacetIdFromVertices(int const* vtcs)
{
  /* ideia: pega todas as células incidentes ao primeiro vértice
   * de vtcs, e procura pela facet em cada célula
   */

  PointT const*const pt = this->MeshT::getNode(*vtcs);
  CellT  const*      cell;
  int trash[FEPIC_MAX_ICELLS], iCs[FEPIC_MAX_ICELLS]; // MAX_ICELLS
  int fC;

  // células incidentes ao primeiro vértice da facet
  int* iCs_end = this->MeshT::vertexStar(pt->PointT::getIncidCell(), pt->PointT::getPosition(), iCs, trash);

  for(int* iC = iCs; iC!=iCs_end; ++iC)
  {
    cell = this->MeshT::getCell(*iC);
    fC = cell->CellT::isFacet(vtcs);

    if (fC<0) continue;

    return cell->CellT::getFacetId(fC);
  }
  return -1;
}

/** Check if the vertices form a corner of this mesh, if so returns corner's id.
* @param vtcs the vertices.
* @return the id of the corner. If vertices do not form a corner, then function returns -1.
* @note vertices form a corner when they are cyclically equal to the corner's vertices.
*/
template<class CT, int SD>
int SMesh<CT,SD>::getCornerIdFromVertices(int const* vtcs)
{
  /* ideia: pega todas as células incidentes ao primeiro vértice
   * de vtcs, e procura pela corner em cada célula
   */

  PointT const*const pt = this->MeshT::getNode(*vtcs);
  CellT  const*      cell;
  int trash[FEPIC_MAX_ICELLS], iCs[FEPIC_MAX_ICELLS]; // MAX_ICELLS
  int bfC;

  // células incidentes ao primeiro vértice da corner
  int* iCs_end = this->MeshT::vertexStar(pt->PointT::getIncidCell(), pt->PointT::getPosition(), iCs, trash);

  for(int* iC = iCs; iC!=iCs_end; ++iC)
  {
    cell = this->MeshT::getCell(*iC);
    bfC = cell->CellT::isCorner(vtcs);

    if (bfC<0) continue;

    return cell->CellT::getCornerId(bfC);
  }
  return -1;
}


/** @brief [DEVELOPER] Assign to the cell the minimal color that is not in the neighborhood.
 *  @param cell pointer to the cell.
 *  @param neibs_beg begin of a vector with the neighbors of the cell.
 *  @param neibs_end end of the vector.
 *  @note Faster version. Update _n_cell_colors attribute of the mesh.
 *  @note if neibs includes the cell itself, then the cell must have color -1.
 */
template<class CT, int SD>
void SMesh<CT,SD>::setColorCell(Cell* cell_, int const* neibs_beg, int const* neibs_end)
{
  CT * cell = static_cast<CT*>(cell_);

  // all colors in compact form
  Uint_color_t all_colors(0);
  Uint_color_t const one(1);
  char color;

  // compacting colors ...
  for (int const*it = neibs_beg; it != neibs_end; ++it)
    if ((color = MeshT::getCell(*it)->CT::getColor()) >= 0)
      all_colors |=  one << color ;


  char new_color;
  for (Uint_color_t c = 1; ; c<<=1)
  {
    FEPIC_CHECK(c, "infinite loop", std::runtime_error);
    // if not found in the range
    if ( c & all_colors )
      continue;

    cell->CT::setColor( new_color = static_cast<char>( log2_i128(c) )  );
    break;
  }


  //if (new_color >= _n_cell_colors)
  //{
  //  _n_cell_colors = new_color;
  //}
}

/** @brief choose a color that is not in the neighborhood and assigns to the cell.
 *  @param id cell (id) to assign a color.
 *  @warning only for colored meshes.
 */
template<class CT, int SD>
void SMesh<CT,SD>::setColorCell(int id)
{
  int neibs[FEPIC_MAX_ICELLS*4/(CT::dim==3 ? 1 : 2)];
  int vops [FEPIC_MAX_ICELLS];

  CT *cell = this->MeshT::getCell(id);
  cell->CT::setColor(-1);

  int *mid = neibs;

  // contém informações redundantes
  for (int i = 0; i < CT::n_vertices; ++i)
    mid = this->MeshT::vertexStar(id,i,mid,vops);

  this->MeshT::setColorCell(cell,neibs,mid);
}

/** @brief check if a cell has invalid color.
 *  @param color the color of the cell.
 *  @param neibs_beg begining of a range that contain the cell's neighbors.
 *  @param neibs_end end of the range.
 *  @note all neighbors must be colored!
 *  @return if cell's color is illegal.
 */
template<class CT, int SD>
bool SMesh<CT,SD>::checkInvalidColorCell(Cell * cell_, int const* neibs_beg, int const* neibs_end) const
{
  CT * cell = static_cast<CT*>(cell_);
  char const this_color = cell->CT::getColor();

  if (this_color < 0)
    return true;

  // all colors in compact form
  Uint_color_t all_colors(0);
  Uint_color_t const one(1);

  // compacting colors ...
  for (int const*it = neibs_beg; it != neibs_end; ++it)
    if ((*it) >= 0)
      all_colors |= ( one << (this->MeshT::getCell(*it)->CT::getColor()) );

  return (one << this_color) & all_colors ;
}

/** @brief cells coloring.
 */
template<class CT, int SD>
void SMesh<CT,SD>::setColorAllCells()
{
  /*****************************************************
   *                 Coloring Cells
   *****************************************************/


  const int n_cells = this->numCellsTotal();

  // tab_neibs = table with neighbors of each cell
  std::vector<std::vector<int> > tab_neibs(n_cells, std::vector<int>(FEPIC_MAX_ICELLS*4/(CT::dim==3 ? 1 : 2)));

  // building tab_neibs
  #pragma omp parallel shared(tab_neibs) default(none)
  {
    CT * cell;
    int *mid;
    int vops[FEPIC_MAX_ICELLS];

    #pragma omp for
    for (int c = 0; c < n_cells; ++c)
    {
      cell = this->MeshT::getCell(c);
      if (cell->CT::disabled())
        continue;

      mid = tab_neibs[c].data();
      for (int i = 0; i < CT::n_vertices; ++i)
        mid = vertexStarWithoutColor(c,i,mid,vops);

      tab_neibs[c].resize(mid - tab_neibs[c].data());

      cell->CT::setColor(-1);
    }
  }


  CT * cell;
  // consume ~ 28% cpu
  for (int c = 0; c < n_cells; ++c)
  {
    cell = this->getCell(c);
    if (cell->disabled())
      continue;

    this->MeshT::setColorCell(cell, tab_neibs[c].data(), tab_neibs[c].data() + tab_neibs[c].size());
  }
  
  _cellL.linkColors();
  
  //// DO NOT ERASE ME ... I'M GOOD FOR TESTING.
  //Cell * cell;

  //for (int c = 0; c < this->numCellsTotal(); ++c)
  //{

    //cell = this->getCell(c);
    //if (cell->disabled())
      //continue;

    //this->setColorCell(c);
  //}


}

/** @brief cells coloring.
 *  @warning must be called after setColorAllCells()
 */
template<class CT, int SD>
void SMesh<CT,SD>::setColorAllPoints()
{

  /*****************************************************
   *                 Coloring Nodes
   *****************************************************/


  const int n_nodes_total = this->numNodesTotal();
  //const int n_cells = this->numCellsTotal();
  //const int n_c_colors = this->numCellColors();


  // tab_neibs = table with neighbors of each node
  std::vector<std::set<int> > tab_neibs(n_nodes_total);

  #pragma omp parallel for
  for (int i = 0; i < n_nodes_total; ++i)
    MeshT::getNode(i)->PointT::setColor(-1);

  //for (int k = 0; k < n_c_colors; ++k)
  //{
  //  #pragma omp parallel shared(tab_neibs) firstprivate(k) default(none)
  //  {
  //    CT const* cell;
  //    int cnodes[CT::n_nodes];
  //
  //    #pragma omp for nowait
  //    for (int c = 0; c < n_cells; ++c)
  //    {
  //      if ((cell = MeshT::getCell(c))->CT::getColor() != k)
  //        continue;
  //      if (cell->CT::disabled())
  //        continue;
  //
  //      cell->CT::getNodesId(cnodes);
  //
  //      for (int i = 0; i < CT::n_nodes; ++i)
  //        for (int j = 0; j < CT::n_nodes; ++j)
  //          tab_neibs[cnodes[i]].insert(cnodes[j]);
  //    }
  //  }
  //}

  #pragma omp parallel
  {
    const int tid      = omp_get_thread_num();
    const int nthreads = omp_get_num_threads();
    const int n_cell_colors = MeshT::numCellColors();
    int cnodes[CT::n_nodes];
    
    for (int k = 0; k < n_cell_colors; ++k)
    {
      // TODO: try to use const iterator. :error: no match for ‘operator!=’ in ‘cell != cell_end’
      CellColorIteratorT cell_end = _cellL.end  (EColor(k), tid, nthreads);
      CellColorIteratorT cell     = _cellL.begin(EColor(k), tid, nthreads);
      
      for (; cell != cell_end; ++cell)
      {
        cell->CT::getNodesId(cnodes);
        
        for (int i = 0; i < CT::n_nodes; ++i)
          for (int j = 0; j < CT::n_nodes; ++j)
            tab_neibs[cnodes[i]].insert(cnodes[j]);
      }
      #pragma omp barrier
    }
  } // end parallel

  #pragma omp parallel for shared(tab_neibs)
  for (int i = 0; i < n_nodes_total; ++i)
    tab_neibs[i].erase(i);


  //// the NodeList iterator must be "random access" for the algorithms that follows
  //// otherwise, the algorithms must be reimplemented.
  //typedef typename PointIteratorT::iterator_category _category;
  //FEP_STATIC_ASSERT_ITERATOR((std::tr1::is_same<_category,std::random_access_iterator_tag>::value));


  std::vector<std::set<int> >::iterator it;
  std::set<int>::iterator at;
  int kk = 0;
  PointT * point;
  int max_color(-1);
  
  for (it = tab_neibs.begin(); it!=tab_neibs.end(); ++it, ++kk)
  {
    point = MeshT::getNode(kk);
    if (point->disabled())
      continue;

    Uint_color_t all_colors(0);
    Uint_color_t const one(1);
    char color;
    
    for (at = (*it).begin(); at != (*it).end(); ++at)
      if ((color = MeshT::getNode(*at)->PointT::getColor()) >= 0)
        all_colors |= one << color;
    
    // find a suitable color
    for (Uint_color_t c = 1;; c<<=1)
    {
      if ( c & all_colors )
        continue;
      point->PointT::setColor(color=log2_i128(c));
      if (max_color < color)
        max_color = color;
      break;
    }
    
  }
//  _n_node_colors = max_color+1;

  _pointL.linkColors();

}


/** @brief facets coloring.
 *  @warning must be called after setColorAllCells()
 */
template<class CT, int SD>
void SMesh<CT,SD>::setColorAllFacets()
{

  /*****************************************************
   *                 Coloring Facets
   *****************************************************/

  const int n_facets_total = this->numFacetsTotal();

  #pragma omp parallel
  {
    FacetT * a;
    CellT const* ic;
    #pragma omp for
    for (int i=0; i<n_facets_total; ++i)
    {
      a = MeshT::getFacet(i);
      if (a->FacetT::disabled())
        continue;
      ic = MeshT::getCell(a->FacetT::getIncidCell());
      a->FacetT::setColor(ic->CT::getColor());
    };
  }

  _facetL.linkColors();
}

/** @brief corners coloring.
 *  @warning must be called after setColorAllCells()
 */
template<class CT, int SD>
void SMesh<CT,SD>::setColorAllCorners()
{

  /*****************************************************
   *                 Coloring Facets
   *****************************************************/

  const int n_corners_total = this->numCornersTotal();

  #pragma omp parallel
  {
    CornerT * a;
    CellT const* ic;
    #pragma omp for
    for (int i=0; i<n_corners_total; ++i)
    {
      a = MeshT::getCorner(i);
      if (a->CornerT::disabled())
        continue;
      ic = MeshT::getCell(a->CornerT::getIncidCell());
      a->CornerT::setColor(ic->CT::getColor());
    };
  }

  _cornerL.linkColors();
}

template<class CT, int SD>
bool SMesh<CT,SD>::inSingleCell(Point const* p) const
{
  int iC  = static_cast<PointT const*>(p)->getIncidCell();
  int viC = static_cast<PointT const*>(p)->getPosition();
  
  CT const* cell = MeshT::getCell(iC);
  
  if (viC<CellT::n_vertices) // is a vertex
  {
    int counter = 0;
    for (int fv = 0; fv<CT::dim; ++fv) // for each incident facet
    {
      int f = CT::table_vC_x_fC[ viC ][ fv ];
      if (cell->getIncidCell(f) < 0)
        ++counter;
    }
    if (counter == CT::dim)
      return true;
    else
      return false;
  }
  else
  if (viC<CellT::n_vertices + CellT::n_corners) // is an edge vertex
  {
    if (CT::dim == 1) return true;
    
    if (CT::dim == 2)
      return cell->getIncidCell(viC - CT::n_vertices) < 0;
    
    if (CT::dim == 3)
    {
      int counter = 0;
      for (int fe = 0; fe < 2; ++fe) // for each incident facet
      {
        const int edge_lid = viC - CT::n_vertices;
        int f = CT::table_bC_x_fC[ edge_lid ][ fe ];
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
  if (viC<CellT::n_vertices + CellT::n_corners + CellT::n_facets) // is an face vertex
  {
    if (CT::dim < 3) return true;
    
    int f = viC - CT::n_vertices - CT::n_corners;
    return cell->getIncidCell(f) < 0;
  }
  else // is a volume point
  {
    return true;
  }
  
}

// only dim = 3
template<class CT, int SD>
bool SMesh<CT,SD>::inSingleCell(Corner const* r) const
{
  if (CT::dim < 3)
    return false;
  
  int iC  = static_cast<CornerT const*>(r)->getIncidCell();
  int eiC = static_cast<CornerT const*>(r)->getPosition();  
  
  CT const* cell = MeshT::getCell(iC);

  int counter = 0;
  for (int fe = 0; fe < 2; ++fe) // for each incident facet
  {
    int f = CT::table_bC_x_fC[ eiC ][ fe ];
    if (cell->getIncidCell(f) < 0)
      ++counter;
  }
  if (counter == 2)
    return true;
  else
    return false;  
  
}
template<class CT, int SD>
bool SMesh<CT,SD>::inSingleCell(Facet const* fa) const
{
  int iC  = static_cast<FacetT const*>(fa)->getIncidCell();
  int fiC = static_cast<FacetT const*>(fa)->getPosition();  
  
  if (MeshT::getCell(iC)->CT::getIncidCell(fiC)<0)
    return true;
  else
    return false;
  
}

// LINEAR CASE ONLY
template<class CT, int SD>
void SMesh<CT,SD>::getCenterCoord(Cell const* cell, Real* Xc) const
{
  int i;
  for (i=0; i<SD; ++i) Xc[i]=0;
  
  for (i = 0; i < CT::n_vertices; ++i)
    for (int j = 0; j < SD; ++j)
      Xc[j] += MeshT::getNode(static_cast<CT const*>(cell)->getNodeId(i))->getCoord(j);
  
  for (i=0; i<SD; ++i)
    Xc[i] /= CT::n_vertices;
  
}
template<class CT, int SD>
void SMesh<CT,SD>::getCenterCoord(Facet const* facet, Real* Xc) const
{
  int ids[CT::n_nodes_per_facet];
  int i;
  for (i=0; i<SD; ++i) Xc[i]=0;
  
  MeshT::getFacetNodesId(facet,ids);
  
  for (int i = 0; i < CT::n_vertices_per_facet; ++i)
    for (int j = 0; j < SD; ++j)
      Xc[j] += MeshT::getNode( ids[i] )->getCoord(j);

  for (i=0; i<SD; ++i)
    Xc[i] /= CT::n_vertices_per_facet;
  
}
template<class CT, int SD>
void SMesh<CT,SD>::getCenterCoord(Corner const* corner, Real* Xc) const
{
  int ids[CT::n_nodes_per_corner];
  int i;
  for (i=0; i<SD; ++i) Xc[i]=0;
  
  MeshT::getCornerNodesId(corner,ids);
  
  for (int i = 0; i < CT::n_vertices_per_corner; ++i)
    for (int j = 0; j < SD; ++j)
      Xc[j] += MeshT::getNode( ids[i] )->getCoord(j);

  for (i=0; i<SD; ++i)
    Xc[i] /= CT::n_vertices_per_corner + (CT::dim==1?1:0);
  
}


template class SMesh<Edge2,3>;
template class SMesh<Edge3,3>;
template class SMesh<Triangle3,3>;
template class SMesh<Triangle6,3>;
template class SMesh<Quadrangle4,3>;
template class SMesh<Quadrangle8,3>;
template class SMesh<Quadrangle9,3>;
template class SMesh<Tetrahedron4,3>;
template class SMesh<Tetrahedron10,3>;
template class SMesh<Hexahedron8,3>;
template class SMesh<Hexahedron20,3>;
template class SMesh<Hexahedron27,3>;

template class SMesh<Edge2,2>;
template class SMesh<Edge3,2>;
template class SMesh<Triangle3,2>;
template class SMesh<Triangle6,2>;
template class SMesh<Quadrangle4,2>;
template class SMesh<Quadrangle8,2>;
template class SMesh<Quadrangle9,2>;

template class SMesh<Edge2,1>;
template class SMesh<Edge3,1>;





namespace _MeshStaticTablesInitializers
{
  typedef Mesh::CreatorMemFunPtr MemFunPtr;


  std::tr1::array<MemFunPtr, 3*N_CELL_TYPES>
  creators_tab()
  {
    std::tr1::array<MemFunPtr, 3*N_CELL_TYPES> tab;

    tab[0*N_CELL_TYPES + log2_i32(EDGE2        )] = &Mesh::create<Edge2,1>;
    tab[0*N_CELL_TYPES + log2_i32(EDGE3        )] = &Mesh::create<Edge3,1>;
    tab[0*N_CELL_TYPES + log2_i32(TRIANGLE3    )] = NULL;
    tab[0*N_CELL_TYPES + log2_i32(TRIANGLE6    )] = NULL;
    tab[0*N_CELL_TYPES + log2_i32(QUADRANGLE4  )] = NULL;
    tab[0*N_CELL_TYPES + log2_i32(QUADRANGLE8  )] = NULL;
    tab[0*N_CELL_TYPES + log2_i32(QUADRANGLE9  )] = NULL;
    tab[0*N_CELL_TYPES + log2_i32(TETRAHEDRON4 )] = NULL;
    tab[0*N_CELL_TYPES + log2_i32(TETRAHEDRON10)] = NULL;
    tab[0*N_CELL_TYPES + log2_i32(HEXAHEDRON8  )] = NULL;
    tab[0*N_CELL_TYPES + log2_i32(HEXAHEDRON20 )] = NULL;
    tab[0*N_CELL_TYPES + log2_i32(HEXAHEDRON27 )] = NULL;

    tab[1*N_CELL_TYPES + log2_i32(EDGE2        ) ] = &Mesh::create<Edge2,2>;
    tab[1*N_CELL_TYPES + log2_i32(EDGE3        ) ] = &Mesh::create<Edge3,2>;
    tab[1*N_CELL_TYPES + log2_i32(TRIANGLE3    ) ] = &Mesh::create<Triangle3,2>;
    tab[1*N_CELL_TYPES + log2_i32(TRIANGLE6    ) ] = &Mesh::create<Triangle6,2>;
    tab[1*N_CELL_TYPES + log2_i32(QUADRANGLE4  ) ] = &Mesh::create<Quadrangle4,2>;
    tab[1*N_CELL_TYPES + log2_i32(QUADRANGLE8  ) ] = &Mesh::create<Quadrangle8,2>;
    tab[1*N_CELL_TYPES + log2_i32(QUADRANGLE9  ) ] = &Mesh::create<Quadrangle9,2>;
    tab[1*N_CELL_TYPES + log2_i32(TETRAHEDRON4 ) ] = NULL;
    tab[1*N_CELL_TYPES + log2_i32(TETRAHEDRON10) ] = NULL;
    tab[1*N_CELL_TYPES + log2_i32(HEXAHEDRON8  ) ] = NULL;
    tab[1*N_CELL_TYPES + log2_i32(HEXAHEDRON20 ) ] = NULL;
    tab[1*N_CELL_TYPES + log2_i32(HEXAHEDRON27 ) ] = NULL;

    tab[2*N_CELL_TYPES + log2_i32(EDGE2        ) ] = &Mesh::create<Edge2,3>;
    tab[2*N_CELL_TYPES + log2_i32(EDGE3        ) ] = &Mesh::create<Edge3,3>;
    tab[2*N_CELL_TYPES + log2_i32(TRIANGLE3    ) ] = &Mesh::create<Triangle3,3>;
    tab[2*N_CELL_TYPES + log2_i32(TRIANGLE6    ) ] = &Mesh::create<Triangle6,3>;
    tab[2*N_CELL_TYPES + log2_i32(QUADRANGLE4  ) ] = &Mesh::create<Quadrangle4,3>;
    tab[2*N_CELL_TYPES + log2_i32(QUADRANGLE8  ) ] = &Mesh::create<Quadrangle8,3>;
    tab[2*N_CELL_TYPES + log2_i32(QUADRANGLE9  ) ] = &Mesh::create<Quadrangle9,3>;
    tab[2*N_CELL_TYPES + log2_i32(TETRAHEDRON4 ) ] = &Mesh::create<Tetrahedron4,3>;
    tab[2*N_CELL_TYPES + log2_i32(TETRAHEDRON10) ] = &Mesh::create<Tetrahedron10,3>;
    tab[2*N_CELL_TYPES + log2_i32(HEXAHEDRON8  ) ] = &Mesh::create<Hexahedron8,3>;
    tab[2*N_CELL_TYPES + log2_i32(HEXAHEDRON20 ) ] = &Mesh::create<Hexahedron20,3>;
    tab[2*N_CELL_TYPES + log2_i32(HEXAHEDRON27 ) ] = &Mesh::create<Hexahedron27,3>;

    return tab;
  }

};



Mesh* Mesh::create(ECellType type, int spacedim)
{

  if (static_cast<unsigned>(spacedim-1)>2)
    spacedim = dimForCtype(type);

  FEPIC_CHECK(spacedim>=dimForCtype(type), "space dimension must be greater than or equal to cell's dimension ", std::invalid_argument);

  static const
  std::tr1::array<CreatorMemFunPtr, 3*N_CELL_TYPES> creators = _MeshStaticTablesInitializers::creators_tab();

  unsigned idx = log2_i32(type);
  if (idx >= N_CELL_TYPES)
  {
    FEPIC_CHECK(false, "invalid or not supported mesh type", std::invalid_argument);
    return NULL;
  }

  return (*(creators[(spacedim-1)*N_CELL_TYPES + idx]))();
};












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



