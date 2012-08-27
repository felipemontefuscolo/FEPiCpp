#include "mesh_tools.hpp"
#include "Fepic/src/mesh/mesh.hpp"
#include <vector>
//#include <type_traits>


//int MeshTools::checkConsistency(Mesh *mesh)
//{
//  cell_iterator cell = mesh->cellBegin();
//  cell_iterator cell_end = mesh->cellEnd();
//  
//  int const nfpc = mesh->numFacetsPerCell();
//  //  int const nnpc = mesh->numNodesPerCell();
//  int const nnpf = mesh->numNodesPerFacet();
//  int const dim  = mesh->spaceDim();
//  
//  int f_nds[nnpf];
//  
//  //  Point *p;
//  Facet *f;
//  Cell  *c;
//  
//  for (; cell != cell_end; ++cell)
//  {
//    int myid = mesh->getCellId(&*cell);
//    for (int i = 0; i < nfpc; ++i)
//    {
//      if (cell->getIncidCell(i) >= 0)
//      {
//        // checks whether the neighbor cell contains this cell.
//        c = mesh->getCell(cell->getIncidCell(i));
//        int pos = cell->getIncidCellPos(i);
//        if (myid != c->getIncidCell(pos))
//          return -1;
//        
//        // checks facets
//        if (!(static_cast<unsigned>(cell->getFacetId(i)) < mesh->numFacetsTotal()))
//          return -2;
//        f = mesh->getFacet(cell->getFacetId(i));
//        int icf = f->getIncidCell();
//        if (!(icf==myid || icf==cell->getIncidCell(i)))
//          return -3;
//        if (icf==myid)
//          EXPECT_TRUE(f->getPosition() == i);
//        else
//        {
//          EXPECT_TRUE(f->getPosition() == pos);
//        }
//      }
//      else // bordo
//      {
//        // verifica a face
//        f = mesh->getFacet(cell->getFacetId(i));
//        int icf = f->getIncidCell();
//        // só pode ser o myid, pq do outro lado não tem ngm
//        EXPECT_TRUE(icf==myid);
//        
//        // verifica se os nós da face estão no contorno
//        mesh->getFacetNodesId(f, f_nds);
//        for (int j = 0; j < nnpf; ++j)
//        {
//          EXPECT_TRUE(mesh->inBoundary(mesh->getNode(f_nds[j])));
//        }
//        
//      }
//    }
//    
//  }
//  
//  for (point_iterator point = mesh->pointBegin(); point != mesh->pointEnd(); ++point)
//  {
//    int myid = mesh->getPointId(&*point);
//    int ic = point->getIncidCell();
//    int pos = point->getPosition();
//    Cell *c = mesh->getCell(ic);
//    
//    EXPECT_TRUE(c->getNodeId(pos) == myid);
//    
//  }
//  
//  
//}
//

// da pra melhorar
void MeshTools::removeCell(Cell * cell, Mesh *mesh)
{
  const int cell_dim = mesh->cellDim();
  const int nds_per_cell = mesh->numNodesPerCell();

  const int n_facets_per_cell = mesh->numFacetsPerCell();
  const int n_corners_per_cell = mesh->numCornersPerCell();

  typedef unsigned char MyBool;



  int id;

  // killing nodes
  for (int i = 0; i < nds_per_cell; ++i)
  {
    id = cell->getNodeId(i);
    if (mesh->inSingleCell(mesh->getNode(id)))
      mesh->disablePoint(id);
  }

  // killing facets
  if (cell_dim > 1)
    for (int i = 0; i < n_facets_per_cell; ++i)
    {
      id = cell->getFacetId(i);
      if (mesh->inSingleCell(mesh->getFacet(id)))
        mesh->disableFacet(id);
    }

  // killing corners
  if (cell_dim > 2)
    for (int i = 0; i < n_corners_per_cell; ++i)
    {
      id = cell->getCornerId(i);
      if (mesh->inSingleCell(mesh->getCorner(id)))
        mesh->disableCorner(id);
    }



  Point  * pt;
  Facet  * ft;
  Corner * cr;

  Cell * neib;
  for (int i = 0; i < n_facets_per_cell; ++i)
  {
    const int neib_id = cell->getIncidCell(i);
    if (neib_id < 0)
      continue;
    neib = mesh->getCell( neib_id );

    // nodes
    for (int j = 0; j < nds_per_cell; ++j)
    {
      pt = mesh->getNode(neib->getNodeId(j));
      pt->setIncidCell(neib_id);
      pt->setPosition(j);
    }

    //facets
    if (cell_dim > 1)
      for (int j = 0; j < n_facets_per_cell; ++j)
      {
        ft = mesh->getFacet(neib->getFacetId(j));
        ft->setIncidCell(neib_id);
        ft->setPosition(j);
      }

    //corners
    if (cell_dim > 2)
      for (int j = 0; j < n_corners_per_cell; ++j)
      {
        cr = mesh->getCorner(neib->getCornerId(j));
        cr->setIncidCell(neib_id);
        cr->setPosition(j);
      }
      
    // neib
    int const ith = cell->getIncidCellPos(i);
    neib->setIncidCell(ith, -1);
    
  }


  // killing cells
  id = mesh->getCellId(cell);
  mesh->disableCell(id);
}


std::pair<bool, Cell *> MeshToolsTri::searchConvexPoint(Real const* x, Cell const* c0, Mesh const* mesh)
{
  Cell const* cell = c0;
  
  double const* x0;
  double const* x1;
  double const* x2;
  double L0, L1, L2, a0, a1, a2, l;
  int iC;

  while(true)
  {
    x0 = mesh->getNode(cell->getNodeId(0))->getCoord();
    x1 = mesh->getNode(cell->getNodeId(1))->getCoord();
    x2 = mesh->getNode(cell->getNodeId(2))->getCoord();
    
    a0 = x1[0]*x2[1] - x2[0]*x1[1];
    a1 = x2[0]*x0[1] - x0[0]*x2[1];
    a2 = x0[0]*x1[1] - x1[0]*x0[1];
    
    L0 = a0 + (x1[1] - x2[1])*x[0] + (x2[0] - x1[0])*x[1];
    L1 = a1 + (x2[1] - x0[1])*x[0] + (x0[0] - x2[0])*x[1];
    L2 = a0+a1+a2 - L0 - L1;
    
    if (L0<0)
    {
      if (L0<L1)
      {
        if (L0<L2)
          l = 1;
        else
          l = 0;
      }
      else
        l = 2;
    }
    else if (L1<0)
    {
      if (L1<L2)
        l = 2;
      else
        l = 0;
    }
    else if (L2<0)
    {
      l=0;
    }
    else
      return std::pair<bool, Cell *>(true, const_cast<Cell *>(cell));
    
    iC = cell->getIncidCell(l);
    if (iC >= 0)
      cell = mesh->getCell(cell->getIncidCell(l));
    else
      return std::pair<bool, Cell *>(false, const_cast<Cell *>(cell));
    
    
  }
  
};

