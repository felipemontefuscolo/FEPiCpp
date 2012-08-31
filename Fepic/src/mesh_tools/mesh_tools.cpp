#include "mesh_tools.hpp"
#include "Fepic/src/mesh/mesh.hpp"
#include <vector>
//#include <type_traits>

/*
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
*/

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


/** flips Tri3 only TODO
 *  @param acell the cell that will be removed.
 *  @param afid local facet's id that will be flipped.
 *  @param mesh mesh context.
 *  @return true if an error occurred, false otherwise.
 */ 
bool MeshTools::flipTri(Cell * acell, int afid, Mesh *mesh)
{
  /*   
      A-reference                          | B-reference    
                                           |
           j                  j            |        k                             
           ^                  ^            |        ^                  ^          
          /|\                / \           |       /|\                / \         
         / | \        afid  /   \          |      / | \              /   \        
        /  |  \            /  A  \         |     /  |  \            /     \       
       /   |   \          /       \        |    /   |   \          /       \      
      <    | A  > i    k <---------> i     | i <  B |    >      i <---------> k   
       \   |   /          \       /        |    \   |   /          \       /      
        \  |  /            \     /         |     \  |  /            \  B  /       
         \ | /              \   /          |      \ | /              \   /  bfid 
          \|/                \ /           |       \|/                \ /         
           v                  v            |        v                  v          
           k                               |        j                  j          
   
  */

  int const bcid = acell->getIncidCell(afid);
  int const bfid = acell->getIncidCellPos(afid);
  int const acid = mesh->getCellId(acell);
  
  FEPIC_CHECK(bcid >= 0 && bfid>=0, "invalid mesh or argument", std::runtime_error);
  
  Cell *bcell = mesh->getCell(bcid);


  /* changing the incidences of the facets */
  Facet *f;
  f = mesh->getFacet(acell->getFacetId((afid+1)%3)); // inferior-direita
  f->setIncidCell(bcid);
  f->setPosition(bfid);
  
  f = mesh->getFacet(bcell->getFacetId((bfid+1)%3)); // superior-esquerda
  f->setIncidCell(acid);
  f->setPosition(afid);
  
  f = mesh->getFacet(acell->getFacetId(afid)); // meio
  f->setIncidCell(acid);
  f->setPosition((afid+1)%3);
  // as outras faces são preservadas
  
  /* changing the incidences of the nodes */
  Point *p;
  p = mesh->getNode(acell->getNodeId(afid)); // j de a
  p->setIncidCell(acid);
  p->setPosition(afid);
  
  p = mesh->getNode(bcell->getNodeId(bfid)); // j de b
  p->setIncidCell(bcid);
  p->setPosition(bfid);
  
  /* changing the incidences of the cells */
  Cell *c;
  // superior-esquerda
  acell->setIncidCell(afid, bcell->getIncidCell((bfid+1)%3));
  acell->setIncidCellPos(afid, bcell->getIncidCellPos((bfid+1)%3));
  c = mesh->getCell(acell->getIncidCell(afid));
  c->setIncidCell(acell->getIncidCellPos(afid), acid);
  c->setIncidCellPos(acell->getIncidCellPos(afid),afid);
  // inferior-direita
  bcell->setIncidCell(bfid, acell->getIncidCell((afid+1)%3));
  bcell->setIncidCellPos(bfid, acell->getIncidCellPos((afid+1)%3));
  c = mesh->getCell(bcell->getIncidCell(bfid));
  c->setIncidCell(bcell->getIncidCellPos(bfid), bcid);
  c->setIncidCellPos(bcell->getIncidCellPos(bfid),bfid);
  // meio
  acell->setIncidCell((afid+1)%3, bcid);
  acell->setIncidCellPos((afid+1)%3, (bfid+1)%3);
  bcell->setIncidCell((bfid+1)%3, acid);
  bcell->setIncidCellPos((bfid+1)%3, (afid+1)%3);  
  
  int f_ie, f_id, f_se, f_sd, f_m;   // i=inferior; s=superior; e=esquerda; d=direita
  
  f_ie = bcell->getFacetId((bfid+2)%3);
  f_id = acell->getFacetId((afid+1)%3);
  f_se = bcell->getFacetId((bfid+1)%3);
  f_sd = acell->getFacetId((afid+2)%3);
  f_m  = acell->getFacetId(afid);
  
  acell->setFacetId((afid+0)%3, f_se);
  acell->setFacetId((afid+1)%3, f_m );
  acell->setFacetId((afid+2)%3, f_sd);
  
  bcell->setFacetId((bfid+0)%3, f_id);
  bcell->setFacetId((bfid+1)%3, f_m );
  bcell->setFacetId((bfid+2)%3, f_ie);
  
  acell->setNode((afid+1)%3, bcell->getNodeId((bfid+2)%3));
  bcell->setNode((bfid+1)%3, acell->getNodeId((afid+2)%3));

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


