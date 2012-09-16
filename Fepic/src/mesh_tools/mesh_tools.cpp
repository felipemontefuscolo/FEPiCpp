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

#include "mesh_tools.hpp"
#include "Fepic/src/mesh/mesh.hpp"
#include <vector>
//#include <type_traits>


// WARNING:
/*
 *  Antes de implementar uma operação sobre a malha, certifique-se
 *  operação mantém a consistência da mesma, verificando os seguintes
 *  itens:
 * 
 *  - tag of all elements.
 * 
 *  - nodes : _icell                   
 *            _icell_pos               
 *            _status                  
 *            _incidences              
 *                                     
 *  - facets : _icell                  
 *             _icell_pos              
 *             _bound_comp_id                 
 * 
 *  - corners : _icell
 *              _icell_pos
 * 
 *  - cells : _icells_pos 
 *            _icells_anchors
 *            _facets
 *            _icells
 *            _corners
 *            _nodes          
 *            _conn_comp_id
 * 
 *  - SMesh : _connected_compL // (connected component vs initial cell id) list
 *            _boundary_compL  // (boundary component vs initialfacet id) list
 * 
 */ 


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

/** Safely removes a cell
 *  @param cell the cell that will be removed.
 *  @param mesh mesh context.
 *  @return nothing.
 */ 
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

/** Read a mesh from input arrays.
 * 
 *  @param n_cells number of cells of the mesh.
 *  
 *  @param n_nodes number of nodes of the mesh.
 * 
 *  @param nodes a vector containing the nodes of the mesh where each
 *          cell c has nnpc nodes, and where nodes[nnpe*c + i] is the i-th
 *          node of the cell c.
 *  @param xyz coordinates of the nodes, where xyz[spacedim*n + i] is
 *          the ith coordinate of the node n.
 *  
 *  @param[out] mesh a pointer to the mesh.
 */ 
void MeshTools::readMesh(int n_nodes, int n_cells, int const* nodes, Real const* xyz, Mesh *mesh)
{
  FEPIC_ASSERT(mesh!=NULL, "invalid argument", std::runtime_error);
  
  int const sdim = mesh->spaceDim();
  int const nnpe = mesh->numNodesPerCell();
  
  mesh->resizePointL(n_nodes);
  mesh->resizeCellL(n_cells);

  FEP_PRAGMA_OMP(paralle for)
  for (int n = 0; n < n_nodes; ++n)
      mesh->getNode(n)->setCoord(&xyz[sdim*n], sdim);
  
  FEP_PRAGMA_OMP(paralle for)
  for (int c = 0; c < n_cells; ++c)
    for (int i = 0; i < nnpe; ++i)
      mesh->getCell(c)->setNodeId(i, nodes[nnpe*c + i]);
  
  if (mesh->qBuildAdjacency())
    mesh->buildAdjacency();
  
}


/** edge flipping for Tri3 and Tri6 cells.
 *  @param acell the cell that will be removed.
 *  @param afid local facet's id that will be flipped.
 *  @param mesh mesh context.
 *  @param move_edge_nds true to move high order nodes, false otherwise.
 *  @return true if an error occurred, false otherwise.
 *  @note none of the ids (of the edge or cells) are changed, only the incidences.
 *  @note for 2D and 3D space.
 *  @note all tags are kept.
 *  @pre the cells acell and acell->getIncidCell(afid) form a convex quadrilateral.
 */ 
bool MeshToolsTri::flipEdge(Cell * acell, int afid, Mesh *mesh, bool move_edge_nds)
{
  /*   
      A-reference                          | B-reference    
                                           |
           j                  j            |        k                             
           ^                  ^            |        ^                  ^          
          /|\         new    / \           |       /|\                / \         
         / | \        afid  /   \          |      / | \              /   \        
        /  |  \            /  A  \         |     /  |  \            /     \       
       /   |   \          /       \        |    /   |   \          /       \      
      <    | A  > i    k <---------> i     | i <  B |    >      i <---------> k   
       \   |   /          \       /        |    \   |   /          \       /      
        \  |  /            \     /         |     \  |  /            \  B  / new   
         \ | /              \   /          |      \ | /              \   /  bfid 
          \|/                \ /           |       \|/                \ /         
           v                  v            |        v                  v          
           k                               |        j                  j          
   
  */

  int const bcid = acell->getIncidCell(afid);
  int const bfid = acell->getIncidCellPos(afid);
  int const acid = mesh->getCellId(acell);
  
  //int const afid_plus0 = (afid+0)%3;
  int const afid_plus1 = (afid+1)%3;
  int const afid_plus2 = (afid+2)%3;
  //int const bfid_plus0 = (bfid+0)%3;
  int const bfid_plus1 = (bfid+1)%3;
  int const bfid_plus2 = (bfid+2)%3;
  
  FEPIC_CHECK(bcid >= 0 && bfid>=0, "invalid mesh or argument", std::runtime_error);
  
  Cell *bcell = mesh->getCell(bcid);

  /* changing the incidences of the facets */
  Facet *f;
  f = mesh->getFacet(acell->getFacetId(afid_plus1)); // inferior-direita
  f->setIncidence(bcid,bfid);
  
  f = mesh->getFacet(bcell->getFacetId(bfid_plus1)); // superior-esquerda
  f->setIncidence(acid, afid);
  
  f = mesh->getFacet(acell->getFacetId(afid)); // meio
  f->setIncidence(acid,afid_plus1);
  // as outras faces são preservadas
  
  /* changing the incidences of the nodes only when is necessary */
  Point *p;
  p = mesh->getNode(acell->getNodeId(afid)); // j de a
  if (p->getIncidCell() == bcid)
    p->setIncidence(acid, afid);

  p = mesh->getNode(bcell->getNodeId(bfid)); // j de b
  if (p->getIncidCell() == acid)
    p->setIncidence(bcid, bfid);
  
  /* changing the incidences of the cells */
  Cell *c;
  // superior-esquerda
  acell->setIncidence(afid, bcell->getIncidCell(bfid_plus1), bcell->getIncidCellPos(bfid_plus1));
  if (acell->getIncidCell(afid) >= 0)
  {
    c = mesh->getCell(acell->getIncidCell(afid));
    c->setIncidence(acell->getIncidCellPos(afid), acid, afid);
  }
  // inferior-direita
  bcell->setIncidence(bfid, acell->getIncidCell(afid_plus1), acell->getIncidCellPos(afid_plus1));
  if (bcell->getIncidCell(bfid) >= 0)
  {
    c = mesh->getCell(bcell->getIncidCell(bfid));
    c->setIncidence(bcell->getIncidCellPos(bfid), bcid, bfid);
  }
  // meio (tem que ser depois mesmo)
  acell->setIncidence(afid_plus1, bcid, bfid_plus1);
  bcell->setIncidence(bfid_plus1, acid, afid_plus1);
  
  //int f_ie, f_id, f_se, f_sd, f_m;   // i=inferior; s=superior; e=esquerda; d=direita
  int f_id, f_se, f_m;   // i=inferior; s=superior; e=esquerda; d=direita
  
  
  //f_ie = bcell->getFacetId(bfid_plus2);
  f_id = acell->getFacetId(afid_plus1);
  f_se = bcell->getFacetId(bfid_plus1);
  //f_sd = acell->getFacetId(afid_plus2);
  f_m  = acell->getFacetId(afid);
  
  acell->setFacetId(afid      , f_se);
  acell->setFacetId(afid_plus1, f_m );
  //acell->setFacetId(afid_plus2, f_sd // isnt necessary
  
  bcell->setFacetId(bfid      , f_id);
  bcell->setFacetId(bfid_plus1, f_m );
  //bcell->setFacetId(bfid_plus2, f_ie); // isnt necessary
  
  acell->setNodeId(afid_plus1, bcell->getNodeId(bfid_plus2));
  bcell->setNodeId(bfid_plus1, acell->getNodeId(afid_plus2));


  // for high-order-cells
  if (mesh->numNodesPerCell() > mesh->numVerticesPerCell())
  {
    // high-order nodes id ... i=inferior; s=superior; e=esquerda; d=direita
    //f_ie = bcell->getNodeId(bfid_plus2+3);
    f_id = acell->getNodeId(afid_plus1+3);
    f_se = bcell->getNodeId(bfid_plus1+3);
    //f_sd = acell->getNodeId(afid_plus2+3);
    f_m  = acell->getNodeId(afid+3);
    
    acell->setNodeId(afid      +3, f_se);
    acell->setNodeId(afid_plus1+3, f_m );
    //acell->setNodeId(afid_plus2+3, f_sd);
    
    bcell->setNodeId(bfid      +3, f_id);
    bcell->setNodeId(bfid_plus1+3, f_m );
    //bcell->setNodeId(bfid_plus2+3, f_ie);
    
    mesh->getNode(f_id)->setIncidence(bcid, bfid+3);       // inferior-direita
    mesh->getNode(f_se)->setIncidence(acid, afid+3);       // superior-esquerda
    mesh->getNode(f_m )->setIncidence(acid, afid_plus1+3); // meio
    
    if (move_edge_nds)
    {
      int const sdim = mesh->spaceDim();
      Real pleft[3], pright[3];
      
      mesh->getNode(acell->getNodeId(afid_plus2))->getCoord(pright,sdim); // right node
      mesh->getNode(bcell->getNodeId(bfid_plus2))->getCoord(pleft,sdim);  // left node
      
      // compute the center of the new edge
      pleft[0] = 0.5*(pleft[0]+pright[0]);
      pleft[1] = 0.5*(pleft[1]+pright[1]);
      if (sdim==3)
        pleft[2] = 0.5*(pleft[2]+pright[2]);
      
      p = mesh->getNode(acell->getNodeId(afid_plus1+3)); // mid node
      p->setCoord(pleft, sdim);
    }
    
  }

  return false;
}

/** Checks if an edge is a Delaunay edge.
 *  @param cell edge's incident cell
 *  @param fid edge's position in the cell.
 *  @param mesh mesh context.
 *  @return true if the edge is Delaunay, false otherwise.
 *  @note for two-dimensinal space only.
 */ 
bool MeshToolsTri::inCircle2d(Cell const* cell, int const fid, Mesh const* mesh)
{
  int const ocid = cell->getIncidCell(fid);
  if (ocid<0)
    return true;

  int const ofid = cell->getIncidCellPos(fid);
  Cell const *ocell = mesh->getCell(ocid);

  Point const *a = mesh->getNode(cell->getNodeId(0));
  Point const *b = mesh->getNode(cell->getNodeId(1));
  Point const *c = mesh->getNode(cell->getNodeId(2));
  Point const *d = mesh->getNode(ocell->getNodeId((ofid+2)%3));
 
  Real const ax = a->getCoord(0), ay = a->getCoord(1), aq = ax*ax + ay*ay;
  Real const bx = b->getCoord(0), by = b->getCoord(1), bq = bx*bx + by*by;
  Real const cx = c->getCoord(0), cy = c->getCoord(1), cq = cx*cx + cy*cy;
  Real const dx = d->getCoord(0), dy = d->getCoord(1), dq = dx*dx + dy*dy;
  /*
     criteria:
     true if
         | 1 ax ay ax^2 + ay^2 |
     det | 1 bx by bx^2 + by^2 | < 0
         | 1 cx cy cx^2 + cy^2 |
         | 1 dx dy dx^2 + dy^2 |
      
         | A  B  C |
     det | D  E  F | = A*(E*I - F*H) - B*(D*I - F*G) + C*(D*H - E*G)
         | G  H  I |
  */
  Real const A = bx-ax, B = by-ay, C = bq-aq;
  Real const D = cx-ax, E = cy-ay, F = cq-aq;
  Real const G = dx-ax, H = dy-ay, I = dq-aq;
  
  Real const det = A*(E*I - F*H) - B*(D*I - F*G) + C*(D*H - E*G);

  if (det<0)
    return false;
  else
    return true;
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
  
}



