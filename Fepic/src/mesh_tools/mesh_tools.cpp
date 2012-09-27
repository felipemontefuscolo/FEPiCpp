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

#include <tr1/memory>
#include "mesh_tools.hpp"
#include "Fepic/src/mesh/mesh.hpp"
#include <vector>
#include <algorithm>
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
//        c = mesh->getCellPtr(cell->getIncidCell(i));
//        int pos = cell->getIncidCellPos(i);
//        if (myid != c->getIncidCell(pos))
//          return -1;
//
//        // checks facets
//        if (!(static_cast<unsigned>(cell->getFacetId(i)) < mesh->numFacetsTotal()))
//          return -2;
//        f = mesh->getFacetPtr(cell->getFacetId(i));
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
//        f = mesh->getFacetPtr(cell->getFacetId(i));
//        int icf = f->getIncidCell();
//        // só pode ser o myid, pq do outro lado não tem ngm
//        EXPECT_TRUE(icf==myid);
//
//        // verifica se os nós da face estão no contorno
//        mesh->getFacetNodesId(f, f_nds);
//        for (int j = 0; j < nnpf; ++j)
//        {
//          EXPECT_TRUE(mesh->inBoundary(mesh->getNodePtr(f_nds[j])));
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
//    Cell *c = mesh->getCellPtr(ic);
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
  const int cell_id_to_remove = mesh->getCellId(cell);
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
    if (mesh->inSingleCell(mesh->getNodePtr(id)))
      mesh->disablePoint(id);
  }

  // killing facets
  if (cell_dim > 1)
    for (int i = 0; i < n_facets_per_cell; ++i)
    {
      id = cell->getFacetId(i);
      if (mesh->inSingleCell(mesh->getFacetPtr(id)))
        mesh->disableFacet(id);
    }

  // killing corners
  if (cell_dim > 2)
    for (int i = 0; i < n_corners_per_cell; ++i)
    {
      id = cell->getCornerId(i);
      if (mesh->inSingleCell(mesh->getCornerPtr(id)))
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
    neib = mesh->getCellPtr( neib_id );

    // nodes
    for (int j = 0; j < nds_per_cell; ++j)
    {
      pt = mesh->getNodePtr(neib->getNodeId(j));
      pt->setIncidCell(neib_id);
      pt->setPosition(j);
    }

    //facets
    if (cell_dim > 1)
      for (int j = 0; j < n_facets_per_cell; ++j)
      {
        ft = mesh->getFacetPtr(neib->getFacetId(j));
        ft->setIncidCell(neib_id);
        ft->setPosition(j);
      }

    //corners
    if (cell_dim > 2)
      for (int j = 0; j < n_corners_per_cell; ++j)
      {
        cr = mesh->getCornerPtr(neib->getCornerId(j));
        cr->setIncidCell(neib_id);
        cr->setPosition(j);
      }

    // neib
    int const ith = cell->getIncidCellPos(i);
    neib->setIncidCell(ith, -1);

  }


  // killing cells
  mesh->disableCell(cell_id_to_remove);
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
  std::tr1::shared_ptr<Point> p_ptr(mesh->createPoint());
  std::tr1::shared_ptr<Cell>  c_ptr(mesh->createCell());

  //mesh->resizePointL(n_nodes);
  //mesh->resizeCellL(n_cells);

  
  for (int n = 0; n < n_nodes; ++n)
  {
    p_ptr->setCoord(&xyz[sdim*n], sdim);
    mesh->pushPoint(p_ptr.get());
  }

  for (int c = 0; c < n_cells; ++c)
  {
    for (int i = 0; i < nnpe; ++i)
      c_ptr->setNodeId(i, nodes[nnpe*c + i]);
    mesh->pushCell(c_ptr.get());
  }

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

  Cell *bcell = mesh->getCellPtr(bcid);

  /* changing the incidences of the facets */
  Facet *f;
  f = mesh->getFacetPtr(acell->getFacetId(afid_plus1)); // inferior-direita
  f->setIncidence(bcid,bfid);

  f = mesh->getFacetPtr(bcell->getFacetId(bfid_plus1)); // superior-esquerda
  f->setIncidence(acid, afid);

  f = mesh->getFacetPtr(acell->getFacetId(afid)); // meio
  f->setIncidence(acid,afid_plus1);
  // as outras faces são preservadas

  /* changing the incidences of the nodes only when is necessary */
  Point *p;
  p = mesh->getNodePtr(acell->getNodeId(afid)); // j de a
  if (p->getIncidCell() == bcid)
    p->setIncidence(acid, afid);

  p = mesh->getNodePtr(bcell->getNodeId(bfid)); // j de b
  if (p->getIncidCell() == acid)
    p->setIncidence(bcid, bfid);

  /* changing the incidences of the cells */
  Cell *c;
  // superior-esquerda
  acell->setIncidence(afid, bcell->getIncidCell(bfid_plus1), bcell->getIncidCellPos(bfid_plus1));
  if (acell->getIncidCell(afid) >= 0)
  {
    c = mesh->getCellPtr(acell->getIncidCell(afid));
    c->setIncidence(acell->getIncidCellPos(afid), acid, afid);
  }
  // inferior-direita
  bcell->setIncidence(bfid, acell->getIncidCell(afid_plus1), acell->getIncidCellPos(afid_plus1));
  if (bcell->getIncidCell(bfid) >= 0)
  {
    c = mesh->getCellPtr(bcell->getIncidCell(bfid));
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

    mesh->getNodePtr(f_id)->setIncidence(bcid, bfid+3);       // inferior-direita
    mesh->getNodePtr(f_se)->setIncidence(acid, afid+3);       // superior-esquerda
    mesh->getNodePtr(f_m )->setIncidence(acid, afid_plus1+3); // meio

    if (move_edge_nds)
    {
      int const sdim = mesh->spaceDim();
      Real pleft[3], pright[3];

      mesh->getNodePtr(acell->getNodeId(afid_plus2))->getCoord(pright,sdim); // right node
      mesh->getNodePtr(bcell->getNodeId(bfid_plus2))->getCoord(pleft,sdim);  // left node

      // compute the center of the new edge
      pleft[0] = 0.5*(pleft[0]+pright[0]);
      pleft[1] = 0.5*(pleft[1]+pright[1]);
      if (sdim==3)
        pleft[2] = 0.5*(pleft[2]+pright[2]);

      p = mesh->getNodePtr(acell->getNodeId(afid_plus1+3)); // mid node
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
  Cell const *ocell = mesh->getCellPtr(ocid);

  Point const *a = mesh->getNodePtr(cell->getNodeId(0));
  Point const *b = mesh->getNodePtr(cell->getNodeId(1));
  Point const *c = mesh->getNodePtr(cell->getNodeId(2));
  Point const *d = mesh->getNodePtr(ocell->getNodeId((ofid+2)%3));

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

/**
 *  Insert a vertex r on an edge defined by two vertices p and q, in direction p -> q,
 *  splitting the edge in two parts.
 *  The vertex r is given by r = t*(p-q) + q, where t is a given parameter in ]0,1[.
 *  If the cell is of type Triangle6, high order nodes are placed in the middle of its edge.
 *  @param cellA a pointer to cell where the edge lives.
 *  @param fidA loca id of the edge on cellA.
 *  @param t a real value between 0 and 1, such that r = t*(p-q) + q.
 *  @param mesh the mesh context.
 *  @return a pointer to the inserted vertex.
 *  @note new vertex tag is same of the splitted edge. The two parts of the splitted edge, has the same
 *        tag of the old edge.
 *
 */
Point* MeshToolsTri::insertVertexOnEdge(Cell *cellA, int fidA, Real t, Mesh *mesh)
{
/* ---------------------------------------------------------------------------------------------------------------
 *                      vtx_t                                                              vtx_t
 *                  ,/|\,                                                              ,/|\
 *                ,/  |  \,                                                          ,/  | `\
 *   cellB      ,/    |    \,    cellA                                  cellB      ,/    |   `\    cellA    <-- edge_t
 *            ,/      |      \,                                                  ,/      |vtx_m`\             = edge_old
 *          ,/        |        \,                    to                        ,/      _-+-_     `\
 *        ,/          | edge_    \,                     \,                   ,/   _,--'  |  '--,_  `\,
 * vtx_l /            | old        \ vtx_r                \,                /_,--'edge_l | edge_r`--,_\
 *       \,           |           ,/            ------------\               \,           |          ,/
 *         \,         |         ,/              -----------,/                 \,         |        ,/
 *           \,       |       ,/                         ,/                     \,       |       ,/         <-- edge_b
 *             \,     |     ,/                          /                         \,     |     ,/
 *               \,   |   ,/                                            cellC       \,   |   ,/    cellD
 *                 \, | ,/                                                            \, | ,/
 *                   \|/                                                                \|/
 *                     vtx_b                                                               vtx_b
 *
 * cellA nodes order:   j                                                 cellC nodes order:     2
 *                      | \ i                                                                0 / |
 *                      | /                                                                    \ |
 *                      k                                                                        1
 *
 * cellB nodes order:    k                                                cellD nodes order:  1
 *                   i / |                                                                    | \ 0
 *                     \ |                                                                    | /
 *                       j                                                                    2
 *
 * nod_tl, nod_tr, nod_bl, nod_br, nod_t*, nod_b, nod_l, nod_r = high order nodes
 *
 * nod_top = node of the old edge
 *
 *
 *///-------------------------------------------------------------------------------------------------------------

  int const  cidA = mesh->getCellId(cellA);
  
  FEPIC_CHECK((t<1 && t>0) && (cidA>=0 && cidA<mesh->numCellsTotal()), "t must be in ]0,1[", std::runtime_error);

  int const  fidA_plus1 = (fidA+1)%3;
  int const  fidA_plus2 = (fidA+2)%3;
  
  int const  cidB = cellA->getIncidCell(fidA);
  Real       xyz_tmp[3];
  int const  sdim = mesh->spaceDim();
  int        nds[6], fts[3], ics[3], ics_pos[3];
  int        pos, tag, stat, bnd;
  const bool edge_in_boundary = cidB < 0;
  const bool high_order = mesh->numNodesPerCell() > mesh->numVerticesPerCell();

  int        fidB = -1;
  //int        fidB_plus1 = -1;
  int        fidB_plus2 = -1;

  int        cidC = -1;
  int        cidD = -1;

  // creating cells to store informations ... the adjacencies
  // of the old cells are kept.
  Cell *cellB = NULL;
  Cell *cellC = NULL;
  Cell *cellD = mesh->pushCell(&cidD); //Cell *cellD = mesh->createCell(); cidD = mesh->pushCell(cellD);

  if (!edge_in_boundary)
  {
    fidB = cellA->getIncidCellPos(fidA);
    //fidB_plus1 = (fidB+1)%3;
    fidB_plus2 = (fidB+2)%3;

    cellB = mesh->getCellPtr(cidB);
    cellC = mesh->pushCell(&cidC); //cellC = mesh->createCell(); cidC = mesh->pushCell(cellC);
  }

  /* ---- First SetUp Lower dimension elements first ------  */
  // nodes First
  int const vtx_t_id = cellA->getNodeId(fidA      );
  int const vtx_b_id = cellA->getNodeId(fidA_plus1);
  int const vtx_r_id = cellA->getNodeId(fidA_plus2);
  int        vtx_l_id = -1;
  int        vtx_m_id;
  int        nod_t_id;
  int        nod_b_id;
  int        nod_r_id;
  int        nod_l_id;
  //int        nod_tr_id;
  //int        nod_tl_id;
  int        nod_br_id;
  int        nod_bl_id;

  Point *vtx_t = mesh->getNodePtr(vtx_t_id);
  Point *vtx_b = mesh->getNodePtr(vtx_b_id);
  Point *vtx_l = NULL;
  Point *vtx_r = mesh->getNodePtr(vtx_r_id);
  Point *vtx_m = mesh->pushPoint(&vtx_m_id);
  // high order nodes
  Point *nod_t = NULL;
  Point *nod_b = NULL;
  Point *nod_l = NULL;
  Point *nod_r = NULL;
  //Point *nod_tr = NULL;
  //Point *nod_tl = NULL;
  Point *nod_br = NULL;
  Point *nod_bl = NULL;
  Real const *vb_xyz = vtx_b->getCoord();
  Real const *vt_xyz = vtx_t->getCoord();
  Real const *vl_xyz = NULL;
  Real const *vr_xyz = vtx_r->getCoord();
  Real const *vm_xyz = NULL;

  if (!edge_in_boundary)
  {
    vtx_l_id = cellB->getNodeId(fidB_plus2);
    vtx_l = mesh->getNodePtr(vtx_l_id);
    vl_xyz = vtx_l->getCoord();
  }
  if (high_order)
  {
    nod_t_id  = cellA->getNodeId(fidA + 3); // node of the old edge
    //nod_tr_id = cellA->getNodeId(fidA_plus2 + 3);
    nod_br_id = cellA->getNodeId(fidA_plus1 + 3);

    nod_t  = mesh->getNodePtr(nod_t_id);
    //nod_tr = mesh->getNodePtr(nod_tr_id);
    nod_br = mesh->getNodePtr(nod_br_id);
    nod_b  = mesh->pushPoint(&nod_b_id);
    nod_r  = mesh->pushPoint(&nod_r_id);

    if (!edge_in_boundary)
    {
      //nod_tl_id = cellB->getNodeId(fidB_plus1 + 3);
      nod_bl_id = cellB->getNodeId(fidB_plus2 + 3);
      nod_l  = mesh->pushPoint(&nod_l_id);
      //nod_tl = mesh->getNodePtr(nod_tl_id);
      nod_bl = mesh->getNodePtr(nod_bl_id);
    }
  }

  /* -------- Edges ------  */

  int const edge_t_id  = cellA->getFacetId(fidA); // = old edge
  int       edge_b_id;
  int       edge_l_id  = -1;
  int       edge_r_id;
  //int const edge_tr_id = cellA->getFacetId(fidA_plus2);
  //int       edge_tl_id = -1;
  int const edge_br_id = cellA->getFacetId(fidA_plus1);
  int       edge_bl_id = -1;

  Facet *edge_t = mesh->getFacetPtr(edge_t_id);
  Facet *edge_b = mesh->pushFacet(&edge_b_id);
  Facet *edge_l = NULL;
  Facet *edge_r = mesh->pushFacet(&edge_r_id);;
  //Facet *edge_tr = mesh->getFacetPtr(edge_tr_id);
  //Facet *edge_tl = NULL;
  Facet *edge_br = mesh->getFacetPtr(edge_br_id);
  Facet *edge_bl = NULL;

  if (!edge_in_boundary)
  {
    //edge_tl_id = cellB->getFacetId(fidB_plus1);
    edge_bl_id = cellB->getFacetId(fidB_plus2);

    edge_l = mesh->pushFacet(&edge_l_id);
    //edge_tl = mesh->getFacetPtr(edge_tl_id);
    edge_bl = mesh->getFacetPtr(edge_bl_id);
  }

  /* -------- All elements created, now setup adjacencies -----------*/

  // remeber: vtx_t, vtx_b, vtx_l, vtx_r, vtx_m, nod_t, nod_b, nod_l, nod_r, nod_tr, nod_tl, nod_br, nod_bl
  // void Point::setAllMembers(Real const* coord, int spacedim, int const* ic, int const* pos, int const* tag,
  //                            int const* flags, int const* stat, std::list<std::pair<int,char> > * incidences)
  // setting vtx_t: nothing to do
  // setting vtx_b:
  if ( ! vtx_b->replacesIncidCell(cidA, cidD, 2))
    vtx_b->replacesIncidCell(cidB, cidD, 2);
  // setting vtx_r: nothing to do
  // setting vtx_m:
  for (int i = 0; i < sdim; ++i)
    xyz_tmp[i] = t*(vt_xyz[i] - vb_xyz[i]) + vb_xyz[i];
  tag = edge_t->getTag();
  stat = edge_in_boundary ? Point::mk_inboundary : 0;
  vtx_m->setAllMembers(xyz_tmp, sdim, &cidA, &fidA_plus1, &tag, NULL/*flags*/, &stat, NULL/*singular list*/);
  vm_xyz = vtx_m->getCoord();
  // setting vtx_l: nothing to do

  if (high_order)
  {
    // setting nod_t: = (old edge node)
    for (int i = 0; i < sdim; ++i)
      xyz_tmp[i] = 0.5 * ( vt_xyz[i] + vm_xyz[i]);
    nod_t->setCoord(xyz_tmp, sdim);

    // setting nod_b:
    for (int i = 0; i < sdim; ++i)
      xyz_tmp[i] = 0.5 * ( vb_xyz[i] + vm_xyz[i]);
    tag = edge_t->getTag();
    pos = 4;
    stat = edge_in_boundary ? Point::mk_inboundary : 0;
    nod_b->setAllMembers(xyz_tmp, sdim, &cidD, &pos, &tag, NULL, &stat, NULL);

    //setting nod_r:
    for (int i = 0; i < sdim; ++i)
      xyz_tmp[i] = 0.5 * ( vr_xyz[i] + vm_xyz[i]);
    tag = cellA->getTag();
    pos = 3;
    stat = 0;
    nod_r->setAllMembers(xyz_tmp, sdim, &cidD, &pos, &tag, NULL, &stat, NULL);

    //setting nod_tr: nothing to do

    //setting nod_br:
    nod_br->setIncidence(cidD, 5);

    if (!edge_in_boundary)
    {
      //setting nod_bl
      nod_bl->setIncidence(cidC, 3);

      //setting nod_tl : nothing to do

      //setting nod_l:
      for (int i = 0; i < sdim; ++i)
        xyz_tmp[i] = 0.5 * ( vl_xyz[i] + vm_xyz[i]);
      tag = cellB->getTag();
      pos = 5;
      stat = 0;
      nod_l->setAllMembers(xyz_tmp, sdim, &cidC, &pos, &tag, NULL, &stat, NULL);

    }
  }

  //  ----- SetUp edges -------------------

  // remeber: edge_t, edge_b, edge_l, edge_r, edge_tr, edge_tl, edge_br, edge_bl
  // void Facet::setAllMembers(int const* ic, int const* pos, int const* tag, int const* flags, int const* bound_comp_id))

  //setting edge_t (= old edge) : nothing to do
  //setting edge_tr: nothing to do
  //setting edge_tl: nothing to do

  //setting edge_b:
  tag = edge_t->getTag();
  pos = 1; // in cidD
  bnd = edge_t->getBoundaryComponentId();
  edge_b->setAllMembers(&cidD, &pos, &tag, NULL, &bnd);

  //setting edge_r:
  tag = cellA->getTag();
  pos = 0; // in cidD
  bnd = -1;
  edge_r->setAllMembers(&cidD, &pos, &tag, NULL, &bnd);

  //setting edge_br:
  edge_br->setIncidence(cidD, 2);

  if (!edge_in_boundary)
  {
    //setting edge_l:
    tag = cellB->getTag();
    pos = 2; // in cidC
    bnd = -1;
    edge_l->setAllMembers(&cidC, &pos, &tag, NULL, &bnd);

    //setting edge_bl
    edge_bl->setIncidence(cidC, 0);
  }


  //  ----- SetUp Neighbors -------------------
  // top-left: nothing to do
  // top-right: nothing to do
  // bot-right
  Cell *neighbor = mesh->getCellPtr(cellA->getIncidCell(fidA_plus1));
  if (neighbor)
    neighbor->setIncidence(cellA->getIncidCellPos(fidA_plus1), cidD, 2);

  if (!edge_in_boundary)
  {
    // bot-left
    neighbor = mesh->getCellPtr(cellB->getIncidCell(fidB_plus2));
    if (neighbor)
      neighbor->setIncidence(cellB->getIncidCellPos(fidB_plus2), cidC, 0);
  }


  //  ----- SetUp Cells -------------------
  // remember: cellA, cellB, cellC, cellD
  //void setAllMembers(int const* nodes, int const* corners, int const* facets, int const* icells, int const* icells_pos,
  //                    int const* icells_ancs, int const* conn_comp_id, int const* tag, int const* flags)

  // set C and D first

  //setting cellD
  nds[0] = vtx_r_id;
  nds[1] = vtx_m_id;
  nds[2] = vtx_b_id;
  if (high_order)
  {
    nds[3] = nod_r_id;
    nds[4] = nod_b_id;
    nds[5] = nod_br_id;
  }
  fts[0] =  edge_r_id;
  fts[1] =  edge_b_id;
  fts[2] =  edge_br_id;
  ics[0] =  cidA;                             ics_pos[0] =  fidA_plus1;
  ics[1] =  cidC;                             ics_pos[1] =  1;
  ics[2] =  cellA->getIncidCell(fidA_plus1);  ics_pos[2] =  cellA->getIncidCellPos(fidA_plus1);
  bnd = cellA->getConnectedComponentId();
  tag = cellA->getTag();
  cellD->setAllMembers(nds, NULL, fts, ics, ics_pos, NULL, &bnd, &tag, NULL);


  if (!edge_in_boundary)
  {
    //setting cellC
    nds[0] = vtx_l_id;
    nds[1] = vtx_b_id;
    nds[2] = vtx_m_id;
    if (high_order)
    {
      nds[3] = nod_bl_id;
      nds[4] = nod_b_id;
      nds[5] = nod_l_id;
    }
    fts[0] =  edge_bl_id;
    fts[1] =  edge_b_id;
    fts[2] =  edge_l_id;
    ics[0] =  cellB->getIncidCell(fidB_plus2);  ics_pos[0] =  cellB->getIncidCellPos(fidB_plus2);
    ics[1] =  cidD;                             ics_pos[1] =  1;
    ics[2] =  cidB;                             ics_pos[2] =  fidB_plus2;
    bnd = cellB->getConnectedComponentId();
    tag = cellB->getTag();
    cellC->setAllMembers(nds, NULL, fts, ics, ics_pos, NULL, &bnd, &tag, NULL);
  }

  //setting cellA
  cellA->setNodeId(fidA_plus1, vtx_m_id);
  if (high_order)
    cellA->setNodeId(fidA_plus1+3, nod_r_id);
  cellA->setFacetId(fidA_plus1, edge_r_id);
  cellA->setIncidence(fidA_plus1, cidD, 0);

  //setting cellB
  if (!edge_in_boundary)
  {
    cellB->setNodeId(fidB, vtx_m_id);
    if (high_order)
      cellB->setNodeId(fidB_plus2+3, nod_l_id);
    cellB->setFacetId(fidB_plus2, edge_l_id);
    cellB->setIncidence(fidB_plus2, cidC, 2);
  }

  return vtx_m;
}

///** Safely removes a triangle (high orders to)
// *  @param cell the cell that will be removed.
// *  @param mesh mesh context.
// *  @return nothing.
// */
//void MeshToolsTri::removeTriCell(cell_handler cell)
//{
//  FEPIC_CHECK(cell.getPtr(), "invalid argument", std::runtime_error);
//
//  const int cell_dim = 2;
//  const int n_facets_per_cell = 3;
//
//  Mesh * mesh = cell.getMesh();
//  
//  const int cell_id_to_remove = cell.getIdx();
//  const int nds_per_cell = mesh->numNodesPerCell();
//  int       facets_id[n_facets_per_cell];
//  
//  // check what type of remove algo will be used .. see doc
//  
//  
//}
//
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
    x0 = mesh->getNodePtr(cell->getNodeId(0))->getCoord();
    x1 = mesh->getNodePtr(cell->getNodeId(1))->getCoord();
    x2 = mesh->getNodePtr(cell->getNodeId(2))->getCoord();

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
      cell = mesh->getCellPtr(cell->getIncidCell(l));
    else
      return std::pair<bool, Cell *>(false, const_cast<Cell *>(cell));


  }

}


































/*
 * 
 * 
 * 
 * 
 * 
 * 
 * backup
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 */














#if 0


Point* MeshToolsTri::insertVertexOnEdge(Cell *cellA, int fidA, Real t, Mesh *mesh)
{
/* ---------------------------------------------------------------------------------------------------------------
 *                      vtx_t                                                              vtx_t
 *                  ,/|\,                                                              ,/|\
 *                ,/  |  \,                                                          ,/  | `\
 *   cellB      ,/    |    \,    cellA                                  cellB      ,/    |   `\    cellA    <-- edge_t
 *            ,/      |      \,                                                  ,/      |vtx_m`\             = edge_old
 *          ,/        |        \,                    to                        ,/      _-+-_     `\
 *        ,/          | edge_    \,                     \,                   ,/   _,--'  |  '--,_  `\,
 * vtx_l /            | old        \ vtx_r                \,                /_,--'edge_l | edge_r`--,_\
 *       \,           |           ,/            ------------\               \,           |          ,/
 *         \,         |         ,/              -----------,/                 \,         |        ,/
 *           \,       |       ,/                         ,/                     \,       |       ,/         <-- edge_b
 *             \,     |     ,/                          /                         \,     |     ,/
 *               \,   |   ,/                                            cellC       \,   |   ,/    cellD
 *                 \, | ,/                                                            \, | ,/
 *                   \|/                                                                \|/
 *                     vtx_b                                                               vtx_b
 *
 * cellA nodes order:   j                                                 cellC nodes order:     2
 *                      | \ i                                                                0 / |
 *                      | /                                                                    \ |
 *                      k                                                                        1
 *
 * cellB nodes order:    k                                                cellD nodes order:  1
 *                   i / |                                                                    | \ 0
 *                     \ |                                                                    | /
 *                       j                                                                    2
 *
 * nod_tl, nod_tr, nod_bl, nod_br, nod_t*, nod_b, nod_l, nod_r = high order nodes
 *
 * nod_top = node of the old edge
 *
 *
 *///-------------------------------------------------------------------------------------------------------------

  int const  cidA = mesh->getCellId(cellA);
  
  FEPIC_CHECK((t<1 && t>0) && (cidA>=0 && cidA<mesh->numCellsTotal()), "t must be in ]0,1[", std::runtime_error);

  int const  fidA_plus1 = (fidA+1)%3;
  int const  fidA_plus2 = (fidA+2)%3;
  
  int const  cidB = cellA->getIncidCell(fidA);
  Real       xyz_tmp[3];
  int const  sdim = mesh->spaceDim();
  int        nds[6], fts[3], ics[3], ics_pos[3];
  int        pos, tag, stat, bnd;
  const bool edge_in_boundary = cidB < 0;
  const bool high_order = mesh->numNodesPerCell() > mesh->numVerticesPerCell();

  int        fidB = -1;
  //int        fidB_plus1 = -1;
  int        fidB_plus2 = -1;

  int        cidC = -1;
  int        cidD = -1;

  // creating cells to store informations ... the adjacencies
  // of the old cells are kept.
  Cell *cellB = NULL;
  Cell *cellC = NULL;
  Cell *cellD = mesh->pushCell(&cidD); //Cell *cellD = mesh->createCell(); cidD = mesh->pushCell(cellD);

  if (!edge_in_boundary)
  {
    fidB = cellA->getIncidCellPos(fidA);
    //fidB_plus1 = (fidB+1)%3;
    fidB_plus2 = (fidB+2)%3;

    cellB = mesh->getCellPtr(cidB);
    cellC = mesh->pushCell(&cidC); //cellC = mesh->createCell(); cidC = mesh->pushCell(cellC);
  }

  /* ---- First SetUp Lower dimension elements first ------  */
  // nodes First
  int const vtx_t_id = cellA->getNodeId(fidA      );
  int const vtx_b_id = cellA->getNodeId(fidA_plus1);
  int const vtx_r_id = cellA->getNodeId(fidA_plus2);
  int        vtx_l_id = -1;
  int        vtx_m_id;
  int        nod_t_id;
  int        nod_b_id;
  int        nod_r_id;
  int        nod_l_id;
  //int        nod_tr_id;
  //int        nod_tl_id;
  int        nod_br_id;
  int        nod_bl_id;

  Point *vtx_t = mesh->getNodePtr(vtx_t_id);
  Point *vtx_b = mesh->getNodePtr(vtx_b_id);
  Point *vtx_l = NULL;
  Point *vtx_r = mesh->getNodePtr(vtx_r_id);
  Point *vtx_m = mesh->pushPoint(&vtx_m_id);
  // high order nodes
  Point *nod_t = NULL;
  Point *nod_b = NULL;
  Point *nod_l = NULL;
  Point *nod_r = NULL;
  //Point *nod_tr = NULL;
  //Point *nod_tl = NULL;
  Point *nod_br = NULL;
  Point *nod_bl = NULL;
  Real const *vb_xyz = vtx_b->getCoord();
  Real const *vt_xyz = vtx_t->getCoord();
  Real const *vl_xyz = NULL;
  Real const *vr_xyz = vtx_r->getCoord();
  Real const *vm_xyz = NULL;

  if (!edge_in_boundary)
  {
    vtx_l_id = cellB->getNodeId(fidB_plus2);
    vtx_l = mesh->getNodePtr(vtx_l_id);
    vl_xyz = vtx_l->getCoord();
  }
  if (high_order)
  {
    nod_t_id  = cellA->getNodeId(fidA + 3); // node of the old edge
    //nod_tr_id = cellA->getNodeId(fidA_plus2 + 3);
    nod_br_id = cellA->getNodeId(fidA_plus1 + 3);

    nod_t  = mesh->getNodePtr(nod_t_id);
    //nod_tr = mesh->getNodePtr(nod_tr_id);
    nod_br = mesh->getNodePtr(nod_br_id);
    nod_b  = mesh->pushPoint(&nod_b_id);
    nod_r  = mesh->pushPoint(&nod_r_id);

    if (!edge_in_boundary)
    {
      //nod_tl_id = cellB->getNodeId(fidB_plus1 + 3);
      nod_bl_id = cellB->getNodeId(fidB_plus2 + 3);
      nod_l  = mesh->pushPoint(&nod_l_id);
      //nod_tl = mesh->getNodePtr(nod_tl_id);
      nod_bl = mesh->getNodePtr(nod_bl_id);
    }
  }

  /* -------- Edges ------  */

  int const edge_t_id  = cellA->getFacetId(fidA); // = old edge
  int       edge_b_id;
  int       edge_l_id  = -1;
  int       edge_r_id;
  //int const edge_tr_id = cellA->getFacetId(fidA_plus2);
  //int       edge_tl_id = -1;
  int const edge_br_id = cellA->getFacetId(fidA_plus1);
  int       edge_bl_id = -1;

  Facet *edge_t = mesh->getFacetPtr(edge_t_id);
  Facet *edge_b = mesh->pushFacet(&edge_b_id);
  Facet *edge_l = NULL;
  Facet *edge_r = mesh->pushFacet(&edge_r_id);;
  //Facet *edge_tr = mesh->getFacetPtr(edge_tr_id);
  //Facet *edge_tl = NULL;
  Facet *edge_br = mesh->getFacetPtr(edge_br_id);
  Facet *edge_bl = NULL;

  if (!edge_in_boundary)
  {
    //edge_tl_id = cellB->getFacetId(fidB_plus1);
    edge_bl_id = cellB->getFacetId(fidB_plus2);

    edge_l = mesh->pushFacet(&edge_l_id);
    //edge_tl = mesh->getFacetPtr(edge_tl_id);
    edge_bl = mesh->getFacetPtr(edge_bl_id);
  }

  /* -------- All elements created, now setup adjacencies -----------*/

  // remeber: vtx_t, vtx_b, vtx_l, vtx_r, vtx_m, nod_t, nod_b, nod_l, nod_r, nod_tr, nod_tl, nod_br, nod_bl
  // void Point::setAllMembers(Real const* coord, int spacedim, int const* ic, int const* pos, int const* tag,
  //                            int const* flags, int const* stat, std::list<std::pair<int,char> > * incidences)
  // setting vtx_t: nothing to do
  // setting vtx_b:
  if ( ! vtx_b->replacesIncidCell(cidA, cidD, 2))
    vtx_b->replacesIncidCell(cidB, cidD, 2);
  // setting vtx_r: nothing to do
  // setting vtx_m:
  for (int i = 0; i < sdim; ++i)
    xyz_tmp[i] = t*(vt_xyz[i] - vb_xyz[i]) + vb_xyz[i];
  tag = edge_t->getTag();
  stat = edge_in_boundary ? Point::mk_inboundary : 0;
  vtx_m->setAllMembers(xyz_tmp, sdim, &cidA, &fidA_plus1, &tag, NULL/*flags*/, &stat, NULL/*singular list*/);
  vm_xyz = vtx_m->getCoord();
  // setting vtx_l: nothing to do

  if (high_order)
  {
    // setting nod_t: = (old edge node)
    for (int i = 0; i < sdim; ++i)
      xyz_tmp[i] = 0.5 * ( vt_xyz[i] + vm_xyz[i]);
    nod_t->setCoord(xyz_tmp, sdim);

    // setting nod_b:
    for (int i = 0; i < sdim; ++i)
      xyz_tmp[i] = 0.5 * ( vb_xyz[i] + vm_xyz[i]);
    tag = edge_t->getTag();
    pos = 4;
    stat = edge_in_boundary ? Point::mk_inboundary : 0;
    nod_b->setAllMembers(xyz_tmp, sdim, &cidD, &pos, &tag, NULL, &stat, NULL);

    //setting nod_r:
    for (int i = 0; i < sdim; ++i)
      xyz_tmp[i] = 0.5 * ( vr_xyz[i] + vm_xyz[i]);
    tag = cellA->getTag();
    pos = 3;
    stat = 0;
    nod_r->setAllMembers(xyz_tmp, sdim, &cidD, &pos, &tag, NULL, &stat, NULL);

    //setting nod_tr: nothing to do

    //setting nod_br:
    nod_br->setIncidence(cidD, 5);

    if (!edge_in_boundary)
    {
      //setting nod_bl
      nod_bl->setIncidence(cidC, 3);

      //setting nod_tl : nothing to do

      //setting nod_l:
      for (int i = 0; i < sdim; ++i)
        xyz_tmp[i] = 0.5 * ( vl_xyz[i] + vm_xyz[i]);
      tag = cellB->getTag();
      pos = 5;
      stat = 0;
      nod_l->setAllMembers(xyz_tmp, sdim, &cidC, &pos, &tag, NULL, &stat, NULL);

    }
  }

  //  ----- SetUp edges -------------------

  // remeber: edge_t, edge_b, edge_l, edge_r, edge_tr, edge_tl, edge_br, edge_bl
  // void Facet::setAllMembers(int const* ic, int const* pos, int const* tag, int const* flags, int const* bound_comp_id))

  //setting edge_t (= old edge) : nothing to do
  //setting edge_tr: nothing to do
  //setting edge_tl: nothing to do

  //setting edge_b:
  tag = edge_t->getTag();
  pos = 1; // in cidD
  bnd = edge_t->getBoundaryComponentId();
  edge_b->setAllMembers(&cidD, &pos, &tag, NULL, &bnd);

  //setting edge_r:
  tag = cellA->getTag();
  pos = 0; // in cidD
  bnd = -1;
  edge_r->setAllMembers(&cidD, &pos, &tag, NULL, &bnd);

  //setting edge_br:
  edge_br->setIncidence(cidD, 2);

  if (!edge_in_boundary)
  {
    //setting edge_l:
    tag = cellB->getTag();
    pos = 2; // in cidC
    bnd = -1;
    edge_l->setAllMembers(&cidC, &pos, &tag, NULL, &bnd);

    //setting edge_bl
    edge_bl->setIncidence(cidC, 0);
  }


  //  ----- SetUp Neighbors -------------------
  // top-left: nothing to do
  // top-right: nothing to do
  // bot-right
  Cell *neighbor = mesh->getCellPtr(cellA->getIncidCell(fidA_plus1));
  neighbor->setIncidence(cellA->getIncidCellPos(fidA_plus1), cidD, 2);

  if (!edge_in_boundary)
  {
    // bot-left
    neighbor = mesh->getCellPtr(cellB->getIncidCell(fidB_plus2));
    neighbor->setIncidence(cellB->getIncidCellPos(fidB_plus2), cidC, 0);
  }


  //  ----- SetUp Cells -------------------
  // remember: cellA, cellB, cellC, cellD
  //void setAllMembers(int const* nodes, int const* corners, int const* facets, int const* icells, int const* icells_pos,
  //                    int const* icells_ancs, int const* conn_comp_id, int const* tag, int const* flags)

  // set C and D first

  //setting cellD
  nds[0] = vtx_r_id;
  nds[1] = vtx_m_id;
  nds[2] = vtx_b_id;
  if (high_order)
  {
    nds[3] = nod_r_id;
    nds[4] = nod_b_id;
    nds[5] = nod_br_id;
  }
  fts[0] =  edge_r_id;
  fts[1] =  edge_b_id;
  fts[2] =  edge_br_id;
  ics[0] =  cidA;                             ics_pos[0] =  fidA_plus1;
  ics[1] =  cidC;                             ics_pos[1] =  1;
  ics[2] =  cellA->getIncidCell(fidA_plus1);  ics_pos[2] =  cellA->getIncidCellPos(fidA_plus1);
  bnd = cellA->getConnectedComponentId();
  tag = cellA->getTag();
  cellD->setAllMembers(nds, NULL, fts, ics, ics_pos, NULL, &bnd, &tag, NULL);


  if (!edge_in_boundary)
  {
    //setting cellC
    nds[0] = vtx_l_id;
    nds[1] = vtx_b_id;
    nds[2] = vtx_m_id;
    if (high_order)
    {
      nds[3] = nod_bl_id;
      nds[4] = nod_b_id;
      nds[5] = nod_l_id;
    }
    fts[0] =  edge_bl_id;
    fts[1] =  edge_b_id;
    fts[2] =  edge_l_id;
    ics[0] =  cellB->getIncidCell(fidB_plus2);  ics_pos[0] =  cellB->getIncidCellPos(fidB_plus2);
    ics[1] =  cidD;                             ics_pos[1] =  1;
    ics[2] =  cidB;                             ics_pos[2] =  fidB_plus2;
    bnd = cellB->getConnectedComponentId();
    tag = cellB->getTag();
    cellC->setAllMembers(nds, NULL, fts, ics, ics_pos, NULL, &bnd, &tag, NULL);
  }

  //setting cellA
  cellA->setNodeId(fidA_plus1, vtx_m_id);
  if (high_order)
    cellA->setNodeId(fidA_plus1+3, nod_r_id);
  cellA->setFacetId(fidA_plus1, edge_r_id);
  cellA->setIncidence(fidA_plus1, cidD, 0);

  //setting cellB
  if (!edge_in_boundary)
  {
    cellB->setNodeId(fidB, vtx_m_id);
    if (high_order)
      cellB->setNodeId(fidB_plus2+3, nod_l_id);
    cellB->setFacetId(fidB_plus2, edge_l_id);
    cellB->setIncidence(fidB_plus2, cidC, 2);
  }

  return vtx_m;
}


#endif

