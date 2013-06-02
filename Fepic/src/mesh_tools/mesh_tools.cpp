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
#include <cmath>

//#include <type_traits>


// WARNING:
/*
 *  Antes de implementar uma operação sobre a malha, certifique-se
 *  operação mantém a consistência da mesma, verificando os seguintes
 *  itens:
 *
 *  - tag of all elements.
 *
 *  - nodes : m_icell
 *            m_icell_pos
 *            m_status
 *            m_incidences
 *
 *  - facets : m_icell
 *             m_icell_pos
 *             m_bound_comp_id
 *
 *  - corners : m_icell
 *              m_icell_pos
 *
 *  - cells : m_icells_pos
 *            m_icells_anchors
 *            m_facets
 *            m_icells
 *            m_corners
 *            m_nodes
 *            m_conn_comp_id
 *
 *  - SMesh : m_connected_comp_l // (connected component vs initial cell id) list
 *            m_boundary_comp_l  // (boundary component vs initialfacet id) list
 *
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
  //Cell const *ocell = mesh->getCellPtr(ocid);

  Point const *a = mesh->getNodePtr(cell->getNodeId(0));
  Point const *b = mesh->getNodePtr(cell->getNodeId(1));
  Point const *c = mesh->getNodePtr(cell->getNodeId(2));
  Point const *d = mesh->getNodePtr(cell->getNodeId((ofid+2)%3));

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
int MeshToolsTri::insertVertexOnEdge(int cell_A_id, int face_Am_id, Real t, Mesh *mesh)
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

  int const  sdim = mesh->spaceDim(); 		   // Pega a dimensão da malha
  
  Cell *const  cell_A = mesh->getCellPtr(cell_A_id);    // Pega o id da celula passada
  int const  cell_B_id = cell_A->getIncidCell(face_Am_id); // Pega id da outra celula que contem a aresta que será dividida
  
  FEPIC_CHECK((t<1 && t>0) && (cell_A_id>=0 && cell_A_id<mesh->numCellsTotal()), "t must be in ]0,1[", std::runtime_error);

  int const  face_Abr_id = (face_Am_id+1)%3; // Aresta seguinte à aresta que será dividida
  int const  face_Atr_id = (face_Am_id+2)%3; // Segunda aresta seguinte à aresta que será dividida
  
  // Variaveis auxiliares
  Real       xyz_tmp[3];                       // 
  
  int        nds[6], 
			 fts[3], 
			 ics[3], 
			 ics_pos[3]; //
  int        pos, 
			 tag, 
			 stat, 
			 bnd;
			 

  // testes para saber se a aresta é de bordo, e se os elementos da malha são de ordem superior
  const bool edge_in_boundary = cell_B_id < 0; // testa se a aresta a ser dividida é de bordo
  const bool highm_order = mesh->numNodesPerCell() > mesh->numVerticesPerCell(); // checa se o elemento é de ordem maior que 1

  int        face_Bm_id = -1;       // Id local da aresta que será dividida
  int        face_Bbl_id = -1; // Id local da segunda aresta seguinte à aresta que será dividida

  int        cell_C_id = -1;   // ids de uma das celulas que serão geradas
  int        cell_D_id = -1;   // id de uma das celulas que serão geradas

  // creating cells to store informations ... the adjacencies
  // of the old cells are kept.
  Cell *cell_B = NULL;   // Ponteiros que não se sabe se serão usadados ainda
  Cell *cell_C = NULL;
  Cell *cell_D = mesh->pushCell(&cell_D_id); // Ponteiro para a celula que será criada
  
  if (!edge_in_boundary)
  {// Se a aresta não for uma celula do bordo 
    face_Bm_id = cell_A->getIncidCellPos(face_Am_id); // Pega id local da aresta na celula oposta
    //fidB_plus1 = (fidB+1)%3;
    face_Bbl_id = (face_Bm_id+2)%3;  // segunda aresta seguinte na celula oposta

    cell_B = mesh->getCellPtr(cell_B_id); // Ponteiro para a celula oposta
    cell_C = mesh->pushCell(&cell_C_id); // Ponteiro para nova celula que será criada
  }

  /* ---- First SetUp Lower dimension elements first ------  */
  // nodes First
  // Os vertices por ordem a partir do primeiro vertice da aresta
  int const vtx_t_id = cell_A->getNodeId(face_Am_id ); 
  int const vtx_b_id = cell_A->getNodeId(face_Abr_id);
  int const vtx_r_id = cell_A->getNodeId(face_Atr_id);
  
  int        vtx_l_id = -1;
  int        vtx_m_id;
  // Armazenaria os nós caso seja um elemento de maior ordem
  int        nod_t_id;
  int        nod_b_id;
  int        nod_r_id;
  int        nod_l_id;
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
  //Point *nod_tr = NULL; // Não são alterados
  //Point *nod_tl = NULL; 
  Point *nod_br = NULL;
  Point *nod_bl = NULL;
  
  // Pega as coordenadas dos vertices
  Real const *vb_xyz = vtx_b->getCoord();
  Real const *vt_xyz = vtx_t->getCoord();
  Real const *vl_xyz = NULL; // vertice não exite caso a aresta seja de bordo
  Real const *vr_xyz = vtx_r->getCoord();
  Real const *vm_xyz = NULL; // vertice que será criado

  if (!edge_in_boundary)
  { // Se a aresta não for de bordo
    vtx_l_id = cell_B->getNodeId(face_Bbl_id); // Pega o Id do vertice oposto à aresta na celula oposta
    vtx_l = mesh->getNodePtr(vtx_l_id);      // Pega um ponteiro o vertice acima
    vl_xyz = vtx_l->getCoord();              // Pega as coordenadas do vertice
  }
  
  if (highm_order) // Se o elemento for de ordem superior
  {
    nod_t_id  = cell_A->getNodeId(face_Am_id + 3); // nó do meio da aresta que será dividida
    //nod_tr_id = cellA->getNodeId(fidA_plus2 + 3); ??????????????????????????? Não mexe pois não precisa
    nod_br_id = cell_A->getNodeId(face_Abr_id + 3); // Nó do meio da aresta entre vtx_b e vtx_r

    nod_t  = mesh->getNodePtr(nod_t_id); // Ponteiro para o nó do meio da aresta que será dividida
    //nod_tr = mesh->getNodePtr(nod_tr_id); ???????????????????????? Não mexe pois não precisa
    nod_br = mesh->getNodePtr(nod_br_id); // Ponteiro para o nó do meio da aresta entre vtx_b e vtx_r
    nod_b  = mesh->pushPoint(&nod_b_id); // 
    nod_r  = mesh->pushPoint(&nod_r_id); // 

    if (!edge_in_boundary)
    {// Se não for uma aresta de bordo
      //nod_tl_id = cellB->getNodeId(fidB_plus1 + 3);
      nod_bl_id = cell_B->getNodeId(face_Bbl_id + 3);
      nod_l  = mesh->pushPoint(&nod_l_id);
      //nod_tl = mesh->getNodePtr(nod_tl_id);
      nod_bl = mesh->getNodePtr(nod_bl_id);
    }
  }

  // ================= Recolhe informação das edges ====================================================================
  
  // Ids das arestas
  int const edge_t_id  = cell_A->getFacetId(face_Am_id); // = old edge
  int       edge_b_id; // Aresta que segue o vertice vtx_b
  int       edge_l_id  = -1; // Aresta que segue o vtx_l
  int       edge_r_id; // Aresta que segue o vertice vtx_r
  int const edge_br_id = cell_A->getFacetId(face_Abr_id); // Aresta entre o vertices vtx_b e vtx_r
  int       edge_bl_id = -1; // Aresta entre os vertices vtx_b e vtx_l

  // Ponteiros para as arestas
  Facet *edge_t = mesh->getFacetPtr(edge_t_id); 	// 
  Facet *edge_b = mesh->pushFacet(&edge_b_id);  	// 
  Facet *edge_l = NULL;								// 
  Facet *edge_r = mesh->pushFacet(&edge_r_id);; 	// 
  Facet *edge_br = mesh->getFacetPtr(edge_br_id);	// 
  Facet *edge_bl = NULL;							// 

  if (!edge_in_boundary)
  { // Se a aresta não for da borda, então o vertice vtx_l existe a as arestas insidentes a este
    edge_bl_id = cell_B->getFacetId(face_Bbl_id);

    edge_l = mesh->pushFacet(&edge_l_id);
    //edge_tl = mesh->getFacetPtr(edge_tl_id);
    edge_bl = mesh->getFacetPtr(edge_bl_id);
  }

  // ================================ All elements created, now setup adjacencies ======================================

  // remeber: vtx_t, vtx_b, vtx_l, vtx_r, vtx_m, nod_t, nod_b, nod_l, nod_r, nod_tr, nod_tl, nod_br, nod_bl
  // void Point::setAllMembers(Real const* coord, int spacedim, int const* ic, int const* pos, int const* tag,
  //                            int const* flags, int const* stat, std::list<std::pair<int,char> > * incidences)
  // setting vtx_t: nothing to do
  // setting vtx_b:
  
  if ( ! vtx_b->replacesIncidCell(cell_A_id, cell_D_id, 2) )
    vtx_b->replacesIncidCell(cell_B_id, cell_D_id, 2);  // !!!!!!!!!!!!!!!!!!!!!!!!!!!! aqui acho que era C no lugar de D
  // setting vtx_r: nothing to do
  // setting vtx_m: 
  
  // Calculo da posição onde deve ser inserido um novo ponto
  for (int i = 0; i < sdim; ++i)
    xyz_tmp[i] = t*(vt_xyz[i] - vb_xyz[i]) + vb_xyz[i];
  
  // Guarda a Tag associada à aresta dividida
  tag = edge_t->getTag();
  
  // Guarda se a aresta é de bordo em uma variavel de estado 
  stat = edge_in_boundary ? Point::mk_inboundary : 0;
  
  // Seta os valores referentes ao vertice criado
  vtx_m->setAllMembers(xyz_tmp, sdim, &cell_A_id, &face_Abr_id, &tag, NULL/*flags*/, &stat, NULL/*singular list*/);
  
  // Pega as coordenadas deste vertice novo
  vm_xyz = vtx_m->getCoord();
  // setting vtx_l: nothing to do

  if (highm_order)
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
    nod_b->setAllMembers(xyz_tmp, sdim, &cell_D_id, &pos, &tag, NULL, &stat, NULL);

    //setting nod_r:
    for (int i = 0; i < sdim; ++i)
      xyz_tmp[i] = 0.5 * ( vr_xyz[i] + vm_xyz[i]);
    tag = cell_A->getTag();
    pos = 3;
    stat = 0;
    nod_r->setAllMembers(xyz_tmp, sdim, &cell_D_id, &pos, &tag, NULL, &stat, NULL);

    //setting nod_tr: nothing to do

    //setting nod_br:
    nod_br->setIncidence(cell_D_id, 5);

    if (!edge_in_boundary)
    {
      //setting nod_bl
      nod_bl->setIncidence(cell_C_id, 3);

      //setting nod_tl : nothing to do

      //setting nod_l:
      for (int i = 0; i < sdim; ++i)
        xyz_tmp[i] = 0.5 * ( vl_xyz[i] + vm_xyz[i]);
      tag = cell_B->getTag();
      pos = 5;
      stat = 0;
      nod_l->setAllMembers(xyz_tmp, sdim, &cell_C_id, &pos, &tag, NULL, &stat, NULL);

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
 
  edge_b->setAllMembers(&cell_D_id, &pos, &tag, NULL, &bnd);

  //setting edge_r:
  tag = cell_A->getTag();
  pos = 0; // in cidD
  bnd = -1;
  
  edge_r->setAllMembers(&cell_D_id, &pos, &tag, NULL, &bnd);

  //setting edge_br:
  edge_br->setIncidence(cell_D_id, 2);

  if (!edge_in_boundary)
  {
    //setting edge_l:
    tag = cell_B->getTag();
    pos = 2; // in cidC
    bnd = -1;
    edge_l->setAllMembers(&cell_C_id, &pos, &tag, NULL, &bnd);

    //setting edge_bl
    edge_bl->setIncidence(cell_C_id, 0);
  }


  //=======================================  ----- SetUp Neighbors --------------  =====================================
  // top-left: nothing to do
  // top-right: nothing to do
  // bot-right
  
  Cell *neighbor = mesh->getCellPtr(cell_A->getIncidCell(face_Abr_id));
  if (neighbor)
    neighbor->setIncidence(cell_A->getIncidCellPos(face_Abr_id), cell_D_id, 2);

  if (!edge_in_boundary)
  {
    // bot-left
    neighbor = mesh->getCellPtr(cell_B->getIncidCell(face_Bbl_id));
    if (neighbor)
      neighbor->setIncidence(cell_B->getIncidCellPos(face_Bbl_id), cell_C_id, 0);
  }

  // ===================================== ----- SetUp Cells ------------------- =======================================
  // remember: cellA, cellB, cellC, cellD
  //void setAllMembers(int const* nodes, int const* corners, int const* facets, int const* icells, int const* icells_pos,
  //                    int const* icells_ancs, int const* conn_comp_id, int const* tag, int const* flags)

  // set C and D first

  // ============ setting cellD
  nds[0] = vtx_r_id;
  nds[1] = vtx_m_id;
  nds[2] = vtx_b_id;
  
  if (highm_order)
  {
    nds[3] = nod_r_id;
    nds[4] = nod_b_id;
    nds[5] = nod_br_id;
  }
  
  fts[0] =  edge_r_id;
  fts[1] =  edge_b_id;
  fts[2] =  edge_br_id;
  
  ics[0] =  cell_A_id;                             ics_pos[0] =  face_Abr_id;
  ics[1] =  cell_C_id;                             ics_pos[1] =  1;
  ics[2] =  cell_A->getIncidCell(face_Abr_id);     ics_pos[2] =  cell_A->getIncidCellPos(face_Abr_id);
  
  bnd = cell_A->getConnectedComponentId();
  tag = cell_A->getTag();
  
  cell_D->setAllMembers(nds, NULL, fts, ics, ics_pos, NULL, &bnd, &tag, NULL);

  if (!edge_in_boundary)
  {
    // ================================ setting cellC ==================================================================
    nds[0] = vtx_l_id;
    nds[1] = vtx_b_id;
    nds[2] = vtx_m_id;
    
    if (highm_order)
    {
      nds[3] = nod_bl_id;
      nds[4] = nod_b_id;
      nds[5] = nod_l_id;
    }
    fts[0] =  edge_bl_id;
    fts[1] =  edge_b_id;
    fts[2] =  edge_l_id;
    ics[0] =  cell_B->getIncidCell(face_Bbl_id);  ics_pos[0] =  cell_B->getIncidCellPos(face_Bbl_id);
    ics[1] =  cell_D_id;                               ics_pos[1] =  1;
    ics[2] =  cell_B_id;                               ics_pos[2] =  face_Bbl_id;
    bnd = cell_B->getConnectedComponentId();
    tag = cell_B->getTag();
    cell_C->setAllMembers(nds, NULL, fts, ics, ics_pos, NULL, &bnd, &tag, NULL);
  }

  // ======================================== setting cellA ============================================================
  cell_A->setNodeId(face_Abr_id, vtx_m_id);
  if (highm_order)
    cell_A->setNodeId(face_Abr_id+3, nod_r_id);
  cell_A->setFacetId(face_Abr_id, edge_r_id);
  cell_A->setIncidence(face_Abr_id, cell_D_id, 0);

  // ========================================= setting cellB ===========================================================
  if (!edge_in_boundary)
  {
    cell_B->setNodeId(face_Bm_id, vtx_m_id);
    if (highm_order)
      cell_B->setNodeId(face_Bbl_id+3, nod_l_id);
    cell_B->setFacetId(face_Bbl_id, edge_l_id);
    cell_B->setIncidence(face_Bbl_id, cell_C_id, 2);
  }

  return vtx_m_id;
}

// =====================================================================================================================


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
//  Mesh * mesh = cell.meshPtr();
//  
//  const int cell_id_to_remove = cell.index();
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


Real MeshToolsTri::edgeSize(Cell const* cell, const int fid, Mesh const* mesh)
{
	  Point const *a = mesh->getNodePtr(cell->getNodeId(fid));
	  Point const *b = mesh->getNodePtr(cell->getNodeId((fid+1)%3));
	  
	  Real const ax = a->getCoord(0), ay = a->getCoord(1);
	  Real const bx = b->getCoord(0), by = b->getCoord(1);
	  
	  Real const size = sqrt( (ax-bx)*(ax-bx) + (ay-by)*(ay-by) );
	  
	  return size;
}


Real MeshToolsTri::meanRatio(Cell const* cell, Mesh const* mesh)
{
	Point const *a = mesh->getNodePtr(cell->getNodeId(0));
	Point const *b = mesh->getNodePtr(cell->getNodeId(1));
	Point const *c = mesh->getNodePtr(cell->getNodeId(2));
	
	Real const ax = a->getCoord(0), ay = a->getCoord(1);
	Real const bx = b->getCoord(0), by = b->getCoord(1);
	Real const cx = c->getCoord(0), cy = c->getCoord(1);
	
	Real const sizee = MeshToolsTri::edgeSize(cell, 0, mesh);
	Real const sizef = MeshToolsTri::edgeSize(cell, 1, mesh);
	Real const sizeg = MeshToolsTri::edgeSize(cell, 2, mesh);
	
	Real At = 0.5 * ( ax*by + ay*cx + bx*cy - by*cx - ay*bx - ax*cy );
	
	
	Real ratio = 4* sqrt(3) * At / ( sizee*sizee + sizef*sizef + sizeg*sizeg );
	
	return ratio;	
}

void MeshToolsTri::centroid(Cell const* cell, Real *coord_ctrd, Mesh const* mesh)
{
	Point const *a = mesh->getNodePtr(cell->getNodeId(0));
	Point const *b = mesh->getNodePtr(cell->getNodeId(1));
	Point const *c = mesh->getNodePtr(cell->getNodeId(2));
	
	Real const ax = a->getCoord(0), ay = a->getCoord(1);
	Real const bx = b->getCoord(0), by = b->getCoord(1);
	Real const cx = c->getCoord(0), cy = c->getCoord(1);
	
	coord_ctrd[0]=(1.0/3.0) * (ax + bx + cx);
	coord_ctrd[1]=(1.0/3.0) * (ay + by + cy);
		
}

Real MeshToolsTri::area(Cell const* cell, Mesh const* mesh)
{
	Point const *a = mesh->getNodePtr(cell->getNodeId(0));
	Point const *b = mesh->getNodePtr(cell->getNodeId(1));
	Point const *c = mesh->getNodePtr(cell->getNodeId(2));
	
	Real const ax = a->getCoord(0), ay = a->getCoord(1);
	Real const bx = b->getCoord(0), by = b->getCoord(1);
	Real const cx = c->getCoord(0), cy = c->getCoord(1);
	
	Real At = 0.5 * ( ax*by + ay*cx + bx*cy - by*cx - ay*bx - ax*cy );
	
	return At;	
}

int MeshToolsTri::collapseEdge2d(int cell_A_id, int face_Am_id, Real t, Mesh *mesh)
{ 
	int const  sdim = mesh->spaceDim();
	
	// ============================== Ids das celulas opostas à aresta =================================================
	Cell const *cell_A = mesh->getCellPtr(cell_A_id);
	int	const cell_B_id = cell_A->getIncidCell(face_Am_id);
	Cell *cell_B = NULL;
	
	// ================================ Teste se a aresta é de bordo ===================================================
	const bool edge_in_boundary = cell_B_id < 0;
	
	// ================================= Ainda não tratei este caso ====================================================
	//const bool highm_order = mesh->numNodesPerCell() > mesh->numVerticesPerCell();
	
	// =============================== Ids locais das arestas das celulas ==============================================
	int face_Abr_id = (face_Am_id+1)%3,
		face_Atr_id = (face_Am_id+2)%3,
		face_Bbl_id = -1,
		face_Btl_id = -1,
		face_Bm_id   = -1;
		
	if ( !edge_in_boundary )
	{ // Trata o caso onde a aresta não é de bordo
		cell_B 		 = mesh->getCellPtr(cell_B_id);  // Ponteiro para a celula B
		face_Bm_id 	 = cell_A->getIncidCellPos(face_Am_id);
		face_Btl_id  = (face_Bm_id+1)%3;
		face_Bbl_id  = (face_Bm_id+2)%3;
	}
	
	// =============================== Ids dos vertices ================================================================
	int const vtx_t_id = cell_A->getNodeId(face_Am_id); 
	int const vtx_b_id = cell_A->getNodeId(face_Abr_id);
	int const vtx_r_id = cell_A->getNodeId(face_Atr_id);
	int 	  vtx_l_id = -1;
	if ( !edge_in_boundary )
	{ // Trata o caso onde a aresta não é de bordo
		vtx_l_id = cell_B->getNodeId(face_Bbl_id);
	}
	
	// ================================ Id global das arestas ==========================================================
	int const edge_m_id  = cell_A->getFacetId(face_Am_id); // = old edge
	int       edge_br_id = cell_A->getFacetId(face_Abr_id); // Aresta que segue o vertice vtx_b
	int       edge_tr_id = cell_A->getFacetId(face_Atr_id); // Aresta que segue o vtx_l
	int       edge_bl_id = -1; // Aresta que segue o vertice vtx_r
	int 	  edge_tl_id = -1; // Aresta entre o vertices vtx_b e vtx_r
	
	if ( !edge_in_boundary )
	{ // Trata o caso onde a aresta não é de bordo
		edge_bl_id = cell_B->getFacetId(face_Bbl_id);
		edge_tl_id = cell_B->getFacetId(face_Btl_id);
	}
		
	
	// ========================== Guada o Id global das celulas vizinhas ===============================================
	int const cell_AT_id = cell_A->getIncidCell(face_Atr_id),
			  cell_AB_id = cell_A->getIncidCell(face_Abr_id);
	int 	  cell_BT_id = -1,
			  cell_BB_id = -1;
	if ( !edge_in_boundary )
	{ // Trata o caso onde a aresta não é de bordo
		cell_BT_id = cell_B->getIncidCell(face_Btl_id);
		cell_BB_id = cell_B->getIncidCell(face_Bbl_id);
	}
	
	// =========================== Guarda o Id local das faces nas celulas vizinhas ====================================
	int const face_ATtr_id = cell_A->getIncidCellPos(face_Atr_id),
			  face_ABbr_id = cell_A->getIncidCellPos(face_Abr_id);
	int 	  face_BTtl_id = -1,
			  face_BBbl_id = -1;
	if ( !edge_in_boundary )
	{ // Trata o caso onde a aresta não é de bordo
		face_BTtl_id = cell_B->getIncidCellPos(face_Btl_id);
		face_BBbl_id = cell_B->getIncidCellPos(face_Bbl_id);	
	}	
	
	// TO DO
    // ================================ Trata as componentes conexa e de bordo =========================================
    // int bdr=mesh->getFacetPtr(edge_tr_id)->getBoundaryComponentId();
    
	// ================================ Mover os vertices vtx_t e vtx_b para a posição de colapso ======================
	Point *vtx_t = mesh->getNodePtr(vtx_t_id),
          *vtx_b = mesh->getNodePtr(vtx_b_id),
          *vtx_r = mesh->getNodePtr(vtx_r_id),
          *vtx_l = NULL;
          
    if( !edge_in_boundary) vtx_l = mesh->getNodePtr(vtx_l_id);
    
    Real const *coords_t = vtx_t->getCoord(), 
			   *coords_b = vtx_b->getCoord();
    Real coords_m[3];
    
    for (int i = 0; i < sdim; ++i)
		coords_m[i] = t*(coords_t[i] - coords_b[i]) + coords_b[i];
    
	vtx_t->setCoord(coords_m,sdim);
	vtx_b->setCoord(coords_m,sdim);
	
	// ================================ Pegar a estrela do vertice que vai ser substituido =============================
    int iCs_t[64],     // Ids das celulas da estrela 
        viCs_t[64];  // Ids locais do vertice nas celulas da estrela
     
    mesh->vertexStar(cell_A_id, face_Am_id , iCs_t, viCs_t);
	
	// Indica a todas as celulas que continham vtx_t que vtx_b está na posição agora
    for( int i=0; iCs_t[i] > -1; i++)
    {
	   Cell *cell = mesh->getCellPtr(iCs_t[i]);
	   cell->setNodeId(viCs_t[i], vtx_b_id); 	
    }
    
    // =========================== Atribuir as adjacencias aos vertices ================================================
    // Vertice vtx_b e vtx_r
    if( cell_AT_id > -1 )
    {
		if(! vtx_b->replacesIncidCell(cell_A_id, cell_AT_id, face_ATtr_id))
		{			
			if( cell_BT_id > -1 )
				vtx_b->replacesIncidCell(cell_B_id, cell_BT_id, (face_BTtl_id+1)%3 );
			else
				vtx_b->replacesIncidCell(cell_B_id, cell_BB_id, (face_BBbl_id) );
		}
		
		vtx_r->replacesIncidCell(cell_A_id, cell_AT_id, (face_ATtr_id+1)%3 );
	}else
	{
		if(! vtx_b->replacesIncidCell(cell_A_id, cell_AB_id, (face_ABbr_id+1)%3))
		{
			if( cell_BT_id > -1 )
				vtx_b->replacesIncidCell(cell_B_id, cell_BT_id, (face_BTtl_id+1)%3);
			else
				vtx_b->replacesIncidCell(cell_B_id, cell_BB_id, (face_BBbl_id)%3 );
		}
		vtx_r->replacesIncidCell(cell_A_id, cell_AB_id, (face_ABbr_id) );
	}	    
    
    if ( !edge_in_boundary ) 
		vtx_l->replacesIncidCell(cell_B_id, cell_BT_id, face_BTtl_id);        
    
    // =========================== Atribuir as adjacencias das arestas =================================================
	
	Facet *edge_br = mesh->getFacetPtr(edge_br_id);
		  //*edge_tr = mesh->getFacetPtr(edge_tr_id) ;	
	Facet *edge_bl = NULL;
		  //*edge_tl = NULL;						
    
    // Aresta edge_br
    if(cell_AT_id > -1 ) edge_br->setIncidence(cell_AT_id, face_ATtr_id );
    else  				 edge_br->setIncidence(cell_AB_id, face_ABbr_id );
    
    // Aresta edge_bl
    if ( !edge_in_boundary ) 
    {
		edge_bl=mesh->getFacetPtr(edge_bl_id);
		// edge_tl=mesh->getFacetPtr(edge_tl_id);
		edge_bl->setIncidence(cell_BT_id, face_BTtl_id);
	}
    
    // ======= Fala para as celulas vizinhas o id das arestas que são suas faces, e quais são suas celulas vizinhas ====
    // =======================================  ----- SetUp Neighbors --------------  ==================================
  
	if ( cell_AT_id > -1 )
	{
		Cell *cell_AT = mesh->getCellPtr(cell_AT_id);
		cell_AT->setIncidence(face_ATtr_id, cell_AB_id, face_ABbr_id);
		cell_AT->setFacetId(face_ATtr_id, edge_br_id);
	}
	if ( cell_AB_id > -1 )
	{
		Cell *cell_AB = mesh->getCellPtr(cell_AB_id);
		cell_AB->setIncidence(face_ABbr_id, cell_AT_id, face_ATtr_id);
	}

	if (!edge_in_boundary)
	{
		if(cell_BT_id > -1 )
		{
			Cell *cell_BT = mesh->getCellPtr(cell_BT_id);
			cell_BT->setFacetId(face_BTtl_id, edge_bl_id);
			cell_BT->setIncidence(face_BTtl_id, cell_BB_id, face_BBbl_id);
		}
		if(cell_BB_id > -1 )
		{
			Cell *cell_BB = mesh->getCellPtr(cell_BB_id);
			cell_BB->setIncidence(face_BBbl_id, cell_BT_id, face_BTtl_id);
		}	
	}
  
  // ========================================== Caso de vtx_t estar no bordo ===========================================
  // Se vtx_t esta no bordo, vtx_b tambem estara
  if (mesh->inBoundary(vtx_t))
    vtx_b->setAsBoundary(true);
    
    // ================================ Deletar os elementos ===========================================================
	mesh->disablePoint(vtx_t_id);
	
	mesh->disableFacet(edge_tr_id);
	mesh->disableFacet(edge_m_id);
	if(!edge_in_boundary)
		mesh->disableFacet(edge_tl_id);
	
	mesh->disableCell(cell_A_id);
	if(!edge_in_boundary) 
		mesh->disableCell(cell_B_id);   
    
  // Casos degenerados
  if( ( cell_AT_id == -1 )&&( cell_AB_id ==-1 ) )
  {
		mesh->disableFacet(edge_br_id);
		mesh->disablePoint(vtx_r_id);
	}
	if( ( cell_BT_id == -1 )&&( cell_BB_id ==-1 )&&(edge_bl_id!=-1))
  {
		mesh->disableFacet(edge_bl_id);
		mesh->disablePoint(vtx_l_id);
	}
			
    return vtx_b_id; 	
}

int MeshToolsTetr::calcAnchors(int cellA_id, int facetA_id, int cellB_id, int facetB_id, Mesh *mesh)
{
	int const nVtxPerFacet = mesh->numVerticesPerFacet();
	int	*idsA=new int[nVtxPerFacet], 
		*idsB=new int[nVtxPerFacet];
	mesh->getCellPtr(cellA_id)->getFacetVerticesId(facetA_id, idsA);
	mesh->getCellPtr(cellB_id)->getFacetVerticesId(facetB_id, idsB);
	
	int anch = calcAnchors( idsA, idsB, mesh);
	delete idsA; delete idsB;
	return anch;	
}

int MeshToolsTetr::calcAnchors(int *idsA, int *idsB, Mesh *mesh)
{
	int n=mesh->numVerticesPerFacet(), 
	    idsB_temp[4], idsA_temp[4], idsA_aux[4];
	
	for( int i = 0; i< n; i++)
	{
		idsB_temp[i] = idsB[n-1-i];
		idsA_temp[i] = idsA[i];
	}		
	
	for( int i = 0; i< n; i++) 
	{
		bool foundAnchor = false;
		for( int j = 0; j< n; j++)
		{
			if( idsB_temp[j] == idsA_temp[j] )
				foundAnchor = true;
			else
			{
				foundAnchor = false;
				break;
			}
		}
		
		if( foundAnchor == true ) 
			return i;
		else
		{ // Rotação a direita
			for( int j = 0; j< n; j++) idsA_aux[j] = idsA_temp[j];
			for( int j = 0; j< n; j++) idsA_temp[(j+1)%n] = idsA_aux[j];
		}
	}
	return -1;
}

int MeshToolsTetr::insertVertexOnEdge3d(int cell_A_id, int corner_Am_pos, Real t, Mesh *mesh)
{
	int sdim=mesh->spaceDim();
	// ====== Testa se a aresta é de bordo e se a celula é de ordem superior ============= ok
	Cell *cell_A = mesh->getCellPtr(cell_A_id);
	int corner_m_id = cell_A->getCornerId(corner_Am_pos);
	
	const bool edge_in_boundary = mesh->inBoundary( mesh->getCornerPtr(cell_A->getCornerId(corner_Am_pos)) ) ; // testa se a aresta a ser dividida é de bordo
    //const bool highm_order = mesh->numNodesPerCell() > mesh->numVerticesPerCell();      // checa se o elemento é de ordem maior que 1
		
	// ============================ Guarda a estrela da aresta =========================== ok
	int idsC[32], 					// Id das células da estrela da aresta
		idsC_e[32], 				// Id local da aresta nas células da estrela
		nStar=0; 					// número de células na estrela.
		
	mesh->edgeStar(mesh->getCornerPtr(corner_m_id)->getIncidCell(), mesh->getCornerPtr(corner_m_id)->getPosition(), idsC, idsC_e);
	
	for( int i=0; idsC[i]!=-1; i++){
		nStar++;
	}
	
	// ================= Cria um nó no meio da aresta ==================================== ok
	int vtx_m_id;
	mesh->pushPoint(&vtx_m_id);
	
	// ======= Para cada célula da estrela cria uma célula que será sua vizinha ========== ok
	int *idsC_neighbor=new int[nStar];	
	
	for( int i=0; i<nStar; i++)
	{
		mesh->pushCell(&idsC_neighbor[i]);
	}
		
	// ========== Guarda os ids dos nós da aresta========================================= ok
	// A ordem destes pontos vai descidir qual o resultado final da adaptação, a células 
	// novas conterão o primeiro nó, e as células mantida conterão o segundo ============= ok
	int nNodesPerCorner = mesh->numNodesPerCorner(), 
	    *edgeNodes=new int[nNodesPerCorner];
	mesh->getCellPtr(idsC[0])->getCornerNodesId(idsC_e[0], edgeNodes);
	
	int vtx_t_id = edgeNodes[0],
		vtx_b_id = edgeNodes[1];

// ========== DEFINE TODAS AS ADJACENCIAS DO PONTO M ===================
	Point *vtx_t = mesh->getNodePtr(vtx_t_id),
		  *vtx_b = mesh->getNodePtr(vtx_b_id),
		  *vtx_m = mesh->getNodePtr(vtx_m_id);
	int tag, stat, pos=0;
		  
	Real coords_t[3], coords_b[3], coords_m[3];
	vtx_t->getCoord(coords_t,mesh->spaceDim());
	vtx_b->getCoord(coords_b,mesh->spaceDim());
	
	for(int i=0; i<3; i++) coords_m[i]=coords_t[i]+t*(coords_b[i]-coords_t[i]);
  
    // Guarda a Tag associada à aresta dividida
    tag = mesh->getCornerPtr(corner_m_id)->getTag();
  
    // Guarda se a aresta é de bordo em uma variavel de estado 
    stat = edge_in_boundary ? Point::mk_inboundary : 0;
  
    // Seta os valores referentes ao vertice criado
    // a posição sera deixada como zero, mas será atualizada quando a primeira célula
    // da estrela for avaliada.
    vtx_m->setAllMembers(coords_m, sdim, &idsC[0], &pos, &tag, NULL/*flags*/, &stat, NULL/*singular list*/);		
  
// ============= Cria aresta que será parte da aresta dividida ====================== ok
	int vertCorner_id; 
	mesh->pushCorner(&vertCorner_id);
	  
// ========== Cria as novas arestas ================================================== ok
	// Se a aresta a ser dividida é de bordo cria ( nStar+1 ) arestas, senão cria ( nStar )
	// arestas.
	int nCorners = edge_in_boundary ? (nStar+1) : (nStar),
		*horiCorners_ids=new int[nCorners];	
	
	for( int i=0 ; i< nCorners; i++){
		 mesh->pushCorner(&horiCorners_ids[i]);	
	}
	
// ========== Cria as novas faces =================================================== ok
	// Se a aresta a ser dividida é de bordo cria ( nStar+1 ) faces verticais, 
	// senão cria ( nStar ) faces verticais. Cria nStar faces horizontais
	int *vertFaces_ids=new int[nCorners], 
	    *horiFaces_ids=new int[nStar];
	    
	for( int i=0 ; i<nCorners; i++){
		 mesh->pushFacet(&vertFaces_ids[i]);	
	}
	
	for( int i=0 ; i<nStar ; i++){
		 mesh->pushFacet(&horiFaces_ids[i]);	
	}	
	
// ***********************************************************************************
// ====== Variaveis que serão usadas nas atribuições ================================= ok
	int corner_m_pos, vtx_t_pos, vtx_b_pos;
		corner_m_pos = idsC_e[0];
		
	mesh->getCellPtr(idsC[0])->getCornerNodesId(corner_m_pos, edgeNodes);
		
	if( edgeNodes[0] == vtx_t_id )
	{
		vtx_t_pos = Tetrahedron4::m_table_bC_x_vC[corner_m_pos][0];
		vtx_b_pos = Tetrahedron4::m_table_bC_x_vC[corner_m_pos][1];
	}else
	{
		vtx_t_pos = Tetrahedron4::m_table_bC_x_vC[corner_m_pos][1];
		vtx_b_pos = Tetrahedron4::m_table_bC_x_vC[corner_m_pos][0];
	}	
	
	int *nds=new int[mesh->numNodesPerCell()],   					// Id dos nós da célula
		*fts=new int[mesh->numFacetsPerCell()],						// Id das faces da célula
		*cns=new int[mesh->numCornersPerCell()],						// Id dos cantos da célula
		*ics=new int[mesh->numFacetsPerCell()] , 					// Id das células vizinhas
		*ics_pos=new int[mesh->numFacetsPerCell()], 					// Posição nas células vizinhas
		*ics_anc=new int[mesh->numFacetsPerCell()]; 					// Ancora nas células vizinhas
	
	int vtx_l_pos, 
		vtx_r_pos, 
	    corner_tl_pos, 
	    corner_tr_pos, 
	    corner_bl_pos, 
	    corner_br_pos,
	    corner_op_pos;
	getOpositVertexsPos(vtx_t_pos,vtx_b_pos, &vtx_r_pos, &vtx_l_pos);
	getOpositCornersPos(vtx_t_pos,vtx_b_pos, &corner_tr_pos,&corner_tl_pos,&corner_br_pos, &corner_bl_pos, &corner_op_pos);
	int facet_r_pos=getOpositFacetPos(vtx_l_pos),
	    facet_l_pos=getOpositFacetPos(vtx_r_pos),
	    facet_t_pos=getOpositFacetPos(vtx_b_pos),
	    facet_b_pos=getOpositFacetPos(vtx_t_pos);
	
// ===================================================================================
// Tratamento especial para o caso de arestas no bordo onde a estrela é dada na ordem
// inversa ===========================================================================

	if( (mesh->getCellPtr(idsC[0])->getIncidCell(facet_l_pos)!=-1) && (mesh->getCellPtr(idsC[0])->getIncidCell(facet_l_pos) == idsC[1]))
	{ // neste caso a estrela será dada em ordem inversa, e para arrumar a orientação trocamos
	  // a orientação da aresta pelos seus vertice t e b
		int aux=vtx_t_id;
		Point *vtx_aux = vtx_t;
		vtx_t_id = vtx_b_id;  	vtx_t = vtx_b;
		vtx_b_id = aux;			vtx_b = vtx_aux;
		std::cerr << "caso da estrela invertida " << mesh->getCellPtr(idsC[0])->getIncidCell(facet_l_pos) << std::endl;
	}

// ITERA SOBRE TODAS AS CÉLULAS DA ESTRELA DA ARESTA==================== ok
	for( int i=0; i<nStar; i++) 
	{
// ==================== Células envolvidas ===================== ok
		Cell *cell_old = mesh->getCellPtr(idsC[i]),
		     *cell_new = mesh->getCellPtr(idsC_neighbor[i]);
		     
		int cell_old_id = idsC[i],
		    cell_new_id = idsC_neighbor[i];
		
// Guarda a posição do canto na célula e dos vértices 
// extremos deste canto. ======================================= ok
		corner_m_pos = idsC_e[i];
		
		cell_old->getCornerNodesId(corner_m_pos, edgeNodes);
		
		if( edgeNodes[0] == vtx_t_id )
		{
			vtx_t_pos = Tetrahedron4::m_table_bC_x_vC[corner_m_pos][0];
			vtx_b_pos = Tetrahedron4::m_table_bC_x_vC[corner_m_pos][1];
		}else
		{
			vtx_t_pos = Tetrahedron4::m_table_bC_x_vC[corner_m_pos][1];
			vtx_b_pos = Tetrahedron4::m_table_bC_x_vC[corner_m_pos][0];
		}	
		
// ==== Determina os valores das posições usadas =============== ok
		getOpositVertexsPos(vtx_t_pos,vtx_b_pos, &vtx_r_pos, &vtx_l_pos);
		getOpositCornersPos(vtx_t_pos,vtx_b_pos, &corner_tr_pos,&corner_tl_pos,&corner_br_pos, &corner_bl_pos, &corner_op_pos);
		facet_r_pos=getOpositFacetPos(vtx_l_pos),
		facet_l_pos=getOpositFacetPos(vtx_r_pos),
		facet_t_pos=getOpositFacetPos(vtx_b_pos),
		facet_b_pos=getOpositFacetPos(vtx_t_pos);
		
// *********************************************************************
// Atualiza as informações dos elementos =============================== ok
//======== Atualiza as adjacencias dos pontos vtx_t ==================== ok
		vtx_t->replacesIncidCell(cell_old_id, cell_new_id, vtx_t_pos);
		// std::cout << "celulas : " << cell_old_id << " " << cell_new_id  << " vtx " << vtx_t_id << " " << vtx_b_id << " -> " << vtx_t->getIncidCell() << std::endl;
		
		if(i == 0)
		{
			vtx_m->setIncidence(cell_old_id, vtx_t_pos);
		}

		// corner_tl, corner_tr
		Corner *corner_tr = mesh->getCornerPtr(cell_old->getCornerId(corner_tr_pos)),
			   *corner_tl = mesh->getCornerPtr(cell_old->getCornerId(corner_tl_pos)),
			   *corner_op = mesh->getCornerPtr(cell_old->getCornerId(corner_op_pos));
			   
		if(corner_tl->getIncidCell() == cell_old_id )
		{ 	
		    corner_tl->setIncidence(cell_new_id,corner_tl->getPosition());
		}
		    
		if(corner_tr->getIncidCell() == cell_old_id )
		{
			corner_tr->setIncidence(cell_new_id,corner_tr->getPosition());
		}
		
		if( (corner_op->getIncidCell()==cell_old_id) && (cell_old->getIncidCell(facet_t_pos)==-1) )
		{
			if( !(cell_old->getIncidCell(facet_b_pos)==-1) )
			{			
				corner_op->setIncidence(cell_new_id, corner_op_pos);
			}
			else
			{
				if( Tetrahedron4::m_table_bC_x_vC[corner_op_pos][0] == vtx_l_pos )
				{
					corner_op->setIncidence(cell_new_id, corner_op_pos);
				}
			}			
		}
		
// Atribui as adjacencias dos cantos horizontais =========================== ok
		// corner_l, corner_r
		Corner *corner_l = mesh->getCornerPtr(horiCorners_ids[i]),
			   *corner_r = mesh->getCornerPtr(horiCorners_ids[(i+1)%nCorners]), 
			   *corner_t = mesh->getCornerPtr(vertCorner_id),
			   *corner_m = mesh->getCornerPtr(corner_m_id);
		
		tag = mesh->getFacetPtr(cell_old->getFacetId(facet_l_pos))->getTag();
		corner_l->setAllMembers(&cell_old_id, &corner_tl_pos, &tag, NULL);
				
		tag = mesh->getFacetPtr(cell_old->getFacetId(facet_r_pos))->getTag();
		corner_r->setAllMembers(&cell_old_id, &corner_tr_pos, &tag, NULL);
		
		if( corner_m->getIncidCell() == cell_old_id )
		{
		      corner_t->setIncidence(cell_new_id, corner_m_pos);
		}	
		
// ======== Atribui as adjacencias da face facet_t ====================== ok
		Facet *facet_tl = mesh->getFacetPtr( vertFaces_ids[i] ),
			  *facet_tr = mesh->getFacetPtr( vertFaces_ids[(i+1)%nCorners] ),
			  *facet_bl = mesh->getFacetPtr( cell_old->getFacetId(facet_l_pos) ),
			  *facet_br = mesh->getFacetPtr( cell_old->getFacetId(facet_r_pos) ),
			  *facet_m = mesh->getFacetPtr( horiFaces_ids[i] ),
			  *facet_t = mesh->getFacetPtr( cell_old->getFacetId(facet_t_pos));
		
		if( facet_t->getIncidCell() == cell_old_id )
		{
		     facet_t->setIncidCell(cell_new_id);
		}
		    
		int bnd;						// Componente conexa e tag da célula
		if( facet_bl->getIncidCell() == cell_old_id )
		{
			tag = facet_bl->getTag();
			bnd = facet_bl->getBoundaryComponentId();
			facet_tl->setAllMembers(&cell_new_id, &facet_l_pos, &tag, NULL, &bnd);	
		}
		
		if( facet_br->getIncidCell() == cell_old_id )
		{
			tag = facet_br->getTag();
			bnd = facet_br->getBoundaryComponentId();
			facet_tr->setAllMembers(&cell_new_id, &facet_r_pos, &tag, NULL, &bnd);	
		}
		
		tag = cell_old->getTag();
		bnd = -1;
		facet_m->setAllMembers(&cell_new_id, &facet_b_pos, &tag, NULL, &bnd);

// **********************************************************************************************
// =========================================================================================== //
//                     ATUALIZA  AS ADJACENCIAS DAS CELULAS 								   //
// =========================================================================================== //
		//std::cout << "\n l " << cell_old->getIncidCell(facet_l_pos) << " c " << cell_old_id << " r " << cell_old->getIncidCell(facet_r_pos) ;
		
		// Célula nova
		// Id dos nós
		cell_old->getNodesId(nds);
		nds[vtx_b_pos]=vtx_m_id;
		
		// Id dos cantos
		cell_old->getCornersId(cns);
		
		cns[corner_bl_pos]=horiCorners_ids[i];
		cns[corner_br_pos]=horiCorners_ids[(i+1)%nCorners];
		cns[corner_m_pos]=vertCorner_id;

		// Id das faces
		fts[facet_t_pos]=cell_old->getFacetId(facet_t_pos);
		fts[facet_b_pos]=horiFaces_ids[i];
		fts[facet_l_pos]=vertFaces_ids[i];
		fts[facet_r_pos]=vertFaces_ids[(i+1)%nCorners];		
		
		// Incidencias 
		// Celula de cima
		ics[facet_t_pos]		= cell_old->getIncidCell(facet_t_pos);
		ics_pos[facet_t_pos] 	= cell_old->getIncidCellPos(facet_t_pos);
		ics_anc[facet_t_pos] 	= cell_old->getIncidCellAnch(facet_t_pos);
		
		// Celula de baixo
		ics[facet_b_pos] = cell_old_id;
		ics_pos[facet_b_pos] = facet_t_pos;
		ics_anc[facet_b_pos] = 0; // Sera atualizada a seguir
		
		// Celula da direita
		if( cell_old->getIncidCell(facet_r_pos) == -1)
		{
			ics[facet_r_pos] = -1;
		} 
		else 
		{
			ics[facet_r_pos] = idsC_neighbor[(i+1)%nStar];
			/*
			if( cell_old->getIncidCell(facet_r_pos) == idsC[(i+1)%nCorners] )
			{
				ics[facet_r_pos] = idsC_neighbor[(i+1)%nCorners];
			}
			else
			{
				ics[facet_r_pos] = idsC_neighbor[(i+nStar-1)%nStar];
				std::cout << "caso estranho "<< edge_in_boundary << " : " << cell_old->getIncidCell(facet_r_pos) << std::endl;
			}*/
		}
		
		ics_pos[facet_r_pos] = cell_old->getIncidCellPos(facet_r_pos);
		ics_anc[facet_r_pos] = cell_old->getIncidCellAnch(facet_r_pos);
		
		// Celula da esquerda
		if( cell_old->getIncidCell(facet_l_pos) == -1)
		{
			ics[facet_l_pos] = -1;
		} 
		else 
		{
			ics[facet_l_pos] = idsC_neighbor[(i+nStar-1)%nStar];
			/*
			if( !(cell_old->getIncidCell(facet_l_pos) == idsC[(i+1)%nCorners]) )
			{
				ics[facet_l_pos] = idsC_neighbor[(i+nStar-1)%nCorners];
			}
			else
			{
				ics[facet_l_pos] = idsC_neighbor[(i+1)%nCorners];
				std::cout << "caso estranho esq "<< edge_in_boundary << " : " << cell_old->getIncidCell(facet_l_pos) << " " << idsC[(i+1)%nCorners] << std::endl;
			}*/
		}
		
		delete idsC_neighbor;
		
		ics_pos[facet_l_pos] = cell_old->getIncidCellPos(facet_l_pos);;
		ics_anc[facet_l_pos] = cell_old->getIncidCellAnch(facet_l_pos);
		
		// Indicadores
		bnd = cell_old->getConnectedComponentId();
		tag = cell_old->getTag();
  
		cell_new->setAllMembers(nds, cns, fts, ics, ics_pos, ics_anc, &bnd, &tag, NULL);	
		
		// Célula velha
		cell_old->setNodeId(vtx_t_pos, vtx_m_id);
		
		cell_old->setCornerId(corner_tl_pos, horiCorners_ids[i]);
		cell_old->setCornerId(corner_tr_pos, horiCorners_ids[(i+1)%nCorners]);
		
		cell_old->setFacetId(facet_t_pos, horiFaces_ids[i]);
		
		cell_old->setIncidCell(facet_t_pos, cell_new_id);
		cell_old->setIncidCellPos(facet_t_pos, facet_b_pos);
		
		// Calculo da ancora das faces
		int nods_m_new[3], nods_m_old[3];
		cell_old->getFacetVerticesId(facet_t_pos, nods_m_old);
		cell_new->getFacetVerticesId(facet_b_pos, nods_m_new);
				
		int anchor = calcAnchors(nods_m_old, nods_m_new, mesh);
		
		cell_old->setIncidCellAnch(facet_t_pos, anchor);
		cell_new->setIncidCellAnch(facet_b_pos, anchor);
		 
// ============= Atualiza a referencia de incidencia na celula incidente de cima ========== ok
		if( cell_new->getIncidCell(facet_t_pos) != -1)
		{
			Cell *cell_t = mesh->getCellPtr(cell_new->getIncidCell(facet_t_pos));
			cell_t->setIncidCell( cell_new->getIncidCellPos(facet_t_pos) , cell_new_id);			
		}
	}
	
	delete nds;
	delete fts;
	delete cns;
	delete ics;
	delete ics_pos;
	delete ics_anc;
	delete edgeNodes;
	delete horiCorners_ids;
	delete vertFaces_ids;
	delete horiFaces_ids;

	return vtx_m_id;
}

int MeshToolsTetr::getOpositFacetPos(int vtx_pos)
{
	int opositFacetPos[4]={3,2,1,0};
	return opositFacetPos[vtx_pos];	
}

void MeshToolsTetr::getOpositVertexsPos(int vtx_t_pos, int vtx_b_pos, int* vtx_r_pos, int* vtx_l_pos)
{
	// OpositVertexsPos[vtx_t_pos][vtx_b_pos][0]  retorna a posição do vertice vtx_r
	// OpositVertexsPos[vtx_t_pos][vtx_b_pos][1]  retorna a posição do vertice vtx_l
	FEPIC_CHECK(vtx_t_pos != vtx_b_pos, "vtx's must be different", std::runtime_error);

	int OpositVertexsPos[4][4][2] = { { {-1,-1} , {2,3} , {3,1} , {1,2} },
									  { {3,2} , {-1,-1} , {0,3} , {2,0} },
									  { {1,3} , {3,0} , {-1,-1} , {0,1} },
									  { {2,1} , {0,2} , {1,0} , {-1,-1} } };
									  
	(*vtx_r_pos) = OpositVertexsPos[vtx_t_pos][vtx_b_pos][0];
	(*vtx_l_pos) = OpositVertexsPos[vtx_t_pos][vtx_b_pos][1];
}

void MeshToolsTetr::getOpositCornersPos(int vtx_t_pos, int vtx_b_pos, int* corner_tr_pos, int* corner_tl_pos, 
																	 int* corner_br_pos, int* corner_bl_pos,
																	 int* corner_op_pos)
{
	// OpositVertexsPos[vtx_t_pos][vtx_b_pos][0]  retorna a posição do canto corner_tr
	// OpositVertexsPos[vtx_t_pos][vtx_b_pos][1]  retorna a posição do canto corner_tl
	// OpositVertexsPos[vtx_t_pos][vtx_b_pos][2]  retorna a posição do canto corner_br
	// OpositVertexsPos[vtx_t_pos][vtx_b_pos][3]  retorna a posição do canto corner_bl
	// OpositVertexsPos[vtx_t_pos][vtx_b_pos][4]  retorna a posição do canto corner_op
	FEPIC_CHECK(vtx_t_pos != vtx_b_pos, "vtx's must be different", std::runtime_error);
	
	int OpositVertexsPos[4][4][5] = { { {-1,-1,-1,-1,-1} , { 2, 3, 1, 5, 4} , { 3, 0, 4, 1, 5} , { 0, 2, 5, 4, 1} },
									  { { 5, 1, 3, 2, 4} , {-1,-1,-1,-1,-1} , { 0, 5, 2, 4, 3} , { 1, 0, 4, 3, 2} },
									  { { 1, 4, 0, 3, 5} , { 4, 2, 5, 0, 3} , {-1,-1,-1,-1,-1} , { 2, 1, 3, 5, 0} },
									  { { 4, 5, 2, 0, 1} , { 3, 4, 0, 1, 2} , { 5, 3, 1, 2, 0} , {-1,-1,-1,-1,-1} } };
									  
	(*corner_tr_pos) = OpositVertexsPos[vtx_t_pos][vtx_b_pos][0];
	(*corner_tl_pos) = OpositVertexsPos[vtx_t_pos][vtx_b_pos][1];
	(*corner_br_pos) = OpositVertexsPos[vtx_t_pos][vtx_b_pos][2];
	(*corner_bl_pos) = OpositVertexsPos[vtx_t_pos][vtx_b_pos][3];
	(*corner_op_pos) = OpositVertexsPos[vtx_t_pos][vtx_b_pos][4];
}

bool MeshToolsTetr::checkMesh(Mesh *mesh)
{
	int nCells = mesh->numCellsTotal();
	
	for(int i=0; i<nCells; i++)
	{
		int nFacetsPerCell=mesh->numFacetsPerCell();
		Cell *cell = mesh->getCellPtr(i);
		
		if( !cell->isDisabled() )
		{
			for(int j=0; j< nFacetsPerCell; j++)
			{
				if( cell->getIncidCell(j) != -1 )
				{
					// Check incidencias and position
					Cell *cell_neig = mesh->getCellPtr(cell->getIncidCell(j));
					if( cell_neig->getIncidCell( cell->getIncidCellPos(j) ) != i ) return false; 
					
					// Check anchor
					if( cell->getIncidCellAnch(j) != cell_neig->getIncidCellAnch( cell->getIncidCellPos(j) ) ) return false;
					
					// Check nodes
				    int anch = cell->getIncidCellAnch(j), 
				        nNodesPerFacet=mesh->numNodesPerFacet();
				    for( int k=0; k< nNodesPerFacet; k++)
				    {
						int pos = ((2*nNodesPerFacet-1)-k-anch)%nNodesPerFacet,
						    nk_inCell = cell->getNodeId(      Tetrahedron4::m_table_fC_x_nC[j][k] ),
						    nk_inNeig = cell_neig->getNodeId( Tetrahedron4::m_table_fC_x_nC[(int)cell->getIncidCellPos(j)][pos] );
						    
						if( nk_inCell != nk_inNeig ) return false;
					}
				    
					// Check corners
					int nCornersPerFacet=mesh->numCornersPerFacet();
					for(int k=0; k<nCornersPerFacet; k++)
					{
						int pos = ((2*nCornersPerFacet-2)-k-anch)%nCornersPerFacet,
						    ck_inCell = cell->getCornerId(      Tetrahedron4::m_table_fC_x_bC[j][k] ),
						    ck_inNeig = cell_neig->getCornerId( Tetrahedron4::m_table_fC_x_bC[(int)cell->getIncidCellPos(j)][pos] );
						    
						if( ck_inCell != ck_inNeig ) return false;						
					}			
					
					// Check facets
					if( cell->getFacetId(j) != cell_neig->getFacetId(cell->getIncidCellPos(j)) ) return false;
				
				}
				else
				{
					// Check nodes
				    int nNodesPerFacet=mesh->numNodesPerFacet();
				    for( int k=0; k< nNodesPerFacet; k++)
				    {
						int nk_inCell = cell->getNodeId(Tetrahedron4::m_table_fC_x_nC[j][k] );
						    
						if( !mesh->inBoundary( mesh->getNodePtr(nk_inCell) ) ) return false;
					}
				    
				    
					// Check corners
					int nCornersPerFacet=mesh->numCornersPerFacet();
					for(int k=0; k<nCornersPerFacet; k++)
					{
						int ck_inCell = cell->getCornerId( Tetrahedron4::m_table_fC_x_bC[j][k] );
						
						if( !mesh->inBoundary( mesh->getCornerPtr(ck_inCell)) ) return false;
					}
					
					// Check facets
					if(  !mesh->inBoundary( mesh->getFacetPtr( cell->getFacetId(j) ) ) ) return false;
					
				}
			}
			
			if( volume(cell, mesh) < 0 ) return false;
		}
		
	}
	
	int nFacets = mesh->numFacetsTotal();
	
	for(int i=0; i<nFacets; i++)
	{
		Facet *facet = mesh->getFacetPtr(i);
		
		if( !facet->isDisabled() )
		{
			Cell *cell = mesh->getCellPtr(facet->getIncidCell());
			int pos = facet->getPosition();
			
			if( cell->getFacetId( pos ) != i ) return false;
		}
	}
	
	
	int nCorners = mesh->numCornersTotal();	
	for(int i=0; i<nCorners; i++)
	{
		Corner *corner = mesh->getCornerPtr(i);
		
		if( !corner->isDisabled() )
		{
			Cell *cell = mesh->getCellPtr(corner->getIncidCell());
			int pos = corner->getPosition();
			
			if( cell->getCornerId( pos ) != i ) return false;
		}
	}
		
	int nNodes = mesh->numNodesTotal();
	for(int i=0; i<nNodes; i++)
	{
		Point *point = mesh->getNodePtr(i);
				
		if( !point->isDisabled() )
		{
			Cell *cell = mesh->getCellPtr(point->getIncidCell());
			int pos = point->getPosition();
			//std::cerr << i << " " << cell->getNodeId( pos ) << " " << point->getIncidCell() << " " << pos << std::endl;
			if( cell->getNodeId( pos ) != i ) return false;
		}
	}// */
	
	return true;
}

Real MeshToolsTetr::volume(Cell const* cell, Mesh const* mesh)
{
  Point const *a = mesh->getNodePtr(cell->getNodeId(0));
  Point const *b = mesh->getNodePtr(cell->getNodeId(1));
  Point const *c = mesh->getNodePtr(cell->getNodeId(2));
  Point const *d = mesh->getNodePtr(cell->getNodeId(3));

  Real const ax = a->getCoord(0), ay = a->getCoord(1), az = a->getCoord(2);
  Real const bx = b->getCoord(0), by = b->getCoord(1), bz = b->getCoord(2);
  Real const cx = c->getCoord(0), cy = c->getCoord(1), cz = c->getCoord(2);
  Real const dx = d->getCoord(0), dy = d->getCoord(1), dz = d->getCoord(2);
  /*
         | 1 ax ay az |
     det | 1 bx by bz | 
         | 1 cx cy cz |
         | 1 dx dy dz |

         | A  B  C |
     det | D  E  F | = A*(E*I - F*H) - B*(D*I - F*G) + C*(D*H - E*G)
         | G  H  I |
  */
  
  Real const A = bx-ax, B = by-ay, C = bz-az;
  Real const D = cx-ax, E = cy-ay, F = cz-az;
  Real const G = dx-ax, H = dy-ay, I = dz-az;

  Real const det = A*(E*I - F*H) - B*(D*I - F*G) + C*(D*H - E*G);

  return det;
}
