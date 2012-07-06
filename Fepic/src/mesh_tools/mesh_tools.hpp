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


#ifndef FEPIC_MESH_TOOLS_HPP
#define FEPIC_MESH_TOOLS_HPP

#include "../util/common.hpp"
#include "../mesh/mesh.hpp"
#include <utility>
#include <vector>

class Mesh;
class Cell;

class MeshTools
{
public:

  /** Safely removes a cell
   *  @param cell the cell that will be removed.
   *  @param mesh mesh context.
   *  @return nothing.
   */ 
  static void removeCell(Cell * cell, Mesh *mesh);
  
  
};

/// specific for triangular meshes
class MeshToolsTri
{
public:

  /** flips Tri3 only TODO
   *  @param cell the cell that will be removed.
   *  @param fid local facet's id that will be flipped.
   *  @param mesh mesh context.
   *  @return true if an error occurred, false otherwise.
   */ 
  static bool flip(Cell * cell, int fid, Mesh *mesh);

  /** Tri3 only
   *  slice a convex part of the mesh with an isocontour of a 
   *  scalar function F(x,y).
   *  @param x0 a point such that F(x0)=0;
   *  @param c0 a initial cell (need not contain the point x0,
   *  but should be close enough to avoid convex issues).
   *  @param F a functor that has an operator () (double *x), where x is
   *  pointer to the coordinates.
   *  @return returns true if it was successful.
   *  @note the mesh has to be convex, otherwise the result is not guaranteed.
   * 
   */
  template<class Functor>
  static bool cutConvexPart(Real *x0, Cell *c0, Functor const& F, Mesh * mesh)
  {
    std::vector<int> cells_id;
    std::vector<Real> slices;
    std::vector<int> mid_nodes(cells_id.size()*6, int(-1));
    
    MeshToolsTri::createPath(x0, c0, F, cells_id, slices,  mesh);
    
    Point pt;
    Cell  ce;
    Facet fa;
    int new_cells[3];
    
    for (int cth = 0; cth < (int)cells_id.size(); ++cth)
    {
      c0 = mesh->getCell(cells_id[cth]);
      
      int I;
      for (int i = 0; i < 3; ++i)
      {
        if (slices[3*cth + i] < 0)
        {
          I = i;
          break;
        }
      }
      int J = (I+1)%3;
      int K = (I+2)%3;

      // cria as 3 células
      ce.setTag(c0->getTag());
      for (int i = 0; i < 3; ++i)
      {
        ce.setIncidCell(i, -1);
        ce.setIncidCellPos(i, -1);
      }
      new_cells[0] = mesh->pushCell(&ce);
      new_cells[1] = mesh->pushCell(&ce);
      new_cells[2] = mesh->pushCell(&ce);

      // para cada lado que é cortado
      for (int kk=0; kk<2; ++kk)
      {
        int f = (I+kk+1)%3;
        int f_nds[2];
        c0->getFacetNodesId(f,f_nds);
        Real xa[2] = {mesh->getNode(f_nds[0])->getCoord(0), mesh->getNode(f_nds[0])->getCoord(1)};
        Real xb[2] = {mesh->getNode(f_nds[1])->getCoord(0), mesh->getNode(f_nds[1])->getCoord(1)};
        Real xmid[2] = {(xb[0]-xa[0])*slices[3*cth+f]+xa[0],(xb[1]-xa[1])*slices[3*cth+f]+xa[1]};
        //Real Fa = F(xa);
        //Real Fb = F(xb);
        
        if (mid_nodes[6*cth + 2*f] < 0) // não foi construído nessa f
        {
          //cria os vértices
          pt.setCoord(xmid);
          pt.setTag(c0->getTag());
          mid_nodes[6*cth + 2*f + 0] = mesh->pushPoint(&pt);  // F(x) < 0
          mid_nodes[6*cth + 2*f + 1] = mesh->pushPoint(&pt);  // F(x) > 0
          
          int oth = c0->getIncidCell(f);
          std::vector<int>::iterator itl = find(cells_id.begin(), cells_id.end(), oth);
          oth = std::distance(itl);
          int oth_f = c0->getIncidCellPos(f);
          if (oth>=0)
          {
            mid_nodes[6*oth + 2*oth_f + 0] = mid_nodes[6*cth + 2*f + 0];  // F(x) < 0
            mid_nodes[6*oth + 2*oth_f + 1] = mid_nodes[6*cth + 2*f + 1];  // F(x) > 0
          }
          
          mesh->getNode(mid_nodes[6*cth + 2*f + 0])->setIncidCell(new_cells[kk+1]);
          mesh->getNode(mid_nodes[6*cth + 2*f + 1])->setIncidCell(new_cells[(kk+2)%3]);
          mesh->getNode(mid_nodes[6*cth + 2*f + 0])->setPosition(new_cells[(kk+2)%3]);
          mesh->getNode(mid_nodes[6*cth + 2*f + 1])->setPosition(new_cells[1]);
          
          int new_efacets[2]; // externo
          fa.setTag(c0->getTag());
          fa.setIncidCell(new_cells[kk+1]);
          fa.setPosition(kk+1);
          new_efacets[0] = mesh->pushFacet(&fa); 
          fa.setIncidCell(new_cells[(kk+2)%3]);
          fa.setPosition(1);
          new_efacets[1] = mesh->pushFacet(&fa); 
          
          mesh->getCell(new_cells[kk+1])->setFacetId(kk+1,new_efacets[0]); 
          mesh->getCell(new_cells[(kk+2)%3])->setFacetId(1,new_efacets[1]); 
          
          mesh->disableFacet(mesh->getFacet(c0->getFacetId(f)));
        }
        else
        {
          for (int bingo = 0; bingo<2; ++bingo)
          {
            Cell *some_cell = mesh->getCell(mesh->getNode(mid_nodes[6*cth + 2*f + bingo])->getIncidCell());
            
            // se some_cell tem o no k, então ele está na face do sem tracejado.
            int fth_node = c0->getNodeId(f);
            for (int pos = 0; pos<3; ++pos)
            {
              if (some_cell->getNodeId(pos) == fth_node)
              {
                mesh->getCell(new_cells[kk+1])->setIncidCell(kk+1, mesh->getCellId(some_cell));
                mesh->getCell(new_cells[kk+1])->setIncidCellPos(kk+1, (pos+2)%3);
                mesh->getCell(new_cells[kk+1])->setFacetId(kk+1, some_cell->getFacetId((pos+2)%3));
                some_cell->setIncidCell((pos+2)%3, new_cells[kk+1]);
                some_cell->setIncidCellPos((pos+2)%3, kk+1);
                break;
              }
            }
            for (int pos = 0; pos<3; ++pos)
            {
              if (some_cell->getNodeId(pos) == (fth_node+1)%3)
              {
                mesh->getCell(new_cells[(kk+2)%3])->setIncidCell(1, mesh->getCellId(some_cell));
                mesh->getCell(new_cells[(kk+2)%3])->setIncidCellPos(1, pos);
                mesh->getCell(new_cells[(kk+2)%3])->setFacetId(1, some_cell->getFacetId(pos));
                some_cell->setIncidCell(pos, new_cells[(kk+2)%3]);
                some_cell->setIncidCellPos(pos, 1);
                break;
              }
            }
          }
        }
        
      } // end kk (cada lado cortado)
      
      // crias as arestas interiores
      int new_ifacets[3];
      //   1
      fa.setTag(c0->getTag());
      fa.setIncidCell(new_cells[2]);
      fa.setPosition(0);
      new_ifacets[0] = mesh->pushFacet(&fa);
      //   2
      fa.setIncidCell(new_cells[1]);
      fa.setPosition(2);
      new_ifacets[1] = mesh->pushFacet(&fa);
      //   3
      fa.setIncidCell(new_cells[0]);
      fa.setPosition(0);
      new_ifacets[2] = mesh->pushFacet(&fa);
      
      
      // seta as adjacências restantes
      mesh->getCell(new_cells[0])->setIncidCell(0, new_cells[1]);
      mesh->getCell(new_cells[0])->setIncidCellPos(0, 0);
      mesh->getCell(new_cells[0])->setIncidCell(2, c0->getIncidCell(I));
      mesh->getCell(new_cells[0])->setIncidCellPos(2, c0->getIncidCellPos(I));
      mesh->getCell(new_cells[0])->setNode(0, c0->getNodeId(J));
      mesh->getCell(new_cells[0])->setNode(1, mid_nodes[6*cth + 2*K + 1]);
      mesh->getCell(new_cells[0])->setNode(2, c0->getNodeId(I));
      mesh->getCell(new_cells[0])->setFacetId(0, new_ifacets[2]);
      mesh->getCell(new_cells[0])->setFacetId(2, c0->getFacetId(I));
      
      mesh->getCell(new_cells[1])->setIncidCell(0, new_cells[0]);
      mesh->getCell(new_cells[1])->setIncidCellPos(0, 0);
      mesh->getCell(new_cells[1])->setIncidCell(2, -1);
      mesh->getCell(new_cells[1])->setIncidCellPos(2, -1);
      mesh->getCell(new_cells[1])->setNode(0, mid_nodes[6*cth + 2*K + 1]);
      mesh->getCell(new_cells[1])->setNode(1, c0->getNodeId(J));
      mesh->getCell(new_cells[1])->setNode(2, mid_nodes[6*cth + 2*J + 0]);
      mesh->getCell(new_cells[1])->setFacetId(0, new_ifacets[2]);
      mesh->getCell(new_cells[1])->setFacetId(2, -1);
      
      mesh->getCell(new_cells[2])->setIncidCell(0, -1);
      mesh->getCell(new_cells[2])->setIncidCellPos(0, -1);
      mesh->getCell(new_cells[2])->setNode(0, mid_nodes[6*cth + 2*K + 0]);
      mesh->getCell(new_cells[2])->setNode(1, mid_nodes[6*cth + 2*J + 0]);
      mesh->getCell(new_cells[2])->setNode(2, c0->getNodeId(K));
      mesh->getCell(new_cells[2])->setFacetId(0, new_ifacets[2]);
      mesh->getCell(new_cells[2])->setFacetId(2, -1);
      
      mesh->getNode(c0->getNodeId(I))->setIncidCell(new_cells[0]);
      mesh->getNode(c0->getNodeId(I))->setPosition(2);
      mesh->getNode(c0->getNodeId(J))->setIncidCell(new_cells[1]);
      mesh->getNode(c0->getNodeId(J))->setPosition(1);
      mesh->getNode(c0->getNodeId(K))->setIncidCell(new_cells[2]);
      mesh->getNode(c0->getNodeId(K))->setPosition(2);

      mesh->disableCell(c0);
      
    } // for celula
    
    
  }
  
  /** Creates a path represented by elements slices.
   *  @param x0 a point such that F(x0)=0;
   *  @param c0 a initial cell (need not contain the point x0,
   *  but should be close enough to avoid convex issues).
   *  @param F a functor that has an operator () (double *x), where x is
   *  @param[out] cells_id vector with cells that are sliced​​.
   *  @param[out] slices vector with the slices. slices(3*c + j) is a value
   *  between 0 and 1 which indicates where the jth facet of the cth cell
   *  were sliced. If it is negative so that facet is not sliced. For example,
   *  [0.1, 0.3, 0.5] means that facet 0, facet 1 and facet 2 are sliced at 0.1, 0.3
   *  and 0.5.
   *  @param mesh mesh context.
   *  pointer to the coordinates.
   *  @note o algoritmo garante que não vai er nenhuma célula singular.
   *  @warning pode perturbar alguns pontos! ... só serve para TRIANGULO!!
   *  @warning usa o flag visited.
   */ 
  template<class Functor>
  static bool createPath(Real const* x0, Cell * c_ini, Functor const& F, std::vector<int> &cells_id, std::vector<Real> &slices, Mesh * mesh)
  {
    FEPIC_ASSERT(x0!=NULL && c_ini!=NULL && mesh!=NULL,"NULL pointers",std::runtime_error);
    
    int f_nodes[2];
    
    cells_id.clear();
    slices.clear();
    cells_id.reserve(32);
    slices.reserve(128);
    
    std::pair<bool, Cell *> S0 = searchConvexPoint(x0,c_ini,mesh);
    
    if (!S0.first)
      return false;
    
    cells_id.push_back(mesh->getCellId(S0.second));
    int cth_cell = 0;
    
    Cell *c0;
    while (cth_cell < cells_id.size())
    {
      c0 = mesh->getCell(cells_id[cth_cell]);
      
      const Real TOL = 1.e-15;
      // perturba os pontos se eles estiverem no level-set
      for (int i=0; i<3; ++i)
      {
        int const node_id = c0->getNodeId(i);
        Real const* xa = mesh->getNode(node_id)->getCoord();
        if (fabs(F(xa))<TOL)
        {
          
          Real const xp[] = {xa[0]+1.e-7, xa[1]+1.e-7};
          mesh->getNode(node_id)->setCoord(xp);
          
          // de novo, só que pra outro lado
          if (fabs(F(xa))<TOL)
          {
            Real const xp[] = {xa[0]-2.e-7, xa[1]};
            mesh->getNode(node_id)->setCoord(xp);
          };
        };
      }
      
      // assume que apenas duas facets são cortadas
      slices.push_back(-1.0);
      slices.push_back(-1.0);
      slices.push_back(-1.0);
      int f_free=-1;
      for (int f = 0; f < 3; ++f)
      {
        c0->getFacetNodesId(f, f_nodes);
        Real const* xa = mesh->getNode(f_nodes[0])->getCoord();
        Real const* xb = mesh->getNode(f_nodes[1])->getCoord();
        
        Real const Fa = F(xa);
        Real const Fb = F(xb);
        // se tem sinais opostos, cruza
        if ( ((Fa<0) && (Fb>0)) || ((Fa>0) && (Fb<0)))
        {
          slices[cth_cell*3 + f] = Fa/(Fa-Fb);
        }
        else
        {
          if (f_free>=0)
          {
            std::cout << "n pode duas faces sem cruzar" << std::endl;
            std::cout << "Fa =" << Fa << " Fb = " << Fb << std::endl;
            std::cout << "At nodes " << f_nodes[0]<<" "<<f_nodes[1] << " " << mesh->getNode(52)->getCoord(0) << " " <<mesh->getNode(52)->getCoord(1) << std::endl;
            std::cout << "f = " << f << " cell = "<<mesh->getCellId(c0)<< std::endl;
            throw; // n pode duas faces sem cruzar
          }
          f_free = f;
        }
        
      }

      // coloca os vizinhos na pilha
      int const c_a = c0->getIncidCell((f_free+1)%3);
      int const c_b = c0->getIncidCell((f_free+2)%3);
      if (c_a>=0 && !mesh->getCell(c_a)->visited())
        cells_id.push_back(c_a);
      if (c_b>=0 && !mesh->getCell(c_b)->visited())
        cells_id.push_back(c_b);
        
      c0->visited(true);
      ++cth_cell;
      
    }
    
    // remove os flags
    for (std::vector<int>::iterator cit = cells_id.begin(); cit != cells_id.end(); ++cit)
    {
      mesh->getCell(*cit)->visited(false);
    }
    
    
  }
  
  
  /** Search a point int the mesh.
   * @param x the point coordinates.
   * @param c0 a first cell to start the search.
   * @param mesh the mesh.
   * @return a "pair" type, where pair.first is a bool that is true if
   * the point was found, and pair.second is a pointer to the cell where
   * the point is. If the point was not found, pair.second returns closest
   * cell.
   * @note the mesh has to be convex, otherwise the result is not guaranteed.
   */ 
  static std::pair<bool, Cell *> searchConvexPoint(Real const* x, Cell const* c0, Mesh const* mesh);
  
  
  
  
  
};




#endif // FEPIC_IMESH_IMPL_HPP
