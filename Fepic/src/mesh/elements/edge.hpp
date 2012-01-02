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

#ifndef FEPIC_EDGE_HPP
#define FEPIC_EDGE_HPP



/**
 *  Esta classe representa segmentos de linha orientados. As Edge podem
 * ter 2 ou mais nós, um nó em cada extremo e o restante no seu interior
 * com espaçamento constante. Uma Edge com N nós tem a seginte forma: \n
 * 
 * *n0_____*n2______*n3____ ... ____*nN-1_____*n1 \n
 * 
 * Uma Edge de N+1 nós é dita Edge de ordem N, e sua orientação é definida
 * como a direção ao longo da Edge de n0 até n1.
 * 
 */ 
 
 // linear edge
template<int N> // N = order
class Edge :  public _CellCore<Edge<N> >
{
public:
  friend class _CellCore<Edge<N> >;

  typedef Polytope<1> PolytopeT;
  typedef Point       Derived;
  
  enum { dim=1,
         n_vertices=2,
         n_nodes=N+1,
         n_facets=2,
         n_corners=0,
         n_vertices_per_facet=1,
         n_vertices_per_corner=0,
         n_nodes_per_facet=1,
         n_nodes_per_corner=0};

  static const bool has_edge_nodes = N==1 ? false : true;
  static const bool has_face_nodes = false;
  static const bool has_volume_nodes = false;
  static const ECellType fep_tag = N==1 ? EDGE2     : (N==2 ? EDGE3     : UNDEFINED_CELLT);
  static const EMshTag msh_tag = N==1 ? MSH_LIN_2 : (N==2 ? MSH_LIN_3 : MSH_UNDEFINED_ELEM);

  /** Verifica se esta Edge está alinha com outra Edge. O alinhamento é verificado
  * comparando-se os nós dos extremos das Edges.
  *  @param e edge na qual se vai verificar o alinhamento.
  *  @return 1 se a edge é paralela, -1 se a edge é anti-paralela, e 0 se não é paralela.
  */
  int isAligned (Edge const& e) const
  {
    if ( (_nodes[0] == e._nodes[0]) && (_nodes[1] == e._nodes[1]) )
      return 1;
    else if ( (_nodes[0] == e._nodes[1]) && (_nodes[1] == e._nodes[0]) )
      return -1;
    else
      return 0;
  }
  
  /** Verifica se esta Edge está alinha com dois nós dados. O alinhamento é verificado
  * comparando-se os nós dos extremos desta Edge com os dois nós.
  *  @param nodes os vetor com os dois nós nos quais se vai verificar o alinhamento.
  *  @return 1 se a edge é paralela, -1 se a edge é anti-paralela, e 0 se não é paralela.
  */
  int isAligned (int* nodes) const
  {
    if ((_nodes[1] == nodes[1]) && (_nodes[0] == nodes[0]) )
      return 1;
    else if ( (_nodes[0] == nodes[1]) && (_nodes[1] == nodes[0]) )
      return -1;
    else
      return 0;
  }       

  
  /** NÃO IMPLEMENTADO
  *  Retorna o comprimento desta Edge.
  */ 
  double lenght() const; // IMPLEMENTAR E LEMBRAR QUE PODER SER LINHA CURVA
  
  virtual ~Edge() {};

private:
  char _icells_pos[2];
  int _icells[2];  // incident cells id
  int _nodes[N+1]; // nodes id
  
  /* não utilizados */
  int _icells_anchors[0];
  int _facets[0];
  int _corners[0];
  
  
public:
  static const int table_fC_x_vC[2][1];
  static const int table_fC_x_nC[2][1];
  static const int table_bC_x_nC[0][0];
  static const int table_bC_x_vC[0][0];
  static const int table_fC_x_bC[0][0];  
  static const int table_bC_x_fC[0][0];
  static const int table_vC_x_fC[2][1];
};





#endif

