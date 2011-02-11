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
template<class _Traits>
class Edge : public _Labelable, public _CellCore<_Traits>
{
public:

  typedef          Polytope<1>     PolytopeT;
  typedef typename _Traits::PointT PointT;
  typedef typename _Traits::HalfT  HalfT;      
  typedef typename _Traits::MeshT  MeshT;  
  
  friend class _CellCore<_Traits>;
  
  enum { dim=1,
         n_vertices=2,
         n_borders=2,
         n_vertices_per_border=1};

  template<class... LabeableArgs>
  Edge(vectori const& nodes, int order, LabeableArgs... args) :
                          _Labelable(args...), _nodes(nodes), _order(order)
  {
    FEPIC_CHECK(nodes.size()==numNodes<Polytope<1>>(order), "", std::invalid_argument);
  }
  
  template<class... LabeableArgs>
  Edge(vectori && nodes, int order, LabeableArgs... args) :
                      _Labelable(args...), _nodes(nodes), _order(order)
  {
    FEPIC_CHECK(nodes.size()==numNodes<Polytope<1>>(order), "", std::invalid_argument);
  }  
  
  Edge() : _nodes({0,0}), _order(1) {};
  Edge(Edge const&) = default;
  ~Edge() = default;

    /** Atualiza esta Edge para a ordem <em>order</em>.
  *  @warning Toda vez que a ordem da malha for alterada, essa função DEVE SER CHAMADA.
  */ 
  void setOrder(int order)
  {
    _nodes.resize(order+1);
    this->_order = order;
  }
  
  /** Retorna o número de nós desta Edge.
  */ 
  int getNumNodes() const
  {
    return _nodes.size();
  }
  
  /** Retorna o ith-ésimo nó desta Edge.
  */ 
  int getNodeIdx(int const ith) const
  {
    return _nodes[ith];
  }
  
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
  int isAligned (vectori const& nodes) const
  {
    if ((_nodes[1] == nodes[1]) && (_nodes[0] == nodes[0]) )
      return 1;
    else if ( (_nodes[0] == nodes[1]) && (_nodes[1] == nodes[0]) )
      return -1;
    else
      return 0;
  }       
  
  /** Atribui o i-ésimo nó da aresta como o n-ésimo nó da malha.
  */ 
  void setNode(int const ith, int const nth)
  {
    _nodes[ith] = nth;
  }
  
  /** NÃO IMPLEMENTADO
  *  Retorna o comprimento desta Edge.
  */ 
  double lenght() const; // IMPLEMENTAR E LEMBRAR QUE PODER SER LINHA CURVA
  
  /** Imprime os nós da aresta. 
  *  @param o o stream aonde se vai escrever.
  *  @param order a ordem da aresta.
  */ 
  void printSelfVtk(std::ostream &o, int order) const
  {
  

          
  }
  
  
  /** Imprime esta Edge no formato State
  */ 
  void printSelfState(std::ostream &o, int order)
  {
    o << getNodeIdx(0);
    for (int i = 0; i < order+1; ++i)
      o << " " << getNodeIdx(i);
  }
  
  /** Retorna a tag de uma Edge no formato Msh.
  *  @param order a ordem da Edge.
  */  
  static int getMshTag(int order)
  {
    /* poderia fazer meta programção aqui, mas acho complicação desnecessária */
    switch (order)
    {
      case 1: return MSH_LIN_2;
      case 2: return MSH_LIN_3;
      case 3: return MSH_LIN_4;
      case 4: return MSH_LIN_5;
      case 5: return MSH_LIN_6;
      default:
      {
        std::cout << "edge order not supported yet" << std::endl;
        throw;
      }
    }
  }

protected:  
  unsigned char _order;
  vectori       _nodes;
};






#endif
