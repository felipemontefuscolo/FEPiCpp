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
 * com espaçamento constante. Uma Edge com N nós tem a seguinte forma: \n
 * 
 * *n0_____*n2______*n3____ ... ____*nN-1_____*n1 \n
 * 
 * Uma Edge de N+1 nós é dita Edge de ordem N, e sua orientação é definida
 * como a direção ao longo da Edge de n0 até n1.
 * 
 */ 
template<class Traits>
class Edge : public _Labelable {
public:
  static const int Dim = 1;

  typedef Simplex<1>  ElmClass;
    
  typedef typename Traits::CellT                         CellT;
  typedef typename VolumeDef<CellT::Dim, CellT>::VolumeT VolumeT;
  typedef typename FaceDef<CellT::Dim, CellT>::FaceT     FaceT; 
    
  typedef typename Traits::PointT  PointT;       
  typedef typename Traits::MeshT   MeshT;        
  typedef typename Traits::PointT  BorderT;  
  typedef UndefElement             BndBorderT;

  /** Construtor.
  * @param nodes vetor com os nós que compõe a edge.
  * @note devem ser passados pelo menos dois nós.
  */ 
  Edge(Fepic::vectorui const& nodes, int label=0) : _Labelable(label), _node(nodes)
  {
#ifdef FEPIC_DEBUG_ON
    if (_node.size() < 2)
    {
      std::cout << "erro: Edge constructor: devem ser passados pelo menos dois nós\n";
      throw;
    }
#endif
  }
    
  /** Construtor.
  */ 
  Edge() : _Labelable(), _node({0,0})
  {
  }
  
  /** Retorna o número de nós de uma Edge dada uma ordem order.
  */ 
  static int getNumCellNodes(int const order)
  {
    return order + 1;
  }
  
  /** Retorna o número de nós desta Edge.
  */ 
  int getNumNodes() const
  {
    return _node.size();
  }
  
  /** Retorna o ith-ésimo nó desta Edge.
  */ 
  uint getNodeIdx(int const ith) const
  {
    return _node[ith];
  }
  
  /** Verifica se esta Edge está alinha com outra Edge. O alinhamento é verificado
  * comparando-se os nós dos extremos das Edges.
  *  @param e edge na qual se vai verificar o alinhamento.
  *  @return 1 se a edge é paralela, -1 se a edge é anti-paralela, e 0 se não é paralela.
  */
  int isAligned (Edge const& e) const
  {
    if ( (_node[0] == e._node[0]) && (_node[1] == e._node[1]) )
      return 1;
    else if ( (_node[0] == e._node[1]) && (_node[1] == e._node[0]) )
      return -1;
    else
      return 0;
  }
  
  /** Verifica se esta Edge está alinha com dois nós dados. O alinhamento é verificado
  * comparando-se os nós dos extremos desta Edge com os dois nós.
  *  @param nodes os vetor com os dois nós nos quais se vai verificar o alinhamento.
  *  @return 1 se a edge é paralela, -1 se a edge é anti-paralela, e 0 se não é paralela.
  */
  int isAligned (Fepic::vectorui const& nodes) const
  {
    if ((_node[1] == nodes[1]) && (_node[0] == nodes[0]) )
      return 1;
    else if ( (_node[0] == nodes[1]) && (_node[1] == nodes[0]) )
      return -1;
    else
      return 0;
  }       
  
  /** Atribui o i-ésimo nó da aresta como o n-ésimo nó da malha.
  */ 
  void setNode(int const ith, uint const nth)
  {
    _node[ith] = nth;
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
  
    if (order<=1)
      o << "2 " << _node[0] << " " << _node[1];
    else
    {
      o << "2 " << _node[0] << " " << _node[2];
      for (int i = 0; i < order-2; i++)
      {   
        o << std::endl;
        o << "2 " << _node[2+i] << " " << _node[3+i];
      }
      o << std::endl;
      o << "2 " << _node[order] << " " << _node[1];
    }
          
  }
  
  /** INCOMPLETO: APENAS LINEAR
  * Retorna a tag desta Edge definida no formato Vtk.
  */ 
  static int getCellTypeVtk()
  {
    /* TRABALHO: apenas linear
    */
    return 3; // Vtk_LINE (3)
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
      case 1: return Msh_LIN_2;
      case 2: return Msh_LIN_3;
      case 3: return Msh_LIN_4;
      case 4: return Msh_LIN_5;
      case 5: return Msh_LIN_6;
      default:
      {
        std::cout << "edge order not supported yet" << std::endl;
        throw;
      }
    }
  }
  
  /** Atualiza esta Edge para a ordem <em>order</em>.
  *  @warning Toda vez que a ordem da malha for alterada, essa função DEVE SER CHAMADA.
  */ 
  void setOrder(int order)
  {
    _node.resize(order+1);
  }
  
  static const int n_borders = 2;
  static const int n_vertices = 2;
  
  /** Destrutor.
  */ 
  ~Edge() {}
  
protected:
  Fepic::vectori _node;
};






#endif
