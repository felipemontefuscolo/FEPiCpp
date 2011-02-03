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

#ifndef FEPIC_POLY2D_HPP
#define FEPIC_POLY2D_HPP



/** A classe _Poly2d representa os polígonos (triângulo, quadrângulos, ...), entidades
 * da malha de dimensão 2. Dependendo da ordem, os polígonos dessa classe podem ter arestas curvas.
 * 
 * @note Esta é uma classe de conceito abstrato (apesar se não ser abstrata sob o ponto
 * de vista da linguagem c++ pois não tem funções virtuais puras), logo não se deve instanciá-la
 * a não ser que se saiba exatamente o que se está fazendo.
 */ 
template<class Traits>
class _Poly2d : public _Labelable
{
public:
  static const int Dim = 2;
  
  /* Se Polygon está sendo instanciado, então CellT é com certeza uma FaceT*/
  typedef typename Traits::CellT    FaceT;
  typedef typename Traits::MeshT        MeshT;
  
  static const int n_borders  = ElementProperties<FaceT, Traits>::n_borders;
  static const int n_vertices = ElementProperties<FaceT, Traits>::n_vertices;
  
  typedef Eigen::Matrix<double, Traits::spacedim, 1> VecT;
  
  /** Construtor.
  */ 
  _Poly2d() : _Labelable()
  {
    _node.resize(n_vertices);
  }
  
  /** Construtor.
  *      @param nodes vetor com os nós que compõe o _Poly2d.
  */
  _Poly2d(Fepic::vectorui const& nodes, int label=0) : _Labelable(label), _node(nodes)
  {
#ifdef FEPIC_DEBUG_ON
    if (_node.size() < static_cast<uint>(n_vertices))
    {
      std::cout << "error:_Poly2d constructor: número insuficiente de nós.\n";
      throw;
    }
#endif
  }

  /** Retorna o número de nós deste polígono.
  */ 
  int getNumNodes() const
  {
    return _node.size();
  }

  /** Retorna o índice do i-ésimo nó do polígono.
  */ 
  uint getNodeIdx(int ith) const
  {   
#ifdef FEPIC_DEBUG_ON
    return _node.at(ith);
#else
    return _node[ith];
#endif
  }

  /** Defini o i-ésimo nó do polígono como nodeid
  */ 
  void setNode(int ith, uint nodeid)
  {
#ifdef FEPIC_DEBUG_ON
    _node.at(ith) = nodeid;
#else
    _node[ith] = nodeid;
#endif
  }
  
  /** Retorna um ponteiro para a i-ésima HalfEdge.
  */ 
  HalfEdge<Traits>* getHalf(int ith)
  {
#ifdef FEPIC_DEBUG_ON
    if (ith >= n_borders)
    {
      std::cout << "getHalf: out of range" << std::endl;
      throw;
    }
#endif
    return &_he[ith];
  }

  /** Retorna se os vertices passados formam uma aresta do polígono.
  * @param[in] vertices um vetor com exatamente 2 vertices.
  * @param[out] ith o índice local da aresta que os pontos formam.
  * @return true se e somente se forma uma aresta.
  * @warning a ordem dos nós É importante.
  */
  bool isAnEdge(Fepic::vectorui const& vertices, int &ith) const
  {
          
#ifdef FEPIC_DEBUG_ON
    if (vertices.size() != 2)
    {
      std::cout << "isAnEdge: first argument has wrong dimension.\n" << std::endl;
      throw;
    }
#endif        
    Fepic::vectorui vtx(2);
                
    for (int i = 0; i != n_borders; ++i)
    {
      vtx[0] = _node[i];
      vtx[1] = _node[(i+1)%n_borders];
      
      if (vtx == vertices)
      {
        ith = i;
        return true;
      }
    }
    return false;
  }
        
  /** Retorna se os vertices passados formam uma aresta do polígono.
  * @param[in] vertices um vetor com exatamente 2 vertices.
  * @param[out] ith o índice local da aresta que os pontos formam.
  * @return true se e somente se forma uma aresta.
  * @warning a ordem dos nós NÃO é importante.
  */
  bool isAnParallelEdge(Fepic::vectorui&& vertices, int &ith) const
  {
#ifdef FEPIC_DEBUG_ON
    if (vertices.size() != 2)
    {
      std::cout << "isAnEdge: first argument has wrong dimension.\n" << std::endl;
      throw;
    }
#endif        
    if (this->isAnEdge(vertices, ith))
      return true;
                        
    std::swap(vertices[0], vertices[1]);
    if (this->isAnEdge(vertices, ith))
    {
      std::swap(vertices[0], vertices[1]);
      return true;
    }
    std::swap(vertices[0], vertices[1]);
    return false;
  }

  /** OBSOLETO: refazer\n
  *  Procura por uma edge em comum com outro polígono: se encontra retorna true, caso contrário
  * retorna false.
  */ 
  bool getCommonBorder(FaceT *other, int &thisborder, int &otherborder) const
  {
#ifdef FEPIC_DEBUG_ON
    if (other==NULL)
    {
      std::cout << "getCommonBorder: null FaceT.\n" << std::endl;
      throw;
    }
#endif
    Fepic::vectorui this_edge_vtx(2);
  
    for (int s = 0; s < n_borders; ++s)
    {
      this_edge_vtx = this->getBorderVertices(s);
      reverse(this_edge_vtx.begin(), this_edge_vtx.end());
      if (other->isAnEdge(this_edge_vtx, otherborder))
      {
        thisborder = s;
        return true;
      }
    }
  
    return false;
  }
    
  /** Atribui em cada nó deste polígono uma HalfEdge compatível.
  * @param mesh a malha em que este polígono está contido.
  */ 
  void propagateHalf(MeshT & mesh) const
  {
#ifdef FEPIC_DEBUG_ON
    if (mesh.edges_local_nodes.empty())
    {
      std::cout << "edges_local_nodes must be initializated." << std::endl;
      throw;
    }
#endif
    for (uint s = 0; s < static_cast<uint>(n_borders); ++s) // loop  nas arestas
    {
      for (uint i = 0, tam=mesh.edges_local_nodes[0].size(); i < tam; ++i)
      {
        mesh.getNode(_node[mesh.edges_local_nodes[s][i]])->setHalf(_he[s]);
      }
    }
  }

  /** Atribui o rótulo deste polígono a seus nós.
  * @param mesh a malha em que este polígono está contido.
  * @param force quando true indica atribuição incondicional, quando false,
  * a atribuição é feita somente se cada nó tem label=0;
  * @warning CONFERIR!
  */ 
  void propagateTag(MeshT & mesh, bool force=false) const
  {
    if (force)
      for (int i = 0; i < _node.size(); i++)
        mesh.getNode(_node[i])->setTag(this->getTag());
    else
      for (int i = 0; i < _node.size(); i++)
        if (mesh.getNode(_node[i])->getTag() == 0)
          mesh.getNode(_node[i])->setTag(this->getTag());
  }
  
  /** Retorna um vetor com os nós da i-ésima edge.
  */  
  Fepic::vectorui getBorderNodes(int ith, MeshT& mesh) const
  {
    Fepic::vectorui nodes;
    
    for (uint i = 0; i < mesh.edges_local_nodes[ith].size(); i++)
      nodes.push_back(_node[ mesh.edges_local_nodes[ith][i] ]);
    return nodes;
  }
  
  /** Retorna um vetor com os vértices da i-ésima edge
  */ 
  Fepic::vectorui getBorderVertices(int ith) const
  {
#ifdef FEPIC_DEBUG_ON
    if (ith >= n_borders)
    {
      std::cout << "getCommonBorder: null FaceT.\n" << std::endl;
      throw;
    }
#endif
    return {_node[ith], _node[(ith+1)%FaceT::n_borders]};
  }
    
  /** Imprime este polígono no formato State (seus nós).
  * @param o a stream onde se vai imprimir, e.g., std::cout.
  */ 
  void printSelfState(std::ostream &o) const
  {
    o << getNodeIdx(0);
    for (uint i = 1; i < _node.size(); ++i)
      o << " " << getNodeIdx(i);
  }

protected:
  Fepic::vectorui     _node;
  HalfEdge<Traits>   _he[n_borders];
};






#endif
