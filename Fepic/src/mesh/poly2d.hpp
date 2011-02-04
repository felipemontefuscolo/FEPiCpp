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
template<class _Traits>
class _Poly2d : public _CellCore<_Traits>
{
#if !defined(THIS) && !defined(CONST_THIS)
  #define THIS static_cast<CellT*>(this)
  #define CONST_THIS static_cast<const CellT*>(this)
#endif  

public:
  
  typedef typename _Traits::CellT  CellT;
  typedef typename _Traits::HalfT  HalfT;
  typedef typename _Traits::HalfLT HalfLT;
  typedef typename _Traits::MeshT  MeshT;
  //typedef typename ElementProperties<CellT, _Traits>::FaceT FaceT;

protected:      
  _Poly2d(_Poly2d const&) {};
  _Poly2d() {};

public:  
  /** Retorna se os vertices passados formam uma aresta do polígono.
  * @param[in] vertices um vetor com exatamente 2 vertices.
  * @param[out] ith o índice local da aresta que os pontos formam.
  * @return true se e somente se forma uma aresta.
  * @warning a ordem dos nós É importante.
  */
  bool isAnEdge(vectorui const& vertices, int &ith) const
  {
    FEPIC_ASSERT(vertices.size()==2, "");
       
    vectorui vtx(2);
                
    for (int i = 0; i != CellT::n_borders; ++i)
    {
      vtx[0] = CONST_THIS->_nodes[i];
      vtx[1] = CONST_THIS->_nodes[(i+1)%CellT::n_borders];
      
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
  bool isAnParallelEdge(vectorui const& vertices, int &ith) const
  {
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

  /** Imprime este polígono no formato State (seus nós).
  * @param o a stream onde se vai imprimir, e.g., std::cout.
  */ 
  void printSelfState(std::ostream &o) const
  {
    o << this->getNodeIdx(0);
    for (uint i = 1; i < CONST_THIS->_nodes.size(); ++i)
      o << " " << this->getNodeIdx(i);
  }

#undef THIS
#undef CONST_THIS 
};






#endif
