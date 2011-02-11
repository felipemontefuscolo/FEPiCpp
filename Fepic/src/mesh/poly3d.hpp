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

#ifndef FEPIC_POLY3D_HPP
#define FEPIC_POLY3D_HPP



template<class _Traits>
class _Poly3d : public _CellCore<_Traits>
{
#if !defined(THIS) && !defined(CONST_THIS)
  #define THIS static_cast<CellT*>(this)
  #define CONST_THIS static_cast<const CellT*>(this)
#endif  

public:
  
  typedef typename _Traits::CellT  CellT;
  typedef typename _Traits::HalfT  HalfT;
  typedef typename _Traits::MeshT  MeshT;
  //typedef typename _MetaCellOf<CellT, _Traits>::FaceT FaceT;

protected:      
  _Poly3d(_Poly3d const&) {};
  _Poly3d() {};

public:  
  /** Retorna se os vértices passados formam uma face nesse poliedro.
  * @param[in] nodes um vetor com os vértices. O tamanho do vetor deve ser exatamente o número
  *            de lados (n_borders) do elemento degenerado deste poliedro.
  * @param[out] O índice local da face.
  * @param[out] a âncora.
  * @return true se forma uma face, false caso contrário.
  * @warning A ordem da numeração dos nós é importante!
  * @note definição para ancora: dada uma face n0,n1,...,nB, se ela é a face F
  * e ancora A de um poliedro, então, o nó B-A é o primeiro nó da face F do poliedro.
  */ 
  bool isAnFace(vectori const& vertices, int &face, int &anchor) const
  {
    static const matrixi faces_vtx(_MetaCellOf<CellT, _Traits>::get_faces_vtx());
  
    FEPIC_CHECK(static_cast<int>(vertices.size())==CellT::n_vertices_per_border, "", std::invalid_argument);
    
    vectori temp(CellT::n_vertices_per_border);
    vectori this_vtx(CellT::n_vertices_per_border);
    
    for (int f = 0; f < CellT::n_borders; ++f) // para cada face
    {
      for (int n = 0; n < CellT::n_vertices_per_border; ++n)
      {
        this_vtx[n] = CONST_THIS->_nodes[ faces_vtx[f][n] ];
      }                       
      for (int a = 0; a < CellT::n_vertices_per_border; ++a) // para cada ancora
      {
        rotate_copy(vertices.begin(), vertices.begin()+a, vertices.end(), temp.begin()); // temp = rot(vertices + a);
        
        if(temp == this_vtx)
        {
          face = f;
          anchor = a;
          return true;
        }
              
      }
    }
    
    return false;
  
  }
  
  /** O mesmo que a função isAnFace, mas aqui a ordem dos vértices podem ser anti-paralelos a face.
  */ 
  bool isAnParallelFace(vectori vertices, int &face, int &anchor) const
  {
    if (this->isAnFace(vertices, face, anchor))
      return true;
            
    reverse(vertices.begin(), vertices.end());
    return (this->isAnFace(vertices, face, anchor));
  }
  
  void printSelfState(std::ostream &o) const
  {
    o << this->getNodeIdx(0);
    for (int i = 1, tam=CONST_THIS->_nodes.size(); i < tam; i++)
      o << " " << this->getNodeIdx(i);
  }       
  
#undef THIS
#undef CONST_THIS 
};







#endif
