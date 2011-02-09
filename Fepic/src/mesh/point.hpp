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

#ifndef FEPIC_POINT_HPP
#define FEPIC_POINT_HPP

/** 
 *  Os objetos desta classe representam pontos no espaço.
 */ 
template<class _Traits>
class Point : public _Labelable {
public:
  enum {dim=0};
  
  typedef typename _Traits::CellT                         CellT;
  typedef typename _MetaHalfOf<typename CellT::PolytopeT, _Traits>::Type    HalfT;
      
  typedef Edge<_Traits>          EdgeT;
  typedef typename _Traits::MeshT MeshT;
  
  typedef Eigen::Matrix<double, _Traits::spacedim, 1> VecT;
    
  /** Construtor.
  *  @param coord um vetor com Dim elementos que armazena a coordenada.
  *  @param label seu rótulo.
  */ 
  template<class T>
  Point(T const& coord, char label=0) : _Labelable(label)
  {
    for (int i = 0; i < _Traits::spacedim; ++i)
      _coord[i] = coord[i];
  }
  
  /** Construtor.
  */ 
  Point() : _Labelable() {}
  // construtor de cópia não necessário, pois não há nenhum ponteiro.
  
  /** @return a dimensão do espaço.
  */  
  int getSpaceDim() const 
  {
    return _Traits::spacedim;
  }
  
  /** Faz com que este ponto tenha a HalfT ha.
  */ 
  void setHalf(HalfT const& ha)
  {
    _half = ha;
  }
  
  /** Retorna um ponteiro para sua HalfT.
  */ 
  HalfT* getHalf()
  {
    return &_half;
  }
  
  /** Retorna um ponteiro para sua HalfT.
  */ 
  const HalfT* getHalf() const
  {
    return &_half;
  }
          
  /** Imprime as coordenadas de um ponto em um stream dado.
  *  @param o stream onde se vai imprimir.
  *  @param space Espaço entre a impressão de cada dimensão da coordenada.
  */ 
  void printSelfVtk(std::ostream &o, int space = 22) const
  {
    switch (_Traits::spacedim) {
      case 1:
      {
        o << std::left << std::setw(space) << _coord[0] << " 0.0 0.0";
        break;
      }
      case 2:
      {
        o << std::left << std::setw(space) << _coord[0]
                       << std::setw(space) << _coord[1] << " 0.0";
        break;
      }
      case 3:
      {
        o << std::left << std::setw(space) << _coord[0]
                       << std::setw(space) << _coord[1]
                       << std::setw(space) << _coord[2];
        break;
      }
      default:
      {
        std::cout << "Error: invalid dimension\n";
        break;
      }
    }
  }
  
  /** Define a coordenada deste ponto.
  *  @param coord um vetor com a coordenada.
  */ 
  template<class Vec>
  void setCoord(Vec const& coord) 
  {
    for (int i = 0; i < _Traits::spacedim; ++i)
      _coord[i] = coord[i];
  }
  
  /** Retorna a coordenada deste ponto em coord.
  *  @param[out] coord a coordenada.
  */ 
  template<class Vec>
  void getCoord(Vec & coord) const 
  {
    for (int i = 0; i < _Traits::spacedim; ++i)
      coord[i] = _coord[i];
  }
  
  /** Retorna a i-ésima componente da coordenada
  */ 
  double getCoord(int i) const
  {
    return _coord[i];
  }
  
  /** Retorna a coordenada deste ponto.
  */
  VecT getCoord() const
  {
    return _coord;
  }
  
  /** Retorna a distância até a coordenada p.
  */     
  double getDistance(VecT const& p) const
  {
    return sqrt( (_coord-p).dot(_coord-p) );
  }
  
  /** Retorna a distância até o ponto p.
  */  
  double getDistance(Point const& p) const
  {
    return sqrt( (_coord-p._coord).dot(_coord-p._coord) );
  }
  
  
  /** Retorna o tag correspondente ao formato de arquivo .msh
  */ 
  static int getMshTag()
  {
    return MSH_PNT;
  }
  
  /** NOT FOR USERS
  */ 
  static void setOrder()
  {
    // nada
  }
  
  /** Destrutor.
  */ 
  ~Point() {}
  
  static const int n_borders = 1;
protected:
  VecT _coord;
  HalfT _half;
};








#endif
