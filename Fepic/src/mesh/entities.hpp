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

#ifndef FEPIC_ENTITIES_HPP    
#define FEPIC_ENTITIES_HPP






/* pre declarações */


/*--------------------------------------------------------------------*/



template<class Traits>
class ElementProperties<Simplex<1>, Traits> {
public:
  typedef Edge<Traits> Type;
};

template<class Traits>
class ElementProperties<Simplex<2>, Traits> {
public:
  typedef Triangle<Traits> Type;
};

template<class Traits>
class ElementProperties<Simplex<3>, Traits> {
public:
  typedef Tetrahedron<Traits> Type;
};

template<class Traits>
class ElementProperties<Hypercube<1>, Traits> {
public:
  typedef Edge<Traits> Type;
};

template<class Traits>
class ElementProperties<Hypercube<2>, Traits> {
public:
  typedef iQuadrangle<Traits> Type;
};

template<class Traits>
class ElementProperties<Hypercube<3>, Traits> {
public:
  typedef iHexahedron<Traits> Type;
};

template<class Traits>
class ElementProperties<Edge<Traits>, Traits> {
public:
  static const int n_borders = 2;
  static const int n_vertices = 2;
};
template<class Traits>
class ElementProperties<Triangle<Traits>, Traits> {
public:
  static const int n_borders = 3;
  static const int n_vertices = 3;
};
template<class Traits>
class ElementProperties<Tetrahedron<Traits>, Traits> {
public: static const int n_borders = 4;
  static const int n_vertices = 4;
  static Fepic::matrixi get_faces_vtx()
  {
    static const Fepic::matrixi temp = { {1,0,2}, {0,1,3}, {3,2,0}, {2,3,1} };
    return temp;
  }
  static Fepic::matrixi get_edges_vtx()
  {
    static const Fepic::matrixi temp = { {0,1}, {1,2}, {2,0}, {3,0}, {3,2}, {3,1} };
    return temp;
  }
  
  typedef Triangle<Traits> FaceT;
};


class UndefElement {
public:
  static int getMshTag()
  {
    return Msh_UNDEFINED_ELEM;
  }
};

/* Define type of volume: Cell::dim < 3 ? UndefVol : Cell::Volume */

template<class Cell>
class VolumeDef<1, Cell> {
public:
  typedef UndefElement VolumeT;
}; 

template<class Cell>
class VolumeDef<2, Cell> {
public:
  typedef UndefElement VolumeT;
};

template<class Cell>
class VolumeDef<3, Cell> {
public:
  typedef typename Cell::VolumeT VolumeT;
};



template<class Cell>
class FaceDef<1, Cell> {
public:
  typedef UndefElement FaceT;
};

template<class Cell>
class FaceDef<2, Cell> {
public:
  typedef typename Cell::FaceT FaceT;
};

template<class Cell>
class FaceDef<3, Cell> {
public:
  typedef typename Cell::FaceT FaceT;
};

//------------------------------------



template<class Traits>
class HalfDef<1, Traits> {
public: 
  typedef HalfEdge<Traits> HalfT;
};

template<class Traits>
class HalfDef<2, Traits> {
public: 
  typedef HalfEdge<Traits> HalfT;
};

template<class Traits>
class HalfDef<3, Traits> {
public: 
  typedef HalfFace<Traits> HalfT;
};


//------------------------------------



template<class Traits>
class HalflDef<1, Traits> {
public: 
  typedef HalfEdgeLab<Traits> HalflT;
};

template<class Traits>
class HalflDef<2, Traits> {
public: 
  typedef HalfEdgeLab<Traits> HalflT;
};

template<class Traits>
class HalflDef<3, Traits> {
public: 
  typedef HalfFaceLab<Traits> HalflT;
};


/** @class DefaultTraits
 * 
 * Default definitions of a traits.
 * 
 * Users can do their own traits.
 * 
 * In custom traits must be defined:
 * - CellT      := type of grid cell
 * - EdgeT
 * - PointT
 * - MeshT
 * - spacedim   := dimension of the space
 * 
 * 
 */
template<int _spacedim, class CellType = Simplex<_spacedim> >
class DefaultTraits {
public:

  DefaultTraits(DefaultTraits const&) = delete; // dont copy me
  ~DefaultTraits() = delete;
  
  typedef DefaultTraits Traits;
  
  typedef typename ElementProperties<CellType, Traits>::Type CellT;
  
  typedef Point<Traits>  PointT;
  
  typedef iMesh<Traits>   MeshT;
  
  static const int spacedim = _spacedim;
};


/** versão estática.\n
 *  Faz um mapeamento de pontos na célula unitário para a célula real
 *  @param list_pts uma lista com pontos no triângulo unitário nos quais se deseja fazer o mapeamento
 *  @param intp_pts a lista de pontos de interpolação da célula; ex, para interpolação linear, passar os vértices
 *  @param Phi funções de interpolação que correspondem aos pontos de interpolação
 *  @return uma lista com as coordenadas dos pontos passados em list_pts na célula real
 *  @warning as funções Phi DEVEM corresponder aos pontos de interpolação
 */
template<class Traits, class ShapeFun, 
         int   sdim = Traits::spacedim,
         int   cdim = Traits::CellT::Dim,
         class VecT = Eigen::Matrix<double, sdim, 1>,   // vetor no espaço da célula real
         class VecU = Eigen::Matrix<double, cdim, 1> >  // vetor na espaço da célula unitária
std::vector<VecT> map2RealCell(std::vector<VecU> const& list_pts,
                               std::vector<VecT> const& intp_pts,
                               ShapeFun          const& Phi)
{
  int tam = static_cast<int>( list_pts.size() );
  std::vector<VecT> ret(tam, VecT::Zero());
  
  for (int i = 0; i < tam; ++i)
    for (int k = 0, size=static_cast<int>( intp_pts.size() ); k < size; ++k)
      ret[i] += Phi(list_pts[i], k)*intp_pts[k];
  
  return ret;
}











#endif





