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




template<int dim>
class Polytope { public:
  typedef Polytope<dim-1> Derived;
  static std::string name()
  {
    return "Polytope"+std::string(itoa(dim));
  }
};

template<int dim>
class Simplex { public:
  typedef Simplex<dim-1> Derived;
  static std::string name()
  {
    return "Simplex"+std::string(itoa(dim));
  }
};

template<int dim>
class Hypercube { public:
  typedef Hypercube<dim-1> Derived;
  static std::string name()
  {
    return "Hypercube"+std::string(itoa(dim));
  }
};



/*--------------------------------------------------------------------*/

template<> class _MetaCellOf<void, void> {}; // for a sweet syntax highlight


#define FEPIC_CELLOF(obj, type) template<class _Traits>                  \
                                class _MetaCellOf<obj, _Traits> {public: \
                                typedef type<_Traits> Type;}

FEPIC_CELLOF(Simplex<1>, Edge);
FEPIC_CELLOF(Simplex<2>, Triangle);
FEPIC_CELLOF(Simplex<3>, Tetrahedron);
FEPIC_CELLOF(Hypercube<1>, Edge);
FEPIC_CELLOF(Hypercube<2>, Quadrangle);
FEPIC_CELLOF(Hypercube<3>, Hexahedron);
FEPIC_CELLOF(Polytope<1>, Edge);

#undef FEPIC_CELLOF

class UndefElement {
public:
  static int getMshTag()
  {
    return MSH_UNDEFINED_ELEM;
  }
};


//------------------------------------

template<> class _MetaHalfOf<void, void> {}; // for a sweet syntax highlight

#define FEPIC_HALFOF(obj, type) template<class _Traits>                    \
                                class _MetaHalfOf<obj, _Traits> { public:  \
                                typedef type<_Traits> Type; }
                                
FEPIC_HALFOF(Polytope<2>, HalfEdge);
FEPIC_HALFOF(Polytope<3>, HalfFace);
FEPIC_HALFOF(Simplex<2>, HalfEdge);
FEPIC_HALFOF(Simplex<3>, HalfFace);
FEPIC_HALFOF(Hypercube<2>, HalfEdge);
FEPIC_HALFOF(Hypercube<3>, HalfFace);

#undef FEPIC_HALFOF

//------------------------------------

template<> class _MetaHalfLabOf<void, void> {}; // for a sweet syntax highlight

#define FEPIC_HALFLABOF(obj, type) template<class _Traits>                       \
                                   class _MetaHalfLabOf<obj, _Traits> { public:  \
                                   typedef type<_Traits> Type; }

FEPIC_HALFLABOF(Polytope<2>, HalfEdgeLab);
FEPIC_HALFLABOF(Polytope<3>, HalfFaceLab);
FEPIC_HALFLABOF(Simplex<2>, HalfEdgeLab);
FEPIC_HALFLABOF(Simplex<3>, HalfFaceLab);
FEPIC_HALFLABOF(Hypercube<2>, HalfEdgeLab);
FEPIC_HALFLABOF(Hypercube<3>, HalfFaceLab);

#undef FEPIC_HALFLABOF

/** @class DefaultTraits
 * 
 * Default definitions of a _Traits.
 * 
 * Users can do their own _Traits.
 * 
 * In custom _Traits must be defined:
 * - CellT      := type of grid cell
 * - EdgeT
 * - PointT
 * - HalfT
 * - HalfLT
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
  
  typedef DefaultTraits _Traits;
  
  typedef typename _MetaCellOf<CellType, _Traits>::Type CellT;
  
  typedef typename _MetaHalfOf<CellType, _Traits>::Type HalfT;
  
  typedef typename _MetaHalfLabOf<CellType, _Traits>::Type HalfLT;
  
  typedef Point<_Traits>  PointT;
  
  typedef iMesh<_Traits>  MeshT;
  
  enum {spacedim = _spacedim};
  
  //static const int spacedim = _spacedim;
};


/** versão estática.\n
 *  Faz um mapeamento de pontos na célula unitário para a célula real
 *  @param list_pts uma lista com pontos no triângulo unitário nos quais se deseja fazer o mapeamento
 *  @param intp_pts a lista de pontos de interpolação da célula; ex, para interpolação linear, passar os vértices
 *  @param Phi funções de interpolação que correspondem aos pontos de interpolação
 *  @return uma lista com as coordenadas dos pontos passados em list_pts na célula real
 *  @warning as funções Phi DEVEM corresponder aos pontos de interpolação
 */
template<class _Traits, class ShapeFun, 
         int   sdim = _Traits::spacedim,
         int   cdim = _Traits::CellT::Dim,
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





