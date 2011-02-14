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



#undef THIS
#undef CONST_THIS
};







#endif
