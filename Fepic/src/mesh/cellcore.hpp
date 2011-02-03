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

#ifndef FEPIC_CELLCORE_HPP
#define FEPIC_CELLCORE_HPP

template<class _Poly>
class _CellCore : public _Labelable
{
#if !defined(THIS) && !defined(CONST_THIS)
  #define THIS static_cast<_Poly*>(this)
  #define CONST_THIS static_cast<const _Poly*>(this)
  #define 
#endif  

public:
  _CellCore(int tag, bool disabled=false, bool wb_disabled=false, bool useless=false) :
                                      _Labelable(tag, disabled, wb_disabled, useless) {}
  _CellCore() = default;
  _CellCore(_CellCore const&) = default;
  ~_CellCore() = default;
  
  int getNumNodes() const
  {
    return THIS->_node.size();
  }

#undef THIS
#undef CONST_THIS  
protected:

  
};



#endif
