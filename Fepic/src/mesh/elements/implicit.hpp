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

#ifndef FEPIC_IMPLICIT_ELEMENT
#define FEPIC_IMPLICIT_ELEMENT

#include "../labelable.hpp"

class ImplicitElement : public _Labelable
{
public:  

  ImplicitElement(int ic=-1,
                  int pos=0,
                  int anch=0,
                  int tag=0,
                  int flags=0) : _Labelable(tag,flags),
                                 _icell_pos(pos),
                                 _icell_anch(anch),
                                 _icell(ic)
                                 {}

  int getIncidCell() const
  {
    return _icell;
  }
  
  int getPosition() const
  {
    return _icell_pos;
  }
  
  int getAnchor() const
  {
    return _icell_anch;
  }

  void setIncidCell(int icell_id)
  {
    _icell = icell_id;
  }

  // facet lid of incident cell
  void setPosition(int pos)
  {
    _icell_pos = pos;
  }
  
  void setAnchor(int anch)
  {
    _icell_anch = anch;
  }

private:
  char  _icell_pos;   // facet lid of incident cell
  char  _icell_anch;
  int   _icell;
  
};


#endif
