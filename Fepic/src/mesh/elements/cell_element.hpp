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


#ifndef FEPIC_CELLELEMENT_ELEMENT
#define FEPIC_CELLELEMENT_ELEMENT

#include "../labelable.hpp"

class CellElement : public _Labelable
{
public:  
  
  CellElement(int icell, int position, int tag, int flags, int status = 0)
      :  _Labelable(tag,flags),
         _status(static_cast<char>(status)),
         _position(static_cast<char>(position)),
         _icell(icell)
         {}
  CellElement() : _Labelable(), _status(0), _position(-1), _icell(-1) {}

  
  int getIncidCell() const
  {
    return _icell;
  }
  
  int getPosition() const
  {
    return static_cast<int>(_position);
  }
  
  void setIncidCell(int icell_id)
  {
    _icell = icell_id;
  }

  // facet lid of incident cell
  void setPosition(int pos)
  {
    _position = static_cast<char>( pos );
  }
  
  /// is the same as doing setIncidCell(icell_id); setPosition(pos);
  void setIncidence(int icell_id, int pos)
  {
    _icell = icell_id;
    _position = static_cast<char>( pos );
  }
  
  void copy(CellElement const* other)
  {
    this->setTag(other->getTag());
    this->setFlags(other->getFlags());
    this->setIncidence(other->getIncidCell(), other->getPosition());
  }

  ~CellElement() {}

protected:
  char  _status;     // ONLY POINT CLASS USE IT .. it was put here to avoid problems with padding.
  char  _position;   // local id of this element on the incident cell
  int   _icell;      // global id of the incident cell
  
};



#endif
