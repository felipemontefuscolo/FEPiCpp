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
  
  CellElement(int tag, int flags) :  _Labelable(tag,flags) {}
  CellElement() : _Labelable() {}

  virtual int getIncidCell() const = 0;
  virtual int getPosition() const = 0;
  virtual void setIncidCell(int icell_id) = 0;
  virtual void setPosition(int pos) = 0;
  virtual void setIncidence(int icell_id, int pos) = 0;
  void copy(CellElement const* other)
  {
    this->setTag(other->getTag());
    this->setFlags(other->getFlags());
    this->setIncidence(other->getIncidCell(), other->getPosition());
  }

  ~CellElement() {}

};



// Facets and corners have this
class _NodeLessElement : public CellElement
{
public:  
  explicit _NodeLessElement(int ic=-1,
                               int pos=-1,
                               int tag=0,
                               int flags=0) : CellElement(tag,flags),
                                              _icell_pos(pos),
                                              _icell(ic)
                                              {}
  
  // --- inherited from CellElement --- ///
  
  int getIncidCell() const
  {
    return _icell;
  }
  
  int getPosition() const
  {
    return _icell_pos;
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
  
  /// is the same as doing setIncidCell(icell_id); setPosition(pos);
  void setIncidence(int icell_id, int pos)
  {
    _icell = icell_id;
    _icell_pos = pos;
  }
  // --------------------------------------- //
  
  
  
  ~_NodeLessElement() {}

protected:
  // manter nessa ordem mesmo, para evitar padding
  char  _icell_pos;   // facet lid of incident cell
  int   _icell;
  
};

#endif
