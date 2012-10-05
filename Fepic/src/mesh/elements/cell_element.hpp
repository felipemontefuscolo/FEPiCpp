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

class CellElement : public Labelable
{
public:  
  
  CellElement(int icell, int position, int tag, int flags, int status = 0)
      :  Labelable(tag,flags),
         m_status(static_cast<char>(status)),
         m_position(static_cast<char>(position)),
         m_icell(icell)
         {}
  CellElement() : Labelable(), m_status(0), m_position(-1), m_icell(-1) {}

  
  int getIncidCell() const
  {
    return m_icell;
  }
  
  int getPosition() const
  {
    return static_cast<int>(m_position);
  }
  
  void setIncidCell(int icell_id)
  {
    m_icell = icell_id;
  }

  // facet lid of incident cell
  void setPosition(int pos)
  {
    m_position = static_cast<char>( pos );
  }
  
  /// is the same as doing setIncidCell(icell_id); setPosition(pos);
  void setIncidence(int icell_id, int pos)
  {
    m_icell = icell_id;
    m_position = static_cast<char>( pos );
  }
  
  void copy(CellElement const* other)
  {
    this->setTag(other->getTag());
    this->setFlags(other->getFlags());
    this->setIncidence(other->getIncidCell(), other->getPosition());
  }

  ~CellElement() {}

protected:
  char  m_status;     // ONLY POINT CLASS USE IT .. it was put here to avoid problems with padding.
  char  m_position;   // local id of this element on the incident cell
  int   m_icell;      // global id of the incident cell
  
};



#endif
