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

#ifndef FEPIC_FACET_HPP
#define FEPIC_FACET_HPP

#include "cell_element.hpp"
#include "../labelable.hpp"

class Facet : public _NodeLessElement
{
public:
  Facet(int ic,
        int pos,
        int tag,
        int flags,
        int bid=-1) : _NodeLessElement(ic,pos,tag,flags), _bound_comp_id(bid)
  {}
  
  void setBoundaryComponentId(int id)
  {
    _bound_comp_id=id;
  }
  int getBoundaryComponentId() const
  {
    return _bound_comp_id;
  }
  
  Facet() : _NodeLessElement(), _bound_comp_id(-1) {}
private:
  int _bound_comp_id;
};

#endif
