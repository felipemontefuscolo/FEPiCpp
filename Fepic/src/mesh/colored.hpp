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

#ifndef FEPIC_COLORED_HPP
#define FEPIC_COLORED_HPP

#include <algorithm>
#include <tr1/cstdint>
#include "../util/misc2.hpp"
#include "../mesh/enums.hpp"

// fwd
class Mesh;
template<class CT, int SD> class SMesh;
template<class,class,class> class SeqList;



class _Colored
{
public:

  friend class Mesh;
  template<class CT, int SD> friend class SMesh;
  template<class,class,class> friend class SeqList;

  explicit
  _Colored(EColor c) : _same_col_next(-1), _same_col_prev(-1), _color((char)c) {}
  _Colored() : _same_col_next(-1), _same_col_prev(-1), _color(-1) {}
  
  EColor getColor() const {
    return static_cast<EColor>(_color);
  }
  
  int sameColorNext() const {
    return _same_col_next;
  }
  
  int sameColorPrev() const
  {
    return _same_col_prev;
  }
  
protected:
  
  void setColor(EColor c)
  {
    _color = static_cast<char>(c);
  }

  void setColor(char c)
  {
    _color = c;
  }

  void setNext(int n)
  {
    _same_col_next = n;
  }
  
  void setPrev(int p)
  {
    _same_col_prev = p;
  }
  
protected:
  int  _same_col_next;
  int  _same_col_prev;
  char _color;
};


#endif


  
  // NÂO APAGUE
  ///** @brief choose a color that is not in the range [first,last) and assing to the cell.
   //*  @return returns which color was chosen.
   //*/ 
  //template<class T>
  //char setColor(T const* first, T const* last)
  //{
  
    //Uint_color_t all_colors(0);
    //for (T const* p = first; p != last; ++p)
      //if ((*p) >= 0)
        //all_colors |= ( static_cast<Uint_color_t>(1) << (*p) );
    
    //for (Uint_color_t c = 1; ; c<<=1)
    //{
      //// if not found in the range
      //if ( c & all_colors )
        //continue;
      
      //_color = static_cast<char>( log2_i128(c) );
      //return _color;
    //}
  //}


