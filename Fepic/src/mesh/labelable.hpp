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

#ifndef FEPIC_LABELABLE_HPP
#define FEPIC_LABELABLE_HPP

#include "../util/assert.hpp"

class _Labelable
{
public:
  enum {tag_size = 256,
       flags_size = 256};

  enum Masks {
    mk_disabled = (1<<0),
    mk_marked   = (1<<1),
    mk_visited  = (1<<2)
  };

protected:
  explicit _Labelable(int tag, int flags=0) : _tag(static_cast<unsigned char>(tag)), _flags(static_cast<unsigned char>(flags))
  {
    FEPIC_CHECK((tag>=0)&&(tag<tag_size), "tag number must be less than "+std::string(itoa(tag_size))+" and greater than 0", std::out_of_range);
    FEPIC_CHECK((flags>=0)&&(flags<flags_size), "wrong flags", std::out_of_range);
  }

  _Labelable() : _tag(0), _flags(0) {};

public:
  
  int getTag() const
  {
    return _tag;
  }

  void setTag(int tag)
  {
    FEPIC_CHECK(unsigned(tag)<tag_size, "tag number must be less or equal "+std::string(itoa(tag_size)), std::out_of_range);
    _tag = static_cast<unsigned char>( tag );
  }

  bool disabled() const
  {
    return _flags & mk_disabled;
  }

  void disabled(bool disable_this)
  {
    _flags = static_cast<unsigned char>(   disable_this ? (_flags | static_cast<unsigned char>(mk_disabled)) : (_flags & (~static_cast<unsigned char>(mk_disabled)))  );
  }

  bool marked() const
  {
    return _flags & mk_marked;
  }

  void marked(bool mark_this)
  {
    _flags = static_cast<unsigned char>(  mark_this ? (_flags | static_cast<unsigned char>(mk_marked)) : (_flags & (~static_cast<unsigned char>(mark_this)))  );
  }

  bool visited() const
  {
    return _flags & mk_visited;
  }

  void visited(bool visit)
  {
    _flags = static_cast<unsigned char>( visit ? (_flags | static_cast<unsigned char>(mk_visited)) : (_flags & (~static_cast<unsigned char>(mk_visited)))  );
  }

  bool getFlag(unsigned flag_no) const
  {
    return static_cast<bool>(_flags & ( 1 << flag_no));
  }

  int getFlags() const
  {
    return _flags;
  }

  void setFlag(int flag_no, bool set=true)
  {
    unsigned char const one = static_cast<unsigned char>(1);
    _flags = static_cast<unsigned char>(  set ? (_flags | (one<<flag_no)) : (_flags & (~(one<<flag_no)))  );
  }

  void setFlags(int flags)
  {
    _flags = static_cast<unsigned char>( flags );
  }

  //inline void printFlags() const
  //{
    //for (unsigned i=0; i<sizeof(_flags)*8; ++i)
    //{
      //std::cout << static_cast<bool>(_flags & (1<<i));
    //}
    //std::cout << std::endl;
  //}
  
protected:
  unsigned char _tag; // 0 a 256
  unsigned char _flags; // 8 flags ...

};


static const int DISABLED = _Labelable::mk_disabled;
static const int MARKED   = _Labelable::mk_marked;
static const int VISITED  = _Labelable::mk_visited;




#endif

