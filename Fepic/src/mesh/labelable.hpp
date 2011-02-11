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


class _Labelable
{
public:
  enum {tag_size = 256,
        flags_size = 256};

  enum Masks {
    mk_disabled=1,
    mk_marked=2
  };

protected:
  _Labelable(int tag, int flags=0) : _tag(static_cast<unsigned char>(tag)), _flags(static_cast<unsigned char>(flags))
  {
    FEPIC_CHECK((tag>=0)&&(tag<tag_size), "tag number must be less than "+std::string(itoa(tag_size))+" and greater than 0", std::out_of_range);
    FEPIC_CHECK((flags>=0)&&(flags<flags_size), "wrong flags", std::out_of_range);
  }
  
  _Labelable() : _tag(0), _flags(0) {};

public:
  int getTag() const
  {
    return static_cast<int>(_tag);
  }

  void setTag(int tag)
  {
    FEPIC_CHECK((tag>=0)&&(tag<tag_size), "tag number must be less or equal "+std::string(itoa(tag_size)), std::out_of_range);
    _tag = static_cast<unsigned char>(tag);
  }
  
  bool disabled() const
  {
    return _flags & mk_disabled;
  }

  void disabled(bool disable_this)
  {
    _flags = disable_this ? (_flags | mk_disabled) : (_flags & (~mk_disabled));
  }
  
  bool marked() const
  {
    return _flags & mk_marked;
  }

  bool marked(bool mark_this)
  {
    _flags = mark_this ? (_flags | mk_marked) : (_flags & (~mark_this));
  }

  bool getFlag(unsigned flag_no) const
  {
    return ((_flags & ( 1 << flag_no))!=0) ? 1 : 0;
  }
  
  int getFlags() const
  {
    return static_cast<int>(_flags);
  }
    
  void setFlag(uint flag_no, bool set=true)
  {
    _flags = set ? (_flags | (1<<flag_no)) : (_flags & (~(1<<flag_no)));
  }
  
  void setFlags(uint flags)
  {
    _flags = static_cast<unsigned char>(flags);
  }
  
protected:
  unsigned char _tag; // 0 a 256
  unsigned char _flags; // 256 flags ...
  
};







#endif

