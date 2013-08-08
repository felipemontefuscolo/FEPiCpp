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

template<class,class> class SeqList;

class Labelable
{
  template<class,class> friend class SeqList;
  
public:
  enum {tag_size = 256,
       flags_size = 256};

  enum Masks {
    mk_disabled = (1<<0),
    mkm_marked   = (1<<1),
    mk_visited  = (1<<2),
    mk_blocked  = (1<<3),
  };

protected:
  explicit Labelable(int tag, int flags=0) : m_tag(static_cast<unsigned char>(tag)), m_flags(static_cast<unsigned char>(flags))
  {
    FEPIC_CHECK((tag>=0)&&(tag<tag_size), "tag number must be less than "+std::string(itoa(tag_size))+" and greater than 0", std::out_of_range);
    FEPIC_CHECK((flags>=0)&&(flags<flags_size), "wrong flags", std::out_of_range);
  }

  Labelable() : m_tag(0), m_flags(0) {};

public:
  
  int getTag() const
  { return m_tag; }

  void setTag(int tag)
  {
    FEPIC_CHECK(unsigned(tag)<tag_size, "tag number must be less or equal "+std::string(itoa(tag_size)), std::out_of_range);
    m_tag = static_cast<unsigned char>( tag );
  }
  
public:

  bool isDisabled() const
  { return m_flags & mk_disabled; }

protected:
  void setDisabledTo(bool disable_this)
  {
    m_flags = disable_this ? (m_flags | mk_disabled) : (m_flags & (~mk_disabled))  ;
  }
  
public:
  bool isMarked() const
  { return m_flags & mkm_marked; }

  void setMarkedTo(bool mark_this)
  { m_flags = mark_this ? (m_flags | mkm_marked) : (m_flags & (~mkm_marked))  ; }

  bool isVisited() const
  { return m_flags & mk_visited; }

  void setVisitedTo(bool visit)
  { m_flags =  visit ? (m_flags | mk_visited) : (m_flags & (~mk_visited))  ; }

  bool isBlocked() const
  { return m_flags & mk_blocked; }

  void setBlockedTo(bool visit)
  { m_flags =  visit ? (m_flags | mk_blocked) : (m_flags & (~mk_blocked))  ; }

  bool getFlag(unsigned flag_no) const
  { return static_cast<bool>(m_flags & ( 1 << flag_no)); }

  int getFlags() const
  { return m_flags; }

protected:
  void setFlag(int flag_no, bool set=true)
  {
    unsigned char const one = static_cast<unsigned char>(1);
    m_flags =  set ? (m_flags | (one<<flag_no)) : (m_flags & (~(one<<flag_no)))  ;
  }

  void setFlags(int flags)
  {  m_flags = static_cast<unsigned char>( flags ); }

  //inline void printFlags() const
  //{
    //for (unsigned i=0; i<sizeof(m_flags)*8; ++i)
    //{
      //std::cout << static_cast<bool>(m_flags & (1<<i));
    //}
    //std::cout << std::endl;
  //}
  
protected:
  // poderia usar bit-field, mas não é portável
  unsigned char m_tag; // 0 a 256
  unsigned char m_flags; // 8 flags ...

};


static const int DISABLED = Labelable::mk_disabled;
static const int MARKED   = Labelable::mkm_marked;
static const int VISITED  = Labelable::mk_visited;




#endif

