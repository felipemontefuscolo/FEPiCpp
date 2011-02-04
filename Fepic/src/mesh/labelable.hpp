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
  enum {tag_size = 32};

protected:
  _Labelable(int tag, bool disabled=false, bool wb_disabled=false, bool useless=false) :
        _tag(static_cast<unsigned char>(tag)), _disabled(disabled), _wb_disabled(wb_disabled), _useless(useless)
  {
    FEPIC_ASSERT((tag>=0)&&(tag<tag_size), "tag number must be less than 32 and greater than 0");
  }
  
  _Labelable() : _tag(0), _disabled(false), _wb_disabled(false), _useless(false) {};

  _Labelable(_Labelable const& L) :
        _tag(L._tag), _disabled(L._disabled), _wb_disabled(L._wb_disabled), _useless(L._useless) {};
  
public:
  int getTag() const
  {
    return static_cast<int>(_tag);
  }
  
  bool disabled() const
  {
    return _disabled;
  }
  
  bool willBeDisabled() const
  {
    return _wb_disabled;
  }
  
  bool useless() const
  {
    return _useless;
  }
  
  void setTag(int tag)
  {
    FEPIC_ASSERT((tag>=0)&&(tag<tag_size), "tag number must be less or equal 31");
    _tag = static_cast<unsigned char>(tag);
  }
  
  void disable(bool dis=true)
  {
    _disabled=dis;
  }
  
  void disableLater(bool dis=true)
  {
    _wb_disabled=dis;
  }
  
  void setUseless(bool ul=true)
  {
    _useless=ul;
  }
  
  void resetFlags()
  {
    _tag=0;
    _disabled = _wb_disabled = _useless = false;
  }

protected:
  unsigned char _tag            : 5; // 0 a 31
  bool          _disabled       : 1;
  bool          _wb_disabled    : 1;
  bool          _useless        : 1;
  
};







#endif

