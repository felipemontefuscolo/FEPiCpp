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

// NO DEFINE HERE!

inline Scalar* begin()  { return &this->operator()(0); }
inline Scalar const* begin() const { return &this->operator()(0); }

inline Scalar* end()  { return this->begin()+this->size(); }
inline Scalar const* end() const  { return this->begin()+this->size(); }

inline void push_back(Scalar const& v)
{
  const unsigned old_size = this->size();
  this->conservativeResize(old_size+1);
  this->operator()(old_size) = v;
}

inline void pop_back()
{
  this->conservativeResize(this->size()-1);
}

inline Scalar& back() const {return this->operator()(this->size()-1);}
inline Scalar& front() const {return this->operator()(0);}

inline void reserve() const {/* nothing */};

inline bool empty() const {return this->size()==0;}

