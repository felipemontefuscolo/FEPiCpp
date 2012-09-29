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


#ifndef FEPIC_BASEHANDLER_HPP
#define FEPIC_BASEHANDLER_HPP

#include "../util/macros.hpp"

template<class> class BaseHandler;
template<class> class ConstBaseHandler;

class Mesh;
class Cell;


// EntityType = Cell, Facet, Corner or Point
template<class EntityType>
class BaseHandler
{
  friend class Mesh;
  
  typedef BaseHandler Self;
public:
  
  //typedef typename long     difference_type;
  typedef  EntityType         value_type;
  typedef  EntityType*        pointer;
  typedef  EntityType const*  const_pointer;
  typedef  EntityType&        reference;
  typedef  EntityType const&  const_reference;
  
  BaseHandler(Mesh * mesh, pointer elem, int id) : _elem_ptr(elem), _mesh_ptr(mesh), _elem_id(id) {};
  BaseHandler(Mesh * mesh, int id);
  
  BaseHandler() : _elem_ptr(0), _mesh_ptr(0), _elem_id(-1) {};
  
  FEP_STRONG_INLINE
  bool isValid() const {return _elem_ptr != 0;}
  
  FEP_STRONG_INLINE
  reference operator*() {return *(getPtr());}
  
  FEP_STRONG_INLINE
  pointer   operator->() {return &*(getPtr());}
    
  FEP_STRONG_INLINE
  bool
  operator==(const Self& x) const
  { return _elem_ptr == x._elem_ptr; }

  FEP_STRONG_INLINE
  bool
  operator!=(const Self& x) const
  { return _elem_ptr != x._elem_ptr; }

  FEP_STRONG_INLINE
  bool
  operator<(const Self& x) const
  { return _elem_id < x._elem_id; }

  FEP_STRONG_INLINE
  bool
  operator==(ConstBaseHandler<EntityType> const& x) const
  { return _elem_ptr == x.getPtr(); }

  FEP_STRONG_INLINE
  bool
  operator!=(ConstBaseHandler<EntityType> const& x) const
  { return _elem_ptr != x.getPtr(); }

  FEP_STRONG_INLINE
  bool
  operator<(ConstBaseHandler<EntityType> const& x) const
  { return _elem_id < x.getIdx(); }

  FEP_STRONG_INLINE
  pointer getPtr()
  {
    return _elem_ptr;
  }

  FEP_STRONG_INLINE
  const_pointer getPtr() const
  {
    return _elem_ptr;
  }

  FEP_STRONG_INLINE
  int getIdx() const
  {
    return _elem_id;
  }

  Mesh const* getMesh() const
  {
    return _mesh_ptr;
  }
  
protected:
  pointer _elem_ptr;
  Mesh * _mesh_ptr;
  int _elem_id;
};


template<class EntityType>
class ConstBaseHandler
{
  friend class Mesh;
  
  typedef ConstBaseHandler Self;
public:
  
  //typedef typename long     difference_type;
  typedef  EntityType         value_type;
  typedef  EntityType*        pointer;
  typedef  EntityType const*  const_pointer;
  typedef  EntityType&        reference;
  typedef  EntityType const&  const_reference;
  
  ConstBaseHandler(Mesh const* mesh, const_pointer elem, int id) : _elem_ptr(elem), _mesh_ptr(mesh), _elem_id(id) {};
  ConstBaseHandler(Mesh const* mesh, int id);
  
  ConstBaseHandler() : _elem_ptr(0), _mesh_ptr(0), _elem_id(-1) {};

  ConstBaseHandler(BaseHandler<EntityType> const& h): _elem_ptr(h.getPtr()), _mesh_ptr(h.getMesh()), _elem_id(h.getIdx()) {};
  
  FEP_STRONG_INLINE
  const_reference operator*() const {return *this->getPrt();}
  
  FEP_STRONG_INLINE
  const_pointer   operator->() const {return &(*this->getPrt());}  
  
  FEP_STRONG_INLINE
  Self&
  operator=(BaseHandler<EntityType> const& x)
  {
    _elem_ptr = x.getPtr();
    _mesh_ptr = x.getMesh();
    _elem_id  = x.getIdx();
    return *this;
  }

  FEP_STRONG_INLINE
  bool
  operator!=(const Self& x) const
  { return _elem_ptr != x._elem_ptr; }

  FEP_STRONG_INLINE
  bool
  operator<(const Self& x) const
  { return _elem_id < x._elem_id; }

  FEP_STRONG_INLINE
  bool
  operator==(BaseHandler<EntityType> const& x) const
  { return _elem_ptr == x.getPtr(); }

  FEP_STRONG_INLINE
  bool
  operator!=(BaseHandler<EntityType> const& x) const
  { return _elem_ptr != x.getPtr(); }

  FEP_STRONG_INLINE
  bool
  operator<(BaseHandler<EntityType> const& x) const
  { return _elem_id < x.getIdx(); }
  

  FEP_STRONG_INLINE
  const_pointer getPtr() const
  {
    return _elem_ptr;
  }

  FEP_STRONG_INLINE
  int getIdx() const
  {
    return _elem_id;
  }

  FEP_STRONG_INLINE
  Mesh const* getMesh() const
  {
    return _mesh_ptr;
  }
  
protected:
  const_pointer _elem_ptr;
  Mesh const* _mesh_ptr;
  int _elem_id;
};




#endif

