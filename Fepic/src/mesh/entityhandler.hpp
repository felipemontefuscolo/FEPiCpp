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


#ifndef FEPIC_ENTITYHANDLER_HPP
#define FEPIC_ENTITYHANDLER_HPP

#include "../util/macros.hpp"

class Mesh;

// EntityType = Cell, Facet, Corner or Point
template<class _Entity, class _Mesh>
class EntityHandler
{
  friend class Mesh;
  
  typedef EntityHandler Self;
  typedef _Mesh*        MeshPtr;
public:
  
  //typedef typename long                  difference_type;
  typedef _Entity   value_type;
  typedef _Entity*  pointer;
  typedef _Entity&  reference;

protected:
  MeshPtr _mesh_ptr;
  pointer _entity_ptr;
  int     _entity_id;

public:
  
  EntityHandler(MeshPtr mesh, pointer elem, int id) : _mesh_ptr(mesh), _entity_ptr(elem), _entity_id(id) {}
  //EntityHandler(MeshPtr mesh, int id) :  _mesh_ptr(mesh), _entity_ptr(mesh->entityPtr<value_type>(id)), _entity_id(id) {}
  EntityHandler(MeshPtr mesh, int id) :  _mesh_ptr(mesh), _entity_ptr(0), _entity_id(id) {}
  
    // Allow iterator to const_iterator conversion
  template<class _Etty, class _Mse>
  EntityHandler(const EntityHandler<_Etty, _Mse> & i) : _mesh_ptr(i.meshPtr()), _entity_ptr(i.ptr()), _entity_id(i.index()) {}
  
  EntityHandler() : _mesh_ptr(0), _entity_ptr(0), _entity_id(-1) {}

  MeshPtr const& meshPtr() const
  { return _mesh_ptr; }

  pointer ptr() const
  {return _entity_ptr; }

  int index() const
  { return _entity_id; }

  reference
  operator*() const
  { return *_entity_ptr; }

  pointer
  operator->() const
  { return _entity_ptr; }

  bool isValid() const {return ptr() != 0;}

};

template<class _EntityL, class _EntityR, class _Mesh>
  inline bool
  operator==(const EntityHandler<_EntityL, _Mesh>& __lhs,
             const EntityHandler<_EntityR, _Mesh>& __rhs)
  { return __lhs.ptr() == __rhs.ptr(); }

template<class _Entity, class _Mesh>
  inline bool
  operator==(const EntityHandler<_Entity, _Mesh>& __lhs,
             const EntityHandler<_Entity, _Mesh>& __rhs)
  { return __lhs.ptr() == __rhs.ptr(); }

template<class _EntityL, class _EntityR, class _Mesh>
  inline bool
  operator!=(const EntityHandler<_EntityL, _Mesh>& __lhs,
             const EntityHandler<_EntityR, _Mesh>& __rhs)
  { return __lhs.ptr() != __rhs.ptr(); }

template<class _Entity, class _Mesh>
  inline bool
  operator!=(const EntityHandler<_Entity, _Mesh>& __lhs,
             const EntityHandler<_Entity, _Mesh>& __rhs)
  { return __lhs.ptr() != __rhs.ptr(); }

// Random access iterator requirements
//template<class _IteratorL, class _IteratorR, class _Container, class _Mesh>
  //inline bool
  //operator<(const EntityHandler<_IteratorL, _Container, _Mesh>& __lhs,
            //const EntityHandler<_IteratorR, _Container, _Mesh>& __rhs)
  //{ return __lhs.index() < __rhs.index(); }
//
//template<class _Iterator, class _Container, class _Mesh>
  //inline bool
  //operator<(const EntityHandler<_Iterator, _Container, _Mesh>& __lhs,
            //const EntityHandler<_Iterator, _Container, _Mesh>& __rhs)
  //{ return __lhs.index() < __rhs.index(); }
//
//template<class _IteratorL, class _IteratorR, class _Container, class _Mesh>
  //inline bool
  //operator>(const EntityHandler<_IteratorL, _Container, _Mesh>& __lhs,
            //const EntityHandler<_IteratorR, _Container, _Mesh>& __rhs)
  //{ return __lhs.index() > __rhs.index(); }
//
//template<class _Iterator, class _Container, class _Mesh>
  //inline bool
  //operator>(const EntityHandler<_Iterator, _Container, _Mesh>& __lhs,
            //const EntityHandler<_Iterator, _Container, _Mesh>& __rhs)
  //{ return __lhs.index() > __rhs.index(); }
//
//template<class _IteratorL, class _IteratorR, class _Container, class _Mesh>
  //inline bool
  //operator<=(const EntityHandler<_IteratorL, _Container, _Mesh>& __lhs,
             //const EntityHandler<_IteratorR, _Container, _Mesh>& __rhs)
  //{ return __lhs.index() <= __rhs.index(); }
//
//template<class _Iterator, class _Container, class _Mesh>
  //inline bool
  //operator<=(const EntityHandler<_Iterator, _Container, _Mesh>& __lhs,
             //const EntityHandler<_Iterator, _Container, _Mesh>& __rhs)
  //{ return __lhs.index() <= __rhs.index(); }
//
//template<class _IteratorL, class _IteratorR, class _Container, class _Mesh>
  //inline bool
  //operator>=(const EntityHandler<_IteratorL, _Container, _Mesh>& __lhs,
             //const EntityHandler<_IteratorR, _Container, _Mesh>& __rhs)
  //{ return __lhs.index() >= __rhs.index(); }
//
//template<class _Iterator, class _Container, class _Mesh>
  //inline bool
  //operator>=(const EntityHandler<_Iterator, _Container, _Mesh>& __lhs,
             //const EntityHandler<_Iterator, _Container, _Mesh>& __rhs)
  //{ return __lhs.index() >= __rhs.index(); }


#endif

