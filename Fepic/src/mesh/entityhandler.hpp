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
template<class Entity_t, class Mesh_t>
class EntityHandler
{
  friend class Mesh;
  
  typedef EntityHandler Self;
  typedef Mesh_t*        MeshPtr;
public:
  
  //typedef typename long                  difference_type;
  typedef Entity_t   value_type;
  typedef Entity_t*  pointer;
  typedef Entity_t&  reference;

protected:
  MeshPtr m_mesh_ptr;
  pointer m_entity_ptr;
  int     m_entity_id;

public:
  
  EntityHandler(MeshPtr mesh, pointer elem, int id) : m_mesh_ptr(mesh), m_entity_ptr(elem), m_entity_id(id) {}
  //EntityHandler(MeshPtr mesh, int id) :  m_mesh_ptr(mesh), m_entity_ptr(mesh->entityPtr<value_type>(id)), m_entity_id(id) {}
  EntityHandler(MeshPtr mesh, int id) :  m_mesh_ptr(mesh), m_entity_ptr(0), m_entity_id(id) {}
  
    // Allow iterator to const_iterator conversion
  template<class Etty_t, class Mse_t>
  EntityHandler(const EntityHandler<Etty_t, Mse_t> & i) : m_mesh_ptr(i.meshPtr()), m_entity_ptr(i.ptr()), m_entity_id(i.index()) {}
  
  EntityHandler() : m_mesh_ptr(0), m_entity_ptr(0), m_entity_id(-1) {}

  MeshPtr const& meshPtr() const
  { return m_mesh_ptr; }

  pointer ptr() const
  {return m_entity_ptr; }

  int index() const
  { return m_entity_id; }

  reference
  operator*() const
  { return *m_entity_ptr; }

  pointer
  operator->() const
  { return m_entity_ptr; }

  bool isValid() const {return ptr() != 0;}

};

template<class EntityL_t, class EntityR_t, class Mesh_t>
  inline bool
  operator==(const EntityHandler<EntityL_t, Mesh_t>& lhs,
             const EntityHandler<EntityR_t, Mesh_t>& rhs)
  { return lhs.ptr() == rhs.ptr(); }

template<class Entity_t, class Mesh_t>
  inline bool
  operator==(const EntityHandler<Entity_t, Mesh_t>& lhs,
             const EntityHandler<Entity_t, Mesh_t>& rhs)
  { return lhs.ptr() == rhs.ptr(); }

template<class EntityL_t, class EntityR_t, class Mesh_t>
  inline bool
  operator!=(const EntityHandler<EntityL_t, Mesh_t>& lhs,
             const EntityHandler<EntityR_t, Mesh_t>& rhs)
  { return lhs.ptr() != rhs.ptr(); }

template<class Entity_t, class Mesh_t>
  inline bool
  operator!=(const EntityHandler<Entity_t, Mesh_t>& lhs,
             const EntityHandler<Entity_t, Mesh_t>& rhs)
  { return lhs.ptr() != rhs.ptr(); }



#endif

