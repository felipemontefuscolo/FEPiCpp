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


#ifndef FEPIC_MESH_ITERATORS_HPP
#define FEPIC_MESH_ITERATORS_HPP

#include <iterator>
#include "basehandler.hpp"

class Mesh;
class Cell;

// EntityType = Cell, Facet, Corner or Point
template<class EntityType>
class _MeshIterator : public BaseHandler<EntityType>
{
  friend class Mesh;
  
  typedef _MeshIterator Self;
public:
  
  //typedef typename long                  difference_type;
  typedef  std::bidirectional_iterator_tag iterator_category;  
  typedef typename BaseHandler<EntityType>::value_type      value_type;
  typedef typename BaseHandler<EntityType>::pointer         pointer;
  typedef typename BaseHandler<EntityType>::const_pointer   const_pointer;
  typedef typename BaseHandler<EntityType>::reference       reference;
  typedef typename BaseHandler<EntityType>::const_reference const_reference;  
  
  _MeshIterator(Mesh * mesh, pointer elem, int id) : BaseHandler<EntityType>(mesh,elem,id) {}
  _MeshIterator(Mesh * mesh, int id) : BaseHandler<EntityType>(mesh,id) {}
  
  _MeshIterator() : BaseHandler<EntityType>() {}

  Self&     operator++();
  Self      operator++(int);
  Self&     operator--();
  Self      operator--(int);

};


// EntityType = Cell, Facet, Corner or Point
template<class EntityType>
class _MeshConstIterator : public ConstBaseHandler<EntityType>
{
  friend class Mesh;
  
  typedef _MeshConstIterator Self;
public:
  
  //typedef typename long                  difference_type;
  typedef  std::bidirectional_iterator_tag iterator_category;  
  typedef  EntityType         value_type;
  typedef  EntityType*        pointer;
  typedef  EntityType const*  const_pointer;
  typedef  EntityType&        reference;
  typedef  EntityType const&  const_reference;
  
  _MeshConstIterator(Mesh const* mesh, const_pointer elem, int id) : ConstBaseHandler<EntityType>(mesh,elem,id) {}
  _MeshConstIterator(Mesh const* mesh, int id) : ConstBaseHandler<EntityType>(mesh,id) {}
  
  _MeshConstIterator() : ConstBaseHandler<EntityType>() {}
  
  Self&           operator++();
  Self            operator++(int);
  Self&           operator--();
  Self            operator--(int);

};




#endif
