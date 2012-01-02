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

template<class,int> class SMesh;
class Mesh;
class Cell;

// VirtualEntity = Cell, Facet, Corner or Point
template<class VirtualEntity>
class _MeshIterator
{
  template<class,int> friend class SMesh;
  friend class Mesh;
  
  typedef _MeshIterator Self;
public:
  
  //typedef typename long                            difference_type;
  typedef  VirtualEntity                   value_type;
  typedef  VirtualEntity*                  pointer;
  typedef  VirtualEntity&                  reference;
  typedef  std::bidirectional_iterator_tag iterator_category;  
  
  explicit
  _MeshIterator(Mesh * mesh, pointer elem) : _elem_ptr(elem), _mesh_ptr(mesh) {};
  
  _MeshIterator() : _elem_ptr(NULL), _mesh_ptr(NULL) {};
  
  reference operator*() const {return *_elem_ptr;}
  pointer   operator->() const {return &(*_elem_ptr);}
  Self&     operator++();
  Self      operator++(int);
  Self&     operator--();
  Self      operator--(int);

  bool
  operator==(const Self& x) const
  { return _elem_ptr == x._elem_ptr; }

  bool
  operator!=(const Self& x) const
  { return _elem_ptr != x._elem_ptr; }
  
private:
  VirtualEntity * _elem_ptr;
  Mesh * _mesh_ptr;
};





// VirtualEntity = Cell, Facet, Corner or Point
template<class VirtualEntity>
class _MeshColorIterator
{
  template<class,int> friend class SMesh;
  friend class Mesh;
  
  typedef _MeshColorIterator Self;
public:
  
  //typedef typename long                            difference_type;
  typedef  VirtualEntity              value_type;
  typedef  VirtualEntity*             pointer;
  typedef  VirtualEntity&             reference;
  typedef  std::forward_iterator_tag  iterator_category;  
  
  explicit
  _MeshColorIterator(Mesh * mesh, pointer elem) : _elem_ptr(elem), _mesh_ptr(mesh) {};
  
  _MeshColorIterator() : _elem_ptr(NULL), _mesh_ptr(NULL) {};
  
  reference operator*() const {return *_elem_ptr;}
  pointer   operator->() const {return &(*_elem_ptr);}
  Self&     operator++();
  Self      operator++(int);

  bool
  operator==(const Self& x) const
  { return _elem_ptr == x._elem_ptr; }

  bool
  operator!=(const Self& x) const
  { return _elem_ptr != x._elem_ptr; }
  
private:
  VirtualEntity * _elem_ptr;
  Mesh * _mesh_ptr;
};

#endif
