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

#include "basehandler.hpp"
#include "mesh.hpp"

// constructors 
template<>
FEP_STRONG_INLINE
BaseHandler<Cell>:: BaseHandler(Mesh * mesh, int id) : _elem_ptr(mesh->getCellPtr(id)), _mesh_ptr(mesh), _elem_id(id) {}

template<>
FEP_STRONG_INLINE
BaseHandler<Point>:: BaseHandler(Mesh * mesh, int id) : _elem_ptr(mesh->getNodePtr(id)), _mesh_ptr(mesh), _elem_id(id) {}

template<>
FEP_STRONG_INLINE
BaseHandler<Facet>:: BaseHandler(Mesh * mesh, int id) : _elem_ptr(mesh->getFacetPtr(id)), _mesh_ptr(mesh), _elem_id(id) {}

template<>
FEP_STRONG_INLINE
BaseHandler<Corner>:: BaseHandler(Mesh * mesh, int id) : _elem_ptr(mesh->getCornerPtr(id)), _mesh_ptr(mesh), _elem_id(id) {}



// constructors 
template<>
FEP_STRONG_INLINE
ConstBaseHandler<Cell>:: ConstBaseHandler(Mesh const* mesh, int id) : _elem_ptr(mesh->getCellPtr(id)), _mesh_ptr(mesh), _elem_id(id) {}

template<>
FEP_STRONG_INLINE
ConstBaseHandler<Point>:: ConstBaseHandler(Mesh const* mesh, int id) : _elem_ptr(mesh->getNodePtr(id)), _mesh_ptr(mesh), _elem_id(id) {}

template<>
FEP_STRONG_INLINE
ConstBaseHandler<Facet>:: ConstBaseHandler(Mesh const* mesh, int id) : _elem_ptr(mesh->getFacetPtr(id)), _mesh_ptr(mesh), _elem_id(id) {}

template<>
FEP_STRONG_INLINE
ConstBaseHandler<Corner>:: ConstBaseHandler(Mesh const* mesh, int id) : _elem_ptr(mesh->getCornerPtr(id)), _mesh_ptr(mesh), _elem_id(id) {}
