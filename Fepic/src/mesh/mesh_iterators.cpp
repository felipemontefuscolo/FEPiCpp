#include "mesh_iterators.hpp"
#include "mesh.hpp"

template<>
typename _MeshIterator<Cell>::Self&
_MeshIterator<Cell>::operator++()
{ _elem_ptr = _mesh_ptr->incEnabledCell(_elem_id); return *this; }

template<>
typename _MeshIterator<Point>::Self&
_MeshIterator<Point>::operator++()
{ _elem_ptr = _mesh_ptr->incEnabledPoint(_elem_id); return *this; }

template<>
typename _MeshIterator<Facet>::Self&
_MeshIterator<Facet>::operator++()
{ _elem_ptr = _mesh_ptr->incEnabledFacet(_elem_id); return *this; }

template<>

typename _MeshIterator<Corner>::Self&
_MeshIterator<Corner>::operator++()
{ _elem_ptr = _mesh_ptr->incEnabledCorner(_elem_id); return *this; }

// =========================================================================

template<>
typename _MeshIterator<Cell>::Self
_MeshIterator<Cell>::operator++(int)
{ Self tmp = *this; _elem_ptr = _mesh_ptr->incEnabledCell(_elem_id); return tmp; }

template<>
typename _MeshIterator<Point>::Self
_MeshIterator<Point>::operator++(int)
{ Self tmp = *this; _elem_ptr = _mesh_ptr->incEnabledPoint(_elem_id); return tmp; }

template<>
typename _MeshIterator<Facet>::Self
_MeshIterator<Facet>::operator++(int)
{ Self tmp = *this; _elem_ptr = _mesh_ptr->incEnabledFacet(_elem_id); return tmp; }

template<>
typename _MeshIterator<Corner>::Self
_MeshIterator<Corner>::operator++(int)
{ Self tmp = *this; _elem_ptr = _mesh_ptr->incEnabledCorner(_elem_id); return tmp; }

// =========================================================================


template<>
typename _MeshIterator<Cell>::Self&
_MeshIterator<Cell>::operator--()
{ _elem_ptr = _mesh_ptr->decEnabledCell(_elem_id); return *this; }

template<>
typename _MeshIterator<Point>::Self&
_MeshIterator<Point>::operator--()
{ _elem_ptr = _mesh_ptr->decEnabledPoint(_elem_id); return *this; }

template<>
typename _MeshIterator<Facet>::Self&
_MeshIterator<Facet>::operator--()
{ _elem_ptr = _mesh_ptr->decEnabledFacet(_elem_id); return *this; }

template<>
typename _MeshIterator<Corner>::Self&
_MeshIterator<Corner>::operator--()
{ _elem_ptr = _mesh_ptr->decEnabledCorner(_elem_id); return *this; }

// =========================================================================

template<> 
typename _MeshIterator<Cell>::Self
_MeshIterator<Cell>::operator--(int)
{ Self tmp = *this; _elem_ptr = _mesh_ptr->decEnabledCell(_elem_id); return tmp; }

template<>
typename _MeshIterator<Point>::Self
_MeshIterator<Point>::operator--(int)
{ Self tmp = *this; _elem_ptr = _mesh_ptr->decEnabledPoint(_elem_id); return tmp; }

template<>
typename _MeshIterator<Facet>::Self
_MeshIterator<Facet>::operator--(int)
{ Self tmp = *this; _elem_ptr = _mesh_ptr->decEnabledFacet(_elem_id); return tmp; }

template<>
typename _MeshIterator<Corner>::Self
_MeshIterator<Corner>::operator--(int)
{ Self tmp = *this; _elem_ptr = _mesh_ptr->decEnabledCorner(_elem_id); return tmp; }


// =========================================================================
// =========================================================================
