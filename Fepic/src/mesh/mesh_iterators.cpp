#include "mesh_iterators.hpp"
#include "mesh.hpp"

template<> 
typename _MeshIterator<Cell>::Self&
_MeshIterator<Cell>::operator++()
{ _elem_ptr = _mesh_ptr->incCell(_elem_ptr); return *this; }

template<> 
typename _MeshIterator<Point>::Self&
_MeshIterator<Point>::operator++()
{ _elem_ptr = _mesh_ptr->incPoint(_elem_ptr); return *this; }

template<> 
typename _MeshIterator<Facet>::Self&
_MeshIterator<Facet>::operator++()
{ _elem_ptr = _mesh_ptr->incFacet(_elem_ptr); return *this; }

template<> 
typename _MeshIterator<Corner>::Self&
_MeshIterator<Corner>::operator++()
{ _elem_ptr = _mesh_ptr->incCorner(_elem_ptr); return *this; }

// =========================================================================

template<> 
typename _MeshIterator<Cell>::Self
_MeshIterator<Cell>::operator++(int)
{ Self tmp = *this; _elem_ptr = _mesh_ptr->incCell(_elem_ptr); return tmp; }

template<> 
typename _MeshIterator<Point>::Self
_MeshIterator<Point>::operator++(int)
{ Self tmp = *this; _elem_ptr = _mesh_ptr->incPoint(_elem_ptr); return tmp; }

template<> 
typename _MeshIterator<Facet>::Self
_MeshIterator<Facet>::operator++(int)
{ Self tmp = *this; _elem_ptr = _mesh_ptr->incFacet(_elem_ptr); return tmp; }

template<> 
typename _MeshIterator<Corner>::Self
_MeshIterator<Corner>::operator++(int)
{ Self tmp = *this; _elem_ptr = _mesh_ptr->incCorner(_elem_ptr); return tmp; }


// =========================================================================


template<> 
typename _MeshIterator<Cell>::Self&
_MeshIterator<Cell>::operator--()
{ _elem_ptr = _mesh_ptr->decCell(_elem_ptr); return *this; }

template<> 
typename _MeshIterator<Point>::Self&
_MeshIterator<Point>::operator--()
{ _elem_ptr = _mesh_ptr->decPoint(_elem_ptr); return *this; }

template<> 
typename _MeshIterator<Facet>::Self&
_MeshIterator<Facet>::operator--()
{ _elem_ptr = _mesh_ptr->decFacet(_elem_ptr); return *this; }

template<> 
typename _MeshIterator<Corner>::Self&
_MeshIterator<Corner>::operator--()
{ _elem_ptr = _mesh_ptr->decCorner(_elem_ptr); return *this; }

// =========================================================================

template<> 
typename _MeshIterator<Cell>::Self
_MeshIterator<Cell>::operator--(int)
{ Self tmp = *this; _elem_ptr = _mesh_ptr->decCell(_elem_ptr); return tmp; }

template<> 
typename _MeshIterator<Point>::Self
_MeshIterator<Point>::operator--(int)
{ Self tmp = *this; _elem_ptr = _mesh_ptr->decPoint(_elem_ptr); return tmp; }

template<> 
typename _MeshIterator<Facet>::Self
_MeshIterator<Facet>::operator--(int)
{ Self tmp = *this; _elem_ptr = _mesh_ptr->decFacet(_elem_ptr); return tmp; }

template<> 
typename _MeshIterator<Corner>::Self
_MeshIterator<Corner>::operator--(int)
{ Self tmp = *this; _elem_ptr = _mesh_ptr->decCorner(_elem_ptr); return tmp; }

// =========================================================================
// =========================================================================




