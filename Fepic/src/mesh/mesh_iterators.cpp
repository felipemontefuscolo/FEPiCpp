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


template<> 
typename _MeshColorIterator<Cell>::Self&
_MeshColorIterator<Cell>::operator++()
{ _elem_ptr = _mesh_ptr->incColorCell(_elem_ptr); return *this; }

template<> 
typename _MeshColorIterator<Point>::Self&
_MeshColorIterator<Point>::operator++()
{ _elem_ptr = _mesh_ptr->incColorPoint(_elem_ptr); return *this; }

template<> 
typename _MeshColorIterator<Facet>::Self&
_MeshColorIterator<Facet>::operator++()
{ _elem_ptr = _mesh_ptr->incColorFacet(_elem_ptr); return *this; }

template<> 
typename _MeshColorIterator<Corner>::Self&
_MeshColorIterator<Corner>::operator++()
{ _elem_ptr = _mesh_ptr->incColorCorner(_elem_ptr); return *this; }

// =========================================================================

template<> 
typename _MeshColorIterator<Cell>::Self
_MeshColorIterator<Cell>::operator++(int)
{ Self tmp = *this; _elem_ptr = _mesh_ptr->incColorCell(_elem_ptr); return tmp; }

template<> 
typename _MeshColorIterator<Point>::Self
_MeshColorIterator<Point>::operator++(int)
{ Self tmp = *this; _elem_ptr = _mesh_ptr->incColorPoint(_elem_ptr); return tmp; }

template<> 
typename _MeshColorIterator<Facet>::Self
_MeshColorIterator<Facet>::operator++(int)
{ Self tmp = *this; _elem_ptr = _mesh_ptr->incColorFacet(_elem_ptr); return tmp; }

template<> 
typename _MeshColorIterator<Corner>::Self
_MeshColorIterator<Corner>::operator++(int)
{ Self tmp = *this; _elem_ptr = _mesh_ptr->incColorCorner(_elem_ptr); return tmp; }
















//template<class T> 
//typename _MeshIterator<T>::Self
//_MeshIterator<T>::operator++(int)
//{
//  Self tmp = *this;
//  const int next = _iter_to_t->sameColorNext();
//  if (next < 0)
//  {
//    _iter_to_t = _ptr_to_seq->_data.end();
//    return tmp;
//  }
//  _iter_to_t = ContainerIterator( &(_ptr_to_seq->_data[next]));
//  return tmp;
//}
//
//template<class T> 
//typename _MeshIterator<T>::Self&
//_MeshIterator<T>::operator--()
//{
//  const int next = _iter_to_t->sameColorNext();
//  if (next < 0)
//  {
//    _iter_to_t = _ptr_to_seq->_data.end();
//    return *this;
//  }
//  _iter_to_t = ContainerIterator( &(_ptr_to_seq->_data[next]));
//  return *this;
//}
//
//template<class T> 
//typename _MeshIterator<T>::Self
//_MeshIterator<T>::operator--(int)
//{
//  Self tmp = *this;
//  const int next = _iter_to_t->sameColorNext();
//  if (next < 0)
//  {
//    _iter_to_t = _ptr_to_seq->_data.end();
//    return tmp;
//  }
//  _iter_to_t = ContainerIterator( &(_ptr_to_seq->_data[next]));
//  return tmp;
//}


