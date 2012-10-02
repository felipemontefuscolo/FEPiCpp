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

class Mesh;

// EntityType = Cell, Facet, Corner or Point
template<class _Iterator, class _Container, class _Mesh>
class MeshIterator
{
protected:  

  friend class Mesh;

  typedef _Iterator      SeqListIter;
  typedef MeshIterator   Self;
  typedef _Mesh*         MeshPtr;
  typedef _Container     ContainerType;
  
public:
  
  //typedef typename long                  difference_type;
  typedef typename SeqListIter::difference_type   difference_type;
  typedef typename SeqListIter::iterator_category iterator_category;  
  typedef typename SeqListIter::value_type        value_type;
  typedef typename SeqListIter::pointer           pointer;
  typedef typename SeqListIter::reference         reference;

protected:
  //MeshPtr     _mesh_ptr; // DO NOT DELETE ME ... I will be useful later.
  SeqListIter _seq_iter;

public:
  
  template<class _Iter>
  MeshIterator(MeshPtr, _Iter it) : _seq_iter(it) {}
  
  MeshIterator() {}

  // Allow iterator to const_iterator conversion
  template<class _Iter>
  MeshIterator(const MeshIterator<_Iter,_Container, _Mesh> & i)
              : _seq_iter(i.base()) { }  

  SeqListIter const& base() const
  { return _seq_iter; }

public:

  difference_type index() const
  { return _seq_iter.index(); }

  reference
  operator*() const
  { return *_seq_iter; }

  pointer
  operator->() const
  { return &operator*(); }

  Self&
  operator++()
  {
    ++_seq_iter;
    return *this;
  }

  Self
  operator++(int)
  { return Self(_seq_iter++); }

  // Bidirectional iterator requirements
  Self&
  operator--()
  {
    --_seq_iter;
    return *this;
  }

  Self
  operator--(int)
  { return Self(_seq_iter--); }

};

// Forward iterator requirements
template<class _IteratorL, class _IteratorR, class _Container, class _Mesh>
  inline bool
  operator==(const MeshIterator<_IteratorL, _Container, _Mesh>& __lhs,
             const MeshIterator<_IteratorR, _Container, _Mesh>& __rhs)
  { return __lhs.base() == __rhs.base(); }

template<class _Iterator, class _Container, class _Mesh>
  inline bool
  operator==(const MeshIterator<_Iterator, _Container, _Mesh>& __lhs,
             const MeshIterator<_Iterator, _Container, _Mesh>& __rhs)
  { return __lhs.base() == __rhs.base(); }

template<class _IteratorL, class _IteratorR, class _Container, class _Mesh>
  inline bool
  operator!=(const MeshIterator<_IteratorL, _Container, _Mesh>& __lhs,
             const MeshIterator<_IteratorR, _Container, _Mesh>& __rhs)
  { return __lhs.base() != __rhs.base(); }

template<class _Iterator, class _Container, class _Mesh>
  inline bool
  operator!=(const MeshIterator<_Iterator, _Container, _Mesh>& __lhs,
             const MeshIterator<_Iterator, _Container, _Mesh>& __rhs)
  { return __lhs.base() != __rhs.base(); }

// Random access iterator requirements
template<class _IteratorL, class _IteratorR, class _Container, class _Mesh>
  inline bool
  operator<(const MeshIterator<_IteratorL, _Container, _Mesh>& __lhs,
            const MeshIterator<_IteratorR, _Container, _Mesh>& __rhs)
  { return __lhs.base() < __rhs.base(); }

template<class _Iterator, class _Container, class _Mesh>
  inline bool
  operator<(const MeshIterator<_Iterator, _Container, _Mesh>& __lhs,
            const MeshIterator<_Iterator, _Container, _Mesh>& __rhs)
  { return __lhs.base() < __rhs.base(); }

template<class _IteratorL, class _IteratorR, class _Container, class _Mesh>
  inline bool
  operator>(const MeshIterator<_IteratorL, _Container, _Mesh>& __lhs,
            const MeshIterator<_IteratorR, _Container, _Mesh>& __rhs)
  { return __lhs.base() > __rhs.base(); }

template<class _Iterator, class _Container, class _Mesh>
  inline bool
  operator>(const MeshIterator<_Iterator, _Container, _Mesh>& __lhs,
            const MeshIterator<_Iterator, _Container, _Mesh>& __rhs)
  { return __lhs.base() > __rhs.base(); }

template<class _IteratorL, class _IteratorR, class _Container, class _Mesh>
  inline bool
  operator<=(const MeshIterator<_IteratorL, _Container, _Mesh>& __lhs,
             const MeshIterator<_IteratorR, _Container, _Mesh>& __rhs)
  { return __lhs.base() <= __rhs.base(); }

template<class _Iterator, class _Container, class _Mesh>
  inline bool
  operator<=(const MeshIterator<_Iterator, _Container, _Mesh>& __lhs,
             const MeshIterator<_Iterator, _Container, _Mesh>& __rhs)
  { return __lhs.base() <= __rhs.base(); }

template<class _IteratorL, class _IteratorR, class _Container, class _Mesh>
  inline bool
  operator>=(const MeshIterator<_IteratorL, _Container, _Mesh>& __lhs,
             const MeshIterator<_IteratorR, _Container, _Mesh>& __rhs)
  { return __lhs.base() >= __rhs.base(); }

template<class _Iterator, class _Container, class _Mesh>
  inline bool
  operator>=(const MeshIterator<_Iterator, _Container, _Mesh>& __lhs,
             const MeshIterator<_Iterator, _Container, _Mesh>& __rhs)
  { return __lhs.base() >= __rhs.base(); }

// _GLIBCXX_RESOLVE_LIB_DEFECTS
// According to the resolution of DR179 not only the various comparison
// operators but also operator- must accept mixed iterator/const_iterator
// parameters.
template<class _IteratorL, class _IteratorR, class _Container, class _Mesh>
#ifdef __GXX_EXPERIMENTAL_CXX0X__
  // DR 685.
  inline auto
  operator-(const MeshIterator<_IteratorL, _Container, _Mesh>& __lhs,
            const MeshIterator<_IteratorR, _Container, _Mesh>& __rhs)
  -> decltype(__lhs.base() - __rhs.base())
#else
  inline class MeshIterator<_IteratorL, _Container, _Mesh>::difference_type
  operator-(const MeshIterator<_IteratorL, _Container, _Mesh>& __lhs,
            const MeshIterator<_IteratorR, _Container, _Mesh>& __rhs)
#endif
  { return __lhs.base() - __rhs.base(); }

template<class _Iterator, class _Container, class _Mesh>
  inline class MeshIterator<_Iterator, _Container, _Mesh>::difference_type
  operator-(const MeshIterator<_Iterator, _Container, _Mesh>& __lhs,
            const MeshIterator<_Iterator, _Container, _Mesh>& __rhs)
  { return __lhs.base() - __rhs.base(); }


#endif
