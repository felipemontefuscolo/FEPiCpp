#ifndef FEPIC_LIST_TYPE_HPP
#define FEPIC_LIST_TYPE_HPP

#include <vector>
#include <map>
#include <set>
#include <utility>
#include <iterator>
#include <tr1/type_traits>
#ifdef FEP_HAS_OPENMP
#  include "omp.h"
#endif
#include "../mesh/enums.hpp"
#include "../util/misc.hpp"
#include "contrib/Loki/set_vector.hpp"

//#include "../mesh/labelable.hpp"  //nao precisa
#include <iostream>
#include "macros.hpp"


#include "boost/static_assert.hpp"
#include "boost/utility/enable_if.hpp"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wundef"
#include "boost/type_traits/remove_pointer.hpp"
#include "boost/type_traits/is_same.hpp"
#pragma GCC diagnostic pop

//
//                      _      _  _       _
//                     | |    | |(_) ___ | |_
//                     | |    | || |/ __|| __|
//                     | |___ | || |\__ \| |_
//                     |_____||_||_||___/ \__|
//
//
//

// fwd
template<class,class> class SeqList_iterator;


/// @brief A container with a specific use for the mesh. The value type of this container is
/// deducted from the container passed as template argument "C". This value type should have two
/// functions: "bool isDisabled() const" and "void setDisabledTo(bool)", see _Labelable class.
/// These functions are used to "delete" containers elements. If these functions don't have
/// public access, make SeqList as a friend class.
///
template<class C,                      ///< A random access data container: std::vector, boost::ptr_vector, etc.
                                       ///< it should have public access to its value type.
                                       /// for store pointers, DO NOT USE std::vector<T*>. Use boost::ptr_vector<T> instead.
         class S  = SetVector<int> >   ///< A sorted container for internal indices. For most
                                       ///< purposes you will never have to change this default.
class SeqList
{
  template<class,class> friend class SeqList_iterator;

  typedef          SeqList<C,S>                          Self;
  typedef typename C::iterator                           DataIterator;
  typedef typename C::const_iterator                     DataConstIterator;
  typedef typename C::reverse_iterator                   DataReverseIterator;

public:

  // stl standard
  typedef          C                                      container_type;
  typedef typename container_type::value_type             value_type;
  typedef typename container_type::allocator_type         allocator_type;
  typedef          S                                      ids_container_type;
  typedef typename container_type::reference              reference;
  typedef typename container_type::const_reference        const_reference;
  typedef typename container_type::pointer                pointer;
  typedef typename container_type::size_type              size_type;
  typedef typename container_type::difference_type        difference_type;
  
  typedef SeqList_iterator<DataIterator,Self>             iterator;
  typedef SeqList_iterator<DataConstIterator,Self const>  const_iterator;

private:
  
  // auxiliary typedefs
  typedef typename boost::remove_pointer<value_type>::type RP_ValueType; // remove pointer value type
  typedef          RP_ValueType &                          RP_Reference;
  typedef          RP_ValueType const&                     RP_ConstReference;
  typedef          RP_ValueType *                          RP_Pointer;
  typedef          RP_ValueType const*                     RP_ConstPointer;

public:

  template<class T> // T = value_type
  struct _InsertFunArgument {
    typedef T const& type;
  };
  
  template<class T> // T* = value_type
  struct _InsertFunArgument<T*> {
    typedef T*  type;
  };

  explicit SeqList()
  {
    _data.clear();
    _actived_beg = _data.begin();
  }

  void clear()
  {
    _data.clear();
    _disabled_idcs.clear();
  }

  //// TODO: only to vector type ...
  //int getDataId(const_pointer v) const
  //{
  //  return static_cast<int>(std::distance(_data.begin(), DataConstIterator(v)));
  //}

  size_type size() const
  {return _data.size() - _disabled_idcs.size();};

  size_type totalSize() const
  {return _data.size();};

  //  SERIAL VERSION 

  iterator begin()
  { return iterator(this, _actived_beg); }

  const_iterator begin() const
  { return const_iterator(this, _actived_beg); }  

  iterator end()
  { return iterator(this, _data.end()); }

  const_iterator end() const
  { return const_iterator(this, _data.end()); }



  //  PARALLEL VERSION

  iterator begin(int tid, int nthreads, int * begin_idx = NULL)
  {
    int const N = size();
    int const start = tid*(N/nthreads) + ((N%nthreads) < tid ? (N%nthreads) : tid);
    if(begin_idx)
      *begin_idx = start;
    iterator s = begin();
    for (int i = 0; i < start; ++i)
      ++s;
    return s;
  }

  iterator end(int tid, int nthreads, int * end_idx = NULL)
  {
    int const N = size();
    int const start = tid*(N/nthreads) + ((N%nthreads) < tid ? (N%nthreads) : tid);
    int const end_   = start + N/nthreads + ((N%nthreads) > tid);
    if (end_idx)
      *end_idx = end_;
    iterator e = begin(tid,nthreads);
    for (int i = 0; i < end_-start; ++i)
      ++e;
    return e;
  }

  const_iterator begin(int tid, int nthreads, int * begin_idx = NULL) const
  {
    int const N = size();
    int const start = tid*(N/nthreads) + ((N%nthreads) < tid ? (N%nthreads) : tid);
    if(begin_idx)
      *begin_idx = start;
    const_iterator s = begin();
    for (int i = 0; i < start; ++i)
      ++s;
    return s;
  }

  const_iterator end(int tid, int nthreads, int * end_idx = NULL) const
  {
    int const N = size();
    int const start = tid*(N/nthreads) + ((N%nthreads) < tid ? (N%nthreads) : tid);
    int const end_   = start + N/nthreads + ((N%nthreads) > tid);
    if (end_idx)
      *end_idx = end_;
    const_iterator e = begin(tid,nthreads);
    for (int i = 0; i < end_-start; ++i)
      ++e;
    return e;
  }


  // modifiers

  void disable(const_iterator it)
  { disable(it.index());};

  void disable(int del_id)
  {
    //DataIterator it = DataIterator(&_data[del_id]);
    RP_Pointer it = &_data.at(del_id);
    if (it->isDisabled())
    {
      //std::cout << "ERRROR: trying to disable a disabled element \n";
      //throw;
      return;
    }
    _disabled_idcs.insert(del_id);
    it->setDisabledTo(true);

    if (it == &*_actived_beg)
      _update_member_beg();
  }

  /** @brief insert an element and return your id.
   */
  int insert(typename _InsertFunArgument<value_type>::type obj)
  {
    //_InsertFunArgument<value_type>::type 
    //    = T const& for non-pointer containers
    //    = T*   for pointer containers
    return insert_impl<value_type>(obj);
  }

  // com dor no coracao
  //void resize(size_type s)
  //{
  //  _data.resize(s);
  //  _update_member_beg();
  //}

  RP_Reference operator[](size_type n)
  { return _data[n]; }

  RP_ConstReference operator[](size_type n) const
  { return _data[n]; }

  int contiguousId(int id) const
  { return id - static_cast<int>( std::distance(_disabled_idcs.begin(), std::lower_bound(_disabled_idcs.begin(), _disabled_idcs.end(), id)));  }
  
  /// A shortcut to get multiples cids at same time.
  /// @param[out] cids contiguous ids.
  ///  
  void contiguousIds(int *ids_beg, int const* ids_end, int *cids) const
  {
    while (ids_beg != ids_end)
    { *cids++ = contiguousId(*ids_beg++); }
  }
    
protected:

  template<class Value_type>
  int insert_impl(const_reference obj, typename EnableIf< ! std::tr1::is_pointer<Value_type>::value >::type * = NULL)
  {
    if (_disabled_idcs.empty())
    {
      // --- push_back ----
      _data.push_back(obj);
      _update_member_beg();
      // -------------------
      return _data.size()-1;
    }

    int const new_id = _disabled_idcs.back();
    _disabled_idcs.pop_back();
    _data[new_id] = obj;

    _update_member_beg();
    return new_id;

  }

  template<class Value_type>
  int insert_impl(value_type obj, typename EnableIf< std::tr1::is_pointer<Value_type>::value >::type * = NULL)
  {
    if (_disabled_idcs.empty())
    {
      // --- push_back ----
      _data.push_back(obj);
      _update_member_beg();
      // -------------------
      return _data.size()-1;
    }

    int const new_id = _disabled_idcs.back();
    _disabled_idcs.pop_back();
    _data[new_id] = *obj;

    _update_member_beg();
    return new_id;
  }

  void _update_member_beg()
  {
    _actived_beg = _data.begin();
    while (_actived_beg != _data.end() && _actived_beg->isDisabled())
      ++_actived_beg;
  }

  container_type      _data;
  ids_container_type  _disabled_idcs; // sorted vector
  DataIterator        _actived_beg;   // iterator to the beginning of valid data

};




// an iterator adaptor
template<typename _Iterator, typename _Container>
class SeqList_iterator
{

  template<class,class> friend class SeqList;
  template<class> friend class SeqList_const_iterator;

  typedef _Container*      SeqListPtr;
  typedef SeqList_iterator Self;
  typedef _Iterator        DataIterator;
public:

  // stl standard
  typedef typename DataIterator::difference_type   difference_type;
  typedef typename DataIterator::value_type        value_type;
  typedef typename DataIterator::pointer           pointer;
  typedef typename DataIterator::reference         reference;
  typedef typename DataIterator::iterator_category iterator_category;

protected:

  // attributes
  DataIterator _data_iter;
  SeqListPtr   _seq_ptr;

public:

  SeqList_iterator(SeqListPtr sq, DataIterator x) : _data_iter(x), _seq_ptr(sq) {}

  SeqList_iterator() : _data_iter(), _seq_ptr(NULL) {}

  // Allow iterator to const_iterator conversion
  template<class _Iter>
  SeqList_iterator(const SeqList_iterator<_Iter,_Container> & __i)
  : _data_iter(__i.base()),  _seq_ptr(__i._seq_ptr) { }

  DataIterator const& base() const
  { return _data_iter; }

  difference_type index() const
  { return base() - _seq_ptr->_data.begin();}

  reference
  operator*() const
  {return *_data_iter; }

  pointer
  operator->() const
  { return &(operator*()); }
  
  Self&
  operator++()
  {
    ++_data_iter;
    while(_data_iter != _seq_ptr->_data.end() && _data_iter->isDisabled())
      ++_data_iter;
    return *this;
  }

  Self
  operator++(int)
  {
    Self tmp = *this;
    ++_data_iter;
    while(_data_iter != _seq_ptr->_data.end() && _data_iter->isDisabled())
      ++_data_iter;
    return tmp;
  }

  Self&
  operator--()
  {
    --_data_iter;
    while(_data_iter != _seq_ptr->_data.begin() && _data_iter->isDisabled())
      --_data_iter;
    return *this;
  }

  Self
  operator--(int)
  {
    Self tmp = *this;
    --_data_iter;
    while(_data_iter != _seq_ptr->_data.begin() && _data_iter->isDisabled())
      --_data_iter;
    return tmp;
  }
 
};

  //
  // inspired by /usr/include/c++/4.6.3/bits/stl_iterator.h
  //

  // Note: In what follows, the left- and right-hand-side iterators are
  // allowed to vary in types (conceptually in cv-qualification) so that
  // comparison between cv-qualified and non-cv-qualified iterators be
  // valid.  However, the greedy and unfriendly operators in std::rel_ops
  // will make overload resolution ambiguous (when in scope) if we don't
  // provide overloads whose operands are of the same type.  Can someone
  // remind me what generic programming is about? -- Gaby

  // Forward iterator requirements
  template<typename _IteratorL, typename _IteratorR, typename _Container>
    inline bool
    operator==(const SeqList_iterator<_IteratorL, _Container>& __lhs,
               const SeqList_iterator<_IteratorR, _Container>& __rhs)
    { return __lhs.base() == __rhs.base(); }

  template<typename _Iterator, typename _Container>
    inline bool
    operator==(const SeqList_iterator<_Iterator, _Container>& __lhs,
               const SeqList_iterator<_Iterator, _Container>& __rhs)
    { return __lhs.base() == __rhs.base(); }

  template<typename _IteratorL, typename _IteratorR, typename _Container>
    inline bool
    operator!=(const SeqList_iterator<_IteratorL, _Container>& __lhs,
               const SeqList_iterator<_IteratorR, _Container>& __rhs)
    { return __lhs.base() != __rhs.base(); }

  template<typename _Iterator, typename _Container>
    inline bool
    operator!=(const SeqList_iterator<_Iterator, _Container>& __lhs,
               const SeqList_iterator<_Iterator, _Container>& __rhs)
    { return __lhs.base() != __rhs.base(); }

  // Random access iterator requirements
  template<typename _IteratorL, typename _IteratorR, typename _Container>
    inline bool
    operator<(const SeqList_iterator<_IteratorL, _Container>& __lhs,
              const SeqList_iterator<_IteratorR, _Container>& __rhs)
    { return __lhs.base() < __rhs.base(); }

  template<typename _Iterator, typename _Container>
    inline bool
    operator<(const SeqList_iterator<_Iterator, _Container>& __lhs,
              const SeqList_iterator<_Iterator, _Container>& __rhs)
    { return __lhs.base() < __rhs.base(); }

  template<typename _IteratorL, typename _IteratorR, typename _Container>
    inline bool
    operator>(const SeqList_iterator<_IteratorL, _Container>& __lhs,
              const SeqList_iterator<_IteratorR, _Container>& __rhs)
    { return __lhs.base() > __rhs.base(); }

  template<typename _Iterator, typename _Container>
    inline bool
    operator>(const SeqList_iterator<_Iterator, _Container>& __lhs,
              const SeqList_iterator<_Iterator, _Container>& __rhs)
    { return __lhs.base() > __rhs.base(); }

  template<typename _IteratorL, typename _IteratorR, typename _Container>
    inline bool
    operator<=(const SeqList_iterator<_IteratorL, _Container>& __lhs,
               const SeqList_iterator<_IteratorR, _Container>& __rhs)
    { return __lhs.base() <= __rhs.base(); }

  template<typename _Iterator, typename _Container>
    inline bool
    operator<=(const SeqList_iterator<_Iterator, _Container>& __lhs,
               const SeqList_iterator<_Iterator, _Container>& __rhs)
    { return __lhs.base() <= __rhs.base(); }

  template<typename _IteratorL, typename _IteratorR, typename _Container>
    inline bool
    operator>=(const SeqList_iterator<_IteratorL, _Container>& __lhs,
               const SeqList_iterator<_IteratorR, _Container>& __rhs)
    { return __lhs.base() >= __rhs.base(); }

  template<typename _Iterator, typename _Container>
    inline bool
    operator>=(const SeqList_iterator<_Iterator, _Container>& __lhs,
               const SeqList_iterator<_Iterator, _Container>& __rhs)
    { return __lhs.base() >= __rhs.base(); }

  // _GLIBCXX_RESOLVE_LIB_DEFECTS
  // According to the resolution of DR179 not only the various comparison
  // operators but also operator- must accept mixed iterator/const_iterator
  // parameters.
  template<typename _IteratorL, typename _IteratorR, typename _Container>
#ifdef __GXX_EXPERIMENTAL_CXX0X__
    // DR 685.
    inline auto
    operator-(const SeqList_iterator<_IteratorL, _Container>& __lhs,
              const SeqList_iterator<_IteratorR, _Container>& __rhs)
    -> decltype(__lhs.base() - __rhs.base())
#else
    inline typename SeqList_iterator<_IteratorL, _Container>::difference_type
    operator-(const SeqList_iterator<_IteratorL, _Container>& __lhs,
              const SeqList_iterator<_IteratorR, _Container>& __rhs)
#endif
    { return __lhs.base() - __rhs.base(); }

  template<typename _Iterator, typename _Container>
    inline typename SeqList_iterator<_Iterator, _Container>::difference_type
    operator-(const SeqList_iterator<_Iterator, _Container>& __lhs,
              const SeqList_iterator<_Iterator, _Container>& __rhs)
    { return __lhs.base() - __rhs.base(); }



#endif


                          
                          
                          
