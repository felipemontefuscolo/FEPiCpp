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
template<class> class SeqList_iterator;
template<class> class SeqList_const_iterator;



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
  template<class> friend class SeqList_iterator;
  template<class> friend class SeqList_const_iterator;

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
  typedef          SeqList_iterator<Self>                 iterator;
  typedef          SeqList_const_iterator<Self>           const_iterator;

private:
  
  // auxiliary typedefs
  typedef typename boost::remove_pointer<value_type>::type RP_ValueType; // remove pointer value type
  typedef          RP_ValueType &                          RP_Reference;
  typedef          RP_ValueType const&                     RP_ConstReference;
  typedef          RP_ValueType *                          RP_Pointer;
  typedef          RP_ValueType const*                     RP_ConstPointer;

public:

  // build TypeHas_reserve singnature checker.
  // suffix,mem_fun_name, mem_fun_return, qualif, mem_fun_args
  //FEP_BUILD_MEM_FUN_CHECKER(reserve,reserve, void, /*empty*/ , typename T::size_type);
  //FEP_BUILD_MEM_FUN_CHECKER(capacity,capacity, typename T::size_type, const, /*empty*/);
  FEP_BUILD_MEM_FUN_CHECKER_1ARG(reserve,reserve, void, ;char none_qualif , typename T::size_type);
  FEP_BUILD_MEM_FUN_CHECKER_0ARG(capacity,capacity, typename T::size_type, const);

  template<class T> // T = value_type
  struct _InsertFunArgument {
    typedef T const& type;
  };
  
  template<class T> // T* = value_type
  struct _InsertFunArgument<T*> {
    typedef T*  type;
  };


  explicit SeqList(float grow_f=0.05f) : _grow_factor(grow_f), _n_reserve_calls(0),
                                         _data(), _disabled_idcs(), _actived_beg(_data.begin())
                                        {}

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
  { disable((int)std::distance((DataConstIterator)_data.begin(), it.get()));};

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

  void setGrowFactor(float factor)
  { _grow_factor = factor; }

  float getGrowFactor() const
  { return _grow_factor; }

  void reserve(size_type amount)
  {
    reserve_impl(_data, amount);
  }

  void resize(size_type s)
  {
    if (s > _data.size())
      reserve_impl(_data, static_cast<size_type>( static_cast<float>(s)*(1.0f + _grow_factor) ) );
    _data.resize(s);
    _update_member_beg();
  }

  size_type capacity() const
  {
    return capacity_impl(_data);
  }

  RP_Reference operator[](size_type n)
  {
    return _data[n];
  }

  RP_ConstReference operator[](size_type n) const
  {
    return _data[n];
  }

  int numReserveCalls() const
  {
    return _n_reserve_calls;
  }
  
  int contiguousId(int id) const
  {
    return id - static_cast<int>( std::distance(_disabled_idcs.begin(), std::lower_bound(_disabled_idcs.begin(), _disabled_idcs.end(), id)));
  }
  
  /// A shortcut to get multiples cids at same time.
  /// @param[out] cids contiguous ids.
  ///  
  void contiguousIds(int *ids_beg, int const* ids_end, int *cids) const
  {
    while (ids_beg != ids_end)
    {
      *cids++ = contiguousId(*ids_beg++);
    }
  }
    
protected:

  template<class V>
  void reserve_impl(V &/*vec*/, typename V::size_type /*amount*/,
              typename boost::enable_if_c<!TypeHas_reserve<V>::value>::type * = NULL) {}

  template<class V>
  void reserve_impl(V &vec, typename V::size_type amount,
              typename boost::enable_if_c<TypeHas_reserve<V>::value>::type * = NULL)
  {
    if ( capacity_impl(vec) < amount )
    {
      ++_n_reserve_calls;
      _update_member_beg();
    }
    vec.reserve(amount);
  }

  template<class V>
  static
  typename V::size_type capacity_impl(V const& /*vec*/,
                        typename boost::enable_if_c<!TypeHas_capacity<V>::value>::type * = NULL)
  { return ~(size_type)0; }

  template<class V>
  static
  typename V::size_type capacity_impl(V const&vec,
                        typename boost::enable_if_c<TypeHas_capacity<V>::value>::type * = NULL)
  {
    return vec.capacity();
  }

  template<class Value_type>
  int insert_impl(const_reference obj, typename boost::enable_if_c< ! std::tr1::is_pointer<Value_type>::value >::type * = NULL)
  {
    if (_disabled_idcs.empty())
    {
      // --- push_back ----
      if (capacity_impl(_data) < _data.size()+1)
        reserve_impl(_data, (size_type)((_data.size()+1)*(1. + _grow_factor))+0.5);
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
  int insert_impl(value_type obj, typename boost::enable_if_c< std::tr1::is_pointer<Value_type>::value >::type * = NULL)
  {
    if (_disabled_idcs.empty())
    {
      // --- push_back ----
      if (capacity_impl(_data) < _data.size()+1)
        reserve_impl(_data, (size_type)((_data.size()+1)*(1. + _grow_factor))+0.5);
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

  float               _grow_factor;
  unsigned            _n_reserve_calls;
  container_type      _data;
  ids_container_type  _disabled_idcs; // sorted vector
  DataIterator        _actived_beg;   // iterator to the beginning of valid data

};




// fwd
template<class SeqListType>
class SeqList_iterator
{

  template<class,class> friend class SeqList;
  template<class> friend class SeqList_const_iterator;

  typedef SeqListType*                            PtrToSeqListType;
  typedef SeqList_iterator<SeqListType>           Self;
  typedef typename SeqListType::DataIterator DataIterator;
public:

  // stl standard
  typedef typename SeqListType::difference_type    difference_type;
  typedef typename SeqListType::value_type         value_type;
  typedef typename SeqListType::pointer            pointer;
  typedef typename SeqListType::reference          reference;
  typedef          std::bidirectional_iterator_tag iterator_category;

private:
  
  // auxiliary typedefs
  typedef typename SeqListType::RP_ValueType      RP_ValueType;
  typedef typename SeqListType::RP_Reference      RP_Reference;
  typedef typename SeqListType::RP_ConstReference RP_ConstReference;
  typedef typename SeqListType::RP_Pointer        RP_Pointer;


public:

  explicit
  SeqList_iterator(SeqListType * sq, DataIterator x) : _iter_to_t(x), _ptr_to_seq(sq) {}
  //explicit
  //SeqList_iterator(SeqListType * sq, pointer p) : _iter_to_t(DataIterator(p)), _ptr_to_seq(sq) {}
  explicit
  SeqList_iterator(SeqListType * sq) : _iter_to_t(), _ptr_to_seq(sq) {}

  SeqList_iterator(Self const& it) : _iter_to_t(it._iter_to_t), _ptr_to_seq(it._ptr_to_seq) {}

  SeqList_iterator() : _iter_to_t(), _ptr_to_seq(NULL) {}

  DataIterator const& get() const
  { return _iter_to_t; }

  FEP_STRONG_INLINE
  Self& operator=(Self const& foo)
  {
    _iter_to_t = foo._iter_to_t;
    _ptr_to_seq = foo._ptr_to_seq;
    return *this;
  }

  FEP_STRONG_INLINE
  RP_Reference
  operator*() const
  {
    return *_iter_to_t;
  }

  FEP_STRONG_INLINE
  RP_Pointer
  operator->() const
  {
    return &(*_iter_to_t);
  }

  
  FEP_STRONG_INLINE
  Self&
  operator++()
  {
    ++_iter_to_t;
    while(_iter_to_t != _ptr_to_seq->_data.end() && _iter_to_t->isDisabled())
      ++_iter_to_t;
    return *this;
  }

  FEP_STRONG_INLINE
  Self
  operator++(int)
  {
    Self tmp = *this;
    ++_iter_to_t;
    while(_iter_to_t != _ptr_to_seq->_data.end() && _iter_to_t->isDisabled())
      ++_iter_to_t;
    return tmp;
  }

  FEP_STRONG_INLINE
  Self&
  operator--()
  {
    --_iter_to_t;
    while(_iter_to_t != _ptr_to_seq->_data.begin() && _iter_to_t->isDisabled())
      --_iter_to_t;
    return *this;
  }

  FEP_STRONG_INLINE
  Self
  operator--(int)
  {
    Self tmp = *this;
    --_iter_to_t;
    while(_iter_to_t != _ptr_to_seq->_data.begin() && _iter_to_t->isDisabled())
      --_iter_to_t;
    return tmp;
  }

  FEP_STRONG_INLINE
  bool
  operator==(const Self& x) const
  { return _iter_to_t == x._iter_to_t; }

  FEP_STRONG_INLINE
  bool
  operator!=(const Self& x) const
  { return _iter_to_t != x._iter_to_t; }

  FEP_STRONG_INLINE
  Self&
  next()
  {
    ++_iter_to_t;
    return *this;
  }

  FEP_STRONG_INLINE
  Self&
  previous()
  {
    --_iter_to_t;
    return *this;
  }

  FEP_STRONG_INLINE
  Self
  plus(difference_type const& n)
  {
    Self tmp = Self(_ptr_to_seq, _iter_to_t + n);
    return tmp;
  }

private:
  DataIterator _iter_to_t;
  PtrToSeqListType  _ptr_to_seq;

};




// fwd
template<class SeqListType>
class SeqList_const_iterator
{
  template<class,class> friend class SeqList;
  template<class> friend class SeqList_iterator;

  typedef SeqListType const*                           PtrToSeqListType;
  typedef SeqList_const_iterator<SeqListType>          Self;
  typedef SeqList_iterator<SeqListType>                SelfNoConst;
  typedef typename SeqListType::DataConstIterator DataIterator;

public:
  typedef typename SeqListType::difference_type    difference_type;
  typedef typename SeqListType::value_type         value_type;
  //typedef typename SeqListType::const_pointer      pointer;
  typedef typename SeqListType::const_reference    reference;
  typedef          std::bidirectional_iterator_tag iterator_category;

private:
  
  // auxiliary typedefs
  typedef typename SeqListType::RP_ValueType      RP_ValueType;
  typedef typename SeqListType::RP_Reference      RP_Reference;
  typedef typename SeqListType::RP_ConstReference RP_ConstReference;
  typedef typename SeqListType::RP_Pointer        RP_Pointer;
  typedef typename SeqListType::RP_ConstPointer   RP_ConstPointer;

public:


  SeqList_const_iterator(SeqListType const* sq, DataIterator const& x) : _iter_to_t(x), _ptr_to_seq(sq) {}
  explicit
  //SeqList_const_iterator(SeqListType const* sq, pointer p) : _iter_to_t(DataIterator(p)), _ptr_to_seq(sq) {}
  //explicit
  SeqList_const_iterator(SeqListType const* sq) : _iter_to_t(), _ptr_to_seq(sq) {}

  SeqList_const_iterator(Self const& it)        : _iter_to_t(it._iter_to_t), _ptr_to_seq(it._ptr_to_seq) {}
  SeqList_const_iterator(SelfNoConst const& it) : _iter_to_t(it._iter_to_t), _ptr_to_seq(it._ptr_to_seq) {}

  SeqList_const_iterator() : _iter_to_t(), _ptr_to_seq(NULL) {}

  FEP_STRONG_INLINE
  DataIterator const& get() const
  { return _iter_to_t; }

  FEP_STRONG_INLINE
  Self& operator=(Self const& foo)
  {
    _iter_to_t = foo._iter_to_t;
    _ptr_to_seq = foo._ptr_to_seq;
    return *this;
  }

  FEP_STRONG_INLINE
  RP_Reference
  operator*() const
  {
    return *_iter_to_t;
  }

  FEP_STRONG_INLINE
  RP_Pointer
  operator->() const
  {
    return &(*_iter_to_t);
  }

  FEP_STRONG_INLINE
  Self&
  operator++()
  {
    ++_iter_to_t;
    while(_iter_to_t != _ptr_to_seq->_data.end() && _iter_to_t->isDisabled())
      ++_iter_to_t;
    return *this;
  }

  FEP_STRONG_INLINE
  Self
  operator++(int)
  {
    Self tmp = *this;
    ++_iter_to_t;
    while(_iter_to_t != _ptr_to_seq->_data.end() && _iter_to_t->isDisabled())
      ++_iter_to_t;
    return tmp;
  }

  FEP_STRONG_INLINE
  Self&
  operator--()
  {
    --_iter_to_t;
    while(_iter_to_t != _ptr_to_seq->_data.begin() && _iter_to_t->isDisabled())
      --_iter_to_t;
    return *this;
  }

  FEP_STRONG_INLINE
  Self
  operator--(int)
  {
    Self tmp = *this;
    --_iter_to_t;
    while(_iter_to_t != _ptr_to_seq->_data.begin() && _iter_to_t->isDisabled())
      --_iter_to_t;
    return tmp;
  }

  FEP_STRONG_INLINE
  bool
  operator==(const Self& x) const
  { return _iter_to_t == x._iter_to_t; }

  FEP_STRONG_INLINE
  bool
  operator!=(const Self& x) const
  { return _iter_to_t != x._iter_to_t; }

  FEP_STRONG_INLINE
  Self&
  next()
  {
    ++_iter_to_t;
    return *this;
  }

  FEP_STRONG_INLINE
  Self&
  previous()
  {
    --_iter_to_t;
    return *this;
  }

  FEP_STRONG_INLINE
  Self
  plus(difference_type const& n)
  {
    Self tmp = Self(_ptr_to_seq, _iter_to_t + n);
    return tmp;
  }

private:
  DataIterator _iter_to_t;
  PtrToSeqListType  _ptr_to_seq;

};



#endif
