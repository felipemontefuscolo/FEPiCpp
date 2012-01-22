#ifndef FEPIC_LIST_TYPE_HPP
#define FEPIC_LIST_TYPE_HPP

#include <vector>
#include <map>
#include <set>
#include <utility>
#include <iterator>
#if (FEP_HAS_OPENMP)
#  include "omp.h"
#endif
#include "../mesh/enums.hpp"
#include "../util/misc.hpp"
#include "contrib/Loki/set_vector.hpp"
#include <tr1/type_traits>
#include "../mesh/labelable.hpp"
#include <iostream>

/*
                        _      _  _       _
                       | |    | |(_) ___ | |_
                       | |    | || |/ __|| __|
                       | |___ | || |\__ \| |_
                       |_____||_||_||___/ \__|


*/

// fwd
template<class> class SeqList_iterator;
template<class> class SeqList_const_iterator;



/** @brief a mix between list and vector. The type T must be an inheritance
 *  of classe _Labelable.
 *
 *  @note Container_t, SmallCont_t : random access containers.
 *  @note dont insert disabled elements.
 */
template<class K,   // value type: class must be inherited of classe Labelable
         class C  = std::vector<K>, // a random access container
         class S  = SetVector<int> >   // a sorted container for indices
class SeqList
{

  template<class> friend class SeqList_iterator;
  template<class> friend class SeqList_const_iterator;


  typedef          SeqList<K,C,S>                        Self;
  typedef typename C::iterator                           ContainerIterator;
  typedef typename C::const_iterator                     ContainerConstIterator;
  typedef typename C::reverse_iterator                   ContainerReverseIterator;

public:
  typedef          K                                      value_type;
  typedef          C                                      container_type;
  typedef typename C::allocator_type                      allocator_type;
  typedef          S                                      ids_container_type;
  typedef typename container_type::reference              reference;
  typedef typename container_type::const_reference        const_reference;
  typedef typename container_type::pointer                pointer;
  typedef typename container_type::const_pointer          const_pointer;
  typedef typename container_type::size_type              size_type;
  typedef typename container_type::difference_type        difference_type;
  typedef          SeqList_iterator<Self>                 iterator;
  typedef          SeqList_const_iterator<Self>           const_iterator;

  // build TypeHas_reserve singnature checker.
  // mem_fun_name, mem_fun_return, qualif, mem_fun_args
  FEP_BUILD_MEM_FUN_CHECKER(reserve,reserve, void, , typename T::size_type);
  FEP_BUILD_MEM_FUN_CHECKER(capacity,capacity, typename T::size_type, const, );

  explicit SeqList(float grow_f=0.05) : _grow_factor(grow_f), _n_reserve_calls(0),
                                        _data(), _disabled_idcs(), _actived_beg(_data.begin())
                                        {}

  void clear()
  {
    _data.clear();
    _disabled_idcs.clear();
  }

  int getDataId(value_type const* v) const
  {
    return static_cast<int>(std::distance(_data.begin(), ContainerConstIterator(v)));
  }

  size_type size() const
  {return _data.size() - _disabled_idcs.size();};

  size_type total_size() const
  {return _data.size();};


  /*  SERIAL VERSION  */

  iterator begin()
  { return iterator(this, _actived_beg); }

  const_iterator begin() const
  { return const_iterator(this, _actived_beg); }  

  iterator end()
  { return iterator(this, _data.end()); }

  const_iterator end() const
  { return const_iterator(this, _data.end()); }



  /*  PARALLEL VERSION  */

  iterator begin(int tid, int nthreads)
  {
    int const N = size();
    int const start = tid*(N/nthreads) + ((N%nthreads) < tid ? (N%nthreads) : tid);
    iterator s = begin();
    for (int i = 0; i < start; ++i)
      ++s;
    return s;
  }

  iterator end(int tid, int nthreads)
  {
    int const N = size();
    int const start = tid*(N/nthreads) + ((N%nthreads) < tid ? (N%nthreads) : tid);
    int const end   = start + N/nthreads + ((N%nthreads) > tid);
    iterator e = begin(tid,nthreads);
    for (int i = 0; i < end-start; ++i)
      ++e;
    return e;
  }

  const_iterator begin(int tid, int nthreads) const
  {
    int const N = size();
    int const start = tid*(N/nthreads) + ((N%nthreads) < tid ? (N%nthreads) : tid);
    const_iterator s = begin();
    for (int i = 0; i < start; ++i)
      ++s;
    return s;
  }

  const_iterator end(int tid, int nthreads) const
  {
    int const N = size();
    int const start = tid*(N/nthreads) + ((N%nthreads) < tid ? (N%nthreads) : tid);
    int const end   = start + N/nthreads + ((N%nthreads) > tid);
    const_iterator e = begin(tid,nthreads);
    for (int i = 0; i < end-start; ++i)
      ++e;
    return e;
  }


  // modifiers
  void push_back(value_type const& obj)
  {

    if (_impl_capacity(_data) < _data.size()+1)
      _impl_reserve(_data, (size_type)((_data.size()+1)*(1. + _grow_factor))+0.5);

    _data.push_back(obj);

    _update_member_beg();
  };

  void disable(const_iterator it)
  { disable((int)std::distance((ContainerConstIterator)_data.begin(), it.get()));};

  void disable(int del_id)
  {
    ContainerIterator it = ContainerIterator(&_data[del_id]);
    if (it->disabled())
      return;
    _disabled_idcs.insert(del_id);
    it->disabled(true);

    if (it == _actived_beg)
      _update_member_beg();
  }

  /** @brief insert an element and return your id.
   *  @warning dont insert disabled elements.
   */
  int insert(value_type const& obj)
  {
    if (_disabled_idcs.empty())
    {
      push_back(obj);
      return _data.size()-1;
    }

    int const new_id = _disabled_idcs.back();
    _disabled_idcs.pop_back();
    _data[new_id] = obj;

    _update_member_beg();
    return new_id;

  }

  void setGrowFactor(float factor)
  { _grow_factor = factor; }

  float getGrowFactor() const
  { return _grow_factor; }

  void reserve(size_type amount)
  {
    _impl_reserve(_data, amount);
  }

  void resize(size_type s)
  {
    if (s > _data.size())
      _impl_reserve(_data, (size_type)s*(1. + _grow_factor));
    _data.resize(s);
    _update_member_beg();
  }

  size_type capacity() const
  {
    return _impl_capacity(_data);
  }

  reference operator[](size_type n)
  {
    return _data[n];
  }

  const_reference operator[](size_type n) const
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
  
  /** @param[out] cids contiguous ids.
   */  
  void contiguousIds(int *ids_beg, int const* ids_end, int *cids) const
  {
    while (ids_beg != ids_end)
    {
      *cids++ = contiguousId(*ids_beg++);
    }
  }
    
protected:

  template<class V>
  void _impl_reserve(V &vec, typename V::size_type amount,
              typename EnableIf<!TypeHas_reserve<V>::value>::type * = NULL) {}

  template<class V>
  void _impl_reserve(V &vec, typename V::size_type amount,
              typename EnableIf<TypeHas_reserve<V>::value>::type * = NULL)
  {
    if ( _impl_capacity(vec) < amount )
    {
      ++_n_reserve_calls;
      _update_member_beg();
    }
    vec.reserve(amount);
  }

  template<class V>
  static
  typename V::size_type _impl_capacity(V const&vec,
                        typename EnableIf<!TypeHas_capacity<V>::value>::type * = NULL)
  { return 0; }

  template<class V>
  static
  typename V::size_type _impl_capacity(V const&vec,
                        typename EnableIf<TypeHas_capacity<V>::value>::type * = NULL)
  {
    return vec.capacity();
  }

  void _update_member_beg()
  {
    _actived_beg = _data.begin();
    while (_actived_beg != _data.end() && _actived_beg->disabled())
      ++_actived_beg;
  }

  float               _grow_factor;
  unsigned            _n_reserve_calls;
  container_type      _data;
  ids_container_type  _disabled_idcs; // sorted vector
  ContainerIterator   _actived_beg; // beginning of valid data

};




// fwd
template<class SeqListType>
class SeqList_iterator
{

  template<class,class,class> friend class SeqList;
  template<class> friend class SeqList_const_iterator;

  typedef SeqListType*                            PtrToSeqListType;
  typedef SeqList_iterator<SeqListType>           Self;
  typedef typename SeqListType::ContainerIterator ContainerIterator;
public:
  typedef typename SeqListType::difference_type    difference_type;
  typedef typename SeqListType::value_type         value_type;
  typedef typename SeqListType::pointer            pointer;
  typedef typename SeqListType::reference          reference;
  typedef          std::bidirectional_iterator_tag iterator_category;

  explicit
  SeqList_iterator(SeqListType * sq, ContainerIterator x) : _iter_to_t(x), _ptr_to_seq(sq) {}
  explicit
  SeqList_iterator(SeqListType * sq, pointer p) : _iter_to_t(ContainerIterator(p)), _ptr_to_seq(sq) {}
  explicit
  SeqList_iterator(SeqListType * sq) : _iter_to_t(), _ptr_to_seq(sq) {}

  SeqList_iterator(Self const& it) : _iter_to_t(it._iter_to_t), _ptr_to_seq(it._ptr_to_seq) {}

  SeqList_iterator() : _iter_to_t(), _ptr_to_seq(NULL) {}

  ContainerIterator const& get() const
  { return _iter_to_t; }

  Self& operator=(Self const& foo)
  {
    _iter_to_t = foo._iter_to_t;
    _ptr_to_seq = foo._ptr_to_seq;
    return *this;
  }

  reference
  operator*() const
  {
    return *_iter_to_t;
  }

  pointer
  operator->() const
  {
    return &(*_iter_to_t);
  }

  Self&
  operator++()
  {
    ++_iter_to_t;
    while(_iter_to_t != _ptr_to_seq->_data.end() && _iter_to_t->disabled())
      ++_iter_to_t;
    return *this;
  }

  Self
  operator++(int)
  {
    Self tmp = *this;
    ++_iter_to_t;
    while(_iter_to_t != _ptr_to_seq->_data.end() && _iter_to_t->disabled())
      ++_iter_to_t;
    return tmp;
  }

  Self&
  operator--()
  {
    --_iter_to_t;
    while(_iter_to_t != _ptr_to_seq->_data.begin() && _iter_to_t->disabled())
      --_iter_to_t;
    return *this;
  }

  Self
  operator--(int)
  {
    Self tmp = *this;
    --_iter_to_t;
    while(_iter_to_t != _ptr_to_seq->_data.begin() && _iter_to_t->disabled())
      --_iter_to_t;
    return tmp;
  }

  bool
  operator==(const Self& x) const
  { return _iter_to_t == x._iter_to_t; }

  bool
  operator!=(const Self& x) const
  { return _iter_to_t != x._iter_to_t; }

private:
  ContainerIterator _iter_to_t;
  PtrToSeqListType  _ptr_to_seq;

};




// fwd
template<class SeqListType>
class SeqList_const_iterator
{
  template<class,class,class> friend class SeqList;
  template<class> friend class SeqList_iterator;

  typedef SeqListType const*                           PtrToSeqListType;
  typedef SeqList_const_iterator<SeqListType>          Self;
  typedef SeqList_iterator<SeqListType>                SelfNoConst;
  typedef typename SeqListType::ContainerConstIterator ContainerIterator;

public:
  typedef typename SeqListType::difference_type    difference_type;
  typedef typename SeqListType::value_type         value_type;
  typedef typename SeqListType::const_pointer      pointer;
  typedef typename SeqListType::const_reference    reference;
  typedef          std::bidirectional_iterator_tag iterator_category;

  SeqList_const_iterator(SeqListType const* sq, ContainerIterator const& x) : _iter_to_t(x), _ptr_to_seq(sq) {}
  explicit
  SeqList_const_iterator(SeqListType const* sq, pointer p) : _iter_to_t(ContainerIterator(p)), _ptr_to_seq(sq) {}
  explicit
  SeqList_const_iterator(SeqListType const* sq) : _iter_to_t(), _ptr_to_seq(sq) {}

  SeqList_const_iterator(Self const& it)        : _iter_to_t(it._iter_to_t), _ptr_to_seq(it._ptr_to_seq) {}
  SeqList_const_iterator(SelfNoConst const& it) : _iter_to_t(it._iter_to_t), _ptr_to_seq(it._ptr_to_seq) {}

  SeqList_const_iterator() : _iter_to_t(), _ptr_to_seq(NULL) {}

  ContainerIterator const& get() const
  { return _iter_to_t; }

  Self& operator=(Self const& foo)
  {
    _iter_to_t = foo._iter_to_t;
    _ptr_to_seq = foo._ptr_to_seq;
    return *this;
  }

  reference
  operator*() const
  {
    return *_iter_to_t;
  }

  pointer
  operator->() const
  {
    return &(*_iter_to_t);
  }

  Self&
  operator++()
  {
    ++_iter_to_t;
    while(_iter_to_t != _ptr_to_seq->_data.end() && _iter_to_t->disabled())
      ++_iter_to_t;
    return *this;
  }

  Self
  operator++(int)
  {
    Self tmp = *this;
    ++_iter_to_t;
    while(_iter_to_t != _ptr_to_seq->_data.end() && _iter_to_t->disabled())
      ++_iter_to_t;
    return tmp;
  }

  Self&
  operator--()
  {
    --_iter_to_t;
    while(_iter_to_t != _ptr_to_seq->_data.begin() && _iter_to_t->disabled())
      --_iter_to_t;
    return *this;
  }

  Self
  operator--(int)
  {
    Self tmp = *this;
    --_iter_to_t;
    while(_iter_to_t != _ptr_to_seq->_data.begin() && _iter_to_t->disabled())
      --_iter_to_t;
    return tmp;
  }

  bool
  operator==(const Self& x) const
  { return _iter_to_t == x._iter_to_t; }

  bool
  operator!=(const Self& x) const
  { return _iter_to_t != x._iter_to_t; }

private:
  ContainerIterator _iter_to_t;
  PtrToSeqListType  _ptr_to_seq;

};



#endif
