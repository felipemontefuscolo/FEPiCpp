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
#include "../mesh/colored.hpp"
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
template<class> class SeqList_color_iterator;
template<class> class SeqList_color_const_iterator;



/** @brief a mix between list and vector. The type T must be an inheritance
 *  of classes _Colored and _Labelable.
 *
 *  @note Container_t, SmallCont_t : random access containers.
 *  @note dont insert disabled elements.
 */
template<class K,   // value type: class must be inherited of classes Colored and Labelable
         class C  = std::vector<K>, // a random access container
         class S  = SetVector<int> >   // a sorted container for indices
class SeqList
{

  template<class> friend class SeqList_iterator;
  template<class> friend class SeqList_const_iterator;
  template<class> friend class SeqList_color_iterator;
  template<class> friend class SeqList_color_const_iterator;

  struct _FirstNode
  {
    _FirstNode(int a=-1, int b=0) : head(a), size(b) {}
    int head;
    int size;
  };

  typedef          SeqList<K,C,S>                        Self;
  typedef          std::map<EColor, _FirstNode>          HeadsContainer; /* <color,(head, size)> */
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
  typedef          SeqList_color_iterator<Self>           color_iterator;
  typedef          SeqList_color_const_iterator<Self>     color_const_iterator;

  // build TypeHas_reserve singnature checker.
  // mem_fun_name, mem_fun_return, qualif, mem_fun_args
  FEP_BUILD_MEM_FUN_CHECKER(reserve,reserve, void, , typename T::size_type);
  FEP_BUILD_MEM_FUN_CHECKER(capacity,capacity, typename T::size_type, const, );

  explicit SeqList(float grow_f=0.05) : _grow_factor(grow_f), _n_reserve_calls(0),
                                        _colors_heads(), _data(), _disabled_idcs(), _actived_beg(_data.begin())
                                        {}

  void clear()
  {
    _data.clear();
    _colors_heads.clear();
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

  int numColors() const
  {
    int c = 0;
    typename HeadsContainer::const_iterator it = _colors_heads.begin();
    for (; it != _colors_heads.end(); ++it)
    {
      c += (*it).second.size > 0;
    }
    return c;
  }

  int maxColor() const
  {
    if (!_colors_heads.empty())
    {
      return (*_colors_heads.rbegin()).first;
    }
    else
      return -1;
  }

  int colorSize(EColor col) const
  {
    typename HeadsContainer::const_iterator it = _colors_heads.find(col);
    
    if (it == _colors_heads.end())
      return 0;
    return (*it).second.size;
  }

  /*  SERIAL VERSION  */

  iterator begin()
  { return iterator(this, _actived_beg); }

  const_iterator begin() const
  { return const_iterator(this, _actived_beg); }  

  iterator end()
  { return iterator(this, _data.end()); }

  const_iterator end() const
  { return const_iterator(this, _data.end()); }

  color_iterator begin(EColor col)
  {
    typename HeadsContainer::iterator i = _colors_heads.find(col);
    if ( i == _colors_heads.end())
      return color_iterator(this, _data.end());
    return color_iterator(this, &_data[(*i).second.head]);
  }

  color_const_iterator begin(EColor col) const
  {
    typename HeadsContainer::const_iterator i = _colors_heads.find(col);
    if ( i == _colors_heads.end())
    {
      return color_const_iterator(this, _data.end());
    }
    return color_const_iterator(this, &_data[(*i).second.head]);
  }
  

  color_iterator end(EColor)
  { return color_iterator(this, _data.end()); }

  color_const_iterator end(EColor) const
  { return color_const_iterator(this, _data.end()); }




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

  // =======================================================================

  color_iterator begin(EColor col, int tid, int nthreads)
  {
    typename HeadsContainer::iterator it = _colors_heads.find(col);
    if ( it == _colors_heads.end())
      return color_iterator(this, _data.end());
    int const N     = (*it).second.size;
    int const start = tid*(N/nthreads) + ((N%nthreads) < tid ? (N%nthreads) : tid);
    color_iterator s = begin(col);
    for (int i = 0; i < start; ++i)
      ++s;
    return s;
  }

  // parallel version
  color_const_iterator begin(EColor col, int tid, int nthreads) const
  {
    typename HeadsContainer::iterator it = _colors_heads.find(col);
    if ( it == _colors_heads.end())
      return color_iterator(this, _data.end());
    int const N     = (*it).second.size;
    int const start = tid*(N/nthreads) + ((N%nthreads) < tid ? (N%nthreads) : tid);
    color_const_iterator s = begin(col);
    for (int i = 0; i < start; ++i)
      ++s;
    return s;
  }


  // parallel version
  color_iterator end(EColor col, int tid, int nthreads)
  {
    typename HeadsContainer::iterator it = _colors_heads.find(col);
    if ( it == _colors_heads.end())
      return color_iterator(this, _data.end());
    int const N     = (*it).second.size;
    int const start = tid*(N/nthreads) + ((N%nthreads) < tid ? (N%nthreads) : tid);
    int const end   = start + N/nthreads + ((N%nthreads) > tid);
    color_iterator e = begin(col,tid,nthreads);
    for (int i = 0; i < end-start; ++i)
      ++e;
    return e;
  }

  // parallel version
  color_const_iterator end(EColor col, int tid, int nthreads) const
  {
    typename HeadsContainer::iterator it = _colors_heads.find(col);
    if ( it == _colors_heads.end())
      return color_iterator(this, _data.end());
    int const N     = (*it).second.size;
    int const start = tid*(N/nthreads) + ((N%nthreads) < tid ? (N%nthreads) : tid);
    int const end   = start + N/nthreads + ((N%nthreads) > tid);
    color_const_iterator e = begin(col,tid,nthreads);
    for (int i = 0; i < end-start; ++i)
      ++e;
    return e;
  }



  // modifiers
  void push_back(value_type const& obj)
  {

    if (_impl_capacity(_data) < _data.size()+1)
      _impl_reserve(_data, (size_type)((_data.size()+1)*(1. + _grow_factor))+0.5);

    char const c = obj.getColor();
    if (obj.disabled() || c<0)
    {
      _data.push_back(obj);
    }
    else
    {
      const int new_id = _data.size();
      std::pair<typename HeadsContainer::iterator, bool> result
                        = _colors_heads.insert(std::make_pair(EColor(c), _FirstNode(new_id, 1)));

      // if already exists
      if (!result.second) {
        typename HeadsContainer::iterator h = result.first;
        int const old_head = (*h).second.head;
        _data[old_head].setPrev(new_id);
        _data.push_back(obj);
        _data.back().setNext(old_head);
        _data.back().setPrev(-1);
        (*h).second.head = new_id;
        ++(*h).second.size;
      }
      else
      {
        _data.push_back(obj);
        _data.back().setNext(-1);
        _data.back().setPrev(-1);
      }
    }

    _update_member_beg();
  };

  void disable(const_iterator it)
  { disable((int)std::distance((ContainerConstIterator)_data.begin(), it.get()));};

  void disable(color_const_iterator it)
  { disable((int)std::distance((ContainerConstIterator)_data.begin(), it.get()));};

  void disable(int del_id)
  {
    ContainerIterator it = ContainerIterator(&_data[del_id]);
    if (it->disabled())
      return;
    _disabled_idcs.insert(del_id);
    const char c = it->getColor();
    if (c<0)
    {
      it->disabled(true);
      if (it == _actived_beg)
        _update_member_beg();
      return;
    }
    it->disabled(true);

    typename HeadsContainer::iterator h = _colors_heads.find(EColor(c));

    int next = it->sameColorNext();
    int prev = it->sameColorPrev();
    if (prev >= 0)
    {
      _data[prev].setNext(next);
      if (next >=0)
        _data[next].setPrev(prev);
      --(*h).second.size;
    }
    else
    if (next >=0)
    {
      _data[next].setPrev(prev);
      (*h).second.head = next;
      --(*h).second.size;
    }
    else
      _colors_heads.erase(h);


    if (it == _actived_beg)
      _update_member_beg();
  }

  /** @brief insert an element and return your id.
   *  @warning dont insert disabled elements.
   *  @note automatically inserted in the color list.
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
    const char c = obj.getColor();

    if (c>=0 && !obj.disabled())
    {
      std::pair<typename HeadsContainer::iterator, bool> result
                        = _colors_heads.insert(std::make_pair(EColor(c), _FirstNode(new_id, 1)));
      // if already exists
      if (!result.second) {
        typename HeadsContainer::iterator h = result.first;
        int const old_head = (*h).second.head;
        _data[old_head].setPrev(new_id);
        _data[new_id].setNext(old_head);
        _data[new_id].setPrev(-1);
        (*h).second.head = new_id;
        ++(*h).second.size;
      }
      else
      {
        _data.back().setNext(-1);
        _data.back().setPrev(-1);
      }
      
    }

    _update_member_beg();
    return new_id;

  }

  void setGrowFactor(float factor)
  { _grow_factor = factor; }

  float getGrowFactor() const
  { return _grow_factor; }

  // links elements of color k and counts how many elements has this color
  // not thread safety in general, but you can call _linkColor(0) and _linkColor(1) for example
  void _linkColor(int const k)
  {
    unsigned color_counter = 0;
    ContainerReverseIterator next_it; // reverse iterator
    int next_id(-1);
    int current_id = _data.size()-1;
    for (ContainerReverseIterator it = _data.rbegin(); it != _data.rend(); ++it, --current_id)
    {
      if (it->value_type::getColor() != k || it->value_type::disabled())
        continue;
      it->value_type::setNext(next_id);
      if (next_id >= 0)
        _data[next_id].setPrev(current_id);
      it->value_type::setPrev(-1); // in case of no previous data
      next_id = current_id;
      ++color_counter;
    }
    typename HeadsContainer::iterator h = _colors_heads.find(EColor(k));
    (*h).second.head  = next_id;
    (*h).second.size  = color_counter;
  }


  /** @brief creates a list for each color;
   */
  void linkColors()
  {
    _colors_heads.clear();
    int max_color = -1;

    // max reduction to search max color
    #pragma omp parallel shared(max_color)
    {
      int max_local_color = -1;

      #pragma omp for
      for (ContainerIterator it = _data.begin(); it < _data.end(); ++it)
        if (max_local_color < it->value_type::getColor())
          max_local_color = it->value_type::getColor();

      #pragma omp critical
      if (max_color < max_local_color)
        max_color = max_local_color;

    } // end parallel

    int const num_colors = max_color + 1;

    for (int k = 0; k < num_colors; ++k)
      _colors_heads.insert( std::make_pair(EColor(k), _FirstNode()) );

    #pragma omp parallel default(none)
    {
      // make link
      #pragma omp for
      for (int k = 0; k < num_colors; ++k)
        _linkColor(k);
        
    } // end parallel
    
    
    // deleting
    // all this work is to avoid gcc warning
    std::vector<int> to_del;
    to_del.reserve(_colors_heads.size());
    typename HeadsContainer::iterator h = _colors_heads.begin();
    for(;h != _colors_heads.end(); ++h)
    {
      if ((*h).second.head==-1)
      {
        to_del.push_back((*h).first);
      }
    }
    
    for (unsigned i = 0; i < to_del.size(); ++i)
      _colors_heads.erase(EColor(to_del[i]));
  }

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
  HeadsContainer      _colors_heads;
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





// fwd
template<class SeqListType>
class SeqList_color_iterator
{
  template<class,class,class> friend class SeqList;
  template<class> friend class SeqList_color_const_iterator;

  typedef SeqListType*                            PtrToSeqListType;
  typedef SeqList_color_iterator<SeqListType>     Self;
  typedef typename SeqListType::ContainerIterator ContainerIterator;
public:
  typedef typename SeqListType::difference_type    difference_type;
  typedef typename SeqListType::value_type         value_type;
  typedef typename SeqListType::pointer            pointer;
  typedef typename SeqListType::reference          reference;
  typedef          std::forward_iterator_tag       iterator_category;

  explicit SeqList_color_iterator(SeqListType * sq, ContainerIterator x)
                                  : _iter_to_t(x), _ptr_to_seq(sq) {}

  explicit SeqList_color_iterator(SeqListType * sq, pointer p)
                                  : _iter_to_t(ContainerIterator(p)), _ptr_to_seq(sq) {}

  explicit SeqList_color_iterator(SeqListType * sq) : _iter_to_t(), _ptr_to_seq(sq) {}

  SeqList_color_iterator(Self const& it)
                         : _iter_to_t(it._iter_to_t), _ptr_to_seq(it._ptr_to_seq) {}

  SeqList_color_iterator() : _iter_to_t(), _ptr_to_seq(NULL) {}

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
    const int next = _iter_to_t->sameColorNext();
    if (next < 0)
    {
      _iter_to_t = _ptr_to_seq->_data.end();
      return *this;
    }
    _iter_to_t = ContainerIterator( &(_ptr_to_seq->_data[next]));
    return *this;
  }

  Self
  operator++(int)
  {
    Self tmp = *this;
    const int next = _iter_to_t->sameColorNext();
    if (next < 0)
    {
      _iter_to_t = _ptr_to_seq->_data.end();
      return tmp;
    }
    _iter_to_t = ContainerIterator( &(_ptr_to_seq->_data[next]));
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
class SeqList_color_const_iterator
{
  template<class,class,class> friend class SeqList;
  template<class> friend class SeqList_color_iterator;

  typedef SeqListType const*                        PtrToSeqListType;
  typedef SeqList_color_const_iterator<SeqListType> Self;
  typedef SeqList_color_iterator<SeqListType>       SelfNoConst;
  typedef typename SeqListType::ContainerConstIterator   ContainerIterator;
public:
  typedef typename SeqListType::difference_type    difference_type;
  typedef typename SeqListType::value_type         value_type;
  typedef typename SeqListType::const_pointer      pointer;
  typedef typename SeqListType::const_reference    reference;
  typedef          std::forward_iterator_tag       iterator_category;

  explicit SeqList_color_const_iterator(SeqListType const* sq, ContainerIterator x)
                                  : _iter_to_t(x), _ptr_to_seq(sq) {}

  explicit SeqList_color_const_iterator(SeqListType const* sq, pointer p)
                                  : _iter_to_t(ContainerIterator(p)), _ptr_to_seq(sq) {}

  explicit SeqList_color_const_iterator(SeqListType const* sq) : _iter_to_t(), _ptr_to_seq(sq) {}

  SeqList_color_const_iterator(Self const& it)
                         : _iter_to_t(it._iter_to_t), _ptr_to_seq(it._ptr_to_seq) {}

  SeqList_color_const_iterator(SelfNoConst const& it)
                         : _iter_to_t(it._iter_to_t), _ptr_to_seq(it._ptr_to_seq) {}
  
  SeqList_color_const_iterator() : _iter_to_t(), _ptr_to_seq(NULL) {}

  ContainerIterator const& get() const
  { return _iter_to_t; }

  Self& operator=(Self const& foo)
  {
    _iter_to_t  = foo._iter_to_t;
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
    const int next = _iter_to_t->sameColorNext();
    if (next < 0)
    {
      _iter_to_t = _ptr_to_seq->_data.end();
      return *this;
    }
    _iter_to_t = ContainerIterator( &(_ptr_to_seq->_data.at(next)));
    return *this;
  }

  Self
  operator++(int)
  {
    Self tmp = *this;
    const int next = _iter_to_t->sameColorNext();
    if (next < 0)
    {
      _iter_to_t = _ptr_to_seq->_data.end();
      return tmp;
    }
    _iter_to_t = ContainerIterator( &(_ptr_to_seq->_data[next]));
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
