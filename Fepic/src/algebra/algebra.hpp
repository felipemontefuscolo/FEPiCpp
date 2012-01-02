#ifndef FEPIC_ALGEBRA_HPP
#define FEPIC_ALGEBRA_HPP

#include <vector>
#include <tr1/array>

// Dim = 1 vector
// Dim = 2 matrix
template<class T, int Dim=1, bool Row_major=true, class A = std::allocator<T> >
class FepArray : private std::vector<T, A>
{
  typedef std::vector<T, A> _Base;  
  
public:

  typedef typename _Base::reference               reference;
  typedef typename _Base::const_reference         const_reference;
  typedef typename _Base::iterator                iterator;
  typedef typename _Base::const_iterator          const_iterator;
  typedef typename _Base::size_type               size_type;
  typedef typename _Base::difference_type         difference_type;
  typedef typename _Base::value_type              value_type;
  typedef typename _Base::allocator_type          allocator_type;
  typedef typename _Base::pointer                 pointer;
  typedef typename _Base::const_pointer           const_pointer;
  typedef typename _Base::reverse_iterator        reverse_iterator;
  typedef typename _Base::const_reverse_iterator  const_reverse_iterator;

private:
  std::tr1::array<size_type,Dim> _dims;
  
public:
  FepArray() : _Base() {};
  explicit FepArray(size_type m) : _Base(m)
  {
    _dims.at(0) = m;
  };
  explicit FepArray(size_type m, size_type n) : _Base(m*n)
  {
    _dims.at(0) = m;
    _dims.at(1) = n;
  };   
  FepArray(const FepArray& x) : _Base(x.size())
  {
    (*this) = x;
  }  
  
  size_type size() const
  { return _Base::size(); }

  void resize ( size_type m)
  { _Base::resize(m); }

  void resize ( size_type m, size_type n)
  { _Base::resize(m*n); }

  void reserve(size_type a)
  { _Base::reserve(a); }
  
  size_type capacity() const
  { return _Base::capacity(); }
  
  reference back()
  { return _Base::back(); }
  
  const_reference back ( ) const
  { return _Base::back(); }
  
  void pop_back()
  { _Base::pop_back(); }
  
  void push_back ( const T& x )
  { _Base::push_back(x); }

  iterator begin() { return _Base::begin(); }
  const_iterator begin() const { return _Base::begin(); }
  iterator end() { return _Base::end(); }
  const_iterator end() const { return _Base::end(); }
  reverse_iterator rbegin() { return _Base::rbegin(); }
  const_reverse_iterator rbegin() const { return _Base::rbegin(); }
  reverse_iterator rend() { return _Base::rend(); }
  const_reverse_iterator rend() const { return _Base::rend(); }

  void clear()
  { _Base::clear(); }

  pointer
  data()
  { return _Base::data(); }
  
  const_pointer
  data() const
  { return _Base::data(); }  
    
  reference operator[] ( size_type n )
  { return _Base::operator[] (n); };
  
  const_reference operator[] ( size_type n ) const
  { return _Base::operator[] (n); };
};



#endif


