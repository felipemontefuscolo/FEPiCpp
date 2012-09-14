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

#ifndef FEPIC_TIMER_HPP
#define FEPIC_TIMER_HPP


#include <cstdio>
#include <cstring>
#include <list>

#ifdef FEP_HAS_OPENMP
#  include <omp.h>
#else
#  include <ctime>
#endif

namespace __FepTimer
{
  struct Item {
    Item(double e, const char* n) : elapsed(e) {
      strcpy(fname, n);
    }

    Item() : elapsed(0) {
      //fname[0] = '\n';
    };
    double elapsed;
    char fname[256];
  };

}

class Timer
{
public:
  typedef __FepTimer::Item Item;
  typedef typename std::list<Item> List;
  typedef typename List::iterator Iter;
  typedef typename List::const_iterator CIter;

  enum {LimitPushBacks = 64};

  explicit Timer()
  {
    #ifdef FEP_HAS_OPENMP
    _method = "OpenMp";
    #else
    _method = "ctime";
    #endif
  }

  void restart()
  {
    #ifdef FEP_HAS_OPENMP
    _elapsed = omp_get_wtime();
    #else
    _temp = clock();
    #endif
  }

  double elapsed(const char* fname = "", bool print_now = false)
  {
    char  buff[256];

    #ifdef FEP_HAS_OPENMP
    _elapsed = omp_get_wtime() - _elapsed;
    #else
    _elapsed = static_cast<double>( clock()  - _temp)/(1.*CLOCKS_PER_SEC);
    #endif
    sprintf(buff, "elapsed %8.3fs in : %s\n", static_cast<float>(_elapsed), fname);

    if (print_now)
      printf("%s", buff);

    if (_list.size() <= LimitPushBacks)
      _list.push_back(Item(_elapsed, buff));
  
    return _elapsed;
  }

  void printTimes() const
  {
    for (CIter it= _list.begin(); it != _list.end(); ++it)
    {
      printf("%s", it->fname);
    }
  }

  void addItems(Timer const& other)
  {
    if (other._list.empty())
      return;
    for (CIter it= other._list.begin(); it != other._list.end(); ++it)
      (this->_list).push_back( *it );
  }
  
  void printMethod() const
  {
    printf("%s\n", _method);
  }

protected:
  double      _elapsed;
  const char* _method;
  List	      _list;
  #ifndef FEP_HAS_OPENMP
  clock_t     _temp;
  #endif


};


#endif

