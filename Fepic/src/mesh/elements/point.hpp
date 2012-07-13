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

#ifndef FEPIC_POINT_HPP
#define FEPIC_POINT_HPP

#include <list>
#include <utility>
#include "Fepic/src/util/assert.hpp"
#include <algorithm>


class Point : public CellElement
{
public:

  enum { dim=0,
         n_vertices=1,
         n_nodes=1,
         n_facets=0,
         n_corners=0,
         n_vertices_per_facet=0,
         n_vertices_per_corner=0,
         n_nodes_per_facet=0,
         n_nodes_per_corner=1};

  virtual int  spaceDim() const = 0;
  virtual void setCoord(Real const*const coord) = 0;
  virtual void getCoord(Real *coord) const = 0;
  virtual Real getCoord(int const i) const = 0;
  virtual Real const* getCoord() const = 0;
  virtual void pushIncidCell(int cid, int pos) = 0;
  virtual int singularity() const = 0;
  virtual void replacesIncidCell(int cid, int cid_subs, int cid_subs_pos) = 0;
  virtual void getIthIncidCell(int ith, int &ic, int &pos) const = 0;

  virtual ~Point(){};
};


template<int sdim>
class PointX : public Point
{
public:
  enum {spacedim=sdim};

  static const EMshTag msh_tag = MSH_PNT;

  //typedef Eigen::Matrix<Real, spacedim, 1> VecT;

  /** Construtor.
  *  @param coord um vetor com Dim elementos que armazena a coordenada.
  *  @param label seu rótulo.
  */
  PointX(double const*const coord)
  {
    for (int i = 0; i < spacedim; ++i)
      _coord[i] = coord[i];
  }

  /** Construtor.
  */
  PointX() {}
  // construtor de cópia não necessário, pois não há nenhum ponteiro.

  /** @return a dimensão do espaço.
  */
  int spaceDim() const
  {
    return spacedim;
  }

  /** Define a coordenada deste ponto.
  *  @param coord um vetor com a coordenada.
  */
  void setCoord(Real const*const coord)
  {
    for (int i = 0; i < spacedim; ++i)
      _coord[i] = coord[i];
  }

  /** Retorna a coordenada deste ponto em coord.
  *  @param[out] coord a coordenada.
  */
  void getCoord(Real *coord) const
  {
    for (int i = 0; i < spacedim; ++i)
      coord[i] = _coord[i];
  }

  Real const* getCoord() const
  {
    return &_coord[0];
  }

  /** Retorna a i-ésima componente da coordenada
  */
  Real getCoord(int i) const
  {
    return _coord[i];
  }

  ///** Retorna a distância até a coordenada p.
  //*/
  //Real getDistance(Eigen::VectorXr const& p) const
  //{
    //return sqrt( (_coord-p).dot(_coord-p) );
  //}

  ///** Retorna a distância até o ponto p.
  //*/
  //Real getDistance(Point const*const p) const
  //{
    //PointX const*const pp = static_cast<PointX const*const>(p);
    //return sqrt( (_coord - pp->_coord).dot(_coord - pp->_coord) );
  //}

  void pushIncidCell(int cid, int pos)
  {
    FEPIC_CHECK(cid>=0 && pos >=0, "Point::pushIncidCell: invalid argument", std::invalid_argument);
    
    if (cid == getIncidCell())
      return;
    
    if (_extra_icells.empty())
    {
      _extra_icells.push_back(std::pair<int,char>(cid,(char)pos));
      return;
    }
    
    std::list<std::pair<int,char> >::iterator it = std::find_if(_extra_icells.begin(), _extra_icells.end(), _RemoveCriteria(cid));
    
    if (it == _extra_icells.end())
      _extra_icells.push_back(std::pair<int,char>(cid,(char)pos));
    
  }

  /** @return number of connected component incident to this node minus one.
   */
  int singularity() const
  {
    return _extra_icells.size();
  }

  class _RemoveCriteria { public:
    _RemoveCriteria(int id) : _id_to_compare(id) {};
    bool operator() (std::pair<int,char> const& p) {return p.first==_id_to_compare;};
    int _id_to_compare;
  };

  //class _PairIntChar { public:
  //  bool operator() (std::pair<int,char> const& p, int const& a) {return p.first<a};
  //  bool operator() (int const& a, std::pair<int,char> const& p) {return a<p.first};
  //  bool operator() (std::pair<int,char> const& p, std::pair<int,char> const& a) {return p.first<a.first};
  //};

  /** replaces (or removes) a incident cell id (singular or not).
   * @param cid id of the incident cell of this node that will be replaced/removed.
   * @param cid_subs id of an incident cell that will replace the cid. if cid_subs<0,
   * then cid is just removed.
   */ 
  void replacesIncidCell(int cid, int cid_subs, int cid_subs_pos)
  {
    if (cid == this->getIncidCell())
    {
      if (cid_subs < 0) // remove
      {
        if (_extra_icells.empty())
        {
          this->setIncidCell(cid_subs);
          this->setPosition(cid_subs_pos);
        }
        else
        {
          this->setIncidCell(_extra_icells.front().first);
          this->setPosition(_extra_icells.front().second);
          _extra_icells.pop_front();
        }
      }
      else
      {
        this->setIncidCell(cid_subs);
        this->setPosition(cid_subs_pos);
      }
    }
    else
    {
      if (cid_subs < 0)
        _extra_icells.remove_if(_RemoveCriteria(cid));
      else
      {
        std::replace_if(_extra_icells.begin(), _extra_icells.end(), _RemoveCriteria(cid), std::pair<int,char>(cid_subs,cid_subs_pos));
      }
    }
    
  } // end replacesIncidCell()


  void getIthIncidCell(int ith, int &ic, int &pos) const
  {
    //FEPIC_CHECK((unsigned)ith <= _extra_icells.size(), "invalid index", std::invalid_argument);
    if (ith == 0)
    {
      ic = getIncidCell();
      pos = getPosition();
      return;
    }
    else
    {
      std::list<std::pair<int,char> >::const_iterator it = _extra_icells.begin();
      std::list<std::pair<int,char> >::const_iterator et = _extra_icells.end();
      int k = 1;
      for (; it != et; ++it)
      {
        if (k==ith)
        {
          ic = (*it).first;
          pos = (*it).second;
          return;
        }
      }
      FEPIC_CHECK((unsigned)ith <= _extra_icells.size(), "invalid index", std::invalid_argument);
    }
  }

  /** Destrutor.
  */
  ~PointX() {}

protected:
  Real _coord[sdim];

  // warning: main icell not included!
    //                 iC  poiC
  std::list<std::pair<int,char> > _extra_icells; // for singular nodes
};








#endif
