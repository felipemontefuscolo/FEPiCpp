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
  virtual void pushIncidCell(int cid, char pos) = 0;
  virtual int  numConnectedComps() const = 0;
  virtual bool isSingular() const = 0;
  virtual void replacesIncidCell(int cid, int cid_subs, char cid_subs_pos) = 0;
  virtual void getIthIncidCell(int ith, int &ic, int &pos) const = 0;
  virtual void clearIncidences() = 0;
  
  // inherited from CellElement
  virtual int getIncidCell() const = 0;
  virtual int getPosition() const = 0;
  virtual void setIncidCell(int icell_id) = 0;
  virtual void setPosition(int pos) = 0;

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
  explicit PointX(double const* coord, int ic=-1, char pos=-1) : _icell(ic), _icell_pos(pos)
  {
    for (int i = 0; i < spacedim; ++i)
      _coord[i] = coord[i];
  }

  /** Construtor.
  */
  PointX()  : _icell(-1), _icell_pos(-1) {}
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

  /** add a singular cell if it does not yet exist.
   *  @param cid cell id.
   *  @param pos position of this node on the cell cid.
   */ 
  void pushIncidCell(int cid, char pos)
  {
    FEPIC_CHECK(cid>=0 && pos >=0, "Point::pushIncidCell: invalid argument", std::invalid_argument);
    
    if (cid==this->_icell) // nothing to do
      return;
    else
    if (this->_icell<0)
    {
      this->_icell = cid;
      this->_icell_pos = pos;
      return;
    }
    else
    {
      std::list<std::pair<int,char> >::iterator it = std::find_if(_incidences.begin(), _incidences.end(), _RemoveCriteria(cid));
        
      if (it == _incidences.end())
        _incidences.push_back(std::make_pair(cid,pos));
    }
    
  }

  /** @return number of connected component incident to this node.
   */
  int numConnectedComps() const
  {
    return _incidences.size()+(_icell>=0);
  }

  bool isSingular() const
  {
    if (this->_icell < 0)
      return _incidences.size() > 1;
    else
      return !_incidences.empty();
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

  /** consistently replaces (or removes) a incident cell id (singular or not).
   * @param cid id of the incident cell of this node that will be replaced/removed.
   * @param cid_subs id of an incident cell that will replace the cid. if cid_subs<0,
   * then cid is just removed.
   */ 
  void replacesIncidCell(int cid, int cid_subs, char cid_subs_pos)
  {
    if (this->_icell == cid)
    {
      if (cid_subs < 0) // remove
      {
        if (_incidences.empty())
        {
          this->_icell = -1;
          return;
        }
        else
        {
          this->_icell     = _incidences.front().first;
          this->_icell_pos = _incidences.front().second;
          _incidences.pop_front();
          return;
        }
      }
      else // replaces
      {
        this->_icell = cid_subs;
        this->_icell_pos = cid_subs_pos;
        return;
      }
    }
    else
    {
      std::replace_if(_incidences.begin(), _incidences.end(), _RemoveCriteria(cid), std::make_pair(cid_subs,cid_subs_pos));
    }
    
  } // end replacesIncidCell()


  void getIthIncidCell(int ith, int &ic, int &pos) const
  {
    //FEPIC_CHECK((unsigned)ith <= _incidences.size(), "invalid index", std::invalid_argument);
    
    if (ith == 0)
    {
      ic = _icell;
      return;
    }
    
    --ith;
    
    std::list<std::pair<int,char> >::const_iterator it = _incidences.begin();
    std::list<std::pair<int,char> >::const_iterator et = _incidences.end();
    
    for (int k=0; k<ith; ++k)
    {
      ++it;
      if (it==et)
      {
        printf("ERROR: getIthIncidCell: invalid index\n");
        printf("ith = %d  _incidences.size()=%d\n", ith, (int)_incidences.size());
        throw;          
      }
    }
    ic = (*it).first;
    pos = (*it).second;
    
  }

  void clearIncidences()
  {
    this->_incidences.clear();
    this->_icell = -1;
  }

  // --- inherited from CellElement ----- //
  
  virtual int getIncidCell() const
  {
    return _icell;
  }
  virtual int getPosition() const
  {
    return _icell_pos;
  }
  virtual void setIncidCell(int icell_id)
  {
    _icell = icell_id;
  }
  virtual void setPosition(int pos)
  {
    _icell_pos = pos;
  }

  // --- inherited from CellElement ----- //



  /// Destrutor.
  ~PointX() {}

protected:
  Real _coord[sdim];

  int   _icell;
  char  _icell_pos;   // facet lid of incident cell
  char  _padd[3]; // supondo sizeof(int)=4 e sizeof(char)=1

  // main icell not included .. this list is for singular nodes only
  // the reason for this is to save memory.
  //                  iC  poiC
  std::list<std::pair<int,char> > _incidences; // for singular nodes
};








#endif
