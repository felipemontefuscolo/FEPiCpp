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
#include <vector>

class Point : public CellElement
{
protected:
  Real  m_coord[3]; // x,y,z
  std::list<std::pair<int,char> > m_incidences; // for singular nodes  
  
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

  static const EMshTag msh_tag = MSH_PNT;

  /** Construtor.
  *  @param coord um vetor com Dim elementos que armazena a coordenada.
  *  @param label seu rótulo.
  */
  template<class VectorT>
  Point(VectorT coord, int spacedim, int ic=-1, int pos=-1, int tag=0, int flags=0, int stat=0) : CellElement(ic, pos, tag, flags, stat)
  {
    for (int i = 0; i < spacedim; ++i)
      m_coord[i] = coord[i];
  }

  /** Construtor.
  */
  Point() {}

  static Point* create()
  {
    return new Point();
  }

  /** Define a coordenada deste ponto.
  *  @param coord um vetor com a coordenada.
  */
  void setCoord(Real const* coord, int spacedim)
  {
    for (int i = 0; i < spacedim; ++i)
      m_coord[i] = coord[i];
  }

  /** Retorna a coordenada deste ponto em coord.
  *  @param[out] coord a coordenada.
  */
  void getCoord(Real *coord, int spacedim) const
  {
    for (int i = 0; i < spacedim; ++i)
      coord[i] = m_coord[i];
  }

  Real const* getCoord() const
  {
    return &m_coord[0];
  }

  /** Retorna a i-ésima componente da coordenada
  */
  Real getCoord(int i) const
  {
    return m_coord[i];
  }

  /** add a singular cell if it does not yet exist.
   *  @param cid cell id.
   *  @param pos position of this node on the cell cid.
   */ 
  void pushIncidCell(int cid, int pos)
  {
    FEPIC_CHECK(cid>=0 && pos >=0, "Point::pushIncidCell: invalid argument", std::invalid_argument);
    
    if (cid==getIncidCell()) // nothing to do
      return;
    else
    if (getIncidCell()<0)
    {
      setIncidCell(cid);
      setPosition(pos);
      return;
    }
    else
    {
      std::list<std::pair<int,char> >::iterator it = std::find_if(m_incidences.begin(), m_incidences.end(), fi_RemoveCriteria(cid));
        
      if (it == m_incidences.end())
        m_incidences.push_back(std::make_pair(cid,static_cast<char>( pos )));
    }
    
  }

  /** @return number of connected component incident to this node.
   */
  int numConnectedComps() const
  {
    return m_incidences.size()+(getIncidCell()>=0);
  }

  bool isSingular() const
  {
    if (getIncidCell() < 0)
      return m_incidences.size() > 1;
    else
      return !m_incidences.empty();
  }

  class fi_RemoveCriteria { public:
    fi_RemoveCriteria(int id) : m_id_to_compare(id) {};
    bool operator() (std::pair<int,char> const& p) {return p.first==m_id_to_compare;};
    bool apply (std::pair<int,char> const& p) {return p.first==m_id_to_compare;};
    int m_id_to_compare;
  };

  /** consistently replaces (or removes) a incident cell id (singular or not).
   * @param cid id of the incident cell of this node that will be replaced/removed. If cid
   *         is not found in node's incidentes cell list, the function return false.
   * @param cid_subs id of an incident cell that will replace the cid. if cid_subs<0,
   * then cid is just removed.
   * @return true if cid was found and removed/replaced, false otherwise.
   */ 
  bool replacesIncidCell(int cid, int cid_subs, int cid_subs_pos)
  {
    if (getIncidCell() == cid)
    {
      if (cid_subs < 0) // remove
      {
        if (m_incidences.empty())
        {
          setIncidCell(-1);
          setPosition(-1);
          return true;
        }
        else
        {
          setIncidCell(m_incidences.front().first);
          setPosition(m_incidences.front().second);
          m_incidences.pop_front();
          return true;
        }
      }
      else // replaces
      {
        setIncidCell(cid_subs);
        setPosition(cid_subs_pos);
        return true;
      }
    }
    else
    {
      std::list<std::pair<int,char> >::iterator first = m_incidences.begin();
      std::list<std::pair<int,char> >::iterator last = m_incidences.end();
      fi_RemoveCriteria pred = fi_RemoveCriteria(cid);
      
      for (; first != last; ++first)
        if (pred.apply(*first))
        {
          *first = std::make_pair(cid_subs, static_cast<char>(cid_subs_pos) );
          return true;
        }
      //std::replace_if(m_incidences.begin(), m_incidences.end(), fi_RemoveCriteria(cid), std::make_pair(cid_subs, static_cast<char>(cid_subs_pos) ));
    }
    return false;
    
  } // end replacesIncidCell()


  void getIthIncidCell(int ith, int &ic, int &pos) const
  {
    //FEPIC_CHECK((unsigned)ith <= m_incidences.size(), "invalid index", std::invalid_argument);
    
    if (ith == 0)
    {
      ic = getIncidCell();
      pos = getPosition();
      return;
    }
    
    --ith;
    
    std::list<std::pair<int,char> >::const_iterator it = m_incidences.begin();
    std::list<std::pair<int,char> >::const_iterator et = m_incidences.end();
    
    for (int k=0; k<ith; ++k)
    {
      ++it;
      if (it==et)
      {
        printf("ERROR: getIthIncidCell: invalid index\n");
        printf("ith = %d  m_incidences.size()=%d\n", ith, (int)m_incidences.size());
        throw;          
      }
    }
    ic = (*it).first;
    pos = (*it).second;
    
  }

  /** @param v the vector where will be places incid cells.
   *         v[2*i]   = i-th incident cell
   *         v[2*i+1] = position in i-th incident cell.
   */ 
  void getAllIncidences(std::vector<int> & v) const
  {
    if (!v.empty())
      v.clear();
    v.push_back(this->getIncidCell());
    v.push_back(this->getPosition());
    
    std::list<std::pair<int,char> >::const_iterator it = m_incidences.begin();
    std::list<std::pair<int,char> >::const_iterator et = m_incidences.end();
    
    for (; it != et; ++it)
    {
      v.push_back((*it).first);
      v.push_back((*it).second);
    }
  }

  void clearIncidences()
  {
    this->m_incidences.clear();
    setIncidCell(-1);
    this->setAsBoundary(false);
  }

  /// @return true if this point is a boundary point, false otherwise.
  bool inBoundary() const
  {
    return CellElement::m_status & mk_inboundary;
  }

  /// @param ib set ib=true(false) to (un)set this point as boundary point.
  void setAsBoundary(bool ib)
  {
    CellElement::m_status = ib ? (m_status | mk_inboundary) : (m_status & (~mk_inboundary));
  }  


  // ----------------------------------------------- inherited from CellElement ----- //

  /// Destrutor.
  ~Point() {}

  void setAllMembers(Real const* coord, int spacedim, int const* ic, int const* pos, int const* tag,
                    int const* flags, int const* stat, std::list<std::pair<int,char> > * incidences)
  {
    if (coord != NULL)
    {
      FEPIC_CHECK(static_cast<unsigned> (spacedim-1) < 3, "invalid argument", std::runtime_error);
      for (int i = 0; i < spacedim; ++i)
        m_coord[i] = coord[i];
    }
    if (ic != NULL)
      setIncidCell(*ic);
    if (pos != NULL)
      setPosition(*pos);
    if (tag != NULL)
      setTag(*tag);
    if (flags != NULL)
      setFlags(*flags);
    if (stat != NULL)
      m_status = *stat;
    if (incidences != NULL)
      m_incidences = *incidences;
  }  

  enum Masks
  {
    mk_inboundary = (1<<0)
  };

};








#endif
