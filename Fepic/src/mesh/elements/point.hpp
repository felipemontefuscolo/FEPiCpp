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


class Point : public _Labelable
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
  
  virtual int spaceDim() const = 0;
  virtual void setCoord(Real const*const coord) = 0;
  virtual void getCoord(Real *coord) const = 0;
  virtual Real getCoord(int const i) const = 0;
  virtual int getIncidCell() const = 0;
  virtual int getPosition() const = 0;
  virtual void setIncidCell(int const icell_id) = 0;
  virtual void setPosition(int const pos) = 0;
  virtual Real const* getCoord() const = 0;
  
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

  int getIncidCell() const
  {
    return _icell;
  }
  
  int getPosition() const
  {
    return _icell_pos;
  }

  void setIncidCell(int const icell_id)
  {
    _icell = icell_id;
  }

  // facet lid of incident cell
  void setPosition(int const pos)
  {
    _icell_pos = pos;
  }


  /** Destrutor.
  */
  ~PointX() {}

protected:
  
  Real _coord[sdim];
  int  _icell;      // incident cells
  char _icell_pos;

};








#endif
