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

#ifndef FEPIC_QUADRATURE_HPP
#define FEPIC_QUADRATURE_HPP

#include "conf/directives.hpp"
#include "../mesh/enums.hpp"

void igetQuadrPtsHypercube(int n, int dim, Real **points, Real *weights, int &nm_points);

class Quadrature
{
public:

  Quadrature()
  {
    
  }

  virtual Real const* point(int qp) const = 0;
  virtual void getPoint(int qp, Real *x) const = 0;
  virtual Real weight(int qp) const = 0;

  virtual void setOrder(int n) = 0;
  virtual int  maxOrder() const = 0;

  int numPoints() const
  {
    return m_nm_points;
  }

  int order() const
  {
    return m_order;
  }

  static
  Quadrature* create(ECellType ct);

  static
  Quadrature* create(ECellClass cc);

  virtual ~Quadrature() {};
protected:
  int m_nm_points;
  int m_order;
};


template<ECellType CType, int Dim_, int MaxNumQPoints_, int MaxOrder_>
class Quadrature_ : public Quadrature
{
public:

  enum{Dim=Dim_, MaxOrder=MaxOrder_, MaxNumQPoints=MaxNumQPoints_};

  Quadrature_(int order_ = 1) : Quadrature()
  {
    this->setOrder(order_);
  }

  Real const* point(int qp) const
  {
    return this->m_points[qp];
  }

  void getPoint(int qp, Real *x) const
  {
    for (int i = 0; i < Dim; ++i)
      x[i] = this->m_points[qp][i];
  }

  Real weight(int qp) const
  {
    return this->m_weights[qp];
  }

  int  maxOrder() const
  {
    return MaxOrder;
  }

  void setOrder(int n);

  ~Quadrature_() {}

protected:
  Real m_points[MaxNumQPoints][Dim!=0?Dim:1];
  Real m_weights[MaxNumQPoints];

};


// =============================================

//typedef Quadrature_<EDGE2,        1, 5  > Quadrature_EDGE;
//typedef Quadrature_<TRIANGLE3,    2, 25 > Quadrature_TRIANGLE;
//typedef Quadrature_<QUADRANGLE4,  2, 25 > Quadrature_QUADRANGLE;
//typedef Quadrature_<TETRAHEDRON4, 3, 46 > Quadrature_TETRAHEDRON;
//typedef Quadrature_<HEXAHEDRON8,  3, 125> Quadrature_HEXAHEDRON;

#define Quadrature_POINT       Quadrature_<POINT1,       0, 1  , 10>
#define Quadrature_EDGE        Quadrature_<EDGE2,        1, 5  , 9 >
#define Quadrature_TRIANGLE    Quadrature_<TRIANGLE3,    2, 25 , 10>
#define Quadrature_QUADRANGLE  Quadrature_<QUADRANGLE4,  2, 25 , 9 >
#define Quadrature_TETRAHEDRON Quadrature_<TETRAHEDRON4, 3, 46 , 8 >
#define Quadrature_HEXAHEDRON  Quadrature_<HEXAHEDRON8,  3, 125, 9 >



#endif


